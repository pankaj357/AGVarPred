import pandas as pd
import numpy as np
import gc, re, joblib, os
from glob import glob

from sklearn.model_selection import GroupKFold   # 🔥 UPDATED
from sklearn.metrics import (
    roc_auc_score, f1_score, precision_score, recall_score,
    roc_curve, average_precision_score,
    brier_score_loss, confusion_matrix, matthews_corrcoef,
    precision_recall_curve
)
from sklearn.dummy import DummyClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.calibration import calibration_curve
from sklearn.preprocessing import RobustScaler

from lightgbm import LGBMClassifier
import optuna
import shap
import matplotlib.pyplot as plt

# =========================
# OUTPUT DIR
# =========================

OUTDIR = "final_model_output"
os.makedirs(OUTDIR, exist_ok=True)

# =========================
# CLASS IMBALANCE
# =========================

SCALE_POS_WEIGHT = None  # computed dynamically after loading data

# =========================
# CLEAN NAME
# =========================

def clean_name(s):
    s = str(s)
    s = s.replace(":", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r'[^A-Za-z0-9_]', '', s)
    s = re.sub(r'_+', '_', s)
    return s

# =========================
# LOAD DATA
# =========================

def load_dataset(path):
    files = glob(f"{path}/*.parquet")
    dfs = []

    for i, f in enumerate(files):
        dfs.append(pd.read_parquet(f))
        if i % 500 == 0:
            print(f"{path}: {i}/{len(files)}")

    df = pd.concat(dfs, ignore_index=True)
    del dfs; gc.collect()

    df.columns = [clean_name(c) for c in df.columns]

    return df

print("🔄 Loading datasets...")

train_df = load_dataset("final_dataset_parts_train")
cal_df   = load_dataset("final_dataset_parts_cal")
test_df  = load_dataset("final_dataset_parts_test")

# =========================
# FEATURES
# =========================

with open("feature_selection_output_nested/selected_features.txt") as f:
    selected_features = [clean_name(x.strip()) for x in f]

selected_features = [f for f in selected_features if f in train_df.columns]

print("Selected features:", len(selected_features))

X_train = train_df[selected_features]
y_train = train_df['label']
SCALE_POS_WEIGHT = (y_train == 0).sum() / (y_train == 1).sum()

# 🔥 ADDED (GENE GROUPS)
groups = train_df["gene"].values

X_cal = cal_df[selected_features]
y_cal = cal_df['label']

X_test = test_df[selected_features]
y_test = test_df['label']

# =========================
# IMPUTER
# =========================

imputer_full = X_train.median()

# =========================
# OPTUNA
# =========================

print("\n🔍 Optuna tuning...")

def objective(trial):

    model = LGBMClassifier(
        n_estimators=trial.suggest_int("n_estimators", 200, 600),
        learning_rate=trial.suggest_float("learning_rate", 0.01, 0.1),
        num_leaves=trial.suggest_int("num_leaves", 20, 150),
        max_depth=trial.suggest_int("max_depth", 3, 12),
        scale_pos_weight=SCALE_POS_WEIGHT,
        random_state=42,
        verbose=-1,
        n_jobs=-1
    )

    kf = GroupKFold(n_splits=5)   # 🔥 UPDATED
    scores = []

    for train_idx, val_idx in kf.split(X_train, y_train, groups):  # 🔥 UPDATED

        X_tr, X_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_tr, y_val = y_train.iloc[train_idx], y_train.iloc[val_idx]

        imputer = X_tr.median()
        X_tr_imp = X_tr.fillna(imputer)
        X_val_imp = X_val.fillna(imputer)

        scaler = RobustScaler()
        X_tr_s = pd.DataFrame(
            scaler.fit_transform(X_tr_imp),
            columns=X_tr_imp.columns,
            index=X_tr_imp.index
        )
        X_val_s = pd.DataFrame(
            scaler.transform(X_val_imp),
            columns=X_val_imp.columns,
            index=X_val_imp.index
        )

        model.fit(X_tr_s, y_tr)
        probs = model.predict_proba(X_val_s)[:,1]

        scores.append(average_precision_score(y_val, probs))

        del X_tr, X_val
        gc.collect()

    return np.mean(scores)

study = optuna.create_study(
    direction="maximize",
    sampler=optuna.samplers.TPESampler(seed=42)
)
study.optimize(objective, n_trials=100)

joblib.dump(study, f"{OUTDIR}/optuna_study.pkl")

best_params = study.best_params

# =========================
# BASELINE
# =========================

dummy = DummyClassifier(strategy="stratified")
dummy.fit(X_train.fillna(imputer_full), y_train)

baseline_probs = dummy.predict_proba(
    X_test.fillna(imputer_full)
)[:, 1]

baseline_auc = roc_auc_score(y_test, baseline_probs)
print(f"Baseline AUC: {baseline_auc:.4f}")

# =========================
# FINAL MODEL
# =========================

X_train_imp = X_train.fillna(imputer_full)
X_cal_imp   = X_cal.fillna(imputer_full)
X_test_imp  = X_test.fillna(imputer_full)

scaler_full = RobustScaler()
X_train_f = pd.DataFrame(
    scaler_full.fit_transform(X_train_imp),
    columns=X_train_imp.columns,
    index=X_train_imp.index
)
X_cal_f = pd.DataFrame(
    scaler_full.transform(X_cal_imp),
    columns=X_cal_imp.columns,
    index=X_cal_imp.index
)
X_test_f = pd.DataFrame(
    scaler_full.transform(X_test_imp),
    columns=X_test_imp.columns,
    index=X_test_imp.index
)

final_model = LGBMClassifier(
    **best_params,
    scale_pos_weight=SCALE_POS_WEIGHT,
    random_state=42,
    verbose=-1,
    n_jobs=-1
)
final_model.fit(X_train_f, y_train)

# =========================
# CALIBRATION
# =========================

probs_cal = final_model.predict_proba(X_cal_f)[:,1]

iso = IsotonicRegression(out_of_bounds='clip')
iso.fit(probs_cal, y_cal)

# =========================
# THRESHOLD
# =========================

probs_cal_calibrated = iso.transform(probs_cal)

thresholds_grid = np.arange(0.05, 0.95, 0.01)
final_thresh = max(
    thresholds_grid,
    key=lambda t: f1_score(y_cal, (probs_cal_calibrated > t).astype(int))
)

print(f"\n✅ Threshold tuned on CAL set: {final_thresh:.2f}")

np.savetxt(f"{OUTDIR}/cv_thresholds.txt", [final_thresh])

# =========================
# TEST
# =========================

probs_test_raw = final_model.predict_proba(X_test_f)[:,1]
probs_test = iso.transform(probs_test_raw)

preds = (probs_test > final_thresh).astype(int)

roc_auc = roc_auc_score(y_test, probs_test)
auprc = average_precision_score(y_test, probs_test)
f1 = f1_score(y_test, preds)
precision = precision_score(y_test, preds)
recall = recall_score(y_test, preds)
mcc = matthews_corrcoef(y_test, preds)
brier = brier_score_loss(y_test, probs_test)

print("\n📊 FINAL TEST PERFORMANCE")
print("ROC-AUC:", roc_auc)

# =========================
# ECE
# =========================

def compute_ece(y_true, y_prob, n_bins=10):
    bins = np.linspace(0, 1, n_bins + 1)
    binids = np.digitize(y_prob, bins) - 1
    ece = 0
    for i in range(n_bins):
        mask = binids == i
        if np.sum(mask) > 0:
            acc = np.mean(y_true[mask])
            conf = np.mean(y_prob[mask])
            ece += np.abs(acc - conf) * np.sum(mask) / len(y_true)
    return ece

ece = compute_ece(y_test.values, probs_test)

# =========================
# SAVE OUTPUTS
# =========================

pd.DataFrame({
    "y_true": y_test.values,
    "raw_prob": probs_test_raw,
    "cal_prob": probs_test
}).to_csv(f"{OUTDIR}/test_predictions.csv", index=False)

np.savetxt(f"{OUTDIR}/confusion_matrix.txt",
           confusion_matrix(y_test, preds), fmt='%d')

# ROC
fpr, tpr, _ = roc_curve(y_test, probs_test)
plt.plot(fpr, tpr)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.savefig(f"{OUTDIR}/roc_curve.png")
plt.close()

# PR
prec, rec, _ = precision_recall_curve(y_test, probs_test)
plt.plot(rec, prec)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig(f"{OUTDIR}/pr_curve.png")
plt.close()

# Calibration
prob_true, prob_pred = calibration_curve(y_test, probs_test, n_bins=10)
plt.plot(prob_pred, prob_true, marker='o')
plt.plot([0,1],[0,1],'--')
plt.xlabel('Mean Predicted Probability')
plt.ylabel('Fraction of Positives')
plt.savefig(f"{OUTDIR}/calibration_curve.png")
plt.close()

# =========================
# SHAP
# =========================

X_sample = X_train_f.sample(5000, random_state=42)

explainer = shap.TreeExplainer(final_model)
shap_values = explainer.shap_values(X_sample)

if isinstance(shap_values, np.ndarray) and shap_values.ndim == 3:
    shap_values = shap_values[:, :, 1]
elif isinstance(shap_values, list):
    shap_values = shap_values[1]

shap.summary_plot(shap_values, X_sample, show=False)
plt.xlabel('SHAP Value')
plt.ylabel('Features')
plt.savefig(f"{OUTDIR}/shap_summary.png", bbox_inches='tight')
plt.close()

pd.DataFrame(shap_values, columns=X_sample.columns)\
  .to_csv(f"{OUTDIR}/shap_values_sample.csv", index=False)

pd.DataFrame({
    "feature": X_sample.columns,
    "mean_abs_shap": np.abs(shap_values).mean(axis=0)
}).sort_values("mean_abs_shap", ascending=False)\
 .to_csv(f"{OUTDIR}/feature_importance_shap.csv", index=False)

# =========================
# CI + METRICS
# =========================

rng = np.random.RandomState(42)
auc_scores = []

for _ in range(200):
    idx = rng.choice(len(y_test), len(y_test), replace=True)
    auc_scores.append(roc_auc_score(y_test.values[idx], probs_test[idx]))

ci_low = np.percentile(auc_scores, 2.5)
ci_high = np.percentile(auc_scores, 97.5)

pd.DataFrame({
    "metric": [
        "AUC","AUPRC","F1","Precision","Recall","MCC",
        "Brier","ECE","Baseline_AUC",
        "AUC_CI_low","AUC_CI_high"
    ],
    "value": [
        roc_auc, auprc, f1, precision, recall, mcc,
        brier, ece, baseline_auc,
        ci_low, ci_high
    ]
}).to_csv(f"{OUTDIR}/final_metrics.csv", index=False)

# =========================
# PIPELINE
# =========================

joblib.dump({
    "model": final_model,
    "calibrator": iso,
    "imputer": imputer_full,
    "scaler": scaler_full,
    "features": selected_features,
    "threshold": final_thresh
}, f"{OUTDIR}/final_pipeline.pkl")

print("\n🎯 FINAL PIPELINE COMPLETE")
