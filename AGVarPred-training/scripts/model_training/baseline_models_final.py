import pandas as pd
import numpy as np
import gc, re, joblib, os, json
from glob import glob

from sklearn.model_selection import GroupKFold
from sklearn.metrics import (roc_auc_score, f1_score, precision_score, recall_score,
average_precision_score, brier_score_loss,
confusion_matrix, matthews_corrcoef,
roc_curve, precision_recall_curve)

from sklearn.calibration import calibration_curve
from sklearn.isotonic import IsotonicRegression
from sklearn.preprocessing import RobustScaler
from sklearn.dummy import DummyClassifier

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

import optuna
from optuna.trial import TrialState
import matplotlib.pyplot as plt

# =========================
# OUTPUT DIR
# =========================

OUT_DIR = "baseline_output"
os.makedirs(OUT_DIR, exist_ok=True)

CHECKPOINT_PATH = os.path.join(OUT_DIR, "checkpoint.json")

def load_checkpoint():
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH, 'r') as f:
            return json.load(f)
    return {}

def save_checkpoint(cp):
    """Atomic write so a power cut doesn't corrupt the checkpoint."""
    tmp = CHECKPOINT_PATH + ".tmp"
    with open(tmp, 'w') as f:
        json.dump(cp, f, indent=2)
    os.replace(tmp, CHECKPOINT_PATH)

checkpoint = load_checkpoint()

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
    del dfs
    gc.collect()
    df.columns = [clean_name(c) for c in df.columns]
    return df

print("Loading datasets...")

train_df = load_dataset("final_dataset_parts_train")
cal_df   = load_dataset("final_dataset_parts_cal")
test_df  = load_dataset("final_dataset_parts_test")

# =========================
# FEATURES (optimal 40)
# =========================

with open("feature_selection_output_nested/selected_features.txt") as f:
    selected_features = [clean_name(x.strip()) for x in f]

selected_features = [f for f in selected_features
                     if f in train_df.columns]

print("Selected features:", len(selected_features))

X_train = train_df[selected_features]
y_train = train_df['label']
groups = train_df["gene"].values

X_cal = cal_df[selected_features]
y_cal = cal_df['label']

X_test = test_df[selected_features]
y_test = test_df['label']

# =========================
# IMPUTER + SCALER (full train)
# =========================

imputer_full = X_train.median()

X_train_imp = X_train.fillna(imputer_full)
X_cal_imp   = X_cal.fillna(imputer_full)
X_test_imp  = X_test.fillna(imputer_full)

scaler_full = RobustScaler()
X_train_s = pd.DataFrame(scaler_full.fit_transform(X_train_imp),
                        columns=X_train_imp.columns,
                        index=X_train_imp.index)
X_cal_s = pd.DataFrame(scaler_full.transform(X_cal_imp),
                      columns=X_cal_imp.columns,
                      index=X_cal_imp.index)
X_test_s = pd.DataFrame(scaler_full.transform(X_test_imp),
                       columns=X_test_imp.columns,
                       index=X_test_imp.index)

scale_pos_weight = (y_train == 0).sum() / (y_train == 1).sum()

# =========================
# OPTUNA TUNING (all models)
# =========================

def tune_model(name, model_cls, param_fn, extra_kwargs=None, n_trials=100):
    print(f"\nTuning {name}...")
    
    best_params_path = os.path.join(OUT_DIR, f"{name}_best_params.json")
    if checkpoint.get(f"{name}_tuned") and os.path.exists(best_params_path):
        print(f"  {name} already tuned. Loading best params.")
        with open(best_params_path, 'r') as f:
            return json.load(f)
    
    storage_path = os.path.abspath(os.path.join(OUT_DIR, "optuna.db"))
    storage_url = f"sqlite:///{storage_path}"
    
    study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=42),
        study_name=name,
        storage=storage_url,
        load_if_exists=True
    )
    
    n_completed = sum(1 for t in study.trials if t.state == TrialState.COMPLETE)
    remaining = max(0, n_trials - n_completed)
    print(f"  {name} study has {n_completed} completed trials; running {remaining} more.")
    
    if remaining > 0:
        def objective(trial):
            params = param_fn(trial)
            if extra_kwargs:
                params.update(extra_kwargs)
            
            model = model_cls(**params)
            kf = GroupKFold(n_splits=5)
            scores = []
            
            for train_idx, val_idx in kf.split(X_train, y_train, groups):
                X_tr, X_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
                y_tr, y_val = y_train.iloc[train_idx], y_train.iloc[val_idx]
                
                imp = X_tr.median()
                X_tr_i = X_tr.fillna(imp)
                X_val_i = X_val.fillna(imp)
                
                sc = RobustScaler()
                X_tr_sc = pd.DataFrame(
                    sc.fit_transform(X_tr_i),
                    columns=X_tr_i.columns,
                    index=X_tr_i.index
                )
                X_val_sc = pd.DataFrame(
                    sc.transform(X_val_i),
                    columns=X_val_i.columns,
                    index=X_val_i.index
                )
                
                model.fit(X_tr_sc, y_tr)
                probs = model.predict_proba(X_val_sc)[:, 1]
                scores.append(average_precision_score(y_val, probs))
                
                del X_tr, X_val
                gc.collect()
            
            return np.mean(scores)
        
        study.optimize(objective, n_trials=remaining)
    
    # Save outputs immediately so they survive a crash
    joblib.dump(study, os.path.join(OUT_DIR, f"{name}_optuna_study.pkl"))
    with open(best_params_path, 'w') as f:
        json.dump(study.best_params, f, indent=2)
    
    checkpoint[f"{name}_tuned"] = True
    save_checkpoint(checkpoint)
    
    print(f"Best {name}:", study.best_params)
    return study.best_params

# --- Random Forest ---

def rf_params(trial):
    return {
        "n_estimators": trial.suggest_int("n_estimators", 200, 600),
        "max_depth": trial.suggest_int("max_depth", 3, 30),
        "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
        "max_features": trial.suggest_categorical("max_features", ["sqrt", "log2"]),
    }

best_rf = tune_model("RandomForest", RandomForestClassifier, rf_params,
                    extra_kwargs={"class_weight": "balanced", "random_state": 42, "n_jobs": -1})

# --- Logistic Regression ---

def lr_params(trial):
    return {
        "C": trial.suggest_float("C", 1e-3, 100.0, log=True),
    }

best_lr = tune_model("LogisticRegression", LogisticRegression, lr_params,
                    extra_kwargs={"max_iter": 1000, "solver": "liblinear", "class_weight": "balanced", "random_state": 42})

# =========================
# FINAL MODELS
# =========================

models = {
    "RandomForest": RandomForestClassifier(**best_rf, class_weight="balanced", random_state=42, n_jobs=-1),
    "LogisticRegression": LogisticRegression(**best_lr, max_iter=1000, solver="liblinear",
                                            class_weight="balanced", random_state=42),
}

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

# =========================
# EVALUATION LOOP
# =========================

results = []

for name, model in list(models.items()):
    print(f"\n=== {name} ===")
    model_path = os.path.join(OUT_DIR, f"{name}_best.pkl")
    results_json_path = os.path.join(OUT_DIR, f"{name}_results.json")
    
    # --- Training (resumable) ---
    if checkpoint.get(f"{name}_trained") and os.path.exists(model_path):
        print(f"{name} already trained. Loading from disk...")
        model = joblib.load(model_path)
        models[name] = model
    else:
        print(f"Training {name} on full training set...")
        model.fit(X_train_s, y_train)
        joblib.dump(model, model_path)
        models[name] = model
        checkpoint[f"{name}_trained"] = True
        save_checkpoint(checkpoint)
    
    # --- Evaluation (resumable) ---
    if checkpoint.get(f"{name}_evaluated") and os.path.exists(results_json_path):
        print(f"{name} already evaluated. Loading results...")
        with open(results_json_path, 'r') as f:
            results.append(json.load(f))
        continue
    
    print(f"Evaluating {name}...")
    
    # --- Calibration ---
    probs_cal = model.predict_proba(X_cal_s)[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(probs_cal, y_cal)
    joblib.dump(iso, os.path.join(OUT_DIR, f"{name}_calibrator.pkl"))
    
    # --- Threshold ---
    probs_cal_cal = iso.transform(probs_cal)
    thresholds_grid = np.arange(0.05, 0.95, 0.01)
    final_thresh = max(
        thresholds_grid,
        key=lambda t: f1_score(y_cal, (probs_cal_cal > t).astype(int))
    )
    print(f"{name} threshold: {final_thresh:.2f}")
    np.savetxt(os.path.join(OUT_DIR, f"{name}_threshold.txt"), [final_thresh])
    
    # --- Test ---
    probs_test_raw = model.predict_proba(X_test_s)[:, 1]
    probs_test = iso.transform(probs_test_raw)
    preds = (probs_test > final_thresh).astype(int)
    
    roc_auc = roc_auc_score(y_test, probs_test)
    auprc = average_precision_score(y_test, probs_test)
    f1 = f1_score(y_test, preds)
    precision = precision_score(y_test, preds)
    recall = recall_score(y_test, preds)
    mcc = matthews_corrcoef(y_test, preds)
    brier = brier_score_loss(y_test, probs_test)
    ece = compute_ece(y_test.values, probs_test)
    
    print(f"{name} AUC: {roc_auc:.4f}")
    print(f"{name} AUPRC: {auprc:.4f}")
    
    # --- CI ---
    rng = np.random.RandomState(42)
    auc_scores = []
    for _ in range(200):
        idx = rng.choice(len(y_test), len(y_test), replace=True)
        auc_scores.append(roc_auc_score(y_test.values[idx], probs_test[idx]))
    ci_low = np.percentile(auc_scores, 2.5)
    ci_high = np.percentile(auc_scores, 97.5)
    
    # --- Save predictions ---
    pd.DataFrame({
        "y_true": y_test.values,
        "raw_prob": probs_test_raw,
        "cal_prob": probs_test
    }).to_csv(os.path.join(OUT_DIR, f"{name}_predictions.csv"), index=False)
    
    # --- Confusion matrix ---
    np.savetxt(
        os.path.join(OUT_DIR, f"{name}_confusion_matrix.txt"),
        confusion_matrix(y_test, preds),
        fmt='%d'
    )
    
    # --- ROC ---
    fpr, tpr, _ = roc_curve(y_test, probs_test)
    plt.figure()
    plt.plot(fpr, tpr)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{name} ROC Curve")
    plt.savefig(os.path.join(OUT_DIR, f"roc_{name}.png"))
    plt.close()
    
    # --- PR ---
    prec_curve, rec_curve, _ = precision_recall_curve(y_test, probs_test)
    plt.figure()
    plt.plot(rec_curve, prec_curve)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{name} Precision-Recall Curve")
    plt.savefig(os.path.join(OUT_DIR, f"pr_{name}.png"))
    plt.close()
    
    # --- Calibration ---
    prob_true, prob_pred = calibration_curve(y_test, probs_test, n_bins=10)
    plt.figure()
    plt.plot(prob_pred, prob_true, marker='o')
    plt.plot([0, 1], [0, 1], '--')
    plt.xlabel("Predicted Probability")
    plt.ylabel("Observed Frequency")
    plt.title(f"{name} Calibration Curve")
    plt.savefig(os.path.join(OUT_DIR, f"calibration_{name}.png"))
    plt.close()
    
    # --- Results ---
    result_row = [
        name, roc_auc, auprc, f1, precision, recall,
        mcc, brier, ece, ci_low, ci_high
    ]
    
    with open(results_json_path, 'w') as f:
        json.dump(result_row, f, indent=2)
    
    checkpoint[f"{name}_evaluated"] = True
    save_checkpoint(checkpoint)
    
    results.append(result_row)

# =========================
# DUMMY BASELINE
# =========================

dummy_path = os.path.join(OUT_DIR, "Dummy_best.pkl")
if checkpoint.get("dummy_done") and os.path.exists(dummy_path):
    print("\nDummy already trained. Loading...")
    dummy = joblib.load(dummy_path)
else:
    print("\nTraining Dummy classifier...")
    dummy = DummyClassifier(strategy="stratified")
    dummy.fit(X_train_s, y_train)
    joblib.dump(dummy, dummy_path)
    checkpoint["dummy_done"] = True
    save_checkpoint(checkpoint)

probs_dummy = dummy.predict_proba(X_test_s)[:, 1]
baseline_auc = roc_auc_score(y_test, probs_dummy)

# =========================
# FINAL RESULTS TABLE
# =========================

if results:
    results_df = pd.DataFrame(results, columns=["Model", "AUC", "AUPRC", "F1", "Precision", "Recall",
                                               "MCC", "Brier", "ECE", "AUC_CI_low", "AUC_CI_high"])
    
    results_df.to_csv(os.path.join(OUT_DIR, "baseline_comparison_final.csv"), index=False)
    checkpoint["results_table_done"] = True
    save_checkpoint(checkpoint)

# =========================
# SAVE MODELS (ensure all present)
# =========================

for name, model in models.items():
    model_path = os.path.join(OUT_DIR, f"{name}_best.pkl")
    if not os.path.exists(model_path):
        joblib.dump(model, model_path)

print("\nBASELINE PIPELINE COMPLETE")
