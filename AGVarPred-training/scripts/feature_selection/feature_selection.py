import pandas as pd
import numpy as np
import gc, re, os
from glob import glob
import joblib

from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import RobustScaler

from lightgbm import LGBMClassifier
import matplotlib.pyplot as plt

# =========================
# OUTPUT DIR
# =========================

OUTDIR_FS = "feature_selection_output_nested"
os.makedirs(OUTDIR_FS, exist_ok=True)

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
# PREPROCESSING HELPER
# =========================

def preprocess_fold(X_tr_raw, X_val_raw, num_cols, cat_cols):
    """
    Impute + scale numerics; impute categoricals with 'missing'.
    Returns processed train/val DataFrames and fitted scaler.
    """
    # --- NUMERIC: median impute + RobustScaler ---
    median_values = X_tr_raw[num_cols].median(numeric_only=True)
    X_tr_num = X_tr_raw[num_cols].fillna(median_values)
    X_val_num = X_val_raw[num_cols].fillna(median_values)

    scaler = RobustScaler()
    X_tr_num_s = pd.DataFrame(
        scaler.fit_transform(X_tr_num),
        columns=num_cols,
        index=X_tr_num.index
    )
    X_val_num_s = pd.DataFrame(
        scaler.transform(X_val_num),
        columns=num_cols,
        index=X_val_num.index
    )

    # --- CATEGORICAL: fill NaN + align categories ---
    X_tr_cat = X_tr_raw[cat_cols].fillna("missing")
    X_val_cat = X_val_raw[cat_cols].fillna("missing")

    for col in cat_cols:
        train_cats = X_tr_cat[col].unique().tolist()
        val_cats = X_val_cat[col].dropna().unique().tolist()
        combined = train_cats + [v for v in val_cats if v not in train_cats]
        cat_dtype = pd.CategoricalDtype(categories=combined)
        X_tr_cat[col] = X_tr_cat[col].astype(cat_dtype)
        X_val_cat[col] = X_val_cat[col].astype(cat_dtype)

    # --- RECOMBINE ---
    X_tr_proc = pd.concat([X_tr_num_s, X_tr_cat], axis=1)
    X_val_proc = pd.concat([X_val_num_s, X_val_cat], axis=1)

    return X_tr_proc, X_val_proc, scaler


def preprocess_global(X_raw, num_cols, cat_cols):
    """
    Global imputation + scaling (fit on full train set).
    Returns processed DataFrame, scaler, and median values.
    """
    global_medians = X_raw[num_cols].median(numeric_only=True)
    X_num = X_raw[num_cols].fillna(global_medians)

    scaler = RobustScaler()
    X_num_s = pd.DataFrame(
        scaler.fit_transform(X_num),
        columns=num_cols,
        index=X_num.index
    )

    X_cat = X_raw[cat_cols].fillna("missing")
    for col in cat_cols:
        train_cats = pd.CategoricalDtype(categories=X_cat[col].unique())
        X_cat[col] = X_cat[col].astype(train_cats)

    X_proc = pd.concat([X_num_s, X_cat], axis=1)
    return X_proc, scaler, global_medians


# =========================
# LOAD TRAIN DATA ONLY
# =========================

print("🔄 Loading TRAIN dataset...")

files = glob("final_dataset_parts_train/*.parquet")

dfs = []
for f in files:
    df = pd.read_parquet(f)
    df.columns = [clean_name(c) for c in df.columns]
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)
del dfs; gc.collect()

print("Dataset:", df.shape)

# =========================
# PREPARE DATA
# =========================

groups = df["gene"].values

drop_cols = [
    'variant_id',
    'label',
    'gene',
    'GeneSymbol'
]

drop_cols = [c for c in drop_cols if c in df.columns]

X = df.drop(columns=drop_cols)
y = df['label']

X = X.loc[:, ~X.columns.duplicated()]

# --- SEPARATE NUMERIC AND CATEGORICAL ---
cat_cols = X.select_dtypes(include=['category', 'object', 'str']).columns.tolist()
num_cols = [c for c in X.columns if c not in cat_cols]
print(f"Numeric columns: {len(num_cols)}, Categorical columns: {len(cat_cols)}")
if cat_cols:
    print(f"Categoricals: {cat_cols}")

SCALE_POS_WEIGHT = (y == 0).sum() / (y == 1).sum()

# =========================
# CONFIG
# =========================

feature_steps = [40, 60, 80, 100, 120, 140, 160, 200, 250, 300, 400, 500, 750, 1000, 1500, 2000]

kf = GroupKFold(n_splits=5)

BASE_PARAMS = dict(
    n_estimators=300,
    n_jobs=-1,
    random_state=42,
    verbose=-1,
    importance_type='gain',
    scale_pos_weight=SCALE_POS_WEIGHT
)

# =========================
# NESTED FEATURE SELECTION
# =========================

mean_auc, std_auc, mean_auprc, std_auprc = [], [], [], []

print("\n🔄 Running nested feature selection...")

for n in feature_steps:

    fold_auc = []
    fold_auprc = []

    for fold, (train_idx, val_idx) in enumerate(kf.split(X, y, groups)):

        # Split raw data
        X_tr_raw, X_val_raw = X.iloc[train_idx], X.iloc[val_idx]
        y_tr, y_val = y.iloc[train_idx], y.iloc[val_idx]

        # Preprocess (no leakage)
        X_tr, X_val, _ = preprocess_fold(X_tr_raw, X_val_raw, num_cols, cat_cols)

        # Categorical feature names for this fold
        fold_cat_cols = [c for c in cat_cols if c in X_tr.columns]

        # Train importance model on ALL features
        model_imp = LGBMClassifier(**BASE_PARAMS)
        model_imp.fit(X_tr, y_tr, categorical_feature=fold_cat_cols)

        importance = pd.Series(
            model_imp.feature_importances_,
            index=X_tr.columns
        ).sort_values(ascending=False)

        top_features = importance.index[:n]

        # Train reduced model on top features
        top_cat_cols = [c for c in fold_cat_cols if c in top_features]
        model = LGBMClassifier(**BASE_PARAMS)
        model.fit(X_tr[top_features], y_tr, categorical_feature=top_cat_cols)

        probs = model.predict_proba(X_val[top_features])[:, 1]

        fold_auc.append(roc_auc_score(y_val, probs))
        fold_auprc.append(average_precision_score(y_val, probs))

        del X_tr_raw, X_val_raw, X_tr, X_val, model, model_imp
        gc.collect()

    mean_auc.append(np.mean(fold_auc))
    std_auc.append(np.std(fold_auc))
    mean_auprc.append(np.mean(fold_auprc))
    std_auprc.append(np.std(fold_auprc))

    print(f"{n} features → AUC: {mean_auc[-1]:.4f} | AUPRC: {mean_auprc[-1]:.4f}")

# =========================
# 1-STD RULE
# =========================

best_idx = np.argmax(mean_auprc)
best_auprc_val = mean_auprc[best_idx]
best_std = std_auprc[best_idx]

threshold = best_auprc_val - best_std

candidates = [n for n, auprc in zip(feature_steps, mean_auprc) if auprc >= threshold]

optimal_n = min(candidates)

print("\n🏆 Optimal number of features (by AUPRC):", optimal_n)

# =========================
# GLOBAL IMPUTATION + SCALER
# =========================

print("🔄 Fitting global imputer + RobustScaler on full training set...")

X_scaled, scaler_global, global_medians = preprocess_global(X, num_cols, cat_cols)

joblib.dump(scaler_global, f"{OUTDIR_FS}/robust_scaler.pkl")
joblib.dump(global_medians, f"{OUTDIR_FS}/global_medians.pkl")

# =========================
# FINAL FEATURE SET
# =========================

print("🔄 Computing final global importance...")

final_model = LGBMClassifier(**BASE_PARAMS)
final_model.fit(X_scaled, y, categorical_feature=cat_cols)

# SAVE MODEL
joblib.dump(final_model, f"{OUTDIR_FS}/feature_selection_model.pkl")

final_importance = pd.Series(
    final_model.feature_importances_,
    index=X_scaled.columns
).sort_values(ascending=False)

selected_features = final_importance.index[:optimal_n]

# =========================
# SAVE OUTPUTS
# =========================

importance_df = pd.DataFrame({
    "feature": final_importance.index,
    "importance": final_importance.values
})
importance_df.to_csv(f"{OUTDIR_FS}/feature_importance_full.csv", index=False)

selected_df = importance_df.head(optimal_n)
selected_df.to_csv(f"{OUTDIR_FS}/selected_features_with_importance.csv", index=False)

with open(f"{OUTDIR_FS}/selected_features.txt", "w") as f:
    for feat in selected_features:
        f.write(feat + "\n")

pd.DataFrame({
    "features": feature_steps,
    "auc": mean_auc,
    "std_auc": std_auc,
    "auprc": mean_auprc,
    "std_auprc": std_auprc
}).to_csv(f"{OUTDIR_FS}/feature_selection_results.csv", index=False)

# =========================
# BASELINE COMPARISON
# =========================

print("\n🔄 Running baseline comparison...")

baseline_auc = []

for train_idx, val_idx in kf.split(X, y, groups):

    X_tr_raw, X_val_raw = X.iloc[train_idx], X.iloc[val_idx]

    X_tr, X_val, _ = preprocess_fold(X_tr_raw, X_val_raw, num_cols, cat_cols)

    fold_cat_cols = [c for c in cat_cols if c in X_tr.columns]

    model = LGBMClassifier(**BASE_PARAMS)
    model.fit(
        X_tr,
        y.iloc[train_idx],
        categorical_feature=fold_cat_cols
    )

    probs = model.predict_proba(X_val)[:, 1]

    baseline_auc.append(
        roc_auc_score(
            y.iloc[val_idx],
            probs
        )
    )

    del X_tr_raw, X_val_raw, X_tr, X_val, model
    gc.collect()

selected_auc = mean_auc[feature_steps.index(optimal_n)]
selected_auprc = mean_auprc[feature_steps.index(optimal_n)]

plt.figure()
plt.bar(['Baseline', 'Selected'], [np.mean(baseline_auc), selected_auc])
plt.ylabel("AUC")
plt.title(f"AUC Comparison (optimal n={optimal_n}, AUPRC={selected_auprc:.4f})")
plt.savefig(f"{OUTDIR_FS}/baseline_comparison.png", dpi=200)
plt.close()

# =========================
# PLOTS
# =========================

plt.figure()
plt.plot(
    feature_steps,
    mean_auc,
    marker='o'
)
plt.axvline(
    x=optimal_n,
    linestyle=':',
    color='red',
    label=f'n={optimal_n}'
)
plt.legend()
plt.xlabel("Features")
plt.ylabel("AUC")

plt.savefig(
    f"{OUTDIR_FS}/feature_selection_curve_auc.png",
    dpi=300
)
plt.close()

plt.figure()
plt.plot(feature_steps, mean_auprc, marker='o')
plt.axvline(x=optimal_n, linestyle=':', color='red')
plt.xlabel("Features")
plt.ylabel("AUPRC")
plt.savefig(f"{OUTDIR_FS}/feature_selection_curve_auprc.png", dpi=300)
plt.close()

top = final_importance.head(30)

plt.figure(figsize=(10,8))
top[::-1].plot(kind='barh')
plt.tight_layout()
plt.savefig(f"{OUTDIR_FS}/top_features.png", dpi=300)
plt.close()

non_zero = final_importance[final_importance > 0]

plt.figure()
plt.hist(non_zero, bins=50)
plt.yscale('log')
plt.savefig(f"{OUTDIR_FS}/importance_distribution.png", dpi=300)
plt.close()

print("\n🎯 DONE — Nested feature selection complete")
