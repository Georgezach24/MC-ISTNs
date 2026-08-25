"""
Model/train_model_geometry_only.py

Second, harder ML pass on top of the Part 1 simulation. train_model.py
gives the classifier ground-truth CandBS_SNR_dB/CandSat_SNR_dB, which
already near-determine the label (ServingType = argmax of the two) - a
pipeline sanity check, not a realistic prediction problem.

This script removes that shortcut: it drops both candidate SNR columns
AND both candidate path-loss columns (PathLoss is a near-affine proxy for
SNR given the fixed per-type EIRP/noise floor in simulateScenario.m, so
keeping it would let the model reconstruct SNR anyway and reintroduce the
same shortcut under a different name). What remains is only what a real
system would know about a link *before* measuring it: geometry
(CandBS_Distance_m, CandSat_Elevation_deg, CandSat_SlantRange_m,
CandSat_Visible) and scenario context (NumBS, NumUsers, ScenarioType).

Because the underlying channel is stochastic (per-link LOS/NLOS draw +
log-normal shadow fading, TR 38.901 SS7.4), geometry alone does not
determine the winner - a lower accuracy here than in train_model.py is
the expected, correct outcome, not a regression. The point is realism,
not the metric.

Usage:
    python train_model_geometry_only.py
"""

import json
import sys
from pathlib import Path

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    RocCurveDisplay,
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

ROOT = Path(__file__).resolve().parent.parent
DATASET_PATH = ROOT / "Dataset" / "dataset.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "results_geometry_only"

# Ίδιο κατώφλι με satParameters.MinElevationDeg στο monteCarloDriver.m
# (TR 38.821 visibility mask) - όχι μια νέα υπόθεση, απλά επαναχρησιμοποίηση.
MIN_ELEVATION_DEG = 10.0

# CandBS_SNR_dB/CandSat_SNR_dB και CandBS_PathLoss_dB/CandSat_PathLoss_dB
# αποκλείονται σκόπιμα (βλ. docstring): δίνουν στο μοντέλο την απάντηση, ή
# ένα σχεδόν-affine ισοδύναμό της. Μένουν μόνο γεωμετρικά/context
# χαρακτηριστικά, διαθέσιμα σε ένα πραγματικό σύστημα πριν τη μέτρηση SNR.
FEATURE_COLUMNS_NUMERIC = [
    # NodeLoad exclude σκόπιμα: είναι συνέπεια του ServingType, όχι
    # ανεξάρτητος predictor (βλ. train_model.py).
    "NumBS", "NumUsers",
    "CandBS_Distance_m",
    "CandSat_Elevation_deg", "CandSat_SlantRange_m",
]
FEATURE_COLUMNS_CATEGORICAL = ["ScenarioType"]
FEATURE_COLUMNS_BOOL = ["CandSat_Visible"]
TARGET_COLUMN = "ServingType"


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # simulateScenario.m πλέον καταγράφει και ServingType="Outage" (κανένας
    # υποψήφιος δεν ξεπερνά το ελάχιστο χρησιμοποιήσιμο SNR) - εξαιρείται
    # εδώ, ίδια λογική με το train_model.py.
    numOutage = int((df["ServingType"] == "Outage").sum())
    if numOutage:
        print(f"Excluding {numOutage} Outage rows (no candidate above minimum usable SNR) "
              f"out of {len(df)} - binary Terrestrial/Satellite target only.")
        df = df[df["ServingType"] != "Outage"].reset_index(drop=True)

    # CandSat_Elevation_deg/CandSat_SlantRange_m είναι πάντα πεπερασμένα
    # (γεωμετρία, όχι SNR/path loss) - το μόνο που χρειάζεται είναι η
    # boolean σημαία ορατότητας.
    df["CandSat_Visible"] = df["CandSat_Elevation_deg"] >= MIN_ELEVATION_DEG
    return df


def build_preprocessor() -> ColumnTransformer:
    return ColumnTransformer([
        ("num", StandardScaler(), FEATURE_COLUMNS_NUMERIC),
        ("cat", OneHotEncoder(drop="if_binary"),
         FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL),
    ])


def group_train_test_split(df: pd.DataFrame, test_size=0.25, seed=42):
    # Split ανά ScenarioID (όχι ανά γραμμή): χρήστες του ίδιου σεναρίου
    # μοιράζονται τις ίδιες θέσεις BS/δορυφόρου, άρα ένα row-level split θα
    # διέρρεε γεωμετρία σεναρίου ανάμεσα σε train/test.
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
    train_idx, test_idx = next(splitter.split(df, groups=df["ScenarioID"]))
    return df.iloc[train_idx].reset_index(drop=True), df.iloc[test_idx].reset_index(drop=True)


def evaluate_model(name, pipeline, X_test, y_test, results):
    y_pred = pipeline.predict(X_test)
    y_proba = pipeline.predict_proba(X_test)[:, list(pipeline.classes_).index("Satellite")]

    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, pos_label="Satellite")
    auc = roc_auc_score((y_test == "Satellite").astype(int), y_proba)
    report = classification_report(y_test, y_pred, output_dict=True)

    print(f"\n=== {name} ===")
    print(f"Accuracy: {acc:.4f}  F1(Satellite): {f1:.4f}  ROC-AUC: {auc:.4f}")
    print(classification_report(y_test, y_pred))

    results[name] = {
        "accuracy": acc,
        "f1_satellite": f1,
        "roc_auc": auc,
        "classification_report": report,
    }

    cm = confusion_matrix(y_test, y_pred, labels=["Terrestrial", "Satellite"])
    disp = ConfusionMatrixDisplay(cm, display_labels=["Terrestrial", "Satellite"])
    fig, ax = plt.subplots(figsize=(4, 4))
    disp.plot(ax=ax, cmap="Blues", colorbar=False)
    ax.set_title(f"{name} - Confusion Matrix (geometry-only)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"confusion_matrix_{name}.png", dpi=150)
    plt.close(fig)

    return y_proba


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not DATASET_PATH.exists():
        raise FileNotFoundError(
            f"{DATASET_PATH} not found - run monteCarloDriver.m in MATLAB first "
            "to generate the dataset (see PROD/monteCarloDriver.m)."
        )

    df = load_dataset(DATASET_PATH)
    train_df, test_df = group_train_test_split(df)

    print(f"Loaded {len(df)} user-rows from {df['ScenarioID'].nunique()} scenarios")
    print(f"Train: {len(train_df)} rows ({train_df['ScenarioID'].nunique()} scenarios)")
    print(f"Test:  {len(test_df)} rows ({test_df['ScenarioID'].nunique()} scenarios)")
    print(f"Class balance (all data): "
          f"{(df[TARGET_COLUMN] == 'Terrestrial').mean():.1%} Terrestrial / "
          f"{(df[TARGET_COLUMN] == 'Satellite').mean():.1%} Satellite")

    feature_cols = FEATURE_COLUMNS_NUMERIC + FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL
    X_train, y_train = train_df[feature_cols], train_df[TARGET_COLUMN]
    X_test, y_test = test_df[feature_cols], test_df[TARGET_COLUMN]

    models = {
        "LogisticRegression": LogisticRegression(max_iter=1000),
        "RandomForest": RandomForestClassifier(n_estimators=300, max_depth=12, random_state=42),
    }

    results = {}
    roc_curves = {}
    for name, clf in models.items():
        pipeline = Pipeline([
            ("preprocess", build_preprocessor()),
            ("clf", clf),
        ])
        pipeline.fit(X_train, y_train)
        y_proba = evaluate_model(name, pipeline, X_test, y_test, results)
        roc_curves[name] = (y_test, y_proba)

        if name == "RandomForest":
            ohe = pipeline.named_steps["preprocess"].named_transformers_["cat"]
            cat_names = list(ohe.get_feature_names_out(FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL))
            all_feature_names = FEATURE_COLUMNS_NUMERIC + cat_names
            importances = pipeline.named_steps["clf"].feature_importances_
            order = np.argsort(importances)[::-1]

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.barh([all_feature_names[i] for i in order][::-1], importances[order][::-1])
            ax.set_title("RandomForest Feature Importance (geometry-only)")
            ax.set_xlabel("Importance")
            fig.tight_layout()
            fig.savefig(RESULTS_DIR / "feature_importance_RandomForest.png", dpi=150)
            plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    for name, (y_true, y_proba) in roc_curves.items():
        RocCurveDisplay.from_predictions((y_true == "Satellite").astype(int), y_proba, name=name, ax=ax)
    ax.set_title("ROC Curve - Predicting Satellite vs Terrestrial (geometry-only)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "roc_curve.png", dpi=150)
    plt.close(fig)

    with open(RESULTS_DIR / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved metrics + plots to {RESULTS_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
