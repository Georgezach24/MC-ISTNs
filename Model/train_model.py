"""
Model/train_model.py

Baseline ML proof-of-concept: predict ServingType (Terrestrial/Satellite/
DualConnectivity) from per-user candidate-link features (best-BS and
satellite SNR/distance/elevation). This is a pipeline sanity check, not the
final Part 2 deliverable - the label is close to a threshold rule on
(CandBS_SNR_dB, CandSat_SNR_dB), so near-perfect accuracy is expected.

Usage:
    python train_model.py
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
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ml_common import CLASS_LABELS, evaluate_model, plot_roc_ovr

ROOT = Path(__file__).resolve().parent.parent
DATASET_PATH = ROOT / "Dataset" / "dataset.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "results"

# Ίδιο κατώφλι με satParameters.MinElevationDeg στο monteCarloDriver.m.
MIN_ELEVATION_DEG = 10.0
SENTINEL_SNR_DB = -50.0
SENTINEL_PATHLOSS_DB = 300.0

FEATURE_COLUMNS_NUMERIC = [
    # BsLoad/SatLoad εξαιρούνται - downstream συνέπεια του ServingType, όχι αιτία.
    "NumBS", "NumUsers",
    "CandBS_SNR_dB", "CandBS_Distance_m", "CandBS_PathLoss_dB",
    "CandSat_SNR_dB", "CandSat_Elevation_deg", "CandSat_SlantRange_m",
    "CandSat_PathLoss_dB",
]
FEATURE_COLUMNS_CATEGORICAL = ["ScenarioType"]
FEATURE_COLUMNS_BOOL = ["CandSat_Visible"]
TARGET_COLUMN = "ServingType"


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Δορυφόρος: -Inf SNR/Inf path loss όταν elevation < MinElevationDeg -> sentinel + flag.
    df["CandSat_Visible"] = df["CandSat_Elevation_deg"] >= MIN_ELEVATION_DEG
    df["CandSat_SNR_dB"] = df["CandSat_SNR_dB"].replace([np.inf, -np.inf], SENTINEL_SNR_DB)
    df["CandSat_PathLoss_dB"] = df["CandSat_PathLoss_dB"].replace([np.inf, -np.inf], SENTINEL_PATHLOSS_DB)

    # Outage (καμία κάλυψη) εξαιρείται - διαφορετικό ερώτημα από ποια ζεύξη χρησιμοποιείται.
    numOutage = int((df["ServingType"] == "Outage").sum())
    if numOutage:
        print(f"Excluding {numOutage} Outage rows (no candidate above minimum usable SNR) "
              f"out of {len(df)} - {'/'.join(CLASS_LABELS)} target only.")
        df = df[df["ServingType"] != "Outage"].reset_index(drop=True)
    return df


def build_preprocessor() -> ColumnTransformer:
    return ColumnTransformer([
        ("num", StandardScaler(), FEATURE_COLUMNS_NUMERIC),
        ("cat", OneHotEncoder(drop="if_binary"),
         FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL),
    ])


def group_train_test_split(df: pd.DataFrame, test_size=0.25, seed=42):
    # Split ανά ScenarioID, όχι ανά γραμμή - αποφυγή διαρροής γεωμετρίας σεναρίου.
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
    train_idx, test_idx = next(splitter.split(df, groups=df["ScenarioID"]))
    return df.iloc[train_idx].reset_index(drop=True), df.iloc[test_idx].reset_index(drop=True)


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
    balance = " / ".join(f"{(df[TARGET_COLUMN] == c).mean():.1%} {c}" for c in CLASS_LABELS)
    print(f"Class balance (all data): {balance}")

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
        y_proba, class_order = evaluate_model(name, pipeline, X_test, y_test, results, RESULTS_DIR)
        roc_curves[name] = (y_test, y_proba, class_order)

        if name == "RandomForest":
            ohe = pipeline.named_steps["preprocess"].named_transformers_["cat"]
            cat_names = list(ohe.get_feature_names_out(FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL))
            all_feature_names = FEATURE_COLUMNS_NUMERIC + cat_names
            importances = pipeline.named_steps["clf"].feature_importances_
            order = np.argsort(importances)[::-1]

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.barh([all_feature_names[i] for i in order][::-1], importances[order][::-1])
            ax.set_title("RandomForest Feature Importance")
            ax.set_xlabel("Importance")
            fig.tight_layout()
            fig.savefig(RESULTS_DIR / "feature_importance_RandomForest.png", dpi=150)
            plt.close(fig)

    plot_roc_ovr(roc_curves, RESULTS_DIR)

    with open(RESULTS_DIR / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved metrics + plots to {RESULTS_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
