"""
Model/train_model_noisy_snr.py

Third ML pass, between train_model.py (exact SNR) and
train_model_geometry_only.py (no SNR): keeps candidate SNRs as features but
adds Gaussian noise before the model sees them, standing in for a stale/
imperfect measurement report (TS 38.331 Event A3/A5). PathLoss columns
still dropped (near-affine SNR proxy).

Noise sigma: per-ScenarioType TR 38.901 v16.1.0 Table 7.4.1-1 NLOS
shadow-fading sigma (UMa=6dB, UMi=7.82dB - NLOS since borderline placement
distances are NLOS-dominated and the dataset has no per-row LOS flag).
Satellite sigma still read from Results/kpi_link_snr_std.csv (~0.07dB,
deterministic channel, no fading model on that side).

Usage:
    python train_model_noisy_snr.py
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
KPI_SUMMARY_PATH = ROOT / "Results" / "kpi_link_snr_std.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "results_noisy_snr"

MIN_ELEVATION_DEG = 10.0
SENTINEL_SNR_DB = -50.0
NOISE_SEED = 42

# 3GPP TR 38.901 v16.1.0, Πίνακας 7.4.1-1, σ_SF NLOS ανά ScenarioType (dB).
BS_SIGMA_SF_NLOS_DB = {"UMa": 6.0, "UMi": 7.82}

FEATURE_COLUMNS_NUMERIC = [
    "NumBS", "NumUsers",
    "CandBS_Distance_m",
    "CandSat_Elevation_deg", "CandSat_SlantRange_m",
    "CandBS_SNR_noisy_dB", "CandSat_SNR_noisy_dB",
]
FEATURE_COLUMNS_CATEGORICAL = ["ScenarioType"]
FEATURE_COLUMNS_BOOL = ["CandSat_Visible"]
TARGET_COLUMN = "ServingType"


def load_sat_noise_sigma(path: Path) -> float:
    kpi = pd.read_csv(path).set_index("LinkType")
    return kpi.loc["Satellite", "std_SNR_dB"]


def load_dataset(path: Path, sigma_bs_map: dict, sigma_sat: float) -> pd.DataFrame:
    df = pd.read_csv(path)

    numOutage = int((df["ServingType"] == "Outage").sum())
    if numOutage:
        print(f"Excluding {numOutage} Outage rows (no candidate above minimum usable SNR) "
              f"out of {len(df)} - {'/'.join(CLASS_LABELS)} target only.")
        df = df[df["ServingType"] != "Outage"].reset_index(drop=True)

    df["CandSat_Visible"] = df["CandSat_Elevation_deg"] >= MIN_ELEVATION_DEG

    # Sentinel πριν τον θόρυβο, ώστε -Inf να μη γίνει πεπερασμένος αλλά παραπλανητικός αριθμός.
    df["CandSat_SNR_dB"] = df["CandSat_SNR_dB"].replace([np.inf, -np.inf], SENTINEL_SNR_DB)

    rng = np.random.default_rng(NOISE_SEED)
    sigma_bs_per_row = df["ScenarioType"].map(sigma_bs_map).to_numpy(dtype=float)
    df["CandBS_SNR_noisy_dB"] = df["CandBS_SNR_dB"] + rng.normal(0.0, sigma_bs_per_row, size=len(df))
    df["CandSat_SNR_noisy_dB"] = df["CandSat_SNR_dB"] + rng.normal(0.0, sigma_sat, size=len(df))
    return df


def build_preprocessor() -> ColumnTransformer:
    return ColumnTransformer([
        ("num", StandardScaler(), FEATURE_COLUMNS_NUMERIC),
        ("cat", OneHotEncoder(drop="if_binary"),
         FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL),
    ])


def group_train_test_split(df: pd.DataFrame, test_size=0.25, seed=42):
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
    if not KPI_SUMMARY_PATH.exists():
        raise FileNotFoundError(
            f"{KPI_SUMMARY_PATH} not found - run kpiRepeatedRuns.m in MATLAB first "
            "to generate the per-link noise-sigma source (see PROD/kpiRepeatedRuns.m)."
        )

    sigma_sat = load_sat_noise_sigma(KPI_SUMMARY_PATH)
    sigmas = {"Terrestrial_by_ScenarioType": BS_SIGMA_SF_NLOS_DB, "Satellite": sigma_sat}
    print(f"Noise sigma: Terrestrial(UMa)={BS_SIGMA_SF_NLOS_DB['UMa']:.2f} dB "
          f"(TR 38.901 Table 7.4.1-1 NLOS), Terrestrial(UMi)={BS_SIGMA_SF_NLOS_DB['UMi']:.2f} dB, "
          f"Satellite={sigma_sat:.4f} dB (from {KPI_SUMMARY_PATH.relative_to(ROOT)})")

    df = load_dataset(DATASET_PATH, BS_SIGMA_SF_NLOS_DB, sigma_sat)
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

    results = {"noise_sigma_dB": sigmas}
    roc_curves = {}
    for name, clf in models.items():
        pipeline = Pipeline([
            ("preprocess", build_preprocessor()),
            ("clf", clf),
        ])
        pipeline.fit(X_train, y_train)
        y_proba, class_order = evaluate_model(name, pipeline, X_test, y_test, results, RESULTS_DIR, " (noisy SNR)")
        roc_curves[name] = (y_test, y_proba, class_order)

        if name == "RandomForest":
            ohe = pipeline.named_steps["preprocess"].named_transformers_["cat"]
            cat_names = list(ohe.get_feature_names_out(FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL))
            all_feature_names = FEATURE_COLUMNS_NUMERIC + cat_names
            importances = pipeline.named_steps["clf"].feature_importances_
            order = np.argsort(importances)[::-1]

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.barh([all_feature_names[i] for i in order][::-1], importances[order][::-1])
            ax.set_title("RandomForest Feature Importance (noisy SNR)")
            ax.set_xlabel("Importance")
            fig.tight_layout()
            fig.savefig(RESULTS_DIR / "feature_importance_RandomForest.png", dpi=150)
            plt.close(fig)

    plot_roc_ovr(roc_curves, RESULTS_DIR, " (noisy SNR)")

    with open(RESULTS_DIR / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved metrics + plots to {RESULTS_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
