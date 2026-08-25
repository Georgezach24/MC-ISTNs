"""
Model/train_model_noisy_snr.py

Third ML pass, sitting between train_model.py (exact, instantaneous
ground-truth SNR - unrealistic, a real system never has this for a link
it hasn't connected to) and train_model_geometry_only.py (no SNR at all -
more pessimistic than a real system too, since real UEs do get *some*
signal-quality estimate for candidate cells via periodic neighbor-cell
measurement reports, e.g. 3GPP TS 38.331 Event A3/A5).

Here the candidate SNRs are kept as features, but each one has Gaussian
noise added before the model ever sees it, standing in for "a slightly
stale / imperfect measurement report" rather than the exact instantaneous
value. The noise magnitude is not an invented fudge factor: it is the
*actual measured* run-to-run SNR standard deviation from 500 repeated
stochastic realizations of the same static topology
(kpiRepeatedRuns.m -> Results/kpi_link_snr_std.csv), i.e. empirically
"how differently would a second, independent measurement of this same
link read" under this simulation's own channel model (TR 38.901 LOS/NLOS
draw + shadow fading for the terrestrial side). CandBS_PathLoss_dB /
CandSat_PathLoss_dB are still dropped (near-affine proxies for the exact
SNR, would reintroduce the shortcut train_model_geometry_only.py removes).

The sigma source is Results/kpi_link_snr_std.csv (std of CandBS_SNR_dB /
CandSat_SNR_dB over usable candidates, i.e. per-LINK variance), not the
older Results/kpi_summary_by_type.csv (std of the WINNER's SNR, grouped by
ServingType): after DualConnectivity was added to simulateScenario.m, the
static 6-user topology no longer produces any "Terrestrial"-only winning
row at all (near users always qualify for DualConnectivity - see
CLAUDE.md, dual-connectivity), so a group-by-ServingType summary has no
"Terrestrial" row to read a sigma from any more. The satellite sigma is
small because this simulation's satellite channel is deterministic given
geometry (no fading model on that side yet - see CLAUDE.md Standards
section); it is used as-is rather than inflated, since inventing a larger
number would not be grounded in anything the project has actually modeled
or measured.

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
SENTINEL_SNR_DB = -50.0  # πρακτικά άχρηστος, ίδιο sentinel με το train_model.py
NOISE_SEED = 42

FEATURE_COLUMNS_NUMERIC = [
    # BsLoad/SatLoad exclude σκόπιμα (βλ. train_model.py) - κυκλικός predictor.
    "NumBS", "NumUsers",
    "CandBS_Distance_m",
    "CandSat_Elevation_deg", "CandSat_SlantRange_m",
    # "Θορυβώδεις" εκδοχές του SNR αντί για το ακριβές (βλ. docstring) -
    # τα CandBS_PathLoss_dB/CandSat_PathLoss_dB παραμένουν εκτός, όπως στο
    # train_model_geometry_only.py.
    "CandBS_SNR_noisy_dB", "CandSat_SNR_noisy_dB",
]
FEATURE_COLUMNS_CATEGORICAL = ["ScenarioType"]
FEATURE_COLUMNS_BOOL = ["CandSat_Visible"]
TARGET_COLUMN = "ServingType"


def load_noise_sigmas(path: Path) -> dict:
    kpi = pd.read_csv(path).set_index("LinkType")
    return {
        "Terrestrial": kpi.loc["Terrestrial", "std_SNR_dB"],
        "Satellite": kpi.loc["Satellite", "std_SNR_dB"],
    }


def load_dataset(path: Path, sigma_bs: float, sigma_sat: float) -> pd.DataFrame:
    df = pd.read_csv(path)

    # simulateScenario.m πλέον καταγράφει και ServingType="Outage" (κανένας
    # υποψήφιος δεν ξεπερνά το ελάχιστο χρησιμοποιήσιμο SNR) - εξαιρείται
    # εδώ, ίδια λογική με το train_model.py.
    numOutage = int((df["ServingType"] == "Outage").sum())
    if numOutage:
        print(f"Excluding {numOutage} Outage rows (no candidate above minimum usable SNR) "
              f"out of {len(df)} - {'/'.join(CLASS_LABELS)} target only.")
        df = df[df["ServingType"] != "Outage"].reset_index(drop=True)

    df["CandSat_Visible"] = df["CandSat_Elevation_deg"] >= MIN_ELEVATION_DEG

    # Ο δορυφόρος έχει -Inf SNR όταν elevation < MinElevationDeg (visibility
    # mask). Sentinel πριν προστεθεί θόρυβος, ώστε ο θόρυβος να μην
    # μετατρέψει ένα -Inf σε έναν πεπερασμένο, παραπλανητικό αριθμό.
    df["CandSat_SNR_dB"] = df["CandSat_SNR_dB"].replace([np.inf, -np.inf], SENTINEL_SNR_DB)

    rng = np.random.default_rng(NOISE_SEED)
    df["CandBS_SNR_noisy_dB"] = df["CandBS_SNR_dB"] + rng.normal(0.0, sigma_bs, size=len(df))
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

    sigmas = load_noise_sigmas(KPI_SUMMARY_PATH)
    print(f"Noise sigma (from {KPI_SUMMARY_PATH.relative_to(ROOT)}): "
          f"Terrestrial={sigmas['Terrestrial']:.4f} dB, Satellite={sigmas['Satellite']:.4f} dB")

    df = load_dataset(DATASET_PATH, sigmas["Terrestrial"], sigmas["Satellite"])
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
