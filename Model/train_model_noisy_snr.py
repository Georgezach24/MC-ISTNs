"""
Model/train_model_noisy_snr.py

Third ML pass, sitting between train_model.py (exact, instantaneous
ground-truth SNR - unrealistic, a real system never has this for a link
it hasn't connected to) and train_model_geometry_only.py (no SNR at all -
more pessimistic than a real system too, since real UEs do get *some*
signal-quality estimate for candidate cells via periodic neighbor-cell
measurement reports, e.g. 3GPP TS 38.331 Event A3/A5).

Here the candidate SINRs are kept as features, but each one has Gaussian
noise added before the model ever sees it, standing in for "a slightly
stale / imperfect measurement report" rather than the exact instantaneous
value. The noise magnitude is not an invented fudge factor: it is the
*actual measured* SNR standard deviation per serving type, computed from
the same unified run that produced the dataset (runSimulation.m), i.e.
empirically "how differently would a second, independent measurement of
this same link read" under this simulation's own channel model (TR 38.901
LOS/NLOS draw + shadow fading terrestrially, shadowed-Rician on the
satellite side). CandBS_PathLoss_dB / CandSat_PathLoss_dB are still
dropped (near-affine proxies for the exact SNR, would reintroduce the
shortcut train_model_geometry_only.py removes).

Both sigmas are printed at run time rather than hard-coded, so they always
match the dataset actually being trained on.

Usage:
    python train_model_noisy_snr.py
"""

import hashlib
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

import joblib
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
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

ROOT = Path(__file__).resolve().parent.parent
DATASET_PATH = ROOT / "Dataset" / "dataset.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "results_noisy_snr"

MIN_ELEVATION_DEG = 20.0
SENTINEL_SINR_DB = -50.0  # πρακτικά άχρηστος, ίδιο sentinel με το train_model.py
NOISE_SEED = 42

# Η υποδειγματοληψία στον χρόνο γίνεται ΣΤΗΝ ΠΗΓΗ: το runSimulation.m
# προχωρά με βήμα 1 s αλλά γράφει μία γραμμή κάθε datasetStride βήματα
# (5 s), γιατί διαδοχικά δείγματα του ίδιου χρήστη απέχουν 0.83 m και είναι
# σχεδόν ταυτόσημα. Εδώ δεν χρειάζεται επιπλέον αραίωση· άλλαξε το μόνο αν
# θέλεις ακόμη πιο αραιό δείγμα.
ML_SAMPLE_STRIDE = 1

FEATURE_COLUMNS_NUMERIC = [
    # NodeLoad exclude σκόπιμα (βλ. train_model.py) - κυκλικός predictor.
    "NumUsers",
    "CandBS_Distance_m",
    "CandSat_Elevation_deg", "CandSat_SlantRange_m",
    # "Θορυβώδεις" εκδοχές του SNR αντί για το ακριβές (βλ. docstring) -
    # τα CandBS_PathLoss_dB/CandSat_PathLoss_dB παραμένουν εκτός, όπως στο
    # train_model_geometry_only.py.
    "CandBS_SINR_noisy_dB", "CandSat_SINR_noisy_dB",
]
FEATURE_COLUMNS_CATEGORICAL = []  # το σενάριο διάδοσης είναι σταθερό (UMa)
FEATURE_COLUMNS_BOOL = ["CandSat_Visible"]
TARGET_COLUMN = "ServingType"

# Παράμετροι διαχωρισμού, ορισμένοι μία φορά ώστε το split.json να μην
# μπορεί να διαφωνήσει με τον διαχωρισμό που έγινε στην πραγματικότητα.
SPLIT_TEST_SIZE = 0.25
SPLIT_SEED = 42


def load_noise_sigmas(df: pd.DataFrame) -> dict:
    """Τυπική απόκλιση του SNR ανά τύπο εξυπηρέτησης, υπολογισμένη από το ίδιο
    το σύνολο δεδομένων. Παλαιότερα διαβαζόταν από το kpi_summary_by_type.csv
    που παρήγαγε χωριστό script· τώρα που η προσομοίωση είναι ενιαία, ο ίδιος
    αριθμός προκύπτει από πολύ μεγαλύτερο δείγμα της ίδιας εκτέλεσης."""
    g = df.groupby("ServingType")["SINR_dB"].std()
    return {
        "Terrestrial": float(g.loc["Terrestrial"]),
        "Satellite": float(g.loc["Satellite"]),
    }


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Υποδειγματοληψία στον χρόνο (βλ. ML_SAMPLE_STRIDE).
    if ML_SAMPLE_STRIDE > 1 and "Step" in df.columns:
        before = len(df)
        df = df[df["Step"] % ML_SAMPLE_STRIDE == 1].reset_index(drop=True)
        print(f"Time decimation: kept {len(df)} of {before} rows "
              f"(1 sample every {ML_SAMPLE_STRIDE} s)")

    # simulateScenario.m πλέον καταγράφει και ServingType="Outage" (κανένας
    # υποψήφιος δεν ξεπερνά το ελάχιστο χρησιμοποιήσιμο SINR) - εξαιρείται
    # εδώ, ίδια λογική με το train_model.py.
    numOutage = int((df["ServingType"] == "Outage").sum())
    if numOutage:
        print(f"Excluding {numOutage} Outage rows (no candidate above minimum usable SINR) "
              f"out of {len(df)} - binary Terrestrial/Satellite target only.")
        df = df[df["ServingType"] != "Outage"].reset_index(drop=True)

    # Ζεύξεις εκτός του πεδίου ισχύος των UMa/UMi (BsUnavailReason =
    # "OutOfModelRange") δεν έχουν υπολογισμένο επίγειο υποψήφιο: τα
    # CandBS_* είναι NaN εξ ορισμού. Πρόκειται για περιορισμό της
    # προσομοίωσης και όχι για φυσική κατάσταση προς πρόβλεψη (βλ. κεφ.
    # μεθοδολογίας, διάκριση αιτίων μη διαθεσιμότητας), οπότε οι γραμμές
    # αυτές εξαιρούνται αντί να τους αποδοθεί τεχνητή τιμή.
    numOutOfRange = int(df["BsUnavailReason"].eq("OutOfModelRange").sum()) if "BsUnavailReason" in df.columns else 0
    if numOutOfRange:
        print(f"Excluding {numOutOfRange} rows with no valid terrestrial candidate "
              f"(outside UMa/UMi validity range) out of {len(df)}.")
        df = df[~df["BsUnavailReason"].eq("OutOfModelRange")].reset_index(drop=True)

    df["CandSat_Visible"] = df["CandSat_Elevation_deg"] >= MIN_ELEVATION_DEG

    # Ο δορυφόρος έχει -Inf SNR όταν elevation < MinElevationDeg (visibility
    # mask). Sentinel πριν προστεθεί θόρυβος, ώστε ο θόρυβος να μην
    # μετατρέψει ένα -Inf σε έναν πεπερασμένο, παραπλανητικό αριθμό.
    df["CandSat_SINR_dB"] = df["CandSat_SINR_dB"].replace([np.inf, -np.inf], SENTINEL_SINR_DB)
    return df


def add_snr_noise(df: pd.DataFrame, sigma_bs: float, sigma_sat: float) -> pd.DataFrame:
    """Προσθέτει γκαουσιανό θόρυβο μέτρησης στα δύο υποψήφια SINR."""
    rng = np.random.default_rng(NOISE_SEED)
    df["CandBS_SINR_noisy_dB"] = df["CandBS_SINR_dB"] + rng.normal(0.0, sigma_bs, size=len(df))
    df["CandSat_SINR_noisy_dB"] = df["CandSat_SINR_dB"] + rng.normal(0.0, sigma_sat, size=len(df))
    return df


def build_preprocessor() -> ColumnTransformer:
    return ColumnTransformer([
        ("num", StandardScaler(), FEATURE_COLUMNS_NUMERIC),
        ("cat", OneHotEncoder(drop="if_binary"),
         FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL),
    ])


def group_train_test_split(df: pd.DataFrame, test_size=SPLIT_TEST_SIZE,
                          seed=SPLIT_SEED):
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
    train_idx, test_idx = next(splitter.split(df, groups=df["PassID"]))
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
    ax.set_title(f"{name} - Confusion Matrix (noisy SNR)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"confusion_matrix_{name}.png", dpi=150)
    plt.close(fig)

    return y_pred, y_proba


def dataset_sha256(path: Path) -> str:
    """SHA-256 του αρχείου δεδομένων. Συμφωνεί με την τιμή που καταγράφει το
    runSimulation.m στο params.txt, οπότε οι μετρικές συνδέονται με το
    συγκεκριμένο αρχείο και όχι απλώς με το όνομά του."""
    return file_sha256(path)


def save_split(train_df: pd.DataFrame, test_df: pd.DataFrame, sha: str) -> dict:
    """Καταγράφει ποιες διελεύσεις πήγαν σε εκπαίδευση και ποιες σε έλεγχο.
    Χωρίς αυτό, ο διαχωρισμός αναπαράγεται μόνο εκτελώντας ξανά τον ίδιο
    κώδικα με την ίδια έκδοση της βιβλιοθήκης."""
    info = {
        "dataset_sha256": sha,
        "group_column": "PassID",
        "test_size": SPLIT_TEST_SIZE,
        "seed": SPLIT_SEED,
        "n_train_rows": int(len(train_df)),
        "n_test_rows": int(len(test_df)),
        "train_passes": sorted(int(p) for p in train_df["PassID"].unique()),
        "test_passes": sorted(int(p) for p in test_df["PassID"].unique()),
    }
    with open(RESULTS_DIR / "split.json", "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)
    return info


def save_predictions(name: str, test_df: pd.DataFrame, y_test, y_pred, y_proba) -> Path:
    """Μία γραμμή ανά δείγμα ελέγχου, με τα αναγνωριστικά του δείγματος ώστε
    κάθε μετρική να επανυπολογίζεται απευθείας από το αρχείο."""
    keys = [c for c in ("PassID", "Step", "UserID") if c in test_df.columns]
    out = test_df[keys].copy()
    out["y_true"] = y_test.to_numpy()
    out["y_pred"] = y_pred
    out["proba_satellite"] = y_proba
    path = RESULTS_DIR / f"predictions_{name}.csv"
    out.to_csv(path, index=False)
    return path


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def save_pipeline(name: str, pipeline) -> Path:
    """Αποθηκεύει ολόκληρο το pipeline, δηλαδή προεπεξεργασία και ταξινομητή
    μαζί. Η αποθήκευση μόνο του ταξινομητή θα άφηνε ανοιχτό το ενδεχόμενο
    διαφορετικού μετασχηματισμού των εισόδων κατά την επαναχρησιμοποίηση."""
    path = RESULTS_DIR / f"pipeline_{name}.joblib"
    joblib.dump(pipeline, path, compress=3)
    return path


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not DATASET_PATH.exists():
        raise FileNotFoundError(
            f"{DATASET_PATH} not found - run runSimulation.m in MATLAB first "
            "to generate the dataset (see PROD/runSimulation.m)."
        )
    df = load_dataset(DATASET_PATH)
    sigmas = load_noise_sigmas(df)
    print(f"Noise sigma (empirical, from the dataset itself): "
          f"Terrestrial={sigmas['Terrestrial']:.4f} dB, Satellite={sigmas['Satellite']:.4f} dB")
    df = add_snr_noise(df, sigmas["Terrestrial"], sigmas["Satellite"])
    train_df, test_df = group_train_test_split(df)
    sha = dataset_sha256(DATASET_PATH)
    split_info = save_split(train_df, test_df, sha)

    print(f"Loaded {len(df)} user-rows from {df['PassID'].nunique()} passes")
    print(f"Train: {len(train_df)} rows ({train_df['PassID'].nunique()} passes)")
    print(f"Test:  {len(test_df)} rows ({test_df['PassID'].nunique()} passes)")
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

    results = {"noise_sigma_dB": sigmas}
    roc_curves = {}
    results["dataset_sha256"] = sha
    results["split"] = {k: v for k, v in split_info.items()
                       if k not in ("train_passes", "test_passes")}
    for name, clf in models.items():
        pipeline = Pipeline([
            ("preprocess", build_preprocessor()),
            ("clf", clf),
        ])
        pipeline.fit(X_train, y_train)
        y_pred, y_proba = evaluate_model(name, pipeline, X_test, y_test, results)
        roc_curves[name] = (y_test, y_proba)
        pred_path = save_predictions(name, test_df, y_test, y_pred, y_proba)
        pipe_path = save_pipeline(name, pipeline)
        results[name]["predictions_file"] = pred_path.name
        results[name]["pipeline_file"] = pipe_path.name
        results[name]["pipeline_sha256"] = file_sha256(pipe_path)
        results[name]["pipeline_bytes"] = pipe_path.stat().st_size

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

            # Σπουδαιότητα χαρακτηριστικών: η impurity-based μετρική είναι
            # μεροληπτική υπέρ συνεχών/συσχετισμένων χαρακτηριστικών, οπότε
            # καταγράφεται και permutation importance πάνω στο σύνολο
            # ελέγχου (10 επαναλήψεις), το οποίο μετρά την πτώση απόδοσης
            # όταν ένα χαρακτηριστικό ανακατευθεί.
            perm = permutation_importance(pipeline, X_test, y_test, n_repeats=10,
                                          random_state=42, scoring="accuracy")
            results[name]["feature_importance_impurity"] = {
                all_feature_names[i]: float(importances[i]) for i in order
            }
            results[name]["feature_importance_permutation"] = {
                feature_cols[i]: {"mean": float(perm.importances_mean[i]),
                                  "std": float(perm.importances_std[i])}
                for i in np.argsort(perm.importances_mean)[::-1]
            }

    fig, ax = plt.subplots(figsize=(5, 5))
    for name, (y_true, y_proba) in roc_curves.items():
        RocCurveDisplay.from_predictions((y_true == "Satellite").astype(int), y_proba, name=name, ax=ax)
    ax.set_title("ROC Curve - Predicting Satellite vs Terrestrial (noisy SNR)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "roc_curve.png", dpi=150)
    plt.close(fig)

    with open(RESULTS_DIR / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved metrics + plots to {RESULTS_DIR.relative_to(ROOT)}")
    print(f"Split: {split_info['n_train_rows']} train rows "
          f"({len(split_info['train_passes'])} passes) / "
          f"{split_info['n_test_rows']} test rows "
          f"({len(split_info['test_passes'])} passes) -> split.json")
    print(f"Dataset SHA-256: {sha}")


if __name__ == "__main__":
    main()
