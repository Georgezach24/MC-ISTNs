"""
Model/train_model.py

Baseline ML proof-of-concept for Part 2 of the thesis: given per-user link
candidate features (best-BS and satellite SNR/distance/elevation, produced by
PROD/runSimulation.m), predict which node TYPE (Terrestrial vs Satellite)
the SNR-greedy baseline in simulateScenario.m would select.

This is a sanity-check model, not the final Part 2 deliverable: since the
label is essentially argmax(CandBS_SNR_dB, CandSat_SNR_dB), a model given both
candidate SNRs is expected to reproduce the rule almost perfectly. The point
of this run is to validate the dataset pipeline end-to-end (MATLAB -> CSV ->
Python -> trained model -> metrics) before tackling harder Part 2 targets
(e.g. predicting from imperfect/estimated SNR, joint/fair allocation, or
multi-KPI objectives).

Usage:
    python train_model.py
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
RESULTS_DIR = Path(__file__).resolve().parent / "results"

# Ίδιο κατώφλι με satParameters.MinElevationDeg στο runSimulation.m
# (TR 38.821 visibility mask) - όχι μια νέα υπόθεση, απλά επαναχρησιμοποίηση.
MIN_ELEVATION_DEG = 20.0
SENTINEL_SNR_DB = -50.0        # "πρακτικά άχρηστος" όταν ο δορυφόρος δεν είναι ορατός
SENTINEL_PATHLOSS_DB = 300.0

# Η υποδειγματοληψία στον χρόνο γίνεται ΣΤΗΝ ΠΗΓΗ: το runSimulation.m
# προχωρά με βήμα 1 s αλλά γράφει μία γραμμή κάθε datasetStride βήματα
# (5 s), γιατί διαδοχικά δείγματα του ίδιου χρήστη απέχουν 0.83 m και είναι
# σχεδόν ταυτόσημα. Εδώ δεν χρειάζεται επιπλέον αραίωση· άλλαξε το μόνο αν
# θέλεις ακόμη πιο αραιό δείγμα.
ML_SAMPLE_STRIDE = 1

FEATURE_COLUMNS_NUMERIC = [
    # NodeLoad is deliberately excluded: it's the count of users sharing the
    # SAME winning node within a scenario, which is a downstream consequence
    # of ServingType for every user in that scenario (satellite scenarios
    # mechanically have larger groups) - a circular predictor, not a cause.
    "NumUsers",
    "CandBS_SNR_dB", "CandBS_Distance_m", "CandBS_PathLoss_dB",
    "CandSat_SNR_dB", "CandSat_Elevation_deg", "CandSat_SlantRange_m",
    "CandSat_PathLoss_dB",
]
FEATURE_COLUMNS_CATEGORICAL = []  # το σενάριο διάδοσης είναι σταθερό (UMa)
FEATURE_COLUMNS_BOOL = ["CandSat_Visible"]
TARGET_COLUMN = "ServingType"

# Παράμετροι διαχωρισμού, ορισμένοι μία φορά ώστε το split.json να μην
# μπορεί να διαφωνήσει με τον διαχωρισμό που έγινε στην πραγματικότητα.
SPLIT_TEST_SIZE = 0.25
SPLIT_SEED = 42


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Υποδειγματοληψία στον χρόνο (βλ. ML_SAMPLE_STRIDE).
    if ML_SAMPLE_STRIDE > 1 and "Step" in df.columns:
        before = len(df)
        df = df[df["Step"] % ML_SAMPLE_STRIDE == 1].reset_index(drop=True)
        print(f"Time decimation: kept {len(df)} of {before} rows "
              f"(1 sample every {ML_SAMPLE_STRIDE} s)")
    # Ο δορυφόρος έχει -Inf SNR / Inf path loss όταν elevation < MinElevationDeg
    # (visibility mask στο simulateScenario.m). Αντικατάσταση με sentinel τιμές
    # + ρητό boolean flag, ώστε το μοντέλο να μη σκάει σε μη-πεπερασμένες τιμές.
    df["CandSat_Visible"] = df["CandSat_Elevation_deg"] >= MIN_ELEVATION_DEG
    df["CandSat_SNR_dB"] = df["CandSat_SNR_dB"].replace([np.inf, -np.inf], SENTINEL_SNR_DB)
    df["CandSat_PathLoss_dB"] = df["CandSat_PathLoss_dB"].replace([np.inf, -np.inf], SENTINEL_PATHLOSS_DB)

    # simulateScenario.m πλέον καταγράφει και ServingType="Outage" (κανένας
    # υποψήφιος δεν ξεπερνά το ελάχιστο χρησιμοποιήσιμο SNR). Εξαιρούνται
    # εδώ: το "ποιος από τους δύο διαθέσιμους κόμβους κερδίζει" είναι
    # διαφορετικό ερώτημα από το "υπάρχει καθόλου κάλυψη" - η ανάμειξή τους
    # θα αλλοίωνε το ήδη καθιερωμένο binary πρόβλημα Terrestrial/Satellite.
    numOutage = int((df["ServingType"] == "Outage").sum())
    if numOutage:
        print(f"Excluding {numOutage} Outage rows (no candidate above minimum usable SNR) "
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
    return df


SNR_MIN_DB = 10 * np.log10(2 ** 0.2344 - 1)   # -7.5346 dB, MCS 0 (TS 38.214 Πίν. 5.1.3.1-1)


def label_identity_check(path: Path) -> dict:
    """Έλεγχος ταυτότητας της ετικέτας (σημείο 15 της αξιολόγησης).

    Εφαρμόζει τον ΙΔΙΟ τον κανόνα δημιουργίας των ετικετών απευθείας στα δύο
    υποψήφια SNR -- μάσκα ορατότητας, argmax με τις ισοπαλίες να πηγαίνουν στον
    επίγειο (στο simulateScenario.m ο δορυφόρος κερδίζει μόνο με `>`), και το
    κατώφλι ελάχιστου χρησιμοποιήσιμου SNR -- και το συγκρίνει με την
    καταγεγραμμένη ετικέτα. Τρέχει στα ΑΚΑΤΕΡΓΑΣΤΑ δεδομένα, πριν από κάθε
    φιλτράρισμα, ώστε να καλύπτει και τις γραμμές εκτός κάλυψης.

    Αν ο κανόνας δεν αναπαράγει τις ετικέτες, υπάρχει ασυνέπεια δεδομένων ή
    υλοποίησης. Αν τις αναπαράγει, τότε η ετικέτα είναι εξ ορισμού συνάρτηση
    δύο χαρακτηριστικών εισόδου, και η ακρίβεια της παραλλαγής με ακριβές SNR
    είναι έλεγχος ροής δεδομένων, όχι αποτέλεσμα πρόβλεψης.
    """
    raw = pd.read_csv(path, usecols=["CandBS_SNR_dB", "CandSat_SNR_dB",
                                     "CandSat_Elevation_deg", "ServingType"])
    bs = raw["CandBS_SNR_dB"].fillna(-np.inf).to_numpy()
    sat = raw["CandSat_SNR_dB"].to_numpy()
    sat = np.where(raw["CandSat_Elevation_deg"].to_numpy() >= MIN_ELEVATION_DEG, sat, -np.inf)
    sat = np.where(np.isnan(sat), -np.inf, sat)

    sat_wins = sat > bs                      # ισοπαλία -> επίγειος, όπως στη MATLAB
    best = np.where(sat_wins, sat, bs)
    predicted = np.where(best < SNR_MIN_DB, "Outage",
                         np.where(sat_wins, "Satellite", "Terrestrial"))

    actual = raw["ServingType"].to_numpy()
    matches = int((predicted == actual).sum())
    total = len(actual)
    print(f"\nLabel identity check (SNR_min = {SNR_MIN_DB:.4f} dB): "
          f"rule reproduces {matches}/{total} labels "
          f"({100*matches/total:.4f}%, {total-matches} mismatches)")
    return {"snr_min_db": float(SNR_MIN_DB), "rows": total,
            "matches": matches, "mismatches": total - matches}


def build_preprocessor() -> ColumnTransformer:
    return ColumnTransformer([
        ("num", StandardScaler(), FEATURE_COLUMNS_NUMERIC),
        ("cat", OneHotEncoder(drop="if_binary"),
         FEATURE_COLUMNS_CATEGORICAL + FEATURE_COLUMNS_BOOL),
    ])


def group_train_test_split(df: pd.DataFrame, test_size=SPLIT_TEST_SIZE,
                          seed=SPLIT_SEED):
    # Split ανά PassID (όχι ανά γραμμή): χρήστες του ίδιου σεναρίου
    # μοιράζονται τις ίδιες θέσεις BS/δορυφόρου, άρα ένα row-level split θα
    # διέρρεε γεωμετρία σεναρίου ανάμεσα σε train/test.
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
    ax.set_title(f"{name} - Confusion Matrix")
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

    label_check = label_identity_check(DATASET_PATH)

    df = load_dataset(DATASET_PATH)
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

    results = {"label_identity_check": label_check}
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
            ax.set_title("RandomForest Feature Importance")
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
    ax.set_title("ROC Curve - Predicting Satellite vs Terrestrial")
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
