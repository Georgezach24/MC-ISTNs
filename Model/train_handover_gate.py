"""
Model/train_handover_gate.py

Απόφαση handover με ML: όταν η αναφορά A3 φτάνει στον σταθμό, προβλέπεται αν
το handover θα καταλήξει σε ping-pong, ώστε ο σταθμός να περιμένει την επόμενη
αναφορά (TS 38.300 §9.2.3.2.1 βήμα 3· TR 37.817 §5.3).

Ετικέτα (TR 36.839 Ορισμός 5)
    Outcome == "HO", NextNode == FromNode, NextCause == "HO",
    StayTo_s <= pingPongMaxStay_s (0.8 s για MTS 1 s και Δt = 0.2 s).

Χωρίς διαρροή
    Ένα δείγμα ανά εντολή handover (επιτυχή ή όχι). Ό,τι συνέβη μετά την
    απόφαση χρησιμοποιείται μόνο για την ετικέτα· τα αναγνωριστικά δεν είναι
    είσοδοι (assert). Διαχωρισμός και CV ανά διέλευση· οι σπόροι εκπαίδευσης
    δεν επικαλύπτονται με την αξιολόγηση (seed 42)· το κατώφλι επιλέγεται μόνο
    στις διελεύσεις εκπαίδευσης.

Δεδομένα
    Dataset/ml_train/seed*/handover_events.mat (PROD/runSimulation.m), μόνο οι
    πολιτικές χωρίς ML (Actual, ActualMLB). Χαρακτηριστικά: handoverFeatures.m
    + FromSat, ToSat.

Μοντέλο
    RandomForest (300 δέντρα, βάθος 12)· min_samples_leaf με group CV
    (average precision). Οι τιμές που δεν ορίζονται παίρνουν SENTINEL, όπως
    και στο PROD/handoverGatePredict.m.

Κανόνας
    αναμονή αν P(ping-pong) > κατώφλι, με κατώφλι μέγιστου F1 στις out-of-fold
    προβλέψεις (Lipton, Elkan & Naryanaswamy, ECML PKDD 2014). Το 0.5 δεν
    ενεργοποιείται ποτέ, αφού το ping-pong είναι ~1% των εντολών.

Έξοδοι: Model/results_handover_gate/
    handover_gate.mat   δέντρα, κατώφλι, δείγματα αυτοελέγχου (loadHandoverGate)
    metrics.json        μετρικές σε διελεύσεις ελέγχου, CV, SHA-256 εισόδων
    *.png               precision-recall, βαθμονόμηση, σημασία χαρακτηριστικών

Χρήση:
    python train_handover_gate.py
    python train_handover_gate.py --events a.mat b.mat
    python train_handover_gate.py --plots-only   # μόνο τα PNG, χωρίς εγγραφή μοντέλου
"""

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

import matplotlib
import numpy as np
import pandas as pd
import scipy
import scipy.io
import sklearn

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.calibration import calibration_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, GroupShuffleSplit

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EVENTS = sorted((ROOT / "Dataset" / "ml_train").glob("seed*/handover_events.mat"))
RESULTS_DIR = Path(__file__).resolve().parent / "results_handover_gate"

# Οι πολιτικές χωρίς απόφαση ML: από αυτές μαθαίνει το μοντέλο.
TRAIN_POLICIES = ("Actual", "ActualMLB")

# Ίδια σειρά και ονόματα με το handoverFeatures.m (+ FromSat, ToSat).
FEATURES = [
    "NumUsers", "ReturnToPrev", "StayFrom_s",
    "SinrServ_dB", "SinrCand_dB", "Gap_dB", "DeltaSinr1s_dB", "DeltaGap1s_dB",
    "AvgSinr3s_dB", "AvgGap3s_dB", "SlopeSinr3s_dBps", "SlopeGap3s_dBps",
    "ElevFrom_deg", "ElevTo_deg", "DistFrom_m", "DistTo_m",
    "LoadFrom", "LoadTo", "UsersOnFrom", "UsersOnTo", "MlbOffset_dB",
    "FromSat", "ToSat",
]
# Ό,τι συνέβη μετά την απόφαση: μόνο για την ετικέτα.
OUTCOME_COLUMNS = {"Outcome", "StayTo_s", "NextNode", "NextCause", "Censored",
                   "TimeToPassEnd_s", "PingPong"}
# Αναγνωριστικά: ποτέ είσοδοι.
ID_COLUMNS = {"PassID", "Step", "Time_s", "UserID", "Seed", "Group", "Policy",
              "FromNode", "ToNode", "PrevNode"}
assert not (set(FEATURES) & (OUTCOME_COLUMNS | ID_COLUMNS)), "leakage: an outcome or ID column is used as a feature"

# Τιμή για «δεν ορίζεται», κάτω από κάθε πραγματική τιμή.
SENTINEL = -1000.0

N_TREES = 300
MAX_DEPTH = 12
LEAF_GRID = [1, 5, 20, 50]
CV_TREES = 100           # για τη σύγκριση ρυθμίσεων
CV_FOLDS = 3
TEST_SIZE = 0.25         # μερίδιο διελεύσεων ελέγχου
SELF_TEST_SAMPLES = 2000
RANDOM_STATE = 42


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_events(path: Path) -> tuple[pd.DataFrame, dict]:
    """handover_events.mat -> DataFrame (μία γραμμή ανά handover) + μεταδεδομένα."""
    m = scipy.io.loadmat(path, squeeze_me=True, chars_as_strings=True)
    meta = {k[len("Meta_"):]: m[k] for k in m if k.startswith("Meta_")}
    cols = {}
    for k, v in m.items():
        if k.startswith("__") or k.startswith("Meta_"):
            continue
        v = np.atleast_1d(v)
        if v.dtype == object:
            v = np.array([str(x) if not isinstance(x, np.ndarray) else "" for x in v])
        cols[k] = v
    df = pd.DataFrame(cols)
    df["Seed"] = int(meta["rngSeed"])
    df["Group"] = df["Seed"].astype(str) + "-" + df["PassID"].astype(int).astype(str)
    return df, meta


def build_frame(paths: list[Path]) -> tuple[pd.DataFrame, float, list[dict]]:
    frames, provenance, max_stay = [], [], None
    for p in paths:
        df, meta = load_events(p)
        stay = float(meta["pingPongMaxStay_s"])
        if max_stay is not None and abs(stay - max_stay) > 1e-12:
            raise ValueError(f"{p}: pingPongMaxStay_s {stay} != {max_stay}")
        max_stay = stay
        frames.append(df)
        provenance.append({
            "path": str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),
            "sha256": file_sha256(p),
            "rngSeed": int(meta["rngSeed"]),
            "events": int(len(df)),
            "passes": int(df["PassID"].nunique()),
        })
    df = pd.concat(frames, ignore_index=True)
    seeds = sorted({d["rngSeed"] for d in provenance})
    if 42 in seeds:
        raise ValueError("The training data contain the evaluation seed (42).")
    return df, max_stay, provenance


def label_and_filter(df: pd.DataFrame, max_stay: float) -> pd.DataFrame:
    df = df[df["Policy"].isin(TRAIN_POLICIES)].copy()
    df["FromSat"] = df["FromNode"].str.startswith("SAT").astype(float)
    df["ToSat"] = df["ToNode"].str.startswith("SAT").astype(float)
    df["ReturnToPrev"] = df["ReturnToPrev"].astype(float)
    ho = (df["Outcome"] == "HO").to_numpy()
    # Η διέλευση τελείωσε πριν κριθεί η ετικέτα: εκτός.
    unknown = ho & df["Censored"].astype(bool).to_numpy() & (df["StayTo_s"] <= max_stay + 1e-9).to_numpy()
    df = df[~unknown].copy()
    # Αποτυχημένη εντολή: όχι ping-pong.
    df["PingPong"] = ((df["Outcome"] == "HO") & (df["NextNode"] == df["FromNode"])
                      & (df["NextCause"] == "HO")
                      & (df["StayTo_s"] <= max_stay + 1e-9)).astype(int)
    return df


def to_matrix(df: pd.DataFrame) -> np.ndarray:
    X = df[FEATURES].to_numpy(dtype=float)
    X[~np.isfinite(X)] = SENTINEL
    return X


def make_forest(n_trees: int, leaf: int) -> RandomForestClassifier:
    return RandomForestClassifier(n_estimators=n_trees, max_depth=MAX_DEPTH,
                                  min_samples_leaf=leaf, n_jobs=-1,
                                  random_state=RANDOM_STATE)


def choose_leaf(X, y, groups) -> tuple[int, dict]:
    cv = GroupKFold(n_splits=CV_FOLDS)
    scores = {}
    for leaf in LEAF_GRID:
        ap = []
        for tr, va in cv.split(X, y, groups):
            rf = make_forest(CV_TREES, leaf).fit(X[tr], y[tr])
            ap.append(average_precision_score(y[va], rf.predict_proba(X[va])[:, 1]))
        scores[leaf] = {"mean_ap": float(np.mean(ap)), "folds": [float(a) for a in ap]}
        print(f"  min_samples_leaf={leaf:3d}: AP = {np.mean(ap):.4f} ({', '.join(f'{a:.4f}' for a in ap)})")
    best = max(scores, key=lambda k: scores[k]["mean_ap"])
    return best, scores


def choose_threshold(X, y, groups, leaf) -> tuple[float, dict]:
    """Κατώφλι μέγιστου F1 στις out-of-fold προβλέψεις των διελεύσεων εκπαίδευσης."""
    oof = np.zeros(len(y))
    for tr, va in GroupKFold(n_splits=CV_FOLDS).split(X, y, groups):
        oof[va] = make_forest(CV_TREES, leaf).fit(X[tr], y[tr]).predict_proba(X[va])[:, 1]
    prec, rec, thr = precision_recall_curve(y, oof)
    f1 = 2 * prec[:-1] * rec[:-1] / np.maximum(prec[:-1] + rec[:-1], 1e-12)
    i = int(np.argmax(f1))
    info = {"rule": "max F1, out-of-fold (training passes)", "threshold": float(thr[i]),
            "f1": float(f1[i]), "precision": float(prec[i]), "recall": float(rec[i])}
    print(f"  threshold = {thr[i]:.4f} (F1 {f1[i]:.3f}, precision {prec[i]:.3f}, recall {rec[i]:.3f})")
    return float(thr[i]), info


def export_forest(rf: RandomForestClassifier) -> dict:
    """Δέντρα σε επίπεδους πίνακες (δείκτες από 1) για το PROD/handoverGatePredict.m."""
    feat, thr, left, right, val, roots = [], [], [], [], [], []
    offset, max_depth = 0, 0
    for est in rf.estimators_:
        t = est.tree_
        n = t.node_count
        is_leaf = t.children_left == -1
        roots.append(offset + 1)
        feat.append(np.where(is_leaf, 0, t.feature + 1))
        thr.append(np.where(is_leaf, 0.0, t.threshold))
        left.append(np.where(is_leaf, 0, t.children_left + offset + 1))
        right.append(np.where(is_leaf, 0, t.children_right + offset + 1))
        v = t.value[:, 0, :]
        val.append(v[:, 1] / v.sum(axis=1))
        offset += n
        max_depth = max(max_depth, t.max_depth)
    return {
        "Roots": np.array(roots, dtype=np.int32),
        "Feature": np.concatenate(feat).astype(np.uint8),
        "Threshold": np.concatenate(thr).astype(np.float64),
        "Left": np.concatenate(left).astype(np.int32),
        "Right": np.concatenate(right).astype(np.int32),
        "Value": np.concatenate(val).astype(np.float64),
        "MaxDepth": float(max_depth),
        "Sentinel": SENTINEL,
    }


# Ύφος ίδιο με τα σχήματα του MATLAB (thesisColors.m, thesisFigure.m).
BLUE, ORANGE, NEUTRAL, TEXT, MUTED = "#2a78d6", "#eb6834", "#8a8984", "#1f1f1f", "#52514e"
CM = 1 / 2.54
plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
    "legend.fontsize": 8, "legend.frameon": False, "axes.edgecolor": TEXT,
    "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT,
    "xtick.direction": "out", "ytick.direction": "out", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": "#cccccc",
    "grid.linewidth": 0.5, "savefig.dpi": 300, "savefig.facecolor": "white",
})


def plot_pr(y, p, ap, threshold, path):
    prec, rec, _ = precision_recall_curve(y, p)
    fig, ax = plt.subplots(figsize=(11 * CM, 9 * CM))
    ax.plot(rec, prec, lw=1.5, color=BLUE, label=f"Random forest (AP = {ap:.3f})")
    ax.axhline(y.mean(), color=NEUTRAL, ls=":", lw=1.2,
               label=f"No-skill classifier (ping-pong prevalence {100*y.mean():.2f}%)")
    pt_p = precision_score(y, p > threshold, zero_division=0)
    pt_r = recall_score(y, p > threshold)
    ax.plot(pt_r, pt_p, "o", color=ORANGE, ms=6,
            label=f"Operating point, threshold {threshold:.3f}")
    ax.annotate(f"precision {pt_p:.2f}, recall {pt_r:.2f}", (pt_r, pt_p), xytext=(8, 6),
                textcoords="offset points", fontsize=7.5, color=TEXT)
    ax.set_xlabel("Recall (share of ping-pong handovers held)")
    ax.set_ylabel("Precision (share of held handovers that were ping-pong)")
    ax.set_title("Ping-pong prediction on held-out passes", fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_reliability(y, p, threshold, path):
    frac, mean_p = calibration_curve(y, p, n_bins=10, strategy="quantile")
    fig, ax = plt.subplots(figsize=(11 * CM, 9 * CM))
    hi = 1.1 * max(mean_p.max(), frac.max())
    ax.plot([0, hi], [0, hi], ":", color=NEUTRAL, lw=1.2, label="Perfect calibration")
    ax.plot(mean_p, frac, "o-", color=BLUE, ms=4, lw=1.5, label="Random forest")
    ax.axvline(threshold, color=ORANGE, ls="--", lw=1, label=f"Decision threshold {threshold:.3f}")
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xlabel("Predicted ping-pong probability (decile mean)")
    ax.set_ylabel("Observed ping-pong frequency")
    ax.set_title("Calibration on held-out passes (prediction deciles)", fontweight="bold")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_importance(imp, path):
    order = np.argsort(imp.importances_mean)
    fig, ax = plt.subplots(figsize=(14 * CM, 13 * CM))
    ax.barh(np.array(FEATURES)[order], imp.importances_mean[order],
            xerr=imp.importances_std[order], color=BLUE, height=0.7,
            error_kw={"ecolor": TEXT, "elinewidth": 0.8, "capsize": 2})
    ax.axvline(0, color=TEXT, lw=0.8)
    ax.set_xlabel("Decrease of average precision when the feature is permuted")
    ax.set_title("Permutation feature importance (held-out passes, 5 repeats)", fontweight="bold")
    ax.grid(True, axis="x")
    ax.grid(False, axis="y")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plots_only(args, Xtr, ytr, Xte, yte):
    """Αναπαραγωγή των PNG από το αποθηκευμένο μοντέλο, χωρίς εγγραφή .mat/.json."""
    with open(args.out / "metrics.json", encoding="utf-8") as f:
        saved = json.load(f)
    leaf = int(saved["model"]["min_samples_leaf"])
    threshold = float(saved["threshold"]["threshold"])
    gate = scipy.io.loadmat(args.out / "handover_gate.mat", squeeze_me=True)
    rf = make_forest(N_TREES, leaf).fit(Xtr, ytr)
    tx = np.array(gate["TestX"], dtype=float)
    tx[~np.isfinite(tx)] = SENTINEL
    err = float(np.max(np.abs(rf.predict_proba(tx)[:, 1] - gate["TestP"])))
    if err > 1e-9:
        sys.exit(f"The refitted forest differs from handover_gate.mat (max |dP| = {err:.3g}).")
    print(f"Refitted forest identical to handover_gate.mat (max |dP| = {err:.1e}, "
          f"min_samples_leaf = {leaf}, threshold = {threshold:.4f})")
    p = rf.predict_proba(Xte)[:, 1]
    plot_pr(yte, p, average_precision_score(yte, p), threshold, args.out / "precision_recall.png")
    plot_reliability(yte, p, threshold, args.out / "reliability.png")
    if args.importance:
        imp = permutation_importance(rf, Xte, yte, scoring="average_precision",
                                     n_repeats=5, random_state=RANDOM_STATE, n_jobs=-1)
        plot_importance(imp, args.out / "feature_importance.png")
    print(f"Figures -> {args.out}")


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--events", nargs="*", type=Path, default=DEFAULT_EVENTS)
    ap_.add_argument("--out", type=Path, default=RESULTS_DIR)
    ap_.add_argument("--cv", action=argparse.BooleanOptionalAction, default=True,
                     help="select min_samples_leaf by group CV (otherwise 1)")
    ap_.add_argument("--plots-only", action="store_true",
                     help="refit the saved configuration and write only the figures")
    ap_.add_argument("--importance", action=argparse.BooleanOptionalAction, default=True)
    args = ap_.parse_args()
    if args.plots_only and args.events == DEFAULT_EVENTS:
        # Τα ίδια αρχεία με το αποθηκευμένο μοντέλο (provenance του metrics.json).
        with open(args.out / "metrics.json", encoding="utf-8") as f:
            args.events = [ROOT / Path(d["path"].replace("\\", "/")) for d in json.load(f)["data"]["files"]]
    if not args.events:
        sys.exit("No handover_events.mat files found (Dataset/ml_train/seed*/).")
    args.out.mkdir(parents=True, exist_ok=True)

    raw, max_stay, provenance = build_frame(args.events)
    df = label_and_filter(raw, max_stay)
    print(f"Handover commands (Actual, ActualMLB): {len(df)} from {df['Group'].nunique()} passes, "
          f"ping-pong {df['PingPong'].sum()} ({100*df['PingPong'].mean():.2f}%)")
    print(df.groupby(["Policy", "Outcome"])["PingPong"].agg(["size", "sum", "mean"]))

    gss = GroupShuffleSplit(n_splits=1, test_size=TEST_SIZE, random_state=RANDOM_STATE)
    tr_idx, te_idx = next(gss.split(df, df["PingPong"], df["Group"]))
    train, test = df.iloc[tr_idx], df.iloc[te_idx]
    Xtr, ytr, gtr = to_matrix(train), train["PingPong"].to_numpy(), train["Group"].to_numpy()
    Xte, yte = to_matrix(test), test["PingPong"].to_numpy()
    print(f"Training: {len(train)} ({train['Group'].nunique()} passes, {ytr.sum()} ping-pong) | "
          f"test: {len(test)} ({test['Group'].nunique()} passes, {yte.sum()} ping-pong)")
    if args.plots_only:
        plots_only(args, Xtr, ytr, Xte, yte)
        return

    if args.cv:
        print(f"Selecting min_samples_leaf ({CV_FOLDS}-fold group CV, {CV_TREES} trees):")
        leaf, cv_scores = choose_leaf(Xtr, ytr, gtr)
    else:
        leaf, cv_scores = 1, {}
    print(f"  -> min_samples_leaf = {leaf}")
    threshold, thr_info = choose_threshold(Xtr, ytr, gtr, leaf)

    rf = make_forest(N_TREES, leaf).fit(Xtr, ytr)
    p = rf.predict_proba(Xte)[:, 1]
    yhat = (p > threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(yte, yhat, labels=[0, 1]).ravel()
    metrics_test = {
        "n": int(len(yte)), "positives": int(yte.sum()),
        "roc_auc": float(roc_auc_score(yte, p)),
        "average_precision": float(average_precision_score(yte, p)),
        "brier": float(brier_score_loss(yte, p)),
        "threshold": threshold,
        "precision": float(precision_score(yte, yhat, zero_division=0)),
        "recall": float(recall_score(yte, yhat)),
        "confusion": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
        "held_fraction": float(yhat.mean()),
    }
    by_pol = {}
    for pol in TRAIN_POLICIES:
        m = (test["Policy"] == pol).to_numpy()
        if m.any() and yte[m].sum() > 0:
            by_pol[pol] = {
                "n": int(m.sum()), "positives": int(yte[m].sum()),
                "average_precision": float(average_precision_score(yte[m], p[m])),
                "precision": float(precision_score(yte[m], yhat[m], zero_division=0)),
                "recall": float(recall_score(yte[m], yhat[m])),
            }
    metrics_test["by_policy"] = by_pol
    print(json.dumps(metrics_test, indent=2, ensure_ascii=False))

    plot_pr(yte, p, metrics_test["average_precision"], threshold, args.out / "precision_recall.png")
    plot_reliability(yte, p, threshold, args.out / "reliability.png")
    importance = None
    if args.importance:
        imp = permutation_importance(rf, Xte, yte, scoring="average_precision",
                                     n_repeats=5, random_state=RANDOM_STATE, n_jobs=-1)
        plot_importance(imp, args.out / "feature_importance.png")
        importance = {f: {"mean": float(m), "std": float(s)} for f, m, s in
                      zip(FEATURES, imp.importances_mean, imp.importances_std)}

    # Δείγματα αυτοελέγχου για τη handoverGatePredict (όλα τα θετικά + τυχαία).
    rng = np.random.default_rng(RANDOM_STATE)
    pos = np.flatnonzero(yte == 1)
    neg = rng.choice(np.flatnonzero(yte == 0), size=max(0, SELF_TEST_SAMPLES - len(pos)), replace=False)
    sel = np.concatenate([pos, neg])[:SELF_TEST_SAMPLES]
    test_x_raw = test[FEATURES].to_numpy(dtype=float)[sel]
    test_p = p[sel]

    forest = export_forest(rf)
    meta = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "sklearn": sklearn.__version__, "numpy": np.__version__, "scipy": scipy.__version__,
        "label": "TR 36.839 Ορισμός 5: NextNode == FromNode & NextCause == HO & StayTo_s <= pingPongMaxStay_s",
        "pingPongMaxStay_s": max_stay,
        "trainPolicies": ", ".join(TRAIN_POLICIES),
        "trainSeeds": ", ".join(str(d["rngSeed"]) for d in provenance),
        "trainEvents": int(len(train)), "trainPingPong": int(ytr.sum()),
        "nTrees": N_TREES, "maxDepth": MAX_DEPTH, "minSamplesLeaf": int(leaf),
        "testAveragePrecision": metrics_test["average_precision"],
        "testRecallAtThreshold": metrics_test["recall"],
        "testPrecisionAtThreshold": metrics_test["precision"],
    }
    gate_path = args.out / "handover_gate.mat"
    scipy.io.savemat(gate_path, {
        **forest,
        "FeatureNames": np.array(FEATURES, dtype=object),
        "Threshold_p": threshold,
        "TestX": test_x_raw, "TestP": test_p,
        "Meta": meta,
    }, do_compression=True)

    out = {
        "data": {"files": provenance, "pingPongMaxStay_s": max_stay,
                 "policies": list(TRAIN_POLICIES),
                 "events": int(len(df)), "passes": int(df["Group"].nunique()),
                 "pingPong": int(df["PingPong"].sum())},
        "split": {"test_size": TEST_SIZE, "by": "pass (seed, PassID)",
                  "train_passes": int(train["Group"].nunique()),
                  "test_passes": int(test["Group"].nunique())},
        "model": {"type": "RandomForestClassifier", "n_estimators": N_TREES,
                  "max_depth": MAX_DEPTH, "min_samples_leaf": int(leaf),
                  "sentinel": SENTINEL, "features": FEATURES,
                  "nodes": int(len(forest["Feature"]))},
        "cv": {str(k): v for k, v in cv_scores.items()},
        "threshold": thr_info,
        "test": metrics_test,
        "permutation_importance": importance,
        "export": {"path": str(gate_path.relative_to(ROOT)), "sha256": file_sha256(gate_path),
                   "bytes": gate_path.stat().st_size, "self_test_samples": int(len(sel))},
        "meta": meta,
    }
    with open(args.out / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"Model: {gate_path} ({gate_path.stat().st_size/1e6:.1f} MB, "
          f"{len(forest['Feature'])} nodes)")


if __name__ == "__main__":
    main()
