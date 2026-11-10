"""
Model/train_handover_gate.py

ML-assisted handover decision (Part 2): at the moment an A3 measurement report
reaches the serving station, predict whether executing the handover would end
in a ping-pong, so that the station can wait for the next report instead.

Why this is the decision point
    TS 38.300 §9.2.3.2.1 step 3: "The source gNB decides to handover the UE,
    based on MeasurementReport and RRM information". The A3 event only
    triggers a report; the decision is left to the station, and TR 37.817 §5.3
    studies AI/ML-based mobility optimisation exactly for reducing ping-pong
    and handover failures, with the inputs listed in §5.3.2.4.

Label (TR 36.839 Definition 5, the simulation's own count)
    ping-pong = the user returns to the source node with a handover after
    staying on the target less than the MTS (1 s):
        Outcome == "HO" and NextNode == FromNode and NextCause == "HO"
        and StayTo_s <= pingPongMaxStay_s
    (pingPongMaxStay_s is read from the events file: 0.8 s at dt = 0.2 s).

No outcome information in the inputs
    One sample per handover COMMAND, successful or not: success is decided
    after the command (true SINR of the target, T310), so keeping only the
    successful ones would select samples by the outcome. Everything that
    happened after the decision (Outcome, StayTo_s, NextNode, NextCause,
    Censored, TimeToPassEnd_s) is used only to build the label; identifiers
    (pass, step, user, seed) are never inputs. Both are checked by an assert.
    Train/test split and CV are by pass; the evaluation run (seed 42) shares
    no pass with the training seeds, and the threshold is chosen on the
    training passes only (out-of-fold), never on the test passes or seed 42.

Data
    Dataset/ml_train/seed*/handover_events.mat, written by PROD/runSimulation.m
    with rngSeed values that do not overlap the evaluation run (seed 42:
    passes use rng(seed + passIdx)). Only the policies WITHOUT the ML decision
    (Actual, ActualMLB) are used, so the model learns what the plain A3 does.
    Features are the ones logged at the decision (handoverFeatures.m) plus the
    node types; their sources are listed in runSimulation.m, logHandoverEvents.

Model
    RandomForestClassifier, as in the Part 2 proof of concept (300 trees,
    depth 12). min_samples_leaf is chosen by grouped cross-validation on the
    training passes (average precision, the summary of the precision-recall
    curve that suits a rare positive class). Missing values (elevation of a
    base station, history before the pass started) get a sentinel below any
    real value, identically in PROD/handoverGatePredict.m.

Decision rule
    wait if P(ping-pong) > threshold. Ping-pong is rare (~1% of commands), so
    the minimum-error rule (0.5) never fires: no prediction reaches it and the
    gate would do nothing. The threshold is the one that maximises F1 on the
    out-of-fold predictions of the training passes (Lipton, Elkan &
    Naryanaswamy, "Optimal Thresholding of Classifiers to Maximize F1
    Measure", ECML PKDD 2014): the usual choice for a rare positive class.
    It is chosen with the hyperparameters, before the test passes are touched,
    and never on the evaluation run.

Output: Model/results_handover_gate/
    handover_gate.mat   trees + threshold + self-test samples, read by
                        runSimulation.m (loadHandoverGate)
    metrics.json        offline metrics on held-out passes, CV results, data
                        provenance (SHA-256 of every input file)
    *.png               precision-recall, reliability, feature importance

Usage:
    python train_handover_gate.py
    python train_handover_gate.py --events path/a.mat path/b.mat
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
# Ό,τι συνέβη ΜΕΤΑ την απόφαση: μόνο για την ετικέτα, ποτέ ως είσοδος. Το
# μοντέλο ξέρει μόνο ό,τι ξέρει ο σταθμός τη στιγμή της αναφοράς.
OUTCOME_COLUMNS = {"Outcome", "StayTo_s", "NextNode", "NextCause", "Censored",
                   "TimeToPassEnd_s", "PingPong"}
# Αναγνωριστικά (διέλευση, βήμα, χρήστης, σπόρος): όχι είσοδοι, θα έδιναν
# στο μοντέλο ταυτότητα δείγματος αντί για κατάσταση.
ID_COLUMNS = {"PassID", "Step", "Time_s", "UserID", "Seed", "Group", "Policy",
              "FromNode", "ToNode", "PrevNode"}
assert not (set(FEATURES) & (OUTCOME_COLUMNS | ID_COLUMNS)), "διαρροή στήλης στις εισόδους"

# Κάτω από κάθε πραγματική τιμή (SINR > -300 dB, ανύψωση >= 0, απόσταση > 0,
# φόρτος >= 0, μετατόπιση >= -24 dB): ένας κόμβος του δέντρου τη στέλνει πάντα
# αριστερά, δηλαδή «δεν ορίζεται» είναι μια ξεχωριστή κατηγορία.
SENTINEL = -1000.0

N_TREES = 300            # όπως το Part 2 (train_model.py)
MAX_DEPTH = 12           # όπως το Part 2 (train_model.py)
LEAF_GRID = [1, 5, 20, 50]
CV_TREES = 100           # αρκεί για τη σύγκριση των ρυθμίσεων
CV_FOLDS = 3
TEST_SIZE = 0.25         # ανά διέλευση, όπως το Part 2
SELF_TEST_SAMPLES = 2000
RANDOM_STATE = 42


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_events(path: Path) -> tuple[pd.DataFrame, dict]:
    """One handover_events.mat -> DataFrame (one row per handover) + metadata."""
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
        raise ValueError("Τα δεδομένα εκπαίδευσης περιέχουν τον σπόρο της αξιολόγησης (42).")
    return df, max_stay, provenance


def label_and_filter(df: pd.DataFrame, max_stay: float) -> pd.DataFrame:
    df = df[df["Policy"].isin(TRAIN_POLICIES)].copy()
    df["FromSat"] = df["FromNode"].str.startswith("SAT").astype(float)
    df["ToSat"] = df["ToNode"].str.startswith("SAT").astype(float)
    df["ReturnToPrev"] = df["ReturnToPrev"].astype(float)
    # Μία γραμμή ανά εντολή handover, επιτυχή ή όχι (runSimulation,
    # logHandoverEvents): το μοντέλο ρωτιέται πριν φανεί η έκβαση, οπότε
    # μαθαίνει από όλες. Αρχεία πριν τη στήλη Outcome είχαν μόνο επιτυχείς.
    if "Outcome" not in df:
        df["Outcome"] = "HO"
    ho = (df["Outcome"] == "HO").to_numpy()
    # Η διέλευση τελείωσε πριν φανεί αν ο χρήστης επιστρέφει μέσα στο MTS:
    # η ετικέτα δεν είναι γνωστή, οπότε το handover δεν χρησιμοποιείται.
    unknown = ho & df["Censored"].astype(bool).to_numpy() & (df["StayTo_s"] <= max_stay + 1e-9).to_numpy()
    df = df[~unknown].copy()
    # Αποτυχημένη εντολή: ο χρήστης δεν έφτασε στον στόχο, άρα δεν υπάρχει
    # επιστροφή από αυτόν - όχι ping-pong.
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
    print(f"  κατώφλι = {thr[i]:.4f} (F1 {f1[i]:.3f}, precision {prec[i]:.3f}, recall {rec[i]:.3f})")
    return float(thr[i]), info


def export_forest(rf: RandomForestClassifier) -> dict:
    """Trees in flat arrays, 1-based, as read by PROD/handoverGatePredict.m."""
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


def plot_pr(y, p, ap, threshold, path):
    prec, rec, _ = precision_recall_curve(y, p)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(rec, prec, lw=1.5, label=f"RandomForest (AP = {ap:.3f})")
    ax.axhline(y.mean(), color="k", ls=":", label=f"τυχαίο (= ποσοστό ping-pong {100*y.mean():.2f}%)")
    pt_p = precision_score(y, p > threshold, zero_division=0)
    pt_r = recall_score(y, p > threshold)
    ax.plot(pt_r, pt_p, "o", color="C3", label=f"κατώφλι {threshold:.3f}")
    ax.set_xlabel("Recall (ping-pong που προλαμβάνονται)")
    ax.set_ylabel("Precision (αναμονές που ήταν όντως ping-pong)")
    ax.set_title("Πρόβλεψη ping-pong σε διελεύσεις εκτός εκπαίδευσης")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_reliability(y, p, threshold, path):
    frac, mean_p = calibration_curve(y, p, n_bins=10, strategy="quantile")
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot([0, 1], [0, 1], "k:", label="τέλεια βαθμονόμηση")
    ax.plot(mean_p, frac, "o-", label="RandomForest")
    ax.axvline(threshold, color="C3", ls="--", lw=1, label=f"κατώφλι {threshold:.3f}")
    ax.set_xlabel("Προβλεπόμενη πιθανότητα ping-pong")
    ax.set_ylabel("Παρατηρούμενο ποσοστό ping-pong")
    ax.set_title("Βαθμονόμηση (δεκατημόρια πρόβλεψης)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_importance(imp, path):
    order = np.argsort(imp.importances_mean)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.barh(np.array(FEATURES)[order], imp.importances_mean[order],
            xerr=imp.importances_std[order], color="C0")
    ax.set_xlabel("Μείωση average precision όταν ανακατεύεται η στήλη")
    ax.set_title("Σημασία χαρακτηριστικών (permutation, διελεύσεις ελέγχου)")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--events", nargs="*", type=Path, default=DEFAULT_EVENTS)
    ap_.add_argument("--out", type=Path, default=RESULTS_DIR)
    ap_.add_argument("--cv", action=argparse.BooleanOptionalAction, default=True,
                     help="επιλογή min_samples_leaf με group CV (αλλιώς 1)")
    ap_.add_argument("--importance", action=argparse.BooleanOptionalAction, default=True)
    args = ap_.parse_args()
    if not args.events:
        sys.exit("Δεν βρέθηκαν αρχεία handover_events.mat (Dataset/ml_train/seed*/).")
    args.out.mkdir(parents=True, exist_ok=True)

    raw, max_stay, provenance = build_frame(args.events)
    df = label_and_filter(raw, max_stay)
    print(f"Εντολές handover (Actual, ActualMLB): {len(df)} από {df['Group'].nunique()} διελεύσεις, "
          f"ping-pong {df['PingPong'].sum()} ({100*df['PingPong'].mean():.2f}%)")
    print(df.groupby(["Policy", "Outcome"])["PingPong"].agg(["size", "sum", "mean"]))

    gss = GroupShuffleSplit(n_splits=1, test_size=TEST_SIZE, random_state=RANDOM_STATE)
    tr_idx, te_idx = next(gss.split(df, df["PingPong"], df["Group"]))
    train, test = df.iloc[tr_idx], df.iloc[te_idx]
    Xtr, ytr, gtr = to_matrix(train), train["PingPong"].to_numpy(), train["Group"].to_numpy()
    Xte, yte = to_matrix(test), test["PingPong"].to_numpy()
    print(f"Εκπαίδευση: {len(train)} ({train['Group'].nunique()} διελεύσεις, {ytr.sum()} ping-pong) | "
          f"έλεγχος: {len(test)} ({test['Group'].nunique()} διελεύσεις, {yte.sum()} ping-pong)")

    if args.cv:
        print(f"Επιλογή min_samples_leaf ({CV_FOLDS}-fold group CV, {CV_TREES} δέντρα):")
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

    # Δείγματα αυτοελέγχου: η handoverGatePredict πρέπει να δώσει τις ίδιες
    # πιθανότητες. Όλα τα θετικά του ελέγχου και τυχαία από τα υπόλοιπα, με
    # τις τιμές NaN όπως έρχονται από την προσομοίωση.
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
    print(f"Μοντέλο: {gate_path} ({gate_path.stat().st_size/1e6:.1f} MB, "
          f"{len(forest['Feature'])} κόμβοι)")


if __name__ == "__main__":
    main()
