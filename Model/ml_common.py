"""
Model/ml_common.py

Shared 3-class (Terrestrial/Satellite/DualConnectivity) evaluation/plotting
helpers for the three train_model*.py scripts, so their metrics stay
comparable. F1 and ROC-AUC (one-vs-rest, manual macro-average) are masked
to classes actually present in y_test, since class balance varies a lot
between experiments and sklearn's roc_auc_score would raise on an absent
class.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    auc as sk_auc,
    classification_report,
    confusion_matrix,
    f1_score,
    roc_auc_score,
    roc_curve,
)

CLASS_LABELS = ["Terrestrial", "Satellite", "DualConnectivity"]


def macro_ovr_auc(y_test, y_proba: np.ndarray, class_order: list) -> tuple[float, list]:
    """One-vs-rest ROC-AUC per class in class_order, macro-averaged over the
    classes actually present (with >=1 positive AND >=1 negative example) in
    y_test. Returns (macro_auc, list_of_classes_used)."""
    aucs = []
    used = []
    for i, cls in enumerate(class_order):
        y_true_bin = (y_test == cls).astype(int)
        if y_true_bin.nunique() < 2:
            continue
        aucs.append(roc_auc_score(y_true_bin, y_proba[:, i]))
        used.append(cls)
    macro = float(np.mean(aucs)) if aucs else float("nan")
    return macro, used


def evaluate_model(name, pipeline, X_test, y_test, results, results_dir: Path, title_suffix=""):
    y_pred = pipeline.predict(X_test)
    y_proba = pipeline.predict_proba(X_test)          # (n_samples, n_classes)
    class_order = list(pipeline.classes_)              # σειρά στηλών του y_proba

    # Μάσκα σε κλάσεις παρούσες στο y_test - μια απούσα κλάση θα τραβούσε
    # τεχνητά κάτω το macro F1 (0 precision/recall by convention).
    present_labels = [c for c in CLASS_LABELS if (y_test == c).any()]

    acc = accuracy_score(y_test, y_pred)
    f1_macro = f1_score(y_test, y_pred, average="macro", labels=present_labels)
    auc_macro, auc_classes_used = macro_ovr_auc(y_test, y_proba, class_order)
    report = classification_report(y_test, y_pred, output_dict=True, labels=CLASS_LABELS, zero_division=0)

    print(f"\n=== {name} ===")
    print(f"Accuracy: {acc:.4f}  F1(macro over {len(present_labels)}/{len(CLASS_LABELS)} classes): {f1_macro:.4f}  "
          f"ROC-AUC(macro, OvR over {len(auc_classes_used)}/{len(CLASS_LABELS)} classes): {auc_macro:.4f}")
    print(classification_report(y_test, y_pred, labels=CLASS_LABELS, zero_division=0))

    results[name] = {
        "accuracy": acc,
        "f1_macro": f1_macro,
        "f1_classes_used": present_labels,
        "roc_auc_macro_ovr": auc_macro,
        "roc_auc_classes_used": auc_classes_used,
        "classification_report": report,
    }

    cm = confusion_matrix(y_test, y_pred, labels=CLASS_LABELS)
    disp = ConfusionMatrixDisplay(cm, display_labels=CLASS_LABELS)
    fig, ax = plt.subplots(figsize=(5, 5))
    disp.plot(ax=ax, cmap="Blues", colorbar=False, xticks_rotation=20)
    ax.set_title(f"{name} - Confusion Matrix{title_suffix}")
    fig.tight_layout()
    fig.savefig(results_dir / f"confusion_matrix_{name}.png", dpi=150)
    plt.close(fig)

    return y_proba, class_order


def plot_roc_ovr(roc_data: dict, results_dir: Path, title_suffix=""):
    """roc_data: {model_name: (y_test, y_proba, class_order)}. One subplot
    per class, each showing every model's one-vs-rest ROC curve for that
    class (classes absent from a given model's y_test are skipped, matching
    macro_ovr_auc)."""
    fig, axes = plt.subplots(1, len(CLASS_LABELS), figsize=(5 * len(CLASS_LABELS), 4.5))
    for ax, cls in zip(axes, CLASS_LABELS):
        plotted_any = False
        for name, (y_test, y_proba, class_order) in roc_data.items():
            if cls not in class_order:
                continue
            col = class_order.index(cls)
            y_true_bin = (y_test == cls).astype(int)
            if y_true_bin.nunique() < 2:
                continue
            fpr, tpr, _ = roc_curve(y_true_bin, y_proba[:, col])
            ax.plot(fpr, tpr, label=f"{name} (AUC={sk_auc(fpr, tpr):.3f})")
            plotted_any = True
        ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
        ax.set_title(f"{cls} vs rest" + ("" if plotted_any else " (no test examples)"))
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        if plotted_any:
            ax.legend(loc="lower right", fontsize=8)
    fig.suptitle(f"One-vs-Rest ROC Curves{title_suffix}")
    fig.tight_layout()
    fig.savefig(results_dir / "roc_curve.png", dpi=150)
    plt.close(fig)
