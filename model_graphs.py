import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Paths ---
logfile = "masif/data/masif_ppi_search/nn_models/sc05/all_feat/model_data/log.txt"

# --- Regex patterns ---
iter_pattern = re.compile(r"Iteration (\d+) validation roc auc: ([0-9\.eE+-]+)")
train_loss_pattern = re.compile(r"training_loss: ([0-9\.eE+-]+)")
mean_pos_pattern = re.compile(r"Mean validation positive score: ([0-9\.eE+-]+)")
mean_neg_pattern = re.compile(r"Mean validation negative score: ([0-9\.eE+-]+)")

# --- Lists to store values ---
iterations = []
val_aucs = []
train_loss = []
mean_pos_scores = []
mean_neg_scores = []

last_train_loss = None

with open(logfile, "r") as f:
    for line in f:
        # Training loss
        match_loss = train_loss_pattern.search(line)
        if match_loss:
            last_train_loss = float(match_loss.group(1))

        # Validation ROC
        match_iter = iter_pattern.search(line)
        if match_iter:
            iter_num = int(match_iter.group(1))
            val_auc_val = float(match_iter.group(2))
            if iter_num != 0:  # skip first iteration
                iterations.append(iter_num)
                val_aucs.append(val_auc_val)
                train_loss.append(last_train_loss if last_train_loss is not None else np.nan)

        # Mean positive score
        match_pos = mean_pos_pattern.search(line)
        if match_pos:
            mean_pos_scores.append(float(match_pos.group(1)))

        # Mean negative score
        match_neg = mean_neg_pattern.search(line)
        if match_neg:
            mean_neg_scores.append(float(match_neg.group(1)))

# --- Align lengths ---
min_len = min(len(iterations), len(val_aucs), len(train_loss),
              len(mean_pos_scores), len(mean_neg_scores))
iterations = iterations[:min_len]
val_aucs = val_aucs[:min_len]
train_loss = train_loss[:min_len]
mean_pos_scores = mean_pos_scores[:min_len]
mean_neg_scores = mean_neg_scores[:min_len]

# --- Create DataFrame ---
df = pd.DataFrame({
    "iteration": iterations,
    "train_loss": train_loss,
    "val_auc": val_aucs,
    "mean_pos_score": mean_pos_scores,
    "mean_neg_score": mean_neg_scores
})

# --- Plot Training Loss ---
plt.figure(figsize=(8,5))
plt.plot(df['iteration'], df['train_loss'], color='orange', label='Training Loss')
plt.xlabel("Iteration")
plt.ylabel("Training Loss")
plt.title("Training Loss over Iterations")
plt.grid(True)
plt.legend()
plt.savefig("training_loss.png", dpi=300)
plt.close()

# --- Plot Validation ROC AUC ---
plt.figure(figsize=(8,5))
plt.plot(df['iteration'], df['val_auc'], color='blue', label='Validation ROC AUC')
plt.xlabel("Iteration")
plt.ylabel("Validation ROC AUC")
plt.title("Validation ROC AUC over Iterations")
plt.grid(True)
plt.legend()
plt.savefig("val_auc.png", dpi=300)
plt.close()

# --- Compute ROC curve manually using mean scores ---
# Flip the scores: smaller distance = more likely positive
y_true = np.array([1]*len(mean_pos_scores) + [0]*len(mean_neg_scores))
y_scores = np.array([-s for s in mean_pos_scores] + [-s for s in mean_neg_scores])  # flip

# Sort scores descending
sorted_idx = np.argsort(y_scores)[::-1]
y_true_sorted = y_true[sorted_idx]

# Compute TPR and FPR
tpr = []
fpr = []
P = np.sum(y_true == 1)
N = np.sum(y_true == 0)
tp = 0
fp = 0
for label in y_true_sorted:
    if label == 1:
        tp += 1
    else:
        fp += 1
    tpr.append(tp / P)
    fpr.append(fp / N)

tpr = np.array(tpr)
fpr = np.array(fpr)
roc_auc = np.trapz(tpr, fpr)

# --- Plot ROC Curve ---
plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, color='blue', label=f'ROC curve (AUC = {roc_auc:.3f})')
plt.plot([0,1], [0,1], 'k--', label='Random')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve from log.txt")
plt.legend()
plt.grid(True)
plt.savefig("roc_curve.png", dpi=300)
plt.close()

print("Plots saved: training_loss.png, val_auc.png, roc_curve.png")
