import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Paths ---
logfile = "masif/data/masif_ppi_search/nn_models/sc05/all_feat/model_data/log.txt"
pos_file = "pos_dists.npy"
neg_file = "neg_dists.npy"

# --- Regex patterns ---
iter_pattern = re.compile(r"Iteration (\d+) validation roc auc: ([0-9\.eE+-]+)")
train_loss_pattern = re.compile(r"training_loss: ([0-9\.eE+-]+)")
mean_pos_pattern = re.compile(r"Mean validation positive score: ([0-9\.eE+-]+)")
mean_neg_pattern = re.compile(r"Mean validation negative score: ([0-9\.eE+-]+)")

# --- Lists to store values ---
iterations = []
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

        # Validation iteration → append loss
        match_iter = iter_pattern.search(line)
        if match_iter:
            iter_num = int(match_iter.group(1))
            if iter_num != 0:
                iterations.append(iter_num)
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
min_len = min(len(iterations), len(train_loss), len(mean_pos_scores), len(mean_neg_scores))
iterations = iterations[:min_len]
train_loss = train_loss[:min_len]
mean_pos_scores = mean_pos_scores[:min_len]
mean_neg_scores = mean_neg_scores[:min_len]

# --- Create DataFrame ---
df = pd.DataFrame({
    "iteration": iterations,
    "train_loss": train_loss,
    "mean_pos_score": mean_pos_scores,
    "mean_neg_score": mean_neg_scores
})

# --- Plot training loss ---
plt.figure(figsize=(8,5))
plt.plot(df['iteration'], df['train_loss'], marker='o', label='Training Loss', color='orange')
plt.xlabel("Iteration")
plt.ylabel("Training Loss")
plt.title("Training Loss over Iterations")
plt.grid(True)
plt.legend()
plt.savefig("training_loss.png", dpi=300)
plt.close()

# --- Plot mean positive vs negative distances ---
plt.figure(figsize=(8,5))
plt.plot(df['iteration'], df['mean_pos_score'], marker='o', label='Mean positive score', color='blue')
plt.plot(df['iteration'], df['mean_neg_score'], marker='o', label='Mean negative score', color='red')
plt.xlabel("Iteration")
plt.ylabel("Mean distance")
plt.title("Mean positive vs negative scores over iterations")
plt.grid(True)
plt.legend()
plt.savefig("mean_scores_over_iterations.png", dpi=300)
plt.close()

# --- Scatter plot ---
plt.figure(figsize=(6,6))
plt.scatter(df['mean_neg_score'], df['mean_pos_score'], c=df['iteration'], cmap='viridis')
plt.xlabel("Mean negative score")
plt.ylabel("Mean positive score")
plt.title("Mean positive vs negative scores (colour = iteration)")
plt.grid(True)
plt.savefig("pos_vs_neg_scatter.png", dpi=300)
plt.close()

print("Saved: training_loss.png, mean_scores_over_iterations.png, pos_vs_neg_scatter.png")

# ====================================================
# === NEW SECTION: Load pos/neg distances & ROC AUC ===
# ====================================================

print("Loading distance files for ROC AUC...")

posd = np.load(pos_file)
negd = np.load(neg_file)

ytrue = np.concatenate([np.ones_like(posd), np.zeros_like(negd)])
ypred = 1.0 / np.concatenate([posd, negd])  # smaller dist = higher score

# --- Manual ROC computation ---
combined = np.vstack([ypred, ytrue]).T
combined = combined[combined[:,0].argsort()[::-1]]  # sort by score descending

tpr = []
fpr = []

P = np.sum(ytrue == 1)
N = np.sum(ytrue == 0)

tp = 0
fp = 0

for _, label in combined:
    if label == 1:
        tp += 1
    else:
        fp += 1
    tpr.append(tp / P)
    fpr.append(fp / N)

tpr = np.array(tpr)
fpr = np.array(fpr)

roc_auc = np.trapz(tpr, fpr)
print(f"ROC AUC = {roc_auc:.4f}")

# --- Plot ROC curve ---
plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.3f}")
plt.plot([0,1], [0,1], 'k--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve (from pos_dists.npy and neg_dists.npy)")
plt.grid(True)
plt.legend()
plt.savefig("roc_auc_curve.png", dpi=300)
plt.close()

print("Saved: roc_auc_curve.png")
