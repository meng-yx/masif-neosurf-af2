#!/usr/bin/env python3
import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import tensorflow as tf
from sklearn.metrics import average_precision_score, roc_auc_score

try:
    from protein_interaction.dataset import (
        ProteinPairDataset,
        batch_iterator,
        read_pairs_csv,
        sanity_check_batch,
    )
    from protein_interaction.models.cross_set_v1 import CrossSetInteractionV1
    from protein_interaction.utils import ensure_dir, load_config, resolve_path
except ModuleNotFoundError:
    # Allow direct execution: python path/to/train.py ...
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from protein_interaction.dataset import (
        ProteinPairDataset,
        batch_iterator,
        read_pairs_csv,
        sanity_check_batch,
    )
    from protein_interaction.models.cross_set_v1 import CrossSetInteractionV1
    from protein_interaction.utils import ensure_dir, load_config, resolve_path


def parse_args():
    parser = argparse.ArgumentParser(description="Train V1 protein interaction model.")
    parser.add_argument("--config", required=True, help="Path to config file.")
    parser.add_argument("--train_pairs", default="", help="Optional override for train CSV.")
    parser.add_argument("--val_pairs", default="", help="Optional override for val CSV.")
    return parser.parse_args()


def set_tf_seed(seed: int):
    # Support both TF2 and TF1 APIs.
    if hasattr(tf.random, "set_seed"):
        tf.random.set_seed(seed)
    elif hasattr(tf.random, "set_random_seed"):
        tf.random.set_random_seed(seed)
    elif hasattr(tf, "set_random_seed"):
        tf.set_random_seed(seed)
    elif hasattr(tf, "compat") and hasattr(tf.compat, "v1"):
        tf.compat.v1.set_random_seed(seed)


def make_adam_optimizer(lr: float):
    # TF2 uses learning_rate; TF1 tf.keras often expects lr.
    try:
        opt = tf.keras.optimizers.Adam(learning_rate=lr)
    except TypeError:
        opt = tf.keras.optimizers.Adam(lr=lr)
    if hasattr(opt, "apply_gradients"):
        return opt
    # TF1 fallback with explicit apply_gradients support.
    return tf.train.AdamOptimizer(learning_rate=lr)


def make_bce_loss():
    # TF2: class-based BCE. TF1 fallback: backend BCE tensor.
    try:
        return tf.keras.losses.BinaryCrossentropy()
    except AttributeError:
        def _bce(y_true, y_pred):
            return tf.reduce_mean(tf.keras.backend.binary_crossentropy(y_true, y_pred))
        return _bce


def tensor_to_numpy(x):
    if isinstance(x, np.ndarray):
        return x
    if hasattr(x, "numpy"):
        return x.numpy()
    try:
        return tf.keras.backend.get_value(x)
    except Exception:
        return np.asarray(x)


def tensor_to_float(x):
    arr = tensor_to_numpy(x)
    return float(np.asarray(arr).reshape(-1)[0])


def to_model_inputs(batch: Dict):
    return {
        "query_desc": tf.convert_to_tensor(batch["query_desc"], dtype=tf.float32),
        "query_xyz": tf.convert_to_tensor(batch["query_xyz"], dtype=tf.float32),
        "query_mask": tf.convert_to_tensor(batch["query_mask"], dtype=tf.float32),
        "matched_desc": tf.convert_to_tensor(batch["matched_desc"], dtype=tf.float32),
        "matched_xyz": tf.convert_to_tensor(batch["matched_xyz"], dtype=tf.float32),
        "matched_mask": tf.convert_to_tensor(batch["matched_mask"], dtype=tf.float32),
    }


def compute_metrics(y_true: np.ndarray, y_prob: np.ndarray) -> Dict:
    if len(np.unique(y_true)) < 2:
        return {"roc_auc": float("nan"), "pr_auc": float("nan")}
    return {
        "roc_auc": float(roc_auc_score(y_true, y_prob)),
        "pr_auc": float(average_precision_score(y_true, y_prob)),
    }


def evaluate(model, dataset, batch_size, descriptor_dim):
    y_true = []
    y_prob = []
    losses = []
    bce = make_bce_loss()
    for batch in batch_iterator(
        dataset, batch_size=batch_size, descriptor_dim=descriptor_dim, shuffle=False
    ):
        sanity_check_batch(batch)
        inputs = to_model_inputs(batch)
        labels = tf.convert_to_tensor(batch["labels"], dtype=tf.float32)
        prob = model(inputs, training=False)
        loss = tensor_to_float(bce(labels, prob))
        losses.append(loss)
        prob_np = tensor_to_numpy(prob).reshape(-1)
        y_true.extend(batch["labels"].reshape(-1).tolist())
        y_prob.extend(prob_np.tolist())
    metrics = compute_metrics(np.asarray(y_true), np.asarray(y_prob))
    metrics["loss"] = float(np.mean(losses)) if losses else float("nan")
    return metrics


def main():
    args = parse_args()
    cfg = load_config(args.config)
    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)
    set_tf_seed(seed)

    data_cfg = cfg["data"]
    model_cfg = cfg["model"]
    train_cfg = cfg["train"]

    pairs_dir = resolve_path(data_cfg["pairs_out_dir"])
    train_pairs_path = (
        resolve_path(args.train_pairs) if args.train_pairs else pairs_dir / "pairs_train.csv"
    )
    val_pairs_path = (
        resolve_path(args.val_pairs) if args.val_pairs else pairs_dir / "pairs_val.csv"
    )

    train_records = read_pairs_csv(str(train_pairs_path))
    val_records = read_pairs_csv(str(val_pairs_path))
    if len(train_records) == 0:
        raise ValueError(f"No training pairs found in {train_pairs_path}")

    train_ds = ProteinPairDataset(cfg, train_records, seed=seed)
    val_ds = ProteinPairDataset(cfg, val_records, seed=seed + 1)

    descriptor_dim = int(model_cfg["descriptor_dim"])
    batch_size = int(train_cfg["batch_size"])
    epochs = int(train_cfg["epochs"])
    lr = float(train_cfg["learning_rate"])
    patience = int(train_cfg.get("early_stop_patience", 6))

    model = CrossSetInteractionV1(
        hidden_dim=int(model_cfg["encoder_hidden_dim"]),
        classifier_hidden_dim=int(model_cfg["classifier_hidden_dim"]),
        dropout=float(model_cfg.get("dropout", 0.2)),
    )
    optimizer = make_adam_optimizer(lr)
    bce = make_bce_loss()

    ckpt_dir = ensure_dir(train_cfg["checkpoint_dir"])
    log_dir = ensure_dir(train_cfg["log_dir"])
    best_path = ckpt_dir / "best"
    history_path = log_dir / "history.json"

    best_val_auc = -1.0
    no_improve = 0
    history = []

    for epoch in range(1, epochs + 1):
        train_losses = []
        for batch in batch_iterator(
            train_ds, batch_size=batch_size, descriptor_dim=descriptor_dim, shuffle=True
        ):
            sanity_check_batch(batch)
            inputs = to_model_inputs(batch)
            labels = tf.convert_to_tensor(batch["labels"], dtype=tf.float32)
            with tf.GradientTape() as tape:
                prob = model(inputs, training=True)
                loss = bce(labels, prob)
            grads = tape.gradient(loss, model.trainable_variables)
            optimizer.apply_gradients(zip(grads, model.trainable_variables))
            train_losses.append(tensor_to_float(loss))

        train_loss = float(np.mean(train_losses)) if train_losses else float("nan")
        val_metrics = evaluate(model, val_ds, batch_size=batch_size, descriptor_dim=descriptor_dim)
        record = {
            "epoch": epoch,
            "train_loss": train_loss,
            "val_loss": val_metrics["loss"],
            "val_roc_auc": val_metrics["roc_auc"],
            "val_pr_auc": val_metrics["pr_auc"],
        }
        history.append(record)
        print(json.dumps(record, sort_keys=True))

        val_auc = val_metrics["roc_auc"]
        if np.isnan(val_auc):
            val_auc = -1.0
        if val_auc > best_val_auc:
            best_val_auc = val_auc
            no_improve = 0
            model.save_weights(str(best_path), save_format="tf")
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"Early stopping at epoch {epoch}.")
                break

    with open(history_path, "w") as f:
        json.dump(
            {
                "config_path": str(Path(args.config)),
                "best_val_roc_auc": best_val_auc,
                "best_checkpoint": str(best_path),
                "history": history,
            },
            f,
            indent=2,
        )
    print(f"Saved training history to {history_path}")
    print(f"Best checkpoint: {best_path}")


if __name__ == "__main__":
    main()

