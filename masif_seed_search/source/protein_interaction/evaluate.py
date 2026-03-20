#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf
from sklearn.metrics import average_precision_score, roc_auc_score

try:
    from protein_interaction.dataset import ProteinPairDataset, batch_iterator, read_pairs_csv, sanity_check_batch
    from protein_interaction.models.cross_set_v1 import CrossSetInteractionV1
    from protein_interaction.utils import ensure_dir, load_config, resolve_path
except ModuleNotFoundError:
    # Allow direct execution: python path/to/evaluate.py ...
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from protein_interaction.dataset import ProteinPairDataset, batch_iterator, read_pairs_csv, sanity_check_batch
    from protein_interaction.models.cross_set_v1 import CrossSetInteractionV1
    from protein_interaction.utils import ensure_dir, load_config, resolve_path


def parse_args():
    parser = argparse.ArgumentParser(description="Evaluate V1 protein interaction model.")
    parser.add_argument("--config", required=True, help="Path to config.")
    parser.add_argument("--pairs_csv", default="", help="Optional override split CSV.")
    parser.add_argument("--checkpoint", default="", help="Optional checkpoint override.")
    parser.add_argument("--split", default="test", choices=["train", "val", "test"])
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


def to_model_inputs(batch):
    return {
        "query_desc": tf.convert_to_tensor(batch["query_desc"], dtype=tf.float32),
        "query_xyz": tf.convert_to_tensor(batch["query_xyz"], dtype=tf.float32),
        "query_mask": tf.convert_to_tensor(batch["query_mask"], dtype=tf.float32),
        "matched_desc": tf.convert_to_tensor(batch["matched_desc"], dtype=tf.float32),
        "matched_xyz": tf.convert_to_tensor(batch["matched_xyz"], dtype=tf.float32),
        "matched_mask": tf.convert_to_tensor(batch["matched_mask"], dtype=tf.float32),
    }


def checkpoint_exists(path_like: str) -> bool:
    path = Path(path_like)
    return path.exists() or Path(str(path) + ".index").exists()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)
    set_tf_seed(seed)

    data_cfg = cfg["data"]
    model_cfg = cfg["model"]
    train_cfg = cfg["train"]
    eval_cfg = cfg["eval"]

    pairs_dir = resolve_path(data_cfg["pairs_out_dir"])
    pairs_path = resolve_path(args.pairs_csv) if args.pairs_csv else (pairs_dir / f"pairs_{args.split}.csv")
    records = read_pairs_csv(str(pairs_path))
    dataset = ProteinPairDataset(cfg, records, seed=seed + 11)

    descriptor_dim = int(model_cfg["descriptor_dim"])
    model = CrossSetInteractionV1(
        hidden_dim=int(model_cfg["encoder_hidden_dim"]),
        classifier_hidden_dim=int(model_cfg["classifier_hidden_dim"]),
        dropout=float(model_cfg.get("dropout", 0.2)),
    )

    ckpt_path = args.checkpoint
    if not ckpt_path:
        ckpt_path = eval_cfg.get("checkpoint_path", "")
    if not ckpt_path:
        ckpt_dir = resolve_path(train_cfg["checkpoint_dir"])
        ckpt_best_tf = str(ckpt_dir / "best")
        ckpt_best_h5 = str(ckpt_dir / "best.weights.h5")
        if checkpoint_exists(ckpt_best_tf):
            ckpt_path = ckpt_best_tf
        elif checkpoint_exists(ckpt_best_h5):
            ckpt_path = ckpt_best_h5
        else:
            ckpt_path = ckpt_best_tf
    ckpt_path = str(resolve_path(ckpt_path))

    # Build model once to initialize variables, then load weights.
    first_batch = next(batch_iterator(dataset, batch_size=1, descriptor_dim=descriptor_dim, shuffle=False))
    _ = model(to_model_inputs(first_batch), training=False)
    model.load_weights(ckpt_path)

    y_true = []
    y_prob = []
    pair_ids = []
    bce = make_bce_loss()
    losses = []
    for batch in batch_iterator(
        dataset,
        batch_size=int(train_cfg["batch_size"]),
        descriptor_dim=descriptor_dim,
        shuffle=False,
    ):
        sanity_check_batch(batch)
        probs = tensor_to_numpy(model(to_model_inputs(batch), training=False)).reshape(-1)
        labels = batch["labels"].reshape(-1)
        labels_tf = tf.convert_to_tensor(labels.reshape(-1, 1), dtype=tf.float32)
        probs_tf = tf.convert_to_tensor(probs.reshape(-1, 1), dtype=tf.float32)
        loss = tensor_to_float(bce(labels_tf, probs_tf))
        losses.append(loss)
        y_true.extend(labels.tolist())
        y_prob.extend(probs.tolist())
        pair_ids.extend(batch["pair_ids"])

    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    metrics = {
        "split": args.split,
        "num_samples": int(len(y_true)),
        "loss": float(np.mean(losses)) if losses else float("nan"),
    }
    if len(np.unique(y_true)) >= 2:
        metrics["roc_auc"] = float(roc_auc_score(y_true, y_prob))
        metrics["pr_auc"] = float(average_precision_score(y_true, y_prob))
    else:
        metrics["roc_auc"] = float("nan")
        metrics["pr_auc"] = float("nan")

    out_dir = ensure_dir(eval_cfg["output_dir"])
    metrics_path = out_dir / f"metrics_{args.split}.json"
    preds_path = out_dir / f"predictions_{args.split}.csv"

    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2, sort_keys=True)
    with open(preds_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["query_id", "matched_id", "label", "probability"])
        for (query_id, matched_id), label, prob in zip(pair_ids, y_true, y_prob):
            writer.writerow([query_id, matched_id, int(label), float(prob)])

    print(json.dumps(metrics, indent=2, sort_keys=True))
    print(f"Saved metrics to {metrics_path}")
    print(f"Saved predictions to {preds_path}")


if __name__ == "__main__":
    main()

