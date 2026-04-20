"""Entry point to train MaSIF-Coherence v1.

Mirrors the call pattern of `masif/source/masif_ppi_search/masif_ppi_search_train.py`:

    python train.py <custom_params_module>

where `<custom_params_module>` is an importable Python module exposing a
`custom_params` dict (e.g. `nn_models.v1.custom_params`). The defaults in
`masif_coherence.config` are overlaid with that dict.

Positives come from `params["training_list"]`; one random non-cognate negative
per positive is generated. Validation AUROC is printed at the end of every
epoch and a checkpoint is saved to `params["model_dir"]`.
"""

from __future__ import annotations

import importlib
import os
import random
import sys
import time
from typing import List

import numpy as np
import torch
from sklearn.metrics import roc_auc_score

# Make the module importable whether it is invoked as a script or as a module.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

from config import get_default_params, merge_custom_params  # noqa: E402
from data import CoherenceDataset, read_pair_list  # noqa: E402
from model import build_model_from_params  # noqa: E402


def _resolve_device(pref: str) -> torch.device:
    if pref == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(pref)


def _set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def _split_train_val(pair_ids: List[str], val_fraction: float, seed: int):
    rng = random.Random(seed)
    shuffled = list(pair_ids)
    rng.shuffle(shuffled)
    n_val = max(1, int(round(val_fraction * len(shuffled)))) if len(shuffled) > 1 else 0
    return shuffled[n_val:], shuffled[:n_val]


def _evaluate(model, dataset, device) -> float:
    model.eval()
    y_true: List[float] = []
    y_pred: List[float] = []
    with torch.no_grad():
        for i in range(len(dataset)):
            prot_A, prot_B, label = dataset[i]
            X_A = prot_A.X.to(device)
            D_A = prot_A.D.to(device)
            X_B = prot_B.X.to(device)
            D_B = prot_B.D.to(device)
            logit = model(X_A, D_A, X_B, D_B)
            y_true.append(float(label.item()))
            y_pred.append(float(torch.sigmoid(logit).item()))
    model.train()
    if len(set(y_true)) < 2:
        return float("nan")
    return float(roc_auc_score(y_true, y_pred))


def main(custom_params_module: str) -> None:
    params = get_default_params()
    custom_mod = importlib.import_module(custom_params_module, package=None)
    params = merge_custom_params(params, custom_mod.custom_params)

    _set_seed(params["seed"])
    device = _resolve_device(params["device"])
    print(f"Using device: {device}")

    training_ids = read_pair_list(params["training_list"])
    train_ids, val_ids = _split_train_val(
        training_ids, params["val_fraction"], params["seed"]
    )
    print(f"Train pairs: {len(train_ids)}  Val pairs: {len(val_ids)}")

    train_ds = CoherenceDataset(
        pair_ids=train_ids,
        precomputation_dir=params["masif_precomputation_dir"],
        ply_chain_dir=params["ply_chain_dir"],
        neg_per_pos=params["neg_per_pos"],
        seed=params["seed"],
    )
    val_ds = CoherenceDataset(
        pair_ids=val_ids,
        precomputation_dir=params["masif_precomputation_dir"],
        ply_chain_dir=params["ply_chain_dir"],
        neg_per_pos=params["neg_per_pos"],
        seed=params["seed"] + 1,
    )

    model = build_model_from_params(params).to(device)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Model parameters: {n_params:,}")

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=params["learning_rate"],
        weight_decay=params["weight_decay"],
    )
    loss_fn = torch.nn.BCEWithLogitsLoss()

    os.makedirs(params["model_dir"], exist_ok=True)

    grad_accum = max(1, int(params["grad_accum_steps"]))
    log_every = max(1, int(params["log_every"]))

    for epoch in range(int(params["num_epochs"])):
        t0 = time.time()
        order = list(range(len(train_ds)))
        random.shuffle(order)

        optimizer.zero_grad(set_to_none=True)
        running_loss = 0.0
        running_count = 0
        micro_step = 0

        for step, idx in enumerate(order):
            prot_A, prot_B, label = train_ds[idx]
            X_A = prot_A.X.to(device)
            D_A = prot_A.D.to(device)
            X_B = prot_B.X.to(device)
            D_B = prot_B.D.to(device)
            label_t = label.to(device)

            logit = model(X_A, D_A, X_B, D_B)
            loss = loss_fn(logit, label_t) / grad_accum
            loss.backward()

            running_loss += float(loss.item()) * grad_accum
            running_count += 1
            micro_step += 1

            if micro_step % grad_accum == 0:
                optimizer.step()
                optimizer.zero_grad(set_to_none=True)

            if (step + 1) % log_every == 0:
                avg = running_loss / max(1, running_count)
                print(f"  epoch {epoch} step {step + 1}/{len(order)} loss {avg:.4f}")

        # Flush any remaining micro-batch.
        if micro_step % grad_accum != 0:
            optimizer.step()
            optimizer.zero_grad(set_to_none=True)

        train_loss = running_loss / max(1, running_count)
        val_auroc = _evaluate(model, val_ds, device) if len(val_ds) > 0 else float("nan")
        dt = time.time() - t0
        print(
            f"epoch {epoch} done: train_loss {train_loss:.4f} "
            f"val_auroc {val_auroc:.4f} time {dt:.1f}s"
        )

        if params["save_every_epoch"]:
            ckpt_path = os.path.join(params["model_dir"], f"epoch_{epoch:03d}.pt")
            torch.save(
                {
                    "epoch": epoch,
                    "model_state": model.state_dict(),
                    "optimizer_state": optimizer.state_dict(),
                    "params": params,
                    "val_auroc": val_auroc,
                },
                ckpt_path,
            )
            print(f"  saved {ckpt_path}")

    final_path = os.path.join(params["model_dir"], "final.pt")
    torch.save({"model_state": model.state_dict(), "params": params}, final_path)
    print(f"Wrote final checkpoint to {final_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python train.py <custom_params_module>")
        sys.exit(1)
    main(sys.argv[1])
