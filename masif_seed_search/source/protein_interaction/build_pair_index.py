#!/usr/bin/env python3
import argparse
import csv
import json
import random
import sys
from pathlib import Path
from typing import Dict, List, Tuple

try:
    from protein_interaction.utils import ensure_dir, load_config, read_id_list, resolve_path
except ModuleNotFoundError:
    # Allow direct execution: python path/to/build_pair_index.py ...
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from protein_interaction.utils import ensure_dir, load_config, read_id_list, resolve_path



def parse_args():
    parser = argparse.ArgumentParser(description="Build protein interaction pair CSV splits.")
    parser.add_argument("--config", required=True, help="Path to YAML/JSON config file.")
    parser.add_argument(
        "--tiny",
        action="store_true",
        help="Enable tiny subset mode regardless of config.",
    )
    return parser.parse_args()


def available_ppi_ids(cfg: Dict) -> List[str]:
    data_cfg = cfg["data"]
    desc_dir = resolve_path(data_cfg["descriptors_dir"])
    query_file = data_cfg["query_desc_file"]
    matched_file = data_cfg["matched_desc_file"]
    ids: List[str] = []
    for path in sorted(desc_dir.iterdir()):
        if not path.is_dir():
            continue
        if (path / query_file).exists() and (path / matched_file).exists():
            ids.append(path.name)
    return ids


def split_ids(cfg: Dict, ids: List[str]) -> Dict[str, List[str]]:
    data_cfg = cfg["data"]
    rng = random.Random(cfg.get("seed", 42))
    train_ids = [x for x in read_id_list(data_cfg["train_list"]) if x in ids]
    test_ids = [x for x in read_id_list(data_cfg["test_list"]) if x in ids]

    rng.shuffle(train_ids)
    n_val = max(1, int(len(train_ids) * float(data_cfg.get("val_fraction_from_train", 0.1))))
    val_ids = train_ids[:n_val]
    train_ids = train_ids[n_val:]

    tiny_cfg = data_cfg.get("tiny_subset", {})
    tiny_enabled = bool(tiny_cfg.get("enabled", False))
    if tiny_enabled:
        max_pos = int(tiny_cfg.get("max_pos_per_split", 20))
        train_ids = train_ids[:max_pos]
        val_ids = val_ids[:max_pos]
        test_ids = test_ids[:max_pos]

    return {"train": train_ids, "val": val_ids, "test": test_ids}


def make_pairs_for_split(
    split_name: str,
    split_ids_list: List[str],
    negatives_per_positive: int,
    rng: random.Random,
) -> List[Tuple[str, str, int, str, str, str]]:
    rows: List[Tuple[str, str, int, str, str, str]] = []
    if not split_ids_list:
        return rows

    for ppi_id in split_ids_list:
        rows.append((ppi_id, ppi_id, 1, split_name, "", "native_pair"))

        neg_pool = [x for x in split_ids_list if x != ppi_id]
        if not neg_pool:
            continue
        n_neg = min(len(neg_pool), negatives_per_positive)
        sampled = rng.sample(neg_pool, k=n_neg)
        for neg_id in sampled:
            rows.append((ppi_id, neg_id, 0, split_name, "random_cross", "synthetic_negative"))
    return rows


def write_split_csv(out_dir: Path, split_name: str, rows: List[Tuple[str, str, int, str, str, str]]):
    outfile = out_dir / f"pairs_{split_name}.csv"
    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["query_id", "matched_id", "label", "split", "negative_type", "source_metadata"]
        )
        writer.writerows(rows)


def summarize_rows(rows: List[Tuple[str, str, int, str, str, str]]) -> Dict:
    positives = sum(1 for row in rows if row[2] == 1)
    negatives = sum(1 for row in rows if row[2] == 0)
    return {"num_rows": len(rows), "num_positive": positives, "num_negative": negatives}


def main():
    args = parse_args()
    cfg = load_config(args.config)
    data_cfg = cfg["data"]

    if args.tiny:
        data_cfg.setdefault("tiny_subset", {})
        data_cfg["tiny_subset"]["enabled"] = True

    ids = available_ppi_ids(cfg)
    splits = split_ids(cfg, ids)
    neg_per_pos = int(data_cfg.get("negatives_per_positive", 3))
    rng = random.Random(cfg.get("seed", 42))

    out_dir = ensure_dir(data_cfg["pairs_out_dir"])
    summary = {"available_ppi_ids": len(ids), "splits": {}}

    for split_name in ("train", "val", "test"):
        rows = make_pairs_for_split(split_name, splits[split_name], neg_per_pos, rng)
        write_split_csv(out_dir, split_name, rows)
        summary["splits"][split_name] = {
            "num_ppi_ids": len(splits[split_name]),
            **summarize_rows(rows),
        }

    with open(out_dir / "pairs_summary.json", "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

