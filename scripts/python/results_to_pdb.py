#!/usr/bin/env python3
"""
results_to_pdb.py — Build assembled target + transformed matched-protein PDBs.

For each row in a Neosurf results CSV, load the target and matched-protein
structures (using paths from the preprocess manifest), keep the relevant PDB
chains, apply the row's flattened 4x4 affine transform to the matched protein,
merge the transformed matched protein into the target structure, strip hydrogens
and water, and write one output PDB per row.

This mirrors the PyMOL command sequence in notebooks/scratch.ipynb but uses
BioPython instead of PyMOL.

Usage:
    python scripts/python/results_to_pdb.py \\
        --results_csv df_results_dedup_picked.csv \\
        --preprocess_metadata_csv data/processing/2_masif_preprocess/df_preprocess_ok.csv \\
        --out_dir assembled_pdbs
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.StructureBuilder import StructureBuilder
from tqdm import tqdm

WATER_RESNAMES = frozenset({"HOH", "WAT", "H2O", "DOD"})
CHAIN_ID_LETTERS = [chr(c) for c in range(ord("A"), ord("Z") + 1)]


def build_target_lookup(df_preprocess_ok: pd.DataFrame) -> dict:
    """Map target id -> manifest fields needed for structure assembly."""
    return (
        df_preprocess_ok.set_index("target")[
            ["pdb_id", "pdb_protein_chain", "pdb_ligand_chain", "pdb_path"]
        ]
        .to_dict("index")
    )


def get_target_info(target_id: str, target_lookup: dict) -> tuple[str, list[str], str | None]:
    """Return (pdb_id, chain_list, pdb_path) for a target id.

    Falls back to splitting on the last underscore for targets not in the
    manifest (e.g. manually specified query targets).
    """
    if target_id in target_lookup:
        info = target_lookup[target_id]
        chains = [info["pdb_protein_chain"]]
        if info["pdb_ligand_chain"] != info["pdb_protein_chain"]:
            chains.append(info["pdb_ligand_chain"])
        return info["pdb_id"], chains, info.get("pdb_path")

    parts = target_id.rsplit("_", 1)
    chains = list(parts[1]) if len(parts) > 1 else []
    return parts[0], chains, None


def apply_transform(structure, transform_str: str) -> None:
    """Apply a flattened 4x4 affine transform to every atom in a structure."""
    transform_flat = np.fromstring(str(transform_str), sep=",")
    if transform_flat.size != 16:
        raise ValueError(f"Expected flattened transform of size 16, got {transform_flat.size}")
    matrix = transform_flat.reshape((4, 4))
    rotation = matrix[:3, :3]
    translation = matrix[:3, 3]
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.set_coord(np.dot(rotation, atom.get_coord()) + translation)


def extract_chains(structure, chain_ids: list[str]):
    """Return a new structure containing only the requested chains."""
    builder = StructureBuilder()
    builder.init_structure("extracted")
    builder.init_model(0)
    for chain_id in chain_ids:
        builder.init_chain(chain_id)

    out = builder.get_structure()
    out_model = out[0]

    for model in structure:
        for chain_id in chain_ids:
            if chain_id not in model:
                continue
            for residue in model[chain_id]:
                out_model[chain_id].add(residue.copy())

    return out


def remove_hydrogens_and_water(structure) -> None:
    """Remove hydrogen atoms and water residues (in-place)."""
    for model in structure:
        for chain in list(model):
            for residue in list(chain):
                resname = residue.get_resname().strip()
                if resname in WATER_RESNAMES:
                    chain.detach_child(residue.id)
                    continue
                for atom in list(residue):
                    element = (atom.element or "").strip().upper()
                    name = atom.get_name().strip()
                    if element == "H" or (not element and name.upper().startswith("H")):
                        residue.detach_child(atom.get_id())


def _next_available_chain_id(used_ids: set[str]) -> str:
    """Return the next unused single-letter chain ID."""
    for chain_id in CHAIN_ID_LETTERS:
        if chain_id not in used_ids:
            return chain_id
    raise ValueError("No available single-letter chain IDs for matched protein")


def merge_structures(target_structure, source_structure) -> None:
    """Append matched-protein chains into target, renaming on chain-id conflict."""
    target_model = target_structure[0]
    source_model = source_structure[0]
    used_chain_ids = {chain.id for chain in target_model}

    for src_chain in source_model:
        chain_copy = copy.deepcopy(src_chain)
        if chain_copy.id in used_chain_ids:
            chain_copy.id = _next_available_chain_id(used_chain_ids)
        used_chain_ids.add(chain_copy.id)
        target_model.add(chain_copy)


def output_pdb_name(row: pd.Series) -> str:
    """Build output filename for one results row."""
    return (
        f"{row['target']}_{row['target_vix']}_to_"
        f"{row['matched_protein']}_{row['matched_patch_id']}.pdb"
    )


def build_assembled_structure(
    row: pd.Series,
    target_lookup: dict,
    parser: PDBParser,
):
    """Load, transform, and merge target + matched protein for one results row."""
    target_id = row["target"]
    matched_id = row["matched_protein"]
    transform = row["flattened_transform"]

    _, target_chains, target_path = get_target_info(target_id, target_lookup)
    _, matched_chains, matched_path = get_target_info(matched_id, target_lookup)

    if not target_path or not Path(target_path).is_file():
        raise FileNotFoundError(f"Target PDB not found for '{target_id}': {target_path}")
    if not matched_path or not Path(matched_path).is_file():
        raise FileNotFoundError(f"Matched protein PDB not found for '{matched_id}': {matched_path}")

    target_struct = parser.get_structure("target", target_path)
    matched_struct = parser.get_structure("matched", matched_path)

    target_struct = extract_chains(target_struct, target_chains)
    matched_struct = extract_chains(matched_struct, matched_chains)
    apply_transform(matched_struct, transform)
    merge_structures(target_struct, matched_struct)
    remove_hydrogens_and_water(target_struct)

    return target_struct


def write_structure(structure, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write assembled target + transformed matched-protein PDB files."
    )
    parser.add_argument("--results_csv", required=True, help="Neosurf results CSV.")
    parser.add_argument(
        "--preprocess_metadata_csv",
        required=True,
        help="Preprocess manifest CSV (df_preprocess_ok.csv).",
    )
    parser.add_argument("--out_dir", required=True, help="Output directory for PDB files.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    results_csv = Path(args.results_csv)
    preprocess_csv = Path(args.preprocess_metadata_csv)
    out_dir = Path(args.out_dir)

    if not results_csv.is_file():
        print(f"Results CSV not found: {results_csv}", file=sys.stderr)
        return 1
    if not preprocess_csv.is_file():
        print(f"Preprocess metadata CSV not found: {preprocess_csv}", file=sys.stderr)
        return 1

    df_results = pd.read_csv(results_csv)
    df_preprocess_ok = pd.read_csv(preprocess_csv)
    target_lookup = build_target_lookup(df_preprocess_ok)

    required_cols = {
        "target",
        "target_vix",
        "matched_protein",
        "matched_patch_id",
        "flattened_transform",
    }
    missing = required_cols - set(df_results.columns)
    if missing:
        print(f"Results CSV missing required columns: {sorted(missing)}", file=sys.stderr)
        return 1

    pdb_parser = PDBParser(QUIET=True)
    failures: list[tuple[str, str]] = []

    for _, row in tqdm(df_results.iterrows(), total=len(df_results), desc="Writing PDBs"):
        out_name = output_pdb_name(row)
        out_path = out_dir / out_name
        try:
            structure = build_assembled_structure(row, target_lookup, pdb_parser)
            write_structure(structure, out_path)
        except Exception as exc:
            failures.append((out_name, str(exc)))

    if failures:
        print(f"\nFailed to write {len(failures)} / {len(df_results)} PDB files:", file=sys.stderr)
        for name, err in failures:
            print(f"  {name}: {err}", file=sys.stderr)
        return 1

    print(f"Wrote {len(df_results)} PDB files to {out_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
