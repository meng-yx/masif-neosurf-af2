#!/usr/bin/env python3
"""
prepare_input.py — CLI for preparing MaSIF input structures.

For each row in --input_csv (columns: pdb_id, protein_chain, ligand_chain,
ligand_code), this script:

  1. Downloads the full structure (PDB or mmCIF) from files.rcsb.org and the
     ligand SDF from models.rcsb.org (with connectivity/bond orders).
  2. Trims to the target protein chain (standard amino acids) plus one ligand
     residue on the ligand chain.
  3. Runs EvoEF2 RepairStructure on the protein-only PDB.
  4. Merges ligand coordinates from the downloaded SDF onto the repaired protein
     as HETATM records using a 3-letter ligand code (ligand_code truncated to
     3 characters).
  5. Appends four manifest columns to every input row:
       pdb_path, target, ligand, ligand_path
     Manifest ligand/ligand_path use the 3-letter code; the input ligand_code
     column is preserved unchanged.
     (empty when a row fails) and writes the result to --out_csv.

Usage:
    python scripts/python/prepare_input.py \
        --input_csv data/processing/1_prep_input/input_subsets/input_1.csv \
        --outdir    data/input \
        --out_csv   data/processing/1_prep_input/input_1_manifest.csv \
        --evoef2_bin EvoEF2/EvoEF2 \
        [--overwrite]

Exit code is non-zero when at least one row failed.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.MMCIFParser import MMCIFParser
from rdkit import Chem
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",
}

REQUIRED_COLUMNS = {"pdb_id", "protein_chain", "ligand_chain", "ligand_code"}

# ---------------------------------------------------------------------------
# BioPython selectors
# ---------------------------------------------------------------------------


def _is_standard_protein_residue(residue):
    return residue.id[0] == " " and residue.get_resname().strip() in STANDARD_AA


class _ProteinOnlySelect(Select):
    def __init__(self, protein_chain):
        self.protein_chain = protein_chain

    def accept_chain(self, chain):
        return chain.id == self.protein_chain

    def accept_residue(self, residue):
        return _is_standard_protein_residue(residue)


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def chain_suffix(protein_chain, ligand_chain):
    """Return deduplicated chain suffix: C+A -> CA, A+A -> A."""
    if protein_chain == ligand_chain:
        return protein_chain
    return protein_chain + ligand_chain


def ligand_code_tla(ligand_code):
    """Return 3-letter ligand code for PDB/SDF filenames and manifest columns."""
    return ligand_code.strip()[:3]


def _pdb_chain_ids(protein_chain, ligand_chain):
    """Map input chain IDs to single-character IDs required by PDB format."""
    protein_chain = str(protein_chain).strip()
    ligand_chain = str(ligand_chain).strip()
    used = set()

    def remap_one(chain_id):
        if len(chain_id) == 1:
            used.add(chain_id)
            return chain_id
        candidate = chain_id[0]
        if candidate not in used:
            used.add(candidate)
            return candidate
        for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            if letter not in used:
                used.add(letter)
                return letter
        raise ValueError(
            f"Cannot remap chains {protein_chain!r}, {ligand_chain!r} to PDB format"
        )

    pdb_protein = remap_one(protein_chain)
    if ligand_chain == protein_chain:
        return pdb_protein, pdb_protein
    return pdb_protein, remap_one(ligand_chain)


def _remap_structure_chains(structure, chain_map):
    """Rename chains in-place so PDB writers only see single-character IDs."""
    for model in structure:
        for chain in model:
            new_id = chain_map.get(chain.id)
            if new_id is not None:
                chain.id = new_id


def expected_output_paths(pdb_id, protein_chain, ligand_chain, ligand_code, outdir):
    """Return the manifest dict without touching the filesystem."""
    outdir = Path(outdir)
    ligand_tla = ligand_code_tla(ligand_code)
    chains = chain_suffix(protein_chain, ligand_chain)
    target = f"{pdb_id}_{chains}"
    ligand_name = f"{ligand_tla}_{ligand_chain}"
    pdb_path = outdir / f"{target}.pdb"
    ligand_path = outdir / f"{pdb_id}_{protein_chain}_{ligand_tla}.sdf"
    return {
        "pdb_path": str(pdb_path.resolve()),
        "target": target,
        "ligand": ligand_name,
        "ligand_path": str(ligand_path.resolve()),
    }


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------


RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core"
RCSB_MODEL_SERVER_API = "https://models.rcsb.org/v1"


def _fetch_json(url):
    with urlopen(url) as response:
        return json.load(response)


def _try_download_url(url, dest_path):
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urlopen(url) as response:
            dest_path.write_bytes(response.read())
        return True
    except (HTTPError, URLError):
        return False


def resolve_ligand_model_server_ids(pdb_id, auth_asym_id, ligand_code):
    """
    Map author chain/residue name to model-server query parameters.

    Input chain IDs are author-assigned (auth_asym_id). The model server expects
    label_asym_id plus auth_seq_id from the deposited mmCIF/PDBx tables.
    """
    pdb_id = pdb_id.strip().upper()
    auth_asym_id = str(auth_asym_id).strip()
    ligand_code = ligand_code.strip().upper()
    ligand_tla = ligand_code_tla(ligand_code)

    entry = _fetch_json(f"{RCSB_DATA_API}/entry/{pdb_id}")
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "non_polymer_entity_ids", []
    )
    if not entity_ids:
        raise ValueError(f"No non-polymer entities found for {pdb_id}")

    label_asym_id = None
    for entity_id in entity_ids:
        entity = _fetch_json(f"{RCSB_DATA_API}/nonpolymer_entity/{pdb_id}/{entity_id}")
        identifiers = entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
        comp_id = identifiers.get("nonpolymer_comp_id", "").strip().upper()
        if comp_id not in {ligand_code, ligand_tla}:
            continue

        asym_ids = identifiers.get("asym_ids", [])
        auth_asym_ids = identifiers.get("auth_asym_ids", [])
        for label_id, auth_id in zip(asym_ids, auth_asym_ids):
            if str(auth_id).strip() == auth_asym_id:
                label_asym_id = str(label_id).strip()
                break
        if label_asym_id is not None:
            break

    if label_asym_id is None:
        raise ValueError(
            f"No non-polymer entity for {pdb_id} with comp_id={ligand_code} "
            f"on auth chain {auth_asym_id}"
        )

    instance = _fetch_json(
        f"{RCSB_DATA_API}/nonpolymer_entity_instance/{pdb_id}/{label_asym_id}"
    )
    instance_ids = instance.get("rcsb_nonpolymer_entity_instance_container_identifiers", {})
    auth_seq_id = str(instance_ids.get("auth_seq_id", "")).strip()
    if not auth_seq_id:
        raise ValueError(
            f"Could not resolve auth_seq_id for {pdb_id} ligand {ligand_code} "
            f"(label_asym_id={label_asym_id})"
        )

    return label_asym_id, auth_seq_id


def build_ligand_sdf_url(pdb_id, auth_asym_id, ligand_code):
    """Build models.rcsb.org ligand SDF download URL."""
    pdb_id_lower = pdb_id.strip().lower()
    ligand_tla = ligand_code_tla(ligand_code)
    label_asym_id, auth_seq_id = resolve_ligand_model_server_ids(
        pdb_id, auth_asym_id, ligand_code
    )
    filename = f"{pdb_id_lower}_{label_asym_id}_{ligand_tla}.sdf"
    query = urlencode(
        {
            "auth_seq_id": auth_seq_id,
            "label_asym_id": label_asym_id,
            "encoding": "sdf",
            "filename": filename,
        }
    )
    return f"{RCSB_MODEL_SERVER_API}/{pdb_id_lower}/ligand?{query}"


def download_ligand_sdf(pdb_id, auth_asym_id, ligand_code, dest_path):
    """Download ligand SDF with connectivity from models.rcsb.org."""
    dest_path = Path(dest_path)
    url = build_ligand_sdf_url(pdb_id, auth_asym_id, ligand_code)
    if not _try_download_url(url, dest_path):
        raise ValueError(f"Could not download ligand SDF from {url}")
    mol = Chem.SDMolSupplier(str(dest_path), sanitize=False, removeHs=False)[0]
    if mol is None:
        dest_path.unlink(missing_ok=True)
        raise ValueError(f"Downloaded ligand SDF is unreadable: {url}")


def download_structure(pdb_id, dest_dir):
    """Download structure from RCSB; try PDB first, then mmCIF. Return path or None."""
    dest_dir = Path(dest_dir)
    pdb_path = dest_dir / f"{pdb_id}.pdb"
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    if _try_download_url(pdb_url, pdb_path):
        return pdb_path
    cif_path = dest_dir / f"{pdb_id}.cif"
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    if _try_download_url(cif_url, cif_path):
        return cif_path
    warnings.warn(
        f"Could not download structure for {pdb_id} (tried {pdb_url} and {cif_url}); skipping.",
        UserWarning,
        stacklevel=2,
    )
    return None


# ---------------------------------------------------------------------------
# Structure helpers
# ---------------------------------------------------------------------------


def _load_structure(pdb_id, structure_path):
    structure_path = Path(structure_path)
    if structure_path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, str(structure_path))


def _format_hetatm_line(serial, atom_name, resname, chain_id, resseq, x, y, z, element):
    return (
        f"HETATM{serial:5d} {atom_name:>4s} {resname:>3s} {chain_id:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}  \n"
    )


def _sdf_to_hetatm_lines(sdf_path, resname, chain_id, resseq=1):
    mol = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not read ligand SDF: {sdf_path}")
    conf = mol.GetConformer()
    element_counts = {}
    lines = []
    for i, atom in enumerate(mol.GetAtoms()):
        element = atom.GetSymbol()
        element_counts[element] = element_counts.get(element, 0) + 1
        atom_name = f"{element}{element_counts[element]:02d}"[:4]
        pos = conf.GetAtomPosition(i)
        lines.append(
            _format_hetatm_line(
                i + 1,
                atom_name,
                resname,
                chain_id,
                resseq,
                pos.x,
                pos.y,
                pos.z,
                element,
            )
        )
    if not lines:
        raise ValueError(f"No atoms found in ligand SDF: {sdf_path}")
    return lines


def _merge_repaired_protein_with_ligand_sdf(
    repaired_pdb, ligand_sdf, ligand_chain, ligand_tla, output_pdb, resseq=1
):
    ligand_lines = _sdf_to_hetatm_lines(ligand_sdf, ligand_tla, ligand_chain, resseq)
    out_lines = []
    with open(repaired_pdb) as handle:
        for line in handle:
            if line.startswith("END"):
                break
            if line.strip():
                out_lines.append(line)
    out_lines.extend(ligand_lines)
    out_lines.append("END\n")
    Path(output_pdb).write_text("".join(out_lines))


def evoef2_repair_structure(input_pdb, output_pdb, evoef2_bin):
    """Run EvoEF2 RepairStructure; intermediate {stem}_Repair.pdb lives in a temp dir."""
    evoef2_bin = str(evoef2_bin)
    if not os.path.isfile(evoef2_bin):
        raise FileNotFoundError(
            f"EvoEF2 binary not found at {evoef2_bin}. Build with: cd EvoEF2 && ./build.sh"
        )
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        stem = "structure"
        work_pdb = tmp / f"{stem}.pdb"
        shutil.copy2(input_pdb, work_pdb)
        subprocess.run(
            [evoef2_bin, "--command=RepairStructure", f"--pdb={stem}.pdb"],
            cwd=tmp,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        repaired = tmp / f"{stem}_Repair.pdb"
        if not repaired.is_file():
            raise FileNotFoundError(f"EvoEF2 did not create {repaired}")
        shutil.copy2(repaired, output_pdb)


# ---------------------------------------------------------------------------
# Core orchestration
# ---------------------------------------------------------------------------


def prepare_input_structures(
    pdb_id, protein_chain, ligand_chain, ligand_code, outdir, evoef2_bin, overwrite=False
):
    """
    Download, trim, EvoEF2-repair, and merge ligand for one complex.

    Returns a manifest dict on success, or None on failure.
    When overwrite=False and both output files already exist, returns the
    expected dict immediately without downloading or reprocessing.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ligand_tla = ligand_code_tla(ligand_code)
    _, pdb_ligand_chain = _pdb_chain_ids(protein_chain, ligand_chain)
    paths = expected_output_paths(pdb_id, protein_chain, ligand_chain, ligand_code, outdir)
    paths["ligand"] = f"{ligand_tla}_{pdb_ligand_chain}"
    pdb_path = Path(paths["pdb_path"])
    ligand_path = Path(paths["ligand_path"])

    if not overwrite and pdb_path.is_file() and ligand_path.is_file():
        return paths

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        structure_path = download_structure(pdb_id, tmp)
        if structure_path is None:
            return None

        download_ligand_sdf(pdb_id, ligand_chain, ligand_code, ligand_path)
        structure = _load_structure(pdb_id, structure_path)

        pdb_protein_chain, pdb_ligand_chain = _pdb_chain_ids(protein_chain, ligand_chain)
        chain_map = {}
        if protein_chain != pdb_protein_chain:
            chain_map[protein_chain] = pdb_protein_chain
        if ligand_chain != pdb_ligand_chain:
            chain_map[ligand_chain] = pdb_ligand_chain
        if chain_map:
            _remap_structure_chains(structure, chain_map)

        protein_only_pdb = tmp / "protein_only.pdb"
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(str(protein_only_pdb), _ProteinOnlySelect(pdb_protein_chain))

        repaired_protein_pdb = tmp / "repaired_protein.pdb"
        evoef2_repair_structure(protein_only_pdb, repaired_protein_pdb, evoef2_bin)
        _merge_repaired_protein_with_ligand_sdf(
            repaired_protein_pdb, ligand_path, pdb_ligand_chain, ligand_tla, pdb_path
        )

    return paths


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Prepare MaSIF input PDB/SDF files from a seed CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_csv", required=True, help="Seed CSV with columns: pdb_id, protein_chain, ligand_chain, ligand_code.")
    parser.add_argument("--outdir", required=True, help="Directory to write .pdb and .sdf files.")
    parser.add_argument("--out_csv", required=True, help="Output CSV path (input rows + manifest columns).")
    parser.add_argument("--evoef2_bin", required=True, help="Path to EvoEF2 binary.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Rebuild even when output files already exist.",
    )
    args = parser.parse_args()

    evoef2_bin = Path(args.evoef2_bin)
    if not evoef2_bin.is_file():
        print(
            f"Error: EvoEF2 binary not found at {evoef2_bin}. "
            "Build it first: cd EvoEF2 && ./build.sh",
            file=sys.stderr,
        )
        sys.exit(1)

    df = pd.read_csv(args.input_csv)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        print(
            f"Error: --input_csv is missing required columns: {sorted(missing)}",
            file=sys.stderr,
        )
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_cols = ["pdb_path", "target", "ligand", "ligand_path"]
    for col in manifest_cols:
        df[col] = None

    n_failed = 0
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Preparing structures"):
        pdb_id = str(row.pdb_id)
        try:
            result = prepare_input_structures(
                pdb_id=pdb_id,
                protein_chain=str(row.protein_chain),
                ligand_chain=str(row.ligand_chain),
                ligand_code=str(row.ligand_code),
                outdir=outdir,
                evoef2_bin=evoef2_bin,
                overwrite=args.overwrite,
            )
        except Exception as exc:
            tqdm.write(f"Error processing {pdb_id}: {exc}", file=sys.stderr)
            result = None

        if result is None:
            n_failed += 1
        else:
            for col in manifest_cols:
                df.at[idx, col] = result[col]

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}  ({len(df) - n_failed}/{len(df)} rows succeeded)")

    if n_failed:
        print(f"{n_failed} row(s) failed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
