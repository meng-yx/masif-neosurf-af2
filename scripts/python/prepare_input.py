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
  4. Merges the ligand onto the EvoEF2-repaired protein: atom coordinates come
     from the downloaded SDF; chain ID, residue number, and residue name come
     from the matching ligand in the downloaded structure. Atom names are
     generated from element symbols (C01, N01, ...) to avoid clashes with
     standard amino-acid nomenclature in downstream preprocessing.
  5. Appends seven manifest columns to every input row:
       pdb_protein_chain, pdb_ligand_chain, ligand_resname,
       pdb_path, target, ligand, ligand_path
     (empty when a row fails) and writes the result to --out_csv.

Identifier scheme
-----------------
Author-assigned identifiers from the seed CSV (protein_chain, ligand_chain,
ligand_code) are preserved unchanged in the output CSV as an audit trail. All
downstream pipeline fields use PDB-normalized identifiers:

  pdb_protein_chain / pdb_ligand_chain
      Single-character chain IDs as written in the output PDB. Multi-letter
      author chains (e.g. "AAA") are remapped to the shortest available letter.

  ligand_resname
      The ≤3-character residue name written in HETATM records, taken from the
      deposited structure. Long ligand_codes (e.g. "A1H3Q") are stored here as
      their deposited abbreviation (e.g. "A1H").

  target = {pdb_id}-{ligand_code}_{pdb_chain_suffix}
      The MaSIF ppi_pair_id / seed_id. Uses the full author ligand_code for
      unambiguous identification, and PDB-normalized chain IDs in the suffix.
      Example: author chains AAA/AAA, ligand M0B → target "6SYP-M0B_A".

  ligand = {ligand_resname}_{pdb_ligand_chain}
      Passed to preprocess as the -l argument; matches PDB HETATM content.

  ligand_path = {outdir}/{pdb_id}_{pdb_protein_chain}_{ligand_resname}.sdf

Usage:
    python scripts/python/prepare_input.py \\
        --input_csv data/processing/1_prep_input/input_subsets/input_1.csv \\
        --outdir    data/input \\
        --out_csv   data/processing/1_prep_input/input_1_manifest.csv \\
        --evoef2_bin EvoEF2/EvoEF2 \\
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
# Path / identifier helpers
# ---------------------------------------------------------------------------


def ligand_code_tla(ligand_code):
    """Return the first 3 characters of a ligand code (for RCSB API and resname matching)."""
    return ligand_code.strip()[:3]


def pdb_resname(resname):
    """Normalize a residue name to the 3-character PDB HETATM field (columns 18-20)."""
    return str(resname).strip().upper()[:3]


def chain_suffix(protein_chain, ligand_chain):
    """Return deduplicated chain suffix from author chain IDs: C+A -> CA, A+A -> A.

    Kept for backward compatibility. Manifest construction uses pdb_chain_suffix()
    with PDB-normalized single-character chain IDs instead.
    """
    if protein_chain == ligand_chain:
        return protein_chain
    return protein_chain + ligand_chain


def pdb_chain_suffix(pdb_protein_chain, pdb_ligand_chain):
    """Return deduplicated chain suffix from PDB-normalized single-char chain IDs."""
    if pdb_protein_chain == pdb_ligand_chain:
        return pdb_protein_chain
    return pdb_protein_chain + pdb_ligand_chain


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


def build_manifest_paths(
    pdb_id, ligand_code, ligand_resname, pdb_protein_chain, pdb_ligand_chain, outdir
):
    """Return the canonical manifest dict using PDB-normalized identifiers.

    All downstream pipeline steps should use these fields rather than re-deriving
    identifiers from the author seed columns.
    """
    outdir = Path(outdir)
    chains = pdb_chain_suffix(pdb_protein_chain, pdb_ligand_chain)
    target = f"{pdb_id}-{ligand_code}_{chains}"
    ligand = f"{ligand_resname}_{pdb_ligand_chain}"
    pdb_path = outdir / f"{target}.pdb"
    ligand_path = outdir / f"{pdb_id}_{pdb_protein_chain}_{ligand_resname}.sdf"
    return {
        "pdb_protein_chain": pdb_protein_chain,
        "pdb_ligand_chain": pdb_ligand_chain,
        "ligand_resname": ligand_resname,
        "pdb_path": str(pdb_path.resolve()),
        "target": target,
        "ligand": ligand,
        "ligand_path": str(ligand_path.resolve()),
    }


def _read_ligand_resname_from_pdb(pdb_path, pdb_ligand_chain):
    """Read the residue name of the first HETATM on pdb_ligand_chain in an existing PDB.

    Used to reconstruct ligand_resname for the skip-if-exists check without
    downloading the structure again.  Returns None if the file is absent or
    contains no matching HETATM record.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        return None
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith("HETATM") and len(line) >= 26:
                if line[21] == pdb_ligand_chain:
                    return line[17:20].strip() or None
    return None


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
    auth_seq_id = None
    for entity_id in entity_ids:
        entity = _fetch_json(f"{RCSB_DATA_API}/nonpolymer_entity/{pdb_id}/{entity_id}")
        identifiers = entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
        comp_id = identifiers.get("nonpolymer_comp_id", "").strip().upper()
        if comp_id not in {ligand_code, ligand_tla}:
            continue

        # asym_ids and auth_asym_ids are independently sorted on the entity
        # endpoint, so zip() does not pair label/instance chains correctly.
        for label_id in identifiers.get("asym_ids", []):
            instance = _fetch_json(
                f"{RCSB_DATA_API}/nonpolymer_entity_instance/{pdb_id}/{label_id}"
            )
            instance_ids = instance.get(
                "rcsb_nonpolymer_entity_instance_container_identifiers", {}
            )
            inst_auth_asym_id = str(instance_ids.get("auth_asym_id", "")).strip()
            if inst_auth_asym_id != auth_asym_id:
                continue
            label_asym_id = str(label_id).strip()
            auth_seq_id = str(instance_ids.get("auth_seq_id", "")).strip()
            break
        if label_asym_id is not None:
            break

    if label_asym_id is None:
        raise ValueError(
            f"No non-polymer entity for {pdb_id} with comp_id={ligand_code} "
            f"on auth chain {auth_asym_id}"
        )

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


def _resname_matches(resname, ligand_code):
    """Return True when a PDB residue name matches the requested ligand code."""
    resname = resname.strip().upper()
    ligand_code = ligand_code.strip().upper()
    return resname == ligand_code or resname == ligand_code_tla(ligand_code)


def _format_pdb_atom_name(element, index_within_element):
    """Return a 4-character PDB atom name with correct element alignment."""
    element = element.strip().upper()
    if len(element) == 2 and element.isalpha():
        stem = f"{element}{index_within_element}"
        return f"{stem:<4}"[:4]
    stem = f"{element}{index_within_element}"
    return f" {stem:<3}"[:4]


def _normalize_pdb_atom_name(atom_name, element):
    """Normalize an atom name to the 4-character PDB format."""
    name = atom_name.strip()
    if len(name) == 4:
        return name
    if len(name) == 3:
        return f" {name}"
    if len(name) == 2:
        return f" {name} "
    if len(name) == 1:
        return f" {name}  "
    return _format_pdb_atom_name(element, 1)


def _get_ligand_placement(structure, pdb_ligand_chain, ligand_code):
    """
    Read chain ID, residue number, and atom metadata for the target ligand.

    Returns a placement dict with chain_id, residue_id, resname, and atoms
    (name/element/occupancy/bfactor) from the downloaded structure.
    """
    ligand_code = ligand_code.strip().upper()
    matching_residues = []
    for model in structure:
        for chain in model:
            if chain.id != pdb_ligand_chain:
                continue
            for residue in chain:
                hetflag, resseq, icode = residue.id
                if hetflag.strip() and _resname_matches(residue.resname, ligand_code):
                    matching_residues.append(residue)

    if not matching_residues:
        raise ValueError(
            f"No HETATM residue matching ligand_code={ligand_code!r} found on chain "
            f"{pdb_ligand_chain!r} in the downloaded structure."
        )

    unique_resnums = {r.id[1] for r in matching_residues}
    if len(unique_resnums) > 1:
        warnings.warn(
            f"Multiple instances of ligand {ligand_code} on chain {pdb_ligand_chain} "
            f"(resseq: {sorted(unique_resnums)}); using the first instance.",
            UserWarning,
            stacklevel=2,
        )

    first_resseq = matching_residues[0].id[1]
    selected = [r for r in matching_residues if r.id[1] == first_resseq]
    template_residue = selected[0]
    atoms = []
    for atom in template_residue.get_atoms():
        element = atom.element.strip() if atom.element else atom.get_id().strip()[:1]
        atoms.append(
            {
                "name": _normalize_pdb_atom_name(atom.get_name(), element),
                "element": element,
                "occupancy": atom.get_occupancy(),
                "bfactor": atom.get_bfactor(),
            }
        )

    if not atoms:
        raise ValueError(
            f"Ligand residue {ligand_code} on chain {pdb_ligand_chain} has no atoms."
        )

    raw_resname = template_residue.resname.strip()
    resname = pdb_resname(raw_resname)
    if len(raw_resname) > 3:
        warnings.warn(
            f"Deposited ligand resname {raw_resname!r} exceeds PDB 3-character limit; "
            f"using {resname!r} in HETATM records.",
            UserWarning,
            stacklevel=2,
        )

    return {
        "chain_id": pdb_ligand_chain,
        "residue_id": template_residue.id,
        "resname": resname,
        "atoms": atoms,
    }


def _format_hetatm_line(
    serial, atom_name, resname, chain_id, resseq, x, y, z, occupancy, bfactor, element
):
    """Write an 80-column HETATM record with the element symbol in columns 77-78."""
    resname = pdb_resname(resname)
    return (
        f"HETATM{serial:5d} {atom_name:>4s} {resname:>3s} {chain_id:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}  \n"
    )


def _last_pdb_serial(pdb_path):
    last_serial = 0
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                last_serial = int(line[6:11])
    return last_serial


def _ligand_hetatm_lines_from_sdf(sdf_path, placement, start_serial):
    """Return formatted HETATM lines for a ligand using SDF coordinates.

    Atom names are generated from element symbols (C01, N01, ...) rather than
    copied from the deposited structure, because ligands such as 3JF use atom
    names (CA, CB, CG, ...) that overlap with standard amino-acid nomenclature
    and break OpenBabel/PDB2PQR during MaSIF preprocessing.
    """
    mol = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not read ligand SDF: {sdf_path}")

    conf = mol.GetConformer()
    sdf_atoms = list(mol.GetAtoms())
    pdb_atoms = placement["atoms"]
    if len(sdf_atoms) != len(pdb_atoms):
        warnings.warn(
            f"SDF atom count ({len(sdf_atoms)}) differs from downloaded PDB ligand "
            f"({len(pdb_atoms)}); using SDF atom order and generated atom names.",
            UserWarning,
            stacklevel=2,
        )

    hetflag, resseq, icode = placement["residue_id"]
    if icode and icode != " ":
        raise ValueError(
            f"Ligand insertion codes are not supported (got {icode!r} for resseq {resseq})."
        )

    chain_id = placement["chain_id"]
    resname = placement["resname"]
    lines = []
    element_counts = {}
    serial = start_serial

    for index, sdf_atom in enumerate(sdf_atoms):
        element = sdf_atom.GetSymbol()
        pos = conf.GetAtomPosition(index)
        element_counts[element] = element_counts.get(element, 0) + 1
        atom_name = f"{element}{element_counts[element]:02d}"[:4]
        if index < len(pdb_atoms):
            occupancy = pdb_atoms[index]["occupancy"]
            bfactor = pdb_atoms[index]["bfactor"]
        else:
            occupancy = 1.0
            bfactor = 0.0

        serial += 1
        lines.append(
            _format_hetatm_line(
                serial,
                atom_name,
                resname,
                chain_id,
                resseq,
                pos.x,
                pos.y,
                pos.z,
                occupancy,
                bfactor,
                element,
            )
        )

    return lines


def _merge_repaired_protein_with_ligand(
    repaired_pdb,
    placement,
    ligand_sdf,
    output_pdb,
):
    """
    Merge EvoEF2-repaired protein with SDF-based ligand HETATM records.

    Ligand coordinates come from the SDF; chain ID, residue number, and residue
    name are inherited from the pre-computed placement dict. Atom names are
    generated from element symbols. Ligand records are written with explicit
    80-column formatting so downstream prody/RDKit parsing sees the element
    symbol in columns 77-78.
    """
    repaired_pdb = Path(repaired_pdb)
    start_serial = _last_pdb_serial(repaired_pdb)
    hetatm_lines = _ligand_hetatm_lines_from_sdf(ligand_sdf, placement, start_serial)

    out_lines = []
    with open(repaired_pdb) as handle:
        for line in handle:
            if line.startswith("END"):
                break
            if line.strip():
                out_lines.append(line if line.endswith("\n") else f"{line}\n")
    out_lines.extend(hetatm_lines)
    out_lines.append("END\n")

    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("".join(out_lines))


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

    Returns a manifest dict on success, or None on failure.  The dict includes
    the seven manifest columns (pdb_protein_chain, pdb_ligand_chain,
    ligand_resname, pdb_path, target, ligand, ligand_path).

    When overwrite=False and both output files already exist, ligand_resname is
    read from the first HETATM in the existing PDB and the manifest is returned
    immediately without re-downloading or reprocessing.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdb_protein_chain, pdb_ligand_chain = _pdb_chain_ids(protein_chain, ligand_chain)

    # The PDB output path is deterministic (no ligand_resname needed).
    chains = pdb_chain_suffix(pdb_protein_chain, pdb_ligand_chain)
    target_base = f"{pdb_id}-{ligand_code}_{chains}"
    candidate_pdb_path = outdir / f"{target_base}.pdb"

    if not overwrite and candidate_pdb_path.is_file():
        ligand_resname = _read_ligand_resname_from_pdb(candidate_pdb_path, pdb_ligand_chain)
        if ligand_resname is not None:
            paths = build_manifest_paths(
                pdb_id, ligand_code, ligand_resname,
                pdb_protein_chain, pdb_ligand_chain, outdir,
            )
            if Path(paths["ligand_path"]).is_file():
                return paths

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        structure_path = download_structure(pdb_id, tmp)
        if structure_path is None:
            return None

        structure = _load_structure(pdb_id, structure_path)

        chain_map = {}
        if protein_chain != pdb_protein_chain:
            chain_map[protein_chain] = pdb_protein_chain
        if ligand_chain != pdb_ligand_chain:
            chain_map[ligand_chain] = pdb_ligand_chain
        if chain_map:
            _remap_structure_chains(structure, chain_map)

        # Resolve ligand_resname from the deposited structure (≤3-char HETATM resname).
        placement = _get_ligand_placement(structure, pdb_ligand_chain, ligand_code)
        ligand_resname = placement["resname"]

        paths = build_manifest_paths(
            pdb_id, ligand_code, ligand_resname,
            pdb_protein_chain, pdb_ligand_chain, outdir,
        )
        pdb_path = Path(paths["pdb_path"])
        ligand_path = Path(paths["ligand_path"])

        download_ligand_sdf(pdb_id, ligand_chain, ligand_code, ligand_path)

        protein_only_pdb = tmp / "protein_only.pdb"
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(str(protein_only_pdb), _ProteinOnlySelect(pdb_protein_chain))

        repaired_protein_pdb = tmp / "repaired_protein.pdb"
        evoef2_repair_structure(protein_only_pdb, repaired_protein_pdb, evoef2_bin)
        _merge_repaired_protein_with_ligand(
            repaired_protein_pdb,
            placement,
            ligand_path,
            pdb_path,
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

    manifest_cols = [
        "pdb_protein_chain", "pdb_ligand_chain", "ligand_resname",
        "pdb_path", "target", "ligand", "ligand_path",
    ]
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
