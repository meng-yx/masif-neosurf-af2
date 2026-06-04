"""
Resolve neosurf ligand anchor residues from preprocessed PDB structures.
"""
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBParser

PDB_PARSER = PDBParser(QUIET=True)

SEED_SKIP_WARNING = (
    "WARNING: skipped {protein_id} with no HETATOM record. "
    "If no ligand is expected, remove --seed-auto-neosurf flag."
)


@dataclass(frozen=True)
class LigandAnchor:
    chain: str
    resid: int
    resname: str
    n_heavy_atoms: int


def _heavy_atom_count(residue):
    return sum(1 for atom in residue.get_atoms() if atom.element != "H")


def find_ligand_anchor(pdb_path):
    """
    Return the hetero residue with the most heavy atoms (BioPython hetero flag).

    Ties are broken by first occurrence in the structure.
    Returns None if no hetero residues are found.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    structure = PDB_PARSER.get_structure(pdb_path.stem, str(pdb_path))
    groups = OrderedDict()
    for residue in structure.get_residues():
        if not residue.id[0].strip():
            continue
        key = (residue.parent.id, residue.id[1], residue.get_resname().strip())
        groups.setdefault(key, []).append(residue)

    if not groups:
        return None

    best_key = None
    best_count = -1
    for key, residues in groups.items():
        count = sum(_heavy_atom_count(r) for r in residues)
        if count > best_count:
            best_count = count
            best_key = key

    chain, resid, resname = best_key
    return LigandAnchor(chain=chain, resid=int(resid), resname=resname, n_heavy_atoms=best_count)


def log_ligand_anchor(protein_id, anchor, pdb_path=None):
    suffix = f" ({pdb_path})" if pdb_path else ""
    print(
        f"INFO: {protein_id} HET anchor chain={anchor.chain} resid={anchor.resid} "
        f"resname={anchor.resname} ({anchor.n_heavy_atoms} heavy atoms){suffix}"
    )
