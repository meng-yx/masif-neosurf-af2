#!/usr/bin/env python
import argparse
import os
import tempfile

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from default_config.masif_opts import masif_opts


class RejectResiduesSelect(Select):
    def __init__(self, reject_keys):
        self.reject_keys = reject_keys

    def accept_residue(self, residue):
        key = (residue.get_parent().get_id(), residue.get_id())
        return 0 if key in self.reject_keys else 1


def remove_backbone_incomplete_residues(pdb_path, required_atoms=("N", "CA", "C", "O")):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("repair_backbone_filter", pdb_path)
    model = next(struct.get_models())
    reject = set()
    removed = []
    total = 0
    for chain in model:
        for residue in chain:
            if residue.get_id()[0] != " " or not is_aa(residue, standard=True):
                continue
            total += 1
            atom_names = {atom.get_name() for atom in residue.get_atoms()}
            missing = [atm for atm in required_atoms if atm not in atom_names]
            if missing:
                reject.add((chain.get_id(), residue.get_id()))
                removed.append(
                    {
                        "chain": chain.get_id(),
                        "resname": residue.get_resname(),
                        "resseq": residue.get_id()[1],
                        "icode": residue.get_id()[2].strip(),
                        "missing_backbone_atoms": missing,
                    }
                )
    if reject:
        io = PDBIO()
        io.set_structure(struct)
        io.save(pdb_path, select=RejectResiduesSelect(reject))
    return {"total_standard_residues": total, "removed_count": len(removed), "removed": removed}


def run_repair(pdb_id):
    raw_path = os.path.join(masif_opts["raw_pdb_dir"], "{}.pdb".format(pdb_id))
    if not os.path.isfile(raw_path):
        raise FileNotFoundError("Raw PDB not found for repair: {}".format(raw_path))

    bb = remove_backbone_incomplete_residues(raw_path)

    fixer = PDBFixer(filename=raw_path)
    fixer.findMissingResidues()
    # Do not add missing sequence residues; only fill atoms in existing residues.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    missing_atom_targets = len(fixer.missingAtoms)
    fixer.addMissingAtoms()

    raw_dir = os.path.dirname(os.path.abspath(raw_path))
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False, dir=raw_dir) as tmp_handle:
        tmp_out = tmp_handle.name
        PDBFile.writeFile(fixer.topology, fixer.positions, tmp_handle, keepIds=True)
    try:
        os.replace(tmp_out, raw_path)
    except Exception:
        raise
    print(
        "Repair summary for {}: removed_backbone_residues={}, sidechain_targets={}".format(
            pdb_id, bb["removed_count"], int(missing_atom_targets)
        )
    )


def main():
    parser = argparse.ArgumentParser(description="Backbone QC + sidechain repair using PDBFixer.")
    parser.add_argument("pdb_id", help="PDB id, e.g. 5JXC-C1-B1-PP")
    args = parser.parse_args()
    run_repair(args.pdb_id)


if __name__ == "__main__":
    main()
