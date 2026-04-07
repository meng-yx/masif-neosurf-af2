#!/usr/bin/python
import argparse
import os
import shutil
import sys

from default_config.masif_opts import masif_opts
# Local includes
from input_output.protonate import protonate

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare raw PDB file from Pinder DB."
    )
    parser.add_argument(
        "complex_id",
        help="Complex identifier in the form PDBID_A_B.",
    )
    parser.add_argument(
        "--pinder-pdb-dir",
        required=True,
        help="Directory containing pre-existing .pdb files, using <complex_id>.pdb naming.",
    )
    return parser.parse_args()


def ensure_dirs():
    if not os.path.exists(masif_opts["raw_pdb_dir"]):
        os.makedirs(masif_opts["raw_pdb_dir"])
    if not os.path.exists(masif_opts["tmp_dir"]):
        os.mkdir(masif_opts["tmp_dir"])


def get_input_pdb(complex_id, pdb_id, pinder_pdb_dir):
    pdb_filename = os.path.join(masif_opts["tmp_dir"], "pdb{}.ent".format(pdb_id.lower()))
    src_pdb = os.path.join(pinder_pdb_dir, "{}.pdb".format(complex_id))
    if not os.path.isfile(src_pdb):
        print("Error: Pinder PDB file not found: {}".format(src_pdb))
        sys.exit(1)
    print("Copying PDB structure from Pinder DB: {}".format(src_pdb))
    shutil.copyfile(src_pdb, pdb_filename)
    return pdb_filename


def main():
    args = parse_args()
    ensure_dirs()

    in_fields = args.complex_id.split("_")
    pdb_id = in_fields[0]
    pdb_filename = get_input_pdb(args.complex_id, pdb_id, args.pinder_pdb_dir)

    # Always protonate; downstream steps can ignore hydrogens if needed.
    protonated_file = masif_opts["raw_pdb_dir"] + "/" + pdb_id + ".pdb"
    protonate(pdb_filename, protonated_file)


if __name__ == "__main__":
    main()

