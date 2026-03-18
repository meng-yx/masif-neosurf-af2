#!/usr/bin/python
import numpy as np
import os
import Bio
import shutil
from Bio.PDB import *
from Bio.PDB.Polypeptide import is_aa
import sys
import importlib
from IPython.core.debugger import set_trace

# Local includes
from default_config.masif_opts import masif_opts
from triangulation.computeMSMS import computeMSMS
from triangulation.fixmesh import fix_mesh
from triangulation.ligand_utils import extract_ligand
import pymesh
from input_output.extractPDB import extractPDB
from input_output.save_ply import save_ply
from input_output.read_ply import read_ply
from input_output.protonate import protonate, get_pdb_conect
from triangulation.computeHydrophobicity import computeHydrophobicity
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from triangulation.computeAPBS import computeAPBS
from triangulation.compute_normal import compute_normal
from sklearn.neighbors import KDTree


def _residue_key(chain_id, residue):
    return (chain_id, residue.get_id())


class _ResidueRejectSelect(Select):
    def __init__(self, reject_keys):
        self.reject_keys = reject_keys

    def accept_residue(self, residue):
        parent_chain = residue.get_parent()
        key = _residue_key(parent_chain.get_id(), residue)
        return 0 if key in self.reject_keys else 1


def filter_incomplete_backbone_residues(
    pdb_path,
    required_backbone_atoms=("N", "CA", "C", "O"),
    fail_if_filtered_fraction_exceeds=None,
):
    """
    Remove standard amino-acid residues missing required backbone heavy atoms.
    Keeps non-protein residues unchanged.
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("filter_backbone", pdb_path)
    model = Selection.unfold_entities(struct, "M")[0]

    reject_keys = set()
    removed_residues = []
    total_protein_residues = 0

    for chain in model:
        for residue in chain:
            hetflag = residue.get_id()[0]
            # Only filter standard polymer amino acids.
            if hetflag != " " or not is_aa(residue, standard=True):
                continue
            total_protein_residues += 1
            atom_names = set(atom.get_name() for atom in residue.get_atoms())
            missing = [atm for atm in required_backbone_atoms if atm not in atom_names]
            if len(missing) > 0:
                key = _residue_key(chain.get_id(), residue)
                reject_keys.add(key)
                removed_residues.append(
                    {
                        "resname": residue.get_resname(),
                        "chain": chain.get_id(),
                        "resseq": residue.get_id()[1],
                        "icode": residue.get_id()[2].strip(),
                        "missing": missing,
                    }
                )

    removed_count = len(removed_residues)
    kept_count = total_protein_residues - removed_count
    filtered_fraction = (
        float(removed_count) / float(total_protein_residues)
        if total_protein_residues > 0
        else 0.0
    )

    # Rewrite the pdb in-place only if anything is removed.
    if removed_count > 0:
        io = PDBIO()
        io.set_structure(struct)
        io.save(pdb_path, select=_ResidueRejectSelect(reject_keys))

    return {
        "total_protein_residues": total_protein_residues,
        "removed_count": removed_count,
        "kept_count": kept_count,
        "filtered_fraction": filtered_fraction,
        "removed_residues": removed_residues,
        "required_backbone_atoms": list(required_backbone_atoms),
    }


if len(sys.argv) <= 1:
    print("Usage: {config} "+sys.argv[0]+" PDBID_A")
    print("A or AB are the chains to include in this surface.")
    sys.exit(1)


# Save the chains as separate files.
in_fields = sys.argv[1].split("_")
pdb_id = in_fields[0]
chain_ids1 = in_fields[1]

if (len(sys.argv)>2) and (sys.argv[2]=='masif_ligand'):
    pdb_filename = os.path.join(masif_opts["ligand"]["assembly_dir"],pdb_id+".pdb")
else:
    pdb_filename = masif_opts['raw_pdb_dir']+pdb_id+".pdb"

# With ligand
if len(sys.argv) >= 3:
    ligand_code, ligand_chain = sys.argv[2].split('_')

    if len(sys.argv) >= 4:
        sdf_file = sys.argv[3]
    else:
        print("No SDF file provided. Connectivity will be inferred from the PDB.")
        sdf_file = None

# Without ligand
else:
    ligand_code, ligand_chain = None, None
    sdf_file = None


mol2_patch = None
if len(sys.argv) >= 5:
    mol2_patch = sys.argv[4]

if ligand_code is not None:
    assert len(ligand_code) <= 3, "Only 3-letter ligand codes are supported."
    print("Including ligand {} in the surface.".format(ligand_code))


tmp_dir = os.environ.get('TMPDIR')
if tmp_dir is None:
    tmp_dir= masif_opts['tmp_dir']

if ligand_code is not None and ligand_code == 'UNL' and sdf_file is not None:
    het_dict = os.path.join(tmp_dir, "{}_{}_conect.txt".format(pdb_id, ligand_code))
    get_pdb_conect(pdb_filename, ligand_code, ligand_chain, sdf_file, save_txt=het_dict)
else:
    het_dict = os.environ.get('REDUCE_HET_DICT')
protonated_file = tmp_dir+"/"+pdb_id+"_protonated.pdb"
protonate(pdb_filename, protonated_file, het_dict=het_dict)
pdb_filename = protonated_file

# Extract chains of interest.
out_filename1 = tmp_dir+"/"+pdb_id+"_"+chain_ids1
extractPDB(pdb_filename, out_filename1+".pdb", chain_ids1, ligand_code, ligand_chain)

if masif_opts.get("filter_incomplete_backbone_residues", True):
    required_backbone_atoms = masif_opts.get(
        "min_backbone_atoms_required", ["N", "CA", "C", "O"]
    )
    fail_if_filtered_fraction_exceeds = masif_opts.get(
        "fail_if_filtered_fraction_exceeds", None
    )
    filter_info = filter_incomplete_backbone_residues(
        out_filename1 + ".pdb",
        required_backbone_atoms=tuple(required_backbone_atoms),
        fail_if_filtered_fraction_exceeds=fail_if_filtered_fraction_exceeds,
    )
    if filter_info["removed_count"] > 0:
        print(
            "Removed {} residues with incomplete backbone (required atoms: {}).".format(
                filter_info["removed_count"], ",".join(filter_info["required_backbone_atoms"])
            )
        )
        for rr in filter_info["removed_residues"]:
            icode = rr["icode"] if rr["icode"] != "" else "-"
            print(
                "  - {} {} {} {} missing [{}]".format(
                    rr["resname"],
                    rr["chain"],
                    rr["resseq"],
                    icode,
                    ",".join(rr["missing"]),
                )
            )
    else:
        print("Backbone filter removed 0 residues.")
    if filter_info["kept_count"] <= 0:
        raise RuntimeError(
            "Backbone filter removed all protein residues for {}. Aborting chain processing.".format(
                out_filename1 + ".pdb"
            )
        )
    if (
        fail_if_filtered_fraction_exceeds is not None
        and filter_info["filtered_fraction"] > float(fail_if_filtered_fraction_exceeds)
    ):
        raise RuntimeError(
            "Backbone filter removed {:.1f}% residues in {} (threshold {:.1f}%).".format(
                100.0 * filter_info["filtered_fraction"],
                out_filename1 + ".pdb",
                100.0 * float(fail_if_filtered_fraction_exceeds),
            )
        )

# Compute MSMS of surface w/hydrogens,
vertices1, faces1, normals1, names1, areas1 = computeMSMS(out_filename1+".pdb", protonate=True, keep_hetatms=[ligand_code])

# Get and RDKit molecule object
if ligand_code is not None and ligand_chain is not None:
    mol2_file = os.path.join(tmp_dir, "{}_{}.mol2".format(ligand_code, ligand_chain))
    rdmol = extract_ligand(pdb_filename, ligand_code, ligand_chain, mol2_file, sdf_template=sdf_file, patched_mol2_file=mol2_patch)
else:
    mol2_file = None
    rdmol = None

# Compute "charged" vertices
if masif_opts['use_hbond']:
    vertex_hbond = computeCharges(out_filename1, vertices1, names1, ligand_code, rdmol)

# For each surface residue, assign the hydrophobicity of its amino acid.
if masif_opts['use_hphob']:
    vertex_hphobicity = computeHydrophobicity(names1, ligand_code, rdmol)

# If protonate = false, recompute MSMS of surface, but without hydrogens (set radius of hydrogens to 0).
vertices2 = vertices1
faces2 = faces1

# Fix the mesh.
mesh = pymesh.form_mesh(vertices2, faces2)
print(f"Fixing mesh...")
regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])
print(f"Fixmesh done!")

# Compute the normals
vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)
# Assign charges on new vertices based on charges of old vertices (nearest
# neighbor)

if masif_opts['use_hbond']:
    vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices1,\
        vertex_hbond, masif_opts)

if masif_opts['use_hphob']:
    vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices1,\
        vertex_hphobicity, masif_opts)

if masif_opts['use_apbs']:
    print(f"Computing APBS...")
    vertex_charges = computeAPBS(regular_mesh.vertices, out_filename1+".pdb", out_filename1, mol2_file)
    print(f"APBS done!")

iface = np.zeros(len(regular_mesh.vertices))
if 'compute_iface' in masif_opts and masif_opts['compute_iface']:
    # Compute the surface of the entire complex and from that compute the interface.
    v3, f3, _, _, _ = computeMSMS(pdb_filename, protonate=True, keep_hetatms=[ligand_code])
    # Regularize the mesh
    mesh = pymesh.form_mesh(v3, f3)
    # I believe It is not necessary to regularize the full mesh. This can speed up things by a lot.
    full_regular_mesh = mesh
    # Find the vertices that are in the iface.
    v3 = full_regular_mesh.vertices
    # Find the distance between every vertex in regular_mesh.vertices and those in the full complex.
    kdt = KDTree(v3)
    d, r = kdt.query(regular_mesh.vertices)
    d = np.square(d) # Square d, because this is how it was in the pyflann version.
    assert(len(d) == len(regular_mesh.vertices))
    iface_v = np.where(d >= 2.0)[0]
    iface[iface_v] = 1.0
    # Convert to ply and save.
    save_ply(out_filename1+".ply", regular_mesh.vertices,\
                        regular_mesh.faces, normals=vertex_normal, charges=vertex_charges,\
                        normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity,\
                        iface=iface)

else:
    # Convert to ply and save.
    save_ply(out_filename1+".ply", regular_mesh.vertices,\
                        regular_mesh.faces, normals=vertex_normal, charges=vertex_charges,\
                        normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity)
if not os.path.exists(masif_opts['ply_chain_dir']):
    os.makedirs(masif_opts['ply_chain_dir'])
if not os.path.exists(masif_opts['pdb_chain_dir']):
    os.makedirs(masif_opts['pdb_chain_dir'])
shutil.copy(out_filename1+'.ply', masif_opts['ply_chain_dir'])
shutil.copy(out_filename1+'.pdb', masif_opts['pdb_chain_dir'])
