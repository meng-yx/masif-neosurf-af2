import os
import sys
import shutil
from pathlib import Path
from argparse import ArgumentParser
from packaging import version

import pymesh
from scipy.spatial import cKDTree
import numpy as np
from Bio.PDB import PDBParser
import open3d as o3d

# import MaSIF modules
masif_neosurf_dir = Path(__file__).resolve().parent
sys.path.append(str(Path(masif_neosurf_dir, 'masif', 'source').resolve()))
sys.path.append(str(Path(masif_neosurf_dir, 'masif_seed_search', 'source').resolve()))
from masif.source.default_config.masif_opts import masif_opts
from alignment_evaluation_nn import AlignmentEvaluationNN
from alignment_utils import get_patch_coords, load_protein_pcd, get_patch_geo, match_descriptors, align_protein, compute_nn_score
from search_output import filter_seeds_for_resume, write_site_seed_hits
from structural_cluster import cluster_search_hits_for_subset
from neosurf_anchor import find_ligand_anchor, log_ligand_anchor


def _grid_subsample_vertices(o3d_mesh, candidates, n_keep, dists_to_anchor):
    """Poisson-disk subsample on full mesh; return n_keep centers closest to anchor (search_grid.py)."""
    mesh_coords = np.asarray(o3d_mesh.vertices)
    n_mesh = len(mesh_coords)
    candidates = np.asarray(candidates)
    if n_keep >= len(candidates):
        return candidates
    n_poisson = int(n_keep * n_mesh / len(candidates))
    sampled_points = o3d_mesh.sample_points_poisson_disk(n_poisson)
    squared_dists = np.sum(
        np.square(
            np.asarray(sampled_points.points).reshape(-1, 1, 3) - mesh_coords.reshape(1, -1, 3)
        ),
        axis=-1,
    )
    mesh_inds = np.argmin(squared_dists, axis=-1)
    order = np.argsort(dists_to_anchor[mesh_inds])
    return mesh_inds[order[:n_keep]]


def _target_dists_to_anchor(mymesh, target_struct, params):
    """Per-vertex Euclidean distance from mesh vertices to target anchor."""
    vertices = np.asarray(mymesh.vertices)
    if 'target_point' in params:
        anchor_coord = np.array(params['target_point']['coord'])
        return np.sqrt(np.sum(np.square(vertices - anchor_coord), axis=1))

    target_chain = params['target_residue']['chain']
    target_resid = [
        x.id for x in target_struct[0][target_chain].get_residues()
        if x.id[1] == params['target_residue']['resid']
    ]
    if len(target_resid) != 1:
        raise RuntimeError(f"Target residue ID not unique: {target_resid}")
    target_resid = target_resid[0]
    print(f"Using residue: {target_resid}")

    if 'target_atom' in params:
        coord = target_struct[0][target_chain][target_resid][params['target_atom']['atom_id']].get_coord()
        return np.sqrt(np.sum(np.square(vertices - coord), axis=1))

    residue_coords = np.stack([
        a.get_coord() for a in target_struct[0][target_chain][target_resid].get_atoms()
        if a.element != 'H'
    ])
    return np.sqrt(np.min(
        np.sum(np.square(vertices.reshape(-1, 1, 3) - residue_coords.reshape(1, -1, 3)), axis=-1),
        axis=1,
    ))


def _select_target_vertices(mymesh, target_struct, target_ply_fn, params):
    """Select target patch center vertex indices."""
    if params.get('target_site_vix') is not None:
        return np.asarray(params['target_site_vix'], dtype=int)

    dists_to_anchor = _target_dists_to_anchor(mymesh, target_struct, params)
    candidates = np.where(dists_to_anchor < params['target_site_cutoff'])[0]
    if len(candidates) == 0:
        raise RuntimeError(
            f"No mesh vertices within {params['target_site_cutoff']} Å of target anchor. "
            "Check --target_chain, --target_resid, --target_site_cutoff."
        )

    n_keep = max(1, round(len(candidates) * params['target_site_sample_ratio']))
    if n_keep < len(candidates):
        np.random.seed(params['random_seed'])
        o3d_mesh = o3d.io.read_triangle_mesh(target_ply_fn)
        if not o3d_mesh.has_vertex_normals():
            o3d_mesh.compute_vertex_normals()
        target_vertices = _grid_subsample_vertices(o3d_mesh, candidates, n_keep, dists_to_anchor)
    else:
        target_vertices = candidates

    if params.get('target_site_max') and len(target_vertices) > params['target_site_max']:
        order = np.argsort(dists_to_anchor[target_vertices])
        target_vertices = target_vertices[order[:params['target_site_max']]]

    print(f"Selected {len(target_vertices)} target site(s).")
    return target_vertices


def _load_seed_subset_lines(subset_path):
    with open(subset_path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def _filter_seeds_for_auto_neosurf(seed_ids, seed_pdb_dir):
    kept = []
    anchor_cache = {}
    for protein_id in seed_ids:
        pdb_path = os.path.join(seed_pdb_dir, f"{protein_id}.pdb")
        anchor = find_ligand_anchor(pdb_path)
        if anchor is None:
            print(
                f"WARNING: skipped {protein_id} with no HETATOM record. "
                "If no ligand is expected, remove --seed-auto-neosurf flag."
            )
            continue
        log_ligand_anchor(protein_id, anchor, pdb_path)
        kept.append(protein_id)
        anchor_cache[protein_id] = anchor
    return kept, anchor_cache


def _resolve_target_auto_neosurf(args):
    pdb_path = os.path.join(args.target_pdb_dir, f"{args.target_name}.pdb")
    anchor = find_ligand_anchor(pdb_path)
    if anchor is None:
        raise SystemExit(
            f"ERROR: no hetero (HET) residue found in target PDB {pdb_path} "
            f"with --target-auto-neosurf."
        )
    log_ligand_anchor(args.target_name, anchor, pdb_path)
    args.target_chain = anchor.chain
    args.target_resid = anchor.resid
    args.target_residue = {'resid': anchor.resid, 'chain': anchor.chain}


def _validate_and_configure_cli(args, parser):
    has_vix = args.target_site_vix is not None or args.target_site_vix_file is not None
    has_residue = args.target_chain is not None and args.target_resid is not None
    has_point = args.target_coord is not None
    has_partial_residue = (args.target_chain is None) != (args.target_resid is None)

    if args.target_auto_neosurf:
        if has_residue or args.target_atom_id is not None:
            parser.error("--target-auto-neosurf cannot be combined with manual --target_chain/--target_resid/--target_atom_id")
        if has_vix:
            parser.error("--target-auto-neosurf cannot be combined with --target_site_vix/--target_site_vix_file")
        if has_point:
            parser.error("--target-auto-neosurf cannot be combined with --target_coord")
        if args.target_site_cutoff is None:
            parser.error("--target_site_cutoff is required with --target-auto-neosurf")

    if args.seed_auto_neosurf and args.seed_atom_id is not None:
        parser.error("--seed-auto-neosurf cannot be combined with --seed_atom_id")

    if has_partial_residue:
        parser.error("--target_chain and --target_resid must be used together")
    if not has_vix and not has_residue and not has_point and not args.target_auto_neosurf:
        parser.error(
            "Provide --target_site_vix / --target_site_vix_file, "
            "--target_chain + --target_resid, --target_coord, or --target-auto-neosurf"
        )
    if has_point and (has_residue or args.target_auto_neosurf):
        parser.error("--target_coord is mutually exclusive with residue-based target anchors")

    if has_residue and args.target_site_cutoff is None:
        parser.error("--target_site_cutoff is required when using --target_chain and --target_resid")

    has_manual_seed_anchor = args.seed_chain is not None and args.seed_resid is not None
    has_partial_seed_anchor = (args.seed_chain is None) != (args.seed_resid is None)
    if has_partial_seed_anchor:
        parser.error("--seed_chain and --seed_resid must be used together")
    if args.seed_auto_neosurf and (has_manual_seed_anchor or has_partial_seed_anchor):
        parser.error("--seed-auto-neosurf cannot be combined with --seed_chain/--seed_resid")
    if has_manual_seed_anchor and args.seed_site_cutoff is None:
        parser.error("--seed_site_cutoff is required with --seed_chain and --seed_resid")
    if args.seed_auto_neosurf and args.seed_site_cutoff is None:
        parser.error("--seed_site_cutoff is required with --seed-auto-neosurf")
    if args.seed_site_cutoff is not None and not args.seed_auto_neosurf and not has_manual_seed_anchor:
        parser.error(
            "Provide --seed_site_cutoff together with --seed-auto-neosurf "
            "or with --seed_chain and --seed_resid"
        )
    if args.seed_score_binder and (has_manual_seed_anchor or args.seed_auto_neosurf):
        parser.error(
            "--seed_score_binder cannot be combined with seed anchor options "
            "(--seed-auto-neosurf or --seed_chain/--seed_resid)"
        )

    if not (0 < args.target_site_sample_ratio <= 1.0):
        parser.error("--target_site_sample_ratio must be in (0, 1]")

    if args.target_coord is not None:
        args.target_point = {'coord': args.target_coord}
    elif has_residue:
        args.target_residue = {'resid': args.target_resid, 'chain': args.target_chain}
        if args.target_atom_id is not None:
            args.target_atom = {'atom_id': args.target_atom_id}

    if args.target_site_vix_file is not None and args.target_site_vix is None:
        args.target_site_vix = [
            int(line) for line in args.target_site_vix_file.read_text().strip().split('\n') if line.strip()
        ]


def score_complex(
    target_name,
    target_vertices,
    source_name,
    target_pcd,
    target_coord,
    target_desc,
    target_pcd_tree,
    target_iface,
    source_paths,
    params,
    nn_score,
    flip_target_normals=True,
):

    # Go through every selected site
    for site_ix, target_vix in enumerate(target_vertices):

        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
            target_pcd, target_coord, target_vix, target_desc,
            flip_normals=flip_target_normals,
            outward_shift=params['surface_outward_shift']
        )

        # Make a ckdtree with the target vertices.
        target_ckdtree = cKDTree(target_patch.points)

        # Load binder
        pid = 'p1'
        if pid == 'p1':
            chain = source_name.split('_')[1]
            chain_number = 1
        else:
            chain = source_name.split('_')[2]
            chain_number = 2
        source_pcd, source_desc, source_iface = load_protein_pcd(source_name, chain_number, source_paths, flipped_features=False, read_mesh=False)


        # Find closest point on binder
        target_vix_coord = np.asarray(target_pcd.points)[target_vix]
        source_points = np.asarray(source_pcd.points)
        dists = np.linalg.norm(target_vix_coord[None, :] - source_points, axis=-1)
        source_vix = np.argmin(dists)

        # Compute descriptor distance
        desc_dist = np.sqrt(np.sum(np.square(source_desc[0, source_vix] - target_desc[0, target_vix]), axis=-1))

        # Compute NN and descriptor distance scores
        source_coord = get_patch_coords(params['seed_precomp_dir'], source_name, pid, cv=[source_vix])
        source_patch, source_patch_descs, source_patch_idx = get_patch_geo(source_pcd, source_coord, source_vix, source_desc, outward_shift=params['surface_outward_shift'])
        d_vi_at, _= target_pcd_tree.query(np.asarray(source_patch.points), k=1)

        align_scores, _ = compute_nn_score(
            target_patch, source_patch, None, target_patch_descs,
            source_patch_descs, target_ckdtree, nn_score, d_vi_at, 1.0
        )

        results = {
            "query_name": target_name,
            "query_site": site_ix,
            "query_vix": target_vix,
            "binder_name": source_name,
            "binder_vix": source_vix,
            "distance_between_center_points": dists[source_vix],
            "query_iface_score": target_iface[target_vix],
            "binder_iface_score": source_iface[source_vix],
            "descriptor_distance": desc_dist,
            "nn_score": align_scores[0],
            "desc_dist_score": align_scores[1],
            "mean_desc_dist_score": align_scores[2],
        }
        print(', '.join(f"{k}: {v}" for k, v in results.items()))


def masif_search(params):

    # Decide whether to run complementarity or similarity search
    if params['similarity_mode']:
        params['surface_outward_shift'] = 0.0
        params['allowed_CA_clashes'] = float('inf')
        params['allowed_heavy_atom_clashes'] = float('inf')
        params['nn_score_cutoff'] = -1.0  # neural network score does not work here
        flip_target_features = False
        flip_target_normals = False
        print("Running surface similarity search. (experimental feature)")
    else:
        flip_target_features = True
        flip_target_normals = True
        print("Running surface complementarity search.")

    # # Load target patches.
    target_ppi_pair_id = params['target_name']
    target_pid = 'p1'
    target_chain_ix = 1

    seed_ppi_pair_ids = params['seed_ppi_pair_ids']

    # Initialize two neural networks - one that does not account for atomic clashes (initial filter) and one with clashes.
    nn_score_atomic = AlignmentEvaluationNN(params['nn_score_atomic_fn'], selected_features=[0,1,2,3], max_npoints=params['max_npoints'])
    nn_score_atomic.restore_model()

    # Go through every 12A patch in the target protein -- get a sorted least in order of the highest iface mean in the patch
    target_ply_fn = os.path.join(params['target_ply_iface_dir'], target_ppi_pair_id + '.ply')
    mymesh = pymesh.load_mesh(target_ply_fn)
    iface = mymesh.get_attribute('vertex_iface')
    target_coord = get_patch_coords(params['target_precomp_dir'], target_ppi_pair_id, target_pid)


    # Define target and source paths (for interface scores, descriptors, ply files)
    target_paths = {}
    target_paths['surf_dir'] = params['target_surf_dir']
    target_paths['iface_dir'] = params['target_iface_dir']
    target_paths['desc_dir'] = params['target_desc_dir']

    source_paths = {}
    source_paths['surf_dir'] = params['seed_surf_dir']
    source_paths['iface_dir'] = params['seed_iface_dir']
    source_paths['desc_dir'] = params['seed_desc_dir']

    # Load the target point cloud, descriptors, interface and mesh.
    target_pcd, target_desc, target_iface, target_mesh = load_protein_pcd(target_ppi_pair_id, target_chain_ix, target_paths, flipped_features=flip_target_features, read_mesh=True)

    # Open the pdb structure of the target, load into point clouds for fast access.
    parser = PDBParser()
    target_pdb_path = os.path.join(params['target_pdb_dir'], '{}.pdb'.format(target_ppi_pair_id))
    target_struct = parser.get_structure(target_pdb_path, target_pdb_path)
    target_atom_coords = [atom.get_coord() for atom in target_struct.get_atoms() if not atom.get_name().startswith('H') ]
    # Create kdtree search trees (for fast comparision).
    target_pcd_tree = cKDTree(np.array(target_atom_coords))

    # NOTE: target CA coords don't seem to be used
    # target_ca_coords = [atom.get_coord() for atom in target_struct.get_atoms() if atom.get_id() == 'CA']
    target_ca_pcd_tree = None # cKDTree(np.array(target_ca_coords))

    target_vertices = _select_target_vertices(mymesh, target_struct, target_ply_fn, params)

    if params.get('seed_score_binder', None) is not None:
        score_complex(
            target_name=target_ppi_pair_id,
            target_vertices=target_vertices,
            source_name=params['seed_score_binder'],
            target_pcd=target_pcd,
            target_coord=target_coord,
            target_desc=target_desc,
            target_pcd_tree=target_pcd_tree,
            target_iface=target_iface,
            source_paths=source_paths,
            params=params,
            nn_score=nn_score_atomic,
            flip_target_normals=True,
        )
        return  # descriptor matching and alignment are not necessary for scoring

    outdir = os.path.join(params['out_dir'], params['target_name'])
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Copy the pdb structure and the ply file of the target
    shutil.copy(target_pdb_path, outdir)
    shutil.copy(target_ply_fn, outdir)

    # Write out the targeted patches
    with open(os.path.join(outdir, 'selected_sites.vert'), 'w+') as out_patch:
        for site_vix in target_vertices:
            point = target_pcd.points[site_vix]
            out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))


    params['target_name'] = target_ppi_pair_id

    # Go through every selected site
    for site_ix, site_vix in enumerate(target_vertices):
        site_outdir = os.path.join(outdir, 'site_{}'.format(site_ix))
        if not os.path.exists(site_outdir):
            os.makedirs(site_outdir, exist_ok=True)

        params['target_site'] = 'site_{}'.format(site_ix)
        params['target_vix'] = int(site_vix)

        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
            target_pcd, target_coord, site_vix, target_desc,
            flip_normals=flip_target_normals,
            outward_shift=params['surface_outward_shift'])

        # Make a ckdtree with the target vertices.
        target_ckdtree = cKDTree(target_patch.points)

        # Write out the patch itself.
        with open(site_outdir + '/target.vert', 'w+') as out_patch:
            for point in target_patch.points:
                out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))

        with open(site_outdir + '/target_info.txt', 'w') as out_info:
            out_info.write(f'name: {target_ppi_pair_id}, site: {site_ix}, vix: {site_vix}')

        seeds_this_site = filter_seeds_for_resume(
            site_outdir, list(seed_ppi_pair_ids), params.get('resume', False))
        if len(seeds_this_site) == 0:
            print(f"Resume: all seeds already completed for {params['target_site']}, skipping site.")
            continue

        # Match the top descriptors in the database based on descriptor distance.
        print('Starting to match target descriptor to descriptors from {} proteins; this may take a while.'.format(len(seeds_this_site)))
        matched_dict, scores_dict = match_descriptors(seeds_this_site, ['p1', 'p2'], target_desc[0][site_vix], params, return_scores=True)

        matched_protein_ids = {name[0] for name in matched_dict.keys()}
        for seed_id in seeds_this_site:
            if seed_id not in matched_protein_ids:
                write_site_seed_hits(site_outdir, seed_id, [])
                print(f"No stage-1 match for {seed_id} at {params['target_site']}; "
                      f"wrote header-only {seed_id}.csv")

        if len(matched_dict.keys()) == 0:
            continue

        print(" ")
        print("Second stage of MaSIF seed search: each matched descriptor is aligned and scored; this may take a while..")
        count_matched_fragments = 0
        for ix, name in enumerate(matched_dict.keys()):
            try:
                align_protein(
                    name, \
                    target_patch, \
                    target_patch_descs, \
                    target_ckdtree, \
                    target_ca_pcd_tree, \
                    target_pcd_tree, \
                    source_paths, \
                    matched_dict,\
                    nn_score_atomic, \
                    site_outdir, \
                    params, \
                    first_stage_scores=scores_dict[name],
                    n_retry_alignment=params['n_retry_alignment'],
                )
            except Exception as e:
                print(f"Error '{e}' while trying to align {name}.")
            if (ix + 1) % 1000 == 0:
                print('So far, MaSIF has aligned {} fragments from {} proteins.'.format(count_matched_fragments, ix + 1))
            count_matched_fragments += len(matched_dict[name])

    print("Done!")


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--target", dest="target_name", type=str, required=True)
    parser.add_argument("--target_dir", dest="masif_target_root", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, default=None)

    parser.add_argument("--target-auto-neosurf", action="store_true",
                        help="Use largest hetero residue in target PDB as anchor.")
    parser.add_argument("--target_chain", type=str, default=None)
    parser.add_argument("--target_resid", type=int, default=None)
    parser.add_argument("--target_atom_id", type=str, default=None)
    parser.add_argument("--target_coord", type=float, nargs="+", default=None)
    parser.add_argument("--target_site_cutoff", type=float, default=None)
    parser.add_argument("--target_site_sample_ratio", type=float, default=1.0)
    parser.add_argument("--target_site_max", type=int, default=None)
    parser.add_argument("--target_site_vix", type=int, nargs='+', default=None)
    parser.add_argument("--target_site_vix_file", type=Path, default=None)

    parser.add_argument("--seed_dir", dest="top_seed_dir", type=Path, required=True)
    parser.add_argument("--seed_subset", dest="database_subset", type=Path, default=None)
    parser.add_argument("--seed-auto-neosurf", action="store_true",
                        help="Use largest hetero residue per seed PDB as anchor for spatial filtering.")
    parser.add_argument("--seed_chain", type=str, default=None)
    parser.add_argument("--seed_resid", type=int, default=None)
    parser.add_argument("--seed_atom_id", type=str, default=None)
    parser.add_argument("--seed_site_cutoff", type=float, default=None)
    parser.add_argument("--seed_desc_dist_cutoff", type=float, default=2.0,
                        help="Recommended values: [1.5-2.0] (lower is stricter)")
    parser.add_argument("--seed_iface_cutoff", type=float, default=0.75,
                        help="Recommended values: [0.75-0.95] (higher is stricter)")
    parser.add_argument("--seed_nn_score_cutoff", type=float, default=0.8,
                        help="Recommended values: [0.8-0.95] (higher is stricter)")
    parser.add_argument("--seed_desc_dist_score_cutoff", type=float, default=0.0,
                        help="Recommended values: [0.0-20.0] (higher is stricter)")
    parser.add_argument("--seed_score_binder", type=str, default=None,
                        help="Score one processed seed protein without database search.")

    parser.add_argument("--allowed_CA_clashes", type=int, default=0)
    parser.add_argument("--allowed_heavy_atom_clashes", type=int, default=5)
    parser.add_argument("--ransac_iter", type=int, default=100000)
    parser.add_argument("--n_retry_alignment", type=int, default=1)
    parser.add_argument("--similarity_mode", action="store_true")
    parser.add_argument("--resume", action="store_true",
                        help="Skip (site, seed) pairs that already have a completed CSV output.")
    parser.add_argument("--random_seed", type=int, default=42)
    args = parser.parse_args()

    _validate_and_configure_cli(args, parser)

    np.random.seed(args.random_seed)
    if version.parse('0.14.1') <= version.parse(o3d.__version__) <= version.parse('0.15.2'):
        args.maybe_seed = {"seed": args.random_seed}
    else:
        args.maybe_seed = {}

    # Keys consumed by masif_search() and alignment_utils (legacy names)
    args.desc_dist_cutoff = args.seed_desc_dist_cutoff
    args.iface_cutoff = args.seed_iface_cutoff
    args.nn_score_cutoff = args.seed_nn_score_cutoff
    args.desc_dist_score_cutoff = args.seed_desc_dist_score_cutoff

    # Database locations
    args.seed_surf_dir = os.path.join(args.top_seed_dir, masif_opts['ply_chain_dir'])
    args.seed_iface_dir = os.path.join(args.top_seed_dir, masif_opts['site']['out_pred_dir'])
    args.seed_ply_iface_dir = os.path.join(args.top_seed_dir, masif_opts['site']['out_surf_dir'])
    args.seed_pdb_dir = os.path.join(args.top_seed_dir, masif_opts['pdb_chain_dir'])
    args.seed_desc_dir = os.path.join(args.top_seed_dir, masif_opts['ppi_search']['desc_dir'])

    args.seed_precomp_dir = os.path.join(args.top_seed_dir, masif_opts['ppi_search']['masif_precomputation_dir'])

    # Target locations
    args.top_target_dir = os.path.join(args.masif_target_root)
    args.target_surf_dir = os.path.join(args.top_target_dir, masif_opts['ply_chain_dir'])
    args.target_iface_dir = os.path.join(args.masif_target_root, masif_opts['site']['out_pred_dir'])
    args.target_ply_iface_dir = os.path.join(args.masif_target_root, masif_opts['site']['out_surf_dir'])
    args.target_pdb_dir = os.path.join(args.top_target_dir, masif_opts['pdb_chain_dir'])
    args.target_desc_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['desc_dir'])

    args.target_precomp_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['masif_precomputation_dir'])

    if args.target_auto_neosurf:
        _resolve_target_auto_neosurf(args)

    # Database subset
    subset_lines = None
    if args.database_subset is not None:
        if not args.database_subset.is_file():
            parser.error(f"--seed_subset not found: {args.database_subset}")
        subset_lines = _load_seed_subset_lines(args.database_subset)
        if len(subset_lines) == 0:
            print(f"NOTE: empty seed subset file {args.database_subset}; nothing to search.")
            sys.exit(0)
        if args.seed_auto_neosurf:
            args.seed_ppi_pair_ids, args.seed_anchor_cache = _filter_seeds_for_auto_neosurf(
                subset_lines, args.seed_pdb_dir)
            if len(args.seed_ppi_pair_ids) == 0:
                raise SystemExit(
                    f"ERROR: all {len(subset_lines)} seed(s) in {args.database_subset} "
                    "were skipped (no hetero residue)."
                )
        else:
            args.seed_ppi_pair_ids = subset_lines
    elif args.seed_auto_neosurf:
        args.seed_ppi_pair_ids = np.array(os.listdir(args.seed_desc_dir))
        args.seed_anchor_cache = {}
    else:
        args.seed_ppi_pair_ids = np.array(os.listdir(args.seed_desc_dir))

    # Some hard-coded parameters
    args.nn_score_atomic_fn = os.path.join(masif_neosurf_dir, "masif_seed_search/data/scoring_nn/models_std/weights_12A_0129")
    args.max_npoints = 200

    if version.parse(o3d.__version__) > version.parse('0.11.0'):
        args.ransac_convergence_kwargs = {'confidence': 0.999}
    args.ransac_radius = 1.5
    args.surface_outward_shift = 0.25

    # Run main search function
    search_params = vars(args)
    masif_search(search_params)

    # Cluster search hits upon completion of search
    if search_params.get("seed_score_binder") is None:
        cluster_search_hits_for_subset(search_params)
