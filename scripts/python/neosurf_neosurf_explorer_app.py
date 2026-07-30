"""
Neosurf-on-Neosurf Results Explorer (Dash)
==========================================

Interactive explorer for the all-vs-all MaSIF-Neosurf search produced by
``notebooks/neosurf_neosurf.ipynb``.

Unlike the single-target legacy app (``results_explorer_app.py``), this workflow
queries *many* ligand-bound "target" complexes against a large library of
ligand-bound "seed" complexes, so every predicted interaction row has a
different target AND a different seed, each carrying its own ligand.

Data sources
  * ``df_preprocess_ok.csv`` -> per-structure metadata (uniprot/gene/name,
    pdb_path, ligand identity), keyed by the unique ``target`` id. Source of
    metadata for BOTH targets and seeds. ``pdb_path`` values are rewritten at
    load time to ``data/preprocess/data_preparation/01-benchmark_pdbs/``.
  * ``df_results_dedup.csv`` -> PRIMARY source (``4_enrich_metrics``). One
    enriched row per (target, matched_protein, cluster_id) cluster; drives the
    tree and the colored scatter points. Includes ``target_mw``,
    ``total_n_patches``, and cluster_size normalisations.
  * ``df_results_all.csv``   -> every individual pose (unique target_vix +
    matched_patch_id within a cluster), also from ``4_enrich_metrics``. Loaded
    once into memory and indexed by cluster key. Poses are populated reactively:
    only when the user expands a cluster in the tree. Poses inherit all metadata
    from their cluster; only the superposition transform / patch identity /
    score differ per pose.

Panels
  * Left sidebar: a lazily-rendered, re-orderable collapsible tree
    (target-protein -> target-id -> matched-protein -> matched-seed -> cluster
    -> pose). Clusters are expandable to reveal individual poses.
  * Plotly scatter. By default only the selected tree branch is plotted (e.g. a
    single target complex). Optionally enable "show all entries in this branch"
    in the sidebar to keep the top-level scope visible with non-selected points
    as grey crosses. When a cluster is pinned, that cluster's individual poses
    are overlaid as grey crosses so they can be clicked and visualized.
  * py3Dmol viewer showing target + transform-superposed seed, both ligands as
    sticks.
  * UniProt entry iframe for the matched (or target) protein.
  * Flag / unflag the currently selected dedup cluster (pose selections map to
    their parent cluster). Flagged rows live in an in-memory store, are marked
    on the scatter, and can be exported to a user-specified CSV with the same
    columns as ``df_results_dedup.csv`` plus a generated ``pymol_cmd`` column
    (fetch / transform / align sequence matching the notebook visualiser).

The legacy app is intentionally left untouched.
"""

import copy
import json
import os
import uuid

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
from Bio.PDB import PDBIO, PDBParser
from dash import ALL, Dash, Input, Output, State, callback_context, dcc, html, no_update
from dash.exceptions import PreventUpdate

# --------------------------------------------------------------------------- #
# Data paths
# --------------------------------------------------------------------------- #
PREPROCESS_CSV = "/scratch/ymeng/Neosurf_Neosurf/data/processing/2_masif_preprocess/df_preprocess_ok.csv"
RESULTS_CSV = "/scratch/ymeng/Neosurf_Neosurf/data/processing/4_enrich_metrics/df_results_dedup.csv"
RESULTS_ALL_CSV = "/scratch/ymeng/Neosurf_Neosurf/data/processing/4_enrich_metrics/df_results_all.csv"
DEFAULT_FLAGGED_CSV = "/scratch/ymeng/Neosurf_Neosurf/data/processing/4_enrich_metrics/df_flagged.csv"
# Preprocess CSV still records the old ``data/input/`` location; PDBs now live here.
BENCHMARK_PDB_DIR = "/scratch/ymeng/Neosurf_Neosurf/data/preprocess/data_preparation/01-benchmark_pdbs"
# Original on-disk columns of RESULTS_CSV (enriched cols like g_* are excluded on save).
DEDUP_CSV_COLS = list(pd.read_csv(RESULTS_CSV, nrows=0).columns)

# Max rows a selected branch's scope may contain before we refuse to plot
# (guardrail so Plotly is never handed tens of thousands of interactive points).
PLOT_ROW_CAP = 4000

# --------------------------------------------------------------------------- #
# Grouping levels for the tree
# --------------------------------------------------------------------------- #
# The four "protein/id" levels are user-reorderable. ``cluster_id`` always
# follows them and is expandable; ``pose`` is the final leaf, sourced lazily
# from df_results_all.
NONLEAF_LEVELS = ["target_protein", "target_id", "matched_protein", "matched_seed_id"]
CLUSTER_LEVEL = "cluster_id"
POSE_LEVEL = "pose"

LEVEL_LABEL = {
    "target_protein": "Target protein",
    "target_id": "Target complex (id)",
    "matched_protein": "Matched protein",
    "matched_seed_id": "Matched seed (id)",
    "cluster_id": "Cluster",
    "pose": "Pose",
}

# Grouping column backing each df-level (string columns built at load time).
# ``pose`` has no df grouping column; it is resolved against df_all by cluster key.
GROUP_COL = {
    "target_protein": "g_target_protein",
    "target_id": "g_target_id",
    "matched_protein": "g_matched_protein",
    "matched_seed_id": "g_matched_seed_id",
    "cluster_id": "g_cluster_id",
}
DF_LEVELS = set(GROUP_COL)

DEFAULT_ORDER = ["target_protein", "target_id", "matched_protein", "matched_seed_id", CLUSTER_LEVEL, POSE_LEVEL]
PRESET_A = ["target_protein", "target_id", "matched_protein", "matched_seed_id"]
PRESET_B = ["target_protein", "matched_protein", "target_id", "matched_seed_id"]

# Metric columns offered as plot axes / colours.
CANDIDATE_NUMERIC = [
    "score", "desc_dist_score", "iface_score", "mean_desc_dist_score", "desc_dist",
    "cluster_size", "cluster_mean_rmsd", "clashing_ca", "clashing_heavy",
    "n_match_ligands", "matched_mw", "target_mw", "total_n_patches",
    "cluster_size_patch_normalized", "cluster_size_mw_normalized",
]


# --------------------------------------------------------------------------- #
# Structure helpers (shared logic with the legacy app)
# --------------------------------------------------------------------------- #
def apply_transform(structure, transform_str):
    """Apply a flattened 4x4 affine transform to every atom of a structure."""
    transform_flat = np.fromstring(transform_str, sep=",")
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


def structure_to_pdbstr(structure):
    """Serialise a Bio.PDB structure to a PDB string."""
    from io import StringIO

    io = PDBIO()
    io.set_structure(structure)
    buffer = StringIO()
    io.save(buffer)
    return buffer.getvalue()


NON_LIGAND_HET_RESNAMES = frozenset(
    {"HOH", "WAT", "H2O", "DOD", "NA", "CL", "K", "MG", "ZN", "CA", "MN", "FE", "CO", "NI", "CU", "CD"}
)


def find_hetatm_residues(structure, exclude_non_ligands=True):
    """Return HETATM residues (chain/resn/resi) - fallback ligand identification."""
    residues, seen = [], set()
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                if hetflag == " ":
                    continue
                resn = residue.resname.strip()
                if exclude_non_ligands and resn in NON_LIGAND_HET_RESNAMES:
                    continue
                key = (chain.id, resn, resseq, icode.strip())
                if key in seen:
                    continue
                seen.add(key)
                residues.append({"chain": chain.id, "resn": resn, "resi": str(resseq)})
    return residues


# --------------------------------------------------------------------------- #
# Load & enrich data
# --------------------------------------------------------------------------- #
def _s(value):
    """NaN-safe string for label building."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value)


def _cluster_id_str(value):
    if pd.isna(value):
        return ""
    return str(int(value)) if float(value).is_integer() else str(value)


def _add_structural(frame, pre_idx):
    """Attach pdb paths + ligand identity for both target and seed sides."""
    def mp(col):
        return pre_idx[col] if col in pre_idx.columns else pd.Series(dtype=object)

    frame["t_uniprot_id"] = frame["target"].map(mp("uniprot_id"))
    frame["t_recommendedName"] = frame["target"].map(mp("recommendedName"))
    frame["target_pdb_id"] = frame["target"].map(mp("pdb_id"))
    frame["matched_pdb_id"] = frame["matched_protein"].map(mp("pdb_id"))
    frame["target_pdb_path"] = frame["target"].map(mp("pdb_path"))
    frame["target_ligand_resname"] = frame["target"].map(mp("ligand_resname"))
    frame["target_pdb_ligand_chain"] = frame["target"].map(mp("pdb_ligand_chain"))
    frame["target_ligand_code"] = frame["target"].map(mp("ligand_code"))

    frame["matched_pdb_path"] = frame["matched_protein"].map(mp("pdb_path"))
    frame["matched_ligand_resname"] = frame["matched_protein"].map(mp("ligand_resname"))
    frame["matched_pdb_ligand_chain"] = frame["matched_protein"].map(mp("pdb_ligand_chain"))
    frame["matched_ligand_code"] = frame["matched_protein"].map(mp("ligand_code"))
    return frame


def _add_grouping_cols(frame):
    frame["g_target_protein"] = frame["t_uniprot_id"].fillna("—").astype(str)
    frame["g_target_id"] = frame["target"].astype(str)
    frame["g_matched_protein"] = frame["matched_uniprot_id"].fillna("—").astype(str)
    frame["g_matched_seed_id"] = frame["matched_protein"].astype(str)
    frame["g_cluster_id"] = frame["cluster_id"].map(_cluster_id_str)
    return frame


def load_data():
    """Load all three tables and build the enriched frames used throughout."""
    pre = pd.read_csv(PREPROCESS_CSV)
    # Rewrite stale ``data/input/`` paths to the current benchmark PDB directory.
    if "pdb_path" in pre.columns:
        pre["pdb_path"] = pre["pdb_path"].map(
            lambda p: os.path.join(BENCHMARK_PDB_DIR, os.path.basename(str(p)))
            if pd.notna(p) and str(p).strip()
            else p
        )
    pre_idx = pre.set_index("target")

    # ---- Primary (dedup) frame: authoritative, already enriched in the CSV ----
    df = pd.read_csv(RESULTS_CSV).reset_index(drop=True)
    df = _add_structural(df, pre_idx)
    df = _add_grouping_cols(df)

    search_cols = [
        "target_gene_name", "t_uniprot_id", "t_recommendedName",
        "matched_gene_name", "matched_uniprot_id", "matched_recommendedName",
        "target", "matched_protein",
    ]
    search_cols = [c for c in search_cols if c in df.columns]
    df["search_blob"] = df[search_cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower()

    # ---- All-poses frame: raw search output, enriched from the manifest ----
    df_all = pd.read_csv(RESULTS_ALL_CSV).reset_index(drop=True)
    df_all = _add_structural(df_all, pre_idx)
    # Metadata that lived only in the (notebook-processed) dedup CSV -> derive.
    df_all["target_gene_name"] = df_all["target"].map(
        pre_idx["gene_name"] if "gene_name" in pre_idx.columns else pd.Series(dtype=object)
    )
    df_all["matched_uniprot_id"] = df_all["matched_protein"].map(
        pre_idx["uniprot_id"] if "uniprot_id" in pre_idx.columns else pd.Series(dtype=object)
    )
    df_all["matched_gene_name"] = df_all["matched_protein"].map(
        pre_idx["gene_name"] if "gene_name" in pre_idx.columns else pd.Series(dtype=object)
    )
    df_all["matched_recommendedName"] = df_all["matched_protein"].map(
        pre_idx["recommendedName"] if "recommendedName" in pre_idx.columns else pd.Series(dtype=object)
    )
    df_all["matched_mw"] = df_all["matched_protein"].map(
        pre_idx["mw"] if "mw" in pre_idx.columns else pd.Series(dtype=object)
    )
    # n_match_ligands is a per gene-pair value computed in the notebook; borrow it
    # from the dedup frame so poses share their cluster's value (grey points only).
    if "n_match_ligands" in df.columns:
        nml = (
            df.dropna(subset=["n_match_ligands"])
            .drop_duplicates(["target_gene_name", "matched_gene_name"])
            .set_index(["target_gene_name", "matched_gene_name"])["n_match_ligands"]
            .to_dict()
        )
        df_all["n_match_ligands"] = [
            nml.get((tg, mg), np.nan) for tg, mg in zip(df_all["target_gene_name"], df_all["matched_gene_name"])
        ]
    df_all["g_cluster_id"] = df_all["cluster_id"].map(_cluster_id_str)

    # Mark the representative pose (max score) of each cluster.
    df_all["is_rep"] = False
    rep_idx = df_all.groupby(["target", "matched_protein", "g_cluster_id"])["score"].idxmax()
    df_all.loc[rep_idx, "is_rep"] = True

    return pre, df, df_all


pre_df, df, df_all = load_data()

# target id -> {pdb_id, pdb_protein_chain, pdb_ligand_chain} for PyMOL cmds.
_TARGET_INFO = pre_df.set_index("target")[
    ["pdb_id", "pdb_protein_chain", "pdb_ligand_chain"]
].to_dict("index")

# Cluster-key -> pose row indices, precomputed once for O(1) lazy expansion.
POSE_GROUPS = {k: list(v) for k, v in df_all.groupby(["target", "matched_protein", "g_cluster_id"]).groups.items()}
POSE_COUNTS = {k: len(v) for k, v in POSE_GROUPS.items()}

NUMERIC_COLS = [c for c in CANDIDATE_NUMERIC if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]
DEFAULT_X = "score" if "score" in NUMERIC_COLS else NUMERIC_COLS[0]
DEFAULT_Y = "cluster_size" if "cluster_size" in NUMERIC_COLS else NUMERIC_COLS[min(1, len(NUMERIC_COLS) - 1)]
DEFAULT_COLOR = "n_match_ligands" if "n_match_ligands" in NUMERIC_COLS else NUMERIC_COLS[0]

_parser = PDBParser(QUIET=True)


# --------------------------------------------------------------------------- #
# Tree / path logic
# --------------------------------------------------------------------------- #
def subframe_by_path(base, path):
    """Filter ``base`` to rows matching every df-level (level, value) pair.

    Any trailing pose level is ignored (poses do not live in the df frames).
    """
    mask = pd.Series(True, index=base.index)
    for level, value in path:
        if level not in GROUP_COL:
            continue
        mask &= base[GROUP_COL[level]] == value
    return base[mask]


def strip_pose(path):
    """Return the path without a trailing pose entry."""
    return [p for p in path if p[0] != POSE_LEVEL]


def cluster_key_from_path(path):
    """Extract (target, matched_protein, cluster_id_str) from a path prefix."""
    values = {level: value for level, value in path}
    return values.get("target_id"), values.get("matched_seed_id"), values.get("cluster_id")


def node_label(level, value, rep):
    """Human label for a df-level node, given a representative row ``rep``."""
    if level == "target_protein":
        return f"{_s(rep.get('t_uniprot_id'))}|{_s(rep.get('target_gene_name'))}|{_s(rep.get('t_recommendedName'))}"
    if level == "target_id":
        code = _s(rep.get("target_ligand_code"))
        return f"target: {value}" + (f"  [{code}]" if code else "")
    if level == "matched_protein":
        return (
            f"{_s(rep.get('matched_uniprot_id'))}|{_s(rep.get('matched_gene_name'))}|"
            f"{_s(rep.get('matched_recommendedName'))}"
        )
    if level == "matched_seed_id":
        code = _s(rep.get("matched_ligand_code"))
        return f"seed: {value}" + (f"  [{code}]" if code else "")
    if level == "cluster_id":
        size, score = rep.get("cluster_size"), rep.get("score")
        size_s = f"size={int(size)}" if pd.notna(size) else "size=?"
        score_s = f"score={score:.3f}" if pd.notna(score) else "score=?"
        return f"cluster {value} ({size_s}, {score_s})"
    return str(value)


def pose_label(prow):
    """Human label for a pose (df_all row)."""
    patch = prow.get("matched_patch_id")
    vix = prow.get("target_vix")
    score = prow.get("score")
    patch_s = int(patch) if pd.notna(patch) else "?"
    vix_s = int(vix) if pd.notna(vix) else "?"
    score_s = f"{score:.3f}" if pd.notna(score) else "?"
    star = "  ★rep" if prow.get("is_rep") else ""
    return f"pose patch={patch_s} vix={vix_s} (score={score_s}){star}"


def _node_style(depth, is_leaf, is_selected):
    return {
        "display": "block",
        "width": "100%",
        "textAlign": "left",
        "border": "none",
        "background": "#eef4ff" if is_selected else "transparent",
        "padding": "2px 4px",
        "marginLeft": f"{depth * 12}px",
        "cursor": "pointer",
        "fontSize": "12px" if not is_leaf else "11px",
        "color": "#333" if not is_leaf else "#555",
        "fontFamily": "monospace",
        "whiteSpace": "nowrap",
        "overflow": "hidden",
        "textOverflow": "ellipsis",
    }


def _pose_nodes(cluster_path, expanded, selected, depth, criteria=None):
    """Leaf nodes for the individual poses of a cluster (from df_all)."""
    key = cluster_key_from_path(cluster_path)
    idxs = POSE_GROUPS.get(key)
    if not idxs:
        return [html.Div("(no individual poses in df_results_all)",
                         style={"marginLeft": f"{depth * 12}px", "fontSize": "11px", "color": "#999"})]
    sub = df_all.loc[idxs]
    if criteria:
        sub = apply_filtering_criteria(sub, criteria)
    if sub.empty:
        return [html.Div("(all poses filtered out)",
                         style={"marginLeft": f"{depth * 12}px", "fontSize": "11px", "color": "#999"})]
    sub = sub.sort_values("score", ascending=False)
    nodes = []
    for all_idx, prow in sub.iterrows():
        child_path = cluster_path + [[POSE_LEVEL, str(all_idx)]]
        path_str = json.dumps(child_path)
        button = html.Button(
            f"• {pose_label(prow)}",
            id={"type": "tree-node", "path": path_str},
            n_clicks=0, title=pose_label(prow),
            style=_node_style(depth, True, path_str == selected),
        )
        nodes.append(html.Div([button]))
    return nodes


def build_tree_nodes(path, order, expanded, selected, base, depth=0, criteria=None):
    """Recursively build tree nodes, descending only into expanded paths.

    Sibling nodes are ordered by their number of immediate child nodes
    (descending), with best score as a tie-breaker.
    """
    level_idx = len(path)
    if level_idx >= len(order):
        return []
    level = order[level_idx]

    # Pose level: children come from df_all, keyed by the parent cluster.
    if level == POSE_LEVEL:
        return _pose_nodes(path, expanded, selected, depth, criteria)

    col = GROUP_COL[level]
    sub = subframe_by_path(base, path)
    if sub.empty:
        return []

    next_level = order[level_idx + 1] if level_idx + 1 < len(order) else None
    next_col = GROUP_COL.get(next_level) if next_level else None

    def child_count(value, ssub):
        """Number of immediate child nodes this node would expand into."""
        if level == CLUSTER_LEVEL:  # children are poses, sourced from df_all
            pidx = POSE_GROUPS.get(cluster_key_from_path(path + [[level, value]]))
            if not pidx:
                return 0
            return len(apply_filtering_criteria(df_all.loc[pidx], criteria)) if criteria else len(pidx)
        if next_col:
            return int(ssub[next_col].nunique())
        return len(ssub)

    # One entry per distinct value at this level, with its child count + best score.
    reps = []
    for value, idx in sub.groupby(col, sort=False).groups.items():
        ssub = sub.loc[idx]
        reps.append((value, child_count(value, ssub), base.loc[ssub["score"].idxmax()], float(ssub["score"].max())))
    reps.sort(key=lambda t: (t[1], t[3]), reverse=True)  # most children first, then best score

    nodes = []
    for value, n_children, rep, best in reps:
        child_path = path + [[level, value]]
        path_str = json.dumps(child_path)
        is_expanded = path_str in expanded
        marker = "▾" if is_expanded else "▸"
        label = node_label(level, value, rep)
        badge = f"  ({n_children} poses)" if level == CLUSTER_LEVEL else f"  ({n_children})"
        button = html.Button(
            f"{marker} {label}{badge}",
            id={"type": "tree-node", "path": path_str},
            n_clicks=0, title=label,
            style=_node_style(depth, False, path_str == selected),
        )
        block = [button]
        if is_expanded:
            block.append(html.Div(build_tree_nodes(child_path, order, expanded, selected, base, depth + 1, criteria)))
        nodes.append(html.Div(block))
    return nodes


def build_order(vals):
    """Normalise the 4 dropdown values into a valid permutation + cluster + pose."""
    seen = []
    for v in vals:
        if v and v not in seen:
            seen.append(v)
    for v in NONLEAF_LEVELS:
        if v not in seen:
            seen.append(v)
    return seen + [CLUSTER_LEVEL, POSE_LEVEL]


# --------------------------------------------------------------------------- #
# Row reference resolution (dedup rows and pose rows share a tagged reference)
# --------------------------------------------------------------------------- #
def resolve_ref(ref):
    """Resolve a {"src","idx"} reference to a row Series, or None."""
    if not isinstance(ref, dict):
        return None
    src, idx = ref.get("src"), ref.get("idx")
    frame = df if src == "dedup" else df_all
    if idx in frame.index:
        return frame.loc[idx]
    return None


def dedup_idx_from_ref(ref):
    """Map a selected-row ref to a dedup ``df`` index, or None.

    Pose refs resolve to the parent cluster's representative dedup row via
    ``(target, matched_protein, g_cluster_id)``.
    """
    if not isinstance(ref, dict):
        return None
    src, idx = ref.get("src"), ref.get("idx")
    if src == "dedup":
        return int(idx) if idx in df.index else None
    if src != "pose" or idx not in df_all.index:
        return None
    prow = df_all.loc[idx]
    hit = df[
        (df["target"] == prow["target"])
        & (df["matched_protein"] == prow["matched_protein"])
        & (df["g_cluster_id"] == prow["g_cluster_id"])
    ]
    if hit.empty:
        return None
    return int(hit.index[0])


def _get_target_info(target_id):
    """Return (pdb_id, chain_list) for a target/seed id via the preprocess manifest.

    Falls back to splitting on the last underscore when the id is missing from
    the manifest (same convention as the notebook visualiser).
    """
    if target_id in _TARGET_INFO:
        info = _TARGET_INFO[target_id]
        chains = [info["pdb_protein_chain"]]
        if info["pdb_ligand_chain"] != info["pdb_protein_chain"]:
            chains.append(info["pdb_ligand_chain"])
        return info["pdb_id"], chains
    parts = str(target_id).rsplit("_", 1)
    chains = list(parts[1]) if len(parts) > 1 else []
    return parts[0], chains


def _fmt_patch_id(value):
    if pd.isna(value):
        return "?"
    f = float(value)
    return str(int(f)) if f.is_integer() else str(value)


def pymol_cmd_for_row(row):
    """Build a single-line PyMOL command sequence for one dedup result row.

    Matches the format in ``data/df_results_dedup_pymol.csv``: fetch target +
    matched chains, apply the flattened transform, merge into one object, then
    clean up and align/orient to ``reference``.
    """
    target = row["target"]
    matched = row["matched_protein"]
    patch = _fmt_patch_id(row["matched_patch_id"])
    transform = row["flattened_transform"]

    target_pdb, target_chains = _get_target_info(target)
    matched_pdb, matched_chains = _get_target_info(matched)
    target_chains_str = "chain " + " chain ".join(str(c) for c in target_chains)
    matched_chains_str = "chain " + " chain ".join(str(c) for c in matched_chains)

    complex_obj = f"{target}_{matched}_{patch}"
    matched_obj = f"{matched}_{patch}"

    return "; ".join([
        f"fetch {target_pdb}, {complex_obj}",
        f"remove {complex_obj} AND (not {target_chains_str})",
        f"fetch {matched_pdb}, {matched_obj}",
        f"remove {matched_obj} AND (not {matched_chains_str})",
        f"apply_transform {matched_obj}, '{transform}'",
        f"copy_to {complex_obj}, {matched_obj}",
        f"delete {matched_obj}",
        "remove (hydro)",
        "remove resn hoh",
        "show lines",
        f"align {complex_obj}, reference",
        f"orient {complex_obj}",
    ])


def flagged_frame(idxs):
    """Return flagged dedup rows with original CSV columns plus ``pymol_cmd``."""
    cols = DEDUP_CSV_COLS + ["pymol_cmd"]
    if not idxs:
        return pd.DataFrame(columns=cols)
    out = df.loc[list(idxs), DEDUP_CSV_COLS].copy()
    out["pymol_cmd"] = out.apply(pymol_cmd_for_row, axis=1)
    return out


# --------------------------------------------------------------------------- #
# 3D viewer (dual ligand)
# --------------------------------------------------------------------------- #
def _ligand_selection(model, resname, chain):
    sel = {"model": model}
    if isinstance(resname, str) and resname.strip():
        sel["resn"] = resname.strip()
    if isinstance(chain, str) and chain.strip():
        sel["chain"] = chain.strip()
    return sel


def _add_seed_model(view, seed_struct, model_idx, resname, chain, opacity):
    """Add one transform-applied seed structure to the view at a given opacity."""
    view.addModel(structure_to_pdbstr(seed_struct), "pdb")
    cartoon = {"color": "lightblue"}
    line = {"colorscheme": "lightblueCarbon"}
    stick = {"colorscheme": "magentaCarbon", "radius": 0.22}
    if opacity < 1.0:
        cartoon["opacity"] = opacity
        line["opacity"] = opacity
        stick["opacity"] = opacity
    view.setStyle({"model": model_idx}, {"cartoon": cartoon, "line": line})
    seed_lig = _ligand_selection(model_idx, resname, chain)
    if "resn" in seed_lig:
        view.setStyle(seed_lig, {"stick": stick})
    else:
        for het in find_hetatm_residues(seed_struct):
            sel = {"model": model_idx, "resn": het["resn"], "resi": het["resi"]}
            if het["chain"]:
                sel["chain"] = het["chain"]
            view.setStyle(sel, {"stick": stick})


def build_viewer_html(row, pose_transforms):
    """Return py3Dmol HTML for the target + one-or-more transform-superposed seed poses.

    ``row`` supplies the target + seed metadata (paths, ligands) shared by every
    pose of a cluster. ``pose_transforms`` is a list of (transform_str, opacity)
    for the seed: opacity 1.0 = the solid main entry, < 1.0 = a faded sibling pose.
    """
    import os

    target_path = row.get("target_pdb_path")
    seed_path = row.get("matched_pdb_path")
    for label, p in (("target", target_path), ("seed", seed_path)):
        if not isinstance(p, str) or not p.strip():
            raise ValueError(f"Missing {label} pdb_path for selected row.")
        if not os.path.exists(p):
            raise FileNotFoundError(f"{label} PDB not found: {p}")
    if not pose_transforms:
        raise ValueError("No pose transforms to render.")

    target_struct = _parser.get_structure("target", target_path)

    view = py3Dmol.view(width=760, height=560)
    view.addModel(structure_to_pdbstr(target_struct), "pdb")  # model 0 = target (solid)
    view.setStyle({"model": 0}, {"cartoon": {"color": "lightgrey"}, "line": {"colorscheme": "lightgreyCarbon"}})

    # Target ligand -> green sticks.
    target_lig = _ligand_selection(0, row.get("target_ligand_resname"), row.get("target_pdb_ligand_chain"))
    if "resn" in target_lig:
        view.setStyle(target_lig, {"stick": {"colorscheme": "greenCarbon", "radius": 0.22}})
    else:
        for het in find_hetatm_residues(target_struct):
            sel = {"model": 0, "resn": het["resn"], "resi": het["resi"]}
            if het["chain"]:
                sel["chain"] = het["chain"]
            view.setStyle(sel, {"stick": {"colorscheme": "greenCarbon", "radius": 0.22}})

    # Seed poses -> models 1..N (main solid, siblings faded).
    for i, (transform, opacity) in enumerate(pose_transforms):
        if not isinstance(transform, str) or not transform.strip():
            continue
        seed_struct = _parser.get_structure(f"seed{i}", seed_path)
        apply_transform(seed_struct, transform)
        _add_seed_model(view, seed_struct, i + 1, row.get("matched_ligand_resname"),
                        row.get("matched_pdb_ligand_chain"), opacity)

    if "resn" in target_lig:
        view.zoomTo(target_lig)
    else:
        view.zoomTo()
    return view._make_html()


def view_pair_html(row):
    """Single-pose (solid) view for the selected row."""
    return build_viewer_html(row, [(row.get("flattened_transform"), 1.0)])


def cluster_key_from_row(row):
    """(target, matched_protein, cluster_id_str) key matching POSE_GROUPS."""
    return (str(row.get("target")), str(row.get("matched_protein")), _cluster_id_str(row.get("cluster_id")))


def cluster_pose_transforms(row, selected_ref):
    """Build (transform, opacity) list for every pose of ``row``'s cluster.

    The main entry (opacity 1.0) is the selected pose if a pose was picked,
    otherwise the cluster's representative (highest MaSIF score); all other poses
    render at 0.5 opacity. Returns None if the cluster has no poses in df_all.
    """
    pidx = POSE_GROUPS.get(cluster_key_from_row(row))
    if not pidx:
        return None
    poses = df_all.loc[pidx]
    main_idx = None
    if selected_ref.get("src") == "pose" and selected_ref.get("idx") in poses.index:
        main_idx = selected_ref["idx"]
    else:
        reps = poses.index[poses["is_rep"]]
        main_idx = reps[0] if len(reps) else poses["score"].idxmax()
    others = poses.drop(main_idx)
    return [(df_all.at[main_idx, "flattened_transform"], 1.0)] + [(t, 0.5) for t in others["flattened_transform"]]


def build_pdb_view(pdb_id, tag):
    """Return an RCSB structure-page iframe block for a PDB id, or a note if invalid."""
    if not isinstance(pdb_id, str) or not pdb_id.strip():
        return html.Div(f"No valid {tag} PDB id for selected row.")
    pid = pdb_id.strip()
    url = f"https://www.rcsb.org/structure/{pid}"
    return html.Div(
        [
            html.Div(
                [f"{tag} ({pid}) — open in new tab: ",
                 html.A(url, href=url, target="_blank", rel="noopener noreferrer")],
                style={"marginBottom": "6px", "fontSize": "12px"},
            ),
            html.Iframe(src=url, style={"width": "100%", "height": "760px", "border": "none"}),
        ]
    )


def empty_fig(message):
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font={"size": 14, "color": "#666"})
    fig.update_layout(
        xaxis={"visible": False}, yaxis={"visible": False},
        height=560, margin={"l": 20, "r": 20, "t": 30, "b": 20},
    )
    return fig


# --------------------------------------------------------------------------- #
# Numeric range filters (restrict plotted rows; e.g. to get a big node under the
# plot cap). Adapted from the legacy app; no JSON/CSV persistence.
# --------------------------------------------------------------------------- #
def apply_filtering_criteria(frame, filtering_criteria):
    """Apply {column: (min_val, max_val)} range criteria to a frame."""
    filtered = frame
    for col, (min_val, max_val) in filtering_criteria.items():
        if col not in filtered.columns:
            continue
        if min_val is not None:
            filtered = filtered[filtered[col] >= min_val]
        if max_val is not None:
            filtered = filtered[filtered[col] <= max_val]
    return filtered


def build_filtering_criteria(filter_rows):
    """Normalise filter-row store state into a {column: (min, max)} dict."""
    criteria = {}
    for row in filter_rows or []:
        col = row.get("column")
        if not col:
            continue
        criteria[col] = (row.get("min"), row.get("max"))
    return criteria


def parse_cap(value):
    """Parse the plot-row-cap UI value, falling back to the default."""
    try:
        v = int(value)
        return v if v >= 1 else PLOT_ROW_CAP
    except (TypeError, ValueError):
        return PLOT_ROW_CAP


def min_cluster_size_for_cap(scope, cap):
    """Lowest integer cluster_size threshold whose ``>=`` filter fits ``scope`` under ``cap``.

    Returns None when the scope is already within the cap (no filter needed).
    """
    sizes = scope["cluster_size"].dropna()
    n = len(sizes)
    if n <= cap:
        return None
    sizes_int = sizes.astype(int)
    max_k = int(sizes_int.max())
    for k in range(1, max_k + 2):  # counts are monotonic-decreasing in k
        if int((sizes_int >= k).sum()) <= cap:
            return k
    return max_k + 1


def _filters_equiv(a, b):
    """Compare two filter-row lists ignoring their ids."""
    def norm(lst):
        return [(r.get("column"), r.get("min"), r.get("max")) for r in (lst or [])]
    return norm(a) == norm(b)


def auto_cluster_size_filter(top_path, filter_store, cap):
    """Return a filter list that auto-applies a cluster_size floor for a top-level scope.

    Existing non-cluster_size filters are preserved; the cluster_size filter is
    recomputed (or removed) so the plotted count fits under ``cap``.
    """
    scope = subframe_by_path(df, top_path)
    others = [f for f in (filter_store or []) if f.get("column") != "cluster_size"]
    scope_other = apply_filtering_criteria(scope, build_filtering_criteria(others))
    k = min_cluster_size_for_cap(scope_other, cap)
    if k is None:
        return others
    return others + [{"id": str(uuid.uuid4()), "column": "cluster_size", "min": k, "max": None}]


def build_filter_row(index, row_state=None):
    """Create a single column/min/max/remove filter control block."""
    row_state = row_state or {}
    return html.Div(
        [
            dcc.Dropdown(
                id={"type": "filter-column", "index": index},
                options=[{"label": c, "value": c} for c in NUMERIC_COLS],
                value=row_state.get("column"), placeholder="Column", clearable=True,
                style={"minWidth": "200px", "fontSize": "11px"},
            ),
            dcc.Input(id={"type": "filter-min", "index": index}, type="number",
                      value=row_state.get("min"), placeholder="Min", style={"width": "90px"}),
            dcc.Input(id={"type": "filter-max", "index": index}, type="number",
                      value=row_state.get("max"), placeholder="Max", style={"width": "90px"}),
            html.Button("✕", id={"type": "filter-remove", "index": index}, n_clicks=0,
                        title="Remove filter", style={"fontSize": "11px"}),
        ],
        style={"display": "flex", "gap": "6px", "alignItems": "center", "marginBottom": "6px"},
    )


def _grey_trace(frame_sub, x, y, src, name, size):
    """Grey-cross trace for non-selected dedup rows or overlaid poses."""
    if frame_sub.empty:
        return None
    custom = [[src, int(i)] for i in frame_sub.index]
    return go.Scatter(
        x=frame_sub[x], y=frame_sub[y], mode="markers",
        marker={"symbol": "x", "size": size, "color": "lightgrey", "line": {"width": 1, "color": "grey"}},
        customdata=custom, name=name,
        hovertemplate="%{x}, %{y}<extra>" + name + "</extra>",
    )


# --------------------------------------------------------------------------- #
# App layout
# --------------------------------------------------------------------------- #
app = Dash(__name__)
app.title = "Neosurf-Neosurf Explorer"

_order_dropdown = lambda i, default: dcc.Dropdown(
    id=f"order-{i}",
    options=[{"label": LEVEL_LABEL[l], "value": l} for l in NONLEAF_LEVELS],
    value=default, clearable=False, style={"fontSize": "11px", "marginBottom": "4px"},
)

# Resizable (drag the right edge) + retractable sidebar.
SIDEBAR_STYLE = {
    "width": "380px", "minWidth": "220px", "maxWidth": "900px", "flex": "0 0 auto",
    "padding": "10px", "borderRight": "1px solid #ddd",
    "resize": "horizontal", "overflow": "auto", "boxSizing": "border-box", "maxHeight": "96vh",
}

sidebar = html.Div(
    [
        html.Div(
            html.Button("◀ Hide", id="sidebar-toggle", n_clicks=0, title="Hide sidebar",
                        style={"fontSize": "11px", "cursor": "pointer"}),
            style={"textAlign": "right", "marginBottom": "6px"},
        ),
        html.H4("Nesting order", style={"marginBottom": "6px"}),
        html.Div(
            [
                html.Button("Preset A: target → seed", id="preset-a", n_clicks=0,
                            style={"fontSize": "10px", "marginRight": "4px"}),
                html.Button("Preset B: by protein", id="preset-b", n_clicks=0, style={"fontSize": "10px"}),
            ],
            style={"marginBottom": "6px"},
        ),
        _order_dropdown(1, PRESET_A[0]),
        _order_dropdown(2, PRESET_A[1]),
        _order_dropdown(3, PRESET_A[2]),
        _order_dropdown(4, PRESET_A[3]),
        html.Div("→ Cluster → Pose (fixed leaves)", style={"fontSize": "10px", "color": "#999", "marginBottom": "6px"}),
        html.Hr(),
        dcc.Input(id="tree-search", type="text", placeholder="Search gene / uniprot / id...",
                  debounce=True, style={"width": "100%", "marginBottom": "6px"}),
        html.Button("Collapse all", id="collapse-all", n_clicks=0, style={"fontSize": "10px", "marginBottom": "6px"}),
        dcc.Checklist(
            id="show-branch-scope-toggle",
            options=[{"label": " Show all entries in this branch", "value": "on"}],
            value=[],
            style={"fontSize": "11px", "marginBottom": "6px"},
        ),
        html.Div(
            id="tree-container",
            style={"height": "70vh", "overflowY": "auto", "overflowX": "auto",
                   "border": "1px solid #ddd", "borderRadius": "6px", "padding": "6px"},
        ),
    ],
    id="sidebar",
    style=SIDEBAR_STYLE,
)

# Thin strip shown when the sidebar is hidden; click to bring it back.
expand_strip = html.Div(
    html.Button("▶", id="sidebar-expand", n_clicks=0, title="Show sidebar",
                style={"fontSize": "13px", "cursor": "pointer", "writingMode": "vertical-rl"}),
    id="sidebar-expand-strip",
    style={"display": "none", "padding": "6px", "borderRight": "1px solid #ddd"},
)

axis_dd = lambda id_, default: dcc.Dropdown(
    id=id_, options=[{"label": c, "value": c} for c in NUMERIC_COLS], value=default, clearable=False,
    style={"fontSize": "12px"},
)

main = html.Div(
    [
        html.H2("MaSIF Neosurf-on-Neosurf Results Explorer", style={"marginTop": "0"}),
        html.Div(id="breadcrumb", style={"fontFamily": "monospace", "fontSize": "12px", "color": "#555",
                                         "marginBottom": "8px", "minHeight": "18px"}),
        html.Div(
            [
                html.Div([html.Label("X axis", style={"fontSize": "11px"}), axis_dd("x-dd", DEFAULT_X)],
                         style={"flex": "1", "marginRight": "8px"}),
                html.Div([html.Label("Y axis", style={"fontSize": "11px"}), axis_dd("y-dd", DEFAULT_Y)],
                         style={"flex": "1", "marginRight": "8px"}),
                html.Div([html.Label("Colour", style={"fontSize": "11px"}), axis_dd("color-dd", DEFAULT_COLOR)],
                         style={"flex": "1"}),
            ],
            style={"display": "flex", "marginBottom": "6px"},
        ),
        html.Details(
            [
                html.Summary("Filters — restrict plotted rows (e.g. to bring a large branch under the cap)",
                             style={"fontSize": "12px", "cursor": "pointer", "fontWeight": "bold"}),
                html.Div(
                    [
                        html.Button("Add filter", id="add-filter-btn", n_clicks=0, style={"fontSize": "11px"}),
                        html.Div(id="filter-rows-container", style={"marginTop": "8px"}),
                    ],
                    style={"marginTop": "6px"},
                ),
            ],
            style={"marginBottom": "8px", "padding": "6px", "border": "1px solid #eee", "borderRadius": "6px"},
        ),
        html.Div(
            [
                html.Label("Plot row cap: ", style={"fontSize": "12px", "marginRight": "6px"}),
                dcc.Input(id="plot-cap-input", type="number", value=PLOT_ROW_CAP, min=1, step=100,
                          debounce=True, style={"width": "110px"}),
                html.Span("  (top-level clicks auto-apply a cluster_size floor to fit this cap)",
                          style={"fontSize": "11px", "color": "#999", "marginLeft": "8px"}),
            ],
            style={"marginBottom": "6px"},
        ),
        html.Div(id="plot-info", style={"fontSize": "12px", "color": "#666", "marginBottom": "4px"}),
        # Plot (left half) + 3D viewer (right half) in the same row.
        html.Div(
            [
                html.Div(
                    dcc.Graph(id="scatter", clear_on_unhover=True, style={"height": "560px"},
                              config={"responsive": True}),
                    style={"flex": "1", "minWidth": "0", "padding": "10px", "border": "1px solid #ddd",
                           "borderRadius": "6px", "marginRight": "10px"},
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H4("3D viewer — target + seed (both ligands as sticks)",
                                        style={"display": "inline-block", "marginRight": "12px"}),
                                dcc.Checklist(
                                    id="show-cluster-toggle",
                                    options=[{"label": " Show cluster", "value": "on"}],
                                    value=[], inline=True,
                                    style={"display": "inline-block", "fontSize": "12px"},
                                ),
                            ]
                        ),
                        html.Div("Target: grey cartoon / green ligand. Seed: blue cartoon / magenta ligand. "
                                 "\"Show cluster\" overlays the cluster's other poses at 50% transparency.",
                                 style={"fontSize": "11px", "color": "#777", "marginBottom": "4px"}),
                        html.Div(id="viewer-container", children="Select a cluster / pose, or click a plot point."),
                    ],
                    style={"flex": "1", "padding": "10px", "border": "1px solid #ddd",
                           "borderRadius": "6px", "minWidth": "0"},
                ),
            ],
            style={"display": "flex", "marginTop": "6px"},
        ),
        # Full-width RCSB PDB structure page.
        html.Div(
            [
                html.Div(
                    [
                        html.H4("PDB structure (RCSB)", style={"display": "inline-block", "marginRight": "12px"}),
                        dcc.RadioItems(
                            id="pdb-side",
                            options=[{"label": "Matched (seed)", "value": "matched"},
                                     {"label": "Target", "value": "target"}],
                            value="matched", inline=True, style={"display": "inline-block", "fontSize": "12px"},
                        ),
                    ]
                ),
                html.Div(id="pdb-container", children="No selection."),
            ],
            style={"padding": "10px", "border": "1px solid #ddd", "borderRadius": "6px", "marginTop": "12px"},
        ),
        html.Div(
            [
                html.H4("Flagged clusters", style={"marginTop": "0", "marginBottom": "8px"}),
                html.Div(
                    [
                        html.Button("Flag", id="flag-toggle-btn", n_clicks=0,
                                    style={"fontSize": "12px", "marginRight": "8px"}),
                        html.Span(id="flag-status", children="Flagged: 0",
                                  style={"fontSize": "12px", "color": "#555"}),
                    ],
                    style={"marginBottom": "8px"},
                ),
                html.Div(
                    [
                        dcc.Input(
                            id="flagged-path-input", type="text", value=DEFAULT_FLAGGED_CSV,
                            debounce=True, style={"flex": "1", "fontSize": "12px", "minWidth": "0"},
                        ),
                        html.Button("Save flagged", id="save-flagged-btn", n_clicks=0,
                                    style={"fontSize": "12px", "marginLeft": "8px"}),
                    ],
                    style={"display": "flex", "alignItems": "center", "marginBottom": "4px"},
                ),
                html.Div(id="flag-save-status", children="",
                         style={"fontSize": "11px", "color": "#666", "minHeight": "16px"}),
            ],
            style={"padding": "10px", "border": "1px solid #ddd", "borderRadius": "6px", "marginTop": "12px"},
        ),
        html.Div([html.H4("Selected row"), html.Div(id="row-detail", children="No point selected.")],
                 style={"padding": "10px", "border": "1px solid #ddd", "borderRadius": "6px", "marginTop": "12px"}),
    ],
    style={"flex": "1", "padding": "12px", "minWidth": "0"},
)

app.layout = html.Div(
    [
        dcc.Store(id="tree-order-store", data=DEFAULT_ORDER),
        dcc.Store(id="filter-store", data=[]),
        dcc.Store(id="expanded-store", data=[]),
        dcc.Store(id="selected-branch-store", data=None),
        dcc.Store(id="selected-row-store", data=None),
        dcc.Store(id="flagged-store", data=[]),
        html.Div([expand_strip, sidebar, main], style={"display": "flex", "alignItems": "flex-start"}),
    ],
)


# --------------------------------------------------------------------------- #
# Callbacks
# --------------------------------------------------------------------------- #
@app.callback(
    Output("order-1", "value"), Output("order-2", "value"),
    Output("order-3", "value"), Output("order-4", "value"),
    Input("preset-a", "n_clicks"), Input("preset-b", "n_clicks"),
    prevent_initial_call=True,
)
def apply_preset(a_clicks, b_clicks):
    preset = PRESET_B if callback_context.triggered_id == "preset-b" else PRESET_A
    return preset[0], preset[1], preset[2], preset[3]


@app.callback(
    Output("sidebar", "style"),
    Output("sidebar-expand-strip", "style"),
    Input("sidebar-toggle", "n_clicks"),
    Input("sidebar-expand", "n_clicks"),
    prevent_initial_call=True,
)
def toggle_sidebar(_hide_clicks, _show_clicks):
    """Retract / restore the sidebar (width is drag-resizable via CSS resize)."""
    if callback_context.triggered_id == "sidebar-toggle":
        return {**SIDEBAR_STYLE, "display": "none"}, {"display": "block", "padding": "6px", "borderRight": "1px solid #ddd"}
    return SIDEBAR_STYLE, {"display": "none"}


@app.callback(
    Output("filter-store", "data"),
    Input("add-filter-btn", "n_clicks"),
    Input({"type": "filter-remove", "index": ALL}, "n_clicks"),
    Input({"type": "filter-column", "index": ALL}, "value"),
    Input({"type": "filter-min", "index": ALL}, "value"),
    Input({"type": "filter-max", "index": ALL}, "value"),
    State("filter-store", "data"),
    prevent_initial_call=True,
)
def update_filter_store(_add_clicks, _remove_clicks, columns, mins, maxs, filter_store):
    """Maintain the normalised filter-row state (add / remove / edit)."""
    store = copy.deepcopy(filter_store or [])
    ctx = callback_context
    if not ctx.triggered:
        return store
    trigger = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger == "add-filter-btn":
        store.append({"id": str(uuid.uuid4()), "column": None, "min": None, "max": None})
        return store

    if trigger.startswith("{"):
        trigger_id = json.loads(trigger)
        if trigger_id.get("type") == "filter-remove":
            target_idx = trigger_id.get("index")
            return [row for i, row in enumerate(store) if i != target_idx]

    for i, row in enumerate(store):
        row["column"] = columns[i] if i < len(columns) else row.get("column")
        row["min"] = mins[i] if i < len(mins) else row.get("min")
        row["max"] = maxs[i] if i < len(maxs) else row.get("max")
    return store


@app.callback(
    Output("filter-rows-container", "children"),
    Input("filter-store", "data"),
)
def render_filter_rows(filter_store):
    rows = filter_store or []
    if not rows:
        return html.Div("No filters. Click 'Add filter' to restrict which rows are plotted.",
                        style={"fontSize": "11px", "color": "#999"})
    return [build_filter_row(i, row_state=row) for i, row in enumerate(rows)]


@app.callback(
    Output("tree-order-store", "data"),
    Output("expanded-store", "data"),
    Output("selected-branch-store", "data"),
    Input("order-1", "value"), Input("order-2", "value"),
    Input("order-3", "value"), Input("order-4", "value"),
)
def update_order(o1, o2, o3, o4):
    # Changing nesting invalidates expanded paths and the selected branch.
    return build_order([o1, o2, o3, o4]), [], None


@app.callback(
    Output("expanded-store", "data", allow_duplicate=True),
    Input("collapse-all", "n_clicks"),
    prevent_initial_call=True,
)
def collapse_all(_n):
    return []


@app.callback(
    Output("tree-container", "children"),
    Input("tree-order-store", "data"),
    Input("expanded-store", "data"),
    Input("selected-branch-store", "data"),
    Input("tree-search", "value"),
    Input("filter-store", "data"),
)
def render_tree(order, expanded, selected, search, filter_store):
    base = df
    if search and search.strip():
        base = base[base["search_blob"].str.contains(search.strip().lower(), regex=False)]
    criteria = build_filtering_criteria(filter_store)
    if criteria:
        base = apply_filtering_criteria(base, criteria)
    if base.empty:
        return html.Div("No entries match the current search / filters.",
                        style={"fontSize": "12px", "color": "#999"})
    nodes = build_tree_nodes([], order, expanded or [], selected, base, criteria=criteria)
    return nodes or html.Div("No data.", style={"fontSize": "12px", "color": "#999"})


@app.callback(
    Output("expanded-store", "data", allow_duplicate=True),
    Output("selected-branch-store", "data", allow_duplicate=True),
    Output("selected-row-store", "data", allow_duplicate=True),
    Output("filter-store", "data", allow_duplicate=True),
    Input({"type": "tree-node", "path": ALL}, "n_clicks"),
    State("expanded-store", "data"),
    State("tree-order-store", "data"),
    State("filter-store", "data"),
    State("plot-cap-input", "value"),
    prevent_initial_call=True,
)
def on_tree_click(_all_clicks, expanded, order, filter_store, cap_value):
    ctx = callback_context
    if not ctx.triggered or not ctx.triggered[0]["value"]:
        raise PreventUpdate  # fired due to node set changing (re-render), not a click
    trig = ctx.triggered_id
    if not isinstance(trig, dict):
        raise PreventUpdate

    path_str = trig["path"]
    path = json.loads(path_str)
    expanded = list(expanded or [])
    level = path[-1][0]

    if level == POSE_LEVEL:
        # Leaf: select this exact pose.
        return expanded, path_str, {"src": "pose", "idx": int(path[-1][1])}, no_update

    # Any df level -> toggle expansion + set as selected branch.
    if path_str in expanded:
        expanded.remove(path_str)
    else:
        expanded.append(path_str)

    # A top-level (target protein) click changes the plotted set: auto-apply a
    # cluster_size floor if the scope would exceed the plot cap.
    filt_out = no_update
    if len(path) == 1:
        new_filters = auto_cluster_size_filter(path, filter_store, parse_cap(cap_value))
        if not _filters_equiv(new_filters, filter_store):
            filt_out = new_filters

    if level == CLUSTER_LEVEL:
        # Auto-select the cluster's representative (dedup) row for the viewer.
        sub = subframe_by_path(df, path)
        row_out = {"src": "dedup", "idx": int(sub.index[0])} if len(sub) else no_update
        return expanded, path_str, row_out, filt_out
    return expanded, path_str, no_update, filt_out


@app.callback(
    Output("scatter", "figure"),
    Output("plot-info", "children"),
    Input("selected-branch-store", "data"),
    Input("x-dd", "value"), Input("y-dd", "value"), Input("color-dd", "value"),
    Input("selected-row-store", "data"),
    Input("filter-store", "data"),
    Input("plot-cap-input", "value"),
    Input("flagged-store", "data"),
    Input("show-branch-scope-toggle", "value"),
)
def update_plot(branch_str, x, y, color, selected_ref, filter_store, cap_value, flagged_idxs, show_branch_scope):
    if not branch_str:
        return empty_fig("Select a node in the tree to plot its predictions."), ""

    cap = parse_cap(cap_value)
    path_df = strip_pose(json.loads(branch_str))
    show_all_in_branch = "on" in (show_branch_scope or [])
    criteria = build_filtering_criteria(filter_store)

    selected_rows = apply_filtering_criteria(subframe_by_path(df, path_df), criteria)
    if show_all_in_branch:
        # Plotted set is fixed by the TOP level only (target protein). Deeper selection
        # re-styles points: selected sub-branch stays coloured, the rest become grey crosses.
        scope_path = path_df[:1]
        scope = apply_filtering_criteria(subframe_by_path(df, scope_path), criteria)
    else:
        scope = selected_rows

    n = len(scope)
    n_unfiltered = len(subframe_by_path(df, path_df[:1] if show_all_in_branch else path_df))
    filt_note = f" (filtered from {n_unfiltered:,})" if criteria and n != n_unfiltered else ""
    if n == 0:
        return empty_fig("No rows in this branch after filtering."), f"0 rows{filt_note}"
    if n > cap:
        return (
            empty_fig(
                f"{n:,} rows in this branch — too many to plot.\n"
                f"Drill deeper, or add filters above to bring it under {cap:,}."
            ),
            f"{n:,} rows{filt_note} (exceeds cap of {cap:,})",
        )

    sel_plot = selected_rows.copy()
    sel_plot["_src"] = "dedup"
    sel_plot["_idx"] = sel_plot.index
    hover = [c for c in ["matched_gene_name", "matched_recommendedName", "cluster_id", "cluster_size"] if c in sel_plot.columns]
    fig = px.scatter(
        sel_plot, x=x, y=y, color=color, custom_data=["_src", "_idx"],
        hover_data=hover or None, color_continuous_scale="viridis", height=560,
    )
    fig.update_traces(marker={"size": 9})
    fig.update_layout(margin={"l": 40, "r": 20, "t": 30, "b": 40}, legend={"orientation": "h", "y": -0.18})

    if show_all_in_branch:
        non_selected = scope.drop(selected_rows.index, errors="ignore")
        grey = _grey_trace(non_selected, x, y, "dedup", "Other clusters in scope", 5)
        if grey is not None:
            fig.add_trace(grey)

    # If a specific cluster is pinned, overlay its individual poses.
    info_extra = ""
    if len(path_df) >= 5:  # all four levels + cluster present
        key = cluster_key_from_path(path_df)
        pose_idxs = POSE_GROUPS.get(key)
        if pose_idxs:
            poses = apply_filtering_criteria(df_all.loc[pose_idxs], criteria)
            pose_trace = _grey_trace(poses, x, y, "pose", "Poses in selected cluster", 5)
            if pose_trace is not None:
                pose_trace.marker.update({"color": "dimgrey", "line": {"width": 1, "color": "black"}})
                fig.add_trace(pose_trace)
            info_extra = f" | {len(poses)} poses overlaid"

    # Gold stars for flagged dedup rows currently in scope.
    flagged_set = set(flagged_idxs or [])
    if flagged_set:
        flagged_in_scope = scope.loc[scope.index.intersection(flagged_set)]
        if not flagged_in_scope.empty and x in flagged_in_scope.columns and y in flagged_in_scope.columns:
            fig.add_trace(go.Scatter(
                x=flagged_in_scope[x], y=flagged_in_scope[y], mode="markers",
                marker={
                    "symbol": "star", "size": 14, "color": "gold",
                    "line": {"width": 1, "color": "darkorange"},
                },
                name="Flagged", hoverinfo="skip", showlegend=True,
            ))
            info_extra += f" | {len(flagged_in_scope)} flagged"

    # Red X for the specific selected row (dedup rep or a pose).
    sel_row = resolve_ref(selected_ref)
    if sel_row is not None and x in sel_row and y in sel_row and pd.notna(sel_row[x]) and pd.notna(sel_row[y]):
        fig.add_trace(go.Scatter(
            x=[sel_row[x]], y=[sel_row[y]], mode="markers",
            marker={"symbol": "x", "size": 18, "color": "red", "line": {"width": 3, "color": "red"}},
            name="Selected", hoverinfo="skip", showlegend=False,
        ))

    return fig, (
        f"{len(selected_rows):,} selected / {n:,} in scope{filt_note}{info_extra}"
        if show_all_in_branch
        else f"{len(selected_rows):,} plotted{filt_note}{info_extra}"
    )


@app.callback(
    Output("flagged-store", "data"),
    Input("flag-toggle-btn", "n_clicks"),
    State("selected-row-store", "data"),
    State("flagged-store", "data"),
    prevent_initial_call=True,
)
def toggle_flag(_n_clicks, selected_ref, flagged_idxs):
    idx = dedup_idx_from_ref(selected_ref)
    if idx is None:
        raise PreventUpdate
    flagged = set(flagged_idxs or [])
    if idx in flagged:
        flagged.remove(idx)
    else:
        flagged.add(idx)
    return sorted(flagged)


@app.callback(
    Output("flag-toggle-btn", "children"),
    Output("flag-status", "children"),
    Input("selected-row-store", "data"),
    Input("flagged-store", "data"),
)
def update_flag_controls(selected_ref, flagged_idxs):
    flagged = flagged_idxs or []
    n = len(flagged)
    idx = dedup_idx_from_ref(selected_ref)
    if idx is None:
        return "Flag", f"Flagged: {n} — select a cluster or pose first"
    if idx in flagged:
        return "Unflag", f"Flagged: {n} — current cluster is flagged"
    return "Flag", f"Flagged: {n} — current cluster not flagged"


@app.callback(
    Output("flag-save-status", "children"),
    Input("save-flagged-btn", "n_clicks"),
    State("flagged-path-input", "value"),
    State("flagged-store", "data"),
    prevent_initial_call=True,
)
def save_flagged(_n_clicks, path, flagged_idxs):
    flagged = flagged_idxs or []
    if not flagged:
        return "Nothing to save."
    path = (path or "").strip()
    if not path:
        return "Enter a CSV path first."
    try:
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        out = flagged_frame(flagged)
        out.to_csv(path, index=False)
        return f"Saved {len(out)} rows to {path}"
    except OSError as exc:
        return f"Save failed: {exc}"


@app.callback(
    Output("selected-row-store", "data", allow_duplicate=True),
    Input("scatter", "clickData"),
    prevent_initial_call=True,
)
def on_point_click(click_data):
    if not click_data or not click_data.get("points"):
        raise PreventUpdate
    custom = click_data["points"][0].get("customdata")
    if not isinstance(custom, (list, tuple)) or len(custom) < 2:
        raise PreventUpdate
    src, idx = custom[0], custom[1]
    try:
        return {"src": src, "idx": int(idx)}
    except (TypeError, ValueError):
        raise PreventUpdate


@app.callback(
    Output("breadcrumb", "children"),
    Input("selected-branch-store", "data"),
)
def update_breadcrumb(branch_str):
    if not branch_str:
        return "No branch selected."
    path = json.loads(branch_str)
    labels = []
    for i in range(1, len(path) + 1):
        level, value = path[i - 1]
        if level == POSE_LEVEL:
            prow = resolve_ref({"src": "pose", "idx": int(value)})
            labels.append(pose_label(prow) if prow is not None else "pose")
            continue
        sub = subframe_by_path(df, path[:i])
        if sub.empty:
            break
        rep = df.loc[sub["score"].idxmax()]
        labels.append(node_label(level, value, rep))
    return " › ".join(labels)


@app.callback(
    Output("viewer-container", "children"),
    Output("row-detail", "children"),
    Output("pdb-container", "children"),
    Input("selected-row-store", "data"),
    Input("pdb-side", "value"),
    Input("show-cluster-toggle", "value"),
)
def update_selection(selected_ref, pdb_side, show_cluster):
    row = resolve_ref(selected_ref)
    if row is None:
        return "Select a cluster / pose, or click a plot point.", "No point selected.", "No selection."

    src = selected_ref.get("src")
    try:
        pose_transforms = None
        if show_cluster:
            pose_transforms = cluster_pose_transforms(row, selected_ref)
        if pose_transforms is None:
            pose_transforms = [(row.get("flattened_transform"), 1.0)]
        viewer_html = build_viewer_html(row, pose_transforms)
        viewer = html.Iframe(srcDoc=viewer_html, style={"width": "100%", "height": "560px", "border": "none"})
    except Exception as exc:  # noqa: BLE001
        viewer = html.Pre(f"Unable to render structure: {exc}", style={"whiteSpace": "pre-wrap"})

    detail_cols = [
        "target", "target_pdb_id", "t_uniprot_id", "target_gene_name", "t_recommendedName", "target_ligand_code",
        "matched_protein", "matched_pdb_id", "matched_uniprot_id", "matched_gene_name", "matched_recommendedName", "matched_ligand_code",
        "cluster_id", "cluster_size", "cluster_mean_rmsd", "target_vix", "matched_patch_id", "is_rep",
        "score", "desc_dist_score", "iface_score", "mean_desc_dist_score",
        "clashing_ca", "clashing_heavy", "n_match_ligands", "matched_mw",
    ]
    detail_cols = [c for c in detail_cols if c in row.index]
    header = f"src={src} | idx={selected_ref.get('idx')}"
    rows_html = [html.Tr([html.Th("row"), html.Td(header)])]
    rows_html += [html.Tr([html.Th(c), html.Td(str(row[c]))]) for c in detail_cols]
    detail = html.Table(rows_html, style={"width": "100%", "borderCollapse": "collapse", "fontSize": "12px"})

    if pdb_side == "target":
        pdb = build_pdb_view(row.get("target_pdb_id"), "Target")
    else:
        pdb = build_pdb_view(row.get("matched_pdb_id"), "Matched (seed)")

    return viewer, detail, pdb


if __name__ == "__main__":
    app.run(debug=True)
