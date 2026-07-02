"""
Neosurf-on-Neosurf Results Explorer (Dash)
==========================================

Interactive explorer for the all-vs-all MaSIF-Neosurf search produced by
``notebooks/neosurf_neosurf.ipynb``.

Unlike the single-target legacy app (``results_explorer_app.py``), this workflow
queries *many* ligand-bound "target" complexes against a large library of
ligand-bound "seed" complexes, so every predicted interaction row has a
different target AND a different seed, each carrying its own ligand.

Two tables drive the app:
  * ``df_preprocess_ok.csv``  -> per-structure metadata (uniprot/gene/name,
    pdb_path, ligand identity), keyed by the unique ``target`` id. This is the
    source of metadata for BOTH targets and seeds.
  * ``df_results_dedup.csv``  -> one row per (target, matched_protein, cluster_id)
    predicted interaction, with metrics and the superposition transform.

Panels
  * Left sidebar: a lazily-rendered, re-orderable collapsible tree that lets you
    drill target-protein -> target-id -> matched-protein -> matched-seed ->
    cluster. Selecting a node plots only that branch (so Plotly never chokes).
  * Plotly scatter of any two metrics for the selected branch.
  * py3Dmol viewer showing target + transform-superposed seed, with BOTH
    ligands rendered as sticks.
  * UniProt entry iframe for the matched (or target) protein.

The legacy app is intentionally left untouched.
"""

import json

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
RESULTS_CSV = "/scratch/ymeng/Neosurf_Neosurf/data/processing/3_masif_search/df_results_dedup.csv"

# Max rows a selected branch may contain before we refuse to plot (guardrail so
# Plotly is never handed tens of thousands of interactive points).
PLOT_ROW_CAP = 4000

# --------------------------------------------------------------------------- #
# Grouping levels for the tree
# --------------------------------------------------------------------------- #
# Non-leaf levels can be re-ordered by the user; ``cluster_id`` is always the
# leaf (a single df_results_dedup row).
LEAF_LEVEL = "cluster_id"
NONLEAF_LEVELS = ["target_protein", "target_id", "matched_protein", "matched_seed_id"]

LEVEL_LABEL = {
    "target_protein": "Target protein",
    "target_id": "Target complex (id)",
    "matched_protein": "Matched protein",
    "matched_seed_id": "Matched seed (id)",
    "cluster_id": "Cluster",
}

# Grouping column backing each level (string columns built at load time).
GROUP_COL = {
    "target_protein": "g_target_protein",
    "target_id": "g_target_id",
    "matched_protein": "g_matched_protein",
    "matched_seed_id": "g_matched_seed_id",
    "cluster_id": "g_cluster_id",
}

DEFAULT_ORDER = ["target_protein", "target_id", "matched_protein", "matched_seed_id", LEAF_LEVEL]
PRESET_A = ["target_protein", "target_id", "matched_protein", "matched_seed_id"]
PRESET_B = ["target_protein", "matched_protein", "target_id", "matched_seed_id"]

# Metric columns offered as plot axes / colours.
CANDIDATE_NUMERIC = [
    "score", "desc_dist_score", "iface_score", "mean_desc_dist_score", "desc_dist",
    "cluster_size", "cluster_mean_rmsd", "clashing_ca", "clashing_heavy",
    "n_match_ligands", "matched_mw",
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


def load_data():
    """Load both tables and build the enriched results frame used everywhere."""
    pre = pd.read_csv(PREPROCESS_CSV)
    res = pd.read_csv(RESULTS_CSV).reset_index(drop=True)

    # Manifest lookup keyed by the unique `target` id (never parse id strings).
    pre_idx = pre.set_index("target")

    def mp(col):
        return pre_idx[col] if col in pre_idx.columns else pd.Series(dtype=object)

    # Target-side metadata (results only carries target_gene_name).
    res["t_uniprot_id"] = res["target"].map(mp("uniprot_id"))
    res["t_recommendedName"] = res["target"].map(mp("recommendedName"))
    res["target_pdb_path"] = res["target"].map(mp("pdb_path"))
    res["target_ligand_resname"] = res["target"].map(mp("ligand_resname"))
    res["target_pdb_ligand_chain"] = res["target"].map(mp("pdb_ligand_chain"))
    res["target_ligand_code"] = res["target"].map(mp("ligand_code"))

    # Matched/seed-side structural metadata (uniprot/gene/name already present).
    res["matched_pdb_path"] = res["matched_protein"].map(mp("pdb_path"))
    res["matched_ligand_resname"] = res["matched_protein"].map(mp("ligand_resname"))
    res["matched_pdb_ligand_chain"] = res["matched_protein"].map(mp("pdb_ligand_chain"))
    res["matched_ligand_code"] = res["matched_protein"].map(mp("ligand_code"))

    # String grouping columns (NaN-safe, JSON-safe path values).
    res["g_target_protein"] = res["t_uniprot_id"].fillna("—").astype(str)
    res["g_target_id"] = res["target"].astype(str)
    res["g_matched_protein"] = res["matched_uniprot_id"].fillna("—").astype(str)
    res["g_matched_seed_id"] = res["matched_protein"].astype(str)
    res["g_cluster_id"] = res["cluster_id"].map(
        lambda v: "" if pd.isna(v) else (str(int(v)) if float(v).is_integer() else str(v))
    )

    # Free-text search blob across all human-facing identifiers.
    search_cols = [
        "target_gene_name", "t_uniprot_id", "t_recommendedName",
        "matched_gene_name", "matched_uniprot_id", "matched_recommendedName",
        "target", "matched_protein",
    ]
    search_cols = [c for c in search_cols if c in res.columns]
    res["search_blob"] = res[search_cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower()

    return pre, res


pre_df, df = load_data()
NUMERIC_COLS = [c for c in CANDIDATE_NUMERIC if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]
DEFAULT_X = "score" if "score" in NUMERIC_COLS else NUMERIC_COLS[0]
DEFAULT_Y = "cluster_size" if "cluster_size" in NUMERIC_COLS else NUMERIC_COLS[min(1, len(NUMERIC_COLS) - 1)]
DEFAULT_COLOR = "n_match_ligands" if "n_match_ligands" in NUMERIC_COLS else NUMERIC_COLS[0]

_parser = PDBParser(QUIET=True)


# --------------------------------------------------------------------------- #
# Tree logic
# --------------------------------------------------------------------------- #
def subframe_by_path(base, path):
    """Filter ``base`` to rows matching every (level, value) pair in ``path``."""
    if not path:
        return base
    mask = pd.Series(True, index=base.index)
    for level, value in path:
        mask &= base[GROUP_COL[level]] == value
    return base[mask]


def node_label(level, value, rep):
    """Human label for a node, given a representative row ``rep``."""
    if level == "target_protein":
        return f"{_s(rep.get('t_uniprot_id'))}|{_s(rep.get('target_gene_name'))}|{_s(rep.get('t_recommendedName'))}"
    if level == "target_id":
        code = _s(rep.get("target_ligand_code"))
        return f"{value}" + (f"  [{code}]" if code else "")
    if level == "matched_protein":
        return (
            f"{_s(rep.get('matched_uniprot_id'))}|{_s(rep.get('matched_gene_name'))}|"
            f"{_s(rep.get('matched_recommendedName'))}"
        )
    if level == "matched_seed_id":
        code = _s(rep.get("matched_ligand_code"))
        return f"{value}" + (f"  [{code}]" if code else "")
    if level == "cluster_id":
        size = rep.get("cluster_size")
        score = rep.get("score")
        size_s = f"size={int(size)}" if pd.notna(size) else "size=?"
        score_s = f"score={score:.3f}" if pd.notna(score) else "score=?"
        return f"cluster {value} ({size_s}, {score_s})"
    return str(value)


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


def build_tree_nodes(path, order, expanded, selected, base, depth=0):
    """Recursively build tree nodes, descending only into expanded paths."""
    level_idx = len(path)
    if level_idx >= len(order):
        return []
    level = order[level_idx]
    col = GROUP_COL[level]
    sub = subframe_by_path(base, path)
    if sub.empty:
        return []

    # One entry per distinct value at this level, with count + best score.
    reps = []
    for value, idx in sub.groupby(col, sort=False).groups.items():
        ssub = sub.loc[idx]
        best_idx = ssub["score"].idxmax()
        reps.append((value, len(idx), df.loc[best_idx], float(ssub["score"].max())))
    reps.sort(key=lambda t: t[3], reverse=True)  # best-scoring siblings first

    is_leaf_level = level_idx == len(order) - 1
    nodes = []
    for value, count, rep, best in reps:
        child_path = path + [[level, value]]
        path_str = json.dumps(child_path)
        is_expanded = path_str in expanded
        marker = "•" if is_leaf_level else ("▾" if is_expanded else "▸")
        label = node_label(level, value, rep)
        badge = "" if is_leaf_level else f"  ({count})"
        button = html.Button(
            f"{marker} {label}{badge}",
            id={"type": "tree-node", "path": path_str},
            n_clicks=0,
            title=label,
            style=_node_style(depth, is_leaf_level, path_str == selected),
        )
        block = [button]
        if is_expanded and not is_leaf_level:
            block.append(
                html.Div(build_tree_nodes(child_path, order, expanded, selected, base, depth + 1))
            )
        nodes.append(html.Div(block))
    return nodes


def build_order(vals):
    """Normalise 4 dropdown values into a valid level permutation + leaf."""
    seen = []
    for v in vals:
        if v and v not in seen:
            seen.append(v)
    for v in NONLEAF_LEVELS:
        if v not in seen:
            seen.append(v)
    return seen + [LEAF_LEVEL]


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


def view_pair_html(row):
    """Return py3Dmol HTML: target + transform-superposed seed, both ligands as sticks."""
    import os

    target_path = row.get("target_pdb_path")
    seed_path = row.get("matched_pdb_path")
    transform = row.get("flattened_transform")
    for label, p in (("target", target_path), ("seed", seed_path)):
        if not isinstance(p, str) or not p.strip():
            raise ValueError(f"Missing {label} pdb_path for selected row.")
        if not os.path.exists(p):
            raise FileNotFoundError(f"{label} PDB not found: {p}")
    if not isinstance(transform, str) or not transform.strip():
        raise ValueError("Missing flattened_transform for selected row.")

    target_struct = _parser.get_structure("target", target_path)
    seed_struct = _parser.get_structure("seed", seed_path)
    apply_transform(seed_struct, transform)

    target_pdb = structure_to_pdbstr(target_struct)
    seed_pdb = structure_to_pdbstr(seed_struct)

    view = py3Dmol.view(width=760, height=560)
    view.addModel(target_pdb, "pdb")  # model 0 = target
    view.setStyle({"model": 0}, {"cartoon": {"color": "lightgrey"}})
    view.addModel(seed_pdb, "pdb")  # model 1 = seed
    view.setStyle({"model": 1}, {"cartoon": {"color": "lightblue"}})

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

    # Seed ligand -> magenta sticks.
    seed_lig = _ligand_selection(1, row.get("matched_ligand_resname"), row.get("matched_pdb_ligand_chain"))
    if "resn" in seed_lig:
        view.setStyle(seed_lig, {"stick": {"colorscheme": "magentaCarbon", "radius": 0.22}})
    else:
        for het in find_hetatm_residues(seed_struct):
            sel = {"model": 1, "resn": het["resn"], "resi": het["resi"]}
            if het["chain"]:
                sel["chain"] = het["chain"]
            view.setStyle(sel, {"stick": {"colorscheme": "magentaCarbon", "radius": 0.22}})

    # Focus on the ligand-overlap region if we could identify the target ligand.
    if "resn" in target_lig:
        view.zoomTo(target_lig)
    else:
        view.zoomTo()
    return view._make_html()


def build_uniprot_view(accession, tag):
    """Return a UniProt iframe block for an accession, or a note if invalid."""
    if not isinstance(accession, str) or not accession.strip() or accession.strip() == "—":
        return html.Div(f"No valid {tag} UniProt accession for selected row.")
    accn = accession.strip()
    url = f"https://www.uniprot.org/uniprotkb/{accn}/entry"
    return html.Div(
        [
            html.Div(
                [f"{tag}: ", html.A(url, href=url, target="_blank", rel="noopener noreferrer")],
                style={"marginBottom": "6px", "fontSize": "12px"},
            ),
            html.Iframe(src=url, style={"width": "100%", "height": "560px", "border": "none"}),
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
# App layout
# --------------------------------------------------------------------------- #
app = Dash(__name__)
app.title = "Neosurf-Neosurf Explorer"

_order_dropdown = lambda i, default: dcc.Dropdown(
    id=f"order-{i}",
    options=[{"label": LEVEL_LABEL[l], "value": l} for l in NONLEAF_LEVELS],
    value=default,
    clearable=False,
    style={"fontSize": "11px", "marginBottom": "4px"},
)

sidebar = html.Div(
    [
        html.H4("Nesting order", style={"marginBottom": "6px"}),
        html.Div(
            [
                html.Button("Preset A: target → seed", id="preset-a", n_clicks=0,
                            style={"fontSize": "10px", "marginRight": "4px"}),
                html.Button("Preset B: by protein", id="preset-b", n_clicks=0,
                            style={"fontSize": "10px"}),
            ],
            style={"marginBottom": "6px"},
        ),
        _order_dropdown(1, PRESET_A[0]),
        _order_dropdown(2, PRESET_A[1]),
        _order_dropdown(3, PRESET_A[2]),
        _order_dropdown(4, PRESET_A[3]),
        html.Hr(),
        dcc.Input(
            id="tree-search", type="text", placeholder="Search gene / uniprot / id...",
            debounce=True, style={"width": "100%", "marginBottom": "6px"},
        ),
        html.Button("Collapse all", id="collapse-all", n_clicks=0, style={"fontSize": "10px", "marginBottom": "6px"}),
        html.Div(
            id="tree-container",
            style={"height": "70vh", "overflowY": "auto", "overflowX": "auto",
                   "border": "1px solid #ddd", "borderRadius": "6px", "padding": "6px"},
        ),
    ],
    style={"width": "360px", "flex": "0 0 360px", "padding": "10px",
           "borderRight": "1px solid #ddd"},
)

axis_dd = lambda id_, default: dcc.Dropdown(
    id=id_, options=[{"label": c, "value": c} for c in NUMERIC_COLS], value=default, clearable=False,
    style={"fontSize": "12px"},
)

main = html.Div(
    [
        html.H2("MaSIF Neosurf-on-Neosurf Results Explorer", style={"marginTop": "0"}),
        html.Div(id="breadcrumb", style={"fontFamily": "monospace", "fontSize": "12px",
                                         "color": "#555", "marginBottom": "8px", "minHeight": "18px"}),
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
        html.Div(id="plot-info", style={"fontSize": "12px", "color": "#666", "marginBottom": "4px"}),
        dcc.Graph(id="scatter", clear_on_unhover=True, style={"height": "560px"},
                  config={"responsive": True}),
        html.Div(
            [
                html.Div(
                    [
                        html.H4("3D viewer — target + seed (both ligands as sticks)"),
                        html.Div(
                            "Target: grey cartoon / green ligand. Seed: blue cartoon / magenta ligand.",
                            style={"fontSize": "11px", "color": "#777", "marginBottom": "4px"},
                        ),
                        html.Div(id="viewer-container", children="Select a cluster (leaf) or plot point."),
                    ],
                    style={"flex": "1", "padding": "10px", "border": "1px solid #ddd",
                           "borderRadius": "6px", "marginRight": "10px", "minWidth": "0"},
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H4("UniProt", style={"display": "inline-block", "marginRight": "12px"}),
                                dcc.RadioItems(
                                    id="uniprot-side",
                                    options=[{"label": "Matched", "value": "matched"},
                                             {"label": "Target", "value": "target"}],
                                    value="matched", inline=True, style={"display": "inline-block", "fontSize": "12px"},
                                ),
                            ]
                        ),
                        html.Div(id="uniprot-container", children="No selection."),
                    ],
                    style={"flex": "1", "padding": "10px", "border": "1px solid #ddd",
                           "borderRadius": "6px", "minWidth": "0"},
                ),
            ],
            style={"display": "flex", "marginTop": "12px"},
        ),
        html.Div(
            [html.H4("Selected row"), html.Div(id="row-detail", children="No point selected.")],
            style={"padding": "10px", "border": "1px solid #ddd", "borderRadius": "6px", "marginTop": "12px"},
        ),
    ],
    style={"flex": "1", "padding": "12px", "minWidth": "0"},
)

app.layout = html.Div(
    [
        dcc.Store(id="tree-order-store", data=DEFAULT_ORDER),
        dcc.Store(id="expanded-store", data=[]),
        dcc.Store(id="selected-branch-store", data=None),
        dcc.Store(id="selected-row-store", data=None),
        html.Div([sidebar, main], style={"display": "flex", "alignItems": "flex-start"}),
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
    trig = callback_context.triggered_id
    preset = PRESET_B if trig == "preset-b" else PRESET_A
    return preset[0], preset[1], preset[2], preset[3]


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
)
def render_tree(order, expanded, selected, search):
    base = df
    if search and search.strip():
        base = df[df["search_blob"].str.contains(search.strip().lower(), regex=False)]
        if base.empty:
            return html.Div("No matches.", style={"fontSize": "12px", "color": "#999"})
    nodes = build_tree_nodes([], order, expanded or [], selected, base)
    return nodes or html.Div("No data.", style={"fontSize": "12px", "color": "#999"})


@app.callback(
    Output("expanded-store", "data", allow_duplicate=True),
    Output("selected-branch-store", "data", allow_duplicate=True),
    Output("selected-row-store", "data", allow_duplicate=True),
    Input({"type": "tree-node", "path": ALL}, "n_clicks"),
    State("expanded-store", "data"),
    State("tree-order-store", "data"),
    prevent_initial_call=True,
)
def on_tree_click(_all_clicks, expanded, order):
    ctx = callback_context
    if not ctx.triggered or not ctx.triggered[0]["value"]:
        # Fired because the node set changed (re-render), not a real click.
        raise PreventUpdate
    trig = ctx.triggered_id
    if not isinstance(trig, dict):
        raise PreventUpdate

    path_str = trig["path"]
    path = json.loads(path_str)
    expanded = list(expanded or [])
    is_leaf = len(path) >= len(order)

    row_out = no_update
    if is_leaf:
        sub = subframe_by_path(df, path)
        row_out = int(sub.index[0]) if len(sub) else no_update
    else:
        if path_str in expanded:
            expanded.remove(path_str)
        else:
            expanded.append(path_str)
    return expanded, path_str, row_out


@app.callback(
    Output("scatter", "figure"),
    Output("plot-info", "children"),
    Input("selected-branch-store", "data"),
    Input("x-dd", "value"), Input("y-dd", "value"), Input("color-dd", "value"),
    Input("selected-row-store", "data"),
)
def update_plot(branch_str, x, y, color, selected_row):
    if not branch_str:
        return empty_fig("Select a node in the tree to plot its predictions."), ""
    path = json.loads(branch_str)
    sub = subframe_by_path(df, path)
    n = len(sub)
    if n == 0:
        return empty_fig("No rows in this branch."), "0 rows"
    if n > PLOT_ROW_CAP:
        return (
            empty_fig(f"{n:,} rows in this branch — too many to plot.\nDrill deeper into the tree."),
            f"{n:,} rows (exceeds cap of {PLOT_ROW_CAP:,})",
        )

    sub = sub.copy()
    sub["_row_idx"] = sub.index
    hover = [c for c in ["matched_gene_name", "matched_recommendedName", "cluster_id", "cluster_size"] if c in sub.columns]
    fig = px.scatter(
        sub, x=x, y=y, color=color, custom_data=["_row_idx"],
        hover_data=hover or None, color_continuous_scale="viridis", height=560,
    )
    fig.update_layout(margin={"l": 40, "r": 20, "t": 30, "b": 40})

    if selected_row is not None and selected_row in sub.index:
        r = sub.loc[selected_row]
        fig.add_trace(
            go.Scatter(
                x=[r[x]], y=[r[y]], mode="markers",
                marker={"symbol": "x", "size": 16, "color": "red", "line": {"width": 2, "color": "red"}},
                name="Selected", hoverinfo="skip", showlegend=False,
            )
        )
    return fig, f"{n:,} predictions in selected branch"


@app.callback(
    Output("selected-row-store", "data", allow_duplicate=True),
    Input("scatter", "clickData"),
    prevent_initial_call=True,
)
def on_point_click(click_data):
    if not click_data or not click_data.get("points"):
        raise PreventUpdate
    custom = click_data["points"][0].get("customdata")
    idx = custom[0] if isinstance(custom, (list, tuple)) and custom else custom
    try:
        return int(idx)
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
        prefix = path[:i]
        sub = subframe_by_path(df, prefix)
        if sub.empty:
            break
        rep = df.loc[sub["score"].idxmax()]
        labels.append(node_label(prefix[-1][0], prefix[-1][1], rep))
    return " › ".join(labels)


@app.callback(
    Output("viewer-container", "children"),
    Output("row-detail", "children"),
    Output("uniprot-container", "children"),
    Input("selected-row-store", "data"),
    Input("uniprot-side", "value"),
)
def update_selection(row_idx, uniprot_side):
    if row_idx is None or row_idx not in df.index:
        return "Select a cluster (leaf) or plot point.", "No point selected.", "No selection."
    row = df.loc[row_idx]

    # 3D viewer
    try:
        viewer_html = view_pair_html(row)
        viewer = html.Iframe(srcDoc=viewer_html, style={"width": "100%", "height": "560px", "border": "none"})
    except Exception as exc:  # noqa: BLE001
        viewer = html.Pre(f"Unable to render structure: {exc}", style={"whiteSpace": "pre-wrap"})

    # Row detail table
    detail_cols = [
        "target", "t_uniprot_id", "target_gene_name", "t_recommendedName", "target_ligand_code",
        "matched_protein", "matched_uniprot_id", "matched_gene_name", "matched_recommendedName", "matched_ligand_code",
        "cluster_id", "cluster_size", "cluster_mean_rmsd",
        "score", "desc_dist_score", "iface_score", "mean_desc_dist_score",
        "clashing_ca", "clashing_heavy", "n_match_ligands", "matched_mw",
    ]
    detail_cols = [c for c in detail_cols if c in df.columns]
    rows_html = [html.Tr([html.Th("index"), html.Td(str(row_idx))])]
    rows_html += [html.Tr([html.Th(c), html.Td(str(row[c]))]) for c in detail_cols]
    detail = html.Table(rows_html, style={"width": "100%", "borderCollapse": "collapse", "fontSize": "12px"})

    # UniProt
    if uniprot_side == "target":
        uniprot = build_uniprot_view(row.get("t_uniprot_id"), "Target")
    else:
        uniprot = build_uniprot_view(row.get("matched_uniprot_id"), "Matched")

    return viewer, detail, uniprot


if __name__ == "__main__":
    app.run(debug=True)
