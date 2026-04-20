# MaSIF-Coherence: a learned protein-level interaction score on top of `masif_ppi_search` descriptors

## 1. What `masif_ppi_search` gives you per protein

Grounding everything in the actual codebase (paths relative to `masif/`):

The `MaSIF_ppi_search` network (`source/masif_modules/MaSIF_ppi_search.py`) takes the polar geodesic patch around every surface vertex — 5 input features (shape index, distance-dependent curvature, H-bond acceptor/donor potential, Poisson–Boltzmann charge, hydrophobicity) sampled on a polar grid of `n_rhos=5` × `n_thetas=16` — and produces an 80-dimensional descriptor per vertex. `masif_ppi_search_comp_desc.py` writes these to disk as `p1_desc_straight.npy` and `p1_desc_flipped.npy` (shape `(V, 80)`, float32); the "flipped" version negates the physico-chemical features and reflects θ → 2π − θ so that a complementary interface can be scored by a Euclidean distance between `desc_straight(A)` and `desc_flipped(B)`.

After preprocessing, for each protein you also have on disk:

- `p{1,2}_rho_wrt_center.npy`, `p{1,2}_theta_wrt_center.npy` — the polar patch coordinates, shape `(V, 200)`.
- `p{1,2}_input_feat.npy` — shape `(V, 200, 5)`, the 5-channel patch features.
- `p{1,2}_mask.npy` — valid-neighbor mask, shape `(V, 200)`.
- `p{1,2}_sc_labels.npy` — per-vertex median shape-complementarity labels (optional).
- The underlying `.ply` mesh in the precomputation directory gives the 3D vertex coordinates and vertex normals.

Typical numbers (`source/default_config/masif_opts.py`): patch radius 12 Å, `max_shape_size = 200` neighbors per patch, per-protein surfaces of order 10³ vertices.

## 2. How MaSIF currently turns descriptors into a protein-level decision

Two stages today (`source/masif_ppi_search/second_stage_alignment_nn.py`):

1. **Descriptor nearest-neighbor**: each target centre vertex is matched to candidate source vertices by raw Euclidean distance in the 80-D descriptor space (`desc_straight` vs `desc_flipped`).
2. **Geometric verification**: RANSAC + ICP estimates a rigid transform between the two patches, and a small `ScoreNN` (`transformation_training_data/score_nn.py`) reads three post-alignment per-point features — `1/spatial_distance`, `1/descriptor_distance`, `normal_dot_product` — over up to 200 points, passes them through 6 Conv1D layers, a global average pool, and a classifier, to output a binder/non-binder probability.

This is exactly what you don't want: the scoring network is conditioned on a rigid alignment, and the "geometric consistency" is a hard pose hypothesis. Your proposed model replaces stage 2 with something that reasons relationally about spatial structure, without ever fitting a transform.

## 3. Design goals in one line

Given two point clouds of surface vertices, each carrying a 80-D MaSIF descriptor and local geometric context, produce a single interaction score that is high when there exists a *subset of mutually descriptor-compatible vertex pairs whose intra-protein spatial relationships are mutually consistent*, and low otherwise — with this notion of "consistency" learned from data, invariant to independent rigid motions of each protein, invariant to vertex ordering, and tractable at thousands of vertices.

## 4. Proposed architecture: MaSIF-Coherence

The model has four stages: a per-protein geometric encoder, a soft cross-protein match layer, a relational match-graph transformer that is the core innovation, and a global readout. SE(3) invariance is achieved by never letting absolute coordinates enter — coordinates only appear via pairwise distances and angles between normals and edge vectors.

### 4.1 Inputs and subsampling

For each protein P you feed in `X_P ∈ R^{V×3}` (coordinates), `N_P ∈ R^{V×3}` (normals, normalised), `F_P ∈ R^{V×5}` (per-vertex physico-chemical features, taken from the centre row of `input_feat`), `D_P^s, D_P^f ∈ R^{V×80}` (straight and flipped MaSIF descriptors).

Because surfaces can have several thousand vertices, run farthest-point sampling on `X_P` to `M = 512` vertices per protein. FPS preserves surface coverage and is permutation-equivariant — the identity of sampled vertices is determined by geometry, and everything downstream is permutation-invariant anyway, so this introduces no order dependence.

Following the MaSIF complementarity convention, feed protein A with its *straight* descriptor and protein B with its *flipped* descriptor. (An ablation would symmetrise this: feed both straight+flipped into both sides and let the model learn which to use; it is cheap to try.)

### 4.2 Intra-protein geometric encoder

Build a directed k-nearest-neighbour graph (k ≈ 16) on each protein's 3D coordinates. Run L ≈ 4 layers of an SE(3)-invariant graph transformer with the following edge featurisation for an edge (i → j):

- Gaussian radial-basis expansion of `‖x_i − x_j‖` into, say, 16 bins.
- `cos(n_i, n_j)`, `cos(n_i, r_ij/‖r_ij‖)`, `cos(n_j, r_ij/‖r_ij‖)` — three rotation-invariant scalars that together capture the local relative pose of two surface points modulo a global rotation.

Initial node feature `h_i^0 = MLP([D_P[i]; F_P[i]])`. Each layer updates nodes via multi-head attention with the edge features added into both the attention logits and the values (GATv2 / Graphormer-style). The output is a set of per-vertex embeddings `h^A ∈ R^{M×d}`, `h^B ∈ R^{M×d}` (d ≈ 128) that now encode *descriptor identity in the context of local surface geometry*. Crucially every feature used is an SE(3) invariant, so independent rigid motion of A or B leaves `h^A, h^B` bit-identical.

### 4.3 Soft cross-protein compatibility

Compute an initial match-score matrix

`S_ij = (h^A_i · W_c · h^B_j) / √d`

and apply a *dual softmax* (softmax over rows times softmax over columns, elementwise). Dual softmax gives a sharp but differentiable approximation to mutual-nearest-neighbour matching: `P_ij ≈ 1` only if `j` is B's best match for `i` *and* `i` is A's best match for `j`. This is the learned analogue of the `desc_straight(A)`–`desc_flipped(B)` nearest-neighbour step that MaSIF does by raw descriptor distance today, but now reasoning in the geometrically-contextualised embedding space.

### 4.4 Candidate match set

Keep the top-K entries of `P` (K ≈ 256; trade-off parameter) and treat each selected pair `(i, j)` as a node in a new *match graph*. Node feature:

`φ_ij = MLP([h^A_i; h^B_j; P_ij; ‖D^s_A[i] − D^f_B[j]‖])`

Top-K can be made differentiable via straight-through, Gumbel-topK, or simply by backpropagating through the dense match-graph and masking — the exact choice is an engineering detail.

### 4.5 Relational match-graph transformer — the core mechanism

This is where "geometric consistency" becomes a learned, not hard-coded, property. On the match graph, run T ≈ 4 layers of attention. The key point: the attention between two match-nodes `n_ij` and `n_i'j'` has access to a rich edge feature built from intra-protein relational quantities that are themselves SE(3) invariants:

- RBF expansion of `d_A := ‖x_A[i] − x_A[i']‖` and `d_B := ‖x_B[j] − x_B[j']‖`, separately.
- RBF expansion of `|d_A − d_B|` — surfaced as a *feature*, not a threshold.
- Normal-geometry relationals: `cos(n_A[i], n_A[i'])` and `cos(n_B[j], n_B[j'])`, plus their difference.
- Optionally a higher-order invariant like `(n_A[i] × r_A) · n_A[i']` compared to the matching quantity in B, which captures chirality / surface orientation.

These features are the raw material from which the classical "two good matches are consistent iff `d_A[i,i'] ≈ d_B[j,j']`" heuristic is built, but here the network sees them as numeric inputs and learns what to do with them. It can discover that small `|d_A − d_B|` is strong evidence, but also that, say, slight mismatches are tolerable when normal-geometry is consistent, or that a small cluster of mutually-consistent matches outweighs a larger cluster of scattered matches — flexibility that a hard threshold cannot express.

Attention mechanics per layer:

`a(n_ij, n_{i'j'}) = softmax_{i'j'} [ (Q φ_ij) · (K φ_{i'j'}) + MLP_b(edge_{ij,i'j'}) ]`

`φ_ij ← φ_ij + Σ_{i'j'} a · (V φ_{i'j'} + MLP_v(edge_{ij,i'j'}))`

Because the only geometric quantities that enter are pairwise distances and inner products of unit normals / unit edge directions, the whole pipeline is manifestly invariant to independent rigid transforms of A and B.

The reason this works as a "learned coherence check": if a match `(i,j)` is real, the network can find many other high-confidence matches `(i',j')` whose intra-protein relationships on A match the corresponding relationships on B; messages arrive along those edges with consistent signals; the node embedding `φ_ij` gets amplified. A random, scattered descriptor match has no such chorus of corroborating matches and is attenuated.

### 4.6 Global readout

Two heads, combined:

- A *gate* `g_ij = σ(MLP_g(φ_ij))` — the probability that this match is part of a coherent interface.
- A *feature pool* `z = Σ g_ij φ_ij / (ε + Σ g_ij)` — the learned summary of the coherent subset.

The scalar interaction score is `ŷ = σ(MLP_o([z; Σ g_ij; Σ P_ij]))`. The `Σ g_ij` component behaves like a learned "interface size" statistic, and exposing it explicitly helps the model calibrate high scores to require a non-trivial number of coherent matches rather than relying on a single very strong one.

## 5. How the model satisfies each requirement

**Permutation invariance.** Every operation is a set operation (attention, sum pool, FPS followed by order-invariant layers, MLPs on individual nodes). Vertex order is never used.

**Independent rigid-motion invariance.** Coordinates enter only through pairwise distances and through cosines between normals and edge directions. An independent rotation and translation of protein A leaves `d_A`, normal cosines, and therefore every feature in the model unchanged — and symmetrically for B. This is invariance *by construction*, not by data augmentation. You should still use independent random SE(3) augmentation at training time as a numerical sanity check.

**Cross-protein relationships.** These live in the match graph. A node is a candidate correspondence; its existence is gated by descriptor + geometric-context compatibility; its embedding is refined by attention across *other* candidate correspondences, weighted by how relationally consistent each other correspondence is. Cross-protein reasoning is therefore explicit and relational, but never aligned.

**Spatial structure encoded invariantly.** Two levels: within each protein, the k-NN graph transformer injects local geometry into `h^A, h^B` using invariant edge features; across proteins, the match-graph transformer compares *intra-protein* distance and normal structure between pairs of matches. No global frame is ever chosen.

**Coherent vs. random matches.** Intuition: the match-graph transformer is a learned second-order correspondence verifier. A random set of descriptor-similar matches has near-zero mutual information in the edge features `|d_A − d_B|`, normal-relational features, etc., so attention provides no consistent amplification signal and the gates close. A real interface produces a dense clique in the match graph where the edge features are systematically small / consistent, giving a chorus of mutually reinforcing attention updates. The gate and the `Σ g_ij` term then both fire.

**No explicit alignment.** No orthogonal Procrustes, no ICP, no RANSAC, no pose hypotheses. The only geometric quantities ever computed are SE(3) invariants. The model is free to treat two descriptor-compatible matches as consistent on the basis of `|d_A − d_B|` ≈ 0, but nothing forces it to; it could equally learn that normal-relation consistency matters more, or that chains of triples matter, etc.

**Scalability.** Per-protein encoder is O(M · k) = O(512 · 16) per layer. Cross match computation is O(M²) = O(2.6 × 10⁵) once. Match-graph transformer is O(K²) per layer = O(6.5 × 10⁴) for K = 256, easily manageable. If you want K larger, you can drop to linear-attention or Performer kernels; the edge-bias MLP can be low-rank factorised.

## 6. Training

Positives come from the existing MaSIF PPI training lists (`masif/data/masif_ppi_search/lists/`): cognate binder chain pairs from co-crystal structures. Negatives: cross-complex pairings that share rough size statistics, plus hard negatives mined by running the trained descriptor module itself (top-scoring non-cognate pairs). Augment by applying independent random SE(3) transforms to each protein every batch. Loss: binary cross-entropy on `ŷ`.

Two useful auxiliary losses:

- A per-match auxiliary BCE on `g_ij` against "is the MaSIF-site interface label positive for both `i` and `j`" when available (`p{1,2}_sc_labels.npy`). This shapes the gate to point at genuine interface residues.
- Contrastive loss on the embeddings `h^A, h^B` after the intra-protein encoder, anchored by descriptor compatibility from true interface pairs, so that the encoder does not degenerate away from MaSIF's descriptor geometry.

## 7. Ablations that will validate the design claims

The design is opinionated; ablate to check each claim:

- Remove the relational edge features from the match-graph transformer (keep only node-node attention). The model should lose the ability to distinguish scattered from coherent matches — it now has descriptor compatibility information but no geometric-consistency signal. Expect a clear drop on the hardest negatives.
- Scramble the 3D coordinates of B while keeping its descriptors. A faithful model should collapse to chance or nearly so; a model leaking through descriptors alone would not. This isolates the contribution of geometry.
- Independently rotate each protein at inference. Predictions must be bit-identical (up to numerical noise). This is the invariance test.
- Replace the match-graph transformer with SuperGlue-style dense cross-attention only (no intra-protein distance edges). Expect intermediate performance — the model has no explicit machinery for relational consistency but can rediscover some of it.
- Swap FPS for random sampling. Expect modest but real degradation on small interfaces where FPS is more likely to land probes on them.

## 8. Concrete implementation path against MaSIF's on-disk data

1. From the precomputation directory, load the `.ply` mesh with `trimesh` or `pymesh`; extract `mesh.vertices` (`X`), `mesh.vertex_normals` (`N`), and the centre-vertex features from `p1_input_feat.npy[:, 0, :]` (`F`).
2. Load `p1_desc_straight.npy` (`D^s`) and `p1_desc_flipped.npy` (`D^f`). Verify shape `(V, 80)`.
3. Do the same for `p2_*`.
4. Build an in-memory `Protein` object with those six tensors. Apply FPS to 512 points; index all tensors accordingly.
5. Training loop consumes pairs `(Protein_A, Protein_B, y)` from the MaSIF training lists, with optional hard-negative sampling.
6. Implementation in PyTorch: the intra-protein encoder is a standard graph transformer (e.g. `torch_geometric` has ready-made layers); dual softmax and the match-graph transformer are ~200 lines of code.

## 9. Why this is the right shape for the problem

The central question you are asking is second-order: not "are these two patches compatible?" (first-order, MaSIF descriptors already answer that) but "are these compatibilities themselves mutually arranged in a way that is consistent with a single interface?" Second-order correspondence problems have a long and productive history in computer vision under various names — spectral matching, graph matching networks, SuperGlue — and in all of them the representation that works is exactly the one proposed here: nodes are candidate correspondences, edges are pairwise-relation consistency features, and attention learns how to weight those consistencies. The innovation specific to MaSIF-Coherence is using SE(3)-invariant relational edge features derived from the two proteins' intra-distance and normal geometry, so that the whole system is invariant by construction and no pose is ever chosen.

## 10. Framework and layout (v1)

The v1 baseline is implemented in **PyTorch**, not in TensorFlow 1 like the existing `masif_ppi_search` stack. The new module only consumes on-disk `.ply` meshes and `.npy` descriptor files produced by the existing MaSIF preprocessing pipeline, so it does not need to share a graph with any TF1 code. PyTorch also matches the assumptions of `minimal_v1_plan.md` and keeps the baseline under ~300 lines without dragging along the legacy TF1 / Python 3.7 environment.

Directory layout:

- Code: `masif_seed_search/source/masif_coherence/` — `config.py`, `data.py`, `model.py`, `train.py`, `check_invariance.py`.
- Data, lists, and trained checkpoints: `masif_seed_search/data/masif_coherence/` — `lists/`, `nn_models/v{n}/custom_params.py`, `nn_models/v{n}/model_data/`, plus a one-line `train_nn.sh` runner mirroring the one under `masif/data/masif_ppi_search/`.
- Training and testing pair lists are initially re-used verbatim from `masif/data/masif_ppi_search/lists/{training,testing}.txt`, so v1 is directly comparable to MaSIF-search on the same split.

Dependencies live in a self-contained `masif_seed_search/source/masif_coherence/requirements.txt` (`torch`, `numpy`, `trimesh`, `scikit-learn`); installing this module does not modify any existing MaSIF environment.
