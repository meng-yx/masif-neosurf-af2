# MaSIF-Coherence v1 — minimal baseline

The hypothesis to test: **"a relational transformer over cross-protein match candidates, with SE(3)-invariant intra-protein distance edges, learns to tell coherent from scattered descriptor matches."** That is the only load-bearing claim. Everything else in the earlier design is there to *improve* the model once the baseline works, not to make the baseline work. v1 cuts every moving part that is not required to test the hypothesis.

## What v1 keeps and what it drops

Kept: (1) raw MaSIF descriptors as per-vertex features, (2) a simple top-K cross-protein match selection by descriptor distance, (3) the match-graph transformer with intra-protein distance edges, (4) a mean-pool + MLP readout. That is four components.

Dropped (moved to ablations): per-protein graph-transformer encoder, dual-softmax matching, gated pooling, auxiliary losses on per-match interface labels, normals and the 5-channel input features, FPS subsampling, chirality / higher-order invariants, independent SE(3) augmentation (invariance is already guaranteed by construction; training with it is belt-and-braces).

The reasoning for each drop: the MaSIF descriptor already integrates local patch geometry via its polar patch encoder, so a per-protein encoder is redundant for a first-pass signal; dual softmax is a refinement of top-K with similar behaviour at small K; gates, auxiliary losses and extra features all help calibration but cannot be the reason the model works or fails — they tune a working model.

## Minimal architecture

Inputs per training example: two proteins A, B and a binary label y ∈ {0, 1}.

Per protein: three tensors on disk, loaded verbatim.

| Name | Source file | Shape | Dtype |
|---|---|---|---|
| `X_A` | `.ply` mesh in precomputation dir, `mesh.vertices` | `(V_A, 3)` | float32 |
| `D_A` | `p1_desc_straight.npy` | `(V_A, 80)` | float32 |
| `X_B` | `.ply` mesh | `(V_B, 3)` | float32 |
| `D_B` | `p2_desc_flipped.npy` | `(V_B, 80)` | float32 |

Note: A uses straight descriptors, B uses flipped — this respects the MaSIF complementarity convention and is the only non-obvious data choice in v1.

Pipeline, step by step.

**Step 1 — descriptor distance matrix.**

`M[i, j] = ‖D_A[i] − D_B[j]‖₂`,    shape `(V_A, V_B)`, float32.

**Step 2 — top-K match selection.**

Flatten `M`, take indices of the `K = 256` smallest entries. Produce three tensors:

| Name | Shape | Notes |
|---|---|---|
| `idx_A` | `(K,)` int64 | row indices into A |
| `idx_B` | `(K,)` int64 | col indices into B |
| `dist` | `(K,)` float32 | descriptor distance for each selected pair |

This step is non-differentiable, but training only needs gradients through *values* at those indices, which flow fine. Top-K is a hard selector in v1; do not attempt Gumbel / straight-through until you have a working baseline to ablate against.

**Step 3 — per-match node features.**

For each of the K candidate matches k:

`φ_k = MLP_node( [ D_A[idx_A[k]] ; D_B[idx_B[k]] ; dist[k] ] )`

Input dim 80 + 80 + 1 = 161. Output dim `d = 64`. MLP is 2 layers, GELU, no dropout in v1.

Stack into `Φ ∈ R^{K×d}`.

**Step 4 — intra-protein relational edge features.**

For every ordered pair of matches (k, l):

`d^A_{kl} = ‖X_A[idx_A[k]] − X_A[idx_A[l]]‖₂`,  `d^B_{kl} = ‖X_B[idx_B[k]] − X_B[idx_B[l]]‖₂`

Gaussian radial basis expansion with 16 bins between 0 and 40 Å for each of `d^A` and `d^B`:

`e_{kl} = [ RBF(d^A_{kl}) ; RBF(d^B_{kl}) ] ∈ R^{32}`

Edge tensor `E ∈ R^{K×K×32}`, float32. Memory for K=256 is 256·256·32·4 ≈ 8 MB per pair — fine.

The important design decision: **the model sees the two distances `d^A` and `d^B` as separate features, not their difference**. The model is free to learn that the difference matters (it will), but we do not pre-bake that into the representation. This is what makes "geometric consistency" a learned concept rather than a hard-coded one.

**Step 5 — match-graph transformer.**

Two layers. Each layer: multi-head attention over the K match-nodes with an additive edge bias, then a 2-layer feed-forward block. Pre-norm. Heads = 4.

For head h at layer ℓ:

`A_{kl}^h = ( (Q^h φ_k) · (K^h φ_l) ) / √d_h  +  b^h · MLP_edge^h(e_{kl})`

`φ_k ← φ_k + Σ_l softmax_l(A_{kl}^h) · ( V^h φ_l + MLP_v^h(e_{kl}) )`

then the standard FFN block. Only 2 layers in v1 — enough to propagate "I'm consistent with several high-quality neighbours" two hops but no more.

**Step 6 — readout.**

`z = mean_k(φ_k) ∈ R^d`
`ŷ = σ( MLP_out(z) )`,  scalar in (0, 1).

Mean pool is used rather than gated pool. Gates are an ablation.

That's the whole model. Six steps, ~250 lines of PyTorch.

## Parameter count sanity

With `d = 64`, 4 heads, 2 layers, the model is well under 1M parameters. It is small on purpose: v1 is a hypothesis test, not a competitive benchmark model.

## Training setup (minimal)

- Optimiser: Adam, lr 1e-4.
- Loss: `BCEWithLogits(ŷ_logits, y)`. No auxiliary losses.
- Batch size: 1 pair per gradient step; accumulate over 16 if you want smoother gradients.
- Positives: cognate binder pairs from `masif/data/masif_ppi_search/lists/training_ppi.txt` (or the equivalent list shipped with the repo).
- Negatives: one random non-cognate pair per positive, drawn uniformly from the same list. No hard-negative mining in v1.
- Epochs: train until validation AUROC plateaus; budget a few hundred gradient steps for a first signal.
- Data augmentation: none. Invariance is structural.

## Invariance by construction — verify numerically before training

Before any training, run this check once:

1. Randomly rotate and translate `X_A` independently; verify `ŷ` is identical to 1e-5.
2. Same for `X_B`.
3. Randomly permute the order of A's vertices (and re-permute `D_A` accordingly); verify `ŷ` is identical.

If any of these fails, there is a bug. This is cheap and catches almost all invariance regressions.

## Smallest useful feature set — and what is deliberately absent

Per match node: `[D_A[i] ; D_B[j] ; ‖D_A[i] − D_B[j]‖]`. 161 scalars. That is all.

Per match edge: `[RBF(d^A_{kl}) ; RBF(d^B_{kl})]`. 32 scalars.

Absent on purpose: normals, per-vertex chemical features, explicit `|d^A − d^B|`, dot products of edge vectors, surface-site labels.

Rationale for exclusion: if the model cannot distinguish coherent from random matches using just descriptor similarity and intra-protein distances, no amount of normal-cosine features will fix it — they would only mask a broken hypothesis. Conversely, if the baseline works, these features are the first things to add.

## Components moved to ablations

Listed in the order I would add them after the baseline is working, each as a single, isolated change vs. v1:

1. **Replace mean pool with gated pool** — `g_k = σ(MLP_gate(φ_k)); z = Σ g_k φ_k / (ε + Σ g_k)` and concatenate `Σ g_k` into the readout input. Expected to help calibration.
2. **Add normal cosines to edges** — append `cos(n_A[i_k], n_A[i_l])` and `cos(n_B[j_k], n_B[j_l])` to `e_{kl}`. Requires loading `mesh.vertex_normals`. Expected to help on surfaces where distance alone is ambiguous.
3. **Add `|d^A − d^B|` RBF as an edge feature** — a prior that biases learning. Expected to improve sample efficiency, not asymptotic accuracy.
4. **Dual-softmax match selection** — replaces hard top-K by a soft row-and-column-mutual-nearest step. Expected to help when genuine matches are not top-K by raw descriptor distance.
5. **Per-protein graph-transformer encoder** — replaces raw `D_A, D_B` with geometry-aware embeddings before Step 1. Expected to help more on harder, flexible-interface cases.
6. **Auxiliary per-match BCE against `p{1,2}_sc_labels.npy`** — shapes the model to light up on true interface residues.
7. **Hard-negative mining** — mine negatives using the v1 model itself.
8. **FPS subsampling of big proteins** — only matters if `V > ~2000` and memory becomes a concern; v1 at K=256 is indifferent to `V` beyond the O(V·80) descriptor load and the one O(V_A · V_B) distance matrix.
9. **More layers or more K** — scale the match-graph transformer to 4 layers or K=512 only if needed.

Each ablation should be evaluated against v1 on the same validation split. If an ablation does not move AUROC by a clear margin, it is not justified.

## The decision v1 commits to

v1 is a bet that the following three choices are individually responsible for most of the signal: using top-K of raw descriptor distance as the match candidate set, feeding raw `D_A`/`D_B`/distance into match nodes, and biasing attention with RBFs of the two intra-protein distances separately. If v1 does not clear a meaningful AUROC over a descriptor-only baseline, the architecture is wrong in some way no ablation will rescue, and the design should be rethought rather than patched. If v1 *does* clear that bar, the ablations above are the disciplined path to the full model.
