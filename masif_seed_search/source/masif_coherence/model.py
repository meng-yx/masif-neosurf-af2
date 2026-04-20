"""MaSIF-Coherence v1 model.

Implements the minimal six-step pipeline defined in `minimal_v1_plan.md`:

    1. Pairwise descriptor distance M[i, j] = ||D_A[i] - D_B[j]||.
    2. Top-K (default 256) smallest-distance matches.
    3. Per-match node feature phi_k = MLP([D_A[i]; D_B[j]; dist_k]).
    4. RBF-expanded intra-protein edge features for every ordered match pair:
       e_{kl} = [RBF(||X_A[i_k] - X_A[i_l]||); RBF(||X_B[j_k] - X_B[j_l]||)].
    5. Two-layer match-graph transformer with additive edge-bias attention.
    6. Mean pool of node features -> MLP -> scalar logit.

The model has ~250 LOC and a little under 1M parameters at the default config.
Invariance to independent rigid motions of A and B, and to vertex permutations,
is by construction: coordinates only ever enter via pairwise distances.
"""

from __future__ import annotations

from typing import Dict, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


# ---------------------------------------------------------------------------
# Small building blocks.
# ---------------------------------------------------------------------------


class GaussianRBF(nn.Module):
    """Non-learnable Gaussian radial-basis expansion over a fixed range."""

    def __init__(self, num_bins: int, d_min: float, d_max: float):
        super().__init__()
        centers = torch.linspace(d_min, d_max, num_bins)
        gamma = (d_max - d_min) / max(num_bins - 1, 1)
        self.register_buffer("centers", centers)
        self.gamma = float(gamma)

    def forward(self, d: torch.Tensor) -> torch.Tensor:
        # d: (...,) -> (..., num_bins)
        d = d.unsqueeze(-1)
        return torch.exp(-((d - self.centers) ** 2) / (self.gamma ** 2 + 1e-12))


class FeedForward(nn.Module):
    def __init__(self, dim: int, mult: int = 2, dropout: float = 0.0):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, dim * mult),
            nn.GELU(),
            nn.Linear(dim * mult, dim),
        )
        self.dropout = nn.Dropout(dropout) if dropout > 0 else nn.Identity()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.net(x))


# ---------------------------------------------------------------------------
# Match-graph transformer layer with additive edge bias.
# ---------------------------------------------------------------------------


class MatchGraphAttention(nn.Module):
    """Multi-head attention over match nodes with an additive edge-feature bias.

    Edge features `E` are expected shape (K, K, edge_dim). They contribute both
    an additive scalar bias to the attention logits and an additive vector to
    the value stream, following the recipe in section 4.5 / step 5 of the docs.
    """

    def __init__(self, dim: int, num_heads: int, edge_dim: int):
        super().__init__()
        if dim % num_heads != 0:
            raise ValueError(f"dim ({dim}) must be divisible by num_heads ({num_heads}).")
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.scale = self.head_dim ** -0.5

        self.qkv = nn.Linear(dim, 3 * dim, bias=False)
        self.out = nn.Linear(dim, dim)

        self.edge_bias = nn.Sequential(
            nn.Linear(edge_dim, num_heads),
        )
        self.edge_value = nn.Sequential(
            nn.Linear(edge_dim, dim),
            nn.GELU(),
            nn.Linear(dim, dim),
        )

    def forward(self, x: torch.Tensor, edge: torch.Tensor) -> torch.Tensor:
        # x: (K, dim); edge: (K, K, edge_dim)
        K = x.shape[0]
        qkv = self.qkv(x)  # (K, 3*dim)
        q, k, v = qkv.chunk(3, dim=-1)
        # Split heads: (K, H, head_dim)
        q = q.view(K, self.num_heads, self.head_dim)
        k = k.view(K, self.num_heads, self.head_dim)
        v = v.view(K, self.num_heads, self.head_dim)

        # Attention logits: (H, K, K)
        logits = torch.einsum("khd,lhd->hkl", q, k) * self.scale
        edge_b = self.edge_bias(edge)  # (K, K, H)
        logits = logits + edge_b.permute(2, 0, 1)  # (H, K, K)

        attn = torch.softmax(logits, dim=-1)
        val_edge = self.edge_value(edge)  # (K, K, dim)
        val_edge = val_edge.view(K, K, self.num_heads, self.head_dim)

        # out[k, h, d] = sum_l attn[h, k, l] * (v[l, h, d] + val_edge[k, l, h, d])
        out_from_v = torch.einsum("hkl,lhd->khd", attn, v)
        out_from_e = torch.einsum("hkl,klhd->khd", attn, val_edge)
        out = (out_from_v + out_from_e).reshape(K, self.num_heads * self.head_dim)
        return self.out(out)


class MatchGraphLayer(nn.Module):
    """Pre-norm transformer block: MGA + FFN, each with residual."""

    def __init__(self, dim: int, num_heads: int, edge_dim: int, ffn_mult: int = 2):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.attn = MatchGraphAttention(dim, num_heads, edge_dim)
        self.norm2 = nn.LayerNorm(dim)
        self.ffn = FeedForward(dim, mult=ffn_mult)

    def forward(self, x: torch.Tensor, edge: torch.Tensor) -> torch.Tensor:
        x = x + self.attn(self.norm1(x), edge)
        x = x + self.ffn(self.norm2(x))
        return x


# ---------------------------------------------------------------------------
# Top-level model.
# ---------------------------------------------------------------------------


class MaSIFCoherence(nn.Module):
    def __init__(
        self,
        descriptor_dim: int = 80,
        top_k: int = 256,
        hidden_dim: int = 64,
        num_heads: int = 4,
        num_layers: int = 2,
        ffn_mult: int = 2,
        rbf_num_bins: int = 16,
        rbf_min_ang: float = 0.0,
        rbf_max_ang: float = 40.0,
    ):
        super().__init__()
        self.top_k = int(top_k)
        self.hidden_dim = int(hidden_dim)

        node_in_dim = 2 * descriptor_dim + 1  # [D_A[i]; D_B[j]; dist]
        self.node_mlp = nn.Sequential(
            nn.Linear(node_in_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

        self.rbf = GaussianRBF(rbf_num_bins, rbf_min_ang, rbf_max_ang)
        edge_dim = 2 * rbf_num_bins  # [RBF(d_A); RBF(d_B)]

        self.layers = nn.ModuleList(
            [
                MatchGraphLayer(hidden_dim, num_heads, edge_dim, ffn_mult=ffn_mult)
                for _ in range(num_layers)
            ]
        )
        self.final_norm = nn.LayerNorm(hidden_dim)

        self.readout = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, 1),
        )

    # ------------------------------------------------------------------
    # Forward.
    # ------------------------------------------------------------------

    def forward(
        self,
        X_A: torch.Tensor,
        D_A: torch.Tensor,
        X_B: torch.Tensor,
        D_B: torch.Tensor,
        return_aux: bool = False,
    ):
        """Run the coherence pipeline on a single pair (batch size 1).

        Args:
            X_A: (V_A, 3) coordinates of protein A.
            D_A: (V_A, 80) MaSIF descriptors of protein A (straight).
            X_B: (V_B, 3) coordinates of protein B.
            D_B: (V_B, 80) MaSIF descriptors of protein B (flipped).
            return_aux: if True, also return selected indices and distances
                (useful for invariance checks and debugging).

        Returns:
            logit: scalar torch.Tensor. Apply torch.sigmoid for probability.
        """
        if X_A.ndim != 2 or X_A.shape[1] != 3:
            raise ValueError(f"X_A must be (V_A, 3), got {tuple(X_A.shape)}")
        if X_B.ndim != 2 or X_B.shape[1] != 3:
            raise ValueError(f"X_B must be (V_B, 3), got {tuple(X_B.shape)}")
        if D_A.shape[0] != X_A.shape[0]:
            raise ValueError(
                f"D_A has {D_A.shape[0]} rows but X_A has {X_A.shape[0]} vertices."
            )
        if D_B.shape[0] != X_B.shape[0]:
            raise ValueError(
                f"D_B has {D_B.shape[0]} rows but X_B has {X_B.shape[0]} vertices."
            )

        # Step 1: pairwise descriptor distance. cdist is deterministic and fine
        # for V_A * V_B up to a few thousand each.
        M = torch.cdist(D_A, D_B, p=2)  # (V_A, V_B)

        # Step 2: top-K smallest entries of the flattened distance matrix.
        K = min(self.top_k, M.numel())
        flat = M.reshape(-1)
        neg_topk = torch.topk(-flat, k=K, sorted=False)
        flat_idx = neg_topk.indices
        dist_k = -neg_topk.values
        V_B = M.shape[1]
        idx_A = torch.div(flat_idx, V_B, rounding_mode="floor")
        idx_B = flat_idx % V_B

        # Step 3: per-match node features.
        da = D_A[idx_A]  # (K, 80)
        db = D_B[idx_B]  # (K, 80)
        node_in = torch.cat([da, db, dist_k.unsqueeze(-1)], dim=-1)  # (K, 161)
        phi = self.node_mlp(node_in)  # (K, hidden_dim)

        # Step 4: intra-protein RBF edge features.
        xa_sel = X_A[idx_A]  # (K, 3)
        xb_sel = X_B[idx_B]  # (K, 3)
        dA = torch.cdist(xa_sel, xa_sel, p=2)  # (K, K)
        dB = torch.cdist(xb_sel, xb_sel, p=2)  # (K, K)
        edge = torch.cat([self.rbf(dA), self.rbf(dB)], dim=-1)  # (K, K, 2*bins)

        # Step 5: match-graph transformer.
        for layer in self.layers:
            phi = layer(phi, edge)
        phi = self.final_norm(phi)

        # Step 6: mean-pool readout to a scalar logit.
        z = phi.mean(dim=0)  # (hidden_dim,)
        logit = self.readout(z).squeeze(-1)

        if return_aux:
            aux: Dict[str, torch.Tensor] = {
                "idx_A": idx_A.detach(),
                "idx_B": idx_B.detach(),
                "dist_k": dist_k.detach(),
                "phi": phi.detach(),
            }
            return logit, aux
        return logit


def build_model_from_params(params: dict) -> MaSIFCoherence:
    """Instantiate a MaSIFCoherence using keys from the config dict."""
    return MaSIFCoherence(
        descriptor_dim=params.get("descriptor_dim", 80),
        top_k=params["top_k"],
        hidden_dim=params["hidden_dim"],
        num_heads=params["num_heads"],
        num_layers=params["num_layers"],
        ffn_mult=params.get("ffn_mult", 2),
        rbf_num_bins=params["rbf_num_bins"],
        rbf_min_ang=params["rbf_min_ang"],
        rbf_max_ang=params["rbf_max_ang"],
    )
