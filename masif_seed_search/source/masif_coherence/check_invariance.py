"""Numerical invariance sanity check for MaSIF-Coherence v1.

Per section 6 of `minimal_v1_plan.md`, before any training is done we must
verify that the model is invariant to:

    (i)   independent random rigid transforms of X_A,
    (ii)  independent random rigid transforms of X_B,
    (iii) random vertex permutations of A (with D_A permuted the same way).

All three assertions should pass to ~1e-5. A failure means a bug, not a
training artefact.

Runs on synthetic random inputs, so it does not require any real MaSIF
precomputation directory.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import torch

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

from config import get_default_params  # noqa: E402
from model import build_model_from_params  # noqa: E402


def _random_rigid(seed: int) -> tuple[torch.Tensor, torch.Tensor]:
    """Return (R, t) with R a random SO(3) rotation and t a random translation."""
    gen = torch.Generator().manual_seed(seed)
    # Random rotation via QR of a Gaussian matrix; fix sign to keep det = +1.
    A = torch.randn(3, 3, generator=gen)
    Q, R = torch.linalg.qr(A)
    d = torch.sign(torch.diagonal(R))
    d[d == 0] = 1.0
    Q = Q * d.unsqueeze(0)
    if torch.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    t = torch.randn(3, generator=gen) * 10.0
    return Q, t


def _apply_rigid(X: torch.Tensor, R: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
    return X @ R.T + t


def _make_inputs(V_A: int, V_B: int, seed: int = 0):
    gen = torch.Generator().manual_seed(seed)
    X_A = torch.randn(V_A, 3, generator=gen) * 10.0
    X_B = torch.randn(V_B, 3, generator=gen) * 10.0
    D_A = torch.randn(V_A, 80, generator=gen)
    D_B = torch.randn(V_B, 80, generator=gen)
    return X_A, D_A, X_B, D_B


def run_checks(tol: float = 1e-5, V_A: int = 400, V_B: int = 420) -> None:
    params = get_default_params()
    model = build_model_from_params(params).eval()

    X_A, D_A, X_B, D_B = _make_inputs(V_A, V_B, seed=0)

    with torch.no_grad():
        base = model(X_A, D_A, X_B, D_B).item()

    def run(X_A2, X_B2, D_A2=D_A, D_B2=D_B):
        with torch.no_grad():
            return model(X_A2, D_A2, X_B2, D_B2).item()

    # (i) rigid transform on A only.
    R_A, t_A = _random_rigid(seed=1)
    y_A = run(_apply_rigid(X_A, R_A, t_A), X_B)
    print(f"Check (i)   A rotated+translated:   base={base:.8f}  new={y_A:.8f}  diff={abs(y_A - base):.2e}")
    assert abs(y_A - base) < tol, f"A rigid-transform invariance failed ({abs(y_A - base):.2e})"

    # (ii) rigid transform on B only.
    R_B, t_B = _random_rigid(seed=2)
    y_B = run(X_A, _apply_rigid(X_B, R_B, t_B))
    print(f"Check (ii)  B rotated+translated:   base={base:.8f}  new={y_B:.8f}  diff={abs(y_B - base):.2e}")
    assert abs(y_B - base) < tol, f"B rigid-transform invariance failed ({abs(y_B - base):.2e})"

    # (i+ii) rigid transform on both.
    y_AB = run(_apply_rigid(X_A, R_A, t_A), _apply_rigid(X_B, R_B, t_B))
    print(f"Check (iii) A and B independently:  base={base:.8f}  new={y_AB:.8f}  diff={abs(y_AB - base):.2e}")
    assert abs(y_AB - base) < tol, f"Joint rigid-transform invariance failed ({abs(y_AB - base):.2e})"

    # (iv) permutation of A's vertex order (and D_A accordingly).
    perm = torch.randperm(V_A, generator=torch.Generator().manual_seed(3))
    y_perm = run(X_A[perm], X_B, D_A2=D_A[perm], D_B2=D_B)
    print(f"Check (iv)  A vertex permutation:   base={base:.8f}  new={y_perm:.8f}  diff={abs(y_perm - base):.2e}")
    assert abs(y_perm - base) < tol, f"A permutation invariance failed ({abs(y_perm - base):.2e})"

    # (v) permutation of B's vertex order.
    permB = torch.randperm(V_B, generator=torch.Generator().manual_seed(4))
    y_permB = run(X_A, X_B[permB], D_B2=D_B[permB])
    print(f"Check (v)   B vertex permutation:   base={base:.8f}  new={y_permB:.8f}  diff={abs(y_permB - base):.2e}")
    assert abs(y_permB - base) < tol, f"B permutation invariance failed ({abs(y_permB - base):.2e})"

    print("All invariance checks passed.")


if __name__ == "__main__":
    # Use double precision for a tighter numerical tolerance.
    torch.set_default_dtype(torch.float64)
    run_checks()
