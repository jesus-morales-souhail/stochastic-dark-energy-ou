#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Seal the tesseract no-go on first-order optical phase space:

  |B4| = 4! * 2^4 = 384 signed permutations of R^4.
  Optical symplectic form: J = blockdiag(J2, J2) on (x,p_x,y,p_y).
  Exactly 32 of the 384 preserve J (M J M^T = J).
  Their joint commutant inside sp(4,R) (dim 10) has dimension 1.

So B4-as-symmetry *kills* 9 of 10 first-order freedoms; it does not
amplify channels. "8 cells → 8 optical channels" is polytope combinatorics,
not the rank of a Maxwell-derived operator.

Run:
  python scripts/b4_symplectic_count.py
"""

from __future__ import annotations

import itertools

import numpy as np

J2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
J = np.block([[J2, np.zeros((2, 2))], [np.zeros((2, 2)), J2]])


def signed_perms():
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((-1.0, 1.0), repeat=4):
            M = np.zeros((4, 4))
            for i, j in enumerate(perm):
                M[i, j] = signs[i]
            yield M


def is_symplectic(M: np.ndarray, tol: float = 1e-10) -> bool:
    return np.allclose(M @ J @ M.T, J, atol=tol)


def sp4_basis() -> list[np.ndarray]:
    """Standard 10-dimensional basis of sp(4,R): A=[[X,Y],[Z,-X^T]], Y=Y^T, Z=Z^T."""
    basis: list[np.ndarray] = []
    for i in range(2):
        for j in range(2):
            X = np.zeros((2, 2))
            X[i, j] = 1.0
            basis.append(np.block([[X, np.zeros((2, 2))], [np.zeros((2, 2)), -X.T]]))
    for i in range(2):
        for j in range(i, 2):
            Y = np.zeros((2, 2))
            Y[i, j] = Y[j, i] = 1.0
            basis.append(
                np.block(
                    [[np.zeros((2, 2)), Y], [np.zeros((2, 2)), np.zeros((2, 2))]]
                )
            )
    for i in range(2):
        for j in range(i, 2):
            Z = np.zeros((2, 2))
            Z[i, j] = Z[j, i] = 1.0
            basis.append(
                np.block(
                    [[np.zeros((2, 2)), np.zeros((2, 2))], [Z, np.zeros((2, 2))]]
                )
            )
    return basis


def main() -> None:
    total = 24 * 16
    Gs = [M for M in signed_perms() if is_symplectic(M)]
    print("=" * 64)
    print("B4 signed permutations vs optical symplectic form")
    print("=" * 64)
    print(f"|B4| = 4! × 2^4 = {total}")
    print(f"Symplectic (M J M^T = J): {len(Gs)}")
    print(f"Expected 32 (4×4×2): {len(Gs) == 32}")

    basis = sp4_basis()
    print(f"dim sp(4,R) basis: {len(basis)}")

    rows = []
    for G in Gs:
        block = np.column_stack([(G @ B - B @ G).ravel() for B in basis])
        rows.append(block)
    Mcond = np.vstack(rows)
    s = np.linalg.svd(Mcond, compute_uv=False)
    rank = int(np.sum(s > 1e-8))
    nullity = Mcond.shape[1] - rank
    print(f"Commutant dim in sp(4,R): {nullity}  (expected 1: {nullity == 1})")
    print()
    print("Reading:")
    print("  B4 geometry, as a *symmetry*, restricts first-order optical maps.")
    print("  It does not add SLM channels. Controllable DOF = N_pix [× N_ω], hardware.")
    print("  Maxwell+paraxial → Fresnel + phase screens → Sp(4,R) for quadratic Φ.")
    print("  Tesseract combinatorics (8 cells) ≠ rank of that operator.")
    print("=" * 64)


if __name__ == "__main__":
    main()
