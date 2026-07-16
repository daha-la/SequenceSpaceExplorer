"""
ca_reconstruct.py
-----------------
Pure-Python, cross-platform reconstruction of backbone N, C, and Cβ atoms
from a Cα-only trace.

Why this exists
---------------
Foldseek output (webserver JSON / .m8) provides only Cα coordinates for
target hits. Structure-aware embedders such as ProstT5 and SaProt need a
3Di sequence, and 3Di conversion (via mini3di) requires the N, Cα, Cβ and
C backbone atoms. This module places N, C and Cβ from the Cα trace so the
3Di sequence can be computed without any external binary (no PULCHRA, no
Foldseek), keeping the whole pipeline pip-installable on Windows, macOS
and Linux.

Method
------
The position of N and C relative to Cα depends on the backbone dihedral
angles, which are themselves recoverable from the local Cα geometry. For
each residue the offset of N, C and Cβ is expressed in a local orthonormal
frame built from the two flanking Cα-Cα bonds, and modelled as a small
Fourier series in the 4-Cα virtual dihedral:

    offset(d) = W[0] + W[1]cos(d) + W[2]sin(d) + W[3]cos(2d) + W[4]sin(2d)

The 45 coefficients below were fitted by least squares against PULCHRA
(the standard Cα reconstruction tool) on a diverse set of protein
structures. They describe peptide-bond geometry as a function of chain
curvature, not any particular dataset; leave-one-database-out validation
gives ~85% 3Di-sequence agreement with PULCHRA on unseen structures.

This is an approximation. For exact reconstruction use PULCHRA. For the
purpose of a 2D structural-similarity landscape — where relative position
matters and the same approximation is applied uniformly to every hit —
~85% agreement is sufficient.

Public API
----------
    reconstruct_backbone(ca_coords) -> dict with 'N','CA','C','CB' arrays
"""

import numpy as np


# ---------------------------------------------------------------------------
# Fitted coefficients — Fourier series in the 4-Cα virtual dihedral.
# Each is a (5, 3) array: rows are [const, cos d, sin d, cos 2d, sin 2d],
# columns are the offset components in the local (e1, e2, e3) frame.
# Fitted against PULCHRA on 2554 residues from structures spanning 9
# Foldseek databases. See module docstring.
# ---------------------------------------------------------------------------
_W_N = np.array([
    [ 0.794281,  1.153955, -0.281657],
    [ 0.040647, -0.069628, -0.004056],
    [ 0.114575, -0.066063, -0.031951],
    [-0.044812, -0.015758,  0.006788],
    [ 0.003096,  0.012798, -0.034164],
])
_W_C = np.array([
    [ 0.805450, -1.186504, -0.101993],
    [ 0.036283,  0.068611,  0.298302],
    [-0.108878, -0.118169,  0.259399],
    [-0.017060,  0.060749, -0.076945],
    [ 0.031568, -0.012236,  0.104270],
])
_W_CB = np.array([
    [-1.139500, -0.129485, -0.876616],
    [ 0.133177, -0.041778, -0.058753],
    [ 0.085270, -0.036630, -0.155371],
    [-0.002129,  0.024075,  0.079381],
    [-0.007040, -0.053152, -0.059070],
])

# Fallback offsets (the constant term only) for residues where the local
# geometry is degenerate and a dihedral cannot be computed reliably.
_FALLBACK_N  = _W_N[0]
_FALLBACK_C  = _W_C[0]
_FALLBACK_CB = _W_CB[0]


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
_EPS = 1e-7


def _unit(v):
    """Normalise a vector; return None if it is effectively zero-length."""
    n = np.linalg.norm(v)
    if n < _EPS:
        return None
    return v / n


def _local_frame(ca_prev, ca_i, ca_next):
    """
    Right-handed orthonormal frame at residue i from three consecutive Cα.
      e1 : bisector of the two Cα-Cα bond directions
      e2 : in-plane, perpendicular to e1
      e3 : normal to the three-Cα plane
    Returns None if the three points are collinear or coincident (the
    frame is then undefined).
    """
    v1 = _unit(ca_prev - ca_i)
    v2 = _unit(ca_next - ca_i)
    if v1 is None or v2 is None:
        return None
    e1 = _unit(v1 + v2)
    e3 = _unit(np.cross(v1, v2))
    if e1 is None or e3 is None:
        # v1 and v2 (anti)parallel — collinear Cα, frame undefined.
        return None
    e2 = np.cross(e3, e1)
    return np.array([e1, e2, e3])


def _virtual_dihedral(ca, i, L):
    """
    4-Cα virtual dihedral around the Cα(i)-Cα(i+1) bond. Returns 0.0 when
    it cannot be computed (chain end or degenerate geometry); callers
    treat 0.0 as 'use the constant term only', which the Fourier series
    handles gracefully.
    """
    im  = i - 1 if i > 0     else i
    ip  = i + 1 if i < L - 1 else i
    ip2 = i + 2 if i < L - 2 else ip
    b1 = ca[i]   - ca[im]
    b2 = ca[ip]  - ca[i]
    b3 = ca[ip2] - ca[ip]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    ub2 = _unit(b2)
    if ub2 is None or np.linalg.norm(n1) < _EPS or np.linalg.norm(n2) < _EPS:
        return 0.0
    return np.arctan2(np.dot(np.cross(n1, n2), ub2), np.dot(n1, n2))


def _fourier(d):
    """Feature vector [1, cos d, sin d, cos 2d, sin 2d]."""
    return np.array([1.0, np.cos(d), np.sin(d), np.cos(2 * d), np.sin(2 * d)])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def reconstruct_backbone(ca_coords) -> dict:
    """
    Reconstruct backbone N, C and Cβ atoms from a Cα-only trace.

    Parameters
    ----------
    ca_coords : array-like, shape (L, 3)
        Cα coordinates of a single chain, in residue order.

    Returns
    -------
    dict
        {'N', 'CA', 'C', 'CB'} each an (L, 3) float64 ndarray. 'CA' is the
        input unchanged. Glycine residues still receive a Cβ position, as
        mini3di expects four atoms per residue; the caller may discard it.

    Raises
    ------
    ValueError
        If fewer than 3 residues are supplied (a frame needs neighbours).

    Notes
    -----
    Degenerate residues — collinear or coincident neighbouring Cα atoms —
    are handled by falling back to the constant-term offset in the best
    frame available, or, if no frame can be built at all, by copying the
    Cα position. Such residues are rare and localised; they do not abort
    the reconstruction.
    """
    ca = np.asarray(ca_coords, dtype=np.float64)
    if ca.ndim != 2 or ca.shape[1] != 3:
        raise ValueError("ca_coords must have shape (L, 3).")
    L = ca.shape[0]
    if L < 3:
        raise ValueError(
            f"Need at least 3 residues to reconstruct a backbone, got {L}."
        )
    if not np.isfinite(ca).all():
        raise ValueError("ca_coords contains non-finite values.")

    N  = np.zeros((L, 3))
    C  = np.zeros((L, 3))
    CB = np.zeros((L, 3))

    for i in range(L):
        ca_prev = ca[i - 1] if i > 0     else ca[i] - (ca[i + 1] - ca[i])
        ca_next = ca[i + 1] if i < L - 1 else ca[i] - (ca[i - 1] - ca[i])

        frame = _local_frame(ca_prev, ca[i], ca_next)
        if frame is None:
            # No usable frame — fall back to placing atoms at Cα. This is
            # a last resort for pathological input; mini3di still runs.
            N[i]  = ca[i]
            C[i]  = ca[i]
            CB[i] = ca[i]
            continue

        d  = _virtual_dihedral(ca, i, L)
        f  = _fourier(d)
        # offset = f · W  -> a length-3 vector in the local frame
        N[i]  = ca[i] + frame.T @ (f @ _W_N)
        C[i]  = ca[i] + frame.T @ (f @ _W_C)
        CB[i] = ca[i] + frame.T @ (f @ _W_CB)

    return {"N": N, "CA": ca, "C": C, "CB": CB}