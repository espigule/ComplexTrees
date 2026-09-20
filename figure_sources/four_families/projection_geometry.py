"""Exact algebraic charts; binary64 evaluation for the four-family figures.

The formulas themselves are exact. Floating-point mesh evaluation is a display
calculation, not a new stability or interval certificate. Stability is inherited
only from the original certified q-cells.
"""
from __future__ import annotations

import numpy as np

FAMILIES = ("DD10", "DD11", "DO10", "DO11")


def endpoint_maps(q, family: str):
    """Return (t1,A1,p1,t2,A2,p2); p=0 direct, p=1 conjugate."""
    if family not in FAMILIES:
        raise ValueError(f"Expected one of {FAMILIES}; got {family!r}")
    q = np.asarray(q, dtype=np.complex128)
    e2 = int(family[-1])
    return q, -q, 0, q + e2 * (1 - q), (-1) ** e2 * (1 - q), int(family[1] == "O")


def canonical_coefficients(q, family: str, *, validate: bool = True):
    """Return a,b for S1(z)=-1+a*z; S2(z)=1+b*J(z).

    The DO opposite coefficient includes the necessary rotation induced by
    canonical translation normalization. q is the endpoint zipper coordinate.
    """
    if family not in FAMILIES:
        raise ValueError(f"Expected one of {FAMILIES}; got {family!r}")
    q = np.asarray(q, dtype=np.complex128)
    if validate and np.any((abs(q) <= 0) | (abs(q) >= 1) | (abs(1-q) <= 0) | (abs(1-q) >= 1)):
        raise ValueError("q must lie in the strict lens |q|<1, |1-q|<1.")
    a = -q
    b = (-1) ** int(family[-1]) * (1 - q)
    if family.startswith("DO"):
        x, y = q.real, q.imag
        s = x * x + y * y
        if family == "DO10":
            N = 3*x - 2*y*y + 1j*y*(2*x - 1)
        else:
            N = 1 + q + s*(2*q - 1)
        b = b * np.conj(N) / N
    return a, b


def exterior_coefficients(q, family: str, *, validate: bool = True):
    """Canonical exterior coordinates c1=1/a,c2=1/b, with both moduli>1."""
    a, b = canonical_coefficients(q, family, validate=validate)
    return 1 / a, 1 / b


def endpoint_exterior(q):
    """Unsigned endpoint coordinate c=1/q, not canonical c1=-1/q."""
    return 1 / np.asarray(q, dtype=np.complex128)


def atlas_from_coefficients(a, b):
    """Return all 4 atlas coordinates and the two globally periodic phases.

    alpha in [0,pi), gamma in [0,2*pi) uses the manuscript phase lattice.
    Its seam is (0,gamma)~(pi,gamma+pi). The 3D display below evaluates the
    coefficient phases directly and therefore has no artificial mesh seam.
    """
    a, b = np.broadcast_arrays(np.asarray(a, dtype=np.complex128), np.asarray(b, dtype=np.complex128))
    r1, r2 = abs(a), abs(b)
    theta1, theta2 = np.angle(a), np.angle(b)
    alpha = (theta1-theta2)/2
    gamma = -(theta1+theta2)/2
    turn = np.floor(alpha / np.pi)
    alpha = alpha - turn*np.pi
    gamma = np.mod(gamma - turn*np.pi, 2*np.pi)
    return {"w": r1/(r1+r2), "rho": 2/(r1+r2), "alpha": alpha,
            "gamma": gamma, "theta1": theta1, "theta2": theta2}


def atlas_coordinates(q, family: str, *, validate: bool = True):
    a, b = canonical_coefficients(q, family, validate=validate)
    return atlas_from_coefficients(a, b)


def coefficient_torus(theta1, theta2, rho, *, major_radius: float = 4.0):
    """Global phase/radius display: T^2 x (0,R) -> R^3.

    X=(R+rho*cos(theta1))*cos(theta2),
    Y=(R+rho*cos(theta1))*sin(theta2), Z=rho*sin(theta1).

    Valid for 0<rho<R. At fixed w this is injective on the atlas phase
    quotient; varying w is projected out and must be separately encoded.
    This uses theta1=alpha-gamma, theta2=-alpha-gamma. It is intentionally
    distinct from the older representative-dependent (alpha,gamma) torus.
    """
    theta1, theta2, rho = np.broadcast_arrays(theta1, theta2, rho)
    if np.any((rho <= 0) | (rho >= major_radius)):
        raise ValueError("The torus chart requires 0<rho<major_radius.")
    radius = major_radius + rho*np.cos(theta1)
    return np.stack((radius*np.cos(theta2), radius*np.sin(theta2), rho*np.sin(theta1)), axis=-1)


def project_q(q, family: str, *, major_radius: float = 4.0, validate: bool = True):
    """Return (...,3) positions together with complete atlas metadata."""
    atlas = atlas_coordinates(q, family, validate=validate)
    return coefficient_torus(atlas["theta1"], atlas["theta2"], atlas["rho"], major_radius=major_radius), atlas


def inverse_coefficient_torus(xyz, *, major_radius: float = 4.0):
    """Recover theta1,theta2,rho; weight is not present in the 3D projection."""
    xyz = np.asarray(xyz, dtype=float)
    X, Y, Z = np.moveaxis(xyz, -1, 0)
    u = np.hypot(X, Y)-major_radius
    return {"theta1": np.arctan2(Z, u), "theta2": np.arctan2(Y, X), "rho": np.hypot(u, Z)}


def from_log_ratio(xi):
    """q=1/(1+exp(xi)); necessary stable strip |Im(xi)|<pi/2."""
    return 1/(1+np.exp(np.asarray(xi, dtype=np.complex128)))


def stable_radial_bounds(w):
    """Strict contraction and necessary stable band in atlas coordinates.

    Stable nondegenerate binary zippers satisfy square_sum < rho <= 2.
    The inequalities are necessary; they do not certify a stable point.
    """
    w = np.asarray(w, dtype=float)
    return {"contraction": 2*np.maximum(w,1-w),
            "square_sum": 2*np.sqrt(w*w+(1-w)*(1-w)), "outer": np.full_like(w, 2)}
