#!/usr/bin/env python3
"""
pcp_thermo_figures_corrected.py
================================
Main script for the PCP (post-Carroll-Pauli) thermodynamic pressure analysis.

Computes the numerical coefficients, evaluates thermodynamic integrals, and
generates three publication-quality figures:

  1. pcp_pressure_comparison_main          – PCP vs UR pressure scaling (main text)
  2. pcp_pressure_crossover_appendix       – Crossover temperature T_x analysis
  3. pcp_pressure_decomposition_appendix   – Pressure decomposition into PCP modes

Usage
-----
    python pcp_thermo_figures_corrected.py [--outdir DIR] [--dpi DPI] [--format FMT]

Options
-------
    --outdir DIR   Output directory for figures (default: figures/)
    --dpi DPI      Resolution in dots per inch (default: 150)
    --format FMT   Comma-separated list of formats, e.g. png,pdf (default: png,pdf)

References
----------
    Post-Carroll-Pauli thermodynamics; see pcp_thermo_package_notes.md for
    detailed derivations and LaTeX fragments.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Sequence

import matplotlib
import numpy as np
from scipy import integrate

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Physical / model constants
# ---------------------------------------------------------------------------

# Degeneracy factors and group-theory coefficients
N_C: float = 3.0          # number of colours
N_F: float = 2.0          # number of light flavours
_PI2: float = np.pi ** 2

# Ideal-gas Stefan-Boltzmann coefficients (in units where k_B = hbar = c = 1)
# P_UR = C_UR * T^4,   C_UR = (pi^2/90) * g_eff
G_EFF_UR: float = (
    2.0 * (N_C ** 2 - 1.0)          # gluons (bosonic)
    + 7.0 / 8.0 * 2.0 * 2.0 * N_C * N_F  # quarks (fermionic, spin + colour + flavour)
)
C_UR: float = (_PI2 / 90.0) * G_EFF_UR

# PCP sector coefficients (leading terms of the low-temperature expansion)
# P_PCP = C4 * T^4 + C7 * T^7 / M^3 + ...
# where M is a heavy-field mass scale (set to 1 in dimensionless units).
# These module-level values are placeholders; they are overwritten by
# compute_pcp_coefficients() when main() is called.
C4: float = 0.0
C7: float = 0.0

# Crossover temperature (dimensionless, M=1)
T_C_PHI: float = 0.0   # characteristic temperature for Phi condensate, set below
T_X: float = 0.0       # crossover temperature P_PCP ~ P_UR, set below


def compute_pcp_coefficients(M: float = 1.0) -> dict[str, float]:
    """Compute the key PCP thermodynamic coefficients numerically.

    The PCP pressure receives contributions from a heavy Phi sector and the
    standard quarks/gluons screened by PCP corrections.  We parameterise:

        P_PCP(T) = C4 * T^4  +  C7 * T^7 / M^3  +  O(T^10)

    where C4 and C7 are computed from one-loop thermal integrals.

    Parameters
    ----------
    M : float
        Heavy PCP field mass scale (default 1, i.e. dimensionless units).

    Returns
    -------
    dict with keys ``C4``, ``C7``, ``T_c_Phi``, ``T_x``.
    """

    # ------------------------------------------------------------------
    # C4: coefficient of T^4 in the PCP pressure
    # Comes from integrating the one-loop thermal integral for massless modes
    # and their PCP-corrected counterparts.
    # For massless bosons: J_B(0) = -pi^4/45 => P = T^4 * (pi^2/90)
    # Here we include the modified dispersion relation from PCP.
    # ------------------------------------------------------------------
    def integrand_c4(u: float) -> float:
        """Dimensionless integrand for C4 (massless PCP limit)."""
        if u <= 0.0:
            return 0.0
        return u ** 2 * np.log(1.0 - np.exp(-u))

    val_c4, _ = integrate.quad(integrand_c4, 0.0, np.inf, limit=200)
    # P_bosonic / T^4 = -1/(2*pi^2) * integral  (factor of -1 for log(1-e^{-u}))
    # We include the gluon and quark degeneracies
    g_boson = 2.0 * (N_C ** 2 - 1.0)
    g_fermion = 2.0 * N_C * N_F  # factor 7/8 absorbed below

    c4_boson = -g_boson / (2.0 * _PI2) * val_c4
    # Fermion integral: J_F(0) = +7*pi^4/720
    def integrand_c4_f(u: float) -> float:
        if u <= 0.0:
            return 0.0
        return u ** 2 * np.log(1.0 + np.exp(-u))

    val_c4_f, _ = integrate.quad(integrand_c4_f, 0.0, np.inf, limit=200)
    c4_fermion = g_fermion / (2.0 * _PI2) * val_c4_f

    c4 = c4_boson + c4_fermion

    # ------------------------------------------------------------------
    # C7: coefficient of T^7/M^3 — leading PCP correction
    # Arises from expanding the heavy-field propagator at T >> M.
    # Schematically: C7 = g_PCP / (M^3) * integral(u^5 * |n_B'(u)|)
    # We use a representative PCP coupling: g_PCP = 1/(16*pi^2).
    # ------------------------------------------------------------------
    g_pcp = 1.0 / (16.0 * _PI2)

    def integrand_c7(u: float) -> float:
        """Dimensionless integrand for C7 correction.

        Uses the identity  exp(u)/(exp(u)-1)^2 = 1/(2*cosh(u)-2)
        and falls back to the asymptotic form -u^5*exp(-u) for u > 20
        to avoid floating-point overflow.
        """
        if u <= 1e-10:
            return 0.0
        if u > 20.0:
            # Asymptotic: 1/(2*cosh(u)-2) ~ exp(-u)
            return -u ** 5 * np.exp(-u)
        return u ** 5 * (-1.0) / (2.0 * np.cosh(u) - 2.0)

    val_c7, _ = integrate.quad(integrand_c7, 1e-10, np.inf, limit=300)
    c7 = g_pcp / (2.0 * _PI2) * (-val_c7) / M ** 3

    # ------------------------------------------------------------------
    # T_c^{Phi}: characteristic condensation temperature of the Phi field
    # Defined as the temperature where the Phi thermal mass m_th(T) ~ M:
    #   m_th^2(T) = M^2 - g * T^2/12  =>  T_c_Phi = sqrt(12) * M / sqrt(g)
    # with coupling g = 1/4 (representative PCP value).
    # ------------------------------------------------------------------
    g_coupling = 0.25
    t_c_phi = np.sqrt(12.0 / g_coupling) * M

    # ------------------------------------------------------------------
    # T_x: crossover temperature where P_PCP(T) ~ P_UR(T)
    # Solve:  C4 * T_x^4 + C7 * T_x^7 / M^3  =  C_UR * T_x^4
    # =>      C7 * T_x^3 / M^3  =  (C_UR - C4)
    # =>      T_x = M * ((C_UR - C4) / C7)^(1/3)   [if C7 != 0]
    # ------------------------------------------------------------------
    delta_c = C_UR - c4
    if abs(c7) > 1e-30:
        t_x = M * (delta_c / c7) ** (1.0 / 3.0) if delta_c / c7 > 0 else np.nan
    else:
        t_x = np.nan

    return {
        "C4": c4,
        "C7": c7,
        "T_c_Phi": t_c_phi,
        "T_x": t_x,
    }


def pcp_pressure(
    T: np.ndarray,
    c4: float,
    c7: float,
    M: float = 1.0,
    order: int = 2,
) -> np.ndarray:
    """Compute the PCP pressure as a polynomial in T.

    Parameters
    ----------
    T : array_like
        Temperature values (dimensionless, M=1 units).
    c4 : float
        Coefficient of T^4.
    c7 : float
        Coefficient of T^7 / M^3.
    M : float
        Heavy-field mass scale.
    order : int
        Number of terms to include (1 => T^4 only, 2 => T^4 + T^7).

    Returns
    -------
    np.ndarray
        PCP pressure array.
    """
    T = np.asarray(T, dtype=float)
    P = c4 * T ** 4
    if order >= 2:
        P = P + c7 * T ** 7 / M ** 3
    return P


def ur_pressure(T: np.ndarray) -> np.ndarray:
    """Ultra-relativistic (Stefan-Boltzmann) pressure P_UR = C_UR * T^4.

    Parameters
    ----------
    T : array_like
        Temperature array.

    Returns
    -------
    np.ndarray
    """
    return C_UR * np.asarray(T, dtype=float) ** 4


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------

def make_figure_comparison(
    T: np.ndarray,
    P_pcp: np.ndarray,
    P_ur: np.ndarray,
    t_x: float,
    outdir: str,
    dpi: int,
    formats: list[str],
) -> None:
    """Figure 1 – PCP vs UR pressure comparison (main text figure).

    Shows normalised pressures P/T^4 vs T on a log scale to highlight the
    PCP corrections and the crossover at T_x.

    Parameters
    ----------
    T : np.ndarray
        Temperature array.
    P_pcp, P_ur : np.ndarray
        Pressure arrays.
    t_x : float
        Crossover temperature (drawn as vertical line).
    outdir : str
        Output directory.
    dpi : int
        Figure resolution.
    formats : list[str]
        File formats to save.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))

    mask = T > 0
    ax.plot(T[mask], P_pcp[mask] / T[mask] ** 4, color="steelblue",
            lw=2, label=r"$P_{\rm PCP}(T)/T^4$")
    ax.plot(T[mask], P_ur[mask] / T[mask] ** 4, color="tomato",
            lw=2, linestyle="--", label=r"$P_{\rm UR}(T)/T^4$")

    if np.isfinite(t_x) and T[0] < t_x < T[-1]:
        ax.axvline(t_x, color="gray", lw=1.2, linestyle=":", label=rf"$T_x = {t_x:.3f}\,M$")

    ax.set_xlabel(r"$T / M$", fontsize=13)
    ax.set_ylabel(r"$P(T) / T^4$", fontsize=13)
    ax.set_title("PCP vs UR Pressure: Normalised Comparison", fontsize=13)
    ax.set_xscale("log")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    _save_figure(fig, "pcp_pressure_comparison_main", outdir, dpi, formats)
    plt.close(fig)


def make_figure_crossover(
    T: np.ndarray,
    P_pcp: np.ndarray,
    P_ur: np.ndarray,
    t_x: float,
    t_c_phi: float,
    outdir: str,
    dpi: int,
    formats: list[str],
) -> None:
    """Figure 2 – Crossover analysis (appendix figure).

    Plots P_PCP / P_UR vs T/M to show where the two descriptions coincide.

    Parameters
    ----------
    T : np.ndarray
        Temperature array.
    P_pcp, P_ur : np.ndarray
        Pressure arrays.
    t_x : float
        Crossover temperature.
    t_c_phi : float
        PCP condensate temperature.
    outdir : str
        Output directory.
    dpi : int
        Figure resolution.
    formats : list[str]
        File formats.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))

    mask = (T > 0) & (P_ur > 0)
    ratio = P_pcp[mask] / P_ur[mask]
    ax.plot(T[mask], ratio, color="mediumseagreen", lw=2)
    ax.axhline(1.0, color="black", lw=1.0, linestyle="--",
               label=r"$P_{\rm PCP}/P_{\rm UR}=1$")

    if np.isfinite(t_x) and T[0] < t_x < T[-1]:
        ax.axvline(t_x, color="gray", lw=1.2, linestyle=":",
                   label=rf"$T_x = {t_x:.3f}\,M$")

    if np.isfinite(t_c_phi) and T[0] < t_c_phi < T[-1]:
        ax.axvline(t_c_phi, color="orchid", lw=1.2, linestyle="-.",
                   label=rf"$T_c^{{(\Phi)}} = {t_c_phi:.3f}\,M$")

    ax.set_xlabel(r"$T / M$", fontsize=13)
    ax.set_ylabel(r"$P_{\rm PCP}(T) / P_{\rm UR}(T)$", fontsize=13)
    ax.set_title("PCP / UR Pressure Ratio and Crossover", fontsize=13)
    ax.set_xscale("log")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    _save_figure(fig, "pcp_pressure_crossover_appendix", outdir, dpi, formats)
    plt.close(fig)


def make_figure_decomposition(
    T: np.ndarray,
    c4: float,
    c7: float,
    M: float,
    outdir: str,
    dpi: int,
    formats: list[str],
) -> None:
    """Figure 3 – PCP pressure decomposition (appendix figure).

    Shows the individual contributions P ~ C4*T^4 and P ~ C7*T^7/M^3 to the
    total PCP pressure, illustrating the dominant terms at different scales.

    Parameters
    ----------
    T : np.ndarray
        Temperature array.
    c4, c7 : float
        Pressure coefficients.
    M : float
        Heavy-field mass scale.
    outdir : str
        Output directory.
    dpi : int
        Resolution.
    formats : list[str]
        File formats.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))

    T_pos = T[T > 0]
    P_t4 = c4 * T_pos ** 4
    P_t7 = c7 * T_pos ** 7 / M ** 3
    P_total = P_t4 + P_t7

    ax.plot(T_pos, np.abs(P_t4), color="steelblue", lw=2,
            label=r"$|C_4|\,T^4$")
    ax.plot(T_pos, np.abs(P_t7), color="darkorange", lw=2, linestyle="--",
            label=r"$|C_7|\,T^7/M^3$")
    ax.plot(T_pos, np.abs(P_total), color="black", lw=1.5, linestyle="-.",
            label=r"$|P_{\rm PCP}|$  (total)")

    ax.set_xlabel(r"$T / M$", fontsize=13)
    ax.set_ylabel(r"$|P|$ (dimensionless)", fontsize=13)
    ax.set_title("PCP Pressure Decomposition by Power", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    _save_figure(fig, "pcp_pressure_decomposition_appendix", outdir, dpi, formats)
    plt.close(fig)


def _save_figure(
    fig: "matplotlib.figure.Figure",
    stem: str,
    outdir: str,
    dpi: int,
    formats: list[str],
) -> None:
    """Save *fig* to *outdir*/<stem>.<fmt> for each format in *formats*.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    stem : str
        Base file name (without extension).
    outdir : str
        Output directory (created if absent).
    dpi : int
        Resolution for raster formats.
    formats : list[str]
        File format extensions, e.g. ``["png", "pdf"]``.
    """
    os.makedirs(outdir, exist_ok=True)
    for fmt in formats:
        path = os.path.join(outdir, f"{stem}.{fmt}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        log.info("Saved %s", path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate PCP thermodynamic pressure figures "
            "(pcp_pressure_comparison_main, crossover_appendix, "
            "decomposition_appendix)."
        )
    )
    parser.add_argument(
        "--outdir",
        default="figures",
        help="Output directory for figures (default: figures/)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Figure resolution in DPI (default: 150)",
    )
    parser.add_argument(
        "--format",
        default="png,pdf",
        dest="fmt",
        help="Comma-separated output formats, e.g. png,pdf (default: png,pdf)",
    )
    parser.add_argument(
        "--T-min",
        type=float,
        default=0.01,
        dest="T_min",
        help="Minimum temperature T/M (default: 0.01)",
    )
    parser.add_argument(
        "--T-max",
        type=float,
        default=10.0,
        dest="T_max",
        help="Maximum temperature T/M (default: 10.0)",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=500,
        dest="n_points",
        help="Number of temperature grid points (default: 500)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Main entry point.

    Parameters
    ----------
    argv : list[str], optional
        Command-line arguments (defaults to sys.argv).

    Returns
    -------
    int
        Exit code (0 = success).
    """
    args = parse_args(argv)
    formats = [f.strip().lower() for f in args.fmt.split(",") if f.strip()]
    if not formats:
        log.error("No output formats specified.")
        return 1

    M = 1.0  # dimensionless units

    # ------------------------------------------------------------------
    # Compute coefficients
    # ------------------------------------------------------------------
    log.info("Computing PCP thermodynamic coefficients (M = %.4g)…", M)
    coeffs = compute_pcp_coefficients(M=M)
    c4_val = coeffs["C4"]
    c7_val = coeffs["C7"]
    t_c_phi = coeffs["T_c_Phi"]
    t_x = coeffs["T_x"]

    # Update module-level globals so tests can access them
    global C4, C7, T_C_PHI, T_X
    C4 = c4_val
    C7 = c7_val
    T_C_PHI = t_c_phi
    T_X = t_x

    log.info("  C4        = %.6g", c4_val)
    log.info("  C7        = %.6g", c7_val)
    log.info("  T_c_Phi   = %.6g M", t_c_phi)
    log.info("  T_x       = %.6g M", t_x)

    # ------------------------------------------------------------------
    # Temperature grid
    # ------------------------------------------------------------------
    T = np.linspace(args.T_min, args.T_max, args.n_points)

    # ------------------------------------------------------------------
    # Pressure arrays
    # ------------------------------------------------------------------
    P_pcp = pcp_pressure(T, c4_val, c7_val, M=M)
    P_ur = ur_pressure(T)

    # ------------------------------------------------------------------
    # Generate figures
    # ------------------------------------------------------------------
    log.info("Generating figures → %s", args.outdir)

    make_figure_comparison(T, P_pcp, P_ur, t_x, args.outdir, args.dpi, formats)
    make_figure_crossover(T, P_pcp, P_ur, t_x, t_c_phi, args.outdir, args.dpi, formats)
    make_figure_decomposition(T, c4_val, c7_val, M, args.outdir, args.dpi, formats)

    log.info("Done.  %d figure(s) × %d format(s) written to '%s'.",
             3, len(formats), args.outdir)
    return 0


if __name__ == "__main__":
    # Use non-interactive backend when running as a script (safe for CI)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401 – imported after backend set
    sys.exit(main())
