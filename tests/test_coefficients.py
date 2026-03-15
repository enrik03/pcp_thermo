"""
tests/test_coefficients.py
==========================
Unit tests for the PCP thermodynamic coefficients and figure generation.

Tests verify:
  - Sign and order-of-magnitude of C4 and C7.
  - T_c_Phi is in the expected range.
  - T_x is finite and positive (when C7 != 0).
  - All three figure output files are created when the main script runs.
"""

from __future__ import annotations

import importlib
import math
import os
import sys
import tempfile

import matplotlib
import pytest

# Use non-interactive backend for all tests
matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Make the package importable regardless of working directory
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import pcp_thermo_figures_corrected as pcp  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def coeffs():
    """Compute PCP coefficients once for all tests in this module."""
    return pcp.compute_pcp_coefficients(M=1.0)


# ---------------------------------------------------------------------------
# Tests: individual coefficients
# ---------------------------------------------------------------------------

class TestC4:
    """Tests for the C4 coefficient (leading T^4 term)."""

    def test_c4_is_positive(self, coeffs):
        """C4 must be positive (pressure is positive)."""
        assert coeffs["C4"] > 0.0, f"Expected C4 > 0, got {coeffs['C4']}"

    def test_c4_order_of_magnitude(self, coeffs):
        """C4 should be O(0.1) for the chosen degeneracy factors."""
        assert 0.05 < coeffs["C4"] < 5.0, (
            f"C4 = {coeffs['C4']} is outside the expected range (0.05, 5.0)"
        )

    def test_c4_less_than_c_ur(self, coeffs):
        """PCP C4 should be less than the full UR coefficient C_UR."""
        assert coeffs["C4"] < pcp.C_UR, (
            f"Expected C4 ({coeffs['C4']:.4g}) < C_UR ({pcp.C_UR:.4g})"
        )


class TestC7:
    """Tests for the C7 coefficient (sub-leading T^7/M^3 term)."""

    def test_c7_is_positive(self, coeffs):
        """C7 should be positive for the representative PCP coupling."""
        assert coeffs["C7"] > 0.0, f"Expected C7 > 0, got {coeffs['C7']}"

    def test_c7_order_of_magnitude(self, coeffs):
        """C7 should be a small positive number (PCP coupling suppression)."""
        assert 1e-5 < coeffs["C7"] < 1.0, (
            f"C7 = {coeffs['C7']:.4g} is outside the expected range (1e-5, 1.0)"
        )

    def test_c7_much_less_than_c4(self, coeffs):
        """C7 must be much smaller than C4 (perturbative regime)."""
        assert coeffs["C7"] < coeffs["C4"] * 1e-1, (
            f"C7 ({coeffs['C7']:.4g}) is not <<  C4 ({coeffs['C4']:.4g})"
        )


class TestTcPhi:
    """Tests for the PCP condensate temperature T_c^{(Phi)}."""

    def test_tc_phi_is_positive(self, coeffs):
        assert coeffs["T_c_Phi"] > 0.0

    def test_tc_phi_range(self, coeffs):
        """With g=0.25 and M=1, T_c_Phi = sqrt(48) ≈ 6.928."""
        expected = math.sqrt(48.0)
        assert abs(coeffs["T_c_Phi"] - expected) < 1e-6, (
            f"T_c_Phi = {coeffs['T_c_Phi']:.6f}, expected {expected:.6f}"
        )


class TestTx:
    """Tests for the crossover temperature T_x."""

    def test_tx_is_finite(self, coeffs):
        """T_x must be a finite float (not NaN) when C7 > 0."""
        assert math.isfinite(coeffs["T_x"]), f"T_x is not finite: {coeffs['T_x']}"

    def test_tx_is_positive(self, coeffs):
        assert coeffs["T_x"] > 0.0, f"Expected T_x > 0, got {coeffs['T_x']}"

    def test_tx_crossover_condition(self, coeffs):
        """Verify that P_PCP(T_x) ≈ P_UR(T_x) (within 5%)."""
        import numpy as np

        t = coeffs["T_x"]
        p_pcp = pcp.pcp_pressure(
            np.array([t]), coeffs["C4"], coeffs["C7"], M=1.0
        )[0]
        p_ur = pcp.ur_pressure(np.array([t]))[0]
        relative_diff = abs(p_pcp - p_ur) / max(abs(p_ur), 1e-30)
        assert relative_diff < 0.05, (
            f"Crossover condition not satisfied: "
            f"P_PCP={p_pcp:.4g}, P_UR={p_ur:.4g}, rel_diff={relative_diff:.4g}"
        )


# ---------------------------------------------------------------------------
# Tests: figure generation
# ---------------------------------------------------------------------------

EXPECTED_STEMS = [
    "pcp_pressure_comparison_main",
    "pcp_pressure_crossover_appendix",
    "pcp_pressure_decomposition_appendix",
]


@pytest.mark.parametrize("fmt", ["png", "pdf"])
def test_figures_created(fmt, tmp_path):
    """Running main() must create all three figure files in the chosen format."""
    outdir = str(tmp_path / "figures")
    exit_code = pcp.main(["--outdir", outdir, "--format", fmt, "--dpi", "72"])
    assert exit_code == 0, f"main() returned exit code {exit_code}"

    for stem in EXPECTED_STEMS:
        expected_path = os.path.join(outdir, f"{stem}.{fmt}")
        assert os.path.isfile(expected_path), (
            f"Expected figure not found: {expected_path}"
        )
        assert os.path.getsize(expected_path) > 0, (
            f"Figure file is empty: {expected_path}"
        )


def test_figures_created_both_formats(tmp_path):
    """Running with --format png,pdf must create files in both formats."""
    outdir = str(tmp_path / "figures")
    exit_code = pcp.main(
        ["--outdir", outdir, "--format", "png,pdf", "--dpi", "72"]
    )
    assert exit_code == 0

    for stem in EXPECTED_STEMS:
        for fmt in ("png", "pdf"):
            path = os.path.join(outdir, f"{stem}.{fmt}")
            assert os.path.isfile(path), f"Missing: {path}"


def test_invalid_format_returns_error(tmp_path):
    """An empty format string should cause main() to return exit code 1."""
    outdir = str(tmp_path / "figures")
    exit_code = pcp.main(["--outdir", outdir, "--format", ""])
    assert exit_code == 1


# ---------------------------------------------------------------------------
# Tests: pressure functions
# ---------------------------------------------------------------------------

class TestPressureFunctions:
    """Sanity checks for pcp_pressure() and ur_pressure()."""

    def test_ur_pressure_positive(self):
        import numpy as np

        T = np.linspace(0.1, 5.0, 20)
        P = pcp.ur_pressure(T)
        assert np.all(P > 0)

    def test_ur_pressure_scales_as_t4(self):
        import numpy as np

        T1, T2 = 1.0, 2.0
        ratio = pcp.ur_pressure(np.array([T2]))[0] / pcp.ur_pressure(np.array([T1]))[0]
        assert abs(ratio - 16.0) < 1e-10, f"Expected ratio 16, got {ratio}"

    def test_pcp_pressure_order1_equals_c4_t4(self, coeffs):
        import numpy as np

        T = np.array([1.0, 2.0, 3.0])
        P = pcp.pcp_pressure(T, coeffs["C4"], coeffs["C7"], order=1)
        expected = coeffs["C4"] * T ** 4
        np.testing.assert_allclose(P, expected)

    def test_pcp_pressure_order2_includes_t7(self, coeffs):
        import numpy as np

        T = np.array([1.0])
        P1 = pcp.pcp_pressure(T, coeffs["C4"], coeffs["C7"], order=1)[0]
        P2 = pcp.pcp_pressure(T, coeffs["C4"], coeffs["C7"], order=2)[0]
        assert P2 != P1, "order=2 should differ from order=1 when C7 != 0"
