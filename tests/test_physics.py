"""Unit tests for VG functions, TDMA solver, and Laio ET stress function.

Physics correctness is verified against analytical solutions and known limits.
These tests run in pure Python — no Fortran kernel is required.

The Python reference implementations below mirror the elemental Fortran
functions in physics.f08 and the TDMA routine in solver_implicit.f08.  Any
discrepancy with the Fortran would surface as a physics bug in both places.
"""

import math

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Python reference implementations
# ---------------------------------------------------------------------------


def _vg_Se(h: float, alpha: float, n: float) -> float:
    """Van Genuchten effective saturation Se(h)."""
    if h >= 0.0:
        return 1.0
    m = 1.0 - 1.0 / n
    return 1.0 / (1.0 + (alpha * abs(h)) ** n) ** m


def _vg_theta(h: float, theta_r: float, theta_s: float, alpha: float, n: float) -> float:
    """Van Genuchten volumetric water content theta(h)."""
    return theta_r + (theta_s - theta_r) * _vg_Se(h, alpha, n)


def _vg_K(h: float, K_sat: float, alpha: float, n: float, lam: float = 0.5) -> float:
    """Mualem–van Genuchten unsaturated conductivity K(h)."""
    if h >= 0.0:
        return K_sat
    m = 1.0 - 1.0 / n
    Se = _vg_Se(h, alpha, n)
    Se_c = max(Se, 1e-12)
    return K_sat * Se_c**lam * (1.0 - (1.0 - Se_c ** (1.0 / m)) ** m) ** 2


def _vg_C(h: float, theta_r: float, theta_s: float, alpha: float, n: float) -> float:
    """Van Genuchten specific moisture capacity C(h) = dtheta/dh.

    At h >= 0 (saturated / at water table), returns Sy = theta_s - theta_r,
    matching the GWSWEX convention for unconfined phreatic storativity.
    """
    Sy = theta_s - theta_r
    if h >= 0.0:
        return Sy
    m = 1.0 - 1.0 / n
    x = (alpha * abs(h)) ** n
    return Sy * m * n * alpha**n * abs(h) ** (n - 1) / (1.0 + x) ** (m + 1)


def _tdma(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Thomas algorithm for the tridiagonal system  a*x_{i-1} + b*x_i + c*x_{i+1} = d.

    a[0] and c[-1] are ignored (boundary sub/super-diagonals).
    """
    n = len(d)
    cp = np.zeros(n)
    dp = np.zeros(n)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        denom = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / denom if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom
    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def _laio_stress_T(s: float, s_star: float, s_w: float) -> float:
    """Laio (2001) piecewise-linear transpiration stress factor in [0, 1]."""
    if s >= s_star:
        return 1.0
    if s > s_w:
        return (s - s_w) / (s_star - s_w)
    return 0.0


def _laio_stress_E(s: float, s_h: float, s_e: float) -> float:
    """Laio (2001) piecewise-linear evaporation stress factor in [0, 1]."""
    if s >= s_e:
        return 1.0
    if s > s_h:
        return (s - s_h) / (s_e - s_h)
    return 0.0


# ---------------------------------------------------------------------------
# Shared soil parameters (van Genuchten topsoil values, matching test fixtures)
# ---------------------------------------------------------------------------
ALPHA = 3.6
VG_N = 1.56
THETA_R = 0.078
THETA_S = 0.43
K_SAT = 1e-5
LAM = 0.5
VG_M = 1.0 - 1.0 / VG_N
SY = THETA_S - THETA_R


# ---------------------------------------------------------------------------
class TestVanGenuchten:
    """Van Genuchten (1980) retention and conductivity function properties."""

    # ---- theta(h) ----------------------------------------------------------

    def test_theta_saturated_at_h_zero(self):
        assert math.isclose(_vg_theta(0.0, THETA_R, THETA_S, ALPHA, VG_N), THETA_S)

    def test_theta_stays_theta_s_for_positive_h(self):
        # h > 0: below water table — Se = 1, theta = theta_s
        assert math.isclose(_vg_theta(0.5, THETA_R, THETA_S, ALPHA, VG_N), THETA_S)

    def test_theta_approaches_theta_r_when_very_dry(self):
        # VG never reaches exactly theta_r; at h=-1000 m (~100 bar), the difference
        # is ~0.003 for these silty loam parameters — physically correct.
        theta = _vg_theta(-1000.0, THETA_R, THETA_S, ALPHA, VG_N)
        assert abs(theta - THETA_R) < 0.01

    def test_theta_monotone_decreasing_with_drying(self):
        thetas = [_vg_theta(h, THETA_R, THETA_S, ALPHA, VG_N) for h in [0.0, -0.1, -1.0, -10.0]]
        assert all(thetas[i] >= thetas[i + 1] for i in range(len(thetas) - 1))

    def test_theta_between_theta_r_and_theta_s(self):
        for h in [-0.01, -0.1, -1.0, -10.0]:
            theta = _vg_theta(h, THETA_R, THETA_S, ALPHA, VG_N)
            assert THETA_R <= theta <= THETA_S, f"h={h}: theta={theta:.4f} out of [theta_r, theta_s]"

    # ---- K(h) --------------------------------------------------------------

    def test_K_equals_Ksat_at_h_zero(self):
        assert math.isclose(_vg_K(0.0, K_SAT, ALPHA, VG_N, LAM), K_SAT)

    def test_K_equals_Ksat_for_positive_h(self):
        assert math.isclose(_vg_K(1.0, K_SAT, ALPHA, VG_N, LAM), K_SAT)

    def test_K_strictly_less_than_Ksat_when_unsaturated(self):
        K = _vg_K(-0.01, K_SAT, ALPHA, VG_N, LAM)
        assert K < K_SAT

    def test_K_monotone_decreasing_with_drying(self):
        heads = [0.0, -0.1, -1.0, -10.0, -100.0]
        Ks = [_vg_K(h, K_SAT, ALPHA, VG_N, LAM) for h in heads]
        assert all(Ks[i] >= Ks[i + 1] for i in range(len(Ks) - 1))

    def test_K_near_zero_when_very_dry(self):
        K = _vg_K(-1000.0, K_SAT, ALPHA, VG_N, LAM)
        assert K < K_SAT * 1e-6

    def test_K_positive_everywhere(self):
        for h in [-0.01, -0.1, -1.0, -10.0, -100.0]:
            assert _vg_K(h, K_SAT, ALPHA, VG_N, LAM) > 0.0

    # ---- C(h) = dtheta/dh --------------------------------------------------

    def test_C_at_saturation_equals_Sy(self):
        """C(h=0) = Sy — the drainable porosity (GWSWEX phreatic convention)."""
        assert math.isclose(_vg_C(0.0, THETA_R, THETA_S, ALPHA, VG_N), SY, rel_tol=1e-9)

    def test_C_positive_for_all_unsaturated_h(self):
        for h in [-0.01, -0.1, -1.0, -10.0]:
            assert _vg_C(h, THETA_R, THETA_S, ALPHA, VG_N) > 0.0, f"C(h={h}) not positive"

    def test_C_consistent_with_finite_difference_of_theta(self):
        """C(h) ≈ Δtheta / Δh via a central finite difference."""
        h0 = -1.0
        dh = 1e-5
        C_ana = _vg_C(h0, THETA_R, THETA_S, ALPHA, VG_N)
        C_fd = (
            _vg_theta(h0 + dh, THETA_R, THETA_S, ALPHA, VG_N) - _vg_theta(h0 - dh, THETA_R, THETA_S, ALPHA, VG_N)
        ) / (2 * dh)
        assert abs(C_ana - C_fd) / max(abs(C_ana), 1e-12) < 1e-4

    def test_C_decreases_toward_zero_when_very_dry(self):
        C_moderate = _vg_C(-1.0, THETA_R, THETA_S, ALPHA, VG_N)
        C_very_dry = _vg_C(-1000.0, THETA_R, THETA_S, ALPHA, VG_N)
        assert C_very_dry < C_moderate


# ---------------------------------------------------------------------------
class TestTDMA:
    """Thomas algorithm (TDMA) for solving tridiagonal linear systems."""

    def test_diagonal_only_system(self):
        """With zero off-diagonals the solution is element-wise d/b."""
        n = 5
        b = np.full(n, 3.0)
        a = np.zeros(n)
        c = np.zeros(n)
        d = np.arange(1.0, n + 1.0)
        x = _tdma(a, b, c, d)
        np.testing.assert_allclose(x, d / b, atol=1e-14)

    def test_3x3_known_system(self):
        """3×3 tridiagonal: [2,-1; -1,2,-1; -1,2] * [1,1,1]^T = [1,0,1]^T."""
        a = np.array([0.0, -1.0, -1.0])
        b = np.array([2.0, 2.0, 2.0])
        c = np.array([-1.0, -1.0, 0.0])
        d = np.array([1.0, 0.0, 1.0])
        np.testing.assert_allclose(_tdma(a, b, c, d), [1.0, 1.0, 1.0], atol=1e-14)

    def test_matches_numpy_linalg_solve(self):
        """General 6×6 diagonally dominant system: TDMA == numpy.linalg.solve."""
        n = 6
        rng = np.random.default_rng(0)
        off = rng.uniform(-1.0, 1.0, n - 1)
        b = (np.abs(off).sum() + 1.0) * np.ones(n)  # strict diagonal dominance
        a = np.concatenate([[0.0], off])
        c = np.concatenate([off, [0.0]])
        rhs = rng.uniform(1.0, 5.0, n)
        A = np.diag(b) + np.diag(off, -1) + np.diag(off, 1)
        x_ref = np.linalg.solve(A, rhs)
        np.testing.assert_allclose(_tdma(a, b, c, rhs), x_ref, atol=1e-12)

    def test_identity_system(self):
        """Identity diagonal (b=1, a=c=0) gives x = d."""
        n = 8
        d = np.linspace(-3.0, 3.0, n)
        np.testing.assert_allclose(_tdma(np.zeros(n), np.ones(n), np.zeros(n), d), d, atol=1e-14)

    def test_larger_random_system(self):
        """10×10 random diagonally dominant system."""
        n = 10
        rng = np.random.default_rng(1)
        off_l = rng.uniform(-2.0, 0.0, n - 1)
        off_u = rng.uniform(-2.0, 0.0, n - 1)
        b = (np.abs(off_l).sum() + np.abs(off_u).sum() + 1.0) * np.ones(n)
        a = np.concatenate([[0.0], off_l])
        c = np.concatenate([off_u, [0.0]])
        rhs = rng.standard_normal(n)
        A = np.diag(b) + np.diag(off_l, -1) + np.diag(off_u, 1)
        x_ref = np.linalg.solve(A, rhs)
        np.testing.assert_allclose(_tdma(a, b, c, rhs), x_ref, atol=1e-11)


# ---------------------------------------------------------------------------
class TestLaioStress:
    """Laio (2001) piecewise-linear ET stress function.

    Reference: Laio, F., Porporato, A., Ridolfi, L., & Rodriguez-Iturbe, I. (2001).
    Plants in water-controlled ecosystems: Active role in hydrologic processes and
    response to water stress. Adv. Water Resour., 24(7), 707–723.
    """

    # ---- Transpiration stress -----------------------------------------------

    def test_T_full_at_s_star(self):
        assert _laio_stress_T(0.5, s_star=0.5, s_w=0.1) == 1.0

    def test_T_full_above_s_star(self):
        assert _laio_stress_T(0.8, s_star=0.5, s_w=0.1) == 1.0

    def test_T_zero_at_wilting_point(self):
        assert _laio_stress_T(0.1, s_star=0.5, s_w=0.1) == 0.0

    def test_T_zero_below_wilting_point(self):
        assert _laio_stress_T(0.05, s_star=0.5, s_w=0.1) == 0.0

    def test_T_midpoint_gives_half_stress(self):
        s_w, s_star = 0.1, 0.5
        s_mid = 0.5 * (s_w + s_star)
        assert math.isclose(_laio_stress_T(s_mid, s_star, s_w), 0.5, abs_tol=1e-14)

    def test_T_quarter_point(self):
        s_w, s_star = 0.1, 0.5
        s = s_w + 0.25 * (s_star - s_w)
        assert math.isclose(_laio_stress_T(s, s_star, s_w), 0.25, abs_tol=1e-14)

    def test_T_monotone_increasing_with_saturation(self):
        s_w, s_star = 0.1, 0.6
        saturations = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8]
        stresses = [_laio_stress_T(s, s_star, s_w) for s in saturations]
        assert all(stresses[i] <= stresses[i + 1] for i in range(len(stresses) - 1))

    # ---- Evaporation stress -------------------------------------------------

    def test_E_full_at_s_e(self):
        assert _laio_stress_E(0.5, s_h=0.05, s_e=0.5) == 1.0

    def test_E_full_above_s_e(self):
        assert _laio_stress_E(0.7, s_h=0.05, s_e=0.5) == 1.0

    def test_E_zero_at_hygroscopic_point(self):
        assert _laio_stress_E(0.05, s_h=0.05, s_e=0.5) == 0.0

    def test_E_zero_below_hygroscopic_point(self):
        assert _laio_stress_E(0.01, s_h=0.05, s_e=0.5) == 0.0

    def test_E_midpoint_gives_half_stress(self):
        s_h, s_e = 0.05, 0.5
        s_mid = 0.5 * (s_h + s_e)
        assert math.isclose(_laio_stress_E(s_mid, s_h, s_e), 0.5, abs_tol=1e-14)

    def test_E_monotone_increasing(self):
        s_h, s_e = 0.05, 0.5
        sats = [0.0, 0.04, 0.05, 0.1, 0.3, 0.5, 0.8]
        stresses = [_laio_stress_E(s, s_h, s_e) for s in sats]
        assert all(stresses[i] <= stresses[i + 1] for i in range(len(stresses) - 1))
