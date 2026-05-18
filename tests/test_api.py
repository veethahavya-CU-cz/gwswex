"""Unit tests for gwswex Python API: Pydantic validators, broadcasting, type guards.

These tests exercise pure-Python logic — no Fortran kernel is loaded.
"""

import numpy as np
import pytest

from gwswex.config import (
    ETStressParams,
    InitialConditions,
    Material,
    SolverConfig,
    SpatialDomain,
    TemporalDomain,
    VanGenuchtenParams,
    Vegetation,
    _broadcast_to_ne,
    _broadcast_to_nl_ne,
)


# ---------------------------------------------------------------------------
# _broadcast_to_ne
# ---------------------------------------------------------------------------
class TestBroadcastToNe:
    def test_scalar(self):
        out = _broadcast_to_ne(3.14, 5, "x")
        assert out.shape == (5,)
        assert np.all(out == pytest.approx(3.14))

    def test_length_one_list(self):
        out = _broadcast_to_ne([7.0], 3, "x")
        assert out.shape == (3,)
        assert np.all(out == pytest.approx(7.0))

    def test_matching_length(self):
        v = [1.0, 2.0, 3.0]
        out = _broadcast_to_ne(v, 3, "x")
        np.testing.assert_array_equal(out, v)

    def test_wrong_length_raises(self):
        with pytest.raises(ValueError):
            _broadcast_to_ne([1.0, 2.0], 3, "x")


# ---------------------------------------------------------------------------
# _broadcast_to_nl_ne
# ---------------------------------------------------------------------------
class TestBroadcastToNlNe:
    def test_scalar(self):
        out = _broadcast_to_nl_ne(1.0, nl=3, ne=4, name="x")
        assert out.shape == (3, 4)
        assert np.all(out == pytest.approx(1.0))

    def test_1d_per_layer(self):
        out = _broadcast_to_nl_ne([1.0, 2.0, 3.0], nl=3, ne=2, name="x")
        assert out.shape == (3, 2)
        np.testing.assert_array_equal(out[:, 0], [1, 2, 3])
        np.testing.assert_array_equal(out[:, 1], [1, 2, 3])

    def test_2d_ne_nl_transposed(self):
        arr = np.array([[1, 2, 3], [4, 5, 6]])  # (ne=2, nl=3)
        out = _broadcast_to_nl_ne(arr, nl=3, ne=2, name="x")
        assert out.shape == (3, 2)
        np.testing.assert_array_equal(out, arr.T)

    def test_2d_nl_ne_passthrough(self):
        arr = np.ones((3, 2))
        out = _broadcast_to_nl_ne(arr, nl=3, ne=2, name="x")
        assert out.shape == (3, 2)

    def test_wrong_shape_raises(self):
        with pytest.raises(ValueError):
            _broadcast_to_nl_ne(np.ones((4, 5)), nl=3, ne=2, name="x")


# ---------------------------------------------------------------------------
# VanGenuchtenParams
# ---------------------------------------------------------------------------
class TestVanGenuchtenParams:
    def test_derived_m(self):
        vg = VanGenuchtenParams(alpha=3.6, n=1.56, theta_r=0.078, theta_s=0.43)
        assert vg.m == pytest.approx(1.0 - 1.0 / 1.56)

    def test_Sy(self):
        vg = VanGenuchtenParams(alpha=3.6, n=1.56, theta_r=0.078, theta_s=0.43)
        assert vg.Sy == pytest.approx(0.43 - 0.078)

    def test_theta_s_gt_theta_r(self):
        with pytest.raises(Exception):
            VanGenuchtenParams(alpha=3.6, n=1.56, theta_r=0.5, theta_s=0.3)

    def test_n_gt_1(self):
        with pytest.raises(Exception):
            VanGenuchtenParams(alpha=3.6, n=0.9, theta_r=0.078, theta_s=0.43)


# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------
class TestMaterial:
    def test_basic(self):
        m = Material(
            id=1, name="topsoil", K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43}
        )
        assert m.K_sat == pytest.approx(1e-5)
        assert m.lam == pytest.approx(0.5)

    def test_k_sat_positive(self):
        with pytest.raises(Exception):
            Material(id=1, K_sat=-1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})


# ---------------------------------------------------------------------------
# SolverConfig
# ---------------------------------------------------------------------------
class TestSolverConfig:
    def test_defaults(self):
        sc = SolverConfig()
        assert sc.solver == "implicit"
        assert sc.omp_threads == 1
        assert sc.picard_tol == pytest.approx(1e-6)

    def test_solver_type_id_explicit(self):
        sc = SolverConfig(solver="explicit")
        assert sc.solver_type_id == 1

    def test_solver_type_id_implicit(self):
        sc = SolverConfig(solver="implicit")
        assert sc.solver_type_id == 2

    def test_invalid_solver_name(self):
        with pytest.raises(Exception):
            SolverConfig(solver="runge-kutta")

    def test_unknown_kwarg_rejected(self):
        """extra='forbid' must catch mistyped parameter names."""
        with pytest.raises(Exception):
            SolverConfig(picard_tolerance=1e-4)  # correct name is picard_tol

    def test_courant_bounded(self):
        with pytest.raises(Exception):
            SolverConfig(courant_number=1.5)

    def test_Ss_not_accepted(self):
        """Ss was removed; passing it should raise with extra='forbid'."""
        with pytest.raises(Exception):
            SolverConfig(Ss=1e-4)


# ---------------------------------------------------------------------------
# Vegetation
# ---------------------------------------------------------------------------
class TestVegetation:
    def test_invalid_growth_model(self):
        with pytest.raises(Exception):
            Vegetation(id=1, root_growth_model="power_law")


# ---------------------------------------------------------------------------
# ETStressParams
# ---------------------------------------------------------------------------
class TestETStressParams:
    def test_defaults(self):
        e = ETStressParams()
        assert 0 <= e.s_h <= e.s_w <= e.s_star <= 1

    def test_custom(self):
        e = ETStressParams(s_star=0.5, s_w=0.1, s_h=0.05, s_e=0.5)
        assert e.s_w == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# TemporalDomain
# ---------------------------------------------------------------------------
class TestTemporalDomain:
    def test_basic(self):
        td = TemporalDomain(n_steps=100, dt=3600.0)
        assert td.n_steps == 100
        assert td.dt == pytest.approx(3600.0)

    def test_steps_property(self):
        td = TemporalDomain(n_steps=5, dt=3600.0)
        assert list(td.steps) == [0, 1, 2, 3, 4]

    def test_timedelta_dt(self):
        from datetime import timedelta

        td = TemporalDomain(n_steps=24, dt=timedelta(hours=1))
        assert td.dt == pytest.approx(3600.0)

    def test_timedelta_dt_min(self):
        from datetime import timedelta

        td = TemporalDomain(n_steps=24, dt=timedelta(hours=1), dt_min=timedelta(minutes=5))
        assert td.dt_min == pytest.approx(300.0)


# ---------------------------------------------------------------------------
# Freezable: mutation after freeze
# ---------------------------------------------------------------------------
class TestFreezable:
    def test_freeze_blocks_mutation(self):
        sc = SolverConfig()
        sc.freeze()
        with pytest.raises((AttributeError, ValueError, Exception)):
            sc.picard_tol = 1e-3

    def test_unfreeze_allows_mutation(self):
        sc = SolverConfig()
        sc.freeze()
        sc.unfreeze()
        sc.picard_tol = 1e-3  # should not raise
        assert sc.picard_tol == pytest.approx(1e-3)


# ---------------------------------------------------------------------------
# RootParams
# ---------------------------------------------------------------------------
from gwswex.config import (
    InitialConditions,
    LateralFluxes,
    ModelParams,
    RootParams,
    RootGrowthModel,
)


class TestRootParams:
    def test_basic(self):
        rp = RootParams(depth=0.5)
        assert rp.depth == pytest.approx(0.5)

    def test_depth_must_be_positive(self):
        with pytest.raises(Exception):
            RootParams(depth=-0.1)

    def test_depth_zero_rejected(self):
        with pytest.raises(Exception):
            RootParams(depth=0.0)


class TestRootGrowthModel:
    def test_default_is_static(self):
        assert RootGrowthModel().model == "static"


# ---------------------------------------------------------------------------
# ModelParams
# ---------------------------------------------------------------------------
class TestModelParams:
    def test_defaults(self):
        mp = ModelParams()
        assert mp.psi_f > 0 and mp.F_min > 0 and 0 <= mp.ICratio_min <= 1

    def test_psi_f_must_be_positive(self):
        with pytest.raises(Exception):
            ModelParams(psi_f=-0.1)

    def test_F_min_must_be_positive(self):
        with pytest.raises(Exception):
            ModelParams(F_min=0.0)

    def test_ICratio_min_bounded(self):
        with pytest.raises(Exception):
            ModelParams(ICratio_min=1.5)
        with pytest.raises(Exception):
            ModelParams(ICratio_min=-0.1)


# ---------------------------------------------------------------------------
# InitialConditions
# ---------------------------------------------------------------------------
class TestInitialConditions:
    def test_basic(self):
        ic = InitialConditions(
            gw=np.array([-0.5]),
            sw=np.array([0.0]),
            uz=np.array([[-999.0], [-999.0], [-999.0]]),
        )
        assert ic.gw.shape == (1,)
        assert ic.uz.shape == (3, 1)


# ---------------------------------------------------------------------------
# LateralFluxes
# ---------------------------------------------------------------------------
class TestLateralFluxes:
    def test_basic(self):
        lf = LateralFluxes(gw=np.zeros(3), sw=np.zeros(3))
        assert lf.gw.shape == (3,) and lf.sw.shape == (3,)


# ---------------------------------------------------------------------------
# SpatialDomain validators
# ---------------------------------------------------------------------------
from gwswex.config import SpatialDomain


class TestSpatialDomain:
    def test_top_bot_assembles_bnds(self):
        sd = SpatialDomain(ne=2, nl=3, top=[1.0, 1.0], bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        assert sd.bnds.shape == (2, 4)
        np.testing.assert_array_equal(sd.bnds[0], [1.0, 0.5, 0.0, -1.0])

    def test_explicit_bnds(self):
        bnds = np.array([[1.0, 0.5, 0.0, -1.0]])
        sd = SpatialDomain(ne=1, nl=3, bnds=bnds, sID=1, vID=1)
        assert sd.bnds.shape == (1, 4)

    def test_bnds_must_be_decreasing(self):
        bnds = np.array([[1.0, 0.5, 0.6, -1.0]])
        with pytest.raises(Exception):
            SpatialDomain(ne=1, nl=3, bnds=bnds, sID=1, vID=1)

    def test_bnds_shape_validated(self):
        with pytest.raises(Exception):
            SpatialDomain(
                ne=1,
                nl=3,
                bnds=np.array([[1.0, 0.5, 0.0]]),
                sID=1,
                vID=1,
            )

    def test_missing_top_or_bot_raises(self):
        with pytest.raises(Exception):
            SpatialDomain(ne=1, nl=3, top=1.0, sID=1, vID=1)

    def test_sID_orientation_ne_nl(self):
        sd = SpatialDomain(
            ne=2,
            nl=3,
            top=1.0,
            bot=[0.5, 0.0, -1.0],
            sID=[[1, 1, 1], [2, 2, 2]],
            vID=1,
        )
        assert sd.sID.shape == (3, 2)
        np.testing.assert_array_equal(sd.sID[:, 0], [1, 1, 1])
        np.testing.assert_array_equal(sd.sID[:, 1], [2, 2, 2])

    def test_sID_orientation_nl_ne(self):
        sd = SpatialDomain(
            ne=2,
            nl=3,
            top=1.0,
            bot=[0.5, 0.0, -1.0],
            sID=np.array([[1, 2], [1, 2], [1, 2]]),
            vID=1,
        )
        assert sd.sID.shape == (3, 2)

    def test_vID_scalar(self):
        sd = SpatialDomain(ne=4, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=2)
        assert sd.vID.shape == (4,)
        assert np.all(sd.vID == 2)

    def test_vID_size_mismatch(self):
        with pytest.raises(Exception):
            SpatialDomain(ne=4, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=[1, 2, 3])

    def test_ne_must_be_positive(self):
        with pytest.raises(Exception):
            SpatialDomain(ne=0, nl=1, top=1.0, bot=[0.0], sID=1, vID=1)

    def test_nl_must_be_positive(self):
        with pytest.raises(Exception):
            SpatialDomain(ne=1, nl=0, top=1.0, bot=[], sID=1, vID=1)


# ---------------------------------------------------------------------------
# TemporalDomain validators
# ---------------------------------------------------------------------------
class TestTemporalDomainValidators:
    def test_dt_must_be_positive(self):
        with pytest.raises(Exception):
            TemporalDomain(n_steps=10, dt=0.0)
        with pytest.raises(Exception):
            TemporalDomain(n_steps=10, dt=-1.0)

    def test_n_steps_must_be_positive(self):
        with pytest.raises(Exception):
            TemporalDomain(n_steps=0, dt=1.0)


# ---------------------------------------------------------------------------
# SolverConfig: additional safety checks
# ---------------------------------------------------------------------------
class TestSolverConfigBounds:
    def test_picard_tol_must_be_positive(self):
        with pytest.raises(Exception):
            SolverConfig(picard_tol=0.0)
        with pytest.raises(Exception):
            SolverConfig(picard_tol=-1e-9)

    def test_picard_max_iter_min_one(self):
        with pytest.raises(Exception):
            SolverConfig(picard_max_iter=0)

    def test_omp_threads_min_one(self):
        with pytest.raises(Exception):
            SolverConfig(omp_threads=0)

    def test_n_trapz_min_four(self):
        with pytest.raises(Exception):
            SolverConfig(n_trapz=3)

    def test_courant_must_be_positive(self):
        with pytest.raises(Exception):
            SolverConfig(courant_number=0.0)

    def test_beta_hyst_bounded(self):
        with pytest.raises(Exception):
            SolverConfig(beta_hyst=0.0)
        with pytest.raises(Exception):
            SolverConfig(beta_hyst=1.5)


# ---------------------------------------------------------------------------
# Material: van Genuchten propagation
# ---------------------------------------------------------------------------
class TestMaterialValidation:
    def test_invalid_vg_propagates(self):
        with pytest.raises(Exception):
            Material(id=1, K_sat=1e-5, vanG={"alpha": -3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})

    def test_id_must_be_positive(self):
        with pytest.raises(Exception):
            Material(id=0, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})


# ---------------------------------------------------------------------------
# Vegetation: ET stress + root configurations
# ---------------------------------------------------------------------------
class TestVegetationValidation:
    def test_id_must_be_positive(self):
        with pytest.raises(Exception):
            Vegetation(id=0)

    def test_default_growth_model_static(self):
        v = Vegetation(id=1)
        assert v.root_growth_model == "static"

    def test_root_depth_initial_must_be_positive(self):
        with pytest.raises(Exception):
            Vegetation(id=1, root_depth_initial=-0.1)


# ---------------------------------------------------------------------------
# GWSWEXmodel: pure-Python lifecycle (no kernel)
# ---------------------------------------------------------------------------
try:
    from gwswex import GWSWEXmodel

    _HAS_GWSWEX = True
except Exception:
    _HAS_GWSWEX = False


@pytest.mark.skipif(not _HAS_GWSWEX, reason="GWSWEXmodel not importable")
class TestGWSWEXmodelLifecycle:
    def _make_model(self):
        m = GWSWEXmodel(name="api_test", T="s", L="m", write_output=False)
        m.init_space(ne=2, nl=3, top=[1.0, 1.0], bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m.add_material(
            id=1,
            K_sat=1e-5,
            vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
        )
        m.add_vegetation(id=1, root={"depth": 2.0})
        m.init_time(n_steps=4, dt=3600.0, dt_min=1.0)
        return m

    def test_invalid_unit_raises(self):
        with pytest.raises(KeyError):
            GWSWEXmodel(T="parsec")
        with pytest.raises(KeyError):
            GWSWEXmodel(L="furlong")

    def test_unit_scales(self):
        m = GWSWEXmodel(T="h", L="cm")
        assert m._T_scale == pytest.approx(3600.0)
        assert m._L_scale == pytest.approx(0.01)

    def test_add_material_before_init_space_raises(self):
        m = GWSWEXmodel()
        with pytest.raises(RuntimeError):
            m.add_material(
                id=1,
                K_sat=1e-5,
                vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
            )

    def test_set_initial_conditions_before_space_raises(self):
        m = GWSWEXmodel()
        with pytest.raises(RuntimeError):
            m.set_initial_conditions(gw=0.0, sw=0.0, uz=-999)

    def test_register_space_unknown_material_raises(self):
        m = GWSWEXmodel()
        m.init_space(ne=1, nl=2, top=1.0, bot=[0.5, 0.0], sID=99, vID=1)
        m.add_material(
            id=1,
            K_sat=1e-5,
            vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
        )
        m.add_vegetation(id=1, root={"depth": 1.5})
        with pytest.raises(ValueError):
            m._register_space()

    def test_register_space_unknown_vegetation_raises(self):
        m = GWSWEXmodel()
        m.init_space(ne=1, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=99)
        m.add_material(
            id=1,
            K_sat=1e-5,
            vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
        )
        m.add_vegetation(id=1, root={"depth": 1.5})
        with pytest.raises(ValueError):
            m._register_space()

    def test_set_forcing_before_register_raises(self):
        m = GWSWEXmodel()
        with pytest.raises(RuntimeError):
            m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)

    def test_set_forcing_scalar_broadcast(self):
        m = self._make_model()
        m.set_forcing(precip=1e-6, pet=0.0, ptt=0.0)
        assert m._forcing["precip"].shape == (4, 2)
        assert np.all(m._forcing["precip"] == pytest.approx(1e-6))

    def test_set_forcing_per_element(self):
        m = self._make_model()
        m.set_forcing(precip=[1e-6, 2e-6], pet=0.0, ptt=0.0)
        assert m._forcing["precip"].shape == (4, 2)
        np.testing.assert_array_equal(m._forcing["precip"][0], [1e-6, 2e-6])

    def test_set_forcing_per_step(self):
        m = self._make_model()
        m.set_forcing(precip=[1e-6, 2e-6, 3e-6, 4e-6], pet=0.0, ptt=0.0)
        assert m._forcing["precip"].shape == (4, 2)
        np.testing.assert_array_equal(m._forcing["precip"][:, 0], [1e-6, 2e-6, 3e-6, 4e-6])

    def test_set_forcing_full_2d(self):
        m = self._make_model()
        arr = np.arange(8.0).reshape(4, 2)
        m.set_forcing(precip=arr, pet=0.0, ptt=0.0)
        np.testing.assert_array_equal(m._forcing["precip"], arr)

    def test_set_forcing_bad_shape_raises(self):
        m = self._make_model()
        with pytest.raises(ValueError):
            m.set_forcing(precip=np.zeros((3, 3)), pet=0.0, ptt=0.0)

    def test_set_solver_after_init_blocked(self):
        m = self._make_model()
        m._is_initialised = True
        with pytest.raises(RuntimeError):
            m.set_solver(solver="explicit")
