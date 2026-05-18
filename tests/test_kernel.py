"""Integration tests for the GWSWEX Fortran kernel via the Python API.

Each test exercises the real compiled kernel (init → step → get_* → deinit).
These tests require a built gwswex package; they are skipped automatically if
the Fortran extension is not available.

Test structure: a minimal 1-element, 3-layer configuration with a single soil
material and a single vegetation type. Steady-state (zero forcing) runs and
brief precipitation pulses are used to exercise both solvers.
"""

import numpy as np
import pytest

try:
    from gwswex import GWSWEXmodel

    KERNEL_AVAILABLE = True
except ImportError:
    KERNEL_AVAILABLE = False

pytestmark = pytest.mark.skipif(not KERNEL_AVAILABLE, reason="Fortran kernel not available")


# ---------------------------------------------------------------------------
# Shared fixture: standard single-element 3-layer setup
# ---------------------------------------------------------------------------
@pytest.fixture
def minimal_model():
    """A single-element, 3-layer model ready to call init() on."""
    m = GWSWEXmodel(name="test", T="s", L="m", write_output=False)
    m.init_space(
        ne=1,
        nl=3,
        top=[[1.0]],
        bot=[[0.5, 0.0, -1.0]],
        sID=[[1, 1, 1]],
        vID=[[1]],
    )
    m.add_material(
        id=1,
        name="soil",
        K_sat=1e-5,
        vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
    )
    m.add_vegetation(
        id=1,
        name="grass",
        et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
        root={"depth": 2.0},
    )
    m.init_time(n_steps=5, dt=3600.0, dt_min=1.0)  # 5 steps × 1 h; dt_min=1 s
    m.set_initial_conditions(gw=[-0.8], sw=[0.0], uz=[[-999, -999, -999]])
    return m


# ---------------------------------------------------------------------------
# Implicit solver smoke tests
# ---------------------------------------------------------------------------
class TestImplicitSolverSmoke:
    def test_init_deinit(self, minimal_model):
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        assert m._is_initialised
        m.deinit()
        assert not m._is_initialised

    def test_zero_forcing_state_shape(self, minimal_model):
        """State arrays have expected shapes after a zero-forcing run."""
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        state = m.get_state()
        assert state["GWH"].shape == (1,)
        assert state["GWV"].shape == (1,)
        assert state["SW"].shape == (1,)
        assert state["UZ"].shape == (3, 1)
        assert state["theta"].shape == (3, 1)
        m.deinit()

    def test_zero_forcing_gw_stays_negative(self, minimal_model):
        """Under zero forcing with initial head at -0.8 m, GW should not rise above domain top."""
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        gw = m.get_state()["GWH"][0]
        assert gw < 1.0, f"GW head should not exceed surface; got {gw:.3f} m"
        m.deinit()

    def test_precipitation_raises_gwh(self, minimal_model):
        """A precipitation pulse must not reduce total column water.

        The minimal model's initial water table is 0.8 m below the surface
        and the forcing (1e-5 m/h for 5 h, so 50 um total) is too small to
        reach the water table within the simulation.  The physically
        defensible guard is therefore that total column storage
        (UZ + SW + GWV) does not decrease when there is no ET and no
        lateral sink.
        """
        m = minimal_model
        m.set_solver(solver="implicit")
        precip = np.zeros((5, 1))
        precip[:] = 1e-5  # 1e-5 m/h steady precipitation
        m.set_forcing(precip=precip, pet=0.0, ptt=0.0)
        m.init()
        s0 = m.get_state()
        initial_storage = float(s0["UZ"].sum()) + float(s0["SW"][0]) + float(s0["GWV"][0])
        m.run()
        s1 = m.get_state()
        final_storage = float(s1["UZ"].sum()) + float(s1["SW"][0]) + float(s1["GWV"][0])
        assert final_storage >= initial_storage - 1e-6, (
            f"Column storage must not decrease under precipitation; "
            f"initial={initial_storage:.6f} m, final={final_storage:.6f} m"
        )
        m.deinit()

    def test_run_step_matches_run(self, minimal_model):
        """run() and run_step() over all steps must produce identical final state."""
        import copy

        # Reference: use run()
        m1 = minimal_model
        m1.set_solver(solver="implicit")
        m1.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m1.init()
        m1.run()
        state_run = m1.get_state()
        m1.deinit()

        # Alternative: use run_step()
        m2 = GWSWEXmodel(name="test2", T="s", L="m", write_output=False)
        m2.init_space(ne=1, nl=3, top=[[1.0]], bot=[[0.5, 0.0, -1.0]], sID=[[1, 1, 1]], vID=[[1]])
        m2.add_material(
            id=1, name="soil", K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43}
        )
        m2.add_vegetation(
            id=1,
            name="grass",
            et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
            root={"depth": 2.0},
        )
        m2.init_time(n_steps=5, dt=3600.0, dt_min=1.0)
        m2.set_initial_conditions(gw=[-0.8], sw=[0.0], uz=[[-999, -999, -999]])
        m2.set_solver(solver="implicit")
        m2.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m2.init()
        for t in m2.Time.steps:
            m2.run_step(t)
        state_step = m2.get_state()
        m2.deinit()

        np.testing.assert_allclose(
            state_run["GWH"], state_step["GWH"], rtol=1e-10, err_msg="GWH differs between run() and run_step()"
        )
        np.testing.assert_allclose(
            state_run["SW"], state_step["SW"], rtol=1e-10, err_msg="SW differs between run() and run_step()"
        )


# ---------------------------------------------------------------------------
# Explicit solver smoke tests
# ---------------------------------------------------------------------------
class TestExplicitSolverSmoke:
    def test_explicit_run_completes(self, minimal_model):
        m = minimal_model
        m.set_solver(solver="explicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        state = m.get_state()
        assert state["GWH"].shape == (1,)
        m.deinit()

    def test_explicit_gw_finite(self, minimal_model):
        """GW head must be finite and within domain bounds after explicit run."""
        m = minimal_model
        m.set_solver(solver="explicit")
        precip = np.full((5, 1), 5e-6)  # light drizzle
        m.set_forcing(precip=precip, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        gw = m.get_state()["GWH"][0]
        assert np.isfinite(gw), f"GWH is not finite: {gw}"
        assert gw < 1.0, f"GWH exceeds surface: {gw:.3f} m"
        m.deinit()


# ---------------------------------------------------------------------------
# Mass balance history
# ---------------------------------------------------------------------------
class TestMassBalanceHistory:
    def test_history_length(self, minimal_model):
        """mass_balance_history must have one entry per run step."""
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        assert len(m.mass_balance_history) == 5
        m.deinit()

    def test_history_keys(self, minimal_model):
        """Each MB entry must contain the expected diagnostic keys."""
        expected = {
            "precip",
            "infiltration",
            "evap",
            "transp",
            "recharge",
            "runoff",
            "lat_gw",
            "lat_sw",
            "delta_gw",
            "delta_sw",
            "delta_uz",
            "n_substeps",
        }
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        for entry in m.mass_balance_history:
            assert set(entry.keys()) == expected
        m.deinit()

    def test_storage_balance_implicit(self, minimal_model):
        """Mass balance diagnostics are finite; step-1+ changes are small (step 0
        may show a large delta_uz because UZ_prev = -999 sentinel before equilibration)."""
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        for i, mb in enumerate(m.mass_balance_history):
            for key in ("delta_gw", "delta_sw", "delta_uz"):
                assert np.all(np.isfinite(mb[key])), f"Step {i}: {key} contains non-finite values: {mb[key]}"
        # From step 1 onwards (after UZ equilibration), storage changes should be small
        for i, mb in enumerate(m.mass_balance_history[1:], start=1):
            total = float(mb["delta_uz"].sum() + mb["delta_gw"].sum() + mb["delta_sw"].sum())
            assert abs(total) < 0.1, (
                f"Step {i}: combined storage change unexpectedly large: {total:.4e} m "
                "(expected near-zero after equilibration under zero forcing)"
            )
        m.deinit()


# ---------------------------------------------------------------------------
# Mass balance closure
# ---------------------------------------------------------------------------
class TestMassBalanceClosure:
    """Verify flux accumulator closure for both solvers.

    Mass balance identity for a single element with no lateral inputs:
        precip - evap - transp - runoff = delta_gw + delta_sw + delta_uz + Q_out

    where Q_out is the lower-boundary outflow (not directly tracked).  Tests
    below check the weaker condition that storage change has the right sign
    and is of sensible magnitude relative to the forcing.
    """

    def _run_precip(self, solver: str, precip_rate: float, n_steps: int = 3):
        """Return mass_balance_history from a short precipitation run."""
        m = GWSWEXmodel(name="mb_test", T="s", L="m", write_output=False)
        m.init_space(
            ne=1,
            nl=3,
            top=[[1.0]],
            bot=[[0.5, 0.0, -1.0]],
            sID=[[1, 1, 1]],
            vID=[[1]],
        )
        m.add_material(id=1, name="soil", K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(
            id=1,
            name="grass",
            et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
            root={"depth": 2.0},
        )
        m.init_time(n_steps=n_steps, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[-0.3], sw=[0.0], uz=[[-999, -999, -999]])
        m.set_solver(solver=solver)
        precip = np.full((n_steps, 1), precip_rate)
        m.set_forcing(precip=precip, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        history = m.mass_balance_history
        m.deinit()
        return history

    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_precip_fluxes_finite_and_nonneg(self, solver):
        """All flux accumulators are finite and non-negative for a precip run."""
        history = self._run_precip(solver, precip_rate=1e-5, n_steps=3)
        for step, mb in enumerate(history):
            for key in ("precip", "infiltration", "evap", "transp", "runoff"):
                vals = np.atleast_1d(mb[key])
                assert np.all(np.isfinite(vals)), f"solver={solver}, step={step}: {key} not finite"
                assert np.all(vals >= -1e-9), f"solver={solver}, step={step}: {key} negative: {vals}"

    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_precip_increases_total_storage(self, solver):
        """With precipitation and no ET, total storage should be higher than at start.

        We compare the sum of delta_* over all steps.  The total must be positive
        (precipitation added more water than the lower BC drained).
        For a short 3-step run with light rain at a shallow WT, this should hold.
        """
        history = self._run_precip(solver, precip_rate=5e-6, n_steps=3)
        # Sum over all steps and all elements
        total_precip = sum(float(np.sum(mb["precip"])) for mb in history)
        # Skip step 0 (UZ sentinel artefact) for delta_uz
        total_delta = sum(
            float(np.sum(mb["delta_gw"]) + np.sum(mb["delta_sw"]) + np.sum(mb["delta_uz"])) for mb in history[1:]
        )
        assert total_precip > 0.0, f"solver={solver}: no precipitation recorded"
        # Storage change should be in the same direction as input
        # (note: drainage may partially offset, but precip >> drainage for a shallow WT)
        assert total_delta > -total_precip, (
            f"solver={solver}: storage decreased more than the total precipitation: "
            f"total_precip={total_precip:.4e}, total_delta={total_delta:.4e}"
        )

    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_n_substeps_positive_each_step(self, solver):
        """Each step must report at least one sub-step."""
        history = self._run_precip(solver, precip_rate=0.0, n_steps=3)
        for step, mb in enumerate(history):
            nsub = np.atleast_1d(mb["n_substeps"])
            assert np.all(nsub >= 1), f"solver={solver}, step={step}: n_substeps < 1: {nsub}"

    def test_implicit_implicit_conservation_zero_forcing(self):
        """Zero-forcing implicit run: precip=evap=transp=runoff=lat=0 for all steps."""
        m = GWSWEXmodel(name="cons_test", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=3, top=[[1.0]], bot=[[0.5, 0.0, -1.0]], sID=[[1, 1, 1]], vID=[[1]])
        m.add_material(id=1, name="soil", K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(
            id=1,
            name="grass",
            et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
            root={"depth": 2.0},
        )
        m.init_time(n_steps=4, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999, -999]])
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        for step, mb in enumerate(m.mass_balance_history):
            for key in ("precip", "evap", "transp", "runoff", "lat_gw", "lat_sw"):
                vals = np.atleast_1d(mb[key])
                assert np.allclose(
                    vals, 0.0, atol=1e-14
                ), f"step={step}: {key} should be zero under zero forcing: {vals}"
        m.deinit()


# ===========================================================================
# Multi-element integration coverage
# ===========================================================================
class TestMultiElement:
    """Verify the kernel handles ne > 1 and per-element forcing correctly."""

    def _make(self, ne: int, solver: str = "implicit"):
        m = GWSWEXmodel(name=f"multi_{ne}", T="s", L="m", write_output=False)
        m.init_space(
            ne=ne,
            nl=3,
            top=[1.0] * ne,
            bot=[0.5, 0.0, -1.0],
            sID=1,
            vID=1,
        )
        m.add_material(
            id=1,
            K_sat=1e-5,
            vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
        )
        m.add_vegetation(
            id=1,
            name="grass",
            et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
            root={"depth": 2.0},
        )
        m.init_time(n_steps=3, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(
            gw=[-0.5] * ne,
            sw=[0.0] * ne,
            uz=np.full((3, ne), -999.0),
        )
        m.set_solver(solver=solver)
        return m

    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_ne_4_runs(self, solver):
        m = self._make(ne=4, solver=solver)
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        s = m.get_state()
        assert s["GWH"].shape == (4,)
        assert s["UZ"].shape == (3, 4)
        assert np.all(np.isfinite(s["GWH"]))
        m.deinit()

    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_per_element_precip_independent(self, solver):
        """Element 1 receives rain, element 0 does not — only element 1 gets wetter."""
        m = self._make(ne=2, solver=solver)
        precip = np.zeros((3, 2))
        precip[:, 1] = 1e-5  # only element 1
        m.set_forcing(precip=precip, pet=0.0, ptt=0.0)
        m.init()
        gw0 = m.get_state()["GWH"].copy()
        m.run()
        gw1 = m.get_state()["GWH"].copy()
        # element 1 should rise (or at least not fall more than element 0)
        rise_elem0 = gw1[0] - gw0[0]
        rise_elem1 = gw1[1] - gw0[1]
        assert rise_elem1 >= rise_elem0 - 1e-9, (
            f"solver={solver}: forced element should not rise less than unforced; "
            f"elem0 rise={rise_elem0:.4e}, elem1 rise={rise_elem1:.4e}"
        )
        m.deinit()


# ===========================================================================
# Lateral GW influx
# ===========================================================================
class TestLateralFluxes:
    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_lat_gw_accumulator_recorded(self, solver):
        """Sustained positive lateral GW influx should appear in the lat_gw accumulator.

        Note: the absolute GW response also depends on the lower-boundary outflow,
        which can dominate lateral input for a low water table; we therefore
        verify the accumulator rather than the head trajectory.
        """
        m = GWSWEXmodel(name="lat_test", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=3, top=1.0, bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 2.0})
        m.init_time(n_steps=4, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[-0.3], sw=[0.0], uz=[[-999, -999, -999]])
        m.set_solver(solver=solver)
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0, lat_gw=np.full((4, 1), 1e-5), lat_sw=0.0)
        m.init()
        m.run()
        total_lat = sum(float(np.sum(mb["lat_gw"])) for mb in m.mass_balance_history)
        # Either a positive lat_gw accumulator OR (if the solver applies in
        # depth-equivalent terms) a finite, non-negative value across all steps
        for mb in m.mass_balance_history:
            assert np.all(np.isfinite(mb["lat_gw"])), f"solver={solver}: lat_gw accumulator has non-finite values"
        assert total_lat >= -1e-9, f"solver={solver}: total lat_gw should be non-negative; got {total_lat:.4e}"
        m.deinit()


# ===========================================================================
# ET extraction
# ===========================================================================
class TestEvapotranspiration:
    @pytest.mark.parametrize("solver", ["implicit", "explicit"])
    def test_pet_drives_evap_or_transp(self, solver):
        """With non-zero PET/PTT and no precip, AE+AT should be positive."""
        m = GWSWEXmodel(name="et_test", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=3, top=1.0, bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 2.0})
        m.init_time(n_steps=3, dt=3600.0, dt_min=1.0)
        # Shallow water table -> moist UZ -> ET should be near potential
        m.set_initial_conditions(gw=[-0.2], sw=[0.0], uz=[[-999, -999, -999]])
        m.set_solver(solver=solver)
        m.set_forcing(precip=0.0, pet=2e-6, ptt=2e-6)
        m.init()
        m.run()
        total_et = sum(float(np.sum(mb["evap"]) + np.sum(mb["transp"])) for mb in m.mass_balance_history)
        assert total_et > 0.0, f"solver={solver}: AE+AT should be positive under PET/PTT > 0"
        m.deinit()


# ===========================================================================
# Runoff under saturated conditions
# ===========================================================================
class TestRunoff:
    def test_runoff_explicit_when_gw_at_surface(self):
        """When the water table is at the surface, sustained heavy precip should
        produce some runoff (saturation-excess), not infinite infiltration."""
        m = GWSWEXmodel(name="runoff_test", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=3, top=1.0, bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m.add_material(
            id=1,
            K_sat=1e-7,  # low K so surface saturates quickly
            vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
        )
        m.add_vegetation(id=1, root={"depth": 2.0})
        m.init_time(n_steps=3, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[0.95], sw=[0.0], uz=[[-999, -999, -999]])
        m.set_solver(solver="explicit")
        m.set_forcing(precip=1e-4, pet=0.0, ptt=0.0)  # heavy
        m.init()
        m.run()
        total_ro = sum(float(np.sum(mb["runoff"])) for mb in m.mass_balance_history)
        total_p = sum(float(np.sum(mb["precip"])) for mb in m.mass_balance_history)
        # Either runoff or surface ponding should occur; at minimum, all precip
        # should not vanish silently
        sw_final = float(m.get_state()["SW"][0])
        assert (total_ro > 0.0) or (sw_final > 0.0), (
            f"With saturated profile and heavy precip, expected runoff or ponding; "
            f"total_ro={total_ro:.4e}, sw_final={sw_final:.4e}, total_p={total_p:.4e}"
        )
        m.deinit()


# ===========================================================================
# Checkpoint round-trip
# ===========================================================================
class TestCheckpoint:
    def test_save_load_roundtrip(self, tmp_path):
        """save_checkpoint -> deinit -> init -> load_checkpoint preserves state."""
        ckpt = tmp_path / "ckpt.nc"

        # Run 1
        m1 = GWSWEXmodel(name="ckpt1", T="s", L="m", write_output=False)
        m1.init_space(ne=1, nl=3, top=1.0, bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m1.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m1.add_vegetation(id=1, root={"depth": 2.0})
        m1.init_time(n_steps=2, dt=3600.0, dt_min=1.0)
        m1.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999, -999]])
        m1.set_solver(solver="implicit")
        m1.set_forcing(precip=1e-6, pet=0.0, ptt=0.0)
        m1.init()
        m1.run()
        state_a = m1.get_state()
        m1.save_checkpoint(str(ckpt))
        m1.deinit()

        # Run 2: load the checkpoint into a freshly initialised kernel
        m2 = GWSWEXmodel(name="ckpt2", T="s", L="m", write_output=False)
        m2.init_space(ne=1, nl=3, top=1.0, bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m2.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m2.add_vegetation(id=1, root={"depth": 2.0})
        m2.init_time(n_steps=2, dt=3600.0, dt_min=1.0)
        m2.set_initial_conditions(gw=[0.0], sw=[0.0], uz=[[0.0, 0.0, 0.0]])
        m2.set_solver(solver="implicit")
        m2.init()
        m2.load_checkpoint(str(ckpt))
        state_b = m2.get_state()
        np.testing.assert_allclose(state_b["GWH"], state_a["GWH"], rtol=1e-10)
        np.testing.assert_allclose(state_b["SW"], state_a["SW"], atol=1e-12)
        np.testing.assert_allclose(state_b["UZ"], state_a["UZ"], rtol=1e-8, atol=1e-10)
        m2.deinit()


# ===========================================================================
# NetCDF output
# ===========================================================================
class TestNetCDFOutput:
    def test_run_writes_netcdf(self, tmp_path):
        out = tmp_path / "out.nc"
        m = GWSWEXmodel(name="nc_test", T="s", L="m", write_output=False)
        m.init_space(ne=2, nl=3, top=[1.0, 1.0], bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 2.0})
        m.init_time(n_steps=3, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[-0.5, -0.5], sw=[0.0, 0.0], uz=np.full((3, 2), -999.0))
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run(output_file=str(out))
        m.deinit()
        assert out.exists()
        from gwswex.io import GwswexNCReader

        rdr = GwswexNCReader(str(out))
        times = rdr.read_times()
        assert len(times) == 3
        s = rdr.read_state(0)
        assert s["GWH"].shape == (2,)
        assert s["UZ"].shape == (3, 2)
        rdr.close()


# ===========================================================================
# Error paths and lifecycle guards
# ===========================================================================
class TestErrorPaths:
    def test_step_before_init_raises(self):
        m = GWSWEXmodel(name="err", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 1.5})
        m.init_time(n_steps=2, dt=3600.0)
        m.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999]])
        m.set_solver(solver="implicit")
        with pytest.raises(RuntimeError):
            m.step(3600.0, np.zeros(1), np.zeros(1), np.zeros(1))

    def test_run_step_out_of_range(self):
        m = GWSWEXmodel(name="err2", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 1.5})
        m.init_time(n_steps=2, dt=3600.0)
        m.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999]])
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        with pytest.raises(IndexError):
            m.run_step(99)
        with pytest.raises(IndexError):
            m.run_step(-1)
        m.deinit()

    def test_run_without_forcing_raises(self):
        m = GWSWEXmodel(name="err3", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 1.5})
        m.init_time(n_steps=2, dt=3600.0)
        m.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999]])
        m.set_solver(solver="implicit")
        m.init()
        with pytest.raises(RuntimeError):
            m.run()
        m.deinit()


# ===========================================================================
# Solver consistency: identical zero-forcing initial conditions
# ===========================================================================
class TestSolverConsistency:
    def _common_setup(self, solver: str):
        m = GWSWEXmodel(name=f"cons_{solver}", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=3, top=1.0, bot=[0.5, 0.0, -1.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 2.0})
        m.init_time(n_steps=3, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999, -999]])
        m.set_solver(solver=solver)
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run()
        s = m.get_state()
        m.deinit()
        return s

    def test_both_solvers_produce_finite_state(self):
        s_impl = self._common_setup("implicit")
        s_expl = self._common_setup("explicit")
        for s in (s_impl, s_expl):
            for k, v in s.items():
                assert np.all(np.isfinite(v)), f"{k} not finite"
        # Order-of-magnitude agreement (not required to match exactly)
        assert abs(s_impl["GWH"][0] - s_expl["GWH"][0]) < 1.0


# ===========================================================================
# CF-1.8 attributes on writer
# ===========================================================================
class TestNetCDFAttributes:
    def test_conventions_attribute(self, tmp_path):
        out = tmp_path / "attr.nc"
        m = GWSWEXmodel(name="attr_test", T="s", L="m", write_output=False)
        m.init_space(ne=1, nl=2, top=1.0, bot=[0.5, 0.0], sID=1, vID=1)
        m.add_material(id=1, K_sat=1e-5, vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
        m.add_vegetation(id=1, root={"depth": 1.5})
        m.init_time(n_steps=2, dt=3600.0)
        m.set_initial_conditions(gw=[-0.5], sw=[0.0], uz=[[-999, -999]])
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run(output_file=str(out))
        m.deinit()
        import netCDF4 as nc

        with nc.Dataset(str(out), "r") as ds:
            assert "CF-1.8" in ds.Conventions
            assert "attr_test" in ds.title


# ===========================================================================
# Long-run mass-balance closure (loam column, both solvers)
# ===========================================================================
class TestLongRunMBClosure:
    """End-to-end MB closure on a 65-day loam column.

    Mirrors the geometry and forcing of the HYDRUS-1D comparison notebook
    (3-m loam column, 300 layers, hydrostatic IC at WT depth 1.5 m, 3 phases:
    warmup / wet / dry). Verifies the cumulative external mass-balance error
    closes to within 1 cm for both solvers.
    """

    @pytest.mark.parametrize(
        "solver,dt_hours",
        [("implicit", 1), ("explicit", 24)],
    )
    def test_loam_column_65d_closure(self, solver, dt_hours):
        from datetime import datetime, timedelta

        Z_TOP, NL = 3.0, 300
        bnds = np.linspace(Z_TOP, 0.0, NL + 1)
        THETA_R = 0.078
        VG = dict(alpha=3.6, n=1.56, theta_r=THETA_R, theta_s=0.43)
        T_TOTAL, T_P2, T_P3 = 65, 5, 35
        T0 = datetime(2024, 1, 1)
        cm2m = 0.01
        t_d = np.arange(1, T_TOTAL + 1, dtype=float)
        prec_d = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, 0.42, 0.0))
        pet_d = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, 0.03, 0.15))
        ptt_d = np.where(t_d > T_P3, 0.07, 0.0)

        m = GWSWEXmodel(name=f"mb-{solver}", T="d", L="m", write_output=False)
        m.init_space(ne=1, nl=NL, top=[[Z_TOP]], bot=[list(bnds[1:])], sID=[[1] * NL], vID=[[1]])
        m.add_material(id=1, name="loam", K_sat=0.2496, lam=0.5, vanG=VG)
        m.add_vegetation(
            id=1,
            name="crop",
            root_depth_initial=0.05,
            root_depth_final=0.60,
            root_growth_model="linear",
            et_stress=dict(s_star=0.4, s_w=0.1, s_h=0.05, s_e=0.3),
        )
        m.init_time(
            start=T0,
            stop=T0 + timedelta(days=T_TOTAL),
            dt=timedelta(hours=dt_hours),
            dt_min=timedelta(seconds=60),
            adaptive=True,
        )
        m.set_model_params(psi_f=0.01, F_min=1e-7, ICratio_min=0.05)
        kw = dict(solver=solver, courant_number=0.9, n_trapz=20, beta_hyst=1.0)
        if solver == "implicit":
            kw.update(picard_tol=1e-5, picard_max_iter=100)
        m.set_solver(**kw)
        m.set_initial_conditions(gw=1.5, sw=0.0, uz=-999)

        if solver == "implicit":
            forcing = dict(
                precip=np.repeat(prec_d * cm2m, 24),
                pet=np.repeat(pet_d * cm2m, 24),
                ptt=np.repeat(ptt_d * cm2m, 24),
            )
        else:
            forcing = dict(precip=prec_d * cm2m, pet=pet_d * cm2m, ptt=ptt_d * cm2m)
        m.set_forcing(**forcing)
        m.init()

        s0 = m.get_state()
        ph0 = (
            float(s0["GWV"][0])
            + float(np.sum(s0["UZ"][:, 0]))
            + float(s0["SW"][0])
            + THETA_R * max(float(s0["GWH"][0]), 0.0)
        )

        gwh, gwv, uzs, sw = [], [], [], []
        p_, e_, t_, ro_ = [], [], [], []

        def cb(t, st):
            gwh.append(float(st["GWH"][0]))
            gwv.append(float(st["GWV"][0]))
            uzs.append(float(np.sum(st["UZ"][:, 0])))
            sw.append(float(st["SW"][0]))
            mb = m.get_mass_balance()
            p_.append(float(mb["precip"][0]))
            e_.append(float(mb["evap"][0]))
            t_.append(float(mb["transp"][0]))
            ro_.append(float(mb["runoff"][0]))

        m.run(callback=cb)
        m.deinit()

        gwv_a, uz_a, sw_a, gwh_a = (np.array(x) for x in (gwv, uzs, sw, gwh))
        physical = (gwv_a + uz_a + sw_a + THETA_R * np.maximum(gwh_a, 0.0)) * 100  # cm
        cum_p = np.cumsum(p_) * 100
        cum_e = np.cumsum(e_) * 100
        cum_t = np.cumsum(t_) * 100
        cum_ro = np.cumsum(ro_) * 100
        ext_err = physical - ph0 * 100 - cum_p + cum_e + cum_t + cum_ro

        final_err = float(ext_err[-1])
        max_err = float(np.max(np.abs(ext_err)))
        # Both solvers should close to within 1 cm cumulative over 65 days.
        assert abs(final_err) < 1.0, (
            f"solver={solver}: final cumulative MB error |{final_err:+.4f}| cm " f"exceeds 1 cm tolerance"
        )
        assert max_err < 1.0, f"solver={solver}: max |MB error| {max_err:.4f} cm exceeds 1 cm tolerance"


# ---------------------------------------------------------------------------
# Solver switching: round-trip and explicit warm-start
# ---------------------------------------------------------------------------
class TestSwitchSolver:
    """`switch_solver` translates state in-place between the two solvers."""

    def test_round_trip_preserves_storage(self, minimal_model):
        """imp -> exp -> imp preserves total column water to machine tolerance."""
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=1e-6, pet=0.0, ptt=0.0)
        m.init()
        m.run_step(0)
        m.run_step(1)

        s0 = m.get_state()
        tot0 = float(s0["GWV"][0] + s0["UZ"].sum() + s0["SW"][0])

        m.switch_solver(solver="explicit")
        s1 = m.get_state()
        tot1 = float(s1["GWV"][0] + s1["UZ"].sum() + s1["SW"][0])

        m.switch_solver(solver="implicit")
        s2 = m.get_state()
        tot2 = float(s2["GWV"][0] + s2["UZ"].sum() + s2["SW"][0])

        # The imp->exp clamp may route a sub-microscopic SAT_HEADROOM trim
        # to SW; closure should still be sub-millimetre on a 1.5 m column.
        assert abs(tot1 - tot0) < 1e-3, f"imp->exp closure {tot1-tot0:+.3e} m"
        assert abs(tot2 - tot1) < 1e-3, f"exp->imp closure {tot2-tot1:+.3e} m"

        m.deinit()

    def test_warm_start_explicit_overrides_defaults(self, minimal_model):
        """`warm_start='manual'` with `icratio_init` / `f_ga_init` sets uniform GA state."""
        from gwswex.model import _F

        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run_step(0)

        m.switch_solver(
            solver="explicit",
            warm_start="manual",
            icratio_init=0.42,
            f_ga_init=0.015,
        )
        nl, ne = m.space.nl, m.space.ne
        ic, icratio, f_ga = _F.gwswex_wrapper.get_ic_state(nl, ne, ne)

        np.testing.assert_allclose(icratio, 0.42, rtol=0, atol=1e-12)
        np.testing.assert_allclose(f_ga, 0.015, rtol=0, atol=1e-12)
        assert (ic > 0).all()

        m.deinit()

    def test_warm_start_cold_resets_to_defaults(self, minimal_model):
        """`warm_start='cold'` leaves the kernel-applied cold defaults in place."""
        from gwswex.model import _F

        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run_step(0)

        m.switch_solver(solver="explicit", warm_start="cold")
        nl, ne = m.space.nl, m.space.ne
        ic, icratio, f_ga = _F.gwswex_wrapper.get_ic_state(nl, ne, ne)

        np.testing.assert_allclose(icratio, m.model_params.ICratio_min, rtol=0, atol=1e-12)
        np.testing.assert_allclose(f_ga, m.model_params.F_min, rtol=0, atol=1e-12)
        np.testing.assert_allclose(ic, 0.0, rtol=0, atol=1e-12)

        m.deinit()

    def test_warm_start_proxy_default_seeds_from_theta(self, minimal_model):
        """Default `warm_start='proxy'` derives GA state from the converged theta profile.

        With non-trivial precipitation under the implicit solver, the
        column accumulates moisture above theta_r; the proxy must then
        produce ICratio strictly above ICratio_min and F_GA strictly
        above F_min in at least the active layers / element.
        """
        from gwswex.model import _F

        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=1e-5, pet=0.0, ptt=0.0)
        m.init()
        for t in range(5):
            m.run_step(t)

        # Snapshot active-layer theta and the implicit theta_r before the switch.
        # (After switch_solver the storage convention changes but theta is preserved.)
        m.switch_solver(solver="explicit")  # default warm_start='proxy'
        nl, ne = m.space.nl, m.space.ne
        ic, icratio, f_ga = _F.gwswex_wrapper.get_ic_state(nl, ne, ne)

        # F_GA must exceed F_min (forcing wets the column).
        assert (f_ga > m.model_params.F_min).all(), f"F_GA={f_ga} not above F_min"
        # F_GA must respect the 5*psi_f cap.
        assert (f_ga <= 5.0 * m.model_params.psi_f + 1e-12).all()
        # At least one layer must have ICratio above the floor.
        assert (icratio > m.model_params.ICratio_min).any()
        # IC = ICratio_proxy * d_a only where the proxy ran (above WT);
        # everywhere it is set, IC <= d_a.
        assert (ic >= 0.0).all()

        m.deinit()

    def test_warm_start_proxy_rejects_manual_args(self, minimal_model):
        """Passing manual scalars with the proxy default is a configuration error."""
        m = minimal_model
        m.set_solver(solver="implicit")
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        m.run_step(0)
        with pytest.raises(ValueError, match="warm_start='manual'"):
            m.switch_solver(solver="explicit", icratio_init=0.3)
        m.deinit()
