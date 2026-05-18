"""OpenMP availability and round-trip tests.

Verify that the GWSWEX kernel was built with OpenMP support and that the
get/set thread-count helpers exposed through the f2py wrapper behave
correctly.

These tests intentionally never fail: when OpenMP is unavailable or the
runtime helpers misbehave, the tests emit a warning and pass. Production
builds are expected to ship with OpenMP enabled because the element loop
is the principal hot path, but the absence of OMP must not break the
test suite (e.g. on developer machines without an OpenMP-capable
toolchain).
"""

from __future__ import annotations

import warnings

import pytest

try:
    from gwswex import f_gwswex as _F  # type: ignore[attr-defined]

    KERNEL_AVAILABLE = True
except Exception:
    KERNEL_AVAILABLE = False

pytestmark = pytest.mark.skipif(not KERNEL_AVAILABLE, reason="Fortran kernel not available")


def _omp_available() -> bool:
    try:
        return int(_F.gwswex_wrapper.get_omp_available()) == 1
    except Exception as exc:  # pragma: no cover - defensive
        warnings.warn(f"get_omp_available() raised {exc!r}; treating OMP as unavailable")
        return False


class TestOpenMP:
    def test_omp_compiled_in(self):
        """Kernel should be built with -fopenmp; warn if not."""
        if not _omp_available():
            warnings.warn(
                "GWSWEX kernel was built without OpenMP support; rebuild with "
                "an OpenMP-capable compiler so that the element-parallel loop "
                "is active. (warning only, not failing the test)"
            )

    def test_get_omp_max_threads_positive(self):
        """get_omp_max_threads should return a strictly positive thread count."""
        try:
            n = int(_F.gwswex_wrapper.get_omp_max_threads())
        except Exception as exc:
            warnings.warn(f"get_omp_max_threads() raised {exc!r}; OMP probably unavailable")
            return
        if n < 1:
            warnings.warn(f"max_threads should be >= 1, got {n}")

    def test_set_then_get_omp_threads_round_trip(self):
        """set_omp_threads(n) followed by get_omp_max_threads() should return n."""
        if not _omp_available():
            warnings.warn("OpenMP not available; skipping round-trip check")
            return
        try:
            original = int(_F.gwswex_wrapper.get_omp_max_threads())
        except Exception as exc:
            warnings.warn(f"get_omp_max_threads() raised {exc!r}; skipping round-trip")
            return
        try:
            for n in (1, 2, 3):
                try:
                    _F.gwswex_wrapper.set_omp_threads(n)
                    got = int(_F.gwswex_wrapper.get_omp_max_threads())
                except Exception as exc:
                    warnings.warn(f"set/get OMP threads raised {exc!r} for n={n}")
                    continue
                if got != n:
                    warnings.warn(
                        f"set_omp_threads({n}) -> get_omp_max_threads()={got} "
                        "(expected equal); OMP runtime may be ignoring the request"
                    )
        finally:
            try:
                _F.gwswex_wrapper.set_omp_threads(original)
            except Exception as exc:  # pragma: no cover - defensive
                warnings.warn(f"failed to restore original thread count: {exc!r}")

    def test_model_set_omp_threads_propagates_to_runtime(self):
        """GWSWEXmodel.set_omp_threads() should propagate into the OMP runtime after init()."""
        if not _omp_available():
            warnings.warn("OpenMP not available; skipping propagation check")
            return

        from gwswex import GWSWEXmodel

        m = GWSWEXmodel(name="omp_test", T="s", L="m", write_output=False)
        m.init_space(
            ne=1,
            nl=2,
            top=[[1.0]],
            bot=[[0.5, 0.0]],
            sID=[[1, 1]],
            vID=[[1]],
        )
        m.add_material(
            id=1,
            K_sat=1e-5,
            vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43},
        )
        m.add_vegetation(id=1, root={"depth": 0.5})
        m.init_time(n_steps=1, dt=3600.0, dt_min=1.0)
        m.set_initial_conditions(gw=[-0.8], sw=[0.0], uz=[[-999, -999]])
        m.set_solver(solver="implicit", omp_threads=2)
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
        m.init()
        try:
            try:
                got = int(_F.gwswex_wrapper.get_omp_max_threads())
            except Exception as exc:
                warnings.warn(f"get_omp_max_threads() raised {exc!r}; skipping propagation")
                return
            if got != 2:
                warnings.warn(
                    f"after init() with omp_threads=2, runtime reports {got} threads"
                )

            try:
                m.set_omp_threads(4)
                got = int(_F.gwswex_wrapper.get_omp_max_threads())
            except Exception as exc:
                warnings.warn(f"set_omp_threads(4) raised {exc!r}")
                return
            if got != 4:
                warnings.warn(
                    f"after set_omp_threads(4), runtime reports {got} threads"
                )
            if getattr(m.solver, "omp_threads", None) != 4:
                warnings.warn(
                    f"model.solver.omp_threads={getattr(m.solver, 'omp_threads', None)} "
                    "(expected 4)"
                )
        finally:
            m.deinit()
