! wrapper.f90 -- f2py interface module for the GWSWEX kernel.
!
!! Thin wrappers around `gwswex_kernel` subroutines. Each subroutine carries
!! explicit `!f2py intent` annotations so that f2py generates the correct
!! Python-callable signatures with dimension arguments hidden from Python.
!!
!! Build pipeline: f2py processes this file to generate C glue code, which is
!! compiled and linked with the static `gwswex_core` library (see meson.build).
module gwswex_wrapper
  use gwswex_kernel
  implicit none

contains

  !> Initialise the kernel: allocate arrays, expand material properties, and
  !! store the rooting mask. Must be called before any other wrapper subroutine.
  subroutine init(ne, nl, nmat, solver_type, bnds, sID, &
                  K_sat, theta_s, theta_r, alpha, vg_n, lam, &
                  is_root, ierr)
    !f2py intent(in)  :: ne, nl, nmat, solver_type, bnds, sID
    !f2py intent(in)  :: K_sat, theta_s, theta_r, alpha, vg_n, lam
    !f2py intent(in)  :: is_root
    !f2py intent(out) :: ierr
    integer,  intent(in)  :: ne, nl, nmat, solver_type
    real(8), intent(in)  :: bnds(nl+1, ne)
    integer,  intent(in)  :: sID(nl, ne)
    real(8), intent(in)  :: K_sat(nmat), theta_s(nmat), theta_r(nmat)
    real(8), intent(in)  :: alpha(nmat), vg_n(nmat), lam(nmat)
    integer,  intent(in)  :: is_root(nl, ne)
    integer,  intent(out) :: ierr
    call kernel_init(ne, nl, nmat, solver_type, bnds, sID, K_sat, theta_s, theta_r, &
                     alpha, vg_n, lam, is_root, ierr)
  end subroutine init

  subroutine set_ic(ne, nl, gw, sw, uz, ierr)
    !f2py intent(in)  :: ne, nl, gw, sw, uz
    !f2py intent(out) :: ierr
    integer,  intent(in)  :: ne, nl
    real(8), intent(in)  :: gw(ne), sw(ne), uz(nl, ne)
    integer,  intent(out) :: ierr
    call kernel_set_ic(gw, sw, uz, ierr)
  end subroutine set_ic

  subroutine step(dt, ne, precip, pet, ptt, lat_gw, lat_sw, ierr)
    !f2py intent(in)  :: dt, ne, precip, pet, ptt, lat_gw, lat_sw
    !f2py intent(out) :: ierr
    real(8), intent(in)  :: dt
    integer,  intent(in)  :: ne
    real(8), intent(in)  :: precip(ne), pet(ne), ptt(ne)
    real(8), intent(in)  :: lat_gw(ne), lat_sw(ne)
    integer,  intent(out) :: ierr
    call kernel_step(dt, precip, pet, ptt, lat_gw, lat_sw, ierr)
  end subroutine step

  subroutine get_gw(out_gw, ne)
    !f2py intent(out) :: out_gw
    !f2py intent(in)  :: ne
    integer,  intent(in)  :: ne
    real(8), intent(out) :: out_gw(ne)
    call kernel_get_gw(out_gw)
  end subroutine get_gw

  subroutine get_gwv(out_gwv, ne)
    !f2py intent(out) :: out_gwv
    !f2py intent(in)  :: ne
    integer,  intent(in)  :: ne
    real(8), intent(out) :: out_gwv(ne)
    call kernel_get_gwv(out_gwv)
  end subroutine get_gwv

  subroutine get_sw(out_sw, ne)
    !f2py intent(out) :: out_sw
    !f2py intent(in)  :: ne
    integer,  intent(in)  :: ne
    real(8), intent(out) :: out_sw(ne)
    call kernel_get_sw(out_sw)
  end subroutine get_sw

  subroutine get_uz(out_uz, nl, ne)
    !f2py intent(out) :: out_uz
    !f2py intent(in)  :: nl, ne
    integer,  intent(in)  :: nl, ne
    real(8), intent(out) :: out_uz(nl, ne)
    call kernel_get_uz(out_uz)
  end subroutine get_uz

  subroutine get_theta(out_theta, nl, ne)
    !f2py intent(out) :: out_theta
    !f2py intent(in)  :: nl, ne
    integer,  intent(in)  :: nl, ne
    real(8), intent(out) :: out_theta(nl, ne)
    call kernel_get_theta(out_theta)
  end subroutine get_theta

  !> Retrieve accumulated flux diagnostics for the last completed macro-step.
  subroutine get_accumulators(ne, acc_p, acc_inf, acc_e, acc_t, &
                              acc_r, acc_ro, acc_lg, acc_ls, &
                              acc_dgw, acc_dsw, acc_duz, nsub)
    !f2py intent(in)  :: ne
    !f2py intent(out) :: acc_p, acc_inf, acc_e, acc_t, acc_r, acc_ro, acc_lg, acc_ls
    !f2py intent(out) :: acc_dgw, acc_dsw, acc_duz, nsub
    integer,  intent(in)  :: ne
    real(8), intent(out) :: acc_p(ne), acc_inf(ne), acc_e(ne), acc_t(ne)
    real(8), intent(out) :: acc_r(ne), acc_ro(ne), acc_lg(ne), acc_ls(ne)
    real(8), intent(out) :: acc_dgw(ne), acc_dsw(ne), acc_duz(ne)
    integer,  intent(out) :: nsub(ne)
    call kernel_get_accumulators(acc_p, acc_inf, acc_e, acc_t, acc_r, acc_ro, acc_lg, acc_ls, &
                                 acc_dgw, acc_dsw, acc_duz, nsub)
  end subroutine get_accumulators

  subroutine set_solver_params(courant, dt_min, beta_h, n_trapz, h_min)
    !f2py intent(in) :: courant, dt_min, beta_h, n_trapz, h_min
    real(8), intent(in) :: courant, dt_min, beta_h, h_min
    integer,  intent(in) :: n_trapz
    call kernel_set_solver_params(courant, dt_min, beta_h, n_trapz, h_min)
  end subroutine set_solver_params

  !> Set implicit Picard iteration parameters.
  subroutine set_picard_params(picard_tol, picard_max_iter)
    !f2py intent(in) :: picard_tol, picard_max_iter
    real(8), intent(in) :: picard_tol
    integer,  intent(in) :: picard_max_iter
    call kernel_set_picard_params(picard_tol, picard_max_iter)
  end subroutine set_picard_params

  subroutine get_h(out_h, nl, ne)
    !f2py intent(out) :: out_h
    !f2py intent(in)  :: nl, ne
    integer,  intent(in)  :: nl, ne
    real(8), intent(out) :: out_h(nl, ne)
    call kernel_get_h(out_h)
  end subroutine get_h

  subroutine set_h(h_in, nl, ne)
    !f2py intent(in) :: h_in, nl, ne
    integer,  intent(in) :: nl, ne
    real(8), intent(in) :: h_in(nl, ne)
    call kernel_set_h(h_in)
  end subroutine set_h

  !> Set the number of OpenMP threads for the element-parallel loop.
  !! Effective only when compiled with -fopenmp; otherwise stores the value
  !! in Model%solver%omp_threads with no runtime effect.
  subroutine set_omp_threads(n)
    !f2py intent(in) :: n
    integer, intent(in) :: n
    call kernel_set_omp_threads(n)
  end subroutine set_omp_threads

  !> Returns the OpenMP runtime's current max-thread count.
  !! Equals 1 when built without -fopenmp.
  subroutine get_omp_max_threads(n)
    !f2py intent(out) :: n
    integer, intent(out) :: n
    n = kernel_get_omp_max_threads()
  end subroutine get_omp_max_threads

  !> Returns 1 if the kernel was built with OpenMP (-fopenmp), 0 otherwise.
  subroutine get_omp_available(flag)
    !f2py intent(out) :: flag
    integer, intent(out) :: flag
    flag = kernel_omp_available()
  end subroutine get_omp_available

  !> Set Green-Ampt infiltration and connectivity parameters.
  !! ET stress thresholds are per-vegetation-type; set via set_vegetation.
  subroutine set_model_params(psi_f, F_min, ICratio_min)
    !f2py intent(in) :: psi_f, F_min, ICratio_min
    real(8), intent(in) :: psi_f, F_min, ICratio_min
    call kernel_set_model_params(psi_f, F_min, ICratio_min)
  end subroutine set_model_params

  subroutine set_is_root(nl, ne, rmask)
    !f2py intent(in) :: nl, ne, rmask
    integer,  intent(in) :: nl, ne
    integer,  intent(in) :: rmask(nl, ne)
    call kernel_set_is_root(rmask)
  end subroutine set_is_root

  subroutine set_ic_state(nl, ne, ic_arr, icratio_arr, ne2, f_ga_arr)
    !f2py intent(in) :: nl, ne, ic_arr, icratio_arr, ne2, f_ga_arr
    integer,  intent(in)  :: nl, ne, ne2
    real(8), intent(in)  :: ic_arr(nl, ne), icratio_arr(nl, ne), f_ga_arr(ne2)
    Model%IC      = ic_arr
    Model%ICratio = icratio_arr
    Model%F_GA    = f_ga_arr
  end subroutine set_ic_state

  subroutine get_ic_state(nl, ne, ic_arr, icratio_arr, ne2, f_ga_arr)
    !f2py intent(in)  :: nl, ne, ne2
    !f2py intent(out) :: ic_arr, icratio_arr, f_ga_arr
    integer,  intent(in)  :: nl, ne, ne2
    real(8), intent(out) :: ic_arr(nl, ne), icratio_arr(nl, ne), f_ga_arr(ne2)
    ic_arr      = Model%IC
    icratio_arr = Model%ICratio
    f_ga_arr    = Model%F_GA
  end subroutine get_ic_state

  !> Store vegetation type library and (re)compute per-element root masks.
  !!
  !! For vegetation types with `root_depth > 0`, Fortran recomputes `is_root`
  !! from the root geometry. For types with `root_depth <= 0`, the mask set
  !! during `init` (or by `set_is_root`) is preserved.
  !!
  !! Must be called after `init` and before the first `step`.
  subroutine set_vegetation(nveg, ne, vID, root_depth, &
                             s_star, s_w, s_h, s_e, ierr)
    !f2py intent(in)  :: nveg, ne, vID, root_depth
    !f2py intent(in)  :: s_star, s_w, s_h, s_e
    !f2py intent(out) :: ierr
    integer,  intent(in)  :: nveg, ne
    integer,  intent(in)  :: vID(ne)
    real(8), intent(in)  :: root_depth(nveg)
    real(8), intent(in)  :: s_star(nveg), s_w(nveg), s_h(nveg), s_e(nveg)
    integer,  intent(out) :: ierr

    call kernel_set_vegetation(nveg, ne, vID, root_depth, &
                               s_star, s_w, s_h, s_e, ierr)
  end subroutine set_vegetation

  subroutine deinit()
    call kernel_deinit()
  end subroutine deinit

  !> Switch the active solver mid-simulation, translating state in-place.
  !! `new_solver` is 1 = explicit, 2 = implicit (matches `init`'s
  !! `solver_type`). No-op when the new solver matches the active one.
  subroutine switch_solver(new_solver, ierr)
    !f2py intent(in)  :: new_solver
    !f2py intent(out) :: ierr
    integer, intent(in)  :: new_solver
    integer, intent(out) :: ierr
    call kernel_switch_solver(new_solver, ierr)
  end subroutine switch_solver

  !> Optional warm start of the explicit solver's Green-Ampt /
  !! connectivity persistent state. Pass a negative value for either
  !! argument to leave that field at its cold-start default. Intended
  !! to be called immediately after `switch_solver(SOLVER_EXPLICIT)`
  !! and before the next `step`. See `kernel_warm_start_explicit` for
  !! the exact semantics.
  subroutine warm_start_explicit(icratio_in, f_ga_in, ierr)
    !f2py intent(in)  :: icratio_in, f_ga_in
    !f2py intent(out) :: ierr
    real(8), intent(in)  :: icratio_in, f_ga_in
    integer, intent(out) :: ierr
    call kernel_warm_start_explicit(icratio_in, f_ga_in, ierr)
  end subroutine warm_start_explicit

  !> Warm-start the explicit solver's GA / connectivity persistent state
  !! from the converged Picard `theta_curr` profile. See
  !! `kernel_warm_start_explicit_proxy` for the exact estimator. Should
  !! be called immediately after `switch_solver(SOLVER_EXPLICIT)` and
  !! before the next `step`.
  subroutine warm_start_explicit_proxy(ierr)
    !f2py intent(out) :: ierr
    integer, intent(out) :: ierr
    call kernel_warm_start_explicit_proxy(ierr)
  end subroutine warm_start_explicit_proxy

end module gwswex_wrapper
