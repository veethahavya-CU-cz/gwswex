!> Model kernel: lifecycle management and main stepping loop.
!!
!! The kernel owns the single `gwswex_model` singleton (`Model`) and exposes
!! Fortran subroutines that are called through the f2py wrapper. Responsibilities:
!!   - `kernel_init`: allocate and populate all model arrays from Python inputs.
!!   - `kernel_step`: run one macro-step across all elements (OMP-parallelised).
!!   - `kernel_get_*`: state and accumulator retrieval.
!!   - `kernel_set_*`: parameter and IC setters.
!!   - `kernel_deinit`: deallocate all arrays.
module gwswex_kernel
  use gwswex_constants,         only: dp, EPS, STAT_OK, STAT_ERROR, SOLVER_EXPLICIT, SOLVER_IMPLICIT
  use gwswex_types,             only: gwswex_model
  use gwswex_geometry,          only: update_geometry
  use gwswex_solver_explicit,   only: solve_element_explicit
  use gwswex_solver_implicit,   only: solve_element_implicit
  use gwswex_mass_balance,      only: output_calc
  use gwswex_physics,           only: V_gw, vg_theta, vg_h_inv
  !$ use omp_lib, only: omp_set_num_threads, omp_get_max_threads
  implicit none

  !> Module-level singleton: one gwswex_model instance per process.
  type(gwswex_model), save, target :: Model

contains

  !> Allocate and initialise the model singleton from Python-supplied arrays.
  subroutine kernel_init(ne, nl, nmat, solver_type, bnds, sID, &
                         K_sat, theta_s, theta_r, alpha, vg_n, lam, &
                         is_root, ierr)
    integer,  intent(in)  :: ne, nl, nmat, solver_type
    real(dp), intent(in)  :: bnds(nl+1, ne)
    integer,  intent(in)  :: sID(nl, ne)
    real(dp), intent(in)  :: K_sat(nmat), theta_s(nmat), theta_r(nmat)
    real(dp), intent(in)  :: alpha(nmat), vg_n(nmat), lam(nmat)
    integer,  intent(in)  :: is_root(nl, ne)
    integer,  intent(out) :: ierr
    integer :: ex, lx

    ierr = STAT_OK
    Model%ne  = ne;  Model%nl = nl;  Model%nmat = nmat
    Model%solver%solver_type = solver_type

    ! allocate all arrays
    allocate(Model%bnds(nl+1, ne))
    allocate(Model%K_sat(nl, ne), Model%theta_s(nl, ne), Model%theta_r(nl, ne))
    allocate(Model%alpha(nl, ne), Model%vg_n(nl, ne), Model%vg_m(nl, ne))
    allocate(Model%lambda(nl, ne), Model%Sy(nl, ne))
    allocate(Model%is_root(nl, ne), Model%n_root(ne))
    allocate(Model%GWH_prev(ne), Model%GWH_curr(ne))
    allocate(Model%GWV_prev(ne), Model%GWV_curr(ne))
    allocate(Model%SW_prev(ne), Model%SW_curr(ne))
    allocate(Model%UZ_prev(nl,ne), Model%UZ_curr(nl,ne))
    allocate(Model%theta_prev(nl,ne), Model%theta_curr(nl,ne))
    allocate(Model%IC(nl, ne), Model%ICratio(nl, ne), Model%F_GA(ne))
    allocate(Model%gw_bnd_idx(ne), Model%d_a(nl, ne), Model%ePV(nl, ne))
    allocate(Model%UZ_eq(nl, ne), Model%V_cum(nl+1, ne))
    allocate(Model%K_unsat(nl, ne), Model%tc(nl, ne))
    allocate(Model%precip_rate(ne), Model%pet_rate(ne), Model%ptt_rate(ne))
    allocate(Model%lat_gw_rate(ne), Model%lat_sw_rate(ne))
    allocate(Model%acc_precip(ne), Model%acc_infiltration(ne))
    allocate(Model%acc_evap(ne), Model%acc_transp(ne))
    allocate(Model%acc_recharge(ne), Model%acc_runoff(ne))
    allocate(Model%acc_lat_gw(ne), Model%acc_lat_sw(ne))
    allocate(Model%acc_delta_gw(ne), Model%acc_delta_sw(ne), Model%acc_delta_uz(ne))
    allocate(Model%n_substeps(ne))

    ! Implicit solver state (head arrays) — only when solver_type == SOLVER_IMPLICIT
    if (solver_type == SOLVER_IMPLICIT) then
      allocate(Model%h_prev(nl, ne), Model%h_curr(nl, ne))
      Model%h_prev = 0.0_dp   ! triggers hydrostatic warm-start on first step
      Model%h_curr = 0.0_dp
    end if

    ! copy quasi-static parameters
    Model%bnds = bnds

    ! expand material properties to per-layer arrays for vectorisation
    do ex = 1, ne
      do lx = 1, nl
        Model%K_sat(lx, ex)   = K_sat(sID(lx, ex))
        Model%theta_s(lx, ex) = theta_s(sID(lx, ex))
        Model%theta_r(lx, ex) = theta_r(sID(lx, ex))
        Model%alpha(lx, ex)   = alpha(sID(lx, ex))
        Model%vg_n(lx, ex)    = vg_n(sID(lx, ex))
        Model%vg_m(lx, ex)    = 1.0_dp - 1.0_dp / vg_n(sID(lx, ex))
        Model%lambda(lx, ex)  = lam(sID(lx, ex))
        Model%Sy(lx, ex)      = theta_s(sID(lx, ex)) - theta_r(sID(lx, ex))
      end do
    end do
    Model%is_root     = is_root
    do ex = 1, ne
      Model%n_root(ex) = sum(Model%is_root(:, ex))
    end do

    ! initialise all accumulators to zero
    Model%acc_precip(:)       = 0.0_dp
    Model%acc_infiltration(:) = 0.0_dp
    Model%acc_evap(:)         = 0.0_dp
    Model%acc_transp(:)       = 0.0_dp
    Model%acc_recharge(:)     = 0.0_dp
    Model%acc_runoff(:)       = 0.0_dp
    Model%acc_lat_gw(:)       = 0.0_dp
    Model%acc_lat_sw(:)       = 0.0_dp
    Model%acc_delta_gw(:)     = 0.0_dp
    Model%acc_delta_sw(:)     = 0.0_dp
    Model%acc_delta_uz(:)     = 0.0_dp
    Model%n_substeps(:)       = 0
  end subroutine kernel_init

  !> Store the vegetation type library and (re)compute per-element root masks.
  !!
  !! Populates `Model%veg`, `Model%vID`, `Model%nveg`, and recalculates
  !! `Model%is_root` and `Model%n_root` from the per-type rooting depth and
  !! the column geometry. Vegetation types with `root_depth <= 0` leave the
  !! `is_root` mask of their elements unchanged (allowing Python-supplied
  !! masks set during `kernel_init` to be preserved).
  !!
  !! @param nveg         Number of vegetation types (must be >= 1).
  !! @param ne_in        Number of spatial elements (must equal `Model%ne`).
  !! @param vID_arr      Per-element vegetation type IDs, 1-based, shape (ne_in).
  !! @param root_depth   Maximum rooting depth per type [L]; <= 0 leaves mask intact.
  !! @param s_star       ET stress param --- incipient stomatal closure [-], shape (nveg).
  !! @param s_w          ET stress param --- wilting point [-], shape (nveg).
  !! @param s_h          ET stress param --- hygroscopic point [-], shape (nveg).
  !! @param s_e          ET stress param --- capillary-continuity threshold [-], shape (nveg).
  !! @param ierr         Return status: 0 on success, -1 on invalid arguments.
  subroutine kernel_set_vegetation(nveg, ne_in, vID_arr, root_depth, &
                                   s_star, s_w, s_h, s_e, ierr)
    use gwswex_constants, only: EPS
    integer,  intent(in)  :: nveg, ne_in
    integer,  intent(in)  :: vID_arr(ne_in)
    real(dp), intent(in)  :: root_depth(nveg)
    real(dp), intent(in)  :: s_star(nveg), s_w(nveg), s_h(nveg), s_e(nveg)
    integer,  intent(out) :: ierr

    integer :: vt, ex
    ierr = STAT_OK

    if (nveg < 1 .or. ne_in /= Model%ne) then
      ierr = STAT_ERROR
      return
    end if

    ! --- Allocate vegetation library ---
    if (allocated(Model%veg)) deallocate(Model%veg)
    if (allocated(Model%vID)) deallocate(Model%vID)
    allocate(Model%veg(nveg), Model%vID(Model%ne))
    Model%nveg = nveg
    Model%vID  = vID_arr

    ! --- Store per-type parameters ---
    do vt = 1, nveg
      Model%veg(vt)%root_depth = root_depth(vt)
      Model%veg(vt)%s_star     = s_star(vt)
      Model%veg(vt)%s_w        = s_w(vt)
      Model%veg(vt)%s_h        = s_h(vt)
      Model%veg(vt)%s_e        = s_e(vt)
    end do

    ! --- Recompute root mask for elements whose vegetation type carries depth ---
    do ex = 1, Model%ne
      vt = vID_arr(ex)
      if (vt < 1 .or. vt > nveg) cycle
      if (root_depth(vt) <= 0.0_dp) cycle   ! leave Python-supplied mask intact
      call compute_root_mask( &
           Model%nl, Model%bnds(1, ex), Model%bnds(:, ex), &
           root_depth(vt), Model%is_root(:, ex))
      Model%n_root(ex) = sum(Model%is_root(:, ex))
    end do

  end subroutine kernel_set_vegetation

  !> Compute `is_root` for a single element column from the maximum rooting depth.
  !!
  !! A layer is flagged as rooted when its midpoint depth below the surface lies
  !! within `[0, root_depth]`. Transpiration demand is partitioned uniformly
  !! across the rooted layers at solve time using `1 / n_root(ex)` weighting,
  !! so no per-layer density profile is stored.
  !!
  !! @param nl            Number of layers.
  !! @param surface       Surface elevation (= bnds(1)) [L].
  !! @param bnds          Layer boundary elevations [L], shape (nl+1), top first.
  !! @param root_depth    Maximum rooting depth below surface [L].
  !! @param is_root       Output: root layer mask, shape (nl).
  subroutine compute_root_mask(nl, surface, bnds, root_depth, is_root)
    integer,  intent(in)  :: nl
    real(dp), intent(in)  :: surface
    real(dp), intent(in)  :: bnds(nl+1)
    real(dp), intent(in)  :: root_depth
    integer,  intent(out) :: is_root(nl)

    integer  :: lx
    real(dp) :: z_mid, depth

    is_root = 0
    do lx = 1, nl
      z_mid = 0.5_dp * (bnds(lx) + bnds(lx+1))
      depth = surface - z_mid
      if (depth >= 0.0_dp .and. depth <= root_depth) is_root(lx) = 1
    end do
  end subroutine compute_root_mask

  !> Set initial conditions (GW head, SW depth, UZ storage per layer).
  subroutine kernel_set_ic(gw, sw, uz, ierr)
    real(dp), intent(in)  :: gw(:), sw(:), uz(:,:)
    integer,  intent(out) :: ierr
    integer :: ex
    ierr = STAT_OK
    Model%GWH_prev = gw;  Model%GWH_curr = gw
    Model%SW_prev  = sw;  Model%SW_curr  = sw
    Model%UZ_prev  = uz;  Model%UZ_curr  = uz
    ! compute initial GWV from GWH
    do ex = 1, Model%ne
      Model%GWV_prev(ex) = V_gw(gw(ex), Model%bnds(:, ex), Model%Sy(:, ex), Model%nl)
      Model%GWV_curr(ex) = Model%GWV_prev(ex)
    end do
    ! Explicit-solver-only: persistent infiltration parameters and initial geometry.
    ! The implicit solver derives all geometry from h_prev at each Picard step
    ! and does not use IC, ICratio, F_GA, or the d_a/ePV/UZ_eq/V_cum arrays.
    if (Model%solver%solver_type == SOLVER_EXPLICIT) then
      Model%IC(:, :)      = 0.0_dp
      Model%ICratio(:, :) = Model%params%ICratio_min
      Model%F_GA(:)       = Model%params%F_min
      ! initial geometry; resolve UZ sentinel (-999) to hydrostatic VG profile
      !$omp parallel do schedule(static)
      do ex = 1, Model%ne
        call update_geometry(Model, ex)
        ! If UZ was requested at equilibrium (sentinel -999), overwrite with the
        ! hydrostatic VG profile computed by update_geometry
        if (Model%UZ_prev(1, ex) < -900.0_dp) then
          Model%UZ_prev(:, ex) = Model%UZ_eq(:, ex)
          Model%UZ_curr(:, ex) = Model%UZ_eq(:, ex)
        end if
      end do
      !$omp end parallel do
    else if (Model%solver%solver_type == SOLVER_IMPLICIT) then
      ! Resolve UZ sentinel (-999) to hydrostatic VG profile for implicit solver.
      ! The implicit solver derives state from h, but UZ_prev/UZ_curr must be
      ! physically consistent at t=0 so that storage diagnostics and the
      ! physical water formula give correct initial values.
      !$omp parallel do schedule(static) private(ex)
      do ex = 1, Model%ne
        if (Model%UZ_prev(1, ex) < -900.0_dp) then
          call resolve_uz_sentinel_implicit(ex)
        end if
      end do
      !$omp end parallel do
    end if

    ! Populate theta_prev/theta_curr from the hydrostatic head profile so
    ! that get_state()/get_theta() return a physically consistent moisture
    ! field immediately after init(), before any solver step has run.  On
    ! the implicit sentinel path these arrays have already been assigned
    ! above; re-running the assignment there is idempotent.
    !$omp parallel do schedule(static) private(ex)
    do ex = 1, Model%ne
      call init_theta_from_hydrostatic(ex)
    end do
    !$omp end parallel do
  end subroutine kernel_set_ic

  !> Compute theta_prev/theta_curr for element `ex` from the hydrostatic head
  !! profile implied by GWH_prev.  Saturated layers get theta_s; unsaturated
  !! layers get theta(h) from the van Genuchten retention curve.  UZ state is
  !! not modified here.
  subroutine init_theta_from_hydrostatic(ex)
    integer, intent(in) :: ex
    integer  :: l
    real(dp) :: z_mid, h_l, theta_l
    do l = 1, Model%nl
      z_mid = 0.5_dp * (Model%bnds(l, ex) + Model%bnds(l+1, ex))
      h_l   = Model%GWH_prev(ex) - z_mid
      if (h_l >= 0.0_dp) then
        theta_l = Model%theta_s(l, ex)
      else
        theta_l = vg_theta(h_l, Model%theta_r(l, ex), Model%theta_s(l, ex), &
                           Model%alpha(l, ex), Model%vg_n(l, ex), Model%vg_m(l, ex))
      end if
      Model%theta_prev(l, ex) = theta_l
      Model%theta_curr(l, ex) = theta_l
    end do
  end subroutine init_theta_from_hydrostatic

  !> Resolve the -999 UZ sentinel for an implicit-solver element.
  !!
  !! Computes the hydrostatic head profile h = GWH - z_mid for each layer,
  !! evaluates θ(h) from the van Genuchten retention curve, and sets
  !! UZ = θ × dz for unsaturated layers and UZ = θ_r × dz for saturated
  !! layers (matching the h_to_state convention in solver_implicit).
  !! Also initialises h_prev to the hydrostatic profile.
  subroutine resolve_uz_sentinel_implicit(ex)
    integer, intent(in) :: ex
    integer  :: l, nl
    real(dp) :: z_mid, dz_l, h_l, theta_l

    nl = Model%nl
    do l = 1, nl
      z_mid = 0.5_dp * (Model%bnds(l, ex) + Model%bnds(l+1, ex))
      dz_l  = Model%bnds(l, ex) - Model%bnds(l+1, ex)
      h_l   = Model%GWH_prev(ex) - z_mid   ! hydrostatic: h = z_wt - z

      ! Store hydrostatic head for warm-start
      Model%h_prev(l, ex) = h_l
      Model%h_curr(l, ex) = h_l

      if (h_l >= 0.0_dp) then
        ! Saturated: UZ stores 0 (drainable part in GWV; residual theta_r
        ! reported diagnostically as theta_r * max(GWH, 0) — see h_to_state
        ! in solver_implicit for the matching convention).
        Model%UZ_prev(l, ex) = 0.0_dp
        Model%UZ_curr(l, ex) = 0.0_dp
        Model%theta_prev(l, ex) = Model%theta_s(l, ex)
        Model%theta_curr(l, ex) = Model%theta_s(l, ex)
      else
        ! Unsaturated: full θ(h) from van Genuchten
        theta_l = vg_theta(h_l, Model%theta_r(l, ex), Model%theta_s(l, ex), &
                           Model%alpha(l, ex), Model%vg_n(l, ex), Model%vg_m(l, ex))
        Model%UZ_prev(l, ex) = theta_l * dz_l
        Model%UZ_curr(l, ex) = theta_l * dz_l
        Model%theta_prev(l, ex) = theta_l
        Model%theta_curr(l, ex) = theta_l
      end if
    end do
  end subroutine resolve_uz_sentinel_implicit

  !> Run one macro-step across all elements (OMP-parallelised over elements).
  !! Zeros all flux accumulators, calls the appropriate solver per element,
  !! then advances the prev/curr state arrays.
  subroutine kernel_step(dt, precip, pet, ptt, lat_gw, lat_sw, ierr)
    real(dp), intent(in)  :: dt
    real(dp), intent(in)  :: precip(:), pet(:), ptt(:)
    real(dp), intent(in)  :: lat_gw(:), lat_sw(:)
    integer,  intent(out) :: ierr
    integer :: ex

    ierr = STAT_OK

    ! set forcings
    Model%precip_rate = precip
    Model%pet_rate    = pet
    Model%ptt_rate    = ptt
    Model%lat_gw_rate = lat_gw
    Model%lat_sw_rate = lat_sw

    ! zero accumulators
    Model%acc_precip(:)       = 0.0_dp
    Model%acc_infiltration(:) = 0.0_dp
    Model%acc_evap(:)         = 0.0_dp
    Model%acc_transp(:)       = 0.0_dp
    Model%acc_recharge(:)     = 0.0_dp
    Model%acc_runoff(:)       = 0.0_dp
    Model%acc_lat_gw(:)       = 0.0_dp
    Model%acc_lat_sw(:)       = 0.0_dp
    Model%acc_delta_gw(:)     = 0.0_dp
    Model%acc_delta_sw(:)     = 0.0_dp
    Model%acc_delta_uz(:)     = 0.0_dp
    Model%n_substeps(:)       = 0

    if (Model%solver%solver_type == SOLVER_EXPLICIT) then
        !$omp parallel do schedule(static)
        do ex = 1, Model%ne
            call solve_element_explicit(Model, ex, dt)
        end do
        !$omp end parallel do
    else if (Model%solver%solver_type == SOLVER_IMPLICIT) then
        !$omp parallel do schedule(static)
        do ex = 1, Model%ne
            call solve_element_implicit(Model, ex, dt)
        end do
        !$omp end parallel do
    end if

    ! Output / accumulator accounting (solver-independent)
    !$omp parallel do schedule(static)
    do ex = 1, Model%ne
      call output_calc(Model, ex, dt)
    end do
    !$omp end parallel do

    ! Compute theta_curr from UZ for the explicit solver
    ! (the implicit solver already sets theta_curr via h_to_state).
    !
    ! Convention: UZ = theta * d_a  (total volumetric water × active thickness,
    ! including the residual theta_r fraction).  This is consistent with
    ! h_to_state ("UZ = theta * dz") and with all explicit-process code that
    ! uses UZ - theta_r*d_a as the drainable volume above residual.
    ! Therefore theta = UZ / d_a, clamped at theta_s.
    if (Model%solver%solver_type == SOLVER_EXPLICIT) then
      where (Model%d_a > EPS)
        Model%theta_curr = min(Model%UZ_curr / Model%d_a, Model%theta_s)
      elsewhere
        Model%theta_curr = Model%theta_s   ! saturated (d_a = 0)
      end where
    end if

    ! copy current → previous for next step
    Model%GWH_prev   = Model%GWH_curr
    Model%GWV_prev   = Model%GWV_curr
    Model%SW_prev    = Model%SW_curr
    Model%UZ_prev    = Model%UZ_curr
    Model%theta_prev = Model%theta_curr
    if (Model%solver%solver_type == SOLVER_IMPLICIT) then
      Model%h_prev = Model%h_curr
    end if
  end subroutine kernel_step

  !> Retrieve current GW table elevations [L] (previous macro-step end).
  subroutine kernel_get_gw(out_gw)
    real(dp), intent(out) :: out_gw(:)
    out_gw = Model%GWH_prev
  end subroutine kernel_get_gw

  subroutine kernel_get_gwv(out_gwv)
    real(dp), intent(out) :: out_gwv(:)
    out_gwv = Model%GWV_prev
  end subroutine kernel_get_gwv

  subroutine kernel_get_sw(out_sw)
    real(dp), intent(out) :: out_sw(:)
    out_sw = Model%SW_prev
  end subroutine kernel_get_sw

  subroutine kernel_get_uz(out_uz)
    real(dp), intent(out) :: out_uz(:,:)
    out_uz = Model%UZ_prev
  end subroutine kernel_get_uz

  subroutine kernel_get_theta(out_theta)
    real(dp), intent(out) :: out_theta(:,:)
    out_theta = Model%theta_prev
  end subroutine kernel_get_theta

  subroutine kernel_get_accumulators(acc_p, acc_inf, acc_e, acc_t, &
                                     acc_r, acc_ro, acc_lg, acc_ls, &
                                     acc_dgw, acc_dsw, acc_duz, nsub)
    real(dp), intent(out) :: acc_p(:), acc_inf(:), acc_e(:), acc_t(:)
    real(dp), intent(out) :: acc_r(:), acc_ro(:), acc_lg(:), acc_ls(:)
    real(dp), intent(out) :: acc_dgw(:), acc_dsw(:), acc_duz(:)
    integer,  intent(out) :: nsub(:)
    acc_p   = Model%acc_precip
    acc_inf = Model%acc_infiltration
    acc_e   = Model%acc_evap
    acc_t   = Model%acc_transp
    acc_r   = Model%acc_recharge
    acc_ro  = Model%acc_runoff
    acc_lg  = Model%acc_lat_gw
    acc_ls  = Model%acc_lat_sw
    acc_dgw = Model%acc_delta_gw
    acc_dsw = Model%acc_delta_sw
    acc_duz = Model%acc_delta_uz
    nsub    = Model%n_substeps
  end subroutine kernel_get_accumulators

  !> Update root mask in-place (used for time-varying vegetation root depth).
  subroutine kernel_set_is_root(rmask)
    integer, intent(in) :: rmask(:,:)
    integer :: ex
    Model%is_root = rmask
    do ex = 1, Model%ne
      Model%n_root(ex) = sum(Model%is_root(:, ex))
    end do
  end subroutine kernel_set_is_root

  !> Set the number of OpenMP threads for the element-parallel loop.
  !! If compiled without OpenMP (-fopenmp absent), this is a no-op (the OMP
  !! sentinel `!$` causes the call to be compiled only when OMP is active).
  subroutine kernel_set_omp_threads(n)
    integer, intent(in) :: n
    Model%solver%omp_threads = n
    !$ call omp_set_num_threads(n)
  end subroutine kernel_set_omp_threads

  !> Return the OpenMP runtime's current max-thread count (i.e. the value
  !! that would size an upcoming parallel region). Equals 1 when built
  !! without -fopenmp.
  function kernel_get_omp_max_threads() result(n)
    integer :: n
    n = 1
    !$ n = omp_get_max_threads()
  end function kernel_get_omp_max_threads

  !> Returns 1 when the kernel was built with -fopenmp (the `_OPENMP`
  !! preprocessor sentinel is in scope), 0 otherwise. Lets Python-side
  !! callers detect, without ambiguity, whether OMP is actually wired in,
  !! independent of the runtime thread count.
  function kernel_omp_available() result(flag)
    integer :: flag
    flag = 0
    !$ flag = 1
  end function kernel_omp_available

  !> Set solver parameters (shared by both solvers).
  subroutine kernel_set_solver_params(courant, dt_min, beta_h, n_trapz, h_min_in)
    real(dp), intent(in) :: courant, dt_min, beta_h, h_min_in
    integer,  intent(in) :: n_trapz
    Model%solver%courant_number = courant
    Model%solver%dt_min         = dt_min
    Model%solver%beta_h         = beta_h
    Model%solver%n_trapz        = n_trapz
    Model%solver%h_min          = h_min_in
  end subroutine kernel_set_solver_params

  !> Set implicit solver Picard iteration parameters.
  subroutine kernel_set_picard_params(picard_tol, picard_max_iter)
    real(dp), intent(in) :: picard_tol
    integer,  intent(in) :: picard_max_iter
    Model%solver%picard_tol      = picard_tol
    Model%solver%picard_max_iter = picard_max_iter
  end subroutine kernel_set_picard_params

  subroutine kernel_get_h(out_h)
    real(dp), intent(out) :: out_h(:,:)
    out_h = Model%h_prev
  end subroutine kernel_get_h

  subroutine kernel_set_h(h_in)
    real(dp), intent(in) :: h_in(:,:)
    Model%h_prev = h_in
    Model%h_curr = h_in
  end subroutine kernel_set_h

  !> Set Green-Ampt infiltration and connectivity parameters.
  !! ET stress thresholds are per-vegetation-type; set via kernel_set_vegetation.
  subroutine kernel_set_model_params(psi_f, F_min, ICratio_min)
    real(dp), intent(in) :: psi_f, F_min, ICratio_min
    Model%params%psi_f       = psi_f
    Model%params%F_min       = F_min
    Model%params%ICratio_min = ICratio_min
  end subroutine kernel_set_model_params

  !> Switch the active solver mid-simulation, translating state into the
  !! representation expected by the new solver.
  !!
  !! Direction-specific work performed (state arrays GWH/GWV/SW/UZ/theta
  !! are already up-to-date when the previous solver returned; the kernel
  !! advances them in-place each macro-step):
  !!
  !!   * EXPLICIT → IMPLICIT
  !!     - Allocate `h_prev` / `h_curr` if not already present.
  !!     - Build a physically consistent matric-head profile per layer:
  !!       hydrostatic positive head for layers below the water table
  !!       (z_mid < GWH_curr); van-Genuchten inverse `vg_h_inv(θ)` for
  !!       unsaturated layers, so that θ implied by `h` matches the
  !!       θ produced by the explicit cascade.  This avoids the
  !!       spurious storage transient that a naïve hydrostatic warm
  !!       start would introduce on the first implicit Picard solve.
  !!
  !!   * IMPLICIT → EXPLICIT
  !!     - Reset Green-Ampt / connectivity persistent state to its
  !!       init() defaults so that the explicit cascade does not inherit
  !!       a stale infiltration front from the implicit run.
  !!     - Refresh per-element geometry (`update_geometry`) to populate
  !!       `gw_bnd_idx`, `d_a`, `ePV`, `UZ_eq`, `V_cum`, which the
  !!       implicit solver does not maintain.
  !!     - Re-derive `UZ_curr = θ * d_a` for active layers so that the
  !!       explicit storage convention (UZ measured against the active
  !!       thickness `d_a`, not the geometric layer thickness `dz`) is
  !!       preserved.
  !!
  !! `h_prev` / `h_curr` are kept allocated even after a switch to the
  !! explicit solver so that subsequent switches back to implicit do not
  !! pay the allocation cost or lose the warm-start hint.
  subroutine kernel_switch_solver(new_solver, ierr)
    integer, intent(in)  :: new_solver
    integer, intent(out) :: ierr
    integer  :: ex, l, nl
    real(dp) :: z_mid, h_l

    ierr = STAT_OK
    if (new_solver /= SOLVER_EXPLICIT .and. new_solver /= SOLVER_IMPLICIT) then
      ierr = STAT_ERROR
      return
    end if
    if (Model%ne <= 0 .or. Model%nl <= 0) then
      ! Kernel never initialised: nothing to translate.
      ierr = STAT_ERROR
      return
    end if
    if (new_solver == Model%solver%solver_type) return

    nl = Model%nl

    if (new_solver == SOLVER_IMPLICIT) then
      if (.not. allocated(Model%h_prev)) allocate(Model%h_prev(nl, Model%ne))
      if (.not. allocated(Model%h_curr)) allocate(Model%h_curr(nl, Model%ne))

      !$omp parallel do schedule(static) private(ex, l, z_mid, h_l)
      do ex = 1, Model%ne
        do l = 1, nl
          z_mid = 0.5_dp * (Model%bnds(l, ex) + Model%bnds(l+1, ex))
          if (z_mid < Model%GWH_curr(ex)) then
            ! Saturated: hydrostatic positive head from the water table
            h_l = Model%GWH_curr(ex) - z_mid
            Model%theta_curr(l, ex) = Model%theta_s(l, ex)
            Model%theta_prev(l, ex) = Model%theta_s(l, ex)
          else
            ! Unsaturated: invert θ(h) so that warm-start θ matches the
            ! explicit-solver θ_curr (no spurious storage transient).
            h_l = vg_h_inv(Model%theta_curr(l, ex), &
                           Model%theta_r(l, ex), Model%theta_s(l, ex), &
                           Model%alpha(l, ex),   Model%vg_n(l, ex), &
                           Model%vg_m(l, ex),    Model%solver%h_min)
          end if
          Model%h_prev(l, ex) = h_l
          Model%h_curr(l, ex) = h_l
        end do
      end do
      !$omp end parallel do

    else  ! new_solver == SOLVER_EXPLICIT
      ! Reset GA/connectivity persistent state to init() defaults
      Model%IC      = 0.0_dp
      Model%ICratio = Model%params%ICratio_min
      Model%F_GA    = Model%params%F_min

      ! Refresh geometry and re-derive UZ in the explicit convention
      !
      ! NB: the implicit solver may legitimately leave a layer fully
      ! saturated (theta = theta_s) above the water table — e.g. capillary-
      ! /tension-saturated zones produced by the Picard solve.  When such a
      ! profile is fed to the explicit cascade, the precipitation CFL
      ! (eval_cfl) collapses to (max(ePV-UZ, EPS)/precip) ≈ EPS/precip,
      ! flooring sub-steps at dt_min and producing an effective hang.
      !
      ! Translation rule: clamp UZ_curr at most to (1 - SAT_HEADROOM) * ePV
      ! per layer, and route the trimmed water to surface storage to keep
      ! the column mass-balance closed (the alternative — raising the WT —
      ! would silently reorganise the storage partition).  The headroom is
      ! small enough to be physically negligible relative to ePV but large
      ! enough that the precip CFL stays well above dt_min.
      block
        real(dp), parameter :: SAT_HEADROOM = 1.0e-6_dp
        real(dp) :: cap, trim_amount
      !$omp parallel do schedule(static) private(ex, l, cap, trim_amount)
      do ex = 1, Model%ne
        call update_geometry(Model, ex)
        do l = 1, nl
          if (Model%gw_bnd_idx(ex) > 0 .and. l <= Model%gw_bnd_idx(ex)) then
            Model%UZ_curr(l, ex) = Model%theta_curr(l, ex) * Model%d_a(l, ex)
            cap = (1.0_dp - SAT_HEADROOM) * Model%ePV(l, ex)
            if (Model%UZ_curr(l, ex) > cap) then
              trim_amount = Model%UZ_curr(l, ex) - cap
              Model%UZ_curr(l, ex) = cap
              !$omp atomic
              Model%SW_curr(ex) = Model%SW_curr(ex) + trim_amount
            end if
          else
            Model%UZ_curr(l, ex) = 0.0_dp
          end if
        end do
      end do
      !$omp end parallel do
      end block
      Model%UZ_prev = Model%UZ_curr
    end if

    Model%solver%solver_type = new_solver
  end subroutine kernel_switch_solver

  !> Override the explicit-solver Green-Ampt / connectivity persistent
  !! state with user-supplied uniform warm-start values.  Intended to be
  !! called immediately after `kernel_switch_solver(SOLVER_EXPLICIT, ...)`
  !! so that the cascade does not have to re-discover its connectivity
  !! and infiltration history from cold defaults over the first few
  !! macro-steps.
  !!
  !! Sentinel: a negative value of either argument means "leave that
  !! field at its current (cold-start default) value".  Both arguments
  !! are uniform across layers and elements; for a heterogeneous warm
  !! start use the existing `set_ic_state` after computing IC, ICratio,
  !! and F_GA externally.
  !!
  !! The per-layer wetting-front depth IC is derived from ICratio as
  !!   IC(l, ex) = ICratio_in * d_a(l, ex)
  !! to match the inverse of the cascade's own update rule
  !!   ICratio = max(IC / d_a, ICratio_min).
  !! d_a is taken from the geometry refreshed inside `kernel_switch_solver`
  !! and is valid as long as no `step` has run between the two calls.
  !!
  !! ICratio_in is clamped from below by `params%ICratio_min` and from
  !! above by 1.0; F_GA_in is clamped from below by `params%F_min`.
  !! ierr is `STAT_ERROR` if the kernel is not initialised or if the
  !! active solver is not explicit at call time.
  subroutine kernel_warm_start_explicit(ICratio_in, F_GA_in, ierr)
    real(dp), intent(in)  :: ICratio_in, F_GA_in
    integer,  intent(out) :: ierr
    real(dp) :: r_clamped, f_clamped
    integer  :: ex, l

    ierr = STAT_OK
    if (Model%ne <= 0 .or. Model%nl <= 0) then
      ierr = STAT_ERROR
      return
    end if
    if (Model%solver%solver_type /= SOLVER_EXPLICIT) then
      ierr = STAT_ERROR
      return
    end if

    if (ICratio_in >= 0.0_dp) then
      r_clamped = min(max(ICratio_in, Model%params%ICratio_min), 1.0_dp)
      !$omp parallel do schedule(static) private(ex, l)
      do ex = 1, Model%ne
        do l = 1, Model%nl
          Model%ICratio(l, ex) = r_clamped
          Model%IC(l, ex)      = r_clamped * Model%d_a(l, ex)
        end do
      end do
      !$omp end parallel do
    end if

    if (F_GA_in >= 0.0_dp) then
      f_clamped = max(F_GA_in, Model%params%F_min)
      Model%F_GA(:) = f_clamped
    end if
  end subroutine kernel_warm_start_explicit

  !> Warm-start the explicit solver's GA / connectivity persistent state
  !! from the live unsaturated profile (`theta_curr`).  Intended to be
  !! called immediately after `kernel_switch_solver(SOLVER_EXPLICIT, ...)`
  !! when the user wants a physics-aware seed rather than the cold
  !! defaults.  Operates per layer / element from the converged Picard
  !! profile that the imp\u2192exp branch has just translated.
  !!
  !! Per-layer wetting-front depth and connectivity:
  !!   Se(l, ex)        = clip( (theta_curr - theta_r) / (theta_s - theta_r), 0, 1 )
  !!   IC(l, ex)        = Se(l, ex) * d_a(l, ex)
  !!   ICratio(l, ex)   = max( Se(l, ex), ICratio_min )
  !! For layers below the water table (`l > gw_bnd_idx`) the explicit
  !! convention sets `UZ = 0`, so the GA tracker is also zeroed there
  !! (`IC = 0`, `ICratio = ICratio_min`).
  !!
  !! Per-element Green-Ampt cumulative infiltration:
  !!   F_GA(ex) = clip( sum_l (theta_curr - theta_r) * d_a, F_min, F_cap )
  !! with the cap `F_cap = 5 * psi_f` keeping the next-step infiltration
  !! capacity `f_cap = K_sat * (1 + psi_f * dtheta / F_GA)` from
  !! collapsing on a profile that the implicit solver has driven near
  !! saturation.  The proxy conflates infiltration with capillary
  !! upflow from a shallow water table, which biases F_GA upward (and
  !! therefore the next-step infiltration capacity slightly downward) \u2014
  !! a deliberately safe direction.
  !!
  !! ierr is `STAT_ERROR` if the kernel is not initialised or if the
  !! active solver is not explicit at call time.
  subroutine kernel_warm_start_explicit_proxy(ierr)
    integer, intent(out) :: ierr
    integer  :: ex, l
    real(dp) :: Se, dtheta_l, F_acc, F_cap, denom

    ierr = STAT_OK
    if (Model%ne <= 0 .or. Model%nl <= 0) then
      ierr = STAT_ERROR
      return
    end if
    if (Model%solver%solver_type /= SOLVER_EXPLICIT) then
      ierr = STAT_ERROR
      return
    end if

    F_cap = 5.0_dp * Model%params%psi_f

    !$omp parallel do schedule(static) private(ex, l, Se, dtheta_l, F_acc, denom)
    do ex = 1, Model%ne
      F_acc = 0.0_dp
      do l = 1, Model%nl
        if (Model%gw_bnd_idx(ex) > 0 .and. l <= Model%gw_bnd_idx(ex)) then
          denom = max(Model%theta_s(l, ex) - Model%theta_r(l, ex), EPS)
          Se    = (Model%theta_curr(l, ex) - Model%theta_r(l, ex)) / denom
          Se    = min(max(Se, 0.0_dp), 1.0_dp)
          Model%IC(l, ex)      = Se * Model%d_a(l, ex)
          Model%ICratio(l, ex) = max(Se, Model%params%ICratio_min)
          dtheta_l = max(Model%theta_curr(l, ex) - Model%theta_r(l, ex), 0.0_dp)
          F_acc    = F_acc + dtheta_l * Model%d_a(l, ex)
        else
          Model%IC(l, ex)      = 0.0_dp
          Model%ICratio(l, ex) = Model%params%ICratio_min
        end if
      end do
      Model%F_GA(ex) = min(max(F_acc, Model%params%F_min), F_cap)
    end do
    !$omp end parallel do
  end subroutine kernel_warm_start_explicit_proxy

  !> Deallocate all model arrays and reset the singleton.
  subroutine kernel_deinit()
    if (allocated(Model%bnds))       deallocate(Model%bnds)
    if (allocated(Model%GWH_prev))   deallocate(Model%GWH_prev)
    if (allocated(Model%GWH_curr))   deallocate(Model%GWH_curr)
    if (allocated(Model%GWV_prev))   deallocate(Model%GWV_prev)
    if (allocated(Model%GWV_curr))   deallocate(Model%GWV_curr)
    if (allocated(Model%SW_prev))    deallocate(Model%SW_prev)
    if (allocated(Model%SW_curr))    deallocate(Model%SW_curr)
    if (allocated(Model%UZ_prev))    deallocate(Model%UZ_prev)
    if (allocated(Model%UZ_curr))    deallocate(Model%UZ_curr)
    if (allocated(Model%theta_prev)) deallocate(Model%theta_prev)
    if (allocated(Model%theta_curr)) deallocate(Model%theta_curr)
    if (allocated(Model%K_sat))      deallocate(Model%K_sat)
    if (allocated(Model%theta_s))    deallocate(Model%theta_s)
    if (allocated(Model%theta_r))    deallocate(Model%theta_r)
    if (allocated(Model%alpha))      deallocate(Model%alpha)
    if (allocated(Model%vg_n))       deallocate(Model%vg_n)
    if (allocated(Model%vg_m))       deallocate(Model%vg_m)
    if (allocated(Model%lambda))     deallocate(Model%lambda)
    if (allocated(Model%Sy))         deallocate(Model%Sy)
    if (allocated(Model%is_root))    deallocate(Model%is_root)
    if (allocated(Model%n_root))     deallocate(Model%n_root)
    if (allocated(Model%IC))         deallocate(Model%IC)
    if (allocated(Model%ICratio))    deallocate(Model%ICratio)
    if (allocated(Model%F_GA))       deallocate(Model%F_GA)
    if (allocated(Model%gw_bnd_idx)) deallocate(Model%gw_bnd_idx)
    if (allocated(Model%d_a))        deallocate(Model%d_a)
    if (allocated(Model%ePV))        deallocate(Model%ePV)
    if (allocated(Model%UZ_eq))      deallocate(Model%UZ_eq)
    if (allocated(Model%V_cum))      deallocate(Model%V_cum)
    if (allocated(Model%K_unsat))    deallocate(Model%K_unsat)
    if (allocated(Model%tc))         deallocate(Model%tc)
    if (allocated(Model%precip_rate))deallocate(Model%precip_rate)
    if (allocated(Model%pet_rate))   deallocate(Model%pet_rate)
    if (allocated(Model%ptt_rate))   deallocate(Model%ptt_rate)
    if (allocated(Model%lat_gw_rate))deallocate(Model%lat_gw_rate)
    if (allocated(Model%lat_sw_rate))deallocate(Model%lat_sw_rate)
    if (allocated(Model%acc_precip))      deallocate(Model%acc_precip)
    if (allocated(Model%acc_infiltration))deallocate(Model%acc_infiltration)
    if (allocated(Model%acc_evap))        deallocate(Model%acc_evap)
    if (allocated(Model%acc_transp))      deallocate(Model%acc_transp)
    if (allocated(Model%acc_recharge))    deallocate(Model%acc_recharge)
    if (allocated(Model%acc_runoff))      deallocate(Model%acc_runoff)
    if (allocated(Model%acc_lat_gw))      deallocate(Model%acc_lat_gw)
    if (allocated(Model%acc_lat_sw))      deallocate(Model%acc_lat_sw)
    if (allocated(Model%acc_delta_gw))    deallocate(Model%acc_delta_gw)
    if (allocated(Model%acc_delta_sw))    deallocate(Model%acc_delta_sw)
    if (allocated(Model%acc_delta_uz))    deallocate(Model%acc_delta_uz)
    if (allocated(Model%n_substeps))      deallocate(Model%n_substeps)
    if (allocated(Model%h_prev))          deallocate(Model%h_prev)
    if (allocated(Model%h_curr))          deallocate(Model%h_curr)
    Model%is_initialised = .false.
  end subroutine kernel_deinit

end module gwswex_kernel
