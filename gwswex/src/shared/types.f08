!> Derived type definitions for the GWSWEX model kernel.
!!
!! Defines all persistent data structures: solver configuration, material
!! properties, model physics parameters, vegetation types, and the top-level
!! `gwswex_model` type that carries the full model state.
module gwswex_types
  use gwswex_constants, only: dp
  implicit none

  !> Numerical solver configuration parameters.
  !! Shared fields apply to both solvers; explicit/implicit sub-blocks are
  !! only used by the corresponding solver.
  type :: solver_config
    ! Shared
    integer  :: solver_type      !! SOLVER_EXPLICIT or SOLVER_IMPLICIT
    integer  :: omp_threads      !! Number of OpenMP threads [-]
    ! Explicit solver parameters
    real(dp) :: courant_number   !! Courant number for CFL stability [-]
    real(dp) :: dt_min           !! Absolute minimum sub-step duration [T]
    real(dp) :: beta_h           !! Capillary hysteresis damping factor [-]
    integer  :: n_trapz          !! Quadrature points for UZ_eq equilibrium integral [-]
    ! Implicit solver parameters
    real(dp) :: picard_tol       !! Convergence tolerance on max |Δh| [L]
    integer  :: picard_max_iter  !! Maximum Picard iterations per step [-]
    ! Safety bounds
    real(dp) :: h_min = -1.0e6_dp  !! Lower bound on matric head [L]; guards against Picard blow-up
  end type solver_config

  !> Single soil material: Van Genuchten + Mualem parameters.
  type :: material
    real(dp) :: K_sat    !! Saturated hydraulic conductivity [LT-1]
    real(dp) :: theta_s  !! Saturated volumetric water content (porosity) [-]
    real(dp) :: theta_r  !! Residual volumetric water content [-]
    real(dp) :: alpha    !! Van Genuchten alpha parameter [L-1]
    real(dp) :: vg_n     !! Van Genuchten n parameter [-]
    real(dp) :: vg_m     !! Derived: m = 1 - 1/n [-]
    real(dp) :: lambda   !! Mualem pore-connectivity exponent [-]
    real(dp) :: Sy       !! Drainable porosity = theta_s - theta_r [-]
  end type material

  !> Global model physics parameters: Green-Ampt infiltration and connectivity.
  !! ET stress thresholds (s_star, s_w, s_h, s_e) are per-vegetation-type
  !! and stored in `vegetation_type`, not here.
  type :: model_params
    real(dp) :: psi_f        !! Green-Ampt suction head fitting parameter [L]
    real(dp) :: F_min        !! Green-Ampt minimum cumulative infiltration [L]
    real(dp) :: ICratio_min  !! Minimum inter-layer connectivity ratio [-]
  end type model_params

  !> Per-vegetation-type root geometry and Laio (2001) ET stress parameters.
  !>
  !! Used by `kernel_set_vegetation` to compute per-element root masks and
  !! to evaluate ET stress during the solver step. Each vegetation type carries
  !! its own stress thresholds, mirroring the per-material soil hydraulic
  !! properties in `material`.
  !!
  !! Transpiration demand is partitioned uniformly across the rooted layers of
  !! each element (1/n_root weighting), with no per-layer density profile.
  type :: vegetation_type
    ! Root structure
    real(dp) :: root_depth        !! Maximum rooting depth below surface [L]

    ! ET stress function parameters
    real(dp) :: s_star            !! Incipient stomatal closure [-]
    real(dp) :: s_w               !! Wilting point [-]
    real(dp) :: s_h               !! Hygroscopic point [-]
    real(dp) :: s_e               !! Capillary-continuity threshold [-]
  end type vegetation_type

  !> Top-level model state container (singleton allocated in `gwswex_kernel`).
  type :: gwswex_model
    ! --- Dimensions ---
    integer :: ne   = 0     ! number of elements
    integer :: nl   = 0     ! number of layers per element
    integer :: nmat = 0     ! number of soil material types in the library
    integer :: nveg = 0     !! Number of vegetation types in the library

    ! --- Config (set at init, immutable during run) ---
    type(solver_config) :: solver
    type(model_params)  :: params
    type(material), allocatable :: mat(:)           ! (nmat) soil material library
    type(vegetation_type), allocatable :: veg(:)   !! (nveg) vegetation library
    integer,  allocatable :: vID(:)                !! (ne) per-element vegetation type ID (1-based)
    real(dp), allocatable :: bnds(:,:)              ! (nl+1, ne) layer boundaries per element
    integer,  allocatable :: is_root(:,:)           ! (nl, ne) root mask, 0/1
    integer,  allocatable :: n_root(:)              ! (ne) number of rooted layers per element

    ! --- Per-layer material lookups (expanded from sID at init for vectorisation) ---
    real(dp), allocatable :: K_sat(:,:)             ! (nl, ne)
    real(dp), allocatable :: theta_s(:,:)           ! (nl, ne)
    real(dp), allocatable :: theta_r(:,:)           ! (nl, ne)
    real(dp), allocatable :: alpha(:,:)             ! (nl, ne)
    real(dp), allocatable :: vg_n(:,:), vg_m(:,:)   ! (nl, ne)
    real(dp), allocatable :: lambda(:,:)            ! (nl, ne)
    real(dp), allocatable :: Sy(:,:)                ! (nl, ne) drainable porosity per layer

    ! --- State (prev/curr pairs) ---
    real(dp), allocatable :: GWH_prev(:), GWH_curr(:)   ! (ne) GW table elevation
    real(dp), allocatable :: GWV_prev(:), GWV_curr(:)   ! (ne) drainable GW volume
    real(dp), allocatable :: SW_prev(:),  SW_curr(:)    ! (ne) surface water depth
    real(dp), allocatable :: UZ_prev(:,:), UZ_curr(:,:) ! (nl, ne) UZ storage
    real(dp), allocatable :: theta_prev(:,:), theta_curr(:,:) ! (nl, ne) derived VWC

    ! --- Implicit solver state (allocated only when solver_type == SOLVER_IMPLICIT) ---
    real(dp), allocatable :: h_prev(:,:)   ! (nl, ne) matric head — previous macro-step
    real(dp), allocatable :: h_curr(:,:)   ! (nl, ne) matric head — current macro-step

    ! --- Persistent parameters (carry memory across steps) ---
    real(dp), allocatable :: IC(:,:)                ! (nl, ne) infiltration-front tracker
    real(dp), allocatable :: ICratio(:,:)           ! (nl, ne) connectivity ratio
    real(dp), allocatable :: F_GA(:)                ! (ne) Green-Ampt cumulative infiltration

    ! --- Ephemeral geometry (recomputed per sub-step) ---
    integer,  allocatable :: gw_bnd_idx(:)          ! (ne)
    real(dp), allocatable :: d_a(:,:)               ! (nl, ne)
    real(dp), allocatable :: ePV(:,:)               ! (nl, ne)
    real(dp), allocatable :: UZ_eq(:,:)             ! (nl, ne)
    real(dp), allocatable :: V_cum(:,:)             ! (nl+1, ne) cumulative drainable volume

    ! --- Ephemeral flux intermediates (per sub-step) ---
    real(dp), allocatable :: K_unsat(:,:)           ! (nl, ne)
    real(dp), allocatable :: tc(:,:)                ! (nl, ne) transfer capacity

    ! --- Forcings (set per macro-step from Python) ---
    real(dp), allocatable :: precip_rate(:)         ! (ne)
    real(dp), allocatable :: pet_rate(:)            ! (ne)
    real(dp), allocatable :: ptt_rate(:)            ! (ne)
    real(dp), allocatable :: lat_gw_rate(:)         ! (ne)
    real(dp), allocatable :: lat_sw_rate(:)         ! (ne)

    ! --- Flux accumulators (summed across sub-steps, read per macro-step) ---
    real(dp), allocatable :: acc_precip(:)          ! (ne)
    real(dp), allocatable :: acc_infiltration(:)    ! (ne)
    real(dp), allocatable :: acc_evap(:)            ! (ne)
    real(dp), allocatable :: acc_transp(:)          ! (ne)
    real(dp), allocatable :: acc_recharge(:)        ! (ne)
    real(dp), allocatable :: acc_runoff(:)          ! (ne)
    real(dp), allocatable :: acc_lat_gw(:)          ! (ne) actual applied lateral GW volume change [L]
    real(dp), allocatable :: acc_lat_sw(:)          ! (ne) actual applied lateral SW depth change [L]
    real(dp), allocatable :: acc_delta_gw(:)        ! (ne) total GW storage change this step [L] (= ΔV_gw)
    real(dp), allocatable :: acc_delta_sw(:)        ! (ne) total SW storage change this step [L]
    real(dp), allocatable :: acc_delta_uz(:)        ! (ne) total UZ storage change this step [L]
    integer,  allocatable :: n_substeps(:)          ! (ne)

    ! --- Status ---
    integer :: status = 0
    logical :: is_initialised = .false.
  end type gwswex_model

end module gwswex_types
