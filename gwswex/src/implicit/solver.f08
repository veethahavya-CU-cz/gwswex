module gwswex_solver_implicit
  !! Full-column implicit Richards solver with Picard iteration and TDMA.
  !!
  !! Solves the 1-D mixed-form Richards equation (Celia et al., 1990):
  !!
  !!   C(h) * ∂h/∂t  =  ∂/∂z [ K(h) * (∂h/∂z + 1) ]  -  sink(z, t)
  !!
  !! using Backward Euler time discretisation, Picard iteration, and the
  !! Thomas algorithm (TDMA) for the tridiagonal linear system.
  !!
  !! Coordinate convention (matching GWSWEX): layer index 1 = surface (top),
  !! layer index nl = domain bottom. z increases upward so z_c(1) > z_c(nl).
  !!
  !! References
  !! ----------
  !! Celia, M. A., Bouloutas, E. T., & Zarba, R. L. (1990).
  !!   A general mass-conservative numerical solution for the unsaturated
  !!   flow equation.  Water Resour. Res., 26(7), 1483–1496.
  use gwswex_constants, only: dp, EPS
  use gwswex_types,     only: gwswex_model
  use gwswex_physics,   only: vg_theta, vg_mualem_kusat, vg_C, harmonic_K, V_gw
  use gwswex_geometry,  only: compute_K_half
  use gwswex_lateral, only: apply_lateral
  implicit none

contains

  ! ===========================================================================
  ! Internal: upstream-weighted inter-node conductance for Richards
  ! ===========================================================================
  !> Upstream (donor-cell) weighting of unsaturated K at layer interfaces.
  !!
  !! For transient infiltration into dry media, the harmonic mean of K(h) at
  !! adjacent cells underestimates the inter-cell conductance when one cell
  !! is wet (h ~ 0) and its neighbour is dry (h << 0, K(h) ~ 0).  This
  !! collapses the wetting front, especially in soils with a steep
  !! retention curve (e.g. sand, alpha >> 1).  Upstream weighting selects
  !! the conductivity of the donor cell (the one supplying flow), which is
  !! the standard remedy in HYDRUS-1D and Forsyth & Kropinski (1997)
  !! "Robust numerical methods for saturated-unsaturated flow with dry
  !! initial conditions in heterogeneous media", Adv. Water Resour. 17(3).
  !!
  !! Sign convention: upper-layer index i is shallower (higher z); positive
  !! Darcy flux is downward (q_{i,i+1} = K * ((h_i - h_{i+1})/dz_int + 1)).
  !! Donor for downward flux is layer i (upper); donor for upward flux is
  !! layer i+1 (lower).  We fall back to harmonic mean when the gradient
  !! sign cannot be classified (transient near-zero gradient).
  pure subroutine compute_K_half_upstream(K, h, z_c, nl, K_half)
    integer,  intent(in)  :: nl
    real(dp), intent(in)  :: K(nl), h(nl), z_c(nl)
    real(dp), intent(out) :: K_half(nl-1)
    integer  :: i
    real(dp) :: dz_int, grad
    do i = 1, nl - 1
      dz_int = z_c(i) - z_c(i+1)
      grad   = (h(i) - h(i+1)) / max(dz_int, EPS) + 1.0_dp
      if (grad >= 0.0_dp) then
        K_half(i) = K(i)            ! downward flux: donor is upper cell
      else
        K_half(i) = K(i+1)          ! upward flux:   donor is lower cell
      end if
    end do
  end subroutine compute_K_half_upstream


  ! ===========================================================================
  ! Public entry point
  ! ===========================================================================

  !> Solve vertical flow implicitly for element ex over macro-step dt.
  !!
  !! Updates Model%GWH_curr(ex), Model%GWV_curr(ex), Model%SW_curr(ex),
  !! Model%UZ_curr(:,ex), Model%theta_curr(:,ex), Model%h_curr(:,ex), and
  !! the flux accumulators.
  subroutine solve_element_implicit(M, ex, dt)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt

    integer  :: nl
    real(dp) :: dz(M%nl), z_c(M%nl)
    real(dp) :: h_old(M%nl), theta_old(M%nl)
    real(dp) :: sink(M%nl), src_lat_gw(M%nl)
    real(dp) :: q_top
    real(dp) :: h_new(M%nl)
    real(dp) :: gw_new, E_act, T_act
    real(dp) :: p_input, p_supply, sw_supply, i_actual
    real(dp) :: K_new(M%nl), C_dummy(M%nl), theta_new(M%nl), K_half_new(M%nl-1)
    real(dp) :: dz_12, q_12, dtheta1, q_top_actual
    logical  :: converged

    nl = M%nl

    ! Carry any ponded water forward from the previous step.
    M%SW_curr(ex) = M%SW_prev(ex)

    ! Apply lateral SW exchange only.  The lateral GW rate is folded into
    ! the Richards solve below as a distributed source over the saturated
    ! portion of the column so that the converged head profile (and hence
    ! `h_to_state`) is mechanically aware of the lateral input.  Mutating
    ! GWH/GWV here would be discarded by `h_to_state` re-deriving the
    ! water table from `h_new`.
    call apply_lateral(M, ex, dt, apply_gw=.false.)

    ! Surface-water supply available to infiltrate during this step.  Ponded
    ! water at the start of the step is offered to the top Neumann flux so
    ! that pre-existing ponding can re-enter the soil.
    sw_supply = max(M%SW_curr(ex), 0.0_dp)

    call compute_dz_zc(M, ex, nl, dz, z_c)

    ! Warm-start head: reuse h_prev if populated, else hydrostatic from
    ! GWH_prev.  h_min enforces a minimum pressure head for numerical
    ! robustness near residual saturation.
    h_old = M%h_prev(:, ex)
    if (all(h_old == 0.0_dp)) h_old = M%GWH_prev(ex) - z_c
    h_old = max(h_old, M%solver%h_min)

    call eval_richards_props(M, ex, nl, h_old, theta_old, &
                             dummy_K=.false., dummy_C=.false.)

    ! Build distributed lateral-GW source [1/T] for the saturated part of
    ! the column.  Flat distribution per saturated thickness; bottom layer
    ! used as the injection point when the column is fully unsaturated.
    call build_lat_gw_source(M, ex, nl, dz, h_old, src_lat_gw)

    ! Initial ET sink at h_old; the Picard loop re-evaluates S and the
    ! actual (E,T) at each iterate to remove the operator-splitting error
    ! of a lagged sink in an otherwise-implicit solve.
    call compute_et_sink(M, ex, nl, dz, theta_old, sink, E_act, T_act, dt)

    ! Top flux: precipitation plus any ponded water, converted to a rate
    ! over dt.  Excess that the soil cannot accept is returned to SW_curr.
    q_top = M%precip_rate(ex) + sw_supply / max(dt, EPS)

    call picard_solve(M, ex, nl, dt, dz, z_c, h_old, theta_old, sink, &
                      src_lat_gw, q_top, h_new, E_act, T_act, converged)

    call h_to_state(M, ex, nl, dz, z_c, h_new, gw_new)

    M%GWH_curr(ex)    = gw_new
    M%GWV_curr(ex)    = V_gw(gw_new, M%bnds(:, ex), M%Sy(:, ex), nl)
    M%h_curr(:, ex)   = h_new

    ! Actual top-boundary water flux from a layer-1 mass balance on the
    ! converged head profile.  Includes the layer-1 lateral-GW source so
    ! that the inferred surface flux remains the actual cross-interface
    ! atmospheric exchange when the column is surface-saturated.
    call eval_richards_props(M, ex, nl, h_new, theta_new, &
                             dummy_K=.false., dummy_C=.true., K=K_new, C=C_dummy)
    if (nl >= 2) then
      call compute_K_half_upstream(K_new, h_new, z_c, nl, K_half_new)
      dz_12 = z_c(1) - z_c(2)
      q_12  = K_half_new(1) * ((h_new(1) - h_new(2)) / max(dz_12, EPS) + 1.0_dp)
    else
      q_12 = 0.0_dp
    end if
    dtheta1      = theta_new(1) - theta_old(1)
    ! `sink(1)` returned from picard_solve is the net layer sink already,
    ! i.e. ET sink at the converged iterate minus src_lat_gw(1).  Hence the
    ! layer-1 budget reads dz1*dtheta1/dt = q_top - q_12 - sink(1)*dz1.
    q_top_actual = q_12 + dz(1) * dtheta1 / dt + sink(1) * dz(1)

    p_input  = M%precip_rate(ex) * dt
    p_supply = p_input + sw_supply

    ! Bounds: 0 <= I_actual <= P_supply.  Residual ponds in SW_curr.
    ! acc_runoff, acc_recharge, and storage deltas are derived centrally
    ! in gwswex_mass_balance.output_calc after the solver returns.
    i_actual      = min(max(q_top_actual * dt, 0.0_dp), p_supply)
    M%SW_curr(ex) = max(p_supply - i_actual, 0.0_dp)

    M%acc_precip(ex)       = M%acc_precip(ex)       + p_input
    M%acc_infiltration(ex) = M%acc_infiltration(ex) + i_actual
    M%acc_evap(ex)         = M%acc_evap(ex)         + E_act
    M%acc_transp(ex)       = M%acc_transp(ex)       + T_act
    ! acc_lat_sw was written by apply_lateral above.  acc_lat_gw is the
    ! prescribed lateral-GW input over the step (the volume the source
    ! injected into the saturated zone); any over-extraction beyond
    ! residual saturation is held off by the constitutive theta-floor and
    ! manifests as a closure residual in output_calc.
    M%acc_lat_gw(ex)       = M%acc_lat_gw(ex)       + M%lat_gw_rate(ex) * dt

    M%n_substeps(ex) = M%n_substeps(ex) + 1
  end subroutine solve_element_implicit

  ! ===========================================================================
  ! Internal: distributed lateral-GW source for Picard
  ! ===========================================================================
  !> Convert the per-element lateral-GW rate `lat_gw_rate(ex)` [LT-1] into a
  !! per-layer volumetric source `src(:)` [1/T] suitable to be subtracted
  !! from the Richards sink.
  !!
  !! Distribution rules (consistent with the explicit solver, where the
  !! lateral GW volume is poured directly into the saturated store and
  !! `V_gw_inv` redistributes head accordingly):
  !!   * If the column has saturated layers (h(l) >= 0), spread the rate
  !!     uniformly across the saturated thickness so that
  !!     sum(src(l) * dz(l)) = lat_gw_rate.
  !!   * If the column is fully unsaturated (deep water table below the
  !!     domain), inject at the bottom layer.
  pure subroutine build_lat_gw_source(M, ex, nl, dz, h, src)
    type(gwswex_model), intent(in)  :: M
    integer,            intent(in)  :: ex, nl
    real(dp),           intent(in)  :: dz(nl), h(nl)
    real(dp),           intent(out) :: src(nl)
    integer  :: l
    real(dp) :: rate, l_sat

    src = 0.0_dp
    rate = M%lat_gw_rate(ex)
    if (abs(rate) < EPS) return

    l_sat = 0.0_dp
    do l = 1, nl
      if (h(l) >= 0.0_dp) l_sat = l_sat + dz(l)
    end do

    if (l_sat > EPS) then
      do l = 1, nl
        if (h(l) >= 0.0_dp) src(l) = rate / l_sat
      end do
    else
      ! Deep water table: inject at the column base.
      if (dz(nl) > EPS) src(nl) = rate / dz(nl)
    end if
  end subroutine build_lat_gw_source

  ! ===========================================================================
  ! Internal: helper geometry
  ! ===========================================================================

  !> Compute layer thicknesses dz(nl) and centre elevations z_c(nl) for element ex.
  pure subroutine compute_dz_zc(M, ex, nl, dz, z_c)
    type(gwswex_model), intent(in)  :: M
    integer,            intent(in)  :: ex, nl
    real(dp),           intent(out) :: dz(nl), z_c(nl)
    integer :: l
    do l = 1, nl
      dz(l)  = M%bnds(l, ex) - M%bnds(l+1, ex)
      z_c(l) = 0.5_dp * (M%bnds(l, ex) + M%bnds(l+1, ex))
    end do
  end subroutine compute_dz_zc

  ! ===========================================================================
  ! Internal: soil property evaluation
  ! ===========================================================================

  !> Evaluate (theta, K, C) at head profile h for element ex.
  !! dummy_K / dummy_C flags skip computing K or C arrays (for init calls).
  subroutine eval_richards_props(M, ex, nl, h, theta, dummy_K, dummy_C, K, C)
    type(gwswex_model), intent(in)  :: M
    integer,            intent(in)  :: ex, nl
    real(dp),           intent(in)  :: h(nl)
    real(dp),           intent(out) :: theta(nl)
    logical,            intent(in)  :: dummy_K, dummy_C
    real(dp), optional, intent(out) :: K(nl), C(nl)
    integer  :: l
    real(dp) :: h_l

    do l = 1, nl
      h_l = min(h(l), 0.0_dp)   ! VG functions expect h <= 0
      if (h(l) >= 0.0_dp) then
        theta(l) = M%theta_s(l, ex)
      else
        theta(l) = vg_theta(h_l, M%theta_r(l,ex), M%theta_s(l,ex), &
                            M%alpha(l,ex), M%vg_n(l,ex), M%vg_m(l,ex))
      end if
      if (present(K) .and. .not. dummy_K) then
        K(l) = vg_mualem_kusat(h(l), M%K_sat(l,ex), M%alpha(l,ex), M%vg_n(l,ex), M%lambda(l,ex))
      end if
      if (present(C) .and. .not. dummy_C) then
        C(l) = vg_C(h(l), M%theta_r(l,ex), M%theta_s(l,ex), M%alpha(l,ex), &
                    M%vg_n(l,ex), M%Sy(l,ex))   ! use drainable porosity Sy, not Ss
      end if
    end do
  end subroutine eval_richards_props

  ! ===========================================================================
  ! Internal: ET sink
  ! ===========================================================================

  !> Compute ET sink vector sink(nl) [1/T] for the implicit solver.
  !! Also returns E_act and T_act [L] (actual fluxes over dt).
  !!
  !! Evaporation: extracted from the surface layer (l=1) at a rate limited
  !! by available soil moisture above theta_h.
  !!
  !! Transpiration: distributed uniformly across rooted layers (1/n_root
  !! weighting), scaled by the Laio stress factor.
  subroutine compute_et_sink(M, ex, nl, dz, theta, sink, E_act, T_act, dt)
    type(gwswex_model), intent(in)  :: M
    integer,            intent(in)  :: ex, nl
    real(dp),           intent(in)  :: dz(nl), theta(nl), dt
    real(dp),           intent(out) :: sink(nl), E_act, T_act

    integer  :: l
    real(dp) :: s_l, stress_e, stress_t, pet_l, ptt_l, w_root
    real(dp) :: theta_s1, theta_r1, dz1

    sink = 0.0_dp
    E_act = 0.0_dp
    T_act = 0.0_dp

    ! Evaporation: surface layer only
    theta_s1 = M%theta_s(1, ex)
    theta_r1 = M%theta_r(1, ex)
    dz1      = dz(1)
    if (dz1 > EPS) then
      s_l = (theta(1) - theta_r1) / max(theta_s1 - theta_r1, EPS)
      s_l = max(0.0_dp, min(s_l, 1.0_dp))
      stress_e = laio_stress(s_l, M%veg(M%vID(ex))%s_h, M%veg(M%vID(ex))%s_e)
      pet_l    = stress_e * M%pet_rate(ex)
      ! limit by available moisture above hygroscopic point
      pet_l = min(pet_l, max((theta(1) - M%theta_r(1,ex)) * dz1, 0.0_dp) / dt)
      sink(1) = sink(1) + pet_l / max(dz1, EPS)
      E_act   = pet_l * dt
    end if

    ! Transpiration: uniformly distributed across rooted layers, Laio stress
    if (M%n_root(ex) <= 0) return
    w_root = 1.0_dp / real(M%n_root(ex), dp)
    do l = 1, nl
      if (M%is_root(l, ex) == 0) cycle
      if (dz(l) <= EPS) cycle
      s_l = (theta(l) - M%theta_r(l,ex)) / max(M%theta_s(l,ex) - M%theta_r(l,ex), EPS)
      s_l = max(0.0_dp, min(s_l, 1.0_dp))
      stress_t = laio_stress(s_l, M%veg(M%vID(ex))%s_w, M%veg(M%vID(ex))%s_star)
      ptt_l    = stress_t * M%ptt_rate(ex) * w_root
      ptt_l    = min(ptt_l, max((theta(l) - M%theta_r(l,ex)) * dz(l), 0.0_dp) / dt)
      sink(l)  = sink(l) + ptt_l / max(dz(l), EPS)
      T_act    = T_act + ptt_l * dt
    end do
  end subroutine compute_et_sink

  !> Laio (2001) piecewise-linear stress function.
  !! Returns 0 at s <= s_lo, 1 at s >= s_hi, linear between.
  elemental real(dp) function laio_stress(s, s_lo, s_hi)
    real(dp), intent(in) :: s, s_lo, s_hi
    if (s <= s_lo) then
      laio_stress = 0.0_dp
    else if (s >= s_hi) then
      laio_stress = 1.0_dp
    else
      laio_stress = (s - s_lo) / max(s_hi - s_lo, EPS)
    end if
  end function laio_stress

  ! ===========================================================================
  ! Internal: Picard iteration
  ! ===========================================================================

  !> Run Picard iteration for one element over one macro-step.
  !!
  !! The ET sink is re-evaluated at each iterate using the current head profile
  !! (theta_k), making the sink treatment consistent with the implicit Richards
  !! solve. Returns the converged sink-derived E_act, T_act.
  subroutine picard_solve(M, ex, nl, dt, dz, z_c, h_old, theta_old, &
                           sink, src_lat_gw, q_top, h_new, E_act, T_act, converged)
    type(gwswex_model), intent(in)  :: M
    integer,            intent(in)  :: ex, nl
    real(dp),           intent(in)  :: dt, dz(nl), z_c(nl)
    real(dp),           intent(in)  :: h_old(nl), theta_old(nl)
    real(dp),           intent(inout) :: sink(nl)
    real(dp),           intent(in)  :: src_lat_gw(nl)
    real(dp),           intent(in)  :: q_top
    real(dp),           intent(out) :: h_new(nl)
    real(dp),           intent(out) :: E_act, T_act
    logical,            intent(out) :: converged

    integer  :: k
    real(dp) :: h_k(nl), h_try(nl), h_prev(nl), dh(nl)
    real(dp) :: theta_k(nl), K_k(nl), C_k(nl)
    real(dp) :: K_half(nl-1)
    real(dp) :: a(nl), b(nl), c_diag(nl), d_rhs(nl)
    real(dp) :: dh_max, dh_cap, damp
    logical  :: surface_ponded

    surface_ponded = .false.
    h_k     = h_old
    converged = .false.
    E_act = 0.0_dp
    T_act = 0.0_dp

    ! Per-iteration update cap.  Near saturation the retention-curve
    ! capacitance C(h) collapses by three orders of magnitude (the vg_C
    ! floor is Sy * 1.0e-3); when the Celia tridiagonal is inverted at an
    ! iterate already close to h = 0 the implied Δh of a small storage
    ! residual can exceed the column depth, which is unphysical.  Bounding
    ! the per-Picard update to a fraction of the thinnest active layer
    ! restores monotone convergence without altering the converged fixed
    ! point (Paniconi & Putti 1994; Lehmann & Ackerer 1998).
    dh_cap = 0.5_dp * minval(dz)

    do k = 1, M%solver%picard_max_iter
      h_prev = h_k

      ! Properties at current iterate
      call eval_richards_props(M, ex, nl, h_k, theta_k, &
                               dummy_K=.false., dummy_C=.false., K=K_k, C=C_k)

      ! Refresh ET sink at current iterate (implicit treatment of S(h))
      call compute_et_sink(M, ex, nl, dz, theta_k, sink, E_act, T_act, dt)

      ! Subtract the lateral-GW source: net layer sink = ET sink - lateral
      ! source (positive sink leaves the layer; positive source adds water).
      sink = sink - src_lat_gw

      ! K_half at layer interfaces (nl-1 values) — upstream weighting for
      ! robust wetting-front advance into dry media (sand, layered profiles)
      call compute_K_half_upstream(K_k, h_k, z_c, nl, K_half)

      ! Build tridiagonal system
      call build_picard_system(nl, dt, dz, z_c, h_old, h_k, &
                               K_k, K_half, C_k, theta_old, theta_k, &
                               sink, q_top, surface_ponded, &
                               a, b, c_diag, d_rhs)

      ! Solve via TDMA
      call solve_tdma(nl, a, b, c_diag, d_rhs, h_try)

      ! Damp oversized Picard updates (capacitance-collapse protection)
      dh     = h_try - h_k
      dh_max = maxval(abs(dh))
      if (dh_max > dh_cap) then
        damp  = dh_cap / dh_max
        h_try = h_k + damp * dh
      end if

      ! Dynamic top BC switch: Neumann → Dirichlet when h[1] > 0
      if (.not. surface_ponded .and. h_try(1) > 0.0_dp) then
        surface_ponded = .true.
        h_k = h_try
        cycle
      end if
      if (surface_ponded .and. h_try(1) < 0.0_dp .and. q_top < EPS) then
        surface_ponded = .false.
        h_k = h_try
        cycle
      end if

      h_k = h_try

      ! Convergence check: max |Δh| < picard_tol
      if (maxval(abs(h_k - h_prev)) < M%solver%picard_tol) then
        converged = .true.
        exit
      end if
    end do

    h_new = h_k
  end subroutine picard_solve

  ! ===========================================================================
  ! Internal: tridiagonal matrix assembly (Celia 1990 mixed form)
  ! ===========================================================================

  !> Assemble tridiagonal arrays (a, b, c, d) for the full soil column.
  !!
  !! Mixed-form linearisation (Celia et al. 1990):
  !!   C^k * (h^{k+1} - h^k) / dt  +  (θ^k - θ^n) / dt
  !!     = ∂/∂z [ K^k * (∂h^{k+1}/∂z + 1) ]  -  sink
  !!
  !! Layer indexing: 1 = surface (top, highest z), nl = domain bottom.
  !! z increases upward: z_c(1) > z_c(2) > ... > z_c(nl).
  !!
  !! a(i): sub-diagonal   — coupling to h(i-1) (shallower, HIGHER z)
  !! b(i): main diagonal
  !! c(i): super-diagonal — coupling to h(i+1) (deeper,    LOWER  z)
  subroutine build_picard_system(nl, dt, dz, z_c, h_old, h_k, &
                                  K, K_half, C, theta_old, theta_k, &
                                  sink, q_top, surface_ponded, &
                                  a, b, c_diag, d_rhs)
    integer,  intent(in)  :: nl
    real(dp), intent(in)  :: dt, dz(nl), z_c(nl)
    real(dp), intent(in)  :: h_old(nl), h_k(nl)
    real(dp), intent(in)  :: K(nl), K_half(nl-1), C(nl)
    real(dp), intent(in)  :: theta_old(nl), theta_k(nl)
    real(dp), intent(in)  :: sink(nl), q_top
    logical,  intent(in)  :: surface_ponded
    real(dp), intent(out) :: a(nl), b(nl), c_diag(nl), d_rhs(nl)

    integer  :: i
    real(dp) :: K_up, K_dn, up_coef, dn_coef, dz_up, dz_dn, gravity_term

    ! h_old and K (cell-centre) are conceptual inputs to the Picard system;
    ! K_half (interface) is used here directly; h_old appears in the full
    ! mixed-form formula but cancels in the Celia (1990) correction form.
    associate(unused_h_old => h_old(1), unused_K => K(1)); end associate

    do i = 1, nl
      ! Upper interface conductance (i=1: top BC)
      if (i == 1) then
        K_up    = 0.0_dp
        up_coef = 0.0_dp
      else
        K_up    = K_half(i-1)
        dz_up   = z_c(i-1) - z_c(i)    ! > 0 (shallower is higher z)
        up_coef = K_up / max(dz(i) * dz_up, EPS)
      end if

      ! Lower interface conductance (i=nl: impermeable base)
      if (i == nl) then
        K_dn    = 0.0_dp
        dn_coef = 0.0_dp
      else
        K_dn    = K_half(i)
        dz_dn   = z_c(i) - z_c(i+1)    ! > 0 (deeper is lower z)
        dn_coef = K_dn / max(dz(i) * dz_dn, EPS)
      end if

      b(i)      = -(C(i) / dt + up_coef + dn_coef)
      a(i)      =  up_coef       ! coupling to h(i-1): shallower
      c_diag(i) =  dn_coef       ! coupling to h(i+1): deeper

      ! Gravity flux contribution (upward gradient of K)
      gravity_term = (K_up - K_dn) / max(dz(i), EPS)

      ! RHS: storage + gravity (negated system, matching the sign convention
      ! b*h[i] + a*h[i-1] + c*h[i+1] = d where b is negative-definite)
      d_rhs(i) = -(C(i) / dt) * h_k(i) &
                 + (theta_k(i) - theta_old(i)) / dt &
                 - gravity_term &
                 + sink(i)
    end do

    ! --- Top boundary condition ---
    if (.not. surface_ponded) then
      ! Neumann: prescribed surface flux q_top [L/T], positive into soil
      d_rhs(1) = d_rhs(1) - q_top / max(dz(1), EPS)
    else
      ! Dirichlet: ponded surface h(1) = 0
      a(1)     = 0.0_dp
      b(1)     = -1.0_dp
      c_diag(1) = 0.0_dp
      d_rhs(1) = 0.0_dp
    end if
  end subroutine build_picard_system

  ! ===========================================================================
  ! Internal: TDMA (Thomas algorithm)
  ! ===========================================================================

  !> Solve tridiagonal system A x = d via the Thomas algorithm O(n).
  !!
  !! a: sub-diagonal  (coupling to x(i-1))   — a(1) is unused
  !! b: main diagonal
  !! c: super-diagonal (coupling to x(i+1)) — c(n) is unused
  !! d: right-hand side
  !! x: solution (overwritten on output, treated as d workspace)
  subroutine solve_tdma(n, a, b, c, d, x)
    integer,  intent(in)    :: n
    real(dp), intent(in)    :: a(n), b(n), c(n)
    real(dp), intent(inout) :: d(n)
    real(dp), intent(out)   :: x(n)

    real(dp) :: cp(n), dp_w(n), denom
    integer  :: i

    ! Forward sweep
    if (abs(b(1)) < EPS) then
      ! Near-zero pivot: set x = 0 and return (degenerate system)
      x = 0.0_dp
      return
    end if
    cp(1)  = c(1) / b(1)
    dp_w(1) = d(1) / b(1)

    do i = 2, n
      denom = b(i) - a(i) * cp(i-1)
      if (abs(denom) < EPS) then
        x = 0.0_dp
        return
      end if
      cp(i)   = c(i) / denom
      dp_w(i) = (d(i) - a(i) * dp_w(i-1)) / denom
    end do

    ! Back substitution
    x(n) = dp_w(n)
    do i = n-1, 1, -1
      x(i) = dp_w(i) - cp(i) * x(i+1)
    end do
  end subroutine solve_tdma

  ! ===========================================================================
  ! Internal: Extract state from converged head profile
  ! ===========================================================================

  !> Convert converged head profile h(nl) back to GWSWEX state variables.
  !! Updates Model%UZ_curr, Model%theta_curr, Model%SW_curr for element ex.
  !! Returns gw_new [L].  The drainable-volume change (V_gw(gw_new) -
  !! V_gw(GWH_prev)) is recovered downstream by gwswex_mass_balance from
  !! the (GWV_curr, GWV_prev) state pair and is therefore no longer
  !! returned here.
  subroutine h_to_state(M, ex, nl, dz, z_c, h, gw_new)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex, nl
    real(dp),           intent(in)    :: dz(nl), z_c(nl), h(nl)
    real(dp),           intent(out)   :: gw_new

    integer  :: l, k, l_gw
    real(dp) :: frac, theta_new(nl)

    ! --- theta from converged h ---
    do l = 1, nl
      if (h(l) >= 0.0_dp) then
        theta_new(l) = M%theta_s(l, ex)
      else
        theta_new(l) = vg_theta(h(l), M%theta_r(l,ex), M%theta_s(l,ex), &
                                M%alpha(l,ex), M%vg_n(l,ex), M%vg_m(l,ex))
      end if
    end do

    ! --- Locate water table ---
    ! The water table is the top of contiguous saturation propagating up
    ! from the bottom of the column: the shallowest layer k such that
    ! h(k), h(k+1), ..., h(nl) are all >= 0.  A positive h(1) without
    ! saturation of the layer immediately below it represents surface
    ! ponding on an unsaturated column (imposed by the Dirichlet top BC,
    ! or produced by excess Neumann inflow on a sub-step that has not yet
    ! conducted the water downward), NOT a rising water table.  The old
    ! test `h(1) >= 0 => l_gw = 0` conflated the two and caused the water
    ! table to snap to the surface the first time the top layer ponded,
    ! which is unphysical and destroys the (UZ, GWV) partition.
    k = nl + 1
    do l = nl, 1, -1
      if (h(l) >= 0.0_dp) then
        k = l
      else
        exit
      end if
    end do

    if (k == 1) then
      ! Column fully saturated from bottom to top: WT at surface.
      gw_new = M%bnds(1, ex)
      l_gw   = 0
    else if (k > nl) then
      ! No saturated layer in domain: WT below the base.
      gw_new = M%bnds(nl+1, ex)
      l_gw   = nl
    else
      ! WT between unsaturated layer k-1 and saturated layer k.
      frac   = -h(k-1) / max(h(k) - h(k-1), EPS)
      gw_new = z_c(k-1) + frac * (z_c(k) - z_c(k-1))
      l_gw   = k - 1
    end if

    ! --- Update UZ_curr per GWSWEX convention ---
    !   l <= l_gw (unsaturated): UZ = theta * dz   (mobile soil moisture)
    !   l >  l_gw (saturated):   UZ = 0
    !
    ! Saturated-layer water is split between drainable storage (GWV = Sy*dz,
    ! tracked separately) and residual storage (theta_r * dz). The residual
    ! component is reported diagnostically as `theta_r * max(GWH, 0)` in the
    ! external storage formula and is therefore intentionally NOT credited to
    ! UZ_curr — doing so would double-count theta_r * sat_thickness and
    ! produce a spurious storage-gain bias of order theta_r * Delta(WT) over
    ! a run with a moving water table.  This convention matches the explicit
    ! solver (Phase C: UZ_curr = 0 for fully saturated layers).
    do l = 1, nl
      if (l_gw == 0) then
        ! Fully saturated column
        M%UZ_curr(l, ex) = 0.0_dp
      else if (l <= l_gw) then
        M%UZ_curr(l, ex) = theta_new(l) * dz(l)
      else
        M%UZ_curr(l, ex) = 0.0_dp
      end if
      M%theta_curr(l, ex) = theta_new(l)
    end do

    ! --- Surface water: updated by solve_element_implicit after this call ---
    ! SW_curr is written by the caller using the actual infiltration deficit
    ! (precip + ponded_in - infiltration).  h_to_state intentionally does
    ! not touch SW_curr.
  end subroutine h_to_state

end module gwswex_solver_implicit
