!> Operator-split process subroutines for the explicit solver (§3.5–§3.10).
!!
!! Each subroutine handles one physical process in the explicit cascade:
!!  - precip_partition:   Green-Ampt precipitation partitioning (§3.5)
!!  - gravity_flow:       gravity-driven UZ drainage (§3.6)
!!  - evaporation /
!!    transpiration:      ET sinks (§3.7–§3.8)
!!  - capillary_redist:   capillary redistribution (§3.9)
!!  - geometry_resolution: GW/UZ state reconciliation after cascade (§3.10)
!!
!! All subroutines operate on a single element `ex` and update `M` in-place.
!! Lateral GW/SW fluxes (§3.3) are handled in `gwswex_lateral` (shared with the
!! implicit solver) and applied by the kernel dispatcher before this cascade.
module gwswex_explicit_processes
  use gwswex_constants, only: dp, EPS
  use gwswex_types,     only: gwswex_model
  use gwswex_physics,   only: calc_Se, calc_Kunsat, ga_capacity, stress_factor, &
                              V_gw, V_gw_inv
  use gwswex_geometry,  only: update_geometry
  implicit none

contains

  !------------------------------------------------------------------
  ! §3.5  Precipitation partitioning
  ! Returns infiltration volume via the out argument.
  !------------------------------------------------------------------
  pure subroutine precip_partition(M, ex, dt, infiltration)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt
    real(dp),           intent(out)   :: infiltration
    real(dp) :: P, sw_avail, f_cap, delta_theta, pore_space, gw_deficit, sink_capacity

    P = M%precip_rate(ex) * dt
    sw_avail = M%SW_curr(ex) + P

    ! Green-Ampt capacity (top layer)
    delta_theta = M%theta_s(1, ex) - M%UZ_prev(1, ex) / max(M%d_a(1, ex), EPS)
    delta_theta = max(delta_theta, 0.0_dp)
    f_cap = ga_capacity(M%K_sat(1, ex), M%params%psi_f, delta_theta, M%F_GA(ex))

    ! Sink capacity = top-layer unsaturated pore space + column-wide GW
    ! saturation deficit (= drainable volume between current GWH and the
    ! ground surface).  Including the GW deficit allows ponded surface water
    ! to infiltrate through a top layer whose unsaturated pore space is
    ! already filled by the capillary fringe (UZ_curr ≈ ePV); the surplus
    ! that exceeds the top-layer pore space is routed to GW by gravity_flow
    ! and the over-saturation cap in geometry_resolution.  Without this
    ! term, ponded water on a near-saturated profile is permanently stalled
    ! at the surface even when the WT lies below ground level.
    pore_space = max(M%ePV(1, ex) - M%UZ_curr(1, ex), 0.0_dp)
    gw_deficit = max(M%bnds(1, ex) - M%GWH_curr(ex), 0.0_dp) * M%Sy(1, ex)
    sink_capacity = pore_space + gw_deficit

    infiltration = min(f_cap * dt, sw_avail, sink_capacity)

    M%SW_curr(ex) = sw_avail - infiltration
    M%acc_infiltration(ex) = M%acc_infiltration(ex) + infiltration
    M%acc_precip(ex)       = M%acc_precip(ex) + P
  end subroutine precip_partition

  !------------------------------------------------------------------
  ! §3.6  Gravity flow (top-down sequential per-layer)
  !------------------------------------------------------------------
  pure subroutine gravity_flow(M, ex, dt, infiltration)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt, infiltration
    integer  :: lx, l_star
    real(dp) :: inflow, Se, exfil, uz_free, exfil_l1

    l_star = M%gw_bnd_idx(ex)
    if (l_star <= 0) return  ! fully saturated
    exfil = 0.0_dp
    exfil_l1 = 0.0_dp

    do lx = 1, l_star

      ! inflow
      if (lx == 1) then
        inflow = infiltration
      else
        inflow = exfil   ! exfiltration from previous layer
      end if

      ! K_unsat: compute pre-inflow coefficient for frozen baseline
      Se = calc_Se(M%UZ_curr(lx, ex), M%d_a(lx, ex), M%theta_r(lx, ex), M%theta_s(lx, ex))
      M%K_unsat(lx, ex) = calc_Kunsat(M%K_sat(lx, ex), Se, M%vg_m(lx, ex), M%lambda(lx, ex))

      ! add inflow
      M%UZ_curr(lx, ex) = M%UZ_curr(lx, ex) + inflow

      ! Post-inflow K-averaging when surface infiltration is active (S3.6.3):
      ! prevents the frozen-coefficient stalling artifact where a previously dry
      ! layer's K≈0 blocks the wetting front for an entire sub-step.
      if (inflow > EPS .and. infiltration > EPS) then
        Se = calc_Se(M%UZ_curr(lx, ex), M%d_a(lx, ex), M%theta_r(lx, ex), M%theta_s(lx, ex))
        M%K_unsat(lx, ex) = 0.5_dp * (M%K_unsat(lx, ex) + &
                            calc_Kunsat(M%K_sat(lx, ex), Se, M%vg_m(lx, ex), M%lambda(lx, ex)))
      end if
      M%tc(lx, ex) = M%K_unsat(lx, ex) * dt

      ! update IC and ICratio
      if (inflow > EPS) then
        M%IC(lx, ex) = min(M%IC(lx, ex) + M%tc(lx, ex), M%d_a(lx, ex))
      else
        M%IC(lx, ex) = max(M%IC(lx, ex) - M%tc(lx, ex), 0.0_dp)
      end if
      M%ICratio(lx, ex) = max(M%IC(lx, ex) / max(M%d_a(lx, ex), EPS), M%params%ICratio_min)

      ! exfiltration (bound by drainable water above residual)
      uz_free = max(M%UZ_curr(lx, ex) - M%UZ_eq(lx, ex), 0.0_dp)
      exfil   = min(M%tc(lx, ex), uz_free) * M%ICratio(lx, ex)
      exfil   = min(exfil, max(M%UZ_curr(lx, ex) - M%theta_r(lx, ex) * M%d_a(lx, ex), 0.0_dp))

      ! deduct from this layer
      M%UZ_curr(lx, ex) = M%UZ_curr(lx, ex) - exfil

      if (lx == 1) exfil_l1 = exfil
    end do

    ! boundary layer exfiltration -> GW recharge (S3.6.6)
    if (exfil > EPS) then
      M%GWV_curr(ex) = M%GWV_curr(ex) + exfil
      M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                              M%Sy(:, ex), M%V_cum(:, ex), M%nl)
    end if

    ! update F_GA (S3.6.7)
    if (infiltration > EPS) then
      M%F_GA(ex) = max(min(M%F_GA(ex), M%params%F_min), M%params%F_min) + infiltration
    else
      M%F_GA(ex) = max(min(M%F_GA(ex), M%params%F_min) - exfil_l1, M%params%F_min)
    end if
  end subroutine gravity_flow

  !------------------------------------------------------------------
  ! §3.7  Evaporation (SW then stress-limited UZ top layer)
  !------------------------------------------------------------------
  pure subroutine evaporation(M, ex, dt)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt
    real(dp) :: E_pot, E_sw, E_residual, s1, E_lim, uz_avail, E_uz

    E_pot = M%pet_rate(ex) * dt
    if (E_pot < EPS) return

    ! extract from SW first
    E_sw = min(E_pot, M%SW_curr(ex))
    M%SW_curr(ex) = M%SW_curr(ex) - E_sw
    E_residual = E_pot - E_sw

    if (E_residual < EPS) then
      M%acc_evap(ex) = M%acc_evap(ex) + E_sw
      return
    end if

    ! stress-limited extraction from top UZ layer
    if (M%d_a(1, ex) < EPS .or. M%ePV(1, ex) < EPS) then
      M%acc_evap(ex) = M%acc_evap(ex) + E_sw
      return
    end if
    s1 = M%UZ_curr(1, ex) / M%ePV(1, ex)
    E_lim = E_residual * stress_factor(s1, M%veg(M%vID(ex))%s_h, M%veg(M%vID(ex))%s_e)

    uz_avail = max(M%UZ_curr(1, ex) - M%theta_r(1, ex) * M%d_a(1, ex), 0.0_dp)
    E_uz = min(E_lim, uz_avail)

    M%UZ_curr(1, ex) = M%UZ_curr(1, ex) - E_uz
    M%acc_evap(ex) = M%acc_evap(ex) + E_sw + E_uz
  end subroutine evaporation

  !------------------------------------------------------------------
  ! §3.8  Transpiration (rooted UZ layers + saturated GW extraction)
  !------------------------------------------------------------------
  pure subroutine transpiration(M, ex, dt)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt
    real(dp) :: T_pot, T_pot_l, w_root, s_l, T_lim, uz_avail, T_uz_l, T_uz_total, T_gw
    integer  :: lx, l_star

    T_pot = M%ptt_rate(ex) * dt
    if (T_pot < EPS) return
    if (M%n_root(ex) <= 0) return
    w_root = 1.0_dp / real(M%n_root(ex), dp)

    l_star     = M%gw_bnd_idx(ex)
    T_uz_total = 0.0_dp

    ! UZ extraction (layers 1..l*)
    do lx = 1, max(l_star, 0)
      if (M%is_root(lx, ex) == 0) cycle
      if (M%d_a(lx, ex) < EPS .or. M%ePV(lx, ex) < EPS) cycle

      T_pot_l = T_pot * w_root
      s_l     = M%UZ_curr(lx, ex) / M%ePV(lx, ex)
      T_lim   = T_pot_l * stress_factor(s_l, M%veg(M%vID(ex))%s_w, M%veg(M%vID(ex))%s_star)

      uz_avail = max(M%UZ_curr(lx, ex) - M%theta_r(lx, ex) * M%d_a(lx, ex), 0.0_dp)
      T_uz_l   = min(T_lim, uz_avail)

      M%UZ_curr(lx, ex) = M%UZ_curr(lx, ex) - T_uz_l
      T_uz_total = T_uz_total + T_uz_l
    end do

    ! GW extraction from fully saturated rooted layers (l*+1..nl)
    ! NOTE: boundary layer l* is handled only by UZ extraction above;
    ! including it here would double-count root demand for that layer.
    T_gw = 0.0_dp
    do lx = l_star + 1, M%nl
      if (M%is_root(lx, ex) == 0) cycle
      T_gw = T_gw + T_pot * w_root
    end do
    T_gw = min(T_gw, M%GWV_curr(ex))
    if (T_gw > EPS) then
      M%GWV_curr(ex) = M%GWV_curr(ex) - T_gw
      M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                              M%Sy(:, ex), M%V_cum(:, ex), M%nl)
    end if

    M%acc_transp(ex) = M%acc_transp(ex) + T_uz_total + T_gw
  end subroutine transpiration

  !------------------------------------------------------------------
  ! §3.9  Capillary redistribution (upward from l* to layer 1)
  !------------------------------------------------------------------
  pure subroutine capillary_redist(M, ex, dt, infiltration)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt, infiltration
    integer  :: lx, l_star
    real(dp) :: cap_deficit, cap_flux, uz_avail_below, Se_lx, K_lx, tc_lx

    ! only execute when infiltration is negligible (S3.9.1)
    if (infiltration > EPS) return

    l_star = M%gw_bnd_idx(ex)
    if (l_star <= 0) return

    ! boundary layer: draw from GWV
    cap_deficit = max(M%UZ_eq(l_star, ex) * M%solver%beta_h - M%UZ_curr(l_star, ex), 0.0_dp)
    cap_flux = min(cap_deficit, M%GWV_curr(ex))
    if (cap_flux > EPS) then
      M%UZ_curr(l_star, ex) = M%UZ_curr(l_star, ex) + cap_flux
      M%GWV_curr(ex) = M%GWV_curr(ex) - cap_flux
      M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                              M%Sy(:, ex), M%V_cum(:, ex), M%nl)
    end if

    ! upper layers: supply from layer immediately below
    do lx = l_star - 1, 1, -1
      cap_deficit = max(M%UZ_eq(lx, ex) * M%solver%beta_h - M%UZ_curr(lx, ex), 0.0_dp)
      if (cap_deficit < EPS) cycle

      ! Evaluate K of donor layer from its current state
      Se_lx = calc_Se(M%UZ_curr(lx+1, ex), M%d_a(lx+1, ex), &
                      M%theta_r(lx+1, ex), M%theta_s(lx+1, ex))
      K_lx  = calc_Kunsat(M%K_sat(lx+1, ex), Se_lx, M%vg_m(lx+1, ex), M%lambda(lx+1, ex))
      tc_lx = K_lx * dt

      uz_avail_below = max(M%UZ_curr(lx+1, ex) - M%theta_r(lx+1, ex) * M%d_a(lx+1, ex), 0.0_dp)
      cap_flux = min(cap_deficit, tc_lx, uz_avail_below)
      M%UZ_curr(lx,   ex) = M%UZ_curr(lx,   ex) + cap_flux
      M%UZ_curr(lx+1, ex) = M%UZ_curr(lx+1, ex) - cap_flux
    end do
  end subroutine capillary_redist

  !------------------------------------------------------------------
  ! §3.10  Geometry resolution (layer state transitions, over-sat cap)
  !
  ! gwh_sub_start: the GWH at the start of this sub-step (post-lateral),
  !   used to compute the old d_a for boundary-layer UZ correction.
  !
  ! Physical-water conservation:
  !   Physical = GWV + UZ_total + SW + θ_r × z_sat
  !   V_gw_inv changes z_sat (= GWH - z_bot) without adjusting UZ,
  !   so θ_r × Δz_sat of physical water is created (rise) or
  !   destroyed (drop) per event.  The universal correction is:
  !     ΔUZ = −θ_r × ΔGWH   (applied at the boundary layer)
  !   For rises: deduct θ_r × dh from UZ (removes the spurious water).
  !   For drops: add θ_r × |dh| to UZ (restores the moisture that was
  !     implicit in the saturated zone's θ_r pool).
  !------------------------------------------------------------------
  subroutine geometry_resolution(M, ex, gwh_sub_start)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: gwh_sub_start
    integer  :: lx, l_star_prev, l_star_new, l_sweep_top
    real(dp) :: excess, d_a_old, d_a_new, dh, gwh_before

    l_star_prev = M%gw_bnd_idx(ex)

    ! recompute geometry from current GW
    call update_geometry(M, ex)
    l_star_new = M%gw_bnd_idx(ex)

    ! ── §3.10  Internal boundary-layer θ_r correction (within-layer GWH change) ──
    if (l_star_new == l_star_prev .and. l_star_new > 0) then
      d_a_old = max(M%bnds(l_star_prev, ex) - gwh_sub_start, 0.0_dp)
      d_a_new = M%d_a(l_star_new, ex)
      dh = d_a_old - d_a_new   ! positive for rise, negative for drop

      if (abs(dh) > EPS) then
        M%UZ_curr(l_star_new, ex) = M%UZ_curr(l_star_new, ex) &
                                  - M%theta_r(l_star_new, ex) * dh
        if (M%UZ_curr(l_star_new, ex) < 0.0_dp) then
          M%GWV_curr(ex) = M%GWV_curr(ex) + M%UZ_curr(l_star_new, ex)
          M%UZ_curr(l_star_new, ex) = 0.0_dp
          M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                                  M%Sy(:, ex), M%V_cum(:, ex), M%nl)
          call update_geometry(M, ex)
          l_star_new = M%gw_bnd_idx(ex)
        end if
      end if
    end if

    ! ── §3.10.3  Layer state transitions ──
    if (l_star_new > l_star_prev .and. l_star_prev > 0) then
      ! GW drop: newly exposed layers gain θ_r residual moisture.
      M%UZ_curr(l_star_prev, ex) = M%UZ_curr(l_star_prev, ex) + &
                              M%theta_r(l_star_prev, ex) * max(gwh_sub_start - M%bnds(l_star_prev+1, ex), 0.0_dp)
      if (M%UZ_curr(l_star_prev, ex) > M%ePV(l_star_prev, ex)) then
        excess = M%UZ_curr(l_star_prev, ex) - M%ePV(l_star_prev, ex)
        M%GWV_curr(ex) = M%GWV_curr(ex) + excess
        M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                                M%Sy(:, ex), M%V_cum(:, ex), M%nl)
        M%UZ_curr(l_star_prev, ex) = M%ePV(l_star_prev, ex)
      end if
      do lx = l_star_prev + 1, l_star_new
        M%UZ_curr(lx, ex) = M%theta_r(lx, ex) * M%d_a(lx, ex)
      end do

    else if (l_star_new < l_star_prev .and. l_star_new >= 0) then
      ! GW rise: submerge layers.
      ! Transfer submerged layers' UZ into GWV before zeroing, so that
      ! V_gw_inv can account for the water already present in those layers.
      ! Without this, Σ(UZ_sub) of water is destroyed each cross-layer rise.
      ! The transfer itself raises GWH (via V_gw_inv), which may submerge
      ! additional layers; iterate until stable.  Convergence ratio is
      ! θ_r/Sy ≈ 0.22, so 2-3 passes suffice.
      l_sweep_top = l_star_prev   ! upper sweep bound (layers above here were unsaturated)
      do
        excess = 0.0_dp
        do lx = l_star_new + 1, l_sweep_top
          if (M%UZ_curr(lx, ex) > EPS) then
            excess = excess + M%UZ_curr(lx, ex)
            M%UZ_curr(lx, ex) = 0.0_dp
          end if
        end do
        if (excess < EPS) exit
        M%GWV_curr(ex) = M%GWV_curr(ex) + excess
        M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                                M%Sy(:, ex), M%V_cum(:, ex), M%nl)
        call update_geometry(M, ex)
        l_sweep_top = l_star_new     ! next iteration only sweeps newly submerged
        l_star_new = M%gw_bnd_idx(ex)
      end do
      ! θ_r correction for the total rise (including any additional rise
      ! from the V_sub transfer above).
      dh = M%GWH_curr(ex) - gwh_sub_start
      if (dh > EPS .and. l_star_new > 0) then
        M%UZ_curr(l_star_new, ex) = M%UZ_curr(l_star_new, ex) &
                                  - M%theta_r(l_star_new, ex) * dh
        if (M%UZ_curr(l_star_new, ex) < 0.0_dp) then
          M%GWV_curr(ex) = M%GWV_curr(ex) + M%UZ_curr(l_star_new, ex)
          M%UZ_curr(l_star_new, ex) = 0.0_dp
          M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                                  M%Sy(:, ex), M%V_cum(:, ex), M%nl)
          call update_geometry(M, ex)
          l_star_new = M%gw_bnd_idx(ex)
        end if
      end if
    end if

    ! ── §3.10.4  Over-saturation cap ──
    call update_geometry(M, ex)
    l_star_new = M%gw_bnd_idx(ex)
    do lx = 1, l_star_new
      if (M%UZ_curr(lx, ex) > M%ePV(lx, ex)) then
        excess = M%UZ_curr(lx, ex) - M%ePV(lx, ex)
        gwh_before = M%GWH_curr(ex)
        M%GWV_curr(ex) = M%GWV_curr(ex) + excess
        M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                                M%Sy(:, ex), M%V_cum(:, ex), M%nl)
        M%UZ_curr(lx, ex) = M%ePV(lx, ex)
        ! θ_r correction for the over-sat cap rise
        dh = M%GWH_curr(ex) - gwh_before
        if (dh > EPS) then
          call update_geometry(M, ex)
          l_star_new = M%gw_bnd_idx(ex)
          if (l_star_new > 0) then
            M%UZ_curr(l_star_new, ex) = M%UZ_curr(l_star_new, ex) &
                                      - M%theta_r(l_star_new, ex) * dh
            if (M%UZ_curr(l_star_new, ex) < 0.0_dp) then
              M%GWV_curr(ex) = M%GWV_curr(ex) + M%UZ_curr(l_star_new, ex)
              M%UZ_curr(l_star_new, ex) = 0.0_dp
              M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                                      M%Sy(:, ex), M%V_cum(:, ex), M%nl)
              call update_geometry(M, ex)
              l_star_new = M%gw_bnd_idx(ex)
            end if
          end if
        else
          call update_geometry(M, ex)
          l_star_new = M%gw_bnd_idx(ex)
        end if
      end if
    end do

    ! ── §3.10.1  Surface boundary check (final cap) ──
    if (M%GWH_curr(ex) >= M%bnds(1, ex)) then
      M%SW_curr(ex) = M%SW_curr(ex) + &
                             (M%GWH_curr(ex) - M%bnds(1, ex)) * M%Sy(1, ex)
      M%GWH_curr(ex) = M%bnds(1, ex)
      M%GWV_curr(ex) = V_gw(M%bnds(1, ex), M%bnds(:, ex), M%Sy(:, ex), M%nl)
      do lx = 1, M%nl
        M%UZ_curr(lx, ex) = 0.0_dp
      end do
      call update_geometry(M, ex)
    end if
  end subroutine geometry_resolution

end module gwswex_explicit_processes
