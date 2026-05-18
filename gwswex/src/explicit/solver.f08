!> Operator-split explicit solver for one macro-step per element.
!!
!! Implements the CFL-adaptive sub-stepping scheme described in model-physics.md
!! §3.3–§3.10.  The two main pathways are:
!!   - UZ-active cascade: precipitation partitioning, gravity flow, ET, capillary
!!     redistribution, and geometry resolution (§3.5–§3.10).
!!   - UZ-inactive: handles the fully saturated case and the re-activation gate (§3.1).
!!
!! Entry point: `solve_element(M, ex, dt_macro)`.
module gwswex_solver_explicit
  use gwswex_constants, only: dp, EPS
  use gwswex_types,     only: gwswex_model
  use gwswex_geometry,  only: update_geometry
  use gwswex_lateral,            only: apply_lateral
  use gwswex_explicit_processes, only: precip_partition, gravity_flow, &
                              evaporation, transpiration, capillary_redist, &
                              geometry_resolution
  use gwswex_time,      only: eval_cfl
  use gwswex_physics,   only: V_gw, V_gw_inv
  implicit none

contains

  !------------------------------------------------------------------
  ! Solve one macro-step for element ex.  Internal sub-stepping is
  ! handled here; the caller never sees sub-steps.
  !------------------------------------------------------------------
  subroutine solve_element_explicit(M, ex, dt_macro)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt_macro
    real(dp) :: t_rem, dt_sub, gwh_sub_start
    integer  :: l_star

    t_rem = dt_macro
    M%n_substeps(ex) = 0

    ! copy t-1 state into working buffer
    M%GWH_curr(ex) = M%GWH_prev(ex)
    M%GWV_curr(ex) = M%GWV_prev(ex)
    M%SW_curr(ex) = M%SW_prev(ex)
    M%UZ_curr(:, ex) = M%UZ_prev(:, ex)

    !================= SUB-STEP LOOP =================
    do while (t_rem > EPS)

      ! --- CFL evaluation ---
      call update_geometry(M, ex)    ! fresh geometry for CFL
      dt_sub = eval_cfl(M, ex)
      dt_sub = min(dt_sub, t_rem)
      t_rem  = t_rem - dt_sub
      M%n_substeps(ex) = M%n_substeps(ex) + 1

      ! --- S3.3: Lateral fluxes ---
      call apply_lateral(M, ex, dt_sub)

      ! --- S3.4: Geometry update (after lateral) ---
      call update_geometry(M, ex)
      l_star = M%gw_bnd_idx(ex)
      gwh_sub_start = M%GWH_curr(ex)   ! snapshot post-lateral GWH

      ! --- Saturation check (S3.1) ---
      if (l_star <= 0) then
        ! ====== UZ-INACTIVE PATHWAY ======
        call uz_inactive_pathway(M, ex, dt_sub, gwh_sub_start)
      else
        ! ====== UZ-ACTIVE PATHWAY ======
        call uz_active_cascade(M, ex, dt_sub, .true., gwh_sub_start)
      end if

    end do  ! sub-step loop
  end subroutine solve_element_explicit

  !------------------------------------------------------------------
  ! UZ-active cascade (S3.5 to S3.10)
  ! with_et: .true. for full cascade, .false. for partial (re-activation)
  !------------------------------------------------------------------
  subroutine uz_active_cascade(M, ex, dt, with_et, gwh_sub_start)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt
    logical,            intent(in)    :: with_et
    real(dp),           intent(in)    :: gwh_sub_start
    real(dp) :: infiltration

    ! S3.5: Precipitation partitioning
    call precip_partition(M, ex, dt, infiltration)

    ! S3.6: Gravity flow
    call gravity_flow(M, ex, dt, infiltration)

    ! S3.7-S3.8: ET (skipped in partial cascade from re-activation)
    if (with_et) then
      call evaporation(M, ex, dt)
      call transpiration(M, ex, dt)
    end if

    ! S3.9: Capillary redistribution
    call capillary_redist(M, ex, dt, infiltration)

    ! S3.10: Geometry resolution
    call geometry_resolution(M, ex, gwh_sub_start)
  end subroutine uz_active_cascade

  !------------------------------------------------------------------
  ! UZ-inactive pathway (S3.1 Phases A-C + re-activation gate)
  !------------------------------------------------------------------
  subroutine uz_inactive_pathway(M, ex, dt, gwh_sub_start)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt
    real(dp),           intent(in)    :: gwh_sub_start
    real(dp) :: P, E_pot, E_act, T_pot, T_gw
    integer  :: lx
    ! gwh_sub_start passed for API symmetry with uz_active_cascade;
    ! reserved for geometry-correction diagnostics.
    associate(unused_gwh => gwh_sub_start); end associate

    ! Phase A: cap excess GW to surface
    if (M%GWH_curr(ex) > M%bnds(1, ex)) then
      M%SW_curr(ex) = M%SW_curr(ex) + &
                             (M%GWH_curr(ex) - M%bnds(1, ex)) * M%Sy(1, ex)
      M%GWH_curr(ex) = M%bnds(1, ex)
      M%GWV_curr(ex) = V_gw(M%bnds(1, ex), M%bnds(:, ex), M%Sy(:, ex), M%nl)
    end if

    ! Phase B1: P -> SW
    P = M%precip_rate(ex) * dt
    M%SW_curr(ex) = M%SW_curr(ex) + P
    M%acc_precip(ex) = M%acc_precip(ex) + P

    ! Phase B2: E from SW (no stress)
    E_pot = M%pet_rate(ex) * dt
    E_act = min(E_pot, M%SW_curr(ex))
    M%SW_curr(ex) = M%SW_curr(ex) - E_act
    M%acc_evap(ex) = M%acc_evap(ex) + E_act

    ! Phase B3: T from GW (roots submerged)
    T_pot = M%ptt_rate(ex) * dt
    T_gw = min(T_pot, M%GWV_curr(ex), M%K_sat(1, ex) * dt)
    if (T_gw > EPS) then
      M%GWV_curr(ex) = M%GWV_curr(ex) - T_gw
      M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), &
                              M%Sy(:, ex), M%V_cum(:, ex), M%nl)
    end if
    M%acc_transp(ex) = M%acc_transp(ex) + T_gw

    ! Phase B4: Re-activation gate
    if (M%GWH_curr(ex) < M%bnds(1, ex)) then
      ! Reverse Phase B1 precip accumulation: precip_partition in the
      ! partial cascade will re-add P with correct infiltration/runoff
      ! partitioning.  Without this reversal P would be double-counted.
      M%SW_curr(ex) = M%SW_curr(ex) - P
      M%acc_precip(ex) = M%acc_precip(ex) - P
      ! re-activate UZ
      call update_geometry(M, ex)
      do lx = 1, M%gw_bnd_idx(ex)
        M%UZ_curr(lx, ex) = M%theta_r(lx, ex) * M%d_a(lx, ex)
      end do
      ! partial cascade (no ET); use current GWH as sub-start reference
      call uz_active_cascade(M, ex, dt, .false., M%GWH_curr(ex))
      return
    end if

    ! Phase C: fully saturated, all UZ = 0
    M%UZ_curr(:, ex) = 0.0_dp
  end subroutine uz_inactive_pathway

end module gwswex_solver_explicit
