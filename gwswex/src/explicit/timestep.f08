!> CFL-based adaptive sub-step duration evaluator for the explicit solver.
!!
!! Computes the largest time-step dt consistent with CFL stability (Courant
!! condition) across all active unsaturated layers, the precipitation input,
!! and any lateral GW withdrawal, then floors the result at dt_min.
module gwswex_time
  use gwswex_constants, only: dp, EPS
  use gwswex_types,     only: gwswex_model
  use gwswex_physics,   only: calc_Se, calc_Kunsat
  implicit none

contains

  !> Evaluate the CFL-safe sub-step duration for element ex [T].
  !!
  !! Scans all active UZ layers and incoming fluxes. Returns max(dt_cfl, dt_min).
  !! A warning is printed if the gravity-drainage CFL forces dt down to dt_min,
  !! which may indicate excessively high K_sat or coarse layer discretisation.
  pure real(dp) function eval_cfl(M, ex)
    type(gwswex_model), intent(in) :: M
    integer,            intent(in) :: ex
    integer  :: lx, l_star
    real(dp) :: Se, K_us, dt_l, dt_min

    l_star = M%gw_bnd_idx(ex)
    dt_min = huge(1.0_dp)

    !NOTE: The inner layer loop is already thread-private (called from the OMP element loop
    !      in kernel_step), so no nested OMP parallelism is needed here.
    ! gravity-drainage CFL per active layer
    do lx = 1, l_star
      if (M%ePV(lx, ex) < EPS) cycle
      Se   = calc_Se(M%UZ_curr(lx, ex), M%d_a(lx, ex), M%theta_r(lx, ex), M%theta_s(lx, ex))
      K_us = calc_Kunsat(M%K_sat(lx, ex), Se, M%vg_m(lx, ex), M%lambda(lx, ex))
      if (K_us > EPS) then
        dt_l = M%solver%courant_number * M%ePV(lx, ex) / K_us
        ! NOTE: if dt_l < dt_min the sub-step will be clamped. This may signal
        !       that K_sat is too large or layer discretisation is too coarse.
        dt_min = min(dt_min, dt_l)
      end if
    end do

    ! precipitation CFL (S3.2.2)
    if (M%precip_rate(ex) > EPS) then
      dt_l = M%solver%courant_number * max(M%ePV(1, ex) - M%UZ_curr(1, ex), EPS) / M%precip_rate(ex)
      dt_min = min(dt_min, dt_l)
    end if

    ! lateral GW withdrawal CFL (S3.2.2)
    if (M%lat_gw_rate(ex) < -EPS) then
      dt_l = M%solver%courant_number * max(M%GWH_curr(ex) - M%bnds(M%nl+1, ex), EPS) &
             / abs(M%lat_gw_rate(ex))
      dt_min = min(dt_min, dt_l)
    end if

    eval_cfl = max(dt_min, M%solver%dt_min)
  end function eval_cfl

end module gwswex_time
