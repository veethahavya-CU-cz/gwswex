!> Lateral GW/SW flux application — shared between explicit and implicit solvers.
!!
!! Provides `apply_lateral`, the per-element Euler-step update for prescribed
!! lateral SW and GW exchange rates.  The subroutine is invoked at the top of
!! each macro-step by both solver dispatchers in `gwswex_kernel`, before the
!! vertical Richards / cascade machinery is exercised.
!!
!! Reference: model-physics.md §3.3 (Lateral fluxes).  Sign convention:
!! positive `lat_*_rate` adds water to the storage; the GW path clamps at the
!! domain bottom and re-derives `GWH_curr` from the (possibly clamped) volume.
module gwswex_lateral
  use gwswex_constants, only: dp, EPS
  use gwswex_types,     only: gwswex_model
  use gwswex_physics,   only: V_gw_inv
  implicit none

  private
  public :: apply_lateral

contains

  !> Apply prescribed lateral SW and GW rates over the macro-step `dt`.
  !!
  !! Updates `M%SW_curr(ex)`, the SW lateral accumulator, and (when
  !! `apply_gw` is true or absent) `M%GWV_curr(ex)`, `M%GWH_curr(ex)`, and
  !! `M%acc_lat_gw(ex)`.  GW volume is clamped at zero; GW head is clamped
  !! at the domain bottom.
  !!
  !! `apply_gw` (optional, default `.true.`) controls whether the GW path
  !! is executed.  The implicit solver passes `apply_gw=.false.` because
  !! it injects `lat_gw_rate` as a distributed source inside the Picard
  !! Richards solve so that the new water table is mechanically aware of
  !! the lateral inflow; a pre-solve mutation of `GWH_curr` would be
  !! erased by `h_to_state` re-deriving the head from the converged head
  !! profile, which knows nothing about the lateral budget.
  pure subroutine apply_lateral(M, ex, dt, apply_gw)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    real(dp),           intent(in)    :: dt
    logical,  intent(in), optional    :: apply_gw
    real(dp) :: gw_vol_pre
    logical  :: do_gw

    do_gw = .true.
    if (present(apply_gw)) do_gw = apply_gw

    M%SW_curr(ex) = M%SW_curr(ex) + M%lat_sw_rate(ex) * dt
    M%SW_curr(ex) = max(M%SW_curr(ex), 0.0_dp)
    M%acc_lat_sw(ex) = M%acc_lat_sw(ex) + M%lat_sw_rate(ex) * dt

    if (.not. do_gw) return

    ! GW: volume-based update using tracked GWV
    gw_vol_pre = M%GWV_curr(ex)
    M%GWV_curr(ex) = M%GWV_curr(ex) + M%lat_gw_rate(ex) * dt
    if (M%GWV_curr(ex) > EPS) then
      M%GWH_curr(ex) = V_gw_inv(M%GWV_curr(ex), M%bnds(:, ex), M%Sy(:, ex), &
                                       M%V_cum(:, ex), M%nl)
    else
      M%GWV_curr(ex) = 0.0_dp
      M%GWH_curr(ex) = M%bnds(M%nl+1, ex)  ! domain bottom
    end if
    M%GWH_curr(ex) = max(M%GWH_curr(ex), M%bnds(M%nl+1, ex))

    ! Store ACTUAL volume change (may differ from prescribed when clamped
    ! at domain bottom), so that output_calc mass balance is correct.
    M%acc_lat_gw(ex) = M%acc_lat_gw(ex) + (M%GWV_curr(ex) - gw_vol_pre)
  end subroutine apply_lateral

end module gwswex_lateral
