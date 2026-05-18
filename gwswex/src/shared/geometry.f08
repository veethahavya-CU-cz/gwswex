!> Geometric pre-computations for the GWSWEX solvers.
!!
!! This module provides two distinct sets of subroutines, each scoped to a
!! specific solver:
!!
!!   §A  Explicit-solver geometry (explicit solver only)
!!        update_geometry — computes l*, d_a, ePV, UZ_eq, and V_cum for one
!!        element.  Called before every CFL sub-step in the explicit cascade
!!        (model-physics.md §3.4) and once at IC setup in kernel_set_ic.
!!        NOT called by the implicit solver during stepping or at any other
!!        point in the implicit solve path.
!!
!!   §B  Implicit-solver geometry (implicit solver only)
!!        compute_K_half — assembles the inter-node harmonic conductance array
!!        K_half(1:nl−1) from the per-node conductivity vector K(nl).
!!        Called once per Picard iteration inside solve_element_implicit
!!        (model-physics.md §4.3).
module gwswex_geometry
  use gwswex_constants, only: dp, EPS
  use gwswex_types,     only: gwswex_model
  use gwswex_physics,   only: vg_integrate, V_gw, harmonic_K
  implicit none

contains

  ! ===========================================================================
  ! §A  Explicit-solver geometry — called only by the explicit solver
  ! ===========================================================================

  !> Compute layer geometry, UZ equilibrium profiles, and V_cum for element ex.
  !!
  !! Updates the following fields of M for the given element ex:
  !!   - M%gw_bnd_idx(ex)    — index l* of the GW boundary layer
  !!   - M%d_a(lx, ex)       — active (unsaturated) thickness per layer [L]
  !!   - M%ePV(lx, ex)       — effective pore volume per layer [L]
  !!   - M%UZ_eq(lx, ex)     — equilibrium UZ storage at hydrostatic state [L]
  !!   - M%V_cum(lx, ex)     — cumulative drainable volume at each boundary [L]
  !!
  !! Must be called before any solver or CFL evaluation within a sub-step.
  !! Called at the start of each CFL sub-step (§3.4) and at geometry
  !! resolution (§3.10) in the explicit cascade, and once at IC setup
  !! (kernel_set_ic, explicit-only path).  NOT called by the implicit solver.
  !!
  !! @param M   Model singleton (intent inout — updates geometric fields).
  !! @param ex  Element index (1-based).
  pure subroutine update_geometry(M, ex)
    type(gwswex_model), intent(inout) :: M
    integer,            intent(in)    :: ex
    integer  :: lx, l_star
    real(dp) :: gw_h, h_top, h_bot

    gw_h   = M%GWH_curr(ex)
    l_star = M%nl   ! default: GW at or below domain bottom

    ! locate the GW boundary layer l*
    do lx = 1, M%nl
      if (M%bnds(lx+1, ex) <= gw_h .and. gw_h < M%bnds(lx, ex)) then
        l_star = lx
        exit
      end if
      if (gw_h >= M%bnds(lx, ex)) then
        l_star = lx - 1
        exit
      end if
    end do
    if (gw_h >= M%bnds(1, ex)) l_star = 0   ! fully saturated
    M%gw_bnd_idx(ex) = l_star

    ! active thickness d_a, effective pore volume ePV, and VG equilibrium storage UZ_eq
    do lx = 1, M%nl
      if (lx < l_star) then
        M%d_a(lx, ex) = M%bnds(lx, ex) - M%bnds(lx+1, ex)
      else if (lx == l_star) then
        M%d_a(lx, ex) = max(0.0_dp, M%bnds(lx, ex) - gw_h)
      else
        M%d_a(lx, ex) = 0.0_dp
      end if
      M%ePV(lx, ex) = M%d_a(lx, ex) * M%theta_s(lx, ex)

      ! UZ_eq via VG trapezoidal integral (model-physics.md §3.4.2)
      if (lx <= l_star .and. M%d_a(lx, ex) > EPS) then
        h_top = gw_h - M%bnds(lx, ex)       ! < 0 (above GW table)
        if (lx < l_star) then
          h_bot = gw_h - M%bnds(lx+1, ex)
        else
          h_bot = 0.0_dp                     ! boundary layer bottom = GW table
        end if
        M%UZ_eq(lx, ex) = vg_integrate(h_top, h_bot, &
          M%theta_r(lx, ex), M%theta_s(lx, ex), M%alpha(lx, ex), &
          M%vg_n(lx, ex), M%vg_m(lx, ex), M%solver%n_trapz)
      else
        M%UZ_eq(lx, ex) = 0.0_dp
      end if
    end do

    ! V_cum(lx) = cumulative drainable volume at boundary lx (used by V_GW_inv)
    M%V_cum(M%nl+1, ex) = 0.0_dp
    do lx = M%nl, 1, -1
      M%V_cum(lx, ex) = M%V_cum(lx+1, ex) + M%Sy(lx, ex) * (M%bnds(lx, ex) - M%bnds(lx+1, ex))
    end do
  end subroutine update_geometry

  ! ===========================================================================
  ! §B  Implicit-solver geometry — called only by solve_element_implicit
  !     (model-physics.md §4.3)
  ! ===========================================================================

  !> Assemble the inter-node harmonic conductance array K_half(1:nl−1).
  !!
  !! K_half(i) = thickness-weighted harmonic mean conductivity at the interface
  !! between layer i (shallower) and layer i+1 (deeper):
  !!
  !!   K_½(i) = (dz_i + dz_{i+1}) / (dz_i/K_i + dz_{i+1}/K_{i+1})
  !!
  !! Used once per Picard iteration to assemble the tridiagonal conductance
  !! matrix (model-physics.md §4.3).  The harmonic mean is the correct
  !! inter-node conductance for 1-D Darcy flow through adjacent layers of
  !! differing thickness and conductivity.
  !!
  !! Not called by the explicit solver, which uses per-layer K(Se) directly
  !! in the CFL evaluation and gravity-flow step without a matrix solve.
  !!
  !! @param K      Per-node unsaturated conductivity at current Picard iterate [LT⁻¹], shape (nl).
  !! @param dz     Layer thicknesses [L], shape (nl).
  !! @param nl     Number of layers.
  !! @param K_half Output: inter-node conductance [LT⁻¹], shape (nl−1).
  pure subroutine compute_K_half(K, dz, nl, K_half)
    integer,  intent(in)  :: nl
    real(dp), intent(in)  :: K(nl), dz(nl)
    real(dp), intent(out) :: K_half(nl-1)
    integer :: i
    do i = 1, nl - 1
      K_half(i) = harmonic_K(K(i), K(i+1), dz(i), dz(i+1))
    end do
  end subroutine compute_K_half

end module gwswex_geometry
