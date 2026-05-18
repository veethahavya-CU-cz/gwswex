!> Soil hydraulic and ET physics functions for the GWSWEX model kernel.
!!
!! Provides elemental and pure functions implementing the constitutive relations
!! used by both the explicit operator-split solver (model-physics.md §3) and
!! the implicit Richards solver (model-physics.md §4).
!!
!! Organisation:
!!   §1  Van Genuchten (1980) retention and saturation            — shared
!!   §2  Mualem–van Genuchten (1976/1980) unsaturated conductivity — shared
!!   §3  Specific moisture capacity C(h)                          — implicit solver
!!   §4  Explicit-solver effective saturation from UZ storage      — explicit solver
!!   §5  Laio et al. (2001) ET stress                             — shared
!!   §6  Green & Ampt (1911) infiltration capacity                — explicit solver
!!   §7  Drainable GW volume V_GW and its inverse                 — shared
!!
!! References (with DOI):
!!   van Genuchten, M. Th. (1980).  A closed-form equation for predicting the
!!     hydraulic conductivity of unsaturated soils.  Soil Sci. Soc. Am. J.,
!!     44, 892–898.  doi:10.2136/sssaj1980.03615995004400050002x
!!   Mualem, Y. (1976).  A new model for predicting the hydraulic conductivity
!!     of unsaturated porous media.  Water Resour. Res., 12, 513–522.
!!     doi:10.1029/WR012i003p00513
!!   Celia, M. A., Bouloutas, E. T., & Zarba, R. L. (1990).  A general
!!     mass-conservative numerical solution for the unsaturated flow equation.
!!     Water Resour. Res., 26(7), 1483–1496.  doi:10.1029/WR026i007p01483
!!   Laio, F., Porporato, A., Ridolfi, L., & Rodriguez-Iturbe, I. (2001).
!!     Plants in water-controlled ecosystems: active role in hydrologic processes.
!!     Adv. Water Resour., 24, 707–723.  doi:10.1016/S0309-1708(01)00005-7
!!   Green, W. H., & Ampt, G. A. (1911).  Studies on soil physics, Part I:
!!     the flow of air and water through soils.  J. Agric. Sci., 4, 1–24.
!!     doi:10.1017/S0021859600001441
module gwswex_physics
  use gwswex_constants, only: dp, EPS
  implicit none

contains

  ! ===========================================================================
  ! §1  Van Genuchten (1980) retention and saturation
  !     Ref: doi:10.2136/sssaj1980.03615995004400050002x
  !     Used by both solvers and the geometry module.
  ! ===========================================================================

  !> Van Genuchten (1980) volumetric water content θ(h).
  !!
  !! θ(h) = θ_r + (θ_s − θ_r) · [1 + (α|h|)^n]^{−m}  for h < 0;  θ_s  for h ≥ 0.
  !!
  !! @param h       Matric pressure head [L]; h ≥ 0 → saturated (returns θ_s).
  !! @param theta_r Residual volumetric water content [-].
  !! @param theta_s Saturated volumetric water content [-].
  !! @param alpha   VG pore-size parameter α [L⁻¹].
  !! @param vg_n   VG pore-size distribution index n [-]; n > 1.
  !! @param vg_m   VG shape parameter m [-]; m = 1 − 1/n (pass pre-computed).
  elemental real(dp) function vg_theta(h, theta_r, theta_s, alpha, vg_n, vg_m)
    real(dp), intent(in) :: h, theta_r, theta_s, alpha, vg_n, vg_m
    if (h >= 0.0_dp) then
      vg_theta = theta_s
    else
      vg_theta = theta_r + (theta_s - theta_r) / (1.0_dp + (alpha * abs(h))**vg_n)**vg_m
    end if
  end function vg_theta

  !> Van Genuchten (1980) inverse retention: matric head h(θ).
  !!
  !! For unsaturated θ ∈ (θ_r, θ_s), inverts θ(h) analytically:
  !!   Se = (θ - θ_r) / (θ_s - θ_r)
  !!   h  = -1/α · (Se^{-1/m} - 1)^{1/n}
  !!
  !! Used to translate the explicit-solver (θ-based) state into a
  !! physically consistent matric-head warm start when switching to the
  !! implicit solver mid-simulation.  Returns 0 for θ ≥ θ_s (saturated)
  !! and `h_min` for θ ≤ θ_r (residual saturation, h → -∞).
  elemental real(dp) function vg_h_inv(theta, theta_r, theta_s, alpha, vg_n, vg_m, h_min)
    real(dp), intent(in) :: theta, theta_r, theta_s, alpha, vg_n, vg_m, h_min
    real(dp) :: se
    if (theta >= theta_s - EPS) then
      vg_h_inv = 0.0_dp
      return
    end if
    se = (theta - theta_r) / max(theta_s - theta_r, EPS)
    if (se <= EPS) then
      vg_h_inv = h_min
      return
    end if
    vg_h_inv = -((se**(-1.0_dp/vg_m) - 1.0_dp)**(1.0_dp/vg_n)) / max(alpha, EPS)
    vg_h_inv = max(vg_h_inv, h_min)
  end function vg_h_inv

  !> Van Genuchten (1980) effective saturation Se(h) from matric head.
  !!
  !! Se(h) = [1 + (α|h|)^n]^{−m}  for h < 0;  1  for h ≥ 0.
  !!
  !! This function is used by the implicit solver (which carries h as primary
  !! state) and by the Python reference implementation.  The explicit solver
  !! does not track h; it derives Se from volumetric storage via calc_Se (§4).
  !!
  !! @param h      Matric pressure head [L]; h ≥ 0 → returns 1.
  !! @param alpha  VG pore-size parameter α [L⁻¹].
  !! @param vg_n  VG pore-size distribution index n [-].
  elemental real(dp) function vg_Se(h, alpha, vg_n)
    real(dp), intent(in) :: h, alpha, vg_n
    real(dp) :: m
    if (h >= 0.0_dp) then
      vg_Se = 1.0_dp
      return
    end if
    m = 1.0_dp - 1.0_dp / vg_n
    vg_Se = 1.0_dp / (1.0_dp + (alpha * abs(h))**vg_n)**m
  end function vg_Se

  !> Composite trapezoidal integration of θ(h) over a layer for equilibrium UZ storage.
  !!
  !! Computes ∫_{h_bot}^{h_top} θ(h) dh by the n_trapz-interval trapezoidal rule.
  !! Used by the geometry module to evaluate UZ_eq for each layer above the GW table.
  !!
  !! @param h_top   Matric head at the top boundary of the layer [L] (< 0 above GW).
  !! @param h_bot   Matric head at the bottom boundary of the layer [L].
  !! @param theta_r Residual volumetric water content [-].
  !! @param theta_s Saturated volumetric water content [-].
  !! @param alpha   VG pore-size parameter α [L⁻¹].
  !! @param vg_n   VG pore-size distribution index n [-].
  !! @param vg_m   VG shape parameter m [-]; m = 1 − 1/n.
  !! @param n_trapz Number of trapezoidal sub-intervals; solver_config%n_trapz.
  pure real(dp) function vg_integrate(h_top, h_bot, theta_r, theta_s, alpha, vg_n, vg_m, n_trapz)
    real(dp), intent(in) :: h_top, h_bot, theta_r, theta_s, alpha, vg_n, vg_m
    integer,  intent(in) :: n_trapz
    real(dp) :: dh, h_i, s
    integer  :: i
    dh = (h_top - h_bot) / real(n_trapz, dp)
    s  = 0.5_dp * (vg_theta(h_bot, theta_r, theta_s, alpha, vg_n, vg_m) &
                  + vg_theta(h_top, theta_r, theta_s, alpha, vg_n, vg_m))
    do i = 1, n_trapz - 1
      h_i = h_bot + real(i, dp) * dh
      s = s + vg_theta(h_i, theta_r, theta_s, alpha, vg_n, vg_m)
    end do
    vg_integrate = s * abs(dh)
  end function vg_integrate

  ! ===========================================================================
  ! §2  Mualem–van Genuchten unsaturated conductivity
  !     Refs: Mualem (1976) doi:10.1029/WR012i003p00513
  !           van Genuchten (1980) doi:10.2136/sssaj1980.03615995004400050002x
  !
  !  Two interface variants exist because the two solvers track different
  !  primary state variables:
  !    vg_mualem_kusat — input is matric head h        (implicit solver)
  !    calc_Kunsat     — input is effective saturation Se (explicit solver)
  !  Both implement the identical Mualem (1976) pore-connectivity formula
  !  K = K_sat · Se^λ · [1 − (1 − Se^{1/m})^m]².
  !  harmonic_K provides the thickness-weighted inter-node conductance needed
  !  for the finite-difference assembly of the implicit solver.
  ! ===========================================================================

  !> Mualem–van Genuchten unsaturated hydraulic conductivity K(h) from matric head.
  !!
  !! K(h) = K_sat · Se(h)^λ · [1 − (1 − Se(h)^{1/m})^m]²  for h < 0;
  !!        K_sat                                             for h ≥ 0.
  !!
  !! Used by the implicit solver, which carries h as its primary state variable.
  !! For the explicit solver (which carries volumetric UZ storage), use calc_Kunsat.
  !!
  !! @param h      Matric pressure head [L]; h ≥ 0 → returns K_sat.
  !! @param K_sat  Saturated hydraulic conductivity [LT⁻¹].
  !! @param alpha  VG pore-size parameter α [L⁻¹].
  !! @param vg_n  VG pore-size distribution index n [-].
  !! @param lam   Mualem pore-connectivity exponent λ [-]; typically 0.5.
  elemental real(dp) function vg_mualem_kusat(h, K_sat, alpha, vg_n, lam)
    real(dp), intent(in) :: h, K_sat, alpha, vg_n, lam
    real(dp) :: Se, vg_m, Se_c
    vg_m = 1.0_dp - 1.0_dp / vg_n
    if (h >= 0.0_dp) then
      vg_mualem_kusat = K_sat
      return
    end if
    Se = 1.0_dp / (1.0_dp + (alpha * abs(h))**vg_n)**vg_m
    Se_c = max(Se, EPS)
    vg_mualem_kusat = K_sat * Se_c**lam * (1.0_dp - (1.0_dp - Se_c**(1.0_dp/vg_m))**vg_m)**2
  end function vg_mualem_kusat

  !> Mualem–van Genuchten unsaturated hydraulic conductivity K(Se) from effective saturation.
  !!
  !! K = K_sat · Se^λ · [1 − (1 − Se^{1/m})^m]².
  !!
  !! Used by the explicit solver and the CFL evaluator, which carry volumetric
  !! UZ storage as their primary state and compute Se = (θ − θ_r)/(θ_s − θ_r)
  !! via calc_Se, without resolving the matric head h.
  !! For the implicit solver (which carries h), use vg_mualem_kusat.
  !!
  !! @param K_sat  Saturated hydraulic conductivity [LT⁻¹].
  !! @param Se     Effective saturation [-]; clamped to [EPS, 1] internally.
  !! @param vg_m  VG shape parameter m [-]; m = 1 − 1/n.
  !! @param lam   Mualem pore-connectivity exponent λ [-]; typically 0.5.
  elemental real(dp) function calc_Kunsat(K_sat, Se, vg_m, lam)
    real(dp), intent(in) :: K_sat, Se, vg_m, lam
    real(dp) :: Se_c
    Se_c = max(Se, EPS)
    calc_Kunsat = K_sat * Se_c**lam * (1.0_dp - (1.0_dp - Se_c**(1.0_dp/vg_m))**vg_m)**2
  end function calc_Kunsat

  !> Thickness-weighted harmonic mean conductivity at the interface between layers i and i+1.
  !!
  !! K_½ = (dz_i + dz_{i+1}) / (dz_i/K_i + dz_{i+1}/K_{i+1}).
  !!
  !! The harmonic mean is the correct inter-node conductance for 1-D Darcy flow
  !! through adjacent layers of differing thickness and conductivity.  Used by
  !! the implicit solver to assemble the tridiagonal conductance matrix
  !! (model-physics.md §4.3).
  !!
  !! @param K_i    Conductivity at node i (shallower layer) [LT⁻¹].
  !! @param K_ip1  Conductivity at node i+1 (deeper layer) [LT⁻¹].
  !! @param dz_i   Thickness of the shallower layer [L].
  !! @param dz_ip1 Thickness of the deeper layer [L].
  elemental real(dp) function harmonic_K(K_i, K_ip1, dz_i, dz_ip1)
    real(dp), intent(in) :: K_i, K_ip1, dz_i, dz_ip1
    real(dp) :: denom
    denom = dz_i / max(K_i, EPS) + dz_ip1 / max(K_ip1, EPS)
    harmonic_K = (dz_i + dz_ip1) / max(denom, EPS)
  end function harmonic_K

  ! ===========================================================================
  ! §3  Specific moisture capacity C(h) — implicit solver only
  !
  !  C(h) = dθ/dh is the analytical derivative of the van Genuchten θ(h) curve.
  !  The formula is explicit in van Genuchten (1980), Eq. 20
  !  (doi:10.2136/sssaj1980.03615995004400050002x).
  !
  !  Its role in the mass-conservative mixed-form Richards discretisation —
  !  evaluating C at the current Picard iterate h^k rather than at θ —
  !  is due to Celia et al. (1990), doi:10.1029/WR026i007p01483.
  ! ===========================================================================

  !> Van Genuchten (1980) specific moisture capacity C(h) = dθ/dh.
  !!
  !! For h < 0 (unsaturated):
  !!   C = (θ_s − θ_r)·m·n·α^n·|h|^{n−1}·Se^{1+1/m}
  !!   — exact analytical derivative of the VG retention curve
  !!     (van Genuchten 1980, Eq. 20; doi:10.2136/sssaj1980.03615995004400050002x).
  !!
  !! For h ≥ 0 (saturated or ponded): returns Sy = θ_s − θ_r (drainable
  !! porosity), not the compressibility-based specific storativity S_s.  For
  !! an unconfined phreatic zone S_y is the correct storage coefficient at the
  !! free surface (Bear 1979, §5.3).
  !!
  !! The unsaturated branch is floored at Sy × 10⁻³ to prevent near-zero pivot
  !! values in the TDMA solve near saturation.
  !!
  !! @param h       Matric pressure head [L]; h ≥ 0 → saturated branch.
  !! @param theta_r Residual volumetric water content [-].
  !! @param theta_s Saturated volumetric water content [-].
  !! @param alpha   VG pore-size parameter α [L⁻¹].
  !! @param vg_n   VG pore-size distribution index n [-]; n > 1.
  !! @param Sy     Drainable porosity (θ_s − θ_r) [-]; returned for h ≥ 0 and
  !!               used as the numerical floor Sy × 10⁻³ in the unsaturated branch.
  elemental real(dp) function vg_C(h, theta_r, theta_s, alpha, vg_n, Sy)
    real(dp), intent(in) :: h, theta_r, theta_s, alpha, vg_n, Sy
    real(dp) :: vg_m, x, Se
    if (h >= 0.0_dp) then
      vg_C = Sy
      return
    end if
    vg_m = 1.0_dp - 1.0_dp / vg_n
    x    = (alpha * abs(h))**vg_n
    Se   = 1.0_dp / (1.0_dp + x)**vg_m
    ! C = (θ_s − θ_r) · m · n · α^n · |h|^(n−1) · Se^(1+1/m)
    vg_C = (theta_s - theta_r) * vg_m * vg_n * (alpha**vg_n) * (abs(h)**(vg_n - 1.0_dp)) &
           * Se**(1.0_dp + 1.0_dp/vg_m)
    vg_C = max(vg_C, Sy * 1.0e-3_dp)   ! numerical floor near saturation
  end function vg_C

  ! ===========================================================================
  ! §4  Effective saturation from UZ storage — explicit solver only
  !
  !  The explicit operator-split solver carries volumetric water storage UZ [L]
  !  (water depth per unit area) as its primary state variable for each layer,
  !  not matric head h.  Effective saturation is therefore computed from
  !  θ = UZ/d_a rather than from h via the VG relation.  This is fundamentally
  !  different from vg_Se (§1), which requires h.
  ! ===========================================================================

  !> Effective saturation Se from volumetric UZ water storage (explicit solver).
  !!
  !! Se = (θ − θ_r) / (θ_s − θ_r),   θ = UZ / d_a,   clamped to [0, 1].
  !!
  !! This is specific to the explicit solver.  The implicit solver uses the
  !! VG relation vg_Se(h, α, n) instead, because it carries h as state.
  !!
  !! @param uz      Volumetric UZ water per unit area in the layer [L].
  !! @param d_a     Active (unsaturated) layer thickness [L]; floored at EPS.
  !! @param theta_r Residual volumetric water content [-].
  !! @param theta_s Saturated volumetric water content [-].
  elemental real(dp) function calc_Se(uz, d_a, theta_r, theta_s)
    real(dp), intent(in) :: uz, d_a, theta_r, theta_s
    real(dp) :: theta
    theta = uz / max(d_a, EPS)
    calc_Se = max(0.0_dp, min((theta - theta_r) / max(theta_s - theta_r, EPS), 1.0_dp))
  end function calc_Se

  ! ===========================================================================
  ! §5  Laio et al. (2001) ET stress — shared
  !     Ref: doi:10.1016/S0309-1708(01)00005-7
  ! ===========================================================================

  !> Laio et al. (2001) piecewise-linear soil-moisture stress function.
  !!
  !! Returns 0 at s ≤ s_lo (water stress fully suppresses flux),
  !!         1 at s ≥ s_hi (no stress; actual flux equals potential flux),
  !!         and a linear ramp in between.
  !!
  !! Used for both evaporation and transpiration stress in both solvers.
  !! Per-vegetation-type stress thresholds (s_h, s_e, s_w, s_star) are passed
  !! from M%veg(M%vID(ex)) at the call site; see model-physics.md §3.7–§3.8.
  !!
  !! @param s     Relative soil saturation [-]; typically Se or θ/θ_s.
  !! @param s_lo  Lower saturation threshold below which the flux is zero [-].
  !! @param s_hi  Upper saturation threshold above which the flux is potential [-].
  elemental real(dp) function stress_factor(s, s_lo, s_hi)
    real(dp), intent(in) :: s, s_lo, s_hi
    if (s <= s_lo) then
      stress_factor = 0.0_dp
    else if (s >= s_hi) then
      stress_factor = 1.0_dp
    else
      stress_factor = (s - s_lo) / (s_hi - s_lo)
    end if
  end function stress_factor

  ! ===========================================================================
  ! §6  Green & Ampt (1911) infiltration capacity — explicit solver
  !     Ref: doi:10.1017/S0021859600001441
  ! ===========================================================================

  !> Green–Ampt (1911) ponded infiltration capacity f [LT⁻¹].
  !!
  !! f = K_sat · (1 + ψ_f · Δθ / F_cum)
  !!
  !! where ψ_f is the effective wetting-front suction head [L] and F_cum is
  !! the cumulative infiltrated depth since ponding onset [L].  As F_cum → 0⁺,
  !! f → ∞ and the rainfall rate limits infiltration; as F_cum grows, f decays
  !! asymptotically toward K_sat.
  !!
  !! @param K_sat       Saturated hydraulic conductivity [LT⁻¹].
  !! @param psi_f       Wetting-front suction head ψ_f [L]; model_params%psi_f.
  !! @param delta_theta Soil moisture deficit = θ_s − θ_initial [-].
  !! @param F_cum       Cumulative infiltrated depth [L]; floored at EPS.
  elemental real(dp) function ga_capacity(K_sat, psi_f, delta_theta, F_cum)
    real(dp), intent(in) :: K_sat, psi_f, delta_theta, F_cum
    ga_capacity = K_sat * (1.0_dp + psi_f * delta_theta / max(F_cum, EPS))
  end function ga_capacity

  ! ===========================================================================
  ! §7  Drainable GW volume V_GW and its inverse — shared
  ! ===========================================================================

  !> Drainable GW volume V_GW(h) by summing Sy·dz over saturated layers [L].
  !!
  !! V_GW = Σ_{l: bnds(l+1) < h} Sy(l) · [min(h, bnds(l)) − bnds(l+1)]
  !!       + Sy(1) · max(h − bnds(1), 0)    [above-surface ponding extension]
  !!
  !! @param h     Current water table elevation [L] above datum.
  !! @param bnds  Layer boundary elevations [L], shape (nl+1), top first.
  !! @param Sy    Drainable porosity per layer [-], shape (nl).
  !! @param nl    Number of layers.
  pure real(dp) function V_gw(h, bnds, Sy, nl)
    real(dp), intent(in) :: h, bnds(:), Sy(:)
    integer,  intent(in) :: nl
    integer  :: l
    real(dp) :: contrib
    V_gw = 0.0_dp
    do l = 1, nl
      contrib = max(0.0_dp, min(h, bnds(l)) - bnds(l+1))
      V_gw = V_gw + Sy(l) * contrib
    end do
    if (h > bnds(1)) V_gw = V_gw + Sy(1) * (h - bnds(1))
  end function V_gw

  !> Inverse drainable GW volume V_GW^{−1}(V): water table elevation from volume.
  !!
  !! Locates the layer containing volume V using the pre-computed cumulative
  !! V_cum array, then linearly interpolates the water table elevation within
  !! that layer.  Returns bnds(nl+1) (domain bottom) when V ≤ 0.
  !!
  !! @param V      Target drainable GW volume [L].
  !! @param bnds   Layer boundary elevations [L], shape (nl+1), top first.
  !! @param Sy     Drainable porosity per layer [-], shape (nl).
  !! @param V_cum  Cumulative drainable volume at each boundary [L], shape (nl+1);
  !!               V_cum(1) = total column volume, V_cum(nl+1) = 0.
  !! @param nl     Number of layers.
  pure real(dp) function V_gw_inv(V, bnds, Sy, V_cum, nl)
    real(dp), intent(in) :: V, bnds(:), Sy(:), V_cum(:)
    integer,  intent(in) :: nl
    integer :: l
    V_gw_inv = -1.0_dp   ! error sentinel (should not be returned under normal operation)
    if (V <= 0.0_dp) then
      V_gw_inv = bnds(nl+1)   ! domain bottom
      return
    end if
    ! above-surface ponding: WT above column top
    if (V > V_cum(1)) then
      V_gw_inv = bnds(1) + (V - V_cum(1)) / max(Sy(1), EPS)
      return
    end if
    ! locate layer by descending through cumulative volumes
    do l = 1, nl
      if (V > V_cum(l+1)) then
        V_gw_inv = bnds(l+1) + (V - V_cum(l+1)) / max(Sy(l), EPS)
        return
      end if
    end do
    ! V ≤ V_cum(nl+1) = 0: WT at domain bottom
    V_gw_inv = bnds(nl+1)
  end function V_gw_inv

end module gwswex_physics
