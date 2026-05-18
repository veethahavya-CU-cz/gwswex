!> Global constants for the GWSWEX model kernel.
!!
!! Defines precision kind `dp`, numerical tolerance `EPS`, return-status codes,
!! solver type enumerators, and solver-specific bounds.
module gwswex_constants
  implicit none

  !> Double precision kind parameter (IEEE-754 64-bit).
  integer, parameter :: dp = selected_real_kind(15, 307)

  !> Numerical zero tolerance. Values smaller than this are treated as zero [-].
  real(dp), parameter :: EPS = 1.0e-12_dp

  ! --- Return status codes ---
  integer, parameter :: STAT_OK    = 0  !! Successful completion
  integer, parameter :: STAT_WARN  = 1  !! Non-fatal warning
  integer, parameter :: STAT_ERROR = 2  !! Fatal error

  ! --- Solver type identifiers ---
  integer, parameter :: SOLVER_EXPLICIT = 1  !! Operator-split explicit solver
  integer, parameter :: SOLVER_IMPLICIT = 2  !! Mixed-Richards implicit Picard solver

end module gwswex_constants
