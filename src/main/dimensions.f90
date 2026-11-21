!> Maximum array sizes in paraest.
!>
!> Each subroutine used to repeat its own "integer, parameter :: NCM = 50, ...";
!> now they are defined once here and each subroutine does "use dimensions".
module dimensions
  implicit none

  integer, parameter :: NCM   = 50     !< components
  integer, parameter :: NGM   = 20     !< dispersive groups
  integer, parameter :: NDM   = 1520   !< experimental data points (in total)
  integer, parameter :: NSetM = 100    !< data sets

  !> Fitted parameters of each type (g or k, q, alpha, association energy
  !> and volume, Tc, dc): how many of each can be fitted.
  integer, parameter :: NPTIPOM = 20

  !> "Identical" parameters (those that copy the value of a fitted one),
  !> adding up all those of the same type.
  integer, parameter :: NIDM = 104

  !> Total number of fitted parameters. There are 7 types; the g's can add
  !> up to 2 temperature-dependence terms (g' and g''):
  !> 3*NPTIPOM + 6*NPTIPOM = 9*20 = 180.
  integer, parameter :: NPARM = 9*NPTIPOM

  !> Length of the ZXSSQ work vector: 5*N + 2*M + N*(N+1)/2,
  !> with N = parameters and M = data points.
  integer, parameter :: LWORK = 5*NPARM + 2*NDM + NPARM*(NPARM + 1)/2

end module dimensions
