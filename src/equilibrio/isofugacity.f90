!--------------------------------------------------------------------------------------------------------------------------
!
!  Departure from isofugacity at experimental conditions (no iteration in
!  T or P: the model is evaluated at the measured T, P and compositions).
!  Fugacities come from the thermodynamic model through ln_phi.
!
!--------------------------------------------------------------------------------------------------------------------------
!
!  Pure component, at its experimental vapor pressure P and temperature T:
!
!     apart = phi_liq/phi_vap - 1
!
!  The model must be loaded with the component alone (nc = 1).
!  ier = 1 if the liquid or the vapor root is not found, or if both roots
!  are the same (|zL - zV| < 1e-3 zV): then apart would be 0 for free.
!
subroutine isofug_pure(ng, nst, T, P, apart, zL, zV, ier)
  implicit none
  integer, intent(in) :: ng, nst
  double precision, intent(in) :: T, P
  double precision, intent(out) :: apart, zL, zV
  integer, intent(out) :: ier

  double precision :: x(1), lnPhiL(1), lnPhiV(1)
  integer :: ierL, ierV

  x(1) = 1D0
  call ln_phi(1, ng, nst, T, P, x, 1, lnPhiL, zL, ierL)
  call ln_phi(1, ng, nst, T, P, x, -1, lnPhiV, zV, ierV)
  apart = 0D0
  ier = 0
  if (ierL /= 0 .or. ierV /= 0 .or. dabs(zL - zV) < 1D-3*zV) then
    ier = 1
    return
  endif
  apart = dexp(lnPhiL(1) - lnPhiV(1)) - 1D0
endsubroutine isofug_pure

!--------------------------------------------------------------------------------------------------------------------------
!
!  Mixture: liquid of composition x and vapor of composition y, at T and P:
!
!     apart(i) = x(i) phi_liq(i) / (y(i) phi_vap(i)) - 1
!
!  The model must be loaded with these nc components.
!  If the vapor root is not found, the y phase is evaluated with the liquid
!  root (as the former paraest routine did: this allows liquid-liquid data).
!  ier = 1 if the liquid root is not found, or if a component is in the
!  liquid but not in the vapor (y(i) = 0, x(i) > 0).
!  A component absent from both phases gets apart(i) = 0.
!
subroutine isofug_mix(nc, ng, nst, T, P, x, y, apart, zL, zV, ier)
  implicit none
  integer, intent(in) :: nc, ng, nst
  double precision, intent(in) :: T, P, x(nc), y(nc)
  double precision, intent(out) :: apart(nc), zL, zV
  integer, intent(out) :: ier

  double precision :: lnPhiL(nc), lnPhiV(nc)
  integer :: i, ierL, ierV

  apart = 0D0
  call ln_phi(nc, ng, nst, T, P, x, 1, lnPhiL, zL, ierL)
  call ln_phi(nc, ng, nst, T, P, y, -1, lnPhiV, zV, ierV)
  if (ierV /= 0) call ln_phi(nc, ng, nst, T, P, y, 1, lnPhiV, zV, ierV)
  ier = 0
  if (ierL /= 0) then
    ier = 1
    return
  endif
  do i = 1, nc
    if (x(i) == 0D0 .and. y(i) == 0D0) cycle
    if (y(i) == 0D0) then           ! present in the liquid only: no finite departure
      apart = 0D0
      ier = 1
      return
    endif
    apart(i) = x(i)*dexp(lnPhiL(i) - lnPhiV(i))/y(i) - 1D0
  enddo
endsubroutine isofug_mix
