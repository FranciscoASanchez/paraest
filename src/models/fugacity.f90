!--------------------------------------------------------------------------------------------------------------------------
!
!  Link between paraest and the thermodynamic model: logarithm of the
!  fugacity coefficients of one phase at given T, P and composition.
!
!  The equilibrium routines written for paraest (isofugacity, ...) get
!  fugacities only through here, so using another model means changing only
!  this routine. It calls THERMO (thermopack) with the model of module
!  thermo_model (1 = GCA-EoS, 4 = RK-PR EoS).
!
!     phase = 1    liquid (dense) root
!            -1    vapor root
!     ier   = 0    the requested root was found
!             1    it was not found (lnPhi and Z are those of the other root)
!
!  T in K, P in atm. The model must already be loaded with nc components.
!
!--------------------------------------------------------------------------------------------------------------------------
subroutine ln_phi(nc, ng, nst, T, P, x, phase, lnPhi, Z, ier)
  use thermo_model, only: model
  implicit none
  integer, intent(in) :: nc, ng, nst, phase
  double precision, intent(in) :: T, P, x(nc)
  double precision, intent(out) :: lnPhi(nc), Z
  integer, intent(out) :: ier

  integer, parameter :: NCM = 30          ! must match thermopack
  double precision :: xn(NCM), phi(NCM), dLPhi(NCM,NCM), dLPhiT(NCM), dLPhiP(NCM), dHn(NCM), H, Cp, B
  integer :: ic

  xn = 0D0
  xn(:nc) = x
  Z = 0D0                                 ! no initial guess (guess = 0)
  call THERMO(model, nc, phase, 1, 0, T, P, xn, Z, phi, dLPhi, dLPhiT, dLPhiP, H, dHn, Cp, B, ic)
  lnPhi = phi(:nc)
  ier = 0
  if (ic /= 1) ier = 1
endsubroutine ln_phi
