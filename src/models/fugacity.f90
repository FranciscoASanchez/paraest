!--------------------------------------------------------------------------------------------------------------------------
!
!  Link between paraest and the thermodynamic model: logarithm of the
!  fugacity coefficients of one phase at given T, P and composition.
!
!  The equilibrium routines written for paraest (isofugacity, ...) get
!  fugacities only through here, so using another model means changing only
!  this routine. For now it calls GCEOS from the GCA-EoS package.
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
  implicit none
  integer, intent(in) :: nc, ng, nst, phase
  double precision, intent(in) :: T, P, x(nc)
  double precision, intent(out) :: lnPhi(nc), Z
  integer, intent(out) :: ier

  integer, parameter :: NCM = 30          ! must match the GCA package
  double precision :: xn(NCM), phi(NCM), dLPhi(NCM,NCM), dLPhiT(NCM), dLPhiP(NCM)
  integer :: ic

  call PARAGC(T, nc, ng, nst, 0)          ! temperature-dependent parameters
  xn = 0D0
  xn(:nc) = x
  ! First derivatives only (iDer = 1), no temperature derivatives (iTemp = 0),
  ! automatic initial guess (IGz = 0).
  call GCEOS(nc, ng, nst, 1, 0, T, P, xn, phi, dLPhi, dLPhiT, dLPhiP, Z, 0, phase, ic)
  lnPhi = phi(:nc)
  ier = 0
  if (ic /= 1) ier = 1
endsubroutine ln_phi
