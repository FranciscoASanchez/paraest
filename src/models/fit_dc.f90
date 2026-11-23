!--------------------------------------------------------------------------------------------------------------------------
!
!  Fits the critical diameter of component ic so that the GCA-EoS reproduces
!  exactly (|P - pFit| < 1e-8) the vapor pressure pFit at temperature T.
!
!  The Newton loop is the one in PARMOL (GCA-v1.13.6.2), with two changes:
!    - it works on the whole mixture, calling PSAT with IC = ic, instead of
!      moving the component to position 1 (tested: same result, same speed);
!    - it has an iteration limit and returns ier instead of stopping.
!
!  The model must already be loaded (PARGR and PARMOL, or fakePARGR and
!  fakePARMOL in paraest) with nc components.
!
!  On exit, dcFit holds the fitted dc and DC(ic) in COMMON /MOL/ is updated.
!  If the fit fails, dcFit and DC(ic) keep their initial value.
!
!     ier = 0   converged
!           1   no convergence in maxIter iterations
!           2   zero derivative dP/dd
!           3   PSAT failed (before, its error was ignored: PSAT leaves P
!               as it was, so F = 0 and a meaningless d was accepted)
!
!--------------------------------------------------------------------------------------------------------------------------
subroutine fit_dc(ic, nc, ng, nst, T, pFit, dcFit, ier)
  implicit none
  integer, intent(in) :: ic, nc, ng, nst
  double precision, intent(in) :: T, pFit
  double precision, intent(out) :: dcFit
  integer, intent(out) :: ier

  integer, parameter :: NCM = 30          ! must match the GCA package
  integer, parameter :: maxIter = 100
  double precision, parameter :: tol = 1D-8, dd = 1D-4

  double precision :: DC, D, DT, HA, HB, TC, PC
  COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
  COMMON/CRIT/TC(NCM),PC(NCM)

  integer :: iter, IGP, IGZ, ierPsat
  double precision :: Z(2), P1, P2, F1, F2, D1, DFDD, DPDD, DEL, DLC, fCrt

  dcFit = DC(ic)
  ier = 0

  call PARAGC(T, nc, ng, nst, 0)          ! D(:) at temperature T
  IGP = 1
  IGZ = 0
  P1 = pFit
  Z = 0D0
  fCrt = 1D0
  iter = 0
  do while (fCrt >= tol)
    iter = iter + 1
    if (iter > maxIter) then
      ier = 1
      return
    endif
    D1 = D(ic)
    call PSAT(T, P1, Z, nc, ng, nst, ic, IGP, IGZ, ierPsat)
    if (ierPsat < 0) then
      ier = 3
      DC(ic) = dcFit
      return
    endif
    if (iter == 3) IGZ = 5
    F1 = P1 - pFit
    fCrt = dabs(F1)
    ! Numerical derivative with respect to d
    P2 = P1
    D(ic) = D1 + dd
    call PSAT(T, P2, Z, nc, ng, nst, ic, IGP, IGZ, ierPsat)
    if (ierPsat < 0) then
      ier = 3
      DC(ic) = dcFit
      return
    endif
    F2 = P2 - pFit
    DFDD = (F2 - F1)/dd
    DPDD = (P2 - P1)/dd
    if (DFDD == 0D0) then
      ier = 2
      return
    endif
    ! Newton step, limited to 5% of d
    DEL = -F1/DFDD
    DLC = DEL/D(ic)
    if (DLC > .05D0) DEL = .05D0*D(ic)
    if (DLC < -.05D0) DEL = -.05D0*D(ic)
    D(ic) = D1 + DEL
    ! Initial guess for the next iteration
    P2 = P1 + DEL*DPDD
    Z(2) = Z(2)/P1*P2
    P1 = P2
  enddo

  ! d at T -> dc
  DC(ic) = D(ic)/(1D0 - .12D0*dexp(-2D0*TC(ic)/3D0/T))/1.065655D0
  dcFit = DC(ic)
endsubroutine fit_dc
