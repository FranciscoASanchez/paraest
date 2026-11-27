!--------------------------------------------------------------------------------------------------------------------------
!
!  Module saturation: saturation pressure of a pure compound, for any model
!  of THERMO (module thermo_model).
!
!  psat_pure is that of GCA-cape-open (src/saturation.f90, which is based on
!  routine PSAT of paraest), with THERMO in place of lnphi_TP. Differences:
!    - P in atm and Z scalar per phase, as THERMO takes them. The algorithm
!      does not depend on the units.
!    - THERMO tells that the requested root does not exist (ic = -1) only
!      for some models (GCA). With RK-PR, both searches converge to the
!      same root. Then the root is classified with the critical point of
!      the model (Tc, Vc): T >= Tc, no L-V equilibrium; otherwise the root
!      is liquid if V < Vc (P too high), vapor if V > Vc (P too low).
!    - When P has been too high and too low, the next P is the midpoint (in
!      ln P) of the two, so steps of 0.5 in ln P cannot jump back and forth
!      over a narrow two-root region (near the critical point).
!    - Single root when |Z_V - Z_L| <= 1e-6 Z_V (1e-10 in the original):
!      the two searches of RK-PR can end at the same root with differences
!      of ~1e-9, and the iteration then "converged" to that trivial solution.
!    - Tolerance 1e-8 in ln phi^L - ln phi^V and in the step in ln P (1e-12
!      and 1e-13 in the original): the volume solver of RK-PR (VCALC) gives
!      ln phi of the liquid to ~1e-9, and the iteration never got below.
!
module saturation
  implicit none
  private
  public :: psat_pure

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Saturation pressure of component i (of the NC loaded in the model) at T:
!  ln phi^L(T,P) = ln phi^V(T,P).
!
!  Newton in ln P with exact derivative: for a pure compound,
!  d ln phi/d ln P = Z - 1, so d(ln phi^L - ln phi^V)/d ln P = Z_L - Z_V.
!  The step is limited to |d ln P| <= 1 only near the critical point
!  (Z_V - Z_L < 0.9); far from it, the full step is p^sat = f^L(T,P)/phi^V.
!  If the liquid (or vapor) root does not exist at the current P, P is
!  increased (or decreased) and the iteration goes on.
!
!  T in K, P in atm, Vc in cm3/mol (V = Z R T/P, R = 82.05 atm cm3/(mol K)).
!  P (inout): initial guess (> 0) and result.
!  Tc, Vc: critical point of the model for component i; 0 if unknown (then
!          a single root means no L-V equilibrium, as in the original).
!  info: 0 converged; 1 not converged; 2 no two distinct phases (T >= Tc
!        of the model, or trivial solution).
!  Z(2): Z of vapor and liquid at the solution.
!  nIt: number of iterations.
!
subroutine psat_pure(NC, i, T, P, Tc, Vc, info, Z, nIt)
  use thermo_model, only: model
  integer, intent(in) :: NC, i
  double precision, intent(in) :: T, Tc, Vc
  double precision, intent(inout) :: P
  integer, intent(out) :: info, nIt
  double precision, intent(out) :: Z(2)

  integer, parameter :: maxIt = 100, NCM = 30     ! NCM: as thermopack
  double precision, parameter :: R = 82.05D0
  double precision, parameter :: tol = 1D-8       ! in ln phi^L - ln phi^V and in ln P
  integer :: it, ICV, ICL, iGues, side
  double precision :: n(NCM), lnphiV(NCM), lnphiL(NCM), ZV, ZL, f, dlnP, lnP, lnPLow, lnPHigh
  double precision :: dLPhi(NCM,NCM), dLPhiT(NCM), dLPhiP(NCM), dHn(NCM), H, Cp, B

  n = 0
  n(i) = 1
  ZV = 0; ZL = 0
  Z = 0
  iGues = 0
  info = 1
  lnP = log(P)
  lnPLow = -huge(1D0)       ! highest ln P with vapor only
  lnPHigh = huge(1D0)       ! lowest ln P with liquid only
  if (Tc > 0 .and. T >= Tc) then
    info = 2
    nIt = 0
    return
  endif
  do it = 1, maxIt

    P = exp(lnP)
    call THERMO(model, NC, -1, 1, iGues, T, P, n, ZV, lnphiV, dLPhi, dLPhiT, dLPhiP, H, dHn, Cp, B, ICV)
    call THERMO(model, NC, +1, 1, iGues, T, P, n, ZL, lnphiL, dLPhi, dLPhiT, dLPhiP, H, dHn, Cp, B, ICL)
    iGues = 1

    ! Missing phase at this P (side: -1 vapor only, +1 liquid only)
    side = 0
    if (ICV == -1 .and. ICL == 1) then
      side = 1
    elseif (ICL == -1 .and. ICV == 1) then
      side = -1
    elseif (.not. (abs(ZV - ZL) > 1D-6*ZV)) then
      ! Single root (or NaN): classified with the critical volume
      if (Vc > 0 .and. ZL > 0) then
        side = 1
        if (ZL*R*T/P > Vc) side = -1
      else
        info = 2
        exit
      endif
    endif
    if (side /= 0) then
      if (side > 0) then
        lnPHigh = min(lnPHigh, lnP)
      else
        lnPLow = max(lnPLow, lnP)
      endif
      if (lnPHigh - lnPLow < 1D-12) then
        info = 2                            ! no two-root region left
        exit
      elseif (lnPHigh < huge(1D0) .and. lnPLow > -huge(1D0)) then
        lnP = (lnPLow + lnPHigh)/2
      else
        lnP = lnP - side*0.5D0
      endif
      iGues = 0
      ZV = 0; ZL = 0
      cycle
    endif

    f = lnphiL(i) - lnphiV(i)
    dlnP = -f/(ZL - ZV)
    ! Far from the critical point (Z_V - Z_L > 0.9) the full step is safe:
    ! with a nearly ideal vapor it is ln P_new = ln phi^L(T,P) + ln P - ln phi^V,
    ! i.e., p^sat = f^L(T,P)/phi^V up to the Poynting correction, which
    ! vanishes at the next iteration (contraction ~ Z_L/Z_V). Near the
    ! critical point, the step is limited to |d ln P| <= 1.
    if (ZV - ZL < 0.9D0 .and. abs(dlnP) > 1) dlnP = sign(1D0, dlnP)
    lnP = lnP + dlnP
    if (abs(f) < tol .or. abs(dlnP) < tol) then
      P = exp(lnP)
      info = 0
      exit
    endif

  enddo
  Z(1) = ZV
  Z(2) = ZL
  nIt = min(it, maxIt)

endsubroutine psat_pure

endmodule saturation
