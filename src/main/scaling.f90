!--------------------------------------------------------------------------------------------------------------------------
!
!  Module scaling: transformation between the model parameters and the
!  variables seen by the optimizer.
!
!    XXSCAL:   parameters -> scaled variables (before optimizing)
!    XXUNSCAL: scaled variables -> parameters (at each evaluation)
!
!  The type of each parameter is in IDXX (set by set_param_types):
!     1  g*  (g_ii)            2  k*  (k_ij)          3  g' (or k' if inIJ)
!     4  Q (area)              5  alpha               6  association energy
!     7  association volume    8  Tc                  9  dc
!    10  g"
!
!  Bounded parameters are scaled with fScld: scl*log((x - vMin)/(vMax - x)).
!
!--------------------------------------------------------------------------------------------------------------------------
module scaling
  use dimensions, only: NPARM
  implicit none
  private
  public :: XXSCAL, XXUNSCAL, fScld, fUnscld, search_box
  public :: inIJ, dc_mean, error_Dc

  ! Scaling factors (scl) of each parameter type.
  ! Note: since the step of ZXSSQ's numerical derivatives is almost fixed in
  ! the scaled variable (~4e-10), these factors set in practice that step in
  ! the transformed variable u (u = x for linear types, u = log(...) for
  ! bounded ones): du ~ 4e-10/scl.
  ! Study of 2026-11-22 (co2-phenols-fit): derivatives of k*, alpha, eps and
  ! kappa are exact for du between ~3e-5 and 1e-8; larger steps cross kinks
  ! of the residuals (with the old factors, du = 4e-3 to 4e-2, the gradient
  ! had the wrong sign near the end and the fit stopped with IER = 39).
  ! Those four use 4e-5 (du ~ 1e-5).
  ! Study of 2026-11-22 (n-alkanes): with fTemp = 1e-7 the derivatives of g'
  ! and g" had errors of 2e-3 to 1.5e-2 (du = 4e-3); with 4e-5, 1e-5. dc
  ! (fDc = 2e-7, du = 2e-3) reached the same minimum but its gradient was
  ! 25-75 % off near the end. Both use 4e-5. g* (du = 4 on ~4e5) is fine.
  ! Not studied yet: Q and Tc.
  double precision, parameter :: fGgg   = 1D-10
  double precision, parameter :: fRkij  = 4D-5
  double precision, parameter :: fTemp  = 4D-5
  double precision, parameter :: fQarea = 1D0
  double precision, parameter :: fAlfa  = 4D-5
  double precision, parameter :: fEnAss = 4D-5
  double precision, parameter :: fRkass = 4D-5
  double precision, parameter :: fTc    = 1D-8
  double precision, parameter :: fDc    = 4D-5

  ! Bounds of the bounded parameters.
  double precision, parameter :: vMinKstr = 0.2D0, vMaxKstr = 30D0   ! k*
  double precision, parameter :: vMaxKp   = 9D0                      ! k' between -vMaxKp and vMaxKp
  double precision, parameter :: vMinGp   = -2D0, vMaxGp = 2D0       ! g' between -vMaxGp and -vMinGp
  double precision, parameter :: vMaxAlfa = 29D0                     ! alpha between -vMaxAlfa and vMaxAlfa
  ! Association volume: between range_kappa(1) and range_kappa(2) (module
  ! config, option of paraest.cfg; default 0 and 20), for LM and PSO.

  ! Minimum relative distance to the bounds: the logarithm never gets 0.
  double precision, parameter :: margin = 1D-10

  ! inIJ(i): parameter i of type 3 is a k' (i-j interaction), not a g'.
  logical :: inIJ(NPARM) = .false.

  ! dc: each one moves within dc_mean(i) +- error_Dc(i) (set by read_input).
  double precision :: dc_mean(NPARM) = 0D0, error_Dc(NPARM) = 0D0

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Transforms the values of XX with the scaling factors.
!
!  boundFrac: fraction of the interval of each bounded parameter left free
!  next to the bounds (see fScld). 0 in prediction, 1e-3 in a fit.
!
subroutine XXSCAL(NPAR, XX, IDXX, boundFrac)
  use config, only: range_kappa
  integer, intent(in) :: NPAR, IDXX(NPAR)
  double precision, intent(inout) :: XX(NPAR)
  double precision, intent(in) :: boundFrac
  integer :: i

  do i = 1, NPAR
    select case (IDXX(i))
    case (1)
      XX(i) = fGgg*XX(i)
    case (2)
      XX(i) = fScld(XX(i), fRkij, vMinKstr, vMaxKstr, boundFrac)
    case (3)
      if (inIJ(i)) then
        XX(i) = fScld(XX(i), fTemp, -vMaxKp, vMaxKp, boundFrac)
      else
        XX(i) = fScld(XX(i), fTemp, -vMaxGp, -vMinGp, boundFrac)
      endif
    case (4)
      XX(i) = fQarea*XX(i)
    case (5)
      XX(i) = fScld(XX(i), fAlfa, -vMaxAlfa, vMaxAlfa, boundFrac)
    case (6)
      ! Logarithmic: it has only a lower bound (0). max avoids the log of 0
      ! or of a negative number. In a fit, at least 1 K (negligible against
      ! typical values of 1000-3000 K): from 0 it could not move.
      if (boundFrac > 0D0) then
        XX(i) = fEnAss*dlog(max(XX(i), 1D0))
      else
        XX(i) = fEnAss*dlog(max(XX(i), tiny(1D0)))
      endif
    case (7)
      XX(i) = fScld(XX(i), fRkass, range_kappa(1), range_kappa(2), boundFrac)
    case (8)
      XX(i) = fTc*XX(i)
    case (9)
      XX(i) = fScld(XX(i), fDc, dc_mean(i) - error_Dc(i), dc_mean(i) + error_Dc(i), boundFrac)
    case (10)
      XX(i) = fTemp*XX(i)       ! g": no bounds for now
    endselect
  enddo
endsubroutine XXSCAL

!--------------------------------------------------------------------------------------------------------------------------
!
!  Inverse transformation of XX, with the result in XXAUX. XX is not
!  changed: FUNC must not alter the values ZXSSQ passes to it.
!
subroutine XXUNSCAL(NPAR, XX, IDXX, XXAUX)
  use config, only: range_kappa
  integer, intent(in) :: NPAR, IDXX(NPAR)
  double precision, intent(in) :: XX(NPAR)
  double precision, intent(out) :: XXAUX(NPAR)
  integer :: i

  do i = 1, NPAR
    select case (IDXX(i))
    case (1)
      XXAUX(i) = XX(i)/fGgg
    case (2)
      XXAUX(i) = fUnscld(XX(i), fRkij, vMinKstr, vMaxKstr)
    case (3)
      if (inIJ(i)) then
        XXAUX(i) = fUnscld(XX(i), fTemp, -vMaxKp, vMaxKp)
      else
        XXAUX(i) = fUnscld(XX(i), fTemp, -vMaxGp, -vMinGp)
      endif
    case (4)
      XXAUX(i) = XX(i)/fQarea
    case (5)
      XXAUX(i) = fUnscld(XX(i), fAlfa, -vMaxAlfa, vMaxAlfa)
      ! An alpha right on the bound is moved away slightly (1e-8 relative).
      if (vMaxAlfa - dabs(XXAUX(i)) < 1D-8) XXAUX(i) = sign((1D0 - 1D-8)*vMaxAlfa, XXAUX(i))
    case (6)
      XXAUX(i) = dexp(min(XX(i)/fEnAss, log(huge(1D0))))   ! min: exp does not overflow
    case (7)
      XXAUX(i) = fUnscld(XX(i), fRkass, range_kappa(1), range_kappa(2))
      ! Same as alpha: not left right on the upper bound.
      if (dabs(range_kappa(2) - XXAUX(i)) < 1D-8*range_kappa(2)) XXAUX(i) = (1D0 - 1D-8)*range_kappa(2)
    case (8)
      XXAUX(i) = XX(i)/fTc
    case (9)
      XXAUX(i) = fUnscld(XX(i), fDc, dc_mean(i) - error_Dc(i), dc_mean(i) + error_Dc(i))
    case (10)
      XXAUX(i) = XX(i)/fTemp
    case default
      XXAUX(i) = XX(i)
    endselect
  enddo
endsubroutine XXUNSCAL

!--------------------------------------------------------------------------------------------------------------------------
!
!  Scaling of a bounded variable.
!
!  vMin and vMax are the minimum and maximum values x can take.
!  scl stretches or compresses the returned value (see the note on the
!  factors).
!
!  So that it never raises an exception (log of 0 or of a negative number),
!  x is first brought inside the interval, at a relative distance "margin"
!  from the bounds. So it does not matter whether the initial value is right
!  on a bound (for instance kappa = 0) or rounding leaves it on the bound or
!  outside.
!
!  frac is the fraction of the interval left free next to each bound, if it
!  is larger than that margin. A fit uses frac = 1e-3: on the bound the
!  derivative is almost zero and the optimizer would not move the
!  parameter. In prediction frac = 0 (the given value is kept).
!
double precision function fScld(x, scl, vMin, vMax, frac)
  double precision, intent(in) :: x, scl, vMin, vMax, frac
  double precision :: d, xin

  d = max(frac, margin)*(vMax - vMin)
  xin = min(max(x, vMin + d), vMax - d)
  fScld = scl*dlog((xin - vMin)/(vMax - xin))
endfunction fScld

!--------------------------------------------------------------------------------------------------------------------------
!
!  Unscaling: inverse of fScld.
!
!  It is computed in two ways depending on the sign of x/scl, so that the
!  exponential always has an argument <= 0: it never overflows (with a very
!  small scl, exp(x/scl) gave infinity, and then infinity/infinity = NaN).
!  The result is always between vMin and vMax.
!
!  fScld moves x away from the bounds by a relative margin; if the result
!  falls within that margin (twice, for rounding), the bound itself is
!  returned. So a value given right on the bound (for instance kappa = 0)
!  comes back unchanged, and the user does not see the margin.
!
double precision function fUnscld(x, scl, vMin, vMax)
  double precision, intent(in) :: x, scl, vMin, vMax
  double precision :: t, ex, d

  t = x/scl
  if (t >= 0D0) then
    ex = dexp(-t)
    fUnscld = (vMin*ex + vMax)/(1D0 + ex)
  else
    ex = dexp(t)
    fUnscld = (vMin + vMax*ex)/(1D0 + ex)
  endif

  d = 2D0*margin*(vMax - vMin)
  if (fUnscld - vMin <= d) fUnscld = vMin
  if (vMax - fUnscld <= d) fUnscld = vMax
endfunction fUnscld

!--------------------------------------------------------------------------------------------------------------------------
!
!  Search box of each parameter for PSO, in model (unscaled) units.
!
!  Parameters with physical bounds are searched between them, leaving free
!  the fraction boundFrac next to each bound (as XXSCAL does). The others,
!  between multiples of the initial value x0. k*, alpha, g' and g" between
!  absolute values (within the physical bounds), widened to include x0:
!  narrower than the bounds used by LM. The ranges are options of
!  paraest.cfg. The association energy keeps XXSCAL's lower
!  limit of 1 K.
!
subroutine search_box(NPAR, IDXX, x0, boundFrac, lo, hi)
  use config, only: range_gstr, range_kstr, range_alpha, range_gp, range_gdd, range_q, range_eps, range_tc, range_kappa
  integer, intent(in) :: NPAR, IDXX(NPAR)
  double precision, intent(in) :: x0(NPAR), boundFrac
  double precision, intent(out) :: lo(NPAR), hi(NPAR)
  integer :: i

  do i = 1, NPAR
    select case (IDXX(i))
    case (1)
      call times(range_gstr)
    case (2)
      call between(vMinKstr, vMaxKstr)
      call absolute(range_kstr)
    case (3)
      if (inIJ(i)) then
        call between(-vMaxKp, vMaxKp)
      else
        call between(-vMaxGp, -vMinGp)
        call absolute(range_gp)
      endif
    case (4)
      call times(range_q)
    case (5)
      call between(-vMaxAlfa, vMaxAlfa)
      call absolute(range_alpha)
    case (6)
      call times(range_eps)
      lo(i) = max(lo(i), 1D0)
      hi(i) = max(hi(i), 1D0)
    case (7)
      call between(range_kappa(1), range_kappa(2))
    case (8)
      call times(range_tc)
    case (9)
      call between(dc_mean(i) - error_Dc(i), dc_mean(i) + error_Dc(i))
    case (10)
      lo(i) = min(range_gdd(1), x0(i))
      hi(i) = max(range_gdd(2), x0(i))
    case default
      lo(i) = x0(i)
      hi(i) = x0(i)
    endselect
  enddo

contains

  subroutine between(vMin, vMax)
    double precision, intent(in) :: vMin, vMax
    double precision :: d
    d = max(boundFrac, margin)*(vMax - vMin)
    lo(i) = vMin + d
    hi(i) = vMax - d
  endsubroutine between

  ! Narrows the physical bounds (already in lo, hi) to the absolute range
  ! r, widened to include x0.
  subroutine absolute(r)
    double precision, intent(in) :: r(2)
    lo(i) = max(lo(i), min(r(1), x0(i)))
    hi(i) = min(hi(i), max(r(2), x0(i)))
  endsubroutine absolute

  subroutine times(r)
    double precision, intent(in) :: r(2)
    lo(i) = min(r(1)*x0(i), r(2)*x0(i))
    hi(i) = max(r(1)*x0(i), r(2)*x0(i))
  endsubroutine times

endsubroutine search_box

endmodule scaling
