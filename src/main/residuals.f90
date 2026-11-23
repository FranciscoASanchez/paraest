!--------------------------------------------------------------------------------------------------------------------------
!
!  Module residuals: the residuals of the objective function.
!
!  FUNC is the routine ZXSSQ calls at each evaluation. It sets the model
!  parameters from the optimizer variables, fits the critical diameters,
!  and computes one residual per experimental point:
!
!     pure components        vapor pressure (IOBJSAT >= 80) or isofugacity (50)
!     data sets (IHEN)       0 phase equilibrium     1 excess properties
!                            2 solid solubility      3 critical points
!                            5 activity coefficients
!
!  Besides the residuals (FCION), each routine stores the calculated values
!  (PCALC, XCALC, ...) and the errors (FCAUX, ...) used later by
!  write_results.
!
!--------------------------------------------------------------------------------------------------------------------------
module residuals
  use dimensions
  use fit_data
  use gca_parameters, PSATu => PSAT   ! PSAT is also a subroutine of the GCA package
  use gca_bridge, only: fakePARGR, fakePARMOL, negative_cohesion
  use config, only: PENALIZ, PENALIZNAN
  implicit none
  private
  public :: FUNC, write_failed_points
  public :: nFailed, nNaN, nFailedFirst, nNaNFirst     ! read by the PSO objective

  ! Arrays passed to the GCA package: dimensioned with its maximum number of components.
  integer, parameter :: NCMgca = 30

  double precision, parameter :: Rgas_SI = 8.314D0      ! J/(mol K)

  ! Set by fakePARGR: number of association sites.
  integer :: NST

  ! Points that could not be calculated (failed, or NaN/Inf residual), in the
  ! last evaluation and in the first one (the initial parameters).
  integer, protected :: nFailed = 0, nNaN = 0, nFailedFirst = 0, nNaNFirst = 0
  integer :: nEval = 0

  ! dc of each component at the first evaluation (the user's values): the dc
  ! fit of fit_critical_diameters always starts from them.
  double precision :: dcStart(NCM) = 0D0

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Residuals of the objective function (called by ZXSSQ). FUNC must not
!  alter XX, NCXP or NPAR (see the IMSL documentation).
!
subroutine FUNC(XX, NCXP, NPAR, FCION)
  use scaling, only: XXUNSCAL
  use utilities, only: rep
  integer, intent(in) :: NCXP, NPAR
  double precision, intent(in) :: XX(NPAR)
  double precision, intent(inout) :: FCION(NCXP)

  double precision :: XXAUX(NPARM)      ! unscaled parameters
  double precision :: SUMQ
  integer :: numExp, ks, NCw, i

  nFailed = 0
  ! Model parameters from the unscaled optimizer variables
  call XXUNSCAL(NPAR, XX, IDXX, XXAUX)
  call set_model_params(XXAUX, TEMPCRIT)

  call fit_critical_diameters()

  ! numExp counts the experimental points already processed
  numExp = 0
  call pure_component_residuals(FCION, numExp)

  do ks = 1, IGLE
    call load_set_components(ks, NCw)
    select case (IHEN(ks))
    case (0)
      call phase_equilibrium_residuals(ks, NCw, FCION, numExp)
    case (1)
      call excess_property_residuals(ks, FCION, numExp)
    case (2)
      call sle_residuals(ks, NCw, FCION, numExp)
    case (3)
      call critical_point_residuals(ks, NCw, FCION, numExp)
    case (5)
      call activity_coefficient_residuals(ks, NCw, FCION, numExp)
    endselect
  enddo

  ! NaN guard (added 2026-11-06): a NaN or Inf residual is replaced by
  ! PENALIZNAN, so it does not spoil the objective function.
  nNaN = 0
  do i = 1, NCXP
    if (FCION(i) /= FCION(i) .or. dabs(FCION(i)) > 1.D30) then
      FCION(i) = PENALIZNAN
      nNaN = nNaN + 1
    endif
  enddo
  if (nNaN > 0) write(*,'(A,I6)') ' WARNING: non-finite FCION at points: ', nNaN

  nEval = nEval + 1
  if (nEval == 1) then
    nFailedFirst = nFailed
    nNaNFirst = nNaN
  endif

  ! Course of the optimization: unscaled parameters and objective function
  SUMQ = sum(FCION(:NCXP)**2)
  open(unit=1, status='OLD', file='paraest.out', form='FORMATTED', position='APPEND')
  if (NPAR > 0) then
    call XXUNSCAL(NPAR, XX, IDXX, XXAUX)
    write(1,'('//rep(NPAR,"F12.4,X")//'E12.5)') (XXAUX(i), i = 1, NPAR), SUMQ
  endif
  close(unit=1)
endsubroutine FUNC

!--------------------------------------------------------------------------------------------------------------------------
!
!  Writes to paraest.out (and to the screen) how many points could not be
!  calculated: in the tables their calculated values are zeros. In a fit,
!  at the initial and at the final parameters; with PENALIZ = 0 a failed
!  point adds nothing to SSQ, so a fit that ends with more failures than it
!  started with may have lowered SSQ just by losing points.
!
subroutine write_failed_points(fit)
  logical, intent(in) :: fit
  character(*), parameter :: fmt = '(X,A,I6)'
  integer :: u, k

  open(unit=1, status='OLD', file='paraest.out', form='FORMATTED', position='APPEND')
  do k = 1, 2
    u = merge(1, 6, k == 1)
    write(u,*)
    if (fit) then
      write(u,fmt) 'Points not calculated (zeros in the tables) at the initial parameters:', nFailedFirst
      write(u,fmt) 'Points not calculated (zeros in the tables) at the final parameters:  ', nFailed
      if (nNaNFirst + nNaN > 0) then
        write(u,fmt) 'Residuals NaN/Inf (replaced by PENALIZNAN) at the initial parameters: ', nNaNFirst
        write(u,fmt) 'Residuals NaN/Inf (replaced by PENALIZNAN) at the final parameters:   ', nNaN
      endif
      if (nFailed > nFailedFirst .and. PENALIZ == 0D0) then
        write(u,'(X,A,I0,A)') 'WARNING: the fit ended with ', nFailed - nFailedFirst, &
          ' more points not calculated than it started with; with PENALIZ = 0 they'
        write(u,'(X,A)') 'add nothing to SSQ, so part of its decrease may come from losing points.'
      endif
    else
      write(u,fmt) 'Points not calculated (zeros in the tables):', nFailed
      if (nNaN > 0) write(u,fmt) 'Residuals NaN/Inf (replaced by PENALIZNAN):', nNaN
    endif
  enddo
  close(unit=1)
endsubroutine write_failed_points

!--------------------------------------------------------------------------------------------------------------------------
!
!  Critical diameters: each dc that is not a fitted parameter is adjusted to
!  reproduce exactly the vapor pressure point (TSAT, PSAT) of its component.
!
!  The Newton iteration always starts from the user's dc (dcStart), not from
!  the result of the previous evaluation: otherwise the result depended on
!  the history. A PSO particle with extreme parameters took dc to another
!  root, and the following evaluations, even at the initial parameters,
!  started from there (n-alkanes-fit-g with pso+lm: SSQ 1.06 at the initial
!  values before the PSO, 1.8e10 after it). If the fit fails, dc keeps the
!  user's value.
!
subroutine fit_critical_diameters()
  integer :: i, iErrDc
  double precision :: dcFit, wDC(NCM)

  if (.not. any(TSAT(:NC) /= 0)) return

  if (nEval == 0) dcStart(:NC) = DC(:NC)
  wDC(:NC) = DC(:NC)
  do i = 1, NC
    if (TSAT(i) /= 0 .and. .not. any(IDDC(:numDC) == i)) wDC(i) = dcStart(i)
  enddo

  call fakePARGR(NC, NG, NGA, NST, NU(:NG,:NC), QAREA(:NG), TSTR(:NG), GSTR(:NG), GDOT(:NG), GDDOT(:NG), &
                 RKIJ(:NG,:NG), RKTIJ(:NG,:NG), ALFAIJ(:NG,:NG), ONUASOC(:NGA,:NC), MASSOC(:NGA),             &
                 ENASS(:2,:NGA,:2,:NGA), RKASS(:2,:NGA,:2,:NGA))
  call fakePARMOL(NC, NG, NST, wDC(:NC), TEMPCRIT(:NC), PRCRIT(:NC))
  do i = 1, NC
    if (TSAT(i) == 0) cycle
    if (any(IDDC(:numDC) == i)) cycle
    call fit_dc(i, NC, NG, NST, TSAT(i), PSATu(i), dcFit, iErrDc)
    if (iErrDc == 0) then
      DC(i) = dcFit
    else
      DC(i) = dcStart(i)
      write(6,'(" dc fit failed for component ",I2,", ier = ",I2)') i, iErrDc
    endif
    write(6,'(" End of dc fit, component No. ",I2)') i
  enddo
endsubroutine fit_critical_diameters

!--------------------------------------------------------------------------------------------------------------------------
!
!  Pure components: vapor pressure (IOBJSAT >= 80; 81 adds the liquid
!  density) or isofugacity (IOBJSAT = 50). Below the melting point the
!  vapor pressure is turned into a sublimation pressure.
!
subroutine pure_component_residuals(FCION, numExp)
  double precision, intent(inout) :: FCION(*)
  integer, intent(inout) :: numExp

  double precision :: WDC(1), WTEMPCRIT(1), WPRCRIT(1), ONUASOCW(12,1), wZ(2)
  integer :: NUW(NGM,1)
  double precision :: feed(NCMgca), lnPhi(NCMgca), dLnPhi(NCMgca,NCMgca), dLnPhiT(NCMgca), dLnPhiP(NCMgca)
  double precision :: wTemp, wPsat, APART, WZFASL, WZFASV, HRV, HRL, TR, Yu, DHSL
  integer :: i, j, n, IERRYAF, iGP, iGz, iErrPSat, iComp, iErr

  do i = 1, NC
    if (NSAT(i) == 0) cycle
    ! The component alone
    WDC(1) = DC(i)
    WTEMPCRIT(1) = TEMPCRIT(i)
    WPRCRIT(1) = PRCRIT(i)
    NUW(:NG,1) = NU(:NG,i)
    ONUASOCW(:NGA,1) = ONUASOC(:NGA,i)

    do j = 1, NSAT(i)
      n = numExp + 1
      wTemp = TEXP(n)
      wPsat = PEXP(n)

      if (IOBJSAT(i) == 50) then
        ! Objective function: isofugacity
        write(6,'("Point ",I3," of component",I3)') j, i
        call load_component(i, WDC, WTEMPCRIT, WPRCRIT, NUW, ONUASOCW)
        if (negative_cohesion(wTemp, 1, NG)) then
          IERRYAF = 1
        else
          call isofug_pure(NG, NST, wTemp, wPsat, APART, WZFASL, WZFASV, IERRYAF)
        endif
        if (IERRYAF == 0) then
          PCALC(n) = wPsat
          ZLCALC(n) = WZFASL
          ZVCALC(n) = WZFASV
          FCAUX(n) = APART
          FCION(n) = W1 * FCAUX(n) * wSat(n)
        else
          PCALC(n) = 0.
          ZLCALC(n) = 0.
          ZVCALC(n) = 0.
          FCAUX(n) = 0.
          FCION(n) = PENALIZ
          nFailed = nFailed + 1
        endif
      endif

      if (IOBJSAT(i) >= 80) then
        ! Objective function: vapor pressure (+ density if IOBJSAT = 81)
        write(6,'("Point ",I3," of component",I3)') j, i
        call load_component(i, WDC, WTEMPCRIT, WPRCRIT, NUW, ONUASOCW)
        call paraGC(wTemp, 1, NG, NST, 1)
        iGP = 1
        iGz = 0       ! could use the result of the dc fit...
        iErrPSat = 0
        iComp = 1
        if (negative_cohesion(wTemp, 1, NG)) then
          iErrPSat = -1
        elseif (wPsat > 1E-10) then
          ! Standard vapor pressure calculation
          call PSAT(wTemp, wPsat, wZ, 1, NG, NST, iComp, iGP, iGz, iErrPSat)
        else
          ! psat(T) = phiL(T)*(1 atm)
          feed(1) = 1
          feed(2:NCMgca) = 0
          wPsat = 1
          call GCEOS(1, NG, NST, 1, 0, wTemp, wPsat, feed, lnPhi, dLnPhi, dLnPhiT, dLnPhiP, wZ(2), 0, 1, iErr)
          wPsat = dexp(lnPhi(1))
          wZ(1) = 1
          wZ(2) = wZ(2)*wPsat
        endif
        if (iErrPSat >= 0) then
          ! Residual enthalpies of vapor and liquid, for the vaporization enthalpy
          feed(1) = 1
          feed(2:NCMgca) = 0
          call NTLPY(1, NG, NST, wTemp, wPsat, feed, HRV, wZ(1), -1, -1, iErr)
          call NTLPY(1, NG, NST, wTemp, wPsat, feed, HRL, wZ(2), 1, 1, iErr)
        endif
        if (iErrPSat >= 0 .and. wTemp < Tfu(i)) then
          ! Sublimation pressure
          TR = (Tfu(i)/wTemp - 1)
          Yu = -DHfu(i)/Rgas_SI/Tfu(i) * TR
          DHSL = DHfu(i)
          wZ(2) = wZ(2) - wPsat*DVfu(i)/82.05/wTemp
          if (DCpfu(i) > 0) then
            ! Change in heat capacity
            Yu = Yu + DCpfu(i)/Rgas_SI * (TR + dlog(wTemp/Tfu(i)))
            DHSL = DHfu(i) + DCpfu(i)*(wTemp - Tfu(i))
          endif
          if (TtrSS(i) > 0 .and. DHss(i) > 0) then
            ! Solid-solid transitions
            if (wTemp < TtrSS(i)) then
              Yu = Yu - DHss(i)/Rgas_SI/TtrSS(i) * (TtrSS(i)/wTemp - 1)
              DHSL = DHSL + DHss(i)
            endif
          endif
          wPsat = wPsat*dexp(Yu)
          HRL = HRL - DHSL/.101325D0    ! HRL is in atm cm3, and it is a solid, not a liquid
        endif

        if (iErrPSat >= 0) then
          PCALC(n) = wPsat
          ZLCALC(n) = wZ(2)
          ZVCALC(n) = wZ(1)
          DHVL(n) = HRV - HRL
          FCAUX(n) = (PCALC(n)/PEXP(n) - 1)
          FCION(n) = W1 * FCAUX(n) * wSat(n)*PEXP(n)
          if (IOBJSAT(i) == 81) FCION(n) = FCION(n) + abs(wPsat/82.05D0/wTemp/wZ(2)/rhoExp(n) - 1)
        else
          PCALC(n) = 0.
          ZLCALC(n) = 0.
          ZVCALC(n) = 0.
          FCAUX(n) = 0.
          FCION(n) = PENALIZ
          nFailed = nFailed + 1
          DHVL(n) = 0
        endif
      endif

      numExp = numExp + 1
    enddo
  enddo
endsubroutine pure_component_residuals

!--------------------------------------------------------------------------------------------------------------------------
!
!  Loads the model with component i alone (fakePARGR and fakePARMOL replace
!  PARGR and PARMOL of the GCA package).
!
subroutine load_component(i, WDC, WTEMPCRIT, WPRCRIT, NUW, ONUASOCW)
  integer, intent(in) :: i
  double precision, intent(in) :: WDC(1), WTEMPCRIT(1), WPRCRIT(1), ONUASOCW(12,1)
  integer, intent(in) :: NUW(NGM,1)

  call fakePARGR(1, NG, NGA, NST, NUW(:NG,1), QAREA(:NG), TSTR(:NG), GSTR(:NG), GDOT(:NG), GDDOT(:NG), &
                 RKIJ(:NG,:NG), RKTIJ(:NG,:NG), ALFAIJ(:NG,:NG), ONUASOCW(:NGA,1), MASSOC(:NGA),           &
                 ENASS(:2,:NGA,:2,:NGA), RKASS(:2,:NGA,:2,:NGA))
  call fakePARMOL(1, NG, NST, WDC, WTEMPCRIT, WPRCRIT)
endsubroutine load_component

!--------------------------------------------------------------------------------------------------------------------------
!
!  Loads the model with the components of data set ks: 2, or 3 for ternary
!  data (NCw).
!
subroutine load_set_components(ks, NCw)
  integer, intent(in) :: ks
  integer, intent(out) :: NCw

  double precision :: WWDC(3), WWTEMPCRIT(3), WWPRCRIT(3), ONUASOCWW(12,3)
  integer :: NUWW(NGM,3), k, ic

  ! mez is read only for some data types (0 for the others)
  NCw = 2
  if (mez(ks) == 3) NCw = 3
  do k = 1, NCw
    ic = IDC(ks,k)
    WWDC(k) = DC(ic)
    WWTEMPCRIT(k) = TEMPCRIT(ic)
    WWPRCRIT(k) = PRCRIT(ic)
    NUWW(:NG,k) = NU(:NG,ic)
    ONUASOCWW(:NGA,k) = ONUASOC(:NGA,ic)
  enddo

  call fakePARGR(NCw, NG, NGA, NST, NUWW(:NG,:NCw), QAREA(:NG), TSTR(:NG), GSTR(:NG), GDOT(:NG), GDDOT(:NG), &
                 RKIJ(:NG,:NG), RKTIJ(:NG,:NG), ALFAIJ(:NG,:NG), ONUASOCWW(:NGA,:NCw), MASSOC(:NGA),             &
                 ENASS(:2,:NGA,:2,:NGA), RKASS(:2,:NGA,:2,:NGA))
  call fakePARMOL(NCw, NG, NST, WWDC(:NCw), WWTEMPCRIT(:NCw), WWPRCRIT(:NCw))
endsubroutine load_set_components

!--------------------------------------------------------------------------------------------------------------------------
!
!  Phase equilibrium data (IHEN = 0). Depending on IFL:
!    -1, -2  isofugacity at the experimental T, P, x, y (with -2 the
!            departure of component 2 is not used in the residual)
!     0      binary flash at T and P
!     1, 2   bubble pressure (T, x) or bubble temperature (P, x)
!     3      ternary flash at T and P
!
!  Errors are absolute, divided by the experimental uncertainty through the
!  weights wBin and wTer (Francisco, 07/08/2013).
!
subroutine phase_equilibrium_residuals(ks, NCw, FCION, numExp)
  integer, intent(in) :: ks, NCw
  double precision, intent(inout) :: FCION(*)
  integer, intent(inout) :: numExp

  double precision :: feed(NCMgca), zFrac(NCMgca,4), zVal(4)
  double precision :: xIsof(2), yIsof(2), apartIsof(2)
  double precision :: TEMP, PRESS, X1, Y1, X2, Y2, ZLL, ZVV, zLiquido, zVapor, ZLiq, ZVap, sq
  integer :: indF(4), iVap(1)
  integer :: i, n, nFailure, iErr, iOrder, iType
  logical :: badG

  nFailure = NOBS(ks)
  STDVA(ks) = 0D0
  STDVB(ks) = 0
  do i = 1, NOBS(ks)
    n = numExp + 1
    TEMP = TEXP(n)
    PRESS = PEXP(n)
    X1 = XEXP(n)
    Y1 = YEXP(n)
    ! g(T) <= 0 at the experimental T: not calculated (for bubble
    ! temperatures, IFL = 2, the check is at the experimental T)
    badG = negative_cohesion(TEMP, NCw, NG)

    if (IFL(ks) == -1 .or. IFL(ks) == -2) then
      ! Isofugacity
      write(6,*) "i=", i
      xIsof(1) = X1
      xIsof(2) = 1D0 - X1
      yIsof(1) = Y1
      yIsof(2) = 1D0 - Y1
      if (badG) then
        iErr = 1
      else
        call isofug_mix(2, NG, NST, TEMP, PRESS, xIsof, yIsof, apartIsof, ZLL, ZVV, iErr)
      endif
      if (iErr == 0) then
        ZLCALC(n) = ZLL
        ZVCALC(n) = ZVV
        FCION(n) = (W4**2.)*(2.+IFL(ks))*(apartIsof(2)**2) + (apartIsof(1)**2.)
        FCION(n) = dsqrt(FCION(n))
        FCAUX(n) = apartIsof(1)
        FCAUX2(n) = apartIsof(2)
      else
        ! A failed point is not counted for the standard deviation
        nFailure = nFailure - 1
        ZLCALC(n) = 0.
        ZVCALC(n) = 0.
        FCION(n) = PENALIZ
        nFailed = nFailed + 1
        FCAUX(n) = 0.
        FCAUX2(n) = 0.
      endif
      STDVA(ks) = STDVA(ks) + FCAUX(n)**2
      STDVB(ks) = STDVB(ks) + FCAUX2(n)**2

    elseif (IFL(ks) >= 0) then
      iErr = 1
      select case (IFL(ks))

      case (0)
        ! Binary multiphase flash. Phase 1 is liquid; phase 2 can be vapor or
        ! liquid, depending on the stability test. The phases are given in the
        ! same order as the experimental data.
        write(6,*) " i = ", i
        feed(1) = (X1 + Y1)/2
        feed(2) = 1D0 - feed(1)
        indF = 0
        zFrac = 0
        zVal = 0
        if (Y1 - X1 < 0) then
          iOrder = 1
        else
          iOrder = 2
        endif
        iErr = 0
        if (badG) then
          iErr = 1
        else
          call GC3Flash(NCw, NG, NST, TEMP, PRESS, feed(:NCMgca), zFrac, zVal, indF, iErr)
        endif
        if (zFrac(1,1) > zFrac(1,2)) then
          X1 = zFrac(1,iOrder)
          zLiquido = zVal(iOrder)
          Y1 = zFrac(1,3-iOrder)
          zVapor = zVal(3-iOrder)
        else
          X1 = zFrac(1,3-iOrder)
          zLiquido = zVal(3-iOrder)
          Y1 = zFrac(1,iOrder)
          zVapor = zVal(iOrder)
        endif
        if (iErr == 0) then
          ZLCALC(n) = zLiquido
          ZVCALC(n) = zVapor
          XCALC(n) = X1
          YCALC(n) = Y1
          FCAUX(n) = abs(XCALC(n) - XEXP(n))
        else
          nFailure = nFailure - 1
          ZLCALC(n) = 0.
          ZVCALC(n) = 0.
          XCALC(n) = 0.
          YCALC(n) = 0.
          FCAUX(n) = 0.
        endif

      case (1, 2)
        ! Bubble point. iType for GCADB: 3 = bubble pressure (IFL = 1),
        ! 1 = bubble temperature (IFL = 2).
        write(6,*) "i = ", i
        if (IFL(ks) == 1) then
          iType = 3
        else
          iType = 1
        endif
        iErr = 0
        if (badG) then
          iErr = 1
        else
          call GCADB(1, NCw, NG, NST, iType, TEMP, PRESS, X1, Y1, iErr, ZLiq, ZVap)
        endif
        if (iErr == 0) then
          ZLCALC(n) = ZLiq
          ZVCALC(n) = ZVap
          YCALC(n) = Y1
          if (iType == 3) then
            PCALC(n) = PRESS
            FCAUX(n) = abs(PCALC(n) - PEXP(n))
          else
            PCALC(n) = TEMP          ! the calculated bubble T is kept in PCALC
            FCAUX(n) = abs(PCALC(n) - TEXP(n))
          endif
        else
          nFailure = nFailure - 1
          ZLCALC(n) = 0.
          ZVCALC(n) = 0.
          YCALC(n) = 0.
          PCALC(n) = 0.
          FCAUX(n) = 0.
        endif

      case (3)
        ! Ternary multiphase flash
        write(6,*) " i = ", i
        feed(1) = (XEXP(n) + YEXP(n))/2
        feed(2) = (xExp2(n) + yExp2(n))/2
        feed(3) = 1D0 - feed(1) - feed(2)
        if (YEXP(n) - XEXP(n) < 0) then
          iOrder = 1
        else
          iOrder = 2
        endif
        iErr = 0
        if (badG) then
          iErr = 1
        else
          call GC3Flash(NCw, NG, NST, TEMP, PRESS, feed(:NCMgca), zFrac, zVal, indF, iErr)
        endif
        if (minval(indF(:4)) > -1) then
          ! LLE: phases in the same order as the experimental data
          if (zFrac(1,1) > zFrac(1,2)) then
            X1 = zFrac(1,iOrder)
            X2 = zFrac(2,iOrder)
            zLiquido = zVal(iOrder)
            Y1 = zFrac(1,3-iOrder)
            Y2 = zFrac(2,3-iOrder)
            zVapor = zVal(3-iOrder)
          else
            X1 = zFrac(1,3-iOrder)
            X2 = zFrac(2,3-iOrder)
            zLiquido = zVal(3-iOrder)
            Y1 = zFrac(1,iOrder)
            Y2 = zFrac(2,iOrder)
            zVapor = zVal(iOrder)
          endif
        else
          ! VLE: minval(indF) = -1, since there must be a vapor phase;
          ! minloc gives the position of the -1.
          iVap(1) = minloc(indF(:4),1)
          Y1 = zFrac(1,iVap(1))
          Y2 = zFrac(2,iVap(1))
          X1 = zFrac(1,3-iVap(1))
          X2 = zFrac(2,3-iVap(1))
          zLiquido = zVal(3-iVap(1))
          zVapor = zVal(iVap(1))
        endif
        if (iErr == 0) then
          ZLCALC(n) = zLiquido
          ZVCALC(n) = zVapor
          XCALC(n) = X1
          xCalc2(n) = X2
          YCALC(n) = Y1
          yCalc2(n) = Y2
          FCAUX(n) = abs(XCALC(n) - XEXP(n))
        else
          nFailure = nFailure - 1
          ZLCALC(n) = 0.
          ZVCALC(n) = 0.
          XCALC(n) = 0.
          xCalc2(n) = 0.
          YCALC(n) = 0.
          yCalc2(n) = 0.
          FCAUX(n) = 0.
        endif
      endselect

      ! Residual. An experimental composition equal to 0 means it was not
      ! measured: its term is left out.
      if (iErr == 0) then
        sq = 0D0
        FCAUX2(n) = 0D0
        select case (IFL(ks))
        case (0)
          if (YEXP(n) > 0D0) then
            sq = W4**2*(YCALC(n) - YEXP(n))**2*wBin(2,n)**2
            FCAUX2(n) = abs(YCALC(n) - YEXP(n))
          endif
          if (XEXP(n) > 0D0) sq = sq + W5**2*(XCALC(n) - XEXP(n))**2*wBin(1,n)**2
        case (1)
          ! Bubble pressure and y
          if (YEXP(n) > 0D0) then
            sq = (W4**2)*(YCALC(n) - YEXP(n))**2*wBin(2,n)**2
            FCAUX2(n) = abs(YCALC(n) - YEXP(n))
          endif
          sq = sq + W2**2*FCAUX(n)**2*wBin(1,n)**2
        case (2)
          ! Bubble temperature and y
          if (YEXP(n) > 0D0) then
            sq = (W4**2)*(YCALC(n) - YEXP(n))**2*wBin(2,n)**2
            FCAUX2(n) = abs(YCALC(n) - YEXP(n))
          endif
          sq = sq + W3**2*FCAUX(n)**2*wBin(1,n)**2
        case (3)
          sq = W5**2*(FCAUX(n))**2*wBin(1,n)**2 &
               + W6**2*(xCalc2(n) - xExp2(n))**2*wTer(1,n)**2
          if (YEXP(n) > 0D0) then
            FCAUX2(n) = abs(YCALC(n) - YEXP(n))
            sq = sq + W4**2*(FCAUX2(n))**2*wBin(2,n)**2 &
                    + W7**2*(yCalc2(n) - yExp2(n))**2*wTer(2,n)**2
          endif
        endselect
        FCION(n) = sqrt(sq)
      else
        FCION(n) = PENALIZ*wBin(1,n)**2
        nFailed = nFailed + 1
        FCAUX2(n) = 0.D0
      endif
      STDVA(ks) = STDVA(ks) + FCAUX(n)
      STDVB(ks) = STDVB(ks) + FCAUX2(n)
    endif

    numExp = numExp + 1
  enddo

  ! Average errors (isofugacity: average of the squares)
  if (nFailure /= 0) then
    STDVA(ks) = STDVA(ks)/nFailure
    STDVB(ks) = STDVB(ks)/nFailure
  else
    STDVA(ks) = 0.D0
    STDVB(ks) = 0.D0
  endif
endsubroutine phase_equilibrium_residuals

!--------------------------------------------------------------------------------------------------------------------------
!
!  Excess properties (IHEN = 1): excess volume (IFL = 1) or excess enthalpy
!  (IFL = 2).
!
subroutine excess_property_residuals(ks, FCION, numExp)
  integer, intent(in) :: ks
  double precision, intent(inout) :: FCION(*)
  integer, intent(inout) :: numExp

  double precision :: feed(NCMgca), Temp, Press
  integer :: i, n, nFailure, iErrEx, numExp1, numExp2

  nFailure = 0
  STDVA(ks) = 0D0
  write(*,'(/,"Set (excess property) No. ",I2)') ks
  do i = 1, NOBS(ks)
    n = numExp + 1
    Temp = TEXP(n)
    Press = PEXP(n)
    feed(1) = XEXP(n)
    feed(2) = 1D0 - feed(1)
    if (feed(2) < 0.D0) feed(2) = 1D-20
    write(*,'("  Point No. ",I2)') i
    ! GCXS does not set iErrEx (its error exits are commented out)
    iErrEx = 0
    if (negative_cohesion(Temp, 2, NG)) then
      iErrEx = 1
    else
      call GCXS(IFL(ks), NG, NST, Temp, Press, feed(:NCMgca), prop_E_calc(n), zLiCalc(n), iErrEx)
    endif
    if (iErrEx == 0) then
      if (prop_E_exp(n) /= 0D0) then
        FCAUX(n) = dabs(prop_E_calc(n) - prop_E_exp(n))
      else
        FCAUX(n) = 0
      endif
      if (IFL(ks) == 1) then
        FCION(n) = W11*FCAUX(n)*wEx(n)
      else
        FCION(n) = W12*FCAUX(n)*wEx(n)
      endif
    else
      nFailure = nFailure + 1
      zLiCalc(n) = 0.
      prop_E_calc(n) = 0.
      FCAUX(n) = 0.
      FCION(n) = PENALIZ*wEx(n)
      nFailed = nFailed + 1
    endif
    numExp = numExp + 1
    if (i == 1) numExp1 = numExp
    if (i == NOBS(ks)) numExp2 = numExp
  enddo
  if (NOBS(ks) - nFailure > 0) STDVA(ks) = sqrt(sum(FCION(numExp1:numExp2)**2)/(NOBS(ks) - nFailure))
endsubroutine excess_property_residuals

!--------------------------------------------------------------------------------------------------------------------------
!
!  Solubility of a solid (component IDC(ks,1)) in a liquid (IHEN = 2).
!
subroutine sle_residuals(ks, NCw, FCION, numExp)
  integer, intent(in) :: ks, NCw
  double precision, intent(inout) :: FCION(*)
  integer, intent(inout) :: numExp

  double precision :: feed(NCMgca), zFrac(NCMgca,2), Temp, Press, CZ
  double precision :: wwDHfu, wwDVfu, wwTfu, wwPfu, wwDCpfu, wwTtrSS, wwDHss
  integer :: i, n, i1, nFailure, iErrSol, numExp1, numExp2

  ! Melting data of the solid (copies: GCSLE could change its arguments)
  i1 = IDC(ks,1)
  wwDHfu = DHfu(i1)
  wwDVfu = DVfu(i1)
  wwTfu = Tfu(i1)
  wwPfu = Pfu(i1)
  wwDCpfu = DCpfu(i1)
  wwTtrSS = TtrSS(i1)
  wwDHss = DHss(i1)
  nFailure = 0
  STDVA(ks) = 0
  write(*,'(/,"Set (SLE) No. ",I2)') ks
  do i = 1, NOBS(ks)
    n = numExp + 1
    Temp = TEXP(n)
    Press = PEXP(n)
    feed(1) = XEXP(n)
    feed(2) = 1D0 - feed(1)
    iErrSol = 0
    ! zFrac(:,1) gets the liquid composition, zFrac(:,2) the fugacity (kept
    ! only because the original program printed it).
    if (negative_cohesion(Temp, NCw, NG)) then
      iErrSol = 1
    else
      call GCSLE(NCw, NG, NST, Temp, Press, feed(:NCMgca), wwDHfu, wwDVfu, wwTfu, wwPfu, wwDCpfu, &
                 wwTtrSS, wwDHss, CZ, zFrac(:NCMgca,1), zFrac(:NCMgca,2), iErrSol)
    endif
    if (iErrSol == 0) then
      zLiCalc(n) = CZ
      XCALC(n) = zFrac(1,1)
      YCALC(n) = zFrac(1,2)
      FCAUX(n) = dabs(zFrac(1,1) - feed(1))
      FCION(n) = wEx(n)*FCAUX(n)*W13
    else
      zLiCalc(n) = 0
      XCALC(n) = 0
      YCALC(n) = 0
      nFailure = nFailure + 1
      FCION(n) = PENALIZ*wEx(n)
      nFailed = nFailed + 1
      FCAUX(n) = 0
    endif
    numExp = numExp + 1
    if (i == 1) numExp1 = numExp
    if (i == NOBS(ks)) numExp2 = numExp
  enddo
  if (NOBS(ks) - nFailure > 0) STDVA(ks) = sum(FCAUX(numExp1:numExp2))/(NOBS(ks) - nFailure)
endsubroutine sle_residuals

!--------------------------------------------------------------------------------------------------------------------------
!
!  Critical points of binary mixtures (IHEN = 3): at the experimental
!  critical T, P and x, the first and second derivatives of ln f1 with
!  respect to z1 must vanish.
!
subroutine critical_point_residuals(ks, NCw, FCION, numExp)
  integer, intent(in) :: ks, NCw
  double precision, intent(inout) :: FCION(*)
  integer, intent(inout) :: numExp

  double precision :: feed(NCMgca), Temp, Press, dlnF_dz1, d2lnF_dz12
  integer :: i, n

  STDVA(ks) = 0D0
  write(*,'("Critical point data ", I2)') ks
  do i = 1, NOBS(ks)
    n = numExp + 1
    Temp = TEXP(n)
    Press = PEXP(n)
    feed(1) = XEXP(n)
    feed(2) = 1D0 - feed(1)
    if (negative_cohesion(Temp, NCw, NG)) then
      zLiCalc(n) = 0
      FCION(n) = PENALIZ*wBin(1,n)**2
      nFailed = nFailed + 1
      FCAUX(n) = 0
      FCAUX2(n) = 0
    else
      call implicitPCrit(NCw, NG, NST, Temp, Press, feed(:NCMgca), dlnF_dz1, d2lnF_dz12, zLiCalc(n))
      FCION(n) = wBin(1,n)**2*dlnF_dz1**2 + wBin(2,n)**2*d2lnF_dz12**2
      FCAUX(n) = dlnF_dz1
      FCAUX2(n) = d2lnF_dz12
    endif
    numExp = numExp + 1
  enddo
endsubroutine critical_point_residuals

!--------------------------------------------------------------------------------------------------------------------------
!
!  Activity coefficients in binary or ternary mixtures (IHEN = 5). The
!  residual uses absolute errors; the statistics (FCAUX), relative ones.
!
subroutine activity_coefficient_residuals(ks, NCw, FCION, numExp)
  integer, intent(in) :: ks, NCw
  double precision, intent(inout) :: FCION(*)
  integer, intent(inout) :: numExp

  double precision :: feed(NCMgca), lnGamma(NCMgca), dlnGammadT(NCMgca), Temp, Press
  integer :: i, n, io, nFailure, iTyp, iGz, jTemp

  nFailure = NOBS(ks)
  STDVA(ks) = 0D0
  STDVB(ks) = 0D0
  STDVC(ks) = 0D0
  write(6,'("Activity coefficient data ",I2)') ks
  do i = 1, NOBS(ks)
    n = numExp + 1
    Temp = TEXP(n)
    Press = PEXP(n)
    feed(1) = XEXP(n)
    if (mez(ks) == 2) then
      feed(2) = 1.D0 - feed(1)
    else
      feed(2) = xExp2(n)
      feed(3) = 1.D0 - feed(1) - feed(2)
    endif
    do io = 1, mez(ks)
      if (feed(io) <= 2.E-16) feed(io) = 2E-16
    enddo
    write(6,*) "OBS=", NOBS(ks)
    iTyp = 1
    iGz = 0
    jTemp = 0
    if (negative_cohesion(Temp, NCw, NG)) then
      ! Not calculated: g(T) <= 0
      nFailure = nFailure - 1
      rGamCalc1(n) = 0
      rGamCalc2(n) = 0
      FCAUX(n) = 0
      FCAUX2(n) = 0
      if (mez(ks) == 3) then
        rGamCalc3(n) = 0
        FCAUX3(n) = 0
      endif
      zLiCalc(n) = 0
      FCION(n) = PENALIZ
      nFailed = nFailed + 1
    else
      call paraGC(Temp, NCw, NG, NST, jTemp)
      call GAMMA(NCw, NG, NST, iTyp, iGz, jTemp, feed(:NCMgca), Temp, Press, lnGamma, dlnGammadT, zLiCalc(n))

      rGamCalc1(n) = dexp(lnGamma(1))
      rGamCalc2(n) = dexp(lnGamma(2))
      FCAUX(n) = rGamCalc1(n)/rGamExp1(n) - 1.D0
      FCAUX2(n) = rGamCalc2(n)/rGamExp2(n) - 1.D0
      FCION(n) = W8**2*(rGamCalc1(n) - rGamExp1(n))**2*wGam(1,n)**2 + &
                 W9**2*(rGamCalc2(n) - rGamExp2(n))**2*wGam(2,n)**2
      if (mez(ks) == 3) then
        rGamCalc3(n) = dexp(lnGamma(3))
        FCAUX3(n) = rGamCalc3(n)/rGamExp3(n) - 1.d0
        FCION(n) = FCION(n) + W10**2*(rGamCalc3(n) - rGamExp3(n))**2*wGam(3,n)**2
      endif
      FCION(n) = sqrt(FCION(n))
    endif

    STDVA(ks) = STDVA(ks) + FCAUX(n)**2
    STDVB(ks) = STDVB(ks) + FCAUX2(n)**2
    STDVC(ks) = STDVC(ks) + FCAUX3(n)**2
    numExp = numExp + 1
  enddo
  if (nFailure /= 0) then
    STDVA(ks) = dsqrt(STDVA(ks)/nFailure)
    STDVB(ks) = dsqrt(STDVB(ks)/nFailure)
    STDVC(ks) = dsqrt(STDVC(ks)/nFailure)
  else
    STDVA(ks) = 0.D0
    STDVB(ks) = 0.0D0
    STDVC(ks) = 0.
  endif
endsubroutine activity_coefficient_residuals

endmodule residuals
