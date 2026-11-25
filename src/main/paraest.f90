!--------------------------------------------------------------------------------------------------------------------------
!
!  paraest: parameter fitting for the GCA-EoS.
!
!  The molecular structure can be defined differently for the dispersive and
!  the association contributions.
!
!  The critical diameter dc of each component can be:
!    - estimated from the critical point (explicit equation, no fit);
!    - fitted to one point of the component's vapor pressure curve;
!    - fitted like any other parameter, i.e. to all the data (binary and/or
!      pure-component) in which the component takes part.
!
!  The input is read from paraest.dat; the fitted parameters are scaled
!  (module scaling) before being passed to the optimizer.
!
!--------------------------------------------------------------------------------------------------------------------------
program paraest
  use dimensions
  use fit_data
  use gca_parameters
  use scaling, only: XXSCAL, XXUNSCAL, search_box
  use config, only: read_options, optimizer, pso_particles, pso_iterations, pso_seed, pso_starts, pso_ssq_factor, &
                    pso_cluster_dist
  use pso_optimizer, only: PSO
  use residuals, only: FUNC, write_failed_points, nFailed, nNaN, nFailedFirst, nNaNFirst   ! residuals of the objective function, for ZXSSQ
  implicit none
!
  ! The data shared with FUNC (system, model parameters, experimental data,
  ! weights, identification of the fitted parameters) are in modules fit_data
  ! and gca_parameters. The meaning of IDXX is in module scaling.
!
  ! Locals:
  INTEGER NSIG, NPRED, NPAR, NCXP, IOPT, MAXFN, IXJAC, INFER, IER
  INTEGER NTEXT(NSetM,65)              ! title of each data set
  DOUBLE PRECISION XX(NPARM)           ! fitted parameters (scaled during the fit)
  DOUBLE PRECISION PARM(4)             ! ZXSSQ parameters
  DOUBLE PRECISION EPS, DELTAZ, SSQ, boundFrac
  DOUBLE PRECISION XJAC(NDM,NPARM), XJTJ(NPARM*(NPARM+1)/2), WORK(LWORK)
  DOUBLE PRECISION FCION(NDM)          ! residuals
  INTEGER iGr
  DOUBLE PRECISION epxOut               ! epx of spline_temperature (not used here)
  ! PSO: search box and every evaluation (model units), for the LM starts
  DOUBLE PRECISION boxLo(NPARM), boxHi(NPARM)
  DOUBLE PRECISION, ALLOCATABLE :: archX(:,:), archF(:)
  INTEGER, ALLOCATABLE :: archN(:)
  INTEGER nArch
!
  ! IDC(J,1): component whose experimental mole fractions are given in data set j.
  ! IDC(J,2): component whose mole fractions are NOT given.
!
  ! Read the input file paraest.dat
  CALL read_input(NSIG,NPRED,NPAR,PARM,NTEXT,XX,NCXP,IOPT)
!
  ! Spline temperatures of the initial parameters (for the output only: the
  ! model computes them in fakePARGR).
  do iGr = 1, NG
    call spline_temperature(TSTR(iGr), GDOT(iGr), GDDOT(iGr), TSPLINE(iGr), epxOut)
  enddo
!
  ! Write the input data and initial parameters to paraest.out
  CALL write_output_header(NSIG,NPRED,NPAR,PARM,NTEXT,XX,IOPT)
!
  ! Options of the optional file paraest.cfg (e.g. the residual of failed points)
  CALL read_options()
!
  ! Type of each value in XX, to assign it the right scaling factor
  CALL set_param_types(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,numTc,NUMDC)
!
  ! Scale the values of XX. In a fit, initial values on a bound (or very
  ! close) are moved 0.1 % of the interval away, so the optimizer can move
  ! them; in prediction they are used as given.
  IF(NPRED.EQ.1)THEN
    boundFrac = 0D0
  ELSE
    boundFrac = 1D-3
  ENDIF
  CALL XXSCAL(NPAR,XX,IDXX,boundFrac)
!
  ! Prediction or fit
  IF(NPRED.EQ.1)THEN
    ! Prediction
    CALL FUNC(XX,NCXP,NPAR,FCION)
    CALL write_failed_points(.FALSE.)
!
  ELSE
    ! Parameter fit. ZXSSQ parameters:
    EPS = 1.D-6
    DELTAZ = 1.D-6
    ! Maximum number of FUNC evaluations. ZXSSQ spends NPAR+1 evaluations per
    ! iteration on the Jacobian alone (finite differences): with a fixed
    ! MAXFN = 100, 8 parameters gave ~10 iterations and the fit ended with
    ! IER = 133. Now it allows ~50 iterations, whatever is fitted.
    MAXFN = 50*(NPAR + 1)
    PARM(3)=2.D2
    PARM(4)=.1D0
    IXJAC = NCXP
!
    ! Optimizer (option of paraest.cfg): PSO explores the search box; then
    ! Levenberg-Marquardt starts from the best point found (or from the
    ! initial values, with 'lm').
    IF(optimizer.NE.'lm') CALL run_pso()
    IF(optimizer.EQ.'lm') THEN
      CALL ZXSSQ(FUNC,NCXP,NPAR,NSIG,EPS,DELTAZ,MAXFN,IOPT,PARM,XX,SSQ,FCION,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
    ELSEIF(optimizer.EQ.'pso+lm') THEN
      CALL multi_start_lm()
    ENDIF
!
    ! One more evaluation at the final parameters: ZXSSQ's last one may be a
    ! rejected trial step or a Jacobian increment, and the calculated values
    ! written by write_results come from the last evaluation.
    CALL FUNC(XX,NCXP,NPAR,FCION)
!
    ! Statistics of the fit (PSO alone has no Jacobian: only SSQ)
    IF(optimizer.NE.'pso') CALL write_fit_stats(SSQ,XJAC,WORK,INFER,IER)
    CALL write_failed_points(.TRUE.)
!
    ! Spline temperatures of the final parameters (for the output)
    do iGr = 1, NG
      call spline_temperature(TSTR(iGr), GDOT(iGr), GDDOT(iGr), TSPLINE(iGr), epxOut)
    enddo
!
  END IF
!
  ! Final parameters, and experimental and calculated values
  CALL write_results(NTEXT)
!
contains

!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Particle swarm (PSO) over the search box of the parameters, in model
  ! units; XX ends with the best point found (scaled) and SSQ with its value.
!
subroutine run_pso()
  double precision :: x0(NPARM), lo(NPARM), hi(NPARM), xBest(NPARM)
  integer :: i, nEval, nLost

  CALL XXUNSCAL(NPAR, XX, IDXX, x0)
  CALL search_box(NPAR, IDXX, x0, boundFrac, lo, hi)
  boxLo = lo
  boxHi = hi
  allocate(archX(NPAR, pso_particles*(pso_iterations + 1)), archF(pso_particles*(pso_iterations + 1)), &
           archN(pso_particles*(pso_iterations + 1)))
  nArch = 0

  OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
  WRITE(1,'(/,X,A,/,X,A)') 'PSO search box:', '  Param. Type      Minimum      Initial      Maximum'
  DO i = 1, NPAR
    WRITE(1,'(X,I6,I6,3ES13.5)') i, IDXX(i), lo(i), x0(i), hi(i)
  ENDDO
  WRITE(1,'(/,X,A)') 'Course of the PSO (one line per evaluation):'
  CLOSE(UNIT=1)

  CALL PSO(pso_objective, NPAR, x0(:NPAR), lo(:NPAR), hi(:NPAR), pso_particles, pso_iterations, pso_seed, &
           xBest(:NPAR), SSQ, nLost, nEval)

  XX(:NPAR) = xBest(:NPAR)
  CALL XXSCAL(NPAR, XX, IDXX, boundFrac)

  OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
  WRITE(1,'(/,X,A,ES20.12,/,X,A,I0,/,X,A,I0)') 'Best SSQ found by PSO: ', SSQ, &
    'Points not calculated beyond those at the initial parameters: ', nLost, 'Number of PSO evaluations: ', nEval
  WRITE(1,'(X,A,/,(X,5ES16.8))') 'Best parameters found by PSO:', xBest(:NPAR)
  CLOSE(UNIT=1)
endsubroutine run_pso

!------------------------------------------------------------------------------------------------------------------------------------
!
!  Objective function for PSO at the parameters x (model units): SSQ, and the
!  number of points not calculated (failed or NaN) beyond those at the
!  initial parameters (the first evaluation, particle 1). Points that already
!  fail at the initial values do not count: the model may not reach them.
!
subroutine pso_objective(n, x, ssqX, nLostX)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: ssqX
  integer, intent(out) :: nLostX
  double precision :: xs(NPARM), f(NDM)

  xs(:n) = x
  CALL XXSCAL(n, xs, IDXX, boundFrac)
  CALL FUNC(xs, NCXP, n, f)
  ssqX = sum(f(:NCXP)**2)
  nLostX = max(0, nFailed + nNaN - nFailedFirst - nNaNFirst)
  if (nArch < size(archF)) then
    nArch = nArch + 1
    archX(:,nArch) = x
    archF(nArch) = ssqX
    archN(nArch) = nLostX
  endif
endsubroutine pso_objective

!------------------------------------------------------------------------------------------------------------------------------------
!
!  Distance between two points (model units): RMS of the differences divided
!  by the width of the PSO box of each parameter (fixed parameters, with
!  zero width, do not count).
!
double precision function box_distance(a, b)
  double precision, intent(in) :: a(:), b(:)
  integer :: j, m
  box_distance = 0D0
  m = 0
  do j = 1, NPAR
    if (boxHi(j) > boxLo(j)) then
      box_distance = box_distance + ((a(j) - b(j))/(boxHi(j) - boxLo(j)))**2
      m = m + 1
    endif
  enddo
  if (m > 0) box_distance = sqrt(box_distance/m)
endfunction box_distance

!------------------------------------------------------------------------------------------------------------------------------------
!
!  Indices of the n values (f, nb) in increasing order of nb, then f (the
!  order of the feasibility rule); a stable insertion sort.
!
subroutine rank(n, f, nb, idx)
  integer, intent(in) :: n, nb(n)
  double precision, intent(in) :: f(n)
  integer, intent(out) :: idx(n)
  integer :: i, j, t
  idx = [(i, i = 1, n)]
  do i = 2, n
    t = idx(i)
    j = i - 1
    do while (j >= 1)
      if (.not. (nb(t) < nb(idx(j)) .or. (nb(t) == nb(idx(j)) .and. f(t) < f(idx(j))))) exit
      idx(j+1) = idx(j)
      j = j - 1
    enddo
    idx(j+1) = t
  enddo
endsubroutine rank

!------------------------------------------------------------------------------------------------------------------------------------
!
!  'pso+lm': Levenberg-Marquardt from up to pso_starts PSO points in distinct
!  regions of the search box (see module config), and the best result
!  (fewest points not calculated, then the lowest SSQ). All the distinct
!  optima are reported in paraest.out. On exit XX, SSQ, XJAC, WORK, INFER and
!  IER are those of the best one (LM is run again from it if it was not the
!  last run, so that the statistics are its own).
!
subroutine multi_start_lm()
  double precision :: f(nArch), fMax, xs(NPARM)
  double precision, allocatable :: start(:,:), xEnd(:,:), xEndS(:,:), ssqEnd(:)
  integer, allocatable :: nEnd(:), inferEnd(:), ierEnd(:), same(:), order(:)
  integer :: idx(nArch), nb(nArch), i, k, m, nS, nMin, best

  ! Candidates, best first
  f = archF(:nArch)
  nb = archN(:nArch)
  where (f /= f) f = huge(1D0)
  call rank(nArch, f, nb, idx)
  nMin = nb(idx(1))
  fMax = pso_ssq_factor*f(idx(1))
  allocate(start(NPAR, pso_starts))
  nS = 0
  do i = 1, nArch
    k = idx(i)
    if (nb(k) > nMin + 1 .or. f(k) > fMax) cycle
    if (any([(box_distance(archX(:,k), start(:,m)) <= pso_cluster_dist, m = 1, nS)])) cycle
    nS = nS + 1
    start(:,nS) = archX(:,k)
    if (nS == pso_starts) exit
  enddo

  allocate(xEnd(NPAR,nS), xEndS(NPAR,nS), ssqEnd(nS), nEnd(nS), inferEnd(nS), ierEnd(nS), same(nS), order(nS))
  do i = 1, nS
    OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
    WRITE(1,'(/,X,A,I0,A,I0,A,/,(X,5ES16.8))') 'Levenberg-Marquardt from PSO point ', i, ' of ', nS, ':', start(:,i)
    CLOSE(UNIT=1)
    XX(:NPAR) = start(:,i)
    CALL XXSCAL(NPAR, XX, IDXX, boundFrac)
    CALL ZXSSQ(FUNC,NCXP,NPAR,NSIG,EPS,DELTAZ,MAXFN,IOPT,PARM,XX,SSQ,FCION,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
    ! Evaluation at the result: SSQ and points not calculated
    CALL FUNC(XX,NCXP,NPAR,FCION)
    ssqEnd(i) = sum(FCION(:NCXP)**2)
    nEnd(i) = nFailed + nNaN
    inferEnd(i) = INFER
    ierEnd(i) = IER
    CALL XXUNSCAL(NPAR, XX, IDXX, xs)
    xEnd(:,i) = xs(:NPAR)
    xEndS(:,i) = XX(:NPAR)            ! scaled: LM is run again from it exactly
  enddo

  ! Distinct optima: in order of quality, a result is the same optimum as a
  ! better one closer than pso_cluster_dist
  call rank(nS, ssqEnd, nEnd, order)
  same = 0
  do i = 2, nS
    do m = 1, i - 1
      if (same(order(m)) == 0 .and. box_distance(xEnd(:,order(i)), xEnd(:,order(m))) <= pso_cluster_dist) then
        same(order(i)) = order(m)
        exit
      endif
    enddo
  enddo
  best = order(1)

  OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
  WRITE(1,'(/,X,A,I0,A,I0,A)') 'Results of Levenberg-Marquardt from ', nS, ' PSO points (', &
    count(same(:nS) == 0), ' distinct optima), best first:'
  WRITE(1,'(X,A)') ' Start  Lost        SSQ       INFER   IER  Same as  Parameters'
  do m = 1, nS
    i = order(m)
    if (same(i) == 0) then
      WRITE(1,'(X,I6,I6,ES16.8,2I6,9X,*(ES15.7))') i, nEnd(i), ssqEnd(i), inferEnd(i), ierEnd(i), xEnd(:,i)
    else
      WRITE(1,'(X,I6,I6,ES16.8,2I6,I9,*(ES15.7))') i, nEnd(i), ssqEnd(i), inferEnd(i), ierEnd(i), same(i), xEnd(:,i)
    endif
  enddo
  WRITE(1,'(X,A,I0)') 'Final parameters: those of start ', best
  CLOSE(UNIT=1)

  ! XX and the statistics of the best result
  if (best /= nS) then
    OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
    WRITE(1,'(/,X,A,I0,A)') 'Levenberg-Marquardt again from the result of start ', best, ' (for its statistics):'
    CLOSE(UNIT=1)
    XX(:NPAR) = xEndS(:,best)
    CALL ZXSSQ(FUNC,NCXP,NPAR,NSIG,EPS,DELTAZ,MAXFN,IOPT,PARM,XX,SSQ,FCION,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
  endif
endsubroutine multi_start_lm

endprogram paraest
!
!-------------------------------------------------------------------------------------------------------------------------
!
  ! Reads the input file paraest.dat. It also echoes what it reads, with
  ! format, to the screen, to help find an error in the input.
!
SUBROUTINE read_input(NSIG,NPRED,NPAR,PARM,NTEXT,XX,NCXP,IOPT)
  use dimensions
  use fit_data, only: NC, W1, W2, W3, W4, W5, w6, W7, W8, W9, W10, W11, W12, w13, ISAT, IGLE, IVE, IUNI, &
&                      PRCRIT, TEMPCRIT, NSAT, IDC, NOBS, iHEN, mez, IFL, rGamExp1, rGamExp2, rGamExp3, &
&                      prop_E_calc, prop_E_exp, TEXP, PEXP, XEXP, xexp2, YEXP, yExp2, rhoExp, IDSAT, IOBJSAT, &
&                      wSat, wBin, wTer, wEx, wGam, DHfu, Tfu, DCpfu, Dvfu, TtrSS, DHss, Pfu
  use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&                      ENASS, RKASS, NNP, NNQ, NALFA, NENASS, NRKASS, TSAT, PSAT, dc, IDG, ITEMP, NIDG, IDIDG, &
&                      iDQ, IDALF, NIDALF, IDIDALF, IDENASS, NIDENASS, IDIDENASS, IDRKASS, NIDRKASS, &
&                      IDIDRKASS, NUMDC, IDDC, numTc, idTc, check_group, check_assoc
  use gca_bridge, only: load_from_pargr, load_from_parmol
  use scaling, only: dc_mean, error_Dc
  use utilities, only: i2s, rep, check_limit
  use model_output, only: keep_pargr_warnings
!
  implicit none
  integer :: I, IOPT, J, K, L, NCXP, NPAR, NPRED, NSIG, NSSA, NSSE, NSSRK
  integer :: NSSS, NST, NSUMK, NSUML
  double precision :: aIOPT
!
  double precision :: PARM(4)               ! ZXSSQ parameters
  integer :: NTEXT(NSetM,65)                ! title of each data set
  double precision :: XX(NPARM)             ! initial values of the fitted parameters
  double precision :: error_Dc1(NCM)
!
  OPEN(UNIT=1,STATUS='OLD',FILE='paraest.dat',FORM='FORMATTED',ACCESS='SEQUENTIAL')
  READ (1, *) W1
  read (1, *) W2, W3, W4, W5, W6, W7
  read (1, *) W8, W9, W10
  read (1, *) W11, W12
  read (1, *) w13
  ! IOPT is read as a real because data files usually have "2." (ifx accepts
  ! it for an integer, gfortran does not).
  read (1, *) PARM(1),PARM(2),aIOPT
  IOPT = nint(aIOPT)
  write(*,'(15F10.4,I2)')W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,w13,PARM(1),PARM(2),IOPT
!
1  FORMAT(8F10.4)
2  FORMAT(*(I3))
!
  READ(1,*) NC,NG,NGA,NSIG
  write(*,'(4I3)')NC,NG,NGA,NSIG
  READ(1,*) ISAT,IGLE,IVE,IUNI,NPRED
  write(*,'(5I3,/)')ISAT,IGLE,IVE,Iuni,NPRED
!
  ! Count the equilibrium data points
  NCXP = 0
  IF(ISAT.NE.0)then

    READ(1,*) (IDSAT(I), I=1,ISAT)
    write(*,'(*(I3,X))') (IDSAT(I), I=1,ISAT)
    READ(1,*) (NSAT(IDSAT(I)), I=1,ISAT)
    write(*,'(*(I3,X))')(NSAT(IDSAT(I)), I=1,ISAT)
    READ(1,*) (IOBJSAT(IDSAT(I)), I=1,ISAT)
    write(*,'('//rep(iSat,"I3,X")//'/)') (IOBJSAT(IDSAT(I)), I=1,ISAT)
!
    DO I=1,ISAT

      NSUMK = 0
      do K=1,(IDSAT(I)-1)

        NSUMK = NSUMK + 1

      end do
      NSUML = 0
      do L=1,NSUMK

        NSUML = NSUML + NSAT(L)

      end do
    !Pure-component equilibrium data: T and P, or T, P and rho.
      DO J= (NSUML +1 ),( NSUML + NSAT(IDSAT(I)) )

        if(iObjSat(idSat(i)).NE.81)then
          READ(1,*) TEXP(J),PEXP(J),wSat(j)
          write(*,'(F7.2,F12.8,I12)') TEXP(J),PEXP(J), j
        else
          READ(1,*) TEXP(J),PEXP(J),rhoExp(j),wSat(j)
          write(*,'(F7.2,2(F12.8,X))') TEXP(J),PEXP(J),rhoExp(j)
      endif

      END DO
!
      NCXP = NCXP + NSAT(IDSAT(I)) ! number of data points read so far
!
    END DO

  endif
  write(*,'(//)')
  ! Attractive and associative parameters: read by PARGR of the GCA package
  ! (v1.14), as in the flash programs: group composition of each component
  ! as pairs (database number, number of groups), parameters from the
  ! database FOR010-FOR015 plus those given by the user. What PARGR writes
  ! goes to a scratch file (unit 31), closed at once; its warnings are kept
  ! and write_output_header copies them to paraest.out.
  open(unit=31, status='scratch', form='formatted')
  call PARGR(NC, NG, NGA, NST, 1, 31)
  call keep_pargr_warnings(31)
  call load_from_pargr(NC)
!-----------------------------------------------------------------------
!
  write(*,'(//)')
  ! Pure-component properties of the GCA: read by PARMOL of the package
  ! (v1.14): unit option, then for each component its name, "Tc Pc omega
  ! Tsat Psat dc" and "M c A B C D E F G" (Peneloux, Passut-Danner). PARMOL
  ! estimates dc from the critical point when it is 0 and fits it to
  ! (Tsat, Psat) when Tsat > 0. What it writes is discarded (unit 32).
  open(unit=32, status='scratch', form='formatted')
  call PARMOL(NC, NG, NST, 1, 32)
  close(32)
  call load_from_parmol(NC, TEMPCRIT(:NC), PRCRIT(:NC), TSAT(:NC), PSAT(:NC), DC(:NC))
  DO I=1,NC
    write(*,'(2(F7.2,2X,F12.6,2X),F9.5,2X,"! Comp. No.",I2)')TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I),DC(I),i
  END DO
  write(*,'(//)')
  call read_fusion(NC)
  write(*,'(//)')
!
  ! (The group composition of each component was read by PARGR.)
!
  ! Data sets (binary or ternary equilibrium, excess properties, ...)
  DO I=1,IGLE
!
    read(1,*)IHEN(I)
    IF(IHEN(I).EQ.0) THEN

      write(*,'(//,I2,5X,"Binary equilibrium data")')iHen(i)
      READ(1,*)IFL(I)
      write(*,'(i2)')iFL(i)
      if(iFL(i).LT.3)then

        mez(i) = 2
        READ(1,*) NOBS(I),IDC(I,1),IDC(I,2)
        write(*,'(3(I2,4X))')NOBS(I),IDC(I,1),IDC(I,2)

      else

        mez(i) = 3
        READ(1,*) NOBS(I),(IDC(I,j),j=1,mez(i))
        write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,mez(i))

      endif

    elseif(iHen(i).EQ.1)then
!
      ! Excess property data
      read(1,*)ifl(i)
      if(ifl(i).EQ.1)then

        write(*,'(//,I2,5X,"Excess volume data")')ifl(i)

      elseif(ifl(i).EQ.2)then

        write(*,'(//,I2,5X,"Excess enthalpy data")')ifl(i)

      endif
      read(1,*)nObs(i), (idC(i,j), j = 1, 2)
      write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)

    elseif(iHen(i).EQ.2)then

      write(*,'(//,I2,5X,"Solid solubility data")')iHen(i)
      read(1,*)nObs(i), (idC(i,j), j = 1, 2)
      write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)

    elseif(iHen(i) == 3)then

      write(*,'(//, I2, 5X, "Critical point data")') iHen(i)
      read(1,*)nObs(i), (idC(i,j), j = 1, 2)
      write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)

    ELSEif(ihen(i).EQ.5)then

      write(*,'(//,I2,5X,"Activity coefficient data")')iHen(i)
      read(1,*) mez(i)
      write(*,'(I2)')mez(i)
      read(1,*) nObs(i), (idC(i,j),j=1,mez(i))
      write(*,'(*(I2,2X))')(idC(i,j),j=1,mez(i)),nObs(i)

    END IF
        READ(1,23) (NTEXT(I,J),J=1,65)
    write(*,23)(NTEXT(I,J),J=1,65)
23    FORMAT(40A2)
!
    write(*,'(/)')
    DO J=( NCXP+1 ),(NCXP + NOBS(I))

      If(ihen(i).EQ.0)then
!
        ! ****************** Phase equilibrium *****************
        if(mez(i).EQ.2)then

          ! Binary
          READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J),wBin(:2,j)
          write(*,'(F7.2,2X,F9.4,2(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),YEXP(J),j

        elseif(mez(i).EQ.3)then

          ! Ternary
          READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J),xExp2(j),yExp2(j),wBin(:2,j),wTer(:2,j)
          write(*,'(F7.2,2X,F9.4,4(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),YEXP(J),xExp2(j),yExp2(j),j

        endif
!
      elseif(iHen(i).EQ.1)then
!
        ! ****************** Excess properties *****************
        READ(1,*) TEXP(J),PEXP(J),XEXP(J),prop_E_exp(J),wEx(j)
        write(*,'(F7.2,2X,F9.4,2X,F7.4,F12.6, "! ", I3)')TEXP(J),PEXP(J),XEXP(J),prop_E_exp(J), j

      elseif(iHen(i).EQ.2)then

        ! ****************** SLE *****************
        read(1,*)Texp(j), Pexp(j), xExp(j), wEx(j)
        write(*,'(F7.2,2X,F9.4,2X,F7.4, "! ", I3)')TEXP(J),PEXP(J),XEXP(J), j

      elseif(iHen(i) == 3)then

        ! ****************** Critical point *****************
        read(1,*)Texp(j), Pexp(j), xExp(j), wBin(:2,j)
        write(*,'(F7.2,2X,F9.4,2X,F7.4, "! ", i3)')TEXP(J), PEXP(J), XEXP(J), j

      elseif(iHen(i).EQ.5)then
!
        ! ****************** Activity coefficients *****************
        if(mez(i).EQ.2)then

          read(1,*) TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),wGam(:2,j)
          write(*,'(F7.2,2X,F9.4,3(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),j

        elseif(mez(i).EQ.3)then

          read(1,*) TEXP(J),PEXP(J),XEXP(J),xexp2(J),rgamexp1(J),rgamexp2(J),rgamexp3(J),wGam(:3,j)
          write(*,'(F7.2,2X,F9.4,4(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),rgamexp3(J),j

        endif
!
      end if

    END DO
    NCXP = NCXP + NOBS(I)
!
    ! End of the loop over data sets
  END DO
!
  write(*,'(//)')
  READ(1,*) NNP,NNQ,NALFA,NENASS,NRKASS,numTc, NUMDC
  write(*,'(7(I2,4X),"! Interaction - Area - Non-random. - Assoc. energy - Assoc. volume - Crit. temp. - Crit. diameters")') NNP,NNQ,NALFA,NENASS,NRKASS,numTc, NUMDC
  call check_limit(NNP,    NPTIPOM, 'g/k interaction parameters (NNP)')
  call check_limit(NNQ,    NPTIPOM, 'area parameters (NNQ)')
  call check_limit(NALFA,  NPTIPOM, 'non-randomness parameters (NALFA)')
  call check_limit(NENASS, NPTIPOM, 'association energies (NENASS)')
  call check_limit(NRKASS, NPTIPOM, 'association volumes (NRKASS)')
  call check_limit(numTc,  NPTIPOM, 'critical temperatures (numTc)')
  call check_limit(NUMDC,  NPTIPOM, 'critical diameters (NUMDC)')
!
  if(NNP.NE.0) then
!
    ! Dispersive interaction parameters to fit: g(i,i) and k(i,j)
    READ(1,*) (IDG(I,1),IDG(I,2),I=1,NNP)
    write(*,'(*(2(I2,X),4X))')(IDG(I,1),IDG(I,2),I=1,NNP)
!
    READ(1,*) (ITEMP(I), I=1,NNP)
    write(*,'(*(I2,6X))')(ITEMP(I), I=1,NNP)
!
    READ(1,*) (NIDG(I), I=1,NNP)
    call check_limit(sum(NIDG(:NNP)), NIDM, 'identical g/k parameters')
    write(*,'(*(I2,6X))')(NIDG(I), I=1,NNP)
!
    NSSS = 0
    ! Identical parameters:
    DO L=1,NNP
!
      if(NIDG(L).NE.0) then
!
        READ(1,*) (IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
        write(*,'(*(I2,X))')(IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
        NSSS = NSSS + NIDG(L)
!
      end if
!
    END DO

  end if
!
  ! Area parameters to fit
  if(NNQ.NE.0)then

    READ(1,*) (IDQ(I), I=1,NNQ)
    write(*,'(*(I2,X))') (IDQ(I), I=1,NNQ)

  endif
!
  ! Non-randomness parameters to fit
  if(NALFA.NE.0) then
!
    READ(1,*) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
    write(*,'(*(2(I2,X),4X))')(IDALF(I,1),IDALF(I,2),I=1,NALFA)
!
    READ(1,*) (NIDALF(I),I=1,NALFA)
    call check_limit(sum(NIDALF(:NALFA)), NIDM, 'identical alpha parameters')
    write(*,'(*(I2,6X))')(NIDALF(I),I=1,NALFA)
    NSSA = 0
    DO L=1,NALFA

      if(NIDALF(L).NE.0) then

        READ(1,*) (IDIDALF(J,1),IDIDALF(J,2),J=(1+NSSA),(NIDALF(L)+NSSA))
        write(*,'(*(2I2,X))')
        NSSA = NSSA + NIDALF(L)

      end if

    END DO

  end if
!
  ! Association energies to fit
  if(NENASS.NE.0) then

    DO I=1,NENASS

      READ(1,*) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
      write(*,'(4I3)')IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)

    END DO
    READ(1,*) (NIDENASS(I), I=1,NENASS )
    call check_limit(sum(NIDENASS(:NENASS)), NIDM, 'identical association energies')
    write(*,'(*(I2,X))')(NIDENASS(I), I=1,NENASS )
!
    ! Identical ones:
    NSSE = 0
    DO L=1,NENASS

      if(NIDENASS(L).NE.0) then

        DO J=(1+NSSE),(NIDENASS(L)+NSSE)

          READ(1,*) IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)
          write(*,'(4I3)') IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)

        END DO
        NSSE = NSSE + NIDENASS(L)

      end if

    END DO
!
  end if
!
  ! Association volumes to fit:
  if(NRKASS.NE.0) then
!
    DO I=1,NRKASS

      READ(1,*) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
      write(*,'(4I3)') IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)

    END DO
    READ(1,*) (NIDRKASS(I), I=1,NRKASS )
    call check_limit(sum(NIDRKASS(:NRKASS)), NIDM, 'identical association volumes')
    write(*,'(*(I2,X))')(NIDRKASS(I), I=1,NRKASS )
!
    NSSRK = 0
    DO L=1,NRKASS

      if(NIDRKASS(L).NE.0) then

        DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

          READ(1,*) IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)
          write(*,'(4I3)')IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)

        END DO
        NSSRK = NSSRK + NIDRKASS(L)

      end if

    END DO

  end if
!
  ! Components whose Tc is fitted:
  IF(numTc.NE.0)then

    READ(1,*) (idTc(I), I=1,numtc)
    write(*,'('//rep(numTc,"I2,X")//'/)')(idTc(I), I=1,numtc)

  endif
!
  ! Components whose dc is fitted (IDDC), and the range allowed around the
  ! initial value of each one (error_Dc1).
  IF(NUMDC.NE.0)then

    READ(1,*) (IDDC(I), I=1,numDC)
    READ(1,*) (error_Dc1(i), I=1,numDC)
    write(*,*) (IDDC(I), I=1,numDC)
    write(*,*) (error_Dc1(i), I=1,numDC)

  endif
!
  ! The fitted parameters refer to the groups by their position in the
  ! system (column Pos of the parameter tables: groups in increasing order of
  ! their database number), 1..NG for attractive groups and 1..NGA for
  ! associating ones.
  do I = 1, NNP
    call check_group(IDG(I,:), 'g/k')
  enddo
  do I = 1, sum(NIDG(:NNP))
    call check_group(IDIDG(I,:), 'identical g/k')
  enddo
  if (NNQ > 0) call check_group(IDQ(:NNQ), 'area')
  do I = 1, NALFA
    call check_group(IDALF(I,:), 'alpha')
  enddo
  do I = 1, sum(NIDALF(:NALFA))
    call check_group(IDIDALF(I,:), 'identical alpha')
  enddo
  ! Association: (site, group, site, group)
  do I = 1, NENASS
    call check_assoc(IDENASS(I,:), 'association energy')
  enddo
  do I = 1, sum(NIDENASS(:NENASS))
    call check_assoc(IDIDENASS(I,:), 'identical association energy')
  enddo
  do I = 1, NRKASS
    call check_assoc(IDRKASS(I,:), 'association volume')
  enddo
  do I = 1, sum(NIDRKASS(:NRKASS))
    call check_assoc(IDIDRKASS(I,:), 'identical association volume')
  enddo
!
  ! Number of parameters
  NPAR = NNP+NNQ+NALFA+NENASS+NRKASS+numTc+NUMDC
  IF(NNP.NE.0) THEN

    DO I = 1,NNP

      NPAR = NPAR + ITEMP(I)

    END DO

  END IF
!
  ! Initial values of the parameters
  if(nPar > 0)READ(1,*)  (XX(I),I=1,NPAR)
!
  ! The dc values are also kept apart, for their scaling between bounds:
  if(numDc > 0)then

    dc_mean(nPar - numDc + 1:nPar) = xx(nPar - numDc + 1:nPar)
    error_Dc(nPar - numDc + 1:nPar) = error_Dc1(1:numDc)

  endif
!
  CLOSE(UNIT=1)
!
ENDsubroutine read_input
!
!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Melting and solid-solid transition properties of each component (any
  ! model), as GCSLE 2.0.0 (option 0) takes them:
  !    Tfus/K  Pfus/atm  DHfus/(J/mol)  DVfus/(cm3/mol)  ACpfus/(J/mol K)  BCpfus/(J/mol K2)  NTSS
  ! followed by NTSS lines "Ttr/K  DHtr/(J/mol)". DCpfus = ACpfus + BCpfus*T.
  ! For now paraest uses DCpfus = ACpfus and only the first transition.
!
subroutine read_fusion(NC)
  use fit_data, only: Tfu, Pfu, DHfu, DVfu, DCpfu, TtrSS, DHss
  implicit none
  integer, intent(in) :: NC
  integer :: i, k, NTSS
  double precision :: BCpfu, Ttr, DHtr

  do i = 1, NC
    read(1,*) Tfu(i), Pfu(i), DHfu(i), DVfu(i), DCpfu(i), BCpfu, NTSS
    write(*,'(F8.2,2X,G12.5,2X,G12.5,2X,G12.5,2X,G12.5,2X,G12.5,I3,2X,"! fusion, comp.",I3)') &
      Tfu(i), Pfu(i), DHfu(i), DVfu(i), DCpfu(i), BCpfu, NTSS, i
    if (BCpfu /= 0) write(*,'(" WARNING: BCpfus of component",I3," is not used yet (DCpfus = ACpfus)")') i
    TtrSS(i) = 0D0
    DHss(i) = 0D0
    do k = 1, NTSS
      read(1,*) Ttr, DHtr
      if (k == 1) then
        TtrSS(i) = Ttr
        DHss(i) = DHtr
      else
        write(*,'(" WARNING: solid-solid transition",I2," of component",I3," is not used yet")') k, i
      endif
    enddo
  enddo
endsubroutine read_fusion
!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Writes to paraest.out the data and parameters read by read_input (not
  ! the results).
!
SUBROUTINE write_output_header(NSIG,NPRED,NPAR,PARM,NTEXT,XX,IOPT)
  use dimensions
  use fit_data, only: NC, W1, W2, W3, W4, W5, w6, W7, W8, W9, W10, W11, W12, w13, ISAT, IGLE, IVE, IUNI, &
&                      PRCRIT, TEMPCRIT, NSAT, IDC, NOBS, IHEN, mez, IFL, TEXP, PEXP, XEXP, xexp2, YEXP, IDSAT
  use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&                      ENASS, RKASS, NNP, NNQ, NALFA, NENASS, NRKASS, TSAT, PSAT, IDG, ITEMP, NIDG, IDIDG, &
&                      IDQ, IDALF, NIDALF, IDIDALF, IDENASS, NIDENASS, IDIDENASS, IDRKASS, NIDRKASS, &
&                      IDIDRKASS, DC, NUMDC, IDDC, TSPLINE, numtc, idtc
  use fit_data, only: DHfu, Tfu, DCpfu, DVfu, TtrSS, DHss, Pfu
  use model_output, only: write_model_params, write_pargr_warnings
  use utilities, only: i2s, rep
!
  implicit none
  integer :: I, IOPT, J, k, L, NPAR, NPRED, NSIG, NSSA, NSSE, NSSRK, NSSS
!
  double precision :: PARM(4)               ! ZXSSQ parameters
  integer :: NTEXT(NSetM,65)                ! title of each data set
  double precision :: XX(NPARM)             ! initial values of the fitted parameters
!
  OPEN(UNIT=1,STATUS='unknown',file='paraest.out',FORM='FORMATTED')
!
  WRITE(1,1) W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,w13,PARM(1),PARM(2),IOPT
!
1  FORMAT(/,X,' Weights in the objective function:',/,                                              &
&         X,'Weight of vapor pressures:                       ',F10.4,/,                            &
&         X,'Weight of bubble pressure of binary mixtures:    ',F10.4,/,                            &
&         X,'Weight of bubble temperature of binary mixtures: ',F10.4,/,                            &
&         X,'Weight of composition of comp. 1 in phase II:    ',F10.4,/,                            &
&         X,'Weight of composition of comp. 1 in phase I:     ',F10.4,/,                            &
&         X,'Weight of composition of comp. 2 in phase I:     ',F10.4,' (ternary mixtures only)',/  &
&         X,'Weight of composition of comp. 2 in phase II:    ',F10.4,' (ternary mixtures only)',/  &
&         X,'Weight of activity coefficient of comp. 1:       ',F10.4,/                             &
&         X,'Weight of activity coefficient of comp. 2:       ',F10.4,/                             &
&         X,'Weight of activity coefficient of comp. 3:       ',F10.4,' (ternary mixtures only)',/  &
&         x,'Weight of excess volume:                         ',F10.4,/                             &
&         x,'Weight of excess enthalpy:                       ',F10.4,/                             &
&         x,'Weight of solid solubility:                      ',F10.4,/                             &
&         X,' PARM(1) = ',F10.4,'     PARM(2) = ',F10.4,/,                                          &
&         X,' Optimization algorithm: IOPT = ',I3)
!
  WRITE(1,'(/,X," Ratio between g*(i,i) and g(i,i) at which the ",/,&
&               X," spline is made = 4 ")')
!
  WRITE(1,'(/,X,"--------------------------------------------------",/,&
&            X,"--------------- INITIAL PARAMETERS ---------------",/,  &
&            X,"--------------------------------------------------")')
  call write_pargr_warnings(1)
  call write_model_params(1, NC)
!
  WRITE(1,18)
18  FORMAT(/,X,'Pure-component data')
  WRITE(1,19)
19  FORMAT (/, 4X, 'Tc/K      Pc/atm    Tsat/K    Psat/atm', 4X, 'dc/(cm/mol^(1/3))', X, "Tfus/K", X, "DHfus/(J/mol)", X, &
&              "DCpfus/(J/molK)", X, "Ttr/K", X, "DHtr/(J/mol)", X, "DVfus/(cm3/mol)", X, "Pfus/atm"/,X,' (initial dc estimated from the exptl. critical point, or given) ',/)
!
  DO I=1,NC

    WRITE(1,500) TEMPCRIT(I), PRCRIT(I), TSAT(I), PSAT(I), DC(I), Tfu(i), DHfu(i), DCpfu(i), TtrSS(i), DHss(i), DVfu(i), Pfu(i)
500    FORMAT(X, F7.2, 3X, F10.4, 2X,F7.2, X, G12.4, 6X, F7.4, 7X, F7.2, X, G12.4, 5X, F6.1, 5X, F7.2, 3(X, G12.4) )

  END DO
!
  WRITE(1,201)
201  FORMAT(/,X,' Parameters to estimate (groups by their position, column Pos of the tables)',/)
!
  if(NNP.NE.0) then
!
    WRITE(1,202)
202    FORMAT(' Parameters gii or kij ')
!
    WRITE(1,203) (IDG(I,1),IDG(I,2),I=1,NNP)
203    FORMAT(X,*(I5))
!
    WRITE(1,204) (ITEMP(I), I=1,NNP)
204    FORMAT(X,' Temperature dependence ',/,X,*(I3))
!
    WRITE(1,206) (NIDG(I), I=1,NNP)
206    FORMAT(X,' Number of identical par. ',/,X,*(I3))
!
    WRITE(1,208)
208    FORMAT(X,' Identical parameters ')
!
    NSSS = 0
!
    DO L=1,NNP

      if(NIDG(L).ne.0) then
!
        WRITE(1,210) (IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
210        FORMAT(X,*(I5))
        NSSS = NSSS + NIDG(L)
!
      end if

    END DO

  end if
!
  if(NNQ.NE.0) WRITE(1,211) (IDQ(I), I=1,NNQ)
211  FORMAT(/,X,' Area parameters  ',/,X,*(I5))
!
  if(NALFA.ne.0) then

    WRITE(1,213) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
213    FORMAT(/,X,' Non-randomness parameters ',/,X,*(I5))
!
    WRITE(1,215) (NIDALF(I),I=1,NALFA)
215    FORMAT(X,' Number of identical par. ',/,X,*(I3))
!
    NSSA = 0
    DO L=1,NALFA

          if(NIDALF(L).ne.0) then

        WRITE(1,217) (IDIDALF(J,1),IDIDALF(J,2),J=(1+NSSA),(NIDALF(L)+NSSA))
217        FORMAT(X,' Identical : ',*(I5))
        NSSA = NSSA + NIDALF(L)

      end if

    END DO

  end if
!
  if(NENASS.ne.0) then

  WRITE(1,219)
219  FORMAT(/,X,' Association energy parameters ')
!
    DO I=1,NENASS

      WRITE(1,221) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
221      FORMAT(X,4I3)

    END DO
    WRITE(1,223) (NIDENASS(I), I=1,NENASS )
223    FORMAT(X,' Number of identical par. ',*(I3))
    NSSE = 0
    DO L=1,NENASS

      if(NIDENASS(L).ne.0) then

        DO J=(1+NSSE),(NIDENASS(L)+NSSE)

          WRITE(1,225) IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)
225          FORMAT(' Identical: ',4I5)

        END DO
        NSSE = NSSE + NIDENASS(L)

      end if

    END DO

  end if
!
  if(NRKASS.ne.0) then

    wRITE(1,503)
503    FORMAT(/,X,' Association volume parameters ')

    DO I=1,NRKASS

      WRITE(1,229) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
229      FORMAT(X,4I3)

    END DO
    WRITE(1,231) (NIDRKASS(I), I=1,NRKASS )
231    FORMAT(X,' Number of identical parameters ',*(I3))
    NSSRK = 0
    DO L=1,NRKASS

      if(NIDRKASS(L).ne.0) then

        DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

          WRITE(1,232) IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)
232          FORMAT(X,' Identical: ',4I3)

        END DO
        NSSRK = NSSRK + NIDRKASS(L)

      end if

    END DO

  end if
  IF(numtc.ne.0) WRITE(1,233) (idtc(I), I=1,numTc)
233  FORMAT(/,X,' Critical temperatures ',/,X,*(I3))
  IF(NUMDC.NE.0) WRITE(1,234) (IDDC(I), I=1,NUMDC)
234   FORMAT(/,X,' Hard-sphere diameters ',/,X,*(I3))
!
  WRITE(1,235)  (XX(I),I=1,NPAR)
 235  FORMAT(/,X,' Initial values of the fitted parameters ',/,X,*(F13.4))
!
  WRITE(1,91)
91  FORMAT(//,X,' Course of the optimization ',X)
  if (NPAR > 0) then
    WRITE(1,"(//,X,' Parameter No. (SCALED). ',/,X,"//rep(NPAR,'I15')//"X,'        OBJ FUN ',/)") (I,I=1,NPAR)
  else   ! in prediction there are no parameters: only the title
    WRITE(1,"(//,X,' Parameter No. (SCALED). ',/)")
  endif
!
  CLOSE(UNIT=1)
!

ENDsubroutine write_output_header
!
!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Sets the type of each fitted parameter (IDXX, see module scaling), in
  ! the order they have in XX.
!
SUBROUTINE set_param_types(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,numTc,NUMDC)
  use dimensions
  use scaling, only: inIJ
!
  implicit none
  integer :: I, NALFA, NENASS, NNP, NNQ, NPAR, NRKASS, NUMDC, numTc
!
  ! Indices of the gii or kij parameters to fit.
  integer :: IDG(NPTIPOM,2), ITEMP(NPTIPOM), IDXX(NPARM)
!
  inIJ(:) = .FALSE.
  NPAR = 0
  IF(NNP.NE.0) THEN

    DO I=1,NNP

      IF(IDG(I,1).EQ.IDG(I,2)) THEN

        IDXX(NPAR+1)=1
        NPAR = NPAR + 1

      ELSE

        IDXX(NPAR+1)=2
        NPAR = NPAR + 1
        inIJ(nPar) = .TRUE.

      END IF
      IF(ITEMP(I).EQ.1)then

        IDXX(NPAR+1) = 3
        NPAR = NPAR + 1
        IF(IDG(I,1).NE.IDG(I,2))inIJ(nPar) = .TRUE.
!
      elseIF(ITEMP(I).EQ.2)then

        IDXX(NPAR+1) = 3
!
        ! g' and g" have different scaling types:
        idXX(nPar+2) = 10
        NPAR = NPAR + 2

      endif
!
    END DO

  END IF
  IF(NNQ.NE.0) THEN

    DO I=1, NNQ

      IDXX(NPAR+1)= 4
      NPAR = NPAR + 1

    eND DO

  END IF
  IF(NALFA.NE.0) THEN

    DO I=1, NALFA

      IDXX(NPAR+1)= 5
      NPAR = NPAR + 1

    END DO

  END IF
  IF(NENASS.NE.0) THEN

    DO I=1,NENASS

      IDXX(NPAR+1)= 6
      NPAR = NPAR + 1

    END DO

  END IF
  IF(NRKASS.NE.0) THEN

    DO I=1,NRKASS

      IDXX(NPAR+1)= 7
      NPAR = NPAR + 1

    END DO

  END IF
  IF(numTc.NE.0) THEN

    DO I=1,numtc

      IDXX(NPAR+1)= 8
      NPAR = NPAR + 1

    END DO

  END IF
  IF(NUMDC.NE.0) THEN

    DO I=1,NUMDC

      IDXX(NPAR+1)= 9
      NPAR = NPAR + 1

    END DO

  END IF
  return
ENDsubroutine set_param_types
!
!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Writes the statistics of the fit to paraest.out.
!
SUBROUTINE write_fit_stats(SSQ,XJAC,WORK,INFER,IER)
  use dimensions
!
  implicit none
  integer :: IER, INFER
  double precision :: SSQ
  double precision :: XJAC(NDM,NPARM), WORK(LWORK)
!
  OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
!
  WRITE(1,95)
95  FORMAT(/,X,' Final value of the objective function = SSQ = ',/)
  WRITE(1,*) SSQ
!
  WRITE(1,97) IER,INFER
 97  FORMAT(/,' (See the manual of routine ZXSSQ (IMSL) ',/,X,' IER = ',I5,/,X,' INFER = ',I5)
!
  WRITE(1,99) WORK(1),WORK(2),WORK(3),WORK(5)
99  FORMAT(/,'  Gradient norm = ',D12.5,/,X,' Number of objective function evaluations = ',F10.2,/, &
&         X,' Number of significant digits in the scaled parameters = ',F10.2,/,X,' Number of iterations = ',F10.2)
  CLOSE(UNIT=1)
  return

END SUBROUTINE write_fit_stats
!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Writes the final parameters, and the experimental data against the
  ! GCA-EoS results (with the given parameters or with the fitted ones), to
  ! paraest.out, tables.csv and graficos.csv.
!
SUBROUTINE write_results(NTEXT)
  use dimensions
  use fit_data, only: NC, ISAT, IGLE, IVE, IUNI, NSAT, IDC, NOBS, IHEN, mez, IFL, TEXP, PEXP, XEXP, xexp2, &
&                      YEXP, yExp2, IDSAT, tempcrit, XCALC, YCALC, xCalc2, yCalc2, ZLCALC, ZVCalc, FCAUX, &
&                      PCALC, IOBJSAT, FCAUX2, fcaux3, STDVA, STDVB, stdvc, rgamcalc1, rgamcalc2, rgamcalc3, &
&                      rgamexp1, rgamexp2, rgamexp3, zLicalc, prop_E_calc, prop_E_exp, DHVL
  use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&                      ENASS, RKASS, TSAT, PSAT, DC, numTc, idTc, NUMDC, IDDC, TSPLINE
  use model_output, only: write_model_params
  use utilities, only: i2s, rep
  use statistics, only: AAD, ARD, ASD, RSD
!
  implicit none
  integer :: I, iDim, iDim1, iDim2, iDim3, iDim4, J, JAMO, k, KS, L, mm
  integer :: N1, NCXP, NCXPy, NFOLI, No
  double precision :: E1, E2, E3, E4, E5, E6, E7, E8, er2, er4, Pmax, Pmin
!
  character(1) :: tab = achar(9)
  integer :: NTEXT(NSetM,65)                ! title of each data set
  double precision :: STDVSAT(NCM)          ! standard deviation of the pure components
!
  OPEN (UNIT = 1, STATUS='OLD', file = 'paraest.out', FORM = 'FORMATTED', POSITION='APPEND')
  open (unit = 2, file = "tables.csv", form = "formatted")
  open (unit = 10, file = "graficos.csv", form = "FORMATTED")
!
!
  WRITE(1,47)
 47  FORMAT(/,X,'--------------------------------------------------',/, &
&           X,' ---------------- FINAL PARAMETERS ----------------',/,   &
&           X,'--------------------------------------------------')
  call write_model_params(1, NC)
!
!
  ! Result of the dc fit
  DO I=1,NC

    DO JAMO=1,NUMDC

      IF(IDDC(JAMO).EQ.I)  THEN

        WRITE(1,1548) I,DC(I)
1548        FORMAT(/,X,' The hard-sphere diameter at the critical point of component No.',I3,/,&
&               ' was a fitted parameter ',/,X,' DC = ', F10.7)
        GO TO 505

      END IF

    END DO
!
    IF(TSAT(I).NE.0) WRITE(1,500) I,TSAT(I),PSAT(I),DC(I)
500    FORMAT(/,X,' The hard-sphere diameter at the critical point of component No.',I3,' was fitted to ',/, &
&           X,' reproduce the point of the vapor pressure curve ',/,x,' [ ',F10.4,' K ;',G10.4,' atm ] ',/,X,' DC = ', F10.7)
505  END DO
   DO I=1,NC
!
    DO mm=1,numtc

      IF(idtc(mm).EQ.I)  THEN

        WRITE(1,1549) I,tempcrit(I)
1549        FORMAT(/,X,'The critical temperature of component No.',I3,/,X,'was a fitted parameter ',/,'Tc = ', F8.2," K")

      END IF

    END DO

  END DO
!
!
  ! **************** Pure-component results *****************
!
  ! Besides paraest.out, the results go to graficos.csv (unit 10), only to make
  ! plots easily, and a summary to tables.csv (unit 2).
  NCXP = 0
 2  FORMAT(*(I3))
!
  IF(ISAT.NE.0) THEN

    WRITE (1, 338)
    write (2, 338)
    write (10, 338)
338    FORMAT(/,X,'***** Results: pure components *****',/)
    write (2, '("Compound,Tc (K),dc (cm/mol^1/3), DTr, ARD(P)%,Reference")')
    DO I=1,NC

      IF(NSAT(I).NE.0) THEN

        IF(IOBJSAT(I).EQ.50) THEN

          WRITE(1,393)
          write(10,393)
393          FORMAT(/,X,'  Comp   Texp(K)    Pexp(atm)       Pexp(atm)        Zliq       Zvap  REL,ERR,%,ISOF',/)

        END IF
        IF(IOBJSAT(I).EQ.80) THEN

          WRITE(1,333)
          write(10,3331)
333          FORMAT(/,X,'  Comp   Texp(K)    Pexp(atm)    Pcalc(atm)     Zliq      Zvap        DP%     DH_LV(J/mol) ',/)
3331          FORMAT(/,X,'Comp,Texp(K),Pexp(atm),Pcalc(atm),Zliq,Zvap,DP%,DH_LV(J/mol)',/)

        END IF
!
        NFOLI=NSAT(I)
        STDVSAT(I)=0D0
        DO J=1,NSAT(I)

          WRITE(1,334) I,TEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100, DHVL(NCXP+1)*.101325D0
          WRITE(10,3341) I,TEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100, DHVL(NCXP+1)*.101325D0
334          FORMAT(X,I3,4X,F10.4,X,G12.5,X,G12.5,X,G11.4,X,G11.4,X,G10.3,G11.4)
3341          FORMAT(I3,",",F10.4,",",G12.5,",",G12.5,",",G11.4,",",G11.4,",",G10.3,",",G11.4)
          NCXP = NCXP + 1

        END DO
!
        ! Average absolute and relative errors and standard deviations
        N1 = NCXP-nSat(i)+1
        E1 = AAD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i)); E2 = 100*ARD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i))
        E3 = ASD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i)); E4 = 100*RSD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i))
        write(*,335) E1, E2
        write(1,335) E1, E2, E3, E4
        write (2, '(I2, ",", F6.2, ",", F7.4, ",", F4.2, "-", F4.2, ",", F7.3  )') i, TempCrit(i), dc(i), Texp(N1)/TempCrit(i),Texp(NCXP)/TempCrit(i), E2
        write(10,335) E1, E2, E3, E4
335        format(/" AAD(ARD%) = ,",F10.5,", atm , (, ",F7.2,"% ,)", /, " AbsStdDev(Rel.%) = ,",F10.5,", atm , (, ",F7.2,"% ,)"/)

      END IF

    END DO

  END IF
!
  ! Mixture results
  IF(IGLE.NE.0) THEN
!
    WRITE(1,340)
    write (2, 340)
    write(10,340)
340    FORMAT(/,X,'********** RESULTS: BINARY MIXTURES *******************',/)
    write (2, '("Compound,-,T (K),P (kPa),AARD%,-,No exp.,Reference")')
    write (2, '("1,2,-,-,P,y1,points,")')
    DO KS=1,IGLE

      WRITE(1,339) (NTEXT(KS,J),J=1,65)
      WRITE(10,339) (NTEXT(KS,J),J=1,65)
339      FORMAT(//,80A2)
!
      ! Width of the number in parentheses when printing x(i): with blanks, Excel
      ! does not split the column headers correctly.
      if(idC(ks,1).LT.10)then

        iDim = 1

      else

        iDim = 2

      endif
      if(idC(ks,2).LT.10)then

        iDim2 = 1

      else

        iDim2 = 2

      endif
      IF(IHEN(KS).EQ.0) THEN

        ! ******************** Binary equilibrium ***********************
!
        IF(IFL(KS).EQ.-2) WRITE(1,343) IDC(KS,2)
        IF(IFL(KS).EQ.-2) WRITE(10,343) IDC(KS,2)
343        FORMAT(X,/,' NOTE: the departure from isofugacity of component ',I3,/, &
&               X,' was not used in the objective function, but',/,' it is reported anyway ',/)
!
          IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then

          WRITE(1,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2)
          WRITE(10,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2)
341          FORMAT(X,'  Temp/K    Pres/ATM  Xexp(',I0,')  Yexp(',I0,') Fug(',I0,')% Fug(',I0,')%   Zliq      Zvap ',/)
          DO I=1,NOBS(KS)

            WRITE(1,345) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),YEXP(NCXP+1),FCAUX(NCXP+1),FCAUX2(NCXP+1),ZLCALC(NCXP+1), &
&            ZVCALC(NCXP+1)
            WRITE(10,345) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),YEXP(NCXP+1),FCAUX(NCXP+1),FCAUX2(NCXP+1),ZLCALC(NCXP+1),&
&            ZVCALC(NCXP+1)
345            FORMAT(6F10.5,2F9.6)
            NCXP = NCXP + 1

          END DO
          WRITE(1,347) STDVA(KS),STDVB(KS)
          WRITE(*,347) STDVA(KS),STDVB(KS)
          WRITE(10,347) STDVA(KS),STDVB(KS)
347          FORMAT(/,X,' Std.Dev. %    ',24X,2f10.4)

        elseIF(IFL(KS).EQ.0) then

          ! Flash at T and P
          WRITE(1,351)tab,tab,tab,IDC(KS,1),tab,IDC(KS,1),tab,tab,IDC(KS,1),tab,IDC(KS,1),tab,tab,tab
351          FORMAT(/,A,'T/K   ',A,'P   /atm',A,'Xexp(',I0,')',A,'Xclc(',I0,')',A,'Err%',A,'Yexp(',I0,')',A,  &
&                  'Ycalc(',I0,')',A,'Err%',A,'Z1',A,'Z2',/)
!
          WRITE(10,3511) IDC(KS,1),IDC(KS,1)

3511          FORMAT(/,X," (x-y)(",I0,")exp,(x-y)(",I0,")calc,T/K,P/atm")
!
          ! Second counter, for graficos.csv
          NCXPy = NCXP
          DO I=1,NOBS(KS)
!
            WRITE(1,349) tab,TEXP(NCXP+1),tab,PEXP(NCXP+1),tab,XEXP(NCXP+1),tab,XCALC(NCXP+1),tab, FCAUX(NCXP+1),tab,YEXP(NCXP+1),tab,YCALC(NCXP+1),tab,FCAUX2(NCXP+1),tab, &
&                         ZLCALC(NCXP+1),tab, ZVCALC(NCXP+1)
349            FORMAT(A,2(F10.5,A),6(G10.4,A),G10.4,A,F10.8)
!
            ! x in graficos.csv
            write(10,3490)XEXP(NCXP+1),XCALC(NCXP+1),TEXP(NCXP+1),PEXP(NCXP+1)
            NCXP = NCXP + 1

          END DO
          write(10,*)
          ! y in graficos.csv
          do i = 1,Nobs(ks)

            write(10,3490)YEXP(NCXPy+1),YCALC(NCXPy+1),TEXP(NCXPy+1),PEXP(NCXPy+1)
            NCXPy = NCXPy + 1

          enddo
3490          format(2(G10.4,","),F7.2,",",F10.4)
          NO = nObs(ks)
          N1 = NCXP - NO + 1
          E1 = AAD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)
          E2 = ARD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)*100
          E3 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
          E4 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
          write(1,353)idC(ks,1),E1,E2, idC(ks,1),E3,E4
          write(*,'(/,X," ARD(x) =",f10.4,"%",20X,"ARD(y) = ",f10.4,"%")')E2,E4
          write(10,353)idC(ks,1),E1,E2, idC(ks,1),E3,E4
353          Format(/" AAD(ARD%): x(",I0,") =, ",E10.4,", (, ",F6.2,"% ,)  ,",20X," y(",I0,") = ,",E10.4,", (, ",F6.2,"% ,)")


        elseIF(IFL(KS).EQ.1) then

          ! Bubble pressure (T and x)
          NCXPy = NCXP
          WRITE(1,355)tab,tab,IDC(KS,1),tab,tab,tab,tab,IDC(KS,1),tab,IDC(KS,1),tab,tab,tab
355          FORMAT(/,A,'T/K',A,'Xexp(',I0,')',A,'Pexp/atm',A,'Pcalc/atm',A,'ERR%',A,'YEXP(',I0,')',A,'YCALC(',I0,')', &
&                   A,'ERR',A,'Z1',A,'Z2'/)
          write(10,3550) IDC(KS,1), IDC(KS,1)
3550          format(/,"T /K,(x-y)(",I0,")exp,Pexp /atm,(xexp-ycalc)(",I0,"),Pcalc /atm")
          DO I = 1, NOBS(KS)

            WRITE (1, 357) tab,TEXP(NCXP+1),tab,XEXP(NCXP+1),tab,PEXP(NCXP+1),tab,PCALC(NCXP+1),tab,FCAUX(NCXP+1),tab,YEXP(NCXP+1), &
&                           tab,YCALC(NCXP+1),tab,FCAUX2(NCXP+1),tab,ZLCALC(NCXP+1),tab,ZVCALC(NCXP+1)
357            FORMAT(A,F10.2,A,G10.4,A,G10.4,A,G10.4,A,F10.5,A,G10.5,A,F10.5,A,F10.5,A,F10.8,A,F10.8)
!
            ! x in graficos.csv
            write(10,3570)TEXP(NCXP+1),XEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),PCALC(NCXP+1)
            NCXP = NCXP + 1

          END DO
          write(10,*)
          do i = 1,NObs(ks)

            write(10,3570)TEXP(NCXPy+1),yEXP(NCXPy+1),PEXP(NCXPy+1),ycalc(NCXPy+1),PCALC(NCXPy+1)
            NCXPy = NCXPy + 1

          enddo
          NO = nObs(ks)
          N1 = NCXP - NO + 1
          E1 = AAD(PExp(N1:NCXP), PCalc(N1:NCXP), NO)
          E2 = ARD(PExp(N1:NCXP), PCalc(N1:NCXP), NO)*100
          E3 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
          E4 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
          write(1,359)E1,E2,idC(ks,1),E3,E4
          write(*,'(/,X," ARD(P) =",f10.4,"%",20X,"ARD(y) = ",f10.4,"%")')E2,E4
          write(10,359)E1,E2,idC(ks,1),E3,E4
359          Format(/" AAD(ARD%): P = ,",E10.4,", atm (, ",F6.2,"% ,)  ,",20X," y(",I0,") = ,",E10.4," (, ",F6.2,"% ,)")
3570          format(F7.2,",",F9.6,",",F10.4,",",F9.6,",",F10.4)
          ! Summary line in tables.csv (decimals chosen by magnitude)
          Pmin = 101.325*minval(Pexp(N1:NCXP)); Pmax = 101.325*maxval(Pexp(N1:NCXP))
          if (Pmin >= 10 )then
            iDim1 = 0
          elseIF (Pmin >= 1) then
            iDim1 = 1
          else
            iDim1 = 2
          endif
          if (Pmax >= 10 )then
            iDim2 = 0
          elseIF (Pmax >= 1) then
            iDim2 = 1
          else
            iDim2 = 2
          endif
          if (E2 >= 10) then
            iDim3 = 0
          elseif (E2 >= 1) then
            iDim3 = 1
          else
            iDim3 = 2
          endif
          if (E4 >= 10) then
            iDim4 = 0
          elseif (E4 >= 1) then
            iDim4 = 1
          else
            iDim4 = 2
          endif
          write (2, '(2(I2, ","), F6.0, "-", F6.0, ",", F6.'//i2s(iDim1)//', "-", F6.'//i2s(iDim2)//', ",", F6.'//i2s(iDim3)//', ",", F6.'//i2s(iDim4)//', ",", I3, ",", 80A2)')  &

              idc(ks,1), idc(ks,2), minval(Texp(N1:NCXP)),maxval(Texp(N1:NCXP)), Pmin, Pmax, E2, E4, nObs(ks), (NTEXT(ks,j), j = 1, 65)

        elseIF(IFL(KS).EQ.2) then

          ! Bubble temperature (P and x)
          NCXPy = NCXP
          WRITE(1,356) IDC(KS,1),IDC(KS,1),IDC(KS,1)
356          FORMAT(/,X,'  Pres/ATM    Xexp(',I0,')    Temp/K     TempCalc       ERR%       Yexp(',I0,')     Ycalc(',I0,')     ERR ',/)
          write(10,'(/,"P /atm,(x-y)(",I0,")exp,T /K,(xExp-yCalc)(",I0,"),Tcalc /K")')IDC(KS,1),IDC(KS,1)
          DO I=1,NOBS(KS)

            WRITE(1,358) PEXP(NCXP+1),XEXP(NCXP+1),TEXP(NCXP+1),PCALC(NCXP+1), FCAUX(NCXP+1),YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&                         ZLCALC(NCXP+1), ZVCALC(NCXP+1)
358            FORMAT(F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,'  ZLIQ = ',F10.8,'      ZVAP = ',F10.8)
!
            ! x in graficos.csv
            write(10,3580)PEXP(NCXP+1),XEXP(NCXP+1),TEXP(NCXP+1),XEXP(NCXP+1),PCALC(NCXP+1)
            NCXP = NCXP + 1

          END DO
          write(10,*)
          do i = 1,NObs(ks)

            ! y in graficos.csv
            write(10,3580)PEXP(NCXPy+1),yEXP(NCXPy+1),TEXP(NCXPy+1),ycalc(NCXPy+1),PCALC(NCXPy+1)
            NCXPy = NCXPy + 1

          enddo
          NO = nObs(ks)
          N1 = NCXP - NO + 1
          E1 = AAD(TExp(N1:NCXP), PCalc(N1:NCXP), NO)
          E2 = ARD(TExp(N1:NCXP), PCalc(N1:NCXP), NO)*100
          E3 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
          E4 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
          write(1,360)E1,E2,idC(ks,1),E3,E4
          write(*,'(/,X," ARD(T) =",f10.4,20X,"ARD(y) = ",f10.4)')E2, E4
          write(10,360)E1,E2,idC(ks,1),E3,E4
360          Format(/" AAD(ARD%): T = ,",E10.4,", K (, ",F6.2,"% ,)  ,",20X," y(",I0,") =, ",E10.4,", (, ",F6.2,"% ,)")
3580          format(F10.4,",",F9.6,",",F7.2,",",F9.6,",",F7.2)

        elseIF(iFL(KS).EQ.3)then

          ! Ternary equilibrium
          WRITE(1,451) IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
451          FORMAT(/,X,"  T/K     P/ATM   xExp(",I0,',I) xClc(',I0,',I)    ERR    xExp(',I0,',I) xClc(',I0,',I)    ERR    xExp(',I0,',II)  xClc(',I0,',II)      ERR     xExp(',I0,',II)  xClc(',I0,',II)     ERR%',/)
!
          WRITE(10,4511) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,2)
4511          FORMAT(/,X,"xExp(",I0,"),xExp(",I0,"),xClc(",I0,"),xClc(",I0,"),T /K,P /atm")
!
          ! Second counter, for graficos.csv
          NCXPy = NCXP
          DO I=1,NOBS(KS)
!
            er2 = abs(XEXP2(NCXP+1) - XCALC2(NCXP+1))
            er4 = abs(YEXP2(NCXP+1) - YCALC2(NCXP+1))
            WRITE(1,449) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),XCALC(NCXP+1),fcAux(ncXp+1),XEXP2(NCXP+1),XCALC2(NCXP+1),er2,YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&                         YEXP2(NCXP+1),YCALC2(NCXP+1),ZLCALC(NCXP+1), ZVCALC(NCXP+1)
449            FORMAT(F7.2,7F10.5,5(X,F10.5,X),2X,' ZI = ',G12.4,'       ZII =  ',G12.4)
!
            ! x in graficos.csv
            write(10,4490)XEXP(NCXP+1),xExp2(NCXP+1),xCalc(ncXP+1),xCalc2(ncXP+1),TEXP(NCXP+1),PEXP(NCXP+1),yExp(ncXp+1),yExp2(ncXp+1), yCalc(ncXP+1),yCalc2(ncXP+1)
            NCXP = NCXP + 1

          END DO
          write(10,*)
          NO = nObs(ks)
          N1 = NCXP - NO + 1
          E1 = AAD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)
          E2 = ARD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)*100
          E3 = AAD(xExp2(N1:NCXP), xCalc2(N1:NCXP), NO)
          E4 = ARD(xExp2(N1:NCXP), xCalc2(N1:NCXP), NO)*100
          E5 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
          E6 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
          E7 = AAD(yExp2(N1:NCXP), yCalc2(N1:NCXP), NO)
          E8 = ARD(yExp2(N1:NCXP), yCalc2(N1:NCXP), NO)*100   ! it used yExp/yCalc before (repeating E6)
          write(1,453)idC(ks,1),E1,E2, idC(ks,2),E3,E4, idC(ks,1),E5,E6, idC(ks,2),E7,E8
          write(*,'(/,X," ARD(x) =",f10.4,20X,"ARD(y) = ",f10.4)')E2,E4,E6,E8
          write(10,453)idC(ks,1),E1,E2, idC(ks,2),E3,E4, idC(ks,1),E5,E6, idC(ks,2),E7,E8
453          Format(/" AAD(ARD%): x(",I0,") = ,",E10.4,", ( ",F7.2,"%, ) ,  x(",I0,") = ,",E10.4," ,( ",F7.2,"%, )",20X,         &
&                    " y(",I0,") = ,",E10.4," ,(, ",F7.2,"% ,) ,  y(",I0,") = ,",E10.4," ( ",F7.2,"% )")
4490          format(4(F9.6,","),F7.2,",",F10.4,/,4(F9.6,","),/)

        END IF
!
      elseIF (IHEN(KS).EQ.1)THEN

        if(ifl(ks).EQ.1)then

          ! Excess volume

          write(1,'(/,"  T(K)       P(atm)    Xexp(",I0,")  Vexp^E(cm3/mol)   Vcalc^E(cm3/mol)       Z",/)')idC(ks,1)
          write(10,'(/,"Temp (K),P (atm),Xexp(",I0,"),Vexp^E (cm3/mol),Vcalc^E (cm3/mol),Z",/)')idC(ks,1)

        elseif(ifl(ks).EQ.2)then

          ! Excess enthalpy
          write(1,'(/,"  T(K)       P(atm)    Xexp(",I0,")  Hexp^E(J/mol)   Hcalc^E(J/mol)       Z",/)')idC(ks,1)
          write(10,'(/,"T (K),P (atm),Xexp(",I0,"),Hexp^E (J/mol),Hcalc^E (J/mol),Z",/)')idC(ks,1)

        endif
        do i = 1,nObs(ks)

          write(1,'(F7.2,2X,F10.5,2X,F10.8,X,2(F10.3,8X),F11.8)')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),prop_E_exp(NCXP+1),prop_E_calc(NCXP+1),zLiCalc(NCXP+1)
          write(10,'(F7.2,",",F10.5,",",F10.8,",",2(F10.2,","),F11.8)')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),prop_E_exp(NCXP+1),prop_E_calc(NCXP+1),zLiCalc(NCXP+1)
          NCXP = NCXP + 1

        enddo
        NO = nObs(ks)
        N1 = NCXP - NO + 1
        E1 = AAD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)  ! average absolute deviation
        E2 = ARD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)*100  ! average relative deviation
        E3 = ASD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)  ! standard deviation
        E4 = RSD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)*100 ! relative standard deviation
        if(ifl(ks).EQ.1)then
          write(1,5000)E1,E2,E3,E4
          write(10,5000)E1,E2,E3,E4
          write(*,'("(ARD(VE) = ",F10.4)')E2
        elseif(ifl(ks).EQ.2)then
          write(1,5001)E1,E2,E3,E4
          write(10,5001)E1,E2,E3,E4
          write(*,'("ARD(HE) = ",F10.4)')E2
        endif
5000        format(/" AAD(ARD%) = ",F8.4," cm3/mol (, ",F8.2,"% ,)",10X," StdDev(%) = ",F8.4," cm3/mol (, ",F8.2,"% ,)")
5001        format(/" AAD(ARD%) = ",F8.4," J/mol ,(, ",F8.2,"% ,)",10X," StdDev(%) = ",F8.4," J/mol ,(, ",F8.2,"% ,)")
!
      elseif(iHen(ks).EQ.2)then

        ! Solubility of solids in liquids
        write(1,'(/,2X,"T(K)",6X,"P(atm)",5X,"xExp(",I0,")",5x,"xCal(",I0,")",8X,"Zliq",7X,"PhiL(",I0,")",/)')idC(ks,1),idC(ks,1),idC(ks,1)
        write(10,'(/,"xExp(",I0,"),xCalc(",I0,"),T(K),P(atm)",/)')idC(ks,1),idC(ks,1)
        do i = 1, nObs(ks)

          write(1,'(F7.2,2X,F10.5,2X,2(F10.8,2X),2(F11.8,2x))')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),xCalc(NCXP+1),ZliCalc(NCXP+1),yCalc(NCXP+1)
          write(10,'(2(F10.8,","),F7.2,",",F10.5)')xExp(NCXP+1),xCalc(NCXP+1),Texp(NCXP+1),Pexp(NCXP+1)
          NCXP = NCXP + 1

        enddo
        NO = nObs(ks)
        N1 = NCXP - NO + 1
        E1 = AAD(xExp(N1:NCXP),xCalc(N1:NCXP),NO)  ! average absolute deviation
        E2 = ARD(xExp(N1:NCXP),xCalc(N1:NCXP),NO)*100  ! average relative deviation
        write(1,5002)idC(ks,1),E1,E2
        write(*,'(/"ARD(xS) = ",F10.4)')E2
        write(10,5002)idC(ks,1),E1,E2
5002        format(/" AAD(ARD%) xS(",I0,") = ,",E10.4," ,(, ",F6.2,"% ,)")

      elseif(iHen(ks) == 3)then

!         Critical points
        write(1,'(/, 2X, "Tc(K)", 5X, "Pc(atm)", 4X, "xcExp(", I0, ")", 8X, "Zc", 4X, "d[ln(f(", I0, "))]/dz(", I0, ")", 2X, "d2[ln(f(", I0, "))]/dz(", I0, ")2",  /)')idC(ks,1), idC(ks,1), idC(ks,1), idC(ks,1), idC(ks,1)
        do i = 1, nObs(ks)

          write(1,'(F7.2, 2X, F10.5, 2X, 2(F10.8, 2X), 2(F11.8, 2x))')Texp(NCXP+1), Pexp(NCXP+1), xExp(NCXP+1),  ZliCalc(NCXP+1), fcAux(NCXP+1), fcAux2(NCXP+1)
          NCXP = NCXP + 1

        enddo

      elseIF (IHEN(KS).EQ.5)THEN
!
        ! Activity coefficients
        If(mez(ks)==2)then

          write(1,1344)IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
1344          FORMAT(X,'Temp/K       Xex(',I0,')   Gex(',I0,')   Gc(',I0,')     Gex(',I0,')   Gc(',I0,')'       ,'Zliq'/)
1444          format('T/K,x',I0,'),ln(gamma(',I0,')exp),ln(gamma(',I0,'))calc,ln(gamma(',I0,')exp),ln(gamma(',I0,')calc),Z')
          write(10,1444)IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
          DO I=1,NOBS(KS)

            WRITE(1,1345)TEXP(NCXP+1),XEXP(NCXP+1),rgamexp1(ncxp+1),rgamcalc1(NCXP+1),rgamexp2(ncxp+1),rgamcalc2(NCXP+1), zlicalc(NCXP+1)
        WRITE(10,1445)TEXP(NCXP+1),XEXP(NCXP+1),dlog(rgamexp1(ncxp+1)),dlog(rgamcalc1(NCXP+1)),dlog(rgamexp2(ncxp+1)),dlog(rgamcalc2(NCXP+1)), zlicalc(NCXP+1)
1345            format(X, 2f10.2, 4f10.3,f10.3)
1445            format(F7.2, ",", 6(G10.4,","))
1347            FORMAT(2F8.4,f7.1,2F8.4,f6.3,f8.4)
             NCXP = NCXP + 1

          END DO
          NO = nObs(ks)
          N1 = NCXP - NO + 1
          E1 = AAD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)
          E2 = ARD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)*100
          E3 = AAD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)
          E4 = ARD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)*100
          write(1,1379)idC(ks,1),E1,E2,idC(ks,2),E3,E4
          write(10,1379)idC(ks,1),E1,E2,idC(ks,2),E3,E4
          write(*,'(" ARD(gamma1) = ",F10.4,20X," ARD(gamma2) = ",F10.4)')E2,E4
1379          format(/" AAD(ARD%) gamma(",I0,") = ",E10.4,", (, ",F6.2"% ,) , AAD(ARD%) gamma(",I0,") = ",E10.4,", (, ",F6.2"% ,)")

        elseif(mez(ks)==3)then

          write(1,1350) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,3),IDC(KS,3)
          write(10,1350) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,3),IDC(KS,3)
          DO I=1,NOBS(KS)

            WRITE(1,1348)TEXP(NCXP+1),XEXP(NCXP+1),xexp2(ncxp+1),rgamexp1(NCXP+1),rgamcalc1(NCXP+1),rgamexp2(NCXP+1),rgamcalc2(NCXP+1), &
&                         rgamexp3(NCXP+1),rgamcalc3(NCXP+1)
1348            FORMAT(X, f6.2 , 2(1X,e10.4), 6(1X,E12.3))
1349            FORMAT(3F9.6)
1350            FORMAT(X,'Temp/K Xex(',I0,') Xex(',I0,') Gex(',I0,') Gc(',I0,') Gex(',I0,') Gc(',I0,'),    &
&                      Gex(',I2,') Gc(',I2,')')
                NCXP = NCXP + 1

          END DO
          NO = nObs(ks)
          N1 = NCXP - NO + 1
          E1 = AAD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)
          E2 = ARD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)*100
          E3 = AAD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)
          E4 = ARD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)*100
          E5 = AAD(rGamExp3(N1:NCXP), rGamCalc3(N1:NCXP), NO)   ! it used component 2 before
          E6 = ARD(rGamExp3(N1:NCXP), rGamCalc3(N1:NCXP), NO)*100
          write(1,1351)idC(ks,1),E1,E2,idC(ks,2),E3,E4,idC(ks,3),E5,E6
          write(10,1351)idC(ks,1),E1,E2,idC(ks,2),E3,E4,idC(ks,3),E5,E6
          write(*,'(" ARD(gamma1) = ",F10.4,20X," ARD(gamma2) = ",F10.4,20X," ARD(gamma3) = ",F10.4)') E2,E4,E6
1351          format(/" AAD(ARD%) gamma(",I0,"): ",E10.4,", (, ",F6.2"%, ) , gamma(",I0,"): ",E10.4," (, ",F6.2,                  &
&                     "% ) , gamma(",I2,"): ",E10.4,", (, ",F6.2"% ,)")

        end if

      END IF

    END DO
!
  END IF
!
  CLOSE(UNIT=1)
  return
ENDsubroutine write_results
