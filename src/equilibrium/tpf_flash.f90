!--------------------------------------------------------------------------------------------------------------------------
!
!  Module tpf_flash: flash at specified vapor fraction (V/F = beta) and
!  temperature or pressure, adapted from TPfractionFlash v0.5.3 (program
!  fractionFlash, subroutines flash, Func, initialization, molar_fractions).
!  In paraest it replaces GCADB for bubble points (beta = 0).
!
!  Iteration variables, scaled logarithmically:
!     xi(i)    = ln K(i), i = 1..NC
!     xi(NC+1) = ln P (T specified) or ln T (P specified)
!  Equations: NC isofugacity conditions and sum(y) - sum(x) = 0, solved by
!  Newton-Raphson with step control [1]. The initial K come from Wilson's
!  equation and four direct-substitution steps.
!
!  Changes with respect to TPfractionFlash v0.5.3:
!  * the model is called directly (PARAGC + GCEOS of the GCA package), not
!    through THERMO;
!  * no stop: every failure returns ierr /= 0 (T above the largest Tc with T
!    specified, both phases missing at the start, no step reduction gives a
!    good evaluation, no convergence in maxIt iterations, trivial solution);
!  * an evaluation fails when the model does not find a phase (v0.5.3: only
!    when both were missing); a phase not found with a guess is tried again
!    without it (with a guess, GCEOS does not compute that phase);
!  * a liquid that ZMAX does not find because the fluid has a single root
!    (its liquid and vapor searches reach the same root, zmax_ier(2) = 199,
!    near a critical point) is accepted: that root is the liquid's;
!  * trivial solution: equal compositions and Z within 0.1% at the end
!    (v0.5.3: only |ZV - ZL| < 1e-6, which let trivial solutions through);
!  * no permanent gases, heavy compounds or traces (paraest gives every
!    component a positive amount), no log or result files.
!
!  [1] ML Michelsen, JM Mollerup, "Thermodynamic Models: Fundamentals &
!      Computational Aspects", 2nd Ed., Tie-Line Publications, 2007,
!      pp. 302-303.
!
!--------------------------------------------------------------------------------------------------------------------------
module tpf_flash
  use zmax_status, only: zmax_ier
  implicit none
  private
  public :: flash_beta

  integer, parameter :: NCM = 30          ! must match the GCA package
  integer, parameter :: V = 1, L = 2
  ! Bound of ln K: K up to 1e43, and sum(y) cannot overflow
  double precision, parameter :: lnKmax = 100D0

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  spec = 'T': T given, P calculated;  'P': P given, T calculated.
!  T and P: the values of the specified variable and an estimate of the other
!  (it is the starting point when spec = 'T'); on return, the result.
!  feed(NC) > 0, beta in [0, 1]. x, y: liquid and vapor mole fractions,
!  ZL, ZV: their compressibility factors.
!  ierr = 0  converged
!         1  T above the largest critical temperature (spec = 'T')
!         2  neither phase found at the initial conditions
!         3  no step reduction gives an evaluation with both phases
!         4  not converged in maxIt iterations
!         5  trivial solution (equal compositions and Z)
!
subroutine flash_beta(NC, NG, NST, spec, beta, T, P, feed, x, y, ZL, ZV, ierr)
  integer, intent(in) :: NC, NG, NST
  character, intent(in) :: spec
  double precision, intent(in) :: beta, feed(NC)
  double precision, intent(inout) :: T, P
  double precision, intent(out) :: x(NC), y(NC), ZL, ZV
  integer, intent(out) :: ierr

  integer, parameter :: maxIt = 50, maxAlphaIt = 10
  double precision, parameter :: tolerance = 1D-8, maxStepTP = 0.25D0

  integer :: iteration, alphaIt, guess, i, NDim, pivot(NC+1)
  double precision :: K(NC), xi(NC+1), xiNew(NC+1), dxi(NC+1), dxiOld(NC+1), F(NC+1), jac(NC+1,NC+1)
  double precision :: w(NCM,2), Z(2), lnPhi(NCM,2), alpha, alphaOld, d, err, errNew, Pold
  logical :: thermoOK

  NDim = NC + 1
  Z = 0D0
  w = 0D0
  call initialization(NC, NG, NST, spec, beta, T, P, feed, K, ierr)
  if (ierr /= 0) return

  xi(:NC) = log(K)
  if (spec == 'P') then
    xi(NDim) = log(T)
  else
    xi(NDim) = log(P)
  endif
  guess = 0
  alphaOld = 1D0
  dxi = 0D0
  call func(NC, NG, NST, NDim, spec, beta, T, P, feed, K, xi, guess, w, Z, lnPhi, F, jac, err, thermoOK)
  if (.not. thermoOK) then
    ierr = 3
    return
  endif

  ierr = 4
  do iteration = 1, maxIt
    if (err < tolerance) then
      ierr = 0
      exit
    endif
    ! Newton step: J dxi = -F
    dxiOld = dxi
    call LUDcmp(jac(:NDim,:NDim), NDim, NDim, pivot, d)
    dxi = -F
    call LUBksb(jac(:NDim,:NDim), NDim, NDim, pivot, dxi)
    if (any(dxi /= dxi)) then             ! singular Jacobian (LUDcmp returns NaN)
      ierr = 3
      return
    endif
    ! Limit the step in ln T or ln P, which go to the model as they are
    if (abs(dxi(NDim)) < maxStepTP) then
      alpha = 1D0
    else
      alpha = maxStepTP/abs(dxi(NDim))
    endif
    ! Halve the step when it nearly cancels the previous one (oscillation)
    if (norm2(dxi*alpha + dxiOld*alphaOld)/err < 0.1D0) alpha = alpha/2
    ! Cut the step until the model finds both phases and the error does not
    ! grow more than tenfold
    alphaIt = 0
    do
      alphaIt = alphaIt + 1
      xiNew = xi + dxi*alpha
      xiNew(:NC) = max(-lnKmax, min(lnKmax, xiNew(:NC)))
      K = exp(xiNew(:NC))
      Pold = P
      if (spec == 'T') then
        P = exp(xiNew(NDim))
      else
        T = exp(xiNew(NDim))
      endif
      ! After 5 iterations, or once the error is small, start the volume
      ! from the previous one (scaled with P when P changes)
      if ((iteration > 5 .or. err < 1D0) .and. thermoOK) then
        guess = 2
        if (spec == 'T') Z(L) = Z(L)*P/Pold
      endif
      call func(NC, NG, NST, NDim, spec, beta, T, P, feed, K, xiNew, guess, w, Z, lnPhi, F, jac, errNew, thermoOK)
      if (thermoOK .and. errNew < 10D0*err) then
        xi = xiNew
        err = errNew
        alphaOld = alpha
        exit
      endif
      if (alphaIt > maxAlphaIt) then
        ierr = 3
        return
      endif
      alpha = alpha/3.1622777D0           ! one order of magnitude every two cuts
      if (.not. thermoOK) guess = 0
    enddo
  enddo
  if (ierr /= 0) return

  ! Trivial solution: both phases the same, in composition (|ln K| < 0.01)
  ! and in Z (within 0.1%). An azeotrope has x = y but different Z; near a
  ! critical point Z are close but the compositions differ.
  if (maxval(abs(xi(:NC))) < 1D-2 .and. abs(Z(V) - Z(L)) <= 1D-3*max(Z(V), Z(L))) then
    ierr = 5
    return
  endif
  x = w(:NC,L)
  y = w(:NC,V)
  ZL = Z(L)
  ZV = Z(V)
endsubroutine flash_beta

!--------------------------------------------------------------------------------------------------------------------------
!
!  Isofugacity conditions and the mole-fraction sum (F), the Jacobian with
!  respect to xi, the phase compositions w and the model results.
!  thermoOK is false when the model does not find a phase or both phases are
!  the same.
!
subroutine func(NC, NG, NST, NDim, spec, beta, T, P, feed, K, xi, guess, w, Z, lnPhi, F, jac, err, thermoOK)
  integer, intent(in) :: NC, NG, NST, NDim, guess
  character, intent(in) :: spec
  double precision, intent(in) :: beta, T, P, feed(NC), K(NC), xi(NDim)
  double precision, intent(inout) :: w(NCM,2), Z(2)
  double precision, intent(out) :: lnPhi(NCM,2), F(NDim), jac(NDim,NDim), err
  logical, intent(out) :: thermoOK

  double precision :: dLnPhi_dn(NCM,NCM,2), dLnPhi_dT(NCM,2), dLnPhi_dP(NCM,2), denom(NC), dxdxi(NC), dydxi(NC)
  integer :: ic(2), i, j

  call molar_fractions(NC, beta, K, feed, w, denom)
  call PARAGC(T, NC, NG, NST, 1)
  call GCEOS(NC, NG, NST, 2, 1, T, P, w(:,L), lnPhi(:,L), dLnPhi_dn(:,:,L), dLnPhi_dT(:,L), dLnPhi_dP(:,L), Z(L), &
             guess, 1, ic(L))
  ! With a guess, GCEOS does not compute a phase it does not find: again
  ! without guess
  if (ic(L) /= 1 .and. guess /= 0) &
    call GCEOS(NC, NG, NST, 2, 1, T, P, w(:,L), lnPhi(:,L), dLnPhi_dn(:,:,L), dLnPhi_dT(:,L), dLnPhi_dP(:,L), Z(L), &
               0, 1, ic(L))
  ! Liquid "not found" because the liquid and vapor searches reached the same
  ! root: the fluid has a single root (near a critical point), and it is the
  ! liquid's
  if (ic(L) /= 1 .and. zmax_ier(2) == 199) ic(L) = 1
  call GCEOS(NC, NG, NST, 2, 1, T, P, w(:,V), lnPhi(:,V), dLnPhi_dn(:,:,V), dLnPhi_dT(:,V), dLnPhi_dP(:,V), Z(V), &
             guess, -1, ic(V))
  thermoOK = all(ic == 1) .and. abs(Z(V) - Z(L)) >= 1D-6
  if (.not. thermoOK) then
    F = 0D0
    jac = 0D0
    err = huge(1D0)
    return
  endif
  lnPhi(:NC,:) = min(lnPhi(:NC,:), 308D0)

  do j = 1, NC
    dxdxi(j) = -w(j,V)*beta/denom(j)
    dydxi(j) = w(j,V) + dxdxi(j)*K(j)
    do i = 1, NC
      ! Derivatives with respect to the amounts (L x(i) and V y(i))
      jac(i,j) = dLnPhi_dn(i,j,V)*dydxi(j) - dLnPhi_dn(i,j,L)*dxdxi(j)
    enddo
    jac(j,j) = jac(j,j) + 1D0
    F(j) = xi(j) + lnPhi(j,V) - lnPhi(j,L)
    if (spec == 'T') then
      jac(j,NDim) = (dLnPhi_dP(j,V) - dLnPhi_dP(j,L))*P
    else
      jac(j,NDim) = (dLnPhi_dT(j,V) - dLnPhi_dT(j,L))*T
    endif
    jac(NDim,j) = dydxi(j) - dxdxi(j)
  enddo
  F(NDim) = sum(w(:NC,V) - w(:NC,L))
  jac(NDim,NDim) = 0D0
  err = norm2(F)
endsubroutine func

!--------------------------------------------------------------------------------------------------------------------------
!
!  Initial K: Wilson's equation, then four direct-substitution steps that
!  also correct P (T specified) or T (P specified; ln K of a reference
!  component taken linear in 1/T).
!
subroutine initialization(NC, NG, NST, spec, beta, T, P, feed, K, ierr)
  integer, intent(in) :: NC, NG, NST
  character, intent(in) :: spec
  double precision, intent(in) :: beta, feed(NC)
  double precision, intent(inout) :: T, P
  double precision, intent(out) :: K(NC)
  integer, intent(out) :: ierr

  double precision :: Tc(NCM), Pc(NCM), omega(NCM)
  common /CRIT/ Tc, Pc, omega

  double precision :: w(NCM,2), Z(2), lnPhi(NCM,2), dLnPhi_dn(NCM,NCM), dLnPhi_dT(NCM,2), dLnPhi_dP(NCM)
  double precision :: denom(NC), A, B, lnKR, best, sumL, sumV
  integer :: ic(2), i, iteration, ref

  ierr = 0
  B = 0D0
  ref = 1
  do i = 1, NC
    K(i) = exp(max(-lnKmax, min(lnKmax, -log(P/Pc(i)) + 5.373D0*(1D0 + omega(i))*(1D0 - Tc(i)/T))))
  enddo
  if (T > maxval(Tc(:NC))) then
    if (spec == 'T') then
      ierr = 1
      return
    endif
    T = 0.9D0*maxval(Tc(:NC))
  endif
  if (spec == 'P') then
    ! Reference component: the most abundant one, or, if it is supercritical,
    ! the one whose ln K is closest to the middle of the range
    ref = maxloc(feed, 1)
    if (Tc(ref) < T) then
      lnKR = (log(maxval(K)) + log(minval(K)))/2
      best = log(maxval(K))
      do i = 1, NC
        if (abs(lnKR - log(K(i))) < best) then
          ref = i
          best = abs(lnKR - log(K(i)))
        endif
      enddo
    endif
  endif

  w = 0D0
  Z = 0D0
  call molar_fractions(NC, beta, K, feed, w, denom)
  do iteration = 1, 4
    call PARAGC(T, NC, NG, NST, 1)
    call GCEOS(NC, NG, NST, 1, 1, T, P, w(:,L), lnPhi(:,L), dLnPhi_dn, dLnPhi_dT(:,L), dLnPhi_dP, Z(L), 0, 1, ic(L))
    call GCEOS(NC, NG, NST, 1, 1, T, P, w(:,V), lnPhi(:,V), dLnPhi_dn, dLnPhi_dT(:,V), dLnPhi_dP, Z(V), 0, -1, ic(V))
    if (ic(L) < 0 .and. ic(V) < 0) then
      ierr = 2
      return
    elseif (ic(L) > 0 .and. ic(V) < 0) then
      lnPhi(:NC,V) = 0D0                  ! ideal gas
    endif
    ! Slope of ln K(ref) = A + B/T
    if (iteration < 3 .and. spec == 'P') B = -T*T*(dLnPhi_dT(ref,L) - dLnPhi_dT(ref,V))
    do i = 1, NC
      K(i) = exp(max(-lnKmax, min(lnKmax, lnPhi(i,L) - lnPhi(i,V))))
    enddo
    call molar_fractions(NC, beta, K, feed, w, denom)
    sumV = log(sum(w(:NC,V)))
    sumL = -log(sum(w(:NC,L)))
    if (spec == 'T') then
      P = exp(log(P) + sumV + sumL)
    else
      lnKR = log(K(ref))
      A = lnKR - B/T
      lnKR = lnKR - (sumV + sumL)
      T = B/(lnKR - A)
    endif
  enddo
  if (.not. (T > 0D0 .and. P > 0D0 .and. T < huge(T) .and. P < huge(P))) ierr = 2
endsubroutine initialization

!--------------------------------------------------------------------------------------------------------------------------
!
!  Phase compositions (not normalized during the iterations) from the feed,
!  beta and K: x = z/(1 - beta + beta K), y = K x.
!
subroutine molar_fractions(NC, beta, K, z, w, denom)
  integer, intent(in) :: NC
  double precision, intent(in) :: beta, K(NC), z(NC)
  double precision, intent(inout) :: w(NCM,2)
  double precision, intent(out) :: denom(NC)
  denom = 1D0 - beta + beta*K
  w(:NC,L) = z/denom
  w(:NC,V) = K*w(:NC,L)
endsubroutine molar_fractions

endmodule tpf_flash
