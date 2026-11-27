!--------------------------------------------------------------------------------------------------------------------------
!
!  Module rkpr_bridge: link between paraest and the RK-PR package (rkpr.f90).
!
!  readcomp reads the parameters of all the components of paraest.dat and
!  leaves them in the package's COMMON blocks. load_from_readcomp keeps a
!  copy of them here (the whole system); rkpr_select puts back in the COMMON
!  blocks those of the components of a data set (or of one component), in
!  the order given, as fakePARGR/fakePARMOL do for the GCA.
!
!  Units: as inside the package (K, atm, cm3, mol). write_rkpr_params writes
!  them as readcomp reads them (bar, L).
!
!  Fitted parameters (read_rkpr_fit, rkpr_set_params): Tc, Pc and delta1 of
!  a component instead of ac and b, which follow from them (as readcomp does
!  for input j = 1, inverted); k; k_ij(inf), with k'_ij and T*_ij when its
!  temperature dependence is fitted too; l_ij. The binary parameters are
!  symmetric: (i,j) also sets (j,i).
!
!--------------------------------------------------------------------------------------------------------------------------
module rkpr_bridge
  use dimensions, only: NPTIPOM, NIDM
  implicit none
  private
  public :: load_from_readcomp, rkpr_select, write_rkpr_params
  public :: read_rkpr_fit, rkpr_fit_bounds, rkpr_set_params, write_rkpr_fit_spec, rkpr_critical

  integer, parameter :: NCM = 30                     ! must match the RK-PR package
  double precision, parameter :: Rgas = 0.08314472D0, Rgas2 = 82.05D0
  double precision, parameter :: atm_to_bar = Rgas/Rgas2*1D3   ! as readcomp converts Pc
  ! Constants of readcomp: k from omega and Zc
  double precision, parameter :: Apar(2) = [0.0017D0, -2.4407D0], Bpar(2) = [1.9681D0, 7.4513D0], &
                                 Cpar(2) = [-2.7238D0, 12.504D0]
  ! LM bounds (module config has the PSO boxes)
  double precision, parameter :: kBounds(2) = [0.1D0, 10D0]    ! k
  double precision, parameter :: kpMax = 1D0                   ! k' between -kpMax and kpMax
  double precision, parameter :: kijBounds(2) = [-29D0, 0.8D0] ! k_ij(inf), l_ij: 1 - (30 to 0.2), as 1 - k* of the GCA

  ! Whole system, as readcomp left it
  integer :: nSys = 0, sysNComb = 0, sysNTDEP = 0
  double precision :: sysAc(NCM), sysB(NCM), sysDel1(NCM), sysK(NCM)
  double precision :: sysKinf(NCM,NCM), sysKp(NCM,NCM), sysTstar(NCM,NCM), sysBij(NCM,NCM)
  double precision :: sysTc(NCM), sysPc(NCM), sysOm(NCM), sysDc(NCM)
  double precision :: sysPM(NCM), sysPen(NCM), sysH(7,NCM)
  character(10) :: sysName(NCM)
  double precision :: sysL(NCM,NCM) = 0D0

  ! Fitted parameters: components of Tc, Pc, delta1 and k; pairs of
  ! k_ij(inf) (iTempK: 0 alone, 1 with k', 2 with k' and T*) and of l_ij,
  ! each with the pairs that copy its values ("identical").
  integer :: nTcF = 0, nPcF = 0, nDelF = 0, nKF = 0, nKijF = 0, nLijF = 0
  integer :: idTcF(NPTIPOM), idPcF(NPTIPOM), idDelF(NPTIPOM), idKF(NPTIPOM)
  integer :: ijK(NPTIPOM,2), iTempK(NPTIPOM), nIdK(NPTIPOM), idIdK(NIDM,2)
  integer :: ijL(NPTIPOM,2), nIdL(NPTIPOM), idIdL(NIDM,2)

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  After readcomp: copies the parameters of the NC components to this
!  module, and returns Tc and Pc (atm).
!
subroutine load_from_readcomp(NC, Tc_out, Pc_out)
  integer, intent(in) :: NC
  double precision, intent(out) :: Tc_out(NC), Pc_out(NC)

  double precision :: ac, b, del1, rk, Kijinf, Kijp, Tstar, bij, TC, PC, om, DCeos
  double precision :: PM, Pen, HHA, HHB, HHC, HHD, HHE, HHF, HHG
  integer :: NTDEP, ncomb, i, j
  character(10) :: fluid
  COMMON /COMPONENTS/ ac(NCM), b(NCM), del1(NCM), rk(NCM), Kijinf(NCM,NCM), NTDEP
  COMMON /rule/ ncomb
  COMMON /Tdep/ Kijp(NCM,NCM), Tstar(NCM,NCM)
  COMMON /bcross/ bij(NCM,NCM)
  COMMON /CRIT/ TC(NCM), PC(NCM), om(NCM), DCeos(NCM)
  COMMON /name/ fluid(NCM)
  common /GCPROM/ PM(NCM), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM), HHG(NCM)

  nSys = NC
  sysNComb = ncomb
  sysNTDEP = NTDEP
  sysAc(:NC) = ac(:NC)
  sysB(:NC) = b(:NC)
  sysDel1(:NC) = del1(:NC)
  sysK(:NC) = rk(:NC)
  sysKinf(:NC,:NC) = Kijinf(:NC,:NC)
  sysKp(:NC,:NC) = Kijp(:NC,:NC)
  sysTstar(:NC,:NC) = Tstar(:NC,:NC)
  sysBij(:NC,:NC) = bij(:NC,:NC)
  sysTc(:NC) = TC(:NC)
  sysPc(:NC) = PC(:NC)
  sysOm(:NC) = om(:NC)
  sysDc(:NC) = DCeos(:NC)
  sysName(:NC) = fluid(:NC)
  sysPM(:NC) = PM(:NC)
  sysPen(:NC) = Pen(:NC)
  sysH(1,:NC) = HHA(:NC)
  sysH(2,:NC) = HHB(:NC)
  sysH(3,:NC) = HHC(:NC)
  sysH(4,:NC) = HHD(:NC)
  sysH(5,:NC) = HHE(:NC)
  sysH(6,:NC) = HHF(:NC)
  sysH(7,:NC) = HHG(:NC)
  sysL = 0D0
  do i = 1, NC
    do j = 1, NC
      if (i == j) cycle
      if (ncomb == 0) then
        sysL(i,j) = 1D0 - 2D0*bij(i,j)/(b(i) + b(j))
      else
        sysL(i,j) = 1D0 - 2D0*bij(i,j)**(1D0/3D0)/(b(i)**(1D0/3D0) + b(j)**(1D0/3D0))
      endif
    enddo
  enddo
  Tc_out = TC(:NC)
  Pc_out = PC(:NC)
endsubroutine load_from_readcomp

!--------------------------------------------------------------------------------------------------------------------------
!
!  Loads the package with the components idx(1..NCw) of the system, in
!  that order.
!
subroutine rkpr_select(NCw, idx)
  integer, intent(in) :: NCw, idx(NCw)

  double precision :: ac, b, del1, rk, Kijinf, Kijp, Tstar, bij, TC, PC, om, DCeos
  double precision :: PM, Pen, HHA, HHB, HHC, HHD, HHE, HHF, HHG
  integer :: NTDEP, ncomb, iEntalp, i, j
  character(10) :: fluid
  COMMON /COMPONENTS/ ac(NCM), b(NCM), del1(NCM), rk(NCM), Kijinf(NCM,NCM), NTDEP
  COMMON /rule/ ncomb
  COMMON /Tdep/ Kijp(NCM,NCM), Tstar(NCM,NCM)
  COMMON /bcross/ bij(NCM,NCM)
  COMMON /CRIT/ TC(NCM), PC(NCM), om(NCM), DCeos(NCM)
  COMMON /name/ fluid(NCM)
  common /GCPROM/ PM(NCM), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM), HHG(NCM)
  common /ENTALP/ iEntalp

  ncomb = sysNComb
  NTDEP = sysNTDEP
  do i = 1, NCw
    ac(i) = sysAc(idx(i))
    b(i) = sysB(idx(i))
    del1(i) = sysDel1(idx(i))
    rk(i) = sysK(idx(i))
    TC(i) = sysTc(idx(i))
    PC(i) = sysPc(idx(i))
    om(i) = sysOm(idx(i))
    DCeos(i) = sysDc(idx(i))
    fluid(i) = sysName(idx(i))
    PM(i) = sysPM(idx(i))
    Pen(i) = sysPen(idx(i))
    HHA(i) = sysH(1,idx(i))
    HHB(i) = sysH(2,idx(i))
    HHC(i) = sysH(3,idx(i))
    HHD(i) = sysH(4,idx(i))
    HHE(i) = sysH(5,idx(i))
    HHF(i) = sysH(6,idx(i))
    HHG(i) = sysH(7,idx(i))
    do j = 1, NCw
      Kijinf(i,j) = sysKinf(idx(i),idx(j))
      Kijp(i,j) = sysKp(idx(i),idx(j))
      Tstar(i,j) = sysTstar(idx(i),idx(j))
      bij(i,j) = sysBij(idx(i),idx(j))
    enddo
  enddo
  ! as readcomp: number of components with ideal-gas enthalpy constants
  iEntalp = count(any(sysH(1:6,idx) /= 0D0, dim=1))
endsubroutine rkpr_select

!--------------------------------------------------------------------------------------------------------------------------
!
!  Writes the parameters of the system to unit u, in the units readcomp
!  reads them (bar, L/mol).
!
subroutine write_rkpr_params(u)
  integer, intent(in) :: u
  integer :: i, j

  write(u,'(/,X,"RK-PR parameters",/,/,3X,"Compound",5X,"Tc/K",6X,"Pc/bar",4X,"omega",3X, &
    & "ac/(bar L2/mol2)",3X,"b/(L/mol)",5X,"delta1",6X,"k")')
  do i = 1, nSys
    write(u,'(3X,A10,X,F8.2,2X,F9.4,2X,F7.4,2X,G15.8,2X,G13.6,2X,F9.5,2X,F9.5)') sysName(i), sysTc(i), &
      sysPc(i)*Rgas/Rgas2*1D3, sysOm(i), sysAc(i)/(Rgas2/Rgas*1D3), sysB(i)/1D3, sysDel1(i), sysK(i)
  enddo
  write(u,'(/,X,"Binary parameters (combining rule ",I1,", T dependence ",I1,")",/,/, &
    & 3X,"i",2X,"j",6X,"k_inf",8X,"k''",9X,"T*/K",8X,"l")') sysNComb, sysNTDEP
  do i = 1, nSys
    do j = i + 1, nSys
      write(u,'(X,I3,X,I3,4(2X,G11.4))') i, j, sysKinf(i,j), sysKp(i,j), sysTstar(i,j), sysL(i,j)
    enddo
  enddo
endsubroutine write_rkpr_params

!--------------------------------------------------------------------------------------------------------------------------
!
!  Reads from unit u the fitted parameters, in the format of those of the
!  GCA (after the data sets):
!
!    NTc NPc Ndelta Nk Nkij Nlij
!    components whose Tc is fitted            (if NTc > 0; the same for
!    ... Pc, delta1, k                          Pc, delta1 and k)
!    i j  i j ...                             pairs of k_ij(inf)  } if
!    T dependence of each: 0, 1 (k'), 2 (k', T*)                  } Nkij
!    number of identical pairs of each                            } > 0
!    i j  i j ...                             one line per pair with identical ones
!    i j  i j ...                             pairs of l_ij       } if
!    number of identical pairs of each                            } Nlij
!    i j  i j ...                             one line per pair with identical ones
!    initial values                           Tc/K, Pc/bar, delta1, k, then for
!                                             each k_ij: k_inf [k' [T*/K]], then l_ij
!
!  The binary parameters are symmetric: (j,i) is the same as (i,j), and
!  need not be listed as identical. An initial value 0 of Tc, Pc, delta1, k
!  or T* takes the one of the parameters read (T* without temperature
!  dependence in the input: min(Tc_i, Tc_j)).
!  Returns the number of parameters, their initial values (XX) and types
!  (IDXX, see module scaling), and sets their bounds (module scaling).
!
subroutine read_rkpr_fit(u, NPAR, XX, IDXX)
  integer, intent(in) :: u
  integer, intent(out) :: NPAR, IDXX(*)
  double precision, intent(inout) :: XX(*)
  integer :: i, j, k, m, n, nId
  double precision :: TcOut(NCM), PcOut(NCM)

  read(u,*) nTcF, nPcF, nDelF, nKF, nKijF, nLijF
  write(*,'(6(I2,4X),"! Tc - Pc - delta1 - k - k_ij - l_ij")') nTcF, nPcF, nDelF, nKF, nKijF, nLijF
  if (any([nTcF, nPcF, nDelF, nKF, nKijF, nLijF] < 0) .or. any([nTcF, nPcF, nDelF, nKF, nKijF, nLijF] > NPTIPOM)) then
    write(*,'(A,I0)') ' ERROR: the number of fitted parameters of each type must be between 0 and ', NPTIPOM
    stop 1
  endif
  call read_components(nTcF, idTcF, 'Tc')
  call read_components(nPcF, idPcF, 'Pc')
  call read_components(nDelF, idDelF, 'delta1')
  call read_components(nKF, idKF, 'k')
  if (nKijF > 0) then
    call read_pairs(nKijF, ijK, 'k_ij')
    read(u,*) iTempK(:nKijF)
    write(*,'(*(I2,6X))') iTempK(:nKijF)
    if (any(iTempK(:nKijF) < 0) .or. any(iTempK(:nKijF) > 2)) then
      write(*,'(A)') " ERROR: temperature dependence of k_ij: 0 (k_inf), 1 (k_inf, k') or 2 (k_inf, k', T*)"
      stop 1
    endif
    call read_identical(nKijF, ijK, nIdK, idIdK, 'k_ij')
  endif
  if (nLijF > 0) then
    call read_pairs(nLijF, ijL, 'l_ij')
    call read_identical(nLijF, ijL, nIdL, idIdL, 'l_ij')
  endif

  ! Types
  NPAR = 0
  call add_type(nTcF, 11)
  call add_type(nPcF, 12)
  call add_type(nDelF, 13)
  call add_type(nKF, 14)
  do m = 1, nKijF
    call add_type(1, 15)
    if (iTempK(m) >= 1) call add_type(1, 16)
    if (iTempK(m) == 2) call add_type(1, 17)
  enddo
  call add_type(nLijF, 18)
  if (NPAR == 0) return

  read(u,*) XX(:NPAR)
  write(*,'(*(G14.6))') XX(:NPAR)

  ! k(T) needs the temperature dependence of the package
  if (any(iTempK(:nKijF) > 0) .and. sysNTDEP == 0) then
    sysNTDEP = 1
    sysKp(:nSys,:nSys) = 0D0
    do i = 1, nSys
      do j = 1, nSys
        sysTstar(i,j) = min(sysTc(i), sysTc(j))
      enddo
    enddo
    write(*,'(A)') " Temperature dependence of k_ij turned on (Tdep = 1): k' = 0, T* = min(Tc_i, Tc_j) for the other pairs"
  endif

  ! Pairs fitted without temperature dependence (flag 0): k' = 0, whatever
  ! the file says (and for their identical pairs); T* only has to be /= 0.
  nId = 0
  do m = 1, nKijF
    if (iTempK(m) == 0) then
      i = ijK(m,1)
      j = ijK(m,2)
      call set_pair(sysKp, i, j, 0D0)
      if (sysTstar(i,j) == 0D0) call set_pair(sysTstar, i, j, min(sysTc(i), sysTc(j)))
      do k = nId + 1, nId + nIdK(m)
        call set_pair(sysKp, idIdK(k,1), idIdK(k,2), 0D0)
        call set_pair(sysTstar, idIdK(k,1), idIdK(k,2), sysTstar(i,j))
      enddo
    endif
    nId = nId + nIdK(m)
  enddo

  ! Initial values 0: those read
  n = 0
  do m = 1, nTcF
    n = n + 1
    if (XX(n) == 0D0) XX(n) = sysTc(idTcF(m))
  enddo
  do m = 1, nPcF
    n = n + 1
    if (XX(n) == 0D0) XX(n) = sysPc(idPcF(m))*atm_to_bar
  enddo
  do m = 1, nDelF
    n = n + 1
    if (XX(n) == 0D0) XX(n) = sysDel1(idDelF(m))
  enddo
  do m = 1, nKF
    n = n + 1
    if (XX(n) == 0D0) XX(n) = sysK(idKF(m))
  enddo
  do m = 1, nKijF
    n = n + 1 + min(iTempK(m), 1)
    if (iTempK(m) == 2) then
      n = n + 1
      if (XX(n) == 0D0) XX(n) = sysTstar(ijK(m,1), ijK(m,2))
      if (XX(n) <= 0D0) then
        write(*,'(A,2(I0,A))') ' ERROR: T* of pair (', ijK(m,1), ',', ijK(m,2), ') must be positive'
        stop 1
      endif
    endif
  enddo

  ! The initial values in the model: the bounds of T* use them (Tc)
  call rkpr_set_params(XX, TcOut, PcOut)

contains

  subroutine read_components(nF, idF, what)
    integer, intent(in) :: nF
    integer, intent(out) :: idF(:)
    character(*), intent(in) :: what
    integer :: k
    if (nF == 0) return
    read(u,*) idF(:nF)
    write(*,'(*(I3))') idF(:nF)
    do k = 1, nF
      if (idF(k) < 1 .or. idF(k) > nSys) then
        write(*,'(3A,I0,A,I0)') ' ERROR in the fitted parameters (', what, '): component ', idF(k), &
          ' is not between 1 and ', nSys
        stop 1
      endif
      if (any(idF(:k-1) == idF(k))) then
        write(*,'(3A,I0,A)') ' ERROR in the fitted parameters (', what, '): component ', idF(k), ' listed twice'
        stop 1
      endif
    enddo
  endsubroutine read_components

  ! Pairs (i,j), i /= j; kept as i < j (symmetric parameters)
  subroutine read_pairs(nF, ij, what)
    integer, intent(in) :: nF
    integer, intent(out) :: ij(:,:)
    character(*), intent(in) :: what
    integer :: k
    read(u,*) (ij(k,1), ij(k,2), k = 1, nF)
    write(*,'(*(2(I2,X),4X))') (ij(k,1), ij(k,2), k = 1, nF)
    do k = 1, nF
      call check_pair(ij(k,:), what)
      if (any(ij(:k-1,1) == ij(k,1) .and. ij(:k-1,2) == ij(k,2))) then
        write(*,'(3A,2(I0,A))') ' ERROR in the fitted parameters (', what, '): pair (', ij(k,1), ',', ij(k,2), &
          ') listed twice (the binary parameters are symmetric)'
        stop 1
      endif
    enddo
  endsubroutine read_pairs

  subroutine check_pair(p, what)
    integer, intent(inout) :: p(2)
    character(*), intent(in) :: what
    if (p(1) < 1 .or. p(1) > nSys .or. p(2) < 1 .or. p(2) > nSys .or. p(1) == p(2)) then
      write(*,'(3A,2(I0,A),I0)') ' ERROR in the fitted parameters (', what, '): pair (', p(1), ',', p(2), &
        ') is not two different components between 1 and ', nSys
      stop 1
    endif
    p = [minval(p), maxval(p)]
  endsubroutine check_pair

  subroutine read_identical(nF, ij, nIdF, idId, what)
    integer, intent(in) :: nF, ij(:,:)
    integer, intent(out) :: nIdF(:), idId(:,:)
    character(*), intent(in) :: what
    integer :: k, l, nS
    read(u,*) nIdF(:nF)
    write(*,'(*(I2,6X))') nIdF(:nF)
    if (any(nIdF(:nF) < 0) .or. sum(nIdF(:nF)) > NIDM) then
      write(*,'(3A,I0)') ' ERROR in the fitted parameters (', what, '): identical pairs must be between 0 and ', NIDM
      stop 1
    endif
    nS = 0
    do k = 1, nF
      if (nIdF(k) == 0) cycle
      read(u,*) (idId(l,1), idId(l,2), l = nS + 1, nS + nIdF(k))
      write(*,'(*(2(I2,X),4X))') (idId(l,1), idId(l,2), l = nS + 1, nS + nIdF(k))
      do l = nS + 1, nS + nIdF(k)
        call check_pair(idId(l,:), what)
        if (all(idId(l,:) == ij(k,:))) write(*,'(3A,2(I0,A))') ' WARNING (', what, '): pair (', ij(k,1), ',', &
          ij(k,2), ') listed as identical to itself (the binary parameters are symmetric): ignored'
      enddo
      nS = nS + nIdF(k)
    enddo
  endsubroutine read_identical

  subroutine add_type(k, t)
    integer, intent(in) :: k, t
    IDXX(NPAR+1:NPAR+k) = t
    NPAR = NPAR + k
  endsubroutine add_type

endsubroutine read_rkpr_fit

!--------------------------------------------------------------------------------------------------------------------------
!
!  LM bounds (pMin, pMax) and PSO box (pLo, pHi) of the parameters read by
!  read_rkpr_fit, from their initial values XX (unscaled). Called after the
!  options of paraest.cfg have been read (they give some of the bounds).
!
subroutine rkpr_fit_bounds(NPAR, XX, IDXX)
  use scaling, only: pMin, pMax, pLo, pHi
  use config, only: range_tc, range_pc, range_k, range_kp, range_delta1, range_kstr, lm_range_tc, lm_range_pc
  integer, intent(in) :: NPAR, IDXX(*)
  double precision, intent(in) :: XX(*)
  integer :: m, n
  double precision :: x0, Tmin, Tmax

  ! Bounds (LM) and PSO box of each parameter
  do n = 1, NPAR
    x0 = XX(n)
    select case (IDXX(n))
    case (11)
      call set_bounds(n, lm_range_tc*x0, range_tc*x0)
    case (12)
      call set_bounds(n, lm_range_pc*x0, range_pc*x0)
    case (13)
      call set_bounds(n, range_delta1, range_delta1)
    case (14)
      call set_bounds(n, kBounds, range_k*x0)
    case (15, 18)
      call set_bounds(n, kijBounds, 1D0 - [range_kstr(2), range_kstr(1)])
    case (16)
      call set_bounds(n, [-kpMax, kpMax], range_kp)
    endselect
  enddo
  n = nTcF + nPcF + nDelF + nKF
  do m = 1, nKijF
    n = n + 1 + min(iTempK(m), 1)
    if (iTempK(m) == 2) then
      n = n + 1
      Tmin = min(sysTc(ijK(m,1)), sysTc(ijK(m,2)))
      Tmax = max(sysTc(ijK(m,1)), sysTc(ijK(m,2)))
      call set_bounds(n, [0.1D0*Tmin, 2D0*Tmax], [0.5D0*Tmin, 2D0*Tmin])
    endif
  enddo

contains

  ! LM bounds b and PSO box r of parameter k. The bounds are widened (by 10 %
  ! of their width) if the initial value is not inside them; a positive
  ! lower bound (Tc, Pc, k, T*) stays positive: at most half the initial value.
  subroutine set_bounds(k, b, r)
    integer, intent(in) :: k
    double precision, intent(in) :: b(2), r(2)
    double precision :: w
    pMin(k) = minval(b)
    pMax(k) = maxval(b)
    w = pMax(k) - pMin(k)
    if (XX(k) <= pMin(k)) then
      if (pMin(k) > 0D0) then
        pMin(k) = max(XX(k) - 0.1D0*w, 0.5D0*XX(k))
      else
        pMin(k) = XX(k) - 0.1D0*w
      endif
    endif
    if (XX(k) >= pMax(k)) pMax(k) = XX(k) + 0.1D0*w
    if (pMin(k) /= minval(b) .or. pMax(k) /= maxval(b)) write(*,'(A,I0,A,G12.5,A,2G12.5)') &
      ' WARNING: initial value of parameter ', k, ' (', XX(k), ') out of its bounds; they become', pMin(k), pMax(k)
    pLo(k) = minval(r)
    pHi(k) = maxval(r)
  endsubroutine set_bounds

endsubroutine rkpr_fit_bounds

!--------------------------------------------------------------------------------------------------------------------------
!
!  Sets the parameters of the system from the fitted ones (unscaled, in
!  the order of read_rkpr_fit), and returns Tc and Pc (atm).
!
subroutine rkpr_set_params(XXAUX, Tc_out, Pc_out)
  double precision, intent(in) :: XXAUX(*)
  double precision, intent(out) :: Tc_out(*), Pc_out(*)
  integer :: i, j, m, n, nId, k
  logical :: crit(NCM)

  n = 0
  crit = .false.
  do m = 1, nTcF
    n = n + 1
    sysTc(idTcF(m)) = XXAUX(n)
    crit(idTcF(m)) = .true.
  enddo
  do m = 1, nPcF
    n = n + 1
    sysPc(idPcF(m)) = XXAUX(n)/atm_to_bar
    crit(idPcF(m)) = .true.
  enddo
  do m = 1, nDelF
    n = n + 1
    sysDel1(idDelF(m)) = XXAUX(n)
    crit(idDelF(m)) = .true.
  enddo
  do m = 1, nKF
    n = n + 1
    sysK(idKF(m)) = XXAUX(n)
    crit(idKF(m)) = .true.
  enddo
  do i = 1, nSys
    if (crit(i)) call critical_to_params(i)
  enddo

  nId = 0
  do m = 1, nKijF
    i = ijK(m,1)
    j = ijK(m,2)
    n = n + 1
    call set_pair(sysKinf, i, j, XXAUX(n))
    if (iTempK(m) >= 1) then
      n = n + 1
      call set_pair(sysKp, i, j, XXAUX(n))
    endif
    if (iTempK(m) == 2) then
      n = n + 1
      call set_pair(sysTstar, i, j, XXAUX(n))
    endif
    do k = nId + 1, nId + nIdK(m)
      call set_pair(sysKinf, idIdK(k,1), idIdK(k,2), sysKinf(i,j))
      call set_pair(sysKp, idIdK(k,1), idIdK(k,2), sysKp(i,j))
      call set_pair(sysTstar, idIdK(k,1), idIdK(k,2), sysTstar(i,j))
    enddo
    nId = nId + nIdK(m)
  enddo
  nId = 0
  do m = 1, nLijF
    i = ijL(m,1)
    j = ijL(m,2)
    n = n + 1
    call set_pair(sysL, i, j, XXAUX(n))
    do k = nId + 1, nId + nIdL(m)
      call set_pair(sysL, idIdL(k,1), idIdL(k,2), sysL(i,j))
    enddo
    nId = nId + nIdL(m)
  enddo

  ! b_ij, as readcomp
  if (any(crit) .or. nLijF > 0) then
    do i = 1, nSys
      do j = 1, nSys
        if (sysNComb == 0) then
          sysBij(i,j) = (1D0 - sysL(i,j))*(sysB(i) + sysB(j))/2D0
        else
          sysBij(i,j) = ((1D0 - sysL(i,j))*(sysB(i)**(1D0/3D0) + sysB(j)**(1D0/3D0))/2D0)**3
        endif
      enddo
    enddo
  endif
  Tc_out(:nSys) = sysTc(:nSys)
  Pc_out(:nSys) = sysPc(:nSys)

endsubroutine rkpr_set_params

!--------------------------------------------------------------------------------------------------------------------------
!
!  Symmetric binary parameter: a(i1,i2) = a(i2,i1) = v.
!
subroutine set_pair(a, i1, i2, v)
  double precision, intent(inout) :: a(NCM,NCM)
  integer, intent(in) :: i1, i2
  double precision, intent(in) :: v
  a(i1,i2) = v
  a(i2,i1) = v
endsubroutine set_pair

!--------------------------------------------------------------------------------------------------------------------------
!
!  ac and b of component i from its Tc, Pc and delta1, and omega from k
!  (for the output), as readcomp does (Pc in atm, cm3: R = Rgas2).
!
subroutine critical_to_params(i)
  integer, intent(in) :: i
  double precision :: d1, y, OmegaA, OmegaB, Zc, AA, BB, CC, disc

  d1 = sysDel1(i) + 1D0
  y = 1D0 + (2D0*d1)**(1D0/3D0) + (4D0/d1)**(1D0/3D0)
  d1 = (1D0 + sysDel1(i)**2)/(1D0 + sysDel1(i))
  OmegaB = 1D0/(3D0*y + d1 - 1D0)
  OmegaA = (3D0*y*y + 3D0*y*d1 + d1*d1 + d1 - 1D0)*OmegaB**2
  sysAc(i) = OmegaA*(Rgas2*sysTc(i))**2/sysPc(i)
  sysB(i) = OmegaB*Rgas2*sysTc(i)/sysPc(i)
  Zc = y*OmegaB
  AA = Apar(2)*Zc + Apar(1)
  BB = Bpar(2)*Zc + Bpar(1)
  CC = Cpar(2)*Zc + Cpar(1) - sysK(i)
  disc = BB*BB - 4D0*AA*CC
  if (disc >= 0D0) sysOm(i) = (-BB + sqrt(disc))/2D0/AA
endsubroutine critical_to_params

!--------------------------------------------------------------------------------------------------------------------------
!
!  Critical point of component i (number in the system) in the model: Tc (K)
!  and Vc = Zc R Tc/Pc = y b (cm3/mol), with the current parameters.
!
subroutine rkpr_critical(i, Tc, Vc)
  integer, intent(in) :: i
  double precision, intent(out) :: Tc, Vc
  double precision :: d1, y

  d1 = sysDel1(i) + 1D0
  y = 1D0 + (2D0*d1)**(1D0/3D0) + (4D0/d1)**(1D0/3D0)
  Tc = sysTc(i)
  Vc = y*sysB(i)
endsubroutine rkpr_critical

!--------------------------------------------------------------------------------------------------------------------------
!
!  Writes to unit u the fitted parameters (what read_rkpr_fit read), one
!  per line with its number.
!
subroutine write_rkpr_fit_spec(u)
  integer, intent(in) :: u
  integer :: m, n, k, nId

  write(u,'(/,X,A)') 'Parameters to estimate (components and pairs by their number in the system):'
  n = 0
  do m = 1, nTcF
    n = n + 1
    write(u,'(X,I4,2X,"Tc/K",6X,A)') n, sysName(idTcF(m))
  enddo
  do m = 1, nPcF
    n = n + 1
    write(u,'(X,I4,2X,"Pc/bar",4X,A)') n, sysName(idPcF(m))
  enddo
  do m = 1, nDelF
    n = n + 1
    write(u,'(X,I4,2X,"delta1",4X,A)') n, sysName(idDelF(m))
  enddo
  do m = 1, nKF
    n = n + 1
    write(u,'(X,I4,2X,"k",9X,A)') n, sysName(idKF(m))
  enddo
  nId = 0
  do m = 1, nKijF
    n = n + 1
    write(u,'(X,I4,2X,"k_inf",5X,I0,"-",I0,:,"   identical:",*(2X,I0,"-",I0))') n, ijK(m,:), &
      (idIdK(k,:), k = nId + 1, nId + nIdK(m))
    if (iTempK(m) >= 1) then
      n = n + 1
      write(u,'(X,I4,2X,"k''",8X,I0,"-",I0)') n, ijK(m,:)
    endif
    if (iTempK(m) == 2) then
      n = n + 1
      write(u,'(X,I4,2X,"T*/K",6X,I0,"-",I0)') n, ijK(m,:)
    endif
    nId = nId + nIdK(m)
  enddo
  nId = 0
  do m = 1, nLijF
    n = n + 1
    write(u,'(X,I4,2X,"l",9X,I0,"-",I0,:,"   identical:",*(2X,I0,"-",I0))') n, ijL(m,:), &
      (idIdL(k,:), k = nId + 1, nId + nIdL(m))
    nId = nId + nIdL(m)
  enddo
  if (n == 0) write(u,'(X,A)') '  none'
endsubroutine write_rkpr_fit_spec

endmodule rkpr_bridge
