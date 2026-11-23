!--------------------------------------------------------------------------------------------------------------------------
!
!  Module gca_bridge: temporary bridge between paraest and the GCA package.
!
!  fakePARGR and fakePARMOL replace PARGR and PARMOL of the package: they
!  put the parameters of paraest in the package's COMMON blocks instead of
!  reading them from files. They disappear when paraest uses GCA 2.0, whose
!  parameters are in modules.
!
!  negative_cohesion tells whether the attractive energy g(T) of a group of
!  the loaded components is not positive at T (the calculation is then
!  meaningless: paraest counts the point as not calculated).
!
!--------------------------------------------------------------------------------------------------------------------------
module gca_bridge
  implicit none
  private
  public :: fakePARGR, fakePARMOL, negative_cohesion, load_from_pargr

contains

!--------------------------------------------------------------------------------------------------------------------------
!
  ! Replaces PARGR of the GCA package: instead of reading the group
  ! parameters from a file, it takes those of paraest and puts them in the
  ! package's COMMON blocks (a dummy argument cannot be in a COMMON, so they
  ! are copied).
!
subroutine fakePARGR(NC,NG,NGA,NST,NUWW,QAREA,TS,GS,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,ONUASOCWW, &
&                 MASSOC,ENASS,RKASS)
  implicit none
  integer, parameter :: NCM=30, NGM=30, NGAM=14, NSM=24     ! must match the GCA package (v1.14)
  integer :: i, J, k, l, m1, m2, NC, NG, NGA, NGROUP, NownII, NownIJ, NST
  integer :: NY, major
  double precision :: A, AKIJ, ALFA, DADT, dDeldT, Delta, eps_R, EPX, G1, G2, GSTR
  double precision :: Q, R, TSPL, TSTR, XKIJ, ZZ
  character*10 AName

  integer :: nuWW(NG,NC)
  double precision :: qArea(NG), Ts(NG), Gs(NG), GDOT(NG), GDDOT(NG), rKij(NG,NG), rKTij(NG,NG), alfaij(NG,NG)
!
  ! Association energy and volume with 4 indices: only in this subroutine.
  integer :: MAssoc(NGA)
  double precision :: enAss(2,NGA,2,NGA), rkAss(2,NGA,2,NGA), onuAsocWW(NGA,NC)
  integer :: sigma
  double precision :: kappa
!
  COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),XKIJ(NGM,NGM),      &
&                AKIJ(NGM,NGM),EPX(NGM)
  COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)

  common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
  common/GrupAs2/sigma(NSM,NCM), major
  common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
  common/NAME/AName(NCM)
  COMMON/COORD/ZZ

  if(NC.GT.NCM)then

    write(*,*)'Error: maximum number of components exceeded in GC3'
    stop

  elseif(NG.GT.NGM)then

    write(*,*)'Error: maximum number of attractive groups exceeded in GC3'
    stop

  elseif(NGA.GT.NGAM)then

    write(*,*)'Error: maximum number of associating groups exceeded in GC3'
    stop

  endif
!
  ! Coordination number
  ZZ=10.D0
!
  ! Gas constant
  R=82.05D0  !atm cm3/(mol K)
  NGROUP=35
  NownII = NG
  NownIJ = ((NG + 1)**2 - NG -1)/2 - NG
!
  ! Dispersive parameters
  do j = 1,NG
!
    q(j) = qArea(j)
    Tstr(j) = Ts(j)
    ny(:NC,j) = nuWW(j,:NC)
    gStr(j) = Gs(j)
    g1(j) = GDOT(j)
    g2(j) = GDDOT(j)
    call spline_temperature(Tstr(j), g1(j), g2(j), Tspl(j), epx(j))
    do k = 1,NG

      xKij(j,k) = rKij(j,k)
      aKij(j,k) = rKTij(j,k)
      alfa(j,k) = alfaij(j,k)

    enddo

  enddo

!-----------------------------------------------------------------------
  ! Association parameters
!-------Conversion to the notation by site "m" (instead of by site "k" of
    ! group "i")
  NST = 0
  if(NGA.GT.0)then

    m1 = 0
    m2 = 0
    do i = 1,NGA
      do k = 1,Massoc(i)

        m2 = m1
        m1 = m1 + 1        ! current site "k" is site "m1"
        do j = i,NGA   ! this is valid only because, so far, there has been only
          if(i.EQ.j)then  ! one type of site in each group
            do l = k,Massoc(j)

              m2 = m2 + 1   ! current site "l" is site "m2"
              eps_R(m1,m2) = enAss(l,j,k,i)
              eps_R(m2,m1) = eps_R(m1,m2)
              kappa(m1,m2) = rkAss(l,j,k,i)
              kappa(m2,m1) = kappa(m1,m2)

            enddo
          else
            do l = 1,Massoc(j)

              m2 = m2 + 1   ! current site "l" is site "m2"
              eps_R(m1,m2) = enAss(l,j,k,i)
              eps_R(m2,m1) = eps_R(m1,m2)
              kappa(m1,m2) = rkAss(l,j,k,i)
              kappa(m2,m1) = kappa(m1,m2)

            enddo
          endif
        enddo

      enddo
    enddo
    NST = sum(Massoc(:NGA))
!
    ! Sites of each component, in the notation 1..NST
    do i = 1,NC

      m1 = 0
      do j = 1,NGA
        do k = 1, MAssoc(j)

          m1 = m1 + 1
          sigma(m1,i) = onuAsocWW(j,i)

        enddo
      enddo

    enddo
!
  endif

  return
endsubroutine fakePARGR
!---------------------------------------------------------------------------------------------------------
!
  ! Same as fakePARGR, for the component parameters of PARMOL (Tc, Pc, dc).
!
subroutine fakePARMOL(NC,NG,NST,WWDC,WWTEMPCRIT,wwPRCRIT)
!
  implicit none
  integer, parameter :: NCM=30, NGM=30, NSM=24     ! must match the GCA package (v1.14)
  integer :: i, iEntalp, NC, NG, NST, NY, major
  double precision :: HHG, omega
  double precision :: A, AKIJ, ALFA, D, DADT, DC, dDeldT, Delta, DT, eps_R, EPX, G1
  double precision :: G2, GSTR, HA, HB, HHA, HHB, hhc, HHD, HHE, HHF, PC, Pen
  double precision :: pmm, Q, R, TC, TSPL, TSTR, XKIJ, ZZ
  character*10 AName
  character*20 :: versGCA
!
  double precision :: wwDc(NC), wwTempCrit(NC), wwPRCRIT(NC)
!
  common/versSUB/versGCA
!
  COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
  COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
  common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
  common/GrupAs2/sigma(NSM,NCM), major
  common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
  COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
  COMMON/CRIT/TC(NCM),PC(NCM),omega(NCM)
!
  common/NAME/AName(NCM)
  common/GCPROM/pmm(ncm),Pen(NCM),HHA(ncm),HHB(ncm),hhc(ncm),HHD(ncm),HHE(NCM),HHF(NCM),HHG(NCM)
  common/ENTALP/iEntalp
  COMMON/COORD/ZZ

  double precision :: kappa
  integer :: sigma

  iEntalp = 0
  do i = 1,NC

    Tc(i) = wwTempCrit(i)
    dc(i) = wwDC(i)
    Pc(i) = wwPRCRIT(i)
    write(AName(i),'("Comp ",I2)')i
    PMM(i) = 0.
    Pen(i) = 0.
    HHA(i) = 0.
    HHB(i) = 0.
    HHC(i) = 0.
    HHD(i) = 0.
    HHE(i) = 0.
    HHF(i) = 0.
    HHG(i) = 0.
    omega(i) = 0.

  enddo
!
  ! iEntalp: number of components with Passut-Danner constants
  iEntalp = NC
  do i = 1,NC
    if((HHA(i).EQ.0).AND.(HHB(i).EQ.0).AND.(HHC(i).EQ.0).AND.(HHD(i).EQ.0).AND.(HHE(i).EQ.0).AND.(HHF(i).EQ.0))then

      iEntalp = iEntalp - 1

    endif
  enddo
  return
endsubroutine fakePARMOL

!---------------------------------------------------------------------------------------------------------
!
!  True if g(T) = g*(1 + g'(T/T* - 1) + g" ln(T/T*)) <= 0 at temperature T
!  for a group present in any of the NC components loaded by fakePARGR.
!  Above the spline temperature PARAGC uses g*/4/(T/Tspl)**epx, which is
!  positive. A group that is not present does not count: its g(T) only
!  enters through k_ij terms that are multiplied by zero.
!
logical function negative_cohesion(T, NC, NG)
  implicit none
  integer, parameter :: NCM=30, NGM=30     ! must match the GCA package (v1.14)
  double precision, intent(in) :: T
  integer, intent(in) :: NC, NG
  integer :: j, NY
  double precision :: A, AKIJ, ALFA, DADT, EPX, G1, G2, GSTR, Q, R, TSPL, TSTR, XKIJ, Tr

  COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),XKIJ(NGM,NGM),      &
&                AKIJ(NGM,NGM),EPX(NGM)
  COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)

  negative_cohesion = .false.
  do j = 1, NG
    if (all(NY(:NC,j) == 0)) cycle
    if (T >= TSPL(j)) cycle
    Tr = T/TSTR(j)
    if (GSTR(j)*(1D0 + G1(j)*(Tr - 1D0) + G2(j)*dlog(Tr)) <= 0D0) then
      negative_cohesion = .true.
      return
    endif
  enddo
endfunction negative_cohesion

!---------------------------------------------------------------------------------------------------------
!
!  After PARGR of the GCA package has read the attractive and associative
!  input (database FOR010-FOR015 plus the user's values), copies what it
!  left in its COMMON blocks and in module pargr_export to paraest's own
!  arrays (module gca_parameters), which the fit changes and fakePARGR puts
!  back in the package at each evaluation.
!
subroutine load_from_pargr(NC)
  use gca_parameters, only: NG, NGA, NU, MASSOC, ONUASOC, QAREA, TSTR, GSTR, GDOT, GDDOT, &
                            RKIJ, RKTIJ, ALFAIJ, ENASS, RKASS, groupID, assocID, groupName, assocName
  use pargr_export
  implicit none
  integer, parameter :: NCM=30, NGM=30, NGAM=14     ! must match the GCA package (v1.14)
  integer, intent(in) :: NC
  integer :: i, j, k, l, NY
  double precision :: A, AKIJ, ALFA, DADT, EPX, G1, G2, GSTR_, Q, R, TSPL, TSTR_, XKIJ
  character(10) :: CName, GName, GAName

  COMMON/GROUP1/GSTR_(NGM),G1(NGM),G2(NGM),TSTR_(NGM),TSPL(NGM),XKIJ(NGM,NGM),      &
&                AKIJ(NGM,NGM),EPX(NGM)
  COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)
  common/NAME/CName(NCM), GName(NGM), GAName(NGAM)

  do j = 1, NG
    groupID(j) = px_groupID(j)
    groupName(j) = GName(j)
    QAREA(j) = Q(j)
    TSTR(j) = TSTR_(j)
    GSTR(j) = GSTR_(j)
    GDOT(j) = G1(j)
    GDDOT(j) = G2(j)
    do i = 1, NC
      NU(j,i) = NY(i,j)
    enddo
    do k = 1, NG
      RKIJ(j,k) = XKIJ(j,k)
      RKTIJ(j,k) = AKIJ(j,k)
      ALFAIJ(j,k) = ALFA(j,k)
    enddo
  enddo
  ! PARGR does not fill the diagonal of k' and alpha when there is only one
  ! group; nor does it matter: PARAGC uses only i < j.
  ENASS = 0D0
  RKASS = 0D0
  do i = 1, NGA
    assocID(i) = px_assocID(i)
    assocName(i) = GAName(i)
    MASSOC(i) = px_sites(i)
    if (MASSOC(i) > 2) then
      write(*,'(A,I0,A)') ' Error: associating group ', assocID(i), &
        ' has more than 2 sites (paraest supports up to 2)'
      stop 1
    endif
    do j = 1, NC
      ONUASOC(i,j) = px_nyAss(j,i)
    enddo
    do j = 1, NGA
      do k = 1, MASSOC(i)
        do l = 1, px_sites(j)
          ENASS(k,i,l,j) = px_enAs(k,i,l,j)
          RKASS(k,i,l,j) = px_volAs(k,i,l,j)
        enddo
      enddo
    enddo
  enddo
endsubroutine load_from_pargr

endmodule gca_bridge
