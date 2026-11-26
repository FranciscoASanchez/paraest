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
!--------------------------------------------------------------------------------------------------------------------------
module rkpr_bridge
  implicit none
  private
  public :: load_from_readcomp, rkpr_select, write_rkpr_params

  integer, parameter :: NCM = 30                     ! must match the RK-PR package
  double precision, parameter :: Rgas = 0.08314472D0, Rgas2 = 82.05D0

  ! Whole system, as readcomp left it
  integer :: nSys = 0, sysNComb = 0, sysNTDEP = 0
  double precision :: sysAc(NCM), sysB(NCM), sysDel1(NCM), sysK(NCM)
  double precision :: sysKinf(NCM,NCM), sysKp(NCM,NCM), sysTstar(NCM,NCM), sysBij(NCM,NCM)
  double precision :: sysTc(NCM), sysPc(NCM), sysOm(NCM), sysDc(NCM)
  double precision :: sysPM(NCM), sysPen(NCM), sysH(7,NCM)
  character(10) :: sysName(NCM)

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
  integer :: NTDEP, ncomb
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
  double precision :: lij

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
      if (sysNComb == 0) then
        lij = 1D0 - 2D0*sysBij(i,j)/(sysB(i) + sysB(j))
      else
        lij = 1D0 - 2D0*sysBij(i,j)**(1D0/3D0)/(sysB(i)**(1D0/3D0) + sysB(j)**(1D0/3D0))
      endif
      write(u,'(X,I3,X,I3,4(2X,G11.4))') i, j, sysKinf(i,j), sysKp(i,j), sysTstar(i,j), lij
    enddo
  enddo
endsubroutine write_rkpr_params

endmodule rkpr_bridge
