!--------------------------------------------------------------------------------------------------------------------------
!
!  Module gca_parameters: structure and parameters of the GCA-EoS model, and
!  the mapping between the fitted parameters and their place in the model.
!
!  Transitional module: when paraest uses the GCA package (with its
!  parameters in modules), this is replaced by the package's module.
!
!--------------------------------------------------------------------------------------------------------------------------
module gca_parameters
  use dimensions, only: NCM, NGM, NPTIPOM, NIDM
  implicit none

  ! --- Structure ----------------------------------------------------------
  integer :: NG, NGA                       ! groups and associating groups
  integer :: NU(NGM,NCM)                   ! groups of each component
  integer :: MASSOC(12)                    ! sites of each associating group
  double precision :: ONUASOC(12,NCM)      ! associating groups of each component

  ! --- Dispersive parameters (per group and per pair of groups) -------------
  double precision :: QAREA(NGM), TSTR(NGM), GSTR(NGM), GDOT(NGM), GDDOT(NGM), TSPLINE(NGM)
  double precision :: RKIJ(NGM,NGM), RKTIJ(NGM,NGM), ALFAIJ(NGM,NGM)

  ! --- Association parameters ---------------------------------------------
  double precision :: ENASS(6,12,6,12), RKASS(6,12,6,12)

  ! --- Critical diameter --------------------------------------------------
  double precision :: DC(NCM)
  ! Vapor pressure point to which dc is fitted
  double precision :: TSAT(NCM), PSAT(NCM)

  ! --- Fitted parameters --------------------------------------------------
  ! Number of each type
  integer :: NNP, NNQ, NALFA, NENASS, NRKASS, numTc, numDC
  ! g_ii / k_ij: indices, how many identical ones each has and which they are
  integer :: IDG(NPTIPOM,2), ITEMP(NPTIPOM), NIDG(NPTIPOM), IDIDG(NIDM,2)
  integer :: IDQ(NPTIPOM)                                                     ! areas
  integer :: IDALF(NPTIPOM,2), NIDALF(NPTIPOM), IDIDALF(NIDM,2)               ! alpha
  integer :: IDENASS(NPTIPOM,4), NIDENASS(NPTIPOM), IDIDENASS(NIDM,4)         ! assoc. energy
  integer :: IDRKASS(NPTIPOM,4), NIDRKASS(NPTIPOM), IDIDRKASS(NIDM,4)         ! assoc. volume
  integer :: idTc(NPTIPOM), IDDC(NPTIPOM)                                     ! Tc and dc

endmodule gca_parameters
