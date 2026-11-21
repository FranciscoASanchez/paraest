!--------------------------------------------------------------------------------------------------------------------------
!
!  Módulo parametros_gca: estructura y parámetros del modelo GCA-EoS, y la
!  correspondencia entre los parámetros a ajustar y su lugar en el modelo.
!
!  Módulo de transición: cuando paraest use el paquete GCA (con sus
!  parámetros en módulos), esto se reemplaza por el módulo del paquete.
!
!--------------------------------------------------------------------------------------------------------------------------
module parametros_gca
  use dimensiones, only: NCM, NGM, NPTIPOM, NIDM
  implicit none

  ! --- Estructura ---------------------------------------------------------
  integer :: NG, NGA                       ! grupos y grupos asociativos
  integer :: NU(NGM,NCM)                   ! grupos de cada componente
  integer :: MASSOC(12)                    ! sitios de cada grupo asociativo
  double precision :: ONUASOC(12,NCM)      ! grupos asociativos de cada componente

  ! --- Parámetros dispersivos (por grupo y por par de grupos) -------------
  double precision :: QAREA(NGM), TSTR(NGM), GSTR(NGM), GDOT(NGM), GDDOT(NGM), TSPLINE(NGM)
  double precision :: RKIJ(NGM,NGM), RKTIJ(NGM,NGM), ALFAIJ(NGM,NGM)

  ! --- Parámetros asociativos ---------------------------------------------
  double precision :: ENASS(6,12,6,12), RKASS(6,12,6,12)

  ! --- Diámetro crítico ---------------------------------------------------
  double precision :: DC(NCM)
  ! Punto de presión de vapor al que se ajusta dc
  double precision :: TSAT(NCM), PSAT(NCM)

  ! --- Parámetros a ajustar -----------------------------------------------
  ! Cantidad de cada tipo
  integer :: NNP, NNQ, NALFA, NENASS, NRKASS, numTc, numDC
  ! g_ii / k_ij: índices, cuántos idénticos tiene cada uno y cuáles son
  integer :: IDG(NPTIPOM,2), ITEMP(NPTIPOM), NIDG(NPTIPOM), IDIDG(NIDM,2)
  integer :: IDQ(NPTIPOM)                                                     ! áreas
  integer :: IDALF(NPTIPOM,2), NIDALF(NPTIPOM), IDIDALF(NIDM,2)               ! alfa
  integer :: IDENASS(NPTIPOM,4), NIDENASS(NPTIPOM), IDIDENASS(NIDM,4)         ! energía asoc.
  integer :: IDRKASS(NPTIPOM,4), NIDRKASS(NPTIPOM), IDIDRKASS(NIDM,4)         ! volumen asoc.
  integer :: idTc(NPTIPOM), IDDC(NPTIPOM)                                     ! Tc y dc

endmodule parametros_gca
