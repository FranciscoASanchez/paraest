!--------------------------------------------------------------------------------------------------------------------------
!
!  Módulo datos_ajuste: los datos que el programa principal lee y que FUNC
!  necesita en cada evaluación. ZXSSQ llama a FUNC solo con (XX, NCXP, NPAR,
!  FCION), así que el resto le llega por acá.
!
!  Reemplaza a los COMMON /FORFUNC/ ... /forFunc6/ e /IDSC/.
!
!--------------------------------------------------------------------------------------------------------------------------
module datos_ajuste
  use dimensiones, only: NCM, NGM, NDM, NSetM, NPTIPOM, NIDM, NPARM
  implicit none

  ! --- Sistema ------------------------------------------------------------
  integer :: NC, NG, NGA                   ! componentes, grupos, grupos asociativos
  integer :: NU(NGM,NCM)                   ! grupos de cada componente

  ! --- Parámetros dispersivos (por grupo y por par de grupos) -------------
  double precision :: QAREA(NGM), TSTR(NGM), GSTR(NGM), GDOT(NGM), GDDOT(NGM), TSPLINE(NGM)
  double precision :: RKIJ(NGM,NGM), RKTIJ(NGM,NGM), ALFAIJ(NGM,NGM)

  ! --- Parámetros asociativos ---------------------------------------------
  integer :: MASSOC(12)
  double precision :: ENASS(6,12,6,12), RKASS(6,12,6,12)
  double precision :: ONUASOC(12,NCM)

  ! --- Componentes puros --------------------------------------------------
  double precision :: DC(NCM), PRCRIT(NCM), TEMPCRIT(NCM)
  double precision :: TSAT(NCM), PSAT(NCM)          ! punto de presión de vapor para ajustar dc
  double precision :: PSTCLC(NCM)                   ! presión de vapor calculada a TSAT
  integer :: NSAT(NCM), IDSAT(NCM), IOBJSAT(NCM)
  ! Datos de fusión y transición sólido-sólido (equilibrio sólido-líquido)
  double precision :: DHfu(NCM), Tfu(NCM), DCpfu(NCM), DVfu(NCM), Pfu(NCM), TtrSS(NCM), DHss(NCM)

  ! --- Pesos de la función objetivo ---------------------------------------
  double precision :: W1, W2, W3, W4, W5, W6, W7, W8, W9, W10, W11, W12, W13
  double precision :: wSat(NDM), wBin(2,NDM), wTer(2,NDM), wEx(NDM), wGam(3,NDM)
  integer :: ISAT, IGLE, IVE, IUNI

  ! --- Conjuntos de datos -------------------------------------------------
  integer :: IDC(NSetM,3), NOBS(NSetM), IHEN(NSetM), mez(NSetM), IFL(NSetM)
  double precision :: STDVA(NSetM), STDVB(NSetM), STDVC(NSetM)

  ! --- Datos experimentales -----------------------------------------------
  double precision :: TEXP(NDM), PEXP(NDM), XEXP(NDM), YEXP(NDM), xExp2(NDM), yExp2(NDM)
  double precision :: rhoExp(NDM), rGamExp1(NDM), rGamExp2(NDM), rGamExp3(NDM), prop_E_exp(NDM)

  ! --- Valores calculados -------------------------------------------------
  double precision :: PCALC(NDM), DHVL(NDM), XCALC(NDM), YCALC(NDM), xCalc2(NDM), yCalc2(NDM)
  double precision :: ZLCALC(NDM), ZVCALC(NDM), zLiCalc(NDM)
  double precision :: rGamCalc1(NDM), rGamCalc2(NDM), rGamCalc3(NDM), prop_E_calc(NDM)
  double precision :: FCAUX(NDM), FCAUX2(NDM), FCAUX3(NDM)

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
  integer :: IDXX(NPARM)                   ! tipo de cada parámetro (ver módulo escalado)

endmodule datos_ajuste
