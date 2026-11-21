!--------------------------------------------------------------------------------------------------------------------------
!
!  Módulo datos_ajuste: los datos que el programa principal lee y que FUNC
!  necesita en cada evaluación. ZXSSQ llama a FUNC solo con (XX, NCXP, NPAR,
!  FCION), así que el resto le llega por acá.
!
!  Aquí va lo que no depende del modelo; lo propio del GCA está en el módulo
!  parametros_gca.
!
!--------------------------------------------------------------------------------------------------------------------------
module datos_ajuste
  use dimensiones, only: NCM, NDM, NSetM, NPARM
  implicit none

  ! --- Componentes puros --------------------------------------------------
  integer :: NC                                     ! número de componentes
  double precision :: PRCRIT(NCM), TEMPCRIT(NCM)
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

  ! --- Tipo de cada parámetro a ajustar (ver módulo escalado) -------------
  integer :: IDXX(NPARM)

endmodule datos_ajuste
