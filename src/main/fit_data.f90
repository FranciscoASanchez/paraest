!--------------------------------------------------------------------------------------------------------------------------
!
!  Module fit_data: the data the main program reads and FUNC needs at each
!  evaluation. ZXSSQ calls FUNC only with (XX, NCXP, NPAR, FCION), so the
!  rest reaches it through here.
!
!  This holds what does not depend on the model; what belongs to the GCA is
!  in module gca_parameters.
!
!--------------------------------------------------------------------------------------------------------------------------
module fit_data
  use dimensions, only: NCM, NDM, NSetM, NPARM
  implicit none

  ! --- Pure components ----------------------------------------------------
  integer :: NC                                     ! number of components
  double precision :: PRCRIT(NCM), TEMPCRIT(NCM)
  integer :: NSAT(NCM), IDSAT(NCM), IOBJSAT(NCM)
  ! Melting and solid-solid transition data (solid-liquid equilibrium):
  ! DCpfus = DCpfu + BCpfu*T; NTSS transitions (Ttr, DHtr), from the
  ! highest temperature down; TtrSS, DHss: the first one.
  integer, parameter :: NTRM = 10
  double precision :: DHfu(NCM), Tfu(NCM), DCpfu(NCM), BCpfu(NCM), DVfu(NCM), Pfu(NCM), TtrSS(NCM), DHss(NCM)
  integer :: NTSS(NCM)
  double precision :: Ttr(NTRM,NCM), DHtr(NTRM,NCM)

  ! --- Weights of the objective function ----------------------------------
  double precision :: W1, W2, W3, W4, W5, W6, W7, W8, W9, W10, W11, W12, W13
  double precision :: wSat(NDM), wBin(2,NDM), wTer(2,NDM), wEx(NDM), wGam(3,NDM)
  integer :: ISAT, IGLE, IVE, IUNI

  ! --- Data sets ----------------------------------------------------------
  integer :: IDC(NSetM,3), NOBS(NSetM), IHEN(NSetM), mez(NSetM), IFL(NSetM)
  double precision :: STDVA(NSetM), STDVB(NSetM), STDVC(NSetM)

  ! --- Experimental data --------------------------------------------------
  double precision :: TEXP(NDM), PEXP(NDM), XEXP(NDM), YEXP(NDM), xExp2(NDM), yExp2(NDM)
  double precision :: vSatExp(NDM), rGamExp1(NDM), rGamExp2(NDM), rGamExp3(NDM), prop_E_exp(NDM)

  ! --- Calculated values --------------------------------------------------
  double precision :: PCALC(NDM), DHVL(NDM), XCALC(NDM), YCALC(NDM), xCalc2(NDM), yCalc2(NDM)
  double precision :: ZLCALC(NDM), ZVCALC(NDM), zLiCalc(NDM)
  double precision :: rGamCalc1(NDM), rGamCalc2(NDM), rGamCalc3(NDM), prop_E_calc(NDM)
  double precision :: FCAUX(NDM), FCAUX2(NDM), FCAUX3(NDM)

  ! --- Type of each fitted parameter (see module scaling) ------------------
  integer :: IDXX(NPARM)

  ! IOBJSAT = 81 (vSatExp above): molar volume (cm3/mol; of the solid below
  ! Tfus), its weight and the calculated one
  double precision :: wVSat(NDM), vSatCalc(NDM)
endmodule fit_data
