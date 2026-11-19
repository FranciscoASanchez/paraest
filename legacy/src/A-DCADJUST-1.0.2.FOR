	SUBROUTINE ASSOCDCADJUSTSUB(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,TEMP,PSAT,NMAX,PSATCALC)
C
C  MODELO ASSOC-GC-EOS
C
C no contiene modulo de inicializacion de dc. Se supone que dc viene
C  inicializado desde el calling program
C
C DC(1):
C
C on input dc contains the initial value of dc
C on output dc contains the final value of dc
C
C	' PSAT, ATM ='
C
C ESTE PROGRAMA PERMITE AJUSTAR EL DIAMETRO DE ESFERAS RIGIDAS DE UN
C DETERMINADO COMPONENTE DE MANERA DE REPRODUCIR UN PUNTO DE LA CURVA
C  DE PRESION DE VAPOR DEL COMPONENTE PURO.
C

C Las coordenadas del punto de la curva de PVAP a reproducir
C     son  (TEMP,PSAT)   en KELVIN y ATM RESPECTIVAMENTE
C
C    PSATCALC = Presion de Vapor calculada (en ATM) calculada con
C               el valor final de DC. Se supone que es muy parecida
C	        a PSAT. (a la temperatura TEMP)
C
C USA NEWTON para la primer
C  iteracion. EN TODAS LA ITERACIONES SUBSIGUIENTES USA UNA PARABOLA
C
C  CONTIENE UNA OPCION QUE PERMITE LIMITAR EL NUMERO DE ITERACIONES
C
C  LLAMA A UNA SUBRUTINA QUE CALCULA PRESIONES DE VAPOR SIN REQUERIR
C   INICIALIZACION EN LA PRESION DE VAPOR
C
C  NMAX = NUMERO MAXIMO DE ITERACIONES PERMITIDO
C        SI NMAX = 1 SOLO SE HACE UN PASO DE NEWTON
C        SI NMAX = 2 O MAYOR ,   SE HACE UN PASO DE NEWTON
C                     Y TODOS LOS SIGUIENTES USAN UNA PARABOLA, EN LUGAR
C                      DE LA RECTA TANGENTE DE NEWTON
C        HASTA HOY (18-AGO-1993) LA EXPERIENCIA MUESTRA QUE CON NMAX = 2
C          LA PRECISION ES BUENA.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	integer, parameter :: NGM = 10
	
	DIMENSION DC(1), D(1), RN(1)
	DIMENSION  TEMPCRIT(1)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
C DEL TERMINO DE DISPERSION
C
	DIMENSION NU(NGM,1), QAREA(NGM), TSTR(NGM), GSTR(NGM)
	DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
	DIMENSION ALFAIJ(NGM,NGM)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
C
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,1)
C
C
C	TEMP= TEMPERATURA DEL SISTEMA EN K
C
C       D(I)= Diametro de esferas rigidas (funcion de la
C                                          temperatura) (un valor por
C                                                              componente)
C
C       DC(I)= Diametro de esferas rigidas en el punto critico
C                                          (un valor por componente)
C
C           D(I) es el diametro de una esfera imaginaria de un volumen
C                igual al de un mol de esferas rigidas. Por eso se da por
C                                                                  mol
C
C       NC  = Numero de componentes (mAX 1)
C
C       RN (I) = NUMERO DE MOLES DEL COMPONENTE I
C
C
C
C
C
C
C CALCULO DE DC, QUE AJUSTA EL PUNTO DE SATURACION
C
	CALL DCCALC(TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,DC,TEMP,PSAT,PSATCALC,NMAX)
C
C	CALL LIB$SHOW_TIMER(,2,,)
C
	ENDsubroutine
C/////////////////////////////////////////////////////////////////////////
C////////////////////////////////////////////////////////////////////
C////////////////////////////////////////////////////////////////////
	SUBROUTINE DCCALC(TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,DC,TEMP,PSAT,PSATCALC,NMAX)
C
C CALCULO DE DC
C
	IMPLICIT REAL*8 (A-H,O-Z)
	integer, parameter :: NGM = 10
	
	DIMENSION DC(1), D(1), RN(1)
	DIMENSION  TEMPCRIT(1)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
C DEL TERMINO DE DISPERSION
C
	DIMENSION NU(NGM,1), QAREA(NGM), TSTR(NGM), GSTR(NGM)
	DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
	DIMENSION ALFAIJ(NGM,NGM)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
C
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,1)
C
C ARREGLOS AUXILIARES
C
	DIMENSION X(2),Y(2)
C
	NITERP=0
	NINFA=0
	MINFA=0
C
C
C GUARDO EL VALOR DE PSAT Y DE DC
C
	PSAT0 = PSAT
C
	DC0 = DC(1)
C
C
 	DC1 = DC0
C
 321	DC(1) = DC1
C
	CALL ASSOCPSATSUB(TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,DC,TEMP,PSAT,ZFASL,ZFASV,IERRPSAT)
C
	PSAT1 = PSAT
	IF(NITERP.GE.NMAX) GO TO 999
	MINFA = MINFA + 1
	NINFA = NINFA + 1
	X(NINFA) = PSAT1
	Y(NINFA) = DC1
C
	DC2 = DC1 + 1.D-5*DC1
C
	DC(1) = DC2
C
	CALL ASSOCPSATSUB(TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,DC,TEMP,PSAT,ZFASL,ZFASV,IERRPSAT)
C
	PSAT2 = PSAT
C
	if (psat2 == psat1) then
	  psat2 = psat1*dc2/dc1
	endif
	YPRIMA2 = (DC2-DC1)/(PSAT2-PSAT1)
C
	IF(MINFA.EQ.1) THEN
C                       Newton
!
!	  Paso:
	  dDc =  YPRIMA2*(PSAT0 - PSAT1)
	  if(dAbs(dDc).GT.0.1)dDc = 0.1*(dDc/dAbs(dDc))
	  dc3 = dc1 + dDc
!	  DC3 = YPRIMA2*(PSAT0 - PSAT1) + DC1
	  NITERP= NITERP+1

	ELSE
!	  Acelero con parabola
	  CALL PARAB(X,Y,YPRIMA2,AA,BB,CC)
	  DC3 = AA*PSAT0**2 + BB*PSAT0 + CC
!
!	  Hay un paso máximo:
        if(dAbs(dc3 - dc1).GT.0.15D0)dc3 = 1.15D0*dc1
!
	  NINFA = 0
	  X(2) = PSAT1
	  Y(2) = DC1
	  NITERP=NITERP+1

	END IF
C
	IF(DABS(DC3-DC1).LT.1.D-8.AND.DABS(PSAT1-PSAT0).LT.1.D-8) THEN
					GO TO 999
				   ELSE
					DC1 = DC3
					GO TO 321
	END IF
C
C
C
C
C
C
 999	PSATCALC = PSAT1
	PSAT = PSAT0
 	ENDSUBROUTINE DCCALC
C////////////////////////////////////////////////////////
	SUBROUTINE PARAB(X,Y,YPRIMA2,A,B,C)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION X(2),Y(2)
C
C CALCULA LOS COEFICIENTES DE UNA PARABOLA PASANTE POR 2 PUNTOS
C   CUYA DERIVADA EN EL PUNTO X(2),Y(2) COINCIDE CON EL VALOR
C     ALIMENTADO EN YPRIMA2
C
C	        2
C	 Y = A X  +  B X  +  C
C
	A21 = Y(2) - Y(1)
	B21 = X(2) - X(1)
	C21 = X(2)**2 - X(1)**2
C
	A = (A21 - YPRIMA2 * B21)/ (C21 - 2 * B21 * X(2))
C
	B = YPRIMA2 - 2 * A * X(2)
C
	C = Y(2) - A * X(2)**2 - B * X(2)
C
C
	ENDSUBROUTINE PARAB
C//////////////////////////////////////////////////////////////
C/////////////////////////////////////////////////////////////////

