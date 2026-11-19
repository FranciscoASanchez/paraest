!
!
!	Calcula un flash a T y P para un sistema binario con el modelo ASS-GC-EOS
!	Utiliza metodo de Newton. Si este diverge aplica sustitucion directa
!
SUBROUTINE ASSocTPflashSUB(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS, &
&                          TEMP,PRESS,X1,Y1,ZLIQUIDO,ZVAPOR,IERRTP)

	IMPLICIT REAL*8 (A-H,O-Z)

	DIMENSION DC(2)
	DIMENSION  TEMPCRIT(2)
!
!	Dimensionamiento de arreglos relacionados con el caracter grupal del termino de dispersion:
	DIMENSION NU(10,2), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
!
!	Dimensionamiento de arreglos relacionados con el termino de asociacion
	DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
	DIMENSION ONUASOC(12,2)
!
!	Dimensionamiento de arreglos propios de un calculo flash
	DIMENSION XLIQ(2), YVAP(2)
!
!	Dimensionamiento de coeficientes de fugacidad en fase liquida (lq) y vapor (vp)
	DIMENSION FFILQ(2)
	DIMENSION FFIVP(2)
!
!	Valores de x y de la funcion que debe hacerse nula
	DIMENSION FISOFC(2)
	DIMENSION FISOF(6)
	DIMENSION XISOF(6)
!
!	FisoF(1)  = Diferencia entre las fugacidades del comp 1 en liq y vapor
!	            (si converge esto termina siendo cero), en el punto inicial (o central)
!
!	FisoF(2)  = IDEM FISOF(1) pero para el componente 2
!	FisoF(3)  = IDEM FISOF(1) pero cuando se perturba Y1 para calcular el jacobiano.
!	FisoF(4)  = IDEM FISOF(2) pero cuando se perturba Y1
!	FisoF(5)  = IDEM FISOF(1) pero cuando se perturba X1
!	FisoF(6)  = IDEM FISOF(2) pero cuando se perturba X1
!
!	XisoF(1)  = Y1 inicial (o central), fraccion molar del componente 1 en fase vapor.
!
!	XISOF(2)  = X1 inicial (o central), fraccion molar del componente 1 en fase líquida.
!
!	Lectura de temperatura, presion  y de los valores inciales de la composicion del vapor y de liquido.
!
	ITERRR=0
	naux=0
	IERRTP = 0
	RUNIV = 82.05D0 !ATM cm3 / (K mol)
!
	XLIQ(1)= X1
	XLIQ(2)= 1.D0 - X1
	YVAP(1)= Y1

	YVAP(2)= 1.D0 - Y1
!
!	Variable auxiliar
909	NINJA = 1
	ITERRR=ITERRR+1
!
!	Almaceno los valores centrales de y1 y x1
 	reservY1 = YVAP(1)
	reservX1 = XLIQ(1)
!
!	Inicializacion del volumen liquido
 192	CALL LIQMIXINICI(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS, &
&	                 TEMP,PRESS,XLIQ,VLIQ,IERRINI)
!
!	Si hay error fin de ejecucion
	IF(IERRINI.NE.0) WRITE(6,*) 'Error liq ini'
	IF(IERRINI.NE.0) GO TO 107
!
!	Compatibilizacion del volumen del liquido con la presion press
	CALL LIQMIXCOMPAT(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS, &
&	                  TEMP,PRESS,XLIQ,VLIQ,FFILQ,IERRLC)

!	Si hay error fin de ejecucion
	IF(IERRLC.NE.0) WRITE(6,*) 'Error liq compat'
	IF(IERRLC.NE.0) GO TO 107
!
!	Almaceno volumen liquido correspondiente al punto central
	IF(NINJA.EQ.1) RESERVVLIQ = VLIQ
!
!	Inicializacion del volumen del vapor
888	If (naux==0) then

	  CALL VAPMIXINICI(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS, &
&	                   TEMP,PRESS,YVAP,VVAP,IERRVI)
!
!	  Si hay error fin de ejecucion
	  IF(IERRVI.NE.0) WRITE(6,*) 'Error vapini'
!
!	  Compatibilizacion del volumen del vapor con la presion press
	  CALL VAPMIXCOMPAT(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS, &
     &                    TEMP,PRESS,YVAP,VVAP,FFIVP,IERRVC)

	end if
!
!	Si hay error fin de ejecucion
	IF(IERRVC.NE.0) WRITE(6,*) 'Error vap compat'
!
!	Almaceno volumen vapor correspondiente al punto central
	IF(NINJA.EQ.1) RESERVVVAP = VVAP
!
!	Calculo de la "isofugacidad" (mas bien, error en la fugacidad)
	FISOFC(1)= FFILQ(1)*XLIQ(1) - FFIVP(1)*YVAP(1)
	FISOFC(2)= FFILQ(2)*XLIQ(2) - FFIVP(2)*YVAP(2)
!
!	Almaceno coeficientes de fugacidad correspondientes al pto. central
	IF(NINJA.EQ.1) THEN

	  RESFI1L = FFILQ(1)
	  RESFI1V = FFIVP(1)
	  RESFI2L = FFILQ(2)
	  RESFI2V = FFIVP(2)

	END IF
!
!	Almacenamiento de valores
	XISOF(2*NINJA-1) = YVAP(1)
	XISOF(2*NINJA) = XLIQ(1)
!
	FISOF(2*NINJA-1) = FISOFC(1)
	FISOF(2*NINJA)   = FISOFC(2)
!
!	PERTURBO YVAP(1)
	IF(NINJA.EQ.1) THEN

	  YVAP(1) = YVAP(1)*(1.001)
	  IF(YVAP(1).GT.1.) YVAP(1) = 2.D0 - YVAP(1)
	  YVAP(2) = 1D0 - YVAP(1)
	  XLIQ(1) = RESERVX1
	  XLIQ(2) = 1D0 - XLIQ(1)
	  NINJA = 2
	  GO TO 888

	END IF
!
!	Perturbo xliq(1)
	IF(NINJA.EQ.2)  THEN

	  XLIQ(1) = XLIQ(1)*(1.001)
	  IF(XLIQ(1).GT.1.) XLIQ(1)= 2.D0 - XLIQ(1)
	  XLIQ(2) = 1D0 - XLIQ(1)
	  YVAP(1) = RESERVY1
	  YVAP(2) = 1D0 - YVAP(1)
	  NINJA = 3
	  GO TO 192

	END IF
!
!	Actualizacion de y1 y x1  aplicando met. de newton (o sust. directa)
	IF(NINJA.EQ.3) THEN
!
	  XLIQ(1) = RESERVX1
	  XLIQ(2) = 1 - XLIQ(1)
	  yVAP(1) = RESERVY1
	  YVAP(2) = 1 - YVAP(1)
	  CALL UPDATE2(XISOF,FISOF,YVAP,XLIQ)
!
!	  Si Newton no funciona ==> Sustitucion directa
	  IF(XLIQ(1).LT.0.OR.XLIQ(1).GT.1.OR.YVAP(1).LT.0.OR.YVAP(1).GT.1) THEN

	    CALL UPDATE5(RESFI1L,RESFI1V,RESFI2L,RESFI2V,RESERVY1,RESERVX1,YVAP,XLIQ)

	  END IF
!
!	  TEST DE CONVERGENCIA
	  IF (ITERRR.GT.200) THEN

	    IERRTP = 1
!----------Salto-------------
	    GOTO 107

	  END IF
!----------------------------Salto---------------------
	  IF(DABS(YVAP(1)-RESERVY1).LT.1.D-9.AND.DABS(XLIQ(1)-RESERVX1).LT.1.D-7) GO TO 191
!---------------------------otro salto-------
	  GO TO 909

	END IF
!
!	Si se cumple el test de convergencia el programa informa todo lo referente al punto central
191	XLIQ(1) = RESERVX1
	XLIQ(2) = 1D0 - XLIQ(1)
	YVAP(1) = RESERVY1
	YVAP(2) = 1D0 - YVAP(1)
	VLIQ = RESERVVLIQ
	VVAP = RESERVVVAP
	ZLIQUIDO = PRESS*VLIQ/RUNIV/TEMP
	ZVAPOR   = PRESS*VVAP/RUNIV/TEMP
	X1 = XLIQ(1)
	Y1 = YVAP(1)
!
!	Las variables de interes estan en , yvap, press, vliq y vvap, al final de la corrida
!
107	IF(IERRLC.NE.0.OR.ierrini2.ne.0.or.IERRINI.NE.0)IERRTP = 1

	return
ENDsubroutine !FIN  ASSocTPflashSUB
!
!-----------------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!	Recalcula yvap y xliq utilizando el metodo de newton
!
!	RJAC es el jacobiano (derivadas parciales por diferencias finitas)
!
SUBROUTINE UPDATE2(XISOF,FISOF,YVAP,XLIQ)

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION FISOF(6)
	DIMENSION XISOF(6)
	DIMENSION XLIQ(2), YVAP(2)
	DIMENSION RJAC(2,2)
!
	RJAC(1,1) = (FISOF(3)-FISOF(1))/ (XISOF(3)-XISOF(1))
	RJAC(1,2) = (FISOF(5)-FISOF(1))/ (XISOF(6)-XISOF(2))
	RJAC(2,1) = (FISOF(4)-FISOF(2))/ (XISOF(3)-XISOF(1))
	RJAC(2,2) = (FISOF(6)-FISOF(2))/ (XISOF(6)-XISOF(2))
!
!	Calculo del determinante
	DETERM = RJAC(1,1)*RJAC(2,2) - RJAC(1,2)*RJAC(2,1)
!
!	Recalculo e y1 y x1
	YVAPAUX = FISOF(2) * RJAC(1,2) - FISOF(1) * RJAC(2,2)
	YVAP(1) = YVAPAUX/DETERM + YVAP(1)
	YVAP(2) = 1D0 - YVAP(1)
!
	XLIQAUX   = FISOF(1) * RJAC(2,1) - FISOF(2) * RJAC(1,1)
	if (determ.NE.0D0) then

	  XLIQ(1) =  XLIQAUX/DETERM + XLIQ(1)
	  XLIQ(2) = 1D0 - XLIQ(1)

	else

	  XLIQ(1)=-.1d0
	  XLIQ(2)=0d0

	end if
	return
ENDsubroutine !FIN  UPDATE2
!
!-----------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
!	Encuentra nuevos valores de xliq e yvap por sustitucion directa (con la ecuacion de isofugacidad)
!
SUBROUTINE UPDATE5(RESFI1L,RESFI1V,RESFI2L,RESFI2V,RESERVY1,RESERVX1,YVAP,XLIQ)

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XLIQ(2), YVAP(2)

	XLIQ(1) = RESERVY1 * RESFI1V / RESFI1L
	YVAP(2) = (1D0 - RESERVX1) * RESFI2L / RESFI2V
 	XLIQ(2) = 1D0 - XLIQ(1)
	YVAP(1) = 1D0 - YVAP(2)

	return
ENDsubroutine !FIN UPDATE5