	SUBROUTINE ASSOCDCADJUSTFASTSUB(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,TEMP,PSAT,
     @NMAX,IERRDC2)
!
!	Este modulo permite ajustar el diametro de esferas rigidas de un
!	determinado componente de manera de reproducir un punto de la curva
!	de presion de vapor del componente puro. a t y p fijas itera sobre
!	dc hasta lograr isofugacidad. Es mucho mas rapido que dcadjust-sub5
!
!
!
!	Se supone que el dc viene inicializado desde el calling program
!
	IMPLICIT REAL*8 (A-H,O-Z)
	integer, parameter :: NGM = 10
	DIMENSION DC(1), D(1), RN(1)
	DIMENSION  TEMPCRIT(1), PRCRIT(1)
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
C       D(I) es el diametro de una esfera imaginaria de un volumen
C       igual al de un mol de esferas rigidas. Por eso se da por mol (en realidad está en mol^(1/3))
C
C
C       NC  = Numero de componentes (mAX 1)
C
C       RN (I) = NUMERO DE MOLES DEL COMPONENTE I
C
C
C
C
C
	IERRDC2=0
C
C
C CALCULO DE DC, QUE AJUSTA EL PUNTO DE SATURACION
C
	CALL DCCALC2(TEMPCRIT,NG,NU,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @RKASS,DC,TEMP,PSAT,NMAX,IERRDC2)
C
C STOP CRONOMETRO
C
C	CALL LIB$SHOW_TIMER(,2,,)
C
	END
C/////////////////////////////////////////////////////////////////////////
C////////////////////////////////////////////////////////////////////
	SUBROUTINE DCCALC2(TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @                   RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,
     @                   RKASS,DC,TEMP,PSAT,NMAX,IERRDC2)
C
C	A T y P fijas itera sobre
C    	 DC hasta lograr ISOFugacidad. Este esquema de calculo es mucho
C	mas rapido que el basado en sucesivas evaluaciones de la presion
C	 de vapor. Sin embargo puede ser menos robusto en la region critica
C
	IMPLICIT REAL*8 (A-H,O-Z)
	integer, parameter :: NGM = 10
	DIMENSION DC(1), D(1), RN(1)
	DIMENSION TEMPCRIT(1), PRCRIT(1)
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

        dimension FFIFV(1),FFIFL(1)
C
C
	IERRDC2=0
C
	CALL CORRECCION(NC,RN)
CTIZ
C	RN(1)= 1.0
CTIZ
C
C
C CALCULO, PROPIAMENTE DICHO, DE DC
C
C
	CALL DCCALC3(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,PSAT,VLIQ,VVAP,ZFASL,ZFASV,FFIFV,FFIFL,
     @IERRDC,NMAX)
C
	IF(IERRDC.NE.0) THEN
			IERRDC2 = 1
	ELSE
			IERRDC2 = 0
	END IF
C
C FIN DE PROGRAMA PRINCIPAL
C
c	CALL LIB$SHOW_TIMER(,2,,)
	END
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////////////////////////////////////////
C//////////////////////////////////////////////////////////////////////
	SUBROUTINE DCCALC3(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,PSAT,VLIQ,VVAP,ZFASL,ZFASV,FFIFV,FFIFL,
     @IERRDC,NMAX)
C
C NMAX = Maximun Number of Newton Iterations

C
C PSAT = PRESION DE VAPOR, EN ATM
C ZFASL = FACTOR DE COMPRESIBILIDAD FASE LIQUIDA
C ZFASV = FACTOR DE COMPRESIBILIDAD FASE VAPOR
C FFIFV(1) = COEF DE FUGACIDAD EN EL VAPOR
C FFIFL(1) = COEF DE FUGACIDAD EN EL LIQUIDO
C
	IMPLICIT REAL*8 (A-H,O-Z)
	integer, parameter :: NGM = 10
C
	DIMENSION DC(1), RN(1), D(1)
	DIMENSION TEMPCRIT(1)
	DIMENSION NU(NGM,1), QAREA(NGM), TSTR(NGM), GSTR(NGM)
	DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
	DIMENSION ALFAIJ(NGM,NGM)
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,1)
C
	DIMENSION FFIFV(1), FFIFL(1)
C
C
C
C    TEMP = TEMPERATURA EN KELVIN
C
C    VVAP = VOLUMEN DEL VAPOR EN CM3
C
C    VLIQ = VOLUMEN DEL LIQUIDO EN CM3
C
C    PSAT = PRESION DE VAPOR EN ATM
C
C	OPEN(UNIT=2,NAME='DIAM.TXT',FORM='FORMATTED',
C     @ ACCESS='APPEND')
C

C
C
	IERRDC=0
C
	NUMIT = 0
c	LECHUZA = 0
	IERRLECH = 0
        ierrit = 0
C
	RUNIV = 82.05 !atm cm3/(K mol)
C
	DC1=DC(1)
C
	NUITER = 1
C
c 542	ERRCOMP=0
C
C
 449 	CONTINUE

	WRITE(*,'(X,A9,X,F8.4)') "dc(1)1 = ",DC(1)
C
C INICIALIZACION VOLUMEN DEL VAPOR
C
 	CALL VAPINICI(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,PSAT,VVAP,FFIFV,IERRVAP)
C	write(2,*) "VVAP=",VVAP
C	write(2,*) "IERRVAP=",IERRVAP

   	IF(IERRVAP.NE.0) GO TO 419
C
c	IF(ERRCOMP.NE.0) GO TO 132
C
C COMPATIBILIZACION DEL VOLUMEN DEL VAPOR CON LA PSAT
C
c	CALL COMPATVVAP(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
c     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
c     @,ENASS,RKASS,PSAT,VVAP,FFIFV,IERRVC)
c   	IF(IERRVC.NE.0) GO TO 419
C
C INICIALIZACION DEL VOLUMEN DEL LIQUIDO
c
c LOS DOS LLAMADOS QUE VIENEN SON PARA LUEGO CALCULAR EL VALOR
c   INICIAL DEL VOLUMEN LIQUIDO
c
c CALCULO DE D(I) (DIAMETROS DE ESFERAS RIGIDAS A TEMPERATURA TEMP
c
 	CALL SOFT2(NC,DC,TEMPCRIT,TEMP,D)
c
c CALCULO DE LOS LAMBDA NECESARIOS PARA LA CONTRIBUCION REPULSIVA
c
	CALL LAMBDA2(NC,RN,D,RLAMBDA1,RLAMBDA2,RLAMBDA3)
c
 	CALL LIQINICI(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,RLAMBDA3,PSAT,VLIQ,FFIFL,IERRLI)
C	write(2,*) "VLIQ=",VLIQ
C	WRITE(2,*) "IERRLI=",IERRLI
   	IF(IERRLI.NE.0) GO TO 419
C
c	IF(ERRCOMP.NE.0) GO TO 135
C
C COMPATIBILIZACION DEL VOL. LIQ. CON LA PRESION DE SATURACION
C
c	CALL COMPATVLIQ(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
c     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
c     @,ENASS,RKASS,RLAMBDA3,PSAT,VLIQ,FFIFL,IERRLC)
c   	IF(IERRLC.NE.0) GO TO 419
C
  	IF(DABS((VLIQ-VVAP)/VVAP).LT.1.D-3) THEN
		   IERRLECH = 1
C	WRITE(2,*)"IERRLECH=",IERRLECH
		   GO TO 419
	ELSE
	END IF
C
C
C
C CALCULO DEL APARTAMIENTO CON RESPECTO A LA ISOFUG
C
	RISOFAP = 1- FFIFL(1)/FFIFV(1)
C	write(2,*) "RISOFAP=",RISOFAP
C
	IF(NUITER.EQ.1) THEN
      !write (6,*) "VLIQ=", VLIQ
	!write (6,*) "VVAP=", VVAP
	!PAUSE

	  RISOFAP1 = RISOFAP
	  IF(DABS(RISOFAP1).LT.1.D-4) GO TO 555

	  NUITER = 2
	  DC2 = DC1 + 1.D-5*DC1
	  DC(1) = DC2
!	  write(2,*) "DC(1)2=",DC(1)
	  GO TO 449

	END IF
C
	IF(NUITER.EQ.2) THEN
			RISOFAP2 = RISOFAP
C	write (6,*) "VLIQ=", VLIQ
C	write (6,*) "VVAP=", VVAP
	!PAUSE
C			Estimacion utilizando Newton
!	  Paso:
	  dDc = -(DC2 - DC1)*RISOFAP1/(RISOFAP2 - RISOFAP1)
!
!	  Chequeo paso máximo:
	  if(dabs(dDc).GT.1D-1)dDc = 1D-1*(dDc/dabs(dDc))
	  DC3 = dDc + DC1

!	  DC3 = (DC2 - DC1)*(-1.D0)*RISOFAP1
!     @                  / (RISOFAP2 - RISOFAP1)   + DC1
	  NUMIT = NUMIT + 1
C
			IF(NUMIT.GT.NMAX) THEN
					IERRIT = 1
					GO TO 419
			 ELSE
			END IF
C
			IF(DABS(DC3-DC1).LT.1.D-4) GO TO 555
C
			DC1 = DC3
			DC(1) = DC1
			NUITER = 1
			GO TO 449
	 ELSE
	END IF
C
C
C
 555	DC(1) = DC1
C
C CALCULO DE FACTORES DE COMPRESIBILIDAD
C
 	ZFASL = PSAT * VLIQ / RUNIV / TEMP
	ZFASV = PSAT * VVAP / RUNIV / TEMP
C
 419	IF(IERRVAP.NE.0.OR.
     @IERRLI.NE.0.OR.IERRIT.NE.0.
     @OR.IERRLECH.NE.0) IERRDC=1
C
	CLOSE(2)
 	END  !FIN
C////////////////////////////////////////////////////////////////
C/////////////////////////////////////////////////////////////////
C////////////////////////////////////////////////////////////////
