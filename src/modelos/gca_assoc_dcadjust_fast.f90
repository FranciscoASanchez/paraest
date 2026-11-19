      SUBROUTINE ASSOCDCADJUSTFASTSUB(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,TEMP,PSAT,&
     &NMAX,IERRDC2)
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
      integer, parameter :: NCM = 50, NGM = 20
      DIMENSION DC(1), D(1), RN(1)
      DIMENSION  TEMPCRIT(1), PRCRIT(1)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
! DEL TERMINO DE DISPERSION
!
      DIMENSION NU(NGM,1), QAREA(NGM), TSTR(NGM), GSTR(NGM)
      DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
      DIMENSION ALFAIJ(NGM,NGM)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
!
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,1)
!
!
!	TEMP= TEMPERATURA DEL SISTEMA EN K
!
!       D(I)= Diametro de esferas rigidas (funcion de la
!                                          temperatura) (un valor por
!                                                              componente)
!
!       DC(I)= Diametro de esferas rigidas en el punto critico
!                                          (un valor por componente)
!
!       D(I) es el diametro de una esfera imaginaria de un volumen
!       igual al de un mol de esferas rigidas. Por eso se da por mol (en realidad está en mol^(1/3))
!
!
!       NC  = Numero de componentes (mAX 1)
!
!       RN (I) = NUMERO DE MOLES DEL COMPONENTE I
!
!
!
!
!
      IERRDC2=0
!
!
! CALCULO DE DC, QUE AJUSTA EL PUNTO DE SATURACION
!
      CALL DCCALC2(TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,DC,TEMP,PSAT,NMAX,IERRDC2)
!
! STOP CRONOMETRO
!
!	CALL LIB$SHOW_TIMER(,2,,)
!
      END
!/////////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////
      SUBROUTINE DCCALC2(TEMPCRIT,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,&
     &                   RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &                   RKASS,DC,TEMP,PSAT,NMAX,IERRDC2)
!
!	A T y P fijas itera sobre
!    	 DC hasta lograr ISOFugacidad. Este esquema de calculo es mucho
!	mas rapido que el basado en sucesivas evaluaciones de la presion
!	 de vapor. Sin embargo puede ser menos robusto en la region critica
!
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
      DIMENSION DC(1), D(1), RN(1)
      DIMENSION TEMPCRIT(1), PRCRIT(1)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
! DEL TERMINO DE DISPERSION
!
      DIMENSION NU(NGM,1), QAREA(NGM), TSTR(NGM), GSTR(NGM)
      DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
      DIMENSION ALFAIJ(NGM,NGM)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
!
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,1)

        dimension FFIFV(1),FFIFL(1)
!
!
      IERRDC2=0
!
      CALL CORRECCION(NC,RN)
!TIZ
!	RN(1)= 1.0
!TIZ
!
!
! CALCULO, PROPIAMENTE DICHO, DE DC
!
!
      CALL DCCALC3(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,PSAT,VLIQ,VVAP,ZFASL,ZFASV,FFIFV,FFIFL,&
     &IERRDC,NMAX)
!
      IF(IERRDC.NE.0) THEN
      		IERRDC2 = 1
      ELSE
      		IERRDC2 = 0
      END IF
!
! FIN DE PROGRAMA PRINCIPAL
!
!	CALL LIB$SHOW_TIMER(,2,,)
      END
!//////////////////////////////////////////////////////////////////////
!//////////////////////////////////////////////////////////////////////
!//////////////////////////////////////////////////////////////////////
!//////////////////////////////////////////////////////////////////////
      SUBROUTINE DCCALC3(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,PSAT,VLIQ,VVAP,ZFASL,ZFASV,FFIFV,FFIFL,&
     &IERRDC,NMAX)
!
! NMAX = Maximun Number of Newton Iterations

!
! PSAT = PRESION DE VAPOR, EN ATM
! ZFASL = FACTOR DE COMPRESIBILIDAD FASE LIQUIDA
! ZFASV = FACTOR DE COMPRESIBILIDAD FASE VAPOR
! FFIFV(1) = COEF DE FUGACIDAD EN EL VAPOR
! FFIFL(1) = COEF DE FUGACIDAD EN EL LIQUIDO
!
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
!
      DIMENSION DC(1), RN(1), D(1)
      DIMENSION TEMPCRIT(1)
      DIMENSION NU(NGM,1), QAREA(NGM), TSTR(NGM), GSTR(NGM)
      DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
      DIMENSION ALFAIJ(NGM,NGM)
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,1)
!
      DIMENSION FFIFV(1), FFIFL(1)
!
!
!
!    TEMP = TEMPERATURA EN KELVIN
!
!    VVAP = VOLUMEN DEL VAPOR EN CM3
!
!    VLIQ = VOLUMEN DEL LIQUIDO EN CM3
!
!    PSAT = PRESION DE VAPOR EN ATM
!
!	OPEN(UNIT=2,NAME='DIAM.TXT',FORM='FORMATTED',
!     @ ACCESS='APPEND')
!

!
!
      IERRDC=0
!
      NUMIT = 0
!	LECHUZA = 0
      IERRLECH = 0
        ierrit = 0
!
      RUNIV = 82.05 !atm cm3/(K mol)
!
      DC1=DC(1)
!
      NUITER = 1
!
! 542	ERRCOMP=0
!
!
449   CONTINUE

      WRITE(*,'(X,A9,X,F8.4)') "dc(1)1 = ",DC(1)
!
! INICIALIZACION VOLUMEN DEL VAPOR
!
      CALL VAPINICI(NC,DC(1),RN(1),TEMPCRIT,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,PSAT,VVAP,FFIFV,IERRVAP)
!	write(2,*) "VVAP=",VVAP
!	write(2,*) "IERRVAP=",IERRVAP

      IF(IERRVAP.NE.0) GO TO 419
!
!	IF(ERRCOMP.NE.0) GO TO 132
!
! COMPATIBILIZACION DEL VOLUMEN DEL VAPOR CON LA PSAT
!
!	CALL COMPATVVAP(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
!     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
!     @,ENASS,RKASS,PSAT,VVAP,FFIFV,IERRVC)
!   	IF(IERRVC.NE.0) GO TO 419
!
! INICIALIZACION DEL VOLUMEN DEL LIQUIDO
!
! LOS DOS LLAMADOS QUE VIENEN SON PARA LUEGO CALCULAR EL VALOR
!   INICIAL DEL VOLUMEN LIQUIDO
!
! CALCULO DE D(I) (DIAMETROS DE ESFERAS RIGIDAS A TEMPERATURA TEMP
!
      CALL SOFT2(NC,DC,TEMPCRIT,TEMP,D)
!
! CALCULO DE LOS LAMBDA NECESARIOS PARA LA CONTRIBUCION REPULSIVA
!
      CALL LAMBDA2(NC,RN,D,RLAMBDA1,RLAMBDA2,RLAMBDA3)
!
      CALL LIQINICI(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,RLAMBDA3,PSAT,VLIQ,FFIFL,IERRLI)
!	write(2,*) "VLIQ=",VLIQ
!	WRITE(2,*) "IERRLI=",IERRLI
      IF(IERRLI.NE.0) GO TO 419
!
!	IF(ERRCOMP.NE.0) GO TO 135
!
! COMPATIBILIZACION DEL VOL. LIQ. CON LA PRESION DE SATURACION
!
!	CALL COMPATVLIQ(NC,DC,RN,TEMPCRIT,TEMP,NG,NU,QAREA
!     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
!     @,ENASS,RKASS,RLAMBDA3,PSAT,VLIQ,FFIFL,IERRLC)
!   	IF(IERRLC.NE.0) GO TO 419
!
      IF(DABS((VLIQ-VVAP)/VVAP).LT.1.D-3) THEN
      	   IERRLECH = 1
!	WRITE(2,*)"IERRLECH=",IERRLECH
      	   GO TO 419
      ELSE
      END IF
!
!
!
! CALCULO DEL APARTAMIENTO CON RESPECTO A LA ISOFUG
!
      RISOFAP = 1- FFIFL(1)/FFIFV(1)
!	write(2,*) "RISOFAP=",RISOFAP
!
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
!
      IF(NUITER.EQ.2) THEN
      		RISOFAP2 = RISOFAP
!	write (6,*) "VLIQ=", VLIQ
!	write (6,*) "VVAP=", VVAP
	!PAUSE
!			Estimacion utilizando Newton
!	  Paso:
        dDc = -(DC2 - DC1)*RISOFAP1/(RISOFAP2 - RISOFAP1)
!
!	  Chequeo paso máximo:
        if(dabs(dDc).GT.1D-1)dDc = 1D-1*(dDc/dabs(dDc))
        DC3 = dDc + DC1

!	  DC3 = (DC2 - DC1)*(-1.D0)*RISOFAP1
!     @                  / (RISOFAP2 - RISOFAP1)   + DC1
        NUMIT = NUMIT + 1
!
      		IF(NUMIT.GT.NMAX) THEN
      				IERRIT = 1
      				GO TO 419
      		 ELSE
      		END IF
!
      		IF(DABS(DC3-DC1).LT.1.D-4) GO TO 555
!
      		DC1 = DC3
      		DC(1) = DC1
      		NUITER = 1
      		GO TO 449
       ELSE
      END IF
!
!
!
555   DC(1) = DC1
!
! CALCULO DE FACTORES DE COMPRESIBILIDAD
!
      ZFASL = PSAT * VLIQ / RUNIV / TEMP
      ZFASV = PSAT * VVAP / RUNIV / TEMP
!
419   IF(IERRVAP.NE.0.OR.&
     &IERRLI.NE.0.OR.IERRIT.NE.0.&
     &OR.IERRLECH.NE.0) IERRDC=1
!
      CLOSE(2)
      END  !FIN
!////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////
