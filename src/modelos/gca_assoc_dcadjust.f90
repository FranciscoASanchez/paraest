      SUBROUTINE ASSOCDCADJUSTSUB(NC,DC,TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,TEMP,PSAT,NMAX,PSATCALC)
!
!  MODELO ASSOC-GC-EOS
!
! no contiene modulo de inicializacion de dc. Se supone que dc viene
!  inicializado desde el calling program
!
! DC(1):
!
! on input dc contains the initial value of dc
! on output dc contains the final value of dc
!
!	' PSAT, ATM ='
!
! ESTE PROGRAMA PERMITE AJUSTAR EL DIAMETRO DE ESFERAS RIGIDAS DE UN
! DETERMINADO COMPONENTE DE MANERA DE REPRODUCIR UN PUNTO DE LA CURVA
!  DE PRESION DE VAPOR DEL COMPONENTE PURO.
!

! Las coordenadas del punto de la curva de PVAP a reproducir
!     son  (TEMP,PSAT)   en KELVIN y ATM RESPECTIVAMENTE
!
!    PSATCALC = Presion de Vapor calculada (en ATM) calculada con
!               el valor final de DC. Se supone que es muy parecida
!	        a PSAT. (a la temperatura TEMP)
!
! USA NEWTON para la primer
!  iteracion. EN TODAS LA ITERACIONES SUBSIGUIENTES USA UNA PARABOLA
!
!  CONTIENE UNA OPCION QUE PERMITE LIMITAR EL NUMERO DE ITERACIONES
!
!  LLAMA A UNA SUBRUTINA QUE CALCULA PRESIONES DE VAPOR SIN REQUERIR
!   INICIALIZACION EN LA PRESION DE VAPOR
!
!  NMAX = NUMERO MAXIMO DE ITERACIONES PERMITIDO
!        SI NMAX = 1 SOLO SE HACE UN PASO DE NEWTON
!        SI NMAX = 2 O MAYOR ,   SE HACE UN PASO DE NEWTON
!                     Y TODOS LOS SIGUIENTES USAN UNA PARABOLA, EN LUGAR
!                      DE LA RECTA TANGENTE DE NEWTON
!        HASTA HOY (18-AGO-1993) LA EXPERIENCIA MUESTRA QUE CON NMAX = 2
!          LA PRECISION ES BUENA.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20

      DIMENSION DC(1), D(1), RN(1)
      DIMENSION  TEMPCRIT(1)
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
!           D(I) es el diametro de una esfera imaginaria de un volumen
!                igual al de un mol de esferas rigidas. Por eso se da por
!                                                                  mol
!
!       NC  = Numero de componentes (mAX 1)
!
!       RN (I) = NUMERO DE MOLES DEL COMPONENTE I
!
!
!
!
!
!
! CALCULO DE DC, QUE AJUSTA EL PUNTO DE SATURACION
!
      CALL DCCALC(TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,DC,TEMP,PSAT,PSATCALC,NMAX)
!
!	CALL LIB$SHOW_TIMER(,2,,)
!
      ENDsubroutine
!/////////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////
      SUBROUTINE DCCALC(TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,DC,TEMP,PSAT,PSATCALC,NMAX)
!
! CALCULO DE DC
!
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20

      DIMENSION DC(1), D(1), RN(1)
      DIMENSION  TEMPCRIT(1)
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
! ARREGLOS AUXILIARES
!
      DIMENSION X(2),Y(2)
!
      NITERP=0
      NINFA=0
      MINFA=0
!
!
! GUARDO EL VALOR DE PSAT Y DE DC
!
      PSAT0 = PSAT
!
      DC0 = DC(1)
!
!
      DC1 = DC0
!
321   DC(1) = DC1
!
      CALL ASSOCPSATSUB(TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,DC,TEMP,PSAT,ZFASL,ZFASV,IERRPSAT)
!
      PSAT1 = PSAT
      IF(NITERP.GE.NMAX) GO TO 999
      MINFA = MINFA + 1
      NINFA = NINFA + 1
      X(NINFA) = PSAT1
      Y(NINFA) = DC1
!
      DC2 = DC1 + 1.D-5*DC1
!
      DC(1) = DC2
!
      CALL ASSOCPSATSUB(TEMPCRIT,NG,NU,QAREA,TSTR,&
     &GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,&
     &RKASS,DC,TEMP,PSAT,ZFASL,ZFASV,IERRPSAT)
!
      PSAT2 = PSAT
!
      if (psat2 == psat1) then
        psat2 = psat1*dc2/dc1
      endif
      YPRIMA2 = (DC2-DC1)/(PSAT2-PSAT1)
!
      IF(MINFA.EQ.1) THEN
!                       Newton
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
!
      IF(DABS(DC3-DC1).LT.1.D-8.AND.DABS(PSAT1-PSAT0).LT.1.D-8) THEN
      				GO TO 999
      			   ELSE
      				DC1 = DC3
      				GO TO 321
      END IF
!
!
!
!
!
!
999   PSATCALC = PSAT1
      PSAT = PSAT0
      ENDSUBROUTINE DCCALC
!////////////////////////////////////////////////////////
      SUBROUTINE PARAB(X,Y,YPRIMA2,A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(2),Y(2)
!
! CALCULA LOS COEFICIENTES DE UNA PARABOLA PASANTE POR 2 PUNTOS
!   CUYA DERIVADA EN EL PUNTO X(2),Y(2) COINCIDE CON EL VALOR
!     ALIMENTADO EN YPRIMA2
!
!	        2
!	 Y = A X  +  B X  +  C
!
      A21 = Y(2) - Y(1)
      B21 = X(2) - X(1)
      C21 = X(2)**2 - X(1)**2
!
      A = (A21 - YPRIMA2 * B21)/ (C21 - 2 * B21 * X(2))
!
      B = YPRIMA2 - 2 * A * X(2)
!
      C = Y(2) - A * X(2)**2 - B * X(2)
!
!
      ENDSUBROUTINE PARAB
!//////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////

