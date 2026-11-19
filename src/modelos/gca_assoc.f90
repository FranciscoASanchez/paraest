!
!	Paquete termodinámico A-GC-EOS versión 3.0.1. Este paquete es el utilizado por el
!	programa Paraest, y es fundamentalmente distinto del utilizado por programas como
!	GC3, GCEnvel, GCDist, THDist, GCExtract, GCXS, GCDB, etc.
!
!
!	Control de cambios:
!	v3.0.0: Versión tal cual lo encontré de entrada al llegar a PLAPIQUI en abril de
!	        2008.
!
!	v3.0.1: Corrección de errores. Se amplio la presición de diversos números presen-
!	        tes, ej: 1.0 -> 1.0D0
!	        Acomodo de identaciones en algunas subroutinas.
!	        Eliminación de parte de código antigo y comentado. Si se desea recuperar,
!	        ver el archivo "A-gc-eso3.for".
!
!
!
!
!
      SUBROUTINE ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,&
     &TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)



	!MODIFICAR CROSSASOC4,DCROSSASOC4,DCROSSASOCVL4
!
!	NCOD es un codigo que especifica que' es lo que se quiere calcular:
!		(se introdujo NCOD para ahorrar tiempo de computo)
!
!     NCOD      SIGNIFICADO
!
!      1         factor de compresibilidad solamente
!      2         coeficientes de fugacidad con factor de compresibilidad
!                   entrado como dato (no es posible calcular el coef.
!                    de fugacidad si no se conoce el factor de compres)
!                    Por supuesto el volumen y la temperatura tb. entran como
!			dato
!      3         coeficientes de fugacidad y el necesario factor de
!	          compresibilidad (calculados ambos dentro de esta rutina)
!      4	 Energia de Helmholtz (util para chequeo de derivadas
!	 	analiticas a traves de derivadas numericas)
!      5         factor de compresibilidad, coeficientes de fugacidad y
!		 energia de Helmholtz (es decir todo)
!
!
!
!
!    ARESTOT      = Energia residual total
!
!                      R
!    ARESTOT      = ( A  / RT )         (las unidades son "moles")
!
!
! NGA = NRO DE GRUPOS DISTINTOS EN EL TERMINO DE ASOCIACION
! NU (J,L)= NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO L (PARA EL
!                                       TERMINO DISPERSIVO)
! ONUASOC(J,L) = NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO L (PARA
!                               EL TERMINO DE ASOCIACION)
!
! CALCULA EL FACTOR DE COMPRESIBILIDAD (ZSIS) Y LOS COEFICIENTES DE
! FUGACIDAD DE CADA COMPONENTE [FFI(I)] y la egia residual total, EN
! FUNCION DEL NUMERO DE MOLES  [RN(I)], LA TEMPERATURA(TEMP) Y EL
! VOLUMEN DEL SISTEMA (V)
!
! EL MODELO ES ASSOC-GC-EOS (NO CONFUNDIR CON ASS-GC-EOS). LA
!  CONTRIBUCIONES repulsiva y dispersiva son identicas a las de la
!  GCEOS. El te'rmino de asociacion es a CONTRIBUCION GRUPAL
!   EL NRO. Y TIPO DE GRUPOS PARA UNA DADA MOLECULA ES DISTINTO
!     SEGUN SE TRATE DE LA CONTRIBUCION DISPERSIVA O DE ASOCIACION
!
!        IMPLICIT REAL*8 (A-H,O-Z)

      IMPLICIT NONE
      integer, parameter :: NCM = 50, NGM = 20

!
! SE manejan:
!               NC componentes (aunque en realidad el maximo es NCM)
!               NGM grupos diferentes desde el punto de vista del
!                                                       termino dispersivo
!                       12  grupos diferentes desde el punto de vista del
!                                                       termino asociativo
!                        6 sitios de  asociacion por grupo
!
! DC tiene dimension NC porque DC es un argumento de la subrutina
! D   es un arreglo interno , por lo que su dimension es (debe ser) fija e
!       igual al maximo nro. de componentes
!
      INTEGER NC
      DOUBLE PRECISION DC(NC), D(NCM), RN(NC)
      DOUBLE PRECISION  TEMPCRIT(NC)
!
! DOUBLE PRECISIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
! DEL TERMINO DE DISPERSION
!
        INTEGER NU(NGM,NC)
        DOUBLE PRECISION  QAREA(NGM), TSTR(NGM), GSTR(NGM)
        DOUBLE PRECISION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM),&
     &                   RKTIJ(NGM,NGM)
        DOUBLE PRECISION ALFAIJ(NGM,NGM)
        DOUBLE PRECISION TITA(NGM), TAU(NGM,NGM)
        DOUBLE PRECISION TSPLINE(NGM)
        DOUBLE PRECISION GGG(NGM,NGM)
        DOUBLE PRECISION VERDKIJ(NGM,NGM)
        DOUBLE PRECISION DLTG(NGM,NGM)
!
! DOUBLE PRECISIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
!
        DOUBLE PRECISION ONUASOC(12,NC)
        DOUBLE PRECISION RNGRUP(12)  !nro. de moles de los grupos
!
! MASSOC(J) = NRO. DE SITIOS DE ASOCIACION ASIGNADOS AL GRUPO J
! ENASS(I,J,K,L)= PARAMETRO DE ENERGIA DE ASOCIACION entre el sitio
!                 I del GRUPO J, y el sitio K del grupo L
! RKASS(I,J,K,L)= PARAMETRO DE VOLUMEN DE ASOCIACION entre el sitio
!                 I del GRUPO J, y el sitio K del grupo L
!
! DELTASS(I,J,K,L) = FUERZA ADOUBLE PRECISIONAL DE ASOCIACION ENTRE  el sitio
!                 I del GRUPO J, y el sitio K del grupo L
!
!
        INTEGER MASSOC(12)
        DOUBLE PRECISION ENASS(6,12,6,12), RKASS(6,12,6,12)
        DOUBLE PRECISION DELTASS(6,12,6,12)
        DOUBLE PRECISION XNOASS(6,12)
        DOUBLE PRECISION FASS(12)
        DOUBLE PRECISION FAASS(6,12)
!
! DOUBLE PRECISIONAMIENTO RELAC CON COEF FUGACIDAD REPULSIVO
!
        DOUBLE PRECISION  DL1N(NCM),DL2N(NCM),DL3N(NCM),DYN(NCM)
        DOUBLE PRECISION FUGREP(NCM)
!
! DOUBLE PRECISIONAMIENTO DE VARIABLES PARA EL CALCULO
!    DEL COEFICIENTE DE FUGACIDAD DE DISPERSIVO
!
        DOUBLE PRECISION PS(NCM,NGM), RMS(NCM), H2(NGM)
        DOUBLE PRECISION H3(NCM,NGM), H4(NGM), H5(NGM)
        DOUBLE PRECISION H6(NGM), H7(NCM,NGM)
        DOUBLE PRECISION FUGDISP(NCM)
!
! DOUBLE PRECISIONAMIENTO RELAC. CON COEF. FUG. ASOCIATIVO
!
        DOUBLE PRECISION DROSSN(12,NCM)
        DOUBLE PRECISION DXNOASSN(6,12,NCM)
        DOUBLE PRECISION DFAASSXNOASS(6,12)
        DOUBLE PRECISION DFASSN(12,NCM)
        DOUBLE PRECISION FUGASS(NCM)
!
! DOUBLE PRECISIONAMIENTO RELAC. CON FACTOR DE COMPRES. ASOCIATIVO
!
        DOUBLE PRECISION DROSSV(12)
        DOUBLE PRECISION DXNOASSV(6,12)
        DOUBLE PRECISION DFASSV(12)
!
        DOUBLE PRECISION FFI(NC)  !DOUBLE PRECISIONAMIENTO COEFICIENTES DE FUGACIDAD
!
!
!
        DOUBLE PRECISION FUGSIS(NCM)

        INTEGER NCOD,IERRFUG,IERRZSIS,IERRREP,NG,NGA,L,NGA2
        DOUBLE PRECISION ARESTOT,ZSIS,ARESASS,ARESDISP,ARESFV,QTOT,&
     &RLAMBDA1,RLAMBDA2,RLAMBDA3,RNTOT,Y,ZASS,ZDISP,ZFV,TEMP,V

        INTEGER ierrassoc, J, KK
        DOUBLE PRECISION TEMP1
!
!
!       TEMP= TEMPERATURA DEL SISTEMA EN K
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
!       V = VOLUMEN DEL SISTEMA EN CM3
!
!       NC  = Numero de componentes (mAX NCM)
!
!       RN (I) = NUMERO DE MOLES DEL COMPONENTE I
!
!
!  V = VOLUMEN DEL SISTEMA
!      V debe tener unidades congruentes con las de D. Es decir si el
!      diametro de esferas rigidas esta dado en cm/mol el volumen V debe
!      estar en cm3.
!
!  Y = VARIABLE INTERMEDIA EN LA EXPRESION DE MANSOORI Y LELAND (1972)
!
!
!  PARA EVITAR ERRORES
!
      IF(NCOD.LT.1.OR.NCOD.GT.5) NCOD = 5  !CALCULA TODO
!
!  AQUI SE OBLIGA A TODOS LOS ARREGLOS QUE NO SON DATO A TENER VALORES
!  NULOS
!
! AQUI SE LIMPIAN TODAS LAS VARIABLES DE SALIDA
!       Excepcion: Si NCOD es 2 ZSIS entra como dato. Es decir que en
!		    ese caso no se limpia ZSIS
!
      CALL ACLAR(NCOD,NC,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,&
     &ierrassoc)
!
! Limpieza de variables locales desde el punto de vista de la rutina principal
! (no de las sub-sub-rutinas
!
        CALL CLEARANCE(NC,NG,NGA,ARESASS,ARESDISP,ARESFV,QTOT,&
     &RLAMBDA1,RLAMBDA2,RLAMBDA3,RNTOT,Y,ZASS,ZDISP,ZFV,&
     &MASSOC,D,DL1N,DL2N,DL3N,DYN,FUGASS,FUGDISP,FUGREP,FUGSIS,&
     &RMS,H3,H7,PS,DFASSN,DXNOASSN,H2,H4,H5,H6,TITA,TSPLINE,&
     &DLTG,TAU,VERDKIJ,GGG,DFASSV,FASS,RNGRUP,DXNOASSV,FAASS,XNOASS,&
     &DFAASSXNOASS,DROSSV,DELTASS,DROSSN)
!
! CALCULO DEL NUMERO TOTAL DE MOLES (RNTOT=NUMERO TOTAL DE MOLES)

!	OPEN(UNIT=4,NAME='DESVIOS.TXT',FORM='FORMATTED')

!	 WRITE(4,*) "TEMP,V",TEMP,V
      if(V.lt.1d-13) then
      	IERRZSIS=1
      	goto 1000
      end if
!
        CALL SUMATORIA(NC,RN,RNTOT)
!
! CALCULO DEL NUMERO TOTAL DE SEGMENTOS DE AREA (PARA CONTRIB. POR DISP)
!
        CALL SURFSEG(NC,NG,NU,RN,QAREA,QTOT)
!
! CALCULO DE LAS FRACCIONES DE SUPERFICIE (PARA CONTRIB. POR DISP.)
!
        CALL FRACAREA(NC,NG,NU,RN,QAREA,QTOT,TITA)
!
! CALCULO DE LA TEMPERATURA A LA QUE SE REALIZA EL SPLINE
!
        CALL SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
!
! CALCULO DE PARAMETROS ENERGETICOS ATRACTIVOS (LIKE)
!
        CALL CALGGG(NG,TEMP,TSTR,GSTR,GDOT,GDDOT,TSPLINE,GGG)
!
! CALCULO DE LOS KIJ (VERDKIJ==>VERDADERO KIJ)
!
        CALL CALCKIJ(NG,RKIJ,RKTIJ,TSTR,VERDKIJ,TEMP)
!
!
! CALCULO DE PARAMETROS ENERGETICOS ATRACTIVOS (unLIKE)
!
        CALL CALGGGIJ(NG,VERDKIJ,GGG)
!
! CALCULO DE DELTAgij
!
        CALL DELTADISP(NG,GGG,DLTG)
!
! CALCULO DE TAU(J,I) (NECESARIO PARA CONTRIB.  FZAS DISPERSION

!
        CALL TAUCALC(NG,TEMP,V,ALFAIJ,QTOT,DLTG,TAU)
!
! CALCULO DE LA ENERGIA RESIDUAL DE HELMHOLTZ DE DISPERSION
!
! SEGUN GC-EOS ORIGINAL
!
        IF(NCOD.GT.3) CALL &
     &AHELMRESDISP(NC,NG,RN,V,TEMP,NU,TITA,TAU,QTOT,&
     &GGG,QAREA,ARESDISP)
!
!
! CALCULO DE D(I) (DIAMETROS DE ESFERAS RIGIDAS A TEMPERATURA TEMP
!
        CALL SOFT(NC,DC,TEMPCRIT,TEMP,D)
!
! CALCULO DE LOS LAMBDA NECESARIOS PARA LA CONTRIBUCION REPULSIVA
!
        CALL LAMBDA(NC,RN,D,RLAMBDA1,RLAMBDA2,RLAMBDA3)
!
! CALCULO DE LA VARIABLE INTERMEDIA Y, NECESARIA PARA LA CONTRIB.
!  RPULSIVA
!
        CALL YCALC(V,RLAMBDA3,Y,IERRREP)
!
! TERMINA LA EJECUCION SI LA FRACCION DE VOLUMEN OCUPADA POR LAS
!  ESFERAS BLANDAS ES MAYOR A LA UNIDAD ( ES DECIR SI EL VOLUMEN
!                                             ES MUY CHICO)
!
        IF(IERRREP.EQ.1) GO TO 1000
!
!  CALCULO DE LA ENERGIA RESIDUAL DE HELMHOLTZ DE FREE VOLUME SEGUN
! GC-EOS ORIGINAL
!
!
        IF(NCOD.GT.3)  CALL &
     &AHELMRESFV(RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ARESFV)
!

      If (nc==1) then

        do l=1,nga
	    !print *, "ONUASOC(1,1)=", ONUASOC(l,1)
          if(ONUASOC(l,1).ne.0) then

            nga2 = nga
            go to 1003

          end if

        end do
        nga2=0
1003    continue

      else

        nga2=nga

      end if



!   CALCULO DEL NUMERO DE MOLES DE CADA GRUPO EN LA CONTRIBUCION DE
!               ASOCIACION
!

         IF(nga2.NE.0) THEN  !AVOID ASSOC 1

        CALL ASSOCNMOL(NC,nga2,RN,ONUASOC,RNGRUP)
!
! Calculo del Delta de asociacion (para nosotros es adimensional
!  porque incluye la densidad de cierto grupo)
!
!
        CALL CALCDELTASS(nga2,TEMP,MASSOC,ENASS,RKASS,DELTASS)
!
! CALCULO DE LAS FRACCIONES DE SITIOS NO ASOCIADOS (XNOASS)
!
        CALL FRACNOAS(nga2,MASSOC,DELTASS,XNOASS, nc, ONUASOC,V,RNGRUP&
     &,ierrassoc)


      TEMP1=RNGRUP(1)



!
! CALCULO DE LA ENERGIA DE HELMHOLTZ RESIDUAL DE ASOCIACION
!
!  FASS(I) Y FAASS(I) SON  VARIABLES INTERMEDIAS DE LAS QUE MAS ADELANTE
!                       SE OBTIENEN DERIVADAS
!
        IF(NCOD.NE.1)  CALL &
     &      AHELMRESASS(nga2,MASSOC,RNGRUP,XNOASS,ARESASS,FASS,FAASS)
!
      END IF  !AVOID ASSOC 1
!
! CALCULO DE LA ENERGIA RESIDUAL TOTAL
!
      IF(NCOD.GT.3)  CALL AHELMRESTOT(ARESFV,ARESDISP,ARESASS,ARESTOT)
!
      IF(NCOD.EQ.4) GO TO 1000
!
! CALCULOS AUXILIARES NECESARIOS PARA EL CALCULO DE LA
!  CONTRIBUCION REPULSIVA AL COEFICIENTE DE FUGACIDAD
!
        IF(NCOD.NE.1)  CALL DERIVFV(NC,V,D,Y,DL1N,DL2N,DL3N,DYN)
!
! CALCULO DE LA CONTRIBUCION REPULSIVA AL COEFICIENTE DE FUGACIDAD
!
        IF(NCOD.NE.1)  CALL &
     &FUGAREP(NC,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,&
     &DL1N,DL2N,DL3N,DYN,FUGREP)

!
! CALCULOS AUXILIARES NECESARIOS PARA EL CALCULO DE LA
!  CONTRIBUCION DISPERSIVA AL COEFICIENTE DE FUGACIDAD
!
        CALL HCALC(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,&
     &TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
!
! CALCULO DE LA CONTRIBUCION DISPERSIVA AL COEFICIENTE DE FUGACIDAD
!
        IF(NCOD.NE.1)  CALL &
     &FUGADISP(NC,NG,TITA,PS,RMS,H2,H3,H4,H5,H6,H7,FUGDISP)

!
!
      IF(nga2.NE.0) THEN  !AVOID ASSOC 2
!
! CALCULO DE DERIVADAS DE LA DENSIDAD CON RESPECTO
!    AL NUMERO DE MOLES
!
        IF(NCOD.NE.1)  CALL &
     &DERIROSSN(NC,nga2,MASSOC,ONUASOC,RNGRUP,DROSSN,V)

!
! CALCULO DE LAS DERIVADAS DE LAS FRACCIONES NO ASOCIADAS CON RESPECTO
!  AL NUMERO DE MOLES entra DROSSN e sai DXNOASSN
!
        IF(NCOD.NE.1)  CALL DERIXNOASS(NC,nga2,MASSOC,DELTASS,XNOASS,&
     &DROSSN,DXNOASSN,ONUASOC,RNGRUP,V)
	!		 write (6,*) "cheguei"

!
!
! CALCULO DE LAS DERIVADAS DE FAASS(I) CON RESPECTO A LAS FRACCIONES NO
!       ASOCIADAS
!
        CALL DERIFAASS(nga2,MASSOC,XNOASS,DFAASSXNOASS)
!
! CALCULO DE LAS DERIVADAS DE FASS(I) CON RESPECTO A LOS NUMEROS DE MOLES
!
        IF(NCOD.NE.1)  CALL &
     &DERIFASSN(NC,nga2,MASSOC,DFAASSXNOASS,DXNOASSN,DFASSN)
!
! CALCULO DE LA CONTRIBUCION ASOCIATIVA AL COEFICIENTE DE FUGACIDAD
!
        IF(NCOD.NE.1)  CALL &
     &FUGAASS(NC,nga2,MASSOC,ONUASOC,RNGRUP,XNOASS,FUGASS)

!
      END IF  !AVOID ASSOC 2
!
! CALCULO DE LA CONTRIBUCION  DE VOLUMEN LIBRE AL FACTOR DE COMPRESIBILIDAD
!
        IF(NCOD.NE.2) CALL &
     &COMPRESFV(V,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ZFV)
!
! CALCULO DE LA CONTRIBUCION DISPERSIVA AL FACTOR DE COMPRESIBILIDAD
!
        IF(NCOD.NE.2)   CALL &
     &COMPRESDISP(NC,NG,RNTOT,RN,NU,QAREA,H2,H4,H5,H6,ZDISP)
!
!
      IF(nga2.NE.0) THEN  !AVOID ASSOC 3
!
! CALCULO DE DERIVADAS DE LA FZA. DE ASOCIACION CON RESPECTO
!    AL VOLUMEN
!
        IF(NCOD.NE.2)CALL DERIROSSV(nga2,V,MASSOC,DELTASS,DROSSV,RNGRUP)
!
! CALCULO DE LAS DERIVADAS DE LAS FRACCIONES NO ASOCIADAS CON RESPECTO
!                       AL VOLUMEN
!
        IF(NCOD.NE.2)  CALL DERIXNOASSVL(nga2,MASSOC,DELTASS,XNOASS,&
     &DROSSV,DXNOASSV, ONUASOC, nc,RNGRUP,V)
!
! CALCULO DE LAS DERIVADAS DE FASS(I) CON RESPECTO AL VOLUMEN
!
        IF(NCOD.NE.2)  CALL &
     &DERIFASSV(nga2,MASSOC,DFAASSXNOASS,DXNOASSV,DFASSV)
!
!
! CALCULO DE LA CONTRIBUCION ASOCIATIVA AL FACTOR DE COMPRESIBILIDAD
!
        IF(NCOD.NE.2)  CALL COMPRESASS(nga2,V,RNTOT,RNGRUP,DFASSV,ZASS)
!
!
       ELSE
        CONTINUE
      END IF  !AVOID ASSOC 3
!
! CALCULO DEL FACTOR DE COMPRESIBILIDAD DEL SISTEMA
!
        IF(NCOD.NE.2)  CALL FCTRCOMPRES(ZFV,ZDISP,ZASS,ZSIS,IERRZSIS)
!	WRITE(4,*)"ZFV,ZDISP,ZASS,ZSIS=",ZFV,ZDISP,ZASS,ZSIS
!
! TERMINA LA EJECUCION SI EL FACTOR DE COMPRESIBILIDAD ES MENOR QUE CERO
!
        IF(IERRZSIS.EQ.1) GO TO 1000
!
! CALCULO DE LOS COEFICIENTES DE FUGACIDAD
!
        IF(NCOD.NE.1)  CALL &
     &FUGASIS(NC,FUGREP,FUGDISP,FUGASS,ZSIS,FUGSIS,FFI,IERRFUG)
!	write(4,*)"rntot=", RNTOT
!	write(4,*)"rnGRUP=", RNGRUP(1)
!	write(4,*)"xnoass(1,1)=", XNOASS(1,1)
!	WRITE(4,*)"FUGREP,FUGDISP,FUGASS,FUGSIS,ZSIS",
!     @FUGREP(1),FUGDISP(1),FUGASS(1),FUGSIS(1),ZSIS
!	WRITE(4,*)FUGREP(2),FUGDISP(2),FUGASS(2),FUGSIS(2),ZSIS

! FINAL DE LA SUBRUTINA ASSOC-GC-EOS
!
1000    END
!
!/////////////////////////////////////////////////////////////////
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
      SUBROUTINE SUMATORIA(NC,RN,RNTOT)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION RN(NC)
!
      rNtot = sum(rN(:NC))
      !RNTOT=0.D0
	!DO I=1,NC
	  !RNTOT=RNTOT+RN(I)
      !END DO
!
      return
      ENDsubroutine
!///////////////////////////////////////////////////////////////////
!//////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////
      SUBROUTINE AHELMRESFV(RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ARESFV)
!
      Implicit real*8 (A-H, O-Z)

!-----------------------------------------------------------------------
!
!       ARESFV= CONTRIBUCION REPULSIVA (O DE VOLUMEN LIBRE) A LA ENERGIA
!                RESIDUAL DE HELMHOLTZ =
!
!                      R
!                 = ( A  / RT )         (las unidades son "moles")
!                              fv
!
!----------------------------------------------------------------
!
!
!  CALCULO DE ARESFV
!
!
        ARESFV=3D0*(RLAMBDA1*RLAMBDA2/RLAMBDA3)*(Y-1D0)
!
        ARESFV= ARESFV + (RLAMBDA2**3/RLAMBDA3**2)*(Y**2-Y-DLOG(Y))
!
        ARESFV= ARESFV + RNTOT*DLOG(Y)
!
      return
      ENDsubroutine
!///////////////////////////////////////////////////////////////////
      SUBROUTINE SOFT(NC,DC,TEMPCRIT,TEMP,D)
!
! CAlCULA LOS DIAMETROS DE ESFERAS RIGIDAS A LA TEMP DEL SISTEMA
!
        Implicit real*8 (A-H, O-Z)
        integer, parameter :: NCM = 50
        DIMENSION DC(NC), TEMPCRIT(NC), D(NCM)
!
        HELPDC = 0d0
!
! CALCULO DE LOS D(I)
!
        DO I=1,NC
                HELPDC=DEXP(-2D0*TEMPCRIT(I)/3D0/TEMP)
                D(I)=1.065655d0*DC(I)*(1D0-0.12D0*HELPDC)
        END DO
!
        END
!///////////////////////////////////////////////////////////////////
        SUBROUTINE LAMBDA(NC,RN,D,RLAMBDA1,RLAMBDA2,RLAMBDA3)
        Implicit real*8 (A-H, O-Z)
        integer, parameter :: NCM = 50
        DIMENSION D(NCM),RN(NC)
!
! CALCULO DE LOS LAMBDA NECESARIOS PARA LA CONTRIBUCION REPULSIVA
!
!
!  RLAMBDA3 ES PROPORCIONAL AL VOLUMEN OCUPADO POR LAS ESFERAS PROPIAMENTE
!   DICHAS, EXCLUYENDO EL ESPACIO INTERSTICIAL
!
      RLAMBDA1=0.D0
      RLAMBDA2=0.D0
      RLAMBDA3=0.D0
!
      DO I=1,NC

        aux = RN(I)*D(I)
        RLAMBDA1= aux + RLAMBDA1
        RLAMBDA2= aux*d(i) + RLAMBDA2
        RLAMBDA3= aux*D(I)**2  +  RLAMBDA3

      END DO
!
      return
      ENDsubroutine
!//////////////////////////////////////////////////////////////////////
        SUBROUTINE YCALC(V,RLAMBDA3,Y,IERRREP)
        Implicit real*8 (A-H, O-Z)
!
        HELPY = 0.D0
!
!  Y = VARIABLE INTERMEDIA EN LA EXPRESION DE MANSOORI Y LELAND (1972)
!----------------------------------------------------------------
        HELPY = 3.141592654D0*RLAMBDA3/(6.D0*V)
!	WRITE(6,*)"V=",V
!	WRITE(6,*)"RLAMDA3=",V

!        WRITE(6,*) ' '
        IF(HELPY.GT.(1.)) THEN
      	WRITE(6,*) ' V ES DEMASIADO CHICO'
!	PAUSE
         END IF
!        WRITE(6,*) ' '
        HELPY = 1 - HELPY
        Y = 1/HELPY
!
        IF(Y.LT.0) THEN
                        IERRREP = 1
                   ELSE
                        IERRREP = 0
        END IF
!
        END
!/////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////
        SUBROUTINE SURFSEG(NC,NG,NU,RN,QAREA,QTOT)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION RN(NC)
!
        DIMENSION QAREA(NGM)
        DIMENSION NU(NGM,NC)
!
        HELPQ = 0.D0
!
! QTOT=NUMERO TOTAL DE SEGMENTOS DE AREA
!
! NU (I,J)= NUMERO DE GRUPOS DEL TIPO I EN EL COMPUESTO DE TIPO J
! NU (J,I)= NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO I
!
        QTOT=0.D0
!
        DO I=1,NC
                HELPQ=0.D0
                  DO J=1,NG
                        HELPQ=HELPQ + NU(J,I)*QAREA(J)
                  END DO
                QTOT=QTOT + RN(I)*HELPQ
        END DO
!
        END
!/////////////////////////////////////////////////////////////////////
        SUBROUTINE FRACAREA(NC,NG,NU,RN,QAREA,QTOT,TITA)
!
! CALCULA LAS FRACCIONES DE AREA DE TODOS LOS GRUPOS
!
!  TITA(J)= FRACCION DE AREA DEL GRUPO J
!
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION RN(NC)
!
        DIMENSION QAREA(NGM)
        DIMENSION NU(NGM,NC)
        DIMENSION TITA(NGM)
!
        HELPFRAC = 0.D0
!
        DO J=1,NG
               HELPFRAC=0.D0
                DO I=1,NC
                  HELPFRAC=HELPFRAC + RN(I)*NU(J,I)
                END DO
               TITA(J)=HELPFRAC*QAREA(J)/QTOT
        END DO
!
        END
!/--------------------------------------------------------------------------------
!
!
!
      SUBROUTINE SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
      DIMENSION TSTR(NGM), GSTR(NGM)
      DIMENSION GDOT(NGM), GDDOT(NGM)
      DIMENSION TSPLINE(NGM)
!
      DIF = 0.D0
      HELPAUX = 0.D0
!
! PRIMERA APROXIMACION A TSPLINE(J)
!	DO J=1,NG

!	  IF(GDOT(J).EQ.0.AND.GDDOT(J).EQ.0) TSPLINE(J) = 10000000.

!	  IF(GDOT(J).EQ.0.AND.GDDOT(J).EQ.0) GO TO 127
!***********************************************************************

!	    if(GDOT(J).EQ.0.AND.GDDOT(J).LT.0)then	!Fran 17/06/09	 !

!	      HELPAUX = -0.75D0/GDDOT(J)
!	      TSPLINE(J) = exp(HELPAUX)*TSTR(J)
!	      go to 1270

!	    else
!***********************************************************************
!	      HELPAUX=(0.25-1.)/GDOT(J) + 1.
!	      TSPLINE(J)=HELPAUX*TSTR(J)

!1270	    end if

!127	END DO
!
! CALCULO FINO DE TSPLINE
!
!	DO J=1,NG
!	  IF(GDDOT(J).NE.0) THEN
!	    CONTINUE
!	  ELSE
!***********************************************************************
!	    IF(GDOT(J).LT.0) THEN		!Fran 17/06/09	 !
!
!***********************************************************************

!15	      HELPAUX=(0.25 - 1. -GDDOT(J)*DLOG(TSPLINE(J)/TSTR(J)))
!	      HELPAUX=HELPAUX/GDOT(J) + 1.
!	      HELPAUX=HELPAUX*TSTR(J)
!	      DIF=DABS(HELPAUX-TSPLINE(J))
!	      IF(DIF.LE.1.D-7) TSPLINE(J)=abs(HELPAUX)
!	      IF(DIF.GT.1.D-7) TSPLINE(J)=abs(HELPAUX)
!	      IF(DIF.GT.1.D-7) GO TO 15
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!	    end if
!************************************************************
!	  END IF
!	END DO

!-----------------------------------------------------------------------
!
!	Halla temperatura de spline con Newton

      do I = 1,NG

        IF(GDOT(I).LT.-1.D-4)then	!GT.-1.D-4) GOTO 18

          TR=1.D0-1.D0/GDOT(I)
          IT=0
          abs_AX = 1.
          do while(abs_AX.GT.1.D-15)

!	      Función a ser igualada a 0 = AX(Tspl)
16          AX=1.D0+Gdot(I)*(TR-1.D0)+Gddot(I)*DLOG(TR)-1.D0/4.D0
!	      Derivada de dicha función: d(AX)/d(Tspl)
            DAX=GDOT(I)+GDDOT(I)/TR
!	      Cociente de Newton:
            DTR=-AX/DAX
            TR=TR+DTR
            IT=IT+1
!	    IF(TR.LT.0.D0.OR.IT.GT.10) GOTO 18
            abs_AX=DABS(AX)
!	      IF(AX.GT.1.D-15) GOTO 16
            if(TR.GE.0D0.AND.it.LE.10)then

              TSPLINE(I)=TR*TSTR(I)
!	        EPX(I)=-4.D0*(GDOT(I)*TR+GDDOT(I))
!	        GOTO 15

            else

!	        Si Ts está dando < 0 o no converge, Ts = 1000000.
              TSPLINE(I)=1.D7
              exit

            endif

          enddo

        else
!
!	   Si g' > 0, no puede haber una Ts > 0 (o difícilmente la haya)
18       TSPLINE(I)=1.D3

        endif
!15	  WRITE(IOUT,106)ID(I),TSTR(I),Q(I),GSTR(I),GDOT(I),GDDOT(I),(NY(J,I),J=1,NC)

      enddo

!
      ENDsubroutine
!//////////////////////////////////////////////////////////////
        SUBROUTINE CALGGG(NG,TEMP,TSTR,GSTR,GDOT,GDDOT,TSPLINE,GGG)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION TSTR(NGM), GSTR(NGM)
        DIMENSION GDOT(NGM), GDDOT(NGM)
        DIMENSION TSPLINE(NGM)
        DIMENSION GGG(NGM,NGM)
!
        ADIM = 0.D0
        EXPONENTE = 0.0D0
!
! GGG(I,J) = PARAMETRO DE ENERGIA ATRACTIVA (DE DISPERSION) ENTRE
!             EL GRUPO I Y EL GRUPO J
!
        DO J=1,NG
!
          ADIM=TEMP/TSTR(J)
!	    write(*,*)'T/Tstr[',J,'] =',ADIM
          IF(TEMP.LE.TSPLINE(J)) THEN

      	  GGG(J,J)=  ADIM - 1.D0
            GGG(J,J)= GGG(J,J)*GDOT(J)
!		  write(*,*)'gp[',J,J,'](T/T* - 1) =',GGG(J,J)
            GGG(J,J)= GGG(J,J) + 1.D0
            GGG(J,J)= GGG(J,J) + GDDOT(J)*DLOG(ADIM)
            GGG(J,J)= GGG(J,J)*GSTR(J)
!		  write(*,*)'g[',J,J,'] = ',GGG(J,J)

          ELSE

!	      write(*,*)'T spline =',J,TSPLINE(J)
            EXPONENTE= GDOT(J)*TSPLINE(J)/TSTR(J)
            EXPONENTE= EXPONENTE + GDDOT(J)
            EXPONENTE= (-4.D0)*EXPONENTE
            GGG(J,J)= 0.25D0*GSTR(J)
            GGG(J,J)= GGG(J,J)/(TEMP/TSPLINE(J))&
     &                **EXPONENTE
!	      write(*,*)'g[',J,J,'] = ',GGG(J,J)

          END IF

        END DO
!
        END
!////////////////////////////////////////////////////////////////
        SUBROUTINE CALCKIJ(NG,RKIJ,RKTIJ,TSTR,VERDKIJ,TEMP)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION TSTR(NGM)
        DIMENSION RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
        DIMENSION VERDKIJ(NGM,NGM)
!
        HELPVERD = 0.D0
!
        DO I=1,NG
                DO J=1,NG
                        HELPVERD=TSTR(I)+TSTR(J)
                        HELPVERD=HELPVERD/2
!	write(*,*)TEMP
                        HELPVERD=DLOG(TEMP/HELPVERD)
                        HELPVERD=RKTIJ(I,J)*HELPVERD + 1.D0
                        VERDKIJ(I,J)=RKIJ(I,J)*HELPVERD
                END DO
        END DO
!
        END
!////////////////////////////////////////////////////////////////
        SUBROUTINE CALGGGIJ(NG,VERDKIJ,GGG)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION GGG(NGM,NGM)
        DIMENSION VERDKIJ(NGM,NGM)
!
!
!
        DO I=1,NG
!
                DO J=1,NG
!					   write(*,*)TEMP
!					   write(*,*)'g(',I,I,') =',GGG(I,I)
!					   write(*,*)'g(',j,j,') =',GGG(j,j)
                         GGG(I,J)=VERDKIJ(I,J)*DSQRT(GGG(I,I)*GGG(J,J))
                END DO
        END DO
!
        END
!////////////////////////////////////////////////////////////////
        SUBROUTINE DELTADISP(NG,GGG,DLTG)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION GGG(NGM,NGM)
        DIMENSION DLTG(NGM,NGM)
!
! DLTG(I,J) = DELTAgij  (NO ES SIMETRICO)
!
        DO J=1,NG
                DO I=1,NG
                        DLTG(J,I)=GGG(J,I)-GGG(I,I)
                END DO
        END DO
!
        END
!///////////////////////////////////////////////////////////
        SUBROUTINE TAUCALC(NG,TEMP,V,ALFAIJ,QTOT,DLTG,TAU)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
!
        DIMENSION ALFAIJ(NGM,NGM)
        DIMENSION DLTG(NGM,NGM)
        DIMENSION TAU(NGM,NGM)
!
        HELPTAU = 0.D0
!
! CTE. UNIVERSAL DE LOS GASES
!
        RUNIV=82.05D0 !atm cm3/(K mol)
!
        DO J=1,NG
                DO I=1,NG
                        HELPTAU=ALFAIJ(J,I)*DLTG(J,I)*QTOT

                        HELPTAU= HELPTAU/(RUNIV*TEMP*V)
!					  write(*,*)HELPTAU
                        TAU(J,I)= DEXP(HELPTAU)
                END DO
        END DO
!
        END
!/////////////////////////////////////////////////////////
        SUBROUTINE AHELMRESDISP(NC,NG,RN,V,TEMP,NU,TITA,TAU,QTOT,&
     &GGG,QAREA,ARESDISP)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50, NGM = 20
        DIMENSION RN(NC)
!
        DIMENSION NU(NGM,NC), QAREA(NGM)
        DIMENSION TITA(NGM), TAU(NGM,NGM)
        DIMENSION GGG(NGM,NGM)
        DIMENSION AUXIL2(NCM,NCM)
        DIMENSION SUM1(NGM),SUM2(NGM),SUM3(NCM)
!
! LIMPIEZA DE VARIABLES
!
        SUMAT = 0.D0
!
        DO I=1,NC+1
                 SUM3(I)=0.D0
!
                DO J=1,NC+1
                        AUXIL2(I,J)=0.D0
                END DO
        END DO
!
!
        DO I=1,NG
                SUM1(I)=0.D0
                SUM2(I)=0.D0
        END DO
!
! CTE. UNIVERSAL DE LOS GASES
!
        RUNIV=82.05D0 !atm cm3/(K mol)
!
        ZETA=10.D0 !NUMERO DE COORDINACION
!
! CALCULO DE LOS SUM1(J)
!
        DO J=1,NG
                SUM1(J)=0.D0
                DO L=1,NG
                        SUM1(J) = SUM1(J) + TITA(L)*TAU(L,J)
                END DO
        END DO
!
!
!
! AUXIL2 --------------------------------------
!
        DO J=1,NG
                DO K=1,NG
                        AUXIL2(K,J)=TITA(K)*(GGG(K,J)*QTOT*TAU(K,J))
                        AUXIL2(K,J)=AUXIL2(K,J)/(RUNIV*TEMP*V)
                END DO
        END DO
!
! CALCULO DE LOS SUM2(J)
!
        DO J=1,NG
                SUM2(J)=0.D0
                 DO K=1,NG
                        SUM2(J) = SUM2(J) + AUXIL2(K,J)
                 END DO
        END DO
!
! CALCULO DE LOS SUM3(I)
!
        DO I=1,NC
               SUM3(I)=0.D0
                DO J=1,NG
                  SUM3(I)= SUM3(I) + NU(J,I)*QAREA(J)*SUM2(J)/SUM1(J)
                END DO
        END DO
!
! SUMATORIA EN  EL NUMERO DE COMPONENTES
!
        SUMAT=0.D0
        DO I=1,NC
                 SUMAT= SUMAT + RN(I)*SUM3(I)
        END DO
!
! CALCULO DE ARESDISP(energ. resid. de dispersion)
        ARESDISP= (-1.D0)*(ZETA/2D0)*SUMAT
!C-----------------------------------------------------------------------
!
!       ARESDISP= CONTRIBUCION DE FZAS. DE DISPERSION A LA ENERGIA
!                RESIDUAL DE HELMHOLTZ =
!
!                      R
!                 = ( A  / RT )         (las unidades son "moles")
!                              disp.
!
!----------------------------------------------------------------
!
      return
        ENDsubroutine
!------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!
      SUBROUTINE CALCDELTASS(NGA,TEMP,MASSOC,ENASS,RKASS,DELTASS)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION DELTASS(6,12,6,12)
!
!	ENASS(I,J,K,L) = PARAMETRO DE ENERGIA DE ASOCIACION ENTRE
!	             EL SITIO I DEL grupo J, Y EL SITIO K
!                  DEL grupo L (EN GRADOS KELVIN)
!
! DELTASS(I,J,K,L) = FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL grupo J, Y EL SITIO K
!                  DEL grupo L (ADIMENSIONAL). NUESTRA DEFINICION
!                  NO COINCIDE CON LA DE RADOSZ.
!     Con esta definicion EL SISTEMA DE ECUACIONES (CUYAS VARIABLES SON
!       LAS FRACCIONES NO ASOCIADAS DE CADA SITIO) QUEDA EN FUNCION DE:
!
!            - FRACCIONES DE SITIOS NO ASOCIADOS y,
!            - FUERZA DE ASOCIACION (ADIMENSIONAL)
      DeltAss(:NGA,:2,:NGA,:2) = 0D0
      DO J=1,NGA

        IF(MASSOC(J).NE.0)then	!EQ.0) GO TO 201

          DO I=1,NGA

            IF(MASSOC(I).NE.0)then	!EQ.0) GO TO 101

              DO L=1,MASSOC(J)

                DO K=1,MASSOC(I)

                  IF(ENASS(K,I,L,J).NE.0)then	!EQ.0) THEN
!                                         CONTINUE
!                                    ELSE
!	               AUXO = 0.D0
!	               AUXO = DEXP(ENASS(K,I,L,J)/TEMP) - 1.D0
!	               AUXO = AUXO * RKASS(K,I,L,J)		!* RNGRUP(J) / V
!	               DELTASS(K,I,L,J) = AUXO
                     DeltAss(k,i,l,j) = rKAss(k,i,l,j)*&
     &                                 (dexp(EnAss(k,i,l,j)/Temp) - 1D0)

                  END IF

                END DO

              END DO

            endif

          END DO

        endif

      END DO

      return
      ENDsubroutine
!--------------------------------------------------------------------------
!
!
!	Antigua subrutina de cálculo de la fracción no asociada. Ahora es
!	puente a la subrutina OptiNewton, utilizada en el paquete GCA-2,
!	presente en GC3, GCDist, GCEnvel, etc.
!	Se eliminó el código antiguo. Si se quiere recuperar, buscarlo en
!	la versión 3.0.0 (A-gc-eos3.for)
!
!	Francisco, durante versión 3.1.0.
!
      SUBROUTINE FRACNOAS(NGA,MASSOC,DELTASS,XNOASS, nc, ONUASOC,V,&
     &                    RNGRUP,iErrAssoc)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION MASSOC(12),sumsi(12)
      DIMENSION DELTASS(6,12,6,12)
      DIMENSION XNOASS(6,12)
      DIMENSION ONUASOC(12,NC)
      DIMENSION RNGRUP(12)
      dimension dif(12,12)
      dimension xnoas_old(6,12)
!
!**************************ALFONSINA******************************************
!	GENERAL
!	Iniciacio
      if(NGA.NE.0) then

        xnoass=1.0d0
        del=1.d0
!	  Iteraciones con tolerancia de 10-9
        do while (del>1.0d-9)

          xnoas_old=xnoass
          do i=1,NGA

            do j=1,2

              sum1=0.D0
              do k=1,NGA

                sum2=0.D0
                do l=1,2

                  sum2=sum2+xnoas_old(l,k)*deltass(l,k,j,i)

                end do
                sum1=sum1+rngrup(k)*sum2

              end do
              xnoass(j,i)=1.0d0/(1.0d0+sum1/v)
              dif(j,i)=dabs((xnoass(j,i)-xnoas_old(j,i))/xnoass(j,i))

            end do

          end do
          del=maxval(dif)

        end do

      end if

!	Fin del Cálculo de la  ffracción no asociada
      return
      endsubroutine

!CCCCCCCCCCCCCCCCC ALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC4(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP&
     &,ierrassoc)

    !alcohol+AMINA PRIMARIA

! CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
! QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION

        IMPLICIT REAL*8 (A-H,O-Z)
!        implicit none
        double precision DELTASS(6,12,6,12)
        double precision XNOASS(6,12)
        double precision RNGRUP(12)
        double precision V,ROOT1(4)
        integer JJ,IDAUT1,IDAUT2, ierrassoc
        DIMENSION ZC(3)
        DOUBLE precision :: qa0,qa1,qa2,qa3,qa4,zq(4)
        DOUBLE PRECISION DELCROSS_2

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
!   UN SITIO (ESTER)

! IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
!   DOS SITIOS (ALCOHOL)
!
! XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
! XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
! XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
!

!	PRIMEIRO CALCULAR XB2=X(2,2)
!	qa4=DELTASS(1,1,1,1)*(RNGRUP(2)/V)**2*DELTASS(1,2,2,2)**2-
!     @2D0*(RNGRUP(2)/V)**2*DELTASS(1,1,2,2)**2*DELTASS(1,2,2,2)



      DELTASS(1,1,2,2)=(DELTASS(1,1,2,1)*DELTASS(1,2,2,2))**0.5D0
      DELCROSS_2=(DELTASS(1,1,2,1)*DELTASS(1,2,2,2))

!	a3=DELTASS(1,1,2,1)*(RNGRUP(1)/V)*(DELTASS(1,1,2,2)-
!     @DELTASS(1,1,2,1))
!	a2=DELTASS(1,1,2,2)-DELTASS(1,1,2,1)-
!     @((RNGRUP(2)/V)+(RNGRUP(1)/V))*DELTASS(1,1,2,1)*
!     @DELTASS(1,1,2,2)
!	a1=-2D0*DELTASS(1,1,2,2)+DELTASS(1,1,2,1)
!	a0=DELTASS(1,1,2,2)


!		call cubicequation(a3,a2,a1,a0,zc)


      qa4=(RNGRUP(1)/V)**2D0*DELTASS(1,1,2,1)*&
     &(DELTASS(1,1,2,1)*DELTASS(1,2,2,2)-DELCROSS_2)
      qa3=(RNGRUP(1)/V)*(2d0*DELTASS(1,1,2,1)*DELTASS(1,2,2,2)-&
     &DELTASS(1,1,2,2)*DELTASS(1,1,2,1)-DELTASS(1,1,2,2)**2D0)
      qa2=DELTASS(1,2,2,2)-2d0*DELTASS(1,1,2,1)*DELTASS(1,2,2,2)*&
     &(RNGRUP(1)/V)+(RNGRUP(1)/V)*DELTASS(1,1,2,2)**2D0-DELTASS(1,1,2,2)&
     &-(RNGRUP(2)/V)*DELTASS(1,1,2,2)**2D0
      qa1=-2D0*DELTASS(1,2,2,2)+DELTASS(1,1,2,2)
      qa0=DELTASS(1,2,2,2)

!		call quarticequation(qa0,qa1,qa2,qa3,qa4,zq)


      	ierrassoc=1

!		DO JJ=1,3
!			If (zc(JJ)<1d0.and.zc(JJ)>0d0)then  !numero real entre 0 e 1
!					root1(jj)=zc(jj)
!					XNOASS(1,1)=root1(jj)

      	DO JJ=1,4
      		If (zq(JJ)<1d0.and.zq(JJ)>0d0)then  !numero real entre 0 e 1
      				root1(jj)=zq(jj)
      				XNOASS(1,1)=root1(jj)

      				If(rngrup(2)==0d0) then
      					XNOASS(2,2)=1.D0/(1d0+RNGRUP(1)/V&
     &				*XNOASS(1,1)*DELTASS(1,1,2,2))

      				else
      					XNOASS(2,2)=(1.D0-XNOASS(1,1)-&
     &				(RNGRUP(1)/V)*DELTASS(1,1,2,1)*XNOASS(1,1)&
     &				**2)/(DELTASS(1,1,2,2)*(RNGRUP(2)/V)&
     &				*XNOASS(1,1))
      				end if

      				If (XNOASS(1,1)<=1d0.and.XNOASS(1,1)>=0d0)then
      					XNOASS(1,2) = XNOASS(2,2)
      					XNOASS(2,1) = XNOASS(1,1)
      					ierrassoc=0
      					goto 1132
      				end if
      			endif
      	END DO


1132  	If (ierrassoc==1) then
      		write(6,*)"No roots found"
      		XNOASS(1,1)=1d0
      		XNOASS(2,1)=1d0
      		XNOASS(1,2)=1d0
      		XNOASS(2,2)=1d0

      	else

      	end if

        END
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC3(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP) !!ester+acid
!
! CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
! QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
!
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION ZC(3)
        DIMENSION RNGRUP(12)
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
!   UN SITIO (ESTER)

! IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
!   DOS SITIOS (ALCOHOL)
!
! XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
! XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
! XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
!

!	PRIMEIRO CALCULAR XB2=X(2,2)

      a3=DELTASS(1,1,1,1)*DELTASS(1,1,1,2)*(RNGRUP(1)/V)**2D0
      a2=(RNGRUP(1)/V)*(DELTASS(1,1,1,2)+DELTASS(1,1,1,1))
      a1=1D0+DELTASS(1,1,1,2)*(RNGRUP(2)/V-RNGRUP(1)/V)
      a0=-1d0

!	If (a3.gt.1e-13) then

!		call cubicequation(a3,a2,a1,a0,zc)

      	do jj=1, 3
      		If(zc(jj)>0d0.and.zc(jj)<=1d0)then
      			root1=zc(jj)
      		end if
      	end do

      	If (root1.ne.0d0) then
      		XNOASS(1,1)=root1
      		XNOASS(1,2) = 1.D0/(1D0+DELTASS(1,1,1,2)*(RNGRUP(1)/V)&
     &		            *XNOASS(1,1))
      	else
      		write(6,*)"No roots found"
      		XNOASS(1,1)=1d0
      		XNOASS(1,2)=1d0
      	end if
!	end if

!	write(2,*)"XNOASS(1,1)=",	XNOASS(1,1)
!	write(2,*)"XNOASS(1,2)=", XNOASS(1,2)
        END
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC2(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP& !alcohol+acid
     &,ierrassoc)

!
! CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
! QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
!
!        IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
        double precision DELTASS(6,12,6,12)
        double precision XNOASS(6,12)
        double precision RNGRUP(12)
        double precision V,ROOT1(4)
        integer JJ,IDAUT1,IDAUT2, ierrassoc
        DOUBLE precision :: qa0,qa1,qa2,qa3,qa4,zq(4)
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
!   UN SITIO (ESTER)

! IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
!   DOS SITIOS (ALCOHOL)
!
! XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
! XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
! XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
!

!	PRIMEIRO CALCULAR XB2=X(2,2)
      qa4=DELTASS(1,1,1,1)*(RNGRUP(2)/V)**2*DELTASS(1,2,2,2)**2-&
     &2D0*(RNGRUP(2)/V)**2*DELTASS(1,1,2,2)**2*DELTASS(1,2,2,2)
      qa3=-DELTASS(1,1,2,2)*(RNGRUP(2)/V)*DELTASS(1,2,2,2)-2D0*&
     &(RNGRUP(2)/V)*DELTASS(1,1,2,2)**2+2D0*DELTASS(1,1,1,1)*&
     &(RNGRUP(2)/V)*DELTASS(1,2,2,2)
      qa2=DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-2D0*DELTASS(1,1,1,1)*&
     &(RNGRUP(2)/V)*DELTASS(1,2,2,2)-(RNGRUP(1)/V)*DELTASS(1,1,2,2)**2+&
     &2D0*(RNGRUP(2)/V)*DELTASS(1,1,2,2)**2
      qa1=-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2)
      qa0=DELTASS(1,1,1,1)

!		call quarticequation(qa0,qa1,qa2,qa3,qa4,zq)

      	ierrassoc=1

      	DO JJ=1,4
      		If (zq(JJ)<1d0.and.zq(JJ)>0d0)then  !numero real entre 0 e 1
      				root1(jj)=zq(jj)
      				XNOASS(2,2)=root1(jj)

      				If(rngrup(1)==0d0) then
      					XNOASS(1,1)=1.D0/(1d0+2d0*RNGRUP(2)/V&
     &				*XNOASS(2,2)*DELTASS(1,1,2,2))
      				else
      					XNOASS(1,1)=(1.D0-XNOASS(2,2)-&
     &				(RNGRUP(2)/V)*DELTASS(1,2,2,2)*XNOASS(2,2)&
     &				**2)/(DELTASS(1,1,2,2)*(RNGRUP(1)/V)&
     &				*XNOASS(2,2))
      				end if

      				If (XNOASS(1,1)<=1d0.and.XNOASS(1,1)>=0d0)then
      					XNOASS(1,2) = XNOASS(2,2)
      					ierrassoc=0
      					goto 1132
      				end if
      			endif
      	END DO


1132  	If (ierrassoc==1) then
      		write(6,*)"No roots found"
      		XNOASS(1,1)=1d0
      		XNOASS(2,2)=1d0
      		XNOASS(1,2)=1d0
      	else
!			write(6,*)"Roots found"
!	write(6,*)"XNOASS(1,1)=",	XNOASS(1,1)
!	write(6,*)"XNOASS(2,2)=", XNOASS(2,2)
!	write(6,*)"XNOASS(1,2)=", XNOASS(1,2)

      	end if

!	write(2,*)"XNOASS(1,1)=",	XNOASS(1,1)
!	write(2,*)"XNOASS(2,2)=", XNOASS(2,2)
!	write(2,*)"XNOASS(1,2)=", XNOASS(1,2)


!
        END


!////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP)!ester+alcohol
!
! CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
! QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
!
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION ZC(3)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
!   UN SITIO (ESTER)

! IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
!   DOS SITIOS (ALCOHOL)
!
! XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
! XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
! XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
!

!	PRIMEIRO CALCULAR XB2=X(2,2)

      a3=DELTASS(1,1,2,2)*DELTASS(1,2,2,2)*(RNGRUP(2)/V)**2D0
      a2=DELTASS(2,2,1,1)*DELTASS(1,2,2,2)*(RNGRUP(2)/V)*(RNGRUP(1)/V)+&
     &DELTASS(1,1,2,2)*(RNGRUP(2)/V)+DELTASS(1,2,2,2)*(RNGRUP(2)/V)
      a1=DELTASS(2,2,1,1)*(RNGRUP(1)/V)-DELTASS(1,1,2,2)*(RNGRUP(2)/V)+1
      a0=-1d0

!	If (a3.gt.1e-13) then

!		call cubicequation(a3,a2,a1,a0,zc)

      	do jj=1, 3
      		If(zc(jj)>0d0.and.zc(jj)<=1d0)then
      			root1=zc(jj)
      		end if
      	end do

      	If (root1.ne.0d0) then
      		XNOASS(2,2)=root1
      		XNOASS(1,1) = 1.D0/(1D0+DELTASS(1,1,2,2)*(RNGRUP(2)/V)&
     &					*XNOASS(2,2))
      		XNOASS(1,2) = 1.D0/(1D0+DELTASS(1,2,2,2)*(RNGRUP(2)/V)&
     &		            *XNOASS(2,2))
      	else
      		write(6,*)"No roots found"
      		XNOASS(1,1)=1d0
      		XNOASS(2,2)=1d0
      		XNOASS(1,2)=1d0
      	end if
!	end if

!	write(2,*)"XNOASS(1,1)=",	XNOASS(1,1)
!	write(2,*)"XNOASS(2,2)=", XNOASS(2,2)
!	write(2,*)"XNOASS(1,2)=", XNOASS(1,2)
        END
!////////////////////////////////////////////////////////////////////////

! A ESTA SUBRUTINA ES NECESARIO AMPLIARLA PARA INCLUIR LA POSIBILIDAD DE
! AUTOASOCIACION CON UN NUMERO DE SITIOS DISTINTO DE DOS, Y LA POSIBILIDAD
! DE ASOCIACION CRUZADA.
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE AUTOASOC(IDAUT,DELTASS,XNOASS,V,RNGRUP)
!
! CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
! QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
!
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION RNGRUP(12)
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
!   DOS SITIOS
!
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo
!                                                       IDAUT
!
        IF(DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

        AUXO2 = 1D0 +  4D0 *DELTASS(1,IDAUT,2,IDAUT)* RNGRUP(IDAUT)/V

!
      	If ((AUXO2.gt.0d0).and.(RNGRUP(IDAUT).gt.0d0)) then
        AUXO2 =  (-1D0) + DSQRT(AUXO2)
!
        AUXO2 = AUXO2 / (2D0 * DELTASS(1,IDAUT,2,IDAUT)*RNGRUP(IDAUT)/V)
        XNOASS(1,IDAUT) = AUXO2
      	else
      	 XNOASS(1,IDAUT) = 1.
      	end if
         ELSE
!
!
                XNOASS(1,IDAUT) = 1.
!
!       ----------------
!
        END IF
!
!       ----------------
!

          IF(DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

      		XNOASS(2,IDAUT) = XNOASS(1,IDAUT)

      	END IF
!	       write(2,*)"XNOASS(1,",IDAUT,")=", XNOASS(1,IDAUT)
!            write(2,*)"XNOASS(2,",IDAUT,")=", XNOASS(2,IDAUT)

!
        END
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE autoasocumsitio(IDAUT,DELTASS,XNOASS,V,RNGRUP)
!
! CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
! QUE ASOCIA, QUE TIENE uno SITIO DE ASOCIACION
!
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION RNGRUP(12)
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
!   uno SITIO
!
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo
!                                                       IDAUT
!
        IF(DELTASS(1,IDAUT,1,IDAUT).NE.0.and.rngrup(1).gt.0d0) THEN

        AUXO2 = 1D0 + 4D0 * DELTASS(1,IDAUT,1,IDAUT)* RNGRUP(IDAUT)/V
!
        AUXO2 =(-1D0) + DSQRT(AUXO2)
!
        AUXO2 = AUXO2 / (2D0* DELTASS(1,IDAUT,1,IDAUT)*RNGRUP(IDAUT)/V)
        XNOASS(1,IDAUT) = AUXO2
      	ELSE
                XNOASS(1,IDAUT) = 1D0
!
!       -----

      	endif
!
!       write(2,*)"XNOASS(1,",IDAUT,")=", XNOASS(1,IDAUT)
         END
!/////////////////////////////////////////////////////////////////////////


        SUBROUTINE AHELMRESASS(NGA,MASSOC,RNGRUP,XNOASS,ARESASS,&
     &FASS,FAASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION RNGRUP(12)
        DIMENSION XNOASS(6,12)
        DIMENSION FASS(12)
        DIMENSION FAASS(6,12)
!
! CALCULO DE LOS FAASS(K,I)
!
        DO I=1,NGA
                IF(MASSOC(I).EQ.0) GO TO 881
                DO K=1,MASSOC(I)
!
                        AUXO3 = 0.D0
      				  If(XNOASS(K,I).gt.1d-13)then
      					AUXO3 = DLOG( XNOASS(K,I) )
      					AUXO3 = AUXO3 - XNOASS(K,I)/2.D0
      					FAASS(K,I) = AUXO3
      				  end if
!
                END DO
881             CONTINUE
        END DO
!
!
! CALCULO DE LAS SUMATORIAS SOBRE LOS SITIOS DE ASOCIACION PARA CADA
!  GRUPO
!
! FASS(I)= {SUMATORIA SOBRE SITIOS DE [ FAASS(K,I) ] }  +
!
!                                + 0.5 MASSOC(I)
!
! SI FASS(I) SE MULTIPLICA POR EL NRO. DE MOLES DEL grupo(I)
!
!               SE OBTIENE LA CONTRIBUCION DE ESE grupo a ARESASS
!
        DO I=1,NGA
!          -----------------------------------------------
           IF(MASSOC(I).NE.0) THEN
!
               AUXO4=0.D0
!
                DO K=1,MASSOC(I)
                  AUXO4 = AUXO4 + FAASS(K,I)
                END DO
!
               FASS(I)= AUXO4 + 0.5D0 * MASSOC(I)
!
           END IF
!          -----------------------------------------------
!
        END DO
!
! CALCULO FINAL DE LA ENERGIA RESIDUAL DE HELMHOLTZ POR ASOCIACION
!
        SUME=0.D0
!
        DO I=1,NGA
          SUME = RNGRUP(I) * FASS(I) + SUME
        END DO
!
        ARESASS = SUME
!-----------------------------------------------------------------------
!
!       ARESASS = CONTRIBUCION DE LA ASOCIACION A LA ENERGIA
!                RESIDUAL DE HELMHOLTZ =
!
!                      R
!                 = ( A  / RT )         (las unidades son "moles")
!                              ass
!
!----------------------------------------------------------------
!
        END
!///////////////////////////////////////////////////////////////////////
        SUBROUTINE AHELMRESTOT(ARESFV,ARESDISP,ARESASS,ARESTOT)
        IMPLICIT REAL*8 (A-H,O-Z)
!
! CALCULO DE LA ENERGIA RESIDUAL TOTAL
!
        ARESTOT = ARESFV + ARESDISP + ARESASS
!
        END
!//////////////////////////////////////////////////////////////////
        SUBROUTINE DERIVFV(NC,V,D,Y,DL1N,DL2N,DL3N,DYN)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
        DIMENSION D(NCM),DL1N(NCM),DL2N(NCM),DL3N(NCM),DYN(NCM)
!
! CALCULA DERIVADAS SEGUN DETALLE
!
! DL1N(I) = DERIVADA DE RLAMBDA1 CON RESPECTO AL NUMERO DE MOLES DEL
!                                       COMPONENTE I
!
! DL2N(I) = DERIVADA DE RLAMBDA2 CON RESPECTO AL NUMERO DE MOLES DEL
!                                       COMPONENTE I
!
!
! DL3N(I) = DERIVADA DE RLAMBDA3 CON RESPECTO AL NUMERO DE MOLES DEL
!                                       COMPONENTE I
!
! DYN(I)  = DERIVADA DE LA VARIABLE Y CON RESPECTO AL NRO. DE MOLES
!                                       DEL COMPONENTE I
!
!
        DO I=1,NC
                DL1N(I) = D(I)
                DL2N(I) = D(I)**2
                DL3N(I) = D(I)**3
                DYN(I)  = Y**2 * 3.141592654D0 * D(I)**3 / (6.D0*V)
        END DO
!
        END
!
!//////////////////////////////////////////////////////////////////
        SUBROUTINE FUGAREP(NC,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,&
     &DL1N,DL2N,DL3N,DYN,FUGREP)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
        DIMENSION DL1N(NCM),DL2N(NCM),DL3N(NCM),DYN(NCM)
        DIMENSION FUGREP(NCM)
!
        BLOQUE1 =0.
        BLOQUE2 =0.
        BLOQUE3 =0.
        BLOQUE4 =0.
        HELPFUR =0.
!
! CALCULO DE LA CONTRIBUCION REPULSIVA A LOS COEFICIENTES DE FUGACIDAD
!
! FUGREP(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD REPULSIVO
!                DEL COMPONENTE I
!
        BLOQUE1= 3D0*(Y-1)*RLAMBDA1*RLAMBDA2/RLAMBDA3
!
        BLOQUE2= RLAMBDA2**3/RLAMBDA3**2
        BLOQUE3= BLOQUE2*( Y**2 - Y - DLOG(Y)  )
!
        BLOQUE4= BLOQUE2*( 2.D0*Y - 1D0 - 1D0/Y  )
!
!
!
        DO I=1,NC
!
         HELPFUR = DL1N(I)/RLAMBDA1+DL2N(I)/RLAMBDA2-DL3N(I)/RLAMBDA3
         HELPFUR = HELPFUR + DYN(I)/(Y-1D0)
         HELPFUR = HELPFUR*BLOQUE1
         HELPFUR = HELPFUR + BLOQUE3 *&
     &( 3.D0 * DL2N(I)/RLAMBDA2 - 2.D0 * DL3N(I)/RLAMBDA3)
         FUGREP(I)= HELPFUR + BLOQUE4*DYN(I) + DLOG(Y) +&
     &(RNTOT/Y)*DYN(I)
!
        END DO
!
        END
!/////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////
        SUBROUTINE HCALC(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,&
     &TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NCM = 50, NGM = 20
!
! NOTAR QUE EN EL CASO DE NU EL PRIMER ARGUMENTO CORRESPONDE A GRUPO
!   Y EL SEGUNDO A COMPONENTE
!
! EN EL CASO DE PS,H3, Y H7 : PRIMER ARG --- COMPONENTE
!                             SEGUNDO ARG--- GRUPO
!
        DIMENSION PS(NCM,NGM), RMS(NCM), H2(NGM)
        DIMENSION H3(NCM,NGM), H4(NGM), H5(NGM)
        DIMENSION H6(NGM), H7(NCM,NGM)
!
        DIMENSION NU(NGM,NC), QAREA(NGM)
        DIMENSION ALFAIJ(NGM,NGM)
        DIMENSION TITA(NGM), TAU(NGM,NGM)
        DIMENSION GGG(NGM,NGM)
        DIMENSION DLTG(NGM,NGM)
!
        AUXH2 = 0.D0
        AUXH3 = 0.D0
        AUXH4 = 0.D0
        AUXH5 = 0.D0
        AUXH6 = 0.D0
        AUXH7 = 0.D0
        AUXRMS = 0.D0
!
        RUNIV=82.05D0 !atm cm3/(K mol)
!
!--------CALCULO DE PS
!
        DO I=1,NC
                DO J=1,NG
                        PS(I,J)= NU(J,I) * QAREA(J)
                END DO
        END DO
!
!--------CALCULO DE RMS
!
        DO I=1,NC
                AUXRMS=0.D0
                 DO J=1,NG
                        AUXRMS= PS(I,J) + AUXRMS
                 END DO
                RMS(I)=AUXRMS
        END DO
!
!--------CALCULO DE H2
!
!
        DO J=1,NG
                AUXH2=0.D0
                 DO K=1,NG
!
          AUXH2 = TITA(K)*TAU(K,J)*GGG(K,J)*QTOT/(V*RUNIV*TEMP)+ AUXH2
!
                 END DO
                H2(J)= AUXH2
        END DO
!
!--------CALCULO DE H3
!
!
        DO I=1,NC
                DO J=1,NG
!
                      AUXH3 =0.D0
                        DO K=1,NG
!
         AUXH3 = PS(I,K)*TAU(K,J)*GGG(K,J)*QTOT/(V*RUNIV*TEMP) + AUXH3
!
                        END DO
                      H3(I,J) = AUXH3
!
                END DO
!
        END DO
!
!
!--------CALCULO DE H4
!
!
        DO J=1,NG
!
               AUXH4 = 0.D0
!
                DO K=1,NG
!
                  AUXH4 = AUXH4 + TITA(K) * TAU(K,J)
!
                END DO
!
               H4(J)= AUXH4
!
        END DO
!
!
!--------CALCULO DE H5
!
!
        DO J=1,NG
!
             AUXH5 = 0.D0
!
                DO K=1,NG
!
                 AUXH5 = AUXH5 +  TITA(K)*TAU(K,J)*(GGG(K,J)*&
!
     &QTOT/(V*RUNIV*TEMP))*ALFAIJ(K,J)*DLTG(K,J)*QTOT/(V*RUNIV*TEMP)
!
                END DO
!
             H5(J)= AUXH5
!
        END DO
!
!
!
!--------CALCULO DE H6
!
!(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,
!     @TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
!
!
        DO J=1,NG
!
             AUXH6 = 0.D0
!
                DO K=1,NG
!
                AUXH6 = AUXH6 + TITA(K)*TAU(K,J)*ALFAIJ(K,J)&
     &*DLTG(K,J)*QTOT/(V*RUNIV*TEMP)
!
                END DO
!
             H6(J) = AUXH6
!
        END DO
!
!
!
!
!
!--------CALCULO DE H7
!
!(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,
!     @TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
!
        DO I=1,NC
                DO J=1,NG
!
                  AUXH7 = 0.D0
!
                        DO K=1,NG
!
                           AUXH7 = AUXH7 + PS(I,K)*TAU(K,J)
!
                        END DO
!
                  H7(I,J) = AUXH7
!
                END DO
!
        END DO
!
!
        END
!///////////////////////////////////////////////////////////
        SUBROUTINE FUGADISP(NC,NG,TITA,PS,RMS,H2,H3,H4,H5,H6,H7,FUGDISP)
!
! calculo de la contribucion dispersiva al coeficiente de fugacidad
!
! FUGDISP(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD DISPERSIVO
!                DEL COMPONENTE I
!
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50, NGM = 20
        DIMENSION PS(NCM,NGM), RMS(NCM), H2(NGM)
        DIMENSION H3(NCM,NGM), H4(NGM), H5(NGM)
        DIMENSION H6(NGM), H7(NCM,NGM)
!
        DIMENSION TITA(NGM)
!
        DIMENSION FUGDISP(NCM)
!
        SUM24 = 0.D0
        SUM2467 = 0.D0
        SUM345 = 0.D0
!
        ZETA=10.D0 !NUMERO DE COORDINACION
!
! CALCULO DE SUMATORIAS SOBRE EL NUMERO DE GRUPOS (p' c/componente)
!
        DO I=1,NC
!
          SUM24 = 0.D0
          SUM345 = 0.D0
          SUM2467 = 0.D0
!
          DO J=1,NG
!
            SUM24 = SUM24 + PS(I,J) * H2(J) / H4(J)
!
            SUM345 = SUM345 + TITA(J)*(H3(I,J) + RMS(I)*H5(J)) / H4(J)
!
                        SUM2467 = SUM2467 + TITA(J) * H2(J) *&
!
     &( H7(I,J) - H4(J) * RMS(I) + H6(J) * RMS(I) ) / H4(J)**2
!
          END DO
          FUGDISP(I) =(-1.D0) * (ZETA/2D0) *&
!
     &( SUM24 + SUM345 - SUM2467 )
!
        END DO
!
        ENDsubroutine FUGADISP
!//////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////
        SUBROUTINE DERIROSSN(NC,NGA,MASSOC,ONUASOC,RNGRUP,&
     &DROSSN,V)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NCM = 50
        DIMENSION MASSOC(12)
        DIMENSION ONUASOC(12,NC)
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION DROSSN(12,NCM)
!
! DELTASS(I,J,K,L) = FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L (ADIMENSIONAL). NUESTRA DEFINICION
!                  NO COINCIDE CON LA DE RADOSZ
!
! DDLTSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M  (ESTO ES LO QUE SE CALCULA EN ESTA
!                  SUBRUTINA)
!
!
        DO J=1,NGA
          IF(MASSOC(J) /= 0) then
!
            DO M=1,NC
!
              DROSSN(J,M) =  ONUASOC(J,M)/V
!
            END DO
201       endif
        END DO
!
!
        END
!///////////////////////////////////////////////////////////////////
      SUBROUTINE DERIXNOASS(NC,NGA,MASSOC,DELTASS,XNOASS,DROSSN,&
     &                      DXNOASSN, ONUASOC,RNGRUP,V)

!
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
      IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
      DIMENSION MASSOC(12), sumsi(12)
      DIMENSION DELTASS(6,12,6,12)
      DIMENSION XNOASS(6,12)
      DIMENSION DROSSN(12,NCM)
      DIMENSION DXNOASSN(6,12,NCM)
      DIMENSION ONUASOC(12,NC)
      DIMENSION RNGRUP(12)
      dimension psin(12)
      real*8  m_lambda(nga*2,nga*2), m_lambda1(nga*2,nga*2)
      integer order,w,s
!	real d1, psin
      dimension indx(NCM)
!
!	Para información sobre el código anterior a las modificaciones de Alfonsina,
!	ver subroutina "FracNoAsoc".
!	Francisco.

!CCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C GENERAL
!	cCalcula los elementos de la matriz deltapq
      w=0
      s=0
      do j=1,2
        do i=1,nga

          w=w+1
          do l=1, 2
            do k=1, nga

              s=s+1
              m_lambda(w,s)=rngrup(k)*deltass(l,k,j,i)*xnoass(j,i)**2/v
              if (w.eq.s) then

                m_lambda(w,s)=m_lambda(w,s)+ 1.0d0

              end if

            end do
          end do
          s=0

        end do
      end do
      order=nga*2
! Calcula los elementos de la matriz [Yp]
      do k=1, NGA

        w=0
        do j=1,2

          do i=1,nga

            sum10=0.0d0
            do l=1,2

              sum10=sum10+ (xnoass(l,k)*deltass(l,k,j,i))

            end do
            w=w+1
            psin(w)=-(xnoass(j,i)**2)*sum10/v

          end do

        end do
        N=order
        NP=order
        m_lambda1=m_lambda
        call  ludcmp(m_lambda1,N,NP,indx,d1)
        call lubksb(m_lambda1,N,NP,indx,psin)
! colectando las derivadas en su correspondiente subíndice
        w=0
        do j=1,2

          do i=1, nga

            w= w+1
            dxnoassn(j,i,k)=psin(w)

          end do

        end do

      end do
      do m=1,nc

        do i=1,nga

          do j=1,2

            do k=1,nga

              if (onuasoc(k,m).ne.0) then

                dxnoassn(j,i,m)=dxnoassn(j,i,k)*onuasoc(k,m)

              end if

            end do

          end do

        end do

      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      return
      ENDsubroutine

!////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOC3 (IDAUT1, IDAUT2, NC,&
     &DELTASS, XNOASS, DROSSN, DXNOASSN, RNGRUP, V)
        IMPLICIT REAL*8 (A-H, O-Z)
        integer, parameter :: NCM = 50
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION RNGRUP(12)

      DO M=1,NC

      	t1 = RNGRUP(1)/V
      	t2 = DROSSN(1,M)  !diff(t1,nq)
      	t3 = t2* DELTASS(1,1,1,2) !delcross
      	t4 = XNOASS(1,1)  !XA1(nq)
      	t5 = t4**2
      	t9 = t1*DELTASS(1,1,1,1)   !DELTA1
      	t15 = RNGRUP(2)/V   !RO2(nq)
      	t16 = DROSSN(2,M)  !diff(t15,nq)
      	t21 = t1*DELTASS(1,1,1,2)  !delcross
      	t26 = t1**2
      DXNOASSN(1,1,M)= (-t3*t5-t2*DELTASS(1,1,1,1)*t5-2D0*t9*t5*t4*&
     &DELTASS(1,1,1,2)*t2-t4*t16*DELTASS(1,1,1,2)+t3*t4)/(1D0+2D0*t21*t4&
     &+2D0*t9*t4+3D0*t26*DELTASS(1,1,1,1)*t5*DELTASS(1,1,1,2)+t15*&
     &DELTASS(1,1,1,2)-t21)



      DXNOASSN(1,2,M)=-1D0*XNOASS(1,2)*DELTASS(1,1,1,2)*&
     &(DROSSN(1,M)*XNOASS(1,1)+&
     & (RNGRUP(1)/V)*DXNOASSN(1,1,M))&
     &/(1D0+(RNGRUP(1)/V)*DELTASS(1,1,1,2)*XNOASS(1,1))

!    WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)
!	WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)


      END DO



      END
!////////////////////////////////////////////////////////////////////////
       SUBROUTINE DCROSSASOC4(IDAUT1,IDAUT2,NC,&
     &DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR




! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL
! NRO. DE MOLES DEL COMPUESTO M

!
      DO M=1,NC
      	t1 = RNGRUP(2)/V
      	t2 = DELTASS(1,1,1,1)*t1
      	t3 = DELTASS(1,2,2,2)**2
      	t4 = DROSSN(2,M)  !diff(t1,nq)
      	t8 = DELTASS(1,1,2,2)**2
      	t9 = t1*t8
      	t14 = XNOASS(2,2) !X(nq)
      	t15 = t14**2
      	t16 = t15**2
      	t21 = 2D0*t4*t8
      	t24 = 2D0*DELTASS(1,1,1,1)*t4*DELTASS(1,2,2,2)
      	t26 = t15*t14
      	t28 = RNGRUP(1)/V
      	t29 = DROSSN(1,M)  !diff(t28,nq)
      	t34 = t1**2
      	t45 = 2D0*t9
      	t47 = 2D0*t2*DELTASS(1,2,2,2)
      DXNOASSN(2,2,M) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-&
     &(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+&
     &t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*DELTASS(1,2,2,2&
     &))*&
     &t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*&
     &(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4&
     &5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


      DXNOASSN(1,1,M)=(-DROSSN(1,M)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2&
     &*XNOASS(1,1)*DROSSN(2,M)*DELTASS(1,1,2,2)*XNOASS(2,2)-2*&
     &XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSN(2,2,M))&
     &/(1+2*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2*(RNGRUP(2)/V)&
     &*DELTASS(1,1,2,2)*XNOASS(2,2))



      DXNOASSN(1,2,M) = DXNOASSN(2,2,M)
!      WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)


      END DO

!	DO M=1,2
!		WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)
!		WRITE(2,*)"DXNOASSN(2,2,",M,")=",DXNOASSN(2,2,M)
!	END DO

      END
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOC2(IDAUT1,IDAUT2,NC,&
     &DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR




! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL
! NRO. DE MOLES DEL COMPUESTO M

!
      DO M=1,NC
      	t1 = RNGRUP(2)/V
      	t2 = DELTASS(1,1,1,1)*t1
      	t3 = DELTASS(1,2,2,2)**2
      	t4 = DROSSN(2,M)  !diff(t1,nq)
      	t8 = DELTASS(1,1,2,2)**2
      	t9 = t1*t8
      	t14 = XNOASS(2,2) !X(nq)
      	t15 = t14**2
      	t16 = t15**2
      	t21 = 2D0*t4*t8
      	t24 = 2D0*DELTASS(1,1,1,1)*t4*DELTASS(1,2,2,2)
      	t26 = t15*t14
      	t28 = RNGRUP(1)/V
      	t29 = DROSSN(1,M)  !diff(t28,nq)
      	t34 = t1**2
      	t45 = 2D0*t9
      	t47 = 2D0*t2*DELTASS(1,2,2,2)
      DXNOASSN(2,2,M) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-&
     &(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+&
     &t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*DELTASS(1,2,2,2&
     &))*&
     &t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*&
     &(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4&
     &5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


      DXNOASSN(1,1,M)=(-DROSSN(1,M)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2D0&
     &*XNOASS(1,1)*DROSSN(2,M)*DELTASS(1,1,2,2)*XNOASS(2,2)-2*&
     &XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSN(2,2,M))&
     &/(1+2*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2*(RNGRUP(2)/V)&
     &*DELTASS(1,1,2,2)*XNOASS(2,2))



      DXNOASSN(1,2,M) = DXNOASSN(2,2,M)
!      WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)


      END DO

!	DO M=1,2
!		WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)
!		WRITE(2,*)"DXNOASSN(2,2,",M,")=",DXNOASSN(2,2,M)
!	END DO

      END



!////////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOC(IDAUT1,IDAUT2,NC,&
     &DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR




! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL
! NRO. DE MOLES DEL COMPUESTO M

!
      DO M=1,NC
      	t1 = RNGRUP(2)/V  !RO2(nq)
      	t2 = DROSSN(2,M) !diff(t1,nq)
      	t4 = XNOASS(2,2)  !XB2(nq)
      	t5 = t4**2
      	t6 = t2*DELTASS(1,2,2,2)*t5  !t2*DELTA2*t5
      	t7 = t2*DELTASS(1,1,2,2)  !t2*DELTA
      	t9 = t1*DELTASS(1,1,2,2)  !t1*DELTA
      	t15 = RNGRUP(1)/V  ! RO1(nq)
      	t16 = DROSSN(1,M) !diff(t15,nq)
      	t20 = t1*DELTASS(1,2,2,2)
      	t23 = t15*DELTASS(1,1,2,2)
      	t27 = t20*t4
      	t31 = t1**2

      DXNOASSN(2,2,M) = (-t6-t7*t5-2*t9*t5*t4*DELTASS(1,2,2,2)*t2-t4*t16&
     &*DELTASS(1,1,2,2)-t16*DELTASS(1,1,2,2)*t20*t5-t23*t6+t7*t4)/(1+2*&
     &t27+2*t9*t4+3*t31*DELTASS(1,1,2,2)*t5*DELTASS(1,2,2,2)+t23+2*t23*t&
     &27-t9)

      DXNOASSN(1,1,M) =-1D0*XNOASS(1,1)*DELTASS(1,1,2,2)*(RNGRUP(2)/V*&
     &DXNOASSN(2,2,M)+XNOASS(2,2)*DROSSN(2,M))/(1D0+RNGRUP(2)/V*&
     &XNOASS(2,2)*DELTASS(1,1,2,2))

      DXNOASSN(1,2,M) = -1D0*XNOASS(1,2)*DELTASS(1,2,2,2)*(RNGRUP(2)/V*&
     &DXNOASSN(2,2,M)+XNOASS(2,2)*DROSSN(2,M))/(1D0+RNGRUP(2)/V*&
     &XNOASS(2,2)*DELTASS(1,2,2,2))

!      WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)


      END DO

!	DO M=1,2
!		WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)
!		WRITE(2,*)"DXNOASSN(2,2,",M,")=",DXNOASSN(2,2,M)
!	END DO

      END
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE DAUTOASOC(IDAUT,NC,&
     &DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
!   QUE AUTOASOCIA (CON DOS SITIOS).
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DDLTSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR

      IF (DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

        BLOQ =0.D0
        HELPI =0.D0
!
!
        BLOQ = 1D0+2D0*DELTASS(1,IDAUT,2,IDAUT)*XNOASS(1,IDAUT)*&
     &RNGRUP(IDAUT)/V
!
        BLOQ = 1D0 / BLOQ
!
! EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
!                               respecto al DELTASS(1,IDAUT,2,IDAUT)
!
!
        BLOQ = (-1D0)  *  (XNOASS(1,IDAUT))**2D0 * BLOQ
!
! CALCULO DE DERIVADAS RESPECTO AL NUMERO DE MOLES
!
        DO M = 1,NC

                HELPI =BLOQ *DELTASS(1,IDAUT,2,IDAUT)*DROSSN(IDAUT,M)
!
                DXNOASSN(1,IDAUT,M) = HELPI
!
        END DO
!
! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL
! NRO. DE MOLES DEL COMPUESTO M
!
        DO M=1,NC
!
                DXNOASSN(2,IDAUT,M) = DXNOASSN(1,IDAUT,M)

!			WRITE(2,*)"DXNOASSN(1,",IDAUT,M,")=",DXNOASSN(1,IDAUT,M)
!			WRITE(2,*)"DXNOASSN(2,",IDAUT,M,")=",DXNOASSN(2,IDAUT,M)

!
        END DO
      END IF
!
!
        END
!/////////////////////////////////////////////////////////////////////////
        SUBROUTINE DAUTOASOCumsitio(IDAUT,NC,&
     &DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
!   QUE AUTOASOCIA (CON DOS SITIOS).
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DDLTSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR
!        if (DELTASS(1,IDAUT,1,IDAUT).NE.0)THEN
        IF(DELTASS(1,IDAUT,1,IDAUT).NE.0) THEN

        BLOQ =0.D0
        HELPI =0.D0
!
!
        BLOQ = 1D0+2D0*DELTASS(1,IDAUT,1,IDAUT)*XNOASS(1,IDAUT)*&
     &RNGRUP(IDAUT)/V
!
        BLOQ = 1D0 / BLOQ
!
! EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
!                               respecto al DELTASS(1,IDAUT,2,IDAUT)
!
!
        BLOQ = (-1D0)  *  (XNOASS(1,IDAUT))**2D0 * BLOQ
!
! CALCULO DE DERIVADAS RESPECTO AL NUMERO DE MOLES
!
        DO M = 1,NC

                HELPI =BLOQ *DELTASS(1,IDAUT,1,IDAUT)*DROSSN(IDAUT,M)
!
                DXNOASSN(1,IDAUT,M) = HELPI
!			  WRITE(2,*)"DXNOASSN(1,",IDAUT,M,")=",DXNOASSN(1,IDAUT,M)
!
        END DO
      endif

        END
!/////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////
        SUBROUTINE DERIFASSN(NC,NGA,MASSOC,DFAASSXNOASS,DXNOASSN,DFASSN)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
        DIMENSION MASSOC(12)
        DIMENSION DFAASSXNOASS(6,12)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION DFASSN(12,NCM)
!
        HELX = 0.0D0
!
! DFASSN(I,J) = DERIVADA DE FASS DEL GRUPO I CON RESPECTO AL NRO. DE
!               MOLES DEL COMPUESTO J
!
        DO I=1,NGA   !grupos
                DO M=1,NC  !componentes
                        IF(MASSOC(I).NE.0) THEN
!
                                HELX = 0.D0
!
                                  DO K=1,MASSOC(I)
!
         HELX = HELX +  DFAASSXNOASS(K,I) * DXNOASSN(K,I,M)
!
                                  END DO
!
                                DFASSN(I,M) = HELX
                          ELSE
                            CONTINUE
                        END IF
                END DO
        END DO
!
        END
!///////////////////////////////////////////////////////////////////
        SUBROUTINE FUGAASS(NC,NGA,MASSOC,ONUASOC,RNGRUP,XNOASS,FUGASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
        DOUBLE PRECISION XNOASS(6,12)
        DIMENSION ONUASOC(12,NC)
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
        DIMENSION MASSOC(12)
        DIMENSION FUGASS(NCM)
        DOUBLE PRECISION SUMAJ
!
        SUMAJ = 0d0
!
! FUGASS(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD
!               ASOCIATIVO PARA EL COMPONENTE I
!
!
! DFASSN(I,M) = DERIVADA DE FASS DEL grupo I CON RESPECTO AL NRO. DE
!               MOLES DEL COMPUESTO M
!


        DO M=1,NC
!
                       SUMAJ = 0d0
!
                        DO I=1,NGA !grupos

                              IF(MASSOC(I).NE.0) THEN
!

!
                                  DO K=1,MASSOC(I)

      					  If(XNOASS(K,I).gt.1d-13)then

!
       	   SUMAJ = SUMAJ + ONUASOC(I,M) * dlog(XNOASS(K,I))

      					  end if
!
                                  END DO
!
                          ELSE
                            CONTINUE
                        END IF

!

!
                        END DO
!
                       FUGASS(M) = SUMAJ
!
        END DO
!
        END
!///////////////////////////////////////////////////////////////////
        SUBROUTINE COMPRESFV(V,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ZFV)
        IMPLICIT REAL*8 (A-H,O-Z)
!
! ZFV = CONTRIBUCION REPULSIVA AL FACTOR DE COMPRESIBILIDAD
!
        ZFV = 3D0*RLAMBDA1*RLAMBDA2/RLAMBDA3 +&
     &      (RLAMBDA2**3/RLAMBDA3**2)*( 2.D0*Y - 1.D0 - 1D0/Y) + RNTOT/Y
!
        ZFV = ZFV * (-1.D0) * ( V / RNTOT)
!
        ZFV = ZFV * (-1.D0) * Y**2 * 3.141592654D0* RLAMBDA3/(6.D0*V**2)
!
        END
!/////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////
        SUBROUTINE COMPRESDISP(NC,NG,RNTOT,RN,NU,QAREA,H2,H4,H5,H6,&
     &ZDISP)
        IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NGM = 20
        DIMENSION RN(NC)
        DIMENSION NU(NGM,NC), QAREA(NGM)
        DIMENSION H2(NGM)
        DIMENSION H4(NGM), H5(NGM)
        DIMENSION H6(NGM)
        DIMENSION RNURAYA(NGM)
!
        SSUM =0.D0
        SUMM =0.D0
!
        DO I=1,NG
                RNURAYA(I)=0.D0
        END DO
!
!
! ZDISP = CONTRIBUCION DISPERSIVA AL FACTOR DE COMPRESIBILIDAD
!
! RNURAYA(J) ES UN ARREGLO AUXILIAR
!
! CALCULO DE  RNURAYA(J)
!
        DO J = 1,NG
!
              SSUM = 0.D0
!
                DO L=1,NC
                        SSUM = SSUM + RN(L) * NU(J,L)
                END DO
!
              RNURAYA(J) = SSUM
!
        END DO
!
! NU (J,L)= NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO L
!
        SUMM = 0.D0
!
        DO J = 1,NG
!
           SUMM  = RNURAYA(J) * QAREA(J) *&
!
     & ( H5(J) + H2(J) - H2(J)*H6(J)/H4(J) ) /&
!
     &  H4(J) + SUMM
!
        END DO
!
!
        ZETA=10.D0 !NUMERO DE COORDINACION
!
        ZDISP  =  (ZETA/2.D0) * ( -1.D0/ RNTOT ) * SUMM
!
        END
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////
        SUBROUTINE DERIROSSV(NGA,V,MASSOC,DELTASS,DROSSV,RNGRUP)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION DROSSV(12)
        DIMENSION RNGRUP(12)

!
! DROSSV(I,J,K,L) = DERIVADA DEL DELTA DE ASOCIACION ENTRE EL SITIO I
!                       DEL grupo J Y EL SITIO K DEL grupo L CON RESPECTO
!                               AL VOLUMEN
!
!
        DO J=1,NGA
              IF(MASSOC(J).EQ.0) GO TO 201
                        DROSSV(J) = (-1D0*RNGRUP(J))/V**2D0
!						write(2,*)"DROSSV(",J,")=",DROSSV(J)
201           CONTINUE
        END DO
!
        END
!///////////////////////////////////////////////////////////////////
        SUBROUTINE DERIFASSV(NGA,MASSOC,DFAASSXNOASS,DXNOASSV,DFASSV)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION DFAASSXNOASS(6,12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION DFASSV(12)
!
        SUMV = 0.D0
!
! DFASSV(I) = DERIVADA DE FASS(I) CON RESPECTO AL VOLUMEN
!
        DO I =1,NGA
!
                SUMV =0.D0
!               -----------------------------------------
                IF(MASSOC(I).NE.0) THEN
!               -----------------------------------------
!
                          DO K=1,MASSOC(I)

                           SUMV = SUMV + DFAASSXNOASS(K,I)

                          END DO
!               -----------------------------------------
!
      			ELSE
                  CONTINUE
                END IF
!               -----------------------------------------
!
                DFASSV(I) = SUMV
        END DO
!
!
        END
!/////////////////////////////////////////////////////////////////////////
        SUBROUTINE COMPRESASS(NGA,V,RNTOT,RNGRUP,DFASSV,ZASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
        DIMENSION DFASSV(12)
!
! ZASS = CONTRIBUCION ASOCIATIVA AL FACTOR DE COMPRESIBILIDAD
!
        SUMZ = 0D0
!
        DO I=1,NGA
!
                SUMZ = SUMZ + RNGRUP(I)*DFASSV(I)
!
        END DO
!
!
        ZASS =-0.5d0 * SUMZ / RNTOT
!
        END
!///////////////////////////////////////////////////////////////////
        SUBROUTINE FCTRCOMPRES(ZFV,ZDISP,ZASS,ZSIS,IERRZSIS)
        IMPLICIT REAL*8 (A-H,O-Z)
!
! ZSIS = FACTOR DE COMPRESIBILIDAD DEL SISTEMA
!
        ZSIS = ZFV + ZDISP + ZASS + 1.D0
!
        IF(ZSIS.LT.0) THEN
                                IERRZSIS = 1D0

                      ELSE
                                IERRZSIS = 0D0

        END IF
!
        END
!////////////////////////////////////////////////////////////////////
        SUBROUTINE FUGASIS(NC,FUGREP,FUGDISP,FUGASS,ZSIS,FUGSIS,&
     &FFI,IERRFUG)
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
        DIMENSION FUGREP(NCM)
        DIMENSION FUGDISP(NCM)
        DIMENSION FUGASS(NCM)
        DIMENSION FUGSIS(NCM)
        DIMENSION FFI(NC)
!
! FUGSIS(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD DEL COMPONENTE
!                                                                   I
! FFI(I)=  COEFICIENTE DE FUGACIDAD DEL COMPONENTE I
!
        DO I = 1,NC
!
          FUGSIS(I)= FUGREP(I) + FUGDISP(I) + FUGASS(I) - DLOG(ZSIS)
!          evitando overflow
           IF(FUGSIS(I).GT.88) THEN
!
                                IERRFUG = 1
!
                              ELSE
!
                                IERRFUG = 0
!							write(*,*)'ln(phi',i,')=',FUGSIS(I)
                                FFI(I)= DEXP(   FUGSIS(I) )
           END IF
!
        END DO
!
        END
!//////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////
        SUBROUTINE ASSOCNMOL(NC,NGA,RN,ONUASOC,RNGRUP)
        IMPLICIT REAL*8 (A-H,O-Z)
!
        DIMENSION RN(NC)
        DIMENSION ONUASOC(12,NC)
        DIMENSION RNGRUP(12)
!
! Calcula el numero de moles de cada grupo considerado en el termino de
!  asociacion. Se realiza el cálculo tanto para los grupos que presentan
!  sitios de asociacion como para los que no los presentan. Tal vez
!   alguna modificacion futura del modelo introduzca en la expresion
!    de la funcion de Helmholtz de asociacion, el efecto de dilucion
!     producido por los grupos que no asocian sobre los que si lo hacen
!
! RNGRUP(I)= NUMERO DE MOLES DEL GRUPO I (TERMINO DE ASOCIACION)
!

       DO I=1,NGA

                S1 = 0.D0
                DO M=1,NC
!
                           S1 = RN(M) * ONUASOC(I,M) + S1

!
                END DO
                RNGRUP(I) = S1
!
        END DO


   !   write(6,*)"nc=", nc
   !  	do i=1, nga
!		write(6,*)"rngrup=", rngrup(i)
!		do k=1,nc
!			write(6,*) "onuasoc(",k,",",i,")=",onuasoc(k,i)
!		end do
!	end do

!
        END
!//////////////////////////////////////////////////////////////////////////
        SUBROUTINE DERIFAASS(NGA,MASSOC,XNOASS,DFAASSXNOASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION XNOASS(6,12)
        DIMENSION FAASS(6,12)
        DIMENSION DFAASSXNOASS(6,12)
!
!
        DO I=1,NGA
             IF(MASSOC(I).EQ.0) GO TO 13
      		  DO K=1,MASSOC(I)
      			If(XNOASS(K,I).gt.1d-13)then
      				DFAASSXNOASS(K,I) = 1.d0 - XNOASS(K,I)
      			END IF
                END DO
13           CONTINUE
        END DO
!
        END
!/////////////////////////////////////////////////////////////////////////////
        SUBROUTINE DERIXNOASSVL(NGA,MASSOC,DELTASS,XNOASS,&
     &DROSSV,DXNOASSV, ONUASOC, nc,RNGRUP,V)

!
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
! DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
!                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
       IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12), sumsi(12)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION ONUASOC(12,NC)
!	dimension m_lambda(nga*2,nga*2), m_lambda1(nga*2,nga*2)
      dimension psiv(12)
      real*8  m_lambda(nga*2,nga*2), m_lambda1(nga*2,nga*2)
      integer order
!	real d1, psiv
      dimension indx(12)
!     	sumsi=0
!	do k=1, nc
!		do i=1, nga
!			If(ONUASOC(i,k).ne.0.)then
!				sumsi(i)=massoc(i)
!			end if
!		end do
!	end do


!	IF ((DELTASS(1,1,1,1)<1E-13).AND.(DELTASS(1,1,2,2)>1E-13)
!     *.AND.(SUMSI(1)==1)
!     *.AND.(DELTASS(1,2,2,2)>1E-13).AND.(SUMSI(2)==2))THEN

	!ASOCIACION CRUZADA ENTRE 1 GRUPO TIPO ESTER E UM GRUPO TIPO ALCOOL
	!COM UM SITIO (ESTER) E UM GRUPO COM 2 SITIOS(ALCOOL)
	!DELTASS(1,1,1,1)=0 (nao se autoassocia: ester ou alcano)
	!DELTASS(1,1,2,2).NE.0  OU  DELTASS(2,2,1,1).NE.0 SE ASOCIA DE FORMA CRUZADA
	!DELTASS(1,2,2,2).NE.0

!		IDAUT1=1
!		IDAUT2=2
!	   CALL DCROSSASOCVL(IDAUT1,IDAUT2,NC,DELTASS,XNOASS,DROSSV,
!     @DXNOASSV,RNGRUP,V)

!	ELSEIF ((DELTASS(1,2,1,2)<1E-13).AND.(DELTASS(1,1,1,2)>1E-13)
!     *.AND.(SUMSI(1)==1)
!    *.AND.(DELTASS(1,1,1,1)>1E-13).AND.(SUMSI(2)==1))THEN

	!ASOCIACION CRUZADA ENTRE 1 GRUPO TIPO ACIDO E UM GRUPO TIPO ESTER
	!COM UM SITIO (ACIDO) E UM GRUPO COM 1 SITIO(ESTER)
	!DELTASS(1,1,1,1)=0 (nao se autoassocia: ester ou alcano)
	!DELTASS(1,1,1,2).NE.0  OU  DELTASS(1,2,1,1).NE.0 SE ASOCIA DE FORMA CRUZADA
	!DELTASS(1,2,2,2).NE.0

!		IDAUT1=1
!		IDAUT2=2
!	   CALL DCROSSASOCVL3(IDAUT1,IDAUT2,NC,DELTASS,XNOASS,DROSSV,
!     @DXNOASSV,RNGRUP,V)



!	ELSEIF((DELTASS(1,1,1,1)>1E-13).AND.(DELTASS(1,1,2,2)>1E-13)
!     *.AND.(SUMSI(1)>1E-13)
!     *.AND.(DELTASS(1,2,2,2)>1E-13).AND.(SUMSI(2)==2))THEN


!		IDAUT1=1
!		IDAUT2=2
!	   CALL DCROSSASOCVL2(IDAUT1,IDAUT2,NC,
!     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V,nga)

!	ELSEIF((DELTASS(1,1,2,1)>1E-13).AND.(DELTASS(1,1,2,2)>1E-13)
!     *.AND.(SUMSI(1)==2)
!     *.AND.(DELTASS(1,2,2,2)>1E-13).AND.(SUMSI(2)==2))THEN


!		IDAUT1=1
!		IDAUT2=2
!	   CALL DCROSSASOCVL4(IDAUT1,IDAUT2,NC,
!     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)


!	ELSE !AUTO ASOCIACIONES EN PARALELO


!		do m=1, nga

!		  IF(SUMSI(m).EQ.2) THEN
!				  DO I=1,NGA
!					 IF(MASSOC(I).EQ.2) THEN
!						  IDAUT=I
!						  CALL DAUTOASOCVL(IDAUT,
!     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
!
!					  ELSE
!					   CONTINUE
!					 END IF
!				  END DO
!			elseif (sumsi(m).eq.1)then
!				  DO I=1,NGA
!					 IF(MASSOC(I).EQ.1) THEN
!						  IDAUT=I
!
!						  CALL DAUTOASOCVLumsitio(IDAUT,
!     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
!
!					  ELSE
!					   CONTINUE
!					 END IF
!				  END DO
!
!			ELSE
!			 CONTINUE
!		  END IF
!		end do
!
!	END IF


!CCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  General

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
! DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
!                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
! Cálculo de la derivada de la fracción no asociada respecto VOLUMEN
!Calcula los elementos de la matriz deltapq

      w=0;
      s=0
      do j=1,2
      do i=1,nga
      w=w+1
      do l=1, 2
      do k=1, nga
      s=s+1
      m_lambda(w,s)=rngrup(k)*deltass(l,k,j,i)*xnoass(j,i)**2/v
      if (w.eq.s) then
      m_lambda(w,s)=m_lambda(w,s)+ 1.0d0
      end if
      end do
      end do
      s=0
      end do
      end do
      order=nga*2
! Calcula los elementos de la matriz [Yp]
      w=0
      do j=1,2
      do i=1,nga
      w=w+1
      psiv(w)=xnoass(j,i)**2*(1.0d0/xnoass(j,i)-1.0d0)/v
      end do
      end do

      N=order
      NP=order
      m_lambda1=m_lambda
      call  ludcmp(m_lambda1,N,NP,indx,d1)
       call lubksb(m_lambda1,N,NP,indx,psiv)

! colectando las derivadas en correspondiente subíndice
      w=0
      do j=1,2
      do i=1, nga
      		w= w+1
      dxnoassv(j,i)=psiv(w)
      end do
      end do


      END
!CCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOCVL3(IDAUT1,IDAUT2,NC,&
     &DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION RNGRUP(12)


      	t1 = RNGRUP(1)/V
      	t2 = DROSSV(1)  !diff(t1,nq)
      	t3 = t2* DELTASS(1,1,1,2) !delcross
      	t4 = XNOASS(1,1)  !XA1(nq)
      	t5 = t4**2
      	t9 = t1*DELTASS(1,1,1,1)   !DELTA1
      	t15 = RNGRUP(2)/V   !RO2(nq)
      	t16 = DROSSV(2)  !diff(t15,nq)
      	t21 = t1*DELTASS(1,1,1,2)  !delcross
      	t26 = t1**2
      DXNOASSV(1,1)= (-t3*t5-t2*DELTASS(1,1,1,1)*t5-2D0*t9*t5*t4*&
     &DELTASS(1,1,1,2)*t2-t4*t16*DELTASS(1,1,1,2)+t3*t4)/(1D0+2D0*t21*t4&
     &+2D0*t9*t4+3D0*t26*DELTASS(1,1,1,1)*t5*DELTASS(1,1,1,2)+t15*&
     &DELTASS(1,1,1,2)-t21)



      DXNOASSV(1,2)=-1D0*XNOASS(1,2)*DELTASS(1,1,1,2)*&
     &(DROSSV(1)*XNOASS(1,1)+&
     & (RNGRUP(1)/V)*DXNOASSV(1,1))&
     &/(1D0+(RNGRUP(1)/V)*DELTASS(1,1,1,2)*XNOASS(1,1))

!     WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)
!	WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)

      END
!////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOCVL4(IDAUT1,IDAUT2,NC,&
     &DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR




! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL
! NRO. DE MOLES DEL COMPUESTO M

!

      	t1 = RNGRUP(2)/V
      	t2 = DELTASS(1,1,1,1)*t1
      	t3 = DELTASS(1,2,2,2)**2
      	t4 = DROSSV(2)  !diff(t1,nq)
      	t8 = DELTASS(1,1,2,2)**2
      	t9 = t1*t8
      	t14 = XNOASS(2,2) !X(nq)
      	t15 = t14**2
      	t16 = t15**2
      	t21 = 2*t4*t8
      	t24 = 2*DELTASS(1,1,1,1)*t4*DELTASS(1,2,2,2)
      	t26 = t15*t14
      	t28 = RNGRUP(1)/V
      	t29 = DROSSV(1)  !diff(t28,nq)
      	t34 = t1**2
      	t45 = 2*t9
      	t47 = 2*t2*DELTASS(1,2,2,2)
      DXNOASSV(2,2) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-&
     &(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+&
     &t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*&
     &DELTASS(1,2,2,2))*&
     &t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*&
     &(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4&
     &5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


      DXNOASSV(1,1)=(-DROSSV(1)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2D0&
     &*XNOASS(1,1)*DROSSV(2)*DELTASS(1,1,2,2)*XNOASS(2,2)-2D0*&
     &XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSV(2,2))&
     &/(1+2D0*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2D0*&
     &(RNGRUP(2)/V)*DELTASS(1,1,2,2)*XNOASS(2,2))



      DXNOASSV(1,2) = DXNOASSV(2,2)
!      WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)

!		WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)
!		WRITE(2,*)"DXNOASSV(2,2)=",DXNOASSV(2,2)


      END
!************************************************************************************
!////////////////////////////////////////////////////////////////////////
       SUBROUTINE DCROSSASOCVL2(IDAUT1,IDAUT2,NC,&
     &DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION RNGRUP(12)

! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
! DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
!                  EL SITIO I DEL GRUPO J, Y EL SITIO K
!                  DEL GRUPO L con respecto al numero de moles
!                  del componente M
!
! DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
!                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
!  CALCULO DE BLOQUE AUXILIAR




! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL
! NRO. DE MOLES DEL COMPUESTO M
!      xnoass(1,1)=0.5873
!     xnoass(1,2)=0.9437
!	xnoass(2,2)=0.9437
!      rngrup(1)=1
!	rngrup(2)=0
!	deltass(1,1,1,1)=81256.6589
!	deltass(1,2,1,1)=6899.99
!	deltass(2,2,1,1)=6899.99
!	deltass(1,1,1,2)=6899.99
!	deltass(2,2,1,2)=585.15
!	deltass(1,1,2,2)=6899.993
!	deltass(1,2,2,2)=585.1586
!	deltass(2,2,2,2)=0
!	v=67937.4025

      	t1 = RNGRUP(2)/V
      	t2 = DELTASS(1,1,1,1)*t1
      	t3 = DELTASS(1,2,2,2)**2
      	t4 = DROSSV(2)  !diff(t1,nq)
      	t8 = DELTASS(1,1,2,2)**2
      	t9 = t1*t8
      	t14 = XNOASS(2,2) !X(nq)
      	t15 = t14**2
      	t16 = t15**2
      	t21 = 2*t4*t8
      	t24 = 2*DELTASS(1,1,1,1)*t4*DELTASS(1,2,2,2)
      	t26 = t15*t14
      	t28 = RNGRUP(1)/V
      	t29 = DROSSV(1)  !diff(t28,nq)
      	t34 = t1**2
      	t45 = 2*t9
      	t47 = 2*t2*DELTASS(1,2,2,2)
      DXNOASSV(2,2) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-&
     &(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+&
     &t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*&
     &DELTASS(1,2,2,2))*&
     &t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*&
     &(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4&
     &5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


      DXNOASSV(1,1)=(-DROSSV(1)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2D0&
     &*XNOASS(1,1)*DROSSV(2)*DELTASS(1,1,2,2)*XNOASS(2,2)-2D0*&
     &XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSV(2,2))&
     &/(1+2D0*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2D0*&
     &(RNGRUP(2)/V)*DELTASS(1,1,2,2)*XNOASS(2,2))



      DXNOASSV(1,2) = DXNOASSV(2,2)
!      WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)

!		WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)
!		WRITE(2,*)"DXNOASSV(2,2)=",DXNOASSV(2,2)


      END




!************************************************************************************
      SUBROUTINE DCROSSASOCVL(IDAUT1,IDAUT2,NC,DELTASS,XNOASS,DROSSV,&
     &DXNOASSV,RNGRUP,V)

        IMPLICIT REAL*8 (A-H,O-Z)
!
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)


      t1 = RNGRUP(2)/V  !RO2(nq)
      t2 = DROSSV(2) !diff(t1,nq)
      t4 = XNOASS(2,2)  !XB2(nq)
      t5 = t4**2
      t6 = t2*DELTASS(1,2,2,2)*t5  !t2*DELTA2*t5
      t7 = t2*DELTASS(1,1,2,2)  !t2*DELTA
      t9 = t1*DELTASS(1,1,2,2)  !t1*DELTA
      t15 = RNGRUP(1)/V  ! RO1(nq)
      t16 = DROSSV(1) !diff(t15,nq)
      t20 = t1*DELTASS(1,2,2,2)
      t23 = t15*DELTASS(1,1,2,2)
      t27 = t20*t4
      t31 = t1**2

      DXNOASSV(2,2) = (-t6-t7*t5-2d0*t9*t5*t4*DELTASS(1,2,2,2)*t2-t4*t16&
     &*DELTASS(1,1,2,2)-t16*DELTASS(1,1,2,2)*t20*t5-t23*t6+t7*t4)/(1d0+&
     &2d0*t27+2d0*t9*t4+3d0*t31*DELTASS(1,1,2,2)*t5*DELTASS(1,2,2,2)+t23&
     &+2d0*t23*t27-t9)

      DXNOASSV(1,1) =-1D0*XNOASS(1,1)*DELTASS(1,1,2,2)*(RNGRUP(2)/V*&
     &DXNOASSV(2,2)+XNOASS(2,2)*DROSSV(2))/(1D0+RNGRUP(2)/V*&
     &XNOASS(2,2)*DELTASS(1,1,2,2))

      DXNOASSV(1,2) = -1D0*XNOASS(1,2)*DELTASS(1,2,2,2)*(RNGRUP(2)/V*&
     &DXNOASSV(2,2)+XNOASS(2,2)*DROSSV(2))/(1D0+RNGRUP(2)/V*&
     &XNOASS(2,2)*DELTASS(1,2,2,2))

!	WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)
!      WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)
!      WRITE(2,*)"DXNOASSV(2,2)=",DXNOASSV(2,2)


      END
!//////////////////////////////////////////////////////////////////////////

!//////////////////////////////////////////////////////////////////////////


        SUBROUTINE DAUTOASOCVL(IDAUT,DELTASS,XNOASS,DROSSV,DXNOASSV,&
     &  RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
!
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)
!
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL VOLUMEN DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
!
!
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
!   QUE AUTOASOCIA (CON DOS SITIOS).
!
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
!
!
! DDLTSSV(I,J,K,L) = DERIVADA DEL DELTA DE ASOCIACION ENTRE EL SITIO I
!                       DEL grupo J Y EL SITIO K DEL grupo L CON RESPECTO
!                               AL VOLUMEN
!
! DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
!                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
!
!  CALCULO DE BLOQUE AUXILIAR

        IF(DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

!
        BLOQ =0.D0
        HELPI =0.D0
!
!
        BLOQ =1D0+2D0*DELTASS(1,IDAUT,2,IDAUT)*XNOASS(1,IDAUT)*&
     &RNGRUP(IDAUT)/V
!
        BLOQ = 1D0 / BLOQ
!
! EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
!                               respecto al DELTASS(1,IDAUT,2,IDAUT)
!
!
        BLOQ=(-1D0)*DELTASS(1,IDAUT,2,IDAUT)*(XNOASS(1,IDAUT))**2D0*BLOQ
!
! CALCULO DE LA DERIVADA RESPECTO VOLUMEN
!
!
        DXNOASSV(1,IDAUT) = BLOQ  * DROSSV(IDAUT)
!
!
! IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
!  IDAUT CON RESPECTO AL VOLUMEN
!
!
        DXNOASSV(2,IDAUT)   = DXNOASSV(1,IDAUT)

!		WRITE(2,*)"DXNOASSV(1,",IDAUT,")=",DXNOASSV(1,IDAUT)
!		WRITE(2,*)"DXNOASSV(2,",IDAUT,")=",DXNOASSV(2,IDAUT)
      END IF
!
!
!
        END
!//////////////////////////////////////////////////////////////////////////
        SUBROUTINE DAUTOASOCVLumsitio(IDAUT,DELTASS,XNOASS,DROSSV,&
     &DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
!
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)
!
!
! NO GENERAL
!
! CALCULA DERIV. RESP. AL VOLUMEN DE LA FRACCION NO ASOCIADA
!  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
! SITIOS)
!
!
! XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
!
!
! IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
!   QUE AUTOASOCIA (CON DOS SITIOS).
!
! XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
!                                                       IDAUT
!
!
! DDLTSSV(I,J,K,L) = DERIVADA DEL DELTA DE ASOCIACION ENTRE EL SITIO I
!                       DEL grupo J Y EL SITIO K DEL grupo L CON RESPECTO
!                               AL VOLUMEN
!
! DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
!                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
!
!  CALCULO DE BLOQUE AUXILIAR
!
        IF(DELTASS(1,IDAUT,1,IDAUT).NE.0) THEN

        BLOQ =0.D0
        HELPI =0.D0
!
!
        BLOQ =1D0+2D0*DELTASS(1,IDAUT,1,IDAUT)*XNOASS(1,IDAUT)*&
     &RNGRUP(IDAUT)/V
!
        BLOQ = 1D0 / BLOQ
!
! EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
!                               respecto al DELTASS(1,IDAUT,2,IDAUT)
!
!
        BLOQ=(-1D0)*DELTASS(1,IDAUT,1,IDAUT)*(XNOASS(1,IDAUT))**2D0*BLOQ
!
! CALCULO DE LA DERIVADA RESPECTO VOLUMEN
!
!
        DXNOASSV(1,IDAUT) = BLOQ  * DROSSV(IDAUT)
! 	WRITE(2,*)"DXNOASSV(1,",IDAUT,")=",DXNOASSV(1,IDAUT)
      END IF
!
        END
!//////////////////////////////////////////////////////////////////////////

        SUBROUTINE &
     &ACLAR(NCOD,NC,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)
!
! PONE EN CERO LAS VARIABLES DE SALIDA CUANDO COMIENZA LA CORRIDA
!
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION FFI(NC)  !DIMENSIONAMIENTO COEFICIENTES DE FUGACIDAD
!
        ARESTOT=0.D0
!
! SI NCOD=2, ZSIS entra como dato, por lo que no se blanquea esa variable
!
        IF(NCOD.NE.2) ZSIS =0.D0
!
        DO I=1,NC
                FFI(I)=0.D0
        END DO
!
        IERRFUG=0
        IERRZSIS=0
        IERRREP=0
      iErrAssoc=0
!
!
        END !FIN ACLAR
!//////////////////////////////////////////////////////////////////
        SUBROUTINE CLEARANCE(NC,NG,NGA,ARESASS,ARESDISP,ARESFV,QTOT,&
     &RLAMBDA1,RLAMBDA2,RLAMBDA3,RNTOT,Y,ZASS,ZDISP,ZFV,&
     &MASSOC,D,DL1N,DL2N,DL3N,DYN,FUGASS,FUGDISP,FUGREP,FUGSIS,&
     &RMS,H3,H7,PS,DFASSN,DXNOASSN,H2,H4,H5,H6,TITA,TSPLINE,&
     &DLTG,TAU,VERDKIJ,GGG,DFASSV,FASS,RNGRUP,DXNOASSV,FAASS,XNOASS,&
     &DFAASSXNOASS,DROSSV,DELTASS,DROSSN)
!
! PONE EN CERO UN NUMERO IMPORTANTE DE VARIABLES
!
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50
!
        DIMENSION D(NCM)
      integer, parameter :: NGM = 20
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
! DEL TERMINO DE DISPERSION
!
        DIMENSION TITA(NGM), TAU(NGM,NGM)
        DIMENSION TSPLINE(NGM)
        DIMENSION GGG(NGM,NGM)
        DIMENSION VERDKIJ(NGM,NGM)
        DIMENSION DLTG(NGM,NGM)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
!
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
!
! MASSOC(J) = NRO. DE SITIOS DE ASOCIACION ASIGNADOS AL GRUPO J
!
! DELTASS(I,J,K,L) = FUERZA ADIMENSIONAL DE ASOCIACION ENTRE  el sitio
!                 I del GRUPO J, y el sitio K del grupo L
!
!
        DIMENSION MASSOC(12)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION FASS(12)
        DIMENSION FAASS(6,12)
!
! DIMENSIONAMIENTO RELAC CON COEF FUGACIDAD REPULSIVO
!
        DIMENSION  DL1N(NCM),DL2N(NCM),DL3N(NCM),DYN(NCM)
        DIMENSION FUGREP(NCM)
!
! DIMENSIONAMIENTO DE VARIABLES PARA EL CALCULO
!    DEL COEFICIENTE DE FUGACIDAD DE DISPERSIVO
!
        DIMENSION PS(NCM,NGM), RMS(NCM), H2(NGM)
        DIMENSION H3(NCM,NGM), H4(NGM), H5(NGM)
        DIMENSION H6(NGM), H7(NCM,NGM)
        DIMENSION FUGDISP(NCM)
!
! DIMENSIONAMIENTO RELAC. CON COEF. FUG. ASOCIATIVO
!
        DIMENSION DROSSN(12,NCM)
        DIMENSION DXNOASSN(6,12,NCM)
        DIMENSION DFAASSXNOASS(6,12)
        DIMENSION DFASSN(12,NCM)
        DIMENSION FUGASS(NCM)
!
! DIMENSIONAMIENTO RELAC. CON FACTOR DE COMPRES. ASOCIATIVO
!
!        DIMENSION DROSSV(6,12,6,12)
        DIMENSION DROSSV(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION DFASSV(12)
!
!
        DIMENSION FUGSIS(NCM)
!
!
! LIMPIEZA DE VARIABLES
!
!
!         VARIABLES SIMPLES
!
!         ARESASS                     =0.
          ARESDISP                    =0.D0
          ARESFV                      =0.D0
          QTOT                        =0.D0
          RLAMBDA1                    =0.D0
          RLAMBDA2                    =0.D0
          RLAMBDA3                    =0.D0
          RNTOT                       =0.D0
          Y                           =0.D0
          ZASS                        =0.D0
          ZDISP                       =0.D0
          ZFV                         =0.D0

!       Arreglos
!
! LA PROXIMA ASIGNACION PERMITE ASEGURAR EL VACIADO DE LOS ARREGLOS
!  CUANDO LA HISTORIA DEL CALCULO FUE TAL , QUE PRIMERO SE REALIZO UN
!    CALCULO' CON CIERTO NUMERO DE COMPONENTES, Y LUEGO CON UN NUMERO
!      MENOR EN UNA UNIDAD ( P.E. primero un flash binario, y luego una
!          presion de vapor de un componente puro)
!
        NCSAVE = NC + 1
!
!
!
        DO I=1,NCSAVE
!
          D(I)=0.D0
          DL1N(I)=0.D0
          DL2N(I)=0.D0
          DL3N(I)=0.D0
          DYN(I)=0.D0
          FUGASS(I)=0.D0
          FUGDISP(I)=0.D0
          FUGREP(I)=0.D0
          FUGSIS(I)=0.D0
          RMS(I)=0.D0
!
           DO J=1,NG
                H3(I,J)=0.D0
                H7(I,J)=0.D0
                PS(I,J)=0.D0
           END DO
!
            DO K=1,NGA
!
                        DFASSN(K,I)=0.D0
!
                        DO L=1,6
!
                                DXNOASSN(L,K,I)=0.D0
!
                        END DO
!
            END DO
!
        END DO
!       --------------------
!
!
        DO I=1,NG
!
          H2(I)=0.D0
          H4(I)=0.D0
          H5(I)=0.D0
          H6(I)=0.D0
          TITA(I)=0.D0
          TSPLINE(I)=0.D0
!
            DO J=1,NG
                DLTG(I,J)=0.D0
                TAU(I,J)=0.D0
                VERDKIJ(I,J)=0.D0
                GGG(I,J)=0.D0
            END DO
!
        END DO
!
!
        DO I=1,NGA
!
         DFASSV(I)=0.D0
         FASS(I)=0.D0
         RNGRUP(I)=0.D0
!
                DO J=1,6
!
                        DXNOASSV(J,I)=0.D0
                        FAASS(J,I)=0.D0
                        XNOASS(J,I)=0.D0
                        DFAASSXNOASS(J,I)=0.D0
!
                END DO
!
        END DO
!
        DROSSV(:12)=0.D0
        DO I=1,NGA
                DO J=1,6
                        DO K=1,NGA
                                DO L=1,6
!
!                                   DROSSV(J,I,L,K)=0.D0
                                   DELTASS(J,I,L,K)=0.D0
!
                                   DO M=1,NCSAVE
!
                                        DROSSN(J,M)=0.D0
!
                                   END DO
                                END DO
                        END DO
                END DO
        END DO
!
        END
