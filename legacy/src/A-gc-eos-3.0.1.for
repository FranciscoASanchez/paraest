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
	SUBROUTINE ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,
     @TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)



	!MODIFICAR CROSSASOC4,DCROSSASOC4,DCROSSASOCVL4
C
C	NCOD es un codigo que especifica que' es lo que se quiere calcular:
C		(se introdujo NCOD para ahorrar tiempo de computo)
C
C     NCOD      SIGNIFICADO
C
C      1         factor de compresibilidad solamente
C      2         coeficientes de fugacidad con factor de compresibilidad
C                   entrado como dato (no es posible calcular el coef.
C                    de fugacidad si no se conoce el factor de compres)
C                    Por supuesto el volumen y la temperatura tb. entran como
C			dato
C      3         coeficientes de fugacidad y el necesario factor de
C	          compresibilidad (calculados ambos dentro de esta rutina)
C      4	 Energia de Helmholtz (util para chequeo de derivadas
C	 	analiticas a traves de derivadas numericas)
C      5         factor de compresibilidad, coeficientes de fugacidad y
C		 energia de Helmholtz (es decir todo)
C
C
C
C
C    ARESTOT      = Energia residual total
C
C                      R
C    ARESTOT      = ( A  / RT )         (las unidades son "moles")
C
C
C NGA = NRO DE GRUPOS DISTINTOS EN EL TERMINO DE ASOCIACION
C NU (J,L)= NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO L (PARA EL
C                                       TERMINO DISPERSIVO)
C ONUASOC(J,L) = NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO L (PARA
C                               EL TERMINO DE ASOCIACION)
C
C CALCULA EL FACTOR DE COMPRESIBILIDAD (ZSIS) Y LOS COEFICIENTES DE
C FUGACIDAD DE CADA COMPONENTE [FFI(I)] y la egia residual total, EN
C FUNCION DEL NUMERO DE MOLES  [RN(I)], LA TEMPERATURA(TEMP) Y EL
C VOLUMEN DEL SISTEMA (V)
C
C EL MODELO ES ASSOC-GC-EOS (NO CONFUNDIR CON ASS-GC-EOS). LA
C  CONTRIBUCIONES repulsiva y dispersiva son identicas a las de la
C  GCEOS. El te'rmino de asociacion es a CONTRIBUCION GRUPAL
C   EL NRO. Y TIPO DE GRUPOS PARA UNA DADA MOLECULA ES DISTINTO
C     SEGUN SE TRATE DE LA CONTRIBUCION DISPERSIVA O DE ASOCIACION
C
c        IMPLICIT REAL*8 (A-H,O-Z)

	IMPLICIT NONE

C
C SE manejan:
C               NC componentes (aunque en realidad el maximo es 20)
C               10 grupos diferentes desde el punto de vista del
C                                                       termino dispersivo
C                       12  grupos diferentes desde el punto de vista del
C                                                       termino asociativo
C                        6 sitios de  asociacion por grupo
C
C DC tiene dimension NC porque DC es un argumento de la subrutina
C D   es un arreglo interno , por lo que su dimension es (debe ser) fija e
C       igual al maximo nro. de componentes
C
        INTEGER NC
	  DOUBLE PRECISION DC(NC), D(20), RN(NC)
        DOUBLE PRECISION  TEMPCRIT(NC)
C
C DOUBLE PRECISIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
C DEL TERMINO DE DISPERSION
C
	  INTEGER NU(10,NC)
        DOUBLE PRECISION  QAREA(10), TSTR(10), GSTR(10)
        DOUBLE PRECISION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
        DOUBLE PRECISION ALFAIJ(10,10)
        DOUBLE PRECISION TITA(10), TAU(10,10)
        DOUBLE PRECISION TSPLINE(10)
        DOUBLE PRECISION GGG(10,10)
        DOUBLE PRECISION VERDKIJ(10,10)
        DOUBLE PRECISION DLTG(10,10)
C
C DOUBLE PRECISIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
C
        DOUBLE PRECISION ONUASOC(12,NC)
        DOUBLE PRECISION RNGRUP(12)  !nro. de moles de los grupos
C
C MASSOC(J) = NRO. DE SITIOS DE ASOCIACION ASIGNADOS AL GRUPO J
C ENASS(I,J,K,L)= PARAMETRO DE ENERGIA DE ASOCIACION entre el sitio
C                 I del GRUPO J, y el sitio K del grupo L
C RKASS(I,J,K,L)= PARAMETRO DE VOLUMEN DE ASOCIACION entre el sitio
C                 I del GRUPO J, y el sitio K del grupo L
C
C DELTASS(I,J,K,L) = FUERZA ADOUBLE PRECISIONAL DE ASOCIACION ENTRE  el sitio
C                 I del GRUPO J, y el sitio K del grupo L
C
C
	  INTEGER MASSOC(12)
        DOUBLE PRECISION ENASS(6,12,6,12), RKASS(6,12,6,12)
        DOUBLE PRECISION DELTASS(6,12,6,12)
        DOUBLE PRECISION XNOASS(6,12)
        DOUBLE PRECISION FASS(12)
        DOUBLE PRECISION FAASS(6,12)
C
C DOUBLE PRECISIONAMIENTO RELAC CON COEF FUGACIDAD REPULSIVO
C
        DOUBLE PRECISION  DL1N(20),DL2N(20),DL3N(20),DYN(20)
        DOUBLE PRECISION FUGREP(20)
C
C DOUBLE PRECISIONAMIENTO DE VARIABLES PARA EL CALCULO
C    DEL COEFICIENTE DE FUGACIDAD DE DISPERSIVO
C
        DOUBLE PRECISION PS(20,10), RMS(20), H2(10)
        DOUBLE PRECISION H3(20,10), H4(10), H5(10)
        DOUBLE PRECISION H6(10), H7(20,10)
        DOUBLE PRECISION FUGDISP(20)
C
C DOUBLE PRECISIONAMIENTO RELAC. CON COEF. FUG. ASOCIATIVO
C
        DOUBLE PRECISION DROSSN(12,20)
        DOUBLE PRECISION DXNOASSN(6,12,20)
        DOUBLE PRECISION DFAASSXNOASS(6,12)
        DOUBLE PRECISION DFASSN(12,20)
        DOUBLE PRECISION FUGASS(20)
C
C DOUBLE PRECISIONAMIENTO RELAC. CON FACTOR DE COMPRES. ASOCIATIVO
C
        DOUBLE PRECISION DROSSV(12)
        DOUBLE PRECISION DXNOASSV(6,12)
        DOUBLE PRECISION DFASSV(12)
C
        DOUBLE PRECISION FFI(NC)  !DOUBLE PRECISIONAMIENTO COEFICIENTES DE FUGACIDAD
C
C
C
        DOUBLE PRECISION FUGSIS(20)

	  INTEGER NCOD,IERRFUG,IERRZSIS,IERRREP,NG,NGA,L,NGA2
	  DOUBLE PRECISION ARESTOT,ZSIS,ARESASS,ARESDISP,ARESFV,QTOT,
     @RLAMBDA1,RLAMBDA2,RLAMBDA3,RNTOT,Y,ZASS,ZDISP,ZFV,TEMP,V

	  INTEGER ierrassoc, J, KK
	  DOUBLE PRECISION TEMP1
C
C
C       TEMP= TEMPERATURA DEL SISTEMA EN K
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
C       V = VOLUMEN DEL SISTEMA EN CM3
C
C       NC  = Numero de componentes (mAX 20)
C
C       RN (I) = NUMERO DE MOLES DEL COMPONENTE I
C
C
C  V = VOLUMEN DEL SISTEMA
C      V debe tener unidades congruentes con las de D. Es decir si el
C      diametro de esferas rigidas esta dado en cm/mol el volumen V debe
C      estar en cm3.
C
C  Y = VARIABLE INTERMEDIA EN LA EXPRESION DE MANSOORI Y LELAND (1972)
C
C
C  PARA EVITAR ERRORES
C
	IF(NCOD.LT.1.OR.NCOD.GT.5) NCOD = 5  !CALCULA TODO
C
C  AQUI SE OBLIGA A TODOS LOS ARREGLOS QUE NO SON DATO A TENER VALORES
C  NULOS
C
C AQUI SE LIMPIAN TODAS LAS VARIABLES DE SALIDA
C       Excepcion: Si NCOD es 2 ZSIS entra como dato. Es decir que en
C		    ese caso no se limpia ZSIS
C
      CALL ACLAR(NCOD,NC,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,
     @ierrassoc)
C
C Limpieza de variables locales desde el punto de vista de la rutina principal
C (no de las sub-sub-rutinas
C
        CALL CLEARANCE(NC,NG,NGA,ARESASS,ARESDISP,ARESFV,QTOT,
     @RLAMBDA1,RLAMBDA2,RLAMBDA3,RNTOT,Y,ZASS,ZDISP,ZFV,
     @MASSOC,D,DL1N,DL2N,DL3N,DYN,FUGASS,FUGDISP,FUGREP,FUGSIS,
     @RMS,H3,H7,PS,DFASSN,DXNOASSN,H2,H4,H5,H6,TITA,TSPLINE,
     @DLTG,TAU,VERDKIJ,GGG,DFASSV,FASS,RNGRUP,DXNOASSV,FAASS,XNOASS,
     @DFAASSXNOASS,DROSSV,DELTASS,DROSSN)
C
C CALCULO DEL NUMERO TOTAL DE MOLES (RNTOT=NUMERO TOTAL DE MOLES)

c	OPEN(UNIT=4,NAME='DESVIOS.TXT',FORM='FORMATTED')

c	 WRITE(4,*) "TEMP,V",TEMP,V
	if(V.lt.1d-13) then
		IERRZSIS=1
		goto 1000
	end if
C
        CALL SUMATORIA(NC,RN,RNTOT)
C
C CALCULO DEL NUMERO TOTAL DE SEGMENTOS DE AREA (PARA CONTRIB. POR DISP)
C
        CALL SURFSEG(NC,NG,NU,RN,QAREA,QTOT)
C
C CALCULO DE LAS FRACCIONES DE SUPERFICIE (PARA CONTRIB. POR DISP.)
C
        CALL FRACAREA(NC,NG,NU,RN,QAREA,QTOT,TITA)
C
C CALCULO DE LA TEMPERATURA A LA QUE SE REALIZA EL SPLINE
C
        CALL SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
C
C CALCULO DE PARAMETROS ENERGETICOS ATRACTIVOS (LIKE)
C
        CALL CALGGG(NG,TEMP,TSTR,GSTR,GDOT,GDDOT,TSPLINE,GGG)
C
C CALCULO DE LOS KIJ (VERDKIJ==>VERDADERO KIJ)
C
        CALL CALCKIJ(NG,RKIJ,RKTIJ,TSTR,VERDKIJ,TEMP)
C
C
C CALCULO DE PARAMETROS ENERGETICOS ATRACTIVOS (unLIKE)
C
        CALL CALGGGIJ(NG,VERDKIJ,GGG)
C
C CALCULO DE DELTAgij
C
        CALL DELTADISP(NG,GGG,DLTG)
C
C CALCULO DE TAU(J,I) (NECESARIO PARA CONTRIB.  FZAS DISPERSION

C
        CALL TAUCALC(NG,TEMP,V,ALFAIJ,QTOT,DLTG,TAU)
C
C CALCULO DE LA ENERGIA RESIDUAL DE HELMHOLTZ DE DISPERSION
C
C SEGUN GC-EOS ORIGINAL
C
        IF(NCOD.GT.3) CALL
     @AHELMRESDISP(NC,NG,RN,V,TEMP,NU,TITA,TAU,QTOT,
     @GGG,QAREA,ARESDISP)
C
C
C CALCULO DE D(I) (DIAMETROS DE ESFERAS RIGIDAS A TEMPERATURA TEMP
C
        CALL SOFT(NC,DC,TEMPCRIT,TEMP,D)
C
C CALCULO DE LOS LAMBDA NECESARIOS PARA LA CONTRIBUCION REPULSIVA
C
        CALL LAMBDA(NC,RN,D,RLAMBDA1,RLAMBDA2,RLAMBDA3)
C
C CALCULO DE LA VARIABLE INTERMEDIA Y, NECESARIA PARA LA CONTRIB.
C  RPULSIVA
C
        CALL YCALC(V,RLAMBDA3,Y,IERRREP)
C
C TERMINA LA EJECUCION SI LA FRACCION DE VOLUMEN OCUPADA POR LAS
C  ESFERAS BLANDAS ES MAYOR A LA UNIDAD ( ES DECIR SI EL VOLUMEN
C                                             ES MUY CHICO)
C
        IF(IERRREP.EQ.1) GO TO 1000
C
C  CALCULO DE LA ENERGIA RESIDUAL DE HELMHOLTZ DE FREE VOLUME SEGUN
C GC-EOS ORIGINAL
C
C
        IF(NCOD.GT.3)  CALL
     @AHELMRESFV(RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ARESFV)
C

      If (nc==1) then

	  do l=1,nga
	    !print *, "ONUASOC(1,1)=", ONUASOC(l,1)
	    if(ONUASOC(l,1).ne.0) then

	      nga2 = nga
	      go to 1003

	    end if

	  end do
	  nga2=0
 1003   continue

	else

	  nga2=nga

      end if



C   CALCULO DEL NUMERO DE MOLES DE CADA GRUPO EN LA CONTRIBUCION DE
C               ASOCIACION
C

	   IF(nga2.NE.0) THEN  !AVOID ASSOC 1

        CALL ASSOCNMOL(NC,nga2,RN,ONUASOC,RNGRUP)
C
C Calculo del Delta de asociacion (para nosotros es adimensional
C  porque incluye la densidad de cierto grupo)
C
C
        CALL CALCDELTASS(nga2,TEMP,MASSOC,ENASS,RKASS,DELTASS)
C
C CALCULO DE LAS FRACCIONES DE SITIOS NO ASOCIADOS (XNOASS)
C
        CALL FRACNOAS(nga2,MASSOC,DELTASS,XNOASS, nc, ONUASOC,V,RNGRUP
     @,ierrassoc)


	TEMP1=RNGRUP(1)



C
C CALCULO DE LA ENERGIA DE HELMHOLTZ RESIDUAL DE ASOCIACION
C
C  FASS(I) Y FAASS(I) SON  VARIABLES INTERMEDIAS DE LAS QUE MAS ADELANTE
C                       SE OBTIENEN DERIVADAS
C
        IF(NCOD.NE.1)  CALL
     @      AHELMRESASS(nga2,MASSOC,RNGRUP,XNOASS,ARESASS,FASS,FAASS)
C
	END IF  !AVOID ASSOC 1
C
C CALCULO DE LA ENERGIA RESIDUAL TOTAL
C
	IF(NCOD.GT.3)  CALL AHELMRESTOT(ARESFV,ARESDISP,ARESASS,ARESTOT)
C
 	IF(NCOD.EQ.4) GO TO 1000
C
C CALCULOS AUXILIARES NECESARIOS PARA EL CALCULO DE LA
C  CONTRIBUCION REPULSIVA AL COEFICIENTE DE FUGACIDAD
C
        IF(NCOD.NE.1)  CALL DERIVFV(NC,V,D,Y,DL1N,DL2N,DL3N,DYN)
C
C CALCULO DE LA CONTRIBUCION REPULSIVA AL COEFICIENTE DE FUGACIDAD
C
        IF(NCOD.NE.1)  CALL
     @FUGAREP(NC,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,
     @DL1N,DL2N,DL3N,DYN,FUGREP)

C
C CALCULOS AUXILIARES NECESARIOS PARA EL CALCULO DE LA
C  CONTRIBUCION DISPERSIVA AL COEFICIENTE DE FUGACIDAD
C
        CALL HCALC(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,
     @TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
C
C CALCULO DE LA CONTRIBUCION DISPERSIVA AL COEFICIENTE DE FUGACIDAD
C
        IF(NCOD.NE.1)  CALL
     @FUGADISP(NC,NG,TITA,PS,RMS,H2,H3,H4,H5,H6,H7,FUGDISP)

C
C
	IF(nga2.NE.0) THEN  !AVOID ASSOC 2
C
C CALCULO DE DERIVADAS DE LA DENSIDAD CON RESPECTO
C    AL NUMERO DE MOLES
C
        IF(NCOD.NE.1)  CALL
     @DERIROSSN(NC,nga2,MASSOC,ONUASOC,RNGRUP,DROSSN,V)

C
C CALCULO DE LAS DERIVADAS DE LAS FRACCIONES NO ASOCIADAS CON RESPECTO
C  AL NUMERO DE MOLES entra DROSSN e sai DXNOASSN
C
        IF(NCOD.NE.1)  CALL DERIXNOASS(NC,nga2,MASSOC,DELTASS,XNOASS,
     @DROSSN,DXNOASSN,ONUASOC,RNGRUP,V)
	!		 write (6,*) "cheguei"

C
C
C CALCULO DE LAS DERIVADAS DE FAASS(I) CON RESPECTO A LAS FRACCIONES NO
C       ASOCIADAS
C
        CALL DERIFAASS(nga2,MASSOC,XNOASS,DFAASSXNOASS)
C
C CALCULO DE LAS DERIVADAS DE FASS(I) CON RESPECTO A LOS NUMEROS DE MOLES
C
        IF(NCOD.NE.1)  CALL
     @DERIFASSN(NC,nga2,MASSOC,DFAASSXNOASS,DXNOASSN,DFASSN)
C
C CALCULO DE LA CONTRIBUCION ASOCIATIVA AL COEFICIENTE DE FUGACIDAD
C
        IF(NCOD.NE.1)  CALL
     @FUGAASS(NC,nga2,MASSOC,ONUASOC,RNGRUP,XNOASS,FUGASS)

C
	END IF  !AVOID ASSOC 2
C
C CALCULO DE LA CONTRIBUCION  DE VOLUMEN LIBRE AL FACTOR DE COMPRESIBILIDAD
C
        IF(NCOD.NE.2) CALL
     @COMPRESFV(V,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ZFV)
C
C CALCULO DE LA CONTRIBUCION DISPERSIVA AL FACTOR DE COMPRESIBILIDAD
C
        IF(NCOD.NE.2)   CALL
     @COMPRESDISP(NC,NG,RNTOT,RN,NU,QAREA,H2,H4,H5,H6,ZDISP)
C
C
	IF(nga2.NE.0) THEN  !AVOID ASSOC 3
C
C CALCULO DE DERIVADAS DE LA FZA. DE ASOCIACION CON RESPECTO
C    AL VOLUMEN
C
        IF(NCOD.NE.2)CALL DERIROSSV(nga2,V,MASSOC,DELTASS,DROSSV,RNGRUP)
C
C CALCULO DE LAS DERIVADAS DE LAS FRACCIONES NO ASOCIADAS CON RESPECTO
C                       AL VOLUMEN
C
        IF(NCOD.NE.2)  CALL DERIXNOASSVL(nga2,MASSOC,DELTASS,XNOASS,
     @DROSSV,DXNOASSV, ONUASOC, nc,RNGRUP,V)
C
C CALCULO DE LAS DERIVADAS DE FASS(I) CON RESPECTO AL VOLUMEN
C
        IF(NCOD.NE.2)  CALL
     @DERIFASSV(nga2,MASSOC,DFAASSXNOASS,DXNOASSV,DFASSV)
C
C
C CALCULO DE LA CONTRIBUCION ASOCIATIVA AL FACTOR DE COMPRESIBILIDAD
C
        IF(NCOD.NE.2)  CALL COMPRESASS(nga2,V,RNTOT,RNGRUP,DFASSV,ZASS)
C
C
	 ELSE
	  CONTINUE
	END IF  !AVOID ASSOC 3
C
C CALCULO DEL FACTOR DE COMPRESIBILIDAD DEL SISTEMA
C
        IF(NCOD.NE.2)  CALL FCTRCOMPRES(ZFV,ZDISP,ZASS,ZSIS,IERRZSIS)
c	WRITE(4,*)"ZFV,ZDISP,ZASS,ZSIS=",ZFV,ZDISP,ZASS,ZSIS
C
C TERMINA LA EJECUCION SI EL FACTOR DE COMPRESIBILIDAD ES MENOR QUE CERO
C
        IF(IERRZSIS.EQ.1) GO TO 1000
C
C CALCULO DE LOS COEFICIENTES DE FUGACIDAD
C
        IF(NCOD.NE.1)  CALL
     @FUGASIS(NC,FUGREP,FUGDISP,FUGASS,ZSIS,FUGSIS,FFI,IERRFUG)
C	write(4,*)"rntot=", RNTOT
C	write(4,*)"rnGRUP=", RNGRUP(1)
C	write(4,*)"xnoass(1,1)=", XNOASS(1,1)
C	WRITE(4,*)"FUGREP,FUGDISP,FUGASS,FUGSIS,ZSIS",
C     @FUGREP(1),FUGDISP(1),FUGASS(1),FUGSIS(1),ZSIS
C	WRITE(4,*)FUGREP(2),FUGDISP(2),FUGASS(2),FUGSIS(2),ZSIS

C FINAL DE LA SUBRUTINA ASSOC-GC-EOS
C
 1000   END
C
C/////////////////////////////////////////////////////////////////
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
	SUBROUTINE SUMATORIA(NC,RN,RNTOT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION RN(NC)
C
	rNtot = sum(rN(:NC))
      !RNTOT=0.D0
	!DO I=1,NC
	  !RNTOT=RNTOT+RN(I)
      !END DO
C
	return
	ENDsubroutine
C///////////////////////////////////////////////////////////////////
C//////////////////////////////////////////////////////
C///////////////////////////////////////////////////////////////
	SUBROUTINE AHELMRESFV(RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ARESFV)
C
      Implicit real*8 (A-H, O-Z)

C-----------------------------------------------------------------------
C
C       ARESFV= CONTRIBUCION REPULSIVA (O DE VOLUMEN LIBRE) A LA ENERGIA
C                RESIDUAL DE HELMHOLTZ =
C
C                      R
C                 = ( A  / RT )         (las unidades son "moles")
C                              fv
C
C----------------------------------------------------------------
C
C
C  CALCULO DE ARESFV
C
C
        ARESFV=3D0*(RLAMBDA1*RLAMBDA2/RLAMBDA3)*(Y-1D0)
C
        ARESFV= ARESFV + (RLAMBDA2**3/RLAMBDA3**2)*(Y**2-Y-DLOG(Y))
C
        ARESFV= ARESFV + RNTOT*DLOG(Y)
C
	return
      ENDsubroutine
C///////////////////////////////////////////////////////////////////
	SUBROUTINE SOFT(NC,DC,TEMPCRIT,TEMP,D)
C
C CAlCULA LOS DIAMETROS DE ESFERAS RIGIDAS A LA TEMP DEL SISTEMA
C
        Implicit real*8 (A-H, O-Z)
        DIMENSION DC(NC), TEMPCRIT(NC), D(20)
C
        HELPDC = 0d0
C
C CALCULO DE LOS D(I)
C
        DO I=1,NC
                HELPDC=DEXP(-2D0*TEMPCRIT(I)/3D0/TEMP)
                D(I)=1.065655d0*DC(I)*(1D0-0.12D0*HELPDC)
        END DO
C
        END
C///////////////////////////////////////////////////////////////////
        SUBROUTINE LAMBDA(NC,RN,D,RLAMBDA1,RLAMBDA2,RLAMBDA3)
        Implicit real*8 (A-H, O-Z)
        DIMENSION D(20),RN(NC)
C
C CALCULO DE LOS LAMBDA NECESARIOS PARA LA CONTRIBUCION REPULSIVA
C
C
C  RLAMBDA3 ES PROPORCIONAL AL VOLUMEN OCUPADO POR LAS ESFERAS PROPIAMENTE
C   DICHAS, EXCLUYENDO EL ESPACIO INTERSTICIAL
C
	RLAMBDA1=0.D0
	RLAMBDA2=0.D0
	RLAMBDA3=0.D0
C
	DO I=1,NC

	  aux = RN(I)*D(I)
	  RLAMBDA1= aux + RLAMBDA1
	  RLAMBDA2= aux*d(i) + RLAMBDA2
	  RLAMBDA3= aux*D(I)**2  +  RLAMBDA3

	END DO
C
	return
	ENDsubroutine
C//////////////////////////////////////////////////////////////////////
        SUBROUTINE YCALC(V,RLAMBDA3,Y,IERRREP)
        Implicit real*8 (A-H, O-Z)
C
        HELPY = 0.D0
C
C  Y = VARIABLE INTERMEDIA EN LA EXPRESION DE MANSOORI Y LELAND (1972)
C----------------------------------------------------------------
        HELPY = 3.141592654D0*RLAMBDA3/(6.D0*V)
C	WRITE(6,*)"V=",V
C	WRITE(6,*)"RLAMDA3=",V

C        WRITE(6,*) ' '
        IF(HELPY.GT.(1.)) THEN
		WRITE(6,*) ' V ES DEMASIADO CHICO'
C	PAUSE
	   END IF
C        WRITE(6,*) ' '
        HELPY = 1 - HELPY
        Y = 1/HELPY
C
        IF(Y.LT.0) THEN
                        IERRREP = 1
                   ELSE
                        IERRREP = 0
        END IF
C
        END
C/////////////////////////////////////////////////////////////////////
C///////////////////////////////////////////////////////////////
C///////////////////////////////////////////////////////////////
        SUBROUTINE SURFSEG(NC,NG,NU,RN,QAREA,QTOT)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RN(NC)
C
        DIMENSION QAREA(10)
        DIMENSION NU(10,NC)
C
        HELPQ = 0.D0
C
C QTOT=NUMERO TOTAL DE SEGMENTOS DE AREA
C
C NU (I,J)= NUMERO DE GRUPOS DEL TIPO I EN EL COMPUESTO DE TIPO J
C NU (J,I)= NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO I
C
        QTOT=0.D0
C
        DO I=1,NC
                HELPQ=0.D0
                  DO J=1,NG
                        HELPQ=HELPQ + NU(J,I)*QAREA(J)
                  END DO
                QTOT=QTOT + RN(I)*HELPQ
        END DO
C
        END
C/////////////////////////////////////////////////////////////////////
        SUBROUTINE FRACAREA(NC,NG,NU,RN,QAREA,QTOT,TITA)
C
C CALCULA LAS FRACCIONES DE AREA DE TODOS LOS GRUPOS
C
C  TITA(J)= FRACCION DE AREA DEL GRUPO J
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RN(NC)
C
        DIMENSION QAREA(10)
        DIMENSION NU(10,NC)
        DIMENSION TITA(10)
C
        HELPFRAC = 0.D0
C
        DO J=1,NG
               HELPFRAC=0.D0
                DO I=1,NC
                  HELPFRAC=HELPFRAC + RN(I)*NU(J,I)
                END DO
               TITA(J)=HELPFRAC*QAREA(J)/QTOT
        END DO
C
        END
C//////////////////////////////////////////////////////////////////////////
	SUBROUTINE SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10)
	DIMENSION TSPLINE(10)
C
	DIF = 0.D0
	HELPAUX = 0.D0
C
C PRIMERA APROXIMACION A TSPLINE(J)
C
	DO J=1,NG

	  IF(GDOT(J).EQ.0.AND.GDDOT(J).EQ.0) TSPLINE(J) = 10000000.
	  IF(GDOT(J).EQ.0.AND.GDDOT(J).EQ.0) GO TO 127
C***********************************************************************

	    if(GDOT(J).EQ.0.AND.GDDOT(J).LT.0)then	!Fran 17/06/09	 !

	      HELPAUX = -0.75D0/GDDOT(J)							 !
	      TSPLINE(J) = exp(HELPAUX)*TSTR(J)
!	      write(*,*)'Ts[',J,']=',TSPLINE(J)					 !
!	      write(*,*)'Tstr[',J,']=',TSTR(J)
!	      write(*,*)'gp[',J,']=',GDOT(J)
	      go to 1270											 !

	    else
C***********************************************************************
	      HELPAUX=(0.25-1.)/GDOT(J) + 1.
	      TSPLINE(J)=HELPAUX*TSTR(J)
C	      write(*,*)'Ts[',J,'] =',TSPLINE(J)

1270	    end if

127	END DO
C
C CALCULO FINO DE TSPLINE
C
        DO J=1,NG
           IF(GDDOT(J).EQ.0) THEN
				CONTINUE
           ELSE
C***********************************************************************
		    IF(GDOT(J).LT.0) THEN		!Fran 17/06/09	 !
C				write(*,*)'gp[',J,']=',GDOT(J)
C***********************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 15               HELPAUX=(0.25 - 1. - GDDOT(J)*
     @DLOG(TSPLINE(J)/TSTR(J)))
                  HELPAUX=HELPAUX/GDOT(J) + 1.
                  HELPAUX=HELPAUX*TSTR(J)
                  DIF=DABS(HELPAUX-TSPLINE(J))
                  IF(DIF.LE.1.D-7) TSPLINE(J)=HELPAUX
                  IF(DIF.GT.1.D-7) TSPLINE(J)=HELPAUX
                  IF(DIF.GT.1.D-7) GO TO 15
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			end if
*************************************************************
           END IF
        END DO
C
        END
C//////////////////////////////////////////////////////////////
        SUBROUTINE CALGGG(NG,TEMP,TSTR,GSTR,GDOT,GDDOT,TSPLINE,GGG)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION TSTR(10), GSTR(10)
        DIMENSION GDOT(10), GDDOT(10)
        DIMENSION TSPLINE(10)
        DIMENSION GGG(10,10)
C
        ADIM = 0.D0
        EXPONENTE = 0.0D0
C
C GGG(I,J) = PARAMETRO DE ENERGIA ATRACTIVA (DE DISPERSION) ENTRE
C             EL GRUPO I Y EL GRUPO J
C
        DO J=1,NG
C
          ADIM=TEMP/TSTR(J)
C	    write(*,*)'T/Tstr[',J,'] =',ADIM
          IF(TEMP.LE.TSPLINE(J)) THEN

		  GGG(J,J)=  ADIM - 1.D0
            GGG(J,J)= GGG(J,J)*GDOT(J)
C		  write(*,*)'gp[',J,J,'](T/T* - 1) =',GGG(J,J)
            GGG(J,J)= GGG(J,J) + 1.D0
            GGG(J,J)= GGG(J,J) + GDDOT(J)*DLOG(ADIM)
            GGG(J,J)= GGG(J,J)*GSTR(J)
C		  write(*,*)'g[',J,J,'] = ',GGG(J,J)

	    ELSE

C	      write(*,*)'T spline =',J,TSPLINE(J)
            EXPONENTE= GDOT(J)*TSPLINE(J)/TSTR(J)
            EXPONENTE= EXPONENTE + GDDOT(J)
            EXPONENTE= (-4.D0)*EXPONENTE
            GGG(J,J)= 0.25D0*GSTR(J)
            GGG(J,J)= GGG(J,J)/(TEMP/TSPLINE(J))
     @                **EXPONENTE
C	      write(*,*)'g[',J,J,'] = ',GGG(J,J)

	    END IF

        END DO
C
        END
C////////////////////////////////////////////////////////////////
        SUBROUTINE CALCKIJ(NG,RKIJ,RKTIJ,TSTR,VERDKIJ,TEMP)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION TSTR(10)
        DIMENSION RKIJ(10,10), RKTIJ(10,10)
        DIMENSION VERDKIJ(10,10)
C
        HELPVERD = 0.D0
C
        DO I=1,NG
                DO J=1,NG
                        HELPVERD=TSTR(I)+TSTR(J)
                        HELPVERD=HELPVERD/2
c	write(*,*)TEMP
                        HELPVERD=DLOG(TEMP/HELPVERD)
                        HELPVERD=RKTIJ(I,J)*HELPVERD + 1.D0
                        VERDKIJ(I,J)=RKIJ(I,J)*HELPVERD
                END DO
        END DO
C
        END
C////////////////////////////////////////////////////////////////
        SUBROUTINE CALGGGIJ(NG,VERDKIJ,GGG)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION GGG(10,10)
        DIMENSION VERDKIJ(10,10)
C
C
C
        DO I=1,NG
C
                DO J=1,NG
C					   write(*,*)TEMP
C					   write(*,*)'g(',I,I,') =',GGG(I,I)
C					   write(*,*)'g(',j,j,') =',GGG(j,j)
                         GGG(I,J)=VERDKIJ(I,J)*DSQRT(GGG(I,I)*GGG(J,J))
                END DO
        END DO
C
        END
C////////////////////////////////////////////////////////////////
        SUBROUTINE DELTADISP(NG,GGG,DLTG)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION GGG(10,10)
        DIMENSION DLTG(10,10)
C
C DLTG(I,J) = DELTAgij  (NO ES SIMETRICO)
C
        DO J=1,NG
                DO I=1,NG
                        DLTG(J,I)=GGG(J,I)-GGG(I,I)
                END DO
        END DO
C
        END
C///////////////////////////////////////////////////////////
        SUBROUTINE TAUCALC(NG,TEMP,V,ALFAIJ,QTOT,DLTG,TAU)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        DIMENSION ALFAIJ(10,10)
        DIMENSION DLTG(10,10)
        DIMENSION TAU(10,10)
C
        HELPTAU = 0.D0
C
C CTE. UNIVERSAL DE LOS GASES
C
        RUNIV=82.05D0 !atm cm3/(K mol)
C
        DO J=1,NG
                DO I=1,NG
                        HELPTAU=ALFAIJ(J,I)*DLTG(J,I)*QTOT

                        HELPTAU= HELPTAU/(RUNIV*TEMP*V)
C					  write(*,*)HELPTAU
                        TAU(J,I)= DEXP(HELPTAU)
                END DO
        END DO
C
        END
C/////////////////////////////////////////////////////////
        SUBROUTINE AHELMRESDISP(NC,NG,RN,V,TEMP,NU,TITA,TAU,QTOT,
     @GGG,QAREA,ARESDISP)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RN(NC)
C
        DIMENSION NU(10,NC), QAREA(10)
        DIMENSION TITA(10), TAU(10,10)
        DIMENSION GGG(10,10)
        DIMENSION AUXIL2(20,20)
        DIMENSION SUM1(10),SUM2(10),SUM3(20)
C
C LIMPIEZA DE VARIABLES
C
        SUMAT = 0.D0
C
        DO I=1,NC+1
                 SUM3(I)=0.D0
C
                DO J=1,NC+1
                        AUXIL2(I,J)=0.D0
                END DO
        END DO
C
C
        DO I=1,NG
                SUM1(I)=0.D0
                SUM2(I)=0.D0
        END DO
C
C CTE. UNIVERSAL DE LOS GASES
C
        RUNIV=82.05D0 !atm cm3/(K mol)
C
        ZETA=10.D0 !NUMERO DE COORDINACION
C
C CALCULO DE LOS SUM1(J)
C
        DO J=1,NG
                SUM1(J)=0.D0
                DO L=1,NG
                        SUM1(J) = SUM1(J) + TITA(L)*TAU(L,J)
                END DO
        END DO
C
C
C
C AUXIL2 --------------------------------------
C
        DO J=1,NG
                DO K=1,NG
                        AUXIL2(K,J)=TITA(K)*(GGG(K,J)*QTOT*TAU(K,J))
                        AUXIL2(K,J)=AUXIL2(K,J)/(RUNIV*TEMP*V)
                END DO
        END DO
C
C CALCULO DE LOS SUM2(J)
C
        DO J=1,NG
                SUM2(J)=0.D0
                 DO K=1,NG
                        SUM2(J) = SUM2(J) + AUXIL2(K,J)
                 END DO
        END DO
C
C CALCULO DE LOS SUM3(I)
C
        DO I=1,NC
               SUM3(I)=0.D0
                DO J=1,NG
                  SUM3(I)= SUM3(I) + NU(J,I)*QAREA(J)*SUM2(J)/SUM1(J)
                END DO
        END DO
C
C SUMATORIA EN  EL NUMERO DE COMPONENTES
C
        SUMAT=0.D0
        DO I=1,NC
                 SUMAT= SUMAT + RN(I)*SUM3(I)
        END DO
C
C CALCULO DE ARESDISP(energ. resid. de dispersion)
        ARESDISP= (-1.D0)*(ZETA/2D0)*SUMAT
CC-----------------------------------------------------------------------
C
C       ARESDISP= CONTRIBUCION DE FZAS. DE DISPERSION A LA ENERGIA
C                RESIDUAL DE HELMHOLTZ =
C
C                      R
C                 = ( A  / RT )         (las unidades son "moles")
C                              disp.
C
C----------------------------------------------------------------
C
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
C                  DEL grupo L (EN GRADOS KELVIN)
C
C DELTASS(I,J,K,L) = FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL grupo J, Y EL SITIO K
C                  DEL grupo L (ADIMENSIONAL). NUESTRA DEFINICION
C                  NO COINCIDE CON LA DE RADOSZ.
C     Con esta definicion EL SISTEMA DE ECUACIONES (CUYAS VARIABLES SON
C       LAS FRACCIONES NO ASOCIADAS DE CADA SITIO) QUEDA EN FUNCION DE:
C
C            - FRACCIONES DE SITIOS NO ASOCIADOS y,
C            - FUERZA DE ASOCIACION (ADIMENSIONAL)
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
	               DeltAss(k,i,l,j) = rKAss(k,i,l,j)*
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
	SUBROUTINE FRACNOAS(NGA,MASSOC,DELTASS,XNOASS, nc, ONUASOC,V,
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

CCCCCCCCCCCCCCCCCC ALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC4(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP
     @,ierrassoc)

    !alcohol+AMINA PRIMARIA

C CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
C QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION

        IMPLICIT REAL*8 (A-H,O-Z)
C        implicit none
        double precision DELTASS(6,12,6,12)
        double precision XNOASS(6,12)
	  double precision RNGRUP(12)
	  double precision V,ROOT1(4)
	  integer JJ,IDAUT1,IDAUT2, ierrassoc
  	  DIMENSION ZC(3)
	  DOUBLE precision :: qa0,qa1,qa2,qa3,qa4,zq(4)
	  DOUBLE PRECISION DELCROSS_2

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
C   UN SITIO (ESTER)

C IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
C   DOS SITIOS (ALCOHOL)
C
C XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
C XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
C XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
C

C	PRIMEIRO CALCULAR XB2=X(2,2)
c	qa4=DELTASS(1,1,1,1)*(RNGRUP(2)/V)**2*DELTASS(1,2,2,2)**2-
c     @2D0*(RNGRUP(2)/V)**2*DELTASS(1,1,2,2)**2*DELTASS(1,2,2,2)



	DELTASS(1,1,2,2)=(DELTASS(1,1,2,1)*DELTASS(1,2,2,2))**0.5D0
	DELCROSS_2=(DELTASS(1,1,2,1)*DELTASS(1,2,2,2))

c	a3=DELTASS(1,1,2,1)*(RNGRUP(1)/V)*(DELTASS(1,1,2,2)-
c     @DELTASS(1,1,2,1))
c	a2=DELTASS(1,1,2,2)-DELTASS(1,1,2,1)-
c     @((RNGRUP(2)/V)+(RNGRUP(1)/V))*DELTASS(1,1,2,1)*
c     @DELTASS(1,1,2,2)
c	a1=-2D0*DELTASS(1,1,2,2)+DELTASS(1,1,2,1)
c	a0=DELTASS(1,1,2,2)


c		call cubicequation(a3,a2,a1,a0,zc)


      qa4=(RNGRUP(1)/V)**2D0*DELTASS(1,1,2,1)*
     @(DELTASS(1,1,2,1)*DELTASS(1,2,2,2)-DELCROSS_2)
	qa3=(RNGRUP(1)/V)*(2d0*DELTASS(1,1,2,1)*DELTASS(1,2,2,2)-
     @DELTASS(1,1,2,2)*DELTASS(1,1,2,1)-DELTASS(1,1,2,2)**2D0)
      qa2=DELTASS(1,2,2,2)-2d0*DELTASS(1,1,2,1)*DELTASS(1,2,2,2)*
     @(RNGRUP(1)/V)+(RNGRUP(1)/V)*DELTASS(1,1,2,2)**2D0-DELTASS(1,1,2,2)
     @-(RNGRUP(2)/V)*DELTASS(1,1,2,2)**2D0
	qa1=-2D0*DELTASS(1,2,2,2)+DELTASS(1,1,2,2)
	qa0=DELTASS(1,2,2,2)

c		call quarticequation(qa0,qa1,qa2,qa3,qa4,zq)


		ierrassoc=1

C		DO JJ=1,3
C			If (zc(JJ)<1d0.and.zc(JJ)>0d0)then  !numero real entre 0 e 1
C					root1(jj)=zc(jj)
C					XNOASS(1,1)=root1(jj)

		DO JJ=1,4
			If (zq(JJ)<1d0.and.zq(JJ)>0d0)then  !numero real entre 0 e 1
					root1(jj)=zq(jj)
					XNOASS(1,1)=root1(jj)

					If(rngrup(2)==0d0) then
						XNOASS(2,2)=1.D0/(1d0+RNGRUP(1)/V
     @					*XNOASS(1,1)*DELTASS(1,1,2,2))

					else
						XNOASS(2,2)=(1.D0-XNOASS(1,1)-
     @					(RNGRUP(1)/V)*DELTASS(1,1,2,1)*XNOASS(1,1)
     @					**2)/(DELTASS(1,1,2,2)*(RNGRUP(2)/V)
     @					*XNOASS(1,1))
					end if

					If (XNOASS(1,1)<=1d0.and.XNOASS(1,1)>=0d0)then
						XNOASS(1,2) = XNOASS(2,2)
						XNOASS(2,1) = XNOASS(1,1)
						ierrassoc=0
						goto 1132
					end if
				endif
     		END DO


1132		If (ierrassoc==1) then
			write(6,*)"No roots found"
			XNOASS(1,1)=1d0
			XNOASS(2,1)=1d0
			XNOASS(1,2)=1d0
			XNOASS(2,2)=1d0

		else

		end if

        END
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC3(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP) !!ester+acid
C
C CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
C QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
	  DIMENSION ZC(3)
	  DIMENSION RNGRUP(12)
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
C   UN SITIO (ESTER)

C IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
C   DOS SITIOS (ALCOHOL)
C
C XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
C XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
C XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
C

C	PRIMEIRO CALCULAR XB2=X(2,2)

	a3=DELTASS(1,1,1,1)*DELTASS(1,1,1,2)*(RNGRUP(1)/V)**2D0
	a2=(RNGRUP(1)/V)*(DELTASS(1,1,1,2)+DELTASS(1,1,1,1))
	a1=1D0+DELTASS(1,1,1,2)*(RNGRUP(2)/V-RNGRUP(1)/V)
	a0=-1d0

C	If (a3.gt.1e-13) then

c		call cubicequation(a3,a2,a1,a0,zc)

		do jj=1, 3
			If(zc(jj)>0d0.and.zc(jj)<=1d0)then
				root1=zc(jj)
			end if
		end do

		If (root1.ne.0d0) then
			XNOASS(1,1)=root1
			XNOASS(1,2) = 1.D0/(1D0+DELTASS(1,1,1,2)*(RNGRUP(1)/V)
     @			            *XNOASS(1,1))
		else
			write(6,*)"No roots found"
			XNOASS(1,1)=1d0
			XNOASS(1,2)=1d0
		end if
C	end if

c	write(2,*)"XNOASS(1,1)=",	XNOASS(1,1)
c	write(2,*)"XNOASS(1,2)=", XNOASS(1,2)
        END
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC2(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP !alcohol+acid
     @,ierrassoc)

C
C CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
C QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
C
c        IMPLICIT REAL*8 (A-H,O-Z)
	implicit none
        double precision DELTASS(6,12,6,12)
        double precision XNOASS(6,12)
	  double precision RNGRUP(12)
	  double precision V,ROOT1(4)
	  integer JJ,IDAUT1,IDAUT2, ierrassoc
	  DOUBLE precision :: qa0,qa1,qa2,qa3,qa4,zq(4)
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
C   UN SITIO (ESTER)

C IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
C   DOS SITIOS (ALCOHOL)
C
C XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
C XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
C XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
C

C	PRIMEIRO CALCULAR XB2=X(2,2)
	qa4=DELTASS(1,1,1,1)*(RNGRUP(2)/V)**2*DELTASS(1,2,2,2)**2-
     @2D0*(RNGRUP(2)/V)**2*DELTASS(1,1,2,2)**2*DELTASS(1,2,2,2)
	qa3=-DELTASS(1,1,2,2)*(RNGRUP(2)/V)*DELTASS(1,2,2,2)-2D0*
     @(RNGRUP(2)/V)*DELTASS(1,1,2,2)**2+2D0*DELTASS(1,1,1,1)*
     @(RNGRUP(2)/V)*DELTASS(1,2,2,2)
	qa2=DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-2D0*DELTASS(1,1,1,1)*
     @(RNGRUP(2)/V)*DELTASS(1,2,2,2)-(RNGRUP(1)/V)*DELTASS(1,1,2,2)**2+
     @2D0*(RNGRUP(2)/V)*DELTASS(1,1,2,2)**2
	qa1=-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2)
	qa0=DELTASS(1,1,1,1)

C		call quarticequation(qa0,qa1,qa2,qa3,qa4,zq)

		ierrassoc=1

		DO JJ=1,4
			If (zq(JJ)<1d0.and.zq(JJ)>0d0)then  !numero real entre 0 e 1
					root1(jj)=zq(jj)
					XNOASS(2,2)=root1(jj)

					If(rngrup(1)==0d0) then
						XNOASS(1,1)=1.D0/(1d0+2d0*RNGRUP(2)/V
     @					*XNOASS(2,2)*DELTASS(1,1,2,2))
					else
						XNOASS(1,1)=(1.D0-XNOASS(2,2)-
     @					(RNGRUP(2)/V)*DELTASS(1,2,2,2)*XNOASS(2,2)
     @					**2)/(DELTASS(1,1,2,2)*(RNGRUP(1)/V)
     @					*XNOASS(2,2))
					end if

					If (XNOASS(1,1)<=1d0.and.XNOASS(1,1)>=0d0)then
						XNOASS(1,2) = XNOASS(2,2)
						ierrassoc=0
						goto 1132
					end if
				endif
     		END DO


1132		If (ierrassoc==1) then
			write(6,*)"No roots found"
			XNOASS(1,1)=1d0
			XNOASS(2,2)=1d0
			XNOASS(1,2)=1d0
		else
c			write(6,*)"Roots found"
c	write(6,*)"XNOASS(1,1)=",	XNOASS(1,1)
c	write(6,*)"XNOASS(2,2)=", XNOASS(2,2)
c	write(6,*)"XNOASS(1,2)=", XNOASS(1,2)

		end if

c	write(2,*)"XNOASS(1,1)=",	XNOASS(1,1)
c	write(2,*)"XNOASS(2,2)=", XNOASS(2,2)
c	write(2,*)"XNOASS(1,2)=", XNOASS(1,2)


C
        END


C////////////////////////////////////////////////////////////////////////
        SUBROUTINE CROSSASOC(IDAUT1,IDAUT2,DELTASS,XNOASS,V,RNGRUP)!ester+alcohol
C
C CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
C QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
	  DIMENSION ZC(3)
	  DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT1 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE SO ASOCIA DE FORMA CRUZADA CON
C   UN SITIO (ESTER)

C IDAUT2 ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
C   DOS SITIOS (ALCOHOL)
C
C XNOASS(1,1)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 1
C XNOASS(1,2)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo 2
C XNOASS(2,2)= FRACCION NO ASOC. CORRESP. AL SITIO 2 DEL grupo 2
C

C	PRIMEIRO CALCULAR XB2=X(2,2)

	a3=DELTASS(1,1,2,2)*DELTASS(1,2,2,2)*(RNGRUP(2)/V)**2D0
	a2=DELTASS(2,2,1,1)*DELTASS(1,2,2,2)*(RNGRUP(2)/V)*(RNGRUP(1)/V)+
     @DELTASS(1,1,2,2)*(RNGRUP(2)/V)+DELTASS(1,2,2,2)*(RNGRUP(2)/V)
	a1=DELTASS(2,2,1,1)*(RNGRUP(1)/V)-DELTASS(1,1,2,2)*(RNGRUP(2)/V)+1
	a0=-1d0

C	If (a3.gt.1e-13) then

c		call cubicequation(a3,a2,a1,a0,zc)

		do jj=1, 3
			If(zc(jj)>0d0.and.zc(jj)<=1d0)then
				root1=zc(jj)
			end if
		end do

		If (root1.ne.0d0) then
			XNOASS(2,2)=root1
			XNOASS(1,1) = 1.D0/(1D0+DELTASS(1,1,2,2)*(RNGRUP(2)/V)
     @						*XNOASS(2,2))
			XNOASS(1,2) = 1.D0/(1D0+DELTASS(1,2,2,2)*(RNGRUP(2)/V)
     @			            *XNOASS(2,2))
		else
			write(6,*)"No roots found"
			XNOASS(1,1)=1d0
			XNOASS(2,2)=1d0
			XNOASS(1,2)=1d0
		end if
C	end if

c	write(2,*)"XNOASS(1,1)=",	XNOASS(1,1)
c	write(2,*)"XNOASS(2,2)=", XNOASS(2,2)
c	write(2,*)"XNOASS(1,2)=", XNOASS(1,2)
        END
C////////////////////////////////////////////////////////////////////////

C A ESTA SUBRUTINA ES NECESARIO AMPLIARLA PARA INCLUIR LA POSIBILIDAD DE
C AUTOASOCIACION CON UN NUMERO DE SITIOS DISTINTO DE DOS, Y LA POSIBILIDAD
C DE ASOCIACION CRUZADA.
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE AUTOASOC(IDAUT,DELTASS,XNOASS,V,RNGRUP)
C
C CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
C QUE ASOCIA, QUE TIENE DOS SITIOS DE ASOCIACION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
	  DIMENSION RNGRUP(12)
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
C   DOS SITIOS
C
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo
C                                                       IDAUT
C
        IF(DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

        AUXO2 = 1D0 +  4D0 *DELTASS(1,IDAUT,2,IDAUT)* RNGRUP(IDAUT)/V

C
		If ((AUXO2.gt.0d0).and.(RNGRUP(IDAUT).gt.0d0)) then
        AUXO2 =  (-1D0) + DSQRT(AUXO2)
C
        AUXO2 = AUXO2 / (2D0 * DELTASS(1,IDAUT,2,IDAUT)*RNGRUP(IDAUT)/V)
        XNOASS(1,IDAUT) = AUXO2
		else
		 XNOASS(1,IDAUT) = 1.
		end if
         ELSE
C
C
                XNOASS(1,IDAUT) = 1.
C
C       ----------------
C
        END IF
C
C       ----------------
C

          IF(DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

			XNOASS(2,IDAUT) = XNOASS(1,IDAUT)

		END IF
c	       write(2,*)"XNOASS(1,",IDAUT,")=", XNOASS(1,IDAUT)
c            write(2,*)"XNOASS(2,",IDAUT,")=", XNOASS(2,IDAUT)

C
        END
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE autoasocumsitio(IDAUT,DELTASS,XNOASS,V,RNGRUP)
C
C CALCULA LA FRACCION NO ASOCIADA SOLO EN EL CASO DE UN grupo
C QUE ASOCIA, QUE TIENE uno SITIO DE ASOCIACION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
	  DIMENSION RNGRUP(12)
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT ES IGUAL AL NUMERO ASIGNADO AL grupo QUE AUTOASOCIA CON
C   uno SITIO
C
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL grupo
C                                                       IDAUT
C
        IF(DELTASS(1,IDAUT,1,IDAUT).NE.0.and.rngrup(1).gt.0d0) THEN

        AUXO2 = 1D0 + 4D0 * DELTASS(1,IDAUT,1,IDAUT)* RNGRUP(IDAUT)/V
C
        AUXO2 =(-1D0) + DSQRT(AUXO2)
C
        AUXO2 = AUXO2 / (2D0* DELTASS(1,IDAUT,1,IDAUT)*RNGRUP(IDAUT)/V)
        XNOASS(1,IDAUT) = AUXO2
		ELSE
                XNOASS(1,IDAUT) = 1D0
C
C       -----

		endif
C
c       write(2,*)"XNOASS(1,",IDAUT,")=", XNOASS(1,IDAUT)
	   END
C/////////////////////////////////////////////////////////////////////////


        SUBROUTINE AHELMRESASS(NGA,MASSOC,RNGRUP,XNOASS,ARESASS,
     @FASS,FAASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION RNGRUP(12)
        DIMENSION XNOASS(6,12)
        DIMENSION FASS(12)
        DIMENSION FAASS(6,12)
C
C CALCULO DE LOS FAASS(K,I)
C
        DO I=1,NGA
                IF(MASSOC(I).EQ.0) GO TO 881
                DO K=1,MASSOC(I)
C
                        AUXO3 = 0.D0
					  If(XNOASS(K,I).gt.1d-13)then
						AUXO3 = DLOG( XNOASS(K,I) )
						AUXO3 = AUXO3 - XNOASS(K,I)/2.D0
						FAASS(K,I) = AUXO3
					  end if
C
                END DO
 881            CONTINUE
        END DO
C
C
C CALCULO DE LAS SUMATORIAS SOBRE LOS SITIOS DE ASOCIACION PARA CADA
C  GRUPO
C
C FASS(I)= {SUMATORIA SOBRE SITIOS DE [ FAASS(K,I) ] }  +
C
C                                + 0.5 MASSOC(I)
C
C SI FASS(I) SE MULTIPLICA POR EL NRO. DE MOLES DEL grupo(I)
C
C               SE OBTIENE LA CONTRIBUCION DE ESE grupo a ARESASS
C
        DO I=1,NGA
C          -----------------------------------------------
           IF(MASSOC(I).NE.0) THEN
C
               AUXO4=0.D0
C
                DO K=1,MASSOC(I)
                  AUXO4 = AUXO4 + FAASS(K,I)
                END DO
C
               FASS(I)= AUXO4 + 0.5D0 * MASSOC(I)
C
           END IF
C          -----------------------------------------------
C
        END DO
C
C CALCULO FINAL DE LA ENERGIA RESIDUAL DE HELMHOLTZ POR ASOCIACION
C
        SUME=0.D0
C
        DO I=1,NGA
          SUME = RNGRUP(I) * FASS(I) + SUME
        END DO
C
        ARESASS = SUME
C-----------------------------------------------------------------------
C
C       ARESASS = CONTRIBUCION DE LA ASOCIACION A LA ENERGIA
C                RESIDUAL DE HELMHOLTZ =
C
C                      R
C                 = ( A  / RT )         (las unidades son "moles")
C                              ass
C
C----------------------------------------------------------------
C
        END
C///////////////////////////////////////////////////////////////////////
        SUBROUTINE AHELMRESTOT(ARESFV,ARESDISP,ARESASS,ARESTOT)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C CALCULO DE LA ENERGIA RESIDUAL TOTAL
C
        ARESTOT = ARESFV + ARESDISP + ARESASS
C
        END
C//////////////////////////////////////////////////////////////////
        SUBROUTINE DERIVFV(NC,V,D,Y,DL1N,DL2N,DL3N,DYN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION D(20),DL1N(20),DL2N(20),DL3N(20),DYN(20)
C
C CALCULA DERIVADAS SEGUN DETALLE
C
C DL1N(I) = DERIVADA DE RLAMBDA1 CON RESPECTO AL NUMERO DE MOLES DEL
C                                       COMPONENTE I
C
C DL2N(I) = DERIVADA DE RLAMBDA2 CON RESPECTO AL NUMERO DE MOLES DEL
C                                       COMPONENTE I
C
C
C DL3N(I) = DERIVADA DE RLAMBDA3 CON RESPECTO AL NUMERO DE MOLES DEL
C                                       COMPONENTE I
C
C DYN(I)  = DERIVADA DE LA VARIABLE Y CON RESPECTO AL NRO. DE MOLES
C                                       DEL COMPONENTE I
C
C
        DO I=1,NC
                DL1N(I) = D(I)
                DL2N(I) = D(I)**2
                DL3N(I) = D(I)**3
                DYN(I)  = Y**2 * 3.141592654D0 * D(I)**3 / (6.D0*V)
        END DO
C
        END
C
C//////////////////////////////////////////////////////////////////
        SUBROUTINE FUGAREP(NC,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,
     @DL1N,DL2N,DL3N,DYN,FUGREP)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION DL1N(20),DL2N(20),DL3N(20),DYN(20)
        DIMENSION FUGREP(20)
C
        BLOQUE1 =0.
        BLOQUE2 =0.
        BLOQUE3 =0.
        BLOQUE4 =0.
        HELPFUR =0.
C
C CALCULO DE LA CONTRIBUCION REPULSIVA A LOS COEFICIENTES DE FUGACIDAD
C
C FUGREP(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD REPULSIVO
C                DEL COMPONENTE I
C
        BLOQUE1= 3D0*(Y-1)*RLAMBDA1*RLAMBDA2/RLAMBDA3
C
        BLOQUE2= RLAMBDA2**3/RLAMBDA3**2
        BLOQUE3= BLOQUE2*( Y**2 - Y - DLOG(Y)  )
C
        BLOQUE4= BLOQUE2*( 2.D0*Y - 1D0 - 1D0/Y  )
C
C
C
        DO I=1,NC
C
         HELPFUR = DL1N(I)/RLAMBDA1+DL2N(I)/RLAMBDA2-DL3N(I)/RLAMBDA3
         HELPFUR = HELPFUR + DYN(I)/(Y-1D0)
         HELPFUR = HELPFUR*BLOQUE1
         HELPFUR = HELPFUR + BLOQUE3 *
     @( 3.D0 * DL2N(I)/RLAMBDA2 - 2.D0 * DL3N(I)/RLAMBDA3)
         FUGREP(I)= HELPFUR + BLOQUE4*DYN(I) + DLOG(Y) +
     @(RNTOT/Y)*DYN(I)
C
        END DO
C
        END
C/////////////////////////////////////////////////////////
C/////////////////////////////////////////////////////////
        SUBROUTINE HCALC(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,
     @TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NOTAR QUE EN EL CASO DE NU EL PRIMER ARGUMENTO CORRESPONDE A GRUPO
C   Y EL SEGUNDO A COMPONENTE
C
C EN EL CASO DE PS,H3, Y H7 : PRIMER ARG --- COMPONENTE
C                             SEGUNDO ARG--- GRUPO
C
        DIMENSION PS(20,10), RMS(20), H2(10)
        DIMENSION H3(20,10), H4(10), H5(10)
        DIMENSION H6(10), H7(20,10)
C
        DIMENSION NU(10,NC), QAREA(10)
        DIMENSION ALFAIJ(10,10)
        DIMENSION TITA(10), TAU(10,10)
        DIMENSION GGG(10,10)
        DIMENSION DLTG(10,10)
C
        AUXH2 = 0.D0
        AUXH3 = 0.D0
        AUXH4 = 0.D0
        AUXH5 = 0.D0
        AUXH6 = 0.D0
        AUXH7 = 0.D0
        AUXRMS = 0.D0
C
        RUNIV=82.05D0 !atm cm3/(K mol)
C
C--------CALCULO DE PS
C
        DO I=1,NC
                DO J=1,NG
                        PS(I,J)= NU(J,I) * QAREA(J)
                END DO
        END DO
C
C--------CALCULO DE RMS
C
        DO I=1,NC
                AUXRMS=0.D0
                 DO J=1,NG
                        AUXRMS= PS(I,J) + AUXRMS
                 END DO
                RMS(I)=AUXRMS
        END DO
C
C--------CALCULO DE H2
C
C
        DO J=1,NG
                AUXH2=0.D0
                 DO K=1,NG
C
          AUXH2 = TITA(K)*TAU(K,J)*GGG(K,J)*QTOT/(V*RUNIV*TEMP)+ AUXH2
C
                 END DO
                H2(J)= AUXH2
        END DO
C
C--------CALCULO DE H3
C
C
        DO I=1,NC
                DO J=1,NG
C
                      AUXH3 =0.D0
                        DO K=1,NG
C
         AUXH3 = PS(I,K)*TAU(K,J)*GGG(K,J)*QTOT/(V*RUNIV*TEMP) + AUXH3
C
                        END DO
                      H3(I,J) = AUXH3
C
                END DO
C
        END DO
C
C
C--------CALCULO DE H4
C
C
        DO J=1,NG
C
               AUXH4 = 0.D0
C
                DO K=1,NG
C
                  AUXH4 = AUXH4 + TITA(K) * TAU(K,J)
C
                END DO
C
               H4(J)= AUXH4
C
        END DO
C
C
C--------CALCULO DE H5
C
C
        DO J=1,NG
C
             AUXH5 = 0.D0
C
                DO K=1,NG
C
                 AUXH5 = AUXH5 +  TITA(K)*TAU(K,J)*(GGG(K,J)*
C
     @QTOT/(V*RUNIV*TEMP))*ALFAIJ(K,J)*DLTG(K,J)*QTOT/(V*RUNIV*TEMP)
C
                END DO
C
             H5(J)= AUXH5
C
        END DO
C
C
C
C--------CALCULO DE H6
C
C(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,
C     @TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
C
C
        DO J=1,NG
C
             AUXH6 = 0.D0
C
                DO K=1,NG
C
                AUXH6 = AUXH6 + TITA(K)*TAU(K,J)*ALFAIJ(K,J)
     @*DLTG(K,J)*QTOT/(V*RUNIV*TEMP)
C
                END DO
C
             H6(J) = AUXH6
C
        END DO
C
C
C
C
C
C--------CALCULO DE H7
C
C(NC,NG,NU,QAREA,ALFAIJ,QTOT,TITA,GGG,DLTG,TAU,
C     @TEMP,V,PS,RMS,H2,H3,H4,H5,H6,H7)
C
        DO I=1,NC
                DO J=1,NG
C
                  AUXH7 = 0.D0
C
                        DO K=1,NG
C
                           AUXH7 = AUXH7 + PS(I,K)*TAU(K,J)
C
                        END DO
C
                  H7(I,J) = AUXH7
C
                END DO
C
        END DO
C
C
        END
C///////////////////////////////////////////////////////////
        SUBROUTINE FUGADISP(NC,NG,TITA,PS,RMS,H2,H3,H4,H5,H6,H7,FUGDISP)
C
C calculo de la contribucion dispersiva al coeficiente de fugacidad
C
C FUGDISP(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD DISPERSIVO
C                DEL COMPONENTE I
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION PS(20,10), RMS(20), H2(10)
        DIMENSION H3(20,10), H4(10), H5(10)
        DIMENSION H6(10), H7(20,10)
C
        DIMENSION TITA(10)
C
        DIMENSION FUGDISP(20)
C
        SUM24 = 0.D0
        SUM2467 = 0.D0
        SUM345 = 0.D0
C
        ZETA=10.D0 !NUMERO DE COORDINACION
C
C CALCULO DE SUMATORIAS SOBRE EL NUMERO DE GRUPOS (p' c/componente)
C
        DO I=1,NC
C
                SUM24 = 0.D0
                SUM345 = 0.D0
                SUM2467 = 0.D0
C
                  DO J=1,NG
C
                        SUM24 = SUM24 + PS(I,J) * H2(J) / H4(J)
C
                        SUM345 = SUM345 + TITA(J) *
C
     @( H3(I,J) + RMS(I) * H5(J) ) / H4(J)
C
                        SUM2467 = SUM2467 + TITA(J) * H2(J) *
C
     @( H7(I,J) - H4(J) * RMS(I) + H6(J) * RMS(I) ) / H4(J)**2
C
                  END DO
C
C
                FUGDISP(I) =(-1.D0) * (ZETA/2D0) *
C
     @( SUM24 + SUM345 - SUM2467 )
C
        END DO
C
        END
C//////////////////////////////////////////////////////////////////////
C////////////////////////////////////////////////////////////////
        SUBROUTINE DERIROSSN(NC,NGA,MASSOC,ONUASOC,RNGRUP,
     @DROSSN,V)

        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION ONUASOC(12,NC)
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION DROSSN(12,20)
C
C DELTASS(I,J,K,L) = FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L (ADIMENSIONAL). NUESTRA DEFINICION
C                  NO COINCIDE CON LA DE RADOSZ
C
C DDLTSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M  (ESTO ES LO QUE SE CALCULA EN ESTA
C                  SUBRUTINA)
C
C
        DO J=1,NGA
              IF(MASSOC(J).EQ.0) GO TO 201
C
                 DO M=1,NC
C
					DROSSN(J,M) =  ONUASOC(J,M)/V
C
                 END DO
 201          CONTINUE
        END DO
C
C
        END
C///////////////////////////////////////////////////////////////////
	SUBROUTINE DERIXNOASS(NC,NGA,MASSOC,DELTASS,XNOASS,DROSSN,
     &                      DXNOASSN, ONUASOC,RNGRUP,V)

C
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION MASSOC(12), sumsi(12)
	DIMENSION DELTASS(6,12,6,12)
	DIMENSION XNOASS(6,12)
	DIMENSION DROSSN(12,20)
	DIMENSION DXNOASSN(6,12,20)
	DIMENSION ONUASOC(12,NC)
	DIMENSION RNGRUP(12)
	dimension psin(12)
      real*8  m_lambda(nga*2,nga*2), m_lambda1(nga*2,nga*2)
      integer order,w,s
C	real d1, psin
	dimension indx(20)
!
!	Para información sobre el código anterior a las modificaciones de Alfonsina,
!	ver subroutina "FracNoAsoc".
!	Francisco.

CCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC GENERAL
c	cCalcula los elementos de la matriz deltapq
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
c Calcula los elementos de la matriz [Yp]
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
c colectando las derivadas en su correspondiente subíndice
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
CCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	return
      ENDsubroutine

C////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOC3(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
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
	DXNOASSN(1,1,M)= (-t3*t5-t2*DELTASS(1,1,1,1)*t5-2D0*t9*t5*t4*
     #DELTASS(1,1,1,2)*t2-t4*t16*DELTASS(1,1,1,2)+t3*t4)/(1D0+2D0*t21*t4
     #+2D0*t9*t4+3D0*t26*DELTASS(1,1,1,1)*t5*DELTASS(1,1,1,2)+t15*
     #DELTASS(1,1,1,2)-t21)



	DXNOASSN(1,2,M)=-1D0*XNOASS(1,2)*DELTASS(1,1,1,2)*
     #(DROSSN(1,M)*XNOASS(1,1)+
     # (RNGRUP(1)/V)*DXNOASSN(1,1,M))
     #/(1D0+(RNGRUP(1)/V)*DELTASS(1,1,1,2)*XNOASS(1,1))

c    WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)
c	WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)


	END DO



      END
C////////////////////////////////////////////////////////////////////////
       SUBROUTINE DCROSSASOC4(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR




C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL
C NRO. DE MOLES DEL COMPUESTO M

C
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
	DXNOASSN(2,2,M) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-
     #(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+
     #t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*DELTASS(1,2,2,2
     #))*
     #t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*
     #(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4
     #5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


	DXNOASSN(1,1,M)=(-DROSSN(1,M)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2
     #*XNOASS(1,1)*DROSSN(2,M)*DELTASS(1,1,2,2)*XNOASS(2,2)-2*
     #XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSN(2,2,M))
     #/(1+2*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2*(RNGRUP(2)/V)
     #*DELTASS(1,1,2,2)*XNOASS(2,2))



	DXNOASSN(1,2,M) = DXNOASSN(2,2,M)
c      WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)


	END DO

c	DO M=1,2
c		WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)
c		WRITE(2,*)"DXNOASSN(2,2,",M,")=",DXNOASSN(2,2,M)
c	END DO

      END
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOC2(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR




C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL
C NRO. DE MOLES DEL COMPUESTO M

C
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
	DXNOASSN(2,2,M) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-
     #(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+
     #t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*DELTASS(1,2,2,2
     #))*
     #t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*
     #(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4
     #5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


	DXNOASSN(1,1,M)=(-DROSSN(1,M)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2D0
     #*XNOASS(1,1)*DROSSN(2,M)*DELTASS(1,1,2,2)*XNOASS(2,2)-2*
     #XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSN(2,2,M))
     #/(1+2*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2*(RNGRUP(2)/V)
     #*DELTASS(1,1,2,2)*XNOASS(2,2))



	DXNOASSN(1,2,M) = DXNOASSN(2,2,M)
c      WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)


	END DO

c	DO M=1,2
c		WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)
c		WRITE(2,*)"DXNOASSN(2,2,",M,")=",DXNOASSN(2,2,M)
c	END DO

      END



C////////////////////////////////////////////////////////////////////////
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOC(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR




C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL
C NRO. DE MOLES DEL COMPUESTO M

C
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

	DXNOASSN(2,2,M) = (-t6-t7*t5-2*t9*t5*t4*DELTASS(1,2,2,2)*t2-t4*t16
     #*DELTASS(1,1,2,2)-t16*DELTASS(1,1,2,2)*t20*t5-t23*t6+t7*t4)/(1+2*
     #t27+2*t9*t4+3*t31*DELTASS(1,1,2,2)*t5*DELTASS(1,2,2,2)+t23+2*t23*t
     #27-t9)

	DXNOASSN(1,1,M) =-1D0*XNOASS(1,1)*DELTASS(1,1,2,2)*(RNGRUP(2)/V*
     #DXNOASSN(2,2,M)+XNOASS(2,2)*DROSSN(2,M))/(1D0+RNGRUP(2)/V*
     #XNOASS(2,2)*DELTASS(1,1,2,2))

	DXNOASSN(1,2,M) = -1D0*XNOASS(1,2)*DELTASS(1,2,2,2)*(RNGRUP(2)/V*
     #DXNOASSN(2,2,M)+XNOASS(2,2)*DROSSN(2,M))/(1D0+RNGRUP(2)/V*
     #XNOASS(2,2)*DELTASS(1,2,2,2))

c      WRITE(2,*)"DXNOASSN(1,1,",M,")=",DXNOASSN(1,1,M)


	END DO

c	DO M=1,2
c		WRITE(2,*)"DXNOASSN(1,2,",M,")=",DXNOASSN(1,2,M)
c		WRITE(2,*)"DXNOASSN(2,2,",M,")=",DXNOASSN(2,2,M)
c	END DO

      END
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE DAUTOASOC(IDAUT,NC,
     @DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
C   QUE AUTOASOCIA (CON DOS SITIOS).
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DDLTSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR

	IF (DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

	  BLOQ =0.D0
        HELPI =0.D0
C
C
        BLOQ = 1D0+2D0*DELTASS(1,IDAUT,2,IDAUT)*XNOASS(1,IDAUT)*
     @RNGRUP(IDAUT)/V
C
        BLOQ = 1D0 / BLOQ
C
C EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
C                               respecto al DELTASS(1,IDAUT,2,IDAUT)
C
C
        BLOQ = (-1D0)  *  (XNOASS(1,IDAUT))**2D0 * BLOQ
C
C CALCULO DE DERIVADAS RESPECTO AL NUMERO DE MOLES
C
        DO M = 1,NC

                HELPI =BLOQ *DELTASS(1,IDAUT,2,IDAUT)*DROSSN(IDAUT,M)
C
                DXNOASSN(1,IDAUT,M) = HELPI
C
        END DO
C
C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL
C NRO. DE MOLES DEL COMPUESTO M
C
        DO M=1,NC
C
                DXNOASSN(2,IDAUT,M) = DXNOASSN(1,IDAUT,M)

c			WRITE(2,*)"DXNOASSN(1,",IDAUT,M,")=",DXNOASSN(1,IDAUT,M)
c			WRITE(2,*)"DXNOASSN(2,",IDAUT,M,")=",DXNOASSN(2,IDAUT,M)

C
        END DO
	END IF
C
C
        END
C/////////////////////////////////////////////////////////////////////////
        SUBROUTINE DAUTOASOCumsitio(IDAUT,NC,
     @DELTASS,XNOASS,DROSSN,DXNOASSN,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
C   QUE AUTOASOCIA (CON DOS SITIOS).
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DDLTSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR
C        if (DELTASS(1,IDAUT,1,IDAUT).NE.0)THEN
        IF(DELTASS(1,IDAUT,1,IDAUT).NE.0) THEN

	  BLOQ =0.D0
        HELPI =0.D0
C
C
        BLOQ = 1D0+2D0*DELTASS(1,IDAUT,1,IDAUT)*XNOASS(1,IDAUT)*
     @RNGRUP(IDAUT)/V
C
        BLOQ = 1D0 / BLOQ
C
C EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
C                               respecto al DELTASS(1,IDAUT,2,IDAUT)
C
C
        BLOQ = (-1D0)  *  (XNOASS(1,IDAUT))**2D0 * BLOQ
C
C CALCULO DE DERIVADAS RESPECTO AL NUMERO DE MOLES
C
        DO M = 1,NC

                HELPI =BLOQ *DELTASS(1,IDAUT,1,IDAUT)*DROSSN(IDAUT,M)
C
                DXNOASSN(1,IDAUT,M) = HELPI
c			  WRITE(2,*)"DXNOASSN(1,",IDAUT,M,")=",DXNOASSN(1,IDAUT,M)
C
        END DO
	endif

        END
C/////////////////////////////////////////////////////////////////////////
C/////////////////////////////////////////////////////////
        SUBROUTINE DERIFASSN(NC,NGA,MASSOC,DFAASSXNOASS,DXNOASSN,DFASSN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION DFAASSXNOASS(6,12)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION DFASSN(12,20)
C
        HELX = 0.0D0
C
C DFASSN(I,J) = DERIVADA DE FASS DEL GRUPO I CON RESPECTO AL NRO. DE
C               MOLES DEL COMPUESTO J
C
        DO I=1,NGA   !grupos
                DO M=1,NC  !componentes
                        IF(MASSOC(I).NE.0) THEN
C
                                HELX = 0.D0
C
                                  DO K=1,MASSOC(I)
C
         HELX = HELX +  DFAASSXNOASS(K,I) * DXNOASSN(K,I,M)
C
                                  END DO
C
                                DFASSN(I,M) = HELX
                          ELSE
                            CONTINUE
                        END IF
                END DO
        END DO
C
        END
C///////////////////////////////////////////////////////////////////
        SUBROUTINE FUGAASS(NC,NGA,MASSOC,ONUASOC,RNGRUP,XNOASS,FUGASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DOUBLE PRECISION XNOASS(6,12)
        DIMENSION ONUASOC(12,NC)
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
        DIMENSION MASSOC(12)
        DIMENSION FUGASS(20)
	  DOUBLE PRECISION SUMAJ
C
        SUMAJ = 0d0
C
C FUGASS(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD
C               ASOCIATIVO PARA EL COMPONENTE I
C
C
C DFASSN(I,M) = DERIVADA DE FASS DEL grupo I CON RESPECTO AL NRO. DE
C               MOLES DEL COMPUESTO M
C


        DO M=1,NC
C
                       SUMAJ = 0d0
C
                        DO I=1,NGA !grupos

	                        IF(MASSOC(I).NE.0) THEN
C

C
                                  DO K=1,MASSOC(I)

						  If(XNOASS(K,I).gt.1d-13)then

C
       	   SUMAJ = SUMAJ + ONUASOC(I,M) * dlog(XNOASS(K,I))

						  end if
C
                                  END DO
C
                          ELSE
                            CONTINUE
                        END IF

C

C
                        END DO
C
                       FUGASS(M) = SUMAJ
C
        END DO
C
        END
C///////////////////////////////////////////////////////////////////
        SUBROUTINE COMPRESFV(V,RNTOT,RLAMBDA1,RLAMBDA2,RLAMBDA3,Y,ZFV)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C ZFV = CONTRIBUCION REPULSIVA AL FACTOR DE COMPRESIBILIDAD
C
        ZFV = 3D0*RLAMBDA1*RLAMBDA2/RLAMBDA3 +
     @      (RLAMBDA2**3/RLAMBDA3**2)*( 2.D0*Y - 1.D0 - 1D0/Y) + RNTOT/Y
C
        ZFV = ZFV * (-1.D0) * ( V / RNTOT)
C
        ZFV = ZFV * (-1.D0) * Y**2 * 3.141592654D0* RLAMBDA3/(6.D0*V**2)
C
        END
C/////////////////////////////////////////////////////////////////////
C/////////////////////////////////////////////////////////////////////
        SUBROUTINE COMPRESDISP(NC,NG,RNTOT,RN,NU,QAREA,H2,H4,H5,H6,
     @ZDISP)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RN(NC)
        DIMENSION NU(10,NC), QAREA(10)
        DIMENSION H2(10)
        DIMENSION H4(10), H5(10)
        DIMENSION H6(10)
        DIMENSION RNURAYA(10)
C
        SSUM =0.D0
        SUMM =0.D0
C
        DO I=1,NG
                RNURAYA(I)=0.D0
        END DO
C
C
C ZDISP = CONTRIBUCION DISPERSIVA AL FACTOR DE COMPRESIBILIDAD
C
C RNURAYA(J) ES UN ARREGLO AUXILIAR
C
C CALCULO DE  RNURAYA(J)
C
        DO J = 1,NG
C
              SSUM = 0.D0
C
                DO L=1,NC
                        SSUM = SSUM + RN(L) * NU(J,L)
                END DO
C
              RNURAYA(J) = SSUM
C
        END DO
C
C NU (J,L)= NUMERO DE GRUPOS DEL TIPO J EN EL COMPUESTO DE TIPO L
C
        SUMM = 0.D0
C
        DO J = 1,NG
C
           SUMM  = RNURAYA(J) * QAREA(J) *
C
     @ ( H5(J) + H2(J) - H2(J)*H6(J)/H4(J) ) /
C
     @  H4(J) + SUMM
C
        END DO
C
C
        ZETA=10.D0 !NUMERO DE COORDINACION
C
        ZDISP  =  (ZETA/2.D0) * ( -1.D0/ RNTOT ) * SUMM
C
        END
C///////////////////////////////////////////////////////////////////////
C///////////////////////////////////////////////////////////////////////////
        SUBROUTINE DERIROSSV(NGA,V,MASSOC,DELTASS,DROSSV,RNGRUP)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION DROSSV(12)
        DIMENSION RNGRUP(12)

C
C DROSSV(I,J,K,L) = DERIVADA DEL DELTA DE ASOCIACION ENTRE EL SITIO I
C                       DEL grupo J Y EL SITIO K DEL grupo L CON RESPECTO
C                               AL VOLUMEN
C
C
        DO J=1,NGA
              IF(MASSOC(J).EQ.0) GO TO 201
                        DROSSV(J) = (-1D0*RNGRUP(J))/V**2D0
c						write(2,*)"DROSSV(",J,")=",DROSSV(J)
 201          CONTINUE
        END DO
C
        END
C///////////////////////////////////////////////////////////////////
        SUBROUTINE DERIFASSV(NGA,MASSOC,DFAASSXNOASS,DXNOASSV,DFASSV)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION DFAASSXNOASS(6,12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION DFASSV(12)
C
        SUMV = 0.D0
C
C DFASSV(I) = DERIVADA DE FASS(I) CON RESPECTO AL VOLUMEN
C
        DO I =1,NGA
C
                SUMV =0.D0
C               -----------------------------------------
                IF(MASSOC(I).NE.0) THEN
C               -----------------------------------------
C
                          DO K=1,MASSOC(I)

                           SUMV = SUMV + DFAASSXNOASS(K,I)

                          END DO
C               -----------------------------------------
C
				ELSE
                  CONTINUE
                END IF
C               -----------------------------------------
C
                DFASSV(I) = SUMV
        END DO
C
C
        END
C/////////////////////////////////////////////////////////////////////////
        SUBROUTINE COMPRESASS(NGA,V,RNTOT,RNGRUP,DFASSV,ZASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
        DIMENSION DFASSV(12)
C
C ZASS = CONTRIBUCION ASOCIATIVA AL FACTOR DE COMPRESIBILIDAD
C
        SUMZ = 0D0
C
        DO I=1,NGA
C
                SUMZ = SUMZ + RNGRUP(I)*DFASSV(I)
C
        END DO
C
C
        ZASS =-0.5d0 * SUMZ / RNTOT
C
        END
C///////////////////////////////////////////////////////////////////
        SUBROUTINE FCTRCOMPRES(ZFV,ZDISP,ZASS,ZSIS,IERRZSIS)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C ZSIS = FACTOR DE COMPRESIBILIDAD DEL SISTEMA
C
        ZSIS = ZFV + ZDISP + ZASS + 1.D0
C
        IF(ZSIS.LT.0) THEN
                                IERRZSIS = 1D0

                      ELSE
                                IERRZSIS = 0D0

        END IF
C
        END
C////////////////////////////////////////////////////////////////////
        SUBROUTINE FUGASIS(NC,FUGREP,FUGDISP,FUGASS,ZSIS,FUGSIS,
     @FFI,IERRFUG)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION FUGREP(20)
        DIMENSION FUGDISP(20)
        DIMENSION FUGASS(20)
        DIMENSION FUGSIS(20)
        DIMENSION FFI(NC)
C
C FUGSIS(I) = LOGARITMO NATURAL DEL COEFICIENTE DE FUGACIDAD DEL COMPONENTE
C                                                                   I
C FFI(I)=  COEFICIENTE DE FUGACIDAD DEL COMPONENTE I
C
        DO I = 1,NC
C
          FUGSIS(I)= FUGREP(I) + FUGDISP(I) + FUGASS(I) - DLOG(ZSIS)
C          evitando overflow
           IF(FUGSIS(I).GT.88) THEN
C
                                IERRFUG = 1
C
                              ELSE
C
                                IERRFUG = 0
C							write(*,*)'ln(phi',i,')=',FUGSIS(I)
                                FFI(I)= DEXP(   FUGSIS(I) )
           END IF
C
        END DO
C
        END
C//////////////////////////////////////////////////////////////////
C/////////////////////////////////////////////////
        SUBROUTINE ASSOCNMOL(NC,NGA,RN,ONUASOC,RNGRUP)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        DIMENSION RN(NC)
        DIMENSION ONUASOC(12,NC)
        DIMENSION RNGRUP(12)
C
C Calcula el numero de moles de cada grupo considerado en el termino de
C  asociacion. Se realiza el cálculo tanto para los grupos que presentan
C  sitios de asociacion como para los que no los presentan. Tal vez
C   alguna modificacion futura del modelo introduzca en la expresion
C    de la funcion de Helmholtz de asociacion, el efecto de dilucion
C     producido por los grupos que no asocian sobre los que si lo hacen
C
C RNGRUP(I)= NUMERO DE MOLES DEL GRUPO I (TERMINO DE ASOCIACION)
C

 	 DO I=1,NGA

                S1 = 0.D0
                DO M=1,NC
C
                           S1 = RN(M) * ONUASOC(I,M) + S1

C
                END DO
                RNGRUP(I) = S1
C
        END DO


   !   write(6,*)"nc=", nc
   !  	do i=1, nga
!		write(6,*)"rngrup=", rngrup(i)
!		do k=1,nc
!			write(6,*) "onuasoc(",k,",",i,")=",onuasoc(k,i)
!		end do
!	end do

C
        END
C//////////////////////////////////////////////////////////////////////////
        SUBROUTINE DERIFAASS(NGA,MASSOC,XNOASS,DFAASSXNOASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12)
        DIMENSION XNOASS(6,12)
        DIMENSION FAASS(6,12)
        DIMENSION DFAASSXNOASS(6,12)
C
C
        DO I=1,NGA
             IF(MASSOC(I).EQ.0) GO TO 13
			  DO K=1,MASSOC(I)
				If(XNOASS(K,I).gt.1d-13)then
					DFAASSXNOASS(K,I) = 1.d0 - XNOASS(K,I)
				END IF
                END DO
 13          CONTINUE
        END DO
C
        END
C/////////////////////////////////////////////////////////////////////////////
        SUBROUTINE DERIXNOASSVL(NGA,MASSOC,DELTASS,XNOASS,
     @DROSSV,DXNOASSV, ONUASOC, nc,RNGRUP,V)

C
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
C                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
       IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION MASSOC(12), sumsi(12)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
	  DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)
	  DIMENSION ONUASOC(12,NC)
C	dimension m_lambda(nga*2,nga*2), m_lambda1(nga*2,nga*2)
	dimension psiv(12)
      real*8  m_lambda(nga*2,nga*2), m_lambda1(nga*2,nga*2)
      integer order
C	real d1, psiv
	dimension indx(12)
c     	sumsi=0
c	do k=1, nc
c		do i=1, nga
c			If(ONUASOC(i,k).ne.0.)then
c				sumsi(i)=massoc(i)
c			end if
c		end do
c	end do


c	IF ((DELTASS(1,1,1,1)<1E-13).AND.(DELTASS(1,1,2,2)>1E-13)
c     *.AND.(SUMSI(1)==1)
c     *.AND.(DELTASS(1,2,2,2)>1E-13).AND.(SUMSI(2)==2))THEN

	!ASOCIACION CRUZADA ENTRE 1 GRUPO TIPO ESTER E UM GRUPO TIPO ALCOOL
	!COM UM SITIO (ESTER) E UM GRUPO COM 2 SITIOS(ALCOOL)
	!DELTASS(1,1,1,1)=0 (nao se autoassocia: ester ou alcano)
	!DELTASS(1,1,2,2).NE.0  OU  DELTASS(2,2,1,1).NE.0 SE ASOCIA DE FORMA CRUZADA
	!DELTASS(1,2,2,2).NE.0

c		IDAUT1=1
c		IDAUT2=2
c	   CALL DCROSSASOCVL(IDAUT1,IDAUT2,NC,DELTASS,XNOASS,DROSSV,
c     @DXNOASSV,RNGRUP,V)

c	ELSEIF ((DELTASS(1,2,1,2)<1E-13).AND.(DELTASS(1,1,1,2)>1E-13)
c     *.AND.(SUMSI(1)==1)
c    *.AND.(DELTASS(1,1,1,1)>1E-13).AND.(SUMSI(2)==1))THEN

	!ASOCIACION CRUZADA ENTRE 1 GRUPO TIPO ACIDO E UM GRUPO TIPO ESTER
	!COM UM SITIO (ACIDO) E UM GRUPO COM 1 SITIO(ESTER)
	!DELTASS(1,1,1,1)=0 (nao se autoassocia: ester ou alcano)
	!DELTASS(1,1,1,2).NE.0  OU  DELTASS(1,2,1,1).NE.0 SE ASOCIA DE FORMA CRUZADA
	!DELTASS(1,2,2,2).NE.0

c		IDAUT1=1
c		IDAUT2=2
c	   CALL DCROSSASOCVL3(IDAUT1,IDAUT2,NC,DELTASS,XNOASS,DROSSV,
c     @DXNOASSV,RNGRUP,V)



c	ELSEIF((DELTASS(1,1,1,1)>1E-13).AND.(DELTASS(1,1,2,2)>1E-13)
c     *.AND.(SUMSI(1)>1E-13)
c     *.AND.(DELTASS(1,2,2,2)>1E-13).AND.(SUMSI(2)==2))THEN


c		IDAUT1=1
c		IDAUT2=2
c	   CALL DCROSSASOCVL2(IDAUT1,IDAUT2,NC,
c     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V,nga)

c	ELSEIF((DELTASS(1,1,2,1)>1E-13).AND.(DELTASS(1,1,2,2)>1E-13)
c     *.AND.(SUMSI(1)==2)
c     *.AND.(DELTASS(1,2,2,2)>1E-13).AND.(SUMSI(2)==2))THEN


c		IDAUT1=1
c		IDAUT2=2
c	   CALL DCROSSASOCVL4(IDAUT1,IDAUT2,NC,
c     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)


c	ELSE !AUTO ASOCIACIONES EN PARALELO


c		do m=1, nga

c		  IF(SUMSI(m).EQ.2) THEN
c				  DO I=1,NGA
c					 IF(MASSOC(I).EQ.2) THEN
c						  IDAUT=I
c						  CALL DAUTOASOCVL(IDAUT,
c     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
c
c					  ELSE
c					   CONTINUE
c					 END IF
c				  END DO
c			elseif (sumsi(m).eq.1)then
c				  DO I=1,NGA
c					 IF(MASSOC(I).EQ.1) THEN
c						  IDAUT=I
c
c						  CALL DAUTOASOCVLumsitio(IDAUT,
c     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
c
c					  ELSE
c					   CONTINUE
c					 END IF
c				  END DO
c
c			ELSE
c			 CONTINUE
c		  END IF
c		end do
c
c	END IF


CCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  General

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
C                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
c Cálculo de la derivada de la fracción no asociada respecto VOLUMEN
cCalcula los elementos de la matriz deltapq

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
c Calcula los elementos de la matriz [Yp]
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

c colectando las derivadas en correspondiente subíndice
      w=0
	do j=1,2
	do i=1, nga
			w= w+1
	dxnoassv(j,i)=psiv(w)
	end do
	end do


	END
CCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOCVL3(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
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
	DXNOASSV(1,1)= (-t3*t5-t2*DELTASS(1,1,1,1)*t5-2D0*t9*t5*t4*
     #DELTASS(1,1,1,2)*t2-t4*t16*DELTASS(1,1,1,2)+t3*t4)/(1D0+2D0*t21*t4
     #+2D0*t9*t4+3D0*t26*DELTASS(1,1,1,1)*t5*DELTASS(1,1,1,2)+t15*
     #DELTASS(1,1,1,2)-t21)



	DXNOASSV(1,2)=-1D0*XNOASS(1,2)*DELTASS(1,1,1,2)*
     #(DROSSV(1)*XNOASS(1,1)+
     # (RNGRUP(1)/V)*DXNOASSV(1,1))
     #/(1D0+(RNGRUP(1)/V)*DELTASS(1,1,1,2)*XNOASS(1,1))

c     WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)
c	WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)

      END
C////////////////////////////////////////////////////////////////////////
        SUBROUTINE DCROSSASOCVL4(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR




C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL
C NRO. DE MOLES DEL COMPUESTO M

C

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
	DXNOASSV(2,2) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-
     #(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+
     #t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*
     #DELTASS(1,2,2,2))*
     #t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*
     #(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4
     #5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


	DXNOASSV(1,1)=(-DROSSV(1)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2D0
     #*XNOASS(1,1)*DROSSV(2)*DELTASS(1,1,2,2)*XNOASS(2,2)-2D0*
     #XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSV(2,2))
     #/(1+2D0*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2D0*
     #(RNGRUP(2)/V)*DELTASS(1,1,2,2)*XNOASS(2,2))



	DXNOASSV(1,2) = DXNOASSV(2,2)
c      WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)

c		WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)
c		WRITE(2,*)"DXNOASSV(2,2)=",DXNOASSV(2,2)


      END
C************************************************************************************
C////////////////////////////////////////////////////////////////////////
       SUBROUTINE DCROSSASOCVL2(IDAUT1,IDAUT2,NC,
     @DELTASS,XNOASS,DROSSV,DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL NRO. DE MOLES DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION RNGRUP(12)

C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C DROSSN(I,J,K,L,M) = DERIVADA DE LA FUERZA DE ASOCIACION ENTRE
C                  EL SITIO I DEL GRUPO J, Y EL SITIO K
C                  DEL GRUPO L con respecto al numero de moles
C                  del componente M
C
C DXNOASSN(I,J,M) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL NRO. DE
C                    MOLES DEL COMPUESTO M (ESTO ES LO QUE SE CALCULA)
C  CALCULO DE BLOQUE AUXILIAR




C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL
C NRO. DE MOLES DEL COMPUESTO M
c      xnoass(1,1)=0.5873
c     xnoass(1,2)=0.9437
c	xnoass(2,2)=0.9437
c      rngrup(1)=1
c	rngrup(2)=0
c	deltass(1,1,1,1)=81256.6589
c	deltass(1,2,1,1)=6899.99
c	deltass(2,2,1,1)=6899.99
c	deltass(1,1,1,2)=6899.99
c	deltass(2,2,1,2)=585.15
c	deltass(1,1,2,2)=6899.993
c	deltass(1,2,2,2)=585.1586
c	deltass(2,2,2,2)=0
c	v=67937.4025

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
	DXNOASSV(2,2) = (-(2D0*t2*t3*t4-4D0*t9*DELTASS(1,2,2,2)*t4)*t16-
     #(-DELTASS(1,1,2,2)*t4*DELTASS(1,2,2,2)-t21+t24)*t26-(-t24-t29*t8+
     #t21)*t15)/(4D0*(DELTASS(1,1,1,1)*t34*t3-2D0*t34*t8*
     #DELTASS(1,2,2,2))*
     #t26+3D0*(-DELTASS(1,1,2,2)*t1*DELTASS(1,2,2,2)-t45+t47)*t15+2D0*
     #(DELTASS(1,1,1,1)-DELTASS(1,1,2,2)-t47-t28*t8+t4
     #5)*t14-2D0*DELTASS(1,1,1,1)+DELTASS(1,1,2,2))


	DXNOASSV(1,1)=(-DROSSV(1)*DELTASS(1,1,1,1)*XNOASS(1,1)**2-2D0
     #*XNOASS(1,1)*DROSSV(2)*DELTASS(1,1,2,2)*XNOASS(2,2)-2D0*
     #XNOASS(1,1)*(RNGRUP(2)/V)*DELTASS(1,1,2,2)*DXNOASSV(2,2))
     #/(1+2D0*(RNGRUP(1)/V)*DELTASS(1,1,1,1)*XNOASS(1,1)+2D0*
     #(RNGRUP(2)/V)*DELTASS(1,1,2,2)*XNOASS(2,2))



	DXNOASSV(1,2) = DXNOASSV(2,2)
c      WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)

c		WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)
c		WRITE(2,*)"DXNOASSV(2,2)=",DXNOASSV(2,2)


      END




C************************************************************************************
	SUBROUTINE DCROSSASOCVL(IDAUT1,IDAUT2,NC,DELTASS,XNOASS,DROSSV,
     @DXNOASSV,RNGRUP,V)

        IMPLICIT REAL*8 (A-H,O-Z)
C
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

	DXNOASSV(2,2) = (-t6-t7*t5-2d0*t9*t5*t4*DELTASS(1,2,2,2)*t2-t4*t16
     #*DELTASS(1,1,2,2)-t16*DELTASS(1,1,2,2)*t20*t5-t23*t6+t7*t4)/(1d0+
     #2d0*t27+2d0*t9*t4+3d0*t31*DELTASS(1,1,2,2)*t5*DELTASS(1,2,2,2)+t23
     #+2d0*t23*t27-t9)

	DXNOASSV(1,1) =-1D0*XNOASS(1,1)*DELTASS(1,1,2,2)*(RNGRUP(2)/V*
     #DXNOASSV(2,2)+XNOASS(2,2)*DROSSV(2))/(1D0+RNGRUP(2)/V*
     #XNOASS(2,2)*DELTASS(1,1,2,2))

	DXNOASSV(1,2) = -1D0*XNOASS(1,2)*DELTASS(1,2,2,2)*(RNGRUP(2)/V*
     #DXNOASSV(2,2)+XNOASS(2,2)*DROSSV(2))/(1D0+RNGRUP(2)/V*
     #XNOASS(2,2)*DELTASS(1,2,2,2))

c	WRITE(2,*)"DXNOASSV(1,1)=",DXNOASSV(1,1)
c      WRITE(2,*)"DXNOASSV(1,2)=",DXNOASSV(1,2)
c      WRITE(2,*)"DXNOASSV(2,2)=",DXNOASSV(2,2)


	END
C//////////////////////////////////////////////////////////////////////////

C//////////////////////////////////////////////////////////////////////////


        SUBROUTINE DAUTOASOCVL(IDAUT,DELTASS,XNOASS,DROSSV,DXNOASSV,
     @  RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
	  DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)
C
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL VOLUMEN DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
C
C
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
C   QUE AUTOASOCIA (CON DOS SITIOS).
C
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C
C
C DDLTSSV(I,J,K,L) = DERIVADA DEL DELTA DE ASOCIACION ENTRE EL SITIO I
C                       DEL grupo J Y EL SITIO K DEL grupo L CON RESPECTO
C                               AL VOLUMEN
C
C DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
C                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
C
C  CALCULO DE BLOQUE AUXILIAR

        IF(DELTASS(1,IDAUT,2,IDAUT).NE.0) THEN

C
        BLOQ =0.D0
        HELPI =0.D0
C
C
        BLOQ =1D0+2D0*DELTASS(1,IDAUT,2,IDAUT)*XNOASS(1,IDAUT)*
     @RNGRUP(IDAUT)/V
C
        BLOQ = 1D0 / BLOQ
C
C EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
C                               respecto al DELTASS(1,IDAUT,2,IDAUT)
C
C
        BLOQ=(-1D0)*DELTASS(1,IDAUT,2,IDAUT)*(XNOASS(1,IDAUT))**2D0*BLOQ
C
C CALCULO DE LA DERIVADA RESPECTO VOLUMEN
C
C
        DXNOASSV(1,IDAUT) = BLOQ  * DROSSV(IDAUT)
C
C
C IGUALACION CON LA DERIV. DE LA FRACC NO ASS DEL SITIO 2 DEL GRUPO
C  IDAUT CON RESPECTO AL VOLUMEN
C
C
        DXNOASSV(2,IDAUT)   = DXNOASSV(1,IDAUT)

c		WRITE(2,*)"DXNOASSV(1,",IDAUT,")=",DXNOASSV(1,IDAUT)
c		WRITE(2,*)"DXNOASSV(2,",IDAUT,")=",DXNOASSV(2,IDAUT)
	END IF
C
C
C
        END
C//////////////////////////////////////////////////////////////////////////
        SUBROUTINE DAUTOASOCVLumsitio(IDAUT,DELTASS,XNOASS,DROSSV,
     @DXNOASSV,RNGRUP,V)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION DROSSV(12)
	  DIMENSION RNGRUP(12)
        DIMENSION DXNOASSV(6,12)
C
C
C NO GENERAL
C
C CALCULA DERIV. RESP. AL VOLUMEN DE LA FRACCION NO ASOCIADA
C  SOLO EN EL CASO DE QUE ASOCIE UN UNICO GRUPO (HACIENDOLO CON 2
C SITIOS)
C
C
C XNOASS(I,J) = FRACCION NO ASOCIADA DEL SITIO I DEL GRUPO J
C
C
C IDAUT ES IGUAL AL NUMERO DE IDENTIFICACION ASIGNADO AL GRUPO
C   QUE AUTOASOCIA (CON DOS SITIOS).
C
C XNOASS(1,IDAUT)= FRACCION NO ASOC. CORRESP. AL SITIO 1 DEL GRUPO
C                                                       IDAUT
C
C
C DDLTSSV(I,J,K,L) = DERIVADA DEL DELTA DE ASOCIACION ENTRE EL SITIO I
C                       DEL grupo J Y EL SITIO K DEL grupo L CON RESPECTO
C                               AL VOLUMEN
C
C DXNOASSV(I,J) = DERIVADA DE LA XNOASS(I,J) CON RESPECTO AL VOLUMEN
C                   DEL SISTEMA (ESTO ES LO QUE SE CALCULA)
C
C  CALCULO DE BLOQUE AUXILIAR
C
        IF(DELTASS(1,IDAUT,1,IDAUT).NE.0) THEN

        BLOQ =0.D0
        HELPI =0.D0
C
C
        BLOQ =1D0+2D0*DELTASS(1,IDAUT,1,IDAUT)*XNOASS(1,IDAUT)*
     @RNGRUP(IDAUT)/V
C
        BLOQ = 1D0 / BLOQ
C
C EN LA LINEA DE ABAJO    BLOQ  (IZQ) es la derivada de XNOASS(1,IDAUT) con
C                               respecto al DELTASS(1,IDAUT,2,IDAUT)
C
C
        BLOQ=(-1D0)*DELTASS(1,IDAUT,1,IDAUT)*(XNOASS(1,IDAUT))**2D0*BLOQ
C
C CALCULO DE LA DERIVADA RESPECTO VOLUMEN
C
C
        DXNOASSV(1,IDAUT) = BLOQ  * DROSSV(IDAUT)
c 	WRITE(2,*)"DXNOASSV(1,",IDAUT,")=",DXNOASSV(1,IDAUT)
	END IF
C
        END
C//////////////////////////////////////////////////////////////////////////

        SUBROUTINE
     @ACLAR(NCOD,NC,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)
C
C PONE EN CERO LAS VARIABLES DE SALIDA CUANDO COMIENZA LA CORRIDA
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION FFI(NC)  !DIMENSIONAMIENTO COEFICIENTES DE FUGACIDAD
C
        ARESTOT=0.D0
C
C SI NCOD=2, ZSIS entra como dato, por lo que no se blanquea esa variable
C
        IF(NCOD.NE.2) ZSIS =0.D0
C
        DO I=1,NC
                FFI(I)=0.D0
        END DO
C
        IERRFUG=0
        IERRZSIS=0
        IERRREP=0
	iErrAssoc=0
C
C
        END !FIN ACLAR
C//////////////////////////////////////////////////////////////////
        SUBROUTINE CLEARANCE(NC,NG,NGA,ARESASS,ARESDISP,ARESFV,QTOT,
     @RLAMBDA1,RLAMBDA2,RLAMBDA3,RNTOT,Y,ZASS,ZDISP,ZFV,
     @MASSOC,D,DL1N,DL2N,DL3N,DYN,FUGASS,FUGDISP,FUGREP,FUGSIS,
     @RMS,H3,H7,PS,DFASSN,DXNOASSN,H2,H4,H5,H6,TITA,TSPLINE,
     @DLTG,TAU,VERDKIJ,GGG,DFASSV,FASS,RNGRUP,DXNOASSV,FAASS,XNOASS,
     @DFAASSXNOASS,DROSSV,DELTASS,DROSSN)
C
C PONE EN CERO UN NUMERO IMPORTANTE DE VARIABLES
C
        IMPLICIT REAL*8 (A-H,O-Z)
C
        DIMENSION D(20)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
C DEL TERMINO DE DISPERSION
C
        DIMENSION TITA(10), TAU(10,10)
        DIMENSION TSPLINE(10)
        DIMENSION GGG(10,10)
        DIMENSION VERDKIJ(10,10)
        DIMENSION DLTG(10,10)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
C
        DIMENSION RNGRUP(12)  !nro. de moles de los grupos
C
C MASSOC(J) = NRO. DE SITIOS DE ASOCIACION ASIGNADOS AL GRUPO J
C
C DELTASS(I,J,K,L) = FUERZA ADIMENSIONAL DE ASOCIACION ENTRE  el sitio
C                 I del GRUPO J, y el sitio K del grupo L
C
C
        DIMENSION MASSOC(12)
        DIMENSION DELTASS(6,12,6,12)
        DIMENSION XNOASS(6,12)
        DIMENSION FASS(12)
        DIMENSION FAASS(6,12)
C
C DIMENSIONAMIENTO RELAC CON COEF FUGACIDAD REPULSIVO
C
        DIMENSION  DL1N(20),DL2N(20),DL3N(20),DYN(20)
        DIMENSION FUGREP(20)
C
C DIMENSIONAMIENTO DE VARIABLES PARA EL CALCULO
C    DEL COEFICIENTE DE FUGACIDAD DE DISPERSIVO
C
        DIMENSION PS(20,10), RMS(20), H2(10)
        DIMENSION H3(20,10), H4(10), H5(10)
        DIMENSION H6(10), H7(20,10)
        DIMENSION FUGDISP(20)
C
C DIMENSIONAMIENTO RELAC. CON COEF. FUG. ASOCIATIVO
C
        DIMENSION DROSSN(12,20)
        DIMENSION DXNOASSN(6,12,20)
        DIMENSION DFAASSXNOASS(6,12)
        DIMENSION DFASSN(12,20)
        DIMENSION FUGASS(20)
C
C DIMENSIONAMIENTO RELAC. CON FACTOR DE COMPRES. ASOCIATIVO
C
        DIMENSION DROSSV(6,12,6,12)
        DIMENSION DXNOASSV(6,12)
        DIMENSION DFASSV(12)
C
C
        DIMENSION FUGSIS(20)
C
C
C LIMPIEZA DE VARIABLES
C
C
C         VARIABLES SIMPLES
C
C         ARESASS                     =0.
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

C       Arreglos
C
C LA PROXIMA ASIGNACION PERMITE ASEGURAR EL VACIADO DE LOS ARREGLOS
C  CUANDO LA HISTORIA DEL CALCULO FUE TAL , QUE PRIMERO SE REALIZO UN
C    CALCULO' CON CIERTO NUMERO DE COMPONENTES, Y LUEGO CON UN NUMERO
C      MENOR EN UNA UNIDAD ( P.E. primero un flash binario, y luego una
C          presion de vapor de un componente puro)
C
        NCSAVE = NC + 1
C
C
C
        DO I=1,NCSAVE
C
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
C
           DO J=1,NG
                H3(I,J)=0.D0
                H7(I,J)=0.D0
                PS(I,J)=0.D0
           END DO
C
            DO K=1,NGA
C
                        DFASSN(K,I)=0.D0
C
                        DO L=1,6
C
                                DXNOASSN(L,K,I)=0.D0
C
                        END DO
C
            END DO
C
        END DO
C       --------------------
C
C
        DO I=1,NG
C
          H2(I)=0.D0
          H4(I)=0.D0
          H5(I)=0.D0
          H6(I)=0.D0
          TITA(I)=0.D0
          TSPLINE(I)=0.D0
C
            DO J=1,NG
                DLTG(I,J)=0.D0
                TAU(I,J)=0.D0
                VERDKIJ(I,J)=0.D0
                GGG(I,J)=0.D0
            END DO
C
        END DO
C
C
        DO I=1,NGA
C
         DFASSV(I)=0.D0
         FASS(I)=0.D0
         RNGRUP(I)=0.D0
C
                DO J=1,6
C
                        DXNOASSV(J,I)=0.D0
                        FAASS(J,I)=0.D0
                        XNOASS(J,I)=0.D0
                        DFAASSXNOASS(J,I)=0.D0
C
                END DO
C
        END DO
C
        DO I=1,NGA
                DO J=1,6
                        DO K=1,NGA
                                DO L=1,6
C
                                   DROSSV(J,I,L,K)=0.D0
                                   DELTASS(J,I,L,K)=0.D0
C
                                   DO M=1,NCSAVE
C
                                        DROSSN(J,M)=0.D0
C
                                   END DO
                                END DO
                        END DO
                END DO
        END DO
C
        END