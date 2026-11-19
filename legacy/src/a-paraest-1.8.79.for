	PROGRAM ASSOCPARAEST

CCCCCCCCCCCCCCCCCCCCCCCCC !!ALFONSINA!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Modificado en 1 de Diciembre del 2005 en el cálculo de la fracción no asociada
c y en el cálculo de las derivadas de las fracciones no asociadas respecto al volumen
c y al número de moles del componente, para que todos los casos particulares se puedan
c resolver con un único procedimiento basado en los paper de: Michelsen, et al.
c (Fluid Phase Equilibria, 180(2001)165-174 )y del paper de Tan, et al.
c (Ind. Eng. Chem. Res, 2004,43,203-208).

c importante!!! Tener cuidado con la entrada de datos de los enass y rkass.
c Como por ejemplo ahora la entrada de datos del ácido con el agua debe ser:
c 1  2
c 6300.  ! 1 1 1 1  AUTO ASOC. COOH  ENERGIA DE ASOCIACION,,
c 4500.  ! 1 1 1 2  CRUZ.
c 4500.  ! 1 1 2 2  CRUZ.
c 4500.  ! 1 2 1 1  CRUZ.
c 4500.  ! 2 2 1 1  CRUZ.
c 0.     ! 1 2 1 2  AUTO ASOC.  OH
c 2700.  ! 1 2 2 2
c 2700.  ! 2 2 1 2
c 0.     ! 2 2 2 2
c 0.02   ! 1 1 1 1  AUTO ASOC. COOH  VOLUMEN DE ASOCIACION
c 0.1313 ! 1 1 1 2  CRUZ.
c 0.1313 ! 1 1 2 2  CRUZ.
c 0.1313 ! 1 2 1 1   CRUZ.
c 0.1313 ! 2 2 1 1   CRUZ.
c 0.     ! 1 2 1 2  AUTO ASOC.  OH
c .8621  ! 1 2 2 2
c .8621  ! 2 2 1 2
c 0.     ! 2 2 2 2
c
c y no como era antes:
c 1  2
c 6300.  ! 1 1 1 1  AUTO ASOC. COOH  ENERGIA DE ASOCIACION
c 0.     ! 1 1 1 2  CRUZ.
c 4500.  ! 1 1 2 2  CRUZ.
c 0.     ! 1 2 1 1  CRUZ.
c 4500.  ! 2 2 1 1  CRUZ.
c 0.     ! 1 2 1 2  AUTO ASOC.  OH
c 2700.  ! 1 2 2 2
c 2700.  ! 2 2 1 2
c 0.     ! 2 2 2 2
c 0.02   ! 1 1 1 1  AUTO ASOC. COOH  VOLUMEN DE ASOCIACION
c 0.0    ! 1 1 1 2  CRUZ.
c 0.1313 ! 1 1 2 2  CRUZ.
c 0.0    ! 1 2 1 1   CRUZ.
c 0.1313 ! 2 2 1 1   CRUZ.
c 0.     ! 1 2 1 2  AUTO ASOC.  OH
c .8621  ! 1 2 2 2
c .8621  ! 2 2 1 2
c 0.     ! 2 2 2 2
c Esto se debe a que esta entrada de datos lleva a las ecuaciones cuárticas y cúbicas
c que se usaban para los casos particulares.

CCCCCCCCCCCCCCCCCCCCCCCCC !!ALFONSINA!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                          !!ALFONSINA!!


C
C AJUSTE DE PARAMETROS CON MODELO ASSOC-GC-EOS: La estructura molecular
C  puede ser definida de distinta manera según se trate de cuantificar
C  la contribucion dispersiva o de asociacion
C
C
C ADMITE A DC COMO UN PARAMETRO AJUSTABLE FRENTE A DATOS DE EQUILIBRIO
C
C  LAS OPCIONES PARA DC SON LAS SIGUIENTES
C
C    - ESTIMARLO EN BASE AL PUNTO CRITICO (ECUACION EXPLICITA, NO HAY
C					   AJUSTE DE DC)
C
C    - AJUSTARLO FRENTE A UN DATO DETERMINADO DE LA CURVA DE PRESION
C       DE VAPOR DEL COMPONENTE EN CUESTION
C
C    - CONSIDERARLO COMO CUALQUIER OTRO PARAMETRO AJUSTABLE. ESTO SIGNIFICA
C       QUE PODRA SER AJUSTADO FRENTE A TODOS LOS DATOS QUE SE ALIMENTEN
C         EN QUE TENGA PARTICIPACION EL COMPONENTE EN CUESTION, SEAN
C           ESTOS BINARIOS Y/O DE COMPONENTE PURO.
C
C CON PENALIZACION
C
C
C  LA LECTURA DE DATOS PUEDE HACERSE POR PANTALLA O DESDE UN ARCHIVO DE
C   DATOS
C
C  CONTIENE MODULO DE ESCALADO DE LOS PARAMETROS A AJUSTAR
C
c	IMPLICIT REAL*8 (A-H,O-Z)
	Implicit none
C
C LA DECLARACION QUE VIENE LA EXIGE LA RUTINA ZXSSQ DE LA IMSL
C
	EXTERNAL FUNC
C
C DIM DISPERSION
C
	INTEGER NU(10,20)
	DOUBLE PRECISION QAREA(10), TSTR(10), GSTR(10)
	DOUBLE PRECISION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DOUBLE PRECISION ALFAIJ(10,10)
	DOUBLE PRECISION TSPLINE(10)
C
C DIM ASOCIACION
C
	INTEGER MASSOC(12)
      DOUBLE PRECISION ENASS(6,12,6,12), RKASS(6,12,6,12)
      DOUBLE PRECISION ONUASOC(12,20)
C
C DOBLE PRECISIONAMIENTO DIAMETRO DE ESFERAS RIGIDAS EN EL PTO CRITICO
C
	DOUBLE PRECISION DC(20)
C
C DIM AJUSTE DE PARAM
C
	DOUBLE PRECISION PARM(4)
C
	INTEGER IDSAT(20)
C
	DOUBLE PRECISION TEMPCRIT(20)
	DOUBLE PRECISION PRCRIT(20), TSAT(20), PSAT(20)
	DOUBLE PRECISION TEXP(260), PEXP(260), XEXP(260), YEXP(260), xexp2(260)
	DOUBLE PRECISION rgamexp1(260),rgamexp2(260),rgamexp3(260)
C
C DOUBLE PRECISIONAMIENTO  NRO. DE DATOS DE PVAP COMP PUROS
C
	INTEGER NSAT(20)
	INTEGER IOBJSAT(20) !Funcion objetivo para el comp. puro
C
C DOUBLE PRECISIONAMIENTOS QUE TIENEN QUE VER CON LOS DATASETS
C
	INTEGER IDC(25,3), NOBS(25), IHEN(25), mez(25)
	INTEGER NTEXT(25,65)
	INTEGER IFL(25)
C
C ID ==> identificacion o identicos
C IDID ==> identificacion de identicos
C
C Los parametros se identifican a traves de su naturaleza y de sus
C	indices
C
	INTEGER ITEMP(8)
	INTEGER IDG(8,2)		!identificacion (indices) de param gii o kij a
C							!ajust.
	INTEGER NIDG(8)			!nro. de param. identicos a c/gii o kij
	INTEGER IDIDG(104,2)	!identificacion (indices) de IDENTicos a gii ó
C							!kij
C
C                8 X 13 = 104
C
	INTEGER IDQ(8)				!IDENTif param de area
C
	INTEGER IDALF(8,2)			!identif (indices) de parametros ALFAij
	INTEGER NIDALF(8)			!nro. de param. IDENTicos a c/alfaij
	INTEGER IDIDALF(104,2)		!identificacion (indices) de IDENTicos
C
	INTEGER IDENASS(8,4)		!IDENTif (indices) de parametros ENASSijkl
	INTEGER NIDENASS(8)			!Nro. de param. IDENTicos  a c/ENASSijkl
	INTEGER IDIDENASS(104,4)	!IDENTif (indices) de IDENTicos a ENASSijkl
C
	INTEGER IDRKASS(8,4)		!Identif (indices) de param. RKASSijkl
	INTEGER NIDRKASS(8)			!Nro. de param. identicos a c/RKASSijkl
	INTEGER IDIDRKASS(104,4)	!Identif (indices) de identicos a c/RKASSijkl
C
	INTEGER IDDC(8)				!Índices para dc
C
      INTEGER idtc(8)				!indice de los parametros temperatura critica
c
	DOUBLE PRECISION XX(8)
	INTEGER IDXX(8)
C
C	IDXX(I) =	identificación del tipo de parámetro a ser ajustado.
C				IDXX se utiliza para saber que factor de escalado
C				debe aplicarse al parámetro que se esta ajustando.
C
C	IDXX(i) =	1   g*
C				2	k*ij
C				3	g'ij / g''ij / k'ij   coeficientes var. Con temp
C 				4	q (param. de área)
C 				5	alfa  (param. de no azar)
C 				6	param. de asoc. Energetico
C 				7	param. de asoc. Volumetrico
c				8	diámetro critico
c				9	temperatura crítica
C
	DOUBLE PRECISION XJAC(260,8), XJTJ(36), WORK(436)
	DOUBLE PRECISION FCION(260)
	DOUBLE PRECISION FCAUX(260)
	DOUBLE PRECISION FCAUX2(260)
	DOUBLE PRECISION FCAUX3(260)
C
	DOUBLE PRECISION PSTCLC(20)		!PRESION DE VAPOR(ATM) CALCULADA A LA
C									!TEMPERATURA TSAT(I)
C
C DOBLE PRECISIONAMIENTO DE VARIABLES CALCULADAS CON EL MODELO
C
	DOUBLE PRECISION PCALC(260), XCALC(260), YCALC(260)
	DOUBLE PRECISION ZLCALC(260), ZVCALC(260)
	DOUBLE PRECISION rgamcalc1(260),rgamcalc2(260),rgamcalc3(260)
	DOUBLE PRECISION zlicalc(260)

C
C DOBLE PRECISIONAMIENTO DESVIACIONES STANDAR PARA MEZCLAS
C
	DOUBLE PRECISION STDVA(25)
	DOUBLE PRECISION STDVB(25)
	DOUBLE PRECISION STDVC(25)
	INTEGER NC,NG,NGA,nga2,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA
     @,NENASS,NRKASS,NPAR,NCXP,IOPT,NUMDC,NUMTC,OPDER,NOPCION,MAXFN,
     @IXJAC,INFER,IER
	double precision W1,W2,W3,W4,W5,W6,W7,W8,EPS,DELTAZ,SSQ

C

C COMMON PARA SUBRUTINA FUNC REQUERIDA POR IMSL
C
C

	COMMON/FORFUNC/NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,MASSOC,ENASS,RKASS,W1,W2,W3,W4,W5,
     @w6,W7,W8,ISAT,IGLE,IVE,IUNI,NNP,NNQ,NALFA,NENASS,NRKASS,
     @  NSAT,IDC,NOBS,IHEN,mez,IFL,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @	IDRKASS,NIDRKASS,IDIDRKASS,DC,PSTCLC,XCALC,YCALC,ZLCALC,
     @  ZVCALC,FCAUX,PCALC,IOBJSAT,FCAUX2,fcaux3,STDVA,STDVB,stdvc,
     @numtc,idtc,NUMDC,IDDC,NGA
	COMMON/FORFUNC2/rgamcalc1, rgamcalc2,rgamcalc3,
     @ rgamexp1, rgamexp2, rgamexp3
	COMMON/FORFUNC3/ONUASOC,zlicalc,PRCRIT,TEMPCRIT,TSAT,PSAT
C
C
	COMMON/IDSC/IDXX
C
C LIMITES
C
C 10 GRUPOS DIFERENTES EN LA CONTRIBUCION DISPERSIVA
C 12 GRUPOS DIFERENTES EN LA CONTRIBUCION DE ASOCIACION
C 20 COMPONENTES
C 180 DATOS EXPERIMENTALES DE CUALQUIER TIPO, SEA  DE PRESION DE
C     VAPOR DE COMP PUROS O CUALQUIER OTRO TIPO DE DATO DE EQUILIBRIO
C 25 DATA SETS DE ELV
C  8  (Ajuste de a lo sumo 8 parametros diferentes)
C 13  (A lo sumo 13 parametros IDENTicos a alguno de los ajustados)
C
C        IDC(J,1)  = CONTIENE EL NUMERO QUE IDENTIFICA AL COMPONENTE
C	                                                                                                    	      PARA EL QUE SE DAN LAS FRACCIONES MOLARES EXPERI
C		       MENTALES, EN EL DATASET NUMERO J
C
C	 IDC(J,2)  =  CONTIENE EL NUMERO QUE IDENTIFICA AL COMPONENTE
C		      PARA EL QUE no SE DAN LAS FRACCIONES MOLARES EXPERI
C		       MENTALES, EN EL DATASET NUMERO J
C


C		Lectura desde el archivo de entrada ASSOCPARA.DAT
C
		CALL READ(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,W1,W2,W3,W4,W5,
     @  w6,W7,W8,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,
     @NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,
     @IFL,NTEXT,rgamexp1, rgamexp2, rgamexp3,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,XX,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @	IDRKASS,NIDRKASS,IDIDRKASS,NCXP,IOPT,IOBJSAT,NUMDC,
     @  IDDC,numtc, idtc)
C
C  . .  . . . . . . . . . . . . . . . . . . . . .. .  . . ..  ..  .

C Inicializacion de DC en base al punto critico, para todo componente
C
	CALL INICDC(NC,PRCRIT,TEMPCRIT,DC)
C
C CUANDO NUMDC NO ES CERO: CARGA DE LOS VALORES DE DC INICIALES EN XX
C
	CALL FORDC(NUMDC,IDDC,DC,NPAR,XX)
C
C   Calculo de las temperaturas de spline para los valores iniciales de
C    los parametros dispersivos
C
      CALL SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
C
C   Escribe  parte de los datos en archivo de salida paraest.out
C
	CALL WRITESAL1(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,W1,W2,W3,W4,W5,
     @ w6,W7,W8, NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,
     @NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,
     @IFL,NTEXT,TEXP,PEXP,XEXP,xexp2,YEXP,XX,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,IDRKASS,
     @NIDRKASS,IDIDRKASS,DC,IOPT,NUMDC,IDDC,TSPLINE,numtc,idtc)
C
*	CALL LIB$INIT_TIMER
C
C	IDENTIFICACION DE LOS VALORES CONTENIDOS EN XX, PARA LUEGO
C  PODER ASIGNARLES FACTORES DE ESCALADO ADECUADOS
C

	CALL FORSCAL(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,
     @numtc,NUMDC)
C
C ESCALADO DE LOS VALORES DE XX
C
	CALL XXSCAL(NPAR,XX,IDXX)
C
C
C       Calculos de simulacion o Ajuste
C
C       ----------------------------------------------------
	IF(NPRED.EQ.1)	THEN   !Prediccion (o simulacion)
C
			 IF (OPDER.EQ.'S')  NOPCION = 1 !chequea deriv
			 IF (OPDER.EQ.'N')  NOPCION = 0 ! no chequea deriv
C
		         CALL FUNC(XX,NCXP,NPAR,FCION)
C
C       ----------------------------------------------------
			ELSE   !Ajuste de Parametros
C       ----------------------------------------------------
C
C
C 			 OPTIMIZACION
C
C
C   			 Parametros para la rutina ZXSSQ
C
		 		EPS = 1.D-6
				DELTAZ = 1.D-6
				MAXFN = 100
			      	PARM(3)=2.D2
			      	PARM(4)=.1D0
C		comentario      IXJAC = 180
				IXJAC = NCXP
C
			 NOPCION = 0   !No chequeo de deriv.
C
C
			 CALL ZXSSQ(FUNC,NCXP,NPAR,NSIG,EPS,
     @DELTAZ,MAXFN,IOPT,PARM,XX,SSQ,FCION,XJAC,IXJAC,
     @XJTJ,WORK,INFER,IER)
C
C 			 ESCRIBE ESTADISTICAS DE LA OPTIMIZACION
C
			 CALL WRITESTAT(SSQ,XJAC,WORK,INFER,IER)
C
C   			 Calculo de las temperaturas de spline para
C			 los valores FINALES de los parametros dispersivos
C
C
			 CALL SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
C
C
			 IF (OPDER.EQ.'S')  NOPCION = 1
C
C			 Simulacion con param. finales y chequeo de derivadas
C
			 IF (OPDER.EQ.'S')  CALL FUNC(XX,NCXP,NPAR,FCION)
C
C
	END IF
C       ----------------------------------------------------
C
C ESCRIBE PARAMETROS FINALES,Y VALORES EXPERIMENTALES Y CALCULADOS
C
	CALL WRITESAL2(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,ISAT,
     @  IGLE,IVE,IUNI,
     @  TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,
     @  NTEXT,
     @TEXP,PEXP,XEXP,xexp2, YEXP,IDSAT,DC,tempcrit,XCALC,YCALC,ZLCALC,
     @  ZVCALC,FCAUX,PCALC,IOBJSAT,FCAUX2,fcaux3,STDVA,STDVB,stdvc,
     @ rgamcalc1, rgamcalc2,rgamcalc3,numtc,idtc,NUMDC,IDDC,TSPLINE,
     @rgamexp1, rgamexp2,rgamexp3,zlicalc)

C
C Informa si hubo chequeo de derivadas
C
	END !FIN ASSOCPARAEST
C     /////////////////////////////////////////////////////////////////////
	SUBROUTINE READ(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,W1,W2,W3,W4,W5,
     @  w6,W7,W8,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,
     @NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,
     @IFL,NTEXT,rgamexp1, rgamexp2, rgamexp3,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,XX,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @	IDRKASS,NIDRKASS,IDIDRKASS,NCXP,IOPT,IOBJSAT,NUMDC,
     @  IDDC,numtc, idtc)
C
C LEE TODO DE UN ARCHIVO DE SALIDA
C
	IMPLICIT REAL*8 (A-H,O-Z)
C
C  DIM DISPERSION
C
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
C
C  DIM ASOCIACION
C
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,20)
C
C  DIM AJUSTE DE PARAM
C
	DIMENSION PARM(4)
C
	DIMENSION IDSAT(20)
C
	DIMENSION  TEMPCRIT(20)
	DIMENSION  PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260),xexp2(260),YEXP(260)
	dimension rgamexp1(260),rgamexp2(260),rgamexp3(260)
***************************************************************************
*	double precision  DC(20)
***************************************************************************
C
C DIMENSIONAMIENTO  NRO. DE DATOS DE PVAP COMP PUROS
C
	DIMENSION NSAT(20)
	DIMENSION IOBJSAT(20)
C
C DIMENSIONAMIENTOS QUE TIENEN QUE VER CON LOS DATASETS
C
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
C
C ID ==> identificacion o identicos
C IDID ==> identificacion de identicos
C
C   los parametros se identifican a traves de su naturaleza y de sus
C	inidices
C
	DIMENSION ITEMP(8)
	DIMENSION IDG(8,2) ! identificacion (indices) de param gii o kij a
C									ajust.
	DIMENSION NIDG(8)   !nro. de param. IDENTicos a c/gii o kij
	DIMENSION IDIDG(104,2) !identificacion (indices) de IDENTicos a gii o
C									  kij
C
C                8 X 13 = 104
C

	DIMENSION IDQ(8)  !IDENTif param de area
C
	DIMENSION IDALF(8,2) !identif (indices) de parametros ALFAij
	DIMENSION NIDALF(8) !nro. de param. IDENTicos a c/alfaij
	DIMENSION IDIDALF(104,2) !identificacion (indices) de IDENTicos
C
	DIMENSION IDENASS(8,4) !IDENTif (indices) de parametros ENASSijkl
	DIMENSION NIDENASS(8) !Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION IDIDENASS(104,4) !IDENTif (indices) de IDENTicos a ENASSijkl
C
	DIMENSION IDRKASS(8,4) !Identif (indices) de param. RKASSijkl
	DIMENSION NIDRKASS(8) !Nro. de param. identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4) !Identif (indices) de identicos a c/RKASSijkl
C
	DIMENSION IDDC(8)
c
	dimension idtc(8)
C
	DIMENSION XX(8)
C
	OPEN(UNIT=1,TYPE='OLD',NAME='paraest.dat',FORM='FORMATTED',
     @ ACCESS='SEQUENTIAL')
c		OPEN(UNIT=4,NAME='XNOASS.TXT',FORM='FORMATTED')


C
C
C
c	READ(1,222) W1,W2,W3,W4,W5,PARM(1),PARM(2),IOPT
c 222	FORMAT(7F10.4,I3)
C	write(*,*)'ssdfaryasd a h st hs th a et a r gry n'
	READ(1,*) W1,W2,W3,W4,W5,W6,W7,W8,PARM(1),PARM(2),IOPT
	write(*,*)W1,W2,W3,W4,W5,W6,W7,W8,PARM(1),PARM(2),IOPT
C
 1	FORMAT(8F10.4)
C
	READ(1,*) NC,NG,NGA,NSIG
C
 2	FORMAT(20I3)
C
	READ(1,*) ISAT,IGLE,IVE,IUNI,NPRED
C
	IF(ISAT.NE.0) READ(1,*) (IDSAT(I), I=1,ISAT)
C	write(*,*)1
C
C
	IF(ISAT.NE.0) READ(1,*) (NSAT(IDSAT(I)), I=1,ISAT)
	IF(ISAT.NE.0) READ(1,*) (IOBJSAT(IDSAT(I)), I=1,ISAT)
C
	NCXP = 0 !Cuenta el numero de datos de equilibrio
C
C
	DO I=1,ISAT
C
	  NSUMK = 0
	  do K=1,(IDSAT(I)-1)
	    NSUMK = NSUMK + 1
	  end do
C
	  NSUML = 0
	  do L=1,NSUMK
	    NSUML = NSUML + NSAT(L)
	  end do
C
	  DO J= (NSUML +1 ),( NSUML + NSAT(IDSAT(I)) )
	    READ(1,*) TEXP(J),PEXP(J)
	    write(*,*) TEXP(J),PEXP(J)
169	    FORMAT(F10.4,F16.10)
	  END DO
C
	  NCXP = NCXP + NSAT(IDSAT(I)) !Nro. de datos YA escritos
C
	END DO
C
C	write(*,*)2
	DO I=1,NG
	  READ(1,*) QAREA(I), TSTR(I)
	END DO
C
	READ(1,*) (GSTR(I), I=1,NG)
 138	FORMAT(<NG>F10.1)
C
	READ(1,*) (GDOT(I), I=1,NG)
C
	READ(1,*) (GDDOT(I), I=1,NG)
C
	DO I=1,NG
	  READ(1,*) (RKIJ(I,J), J=1,NG)
	END DO
C
	DO I=1,NG
	  READ(1,*) (RKTIJ(I,J), J=1,NG)
	END DO
C
C
	DO I=1,NG
	  READ(1,*) (ALFAIJ(I,J), J=1,NG)
	END DO
C
C
C	write(*,*)3
	IF(NGA.GT.0)  READ(1,*) (MASSOC(I),I=1,NGA)
C
	IF(NGA.GT.0)THEN
C
C LECTURA DE PARAMETROS ENERGETICOS DE ASOCIACION
C
C	DO J=1,NGA
C	      IF(MASSOC(J).EQ.0) GO TO 20
C		DO L=1,NGA
C                  IF(MASSOC(L).EQ.0) GO TO 10
C			DO I=1,MASSOC(J)
C				DO K=1,MASSOC(L)
C		                  READ(1,*) ENASS(I,J,K,L)
C				END DO
C			END DO
C 10	        CONTINUE
C	        END DO
C 20	      CONTINUE
C	END DO
C
C
C
C LECTURA DE PARAMETROS VOLUMETRICOS DE ASOCIACION
C
C	DO J=1,NGA
C	      IF(MASSOC(J).EQ.0) GO TO 30
C		DO L=1,NGA
C                 IF(MASSOC(L).EQ.0) GO TO 50
C			DO I=1,MASSOC(J)
C				DO K=1,MASSOC(L)
c				       READ(1,179) RKASS(I,J,K,L)
C				       READ(1,*) RKASS(I,J,K,L)
C 179 				       FORMAT(F16.10)
C				END DO
C			END DO
C 50	        CONTINUE
C	    END DO
C30	      CONTINUE
C       END DO
C
C       ----------------------------------------IF END
C	 ELSE
C	  CONTINUE
C	END IF
C       ----------------------------------------IF END
C
C-----------------------------------------------Francisco-----------------------------------
C
C	Lectura simplificada.
C

C	Lectura de energía y volumen de asociación en una misma línea
	  do i = 1,NGA
	    do j = i,NGA
	      do k = 1,Massoc(i)
	        do l = k,Massoc(j)
C
	          read(1,*)EnAss(k,i,l,j),rKAss(k,i,l,j)
C---------------Simetrías
	          EnAss(l,i,k,j) = EnAss(k,i,l,j)
	          EnAss(l,j,k,i) = EnAss(k,i,l,j)
	          EnAss(k,j,l,i) = EnAss(k,i,l,j)
	          rKAss(l,i,k,j) = rKAss(k,i,l,j)
	          rKAss(l,j,k,i) = rKAss(k,i,l,j)
		    rKAss(k,j,l,i) = rKAss(k,i,l,j)
C
	        enddo
	      enddo
	    enddo
	  enddo
	endif
C	WRITE(*,*)'FGTEGST   RTJUETYHETRUR  HR6R'
C-----------------------------------------------------------------------
C
	DO I=1,NC
c		READ(1,1) TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I)
	  READ(1,*) TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I)		!,DC(I)
c		write(*,*) TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I)
	END DO
C
	DO I=1,NG
	  READ(1,*) (NU(I,J), J=1,NC)
C
C NU (I,J)= NUMERO DE GRUPOS DEL TIPO I EN EL COMPUESTO DE TIPO J (DISPERSION
C
	END DO
C
C
C       ...........................................................
	IF(NGA.GT.0) THEN
C
	  DO I=1,NGA
	    READ(1,*) (ONUASOC(I,J), J=1,NC)
C 			ONUASOC(I,J)= NUMERO DE GRUPOS DEL TIPO I EN EL
C					COMPUESTO DE TIPO J (ASOCIACION
C
	  END DO
C	ELSE
C	  CONTINUE
	END IF
C       ...........................................................
C
C
	DO I=1,IGLE
	  write(*,*)'------------------------------------'
	  write(*,*)I
C
	  read(1,*)IHEN(I)
	  IF(IHEN(I).EQ.0) THEN
	    READ(1,*) NOBS(I),IDC(I,1),IDC(I,2)
	    READ(1,*) IFL(I)
	  ELSEif(ihen(i).eq.5)then
	    read(1,*) mez(i)
	    If (mez(i)==2) then
	      READ(1,*) IDC(I,1),IDC(I,2),NOBS(I)
	    elseif(mez(i)==3)then
	      READ(1,*) IDC(I,1),IDC(I,2),IDC(i,3),NOBS(I)
	    end if
	  END IF
C
        READ(1,23) (NTEXT(I,J),J=1,65)
23	  FORMAT(40A2)
C
	  DO J=( NCXP+1 ),(NCXP + NOBS(I))
	    If(ihen(i).eq.0)then
C
	      READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J)
	      write(*,9010) J,TEXP(J),PEXP(J),XEXP(J),YEXP(J)			!Impresión de datos en pantalla
9010	      format(I5,F7.2,3F7.4)
C
	    elseif(ihen(i).eq.5)then
C
	      if (mez(i)==2)then
	        read(1,*) TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J)
	      elseif(mez(i)==3)then
		  read(1,*) TEXP(J),PEXP(J),XEXP(J),xexp2(J),rgamexp1(J)
     @                  ,rgamexp2(J), rgamexp3(J)
	      endif
C
	    end if
	  END DO
	  NCXP = NCXP + NOBS(I)
C
	END DO !Fin del bucle sobre los DATASETS
C

	READ(1,*) NNP,NNQ,NALFA,NENASS,NRKASS,numtc, NUMDC
C	write(*,*)'NNP,NNQ,NALFA,NENASS,NRKASS,numtc, NUMDC'
C	write(*,*) NNP,NNQ,NALFA,NENASS,NRKASS,numtc, NUMDC

C
C
	if(NNP.ne.0) then
C
		READ(1,*) (IDG(I,1),IDG(I,2),I=1,NNP)
		write(*,*)(IDG(I,1),IDG(I,2),I=1,NNP)
C
		READ(1,*) (ITEMP(I), I=1,NNP)
		write(*,*)(ITEMP(I), I=1,NNP)
C
		READ(1,*) (NIDG(I), I=1,NNP)
		write(*,*)(NIDG(I), I=1,NNP)
C
		NSSS = 0
C
		DO L=1,NNP
C                  ----------------------------
		   if(NIDG(L).ne.0) then
C
			READ(1,*) (IDIDG(J,1), IDIDG(J,2),
     @J=(1+NSSS),(NIDG(L)+NSSS))
C			write(*,*)(IDIDG(J,1), IDIDG(J,2),
C    @J=(1+NSSS),(NIDG(L)+NSSS))
			NSSS = NSSS + NIDG(L)
C
				    else
					continue
		   end if
C                  -----------------------------
		END DO
C
	             else
				continue
	end if
C
	if(NNQ.ne.0)  READ(1,*) (IDQ(I), I=1,NNQ)
C
	if(NALFA.ne.0) then
C
		READ(1,*) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
C
		READ(1,*) (NIDALF(I),I=1,NALFA)
C
		NSSA = 0
C
		DO L=1,NALFA
C                   ------------------------------------------
                    if(NIDALF(L).ne.0) then
			READ(1,*) (IDIDALF(J,1),IDIDALF(J,2),
     @J=(1+NSSA),(NIDALF(L)+NSSA))
			NSSA = NSSA + NIDALF(L)
		  		       else
						continue
		    end if
C                   ------------------------------------------
		END DO
C
		       else
				continue
	end if
C
C - - - -  - - -  - - - - - - - - - - - - - - - - -
	if(NENASS.ne.0) then
C
	DO I=1,NENASS
C
	READ(1,*) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
C
	END DO
C
C
	READ(1,*) (NIDENASS(I), I=1,NENASS )
C
C
	NSSE = 0
	DO L=1,NENASS
C           ----------------------------------------------------
	    if(NIDENASS(L).ne.0) then
		DO J=(1+NSSE),(NIDENASS(L)+NSSE)
			READ(1,*) IDIDENASS(J,1), IDIDENASS(J,2),
     @ IDIDENASS(J,3), IDIDENASS(J,4)
 3			FORMAT(4I5)
		END DO
			NSSE = NSSE + NIDENASS(L)
				 else
					continue
	    end if
C           ----------------------------------------------------
	END DO
C
C
			else
				continue
	end if
C - - - -  - - -  - - - - - - - - - - - - - - - - -
C
C -  - - - - - - - - - - - - - - - - - - - - - -
C
	if(NRKASS.ne.0) then
C
	DO I=1,NRKASS
C
	READ(1,*) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
C
	END DO
C
C
	READ(1,*) (NIDRKASS(I), I=1,NRKASS )
C
	NSSRK = 0
	DO L=1,NRKASS
C           -----------------------------------------------------
	      if(NIDRKASS(L).ne.0) then
		DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)
			READ(1,*) IDIDRKASS(J,1), IDIDRKASS(J,2),
     @ IDIDRKASS(J,3), IDIDRKASS(J,4)
		END DO
			NSSRK = NSSRK + NIDRKASS(L)
				   else
					continue
	      end if
C           -----------------------------------------------------
	END DO
C
			else
				continue
	end if
C -  - - - - - - - - - - - - - - - - - - - - - -
C
C

	IF(numtc.NE.0) READ(1,*) (idtc(I), I=1,numtc)
c
	IF(NUMDC.NE.0) READ(1,*) (IDDC(I), I=1,NUMDC)
C
C NUMERO DE PARAMETROS
C
	NPAR = NNP+NNQ+NALFA+NENASS+NRKASS+numtc+NUMDC
	IF(NNP.NE.0) THEN
		DO I = 1,NNP
	  	 NPAR = NPAR + ITEMP(I)
		END DO
	 ELSE
	  CONTINUE
	END IF
C
C LECTURA DE LOS VALORES INICIALES DE LOS PARAMETROS
C
	READ(1,*)  (XX(I),I=1,NPAR)

C
	CLOSE(UNIT=1)
C
	END !FIN READ
C/////////////////////////////////////////////////////////////////
C
C
	SUBROUTINE WRITESAL1(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,W1,W2,W3,W4,W5,
     @ w6,W7,W8, NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,
     @NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,
     @IFL,NTEXT,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,XX,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @IDRKASS,NIDRKASS,IDIDRKASS,DC,IOPT,NUMDC,IDDC,TSPLINE,numtc,idtc)
C
C
C
C ESCRIBE DATOS ARCHIVO DE SALIDA paraest.out
C
	IMPLICIT REAL*8 (A-H,O-Z)
C
C  DIM DISPERSION
C
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
	DIMENSION TSPLINE(10)
C
C  DIM ASOCIACION
C
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,20)
C
C  DIM DIAMETRO DE ESFERAS RIGIDAS EN EL PUNTO CRITICO
C
	DIMENSION DC(20)
C
C  DIM AJUSTE DE PARAM
C
	DIMENSION PARM(4)
C
	DIMENSION IDSAT(20)
C
	DIMENSION  TEMPCRIT(20)
	DIMENSION  PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260), xexp2(260)
C
C DIMENSIONAMIENTO  NRO. DE DATOS DE PVAP COMP PUROS
C
	DIMENSION NSAT(20)
C
C DIMENSIONAMIENTOS QUE TIENEN QUE VER CON LOS DATASETS
C
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
C
C ID ==> identificacion o identicos
C IDID ==> identificacion de identicos
C
C   los parametros se identifican a traves de su naturaleza y de sus
C	inidices
C
	DIMENSION ITEMP(8)
	DIMENSION IDG(8,2) ! identificacion (indices) de param gii o kij a
C									ajust.
	DIMENSION NIDG(8)   !nro. de param. IDENTicos a c/gii o kij
	DIMENSION IDIDG(104,2) !identificacion (indices) de IDENTicos a gii o
C									  kij
C
C                8 X 13 = 104
C
	DIMENSION IDQ(8)  !IDENTif param de area
C
	DIMENSION IDALF(8,2) !identif (indices) de parametros ALFAij
	DIMENSION NIDALF(8) !nro. de param. IDENTicos a c/alfaij
	DIMENSION IDIDALF(104,2) !identificacion (indices) de IDENTicos
C
	DIMENSION IDENASS(8,4) !IDENTif (indices) de parametros ENASSijkl
	DIMENSION NIDENASS(8) !Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION IDIDENASS(104,4) !IDENTif (indices) de IDENTicos a ENASSijkl
C
	DIMENSION IDRKASS(8,4) !Identif (indices) de param. RKASSijkl
	DIMENSION NIDRKASS(8) !Nro. de param. identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4) !Identif (indices) de identicos a c/RKASSijkl
C
	DIMENSION IDDC(8)
c
	dimension idtc(8)
C
	DIMENSION XX(8)
C
	CHARACTER*9 FECHA,HORA
C	CALL DATE(FECHA)
C	CALL TIME(HORA)
C
c	OPEN(UNIT=1,TYPE='OLD',NAME='ASSOCPARA.RES',FORM='FORMATTED',
c     @ ACCESS='SEQUENTIAL')
	OPEN(UNIT=1,STATUS='unknown',file='paraest.out',FORM='FORMATTED'
     *)
C
C
C	WRITE(1,415) FECHA,HORA
C 415	FORMAT(/,X,'                           FECHA =',A,/,
C     @X,'                           HORA  =',A,/,X,
C     @' ESTIMACION DE PARAMAETROS CON EL MODELO ASSOC-GC-EOS')
C
	WRITE(1,1) W1,W2,W3,W4,W5,w6,W7,W8,PARM(1),PARM(2),cc
C
 1	FORMAT(/,X,' FACTORES DE PESO UTILIZADOS EN LA FUNCION
     @OBJETIVO ',
     @/,X,' WSAT,WHEN,WACT,WVAP o Wcomp2,WVE: ',/,X,8F10.4,
     @/,/,X,' PARM(1) = ',F10.4,'     PARM(2) = ',F10.4,/,X,
     @' ALGORITMO DE OPTIMIZACION: IOPT = ',I3)
C
	WRITE(1,2)
 2	FORMAT(/,X,' RELACION ENTRE g*ii y gii A LA QUE SE REALIZA ',/,X,
     @' EL SPLINE = 4 ')
C
	WRITE(1,47)
 47	FORMAT(/,X,'--------------------------------------------------',
     @/,X,'------------ PARAMETROS INICIALES ----------------',
     @/,X,'--------------------------------------------------',
     @/,X,' PARAMETROS DE GRUPOS PUROS (termino  dispersivo)',/)
C
	WRITE(1,3) (I, I=1,NC)
 3	FORMAT(40X,'COMPONENTES',/,X,
     @' Gr.Nro.:  T*        Q        TSPLINE ',<NC>I3,/)
C
	DO I=1,NG
		WRITE(1,4) I,TSTR(I),QAREA(I),TSPLINE(I),(NU(I,J), J=1,NC)
 4		FORMAT(X,I3,5X,2F10.4,X,E10.4,<NC>I3)
	END DO
C
C
C NU (I,J)= NUMERO DE GRUPOS DEL TIPO I EN EL COMPUESTO DE TIPO J
C
	WRITE(1,5)
 5	FORMAT(/,X,' ENERGIAS DE INTERACCION gjj ',
     @/,X,'(la forma funcional cambia a altas temperaturas)',/)
C
	DO I=1,NG
C
	WRITE(1,6) I,I,GSTR(I), GDOT(I), GDDOT(I)
 6	FORMAT(X,' g(',I3,I3,') = ',
     @F10.1,' * [ 1 + ',F10.4,' *(T/T* - 1) + ',
     @F10.4, ' *LN(T/T*)] ',/)
C
	END DO
C
	WRITE(1,7)
 7	FORMAT(/,X,' PARAMETROS BINARIOS k*ij ',/)
C
C
	DO I=1,NG
		WRITE(1,8) (RKIJ(I,J), J=1,NG)
 8		FORMAT(X,<NG>F10.4)
	END DO
C
	WRITE(1,9)
 9	FORMAT(/,X,' DEPENDENCIA RESP. A LA TEMP. : k´ij ',/)
C
C
	DO I=1,NG
		WRITE(1,8) (RKTIJ(I,J), J=1,NG)
	END DO
C
	WRITE(1,10)
 10	FORMAT(/,X,' PARAMETROS BINARIOS DE NO AZAR alfaij ',/)
C
C
	DO I=1,NG
		WRITE(1,8) (ALFAIJ(I,J), J=1,NG)
	END DO
C
	WRITE(1,186) NGA
 186	FORMAT(/,X,' Numero de grupos considerados en la contribucion
     @de asociacion = ',I3)
C
C
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IF(NGA.GT.0) THEN
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
	WRITE(1,11) (I, I=1,NC)
 11	FORMAT(/,X,' PARAMETROS DE ASOCIACION ',/,23X,'COMPONENTES',/,
     @' GRUPO 	  # SIT ASOC  ',<NC>I5,/)
C
	DO I=1,NGA
	WRITE(1,12) I, MASSOC(I), (ONUASOC(I,J), J=1,NC)
 12	FORMAT(3X,I3,9X,I3,6X,<NC>f5.2)
	END DO
C
	WRITE(1,13)
 13 	FORMAT(/,X,'PARAMETROS DE ENERGIA DE ASOCIACION (Kelvin)  ',/)
C
C
C ESCRITURA DE PARAMETROS ENERGETICOS DE ASOCIACION
C
	DO J=1,NGA
	      IF(MASSOC(J).EQ.0) GO TO 20
		DO L=1,NGA
                  IF(MASSOC(L).EQ.0) GO TO 101
			DO I=1,MASSOC(J)
				DO K=1,MASSOC(L)
		                  WRITE(1,14) I,J,K,L,ENASS(I,J,K,L)
 14	FORMAT(X,' ENASS( ',I3,I3,I3,I3,' )  = ',F10.4)
				END DO
			END DO
 101	        CONTINUE
	        END DO
 20	      CONTINUE
	END DO
C
C
C ESCRITURA DE PARAMETROS VOLUMETRICOS DE ASOCIACION (cm3/mol)
C
	WRITE(1,15)
 15	FORMAT(/,X,'PARAMETROS DE VOLUMEN DE ASOCIACION (cm3/mol) ',/)
C
	DO J=1,NGA
	      IF(MASSOC(J).EQ.0) GO TO 30
		DO L=1,NGA
                  IF(MASSOC(L).EQ.0) GO TO 50
			DO I=1,MASSOC(J)
				DO K=1,MASSOC(L)
				    WRITE(1,16) I,J,K,L,RKASS(I,J,K,L)
 16	FORMAT(X,' RKASS( ',I3,I3,I3,I3,' ) = ',F16.10)
				END DO
			END DO
 50	        CONTINUE
	        END DO
 30	      CONTINUE
	END DO
C
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 ELSE
	  CONTINUE
	END IF
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
	WRITE(1,18)
 18	FORMAT(/,X,'DATOS SOBRE LOS COMPUESTOS PUROS')
C
	WRITE(1,19)
 19	FORMAT(/,4X,'Tc/K      Pc/atm    Tsat/K    Psat/atm    DC'
     @,/,X,' (DC inicial estimado en base al pto critico exptl.) ',/)

C
	DO I=1,NC
		WRITE(1,500) TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I),
     @DC(I)
 500		FORMAT(X,4F10.4,F10.7)
	END DO
C
	WRITE(1,201)
 201	FORMAT(/,X,' PARAMETROS A SER ESTIMADOS ',/)
C
C
	if(NNP.ne.0) then
C
		WRITE(1,202)
 202		FORMAT(' Parametros gii o kij ')
C
		WRITE(1,203) (IDG(I,1),IDG(I,2),I=1,NNP)
 203		FORMAT(X,<2*NNP>I3)
C
		WRITE(1,204) (ITEMP(I), I=1,NNP)
 204		FORMAT(X,' Dependencia respecto a la temperatura ',
     @/,X,<NNP>I3)
C
		WRITE(1,206) (NIDG(I), I=1,NNP)
 206		FORMAT(X,' Numero de par. identicos ',/,X,<NNP>I3)
C
		WRITE(1,208)
 208		FORMAT(X,' Parametros Identicos ')
C
		NSSS = 0
C
		DO L=1,NNP
C                  ----------------------------
		   if(NIDG(L).ne.0) then
C
			WRITE(1,210) (IDIDG(J,1), IDIDG(J,2),
     @J=(1+NSSS),(NIDG(L)+NSSS))
 210			FORMAT(X,<2*NIDG(L)>I3)
			NSSS = NSSS + NIDG(L)
C
				    else
					continue
		   end if
C                  -----------------------------
		END DO
C
	             else
				continue
	end if
C
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if(NNQ.ne.0)  WRITE(1,211) (IDQ(I), I=1,NNQ)
 211	FORMAT(/,X,' Parametros de Area  ',/,X,<NNQ>I3)
C
	if(NALFA.ne.0) then

C
		WRITE(1,213) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
 213		FORMAT(/,X,' Parametros de no azar ',/,X,<2*NALFA>I3)
C
		WRITE(1,215) (NIDALF(I),I=1,NALFA)
 215		FORMAT(X,' Numero de par. identicos ',/,X,<NALFA>I3)
C
		NSSA = 0
C
		DO L=1,NALFA
C                   ------------------------------------------
                    if(NIDALF(L).ne.0) then
			WRITE(1,217) (IDIDALF(J,1),IDIDALF(J,2),
     @J=(1+NSSA),(NIDALF(L)+NSSA))
 217			FORMAT(X,' Identicos : ',<2*NIDALF(L)>I3)
			NSSA = NSSA + NIDALF(L)
		  		       else
						continue
		    end if
C                   ------------------------------------------
		END DO
C
		       else
				continue
	end if
C
C - - - -  - - -  - - - - - - - - - - - - - - - - -
	if(NENASS.ne.0) then
C
	WRITE(1,219)
 219	FORMAT(/,X,' Parametros de energia de asociacion ')
C
	DO I=1,NENASS
C
	WRITE(1,221) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
 221	FORMAT(X,4I3)
C
	END DO
C
C
	WRITE(1,223) (NIDENASS(I), I=1,NENASS )
 223	FORMAT(X,' Numero de par. identicos ',<NENASS>I3)
C
C
	NSSE = 0
	DO L=1,NENASS
C           ----------------------------------------------------
	    if(NIDENASS(L).ne.0) then
		DO J=(1+NSSE),(NIDENASS(L)+NSSE)
			WRITE(1,225) IDIDENASS(J,1), IDIDENASS(J,2),
     @ IDIDENASS(J,3), IDIDENASS(J,4)
 225			FORMAT(' Identicos: ',4I5)
		END DO
			NSSE = NSSE + NIDENASS(L)
				 else
					continue
	    end if
C           ----------------------------------------------------
	END DO
C
C
			else
				continue
	end if
C - - - -  - - -  - - - - - - - - - - - - - - - - -
C
C -  - - - - - - - - - - - - - - - - - - - - - -
C
	if(NRKASS.ne.0) then
C
C
	WRITE(1,503)
 503	FORMAT(/,X,' Parametros de Volumen de Asociacion ')
C
	DO I=1,NRKASS
C
	WRITE(1,229) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
 229	FORMAT(X,4I3)
C
	END DO
C
C
	WRITE(1,231) (NIDRKASS(I), I=1,NRKASS )
 231	FORMAT(X,' Numero de parametros identicos ',<NRKASS>I3)
C
	NSSRK = 0
	DO L=1,NRKASS
C           -----------------------------------------------------
	      if(NIDRKASS(L).ne.0) then
		DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)
			WRITE(1,232) IDIDRKASS(J,1), IDIDRKASS(J,2),
     @ IDIDRKASS(J,3), IDIDRKASS(J,4)
 232		FORMAT(X,' Identicos: ',4I3)
		END DO
			NSSRK = NSSRK + NIDRKASS(L)
				   else
					continue
	      end if
C           -----------------------------------------------------
	END DO
C
			else
				continue
	end if
C -  - - - - - - - - - - - - - - - - - - - - - -

	IF(numtc.ne.0) WRITE(1,233) (idtc(I), I=1,numtc)
  233   FORMAT(/,X,' temperaturas criticas ',/,X,
     @<numtc>I3)

	IF(NUMDC.NE.0) WRITE(1,234) (IDDC(I), I=1,NUMDC)
  234   FORMAT(/,X,' Diametros de Esferas Rigidas ',/,X,
     @<NUMDC>I3)
C
C
	WRITE(1,235)  (XX(I),I=1,NPAR)
 235	FORMAT(/,X,' Valores iniciales de los parametros a ajustar '
     @,/,X,<NPAR>F13.4)
C
	WRITE(1,91)
 91	FORMAT(//,X,' Curso de la optimizacion ',X)
	WRITE(1,93) (I,I=1,NPAR)
 93	FORMAT(//,X,' PARAMETRO NRO (ESCALADO). ',/,X,<NPAR>I15,X,'
     @         FC OBJ ',/)
C
	CLOSE(UNIT=1)
C
	END !FIN WRITESAL1
C/////////////////////////////////////////////////////////////////
C
C
C
C
C
C
C	DO I=1,IGLE
C
C		WRITE(1,2) IDC(I,1),IDC(I,2),NOBS(I),IHEN(I)
C
C		IF(IHEN(I).EQ.0) THEN
C					WRITE(1,2) IFL(I)
C				 ELSE
C	 	END IF
C
C       		WRITE(1,23) (NTEXT(I,J),J=1,65)
C 23 		FORMAT(40A2)
C
C		DO J=( NCXP+1 ),(NCXP + NOBS(I))
C			WRITE(1,1) TEXP(J),PEXP(J),XEXP(J),YEXP(J)
C		END DO
C		NCXP = NCXP + NOBS(I)
C
C	END DO !Fin del bucle sobre los DATASETS
C
C	WRITE(1,2) NNP,NNQ,NALFA,NENASS,NRKASS
C
C	WRITE(1,2) NPAR
C
C FIN RESERVA
C///////////////////////////////////////////////////////////////////////////
	SUBROUTINE MODIFPAR(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,
     @  W1,W2,W3,W4,W5, w6,W7,W8,
     @  NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,NRKASS,
     @  NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,
     @  NTEXT,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,XXAUX,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @	IDRKASS,NIDRKASS,IDIDRKASS,numtc,idtc,NUMDC,
     @  IDDC,DC)
C
C MODIFICA LOS VALORES DE LOS PARAMETROS DEL MODELO EN BASE A LOS
C  NUEVOS VALORES DE LOS PARAMETROS QUE SE AJUSTAN. ESTOS ULTIMOS
C    ESTAN CONTENIDOS EN XXAUX, LIBRES DE ESCALADO
C
C  A ESTA ALTURA NPAR ES DATO
C
	IMPLICIT REAL*8 (A-H,O-Z)
C
C  DIM DISPERSION
C
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
C
C  DIM ASOCIACION
C
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,20)
C
C  DIM AJUSTE DE PARAM
C
	DIMENSION PARM(4)
C
	DIMENSION IDSAT(20)
C
	DIMENSION TEMPCRIT(20)
	DIMENSION PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260),xexp2(260)
C
C DIMENSIONAMIENTO  NRO. DE DATOS DE PVAP COMP PUROS
C
	DIMENSION NSAT(20)
C
C DIMENSIONAMIENTOS QUE TIENEN QUE VER CON LOS DATASETS
C
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
C
C ID ==> identificacion o identicos
C IDID ==> identificacion de identicos
C
C   los parametros se identifican a traves de su naturaleza y de sus
C	inidices
C
	DIMENSION ITEMP(8)
	DIMENSION IDG(8,2) ! identificacion (indices) de param gii o kij a
C									ajust.
	DIMENSION NIDG(8)   !nro. de param. IDENTicos a c/gii o kij
	DIMENSION IDIDG(104,2) !identificacion (indices) de IDENTicos a gii o
C									  kij
C
C                8 X 13 = 104
C
	DIMENSION IDQ(8)  !IDENTif param de area
C
	DIMENSION IDALF(8,2) !identif (indices) de parametros ALFAij
	DIMENSION NIDALF(8) !nro. de param. IDENTicos a c/alfaij
	DIMENSION IDIDALF(104,2) !identificacion (indices) de IDENTicos
C
	DIMENSION IDENASS(8,4) !IDENTif (indices) de parametros ENASSijkl
	DIMENSION NIDENASS(8) !Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION IDIDENASS(104,4) !IDENTif (indices) de IDENTicos a ENASSijkl
C
	DIMENSION IDRKASS(8,4) !Identif (indices) de param. RKASSijkl
	DIMENSION NIDRKASS(8) !Nro. de param. identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4) !Identif (indices) de identicos a c/RKASSijkl
C
C
C DIMENSIONAMIENTO DIAMETRO DE ESFERAS RIGIDAS EN EL PUNTO CRIT
C
	dimension idtc(8)

	DIMENSION DC(20)
	DIMENSION IDDC(8)
C
	DIMENSION XXAUX(8)
C
C NPAR2 ES UNA VARIABLE QUE TIENE EN CUENTA EL NUMERO DE PARAMETROS YA
C  ACTUALIZADOS. UNA VEZ QUE TERMINA LA EJECUCION DE ESTE MODULO NPAR2 DEBE
C   SER IGUAL A NPAR
C
	NPAR2 = 0
C
C ASIGNACION DE VALORES A PARAMETROS DEL TIPO gii o kij
C
	IF(NNP.ne.0) then !A
C
		DO I=1,NNP
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C
	           IF(IDG(I,1).EQ.IDG(I,2)) THEN
C
			GSTR(IDG(I,1)) = XXAUX(NPAR2+1)
	           	NPAR2 = NPAR2 + 1
			IF(ITEMP(I).EQ.0) GDOT(IDG(I,1)) = 0.
			IF(ITEMP(I).EQ.0) GDDOT(IDG(I,1)) = 0.
		   	IF(ITEMP(I).GE.1) GDOT(IDG(I,1)) = XXAUX(NPAR2+1)
		   	IF(ITEMP(I).GE.1) NPAR2 = NPAR2 + 1
		   	IF(ITEMP(I).EQ.1) GDDOT(IDG(I,1)) = 0.
			IF(ITEMP(I).EQ.2) GDDOT(IDG(I,1)) = XXAUX(NPAR2+1)
			IF(ITEMP(I).EQ.2) NPAR2 = NPAR2 + 1
C
					    ELSE
C
			RKIJ(IDG(I,1),IDG(I,2)) = XXAUX(NPAR2+1)
		        NPAR2 = NPAR2 + 1
			IF(ITEMP(I).EQ.0) RKTIJ(IDG(I,1),IDG(I,2))= 0.
 1001			IF(ITEMP(I).EQ.1) RKTIJ(IDG(I,1),IDG(I,2))=
     @							XXAUX(NPAR2+1)
			IF(ITEMP(I).EQ.1) NPAR2 = NPAR2 + 1
			IF(ITEMP(I).EQ.2) THEN  !avoid human error
						ITEMP(I) = 1
						GO TO 1001
			 ELSE
			END IF
C
		   END IF
C
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		END DO
	else !A
 				continue
	end if !A
C
C IDENTICOS DEL DEL TIPO gii o kij
C
C
	IF(NNP.NE.0) THEN
		NSSS = 0
C
		DO I=1,NNP
C                  ----------------------------
		   if(NIDG(I).ne.0) then
C                    ------------------------------
	             IF(IDG(I,1).EQ.IDG(I,2)) THEN
C
			  DO J=(1+NSSS),(NIDG(I)+NSSS)
			    GSTR(IDIDG(J,1))  = GSTR(IDG(I,1))
		   	    GDOT(IDIDG(J,1))  = GDOT(IDG(I,1))
			    GDDOT(IDIDG(J,1)) = GDDOT(IDG(I,1))
			  END DO
			  NSSS = NSSS + NIDG(I)
C
		     ELSE
C
			  DO J=(1+NSSS),(NIDG(I)+NSSS)
			    RKIJ(IDIDG(J,1),IDIDG(J,2)) =
     @ 					RKIJ(IDG(I,1),IDG(I,2))
			    RKTIJ(IDIDG(J,1),IDIDG(J,2))=
     @ 					RKTIJ(IDG(I,1),IDG(I,2))
			  END DO
			  NSSS = NSSS + NIDG(I)
C
		     END IF
C                    ------------------------------
C
	           else
			continue
		   end if
C                  -----------------------------
		END DO
C
	ELSE
	  CONTINUE
	END IF
C
C
C ASIGNACION DE VALORES A PARAMETROS DE AREA
C
C
C
	IF(NNQ.NE.0) THEN
	   DO I=1,NNQ
		QAREA(IDQ(I)) = XXAUX(NPAR2+1)
		NPAR2 = NPAR2 + 1
	   END DO
	ELSE
 	   CONTINUE
	END IF
C
C ASIGNACION DE VALORES A PARAMETROS ALFA
C
	IF(NALFA.NE.0) THEN
		DO I=1,NALFA
			ALFAIJ(IDALF(I,1),IDALF(I,2)) = XXAUX(NPAR2+1)
			NPAR2 = NPAR2 + 1
		END DO
C
	ELSE
          CONTINUE
	END IF
C
C ASIGNACION A IDENTICOS A LOS ALFA
C
	IF(NALFA.NE.0) THEN
C
		NSSA = 0
C
		DO I=1,NALFA
C                   ------------------------------------------
                    if(NIDALF(I).ne.0) then
C
			DO J=(1+NSSA),(NIDALF(I)+NSSA)
			 ALFAIJ(IDIDALF(J,1),IDIDALF(J,2))=
     @ALFAIJ(IDALF(I,1),IDALF(I,2))
			END DO
			NSSA = NSSA + NIDALF(I)
C
		    else
		     continue
		    end if
C                   ------------------------------------------
		END DO

	ELSE
	  CONTINUE
	END IF
C
C ASIGNACION DE VALORES A PARAMETROS DE ENERGIA DE ASOCIACION
C
C
C
	IF(NENASS.NE.0) THEN
C
		DO I=1,NENASS
C
		  ENASS(IDENASS(I,1),IDENASS(I,2),IDENASS(I,3)
     @,IDENASS(I,4))
     @		  = XXAUX(NPAR2+1)
		  NPAR2 = NPAR2 + 1
C
 		END DO
C
	ELSE
	 	CONTINUE
	END IF
C
C ASIGNACION A IDENTICOS TIPO ENERGIA DE ASOCIACION
C
C
	IF(NENASS.NE.0) THEN
C
	NSSE = 0
	DO L=1,NENASS
C           ----------------------------------------------------
	    if(NIDENASS(L).ne.0) then
		DO J=(1+NSSE),(NIDENASS(L)+NSSE)
		        ENASS(IDIDENASS(J,1),IDIDENASS(J,2)
     @                  ,IDIDENASS(J,3),IDIDENASS(J,4)) =
     @ 	                ENASS(IDENASS(L,1),IDENASS(L,2)
     @                  ,IDENASS(L,3),IDENASS(L,4))
		END DO
			NSSE = NSSE + NIDENASS(L)
				 else
					continue
	    end if
C           ----------------------------------------------------
	END DO
C
	ELSE
	 	CONTINUE
	END IF
C
	if(NRKASS.ne.0) then
C
		DO I=1,NRKASS
C
		RKASS(IDRKASS(I,1),IDRKASS(I,2),
     @IDRKASS(I,3),IDRKASS(I,4)) = XXAUX(NPAR2+1)
		NPAR2 = NPAR2 + 1
C
		END DO
C
	else
	 continue
	end if
C
C ASIGNACION DE VALORES DE PARAMETROS DE VOLUMEN DE ASOCIACION
C  IDENTICOS A LOS YA ESTABLECLIDOS
C
	IF(NRKASS.ne.0) THEN
C
	NSSRK = 0
	 	DO L=1,NRKASS
C           -----------------------------------------------------
		      if(NIDRKASS(L).ne.0) then
			DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)
				RKASS(IDIDRKASS(J,1),IDIDRKASS(J,2),
     @ IDIDRKASS(J,3),IDIDRKASS(J,4)) = RKASS(IDRKASS(L,1),IDRKASS(L,2),
     @IDRKASS(L,3),IDRKASS(L,4))
C
			END DO
			NSSRK = NSSRK + NIDRKASS(L)
					   else
						continue
		      end if
C           -----------------------------------------------------
		END DO
	ELSE
	 CONTINUE
	END IF
C
C ASIGNACION DE VALORES A LOS PARAMETROS DC, CUANDO DC SE USA COMO
C  PARAMETRO AJUSTABLE (3RA OPCION DEL PRIMER COMENTARIO DEL PROG.PPAL)
C
C
    	IF(numtc.NE.0) THEN
	   DO I=1,numtc
		tempcrit(idtc(I)) = XXAUX(NPAR2+1)
		NPAR2 = NPAR2 + 1
	   END DO
	 ELSE
 	   CONTINUE
	END IF
C

	IF(NUMDC.NE.0) THEN
	   DO I=1,NUMDC
		DC(IDDC(I)) = XXAUX(NPAR2+1)
		NPAR2 = NPAR2 + 1
	   END DO
	 ELSE
 	   CONTINUE
	END IF
C
	END !FIN MODIFPAR
C
C////////////////////////////////////////////////////////////////////////////
	SUBROUTINE FUNC(XX,NCXP,NPAR,FCION)
C
C  EN FCION APARECEN LOS VALORES A PARTIR DE LOS CUALES SE ARMA LA
C SUMATORIA DE CUADRADOS
C
	IMPLICIT REAL*8 (A-H,O-Z)
C
C  DIM DISPERSION
C
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
C
C  DIM ASOCIACION
C
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,20)
C
C DIMENSIONAMIENTO DIAMETRO DE ESFERAS RIGIDAS EN EL PUNTO CRIT
C
	DIMENSION DC(20)
C
C  DIM AJUSTE DE PARAM
C
C
	DIMENSION IDSAT(20)
C
	DIMENSION TEMPCRIT(20)
	DIMENSION PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260), xexp2(260)
	dimension rgamexp1(260),rgamexp2(260),rgamexp3(260)
	DIMENSION PCALC(260), XCALC(260), YCALC(260)
	dimension rgamcalc1(260),rgamcalc2(260),rgamcalc3(260)
	dimension zlicalc(260)
	DIMENSION ZLCALC(260), ZVCALC(260)
C
C VIENE IMPORTANTE DIMENSION VARIABLE
C
	DIMENSION FCION(NCXP)
	DIMENSION FCAUX(260)
	DIMENSION FCAUX2(260)
	DIMENSION FCAUX3(260)
C
C
C DIMENSIONAMIENTO  NRO. DE DATOS DE PVAP COMP PUROS
C
	DIMENSION NSAT(20)
	DIMENSION IOBJSAT(20)
C
C DIMENSIONAMIENTOS QUE TIENEN QUE VER CON LOS DATASETS
C
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION IFL(25)
C
C ID ==> identificacion o identicos
C IDID ==> identificacion de identicos
C
C   los parametros se identifican a traves de su naturaleza y de sus
C	indices
C
	DIMENSION ITEMP(8)
	DIMENSION IDG(8,2) ! identificacion (indices) de param gii o kij a
C									ajust.
	DIMENSION NIDG(8)   !nro. de param. IDENTicos a c/gii o kij
	DIMENSION IDIDG(104,2) !identificacion (indices) de IDENTicos a gii o
C									  kij
C
C                8 X 13 = 104
C
	DIMENSION IDQ(8)  !IDENTif param de area
C
	DIMENSION IDALF(8,2) !identif (indices) de parametros ALFAij
	DIMENSION NIDALF(8) !nro. de param. IDENTicos a c/alfaij
	DIMENSION IDIDALF(104,2) !identificacion (indices) de IDENTicos
C
	DIMENSION IDENASS(8,4) !IDENTif (indices) de parametros ENASSijkl
	DIMENSION NIDENASS(8) !Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION IDIDENASS(104,4) !IDENTif (indices) de IDENTicos a ENASSijkl
C
	DIMENSION IDRKASS(8,4) !Identif (indices) de param. RKASSijkl
	DIMENSION NIDRKASS(8) !Nro. de param. identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4) !Identif (indices) de identicos a c/RKASSijkl
C
	DIMENSION idtc(8)

      DIMENSION IDDC(8)

C
C LA DIMENSION VARIABLE QUE VIENE ES MUY IMPORTANTE QUE SEA VARIABLE
C
	DIMENSION XX(NPAR), IDXX(8) , XXAUX(8)
C
C DIMENSIONAMIENTO DE ARREGLOS DE TRABAJO NECESARIOS EN EL CALCULO DE
C  PROPIEDADES DE COMPUESTOS PUROS
C
	DIMENSION WDC(1)
	DIMENSION  WTEMPCRIT(1), WPRCRIT(1)
	DIMENSION NUW(10,1)
      DIMENSION ONUASOCW(12,1)
C
	DIMENSION PSTCLC(20)    !PRESION DE VAPOR(ATM) CALCULADA A LA
C				  TEMPERATURA TSAT(I)
C
C DIMENSIONAMIENTO DE ARREGLOS DE TRABAJO NECESARIOS EN EL CALCULO DE
C  PROPIEDADES DE MEZCLAS BINARIAS
C
	DIMENSION WWDC(3)
	DIMENSION  WWTEMPCRIT(3), WWPRCRIT(3)
	DIMENSION NUWW(10,3)
      DIMENSION ONUASOCWW(12,3)
C
C
C DIMENSIONAMIENTO DESVIACIONES STANDAR PARA MEZCLAS
C
	DIMENSION STDVA(25)
	DIMENSION STDVB(25)
	DIMENSION STDVc(25)
	real :: rgam1, rgam2, rgam3
	DIMENSION  :: PM(20),WWPM(2)

	COMMON/FORFUNC/NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,MASSOC,ENASS,RKASS,W1,W2,W3,W4,W5,
     @w6,W7,W8,ISAT,IGLE,IVE,IUNI,NNP,NNQ,NALFA,NENASS,NRKASS,
     @  NSAT,IDC,NOBS,IHEN,mez,IFL,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @	IDRKASS,NIDRKASS,IDIDRKASS,DC,PSTCLC,XCALC,YCALC,ZLCALC,
     @  ZVCALC,FCAUX,PCALC,IOBJSAT,FCAUX2,fcaux3,STDVA,STDVB,stdvc,
     @numtc,idtc,NUMDC,IDDC,NGA
	COMMON/FORFUNC2/rgamcalc1, rgamcalc2,rgamcalc3,
     @ rgamexp1, rgamexp2, rgamexp3
	COMMON/FORFUNC3/ONUASOC,zlicalc,PRCRIT,TEMPCRIT,TSAT,PSAT

	COMMON/IDSC/IDXX
C
	COMMON/ESCALADO/FGGG,FRKIJ,FTEMP,FQAREA,
     @FALFA,FENASS,FRKASS,ftc,FDC
C
	PENALIZ = 0.
C
C
C
C DES-ESCALADO DE LOS VALORES DE XX PARA PODER REALIZAR EL CALCULO
C   DE LOS TERMINOS DE LA FUNCION OBJETIVO. EL RESULTADO DEL DESESCA
C     -LADO SE GUARDA EN XXAUX, DADO QUE LA RUTINA FUNC NO DEBE ALTERAR
C       LOS VALORES DE XX (VER INSTRUCCIONES DE LA LIBRERIA IMSL)
C
	CALL XXDESSCAL(NPAR,XX,IDXX,FGGG,FRKIJ,FTEMP,FQAREA,
     @FALFA,FENASS,FRKASS,XXAUX,ftc,FDC)
C
C MODIFICACION DE LOS PARAMETROS EN BASE A LOS NUEVOS VALORES
C   DE LOS MISMOS, QUE ENTRARON ESCALADOS EN EN EL VECTOR XX, Y SE
C    SE ENCUENTRAN LIBRES DE ESCALAMIENTO EN EL VECTOR XXAUX
C
 	CALL MODIFPAR(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,
     @  W1,W2,W3,W4,W5,w6,W7,W8,
     @  NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,NRKASS,
     @  NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,
     @  NTEXT,
     @	TEXP,PEXP,XEXP,xexp2,YEXP,XXAUX,IDSAT,IDG,ITEMP,NIDG,IDIDG,
     @	IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,
     @	IDRKASS,NIDRKASS,IDIDRKASS,numtc,idtc,NUMDC,
     @  IDDC,DC)
C
C	COMIENZA EVALUACION DE FCION
C
C CALCULO DE DC (de manera de ajustar exactamente un punto de
C  la curva de presion de vapor, para cada componente
C
C Correccion de dc de manera de ajustar un punto de la curva de pvapor
C
c		OPEN(UNIT=5,TYPE='OLD',NAME='DIAM.TXT')

      DO I=1,NC
*****************************
*	  if(DC(I).EQ.0) then ! Fran. 26/05/09
*****************************
 	  IF(TSAT(I).NE.0) THEN
C         ------------------------------------------------------------
C
C	    Verifico si el DC(I) no ha sido ya considerado como una variable
C	     de optimizacion
C
	    DO JAM=1,NUMDC
		IF(IDDC(JAM).EQ.I)  GO TO 1325
	    END DO
C
C	    Carga de las variables de trabajo
	    WDC(1)= DC(I)
	    WTEMPCRIT(1)= TEMPCRIT(I)
 	    WPRCRIT(1)  = PRCRIT(I)
	     DO J=1,NG
		NUW(J,1) = NU(J,I)
	     END DO
C
	     DO K=1,NGA
		ONUASOCW(K,1) = ONUASOC(K,I)
	     END DO
C
	     WTEMP = TSAT(I)
	     WPSAT = PSAT(I)
C
C            Ajuste de DC
C
C
C
	     CALL ASSOCDCADJUSTFASTSUB(1,WDC,WTEMPCRIT,NG,NUW,QAREA,TSTR,
     @     GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCW,MASSOC,ENASS,
     @       RKASS,
     @       WTEMP,WPSAT,20,IERRDC2)
C
C			 20 es el numero maximo de iteraciones
C			  en ASSOCDCADJUSTFASTSUB
C
C            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^ ^^^^
C
	     IF(IERRDC2.NE.0) THEN
C
	       WDC(1)= DC(I)
C
C              Ajuste de DC por un metodo un poco mas robusto pero
C		unas 5 veces mas lento
C
C
	       CALL ASSOCDCADJUSTSUB(1,WDC,WTEMPCRIT,NG,NUW,QAREA,TSTR,
     @	       GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCW,
     @	       MASSOC,ENASS,RKASS,
     @         WTEMP,WPSAT,2,PSATCALC)
C
	       PSTCLC(I) = PSATCALC
C
C                        2  es el numero maximo de iteraciones
C			  en ASSOCDCADJUSTSUB
C
	      ELSE
	       CONTINUE
	     END IF
C
C            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^ ^^^^
C	retirar
c      dc(3)=        5.3796  !PROPIL ACETATO

c	dc(4)=      5.6326  !METILBUTANOATO


c	dc(5)= 4.7525  !metil acetato

c	dc(4)= 4.7525  !metil acetato
c	dc(6)=5.0991   !etil acetato

c	dc(3)=        5.3796  !PROPIL ACETATO

c	dc(4)=      5.6326  !METILBUTANOATO


c	dc(5)= 4.7525  !metil acetato

cc	dc(4)= 4.7525  !metil acetato
c	dc(7)=5.0991   !etil acetato

c	dc(6)= 4.7525  !metil acetato
c	dc(15)= 5.6796


	 DC(I) = WDC(1)

	     WRITE(6,*) ' FIN AJUSTE DE DC COMP NRO.= ',I

C

C
C         ------------------------------------------------------------
	  ELSE
 1325	   CONTINUE
	  END IF
********************************
*	  end if				 ! Fran. 26/05/09
********************************
	END DO
c	dc(6)=5.0991   !etil acetato

!	dc(4)=5.547
!	wdc(4)=5.547
!	dc(2)=5.87546  !decano
!	dc(3)=7.06345  !hexadecano

C
C  CALCULO DE PRESIONES DE VAPOR O DE ISOFUGACIDAD DE PUROS
C

c	dc(3)=7.497845663  !oleic acid
c	dc(11)=7.719153495  !eicosane

	NUMEXP=0  !Numero de Datos experimentales
C
	DO I=1,NC

		IF(NSAT(I).NE.0) THEN
C               ------------------------------------------------
	    	WDC(1)= DC(I)
	    	WTEMPCRIT(1)= TEMPCRIT(I)
 	    	WPRCRIT(1)  = PRCRIT(I)
	     	DO J=1,NG
			NUW(J,1) = NU(J,I)
	     	END DO
C
	     	DO K=1,NGA
			ONUASOCW(K,1) = ONUASOC(K,I)
	     	END DO
C
C               Calculo de Presiones de Vapor O ISOFUGAC. DE PUROS
C
C               ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! Bucle Nsat ! ! ! ! !
C
		DO J=1,NSAT(I)
C
	        WTEMP=TEXP(NUMEXP+1)
	        WPSAT=PEXP(NUMEXP+1)
C
C                  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
		   IF(IOBJSAT(I).EQ.50) THEN
C
		CALL ASSOCISOFPUROSSUB(WTEMPCRIT,NG,NUW,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCW,
     @ MASSOC,ENASS,
     @RKASS,WDC,WTEMP,WPSAT,WZFASL,WZFASV,APART,IERRYAF)
C
			IF(IERRYAF.EQ.0) THEN
C
			 PCALC(NUMEXP+1)=WPSAT   !AL DOPE
			 ZLCALC(NUMEXP+1)= WZFASL
			 ZVCALC(NUMEXP+1)= WZFASV
	        	 FCAUX(NUMEXP+1)= APART
			 FCION(NUMEXP+1)= W1 * FCAUX(NUMEXP+1)
			ELSE
			 PCALC(NUMEXP+1)= 0.
			 ZLCALC(NUMEXP+1)= 0.
			 ZVCALC(NUMEXP+1)= 0.
	        	 FCAUX(NUMEXP+1)= 0.
			 FCION(NUMEXP+1)= PENALIZ
			END IF
C
		    ELSE
		     CONTINUE
		   END IF
C                  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
C
C
C                  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
		   IF(IOBJSAT(I).EQ.80) THEN
C
	write(6,*)"i=",i
		CALL ASSOCPSATSUB(WTEMPCRIT,NG,NUW,QAREA,TSTR,
     @GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCW,
     @MASSOC,ENASS,
     @RKASS,WDC,WTEMP,WPSAT,WZFASL,WZFASV,IERRPSAT)
C
			IF(IERRPSAT.EQ.0) THEN
C
			 PCALC(NUMEXP+1)=WPSAT
			 ZLCALC(NUMEXP+1)= WZFASL
			 ZVCALC(NUMEXP+1)= WZFASV
 	        	 FCAUX(NUMEXP+1)= (PCALC(NUMEXP+1)
     @/PEXP(NUMEXP+1)-1)
			 FCION(NUMEXP+1)= W1 * FCAUX(NUMEXP+1)
			ELSE
C
			 PCALC(NUMEXP+1)= 0.
			 ZLCALC(NUMEXP+1)= 0.
			 ZVCALC(NUMEXP+1)= 0.
 	        	 FCAUX(NUMEXP+1)= 0.
			 FCION(NUMEXP+1)= PENALIZ
C
			END IF
C
		    ELSE
		     CONTINUE
		   END IF
C                  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
C
C
		NUMEXP = NUMEXP + 1
C
		END DO
C
C               ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! Fin Bucle Nsat ! ! ! ! !
C
C               ------------------------------------------------
C
		ELSE
		 CONTINUE
		END IF
	END DO

	!RETIRAR
!	PM(1) = 4.40E+01
!	PM(2) = 1.42285E+02   !1.16160E+02
!	PM(3) = 2.00321E+02   !2.26E+02
!	PM(4) = 2.56429E+02
!	PM(5) = 7.40794E+01
!	PM(6) = 6.00526E+01
!	PM(7) = 8.81063E+01
!	PM(8) = 1.02133E+02
!	PM(9) = 2.82467E+02


C
C COMIENZA BUCLE SOBRE LOS DATASETS > > > >>> >>  >  >  >  >  >  >  >  >  >
C
	IF(IGLE.NE.0) THEN
C
C	````````````````````````````````````````````````````
C	````````````````DATASETS````````````````````````````
	  DO KS=1,IGLE
C	    Carga de Datos
C
	    WWDC(1) = DC(IDC(KS,1))
	    WWDC(2) = DC(IDC(KS,2))
	    WWTEMPCRIT(1) = TEMPCRIT(IDC(KS,1))
	    WWTEMPCRIT(2) = TEMPCRIT(IDC(KS,2))
	    WWPRCRIT(1) = PRCRIT(IDC(KS,1))
	    WWPRCRIT(2) = PRCRIT(IDC(KS,2))
	    WWPM(1)=PM(IDC(KS,1))
	    WWPM(2)=PM(IDC(KS,2))

	    DO J=1,NG

	      NUWW(J,1) = NU(J,IDC(KS,1))
	      NUWW(J,2) = NU(J,IDC(KS,2))

	    END DO
	    DO J=1,NGA

	      ONUASOCWW(J,1) = ONUASOC(J,IDC(KS,1))
	      ONUASOCWW(J,2) = ONUASOC(J,IDC(KS,2))

	    END DO
	    If(mez(ks)==3)then

	      WWDC(3) = DC(IDC(KS,3))
	      WWTEMPCRIT(3) = TEMPCRIT(IDC(KS,3))
	      WWPRCRIT(3) = PRCRIT(IDC(KS,3))
	      DO J=1,NG

	        NUWW(J,3) = NU(J,IDC(KS,3))

	      END DO
	      DO J=1,NGA

	        ONUASOCWW(J,3) = ONUASOC(J,IDC(KS,3))

	      END DO

	    end if
C
	    IF(IHEN(KS).EQ.0) THEN   ! TPXY (IHEN=0)
c
	      NFAILURE = NOBS(KS)
	      STDVA(KS) = 0D0
	      STDVB(KS) = 0
	      DO I=1,NOBS(KS)

	        TEMP  = TEXP(NUMEXP+1)
	        PRESS = PEXP(NUMEXP+1)
	        X1 = XEXP(NUMEXP+1)
	        Y1 = YEXP(NUMEXP+1)
!	        Nota: ifl = -1 o ifl = -2
!
	        IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then
C
!	          WWDC(2) WWTEMPCRIT(2), WWPRCRIT(2)	NUWW(10,2) ONUASOCWW(12,2)
	          write(6,*) "i=",i
!	          Inicio de la variable de error. Si no se inicia puede
!	          dar cualquier cosa.
	          iErrMixIsof = 0
	          CALL ASSocISOFmezclasSUB(2,WWDC,WWTEMPCRIT,NG,NUWW,
     @               QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,
     &               ONUASOCWW,MASSOC,ENASS,RKASS,TEMP,PRESS,X1,Y1,
     &               APART1,APART2,ZLL,ZVV,IERRMIXISOF)
!
	          IF(IERRMIXISOF.EQ.0) THEN

	            ZLCALC(NUMEXP+1)= ZLL
	            ZVCALC(NUMEXP+1)= ZVV
	            FCION(NUMEXP+1) = (W4**2.)*(2.+IFL(KS))*(APART2**2) +
     @                              (APART1**2.)
	            FCION(NUMEXP+1)= DSQRT(FCION(NUMEXP+1))
	            FCAUX(NUMEXP+1)= APART1
	            FCAUX2(NUMEXP+1)= APART2

	          ELSE
C	            Si el calculo fracasa hay que descontar uno al
C	            numero de observaciones al efecto de calcular
C	            la desviacion standar
	            NFAILURE = NFAILURE - 1
	            ZLCALC(NUMEXP+1)= 0.
	            ZVCALC(NUMEXP+1)= 0.
	            FCION(NUMEXP+1)= PENALIZ
	            FCAUX(NUMEXP+1)= 0.
	            FCAUX2(NUMEXP+1)= 0.

	          END IF
	          STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)**2
	          STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)**2
C
!	        else
!	          continue
	        end if
C	        Nota: FIN ifl = -1 o ifl = -2
C
C	        Nota: ifl = 0 o ifl=1
!	        IFL(KS)=1 !RETIRAR

	        IF(IFL(KS).EQ.0.OR.IFL(KS).EQ.1) then
C
	          IF(IFL(KS).EQ.0) THEN   !  ifl=0
C
	            write(6,*) "i=",i
!	            Inicio de la variable de error. Si no se inicia puede
!	            dar cualquier cosa.
	            iErrTP = 0
	            CALL ASSocTPflashSUB(2,WWDC,WWTEMPCRIT,NG,NUWW,
     @                                 QAREA,TSTR,GSTR,GDOT,GDDOT,
     @                                 RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCWW,
     @                                 MASSOC,ENASS,RKASS,TEMP,PRESS,
     @                                 X1,Y1,ZLIQUIDO,ZVAPOR,IERRTP)
!
!	            Para flashes LL (tampoco sé qué hace acá).
!	            CALL TPflashELL(2,WWDC,WWTEMPCRIT,NG,NUWW,QAREA,TSTR,
!     @                           GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,
!     @                           ONUASOCWW,MASSOC,ENASS,RKASS,TEMP,
!     @                           PRESS,X1,Y1,ZLIQUIDO,ZVAPOR,IERRTP)

	            IF(IERRTP.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQUIDO
	              ZVCALC(NUMEXP+1) = ZVAPOR
	              XCALC(NUMEXP+1) = X1
	              YCALC(NUMEXP+1) = Y1
	              FCAUX(NUMEXP+1) = ABS(XCALC(NUMEXP+1)
     @                                                  -XEXP(NUMEXP+1))

	            ELSE
C
	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              XCALC(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF
C
!	          ELSE
!	            CONTINUE
	          END IF   ! fin ifl = 0
C
	          IF(IFL(KS).EQ.1) THEN   !  ifl=1
C
	            write(6,*) "i=",i
!	            Inicio de la variable de error. Si no se inicia puede
!	            dar cualquier cosa.
	            iErrTx = 0
	            CALL ASSocTXFLASHSUB(2,WWDC,WWTEMPCRIT,NG,NUWW,
     @                                 QAREA,TSTR,GSTR,GDOT,GDDOT,
     @                                 RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCWW,
     @                                 MASSOC,ENASS,RKASS,TEMP,PRESS,
     @                                 X1,Y1,ZLIQLIQ,ZVAPVAP,IERRTX)
C
	            IF(IERRTX.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQLIQ
	              ZVCALC(NUMEXP+1) = ZVAPVAP
	              YCALC(NUMEXP+1) = Y1
	              PCALC(NUMEXP+1) = PRESS
	              FCAUX(NUMEXP+1) = ABS(PCALC(NUMEXP+1)/
     @                                          PEXP(NUMEXP+1)-1)

	            ELSE
C
	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              PCALC(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF
C
!	          ELSE
!	            CONTINUE
	          END IF   ! fin ifl = 1
C
C                         #  #  #  #  #  #  #  #  #  #  #  #  #  #
	          IF(IFL(KS).EQ.0.AND.IERRTP.EQ.0.OR.
     @                               IFL(KS).EQ.1.AND.IERRTX.EQ.0.) THEN

!			  FCION(NUMEXP+1) = (W4**2.)*
c     @(    YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.     )**2. +
c     @IFL(KS)*(PCALC(NUMEXP+1)/PEXP(NUMEXP+1)-1.)**2.+
c     @(1-IFL(KS))*((XCALC(NUMEXP+1)*WWPM(1)/(XCALC(NUMEXP+1)*WWPM(1)+(1-
c     @XCALC(NUMEXP+1))*WWPM(2)))/(XEXP(NUMEXP+1)*WWPM(1)/(XEXP(NUMEXP+1)
c     @*WWPM(1)+(1-XEXP(NUMEXP+1))*WWPM(2)))  - 1. )**2.
	            FCION(NUMEXP+1) = (W4**2)
     @                     *(YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.D0)**2 +
     @             IFL(KS)*(PCALC(NUMEXP+1)/PEXP(NUMEXP+1) - 1.D0)**2 +
     @        (1-IFL(KS))*(XCALC(NUMEXP+1)/(XEXP(NUMEXP+1))  - 1.D0)**2

c			  FCION(NUMEXP+1) = (W4**2.)*
c     @(    YCALC(NUMEXP+1)-YEXP(NUMEXP+1)   )**2. +
c     @ IFL(KS) * (  PCALC(NUMEXP+1)-PEXP(NUMEXP+1)    )**2. +
c     @ (1-IFL(KS)) * (XCALC(NUMEXP+1)-(XEXP(NUMEXP+1))  )**2.

	            FCION(NUMEXP+1)=   DSQRT(FCION(NUMEXP+1))
	            FCAUX2(NUMEXP+1) =ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	          ELSE

	            FCION(NUMEXP+1)= PENALIZ
	            FCAUX2(NUMEXP+1) = 0.D0

	          END IF
C                         #  #  #  #  #  #  #  #  #  #  #  #  #  #
	          STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)
		    STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)
C
!	        else
!	          continue
	        end if
C
C	-	-	 Nota: FIN ifl = 0 o ifl=1
C
C
			 NUMEXP = NUMEXP + 1
C
			END DO
C                       desviaciones standar
C
			IF(NFAILURE.NE.0) THEN
			  STDVA(KS) = STDVA(KS)/NFAILURE
			  STDVB(KS) = STDVB(KS)/NFAILURE
			 ELSE
			  STDVA(KS) = 0.
			  STDVB(KS) = 0.
			END IF
C
C                       ddddddddddddddddddddddddddddddddd
		 ELSEIF (ihen(ks).eq.5) then

		    nfailure=nobs(ks)
		    stdva(ks) = 0
			stdvb(ks) = 0
			stdvc(ks) = 0
			WRITE(6,*)"DADOS=", KS
			do i=1,nobs(ks)
			 temp  = texp(numexp+1)
			 press = pexp(numexp+1)
			 x1 = xexp(numexp+1)
			 x2 = xexp2(numexp+1)
	!		 rgam1 = rgamexp1(numexp+1)
	!		 rgam2 = rgamexp2(numexp+1)
	!		 rgam3 = rgamexp3(numexp+1)

			 WRITE(6,*)"OBS=", NOBS(KS)
			iaux=mez(ks)

!	      Iniciación de iErrGam. Sino, el compilador puede hacer lo que se le cante.
	      iErrGam = 0

		call rgamasub(iaux,wwdc, wwtempcrit,ng,nuww,qarea,tstr, gstr,
     *gdot,gddot,rkij, rktij,alfaij,nga, onuasocww, massoc,enass,rkass,
     *temp, press, x1,x2,rgam1,rgam2,rgam3,zli,ierrgam)

       			IF(ierrgam==0) THEN
			  	zlicalc(numexp+1) = zli
			  	rgamcalc1(numexp+1) = rgam1
				rgamcalc2(numexp+1) = rgam2
		    fcaux(numexp+1)=rgamcalc1(numexp+1)/rgamexp1(numexp+1)-1.
			fcaux2(numexp+1)=rgamcalc2(numexp+1)/rgamexp2(numexp+1)-1.
	fcion(numexp+1)=(w6)**2.*(rgamcalc1(numexp+1)/
     @rgamexp1(numexp+1)-1.)**2.+
     @w7**2*(rgamcalc2(numexp+1)/rgamexp2(numexp+1)-1.)**2.
		If(mez(ks)==3)then
			rgamcalc3(numexp+1) = rgam3
			fcaux3(numexp+1)=rgamcalc3(numexp+1)/rgamexp3(numexp+1)-1.
			fcion(numexp+1)=fcion(numexp+1)+
     @		w8*(rgamcalc3(numexp+1)/rgamexp3(numexp+1)-1.)**2.
		end if
				fcion(numexp+1)=sqrt(fcion(numexp+1))

				else
C
			  	nfailure = nfailure - 1
C
  				zlicalc(numexp+1) = 0.
			    rgamcalc1(numexp+1) = 0.
				rgamcalc2(numexp+1) = 0.
				rgamcalc3(numexp+1) = 0.
	           	fcaux(numexp+1) = 0.
		       	fcaux2(numexp+1) = 0.
				fcaux3(numexp+1) = 0.
				fcion(numexp+1)= penaliz
				end if

				STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)**2
			    STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)**2
				STDVC(KS) = STDVC(KS) + FCAUX3(NUMEXP+1)**2

			    numexp=numexp+1

			 end do

			  IF(NFAILURE.NE.0) THEN
			  STDVA(KS) = DSQRT(STDVA(KS)/NFAILURE)
			  STDVB(KS) = DSQRT(STDVB(KS)/NFAILURE)
			  STDVC(KS) = DSQRT(STDVC(KS)/NFAILURE)
			  ELSE
			  STDVA(KS) = 0.
			  STDVB(KS) = 0.
			  STDVC(KS) = 0.
			  END IF

		 ELSE
		CONTINUE
		END IF	! TPXY (IHEN=0)
C

		END DO
C		````````````````DATASETS````````````````````````````
C		````````````````````````````````````````````````````
C
	 ELSE
	  CONTINUE
	END IF
C
C TERMINA BUCLE SOBRE LOS DATASETS > > > >>> >>  >  >  >  >  >  >  >  >  >
C
	SUMQ = 0
	DO I=1,NCXP
		SUMQ = SUMQ + FCION(I)**2
	END DO
C
C
	OPEN(UNIT=1,TYPE='OLD',file='paraest.out',FORM='FORMATTED',
     @ ACCESS='APPEND')
C
	WRITE(1,*) (XX(I),I=1,NPAR),SUMQ
c	WRITE(1,*) ' DC, I=1,NC ',(DC(I),I=1,NC)
c	WRITE(1,*) ' tc, I=1,NC ',(tempcrit(I),I=1,NC)
C
	CLOSE(UNIT=1)
C
	RETURN
	END !FIN FUNC
C//////////////////////////////////////////////////////////////////////////////
	SUBROUTINE INICDC(NC,PRCRIT,TEMPCRIT,DC)
C
C	Estima DC solo con informacion de los puntos criticos de los compuestos
C	puros.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION DC(20)
	DIMENSION  TEMPCRIT(20)
	DIMENSION  PRCRIT(20)
C
	R=82.05
C
C R EN ATM CM3/(K MOL)
C
	DO I=1,NC
C
******************************
*		if(DC(I).EQ.0)then ! Fran. 26/05/09
******************************
*
			DC(I) = 0.4471866767*((R*TEMPCRIT(I)/PRCRIT(I)))**(1./3.)
*
******************************
*		end if              ! Fran. 26/05/09
******************************
C
	END DO
C
	END !FIN INICDC
C//////////////////////////////////////////////////////////////////////////
	SUBROUTINE FORSCAL(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,
     @numtc,NUMDC)
C
C RETORNA IDXX
C
	IMPLICIT REAL*8 (A-H,O-Z)
C
	DIMENSION IDG(8,2) ! identificacion (indices) de param gii o kij a
	DIMENSION ITEMP(8)
	DIMENSION IDXX(8)
C
	NPAR = 0
C - - - - - - - - - - - - - - - - - - - - - - - - - -1
	IF(NNP.NE.0) THEN
C
		DO I=1,NNP
C
		    IF(IDG(I,1).EQ.IDG(I,2)) THEN
			IDXX(NPAR+1)=1
			NPAR = NPAR + 1
		     ELSE
			IDXX(NPAR+1)=2
			NPAR = NPAR + 1
		    END IF
C
		    IF(ITEMP(I).EQ.1) IDXX(NPAR+1) = 3
		    IF(ITEMP(I).EQ.1) NPAR = NPAR + 1
C
		    IF(ITEMP(I).EQ.2) IDXX(NPAR+1) = 3
		    IF(ITEMP(I).EQ.2) NPAR = NPAR + 1
C
		END DO
C
	 ELSE
	  CONTINUE
	END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - -1
C
	IF(NNQ.NE.0) THEN
		DO I=1, NNQ
			IDXX(NPAR+1)= 4
			NPAR = NPAR + 1
		END DO
	 ELSE
	  CONTINUE
	END IF
C
	IF(NALFA.NE.0) THEN
		DO I=1, NALFA
			IDXX(NPAR+1)= 5
			NPAR = NPAR + 1
		END DO
	 ELSE
	  CONTINUE
	END IF
C
	IF(NENASS.NE.0) THEN
		DO I=1,NENASS
			IDXX(NPAR+1)= 6
			NPAR = NPAR + 1
		END DO
	 ELSE
	  CONTINUE
	END IF
C
	IF(NRKASS.NE.0) THEN
		DO I=1,NRKASS
			IDXX(NPAR+1)= 7
			NPAR = NPAR + 1
		END DO
	 ELSE
	  CONTINUE
	END IF
C
C
	IF(numtc.NE.0) THEN
		DO I=1,numtc
			IDXX(NPAR+1)= 8
			NPAR = NPAR + 1
		END DO
	 ELSE
	  CONTINUE
	END IF

	IF(NUMDC.NE.0) THEN
		DO I=1,NUMDC
			IDXX(NPAR+1)= 9
			NPAR = NPAR + 1
		END DO
	 ELSE
	  CONTINUE
	END IF
C
	END !FIN FORSCAL
C////////////////////////////////////////////////////////////////////
	SUBROUTINE XXSCAL(NPAR,XX,IDXX)
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XX(8), IDXX(8)
C
	COMMON/ESCALADO/FGGG,FRKIJ,FTEMP,FQAREA,
     @FALFA,FENASS,FRKASS,ftc,FDC
C
C TRANSFORMA LOS VALORES DE XX SEGUN LOS SIGUIENTES FACTORES DE ESCALADO
C
	FGGG  = 1.D-14
	FRKIJ = 1.D-8
	FTEMP = 1.D-7
	FQAREA = 1.D-5
	FALFA = 1.D-10
	FENASS = 1.D-12    !TENTATIVO
	FRKASS = 1.D-7     !TENTATIVO
	ftc = 1.d-8	       !TENTATIVO
	FDC = 1.D-6        !TENTATIVO

c	FGGG  = 1.D-10
c	FRKIJ = 1.D-5
c	FTEMP = 1.D-6
c	FQAREA = 1.D-5
c	FALFA = 1.D-6
c	FENASS = 1.D-9     !TENTATIVO
c	FRKASS = 1.D-4     !TENTATIVO
c	FDC = 1.D-6        !TENTATIVO

C
	DO I=1,NPAR
		IF(IDXX(I).EQ.1) XX(I) = FGGG   * XX(I)
		IF(IDXX(I).EQ.2) XX(I) = FRKIJ  * XX(I)
		IF(IDXX(I).EQ.3) XX(I) = FTEMP  * XX(I)
		IF(IDXX(I).EQ.4) XX(I) = FQAREA * XX(I)
		IF(IDXX(I).EQ.5) XX(I) = FALFA  * XX(I)
		IF(IDXX(I).EQ.6) XX(I) = FENASS * XX(I)
		IF(IDXX(I).EQ.7) XX(I) = FRKASS * XX(I)
		IF(IDXX(I).EQ.8) XX(I) = ftc    * XX(I)
		IF(IDXX(I).EQ.9) XX(I) = FDC    * XX(I)
	END DO
C
	END !FIN XXSCAL
C//////////////////////////////////////////////////////////////////
	SUBROUTINE XXDESSCAL(NPAR,XX,IDXX,FGGG,FRKIJ,FTEMP,FQAREA,
     @FALFA,FENASS,FRKASS,XXAUX,ftc,FDC)
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XX(8), IDXX(8), XXAUX(8)
C
C	ANTITRANSFORMA LOS VALORES DE XX SEGUN LOS FACTORES DE ESCALADO
C  ESTABLECIDOS EN LA RUTINA XXSCAL
C
C
	DO I=1,NPAR
		XXAUX(I) = XX(I)
		IF(IDXX(I).EQ.1) XXAUX(I) = XXAUX(I)/FGGG
		IF(IDXX(I).EQ.2) XXAUX(I) = XXAUX(I)/FRKIJ
		IF(IDXX(I).EQ.3) XXAUX(I) = XXAUX(I)/FTEMP
		IF(IDXX(I).EQ.4) XXAUX(I) = XXAUX(I)/FQAREA
		IF(IDXX(I).EQ.5) XXAUX(I) = XXAUX(I)/FALFA
		IF(IDXX(I).EQ.6) XXAUX(I) = XXAUX(I)/FENASS
		IF(IDXX(I).EQ.7) XXAUX(I) = XXAUX(I)/FRKASS
		IF(IDXX(I).EQ.8) XXAUX(I) = XXAUX(I)/ftc
		IF(IDXX(I).EQ.9) XXAUX(I) = XXAUX(I)/FDC
	END DO
C
	END !FIN XXDESSCAL
C//////////////////////////////////////////////////////////////////
	SUBROUTINE WRITESTAT(SSQ,XJAC,WORK,INFER,IER)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XJAC(260,8), WORK(436)
C
	OPEN(UNIT=1,TYPE='OLD',file='paraest.out',FORM='FORMATTED',
     @ ACCESS='APPEND')
C
	WRITE(1,95)
 95	FORMAT(/,X,' VALOR FINAL FCION OBJETIVO = SSQ = ',/)
	WRITE(1,*) SSQ
C
	WRITE(1,97) IER,INFER
 97	FORMAT(/,' (VER MANUAL DE LA RUTINA ZXSSQ (Lib IMSL) ',/,
     @X,' IER = ',I5,/,X,' INFER = ',I5)
C
	WRITE(1,99) WORK(1),WORK(2),WORK(3),WORK(5)
 99	FORMAT(/,'  NORMA DEL GRADIENTE = ',D12.5,/,
     @X,' NUMERO DE EVALUACIONES DE LA FC OBJ. = ',F10.2,/
     @X,' NUMERO DE DIGITOS SIGNIFICATIVOS EN LOS PARAM. ESCALADOS = '
     @,F10.2,/,X,' NUMERO DE ITERACIONES = ',F10.2)
	CLOSE(UNIT=1)
C
	END !FIN WRITESTAT
C//////////////////////////////////////////////////////////////////
	SUBROUTINE WRITESAL2(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,
     @  RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,ISAT,
     @  IGLE,IVE,IUNI,
     @  TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,
     @  NTEXT,
     @TEXP,PEXP,XEXP,xexp2, YEXP,IDSAT,DC,tempcrit,XCALC,YCALC,ZLCALC,
     @  ZVCALC,FCAUX,PCALC,IOBJSAT,FCAUX2,fcaux3,STDVA,STDVB,stdvc,
     @ rgamcalc1, rgamcalc2,rgamcalc3,numtc,idtc,NUMDC,IDDC,TSPLINE,
     @rgamexp1, rgamexp2,rgamexp3,zlicalc)

c
C ESCRIBE DATOS ARCHIVO DE SALIDA paraest.out
C
	IMPLICIT REAL*8 (A-H,O-Z)
C
C  DIM DISPERSION
C
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
	DIMENSION TSPLINE(10)
C
C  DIM ASOCIACION
C
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
        DIMENSION ONUASOC(12,20)
C
C  DIM DIAMETRO DE ESFERAS RIGIDAS EN EL PUNTO CRITICO
C
	DIMENSION DC(20)
C
C
	DIMENSION IDSAT(20)
C
	DIMENSION  TSAT(20), PSAT(20),tempcrit(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260), xexp2(260)
	dimension rgamexp1(260),rgamexp2(260),rgamexp3(260)
C
C DIMENSIONAMIENTO  NRO. DE DATOS DE PVAP COMP PUROS
C
	DIMENSION NSAT(20)
	DIMENSION IOBJSAT(20)
C
C DIMENSIONAMIENTOS QUE TIENEN QUE VER CON LOS DATASETS
C
	DIMENSION IDC(25,3),NOBS(25),IHEN(25),mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
C
	DIMENSION PCALC(260), XCALC(260), YCALC(260)
	DIMENSION  ZLCALC(260), ZVCALC(260)
	dimension rgamcalc1(260),rgamcalc2(260),rgamcalc3(260)
 	dimension zlicalc(260)

C
	DIMENSION FCAUX(260)
	DIMENSION FCAUX2(260)
	DIMENSION FCAUX3(260)

C
C DIMENSIONAMIENTO DE LA DESVIACION STANDAR PARA LAS COMP PUROS
C
	DIMENSION STDVSAT(20)
C
C DIMENSIONAMIENTO DESVIACIONES STANDAR PARA MEZCLAS
C
	DIMENSION STDVA(25)
	DIMENSION STDVB(25)
	DIMENSION STDVC(25)

C
C DIMENSIONAMIENTO ARREGLO DE IDENTIFICACION DE PARAM. DC CONSIDERADOS
C   COMO SIMPLES VARIABLES DE OPTIMIZCION
C
	DIMENSION IDDC(8)

	dimension idtc(8)
C

	OPEN(UNIT=1,TYPE='OLD',file='paraest.out',FORM='FORMATTED',
     @ ACCESS='APPEND')
c	OPEN(UNIT=6,NAME='DESVIOS.TXT',FORM='FORMATTED')

C
C
	WRITE(1,47)
 47	FORMAT(/,X,'--------------------------------------------------',
     @/,X,'------------ PARAMETROS FINALES ----------------',
     @/,X,'--------------------------------------------------',
     @/,X,' PARAMETROS DE GRUPOS PUROS (termino dispersivo)',/)
C
	WRITE(1,3) (I, I=1,NC)
 3	FORMAT(40X,'COMPONENTES',/,X,
     @' Gr.Nro.:  T*        Q        TSPLINE ',<NC>I3,/)
C
	DO I=1,NG
		WRITE(1,4) I,TSTR(I),QAREA(I),TSPLINE(I), (NU(I,J), J=1,NC)
 4		FORMAT(X,I3,5X,2F10.4,X,E10.4,<NC>I3)
	END DO
C
C
C NU (I,J)= NUMERO DE GRUPOS DEL TIPO I EN EL COMPUESTO DE TIPO J
C
	WRITE(1,5)
 5	FORMAT(/,X,' ENERGIAS DE INTERACCION gjj ',
     @/,X,'(la forma funcional cambia a altas temperaturas)',/)
C
	DO I=1,NG
C
	WRITE(1,6) I,I,GSTR(I), GDOT(I), GDDOT(I)
 6	FORMAT(X,' g(',I3,I3,') = ',
     @F10.1,' * [ 1 + ',F10.4,' *(T/T* - 1) + ',
     @F10.4, ' *LN(T/T*)] ',/)
C
	END DO
C
	WRITE(1,7)
 7	FORMAT(/,X,' PARAMETROS BINARIOS k*ij ',/)
C
C
	DO I=1,NG
		WRITE(1,8) (RKIJ(I,J), J=1,NG)
 8		FORMAT(X,<NG>F10.4)
	END DO
C
	WRITE(1,9)
 9	FORMAT(/,X,' DEPENDENCIA RESP. A LA TEMP. : kprimaij ',/)
C
C
	DO I=1,NG
		WRITE(1,8) (RKTIJ(I,J), J=1,NG)
	END DO
C
	WRITE(1,10)
 10	FORMAT(/,X,' PARAMETROS BINARIOS DE NO AZAR alfaij ',/)
C
C
	DO I=1,NG
		WRITE(1,8) (ALFAIJ(I,J), J=1,NG)
	END DO
C
C
	WRITE(1,186) NGA
 186	FORMAT(/,X,' Numero de grupos considerados en la contribucion
     @de asociacion = ',I3)
C
C
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IF(NGA.GT.0) THEN
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
	WRITE(1,11) (I, I=1,NC)
 11	FORMAT(/,X,' PARAMETROS DE ASOCIACION ',/,23X,'COMPONENTES',/,
     @' GRUPO 	  # SIT ASOC  ',<NC>I3,/)
C
	DO I=1,NGA
	WRITE(1,12) I, MASSOC(I), (ONUASOC(I,J), J=1,NC)
 12	FORMAT(3X,I3,9X,I3,6X,<NC>f5.2)
	END DO
C
	WRITE(1,13)
 13 	FORMAT(/,X,'PARAMETROS DE ENERGIA DE ASOCIACION (Kelvin)  ',/)
C
C
C ESCRITURA DE PARAMETROS ENERGETICOS DE ASOCIACION
C
	DO J=1,NGA
	      IF(MASSOC(J).EQ.0) GO TO 20
		DO L=1,NGA
                  IF(MASSOC(L).EQ.0) GO TO 101
			DO I=1,MASSOC(J)
				DO K=1,MASSOC(L)
		                  WRITE(1,14) I,J,K,L,ENASS(I,J,K,L)
 14	FORMAT(X,' ENASS( ',I3,I3,I3,I3,' )  = ',F10.4)
				END DO
			END DO
 101	        CONTINUE
	        END DO
 20	      CONTINUE
	END DO
C
C
C ESCRITURA DE PARAMETROS VOLUMETRICOS DE ASOCIACION (cm3/mol)
C
	WRITE(1,15)
 15	FORMAT(/,X,'PARAMETROS DE VOLUMEN DE ASOCIACION (cm3/mol) ',/)
C
	DO J=1,NGA
	      IF(MASSOC(J).EQ.0) GO TO 30
		DO L=1,NGA
                  IF(MASSOC(L).EQ.0) GO TO 50
			DO I=1,MASSOC(J)
				DO K=1,MASSOC(L)
				    WRITE(1,16) I,J,K,L,RKASS(I,J,K,L)
 16	FORMAT(X,' RKASS( ',I3,I3,I3,I3,' ) = ',F16.10)
				END DO
			END DO
 50	        CONTINUE
	        END DO
 30	      CONTINUE
	END DO
C
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 ELSE
	  CONTINUE
	END IF
C      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
C
C
C ESCRITURA DEL RESULTADO DEL AJUSTE DE DC
C
	DO I=1,NC
C
	    DO JAMO=1,NUMDC
		IF(IDDC(JAMO).EQ.I)  THEN
C
		  WRITE(1,1548) I,DC(I)
 1548		  FORMAT(/,X,' El diametro de esferas rigidas en el punto',
     @/,X,' critico del componente Nro.',I3,/' fue considerado otra
     @variable de optimizacion ',/,X,' DC = ', F10.7)
C
		  GO TO 505
C
		 ELSE
		  CONTINUE
		END IF
	    END DO
C
		IF(TSAT(I).NE.0) WRITE(1,500) I,TSAT(I),PSAT(I),DC(I)
 500		FORMAT(/,X,' El diametro de esferas rigidas en el punto',
     @/,X,' critico del componente Nro.',I3,' fue reajustado para ',
     @/,x,' reproducir el punto de la curva de presion de vapor ',
     @/,x,' [ ',F10.4,' K ;',F10.4,' atm ] ',/,X,' DC = ', F10.7)
 505	END DO

 	DO I=1,NC
C
	    DO mm=1,numtc
		IF(idtc(mm).EQ.I)  THEN
C
		  WRITE(1,1549) I,tempcrit(I)
 1549	  FORMAT(/,X,'La temperatura critica del componente Nro.',I3,/,X,
     @ 'fue considerado otra variable de optimizacion ',/,'tc= ', F8.2)
C
C
		 ELSE
		  CONTINUE
		 END IF
		 END DO
	END DO

C

C
C ESCRITURA DEL RESULTADO DEL AJUSTE DE DATOS DE COMPUESTOS PUROS
C
C
	NCXP = 0 !Cuenta el numero de datos de equilibrio
C
 2	FORMAT(20I3)
C
C       ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
	IF(ISAT.NE.0) THEN
		      WRITE(1,338)
 338	FORMAT(/,X,'***** RESULTADOS: COMPUESTOS PUROS *****',/)
C                     +++++++++++++++++++++++++++++++++++++++++++++++
			DO I=1,NC
C                        ---------------------------------------------
			 IF(NSAT(I).NE.0) THEN
			   IF(IOBJSAT(I).EQ.50) THEN
		            WRITE(1,393)
 393			    FORMAT(/,X,'  COMP   TEXP/K     PEXP/atm
     @         PEXP/atm          ZLIQ       ZVAP
     @  ERR,REL,%,ISOF',/)

			    ELSE
			     CONTINUE
			   END IF
C
C                          '  '  '  '  '  '   '   '  '  '  '  ' '  '  '  '
C
			   IF(IOBJSAT(I).EQ.80) THEN
		            WRITE(1,333)
 333			    FORMAT(/,X,'  COMP   TEXP/K     PEXP/atm
     @        PCALC/atm        ZLIQ        ZVAP
     @  ERR,REL,%,P ',/)
			    ELSE
			     CONTINUE
			   END IF
C
C                          '  '  '  '  '  '   '   '  '  '  '  ' '  '  '  '
C
C
C                         de paso calculo desviacion standar para cada
C                          componente
			   NFOLI=NSAT(I)
			   STDVSAT(I)=0
C
C                          ...................................
			   DO J=1,NSAT(I)
C
 	WRITE(1,334) I,TEXP(NCXP+1),PEXP(NCXP+1),
     @PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100
 334	FORMAT(X,I3,4X,F10.4,X,F16.10,X,F16.10,X,F11.8,X,F10.7,X,F10.4)
C
			    STDVSAT(I) = STDVSAT(I) + FCAUX(NCXP+1)**2
			    if(FCAUX(NCXP+1).eq.0) NFOLI = NFOLI - 1
C
			    NCXP = NCXP + 1
C
			   END DO
C                          ...................................
			   IF(NFOLI.NE.0) THEN
			     STDVSAT(I) = STDVSAT(I) / NFOLI
			    ELSE
			     STDVSAT(I) = 0.
			   END IF
C
			   STDVSAT(I) = 100*DSQRT(STDVSAT(I))
		           WRITE(1,335) STDVSAT(I)
					WRITE(6,*) STDVSAT(I)
 335	FORMAT(/,X,' Desviacion Standar %',41X,f10.4,/)

C
			  ELSE
			   CONTINUE
			 END IF
C                        ---------------------------------------------
			END DO
C                     +++++++++++++++++++++++++++++++++++++++++++++++
	 ELSE
	  CONTINUE
	END IF
C       ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
C
C ESCRITURA DEL RESULTADO DEL AJUSTE DE DATOS SOBRE MEZCLAS BINARIAS
C
	IF(IGLE.NE.0) THEN
		WRITE(1,340)
 340		FORMAT(/,X,'********** RESULTADOS: MEZCLAS BINARIAS ***
     @****************',/)
C
C		````````````````````````````````````````````````````
		DO KS=1,IGLE
       		WRITE(1,339) (NTEXT(KS,J),J=1,65)
 339 		FORMAT(/,40A2)
c
C
		IF(IHEN(KS).EQ.0) THEN   !ZIMBAWE
C
			IF(IFL(KS).EQ.-2) WRITE(1,343) IDC(KS,2)
 343			FORMAT(X,/,' NOTA: El apartamiento con respecto
     @ a la isofugacidad del componente ',I3,/,X,
     @' no fue utilizado en la evaluacion de la funcion objetivo,
     @aunque',/,X,'
     @ igualmente se informa ',/)
C
C                       \  \  \  \  \  \  \  \ \ \  \  \  \  \  \  \
			IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then
			WRITE(1,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),
     @IDC(KS,2)
 341			FORMAT(X,'  TEMP/K    PRES/ATM
     @  XEXP(',I2,')  YEXP(',I2,') FUG(',I2,')%
     @ FUG(',I2,')%   ZLIQ      ZVAP ',/)
C                       ddddddddddddddddddddddddddddddddd
			DO I=1,NOBS(KS)
C
 			  WRITE(1,345) TEXP(NCXP+1),PEXP(NCXP+1),
     @			XEXP(NCXP+1),
     @			YEXP(NCXP+1),FCAUX(NCXP+1),
     @			FCAUX2(NCXP+1),
     @			ZLCALC(NCXP+1), ZVCALC(NCXP+1)
 345	FORMAT(6F10.5,2F9.6)
C
			    NCXP = NCXP + 1
C
C
			END DO
C                       ddddddddddddddddddddddddddddddddd

			WRITE(1,347) STDVA(KS),STDVB(KS)
			WRITE(6,*) STDVA(KS),STDVB(KS)

 347			FORMAT(/,X,' Desv.Standar %',24X,2f10.4,/)
C
			 ELSE
			  CONTINUE
			END IF
C                       \  \  \  \  \  \  \  \ \ \  \  \  \  \  \  \
C                       \  \  \  \  \  \  \  \ \ \  \  \  \  \  \  \
			IF(IFL(KS).EQ.0) then
			WRITE(1,351) IDC(KS,1),IDC(KS,1),IDC(KS,1),
     @IDC(KS,1)
 351			FORMAT(X,'  TEMP/K   PRES/ATM  XEXP(',I2,')  XCLC(',I2,')
     @ ERR%    YEXP(',I2,') YCALC(',I2,')
     @ERR%',/)

C                       ddddddddddddddddddddddddddddddddd
			DO I=1,NOBS(KS)
C
 			  WRITE(1,349) TEXP(NCXP+1),PEXP(NCXP+1),
     @			XEXP(NCXP+1),XCALC(NCXP+1), FCAUX(NCXP+1),
     @			YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1),
     @			ZLCALC(NCXP+1), ZVCALC(NCXP+1)
 349	FORMAT(8F10.5,2X,' ZLIQ = ',F10.8,'
     @       ZVAPOR =  ',
     @F10.8)
C
			    NCXP = NCXP + 1
C
C
			END DO
C                       ddddddddddddddddddddddddddddddddd

			WRITE(1,353) STDVA(KS),STDVB(KS)
			WRITE(6,*) STDVA(KS),STDVB(KS)

 353			FORMAT(/,X,' Desv.Standar %',24X,
     @			f10.4,20X,f10.4,/)
C
			 ELSE
			  CONTINUE
			END IF
C                       \  \  \  \  \  \  \  \ \ \  \  \  \  \  \  \
C                       \  \  \  \  \  \  \  \ \ \  \  \  \  \  \  \
			IF(IFL(KS).EQ.1) then
			WRITE(1,355) IDC(KS,1),IDC(KS,1),IDC(KS,1)
C
 355			FORMAT(X,'  TEMP/K   XEXP(',I2,') PRES/ATM  PREScalc  ERR%
     @YEXP(',I2,')  YCALC(',I2,') @ERR% ',/)


C                       ddddddddddddddddddddddddddddddddd
			DO I=1,NOBS(KS)
C
 				  WRITE(1,357) TEXP(NCXP+1),XEXP(NCXP+1),
     @ 			PEXP(NCXP+1),PCALC(NCXP+1), FCAUX(NCXP+1),
     @			YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1),
     @			ZLCALC(NCXP+1), ZVCALC(NCXP+1)
 357	  FORMAT(F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X
     @,F10.5,2X,F10.5,2X,'  ZLIQ = ',F10.8,'
     @      ZVAP = ',F10.8)



C
			    NCXP = NCXP + 1
C
C
			END DO
C                       ddddddddddddddddddddddddddddddddd

			WRITE(1,359) STDVA(KS),STDVB(KS)
			WRITE(6,*) STDVA(KS),STDVB(KS)

 359			FORMAT(/,X,' Desv.Standar %',24X,
     @			f10.4,20X,f10.4,/)
C
			 ELSE
			  CONTINUE
			END IF
C                       \  \  \  \  \  \  \  \ \ \  \  \  \  \  \  \
C
		 ELSE

CTIZ       aqui se agrega escritura de resultados si IHEN=5 (coef. de actividad)

		   CONTINUE
	     END IF
C
		  IF (IHEN(KS).EQ.5)THEN

				If(mez(ks)==2)then
	write(1,1344)IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2),
     @						 IDC(KS,2)

			  	DO I=1,NOBS(KS)

				WRITE(1,1345)TEXP(NCXP+1),XEXP(NCXP+1),
     @rgamexp1(ncxp+1),rgamcalc1(NCXP+1),rgamexp2(ncxp+1),
     @rgamcalc2(NCXP+1), zlicalc(NCXP+1)


 1344	FORMAT(X,'TEMP/K       XEX(',I2,')   Gex(',I2,')   Gc(',I2,')
     @Gex(',I2,')   Gc(',I2,')'       ,'zliq')

 1345			format(X, 2f10.2, 4f10.3,f10.3)


 1347			FORMAT(2F8.4,f7.1,2F8.4,f6.3,f8.4)
 				     NCXP = NCXP + 1
				END DO
				WRITE(1,1379) STDVA(KS),STDVB(KS)
				WRITE(6,*) STDVA(KS),STDVB(KS)

 1379			FORMAT(/,X,' Desv.Standar %',2X,
     @			f10.4,7X,f10.4,/)


				elseif(mez(ks)==3)then
		write(1,1350) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,1),
     @				 IDC(KS,2),IDC(KS,2),IDC(KS,3),IDC(KS,3)

				DO I=1,NOBS(KS)
			  	WRITE(1,1348)TEXP(NCXP+1),XEXP(NCXP+1),xexp2(ncxp+1),
     @rgamexp1(NCXP+1),rgamcalc1(NCXP+1),rgamexp2(NCXP+1),
     @rgamcalc2(NCXP+1),rgamexp3(NCXP+1),rgamcalc3(NCXP+1)

 1348					FORMAT(X, f6.2 , 2(1X,e10.4), 6(1X,E12.3))
 1349					FORMAT(3F9.6)
 1350					FORMAT(X,'TEMP/K XEX(',I2,') XEX(',I2,') Gex(',I
     @2,') Gc(',I2,') Gex(',I2,') Gc(',I2,'),Gex(',I2,') Gc(',I2,')')

                    NCXP = NCXP + 1
				END DO
			WRITE(1,1359) STDVA(KS)*100.,STDVB(KS)*100.,STDVC(KS)*100.
			WRITE(6,*) STDVA(KS)*100.,STDVB(KS)*100.,STDVC(KS)*100.
 1359			FORMAT(/,X,' Desv.Standar %',11X,
     @			f10.4,6X,f10.4,7X,f10.4,/)

				end if


	      ELSE
            CONTINUE
	      END IF
CTIZ       fin de agregado
CTIZ		END IF	!ZIMBAWE
C

		END DO
C		````````````````````````````````````````````````````
C
	 ELSE
	  CONTINUE
	END IF
C
	CLOSE(UNIT=1)
c	CLOSE(UNIT=6)
	CLOSE(UNIT=4)


C
	END !FIN WRITESAL2
C///////////////////////////////////////////////////////////
	SUBROUTINE FORDC(NUMDC,IDDC,DC,NPAR,XX)
C
C CARGA EN LOS LUGARES ADECUADOS DE XX LOS VALORES INICIALES DE
C   LOS DC QUE SE USARAN
C  COMO CUALQUIER OTRO PARAMETRO AJUSTABLE
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION IDDC(8)
	DIMENSION DC(20)
	DIMENSION XX(8)
C
	IF(NUMDC.EQ.0) GO TO 1234
C
	DO I= (NPAR-NUMDC+1),NPAR
		    XX(I)=DC(IDDC(I-npar+numdc))
	END DO
C
 1234	END !FIN FORDC
