PROGRAM ASSOCPARAEST
!
!
!	Ajuste de parametros con modelo GCA-EoS: La estructura molecular
!	puede ser definida de distinta manera según se trate de cuantificar
!	la contribucion dispersiva o de asociacion.
!
!
!	Admite a dc como un parametro ajustable frente a datos de equilibrio
!
!	Las opciones para dc son las siguientes
!
!	- Estimarlo en base al punto critico (ecuacion explicita, no hay
!	  ajuste de dc)
!
!	- Ajustarlo frente a un dato determinado de la curva de presion
!	  de vapor del componente en cuestion.
!
!	- Considerarlo como cualquier otro parametro ajustable. esto significa
!	  que podra ser ajustado frente a todos los datos que se alimenten
!	  en que tenga participacion el componente en cuestion, sean
!	  estos binarios y/o de componente puro.
!
!
!	Con penalizacion.
!
!
!	La lectura de datos puede hacerse por pantalla o desde un archivo de
!	datos
!
!	Contiene modulo de escalado de los parametros a ajustar
!
	IMPLICIT REAL*8 (A-H,O-Z)
!	Implicit none
	parameter(NDM=260,NCM=20 )
!
!	La declaracion que viene la exige la rutina ZXSSQ de la IMSL
	EXTERNAL FUNC
!
!	Dim dispersion
	INTEGER NU(10,20)
	DOUBLE PRECISION QAREA(10), TSTR(10), GSTR(10)
	DOUBLE PRECISION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DOUBLE PRECISION ALFAIJ(10,10)
	DOUBLE PRECISION TSPLINE(10)
!
!	Dim asociacion
	INTEGER MASSOC(12)
      DOUBLE PRECISION ENASS(6,12,6,12), RKASS(6,12,6,12)
      DOUBLE PRECISION ONUASOC(12,20)
!
!	Doble precisionamiento diametro de esferas rigidas en el pto critico.
	DOUBLE PRECISION DC(20)
!
!	DIM ajuste de parametros?
	DOUBLE PRECISION PARM(4)
!
	INTEGER IDSAT(20)
!
	DOUBLE PRECISION TEMPCRIT(20)
	DOUBLE PRECISION PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(NDM), PEXP(NDM), XEXP(NDM), YEXP(NDM),xExp2(NDM),yExp2(NDM),rhoExp(NDM)
	dimension rgamexp1(NDM),rgamexp2(NDM),rgamexp3(NDM)
!
!
!
!	Double precisionamiento  nro. de datos de pvap comp puros ¿por qué es entero el número de datos?
	INTEGER NSAT(20)
!
!	Funcion objetivo para el comp. puro
	INTEGER IOBJSAT(20)
!
!	Doble precision que tiene que ver con los sets de datos:
	INTEGER IDC(25,3), NOBS(25), IHEN(25), mez(25)
	INTEGER NTEXT(25,65)
	INTEGER IFL(25)
!
!		ID 	==> identificacion o identicos
!		IDID	==> identificacion de identicos
!
!	Los parametros se identifican a traves de su naturaleza y de sus indices
	INTEGER ITEMP(8)
!
!	Identificacion (indices) de param gii o kij a ajustar.
	INTEGER IDG(8,2)
!
!	Nro. de param. identicos a c/gii o kij
	INTEGER NIDG(8)
!	Identificacion (indices) de IDENTicos a gii ó kij
	INTEGER IDIDG(104,2)
!
!
!	           8 X 13 = 104
!
!	Identificación parámetros de área
	INTEGER IDQ(8)
!
!	Identif (indices) de parametros ALFAij
	INTEGER IDALF(8,2)
!	Nro. de param. identicos a c/alfaij
	INTEGER NIDALF(8)
!	Identificacion (indices) de IDENTicos
	INTEGER IDIDALF(104,2)
!
!	IDENTif (indices) de parametros ENASSijkl
	INTEGER IDENASS(8,4)
!	Nro. de param. IDENTicos  a c/ENASSijkl
	INTEGER NIDENASS(8)
!	IDENTif (indices) de IDENTicos a ENASSijkl
	INTEGER IDIDENASS(104,4)
!
!	Identif (indices) de param. RKASSijkl
	INTEGER IDRKASS(8,4)
!	Nro. de param. identicos a c/RKASSijkl
	INTEGER NIDRKASS(8)
!	Identif (indices) de identicos a c/RKASSijkl
	INTEGER IDIDRKASS(104,4)
!
!	Índices para dc
	INTEGER IDDC(8)
!	Indice de los parametros temperatura critica
      INTEGER idtc(8)
!
	DOUBLE PRECISION XX(8)
	INTEGER IDXX(8)
!
!	IDXX(i) =	identificación del tipo de parámetro a ser ajustado.
!	            IDXX se utiliza para saber que factor de escalado
!	            debe aplicarse al parámetro que se esta ajustando.
!
!	IDXX(i) = 1   g*
!	          2   k*ij
!	          3   g'ij / g''ij / k'ij   coeficientes var. Con temp
!	          4   q (param. de área)
!	          5   alfa  (param. de no azar)
!	          6   épsilon (energ. asoc.)
!	          7   kappa (vol. asoc.)
!	          8   dc (diámetro critico)
!	          9   Tc (temperatura crítica)
!
	DOUBLE PRECISION XJAC(260,8), XJTJ(36), WORK(500)
	DOUBLE PRECISION FCION(260)
	DOUBLE PRECISION FCAUX(260)
	DOUBLE PRECISION FCAUX2(260)
	DOUBLE PRECISION FCAUX3(260)
!
!	Presion de vapor(atm) calculada a la temperatura Tsat(i)
	DOUBLE PRECISION PsTclc(20)
!
!	Doble precisionamiento de variables calculadas con el modelo
	DIMENSION PCALC(NDM), XCALC(NDM), YCALC(NDM),xCalc2(NDM),yCalc2(NDM)
	DIMENSION ZLCALC(NDM), ZVCALC(NDM)
	dimension rgamcalc1(NDM),rgamcalc2(NDM),rgamcalc3(NDM)
	DOUBLE PRECISION ZLicalc(NDM)
!
!	Doble precisionamiento desviaciones stándar para mezclas
	dimension STDVA(25),STDVB(25),STDVC(25),stDvD(25)
!
!	A partir de esta versión, A-Paraest está en Fortran 90, lo que me permite + de 73 columnas.
!	Aún así separo en renglones cerca de la col. 100 para no restar legibilidad.
!	La separación como se ve, es con un "&" al final, y otro "&" al principio.
	INTEGER NC,NG,NGA,nga2,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,NRKASS,NPAR,NCXP,IOPT,&

&	        NUMDC,NUMTC,NOPCION,MAXFN,IXJAC,INFER,IER

	character(1) OPDER
!
	double precision W1,W2,W3,W4,W5,W6,W7,W8,EPS,DELTAZ,SSQ
!
!	COMMON Para subrutina FUNC requerida por IMSL
	COMMON/FORFUNC/NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,Tspline,RKIJ,RKTIJ,ALFAIJ,MASSOC,ENASS,RKASS,&

&	               W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,ISAT,IGLE,IVE,IUNI,NNP,NNQ,NALFA,NENASS,NRKASS,NSAT,     &

&	               IDC,NOBS,IHEN,mez,IFL,Texp,Pexp,xExp,xExp2,YExp,yExp2,IDSAT,IDG,ITEMP,NIDG,IDIDG,&

&	               IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,IDRKASS,NIDRKASS,IDIDRKASS,  &

&	               DC,PSTCLC,XCALC,YCALC,xCalc2,yCalc2,ZLCALC,ZVCALC,FCAUX,PCALC,IOBJSAT,FCAUX2,fcaux3,&

&	               STDVA,STDVB,stdvc,NumTc,idTc,numDC,IDDC,NGA
!
	COMMON/FORFUNC2/rgamcalc1,rgamcalc2,rgamcalc3,rgamexp1, rgamexp2, rgamexp3,rhoExp
	COMMON/FORFUNC3/ONUASOC,zlicalc,PRCRIT,TEMPCRIT,TSAT,PSAT
	common/forfunc4/prop_E_calc(NDM),prop_E_exp(NDM)
!
	COMMON/IDSC/IDXX
!
!	Límites:
!
!		10 Grupos diferentes en la contribucion dispersiva
!		12 Grupos diferentes en la contribucion de asociacion
!		20 Componentes
!		180 Datos experimentales de cualquier tipo, sea  de presion de vapor de comp puros o cualquier
!		otro tipo de dato de equilibrio
!		25 Data sets de ELV
!		8  (Ajuste de a lo sumo 8 paráetros diferentes)
!		13  (A lo sumo 13 parametros IDENTicos a alguno de los ajustados)
!
!
!	IDC(J,1)  = Contiene el numero que identifica al componente para el que se dan las fracciones molares ex-
!	            perimentales, en el dataset numero j.
!
!	IDC(J,2) = Contiene el numero que identifica al componente para el que NO se dan las fracciones molares
!	           experimentales, en el dataset numero j.
!
!
	opder = "N"


!	Lectura desde el archivo de entrada "paraest.dat"
	CALL READ(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&

&	          W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,                               &

&	          NNP,NNQ,NALFA,NENASS,NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,                     &

&	          NSAT,IDC,NOBS,iHEN,mez,IFL,NTEXT,rGamExp1,rGamExp2,rGamExp3,prop_E_calc,prop_E_exp,   &

&	          TEXP,PEXP,XEXP,xexp2,YEXP,yExp2,rhoExp,XX,IDSAT,IDG,ITEMP,NIDG,IDIDG,iDQ,IDALF,NIDALF,IDIDALF,    &
&	          IDENASS,NIDENASS,IDIDENASS,IDRKASS,NIDRKASS,IDIDRKASS,                               &

&               NCXP,IOPT,IOBJSAT,NUMDC,IDDC,numTc, idTc)
!
!
!	Inicio de DC en base al punto critico predicho por la GC-EoS, para todo componente (no tiene en cuenta asociación)
	CALL INICDC(NC,PRCRIT,TEMPCRIT,DC)
!
!	De ahora en más, los dc iniciales (cuando se ajustan) se cargan abajo.
!
!	Cuando numDc no es cero: carga de los valores de dc iniciales en XX
!	CALL FORDC(NUMDC,IDDC,DC,NPAR,XX)ç
!
!	Calculo de las temperaturas de spline para los valores iniciales de los parametros dispersivos.
      CALL SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
!
!	Escribe PARTE de los datos en archivo de salida paraest.out
	CALL WRITESAL1(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&

&	               W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,                &

&	               NNP,NNQ,NALFA,NENASS,NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,                     &

&                    NSAT,IDC,NOBS,IHEN,mez,IFL,NTEXT,TEXP,PEXP,XEXP,xexp2,YEXP,XX,                       &

&                    IDSAT,IDG,ITEMP,NIDG,IDIDG,IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,      &
&                    IDRKASS,NIDRKASS,IDIDRKASS,DC,IOPT,NUMDC,IDDC,TSPLINE,numTc,idTc)
!
!	Identificacion de los valores contenidos en XX, para luego poder asignarles factores de escalado adecuados
	CALL FORSCAL(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,numTc,NUMDC)
!
!	Escalado de los valores de xx
	CALL XXSCAL(NPAR,XX,IDXX)
!
!	Cálculos de simulacion o Ajuste
!	----------------------------------------------------
	IF(NPRED.EQ.1)THEN
!	  Prediccion (o simulacion)
!
!	  Chequea derivadas:
	  IF (OPDER.EQ.'S')NOPCION = 1
!	  No chequea deriv
	  IF (OPDER.EQ.'N')  NOPCION = 0
!
	  CALL FUNC(XX,NCXP,NPAR,FCION)
!
	ELSE
!	  Ajuste de Parametros
!
!	  ***********************************Optimizacion******************************
!
!       Parametros para la rutina ZXSSQ
	  EPS = 1.D-6
	  DELTAZ = 1.D-6
	  MAXFN = 100
	  PARM(3)=2.D2
	  PARM(4)=.1D0
!	  Comentario      IXJAC = 180
	  IXJAC = NCXP
!       No chequeo de deriv.
	  NOPCION = 0
!
!	  Llamado a la subrutina de optimización.
	  CALL ZXSSQ(FUNC,NCXP,NPAR,NSIG,EPS,DELTAZ,MAXFN,IOPT,PARM,XX,SSQ,FCION,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
!
!	  Escribe estadisticas de la optimizacion
	  CALL WRITESTAT(SSQ,XJAC,WORK,INFER,IER)
!
!	  Cálculo de las temperaturas de spline para los valores FINALES de los parametros dispersivos
	  CALL SPLINECALC(NG,TSTR,GSTR,GDOT,GDDOT,TSPLINE)
	  IF (OPDER.EQ.'S')NOPCION = 1
!
!	  Simulacion con param. finales y chequeo de derivadas
	  IF (OPDER.EQ.'S')  CALL FUNC(XX,NCXP,NPAR,FCION)
!
	END IF
!
!
!	Escribe parametros finales, y valores experimentales y calculados.
	CALL WRITESAL2(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&

&	               ISAT,IGLE,IVE,IUNI,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,NTEXT,                       &

&	               TEXP,PEXP,XEXP,xexp2,YEXP,yExp2,IDSAT,DC,tempcrit,XCALC,YCALC,xCalc2,yCalc2,ZLCALC,ZVCalc,FCAUX,PCALC,  &

&                    IOBJSAT,FCAUX2,fcaux3,STDVA,STDVB,stDvC,stDvD,rgamcalc1, rgamcalc2,rgamcalc3,numTc,idTc,   &
&	               NUMDC,IDDC,TSPLINE,rgamexp1, rgamexp2,rgamexp3,zLicalc,prop_E_calc,prop_E_exp)
!
ENDprogram
!
!	Fin del Paraest.
!
!-------------------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!	Subrutina de lectura del programa Paraest. Lee la información del archivo de entrada paraest.dat. Al mismo
!	tiempo, imprime lo que va leyendo con formato, para -ayudar a- identificar un error en la entrada, de
!	haberlo.
!
SUBROUTINE READ(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&

&	                W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,                &

&	                NNP,NNQ,NALFA,NENASS,NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,                     &

&	                NSAT,IDC,NOBS,iHEN,mez,IFL,NTEXT,rGamExp1,rGamExp2,rGamExp3,prop_E_calc,prop_E_exp,  &

&	                TEXP,PEXP,XEXP,xexp2,YEXP,yExp2,rhoExp,XX,IDSAT,IDG,ITEMP,NIDG,IDIDG,iDQ,IDALF,NIDALF,IDIDALF,    &
&	                IDENASS,NIDENASS,IDIDENASS,IDRKASS,NIDRKASS,IDIDRKASS,                               &

&                     NCXP,IOPT,IOBJSAT,NUMDC,IDDC,numTc, idTc)
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
	parameter(NCM=20,NDM=260)
!
!	DIM DISPERSION
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
!
!	DIM ASOCIACION
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,20)
!
!	DIM AJUSTE DE PARAM
	DIMENSION PARM(4)
!
	DIMENSION IDSAT(20)
!
	DIMENSION  TEMPCRIT(20),PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260),xexp2(260),YEXP(260),yExp2(260),rhoExp(NDM)
	dimension rgamexp1(260),rgamexp2(260),rgamexp3(260)
	dimension prop_E_calc(NDM),prop_E_exp(NDM)
!
!	double precision  DC(20)
!
!
!	Dimensionamiento  nro. de datos de pvap comp puros
	DIMENSION NSAT(20)
	DIMENSION IOBJSAT(20)
!
!	Dimensionamientos que tienen que ver con los datasets
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
!
!		ID 	==> identificacion o identicos
!		IDID	==> identificacion de identicos
!
!	Los parametros se identifican a traves de su naturaleza y de sus inidices
	DIMENSION ITEMP(8)
!
!	Identificacion (indices) de param gii o kij a ajustar.
	DIMENSION IDG(8,2)
!	Nro. de param. IDENTicos a c/gii o kij
	DIMENSION NIDG(8)
!
!	Identificacion (indices) de IDENTicos a gii o kij
	DIMENSION IDIDG(104,2)
!
!                8 X 13 = 104
!
!	IDENTif param de area
	DIMENSION IDQ(8)
!
!	Identif (indices) de parametros ALFAij
	DIMENSION IDALF(8,2)
!	Nro. de param. IDENTicos a c/alfaij
	DIMENSION NIDALF(8)
!	Identificacion (indices) de IDENTicos
	DIMENSION IDIDALF(104,2)
!
!	IDENTif (indices) de parametros ENASSijkl
	DIMENSION IDENASS(8,4)
!	Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION NIDENASS(8)
!	IDENTif (indices) de IDENTicos a ENASSijkl
	DIMENSION IDIDENASS(104,4)
!
!	Identif (indices) de param. RKASSijkl
	DIMENSION IDRKASS(8,4)
!	Nro. de param. identicos a c/RKASSijkl
	DIMENSION NIDRKASS(8)
!	Identif (indices) de identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4),IDDC(8),idTc(8)
!
	DIMENSION XX(8)
	common/escMax/vMaxAlfa,vMaxkp,vMaxKstr,vMinGp,vMaxGp,error_Dc(NCM),dc_medio(NCM)
	dimension error_Dc1(NCM)
!
	OPEN(UNIT=1,TYPE='OLD',NAME='paraest.dat',FORM='FORMATTED',ACCESS='SEQUENTIAL')
	READ(1,*) W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,PARM(1),PARM(2),IOPT
	write(*,'(8F10.4,3I2)')W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,PARM(1),PARM(2),IOPT
!
1	FORMAT(8F10.4)
2	FORMAT(20I3)
!
	READ(1,*) NC,NG,NGA,NSIG
	write(*,'(4I3)')NC,NG,NGA,NSIG
	READ(1,*) ISAT,IGLE,IVE,IUNI,NPRED
	write(*,'(5I3,/)')ISAT,IGLE,IVE,Iuni,NPRED
!
!	Cuenta el numero de datos de equilibrio
	NCXP = 0
	IF(ISAT.NE.0)then

	  READ(1,*) (IDSAT(I), I=1,ISAT)
	  write(*,'(<iSat>(I3,X))') (IDSAT(I), I=1,ISAT)
	  READ(1,*) (NSAT(IDSAT(I)), I=1,ISAT)
	  write(*,'(<iSat>(I3,X))')(NSAT(IDSAT(I)), I=1,ISAT)
	  READ(1,*) (IOBJSAT(IDSAT(I)), I=1,ISAT)
	  write(*,'(<iSat>(I3,X),/)') (IOBJSAT(IDSAT(I)), I=1,ISAT)
!
	  DO I=1,ISAT

	    NSUMK = 0
	    do K=1,(IDSAT(I)-1)

	      NSUMK = NSUMK + 1

	    end do
	    NSUML = 0
	    do L=1,NSUMK

	      NSUML = NSUML + NSAT(L)

	    end do
	    DO J= (NSUML +1 ),( NSUML + NSAT(IDSAT(I)) )

	      if(iObjSat(idSat(i)).NE.81)then
	        READ(1,*) TEXP(J),PEXP(J)
	        write(*,'(F7.2,F12.8)') TEXP(J),PEXP(J)
	      else
	        READ(1,*) TEXP(J),PEXP(J),rhoExp(j)
	        write(*,'(F7.2,F12.8)') TEXP(J),PEXP(J),rhoExp(j)
		endif

	    END DO
!
	    NCXP = NCXP + NSAT(IDSAT(I)) !Nro. de datos YA escritos
!
	  END DO

	endif
	write(*,'(//)')
!	Lectura de áreas y T de referencia de cada grupo dispersivo
	DO I=1,NG

	  READ(1,*) QAREA(I), TSTR(I)
	  write(*,'(F6.3,F7.2)')qArea(i),Tstr(i)

	END DO
!	Lectura de energía de segmento de grupo
	READ(1,*) (GSTR(I), I=1,NG)
	write(*,'(/,<NG>(F7.0,2X),"   !g*")')(GSTR(I), I=1,NG)
!
	READ(1,*) (GDOT(I), I=1,NG)
	write(*,'(<NG>(F7.4,2X),"   !g`")')(GDOT(I), I=1,NG)
!
	READ(1,*) (GDDOT(I), I=1,NG)
	write(*,'(<NG>(F7.4,2X),"   !g``",//)')(GDdOT(I), I=1,NG)
!
!	Parámetros de interacción binaria: k*
	DO I=1,NG

	  READ(1,*) (RKIJ(I,J), J=1,NG)
	  if(i.EQ.1)then

	    write(*,'(/,<NG>(F6.4,2X),"  !k*")')(RKIJ(I,J), J=1,NG)

	  else

	    write(*,'(<NG>(F6.4,2X))')(RKIJ(I,J), J=1,NG)

	  endif


	END DO
!	k': dependencia con la T de k*
	DO I=1,NG

	  READ(1,*) (RKTIJ(I,J), J=1,NG)
	  if(i.EQ.1)then

	    write(*,'(/,<NG>(F6.4,2X),"  !k`")')(RKTIJ(I,J), J=1,NG)

	  else

	    write(*,'(<NG>(F6.4,2X))')(RKTIJ(I,J), J=1,NG)

	  endif

	END DO
!	alfa: parámetros de "no-azar"
	DO I=1,NG

	  READ(1,*) (ALFAIJ(I,J), J=1,NG)
	  if(i.EQ.1)then

	    write(*,'(/,<NG>(F6.4,2X),"  !alfa")')(alfaIJ(I,J), J=1,NG)

	  else

	    write(*,'(<NG>(F6.4,2X))')(alfaIJ(I,J), J=1,NG)

	  endif

	END DO
!
	IF(NGA.GT.0)THEN
!
!	  Contribución asociativa
!	  Lectura de vector de sitios de asociación:
	  READ(1,*) (MASSOC(I),I=1,NGA)
	  write(*,'(//,<NGA>(I2,2X),5X,"! Vector de sitios de asociación",/)') (MASSOC(I),I=1,NGA)

!	  Lectura de energía y volumen de asociación en una misma línea
	  NST = 0
	  EnAss(:6,:12,:6,:12) = 0D0
	  rKAss(:6,:12,:6,:12) = 0D0
	  do i = 1,NGA
	    do j = i,NGA
	      do k = 1,Massoc(i)
	        do l = k,Massoc(j)
!
	          NST = NST + 1
	          read(1,*)EnAss(k,i,l,j),rKAss(k,i,l,j)
		    if(NST.EQ.1)then

	            write(*,'(F6.1,4X,F8.5,4X,"! ",2I2,3X,2I2,4X,"Energía y volúmen de asociación")')EnAss(k,i,l,j),rKAss(k,i,l,j)

	          else

	            write(*,'(F6.1,4X,F8.5,4X,"! ",2I2,3X,2I2,4X)')EnAss(k,i,l,j),rKAss(k,i,l,j)

	          endif
!---------------Simetrías
	          EnAss(l,i,k,j) = EnAss(k,i,l,j)
	          EnAss(l,j,k,i) = EnAss(k,i,l,j)
	          EnAss(k,j,l,i) = EnAss(k,i,l,j)
	          rKAss(l,i,k,j) = rKAss(k,i,l,j)
	          rKAss(l,j,k,i) = rKAss(k,i,l,j)
	          rKAss(k,j,l,i) = rKAss(k,i,l,j)
!
	        enddo
	      enddo
	      write(*,'(/)')
	    enddo
	  enddo
	endif
!-----------------------------------------------------------------------
!
	write(*,'(//)')
!	Lectura de propiedades de los compuestos ingresados:
	DO I=1,NC

	  READ(1,*) TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I)		!,DC(I)
	  write(*,'(2(F7.2,2X,F12.6,2X),"! Comp. Nro.",I2)')TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I),i	!,DC(I)

	END DO
	write(*,'(//)')
!
!	Configuración grupal de fuerzas dispersivas:
	DO I=1,NG
!
!	  NU (I,J)= numero de grupos del tipo i en el compuesto de tipo j.
	  READ(1,*) (NU(I,J), J=1,NC)
	  if(i.EQ.1)then

	    write(*,'(//,<NC>(I2,4X),"!  Grupo disp. Nro.",I2,"   Configuración grupal dispersiva")')(NU(I,J), J=1,NC),I

	  else

	    write(*,'(<NC>(I2,4X),"!  Grupo disp. Nro.",I2)')(NU(I,J), J=1,NC),I

	  endif
!
	END DO
!
!
	IF(NGA.GT.0) THEN
!
!	  Configuración grupal asociativa, si corresponde:
	  DO I=1,NGA
!
!	    ONUASOC(I,J) = número de grupos del tipo i en el compuesto de tipo j (asociación).
	    READ(1,*) (ONUASOC(I,J), J=1,NC)
	    if(i.EQ.1)then

	      write(*,'(/,<NC>(F2.0,4X),"!  Grupo asoc. Nro.",I2,"   Configuración grupal asocitiva")')(oNUAsoc(I,J), J=1,NC),I

	    else

	      write(*,'(<NC>(F2.0,4X),"!  Grupo asoc. Nro.",I2)')(oNUAsoc(I,J), J=1,NC),I

	    endif

	  END DO

	END IF
!	write(*,'(//)')
!
!	Lectura de sets de datos de equilibrio binario o ternario
	DO I=1,IGLE
!
	  read(1,*)IHEN(I)
	  IF(IHEN(I).EQ.0) THEN

	    write(*,'(//,I2,5X,"Datos de equilibrio binario")')iHen(i)
	    READ(1,*)IFL(I)
	    write(*,'(i2)')iFL(i)
	    if(iFL(i).LT.3)then

	      mez(i) = 2
	      READ(1,*) NOBS(I),IDC(I,1),IDC(I,2)
	      write(*,'(3(I2,4X))')NOBS(I),IDC(I,1),IDC(I,2)

	    else

	      mez(i) = 3
	      READ(1,*) NOBS(I),(IDC(I,j),j=1,mez(i))
	      write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,mez(i))

	    endif

	  elseif(iHen(i).EQ.1)then
!
!	    Datos de propiedades en exceso
	    read(1,*)ifl(i)
	    if(ifl(i).EQ.1)then

	      write(*,'(//,I2,5X,"Datos de volumen en exceso")')ifl(i)

	    elseif(ifl(i).EQ.2)then

	      write(*,'(//,I2,5X,"Datos de entalpía en exceso")')ifl(i)

	    endif
	    read(1,*)nObs(i),(idC(i,j),j=1,2)
	    write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)

	  ELSEif(ihen(i).EQ.5)then

	    write(*,'(//,I2,5X,"Datos coeficientes de actividad")')iHen(i)
	    read(1,*) mez(i)
	    write(*,'(I2)')mez(i)
	    read(1,*)(idC(i,j),j=1,mez(i)),nObs(i)
!	    If (mez(i)==2) then
!	      READ(1,*) IDC(I,1),IDC(I,2),NOBS(I)
!	    elseif(mez(i)==3)then
!	      READ(1,*) IDC(I,1),IDC(I,2),IDC(i,3),NOBS(I)
!	    end if
	    write(*,'(<mez(i)+1>(I2,2X))')(idC(i,j),j=1,mez(i)),nObs(i)

	  END IF
        READ(1,23) (NTEXT(I,J),J=1,65)
	  write(*,23)(NTEXT(I,J),J=1,65)
23	  FORMAT(40A2)
!
	  write(*,'(/)')
	  DO J=( NCXP+1 ),(NCXP + NOBS(I))

	    If(ihen(i).EQ.0)then
!
!	      ******************Datos de equilibrio*****************
!
	      if(mez(i).EQ.2)then

!	        Binario
	        READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J)
	        write(*,'(F7.2,3(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),YEXP(J),j

	      elseif(mez(i).EQ.3)then

!	        Ternario
	        READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J),xExp2(j),yExp2(j)
	        write(*,'(F7.2,5(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),YEXP(J),xExp2(j),yExp2(j),j

	      endif
!
	    elseif(iHen(i).EQ.1)then
!
!	      ******************Propiedades en exceso*****************
!
!	      if(ifl(i).EQ.2)READ(1,*) TEXP(J),PEXP(J),XEXP(J),VEexp(J)
!	      if(ifl(i).EQ.2)READ(1,*) TEXP(J),PEXP(J),XEXP(J),HEexp(J)
	      READ(1,*) TEXP(J),PEXP(J),XEXP(J),prop_E_exp(J)

	    elseif(iHen(i).EQ.5)then
!
!	      ******************Coeficientes de actividad*****************
!
	      if(mez(i).EQ.2)then

	        read(1,*) TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J)
	        write(*,'(F7.2,4(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),j

	      elseif(mez(i).EQ.3)then

	        read(1,*) TEXP(J),PEXP(J),XEXP(J),xexp2(J),rgamexp1(J),rgamexp2(J),rgamexp3(J)
	        write(*,'(F7.2,5(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),rgamexp3(J),j

	      endif
!
	    end if

	  END DO
	  NCXP = NCXP + NOBS(I)
!
!	  Fin del bucle sobre los DATASETS
	END DO



!
	write(*,'(//)')
	READ(1,*) NNP,NNQ,NALFA,NENASS,NRKASS,numTc, NUMDC
	write(*,'(7(I2,4X),"! Interaccion - Area - No-azar - Energ. asoc. - Volumen asoc. - Temp. crit. - diametros  crit.")') NNP,NNQ,NALFA,NENASS,NRKASS,numTc, NUMDC
!
	if(NNP.NE.0) then
!
!	  Lectura de identificación de parámetros de interacción dispertiva a ser ajustados: g(i,i) y k(i,j)
	  READ(1,*) (IDG(I,1),IDG(I,2),I=1,NNP)
	  write(*,'(<NNP>(2(I2,X),4X))')(IDG(I,1),IDG(I,2),I=1,NNP)
!
	  READ(1,*) (ITEMP(I), I=1,NNP)
	  write(*,'(<NNP>(I2,6X))')(ITEMP(I), I=1,NNP)
!
	  READ(1,*) (NIDG(I), I=1,NNP)
	  write(*,'(<NNP>(I2,6X))')(NIDG(I), I=1,NNP)
!
	  NSSS = 0
!	  Lectura de parámetros idénticos:
	  DO L=1,NNP
!
	    if(NIDG(L).NE.0) then
!
	      READ(1,*) (IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
!	      NumIde = (NIDG(L)+NSSS) - (1+NSSS)
	      write(*,'(20(I2,X))')(IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
	      NSSS = NSSS + NIDG(L)
!
	    end if
!
	  END DO

	end if
!
!	Lectura de identificación de parámetros de área a ser ajustados
	if(NNQ.NE.0)then

	  READ(1,*) (IDQ(I), I=1,NNQ)
	  write(*,'(<NNQ>(I2,X))') (IDQ(I), I=1,NNQ)

	endif
!
!	Lectura de identificación de parámetros de no-azar a ser ajustados
	if(NALFA.NE.0) then
!
	  READ(1,*) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
	  write(*,'(<nAlfa>(2(I2,X),4X))')(IDALF(I,1),IDALF(I,2),I=1,NALFA)
!
	  READ(1,*) (NIDALF(I),I=1,NALFA)
	  write(*,'(<nAlfa>(I2,6X))')(NIDALF(I),I=1,NALFA)
	  NSSA = 0
	  DO L=1,NALFA

	    if(NIDALF(L).NE.0) then

	      READ(1,*) (IDIDALF(J,1),IDIDALF(J,2),J=(1+NSSA),(NIDALF(L)+NSSA))
!	      NumIde = (NIDALF(L)+NSSA) - (1+NSSA)
	      write(*,'(20(2I2,X))')
	      NSSA = NSSA + NIDALF(L)

	    end if

	  END DO

	end if
!
!	Lectura de identiciación de parámetros de energía de asociación a ser ajustados
	if(NENASS.NE.0) then

	  DO I=1,NENASS

	    READ(1,*) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
	    write(*,'(4I3)')IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)

	  END DO
	  READ(1,*) (NIDENASS(I), I=1,NENASS )
	  write(*,'(<nEnAss>(I2,X))')(NIDENASS(I), I=1,NENASS )
!
!	  Lectura de idénticos:
	  NSSE = 0
	  DO L=1,NENASS

	    if(NIDENASS(L).NE.0) then

	      DO J=(1+NSSE),(NIDENASS(L)+NSSE)

	        READ(1,*) IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)
	        write(*,'(4I3)')

	      END DO
	      NSSE = NSSE + NIDENASS(L)

	    end if

	  END DO
!
	end if
!
!	Lectura de identiciación de paráetros de volumen de asociación a ser ajustados:
	if(NRKASS.NE.0) then
!
	  DO I=1,NRKASS

	    READ(1,*) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
	    write(*,'(4I3)') IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)

	  END DO
	  READ(1,*) (NIDRKASS(I), I=1,NRKASS )
	  write(*,'(<nrKass>(I2,X))')(NIDRKASS(I), I=1,NRKASS )
!
	  NSSRK = 0
	  DO L=1,NRKASS

	    if(NIDRKASS(L).NE.0) then

	      DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

	        READ(1,*) IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)
	        write(*,'(4I3)')IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)

	      END DO
	      NSSRK = NSSRK + NIDRKASS(L)

	    end if

	  END DO

	end if
!
!	Lectura de identiciación de compuestos a los cuales se les ajustará su Tc:
	IF(numTc.NE.0)then

	  READ(1,*) (idTc(I), I=1,numtc)
	  write(*,'(<numTc>(I2,X),/)')(idTc(I), I=1,numtc)

	endif
!
!	Ahora, para ajustar dc, lee la identificación del componente (idDc(i)), sumado al error asociado al dc de ese
!	componente (vMaxDc(i).
	IF(NUMDC.NE.0)then

	  READ(1,*) (IDDC(I), I=1,numDC)
	  READ(1,*) (error_Dc1(i), I=1,numDC)

	endif
!
!	Número de parametros
	NPAR = NNP+NNQ+NALFA+NENASS+NRKASS+numTc+NUMDC
	IF(NNP.NE.0) THEN

	  DO I = 1,NNP

	    NPAR = NPAR + ITEMP(I)

	  END DO
!	  nPar = nPar + sum(iTemp(:NPP))

	END IF
!
!	Lectura de los valores iniciales de los parametros
	READ(1,*)  (XX(I),I=1,NPAR)
!
!	Guardo los valores que sean de dc en una variable aparte, para el escalado con cota máx y mín.:
	dc_medio(nPar - numDc + 1:nPar) = xx(nPar - numDc + 1:nPar)
	error_Dc(nPar - numDc + 1:nPar) = error_Dc1(1:numDc)
!
	CLOSE(UNIT=1)
!
ENDsubroutine !FIN READ
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!
!	Esta subrutina imprime en paraest.out la lectura de los disintos parámetros alimentados en la subrutina READ
!	mas no los resultados de optimización
!
SUBROUTINE WRITESAL1(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&

&	               W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,                  &

&	               NNP,NNQ,NALFA,NENASS,NRKASS,NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,                       &

&                    NSAT,IDC,NOBS,IHEN,mez,IFL,NTEXT,TEXP,PEXP,XEXP,xexp2,YEXP,XX,                       &

&                    IDSAT,IDG,ITEMP,NIDG,IDIDG,IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,      &
&                    IDRKASS,NIDRKASS,IDIDRKASS,DC,IOPT,NUMDC,IDDC,TSPLINE,numtc,idtc)
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
!	Dim. dispersion
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
	DIMENSION TSPLINE(10)
!
!	DIm. asociacion
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,20)
!
!	Dim diametro de esferas rigidas en el punto critico
	DIMENSION DC(20)
!
!	Dim. ajuste de param
	DIMENSION PARM(4)
!
	DIMENSION IDSAT(20)
!
	DIMENSION  TEMPCRIT(20)
	DIMENSION  PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260), xexp2(260)
!
!	Dimensionamiento  nro. de datos de pvap comp puros
	DIMENSION NSAT(20)
!
!	Dimensionamientos que tienen que ver con los datasets
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
!
!	ID ==> identificacion o identicos
!	IDID ==> identificacion de identicos
!
!	Los parametros se identifican a traves de su naturaleza y de sus inidices
	DIMENSION ITEMP(8)
!	Identificacion (indices) de param gii o kij a ajust.
	DIMENSION IDG(8,2)
!	nro. de param. IDENTicos a c/gii o kij
	DIMENSION NIDG(8)
!	identificacion (indices) de IDENTicos a gii o kij
	DIMENSION IDIDG(104,2) !
!
!	           8 X 13 = 104
!
!	IDENTif param de area
	DIMENSION IDQ(8)
!
!	Identif (indices) de parametros ALFAij
	DIMENSION IDALF(8,2)
!	Nro. de param. IDENTicos a c/alfaij
	DIMENSION NIDALF(8)
!	Identificacion (indices) de IDENTicos
	DIMENSION IDIDALF(104,2)
!
!	IDENTif (indices) de parametros ENASSijkl
	DIMENSION IDENASS(8,4)
!	Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION NIDENASS(8)
!	IDENTif (indices) de IDENTicos a ENASSijkl
	DIMENSION IDIDENASS(104,4)
!	Identif (indices) de param. RKASSijkl
	DIMENSION IDRKASS(8,4)
!	Nro. de param. identicos a c/RKASSijkl
	DIMENSION NIDRKASS(8)
!	Identif (indices) de identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4)
!
	DIMENSION IDDC(8)
!
	dimension idtc(8)
!
	DIMENSION XX(8)
!
	CHARACTER*9 FECHA,HORA
	OPEN(UNIT=1,STATUS='unknown',file='paraest.out',FORM='FORMATTED')
!
	WRITE(1,1) W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,PARM(1),PARM(2)
!
1	FORMAT(/,X,' Factores de peso utilizados en la funcion objetivo:',/,                             &
&	       X,'Peso a las presiones de vapor:                   ',F10.4,/,                            &
&	       X,'Peso a la presión de burbuja de mezclas binarias:',F10.4,/,                            &
&	       X,'Peso a la temp. de burbuja de mezclas binarias:  ',F10.4,/,                            &
&	       X,'Peso a la composición del comp. 1 en la fase I:  ',F10.4,/,                            &
&	       X,'Peso a la composición del comp. 1 en la fase II: ',F10.4,/,                            &
&	       X,'Peso a la composición del comp. 2 en la fase I:  ',F10.4,' (sólo mezclas ternarias)',/ &
&	       X,'Peso a la composición del comp. 2 en la fase II: ',F10.4,' (sólo mezclas ternarias)',/ &
&	       X,'Peso al coef. de actividad del comp. 1:          ',F10.4,/                             &
&	       X,'Peso al coef. de actividad del comp. 2:          ',F10.4,/                             &
&	       X,'Peso al coef. de actividad del comp. 3:          ',F10.4,' (sólo mezclas ternarias)',/ &
&	       x,'Peso al volumen en exceso:                       ',F10.4,/                             &
&	       x,'Peso a la entalpía en exceso:                    ',F10.4,/                             &
&	       X,' PARM(1) = ',F10.4,'     PARM(2) = ',F10.4,/,                                          &
&	       X,' Algoritmo de optimizacion: IOPT = ',I3)
!
	WRITE(1,'(/,X," Relacion entre g*(i,i) y g(i,i) a la que se realiza ",/,&
&               X," el spline = 4 ")')
!
	WRITE(1,'(/,X,"--------------------------------------------------",/,&
&	          X,"------------ PARAMETROS INICIALES ----------------",/,  &
&	          X,"--------------------------------------------------",/,  &
&	          X," Parametros de grupos puros (termino  dispersivo)",/)')

	WRITE(1,3) (I, I=1,NC)
3	FORMAT(40X,'COMPONENTES',/,X,' Gr.Nro.:  T*        Q        TSPLINE ',<NC>I3,/)
!
	DO I=1,NG

	  WRITE(1,4) I,TSTR(I),QAREA(I),TSPLINE(I),(NU(I,J), J=1,NC)
4	  FORMAT(X,I3,5X,2F10.4,X,E10.4,<NC>I3)

	END DO
!	NU (I,J)= Numero de grupos del tipo i en el compuesto de tipo j
	WRITE(1,5)
 5	FORMAT(/,X,' Energías de interaccion gjj ',/,X,'(la forma funcional cambia a altas temperaturas)',/)
!
	DO I=1,NG
!
	  WRITE(1,6) I,I,GSTR(I), GDOT(I), GDDOT(I)
6	  FORMAT(X,' g(',I2,I2,') = ',F10.1,' * [ 1 + ',F10.4,' *(T/T* - 1) + ',F10.4, ' *LN(T/T*)] ',/)

	END DO
	WRITE(1,7)
 7	FORMAT(/,X,' Parametros binarios: [k*] ',/)
	DO I=1,NG

	  WRITE(1,8) (RKIJ(I,J), J=1,NG)
8	  FORMAT(X,<NG>F10.4)
	END DO
!
	WRITE(1,9)
 9	FORMAT(/,X,' Dependencia resp. a la temp.: [k`] ',/)
	DO I=1,NG

	  WRITE(1,8) (RKTIJ(I,J), J=1,NG)

	END DO
	WRITE(1,10)
 10	FORMAT(/,X,' Parametros binarios de no-azar: [alfa] ',/)
	DO I=1,NG

	  WRITE(1,8) (ALFAIJ(I,J), J=1,NG)

	END DO
!
	WRITE(1,186) NGA
 186	FORMAT(/,X,' Numero de grupos considerados en la contribucion de asociación = ',I3)
!
	IF(NGA.GT.0) THEN

	  WRITE(1,11) (I, I=1,NC)
11	  FORMAT(/,X,' Parametros de asociacion ',/,23X,'Componentes',/,' Grupo 	  # Sit asoc  ',<NC>I5,/)
	  DO I=1,NGA

	    WRITE(1,12) I, MASSOC(I), (ONUASOC(I,J), J=1,NC)
12	    FORMAT(3X,I3,9X,I3,6X,<NC>F5.2)

	  END DO
	  WRITE(1,13)
13	  FORMAT(/,X,'Parametros de energia de asociacion (K)  ',/)
!
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then

	      DO L=1,NGA

              IF(MASSOC(L).NE.0)then

	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,14) I,J,K,L,ENASS(I,J,K,L)
14	              FORMAT(X,' ENASS(',I2,I2,X,I2,I2,')  = ',F10.4)

	            END DO

	          END DO

	        endif

	      END DO

	    endif

	  END DO
!
!	  Escritura de parametros volumetricos de asociacion (cm3/mol)
	  WRITE(1,15)
15	  FORMAT(/,X,'Parametros de volumen de asociacion (cm3/mol) ',/)
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then

	      DO L=1,NGA

	        IF(MASSOC(L).NE.0)then

	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,16) I,J,K,L,RKASS(I,J,K,L)
16	              FORMAT(X,' RKASS(',I2,I2,X,I2,I2,') = ',F16.10)

	            END DO

	          END DO

	        endif

	      END DO

	    endif

	  END DO

	END IF
	WRITE(1,18)
18	FORMAT(/,X,'Datos sobre los compuestos puros')
	WRITE(1,19)
19	FORMAT(/,4X,'Tc/K      Pc/atm    Tsat/K    Psat/atm    DC',/,X,' (DC inicial estimado en base al pto critico exptl.) ',/)
!
	DO I=1,NC

	  WRITE(1,500) TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I),DC(I)
500	  FORMAT(X,4F10.4,F10.7)

	END DO
!
	WRITE(1,201)
201	FORMAT(/,X,' Parametros a ser estimados ',/)
!
	if(NNP.NE.0) then
!
	  WRITE(1,202)
202	  FORMAT(' Parametros gii o kij ')
!
	  WRITE(1,203) (IDG(I,1),IDG(I,2),I=1,NNP)
203	  FORMAT(X,<2*NNP>I3)
!
	  WRITE(1,204) (ITEMP(I), I=1,NNP)
204	  FORMAT(X,' Dependencia respecto a la temperatura ',/,X,<NNP>I3)
!
	  WRITE(1,206) (NIDG(I), I=1,NNP)
206	  FORMAT(X,' Numero de par. identicos ',/,X,<NNP>I3)
!
	  WRITE(1,208)
208	  FORMAT(X,' Parametros Identicos ')
!
	  NSSS = 0
!
	  DO L=1,NNP

	    if(NIDG(L).ne.0) then
!
	      WRITE(1,210) (IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
210	      FORMAT(X,<2*NIDG(L)>I3)
	      NSSS = NSSS + NIDG(L)
!
	    end if

	  END DO

	end if
!
	if(NNQ.NE.0) WRITE(1,211) (IDQ(I), I=1,NNQ)
211	FORMAT(/,X,' Parámetros de Area  ',/,X,<NNQ>I3)
!
	if(NALFA.ne.0) then

	  WRITE(1,213) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
213	  FORMAT(/,X,' Parametros de no azar ',/,X,<2*NALFA>I3)
!
	  WRITE(1,215) (NIDALF(I),I=1,NALFA)
215	  FORMAT(X,' Numero de par. identicos ',/,X,<NALFA>I3)
!
	  NSSA = 0
	  DO L=1,NALFA

          if(NIDALF(L).ne.0) then

	      WRITE(1,217) (IDIDALF(J,1),IDIDALF(J,2),J=(1+NSSA),(NIDALF(L)+NSSA))
217	      FORMAT(X,' Identicos : ',<2*NIDALF(L)>I3)
	      NSSA = NSSA + NIDALF(L)

	    end if

	  END DO

	end if
!
	if(NENASS.ne.0) then

	WRITE(1,219)
219	FORMAT(/,X,' Parametros de energia de asociacion ')
!
	  DO I=1,NENASS

	    WRITE(1,221) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
221	    FORMAT(X,4I3)

	  END DO
	  WRITE(1,223) (NIDENASS(I), I=1,NENASS )
223	  FORMAT(X,' Numero de par. identicos ',<NENASS>I3)
	  NSSE = 0
	  DO L=1,NENASS

	    if(NIDENASS(L).ne.0) then

	      DO J=(1+NSSE),(NIDENASS(L)+NSSE)

	        WRITE(1,225) IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)
225	        FORMAT(' Identicos: ',4I5)

	      END DO
	      NSSE = NSSE + NIDENASS(L)

	    end if

	  END DO

	end if
!
	if(NRKASS.ne.0) then

	  wRITE(1,503)
503	  FORMAT(/,X,' Parametros de Volumen de Asociacion ')

	  DO I=1,NRKASS

	    WRITE(1,229) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
229	    FORMAT(X,4I3)

	  END DO
	  WRITE(1,231) (NIDRKASS(I), I=1,NRKASS )
231	  FORMAT(X,' Numero de parametros identicos ',<NRKASS>I3)
	  NSSRK = 0
	  DO L=1,NRKASS

	    if(NIDRKASS(L).ne.0) then

	      DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

	        WRITE(1,232) IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)
232	        FORMAT(X,' Identicos: ',4I3)

	      END DO
	      NSSRK = NSSRK + NIDRKASS(L)

	    end if

	  END DO

	end if
	IF(numtc.ne.0) WRITE(1,233) (idtc(I), I=1,numTc)
233	FORMAT(/,X,' temperaturas criticas ',/,X,<numTc>I3)
	IF(NUMDC.NE.0) WRITE(1,234) (IDDC(I), I=1,NUMDC)
234   FORMAT(/,X,' Diametros de Esferas Rigidas ',/,X,<NUMDC>I3)
!
	WRITE(1,235)  (XX(I),I=1,NPAR)
 235	FORMAT(/,X,' Valores iniciales de los parametros a ajustar ',/,X,<NPAR>F13.4)
!
	WRITE(1,91)
91	FORMAT(//,X,' Curso de la optimización ',X)
	WRITE(1,93) (I,I=1,NPAR)
 93	FORMAT(//,X,' Parametro nro (ESCALADO). ',/,X,<NPAR>I15,X,'         FC OBJ ',/)
!
	CLOSE(UNIT=1)
!
	END !FIN WRITESAL1
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!
!	MODIFPAR modifica los valores de los parametros del modelo en base a los nuevos valores de los
!	parametros que se ajustan. Estos ultimos estan contenidos en XXAUX, libres de escalado.
!
!	A esta altura npar es dato
!
SUBROUTINE MODIFPAR(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&
&	                    W1,W2,W3,W4,W5,w6,W7,W8,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,NRKASS,   &
&	                    NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,                      &
&	                    NTEXT,TEXP,PEXP,XEXP,xexp2,YEXP,XXAUX,IDSAT,IDG,ITEMP,NIDG,IDIDG,                    &
&	                    iDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,IDRKASS,NIDRKASS,IDIDRKASS,numtc,&
&	                    idTc,NUMDC,IDDC,DC)
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
!	Dim. dispersion
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
!
!	Dim. asociacion
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,20)
!
!	Dim. ajuste de param.
	DIMENSION PARM(4)
	DIMENSION IDSAT(20)
	DIMENSION TEMPCRIT(20)
	DIMENSION PRCRIT(20), TSAT(20), PSAT(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260),xexp2(260)
!
!	Dimensionamiento  nro. de datos de pvap comp puros
	DIMENSION NSAT(20)
!
!	Dimensionamientos que tienen que ver con los datasets
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
!
!	ID ==> identificacion o identicos
!	IDID ==> identificacion de identicos
!
!	Los parametros se identifican a traves de su naturaleza y de sus indices
	DIMENSION ITEMP(8)
!	Identificacion (indices) de param gii o kij a ajust.
	DIMENSION IDG(8,2)
!
!	Nro. de param. IDENTicos a c/gii o kij
	DIMENSION NIDG(8)
!
!	Identificacion (indices) de IDENTicos a gii o kij
	DIMENSION IDIDG(104,2)
!
!                8 X 13 = 104
!
!	IDENTif param de area
	DIMENSION IDQ(8)
!
	DIMENSION IDALF(8,2) !identif (indices) de parametros ALFAij
	DIMENSION NIDALF(8) !nro. de param. IDENTicos a c/alfaij
	DIMENSION IDIDALF(104,2) !identificacion (indices) de IDENTicos
!
	DIMENSION IDENASS(8,4) !IDENTif (indices) de parametros ENASSijkl
	DIMENSION NIDENASS(8) !Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION IDIDENASS(104,4) !IDENTif (indices) de IDENTicos a ENASSijkl
!
	DIMENSION IDRKASS(8,4) !Identif (indices) de param. RKASSijkl
	DIMENSION NIDRKASS(8) !Nro. de param. identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4) !Identif (indices) de identicos a c/RKASSijkl
!
!
!	Dimensionamiento diametro de esferas rigidas en el punto crit.
	dimension idTc(8)
	DIMENSION DC(20)
	DIMENSION IDDC(8)
!
!	Parámetros siendo ajustados:
	DIMENSION XXAUX(8)
!
!	NPAR2 es una variable que tiene en cuenta el numero de parametros ya actualizados. Una vez que termina la ejecucion
!	de este modulo NPAR2 debe ser igual a npar
!
	NPAR2 = 0
!
!	Asignacion de valores a parametros del tipo gii o kij
	IF(NNP.NE.0) then !A

	  DO I=1,NNP

	    IF(IDG(I,1).EQ.IDG(I,2)) THEN

	      GSTR(IDG(I,1)) = XXAUX(NPAR2+1)
	      NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.0)then

	        GDOT(IDG(I,1)) = 0.D0
	        GDDOT(IDG(I,1)) = 0.D0

	      endif
	      IF(ITEMP(I).GE.1) GDOT(IDG(I,1)) = XXAUX(NPAR2+1)
	      IF(ITEMP(I).GE.1) NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.1) GDDOT(IDG(I,1)) = 0.D0
	      IF(ITEMP(I).EQ.2) GDDOT(IDG(I,1)) = XXAUX(NPAR2+1)
	      IF(ITEMP(I).EQ.2) NPAR2 = NPAR2 + 1
!
	    ELSE

	      RKIJ(IDG(I,1),IDG(I,2)) = XXAUX(NPAR2+1)
	      NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.0) RKTIJ(IDG(I,1),IDG(I,2))= 0.D0
1001	      IF(ITEMP(I).EQ.1) RKTIJ(IDG(I,1),IDG(I,2))=XXAUX(NPAR2+1)
	      IF(ITEMP(I).EQ.1) NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.2) THEN

!	        Avoid human error
	        ITEMP(I) = 1
	        GO TO 1001

		END IF
!
	    END IF

	  END DO

	endif
!
!	Identicos del del tipo gii o kij
	IF(NNP.NE.0) THEN

	  NSSS = 0
	  DO I=1,NNP

	    if(NIDG(I).ne.0) then

	      IF(IDG(I,1).EQ.IDG(I,2)) THEN

	        DO J=(1+NSSS),(NIDG(I)+NSSS)

	          GSTR(IDIDG(J,1))  = GSTR(IDG(I,1))
	          GDOT(IDIDG(J,1))  = GDOT(IDG(I,1))
	          GDDOT(IDIDG(J,1)) = GDDOT(IDG(I,1))

 	        END DO
	        NSSS = NSSS + NIDG(I)

	      ELSE

	        DO J=(1+NSSS),(NIDG(I)+NSSS)

	          RKIJ(IDIDG(J,1),IDIDG(J,2)) = RKIJ(IDG(I,1),IDG(I,2))
	          RKTIJ(IDIDG(J,1),IDIDG(J,2))= RKTIJ(IDG(I,1),IDG(I,2))

	        END DO
	        NSSS = NSSS + NIDG(I)
!
	      END IF

	    end if

	  END DO

	END IF
!
!	Asignacion de valores a parametros de area
	IF(NNQ.NE.0) THEN

	  DO I=1,NNQ

	    QAREA(IDQ(I)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
!	Asignacion de valores a parametros alfa
	IF(NALFA.NE.0) THEN

	  DO I=1,NALFA

	    ALFAIJ(IDALF(I,1),IDALF(I,2)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
!	Asignacion a identicos a los alfa
	IF(NALFA.NE.0) THEN

	  NSSA = 0
	  DO I=1,NALFA

	    if(NIDALF(I).ne.0) then

	      DO J=(1+NSSA),(NIDALF(I)+NSSA)

	        ALFAIJ(IDIDALF(J,1),IDIDALF(J,2))= ALFAIJ(IDALF(I,1),IDALF(I,2))

	      END DO
	      NSSA = NSSA + NIDALF(I)

	    end if

	  END DO

	END IF
!
!	Asignacion de valores a parametros de energia de asociacion
	IF(NENASS.NE.0) THEN
!
	  DO I=1,NENASS

	    ENASS(IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
!	Asignacion a identicos tipo energia de asociacion
	IF(NENASS.NE.0) THEN

	  NSSE = 0
	  DO L=1,NENASS

	    if(NIDENASS(L).ne.0) then

	      DO J=(1+NSSE),(NIDENASS(L)+NSSE)

	        ENASS(IDIDENASS(J,1),IDIDENASS(J,2),IDIDENASS(J,3),IDIDENASS(J,4)) = ENASS(IDENASS(L,1),IDENASS(L,2) ,IDENASS(L,3),IDENASS(L,4))

	      END DO
	      NSSE = NSSE + NIDENASS(L)

	    end if

	  END DO

	END IF
	if(NRKASS.ne.0) then

	  DO I=1,NRKASS

	    RKASS(IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	end if
!
!	Asignacion de valores de parametros de volumen de asociacion identicos a los ya estableclidos
	IF(NRKASS.ne.0) THEN

	  NSSRK = 0
	  DO L=1,NRKASS

	    if(NIDRKASS(L).ne.0) then

	      DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

	        RKASS(IDIDRKASS(J,1),IDIDRKASS(J,2),IDIDRKASS(J,3),IDIDRKASS(J,4)) = RKASS(IDRKASS(L,1),IDRKASS(L,2),IDRKASS(L,3),IDRKASS(L,4))

	      END DO
	      NSSRK = NSSRK + NIDRKASS(L)

	    end if

	  END DO

	END IF
!
!	Asignacion de valores a los parametros dc, cuando dc se usa como parametro ajustable (3ra opcion del primer comentario del prog.ppal)
	IF(numtc.NE.0) THEN

	  DO I=1,numtc

	    tempcrit(idtc(I)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF

	IF(NUMDC.NE.0) THEN

	  DO I=1,NUMDC

	    DC(IDDC(I)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
ENDsubroutine
!
!-------------------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!	En fcion aparecen los valores a partir de los cuales se arma la sumatoria de cuadrados
!
SUBROUTINE FUNC(XX,NCXP,NPAR,FCION)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	parameter(NDM=260)
!
!	Dim dispersion
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10),Tspline(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
!
!	Dim. asociacion
	DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
	DIMENSION ONUASOC(12,20)
!
!	Dimensionamiento diametro de esferas rigidas en el punto crit.
	DIMENSION DC(20)
!
!	Dim. ajuste de param
	DIMENSION IDSAT(20)
!
	DIMENSION TEMPCRIT(20)
	DIMENSION PRCRIT(20), TSAT(20), PSATu(20)
	DIMENSION TEXP(NDM), PEXP(NDM), XEXP(NDM), YEXP(NDM),xexp2(NDM),yExp2(NDM),rhoExp(NDM)
	dimension rgamexp1(NDM),rgamexp2(NDM),rgamexp3(NDM)
	DIMENSION PCALC(NDM), XCALC(NDM), YCALC(NDM),xCalc2(NDM),yCalc2(NDM)
	dimension rgamcalc1(NDM),rgamcalc2(NDM),rgamcalc3(NDM)
	dimension zlicalc(NDM)
	DIMENSION ZLCALC(NDM), ZVCALC(NDM)
!
!
	dimension feed(15),zFrac(15,4),zVal(4),iVap(1),indF(4)

!
!	***********************Viene importante dimension variable**************************
	DIMENSION FCION(NCXP),FCAUX(NDM),FCAUX2(NDM),FCAUX3(NDM)
!
!
!	Dimensionamiento  nro. de datos de pvap comp puros
	DIMENSION NSAT(20),IOBJSAT(20)
!
!	Dimensionamientos que tienen que ver con los datasets
!	Dim. ajuste de param.
	DIMENSION PARM(4),NTEXT(25,65)
	DIMENSION IDC(25,3),NOBS(25),IHEN(25), mez(25),IFL(25)
!
!	ID ==> identificacion o identicos
!	IDID ==> identificacion de identicos
!
!	Los parametros se identifican a traves de su naturaleza y de sus indices
	DIMENSION ITEMP(8)
	DIMENSION IDG(8,2) ! identificacion (indices) de param gii o kij a ajust.
	DIMENSION NIDG(8)   !nro. de param. IDENTicos a c/gii o kij
	DIMENSION IDIDG(104,2) !identificacion (indices) de IDENTicos a gii o kij
!
!                8 X 13 = 104
!
	DIMENSION IDQ(8)  !IDENTif param de area
!
	DIMENSION IDALF(8,2) !identif (indices) de parametros ALFAij
	DIMENSION NIDALF(8) !nro. de param. IDENTicos a c/alfaij
	DIMENSION IDIDALF(104,2) !identificacion (indices) de IDENTicos
!
	DIMENSION IDENASS(8,4) !IDENTif (indices) de parametros ENASSijkl
	DIMENSION NIDENASS(8) !Nro. de param. IDENTicos  a c/ENASSijkl
	DIMENSION IDIDENASS(104,4) !IDENTif (indices) de IDENTicos a ENASSijkl
!
	DIMENSION IDRKASS(8,4) !Identif (indices) de param. RKASSijkl
	DIMENSION NIDRKASS(8) !Nro. de param. identicos a c/RKASSijkl
	DIMENSION IDIDRKASS(104,4) !Identif (indices) de identicos a c/RKASSijkl
!
	DIMENSION idTc(8), IDDC(8)
!
!
!	*************La dimension variable que viene es muy importante que sea variable************
!
	DIMENSION XX(NPAR), IDXX(8) , XXAUX(8)
!
!	Dimensionamiento de arreglos de trabajo necesarios en el calculo de propiedades de compuestos puros.
	DIMENSION WDC(1), WTEMPCRIT(1), WPRCRIT(1),NUW(10,1),ONUASOCW(12,1),wZ(2)
!
	DIMENSION PSTCLC(20)    !Presion de vapor(atm) calculada a la temperatura TSAT(I)
!
!	Dimensionamiento de arreglos de trabajo necesarios en el calculo de propiedades de mezclas binarias.
	DIMENSION WWDC(3),WWTEMPCRIT(3), WWPRCRIT(3),NUWW(10,3),ONUASOCWW(12,3)
!
!	Dimensionamiento desviaciones standar para mezclas
	DIMENSION STDVA(25),STDVB(25),STDVc(25),stDvD(25)

	DIMENSION  :: PM(20),WWPM(2)
!
!	COMMON Para subrutina func requerida por IMSL
	COMMON/FORFUNC/NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,Tspline,RKIJ,RKTIJ,ALFAIJ,MASSOC,ENASS,RKASS,      &

&	               W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,ISAT,IGLE,IVE,IUNI,NNP,NNQ,NALFA,NENASS,NRKASS,NSAT,           &

&	               IDC,NOBS,IHEN,mez,IFL,Texp,Pexp,xExp,xExp2,YExp,yExp2,IDSAT,IDG,ITEMP,NIDG,IDIDG,      &

&	               IDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,IDRKASS,NIDRKASS,IDIDRKASS,        &

&	               DC,PSTCLC,XCALC,YCALC,xCalc2,yCalc2,ZLCALC,ZVCALC,FCAUX,PCALC,IOBJSAT,FCAUX2,fcaux3,   &

&	               STDVA,STDVB,stdvc,NumTc,idTc,numDC,IDDC,NGA
!
	COMMON/FORFUNC2/rGamCalc1, rGamCalc2,rGamCalc3,rGamExp1, rGamExp2, rGamExp3,rhoExp
!
	COMMON/FORFUNC3/ONUASOC,zLicalc,PRCRIT,TEMPCRIT,TSAT,PSATu
!
	common/forfunc4/prop_E_calc(NDM),prop_E_exp(NDM)

	COMMON/IDSC/IDXX
!
	COMMON/ESCALADO/FGGG,FRKIJ,FTEMP,FQAREA,FALFA,FENASS,FRKASS,FTc,FDc
!
	PENALIZ = 1
!
!	Desescalado de los valores de xx para poder realizar el calculo de los terminos de la funcion objetivo.
!	El resultado del desescalado se guarda en XXAUX, dado que la rutina func no debe alterar los valores de
!	XX (ver instrucciones de la libreria IMSL)
	CALL XXDESSCAL(NPAR,XX,IDXX,FGGG,FRKIJ,FTEMP,FQAREA,FALFA,FENASS,FRKASS,XXAUX,ftc,FDC)
!
!	Modificacion de los parametros en base a los nuevos valores de los mismos, que entraron escalados en en
!	el vector XX, y se se encuentran libres de escalamiento en el vector XXAUX
 	CALL MODIFPAR(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&
&	              W1,W2,W3,W4,W5, w6,W7,W8,NSIG,ISAT,IGLE,IVE,IUNI,NPRED,NNP,NNQ,NALFA,NENASS,NRKASS,   &
&	              NPAR,PARM,PRCRIT,TEMPCRIT,TSAT,PSATu,NSAT,IDC,NOBS,IHEN,mez,IFL,                     &
&	              NTEXT,TEXP,PEXP,XEXP,xexp2,YEXP,XXAUX,IDSAT,IDG,ITEMP,NIDG,IDIDG,                    &
&	              iDQ,IDALF,NIDALF,IDIDALF,IDENASS,NIDENASS,IDIDENASS,IDRKASS,NIDRKASS,IDIDRKASS,numtc,&
&	              idTc,NUMDC,IDDC,DC)
!
!
!	********************************Comienza evaluacion de FCION****************************************
!
!	Calculo de dc (de manera de ajustar exactamente un punto de la curva de presion de vapor, para cada com-
!	ponente.
!

    DO I=1,NC

	  iErrorDc = 0
!****************************
!	  if(DC(I).EQ.0) then ! Fran. 26/05/09
!****************************
 	  IF(TSAT(I).NE.0) THEN
!
!	    Verifico si el DC(I) no ha sido ya considerado como una variable de optimizacion:
	    DO JAM=1,NUMDC

	      IF(IDDC(JAM).EQ.I)iErrorDc=1	!GO TO 1325

	    END DO
	    if(iErrorDc.EQ.0)then
!
!	      Carga de las variables de trabajo
	      WDC(1)= DC(I)
	      WTEMPCRIT(1)= TEMPCRIT(I)
 	      WPRCRIT(1)  = PRCRIT(I)
	      DO J=1,NG

	        NUW(J,1) = NU(J,I)

	      END DO
	      DO K=1,NGA

	        ONUASOCW(K,1) = ONUASOC(K,I)

	      END DO
	      WTEMP = TSAT(I)
	      WPSAT = PSATu(I)
!
!	      Ajuste de DC
	      CALL ASSOCDCADJUSTFASTSUB(1,WDC,WTEMPCRIT,NG,NUW,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,      &
&	                                NGA,ONUASOCW,MASSOC,ENASS,RKASS,WTEMP,WPSAT,20,IERRDC2)
!
!	      20 es el numero maximo de iteraciones en ASSOCDCADJUSTFASTSUB
!
!	      Si hubo un error en el ajuste:
	      IF(IERRDC2.NE.0) THEN

	        WDC(1)= DC(I)

!	        Ajuste de DC por un metodo un poco mas robusto pero unas 5 veces mas lento
	        CALL AssocDcAdjustSUB(1,WDC,WTEMPCRIT,NG,NUW,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,         &
&	                              NGA,ONUASOCW,MASSOC,ENASS,RKASS,WTEMP,WPSAT,2,PSATCALC)
!	        2  es el numero maximo de iteraciones ASSOCDCADJUSTSUB
	        PSTCLC(I) = PSATCALC

	      END IF
	      DC(I) = WDC(1)
	      WRITE(6,'(" Fin ajuste de dc comp nro. = ",I2)')I

	    endif

	  END IF
!********************************
!*	  end if		 ! Fran. 26/05/09
!********************************
	END DO
!
!
!	Calculo de presiones de vapor o de isofugacidad de puros
!	Numero de Datos experimentales
	NUMEXP=0
	DO I=1,NC

	  IF(NSAT(I).NE.0) THEN
!               ------------------------------------------------
	    WDC(1)= DC(I)
	    WTEMPCRIT(1)= TEMPCRIT(I)
 	    WPRCRIT(1)  = PRCRIT(I)
	    DO J=1,NG

	      NUW(J,1) = NU(J,I)

	    END DO
	    DO K=1,NGA

	      ONUASOCW(K,1) = ONUASOC(K,I)

	    END DO
!
!	    Calculo de Presiones de Vapor O ISOFUGAC. DE PUROS
!
!         **************************Bucle Nsat***********************
!
		DO J=1,NSAT(I)
!
	        WTEMP=TEXP(NUMEXP+1)
	        WPSAT=PEXP(NUMEXP+1)
	        IF(IOBJSAT(I).EQ.50) THEN

!	          Función objetivo: isofugacidad
	          write(6,'("Punto = ",I3," del compuesto",I3)')j,i
	          CALL assocIsofPurosSub(WTEMPCRIT,NG,NUW,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOCW, &
&	                                 MASSOC,ENASS,RKASS,WDC,WTEMP,WPSAT,WZFASL,WZFASV,APART,IERRYAF)
!
	          IF(IERRYAF.EQ.0) THEN

	            PCALC(NUMEXP+1)=WPSAT
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
!

	        END IF
	        IF(IOBJSAT(I).GE.80) THEN
!
!	          Función objetivo: presión de vapor.
	          write(6,'("Punto = ",I3," del compuesto",I3)')j,i
	          call fakePARGR(1,NG,NGA,NST,IOUt,NUW(:NG,1),QAREA(:NG),TStr(:NG),gStr(:NG),GDOT(:NG),GDDOT(:NG), &
&			                 RKIJ(:NG,:NG),RKTIJ(:NG,:NG),ALFAIJ(:NG,:NG),ONUASOCW(:NGA,1),MASSOC(:NGA),         &
&	                         ENASS(:2,:NGA,:2,:NGA),RKASS(:2,:NGA,:2,:NGA),Tspline(:NG))
	          call fakePARMOL(1,NG,NST,IOUt,WDC,WTEMPCRIT,wPRCRIT)
	          call paraGC(wTemp,1,NG,NST,0)
	          iGP = 1
	          iGz = 0 !debería usar el resultado de ajuste del dc...
	          iErrPSat = 0
			  iComp = 1
	          call PSAT(wTemp,wPSat,wZ,1,NG,NST,iComp,IGP,IGZ,iErrPSat)
!	          CALL ASSOCPSATSUB(WTEMPCRIT,NG,NUW,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,      &
!&	                            NGA,ONUASOCW,MASSOC,ENASS,RKASS,WDC,WTEMP,WPSAT,WZFASL,WZFASV,IERRPSAT)
!
	          IF(IERRPSAT.GE.0) THEN
!
	            PCALC(NUMEXP+1)=WPSAT
	            ZLCALC(NUMEXP+1)= WZ(2)
	            ZVCALC(NUMEXP+1)= WZ(1)
 	            FCAUX(NUMEXP+1)= (PCALC(NUMEXP+1)/PEXP(NUMEXP+1)-1)
	            FCION(NUMEXP+1)= W1 * FCAUX(NUMEXP+1)
	            if(iObjSat(i).EQ.81)fcion(numExp+1) =fcion(numExp+1) + abs(wPsat/82.05/wTemp/wZ(2)/rhoExp(numExp+1) - 1)

	          ELSE
!
	            PCALC(NUMEXP+1)= 0.
	            ZLCALC(NUMEXP+1)= 0.
	            ZVCALC(NUMEXP+1)= 0.
 	            FCAUX(NUMEXP+1)= 0.
	            FCION(NUMEXP+1)= PENALIZ
!
	          END IF

	        END IF
	        NUMEXP = NUMEXP + 1
!
	      END DO
!
!         **********************Fin bucle NSat***********************
	  END IF

	END DO
!
!
!	******************Comienza bucle sobre los datasets************************
!
	IF(IGLE.NE.0) THEN
!
!	************************Datasets*****************
	  DO KS=1,IGLE

!	    Carga de Datos para los 2 componentes del binario (3 en el caso de cálculo de gammas)
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
	    NCw = 2
	    If(mez(ks)==3)then

	      NCw = 3
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
!
!---------------------------------------------------------------------------------->
!
!
!	    Esta subrutina falsea la versión de PARGR del paquete GCA, teniendo dentro de sí los COMMON's
!	    necesarios para los parámetros grupales y de compuesto puro que necesita el paquete.
	    call fakePARGR(NCw,NG,NGA,NST,IOUt,NUWW(:NG,:NCw),QAREA(:NG),TStr(:NG),GStr(:NG),GDOT(:NG),GDDOT(:NG), &
&	                   RKIJ(:NG,:NG),RKTIJ(:NG,:NG),ALFAIJ(:NG,:NG),ONUASOCWW(:NGA,:NCw),MASSOC(:NGA),         &
&	                   ENASS(:2,:NGA,:2,:NGA),RKASS(:2,:NGA,:2,:NGA),Tspline(:NG))
!
!	    Esta hace lo mismo que la anterior, sólo que para los parámetros por componente que define
!	    originalmente PARMOL:
	    call fakePARMOL(NCw,NG,NST,IOUt,WWDC(:NCw),WWTEMPCRIT(:NCw),wwPRCRIT(:NCw))
!
!<--------------------------------------------------------------------------------
!
	    IF(IHEN(KS).EQ.0) THEN   ! TPXY (IHEN=0)
!
	      NFAILURE = NOBS(KS)
	      STDVA(KS) = 0D0
	      STDVB(KS) = 0
	      DO I=1,NOBS(KS)

	        TEMP  = TEXP(NUMEXP+1)
	        PRESS = PEXP(NUMEXP+1)
	        X1 = XEXP(NUMEXP+1)
	        Y1 = YEXP(NUMEXP+1)
!	        Nota: ifl = -1 o ifl = -2
	        IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then
!
	          write(6,*) "i=",i
!	          Inicio de la variable de error. Si no se inicia puede dar cualquier cosa.
	          iErrMixIsof = 0
	          CALL ASSocISOFmezclasSUB(2,WWDC,WWTEMPCRIT,NG,NUWW,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,  &
&	                                   NGA,ONUASOCWW,MASSOC,ENASS,RKASS,TEMP,PRESS,X1,Y1,                       &
&	                                   APART1,APART2,ZLL,ZVV,IERRMIXISOF)
!
	          IF(IERRMIXISOF.EQ.0) THEN

	            ZLCALC(NUMEXP+1)= ZLL
	            ZVCALC(NUMEXP+1)= ZVV
	            FCION(NUMEXP+1) = (W4**2.)*(2.+IFL(KS))*(APART2**2) + (APART1**2.)
	            FCION(NUMEXP+1)= DSQRT(FCION(NUMEXP+1))
	            FCAUX(NUMEXP+1)= APART1
	            FCAUX2(NUMEXP+1)= APART2

	          ELSE
!
!	            Si el calculo fracasa hay que descontar uno al numero de observaciones al efecto de calcular la desviacion standar
	            NFAILURE = NFAILURE - 1
	            ZLCALC(NUMEXP+1)= 0.
	            ZVCALC(NUMEXP+1)= 0.
	            FCION(NUMEXP+1)= PENALIZ
	            FCAUX(NUMEXP+1)= 0.
	            FCAUX2(NUMEXP+1)= 0.

	          END IF
	          STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)**2
	          STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)**2

	        end if
!	        Nota: FIN ifl = -1 o ifl = -2
!
!
!	        IF(IFL(KS).EQ.0.OR.IFL(KS).EQ.1) then
	        IF(IFL(KS).GE.0) then
!
	          IF(IFL(KS).EQ.0) THEN   !  ifl=0
!
	            write(6,*) " i = ",i
!	            Inicio de la variable de error.
	            iErrTP = 0
	            !
!	            Este es el respectivo llamado al flash multifásico, en este caso binario. ZLIQUIDO se refiere
!	            siempre a la fase 1, que será líquida. ZVAPOR puede ser vapor o líquido, dependiendo del resul-
!	            tado del test de estabilidad. En el caso de ELL, se evalúa qué composición es mayor: X1 o Y1.
!	            El resultado del flash multifásico se ordena en base a eso.
!
	            feed(1) = (X1 + Y1)/2
	            feed(2) = 1D0-feed(1)
	            indF = 0
	            zFrac = 0
	            zVal = 0
!
!	            Organización de fases (en el caso de ELL)
	            orden = Y1 - X1
	            if(orden.LT.0)then

	              iOrden = 1

	            else

	              iOrden = 2

	            endif
	            call GC3Flash(NCw,NG,NST,TEMP,PRESS,feed(:15),zFrac,zVal,indF,IERRTP)
!
!	            if(minval(indF(:4)).GT.-1)then

!	              ELL: se procede a organizar las fases de la manera que se ingresaron:
	              if(zFrac(1,1).GT.zFrac(1,2))then

	                X1 = zFrac(1,iOrden)
	                zLiquido = zVal(iOrden)
	                Y1 = zFrac(1,3-iOrden)
	                zVapor = zVal(3-iOrden)

	              else

	                X1 = zFrac(1,3-iOrden)
	                zLiquido = zVal(3-iOrden)
	                Y1 = zFrac(1,iOrden)
	                zVapor = zVal(iOrden)

	              endif

!	            else

!	              ELV: se procede a asignar cada fase, como se ingresaron:
!	              minval(indF) = -1, porque a esta altura, debe haber haber vapor
!	              "minLoc" identifica en qué posición está el "-1".
!	              iVap(1) = minLoc(indF(:4),1)
!	              Y1 = zFrac(1,iVap(1))
!	              X1 = zFrac(1,3-iVap(1))

!	            endif
	            IF(IERRTP.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQUIDO
	              ZVCALC(NUMEXP+1) = ZVAPOR
	              XCALC(NUMEXP+1) = X1
	              YCALC(NUMEXP+1) = Y1
	              FCAUX(NUMEXP+1) = ABS(XCALC(NUMEXP+1)-XEXP(NUMEXP+1))

	            ELSE

	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              XCALC(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF

	          END IF   ! fin ifl = 0
!
!	          IF(IFL(KS).EQ.1) THEN   !  ifl=1
!
	          IF(IFL(KS).GE.1.AND.iFL(ks).LE.2) THEN   !  2>ifl>0

!
	            write(6,*) "i = ",i
!	            Inicio de la variable de error. Si no se inicia puede
!	            dar cualquier cosa.
	            iErrBurb = 0
!	            Antiguo llamado a Pburb:
!	            CALL ASSocTXFLASHSUB(2,WWDC,WWTEMPCRIT,NG,NUWW,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,  &
!&	                                 NGA,ONUASOCWW,MASSOC,ENASS,RKASS,TEMP,PRESS,X1,Y1,Z1,ZVAPVAP,IERRTX)


!	            La variable iType está relacionada con IFL.IFL = 1 implica Pburb, mientras que corresponde a
!	            iType = 3. Ahora existirá IFL = 2, que corresponde a Tburb.
	            if(IFL(KS).EQ.1)iType = 3
	            if(IFL(KS).EQ.2)iType = 1

	            call GCADB(1,NCw,NG,NST,iType,TEMP,PRESS,X1,Y1,iErrBurb,ZLiq,ZVap)
!
!	            IF(IERRTX.EQ.0) THEN	!ahora son cálculos de burbuja.
	            IF(IERRBurb.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQ
	              ZVCALC(NUMEXP+1) = ZVAP
	              YCALC(NUMEXP+1) = Y1
	              if(iType.EQ.3)then

	                PCALC(NUMEXP+1) = PRESS
	                FCAUX(NUMEXP+1) = ABS(PCALC(NUMEXP+1)/PEXP(NUMEXP+1)-1)

	              elseif(iType.EQ.1)then

	                Pcalc(numexp+1) = Temp !medio de los pelos, no?
	                FCAUX(NUMEXP+1) = ABS(PCALC(NUMEXP+1)/TEXP(NUMEXP+1)-1)

	              endif

	            ELSE
!
	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              PCALC(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF
	          END IF   ! fin ifl = 1

!
!	          Equilibrio ternario a TP. Si se estuvieran haciendo las cosas bien, esto iría
!	          todo junto con el flash binario a TP.
	          IF(IFL(KS).EQ.3) THEN   !  ifl=3
!
	            write(6,*) " i = ",i
!	            Inicio de la variable de error.
	            iErrTP = 0
	            !
!	            Este es el respectivo llamado al flash multifásico, en este caso ternario
!
	            feed(1) = (xExp(numExp+1) + yExp(numExp+1))/2
	            feed(2) = (xExp2(numExp+1) + yExp2(numExp+1))/2
	            feed(3) = 1D0-feed(1)-feed(2)
!
!	            Organización de fases (en el caso de ELL)
	            orden = yExp(numExp+1) - xExp(numExp+1)
	            if(orden.LT.0)then

	              iOrden = 1

	            else

	              iOrden = 2

	            endif
	            call GC3Flash(NCw,NG,NST,TEMP,PRESS,feed(:15),zFrac,zVal,indF,IERRTP)
!
	            if(minval(indF(:4)).GT.-1)then

!	              ELL: se procede a organizar las fases de la manera que se ingresaron:
	              if(zFrac(1,1).GT.zFrac(1,2))then

	                X1 = zFrac(1,iOrden)
	                X2 = zFrac(2,iOrden)
	                zLiquido = zVal(iOrden)
	                Y1 = zFrac(1,3-iOrden)
	                Y2 = zFrac(2,3-iOrden)
	                zVapor = zVal(3-iOrden)

	              else

	                X1 = zFrac(1,3-iOrden)
	                X2 = zFrac(2,3-iOrden)
	                zLiquido = zVal(3-iOrden)
	                Y1 = zFrac(1,iOrden)
	                Y2 = zFrac(2,iOrden)
	                zVapor = zVal(iOrden)

	              endif

	            else

!	              ELV: se procede a asignar cada fase, como se ingresaron:
!	              minval(indF) = -1, porque a esta altura, debe haber haber vapor
!	              "minLoc" identifica en qué posición está el "-1".
	              iVap(1) = minLoc(indF(:4),1)
	              Y1 = zFrac(1,iVap(1))
	              Y2 = zFrac(2,iVap(1))
	              X1 = zFrac(1,3-iVap(1))
	              X2 = zFrac(2,3-iVap(1))
	              zLiquido = zVal(3-iVap(1))
	              zVapor = zVal(iVap(1))

	            endif
	            IF(IERRTP.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQUIDO
	              ZVCALC(NUMEXP+1) = ZVAPOR
	              XCALC(NUMEXP+1) = X1
	              XCALC2(NUMEXP+1) = X2
	              YCALC(NUMEXP+1) = Y1
	              YCALC2(NUMEXP+1) = Y2
	              FCAUX(NUMEXP+1) = ABS(XCALC(NUMEXP+1)-XEXP(NUMEXP+1))

	            ELSE

	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              XCALC(NUMEXP+1) = 0.
	              XCALC2(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              YCALC2(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF

	          END IF   ! fin ifl = 0


!	          IF(IFL(KS).EQ.0.AND.IERRTP.EQ.0.OR.IFL(KS).EQ.1.AND.IERRTX.EQ.0.) THEN

!	            FCION(NUMEXP+1) = (W4**2)*(YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.D0)**2 + &
!&	                              IFL(KS)*(PCALC(NUMEXP+1)/PEXP(NUMEXP+1) - 1.D0)**2 + &
!&	                              (1-IFL(KS))*(XCALC(NUMEXP+1)/(XEXP(NUMEXP+1))  - 1.D0)**2
	          IF(IFL(KS).EQ.0.AND.IERRTP.EQ.0)then

!	            Para optimizar orgánicos+agua (sol baja de org en agua)
	            if(yExp(numExp+1) .GT. 0D0 .AND. xExp(numExp+1) .GT. 0D0)then

	              FCION(NUMEXP+1) = SQRT(w4**2*(YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.D0)**2 + &
&	                                     w5**2*(XCALC(NUMEXP+1)/(XEXP(NUMEXP+1))  - 1.D0)**2)
!
!	              Optimización general de LL:
!	              FCION(NUMEXP+1) = SQRT((YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.D0)**2 + &
!&	                                     (XCALC(NUMEXP+1)/(XEXP(NUMEXP+1))  - 1.D0)**2)

!	              FCION(NUMEXP+1)= DSQRT(FCION(NUMEXP+1))
	              FCAUX2(NUMEXP+1) = ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	            endif

	          elseif(IFL(KS).EQ.1.AND.IERRBurb.EQ.0.)THEN

	            if(yExp(numExp+1) .GT. 0D0)then

	              FCION(NUMEXP+1) = SQRT((W4**2)*(YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.D0)**2 + &
&	                                w2**2*fcAux(numexp+1)**2)
	              FCAUX2(NUMEXP+1) = ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	            endif

	          elseif(IFL(KS).EQ.2.AND.IERRBurb.EQ.0.)THEN

	            if(yExp(numExp+1) .GT. 0D0)then

	              FCION(NUMEXP+1) = SQRT((W4**2)*(YCALC(NUMEXP+1)/YEXP(NUMEXP+1) - 1.D0)**2 + &
&	                                w3**2*fcAux(numexp+1)**2)
	              FCAUX2(NUMEXP+1) = ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	            endif

	          elseif(iFL(KS).EQ.3.AND.iErrTP.EQ.0)then

	            if(yExp(numExp+1) .GT. 0D0)then

	              fcAux2(numExp+1) = abs(yCalc(numExp+1) - yExp(numExp+1))
	              fcion(numExp+1) = SQRT(w5**2*(FCAUX(NUMEXP+1)/xExp(numExp + 1))**2 &
&	                              + w6**2*(xCalc2(numExp+1)/xExp2(numExp+1) - 1)**2 &
&	                              + w4**2*(fcAux2(numExp+1)/yExp(numExp+1))**2 &
&	                              + w7**2*(yCalc2(numExp+1)/yExp2(numExp+1) - 1)**2)
	            endif

	          ELSE

	            FCION(NUMEXP+1)= PENALIZ
	            FCAUX2(NUMEXP+1) = 0.D0

	          END IF
	          STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)
	          STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)
	          if(mez(ks).EQ.3)then

	            STDVC(KS) = STDVC(KS) + ABS(XCALC2(NUMEXP+1) - XEXP2(NUMEXP+1))
	            STDVD(KS) = STDVD(KS) + abs(yCalc2(numExp+1) - yExp2(numExp+1))

	          endif
!
	        end if
!
	        NUMEXP = NUMEXP + 1

	      END DO

!           Desviaciones standar
	      IF(NFAILURE.NE.0) THEN

	        STDVA(KS) = STDVA(KS)/NFAILURE
	        STDVB(KS) = STDVB(KS)/NFAILURE
	        if(mez(ks).EQ.3)then

	          stDvC(ks) = stDvC(ks)/nFailure
	          stDvD(ks) = stDvD(ks)/nFailure

	        endif

	      ELSE

	        STDVA(KS) = 0.D0
	        STDVB(KS) = 0.D0

	      END IF
!
!
!	    Datos de propiedades en exceso
	    elseIF(iHen(ks).EQ.1)then

	      nFailure= 0
	      iErrEx = 0
	      stdva(ks) = 0D0
	      write(*,'(/,"Set Nro. ",I2)')ks
	      do i =1,nObs(ks)

	        Temp = Texp(numExp+1)
	        Press = Pexp(numexp+1)
	        feed(1) = xExp(numexp+1)
	        feed(2) = 1- feed(1)
	        write(*,'("  Pto. nro. ",I2)')i
	        call GCXS(ifl(ks),NG,NST,Temp,Press,feed(:15),prop_E_calc(numExp+1),zLiCalc(numexp+1),iErrEx)
	        if(iErrEx.EQ.0)then

	          if(prop_E_exp(numExp+1).NE.0D0)then

	            fcAux(numExp+1) = 1D0 - prop_E_calc(numExp+1)/prop_E_exp(numExp+1)

	          else

	            fcAux(numExp+1) = 0

	          endif
			  if(ifl(ks).EQ.1)then

	            fcion(numExp+1) = w11*fcAux(numExp+1)

			  else

			    fcion(numExp+1) = w12*fcAux(numExp+1)

			  endif

	        else

	          nFailure = nFailure + 1
	          zLiCalc = 0.
	          prop_E_calc(numExp) = 0.
	          fcAux(numExp+1) = 0.
	          fcion(numExp+1)= penaliz

	        endif
	        NUMEXP = NUMEXP + 1
	        if(i.EQ.1)numExp1 = numExp
	        nObsKS = nObs(ks)
	        if(i.EQ.nObsKS)numExp2 = numExp

	      enddo
	      if(nObs(ks)-nFailure.GT.0)stDvA(ks) = sqrt(sum(fcion(numExp1:numExp2)**2)/(nObs(ks)-nFailure))
!
!
!	    Datos de coeficientes de actividad
	    ELSEIF (iHen(ks).EQ.5) then

	      nfailure=nobs(ks)
	      stdva(ks) = 0D0
	      stdvb(ks) = 0D0
	      stdvc(ks) = 0D0
	      WRITE(6,*)"DADOS=", KS
	      do i=1,nobs(ks)

	        Temp  = texp(numexp+1)
	        Press = pexp(numexp+1)
	        x1 = xexp(numexp+1)
	        x2 = xexp2(numexp+1)
	        WRITE(6,*)"OBS=", NOBS(KS)
	        iaux=mez(ks)
!	        Inicio de iErrGam. Sino, el compilador puede hacer lo que se le cante.
	        iErrGam = 0

	        call rGamasub(iaux,wwdc, wwtempcrit,ng,nuww,qarea,tstr, gstr,gdot,gddot,rkij, rktij,alfaij,  &
&	                      nga, onuasocww, massoc,enass,rkass,temp, press, x1,x2,rgam1,rgam2,rgam3,zli,ierrgam)
	        IF(ierrgam==0)THEN

	          zlicalc(numexp+1) = zli
	          rgamcalc1(numexp+1) = rgam1
	          rgamcalc2(numexp+1) = rgam2
	          fcaux(numexp+1)=rgamcalc1(numexp+1)/rgamexp1(numexp+1)-1.D0
	          fcaux2(numexp+1)=rgamcalc2(numexp+1)/rgamexp2(numexp+1)-1.D0

	          fcion(numexp+1) = w8**2.*(rgamcalc1(numexp+1)/rgamexp1(numexp+1)-1.)**2. + &
&	                            w9**2*(rgamcalc2(numexp+1)/rgamexp2(numexp+1)-1.)**2.
	          If(mez(ks)==3)then

	            rgamcalc3(numexp+1) = rgam3
	            fcaux3(numexp+1)=rgamcalc3(numexp+1)/rgamexp3(numexp+1)-1.d0

	            fcion(numexp+1) = fcion(numexp+1)+w10*(rgamcalc3(numexp+1)/rgamexp3(numexp+1)-1.)**2.

	          end if
!
	          fcion(numexp+1) = SQRT(fcion(numexp+1))

	        else

	          nfailure = nfailure - 1
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

	        STDVA(KS) = 0.D0
	        STDVB(KS) = 0.0D0
	        STDVC(KS) = 0.

		END IF


	    END IF

	  END DO
!****************************Fin Datasets******************************

	END IF
!
	SUMQ = 0D0
	DO I=1,NCXP

	  SUMQ = SUMQ + FCION(I)**2

	END DO
!	sumQ = sum(fcion(:NCXP)**2)
!
	OPEN(UNIT=1,TYPE='OLD',file='paraest.out',FORM='FORMATTED',ACCESS='APPEND')
!

!	Descalo los valores antes de imprimirlos:
	CALL XXDESSCAL(NPAR,XX,IDXX,FGGG,FRKIJ,FTEMP,FQAREA,FALFA,FENASS,FRKASS,XXAUX,fTc,FDC)
	WRITE(1,'(<NPAR+1>(F12.4,X))') (XXAUX(I),I=1,NPAR),SUMQ

	CLOSE(UNIT=1)
	RETURN

ENDsubroutine
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
!	Estima DC solo con informacion de los puntos criticos de los compuestos puros.
!
SUBROUTINE INICDC(NC,PRCRIT,TEMPCRIT,DC)

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION DC(20),TEMPCRIT(20),PRCRIT(20)
!
	R=82.05	!atm cm3/(mol K)
	DO I=1,NC
!
!******************************
!	  if(DC(I).EQ.0)then ! Fran. 26/05/09
!******************************
	  DC(I) = 0.4471866767*((R*TEMPCRIT(I)/PRCRIT(I)))**(1./3.)
!******************************
!	  end if              ! Fran. 26/05/09
!*****************************
!
	END DO
!
	return
ENDsubroutine
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
SUBROUTINE FORSCAL(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,numTc,NUMDC)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	parameter(NCM=20)
	logical inIJ(8)
	common/IJ/inIJ
!
!	Identificacion (indices) de param gii o kij a ajustar.
	DIMENSION IDG(8,2),ITEMP(8),IDXX(8),xx(8)
!
	inIJ(:8) = .FALSE.
	NPAR = 0
	IF(NNP.NE.0) THEN

	  DO I=1,NNP

	    IF(IDG(I,1).EQ.IDG(I,2)) THEN

	      IDXX(NPAR+1)=1
	      NPAR = NPAR + 1

	    ELSE

	      IDXX(NPAR+1)=2
	      NPAR = NPAR + 1
	      inIJ(nPar) = .TRUE.

	    END IF
	    IF(ITEMP(I).EQ.1)then

	      IDXX(NPAR+1) = 3
	      NPAR = NPAR + 1
	      IF(IDG(I,1).NE.IDG(I,2))inIJ(nPar) = .TRUE.
!
	    elseIF(ITEMP(I).EQ.2)then

	      IDXX(NPAR+1) = 3
!
!	      Ahora hay un tipo distinto de escalado para g' y g":
	      idXX(nPar+2) = 10
	      NPAR = NPAR + 2

	    endif
!
	  END DO

	END IF
	IF(NNQ.NE.0) THEN

	  DO I=1, NNQ

	    IDXX(NPAR+1)= 4
	    NPAR = NPAR + 1

	  eND DO

	END IF
	IF(NALFA.NE.0) THEN

	  DO I=1, NALFA

	    IDXX(NPAR+1)= 5
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(NENASS.NE.0) THEN

	  DO I=1,NENASS

	    IDXX(NPAR+1)= 6
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(NRKASS.NE.0) THEN

	  DO I=1,NRKASS

	    IDXX(NPAR+1)= 7
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(numTc.NE.0) THEN

	  DO I=1,numtc

	    IDXX(NPAR+1)= 8
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(NUMDC.NE.0) THEN

	  DO I=1,NUMDC

	    IDXX(NPAR+1)= 9
	    NPAR = NPAR + 1

	  END DO

	END IF
	return
ENDsubroutine !FIN FORSCAL
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
SUBROUTINE XXSCAL(NPAR,XX,IDXX)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	parameter(NCM=20)
	parameter(pi=3.141592654D0)
	DIMENSION XX(8), IDXX(8)
	logical inIJ(8)
!
	COMMON/ESCALADO/FGGG,FRKIJ,FTEMP,FQAREA,FALFA,FENASS,FRKASS,ftc,FDC
	common/escMax/vMaxAlfa,vMaxkp,vMaxKstr,vMinGp,vMaxGp,error_Dc(NCM),dc_medio(NCM)
	common/IJ/inIJ
!
!	Transforma los valores de xx segun los siguientes factores de escalado
	FGGG  = 1.D-10
	FRKIJ = 1.D-7
	FTEMP = 1.D-7
	FQAREA = 1.D-0
	FALFA = 1.D-7
	FENASS = 1.D-7    !TENTATIVO
	FRKASS = 1.D-8     !TENTATIVO
	fTc = 1.D-8	       !TENTATIVO
	FDC = 1.D-8        !TENTATIVO

!	FGGG  = 1.D-10
!	FRKIJ = 1.D-5
!	FTEMP = 1.D-6
!	FQAREA = 1.D-5
!	FALFA = 1.D-6
!	FENASS = 1.D-9     !TENTATIVO
!	FRKASS = 1.D-4     !TENTATIVO
!	FDC = 1.D-6        !TENTATIVO
!
!	Valor máximo kij* (0 y vMaxKstr):
	vMaxKstr = 2D0
!	Valor máximo para los kij' (entre -vMaxKp y vMaxKp):
	vMaxkp = 0.151D0
!	Valor máximo para los alfas (entre -vMaxAlfa y vMaxAlfa):
	vMaxAlfa = 135D0
!	Valores mínimo y máximo de g'
	vMaxGp = 2.D0
	vMinGp = 1D-3
!
!
	DO I=1,NPAR
!
!	  Escalado de g*
	  IF(IDXX(I).EQ.1)then

	    XX(I) = FGGG*XX(I)
!
!	  Escalado de k*
	  elseIF(IDXX(I).EQ.2)then

!	      DTan(x) es la "T"angente de "x" con "D"oble precisión.
	      xx(i) = dtan(pi*(xx(i)/vMaxKstr - 1D0)/2D0)*frKij

!	  Escalado de g' y k':
	  elseIF(IDXX(I).EQ.3)then

	    if(.NOT.inIJ(i))then

!	      g' es siempre < -0.001
!	      XX(I) = FTEMP*dlog(-(XX(I)+ vMinGp))
	      xx(i) = fTemp*dtan(pi*((xx(i) + vMinGp)/vMaxGp - 1D0)/2D0)

	    else

!	      Sí para los k':
	      xx(i) = dtan(xx(i)*pi/2D0/vMaxkp)*fTemp

	    endif

	  elseif(idXX(i).EQ.10)then

!	    De momento, g" no tiene ninguna limitación...
	    xx(i) = fTemp*xx(i)

	  elseIF(IDXX(I).EQ.4)then

	    XX(I) = FQAREA * XX(I)

	  elseif(idXX(i).EQ.5)then

!	    Escalado con límites para alfas:
	    xx(i) = dtan(xx(i)*pi/2D0/vMaxAlfa)*fAlfa

	  elseIF(IDXX(I).EQ.6)then

	    XX(I) = FENASS*XX(I)

	  elseIF(IDXX(I).EQ.7)then

	    XX(I) = FRKASS*XX(I)

	  elseIF(IDXX(I).EQ.8)then

	    XX(I) = fTc*XX(I)

	  elseIF(IDXX(I).EQ.9)then

!	    XX(I) = FDC*XX(I)
!
!	    Valores máximo y mínimo para dc:
	    xx(i) = fDc*dtan(pi*(xx(i) - dc_medio(i))/2.D0/error_Dc(i))

	  endif

	END DO
!
	return
!
ENDsubroutine !FIN XXSCAL
!
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!	Antitransforma los valores de xx segun los factores de escalado establecidos en la rutina XXSCAL
!
SUBROUTINE XXDESSCAL(NPAR,XX,IDXX,FGGG,FRKIJ,FTEMP,FQAREA,FALFA,FENASS,FRKASS,XXAUX,fTc,FDC)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XX(8), IDXX(8), XXAUX(8)
	parameter(pi=3.141592654D0,NCM=20)
	logical inIJ(8)
!
	common/escMax/vMaxAlfa,vMaxKp,vMaxKstr,vMinGp,vMaxGp,error_Dc(NCM),dc_medio(NCM)
	common/IJ/inIJ
!
	DO I=1,NPAR

	  XXAUX(I) = XX(I)
	  IF(IDXX(I).EQ.1)then

	    XXAUX(I) = XXAUX(I)/FGGG

	  elseIF(IDXX(I).EQ.2)then

!	    "DATan(x)" es el "A"rco "T"angente de "x" con "D"oble precisión.
	    xxAux(i) = vMaxKstr*(1D0 + daTan(xxAux(i)/frKij)*2D0/pi)

	  elseIF(IDXX(I).EQ.3)then

	    if(.NOT.inIJ(i))then

!	      XXAUX(I) = -dexp(XXAUX(I)/FTEMP) - vMinGp
	      xxAux(i) = vMaxGp*(2D0*daTan(xx(i)/fTemp)/pi - 1D0) - vMinGp

	    else

	      xxAux(i) = datan(xxAux(i)/fTemp)*vMaxKp*2D0/pi

	    endif

	  elseif(idXX(i).EQ.10)then

	    xxAux(i) = xxAux(i)/fTemp

	  elseIF(IDXX(I).EQ.4)then

	    XXAUX(I) = XXAUX(I)/FQAREA


	  elseif(idXX(i).EQ.5)then

	    xxAux(i) = datan(xxAux(i)/fAlfa)*vMaxAlfa*2D0/pi

	  endif
	  IF(IDXX(I).EQ.6) XXAUX(I) = XXAUX(I)/FENASS
	  IF(IDXX(I).EQ.7) XXAUX(I) = XXAUX(I)/FRKASS
	  IF(IDXX(I).EQ.8) XXAUX(I) = XXAUX(I)/fTc
	  IF(IDXX(I).EQ.9)then

	    XXAUX(I) = dc_medio(i) + datan(xxAux(i)/fDc)*error_Dc(i)*2D0/pi

	  endif

	END DO
!
	return
ENDsubroutine !FIN XXDESSCAL
!
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!
!
SUBROUTINE WRITESTAT(SSQ,XJAC,WORK,INFER,IER)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XJAC(260,8), WORK(436)
!
	OPEN(UNIT=1,TYPE='OLD',file='paraest.out',FORM='FORMATTED',ACCESS='APPEND')
!
	WRITE(1,95)
95	FORMAT(/,X,' Valor final fcion objetivo = SSQ = ',/)
	WRITE(1,*) SSQ
!
	WRITE(1,97) IER,INFER
 97	FORMAT(/,' (Ver manual de la rutina ZXSSQ (Lib IMSL) ',/,X,' IER = ',I5,/,X,' INFER = ',I5)
!
	WRITE(1,99) WORK(1),WORK(2),WORK(3),WORK(5)
99	FORMAT(/,'  Norma del gradiente = ',D12.5,/,X,' Numero de evaluaciones de la Fc Obj. = ',F10.2,/, &
&	       X,' Numero de digitos significativos en los param. escalados = ',F10.2,/,X,' Numero de iteraciones = ',F10.2)
	CLOSE(UNIT=1)
	return

END !FIN WRITESTAT
!------------------------------------------------------------------------------------------------------------------------------------
!
!
!	Esta subrutina escribe los resultados de evaluación de los datos experimentales alimentados frente a la predicción del
!	modelo GCA-EoS, sea con los parámetros ingresados o el resultado final del ajuste.
!	Escribe datos archivo de salida paraest.out
!
SUBROUTINE WRITESAL2(NC,NG,NU,QAREA,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC,ENASS,RKASS,&

&	               ISAT,IGLE,IVE,IUNI,TSAT,PSAT,NSAT,IDC,NOBS,IHEN,mez,IFL,NTEXT,                       &

&	               TEXP,PEXP,XEXP,xexp2,YEXP,yExp2,IDSAT,DC,tempcrit,XCALC,YCALC,xCalc2,yCalc2,ZLCALC,ZVCalc,FCAUX,PCALC,  &

&                    IOBJSAT,FCAUX2,fcaux3,STDVA,STDVB,stdvc,stDvD,rgamcalc1, rgamcalc2,rgamcalc3,numTc,idTc,   &
&	               NUMDC,IDDC,TSPLINE,rgamexp1, rgamexp2,rgamexp3,zLicalc,prop_E_calc,prop_E_exp)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	parameter(NDM=260)
!
!	Dim. dispersion
	DIMENSION NU(10,20), QAREA(10), TSTR(10), GSTR(10)
	DIMENSION GDOT(10), GDDOT(10), RKIJ(10,10), RKTIJ(10,10)
	DIMENSION ALFAIJ(10,10)
	DIMENSION TSPLINE(10)
!
!	Dim. asociacion
      DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
      DIMENSION ONUASOC(12,20)
!
!	Dim diametro de esferas rigidas en el punto critico
	DIMENSION DC(20)
!
	DIMENSION IDSAT(20)
!
	DIMENSION  TSAT(20), PSAT(20),tempcrit(20)
	DIMENSION TEXP(260), PEXP(260), XEXP(260), YEXP(260), xexp2(260),yExp2(NDM)
	dimension rgamexp1(260),rgamexp2(260),rgamexp3(260)
!
!	Dimensionamiento  nro. de datos de pvap comp puros
	DIMENSION NSAT(20)
	DIMENSION IOBJSAT(20)
!
!	Dimensionamientos que tienen que ver con los datasets
	DIMENSION IDC(25,3),NOBS(25),IHEN(25),mez(25)
	DIMENSION NTEXT(25,65)
	DIMENSION IFL(25)
!
	DIMENSION PCALC(260), XCALC(260), YCALC(260),xCalc2(260),yCalc2(260)
	DIMENSION  ZLCALC(260), ZVCALC(260)
	dimension rgamcalc1(260),rgamcalc2(260),rgamcalc3(260)
 	dimension zlicalc(260)
	dimension prop_E_calc(NDM),prop_E_exp(NDM)
!
	DIMENSION FCAUX(260)
	DIMENSION FCAUX2(260)
	DIMENSION FCAUX3(260)
!
!	Dimensionamiento de la desviacion standar para las comp puros
	DIMENSION STDVSAT(20)
!
!	Dimensionamiento desviaciones standar para mezclas
	DIMENSION STDVA(25),STDVB(25),STDVC(25),stDvD(25)
!
!	Dimensionamiento arreglo de identificacion de param. dc considerados como simples variables de optimizcion
	DIMENSION IDDC(8),idTc(8)
!
	OPEN(UNIT=1,TYPE='OLD',file='paraest.out',FORM='FORMATTED',ACCESS='APPEND')
!
!
	WRITE(1,47)
 47	FORMAT(/,X,'--------------------------------------------------',/, &
&	       X,'------------ PARAMETROS FINALES ----------------',/,     &
&	       X,'--------------------------------------------------',/,   &
&	       X,' PARAMETROS DE GRUPOS PUROS (termino dispersivo)',/)
!
	WRITE(1,3) (I, I=1,NC)
3	FORMAT(40X,'COMPONENTES',/,X,' Gr.Nro.:  T*        Q        TSPLINE ',<NC>I3,/)
!
	DO I=1,NG

	  WRITE(1,4) I,TSTR(I),QAREA(I),TSPLINE(I), (NU(I,J), J=1,NC)
4	  FORMAT(X,I3,5X,2F10.4,X,E10.4,<NC>I3)

	END DO
!
	WRITE(1,5)
5	FORMAT(/,X,' Energias de interacción gjj ',/,X,'(la forma funcional cambia a altas temperaturas)',/)
!
	DO I=1,NG
!
	  WRITE(1,6) I,I,GSTR(I), GDOT(I), GDDOT(I)
6	  FORMAT(X,' g(',I3,I3,') = ',F10.1,' * [ 1 + ',F10.4,' *(T/T* - 1) + ',F10.4, ' *LN(T/T*)] ',/)

	END DO

	WRITE(1,7)
 7	FORMAT(/,X,' Parametros binarios: [k*] ',/)

	DO I=1,NG

	  WRITE(1,8) (RKIJ(I,J), J=1,NG)
8	  FORMAT(X,<NG>F10.4)

	END DO
!
	WRITE(1,9)
9	FORMAT(/,X,' Dependencia resp. a la temp. : [k`]',/)
	DO I=1,NG

	  WRITE(1,8) (RKTIJ(I,J), J=1,NG)

	END DO
!
	WRITE(1,10)
 10	FORMAT(/,X,' Parametros binarios de no azar: [ALFA]',/)
	DO I=1,NG

	  WRITE(1,8) (ALFAIJ(I,J), J=1,NG)

	END DO
	WRITE(1,186) NGA
 186	FORMAT(/,X,' Numero de grupos considerados en la contribucion de asociacion = ',I3)
!
	IF(NGA.GT.0) THEN

	  WRITE(1,11) (I, I=1,NC)
11	  FORMAT(/,X,' PARAMETROS DE ASOCIACION ',/,23X,'COMPONENTES',/,' GRUPO 	  # SIT ASOC  ',<NC>I3,/)

	  DO I=1,NGA

	    WRITE(1,12) I, MASSOC(I), (ONUASOC(I,J), J=1,NC)
12	    FORMAT(3X,I3,9X,I3,6X,<NC>f5.2)

	  END DO
	  WRITE(1,13)
13 	  FORMAT(/,X,'Parametros de energia de asociacion (K)  ',/)
!
!
!	  Escritura de parametros energeticos de asociacion
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then

	      DO L=1,NGA

              IF(MASSOC(L).NE.0)then

	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,14) I,J,K,L,ENASS(I,J,K,L)
14	              FORMAT(X,' ENASS(',I2,I2,X,I2,I2,')  = ',F10.4)

	            END DO

	          END DO

	        endif

	      END DO

	    endif

	  END DO
!
!	  Escritura de parametros volumetricos de asociacion (cm3/mol)
	  WRITE(1,15)
15	  FORMAT(/,X,'Parametros de volumen de asociacion (cm3/mol) ',/)
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then

	      DO L=1,NGA

	        IF(MASSOC(L).NE.0)then

	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,16) I,J,K,L,RKASS(I,J,K,L)
16	              FORMAT(X,' RKASS(',I2,I2,X,I2,I2,') = ',F16.10)

	            END DO

	          END DO

	        endif

	      END DO

	    endif

	  END DO

	endif
!
!
!	Escritura del resultado del ajuste de dc
	DO I=1,NC

	  DO JAMO=1,NUMDC

	    IF(IDDC(JAMO).EQ.I)  THEN

	      WRITE(1,1548) I,DC(I)
1548	      FORMAT(/,X,' El diametro de esferas rigidas en el punto critico del componente Nro.',I3,/,&
&	             ' fue considerado otra variable de optimizacion ',/,X,' DC = ', F10.7)
	      GO TO 505

	    END IF

	  END DO
!
	  IF(TSAT(I).NE.0) WRITE(1,500) I,TSAT(I),PSAT(I),DC(I)
500	  FORMAT(/,X,' El diametro de esferas rigidas en el punto critico del componente Nro.',I3,' fue reajustado para ',/, &
&	         X,' reproducir el punto de la curva de presion de vapor ',/,x,' [ ',F10.4,' K ;',F10.4,' atm ] ',/,X,' DC = ', F10.7)
505	END DO
 	DO I=1,NC
!
	  DO mm=1,numtc

	    IF(idtc(mm).EQ.I)  THEN

	      WRITE(1,1549) I,tempcrit(I)
1549	      FORMAT(/,X,'La temperatura critica del componente Nro.',I3,/,X,'fue considerado otra variable de optimizacion ',/,'Tc= ', F8.2," K")

	    END IF

	  END DO

	END DO
!
!
!	****************Escritura del resultado del ajuste de datos de compuestos puros*****************
!
!	Cuenta el numero de datos de equilibrio
	NCXP = 0
 2	FORMAT(20I3)
!
!	Nueva unidad de escritura, sólo para graficar con facilidad. No tiene intención de reemplazar
!	al archivo "paraest.out". La extensión "CSV" es "texto delimitado por comas".
	open(unit=10,file="graficos.csv",form="FORMATTED")
!
	IF(ISAT.NE.0) THEN

	  WRITE(1,338)
	  write(10,338)
338	  FORMAT(/,X,'***** Resultados: compuestos puros *****',/)
	  DO I=1,NC

	    IF(NSAT(I).NE.0) THEN

	      IF(IOBJSAT(I).EQ.50) THEN

	        WRITE(1,393)
		  write(10.393)
393	        FORMAT(/,X,'  Comp   Texp/K     Pexp/atm         Pexp/atm          Zliq       Zvap  ERR,REL,%,ISOF',/)

		END IF
	      IF(IOBJSAT(I).EQ.80) THEN

	        WRITE(1,333)
		  write(10,333)
333	        FORMAT(/,X,'  Comp   Texp/K     Pexp/atm        PcalC/atm        Zliq        Zvap  ERR,REL,%,P ',/)

	      END IF
!           De paso calculo desviacion standar para cada componente.
	      NFOLI=NSAT(I)
	      STDVSAT(I)=0D0
	      DO J=1,NSAT(I)

	        WRITE(1,334) I,TEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100
	        WRITE(10,334) I,TEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100
334	        FORMAT(X,I3,4X,F10.4,X,F16.10,X,F16.10,X,F11.8,X,F10.7,X,F10.4)
		  STDVSAT(I) = STDVSAT(I) + FCAUX(NCXP+1)**2
	        if(FCAUX(NCXP+1).EQ.0) NFOLI = NFOLI - 1
	        NCXP = NCXP + 1

	      END DO
	      IF(NFOLI.NE.0) THEN

	        STDVSAT(I) = STDVSAT(I) / NFOLI

	      ELSE

	        STDVSAT(I) = 0.D0

	      END IF
	      STDVSAT(I) = 100*DSQRT(STDVSAT(I))
	      WRITE(1,335) STDVSAT(I)
		write(10,335)STDVSAT(I)
	      WRITE(*,335) STDVSAT(I)
335	      FORMAT(/,X,' Desviación estándar %',41X,f10.4,/)

	    END IF

	  END DO

	END IF
!
!	Escritura del resultado del ajuste de datos sobre mezclas binarias
	IF(IGLE.NE.0) THEN

	  WRITE(1,340)
	  write(10,340)
340	  FORMAT(/,X,'********** RESULTADOS: MEZCLAS BINARIAS *******************',/)
	  DO KS=1,IGLE

	    WRITE(1,339) (NTEXT(KS,J),J=1,65)
	    WRITE(10,339) (NTEXT(KS,J),J=1,65)
339	    FORMAT(//,80A2)
!
!	    Dimensión del número entre paréntesis al imprimir x(i). Si hay espacios, Excel no separa
!	    Correctamente la cabecera de las columnas:
	    if(idC(ks,1).LT.10)then

	      iDim = 1

	    else

	      iDim = 2

	    endif
	    if(idC(ks,2).LT.10)then

	      iDim2 = 1

	    else

	      iDim2 = 2

	    endif
	    IF(IHEN(KS).EQ.0) THEN

!	      ********************Equilibrio binario***********************
!
!
	      IF(IFL(KS).EQ.-2) WRITE(1,343) IDC(KS,2)
	      IF(IFL(KS).EQ.-2) WRITE(10,343) IDC(KS,2)
343	      FORMAT(X,/,' NOTA: El apartamiento con respecto a la isofugacidad del componente ',I3,/, &
&	             X,' no fue utilizado en la evaluacion de la funcion objetivo,aunque',/,           &
&	             X,' igualmente se informa ',/)
	      IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then

	        WRITE(1,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2)
	        WRITE(10,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2)
341	        FORMAT(X,'  Temp/K    Pres/ATM  Xexp(',I<iDim>,')  Yexp(',I<iDim>,') Fug(',I<iDim>,')% Fug(',I<iDim>,')%   Zliq      Zvap ',/)
	        DO I=1,NOBS(KS)

	          WRITE(1,345) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),YEXP(NCXP+1),FCAUX(NCXP+1),FCAUX2(NCXP+1),ZLCALC(NCXP+1), ZVCALC(NCXP+1)
	          WRITE(10,345) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),YEXP(NCXP+1),FCAUX(NCXP+1),FCAUX2(NCXP+1),ZLCALC(NCXP+1), ZVCALC(NCXP+1)
345	          FORMAT(6F10.5,2F9.6)
	          NCXP = NCXP + 1

	        END DO
	        WRITE(1,347) STDVA(KS),STDVB(KS)
	        WRITE(6,347) STDVA(KS),STDVB(KS)
	        WRITE(10,347) STDVA(KS),STDVB(KS)
347	        FORMAT(/,X,' Desv.Standar %',24X,2f10.4)

            elseIF(IFL(KS).EQ.0) then

!	        Resultados de flash a T y P:
	        WRITE(1,351) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,1)
351	        FORMAT(/,X,'  Temp/K   Pres/ATM  Xexp(',I<iDim>,')  Xclc(',I<iDim>,') Err%    Yexp(',I<iDim>,') Ycalc(',I<iDim>,')     Err%',/)
!
	        WRITE(10,3511) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,1)
3511	        FORMAT(/,X," (x-y)(",I<iDim>,")exp,(x-y)(",I<iDim>,")calc,T/K,P/atm",/)
!
!	        Recontador para "gráficos"
	        NCXPy = NCXP
	        DO I=1,NOBS(KS)
!
	          WRITE(1,349) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),XCALC(NCXP+1), FCAUX(NCXP+1),YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&	                       ZLCALC(NCXP+1), ZVCALC(NCXP+1)
349	          FORMAT(8(F10.5,X),2X,' Zliq = ',F10.8,'       Zvapor =  ',F10.8)
!
!	          "x" en "gráficos.cvs":
	          write(10,3490)XEXP(NCXP+1),XCALC(NCXP+1),TEXP(NCXP+1),PEXP(NCXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
!	        "y" en "gráficos.cvs":
	        do i = 1,Nobs(ks)

	          write(10,3490)YEXP(NCXPy+1),YCALC(NCXPy+1),TEXP(NCXPy+1),PEXP(NCXPy+1)
	          NCXPy = NCXPy + 1

	        enddo
3490	        format(2(F9.6,","),F7.2,",",F10.4)
	        WRITE(1,353) STDVA(KS),STDVB(KS)
	        WRITE(6,353) STDVA(KS),STDVB(KS)
	        WRITE(10,353) STDVA(KS),STDVB(KS)
353	        FORMAT(/,X,' AAD(x) =',f10.4,20X,"AAD(y) = ",f10.4)

	      elseIF(IFL(KS).EQ.1) then

!	        Resultados de flash a T y x:
	        NCXPy = NCXP
	        WRITE(1,355)IDC(KS,1), IDC(KS,1), IDC(KS,1)
355	        FORMAT(/,X,'  Temp/K      Xexp(',I<iDim>,')    Pres/ATM    PREScalc      ERR%       YEXP(',I<iDim>,')    YCALC(',I<iDim>,')    ERR ',/)
	        write(10,3550) IDC(KS,1), IDC(KS,1)
3550	        format(/,"T /K,(x-y)(",I<iDim>,")exp,Pexp /atm,(xexp-ycalc)(",I<iDim>,"),Pcalc /atm")
	        DO I=1,NOBS(KS)

	          WRITE(1,357) TEXP(NCXP+1),XEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), FCAUX(NCXP+1),YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&	                       ZLCALC(NCXP+1), ZVCALC(NCXP+1)
357	          FORMAT(F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,'  ZLIQ = ',F10.8,'      ZVAP = ',F10.8)
!
!	          "x" en "gráficos.cvs":
	          write(10,3570)TEXP(NCXP+1),XEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),PCALC(NCXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
	        do i = 1,NObs(ks)

	          write(10,3570)TEXP(NCXPy+1),yEXP(NCXPy+1),PEXP(NCXPy+1),ycalc(NCXPy+1),PCALC(NCXPy+1)
	          NCXPy = NCXPy + 1

	        enddo
	        WRITE(1,359) 100*STDVA(KS),STDVB(KS)
	        WRITE(6,359) 100*STDVA(KS),STDVB(KS)
	        WRITE(10,359) 100*STDVA(KS),STDVB(KS)
3570	        format(F7.2,",",F9.6,",",F10.4,",",F9.6,",",F10.4)
359	        FORMAT(/,X,' ARD(P) =',F8.2," %",20X,"AAD(y) = ",f10.4)

	      elseIF(IFL(KS).EQ.2) then

!	        Resultados de flash a P y x:
	        NCXPy = NCXP
	        WRITE(1,356) IDC(KS,1),IDC(KS,1),IDC(KS,1)
356	        FORMAT(/,X,'  Pres/ATM    Xexp(',I<iDim>,')    Temp/K     TempCalc       ERR%       Yexp(',I<iDim>,')     Ycalc(',I<iDim>,')     ERR ',/)
	        write(10,'(/,"P /atm,(x-y)(",I<iDim>,")exp,T /K,(xexp-ycalc)(",I<iDim>,"),Tcalc /K")')IDC(KS,1),IDC(KS,1)
	        DO I=1,NOBS(KS)

	          WRITE(1,358) PEXP(NCXP+1),XEXP(NCXP+1),TEXP(NCXP+1),PCALC(NCXP+1), FCAUX(NCXP+1),YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&	                       ZLCALC(NCXP+1), ZVCALC(NCXP+1)
358	          FORMAT(F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,'  ZLIQ = ',F10.8,'      ZVAP = ',F10.8)
!
!	          "x" en "gráficos.cvs":
	          write(10,3580)PEXP(NCXP+1),XEXP(NCXP+1),TEXP(NCXP+1),XEXP(NCXP+1),PCALC(NCXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
	        do i = 1,NObs(ks)

!	          "y" en "gráficos.cvs":
	          write(10,3580)PEXP(NCXPy+1),yEXP(NCXPy+1),TEXP(NCXPy+1),ycalc(NCXPy+1),PCALC(NCXPy+1)
	          NCXPy = NCXPy + 1

	        enddo
3580	        format(F10.4,",",F9.6,",",F7.2,",",F9.6,",",F7.2)
	        WRITE(1,360) 100*STDVA(KS),STDVB(KS)
	        WRITE(6,360) 100*STDVA(KS),STDVB(KS)
	        WRITE(10,360) 100*STDVA(KS),STDVB(KS)
360	        FORMAT(/,X,' ARD(T) =',F8.2," %",20X,"AAD(y) = ",f10.4)

	      elseIF(iFL(KS).EQ.3)then

!	        Resultados de equilibrio ternario
	        WRITE(1,451) IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
451	        FORMAT(/,X,"  T/K     P/ATM   xExp(",I<iDim>,',I) xClc(',I<iDim>,',I)    ERR    xExp(',I<iDim2>,',I) xClc(',I<iDim2>,',I)    ERR    xExp(',I<iDim>,',II)  xClc(',I<iDim>,',II)      ERR     xExp(',I<iDim2>,',II)  xClc(',I<iDim2>,',II)     ERR%',/)
!
	        WRITE(10,4511) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,2)
4511	        FORMAT(/,X,"xExp(",I<iDim>,"),xExp(",I<iDim>,"),xClc(",I<iDim>,"),xClc(",I<iDim>,"),T /K,P /atm")
!
!	        Recontador para "gráficos"
	        NCXPy = NCXP
	        DO I=1,NOBS(KS)
!
	          er2 = abs(XEXP2(NCXP+1) - XCALC2(NCXP+1))
	          er4 = abs(YEXP2(NCXP+1) - YCALC2(NCXP+1))
	          WRITE(1,449) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),XCALC(NCXP+1),fcAux(ncXp+1),XEXP2(NCXP+1),XCALC2(NCXP+1),er2,YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&	                       YEXP2(NCXP+1),YCALC2(NCXP+1),ZLCALC(NCXP+1), ZVCALC(NCXP+1)
449	          FORMAT(F7.2,7F10.5,6(X,F10.5,X),2X,' ZI = ',F10.8,'       ZII =  ',F10.8)
!
!	          "x" en "gráficos.cvs":
	          write(10,4490)XEXP(NCXP+1),xExp2(NCXP+1),xCalc(ncXP+1),xCalc2(ncXP+1),TEXP(NCXP+1),PEXP(NCXP+1),yExp(ncXp+1),yExp2(ncXp+1),yCalc(ncXP+1),yCalc2(ncXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
4490	        format(4(F9.6,","),F7.2,",",F10.4,/,4(F9.6,","),/)
	        WRITE(1,453) STDVA(KS),STDVB(KS),stDvC(ks),stDvD(ks)
	        WRITE(6,453) STDVA(KS),STDVB(KS),stDvC(ks),stDvD(ks)
	        WRITE(10,453) STDVA(KS),STDVB(KS),stDvC(ks),stDvD(ks)
453	        FORMAT(/,X,' AAD(x) =',10X,f10.4,20X,f10.4,25X,f10.4,25X,F10.4)

	      END IF

	    END IF
!
	    IF (IHEN(KS).EQ.1)THEN

	      if(ifl(ks).EQ.1)then

!	        Resultados de volumen en exceso:

	        write(1,'(/,"  T(K)       P(atm)    Xexp(",I<iDim>,")  Vexp^E(cm3/mol)   Vcalc^E(cm3/mol)       Z",/)')idC(ks,1)
		  write(10,'(/,"Temp (K),P (atm),Xexp(",I<iDim>,"),Vexp^E (cm3/mol),Vcalc^E (cm3/mol),Z",/)')idC(ks,1)

	      elseif(ifl(ks).EQ.2)then

!	        Resultados de entalpía en exceso:
	        write(1,'(/,"  T(K)       P(atm)    Xexp(",I<iDim>,")  Hexp^E(cm3/mol)   Hcalc^E(cm3/mol)       Z",/)')idC(ks,1)
	        write(10,'(/,"T (K),P (atm),Xexp(",I<iDim>,"),Hexp^E (J/mol),Hcalc^E (J/mol),Z",/)')idC(ks,1)

	      endif
	      do i = 1,nObs(ks)

	        write(1,'(F7.2,2X,F10.5,2X,F10.8,X,2(F10.3,8X),F11.8)')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),prop_E_exp(NCXP+1),prop_E_calc(NCXP+1),zLiCalc(NCXP+1)
	        write(10,'(F7.2,",",F10.5,",",F10.8,",",2(F10.2,","),F11.8)')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),prop_E_exp(NCXP+1),prop_E_calc(NCXP+1),zLiCalc(NCXP+1)
	        NCXP = NCXP + 1

	      enddo
	      write(1,5000)stDvA(ks)*100.
	      write(10,5000)stDvA(ks)*100.
5000	      format(/,"  Desv.Std.:",F10.2," %")

	    endif
	    IF (IHEN(KS).EQ.5)THEN
!
!	      Resultados de coeficientes de actividad:
	      If(mez(ks)==2)then

	        write(1,1344)IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
	        DO I=1,NOBS(KS)

	          WRITE(1,1345)TEXP(NCXP+1),XEXP(NCXP+1),rgamexp1(ncxp+1),rgamcalc1(NCXP+1),rgamexp2(ncxp+1),rgamcalc2(NCXP+1), zlicalc(NCXP+1)
1344	          FORMAT(X,'Temp/K       Xex(',I<iDim>,')   Gex(',I<iDim>,')   Gc(',I<iDim>,')     Gex(',I<iDim2>,')   Gc(',I<iDim2>,')'       ,'Zliq')
1345	          format(X, 2f10.2, 4f10.3,f10.3)
1347	          FORMAT(2F8.4,f7.1,2F8.4,f6.3,f8.4)
 	          NCXP = NCXP + 1

	        END DO
	        WRITE(1,1379) STDVA(KS),STDVB(KS)
	        WRITE(6,1379) STDVA(KS),STDVB(KS)
1379	        FORMAT(/,X,' Desv.Std. %',2X,f10.4,7X,f10.4,/)

	      elseif(mez(ks)==3)then

	        write(1,1350) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,3),IDC(KS,3)
	        DO I=1,NOBS(KS)

	          wRITE(1,1348)TEXP(NCXP+1),XEXP(NCXP+1),xexp2(ncxp+1),rgamexp1(NCXP+1),rgamcalc1(NCXP+1),rgamexp2(NCXP+1),rgamcalc2(NCXP+1), &
&	                       rgamexp3(NCXP+1),rgamcalc3(NCXP+1)
1348	          FORMAT(X, f6.2 , 2(1X,e10.4), 6(1X,E12.3))
1349	          FORMAT(3F9.6)
1350	          FORMAT(X,'Temp/K Xex(',I<iDim>,') Xex(',I<iDim2>,') Gex(',I<iDim>,') Gc(',I<iDim>,') Gex(',I<iDim2>,') Gc(',I<iDim2>,'),Gex(',I2,') Gc(',I2,')')
                NCXP = NCXP + 1

	        END DO
	        WRITE(1,1359) STDVA(KS)*100.,STDVB(KS)*100.,STDVC(KS)*100.
	        WRITE(6,*) STDVA(KS)*100.,STDVB(KS)*100.,STDVC(KS)*100.
1359	        FORMAT(/,X,' Desv.Std. %',11X,F10.4,6X,f10.4,7X,f10.4,/)

	      end if

	    END IF

	  END DO
!
	END IF
!
	CLOSE(UNIT=1)
	return
ENDsubroutine !FIN WRITESAL2
!
!------------------------------------------------------------------------------------------------------------------
!
!
!
!	Carga en los lugares adecuados de XX los valores iniciales de los dc que se usaran como cualquier otro parametro ajustable
!
SUBROUTINE FORDC(NUMDC,IDDC,DC,NPAR,XX)
!
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION IDDC(8)
	DIMENSION DC(20)
	DIMENSION XX(8)
!
	IF(NUMDC.NE.0)then
!
	  DO I= (NPAR-NUMDC+1),NPAR

	    XX(I)=DC(IDDC(I-npar+numdc))

	  END DO

	endif
	return
ENDsubroutine !FIN FORDC
!--------------------------------------------------------------------------------------------------------------
!
!	Falsa subrutina de lectura, para traspasar al paquete termodinámico
!	que usa GC3 las variables leídas desde ParaEst.
!
!	Una variable de entrada no puede estar en un COMMON, por lo tanto,
!	hay que asignar la entrada de esta subrutina a las variables que
!	van en dicho COMMON.
!
subroutine fakePARGR(NC,NG,NGA,NST,IOUt,NUWW,QAREA,TS,GS,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,ONUASOCWW, &
&	               MASSOC,ENASS,RKASS,Tspline)

	IMPLICIT REAL*8(A-H,O-Z)
	parameter(NCM=15,NGM=15,NGAM=12,NSM=24)

	DIMENSION qArea(NG),Ts(NG),nuWW(NG,NC),Gs(NG),GDOT(NG),GDDOT(NG),rKij(NG,NG),rKTij(NG,NG), &
&	          alfaij(NG,NG),Tspline(NG)
!
!	Energía y volumen de asociación con 4 subíndices: sólo están en
!	esta subrutina.
	dimension enAss(2,NGA,2,NGA),rkAss(2,NGA,2,NGA),onuAsocWW(NGA,NC),MAssoc(NGA)
	integer sigma
	real*8 kappa
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),XKIJ(NGM,NGM),      &
&	              AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)

	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM)
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
	common/NAME/AName(NCM)
	COMMON/COORD/ZZ

	if(NC.GT.NCM)then

	  write(*,*)'Error: número de compuestos máximo sobrepasado en GC3'
	  stop

	elseif(NG.GT.NGM)then

	  write(*,*)'Error: número máximo de grupos atractivos sobrepasado en GC3'
	  stop

	elseif(NGA.GT.NGAM)then

	  write(*,*)'Error: número máximo de grupos asociativos sobrepasado en GC3'
	  stop

	endif
	iOut = 20
!
!	Número de coordinación:
	ZZ=10.D0
!
!	Constante de los gases ideales
	R=82.05D0	!atm cm3/(mol K)
	NGROUP=35
	NownII = NG
	NownIJ = ((NG + 1)**2 - NG -1)/2 - NG
!
!	Variables equivalentes dispersivas
12	WRITE(iOut,105) (J,J=1,NC)
	WRITE(iOut,120)
	do j = 1,NG
!
	  q(j) = qArea(j)
	  Tstr(j) = Ts(j)
	  ny(:NC,j) = nuWW(j,:NC)
	  gStr(j) = Gs(j)
	  g1(j) = GDOT(j)
	  g2(j) = GDDOT(j)
	  Tspl(j) = Tspline(j)
	  epx(j) = -4.D0*(g1(j)*Tspl(j)/Tstr(j) + g2(j))
15	  WRITE(IOUT,106)j,TSTR(j),Q(j),GSTR(j),G1(j),G2(j),(NY(i,j),i=1,NC)
	  do k = 1,NG

	    xKij(j,k) = rKij(j,k)
	    aKij(j,k) = rKTij(j,k)
	    alfa(j,k) = alfaij(j,k)

	  enddo

	enddo

!-----------------------------------------------------------------------
!
!	Impresión de parámetros binarios
33	WRITE(iOut,115)
	do i = 1,NG

	  WRITE(IOUT,116) (XKIJ(I,J),J=1,NG)

	enddo
      WRITE(iOut,117)
	do i = 1,NG

	  WRITE(IOUT,116) (AKIJ(I,J),J=1,NG)

	enddo
	WRITE(IOUT,118)
	DO 23 I=1,NG

	  WRITE(IOUT,116) (ALFA(I,J),J=1,NG)

23	enddo
!-----------------------------------------------------------------------
!	Variables equivalentes asociativas
!-------Conversion de variables a notacion por sitio "m" (y no, por sitio
!	  "k" en el grupo "i")
	if(NGA.GT.0)then

	  m1 = 0
	  m2 = 0
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      m2 = m1
	      m1 = m1 + 1      	!Defino actual sitio "k" como sitio "m1"
	      do j = i,NGA 	!esta definicion solo es valida gracias a que hasta ahora, solo ha habido
	        if(i.EQ.j)then	!un tipo de sitio en cada grupo.
	          do l = k,Massoc(j)

	            m2 = m2 + 1 	!defino actual sitio "l" como sitio "m2"
	            eps_R(m1,m2) = enAss(l,j,k,i)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = rkAss(l,j,k,i)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        else
	          do l = 1,Massoc(j)

	            m2 = m2 + 1 	!defino actual sitio "l" como sitio "m2"
	            eps_R(m1,m2) = enAss(l,j,k,i)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = rkAss(l,j,k,i)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        endif
	      enddo

	    enddo
	  enddo
	  NST = sum(Massoc(:NGA))
!
!	  Reasignación en notación 1..NST:
	  do i = 1,NC

	    m1 = 0
	    do j = 1,NGA
	      do k = 1, MAssoc(j)

	        m1 = m1 + 1
	        sigma(m1,i) = onuAsocWW(j,i)

	      enddo
	    enddo

	  enddo
!
!-------Impresion parametros asociativos
	  WRITE(IOUT,217) NGA
217	  FORMAT(//,X,"Número de grupos considerados en la asociación =",I3)
!
	  WRITE(IOUT,218)
218	  FORMAT(/,X,' Configuración grupal asociativa ',/,30X,'Compuesto Nro.')

	  write(iOut,'(" Grupo    Nro. de sitios  ",<NC>(I5,X))')((I), I=1,NC)
!
	  write(iOut,120)
	  DO I=1,NGA

	    WRITE(IOUT,219) I,MAssoc(I), (onuAsocWW(I,J), J=1,NC)
219	    FORMAT(2X,I3,16X,I3,2X,<NC>(F5.0,X))

	  END DO
!
	  WRITE(iout,220)
 220	  FORMAT(/,X,'Energía de asociación [épsilon/R] (K)',/)

!	  Escritura de parámetros energéticos de asociación en forma de matriz
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write(IOUT,2210)((enAss(k,i,l,j),l = 1,MAssoc(j)),j=1,NGA)

	    enddo
	  enddo
2210	  format(3X,'|',<NST>(X,F7.1),'|')
!	  Escritura del volumen de asociación
	  write(IOUT,'(/," Volumen de asociación [kappa] (cm3/mol)",/)')
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write(IOUT,2220)((rKAss(k,i,l,j),l = 1,Massoc(j)),j=1,NGA)

	    enddo
	  enddo
2220	  format(3X,'|',<NST>(X,F7.5),'|')

	endif


100	FORMAT(20I3)
101	FORMAT(5F10.2)
102	FORMAT(I3,5F10.2)
111	FORMAT(8F7.2)
116	FORMAT(3X,'|',<NG>F8.3,'|')
115	FORMAT(/,1X,'Parámetros de interacción binaria, [k*]',/)
117	FORMAT(/,1X,'Dependencia con la temperatura, [k`]',/)
118	FORMAT(/,1X,'Parámetros de composición local, [alfa]',/)
119	FORMAT(2I3,4F10.2)
105	FORMAT(///,2X,'Grupos atractivos:',33X,'Compuesto Nro.',/,3X,'GN',4X,'T*',6X,'q',5X,'g*',10X,'g`',8X,'g"',3X,15I4,/)
106	FORMAT(1X,I4,F7.1,F8.4,F10.1,2F10.5,15I4)
120	FORMAT(1X,90("-"))
130	FORMAT(/,1X,'Los parámetros de interacción binaria entre los ','grupos',I2,' y ',I2,' no están disponibles.',/,1X,'Se usarán', &
&	        'valores estándar a menos que estos sean ingresador por el usuario')
	return
endsubroutine
!---------------------------------------------------------------------------------------------------------
!
!
!
!	Análoga de la subrutina anterior pero para PARMOL
!
!
subroutine fakePARMOL(NC,NG,NST,IOUt,WWDC,WWTEMPCRIT,wwPRCRIT)
!
	implicit real*8(A-H,O-Z)
	parameter(NCM=15,NGM=15,NGAM=12,NSM=24)
	character*10 AName
!
	DIMENSION TSAT(NCM),PS(NCM),Z(2),OMEGA(NCM),NYOLD(NGM),EXPOL(NGM)
	DIMENSION wwDc(NC),wwTempCrit(NC),wwPRCRIT(NC)
!
	common/versSUB/versGCA
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM)
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
!
	common/NAME/AName(NCM)
	common/GCPROM/pmm(ncm),Pen(NCM),HHA(ncm),HHB(ncm),hhc(ncm),HHD(ncm),HHE(NCM),HHF(NCM)
	common/ENTALP/iEntalp
	COMMON/COORD/ZZ

	dimension eps_R_bak(NST,NST)
	real*8 kappa_bak(NST,NST),kappa
	integer sigma,sigma_bak(NST,NC)

	write(iOut,'(//,"Diámetros críticos",/)')
	iEntalp = 0
	do i = 1,NC

	  Tc(i) = wwTempCrit(i)
	  dc(i) = wwDC(i)
	  Pc(i) = wwPRCRIT(i)
	  Tsat(i) = 0.
	  Ps(i) = 0.
	  omega(i) = 0.
	  write(AName(i),'("Comp ",I2)')i
	  PMM(i) = 0.
	  Pen(i) = 0.
	  HHA(i) = 0.
	  HHB(i) = 0.
	  HHC(i) = 0.
	  HHD(i) = 0.
	  HHE(i) = 0.
	  HHF(i) = 0.


	enddo
	WRITE(IOUT,10)
	WRITE(iOut,13)
	do i = 1,NC

	  WRITE(IOUT,11) AName(i),TC(I),PC(I),OMEGA(I),TSAT(I),PS(I),DC(I),PMM(i),Pen(i)

	enddo
!
!	Constantes de Passut-Danner: impresión sólo si al menos se le in-
!	gresa valores a un compuesto:
	iEntalp = NC
	do i = 1,NC
	  if((HHA(i).EQ.0).AND.(HHB(i).EQ.0).AND.(HHC(i).EQ.0).AND.(HHD(i).EQ.0).AND.(HHE(i).EQ.0).AND.(HHF(i).EQ.0))then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	if(iEntalp.NE.0)then

	  write(iOut,'(/,"Constantes de Passut-Danner para entalpía del ","gas ideal",/,"   Compuesto  A(BTU/lb)  B(BTU/(lb ºR))  ", &
     &             "C(BTU/(lb ºR2)·1E3  D(BTU/(lb ºR3))·1E6  E(BTU/(lb ºR4))·1E10  F(BTU/(lb ºR5))·1E14",/,128("-"))')
	  do i = 1,NC

	    write(iOut,14)AName(i),HHA(i),HHB(i),1.D3*HHC(i),1.D6*HHD(i),1D10*HHE(i),1D14*HHF(i)

	  enddo

	endif
91	FORMAT(8F10.4)
10	FORMAT(//,2X,'Propiedades de compuesto puro',//,3X,'Compuesto',5X,'Tc(K)',1X,'Pc(atm)',3X,'omega',2X,'Tsat(K)',2X,'Psat(atm)',2X, &
&	       'dc(cm/mol^1/3)',3X,"M(g/mol)",3X,"c(cm3/mol)")
11	FORMAT(4X,A10,1X,2F7.1,3X,F6.4,F9.2,2X,F8.4,5X,F8.4,7X,F6.2,8X,F6.2)
12	FORMAT(//,2X,'Propiedades de compuesto puro',/,3X,'Nro.',2X,'Tc(K)',2X,'Pc(MPa)',3X,'omega',2X,'Tsat(K)',2X,'Psat(MPa)',2X,  &
&	       'dc(cm/mol^1/3)',X,"M(g/mol)",X,"c(cm3/mol)")
13	format("  ",100("-"))
14	format(3X,A10,X,F9.6,3X,F9.6,10X,F9.6,12X,F9.6,13X,F9.6,13X,F9.6)
	return
endsubroutine