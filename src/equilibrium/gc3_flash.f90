!
!	paraest copy of GC3-1.5.3.1.f90 (GC3/src), made by tools/port_gc3.py.
!	Program GC3 turned into
!	  subroutine GC3Flash (NCOMP, TEMP, PRES, feed, xFrac, Zcomp, indFase, iErrTP)
!	  feed(NCOMP) at TEMP (K), PRES (atm) -> phases xFrac(:,j), Zcomp(j),
!	  indFase(j) (THERMO phase type), j = 1..NF; iErrTP = 1 when the result
!	  is a single phase, after more than maxIt phase splits or more than 4
!	  phases. Changes:
!	  - no GCTHREE.DAT, datter.dat or Compres.dat; GCTHREE.SAL on unit 33.
!	  - model from module thermo_model; no solid phases (solidPhase = 0).
!	  - INDATA does not call ReadParameters (paraest has loaded the model).
!	The other routines are those of the original.
!
!
!       Multiphase P/T Flash program, version 3, may 1983
!       Copyright: Michael L. Michelsen, Instituttet for Kemiteknik,
!                  Bygning 229, DTH, Lundtoftevej 100, DK 2800 Lyngby.
!                  DENMARK;   TEL 02-88 32 88
!
!                  THERMODYNAMIC SUBROUTINES BY STEEN SKJOLD-JOERGENSEN
!
!     M a i n     p r o g r a m  for multiphase flash
!     ************************************************
!
!     Common areas in program complex do not interfere with common
!     areas in thermodynamic part except in /STYR2/ In the present
!     interface routine SRK
!
!		-Esto ya no es así, a partir de la versión 1.2.0 de la GC3: se impri-
!	       men tanto la versión del paquete termodinámico (leído dentro del pa-
!	       quete) como la versión de esta subrutina. Implica un nuevo COMMON.
!			Francisco, 20/05/2011, durante v1.2.0
!
!		-Además, ahora hay un common en "iEntalp" para ver si calcula o no entalpías
!	-------------------------------------------------------------
!
!	Aprovechando las modificaciones introducidas en la versión 1.3.3 del paque-
!	te termodinámico GCA, se modifica esta subrutina para poder imprimir valo-
!	res de Z (coef de compr) afectados por una corrección de Peneloux.
!
!	El usuario debe tener en cuenta que esta corrección no es "mágica", se pasa
!	de tener una predicción pésima, a una simplemente mala.
!
!	La corrección de Peneloux ingresada es siempre negativa (al menos en la ver-
!	sión 1.3.3), ya que el V predicho por la GC es demasiado grande. La línea de
!	ingreso es la siguiente a las propiedades críticas.
!
!	Ej: líneas correspondientes al compuesto "H2O":
!
!	H2O                                                               !nombre
!	647.13    217.67    0.345     373.15     1.0000     2.5763        !Tc  Pc  w  Ts  Ps   dc
!	18	-10	0	0	0	0	0	0				!PM  Peneloux   coeficientes de H(gas ideal) de Passut-Danner
!
!
!	v1.2.90 (1.3.0 beta0): comienzo de soporte para GCA-v2.0
!	     91                cambiada dimensión de versSUB
!
!	v1.2.92 	           Limpieza de algunas impresiones
!	                       Cálculo de entalpías
!
!	v1.2.93: limpieza de código.
!	         algunas identaciones
!	         agregado parámetros NCM
!
!	v1.2.96: PARAGC llamado con nTemp = 1, pero no así GCEOS.
!
!	v1.2.100: las versiones 98-99 están falladas. Esta comienza siendo un clon de la 96, incluyendo modificaciiones
!	          hechas en 97 y 98.
!
!	v1.2.101. Rev. 3. En la subrutina ONEDIR, luego del llamado al paquete termodinámico,
!	                  Si el elemento d(ln phi(i))/dn(i) = 0, le asigno un valor pequeño,
!	                  pero no 0. Esto se debe a que en ciertos casos, la variable SSQ puede
!	                  tomar un valor nulo, generando dps un 1/0.
!
!	v1.3.0. Inclusión del paquete thermopack. 
!
	subroutine GC3Flash (NCOMP, TEMP, PRES, feed, xFrac, Zcomp, indFase, iErrTP)
!
	use thermo_model, only: model
	IMPLICIT REAL*8(A-H,O-Z)
	
	integer, parameter      :: inputFile = 5, NCM = 30, NFM = 4, outputFile = 33, maxIt = 50
	real(8), parameter      :: R2 = 82.05D0
	character, parameter    :: TAB = achar(9)
	dimension                  feed(NCM), xFrac(NCM,NFM), Zcomp(NFM), indFase(NFM)
	
	integer                    solidPhase
	DIMENSION                  Y(NCM), yp(NCM), yt(NCM), DH_n(NCM), YVAL(NCM*NFM)
                                
	COMMON /DERIV/             FX(NCM,NCM)
	COMMON /KEEP/              T,P,GNUL,Z(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /SOLID/             FUGAR(NCM), ZLOG(NCM), Zcomp_pure(NCM), ISOLD(NCM), ISTAT, ID_S
	COMMON /STORE/             GRAD(NCM*(NFM-1)),XMAT(NCM*(NFM-1),NCM*(NFM-1)),WORK(NCM*(NFM-1),NFM+1)
	COMMON /DELST/             YV1
	COMMON /AKTIV/             IVAL(NCM),IVAR(NCM*(NFM-1)),NACT
	COMMON /LIMIT/             EPSV
	COMMON /NTEST/             GNORM,GLIM, NFA,NDUM
	COMMON /MOLFRC/            ZFRAC(NCM,NFM)
! 	COMMON /OLDZ/              ZOLD1,ZOLD2
	COMMON /OLDT/              TOLD
	common /OUTUNIT/           iOut
	COMMON /SOLPRP/            TMELT(NCM), HMELT(NCM), ACpMelt(NCM), BCpMelt(NCM), vMelt(NCM), Ptr(NCM)
!
!	SUbroutine arguments for Murray routine
	EXTERNAL                   FSTAB,PHASN
!
!	Nuevas declaraciones incluidas en la versión 1.2.0------------------------
	real(8), dimension(NCM) :: HHA, HHB, HHC, HHD, HHE, HHF, HHG, Hgi, Pen, PMM
	character(10)           :: aname(NCM)
! 	character(20)           :: versGCA
	character(20)           :: version='1.5.3.0'
! 	common /versSUB/           versGCA
	common /NAME/              AName
	common /GCPROM/            PMM, Pen, HHA, HHB, HHC, HHD, HHE, HHF, HHG
	common /ENTALP/            iEntalp
!
!	Versión asociativa del common "ZZZ"
	real(8)                    Kpart(NCM)
!	--------------------------------------------Francisco, 20/05/2011---------
!
	npunto = 0    !selva
	TOLD=0.D0
! 	GLim = 1D-8
	GLIM=1.D-16
10	FORMAT(8F10.5)
500	format(/,72("-"))
!------INPUT NUMBER OF COMPONENTS AND TOTAL NUMBER OF GROUPS IN MIXTURE
!
	IOUT = outputFile
	OPEN (UNIT = IOUT, FILE='GCTHREE.SAL', STATUS = 'unknown',FORM = 'FORMATTED')
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
	solidPhase = 0	! paraest: no solid phases (yet)
	iErrTP = 0
	nIt = 0

20	FORMAT(16I3)
!
!     Indat sets up tables and common areas for termodynamic subroutine
	CALL INDATA (model, NCOMP, IOUT, solidPhase)
	N=NCOMP
!
!	Read mixture composition, normalization not needed
	Z(:N) = feed(:N)
	T = TEMP
	P = PRES
!260	READ(5,*)(Z(I),I=1,N),T,P
	suma = 0
	indF(:NF) = 0	!restart
	do  i = 1, N

!	  Zero mole fraction not permitted (log-terms)
	  z(I) = z(I) + 1.D-20
	  suma = suma + z(i)
	  
	enddo
	if (suma.LT.1.D-10) GOTO 250
!	NORMALIZE:
	z(:N) = z(:N)/suma
!
!	Input temperature and pressure
!	Units must correspond to (absolute) units in indat
	iSTAT=0
	ND=0
	IF (T.LE.0.) GOTO 250
!
	CALL SRK (3, 0, IC, T, P, ZX, Z, Y, Yp, FX)
!	
!	"y" es ln(phi), por más que aparezca 2 veces.
!	"iC" es el tipo de fase que se encontró para una misma TP, con menor G (-1 para vapor, 1 para líquido).
!	"Z" la alimentación
!	"ZX" el factor de compresibilidad
!	"FX" es la matriz d[ln(phi)]/dn
!
	iTypOrig = iC	!almaceno el tipo de fase original.
! 	Kpart(:N) = exp(y(:N))
!
!	EPSV is tolerance for trace components, it is modified according
!	to the result of the stability test
	EPSV=1.D-4
	WRITE (6, 60) T, P, ZX
	WRITE (IOUT, 60) T, P, ZX
! 	WRITE(9,60) T,P,ZX, iTypOrig, Xs(1,:3)				!Fran 07/08/09 -> JA!!
! 60	FORMAT(1H1,///,' PT-Flash, T = ',F8.2,' K, P = ',F8.2,' atm, Z = ',F10.8, X, i2, 3(2X,F10.8))
60	FORMAT(1H1,///,' PT-Flash, T = ',F8.2,' K, P = ',F8.2,' atm, Z = ',F10.8)
	GNUL=0.D0
	SFAS(1)=1.D0	!fracción de fase 1
	DO 70 I=1,N

	  xVL(I,1)=1.D0
!	  Zero mole fraction not permitted (log-terms)
	  z(I) = z(I) + 1.D-20
	  zFrac(I,1) = z(I)
	  zLog(I) = dlog(z(I))
	  fr(I) = zLog(I) + Y(I)	! fr(i) = d(i), siguiendo definición en pág. 232 del libro.
	  Gnul = Gnul + z(I)*y(I)

70	enddo
!
	WRITE (6, 80) (I,Z(I),FR(I),I=1,N)
	WRITE (IOUT, 80) (I,Z(I),FR(I),I=1,N)
80	FORMAT(/, 8X, 'Composition  ln z(i)*phi(i) ', /, 36("-"), /, (I3, 2X, F12.6, F14.6))
	NF=1
!
!	Initiate stability analysis
90	CALL PURE (Y, N, NONHY, FX)  !<---------------------------------------------<
!	                                                                        |
	IF (YV1.GE.-.01) GOTO 110
!	                                                                        |
!	  Sufficiently unstable, go to phase split
	  write (*, 190)
	  WRITE(IOUT,190)
	  yVal(:N) = y(:N)/z(:N)
100	  continue
!	                                                                        |
	  GOTO 210
!	                                                                        |  
!     Otherwise locate minimum by Murray method
! 110	DO 120 I=1,N
! !	                                                                        |
! 120	  YVAL(I)=DLOG(Y(I))
110	yVal(:N) = dlog(y(:N))
!	                                                                        |
      xLAM=1.D0
      write (*, 130)
      WRITE(IOUT,130)
130	FORMAT(' ')
	DO 140 I=1,N
!	                                                                        |
	  IVAR(I)=I
	  IVAL(I)=1
!	                                                                        |
140	enddo
      NACT=N
!	                                                                        |
!	Murray search for stability minimum
	CALL MURML (NCM*(NFM-1), N, 1, 10, xLAM, 1.D-12, FUN, yVal, GRAD, XMAT, WORK, FSTAB)
	IF (DABS(FUN) .LT. .01) EPSV=1.D-3*DSQRT(DABS(FUN))
!	                                                                        |	
	write (*, 150) NACT
	WRITE(IOUT,150) NACT
150	FORMAT(/,'  Variables active ',I4/)
!	                                                                        |
!	Check result of search
! 160	IF (FUN.LT. -1D-7) GO TO 180		!Valor original
! 160	IF (FUN.LT.-1D-10) GO TO 180		!Fran 19/07/09
160	IF (FUN.LT.-1D-5) GO TO 180		!Fran 19/07/09
!
!	Try onedimensional search if pure is negative
!	Tolerance set in accordance with convergence criterion.
!	Use of a smaller value requires a smaller Glim
	iTry=0
	DO 410 I=1,N
!	                                                                        |
!	  IF (ISOLD(I) == 0) GOTO 410
	  if (iSold(i) /= 0) then
!	                                                                        |
	    dif = fugaR(i) - fr(i)
	    if (dif.LE.-1D-7) then
!	                                                                        |
	      iTry = 1
	      iSold(i) = 2
!	                                                                        |
	    endif
!	                                                                        |
	  endif
!	                                                                        |
410	enddo
!	                                                                        |
!	DIF=FUGAR(I)-FR(I)
!	IF (DIF > -1.D-7) GOTO 410
!	                                                                        |
!	  ITRY=1
!	  ISOLD(I)=2
!	                                                                        |
!410	CONTINUE
!	                                                                        |
	IF (ITRY == 0) GOTO 420
!	                                                                        |
	  ISTAT=ISTAT+ITRY
	  write (*, 403)
	  WRITE(IOUT,403)
403	  FORMAT(/,'   Solid phase formation  ')
	  GOTO 218
!	                                                                        |
420	CONTINUE
!	                                                                        |
!	One dimensional search for close minima as last resort
	CALL ONEDIR (IOUT, N, YV1, YVAL, FX)
!	                                                                        |
	IF (YV1 .LT. -1.D-7) GOTO 180
!	                                                                        |	
	write (*, 170)
	WRITE(IOUT,170)
170	FORMAT(/,'  Solution is now stable * * * ')
	nPunto = nPunto + 1
!	paraest: the phases are returned; a single phase is an error.
	if (NF > 1) then
	  xFrac = zFrac
	  Zcomp = zVal
	  indFase = indF
	else
	  iErrTP = 1
	endif
	return
!	                                                                        |
!	No more phases could be found
!	                                                                        |	                                                                        |
!	GOTO 50
	goto 260
! 	
! 
! -----------------------------------------------------------------------------------------------------------
! 
180	WRITE (IOUT, 190)
190	FORMAT (/,' System is still unstable, phase split performed '/)
	nIt = nIt + 1
	if (nIt > maxIt) GOTO 250
!	                                                                        |
	yVal(:N)=DEXP(yVal(:N))/z(:N)  !yVal(i) = phi(i)
200	continue
210	NF = NF + 1
!	                                                                        |
!	Add new phase and proceed with phase split
	IF (NF > 4) GOTO 250
!	                                                                        |
!	Four phases maximum for present version
	DO 215 J = 2, NF
!	                                                                        |
	  JN = NF + 2 - J
	  SFAS(JN) = SFAS(JN-1)
!	                                                                        |
	  DO 215 I = 1, N
!	                                                                        |
215 	    XVL(I,JN) = XVL(I,JN-1)
!	                                                                        |
	DO 216 I=1,N
!	                                                                        |
	  XVL(I,1)=1.D-16*yVal(I)
	  DO 216 J=1,NF
!	                                                                        |
216	    yVal(I+(J-1)*N)=xVL(I,J)
!	                                                                        |
	SFAS(1)=1.D-16
	GOTO 419
218	IF (ISTAT == 1) GOTO 414
!	                                                                        |
	  ISTAT=1
	  GOTO 419
!	                                                                        |
!	Usa a NF luego, como fase de referencia
414	NF=NF+1
	SFAS(NF)=0.
	ZVAL(NF)=0.
	DO 416 I=1,N
!	                                                                        |
416	  XVL(I,NF)=0.
!	                                                                        |
419	CONTINUE
	XLAM=.2
	IF (NF == 2) XLAM=.5
	M=(NF-1)*N
!	                                                                        |
	write(*,*)NF,N,M
!	                                                                        |
	DO 220 I=1,M
!	                                                                        |
	  IF (I.LE.N) IVAL(I)=1
220	  IVAR(I)=I
!	                                                                        |
	NACT=M
!	                                                                        |
!	Phase split by direct substitution/Gdem
	CALL DSUB (N, FX, YVAL)
	do i = 1, NF
!	  avoid zero-amount phases... 
	  if (sFas(i) <= 0D0) sFas(i) = 1D-20
      enddo
!	                                                                        |	
!	Check for reduction of number of phases
	NF=NFA
	M=(NF-1)*N
	IF (Gnorm.LT. Glim) GO TO 230  !----------V                             |
!	                                          |                             |
!	  Phase split by murray method
        CALL MURML (NCM*(NFM-1), M, 1, 10, XLAM, GLIM, FUN, YVAL, GRAD, XMAT, WORK, PHASN)
!	                                          |                             |
230	EPSV=1.D-4  !<----------------------------<                             |
	write (*, 150) NACT
	WRITE(IOUT,150) NACT
!	                                                                        |	
240	CALL OUT (YVAL, N, FX)
!	                                                                        |
	GO TO 90  !>----------Back to stability test----------------------------^
!
!
250	iErrTP = 1
260	return
	END subroutine GC3Flash
!
!	END OF PROGRAM


!***********************************************************************
!
!     INDATA SETS UP DATA TABLES FOR ACTUAL MIXTURE
!
!***********************************************************************
!
	SUBROUTINE INDATA (model, NC, IOU, solidPhase)
	IMPLICIT REAL*8 (A-H,O-Z)
	
	integer, parameter      :: inputFile = 5, NCM = 30
	real(8), parameter      :: R1 = 8.314, R2 = 82.05
	
	integer                 :: i, NCS, solid_ID(NCM), solidPhase
	real(8), dimension(NCM) :: ACpMelt, aux, BCpMelt, Hmelt, Ptr, Tmelt, vMelt
	
	COMMON /SOLPRP/            TMELT, HMELT, ACpMelt, BCpMelt, vMelt, Ptr
	COMMON /CRIT/              TCR(NCM), PCR(NCM)
	COMMON /STYR2/             model_common, N, NDER, NTEMP
	COMMON /SKIP/              NOTRY(NCM)
	                      
	CHARACTER(10)           :: compoundName(NCM)
	COMMON /NAME/              compoundName

	N = NC
	DO I = 1, NC
    
	  TMELT(I) = 1
	  HMELT(I) = 0
	  vMelt(i) = 0
	  ACpMelt(i) = 0
	  BCpMelt(i) = 0

	enddo	
	noTry(:NC) = 1
!  	noTry(:NC) = 0	!activate this line to avoid pure compounds in stability test (save time fo large NC)
	iOut = iOu
!	paraest: the parameters are already loaded (ReadParameters not called)
	if (solidPhase > 0) then	
	  
	  write (iOu, '(/, " Solid phase enabled", /)')
	  read (inputFile, *) NSC
	  read (inputFile, *) solid_ID(:NSC)
	  read (inputFile, *) aux(:NSC)
	  do i = 1, NSC
	    TMelt(solid_ID(i)) = aux(i)
	  enddo
	  read (inputFile, *) aux(:NSC)
	  do i = 1, NSC
	    Ptr(solid_ID(i)) = aux(i)
	  enddo	  
	  read (inputFile, *) aux(:NSC)
	  do i = 1, NSC
	    HMelt(solid_ID(i)) = aux(i)
	  enddo
	  read (inputFile, *) aux(:NSC)
	  do i = 1, NSC
	    vMelt(solid_ID(i)) = aux(i)
	  enddo
	  read (inputFile, *) aux(:NSC)
	  do i = 1, NSC
	    ACpMelt(solid_ID(i)) = aux(i)
	  enddo
	  read (inputFile, *) aux(:NSC)
	  do i = 1, NSC
	    BCpMelt(solid_ID(i)) = aux(i)
	  enddo
	  write (iOu, '(" Compound", 6X, "Tfus(K)", 5X, "Pfus(atm)", 11X, "DHfus(J/mol)", 5X, "Dvfus(cm3/mol)", 5X,  &
&	                "ACpfus(J/molK)", 5X, "BCpfus(J/molK2)", /, " ",116("-"))')
	  do i = 1, NSC
	    	    
	    write (iOu, '(X, A10, 4X, F7.2, 4X, G12.4, 9X, G12.4, 9X, F7.2, 11X, G12.4, 8X, G12.4)') compoundName(solid_ID(i)), TMelt(solid_ID(i)), Ptr(solid_ID(i)), HMELT(solid_ID(i)), vMelt(solid_ID(i)), ACpMelt(solid_ID(i)), BCpMelt(solid_ID(i))
	    HMELT(solid_ID(i)) = HMELT(solid_ID(i))/R1/Tmelt(solid_ID(i))
	    vMelt(solid_ID(i)) = vMelt(solid_ID(i))/R2
	    ACpMelt(solid_ID(i)) = ACpMelt(solid_ID(i))/R1
	    BCpMelt(solid_ID(i)) = BCpMelt(solid_ID(i))/R1	    
	    
	  enddo
	  write (iOu, *)
	
	endif
	model_common = model
	RETURN
	ENDsubroutine INDATA
!	------------------------------------------------------------------
!
!
!     THIS IS THE END OF THERMODYNAMIC SUBROUTINES
!
!     START GENERAL EQUATION SOLVING ROUTINES AND MINIMIZER
!
!     LU PERFORMS TRIANGULAR DECOMPOSITION OF A FOR SOL'N OF LINEAR EQ.
!
      SUBROUTINE LU(ND,N,IPIV,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IPIV(ND),A(ND,ND)
      IPIV(N)=N
      N1=N-1
      IF (N1 == 0) GO TO 50
      DO 40 I=1,N1
      X=A(I,I)
      IF (X.LT.0.) X=-X
      IPIV(I)=I
      I1=I+1
      DO 10 J=I1,N
      Y=A(J,I)
      IF (Y.LT.0.) Y=-Y
      IF (Y.LE.X) GO TO 10
      X=Y
      IPIV(I)=J
   10 CONTINUE
      IF (IPIV(I) == I) GO TO 30
      K=IPIV(I)
      DO 20 J=I,N
      X=A(I,J)
      A(I,J)=A(K,J)
   20 A(K,J)=X
   30 DO 40 J=I1,N
   
        if (A(i,i) == 0) A(i,i) = 1.D-20
	    
        X=-A(J,I)/A(I,I)
        A(J,I)=X
        DO 40 K=I1,N
   40   A(J,K)=A(J,K)+X*A(I,K)
   
   50 RETURN
      ENDsubroutine LU
!	------------------------------------------------------------------
!
!
!     BACK SOLVES LINEAR EQUATIONS A*X=V BY BACK SUBSTITUTION
!     RESULT IN V
!
      SUBROUTINE BACK(ND,N,IPIV,A,V)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IPIV(ND),A(ND,ND),V(ND)
      N1=N-1
      IF (N1 == 0) GO TO 30
      DO 20 I=1,N1
      I1=I+1
      K=IPIV(I)
      IF (K == I) GO TO 10
      X=V(I)
      V(I)=V(K)
      V(K)=X
   10 DO 20 J=I1,N
   20 V(J)=V(J)+A(J,I)*V(I)
   30 V(N)=V(N)/A(N,N)
      IF (N1 == 0) GO TO 60
      DO 50 II=2,N
      I=N+1-II
      I1=I+1
      DO 40 J=I1,N
   40 V(I)=V(I)-A(I,J)*V(J)
   50 V(I)=V(I)/A(I,I)
   60 RETURN
      ENDsubroutine BACK
!-----------------------------------------------------------------------
!
!
!	MURML performs the minimization using the murray modified Cholesky
!	decomposition of the hessian matrix
!
	SUBROUTINE MURML(ND,N,IPR,NMAX,XLAM,GLIM,F,X,GD,G,W,FUNC)
	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4

	DIMENSION             X(NCM*NFM),GD(ND),G(ND,ND),W(ND,5)
	COMMON /AKTIV/        IVAL(NCM),IVAR(NCM*(NFM-1)),NACT
	common /OUTUNIT/      iOut	
	
	NTM=0
	DEL=1.D-7
!
!	The end of the file jumps here in this way:
!	if (NTM.LT.NMax .AND. GPS > Glim) then
!	
10	NTM=NTM+1	!<-------------------------------------------<
!	                                                         |
	GPS=0.
!	                                                         |
!     FUNC is the subroutine call to the function to be minimized
!	                                                         |
	CALL FUNC(N,2,F,GD,G,X)
	DO 20 I=1,N
!	                                                         |	
	  W(I,4)=X(I)
	  W(I,5)=GD(I)
	  IF (I > NACT) GOTO 20
!	                                                         |	  
	    IV=IVAR(I)
	    GPS=GPS+GD(IV)**2
!	                                                         |   
20	CONTINUE
!	                                                         |
!	GLIM is the tolerance for the squared gradient norm
!     If result is less than tolerance, exit
!
	IF (GPS.LT.GLIM) GOTO 170
	BETA=0.
!	                                                         |
!     Set up hessian matrix for decomposition
!	                                                         |
	DO 30 I=1,NACT
!	                                                         |      
	  IV=IVAR(I)
	  DO 30 J=1,I
!	                                                         |      
	    JV=IVAR(J)
	    G(I,J)=G(IV,JV)
	    S=DABS(G(I,J))
	    IF (I == J) S=S*NACT
	    IF (S > BETA) BETA=S
30	    W(I,1)=0.
!	                                                         |
	BETA=DSQRT(BETA/NACT)
!	                                                         |
!	Subroutine split performs decomposition
!	                                                         |
	CALL SPLIT(ND,N,IDI,BETA,DEL,G,W)
	XLM=0.
	NTS=0
!	                                                         |
!	                                                         |
!	There's a GOTO loop, which starts here. It repeats 
!	optimization without modify the hessian if results
!	are "promising".
!	                                                         |
40	NTS = NTS + 1  !<-----------------------------------<    |
!	                                                    |    |
!	LINE solves the linear equations.
!	Using XLM>0 gives a Marquardt-like procedure (not used here)
!	                                                    |    | 
	CALL LINE(ND,N,XLM,GD,G,W)
!	                                                    |    | 
	SMAX=0.
	GP=0.
	DP=0.
!	                                                    |    | 
	GPS=0.
	GGG=0.
      DO 60 I=1,N
!	                                                    |    | 
	  IF (I > NACT) GOTO 50
!	                                                    |    | 
	    IV=IVAR(I)
	    GGG=GGG+GD(IV)**2
!	                                                    |    | 
50	  CONTINUE
	  S=W(I,3)
!	                                                    ^    ^
	  IF (DABS(S) > SMAX) SMAX=DABS(S)
!	                                                    |    | 
	  GP=GP+S*GD(I)
!	                                                    |    | 
!	  Gradient norm?
	  GPS=GPS+GD(I)**2
60	  DP=DP+S*S*W(I,1)
!	                                                    |    | 
	FAC=1.
!	                                                    |    |     
	IF (SMAX > XLAM) FAC=XLAM/SMAX
!	                                                    |    | 
!     XLAM is the maximal permitted step, if larger, reduce
!	                                                    |    | 
	DER=FAC*GP
	ALFA=1.
!	                                                    |    | 
!	                                                    |    | 
!	"Update X" loop, as a function of alfa factor, and 
!	change in the objetive function
!	------------------------------------
70	FF=ALFA*FAC	!<--------------------------------^     |    |
!	                                              |     |    |
	DO 80 I=1,N
!	                                              |     |    |
!	  x(i) = w(i,4) + alfa*fac*w(i,3)
	  X(I)=W(I,4)+FF*W(I,3)
!	                                              |     |    |
80	CONTINUE
!	                                              |     |    |
!	Reduction less than E-4: change sign
!	Reduction less than E-5: accept and exit
!	                                              |     |    |
	IF (DABS(ALFA).LT. 1.D-4) ALFA=-ALFA
!	                                              |     |    |
!	                                              ^     ^    ^
!	Exit of the loop, if alfa is lesser than 1E-5
	IF (DABS(ALFA).LT. 1.D-5) GO TO 90
!	                                              |     |    |
	CALL FUNC(N,1,FNEW,GD,G,X)
	DELS=FNEW-F
!	                                              |     |    |
!	Test for better result (worse by E-10 accepted)
!	                                              |     |    |
	IF (DELS.LT.1.D-10) GO TO 90  !>-----V        |     |    |
!	                                     |        |     |    |
!	                                     |        |     |    |
!	                                     |        |     |    |
!	If it has already tryed without changing H, exits this loop,
!	skiping alfa reduction.
	IF (NTS > 1) GO TO 120       !>---------V    |     |    |
!	                                     |   |    |     |    |
!	                                     V   V    ^     ^    ^
!	Energy increased, try with smaller step  
!	                                     |   |    |     |    |
	ALFA=ALFA/3
	GO TO 70	!---------------------------------^     |    |
!	----------
!	                                     |   |          |    |
!	                                     |   |          |    |
!	Calculate predicted improvement (quadratic approximation)
!	                                     |   |          |    |
90	PRED = FF*GP - FF**2/2*(GP + DP)  !<-<   |          |    |
!	                                         |          |    |
	CALPRE = DELS/(PRED + 1.D-19)
	GNW=0.D0
	DO 100 I=1,N
!	                                         |          |    |
	  IF (I > NACT) GOTO 100
!	                                         |          |    |
	    IV=IVAR(I)
	    GNW=GNW+GD(IV)**2
!	                                         |          |    |
100	CONTINUE
	RG=GNW/GPS
	GPS=GNW
	IF (GPS.LT. GLIM) GO TO 170
!	                                         |          |    |
	  F=FNEW
	  DO 110 I=1,N
!	                                         |          |    |
110	    W(I,4)=X(I)
!	                                         |          |    |
	    IF (NTS.GE.5) GOTO 120	!>---------V          ^    ^
!	                                         |          |    |
!	                                         |          |    |
!	Repeat with unmodified hessian up to 5 times, if results
!	are sufficiently promising
!	                                         |          |    |
	    IF (RG.LT.0.02 .AND. CALPRE > 0.3) GOTO 40 !>--^    |
!	-------------------------------
!	                                         |               |
!                                              |               |
120	    DO 130 I=1,N   !<--------------------<               |
!	                                                         |
130	      X(I)=W(I,4)
!	                                                         |
	    IF (IPR > 1) write (*, 160)
!	                                                         |
! 	    iOut = 0
	    IF (IPR > 1) WRITE(IOUT,160)
160	    FORMAT(' ')
!	                                                         |
	IF (NTM.LT. NMAX.AND. GPS > GLIM) GO TO 10	!>-------^
!	----------------------------------------------
!
!
	XLAM=GPS/GLIM
170	RETURN
	ENDsubroutine
!	------------------------------------------------------------------
!
!
!     FSTAB EVALUATES FUNCTION, GRADIENT AND HESSIAN FOR STABILITY
!     ANALYSIS USING TANGENT PLANE DISTANCE
!
	SUBROUTINE FSTAB(N,IV,FUN,GRAD,XMAT,Y)
	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4

	DIMENSION             GRAD(NCM*(NFM-1)),XMAT(NCM*(NFM-1),NCM*(NFM-1)),Y(NCM*NFM),YEX(NCM), &
	                      FG(NCM),FP(NCM)
	COMMON /DERIV/        FX(NCM,NCM)
	COMMON /KEEP/         T,P,GNUL,ZZZ(NCM),FR(NCM),VL(NCM,NFM),SFS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /AKTIV/        IVAL(NCM),IVAR(NCM*(NFM-1)),NACT

      DO 10 I=1,N
!	TRACE COMPONENTS BY DIRECT SUBSTITUTION
      IF (IVAL(I) == 0) Y(I)=FR(I)-FG(I)
   10 CONTINUE
      IF (IV /= 2) GO TO 25
      NACT=0
      IS=0
!
!     ELIMINATE TRACE COMPONENTS
!
      DO 20 I=1,N
      IVAL(I)=1
      IF (Y(I).LT.-15.) IVAL(I)=0
      IF (IVAL(I) == 0) GOTO 20
      IS=IS+1
      IVAR(IS)=I
   20 NACT=NACT+IVAL(I)
   25 SUM=0.
      DO 30 I=1,N
      YEX(I)=0.
      IF (Y(I) > -20.) YEX(I)=DEXP(Y(I))
   30 SUM=SUM+YEX(I)
      JD=1
      IF (IV == 2) JD=3
      CALL SRK (JD,0,IC,T,P,Z,YEX,FG,FP,FX)		!equilibrio
      FUN=1.
      DO 40 I=1,N
      S=Y(I)+FG(I)-FR(I)
      IF (IV == 0) GOTO 40
      GRAD(I)=YEX(I)*S
   40 FUN=FUN+YEX(I)*(S-1.)
      IF (IV.LT.2) GOTO 70
      DO 60 I=1,N
      S=YEX(I)/SUM
      DO 50 J=1,I
      XMAT(I,J)=S*YEX(J)*FX(I,J)
   50 XMAT(J,I)=XMAT(I,J)
   60 XMAT(I,I)=XMAT(I,I)+YEX(I)
   70 RETURN
      END
!	------------------------------------------------------------------
!
!
!     SPLIT FERFORMS THE DECOMPOSITION OF THE HESSIAN FOR MURML
!     IDI IS THE NUMBER OF DIAGONAL ELEMENT CHANGES
!     POSITIVE IDI SHOWS AN INDEFINITE HESSIAN
!
	SUBROUTINE SPLIT(ND,N,IDI,BETA,DEL,G,W)
 	IMPLICIT REAL*8(A-H,O-Z)
	integer, parameter :: NCM = 30, NFM = 4
 	DIMENSION             G(ND,ND),W(ND,NFM+1)
 	COMMON /AKTIV/        IVAL(NCM),IVAR(NCM*(NFM-1)),NACT
!
      IDI=0
      DO 90 I=1,NACT
      I1=I-1
      I2=I+1
      ID=0
      TV=G(I,I)
      SV=TV
      IF (I1 == 0) GOTO 20
      DO 10 J=1,I1
   10 SV=SV-G(I,J)**2
   20 if (SV.LT.DEL*DEL) ID=1
      SVR=DSQRT(DABS(SV))
      if (SVR.LT.DEL)SVR=DEL
      XM=0.
      IF (I == NACT) GOTO 60
      DO 50 J=I2,NACT
      S=G(J,I)
      IF (I1 == 0) GOTO 40
      DO 30 K=1,I1
   30 S=S-G(I,K)*G(J,K)
   40 S=S/SVR
      IF (DABS(S) > XM) XM=DABS(S)
   50 G(J,I)=S
   60 IF (XM.LT.BETA) GOTO 80
      ID=1
      XM=XM/BETA
      DO 70 J=I,NACT
   70 G(J,I)=G(J,I)/XM
      SVR=SVR*XM
   80 IF (ID == 1) W(I,1)=SVR**2-SV
      G(I,I)=SVR
      IDI=IDI+ID
   90 CONTINUE
      DO 100 I=1,N
  100 W(I,2)=0.
      DO 110 I=1,NACT
      IV=IVAR(I)
  110 W(IV,2)=W(I,1)
      DO 120 I=1,N
  120 W(I,1)=W(I,2)
      RETURN
      END
!	------------------------------------------------------------------
!
!
!     LINE CALCULATES THE CORRECTION BY BACKSUBSTITUTION
!
	SUBROUTINE LINE(ND,N,XLAM,GD,G,W)
	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4

	DIMENSION             GD(ND),G(ND,ND),W(ND,NFM+1)
	COMMON /AKTIV/        IVAL(NCM),IVAR(NCM*(NFM-1)),NACT
!
!
      DO 10 I=1,NACT
      IV=IVAR(I)
   10 W(I,5)=GD(IV)
!	 PREPARED FOR MARQUARDT-LIKE VARIANT,BUT HERE XLAM=0
   20 W(1,2)=-W(1,5)/(G(1,1)+XLAM)
      IF (NACT == 1) GOTO 50
      DO 40 I=2,NACT
      I1=I-1
      S=-W(I,5)
      DO 30 J=1,I1
   30 S=S-G(I,J)*W(J,2)
   40 W(I,2)=S/(G(I,I)+XLAM)
   50 W(NACT,3)=W(NACT,2)/(G(NACT,NACT)+XLAM)
      IF (NACT == 1) GOTO 80
      DO 70 II=2,NACT
      I=NACT+1-II
      I1=I+1
      S=W(I,2)
      DO 60 J=I1,NACT
   60 S=S-G(J,I)*W(J,3)
   70 W(I,3)=S/(G(I,I)+XLAM)
   80 DO 90 I=1,N
   90 W(I,2)=0.
      DO 100 I=1,NACT
      IV=IVAR(I)
  100 W(IV,2)=W(I,3)
      DO 110 I=1,N
  110 W(I,3)=W(I,2)
      RETURN
      END
!-----------------------------------------------------------------------
!
!
!	SRK is the interface for the thermodynamic package
!	replacement of srk permits interfacing with different
!	thermodynamic packages, provided the arguments conform.
!	Numerical derivatives can be used for testing of new thermo-
!	dynamics. For production runs, they should be replaced with
!	analytic derivatives to avoid an unnecessary waste of time
!
	SUBROUTINE SRK (INDIC, MTYP, IC, T, P, Z, X, FUG, FUGP, FUGX)
	IMPLICIT REAL*8 (A-H,O-Z)
	
	integer, parameter          :: NCM = 30
	                            
	integer                     :: calc
	                            
	real(8), dimension(NCM)     :: X, FUGT, FUG, FUGP, DH_N
	real(8), dimension(NCM,NCM) :: FUGX
	COMMON /STYR2/                 model, N, NDER, NTEMP
! 	COMMON /OLDZ/                  ZOLD1,ZOLD2
	COMMON /OLDT/                  TOLD
!
!	Note that temperature derivatives are not evaluated. It
!	may be of interest for the user to evaluate the temperature
!	derivatives of the fugacity coefficients (ntemp=1), since
!	this permits a calculation of other residual properties di-
!	rectly, e.g. enthalpy and entropy
!
!	Se modificó nTemp
!
!
	NTEMP=0
	NDER=1
	calc = 1 ! only fugacity coefficients
	IF (INDIC > 2) then
	
	  NDER=2
	  calc = 5 !add derivatives wrt composition
	  
	endif	  
! 	TCR=DABS(TOLD-T)
	IGZ=0
	call THERMO (model, N, mTyp, calc, iGz, T, P, x, Z, FUG, FUGX, FUGT, FUGP, H, DH_N, Cp, B, ic)	
	
	RETURN

	ENDsubroutine SRK
!-----------------------------------------------------------------------
!
!
!	O N E  D I M E N S I O N A L  E I G E N V E C T O R  
!	S E A R C H   R O U T I N E
!	*********************************************************
!	This subroutine is a very powerfull tool for locating close
!	minima (i.e. critical behaviour), but the routine has not yet
!	been extensively tested. except for time consumption, it will
!	however never do any harm. it can be decoupled in a simple
!	manner
!
	SUBROUTINE ONEDIR (IOUT, N, RESV, YBEST, FX)
	
	IMPLICIT REAL*8 (A-H,O-Z)
	
	integer, parameter      :: NCM = 30, NFM = 4
	
	DIMENSION                  YBEST(NCM*NFM), FX(NCM,NCM), FX_null(NCM,NCM)
	real(8), dimension(NCM) :: Z, ZR, Y, YD1, YD2, YV, YD, FG, F1, FP, Fpres, yDer, yNew, Corr,            &
	                           yVal, GRD, GDG
	
	COMMON /STORE/             Z, ZR, Y, YD1, YD2, YV, YD, FG, F1, FP, YDER, YNEW, CORR, DIST(NFM), &
	                           YVAL, GRD, GDG, IP(NCM)
	COMMON /KEEP/              T,P,G0,ZFEED(NCM),FR(NCM),XVL(NCM,NFM),SF(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /MOLFRC/            ZFRAC(NCM,NFM)
	COMMON /SOLID/             FUGAR(NCM), ZL(NCM), Zcomp_pure(NCM), IS(NCM), ISTAT
!      
	RESV=0.
	NFF=NF
	IF (ISTAT == 1) NFF=NF-1
	IF (NFF == 0) RETURN
	DO 410  KKK=1,NFF
!
!	Set up phase composition, solid not tested
!
	  DO 20 I=1,N
	  
	    Z(I)=ZFRAC(I,KKK)
20	    ZR(I)=DSQRT(Z(I))
!
!	Evaluate partial derivatives
!
	  CALL SRK (3,0,IC,T,P,ZV,Z,FG,Fpres,FX)
!	  do io = 1,N
!	    write(*,'(/,<N>("| ",<N>(F15.10,2X), "|"),/)')(FX(io,ipo),ipo=1,N)
!	  enddo
	  DO 30 I=1,N
!
!	  Para evitar elementos nulos en la diagonal.
!	  Francisco, 05/08/2012
	    if (fx(i,i) == 0D0)fx(i,i) = 1e-16
!

	    F1(I)=FR(I)
	    RVAL=ZR(I)
	    DO 30 J=I,N
	    
	      FX(I,J)=RVAL*FX(I,J)*ZR(J)
30	      FX(J,I)=FX(I,J)
!
!	At this point, FX(i,j) = sqrt( z(i)*z(j) )*d[ ln(phiZ(i)) ]/dz(j)
!
!
!	...possible, for calculation of d2(tpd)/ds2 (pag 246), being FX the
!	M matriz.
!
!	But maybe it's equation 22 on page 236, though that one is for W moles 
!	(testing phase) and not Z (feed) as in page 246.
!
!
!	Estimate eigenvalue by direct iteration
!
	  XMAX=0.
	  DO 40 I=1,N

	    XMT=DABS(FX(I,I))
	    IF (XMT.LT.XMAX) GO TO 40

            XMAT=XMT	!search of the máximum diagonal element
            IMAX=I	!position of maximum element

40	    CONTINUE

	  SSQ=0.D0
	  DO 50 I=1,N
	  
	    YD(I)=FX(I,IMAX)
50	    SSQ = SSQ + YD(I)**2	!
!	  Ssq = sum of the squared of the raw with the maximum diag element

	  SSQS=DSQRT(SSQ)
	  DO 60 I=1,N
	  
60	    YD(I)=YD(I)/SSQS

	  DO 90 M=1,3
	  
	    PR1=0.D0
	    PR2=0.D0
	    DO 80 I=1,N
	    
	      TERM=0.D0
	      DO 70 J=1,N
	      
70	        TERM=TERM+FX(I,J)*YD(J)

	      PR1=PR1+TERM*YD(I)
	      PR2=PR2+TERM*TERM
80	      YV(I)=TERM

	    PR2S=DSQRT(PR2)
	    DO 90 I=1,N
	    
90	      YD(I)=YV(I)/PR2S


	  RAT=PR2/PR1
	  RAT1 = RAT + 1.D0
	  write (*, 100) KKK,RAT1

	  WRITE(IOUT,100) KKK,RAT1
100	  FORMAT(/,'  Phase investigated ',I2,' Eigenvalueestim.',F9.3)
!
!	Stop testing if estimate is not promising
!	otherwise refine by inverse iteration
!
	  IF (RAT > -.7) GO TO 410
	  
	    DO 110 I=1,N
	    
110	      FX(I,I)=FX(I,I)+1
!	  Definitely, FX is the hessian matriz, evaluated at Z.
!	  But YV (aka YD, aka YV/PR2S, aka a function of maximum diag element)
!	  doesn't seem to be the gradient.

	    CALL LU(NCM,N,IP,FX)
	    RATO=0.D0
	    DO 150 K=1,10
	    
	      DO 120 I=1,N
	      
120	        YV(I)=YD(I)
      
	      CALL BACK(NCM,N,IP,FX,YV)
	      SP1=0.D0
	      SP2=0.D0
	      DO 130 I=1,N
      
	        YV(I)=YV(I)-YD(I)
	        SP1=SP1+YV(I)*YD(I) 
130	        SP2=SP2+YV(I)**2

	      RAT=SP2/SP1  ! RAT seems to be an eigenvalue (eq 33 on page 260)
	                   ! (or an estimate...)
	      RTEST=1/(RAT+1)
	      SP2=DSQRT(SP2)
	      DO 140 I=1,N
	      
140	        YD(I)=YV(I)/SP2
!
!	      if (dabs(Rtest - Rat0) .LT. 1D-6)exit
	      if ( DABS(RTEST-RATO) .LT. 1.D-6) GO TO 160  !>--V
!	                                                      |
	      RATO=RTEST
!	                                                      |
150	    CONTINUE  !<--------------------------------------<
160	    RAT=1D0/(RAT+1D0)
!
!	Start linesearch
!
	    YVS=0.D0
	    YVSQ=0.D0
	    DO 170 I=1,N
	  
	      YV(I)=YD(I)/ZR(I)*.5
170	      YVSQ=YVSQ+YV(I)**2

	    YVS=0.D0
	    IPOS=0
!
!	Search for 'close' minimum limited to distance of 1
!	along eigenvector in direction of negative third derivative
!
	    INEG=0
	    FUNMIN=0.D0
	    GYOLD=1.D0
	    NCN=0
	    SMAX=0.98D0
!
!	Sequence terminates at 0.985, if nothing is found
!	Minimum distance 0.02 may miss a very,very close minimum
!
	    SNUL=.02D0	!Fran. 
!      
	    DS = 1.5*SNUL
	    IG=0
	    IST=0
	    S=SNUL
	    GO TO  190  !>--------------------V
!	                                      |
180	      SOLD=S  !<----------------------------------------------------<
!	                                      |                             |
	      FOLD=FUNM
	      GYOLD=GY
!	                                      |                             |
	      IF (S == SMAX) GO TO 370  !>----------V
!	                                      |     |                       |
	      S=S+DS
	      IF (S > SMAX) S=SMAX
	      IG=1
!	                                      |     |                       |
	      IF (S == 0) GO TO 380  !>-------------------V
!	                                      |     |     |                 |
!	    There's a GOTO loop to which begins here...
190	    LOC=0   !<------------------------<-----------------<           |
!	                                            |     |     |           |
!	To check when calculation and refinement are of opposite sign
195	      LOC=LOC+1  !<-------------------------------------------<     |
!	                                            |     |     |     |     |
	      SYVAL=0.
	      IF (IG == 0) IST=IST+1
	      DO 210 I=1,N
!	                                            |     |     |     |     |
	        SY=S*YV(I)		
	        IF (SY.LT.0.) GO TO 200  !>---V     |     |     |     |     |
!	                                      |     |     |     |     |     |
	          YNEW(I)=(ZR(I)*(1+SY))
!	          exit                        |     |     |     |     |     |
	          GO TO 210 !>----------V     V     V     V     ^     ^     ^
!	                                |     |     |     |     |     |     |
!	Avoid negative mole fractions at large distances
200	        FAC=(1-SY)*(1+SY*SY)  !<------<     |     |     |     |     |
	        YNEW(I)=ZR(I)/FAC
!	                                |           |     |
210	      CONTINUE  !<--------------<           |     |     |     |     |
	      DO 220 I=1,N
!	                                            |     |     |     |     |
	        IF (IG == 0) GO TO 220  !>----V     |     |     |     |     |
!	                                      |     |     |     |     |     |
!	Corr evaluates a quadratic correction to search direction
	          TTT=S*S*CORR(I)
	          YN=YNEW(I)+TTT
	          YNEW(I)=YN
!	                                      V     V     V     ^     ^     ^
!	Replace very small quantities, using direct substitution
	          IF (YN.LE.0. .OR. DABS(YN/TTT).LT..5) then
	            YNEW(I)=DEXP((F1(I)-FP(I))/2)
	          endif
!	                                      |     |     |     |     |     |
220	          YVAL(I)=YNEW(I)**2  !<------<     |     |     |     |     |
!	                                            |     |     |     |     |
	      DO 221 L=1,NF
!	                                            |     |     |     |     |
	        DIS=0.
	        DO 222 I=1,N
!	                                            |     |     |     |     |
222	          DIS=DIS + (YVAL(I)-ZFRAC(I,L))**2
!	                                            V     V     ^     ^     ^
221	        DIST(L)=DIS
!	                                            |     |     |     |     |
!	Test distance in composition space to solutions located allready
!	limit search to avoid overlap
!	                                            |     |     |     |     |
	      DO 223 L=1,NF
!	                                            |     |     |     |     |
	        IF (DIST(L)/DIST(KKK) .LT. .5) SMAX=S
!	                                            |     |     |     |     |
223	      CONTINUE
! 	      CALL SRK (1,0,IC,T,P,ZV,YVAL,FP,FP,FX)
	      CALL SRK (1,0,IC,T,P,ZV,YVAL,FP,Fpres,FX_null)
	      XPR=0.
	      FUN=1.
	      YDP=0.
	      GDP=0.
	      DO 250 I=1,N
!	                                            |     |     |     |     |
	        DDS=F1(I)-FP(I)
	        IF (DDS >  -14.) GO TO 230!>-V
!	                                      |     |     |     |     |     |
	          YNEW(I)=DEXP(DDS/2)
	          YVAL(I)=YNEW(I)**2
	          TM=0
!	          exit                        |     |     |     |     |     |
	          GO TO 240 !>----------V     |     |     |     |     |     |
!	                                |     |     |     |     |     |     |
230	        TM=DLOG(YVAL(I))-DDS  !<------<     |     |     |     |     |
!	                                |           |     |     |     |     |
240	        FUN=FUN+YVAL(I)*(TM-1) !<           V     V     ^     ^     ^
!	                                            |     |     |     |     |
	        TM=DLOG(YVAL(I))+FP(I)-F1(I)	!Modified tpd
	        QQ=YNEW(I)*TM
	        GRD(I)=QQ
	        GDG(I)=QQ
250	        XPR=XPR+YD(I)*QQ
!	                                            |     |     |     |     |
	      GY=XPR/S/RAT
!	                                            |     |     |     |     |
!	Investigate sign of third derivative at origo
!	change direction, if sign is wrong
!	                                            |     |     |     |     |
	    TN=(GY-1)*RAT/S
	    IF (IG == 1) GO TO 270  !>--------V     V     V     ^     ^     ^
!	                                      |     |     |     |     |     |
	    IF (TN*S .LT. 0.) GO TO 270 !>----V     |     |     |     |     |
!	                                      |     |     |     |     |     |
	      IF (IST.GE.2) GO TO 270 !>------V     |     |     |     |     |
!	                                      |     |     |     |     |     |
	        DO 260 I=1,N
!	                                      |     |     |     |     |     |
	          YV(I)=-YV(I)
260	          YD(I)=-YD(I)
!	                                      |     |     |     |     |     |
	  GO TO 190  !>-----------------------------------------^     |     |
!	                                      |     |     |           |     |
270	  CONTINUE  !<------------------------<     |     |           |     |
!	                                            |     |           |     |
	  DO 280 I=1,N
!	                                            |     |           |     |
280	    GDG(I)=GDG(I)-XPR*YD(I)
!	                                            V     V           ^     ^
	  CALL BACK(NCM,N,IP,FX,GDG)
	  GAIN=0.D0
	  DO 290 I=1,N
!	                                            |     |           |     |
290	    GAIN=GAIN+GDG(I)*(GRD(I)+XPR*YD(I))
!	                                            |     |           |     |
	  GAIN=GAIN/2
!	                                            |     |           |     |
!	Include correction term for next step
!	                                            |     |           |     |
	  DO 300 I=1,N
!	                                            |     |           |     |
	    CRN=YNEW(I)-.5*GDG(I)
	    IF (IG == 0 ) CORR(I)=-GDG(I)/SNUL**2/2
	    IF (IG > 0) CORR(I)=CORR(I)-GDG(I)/S/S/2
!	                                            |     |           |     |
300	  CONTINUE
	  FUNM=FUN-GAIN
!	                                            V     V           ^     ^
!	Only accept verified, not suggested results, but repeat
!	the calculation once with the corrected value, if the sugges-
!	ted result is negative
!	                                            |     |           |     |
        IF (LOC == 1 .AND. FUNM*FUN .LT. 0.) GO TO 195 !>-----------^
!	                                            |     |                 |
	  IF (FUN.GE.FUNMIN) GO TO 320 !>-----V
!	                                      |     |     |                 |
	    FUNMIN=FUN
	    DO 310 I=1,N
!	                                      V     V     V                 ^
310	      YBEST(I)=YVAL(I)
!	                                      |     |     |                 |
320	  CONTINUE  !<------------------------<
!	GY is the procected gradient,which determines the search
	  IF (RAT.LT.0.) GY=-GY
!	                                            |     |                 |      
	  IF (IG == 0) GO TO 180  !-----------------------------------------^
!	                                            |     |                 | 
	  IF (INEG == 1) GO TO 340  !>--------V     |     |                 | 
!	                                      |     |     |                 | 
	    IF (GY.LT.0.) GO TO 330 !>--V     V     V     V                 ^
!	                                |     |     |     |                 | 
	      IF (GY > GYOLD) GO TO 380  !>--------------V                 | 
!	                                |     |     |     |                 | 
	        DS=1.5*DS
	        GO TO 180
!	                                |     |     |     |                 | 
330	    INEG=1  !<------------------<-----------------------<                     
!	   3                                  |     |     |     |           | 
	    GYNEG=GY
	    SNEG=S
	    DS=1.5*DS
	    FNEG=FUNM
!	                                      V     V     V     ^           ^
	    GO TO 180
!	                                      |     |     |     |           | 
!340	  IF (GY.LT.0. .AND. IPOS == 0) GO TO 330
340	  continue !<-------------------------<     |     |     |           | 
	  IF (GY.LT.0. .AND. IPOS == 0) GO TO 330 !>------------^           | 
	  IPOS=1
!	                                            |     |                 | 
!	Stationary point bracketed; only 3 further steps taken.
!	                                            |     |                 | 
	  NCN=NCN+1
	  IF (NCN > 3) GO TO 380  !>---------------------V                 | 
!	                                            |     |                 | 
	    IF (GY > 0.) GO TO 350 !>--------V     |     |                 | 
!	                                      |     |     |                 | 
	      GYNEG=GY
	      SNEG=S
	      FNEG=FUNM
!	                                      |     |     |                 | 
	      GO TO 360  !>-------------------------------------V           | 
!	                                      |     |     |     |           | 
350	    GYPOS=GY  !<----------------------<     |     |     |           | 
	    SPOS=S
	    FPOS=FUNM
!	                                            |     |     |           | 
360	    RV=(GYPOS-GYNEG)/(SPOS-SNEG)  !<--------------------<           | 
!	Inverse linear interpolation used in search
	    S=SNEG-GYNEG/RV
	    DS=0
!	                                            |     |                 | 
	    GO TO 180  !>---------------------------------------------------^
!	                                            |     |
370	    CONTINUE  !<----------------------------<
!	                                                  |
380	  CONTINUE  !<------------------------------------<
	  IF (FUNMIN .GE. 0.) GO TO 410  !>---------------------V
!	                                                        |
	    DO 400 I=1,N
!	                                                        |
400	      YBEST(I)=DLOG(YBEST(I))
!	                                                        |
	    write (*, 390 ) FUNMIN
	    WRITE(IOUT,390 ) FUNMIN
390	    FORMAT(/,'  Search conclusive, FMIN = ',D10.2)
	    RESV=FUNMIN
!	                                                        |                                                        |
!	    IF (FUNMIN.LE.-1.D-5) RETURN		!alfo?
!	    IF (FUNMIN.LE.-1.D-7) RETURN
	    IF (FUNMIN.LE.-1.D-10) RETURN		!fran
!	                                                        |    
410	CONTINUE  !<--------------------------------------------<

	RETURN
	END
!-----------------------------------------------------------------------
!
!
!	PITON evaluates phase distribution, given the fugacity coef-
!	ficients, including search for appearing or vanishing solid
!	and fluid phases
!
	SUBROUTINE PITON (N,H,VEC,ALFA,YVAL,XJAC,IPX)
	IMPLICIT REAL*8 (A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4
	DIMENSION             ALFA(NCM,NFM), BETEX(NCM), IPX(NFM+2),G(NFM+2),H(NCM),XJAC(NFM+2,NFM+2)
	DIMENSION             YVAL(NCM*NFM),VEC(NCM*(NFM-1)),BOLD(NFM+2),BTOLD(NCM)
	LOGICAL               NOSOL
	COMMON /SOLID/        FUGAR(NCM),ZLOG(NCM), Zcomp_pure(NCM), IS(NCM),ISTAT
	COMMON /KEEP/         T,P,GNUL,Z(NCM),FR(NCM),XVL(NCM,NFM),BETA(NFM),ZV(NFM),indF(NFM),NF,ND

!	CALCULATION OF PHASE DISTRIBUTION, SOLID PHASES
      IOUT=0
      ITER=0
      NOSOL=(ISTAT == 0)
      NR=NF-1
      IF (NOSOL) NR=NF
      IF (NOSOL) GOTO 205
!	-  K-FACTORS, SOLIDS PRESENT
      IF (NR > 1) GOTO 205
!	-  ONLY SINGLE PHASE;
      DO 202 I=1,N
      DD=0.
      ALFA(I,1)=1.
      IF (IS(I) == 0) GOTO 202
      DIF=VEC(I)
      IF (DIF.LT. -60.) DIF=-60.
      DD=DEXP(DIF)
  202 ALFA(I,2)=DD
      GOTO 300
  205 NG=NR-1
      DO 210 I=1,N
      SMAX=0.
      IP=I-N
      DO 212 J=1,NG
      IP=IP+N
      SS=VEC(IP)
      IF (SS.LE. SMAX) GOTO 212
      SMAX=SS
  212 CONTINUE
      IP=I-N
      DO 214 J=1,NR
      IP=IP+N
      IF (J.LT.NR) SS=VEC(IP)
      IF (J == NR) SS=0.
      SS=SS-SMAX
      IF (SS.LT.-60.) SS=-60.
      IF (SS == 0.) ALFA(I,J)=1.
      IF (SS /= 0.) ALFA(I,J)=DEXP(SS)
  214 CONTINUE
      IF (NOSOL) GOTO 210
      IF (IS(I) > 0.) GOTO 216
      ALFA(I,NF)=0.
      GOTO 210
  216 SS=VEC(I+(NR-1)*N)
      SS=SS-SMAX
      IF (SS.LT. -60.) SS=-60.
      ALFA(I,NF)=DEXP(SS)
  210 CONTINUE
!	-START PHASE ITERATION
  300 IF (NOSOL) GOTO 100
!	--EVALUATE CURRENT STATE OF SOLIDS
      DO 310 I=1,N
      IF (IS(I) == 0) GOTO 310
      HI=0.
      DO 320 J=1,NR
  320 HI=HI+ALFA(I,J)*BETA(J)
      BETX=Z(I)-HI/ALFA(I,NF)
      IS(I)=1
      IF (BETX > 0.) IS(I)=2
      BTOLD(I)=BETX
  310 CONTINUE
  100 ITER=ITER+1
      DO 10 I=1,N
      HI=0.
      IF (IS(I).LT.2 .OR. NOSOL) GOTO 15
!	--MOLE FRACTIONS FIXED BY SOLID EQUIL.
      HI=Z(I)*ALFA(I,NF)
      GOTO 10
   15 DO 20 J=1,NR
   20 HI=HI+ALFA(I,J)*BETA(J)
      IF (HI == 0.) HI=1.
   10 H(I)=1/HI
!	-CONVERGENCE TEST
      IF (IOUT == 1) GOTO 350
!	INITIATE JACOBIAN
      DO 30 I=1,NR
      G(I)=1.
      DO 30 J=1,I
   30 XJAC(I,J)=0.
!	MAIN LOOP, GRADIENT AND JACOBIAN
      DO 40 I=1,N
      ZI=Z(I)
      HI=H(I)
      IST=IS(I)
      DO 45 J=1,NR
      ZL=ALFA(I,J)*HI
      PR=ZL*ZI
      G(J)=G(J)-PR
      IF (BETA(J) == 0.) GOTO 45
      IF (IST == 2 .AND. ISTAT > 0) GOTO 45
      DO 50 K=1,J
      ZK=ALFA(I,K)*HI
   50 XJAC(J,K)=XJAC(J,K)+ZK*PR
   45 CONTINUE
   40 CONTINUE
!	-CHECK GRADIENT AND COPY JACOBIAN
      DO 60 J=1,NR
      IF (BETA(J) > 0.) GOTO 70
      IF (G(J).LT.0.) GOTO 65
!	GRADIENT POSITIVE; SKIP OKAY
      DO 61 K=1,NF
      XJAC(J,K)=0.
   61 XJAC(K,J)=0.
      XJAC(J,J)=1.
      G(J)=0.
      GOTO 70
!	-RE-INITIATE PHASE, GRADIENT NEGATIVE
   65 BETA(J)=1.D-16
      GOTO 100
!	-SAFEGUARD AGAINST ZERO DETERMINANT OF JACOBIAN
   70 XJAC(J,J)=XJAC(J,J)+1.D-7
      DO 60 K=1,J
   60 XJAC(K,J)=XJAC(J,K)
!	- CALCULATE NEWTON STEP
      CALL LU (NFM+2, NR, IPX, XJAC)
      CALL BACK (NFM+2, NR, IPX, XJAC, G)
!	--MONITOR STEPLENGTH TO AVOID NEGATIVE PHASE FRACTIONS
      FACL=0.
      FACR=1.
      JZM=0
      DO 80 J=1,NR
      DA=G(J)
      DAB=DABS(DA)
      IF (DAB > FACL) FACL=DAB
      GNEW=BETA(J)-DA
      IF (GNEW.GE.0.) GOTO 80
      RAT=BETA(J)/DAB
      IF (RAT > FACR) GOTO 80
      FACR=RAT
      JZM=J
   80 CONTINUE
!	-EVALUATE REDUCTION FACTOR
      RANG=5*FACL
      DEC=1.
      IF (RANG > 1.) DEC=1/RANG
      IF (JZM == 0) GOTO 90
      IF (FACR.LT. DEC) GOTO 85
      JZM=0
      GOTO 90
   85 DEC=FACR
   90 DO 110 J=1,NR
      BOLD(J)=BETA(J)
      G(J)=G(J)*DEC
      BETA(J)=BETA(J)-G(J)
      IF (J == JZM) BETA(J)=0.
  110 CONTINUE
!	-END FLUID PHASES; CHECK STATUS OF SOLIDS; FIND REDUCTION
!	-CORRESPONDING TO ONE CHANGE
  350 IF (NOSOL) GOTO 400
      ICHAN=0
      RMIN=1.
      DO 120 I=1,N
      IST=IS(I)
      IF (IST == 0) GOTO 120
      HI=0.
      DO 125 J=1,NR
  125 HI=HI+BETA(J)*ALFA(I,J)
      BETX=Z(I)-HI/ALFA(I,NF)
      IF (IOUT == 1) GOTO 130
      BOLDI=BTOLD(I)
      IF (BETX*BOLDI  > 0.) GOTO 130
      ICHAN=ICHAN+1
      DIF=BETX-BOLDI
      RAT=DABS(BOLDI/DIF)
      IF (RAT.LE. RMIN) RMIN=RAT
  130 BETEX(I)=BETX
  120 CONTINUE
      IF (IOUT == 1) GOTO 400
      IF (ICHAN > 0) GOTO 140
      DO 135 I=1,N
      IF (IS(I) == 0) GOTO 135
      BTOLD(I)=BETEX(I)
  135 CONTINUE
      GOTO 400
  140 RT1=1.01*RMIN
      IF (RT1 >  .99) RT1=(1+RAT)/2
      RT2=1-RT1
      DO 141 J=1,NR
  141 BETA(J)=RT1*BETA(J)+RT2*BOLD(J)
      DO 145 I=1,N
      IF (IS(I) == 0) GOTO 145
      BETEX(I)=RT1*BETEX(I)+RT2*BTOLD(I)
      IF (BETEX(I)*BTOLD(I) .LT.0.) IS(I)=3-IS(I)
      BTOLD(I)=BETEX(I)
  145 CONTINUE
  400 IF (IOUT == 1) GOTO 150
      IF (ITER.GE.20 .OR. FACL.LT. 1.D-8) IOUT=1
      GOTO 100
!	CONVERGENCE;
  150 IF (ISTAT > 0.) BETA(NF)=0.
      DO 160 I=1,N
      HI=H(I)
      DO 165 J=1,NR
      BET=BETA(J)
      IF (BET.LE.0.) BET=-1.
  165 XVL(I,J)=BET*ALFA(I,J)*HI
      IF (NOSOL) GOTO 160
      IF (IS(I).LT.2) XVL(I,NF)=0.
      IF (IS(I) == 2) XVL(I,NF)=BETEX(I)/Z(I)
      BETA(NF)=BETA(NF)+XVL(I,NF)
  160 CONTINUE
      RETURN
      END
!	------------------------------------------------------------------
!
!
!       PHASN EVALUATES FUNCTION, GRADIENT AND HESSIAN FOR PHASE SPLIT
!       CALCULATIONS. PHASN IS AN ARGUMENT OF THE MURRAY ROUTINE
!
	SUBROUTINE PHASN(NARG,IVAR,FUN,GRAD,XMAT,YVAL)
	
	IMPLICIT REAL*8 (A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4
	DIMENSION             INDEX(NCM),NCON(NCM,NFM),ENER(NFM)
	DIMENSION             YVAL(NCM*NFM),GRAD(NCM*(NFM-1)),X(NCM),TM(NCM),FG(NCM),FP(NCM)
	DIMENSION             XMAT(NCM*(NFM-1),NCM*(NFM-1)),XMIN(NCM),DELEN(NFM)
	COMMON /MOLFRC/       FILD(NCM,NFM)
	COMMON /DERIV/        FX(NCM,NCM)
	COMMON /KEEP/         T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /SOLID/        FUGAR(NCM),ZLOG(NCM),Zcomp_pure(NCM),ISol(NCM),IStat
	COMMON /OPTIM1/       VBEST(NCM*NFM),GMIN
	COMMON /AKTIV/        IVAL(NCM),IPOS(NCM*(NFM-1)),NACT
	COMMON /LIMIT/        EPSV
	common /OUTUNIT/      IOUT	

      JD=1
      N=NARG/(NF-1)
   10 NG=NF-1
      NFE=NF-ISTAT
      NARG=N*NG
      if (IVAR.LT.2) GOTO 110
      JD=3
      SMIN=1.
!
!     SECOND DERIVATIVES NEEDED.
!     IDENTIFY TRACE COMPONENTS AND DEPENDENT VARIABLES
!
      DO 20 J=1,NF
      IF (SFAS(J) > SMIN) GO TO 20
      SMIN=SFAS(J)
      JMIN=J
   20 CONTINUE
!
!     ARRAY NCON: 1: INDEPENDENT; 0: DEPENDENT; -1: TRACE
!     INDEX KEEPS TRACE OF DEPENDENT VARIABLES
!
      DO 70 I=1,N
      IVAL(I)=1
      XMAX=0
      XM1=0.
      DO 30 J=1,NF
      IV=I+(J-1)*N
      YVAL(IV)=VBEST(IV)
      XVL(I,J)=VBEST(IV)
      XM=YVAL(IV)*ZZ(I)/SFAS(J)
      if (XM > XM1)XM1=XM
      NCON(I,J)=1
      if (XM.LT.EPSV)NCON(I,J)=-1
      if (J == JMIN)XMIN(I)=XM
      if (YVAL(IV).LT.XMAX) GOTO 30
      XMAX=YVAL(IV)
      JMAX=J
   30 CONTINUE
      INDEX(I)=JMAX
      NCON(I,JMAX)=0
      IG=0
      IF (XM1 > EPSV) GOTO 40
      IVAL(I)=0
      NCON(I,JMAX)=-1
      GO TO 70
   40 DO 60 J=1,NF
      IF (NCON(I,J).GE.0) GO TO 50
      IF (XVL(I,J)*ZZ(I)/SFAS(JMAX) >  EPSV) NCON(I,J)=1
   50 IF (NCON(I,J) >  0) IG=1
   60 CONTINUE
      IF (IG == 0) IVAL(I)=0
   70 CONTINUE
!
!     ORDER ACTUAL INDEPENDENT VARIABLES
!
      DO 80 I=1,N
      IV=INDEX(I)
      S=YVAL(I+(IV-1)*N)
      YVAL(I+(IV-1)*N)=YVAL(I+(NF-1)*N)
   80 YVAL(I+(NF-1)*N)=S
      NACT=0
      IPOS(1)=1
      DO 100 J=1,NG
      JN=(J-1)*N
      DO 90 I=1,N
      IF (IVAL(I) == 0) GOTO 90
      IJ=I+JN
!	 NOTE CHANGE FROM PREVIOUS INCORRECT VERSION
      JJ=J
      IF (J == INDEX(I)) JJ=NF
      IF (NCON(I,JJ).LT.0) GOTO 90
      NACT=NACT+1
      IPOS(NACT)=IJ
   90 CONTINUE
  100 CONTINUE
      IF (NACT == 0) NACT=1
      GOTO 200
!
!     REARRANGE VARIABLES
!
  110 DO 120 I=1,N
      IV=INDEX(I)
      S=YVAL(I+(IV-1)*N)
      YVAL(I+(IV-1)*N)=YVAL(I+(NF-1)*N)
  120 YVAL(I+(NF-1)*N)=S
      DO 190 I=1,N
      IV=INDEX(I)
      XVLV=1.
      RSM=0.
      DO 130 J=1,NF
      if (NCON(I,J) == -1) GOTO 140
  130 CONTINUE
      GOTO 160
 140  BFAS=SFAS(IV)
      IF (IV > NFE) BFAS=XVL(I,NF)*ZZ(I)
      DO 150 J=1,NF
      TMJ=1
      IF (J == IV) GOTO 155
      SFJ=SFAS(J)
      IF (J > NFE) SFJ=XVL(I,J)
      IF (J.LE.NFE .OR. ISOL(I) == 2) GOTO 152
      TMJ=0
      GOTO 155
  152 DSN=FILD(I,IV)-FILD(I,J)
      IF (DSN.LT. -80.) DSN=-80.
      TMJ=SFJ/BFAS*DEXP(DSN)
  155 TM(J)=TMJ
  150 RSM=RSM+TMJ
  160 SM=0.
!
!     CORRECT TRACE COMPONENTS, USING OLD K-FACTORS AND PHASE AMOUNTS
!
      DO 170 J=1,NF
      IF (J == IV) GOTO 170
      XVL(I,J)=YVAL(I+(J-1)*N)
      IF (XVL(I,J) == 0.) GOTO 175
      if (NCON(I,J) == -1)XVL(I,J)=TM(J)/RSM
  175 XVLV=XVLV-XVL(I,J)
  170 CONTINUE
      XVL(I,IV)=XVLV
      DO 180 J=1,NF
      SLIMT=1.D-8
      IF (J == NF .AND. ISTAT == 1) SLIMT=0.
      IF (XVL(I,J).LT.0.) XVL(I,J)=SLIMT
  180 SM=SM+XVL(I,J)
      DO 190 J=1,NF
      XVL(I,J)=XVL(I,J)/SM
  190 YVAL(I+(J-1)*N)=XVL(I,J)
!	-INITIATE OBJECTIVE FUNCTION
  200 FUN=-GNUL
      IF (IVAR == 0) GOTO 230
!	-INITIATE GRADIENT

	write(*,*)NARG


      DO 220 I=1,NARG
      GRAD(I)=0.
      IF (IVAR.LT.2) GOTO 220
!	-INITIATE HESSIAN
      DO 210 J=1,NARG
  210 XMAT(I,J)=0.
  220 CONTINUE
!
!     START CALCULATING CONTRIBUTION FROM INDIVIDUAL PHASES
!
  230 DO 350 K=1,NF
      SFAS(K)=0.
      ENER(K)=0.
      DO 240 I=1,N
      X(I)=XVL(I,K)*ZZ(I)
  240 SFAS(K)=SFAS(K)+X(I)
      IF (K == NF .AND. ISTAT == 1) GOTO 245
      CALL SRK (JD,0,IC,T,P,ZVAL(K),X,FG,FP,FX)
!
	indF(k) = iC
!
  245 DELEN(K)=0
      DO 250 I=1,N
      IF (K /= NF .OR. ISTAT /= 1) GOTO 255
      IF (ISOL(I).LT. 2) FILD(I,K)=100.
      IF (ISOL(I) == 2) FILD(I,K)=FUGAR(I)
      TM(I)=0.
      TMI=0
      IF (ISOL(I) == 2) TMI=FUGAR(I)
      TMI=TMI-ZLOG(I)
      TM(I)=TMI
      GOTO 248
  255 FILD(I,K)=FG(I)
      TM(I)=FG(I)+DLOG(XVL(I,K)/SFAS(K))
      TMI=TM(I)
  248 DELEN(K)=DELEN(K)+XMIN(I)*TMI
  250 ENER(K)=ENER(K)+TMI*X(I)
      FUN=FUN+ENER(K)
!
!     END GIBBS ENERGY, START GRADIENT
!
      IF (IVAR == 0) GOTO 350
      DO 280 I=1,N
      S=TM(I)*ZZ(I)
      IPL=K
      IV=INDEX(I)
      IF (K == NF) IPL=IV
      IF (K == IV) GOTO 260
      IPS=I+(IPL-1)*N
      GRAD(IPS)=GRAD(IPS)+S
      GOTO 280
  260 DO 270 M=1,NG
      IPS=I+(M-1)*N
  270 GRAD(IPS)=GRAD(IPS)-S
  280 CONTINUE
!
!     END GRADIENT, START HESSIAN MATRIX
!
      IF (IVAR.LT.2) GOTO 350
      IF (K == NF .AND. ISTAT == 1) GOTO 350
      DO 300 I=1,N
      S=ZZ(I)/SFAS(K)
      DO 290 J=1,I
      FX(I,J)=S*(FX(I,J)-1.)*ZZ(J)
  290 FX(J,I)=FX(I,J)
  300 FX(I,I)=FX(I,I)+ZZ(I)/XVL(I,K)
      DO 340 I=1,N
      IV=INDEX(I)
      IPL=K
      IF (K == NF) IPL=IV
      DO 330 J=1,N
      JV=INDEX(J)
      JPL=K
      IF (K == NF) JPL=JV
      S=FX(I,J)
      M11=IPL
      M12=IPL
      M21=JPL
      M22=JPL
      IF (IV /= K) GOTO 310
      S=-S
      M11=1
      M12=NG
  310 IF (JV /= K) GOTO 320
      M21=1
      M22=NG
      S=-S
  320 DO 330 M1=M11,M12
      DO 330 M2=M21,M22
      IPS=I+(M1-1)*N
      JPS=J+(M2-1)*N
  330 XMAT(IPS,JPS)=XMAT(IPS,JPS)+S
  340 CONTINUE
  350 CONTINUE
      IF (IVAR /= 1) GOTO 370
!
!     CHANGE BACK TO NATURAL ORDER
!
      DO 360 I=1,N
      IV=INDEX(I)
      S=YVAL(I+(IV-1)*N)
      YVAL(I+(IV-1)*N)=YVAL(I+(NF-1)*N)
      YVAL(I+(NF-1)*N)=S
  360 CONTINUE
  370 CONTINUE
      IF (IVAR.LT.2.OR.NF.LE.2) GOTO 440
!
!     CHECK FOR REMOVAL OF A PHASE
!
      DMIN=1.
      DO 380 J=1,NF
      IF (J == JMIN) GOTO 380
      DELT=DELEN(J)-DELEN(JMIN)
      IF (DELT > DMIN) GOTO 380
      K=J
      DMIN=DELT
  380 CONTINUE
      IF (DMIN == 1.) GOTO 440
      J=JMIN
      DO 390 I=1,N
  390 X(I)=ZZ(I)*(XVL(I,J)+XVL(I,K))
      SFT=SFAS(J)+SFAS(K)
      CALL SRK (1,0,IC,T,P,ZNEW,X,FG,FP,FX)
      ENGO=ENER(J)+ENER(K)
      DO 400 I=1,N
  400 ENGO=ENGO-X(I)*(FG(I)+DLOG(X(I)/SFT/ZZ(I)))
      IF (ENGO.LT.-1.D-10) GOTO 440
      IF (J == NF) ISTAT=0
      IF (K.LT.J) GOTO 410
      JJJ=J
      J=K
      K=JJJ
  410 CONTINUE
!
!     STORE BEST RESULT SO FAR IN VBEST
!
      DO 420 I=1,N
      SS=XVL(I,K)+XVL(I,J)
      XVL(I,K)=SS
      XVL(I,J)=XVL(I,NF)
      VBEST(I+(K-1)*N)=SS
  420 VBEST(I+(J-1)*N)=VBEST(I+(NF-1)*N)
      NF=NF-1
      write (*, 430)
! 	iOut = 0
      WRITE(IOUT,430)
  430 FORMAT(/,'  Number of phases reduced ')
      GOTO 10
  440 IF (FUN > GMIN+1.D-10) GOTO 460
!
!     NO REDUCTION, STORE BEST RESULT
!
      GMIN=FUN
      DO 450 I=1,N
      DO 450 J=1,NF
      IV=I+(J-1)*N
  450 VBEST(IV)=XVL(I,J)
  460 RETURN
      END
!-----------------------------------------------------------------------
!
!	DSUB calculates phase split by direct substitution, using
!	the general dominant eigenvalue method of Crowe and Nishio
!	for acceleration
!
	SUBROUTINE DSUB(N,FX,YVAL)
	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4
	DIMENSION             YVAL(NCM*NFM), FX(NCM,NCM), IP(NFM+2)
!
!     Common area store is used for internal variables in the DSUB
!     routines to save space
!
	COMMON /KEEP/         T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /STORE/        GRAD(NCM*(NFM-1)),WORK(NCM*(NFM-1),NFM+1),SFOLD(NFM),B(NFM+2,NFM+2),  &
	                      Y(NCM*(NFM-1)),H(NCM),FG(NCM),YOLD(NCM*(NFM-1)),YNEW(NCM*(NFM-1)), &
	                      CR(NCM,NFM),BAU(NFM+2),YDAR(NCM,NFM)
	COMMON /OPTIM1/       VBEST(NCM*NFM),GMIN
	COMMON /NTEST/        GNORM,GLIM, NFA,NDUM
	COMMON /SOLID/        FUGAR(NCM),ZLOG(NCM),Zcomp_pure(NCM), IS(NCM),ISTAT
	COMMON /NEWSOL/       NEWVAL
	common /OUTUNIT/      iOut
!
!	Column dimension of CR must be increased to MAFX-1 with more
!     than 4 phases
!
!     GDEM used is of order M-2, where M must be < 6, here 4
!
! 	iOut = 0
	M=4
	LMAX=4
      
5	NFA=NF  !<-----------------------------------------------------<
!	                                                               |    
	NG=NF-1
	IF (NG > 0) GOTO 6  !>-----V                                  |
!	                            |                                  |
	  GNORM=0.D0
!	                            |                                  |
!	                            |                                  |
	  GOTO 440	
!	                            |                                  |
6	IV=1  !<--------------------<                                  |
	NT=N*NG
!	                                                               |
!	Subroutine NTOPI finds fugacity coefficients from compositions.
!	Subroutine PITON finds phase amounts (and composition) from
!	K-factors
!	   
	GMIN = 0                                                        
	CALL NTOPI(N,GTOT,FG,H,FX,YOLD,YDAR)
	GMIN=GTOT
	GBAS=GMIN
	DO 20 I=1,NT
!	                                                               |
20	  YNEW(I)=YOLD(I)
!	                                                               |
	DO 230 L=1,LMAX
!	                                                               |
!	ISO AND ICOR ARE FAILURE COUNTERS
!	                                                               |
	  ISO=0
	  ICOR=0
!	                                                               |	  
40	  CONTINUE  !<-----------------------------<                   |
!	                                           |                   |
	  ISO=ISO+1
	  M2=M+IV
!	                                           |                   |
!	  Following loop has GOTOs which exits to above lines...
	  DO 140 KK=1,M2
!	                                           |                   |
	    KV=KK-IV
	    DO 50 J=1,NF
!	                                           |                   |
50	      SFOLD(J)=SFAS(J)
!	                                           |                   |
	    CALL PITON(N,H,YOLD,CR,YVAL,B,IP)
!	                                           |                   |
	    GOTO 100  !>--------------------------------V
!	                                           |    |              |
!	    This loop is performed until KK<>1 or Gmin<1D-12
70	    DO 80 J=1,NF  !<---------------------------------<
!	                                           |    |    |         |
80	      SFAS(J)=SFOLD(J)
!	                                           |    |    |         |
	    DO 90 I=1,NT
!	      Restart from last value, if extrapolate is worse
	      IF (ISO > 2.OR.L == LMAX) YOLD(I)=YNEW(I)
90	      YOLD(I)=(YNEW(I)+YOLD(I))/2
!	                                           |    |    |         |
	    ICOR=ICOR+1
!	                                           |    |    |         |
	    IF (L == LMAX.OR.ICOR > 4) GOTO 230 !>---------------V
!	                                           |    |    |    |    |
	      GOTO 40  !>--------------------------^    |    |    |    |
!	                                                |    |    |    |
100	      CONTINUE  !<------------------------------<    |    |    |
!	                                                     |    |    |
	      CALL NTOPI(N,GTOT,FG,H,FX,GRAD,YDAR)
!	                                                     |    |    |
!	      Check for improvement
!	                                                     |    |    |
	      if (KK /= 1 .OR. GTOT.LE.GMIN+1.D-12) GOTO 110
!	                                                     |    |    |
	        GOTO 70  !>----------------------------------^    |    |
!	                                                          |    |
110	      IF (GTOT > GMIN+1.D-12) GOTO 130  !>----------V
!	                                                     |    |    |
	        GMIN=GTOT
	        DO 120 I=1,N
	          DO 120 J=1,NF
120	            YVAL(I+(J-1)*N)=XVL(I,J)
!	                                                     |    |    |
!	        Test for solids after extrapolation
	        IF (KK /= 1 .OR. ISTAT == 1) GOTO 125 !>V
!	                                                |    |    |    |
	        IF (NEWVAL == 0) GOTO 125 !>------------V    |    |    |
!	                                                |    |    |    |
!	          Solid phase found
	          write (*, 126)
	          WRITE(IOUT,126)
126	          FORMAT(/,'  Solid phase added in DSUB ',/)
!	                                                |    |    |    | 
	          NF=NF+1
	          ISTAT=1
	          SFAS(NF)=0.
	          indF(NF) = 0
	          DO 128 I=1,N
!	                                                |    |    |    | 	          
128	            XVL(I,NF)=0.
!	                                                |    |    |    |
	          GOTO 5  !--------------------------------------------^
!	                                                |    |    |    |   
125	        CONTINUE  !<----------------------------<    |    |    |
	        IRFAS=0
!	                                                     |    |    |        
115	        JOUT=0  !<----------------------------<
!	                                              |      |    |    |
	          DO 111 J=1,NF
!	A value in this range indicates that the phase can be removed
	            IF (SFAS(J) == 0.) JOUT=J
!	                                              |      |    |    |
111	          CONTINUE
	          IF (JOUT /= 0.) IRFAS=1
!	                                              |      |    |    |
!	          Exit from loop of line 115
	          IF (IRFAS == 0) GO TO 119 !>------------------------------V
!	                                              |      |    |    |    |
	          IF (JOUT == 0 .AND. IRFAS == 1) GOTO 5  !>-----------^
	          NF=NF-1
!	                                              |      |    |    |    |
	          write (*, 116)
	          WRITE(IOUT,116)
116	          FORMAT(/,'  Phase removed in DSUB ',/)
	          IF (JOUT > NF) ISTAT=0
!	                                              |      |    |    |    |
	          IF (JOUT > NF) GO TO 115  !>-------^
!	Reset yield arrays to smaller number of phases
	          DO 112 J=JOUT,NF
!	                                              |      |    |    |    |
	            SFAS(J)=SFAS(J+1)
	            indF(J) = indF(J+1)
	            DO 112 I=1,N
!	                                              |      |    |    |    |      
	              XVL(I,J)=XVL(I,J+1)
112	              YVAL(I+(J-1)*N)=YVAL(I+J*N)
!	                                              |      |    |    |    |
	        GO TO 115  !>-------------------------^      |    |    |    |
!	                                                     |    |    |    |
119	        CONTINUE  !<------------------------------------------------<
!	                                                     |    |     
	        IF (GNORM.LT.GLIM) GOTO 440  !>-----------------------------V
!	                                                     |    |    |    |
	          IF (L == LMAX) GOTO 230  !>---------------------V         |
!	                                                     |    |    |    |
!     Save result of direct substitution
!	                                                     |    |    |    |
130	            DO 140 I=1,NT  !<------------------------<    |    |    |
!	                                                          |    |    |
	              IF (KV > 0) WORK(I,KV)=GRAD(I)
	              Y(I)=GRAD(I)-YOLD(I)
140	              YOLD(I)=GRAD(I)
!	                                                          |    |    |
	  M1=M-1
	  DO 160 I=1,M1
	    DO 160 J=I,M1
!	                                                          |    |    |   
	      B(I,J)=0.
	      DO 150 K=1,NT
!	                                                          |    |    |
	        SS=(WORK(K,I+1)-WORK(K,I))*(WORK(K,J+1)-WORK(K,J))
150	        B(I,J)=B(I,J)+SS
!	                                                          |    |    |
160	      B(J,I)=B(I,J)
!	                                                          |    |    |
!	Convergence; danger of underflow in extrapolation
!	                                                          |    |    |
	  IF (B(M1,M1).LT.1.D-20) GOTO 440  !>------------------------------V
!	                                                          |    |    |
	  DO 170 J=1,M1
!	                                                          |    |    |
	    B(M1,J)=0.
170	    BAU(J)=0.
!	                                                          |    |    |
	  BAU(M1)=1.D0
	  B(M1,M1)=1.D0
!	                                                          |    |    |
!	Solve linear equations for extrapolation
!	                                                          |    |    |
!	  write(*,'(6(6(G16.7,X),/))')B
	  CALL LU(6,M1,IP,B)
	  CALL BACK(6,M1,IP,B,BAU)
!	                                                          |    |    |
	  SNY=0.
	  DO 180 J=1,M1
!	                                                          |    |    |
180	    SNY=SNY+BAU(J)
!	                                                          |    |    |
	  SLIM=.5
	  CPR=0.
	  SMAX=SLIM
!	                                                          |    |    |
!	Check that extrapolate is in correct direction
!
	DO 200 I=1,NT
!	                                                          |    |    |
	  YNEW(I)=YOLD(I)
	  YOLD(I)=0.
	  DO 190 J=1,M1
!	                                                          |    |    |
190	    YOLD(I)=YOLD(I)+BAU(J)*WORK(I,J+1)
!	                                                          |    |    |
	  YOLD(I)=YOLD(I)/SNY
	  DS=YOLD(I)-YNEW(I)
	  IF (DABS(DS) > SMAX) SMAX=DABS(DS)
200	  CPR=CPR+DS*Y(I)
!	                                                          |    |    |
!	  Output energy and coefficients of characteristic equation
	  write (*, 210) GMIN,(BAU(J),J=1,M1)
	  WRITE(IOUT,210) GMIN,(BAU(J),J=1,M1)
210	  FORMAT(5X,D15.5,3X,6F12.7)
!	                                                          |    |    |
	  IF (CPR > 0..AND.SMAX.LE.SLIM) GOTO 230  !>------------V    |    |
!	                                                          |    |    |
!	Limit extrapolation factor to improve stability
!	                                                          |    |    |
	  FAC=SLIM/SMAX
!	Change sign, if in wrong direction
	  IF (CPR.LT.0.) FAC=-FAC
	  DO 220 I=1,NT
!	                                                          |    |    |
220	    YOLD(I)=YNEW(I)+FAC*(YOLD(I)-YNEW(I))
!	                                                          |    |    |
230	  CONTINUE  !<--------------------------------------------<    |    |
!	                                                               |    |
240	CONTINUE
!	                                                               |    |
!	And this is where the Broyden Method used to be
!	                                                               |    |
!	                                                               |    |
!	Save best results seen so far
!440 	IF (NEWVAL.LE.ISTAT) GOTO 445
!	                                                               |    |
440	continue  !<--------------------------------------------------------<
	IF (NEWVAL.LE.ISTAT) GOTO 445  !>-------------------------V    |
!	                                                          |    |
	  write (*, 126)
! 	  iOut = 0
	  WRITE(IOUT,126)
!	                                                          |    |	  
	  NF=NF+1
	  ISTAT=1
	  DO 442 I=1,N
!	                                                          |    |	  
442	    XVL(I,NF)=0.
!	                                                          |    |
	  SFAS(NF)=0.
!	                                                          |    |	  
	  GOTO 5  !>---------------------------------------------------^
!	                                                          |
445	DO 450 I=1,N  !<------------------------------------------<

	  DO 450 J=1,NF
	  
	    VBEST(I+(J-1)*N)=YVAL(I+(J-1)*N)
450	    XVL(I,J)=YVAL(I+(J-1)*N)

	write (*, 460)GMIN,GNORM
	WRITE(IOUT,460)GMIN,GNORM
460	FORMAT(/,'  END DSUB  GMIN =',D25.15,' GNORM =',D15.3/)

      RETURN
      END
!	------------------------------------------------------------------
!
!
!     NTOPI EVALUATES FUGACITY COEFFICIENTS FROM PHASE COMPOSITIONS
!     AND ALSO TESTS FOR CONVERGENCE IN DSUB
!
	SUBROUTINE NTOPI(N,GTOT,FG,H,FX,VECT,YR)
	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4

	DIMENSION             FG(NCM), FP(NCM), H(NCM),FX(NCM,NCM),VECT(NCM*(NFM-1)),X(NCM), YR(NCM,NFM)
	COMMON /SOLID/        FUGAR(NCM), ZLOG(NCM), Zcomp_pure(NCM), ISOLC(NCM), ISOLPC
	COMMON /KEEP/         T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON/NTEST/         GNORM,GLIM, NFA,NDUM
	COMMON/NEWSOL/        NEWFS

      NEWFS=0
      GNORM=0.
      GTOT=-GNUL
      NREF=NF
!	- PHASE NREF IS REFERENCE PHASE FOR K-FACTORS
      IF (ISOLPC /= 0) NREF=NF-1
      NG=NREF-1
      SF=0.
!	--PROPERTIES OF REFERENCE PHASE
      DO 10 I=1,N
      XI=ZZ(I)*XVL(I,NREF)
      IF (SFAS(NREF) == 0.) XVL(I,NREF)=0.
      X(I)=XI
   10 SF=SF+XI
      IF (SF > 0.) SFAS(NREF)=SF
!	-OBJECTIVE FUNCTION EVALUATION
      CALL SRK (1,0,IC,T,P,ZVAL(NREF),X,FG,FP,FX)
!
	indF(nRef) = iC
!
      DO 20 I=1,N
      YI=FG(I)+DLOG(X(I)/ZZ(I)/SF)
      YR(I,NREF)=YI
      IF (SF.LE.0.) GOTO 20
      GTOT=GTOT+X(I)*YI
   20 CONTINUE
      IF (NG == 0) GOTO 50
      DO 40 J=1,NG
      SF=0.
      DO 30 I=1,N
      XI=XVL(I,J)*ZZ(I)
      IF (SFAS(J) == 0.) XVL(I,J)=0.
      X(I)=XI
   30 SF=SF+XI
      IF (SF > 0.) SFAS(J)=SF
	  CALL SRK (1,0,IC,T,P,ZVAL(J),X,H,FP,FX)
!
	indF(j) = iC
!
      DO 40 I=1,N
      HI=H(I)
      FAC=HI+DLOG(X(I)/SF/ZZ(I))
      IF (SFAS(J) == 0.) GOTO 41
      GTOT=GTOT+X(I)*FAC
      YR(I,J)=FAC
      IF (X(I).LT. 1.D-15) GOTO 41
      GNORM=GNORM+(ZZ(I)*(FAC-YR(I,NREF)))**2
!
!     LOG K-FACTORS STORED IN VECT
!
   41 CONTINUE
   40 VECT(I+(J-1)*N)=FG(I)-HI
   50 IF (ISOLPC == 0) GOTO 80
      SFAS(NF)=0.
      DO 60 I=1,N
      IPS=(NF-2)*N+I
      VECT(IPS)=0.
      IF (ISOLC(I) == 0) GOTO 60
      ZZI=ZZ(I)
      XI=XVL(I,NF)*ZZI
      SFAS(NF)=SFAS(NF)+XI
      HI=FUGAR(I)
      VECT(IPS)=FG(I)-HI
      FAC=HI-ZLOG(I)
      GTOT=GTOT+XI*FAC
      IF (XI.LT.1.D-15) GOTO 60
      IF (ISOLC(I) == 2) GNORM=GNORM+(ZZI*(FAC-YR(I,NREF)))**2
   60 CONTINUE
      RETURN
   80 IF (GNORM >  1.D-8) RETURN
      DO 90 I=1,N
      IF (ISOLC(I) == 0) GOTO 90
      VALUE=FUGAR(I)-ZLOG(I)
      JM=0
      DO 100 J=1,NF
      IF (YR(I,J).GE. VALUE) GOTO 100
      JM=1
  100 CONTINUE
      IF (JM > 0) GOTO 90
      ISOLC(I)=2
      NEWFS=1
   90 CONTINUE
	RETURN
	END
!------------------------------------------------------------------------
!
!	SRKPUR evaluates pure component properties and solid fugacity
!	coefficients for use in stability analysis
!
	SUBROUTINE SRKPUR(JPUR,OLDY,FX,N)
	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 30, NFM = 4

	DIMENSION             A(NCM),OLDY(NCM+2,NCM),FX(NCM,NCM),FG(NCM),FP(NCM),FW(NCM,NCM)
	
	COMMON /SOLPRP/       TMELT(NCM),HMELT(NCM), ACpMelt(NCM), BCpMelt(NCM), vMelt(NCM), Ptr(NCM)
	COMMON /SKIP/         NOTRY(NCM)
	COMMON /SOLID/        FUGAR(NCM), ZLOG(NCM), Zcomp_pure(NCM), ISOLC(NCM), ISOLPC
	COMMON /KEEP/         T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	
	IPT=0
	FMAX=0.
	DO 30 I=1,N
	
	  ISOLC(I)=1
	  
! 	  IF (T.GE. TMELT(I)) ISOLC(I)=0.
	  TRED = TMELT(I)/T
	  Y = -HMELT(I)*(TRED-1)
	  if (ACpMelt(i) > 0) then
	    Y = Y - ACpMelt(i)*dlog(TRED) + (ACpMelt(i) + BCpMelt(i)*TMelt(i)/2)*(TRED - 1) + BCpMelt(i)/2 * (T - TMELT(i))
! 	   
	  endif
	  if (vMelt(i) > 0) then
          Pfus = Y*T/vMelt(i) + Ptr(i)
	  endif
	  IF ((T >= TMELT(I) .AND. vMelt(i) == 0) .OR. (vMelt(i)*(P - Pfus) < 0)) ISOLC(I)=0.	  
	  IF (ISOLC(I) > 0) IPT=1
! 	  IF (ISOLC(I) == 0 .AND. NOTRY(I) == 0) GOTO 30
	  if (ISOLC(i) /= 0 .OR. NOTRY(i) /= 0) then
	    
	    DO 10 J=1,N
10	      A(J)=0.
!
	    A(I)=1.
!
!	    Pure component properties are here calculated from full N-compo-
!	    nent thermodynamics. it is, of course more efficient so substi-
!	    tute a separate routine (easy for cubic equations)
!         
	    IC = 0
	    CALL SRK (1,1,IC,T,P,ZZZ,A,FG,FP,FW)
	    Zcomp_pure(i) = ZZZ
!         
!	    log fugacity (coefficient) difference stored
!         
	    IF (P.LT. 1. .AND. IC /= 1) ISOLC(I)=0.
!         
!	    Evaluate properties of solid phase for below-melting comp.
!         
	    TRED = TMELT(I)/T
	    FUGAR(I) = FG(I)-HMELT(I)*(TRED-1)
	    if (vMelt(i) > 0) then
	      FUGAR(i) = FUGAR(i) - vMelt(i)*(P - Ptr(i))/T
	      Zcomp_pure(i) = Zcomp_pure(i) - vMelt(i)*P/T
	    endif
	    if (ACpMelt(i) > 0) FUGAR(i) = FUGAR(i) - ACpMelt(i)*dlog(TRED) + (ACpMelt(i) + BCpMelt(i)*TMelt(i)/2)*(TRED - 1) + BCpMelt(i)/2 * (T - TMELT(i))
	    DO 20 J=1,N
20	      FX(I,J)=FG(J)

	  endif

30	CONTINUE
	  IF (IPT > 0) GOTO 36

	    GOTO 48

36	  write (*, 35)

35	  FORMAT(/,'  Solid formation investigation for components ',//,'       T MELT   LN (FUG/P)    ',/)
	  DO 40 I=1,N
	    IF (ISOLC(I) == 0) GOTO 40
	      write (*, 45) I,TMELT(I),FUGAR(I)
45	      FORMAT(1X,I3,2X,F7.2,2X,F9.5)
40      CONTINUE
48	  SMAX=0.D0

!	  NOTE revision from earlier version
	  DO 50 I=1,N	
	    IF (NOTRY(I) == 0) GOTO 50
	  
	      DO 55 J=1,N	    
55	        OLDY(I,J)=FR(J)-FX(I,J)

	      IF (OLDY(I,I).LT. SMAX) GOTO 50
      
	        SMAX=OLDY(I,I)
	        JPUR=I
	      
50	CONTINUE


	IF (JPUR == 0) RETURN	! En el COMIENZO, se llama a esta subrutina
	                        ! solamente por la capacidad interna de tra-
	                        ! tar con fases sólidas. Sin esto, entra y 
	                        ! sale sin hacer nada.
	
	
	SMAX=SMAX+3
	DO 60 I=1,N
	  IF (NOTRY(I) == 0) GOTO 60
	  
	    IF (I /=  JPUR) OLDY(I,JPUR)=OLDY(I,JPUR)-SMAX
	    
60	CONTINUE

      RETURN
      END
!-----------------------------------------------------------------------
!
!
!	Subroutine PURE initiates the stability analysis, scanning trial
!	phases based on pure component properties for instability
!
	SUBROUTINE PURE (Y,N,NH,FX)
	IMPLICIT REAL*8(A-H,O-Z)
	integer, parameter :: NCM = 30, NFM = 4

	DIMENSION             LIST(NCM+2),RTOLD(NCM+2),OLDY(NCM+2,NCM)
	DIMENSION             Y(NCM),YTEST(NCM),FA(NCM),FP(NCM),DS(NCM),DSG(NCM),FX(NCM,NCM)
	COMMON /KEEP/         T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /DELST/        YV1
	COMMON /SKIP/         NOTRY(NCM)
	COMMON /SOLID/        FUGAR(NCM), ZLOG(NCM), Zcomp_pure(NCM), ISOLC(NCM), ISOLPC
	COMMON /MOLFRC/       ZFRAC(NCM,NFM)
	common /OUTUNIT/      iOut

! 	iOut = 0

	JPUR=0
	JMIN=0
	KMIN=1
	NT=N+1
!
!	List identifies elements, where search is skipped
!
	DO 10 I=1,NT
	  IF (NF == 1) LIST(I)=1
	  IF (LIST(I) /= 0) LIST(I)=1
10	CONTINUE
!	--------------------
	IF (NF == 1) GOTO 40  !--------------------------V
!	                                                 |      
	  NT=N+2
	  ZMAX=0.
	  JMAX=0
	  DO 20 J=1,NF
!	                                                 |	  
	    IF (ZVAL(J).LT. ZMAX) GOTO 20
!	                                                 |	    
	      ZMAX=ZVAL(J)
	      JMAX=J
20	  CONTINUE
!	                                                 |
!     Select additional primer phase in multiphase analysis
!	                                                 |
	  DO 30 I=1,N
!	                                                 |	  
	    YCM=ZFRAC(I,JMAX)
	    YPN=(3*ZZ(I) + 2*YCM)/5
30	    OLDY(NT-1,I)=DLOG(YPN)
!	                                                 |
	  LIST(NT-1)=1
!	                                                 |	  
40	  LIST(NT)=1  !<---------------------------------<
!	-------------------

	IF (NF /= 1) GOTO 45  !----V
!                                |      
	  IF (ISOLPC == 0) CALL SRKPUR(JPUR,OLDY,FX,N)
!                                |	  
!	Ideal gas trial phase
45	NA=3  !<-------------------<
	DO 50 I=1,N
	  IF (NOTRY(I) == 0) LIST(I)=0
50	  OLDY(NT,I)=FR(I)                ! oldY(NT,i) = fr(i) = ln( z(i)*phiZ(i) )
!	-------------------------V
	IF (JPUR.LE.0) GOTO 55
	
	  DO 56 I=1,N
56	    Y(I)=1.D-10

	  Y(JPUR)=1.
	  YV1=-1
	  RETURN
!	-------------------------^
55	DO 140 KK=1,NA
      
	  KMIN=3
	  RTMIN=1.D20
	  DO 130 K=1,NT
!
!	Search phase or not ?

	    IF (LIST(K) == 0 .OR. LIST(K).GE.3) GOTO 130 !>--V
!	                                                     |   
!     N O T E:  Trial phases 1 TO N are N O T restored
!     ********  It may (or may not) be better to restore
!	                                                     |
!	Specify type of primer phase
!	                                                     |
	      MV=1
	      IF (K == NT) MV=-1
	      SUM=0.D0
	      DO 60 I=1,N
!	                                                     |	  
	        Z=1.D-25
	        IF (OLDY(K,I) > -80.) Z=DEXP(OLDY(K,I))	! This suposes that phiNew = 1,  so isofug.
	        YTEST(I)=Z	                              ! is simply yTest(i)*1*P = z(i)*phiZ(i)*P
60	        SUM=SUM+Z
!	                                                     |
	      DLSUM=DLOG(SUM)	!ln de la suma de los nros de moles de la fase de testeo?
!	                                                     |	  
	      CALL SRK (1,MV,IC,T,P,Z,YTEST,FA,FP,FX)
	      DMIN=N
!	  Guarda en FA los ln(phi(yTest))  
!	                                                     |
!	Find closest phase for comparation
!	                                                     |
	      NFF=NF
	      IF (ISOLPC == 1) NFF=NF-1
!	                                                     |	      
	      IF (NFF == 0) NFF=1
	      DO 90 J=1,NFF
!	                                                     |     
	        DIST=0.D0
	        DO 70 I=1,N
!	                                                     |      
	          Z=ZFRAC(I,J)
!	                                                     |     
	          IF (J == 1) YTEST(I)=YTEST(I)/SUM	! Normaliza las fracciones de la
	                                                ! fase de testeo
	          DS(I)=YTEST(I)-Z	                  
70	          DIST=DIST+DS(I)**2	                  ! Desviación cuadrada
!	                                                     |
!	      --------------------------------V
	        IF (DIST > DMIN) GOTO 90
!	                                                     |	 
	          DMIN=DIST
	          MINVL=J
	          DO 80 I=1,N
!	                                                     |	    
80	            DSG(I)=DS(I)
!	      --------------------------------^
90	      CONTINUE
!	                                                     |
	      VAL=0.D0
	      GRAD=0.D0
	      DO 100 I=1,N
!	                                                     |      
	        DNEW=FR(I)-FA(I)	!fr(i) = ln(z(i)*phi(z(i)))
	                        	!fa(i) = ln(phi(yTest(i)))
	                        	!this is eq. (16) in page 235. of Michelsen&Mollerup book.
!	                                                     |
	        GD=OLDY(K,I)-DNEW-DLSUM
	        OLDY(K,I)=DNEW
	        VAL=VAL+YTEST(I)*GD
100	        GRAD=GRAD+GD*DSG(I)
!	                                                     |
!	Evaluate result; check for negative value.
!	                                                     |
	      RT=1+VAL	!RT seems to be the tangent plane distance
	      IF (GRAD > 0.) RT=2*VAL/GRAD
!	                                                     |     
	      IF (VAL.LT.0.) RT=VAL
!	                                                     |      
	      IF (DMIN.LT.2.D-3.AND.RT > .8) LIST(K)=0
!	                                                     |      
	      IF (KK == 1) GOTO 110  !------------------V    |
!	                                                |    |      
	        IF (RT.LT.0.) GOTO 110 !----------------V    |
!	                                                |    |       
	          IF (LIST(K) == 0) GOTO 130  !--------------V
!	                                                |    |
	            IF (RT.LT.RTOLD(K)) GOTO 110  !-----V    |
!	                                                |    |      
	            LIST(K)=LIST(K)+1
!	                                                |    |	        
	            IF (RT > 1.) LIST(K)=LIST(K)+1
!	                                                |    |	        
110	      RTOLD(K)=RT  !<---------------------------<
!	                                                     |        
	      IF (LIST(K).GE.3) GOTO 130  !------------------V
	        IF (RT > RTMIN) GOTO 130  !-----------------V
!	                                                     |              
	        KMIN=LIST(K)
	        RTMIN=RT
!	                                                     |
	        JMIN=MINVL
	        DO 120 I=1,N
120	          Y(I)=YTEST(I)
!	                                                     |              
130	  CONTINUE  !<---------------------------------------<
!
!	Search conclusive for negative RT
!
	  IF (RTMIN.LT.0.) GOTO 150 !>---V	!(If there's already a negative  tangent plane distance)
!	                                 |
140	CONTINUE
!	                                 |
150	YV1=0. !<------------------------<

	IF (RTMIN.LT.0.) YV1=RTMIN
	
	IF (YV1.LT.-1.D-7) GOTO 240    !>-----------------------------v
!	                                                              |
!	Search inconclusive. Try also reflection with respect to the
!	closest "trivial" solution
!
!	(This section arises to see if there's a numerical error)
!	  Francisco, 27/02/14
!	                                                              |
	  DO 160 I=1,N
!	                                                              |	
	    DS(I)=ZFRAC(I,JMIN)
	    Z=DS(I)/Y(I)
	    IF (Z > 2.) Z=2.
	    IF (Z.LT..5) Z=.5
160	    YTEST(I)=DS(I)*Z
!	                                                              |
	  DO 180 K=1,NA
!	                                                              |	
	    CALL SRK (1,1,IC,T,P,Z,YTEST,FA,FP,FX)
	    IF (K == NA) GOTO 180
!	                                                              |	    
	      DO 170 I=1,N
	  
	        Z=FR(I)-FA(I)
	        IF (Z.LT.-80.) Z=-80.
170	        YTEST(I)=DEXP(Z)

180	  CONTINUE
!	                                                              | 
	  SUM=0.D0
	  DO 190 I=1,N
190	    SUM=SUM+YTEST(I)
!	                                                              |
	  VAL=0.D0
	  GRAD=0.D0
	  DO 200 I=1,N
!	                                                              |
	    YTEST(I)=YTEST(I)/SUM
	    GD=DLOG(YTEST(I))+FA(I)-FR(I)
	    VAL=VAL+YTEST(I)*GD
200	    GRAD=GRAD+GD*(YTEST(I)-DS(I))
!	                                                              |
	  IF (VAL.LT.0.) GOTO 210
!	                                                              |      
	    RT=1+VAL
	    IF (GRAD > 0.) RT=2*VAL/GRAD
	    IF (KMIN == 0.OR.KMIN == 3) GOTO 210
!	                                                              |	  
	      RTMIN=(1+KMIN)*RTMIN/2
	      IF (RT.LT.RTMIN) GOTO 210
!	                                                              |
!	        EXIT!!!   
	        GOTO 240	
!	                                                              |  
210	      IF (VAL.LT.0.) YV1=VAL
	  DO 230 I=1,N
230	    Y(I)=YTEST(I)
! 	  iOut = 0
!	                                                              |	  
240	WRITE(IOUT,250)(Y(I),I=1,N)  !<-------------------------------<

250	FORMAT(/,' Composition for stability test ',/,(5X,12F9.6))
	RETURN
	ENDSUBROUTINE PURE
!-----------------------------------------------------------------------
!
!
!       OUT PRINTS RESULT OF (INTERMEDIATE) PHASE SPLIT CALCULATION
!       NOTE THAT SOLID PHASES ARE COMBINED IN THE LAST PHASE, IF
!       PRESENT, AND THAT THEIR "MOLE FRACTIONS" INDICATE THEIR RELA-
!       TIVE AMOUNTS
!
	subroutine OUT (Yval, N, Fx)
!
	IMPLICIT REAL*8(A-H,O-Z)
	
	integer, parameter :: NCM = 30, NFM = 4
	DIMENSION             YVAL(NCM*NFM), DLX(NCM), FP(NCM), FX(NCM,NCM)
	COMMON /SOLID/        FUGAR(NCM), ZLOG(NCM), Zcomp_pure(NCM),ISOLD(NCM), ISTAT, ID_S
	COMMON /KEEP/         T,P,GNUL,Z(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),ZVAL(NFM),indF(NFM),NF,ND
	COMMON /MOLFRC/       ZFRAC(NCM,NFM)
	common /OUTUNIT/      iOut
 
! 	NT = NF*N
! 	NB = NT - N
	if ( ISTAT /= 0 ) then
 
	  write (*, 11) NF
	  write (IOUt, 11) NF
11      FORMAT(/,'  NOTE:   PHASE  ',I2,'  IS SOLID ',/)
	  zFrac(:N,NF) = XVL(:N,NF)*Z(:N)/sFas(NF)
	  zVal(NF) = dot_product(zFrac(:N,NF), Zcomp_pure(:N)) !average compressibility factor of solid phase, 
	  
	endif
!
!	Change output formats for more than 4 phases
!
	write (IOUT, 20) NF
20	format (/,'  FLASH RESULT ', I3,' PHASES * * * * ')
	write (*,30) (j, sFas(j), j = 1, NF)
	write (IOUt,30) (j, sFas(j), j = 1, NF)
30    format(/,'  TOTAL AMOUNT  ', 4(2X, I2, F12.8))

	write (*, 40)(J,ZVAL(J),J=1,NF)
	write (IOUT, 40) (J,ZVAL(J),J=1,NF)
40  	format('  COMPRESS      ',4(2X, I2, F12.8))
 
	write (*, 50) (j, j = 1, NF)
	write (IOUt, 50) (j, j = 1, NF)
50    format(/,'  COMPOSITION   ',4(5X, I2, 9X))

	sum = 0
	NFF = NF
	if ( ISTAT == 1 ) NFF = NF - 1           !Number of Fluid Phases
	  
	if ( NFF /= 0 ) then
 
	  do i = 1, N
 
	    dlx(i) = xVL(i,NFF)*z(i)/sFAS(NFF)   !composition of last fluid phase
	    sum = sum + dlx(i)	                 !total amount of last fluid phase
!
	  enddo 
	  sum = dlog(sum)
	  call SRK (1, 0, ic, T, P, ZZ, dlx, FR, fp, Fx)
	  do i = 1, N
 
	    dlx(i) = dlog(dlx(i))
	    FR(i) = dlx(i) + FR(i) - sum         !ln(n(i)) + ln(phi(i)) - ln(Ntot) = ln(f(i)/P) = fugacity of a reference phase
	    
	  enddo
	  
	else
 
	  FR(:N) = FUGAR(:N)                     !Reference fugacity is that of the solid phases 
	  
	endif
	
	do i = 1, N
 
	  do j = 1, NFF

	    zFrac(i,j) = XVL(i,j)*Z(i)/sFas(j)
	    
	  enddo
 
	  write (IOUt,'("  X(",I2,")         ", 5(1X, G15.7))') i, (zFrac(i,j), j = 1, NF)

	enddo
	return
	endsubroutine OUT
!**************************************************************
