C
C       Multiphase P/T Flash program, version 3, may 1983
C       Copyright: Michael L. Michelsen, Instituttet for Kemiteknik,
C                  Bygning 229, DTH, Lundtoftevej 100, DK 2800 Lyngby.
C                  DENMARK;   TEL 02-88 32 88
C
C                  THERMODYNAMIC SUBROUTINES BY STEEN SKJOLD-JOERGENSEN
C
C     M a i n     p r o g r a m  for multiphase flash
C     ************************************************
C
C     Common areas in program complex do not interfere with common
C     areas in thermodynamic part except in /STYR/ In the present
C     interface routine SRK
C
C		-Esto ya no es así, a partir de la versión 1.2.0 de la GC3: se impri-
C	       men tanto la versión del paquete termodinámico (leído dentro del pa-
C	       quete) como la versión de esta subrutina. Implica un nuevo COMMON.
C			Francisco, 20/05/2011, durante v1.2.0
C
C		-Además, ahora hay un common en "iEntalp" para ver si calcula o no entalpías
C------------------------------------------------------------------
C
C	Aprovechando las modificaciones introducidas en la versión 1.3.3 del paque-
C	te termodinámico GCA, se modifica esta subrutina para poder imprimir valo-
C	res de Z (coef de compr) afectados por una corrección de Peneloux.
C
C	El usuario debe tener en cuenta que esta corrección no es "mágica", se pasa
C	de tener una predicción pésima, a una simplemente mala.
C
C	La corrección de Peneloux ingresada es siempre negativa (al menos en la ver-
C	sión 1.3.3), ya que el V predicho por la GC es demasiado grande. La línea de
C	ingreso es la siguiente a las propiedades críticas.
C
C	Ej: líneas correspondientes al compuesto "H2O":
C
C	H2O                                                               !nombre
C	647.13    217.67    0.345     373.15     1.0000     2.5763        !Tc  Pc  w  Ts  Ps   dc
C	18	-10	0	0	0	0	0	0				!PM  Peneloux   coeficientes de H(gas ideal) de Passut-Danner
C
C
C	v1.2.90 (1.3.0 beta0): comienzo de soporte para GCA-v2.0
C	     91                cambiada dimensión de versSUB
C
C	v1.2.92 	           Limpieza de algunas impresiones
C	                       Cálculo de entalpías
C
C	v1.2.93: limpieza de código.
C	         algunas identaciones
C	         agregado parámetros NCM
C
C	v1.2.96: PARAGC llamado con nTemp = 1, pero no así GCEOS.
!
!	v1.2.100: las versiones 98-99 están falladas. Esta comienza siendo un clon de la 96, incluyendo modificaciiones
!	          hechas en 97 y 98.
!
!	v-1.2.100.0-paraest: 	Modificación para usar el programa como subrutina en el ParaEst. La versión
!	anterior estaba hecha con las "97" y "98", que estaban falladas en algún punto.
C
	subroutine GC3Flash(NComp,NGrup,NSit,TEMP,PRES,feed,xFrac,
     &                    Zcomp,indFase,IERRTP)
C
      IMPLICIT REAL*8(A-H,O-Z)
	parameter(NCM=15,NFM=4)
      DIMENSION Y(NCM)
      DIMENSION YVAL(60)
      COMMON/DERIV/FX(NCM,NCM)
      COMMON/KEEP/NF,ND,T,P,GNUL,Z(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),
     *            ZVAL(NFM),indF(NFM)
      COMMON/SOLID/FUGAR(NCM),ZLOG(NCM),ISOLD(NCM),ISTAT
      COMMON/STORE/GRAD(45),XMAT(45,45),WORK(45,5)
      COMMON/DELST/YV1
      COMMON/AKTIV/IVAL(NCM),IVAR(NCM),NACT
      COMMON/LIMIT/EPSV
      COMMON/NTEST/NFA,NDUM,GNORM,GLIM
      COMMON/MOLFRC/ZFRAC(NCM,NFM)
      COMMON/OLDZ/ZOLD1,ZOLD2
      COMMON/OLDT/TOLD
C
C-----SUbroutine arguments for Murray routine
      EXTERNAL FSTAB,PHASN
C
C-----Nuevas declaraciones incluidas en la versión 1.2.0------------------------
	dimension Hgi(NCM)
	character*10 aname
	character*20 versGCA
	character*20::version='1.2.100.0-paraest'
	common/versSUB/versGCA
	common/NAME/AName(15)
	common/GCPROM/pmm(ncm),Pen(NCM),hha(ncm),hhb(ncm),hhc(ncm),hhd(ncm),
     &              HHE(NCM),HHF(NCM)
	common/ENTALP/iEntalp
!
!	Declaraciones agregadas al incorporar a ParaEst:
	dimension beta(NFM),Zcomp(NFM),feed(NCM),xFrac(NCM,NFM)
	dimension indFase(NFM),iVap(1)
      COMMON/STYR/N,NG,NST,NDER,NTEMP
C-------------------------------------------------Francisco, 20/05/2011---------
C
	npunto = 0    !selva
      TOLD=0.D0
C	GLim = 1e-8
	GLIM=1.D-16
   10 FORMAT(8F10.5)
500	format(/,72("-"))
C------INPUT NUMBER OF COMPONENTS AND TOTAL NUMBER OF GROUPS IN MIXTURE
C
	IOUT=0
!	OPEN(UNIT=5,FILE='GCTHREE.DAT',STATUS='OLD')
	OPEN(UNIT=IOUT,FILE='GCTHREE.SAL',STATUS='unknown',FORM='FORMATTED
     *')
!	open(unit=1,file='datter.dat',status='unknown',form='formatted')
!	open(unit=9,file='Compres.dat',form='formatted')
!------------------------------------------------------------------------------------
!	write(1,500)
!	write(1,*)'Programa de cálculo flash multifásico GC3, versión ',
!     &version,'usando el paquete termodinámico GCA-v',versGCA
!	write(1,*)
!	write(1,500)
!	write(iOut,500)
!	write(iOut,*)'Programa de cálculo flash multifásico GC3, versión ',
!     &version,'usando el paquete termodinámico GCA-v',versGCA
!	write(iOut,*)
!	write(iOut,500)
!------------------------------------------------------------------------------------
!	READ(5,*)NCOMP,NG,NGA,iOut
20	FORMAT(16I3)
!
!	ndat sets up tables and common areas for termodynamic subroutine
!      CALL INDAT(NCOMP,NG,NGA,NST,IOUT)
      N=NCOMP
      NG = NGrup
      NST = NSit
	T = Temp
	P = Pres
!
!	write(1,*)'Fases     Z          P(atm)    T(K)   Frac   H(J/mol)',
!     &	 '     Composición'	       !Fran, durante v1.2.0
!	write(1,15)AName(:N)
!15	format(' ----- ---------   ----------  ------ ----- ------------'
!     *,2X,16(A10,X))
!
!	Read mixture composition, normalization not needed
!260	READ(5,*,ERR=250,END=250)(Z(I),I=1,N),T,P
!260	READ(5,*)(Z(I),I=1,N),T,P
	z = feed
	suma = sum(z(:N))
      IF(suma.LT.1.D-10) GOTO 250
C-----NORMALIZE:
	z(:N) = z(:N)/suma
C
C     Input temperature and pressure
C     Units must correspond to (absolute) units in indat
      iSTAT=0
      ND=0
!     IF (T.EQ.0.) GOTO 260
	if(T.EQ.0)return
C
      CALL SRK(1,0,IC,T,P,ZX,Z,Y,Y,FX)
	iTypOrig = iC
C
C     EPSV is tolerance for trace components, it is modified according
C     to the result of the stability test
      EPSV=1.D-4
      WRITE(6,60) T,P,ZX
      WRITE(IOUT,60) T,P,ZX
	WRITE(9,60) T,P,ZX				!Fran 07/08/09 -> JA!!
60	FORMAT(1H1,///,' PT-Flash, T = ',F8.2,' K, P = ',F8.2,' atm, Z = '
     *       ,F10.8)
      GNUL=0.D0
      SFAS(1)=1.D0
	DO 70 I=1,N

!	  K(i,1):
	  xVL(I,1)=1.D0
C	  Zero mole fraction not permitted (log-terms)
	  z(I) = z(I) + 1.D-20
	  zFrac(I,1) = z(I)
	  zLog(I) = dlog(z(I))
	  fr(I) = zLog(I) + Y(I)
	  Gnul = Gnul+z(I)*y(I)

70	enddo
C
!	WRITE(6,80)(I,Z(I),FR(I),I=1,N)
	WRITE(IOUT,80)(I,Z(I),FR(I),I=1,N)
80	FORMAT(/,8X,'Composition  ln z(i)*phi(i) ',/,36("-"),/,
     *       (I3,2X,F12.6,F14.6))
	NF=1
C
C-----Initiate stability analysis
90	CALL PURE(Y,N,NONHY,FX)
C
C     Sufficiently unstable, go to phase split
	IF (YV1.GE.-.01) GOTO 110
C
C       Otherwise locate minimum by murray method
	  WRITE(6,190)
	  WRITE(IOUT,190)
	  yVal(:N) = y(:N)/z(:N)
100	  continue
C
	  GOTO 210

110	DO 120 I=1,N

120	  YVAL(I)=DLOG(Y(I))
C
      xLAM=1.D0
      WRITE(6,130)
      WRITE(IOUT,130)
130	FORMAT(' ')
	DO 140 I=1,N

	  IVAR(I)=I
C  140   IVAL(I)=1
	  IVAL(I)=1

140	enddo
      NACT=N
C
C     Murray search for stability minimum
      CALL MURML(45,N,1,10,xLAM,1.D-12,FUN,yVal,GRAD,XMAT,WORK,FSTAB)
      IF (DABS(FUN) .LT. .01) EPSV=1.D-3*DSQRT(DABS(FUN))
      WRITE(6,150) NACT
      WRITE(IOUT,150) NACT
150	FORMAT(/,'  Variables active ',I4/)
C
C     Check result of search
160	IF (FUN.LT.-1D-7) GO TO 180		!Fran 19/07/09
C
C     Try onedimensional search if pure is negative
C     Tolerance set in accordance with convergence criterion.
C     Use of a smaller value requires a smaller Glim
      iTry=0
	DO 410 I=1,N

C	  IF (ISOLD(I).EQ.0) GOTO 410
	  if(iSold(i).NE.0)then

	    dif = fugaR(i) - fr(i)
	    if(dif.LE.-1D-7)then

	      iTry = 1
	      iSold(i) = 2

	    endif

	  endif

410	enddo

C	    DIF=FUGAR(I)-FR(I)
C	    IF (DIF .GT. -1.D-7) GOTO 410

C	      ITRY=1
C	      ISOLD(I)=2

C  410 CONTINUE

      IF (ITRY.EQ.0) GOTO 420

	  ISTAT=ISTAT+ITRY
	  WRITE(6,403)
	  WRITE(IOUT,403)
  403   FORMAT(/,'   Solid phase formation  ')
	  GOTO 218

  420 CONTINUE
C
C-----One dimensional search for close minima as last resort
      CALL ONEDIR(N,YV1,YVAL,FX)
C
      IF (YV1 .LT. -1.D-7) GOTO 180
      WRITE(6,170)
      WRITE(IOUT,170)
  170 FORMAT(/,'  Solution is now stable * * * ')
      nPunto = nPunto + 1
C
C	Entalpía del gas ideal, si se ingresaron coef de Passut&Danner de
C	al menos 1 componente
!	if(iEntalp.GT.0)call NTLPY_GI(N,T,Hgi)
	NFases = NF
	if(NF.GT.1)then
!	  do j=1,NF

C---------Corrección de Peneloux:
!	    CorrPen = dot_product(Pen(:N),zFrac(:N,j))
!	    Z_corr = zval(j)+CorrPen*P/82.05D0/T
!	    HR = 0D0
!	    H = 0D0
!	    if(iEntalp.GT.0)then

!	      Cálculo de entalpías:
!	      La variable indF guarda el resultado de los llamados a la sub-
!	      rutina termodinámica, evaluando la fase (L o V) con menor G
!	      iGz = -1
!	      if(indF(j).EQ.1)iGz = 2
!	      call NTLPY(N,NG,NST,T,P,zfrac(:NCM,j),HR,zVal(j),iGz,
!     *                 indF(j),Ie)
!
!	      H = HR + dot_product(Hgi(:N),zFrac(:N,j))

!	    endif
!	    write(1,171)NF,Z_corr,P,T,sFas(j),H*.101325,
!     *                (zFrac(i,j),i=1,N),nPunto

!	  end do
	  xFrac = zFrac
	  Zcomp = zVal
	  indFase = indF

	else

!	  Impresión a partir del cálculo original.
!	  CorrPen = dot_product(Pen(:N),z(:N))
!	  Z_corr = Zx + CorrPen*P/82.05D0/T
!	  HR = 0D0
!	  H = 0D0
!	  if(iEntalp.GT.0)then

!	    Cálculo de entalpías:
!	    La variable indF guarda el resultado de los llamados a la sub-
!	    rutina termodinámica, evaluando la fase (L o V) con menor G
!	    iGz = -1
!	    if(iTypOrig.EQ.1)iGz = 2

!	    call NTLPY(N,NG,NST,T,P,z(:NCM),HR,Zx,iGz,
!     *               iTypOrig,Ie)

!	    H = HR + dot_product(Hgi(:N),z(:N))

!	  endif
!	  write(1,171)NF,Z_corr,P,T,1.0,H*.101325,(z(i),i=1,N),nPunto

!	  Si el resultado es 1 fase en vez de 2, se vuelve con error:
	  iErrTP = 1
	  return

	end if
!	write(1,*)
171	format(2XI2,1X,F11.7,1X,F12.6,1X,F7.2,1X,F5.3,X,F12.3,X,
     *       <N>(1X,F10.8),1X,I5)
C
C
C------------------------------------------------------------------------------------
C
C     No more phases could be found
C
c	GOTO 50
	goto 260
180	WRITE(IOUT,190)
190	FORMAT(/,' System is still unstable, phase split performed '/)
	yVal(:N)=DEXP(yVal(:N))/z(:N)  !yVal(i) = phi(i)
200	continue
210	NF=NF+1
C
C     Add new phase and proceed with phase split
C	IF (NF.GT.4) GOTO 50		!modifique cuando cambie el modo de lectura de condiciones alim. flash
	IF (NF.GT.4) GOTO 260
C
C     Four phases maximum for present version
	DO 215 J=2,NF

	  JN=NF+2-J
	  SFAS(JN)=SFAS(JN-1)

	  DO 215 I=1,N

215 	    XVL(I,JN)=XVL(I,JN-1)

	DO 216 I=1,N

	  XVL(I,1)=1.D-16*yVal(I)
	  DO 216 J=1,NF

216	    yVal(I+(J-1)*N)=xVL(I,J)

	SFAS(1)=1.D-16
	GOTO 419
218	IF (ISTAT.EQ.1) GOTO 414

	  ISTAT=1
	  GOTO 419

C	Usa a NF luego, como fase de referencia
414	NF=NF+1
	SFAS(NF)=0.
	ZVAL(NF)=0.
	DO 416 I=1,N

416	  XVL(I,NF)=0.

419	CONTINUE
	XLAM=.2
	IF (NF.EQ.2) XLAM=.5
	M=(NF-1)*N
	DO 220 I=1,M

	  IF (I.LE.N) IVAL(I)=1
220	  IVAR(I)=I

	NACT=M
C
C     Phase split by direct substitution/Gdem
	CALL DSUB(N,FX,YVAL)
C-----Check for reduction of number of phases
	NF=NFA
	M=(NF-1)*N
	IF (Gnorm.LT. Glim) GO TO 230
C
C     Phase split by murray method
      CALL MURML(45,M,1,10,XLAM,GLIM,FUN,YVAL,GRAD,XMAT,WORK,PHASN)
  230 EPSV=1.D-4
      WRITE(6,150) NACT
      WRITE(IOUT,150) NACT
  240 CALL OUT(YVAL,N,FX)
      GO TO 90
  250 CLOSE(UNIT=5)
     	close(unit=1)   !selva
!
260	return
C
!	CLOSE(UNIT=6)
!	STOP
      END
C
C     END OF PROGRAM


C***********************************************************************
C
C     INDAT SETS UP DATA TABLES FOR ACTUAL MIXTURE
C
C***********************************************************************
C
      SUBROUTINE INDAT(NC,NGG,NGAa,NSitio,IOU)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SOLPRP/TMELT(15),HMELT(15)
      COMMON/CRIT/TCR(15),PCR(15)
      COMMON/STYR/N,NG,NST,NDER,NTEMP
      COMMON/SKIP/NOTRY(15)
	nga=ngaa
      NG=NGG
      N=NC
C
      DO 1 I=1,NC
      TMELT(I)=1.0D0
      HMELT(I)=0.D0
    1 NOTRY(I)=0
	iOut = iOu
C	write(*,*)'iOut =',iOut
      CALL PARGR(NC,NG,NGA,NST,IOUt)
	NSitio = NST
C	write(*,*)'iOut =',iOut
      CALL PARMOL(NC,NG,NST,IOUt)
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C     THIS IS THE END OF THERMODYNAMIC SUBROUTINES
C
C     START GENERAL EQUATION SOLVING ROUTINES AND MINIMIZER
C
C     LU PERFORMS TRIANGULAR DECOMPOSITION OF A FOR SOL'N OF LINEAR EQ.
C
      SUBROUTINE LU(ND,N,IPIV,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IPIV(ND),A(ND,ND)
      IPIV(N)=N
      N1=N-1
      IF (N1.EQ.0) GO TO 50
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
      IF (IPIV(I).EQ.I) GO TO 30
      K=IPIV(I)
      DO 20 J=I,N
      X=A(I,J)
      A(I,J)=A(K,J)
   20 A(K,J)=X
   30 DO 40 J=I1,N
      X=-A(J,I)/A(I,I)
      A(J,I)=X
      DO 40 K=I1,N
   40 A(J,K)=A(J,K)+X*A(I,K)
   50 RETURN
      END
C-----------------------------------------------------------------------
C
C
C     BACK SOLVES LINEAR EQUATIONS A*X=V BY BACK SUBSTITUTION
C     RESULT IN V
C
      SUBROUTINE BACK(ND,N,IPIV,A,V)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IPIV(ND),A(ND,ND),V(ND)
      N1=N-1
      IF (N1.EQ.0) GO TO 30
      DO 20 I=1,N1
      I1=I+1
      K=IPIV(I)
      IF (K.EQ.I) GO TO 10
      X=V(I)
      V(I)=V(K)
      V(K)=X
   10 DO 20 J=I1,N
   20 V(J)=V(J)+A(J,I)*V(I)
   30 V(N)=V(N)/A(N,N)
      IF (N1.EQ.0) GO TO 60
      DO 50 II=2,N
      I=N+1-II
      I1=I+1
      DO 40 J=I1,N
   40 V(I)=V(I)-A(I,J)*V(J)
   50 V(I)=V(I)/A(I,I)
   60 RETURN
      END
C-----------------------------------------------------------------------
C
C
C     MURML PERFORMS THE MINIMIZATION USING THE MURRAY MODIFIED CHOLESKY
C     DECOMPOSITION OF THE HESSIAN MATRIX
C
      SUBROUTINE MURML(ND,N,IPR,NMAX,XLAM,GLIM,F,X,GD,G,W,FUNC)
      IMPLICIT REAL*8(A-H,O-Z)

	parameter(NCM=15,NFM=4)

      DIMENSION X(NCM*NFM),GD(ND),G(ND,ND),W(ND,5)
      COMMON/AKTIV/IVAL(NCM),IVAR(45),NACT
      NTM=0
      DEL=1.D-7
   10 NTM=NTM+1
      GPS=0.
C
C     FUNC IS THE SUBROUTINE CALL TO THE FUNCTION TO BE MINIMIZED
C
      CALL FUNC(N,2,F,GD,G,X)
      DO 20 I=1,N
      W(I,4)=X(I)
      W(I,5)=GD(I)
      IF (I.GT.NACT) GOTO 20
      IV=IVAR(I)
      GPS=GPS+GD(IV)**2
   20 CONTINUE
C
C     GLIM IS THE TOLERANCE FOR THE SQUARED GRADIENT NORM
C     IF RESULT IS LESS THAN TOLERANCE, EXIT
C
      IF (GPS.LT.GLIM) GOTO 170
      BETA=0.
C
C     SET UP HESSIAN MATRIX FOR DECOMPOSITION
C
      DO 30 I=1,NACT
      IV=IVAR(I)
      DO 30 J=1,I
      JV=IVAR(J)
      G(I,J)=G(IV,JV)
      S=DABS(G(I,J))
      IF (I.EQ.J) S=S*NACT
      IF (S.GT.BETA) BETA=S
   30 W(I,1)=0.
      BETA=DSQRT(BETA/NACT)
C
C     SUBROUTINE SPLIT PERFORMS DECOMPOSITION
C
      CALL SPLIT(ND,N,IDI,BETA,DEL,G,W)
      XLM=0.
      NTS=0
   40 NTS=NTS+1
C
C     LINE SOLVES THE LINEAR EQUATIONS.
C     USING XLM>0 GIVES A MARQUARDT-LIKE PROCEDURE (NOT USED HERE)
C
      CALL LINE(ND,N,XLM,GD,G,W)
      SMAX=0.
      GP=0.
      DP=0.
      GPS=0.
      GGG=0.
      DO 60 I=1,N
      IF (I.GT.NACT) GOTO 50
      IV=IVAR(I)
      GGG=GGG+GD(IV)**2
   50 CONTINUE
      S=W(I,3)
      IF (DABS(S).GT.SMAX) SMAX=DABS(S)
      GP=GP+S*GD(I)
      GPS=GPS+GD(I)**2
   60 DP=DP+S*S*W(I,1)
      FAC=1.
      IF (SMAX.GT.XLAM) FAC=XLAM/SMAX
C
C     XLAM IS THE MAXIMAL PERMITTED STEP, IF LARGER, REDUCE
C
      DER=FAC*GP
      ALFA=1.
   70 FF=ALFA*FAC
      DO 80 I=1,N
      X(I)=W(I,4)+FF*W(I,3)
   80 CONTINUE
C
C     REDUCTION LESS THAN E-4: CHANGE SIGN
C     REDUCTION LESS THAN E-5: ACCEPT AND EXIT
C
      IF (DABS(ALFA).LT. 1.D-4) ALFA=-ALFA
      IF (DABS(ALFA).LT. 1.D-5) GO TO 90
      CALL FUNC(N,1,FNEW,GD,G,X)
      DELS=FNEW-F
C
C     TEST FOR BETTER RESULT (WORSE BY E-10 ACCEPTED)
C
      IF (DELS.LT.1.D-10) GO TO 90
      IF (NTS.GT.1) GO TO 120
C
C     ENERGY INCREASED, TRY WITH SMALLER STEP
C
      ALFA=ALFA/3
      GO TO 70
C
C     CALCULATE PREDICTED IMPROVEMENT (QUADRATIC APPROXIMATION)
C
   90 PRED=FF*GP-FF**2/2*(GP+DP)
      CALPRE=DELS/(PRED+1.D-19)
      GNW=0.
      DO 100 I=1,N
      IF (I.GT.NACT) GOTO 100
      IV=IVAR(I)
      GNW=GNW+GD(IV)**2
  100 CONTINUE
      RG=GNW/GPS
      GPS=GNW
      IF (GPS.LT. GLIM) GO TO 170
      F=FNEW
      DO 110 I=1,N
  110 W(I,4)=X(I)
      IF (NTS.GE.5) GOTO 120
C
C     REPEAT WITH UNMODIFIED HESSIAN UP TO 5 TIMES, IF RESULTS
C     ARE SUFFICIENTLY PROMISING
C
      IF (RG.LT..02.AND. CALPRE.GT..3) GOTO 40
  120 DO 130 I=1,N
  130 X(I)=W(I,4)
      IF (IPR.GT.1) WRITE(6,160)
	iOut = 0
      IF (IPR.GT.1) WRITE(IOUT,160)
  160 FORMAT(' ')
      IF (NTM.LT. NMAX.AND. GPS.GT.GLIM) GO TO 10
      XLAM=GPS/GLIM
170   RETURN
      END
C-----------------------------------------------------------------------
C
C
C     FSTAB EVALUATES FUNCTION, GRADIENT AND HESSIAN FOR STABILITY
C     ANALYSIS USING TANGENT PLANE DISTANCE
C
      SUBROUTINE FSTAB(N,IV,FUN,GRAD,XMAT,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GRAD(45),XMAT(45,45),Y(60),YEX(15),FG(15)
      COMMON/DERIV/FX(15,15)
      COMMON/KEEP/NF,ND,T,P,GNUL,ZZZ(15),FR(15),VL(15,4),SFS(4),ZVAL(4)
      COMMON/AKTIV/IVAL(15),IVAR(45),NACT
      DO 10 I=1,N
C-----TRACE COMPONENTS BY DIRECT SUBSTITUTION
      IF (IVAL(I).EQ.0) Y(I)=FR(I)-FG(I)
   10 CONTINUE
      IF (IV.NE.2) GO TO 25
      NACT=0
      IS=0
C
C     ELIMINATE TRACE COMPONENTS
C
      DO 20 I=1,N
      IVAL(I)=1
      IF (Y(I).LT.-15.) IVAL(I)=0
      IF (IVAL(I).EQ.0) GOTO 20
      IS=IS+1
      IVAR(IS)=I
   20 NACT=NACT+IVAL(I)
   25 SUM=0.
      DO 30 I=1,N
      YEX(I)=0.
      IF (Y(I).GT.-20.) YEX(I)=DEXP(Y(I))
   30 SUM=SUM+YEX(I)
      JD=1
      IF (IV.EQ.2) JD=3
      CALL SRK(JD,0,IC,T,P,Z,YEX,FG,FG,FX)		!equilibrio
      FUN=1.
      DO 40 I=1,N
      S=Y(I)+FG(I)-FR(I)
      IF (IV.EQ.0) GOTO 40
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
C-----------------------------------------------------------------------
C
C
C     SPLIT FERFORMS THE DECOMPOSITION OF THE HESSIAN FOR MURML
C     IDI IS THE NUMBER OF DIAGONAL ELEMENT CHANGES
C     POSITIVE IDI SHOWS AN INDEFINITE HESSIAN
C
      SUBROUTINE SPLIT(ND,N,IDI,BETA,DEL,G,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(ND,ND),W(ND,5)
      COMMON/AKTIV/IVAL(15),IVAR(45),NACT
      IDI=0
      DO 90 I=1,NACT
      I1=I-1
      I2=I+1
      ID=0
      TV=G(I,I)
      SV=TV
      IF (I1.EQ.0) GOTO 20
      DO 10 J=1,I1
   10 SV=SV-G(I,J)**2
   20 IF(SV.LT.DEL*DEL) ID=1
      SVR=DSQRT(DABS(SV))
      IF(SVR.LT.DEL)SVR=DEL
      XM=0.
      IF (I.EQ.NACT) GOTO 60
      DO 50 J=I2,NACT
      S=G(J,I)
      IF (I1.EQ.0) GOTO 40
      DO 30 K=1,I1
   30 S=S-G(I,K)*G(J,K)
   40 S=S/SVR
      IF (DABS(S).GT.XM) XM=DABS(S)
   50 G(J,I)=S
   60 IF (XM.LT.BETA) GOTO 80
      ID=1
      XM=XM/BETA
      DO 70 J=I,NACT
   70 G(J,I)=G(J,I)/XM
      SVR=SVR*XM
   80 IF (ID.EQ.1) W(I,1)=SVR**2-SV
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
C-----------------------------------------------------------------------
C
C
C     LINE CALCULATES THE CORRECTION BY BACKSUBSTITUTION
C
      SUBROUTINE LINE(ND,N,XLAM,GD,G,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GD(ND),G(ND,ND),W(ND,5)
      COMMON/AKTIV/IVAL(15),IVAR(45),NACT
      DO 10 I=1,NACT
      IV=IVAR(I)
   10 W(I,5)=GD(IV)
C----- PREPARED FOR MARQUARDT-LIKE VARIANT,BUT HERE XLAM=0
   20 W(1,2)=-W(1,5)/(G(1,1)+XLAM)
      IF (NACT.EQ.1) GOTO 50
      DO 40 I=2,NACT
      I1=I-1
      S=-W(I,5)
      DO 30 J=1,I1
   30 S=S-G(I,J)*W(J,2)
   40 W(I,2)=S/(G(I,I)+XLAM)
   50 W(NACT,3)=W(NACT,2)/(G(NACT,NACT)+XLAM)
      IF (NACT.EQ.1) GOTO 80
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
C-----------------------------------------------------------------------
C
C
C     SRK IS THE INTERFACE FOR THE THERMODYNAMIC PACKAGE
C     REPLACEMENT OF SRK PERMITS INTERFACING WITH DIFFERENT
C     THERMODYNAMIC PACKAGES, PROVIDED THE ARGUMENTS CONFORM.
C     NUMERICAL DERIVATIVES CAN BE USED FOR TESTING OF NEW THERMO-
C     DYNAMICS. FOR PRODUCTION RUNS, THEY SHOULD BE REPLACED WITH
C     ANALYTIC DERIVATIVES TO AVOID AN UNNECESSARY WASTE OF TIME
C
      SUBROUTINE SRK(INDIC,MTYP,IC,T,P,Z,X,FUG,FUGP,FUGX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(15),FUGT(15),FUG(15),FUGP(15),FUGX(15,15)
      COMMON/STYR/ N,NG,NST,NDER,NTEMP
      COMMON/OLDZ/ZOLD1,ZOLD2
      COMMON/OLDT/TOLD
C
C     NOTE THAT TEMPERATURE DERIVATIVES ARE NOT EVALUATED. IT
C     MAY BE OF INTEREST FOR THE USER TO EVALUATE THE TEMPERATURE
C     DERIVATIVES OF THE FUGACITY COEFFICIENTS (NTEMP=1), SINCE
C     THIS PERMITS A CALCULATION OF OTHER RESIDUAL PROPERTIES DI-
C     RECTLY, E.G. ENTHALPY AND ENTROPY
C
C	Se modificó nTemp
C
C
      NTEMP=0
      NDER=1
      IF (INDIC.GT.2) NDER=2
      TCR=DABS(TOLD-T)
      IF(TCR.GT.1.D-6) CALL  PARAGC(T,N,NG,NST,1)
      TOLD=T
      IGZ=0
      CALL GCEOS(N,NG,NST,NDER,NTEMP,T,P,X,FUG,FUGX,FUGP,FUGP,Z,IGZ,
     @MTYP,IC)
      RETURN

      END
C-----------------------------------------------------------------------
C
C
C     ONE DIMENSIONAL EIGENVECTOR SEARCH ROUTINE
C     ******************************************
C     THIS SUBROUTINE IS A VERY POWERFULL TOOL FOR LOCATING CLOSE
C     MINIMA (I.E. CRITICAL BEHAVIOUR), BUT THE ROUTINE HAS NOT YET
C     BEEN EXTENSIVELY TESTED. EXCEPT FOR TIME CONSUMPTION, IT WILL
C     HOWEVER NEVER DO ANY HARM. IT CAN BE DECOUPLED IN A SIMPLE
C     MANNER
C
      SUBROUTINE ONEDIR(N,RESV,YBEST,FX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YBEST(60),FX(15,15)
      COMMON/STORE/Z(15),ZR(15),Y(15),YD1(15),YD2(15),YV(15),
     *YD(15),FG(15),F1(15),FP(15),YDER(15),YNEW(15),CORR(15),DIST(4),
     *YVAL(15),GRD(15),GDG(15),IP(15)
      COMMON/KEEP/NF,ND,T,P,G0,ZFEED(15),FR(15),XVL(15,4),SF(4),ZVAL(4)
      COMMON/MOLFRC/ZFRAC(15,4)
      COMMON/SOLID/FUGAR(15),ZL(15),IS(15),ISTAT
      RESV=0.
      NFF=NF
      IF (ISTAT.EQ.1) NFF=NF-1
      IF (NFF.EQ.0) RETURN
      DO 410  KKK=1,NFF
C
C     SET UP PHASE COMPOSITION, SOLID NOT TESTED
C
      DO 20 I=1,N
      Z(I)=ZFRAC(I,KKK)
   20 ZR(I)=DSQRT(Z(I))
C
C     EVALUATE PARTIAL DERIVATIVES
C
      CALL SRK(3,0,IC,T,P,ZV,Z,FG,FG,FX)
      DO 30 I=1,N
      F1(I)=FR(I)
      RVAL=ZR(I)
      DO 30 J=I,N
      FX(I,J)=RVAL*FX(I,J)*ZR(J)
   30 FX(J,I)=FX(I,J)
C
C     ESTIMATE EIGENVALUE BY DIRECT ITERATION
C
      XMAX=0.
      DO 40 I=1,N
      XMT=DABS(FX(I,I))
      IF (XMT.LT.XMAX) GO TO 40
      XMAT=XMT
      IMAX=I
   40 CONTINUE
      SSQ=0.
      DO 50 I=1,N
      YD(I)=FX(I,IMAX)
   50 SSQ=SSQ+YD(I)**2
      SSQS=DSQRT(SSQ)
      DO 60 I=1,N
   60 YD(I)=YD(I)/SSQS
      DO 90 M=1,3
      PR1=0.
      PR2=0.
      DO 80 I=1,N
      TERM=0.
      DO 70 J=1,N
   70 TERM=TERM+FX(I,J)*YD(J)
      PR1=PR1+TERM*YD(I)
      PR2=PR2+TERM*TERM
   80 YV(I)=TERM
      PR2S=DSQRT(PR2)
      DO 90 I=1,N
   90 YD(I)=YV(I)/PR2S
      RAT=PR2/PR1
      RAT1=RAT+1
      WRITE(6,100) KKK,RAT1

      iOut = 0

      WRITE(IOUT,100) KKK,RAT1
  100 FORMAT(/,'  PHASE INVESTIGATED ',I3,' EIGENVALUEESTIM.',F9.3)
C
C     STOP TESTING IF ESTIMATE IS NOT PROMISING
C     OTHERWISE REFINE BY INVERSE ITERATION
C
      IF (RAT .GT. -.7) GO TO 410
      DO 110 I=1,N
  110 FX(I,I)=FX(I,I)+1
      CALL LU(15,N,IP,FX)
      RATO=0.
      DO 150 K=1,10
      DO 120 I=1,N
  120 YV(I)=YD(I)
      CALL BACK(15,N,IP,FX,YV)
      SP1=0.
      SP2=0.
      DO 130 I=1,N
      YV(I)=YV(I)-YD(I)
      SP1=SP1+YV(I)*YD(I)
  130 SP2=SP2+YV(I)**2
      RAT=SP2/SP1
      RTEST=1/(RAT+1)
      SP2=DSQRT(SP2)
      DO 140 I=1,N
  140 YD(I)=YV(I)/SP2
      IF( DABS(RTEST-RATO) .LT. 1.D-6) GO TO 160
      RATO=RTEST
  150 CONTINUE
  160 RAT=1/(RAT+1)
C
C     START LINESEARCH
C
      YVS=0.
      YVSQ=0.
      DO 170 I=1,N
      YV(I)=YD(I)/ZR(I)*.5
  170 YVSQ=YVSQ+YV(I)**2
      YVS=0.
      IPOS=0
C
C     SEARCH FOR 'CLOSE' MINIMUM LIMITED TO DISTANCE OF 1
C     ALONG EIGENVECTOR IN DIRECTION OF NEGATIVE THIRD DERIVATIVE
C
      INEG=0
      FUNMIN=0.
      GYOLD=1.
      NCN=0
      SMAX=0.98
C------SEQUENCE TERMINATES AT 0.985, IF NOTHING IS FOUND
C------MINIMUM DISTANCE 0.02 MAY MISS A VERY,VERY CLOSE MINIMUM
      SNUL=.02D0
      DS=1.5*SNUL
      IG=0
      IST=0
      S=SNUL
      GO TO  190
  180 SOLD=S
      FOLD=FUNM
      GYOLD=GY
      IF (S.EQ.SMAX) GO TO 370
      S=S+DS
      IF (S.GT.SMAX) S=SMAX
      IG=1
      IF (S.EQ.0) GO TO 380
  190 LOC=0.
C-----TO CHECK WHEN CALCULATION AND REFINEMENT ARE OF OPPOSITE SIGN
  195 LOC=LOC+1
      SYVAL=0.
      IF (IG.EQ.0) IST=IST+1
      DO 210 I=1,N
      SY=S*YV(I)
      IF (SY.LT.0.) GO TO 200
      YNEW(I)=(ZR(I)*(1+SY))
      GO TO 210
C------AVOID NEGATIVE MOLE FRACTIONS AT LARGE DISTANCES
  200 FAC=(1-SY)*(1+SY*SY)
      YNEW(I)=ZR(I)/FAC
  210 CONTINUE
      DO 220 I=1,N
      IF (IG.EQ.0) GO TO 220
C------CORR EVALUATES A QUADRATIC CORRECTION TO SEARCH DIRECTION
      TTT=S*S*CORR(I)
      YN=YNEW(I)+TTT
      YNEW(I)=YN
C------REPLACE VERY SMALL QUANTITIES,USING DIRECT SUBSTITUTION
      IF (YN.LE.0. .OR.DABS(YN/TTT).LT..5) YNEW(I)=DEXP(
     *(F1(I)-FP(I))/2)
  220 YVAL(I)=YNEW(I)**2
      DO 221 L=1,NF
      DIS=0.
      DO 222 I=1,N
  222 DIS=DIS + (YVAL(I)-ZFRAC(I,L))**2
  221 DIST(L)=DIS
C
C     TEST DISTANCE IN COMPOSITION SPACE TO SOLUTIONS LOCATED ALLREADY
C     LIMIT SEARCH TO AVOID OVERLAP
C
      DO 223 L=1,NF
      IF (DIST(L)/DIST(KKK) .LT. .5) SMAX=S
  223 CONTINUE
      CALL SRK(1,0,IC,T,P,ZV,YVAL,FP,FP,FX)
      XPR=0.
      FUN=1.
      YDP=0.
      GDP=0.
      DO 250 I=1,N
      DDS=F1(I)-FP(I)
      IF (DDS.GT. -14.) GO TO 230
      YNEW(I)=DEXP(DDS/2)
      YVAL(I)=YNEW(I)**2
      TM=0
      GO TO 240
  230 TM=DLOG(YVAL(I))-DDS
  240 FUN=FUN+YVAL(I)*(TM-1)
      TM=DLOG(YVAL(I))+FP(I)-F1(I)
      QQ=YNEW(I)*TM
      GRD(I)=QQ
      GDG(I)=QQ
  250 XPR=XPR+YD(I)*QQ
      GY=XPR/S/RAT
C
C     INVESTIGATE SIGN OF THIRD DERIVATIVE AT ORIGO
C     CHANGE DIRECTION, IF SIGN IS WRONG
C
      TN=(GY-1)*RAT/S
      IF (IG.EQ.1) GO TO 270
      IF (TN*S .LT. 0.) GO TO 270
      IF (IST.GE.2) GO TO 270
      DO 260 I=1,N
      YV(I)=-YV(I)
  260 YD(I)=-YD(I)
      GO TO 190
  270 CONTINUE
      DO 280 I=1,N
  280 GDG(I)=GDG(I)-XPR*YD(I)
      CALL BACK(15,N,IP,FX,GDG)
      GAIN=0.
      DO 290 I=1,N
  290 GAIN=GAIN+GDG(I)*(GRD(I)+XPR*YD(I))
      GAIN=GAIN/2
C
C     INCLUDE CORRECTION TERM FOR NEXT STEP
C
      DO 300 I=1,N
      CRN=YNEW(I)-.5*GDG(I)
      IF (IG.EQ.0 ) CORR(I)=-GDG(I)/SNUL**2/2
      IF (IG.GT.0) CORR(I)=CORR(I)-GDG(I)/S/S/2
  300 CONTINUE
      FUNM=FUN-GAIN
C
C     ONLY ACCEPT VERIFIED, NOT SUGGESTED RESULTS, BUT REPEAT
C     THE CALCULATION ONCE WITH THE CORRECTED VALUE, IF THE SUGGES-
C     TED RESULT IS NEGATIVE
C
      IF (LOC.EQ.1 .AND. FUNM*FUN .LT. 0.) GO TO 195
      IF (FUN.GE.FUNMIN) GO TO 320
      FUNMIN=FUN
      DO 310 I=1,N
  310 YBEST(I)=YVAL(I)
  320 CONTINUE
C------GY IS THE PROCECTED GRADIENT,WHICH DETERMINES THE SEARCH
      IF (RAT.LT.0.) GY=-GY
      IF (IG.EQ.0) GO TO 180
      IF (INEG.EQ.1) GO TO 340
      IF (GY.LT.0.) GO TO 330
      IF (GY.GT.GYOLD) GO TO 380
      DS=1.5*DS
      GO TO 180
  330 INEG=1
      GYNEG=GY
      SNEG=S
      DS=1.5*DS
      FNEG=FUNM
      GO TO 180
  340 IF (GY.LT.0. .AND. IPOS.EQ.0) GO TO 330
      IPOS=1
C
C     STATIONARY POINT BRACKETED;ONLY 3 FURTHER STEPS TAKEN.
C
      NCN=NCN+1
      IF (NCN.GT.3) GO TO 380
      IF (GY.GT.0.) GO TO 350
      GYNEG=GY
      SNEG=S
      FNEG=FUNM
      GO TO 360
  350 GYPOS=GY
      SPOS=S
      FPOS=FUNM
  360 RV=(GYPOS-GYNEG)/(SPOS-SNEG)
C-----INVERSE LINEAR INTERPOLATION USED IN SEARCH
      S=SNEG-GYNEG/RV
      DS=0
      GO TO 180
  370 CONTINUE
  380 CONTINUE
  390 FORMAT(/,'  SEARCH CONCLUSIVE, FMIN = ',D10.2)
      IF (FUNMIN .GE. 0.) GO TO 410
      DO 400 I=1,N
  400 YBEST(I)=DLOG(YBEST(I))
      WRITE(6,390 ) FUNMIN
      WRITE(IOUT,390 ) FUNMIN
      RESV=FUNMIN
C	IF (FUNMIN.LE.-1.D-5) RETURN		!alfo?
      IF (FUNMIN.LE.-1.D-7) RETURN
  410 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C       PITON EVALUATES PHASE DISTRIBUTION, GIVEN THE FUGACITY COEF-
C       FICIENTS, INCLUDING SEARCH FOR APPEARING OR VANISHING SOLID
C       AND FLUID PHASES
C
      SUBROUTINE PITON(N,H,VEC,ALFA,YVAL,XJAC,IPX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALFA(15,4),BETEX(15),IPX(6),G(6),H(15),XJAC(6,6)
      DIMENSION YVAL(60),VEC(45),BOLD(6),BTOLD(15)
      LOGICAL NOSOL
      COMMON/SOLID/FUGAR(15),ZLOG(15),IS(15),ISTAT
      COMMON/KEEP/NF,ND,T,P,GNUL,Z(15),FR(15),XVL(15,4),BETA(4),ZV(4)
C-----CALCULATION OF PHASE DISTRIBUTION, SOLID PHASES
      IOUT=0
      ITER=0
      NOSOL=(ISTAT.EQ.0)
      NR=NF-1
      IF (NOSOL) NR=NF
      IF (NOSOL) GOTO 205
C------  K-FACTORS, SOLIDS PRESENT
      IF (NR.GT.1) GOTO 205
C------  ONLY SINGLE PHASE;
      DO 202 I=1,N
      DD=0.
      ALFA(I,1)=1.
      IF (IS(I).EQ.0) GOTO 202
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
      IF (J.EQ.NR) SS=0.
      SS=SS-SMAX
      IF (SS.LT.-60.) SS=-60.
      IF (SS.EQ.0.) ALFA(I,J)=1.
      IF (SS.NE.0.) ALFA(I,J)=DEXP(SS)
  214 CONTINUE
      IF (NOSOL) GOTO 210
      IF (IS(I).GT.0.) GOTO 216
      ALFA(I,NF)=0.
      GOTO 210
  216 SS=VEC(I+(NR-1)*N)
      SS=SS-SMAX
      IF (SS.LT. -60.) SS=-60.
      ALFA(I,NF)=DEXP(SS)
  210 CONTINUE
C------START PHASE ITERATION
  300 IF (NOSOL) GOTO 100
C-------EVALUATE CURRENT STATE OF SOLIDS
      DO 310 I=1,N
      IF (IS(I).EQ.0) GOTO 310
      HI=0.
      DO 320 J=1,NR
  320 HI=HI+ALFA(I,J)*BETA(J)
      BETX=Z(I)-HI/ALFA(I,NF)
      IS(I)=1
      IF (BETX.GT.0.) IS(I)=2
      BTOLD(I)=BETX
  310 CONTINUE
  100 ITER=ITER+1
      DO 10 I=1,N
      HI=0.
      IF (IS(I).LT.2 .OR. NOSOL) GOTO 15
C-------MOLE FRACTIONS FIXED BY SOLID EQUIL.
      HI=Z(I)*ALFA(I,NF)
      GOTO 10
   15 DO 20 J=1,NR
   20 HI=HI+ALFA(I,J)*BETA(J)
      IF (HI.EQ.0.) HI=1.
   10 H(I)=1/HI
C------CONVERGENCE TEST
      IF (IOUT.EQ.1) GOTO 350
C-----INITIATE JACOBIAN
      DO 30 I=1,NR
      G(I)=1.
      DO 30 J=1,I
   30 XJAC(I,J)=0.
C-----MAIN LOOP, GRADIENT AND JACOBIAN
      DO 40 I=1,N
      ZI=Z(I)
      HI=H(I)
      IST=IS(I)
      DO 45 J=1,NR
      ZL=ALFA(I,J)*HI
      PR=ZL*ZI
      G(J)=G(J)-PR
      IF (BETA(J).EQ.0.) GOTO 45
      IF (IST.EQ.2 .AND. ISTAT.GT.0) GOTO 45
      DO 50 K=1,J
      ZK=ALFA(I,K)*HI
   50 XJAC(J,K)=XJAC(J,K)+ZK*PR
   45 CONTINUE
   40 CONTINUE
C------CHECK GRADIENT AND COPY JACOBIAN
      DO 60 J=1,NR
      IF (BETA(J).GT.0.) GOTO 70
      IF (G(J).LT.0.) GOTO 65
C-----GRADIENT POSITIVE; SKIP OKAY
      DO 61 K=1,NF
      XJAC(J,K)=0.
   61 XJAC(K,J)=0.
      XJAC(J,J)=1.
      G(J)=0.
      GOTO 70
C------RE-INITIATE PHASE, GRADIENT NEGATIVE
   65 BETA(J)=1.D-16
      GOTO 100
C------SAFEGUARD AGAINST ZERO DETERMINANT OF JACOBIAN
   70 XJAC(J,J)=XJAC(J,J)+1.D-7
      DO 60 K=1,J
   60 XJAC(K,J)=XJAC(J,K)
C------ CALCULATE NEWTON STEP
      CALL LU(6,NR,IPX,XJAC)
      CALL BACK(6,NR,IPX,XJAC,G)
C-------MONITOR STEPLENGTH TO AVOID NEGATIVE PHASE FRACTIONS
      FACL=0.
      FACR=1.
      JZM=0
      DO 80 J=1,NR
      DA=G(J)
      DAB=DABS(DA)
      IF (DAB.GT.FACL) FACL=DAB
      GNEW=BETA(J)-DA
      IF (GNEW.GE.0.) GOTO 80
      RAT=BETA(J)/DAB
      IF (RAT.GT.FACR) GOTO 80
      FACR=RAT
      JZM=J
   80 CONTINUE
C------EVALUATE REDUCTION FACTOR
      RANG=5*FACL
      DEC=1.
      IF (RANG.GT.1.) DEC=1/RANG
      IF (JZM.EQ.0) GOTO 90
      IF (FACR.LT. DEC) GOTO 85
      JZM=0
      GOTO 90
   85 DEC=FACR
   90 DO 110 J=1,NR
      BOLD(J)=BETA(J)
      G(J)=G(J)*DEC
      BETA(J)=BETA(J)-G(J)
      IF (J.EQ.JZM) BETA(J)=0.
  110 CONTINUE
C------END FLUID PHASES; CHECK STATUS OF SOLIDS; FIND REDUCTION
C------CORRESPONDING TO ONE CHANGE
  350 IF (NOSOL) GOTO 400
      ICHAN=0
      RMIN=1.
      DO 120 I=1,N
      IST=IS(I)
      IF (IST.EQ.0) GOTO 120
      HI=0.
      DO 125 J=1,NR
  125 HI=HI+BETA(J)*ALFA(I,J)
      BETX=Z(I)-HI/ALFA(I,NF)
      IF (IOUT.EQ.1) GOTO 130
      BOLDI=BTOLD(I)
      IF (BETX*BOLDI .GT.0.) GOTO 130
      ICHAN=ICHAN+1
      DIF=BETX-BOLDI
      RAT=DABS(BOLDI/DIF)
      IF (RAT.LE. RMIN) RMIN=RAT
  130 BETEX(I)=BETX
  120 CONTINUE
      IF (IOUT.EQ.1) GOTO 400
      IF (ICHAN.GT.0) GOTO 140
      DO 135 I=1,N
      IF (IS(I).EQ.0) GOTO 135
      BTOLD(I)=BETEX(I)
  135 CONTINUE
      GOTO 400
  140 RT1=1.01*RMIN
      IF (RT1.GT. .99) RT1=(1+RAT)/2
      RT2=1-RT1
      DO 141 J=1,NR
  141 BETA(J)=RT1*BETA(J)+RT2*BOLD(J)
      DO 145 I=1,N
      IF (IS(I).EQ.0) GOTO 145
      BETEX(I)=RT1*BETEX(I)+RT2*BTOLD(I)
      IF (BETEX(I)*BTOLD(I) .LT.0.) IS(I)=3-IS(I)
      BTOLD(I)=BETEX(I)
  145 CONTINUE
  400 IF (IOUT.EQ.1) GOTO 150
      IF (ITER.GE.20 .OR. FACL.LT. 1.D-8) IOUT=1
      GOTO 100
C-----CONVERGENCE;
  150 IF (ISTAT.GT.0.) BETA(NF)=0.
      DO 160 I=1,N
      HI=H(I)
      DO 165 J=1,NR
      BET=BETA(J)
      IF (BET.LE.0.) BET=-1.
  165 XVL(I,J)=BET*ALFA(I,J)*HI
      IF (NOSOL) GOTO 160
      IF (IS(I).LT.2) XVL(I,NF)=0.
      IF (IS(I).EQ.2) XVL(I,NF)=BETEX(I)/Z(I)
      BETA(NF)=BETA(NF)+XVL(I,NF)
  160 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C       PHASN EVALUATES FUNCTION, GRADIENT AND HESSIAN FOR PHASE SPLIT
C       CALCULATIONS. PHASN IS AN ARGUMENT OF THE MURRAY ROUTINE
C
      SUBROUTINE PHASN(NARG,IVAR,FUN,GRAD,XMAT,YVAL)

      IMPLICIT REAL*8(A-H,O-Z)

	parameter(NCM=15,NFM=4)
      DIMENSION INDEX(NCM),NCON(NCM,NFM),ENER(NFM)
      DIMENSION YVAL(NCM*NFM),GRAD(NCM),X(NCM),TM(NCM),FG(NCM)
      DIMENSION XMAT(45,45),XMIN(NCM),DELEN(NFM)
      COMMON/MOLFRC/FILD(NCM,NFM)
      COMMON/DERIV/FX(NCM,NCM)
      COMMON/KEEP/NF,ND,T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),
     *            ZVAL(NFM),indF(NFM)
      COMMON/SOLID/FUGAR(NCM),ZLOG(NCM),ISol(NCM),IStat
      COMMON/OPTIM/VBEST(NCM*NFM),GMIN
      COMMON/AKTIV/IVAL(NCM),IPOS(45),NACT
      COMMON/LIMIT/EPSV
      JD=1
      N=NARG/(NF-1)
   10 NG=NF-1
      NFE=NF-ISTAT
      NARG=N*NG
      IF(IVAR.LT.2) GOTO 110
      JD=3
      SMIN=1.
C
C     SECOND DERIVATIVES NEEDED.
C     IDENTIFY TRACE COMPONENTS AND DEPENDENT VARIABLES
C
      DO 20 J=1,NF
      IF (SFAS(J).GT.SMIN) GO TO 20
      SMIN=SFAS(J)
      JMIN=J
   20 CONTINUE
C
C     ARRAY NCON: 1: INDEPENDENT; 0: DEPENDENT; -1: TRACE
C     INDEX KEEPS TRACE OF DEPENDENT VARIABLES
C
      DO 70 I=1,N
      IVAL(I)=1
      XMAX=0
      XM1=0.
      DO 30 J=1,NF
      IV=I+(J-1)*N
      YVAL(IV)=VBEST(IV)
      XVL(I,J)=VBEST(IV)
      XM=YVAL(IV)*ZZ(I)/SFAS(J)
      IF(XM.GT.XM1)XM1=XM
      NCON(I,J)=1
      IF(XM.LT.EPSV)NCON(I,J)=-1
      IF(J.EQ.JMIN)XMIN(I)=XM
      IF(YVAL(IV).LT.XMAX) GOTO 30
      XMAX=YVAL(IV)
      JMAX=J
   30 CONTINUE
      INDEX(I)=JMAX
      NCON(I,JMAX)=0
      IG=0
      IF (XM1.GT.EPSV) GOTO 40
      IVAL(I)=0
      NCON(I,JMAX)=-1
      GO TO 70
   40 DO 60 J=1,NF
      IF (NCON(I,J).GE.0) GO TO 50
      IF (XVL(I,J)*ZZ(I)/SFAS(JMAX).GT. EPSV) NCON(I,J)=1
   50 IF (NCON(I,J).GT. 0) IG=1
   60 CONTINUE
      IF (IG.EQ.0) IVAL(I)=0
   70 CONTINUE
C
C     ORDER ACTUAL INDEPENDENT VARIABLES
C
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
      IF (IVAL(I).EQ.0) GOTO 90
      IJ=I+JN
C----- NOTE CHANGE FROM PREVIOUS INCORRECT VERSION
      JJ=J
      IF (J.EQ.INDEX(I)) JJ=NF
      IF (NCON(I,JJ).LT.0) GOTO 90
      NACT=NACT+1
      IPOS(NACT)=IJ
   90 CONTINUE
  100 CONTINUE
      IF (NACT.EQ.0) NACT=1
      GOTO 200
C
C     REARRANGE VARIABLES
C
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
      IF(NCON(I,J).EQ.-1) GOTO 140
  130 CONTINUE
      GOTO 160
 140  BFAS=SFAS(IV)
      IF (IV.GT.NFE) BFAS=XVL(I,NF)*ZZ(I)
      DO 150 J=1,NF
      TMJ=1
      IF (J.EQ.IV) GOTO 155
      SFJ=SFAS(J)
      IF (J.GT.NFE) SFJ=XVL(I,J)
      IF (J.LE.NFE .OR. ISOL(I).EQ.2) GOTO 152
      TMJ=0
      GOTO 155
  152 DSN=FILD(I,IV)-FILD(I,J)
      IF (DSN.LT. -80.) DSN=-80.
      TMJ=SFJ/BFAS*DEXP(DSN)
  155 TM(J)=TMJ
  150 RSM=RSM+TMJ
  160 SM=0.
C
C     CORRECT TRACE COMPONENTS, USING OLD K-FACTORS AND PHASE AMOUNTS
C
      DO 170 J=1,NF
      IF (J.EQ.IV) GOTO 170
      XVL(I,J)=YVAL(I+(J-1)*N)
      IF (XVL(I,J).EQ.0.) GOTO 175
      IF(NCON(I,J).EQ.-1)XVL(I,J)=TM(J)/RSM
  175 XVLV=XVLV-XVL(I,J)
  170 CONTINUE
      XVL(I,IV)=XVLV
      DO 180 J=1,NF
      SLIMT=1.D-8
      IF (J.EQ.NF .AND. ISTAT.EQ.1) SLIMT=0.
      IF (XVL(I,J).LT.0.) XVL(I,J)=SLIMT
  180 SM=SM+XVL(I,J)
      DO 190 J=1,NF
      XVL(I,J)=XVL(I,J)/SM
  190 YVAL(I+(J-1)*N)=XVL(I,J)
C------INITIATE OBJECTIVE FUNCTION
  200 FUN=-GNUL
      IF (IVAR.EQ.0) GOTO 230
C------INITIATE GRADIENT
      DO 220 I=1,NARG
      GRAD(I)=0.
      IF (IVAR.LT.2) GOTO 220
C------INITIATE HESSIAN
      DO 210 J=1,NARG
  210 XMAT(I,J)=0.
  220 CONTINUE
C
C     START CALCULATING CONTRIBUTION FROM INDIVIDUAL PHASES
C
  230 DO 350 K=1,NF
      SFAS(K)=0.
      ENER(K)=0.
      DO 240 I=1,N
      X(I)=XVL(I,K)*ZZ(I)
  240 SFAS(K)=SFAS(K)+X(I)
      IF (K.EQ.NF .AND. ISTAT.EQ.1) GOTO 245
      CALL SRK(JD,0,IC,T,P,ZVAL(K),X,FG,FG,FX)
C
	indF(k) = iC
C
  245 DELEN(K)=0
      DO 250 I=1,N
      IF (K.NE.NF .OR. ISTAT.NE.1) GOTO 255
      IF (ISOL(I).LT. 2) FILD(I,K)=100.
      IF (ISOL(I).EQ.2) FILD(I,K)=FUGAR(I)
      TM(I)=0.
      TMI=0
      IF (ISOL(I).EQ.2) TMI=FUGAR(I)
      TMI=TMI-ZLOG(I)
      TM(I)=TMI
      GOTO 248
  255 FILD(I,K)=FG(I)
      TM(I)=FG(I)+DLOG(XVL(I,K)/SFAS(K))
      TMI=TM(I)
  248 DELEN(K)=DELEN(K)+XMIN(I)*TMI
  250 ENER(K)=ENER(K)+TMI*X(I)
      FUN=FUN+ENER(K)
C
C     END GIBBS ENERGY, START GRADIENT
C
      IF (IVAR.EQ.0) GOTO 350
      DO 280 I=1,N
      S=TM(I)*ZZ(I)
      IPL=K
      IV=INDEX(I)
      IF (K.EQ.NF) IPL=IV
      IF (K.EQ.IV) GOTO 260
      IPS=I+(IPL-1)*N
      GRAD(IPS)=GRAD(IPS)+S
      GOTO 280
  260 DO 270 M=1,NG
      IPS=I+(M-1)*N
  270 GRAD(IPS)=GRAD(IPS)-S
  280 CONTINUE
C
C     END GRADIENT, START HESSIAN MATRIX
C
      IF (IVAR.LT.2) GOTO 350
      IF (K.EQ.NF .AND. ISTAT.EQ.1) GOTO 350
      DO 300 I=1,N
      S=ZZ(I)/SFAS(K)
      DO 290 J=1,I
      FX(I,J)=S*(FX(I,J)-1.)*ZZ(J)
  290 FX(J,I)=FX(I,J)
  300 FX(I,I)=FX(I,I)+ZZ(I)/XVL(I,K)
      DO 340 I=1,N
      IV=INDEX(I)
      IPL=K
      IF (K.EQ.NF) IPL=IV
      DO 330 J=1,N
      JV=INDEX(J)
      JPL=K
      IF (K.EQ.NF) JPL=JV
      S=FX(I,J)
      M11=IPL
      M12=IPL
      M21=JPL
      M22=JPL
      IF (IV.NE.K) GOTO 310
      S=-S
      M11=1
      M12=NG
  310 IF (JV.NE.K) GOTO 320
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
      IF (IVAR.NE.1) GOTO 370
C
C     CHANGE BACK TO NATURAL ORDER
C
      DO 360 I=1,N
      IV=INDEX(I)
      S=YVAL(I+(IV-1)*N)
      YVAL(I+(IV-1)*N)=YVAL(I+(NF-1)*N)
      YVAL(I+(NF-1)*N)=S
  360 CONTINUE
  370 CONTINUE
      IF (IVAR.LT.2.OR.NF.LE.2) GOTO 440
C
C     CHECK FOR REMOVAL OF A PHASE
C
      DMIN=1.
      DO 380 J=1,NF
      IF (J.EQ.JMIN) GOTO 380
      DELT=DELEN(J)-DELEN(JMIN)
      IF (DELT.GT.DMIN) GOTO 380
      K=J
      DMIN=DELT
  380 CONTINUE
      IF (DMIN.EQ.1.) GOTO 440
      J=JMIN
      DO 390 I=1,N
  390 X(I)=ZZ(I)*(XVL(I,J)+XVL(I,K))
      SFT=SFAS(J)+SFAS(K)
      CALL SRK(1,0,IC,T,P,ZNEW,X,FG,FG,FX)
      ENGO=ENER(J)+ENER(K)
      DO 400 I=1,N
  400 ENGO=ENGO-X(I)*(FG(I)+DLOG(X(I)/SFT/ZZ(I)))
      IF (ENGO.LT.-1.D-10) GOTO 440
      IF (J.EQ.NF) ISTAT=0
      IF (K.LT.J) GOTO 410
      JJJ=J
      J=K
      K=JJJ
  410 CONTINUE
C
C     STORE BEST RESULT SO FAR IN VBEST
C
      DO 420 I=1,N
      SS=XVL(I,K)+XVL(I,J)
      XVL(I,K)=SS
      XVL(I,J)=XVL(I,NF)
      VBEST(I+(K-1)*N)=SS
  420 VBEST(I+(J-1)*N)=VBEST(I+(NF-1)*N)
      NF=NF-1
      WRITE(6,430)
	iOut = 0
      WRITE(IOUT,430)
  430 FORMAT(/,'  NUMBER OF PHASES REDUCED ')
      GOTO 10
  440 IF (FUN.GT.GMIN+1.D-10) GOTO 460
C
C     NO REDUCTION, STORE BEST RESULT
C
      GMIN=FUN
      DO 450 I=1,N
      DO 450 J=1,NF
      IV=I+(J-1)*N
  450 VBEST(IV)=XVL(I,J)
  460 RETURN
      END
C-----------------------------------------------------------------------
C
C       DSUB CALCULATES PHASE SPLIT BY DIRECT SUBSTITUTION, USING
C       THE GENERAL DOMINANT EIGENVALUE METHOD OF CROWE AND NISHIO
C       FOR ACCELERATION
C
      SUBROUTINE DSUB(N,FX,YVAL)
      IMPLICIT REAL*8(A-H,O-Z)
	parameter(NCM=15,NFM=4)
      DIMENSION YVAL(NCM*NFM),FX(NCM,NCM),IP(6)
C
C     COMMON AREA STORE IS USED FOR INTERNAL VARIABLES IN THE DSUB
C     ROUTINES TO SAVE SPACE
C
      COMMON/KEEP/NF,ND,T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),
     *            ZVAL(NFM),indF(NFM)
      COMMON/STORE/GRAD(45),WORK(45,5),SFOLD(NFM),B(6,6),Y(45),H(NCM),
     *             FG(NCM),YOLD(45),YNEW(45),CR(NCM,NFM),BAU(6),
     *             YDAR(NCM,NFM)
      COMMON/OPTIM/VBEST(NCM*NFM),GMIN
      COMMON/NTEST/NFA,NDUM,GNORM,GLIM
      COMMON/SOLID/FUGAR(NCM),ZLOG(NCM),IS(NCM),ISTAT
      COMMON/NEWSOL/NEWVAL
C
C     COLUMN DIMENSION OF CR MUST BE INCREASED TO MAFX-1 WITH MORE
C     THAN 4 PHASES
C
C     GDEM USED IS OF ORDER M-2, WHERE M MUST BE < 6, HERE 4
C
	iOut = 0
      M=4
      LMAX=4
    5 NFA=NF
      NG=NF-1
      IF (NG.GT.0) GOTO 6
      GNORM=0.
      GOTO 440
    6 IV=1
      NT=N*NG
C
C     SUBROUTINE NTOPI FINDS FUGACITY COEFFICIENTS FROM COMPOSITIONS.
C     SUBROUTINE PITON FINDS PHASE AMOUNTS (AND COMPOSITION) FROM
C     K-FACTORS
C
      CALL NTOPI(N,GTOT,FG,H,FX,YOLD,YDAR)
      GMIN=GTOT
      GBAS=GMIN
      DO 20 I=1,NT
   20 YNEW(I)=YOLD(I)
      DO 230 L=1,LMAX
C
C      ISO AND ICOR ARE FAILURE COUNTERS
C
      ISO=0
      ICOR=0
   40 CONTINUE
      ISO=ISO+1
      M2=M+IV
      DO 140 KK=1,M2
      KV=KK-IV
      DO 50 J=1,NF
   50 SFOLD(J)=SFAS(J)
      CALL PITON(N,H,YOLD,CR,YVAL,B,IP)
      GOTO 100
   70 DO 80 J=1,NF
   80 SFAS(J)=SFOLD(J)
      DO 90 I=1,NT
C-------RESTART FROM LAST VALUE, IF EXTRAPOLATE IS WORSE
      IF (ISO.GT.2.OR.L.EQ.LMAX) YOLD(I)=YNEW(I)
   90 YOLD(I)=(YNEW(I)+YOLD(I))/2
      ICOR=ICOR+1
      IF (L.EQ.LMAX.OR.ICOR.GT.4) GOTO 230
      GOTO 40
  100 CONTINUE
      CALL NTOPI(N,GTOT,FG,H,FX,GRAD,YDAR)
C
C     CHECK FOR IMPROVEMENT
C
      IF (KK.NE.1.OR.GTOT.LE.GMIN+1.D-12) GOTO 110
      GOTO 70
  110 IF (GTOT.GT.GMIN+1.D-12) GOTO 130
      GMIN=GTOT
      DO 120 I=1,N
      DO 120 J=1,NF
  120 YVAL(I+(J-1)*N)=XVL(I,J)
C-------TEST FOR SOLIDS AFTER EXTRAPOLATION
      IF (KK.NE.1 .OR. ISTAT.EQ.1) GOTO 125
      IF (NEWVAL.EQ.0) GOTO 125
C-------SOLID PHASE FOUND
      WRITE(6,126)
      WRITE(IOUT,126)
  126 FORMAT(/,'  SOLID PHASE ADDED IN DSUB ',/)
      NF=NF+1
      ISTAT=1
      SFAS(NF)=0.
      DO 128 I=1,N
  128 XVL(I,NF)=0.
      GOTO 5
  125 CONTINUE
      IRFAS=0
  115 JOUT=0
      DO 111 J=1,NF
C-----A VALUE IN THIS RANGE INDICATES THAT THE PHASE CAN BE REMOVED
      IF (SFAS(J).EQ.0.) JOUT=J
  111 CONTINUE
      IF (JOUT.NE.0.) IRFAS=1
      IF (IRFAS.EQ.0) GO TO 119
      IF (JOUT.EQ.0 .AND. IRFAS.EQ.1) GOTO 5
      NF=NF-1
      WRITE(6,116)
      WRITE(IOUT,116)
      IF (JOUT.GT.NF) ISTAT=0
      IF (JOUT.GT.NF) GO TO 115
C-----RESET YIELD ARRAYS TO SMALLER NUMBER OF PHASES
      DO 112 J=JOUT,NF
      SFAS(J)=SFAS(J+1)
      DO 112 I=1,N
      XVL(I,J)=XVL(I,J+1)
  112 YVAL(I+(J-1)*N)=YVAL(I+J*N)
  116 FORMAT(/,'  PHASE REMOVED IN DSUB ',/)
      GO TO 115
  119 CONTINUE
      IF (GNORM.LT.GLIM) GOTO 440
      IF (L.EQ.LMAX) GOTO 230
C
C     SAVE RESULT OF DIRECT SUBSTITUTION
C
  130 DO 140 I=1,NT
      IF (KV.GT.0) WORK(I,KV)=GRAD(I)
      Y(I)=GRAD(I)-YOLD(I)
  140 YOLD(I)=GRAD(I)
      M1=M-1
      DO 160 I=1,M1
      DO 160 J=I,M1
      B(I,J)=0.
      DO 150 K=1,NT
      SS=(WORK(K,I+1)-WORK(K,I))*(WORK(K,J+1)-WORK(K,J))
  150 B(I,J)=B(I,J)+SS
  160 B(J,I)=B(I,J)
C------CONVERGENCE; DANGER OF UNDERFLOW IN EXTRAPOLATION
      IF (B(M1,M1).LT.1.D-20) GOTO 440
      DO 170 J=1,M1
      B(M1,J)=0.
  170 BAU(J)=0.
      BAU(M1)=1.
      B(M1,M1)=1.
C
C     SOLVE LINEAR EQUATIONS FOR EXTRAPOLATION
C
      CALL LU(6,M1,IP,B)
      CALL BACK(6,M1,IP,B,BAU)
      SNY=0.
      DO 180 J=1,M1
  180 SNY=SNY+BAU(J)
      SLIM=.5
      CPR=0.
      SMAX=SLIM
C
C     CHECK THAT EXTRAPOLATE IS IN CORRECT DIRECTION
C
      DO 200 I=1,NT
      YNEW(I)=YOLD(I)
      YOLD(I)=0.
      DO 190 J=1,M1
  190 YOLD(I)=YOLD(I)+BAU(J)*WORK(I,J+1)
      YOLD(I)=YOLD(I)/SNY
      DS=YOLD(I)-YNEW(I)
      IF (DABS(DS).GT.SMAX) SMAX=DABS(DS)
  200 CPR=CPR+DS*Y(I)
C-----OUTPUT ENERGY AND COEFFICIENTS OF CHARACTERISTIC EQUATION
      WRITE(6,210) GMIN,(BAU(J),J=1,M1)
      WRITE(IOUT,210) GMIN,(BAU(J),J=1,M1)
  210 FORMAT(5X,D15.5,3X,6F12.7)
      IF (CPR.GT.0..AND.SMAX.LE.SLIM) GOTO 230
C
C     LIMIT EXTRAPOLATION FACTOR TO IMPROVE STABILITY
C
      FAC=SLIM/SMAX
C------CHANGE SIGN, IF IN WRONG DIRECTION
      IF (CPR.LT.0.) FAC=-FAC
      DO 220 I=1,NT
  220 YOLD(I)=YNEW(I)+FAC*(YOLD(I)-YNEW(I))
  230 CONTINUE
  240 CONTINUE
C
C     AND THIS IS WHERE THE BROYDEN METHOD USED TO BE
C
C------SAVE BEST RESULTS SEEN SO FAR
  440 IF (NEWVAL.LE.ISTAT) GOTO 445
      WRITE(6,126)
	iOut = 0
      WRITE(IOUT,126)
      NF=NF+1
      ISTAT=1
      DO 442 I=1,N
  442 XVL(I,NF)=0.
      SFAS(NF)=0.
      GOTO 5
  445 DO 450 I=1,N
      DO 450 J=1,NF
      VBEST(I+(J-1)*N)=YVAL(I+(J-1)*N)
  450 XVL(I,J)=YVAL(I+(J-1)*N)
      WRITE(6,460)GMIN,GNORM
      WRITE(IOUT,460)GMIN,GNORM
  460 FORMAT(/,'  END DSUB  GMIN =',D25.15,' GNORM =',D15.3/)
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C     NTOPI EVALUATES FUGACITY COEFFICIENTS FROM PHASE COMPOSITIONS
C     AND ALSO TESTS FOR CONVERGENCE IN DSUB
C
      SUBROUTINE NTOPI(N,GTOT,FG,H,FX,VECT,YR)
      IMPLICIT REAL*8(A-H,O-Z)

	parameter(NCM=15,NFM=4)

      DIMENSION FG(NCM),H(NCM),FX(NCM,NCM),VECT(45),X(NCM),YR(NCM,NFM)
      COMMON/SOLID/FUGAR(NCM),ZLOG(NCM),ISOLC(NCM),ISOLPC
      COMMON/KEEP/NF,ND,T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),
     *            ZVAL(NFM),indF(NFM)
      COMMON/NTEST/NFA,NDUM,GNORM,GLIM
      COMMON/NEWSOL/NEWFS

      NEWFS=0
      GNORM=0.
      GTOT=-GNUL
      NREF=NF
C------ PHASE NREF IS REFERENCE PHASE FOR K-FACTORS
      IF (ISOLPC.NE.0) NREF=NF-1
      NG=NREF-1
      SF=0.
C-------PROPERTIES OF REFERENCE PHASE
      DO 10 I=1,N
      XI=ZZ(I)*XVL(I,NREF)
      IF (SFAS(NREF).EQ.0.) XVL(I,NREF)=0.
      X(I)=XI
   10 SF=SF+XI
      IF (SF.GT.0.) SFAS(NREF)=SF
C------OBJECTIVE FUNCTION EVALUATION
      CALL SRK(1,0,IC,T,P,ZVAL(NREF),X,FG,FG,FX)
C
	indF(nRef) = iC
C
      DO 20 I=1,N
      YI=FG(I)+DLOG(X(I)/ZZ(I)/SF)
      YR(I,NREF)=YI
      IF (SF.LE.0.) GOTO 20
      GTOT=GTOT+X(I)*YI
   20 CONTINUE
      IF (NG.EQ.0) GOTO 50
      DO 40 J=1,NG
      SF=0.
      DO 30 I=1,N
      XI=XVL(I,J)*ZZ(I)
      IF (SFAS(J).EQ.0.) XVL(I,J)=0.
      X(I)=XI
   30 SF=SF+XI
      IF (SF.GT.0.) SFAS(J)=SF
      CALL SRK(1,0,IC,T,P,ZVAL(J),X,H,H,FX)
C
	indF(j) = iC
C
      DO 40 I=1,N
      HI=H(I)
      FAC=HI+DLOG(X(I)/SF/ZZ(I))
      IF (SFAS(J).EQ.0.) GOTO 41
      GTOT=GTOT+X(I)*FAC
      YR(I,J)=FAC
      IF (X(I).LT. 1.D-15) GOTO 41
      GNORM=GNORM+(ZZ(I)*(FAC-YR(I,NREF)))**2
C
C     LOG K-FACTORS STORED IN VECT
C
   41 CONTINUE
   40 VECT(I+(J-1)*N)=FG(I)-HI
   50 IF (ISOLPC.EQ.0) GOTO 80
      SFAS(NF)=0.
      DO 60 I=1,N
      IPS=(NF-2)*N+I
      VECT(IPS)=0.
      IF (ISOLC(I).EQ.0) GOTO 60
      ZZI=ZZ(I)
      XI=XVL(I,NF)*ZZI
      SFAS(NF)=SFAS(NF)+XI
      HI=FUGAR(I)
      VECT(IPS)=FG(I)-HI
      FAC=HI-ZLOG(I)
      GTOT=GTOT+XI*FAC
      IF (XI.LT.1.D-15) GOTO 60
      IF (ISOLC(I).EQ.2) GNORM=GNORM+(ZZI*(FAC-YR(I,NREF)))**2
   60 CONTINUE
      RETURN
   80 IF (GNORM.GT. 1.D-8) RETURN
      DO 90 I=1,N
      IF (ISOLC(I).EQ.0) GOTO 90
      VALUE=FUGAR(I)-ZLOG(I)
      JM=0
      DO 100 J=1,NF
      IF (YR(I,J).GE. VALUE) GOTO 100
      JM=1
  100 CONTINUE
      IF (JM.GT.0) GOTO 90
      ISOLC(I)=2
      NEWFS=1
   90 CONTINUE
      RETURN
      END
C
C     SRKPUR EVALUATES PURE COMPONENT PROPERTIES AND SOLID FUGACITY
C     COEFFICIENTS FOR USE IN STABILITY ANALYSIS
C
      SUBROUTINE SRKPUR(JPUR,OLDY,FX,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(15),OLDY(17,15),FX(15,15),FG(15)
      COMMON/SOLPRP/TMELT(15),HMELT(15)
      COMMON/SKIP/NOTRY(15)
      COMMON/SOLID/FUGAR(15),ZLOG(15),ISOLC(15),ISOLPC
      COMMON/KEEP/NF,ND,T,P,GNUL,ZZ(15),FR(15),XVL(15,4),SFAS(4),ZVAL(4)
      IPT=0
      FMAX=0.
      DO 30 I=1,N
      ISOLC(I)=1
      IF (T.GE. TMELT(I)) ISOLC(I)=0.
      IF (ISOLC(I).GT.0) IPT=1
      IF (ISOLC(I).EQ.0 .AND. NOTRY(I).EQ.0) GOTO 30
      DO 10 J=1,N
   10 A(J)=0.
C
      A(I)=1.
C
C     PURE COMPONENT PROPERTIES ARE HERE CALCULATED FROM FULL N-COMPO-
C     NENT THERMODYNAMICS. IT IS, OF COURSE MORE EFFICIENT SO SUBSTI-
C     TUTE A SEPARATE ROUTINE (EASY FOR CUBIC EQUATIONS)
C
      CALL SRK(1,1,IC,T,P,ZZZ,A,FG,FG,FX)
C
C     LOG FUGACITY (COEFFICIENT) DIFFERENCE STORED
C
      IF (P.LT. 1. .AND. IC.NE.1) ISOLC(I)=0.
C
C-----EVALUATE PROPERTIES OF SOLID PHASE FOR BELOW-MELTING COMP.
C
      TRED=TMELT(I)/T
      FUGAR(I)=FG(I)-HMELT(I)*(TRED-1)
      DO 20 J=1,N
   20 FX(I,J)=FG(J)
   30 CONTINUE
C     IF (IPT.GT.0) GOTO 36
C     GOTO 48
C  36 WRITE(6,35)
C  35 FORMAT(/,'  SOLID FORMATION INVESTIGATION FOR COMPONENTS ',//,
C    *'       T MELT   LN (FUG/P)    ',/)
C     DO 40 I=1,N
C     IF (ISOLC(I).EQ.0) GOTO 40
C     WRITE(6,45) I,TMELT(I),FUGAR(I)
C  45 FORMAT(1X,I3,2X,F7.2,2X,F9.5)
C  40 CONTINUE
   48 SMAX=0.
C------ NOTE REVISION FROM EARLIER VERSION
      DO 50 I=1,N
      IF (NOTRY(I).EQ.0) GOTO 50
      DO 55 J=1,N
   55 OLDY(I,J)=FR(J)-FX(I,J)
      IF (OLDY(I,I).LT. SMAX) GOTO 50
      SMAX=OLDY(I,I)
      JPUR=I
   50 CONTINUE
      IF (JPUR.EQ.0) RETURN
      SMAX=SMAX+3
      DO 60 I=1,N
      IF (NOTRY(I).EQ.0) GOTO 60
      IF (I.NE. JPUR) OLDY(I,JPUR)=OLDY(I,JPUR)-SMAX
   60 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C     SUBROUTINE PURE INITIATES THE STABILITY ANALYSIS, SCANNING TRIAL
C     PHASES BASED ON PURE COMPONENT PROPERTIES FOR INSTABILITY
C
      SUBROUTINE PURE(Y,N,NH,FX)
      IMPLICIT REAL*8(A-H,O-Z)
	parameter(NCM=15,NFM=4)

      DIMENSION LIST(NCM+2),RTOLD(NCM+2),OLDY(NCM+2,NCM)
      DIMENSION Y(NCM),YTEST(NCM),FA(NCM),DS(NCM),DSG(NCM),FX(NCM,NCM)
      COMMON/KEEP/NF,ND,T,P,GNUL,ZZ(NCM),FR(NCM),XVL(NCM,NFM),SFAS(NFM),
     *            ZVAL(NFM),indF(NFM)
      COMMON/DELST/YV1
      COMMON/SKIP/NOTRY(NCM)
      COMMON/SOLID/FUGAR(NCM),ZLOG(NCM),ISOLC(NCM),ISOLPC
      COMMON/MOLFRC/ZFRAC(NCM,NFM)

	iOut = 0

      JPUR=0
      JMIN=0
      KMIN=1
      NT=N+1
C
C     LIST IDENTIFIES ELEMENTS, WHERE SEARCH IS SKIPPED
C
      DO 10 I=1,NT
      IF (NF.EQ.1) LIST(I)=1
      IF (LIST(I).NE.0) LIST(I)=1
   10 CONTINUE
      IF (NF.EQ.1) GOTO 40
      NT=N+2
      ZMAX=0.
      JMAX=0
      DO 20 J=1,NF
      IF (ZVAL(J).LT. ZMAX) GOTO 20
      ZMAX=ZVAL(J)
      JMAX=J
   20 CONTINUE
C
C     SELECT ADDITIONAL PRIMER PHASE IN MULTIPHASE ANALYSIS
C
      DO 30 I=1,N
      YCM=ZFRAC(I,JMAX)
      YPN=(3*ZZ(I)+2*YCM)/5
   30 OLDY(NT-1,I)=DLOG(YPN)
      LIST(NT-1)=1
   40 LIST(NT)=1
      IF (NF.NE.1) GOTO 45
      IF (ISOLPC.EQ.0) CALL SRKPUR(JPUR,OLDY,FX,N)
C------IDEAL GAS TRIAL PHASE
   45 NA=3
      DO 50 I=1,N
      IF (NOTRY(I).EQ.0) LIST(I)=0
   50 OLDY(NT,I)=FR(I)
      IF (JPUR.LE.0) GOTO 55
      DO 56 I=1,N
   56 Y(I)=1.D-10
      Y(JPUR)=1.
      YV1=-1
      RETURN
   55 DO 140 KK=1,NA
      KMIN=3
      RTMIN=1.D20
      DO 130 K=1,NT
C
C     SEARCH PHASE OR NOT ?
C
      IF (LIST(K).EQ.0 .OR. LIST(K).GE.3) GOTO 130
      IF (LIST(K).EQ.0.OR.LIST(K).GE.3) GOTO 130
C
C     N O T E:  TRIAL PHASES 1 TO N ARE N O T RESTORED
C     ********  IT MAY (OR MAY NOT) BE BETTER TO RESTORE
C
C------SPECIFY TYPE OF PRIMER PHASE
      MV=1
      IF (K.EQ.NT) MV=-1
      SUM=0.
      DO 60 I=1,N
      Z=1.D-25
      IF (OLDY(K,I).GT.-80.) Z=DEXP(OLDY(K,I))
      YTEST(I)=Z
   60 SUM=SUM+Z
      DLSUM=DLOG(SUM)
      CALL SRK(1,MV,IC,T,P,Z,YTEST,FA,FA,FX)
      DMIN=N
C
C     FIND CLOSEST PHASE FOR COMPARATION
C
      NFF=NF
      IF (ISOLPC.EQ.1) NFF=NF-1
      IF (NFF.EQ.0) NFF=1
      DO 90 J=1,NFF
      DIST=0.
      DO 70 I=1,N
      Z=ZFRAC(I,J)
      IF (J.EQ.1) YTEST(I)=YTEST(I)/SUM
      DS(I)=YTEST(I)-Z
   70 DIST=DIST+DS(I)**2
      IF (DIST.GT.DMIN) GOTO 90
      DMIN=DIST
      MINVL=J
      DO 80 I=1,N
   80 DSG(I)=DS(I)
   90 CONTINUE
      VAL=0.
      GRAD=0.
      DO 100 I=1,N
      DNEW=FR(I)-FA(I)
      GD=OLDY(K,I)-DNEW-DLSUM
      OLDY(K,I)=DNEW
      VAL=VAL+YTEST(I)*GD
  100 GRAD=GRAD+GD*DSG(I)
C
C     EVALUATE RESULT; CHECK FOR NEGATIVE VALUE.
C
      RT=1+VAL
      IF (GRAD.GT.0.) RT=2*VAL/GRAD
      IF (VAL.LT.0.) RT=VAL
      IF (DMIN.LT.2.D-3.AND.RT.GT..8) LIST(K)=0
      IF (KK.EQ.1) GOTO 110
      IF (RT.LT.0.) GOTO 110
      IF (LIST(K).EQ.0) GOTO 130
      IF (RT.LT.RTOLD(K)) GOTO 110
      LIST(K)=LIST(K)+1
      IF (RT.GT.1.) LIST(K)=LIST(K)+1
  110 RTOLD(K)=RT
      IF (LIST(K).GE.3) GOTO 130
      IF (RT.GT.RTMIN) GOTO 130
      KMIN=LIST(K)
      RTMIN=RT
      JMIN=MINVL
      DO 120 I=1,N
  120 Y(I)=YTEST(I)
  130 CONTINUE
C
C      SEARCH CONCLUSIVE FOR NEGATIVE RT
C
      IF (RTMIN.LT.0.) GOTO 150
  140 CONTINUE
  150 YV1=0.
      IF (RTMIN.LT.0.) YV1=RTMIN
      IF (YV1.LT.-1.D-7) GOTO 240
C
C     SEARCH INCONCLUSIVE. TRY ALSO REFLECTION WITH RESPECT TO THE
C     CLOSEST "TRIVIAL" SOLUTION
C
      DO 160 I=1,N
      DS(I)=ZFRAC(I,JMIN)
      Z=DS(I)/Y(I)
      IF (Z.GT.2.) Z=2.
      IF (Z.LT..5) Z=.5
  160 YTEST(I)=DS(I)*Z
      DO 180 K=1,NA
      CALL SRK(1,1,IC,T,P,Z,YTEST,FA,FA,FX)
      IF (K.EQ.NA) GOTO 180
      DO 170 I=1,N
      Z=FR(I)-FA(I)
      IF (Z.LT.-80.) Z=-80.
  170 YTEST(I)=DEXP(Z)
  180 CONTINUE
      SUM=0.
      DO 190 I=1,N
  190 SUM=SUM+YTEST(I)
      VAL=0.
      GRAD=0.
      DO 200 I=1,N
      YTEST(I)=YTEST(I)/SUM
      GD=DLOG(YTEST(I))+FA(I)-FR(I)
      VAL=VAL+YTEST(I)*GD
  200 GRAD=GRAD+GD*(YTEST(I)-DS(I))
      IF (VAL.LT.0.) GOTO 210
      RT=1+VAL
      IF (GRAD.GT.0.) RT=2*VAL/GRAD
      IF (KMIN.EQ.0.OR.KMIN.EQ.3) GOTO 210
      RTMIN=(1+KMIN)*RTMIN/2
      IF (RT.LT.RTMIN) GOTO 210
      GOTO 240
  210 IF (VAL.LT.0.) YV1=VAL
      DO 230 I=1,N
  230 Y(I)=YTEST(I)
	iOut = 0
  240 WRITE(IOUT,250)(Y(I),I=1,N)
  250 FORMAT(/,' COMPOSITION FOR STABILITY TEST ',/,(5X,12F9.6))
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C       OUT PRINTS RESULT OF (INTERMEDIATE) PHASE SPLIT CALCULATION
C       NOTE THAT SOLID PHASES ARE COMBINED IN THE LAST PHASE, IF
C       PRESENT, AND THAT THEIR "MOLE FRACTIONS" INDICATE THEIR RELA-
C       TIVE AMOUNTS
C
      SUBROUTINE OUT(YVAL,N,FX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION YVAL(60),DLX(15),FX(15,15)
      COMMON/SOLID/FUGAR(15),ZLOG(15),ISOLD(15),ISTAT
      COMMON/KEEP/NF,ND,T,P,GNUL,Z(15),FR(15),XVL(15,4),SFAS(4),ZVAL(4)
      COMMON/MOLFRC/ZFRAC(15,4)
	iOut = 0
      NT=NF*N
      NB=NT-N
      IF (ISTAT.EQ.0) GOTO 15
      WRITE(6,11) NF
      WRITE(IOUT,11) NF
      ZVAL(NF)=0.
   11 FORMAT(/,'  NOTE:   PHASE  ',I2,'  IS SOLID ',/)
C
C     CHANGE OUTPUT FORMATS FOR MORE THAN 4 PHASES
C

   15 WRITE(IOUT,20)NF
   20 FORMAT(/,'  FLASH RESULT ',I3,' PHASES * * * * ')
      WRITE(6,30)(J,SFAS(J),J=1,NF)
      WRITE(IOUT,30)(J,SFAS(J),J=1,NF)
   30 FORMAT(/,'  TOTAL AMOUNT  ',4(2X,I2,F12.8))
      WRITE(6,40)(J,ZVAL(J),J=1,NF)
      WRITE(IOUT,40)(J,ZVAL(J),J=1,NF)
   40 FORMAT('  COMPRESS      ',4(2X,I2,F12.8))
      WRITE(6,50)(J,J=1,NF)
      WRITE(IOUT,50)(J,J=1,NF)
   50 FORMAT(/,'  COMPOSITION   ',4(5X,I2,9X))
      SUM=0.
      NFF=NF
      IF (ISTAT.EQ.1) NFF=NF-1
      IF (NFF.NE.0) GOTO 52
      DO 51 I=1,N
   51 FR(I)=FUGAR(I)
      GOTO 75
   52 DO 60 I=1,N
      DLX(I)=XVL(I,NFF)*Z(I)/SFAS(NFF)
   60 SUM=SUM+DLX(I)
      SUM=DLOG(SUM)
      CALL SRK(1,0,IC,T,P,ZZ,DLX,FR,FR,FX)
      DO 70 I=1,N
      DLX(I)=DLOG(DLX(I))
   70 FR(I)=DLX(I)+FR(I)-SUM
   75 DO 90 I=1,N
      DO 80 J=1,NF
   80 ZFRAC(I,J)=XVL(I,J)*Z(I)/SFAS(J)
   90 WRITE(IOUT,100)I,(ZFRAC(I,J),J=1,NF)

  100 FORMAT('  X(',I2,')         ',5(1X,G15.7))

      RETURN
      END
C**************************************************************
