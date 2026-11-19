!------------------------------------------------------------------------------------------
!
!
!
!
!
!
!	GC DB con  Asociación cruzada y de autoasociación (Abril 2006) -ALFONSINA-
!
!
!	Modificado para incluir paquete GCA 2.0 (calc. asoc. propuestos por Michelsen
!		Francisco, 2011-2012
!--------------------------------------------------------------------------------------
!	Versión adaptada para ser incluida en ParaEst, para reemplazar a ATXFLASH()
!	No estoy eliminando código que no hace nada. Solo adaptando para que fun-
!	cione.
!		Francisco, 03/2012.
!
      subroutine GCADB(MODEL,NC,NG,NST,iTy,T,P,xx,yy,iErrBurb,ZLiq,ZVap)
!
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 K(15),X(15),Y(15),GG(15,16),Z(15)
      DIMENSION NAME(36)
!
      COMMON/CRIT/TC(15),PC(15)
      COMMON/CONTL/ERROR,IPRINT,IEXPON,INIT,ITYPE,MAXIT,ITHERM,IT
      COMMON/COUT/IOUT

!
      COMMON/Fact_comp/ZL,ZV
!
      MAXIT=20
5001  FORMAT(20I3)
5003  FORMAT(36A2)
6003  FORMAT(///,1X,36A2)
!	OPEN(UNIT=5,STATUS='OLD',FILE='Gcdba.dat')

!	Alfonsina
!	open(unit=8,file='datter.dat',status='unknown',form='formatted')
!	write(8,'("   T(K)        P(atm)          y(i)")')
!
!	READ(5,*) ISW,IPRINT,IEXPON,INIT
      iSw = 0
      iPrint = 0
      iExpon = 8
      init = 0
      iType = iTy
!	IF (ISW.EQ.1) THEN

!	  OPEN(UNIT=3,FILE='GCDBA.sal')

!	ENDIF
      IF(IEXPON.EQ.0) IEXPON=6
!	READ(5,*) MODEL,NC,NG,NGA
      WRITE(6,6001) ISW
      IF(ISW.EQ.1) WRITE (3,6001)ISW
6001  FORMAT(1H1,///,1X,'IOUT =',I2,/,'If IOUT = 0: output only on unit'&
     &       ' 6',/,' If IOUT = 1: output on both unit 6 and .sal',//)
!
!	READ(5,5003) NAME
      N=NC
      NOCOMP=NC
      IOUT=ISW
!	WRITE(6,6003) NAME
!	IF(ISW.EQ.1) WRITE(3,6003) NAME
!20	IF(MODEL.GT.1) GO TO 25
20    if(model.EQ.1)then

        ITHERM=1
!	  WRITE(6,6007)
!	  IF(IOUT.EQ.1) WRITE(3,6007)
!6007	  FORMAT(//,' Model used for calculations: GCA-EoS',//)
!	  CALL PARGR(NC,NG,NGA,NST,iOut)
!	  CALL PARMOL(NC,NG,NST,IOUT)
!	  GO TO 30

      elseif(model.EQ.2)then
!25	IF(MODEL.EQ.3) GO TO 27

        WRITE(6,6008)
        IF(IOUT.EQ.1) WRITE(3,6008)
6008    FORMAT(//,' Model used for calculations: SRK-EoS',//)
!	  GO TO 29

      elseif(model.EQ.3)then

27      WRITE(6,6009)
        IF(IOUT.EQ.1) WRITE(3,6009)
6009    FORMAT(//,' Model used for calculations: PR-EoS',//)
29      CALL INDAT_cubica(N,NH)
        ITHERM=2

      endif
30    NCALC=0
      WRITE(6,3050)
      IF(IOUT.EQ.1) WRITE(3,3050)
!
!40	READ(5,*,ERR=700,END=700) (Z(I),I=1,NC)                           00005600
!
!	Lectura en una sola línea, al estilo GC3
40    continue
!	if(init.EQ.0)then

!	  READ(5,*,ERR=700,END=700)(Z(I),I=1,NC),T,P,iType

!	else

!	  READ(5,*,ERR=700,END=700)(Z(I),I=1,NC),T,P,iType,(K(I),I=1,NC)

!	endif
!
!      IF(Z(1).LE.0.0) GO TO 700                                         00005610
      K(:NC) = 0.
      Z(1) = XX
      Z(2) = 1.D0 - XX
      SumZ = sum(z(:NC))
      if(SumZ.LE.0D0)goto 700
      do i = 1,NC

        if(z(i).LE.0.D0)z(i) = 1.D-8

      enddo
!45	READ(5,*,ERR=700,END=700) ITYPE,P,T
!	IF(INIT.GT.0) READ(5,ERR=700,END=700)(K(I),I=1,NC)
      WRITE(6,3000)
      IF(IOUT.EQ.1) WRITE(3,3000)
3000  FORMAT(//,80('-'),/)

!	IF(ITYPE.EQ.0) GO TO 40
      NCALC=NCALC+1
!
!	WRITE(6,3100) NCALC,ITYPE
!	IF(IOUT.EQ.1) WRITE(3,3100) NCALC,ITYPE
!
      CALL SUMA(N,Z,SZ,1)
      IF(INIT.EQ.0)then		!NE.0) GO TO 80

        CALL KINIT(NC,X,Y,Z,T,P)
!	  DO 75 I=1,NC
!75	  X(I)=Z(I)
        x(:NC) = z(:NC)

      endif
80    IF(INIT.NE.0)then		!EQ.0) GO TO 100

        IF(INIT.EQ.1) CALL SUMA(N,K,SK,1)
        DO 90 I=1,NC

          X(I)=Z(I)
          IF(INIT.EQ.1) Y(I)=K(I)

90      enddo

100   endif
!
!	Si se calculan P, calcula propiedades dependiente de T (=cte)
!	ahora. Si se calculan T, como van variando hasta converger, se
!	recalculan en cada iteración y no acá.
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) CALL PARAGC(T,NC,NG,NST,1)
!
!	Llamado al punto de saturación. iType entra por COMMON aparentemente.
      NC1=NC+1
      CALL DBTP(NC,NC1,NG,NST,T,P,X,Y,K,GG,ICH)
      ZLiq = ZL
      ZVap = ZV
!
!	Variable de error:
      IF(ICH.EQ.99)then
!
!	  La siguiente modificación nace para tener una entrada de una sola línea por
!	  punto de composición corrido:
!	  GO TO 45
500     WRITE(6,3600)
        IF(IOUT.EQ.1) WRITE(3,3600)
3600    FORMAT(//,' Calculation stopped: initial estimate out off ',&
     &         'convergence region',//)
!	  goto 40

        iErrBurb = 1

      endif

!	IF(ITYPE.EQ.2.OR.iType.EQ.4)then		!EQ.1.OR.ITYPE.EQ.3) GO TO 300

!       Se calcularon putos de rocío. Por lo tanto, la salida (Y) es líquido
!	  y la alimentación (Z) es vapor.
!	  DO 200 I=1,NC

!	    X(I)=Y(I)
!	    Y(I)=Z(I)

!200	  enddo

!	endif

      YY = Y(1)
!
300   WRITE(6,3200) IT,P,T,ERROR
      IF(IOUT.EQ.1) WRITE(3,3200) IT,P,T,ERROR
3200  FORMAT(/,' ITER=',I3,2X,' Pressure(ATM)=',F10.4,2X,' Temperature='&
     &,F10.4,' Error=',D12.4,/)
!
      WRITE(6,3300)
!
!	Alfonsina
!	write(8,16)	T,P ,( y(i), i=1,nc )
!16	format((F10.4),1X,(F10.4), 1x,(F10.6), 1x,(F10.6), 1x,(F10.6) )! Alfonsina
16    format((F10.4),1X,(F10.4),<NC>(X,(F10.6)))	  ! Francisco
!
      IF(IOUT.EQ.1) WRITE(3,3300)
!	DO 400 I=1,NC

!       Si se calcularon puntos de burbuja, la salida (Y), es líquido.
!	  IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) X(I)=Y(I)
!	  WRITE(6,3400) I,Z(I),X(I),K(I)
!	  IF(IOUT.EQ.1) WRITE(3,3400) I,Z(I),X(I),K(I)

!400	enddo
      IF(ICH.NE.0) WRITE(6,3500) ICH
      IF(ICH.NE.0.AND.IOUT.EQ.1) WRITE(3,3500) ICH
3500  FORMAT(//,' Error in calcuation:  IER = ',I5,/,' IER = 1 : Trivial&
     & sollution y=x',/,' IER = 10: No convergence in 10 iterations',//)
!      GO TO 45
!
!	goto 40
!
1100  FORMAT(8F10.4)
1200  FORMAT(I3,5F10.4)

3050  FORMAT(//,' ITYPE=1 : Bubble point temperature',/,' ITYPE=2 : DEW&
     & POINT TEMPERATURE',/,' ITYPE=3: Bubble point pressure ',/,' ITYPE&
     &=4: Dew point pressure',/)
3100  FORMAT(/,' Calculation number: ',I5,' , ITYPE= ',I5,//)

3300  FORMAT(/,' Comp.      Feed    Incipient phase     ln(K)',/)
3400  FORMAT(I3,7X,F10.7,2X,F10.7,6X,F10.6)
!
700   CLOSE (UNIT=5)
      IF(IOUT.EQ.1) CLOSE(UNIT=3)
!
      return
      ENDsubroutine
!-------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!
!
!	                                                                  00000010
      SUBROUTINE DBTP(NC,NC1,NG,NST,T,P,X,Y,Z,DG,IER)
!     ******************************************************************00000030
!                                                                       00000040
!     ERROR MESSAGES (IER VALUES)                                       00000050
!     0     : NO ERRORS                                                 00000060
!     1     : THE TRIVIAL SOLUTION  X=Y                                 00000070
!     10    : THE EQUILIBRIUM CALCULATIONS DID NOT CONVERGE IN MAXIT    00000080
!             ITERATIONS                                                00000090
!     11    : ERROR ENCOUNTERED DURING CALCULATION OF DYDX AND DPDX.    00000100
!             THE EQUILIBRIUM CALCULATION IS ALL RIGHT.                 00000110
!                                                                       00000120
!     99,100                                                            00000130
!     199,200: SEE ZMAX                                                 00000140
!                                                                       00000150
!     ******************************************************************00000160
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(15),Y(15),Z(15),G(15),DG(NC,NC1),DEL(15),PHIL(15)
      DIMENSION PHIV(15),FX(15,15),FY(15,15),FPX(15),FPY(15),FTX(15)
      DIMENSION FTY(15)
      COMMON/CONTL/ERROR,IPRINT,NDIG,INIT,ITYPE,MAXIT,ITHERM,NIT
      COMMON/COUT/IOUT

      COMMON/Fact_comp/ZL,ZV
      real(8), parameter :: macheps = epsilon(1D0)

      EPS=10.D0**(-NDIG)
      DXYMAX=0.75D0
      DPMAX=0.75D0
      DTMAX=0.75D0
      NIT=0
      IGZ=0
      IER=0
      IGOLD=IGZ
      NCC=NC-1
      IXY=1
      IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) IXY=2
3     GO TO(100,200,100,200),ITYPE
100   IF(ITHERM.EQ.1.AND.ITYPE.LE.2) CALL PARAGC(T,NC,NG,NST,1)
      IF(ITHERM.EQ.1)CALL GCEOS(NC,NG,NST,2,1,T,P,X,PHIL,FX,FTX,FPX,&
     &ZL,IGZ,1,ICL)
      IF(ITHERM.EQ.1) CALL GCEOS(NC,NG,NST,2,1,T,P,Y,PHIV,FY,FTY,FPY,&
     &ZV,IGZ,-1,ICV)
      IF(ITHERM.EQ.2) CALL SRK1(2,1,ICL,T,P,ZL,X,PHIL,FTX,FPX,FX)
      IF(ITHERM.EQ.2) CALL SRK1(2,-1,ICV,T,P,ZV,Y,PHIV,FTY,FPY,FY)
      GO TO 300
200   IF(ITHERM.EQ.1.AND.ITYPE.LE.2) CALL PARAGC(T,NC,NG,NST,1)
      IF(ITHERM.EQ.1) CALL GCEOS(NC,NG,NST,2,1,T,P,Y,PHIL,FX,FTX,FPX,&
     &ZL,IGZ,1,ICL)
      IF(ITHERM.EQ.1) CALL GCEOS(NC,NG,NST,2,1,T,P,X,PHIV,FY,FTY,FPY,&
     &ZV,IGZ,-1,ICV)
      IF(ITHERM.EQ.2) CALL SRK1(2,1,ICL,T,P,ZL,Y,PHIL,FTX,FPX,FX)
      IF(ITHERM.EQ.2) CALL SRK1(2,-1,ICV,T,P,ZV,X,PHIV,FTY,FPY,FY)
300   IF(ICL.EQ.1.AND.ICV.EQ.1) GOTO 400
      IER=99
      GO TO 99
400   ZCRT=DABS(ZV-ZL)
      IF(ZCRT.LT.1.D-6) IER=1
      IF (IER.GT.0) GOTO 2

      GNORM=0.D0
      DO 4 I=1,NC
        xi = X(I)+macheps
        yi = Y(I)+macheps
        IF(IXY.EQ.1) G(I)=DLOG(xi/yi)-PHIV(I)+PHIL(I)
        IF(IXY.EQ.2) G(I)=DLOG(yi/xi)-PHIV(I)+PHIL(I)
4       GNORM=GNORM+G(I)**2

      GNORM=DSQRT(GNORM)
      ERROR=GNORM
      IF(IPRINT.EQ.0) GO TO 40
      WRITE(6,1000) NIT,P,T,ERROR
      IF(IOUT.EQ.1) WRITE(3,1000) NIT,P,T,ERROR
40    DO 6 I=1,NC
      IF(ITYPE.GT.2) DG(I,NC)=-FPY(I)+FPX(I)
      IF(ITYPE.LE.2) DG(I,NC)=-FTY(I)+FTX(I)
      DO 5 J=1,NCC
      IF(IXY.EQ.1) DG(I,J)=-FY(I,J)+FY(I,NC)
      IF(IXY.EQ.2) DG(I,J)=FX(I,J)-FX(I,NC)
      IF(I.EQ.J.AND.IXY.EQ.1) DG(I,J)=DG(I,J)-1.D0/Y(I)
      IF(I.EQ.J.AND.IXY.EQ.2) DG(I,J)=DG(I,J)+1.D0/Y(I)
      IF(I.EQ.NC.AND.IXY.EQ.1) DG(I,J)=DG(I,J)+1.D0/Y(I)
5     IF(I.EQ.NC.AND.IXY.EQ.2) DG(I,J)=DG(I,J)-1.D0/Y(I)
6     DG(I,NC1)=-G(I)
      IF(GNORM.LT.EPS) GOTO 2
      CALL GAUSL(NC,NC1,NC,1,DG)
      RATMIN=1.D0
      DO 7 I=1,NCC
      CRITY=DXYMAX*Y(I)
      IF(Y(I).GT..5D0) CRITY=.5D0*(1.D0-Y(I))
      DYABS=DABS(DG(I,NC1))
      IF(DYABS.LT.CRITY) GOTO 7
      RAT=CRITY/DYABS
      IF(RAT.LT.RATMIN) RATMIN=RAT
7     CONTINUE
      PCT=DABS(DG(NC,NC1))
      IF(ITYPE.GT.2) PCRIT=DPMAX*P
      IF(ITYPE.LE.2) PCRIT=DTMAX*T
      RAT=1.D0
      IF (PCT.GT.PCRIT) RAT=PCRIT/PCT
      IF(RAT.LT.RATMIN) RATMIN=RAT
      YSUM=0.D0
      DO 9 I=1,NCC
9     YSUM=YSUM+Y(I)+DG(I,NC1)*RATMIN
      CRITY=.5D0*Y(NC)
      YTEST=1.D0-YSUM-Y(NC)
      YTEST=DABS(YTEST)
      IF(YTEST.GT.CRITY) RAT=CRITY/YTEST
      IF(RAT.LT.RATMIN) RATMIN=RAT
      IF(ITYPE.GT.2) P=P+DG(NC,NC1)*RATMIN
      IF(ITYPE.LE.2) T=T+DG(NC,NC1)*RATMIN
      YSUM=0.D0
      DO 8 I=1,NCC
      Y(I)=Y(I)+DG(I,NC1)*RATMIN
8     YSUM=YSUM+Y(I)
      Y(NC)=1.D0-YSUM
      NIT=NIT+1
      IF(NIT.EQ.5) IGZ=IGZ+2
      IF(NIT.GT.MAXIT) IER=10
      IF(IER.NE.0) GOTO 2
      GOTO 3
2     DO 98 I=1,NC
      Z(I)=PHIL(I)-PHIV(I)
98    CONTINUE
99    IGZ=IGOLD
1000  FORMAT(' ITER=',I5,' VARIABLE1=',F10.4,' VARIABLE2=',F10.4,' RESID&
     &UE=',D14.4)
!
      RETURN
      ENDsubroutine
!------------------------------------------------------------------------------------------
!
!
!
!
!	Sum normalizes input mole numbers to sum 1
!
      SUBROUTINE SUMA(N,Z,SZ,IC)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Z(15)

      Sz = sum(z(:N))
      IF (IC.EQ.0 .OR. SZ.EQ.0.) RETURN
      z(:N) = z(:N)/Sz
!
      RETURN
      ENDsubroutine
!-------------------------------------------------------------------------------------------
!
!
!
!
!     SRK1 IS THE INTERFACE ROUTINE TO CUBGEN                           00001350
!                                                                       00001360
      SUBROUTINE SRK1(INDIC,MTYP,IC,T,P,Z,X,FUG,FUGT,FUGP,FUGX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(15),FUGT(15),FUG(15),FUGP(15),FUGX(15,15)
      COMMON/STYR/ N,NH,NDER,NTEMP
      NTEMP=0
      NDER=1
      IF (INDIC.GT.1) NDER=2
      IF (INDIC.EQ.2 .OR. INDIC.EQ.4) NTEMP=1
      CALL CUBGEN(MTYP,IC,T,P,Z,X,FUG,FUGT,FUGP,FUGX)
      RETURN
      END
!                                                                       00001480
!     SUBROUTINE ANEW CALCULATES PARAMETERS A AND B AND THEIR           00001490
!     PARTIAL COMPOSITION AND TEMPERATURE DERIVATIVES                   00001500
!                                                                       00001510
      SUBROUTINE ANEW(X,A0,AT,ATT,B0,AD1,AD2,ADT,BD1)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AD1(15),AD2(120),ADT(15),BD1(15),&
     &X(15),AX0(15),AX1(15),AS1(15),AST1(15)
      COMMON/STYR/N,NH,NDER,NTEMP
      COMMON/PAR/CK(120),BC(15),AC0(15),AC1(15)
      S0=0.
      S1=0.
      DO 10 I=1,N
      XI=X(I)
      AS1(I)=0.
      AX0(I)=XI*AC0(I)
      IF (NTEMP.EQ.0) GO TO 10
      AST1(I)=0.
      AX1(I)=XI*AC1(I)
      S1=S1+AX1(I)
10    S0=S0+AX0(I)
      A0=S0*S0
      NADR=0
15    DO 20 I=1,N
      FAC=2*AC0(I)
      AD1(I)=S0*FAC
      IF (NTEMP.NE.0) ADT(I)=FAC*S1+2*S0*AC1(I)
      IF (NDER.LT.2) GO TO 20
      DO 25 J=I,N
      NADR=NADR+1
25    AD2(NADR)=FAC*AC0(J)
20    CONTINUE
      NADR=0
!                                                                       00001810
!     START MAIN LOOP FOR A-PARAMETER                                   00001820
!                                                                       00001830
      DO 50 I=1,N
      IF(I.EQ.N) GO TO 50
      K1=I+1
      NADR=NADR+1
      XI=X(I)
      AX0I=AX0(I)
      IF (NTEMP.NE.0.) AX1I=AX1(I)
      FAC=2*AC0(I)
      DO 60 K=K1,N
!                                                                       00001930
!     START INNER LOOP                                                  00001940
!                                                                       00001950
      NADR=NADR+1
      EF=CK(NADR)
      IF(EF.EQ.0.) GO TO 60
      AX0K=AX0(K)
      AS1(I)=AS1(I)+AX0K*EF
      AS1(K)=AS1(K)+AX0I*EF
      IF (NTEMP.EQ.0) GO TO 65
      AST1(I)=AST1(I)+AX1(K)*EF
      AST1(K)=AST1(K)+AX1I*EF
65    IF (NDER.LT.2) GO TO 60
      FAC1=AC0(K)*FAC*EF
      AD2(NADR)=AD2(NADR)-FAC1
60    CONTINUE
50    CONTINUE
      A0=0
      AT=0.
      DO 70 I=1,N
      FAC=2*AC0(I)
      AD1(I)=AD1(I)-FAC*AS1(I)
      A0=A0+AD1(I)*X(I)/2
      IF (NTEMP.EQ.0) GO TO 70
      ADT(I)=ADT(I)-2*AC1(I)*AS1(I)-FAC*AST1(I)
      AT=AT+X(I)*ADT(I)
70    CONTINUE
      AT=AT/2
!                                                                       00002210
!     START B-CALCULATION                                               00002220
!                                                                       00002230
      B0=0.
      DO 120 I=1,N
      B0=B0+X(I)*BC(I)
      BD1(I)=BC(I)
120   CONTINUE
      RETURN
      END
!                                                                       00002310
!     CUBIC SOLVES FOR THE COMPRESSIBILITY FACTOR                       00002320
!     MTYP=1: LIQUID; MTYP=-1: VAPOUR; MTYP= 0: MIN. GIBBS ENERGY       00002330
!                                                                       00002340
      SUBROUTINE CUBIC(MTYP,A,B,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CUB/C,C1,C2
      BC=B*C
      X2=1-BC
      X1=A-B*(1+B+C+2*BC)
      X0=B*(A-BC*(1+B))
      Z=X2/3
      F=((Z-X2)*Z+X1)*Z-X0
      Z=B
      IF(F.LT.0..OR. B.GT. X2/3)Z=Z+1
10    DF2=3*Z-X2
      DF1=Z*(DF2-X2)+X1
      F=((Z-X2)*Z+X1)*Z-X0
      DZ=F/DF1
      DZ=DZ*(1+DZ*DF2/DF1)
      Z=Z-DZ
      IF (DABS(DZ).GT. 1.D-10) GO TO 10
      IF (MTYP*DF2 .LT.0.) GO TO 100
      E1=Z-X2
      E0=Z*E1+X1
      D=E1*E1-4*E0
      IF(D.LT.0.) GO TO 100
      Z1=(DABS(E1)+DSQRT(D))/2
      IF (E1.GT.0.) Z1=-Z1
      IF (Z.GT.Z1) Z1=E0/Z1
      IF (Z1.LT.B) GO TO 100
      IF (MTYP.EQ.0) GO TO 50
      IF (MTYP.EQ.1) GO TO 60
      IF (Z1.GT.Z) Z=Z1
      GO TO 100
60    IF (Z1.LT.Z) Z=Z1
      GO TO 100
50    F=DLOG((Z-B)/(Z1-B))+A/B/(C2-C1)*DLOG((Z+C2*B)*(Z1+C1*B)/(Z+C1*B)&
     &/(Z1+C2*B))-(Z-Z1)
      IF (F.LT.0.) Z=Z1
100   RETURN
      END
!------------------------------------------------------------------------
!
!
!	Indat reads in components and selects proper parameters from parameter
!	tables to change from 15 to k components:
!	Replace dim. 15 by K and dim. 120 by  K(K+1)/2
!
!
!	Cambié nombre por conflictos con INDAT en GC3
!		Francisco, durante "1.3.1.paraest"
!                                                                       00002780
      SUBROUTINE INDAT_cubica(NC,NONHY)

      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 T(15),P(15),OMEGA(15),CKV(15,4),O(15)
      REAL*4 ACEN(15)
      CHARACTER*4 NAME(15),NEW(15)
      DIMENSION IV(15)
      COMMON/CUB/C,C1,C2
      COMMON/CRIT/TCR(15),PCR(15)
      COMMON/STYR/N,NH,NDER,NTEMP
      COMMON/OVER/AC(15),Q(15),TSQR(15)
      COMMON/PAR/CK(120),BC(15),AC0(15),AC1(15)
      COMMON/COUT/IOUT
!
!     COMPONENTS LIST                                                   00002910
!                                                                       00002920
      DATA NAME/' C1 ',' C2 ',' C3 ',' IC4',' C4 ',' IC5',' C5 ',&
     &' C6 ',' C7 ',' C8 ',' C9 ',' H20',' N2 ',' C02',' H2S'/
!                                                                       00002950
!     CRITICAL TEMPERATURES AND PRESSURES                               00002960
!                                                                       00002970
      DATA T,P/190.6,305.4,369.8,408.1,425.2,460.4,469.6,507.4,&
     &540.2,568.8,594.6,647.3,126.2,304.2,373.2, 45.4,48.2,41.9,36.,&
     &37.5,33.4,33.3,29.3,27.,24.5,22.8,218.,33.5,72.8,88.2/
!                                                                       00003010
!     ACENTRIC FACTORS                                                  00003020
!                                                                       00003030
      DATA OMEGA/.008,.098,.152,.176,.193,.227,.251,.296,.351,&
     &.394,.440,.344,.04,.225,.1/
!                                                                       00003060
!     BINARY INTERACTION COEFFICIENTS                                   00003070
!                                                                       00003080
      DATA CKV  /2*.45,.53,2*.52,6*.5,4*0.,.02,.06,10*.08,3*0.,.12,&
     &11*.15,2*0.,.12,.08,2*.07,4*.06,.05,2*.04,2*.03,0.,.12,0./
!                                                                       00003110
!     SPECIFY COMPONENT LIST, PRINTING OPTION AND EQUATION              00003120
!                                                                       00003130
      READ(5,1) (IV(I),I=1,NC),IPRT,ICEQ
      C=ICEQ
!                                                                       00003160
!     EQUATION OF STATES USED IS:                                       00003170
!                                                                       00003180
!     P = RT/(V-B) - A/(V(V+B)+C(V-B))                                  00003190
!                                                                       00003200
!     (C=0) CORRESPONDS TO SRK. FOR PR, USE C=1                         00003210
!                                                                       00003220
!                                                                       00003230
      IF(ICEQ.EQ.0) GO TO 7
      WRITE(6,6)
      IF(IOUT.EQ.1) WRITE(3,6)
6     FORMAT(//,' PENG-ROBINSON EQUATION',//)
      GO TO 9
7     WRITE(6,8)
      IF(IOUT.EQ.1)WRITE(3,8)
8     FORMAT(//,' SOAVE-REDLICH-KWONG EQUATION',//)
9     CALL CONST(C,C1,C2,CONA,CONB)
      N = NC
1     FORMAT(16I5)
      NH=0
      DO 10 I=1,N
      L=IV(I)
      IF(L.GT.15) GO TO 2
      IF (NH.EQ.0 .AND. L.GT.11) NH=I
      NEW(I)=NAME(L)
      TCR(I)=T(L)
      PCR(I)=P(L)
      OM=OMEGA(L)
      IF(L.LT.16) GO TO 4
2     READ(5,3) NEW(I),TCR(I),PCR(I),ACEN(I)
      OM=ACEN(I)
3     FORMAT(A4,3F10.4)
4     TSQR(I)=1./DSQRT(TCR(I))
      BC(I)=CONB*TCR(I)/PCR(I)
      AC(I)=CONA*TCR(I)/DSQRT(PCR(I))
!                                                                       00003510
!     TEMPERATURE DEPENDENCE SELECTED ACCORDING TO SRK                  00003520
!                                                                       00003530
      Q(I)=.48+OM*(1.574-.176*OM)
      IF (ICEQ.EQ.1) Q(I)=.37464+OM*(1.54226-.26992*OM)
10    CONTINUE
      IF (IPRT.LT.0) GO TO 15
      WRITE(6,11)
      IF(IOUT.EQ.1) WRITE(3,11)
11    FORMAT('  COMP   TCRIT     PCRIT         M-VALUE'/)
      DO 111 I=1,N
      WRITE(6,12) NEW(I),TCR(I),PCR(I),Q(I)
111   IF(IOUT.EQ.1) WRITE(3,12) NEW(I),TCR(I),PCR(I),Q(I)
12    FORMAT(2X,A4,2F10.2,5X,F12.4)
15    NAD=0
      DO 20 I=1,N
      DO 20 J=I,N
      NAD=NAD+1
      EL=0.
      L1=IV(I)
      L2=IV(J)
      IF(L1.GT.15.OR.L2.GT.15) GO TO 20
      IF (NH.NE.0 .AND. J.GE.NH) EL=CKV(L1,L2-11)
20    CK(NAD)=EL
      IF (IPRT.LT.0 .OR. NH.EQ.0) RETURN
      WRITE(6,30) (NEW(K),K=1,N)
      IF(IOUT.EQ.1) WRITE(3,30) (NEW(K),K=1,N)
30    FORMAT(/,'  NONZERO BINARY INTERACTION COEFFICIENTS ',//,&
     &5X,15(A4,1X))
      DO 120 J=NH,N
      L1=IV(J)-11
      DO 125 I=1,N
      L=IV(I)
      IF(L1.GT.4) GO TO 125
      O(I)=CKV(L,L1)
125   CONTINUE
      IF(IOUT.EQ.1) WRITE(3,40) NEW(J),(O(I),I=1,N)
120   WRITE(6,40) NEW(J),(O(I),I=1,N)
40    FORMAT(1X,A4,15F5.2)
      RETURN
      END
!------------------------------------------------------------------------------------------
!
!
!
!
!
!     CALCULATE EQUATION CONSTANTS WHICH REPRODUCE CORRECT TC AND PC    00003930
!                                                                       00003940
      SUBROUTINE CONST(C,C1,C2,AC,BC)
      IMPLICIT REAL*8(A-H,O-Z)
      U=1+C
      W=-C
      D=U*U-4*W
      C2=(U+DSQRT(D))/2
      C1=W/C2
      S1=1+((U+W)*(U+3)-W)/2
      S2=1+U+W
      R=DSQRT(S1*S1-S2*S2*S2)
      A=(S1+R)**(1.D0/3)
      B=S2/A
      X=A+B+1
      BC=1./(3*X+U-1)
      AC=BC*DSQRT(X*(X*X-3*W)-U*W)
      RETURN
      END
!                                                                       00004120
!     CALCULATE ACTUAL VALUES (AND DERIVATIVES) OF TEMPERATURE          00004130
!     DEPENDENT PROPERTIES                                              00004140
!                                                                       00004150
      SUBROUTINE TEMSET(T)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRIT/TCR(15),PCR(15)
      COMMON/STYR/N,NH,NDER,NTEMP
      COMMON/OVER/AC(15),Q(15),TSQR(15)
      COMMON/PAR/CK(120),BC(15),AC0(15),AC1(15)
      DATA TOLD/0.D0/
      IF (T.EQ.TOLD) RETURN
      TOLD=T
      SQT=DSQRT(T)
      T2R=1.D0/T/2
      DO 10 I=1,N
      ALF=TSQR(I)
      X=1-ALF*SQT
      Q1=AC(I)*(1+Q(I))/SQT
      AC0(I)=Q1-AC(I)*Q(I)*ALF
10    AC1(I)=-Q1*T2R
      RETURN
      END
!                                                                       00004350
!     CUBGEN CALCULATES FUGACITY COEFFICIENTS AND THEIR DERIVATIVES,    00004360
!     USING MIXTURE PARAMETER ROUTINE ANEW.                             00004370
!     OTHER MODIFICATIONS ARE POSSIBLE BY REPLACING ANEW                00004380
!                                                                       00004390
      SUBROUTINE CUBGEN(MT,IC,T,P,Z,XX,FG,FT,FP,FX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XX(15),X(15),FG(15),FT(15),FP(15),FX(15,15)
      DIMENSION PD(15),AD1(15),BD1(15),ADT(15),AD2(120)
      COMMON/STYR/N,NH,NDER,NTEMP
      COMMON/CUB/C,C1,C2
      CALL TEMSET(T)
      SUM=0.
      DO 5 I=1,N
5     SUM=SUM+XX(I)
      DO 10 I=1,N
10    X(I)=XX(I)/SUM
      CALL ANEW(X,A,AT,ATT,B,AD1,AD2,ADT,BD1)
      APT=A*P/T
      BPT=B*P/T
      CALL CUBIC(MT,APT,BPT,Z)
      IC=1
      TL=1-BPT*C-3*Z
      IF (MT*TL .LT. 0.) IC=-1
      V=Z*T/P
      S1=1./(V+C1*B)
      S2=1./(V+C2*B)
      P1=P/T+A*S1*S2
      PA=-S1*S2
      PN=P1
      FAC=C1*S1+C2*S2
      P2=A*PA
      PB=P1*P1-FAC*P2
      DPDV=-P1*P1-P2*(S1+S2)
      XL1=DLOG(V/(V-B))
      FN=XL1
      XL2=DLOG(S1/S2)/(C2-C1)
      FA=-XL2/B
      F2=-A*FA
      FF=FN-F2
      GB=-V*PA
      F2B=(A*GB-F2)/B
      FB=P1-F2B
      FNB=P1
      FAB=-F2B/A
      GBB=-GB*FAC
      F2BB=(A*GBB-2*F2B)/B
      FBB=P1*P1-F2BB
      XLZ=DLOG(Z)
      FNP=FN-XLZ
      DO 15 I=1,N
15    FG(I)=FNP+FA*AD1(I)+FB*BD1(I)
      IF (NDER.LT.2 .AND. NTEMP.EQ.0) GO TO 100
      DO 20 I=1,N
20    PD(I)=(PN+PA*AD1(I)+PB*BD1(I))/DPDV
      IF (NDER.LT.2) GO TO 50
      NADR=0
      DO 30 I=1,N
      PDI=PD(I)
      CC=1+FNB*BD1(I)+PN*PDI
      CA=FAB*BD1(I)+PA*PDI
      CB=FNB+FAB*AD1(I)+FBB*BD1(I)+PB*PDI
      DO 30 J=I,N
      NADR=NADR+1
      FX(I,J)=CC+CA*AD1(J)+CB*BD1(J)+FA*AD2(NADR)
30    FX(J,I)=FX(I,J)
50    IF (NTEMP.EQ.0) GO TO 100
      DFT=FA*AT
      PTR=PA*AT
      DPDT=P/T+T*PTR
      CC=1./T
      CB=FAB*AT
      DO 70 I=1,N
      FP(I)=-PD(I)/T-1./P
70    FT(I)=CC+CB*BD1(I)+FA*ADT(I)+DPDT*PD(I)/T
100   RETURN
      END
!-----------------------------------------------------------------------------------------
!
!
!
!
!
!	Inicio de K(i) con una relación entre Tr y Pr de los puros, similar a una pv.
!
      SUBROUTINE KInit(NC,X,Y,Z,T,P)

      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 K
      DIMENSION X(15),Y(15),Z(15)
      COMMON/CRIT/TC(15),PC(15)
!
      DO 60 I=1,NC

        DIST=PC(I)/P*DEXP(5.42*(1.-TC(I)/T))
50      X(I)=Z(I)/DIST
!60	  Y(I)=Z(I)*DIST
        Y(I)=Z(I)*DIST

60    enddo
      CALL SUMA(NC,X,SX,1)
      CALL SUMA(NC,Y,SY,1)
!
      RETURN
      END
!-------------------------------------------------------------------------------------
!
!
!
!
!
      SUBROUTINE GAUSL(ND,NCOL,N,NS,A)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(ND,NCOL)
      N1=N+1
      NT=N+NS
      IF (N .EQ. 1) GO TO 50
!      START ELIMINATION
!
!
      DO 10 I=2,N
      IP=I-1
      I1=IP
      X=DABS(A(I1,I1))
      DO 11 J=I,N
      IF (DABS(A(J,I1)) .LT. X) GO TO 11
      X=DABS(A(J,I1))
      IP=J
11    CONTINUE
      IF (IP .EQ. I1) GO TO 13
!
!     ROW INTERCHANGE
!
      DO 12 J=I1,NT
      X=A(I1,J)
      A(I1,J)=A(IP,J)
12    A(IP,J)=X
13    DO 10 J=I,N
      X=A(J,I1)/A(I1,I1)
      DO 10 K=I,NT
10    A(J,K)=A(J,K) - X*A(I1,K)
!
!      ELIMINATION FINISHED, NOW BACKSUBSTITUTION
!
50    DO 20 IP=1,N
      I=N1-IP
      DO 20 K=N1,NT
      A(I,K) = A(I,K)/A(I,I)
      IF (I .EQ. 1) GO TO 20
      I1=I-1
      DO 25 J=1,I1
25    A(J,K) = A(J,K) - A(I,K)*A(J,I)
20    CONTINUE
      RETURN
      END
