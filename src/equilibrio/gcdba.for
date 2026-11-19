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
5001	FORMAT(20I3)
5003	FORMAT(36A2)
6003	FORMAT(///,1X,36A2)
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
6001	FORMAT(1H1,///,1X,'IOUT =',I2,/,'If IOUT = 0: output only on unit'
     *       ' 6',/,' If IOUT = 1: output on both unit 6 and .sal',//)
!
!	READ(5,5003) NAME
	N=NC
	NOCOMP=NC
	IOUT=ISW
!	WRITE(6,6003) NAME
!	IF(ISW.EQ.1) WRITE(3,6003) NAME
!20	IF(MODEL.GT.1) GO TO 25
20	if(model.EQ.1)then

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
6008	  FORMAT(//,' Model used for calculations: SRK-EoS',//)
!	  GO TO 29

	elseif(model.EQ.3)then

27	  WRITE(6,6009)
	  IF(IOUT.EQ.1) WRITE(3,6009)
6009	  FORMAT(//,' Model used for calculations: PR-EoS',//)
29	  CALL INDAT_cubica(N,NH)
	  ITHERM=2

      endif
30	NCALC=0
	WRITE(6,3050)
	IF(IOUT.EQ.1) WRITE(3,3050)
!
!40	READ(5,*,ERR=700,END=700) (Z(I),I=1,NC)                           00005600
!
!	Lectura en una sola línea, al estilo GC3
40	continue
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
3000	FORMAT(//,80('-'),/)

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
80	IF(INIT.NE.0)then		!EQ.0) GO TO 100

	  IF(INIT.EQ.1) CALL SUMA(N,K,SK,1)
	  DO 90 I=1,NC

	    X(I)=Z(I)
	    IF(INIT.EQ.1) Y(I)=K(I)

90	  enddo

100	endif
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
500	  WRITE(6,3600)
	  IF(IOUT.EQ.1) WRITE(3,3600)
3600	  FORMAT(//,' Calculation stopped: initial estimate out off ',
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
300	WRITE(6,3200) IT,P,T,ERROR
	IF(IOUT.EQ.1) WRITE(3,3200) IT,P,T,ERROR
3200	FORMAT(/,' ITER=',I3,2X,' Pressure(ATM)=',F10.4,2X,' Temperature='
     *,F10.4,' Error=',D12.4,/)
!
	WRITE(6,3300)
!
!	Alfonsina
!	write(8,16)	T,P ,( y(i), i=1,nc )
!16	format((F10.4),1X,(F10.4), 1x,(F10.6), 1x,(F10.6), 1x,(F10.6) )! Alfonsina
16	format((F10.4),1X,(F10.4),<NC>(X,(F10.6)))	  ! Francisco
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
3500	FORMAT(//,' Error in calcuation:  IER = ',I5,/,' IER = 1 : Trivial
     * sollution y=x',/,' IER = 10: No convergence in 10 iterations',//)
!      GO TO 45
!
!	goto 40
!
 1100 FORMAT(8F10.4)
 1200 FORMAT(I3,5F10.4)

 3050 FORMAT(//,' ITYPE=1 : Bubble point temperature',/,' ITYPE=2 : DEW
     * POINT TEMPERATURE',/,' ITYPE=3: Bubble point pressure ',/,' ITYPE
     *=4: Dew point pressure',/)
 3100 FORMAT(/,' Calculation number: ',I5,' , ITYPE= ',I5,//)

 3300 FORMAT(/,' Comp.      Feed    Incipient phase     ln(K)',/)
 3400 FORMAT(I3,7X,F10.7,2X,F10.7,6X,F10.6)                             00006200
!
  700 CLOSE (UNIT=5)
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
      SUBROUTINE DBTP(NC,NC1,NG,NST,T,P,X,Y,Z,DG,IER)                   00000020
C     ******************************************************************00000030
C                                                                       00000040
C     ERROR MESSAGES (IER VALUES)                                       00000050
C     0     : NO ERRORS                                                 00000060
C     1     : THE TRIVIAL SOLUTION  X=Y                                 00000070
C     10    : THE EQUILIBRIUM CALCULATIONS DID NOT CONVERGE IN MAXIT    00000080
C             ITERATIONS                                                00000090
C     11    : ERROR ENCOUNTERED DURING CALCULATION OF DYDX AND DPDX.    00000100
C             THE EQUILIBRIUM CALCULATION IS ALL RIGHT.                 00000110
C                                                                       00000120
C     99,100                                                            00000130
C     199,200: SEE ZMAX                                                 00000140
C                                                                       00000150
C     ******************************************************************00000160
      IMPLICIT REAL*8(A-H,O-Z)                                          00000170
      DIMENSION X(15),Y(15),Z(15),G(15),DG(NC,NC1),DEL(15),PHIL(15)     00000180
      DIMENSION PHIV(15),FX(15,15),FY(15,15),FPX(15),FPY(15),FTX(15)    00000190
      DIMENSION FTY(15)                                                 00000200
      COMMON/CONTL/ERROR,IPRINT,NDIG,INIT,ITYPE,MAXIT,ITHERM,NIT        00000210
      COMMON/COUT/IOUT

	COMMON/Fact_comp/ZL,ZV
	real(8), parameter :: macheps = epsilon(1D0)

      EPS=10.D0**(-NDIG)                                                00000230
      DXYMAX=0.75D0                                                     00000240
      DPMAX=0.75D0                                                      00000250
      DTMAX=0.75D0                                                      00000260
      NIT=0                                                             00000270
      IGZ=0                                                             00000280
      IER=0                                                             00000290
      IGOLD=IGZ                                                         00000300
      NCC=NC-1                                                          00000310
      IXY=1
      IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) IXY=2                                00000330
    3 GO TO(100,200,100,200),ITYPE                                      00000340
  100 IF(ITHERM.EQ.1.AND.ITYPE.LE.2) CALL PARAGC(T,NC,NG,NST,1)         00000350
      IF(ITHERM.EQ.1)CALL GCEOS(NC,NG,NST,2,1,T,P,X,PHIL,FX,FTX,FPX,
     *ZL,IGZ,1,ICL)                                                     00000370
      IF(ITHERM.EQ.1) CALL GCEOS(NC,NG,NST,2,1,T,P,Y,PHIV,FY,FTY,FPY,
     *ZV,IGZ,-1,ICV)                                                    00000390
      IF(ITHERM.EQ.2) CALL SRK1(2,1,ICL,T,P,ZL,X,PHIL,FTX,FPX,FX)       00000400
      IF(ITHERM.EQ.2) CALL SRK1(2,-1,ICV,T,P,ZV,Y,PHIV,FTY,FPY,FY)      00000410
      GO TO 300                                                         00000420
  200 IF(ITHERM.EQ.1.AND.ITYPE.LE.2) CALL PARAGC(T,NC,NG,NST,1)         00000430
      IF(ITHERM.EQ.1) CALL GCEOS(NC,NG,NST,2,1,T,P,Y,PHIL,FX,FTX,FPX,
     *ZL,IGZ,1,ICL)                                                     00000450
      IF(ITHERM.EQ.1) CALL GCEOS(NC,NG,NST,2,1,T,P,X,PHIV,FY,FTY,FPY,
     *ZV,IGZ,-1,ICV)                                                    00000470
      IF(ITHERM.EQ.2) CALL SRK1(2,1,ICL,T,P,ZL,Y,PHIL,FTX,FPX,FX)       00000480
      IF(ITHERM.EQ.2) CALL SRK1(2,-1,ICV,T,P,ZV,X,PHIV,FTY,FPY,FY)      00000490
  300 IF(ICL.EQ.1.AND.ICV.EQ.1) GOTO 400                                00000500
      IER=99                                                            00000510
      GO TO 99                                                          00000520
  400 ZCRT=DABS(ZV-ZL)                                                  00000530
      IF(ZCRT.LT.1.D-6) IER=1                                           00000540
      IF (IER.GT.0) GOTO 2                                              00000550
      
      GNORM=0.D0                                                        00000560
      DO 4 I=1,NC  
        xi = X(I)+macheps
        yi = Y(I)+macheps
        IF(IXY.EQ.1) G(I)=DLOG(xi/yi)-PHIV(I)+PHIL(I)        
        IF(IXY.EQ.2) G(I)=DLOG(yi/xi)-PHIV(I)+PHIL(I)             
    4   GNORM=GNORM+G(I)**2                                             
    
      GNORM=DSQRT(GNORM)                                                00000610
      ERROR=GNORM                                                       00000620
      IF(IPRINT.EQ.0) GO TO 40                                          00000630
      WRITE(6,1000) NIT,P,T,ERROR                                       00000640
      IF(IOUT.EQ.1) WRITE(3,1000) NIT,P,T,ERROR                         00000650
   40 DO 6 I=1,NC                                                       00000660
      IF(ITYPE.GT.2) DG(I,NC)=-FPY(I)+FPX(I)                            00000670
      IF(ITYPE.LE.2) DG(I,NC)=-FTY(I)+FTX(I)                            00000680
      DO 5 J=1,NCC                                                      00000690
      IF(IXY.EQ.1) DG(I,J)=-FY(I,J)+FY(I,NC)                            00000700
      IF(IXY.EQ.2) DG(I,J)=FX(I,J)-FX(I,NC)                             00000710
      IF(I.EQ.J.AND.IXY.EQ.1) DG(I,J)=DG(I,J)-1.D0/Y(I)                 00000720
      IF(I.EQ.J.AND.IXY.EQ.2) DG(I,J)=DG(I,J)+1.D0/Y(I)                 00000730
      IF(I.EQ.NC.AND.IXY.EQ.1) DG(I,J)=DG(I,J)+1.D0/Y(I)                00000740
    5 IF(I.EQ.NC.AND.IXY.EQ.2) DG(I,J)=DG(I,J)-1.D0/Y(I)                00000750
    6 DG(I,NC1)=-G(I)                                                   00000760
      IF(GNORM.LT.EPS) GOTO 2                                           00000770
      CALL GAUSL(NC,NC1,NC,1,DG)                                        00000780
      RATMIN=1.D0                                                       00000790
      DO 7 I=1,NCC                                                      00000800
      CRITY=DXYMAX*Y(I)                                                 00000810
      IF(Y(I).GT..5D0) CRITY=.5D0*(1.D0-Y(I))                           00000820
      DYABS=DABS(DG(I,NC1))                                             00000830
      IF(DYABS.LT.CRITY) GOTO 7                                         00000840
      RAT=CRITY/DYABS                                                   00000850
      IF(RAT.LT.RATMIN) RATMIN=RAT                                      00000860
    7 CONTINUE                                                          00000870
      PCT=DABS(DG(NC,NC1))                                              00000880
      IF(ITYPE.GT.2) PCRIT=DPMAX*P                                      00000890
      IF(ITYPE.LE.2) PCRIT=DTMAX*T                                      00000900
      RAT=1.D0                                                          00000910
      IF (PCT.GT.PCRIT) RAT=PCRIT/PCT                                   00000920
      IF(RAT.LT.RATMIN) RATMIN=RAT                                      00000930
      YSUM=0.D0                                                         00000940
      DO 9 I=1,NCC                                                      00000950
    9 YSUM=YSUM+Y(I)+DG(I,NC1)*RATMIN                                   00000960
      CRITY=.5D0*Y(NC)                                                  00000970
      YTEST=1.D0-YSUM-Y(NC)                                             00000980
      YTEST=DABS(YTEST)                                                 00000990
      IF(YTEST.GT.CRITY) RAT=CRITY/YTEST                                00001000
      IF(RAT.LT.RATMIN) RATMIN=RAT                                      00001010
      IF(ITYPE.GT.2) P=P+DG(NC,NC1)*RATMIN                              00001020
      IF(ITYPE.LE.2) T=T+DG(NC,NC1)*RATMIN                              00001030
      YSUM=0.D0                                                         00001040
      DO 8 I=1,NCC                                                      00001050
      Y(I)=Y(I)+DG(I,NC1)*RATMIN                                        00001060
    8 YSUM=YSUM+Y(I)                                                    00001070
      Y(NC)=1.D0-YSUM                                                   00001080
      NIT=NIT+1                                                         00001090
      IF(NIT.EQ.5) IGZ=IGZ+2                                            00001100
      IF(NIT.GT.MAXIT) IER=10                                           00001110
      IF(IER.NE.0) GOTO 2                                               00001120
      GOTO 3                                                            00001130
    2 DO 98 I=1,NC                                                      00001140
      Z(I)=PHIL(I)-PHIV(I)                                              00001150
   98 CONTINUE                                                          00001160
   99 IGZ=IGOLD                                                         00001170
 1000 FORMAT(' ITER=',I5,' VARIABLE1=',F10.4,' VARIABLE2=',F10.4,' RESID00001180
     *UE=',D14.4)
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
C     SRK1 IS THE INTERFACE ROUTINE TO CUBGEN                           00001350
C                                                                       00001360
      SUBROUTINE SRK1(INDIC,MTYP,IC,T,P,Z,X,FUG,FUGT,FUGP,FUGX)         00001370
      IMPLICIT REAL*8 (A-H,O-Z)                                         00001380
      DIMENSION X(15),FUGT(15),FUG(15),FUGP(15),FUGX(15,15)             00001390
      COMMON/STYR/ N,NH,NDER,NTEMP                                      00001400
      NTEMP=0                                                           00001410
      NDER=1                                                            00001420
      IF (INDIC.GT.1) NDER=2                                            00001430
      IF (INDIC.EQ.2 .OR. INDIC.EQ.4) NTEMP=1                           00001440
      CALL CUBGEN(MTYP,IC,T,P,Z,X,FUG,FUGT,FUGP,FUGX)                   00001450
      RETURN                                                            00001460
      END                                                               00001470
C                                                                       00001480
C     SUBROUTINE ANEW CALCULATES PARAMETERS A AND B AND THEIR           00001490
C     PARTIAL COMPOSITION AND TEMPERATURE DERIVATIVES                   00001500
C                                                                       00001510
      SUBROUTINE ANEW(X,A0,AT,ATT,B0,AD1,AD2,ADT,BD1)                   00001520
      IMPLICIT REAL*8(A-H,O-Z)                                          00001530
      DIMENSION AD1(15),AD2(120),ADT(15),BD1(15),                       00001540
     *X(15),AX0(15),AX1(15),AS1(15),AST1(15)                            00001550
      COMMON/STYR/N,NH,NDER,NTEMP                                       00001560
      COMMON/PAR/CK(120),BC(15),AC0(15),AC1(15)                         00001570
      S0=0.                                                             00001580
      S1=0.                                                             00001590
      DO 10 I=1,N                                                       00001600
      XI=X(I)                                                           00001610
      AS1(I)=0.                                                         00001620
      AX0(I)=XI*AC0(I)                                                  00001630
      IF (NTEMP.EQ.0) GO TO 10                                          00001640
      AST1(I)=0.                                                        00001650
      AX1(I)=XI*AC1(I)                                                  00001660
      S1=S1+AX1(I)                                                      00001670
  10  S0=S0+AX0(I)                                                      00001680
      A0=S0*S0                                                          00001690
      NADR=0                                                            00001700
   15 DO 20 I=1,N                                                       00001710
      FAC=2*AC0(I)                                                      00001720
      AD1(I)=S0*FAC                                                     00001730
      IF (NTEMP.NE.0) ADT(I)=FAC*S1+2*S0*AC1(I)                         00001740
      IF (NDER.LT.2) GO TO 20                                           00001750
      DO 25 J=I,N                                                       00001760
      NADR=NADR+1                                                       00001770
   25 AD2(NADR)=FAC*AC0(J)                                              00001780
   20 CONTINUE                                                          00001790
      NADR=0                                                            00001800
C                                                                       00001810
C     START MAIN LOOP FOR A-PARAMETER                                   00001820
C                                                                       00001830
      DO 50 I=1,N                                                       00001840
      IF(I.EQ.N) GO TO 50                                               00001850
      K1=I+1                                                            00001860
      NADR=NADR+1                                                       00001870
      XI=X(I)                                                           00001880
      AX0I=AX0(I)                                                       00001890
      IF (NTEMP.NE.0.) AX1I=AX1(I)                                      00001900
      FAC=2*AC0(I)                                                      00001910
      DO 60 K=K1,N                                                      00001920
C                                                                       00001930
C     START INNER LOOP                                                  00001940
C                                                                       00001950
      NADR=NADR+1                                                       00001960
      EF=CK(NADR)                                                       00001970
      IF(EF.EQ.0.) GO TO 60                                             00001980
      AX0K=AX0(K)                                                       00001990
      AS1(I)=AS1(I)+AX0K*EF                                             00002000
      AS1(K)=AS1(K)+AX0I*EF                                             00002010
      IF (NTEMP.EQ.0) GO TO 65                                          00002020
      AST1(I)=AST1(I)+AX1(K)*EF                                         00002030
      AST1(K)=AST1(K)+AX1I*EF                                           00002040
   65 IF (NDER.LT.2) GO TO 60                                           00002050
      FAC1=AC0(K)*FAC*EF                                                00002060
      AD2(NADR)=AD2(NADR)-FAC1                                          00002070
   60 CONTINUE                                                          00002080
   50 CONTINUE                                                          00002090
      A0=0                                                              00002100
      AT=0.                                                             00002110
      DO 70 I=1,N                                                       00002120
      FAC=2*AC0(I)                                                      00002130
      AD1(I)=AD1(I)-FAC*AS1(I)                                          00002140
      A0=A0+AD1(I)*X(I)/2                                               00002150
      IF (NTEMP.EQ.0) GO TO 70                                          00002160
      ADT(I)=ADT(I)-2*AC1(I)*AS1(I)-FAC*AST1(I)                         00002170
      AT=AT+X(I)*ADT(I)                                                 00002180
   70 CONTINUE                                                          00002190
      AT=AT/2                                                           00002200
C                                                                       00002210
C     START B-CALCULATION                                               00002220
C                                                                       00002230
      B0=0.                                                             00002240
      DO 120 I=1,N                                                      00002250
      B0=B0+X(I)*BC(I)                                                  00002260
      BD1(I)=BC(I)                                                      00002270
  120 CONTINUE                                                          00002280
      RETURN                                                            00002290
      END                                                               00002300
C                                                                       00002310
C     CUBIC SOLVES FOR THE COMPRESSIBILITY FACTOR                       00002320
C     MTYP=1: LIQUID; MTYP=-1: VAPOUR; MTYP= 0: MIN. GIBBS ENERGY       00002330
C                                                                       00002340
      SUBROUTINE CUBIC(MTYP,A,B,Z)                                      00002350
      IMPLICIT REAL*8(A-H,O-Z)                                          00002360
      COMMON/CUB/C,C1,C2                                                00002370
      BC=B*C                                                            00002380
      X2=1-BC                                                           00002390
      X1=A-B*(1+B+C+2*BC)                                               00002400
      X0=B*(A-BC*(1+B))                                                 00002410
      Z=X2/3                                                            00002420
      F=((Z-X2)*Z+X1)*Z-X0                                              00002430
      Z=B                                                               00002440
      IF(F.LT.0..OR. B.GT. X2/3)Z=Z+1                                   00002450
   10 DF2=3*Z-X2                                                        00002460
      DF1=Z*(DF2-X2)+X1                                                 00002470
      F=((Z-X2)*Z+X1)*Z-X0                                              00002480
      DZ=F/DF1                                                          00002490
      DZ=DZ*(1+DZ*DF2/DF1)                                              00002500
      Z=Z-DZ                                                            00002510
      IF (DABS(DZ).GT. 1.D-10) GO TO 10                                 00002520
      IF (MTYP*DF2 .LT.0.) GO TO 100                                    00002530
      E1=Z-X2                                                           00002540
      E0=Z*E1+X1                                                        00002550
      D=E1*E1-4*E0                                                      00002560
      IF(D.LT.0.) GO TO 100                                             00002570
      Z1=(DABS(E1)+DSQRT(D))/2                                          00002580
      IF (E1.GT.0.) Z1=-Z1                                              00002590
      IF (Z.GT.Z1) Z1=E0/Z1                                             00002600
      IF (Z1.LT.B) GO TO 100                                            00002610
      IF (MTYP.EQ.0) GO TO 50                                           00002620
      IF (MTYP.EQ.1) GO TO 60                                           00002630
      IF (Z1.GT.Z) Z=Z1                                                 00002640
      GO TO 100                                                         00002650
   60 IF (Z1.LT.Z) Z=Z1                                                 00002660
      GO TO 100                                                         00002670
   50 F=DLOG((Z-B)/(Z1-B))+A/B/(C2-C1)*DLOG((Z+C2*B)*(Z1+C1*B)/(Z+C1*B) 00002680
     1/(Z1+C2*B))-(Z-Z1)                                                00002690
      IF (F.LT.0.) Z=Z1                                                 00002700
  100 RETURN                                                            00002710
      END                                                               00002720
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
C     COMPONENTS LIST                                                   00002910
C                                                                       00002920
      DATA NAME/' C1 ',' C2 ',' C3 ',' IC4',' C4 ',' IC5',' C5 ',       00002930
     *' C6 ',' C7 ',' C8 ',' C9 ',' H20',' N2 ',' C02',' H2S'/          00002940
C                                                                       00002950
C     CRITICAL TEMPERATURES AND PRESSURES                               00002960
C                                                                       00002970
      DATA T,P/190.6,305.4,369.8,408.1,425.2,460.4,469.6,507.4,         00002980
     *540.2,568.8,594.6,647.3,126.2,304.2,373.2, 45.4,48.2,41.9,36.,    00002990
     *37.5,33.4,33.3,29.3,27.,24.5,22.8,218.,33.5,72.8,88.2/            00003000
C                                                                       00003010
C     ACENTRIC FACTORS                                                  00003020
C                                                                       00003030
      DATA OMEGA/.008,.098,.152,.176,.193,.227,.251,.296,.351,          00003040
     *.394,.440,.344,.04,.225,.1/                                       00003050
C                                                                       00003060
C     BINARY INTERACTION COEFFICIENTS                                   00003070
C                                                                       00003080
      DATA CKV  /2*.45,.53,2*.52,6*.5,4*0.,.02,.06,10*.08,3*0.,.12,     00003090
     *11*.15,2*0.,.12,.08,2*.07,4*.06,.05,2*.04,2*.03,0.,.12,0./        00003100
C                                                                       00003110
C     SPECIFY COMPONENT LIST, PRINTING OPTION AND EQUATION              00003120
C                                                                       00003130
      READ(5,1) (IV(I),I=1,NC),IPRT,ICEQ                                00003140
      C=ICEQ                                                            00003150
C                                                                       00003160
C     EQUATION OF STATES USED IS:                                       00003170
C                                                                       00003180
C     P = RT/(V-B) - A/(V(V+B)+C(V-B))                                  00003190
C                                                                       00003200
C     (C=0) CORRESPONDS TO SRK. FOR PR, USE C=1                         00003210
C                                                                       00003220
C                                                                       00003230
      IF(ICEQ.EQ.0) GO TO 7                                             00003240
      WRITE(6,6)                                                        00003250
      IF(IOUT.EQ.1) WRITE(3,6)                                          00003260
    6 FORMAT(//,' PENG-ROBINSON EQUATION',//)                           00003270
      GO TO 9                                                           00003280
    7 WRITE(6,8)                                                        00003290
      IF(IOUT.EQ.1)WRITE(3,8)                                           00003300
    8 FORMAT(//,' SOAVE-REDLICH-KWONG EQUATION',//)                     00003310
    9 CALL CONST(C,C1,C2,CONA,CONB)                                     00003320
      N = NC                                                            00003330
    1 FORMAT(16I5)                                                      00003340
      NH=0                                                              00003350
      DO 10 I=1,N                                                       00003360
      L=IV(I)                                                           00003370
      IF(L.GT.15) GO TO 2                                               00003380
      IF (NH.EQ.0 .AND. L.GT.11) NH=I                                   00003390
      NEW(I)=NAME(L)                                                    00003400
      TCR(I)=T(L)                                                       00003410
      PCR(I)=P(L)                                                       00003420
      OM=OMEGA(L)                                                       00003430
      IF(L.LT.16) GO TO 4                                               00003440
    2 READ(5,3) NEW(I),TCR(I),PCR(I),ACEN(I)                            00003450
      OM=ACEN(I)                                                        00003460
    3 FORMAT(A4,3F10.4)                                                 00003470
    4 TSQR(I)=1./DSQRT(TCR(I))                                          00003480
      BC(I)=CONB*TCR(I)/PCR(I)                                          00003490
      AC(I)=CONA*TCR(I)/DSQRT(PCR(I))                                   00003500
C                                                                       00003510
C     TEMPERATURE DEPENDENCE SELECTED ACCORDING TO SRK                  00003520
C                                                                       00003530
      Q(I)=.48+OM*(1.574-.176*OM)                                       00003540
      IF (ICEQ.EQ.1) Q(I)=.37464+OM*(1.54226-.26992*OM)                 00003550
   10 CONTINUE                                                          00003560
      IF (IPRT.LT.0) GO TO 15                                           00003570
      WRITE(6,11)                                                       00003580
      IF(IOUT.EQ.1) WRITE(3,11)                                         00003590
   11 FORMAT('  COMP   TCRIT     PCRIT         M-VALUE'/)               00003600
      DO 111 I=1,N                                                      00003610
      WRITE(6,12) NEW(I),TCR(I),PCR(I),Q(I)                             00003620
111   IF(IOUT.EQ.1) WRITE(3,12) NEW(I),TCR(I),PCR(I),Q(I)               00003630
   12 FORMAT(2X,A4,2F10.2,5X,F12.4)                                     00003640
   15 NAD=0                                                             00003650
      DO 20 I=1,N                                                       00003660
      DO 20 J=I,N                                                       00003670
      NAD=NAD+1                                                         00003680
      EL=0.                                                             00003690
      L1=IV(I)                                                          00003700
      L2=IV(J)                                                          00003710
      IF(L1.GT.15.OR.L2.GT.15) GO TO 20                                 00003720
      IF (NH.NE.0 .AND. J.GE.NH) EL=CKV(L1,L2-11)                       00003730
   20 CK(NAD)=EL                                                        00003740
      IF (IPRT.LT.0 .OR. NH.EQ.0) RETURN                                00003750
      WRITE(6,30) (NEW(K),K=1,N)                                        00003760
      IF(IOUT.EQ.1) WRITE(3,30) (NEW(K),K=1,N)                          00003770
   30 FORMAT(/,'  NONZERO BINARY INTERACTION COEFFICIENTS ',//,         00003780
     *5X,15(A4,1X))                                                     00003790
      DO 120 J=NH,N                                                     00003800
      L1=IV(J)-11                                                       00003810
      DO 125 I=1,N                                                      00003820
      L=IV(I)                                                           00003830
      IF(L1.GT.4) GO TO 125                                             00003840
      O(I)=CKV(L,L1)                                                    00003850
  125 CONTINUE                                                          00003860
      IF(IOUT.EQ.1) WRITE(3,40) NEW(J),(O(I),I=1,N)                     00003870
  120 WRITE(6,40) NEW(J),(O(I),I=1,N)                                   00003880
   40 FORMAT(1X,A4,15F5.2)                                              00003890
      RETURN                                                            00003900
      END                                                               00003910
!------------------------------------------------------------------------------------------
!
!
!
!
!
C     CALCULATE EQUATION CONSTANTS WHICH REPRODUCE CORRECT TC AND PC    00003930
C                                                                       00003940
      SUBROUTINE CONST(C,C1,C2,AC,BC)                                   00003950
      IMPLICIT REAL*8(A-H,O-Z)                                          00003960
      U=1+C                                                             00003970
      W=-C                                                              00003980
      D=U*U-4*W                                                         00003990
      C2=(U+DSQRT(D))/2                                                 00004000
      C1=W/C2                                                           00004010
      S1=1+((U+W)*(U+3)-W)/2                                            00004020
      S2=1+U+W                                                          00004030
      R=DSQRT(S1*S1-S2*S2*S2)                                           00004040
      A=(S1+R)**(1.D0/3)                                                00004050
      B=S2/A                                                            00004060
      X=A+B+1                                                           00004070
      BC=1./(3*X+U-1)                                                   00004080
      AC=BC*DSQRT(X*(X*X-3*W)-U*W)                                      00004090
      RETURN                                                            00004100
      END                                                               00004110
C                                                                       00004120
C     CALCULATE ACTUAL VALUES (AND DERIVATIVES) OF TEMPERATURE          00004130
C     DEPENDENT PROPERTIES                                              00004140
C                                                                       00004150
      SUBROUTINE TEMSET(T)                                              00004160
      IMPLICIT REAL*8(A-H,O-Z)                                          00004170
      COMMON/CRIT/TCR(15),PCR(15)                                       00004180
      COMMON/STYR/N,NH,NDER,NTEMP                                       00004190
      COMMON/OVER/AC(15),Q(15),TSQR(15)                                 00004200
      COMMON/PAR/CK(120),BC(15),AC0(15),AC1(15)                         00004210
      DATA TOLD/0.D0/                                                   00004220
      IF (T.EQ.TOLD) RETURN                                             00004230
      TOLD=T                                                            00004240
      SQT=DSQRT(T)                                                      00004250
      T2R=1.D0/T/2                                                      00004260
      DO 10 I=1,N                                                       00004270
      ALF=TSQR(I)                                                       00004280
      X=1-ALF*SQT                                                       00004290
      Q1=AC(I)*(1+Q(I))/SQT                                             00004300
      AC0(I)=Q1-AC(I)*Q(I)*ALF                                          00004310
   10 AC1(I)=-Q1*T2R                                                    00004320
      RETURN                                                            00004330
      END                                                               00004340
C                                                                       00004350
C     CUBGEN CALCULATES FUGACITY COEFFICIENTS AND THEIR DERIVATIVES,    00004360
C     USING MIXTURE PARAMETER ROUTINE ANEW.                             00004370
C     OTHER MODIFICATIONS ARE POSSIBLE BY REPLACING ANEW                00004380
C                                                                       00004390
      SUBROUTINE CUBGEN(MT,IC,T,P,Z,XX,FG,FT,FP,FX)                     00004400
      IMPLICIT REAL*8(A-H,O-Z)                                          00004410
      DIMENSION XX(15),X(15),FG(15),FT(15),FP(15),FX(15,15)             00004420
      DIMENSION PD(15),AD1(15),BD1(15),ADT(15),AD2(120)                 00004430
      COMMON/STYR/N,NH,NDER,NTEMP                                       00004440
      COMMON/CUB/C,C1,C2                                                00004450
      CALL TEMSET(T)                                                    00004460
      SUM=0.                                                            00004470
      DO 5 I=1,N                                                        00004480
    5 SUM=SUM+XX(I)                                                     00004490
      DO 10 I=1,N                                                       00004500
   10 X(I)=XX(I)/SUM                                                    00004510
      CALL ANEW(X,A,AT,ATT,B,AD1,AD2,ADT,BD1)                           00004520
      APT=A*P/T                                                         00004530
      BPT=B*P/T                                                         00004540
      CALL CUBIC(MT,APT,BPT,Z)                                          00004550
      IC=1                                                              00004560
      TL=1-BPT*C-3*Z                                                    00004570
      IF (MT*TL .LT. 0.) IC=-1                                          00004580
      V=Z*T/P                                                           00004590
      S1=1./(V+C1*B)                                                    00004600
      S2=1./(V+C2*B)                                                    00004610
      P1=P/T+A*S1*S2                                                    00004620
      PA=-S1*S2                                                         00004630
      PN=P1                                                             00004640
      FAC=C1*S1+C2*S2                                                   00004650
      P2=A*PA                                                           00004660
      PB=P1*P1-FAC*P2                                                   00004670
      DPDV=-P1*P1-P2*(S1+S2)                                            00004680
      XL1=DLOG(V/(V-B))                                                 00004690
      FN=XL1                                                            00004700
      XL2=DLOG(S1/S2)/(C2-C1)                                           00004710
      FA=-XL2/B                                                         00004720
      F2=-A*FA                                                          00004730
      FF=FN-F2                                                          00004740
      GB=-V*PA                                                          00004750
      F2B=(A*GB-F2)/B                                                   00004760
      FB=P1-F2B                                                         00004770
      FNB=P1                                                            00004780
      FAB=-F2B/A                                                        00004790
      GBB=-GB*FAC                                                       00004800
      F2BB=(A*GBB-2*F2B)/B                                              00004810
      FBB=P1*P1-F2BB                                                    00004820
      XLZ=DLOG(Z)                                                       00004830
      FNP=FN-XLZ                                                        00004840
      DO 15 I=1,N                                                       00004850
   15 FG(I)=FNP+FA*AD1(I)+FB*BD1(I)                                     00004860
      IF (NDER.LT.2 .AND. NTEMP.EQ.0) GO TO 100                         00004870
      DO 20 I=1,N                                                       00004880
   20 PD(I)=(PN+PA*AD1(I)+PB*BD1(I))/DPDV                               00004890
      IF (NDER.LT.2) GO TO 50                                           00004900
      NADR=0                                                            00004910
      DO 30 I=1,N                                                       00004920
      PDI=PD(I)                                                         00004930
      CC=1+FNB*BD1(I)+PN*PDI                                            00004940
      CA=FAB*BD1(I)+PA*PDI                                              00004950
      CB=FNB+FAB*AD1(I)+FBB*BD1(I)+PB*PDI                               00004960
      DO 30 J=I,N                                                       00004970
      NADR=NADR+1                                                       00004980
      FX(I,J)=CC+CA*AD1(J)+CB*BD1(J)+FA*AD2(NADR)                       00004990
   30 FX(J,I)=FX(I,J)                                                   00005000
   50 IF (NTEMP.EQ.0) GO TO 100                                         00005010
      DFT=FA*AT                                                         00005020
      PTR=PA*AT                                                         00005030
      DPDT=P/T+T*PTR                                                    00005040
      CC=1./T                                                           00005050
      CB=FAB*AT                                                         00005060
      DO 70 I=1,N                                                       00005070
      FP(I)=-PD(I)/T-1./P                                               00005080
   70 FT(I)=CC+CB*BD1(I)+FA*ADT(I)+DPDT*PD(I)/T                         00005090
  100 RETURN                                                            00005100
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
50	  X(I)=Z(I)/DIST
!60	  Y(I)=Z(I)*DIST
	  Y(I)=Z(I)*DIST

60	enddo
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
C      START ELIMINATION
C
C
      DO 10 I=2,N
      IP=I-1
      I1=IP
      X=DABS(A(I1,I1))
      DO 11 J=I,N
      IF (DABS(A(J,I1)) .LT. X) GO TO 11
      X=DABS(A(J,I1))
      IP=J
   11 CONTINUE
      IF (IP .EQ. I1) GO TO 13
C
C     ROW INTERCHANGE
C
      DO 12 J=I1,NT
      X=A(I1,J)
      A(I1,J)=A(IP,J)
   12 A(IP,J)=X
   13 DO 10 J=I,N
      X=A(J,I1)/A(I1,I1)
      DO 10 K=I,NT
   10 A(J,K)=A(J,K) - X*A(I1,K)
C
C      ELIMINATION FINISHED, NOW BACKSUBSTITUTION
C
   50 DO 20 IP=1,N
      I=N1-IP
      DO 20 K=N1,NT
      A(I,K) = A(I,K)/A(I,I)
      IF (I .EQ. 1) GO TO 20
      I1=I-1
      DO 25 J=1,I1
   25 A(J,K) = A(J,K) - A(I,K)*A(J,I)
   20 CONTINUE
      RETURN
      END
