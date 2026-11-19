!***********************************************************************
!
!      SOLID SOLUBILITY IN SCFs USING THE GCAEOS
!                  con efectos asociativos
!         (implementación de mecanismo generalizado)
!
!                      Susana Bottini
!                      NOVEMBER 2003

! Agregado de la asociación Alfonsina Ester Andreatta 10/08/07
!***********************************************************************
!
program GCASolub

	IMPLICIT REAL*8(A-H,O-Z)

	parameter(NCM = 15)

	DIMENSION XN(NCM),XN0(NCM),PHI(NCM),PHIV(NCM),PHILH(NCM),PHIL(NCM),DLPHI(NCM),DLPHIT(NCM),DLPHIP(NCM)
	character(2)txtiOut
!
!	COMMON/COORD/ZZ
!	COMMON/GROUP1/GSTR(15),G1(15),G2(15),TSTR(15),TSPL(15),XKIJ(15,15),AKIJ(15,15),EPX(15)
!	COMMON/GROUP2/Q(15),A(15,15),DADT(15,15),ALFA(15,15),R,NY(15,15)
!	COMMON/MOL/DC(15),D(15),DT(15),HA(15),HB(15)
!	COMMON/CRIT/TC(15),PC(15)
!	COMMON/ZZZ/HELP2(2,15),HELP4(2,15),HELP5(2,15),XNYT(15),QT,TETA(15),HELP12(2,15),HELP15(2,15),HELP7(2,15),E(2,15,15), &
!&	           PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3
!
	OPEN(UNIT=5,FILE='Soluba.DAT',STATUS='OLD',FORM='FORMATTED')
!
	READ(5,*) IOUT
!	IOUT is indicator for duplicate ouptut
	write(txtiOut,'(I2)')iOut
	if(iOut.GT.0)open(unit=iOut,file="Solub"//txtiOut//".out")

!	Input component no.,group no.
	READ(5,*) NCOMP,NGROUP,NGA
	NC=NCOMP
	NG=NGROUP
!	Input GCAEOS parameters
	CALL PARGR(NC,NG,NGA,NST, IOUT)
	CALL PARMOL(NC,NG,NST,IOUT)

!	Input solid vapor pressure	-> Aparentemente ya no es así, sólo lee DH fus y Tm
	READ(5,*) DHFU, TFU
	WRITE(IOUT,89) DHFU,TFU
89	FORMAT(/,' Parámetros del soluto:',//,'DH fusión (J/mol)=',F10.3,3X,' T fusión (K)=', F10.3,3X,/)
!
!	Input solid molar volume [cm3/mol]
	WRITE(6,6)
	IF(IOUT.NE.0) WRITE(IOUT,6)
6	FORMAT(/,' ***************************************************')
!
    	WRITE(IOUT,90)
90	FORMAT(//,2X,'T[K]',5X,'P[BAR]',5X,'Z',4X,'phi(sólido)',5X,'Solubilitdad',/)
!	Input temperature [K] and pressure [bar]
50	READ(5,*) (XN0(i),I=2,NC), T,P
!
!	Lectura de composiciones iniciales
!	READ(5,*) (XN0(i),I=2,NC)
	P=P/1.01325
	IF (T.LT.1.0E-3) GO TO 100
!	PSAT=HARA-HARB/(T)
!	PSAT=10**(PSAT)/1.01325D05
!	ARGP=(P-PSAT)*VOL/(82.054*T)
!	POY=DEXP(ARGP)
!	SID = PSAT*POY/P
!	SID = 10E-6
!	XN(1)=SID
!	do 125 i=2, NC
!125	  XN(i)=XN0(i)*(1.0-XN(1))
!
      NIT = 0
!	Variables de GCEoS
	IGZ=0
	IDER=1
	ITEMP=0
	MTYP=1

	SOLUTO=DEXP(-DHFU/8.314D0*(1D0/T-1D0/TFU))

!	DO 1 I=1,NC

!	  ARG=-2.D0*TC(I)/T/3.D0
!	  DARG=DEXP(ARG)
!1	  D(I)=DC(I)*(1.D0-.12D0*DARG)*1.065655D0


      xn(1)=1D0 - dfloat(NC-1)*1D-8
	do i=2, NC

	  XN(i) = 1D-8

	enddo
!
!	Parámetros dependientes de T:
	CALL PARAGC(T,NC,NG,NST,0)
!
!	Calculo de la fugacidad de 1 "puro":
	iTemp = 0
55	CALL GCEOS(NC,NG,NST,IDER,ITEMP,T,P,XN,PHI, DLPHI,DLPHIT,DLPHIP,CZ,IGZ,MTYP,IC)
!
!	Fugacidad del componente 1:
	PHILH(1)=DEXP(PHI(1))
!
!	Fugacidades de la mezcla:
	SID = SOLUTO
	XN(1)=SID
!	SID = 1.0E-6
!	XN(1)=SID
	do i=2, NC

	  XN(i)=XN0(i)*(1.0-XN(1))

	enddo
!
!-------------------------{Comienzo lazo GO TO}---------------------------------------
59	  CALL GCEOS(NC,NG,NST,IDER,ITEMP,T,P,XN,PHI, DLPHI,DLPHIT,DLPHIP,CZ,IGZ,MTYP,IC)
!
!	  Fugacidad del componente 1 en la mezcla
	  PHIL(1)=DEXP(PHI(1))
!	  SOLUB = SID/DEXP(PHI(1))
!
!	  xs = f(i)/fº(i)*DEXP(-DHFU/8.314D0*(1D0/T-1D0/TFU)) = gamma(i)*DEXP(-DHFU/8.314D0*(1D0/T-1D0/TFU))
        SOLUB = (PHILH(1))/(PHIL(1))*SOLUTO
!
!
	  ERROR=SOLUB*1.0D-03
!
	  DIF=SOLUB-XN(1)

	  IF (DABS(DIF).LT.ERROR) GOTO 56

	    XN(1)=SOLUB
	    do i=2, NC

	      XN(i)=XN0(i)*(1.0-XN(1))

	    enddo
	    if(NIT.EQ.0) iGz = 1	!si ya se realizó al menos una iteración, aprovecho el valor calculado de Z_liq.
	    NIT = NIT+1

	    IF (NIT.GT.50)then		!GO TO 100

	      WRITE(IOUT,*) 'No convergence after 50 iterations'
	      CLOSE(UNIT=5)
	      stop

	    endif

	GOTO 59
!---------------------------------------------------------------------------------------
56	P=1.01325*P
!	PPHI=DEXP(PHI(1))
	WRITE(IOUT,91) T,P,CZ,PHIL(1),SOLUB	!, xn(3)
91	FORMAT(F7.2,3X,F7.2,3X,F8.6,3X,E9.3,6X,E10.4, 6X,F10.4 )
	GO TO 50

100	CONTINUE
	CLOSE(UNIT=5)
	STOP

END