!***********************************************************************
!
!      SOLID SOLUBILITY IN SCFs USING THE GCAEOS
!                  con efectos asociativos
!         (implementación de mecanismo generalizado)
!
!                      Susana Bottini
!                      NOVEMBER 2003

! Agregado de la asociación Alfonsina Ester Andreatta 10/08/07
!
! Modificado para insertar en ParaEst
!	Francisco, 03/08/13
!
!***********************************************************************
!
subroutine GCASolub(NC,NG,NST,T,P,XN0,DHfu,Tfu,CZ,XN,phiL,ErrorSol)
!
!	I M P O R T A N T E:
!	                      El soluto es siempre el compuesto 1. Si se desea calcular toda la línea de ESL se debe
!	                      cambiar el orden de los compuestos.

	IMPLICIT REAL*8(A-H,O-Z)

	parameter(NCM = 15, maxIt = 50)

	DIMENSION XN(NCM),XN0(NCM),PHI(NCM),PHIV(NCM),PHILH(NCM),PHIL(NCM),DLPHI(NCM),DLPHIT(NCM),DLPHIP(NCM)
	common/NAME/AName(NCM)
	common/GCPROM/PMM(ncm),Pen(NCM),HHA(ncm),HHB(ncm),HHC(ncm),HHD(ncm),HHE(NCM),HHF(NCM)
!
!	OPEN(UNIT=5,FILE='Soluba.DAT',STATUS='OLD',FORM='FORMATTED')
!
!	READ(5,*) IOUT
!	IOUT is indicator for duplicate ouptut
!	write(txtiOut,'(I2)')iOut
!	if(iOut.GT.0)open(unit=iOut,file="Solub"//txtiOut//".out")

!	Input component no.,group no.
!	READ(5,*) NCOMP,NGROUP,NGA
!	NC=NCOMP
!	NG=NGROUP
!	Input GCAEOS parameters
!	CALL PARGR(NC,NG,NGA,NST, IOUT)
!	CALL PARMOL(NC,NG,NST,IOUT)

!	Lee DH fus y Tm
!	READ(5,*) DHFU, TFU
!	WRITE(IOUT,89) DHFU,TFU
89	FORMAT(/,' Parámetros del soluto:',//,'DH fusión (J/mol)=',F10.3,3X,' T fusión (K)=', F10.3,3X,/)
!
!	WRITE(6,6)
	IF(IOUT.NE.0) WRITE(IOUT,6)
6	FORMAT(/,' ***************************************************')
!
!    	WRITE(IOUT,90)
90	FORMAT(//,2X,'T[K]',5X,'P[BAR]',5X,'Z',4X,'phi(sólido)',5X,'Solubilitdad',/)
!	Input temperature [K] and pressure [bar]
!50	READ(5,*) (XN0(i),I=2,NC), T,P
!
!	Lectura de composiciones iniciales
!	READ(5,*) (XN0(i),I=2,NC)
!	P=P/1.01325
	IF (T.GT.1.0E-3)then	! GO TO 100
!
	  NIT = 0
!	  Variables de GCEoS
	  IGZ=0
	  IDER=1
	  ITEMP=0
	  MTYP=1

	  SOLUTO=DEXP(-DHFU/8.314D0*(1D0/T-1D0/TFU))
	  XN(1)=1D0 - dfloat(NC-1)*1D-8
	  XN(2:NC) = 1D-8
!
!	  Parámetros dependientes de T:
	  CALL PARAGC(T,NC,NG,NST,0)
!
!	  Calculo de la fugacidad de 1 "puro":
	  iTemp = 0
55	  CALL GCEOS(NC,NG,NST,IDER,ITEMP,T,P,XN,PHI, DLPHI,DLPHIT,DLPHIP,CZ,IGZ,MTYP,IC)
!
!	  Fugacidad del componente 1:
	  PHILH(1)=DEXP(PHI(1))
!
!	  Fugacidades de la mezcla:
	  SID = SOLUTO
	  XN(1)=SID
	  XN(2:NC)=XN0(2:NC)*(1.0-XN(1))
!
!	  Comienzo lazo iterativo
!---------------------------------------------------------------------------------------------------------------------
	  do i = 1, maxIt
	    
59	    CALL GCEOS(NC,NG,NST,IDER,ITEMP,T,P,XN,PHI, DLPHI,DLPHIT,DLPHIP,CZ,IGZ,MTYP,IC)
!
!	    Fugacidad del componente 1 en la mezcla
	    PHIL(1)=DEXP(PHI(1))
!
!	    xs = f(i)/fº(i)*DEXP(-DHFU/8.314D0*(1D0/T-1D0/TFU)) = gamma(i)*DEXP(-DHFU/8.314D0*(1D0/T-1D0/TFU))
	    SOLUB = (PHILH(1))/(PHIL(1))*SOLUTO
!
	    ERROR=SOLUB*1.0D-03	!Error relativo = 10^-3
	    DIF=SOLUB-XN(1)
	    IF (DABS(DIF).LT.ERROR)exit	!	 GOTO 56

	    XN(1)=SOLUB
	    XN(2:NC)=XN0(2:NC)*(1.0-XN(1))
	    if(i.EQ.0) iGz = 1	!si ya se realizó al menos una iteración, aprovecho el valor calculado de Z_liq.
	    IF (i.GT.maxIt)then		!GO TO 100

	      WRITE(*,'("No se logró converger en ",I2," iteraciones")')maxIt
!	      CLOSE(UNIT=5)
	      ErrorSol = 1
	      return

	    endif
	    
	  enddo

!	GOTO 59
!---------------------------------------------------------------------------------------
!56	    P=1.01325*P
!	    PPHI=DEXP(PHI(1))
!	    WRITE(IOUT,91) T,P,CZ,PHIL(1),SOLUB	!, xn(3)
91	  FORMAT(F7.2,3X,F7.2,3X,F8.6,3X,E9.3,6X,E10.4, 6X,F10.4 )
!	  GO TO 50

100	endif	   !CONTINUE
!	CLOSE(UNIT=5)
	return

ENDsubroutine