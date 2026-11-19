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
subroutine GCSLE (NC, NG, NST, T, P, XN0, DHfu, DVfu, Tfu, DCpfu, TtrSS, DHss, CZ, XN, phiL, iErrorSol)
!
!	I M P O R T A N T E:
!	                      El soluto es siempre el compuesto 1. Si se desea calcular toda la línea de ESL se debe
!	                      cambiar el orden de los compuestos.

	IMPLICIT REAL*8(A-H,O-Z)

	integer, parameter :: NCM = 15, maxIt = 50
	real(8), parameter :: Rgas_SI = 8.314D0, Rgas_GCA = 82.05D0

	DIMENSION XN(NCM),XN0(NCM),PHI(NCM),PHIV(NCM),PHILH(NCM),PHIL(NCM),DLPHI(NCM,NCM),DLPHIT(NCM),DLPHIP(NCM)
	common/NAME/AName(NCM)
	common/GCPROM/PMM(ncm),Pen(NCM),HHA(ncm),HHB(ncm),HHC(ncm),HHD(ncm),HHE(NCM),HHF(NCM)
!


	IF (T > 0)then	! GO TO 100
!
	  NIT = 0
!	  Variables de GCEoS
	  IGZ=0
	  IDER=1
	  ITEMP=0
	  MTYP=1

!       Pure compound change in chemical potencial:
	  TR = (Tfu/T - 1)
	  Yu = -DHFU/Rgas_SI/Tfu * TR
	  if (DVfu > 0) then
	    
	    Yu = Yu - DVfu/Rgas_GCA/T*(P - Pfu) ! Estoy mezclando unidades, pero si se respetan las unidades de la entrada
! 	                                                                    ! el resultado es correcto (adimensional).
	  endif
	  if (DCpfu > 0) then
! 	  
! 	    Take into account change in heat capacity
	    Yu = Yu + DCpfu/Rgas_SI*(TR + dlog(T/TFU))
	    DHSL = DHFU + DCpfu*(T - Tfu)
	    
	  endif 
	  if (TtrSS > 0 .AND. DHss > 0) then
	  
! 	    Take into account SS transitions
	    if (T < TtrSS) then
	    
	      Yu = Yu - DHss/Rgas_SI/TtrSS * (TtrSS/T - 1)

	    endif
	    
	  endif  	  		  
	  SOLUTO=DEXP(Yu)
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
!
	  XN(2) = XN0(2)
	  XN(1) = 1.D0 - XN(2)

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
	    IF (DABS(DIF) < ERROR)exit	!	 GOTO 56

	    XN(1)=SOLUB
! 	    XN(2:NC)=XN0(2:NC)*(1.0-XN(1))
	    XN(2) = 1D0 - xn(1)
	    
	    if(i.EQ.0) iGz = 1	!si ya se realizó al menos una iteración, aprovecho el valor calculado de Z_liq.
	    IF (i+1 > maxIt)then		!GO TO 100

	      WRITE(*,'("No se logró converger en ",I2," iteraciones")')maxIt
!	      CLOSE(UNIT=5)
	      iErrorSol = 1
	      return

	    endif
	    do io = 1, NC
	    
	      if(xn(io) <= 0)xn(io) = -xn(io)
	     
	    enddo
	    xn(:NC) = xn(:NC)/sum(xn(:NC))
! 	    write(*,*)xn(:NC)
	    
	  enddo

!	GOTO 59
!---------------------------------------------------------------------------------------
!56	    P=1.01325*P
!	    PPHI=DEXP(PHI(1))
!	    WRITE(IOUT,91) T,P,CZ,PHIL(1),SOLUB	!, xn(3)
91	  FORMAT(F7.2,3X,F7.2,3X,F8.6,3X,E9.3,6X,E10.4, 6X,F10.4 )
!	  GO TO 50

100	endif	   !CONTINUE

!  Storing activity coef. instead of fugacity coef.
	phiL(1) = phiL(1)/phiLH(1)
	return

ENDsubroutine
