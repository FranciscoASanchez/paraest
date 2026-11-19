! Cálculo de valores de equilibrio trifásico SLV con la GCA
! Se supone fase sólida como compuesto puro, y se necesita el
! valor exacto de la presión de vapor predicha por el modelo
! en el punto triple de sólido.
! 
! A diferencia de los otros cálculos, no se evalúa la desvia-
! ción punto por punto sino que se calcula la línea SLV pre-
! dicha hasta Pmáx
! 
subroutine SLV(NData)

	implicit none
	integer, parameter ::             NCM=15, NGM=15, NGAM=24, iTemp = 1
	real(8), parameter ::             delta = 1.D-5
	                                  
	integer ::                        NData, NC, NG, NST, iErrPC
	                                  
	real(8) ::                        T, P, Z, dlnF_dz1, d2lnF_dz12, dlnFper1, dlnFper2
	real(8), dimension(NData) ::      T, P, x, y
	real(8), dimension(NCM) ::        zn, aux, aux2	
	

	
	
! Cálculo de cada punto	

	
! Salida si criterio de convergencia
	
! Cálculo de nueva variable especificada:
	
	return
endsubroutine SLV
! -------------------------------------------------------------------------
! 
!  Cálculo de la función a resolver por punto y el jacoviano de dicha función
! 
subroutine FyJacobianoSLV(NC, NG, NST, T, P, x, Z, o, spec)

	implicit none
	integer, parameter ::             NCM=15, NGM=15, NGAM=24, iTemp = 1, L = 2, V = 1, S = 3
	                                  
	integer ::                        i, NC, NG, NST
	                                  
	real(8) ::                        o(NC+2), P, spec, T
	real(8), dimension(3,NCM) ::      x, lnX, lnPhi, dLnPhi_T, dLnPhi_P, Z
	real(8), dimension(NCM,NCM) ::    dLnPhi_n
	                                  
	real(8), dimension(NC+2) ::       F, o, W
	real(8), dimension(NC+2,NC+2)::   Jac

	x(S,2) = 1.D0
	x(S,1) = 0.D0

	
	W(1) = dlog(T)
	W(2) = dlog(P)
	lnX(:3,:NC) = dlog( x(:3,:NC) )
	W(2+L) = lnX(L,1)
	W(2+V) = lnX(V,2)
	
	call fugacitySLV()
	
	do i = 1, NC
	
	  F(i) = lnPhi(L,i) + lnX(L,i) - lnPhi(V,i) - lnX(V,i)
!  dFi/d(ln T)	  
	  Jac(i,1) = T*(dLnPhi_T(L,i) - dLnPhi_T(V,i))
!  dFi/d(ln P) =	  
	  Jac(i,2) = P*(dLnPhi_P(L,i) - dLnPhi_P(V,i))
!  dFi/d(ln x1) (con n1 + n2 = 1)
	  Jac(i,3) = x(L,1)*(dLnPhi_n(L,i,1) - 1.D0/x(L,1) - dLnPhi_n(L,i,2))
!  dFi/d(ln y2)
	  Jac(i,4) = x(V,1)*(dLnPhi_n(V,i,1) - 1.D0/x(V,2) - dLnPhi_n(V,i,2))
	  
	enddo
	F(1+NC) = lnPhi(S,2) - lnPhi(L,2) - lnX(L,2)
	F(2+NC) = dot_product(o, W) - spec
	
	return
	
endsubroutine FyJacobianoSLV
! -------------------------------------------------------------------------
! 
!  Cálculo del logaritmo natural de la  F U G A C I D A D  de cada componen-
!  te en cada fase, además de las derivadas respecto del número de moles, tem-
!  peratura y presión.
! 
subroutine fugacitySLV(NC, NG, NST, T, P, x, )
	implicit none

!  Cálculo de los coeficients de fugacidad y de sus derivadas respecto de P y T:	
	iGz = 0; iC = 0; iType = 0	
	call paragc(T, NC, NG, NST, iTemp)
	call gceos()
	
	return
endsubroutine fugacitySLV