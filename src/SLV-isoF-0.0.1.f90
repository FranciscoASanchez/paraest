
! References
! [1] S.B. Rodriguez-Reartes, M. Cismondi Duarte, M.S. Zabaloy, Computation of solid?fluid?fluid 
!     equilibria for binary asymmetric mixtures in wide ranges of conditions. J. Supercrit Fluids
!     57 (2011) 9-24.

subroutine SLV_isofugacidad()

	implicit none
	
	integer, parameter ::             NCM=15, NGM=15, NGAM=24, iTemp = 1
	real(8), parameter ::             delta = 1.D-5
	                                  
	integer ::                        NData, NC, NG, NST, iErrPC
	                                  
	real(8) ::                        T, P, Z, dlnF_dz1, d2lnF_dz12, dlnFper1, dlnFper2
	real(8), dimension(NData) ::      T, P, x, y
	real(8), dimension(NCM) ::        zn, aux, aux2		
	
	x(
	
	
	return
	
endsubroutine
! --------------------------------------------------------------------------------------
!
! This subroutine evaluates S, L and V fugacity at T and P.

subroutine fugacitySLV(NC, NG, NST, T, P, x, Dv_fus, C1, C2, C3)
	implicit none
	integer, parameter             :: NCM = 15, NGM = 15, NGAM = 24, NCM=15, NGM=15, NGAM=24, L = 2, V = 1, S = 3
	integer ::                        iGz, NC, NG, NST, status
	                                  
	real(8)                        :: P, T
	real(8), dimension(3,NCM) ::      x, lnX, lnPhi, Z
	real(8), dimension(NCM,NCM) ::    dLnPhi_n, dLnPhi_T, dLnPhi_P
	
! Vapor phase fugacity
	iGz = 0; iC = 0; iType = V - 2, iDer = 1
	call paragc(T, NC, NG, NST, iTemp)
	call gceos(NC, NG, NST, iDer, iTemp, T, P, x(V,:NCM), lnPhi(V,:NCM), dLnPhi_n, dLnPhi_T, dLnPhi_P, Z(V), IGz, iType, status)
	
! Liquid phase fugacity
	iGz = 0; iType = L - 1, iDer = 1
	call paragc(T, NC, NG, NST, iTemp)
	call gceos(NC, NG, NST, iDer, iTemp, T, P, x(L,:NCM), lnPhi(L,:NCM), dLnPhi_n, dLnPhi_T, dLnPhi_P, Z(L), IGz, iType, status)
	

! Pure 1 solid phase fugacity:
	iGz = 0; iType = S - 2, iDer = 1
	call paragc(T, NC, NG, NST, iTemp)
	call gceos(NC, NG, NST, iDer, iTemp, T, P, x(S,:NCM), lnPhi(S,:NCM), dLnPhi_n, dLnPhi_T, dLnPhi_P, Z(S), IGz, iType, status)
	
! Correct solid phase properties:
	Z(S) = Z(S) - P*Dv_fus/R/T
	Tr = T/T_fus
! Se Eq. (A-11) of the supplementary material of reference [1]
	lnPhi(S,1) = lnPhi(S,1) - Dv_fus/R/T_fus*(C1*(1.D0 - Tr) + C2*(Tr - 1.D0 + dlog(Tr)) + C3*(Tr/2.D0 - 1.D0 + 0.5D0/Tr) + (P - P_fus)/Tr)
	
	
	
	
	
	return
endsubroutine fugacitySLV