! Cálculo de un punto de desviación de función objetivo para ajuste de puntos
! críticos a través de un enfoque semi-implícito.[1]
!
! [1] J.I. Ramello, Aspectos avanzados en el modelados del equilibrio entre fases 
!     de mezclas en sistemas asimétricos y a altas presiones. Tesis doctoral, UNS
!     2014, pp. 50-52.
! 
! dlnF_dz1   = derivada del ln(f1) con respecto a z1 a T y P constantes, calculada
!              analíticamente.
! d2lnF_dz12 = derivada segunda de ln(f1) respecto de z1 a T y P constantes, calcu-
!              lada por diferencias finitas centradas.
! 
subroutine implicitPCrit(NC, NG, NST, T, P, zn, dlnF_dz1, d2lnF_dz12, Z)

	implicit none
	integer, parameter ::          NCM=15, NGM=15, NGAM=24
	real(8), parameter ::          delta = 1.D-5
	
	integer ::                     NC, NG, NST, iErrPC
	
	real(8) ::                     T, P, Z, dlnF_dz1, d2lnF_dz12, dlnFper1, dlnFper2
	real(8), dimension(NCM) ::     zn, aux, aux2
	
! Cálculo de dlnF_dz1 a T, P y z:
	call dFug_z1(NC, NG, NST, T, P, zn, dlnF_dz1, Z)
! Cálculo de d2lnF_dz12 por diferencias finitas centradas:
	aux(1) = zn(1)*(1.D0 - delta)
	aux(2) = 1.D0 - aux(1)
	call dFug_z1(NC, NG, NST, T, P, aux, dlnFper1, aux2(1))
	
	aux(1) = zn(1)*(1.D0 + delta)
	aux(2) = 1.D0 - aux(1)
	call dFug_z1(NC, NG, NST, T, P, aux, dlnFper2, aux2(2))
	
	d2lnF_dz12 = (dlnFper2 - dlnFper1)/2.D0/delta
	
	return
	
endsubroutine implicitPCrit

subroutine dFug_z1(NC, NG, NST, T, P, zn, dlnF_dz1, Z)

	implicit none
	integer, parameter ::          NCM=15, NGM=15, NGAM=24

	integer ::                     NC, NG, NST, iErrPC, iC, Igz

	real(8) ::                     T, P, Z, dlnF_dz1, d2lnF_dz12
	real(8), dimension(NCM) ::     zn, lnPhi, dLnPhi_T, dLnPhi_P
	real(8), dimension(NCM,NCM) :: dLnPhi_n
	
	iGz = 0
! 	iTemp = 0
! 	iDer = 2  !Necesito derivadas respecto del nro de moles
	call PARAGC(T, NC, NG, NST, 0)
! 	write(*,*)NC, NG, NST, T, P, zn(:NC)
	call GCEOS(NC, NG, NST, 2, 0, T, P, zn, lnPhi, dLnPhi_n, dLnPhi_T, DLnPHI_P, Z, iGz, 0, IC)
	dlnf_dz1 = dLnPhi_n(1,1) + 1.D0/zn(1) - dLnPhi_n(1,2)	
! 	write(*,*)dlnf_dz1	
	
	return
endsubroutine dFug_z1