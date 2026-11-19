! Cálculo de un punto de desviación de función objetivo para ajuste de puntos
! críticos a través de un enfoque semi-implícito.[1]
!
! Si el NC = 1, el programa calcula dP/dv y d2P/dv2 para incluir en el ajuste el 
! punto crítico del puro.
!
! Si el NC = 2, se calculan las derivadas de la fugacidad respecto de la composi-
! ción, para ajustar locus crítico de sistema binario.
!
! 
! dP_dv      = derivada de P respecto de v a T contante para un componente puro.
! d2P_dv2    = derivada segunda de P respecto de v dos veces a T constante, para
!              un componente puro.
! dlnF_dz1   = derivada del ln(f1) con respecto a z1 a T y P constantes, calculada
!              analíticamente.
! d2lnF_dz12 = derivada segunda de ln(f1) respecto de z1 a T y P constantes, calcu-
!              lada por diferencias finitas centradas.
! D1, D2     = Asignación de las derivadas dependienco de si NC = 1 ó 2.
!
! Estas variables son luego asignadas a las de salida, dependiendo de NC = 1 ó 2.
!
!
! Referencias:
!
! [1] J.I. Ramello, Aspectos avanzados en el modelados del equilibrio entre fases 
!     de mezclas en sistemas asimétricos y a altas presiones. Tesis doctoral, UNS
!     2014, pp. 50-52.
! 
! [2] Derivada de la función inversa. 
!     https://es.wikipedia.org/wiki/Derivada_de_la_función_inversa#Derivadas_superiores
! 
subroutine implicitPCrit(NC, NG, NST, T, P, zn, dlnF_dz1, d2lnF_dz12, Z)

	implicit none
	integer, parameter ::            NCM=15, NGM=15, NGAM=24
	real(8), parameter ::            delta = 1.D-5
	                                 
	integer ::                       NC, NG, NST, iErrPC, phase_index, phase_type = 0, guess = 0
	                                 
	real(8) ::                       T, P, Z, dlnF_dz1, d2lnF_dz12, dlnFper1, dlnFper2, Z2(2)
	real(8), dimension(NCM) ::       zn, aux, aux2
	                                 
!	Variables incluidas en vectores COMMON:
	real(8) ::                       dP_dv(2), lambda1, lambda2, lambda3, Prep(2), q_total
	real(8), dimension(NGM) ::       ny_total, theta
	real(8), dimension(2,NGM) ::     Help2, Help4, Help5, Help7, Help12
	real(8), dimension(2,NGM,NGM) :: tau
	
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON /ZZZ/                     Help2, Help4, Help5, ny_total, q_total,theta, Help12, Help15, &
						   Help7, tau, Prep, dPdV, lambda1, lambda2, lambda3
	
	D1 = 0
	D2 = 0
	iErrPC = 0
	if (NC == 1) then
	
! 	  Cálculo de las derivadas respecto de volumen
!
!	  La derivada primera se obtiene analíticamente del calculo del volumen:
	  call ZMAX (Z2, phase_type, T, P, zn, guess, NC, NG, NST, phase_index)
	  phase_numer = 1+(phase_index+1)/2
	  Z = Z2(phase_number)
	  D1 = dP_dv(phase_number)
	  
! 	  La derivada 2da por diferencias finitas centradas:

!	  Es necesario utilizar el teorema de la función inversa, aplicado 
!       a derivdas de orden superior.
!	  En particular, en nuestro caso se tiene que:
!			
!		d2y/dx2 = - d2x/dy2 * (dy/dx)^3   [2]
!
!	  Reempladando y = P, y x = v:
!	    
! 		d2P/dv2 = - d2v/dP2 * (dP/dv)^3 
!
!	  Mientras que la d2v/dP2 nuevamente gracias al teorema de la función 
!       inversa, por diferencias finitas:
!
! 		d2v/dP2 = (dv1/dP - dv2/dP)/2/delta/P

!	  donde
!
!	 	dv/dP = 1/(dP/dv)

	  guess = 2
	  phase_type = phase_index
	  P_pert = P*(1D0 + delta)
	  call ZMAX (Z2, phase_type, T, P, zn, guess, NC, NG, NST, phase_index)
	  dv1_dP = 1D0/dP_dv(phase_number)
	  P_pert = P*(1D0 - delta)
	  call ZMAX (Z2, phase_type, T, P, zn, guess, NC, NG, NST, phase_index)
	  dv2_dP =  1D0/dP_dv(phase_number)
	  
	  D2 = - (dv1_dP - dv2_dP)/2/delta/P * D1**3
	  
	elseif (NC == 2) then
	
! 	  Cálculo de dlnF_dz1 a T, P y z:
	  call dFug_z1 (NC, NG, NST, T, P, zn, dlnF_dz1, Z)
  !	  Cálculo de d2lnF_dz12 por diferencias finitas centradas:
	  aux(1) = zn(1)*(1.D0 - delta)
	  aux(2) = 1.D0 - aux(1)
	  call dFug_z1 (NC, NG, NST, T, P, aux, dlnFper1, aux2(1))
	  
	  aux(1) = zn(1)*(1.D0 + delta)
	  aux(2) = 1.D0 - aux(1)
	  call dFug_z1 (NC, NG, NST, T, P, aux, dlnFper2, aux2(2))
	  
	  d2lnF_dz12 = (dlnFper2 - dlnFper1)/2.D0/delta
	
	  D1 = dlnF_dz1
	  D2 = d2lnF_dz12
	  
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