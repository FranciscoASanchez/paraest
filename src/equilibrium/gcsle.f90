!
!	paraest copy of GCSLE-2.0.0.f90 (SolubSol23/src), made by
!	tools/port_gcsle.py. Program SLE and readSoluteProperties replaced by
!	subroutine GCSLE (one point, option 0, at the end of this file).
!	Changes in the routines of the original:
!	  - root_from_newton: second_order reset at each call (initialized in
!	    the declaration, it kept the value of the previous point).
!	  - format (4G) -> (4G15.7) (gfortran needs the width).
!
!--------------------------------------------------------------------------------------------------------------------------
!
! Calculate pure-solute properties
!
  subroutine pure_solute (solid_fugacity_option, solid_ID, B_ID, nu, NTSS, model, NC, T, P, T_fus, P_fus, DH_fus, DH_subl, &
                        vS, DCp_fus, DB, Dv_fus, Tss, DHss, Y, lnPhisL, ZsL, PfusT)
	
	implicit none

	integer, parameter                     :: NCM = 30
	real(8), parameter                     :: R1 = 82.05D0, R2 = 8.314D0

	integer                                :: calc, guess = 0, i, ic, phase_type = 0
	real(8)                                :: B = 0, Cp = 0, DLnPhi_N(NCM,NCM), H = 0, lnPhiBL, ln_psubl, psubl, ZAL, ZBL
	real(8), dimension(NCM)                :: aux, lnPhi, DH_N, DLnPhi_P, DLnPhi_T
	                                   
	integer, intent(in)                    :: B_ID, model, NC, NTSS, nu(2), solid_fugacity_option, solid_ID
	real(8), intent(in)                    :: DB, DCp_fus, DH_fus,DH_subl, Dv_fus, P, P_fus, T, T_fus, vS
	real(8), intent(in), dimension(NTSS+1) :: DHss, Tss 
	
	real(8), intent(out)                   :: lnPhisL, Y, ZsL, PfusT

	aux = epsilon(1.D0)/1d3
	aux(solid_ID) = 1.D0
	calc = 1 ! fugacity only
	guess = 0 !start with default initial value
	
!	Even though I'm mixing different units of the R constant, there shouldn't be any problem
!	if units of H, T, P and v are respected since the U is non-dimensional.

	if (solid_fugacity_option == 0 .OR. solid_fugacity_option == 2) then
	  
	  Y = -DH_fus/R2*(1D0/T - 1D0/T_fus)
	  if (DCp_fus /= 0) then
	  
	    Y = Y + DCp_fus*dlog(T/T_fus)/R2 + (DCp_fus*T_fus + DB*T_fus**2/2)/R2*(1D0/T - 1D0/T_fus) + DB/2/R2*(T - T_fus)
	  
	  endif
	  if (Dv_fus /= 0) then
	    
	    PfusT = Y*R1*T/Dv_fus + P_fus
	    Y = Y - Dv_fus/R1/T*(P - P_fus)
	    
	  endif
!	  Addition of solid-solid transition enthalpies:

	  if (NTSS > 0) then	
	    do i = 1, NTSS
	      if(T < Tss(i))then
	    
	        Y = Y - DHss(i)/R2*(1D0/T - 1D0/Tss(i))
	  
	      else
	  
	        exit
	  
	      endif
	    enddo
	  endif
	  
!	  Calculation of hypothetic pure liquid fugacity:

	  phase_type = 0 !whatever phase is the most stable
	  call THERMO (model, NC, phase_type, calc, guess, T, P, aux, ZsL, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)
	  lnPhisL = lnPhi(solid_ID)
	  if (solid_fugacity_option == 2) then

	    phase_type = ic
	    guess = 2
	    ZAL = ZsL
	    call THERMO (model, NC, phase_type, calc, guess, T_fus, P, aux, ZsL, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)
	    lnPhisL = lnPhisL - lnPhi(solid_ID) !change in phiL of A due to DT
	    
	    guess = 0
	    aux(:NC) = 0
	    aux(B_ID) = 1
	    call THERMO (model, NC, phase_type, calc, guess, T, P, aux, ZBL, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)	    
	    lnPhiBL = lnPhi(B_ID)
	    ZsL = nu(1)*ZAL + nu(2)*ZBL !no change in molar volume during the formation of the complex as liquid. The volume ocupied by thee 
	    
	    guess = 2
	    call THERMO (model, NC, phase_type, calc, guess, T_fus, P, aux, ZBL, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)
	    lnPhiBL = lnPhiBL - lnPhi(B_ID) !change in phiL of B due to DT
	    	    
	    aux(solid_ID) = dfloat(nu(1))/(nu(1) + nu(2))	!stoichiometric composition at complex formation
	    aux(B_ID) = 1 - aux(solid_ID)
	    ZBL = ZsL/(nu(1) + nu(2))	!initial compressibility factor of the stoichiometric mixture
	    call THERMO (model, NC, phase_type, calc, guess, T_fus, P, aux, ZBL, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)
	    
! 	    Fugacity of the pure complex.
!	    Take into account that stoichiometric coeff are negative, but the input values are positive.
	    lnPhisL = nu(1)*(dlog(aux(solid_ID)) + lnPhi(solid_ID) + lnPhisL) + nu(2)*(dlog(aux(B_ID)) + lnPhi(B_ID) + lnPhiBL)
	    
	  endif	  
	  
	elseif (solid_fugacity_option == 1) then
	  
	  ln_psubl = dlog(P_fus) - DH_subl/R2*(1D0/T - 1D0/T_fus)
	  psubl = dexp(ln_psubl)
	  Y = vS*(P - psubl)/R1/T + ln_psubl - dlog(P)
	  
!	  Calculation of saturated solute fugacity coefficient 

	  phase_type = -1 !gas
	  call THERMO (model, NC, phase_type, calc, guess, T, psubl, aux, ZsL, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)
	  lnPhisL = lnPhi(solid_ID)
	  
	endif
	
	return
endsubroutine
!--------------------------------------------------------------------------------------------------------------------------
! 
!	Solve the isofugacity of the solute using a damped Newton procedure.
!
subroutine root_from_newton (maxIteration, model, NC, solid_fugacity_option, solid_ID, B_ID, nu, tolerance, lnFS, T, P, Z, x_free, x, lnPhi, success, F_error)

	implicit none

	logical                              :: second_order = .FALSE.
	
	integer, parameter                   :: NCM = 30
	real(8), parameter                   :: deltaMax = .25d0
	
	integer                              :: guess = 0, iteration, maxIteration
	real(8)                              :: deltaW, deltaW_Newton, deltaW_old, dFdw, dFdx, F, fDscld, fScld, &
	                                        test_oscilation, w, w_old, w_ref
	                                     
	integer, intent(in)                  :: B_ID, model, NC, nu(2), solid_fugacity_option, solid_ID	                                     
	real(8), intent(in)                  :: lnFS, P, T, tolerance
	real(8), dimension(NC-1), intent(in) :: x_free
	real(8), intent(inout)               :: x(NCM), Z
	
	logical, intent(inout)               :: F_error, success
	real(8), dimension(NCM), intent(out) :: lnPhi
	
	
!	First value of scaled molar fraction:	  
! 
! 	Variable scaling
	second_order = .FALSE.
	w_old = fScld (x(solid_ID), 1.D0, 0.D0, 1.D0)
	deltaW_old = 10
	do iteration = 1, maxIteration
  
	  ! Activate second order after 3 iterations.
	  if (iteration > 2) second_order = .TRUE.
	  call solute_isofugacity (model, guess, NC, solid_fugacity_option, solid_ID, B_ID, nu, lnFS, T, P, x_free, second_order, F, dFdx, Z, x, lnPhi, F_error)
	  if (F_error) exit
	    
	  if (second_order) then
	    
	    dFdw = dFdx*x(solid_ID)*(1D0 - x(solid_ID))
	    deltaW = -F/dFdw !Newton step
	    
	  else
	    
	    deltaW = F !Direct-substitution step
	    
	  endif
	  if (dabs(deltaW) < tolerance) then
	    
	    success = .TRUE.
	    exit
	    
	  endif    
	  w_ref = max(1D0, dabs(w_old)) ! the scaled molar fraction has a single problem: it can be 0 if x = 0.5.
	                                ! taking 1 as reference considers w = 1 the interval x = 0.27-0.73
	  deltaW = min(dabs(deltaW), deltaMax*w_ref) * deltaW/dabs(deltaW)
	  test_oscilation = dabs(deltaW + deltaW_old)
	  if (test_oscilation < tolerance) deltaW = deltaW/3
	  write (*, '(4G15.7)') x(solid_ID), F, deltaW, test_oscilation
	  w = w_old + deltaW
	  x(solid_ID) = fDscld (w, 1.D0, 0.D0, 1.D0)
	  call normal (NC, solid_ID, x(:NC), x_free)
	  w_old = w
	  deltaW_old = deltaW
  
	enddo
      
	return

endsubroutine root_from_newton

!--------------------------------------------------------------------------------------------------------------------------
!
! This subroutine aims to find a good initial estimate if the conventional has failed. It will be called whenever
! the Newton procedure fails.
! Screening searchs for a root by trial and error, calculating many points in an inverval set by the user. If many
! possible roots are found, the subroutine will work among the closer to the value initially proposed.
!
subroutine screening_SLE (n_screen, model, NC, solid_fugacity_option, solid_ID, B_ID, nu, lnFS, x0, T, P, x, x_free)
	implicit none
	
	integer, parameter                    :: NCM = 30
	                                      
	logical                               :: deriv = .FALSE., F_error
	integer                               :: guess = 0, i, j, n_candidates = 0, phase_type = 0
	                                      
	real(8)                               :: B = 0, Cp = 0, deltaW = 0, dFdx = 0, F = 0, fDscld, fScld, &
&	                                         w = 0, w0 = 0, Z = 0
	real(8), dimension(NCM)               :: lnPhi = 0
	                                      
	integer, dimension(:), allocatable    :: candidate
	real(8), dimension(:), allocatable    :: F_screen, w_screen
	                                      
	integer, intent(in)                   :: B_ID, model, n_screen, NC, nu(2), solid_fugacity_option, solid_ID
	real(8), intent(in)                   :: lnFS, T, P, x0
	real(8), intent(in), dimension(NC-1)  :: x_free
	
	real(8), intent(inout), dimension(NCM):: x

	w0 = fScld (x0, 1D0, 0D0, 1D0)
	allocate (w_screen(n_screen+1), F_screen(n_screen+1), candidate(n_screen+1))
	n_candidates = 0
	candidate(:n_screen+1) = 0
	
! Set interval

	x(solid_ID) = 0
	call normal (NC, solid_ID, x(:NC), x_free)
	call solute_isofugacity (model, guess, NC, solid_fugacity_option, solid_ID, B_ID, nu, lnFS, T, P, x_free, deriv, F_screen(1), dFdx, Z, x, lnPhi, F_error)   
	w_screen(1) = fScld(x(solid_ID), 1D0, 0.D0, 1.D0)
	write (*, *) x(solid_ID), F_screen(1)	

	deltaW = (36 - w_screen(1))/dble(n_screen)
	
! Start screening with logarithmic scaling, within the interval (0,1).

	do i = 1, n_screen
  
	  w_screen(i+1) = w_screen(i) + deltaW
	  x(solid_ID) = fDscld(w_screen(i+1), 1D0, 0.D0, 1.D0)
	  call normal (NC, solid_ID, x(:NC), x_free)
	  call solute_isofugacity (model, guess, NC, solid_fugacity_option, solid_ID, B_ID, nu, lnFS, T, P, x_free, deriv, F_screen(i+1), dFdx, Z, x, lnPhi, F_error)   	  
	  write (*, *) x(solid_ID), F_screen(i+1)
	  
! Look for a candidate:

	  if (F_screen(i)*F_screen(i+1) < 0D0) then
    
	    n_candidates = n_candidates + 1
	    candidate(n_candidates) = i
    
	    if (n_candidates > 1) then

! 	      check if we are moveing away from the initial value originally applied.
	      if (dabs(w0 - w_screen(i+1)) > dabs(w0 - w_screen(candidate(n_candidates-1))) ) then
	        n_candidates = n_candidates - 1
	        exit !we've started to move away from x0
	      endif
	    endif
	    
	  endif
	  
	enddo
	if (n_candidates == 0) then
  
! 	  No candidates have been found :-(
	  return
	  
	else
  
!	  Screening finished, start Bisection within the interval with candidate closer to x0
	  i = candidate(n_candidates)
	  do j = 1, n_screen/4
	
	    w = (w_screen(i) + w_screen(i+1))/2
	    x(solid_ID) = fDscld(w, 1D0, 0.D0, 1.D0)
	    call normal (NC, solid_ID, x(:NC), x_free)
	    call solute_isofugacity (model, guess, NC, solid_fugacity_option, solid_ID, B_ID, nu, lnFS, T, P, x_free, deriv, F, dFdx, Z, x, lnPhi, F_error)   
	    if (F*F_screen(i+1) < 0) then	    
	      w_screen(i) = w
	      F_screen(i) = F
	    else   
	      w_screen(i+1) = w
	      F_screen(i+1) = F	          
	    endif	  
	    
	  enddo
	
	endif
	
	return
endsubroutine screening_SLE
!--------------------------------------------------------------------------------------------------------------------------
subroutine normal (NC, solid_ID, x, x_free)

	implicit none
	real(8), parameter                   :: macheps = epsilon(1.D0)
	integer                              :: i
	
	integer, intent(in)                  :: NC, solid_ID
	real(8), dimension(NC-1), intent(in) :: x_free
	real(8), dimension(NC), intent(inout):: x
	

	do i = 1, NC - 1
	  if (i < solid_ID) then
	    x(i) = x_free(i)*(1 - x(solid_ID))
	  elseif (i >= solid_ID) then
	    x(i+1) = x_free(i)*(1 - x(solid_ID))
	  endif
	enddo
	do i = 1, NC
	  if (x(i) >= 1.D0) then

	    x(i) = 0.999999999999999D0

	  elseif (x(i) <= 0) then

	    x(i) = macheps/1000

	  endif
	enddo
	return
endsubroutine normal
!--------------------------------------------------------------------------------------------------------------------------
!
!
subroutine solute_isofugacity (model, guess, NC, solid_fugacity_option, solid_ID, B_ID, nu, lnFS, T, P, x_free, deriv, F, dFdx, Z, x, lnPhi, F_error)

	implicit none

	integer, parameter                    :: NCM = 30                                     
	                                      
	integer                               :: calc, check, phase_out = 0, guess, i, model, phase_type = 0
	real(8)                               :: B = 0, Cp = 0, H = 0
	real(8), dimension(NCM)               :: DLnPhi_N(NCM,NCM) = 0, DLnPhi_P = 0, DLnPhi_T = 0, DH_N = 0

	real(8), intent(inout)                :: x(NCM), lnPhi(NCM), Z
	                                     
	integer, intent(in)                   :: B_ID, NC, nu(2), solid_fugacity_option, solid_ID                                     
	logical, intent(in)                   :: deriv
	real(8), intent(in)                   :: lnFS, T, P
	real(8), dimension(NC-1), intent(in)  :: x_free
	                                      
	logical, intent(out)                  :: F_error
	real(8), intent(out)                  :: F, dFdx
	
	calc = 1
	if (deriv) then
	  if (calc == 5) guess = 2
	  calc = 5
	endif
	F_error = .FALSE.
	call THERMO (model, NC, phase_type, calc, guess, T, P, x, Z, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, phase_out)
	F = lnFS - nu(1)*(dlog(x(solid_ID)) + lnPhi(solid_ID))
	if (solid_fugacity_option == 2) then
	  
	  F = F - nu(2)*(dlog(x(B_ID)) + lnPhi(B_ID))
	  
	endif
	if (deriv) then
	  
	  dFdx = -dfloat(nu(1))/x(solid_ID)
	  if (solid_fugacity_option == 2) dFdx = dFdx + dfloat(nu(2))/(1 - x(solid_ID)) !sign is inverted here.
	  do i = 1, NC
	    
	    if (i < solid_ID) then
		
	      dFdx = dFdx + nu(1)*x_free(i)*DLnPhi_N(i,solid_ID)
	      if (solid_fugacity_option == 2) dFdx = dFdx + nu(2)*x_free(i)*DLnPhi_N(i,B_ID)
		
	    elseif (i == solid_ID) then
		
	      dFdx = dFdx - nu(1)*DLnPhi_N(solid_ID,solid_ID)
	      if (solid_fugacity_option == 2) dFdx = dFdx - nu(2)*DLnPhi_N(solid_ID,B_ID)				
		
	    elseif (i > solid_ID) then
		
	      dFdx = dFdx + nu(1)*x_free(i-1)*DLnPhi_N(i,solid_ID)
	      if (solid_fugacity_option == 2) dFdx = dFdx + nu(2)*x_free(i-1)*DLnPhi_N(i,B_ID)		
	    
	    endif
	  enddo
	  
	endif
	return
endsubroutine solute_isofugacity
!--------------------------------------------------------------------------------------------------------------------------
!
!	Función de escalado de variables: útil para variables con topes máximo y mínimo.
!
!	vMin y vMax son las cotas mínima y máxima que puede tomar "x". 
!
!	scl es la variable de escalado, que estira o comprime el dominio.
!
! double precision function fScld (x, scl, vMin, vMax)
! 	implicit none
! 	real(8) scl, vMax, vMin, x
! 	
! 	fScld = scl*dlog((x - vMin)/(vMax - x))
! 
! 	return
! endfunction
! !--------------------------------------------------------------------------------------------------------------------------
! !
! !	Función de desescalado de variables: inversa de la función fScld
! !
! double precision function fDscld (x, scl, vMin, vMax)
! 	implicit none
! 	real(8) ex, scl, vMax, vMin, x
! 
! 	ex = dexp(x/scl)
! 	fDscld = (vMin + vMax*ex)/(1.d0 + ex)
! 
! 	return
! endfunction

!--------------------------------------------------------------------------------------------------------------------------
!
!  paraest: solubility of the solid solute solid_ID at T (K), P (atm), in a
!  liquid whose solute-free composition is x_free(NC-1), as the loop of
!  program SLE does for one point, option 0 (subcooled liquid reference).
!  The program loops between Newton and screening until Newton converges;
!  here screening is tried once (it is deterministic, a second one would
!  give the same start). x0: initial estimate of the solute mole fraction
!  (<= 0: ideal solution). Melting data as read by paraest: J/mol, K, atm,
!  cm3/mol; NTSS solid-solid transitions (Tss, DHss).
!  Returns x (liquid composition), Z of the liquid and gamma, the activity
!  coefficient of the solute; complete miscibility (no solid at T, P): x = 1
!  for the solute. iErr = 1 if it did not converge.
!
subroutine GCSLE (NC, solid_ID, T, P, x_free, x0, DH_fus, ADCp_fus, BDCp_fus, Dv_fus, T_fus, P_fus, NTSS, Tss, DHss, &
                  Z, x, gamma, iErr)

	use thermo_model, only: model
	implicit none

	integer, parameter                   :: maxIteration = 10, NCM = 30
	real(8), parameter                   :: tolerance = 1.E-7

	integer, intent(in)                  :: NC, solid_ID, NTSS
	real(8), intent(in)                  :: T, P, x_free(NC-1), x0, DH_fus, ADCp_fus, BDCp_fus, Dv_fus, T_fus, P_fus, &
	                                        Tss(*), DHss(*)
	real(8), intent(out)                 :: Z, x(NCM), gamma
	integer, intent(out)                 :: iErr

	logical                              :: F_error, success
	integer                              :: B_ID, nu(2)
	real(8)                              :: lnFS, lnPhisL, PfusT, x_start, Y, ZsL, Tss_w(NTSS+1), DHss_w(NTSS+1)
	real(8)                              :: lnPhi(NCM)

	iErr = 0
	nu(1) = 1
	nu(2) = 0
	B_ID = solid_ID
	Tss_w = 0
	DHss_w = 0
	Tss_w(:NTSS) = Tss(:NTSS)
	DHss_w(:NTSS) = DHss(:NTSS)
	PfusT = P_fus
	x = 0
	lnPhi = 0
	Z = 0
	gamma = 0
	call pure_solute (0, solid_ID, B_ID, nu, NTSS, model, NC, T, P, T_fus, P_fus, DH_fus, &
	                  0D0, 0D0, ADCp_fus, BDCp_fus, Dv_fus, Tss_w, DHss_w, Y, lnPhisL, ZsL, PfusT)
	if ((T < T_fus .AND. Dv_fus == 0) .OR. (Dv_fus*(P - PfusT) > 0)) then

	  ! Solid fugacity divided by pressure (ln (fs/P))
	  lnFS = lnPhisL + Y
	  x_start = min(x0, 1D0)
	  if (x_start <= 0) x_start = dexp(Y)
	  x(solid_ID) = x_start
	  call normal (NC, solid_ID, x(:NC), x_free)
	  success = .FALSE.
	  F_error = .FALSE.
	  call root_from_newton (maxIteration, model, NC, 0, solid_ID, B_ID, nu, tolerance, lnFS, T, P, Z, x_free, x, lnPhi, &
	                         success, F_error)
	  if (.NOT. success .AND. .NOT. F_error) then
	    call screening_SLE (100, model, NC, 0, solid_ID, B_ID, nu, lnFS, x_start, T, P, x, x_free)
	    call root_from_newton (maxIteration, model, NC, 0, solid_ID, B_ID, nu, tolerance, lnFS, T, P, Z, x_free, x, lnPhi, &
	                           success, F_error)
	  endif
	  if (success) then
	    gamma = exp(lnPhi(solid_ID) - lnPhisL)
	  else
	    iErr = 1
	  endif

	else

	  ! complete miscibility
	  x(solid_ID) = 1
	  Z = ZsL
	  gamma = 1

	endif
	return
endsubroutine GCSLE
