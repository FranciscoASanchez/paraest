!--------------------------------------------------------------------------------------------------------------------------
!
!  paraest: parameter fitting for the GCA-EoS.
!
!  The molecular structure can be defined differently for the dispersive and
!  the association contributions.
!
!  The critical diameter dc of each component can be:
!    - estimated from the critical point (explicit equation, no fit);
!    - fitted to one point of the component's vapor pressure curve;
!    - fitted like any other parameter, i.e. to all the data (binary and/or
!      pure-component) in which the component takes part.
!
!  The input is read from paraest.dat; the fitted parameters are scaled
!  (module scaling) before being passed to the optimizer.
!
!--------------------------------------------------------------------------------------------------------------------------
program paraest
	use dimensions
	use fit_data
	use gca_parameters
	use scaling, only: XXSCAL
	implicit none
!
!	Required by the IMSL routine ZXSSQ
	EXTERNAL FUNC
!
!	The data shared with FUNC (system, model parameters, experimental data,
!	weights, identification of the fitted parameters) are in modules fit_data
!	and gca_parameters. The meaning of IDXX is in module scaling.
!
!	Locals:
	INTEGER NSIG, NPRED, NPAR, NCXP, IOPT, MAXFN, IXJAC, INFER, IER
	INTEGER NTEXT(NSetM,65)              ! title of each data set
	DOUBLE PRECISION XX(NPARM)           ! fitted parameters (scaled during the fit)
	DOUBLE PRECISION PARM(4)             ! ZXSSQ parameters
	DOUBLE PRECISION EPS, DELTAZ, SSQ, boundFrac
	DOUBLE PRECISION XJAC(NDM,NPARM), XJTJ(NPARM*(NPARM+1)/2), WORK(LWORK)
	DOUBLE PRECISION FCION(NDM)          ! residuals
	DOUBLE PRECISION stDvD(NSetM)
	INTEGER iGr
	DOUBLE PRECISION epxOut               ! epx of spline_temperature (not used here)
!
!	IDC(J,1): component whose experimental mole fractions are given in data set j.
!	IDC(J,2): component whose mole fractions are NOT given.
!
!	Read the input file paraest.dat
	CALL read_input(NSIG,NPRED,NPAR,PARM,NTEXT,XX,NCXP,IOPT)
!
!	Initial dc from the critical point predicted by the GC-EoS (without
!	association), for every component whose dc was not given by the user.
	CALL init_dc(NC,PRCRIT,TEMPCRIT,DC)
!
!	Spline temperatures of the initial parameters (for the output only: the
!	model computes them in fakePARGR).
	do iGr = 1, NG
	  call spline_temperature(TSTR(iGr), GDOT(iGr), GDDOT(iGr), TSPLINE(iGr), epxOut)
	enddo
!
!	Write the input data and initial parameters to paraest.out
	CALL write_output_header(NSIG,NPRED,NPAR,PARM,NTEXT,XX,IOPT)
!
!	Type of each value in XX, to assign it the right scaling factor
	CALL set_param_types(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,numTc,NUMDC)
!
!	Scale the values of XX. In a fit, initial values on a bound (or very
!	close) are moved 0.1 % of the interval away, so the optimizer can move
!	them; in prediction they are used as given.
	IF(NPRED.EQ.1)THEN
	  boundFrac = 0D0
	ELSE
	  boundFrac = 1D-3
	ENDIF
	CALL XXSCAL(NPAR,XX,IDXX,boundFrac)
!
!	Prediction or fit
	IF(NPRED.EQ.1)THEN
!	  Prediction
	  CALL FUNC(XX,NCXP,NPAR,FCION)
!
	ELSE
!	  Parameter fit. ZXSSQ parameters:
	  EPS = 1.D-6
	  DELTAZ = 1.D-6
!	  Maximum number of FUNC evaluations. ZXSSQ spends NPAR+1 evaluations per
!	  iteration on the Jacobian alone (finite differences): with a fixed
!	  MAXFN = 100, 8 parameters gave ~10 iterations and the fit ended with
!	  IER = 133. Now it allows ~50 iterations, whatever is fitted.
	  MAXFN = 50*(NPAR + 1)
	  PARM(3)=2.D2
	  PARM(4)=.1D0
	  IXJAC = NCXP
!
	  CALL ZXSSQ(FUNC,NCXP,NPAR,NSIG,EPS,DELTAZ,MAXFN,IOPT,PARM,XX,SSQ,FCION,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
!
!	  Statistics of the fit
	  CALL write_fit_stats(SSQ,XJAC,WORK,INFER,IER)
!
!	  Spline temperatures of the final parameters (for the output)
	  do iGr = 1, NG
	    call spline_temperature(TSTR(iGr), GDOT(iGr), GDDOT(iGr), TSPLINE(iGr), epxOut)
	  enddo
!
	END IF
!
!	Final parameters, and experimental and calculated values
	CALL write_results(NTEXT,stDvD)
!
endprogram paraest
!
!-------------------------------------------------------------------------------------------------------------------------
!
!	Reads the input file paraest.dat. It also echoes what it reads, with
!	format, to the screen, to help find an error in the input.
!
SUBROUTINE read_input(NSIG,NPRED,NPAR,PARM,NTEXT,XX,NCXP,IOPT)
	use dimensions
	use fit_data, only: NC, W1, W2, W3, W4, W5, w6, W7, W8, W9, W10, W11, W12, w13, ISAT, IGLE, IVE, IUNI, &
&	                    PRCRIT, TEMPCRIT, NSAT, IDC, NOBS, iHEN, mez, IFL, rGamExp1, rGamExp2, rGamExp3, &
&	                    prop_E_calc, prop_E_exp, TEXP, PEXP, XEXP, xexp2, YEXP, yExp2, rhoExp, IDSAT, IOBJSAT, &
&	                    wSat, wBin, wTer, wEx, wGam, DHfu, Tfu, DCpfu, Dvfu, TtrSS, DHss, Pfu
	use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&	                    ENASS, RKASS, NNP, NNQ, NALFA, NENASS, NRKASS, TSAT, PSAT, dc, IDG, ITEMP, NIDG, IDIDG, &
&	                    iDQ, IDALF, NIDALF, IDIDALF, IDENASS, NIDENASS, IDIDENASS, IDRKASS, NIDRKASS, &
&	                    IDIDRKASS, NUMDC, IDDC, numTc, idTc
	use scaling, only: dc_mean, error_Dc
	use utilities, only: i2s, rep, check_limit
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
	DIMENSION PARM(4)                ! ZXSSQ parameters
	DIMENSION NTEXT(NSetM,65)        ! title of each data set
	DIMENSION XX(NPARM)              ! initial values of the fitted parameters
	dimension error_Dc1(NCM)
!
	OPEN(UNIT=1,STATUS='OLD',FILE='paraest.dat',FORM='FORMATTED',ACCESS='SEQUENTIAL')
	READ (1, *) W1
	read (1, *) W2, W3, W4, W5, W6, W7
	read (1, *) W8, W9, W10
	read (1, *) W11, W12
	read (1, *) w13
	! IOPT is read as a real because data files usually have "2." (ifx accepts
	! it for an integer, gfortran does not).
	read (1, *) PARM(1),PARM(2),aIOPT
	IOPT = nint(aIOPT)
	write(*,'(15F10.4,I2)')W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,w13,PARM(1),PARM(2),IOPT
!
1	FORMAT(8F10.4)
2	FORMAT(*(I3))
!
	READ(1,*) NC,NG,NGA,NSIG
	write(*,'(4I3)')NC,NG,NGA,NSIG
	READ(1,*) ISAT,IGLE,IVE,IUNI,NPRED
	write(*,'(5I3,/)')ISAT,IGLE,IVE,Iuni,NPRED
!
!	Count the equilibrium data points
	NCXP = 0
	IF(ISAT.NE.0)then

	  READ(1,*) (IDSAT(I), I=1,ISAT)
	  write(*,'(*(I3,X))') (IDSAT(I), I=1,ISAT)
	  READ(1,*) (NSAT(IDSAT(I)), I=1,ISAT)
	  write(*,'(*(I3,X))')(NSAT(IDSAT(I)), I=1,ISAT)
	  READ(1,*) (IOBJSAT(IDSAT(I)), I=1,ISAT)
	  write(*,'('//rep(iSat,"I3,X")//'/)') (IOBJSAT(IDSAT(I)), I=1,ISAT)
!
	  DO I=1,ISAT

	    NSUMK = 0
	    do K=1,(IDSAT(I)-1)

	      NSUMK = NSUMK + 1

	    end do
	    NSUML = 0
	    do L=1,NSUMK

	      NSUML = NSUML + NSAT(L)

	    end do
		!Pure-component equilibrium data: T and P, or T, P and rho.
	    DO J= (NSUML +1 ),( NSUML + NSAT(IDSAT(I)) )

	      if(iObjSat(idSat(i)).NE.81)then
	        READ(1,*) TEXP(J),PEXP(J),wSat(j)
	        write(*,'(F7.2,F12.8,I12)') TEXP(J),PEXP(J), j
	      else
	        READ(1,*) TEXP(J),PEXP(J),rhoExp(j),wSat(j)
	        write(*,'(F7.2,2(F12.8,X))') TEXP(J),PEXP(J),rhoExp(j)
		  endif

	    END DO
!
	    NCXP = NCXP + NSAT(IDSAT(I)) ! number of data points read so far
!
	  END DO

	endif
	write(*,'(//)')
!	Area and reference temperature of each dispersive group
	DO I=1,NG

	  READ(1,*) QAREA(I), TSTR(I)
	  write(*,'(F6.3,F7.2)')qArea(i),Tstr(i)

	END DO
!	Segment energy of each group
	READ(1,*) (GSTR(I), I=1,NG)
	write(*,'(/,'//rep(NG,"F8.0,2X")//'"   !g*")')(GSTR(I), I=1,NG)
!
	READ(1,*) (GDOT(I), I=1,NG)
	write(*,'('//rep(NG,"F7.4,2X")//'"   !g`")')(GDOT(I), I=1,NG)
!
	READ(1,*) (GDDOT(I), I=1,NG)
	write(*,'('//rep(NG,"F7.4,2X")//'"   !g``",//)')(GDdOT(I), I=1,NG)
!
!	Binary interaction parameters: k*
	DO I=1,NG

	  READ(1,*) (RKIJ(I,J), J=1,NG)
	  if(i.EQ.1)then

	    write(*,'(/,'//rep(NG,"F6.4,2X")//'"  !k*")')(RKIJ(I,J), J=1,NG)

	  else

	    write(*,'(*(F6.4,2X))')(RKIJ(I,J), J=1,NG)

	  endif


	END DO
!	k': temperature dependence of k*
	DO I=1,NG

	  READ(1,*) (RKTIJ(I,J), J=1,NG)
	  if(i.EQ.1)then

	    write(*,'(/,'//rep(NG,"F6.4,2X")//'"  !k`")')(RKTIJ(I,J), J=1,NG)

	  else

	    write(*,'(*(F6.4,2X))')(RKTIJ(I,J), J=1,NG)

	  endif

	END DO
!	alpha: non-randomness parameters
	DO I=1,NG

	  READ(1,*) (ALFAIJ(I,J), J=1,NG)
	  if(i.EQ.1)then

	    write(*,'(/,'//rep(NG,"F8.4,2X")//'"  !alfa")')(alfaIJ(I,J), J=1,NG)

	  else

	    write(*,'(*(F8.4,2X))')(alfaIJ(I,J), J=1,NG)

	  endif

	END DO
!
	IF(NGA.GT.0)THEN
!
!	  Association contribution
!	  Number of association sites of each group:
	  READ(1,*) (MASSOC(I),I=1,NGA)
	  write(*,'(//,'//rep(NGA,"I2,2X")//'5X,"! Association sites per group",/)') (MASSOC(I),I=1,NGA)

!	  Association energy and volume, on the same line
	  NST = 0
	  EnAss(:6,:12,:6,:12) = 0D0
	  rKAss(:6,:12,:6,:12) = 0D0
	  do i = 1,NGA
	    do j = i,NGA
	      do k = 1,Massoc(i)
	        do l = k,Massoc(j)
!
	          NST = NST + 1
	          read(1,*)EnAss(k,i,l,j),rKAss(k,i,l,j)
		    if(NST.EQ.1)then

	            write(*,'(F6.1,4X,F8.5,4X,"! ",2I2,3X,2I2,4X,"Association energy and volume")')EnAss(k,i,l,j),rKAss(k,i,l,j)

	          else

	            write(*,'(F6.1,4X,F8.5,4X,"! ",2I2,3X,2I2,4X)')EnAss(k,i,l,j),rKAss(k,i,l,j)

	          endif
!	          Symmetries
	          EnAss(l,i,k,j) = EnAss(k,i,l,j)
	          EnAss(l,j,k,i) = EnAss(k,i,l,j)
	          EnAss(k,j,l,i) = EnAss(k,i,l,j)
	          rKAss(l,i,k,j) = rKAss(k,i,l,j)
	          rKAss(l,j,k,i) = rKAss(k,i,l,j)
	          rKAss(k,j,l,i) = rKAss(k,i,l,j)
!
	        enddo
	      enddo
	      write(*,'(/)')
	    enddo
	  enddo
	endif
!-----------------------------------------------------------------------
!
	write(*,'(//)')
!	Pure-component properties:
	DO I=1,NC

	  READ (1,*) TEMPCRIT(I), PRCRIT(I), TSAT(I), PSAT(I), DC(I), Tfu(i), DHfu(i), DCpfu(i), TtrSS(i), DHss(i), DVfu(i), Pfu(i)
	  write(*,'(2(F7.2,2X,F12.6,2X),F9.5,2X,"! Comp. No.",I2)')TEMPCRIT(I),PRCRIT(I),TSAT(I),PSAT(I),DC(I),i

	END DO
	write(*,'(//)')
!
!	Group composition, dispersive contribution:
	DO I=1,NG
!
!	  NU(I,J) = number of groups of type i in component j.
	  READ(1,*) (NU(I,J), J=1,NC)
	  if(i.EQ.1)then

	    write(*,'(//,'//rep(NC,"I2,4X")//'"!  Disp. group No.",I2,"   Dispersive group composition")')(NU(I,J), J=1,NC),I

	  else

	    write(*,'('//rep(NC,"I2,4X")//'"!  Disp. group No.",I2)')(NU(I,J), J=1,NC),I

	  endif
!
	END DO
!
!
	IF(NGA.GT.0) THEN
!
!	  Group composition, association contribution:
	  DO I=1,NGA
!
!	    ONUASOC(I,J) = number of associating groups of type i in component j.
	    READ(1,*) (ONUASOC(I,J), J=1,NC)
	    if(i.EQ.1)then

	      write(*,'(/,'//rep(NC,"F2.0,4X")//'"!  Assoc. group No.",I2,"   Associating group composition")')(oNUAsoc(I,J), J=1,NC),I

	    else

	      write(*,'('//rep(NC,"F2.0,4X")//'"!  Assoc. group No.",I2)')(oNUAsoc(I,J), J=1,NC),I

	    endif

	  END DO

	END IF
!
!
!	Data sets (binary or ternary equilibrium, excess properties, ...)
	DO I=1,IGLE
!
	  read(1,*)IHEN(I)
	  IF(IHEN(I).EQ.0) THEN

	    write(*,'(//,I2,5X,"Binary equilibrium data")')iHen(i)
	    READ(1,*)IFL(I)
	    write(*,'(i2)')iFL(i)
	    if(iFL(i).LT.3)then

	      mez(i) = 2
	      READ(1,*) NOBS(I),IDC(I,1),IDC(I,2)
	      write(*,'(3(I2,4X))')NOBS(I),IDC(I,1),IDC(I,2)

	    else

	      mez(i) = 3
	      READ(1,*) NOBS(I),(IDC(I,j),j=1,mez(i))
	      write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,mez(i))

	    endif

	  elseif(iHen(i).EQ.1)then
!
!	    Excess property data
	    read(1,*)ifl(i)
	    if(ifl(i).EQ.1)then

	      write(*,'(//,I2,5X,"Excess volume data")')ifl(i)

	    elseif(ifl(i).EQ.2)then

	      write(*,'(//,I2,5X,"Excess enthalpy data")')ifl(i)

	    endif
	    read(1,*)nObs(i), (idC(i,j), j = 1, 2)
	    write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)

	  elseif(iHen(i).EQ.2)then
	    
	    write(*,'(//,I2,5X,"Solid solubility data")')iHen(i)
	    read(1,*)nObs(i), (idC(i,j), j = 1, 2)
	    write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)
	    
	  elseif(iHen(i) == 3)then
	  
	    write(*,'(//, I2, 5X, "Critical point data")') iHen(i)
	    read(1,*)nObs(i), (idC(i,j), j = 1, 2)
	    write(*,'(4(I2,4X))')NOBS(I),(IDC(I,j),j=1,2)	    
	    
	  ELSEif(ihen(i).EQ.5)then

	    write(*,'(//,I2,5X,"Activity coefficient data")')iHen(i)
	    read(1,*) mez(i)
	    write(*,'(I2)')mez(i)
	    read(1,*) nObs(i), (idC(i,j),j=1,mez(i))
	    write(*,'(*(I2,2X))')(idC(i,j),j=1,mez(i)),nObs(i)

	  END IF
        READ(1,23) (NTEXT(I,J),J=1,65)
	  write(*,23)(NTEXT(I,J),J=1,65)
23	  FORMAT(40A2)
!
	  write(*,'(/)')
	  DO J=( NCXP+1 ),(NCXP + NOBS(I))

	    If(ihen(i).EQ.0)then
!
!	      ****************** Phase equilibrium *****************
	      if(mez(i).EQ.2)then

!	        Binary
	        READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J),wBin(:2,j)
	        write(*,'(F7.2,2X,F9.4,2(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),YEXP(J),j

	      elseif(mez(i).EQ.3)then

!	        Ternary
	        READ(1,*) TEXP(J),PEXP(J),XEXP(J),YEXP(J),xExp2(j),yExp2(j),wBin(:2,j),wTer(:2,j)
	        write(*,'(F7.2,2X,F9.4,4(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),YEXP(J),xExp2(j),yExp2(j),j

	      endif
!
	    elseif(iHen(i).EQ.1)then
!
!	      ****************** Excess properties *****************
	      READ(1,*) TEXP(J),PEXP(J),XEXP(J),prop_E_exp(J),wEx(j)
	      write(*,'(F7.2,2X,F9.4,2X,F7.4,F12.6, "! ", I3)')TEXP(J),PEXP(J),XEXP(J),prop_E_exp(J), j
		
	    elseif(iHen(i).EQ.2)then
		
!	      ****************** SLE *****************
	      read(1,*)Texp(j), Pexp(j), xExp(j), wEx(j)
	      write(*,'(F7.2,2X,F9.4,2X,F7.4, "! ", I3)')TEXP(J),PEXP(J),XEXP(J), j
	      
	    elseif(iHen(i) == 3)then
		
!	      ****************** Critical point *****************
	      read(1,*)Texp(j), Pexp(j), xExp(j), wBin(:2,j)
	      write(*,'(F7.2,2X,F9.4,2X,F7.4, "! ", i3)')TEXP(J), PEXP(J), XEXP(J), j     

	    elseif(iHen(i).EQ.5)then
!
!	      ****************** Activity coefficients *****************
	      if(mez(i).EQ.2)then

	        read(1,*) TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),wGam(:2,j)
	        write(*,'(F7.2,2X,F9.4,3(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),j

	      elseif(mez(i).EQ.3)then

	        read(1,*) TEXP(J),PEXP(J),XEXP(J),xexp2(J),rgamexp1(J),rgamexp2(J),rgamexp3(J),wGam(:3,j)
	        write(*,'(F7.2,2X,F9.4,4(2X,F7.4),3X,"! ",I3)')TEXP(J),PEXP(J),XEXP(J),rgamexp1(J),rgamexp2(J),rgamexp3(J),j

	      endif
!
	    end if

	  END DO
	  NCXP = NCXP + NOBS(I)
!
!	  End of the loop over data sets
	END DO
!
	write(*,'(//)')
	READ(1,*) NNP,NNQ,NALFA,NENASS,NRKASS,numTc, NUMDC
	write(*,'(7(I2,4X),"! Interaction - Area - Non-random. - Assoc. energy - Assoc. volume - Crit. temp. - Crit. diameters")') NNP,NNQ,NALFA,NENASS,NRKASS,numTc, NUMDC
	call check_limit(NNP,    NPTIPOM, 'g/k interaction parameters (NNP)')
	call check_limit(NNQ,    NPTIPOM, 'area parameters (NNQ)')
	call check_limit(NALFA,  NPTIPOM, 'non-randomness parameters (NALFA)')
	call check_limit(NENASS, NPTIPOM, 'association energies (NENASS)')
	call check_limit(NRKASS, NPTIPOM, 'association volumes (NRKASS)')
	call check_limit(numTc,  NPTIPOM, 'critical temperatures (numTc)')
	call check_limit(NUMDC,  NPTIPOM, 'critical diameters (NUMDC)')
!
	if(NNP.NE.0) then
!
!	  Dispersive interaction parameters to fit: g(i,i) and k(i,j)
	  READ(1,*) (IDG(I,1),IDG(I,2),I=1,NNP)
	  write(*,'(*(2(I2,X),4X))')(IDG(I,1),IDG(I,2),I=1,NNP)
!
	  READ(1,*) (ITEMP(I), I=1,NNP)
	  write(*,'(*(I2,6X))')(ITEMP(I), I=1,NNP)
!
	  READ(1,*) (NIDG(I), I=1,NNP)
	  call check_limit(sum(NIDG(:NNP)), NIDM, 'identical g/k parameters')
	  write(*,'(*(I2,6X))')(NIDG(I), I=1,NNP)
!
	  NSSS = 0
!	  Identical parameters:
	  DO L=1,NNP
!
	    if(NIDG(L).NE.0) then
!
	      READ(1,*) (IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
	      write(*,'(*(I2,X))')(IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
	      NSSS = NSSS + NIDG(L)
!
	    end if
!
	  END DO

	end if
!
!	Area parameters to fit
	if(NNQ.NE.0)then

	  READ(1,*) (IDQ(I), I=1,NNQ)
	  write(*,'(*(I2,X))') (IDQ(I), I=1,NNQ)

	endif
!
!	Non-randomness parameters to fit
	if(NALFA.NE.0) then
!
	  READ(1,*) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
	  write(*,'(*(2(I2,X),4X))')(IDALF(I,1),IDALF(I,2),I=1,NALFA)
!
	  READ(1,*) (NIDALF(I),I=1,NALFA)
	  call check_limit(sum(NIDALF(:NALFA)), NIDM, 'identical alpha parameters')
	  write(*,'(*(I2,6X))')(NIDALF(I),I=1,NALFA)
	  NSSA = 0
	  DO L=1,NALFA

	    if(NIDALF(L).NE.0) then

	      READ(1,*) (IDIDALF(J,1),IDIDALF(J,2),J=(1+NSSA),(NIDALF(L)+NSSA))
	      write(*,'(*(2I2,X))')
	      NSSA = NSSA + NIDALF(L)

	    end if

	  END DO

	end if
!
!	Association energies to fit
	if(NENASS.NE.0) then

	  DO I=1,NENASS

	    READ(1,*) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
	    write(*,'(4I3)')IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)

	  END DO
	  READ(1,*) (NIDENASS(I), I=1,NENASS )
	  call check_limit(sum(NIDENASS(:NENASS)), NIDM, 'identical association energies')
	  write(*,'(*(I2,X))')(NIDENASS(I), I=1,NENASS )
!
!	  Identical ones:
	  NSSE = 0
	  DO L=1,NENASS

	    if(NIDENASS(L).NE.0) then

	      DO J=(1+NSSE),(NIDENASS(L)+NSSE)

	        READ(1,*) IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)
	        write(*,'(4I3)') IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)

	      END DO
	      NSSE = NSSE + NIDENASS(L)

	    end if

	  END DO
!
	end if
!
!	Association volumes to fit:
	if(NRKASS.NE.0) then
!
	  DO I=1,NRKASS

	    READ(1,*) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
	    write(*,'(4I3)') IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)

	  END DO
	  READ(1,*) (NIDRKASS(I), I=1,NRKASS )
	  call check_limit(sum(NIDRKASS(:NRKASS)), NIDM, 'identical association volumes')
	  write(*,'(*(I2,X))')(NIDRKASS(I), I=1,NRKASS )
!
	  NSSRK = 0
	  DO L=1,NRKASS

	    if(NIDRKASS(L).NE.0) then

	      DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

	        READ(1,*) IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)
	        write(*,'(4I3)')IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)

	      END DO
	      NSSRK = NSSRK + NIDRKASS(L)

	    end if

	  END DO

	end if
!
!	Components whose Tc is fitted:
	IF(numTc.NE.0)then

	  READ(1,*) (idTc(I), I=1,numtc)
	  write(*,'('//rep(numTc,"I2,X")//'/)')(idTc(I), I=1,numtc)

	endif
!
!	Components whose dc is fitted (IDDC), and the range allowed around the
!	initial value of each one (error_Dc1).
	IF(NUMDC.NE.0)then

	  READ(1,*) (IDDC(I), I=1,numDC)
	  READ(1,*) (error_Dc1(i), I=1,numDC)
	  write(*,*) (IDDC(I), I=1,numDC)
	  write(*,*) (error_Dc1(i), I=1,numDC)

	endif
!
!	Number of parameters
	NPAR = NNP+NNQ+NALFA+NENASS+NRKASS+numTc+NUMDC
	IF(NNP.NE.0) THEN

	  DO I = 1,NNP

	    NPAR = NPAR + ITEMP(I)

	  END DO

	END IF
!
!	Initial values of the parameters
	if(nPar > 0)READ(1,*)  (XX(I),I=1,NPAR)
!
!	The dc values are also kept apart, for their scaling between bounds:
	if(numDc > 0)then
	
	  dc_mean(nPar - numDc + 1:nPar) = xx(nPar - numDc + 1:nPar)
	  error_Dc(nPar - numDc + 1:nPar) = error_Dc1(1:numDc)
	  
	endif
!
	CLOSE(UNIT=1)
!
ENDsubroutine read_input
!------------------------------------------------------------------------------------------------------------------------------------
!
!	Writes to paraest.out the data and parameters read by read_input (not
!	the results).
!
SUBROUTINE write_output_header(NSIG,NPRED,NPAR,PARM,NTEXT,XX,IOPT)
	use dimensions
	use fit_data, only: NC, W1, W2, W3, W4, W5, w6, W7, W8, W9, W10, W11, W12, w13, ISAT, IGLE, IVE, IUNI, &
&	                    PRCRIT, TEMPCRIT, NSAT, IDC, NOBS, IHEN, mez, IFL, TEXP, PEXP, XEXP, xexp2, YEXP, IDSAT
	use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&	                    ENASS, RKASS, NNP, NNQ, NALFA, NENASS, NRKASS, TSAT, PSAT, IDG, ITEMP, NIDG, IDIDG, &
&	                    IDQ, IDALF, NIDALF, IDIDALF, IDENASS, NIDENASS, IDIDENASS, IDRKASS, NIDRKASS, &
&	                    IDIDRKASS, DC, NUMDC, IDDC, TSPLINE, numtc, idtc
	use fit_data, only: DHfu, Tfu, DCpfu, DVfu, TtrSS, DHss, Pfu
	use utilities, only: i2s, rep
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
	DIMENSION PARM(4)                ! ZXSSQ parameters
	DIMENSION NTEXT(NSetM,65)        ! title of each data set
	DIMENSION XX(NPARM)              ! initial values of the fitted parameters
!
	OPEN(UNIT=1,STATUS='unknown',file='paraest.out',FORM='FORMATTED')
!
	WRITE(1,1) W1,W2,W3,W4,W5,w6,W7,W8,W9,W10,W11,W12,w13,PARM(1),PARM(2),IOPT
!
1	FORMAT(/,X,' Weights in the objective function:',/,                                              &
&	       X,'Weight of vapor pressures:                       ',F10.4,/,                            &
&	       X,'Weight of bubble pressure of binary mixtures:    ',F10.4,/,                            &
&	       X,'Weight of bubble temperature of binary mixtures: ',F10.4,/,                            &
&	       X,'Weight of composition of comp. 1 in phase II:    ',F10.4,/,                            &
&	       X,'Weight of composition of comp. 1 in phase I:     ',F10.4,/,                            &
&	       X,'Weight of composition of comp. 2 in phase I:     ',F10.4,' (ternary mixtures only)',/  &
&	       X,'Weight of composition of comp. 2 in phase II:    ',F10.4,' (ternary mixtures only)',/  &
&	       X,'Weight of activity coefficient of comp. 1:       ',F10.4,/                             &
&	       X,'Weight of activity coefficient of comp. 2:       ',F10.4,/                             &
&	       X,'Weight of activity coefficient of comp. 3:       ',F10.4,' (ternary mixtures only)',/  &
&	       x,'Weight of excess volume:                         ',F10.4,/                             &
&	       x,'Weight of excess enthalpy:                       ',F10.4,/                             &
&	       x,'Weight of solid solubility:                      ',F10.4,/                             &
&	       X,' PARM(1) = ',F10.4,'     PARM(2) = ',F10.4,/,                                          &
&	       X,' Optimization algorithm: IOPT = ',I3)
!
	WRITE(1,'(/,X," Ratio between g*(i,i) and g(i,i) at which the ",/,&
&               X," spline is made = 4 ")')
!
	WRITE(1,'(/,X,"--------------------------------------------------",/,&
&	          X,"--------------- INITIAL PARAMETERS ---------------",/,  &
&	          X,"--------------------------------------------------",/,  &
&	          X," Pure-group parameters (dispersive term)",/)')

	WRITE(1,"(40X,'COMPONENTS',/,X,' Gr.No.:   T*        Q        TSPLINE ',"//rep(NC,'I3')//"/)") (I, I=1,NC)
!
	DO I=1,NG

	  WRITE(1,4) I,TSTR(I),QAREA(I),TSPLINE(I),(NU(I,J), J=1,NC)
4	  FORMAT(X,I3,5X,2F10.4,X,E10.4,*(I3))

	END DO
!	NU(I,J) = number of groups of type i in component j
	WRITE(1,5)
 5	FORMAT(/,X,' Interaction energies gjj ',/,X,'(the functional form changes at high temperatures)',/)
!
	DO I=1,NG
!
	  WRITE(1,6) I,I,GSTR(I), GDOT(I), GDDOT(I)
6	  FORMAT(X,' g(',I2,I2,') = ',F10.1,' * [ 1 + ',F10.4,' *(T/T* - 1) + ',F10.4, ' *LN(T/T*)] ',/)

	END DO
	WRITE(1,7)
 7	FORMAT(/,X,' Binary parameters: [k*] ',/)
	DO I=1,NG

	  WRITE(1,8) (RKIJ(I,J), J=1,NG)
8	  FORMAT(X,*(F10.4))
	END DO
!
	WRITE(1,9)
 9	FORMAT(/,X,' Temperature dependence: [k`] ',/)
	DO I=1,NG

	  WRITE(1,8) (RKTIJ(I,J), J=1,NG)

	END DO
	WRITE(1,10)
 10	FORMAT(/,X,' Binary non-randomness parameters: [alpha] ',/)
	DO I=1,NG

	  WRITE(1,8) (ALFAIJ(I,J), J=1,NG)

	END DO
!
	WRITE(1,186) NGA
 186	FORMAT(/,X,' Number of groups in the association contribution = ',I3)
!
	IF(NGA.GT.0) THEN

	  WRITE(1,"(/,X,' Association parameters ',/,23X,'Components',/,' Group 	  # Assoc sites',"//rep(NC,'I5')//"/)") (I, I=1,NC)
	  DO I=1,NGA

	    WRITE(1,12) I, MASSOC(I), (ONUASOC(I,J), J=1,NC)
12	    FORMAT(3X,I3,9X,I3,6X,*(F5.2))

	  END DO
	  WRITE(1,13)
13	  FORMAT(/,X,'Association energy parameters (K)  ',/)
!
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then
	      DO L=1,NGA

              IF(MASSOC(L).NE.0)then
	            DO I=1,MASSOC(J)

	              DO K=1,MASSOC(L)

	                WRITE(1,14) I,J,K,L,ENASS(I,J,K,L)
14	                FORMAT(X,' ENASS(',I2,I2,X,I2,I2,')  = ',F10.4)

	              END DO

	            END DO
	          endif

	      END DO
	    endif

	  END DO
!
!	  Association volume parameters (cm3/mol)
	  WRITE(1,15)
15	  FORMAT(/,X,'Association volume parameters (cm3/mol) ',/)
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then
	      DO L=1,NGA

	        IF(MASSOC(L).NE.0)then
	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,16) I,J,K,L,RKASS(I,J,K,L)
16	              FORMAT(X,' RKASS(',I2,I2,X,I2,I2,') = ',F16.10)

	            END DO

	          END DO
	        endif

	      END DO
	    endif

	  END DO

	END IF
	WRITE(1,18)
18	FORMAT(/,X,'Pure-component data')
	WRITE(1,19)
19	FORMAT (/, 4X, 'Tc/K      Pc/atm    Tsat/K    Psat/atm', 4X, 'dc/(cm/mol^(1/3))', X, "Tfus/K", X, "DHfus/(J/mol)", X, &
&              "DCpfus/(J/molK)", X, "Ttr/K", X, "DHtr/(J/mol)", X, "DVfus/(cm3/mol)", X, "Pfus/atm"/,X,' (initial dc estimated from the exptl. critical point, or given) ',/)
!
	DO I=1,NC

	  WRITE(1,500) TEMPCRIT(I), PRCRIT(I), TSAT(I), PSAT(I), DC(I), Tfu(i), DHfu(i), DCpfu(i), TtrSS(i), DHss(i), DVfu(i), Pfu(i)
500	  FORMAT(X, F7.2, 3X, F10.4, 2X,F7.2, X, G12.4, 6X, F7.4, 7X, F7.2, X, G12.4, 5X, F6.1, 5X, F7.2, 3(X, G12.4) )

	END DO
!
	WRITE(1,201)
201	FORMAT(/,X,' Parameters to estimate ',/)
!
	if(NNP.NE.0) then
!
	  WRITE(1,202)
202	  FORMAT(' Parameters gii or kij ')
!
	  WRITE(1,203) (IDG(I,1),IDG(I,2),I=1,NNP)
203	  FORMAT(X,*(I3))
!
	  WRITE(1,204) (ITEMP(I), I=1,NNP)
204	  FORMAT(X,' Temperature dependence ',/,X,*(I3))
!
	  WRITE(1,206) (NIDG(I), I=1,NNP)
206	  FORMAT(X,' Number of identical par. ',/,X,*(I3))
!
	  WRITE(1,208)
208	  FORMAT(X,' Identical parameters ')
!
	  NSSS = 0
!
	  DO L=1,NNP

	    if(NIDG(L).ne.0) then
!
	      WRITE(1,210) (IDIDG(J,1), IDIDG(J,2),J=(1+NSSS),(NIDG(L)+NSSS))
210	      FORMAT(X,*(I3))
	      NSSS = NSSS + NIDG(L)
!
	    end if

	  END DO

	end if
!
	if(NNQ.NE.0) WRITE(1,211) (IDQ(I), I=1,NNQ)
211	FORMAT(/,X,' Area parameters  ',/,X,*(I3))
!
	if(NALFA.ne.0) then

	  WRITE(1,213) (IDALF(I,1),IDALF(I,2),I=1,NALFA)
213	  FORMAT(/,X,' Non-randomness parameters ',/,X,*(I3))
!
	  WRITE(1,215) (NIDALF(I),I=1,NALFA)
215	  FORMAT(X,' Number of identical par. ',/,X,*(I3))
!
	  NSSA = 0
	  DO L=1,NALFA

          if(NIDALF(L).ne.0) then

	      WRITE(1,217) (IDIDALF(J,1),IDIDALF(J,2),J=(1+NSSA),(NIDALF(L)+NSSA))
217	      FORMAT(X,' Identical : ',*(I3))
	      NSSA = NSSA + NIDALF(L)

	    end if

	  END DO

	end if
!
	if(NENASS.ne.0) then

	WRITE(1,219)
219	FORMAT(/,X,' Association energy parameters ')
!
	  DO I=1,NENASS

	    WRITE(1,221) IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)
221	    FORMAT(X,4I3)

	  END DO
	  WRITE(1,223) (NIDENASS(I), I=1,NENASS )
223	  FORMAT(X,' Number of identical par. ',*(I3))
	  NSSE = 0
	  DO L=1,NENASS

	    if(NIDENASS(L).ne.0) then

	      DO J=(1+NSSE),(NIDENASS(L)+NSSE)

	        WRITE(1,225) IDIDENASS(J,1), IDIDENASS(J,2),IDIDENASS(J,3), IDIDENASS(J,4)
225	        FORMAT(' Identical: ',4I5)

	      END DO
	      NSSE = NSSE + NIDENASS(L)

	    end if

	  END DO

	end if
!
	if(NRKASS.ne.0) then

	  wRITE(1,503)
503	  FORMAT(/,X,' Association volume parameters ')

	  DO I=1,NRKASS

	    WRITE(1,229) IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)
229	    FORMAT(X,4I3)

	  END DO
	  WRITE(1,231) (NIDRKASS(I), I=1,NRKASS )
231	  FORMAT(X,' Number of identical parameters ',*(I3))
	  NSSRK = 0
	  DO L=1,NRKASS

	    if(NIDRKASS(L).ne.0) then

	      DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

	        WRITE(1,232) IDIDRKASS(J,1), IDIDRKASS(J,2),IDIDRKASS(J,3), IDIDRKASS(J,4)
232	        FORMAT(X,' Identical: ',4I3)

	      END DO
	      NSSRK = NSSRK + NIDRKASS(L)

	    end if

	  END DO

	end if
	IF(numtc.ne.0) WRITE(1,233) (idtc(I), I=1,numTc)
233	FORMAT(/,X,' Critical temperatures ',/,X,*(I3))
	IF(NUMDC.NE.0) WRITE(1,234) (IDDC(I), I=1,NUMDC)
234   FORMAT(/,X,' Hard-sphere diameters ',/,X,*(I3))
!
	WRITE(1,235)  (XX(I),I=1,NPAR)
 235	FORMAT(/,X,' Initial values of the fitted parameters ',/,X,*(F13.4))
!
	WRITE(1,91)
91	FORMAT(//,X,' Course of the optimization ',X)
	if (NPAR > 0) then
	  WRITE(1,"(//,X,' Parameter No. (SCALED). ',/,X,"//rep(NPAR,'I15')//"X,'        OBJ FUN ',/)") (I,I=1,NPAR)
	else   ! in prediction there are no parameters: only the title
	  WRITE(1,"(//,X,' Parameter No. (SCALED). ',/)")
	endif
!
	CLOSE(UNIT=1)
!
ENDsubroutine write_output_header
!------------------------------------------------------------------------------------------------------------------------------------
!
!	Sets the model parameters from the new values of the fitted parameters,
!	which come unscaled in XXAUX.
!
SUBROUTINE set_model_params(XXAUX)
	use dimensions
	use fit_data, only: NC, W1, W2, W3, W4, W5, w6, W7, W8, ISAT, IGLE, IVE, IUNI, PRCRIT, TEMPCRIT, NSAT, &
&	                    IDC, NOBS, IHEN, mez, IFL, TEXP, PEXP, XEXP, xexp2, YEXP, IDSAT
	use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&	                    ENASS, RKASS, NNP, NNQ, NALFA, NENASS, NRKASS, TSAT, PSAT, IDG, ITEMP, NIDG, IDIDG, &
&	                    iDQ, IDALF, NIDALF, IDIDALF, IDENASS, NIDENASS, IDIDENASS, IDRKASS, NIDRKASS, &
&	                    IDIDRKASS, numtc, idTc, NUMDC, IDDC, DC
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
	DIMENSION XXAUX(NPARM)               ! values of the fitted parameters
!
	NPAR2 = 0
!
!	Parameters of type gii or kij
	IF(NNP.NE.0) then !A

	  DO I=1,NNP

	    IF(IDG(I,1).EQ.IDG(I,2)) THEN

	      GSTR(IDG(I,1)) = XXAUX(NPAR2+1)
	      NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.0)then

	        GDOT(IDG(I,1)) = 0.D0
	        GDDOT(IDG(I,1)) = 0.D0

	      endif
	      IF(ITEMP(I).GE.1) GDOT(IDG(I,1)) = XXAUX(NPAR2+1)
	      IF(ITEMP(I).GE.1) NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.1) GDDOT(IDG(I,1)) = 0.D0
	      IF(ITEMP(I).EQ.2) GDDOT(IDG(I,1)) = XXAUX(NPAR2+1)
	      IF(ITEMP(I).EQ.2) NPAR2 = NPAR2 + 1
!
	    ELSE

	      RKIJ(IDG(I,1),IDG(I,2)) = XXAUX(NPAR2+1)
	      NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.0) RKTIJ(IDG(I,1),IDG(I,2))= 0.D0
1001	      IF(ITEMP(I).EQ.1) RKTIJ(IDG(I,1),IDG(I,2))=XXAUX(NPAR2+1)
	      IF(ITEMP(I).EQ.1) NPAR2 = NPAR2 + 1
	      IF(ITEMP(I).EQ.2) THEN

!	        Avoid human error
	        ITEMP(I) = 1
	        GO TO 1001

	      END IF
!
	    END IF

	  END DO

	endif
!
!	Identical ones of type gii or kij
	IF(NNP.NE.0) THEN

	  NSSS = 0
	  DO I=1,NNP

	    if(NIDG(I).ne.0) then

	      IF(IDG(I,1).EQ.IDG(I,2)) THEN

	        DO J=(1+NSSS),(NIDG(I)+NSSS)

	          GSTR(IDIDG(J,1))  = GSTR(IDG(I,1))
	          GDOT(IDIDG(J,1))  = GDOT(IDG(I,1))
	          GDDOT(IDIDG(J,1)) = GDDOT(IDG(I,1))

 	        END DO
	        NSSS = NSSS + NIDG(I)

	      ELSE

	        DO J=(1+NSSS),(NIDG(I)+NSSS)

	          RKIJ(IDIDG(J,1),IDIDG(J,2)) = RKIJ(IDG(I,1),IDG(I,2))
	          RKTIJ(IDIDG(J,1),IDIDG(J,2))= RKTIJ(IDG(I,1),IDG(I,2))

	        END DO
	        NSSS = NSSS + NIDG(I)
!
	      END IF

	    end if

	  END DO

	END IF
!
!	Area parameters
	IF(NNQ.NE.0) THEN

	  DO I=1,NNQ

	    QAREA(IDQ(I)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
!	Alpha parameters
	IF(NALFA.NE.0) THEN

	  DO I=1,NALFA

	    ALFAIJ(IDALF(I,1),IDALF(I,2)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
!	Identical alphas
	IF(NALFA.NE.0) THEN

	  NSSA = 0
	  DO I=1,NALFA

	    if(NIDALF(I).ne.0) then

	      DO J=(1+NSSA),(NIDALF(I)+NSSA)

	        ALFAIJ(IDIDALF(J,1),IDIDALF(J,2))= ALFAIJ(IDALF(I,1),IDALF(I,2))

	      END DO
	      NSSA = NSSA + NIDALF(I)

	    end if

	  END DO

	END IF
!
!	Association energy parameters
	IF(NENASS.NE.0) THEN
!
	  DO I=1,NENASS

	    ENASS(IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
!	Identical association energies
	IF(NENASS.NE.0) THEN

	  NSSE = 0
	  DO L=1,NENASS

	    if(NIDENASS(L).ne.0) then

	      DO J=(1+NSSE),(NIDENASS(L)+NSSE)

	        ENASS(IDIDENASS(J,1),IDIDENASS(J,2),IDIDENASS(J,3),IDIDENASS(J,4)) = ENASS(IDENASS(L,1),IDENASS(L,2) ,IDENASS(L,3),IDENASS(L,4))

	      END DO
	      NSSE = NSSE + NIDENASS(L)

	    end if

	  END DO

	END IF
!
!	Association volume parameters
	if(NRKASS.ne.0) then

	  DO I=1,NRKASS

	    RKASS(IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	end if
!
!	Identical association volumes
	IF(NRKASS.ne.0) THEN

	  NSSRK = 0
	  DO L=1,NRKASS

	    if(NIDRKASS(L).ne.0) then

	      DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

	        RKASS(IDIDRKASS(J,1),IDIDRKASS(J,2),IDIDRKASS(J,3),IDIDRKASS(J,4)) = RKASS(IDRKASS(L,1),IDRKASS(L,2),IDRKASS(L,3),IDRKASS(L,4))

	      END DO
	      NSSRK = NSSRK + NIDRKASS(L)

	    end if

	  END DO

	END IF
!
!	Critical temperatures and diameters, when they are fitted parameters
	IF(numtc.NE.0) THEN

	  DO I=1,numtc

	    tempCrit(idtc(I)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF

	IF(NUMDC.NE.0) THEN

	  DO I=1,NUMDC

	    DC(IDDC(I)) = XXAUX(NPAR2+1)
	    NPAR2 = NPAR2 + 1

	  END DO

	END IF
!
ENDsubroutine set_model_params
!
!-------------------------------------------------------------------------------------------------------------------------
!
!	Residuals of the objective function (called by ZXSSQ).
!
SUBROUTINE FUNC(XX,NCXP,NPAR,FCION)
	use dimensions
	use fit_data
	use gca_parameters, PSATu => PSAT   ! PSAT is also a subroutine called here
	use utilities, only: i2s, rep
	use scaling, only: XXUNSCAL
!
	IMPLICIT REAL*8 (A-H,O-Z)
	real(8), parameter :: Rgas_SI = 8.314D0
!
!	The system and experimental data and the identification of the
!	parameters come through modules fit_data and gca_parameters (ZXSSQ passes
!	only XX, NCXP, NPAR and FCION).
!
	DIMENSION XX(NPAR), FCION(NCXP)
	DIMENSION XXAUX(NPARM)                 ! unscaled parameters
	DIMENSION stDvD(NSetM)
	DIMENSION PM(NCM), WWPM(2)
!	Work arrays for pure-component properties
	DIMENSION WDC(1), WTEMPCRIT(1), WPRCRIT(1), NUW(NGM,1), ONUASOCW(12,1), wZ(2)
!	Work arrays for binary mixture properties
	DIMENSION WWDC(3), WWTEMPCRIT(3), WWPRCRIT(3), NUWW(NGM,3), ONUASOCWW(12,3)
!	Arrays passed to the GCA package: dimensioned with its maximum number of components
	integer, parameter :: NCMgca = 30
!	Isofugacity of binary data
	DIMENSION xIsof(2), yIsof(2), apartIsof(2)
	DIMENSION feed(NCMgca), zFrac(NCMgca,4), help(NCMgca,NCMgca), zVal(4), iVap(1), indF(4)
!
	PENALIZ = 0	!1e-0
!	NaN guard (added 2026-11-06): NaN/Inf residuals are replaced by
!	PENALIZNAN so they do not spoil the objective function.
	PENALIZNAN = 1.D0
!
!	Unscale XX into XXAUX: FUNC must not alter the values of XX (see the
!	IMSL documentation).
	CALL XXUNSCAL(NPAR,XX,IDXX,XXAUX)
!
!	Set the model parameters from the unscaled values
	CALL set_model_params(XXAUX)
!
!	******************************** Residuals ****************************************
!
!	Critical diameters: each dc that is not a fitted parameter is adjusted to
!	reproduce exactly the vapor pressure point (TSAT, PSAT) of its component.
	if (any(TSAT(:NC) /= 0)) then
	  call fakePARGR(NC,NG,NGA,NST,IOUt,NU(:NG,:NC),QAREA(:NG),TStr(:NG),gStr(:NG),GDOT(:NG),GDDOT(:NG),   &
&	                 RKIJ(:NG,:NG),RKTIJ(:NG,:NG),ALFAIJ(:NG,:NG),ONUASOC(:NGA,:NC),MASSOC(:NGA),           &
&	                 ENASS(:2,:NGA,:2,:NGA),RKASS(:2,:NGA,:2,:NGA))
	  call fakePARMOL(NC,NG,NST,IOUt,DC(:NC),TEMPCRIT(:NC),PRCRIT(:NC))
	  DO I=1,NC
	    IF (TSAT(I) == 0) CYCLE
	    IF (any(IDDC(:NUMDC) == I)) CYCLE
	    call fit_dc(I,NC,NG,NST,TSAT(I),PSATu(I),dcFit,iErrDc)
	    IF (iErrDc == 0) THEN
	      DC(I) = dcFit
	    ELSE
	      WRITE(6,'(" dc fit failed for component ",I2,", ier = ",I2)') I, iErrDc
	    ENDIF
	    WRITE(6,'(" End of dc fit, component No. ",I2)')I
	  END DO
	endif
!
!
!	Pure components: vapor pressure or isofugacity
!	NUMEXP counts the experimental data points
	NUMEXP=0
	DO I=1,NC

	  IF(NSAT(I).NE.0) THEN
	    WDC(1)= DC(I)
	    WTEMPCRIT(1)= TEMPCRIT(I)
 	    WPRCRIT(1)  = PRCRIT(I)
	    DO J=1,NG

	      NUW(J,1) = NU(J,I)

	    END DO
	    DO K=1,NGA

	      ONUASOCW(K,1) = ONUASOC(K,I)

	    END DO
!
!       **************************** Loop over the NSAT points ****************************
!
		DO J=1,NSAT(I)
!
	        WTEMP=TEXP(NUMEXP+1)
	        WPSAT=PEXP(NUMEXP+1)
	        IF(IOBJSAT(I).EQ.50) THEN

!	          Objective function: isofugacity
	          write(6,'("Point ",I3," of component",I3)')j,i
	          call fakePARGR(1,NG,NGA,NST,IOUt,NUW(:NG,1),QAREA(:NG),TStr(:NG),gStr(:NG),GDOT(:NG),GDDOT(:NG), &
&	                         RKIJ(:NG,:NG),RKTIJ(:NG,:NG),ALFAIJ(:NG,:NG),ONUASOCW(:NGA,1),MASSOC(:NGA),       &
&	                         ENASS(:2,:NGA,:2,:NGA),RKASS(:2,:NGA,:2,:NGA))
	          call fakePARMOL(1,NG,NST,IOUt,WDC,WTEMPCRIT,WPRCRIT)
	          call isofug_pure(NG,NST,WTEMP,WPSAT,APART,WZFASL,WZFASV,IERRYAF)
!
	          IF(IERRYAF.EQ.0) THEN

	            PCALC(NUMEXP+1)=WPSAT
	            ZLCALC(NUMEXP+1)= WZFASL
	            ZVCALC(NUMEXP+1)= WZFASV
	            FCAUX(NUMEXP+1)= APART
	            FCION(NUMEXP+1)= W1 * FCAUX(NUMEXP+1) * wSat(numExp+1)

	          ELSE

	            PCALC(NUMEXP+1)= 0.
	            ZLCALC(NUMEXP+1)= 0.
	            ZVCALC(NUMEXP+1)= 0.
	            FCAUX(NUMEXP+1)= 0.
	            FCION(NUMEXP+1)= PENALIZ

	          END IF
!

	        END IF
	        IF(IOBJSAT(I).GE.80) THEN
!
!	          Objective function: vapor pressure (+ density if IOBJSAT = 81)
	          write(6,'("Point ",I3," of component",I3)')j,i
	          call fakePARGR(1,NG,NGA,NST,IOUt,NUW(:NG,1),QAREA(:NG),TStr(:NG),gStr(:NG),GDOT(:NG),GDDOT(:NG), &
&			                 RKIJ(:NG,:NG),RKTIJ(:NG,:NG),ALFAIJ(:NG,:NG),ONUASOCW(:NGA,1),MASSOC(:NGA),         &
&	                         ENASS(:2,:NGA,:2,:NGA),RKASS(:2,:NGA,:2,:NGA))
	          call fakePARMOL(1,NG,NST,IOUt,WDC,WTEMPCRIT,wPRCRIT)
	          call paraGC(wTemp, 1, NG, NST, 1)
	          iGP = 1
	          iGz = 0 ! could use the result of the dc fit...
	          iErrPSat = 0
			  iComp = 1
	          if (wPsat > 1E-10) then
! 	            Standard vapor pressure calculation
	            call PSAT(wTemp,wPSat,wZ,1,NG,NST,iComp,IGP,IGZ,iErrPSat)
	          else 
	          
! 	            Eval psat(T) = phiL(T)*(1 atm)
	            feed(1) = 1
	            feed(2:NCMgca) = 0
	            wPsat = 1
	            call GCEOS(1, NG, NST, 1, 0, wTemp, wPsat, feed, zFrac(:NCMgca,1), help(:NCMgca,:NCMgca), zfrac(:NCMgca,2), zfrac(:NCMgca,3), wZ(2), 0, 1, iErr)
	            wPsat = dexp(zFrac(1,1))	            
	            wZ(1) = 1
	            wZ(2) = wZ(2)*wPsat
	          
	          endif
	          if (iErrPSat >= 0) then
! 	            Vaporization enthalpy calculation
	            feed(1) = 1
	            feed(2:NCMgca) = 0
	            call NTLPY (1, NG, NST, wTemp, wPSat, feed, HRV, wZ(1), -1, -1, iErr)
	            call NTLPY (1, NG, NST, wTemp, wPSat, feed, HRL, wZ(2), 1, 1, iErr)
			
	          endif	            
	          if (wTemp < Tfu(i)) then
!
!                 Sublimation pressure
	            TR = (TFU(i)/wTemp - 1)
	            Yu = -DHFU(i)/Rgas_SI/TFU(i) * TR
	            DHSL = DHFU(i)
	            wZ(2) = wZ(2) - wPSat*DVfu(i)/82.05/wTemp
	            if (DCpfu(i) > 0) then
! 	            
! 	              Take into account change in heat capacity
	              Yu = Yu + DCpfu(i)/Rgas_SI * (TR + dlog(wTemp/TFU(i)))
	              DHSL = DHFU(i) + DCpfu(i)*(wTemp - Tfu(i))
	              
	            endif 
	            if (TtrSS(i) > 0 .AND. DHss(i) > 0) then
	            
! 	              Take into account SS transitions
	              if (wTemp < TtrSS(i)) then
	                Yu = Yu - DHss(i)/Rgas_SI/TtrSS(i) * (TtrSS(i)/wTemp - 1)
	                DHSL = DHSL + DHss(i)
	              endif
	              
	            endif  
	            WPSAT = WPSAT*dexp(Yu)
	            HRL = HRL - DHSL/.101325D0 !HRL is in [atm cm3], and it's not a liquid but a solid.

	          endif  
!
	          IF (IERRPSAT >= 0) THEN
!
	            PCALC(NUMEXP+1)=WPSAT
	            ZLCALC(NUMEXP+1)= WZ(2)
	            ZVCALC(NUMEXP+1)= WZ(1)
	            DHVL(numExp+1) = HRV - HRL
 	            FCAUX(NUMEXP+1)= (PCALC(NUMEXP+1)/PEXP(NUMEXP+1)-1)
	            FCION(NUMEXP+1)= W1 * FCAUX(NUMEXP+1) * wSat(numExp+1)*PEXP(NUMEXP+1)
	            if(iObjSat(i).EQ.81)fcion(numExp+1) =fcion(numExp+1) + abs(wPsat/82.05D0/wTemp/wZ(2)/rhoExp(numExp+1) - 1)

	          ELSE
!
	            PCALC(NUMEXP+1)= 0.
	            ZLCALC(NUMEXP+1)= 0.
	            ZVCALC(NUMEXP+1)= 0.
 	            FCAUX(NUMEXP+1)= 0.
	            FCION(NUMEXP+1)= PENALIZ
	            DHVL(numExp+1) = 0
!
	          END IF

	        END IF
	        NUMEXP = NUMEXP + 1
!
	      END DO
!
!         ********************** End of the loop over the NSAT points ***********************
	  END IF

	END DO
!
!
!	****************** Loop over the data sets ************************
!
	IF(IGLE.NE.0) THEN
	  DO KS=1,IGLE

!	    Data of the 2 components of the binary (3 for ternary data)
	    I1 = IDC(KS,1)
	    I2 = IDC(KS,2)
	    WWDC(1) = DC(I1)
	    WWDC(2) = DC(I2)
	    WWTEMPCRIT(1) = TEMPCRIT(I1)
	    WWTEMPCRIT(2) = TEMPCRIT(I2)
	    WWPRCRIT(1) = PRCRIT(I1)
	    WWPRCRIT(2) = PRCRIT(I2)
	    WWPM(1)=PM(I1)
	    WWPM(2)=PM(I2)
	    wwTfu = Tfu(I1)
	    wwDHfu = DHfu(I1)
	    wwDCpfu = DCpfu(I1)
	    wwDHss = DHss(I1)
	    wwTtrSS = TtrSS(I1)
	    wwDVfu = DVfu(I1)
	    wwPfu = Pfu(I1)

	    DO J=1,NG

	      NUWW(J,1) = NU(J,I1)
	      NUWW(J,2) = NU(J,I2)

	    END DO
	    DO J=1,NGA

	      ONUASOCWW(J,1) = ONUASOC(J,I1)
	      ONUASOCWW(J,2) = ONUASOC(J,I2)

	    END DO
	    NCw = 2
	    If(mez(ks)==3)then

	      NCw = 3
	      WWDC(3) = DC(IDC(KS,3))
	      WWTEMPCRIT(3) = TEMPCRIT(IDC(KS,3))
	      WWPRCRIT(3) = PRCRIT(IDC(KS,3))
	      DO J=1,NG

	        NUWW(J,3) = NU(J,IDC(KS,3))

	      END DO
	      DO J=1,NGA

	        ONUASOCWW(J,3) = ONUASOC(J,IDC(KS,3))

	      END DO

	    end if
!
!	    Load the model with these components: fakePARGR and fakePARMOL replace
!	    PARGR and PARMOL of the GCA package (group and component parameters).
	    call fakePARGR(NCw,NG,NGA,NST,IOUt,NUWW(:NG,:NCw),QAREA(:NG),TStr(:NG),GStr(:NG),GDOT(:NG),GDDOT(:NG), &
&	                   RKIJ(:NG,:NG),RKTIJ(:NG,:NG),ALFAIJ(:NG,:NG),ONUASOCWW(:NGA,:NCw),MASSOC(:NGA),         &
&	                   ENASS(:2,:NGA,:2,:NGA),RKASS(:2,:NGA,:2,:NGA))
	    call fakePARMOL(NCw,NG,NST,IOUt,WWDC(:NCw),WWTEMPCRIT(:NCw),wwPRCRIT(:NCw))
!
!---------Fluid phase equilibrium: VLE and LLE
	    IF(IHEN(KS).EQ.0) THEN   ! TPXY (IHEN=0)
!
	      NFAILURE = NOBS(KS)
	      STDVA(KS) = 0D0
	      STDVB(KS) = 0
	      DO I=1,NOBS(KS)

	        TEMP  = TEXP(NUMEXP+1)
	        PRESS = PEXP(NUMEXP+1)
	        X1 = XEXP(NUMEXP+1)
	        Y1 = YEXP(NUMEXP+1)
!	        IFL = -1 or -2: isofugacity
	        IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then
!
	          write(6,*) "i=",i
	          iErrMixIsof = 0
	          xIsof(1) = X1
	          xIsof(2) = 1D0 - X1
	          yIsof(1) = Y1
	          yIsof(2) = 1D0 - Y1
	          call isofug_mix(2,NG,NST,TEMP,PRESS,xIsof,yIsof,apartIsof,ZLL,ZVV,IERRMIXISOF)
	          APART1 = apartIsof(1)
	          APART2 = apartIsof(2)
!
	          IF(IERRMIXISOF.EQ.0) THEN

	            ZLCALC(NUMEXP+1)= ZLL
	            ZVCALC(NUMEXP+1)= ZVV
	            FCION(NUMEXP+1) = (W4**2.)*(2.+IFL(KS))*(APART2**2) + (APART1**2.)
	            FCION(NUMEXP+1)= DSQRT(FCION(NUMEXP+1))
	            FCAUX(NUMEXP+1)= APART1
	            FCAUX2(NUMEXP+1)= APART2

	          ELSE
!
!	            A failed point is not counted for the standard deviation
	            NFAILURE = NFAILURE - 1
	            ZLCALC(NUMEXP+1)= 0.
	            ZVCALC(NUMEXP+1)= 0.
	            FCION(NUMEXP+1)= PENALIZ
	            FCAUX(NUMEXP+1)= 0.
	            FCAUX2(NUMEXP+1)= 0.

	          END IF
	          STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)**2
	          STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)**2

	        end if
!
	        IF(IFL(KS).GE.0) then
!
	          IF(IFL(KS).EQ.0) THEN   !  ifl=0
!
	            write(6,*) " i = ",i
	            iErrTP = 0
!
!	            Binary multiphase flash. ZLIQUIDO always refers to phase 1, which is
!	            liquid. ZVAPOR can be vapor or liquid, depending on the result of the
!	            stability test. For LLE, the larger composition (X1 or Y1) decides how
!	            the flash result is ordered.
!
	            feed(1) = (X1 + Y1)/2
	            feed(2) = 1D0-feed(1)
	            indF = 0
	            zFrac = 0
	            zVal = 0
!
!	            Phase order (for LLE)
	            orden = Y1 - X1
	            if(orden.LT.0)then

	              iOrden = 1

	            else

	              iOrden = 2

	            endif
	            call GC3Flash(NCw,NG,NST,TEMP,PRESS,feed(:NCMgca),zFrac,zVal,indF,IERRTP)
!
!	            Phases in the same order as the experimental data:
	              if(zFrac(1,1).GT.zFrac(1,2))then

	                X1 = zFrac(1,iOrden)
	                zLiquido = zVal(iOrden)
	                Y1 = zFrac(1,3-iOrden)
	                zVapor = zVal(3-iOrden)

	              else

	                X1 = zFrac(1,3-iOrden)
	                zLiquido = zVal(3-iOrden)
	                Y1 = zFrac(1,iOrden)
	                zVapor = zVal(iOrden)

	              endif

	            IF(IERRTP.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQUIDO
	              ZVCALC(NUMEXP+1) = ZVAPOR
	              XCALC(NUMEXP+1) = X1
	              YCALC(NUMEXP+1) = Y1
	              FCAUX(NUMEXP+1) = ABS(XCALC(NUMEXP+1)-XEXP(NUMEXP+1))

	            ELSE

	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              XCALC(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF

	          END IF   ! end ifl = 0
!
	          IF(IFL(KS).GE.1.AND.iFL(ks).LE.2) THEN   !  ifl = 1 or 2: bubble point
!
	            write(6,*) "i = ",i
	            iErrBurb = 0
!
!	            iType for GCADB: IFL = 1 is a bubble pressure (iType = 3), IFL = 2 a
!	            bubble temperature (iType = 1).
	            if(IFL(KS).EQ.1)iType = 3
	            if(IFL(KS).EQ.2)iType = 1

	            call GCADB(1,NCw,NG,NST,iType,TEMP,PRESS,X1,Y1,iErrBurb,ZLiq,ZVap)
!
	            IF(IERRBurb.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQ
	              ZVCALC(NUMEXP+1) = ZVAP
	              YCALC(NUMEXP+1) = Y1
	              if(iType.EQ.3)then

	                PCALC(NUMEXP+1) = PRESS
!	                Absolute error
	                FCAUX(NUMEXP+1) = ABS(PCALC(NUMEXP+1) - PEXP(NUMEXP+1))

	              elseif(iType.EQ.1)then

	                Pcalc(numexp+1) = Temp ! the calculated bubble T is kept in PCALC
!	                Absolute error
	                FCAUX(NUMEXP+1) = ABS(PCALC(NUMEXP+1) - TEXP(NUMEXP+1))

	              endif

	            ELSE
!
	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              PCALC(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF
	          END IF   ! end ifl = 1, 2
!
!	          Ternary equilibrium at T, P (it could go together with the binary flash)
	          IF(IFL(KS).EQ.3) THEN   !  ifl=3
!
	            write(6,*) " i = ",i
	            iErrTP = 0
!
!	            Ternary multiphase flash
!
	            feed(1) = (xExp(numExp+1) + yExp(numExp+1))/2
	            feed(2) = (xExp2(numExp+1) + yExp2(numExp+1))/2
	            feed(3) = 1D0-feed(1)-feed(2)
!
!	            Phase order (for LLE)
	            orden = yExp(numExp+1) - xExp(numExp+1)
	            if(orden.LT.0)then

	              iOrden = 1

	            else

	              iOrden = 2

	            endif
	            call GC3Flash(NCw,NG,NST,TEMP,PRESS,feed(:NCMgca),zFrac,zVal,indF,IERRTP)
!
	            if(minval(indF(:4)).GT.-1)then

!	              LLE: phases in the same order as the experimental data
	              if(zFrac(1,1).GT.zFrac(1,2))then

	                X1 = zFrac(1,iOrden)
	                X2 = zFrac(2,iOrden)
	                zLiquido = zVal(iOrden)
	                Y1 = zFrac(1,3-iOrden)
	                Y2 = zFrac(2,3-iOrden)
	                zVapor = zVal(3-iOrden)

	              else

	                X1 = zFrac(1,3-iOrden)
	                X2 = zFrac(2,3-iOrden)
	                zLiquido = zVal(3-iOrden)
	                Y1 = zFrac(1,iOrden)
	                Y2 = zFrac(2,iOrden)
	                zVapor = zVal(iOrden)

	              endif

	            else

!	              VLE: minval(indF) = -1, since there must be a vapor phase;
!	              minLoc gives the position of the -1.
	              iVap(1) = minLoc(indF(:4),1)
	              Y1 = zFrac(1,iVap(1))
	              Y2 = zFrac(2,iVap(1))
	              X1 = zFrac(1,3-iVap(1))
	              X2 = zFrac(2,3-iVap(1))
	              zLiquido = zVal(3-iVap(1))
	              zVapor = zVal(iVap(1))

	            endif
	            IF(IERRTP.EQ.0) THEN

	              ZLCALC(NUMEXP+1) = ZLIQUIDO
	              ZVCALC(NUMEXP+1) = ZVAPOR
	              XCALC(NUMEXP+1) = X1
	              XCALC2(NUMEXP+1) = X2
	              YCALC(NUMEXP+1) = Y1
	              YCALC2(NUMEXP+1) = Y2
	              FCAUX(NUMEXP+1) = ABS(XCALC(NUMEXP+1)-XEXP(NUMEXP+1))

	            ELSE

	              NFAILURE = NFAILURE - 1
	              ZLCALC(NUMEXP+1) = 0.
	              ZVCALC(NUMEXP+1) = 0.
	              XCALC(NUMEXP+1) = 0.
	              XCALC2(NUMEXP+1) = 0.
	              YCALC(NUMEXP+1) = 0.
	              YCALC2(NUMEXP+1) = 0.
	              FCAUX(NUMEXP+1) = 0.

	            END IF

	          END IF   ! end ifl = 3
!
!	          Residuals. Errors are absolute, divided by the experimental
!	          uncertainty through the weights wBin and wTer (Francisco, 07/08/2013).
	          IF(IFL(KS).EQ.0.AND.IERRTP.EQ.0)then

	            if(yExp(numExp+1) .GT. 0D0 .AND. xExp(numExp+1) .GT. 0D0)then

	              FCION(NUMEXP+1) = SQRT(w4**2*(YCALC(NUMEXP+1) - YEXP(NUMEXP+1))**2*wBin(2,numExp+1)**2 + &
&	                                     w5**2*(XCALC(NUMEXP+1) - XEXP(NUMEXP+1))**2*wBin(1,numExp+1)**2)
	              FCAUX2(NUMEXP+1) = ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	            endif

	          elseif(IFL(KS).EQ.1.AND.IERRBurb.EQ.0.)THEN
!
!	            Bubble pressure and y
	            if(yExp(numExp+1) .GT. 0D0)then
	              FCION(NUMEXP+1) = SQRT((W4**2)*(YCALC(NUMEXP+1) - YEXP(NUMEXP+1))**2*wBin(2,numExp+1)**2 + &
&	                                w2**2*fcAux(numexp+1)**2*wBin(1,numExp+1)**2)
	              FCAUX2(NUMEXP+1) = ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	            endif

	          elseif(IFL(KS).EQ.2.AND.IERRBurb.EQ.0.)THEN
!	            Bubble temperature and y
	            if(yExp(numExp+1) .GT. 0D0)then
	              FCION(NUMEXP+1) = SQRT((W4**2)*(YCALC(NUMEXP+1) - YEXP(NUMEXP+1))**2*wBin(2,numExp+1)**2 + &
&	                                w3**2*fcAux(numexp+1)**2*wBin(1,numExp+1)**2)
	              FCAUX2(NUMEXP+1) = ABS(YCALC(NUMEXP+1)-YEXP(NUMEXP+1))

	            endif

	          elseif(iFL(KS).EQ.3.AND.iErrTP.EQ.0)then

	            if(yExp(numExp+1) .GT. 0D0)then

	              fcAux2(numExp+1) = abs(yCalc(numExp+1) - yExp(numExp+1))
	              fcion(numExp+1) = SQRT(w5**2*(FCAUX(NUMEXP+1))**2*wBin(1,numExp+1)**2 &
&	                              + w6**2*(xCalc2(numExp+1) - xExp2(numExp+1))**2*wTer(1,numExp+1)**2 &
&	                              + w4**2*(fcAux2(numExp+1))**2*wBin(2,numExp+1)**2 &
&	                              + w7**2*(yCalc2(numExp+1) - yExp2(numExp+1))**2*wTer(2,numExp+1)**2)
	            endif

	          ELSE

	            FCION(NUMEXP+1)= PENALIZ*wBin(1,numExp+1)**2
	            FCAUX2(NUMEXP+1) = 0.D0

	          END IF
	          STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)
	          STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)
	          if(mez(ks).EQ.3)then

	            STDVC(KS) = STDVC(KS) + ABS(XCALC2(NUMEXP+1) - XEXP2(NUMEXP+1))
	            STDVD(KS) = STDVD(KS) + abs(yCalc2(numExp+1) - yExp2(numExp+1))

	          endif
!
	        end if
!
	        NUMEXP = NUMEXP + 1

	      END DO

!           Standard deviations
	      IF(NFAILURE.NE.0) THEN

	        STDVA(KS) = STDVA(KS)/NFAILURE
	        STDVB(KS) = STDVB(KS)/NFAILURE
	        if(mez(ks).EQ.3)then

	          stDvC(ks) = stDvC(ks)/nFailure
	          stDvD(ks) = stDvD(ks)/nFailure

	        endif

	      ELSE

	        STDVA(KS) = 0.D0
	        STDVB(KS) = 0.D0

	      END IF
!
!
!---------Excess property data------------------------------------------------------------------------------------
	    elseIF(iHen(ks).EQ.1)then

	      nFailure= 0
	      iErrEx = 0
	      stdva(ks) = 0D0
	      write(*,'(/,"Set (excess property) No. ",I2)')ks
	      do i =1,nObs(ks)

	        Temp = Texp(numExp+1)
	        Press = Pexp(numexp+1)
	        feed(1) = xExp(numexp+1)
	        feed(2) = 1D0 - feed(1)
	        if(feed(2).LT.0.D0)feed(2) = 1D-20
	        write(*,'("  Point No. ",I2)')i
	        call GCXS(ifl(ks),NG,NST,Temp,Press,feed(:NCMgca),prop_E_calc(numExp+1),zLiCalc(numexp+1),iErrEx)
	        if(iErrEx.EQ.0)then

	          if(prop_E_exp(numExp+1).NE.0D0)then

!	            Absolute error
	            fcAux(numExp+1) = dabs(prop_E_calc(numExp+1) - prop_E_exp(numExp+1))

	          else

	            fcAux(numExp+1) = 0

	          endif
	          if(ifl(ks).EQ.1)then

	            fcion(numExp+1) = w11*fcAux(numExp+1)*wEx(numExp+1)

	          else

	            fcion(numExp+1) = w12*fcAux(numExp+1)*wEx(numExp+1)

	          endif

	        else

	          nFailure = nFailure + 1
	          zLiCalc = 0.
	          prop_E_calc(numExp) = 0.
	          fcAux(numExp+1) = 0.
	          fcion(numExp+1)= penaliz*wEx(numExp+1)

	        endif
	        NUMEXP = NUMEXP + 1
	        if(i.EQ.1)numExp1 = numExp
	        nObsKS = nObs(ks)
	        if(i.EQ.nObsKS)numExp2 = numExp

	      enddo
	      if(nObs(ks)-nFailure.GT.0)stDvA(ks) = sqrt(sum(fcion(numExp1:numExp2)**2)/(nObs(ks)-nFailure))
	      
!  
!
!
!---------SLE data-----------------------------------------------------------------------------------------------
	    elseif(iHen(ks).EQ.2)then
	
	      nFailure = 0

	      StDva(ks) = 0
	      write(*,'(/,"Set (SLE) No. ",I2)')ks
	      do i =1,nObs(ks)

	        Temp = Texp(numExp+1)
	        Press = Pexp(numexp+1)
	        feed(1) = xExp(numexp+1)
	        feed(2) = 1D0 - feed(1)
	        iErrSol = 0
!	        zFrac(:,1) gets the liquid composition, zFrac(:,2) the fugacity (kept
!	        only because the original program printed it).
	        call GCSLE (NCw, NG, NST, Temp, Press, feed(:NCMgca), wwDHfu, wwDVfu, wwTfu, wwPfu, wwDCpfu, wwTtrSS, wwDHss, &
	                    CZ, zFrac(:NCMgca,1), zFrac(:NCMgca,2),iErrSol)	  
	        if(iErrSol.EQ.0)then
	  
	          zLiCalc(numExp+1) = cZ
	          xCalc(numExp+1) = zFrac(1,1)
	          yCalc(numExp+1) = zFrac(1,2)
!	          Absolute error
	          fcAux(numExp+1) = dabs(zFrac(1,1) - feed(1))
	          fcion(numExp+1) = wEx(numExp+1)*fcAux(numExp+1)*w13
	        
	        else
	        
	          zLiCalc(numExp+1) = 0
	          xCalc(numExp+1) = 0
	          yCalc(numExp+1) = 0
	          nFailure = nFailure + 1
	          fcion(numExp+1)= penaliz*wEx(numExp+1)
	          fcAux(numExp+1) = 0
	        
	        endif
	        NUMEXP = NUMEXP + 1
	        if(i.EQ.1)numExp1 = numExp
	        nObsKS = nObs(ks)
	        if(i.EQ.nObsKS)numExp2 = numExp	      
	    
	      enddo
	      if(nObs(ks)-nFailure.GT.0)stDvA(ks) = sum(fcAux(numExp1:numExp2))/(nObsKS - nFailure)
!  
!
!
!---------Critical point data--------------------------------------------------------------------------------------
	    elseif (iHen(ks) == 3) then
	    
	      nFailure = nObs(ks)
	      StDva(ks) = 0D0
	      write(*,'("Critical point data ", I2)')ks
	      do i = 1, nObs(ks)
	      
	        Temp = Texp(numExp+1)
	        Press = Pexp(numexp+1)
	        feed(1) = xExp(numexp+1)
	        feed(2) = 1D0 - feed(1)
	        iErrPC = 0
	        
	        call implicitPCrit(NCw, NG, NST, Temp, Press, feed(:NCMgca), dlnF_dz1, d2lnF_dz12, zLiCalc(numExp+1))
	        fcion(numExp+1) = wBin(1,numExp+1)**2*dlnF_dz1**2 + wBin(2,numExp+1)**2*d2lnF_dz12**2
	        fcAux(numExp+1) = dlnF_dz1
	        fcAux2(numExp+1) = d2lnF_dz12	      
	        NUMEXP = NUMEXP + 1
	        
	      enddo
!  
!
!
!---------Activity coefficient data--------------------------------------------------------------------------------
	    ELSEIF (iHen(ks).EQ.5) then

	      nfailure=nobs(ks)
	      stdva(ks) = 0D0
	      stdvb(ks) = 0D0
	      stdvc(ks) = 0D0
	      WRITE(6,'("Activity coefficient data ",I2)')KS
	      do i=1,nobs(ks)

	        Temp  = texp(numexp+1)
	        Press = pexp(numexp+1)
	        feed(1) = xExp(numExp+1)
	        if(mez(ks) == 2)then
	        
	          feed(2) = 1.D0 - feed(1)
	          
	        else
	        
	          feed(2) = xExp2(numExp+1)
	          feed(3) = 1.D0 - feed(1) - feed(2)
	          
	        endif
	        do io = 1, mez(ks)
	          if(feed(io) <= 2.E-16)feed(io) = 2E-16
	        enddo
	        WRITE(6,*)"OBS=", NOBS(KS)
	        iaux=mez(ks)
	        iErrGam = 0
	        iTyp = 1; iGz = 0; jTemp = 0
! 	        zFrac gets ln(gamma)
	        call paragc(Temp, NCw, NG, NST, jTemp)
	        call GAMMA(NCw, NG, NST, iTyp, iGz, jTemp, feed(:NCMgca), Temp, Press, zFrac(:NCMgca,1), zFrac(:NCMgca,2), zLiCalc(numExp+1))
	        IF(iErrGam==0)THEN

	          rgamcalc1(numexp+1) = dexp(zFrac(1,1))
	          rgamcalc2(numexp+1) = dexp(zFrac(2,1))
!	          Relative errors (for the statistics)
	          fcaux(numexp+1)=rgamcalc1(numexp+1)/rgamexp1(numexp+1)-1.D0
	          fcaux2(numexp+1)=rgamcalc2(numexp+1)/rgamexp2(numexp+1)-1.D0
!	          Absolute errors (for the residual)
	          fcion(numexp+1) = w8**2*(rgamcalc1(numexp+1) - rgamexp1(numexp+1))**2*wGam(1,numExp+1)**2 + &
&	                            w9**2*(rgamcalc2(numexp+1) - rgamexp2(numexp+1))**2*wGam(2,numExp+1)**2
	          If(mez(ks) == 3)then

	            rGamCalc3(numexp+1) = dexp(zFrac(3,1))
	            fcaux3(numexp+1)=rgamcalc3(numexp+1)/rgamexp3(numexp+1)-1.d0
	            fcion(numexp+1) = fcion(numexp+1) + w10**2*(rGamcalc3(numexp+1) - rGamexp3(numexp+1))**2*wGam(3,numExp+1)**2

	          end if
!
	          fcion(numexp+1) = SQRT(fcion(numexp+1))

	        else

	          nfailure = nfailure - 1
	          zlicalc(numexp+1) = 0.
	          rgamcalc1(numexp+1) = 0.
	          rgamcalc2(numexp+1) = 0.
	          rgamcalc3(numexp+1) = 0.
	          fcaux(numexp+1) = 0.
	          fcaux2(numexp+1) = 0.
	          fcaux3(numexp+1) = 0.
	          fcion(numexp+1)= penaliz*wGam(1,numExp+1)

	        end if
	        STDVA(KS) = STDVA(KS) + FCAUX(NUMEXP+1)**2
	        STDVB(KS) = STDVB(KS) + FCAUX2(NUMEXP+1)**2
	        STDVC(KS) = STDVC(KS) + FCAUX3(NUMEXP+1)**2
	        numexp=numexp+1

	      end do
	      IF(NFAILURE.NE.0) THEN

	        STDVA(KS) = DSQRT(STDVA(KS)/NFAILURE)
	        STDVB(KS) = DSQRT(STDVB(KS)/NFAILURE)
	        STDVC(KS) = DSQRT(STDVC(KS)/NFAILURE)

	      ELSE

	        STDVA(KS) = 0.D0
	        STDVB(KS) = 0.0D0
	        STDVC(KS) = 0.

		END IF


	    END IF

	  END DO
!**************************** End of the data sets ******************************

	END IF
!
!	NaN guard (added 2026-11-06): a NaN or Inf FCION(i) is replaced by PENALIZNAN.
	NCOUNTNAN = 0
	DO I=1,NCXP
	  IF (FCION(I).NE.FCION(I) .OR. DABS(FCION(I)).GT.1.D30) THEN
	    FCION(I) = PENALIZNAN
	    NCOUNTNAN = NCOUNTNAN + 1
	  END IF
	END DO
	IF (NCOUNTNAN.GT.0) WRITE(*,'(A,I6)') ' WARNING: non-finite FCION at points: ',NCOUNTNAN
!
	SUMQ = 0D0
	DO I=1,NCXP

	  SUMQ = SUMQ + FCION(I)**2

	END DO
	sumQ = sum(fcion(:NCXP)**2)
!
	OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
!
!	Course of the optimization: unscaled parameters and objective function
	if(nPar.GT.0)then
	  CALL XXUNSCAL(NPAR,XX,IDXX,XXAUX)
	  WRITE(1,'('//rep(NPAR,"F12.4,X")//'E12.5)') (XXAUX(I),I=1,NPAR),SUMQ
	endif

	CLOSE(UNIT=1)
	RETURN

ENDsubroutine FUNC
!------------------------------------------------------------------------------------------------------------------------------------
!
!	Estimates dc from the critical point of each pure component (only for
!	the components whose dc was not given).
!
SUBROUTINE init_dc(NC,PRCRIT,TEMPCRIT,DC)
	use dimensions

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION DC(NCM),TEMPCRIT(NCM),PRCRIT(NCM)
!
	R=82.05	!atm cm3/(mol K)
	DO I=1,NC
	  if(DC(I).EQ.0)then ! Fran. 26/05/09
	    DC(I) = 0.4471866767D0*((R*TEMPCRIT(I)/PRCRIT(I)))**(1./3.)
	  end if
	END DO
!
	return
ENDsubroutine init_dc
!------------------------------------------------------------------------------------------------------------------------------------
!
!	Sets the type of each fitted parameter (IDXX, see module scaling), in
!	the order they have in XX.
!
SUBROUTINE set_param_types(NNP,IDG,ITEMP,NNQ,NALFA,NENASS,NRKASS,IDXX,numTc,NUMDC)
	use dimensions
	use scaling, only: inIJ
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
!	Indices of the gii or kij parameters to fit.
	DIMENSION IDG(NPTIPOM,2),ITEMP(NPTIPOM),IDXX(NPARM),xx(NPARM)
!
	inIJ(:) = .FALSE.
	NPAR = 0
	IF(NNP.NE.0) THEN

	  DO I=1,NNP

	    IF(IDG(I,1).EQ.IDG(I,2)) THEN

	      IDXX(NPAR+1)=1
	      NPAR = NPAR + 1

	    ELSE

	      IDXX(NPAR+1)=2
	      NPAR = NPAR + 1
	      inIJ(nPar) = .TRUE.

	    END IF
	    IF(ITEMP(I).EQ.1)then

	      IDXX(NPAR+1) = 3
	      NPAR = NPAR + 1
	      IF(IDG(I,1).NE.IDG(I,2))inIJ(nPar) = .TRUE.
!
	    elseIF(ITEMP(I).EQ.2)then

	      IDXX(NPAR+1) = 3
!
!	      g' and g" have different scaling types:
	      idXX(nPar+2) = 10
	      NPAR = NPAR + 2

	    endif
!
	  END DO

	END IF
	IF(NNQ.NE.0) THEN

	  DO I=1, NNQ

	    IDXX(NPAR+1)= 4
	    NPAR = NPAR + 1

	  eND DO

	END IF
	IF(NALFA.NE.0) THEN

	  DO I=1, NALFA

	    IDXX(NPAR+1)= 5
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(NENASS.NE.0) THEN

	  DO I=1,NENASS

	    IDXX(NPAR+1)= 6
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(NRKASS.NE.0) THEN

	  DO I=1,NRKASS

	    IDXX(NPAR+1)= 7
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(numTc.NE.0) THEN

	  DO I=1,numtc

	    IDXX(NPAR+1)= 8
	    NPAR = NPAR + 1

	  END DO

	END IF
	IF(NUMDC.NE.0) THEN

	  DO I=1,NUMDC

	    IDXX(NPAR+1)= 9
	    NPAR = NPAR + 1

	  END DO

	END IF
	return
ENDsubroutine set_param_types
!
!------------------------------------------------------------------------------------------------------------------------------------
!
!	Writes the statistics of the fit to paraest.out.
!
SUBROUTINE write_fit_stats(SSQ,XJAC,WORK,INFER,IER)
	use dimensions
!
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XJAC(NDM,NPARM), WORK(LWORK)
!
	OPEN(UNIT=1,STATUS='OLD',file='paraest.out',FORM='FORMATTED',POSITION='APPEND')
!
	WRITE(1,95)
95	FORMAT(/,X,' Final value of the objective function = SSQ = ',/)
	WRITE(1,*) SSQ
!
	WRITE(1,97) IER,INFER
 97	FORMAT(/,' (See the manual of routine ZXSSQ (IMSL) ',/,X,' IER = ',I5,/,X,' INFER = ',I5)
!
	WRITE(1,99) WORK(1),WORK(2),WORK(3),WORK(5)
99	FORMAT(/,'  Gradient norm = ',D12.5,/,X,' Number of objective function evaluations = ',F10.2,/, &
&	       X,' Number of significant digits in the scaled parameters = ',F10.2,/,X,' Number of iterations = ',F10.2)
	CLOSE(UNIT=1)
	return

END SUBROUTINE write_fit_stats
!------------------------------------------------------------------------------------------------------------------------------------
!
!	Writes the final parameters, and the experimental data against the
!	GCA-EoS results (with the given parameters or with the fitted ones), to
!	paraest.out, tables.csv and graficos.csv.
!
SUBROUTINE write_results(NTEXT,stDvD)
	use dimensions
	use fit_data, only: NC, ISAT, IGLE, IVE, IUNI, NSAT, IDC, NOBS, IHEN, mez, IFL, TEXP, PEXP, XEXP, xexp2, &
&	                    YEXP, yExp2, IDSAT, tempcrit, XCALC, YCALC, xCalc2, yCalc2, ZLCALC, ZVCalc, FCAUX, &
&	                    PCALC, IOBJSAT, FCAUX2, fcaux3, STDVA, STDVB, stdvc, rgamcalc1, rgamcalc2, rgamcalc3, &
&	                    rgamexp1, rgamexp2, rgamexp3, zLicalc, prop_E_calc, prop_E_exp, DHVL
	use gca_parameters, only: NG, NU, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, NGA, ONUASOC, MASSOC, &
&	                    ENASS, RKASS, TSAT, PSAT, DC, numTc, idTc, NUMDC, IDDC, TSPLINE
	use utilities, only: i2s, rep
	use statistics, only: AAD, ARD, ASD, RSD
!
	IMPLICIT REAL*8 (A-H,O-Z)
!
	character(1) :: tab = achar(9)
	character(2) :: binary_number
	DIMENSION NTEXT(NSetM,65)        ! title of each data set
	DIMENSION STDVSAT(NCM)           ! standard deviation of the pure components
	DIMENSION stDvD(NSetM)           ! standard deviations of the mixtures
!
	OPEN (UNIT = 1, STATUS='OLD', file = 'paraest.out', FORM = 'FORMATTED', POSITION='APPEND')
	open (unit = 2, file = "tables.csv", form = "formatted")
	open (unit = 10, file = "graficos.csv", form = "FORMATTED")	
!
!
	WRITE(1,47)
 47	FORMAT(/,X,'--------------------------------------------------',/, &
&	         X,' ---------------- FINAL PARAMETERS ----------------',/,   &
&	         X,'--------------------------------------------------',/,   &
&	        X,' PURE-GROUP PARAMETERS (dispersive term)',/)
!
	WRITE(1,"(40X,'COMPONENTS',/,X,' Gr.No.:   T*        Q        TSPLINE ',"//rep(NC,'I3')//"/)") (I, I=1,NC)
!
	DO I=1,NG

	  WRITE(1,4) I,TSTR(I),QAREA(I),TSPLINE(I), (NU(I,J), J=1,NC)
4	  FORMAT(X,I3,5X,2F10.4,X,E10.4,*(I3))

	END DO
!
	WRITE(1,5)
5	FORMAT(/,X,' Interaction energies gjj ',/,X,'(the functional form changes at high temperatures)',/)
!
	DO I=1,NG
!
	  WRITE(1,6) I,I,GSTR(I), GDOT(I), GDDOT(I)
6	  FORMAT(X,' g(',I3,I3,') = ',F10.1,' * [ 1 + ',F10.4,' *(T/T* - 1) + ',F10.4, ' *LN(T/T*)] ',/)

	END DO

	WRITE(1,7)
 7	FORMAT(/,X,' Binary parameters: [k*] ',/)

	DO I=1,NG

	  WRITE(1,8) (RKIJ(I,J), J=1,NG)
8	  FORMAT(X,*(F10.4))

	END DO
!
	WRITE(1,9)
9	FORMAT(/,X,' Temperature dependence: [k`]',/)
	DO I=1,NG

	  WRITE(1,8) (RKTIJ(I,J), J=1,NG)

	END DO
!
	WRITE(1,10)
 10	FORMAT(/,X,' Binary non-randomness parameters: [ALPHA]',/)
	DO I=1,NG

	  WRITE(1,8) (ALFAIJ(I,J), J=1,NG)

	END DO
	WRITE(1,186) NGA
 186	FORMAT(/,X,' Number of groups in the association contribution = ',I3)
!
	IF(NGA.GT.0) THEN

	  WRITE(1,"(/,X,' Association parameters ',/,23X,'Components',/,' Group 	  # Assoc sites',"//rep(NC,'I3')//"/)") (I, I=1,NC)

	  DO I=1,NGA

	    WRITE(1,12) I, MASSOC(I), (ONUASOC(I,J), J=1,NC)
12	    FORMAT(3X,I3,9X,I3,6X,*(f5.2))

	  END DO
	  WRITE(1,13)
13 	  FORMAT(/,X,'Association energy parameters (K)  ',/)
!
!	  Association energies
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then

	      DO L=1,NGA

              IF(MASSOC(L).NE.0)then

	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,14) I,J,K,L,ENASS(I,J,K,L)
14	              FORMAT(X,' ENASS(',I2,I2,X,I2,I2,')  = ',F10.4)

	            END DO

	          END DO

	        endif

	      END DO

	    endif

	  END DO
!
!	  Association volume parameters (cm3/mol)
	  WRITE(1,15)
15	  FORMAT(/,X,'Association volume parameters (cm3/mol) ',/)
	  DO J=1,NGA

	    IF(MASSOC(J).NE.0)then

	      DO L=1,NGA

	        IF(MASSOC(L).NE.0)then

	          DO I=1,MASSOC(J)

	            DO K=1,MASSOC(L)

	              WRITE(1,16) I,J,K,L,RKASS(I,J,K,L)
16	              FORMAT(X,' RKASS(',I2,I2,X,I2,I2,') = ',F16.10)

	            END DO

	          END DO

	        endif

	      END DO

	    endif

	  END DO

	endif
!
!
!	Result of the dc fit
	DO I=1,NC

	  DO JAMO=1,NUMDC

	    IF(IDDC(JAMO).EQ.I)  THEN

	      WRITE(1,1548) I,DC(I)
1548	      FORMAT(/,X,' The hard-sphere diameter at the critical point of component No.',I3,/,&
&	             ' was a fitted parameter ',/,X,' DC = ', F10.7)
	      GO TO 505

	    END IF

	  END DO
!
	  IF(TSAT(I).NE.0) WRITE(1,500) I,TSAT(I),PSAT(I),DC(I)
500	  FORMAT(/,X,' The hard-sphere diameter at the critical point of component No.',I3,' was fitted to ',/, &
&	         X,' reproduce the point of the vapor pressure curve ',/,x,' [ ',F10.4,' K ;',G10.4,' atm ] ',/,X,' DC = ', F10.7)
505	END DO
 	DO I=1,NC
!
	  DO mm=1,numtc

	    IF(idtc(mm).EQ.I)  THEN

	      WRITE(1,1549) I,tempcrit(I)
1549	      FORMAT(/,X,'The critical temperature of component No.',I3,/,X,'was a fitted parameter ',/,'Tc = ', F8.2," K")

	    END IF

	  END DO

	END DO
!
!
!	**************** Pure-component results *****************
!
!	Besides paraest.out, the results go to graficos.csv (unit 10), only to make
!	plots easily, and a summary to tables.csv (unit 2).
	NCXP = 0
 2	FORMAT(*(I3))
!
	IF(ISAT.NE.0) THEN

	  WRITE (1, 338)
	  write (2, 338)
	  write (10, 338)
338	  FORMAT(/,X,'***** Results: pure components *****',/)
	  write (2, '("Compound,Tc (K),dc (cm/mol^1/3), DTr, ARD(P)%,Reference")')
	  DO I=1,NC

	    IF(NSAT(I).NE.0) THEN

	      IF(IOBJSAT(I).EQ.50) THEN

	        WRITE(1,393)
	        write(10,393)
393	        FORMAT(/,X,'  Comp   Texp(K)    Pexp(atm)       Pexp(atm)        Zliq       Zvap  REL,ERR,%,ISOF',/)

	      END IF
	      IF(IOBJSAT(I).EQ.80) THEN

	        WRITE(1,333)
	        write(10,3331)
333	        FORMAT(/,X,'  Comp   Texp(K)    Pexp(atm)    Pcalc(atm)     Zliq      Zvap        DP%     DH_LV(J/mol) ',/)
3331	        FORMAT(/,X,'Comp,Texp(K),Pexp(atm),Pcalc(atm),Zliq,Zvap,DP%,DH_LV(J/mol)',/)

	      END IF
!
	      NFOLI=NSAT(I)
	      STDVSAT(I)=0D0
	      DO J=1,NSAT(I)

	        WRITE(1,334) I,TEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100, DHVL(NCXP+1)*.101325D0
	        WRITE(10,3341) I,TEXP(NCXP+1),PEXP(NCXP+1),PCALC(NCXP+1), ZLCALC(NCXP+1), ZVCALC(NCXP+1),FCAUX(NCXP+1)*100, DHVL(NCXP+1)*.101325D0
334	        FORMAT(X,I3,4X,F10.4,X,G12.5,X,G12.5,X,G11.4,X,G11.4,X,G10.3,G11.4)
3341	        FORMAT(I3,",",F10.4,",",G12.5,",",G12.5,",",G11.4,",",G11.4,",",G10.3,",",G11.4)
	        NCXP = NCXP + 1

	      END DO
!
!	      Average absolute and relative errors and standard deviations
	      N1 = NCXP-nSat(i)+1
	      E1 = AAD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i)); E2 = 100*ARD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i))
	      E3 = ASD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i)); E4 = 100*RSD(Pexp(N1:NCXP),Pcalc(N1:NCXP),nSat(i))
	      write(*,335) E1, E2
	      write(1,335) E1, E2, E3, E4
	      write (2, '(I2, ",", F6.2, ",", F7.4, ",", F4.2, "-", F4.2, ",", F7.3  )') i, TempCrit(i), dc(i), Texp(N1)/TempCrit(i),Texp(NCXP)/TempCrit(i), E2
	      write(10,335) E1, E2, E3, E4
335	      format(/" AAD(ARD%) = ,",F10.5,", atm , (, ",F7.2,"% ,)", /, " AbsStdDev(Rel.%) = ,",F10.5,", atm , (, ",F7.2,"% ,)"/)

	    END IF

	  END DO

	END IF
!
!	Mixture results
	IF(IGLE.NE.0) THEN
!
	  WRITE(1,340)
	  write (2, 340)
	  write(10,340)
340	  FORMAT(/,X,'********** RESULTS: BINARY MIXTURES *******************',/)
	  write (2, '("Compound,-,T (K),P (kPa),AARD%,-,No exp.,Reference")')
	  write (2, '("1,2,-,-,P,y1,points,")')
	  DO KS=1,IGLE

	    WRITE(1,339) (NTEXT(KS,J),J=1,65)
	    WRITE(10,339) (NTEXT(KS,J),J=1,65)
339	    FORMAT(//,80A2)
!
!	    Width of the number in parentheses when printing x(i): with blanks, Excel
!	    does not split the column headers correctly.
	    if(idC(ks,1).LT.10)then

	      iDim = 1

	    else

	      iDim = 2

	    endif
	    if(idC(ks,2).LT.10)then

	      iDim2 = 1

	    else

	      iDim2 = 2

	    endif
	    IF(IHEN(KS).EQ.0) THEN

!	      ******************** Binary equilibrium ***********************
!
	      IF(IFL(KS).EQ.-2) WRITE(1,343) IDC(KS,2)
	      IF(IFL(KS).EQ.-2) WRITE(10,343) IDC(KS,2)
343	      FORMAT(X,/,' NOTE: the departure from isofugacity of component ',I3,/, &
&	             X,' was not used in the objective function, but',/,' it is reported anyway ',/)
!
          IF(IFL(KS).EQ.-1.OR.IFL(KS).EQ.-2) then

	        WRITE(1,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2)
	        WRITE(10,341) IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2)
341	        FORMAT(X,'  Temp/K    Pres/ATM  Xexp(',I0,')  Yexp(',I0,') Fug(',I0,')% Fug(',I0,')%   Zliq      Zvap ',/)
	        DO I=1,NOBS(KS)

	          WRITE(1,345) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),YEXP(NCXP+1),FCAUX(NCXP+1),FCAUX2(NCXP+1),ZLCALC(NCXP+1), &
&	          ZVCALC(NCXP+1)
	          WRITE(10,345) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),YEXP(NCXP+1),FCAUX(NCXP+1),FCAUX2(NCXP+1),ZLCALC(NCXP+1),&
&	          ZVCALC(NCXP+1)
345	          FORMAT(6F10.5,2F9.6)
	          NCXP = NCXP + 1

	        END DO
	        WRITE(1,347) STDVA(KS),STDVB(KS)
	        WRITE(*,347) STDVA(KS),STDVB(KS)
	        WRITE(10,347) STDVA(KS),STDVB(KS)
347	        FORMAT(/,X,' Std.Dev. %    ',24X,2f10.4)

	      elseIF(IFL(KS).EQ.0) then

!	        Flash at T and P
	        WRITE(1,351)tab,tab,tab,IDC(KS,1),tab,IDC(KS,1),tab,tab,IDC(KS,1),tab,IDC(KS,1),tab,tab,tab
351	        FORMAT(/,A,'T/K   ',A,'P   /atm',A,'Xexp(',I0,')',A,'Xclc(',I0,')',A,'Err%',A,'Yexp(',I0,')',A,  &
&	                'Ycalc(',I0,')',A,'Err%',A,'Z1',A,'Z2',/)
!
	        WRITE(10,3511) IDC(KS,1),IDC(KS,1)

3511	        FORMAT(/,X," (x-y)(",I0,")exp,(x-y)(",I0,")calc,T/K,P/atm")
!
!	        Second counter, for graficos.csv
	        NCXPy = NCXP
	        DO I=1,NOBS(KS)
!
	          WRITE(1,349) tab,TEXP(NCXP+1),tab,PEXP(NCXP+1),tab,XEXP(NCXP+1),tab,XCALC(NCXP+1),tab, FCAUX(NCXP+1),tab,YEXP(NCXP+1),tab,YCALC(NCXP+1),tab,FCAUX2(NCXP+1),tab, &
&	                       ZLCALC(NCXP+1),tab, ZVCALC(NCXP+1)
349	          FORMAT(A,2(F10.5,A),6(G10.4,A),G10.4,A,F10.8)
!
!	          x in graficos.csv
	          write(10,3490)XEXP(NCXP+1),XCALC(NCXP+1),TEXP(NCXP+1),PEXP(NCXP+1)
	          write(11,'(F7.2,3(",",G10.4),A,F7.2,3(",",G10.4))')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),yExp(NCXP+1),",|,",Texp(NCXP+1),Pexp(NCXP+1),xCalc(NCXP+1),yCalc(NCXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
!	        y in graficos.csv
	        do i = 1,Nobs(ks)

	          write(10,3490)YEXP(NCXPy+1),YCALC(NCXPy+1),TEXP(NCXPy+1),PEXP(NCXPy+1)
	          NCXPy = NCXPy + 1

	        enddo
3490	        format(2(G10.4,","),F7.2,",",F10.4)
	        NO = nObs(ks)
	        N1 = NCXP - NO + 1
	        E1 = AAD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)
	        E2 = ARD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)*100
	        E3 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
	        E4 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
	        write(1,353)idC(ks,1),E1,E2, idC(ks,1),E3,E4
	        write(*,'(/,X," ARD(x) =",f10.4,"%",20X,"ARD(y) = ",f10.4,"%")')E2,E4
	        write(10,353)idC(ks,1),E1,E2, idC(ks,1),E3,E4
353	        Format(/" AAD(ARD%): x(",I0,") =, ",E10.4,", (, ",F6.2,"% ,)  ,",20X," y(",I0,") = ,",E10.4,", (, ",F6.2,"% ,)")


	      elseIF(IFL(KS).EQ.1) then

!	        Bubble pressure (T and x)
	        NCXPy = NCXP
	        WRITE(1,355)tab,tab,IDC(KS,1),tab,tab,tab,tab,IDC(KS,1),tab,IDC(KS,1),tab,tab,tab
355	        FORMAT(/,A,'T/K',A,'Xexp(',I0,')',A,'Pexp/atm',A,'Pcalc/atm',A,'ERR%',A,'YEXP(',I0,')',A,'YCALC(',I0,')', &
&	                 A,'ERR',A,'Z1',A,'Z2'/)
	        write(10,3550) IDC(KS,1), IDC(KS,1)
3550	        format(/,"T /K,(x-y)(",I0,")exp,Pexp /atm,(xexp-ycalc)(",I0,"),Pcalc /atm")
	        DO I = 1, NOBS(KS)

	          WRITE (1, 357) tab,TEXP(NCXP+1),tab,XEXP(NCXP+1),tab,PEXP(NCXP+1),tab,PCALC(NCXP+1),tab,FCAUX(NCXP+1),tab,YEXP(NCXP+1), &
&	                         tab,YCALC(NCXP+1),tab,FCAUX2(NCXP+1),tab,ZLCALC(NCXP+1),tab,ZVCALC(NCXP+1)
357	          FORMAT(A,F10.2,A,G10.4,A,G10.4,A,G10.4,A,F10.5,A,G10.5,A,F10.5,A,F10.5,A,F10.8,A,F10.8)
!
!	          x in graficos.csv
	          write(10,3570)TEXP(NCXP+1),XEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),PCALC(NCXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
	        do i = 1,NObs(ks)

	          write(10,3570)TEXP(NCXPy+1),yEXP(NCXPy+1),PEXP(NCXPy+1),ycalc(NCXPy+1),PCALC(NCXPy+1)
	          NCXPy = NCXPy + 1

	        enddo
	        NO = nObs(ks)
	        N1 = NCXP - NO + 1
	        E1 = AAD(PExp(N1:NCXP), PCalc(N1:NCXP), NO)
	        E2 = ARD(PExp(N1:NCXP), PCalc(N1:NCXP), NO)*100
	        E3 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
	        E4 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
	        write(1,359)E1,E2,idC(ks,1),E3,E4
	        write(*,'(/,X," ARD(P) =",f10.4,"%",20X,"ARD(y) = ",f10.4,"%")')E2,E4
	        write(10,359)E1,E2,idC(ks,1),E3,E4
359	        Format(/" AAD(ARD%): P = ,",E10.4,", atm (, ",F6.2,"% ,)  ,",20X," y(",I0,") = ,",E10.4," (, ",F6.2,"% ,)")
3570	        format(F7.2,",",F9.6,",",F10.4,",",F9.6,",",F10.4)
!	        Summary line in tables.csv (decimals chosen by magnitude)
	        Pmin = 101.325*minval(Pexp(N1:NCXP)); Pmax = 101.325*maxval(Pexp(N1:NCXP))
	        if (Pmin >= 10 )then
	          iDim1 = 0
	        elseIF (Pmin >= 1) then
	          iDim1 = 1
	        else
	          iDim1 = 2
	        endif
	        if (Pmax >= 10 )then
	          iDim2 = 0
	        elseIF (Pmax >= 1) then
	          iDim2 = 1
	        else
	          iDim2 = 2
	        endif
	        if (E2 >= 10) then
	          iDim3 = 0
	        elseif (E2 >= 1) then
	          iDim3 = 1
	        else
	          iDim3 = 2
	        endif
	        if (E4 >= 10) then
	          iDim4 = 0
	        elseif (E4 >= 1) then
	          iDim4 = 1
	        else
	          iDim4 = 2
	        endif
	        write (2, '(2(I2, ","), F6.0, "-", F6.0, ",", F6.'//i2s(iDim1)//', "-", F6.'//i2s(iDim2)//', ",", F6.'//i2s(iDim3)//', ",", F6.'//i2s(iDim4)//', ",", I3, ",", 80A2)')  &
	        
	            idc(ks,1), idc(ks,2), minval(Texp(N1:NCXP)),maxval(Texp(N1:NCXP)), Pmin, Pmax, E2, E4, nObs(ks), (NTEXT(ks,j), j = 1, 65)

	      elseIF(IFL(KS).EQ.2) then

!	        Bubble temperature (P and x)
	        NCXPy = NCXP
	        WRITE(1,356) IDC(KS,1),IDC(KS,1),IDC(KS,1)
356	        FORMAT(/,X,'  Pres/ATM    Xexp(',I0,')    Temp/K     TempCalc       ERR%       Yexp(',I0,')     Ycalc(',I0,')     ERR ',/)
	        write(10,'(/,"P /atm,(x-y)(",I0,")exp,T /K,(xExp-yCalc)(",I0,"),Tcalc /K")')IDC(KS,1),IDC(KS,1)
	        DO I=1,NOBS(KS)

	          WRITE(1,358) PEXP(NCXP+1),XEXP(NCXP+1),TEXP(NCXP+1),PCALC(NCXP+1), FCAUX(NCXP+1),YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&	                       ZLCALC(NCXP+1), ZVCALC(NCXP+1)
358	          FORMAT(F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,'  ZLIQ = ',F10.8,'      ZVAP = ',F10.8)
!
!	          x in graficos.csv
	          write(10,3580)PEXP(NCXP+1),XEXP(NCXP+1),TEXP(NCXP+1),XEXP(NCXP+1),PCALC(NCXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
	        do i = 1,NObs(ks)

!	          y in graficos.csv
	          write(10,3580)PEXP(NCXPy+1),yEXP(NCXPy+1),TEXP(NCXPy+1),ycalc(NCXPy+1),PCALC(NCXPy+1)
	          NCXPy = NCXPy + 1

	        enddo
	        NO = nObs(ks)
	        N1 = NCXP - NO + 1
	        E1 = AAD(TExp(N1:NCXP), PCalc(N1:NCXP), NO)
	        E2 = ARD(TExp(N1:NCXP), PCalc(N1:NCXP), NO)*100
	        E3 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
	        E4 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
	        write(1,360)E1,E2,idC(ks,1),E3,E4
	        write(*,'(/,X," ARD(T) =",f10.4,20X,"ARD(y) = ",f10.4)')E2, E4
	        write(10,360)E1,E2,idC(ks,1),E3,E4
360	        Format(/" AAD(ARD%): T = ,",E10.4,", K (, ",F6.2,"% ,)  ,",20X," y(",I0,") =, ",E10.4,", (, ",F6.2,"% ,)")        
3580	        format(F10.4,",",F9.6,",",F7.2,",",F9.6,",",F7.2)

	      elseIF(iFL(KS).EQ.3)then

!	        Ternary equilibrium
	        WRITE(1,451) IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
451	        FORMAT(/,X,"  T/K     P/ATM   xExp(",I0,',I) xClc(',I0,',I)    ERR    xExp(',I0,',I) xClc(',I0,',I)    ERR    xExp(',I0,',II)  xClc(',I0,',II)      ERR     xExp(',I0,',II)  xClc(',I0,',II)     ERR%',/)
!
	        WRITE(10,4511) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,2)
4511	        FORMAT(/,X,"xExp(",I0,"),xExp(",I0,"),xClc(",I0,"),xClc(",I0,"),T /K,P /atm")
!
!	        Second counter, for graficos.csv
	        NCXPy = NCXP
	        DO I=1,NOBS(KS)
!
	          er2 = abs(XEXP2(NCXP+1) - XCALC2(NCXP+1))
	          er4 = abs(YEXP2(NCXP+1) - YCALC2(NCXP+1))
	          WRITE(1,449) TEXP(NCXP+1),PEXP(NCXP+1),XEXP(NCXP+1),XCALC(NCXP+1),fcAux(ncXp+1),XEXP2(NCXP+1),XCALC2(NCXP+1),er2,YEXP(NCXP+1),YCALC(NCXP+1),FCAUX2(NCXP+1), &
&	                       YEXP2(NCXP+1),YCALC2(NCXP+1),ZLCALC(NCXP+1), ZVCALC(NCXP+1)
449	          FORMAT(F7.2,7F10.5,5(X,F10.5,X),2X,' ZI = ',G12.4,'       ZII =  ',G12.4)
!
!	          x in graficos.csv
	          write(10,4490)XEXP(NCXP+1),xExp2(NCXP+1),xCalc(ncXP+1),xCalc2(ncXP+1),TEXP(NCXP+1),PEXP(NCXP+1),yExp(ncXp+1),yExp2(ncXp+1), yCalc(ncXP+1),yCalc2(ncXP+1)
	          NCXP = NCXP + 1

	        END DO
	        write(10,*)
	        NO = nObs(ks)
	        N1 = NCXP - NO + 1
	        E1 = AAD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)
	        E2 = ARD(xExp(N1:NCXP), xCalc(N1:NCXP), NO)*100
	        E3 = AAD(xExp2(N1:NCXP), xCalc2(N1:NCXP), NO)
	        E4 = ARD(xExp2(N1:NCXP), xCalc2(N1:NCXP), NO)*100
	        E5 = AAD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)
	        E6 = ARD(yExp(N1:NCXP), yCalc(N1:NCXP), NO)*100
	        E7 = AAD(yExp2(N1:NCXP), yCalc2(N1:NCXP), NO)
	        E8 = ARD(yExp2(N1:NCXP), yCalc2(N1:NCXP), NO)*100   ! it used yExp/yCalc before (repeating E6)
	        write(1,453)idC(ks,1),E1,E2, idC(ks,2),E3,E4, idC(ks,1),E5,E6, idC(ks,2),E7,E8
	        write(*,'(/,X," ARD(x) =",f10.4,20X,"ARD(y) = ",f10.4)')E2,E4,E6,E8
	        write(10,453)idC(ks,1),E1,E2, idC(ks,2),E3,E4, idC(ks,1),E5,E6, idC(ks,2),E7,E8
453	        Format(/" AAD(ARD%): x(",I0,") = ,",E10.4,", ( ",F7.2,"%, ) ,  x(",I0,") = ,",E10.4," ,( ",F7.2,"%, )",20X,         &
&                    " y(",I0,") = ,",E10.4," ,(, ",F7.2,"% ,) ,  y(",I0,") = ,",E10.4," ( ",F7.2,"% )")
4490	        format(4(F9.6,","),F7.2,",",F10.4,/,4(F9.6,","),/)

	      END IF
!
	    elseIF (IHEN(KS).EQ.1)THEN

	      if(ifl(ks).EQ.1)then

!	        Excess volume

	        write(1,'(/,"  T(K)       P(atm)    Xexp(",I0,")  Vexp^E(cm3/mol)   Vcalc^E(cm3/mol)       Z",/)')idC(ks,1)
	        write(10,'(/,"Temp (K),P (atm),Xexp(",I0,"),Vexp^E (cm3/mol),Vcalc^E (cm3/mol),Z",/)')idC(ks,1)

	      elseif(ifl(ks).EQ.2)then

!	        Excess enthalpy
	        write(1,'(/,"  T(K)       P(atm)    Xexp(",I0,")  Hexp^E(J/mol)   Hcalc^E(J/mol)       Z",/)')idC(ks,1)
	        write(10,'(/,"T (K),P (atm),Xexp(",I0,"),Hexp^E (J/mol),Hcalc^E (J/mol),Z",/)')idC(ks,1)

	      endif
	      do i = 1,nObs(ks)

	        write(1,'(F7.2,2X,F10.5,2X,F10.8,X,2(F10.3,8X),F11.8)')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),prop_E_exp(NCXP+1),prop_E_calc(NCXP+1),zLiCalc(NCXP+1)
	        write(10,'(F7.2,",",F10.5,",",F10.8,",",2(F10.2,","),F11.8)')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),prop_E_exp(NCXP+1),prop_E_calc(NCXP+1),zLiCalc(NCXP+1)
	        NCXP = NCXP + 1

	      enddo
	      NO = nObs(ks)
	      N1 = NCXP - NO + 1
	      E1 = AAD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)	! average absolute deviation
	      E2 = ARD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)*100	! average relative deviation
	      E3 = ASD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)	! standard deviation
	      E4 = RSD(prop_E_exp(N1:NCXP), prop_E_calc(N1:NCXP), NO)*100 ! relative standard deviation
	      if(ifl(ks).EQ.1)then
	        write(1,5000)E1,E2,E3,E4
	        write(10,5000)E1,E2,E3,E4
	        write(*,'("(ARD(VE) = ",F10.4)')E2
	      elseif(ifl(ks).EQ.2)then
	        write(1,5001)E1,E2,E3,E4
	        write(10,5001)E1,E2,E3,E4
	        write(*,'("ARD(HE) = ",F10.4)')E2
	      endif
5000	      format(/" AAD(ARD%) = ",F8.4," cm3/mol (, ",F8.2,"% ,)",10X," StdDev(%) = ",F8.4," cm3/mol (, ",F8.2,"% ,)")
5001	      format(/" AAD(ARD%) = ",F8.4," J/mol ,(, ",F8.2,"% ,)",10X," StdDev(%) = ",F8.4," J/mol ,(, ",F8.2,"% ,)")
!
	    elseif(iHen(ks).EQ.2)then

!	      Solubility of solids in liquids
	      write(1,'(/,2X,"T(K)",6X,"P(atm)",5X,"xExp(",I0,")",5x,"xCal(",I0,")",8X,"Zliq",7X,"PhiL(",I0,")",/)')idC(ks,1),idC(ks,1),idC(ks,1)
	      write(10,'(/,"xExp(",I0,"),xCalc(",I0,"),T(K),P(atm)",/)')idC(ks,1),idC(ks,1)
	      do i = 1, nObs(ks)
		  
	        write(1,'(F7.2,2X,F10.5,2X,2(F10.8,2X),2(F11.8,2x))')Texp(NCXP+1),Pexp(NCXP+1),xExp(NCXP+1),xCalc(NCXP+1),ZliCalc(NCXP+1),yCalc(NCXP+1)
	        write(10,'(2(F10.8,","),F7.2,",",F10.5)')xExp(NCXP+1),xCalc(NCXP+1),Texp(NCXP+1),Pexp(NCXP+1)
	        NCXP = NCXP + 1
		  
	      enddo
	      NO = nObs(ks)
	      N1 = NCXP - NO + 1
	      E1 = AAD(xExp(N1:NCXP),xCalc(N1:NCXP),NO)	! average absolute deviation
	      E2 = ARD(xExp(N1:NCXP),xCalc(N1:NCXP),NO)*100	! average relative deviation
	      write(1,5002)idC(ks,1),E1,E2
	      write(*,'(/"ARD(xS) = ",F10.4)')E2
	      write(10,5002)idC(ks,1),E1,E2
5002	      format(/" AAD(ARD%) xS(",I0,") = ,",E10.4," ,(, ",F6.2,"% ,)")
	
	    elseif(iHen(ks) == 3)then
	    
! 	      Critical points
	      write(1,'(/, 2X, "Tc(K)", 5X, "Pc(atm)", 4X, "xcExp(", I0, ")", 8X, "Zc", 4X, "d[ln(f(", I0, "))]/dz(", I0, ")", 2X, "d2[ln(f(", I0, "))]/dz(", I0, ")2",  /)')idC(ks,1), idC(ks,1), idC(ks,1), idC(ks,1), idC(ks,1)
	      do i = 1, nObs(ks)
	      
	        write(1,'(F7.2, 2X, F10.5, 2X, 2(F10.8, 2X), 2(F11.8, 2x))')Texp(NCXP+1), Pexp(NCXP+1), xExp(NCXP+1),  ZliCalc(NCXP+1), fcAux(NCXP+1), fcAux2(NCXP+1)
	        NCXP = NCXP + 1
	        
	      enddo
	    
	    elseIF (IHEN(KS).EQ.5)THEN
!
!	      Activity coefficients
	      If(mez(ks)==2)then

	        write(1,1344)IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
1344	        FORMAT(X,'Temp/K       Xex(',I0,')   Gex(',I0,')   Gc(',I0,')     Gex(',I0,')   Gc(',I0,')'       ,'Zliq'/)
1444	        format('T/K,x',I0,'),ln(gamma(',I0,')exp),ln(gamma(',I0,'))calc,ln(gamma(',I0,')exp),ln(gamma(',I0,')calc),Z')
	        write(10,1444)IDC(KS,1),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2)
	        DO I=1,NOBS(KS)

	          WRITE(1,1345)TEXP(NCXP+1),XEXP(NCXP+1),rgamexp1(ncxp+1),rgamcalc1(NCXP+1),rgamexp2(ncxp+1),rgamcalc2(NCXP+1), zlicalc(NCXP+1)
		    WRITE(10,1445)TEXP(NCXP+1),XEXP(NCXP+1),dlog(rgamexp1(ncxp+1)),dlog(rgamcalc1(NCXP+1)),dlog(rgamexp2(ncxp+1)),dlog(rgamcalc2(NCXP+1)), zlicalc(NCXP+1)
1345	          format(X, 2f10.2, 4f10.3,f10.3)
1445	          format(F7.2, ",", 6(G10.4,","))
1347	          FORMAT(2F8.4,f7.1,2F8.4,f6.3,f8.4)
 	          NCXP = NCXP + 1

	        END DO
	        NO = nObs(ks)
	        N1 = NCXP - NO + 1
	        E1 = AAD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)
	        E2 = ARD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)*100
	        E3 = AAD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)
	        E4 = ARD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)*100
	        write(1,1379)idC(ks,1),E1,E2,idC(ks,2),E3,E4
	        write(10,1379)idC(ks,1),E1,E2,idC(ks,2),E3,E4
	        write(*,'(" ARD(gamma1) = ",F10.4,20X," ARD(gamma2) = ",F10.4)')E2,E4
1379	        format(/" AAD(ARD%) gamma(",I0,") = ",E10.4,", (, ",F6.2"% ,) , AAD(ARD%) gamma(",I0,") = ",E10.4,", (, ",F6.2"% ,)")

	      elseif(mez(ks)==3)then

	        write(1,1350) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,3),IDC(KS,3)
	        write(10,1350) IDC(KS,1),IDC(KS,2),IDC(KS,1),IDC(KS,1),IDC(KS,2),IDC(KS,2),IDC(KS,3),IDC(KS,3)
	        DO I=1,NOBS(KS)

	          WRITE(1,1348)TEXP(NCXP+1),XEXP(NCXP+1),xexp2(ncxp+1),rgamexp1(NCXP+1),rgamcalc1(NCXP+1),rgamexp2(NCXP+1),rgamcalc2(NCXP+1), &
&	                       rgamexp3(NCXP+1),rgamcalc3(NCXP+1)
1348	          FORMAT(X, f6.2 , 2(1X,e10.4), 6(1X,E12.3))
1349	          FORMAT(3F9.6)
1350	          FORMAT(X,'Temp/K Xex(',I0,') Xex(',I0,') Gex(',I0,') Gc(',I0,') Gex(',I0,') Gc(',I0,'),    &
&                      Gex(',I2,') Gc(',I2,')')
                NCXP = NCXP + 1

	        END DO
	        NO = nObs(ks)
	        N1 = NCXP - NO + 1
	        E1 = AAD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)
	        E2 = ARD(rGamExp1(N1:NCXP), rGamCalc1(N1:NCXP), NO)*100
	        E3 = AAD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)
	        E4 = ARD(rGamExp2(N1:NCXP), rGamCalc2(N1:NCXP), NO)*100
	        E5 = AAD(rGamExp3(N1:NCXP), rGamCalc3(N1:NCXP), NO)   ! it used component 2 before
	        E6 = ARD(rGamExp3(N1:NCXP), rGamCalc3(N1:NCXP), NO)*100
	        write(1,1351)idC(ks,1),E1,E2,idC(ks,2),E3,E4,idC(ks,3),E5,E6
	        write(10,1351)idC(ks,1),E1,E2,idC(ks,2),E3,E4,idC(ks,3),E5,E6
	        write(*,'(" ARD(gamma1) = ",F10.4,20X," ARD(gamma2) = ",F10.4,20X," ARD(gamma3) = ",F10.4)') E2,E4,E6
1351	        format(/" AAD(ARD%) gamma(",I0,"): ",E10.4,", (, ",F6.2"%, ) , gamma(",I0,"): ",E10.4," (, ",F6.2,                  &
&                     "% ) , gamma(",I2,"): ",E10.4,", (, ",F6.2"% ,)")

	      end if

	    END IF

	  END DO
!
	END IF
!
	CLOSE(UNIT=1)
	return
ENDsubroutine write_results
!--------------------------------------------------------------------------------------------------------------
!
!	Replaces PARGR of the GCA package: instead of reading the group
!	parameters from a file, it takes those of paraest and puts them in the
!	package's COMMON blocks (a dummy argument cannot be in a COMMON, so they
!	are copied).
!
subroutine fakePARGR(NC,NG,NGA,NST,IOUt,NUWW,QAREA,TS,GS,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,ONUASOCWW, &
&	               MASSOC,ENASS,RKASS)
	use utilities, only: i2s, rep

	IMPLICIT REAL*8(A-H,O-Z)
	parameter(NCM=30,NGM=20,NGAM=12,NSM=24)

	DIMENSION qArea(NG),Ts(NG),nuWW(NG,NC),Gs(NG),GDOT(NG),GDDOT(NG),rKij(NG,NG),rKTij(NG,NG), &
&	          alfaij(NG,NG)
!
!	Association energy and volume with 4 indices: only in this subroutine.
	dimension enAss(2,NGA,2,NGA),rkAss(2,NGA,2,NGA),onuAsocWW(NGA,NC),MAssoc(NGA)
	integer sigma
	real*8 kappa
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),XKIJ(NGM,NGM),      &
&	              AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)

	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM)
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
	common/NAME/AName(NCM)
	COMMON/COORD/ZZ

	if(NC.GT.NCM)then

	  write(*,*)'Error: maximum number of components exceeded in GC3'
	  stop

	elseif(NG.GT.NGM)then

	  write(*,*)'Error: maximum number of attractive groups exceeded in GC3'
	  stop

	elseif(NGA.GT.NGAM)then

	  write(*,*)'Error: maximum number of associating groups exceeded in GC3'
	  stop

	endif
	iOut = 20
!
!	Coordination number
	ZZ=10.D0
!
!	Gas constant
	R=82.05D0	!atm cm3/(mol K)
	NGROUP=35
	NownII = NG
	NownIJ = ((NG + 1)**2 - NG -1)/2 - NG
!
!	Dispersive parameters
12	WRITE(iOut,105) (J,J=1,NC)
	WRITE(iOut,'()')
	WRITE(iOut,120)
	do j = 1,NG
!
	  q(j) = qArea(j)
	  Tstr(j) = Ts(j)
	  ny(:NC,j) = nuWW(j,:NC)
	  gStr(j) = Gs(j)
	  g1(j) = GDOT(j)
	  g2(j) = GDDOT(j)
	  call spline_temperature(Tstr(j), g1(j), g2(j), Tspl(j), epx(j))
15	  WRITE(IOUT,106)j,TSTR(j),Q(j),GSTR(j),G1(j),G2(j),(NY(i,j),i=1,NC)
	  do k = 1,NG

	    xKij(j,k) = rKij(j,k)
	    aKij(j,k) = rKTij(j,k)
	    alfa(j,k) = alfaij(j,k)

	  enddo

	enddo

!-----------------------------------------------------------------------
!
!	Binary parameters
33	WRITE(iOut,115)
	do i = 1,NG

	  WRITE(IOUT,"(3X,'|',"//rep(NG,'F8.3')//"'|')") (XKIJ(I,J),J=1,NG)

	enddo
      WRITE(iOut,117)
	do i = 1,NG

	  WRITE(IOUT,"(3X,'|',"//rep(NG,'F8.3')//"'|')") (AKIJ(I,J),J=1,NG)

	enddo
	WRITE(IOUT,118)
	DO 23 I=1,NG

	  WRITE(IOUT,"(3X,'|',"//rep(NG,'F8.3')//"'|')") (ALFA(I,J),J=1,NG)

23	enddo
!-----------------------------------------------------------------------
!	Association parameters
!-------Conversion to the notation by site "m" (instead of by site "k" of
!	  group "i")
	NST = 0
	if(NGA.GT.0)then

	  m1 = 0
	  m2 = 0
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      m2 = m1
	      m1 = m1 + 1      	! current site "k" is site "m1"
	      do j = i,NGA 	! this is valid only because, so far, there has been only
	        if(i.EQ.j)then	! one type of site in each group
	          do l = k,Massoc(j)

	            m2 = m2 + 1 	! current site "l" is site "m2"
	            eps_R(m1,m2) = enAss(l,j,k,i)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = rkAss(l,j,k,i)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        else
	          do l = 1,Massoc(j)

	            m2 = m2 + 1 	! current site "l" is site "m2"
	            eps_R(m1,m2) = enAss(l,j,k,i)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = rkAss(l,j,k,i)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        endif
	      enddo

	    enddo
	  enddo
	  NST = sum(Massoc(:NGA))
!
!	  Sites of each component, in the notation 1..NST
	  do i = 1,NC

	    m1 = 0
	    do j = 1,NGA
	      do k = 1, MAssoc(j)

	        m1 = m1 + 1
	        sigma(m1,i) = onuAsocWW(j,i)

	      enddo
	    enddo

	  enddo
!
!-------Print the association parameters
	  WRITE(IOUT,217) NGA
217	  FORMAT(//,X,"Number of groups in the association term =",I3)
!
	  WRITE(IOUT,218)
218	  FORMAT(/,X,' Associating group composition ',/,30X,'Component No.')

	  write(iOut,'(" Group    No. of sites    ",*(I5,X))')((I), I=1,NC)
!
	  write(iOut,120)
	  DO I=1,NGA

	    WRITE(IOUT,219) I,MAssoc(I), (onuAsocWW(I,J), J=1,NC)
219	    FORMAT(2X,I3,16X,I3,2X,*(F5.0,X))

	  END DO
!
	  WRITE(iout,220)
 220	  FORMAT(/,X,'Association energy [epsilon/R] (K)',/)

!	  Association energies, as a matrix
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write(IOUT,"(3X,'|',"//rep(NST,'X,F7.1')//"'|')")((enAss(k,i,l,j),l = 1,MAssoc(j)),j=1,NGA)

	    enddo
	  enddo
!	  Association volumes
	  write(IOUT,'(/," Association volume [kappa] (cm3/mol)",/)')
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write(IOUT,"(3X,'|',"//rep(NST,'X,F7.5')//"'|')")((rKAss(k,i,l,j),l = 1,Massoc(j)),j=1,NGA)

	    enddo
	  enddo

	endif

115	FORMAT(/,1X,'Binary interaction parameters, [k*]',/)
117	FORMAT(/,1X,'Temperature dependence, [k`]',/)
118	FORMAT(/,1X,'Local composition parameters, [alpha]',/)
105	FORMAT(///,2X,'Attractive groups:',33X,'Component No.',/,3X,'GN',4X,'T*',6X,'q',5X,'g*',10X,'g`',8X,'g"',3X,*(I4))
106	FORMAT(1X,I4,F7.1,F8.4,F10.1,2F10.5,*(I4))
120	FORMAT(1X,90("-"))
	return
endsubroutine fakePARGR
!---------------------------------------------------------------------------------------------------------
!
!	Same as fakePARGR, for the component parameters of PARMOL (Tc, Pc, dc).
!
subroutine fakePARMOL(NC,NG,NST,IOUt,WWDC,WWTEMPCRIT,wwPRCRIT)
!
	implicit real*8(A-H,O-Z)
	parameter(NCM=30,NGM=20,NGAM=12,NSM=24)
	character*10 AName
	character*20 :: versGCA
!
	DIMENSION TSAT(NCM),PS(NCM),Z(2),OMEGA(NCM),NYOLD(NGM),EXPOL(NGM)
	DIMENSION wwDc(NC),wwTempCrit(NC),wwPRCRIT(NC)
!
	common/versSUB/versGCA
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM)
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
!
	common/NAME/AName(NCM)
	common/GCPROM/pmm(ncm),Pen(NCM),HHA(ncm),HHB(ncm),hhc(ncm),HHD(ncm),HHE(NCM),HHF(NCM)
	common/ENTALP/iEntalp
	COMMON/COORD/ZZ

	dimension eps_R_bak(NST,NST)
	real*8 kappa_bak(NST,NST),kappa
	integer sigma,sigma_bak(NST,NC)

	write(iOut,'(//,"Critical diameters",/)')
	iEntalp = 0
	do i = 1,NC

	  Tc(i) = wwTempCrit(i)
	  dc(i) = wwDC(i)
	  Pc(i) = wwPRCRIT(i)
	  Tsat(i) = 0.
	  Ps(i) = 0.
	  omega(i) = 0.
	  write(AName(i),'("Comp ",I2)')i
	  PMM(i) = 0.
	  Pen(i) = 0.
	  HHA(i) = 0.
	  HHB(i) = 0.
	  HHC(i) = 0.
	  HHD(i) = 0.
	  HHE(i) = 0.
	  HHF(i) = 0.

	enddo
!
!	iEntalp: number of components with Passut-Danner constants
	iEntalp = NC
	do i = 1,NC
	  if((HHA(i).EQ.0).AND.(HHB(i).EQ.0).AND.(HHC(i).EQ.0).AND.(HHD(i).EQ.0).AND.(HHE(i).EQ.0).AND.(HHF(i).EQ.0))then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	return
endsubroutine fakePARMOL
