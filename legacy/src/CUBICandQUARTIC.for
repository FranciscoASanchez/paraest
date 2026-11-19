			!*********************************************************************************!
	subroutine quarticequation(qa0,qa1,qa2,qa3,qa4,zq)
	implicit none
	integer jj,i
	DOUBLE PRECISION :: root1,zc(3),zq(4)
	DOUBLE precision :: qa0,qa1,qA2,qa3,qa4
	double complex :: cqa0,cqa1,cqA2,cqa3,cqa4
	DOUBLE precision ca0,ca1,ca2,ca3
	DOUBLE COMPLEX Q1,R1,E1,D1
	DOUBLE COMPLEX czq(4)
	
	zq=0D0
	root1=0D0
	czq=0d0
	zc=0d0

	If(qa4.eq.0d0.and.qa3.eq.0d0)then !quadratic 
		If((qa1**2-4d0*qa2*qa0)>=0d0)then
			zq(1)=(-qa1+dsqrt(qa1**2-4d0*qa2*qa0))/(2d0*qa2)
			zq(2)=(-qa1-dsqrt(qa1**2-4d0*qa2*qa0))/(2d0*qa2)
		end if
	elseif (qa4.eq.0d0.and.qa3.ne.0d0)then  !cubic 

	call cubicequation(qa3,qa2,qa1,qa0,zc)
		do i=1,3
			zq(i)=zc(i)
		end do
	else !quartic
	qa0=qa0/qa4
	qa1=qa1/qa4
	qa2=qa2/qa4
	qa3=qa3/qa4
	qa4=qa4/qa4

	ca3=1d0
	ca2=-1D0*qa2
	ca1=qa1*qa3-4.D0*qa0
	ca0=4.D0*qa2*qa0-qa1**2-qa3**2*qa0
	
	call cubicequation(ca3,ca2,ca1,ca0,zc)
	!select a real root 
	DO JJ=1,3
		If (zc(jj).ne.0d0)then  !real number
			root1=zc(jj)
		endif
	END DO

	cqa0=DCMPLX(qa0)
	cqa1=DCMPLX(qa1)
	cqa2=DCMPLX(qa2)
	cqa3=DCMPLX(qa3)
	cqa4=DCMPLX(qa4)

	R1=cdsqrt(0.25D0*cqa3**2-cqa2+root1)

	D1 = cdsqrt(0.5D0*cqa3**2-cqa2-root1+(4.D0*cqa3*cqa2-8.D0*cqa1
     @-cqa3**3)/cdsqrt(0.25D0*cqa3**2-cqa2+root1)/4.D0)
	E1=cdsqrt(0.5D0*cqa3**2-cqa2-root1-(4.D0*cqa3*cqa2-8.D0*cqa1
     @-cqa3**3)/cdsqrt(0.25D0*cqa3**2-cqa2+root1)/4.D0)


	!roots
	czq(1)=-0.25D0*cqa3+0.5D0*cdsqrt(0.25D0*cqa3**2-cqa2+root1)+0.5D0*
     @cdsqrt(0.5D0*cqa3**2-cqa2-root1+0.25D0*(4.D0*cqa3*cqa2-8.D0*cqa1
     @-cqa3**3)/cdsqrt(0.25D0*cqa3**2-cqa2+root1))

	czq(2)=-0.25D0*cqa3+0.5D0*cdsqrt(0.25D0*cqa3**2-cqa2+root1)-0.5D0*
     @cdsqrt(0.5D0*cqa3**2-cqa2-root1+0.25D0*(4.D0*cqa3*cqa2-8.D0*cqa1
     @-cqa3**3)/cdsqrt(0.25D0*cqa3**2-cqa2+root1))

	czq(3)=-0.25D0*cqa3-0.5D0*cdsqrt(0.25D0*cqa3**2-cqa2+root1)+0.5D0*
     @cdsqrt(0.5D0*cqa3**2-cqa2-root1-0.25D0*(4.D0*cqa3*cqa2-8.D0*cqa1
     @-cqa3**3)/cdsqrt(0.25D0*cqa3**2-cqa2+root1))

	czq(4)=-0.25D0*cqa3-0.5D0*cdsqrt(0.25D0*cqa3**2-cqa2+root1)-0.5D0*
     @cdsqrt(0.5D0*cqa3**2-cqa2-root1-0.25D0*(4.D0*cqa3*cqa2-8.D0*cqa1
     @-cqa3**3)/cdsqrt(0.25D0*cqa3**2-cqa2+root1))

	do jj=1, 4
		If (dimag(czq(jj)).eq.0d0)then
			zq(jj)=dreal(czq(jj))
		end if
	end do
	
	end if

	end subroutine
		!*********************************************************************************!
	subroutine cubicequation(a3,a2,a1,a0,z1)
	implicit none
	integer jj
	double precision a0,a1,a2,a3
	double precision Q1,R1,theta,pi,D1
	double precision t1,t4,t5,t7,t14,t19,t22,t31
	double precision z1(3),root1,e
	PARAMETER (pi=3.141592653589793238D0)
	
	z1=0D0
	root1=0D0
	If (a3.ne.0d0)then
		a0=a0/a3
		a1=a1/a3
		a2=a2/a3
		a3=a3/a3
		Q1 = (-a2 * a2 + 3d0 * a1) / 9d0
		R1=0D0
		R1 = (-2d0 * a2 * a2 * a2)
		R1=R1+ 9d0 * a2 * a1 
		R1=R1- 27d0 * a0
		R1=R1 / 54d0
		D1 = -Q1**3d0 - R1 * R1

		if (D1 >= 0) then  ! Three real roots
			 theta = acos(R1 / dsqrt(-Q1**3))
			
		z1(1)=2d0*dsqrt(-Q1)*cos(theta/3d0)-a2/3d0
		z1(2)=2d0*dsqrt(-Q1)*cos((theta+2d0*pi)/3d0)-a2/3d0
		z1(3)=2d0*dsqrt(-Q1)*cos((theta+4d0*pi)/3d0)-a2/3d0
		else   ! One real root
			t1 = a1*a2      
			t4 = a2*a2
			t5 = t4*a2
			t7 = a1*a1
			t14 = a0*a0
	
			t19 = dsqrt(12d0*t7*a1-3d0*t7*t4-54d0*t1*a0+81d0*t14+
     @			12d0*a0*t5)
			t22 = (36d0*t1-108d0*a0-8d0*t5+12d0*t19)**(1.D0/3.D0)
			z1(1) = t22/6d0-6d0*(a1/3d0-t4/9d0)/t22-a2/3d0

		end if

	ELSEIF(a3==0d0.AND.a2==0d0)then

		z1(1)=-a0/a1

	else	!quadratic equation
		z1(1)=(-a1+dsqrt(a1**2-4d0*a2*a0))/(2d0*a2)
		z1(2)=(-a1-dsqrt(a1**2-4d0*a2*a0))/(2d0*a2)


	end if

		WRITE(4,*) "Z1=",z1

	
	end subroutine