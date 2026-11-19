
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
	    enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
	enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
	  enddo
        b(i)=sum/a(i,i)
	enddo
      return
      END

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d = 1.
      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
	  enddo
C        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
	  if(aamax.EQ.0.)then
	    write(*,*)
	    write(*,*)"Matriz singular en LUDcmp"
	    write(*,*)
	    stop
	  endif
        vv(i)=1./aamax
	enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
	    enddo
          a(i,j)=sum
	  enddo
        aamax=0.
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum = sum - a(i,k)*a(k,j)
	    enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
	  enddo
        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
	    enddo
          d1=-d1
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)

          do i=j+1,n
            a(i,j)=a(i,j)*dum
	    enddo
        endif
	enddo
      return
      END
