!--------------------------------------------------------------------------------------------------------------------------
!
!  Spline temperature of a group: above Tspl the temperature dependence of
!  g_ii changes form. Same calculation as in PARGR (GCA package).
!
!     Tstr     reference temperature T* of the group (K)
!     g1, g2   g' and g" of the group
!     Tspl     spline temperature (K); 1000 K if g' >= -1e-4 or if Newton
!              does not converge
!     epx      coefficient of the high-temperature form (0 in those cases)
!
!--------------------------------------------------------------------------------------------------------------------------
subroutine spline_temperature(Tstr, g1, g2, Tspl, epx)
  implicit none
  double precision, intent(in) :: Tstr, g1, g2
  double precision, intent(out) :: Tspl, epx
  double precision :: TR, AX, DAX, absAX
  integer :: it

  Tspl = 1D3
  epx = 0D0
  if (g1 >= -1D-4) return               ! g' > 0: (hardly) any Tspl > 0

  ! Newton on Tr = Tspl/T*:  1 + g'(Tr - 1) + g" ln(Tr) - 1/4 = 0
  TR = 1D0 - 1D0/g1
  it = 0
  absAX = 1D0
  do while (absAX > 1D-15)
    AX = 1D0 + g1*(TR - 1D0) + g2*dlog(TR) - 0.25D0
    DAX = g1 + g2/TR
    TR = TR - AX/DAX
    it = it + 1
    absAX = dabs(AX)
    if (TR < 0D0 .or. it > 10) then     ! Tr < 0 or no convergence
      Tspl = 1D3
      epx = 0D0
      return
    endif
    Tspl = TR*Tstr
    epx = -4D0*(g1*TR + g2)
  enddo
endsubroutine spline_temperature
