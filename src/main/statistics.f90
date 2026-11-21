!--------------------------------------------------------------------------------------------------------------------------
!
!  Module statistics: deviations between experimental and calculated data.
!
!  All of them take the experimental and calculated vectors and the number of
!  points N. Points that did not converge (calculated value zero) are left
!  out of the average: if none is left, they return 0.
!
!--------------------------------------------------------------------------------------------------------------------------
module statistics
  implicit none
  private
  public :: AAD, ARD, ASD, RSD

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Average absolute deviation. Keeps the units of the experimental data.
!
double precision function AAD(Xexp, Xcalc, N)
  integer, intent(in) :: N
  double precision, intent(in) :: Xexp(N), Xcalc(N)
  integer :: i, nReal

  AAD = 0D0
  nReal = 0
  do i = 1, N
    if (Xcalc(i) > 1D-20) then
      AAD = AAD + dabs(Xexp(i) - Xcalc(i))
      nReal = nReal + 1
    endif
  enddo
  if (nReal > 0) AAD = AAD/nReal
endfunction AAD

!--------------------------------------------------------------------------------------------------------------------------
!
!  Average relative deviation (dimensionless). Zero experimental values are
!  also left out.
!
double precision function ARD(Xexp, Xcalc, N)
  integer, intent(in) :: N
  double precision, intent(in) :: Xexp(N), Xcalc(N)
  integer :: i, nReal

  ARD = 0D0
  nReal = 0
  do i = 1, N
    if (dabs(Xcalc(i)) > 1D-16 .and. dabs(Xexp(i)) > 1D-16) then
      ARD = ARD + dabs(1D0 - Xcalc(i)/Xexp(i))
      nReal = nReal + 1
    endif
  enddo
  if (nReal > 0) ARD = ARD/nReal
endfunction ARD

!--------------------------------------------------------------------------------------------------------------------------
!
!  Absolute standard deviation. Keeps the units of the experimental data.
!
double precision function ASD(Xexp, Xcalc, N)
  integer, intent(in) :: N
  double precision, intent(in) :: Xexp(N), Xcalc(N)
  integer :: i, nReal

  ASD = 0D0
  nReal = 0
  do i = 1, N
    if (Xcalc(i) > 1D-20) then
      ASD = ASD + (Xexp(i) - Xcalc(i))**2
      nReal = nReal + 1
    endif
  enddo
  if (nReal > 0) ASD = dsqrt(ASD/nReal)
endfunction ASD

!--------------------------------------------------------------------------------------------------------------------------
!
!  Relative standard deviation (dimensionless). Here the points with a zero
!  experimental value are left out.
!
double precision function RSD(Xexp, Xcalc, N)
  integer, intent(in) :: N
  double precision, intent(in) :: Xexp(N), Xcalc(N)
  integer :: i, nReal

  RSD = 0D0
  nReal = 0
  do i = 1, N
    if (Xexp(i) > 1D-20) then
      RSD = RSD + (1D0 - Xcalc(i)/Xexp(i))**2
      nReal = nReal + 1
    endif
  enddo
  if (nReal > 0) RSD = dsqrt(RSD/nReal)
endfunction RSD

endmodule statistics
