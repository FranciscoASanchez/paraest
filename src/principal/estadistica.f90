!--------------------------------------------------------------------------------------------------------------------------
!
!  Módulo estadistica: desviaciones entre datos experimentales y calculados.
!
!  Todas reciben los vectores experimental y calculado y el número de datos N.
!  Los puntos que no convergieron (valor calculado nulo) se excluyen del
!  promedio: si no queda ninguno, devuelven 0.
!
!--------------------------------------------------------------------------------------------------------------------------
module estadistica
  implicit none
  private
  public :: AAD, ARD, ASD, RSD

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Desviación absoluta promedio. Mantiene las unidades del dato experimental.
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
!  Desviación relativa promedio (sin unidades). Se excluyen también los
!  datos experimentales nulos.
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
!  Desviación estándar absoluta. Mantiene las unidades del dato experimental.
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
!  Desviación estándar relativa (sin unidades). Aquí se excluyen los puntos
!  con dato experimental nulo.
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

endmodule estadistica
