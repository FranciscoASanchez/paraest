!> Utilidades generales de paraest.
module utilidades
  implicit none

contains

  !> Convierte un entero en texto, sin blancos: 12 -> '12'.
  !>
  !> Sirve para armar formatos cuyo número de repeticiones se conoce recién
  !> al ejecutar el programa (lo que en Intel se escribía como <n>):
  !>
  !>     write(*,'(/,'//i2s(ng)//'(F7.4,2X),"  !g*")') (g(i), i = 1, ng)
  !>
  !> El resultado es un texto "de largo diferido" (character(:), allocatable):
  !> su largo se ajusta solo a la cantidad de dígitos.
  function i2s(n) result(texto)
    integer, intent(in) :: n
    character(:), allocatable :: texto
    character(20) :: aux

    write(aux, '(I0)') n
    texto = trim(aux)
  end function i2s

end module utilidades
