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

  !> Arma un grupo repetido de un formato, seguido de una coma:
  !>
  !>     rep(3, 'F7.4,2X')  ->  '3(F7.4,2X),'
  !>     rep(0, 'F7.4,2X')  ->  ''            (el grupo desaparece)
  !>
  !> Uso típico:
  !>
  !>     write(*,'(/,'//rep(ng,'F7.4,2X')//'"  !g*")') (g(i), i = 1, ng)
  !>
  !> Con n = 0 no se puede escribir '0(...)' en un formato; Intel lo
  !> aceptaba en su forma <n>, así que se devuelve el texto vacío para que
  !> el resto del formato se escriba igual.
  function rep(n, descriptor) result(texto)
    integer, intent(in) :: n
    character(*), intent(in) :: descriptor
    character(:), allocatable :: texto

    if (n > 0) then
      texto = i2s(n)//'('//descriptor//'),'
    else
      texto = ''
    end if
  end function rep

end module utilidades
