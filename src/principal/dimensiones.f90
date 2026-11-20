!> Tamaños máximos de los arreglos de paraest.
!>
!> Antes cada subrutina repetía su propio "integer, parameter :: NCM = 50, ...";
!> ahora se definen una sola vez acá y cada subrutina hace "use dimensiones".
module dimensiones
  implicit none

  integer, parameter :: NCM   = 50     !< componentes
  integer, parameter :: NGM   = 20     !< grupos dispersivos
  integer, parameter :: NDM   = 1520   !< datos experimentales (en total)
  integer, parameter :: NSetM = 100    !< conjuntos de datos

  !> Parámetros ajustables (en total). Hay 7 tipos (g, q, alfa, energía y
  !> volumen de asociación, Tc, dc) con hasta 8 de cada uno; los g pueden
  !> sumar además hasta 2 términos de dependencia con T: 8*3 + 6*8 = 72.
  integer, parameter :: NPARM = 72

  !> Largo del vector de trabajo de ZXSSQ: 5*N + 2*M + N*(N+1)/2,
  !> con N = parámetros y M = datos.
  integer, parameter :: LWORK = 5*NPARM + 2*NDM + NPARM*(NPARM + 1)/2

end module dimensiones
