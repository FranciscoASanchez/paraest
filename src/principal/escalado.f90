!--------------------------------------------------------------------------------------------------------------------------
!
!  Módulo escalado: transformación entre los parámetros del modelo y las
!  variables que ve el optimizador.
!
!    XXSCAL:    parámetros -> variables escaladas (antes de optimizar)
!    XXDESSCAL: variables escaladas -> parámetros (en cada evaluación)
!
!  El tipo de cada parámetro está en IDXX (lo arma FORSCAL):
!     1  g*  (g_ii)            2  k*  (k_ij)          3  g' (o k' si inIJ)
!     4  Q (área)              5  alfa                6  energía de asociación
!     7  volumen de asociación 8  Tc                  9  dc
!    10  g"
!
!  Los parámetros acotados se escalan con fScld: scl*log((x - vMin)/(vMax - x)).
!
!--------------------------------------------------------------------------------------------------------------------------
module escalado
  use dimensiones, only: NPARM
  implicit none
  private
  public :: XXSCAL, XXDESSCAL, fScld, fDscld
  public :: inIJ, dc_medio, error_Dc

  ! Factores de escalado (scl) de cada tipo de parámetro.
  ! Ojo: como el paso de las derivadas numéricas de ZXSSQ es casi fijo en
  ! la variable escalada, estos factores controlan en la práctica ese paso
  ! en la variable original (paso en log(...) ~ 4e-10/scl).
  double precision, parameter :: fGgg   = 1D-10
  double precision, parameter :: fRkij  = 1D-7
  double precision, parameter :: fTemp  = 1D-7
  double precision, parameter :: fQarea = 1D0
  double precision, parameter :: fAlfa  = 1D-7
  double precision, parameter :: fEnAss = 1D-6
  double precision, parameter :: fRkass = 1D-8
  double precision, parameter :: fTc    = 1D-8
  double precision, parameter :: fDc    = 2D-7

  ! Cotas de los parámetros acotados.
  double precision, parameter :: vMinKstr = 0.2D0, vMaxKstr = 30D0   ! k*
  double precision, parameter :: vMaxKp   = 9D0                      ! k' entre -vMaxKp y vMaxKp
  double precision, parameter :: vMinGp   = -2D0, vMaxGp = 2D0       ! g' entre -vMaxGp y -vMinGp
  double precision, parameter :: vMaxAlfa = 29D0                     ! alfa entre -vMaxAlfa y vMaxAlfa
  double precision, parameter :: vMaxRkass = 20D0                    ! volumen de asoc. entre 0 y vMaxRkass

  ! Distancia relativa mínima a las cotas: el logaritmo nunca recibe 0.
  double precision, parameter :: margen = 1D-10

  ! inIJ(i): el parámetro i de tipo 3 es un k' (interacción i-j) y no un g'.
  logical :: inIJ(NPARM) = .false.

  ! dc: cada uno se mueve en dc_medio(i) +- error_Dc(i) (los fija READ).
  double precision :: dc_medio(NPARM) = 0D0, error_Dc(NPARM) = 0D0

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Transforma los valores de XX según los factores de escalado.
!
!  fracCota: fracción del intervalo de cada parámetro acotado que se deja
!  libre junto a las cotas (ver fScld). 0 en predicción, 1e-3 en ajuste.
!
subroutine XXSCAL(NPAR, XX, IDXX, fracCota)
  integer, intent(in) :: NPAR, IDXX(NPAR)
  double precision, intent(inout) :: XX(NPAR)
  double precision, intent(in) :: fracCota
  integer :: i

  do i = 1, NPAR
    select case (IDXX(i))
    case (1)
      XX(i) = fGgg*XX(i)
    case (2)
      XX(i) = fScld(XX(i), fRkij, vMinKstr, vMaxKstr, fracCota)
    case (3)
      if (inIJ(i)) then
        XX(i) = fScld(XX(i), fTemp, -vMaxKp, vMaxKp, fracCota)
      else
        XX(i) = fScld(XX(i), fTemp, -vMaxGp, -vMinGp, fracCota)
      endif
    case (4)
      XX(i) = fQarea*XX(i)
    case (5)
      XX(i) = fScld(XX(i), fAlfa, -vMaxAlfa, vMaxAlfa, fracCota)
    case (6)
      ! Logarítmico: solo tiene cota inferior (0). max evita el log de 0 o
      ! de un negativo. En un ajuste, al menos 1 K (despreciable frente a
      ! valores típicos de 1000-3000 K): desde 0 no se podría mover.
      if (fracCota > 0D0) then
        XX(i) = fEnAss*dlog(max(XX(i), 1D0))
      else
        XX(i) = fEnAss*dlog(max(XX(i), tiny(1D0)))
      endif
    case (7)
      XX(i) = fScld(XX(i), fRkass, 0D0, vMaxRkass, fracCota)
    case (8)
      XX(i) = fTc*XX(i)
    case (9)
      XX(i) = fScld(XX(i), fDc, dc_medio(i) - error_Dc(i), dc_medio(i) + error_Dc(i), fracCota)
    case (10)
      XX(i) = fTemp*XX(i)       ! g": de momento sin cotas
    endselect
  enddo
endsubroutine XXSCAL

!--------------------------------------------------------------------------------------------------------------------------
!
!  Antitransforma XX y deja el resultado en XXAUX. XX no se toca: FUNC no
!  debe alterar los valores que le pasa ZXSSQ.
!
subroutine XXDESSCAL(NPAR, XX, IDXX, XXAUX)
  integer, intent(in) :: NPAR, IDXX(NPAR)
  double precision, intent(in) :: XX(NPAR)
  double precision, intent(out) :: XXAUX(NPAR)
  integer :: i

  do i = 1, NPAR
    select case (IDXX(i))
    case (1)
      XXAUX(i) = XX(i)/fGgg
    case (2)
      XXAUX(i) = fDscld(XX(i), fRkij, vMinKstr, vMaxKstr)
    case (3)
      if (inIJ(i)) then
        XXAUX(i) = fDscld(XX(i), fTemp, -vMaxKp, vMaxKp)
      else
        XXAUX(i) = fDscld(XX(i), fTemp, -vMaxGp, -vMinGp)
      endif
    case (4)
      XXAUX(i) = XX(i)/fQarea
    case (5)
      XXAUX(i) = fDscld(XX(i), fAlfa, -vMaxAlfa, vMaxAlfa)
      ! Un alfa justo en la cota se aparta apenas (1e-8 relativo).
      if (vMaxAlfa - dabs(XXAUX(i)) < 1D-8) XXAUX(i) = sign((1D0 - 1D-8)*vMaxAlfa, XXAUX(i))
    case (6)
      XXAUX(i) = dexp(min(XX(i)/fEnAss, log(huge(1D0))))   ! min: exp no desborda
    case (7)
      XXAUX(i) = fDscld(XX(i), fRkass, 0D0, vMaxRkass)
      ! Igual que alfa: no se deja justo en la cota superior.
      if (dabs(vMaxRkass - XXAUX(i)) < 1D-8) XXAUX(i) = (1D0 - 1D-8)*vMaxRkass
    case (8)
      XXAUX(i) = XX(i)/fTc
    case (9)
      XXAUX(i) = fDscld(XX(i), fDc, dc_medio(i) - error_Dc(i), dc_medio(i) + error_Dc(i))
    case (10)
      XXAUX(i) = XX(i)/fTemp
    case default
      XXAUX(i) = XX(i)
    endselect
  enddo
endsubroutine XXDESSCAL

!--------------------------------------------------------------------------------------------------------------------------
!
!  Escalado de una variable acotada.
!
!  vMin y vMax son las cotas mínima y máxima que puede tomar x.
!  scl estira o comprime el valor devuelto (ver la nota de los factores).
!
!  Para que nunca dé una excepción (log de 0 o de un negativo), x se lleva
!  primero al interior del intervalo, a una distancia relativa "margen" de
!  las cotas. Así no importa si el valor inicial está justo en una cota
!  (por ejemplo kappa = 0) o si el redondeo lo deja en la cota o afuera.
!
!  frac es la fracción del intervalo que se deja libre junto a cada cota,
!  si es mayor que ese margen. En un ajuste se usa frac = 1e-3: en la cota
!  la derivada es casi nula y el optimizador no movería el parámetro.
!  En predicción, frac = 0 (se respeta el valor dado).
!
double precision function fScld(x, scl, vMin, vMax, frac)
  double precision, intent(in) :: x, scl, vMin, vMax, frac
  double precision :: d, xin

  d = max(frac, margen)*(vMax - vMin)
  xin = min(max(x, vMin + d), vMax - d)
  fScld = scl*dlog((xin - vMin)/(vMax - xin))
endfunction fScld

!--------------------------------------------------------------------------------------------------------------------------
!
!  Desescalado: inversa de fScld.
!
!  Se calcula de dos maneras según el signo de x/scl, para que la
!  exponencial siempre tenga argumento <= 0: nunca desborda (con scl muy
!  chico, exp(x/scl) daba infinito, y luego infinito/infinito = NaN).
!  El resultado queda siempre entre vMin y vMax.
!
!  fScld aparta x de las cotas un margen relativo; si el resultado cae
!  dentro de ese margen (el doble, por redondeo), se devuelve la cota
!  misma. Así un valor dado justo en la cota (por ejemplo kappa = 0) vuelve
!  igual, y el usuario no ve el margen.
!
double precision function fDscld(x, scl, vMin, vMax)
  double precision, intent(in) :: x, scl, vMin, vMax
  double precision :: t, ex, d

  t = x/scl
  if (t >= 0D0) then
    ex = dexp(-t)
    fDscld = (vMin*ex + vMax)/(1D0 + ex)
  else
    ex = dexp(t)
    fDscld = (vMin + vMax*ex)/(1D0 + ex)
  endif

  d = 2D0*margen*(vMax - vMin)
  if (fDscld - vMin <= d) fDscld = vMin
  if (vMax - fDscld <= d) fDscld = vMax
endfunction fDscld

endmodule escalado
