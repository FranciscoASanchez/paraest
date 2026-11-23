!--------------------------------------------------------------------------------------------------------------------------
!
!  Module config: settings that can be changed in the optional file
!  paraest.cfg, a namelist such as
!
!     &options
!       penaliz   = 1.0
!       optimizer = 'pso+lm'
!     /
!
!  Variables left out keep their default values; without the file paraest
!  runs as always (Levenberg-Marquardt, PENALIZ = 0).
!
!--------------------------------------------------------------------------------------------------------------------------
module config
  implicit none
  private
  public :: read_options
  public :: PENALIZ, PENALIZNAN
  public :: optimizer, pso_particles, pso_iterations, pso_seed
  public :: range_gstr, range_gp, range_gdd, range_q, range_eps, range_tc

  ! --- Objective function ---------------------------------------------------
  ! With PENALIZ = 0 a point whose calculation fails adds nothing to the
  ! objective function, so a fit can lower it by driving points into
  ! failure; a positive value avoids that (but Levenberg-Marquardt copes
  ! badly with the jump; PSO does not mind).
  double precision :: PENALIZ = 0D0      ! residual of a point whose calculation failed
  double precision :: PENALIZNAN = 1D0   ! replaces NaN/Inf residuals

  ! --- Optimizer --------------------------------------------------------------
  !   'lm'      Levenberg-Marquardt (ZXSSQ) from the initial values
  !   'pso'     particle swarm (PSO) only: no statistics of the parameters
  !   'pso+lm'  PSO, then Levenberg-Marquardt from the best point found
  character(8) :: optimizer = 'lm'
  integer :: pso_particles = 35
  integer :: pso_iterations = 50
  integer :: pso_seed = 12345            ! same seed, same run

  ! --- PSO search box ---------------------------------------------------------
  ! Parameters with physical bounds (k*, k', alpha, association volume, dc)
  ! are searched between those bounds. g* , Q, eps and Tc between these
  ! multiples of the initial value given by the user. g' and g" between
  ! absolute values: their physical bounds (+-2) give g(T) < 0 (impossible,
  ! the calculation fails) at many temperatures; the box is widened to
  ! include the initial value if needed.
  double precision :: range_gstr(2) = [0.5D0, 1.5D0]   ! g*
  double precision :: range_q(2)    = [0.5D0, 2.0D0]   ! Q (area)
  double precision :: range_eps(2)  = [0.5D0, 2.0D0]   ! association energy
  double precision :: range_tc(2)   = [0.9D0, 1.1D0]   ! Tc
  double precision :: range_gp(2)   = [-1D0, 1D0]      ! g' (absolute)
  double precision :: range_gdd(2)  = [-1D0, 1D0]      ! g" (absolute)

  namelist /options/ PENALIZ, PENALIZNAN, optimizer, pso_particles, pso_iterations, pso_seed, &
                     range_gstr, range_q, range_eps, range_tc, range_gp, range_gdd

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Reads paraest.cfg, if it exists, and writes the values used to
!  paraest.out (which must already exist).
!
subroutine read_options()
  logical :: exists
  integer :: ios

  inquire(file='paraest.cfg', exist=exists)
  if (.not. exists) return
  open(unit=1, status='OLD', file='paraest.cfg', form='FORMATTED')
  read(1, nml=options, iostat=ios)
  close(unit=1)
  if (ios /= 0) then
    write(*,*) 'Error reading the namelist &options of paraest.cfg'
    stop 1
  endif
  if (optimizer /= 'lm' .and. optimizer /= 'pso' .and. optimizer /= 'pso+lm') then
    write(*,*) 'paraest.cfg: optimizer must be ''lm'', ''pso'' or ''pso+lm'', not ''', trim(optimizer), ''''
    stop 1
  endif
  if (pso_particles < 2 .or. pso_iterations < 1 .or. pso_seed <= 0) then
    write(*,*) 'paraest.cfg: pso_particles >= 2, pso_iterations >= 1 and pso_seed > 0 are required'
    stop 1
  endif

  open(unit=1, status='OLD', file='paraest.out', form='FORMATTED', position='APPEND')
  write(1,'(/,X,A,/,X,A,ES12.4,/,X,A,ES12.4)') 'Options read from paraest.cfg:', &
    'Residual of a failed point (PENALIZ):    ', PENALIZ, &
    'Residual replacing NaN/Inf (PENALIZNAN): ', PENALIZNAN
  if (optimizer /= 'lm') then
    write(1,'(X,2A)') 'Optimizer: ', trim(optimizer)
    write(1,'(X,A,I0,A,I0,A,I0)') 'PSO: particles = ', pso_particles, ', iterations = ', pso_iterations, &
      ', seed = ', pso_seed
    write(1,'(X,A,6(/,X,A,2F10.4))') 'PSO search box of the parameters without physical bounds:', &
      '  g*  (times the initial value): ', range_gstr, &
      '  Q   (times the initial value): ', range_q, &
      '  eps (times the initial value): ', range_eps, &
      '  Tc  (times the initial value): ', range_tc, &
      "  g'  (absolute):                ", range_gp, &
      '  g"  (absolute):                ', range_gdd
  endif
  close(unit=1)
endsubroutine read_options

endmodule config
