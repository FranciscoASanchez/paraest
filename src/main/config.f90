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
  public :: optimizer, pso_particles, pso_iterations, pso_seed, pso_starts, pso_ssq_factor, pso_cluster_dist
  public :: range_gstr, range_kstr, range_alpha, range_gp, range_gdd, range_q, range_eps, range_tc
  public :: range_kappa
  public :: bubble_solver

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
  ! 'pso+lm': Levenberg-Marquardt from up to pso_starts points of the PSO in
  ! distinct regions of the box. Candidates: the evaluations with at most one
  ! point lost more than the best, and SSQ <= pso_ssq_factor times the best.
  ! Taken from best to worst, a candidate starts LM if its distance to every
  ! start already chosen is larger than pso_cluster_dist (RMS of the
  ! differences divided by the width of the box of each parameter). The LM
  ! results closer than that are the same optimum.
  integer :: pso_starts = 5
  double precision :: pso_ssq_factor = 2D0
  double precision :: pso_cluster_dist = 0.1D0

  ! --- PSO search box ---------------------------------------------------------
  ! g*, Q, eps and Tc: between these multiples of the initial value given
  ! by the user. k*, alpha, g' and g": between absolute values, within their
  ! physical bounds and widened to include the initial value if needed (the
  ! bounds used by LM are wider: k* 0.2-30, alpha +-29, g' +-2; k* should
  ! statistically be near 1, and g' or g" near +-2 give g(T) < 0 at many
  ! temperatures). k', association volume and dc: their physical bounds.
  double precision :: range_gstr(2) = [0.5D0, 1.5D0]   ! g*
  double precision :: range_q(2)    = [0.5D0, 2.0D0]   ! Q (area)
  double precision :: range_eps(2)  = [0.5D0, 2.0D0]   ! association energy
  double precision :: range_tc(2)   = [0.9D0, 1.1D0]   ! Tc
  double precision :: range_kstr(2) = [0.8D0, 1.2D0]   ! k* (absolute)
  double precision :: range_alpha(2) = [-10D0, 10D0]   ! alpha (absolute)
  double precision :: range_gp(2)   = [-1D0, 1D0]      ! g' (absolute)
  double precision :: range_gdd(2)  = [-1D0, 1D0]      ! g" (absolute)

  ! --- Bounds of the association volume (kappa, cm3/mol), for LM and PSO -----
  ! (not only the PSO box: LM cannot leave them either)
  double precision :: range_kappa(2) = [0D0, 20D0]

  ! --- Bubble points of binary VLE data ----------------------------------------
  !   'tpf'    Newton-Raphson in ln K and ln P (or ln T), from TPfractionFlash
  !   'gcadb'  the old GCADB routine
  character(8) :: bubble_solver = 'tpf'

  namelist /options/ PENALIZ, PENALIZNAN, optimizer, pso_particles, pso_iterations, pso_seed, &
                     pso_starts, pso_ssq_factor, pso_cluster_dist, &
                     range_gstr, range_q, range_eps, range_tc, range_kstr, range_alpha, range_gp, range_gdd, &
                     range_kappa, bubble_solver

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
  if (bubble_solver /= 'tpf' .and. bubble_solver /= 'gcadb') then
    write(*,*) 'paraest.cfg: bubble_solver must be ''tpf'' or ''gcadb'', not ''', trim(bubble_solver), ''''
    stop 1
  endif
  if (pso_particles < 2 .or. pso_iterations < 1 .or. pso_seed <= 0) then
    write(*,*) 'paraest.cfg: pso_particles >= 2, pso_iterations >= 1 and pso_seed > 0 are required'
    stop 1
  endif
  if (range_kappa(1) < 0D0 .or. range_kappa(2) <= range_kappa(1)) then
    write(*,*) 'paraest.cfg: range_kappa needs 0 <= minimum < maximum'
    stop 1
  endif
  if (pso_starts < 1 .or. pso_ssq_factor < 1D0 .or. pso_cluster_dist < 0D0) then
    write(*,*) 'paraest.cfg: pso_starts >= 1, pso_ssq_factor >= 1 and pso_cluster_dist >= 0 are required'
    stop 1
  endif

  open(unit=1, status='OLD', file='paraest.out', form='FORMATTED', position='APPEND')
  write(1,'(/,X,A,/,X,A,ES12.4,/,X,A,ES12.4)') 'Options read from paraest.cfg:', &
    'Residual of a failed point (PENALIZ):    ', PENALIZ, &
    'Residual replacing NaN/Inf (PENALIZNAN): ', PENALIZNAN
  if (any(range_kappa /= [0D0, 20D0])) &
    write(1,'(X,A,2F10.4)') 'Bounds of the association volume (range_kappa):', range_kappa
  if (bubble_solver /= 'tpf') write(1,'(X,2A)') 'Bubble points calculated with: ', trim(bubble_solver)
  if (optimizer /= 'lm') then
    write(1,'(X,2A)') 'Optimizer: ', trim(optimizer)
    write(1,'(X,A,I0,A,I0,A,I0)') 'PSO: particles = ', pso_particles, ', iterations = ', pso_iterations, &
      ', seed = ', pso_seed
    if (optimizer == 'pso+lm') write(1,'(X,A,I0,A,F6.2,A,F6.3)') 'LM starts: up to ', pso_starts, &
      ', SSQ <= ', pso_ssq_factor, ' x best, distance > ', pso_cluster_dist
    write(1,'(X,A,8(/,X,A,2F10.4))') "PSO search box (k', association volume and dc: their bounds):", &
      '  g*  (times the initial value): ', range_gstr, &
      '  Q   (times the initial value): ', range_q, &
      '  eps (times the initial value): ', range_eps, &
      '  Tc  (times the initial value): ', range_tc, &
      '  k*  (absolute):                ', range_kstr, &
      '  alpha (absolute):              ', range_alpha, &
      "  g'  (absolute):                ", range_gp, &
      '  g"  (absolute):                ', range_gdd
  endif
  close(unit=1)
endsubroutine read_options

endmodule config
