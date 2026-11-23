!--------------------------------------------------------------------------------------------------------------------------
!
!  Module pso_optimizer: particle swarm optimization (PSO), global best.
!
!  Based on the PSO of 2019 (DriverPSO.for, "version 0"): inertia W = 1,
!  C1 = C2 = 2, and a maximum velocity that shrinks 10 % per iteration,
!  which is what makes the swarm converge. Corrected with respect to it:
!    - random numbers: one per particle and per dimension (it used the same
!      two for the whole swarm), from one seeded generator (it reseeded with
!      the clock at every iteration, and repeated sequences);
!    - the global best keeps the value of the personal best it comes from
!      (it took the current value of the particle);
!    - initial velocities between -VMAX and VMAX (they were all negative);
!    - the result is the best point found in the whole run (it returned the
!      best of the final swarm).
!  Also, particle 1 starts at the initial values, so the result is never
!  worse than the starting point.
!
!  Comparison of two points (feasibility rule): first the number of points
!  that could not be calculated (nBad, returned by OBJ), then the objective
!  function. Without it the swarm finds parameters where most points fail,
!  which give a small SSQ whatever the penalty (co2-phenols-fit: 118 points
!  lost, SSQ = 0). PSO only needs to know which point is better, so no
!  magnitude of the penalty has to be chosen.
!
!  The random numbers come from a generator of its own (Park-Miller), so a
!  run with the same seed gives the same result with any compiler.
!
!--------------------------------------------------------------------------------------------------------------------------
module pso_optimizer
  implicit none
  private
  public :: PSO

  double precision, parameter :: W = 1D0, C1 = 2D0, C2 = 2D0, alfaDecay = 0.9D0
  integer(8), parameter :: rngA = 16807_8, rngM = 2147483647_8
  integer(8) :: rngState = 1_8

contains

!--------------------------------------------------------------------------------------------------------------------------
!
!  Minimizes OBJ for x between lo and hi (a parameter with lo = hi stays
!  fixed). OBJ is a subroutine OBJ(N, x, f, nBad): f is the objective
!  function and nBad the number of points that could not be calculated.
!
!    x0               starting point (particle 1); brought inside the box
!    nPart, nIter     particles and iterations: nPart*(nIter + 1) evaluations
!    seed             seed of the random numbers (> 0)
!    xBest, fBest     best point found and its value
!    nBest            its number of points not calculated
!    nEval            number of evaluations
!
subroutine PSO(OBJ, N, x0, lo, hi, nPart, nIter, seed, xBest, fBest, nBest, nEval)
  external :: OBJ
  integer, intent(in) :: N, nPart, nIter, seed
  double precision, intent(in) :: x0(N), lo(N), hi(N)
  double precision, intent(out) :: xBest(N), fBest
  integer, intent(out) :: nBest, nEval

  double precision :: x(N,nPart), v(N,nPart), xP(N,nPart), fP(nPart), f, vMax(N), alfa
  integer :: nP(nPart), nBad, i, j, k

  rngState = mod(int(seed, 8), rngM - 1_8) + 1_8

  ! Initial swarm: particle 1 at x0, the others at random in the box
  x(:,1) = min(max(x0, lo), hi)
  do i = 2, nPart
    do j = 1, N
      x(j,i) = lo(j) + rnd()*(hi(j) - lo(j))
    enddo
  enddo
  alfa = 1D0
  vMax = alfa*(hi - lo)
  do i = 1, nPart
    do j = 1, N
      v(j,i) = (2D0*rnd() - 1D0)*vMax(j)
    enddo
  enddo

  fP = huge(1D0)
  nP = huge(1)
  fBest = huge(1D0)
  nBest = huge(1)
  xBest = x(:,1)
  nEval = 0
  do k = 0, nIter
    ! Evaluation, personal and global memory
    do i = 1, nPart
      call OBJ(N, x(:,i), f, nBad)
      nEval = nEval + 1
      if (f /= f) f = huge(1D0)          ! NaN: worst possible
      if (better(f, nBad, fP(i), nP(i))) then
        fP(i) = f
        nP(i) = nBad
        xP(:,i) = x(:,i)
      endif
      if (better(fP(i), nP(i), fBest, nBest)) then
        fBest = fP(i)
        nBest = nP(i)
        xBest = xP(:,i)
      endif
    enddo
    write(*,'(A,I4,A,ES14.6,A,I0)') ' PSO iteration ', k, ': best SSQ = ', fBest, ', points lost = ', nBest
    if (k == nIter) exit

    ! New velocities, limited to vMax (which shrinks), and positions
    vMax = alfa*(hi - lo)
    alfa = alfaDecay*alfa
    do i = 1, nPart
      do j = 1, N
        v(j,i) = W*v(j,i) + C1*rnd()*(xP(j,i) - x(j,i)) + C2*rnd()*(xBest(j) - x(j,i))
        v(j,i) = min(max(v(j,i), -vMax(j)), vMax(j))
        x(j,i) = min(max(x(j,i) + v(j,i), lo(j)), hi(j))
      enddo
    enddo
  enddo
endsubroutine PSO

!--------------------------------------------------------------------------------------------------------------------------
!
!  True if (f1, n1) is better than (f2, n2): fewer points not calculated,
!  or as many and a smaller objective function.
!
logical function better(f1, n1, f2, n2)
  double precision, intent(in) :: f1, f2
  integer, intent(in) :: n1, n2
  better = n1 < n2 .or. (n1 == n2 .and. f1 < f2)
endfunction better

!--------------------------------------------------------------------------------------------------------------------------
!
!  Uniform random number in (0, 1): Park-Miller minimal standard generator
!  (the product fits in 64 bits).
!
double precision function rnd()
  rngState = mod(rngA*rngState, rngM)
  rnd = dble(rngState)/dble(rngM)
endfunction rnd

endmodule pso_optimizer
