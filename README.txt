paraest — Parameter fitting of thermodynamic models (GCA-EoS)
=============================================================

Build
-----
    make                  with Intel ifx (default)
    make FC=gfortran      with GNU gfortran
    make debug            build with run-time checks -> debug.paraest.x
    make clean            remove objects and executables

Before using ifx, load the Intel environment:
    source /opt/intel/oneapi/setvars.sh

Run
---
    ./paraest.x           reads paraest.dat (and FOR010-013.DAT) from the current directory

Options (optional file paraest.cfg in the current directory, a namelist):
    &options
      optimizer = 'pso+lm'
      penaliz   = 1.0
    /
  penaliz         residual of a point whose calculation failed (default 0:
                  the point adds nothing to the objective function, so a fit
                  can lower it by driving points into failure)
  penaliznan      residual that replaces a NaN/Inf residual (default 1)
  optimizer       'lm'     Levenberg-Marquardt (ZXSSQ), default
                  'pso'    particle swarm only (no parameter statistics)
                  'pso+lm' particle swarm, then Levenberg-Marquardt from the
                           best point found
  pso_particles   default 35     } pso_particles*(pso_iterations+1)
  pso_iterations  default 50     } evaluations
  pso_seed        default 12345  (same seed, same run, with any compiler)
  PSO search box: g*, Q, eps and Tc between multiples of the initial value;
  k*, alpha, g' and g" between absolute values, narrower than the bounds
  used by Levenberg-Marquardt (k* 0.2-30, alpha +-29, g' +-2), widened to
  include the initial value; k', association volume and dc between their
  physical bounds:
  range_gstr      g*                    default 0.5, 1.5
  range_q         Q                     default 0.5, 2.0
  range_eps       association energy    default 0.5, 2.0
  range_tc        Tc                    default 0.9, 1.1
  range_kstr      k* (absolute values)  default 0.8, 1.2
  range_alpha     alpha (absolute)      default -10.0, 10.0
  range_gp        g' (absolute values)  default -1.0, 1.0
  range_gdd       g" (absolute values)  default -1.0, 1.0
  PSO compares two parameter sets first by the number of points that could
  not be calculated (beyond those that already fail at the initial values),
  then by SSQ: otherwise it finds sets where most points fail and SSQ is
  small. Levenberg-Marquardt afterwards uses penaliz as usual.
Source: src/main/config.f90. Without the file the defaults are used; with
it, the values used are written to paraest.out.

Regression tests
----------------
    make test

Runs every case in tests/cases/ in tests/tmp/ and compares paraest.out,
graficos.csv and tables.csv with the case's reference/ directory (relative
tolerance 1e-6). Run them after every change: if they end with "Regression
FAILED", the change altered the results. With NUMBERS_ONLY=1 only the
numbers are compared (useful when only the wording of the output changed):

    NUMBERS_ONLY=1 tests/regression.sh ./paraest.x

Layout
------
    src/main/          main program, input/output, FUNC, scaling, statistics
    src/optimization/  optimizers (ZXSSQ: Levenberg-Marquardt)
    src/models/        thermodynamic model (GCA package) and its link with paraest
    src/equilibrium/   equilibrium calculations (flash, bubble point, isofugacity...)
    tests/             regression tests and reference outputs
    Archivos/          experimental data sets
    legacy/            old versions of the sources and makefiles

Versions are kept with git (git log, git tag), not in file names.
