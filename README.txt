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
    ./paraest.x           reads paraest.dat and the database FOR010-FOR015.DAT
                          from the current directory

Model parameters in paraest.dat (since 2026-11-23, as in the flash programs)
-----------------------------------------------------------------------------
The attractive and associative parameters are read by PARGR of the GCA
package (v1.14): the groups are identified by their number in the database
(FOR010.DAT: attractive groups; FOR014.DAT: associating groups) and the
database values can be replaced by the user's ("own") values. The block goes
where the old parameter block was (after the vapor-pressure data, before the
pure-component properties); the old group-composition lines (NU, ONUASOC)
are gone.

    NownII  NownIJ                      own pure-group lines, own binary lines
    id n  id n ...  k*0                 one line per component: (group, number)
                                        pairs, 2*NG values in all (pad with
                                        zeros); the ids in increasing order
    id  Q  T*  g*  g'  g"  name         NownII lines
    i  j  k*  k'  alpha_ij  alpha_ji    NownIJ lines
  only if NGA > 0:
    NownII  NownIJ                      own associating groups, own pairs
    id n  id n ...  k*0                 one line per component, 2*NGA values
    id  sites  name                     NownII lines
    k  i  l  j  eps/R  kappa            NownIJ lines: site k of group i with
                                        site l of group j (i <= j)

A group that is not in the database: use the number of a group that the
system does not use and give all its values as own values. A group listed
with 0 groups in a component is part of the system without being in it.

The fitted parameters (the lines after the data sets) refer to the groups by
their position in the system, not by their database number: PARGR orders
the groups by increasing database number, so position 1 is the group with
the lowest number, and so on (attractive and associating groups separately).
The parameter tables of paraest.out show it in column Pos, next to the
database number (ID). The association parameters are (site, group position,
site, group position). k*, k', and the association energy and volume are
symmetric: fitting (i,j) also sets (j,i).

Old files (all the parameters written out, groups numbered 1..NG) are
converted with
    tools/convert_input.py old.dat new.dat [directory with FOR010-FOR015]
It finds each group in FOR010 by (Q, T*), writes as own values what differs
from the database, and tells which groups it could not find.

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
                  'pso+lm' particle swarm, then Levenberg-Marquardt from
                           PSO points in distinct regions (see below)
  pso_particles   default 35     } pso_particles*(pso_iterations+1)
  pso_iterations  default 50     } evaluations
  pso_seed        default 12345  (same seed, same run, with any compiler)
  With 'pso+lm', Levenberg-Marquardt starts from up to pso_starts PSO
  points that are good and far from each other:
  pso_starts       default 5    (1: only from the best point)
  pso_ssq_factor   default 2.0  candidates: SSQ <= this times the best one,
                                and at most one point lost more than it
  pso_cluster_dist default 0.1  distance between two points: RMS of the
                                differences divided by the width of the
                                PSO box of each parameter
  The candidates are taken from best to worst; one becomes a start if it is
  farther than pso_cluster_dist from every start already chosen. The LM
  results are grouped the same way (closer than pso_cluster_dist: the same
  optimum) and listed in paraest.out, best first (fewest points lost, then
  SSQ); the final parameters and statistics are those of the best.
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
