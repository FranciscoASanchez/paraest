paraest — Parameter fitting of thermodynamic models (GCA-EoS, RK-PR)
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

Model and parameters in paraest.dat (since 2026-11-26)
------------------------------------------------------
The line after the weights is "NC NSIG" (before: "NC NG NGA NSIG"). After
the vapor-pressure data comes the model, followed by what ReadParameters
of thermopack reads for it, as in the flash programs:

    model                               1: GCA, 4: RK-PR
  model 1:
    NG  NGA
    (PARGR block, below)
    (PARMOL block, below)
  model 4: what readcomp of RKPR reads:
    NComb  Tdep                         0: van der Waals; Tdep 1: k(T)
    name  j                             for each component; j = 0: next line
    Tc Pc omega vc Zrat  (j = 0)          is Tc/K Pc/bar omega vc/(L/mol)
    ac b delta1 k        (j = 1)          Zrat; j = 1: ac/(bar L2/mol2)
    M c A B C D E F G                     b/(L/mol) delta1 k
    NownIJ
    i j l_ij k_ij k'_ij T*_ij           NownIJ lines (Tdep = 0: i j l_ij k_ij)
  Then, for any model, the melting block (below).

With RK-PR only data calculated through THERMO can be used: pure-component
isofugacity (IOBJSAT = 50) and vapor pressure (IOBJSAT 80, 81 with liquid
density), phase equilibrium (IHEN = 0; bubble points with the 'tpf' solver)
and solid solubility (IHEN = 2). Other data stop the program.
The vapor pressure is calculated by psat_pure (src/equilibrium/saturation.f90,
from GCA-cape-open, with THERMO in place of lnphi_TP); the experimental
value is the initial guess. The GCA still uses its PSAT.
Example: tests/cases/meoh-water-rkpr-psat (k of water fitted to its vapor
pressure, with the methanol + water data of meoh-water-rkpr).

Pure-component data lines (any model): IOBJSAT = 50 and 80: "T  P  wP";
IOBJSAT = 81 (since 2026-11-27): "T  P  v  wP  wv", with the molar volume
v in cm3/mol (that of the solid below Tfus: v_S = v_L - DVfus) and the
weights of P and v: residual wP (Pcalc - Pexp) and wv (vcalc - vexp), both
times W1, combined as the square root of their sum of squares when both
are nonzero. With wP = 0, P is only the initial guess (e.g. an extrapolated
sublimation pressure). Earlier, IOBJSAT 81 read "T  P  rho  w" (rho in
mol/cm3, unweighted).
Example: tests/cases/rkpr-pure-fit (vanillin and ethyl levulinate: delta1 = 1 + sqrt(2)
fixed, k fitted to the vapor pressure; k and delta1 fitted to the vapor
pressure and the molar volume, of the solid for vanillin and of the liquid
for ethyl levulinate; and the published parameters of ethyl levulinate).
tests/cases/rkpr-pure-fit-tcpc: the same, with the estimated Tc and Pc of
vanillin, fitted too (LM bounds lm_range_tc 0.9-1.1, lm_range_pc 0.75-1.25
in paraest.cfg). With the vapor pressure only (VAN-PR), Tc, Pc and
k are not determined: Pc ends at its bound.

RK-PR fitted parameters (since 2026-11-27): after the data sets, in the
format of those of the GCA:

    NTc NPc Ndelta Nk Nkij Nlij
    components whose Tc is fitted       (only if NTc > 0; the same, one line
    ... Pc, delta1, k                    each, for Pc, delta1 and k)
    i j  i j ...                        pairs of k_ij               } only if
    T dependence of each pair:          0: k_inf; 1: k_inf, k';     } Nkij > 0
                                        2: k_inf, k', T*            }
    number of identical pairs of each   (pairs that copy its values)}
    i j  i j ...                        one line per pair with identical ones
    i j  i j ...                        pairs of l_ij               } only if
    number of identical pairs of each                               } Nlij > 0
    i j  i j ...                        one line per pair with identical ones
    initial values                      Tc/K, Pc/bar, delta1, k (in the order
                                        of the lists), then k_inf [k' [T*/K]]
                                        of each pair, then l_ij

  k_ij(T) = k_inf + k' exp(-T/T*). The pure-component parameters fitted
  are Tc, Pc and delta1 (and k), not ac and b: these follow from them, as
  readcomp does with input j = 0 (with delta1 fixed, ac and b only move the
  critical point). Tc and Pc are those of the EoS (with j = 1 input, those
  calculated from ac and b, written in paraest.out). An initial value 0 of
  Tc, Pc, delta1, k or T* takes the value read with the model. A pair
  fitted with T dependence 0 has k' = 0, whatever the input says. Fitting k'
  turns on the temperature dependence (Tdep = 1) if the input did not have
  it. The binary parameters are symmetric: (j,i) is the same parameter as
  (i,j), it need not be listed as identical, and listing both as fitted
  parameters is an error.
  Bounds (Levenberg-Marquardt) and PSO box:
    Tc, Pc       lm_range_tc, lm_range_pc (default 0.5-2 times the initial value); PSO range_tc, range_pc
    delta1       range_delta1 (default 1-6), both
    k            0.1-10; PSO range_k (times the initial value)
    k_inf, l_ij  1 - (30 to 0.2), i.e. -29 to 0.8, as 1 - k* of the GCA;
                 PSO 1 - range_kstr (-0.2 to 0.2)
    k'           -1 to 1; PSO range_kp (default -0.15 to 0.15)
    T*           0.1 min(Tc_i, Tc_j) to 2 max(Tc_i, Tc_j); PSO 0.5-2 times
                 min(Tc_i, Tc_j)
  An initial value outside its bounds widens them (with a warning).
  Example: tests/cases/meoh-water-rkpr-fit.

Files in the earlier format are converted with
    tools/convert_input.py --model file.dat new.dat
(the conversions below already write the new format).

GCA parameters (since 2026-11-23, as in the flash programs)
-----------------------------------------------------------
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
symmetric: fitting (i,j) also sets (j,i) (for the association parameters,
(k,i,l,j) also sets (l,j,k,i)), so (j,i) need not be listed as identical
(older files that do it give the same results).

Old files (all the parameters written out, groups numbered 1..NG) are
converted with
    tools/convert_input.py old.dat new.dat [directory with FOR010-FOR015]
It finds each group in FOR010 by (Q, T*), and each associating group by the
name in the comment of its composition line ("!asoc AOH": short name or
description in FOR014) or else by its association energies and volumes in
FOR015; it writes as own values what differs from the database, and tells
which groups it could not find (they take free numbers, with own values).

Pure-component properties (since 2026-11-25) are read by PARMOL of the GCA
package, right after the PARGR block, and are followed by the melting block
(independent of the model):

    nUnit                               1: pressures in atm, 0: in MPa
  for each component:
    name
    Tc  Pc  omega  Tsat  Psat  dc       dc = 0: estimated from Tc, Pc; Tsat > 0:
                                        dc is fitted to (Tsat, Psat) at every
                                        evaluation (unless dc is a fitted
                                        parameter); Tsat = 0 and omega > 0:
                                        PARMOL fits it to a point from omega
    M  c  A  B  C  D  E  F  G           molar mass, Peneloux, Passut-Danner
  melting block, for each component:
    Tfus  Pfus  DHfus  DVfus  ACpfus  BCpfus  NTSS
                                        K, atm, J/mol, cm3/mol, J/(mol K),
                                        J/(mol K2); DCpfus = ACpfus + BCpfus*T
    Ttr  DHtr                           NTSS lines (solid-solid transitions)
  The transitions go from the highest Ttr down. The solubility of solids
  (GCSLE) uses all of it; the sublimation pressure of a pure component
  (vapor-pressure data below Tfus) only ACpfus and the first transition
  (it warns otherwise).

Files with the earlier component line "Tc Pc Tsat Psat dc Tfu DHfu DCpfu
TtrSS DHss DVfu Pfu" are converted with
    tools/convert_input.py --components file.dat new.dat
(tools/convert_input.py old.dat new.dat already writes the new format).

Data sets: the title line is read with FORMAT(40A2) for 65 items, so it
takes two records: leave a blank line after it.

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
  used by Levenberg-Marquardt (k* 0.2-30, k' +-0.3, alpha +-29, g' +-2),
  widened to include the initial value; association volume and dc between
  their physical bounds:
  range_gstr      g*                    default 0.5, 1.5
  range_q         Q                     default 0.5, 2.0
  range_eps       association energy    default 0.5, 2.0
  range_tc        Tc                    default 0.9, 1.1
  range_kstr      k* (absolute values)  default 0.8, 1.2
  range_alpha     alpha (absolute)      default -10.0, 10.0
  range_gp        g' (absolute values)  default -1.0, 1.0
  range_gdd       g" (absolute values)  default -1.0, 1.0
  range_kappa     association volume    default 0.0, 20.0 (cm3/mol): bounds
                  for both PSO and Levenberg-Marquardt, not only the box
  range_kp        k' (absolute; GCA and RK-PR)  default -0.15, 0.15
  RK-PR:
  range_pc        Pc (times initial)    default 0.8, 1.25
  range_k         k (times initial)     default 0.7, 1.3
  range_delta1    delta1                default 1.0, 6.0: bounds for both
                  PSO and Levenberg-Marquardt
  lm_range_tc     Tc (times initial)    default 0.5, 2.0: Levenberg-Marquardt
  lm_range_pc     Pc (times initial)    default 0.5, 2.0  bounds (not the box)
  bubble_solver   bubble points of binary VLE data (IFL = 1, 2):
                  'tpf'    Newton-Raphson in ln K and ln P (or ln T), adapted
                           from TPfractionFlash v0.5.3 (default)
                  'gcadb'  the old GCADB routine
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
    src/equilibrium/   equilibrium calculations (flash, bubble point, isofugacity...);
                       tpf_flash.f90: flash at given V/F (bubble points)
                       gc3_flash.f90, gcsle.f90: made from GC3 and GCSLE by
                       tools/port_gc3.py and tools/port_gcsle.py
    src/models/rkpr.f90, thermopack.f90: made by tools/port_rkpr.py
    (to update one of these, run its script with the new version of the
    original; do not edit the generated file)
    tests/             regression tests and reference outputs
    Archivos/          experimental data sets
    legacy/            old versions of the sources and makefiles

Versions are kept with git (git log, git tag), not in file names.
