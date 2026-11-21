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
