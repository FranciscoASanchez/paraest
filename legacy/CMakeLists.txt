cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

project(SLE)

SET(LIBRARY_OUTPUT_PATH  "/opt/intel/lib/intel64/")

add_executable(paraest.x src/a-paraest-1.11.4.f90 src/A-DCADJUST-1.0.1.FOR src/A-dcadjust-fast-1.0.2.for src/A-gamma-1.0.1.for src/A-gc-eos-3.0.3.for src/A-isof-mezclas-1.0.1.for src/A-isof-puros.for src/A-psat-1.0.1.for src/A-rutin-mezclas.for src/A-Rutin-puros-1.02.for src/ASSOC-GC-EOSFev.FOR src/ZxssqFev-2.0.0.f90 src/GC3-1.2.101.0-paraest.for src/GCA-v1.9.27.0.f90 src/GCDBA-v1.3.3.paraest.for src/GCXS-v0.5.0.2.paraest.f90 src/GCSLE-v0.2.0.paraest.f90 src/semiImpCrit-v0.0.1.f90) 
