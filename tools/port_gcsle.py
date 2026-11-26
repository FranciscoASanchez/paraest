#!/usr/bin/env python3
"""Port GCSLE (SolubSol23) into paraest as subroutine GCSLE.

Usage:  python3 tools/port_gcsle.py [SolubSol23/src] [version]
        (defaults: ../SolubSol23/src 2.0.0)

Writes src/equilibrium/gcsle.f90 (UTF-8, LF): the routines of the original
(pure_solute, root_from_newton, screening_SLE, normal, solute_isofugacity),
and, instead of program SLE and readSoluteProperties, subroutine GCSLE,
which solves one point as the loop of the program does (option 0: subcooled
liquid reference). fScld and fDscld are those of thermopack (as in the
original, which has them commented out). Every change is checked: the
script stops if a piece of text it has to change is not found.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT.parent / "SolubSol23" / "src"
version = sys.argv[2] if len(sys.argv) > 2 else "2.0.0"


def replace(text, old, new, count=1):
    n = text.count(old)
    if n != count:
        sys.exit(f"port_gcsle: expected {count} occurrence(s), found {n}:\n{old}")
    return text.replace(old, new)


def resub(text, pattern, new, count=1):
    out, n = re.subn(pattern, lambda m: new, text, flags=re.S | re.M)
    if n != count:
        sys.exit(f"port_gcsle: expected {count} match(es), found {n}:\n{pattern}")
    return out


sle = (src / f"GCSLE-{version}.f90").read_bytes().decode("latin-1").replace("\r\n", "\n")

# Module IOFILES and program SLE: replaced by subroutine GCSLE (below).
sle = resub(sle, r"\Amodule IOFILES\n.*?^endprogram SLE\n", "")
# readSoluteProperties: paraest reads the melting block itself.
sle = resub(sle, r"^!-+\n!\n!\tRead-write subroutine for solute properties\n!\n"
                 r"subroutine readSoluteProperties .*?^endsubroutine readSoluteProperties\n", "")
# root_from_newton: second_order is initialized in its declaration (implicit
# SAVE), so after the first point every call started with Newton steps and
# skipped the direct-substitution ones. Reset it at each call.
sle = replace(sle, "\tw_old = fScld (x(solid_ID), 1.D0, 0.D0, 1.D0)\n",
              "\tsecond_order = .FALSE.\n\tw_old = fScld (x(solid_ID), 1.D0, 0.D0, 1.D0)\n")
# gfortran: G edit descriptor needs a width.
sle = replace(sle, "write (*, '(4G)')", "write (*, '(4G15.7)')")

for bad in (r"\bprogram\b", r"\bIOFILES\b", r"\bthermoOut\b", r"readSoluteProperties"):
    m = re.search(bad, "\n".join(l for l in sle.split("\n") if not l.lstrip().startswith("!")), re.I)
    if m:
        sys.exit(f"port_gcsle: '{m.group(0)}' left in the code")

driver = """!--------------------------------------------------------------------------------------------------------------------------
!
!  paraest: solubility of the solid solute solid_ID at T (K), P (atm), in a
!  liquid whose solute-free composition is x_free(NC-1), as the loop of
!  program SLE does for one point, option 0 (subcooled liquid reference).
!  The program loops between Newton and screening until Newton converges;
!  here screening is tried once (it is deterministic, a second one would
!  give the same start). x0: initial estimate of the solute mole fraction
!  (<= 0: ideal solution). Melting data as read by paraest: J/mol, K, atm,
!  cm3/mol; NTSS solid-solid transitions (Tss, DHss).
!  Returns x (liquid composition), Z of the liquid and gamma, the activity
!  coefficient of the solute; complete miscibility (no solid at T, P): x = 1
!  for the solute. iErr = 1 if it did not converge.
!
subroutine GCSLE (NC, solid_ID, T, P, x_free, x0, DH_fus, ADCp_fus, BDCp_fus, Dv_fus, T_fus, P_fus, NTSS, Tss, DHss, &
                  Z, x, gamma, iErr)

	use thermo_model, only: model
	implicit none

	integer, parameter                   :: maxIteration = 10, NCM = 30
	real(8), parameter                   :: tolerance = 1.E-7

	integer, intent(in)                  :: NC, solid_ID, NTSS
	real(8), intent(in)                  :: T, P, x_free(NC-1), x0, DH_fus, ADCp_fus, BDCp_fus, Dv_fus, T_fus, P_fus, &
	                                        Tss(*), DHss(*)
	real(8), intent(out)                 :: Z, x(NCM), gamma
	integer, intent(out)                 :: iErr

	logical                              :: F_error, success
	integer                              :: B_ID, nu(2)
	real(8)                              :: lnFS, lnPhisL, PfusT, x_start, Y, ZsL, Tss_w(NTSS+1), DHss_w(NTSS+1)
	real(8)                              :: lnPhi(NCM)

	iErr = 0
	nu(1) = 1
	nu(2) = 0
	B_ID = solid_ID
	Tss_w = 0
	DHss_w = 0
	Tss_w(:NTSS) = Tss(:NTSS)
	DHss_w(:NTSS) = DHss(:NTSS)
	PfusT = P_fus
	x = 0
	lnPhi = 0
	Z = 0
	gamma = 0
	call pure_solute (0, solid_ID, B_ID, nu, NTSS, model, NC, T, P, T_fus, P_fus, DH_fus, &
	                  0D0, 0D0, ADCp_fus, BDCp_fus, Dv_fus, Tss_w, DHss_w, Y, lnPhisL, ZsL, PfusT)
	if ((T < T_fus .AND. Dv_fus == 0) .OR. (Dv_fus*(P - PfusT) > 0)) then

	  ! Solid fugacity divided by pressure (ln (fs/P))
	  lnFS = lnPhisL + Y
	  x_start = min(x0, 1D0)
	  if (x_start <= 0) x_start = dexp(Y)
	  x(solid_ID) = x_start
	  call normal (NC, solid_ID, x(:NC), x_free)
	  success = .FALSE.
	  F_error = .FALSE.
	  call root_from_newton (maxIteration, model, NC, 0, solid_ID, B_ID, nu, tolerance, lnFS, T, P, Z, x_free, x, lnPhi, &
	                         success, F_error)
	  if (.NOT. success .AND. .NOT. F_error) then
	    call screening_SLE (100, model, NC, 0, solid_ID, B_ID, nu, lnFS, x_start, T, P, x, x_free)
	    call root_from_newton (maxIteration, model, NC, 0, solid_ID, B_ID, nu, tolerance, lnFS, T, P, Z, x_free, x, lnPhi, &
	                           success, F_error)
	  endif
	  if (success) then
	    gamma = exp(lnPhi(solid_ID) - lnPhisL)
	  else
	    iErr = 1
	  endif

	else

	  ! complete miscibility
	  x(solid_ID) = 1
	  Z = ZsL
	  gamma = 1

	endif
	return
endsubroutine GCSLE
"""

out = ROOT / "src/equilibrium/gcsle.f90"
out.write_text(f"""!
!	paraest copy of GCSLE-{version}.f90 (SolubSol23/src), made by
!	tools/port_gcsle.py. Program SLE and readSoluteProperties replaced by
!	subroutine GCSLE (one point, option 0, at the end of this file).
!	Changes in the routines of the original:
!	  - root_from_newton: second_order reset at each call (initialized in
!	    the declaration, it kept the value of the previous point).
!	  - format (4G) -> (4G15.7) (gfortran needs the width).
!
""" + sle.lstrip("\n") + driver, encoding="utf-8")
print(f"port_gcsle: {out.relative_to(ROOT)}")
