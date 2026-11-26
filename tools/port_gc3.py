#!/usr/bin/env python3
"""Port the GC3 multiphase flash program into paraest as subroutine GC3Flash.

Usage:  python3 tools/port_gc3.py [GC3/src] [version]
        (defaults: ../GC3/src 1.5.3.1)

Writes src/equilibrium/gc3_flash.f90 (UTF-8, LF). Same recipe as the earlier
paraest copy (GC3 1.2.101.0-paraest): the main program becomes a subroutine
that takes the feed, T and P and returns the phases instead of reading
GCTHREE.DAT and writing datter.dat/Compres.dat; INDATA does not read the
parameters (paraest has already loaded the model). The other routines are
kept as in the original. Every change is checked: the script stops if a
piece of text it has to change is not found (e.g. after an update of GC3).
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT.parent / "GC3" / "src"
version = sys.argv[2] if len(sys.argv) > 2 else "1.5.3.1"


def replace(text, old, new, count=1):
    n = text.count(old)
    if n != count:
        sys.exit(f"port_gc3: expected {count} occurrence(s), found {n}:\n{old}")
    return text.replace(old, new)


def resub(text, pattern, new, count=1):
    out, n = re.subn(pattern, lambda m: new, text, flags=re.S | re.M)
    if n != count:
        sys.exit(f"port_gc3: expected {count} match(es), found {n}:\n{pattern}")
    return out


gc3 = (src / f"GC3-{version}.f90").read_bytes().decode("latin-1").replace("\r\n", "\n")

# ----------------------------------------------------------------------------
# Main program -> subroutine GC3Flash
# ----------------------------------------------------------------------------
gc3 = replace(gc3, "\tprogram GC3\n!\n\tIMPLICIT REAL*8(A-H,O-Z)\n",
              "\tsubroutine GC3Flash (NCOMP, TEMP, PRES, feed, xFrac, Zcomp, indFase, iErrTP)\n"
              "!\n\tuse thermo_model, only: model\n\tIMPLICIT REAL*8(A-H,O-Z)\n")
# GCTHREE.SAL on a unit that paraest does not use (2 is tables.csv);
# maximum number of phase splits, as in the earlier paraest copy.
gc3 = replace(gc3, "\tinteger, parameter      :: inputFile = 5, NCM = 30, NFM = 4, outputFile = 2\n",
              "\tinteger, parameter      :: inputFile = 5, NCM = 30, NFM = 4, outputFile = 33, maxIt = 50\n")
gc3 = replace(gc3, "\tcharacter, parameter    :: TAB = achar(9)\n",
              "\tcharacter, parameter    :: TAB = achar(9)\n"
              "\tdimension                  feed(NCM), xFrac(NCM,NFM), Zcomp(NFM), indFase(NFM)\n")
# only used by the debug output to Compres.dat (removed below)
gc3 = resub(gc3, r"^\tcommon /ZZZAs/[^\n]*\n", "")
gc3 = replace(gc3,
              "\tOPEN (UNIT = inputFile, FILE = 'GCTHREE.DAT', STATUS = 'OLD')\n"
              "\tOPEN (UNIT = IOUT, FILE='GCTHREE.SAL', STATUS = 'unknown',FORM = 'FORMATTED')\n"
              "\topen (unit = 1, file= 'datter.dat', status = 'unknown', form = 'formatted')\n"
              "\topen (unit = 9, file= 'Compres.dat', form = 'formatted')\n",
              "\tOPEN (UNIT = IOUT, FILE='GCTHREE.SAL', STATUS = 'unknown',FORM = 'FORMATTED')\n")
# program title (datter.dat and GCTHREE.SAL)
gc3 = resub(gc3, r"^\twrite \(1,500\)\n.*?^\twrite \(iOut,500\)\n(?=!-+\n)", "")
gc3 = replace(gc3, "\tREAD (5,*) model, NCOMP, solidPhase\n",
              "\tsolidPhase = 0\t! paraest: no solid phases (yet)\n\tiErrTP = 0\n\tnIt = 0\n")
# column titles of datter.dat
gc3 = resub(gc3, r"^\twrite\(1,'\(\"Phases\"[^\n]*\n\twrite\(1,15\)[^\n]*\n15\tformat[^\n]*\n", "")
gc3 = replace(gc3, "260\tREAD (5, *, ERR=250, END=250) (Z(I), I=1,N), T, P\n",
              "\tZ(:N) = feed(:N)\n\tT = TEMP\n\tP = PRES\n")
gc3 = replace(gc3, "\tIF (T.LE.0.) GOTO 260\n", "\tIF (T.LE.0.) GOTO 250\n")
# Compres.dat (format with <N>, an Intel extension); format 60 is also used
# for the screen and GCTHREE.SAL: the one of the earlier versions.
gc3 = resub(gc3, r"^\twrite \(9, 60\)[^\n]*\n60\tformat[^\n]*\n",
            "60\tFORMAT(1H1,///,' PT-Flash, T = ',F8.2,' K, P = ',F8.2,' atm, Z = ',F10.8)\n")
# Stable solution: return the phases (one phase: error) instead of writing
# them to datter.dat.
gc3 = resub(gc3, r"^\tnPunto = nPunto \+ 1\n.*?^171\tformat[^\n]*\n",
            "\tnPunto = nPunto + 1\n"
            "!\tparaest: the phases are returned; a single phase is an error.\n"
            "\tif (NF > 1) then\n"
            "\t  xFrac = zFrac\n"
            "\t  Zcomp = zVal\n"
            "\t  indFase = indF\n"
            "\telse\n"
            "\t  iErrTP = 1\n"
            "\tendif\n"
            "\treturn\n")
gc3 = replace(gc3, "190\tFORMAT (/,' System is still unstable, phase split performed '/)\n",
              "190\tFORMAT (/,' System is still unstable, phase split performed '/)\n"
              "\tnIt = nIt + 1\n"
              "\tif (nIt > maxIt) GOTO 250\n")
gc3 = replace(gc3, "\tIF (NF > 4) GOTO 260\n", "\tIF (NF > 4) GOTO 250\n")
gc3 = replace(gc3, "250\tCLOSE (UNIT = 5); close(unit = 1)\n!\n\tSTOP\n\tEND\n",
              "250\tiErrTP = 1\n260\treturn\n\tEND subroutine GC3Flash\n")

# ----------------------------------------------------------------------------
# INDATA: paraest has already loaded the model parameters.
# ----------------------------------------------------------------------------
gc3 = replace(gc3, "\tcall ReadParameters (inputFile, iOut, model, N)\n",
              "!\tparaest: the parameters are already loaded (ReadParameters not called)\n")

# ----------------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------------
code = "\n".join(l for l in gc3.split("\n") if not l.lstrip().startswith("!"))
for bad in (r"<N[-+0-9]*>", r"\bwrite\s*\(\s*(1|9)\s*,", r"\bstop\b", r"\bprogram\b"):
    m = re.search(bad, code, re.I)
    if m:
        sys.exit(f"port_gc3: '{m.group(0)}' left in the code")

out = ROOT / "src/equilibrium/gc3_flash.f90"
out.write_text(f"""!
!	paraest copy of GC3-{version}.f90 (GC3/src), made by tools/port_gc3.py.
!	Program GC3 turned into
!	  subroutine GC3Flash (NCOMP, TEMP, PRES, feed, xFrac, Zcomp, indFase, iErrTP)
!	  feed(NCOMP) at TEMP (K), PRES (atm) -> phases xFrac(:,j), Zcomp(j),
!	  indFase(j) (THERMO phase type), j = 1..NF; iErrTP = 1 when the result
!	  is a single phase, after more than maxIt phase splits or more than 4
!	  phases. Changes:
!	  - no GCTHREE.DAT, datter.dat or Compres.dat; GCTHREE.SAL on unit 33.
!	  - model from module thermo_model; no solid phases (solidPhase = 0).
!	  - INDATA does not call ReadParameters (paraest has loaded the model).
!	The other routines are those of the original.
!
""" + gc3, encoding="utf-8")
print(f"port_gc3: {out.relative_to(ROOT)} ({gc3.count(chr(10))} lines)")
