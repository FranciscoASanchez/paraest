#!/usr/bin/env python3
"""Port RK-PR and thermopack from TPfractionFlash into paraest.

Usage:  python3 tools/port_rkpr.py [TPfractionFlash/src] [vRKPR] [vTPack]
        (defaults: ../TPfractionFlash/src 0.3.6 0.4.1)

Writes src/models/rkpr.f90 and src/models/thermopack.f90 (UTF-8, LF).
The routines are kept as in the originals; only the changes listed in the
headers are applied, each one checked (the script stops if a piece of text
it has to change is not found, e.g. after an update of the originals).
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT.parent / "TPfractionFlash" / "src"
vRKPR = sys.argv[2] if len(sys.argv) > 2 else "0.3.6"
vTPack = sys.argv[3] if len(sys.argv) > 3 else "0.4.1"


def load(name):
    return (src / name).read_bytes().decode("latin-1").replace("\r\n", "\n")


def replace(text, old, new, count=1):
    n = text.count(old)
    if n != count:
        sys.exit(f"port_rkpr: expected {count} occurrence(s), found {n}:\n{old}")
    return text.replace(old, new)


def routine(text, kind, name):
    """Text of one routine, from its first line to its end line (inclusive),
    preceded by the comment block right above it."""
    m = re.search(rf"^[ \t]*(?:double precision[ \t]+)?{kind}[ \t]+{name}[ \t]*\(.*?^[ \t]*end[ \t]*{kind}\b[^\n]*\n",
                  text, re.S | re.M | re.I)
    if not m:
        sys.exit(f"port_rkpr: {kind} {name} not found")
    start = m.start()
    # comment lines immediately above
    lines = text[:start].split("\n")[:-1]
    k = len(lines)
    while k > 0 and lines[k - 1].lstrip().startswith("!"):
        k -= 1
    head = "\n".join(lines[k:])
    return (head + "\n" if head else "") + m.group(0)


def write(path, header, body):
    path.write_text(header + body, encoding="utf-8")
    print(f"port_rkpr: {path.relative_to(ROOT)} ({body.count(chr(10))} lines)")


# ----------------------------------------------------------------------------
# RK-PR: whole file
# ----------------------------------------------------------------------------
rk = load(f"RKPR-v{vRKPR}.f90")
# gfortran: no variable format expressions <NC> (Intel extension).
rk = replace(rk, "18X, <NC>A10)", "18X, 30A10)", 4)
rk = replace(rk, "'|',<NC>(F8.4, 2X),'|')", "'|',30(F8.4, 2X),'|')")
# readcomp: thermoName is not declared (implicit REAL) and never used;
# gfortran rejects assigning it a string.
rk = replace(rk, '\tthermoName = "RKPR EOS"\n', "")
if "<NC>" in rk:
    sys.exit("port_rkpr: <NC> left in RKPR")

write(ROOT / "src/models/rkpr.f90", f"""!
!	paraest copy of RKPR-v{vRKPR}.f90 (TPfractionFlash/src), made by
!	tools/port_rkpr.py, with:
!	  - gfortran: fixed repeat counts (30) instead of <NC> (Intel extension)
!	    in the formats of the binary-parameter matrices of readcomp.
!	  - readcomp: unused assignment thermoName = "RKPR EOS" removed
!	    (undeclared, implicitly REAL; gfortran rejects it).
!
""", rk)

# ----------------------------------------------------------------------------
# thermopack: only what GCA (model 1) and RK-PR (model 4) need
# ----------------------------------------------------------------------------
tp = load(f"thermopack-v{vTPack}.f90")
keep = [("subroutine", "THERMO"), ("subroutine", "ReadParameters"),
        ("SUBROUTINE", "GCA_EoS"), ("SUBROUTINE", "RKPR_EoS"),
        ("subroutine", "Cp_GI"), ("SUBROUTINE", "TERMO"),
        ("SUBROUTINE", "ArVnder"), ("SUBROUTINE", "VCALC"),
        ("SUBROUTINE", "Bcalc"), ("subroutine", "screening"),
        ("subroutine", "deltapressure"),
        ("function", "fScld"), ("function", "fDscld")]
parts = {name: routine(tp, kind, name) for kind, name in keep}

# THERMO: models 1 and 4 only; any other model returns ic = -99.
t = parts["THERMO"]
t = re.sub(r"\n\telseif \(model == 2 \.OR\. model == 3\) then.*?(?=\n\telseif \(model == 4\) then)",
           "", t, flags=re.S)
t = re.sub(r"\n\telseif \(model == 5\) then.*?(?=\n\tendif\n)",
           "\n\telse\n\n\t  ic = -99\t! model not available in paraest\n", t, flags=re.S)
for m in ("SRKPR_EoS", "PCSAFT_EoS", "SanchezLacombe_EoS", "MTC_EoS", "gamma_phi"):
    if m in t:
        sys.exit(f"port_rkpr: {m} left in THERMO")
parts["THERMO"] = t

# GCA_EoS: n_deriv and temp_deriv are initialized in their declaration
# (implicit SAVE), so after one call with calc = 5/6/9 (or 4/6/7/8/9) every
# later call computed second (or temperature) derivatives too. Reset them at
# each call, as RKPR_EoS does.
g = parts["GCA_EoS"]
g = replace(g, "\tcommon /GCA/                    NG, NGA, NST\n!\n",
            "\tcommon /GCA/                    NG, NGA, NST\n!\n\tn_deriv = 1\n\ttemp_deriv = 0\n")
parts["GCA_EoS"] = g

# ReadParameters: models 1 and 4 only. The COMMON /cubicEq/ with initialized
# variables is not valid outside BLOCK DATA (gfortran) and only served SRK/PR.
r = parts["ReadParameters"]
r = replace(r, "\tinteger                 :: iCEq = 0, NONHY = 0\n", "")
r = replace(r, "\tcommon /cubicEq/           iCEq, NONHY\n", "")
r = re.sub(r"\n\telseif \(model == 2 \.OR\. model == 3\) then.*?(?=\n\telseif \(model == 4\) then)",
           "", r, flags=re.S)
r = re.sub(r"\n\telseif \(model == 5\) then.*?(?=\n\tendif\n)",
           "\n\telse\n\n\t  write (outputFile, *) \"Model\", model, \"is not available in paraest\"\n\t  stop\n", r, flags=re.S)
parts["ReadParameters"] = r

# ArVnder: RK-PR only (GCA does not use VCALC/TERMO).
a = parts["ArVnder"]
a = re.sub(r"\n \tELSE IF \(n_model == 5 \) THEN.*?(?=\n\tEND IF\n)", "", a, flags=re.S)
for m in ("HelmPCSAFT", "Helmhotz_SL", "Helmholtz_MTC"):
    if re.search(rf"^[^!\n]*call +{m}", a, re.M | re.I):
        sys.exit(f"port_rkpr: {m} left in ArVnder")
parts["ArVnder"] = a

# Bcalc: COMMON /MOL/ declared with 2 components; with more, D(i) was read
# from the wrong place. PC-SAFT and MTC branches removed.
b = parts["Bcalc"]
b = replace(b, "COMMON /MOL/     DC(2),D(2),DT(2),HA(2),HB(2)",
            "COMMON /MOL/     DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)")
b = re.sub(r"\n\telseif \(n_model == 5\) then.*?(?=\n\telse\n)", "", b, flags=re.S)
for m in ("DIAMET", "MTC_covolume"):
    if re.search(rf"^[^!\n]*call +{m}", b, re.M | re.I):
        sys.exit(f"port_rkpr: {m} left in Bcalc")
parts["Bcalc"] = b

# VCALC: (1) a guess volume V0 <= B (zeta = B/V0 >= 1; e.g. the volume of
# another composition, or NaN) made ArVnder take the log of a negative number
# and, with NaN, (2) the loop never ended: every 30 iterations it screened
# again. A starting zeta outside (0, 1) is replaced by 0.5, and after 300
# iterations in all V is NaN (the caller's result is NaN: in paraest a point
# that could not be calculated).
v = parts["VCALC"]
v = replace(v, "\tIMPLICIT DOUBLE PRECISION (A-H,O-Z)\n",
            "\tuse, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan\n"
            "\tIMPLICIT DOUBLE PRECISION (A-H,O-Z)\n")
v = replace(v, "\tinteger                             :: n_deriv, NC_common, t_deriv\n",
            "\tinteger                             :: n_deriv, NC_common, t_deriv\n"
            "\tinteger                             :: iterTot\n")
v = replace(v, "\t  ZETA = B/V0\n\t  IC = ITYP\n\t  \n\tENDIF\n",
            "\t  ZETA = B/V0\n\t  IC = ITYP\n\t  \n\tENDIF\n"
            "\tIF (.NOT. (ZETA > 0D0 .AND. ZETA < 1D0)) ZETA = .5D0\t! paraest: invalid guess (V0 <= B, NaN)\n"
            "\titerTot = 0\n")
v = replace(v, "\t  V = CPV/ZETA\n\t  ITER = ITER + 1\n",
            "\t  V = CPV/ZETA\n\t  ITER = ITER + 1\n"
            "\t  iterTot = iterTot + 1\n"
            "\t  if (iterTot > 300) then\t! paraest: no convergence\n"
            "\t    V = ieee_value(V, ieee_quiet_nan)\n"
            "\t    IC = ITYP\n"
            "\t    return\n"
            "\t  endif\n")
parts["VCALC"] = v

body = "\n!" + "-" * 90 + "\n".join([""] + [parts[n] for _, n in keep])
write(ROOT / "src/models/thermopack.f90", f"""!
!	paraest copy of thermopack-v{vTPack}.f90 (TPfractionFlash/src), made by
!	tools/port_rkpr.py. Only the routines that GCA (model 1) and RK-PR
!	(model 4) need: THERMO, ReadParameters, GCA_EoS, RKPR_EoS, Cp_GI, TERMO,
!	ArVnder, VCALC, Bcalc, screening, deltapressure, and the external
!	functions fScld, fDscld used by screening (not paraest's module scaling,
!	whose functions are module procedures). Changes:
!	  - THERMO, ReadParameters, ArVnder, Bcalc: branches of the other models
!	    removed (THERMO returns ic = -99 for them; ReadParameters stops).
!	  - ReadParameters: COMMON /cubicEq/ removed (initialized variables in a
!	    COMMON outside BLOCK DATA; it only served SRK/PR).
!	  - Bcalc: COMMON /MOL/ with NCM components instead of 2.
!	  - GCA_EoS: n_deriv and temp_deriv reset at each call (initialized in
!	    the declaration, they kept the value of a previous call).
!	  - VCALC: a starting zeta = B/V0 outside (0, 1) is replaced by 0.5, and
!	    after 300 iterations it returns V = NaN (it looped forever).
!
""", body)
