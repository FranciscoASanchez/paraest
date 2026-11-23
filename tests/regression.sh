#!/usr/bin/env bash
# =============================================================================
#  paraest regression tests
#
#  Each case is a directory tests/cases/<name>/ with:
#     paraest.dat        input file
#     reference/         expected output (paraest.out, graficos.csv, tables.csv)
#  FOR010-FOR015.DAT and fort.3 are taken from the case if it has them, otherwise
#  from the root directory. paraest.cfg (optional) is taken only from the case.
#
#  Each case runs in tests/tmp/<name>/ (all at once, in parallel) and its
#  output is compared with the reference. Numbers are compared with relative
#  tolerance TOL (default 1e-6), so that compiler rounding differences are not
#  taken for real changes.
#
#  Usage:   tests/regression.sh ./paraest.x [TOL] [case1 case2 ...]
#           (no cases: run them all)
#
#  With NUMBERS_ONLY=1 only the numbers of each line are compared and the
#  text is ignored (useful when only the wording of the output changed).
# =============================================================================
set -u
# "C" locale: otherwise, with a Spanish locale awk uses decimal commas and
# reads "5.537" as 5.
export LC_ALL=C
EXE=$(realpath "${1:-./paraest.x}")
TOL=${2:-1e-6}
shift $(( $# >= 2 ? 2 : $# ))
ROOT=$(cd "$(dirname "$0")/.." && pwd)
CASES=$ROOT/tests/cases
TMP=$ROOT/tests/tmp
LIMIT=${LIMIT:-1800}              # maximum seconds per case
NUMBERS_ONLY=${NUMBERS_ONLY:-0}

if [ $# -gt 0 ]; then list=("$@"); else list=($(ls "$CASES")); fi

# Compares two files word by word; numbers with relative tolerance.
compare() {   # compare <reference> <new> <name to show>
    awk -v tol="$TOL" -v name="$3" -v numonly="$NUMBERS_ONLY" '
    function isnum(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([EeDd][-+]?[0-9]+)?$/ }
    function num(s)   { gsub(/[Dd]/, "E", s); return s + 0 }
    function abs(x)   { return x < 0 ? -x : x }
    # Splits a line into words; with numonly, keeps only the numbers.
    function words(line, w,    n, k, i, t) {
        n = split(line, t, /[ \t,;()=]+/); k = 0
        for (i = 1; i <= n; i++) if (!numonly || isnum(t[i])) w[++k] = t[i]
        return k
    }
    # Non-ASCII characters (accents) are replaced by "?": the output mixes
    # program text (UTF-8) with titles copied from paraest.dat (Latin-1).
    # Lines with only blanks count as empty: write(*,*) with no data writes
    # " " with ifx and "" with gfortran.
    NR == FNR { gsub(/\r/, ""); gsub(/[\200-\377]+/, "?"); sub(/^[ \t]+$/, ""); ref[FNR] = $0; nref = FNR; next }
    {
        gsub(/\r/, ""); gsub(/[\200-\377]+/, "?"); sub(/^[ \t]+$/, ""); line = FNR
        if (!(line in ref)) { bad++; if (bad <= 10) print "    " name ": extra line " line ": " $0; next }
        delete a; delete b
        n1 = words(ref[line], a); n2 = words($0, b)
        differ = (n1 != n2)
        for (i = 1; i <= n1 && !differ; i++) {
            if (a[i] == b[i]) continue
            if (isnum(a[i]) && isnum(b[i])) {
                x = num(a[i]); y = num(b[i]); m = abs(x) > abs(y) ? abs(x) : abs(y)
                if (abs(x - y) > tol * (m > 1e-30 ? m : 1)) differ = 1
            } else differ = 1
        }
        if (differ) {
            bad++
            if (bad <= 10) { print "    " name ": line " line; print "       ref: " ref[line]; print "       new: " $0 }
        }
    }
    END {
        if (FNR < nref) { print "    " name ": " nref - FNR " lines missing"; bad++ }
        if (bad > 10) print "    " name ": ... " bad " different lines in total"
        exit (bad > 0)
    }' "$1" "$2"
}

# --- Run all the cases in parallel -------------------------------------------
echo "Running ${#list[@]} case(s) with $EXE ..."
for case in "${list[@]}"; do
    d=$TMP/$case
    rm -rf "$d" && mkdir -p "$d"
    cp "$CASES/$case/paraest.dat" "$d/"
    for f in FOR010.DAT FOR011.DAT FOR012.DAT FOR013.DAT FOR014.DAT FOR015.DAT fort.3; do
        if [ -e "$CASES/$case/$f" ]; then cp "$CASES/$case/$f" "$d/"; else cp "$ROOT/$f" "$d/"; fi
    done
    if [ -e "$CASES/$case/paraest.cfg" ]; then cp "$CASES/$case/paraest.cfg" "$d/"; fi
    ( cd "$d" && timeout "$LIMIT" "$EXE" < /dev/null > screen.txt 2>&1; echo $? > .status ) &
done
wait

# --- Compare -------------------------------------------------------------------
status=0
for case in "${list[@]}"; do
    d=$TMP/$case
    rc=$(cat "$d/.status")
    if [ "$rc" != 0 ]; then
        echo "  ERROR   $case: the program ended with code $rc (see $d/screen.txt)"
        status=1; continue
    fi
    # Details are collected first and printed after the case line.
    bad=""; details=""
    for ref in "$CASES/$case"/reference/*; do
        f=$(basename "$ref")
        if [ ! -e "$d/$f" ]; then
            details="$details    $f: not generated"$'\n'; bad="$bad $f"
        elif ! out=$(compare "$ref" "$d/$f" "$f"); then
            details="$details$out"$'\n'; bad="$bad $f"
        fi
    done
    if [ -z "$bad" ]; then
        echo "  OK      $case"
    else
        echo "  DIFFERENT $case:$bad"; printf "%s" "$details"; status=1
    fi
done
[ $status -eq 0 ] && echo "Regression OK (relative tolerance $TOL)" || echo "Regression FAILED"
exit $status
