#!/usr/bin/env bash
# =============================================================================
#  Prueba de regresión de paraest
#
#  Corre el ejecutable con los datos de entrada del proyecto en una carpeta
#  temporal (pruebas/tmp/) y compara las salidas con pruebas/referencia/.
#
#  Los números se comparan con tolerancia relativa TOL (por defecto 1e-6),
#  para no confundir diferencias de redondeo del compilador con cambios reales.
#
#  Uso:   pruebas/regresion.sh ./paraest.x  [TOL]
# =============================================================================
set -u
# Locale "C": si no, con es_AR awk usa coma decimal y lee "5.537" como 5.
export LC_ALL=C
EXE=$(realpath "${1:-./paraest.x}")
TOL=${2:-1e-6}
RAIZ=$(cd "$(dirname "$0")/.." && pwd)
REF=$RAIZ/pruebas/referencia
TMP=$RAIZ/pruebas/tmp

rm -rf "$TMP" && mkdir -p "$TMP"
cp "$RAIZ"/paraest.dat "$RAIZ"/FOR01?.DAT "$RAIZ"/fort.3 "$TMP"/
cd "$TMP" || exit 2

echo "Corriendo $EXE en $TMP ..."
if ! "$EXE" < /dev/null > pantalla.txt 2>&1; then
    echo "ERROR: el programa terminó con error. Ver $TMP/pantalla.txt"
    exit 1
fi

# Compara dos archivos palabra por palabra; los números con tolerancia relativa.
compara() {
    awk -v tol="$TOL" -v nombre="$1" '
    function esnum(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([EeDd][-+]?[0-9]+)?$/ }
    function num(s)   { gsub(/[Dd]/, "E", s); return s + 0 }
    function abs(x)   { return x < 0 ? -x : x }
    NR == FNR { gsub(/\r/, ""); ref[FNR] = $0; nref = FNR; next }
    {
        gsub(/\r/, ""); linea = FNR
        if (!(linea in ref)) { print nombre ": línea " linea " sobra: " $0; malas++; next }
        n1 = split(ref[linea], a, /[ \t,;()]+/); n2 = split($0, b, /[ \t,;()]+/)
        distinta = (n1 != n2)
        for (i = 1; i <= n1 && !distinta; i++) {
            if (a[i] == b[i]) continue
            if (esnum(a[i]) && esnum(b[i])) {
                x = num(a[i]); y = num(b[i]); m = abs(x) > abs(y) ? abs(x) : abs(y)
                if (abs(x - y) > tol * (m > 1e-30 ? m : 1)) distinta = 1
            } else distinta = 1
        }
        if (distinta) {
            malas++
            if (malas <= 10) { print nombre ": línea " linea; print "   ref: " ref[linea]; print "   new: " $0 }
        }
    }
    END {
        if (FNR < nref) { print nombre ": faltan " nref - FNR " líneas"; malas++ }
        if (malas > 10) print nombre ": ... " malas " líneas distintas en total"
        exit (malas > 0)
    }' "$REF/$1" "$TMP/$1"
}

estado=0
for f in paraest.out graficos.csv tables.csv; do
    if compara "$f"; then echo "  OK      $f"; else echo "  DISTINTO $f"; estado=1; fi
done
[ $estado -eq 0 ] && echo "Regresión OK (tolerancia relativa $TOL)" || echo "Regresión FALLÓ"
exit $estado
