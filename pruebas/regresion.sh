#!/usr/bin/env bash
# =============================================================================
#  Prueba de regresión de paraest
#
#  Cada caso es una carpeta pruebas/casos/<nombre>/ con:
#     paraest.dat        archivo de entrada
#     referencia/        salidas esperadas (paraest.out, graficos.csv, tables.csv)
#  Los FOR01?.DAT y fort.3 se toman del caso si los tiene; si no, de la raíz.
#
#  Cada caso se corre en pruebas/tmp/<nombre>/ (todos a la vez, en paralelo) y
#  sus salidas se comparan con la referencia. Los números se comparan con
#  tolerancia relativa TOL (por defecto 1e-6), para no confundir diferencias
#  de redondeo del compilador con cambios reales.
#
#  Uso:   pruebas/regresion.sh ./paraest.x [TOL] [caso1 caso2 ...]
#         (sin casos: corre todos)
# =============================================================================
set -u
# Locale "C": si no, con es_AR awk usa coma decimal y lee "5.537" como 5.
export LC_ALL=C
EXE=$(realpath "${1:-./paraest.x}")
TOL=${2:-1e-6}
shift $(( $# >= 2 ? 2 : $# ))
RAIZ=$(cd "$(dirname "$0")/.." && pwd)
CASOS=$RAIZ/pruebas/casos
TMP=$RAIZ/pruebas/tmp
LIMITE=${LIMITE:-1800}      # segundos máximos por caso

if [ $# -gt 0 ]; then lista=("$@"); else lista=($(ls "$CASOS")); fi

# Compara dos archivos palabra por palabra; los números con tolerancia relativa.
compara() {   # compara <referencia> <nuevo> <nombre para mostrar>
    awk -v tol="$TOL" -v nombre="$3" '
    function esnum(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([EeDd][-+]?[0-9]+)?$/ }
    function num(s)   { gsub(/[Dd]/, "E", s); return s + 0 }
    function abs(x)   { return x < 0 ? -x : x }
    # Los caracteres no-ASCII (acentos) se reemplazan por "?": la salida mezcla
    # textos del programa (UTF-8) con títulos copiados de paraest.dat (Latin-1).
    # Las líneas con solo blancos se toman como vacías: write(*,*) sin datos
    # escribe " " con ifx y "" con gfortran.
    NR == FNR { gsub(/\r/, ""); gsub(/[\200-\377]+/, "?"); sub(/^[ \t]+$/, ""); ref[FNR] = $0; nref = FNR; next }
    {
        gsub(/\r/, ""); gsub(/[\200-\377]+/, "?"); sub(/^[ \t]+$/, ""); linea = FNR
        if (!(linea in ref)) { malas++; if (malas <= 10) print "    " nombre ": línea " linea " sobra: " $0; next }
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
            if (malas <= 10) { print "    " nombre ": línea " linea; print "       ref: " ref[linea]; print "       new: " $0 }
        }
    }
    END {
        if (FNR < nref) { print "    " nombre ": faltan " nref - FNR " líneas"; malas++ }
        if (malas > 10) print "    " nombre ": ... " malas " líneas distintas en total"
        exit (malas > 0)
    }' "$1" "$2"
}

# --- Correr todos los casos en paralelo -------------------------------------
echo "Corriendo ${#lista[@]} caso(s) con $EXE ..."
for caso in "${lista[@]}"; do
    d=$TMP/$caso
    rm -rf "$d" && mkdir -p "$d"
    cp "$CASOS/$caso/paraest.dat" "$d/"
    for f in FOR010.DAT FOR011.DAT FOR012.DAT FOR013.DAT fort.3; do
        if [ -e "$CASOS/$caso/$f" ]; then cp "$CASOS/$caso/$f" "$d/"; else cp "$RAIZ/$f" "$d/"; fi
    done
    ( cd "$d" && timeout "$LIMITE" "$EXE" < /dev/null > pantalla.txt 2>&1; echo $? > .estado ) &
done
wait

# --- Comparar ----------------------------------------------------------------
estado=0
for caso in "${lista[@]}"; do
    d=$TMP/$caso
    rc=$(cat "$d/.estado")
    if [ "$rc" != 0 ]; then
        echo "  ERROR   $caso: el programa terminó con código $rc (ver $d/pantalla.txt)"
        estado=1; continue
    fi
    malos=""
    for ref in "$CASOS/$caso"/referencia/*; do
        f=$(basename "$ref")
        if [ ! -e "$d/$f" ]; then
            echo "    $f: no se generó"; malos="$malos $f"
        elif ! compara "$ref" "$d/$f" "$f"; then
            malos="$malos $f"
        fi
    done
    if [ -z "$malos" ]; then echo "  OK      $caso"; else echo "  DISTINTO $caso:$malos"; estado=1; fi
done
[ $estado -eq 0 ] && echo "Regresión OK (tolerancia relativa $TOL)" || echo "Regresión FALLÓ"
exit $estado
