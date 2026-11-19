# =============================================================================
#  Makefile de paraest
#
#  Uso:
#     make                    compila optimizado con ifx  -> paraest.x
#     make FC=gfortran        compila optimizado con gfortran
#     make debug              compila con verificaciones -> debug.paraest.x
#     make prueba             compila, corre paraest.dat y compara con la referencia
#     make clean              borra objetos y ejecutables
#
#  Los objetos se guardan en obj/<compilador>/ (o debug/<compilador>/), así
#  ifx y gfortran no se mezclan.
# =============================================================================

# make trae FC=f77 predefinido; si nadie eligió compilador, usamos ifx.
ifeq ($(origin FC),default)
  FC = ifx
endif

# --- Archivos fuente, por capa ----------------------------------------------
PRINCIPAL    = src/principal/utilidades.f90              \
               src/principal/paraest.f90

OPTIMIZACION = src/optimizacion/zxssq.f90

MODELOS      = src/modelos/gca.f90                       \
               src/modelos/gca_assoc.f90                 \
               src/modelos/gca_assoc_fev.f90             \
               src/modelos/gca_assoc_rutin_puros.f90     \
               src/modelos/gca_assoc_rutin_mezclas.f90   \
               src/modelos/gca_assoc_dcadjust.f90        \
               src/modelos/gca_assoc_dcadjust_fast.f90   \
               src/modelos/gca_assoc_gamma.f90

EQUILIBRIO   = src/equilibrio/gc3_flash.f90              \
               src/equilibrio/gcdba.f90                  \
               src/equilibrio/gcxs.f90                   \
               src/equilibrio/gcsle.f90                  \
               src/equilibrio/pcrit_semiimplicito.f90    \
               src/equilibrio/psat_assoc.f90             \
               src/equilibrio/isof_puros_assoc.f90       \
               src/equilibrio/isof_mezclas_assoc.f90

# Orden de compilación: primero lo que define módulos, al final el programa.
FUENTES = $(OPTIMIZACION) $(MODELOS) $(EQUILIBRIO) $(PRINCIPAL)

# --- Opciones de cada compilador --------------------------------------------
ifeq ($(FC),gfortran)
  OPT   = -O2
  DEBUG = -g -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall
  COMUN = -std=legacy -ffree-line-length-none
  MODS  = -J
else
  OPT   = -O3
  DEBUG = -g -O0 -traceback -debug extended -fpe0 -check all -warn
  COMUN = -diag-disable=remark
  MODS  = -module
endif

# --- Derivados (no hace falta tocar de acá para abajo) ----------------------
COMP     = $(notdir $(FC))
OBJDIR   = obj/$(COMP)
DBGDIR   = debug/$(COMP)
objs     = $(addprefix $(1)/,$(addsuffix .o,$(basename $(notdir $(FUENTES)))))
OBJETOS  = $(call objs,$(OBJDIR))
OBJDEBUG = $(call objs,$(DBGDIR))

vpath %.f90 src/principal src/optimizacion src/modelos src/equilibrio

.PHONY: all debug prueba clean
all: paraest.x

paraest.x: $(OBJETOS)
	$(FC) $(OPT) -o $@ $^

debug: debug.paraest.x
debug.paraest.x: $(OBJDEBUG)
	$(FC) $(DEBUG) -o $@ $^

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(COMUN) $(OPT) $(MODS) $(OBJDIR) -c $< -o $@
$(DBGDIR)/%.o: %.f90 | $(DBGDIR)
	$(FC) $(COMUN) $(DEBUG) $(MODS) $(DBGDIR) -c $< -o $@

$(OBJDIR) $(DBGDIR):
	mkdir -p $@

# Prueba de regresión: corre el caso de paraest.dat en pruebas/tmp/
# y compara la salida con pruebas/referencia/.
prueba: paraest.x
	./pruebas/regresion.sh ./paraest.x

clean:
	rm -rf obj debug paraest.x debug.paraest.x
