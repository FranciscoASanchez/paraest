# =============================================================================
#  paraest makefile
#
#  Usage:
#     make                    optimized build with ifx       -> paraest.x
#     make FC=gfortran        optimized build with gfortran
#     make debug              build with run-time checks     -> debug.paraest.x
#     make test               build and run the regression tests
#     make clean              remove objects and executables
#
#  Objects go to obj/<compiler>/ (or debug/<compiler>/), so ifx and
#  gfortran objects never mix.
# =============================================================================

# make predefines FC=f77; if nobody chose a compiler, use ifx.
ifeq ($(origin FC),default)
  FC = ifx
endif

# --- Source files, by layer -------------------------------------------------
MAIN         = src/main/dimensions.f90                   \
               src/main/config.f90                       \
               src/main/fit_data.f90                     \
               src/main/statistics.f90                   \
               src/main/scaling.f90                      \
               src/main/utilities.f90                    \
               src/main/model_output.f90                 \
               src/main/residuals.f90                    \
               src/main/paraest.f90

OPTIMIZATION = src/optimization/zxssq.f90                 \
               src/optimization/pso.f90

MODELS       = src/models/thermo_model.f90               \
               src/models/gca_parameters.f90             \
               src/models/gca.f90                        \
               src/models/rkpr.f90                       \
               src/models/thermopack.f90                 \
               src/models/spline_temperature.f90         \
               src/models/gca_bridge.f90                 \
               src/models/fugacity.f90                   \
               src/models/fit_dc.f90

EQUILIBRIUM  = src/equilibrium/gc3_flash.f90             \
               src/equilibrium/isofugacity.f90           \
               src/equilibrium/gcdba.f90                 \
               src/equilibrium/tpf_flash.f90             \
               src/equilibrium/gcxs.f90                  \
               src/equilibrium/gcsle.f90                 \
               src/equilibrium/pcrit_semi_implicit.f90

# Build order: modules first, the main program last. Within each layer, the
# routines called by others go first (see the note on *__genmod* below).
SOURCES = $(OPTIMIZATION) $(MODELS) $(EQUILIBRIUM) $(MAIN)

# --- Compiler options -------------------------------------------------------
ifeq ($(FC),gfortran)
  OPT    = -O2
  DEBUG  = -g -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall
  COMMON = -std=legacy -ffree-line-length-none
  MODS   = -J
else
  OPT    = -O3
  DEBUG  = -g -O0 -traceback -debug extended -fpe0 -check all,nouninit,noarg_temp_created -warn
  COMMON = -diag-disable=remark
  MODS   = -module
endif

# --- Derived (nothing to change below this line) ----------------------------
COMP     = $(notdir $(FC))
OBJDIR   = obj/$(COMP)
DBGDIR   = debug/$(COMP)
objs     = $(addprefix $(1)/,$(addsuffix .o,$(basename $(notdir $(SOURCES)))))
OBJECTS  = $(call objs,$(OBJDIR))
OBJDEBUG = $(call objs,$(DBGDIR))

vpath %.f90 src/main src/optimization src/models src/equilibrium

.PHONY: all debug test clean
all: paraest.x

# The stamp obj/.compiler-<comp> remembers which compiler built the
# executables: when the compiler changes, the old stamp is removed and a new
# one, newer than paraest.x, is created, so the executable is relinked.
STAMP    = obj/.compiler-$(COMP)

paraest.x: $(OBJECTS) $(STAMP)
	$(FC) $(OPT) -o $@ $(OBJECTS)

debug: debug.paraest.x
debug.paraest.x: $(OBJDEBUG) $(STAMP)
	$(FC) $(DEBUG) -o $@ $(OBJDEBUG)

$(STAMP):
	mkdir -p obj
	rm -f obj/.compiler-* obj/.compilador-*
	touch $@

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(COMMON) $(OPT) $(MODS) $(OBJDIR) -c $< -o $@
$(DBGDIR)/%.o: %.f90 | $(DBGDIR)
	$(FC) $(COMMON) $(DEBUG) $(MODS) $(DBGDIR) -c $< -o $@

$(OBJDIR) $(DBGDIR):
	mkdir -p $@

# Module dependencies: when a module changes, whatever uses it is rebuilt
# (and make builds the module first).
MODULES = dimensions config fit_data gca_parameters statistics scaling utilities residuals pso model_output gca_bridge
dep = $(addprefix $(1)/,$(addsuffix .o,$(2)))
$(OBJDIR)/paraest.o: $(call dep,$(OBJDIR),$(MODULES))
$(DBGDIR)/paraest.o: $(call dep,$(DBGDIR),$(MODULES))
$(OBJDIR)/fit_data.o $(OBJDIR)/gca_parameters.o $(OBJDIR)/scaling.o: $(OBJDIR)/dimensions.o
$(OBJDIR)/scaling.o: $(OBJDIR)/config.o
$(DBGDIR)/fit_data.o $(DBGDIR)/gca_parameters.o $(DBGDIR)/scaling.o: $(DBGDIR)/dimensions.o
$(DBGDIR)/scaling.o: $(DBGDIR)/config.o
$(OBJDIR)/gca_bridge.o $(OBJDIR)/model_output.o: $(OBJDIR)/gca_parameters.o
$(DBGDIR)/gca_bridge.o $(DBGDIR)/model_output.o: $(DBGDIR)/gca_parameters.o
# gca.f90 holds module pargr_export
$(OBJDIR)/gca_bridge.o: $(OBJDIR)/gca.o
$(DBGDIR)/gca_bridge.o: $(DBGDIR)/gca.o
# gca.f90 also holds module zmax_status
$(OBJDIR)/tpf_flash.o: $(OBJDIR)/gca.o
$(DBGDIR)/tpf_flash.o: $(DBGDIR)/gca.o
$(OBJDIR)/fugacity.o: $(OBJDIR)/thermo_model.o
$(DBGDIR)/fugacity.o: $(DBGDIR)/thermo_model.o
RESDEPS = dimensions config fit_data gca_parameters gca_bridge scaling utilities tpf_flash
$(OBJDIR)/residuals.o: $(call dep,$(OBJDIR),$(RESDEPS))
$(DBGDIR)/residuals.o: $(call dep,$(DBGDIR),$(RESDEPS))

# Regression tests: every case in tests/cases/ runs in tests/tmp/ and its
# output is compared with the case's reference/ directory.
test: paraest.x
	./tests/regression.sh ./paraest.x

clean:
	rm -rf obj debug paraest.x debug.paraest.x

install:
	cp *.x ~/.local/bin
