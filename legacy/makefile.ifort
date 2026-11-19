# Makefile para GPEC en fortran, compilado con Intel Fortran 11
#
# Definición de variables
objetos = obj/a-paraest-$(vapara).o \
	  obj/A-DCADJUST-1.1.0.o \
	  obj/A-dcadjust-fast-$(vADCAdfast).o \
	  obj/A-gamma-$(vGAMMA).o \
	  obj/A-gc-eos-3.1.1.o \
	  obj/A-isof-mezclas-$(vIMez).o \
	  obj/A-isof-puros-1.1.0.o \
	  obj/A-psat-$(vPSAT).o \
	  obj/A-rutin-mezclas-$(vMxs).o \
	  obj/A-Rutin-puros-$(vRPur).o \
	  obj/ASSOC-GC-EOSFev.o \
	  obj/ZxssqFev-$(vZXSSQ).o \
	  obj/GC3-$(vGC3)-paraest.o \
	  obj/GCA-v$(vGCA).o \
	  obj/GCDBA-v$(vGCADB).paraest.o \
	  obj/GCXS-v$(vGCXS).paraest.o \
	  obj/GCSLE-v$(vGCSLE).paraest.o \
	  obj/semiImpCrit-v$(vPC).o
#
# Objetos con información para debuggear
objDebug = debug/a-paraest-$(vapara).o \
	   debug/A-DCADJUST-1.1.0.o \
	   debug/A-dcadjust-fast-$(vADCAdfast).o \
	   debug/A-gamma-$(vGAMMA).o \
	   debug/A-gc-eos-3.1.1.o \
	   debug/A-isof-mezclas-$(vIMez).o \
	   debug/A-isof-puros-1.1.0.o \
	   debug/A-psat-$(vPSAT).o \
	   debug/A-rutin-mezclas-$(vMxs).o \
	   debug/A-Rutin-puros-$(vRPur).o \
	   debug/ASSOC-GC-EOSFev.o \
	   debug/ZxssqFev-$(vZXSSQ).o \
	   debug/GC3-$(vGC3)-paraest.o \
	   debug/GCA-v$(vGCA).o \
	   debug/GCDBA-v$(vGCADB).paraest.o \
	   debug/GCXS-v$(vGCXS).paraest.o \
	   debug/GCSLE-v$(vGCSLE).paraest.o \
	   debug/semiImpCrit-v$(vPC).o	   
	  
f90comp = ifort
#
# Debug con traza
switchDebug = -g -traceback -diag-disable remark
#
# Riguroso
switchDebug = -g -traceback -debug extended -fpe0 -check all -warn
switchDebug = -g -traceback -debug extended -fpe0 -diag-disable remark
#switchDebug = -g -traceback
#
vADCAdfast = 1.3.0
#vapara = 1.8.95.0
vapara = 1.9.04.1.densidad
vapara = 1.12.0
vGAMMA = 1.1.0
#vGCA = 1.9.26.6.H2O
vGCA = 1.9.28.1
#1.9.25.1
vGC3 = 1.2.101.0
vGCADB = 1.3.3
vGCXS = 0.5.0.2
vPSAT = 1.0.2
vIMez = 1.1.0
vRPur = 1.04
vZXSSQ = 1.0.3
vZXSSQ = 2.0.0
vGCSLE = 0.2.1
vPC = 0.0.1
vMxs = 1.1.0


# Optimización
#switch = -O3 -parallel -static-intel -guide-par
switch = -O3 -mp1

# Makefile
paraest: $(objetos)
	$(f90comp) -o paraest.x $(switch) $(objetos)

obj/GCDBA-v$(vGCADB).paraest.o: src/GCDBA-v$(vGCADB).paraest.for
	$(f90comp) -c $(switch) src/GCDBA-v$(vGCADB).paraest.for -o obj/GCDBA-v$(vGCADB).paraest.o
	
obj/GC3-$(vGC3)-paraest.o: src/GC3-$(vGC3)-paraest.for
	$(f90comp) -c $(switch) src/GC3-$(vGC3)-paraest.for -o obj/GC3-$(vGC3)-paraest.o
	
obj/GCA-v$(vGCA).o: src/GCA-v$(vGCA).f90
	$(f90comp) -c $(switch) src/GCA-v$(vGCA).f90 -o obj/GCA-v$(vGCA).o
	
obj/a-paraest-$(vapara).o: src/a-paraest-$(vapara).f90
	$(f90comp) -c $(switch) src/a-paraest-$(vapara).f90 -o obj/a-paraest-$(vapara).o
	
obj/A-DCADJUST-1.1.0.o: src/A-DCADJUST-1.1.0.FOR
	$(f90comp) -c $(switch) src/A-DCADJUST-1.1.0.FOR -o obj/A-DCADJUST-1.1.0.o
	
obj/A-dcadjust-fast-$(vADCAdfast).o: src/A-dcadjust-fast-$(vADCAdfast).for
	$(f90comp) -c $(switch) src/A-dcadjust-fast-$(vADCAdfast).for -o obj/A-dcadjust-fast-$(vADCAdfast).o
	
obj/A-gamma-$(vGAMMA).o: src/A-gamma-$(vGAMMA).for
	$(f90comp) -c $(switch) src/A-gamma-$(vGAMMA).for -o obj/A-gamma-$(vGAMMA).o
	
obj/A-gc-eos-3.1.1.o: src/A-gc-eos-3.1.1.for
	$(f90comp) -c $(switch) src/A-gc-eos-3.1.1.for -o obj/A-gc-eos-3.1.1.o
	
obj/A-isof-mezclas-$(vIMez).o: src/A-isof-mezclas-$(vIMez).for
	$(f90comp) -c $(switch) src/A-isof-mezclas-$(vIMez).for -o obj/A-isof-mezclas-$(vIMez).o
	
obj/A-isof-puros-1.1.0.o: src/A-isof-puros-1.1.0.for
	$(f90comp) -c $(switch) src/A-isof-puros-1.1.0.for -o obj/A-isof-puros-1.1.0.o
	
obj/A-psat-$(vPSAT).o: src/A-psat-$(vPSAT).for
	$(f90comp) -c $(switch) src/A-psat-$(vPSAT).for -o obj/A-psat-$(vPSAT).o
	
obj/A-rutin-mezclas-$(vMxs).o: src/A-rutin-mezclas-$(vMxs).for
	$(f90comp) -c $(switch) src/A-rutin-mezclas-$(vMxs).for -o obj/A-rutin-mezclas-$(vMxs).o
	
obj/A-Rutin-puros-$(vRPur).o: src/A-Rutin-puros-$(vRPur).for
	$(f90comp) -c $(switch) src/A-Rutin-puros-$(vRPur).for -o obj/A-Rutin-puros-$(vRPur).o
	
#obj/A-TP-FLASH.o: src/A-TP-FLASH.f90
#	$(f90comp) -c $(switch) src/A-TP-FLASH.f90 -o obj/A-TP-FLASH.o

#obj/A-TX-FLASH.o: src/A-TX-FLASH.FOR
#	$(f90comp) -c $(switch) src/A-TX-FLASH.FOR -o obj/A-TX-FLASH.o

obj/ASSOC-GC-EOSFev.o: src/ASSOC-GC-EOSFev.FOR
	$(f90comp) -c $(switch) src/ASSOC-GC-EOSFev.FOR -o obj/ASSOC-GC-EOSFev.o
	
#obj/sistemaeca.o: src/sistemaeca.for
#	$(f90comp) -c $(switch) src/sistemaeca.for -o obj/sistemaeca.o
# obj/ZxssqFev-$(vZXSSQ).o: src/ZxssqFev-$(vZXSSQ).for
# 	$(f90comp) -c $(switch) src/ZxssqFev-$(vZXSSQ).for -o obj/ZxssqFev-$(vZXSSQ).o
obj/ZxssqFev-$(vZXSSQ).o: src/ZxssqFev-$(vZXSSQ).f90
	$(f90comp) -c $(switch) src/ZxssqFev-$(vZXSSQ).f90 -o obj/ZxssqFev-$(vZXSSQ).o	
	
obj/GCXS-v$(vGCXS).paraest.o: src/GCXS-v$(vGCXS).paraest.f90
	$(f90comp) -c $(switch) src/GCXS-v$(vGCXS).paraest.f90 -o obj/GCXS-v$(vGCXS).paraest.o

obj/GCSLE-v$(vGCSLE).paraest.o: src/GCSLE-v$(vGCSLE).paraest.f90
	$(f90comp) -c $(switch) src/GCSLE-v$(vGCSLE).paraest.f90 -o obj/GCSLE-v$(vGCSLE).paraest.o
	
obj/semiImpCrit-v$(vPC).o: src/semiImpCrit-v$(vPC).f90
	$(f90comp) -c $(switch) src/semiImpCrit-v$(vPC).f90 -o obj/semiImpCrit-v$(vPC).o
#
#--------------------------------------------------------------------------------------------------------------
#
debug: $(objDebug)
	$(f90comp) -o debug.paraest.x $(switchDebug) $(objDebug)

debug/GCDBA-v$(vGCADB).paraest.o: src/GCDBA-v$(vGCADB).paraest.for
	$(f90comp) -c $(switchDebug) src/GCDBA-v$(vGCADB).paraest.for -o debug/GCDBA-v$(vGCADB).paraest.o
	
debug/GC3-$(vGC3)-paraest.o: src/GC3-$(vGC3)-paraest.for
	$(f90comp) -c $(switchDebug) src/GC3-$(vGC3)-paraest.for -o debug/GC3-$(vGC3)-paraest.o
	
debug/GCA-v$(vGCA).o: src/GCA-v$(vGCA).f90
	$(f90comp) -c $(switchDebug) src/GCA-v$(vGCA).f90 -o debug/GCA-v$(vGCA).o
	
debug/a-paraest-$(vapara).o: src/a-paraest-$(vapara).f90
	$(f90comp) -c $(switchDebug) src/a-paraest-$(vapara).f90 -o debug/a-paraest-$(vapara).o
	
debug/A-DCADJUST-1.1.0.o: src/A-DCADJUST-1.1.0.FOR
	$(f90comp) -c $(switchDebug) src/A-DCADJUST-1.1.0.FOR -o debug/A-DCADJUST-1.1.0.o
	
debug/A-dcadjust-fast-$(vADCAdfast).o: src/A-dcadjust-fast-$(vADCAdfast).for
	$(f90comp) -c $(switchDebug) src/A-dcadjust-fast-$(vADCAdfast).for -o debug/A-dcadjust-fast-$(vADCAdfast).o
	
debug/A-gamma-$(vGAMMA).o: src/A-gamma-$(vGAMMA).for
	$(f90comp) -c $(switchDebug) src/A-gamma-$(vGAMMA).for -o debug/A-gamma-$(vGAMMA).o
	
debug/A-gc-eos-3.1.1.o: src/A-gc-eos-3.1.1.for
	$(f90comp) -c $(switchDebug) src/A-gc-eos-3.1.1.for -o debug/A-gc-eos-3.1.1.o
	
debug/A-isof-mezclas-$(vIMez).o: src/A-isof-mezclas-$(vIMez).for
	$(f90comp) -c $(switchDebug) src/A-isof-mezclas-$(vIMez).for -o debug/A-isof-mezclas-$(vIMez).o
	
debug/A-isof-puros-1.1.0.o: src/A-isof-puros-1.1.0.for
	$(f90comp) -c $(switchDebug) src/A-isof-puros-1.1.0.for -o debug/A-isof-puros-1.1.0.o
	
debug/A-psat-$(vPSAT).o: src/A-psat-$(vPSAT).for
	$(f90comp) -c $(switchDebug) src/A-psat-$(vPSAT).for -o debug/A-psat-$(vPSAT).o
	
debug/A-rutin-mezclas-$(vMxs).o: src/A-rutin-mezclas-$(vMxs).for
	$(f90comp) -c $(switchDebug) src/A-rutin-mezclas-$(vMxs).for -o debug/A-rutin-mezclas-$(vMxs).o
	
debug/A-Rutin-puros-$(vRPur).o: src/A-Rutin-puros-$(vRPur).for
	$(f90comp) -c $(switchDebug) src/A-Rutin-puros-$(vRPur).for -o debug/A-Rutin-puros-$(vRPur).o
	
#debug/A-TP-FLASH.o: src/A-TP-FLASH.f90
#	$(f90comp) -c $(switchDebug) src/A-TP-FLASH.f90 -o debug/A-TP-FLASH.o

#debug/A-TX-FLASH.o: src/A-TX-FLASH.FOR
#	$(f90comp) -c $(switchDebug) src/A-TX-FLASH.FOR -o debug/A-TX-FLASH.o

debug/ASSOC-GC-EOSFev.o: src/ASSOC-GC-EOSFev.FOR
	$(f90comp) -c $(switchDebug) src/ASSOC-GC-EOSFev.FOR -o debug/ASSOC-GC-EOSFev.o
	
#debug/sistemaeca.o: src/sistemaeca.for
#	$(f90comp) -c $(switchDebug) src/sistemaeca.for -o debug/sistemaeca.o
# debug/ZxssqFev-$(vZXSSQ).o: src/ZxssqFev-$(vZXSSQ).for
# 	$(f90comp) -c $(switchDebug) src/ZxssqFev-$(vZXSSQ).for -o debug/ZxssqFev-$(vZXSSQ).o
debug/ZxssqFev-$(vZXSSQ).o: src/ZxssqFev-$(vZXSSQ).f90
	$(f90comp) -c $(switchDebug) src/ZxssqFev-$(vZXSSQ).f90 -o debug/ZxssqFev-$(vZXSSQ).o
	
debug/GCXS-v$(vGCXS).paraest.o: src/GCXS-v$(vGCXS).paraest.f90
	$(f90comp) -c $(switchDebug) src/GCXS-v$(vGCXS).paraest.f90 -o debug/GCXS-v$(vGCXS).paraest.o

debug/GCSLE-v$(vGCSLE).paraest.o: src/GCSLE-v$(vGCSLE).paraest.f90
	$(f90comp) -c $(switchDebug) src/GCSLE-v$(vGCSLE).paraest.f90 -o debug/GCSLE-v$(vGCSLE).paraest.o

debug/semiImpCrit-v$(vPC).o: src/semiImpCrit-v$(vPC).f90
	$(f90comp) -c $(switchDebug) src/semiImpCrit-v$(vPC).f90 -o debug/semiImpCrit-v$(vPC).o
	
# Limpiando "objetos"
clean:
	rm obj/*.o
cdebug:
	rm debug/*.o
# Limpiar backups
cleanbak:
	rm src/*.bak *.bak
# Limpiar todo:
veryclean:
	rm obj/*.o debug/*.o *.bak src/*.bak *~ *genmod*
# Instalando localmente
install:
	cp *.x $(HOME)/.local/bin
# Fin del makefile
