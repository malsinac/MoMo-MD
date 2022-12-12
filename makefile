MAKEFLAGS+="-j 4"

FC=gfortran

EFLAGS=-pedantic -Wall -Wextra -Wdo-subscript -Waliasing -Winteger-division -Wsurprising -Wuse-without-only # Error flags
FFLAGS=-ffree-line-length-none -ffree-form -std=gnu # Fortran-language flags
DFLAGS=-g -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wrealloc-lhs # Debugging flags
PFLAGS=-march=native -O3 # Performance flags


all: main.x

# ~~ LINKING ~~
main.x: main.o analysis_m.o init_conditions_mod.o integrators_mod.o mathutils_m.o thermodynamics_mod.o thermostats_m.o writers_mod.o interfaces_m.o
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) $^ -o $@

# ~~ COMPILING ~~
main.o: main.F90 
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

analysis_m.o: analysis_m.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

init_conditions_mod.o: init_conditions_mod.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

integrators_mod.o: integrators_mod.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

mathutils_m.o: mathutils_m.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

thermodynamics_mod.o: thermodynamics_mod.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

thermostats_m.o: thermostats_m.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

writers_mod.o: writers_mod.F90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

interfaces_m.o: interfaces_m.f90
	$(FC) $(EFLAGS) $(FFLAGS) $(DFLAGS) $(PFLAGS) -c $^

.PHONY: clean
clean:
	rm *.o