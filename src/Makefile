#===========================================================
# system parameters
F90 = mpif90
FCOPTS = -O2 -cpp -DMPI -fexternal-blas -fbounds-check -Wall -Wextra #-g -fbacktrace 
LDOPTS = -O2
LIBS = -lopenblas
#==========================================================
# define all obj files
OBJS = mod_randomlib.o mod_matrixlib.o mod_measpool.o mod_dqmc_core.o dqmc.o \
	dqmc_init.o dqmc_measurement.o dqmc_postprocess.o dqmc_tmpout.o
#==========================================================
# link all to generate exe file
all: exe lib
exe: $(OBJS) dqmc_external.o
	$(F90) $(LDOPTS) $(OBJS) dqmc_external.o $(LIBS) -o dqmc.x
lib: $(OBJS)
	ar cr libdqmc.a $(OBJS)
#==========================================================
clean:
	rm *.o *.mod
#==========================================================
# generate every obj and module files
mod_matrixlib.o: mod_matrixlib.f90
	$(F90) $(FCOPTS) -c mod_matrixlib.f90

mod_randomlib.o: mod_randomlib.f90
	$(F90) $(FCOPTS) -c mod_randomlib.f90

mod_measpool.o: mod_measpool.f90
	$(F90) $(FCOPTS) -c mod_measpool.f90

mod_dqmc_core.o: mod_dqmc_core.f90 mod_matrixlib.o mod_randomlib.o mod_measpool.o
	$(F90) $(FCOPTS) -c mod_dqmc_core.f90

dqmc.o: dqmc.f90 mod_dqmc_core.o
	$(F90) $(FCOPTS) -c dqmc.f90

dqmc_init.o: dqmc_init.f90 dqmc.o
	$(F90) $(FCOPTS) -c dqmc_init.f90

dqmc_measurement.o: dqmc_measurement.f90 dqmc.o
	$(F90) $(FCOPTS) -c dqmc_measurement.f90

dqmc_tmpout.o: dqmc_tmpout.f90 dqmc.o
	$(F90) $(FCOPTS) -c dqmc_tmpout.f90

dqmc_postprocess.o: dqmc_postprocess.f90 dqmc.o
	$(F90) $(FCOPTS) -c dqmc_postprocess.f90

dqmc_external.o: dqmc_external.f90 dqmc.o
	$(F90) $(FCOPTS) -c dqmc_external.f90
