#===========================================================
# system parameters
F90 = mpif90
FCOPTS = -O2 -cpp -DMPI -fexternal-blas -fbounds-check -Wall -Wextra #-g -fbacktrace 
LDOPTS = -O2
LIBS = -lopenblas
#==========================================================
# define all obj files
OBJS = mod_randomlib.o mod_matrixlib.o mod_measpool.o mod_dqmc_core.o dqmc.o \
	init.o measurement.o postprocess.o hubhol.o tmpout.o
#==========================================================
# link all to generate exe file
dqmc.x: $(OBJS)
	$(F90) $(LDOPTS) $(OBJS) $(LIBS) -o dqmc.x
#==========================================================
clean:
	rm *.o *.mod
#==========================================================
# generate every obj and module files
mod_matrixlib.o: ../src/mod_matrixlib.f90
	$(F90) $(FCOPTS) -c ../src/mod_matrixlib.f90

mod_randomlib.o: ../src/mod_randomlib.f90
	$(F90) $(FCOPTS) -c ../src/mod_randomlib.f90

mod_measpool.o: ../src/mod_measpool.f90
	$(F90) $(FCOPTS) -c ../src/mod_measpool.f90

mod_dqmc_core.o: ../src/mod_dqmc_core.f90 mod_matrixlib.o mod_randomlib.o mod_measpool.o
	$(F90) $(FCOPTS) -c ../src/mod_dqmc_core.f90

dqmc.o: ../src/dqmc.f90 mod_dqmc_core.o
	$(F90) $(FCOPTS) -c ../src/dqmc.f90

hubhol.o: hubhol.f90 dqmc.o
	$(F90) $(FCOPTS) -c hubhol.f90

init.o: init.f90 dqmc.o hubhol.o
	$(F90) $(FCOPTS) -c init.f90

measurement.o: measurement.f90 dqmc.o hubhol.o
	$(F90) $(FCOPTS) -c measurement.f90

postprocess.o: postprocess.f90 dqmc.o hubhol.o
	$(F90) $(FCOPTS) -c postprocess.f90

tmpout.o: tmpout.f90 dqmc.o hubhol.o
	$(F90) $(FCOPTS) -c tmpout.f90
