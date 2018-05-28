#===========================================================
# system parameters
F90 = mpif90
FCOPTS = -O2 -cpp -DMPI -fexternal-blas -g -fbacktrace -fbounds-check -Wall -Wextra
LDOPTS = -O2
LIBS = -lopenblas #-mkl=sequential -limf 
#==========================================================
# define all obj files
OBJS = mod_randomlib.o mod_matrixlib.o mod_measpool.o mod_dqmc_complex.o mod_model_complex.o \
	init.o measurement.o postprocess.o mod_hubhol.o tmpout.o
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

mod_dqmc_complex.o: ../src/mod_dqmc_complex.f90 mod_matrixlib.o mod_randomlib.o mod_measpool.o
	$(F90) $(FCOPTS) -c ../src/mod_dqmc_complex.f90

mod_model_complex.o: ../src/mod_model_complex.f90 mod_dqmc_complex.o
	$(F90) $(FCOPTS) -c ../src/mod_model_complex.f90

mod_hubhol.o: mod_hubhol.f90 mod_model_complex.o mod_dqmc_complex.o
	$(F90) $(FCOPTS) -c mod_hubhol.f90

init.o: init.f90 mod_model_complex.o mod_dqmc_complex.o
	$(F90) $(FCOPTS) -c init.f90

measurement.o: measurement.f90 mod_model_complex.o mod_dqmc_complex.o
	$(F90) $(FCOPTS) -c measurement.f90

postprocess.o: postprocess.f90 mod_model_complex.o mod_dqmc_complex.o
	$(F90) $(FCOPTS) -c postprocess.f90

tmpout.o: tmpout.f90 mod_dqmc_complex.o mod_model_complex.o
	$(F90) $(FCOPTS) -c tmpout.f90
