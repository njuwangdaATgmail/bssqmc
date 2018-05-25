#===========================================================
# system parameters
F90 = mpif90 #ifort
FCOPTS = -O2 -cpp -DMPI -fexternal-blas #-g -Wall #-Wextra -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow #-CB -CU -g -trace -W1 
LDOPTS = -O2
LIBS = -lopenblas #-mkl=sequential -limf #-static
#==========================================================
# define all obj files
OBJS = mod_matrixlib.o test_env.o 
#==========================================================
# link all to generate exe file
test_env.x: $(OBJS)
	$(F90) $(LDOPTS) $(OBJS) $(LIBS) -o test_env.x
#==========================================================
clean:
	rm *.o *.mod
#==========================================================
# generate every obj and module files
mod_matrixlib.o: ../src/mod_matrixlib.f90
	$(F90) $(FCOPTS) -c ../src/mod_matrixlib.f90

test_env.o: test_env.f90
	$(F90) $(FCOPTS) -c test_env.f90
