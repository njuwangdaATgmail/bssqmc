program test_env
! a test of the environment: BLAS + LAPACK + MPI
#ifdef MPI
  use mpi
#endif
  use matrixlib
  implicit none
  integer,parameter :: n=2000
  integer :: id=0
  integer ierr
  real(8) a(n,n),b(n,n),c(n,n)
  real(8) t0,t1

#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,id,ierr)
#endif

  call random_number(a)
  call random_number(b)

!------------------------------------------
#ifdef MPI
  t0=mpi_wtime()
#else
  call cpu_time(t0)
#endif

  c=matmul(a,b)
#ifdef MPI
  t1=mpi_wtime()
#else
  call cpu_time(t1)
#endif  
  print*,'matmul (with OR without --external-blas):',t1-t0,'s,  on core',id
  
!-----------------------------------------  
  t0=t1
  call dgemm('n','n',n,n,n,1d0,a,n,b,n,0d0,c,n)
#ifdef MPI
  t1=mpi_wtime()
#else
  call cpu_time(t1)
#endif  
  print*,'dgemm:',t1-t0,'s,  on core',id

!---------------------------------------
  t0=t1
  call inverse(n,c)
#ifdef MPI
  t1=mpi_wtime()
#else
  call cpu_time(t1)
#endif
  print*,'inverse:',t1-t0,'s,  on core',id

!----------------------------------------
#ifdef MPI
  call mpi_finalize(ierr)
#endif

end
