program test_env
! a test of the environment: BLAS + LAPACK + MPI
#ifdef MPI
  use mpi
#endif
  use matrixlib
  implicit none
  integer,parameter :: n=2000
  integer,parameter :: s=2
  integer :: i,j,info
  integer :: id=0
  integer ierr
  real(8) a(n,n),b(n,n),c(n,n)
  real(8) t0,t1
  real(8) d(s,s),work(64*s),v(s),r(s,s),m(s,s),g(s,s)
  open(8,file='t.txt')
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
 do i=1,s
   do j=1,s
      d(i,j)=i+j
      end do
  end do    
 do i=1,s
   do j=1,s
      g(i,j)=d(i,j)
    end do
  end do  
  call dsyev('V','U',s,d,s,v,work,64*s,info)
  do i=1,s
   do j=1,s
     m(i,j)=d(i,j)
     end do
  end do
  write(*,*)'v:'
  write(*,*)v

  call inverse(s,m)
  r=matmul(d,g)
  r=matmul(r,m)
  write(*,*)'r:'
  write(*,*) r
  do i=1,s
    write(8,*) v(i)
  
  end do
 close(8)
 open(8,file='t.txt')
 do i=1,s
  do j=1,s
   if(i==j) then
   read(8,*)d(i,j)
   else 
           d(i,j)=0
   end if
  end do
 
 end do 
 do i=1,s
  do j=1,s
   if (abs(d(i,j)-r(i,j))<0.00001) then
   write(*,*)'T'
   else 
   write(*,*)'f'
   end if
  end do
end do 
 close(8)
end
