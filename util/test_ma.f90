program test_env
! a test of the environment: BLAS + LAPACK + MPI
#ifdef MPI
  use mpi
#endif
  use matrixlib
  implicit none
  integer,parameter ::n=4
  integer,parameter ::m=2
  integer,parameter ::s=4
  integer :: i,j,info
  real(8) r(n,n)
  real(8) a(n,n)
  real(8) c(n,n)
  real(8) d(n,n)
  real(8) work(64*s),v(s)
  call creat_m(n,n,a)
  call print_m(n,n,a)
  data c /4,6,5,7,8,9,1,0,2,3,4,5,7,8,9,1/
  data d /4,6,5,7,8,9,1,0,2,3,4,5,7,8,9,1/

  write(*,*)"trace:"
  write(*,*)trace(n,a)
  write(*,*)"eigen"
  call eigen(n,a,v)
  do i=1,n
    write(*,*)v(i)
   end do 
  write(*,*)"v:"
  write(*,*) a(n,n)
  call print_m(n,n,a)
  call creat_m(n,n,a)
  call print_m(n,n,a)
  call inverse(n,c)
  call print_m(n,n,c)
  write(*,*)"det:"
  write(*,*)det(n,a)
  call print_m(n,n,d)
  call qr(n,n,d,r)
  call print_m(n,n,d)
  call print_m(n,n,r)
  call print_m(n,n,matmul(d,r))
  
 














   end
   subroutine creat_m(n,m,a)
          implicit none
          integer :: n
          integer :: m
          real(8) a(n,m)
          integer :: i,j
          do i=1,n
             do j=1,m
               a(i,j)=i+j
               end do
           end do
           end
  subroutine print_m(n,m,a)
          implicit none
          integer :: n
          integer :: m
          real(8) a(n,m)
          integer :: i,j
          write(*,*) "matrix:"
          write(*,"(/)")
          do i=1,n
           do j=1,m
           if(j/=m)then
               write(*,*) a(i,j)
           else
                   write(*,*) a(i,j)
                   write(*,"(/)")
           end if
            end do
           end do
           
end
