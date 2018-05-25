!> A pool to store measured observables during the Monte Carlo simulation.
!! It is designed to help Monte Carlo users easily measure observables and calculate error bars.
!! The user only need to perform two kinds of operations: put and get, i.e. put in different
!! measured observables and finnaly read out their mean values and errorbars (or different moments) directly.
!! Some other functions (including initialization, save and load, etc.) can be performed
!! by the MODULE dqmc_complex and The users don't need to care about them.
!!
!! For users who don't care about the qmc details but only need measure some observables:
!!    put_pool(a)
!!    put_pool(n,a)
!!    get_pool(a)
!!    get_pool(a,a_err)
!!    get_pool(n,a)
!!    get_pool(n,a,a_err)
!! The above can be used for scaler and vector variables.
!! The following can be called to measure array (rank>=2):
!! (If we use interface put_pool/get_pool directly, the compiler may report error.)
!!    xput_array(n,a)
!!    xget_array(n,a)
!!    xget_array_err(n,a,aerr)
!! where x=r,z.
!!
!! **IMPORTANT**: the first element in the pool must be the sign of the current configuration.
!!
!! For QMC developers:
!!    set_pool(nbin_,size_r_,size_z_) size_r_ and size_z_ should be given by users
!!    save_pool(): save the pool temporarily
!!    load_pool(): load previous pool
!!    move_pool(): move to the next bin
!!    add_pool(): measurement times +1, reset head_r=head_z=1
!!    average_pool(): do average on each core and then among different cores
!!
MODULE measpool

#ifdef MPI
  USE mpi
#endif

  INTEGER, PRIVATE :: nbin=1
  INTEGER, PRIVATE :: bin=1
  INTEGER, PRIVATE :: size_r, size_z      ! size of the pool
  INTEGER, PRIVATE :: head_r, head_z      ! current position of the pool
  INTEGER, ALLOCATABLE, PRIVATE :: ndata(:)       ! (nbin)
  REAL(8), ALLOCATABLE, PRIVATE :: work_r(:,:)    ! (size_r,nbin)
  COMPLEX(8), ALLOCATABLE, PRIVATE :: work_z(:,:) ! (size_z,nbin)
  REAL(8), ALLOCATABLE, PRIVATE :: work2_r(:)    ! (size_r)
  COMPLEX(8), ALLOCATABLE, PRIVATE :: work2_z(:) ! (size_z)

  INTERFACE add_pool
    MODULE PROCEDURE add_pool_r,add_pool_z
  END INTERFACE

  INTERFACE put_pool
    MODULE PROCEDURE rput_number,rput_array
    MODULE PROCEDURE zput_number,zput_array
  END INTERFACE

  INTERFACE get_pool
    MODULE PROCEDURE rget_number,rget_number_err,rget_array,rget_array_err
    MODULE PROCEDURE zget_number,zget_number_err,zget_array,zget_array_err
  END INTERFACE

#ifdef MPI
  INTERFACE dmpi_mean
    MODULE PROCEDURE dmpi_mean,dmpi_mean_scalar
  END INTERFACE

  INTERFACE zmpi_mean
    MODULE PROCEDURE zmpi_mean,zmpi_mean_scalar
  END INTERFACE

  INTERFACE mpi_mean
    MODULE PROCEDURE dmpi_mean,dmpi_mean_scalar,zmpi_mean,zmpi_mean_scalar
  END INTERFACE
#endif

CONTAINS

  SUBROUTINE set_pool(nbin_,size_r_,size_z_)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbin_,size_r_,size_z_
    nbin=nbin_
    size_r=size_r_
    size_z=size_z_
    bin=1
    ALLOCATE(ndata(nbin))
    ndata(:)=0
    IF(size_r>0)THEN
      ALLOCATE(work_r(size_r,nbin))
      work_r(:,:)=0d0
      head_r=1
    END IF
    IF(size_z>0)THEN
      ALLOCATE(work_z(size_z,nbin))
      work_z(:,:)=(0d0,0d0)
      head_z=1
    END IF
  END SUBROUTINE

  SUBROUTINE load_pool()
    IMPLICIT NONE
    INTEGER id,ierr
    CHARACTER*4 cid
#ifdef MPI
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
#else
    id=0
#endif
    WRITE(cid,'(1i4)') id
    OPEN(81,FILE='pool.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    READ(81) nbin,bin,size_r,size_z,head_r,head_z
    IF(.not.allocated(ndata))ALLOCATE(ndata(nbin))
    READ(81) ndata(:)
    IF(size_r>0)THEN
      IF(.not.allocated(work_r))ALLOCATE(work_r(size_r,nbin))
      READ(81) work_r(:,:)
    END IF
    IF(size_z>0)THEN
      IF(.not.allocated(work_z))ALLOCATE(work_z(size_z,nbin))
      READ(81) work_z(:,:)
    END IF
    CLOSE(81)
  END SUBROUTINE

  SUBROUTINE save_pool()
    IMPLICIT NONE
    INTEGER id,ierr
    CHARACTER*4 cid
#ifdef MPI
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
#else
    id=0
#endif
    WRITE(cid,'(1i4)') id
    OPEN(81,FILE='pool.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    WRITE(81) nbin,bin,size_r,size_z,head_r,head_z
    WRITE(81) ndata(:)
    IF(size_r>0) WRITE(81) work_r(:,:)
    IF(size_z>0) WRITE(81) work_z(:,:)
    CLOSE(81)
  END SUBROUTINE

  SUBROUTINE move_pool()
    IMPLICIT NONE
    bin=bin+1
    head_r=1
    head_z=1
  END SUBROUTINE

  SUBROUTINE add_pool_r(sgn)
    IMPLICIT NONE
    REAL(8) sgn
    ndata(bin)=ndata(bin)+1
    work_r(1,bin)=work_r(1,bin)+sgn
    head_r=2
  END SUBROUTINE

  SUBROUTINE add_pool_z(sgn)
    IMPLICIT NONE
    COMPLEX(8) sgn
    ndata(bin)=ndata(bin)+1
    work_z(1,bin)=work_z(1,bin)+sgn
    head_z=2
  END SUBROUTINE

  SUBROUTINE rput_number(a)
    IMPLICIT NONE
    REAL(8) a
    work_r(head_r,bin)=work_r(head_r,bin)+a
    head_r=head_r+1
  END SUBROUTINE

  SUBROUTINE rput_array(n,a)
    IMPLICIT NONE
    INTEGER n
    REAL(8) a(*)
    work_r(head_r:head_r+n-1,bin)=work_r(head_r:head_r+n-1,bin)+a(1:n)
    head_r=head_r+n
  END SUBROUTINE

  SUBROUTINE rget_number(a)
    IMPLICIT NONE
    REAL(8) a
    a=work_r(head_r,1)
    head_r=head_r+1
  END SUBROUTINE

  SUBROUTINE rget_number_err(a,a_err)
    IMPLICIT NONE
    REAL(8) a,a_err
    a=work_r(head_r,1)
    a_err=work2_r(head_r)
    head_r=head_r+1
  END SUBROUTINE

  SUBROUTINE rget_array(n,a)
    IMPLICIT NONE
    INTEGER n
    REAL(8) a(*)
    a(1:n)=work_r(head_r:head_r+n-1,1)
    head_r=head_r+n
  END SUBROUTINE

  SUBROUTINE rget_array_err(n,a,a_err)
    IMPLICIT NONE
    INTEGER n
    REAL(8) a(*),a_err(*)
    a(1:n)=work_r(head_r:head_r+n-1,1)
    a_err(1:n)=work2_r(head_r:head_r+n-1)
    head_r=head_r+n
  END SUBROUTINE

  SUBROUTINE zput_number(a)
    IMPLICIT NONE
    COMPLEX(8) a
    work_z(head_z,bin)=work_z(head_z,bin)+a
    head_z=head_z+1
  END SUBROUTINE

  SUBROUTINE zput_array(n,a)
    IMPLICIT NONE
    INTEGER n
    COMPLEX(8) a(*)
    work_z(head_z:head_z+n-1,bin)=work_z(head_z:head_z+n-1,bin)+a(1:n)
    head_z=head_z+n
  END SUBROUTINE

  SUBROUTINE zget_number(a)
    IMPLICIT NONE
    COMPLEX(8) a
    a=work_z(head_z,1)
    head_z=head_z+1
  END SUBROUTINE

  SUBROUTINE zget_number_err(a,a_err)
    IMPLICIT NONE
    COMPLEX(8) a,a_err
    a=work_z(head_z,1)
    a_err=work2_z(head_z)
    head_z=head_z+1
  END SUBROUTINE

  SUBROUTINE zget_array(n,a)
    IMPLICIT NONE
    INTEGER n
    COMPLEX(8) a(*)
    a(1:n)=work_z(head_z:head_z+n-1,1)
    head_z=head_z+n
  END SUBROUTINE

  SUBROUTINE zget_array_err(n,a,a_err)
    IMPLICIT NONE
    INTEGER n
    COMPLEX(8) a(*),a_err(*)
    a(1:n)=work_z(head_z:head_z+n-1,1)
    a_err(1:n)=work2_z(head_z:head_z+n-1)
    head_z=head_z+n
  END SUBROUTINE

  SUBROUTINE average_pool()
    IMPLICIT NONE
    INTEGER b,nd,ierr

    ! do average for each bin
    DO b=1,nbin
      IF(size_r>0)THEN
        work_r(2:size_r,b)=work_r(2:size_r,b)/work_r(1,b)
        ! For sign average, only its absolute value has meaning because different cores have different initial boson fields with
        ! initial currentphase=(1.0,0.0)
        work_r(1,b)=abs(work_r(1,b)/ndata(b))
      END IF
      IF(size_z>0)THEN
        work_z(2:size_z,b)=work_z(2:size_z,b)/work_z(1,b)
        ! For sign average, only its absolute value has meaning because different cores have different initial boson fields with
        ! initial currentphase=(1.0,0.0)
        work_z(1,b)=abs(work_z(1,b)/ndata(b))
      END IF
    END DO

    ! do average between bins within each core
    IF(size_r>0)THEN
      ALLOCATE(work2_r(size_r))
      work2_r(:)=0d0
      DO b=1,nbin
        work2_r(:)=work2_r(:)+work_r(:,b)**2
      END DO
      work2_r(:)=work2_r(:)/nbin
      work_r(:,1)=sum(work_r(:,:),2)/nbin
    END IF
    IF(size_z>0)THEN
      ALLOCATE(work2_z(size_z))
      work2_z(:)=(0d0,0d0)
      DO b=1,nbin
        work2_z(:)=work2_z(:)+dcmplx(real(work_z(:,b))**2,aimag(work_z(:,b))**2)
      END DO
      work2_z(:)=work2_z(:)/nbin
      work_z(:,1)=sum(work_z(:,:),2)/nbin
    END IF

#ifdef MPI
    ! do average among all cores
    IF(size_r>0)THEN
      CALL dmpi_mean(work_r(:,1),size_r)
      CALL dmpi_mean(work2_r(:),size_r)
    END IF
    IF(size_z>0)THEN
      CALL zmpi_mean(work_z(:,1),size_z)
      CALL zmpi_mean(work2_z(:),size_z)
    END IF
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
#else
    nd=1
#endif

    ! get error bars
    IF(size_r>0)THEN
      work2_r(:)=sqrt((work2_r(:)-work_r(:,1)**2)/(nd*nbin-1))
    END IF
    IF(size_z>0)THEN
      work2_z(:)=dcmplx(sqrt((real(work2_z(:))-real(work_z(:,1))**2)/(nd*nbin-1)), &
      &                 sqrt((aimag(work2_z(:))-aimag(work_z(:,1))**2)/(nd*nbin-1)))
    END IF

    bin=1
    head_r=1
    head_z=1

  END SUBROUTINE

#ifdef MPI

  ! do average over different nodes. This is for real(8) data.
  SUBROUTINE dmpi_mean(x,n)  ! add x(n) to master node
    USE mpi
    IMPLICIT NONE
    INTEGER n,ierr,id,nd
    REAL(8) x(n),y(n)
    CALL mpi_barrier(mpi_comm_world,ierr)
    ! get summation of all nodes
    CALL mpi_reduce(x,y,n,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
    IF(id==0)x=y/nd  ! get the mean value
  END SUBROUTINE

  ! do average over different nodes. This is for complex(8) data.
  SUBROUTINE zmpi_mean(x,n)
    USE mpi
    IMPLICIT NONE
    INTEGER n,ierr,id,nd
    COMPLEX(8) x(n),y(n)
    CALL mpi_barrier(mpi_comm_world,ierr)
    CALL mpi_reduce(x,y,n,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
    IF(id==0)x=y/nd
  END SUBROUTINE

  ! do average over different nodes. This is for real(8) data.
  SUBROUTINE dmpi_mean_scalar(x)  ! add x(n) to master node
    USE mpi
    IMPLICIT NONE
    INTEGER n,ierr,id,nd
    REAL(8) x,y
    n=1
    CALL mpi_barrier(mpi_comm_world,ierr)
    ! get summation of all nodes
    CALL mpi_reduce(x,y,n,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
    IF(id==0)x=y/nd  ! get the mean value
  END SUBROUTINE

  ! do average over different nodes. This is for complex(8) data.
  SUBROUTINE zmpi_mean_scalar(x)
    USE mpi
    IMPLICIT NONE
    INTEGER n,ierr,id,nd
    COMPLEX(8) x,y
    n=1
    CALL mpi_barrier(mpi_comm_world,ierr)
    CALL mpi_reduce(x,y,n,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
    IF(id==0)x=y/nd
  END SUBROUTINE

#endif

END MODULE
