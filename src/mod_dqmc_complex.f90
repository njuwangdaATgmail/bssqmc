!-------------------------------------------------------------------
! This module is designed as the core of a general auxiliary field
! DQMC program. It can be used by externally calling
!   dqmc_driver().
! 
!-------------------------------------------------------------------
! List of external subroutines which should be specified by users:
!
! init()
! - the following variables are in need:
!   LOGICAL restart
!   LOGICAL proj
!   INTEGER nflv
!   INTEGER nsite
!   INTEGER nelec
!   INTEGER ntime
!   INTEGER nsp
!   INTEGER nbin
!   INTEGER nwarmup
!   INTEGER nmeasure
!   INTEGER ninterval
!   INTEGER ntmpout
!   INTEGER ngroup
!   INTEGER randomseed
!   INTEGER nfield
!   REAL(8) newMetro
!   COMPLEX(8), DIMENSION(:,:,:)  :: expk, expk_half, inv_expk, inv_expk_half
!   COMPLEX(8), DIMENSION(:,:,:)  :: slater, slater_Q, slater_R
!   COMPLEX(8), DIMENSION(:,:)    :: slater_D
!   COMPLEX(8), DIMENSION(:)      :: field
!   INTEGER,    DIMENSION(:)      :: ndim_field
!   LOGICAL,    DIMENSION(:)      :: mask_field
!   LOGICAL,    DIMENSION(:,:)    :: mask_field_site
!   INTEGER,    DIMENSION(:,:,:)  :: nb_field
!   INTEGER,    DIMENSION(:)      :: ninterval_global
! - set the measure pool by calling set_pool(nbin,poolsize_r,poolsize_z)
! - other initial settings can also be put here
!
! tmpout()
! - output temporarily
!
! measurement(time)
! - perform measurement
! - observable x can be measured by calling
!     put_pool(x)
!
! postprocess()
! - MC averaged observable x and its errbar can be read out by calling
!     get_pool(x,x_err)
! 
! generate_newfield(newfield,site,time,delta,ifield)
! - generate newfield at (site,time) and return
!     delta = exp((newfield-oldfield)*F) - 1
!   where F is the form factor of the auxiliary field
! 
! generate_newfield_global(ifield)
! - generate newfield globally
! 
! acceptprob_local(ratio,newfield,site,time,ifield,rtot)
! - evaluate acceptance ratio locally
!     rtot = ratio(1)*...*ratio(nflv) * W
!   where W accounts for free boson part
! 
! acceptprob_global(ratio,newfield,ifield,rtot)
! - the same as acceptprob_local but for global update
! 
! get_expV(site,time,ifield,flv,inv,expV)
! - return exp(V) = exp(field(site,time,flv)*F)
!   where F is the form factor of the auxiliary field
!
!------------------------------------------------------------------------

MODULE dqmc_complex

#ifdef MPI
  USE mpi
#endif
  USE randomlib
  USE matrixlib
  USE measpool

  !---------------------------------
  ! public scaler variables
  !--------------------------------

  ! restart mode or not
  LOGICAL :: restart                = .false.

  ! projector mode (T=0) or normal mode (T>0)
  LOGICAL :: proj                   = .false.

  ! number of fermion flavors
  INTEGER :: nflv                   = 1

  ! number of sites
  INTEGER :: nsite                  = -1

  ! number of filled electrons, only used for T=0
  INTEGER :: nelec                  = -1

  ! number of time slices
  INTEGER :: ntime                  = -1

  ! For T=0, do measurements between [ntime/2+1-nsp,ntime/2+1+nsp]
  ! For T>0, do measurements between [1,nsp]
  INTEGER :: nsp                    = 1

  ! number of data bins. The covariance is calculated between different bins.
  INTEGER :: nbin                   = 1

  ! warmup steps, in unit of Monte Carlo steps (MCS)
  INTEGER :: nwarmup                = 1000

  ! MCS in each bin
  INTEGER :: nmeasure               = 10000

  ! do measurement every ninterval MCSs
  INTEGER :: ninterval              = 1

  ! output temporarily every ntmpout MCSs
  INTEGER :: ntmpout                = 100

  ! update G from scratch every nscratch steps
  INTEGER :: nscratch               = 10

  ! every ngroup matrices are producted directly before performing QDR
  INTEGER :: ngroup                 = 10

  ! seed of random number generator
  INTEGER :: randomseed             = 20160701

  ! number of auxiliary fields
  INTEGER :: nfield                 = 0

  ! correction of Metropolis ratio
  REAL(8)    :: newMetro            = 0d0

  ! the phase of the current configuration
  COMPLEX(8) :: currentphase        = (1d0,0d0)

  ! the index of the current node
  INTEGER :: id                     = 0

  ! total number of all nodes
  INTEGER :: nd                     = 1

  !---------------------------------
  ! public arrays
  !--------------------------------

  ! exp(-dtau*K), exp(-dtau*K/2), exp(dtau*K), exp(dtau*K/2)
  ! dimension (nsite, nsite, nflv)
  COMPLEX(8), ALLOCATABLE :: expk(:,:,:), expk_half(:,:,:), inv_expk(:,:,:), inv_expk_half(:,:,:)

  ! trial wave function and its QDR decomposition
  ! slater = slater_Q * slater_D *slater_R
  ! dimension (nsite, nelec, nflv) for slater and slater_Q
  ! dimension (nelec, nflv) for slater_D
  ! dimension (nelec, nelec, nflv) for slater_R
  COMPLEX(8), ALLOCATABLE :: slater(:,:,:), slater_Q(:,:,:), slater_D(:,:), slater_R(:,:,:)

  ! equal-time Green's function
  ! dimension (nsite, nsite, nflv)
  COMPLEX(8), ALLOCATABLE :: g(:,:,:)

  ! configuration of auxiliary fields
  ! dimension (nsite, ntime, nfield)
  COMPLEX, ALLOCATABLE :: field(:,:,:)

  ! subspace dimension of each auxiliary field
  ! dimension (nfield)
  INTEGER, ALLOCATABLE :: ndim_field(:)

  ! whether the field exists
  ! dimension (nfield)
  LOGICAL, ALLOCATABLE :: mask_field(:)

  ! whether a site belongs to a given field
  ! dimension (nsite,nfield)
  LOGICAL, ALLOCATABLE :: mask_field_site(:,:)

  ! sites belong to a given field
  ! dimension (nsite, maxval(ndim_field), nfield)
  INTEGER, ALLOCATABLE :: nb_field(:,:,:)

  ! MCS distance to do global update
  ! dimension (nfield)
  INTEGER, ALLOCATABLE :: ninterval_global(:)

  !--------------------------------
  ! private scaler variables
  !--------------------------------

  ! difference between scratch and fast-update
  REAL(8), PRIVATE :: err_fast      = 0d0

  ! used to save starting time, in order to obtain running time
  REAL(8), PRIVATE :: t0            = maxexponent(1d0)

  ! start from a given bin
  INTEGER, PRIVATE :: ith_start              = 0

  ! start from a given MC step
  INTEGER, PRIVATE :: mcs_start              = 1

  ! whether update_scratch() has been finished in global update
  lOGICAL, PRIVATE :: scratch_global_useful  = .false.

  !-----------------------------------
  ! private arrays
  !----------------------------------

  ! string of B*B*...*B, see update_scratch() for details
  ! dimension (nsite, nsite/nelec, nblock, nflv) for Bstring_Q
  ! dimension (nsite/nelec, nblock, nflv) for Bstring_D
  ! dimension (nsite/nelec, nsite/nelec, nblock, nflv) for Bstring_T
  COMPLEX(8), ALLOCATABLE, PRIVATE :: Bstring_Q(:,:,:,:), Bstring_D(:,:,:), Bstring_T(:,:,:,:)

  ! total number of trials
  ! dimension (nfield)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_field(:), Ntotal_field_global(:)

  ! accepted number of trials
  ! dimension (nfield)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_field(:), Naccept_field_global(:)

  !----------------------------------
  ! list of public subroutines
  PUBLIC :: dqmc_driver, HS1, HS2, HSgeneral, evolve_left_K, evolve_right_K,  &
    &       evolve_left_V, evolve_right_V, evolve_left, evolve_right,         &
    &       evolve_left_2nd, evolve_right_2nd

  !----------------------------------
  ! list of private subroutines
  PRIVATE :: dqmc_update_local, dqmc_update_global, update_scratch_T0,        &
    &        update_scratch, init_, runtime_, tmpout_, loadtmp_, sort_

CONTAINS

  ! This is the entrance of this module. It should be called from outside.
  SUBROUTINE dqmc_driver()
    IMPLICIT NONE
    INTEGER ith,mcs,ierr,ifield

#ifdef MPI
    ! set MPI mode
    CALL mpi_init(ierr)
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
    IF(id==0)t0=mpi_wtime()  ! It seems mpi_wtime() may fail on very few machines.
#else
    ! set sequential mode
    nd=1
    id=0
    CALL cpu_time(t0)
#endif

    ! external initialization subroutine, set required parameters listed above
    CALL init()

    ! internal initialization, set remaining parameters
    CALL init_()

    DO ith=ith_start,nbin

      IF(ith==0)THEN

        ! do warmup
        DO mcs=mcs_start,nwarmup

          ! do global update
          DO ifield=1,nfield; IF(.not.mask_field(ifield))CYCLE
            IF(mod(mcs,ninterval_global(ifield))==0) CALL dqmc_update_global(ifield)
          END DO

          ! do local update without doing measurement
          CALL dqmc_update_local(.false.)

          ! output information to screen and save work sapce into hardware
          IF(mod(mcs,ntmpout)==0)THEN
            CALL tmpout_(ith,mcs)
            CALL tmpout()
          END IF

        END DO ! DO mcs=mcs_start,nwarmup

        ! reset mcs_start=1 in case the program starts from a break point with mcs_start>1
        mcs_start=1

      ELSE  ! IF(ith==0)

        ! do measurement
        DO mcs=mcs_start,nmeasure*ninterval

          ! do global update
          DO ifield=1,nfield; IF(.not.mask_field(ifield))CYCLE
            IF(mod(mcs,ninterval_global(ifield))==0) CALL dqmc_update_global(ifield)
          END DO

          ! only do measurement every ninterval space-time sweeps
          IF(mod(mcs,ninterval)==0)THEN
            CALL dqmc_update_local(.true.)
          ELSE
            CALL dqmc_update_local(.false.)
          END IF

          IF(mod(mcs,ntmpout)==0)THEN
            CALL tmpout_(ith,mcs)
            CALL tmpout()
          END IF

        END DO  ! DO mcs=mcs_start,nmeasure*ninterval

        ! reset mcs_start=1 in case the program starts from a break point with mcs_start>1
        mcs_start=1

        ! prepare the pool for the next data bin: move bin to bin+1 and reset head
        CALL move_pool()

      END IF  ! IF(ith==0)

    END DO  ! DO ith=ith_start,nbin

    ! do average of the measured observables in the pool
    CALL average_pool()

    ! external subroutine, data analysis and output
    CALL postprocess()

#ifdef MPI
    ! close MPI
    CALL mpi_barrier(mpi_comm_world,ierr)
    CALL mpi_finalize(ierr)
#endif

    IF(id==0)THEN
      PRINT*,'Job is done with',nd,' cores.'
    END IF

  END SUBROUTINE dqmc_driver

  ! internal initial subroutine, to set remaining parameters
  SUBROUTINE init_()
    IMPLICIT NONE
    INTEGER nblock

    ! set up the random number seed on different cores
    randomseed=randomseed+701703*id
    CALL init_rng(randomseed)

    ALLOCATE(g(nsite,nsite,nflv))

    ALLOCATE(Ntotal_field(nfield),Naccept_field(nfield))
    ALLOCATE(Ntotal_field_global(nfield),Naccept_field_global(nfield))
    Ntotal_field(:)=0
    Naccept_field(:)=0
    Ntotal_field_global(:)=0
    Naccept_field_global(:)=0

    nblock=ntime/nscratch
    IF(mod(ntime,nscratch).ne.0) nblock=nblock+1
    IF(proj)THEN
      ALLOCATE(Bstring_Q(nsite,nelec,nblock,nflv),Bstring_D(nelec,nblock,nflv),Bstring_T(nelec,nelec,nblock,nflv))
    ELSE
      ALLOCATE(Bstring_Q(nsite,nsite,nblock,nflv),Bstring_D(nsite,nblock,nflv),Bstring_T(nsite,nsite,nblock,nflv))
    END IF

    IF(restart) call loadtmp_()

  END SUBROUTINE init_

  ! internal subroutine to output running status to screen
  ! and save variables temporarily in case of unexpected breakdown (e.g. power off)
  SUBROUTINE tmpout_(ith,mcs)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ith,mcs
    INTEGER ifield
    CHARACTER*4 cid

    ! temporirally output the running status
    IF(id==0)THEN
      PRINT'(2i8,1a,1e13.6,1a,1e13.6,1a,1f15.2,1A,1e13.6)',ith,mcs,'     (', &
      &  real(currentphase),',',aimag(currentphase),')',runtime_(),'s     ',err_fast
      err_fast=0d0
      PRINT'(1a50,100f7.4)','acceptance rate of local update:', &
      & Naccept_field(:)*1d0/(Ntotal_field(:)+1d-30)
      PRINT'(1a50,100f7.4)','acceptance rate of global update:', &
      & Naccept_field_global(:)*1d0/(Ntotal_field_global(:)+1d-30)
      PRINT*,'field range:'
      DO ifield=1,nfield
        IF(mask_field(ifield)) PRINT'(1i4,1a46,4e15.4)',ifield,'  max(Re), min(Re), max(Im), min(Im):', &
        & maxval(real(field(:,:,ifield))),minval(real(field(:,:,ifield))), &
        & maxval(aimag(field(:,:,ifield))),minval(aimag(field(:,:,ifield)))
      END DO
      PRINT*,'-----------------------------------'
    END IF

    ! save QMC variables to files
    WRITE(cid,'(1i4)') id
    OPEN(81,FILE='tmpout_.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    WRITE(81) ith,mcs
    WRITE(81) currentphase
    WRITE(81) g
    WRITE(81) field,Ntotal_field,Naccept_field,Ntotal_field_global,Naccept_field_global
    CLOSE(81)

    IF(ith>0)CALL save_pool()

  END SUBROUTINE tmpout_

  ! internal subroutine to load variables from previous break point.
  SUBROUTINE loadtmp_()
    IMPLICIT NONE
    CHARACTER*4 cid

    WRITE(cid,'(1i4)') id
    OPEN(81,FILE='tmpout_.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    READ(81) ith_start,mcs_start
    READ(81) currentphase
    READ(81) g
    READ(81) field,Ntotal_field,Naccept_field,Ntotal_field_global,Naccept_field_global
    CLOSE(81)

    IF(ith_start>0)THEN
      CALL load_pool()
    END IF

    mcs_start=mcs_start+1

  END SUBROUTINE loadtmp_

  ! internal function to get the running time
  FUNCTION runtime_()
    IMPLICIT NONE
    REAL(8) runtime_,t1
#ifdef MPI
    IF(id==0)t1=mpi_wtime()
#else
    CALL cpu_time(t1)
#endif
    runtime_=t1-t0
  END FUNCTION runtime_

  ! do local update trying to update all auxiliary fields at each site and time
  SUBROUTINE dqmc_update_local(do_meas)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: do_meas
    INTEGER time,site,ifield,ndim,j,a,b,sitea,siteb,flv
    REAL(8) paccept
    COMPLEX(8) ratio(nflv),rtot,newfield,oldfield,delta_ja,gtmp(nsite,nsite)
    COMPLEX(8) delta(maxval(ndim_field),maxval(ndim_field),nflv)
    COMPLEX(8) ratiomat(maxval(ndim_field),maxval(ndim_field),nflv)

    DO time=1,ntime

      IF(mod(time-1,nscratch)==0)THEN

        ! update Green's function from scratch every nscratch steps
        IF(proj)THEN
          CALL update_scratch_T0(time)
        ELSE
          CALL update_scratch(time)
        END IF

      ELSE

        ! update Green's function by fast-update algorithm
        DO flv=1,nflv
          CALL evolve_left_K(g(:,:,flv),nsite,flv,.false.,.false.)
          CALL evolve_right_K(g(:,:,flv),nsite,flv,.true.,.false.)
          !g(:,:,flv)=matmul(matmul(expk(:,:,flv),g(:,:,flv)),inv_expk(:,:,flv))
        END DO

      END IF

      IF(proj)THEN

        ! at T=0, do measurement only during [ntime/2+1-nsp,ntime/2+1+nsp]
        IF(do_meas.and.abs(time-ntime/2-1)<=nsp)THEN
          ! prepare the pool for the following measurement
          CALL add_pool(currentphase)
          ! do measurement by calling external subroutine
          CALL measurement(time)
        END IF

      ELSE

        ! at T>0, do measurement only when time<=nsp
        IF(do_meas.and.time<=nsp)THEN
          ! prepare the pool for the following measurement
          CALL add_pool(currentphase)
          ! do measurement by calling external subroutine
          CALL measurement(time)
        END IF

      ENDIF

      ! update the i-th field locally
      DO ifield=1,nfield; IF(.not.mask_field(ifield))CYCLE

        ndim=ndim_field(ifield)

        ! try to update i-th field on each (block-)site
        DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE

          ! record old configuration
          oldfield=field(site,time,ifield)

          ! generate new configuration by calling an external subroutine which also returns delta=exp(V'-V)-1
          CALL generate_newfield(newfield,site,time,delta(1:ndim,1:ndim,1:nflv),ifield)

          ! calculate determinant ratio for each flavor
          IF(ndim==1)THEN
            DO flv=1,nflv
              ratiomat(1,1,flv)=1d0+(1d0-g(site,site,flv))*delta(1,1,flv)
              ratio(flv)=ratiomat(1,1,flv)
            END DO
          ELSE
            ratio=1d0
            DO flv=1,nflv
              DO a=1,ndim; sitea=nb_field(site,a,ifield)
                DO b=1,ndim; siteb=nb_field(site,b,ifield)
                  gtmp(a,b)=g(sitea,siteb,flv)
                END DO
              END DO
              ratiomat(1:ndim,1:ndim,flv)=delta(1:ndim,1:ndim,flv)-matmul(gtmp(1:ndim,1:ndim),delta(1:ndim,1:ndim,flv))
              DO a=1,ndim
                ratiomat(a,a,flv)=ratiomat(a,a,flv)+1d0
              END DO
              ratio(flv)=det(ndim,ratiomat(1:ndim,1:ndim,flv))
            END DO
          END IF

          ! calculate total accept probability externally
          CALL acceptprob_local(ratio,newfield,site,time,ifield,rtot)

          ! correction of Methopolis ratio
          IF(abs(rtot)<1d0)THEN
            paccept=abs(rtot)/(1d0+newMetro*abs(rtot))
          ELSE
            paccept=abs(rtot)/(newMetro+abs(rtot))
          END IF

          ! record the total number of tryings
          Ntotal_field(ifield)=Ntotal_field(ifield)+1

          ! accept the new configuration with probability paccept
          IF(drand()>paccept)CYCLE

          ! record the number of accepted tryings
          Naccept_field(ifield)=Naccept_field(ifield)+1

          ! obtain the phase of the current configuration
          currentphase=currentphase*rtot/abs(rtot)

          ! update Green's function based on Dyson equation
          IF(ndim==1)THEN
            DO flv=1,nflv
              gtmp=g(:,:,flv)
              DO j=1,nsite
                delta_ja=0d0
                IF(j==site)delta_ja=1d0
                gtmp(j,:)=gtmp(j,:)+(g(j,site,flv)-delta_ja)*delta(1,1,flv)/ratiomat(1,1,flv)*g(site,:,flv)
              END DO
              g(:,:,flv)=gtmp
            END DO
          ELSE
            DO flv=1,nflv
              CALL inverse(ndim,ratiomat(1:ndim,1:ndim,flv))
              delta(1:ndim,1:ndim,flv)=matmul(delta(1:ndim,1:ndim,flv),ratiomat(1:ndim,1:ndim,flv))
              gtmp=g(:,:,flv)
              DO j=1,nsite
                DO a=1,ndim; sitea=nb_field(site,a,ifield)
                  delta_ja=0d0; IF(j==sitea) delta_ja=1d0
                  DO b=1,ndim; siteb=nb_field(site,b,ifield)
                    gtmp(j,:)=gtmp(j,:)+(g(j,sitea,flv)-delta_ja)*delta(a,b,flv)*g(siteb,:,flv)
                  END DO
                END DO
              END DO
              g(:,:,flv)=gtmp
            END DO
          END IF

          ! update field
          field(site,time,ifield)=newfield

        END DO ! site-loop

        ! udpate Green's function by exp(phi*F)*G*exp(-phi*F)
        DO flv=1,nflv
          CALL evolve_left_V(time,ifield,g(:,:,flv),nsite,flv,.false.)
          CALL evolve_right_V(time,ifield,g(:,:,flv),nsite,flv,.true.)
        END DO

      END DO ! ifield-loop

    END DO ! time-loop

  END SUBROUTINE dqmc_update_local

  ! do global update
  SUBROUTINE dqmc_update_global(ifield)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifield
    INTEGER flv,i
    INTEGER, ALLOCATABLE :: ising_old(:,:),ising_new(:,:)
    REAL(8) paccept
    COMPLEX(8) oldfield(nsite,ntime),newfield(nsite,ntime)
    COMPLEX(8) ratio(nflv),rtot
    COMPLEX(8), ALLOCATABLE :: dvec(:,:),dvec_old(:,:)


    ! calculate old determinant from scratch
    IF(proj)THEN
      CALL update_scratch_T0(1)
    ELSE
      CALL update_scratch(1)
    END IF

    ! save eigenvalues of the old determinant
    DO flv=1,nflv

      IF(proj)THEN

        ratio(flv)=1d0/det(nelec,Bstring_Q(1:nelec,1:nelec,1,flv))
        dvec_old(1:nelec,flv)=Bstring_D(1:nelec,1,flv)
        CALL sort_(nelec,dvec_old(1:nelec,flv))

      ELSE

        ratio(flv)=1d0/det(nsite,Bstring_Q(1:nsite,1:nsite,1,flv))
        dvec_old(1:nsite,flv)=Bstring_D(1:nsite,1,flv)
        CALL sort_(nsite,dvec_old(1:nsite,flv))

      END IF

    END DO

    ! save old field
    oldfield=field(:,:,ifield)

    ! generate new field
    CALL generate_newfield_global(ifield)

    ! calculate new determinant from scratch
    IF(proj)THEN
      CALL update_scratch_T0(1)
    ELSE
      CALL update_scratch(1)
    END IF

    ! get eigenvalues of the new determinant
    DO flv=1,nflv

      IF(proj)THEN

        ratio(flv)=ratio(flv)*det(nelec,Bstring_Q(1:nelec,1:nelec,1,flv))
        dvec(1:nelec,flv)=Bstring_D(1:nelec,1,flv)
        CALL sort_(nelec,dvec(1:nelec,flv))

        DO i=1,nelec
          ratio(flv)=ratio(flv)*dvec(i,flv)/dvec_old(i,flv)
        END DO

      ELSE

        ratio(flv)=ratio(flv)*det(nsite,Bstring_Q(1:nsite,1:nsite,1,flv))
        dvec(1:nsite,flv)=Bstring_D(1:nsite,1,flv)
        CALL sort_(nsite,dvec(1:nsite,flv))

        DO i=1,nsite
          ratio(flv)=ratio(flv)*dvec(i,flv)/dvec_old(i,flv)
        END DO

      END IF

    END DO

    ! get acceptance ratio
    newfield=field(:,:,ifield)
    field(:,:,ifield)=oldfield
    CALL acceptprob_global(ratio,newfield,ifield,rtot)

    ! correction of Metropolis ratio
    IF(abs(rtot)<1d0)THEN
      paccept=abs(rtot)/(1d0+newMetro*abs(rtot))
    ELSE
      paccept=abs(rtot)/(newMetro+abs(rtot))
    END IF

    ! record the total number of tryings
    Ntotal_field_global(ifield)=Ntotal_field_global(ifield)+1

    ! whether update_scratch is useful
    scratch_global_useful=.false.

    ! accept the new configuration with probability paccept
    IF(drand()>paccept)RETURN

    ! record the number of accepted tryings
    Naccept_field_global(ifield)=Naccept_field_global(ifield)+1

    ! obtain the phase of the current configuration
    currentphase=currentphase*rtot/abs(rtot)

    ! update field
    field(:,:,ifield)=newfield

    ! whether update_scratch useful
    scratch_global_useful=.true.

  END SUBROUTINE dqmc_update_global

  ! sort the array a(n) from large to small absulute values
  SUBROUTINE sort_(n,a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER i,j
    COMPLEX(8) a(n),tmp
    DO j=n-1,1,-1
      DO i=1,j
        IF(abs(a(i))<abs(a(i+1)))THEN
          tmp=a(i);a(i)=a(i+1);a(i+1)=tmp
        END IF
      END DO
    END DO
  END SUBROUTINE sort_

  !------------------------------------------------------------------------------------------------
  ! evaluate Green's function from definition directly
  ! Bstring technique is used to accelerate (a cheap realization, may not be the most efficient)
  !------------------------------------------------------------------------------------------------
  ! For T>0 algorithm:
  ! define Bstring(nsite,nsite,nblock)
  ! when time=1 or block=1,
  !   used  :   none
  !   update:   Bstring(:,:,nblock)=BBB(nblock)
  !             Bstring(:,:,nblock-1)=BBB(nblock)*BBB(nblock-1)
  !             ...
  !             Bstring(:,:,2)=BBB(nblock)*...*BBB(2)
  !             Bstring(:,:,1)=BBB(nblock)*...*BBB(1)+1 which further gives the determinant
  ! when time=nscratch+1 or block=2,
  !   used  :   Bstring(:,:,block)
  !   update:   Bstring(:,:,1)=BBB(1)
  !
  ! when nblock>block>2,
  !   used  :   Bstring(:,:,block)
  !             Bstring(:,:,block-2)
  !   update:   Bstring(:,:,block-1)=BBB(block-1)*Bstring(block-2)
  ! ...
  ! when block=nblock
  !   used  :   Bstring(:,:,nblock)
  !             Bstring(:,:,nblock-2)
  !   update:   Bstring(:,:,1)=BBB+1
  !-----------------------------------------------------------------------------------------------
  ! For T=0 algorithm:
  ! define Bstring(nsite,nelec,nblock)
  ! when time=1 or block=1,
  !   used  :   none
  !   update:   Bstring(:,:,nblock)=transpose( P'*BBB(nblock) )
  !             Bstring(:,:,nblock-1)=transpose( P'*BBB(nblock)*BBB(nblock-1) )
  !             ...
  !             Bstring(:,:,2)=transpose( P'*BBB(nblock)*...*BBB(2) )
  !             Bstring(1:nelec,:,1)=P'*BBB(nblock)*...*BBB(1)*P which further gives the determinant
  ! when time=nscratch+1 or block=2,
  !   used  :   Bstring(:,:,block)
  !   update:   Bstring(:,:,1)=BBB(1)*P
  !
  ! when nblock>block>2,
  !   used  :   Bstring(:,:,block)
  !             Bstring(:,:,block-2)
  !   update:   Bstring(:,:,block-1)=BBB(block-1)*Bstring(block-2)
  ! ...
  ! when block=nblock
  !   used  :   Bstring(:,:,nblock)
  !             Bstring(:,:,nblock-2)
  !   update:   Bstring(:,:,1)=P'BBB*P
  !-------------------------------------------------------------------------------------------------

  ! calculate T=0 Green's function from definition, using QDR decomposition stabilization algorithm.
  SUBROUTINE update_scratch_T0(time)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time
    INTEGER block,pblock,flv,flag,p,i
    COMPLEX(8) dvec(nelec),tri(nelec,nelec),tmat(nelec,nelec),gfast(nsite,nsite)
    COMPLEX(8), ALLOCATABLE :: qmat(:,:)

    IF(scratch_global_useful)THEN
      scratch_global_useful=.false.
      RETURN
    END IF

    ! get the block index: [1,nscratch]->1, [nscratch+1,2*nscratch]->2,...
    ! we come here only when mod(time-1,nscratch)==0
    block=(time-1)/nscratch+1

    DO flv=1,nflv

      ! save g to evaluate the error of fast-update
      gfast=g(:,:,flv)

      IF(block==1)THEN

        ALLOCATE(qmat(nelec,nsite))

        ! tmat*dvec*qmat=P'
        qmat=conjg(transpose(slater_Q(:,:,flv)))
        dvec=conjg(slater_D(:,flv))
        tmat=conjg(transpose(slater_R(:,:,flv)))

        flag=0
        DO p=ntime,1,-1
          CALL evolve_right(p,qmat,nelec,flv,.false.)
          flag=flag+1
          IF(flag==ngroup.or.mod(p-1,nscratch)==0)THEN
            DO i=1,nsite
              qmat(:,i)=dvec(:)*qmat(:,i)
            END DO
            CALL ldq(nelec,nsite,qmat,tri,dvec)
            tmat=matmul(tmat,tri)
            flag=0
          END IF
          IF(mod(p-1,nscratch)==0.and.p>1)THEN
            pblock=(p-1)/nscratch+1
            Bstring_Q(:,:,pblock,flv)=transpose(qmat)
            Bstring_D(:,pblock,flv)=dvec
            Bstring_T(:,:,pblock,flv)=tmat
          END IF
        END DO

        !   g = 1 - Q_R * (Q_L*Q_R)^{-1} * Q_L
        tri=matmul(qmat,slater_Q(:,:,flv))
        CALL inverse(nelec,tri)
        g(:,:,flv)=-matmul(matmul(slater_Q(:,:,flv),tri),qmat)
        DO i=1,nsite
          g(i,i,flv)=g(i,i,flv)+1d0
        END DO

        ! Bstring_T(1)*Bstring_D(1)*Bstring_Q(1) = P'*BBB*P
        DO i=1,nsite
          qmat(:,i)=dvec(:)*qmat(:,i)
        END DO
        qmat(1:nelec,1:nelec)=matmul(qmat,slater(:,:,flv))
        CALL ldq(nelec,nelec,qmat(1:nelec,1:nelec),tri,dvec)
        Bstring_Q(1:nelec,1:nelec,1,flv)=qmat(1:nelec,1:nelec)
        Bstring_D(:,1,flv)=dvec
        Bstring_T(:,:,1,flv)=matmul(tmat,tri)

        DEALLOCATE(qmat)

      ELSE

        ALLOCATE(qmat(nsite,nelec))

        IF(block==2)THEN
          ! qmat*dvec*tmat = slater
          qmat=slater_Q(:,:,flv)
          dvec=slater_D(:,flv)
          tmat=slater_R(:,:,flv)
        ELSE
          ! qmat*dvec*tmat = Bstring(block-2)
          qmat=Bstring_Q(:,:,block-2,flv)
          dvec=Bstring_D(:,block-2,flv)
          tmat=Bstring_T(:,:,block-2,flv)
        END IF

        flag=0
        DO p=(block-2)*nscratch+1,(block-1)*nscratch
          CALL evolve_left(p,qmat,nelec,flv,.false.)
          flag=flag+1
          IF(flag==ngroup.or.mod(p,nscratch)==0)THEN
            DO i=1,nelec
              qmat(:,i)=qmat(:,i)*dvec(i)
            END DO
            CALL qdr(nsite,nelec,qmat,tri,dvec)
            tmat=matmul(tri,tmat)
            flag=0
          END IF
        END DO
        Bstring_Q(:,:,block-1,flv)=qmat
        Bstring_D(:,block-1,flv)=dvec
        Bstring_T(:,:,block-1,flv)=tmat

        ! g = 1 - qmat*dvec*tmat * (BL_T*BL_D*BL_Q*qmat*dvec*tmat)^{-1} * BL_T*BL_D*BL_Q
        !   = 1 - qmat * (BL_Q*qmat)^{-1} * BL_Q

        tri=matmul(transpose(Bstring_Q(:,:,block,flv)),qmat)
        CALL inverse(nelec,tri)
        g(:,:,flv)=-matmul(matmul(qmat,tri),transpose(Bstring_Q(:,:,block,flv)))
        DO i=1,nsite
          g(i,i,flv)=g(i,i,flv)+1d0
        END DO

        DEALLOCATE(qmat)

      END IF

      IF(time>1)THEN
        CALL evolve_left_K(gfast,nsite,flv,.false.,.false.)
        CALL evolve_right_K(gfast,nsite,flv,.true.,.false.)
        err_fast=max(maxval(abs(gfast-g(:,:,flv))),err_fast)
      END IF

    END DO

  END SUBROUTINE update_scratch_T0

  ! calculate T>0 Green's function from definition, using QDR decomposition stabilization algorithm.
  SUBROUTINE update_scratch(time)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time
    INTEGER block,pblock,flv,flag,p,i
    COMPLEX(8) dvec(nsite),qmat(nsite,nsite),tmat(nsite,nsite),tri(nsite,nsite),gfast(nsite,nsite)

    IF(scratch_global_useful)THEN
      scratch_global_useful=.false.
      RETURN
    END IF

    ! get the block index: [1,nscratch]->1, [nscratch+1,2*nscratch]->2,...
    ! we come here only when mod(time-1,nscratch)==0
    block=(time-1)/nscratch+1

    DO flv=1,nflv

      gfast=g(:,:,flv)

      IF(block==1)THEN

        ! tmat*dvec*qmat=1
        qmat=0d0
        dvec=1d0
        tmat=0d0
        DO i=1,nsite
          qmat(i,i)=1d0
          tmat(i,i)=1d0
        END DO

        flag=0
        DO p=ntime,1,-1
          CALL evolve_right(p,qmat,nsite,flv,.false.)
          flag=flag+1
          IF(flag==ngroup.or.mod(p-1,nscratch)==0)THEN
            DO i=1,nsite
              qmat(:,i)=dvec(:)*qmat(:,i)
            END DO
            CALL ldq(nsite,nsite,qmat,tri,dvec)
            tmat=matmul(tmat,tri)
            flag=0
          END IF
          IF(mod(p-1,nscratch)==0.and.p>1)THEN
            pblock=(p-1)/nscratch+1
            Bstring_Q(:,:,pblock,flv)=qmat
            Bstring_D(:,pblock,flv)=dvec
            Bstring_T(:,:,pblock,flv)=tmat
          END IF
        END DO

        ! tmat*dvec*qmat = BBB
        ! g^{-1}= BBB + 1 = ( tmat*dvec + qmat^{-1} ) * qmat

        Bstring_Q(:,:,1,flv)=qmat
        CALL inverse(nsite,qmat)
        g(:,:,flv)=qmat
        DO i=1,nsite
          qmat(:,i)=qmat(:,i)+tmat(:,i)*dvec(i)
        END DO
        CALL ldq(nsite,nsite,qmat,tri,dvec)
        Bstring_Q(:,:,1,flv)=matmul(qmat,Bstring_Q(:,:,1,flv))
        Bstring_D(:,1,flv)=dvec
        Bstring_T(:,:,1,flv)=tri
        CALL inverse(nsite,tri)
        CALL inverse(nsite,qmat)
        DO i=1,nsite
          tri(:,i)=tri(:,i)/dvec(:)
        END DO
        g(:,:,flv)=matmul(matmul(g(:,:,flv),qmat),tri)

      ELSE

        IF(block==2)THEN
          ! tmat*dvec*qmat=1
          qmat=0d0
          dvec=1d0
          tmat=0d0
          DO i=1,nsite
            qmat(i,i)=1d0
            tmat(i,i)=1d0
          END DO
        ELSE
          ! tmat*dvec*qmat=Bstring(block-2)
          tmat=Bstring_T(:,:,block-2,flv)
          qmat=Bstring_Q(:,:,block-2,flv)
          dvec=Bstring_D(:,block-2,flv)
        END IF

        flag=0
        DO p=(block-2)*nscratch+1,(block-1)*nscratch
          CALL evolve_left(p,qmat,nsite,flv,.false.)
          flag=flag+1
          IF(flag==ngroup.or.mod(p,nscratch)==0)THEN
            DO i=1,nsite
              qmat(:,i)=qmat(:,i)*dvec(i)
            END DO
            CALL qdr(nsite,nsite,qmat,tri,dvec)
            tmat=matmul(tri,tmat)
            flag=0
          END IF
        END DO
        Bstring_Q(:,:,block-1,flv)=qmat
        Bstring_D(:,block-1,flv)=dvec
        Bstring_T(:,:,block-1,flv)=tmat

        ! g^{-1} = qmat*dvec*tmat*BL_T*BL_D_BL_Q + 1

        tmat=matmul(tmat,Bstring_T(:,:,block,flv))
        DO i=1,nsite
          tmat(:,i)=dvec(:)*tmat(:,i)*Bstring_D(i,block,flv)
        END DO
        CALL qdr(nsite,nsite,tmat,tri,dvec)
        qmat=matmul(qmat,tmat)
        tri=matmul(tri,Bstring_Q(:,:,block,flv))

        ! g^{-1} = qmat*dvec*tri + 1

        CALL inverse(nsite,qmat)
        g(:,:,flv)=qmat
        DO i=1,nsite
          qmat(:,i)=qmat(:,i)+dvec(:)*tri(:,i)
        END DO
        CALL qdr(nsite,nsite,qmat,tri,dvec)
        CALL inverse(nsite,qmat)
        CALL inverse(nsite,tri)
        DO i=1,nsite
          tri(:,i)=tri(:,i)/dvec(i)
        END DO
        g(:,:,flv)=matmul(tri,matmul(qmat,g(:,:,flv)))

      END IF

      ! when time=1, we may come through a global update, then the following fast update has no meaning
      IF(time>1)THEN
        CALL evolve_left_K(gfast,nsite,flv,.false.,.false.)
        CALL evolve_right_K(gfast,nsite,flv,.true.,.false.)
        err_fast=max(maxval(abs(gfast-g(:,:,flv))),err_fast)
      END IF

    END DO

  END SUBROUTINE update_scratch

  !
  !------------------------------------------------------------------------------
  !  The following 8 subroutines are provided to perform evolutions of a matrix.
  !       evolve_left_K   : exp(K)*matrix
  !       evolve_right_K  : matrix*exp(K)
  !       evolve_left_V   : exp(V)*matrix
  !       evolve_right_V  : matrix*exp(V)
  !       evolve_left     : B*matrix
  !       evolve_right    : matrix*B
  !       evolve_left_2nd : B_2nd*matrix
  !       evolve_right_2nd: matrix*B_2nd
  !------------------------------------------------------------------------------
  !
  ! exp( +- (dtau, dtau/2)*K ) * matrix
  SUBROUTINE evolve_left_K(matrix,d,flv,inv,half)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: d,flv
    COMPLEX(8) matrix(nsite,d)
    LOGICAL, INTENT(IN) :: inv,half
    IF(inv.and.half)THEN
      matrix=matmul(inv_expk_half(:,:,flv),matrix)
    ELSE IF(.not.inv.and.half)THEN
      matrix=matmul(expk_half(:,:,flv),matrix)
    ELSE IF(inv.and..not.half)THEN
      matrix=matmul(inv_expk(:,:,flv),matrix)
    ELSE
      matrix=matmul(expk(:,:,flv),matrix)
    END IF
  END SUBROUTINE evolve_left_K

  ! matrix * exp( +- (dtau, dtau/2)*K )
  SUBROUTINE evolve_right_K(matrix,d,flv,inv,half)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: d,flv
    COMPLEX(8) matrix(d,nsite)
    LOGICAL, INTENT(IN) :: inv,half
    IF(inv.and.half)THEN
      matrix=matmul(matrix,inv_expk_half(:,:,flv))
    ELSE IF(.not.inv.and.half)THEN
      matrix=matmul(matrix,expk_half(:,:,flv))
    ELSE IF(inv.and..not.half)THEN
      matrix=matmul(matrix,inv_expk(:,:,flv))
    ELSE
      matrix=matmul(matrix,expk(:,:,flv))
    END IF
  END SUBROUTINE evolve_right_K

  ! exp( +- V) * matrix
  SUBROUTINE evolve_left_V(time,ifield,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,ifield,d,flv
    LOGICAL, INTENT(IN) :: inv
    INTEGER a,b,site,sitea,siteb
    COMPLEX(8) matrix(nsite,d),expV(ndim_field(ifield),ndim_field(ifield))
    COMPLEX(8) gtmp(ndim_field(ifield),d)
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      CALL get_expV(site,time,ifield,flv,inv,expV)
      IF(ndim_field(ifield)==1)THEN
        matrix(site,:)=expV(1,1)*matrix(site,:)
      ELSE
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          gtmp(a,:)=matrix(sitea,:)
        END DO
        gtmp=matmul(expV,gtmp)
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          matrix(sitea,:)=gtmp(a,:)
        END DO
      END IF
    END DO
  END SUBROUTINE evolve_left_V

  ! matrix * exp( +- V)
  SUBROUTINE evolve_right_V(time,ifield,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,ifield,d,flv
    LOGICAL, INTENT(IN) :: inv
    INTEGER a,b,site,sitea,siteb
    COMPLEX(8) matrix(d,nsite),expV(ndim_field(ifield),ndim_field(ifield))
    COMPLEX(8) gtmp(d,ndim_field(ifield))
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      CALL get_expV(site,time,ifield,flv,inv,expV)
      IF(ndim_field(ifield)==1)THEN
        matrix(:,site)=matrix(:,site)*expV(1,1)
      ELSE
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          gtmp(:,a)=matrix(:,sitea)
        END DO
        gtmp=matmul(gtmp,expV)
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          matrix(:,sitea)=gtmp(:,a)
        END DO
      END IF
    END DO
  END SUBROUTINE evolve_right_V

  ! B(time)*matrix
  SUBROUTINE evolve_left(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(nsite,d)
    INTEGER ifield
    DO ifield=1,nfield
      CALL evolve_left_V(time,ifield,matrix,d,flv,inv)
    END DO
    CALL evolve_left_K(matrix,d,flv,inv,.false.)
  END SUBROUTINE evolve_left

  ! matrix*B(time)
  SUBROUTINE evolve_right(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(d,nsite)
    INTEGER ifield
    CALL evolve_right_K(matrix,d,flv,inv,.false.)
    DO ifield=nfield,1,-1
      CALL evolve_right_V(time,ifield,matrix,d,flv,inv)
    END DO
  END SUBROUTINE evolve_right

  ! B_2nd(time)*matrix
  SUBROUTINE evolve_left_2nd(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(nsite,d)
    INTEGER ifield
    CALL evolve_left_K(matrix,d,flv,inv,.true.)
    DO ifield=1,nfield
      CALL evolve_left_V(time,ifield,matrix,d,flv,inv)
    END DO
    CALL evolve_left_K(matrix,d,flv,inv,.true.)
  END SUBROUTINE evolve_left_2nd

  ! matrix*B_2nd(time)
  SUBROUTINE evolve_right_2nd(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(d,nsite)
    INTEGER ifield
    CALL evolve_right_K(matrix,d,flv,inv,.true.)
    DO ifield=nfield,1,-1
      CALL evolve_right_V(time,ifield,matrix,d,flv,inv)
    END DO
    CALL evolve_right_K(matrix,d,flv,inv,.true.)
  END SUBROUTINE evolve_right_2nd

  !
  !--------------------------------------------------------------------
  ! several usually used discrete Hubbard-Stratonovich transformation
  ! HS1:
  !   exp(-dtau*g*op^2) = ga*exp(lam*op) + ga*exp(-lam*op)
  ! HS2, HSgeneral:
  !   exp(-dtau*g*op^2) = ga1*exp(lam1*op) + ga1*exp(-lam1*op)
  !                     + ga2*exp(lam2*op) + ga2*exp(-lam2*op)
  !--------------------------------------------------------------------
  !
  ! |eig(op)|=0,1, exact
  SUBROUTINE HS1(ga,lam,a)  ! a=exp(-dtau*g)
    IMPLICIT NONE
    COMPLEX(8) ga,lam
    REAL(8) a
    ga=0.5d0
    IF(a>1d0)THEN
      lam=acosh(a)
    ELSE
      lam=dcmplx(0d0,acos(a))
    END IF
  END SUBROUTINE HS1

  ! |eig(op)|=0,1,2,3, exact
  SUBROUTINE HS2(ga1,lam1,ga2,lam2,a)  ! a=exp(-dtau*g)
    IMPLICIT NONE
    COMPLEX(8) ga1,ga2,lam1,lam2
    REAL(8) a,d
    d=sqrt(8d0+a**2*(3d0+a**2)**2)
    ga1=(-a*(3d0+a**2)+d)/(4*d)
    ga2=(a*(3d0+a**2)+d)/(4*d)
    IF(a>1d0)THEN
      lam1=acosh((a+2*a**3+a**5+(a**2-1d0)*d)/4)
      lam2=acosh((a+2*a**3+a**5-(a**2-1d0)*d)/4)
    ELSE
      lam1=DCMPLX(0d0,acos((a+2*a**3+a**5+(a**2-1d0)*d)/4))
      lam2=DCMPLX(0d0,acos((a+2*a**3+a**5-(a**2-1d0)*d)/4))
    END IF
  END SUBROUTINE HS2

  ! |eig(op)|=any value, but approximate.
  SUBROUTINE HSgeneral(ga1,lam1,ga2,lam2,a) ! a=-dtau*g
    IMPLICIT NONE
    COMPLEX(8) ga1,lam1,ga2,lam2
    REAL(8) a
    ga1=(3d0+sqrt(6d0))/12d0
    ga2=(3d0-sqrt(6d0))/12d0
    IF(a>0d0)THEN
      lam1=sqrt(a*2*(3d0-sqrt(6d0)))
      lam2=sqrt(a*2*(3d0+sqrt(6d0)))
    ELSE
      lam1=DCMPLX(0d0,sqrt(-a*2*(3d0-sqrt(6d0))))
      lam2=DCMPLX(0d0,sqrt(-a*2*(3d0+sqrt(6d0))))
    END IF
  END SUBROUTINE HSgeneral

END MODULE dqmc_complex


