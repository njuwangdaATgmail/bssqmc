! =============================================================================
! module of random number generator, it has four interfaces:
!	 FUNCTION irand(n) returns a random integer from 0 to n-1
!  FUNCTION drand() returns a random double precision number from 0 to 1
!  FUNCTION drand_sym() returns a random double precision number from -1 to 1
!  SUBROUTINE init_rng() is used to initialize the generator by system clock
!  SUBROUTINE init_rng(iseed) is used to initialize the generator by input iseed
! =============================================================================
MODULE randomlib

  INTEGER,PARAMETER,PRIVATE :: dp=8
  INTERFACE init_rng
    MODULE PROCEDURE init_rng,init_rng_by_seed
  END INTERFACE

  CONTAINS

  INTEGER FUNCTION irand(n)
    IMPLICIT NONE
    REAL(dp) r
    INTEGER n
    CALL RANDOM_NUMBER(r)
    irand = INT(r*n)
  END FUNCTION

  REAL(dp) FUNCTION drand()
    CALL RANDOM_NUMBER(drand)
  END FUNCTION

  REAL(dp) FUNCTION drand_sym()
    CALL RANDOM_NUMBER(drand_sym)
    drand_sym=2*drand_sym-1d0
  END FUNCTION

  SUBROUTINE  init_rng()
    INTEGER :: i,n,clock
    INTEGER,DIMENSION(:),ALLOCATABLE :: seed
    CALL RANDOM_SEED(size=n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(count=clock)
    seed=clock+37*(/(i-1,i=1,n)/)
    CALL RANDOM_SEED(put=seed)
    DEALLOCATE(seed)
  END SUBROUTINE

  SUBROUTINE init_rng_by_seed(iseed)
    INTEGER iseed,n
    INTEGER,DIMENSION(:),ALLOCATABLE :: seed
    CALL RANDOM_SEED(size=n)
    ALLOCATE(seed(n))
    seed=iseed
    CALL RANDOM_SEED(put=seed)
    DEALLOCATE(seed)
  END SUBROUTINE

END MODULE randomlib
