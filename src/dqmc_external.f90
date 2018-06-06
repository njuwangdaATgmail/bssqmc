! In this file, we define some NULL subroutines which are only needed
! to make an executable program dqmc.x
! Otherwise, this file will not be omitted to make a lib file dqmc.a
! and these subroutines should be given outside.
SUBROUTINE acceptprob_local_external(ratio,newfield,site,time,ifield,rtot)
  USE dqmc
  IMPLICIT NONE
  INTEGER site,time,ifield
  COMPLEX(8) ratio(nflv),newfield,rtot
END SUBROUTINE

SUBROUTINE generate_newfield_global_external(ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield
END SUBROUTINE

SUBROUTINE acceptprob_global_external(ratio,newfield,ifield,rtot)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield
  COMPLEX(8) ratio(nflv),newfield(nsite,ntime),rtot
END SUBROUTINE
