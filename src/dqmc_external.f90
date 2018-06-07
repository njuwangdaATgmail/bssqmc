!--------------------------------------------------------------------
! This file provides some NULL subroutines which are only used
! to make an executable program dqmc.x
! Otherwise, this file will be omitted to make a lib file libdqmc.a
! and the following subroutines must be given outside.
!--------------------------------------------------------------------

! if |type_field|>=100
SUBROUTINE set_field_external(ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield
END SUBROUTINE

! if |type_field|>=100
SUBROUTINE acceptprob_local_external(ratio,newfield,site,time,ifield,rtot)
  USE dqmc
  IMPLICIT NONE
  INTEGER site,time,ifield
  COMPLEX(8) ratio(nflv),newfield,rtot
END SUBROUTINE

! if global_method==100
SUBROUTINE generate_newfield_global_external(ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield
END SUBROUTINE

! if |type_field|>=100
SUBROUTINE acceptprob_global_external(ratio,newfield,ifield,rtot)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield
  COMPLEX(8) ratio(nflv),newfield(nsite,ntime),rtot
END SUBROUTINE

! if do_measure_external=.true.
SUBROUTINE measurement_external(time)
  USE dqmc
  IMPLICIT NONE
  INTEGER time
END SUBROUTINE

! if do_tmpout_external=.true.
SUBROUTINE tmpout_external()
  USE dqmc
  IMPLICIT NONE
END SUBROUTINE

! if do_postprocess_external=.true.
SUBROUTINE postprocess_external()
  USE dqmc
  IMPLICIT NONE
END SUBROUTINE
