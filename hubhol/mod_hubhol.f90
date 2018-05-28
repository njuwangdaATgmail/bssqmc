MODULE hubhol
  
  ! 
  REAL(8) U,V,debye

  !
  REAL(8) gph_x2,gph_p2


END MODULE hubhol

SUBROUTINE generate_newfield(newfield,site,time,delta,ifield)
  USE dqmc_complex, ONLY : ndim_field
  USE model_complex, ONLY : nising,generate_newising_random,generate_newphi_random
  IMPLICIT NONE
  INTEGER site,time,ifield
  COMPLEX(8) newfield,delta(ndim_field(ifield),ndim_field(ifield))
  IF(ifield<=nising)THEN
    CALL generate_newising_random(newising,site,time,delta,ifield)
    newfield=newising
  ELSE
    CALL generate_newphi_random(newfield,site,time,delta,ifield-nising)
  END IF
END SUBROUTINE

SUBROUTINE generate_newfield_global(ifield)
  !USE dqmc_complex
  !USE model_complex
  IMPLICIT NONE
  INTEGER ifield
END SUBROUTINE

SUBROUTINE acceptprob_local(ratio,newfield,site,time,ifield,rtot)
  USE dqmc_complex, ONLY : nflv
  !USE model_complex
  IMPLICIT NONE
  INTEGER site,time,ifield
  COMPLEX(8) ratio(nflv),newfield,rtot
END SUBROUTINE

SUBROUTINE acceptprob_global(ratio,newfield,ifield,rtot)
  USE dqmc_complex, ONLY : nflv,nsite,ntime
  !USE model_complex
  IMPLICIT NONE
  INTEGER ifield
  COMPLEX(8) ratio(nflv),newfield(nsite,ntime),rtot
END SUBROUTINE

SUBROUTINE get_expV(site,time,ifield,flv,inv,expV)
  USE dqmc_complex, ONLY : ndim_field
  !USE model_complex
  IMPLICIT NONE
  INTEGER site,time,ifield,flv
  LOGICAL inv
  COMPLEX(8) expV(ndim_field(ifield),ndim_field(ifield))
END SUBROUTINE

PROGRAM main
  USE dqmc_complex, ONLY : dqmc_driver
  IMPLICIT NONE
  CALL dqmc_driver()
END PROGRAM
