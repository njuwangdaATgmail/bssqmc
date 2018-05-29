MODULE hubhol
  
  !  interactions and phonon frequency
  REAL(8) U,Vph,debye

  ! exp( -gph_x2*phi(site,time)**2 -gph_p2*(phi(site,time)-phi(site,time+1))**2 )
  ! where gph_x2=1/2/dtau/Vph and gph_p2=1/2/dtau**3/Vph/debye**2
  REAL(8), ALLOCATABLE :: gph_x2(:),gph_p2(:)

END MODULE hubhol

! one can defines other types of local update here
SUBROUTINE generate_newfield_local(newfield,site,time,delta,ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: site,time,ifield
  COMPLEX(8), INTENT(OUT) :: newfield,delta(ndim_field(ifield),ndim_field(ifield),nflv)
  INTEGER newising,oldising
  IF(ifield<=nising)THEN
    oldising=nint(real(field(site,time,ifield)))
    CALL generate_newising_local(newising,oldising,delta,ifield)
    newfield=newising
  ELSE
    CALL generate_newphi_local(newfield,field(site,time,ifield),delta,ifield-nising)
  END IF
END SUBROUTINE

SUBROUTINE acceptprob_local(ratio,newfield,site,time,ifield,rtot)
  USE dqmc
  USE hubhol
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: site,time,ifield
  COMPLEX(8), INTENT(IN) :: ratio(nflv),newfield
  COMPLEX(8), INTENT(OUT) :: rtot
  INTEGER newising,oldising,flv
  COMPLEX(8) newphi,oldphi,diffnewphi,diffoldphi

  IF(ifield<=nising)THEN
    newising=nint(real(newfield))
    oldising=nint(real(field(site,time,ifield)))
    rtot=gamma_ising(newising,ifield)/gamma_ising(oldising,ifield)
    DO flv=1,nflv
      rtot=rtot*ratio(flv)
    END DO
  ELSE
    newphi=newfield
    oldphi=field(site,time,ifield-nising)
    rtot=exp(-gph_x2(ifield-nising)*(newphi**2-oldphi**2)) ! potential energy
    diffnewphi=newphi-field(site,mod(time,ntime)+1,ifield-nising)
    diffoldphi=oldphi-field(site,mod(time,ntime)+1,ifield-nising)
    rtot=rtot*exp(-gph_p2(ifield-nising)*(diffnewphi**2-diffoldphi**2)) ! kinetic energy on (time,time+1)
    diffnewphi=newphi-field(site,mod(time-2+ntime,ntime)+1,ifield-nising)
    diffoldphi=oldphi-field(site,mod(time-2+ntime,ntime)+1,ifield-nising)
    rtot=rtot*exp(-gph_p2(ifield-nising)*(diffnewphi**2-diffoldphi**2)) ! kinetic energy on (time,time-1)
    DO flv=1,nflv
      rtot=rtot*ratio(flv)
    END DO
    rtot=rtot*exp(-(newphi-oldphi))  ! ONLY used for PH symmetric case to fix <x>=0
  END IF

END SUBROUTINE

! one can define other types of global update here
SUBROUTINE generate_newfield_global(ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  IF(ifield<=nising)THEN
    CALL init_ising_random(ifield,ifield)
  ELSE
    CALL init_phi_random(ifield-nising,ifield)
  END IF
END SUBROUTINE

SUBROUTINE acceptprob_global(ratio,newfield,ifield,rtot)
  USE dqmc
  USE hubhol
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  COMPLEX(8), INTENT(IN) :: ratio(nflv),newfield(nsite,ntime)
  COMPLEX(8), INTENT(OUT) :: rtot
  INTEGER site,time,flv,newising,oldising
  COMPLEX(8) newphi,oldphi,diffnewphi,diffoldphi
  IF(ifield<=nising)THEN
    rtot=1d0
    DO site=1,nsite; IF(.not.mask_ising_site(site,ifield))CYCLE
      DO time=1,ntime
        newising=nint(real(newfield(site,time)))
        oldising=nint(real(field(site,time,ifield)))
        rtot=rtot*gamma_ising(newising,ifield)/gamma_ising(oldising,ifield)
      END DO
    END DO
    DO flv=1,nflv
      rtot=rtot*ratio(flv)
    END DO
  ELSE
    rtot=1d0
    DO site=1,nsite; IF(.not.mask_phi_site(site,ifield-nising))CYCLE
      DO time=1,ntime
        newphi=newfield(site,time)
        oldphi=field(site,time,ifield)
        diffnewphi=newfield(site,mod(time,ntime)+1)
        diffoldphi=field(site,mod(time,ntime)+1,ifield)
        rtot=rtot*exp( -gph_x2(ifield-nising)*(newphi**2-oldphi**2) ) & ! potential energy
          &        *exp( -gph_p2(ifield-nising)*(diffnewphi**2-diffoldphi**2) ) ! kinetic energy
        rtot=rtot*exp(-(newphi-oldphi))  ! ONLY used for PH symmetric case to fix <x>=0
      END DO
    END DO
    DO flv=1,nflv
      rtot=rtot*ratio(flv)
    END DO
  END IF
END SUBROUTINE

! one should identify whether ifield is ising-type or phi-type
! in default, ifield<=nising with ising-type
SUBROUTINE get_expV(site,time,ifield,flv,inv,expV)
  USE dqmc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: site,time,ifield,flv
  LOGICAL, INTENT(IN) :: inv
  COMPLEX(8), INTENT(OUT) :: expV(ndim_field(ifield),ndim_field(ifield))
  INTEGER newising
  COMPLEX(8) newphi
  IF(ifield<=nising)THEN
    newising=nint(real(field(site,time,ifield)))
    CALL get_expV_ising(newising,ifield,flv,inv,expV)
  ELSE
    newphi=field(site,time,ifield)
    CALL get_expV_phi(newphi,ifield-nising,flv,inv,expV)
  END IF
END SUBROUTINE


