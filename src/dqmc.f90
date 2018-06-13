!---------------------------------------------------------------------------
! This is a wrapper of the dqmc_core. It defines a 3D (1-2D are special cases)
! fermionic lattice coupling to boson fields. The boson fields can be either
! descrete (called ising) or continuous (called phi).
!
! NOTICE: In this module, uniform fermion-fermion and fermion-boson
! interactions are supposed. For e.g. disorder case, we can simply add
! site-dependence to relevant variables, e.g. lam_field, gamma_field, etc.
!
! Several useful subroutines are provided:
! to be added...
!---------------------------------------------------------------------------
! Todo:
! + aobut initialization: realize a general UVJ-phonon-phi4 model and
!   provide relevant subroutines to realize any new boson fields.
! + about measurement: all measurements are correlation functions
!   (1) single-particle:  [(x1,y1,z1,orb1), (x2,y2,z2,orb2)]
!   (2) two-particle:     [(x1,y1,z1,orb1), (x1+dx,y1+dy,z1+dz,orb1')] -->
!       (both PP and PH)  [(x2,y2,z2,orb2), (x2+dx,y2+dy,z2+dz,orb2')]
!                         with form factor fmat_meas(:,:)
!      so the user only need provide the form factor and its subspace
!---------------------------------------------------------------------------
!
! the boson fields are redefined with no need to distinguish ising and phi
! isingmax=0 is used to identify continuous boson fields
MODULE dqmc

  USE dqmc_core

  REAL(8), PARAMETER :: twopi = 2*acos(-1d0)

  !-------------------------------------------------
  ! fermion field parameters
  !-------------------------------------------------

  ! number of fermion flavors
  !INTEGER nflv ! already defined in dqmc_core

  ! number of copies of the flavor
  INTEGER, ALLOCATABLE :: ncopy(:)

  ! lattice size
  INTEGER La, Lb, Lc

  ! number of orbitals
  INTEGER norb

  ! unit cell in real space
  REAL(8) a0r(3), b0r(3), c0r(3)

  ! unit cell in momentum space
  REAL(8) a0k(3), b0k(3), c0k(3)

  ! position of each orbital fermion
  REAL(8), ALLOCATABLE :: rorb(:,:)  ! (3,norb)

  ! hopping range
  INTEGER cuta, cutb, cutc

  ! hopping parameters (including onsite energy)
  COMPLEX(8), ALLOCATABLE :: hop(:,:,:,:,:,:) ! (cuta,cutb,cutc,norb,norb,nflv)

  ! inverse of temperature
  REAL(8) beta

  ! length of each time slice
  REAL(8) dtau

  ! periodic boundary condition or not along each direction
  LOGICAL pbca, pbcb, pbcc

  ! twisted boundary condition
  REAL(8) twista, twistb, twistc

  ! kinetic Hamiltonian
  COMPLEX(8), ALLOCATABLE :: kmat(:,:,:)

  !--------------------------------------------------
  ! boson field parameters
  !   gamma_ising * exp( lam_ising * (fmat-f0) )
  ! ->gamma_ising * exp( lam_ising * fmat )
  ! or
  !   exp( - field * f0 ) * exp( field * fmat )
  !--------------------------------------------------

  ! the type of the boson field
  ! dimension (nfield)
  ! type_field =  1 : HS1
  !               2 : HS2
  !               3 : HSgeneral
  !              -1 : continuous HS
  !              -2 : local phonon
  !          others : undefined yet
  !           >=100 : user-defined ising
  !          <=-100 : user-defined phi
  INTEGER, ALLOCATABLE :: type_field(:)

  ! interaction or other parameters to characterize the boson field
  !   H = g(1)/2 * ( c'*fmat*c - g(2) )^2
  ! for phonon
  !   H = eta * ( c'*fmat*c - g(2) )
  !   with g(3)=debye and g(1)=eta^2/M/debye^2
  ! dimension (3,nfield)
  ! in all default models, at most 3 parameters are enough to describe the boson field
  REAL(8), ALLOCATABLE :: g_field(:,:)

  ! the number of allowed values of the boson field
  ! if isingmax=0, it corresponds to a continuous field
  ! dimension (nfield)
  INTEGER, ALLOCATABLE :: isingmax(:)

  ! local update of ising-type field:
  !   flip of the ising field to a new value
  ! for continuous field, it's useless
  ! dimension (maxval(isingmax)-1,maxval(isingmax))
  INTEGER, ALLOCATABLE :: isingflip(:,:)

  ! local update of phi-type field:
  !   phi -> phi + dphi/dphi_global
  ! for ising-type field, it's useless
  ! dimension (nfield)
  COMPLEX(8), ALLOCATABLE :: dphi(:), dphi_global(:)

  ! indicate the method to do global udpate
  INTEGER, ALLOCATABLE :: global_method(:)

  ! form factor of a given boson field
  ! dimension (maxval(ndim_field),maxval(ndim_field),nfield,nflv)
  COMPLEX(8), ALLOCATABLE :: fmat(:,:,:,:)

  ! dimension (maxval(ndim_field),maxval(ndim_field),nfield,nflv)
  COMPLEX(8), ALLOCATABLE :: fmat_U(:,:,:,:)

  ! dimension (maxval(ndim_field),nfield,nflv)
  REAL(8), ALLOCATABLE :: fmat_expE(:,:,:)

  ! lambda of the ising-type field
  ! dimension (maxval(isingmax),nfield)
  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:)

  ! gamma of the ising-type field
  ! dimension (maxval(isingmax),nfield)
  COMPLEX(8), ALLOCATABLE :: gamma_ising(:,:)

  ! exp(lam*F) and exp(-lam*F)
  ! only used for ising-type field
  ! dimension (maxval(ndim_field),maxval(ndim_field),maxval(isingmax),nfield,nflv)
  COMPLEX(8), ALLOCATABLE :: expflam_ising(:,:,:,:,:), inv_expflam_ising(:,:,:,:,:)

  ! exp((lam'-lam)*F)-1
  ! only used for ising-type field
  ! dimension (maxval(ndim_field),maxval(ndim_field),maxval(isingmax),maxval(isingmax),nfield,nflv)
  COMPLEX(8), ALLOCATABLE :: diff_ef_ising(:,:,:,:,:,:)

  !-------------------------------------
  ! measurement parameters
  ! default measurements are classified into 3 classies:
  !   single particle Green's function <c'(i,t)c(j,0)>
  !   PH channel two-particle Green's function <PH'(i,t)PH(j,0)>
  !   PP channel two-particle Green's function <PP'(i,t)PP(j,0)>
  !-------------------------------------

  !LOGICAL meas_k, meas_r, meas_rr, meas_tau

  INTEGER nk_meas, nr_meas, nrr_meas, ntau_meas

  INTEGER, ALLOCATABLE :: k_array(:,:), r_array(:,:), rr_array(:,:,:)

  COMPLEX(8), ALLOCATABLE :: expikr_array(:,:,:,:)

  INTEGER n_ph_meas, n_pp_meas

  LOGICAL, ALLOCATABLE :: hartree_ph_meas(:), fork_ph_meas(:)
  
  LOGICAL, ALLOCATABLE :: fork14_pp_meas(:), fork13_pp_meas(:)
  
  CHARACTER(20), ALLOCATABLE :: name_ph_meas(:), name_pp_meas(:)

  INTEGER, ALLOCATABLE :: ndim_ph_meas(:), ndim_pp_meas(:)

  INTEGER, ALLOCATABLE :: flv_ph_meas(:,:)

  INTEGER, ALLOCATABLE :: nb_ph_meas(:,:,:,:,:)  !(La,Lb,Lc,maxval(ndim_ph_meas),n_ph_meas)

  INTEGER, ALLOCATABLE :: flv_pp_meas(:,:)

  INTEGER, ALLOCATABLE :: nb_pp_meas(:,:,:,:,:)  !(La,Lb,Lc,maxval(ndim_pp_meas),n_pp_meas)

  COMPLEX(8), ALLOCATABLE :: fmat_ph_meas(:,:,:), fmat_pp_meas(:,:,:)

  INTEGER ncross_ph_meas, ncross_pp_meas

  CHARACTER(8), ALLOCATABLE :: name_cross_ph_meas(:), name_cross_pp_meas(:)

  ! dimension (2,ncross_ph_meas) and (2,ncross_pp_meas)
  INTEGER, ALLOCATABLE :: cross_ph_meas(:,:), cross_pp_meas(:,:)

  LOGICAL FAtech

  LOGICAL do_measure_external

  !---------------------
  ! other controls
  !---------------------

  LOGICAL do_tmpout_external

  LOGICAL do_postprocess_external


CONTAINS

  ! group (a,b,c,orb) to a unique index
  INTEGER FUNCTION label(a,b,c,orb)
    IMPLICIT NONE
    INTEGER a,b,c,orb
    label=(a-1)*Lb*Lc*norb+(b-1)*Lc*norb+(c-1)*norb+orb
  END FUNCTION

  ! set kinetic Hamiltonian matrix
  SUBROUTINE set_kmat()
    IMPLICIT NONE
    REAL(8), PARAMETER :: twopi=acos(-1d0)*2
    INTEGER a,b,c,orb,i,flv
    INTEGER a2,b2,c2,orb2,i2
    INTEGER da,db,dc
    COMPLEX(8) boundary

    IF(.not.allocated(kmat)) ALLOCATE(kmat(La*Lb*Lc*norb,La*Lb*Lc*norb,nflv))

    kmat=0d0

    DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
      DO da=-cuta,cuta; DO db=-cutb,cutb; DO dc=-cutc,cutc; DO orb2=1,norb

        a2=a+da; b2=b+db; c2=c+dc

        IF(a2>=1.and.a2<=La.and.b2>=1.and.b2<=Lb.and.c2>=1.and.c2<=Lc)THEN ! without crossing the boundary

          i2=label(a2,b2,c2,orb2)

          DO flv=1,nflv
            kmat(i,i2,flv)=hop(da,db,dc,orb,orb2,flv)
            !kmat(i2,i,flv)=conjg(kmat(i,i2,flv))
          END DO

        ELSE ! crossing the boundary

          IF(a2<1.and.(.not.pbca))CYCLE
          IF(b2<1.and.(.not.pbcb))CYCLE
          IF(c2<1.and.(.not.pbcc))CYCLE

          IF(a2>La.and.(.not.pbca))CYCLE
          IF(b2>Lb.and.(.not.pbcb))CYCLE
          if(c2>Lc.and.(.not.pbcc))CYCLE

          boundary=1d0

          IF(a2<1)THEN
            boundary=boundary*exp(dcmplx(0d0,-twista*twopi))
            a2=a2+La
          END IF

          IF(a2>La)THEN
            boundary=boundary*exp(dcmplx(0d0,twista*twopi))
            a2=a2-La
          END IF

          IF(b2<1)THEN
            boundary=boundary*exp(dcmplx(0d0,-twistb*twopi))
            b2=b2+Lb
          END IF

          IF(b2>Lb)THEN
            boundary=boundary*exp(dcmplx(0d0,twistb*twopi))
            b2=b2-Lb
          END IF

          IF(c2<1)THEN
            boundary=boundary*exp(dcmplx(0d0,-twistc*twopi))
            c2=c2+Lc
          END IF

          IF(c2>Lc)THEN
            boundary=boundary*exp(dcmplx(0d0,twistc*twopi))
            c2=c2-Lc
          END IF

          IF(a2==a.and.b2==b.and.c2==c)CYCLE

          i2=label(a2,b2,c2,orb2)

          DO flv=1,nflv
            kmat(i,i2,flv)=hop(da,db,dc,orb,orb2,flv)*boundary
            !kmat(i2,i,flv)=conjg(kmat(i,i2))
          END DO

        END IF

      END DO; END DO; END DO; END DO
    END DO; END DO; END DO; END DO

  END SUBROUTINE set_kmat

  ! set exp(K)=exp(-dtau*H0)
  ! NOTICE: before calling this subroutine, kmat must be generated
  !         after this subroutine, kmat is changed and must be set again
  SUBROUTINE set_expk()
    IMPLICIT NONE
    INTEGER i,j,flv
    REAL(8) eval(nsite)
    IF(.not.allocated(expk))ALLOCATE(expk(nsite,nsite,nflv))
    IF(.not.allocated(inv_expk))ALLOCATE(inv_expk(nsite,nsite,nflv))
    IF(.not.allocated(expk_half))ALLOCATE(expk_half(nsite,nsite,nflv))
    IF(.not.allocated(inv_expk_half))ALLOCATE(inv_expk_half(nsite,nsite,nflv))
    DO flv=1,nflv
      CALL eigen(nsite,kmat(:,:,flv),eval)
      DO i=1,nsite
        DO j=1,nsite
          expk(i,j,flv)=sum(kmat(i,:,flv)*exp(-dtau*eval(:))*conjg(kmat(j,:,flv)))
          inv_expk(i,j,flv)=sum(kmat(i,:,flv)*exp(dtau*eval(:))*conjg(kmat(j,:,flv)))
          expk_half(i,j,flv)=sum(kmat(i,:,flv)*exp(-dtau*eval(:)/2)*conjg(kmat(j,:,flv)))
          inv_expk_half(i,j,flv)=sum(kmat(i,:,flv)*exp(dtau*eval(:)/2)*conjg(kmat(j,:,flv)))
        END DO
      END DO
    END DO
  END SUBROUTINE

  ! set slater for T=0 dqmc using kmat
  ! NOTICE: before calling this subroutine, kmat must be generated and expk_half must be set
  !         after this subroutine, kmat is changed and must be set again
  SUBROUTINE set_slater()
    IMPLICIT NONE
    INTEGER flv
    REAL(8) eval(nsite)
    IF(.not.proj)THEN
      PRINT*,'slater is not needed for T>0'
      RETURN
    END IF
    IF(abs(dot_product(expk_half(1,:,1),inv_expk_half(:,1,1))-1d0)>1d-6)THEN
      PRINT*,'expk_half has not been correctly set. It is required for 2nd-order Trotter.'
      CALL exit(0)
    END IF
    IF(.not.allocated(slater))ALLOCATE(slater(nsite,maxval(nelec),nflv))
    IF(.not.allocated(slater_Q))ALLOCATE(slater_Q(nsite,maxval(nelec),nflv))
    IF(.not.allocated(slater_D))ALLOCATE(slater_D(maxval(nelec),nflv))
    IF(.not.allocated(slater_R))ALLOCATE(slater_R(maxval(nelec),maxval(nelec),nflv))
    DO flv=1,nflv
      CALL eigen(nsite,kmat(:,:,flv),eval)
      slater(:,1:nelec(flv),flv)=kmat(:,1:nelec(flv),flv)
      slater(:,1:nelec(flv),flv)=matmul(expk_half(:,:,flv),slater(:,1:nelec(flv),flv))
      slater_Q(:,1:nelec(flv),flv)=slater(:,1:nelec(flv),flv)
      CALL qdr(nsite,nelec(flv),slater_Q(:,1:nelec(flv),flv),slater_R(1:nelec(flv),1:nelec(flv),flv), &
        & slater_D(1:nelec(flv),flv))
    END DO
  END SUBROUTINE

  !
  SUBROUTINE init_field_random(ifield)
    IMPLICIT NONE
    INTEGER ifield,site,time
    field(:,:,ifield)=0d0
    IF(.not.mask_field(ifield))RETURN
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      DO time=1,ntime
        IF(isingmax(ifield)==0)THEN
          field(site,time,ifield)=drand_sym()*dphi(ifield)
        ELSE
          field(site,time,ifield)=irand(isingmax(ifield))+1
        END IF
      END DO
    END DO
  END SUBROUTINE

  ! set isingflip
  SUBROUTINE set_isingflip()
    IMPLICIT NONE
    INTEGER n,oldising,i
    n=maxval(isingmax)
    IF(.not.allocated(isingflip)) ALLOCATE(isingflip(n-1,n))
    DO oldising=1,n
      DO i=1,n
        IF(i<oldising)THEN
          isingflip(i,oldising)=i
        ELSEIF(i>oldising)THEN
          isingflip(i-1,oldising)=i
        END IF
      END DO
    END DO
  END SUBROUTINE

  !
  SUBROUTINE set_expf()
    IMPLICIT NONE
    INTEGER ifield,d,z,a,b,s,s2,flv

    d=maxval(ndim_field)
    z=maxval(isingmax)

    ! allocate fmat_U and fmat_expE
    IF(.not.allocated(fmat_U)) ALLOCATE(fmat_U(d,d,nfield,nflv))
    IF(.not.allocated(fmat_expE)) ALLOCATE(fmat_expE(d,nfield,nflv))

    IF(z>0)THEN
      ! allocate (inv_)expflam_ising and diff_ef_ising used for ising-type fields
      IF(.not.allocated(expflam_ising)) ALLOCATE(expflam_ising(d,d,z,nfield,nflv))
      IF(.not.allocated(inv_expflam_ising)) ALLOCATE(inv_expflam_ising(d,d,z,nfield,nflv))
      IF(.not.allocated(diff_ef_ising)) ALLOCATE(diff_ef_ising(d,d,z,z,nfield,nflv))
    END IF

    DO ifield=1,nfield

      d=ndim_field(ifield)
      z=isingmax(ifield)

      DO flv=1,nflv

        ! set fmat_U and fmat_expE
        fmat_U(1:d,1:d,ifield,flv)=fmat(1:d,1:d,ifield,flv)
        CALL eigen(d,fmat_U(1:d,1:d,ifield,flv),fmat_expE(1:d,ifield,flv))
        fmat_expE(1:d,ifield,flv)=exp(fmat_expE(1:d,ifield,flv))

        ! set expf(lam*fmat) and expf(-lam*fmat)
        DO s=1,z
          DO a=1,d
            DO b=1,d
              expflam_ising(a,b,s,ifield,flv)=sum(fmat_U(a,1:d,ifield,flv) &
                & *fmat_expE(1:d,ifield,flv)**lam_ising(s,ifield)*conjg(fmat_U(b,1:d,ifield,flv)))
              inv_expflam_ising(a,b,s,ifield,flv)=sum(fmat_U(a,1:d,ifield,flv) &
                & /fmat_expE(1:d,ifield,flv)**lam_ising(s,ifield)*conjg(fmat_U(b,1:d,ifield,flv)))
            END DO
          END DO
        END DO

        ! set expf((lam'-lam)*fmat)-1
        DO s=1,z
          DO s2=1,z
            diff_ef_ising(1:d,1:d,s,s2,ifield,flv)=matmul(expflam_ising(1:d,1:d,s,ifield,flv), &
              & inv_expflam_ising(1:d,1:d,s2,ifield,flv))
            DO a=1,d
              diff_ef_ising(a,a,s,s2,ifield,flv)=diff_ef_ising(a,a,s,s2,ifield,flv)-1d0
            END DO
          END DO
        END DO

      END DO

    END DO

  END SUBROUTINE

END MODULE

!
SUBROUTINE generate_newfield_local(newfield,site,time,delta,ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER site,time,ifield,d,a,b,flv,z
  COMPLEX(8) newfield,oldfield,delta(ndim_field(ifield),ndim_field(ifield),nflv)
  d=ndim_field(ifield)
  z=isingmax(ifield)
  oldfield=field(site,time,ifield)
  IF(z==0)THEN
    newfield=oldfield+drand_sym()*dphi(ifield)
    DO flv=1,nflv
      DO a=1,d
        DO b=1,d
          delta(a,b,flv)=sum(fmat_U(a,1:d,ifield,flv)*fmat_expE(1:d,ifield,flv) &
            & **(newfield-oldfield)*conjg(fmat_U(b,1:d,ifield,flv)))
        END DO
        delta(a,a,flv)=delta(a,a,flv)-1d0
      END DO
    END DO
  ELSE
    newfield=isingflip(irand(z-1)+1,nint(real(oldfield)))
    DO flv=1,nflv
      delta(1:d,1:d,flv)=diff_ef_ising(1:d,1:d,nint(real(newfield)),nint(real(oldfield)),ifield,flv)
    END DO
  END IF
END SUBROUTINE

!
SUBROUTINE acceptprob_local(ratio,newfield,site,time,ifield,rtot)
  USE dqmc
  IMPLICIT NONE
  INTEGER site,time,ifield,flv
  COMPLEX(8) ratio(nflv),newfield,rtot,oldfield,gph_x2,gph_p2,f0,diffnew,diffold
  oldfield=field(site,time,ifield)
  IF(type_field(ifield)>=1.and.type_field(ifield)<=3)THEN
    ! ising-type HS field
    ! NOTICE f0 has been absorbed into the definition of gamma_ising
    rtot=gamma_ising(nint(real(newfield)),ifield)/gamma_ising(nint(real(oldfield)),ifield)
    DO flv=1,nflv
      rtot=rtot*ratio(flv)**ncopy(flv)
    END DO
  ELSEIF(type_field(ifield)==-1)THEN
    ! continuous HS field
    gph_x2=0.5d0/g_field(1,ifield)/dtau
    f0=g_field(2,ifield)
    rtot=exp(-gph_x2*(newfield**2-oldfield**2)-f0*(newfield-oldfield))
    DO flv=1,nflv
      rtot=rtot*ratio(flv)**ncopy(flv)
    END DO
  ELSEIF(type_field(ifield)==-2)THEN
    ! local phonon field
    gph_x2=0.5d0/g_field(1,ifield)/dtau
    gph_p2=0.5d0/g_field(1,ifield)/dtau**3/g_field(3,ifield)**2
    f0=g_field(2,ifield)
    rtot=exp(-gph_x2*(newfield**2-oldfield**2)-f0*(newfield-oldfield))
    diffnew=newfield-field(site,mod(time,ntime)+1,ifield)
    diffold=oldfield-field(site,mod(time,ntime)+1,ifield)
    rtot=rtot*exp(-gph_p2*(diffnew**2-diffold**2)) ! kinetic energy on (time,time+1)
    diffnew=newfield-field(site,mod(time-2+ntime,ntime)+1,ifield)
    diffold=oldfield-field(site,mod(time-2+ntime,ntime)+1,ifield)
    rtot=rtot*exp(-gph_p2*(diffnew**2-diffold**2)) ! kinetic energy on (time,time-1)
    DO flv=1,nflv
      rtot=rtot*ratio(flv)**ncopy(flv)
    END DO
  ELSEIF(abs(type_field(ifield))>=100)THEN
    CALL acceptprob_local_external(ratio,newfield,site,time,ifield,rtot)
  ELSE
    IF(id==0) PRINT*,'undefined type_field'
    CALL exit(0)
  END IF
END SUBROUTINE

!
SUBROUTINE generate_newfield_global(ifield)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield,site,time
  COMPLEX(8) shift
  IF(global_method(ifield)==0)THEN
    CALL init_field_random(ifield)
  ELSEIF(global_method(ifield)==100)THEN
    CALL generate_newfield_global_external(ifield)
  ELSE
    IF(isingmax(ifield)==0)THEN
      ! continuous field
      IF(global_method(ifield)==1)THEN
        DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
          DO time=1,ntime
            field(site,time,ifield)=field(site,time,ifield)+drand_sym()*dphi_global(ifield)
          END DO
        END DO
      ELSEIF(global_method(ifield)==2)THEN
        DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
          shift=drand_sym()*dphi_global(ifield)
          field(site,:,ifield)=field(site,:,ifield)+shift
        END DO
      ELSE
        IF(id==0) PRINT*,'undefined global_method for continuous field'
        CALL exit(0)
      END IF
    ELSE
      ! ising-type field
      IF(id==0) PRINT*,'undefined global_method for ising-type field'
      CALL exit(0)
    END IF
  END IF
END SUBROUTINE

SUBROUTINE acceptprob_global(ratio,newfield,ifield,rtot)
  USE dqmc
  IMPLICIT NONE
  INTEGER ifield,site,time,flv
  COMPLEX(8) ratio(nflv),newfield(nsite,ntime),rtot
  COMPLEX(8) gph_x2,gph_p2,f0,newphi,oldphi,diffnew,diffold
  IF(type_field(ifield)>=1.and.type_field(ifield)<=3)THEN
    ! ising-type HS field
    rtot=1d0
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      DO time=1,ntime
        rtot=rtot*gamma_ising(nint(real(newfield(site,time))), ifield) &
          & /gamma_ising(nint(real(field(site,time,ifield))),ifield)
      END DO
    END DO
    DO flv=1,nflv
      rtot=rtot*ratio(flv)**ncopy(flv)
    END DO
  ELSEIF(type_field(ifield)==-1)THEN
    ! continuous HS field
    gph_x2=0.5d0/g_field(1,ifield)/dtau
    f0=g_field(2,ifield)
    rtot=1d0
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      DO time=1,ntime
        newphi=newfield(site,time)
        oldphi=field(site,time,ifield)
        rtot=rtot*exp(-gph_x2*(newphi**2-oldphi**2)-f0*(newphi-oldphi))
      END DO
    END DO
    DO flv=1,nflv
      rtot=rtot*ratio(flv)**ncopy(flv)
    END DO
  ELSEIF(type_field(ifield)==-2)THEN
    ! local phonon field
    gph_x2=0.5d0/g_field(1,ifield)/dtau
    gph_p2=0.5d0/g_field(1,ifield)/dtau**3/g_field(3,ifield)**2
    f0=g_field(2,ifield)
    rtot=1d0
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      DO time=1,ntime
        newphi=newfield(site,time)
        oldphi=field(site,time,ifield)
        rtot=rtot*exp(-gph_x2*(newphi**2-oldphi**2)-f0*(newphi-oldphi))
        diffnew=newphi-newfield(site,mod(time,ntime)+1)
        diffold=oldphi-field(site,mod(time,ntime)+1,ifield)
        rtot=rtot*exp(-gph_p2*(diffnew**2-diffold**2))
      END DO
    END DO
    DO flv=1,nflv
      rtot=rtot*ratio(flv)**ncopy(flv)
    END DO
  ELSEIF(abs(type_field(ifield))>=100)THEN
    CALL acceptprob_global_external(ratio,newfield,ifield,rtot)
  ELSE
    IF(id==0) PRINT*,'undefined type_field'
    CALL exit(0)
  END IF
END SUBROUTINE

!
SUBROUTINE get_expV(site,time,ifield,flv,inv,expV)
  USE dqmc
  IMPLICIT NONE
  INTEGER site,time,ifield,flv,z,d,a,b
  LOGICAL inv
  COMPLEX(8) expV(ndim_field(ifield),ndim_field(ifield)),newfield
  d=ndim_field(ifield)
  z=isingmax(ifield)
  newfield=field(site,time,ifield)
  IF(z==0)THEN
    DO a=1,d
      DO b=1,d
        IF(inv)THEN
          expV(a,b)=sum(fmat_U(a,1:d,ifield,flv)/fmat_expE(1:d,ifield,flv) &
            & **newfield*conjg(fmat_U(b,1:d,ifield,flv)))
        ELSE
          expV(a,b)=sum(fmat_U(a,1:d,ifield,flv)*fmat_expE(1:d,ifield,flv) &
            & **newfield*conjg(fmat_U(b,1:d,ifield,flv)))
        END IF
      END DO
    END DO
  ELSE
    IF(inv)THEN
      expV=inv_expflam_ising(1:d,1:d,nint(real(newfield)),ifield,flv)
    ELSE
      expV=expflam_ising(1:d,1:d,nint(real(newfield)),ifield,flv)
    END IF
  END IF
END SUBROUTINE

PROGRAM main
  USE dqmc
  IMPLICIT NONE
  CALL dqmc_driver()
END PROGRAM
