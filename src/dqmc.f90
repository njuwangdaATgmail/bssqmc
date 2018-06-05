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

  !-------------------------------------------------
  ! fermion field parameters
  !-------------------------------------------------

  ! number of fermion flavors
  !INTEGER nflv ! already defined in dqmc_core

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
  INTEGER, ALLOCATABLE :: type_field(:)
  
  ! the number of allowed values of the boson field
  ! isingmax=0 corresponds to a continuous field
  ! dimension (nfield)
  INTEGER, ALLOCATABLE :: isingmax(:)

  ! flip an ising field to a new value
  ! dimension (maxval(isingmax)-1,maxval(isingmax))
  INTEGER, ALLOCATABLE :: isingflip(:,:)

  ! update phi by changing to phi+dphi/dphi_global
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

  ! lambda of the boson field
  ! dimension (maxval(isingmax),nfield)
  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:)

  ! gamma of the boson field
  ! dimension (maxval(isingmax),nfield)
  COMPLEX(8), ALLOCATABLE :: gamma_ising(:,:)

  ! exp(lam*F) and exp(-lam*F)
  ! dimension (maxval(ndim_field),maxval(ndim_field),maxval(isingmax),nfield,nflv)
  COMPLEX(8), ALLOCATABLE :: expflam_ising(:,:,:,:,:), inv_expflam_ising(:,:,:,:,:)
  
  ! exp((lam'-lam)*F)-1
  ! dimension (maxval(ndim_field),maxval(ndim_field),maxval(isingmax),maxval(isingmax),nfield,nflv)
  COMPLEX(8), ALLOCATABLE :: diff_ef_ising(:,:,:,:,:,:)

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
    IF(.not.allocated(slater))ALLOCATE(slater(nsite,nelec,nflv))
    IF(.not.allocated(slater_Q))ALLOCATE(slater_Q(nsite,nelec,nflv))
    IF(.not.allocated(slater_D))ALLOCATE(slater_D(nelec,nflv))
    IF(.not.allocated(slater_R))ALLOCATE(slater_R(nelec,nelec,nflv))
    DO flv=1,nflv
      CALL eigen(nsite,kmat(:,:,flv),eval)
      slater(:,:,flv)=kmat(:,1:nelec,flv)
      slater(:,:,flv)=matmul(expk_half(:,:,flv),slater(:,:,flv))
      slater_Q(:,:,flv)=slater(:,:,flv)
      CALL qdr(nsite,nelec,slater_Q(:,:,flv),slater_R(:,:,flv),slater_D(:,flv))
    END DO
  END SUBROUTINE

  !
  SUBROUTINE init_ising_random(ifield,ifield_tot)
    IMPLICIT NONE
    INTEGER ifield,ifield_tot,time,site
    field(:,:,ifield_tot)=0
    IF(.not.mask_ising(ifield))RETURN
    DO time=1,ntime
      DO site=1,nsite
        IF(mask_ising_site(site,ifield)) field(site,time,ifield_tot)=irand(isingmax(ifield))+1
      END DO
    END DO
  END SUBROUTINE

  !
  SUBROUTINE init_phi_random(ifield,ifield_tot)
    IMPLICIT NONE
    INTEGER ifield,ifield_tot,time,site
    field(:,:,ifield_tot)=0d0
    IF(.not.mask_phi(ifield))RETURN
    DO time=1,ntime
      DO site=1,nsite
        IF(mask_phi_site(site,ifield)) field(site,time,ifield_tot)=drand_sym()*dphi(ifield)
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
  SUBROUTINE set_expf_ising()
    IMPLICIT NONE
    INTEGER ifield,n,z,d,a,b,s,s2,flv
    
    n=maxval(ndim_ising)
    z=maxval(isingmax)
    
    IF(.not.allocated(fmat_ising_U)) ALLOCATE(fmat_ising_U(n,n,nising,nflv))
    IF(.not.allocated(fmat_ising_expE)) ALLOCATE(fmat_ising_expE(n,nising,nflv))
    IF(.not.allocated(expflam_ising)) ALLOCATE(expflam_ising(n,n,z,nising,nflv))
    IF(.not.allocated(inv_expflam_ising)) ALLOCATE(inv_expflam_ising(n,n,z,nising,nflv))
    IF(.not.allocated(diff_ef_ising)) ALLOCATE(diff_ef_ising(n,n,z,z,nising,nflv))
    
    DO ifield=1,nising
      
      d=ndim_ising(ifield)

      DO flv=1,nflv

        ! set fmat_ising_U and fmat_ising_expE
        fmat_ising_U(1:d,1:d,ifield,flv)=fmat_ising(1:d,1:d,ifield,flv)
        CALL eigen(d,fmat_ising_U(1:d,1:d,ifield,flv),fmat_ising_expE(1:d,ifield,flv))
        fmat_ising_expE(1:d,ifield,flv)=exp(fmat_ising_expE(1:d,ifield,flv))
        
        ! set expf(lam*fmat) and expf(-lam*fmat)
        DO s=1,z
          DO a=1,d
            DO b=1,d
              expflam_ising(a,b,s,ifield,flv)=sum(fmat_ising_U(a,1:d,ifield,flv) &
                & *fmat_ising_expE(1:d,ifield,flv)**lam_ising(s,ifield)*conjg(fmat_ising_U(b,1:d,ifield,flv)))    
              inv_expflam_ising(a,b,s,ifield,flv)=sum(fmat_ising_U(a,1:d,ifield,flv) &
                & /fmat_ising_expE(1:d,ifield,flv)**lam_ising(s,ifield)*conjg(fmat_ising_U(b,1:d,ifield,flv)))    
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
  
  !
  SUBROUTINE set_expf_phi()
    IMPLICIT NONE
    INTEGER ifield,n,d,flv
    
    n=maxval(ndim_phi)
    
    IF(.not.allocated(fmat_phi_U)) ALLOCATE(fmat_phi_U(n,n,nphi,nflv))
    IF(.not.allocated(fmat_phi_expE)) ALLOCATE(fmat_phi_expE(n,nphi,nflv))
    
    DO ifield=1,nphi
      
      d=ndim_phi(ifield)

      DO flv=1,nflv

        ! set fmat_ising_U and fmat_ising_expE
        fmat_phi_U(1:d,1:d,ifield,flv)=fmat_phi(1:d,1:d,ifield,flv)
        CALL eigen(d,fmat_phi_U(1:d,1:d,ifield,flv),fmat_phi_expE(1:d,ifield,flv))
        fmat_phi_expE(1:d,ifield,flv)=exp(fmat_phi_expE(1:d,ifield,flv))

      END DO
      
    END DO

  END SUBROUTINE

  !
  SUBROUTINE generate_newising_local(newising,oldising,delta,ifield)
    IMPLICIT NONE
    INTEGER newising,oldising,ifield,d,flv
    COMPLEX(8) delta(ndim_ising(ifield),ndim_ising(ifield),nflv)
    d=ndim_ising(ifield)
    newising=isingflip(irand(isingmax(ifield)-1)+1,oldising)
    DO flv=1,nflv
      delta(1:d,1:d,flv)=diff_ef_ising(1:d,1:d,newising,oldising,ifield,flv)
    END DO
  END SUBROUTINE

  !
  SUBROUTINE generate_newphi_local(newphi,oldphi,delta,ifield)
    IMPLICIT NONE
    INTEGER ifield,d,a,b,flv
    COMPLEX(8) newphi,oldphi,delta(ndim_phi(ifield),ndim_phi(ifield),nflv)
    d=ndim_phi(ifield)
    newphi=oldphi+drand_sym()*dphi(ifield)
    DO flv=1,nflv
      DO a=1,d
        DO b=1,d
          delta(a,b,flv)=sum(fmat_phi_U(a,1:d,ifield,flv)*fmat_phi_expE(1:d,ifield,flv) &
            & **(newphi-oldphi)*conjg(fmat_phi_U(b,1:d,ifield,flv)))
        END DO
        delta(a,a,flv)=delta(a,a,flv)-1d0
      END DO
    END DO
  END SUBROUTINE

  !
  SUBROUTINE get_expV_ising(newising,ifield,flv,inv,expV)
    IMPLICIT NONE
    INTEGER newising,ifield,flv,d
    LOGICAL inv
    COMPLEX(8) expV(ndim_ising(ifield),ndim_ising(ifield))
    d=ndim_ising(ifield)
    IF(inv)THEN
      expV=inv_expflam_ising(1:d,1:d,newising,ifield,flv)
    ELSE
      expV=expflam_ising(1:d,1:d,newising,ifield,flv)
    END IF
  END SUBROUTINE

  !
  SUBROUTINE get_expV_phi(newphi,ifield,flv,inv,expV)
    IMPLICIT NONE
    INTEGER ifield,flv,a,b,d
    LOGICAL inv
    COMPLEX(8) newphi,expV(ndim_phi(ifield),ndim_phi(ifield))
    d=ndim_phi(ifield)
    DO a=1,d
      DO b=1,d
        IF(inv)THEN
          expV(a,b)=sum(fmat_phi_U(a,1:d,ifield,flv)/fmat_phi_expE(1:d,ifield,flv) &
            & **newphi*conjg(fmat_phi_U(b,1:d,ifield,flv)))
        ELSE
          expV(a,b)=sum(fmat_phi_U(a,1:d,ifield,flv)*fmat_phi_expE(1:d,ifield,flv) &
            & **newphi*conjg(fmat_phi_U(b,1:d,ifield,flv)))
        END IF
      END DO
    END DO
  END SUBROUTINE

END MODULE dqmc

PROGRAM main
  USE dqmc
  IMPLICIT NONE
  CALL dqmc_driver()
END PROGRAM
