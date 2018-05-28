!---------------------------------------------------------------------------
! Tis MODULE is designed for defining a 3D (1-2D are special cases) fermionic 
! lattice coupling to boson fields. The boson fields can be either descrete 
! (called ising) or continuous (called phi).
!
! NOTICE: In this module, uniform fermion-fermion and fermion-boson
! interactions are supposed. For e.g. disorder case, we can simply add
! site-dependence to relevant variables, e.g. lam_field, gamma_field, etc.
!
! Several useful subroutines are provided to help users to do initialization. 
!---------------------------------------------------------------------------
MODULE model_complex

  !-------------------------------------------------
  ! fermion field parameters
  !-------------------------------------------------

  ! number of fermion flavors
  INTEGER nflv

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

  !---------------------------------------------------
  ! Ising-type boson field parameters
  !   gamma(field) * exp( lam(field) * fmat(:,:) )
  ! NOTICE: lam can be viewed as fermion-boson coupling and
  !         gamma can be viewed as free-boson action
  !---------------------------------------------------

  ! number of ising fields
  INTEGER nising

  ! dimension for each ising field
  ! dimension (nising)
  INTEGER, ALLOCATABLE :: ndim_ising(:) 

  ! list of all sites belong to one ising field living on one representative site.
  ! dimension (nsite,maxval(ndim_ising),nfield)
  INTEGER, ALLOCATABLE :: nb_ising(:,:,:)
 
  ! whether an ising field exists?
  ! dimension (nising)
  LOGICAL, ALLOCATABLE :: mask_ising(:)

  ! whether a site lives an ising field
  ! dimension (nsite,nising)
  LOGICAL, ALLOCATABLE :: mask_ising_site(:,:)

  ! the number of allowed values of an ising field
  ! dimension (nising)
  INTEGER, ALLOCATABLE :: isingmax(:)

  ! flip an ising field to a new value
  ! dimension (maxval(isingmax)-1,maxval(isingmax))
  INTEGER, ALLOCATABLE :: isingflip(:,:)

  ! form factor of a given ising field
  ! dimension (maxval(ndim_ising),maxval(ndim_ising),nising)
  COMPLEX(8), ALLOCATABLE :: fmat_ising(:,:,:)

  !
  COMPLEX(8), ALLOCATABLE :: fmat_ising_U(:,:,:)
  
  !
  COMPLEX(8), ALLOCATABLE :: fmat_ising_expE(:,:)

  ! lambda of the ising field
  ! dimension (maxval(isingmax),nising)
  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:)

  ! gamma of the iisng field
  ! dimension (maxval(isingmax),nising)
  COMPLEX(8), ALLOCATABLE :: gamma_ising(:,:)

  ! exp(lam*F) and exp(-lam*F)
  ! dimension (maxval(ndim_ising),maxval(ndim_ising),isingmax,nising)
  COMPLEX(8), ALLOCATABLE :: expflam_ising(:,:,:,:), inv_expflam_ising(:,:,:,:)
  
  ! exp((lam'-lam)*F)-1
  ! dimension (maxval(ndim_ising),maxval(ndim_ising),isingmax,isingmax,nising)
  COMPLEX(8), ALLOCATABLE :: diff_ef_ising(:,:,:,:,:)


  !--------------------------------------------------
  ! continuous boson field parameters
  !   exp( field * fmat(:,:) )
  ! the free boson action is not defined here
  !--------------------------------------------------

  ! number of continuous boson fields
  INTEGER nphi
  
  ! dimension for each continuous field
  INTEGER, ALLOCATABLE :: ndim_phi(:)  ! (nphi)

  ! dimension (nsite,maxval(ndim_phi),nphi)
  INTEGER, ALLOCATABLE :: nb_phi(:,:,:) 

  ! dimension (nphi)
  LOGICAL, ALLOCATABLE :: mask_phi(:)

  ! dimension (nsite,nphi)
  LOGICAL, ALLOCATABLE :: mask_phi_site(:,:)

  ! update phi by changing to phi+dphi/dphi_global
  ! dimension (nphi)
  COMPLEX(8), ALLOCATABLE :: dphi(:), dphi_global(:)

  ! dimension (maxval(ndim_phi),maxval(ndim_phi),nphi)
  COMPLEX(8), ALLOCATABLE :: fmat_phi(:,:,:)

  ! dimension (maxval(ndim_phi),maxval(ndim_phi),nphi)
  COMPLEX(8), ALLOCATABLE :: fmat_phi_U(:,:,:)

  ! dimension (maxval(ndim_phi),nphi)
  REAL(8), ALLOCATABLE :: fmat_phi_expE(:,:)

CONTAINS

  ! group (a,b,c,orb) to a unique index
  INTEGER FUNCTION label(a,b,c,orb)
    IMPLICIT NONE
    INTEGER a,b,orb
    label=(a-1)*Lb*Lc*norb+(b-1)*Lc*norb+(c-1)*norb+orb
  END FUNCTION

  ! set kinetic Hamiltonian matrix
  SUBROUTINE set_kmat()
    IMPLICIT NONE
    REAL(8), PARAMETER :: twopi=acos(-1d0)*2
    INTEGER a,b,c,orb,i
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
    USE dqmc_complex, ONLY : nsite,expk,inv_expk,expk_half,inv_expk_half
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
    USE dqmc_complex, ONLY : proj,nsite,nelec,slater,slater_U,slater_D,slater_R,expk_half
    IMPLICIT NONE
    INTEGER flv
    REAL(8) eval(nsite)
    IF(.not.proj)THEN
      PRINT*,'slater is not needed for T>0'
      RETURN
    END IF
    IF(abs(matmul(expk_half(1,:,1),inv_expk_half(:,1,1))-1d0)>1d-6)THEN
      PRINT*,'expk_half has not been correctly set. It is required for 2nd-order Trotter.'
      CALL exit(0)
    END IF
    IF(.not.allocated(slater))ALLOCATE(slater(nsite,nelec,nflv))
    IF(.not.allocated(slater_U))ALLOCATE(slater_U(nsite,nelec,nflv))
    IF(.not.allocated(slater_D))ALLOCATE(slater_D(nelec,nflv))
    IF(.not.allocated(slater_R))ALLOCATE(slater_R(nelec,nelec,nflv))
    DO flv=1,nflv
      CALL eigen(nsite,kmat(:,:,flv),eval)
      slater(:,:,flv)=kmat(:,1:nelec,flv)
      slater(:,:,flv)=matmul(expk_half(:,:,flv),slater(:,:,flv))
      slater_U(:,:,flv)=slater(:,:,flv)
      CALL dqr(nsite,nelec,slater_U(:,:,flv),slater_R(:,:,flv),slater_D(:,flv))
    END DO
  END SUBROUTINE

  ! initialize auxiliary field configuration randomly
  SUBROUTINE init_field_random()
    USE dqmc_complex, ONLY : nsite,ntime,field
    IMPLICIT NONE
    INTEGER ifield,time,site
    IF(.not.allocated(field)) ALLOCATE(field(nsite,ntime,nising+nphi))
    field=0d0
    DO ifield=1,nising; IF(.not.mask_ising(ifield))CYCLE
      DO time=1,ntime
        DO site=1,nsite
          IF(mask_ising_site(site,ifield)) field(site,time,ifield)=irand(isingmax(ifield))+1
        END DO
      END DO
    END DO
    DO ifield=1,nphi; IF(.not.mask_phi(ifield))CYCLE
      DO time=1,ntime
        DO site=1,nsite
          IF(mask_phi_site(site,ifield)) field(site,time,nising+ifield)=drand_sym()*dphi(ifield)
        END DO
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
    INTEGER ifield,n,z,d,a,b,s,s2
    
    n=maxval(ndim_ising)
    z=maxval(isingmax)
    
    IF(.not.allocated(fmat_ising_U)) ALLOCATE(fmat_ising_U(n,n,nising))
    IF(.not.allocated(fmat_ising_expE)) ALLOCATE(fmat_ising_expE(n,nising))
    IF(.not.allocated(expflam)) ALLOCATE(expflam(n,n,z,nising))
    IF(.not.allocated(inv_expflam)) ALLOCATE(inv_expflam(n,n,z,nising))
    IF(.not.allocated(diff_ef)) ALLOCATE(diff_ef(n,n,z,z,nising))
    
    DO ifield=1,nising
      
      d=ndim_ising(ifield)

      ! set fmat_ising_U and fmat_ising_expE
      fmat_ising_U(1:d,1:d,ifield)=fmat_ising(1:d,1:d,ifield)
      CALL eigen(d,fmat_ising_U(1:d,1:d,ifield),fmat_ising_expE(1:d,ifield))
      fmat_ising_expE(1:d,ifield)=exp(fmat_ising_expE(1:d,ifield))
      
      ! set expf(lam*fmat) and expf(-lam*fmat)
      DO s=1,z
        DO a=1,d
          DO b=1,d
            expflam(a,b,s,ifield)=sum(fmat_ising_U(a,1:d,ifield) &
            & *fmat_ising_expE(1:d,ifield)*conjg(fmat_ising_U(b,1:d,ifield)))    
            inv_expflam(a,b,s,ifield)=sum(fmat_ising_U(a,1:d,ifield) &
            & /fmat_ising_expE(1:d,ifield)*conjg(fmat_ising_U(b,1:d,ifield)))    
          END DO
        END DO
      END DO
      
      ! set expf((lam'-lam)*fmat)-1
      DO s=1,z
        DO s2=1,z
          diff_ef(1:d,1:d,s,s2,ifield)=matmul(expflam(1:d,1:d,s,ifield), &
          & inv_expvlam(1:d,1:d,s2,ifield))
          DO a=1,d
            diff_ef(a,a,s,s2,ifield)=diff_ef(a,a,s,s2,ifield)-1d0
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE
  
  !
  SUBROUTINE set_expf_phi()
    IMPLICIT NONE
    INTEGER ifield,n,d
    
    n=maxval(ndim_phi)
    
    IF(.not.allocated(fmat_ising_U)) ALLOCATE(fmat_ising_U(n,n,nphi))
    IF(.not.allocated(fmat_ising_expE)) ALLOCATE(fmat_ising_expE(n,nphi))
    
    DO ifield=1,nphi
      
      d=ndim_phi(ifield)

      ! set fmat_ising_U and fmat_ising_expE
      fmat_phi_U(1:d,1:d,ifield)=fmat_phi(1:d,1:d,ifield)
      CALL eigen(d,fmat_phi_U(1:d,1:d,ifield),fmat_phi_expE(1:d,ifield))
      fmat_phi_expE(1:d,ifield)=exp(fmat_phi_expE(1:d,ifield))
      
    END DO

  END SUBROUTINE
 
END MODULE model_complex


