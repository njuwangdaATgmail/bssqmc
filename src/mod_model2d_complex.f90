!---------------------------------------------------------------------------
! In this module, uniform fermion-fermion and fermion-boson interactions
! are supposed. For disorder case, we can simply add site-dependence to 
! relevant variables, e.g. lam_field, gamma_field, etc.
!---------------------------------------------------------------------------
MODULE model2d_complex

  !-------------------------------------------------
  ! fermion field parameters
  !-------------------------------------------------

  INTEGER nflv

  INTEGER La, Lb, norb

  REAL(8) a0r(2), b0r(2)

  REAL(8) a0k(2), b0k(2)

  REAL(8), ALLOCATABLE :: rorb(:,:)  ! (2,norb)

  INTEGER cuta, cutb

  COMPLEX(8), ALLOCATABLE :: hop(:,:,:,:,:) ! (cuta,cutb,norb,norb,nflv)

  REAL(8) beta, dtau
  
  LOGICAL pbca, pbcb

  REAL(8) twista, twistb

  COMPLEX(8), ALLOCATABLE :: kmat(:,:,:)

  
  !---------------------------------------------------
  ! Ising-type boson field parameters
  !---------------------------------------------------
  
  INTEGER nising
  
  INTEGER, ALLOCATABLE :: isingmax(:)    

  INTEGER, ALLOCATABLE :: isingflip(:,:) ! (maxval(isingmax)-1,maxval(isingmax))
  
  COMPLEX(8), ALLOCATABLE :: fmat_ising(:,:,:)
  
  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:)

  COMPLEX(8), ALLOCATABLE :: gamma_ising(:,:)
  
  COMPLEX(8), ALLOCATABLE :: expflam_ising(:,:,:,:), inv_expflam_ising(:,:,:,:)
  
  COMPLEX(8), ALLOCATABLE :: diff_ef_ising(:,:,:,:,:)


  !--------------------------------------------------
  ! continuous boson field parameters
  !--------------------------------------------------

  INTEGER nphi

  COMPLEX(8), ALLOCATABLE :: fmat_phi(:,:,:)
  
  COMPLEX(8), ALLOCATABLE :: dphi(:), dphi_global(:)

  COMPLEX(8), ALLOCATABLE :: expf_U_phi(:,:,:), expf_Udag_phi(:,:,:)

  REAL(8), ALLOCATABLE :: expf_E_phi(:,:)

CONTAINS

  INTEGER FUNCTION label(a,b,orb)
    IMPLICIT NONE
    INTEGER a,b,orb
    label=(a-1)*Lb*norb+(b-1)*norb+orb
  END FUNCTION

  SUBROUTINE set_kmat()
    IMPLICIT NONE
    REAL(8), PARAMETER :: twopi=acos(-1d0)*2
    INTEGER a,b,orb,i
    INTEGER a2,b2,orb2,i2
    INTEGER da,db
    COMPLEX(8) boundary

    kmat=0d0

    DO a=1,La; DO b=1,Lb; DO orb=1,norb; i=label(a,b,orb)
      DO da=-cuta,cuta; DO db=-cutb,cutb; DO orb2=1,norb
        
        a2=a+da; b2=b+db

        IF(a2>=1.and.a2<=La.and.b2>=1.and.b2<=Lb)THEN ! without crossing the boundary

          i2=label(a2,b2,orb2)
          kmat(i,i2)=hop(da,db,orb,orb2)
          !kmat(i2,i)=conjg(kmat(i,i2))

        ELSE ! crossing the boundary
          
          IF(a2<1.and.(.not.pbca))CYCLE
          IF(b2<1.and.(.not.pbcb))CYCLE

          IF(a2>La.and.(.not.pbca))CYCLE
          IF(b2>Lb.and.(.not.pbcb))CYCLE

          boundary=1d0
          
          IF(a2<1)THEN
            boundary=exp(dcmplx(0d0,-twista*twopi))
            a2=a2+La
          END IF

          IF(a2>La)THEN
            boundary=exp(dcmplx(0d0,twista*twopi))
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

          IF(a2==a.and.b2==b)CYCLE

          i2=label(a2,b2,orb2)
          kmat(i,i2)=hop(da,db,orb,orb2)*boundary
          !kmat(i2,i)=conjg(kmat(i,i2))

        END IF
 
      END DO; END DO; END DO
    END DO; END DO; END DO

  END SUBROUTINE set_kmat

  !> initialize ising configuration randomly
  SUBROUTINE init_ising_random()
    IMPLICIT NONE
    INTEGER ifield,time,site
    ising=0
    DO ifield=1,nising; IF(.not.mask_ising(ifield))CYCLE
      DO time=1,ntime
        DO site=1,nsite
          IF(mask_form(site,form_ising(ifield))) ising(site,time,ifield)=irand(isingmax(ifield))+1
        END DO
      END DO
    END DO
  END SUBROUTINE

  !> initialize phi configuration randomly
  SUBROUTINE init_phi_random()
    IMPLICIT NONE
    INTEGER ifield,time,site
    phi=0d0
    DO ifield=1,nphi; IF(.not.mask_phi(ifield))CYCLE
      DO time=1,ntime
        DO site=1,nsite
          IF(mask_form(site,form_phi(ifield))) phi(site,time,ifield)=drand_sym()*dphi(ifield)
        END DO
      END DO
    END DO
  END SUBROUTINE

  !> set isingflip
  SUBROUTINE set_isingflip()
    IMPLICIT NONE
    INTEGER n,oldising,i
    n=maxval(isingmax)
    IF(.not.allocated(isingflip))ALLOCATE(isingflip(n-1,n))
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

  !> set expf_E, expf_U, expf_Udag, ie U*exp(E)*Udag=exp(fmat) from fmat
  SUBROUTINE set_expf()
    IMPLICIT NONE
    INTEGER iform,ndim
    IF(.not.allocated(expf_E))ALLOCATE(expf_E(maxval(ndim_form),nform))
    IF(.not.allocated(expf_U))ALLOCATE(expf_U(maxval(ndim_form),maxval(ndim_form),nform))
    IF(.not.allocated(expf_Udag))ALLOCATE(expf_Udag(maxval(ndim_form),maxval(ndim_form),nform))
    DO iform=1,nform
      ndim=ndim_form(iform)
      IF(ndim==1)THEN
        expf_E(1,iform)=exp(1d0)
        expf_U(1,1,iform)=1d0
        expf_Udag(1,1,iform)=1d0
      ELSE
        expf_U(1:ndim,1:ndim,iform)=fmat(1:ndim,1:ndim,iform)
        CALL eigen(ndim,expf_U(1:ndim,1:ndim,iform),expf_E(1:ndim,iform))
        expf_E(1:ndim,iform)=exp(expf_E(1:ndim,iform))
        expf_Udag(1:ndim,1:ndim,iform)=conjg(transpose(expf_U(1:ndim,1:ndim,iform)))
      END IF
    END DO
  END SUBROUTINE

  !> set exp(lam*fmat) and exp(-lam*fmat)
  SUBROUTINE set_expflam()
    IMPLICIT NONE
    INTEGER ifield,iform,ndim,site,s,a,b
    IF(.not.allocated(expflam))THEN
      ALLOCATE(expflam(maxval(ndim_form),maxval(ndim_form),maxval(isingmax),nsite,nising))
    END IF
    IF(.not.allocated(inv_expflam))THEN
      ALLOCATE(inv_expflam(maxval(ndim_form),maxval(ndim_form),maxval(isingmax),nsite,nising))
    END IF
    DO ifield=1,nising
      iform=form_ising(ifield)
      ndim=ndim_form(iform)
      DO site=1,nsite
        DO s=1,isingmax(ifield)
          DO a=1,ndim
            DO b=1,ndim
              expflam(a,b,s,site,ifield)=sum(expf_U(a,1:ndim,iform) &
              & *expf_E(1:ndim,iform)**lam_ising(s,site,ifield)*expf_Udag(1:ndim,b,iform))
              inv_expflam(a,b,s,site,ifield)=sum(expf_U(a,1:ndim,iform) &
              & /expf_E(1:ndim,iform)**lam_ising(s,site,ifield)*expf_Udag(1:ndim,b,iform))
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE
  
  !> set exp(lam'*fmat)/exp(lam*fmat)-1
  SUBROUTINE set_diff_ef()
    IMPLICIT NONE
    INTEGER ifield,ndim,site,a,b,d
    IF(.not.allocated(diff_ef))THEN
      ALLOCATE(diff_ef(maxval(ndim_form),maxval(ndim_form),maxval(isingmax),maxval(isingmax),nsite,nising))
    END IF
    DO ifield=1,nising
      ndim=ndim_form(form_ising(ifield))
      DO site=1,nsite
        DO a=1,isingmax(ifield)
          DO b=1,isingmax(ifield)
            diff_ef(1:ndim,1:ndim,a,b,site,ifield)=matmul(expflam(1:ndim,1:ndim,a,site,ifield),&
            & inv_expflam(1:ndim,1:ndim,b,site,ifield))
            DO d=1,ndim
              diff_ef(d,d,a,b,site,ifield)=diff_ef(d,d,a,b,site,ifield)-1d0
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE
  
  !> set slater for T=0 dqmc
  SUBROUTINE set_slater(kmat)
    IMPLICIT NONE
    COMPLEX(8) kmat(nsite,nsite)
    REAL(8) eval(nsite)
    IF(.not.proj)THEN
      IF(id==0)PRINT*,'Slater is not needed in finite-T DQMC.'
      RETURN
    END IF
    IF(.not.allocated(slater))ALLOCATE(slater(nsite,nelec))
    CALL eigen(nsite,kmat,eval)
    slater(:,:)=kmat(:,1:nelec)
    slater=matmul(expk_half,slater)
  END SUBROUTINE

  !> set exp(K)=exp(-dtau*H0)
  SUBROUTINE set_expk(kmat,dtau)
    IMPLICIT NONE
    COMPLEX(8) kmat(nsite,nsite)
    REAL(8) eval(nsite),dtau
    INTEGER i,j
    IF(.not.allocated(expk))ALLOCATE(expk(nsite,nsite))
    IF(.not.allocated(inv_expk))ALLOCATE(inv_expk(nsite,nsite))
    IF(.not.allocated(expk_half))ALLOCATE(expk_half(nsite,nsite))
    IF(.not.allocated(inv_expk_half))ALLOCATE(inv_expk_half(nsite,nsite))
    CALL eigen(nsite,kmat,eval)
    DO i=1,nsite
      DO j=1,nsite
        expk(i,j)=sum(kmat(i,:)*exp(-dtau*eval(:))*conjg(kmat(j,:)))
        inv_expk(i,j)=sum(kmat(i,:)*exp(dtau*eval(:))*conjg(kmat(j,:)))
        expk_half(i,j)=sum(kmat(i,:)*exp(-dtau*eval(:)/2)*conjg(kmat(j,:)))
        inv_expk_half(i,j)=sum(kmat(i,:)*exp(dtau*eval(:)/2)*conjg(kmat(j,:)))
      END DO
    END DO
  END SUBROUTINE


END MODULE model2d_complex


