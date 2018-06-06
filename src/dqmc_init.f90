! default init()
SUBROUTINE init()
  USE dqmc
  USE hubhol
  IMPLICIT NONE

  INTEGER i,da,db,dc,orb,orb2,nhop
  REAL(8) re,re2,re3
  COMPLEX(8) ga,lam
  
  
  OPEN(10,file='dqmc.in')

  !--------------------------------
  ! read in MC control parameters
  !--------------------------------
  
  READ(10,*) restart
  READ(10,*) proj
  READ(10,*) nflv
  READ(10,*) beta
  READ(10,*) dtau
  
  ntime=ceiling(beta/dtau)
  !IF(mod(ntime,2)==1) ntime=ntime+1
  !IF(ntime<2) ntime=2
  IF(ntime<1) ntime=1
  IF((id==0).and.abs(dtau-beta/ntime)>1d-6) PRINT*, 'dtau is reset to',dtau,' and ntime=',ntime
  dtau=beta/ntime

  READ(10,*) nsp
  IF(proj.and.(nsp>ntime/2-1))THEN
    IF(id==0)PRINT*, 'ERROR: nsp should be less than ntime/2 in T=0 version'
    CALL exit(0)
  END IF

  READ(10,*) nbin
  READ(10,*) nwarmup
  READ(10,*) nmeasure
  READ(10,*) ninterval
  READ(10,*) ntmpout
  READ(10,*) nscratch
  READ(10,*) ngroup
  READ(10,*) randomseed
  IF(randomseed==0)THEN
    CALL init_rng()    ! can different cores have the same time and random number?
  ELSE
    randomseed=randomseed+701703*id
    CALL init_rng(randomseed)
  END IF
  READ(10,*) newMetro
  
  !--------------------------------
  ! read in fermion lattice parameters
  !--------------------------------
  
  READ(10,*) norb
  READ(10,*) La,Lb,Lc; nsite=La*Lb*Lc*norb
  READ(10,*) nelec
  READ(10,*) pbca,pbcb,pbcc
  IF(La<=2)pbca=.false.
  IF(Lb<=2)pbcb=.false.
  IF(Lc<=2)pbcc=.false.
  READ(10,*) twista,twistb,twistc
  READ(10,*) a0r(1:3)
  READ(10,*) b0r(1:3)
  READ(10,*) c0r(1:3)
  ALLOCATE(rorb(3,norb))
  DO orb=1,norb
    READ(10,*) rorb(1:3,orb)
  END DO

  READ(10,*) cuta,cutb,cutc
  ALLOCATE(hop(-cuta:cuta,-cutb:cutb,-cutc:cutc,norb,norb,nflv))
  hop=0d0

  READ(10,*) nhop
  DO i=1,nhop
    READ(10,*) da,db,dc,orb,orb2,re,re2
    hop(da,db,dc,orb,orb2,:)=dcmplx(re,re2)
    hop(-da,-db,-dc,orb2,orb,:)=dcmplx(re,-re2)
  END DO

  ALLOCATE(kmat(nsite,nsite,nflv))

  CALL set_kmat(); CALL set_expk()

  IF(proj)THEN
    re=twista; re2=twistb; re3=twistc
    twista=3d0/17; twistb=3d0/13; twistc=3d0/7
    CALL set_kmat(); CALL set_slater()
    twista=re; twistb=re2; twistc=re3
  END IF

  CALL set_kmat()
  OPEN(20,FILE='kmat.dat'); WRITE(20,'(2f18.10)') kmat; CLOSE(20)

  !--------------------------------
  ! read in boson field parameters
  !--------------------------------
  
  ! In fact, these variables are not necessary. The program should be 
  ! able to get them from the following settings. Right now, they are
  ! read from file directly for simplicity. 
  READ(10,*) n_g, nfield, max_ndim_field, max_isingmax
  
  ALLOCATE( type_field(nfield)                          ); type_field       = 0
  ALLOCATE( g_field(nfield)                             ); g_field          = 0d0
  ALLOCATE( isingmax(nfield)                            ); isingmax         = 0
  ALLOCATE( dphi(nfield)                                ); dphi             = 0d0
  ALLOCATE( dphi_global(nfield)                         ); dphi_global      = 0d0
  ALLOCATE( nglobal_field(nfield)                       ); nglobal_field    = 0
  ALLOCATE( global_method(nfield)                       ); global_method    = 0
  ALLOCATE( ndim_field(nfield)                          ); ndim_field       = 0
  ALLOCATE( lam_ising(max_isingmax,nfield)              ); lam_ising        = 0d0
  ALLOCATE( gamma_ising(max_isingmax,nfield)            ); gamma_ising      = 0d0
  ALLOCATE( fmat(max_isingmax,max_isingmax,nfield,nflv) ); fmat             = 0d0
  ALLOCATE( nb_field(nsite,max_ndim_field,nfield)       ); nb_field         = 0
  ALLOCATE( mask_field_site(nsite,nfield)               ); mask_field_site  = .false.
  ALLOCATE( mask_field(nfield)                          ); mask_field       = .false.
  
  ifield=0

  DO i_g=1,n_g

    ifield=ifield+1
    
    READ(10,*) type_field(ifield), n_checkerboard

    SELECT CASE(type_field(ifield))
    CASE(1)

      READ(10,*) g_field(1:2,ifield)
      IF(abs(g_field(1,ifield))>1d-6) mask_field(ifield)=.true.
      isingmax(ifield)=2
      CALL HS1(ga1,lam1,exp(-dtau*g_field(1,ifield)/2))
      lam_ising(1:2,ifield)=(/lam1,-lam1/)
      gamma_ising(1:2,ifield)=(/ga1*exp(-lam1*g(2,ifield)),ga1*exp(lam1*g(2,ifield))/)

    CASE(2)

      READ(10,*) g_field(1:2,ifield)
      IF(abs(g_field(1,ifield))>1d-6) mask_field(ifield)=.true.
      isingmax(i_ising)=4
      CALL HS2(ga1,lam1,ga2,lam2,exp(-dtau*g_field(1,ifield)/2))
      lam_ising(1:4,ifield)=(/lam1,-lam1,lam2,-lam2/)
      gamma_ising(1:4,ifield)=(/ga1*exp(-lam1*g_field(2,ifield)),ga1*exp(lam1*g_field(2,ifield)), &
        &                       ga2*exp(-lam2*g_field(2,ifield)),ga2*exp(lam2*g_field(2,ifield))/)

    CASE(3)
      
      READ(10,*) g_field(1:2,ifield)
      IF(abs(g_field(1,ifield))>1d-6) mask_field(ifield)=.true.
      isingmax(i_ising)=4
      CALL HSgeneral(ga1,lam1,ga2,lam2,-dtau*g_field(1,ifield)/2)
      lam_ising(1:4,ifield)=(/lam1,-lam1,lam2,-lam2/)
      gamma_ising(1:4,ifield)=(/ga1*exp(-lam1*g_field(2,ifield)),ga1*exp(lam1*g_field(2,ifield)), &
        &                       ga2*exp(-lam2*g_field(2,ifield)),ga2*exp(lam2*g_field(2,ifield))/)
    
    CASE(-1)

      READ(10,*) g_field(1:2,ifield)
      IF(abs(g_field(1,ifield))>1d-6) mask_field(ifield)=.true.
      READ(10,*) dphi(ifield), dphi_global(ifield)

    CASE(-2)

      ! Vph, f0, debye
      READ(10,*) g_field(1:3,ifield)
      IF(abs(g_field(1,ifield))>1d-6) mask_field(ifield)=.true.
      READ(10,*) dphi(ifield), dphi_global(ifield)

    CASE(100:)
      
      READ(10,*) isingmax(ifield)
      CALL set_field_external(type_field(ifield),ifield)

    CASE(:-100)
      READ(10,*) dphi(ifield), dphi_global(ifield)
      CALL set_field_external(type_field(ifield),ifield)

    CASE default

      IF(id==0) PRINT*,'undefined yet'
      CALL exit(0)

    END SELECT

    !
    READ(10,*) nglobal_field(ifield),global_method(ifield)

    !
    READ(10,*) ndim_field(ifield)
  
    ! read in fmat
    DO flv=1,nflv
      DO k=1,ndim_field(ifield)
        READ(10,*) fmat(k,1:ndim_field(ifield),ifield,flv)
      END DO
    END DO

    ! read in the subspace
    DO k=1,ndim_field(ifield)
      READ(10,*) da,db,dc,orb,orb2
      DO a=1,La; DO b=1,Lb; DO c=1,Lc; i=label(da,db,dc,orb)
        a2=a+da; b2=b+db; c2=c+dc

        IF(a2<1 .and.(.not.pbca)) CYCLE
        IF(a2>La.and.(.not.pbca)) CYCLE
        IF(b2<1 .and.(.not.pbcb)) CYCLE
        IF(b2>Lb.and.(.not.pbcb)) CYCLE
        IF(c2<1 .and.(.not.pbcc)) CYCLE
        IF(c2>Lc.and.(.not.pbcc)) CYCLE
      
        IF(a2<1 .and.pbca) a2=a2+La
        IF(a2>La.and.pbca) a2=a2-La
        IF(b2<1 .and.pbcb) b2=b2+Lb
        IF(b2>Lb.and.pbcb) b2=b2-Lb
        IF(c2<1 .and.pbcc) c2=c2+Lc
        IF(c2>Lc.and.pbcc) c2=c2-Lc
        i2=label(a2,b2,c2,orb2)
        nb_field(i,k,ifield)=i2

      END DO; END DO; END DO
    END DO

    ! make n_checkerboard-1 copies of ifield
    DO i=1,n_checkerboard-1
      type_field(ifield+i)    = type_field(ifield)
      g_field(:,ifield+i)     = g_field(:,ifield)
      isingmax(ifield+i)      = isingmax(ifield)
      dphi(ifield+i)          = dphi(ifield)
      dphi_global(ifield+i)   = dphi_global(ifield)
      nglobal_field(ifield+i) = nglobal_field(ifield)
      global_method(ifield+i) = global_method(ifield)
      ndim_field(ifield+i)    = ndim_field(ifield)
      lam_ising(:,ifield+i)   = lam_ising(:,ifield)
      gamma_ising(:,ifield+i) = gamma_ising(:,ifield)
      fmat(:,:,ifield+i,:)    = fmat(:,:,ifield,:)
      nb_field(:,:,ifield+i)  = nb_field(:,:,ifield)
    END DO

    ! the only difference between these n_checkerboard boson fields is mask_field_site
    DO i=1,n_checkerboard
      READ(10,*) ntmp   ! how many conditions to define this HS field
      DO k=1,ntmp
        READ(10,*) ma,moda,mb,modb,mc,modc,orb
        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          IF(mod(a,ma)==moda.and.mod(b,mb)==modb.and.mod(c,mc)==modc)THEN
            mask_field_site(label(a,b,c,orb),ifield+i)=.true.
          END IF
        END DO; END DO; END DO
      END DO
    END DO

    ifield=ifield+n_checkerboard-1

  END DO
  
  CALL set_expf()

  !--------------------------------
  ! read in measurement parameters
  !--------------------------------
  
  CALL set_pool(nbin,0,100)

  CLOSE(10)

  !-------------------------------
  ! other settings
  !-------------------------------
  ALLOCATE(field(nsite,ntime,nfield))
  DO ifield=1,nfield
    CALL init_field_random(ifield)
  END DO

  IF(id==0) CALL print_input()

END SUBROUTINE

SUBROUTINE print_input()
  USE dqmc
  IMPLICIT NONE
  OPEN(10,FILE='input.dat')
  WRITE(10,*) 'This is a place to print input parameters...'
  CLOSE(10)
END SUBROUTINE

  !--------------------------------------------------------------
  ALLOCATE(gph_x2(nphi),gph_p2(nphi))
  gph_x2=1d0/2/Vph/dtau
  gph_p2=1d0/2/Vph/dtau**3/debye**2


