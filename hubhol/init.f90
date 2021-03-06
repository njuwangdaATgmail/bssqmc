SUBROUTINE init()
  USE dqmc
  USE hubhol
  IMPLICIT NONE

  INTEGER i,da,db,dc,orb,orb2,nhop
  REAL(8) re,re2,re3
  COMPLEX(8) ga,lam
  
  nising=1
  nphi=1
  nfield=nising+nphi

  !------------------------------------------------------------------
  ! read parameters from file "job.in"
  OPEN(10,file='job.in')
  READ(10,*) restart
  READ(10,*) proj
  READ(10,*) nflv; nflv=2
  READ(10,*) beta
  READ(10,*) dtau

  ntime=ceiling(beta/dtau)
  IF(mod(ntime,2)==1) ntime=ntime+1
  IF(ntime<2) ntime=2
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
  ALLOCATE(ninterval_global(nfield))
  READ(10,*) ninterval_global(:)
  READ(10,*) ntmpout
  READ(10,*) nscratch
  READ(10,*) ngroup
  READ(10,*) randomseed
  randomseed=randomseed+701703*id; CALL init_rng(randomseed)
  READ(10,*) newMetro

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

  READ(10,*) U,Vph,debye
  ALLOCATE(dphi(nphi),dphi_global(nphi))
  READ(10,*) re; dphi(1)=re
  READ(10,*) re; dphi_global(1)=re

  READ(10,*) cuta,cutb,cutc
  ALLOCATE(hop(-cuta:cuta,-cutb:cutb,-cutc:cutc,norb,norb,nflv))
  hop=0d0

  READ(10,*) nhop
  DO i=1,nhop
    READ(10,*) da,db,dc,orb,orb2,re,re2
    hop(da,db,dc,orb,orb2,:)=dcmplx(re,re2)
    hop(-da,-db,-dc,orb2,orb,:)=dcmplx(re,-re2)
  END DO

  CLOSE(10)
  !-----------------------------------------------------------

  ALLOCATE(kmat(nsite,nsite,nflv))

  CALL set_kmat(); CALL set_expk()

  IF(proj)THEN
    re=twista; re2=twistb; re3=twistc
    twista=3d0/17; twistb=3d0/13; twistc=3d0/7
    CALL set_kmat(); CALL set_slater()
    twista=re; twistb=re2; twistc=re3
  END IF

  CALL set_kmat()
  OPEN(10,FILE='kmat.dat'); WRITE(10,'(2f18.10)') kmat; CLOSE(10)

  !-----------------------------------------------------------------
  ALLOCATE(ndim_ising(nising)); ndim_ising=1
  ALLOCATE(nb_ising(nsite,maxval(ndim_ising),nising)); nb_ising(:,1,1)=(/(i,i=1,nsite)/)
  ALLOCATE(mask_ising(nising))
  IF(abs(U)>1d-6)THEN
    mask_ising=.true.
  ELSE
    mask_ising=.false.
  END IF
  ALLOCATE(mask_ising_site(nsite,nising)); mask_ising_site=.true.
  ALLOCATE(isingmax(nising)); isingmax=2
  CALL set_isingflip()
  ALLOCATE(fmat_ising(maxval(ndim_ising),maxval(ndim_ising),nising,nflv))
  fmat_ising(1,1,1,:)=(/1d0,-1d0/)
  ALLOCATE(gamma_ising(maxval(isingmax),nising),lam_ising(maxval(isingmax),nising))
  CALL HS1(ga,lam,exp(dtau*U/2))
  gamma_ising(:,1)=(/ga,ga/)
  lam_ising(:,1)=(/lam,-lam/)
  CALL set_expf_ising()

  !----------------------------------------------------------------
  ALLOCATE(ndim_phi(nphi)); ndim_phi=1
  ALLOCATE(nb_phi(nsite,maxval(ndim_ising),nphi)); nb_phi(:,1,1)=(/(i,i=1,nsite)/)
  ALLOCATE(mask_phi(nphi))
  IF(abs(Vph)>1d-6)THEN
    mask_phi=.true.
  ELSE
    mask_phi=.false.
  END IF
  ALLOCATE(mask_phi_site(nsite,nphi)); mask_phi_site=.true.
  ALLOCATE(fmat_phi(maxval(ndim_phi),maxval(ndim_phi),nphi,nflv))
  fmat_phi(1,1,1,:)=(/1d0,1d0/)
  CALL set_expf_phi()

  !---------------------------------------------------------------
  ALLOCATE(ndim_field(nfield)); ndim_field=(/1,1/)
  ALLOCATE(mask_field(nfield))
  mask_field(1)=mask_ising(1); mask_field(2)=mask_phi(1)
  ALLOCATE(mask_field_site(nsite,nfield))
  mask_field_site(:,1)=mask_ising_site(:,1)
  mask_field_site(:,2)=mask_phi_site(:,1)
  ALLOCATE(nb_field(nsite,maxval(ndim_field),nfield))
  nb_field(:,1,1)=nb_ising(:,1,1)
  nb_field(:,1,2)=nb_phi(:,1,1)
  
  !--------------------------------------------------------------
  CALL set_pool(nbin,0,100)

  !--------------------------------------------------------------
  ALLOCATE(gph_x2(nphi),gph_p2(nphi))
  gph_x2=1d0/2/Vph/dtau
  gph_p2=1d0/2/Vph/dtau**3/debye**2

  !--------------------------------------------------------------
  ALLOCATE(field(nsite,ntime,nfield))
  CALL init_ising_random(1,1)
  CALL init_phi_random(1,2)

  IF(id==0) CALL print_input()

END SUBROUTINE

SUBROUTINE print_input()
  USE dqmc
  USE hubhol
  IMPLICIT NONE
  OPEN(10,FILE='input.dat')
  WRITE(10,*) 'This is a place to print input parameters...'
  CLOSE(10)
END SUBROUTINE
