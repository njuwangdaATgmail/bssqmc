! This is an interactive engine for the users to prepare the input file 'dqmc.in'.
! This program is standalone.
! to be added: 
! ? add commonly-used interactions smartly
! ? add commonly-used correlation functions smartly
PROGRAM main
  IMPLICIT NONE
  COMPLEX(8), PARAMETER :: one=(1d0,0d0)
  COMPLEX(8), PARAMETER :: zero=(0d0,0d0)
  COMPLEX(8), PARAMETER :: uniti=(0d0,1d0)
  COMPLEX(8), DIMENSION(2,2) :: s0,s1,s2,s3,z0
  COMPLEX(8), DIMENSION(4,4) :: s0t0,s0t1,s0t2,s0t3,s1t0,s1t1,s1t2,s1t3,&
    &                           s2t0,s2t1,s2t2,s2t3,s3t0,s3t1,s3t2,s3t3
  INTEGER i1,i2,i3,i4,i5,i6,nflv,norb,i,j,k,flv,n_g,nfield,ncond,d,n_checkerboard,n2,err,method
  REAL(8) r1,r2,r3,r4
  COMPLEX(8) z1,z2,z3,z4
  LOGICAL l1,l2,l3,l4
  CHARACTER(120) str
  LOGICAL lvec(100)
  INTEGER ivec(100)
  REAL(8) rvec(100)
  INTERFACE saferead
    PROCEDURE saferead_logical, saferead_logical_array
    PROCEDURE saferead_integer, saferead_integer_array
    PROCEDURE saferead_real, saferead_real_array
    PROCEDURE saferead_integer_real_array
    PROCEDURE saferead_integer_character, saferead_character
    PROCEDURE saferead_integer_array_character
    PROCEDURE saferead_logical_integer
  END INTERFACE
  
  z0=0d0
  s0=reshape((/one,zero,zero,one/),(/2,2/))
  s1=reshape((/zero,one,one,zero/),(/2,2/))
  s2=reshape((/zero,uniti,-uniti,zero/),(/2,2/))
  s3=reshape((/one,zero,zero,-one/),(/2,2/))
  
  s0t0=0d0; s0t0(1:2,1:2)=s0;         s0t0(3:4,3:4)=s0
  s0t1=0d0; s0t1(1:2,3:4)=s0;         s0t1(3:4,1:2)=s0
  s0t2=0d0; s0t2(1:2,3:4)=-uniti*s0;  s0t2(3:4,1:2)=uniti*s0
  s0t3=0d0; s0t3(1:2,1:2)=s0;         s0t3(3:4,3:4)=-s0

  s1t0=0d0; s1t0(1:2,1:2)=s1;         s1t0(3:4,3:4)=s1
  s1t1=0d0; s1t1(1:2,3:4)=s1;         s1t1(3:4,1:2)=s1
  s1t2=0d0; s1t2(1:2,3:4)=-uniti*s1;  s1t2(3:4,1:2)=uniti*s1
  s1t3=0d0; s1t3(1:2,1:2)=s1;         s1t3(3:4,3:4)=-s1

  s2t0=0d0; s2t0(1:2,1:2)=s2;         s2t0(3:4,3:4)=s2
  s2t1=0d0; s2t1(1:2,3:4)=s2;         s2t1(3:4,1:2)=s2
  s2t2=0d0; s2t2(1:2,3:4)=-uniti*s2;  s2t2(3:4,1:2)=uniti*s2
  s2t3=0d0; s2t3(1:2,1:2)=s2;         s2t3(3:4,3:4)=-s2

  s3t0=0d0; s3t0(1:2,1:2)=s3;         s3t0(3:4,3:4)=s3
  s3t1=0d0; s3t1(1:2,3:4)=s3;         s3t1(3:4,1:2)=s3
  s3t2=0d0; s3t2(1:2,3:4)=-uniti*s3;  s3t2(3:4,1:2)=uniti*s3
  s3t3=0d0; s3t3(1:2,1:2)=s3;         s3t3(3:4,3:4)=-s3


  OPEN(10,file='dqmc.in')

  WRITE(10,*) '!------------------------------------------------!'
  WRITE(10,*) '!           Monte Carlo control block            !'
  WRITE(10,*) '!------------------------------------------------!'

  PRINT*,'>>restart: restart mode or not'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(1l40,10x,a)') l1,'!restart'
    

  PRINT*,'>>proj: projector mode or not'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(1l40,10x,a)') l1,'!proj'

  PRINT*,'>>nflv: number of flavors'
  i1=1
  CALL saferead(i1); nflv=i1
  WRITE(10,'(1i40,10x,a)') i1,'!nflv'

  PRINT*,'>>beta: inverse of temperature or projecting time'
  r1=1d0
  CALL saferead(r1)
  WRITE(10,'(1f40.6,10x,a)') r1,'!beta'

  PRINT*,'>>dtau: time slice distance'
  r1=0.1d0
  CALL saferead(r1)
  WRITE(10,'(1f40.6,10x,a)') r1,'!dtau'

  PRINT*,'>>nsp: how many measurements during each time-sweep'
  i1=1
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!nsp'

  PRINT*,'>>nbin: how many data bins on each core'
  i1=10
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!nbin'

  PRINT*,'>>nwarmup: warm up steps'
  i1=1000
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!nwarmup'

  PRINT*,'>>nmeasure: measurement steps of each bin'
  i1=1000
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!nmeasure'

  PRINT*,'>>ninterval: perform measurement every ninterval sweeps'
  i1=1
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!ninterval'

  PRINT*,'>>ntmpout: output internal MC information every ntmpout sweeps'
  i1=100
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!ntmpout'

  PRINT*,'>>nscratch: calculate Green function from scratch every nscratch sweeps'
  i1=10
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!nscratch'
  
  PRINT*,'>>ngroup: steps of direct matrix product between QDR/LDQ decompositions'
  i1=10
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!ngroup'

  PRINT*,'>>randomseed: seed of random number generator'
  i1=12321
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!randomseed'

  PRINT*,'>>newMetro: factor to define the modified Metropolis acceptance ratio &
    &(from 0.0 - Metropolis to 1.0 - heat-bath)'
  r1=0d0
  CALL saferead(r1)
  WRITE(10,'(1f40.6,10x,a)') r1,'!newMetro'

  WRITE(10,*)
  WRITE(10,*) '!------------------------------------------------!'
  WRITE(10,*) '!           fermion lattice block                !'
  WRITE(10,*) '!------------------------------------------------!'
  
  PRINT*,'>>norb: number of fermions in each unit cell'
  i1=1
  CALL saferead(i1); norb=i1
  WRITE(10,'(1i40,10x,a)') i1,'!norb'

  PRINT*,'>>La,Lb,Lc: lattice size'
  ivec(1:3)=(/4,4,1/)
  CALL saferead(3,ivec)
  WRITE(10,'(10x,3i10,10x,a)') ivec(1:3),'!La,Lb,Lc'

  PRINT*,'>>nelec: number of filled electrons for each flv'
  ivec(1:nflv)=ivec(1)*ivec(2)*ivec(3)*norb/2
  CALL saferead(nflv,ivec)
  SELECT CASE(nflv)
  CASE(1)
    WRITE(10,'(30x,1i10,10x,a)') ivec(1:nflv),'!nelec'
  CASE(2)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:nflv),'!nelec'
  CASE(3)
    WRITE(10,'(10x,3i10,10x,a)') ivec(1:nflv),'!nelec'
  CASE(4)
    WRITE(10,'(4i10,10x,a)') ivec(1:nflv),'!nelec'
  CASE DEFAULT
    WRITE(10,'(1000i10,10x)',advance='no') ivec(1:nflv)
    WRITE(10,*) '!nelec'
  END SELECT

  PRINT*,'>>ncopy(nflv): number of copies for each flavor'
  ivec(1:nflv)=1
  CALL saferead(nflv,ivec)
  SELECT CASE(nflv)
  CASE(1)
    WRITE(10,'(30x,1i10,10x,a)') ivec(1:nflv),'!ncopy'
  CASE(2)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:nflv),'!ncopy'
  CASE(3)
    WRITE(10,'(10x,3i10,10x,a)') ivec(1:nflv),'!ncopy'
  CASE(4)
    WRITE(10,'(4i10,10x,a)') ivec(1:nflv),'!ncopy'
  CASE DEFAULT
    WRITE(10,'(1000i10,10x)',advance='no') ivec(1:nflv)
    WRITE(10,*) '!ncopy'
  END SELECT

  PRINT*,'>>pbca,pbcb,pbcc: periodic boundary condition or not'
  lvec(1:3)=(/.true.,.true.,.false./)
  CALL saferead(3,lvec)
  WRITE(10,'(10x,3l10,10x,a)') lvec(1:3),'!pbca,pbcb,pbcc'

  PRINT*,'>>twista,twistb,twistc: twisted boundary condition in unit of 2pi'
  rvec(1:3)=0d0
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),'!twista,twistb,twistc in unit of 2pi'

  PRINT*,'>>a0r: unit cell axis a'
  rvec(1:3)=(/1d0,0d0,0d0/)
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),'!a0r'

  PRINT*,'>>b0r: unit cell axis b'
  rvec(1:3)=(/0d0,1d0,0d0/)
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),'!b0r'
  
  PRINT*,'>>c0r: unit cell axis c'
  rvec(1:3)=(/0d0,0d0,1d0/)
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),'!c0r'
  
  DO i=1,norb
    PRINT*,'>>position of the orbital',i
    rvec(1:3)=(/0d0,0d0,0d0/)
    CALL saferead(3,rvec)
    WRITE(10,'(10x,3f10.6,10x,a,1i2)') rvec(1:3),'!rorb-',i
  END DO

  PRINT*,'>>cuta,cutb,cutc: hopping range along 3 directions'
  ivec(1:3)=(/1,1,1/)
  CALL saferead(3,ivec)
  WRITE(10,'(10x,3i10,10x,a)') ivec(1:3),'!cuta,cutb,cutc'

  PRINT*,'>>nhop: number of nonzero hoppings (including onsite energies)'
  i1=0
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!nhop'

  DO i=1,i1
    PRINT*,'>>da,db,dc,orb1,orb2,Re(t),Im(t):'
    ivec(1:5)=(/0,0,0,1,1/)
    rvec(1:2)=(/-1d0,0d0/)
    CALL saferead(5,ivec,2,rvec)
    WRITE(10,'(5i4,2f10.6,10x,a)') ivec(1:5),rvec(1:2),'!da,db,dc,orb1,orb2,Re(t),Im(t)'
  END DO

  WRITE(10,*)
  WRITE(10,*) '!------------------------------------------------!'
  WRITE(10,*) '!           boson field block                    !'
  WRITE(10,*) '!------------------------------------------------!'

  PRINT*,'>>n_g,nfield,max_ndim_field,max_isingmax:'
  PRINT*,'>>  - n_g counts the kinds of interactions/bosons'
  PRINT*,'>>  - nfield is the number of practical boson fields'
  PRINT*,'>>  - max_ndim_field is the max dimension of all boson fields'
  PRINT*,'>>  - max_isingmax is the max value of all ising fields'
  ivec(1:4)=(/1,1,1,2/)
  CALL saferead(4,ivec); n_g=ivec(1); nfield=ivec(2)
  WRITE(10,'(4i10,10x,a)') ivec(1:4),'!n_g,nfield,max_ndim_field,max_isingmax'

  DO i=1,n_g
    
    PRINT*,'>>type_field, n_checkerboard, name for i_g=',i
    PRINT*,'>>  - type_field = 1(HS1), 2(HS2), 3(HSgeneral), -1(HScontinuous), -2(phonon)'
    PRINT*,'>>                 100+(user-defined ising field), -100-(user-defined continuous field)'
    PRINT*,'>>  - n_checkerboard counts the number of independent boson fields belong to this kind of interaction'
    ivec(1:2)=(/1,1/)
    IF(i<=9)THEN
      WRITE(str,'(1a,1i1)') 'interaction_',i
    ELSEIF(i<=99)THEN
      WRITE(str,'(1a,1i2)') 'interaction_',i
    ELSE
      WRITE(str,'(1a,1i3)') 'interaction_',i
    END IF
    CALL saferead(2,ivec,str); n_checkerboard=ivec(2)
    WRITE(10,*)
    WRITE(10,'(50x,a)') '!setting block of '//trim(adjustl(str))
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!type_field, n_checkerboard'
    
    PRINT*,'>>g_field for i_g=',i
    PRINT*,">>  - g_field is defined as H = g(1)/2 * ( c' * fmat *c - g(2) )^2"
    PRINT*,">>  - for phonon, it is defined as H = eta * (c' * fmat * c - g(2) )^2"
    PRINT*,">>    and g(3)=Debye frequency, g(1)=eta^2/M/g(3)^2"
    rvec(1:3)=(/4d0,1d0,0d0/)
    CALL saferead(3,rvec)
    WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),'!g_field'
    
    PRINT*,'>>dphi and dphi_global for i_g=',i
    PRINT*,'>>  - they are the maximal shift of the continuous boson field in local and global updates, respectively'
    rvec(1:4)=(/0d0,0d0,0d0,0d0/)
    CALL saferead(4,rvec)
    WRITE(10,'(4f10.6,10x,a)') rvec(1:4),'!dphi, dphi_global'
    
    PRINT*,'>>ninterval_global and global_method for i_g=',i
    PRINT*,'>>  - do global update every ninterval_global sweeps'
    PRINT*,'>>  - global method = 0(totally random init), 1(totally random shift), 2(random shift at each time), 100(user-defined)'
    ivec(1:2)=(/100000,0/)
    CALL saferead(2,ivec)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!ninterval_global, global_method'
    
    PRINT*,'>>ndim_field: dimension of each boson field'
    i1=1
    CALL saferead(i1); d=i1
    WRITE(10,'(1i40,10x,a)') i1,'!ndim_field'
    
    DO flv=1,nflv
      WRITE(10,'(50x,a,1i4)') '!fmat(:,:) for flv=',flv
      PRINT*,'>>fmat of flv=',flv,'! choose si or sitj (0<=i,j<=3) or input manually'
      str='by hand'
      CALL saferead(str)
      SELECT CASE(trim(str))
      CASE('s0','S0')
        WRITE(10,'(20x,4f5.1)') transpose(s0)
      CASE('s1','S1')
        WRITE(10,'(20x,4f5.1)') transpose(s1)
      CASE('s2','S2')
        WRITE(10,'(20x,4f5.1)') transpose(s2)
      CASE('s3','S3')
        WRITE(10,'(20x,4f5.1)') transpose(s3)
      CASE('s0t0','S0T0')
        WRITE(10,'(8f5.1)') transpose(s0t0)
      CASE('s0t1','S0T1')
        WRITE(10,'(8f5.1)') transpose(s0t1)
      CASE('s0t2','S0T2')
        WRITE(10,'(8f5.1)') transpose(s0t2)
      CASE('s0t3','S0T3')
        WRITE(10,'(8f5.1)') transpose(s0t3)
      CASE('s1t0','S1T0')
        WRITE(10,'(8f5.1)') transpose(s1t0)
      CASE('s1t1','S1T1')
        WRITE(10,'(8f5.1)') transpose(s1t1)
      CASE('s1t2','S1T2')
        WRITE(10,'(8f5.1)') transpose(s1t2)
      CASE('s1t3','S1T3')
        WRITE(10,'(8f5.1)') transpose(s1t3)
      CASE('s2t0','S2T0')
        WRITE(10,'(8f5.1)') transpose(s2t0)
      CASE('s2t1','S2T1')
        WRITE(10,'(8f5.1)') transpose(s2t1)
      CASE('s2t2','S2T2')
        WRITE(10,'(8f5.1)') transpose(s2t2)
      CASE('s2t3','S2T3')
        WRITE(10,'(8f5.1)') transpose(s2t3)
      CASE('s3t0','S3T0')
        WRITE(10,'(8f5.1)') transpose(s3t0)
      CASE('s3t1','S3T1')
        WRITE(10,'(8f5.1)') transpose(s3t1)
      CASE('s3t2','S3T2')
        WRITE(10,'(8f5.1)') transpose(s3t2)
      CASE('s3t3','S3T3')
        WRITE(10,'(8f5.1)') transpose(s3t3)
      CASE DEFAULT
        DO k=1,d
          PRINT*,'>>fmat(k,:,flv): k=',k,' flv=',flv
          rvec(1:2*d)=0d0; rvec(2*d-1)=1d0
          CALL saferead(2*d,rvec)
          SELECT CASE(d)
          CASE(1)
            WRITE(10,'(30x,2f5.1)') rvec(1:2*d)
          CASE(2)
            WRITE(10,'(20x,4f5.1)') rvec(1:2*d)
          CASE(3)
            WRITE(10,'(10x,6f5.1)') rvec(1:2*d)
          CASE(4)
            WRITE(10,'(8f5.1)') rvec(1:2*d)
          CASE DEFAULT
            WRITE(10,'(1000f5.1)') rvec(1:2*d)
          END SELECT
        END DO
      END SELECT
    END DO


    DO k=1,n_checkerboard

      PRINT*,'>>for checkboard',k
      WRITE(10,'(50x,a,1i4)'),'!for checkboard',k
    
      DO j=1,d
        PRINT*,'>>da,db,dc,orb for basis',j
        ivec(1:4)=(/0,0,0,1/)
        CALL saferead(4,ivec)
        WRITE(10,'(4i10,10x,a,1i4)') ivec(1:4),'!da,db,dc,orb for basis',j
      END DO

      PRINT*,'>>how many conditions to define this checkerboard?'
      i1=1
      CALL saferead(i1); ncond=i1
      WRITE(10,'(1i40,10x,a)') i1,'!n_cond'
      PRINT*,'>>for each condition, provide: ma,moda,mb,modb,mc,modc'
      PRINT*,'>>where the boson field lives only when mod(a,ma)==moda, mod(b,mb)==modb, mod(c,mc)==modc'
      DO j=1,ncond
        ivec(1:6)=(/1,0,1,0,1,0/)
        CALL saferead(6,ivec)
        IF(j==1)THEN
          WRITE(10,'(10x,6i5,10x,a)') ivec(1:6),'!ma,moda,mb,modb,mc,modc'
        ELSE
          WRITE(10,'(10x,6i5,10x)') ivec(1:6)
        END IF
      END DO

    END DO

  END DO

  WRITE(10,*)
  WRITE(10,*) '!------------------------------------------------!'
  WRITE(10,*) '!           measurement block                    !'
  WRITE(10,*) '!------------------------------------------------!'

  PRINT*,'>>nk_meas,k_method: how many k-points to measure and choose input method'
  PRINT*,'>>  - method = 0: input a,b,c, while ka=2*pi*a/La, kb=2*pi*b/Lb, kc=2*pi*c/Lc'
  PRINT*,'>>          /= 0: input ka,kb,kc in unit of 2*pi'
  ivec(1:2)=0
  CALL saferead(2,ivec); k=ivec(1); method=ivec(2)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!nk_meas,k_method'
  
  IF(k>0) PRINT*,'>>input the k_array below:'
  DO i=1,k
    IF(method==0)THEN
      ivec(1:3)=(/0,0,0/)
      CALL saferead(3,ivec)
      WRITE(10,'(10x,3i10,10x)') ivec(1:3)
    ELSE
      rvec(1:3)=(/0d0,0d0,0d0/)
      CALL saferead(3,rvec)
      WRITE(10,'(10X,3f10.4,10x)') rvec(1:3)
    END IF
  END DO
  
  WRITE(10,*)
  PRINT*,'>>nr_meas,r_method: how many (relative) r-points to measure and choose input method'
  PRINT*,'>>  - method = 0: input a,b,c directly'
  PRINT*,'>>          /= 0: input a/La, b/Lb, c/Lc'
  ivec(1:2)=0
  CALL saferead(2,ivec); k=ivec(1); method=ivec(2)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!nr_meas,r_method'
  
  IF(k>0) PRINT*,">>input the r_array below:"
  DO i=1,k
    IF(method==0)THEN
      ivec(1:3)=(/0,0,0/)
      CALL saferead(3,ivec)
      WRITE(10,'(10x,3i10,10x)') ivec(1:3)
    ELSE
      rvec(1:3)=(/0d0,0d0,0d0/)
      CALL saferead(3,rvec)
      WRITE(10,'(10x,3f10.4,10x)') rvec(1:3)
    END IF
  END DO

  WRITE(10,*)
  PRINT*,'>>nrr_meas,rr_method: how many rr-points to measure and choose input method'
  PRINT*,'>>                    (useful in models without translation symmetry)'
  PRINT*,'>>  - method = 0: input a,b,c,a2,b2,c2 directly'
  PRINT*,'>>          /= 0: input a/La,b/Lb,c/Lc,a2/La,b2/Lb,c2/Lc'
  ivec(1:2)=0
  CALL saferead(2,ivec); k=ivec(1); method=ivec(2)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!nrr_meas,rr_method'
  
  IF(k>0) PRINT*,">>input the rr_array below:"
  DO i=1,k
    IF(method==0)THEN
      ivec(1:6)=(/1,1,1,1,1,1/)
      CALL saferead(6,ivec)
      WRITE(10,'(10x,6i5)') ivec(1:6)
    ELSE
      rvec(1:6)=0d0
      CALL saferead(6,rvec)
      WRITE(10,'(4x,6f6.3,10x)') rvec(1:6)
    END IF
  END DO

  WRITE(10,*)
  PRINT*,'>>ntau_meas: time displaced Green functions are measured during [-ntau_meas,ntau_meas]'
  i1=0
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,'!ntau_meas'

  WRITE(10,*)
  PRINT*,'>>n_ph_meas, max_ndim_ph_meas: how many PH-channel two-particle Green functions to measure'
  PRINT*,'>>                             and the maximal dimension of fmat given below'
  ivec(1:2)=(/0,0/)
  CALL saferead(2,ivec); n2=ivec(1)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!n_ph_meas,max_ndim_ph_meas'

  DO i=1,n2
    
    WRITE(10,'(50x,a,1i4)') '!setting block for PH-',i

    PRINT*,'>>hartree_ph_meas, fork_ph_meas for PH-', i
    lvec(1:2)=(/.true.,.true./)
    CALL saferead(2,lvec)
    WRITE(10,'(20x,2l10,10x,a)') lvec(1:2),'!hartree_ph_meas,fork_ph_meas'

    PRINT*,'>>ndim_ph_meas and name_ph_meas for PH-', i
    i1=nflv; str='default_ph'
    CALL saferead(i1,str); d=i1
    WRITE(10,'(10x,1i10,1a20,10x,a)') i1,trim(adjustl(str)),'!ndim_ph_meas, name_ph_meas'
    
    PRINT*,'>>for each basis, provide: da,db,dc,orb,flv'
    DO k=1,d
      ivec(1:5)=(/0,0,0,1,1/)
      CALL saferead(5,ivec)
      IF(k==1)THEN
        WRITE(10,'(15x,5i5,10x,a)') ivec(1:5),'!da,db,dc,orb,flv'
      ELSE
        WRITE(10,'(15x,5i5,10x)') ivec(1:5)
      END IF
    END DO

    WRITE(10,'(50x,a,1i4)') '!fmat_ph_meas(:,:) for PH-',i
    PRINT*,'>>fmat of PH-',i,'! choose si or sitj (0<=i,j<=3) or input manually'
    str='by hand'
    CALL saferead(str)
    SELECT CASE(trim(str))
    CASE('s0','S0')
      WRITE(10,'(20x,4f5.1)') transpose(s0)
    CASE('s1','S1')
      WRITE(10,'(20x,4f5.1)') transpose(s1)
    CASE('s2','S2')
      WRITE(10,'(20x,4f5.1)') transpose(s2)
    CASE('s3','S3')
      WRITE(10,'(20x,4f5.1)') transpose(s3)
    CASE('s0t0','S0T0')
      WRITE(10,'(8f5.1)') transpose(s0t0)
    CASE('s0t1','S0T1')
      WRITE(10,'(8f5.1)') transpose(s0t1)
    CASE('s0t2','S0T2')
      WRITE(10,'(8f5.1)') transpose(s0t2)
    CASE('s0t3','S0T3')
      WRITE(10,'(8f5.1)') transpose(s0t3)
    CASE('s1t0','S1T0')
      WRITE(10,'(8f5.1)') transpose(s1t0)
    CASE('s1t1','S1T1')
      WRITE(10,'(8f5.1)') transpose(s1t1)
    CASE('s1t2','S1T2')
      WRITE(10,'(8f5.1)') transpose(s1t2)
    CASE('s1t3','S1T3')
      WRITE(10,'(8f5.1)') transpose(s1t3)
    CASE('s2t0','S2T0')
      WRITE(10,'(8f5.1)') transpose(s2t0)
    CASE('s2t1','S2T1')
      WRITE(10,'(8f5.1)') transpose(s2t1)
    CASE('s2t2','S2T2')
      WRITE(10,'(8f5.1)') transpose(s2t2)
    CASE('s2t3','S2T3')
      WRITE(10,'(8f5.1)') transpose(s2t3)
    CASE('s3t0','S3T0')
      WRITE(10,'(8f5.1)') transpose(s3t0)
    CASE('s3t1','S3T1')
      WRITE(10,'(8f5.1)') transpose(s3t1)
    CASE('s3t2','S3T2')
      WRITE(10,'(8f5.1)') transpose(s3t2)
    CASE('s3t3','S3T3')
      WRITE(10,'(8f5.1)') transpose(s3t3)
    CASE DEFAULT
      DO k=1,d
        PRINT*,'>>fmat_ph_meas(k,:): k=',k
        rvec(1:2*d)=0d0; rvec(2*k-1)=1d0
        CALL saferead(2*d,rvec)
        SELECT CASE(d)
        CASE(1)
          WRITE(10,'(30x,2f5.1)') rvec(1:2*d)
        CASE(2)
          WRITE(10,'(20x,4f5.1)') rvec(1:2*d)
        CASE(3)
          WRITE(10,'(10x,6f5.1)') rvec(1:2*d)
        CASE(4)
          WRITE(10,'(8f5.1)') rvec(1:2*d)
        CASE DEFAULT
          WRITE(10,'(1000f5.1)') rvec(1:2*d)
        END SELECT
      END DO
    END SELECT
  END DO
  
  WRITE(10,*)
  PRINT*,'>>n_pp_meas, max_ndim_pp_meas: how many PP-channel two-particle Green functions to measure'
  PRINT*,'>>                             and the maximal dimension of fmat given below'
  ivec(1:2)=(/0,0/)
  CALL saferead(2,ivec); n2=ivec(1)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!n_pp_meas,max_ndim_pp_meas'

  DO i=1,n2
    
    WRITE(10,'(50x,a,1i4)') '!setting block for PP-',i

    PRINT*,'>>fork13_pp_meas, fork14_pp_meas for PP-', i
    lvec(1:2)=(/.true.,.true./)
    CALL saferead(2,lvec)
    WRITE(10,'(20x,2l10,10x,a)') lvec(1:2),'!fork13_pp_meas,fork14_pp_meas'

    PRINT*,'>>ndim_pp_meas and name_pp_meas for PP-', i
    i1=nflv; str='default_pp'
    CALL saferead(i1,str); d=i1
    WRITE(10,'(10x,1i10,1a20,10x,a)') i1,trim(adjustl(str)),'!ndim_pp_meas, name_pp_meas'
    
    PRINT*,'>>for each basis, provide: da,db,dc,orb,flv'
    DO k=1,d
      ivec(1:5)=(/0,0,0,1,1/)
      CALL saferead(5,ivec)
      IF(k==1)THEN
        WRITE(10,'(15x,5i5,10x,a)') ivec(1:5),'!da,db,dc,orb,flv'
      ELSE
        WRITE(10,'(15x,5i5,10x)') ivec(1:5)
      END IF
    END DO

    WRITE(10,'(50x,a,1i4)') '!fmat_pp_meas(:,:) for PP-',i
    PRINT*,'>>fmat of PP-',i,'! choose si or sitj (0<=i,j<=3) or input manually'
    str='by hand'
    CALL saferead(str)
    SELECT CASE(trim(str))
    CASE('s0','S0')
      WRITE(10,'(20x,4f5.1)') transpose(s0)
    CASE('s1','S1')
      WRITE(10,'(20x,4f5.1)') transpose(s1)
    CASE('s2','S2')
      WRITE(10,'(20x,4f5.1)') transpose(s2)
    CASE('s3','S3')
      WRITE(10,'(20x,4f5.1)') transpose(s3)
    CASE('s0t0','S0T0')
      WRITE(10,'(8f5.1)') transpose(s0t0)
    CASE('s0t1','S0T1')
      WRITE(10,'(8f5.1)') transpose(s0t1)
    CASE('s0t2','S0T2')
      WRITE(10,'(8f5.1)') transpose(s0t2)
    CASE('s0t3','S0T3')
      WRITE(10,'(8f5.1)') transpose(s0t3)
    CASE('s1t0','S1T0')
      WRITE(10,'(8f5.1)') transpose(s1t0)
    CASE('s1t1','S1T1')
      WRITE(10,'(8f5.1)') transpose(s1t1)
    CASE('s1t2','S1T2')
      WRITE(10,'(8f5.1)') transpose(s1t2)
    CASE('s1t3','S1T3')
      WRITE(10,'(8f5.1)') transpose(s1t3)
    CASE('s2t0','S2T0')
      WRITE(10,'(8f5.1)') transpose(s2t0)
    CASE('s2t1','S2T1')
      WRITE(10,'(8f5.1)') transpose(s2t1)
    CASE('s2t2','S2T2')
      WRITE(10,'(8f5.1)') transpose(s2t2)
    CASE('s2t3','S2T3')
      WRITE(10,'(8f5.1)') transpose(s2t3)
    CASE('s3t0','S3T0')
      WRITE(10,'(8f5.1)') transpose(s3t0)
    CASE('s3t1','S3T1')
      WRITE(10,'(8f5.1)') transpose(s3t1)
    CASE('s3t2','S3T2')
      WRITE(10,'(8f5.1)') transpose(s3t2)
    CASE('s3t3','S3T3')
      WRITE(10,'(8f5.1)') transpose(s3t3)
    CASE DEFAULT
      DO k=1,d
        PRINT*,'>>fmat_pp_meas(k,:): k=',k
        rvec(1:2*d)=0d0; rvec(2*k-1)=1d0
        CALL saferead(2*d,rvec)
        SELECT CASE(d)
        CASE(1)
          WRITE(10,'(30x,2f5.1)') rvec(1:2*d)
        CASE(2)
          WRITE(10,'(20x,4f5.1)') rvec(1:2*d)
        CASE(3)
          WRITE(10,'(10x,6f5.1)') rvec(1:2*d)
        CASE(4)
          WRITE(10,'(8f5.1)') rvec(1:2*d)
        CASE DEFAULT
          WRITE(10,'(1000f5.1)') rvec(1:2*d)
        END SELECT
      END DO
    END SELECT
  END DO

  WRITE(10,*)
  PRINT*,">>ncross_ph_meas: how many crossing-PH-channel measurements G(ph1,ph2)=<O(ph1)O'(ph2)>"
  i1=0
  CALL saferead(i1); n2=i1
  WRITE(10,'(30x,1i10,10x,a)') i1,'!ncross_ph_meas'

  DO i=1,n2
    PRINT*,'>>name of the channel'
    str='default_cross_ph'
    CALL saferead(str)
    WRITE(10,'(20x,1a20,10x,a,1i4)') trim(adjustl(str)),'!name of the channel',i
    PRINT*,'>>which two channels are correlated?'
    ivec(1:2)=(/1,2/)
    CALL saferead(2,ivec)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!cross_ph_meas'
  END DO
  
  WRITE(10,*)
  PRINT*,">>ncross_pp_meas: how many crossing-PP-channel measurements G(pp1,pp2)=<O(pp1)O'(pp2)>"
  i1=0
  CALL saferead(i1); n2=i1
  WRITE(10,'(30x,1i10,10x,a)') i1,'!ncross_pp_meas'

  DO i=1,n2
    PRINT*,'>>name of the channel'
    str='default_cross_pp'
    CALL saferead(str)
    WRITE(10,'(20x,1a20,10x,a,1i4)') trim(adjustl(str)),'!name of the channel',i
    PRINT*,'>>which two channels are correlated?'
    ivec(1:2)=(/1,2/)
    CALL saferead(2,ivec)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),'!cross_pp_meas'
  END DO


  WRITE(10,*)
  PRINT*,'>>FAtech: whether to use Feldbacher-Assaad stablization algorithm for T=0'
  l1=.true.
  CALL saferead(l1)
  WRITE(10,'(1l40,10x,a)') l1,'!FAtech'

  PRINT*,'>>do_measurement_external, n_meas_external: whether or not to do and how many external measurements'
  l1=.false.; i1=0
  CALL saferead(l1,i1)
  WRITE(10,'(20x,1l10,1i10,10x,a)') l1,i1,'!do_measurement_external, n_meas_external'


  PRINT*,'>>do_tmpout_external: whether or not to call tmpout externally'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(30x,1l10,10x,a)') l1,'!do_tmpout_external'


  PRINT*,'>>do_postprocess_external: whether or not to call postprocess externally'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(30x,1l10,10x,a)') l1,'!do_postprocess_external'

  CLOSE(10)

CONTAINS

SUBROUTINE saferead_logical(x)
  IMPLICIT NONE
  LOGICAL x
  CHARACTER(120) line
  INTEGER err
  PRINT'(1x,1a,1l6)','>>press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_logical_array(n,x)
  IMPLICIT NONE
  INTEGER n
  LOGICAL x(n)
  CHARACTER(120) line
  INTEGER err
  PRINT'(1x,1a,100l6)','>>press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_integer(x)
  IMPLICIT NONE
  INTEGER x
  CHARACTER(120) line
  INTEGER err
  PRINT'(1x,1a,1i6)','>>press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_integer_array(n,x)
  IMPLICIT NONE
  INTEGER n
  INTEGER x(n)
  CHARACTER(120) line
  INTEGER err
  PRINT'(1x,1a,100i6)','>>press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_real(x)
  IMPLICIT NONE
  REAL(8) x
  CHARACTER(120) line
  INTEGER err
  PRINT'(1x,1a,1f12.6)','>>press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_real_array(n,x)
  IMPLICIT NONE
  INTEGER n
  REAL(8) x(n)
  CHARACTER(120) line
  INTEGER err
  PRINT'(1x,1a,100f12.6)','>>press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_integer_real_array(n,x,m,y)
  IMPLICIT NONE
  INTEGER n,m
  INTEGER x(n)
  REAL(8) y(m)
  CHARACTER(120) line
  INTEGER err
  PRINT*,'>>press ENTER for default value:',x,y
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE


SUBROUTINE saferead_integer_character(x,y)
  IMPLICIT NONE
  INTEGER x
  CHARACTER(120) y
  CHARACTER(120) line
  INTEGER err
  PRINT*,'>>press ENTER for default value:',x,trim(adjustl(y))
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_character(y)
  IMPLICIT NONE
  CHARACTER(120) y
  CHARACTER(120) line
  INTEGER err
  PRINT*,'>>press ENTER for default value:',trim(adjustl(y))
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) y
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_integer_array_character(n,x,y)
  IMPLICIT NONE
  INTEGER n
  INTEGER x(n)
  CHARACTER(120) y
  CHARACTER(120) line
  INTEGER err
  PRINT*,'>>press ENTER for default value:',x,trim(adjustl(y))
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE

SUBROUTINE saferead_logical_integer(x,y)
  IMPLICIT NONE
  LOGICAL x
  INTEGER y
  CHARACTER(120) line
  INTEGER err
  PRINT*,'>>press ENTER for default value:',x,y
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'>>input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE


END PROGRAM
