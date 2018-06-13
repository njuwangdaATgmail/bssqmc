! This is an interactive engine for the users to prepare the input file 'dqmc.in'.
! This program is standalone.
! to be added: 
! expand fmat fmat_ph_meas, fmat_pp_meas using
! + 3 Pauli matrices
! + 8 Gelmann matrices
! + 15 gamma matrices
! + ...
PROGRAM main
  IMPLICIT NONE
  INTEGER i1,i2,i3,i4,i5,i6,nflv,norb,i,j,k,flv,n_g,nfield,ncond,d,n_checkerboard,n2,err
  REAL(8) r1,r2,r3,r4
  COMPLEX(8) z1,z2,z3,z4
  LOGICAL l1,l2,l3,l4
  CHARACTER(120) str,line
  LOGICAL lvec(100)
  INTEGER ivec(100)
  REAL(8) rvec(100)
  INTERFACE saferead
    PROCEDURE saferead_logical, saferead_logical_array
    PROCEDURE saferead_integer, saferead_integer_array
    PROCEDURE saferead_real, saferead_real_array
    PROCEDURE saferead_integer_real_array
    PROCEDURE saferead_integer_character, saferead_character
    PROCEDURE saferead_logical_integer
  END INTERFACE

  OPEN(10,file='dqmc.in')

  WRITE(10,*) '#------------------------------------------------'
  WRITE(10,*) '#           Monte Carlo control block          '
  WRITE(10,*) '#------------------------------------------------'

  PRINT*,'restart: restart mode or not'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(1l40,10x,a)') l1,':restart'
    

  PRINT*,'proj: projector mode or not'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(1l40,10x,a)') l1,':proj'

  PRINT*,'nflv: number of flavors'
  i1=1
  CALL saferead(i1); nflv=i1
  WRITE(10,'(1i40,10x,a)') i1,':nflv'

  PRINT*,'beta: inverse of temperature or projecting time'
  r1=1d0
  CALL saferead(r1)
  WRITE(10,'(1f40.6,10x,a)') r1,':beta'

  PRINT*,'dtau: time slice distance'
  r1=0.1d0
  CALL saferead(r1)
  WRITE(10,'(1f40.6,10x,a)') r1,':dtau'

  PRINT*,'nsp: how many measurements during each time-sweep'
  i1=1
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':nsp'

  PRINT*,'nbin: how many data bins on each core'
  i1=10
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':nbin'

  PRINT*,'nwarmup: warm up steps'
  i1=1000
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':nwarmup'

  PRINT*,'nmeasure: measurement steps of each bin'
  i1=1000
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':nmeasure'

  PRINT*,'ninterval: perform measurement every ninterval sweeps'
  i1=1
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':ninterval'

  PRINT*,'ntmpout: output internal MC information every ntmpout sweeps'
  i1=100
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':ntmpout'

  PRINT*,'nscratch: calculate Green function from scratch every nscratch sweeps'
  i1=10
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':nscratch'
  
  PRINT*,'ngroup: steps of direct matrix product between QDR/LDQ decompositions'
  i1=10
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':ngroup'

  PRINT*,'randomseed: seed of random number generator'
  i1=12321
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':randomseed'

  PRINT*,'newMetro: factor to define the modified Metropolis acceptance ratio &
    &(from 0.0 - Metropolis to 1.0 - heat-bath)'
  r1=0d0
  CALL saferead(r1)
  WRITE(10,'(1f40.6,10x,a)') r1,':newMetro'

  WRITE(10,*) '#------------------------------------------------'
  WRITE(10,*) '#           fermion lattice block          '
  WRITE(10,*) '#------------------------------------------------'
  
  PRINT*,'norb: number of fermions in each unit cell'
  i1=1
  CALL saferead(i1); norb=i1
  WRITE(10,'(1i40,10x,a)') i1,':norb'

  PRINT*,'La,Lb,Lc: lattice size'
  ivec(1:3)=(/4,4,1/)
  CALL saferead(3,ivec)
  WRITE(10,'(10x,3i10,10x,a)') ivec(1:3),':La,Lb,Lc'

  PRINT*,'nelec: number of filled electrons for each flv'
  ivec(1:nflv)=ivec(1)*ivec(2)*ivec(3)*norb/2
  CALL saferead(nflv,ivec)
  WRITE(10,'(100i10)',advance='no') ivec(1:nflv)
  WRITE(10,'(10x,a)') ':nelec'

  PRINT*,'ncopy(nflv): number of copies for each flavor'
  ivec(1:nflv)=1
  CALL saferead(nflv,ivec)
  WRITE(10,'(100i10)',advance='no') ivec(1:nflv)
  WRITE(10,'(10x,a)') ':ncopy'

  PRINT*,'pbca,pbcb,pbcc: periodic boundary condition or not'
  lvec(1:3)=(/.true.,.true.,.false./)
  CALL saferead(3,lvec)
  WRITE(10,'(10x,3l10,10x,a)') lvec(1:3),':pbca,pbcb,pbcc'

  PRINT*,'twista,twistb,twistc: twisted boundary condition in unit of 2pi'
  rvec(1:3)=0d0
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),':twista,twistb,twistc in unit of 2pi'

  PRINT*,'a0r: unit cell axis a'
  rvec(1:3)=(/1d0,0d0,0d0/)
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),':a0r'

  PRINT*,'b0r: unit cell axis b'
  rvec(1:3)=(/0d0,1d0,0d0/)
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),':b0r'
  
  PRINT*,'c0r: unit cell axis c'
  rvec(1:3)=(/0d0,0d0,1d0/)
  CALL saferead(3,rvec)
  WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),':c0r'
  
  DO i=1,norb
    PRINT*,'position of the orbital',i
    rvec(1:3)=(/0d0,0d0,0d0/)
    CALL saferead(3,rvec)
    WRITE(10,'(10x,3f10.6,10x,a,1i2)') rvec(1:3),':rorb-',i
  END DO

  PRINT*,'cuta,cutb,cutc: hopping range along 3 directions'
  ivec(1:3)=(/1,1,1/)
  CALL saferead(3,ivec)
  WRITE(10,'(10x,3i10,10x,a)') ivec(1:3),':cuta,cutb,cutc'

  PRINT*,'nhop: number of nonzero hoppings (including onsite energies)'
  i1=0
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':nhop'

  DO i=1,i1
    PRINT*,'da,db,dc,orb1,orb2,Re(t),Im(t):'
    ivec(1:5)=(/0,0,0,1,1/)
    rvec(1:2)=(/-1d0,0d0/)
    CALL saferead(5,ivec,2,rvec)
    WRITE(10,'(5i4,2f10.6,10x,a)') ivec(1:5),rvec(1:2),':da,db,dc,orb1,orb2,Re(t),Im(t)'
  END DO

  WRITE(10,*) '#------------------------------------------------'
  WRITE(10,*) '#           boson field block          '
  WRITE(10,*) '#------------------------------------------------'

  PRINT*,'n_g,nfield,max_ndim_field,max_isingmax:'
  PRINT*,'  - n_g counts the kinds of interactions/bosons'
  PRINT*,'  - nfield is the number of practical boson fields'
  PRINT*,'  - max_ndim_field is the max dimension of all boson fields'
  PRINT*,'  - max_isingmax is the max value of all ising fields'
  ivec(1:4)=(/1,1,1,2/)
  CALL saferead(4,ivec); n_g=ivec(1); nfield=ivec(2)
  WRITE(10,'(4i10,10x,a)') ivec(1:4),':n_g,nfield,max_ndim_field,max_isingmax'

  DO i=1,n_g
    
    PRINT*,'type_field, n_checkerboard for i_g=',i
    PRINT*,'  - type_field = 1(HS1), 2(HS2), 3(HSgeneral), -1(HScontinuous), -2(phonon)'
    PRINT*,'                 100+(user-defined ising field), -100-(user-defined continuous field)'
    PRINT*,'  - n_checkerboard counts the number of independent boson fields belong to this kind of interaction'
    ivec(1:2)=(/1,1/)
    CALL saferead(2,ivec); n_checkerboard=ivec(2)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),':type_field, n_checkerboard'
    
    PRINT*,'g_field for i_g=',i
    PRINT*,"  - g_field is defined as H = g(1)/2 * ( c' * fmat *c - g(2) )^2"
    PRINT*,"  - for phonon, it is defined as H = eta * (c' * fmat * c - g(2) )^2"
    PRINT*,"    and g(3)=Debye frequency, g(1)=eta^2/M/g(3)^2"
    rvec(1:3)=(/4d0,1d0,0d0/)
    CALL saferead(3,rvec)
    WRITE(10,'(10x,3f10.6,10x,a)') rvec(1:3),':g_field'
    
    PRINT*,'dphi and dphi_global for i_g=',i
    PRINT*,'  - they are the maximal shift of the continuous boson field in local and global updates, respectively'
    rvec(1:4)=(/1d0,0d0,1d0,0d0/)
    CALL saferead(4,rvec)
    WRITE(10,'(4f10.6,10x,a)') rvec(1:4),':dphi, dphi_global'
    
    PRINT*,'ninterval_global and global_method for i_g=',i
    PRINT*,'  - do global update every ninterval_global sweeps'
    PRINT*,'  - global method = 0(totally random init), 1(totally random shift), 2(random shift at each time), 100(user-defined)'
    ivec(1:2)=(/100000,0/)
    CALL saferead(2,ivec)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),':ninterval_global, global_method'
    
    PRINT*,'ndim_field: dimension of each boson field'
    i1=1
    CALL saferead(i1); d=i1
    WRITE(10,'(1i40,10x,a)') i1,':ndim_field'
    
    DO flv=1,nflv
      DO k=1,d
        PRINT*,'fmat(k,:,flv): k=',k,' flv=',flv
        rvec(1:2*d)=0d0; rvec(2*d-1)=1d0
        CALL saferead(2*d,rvec)
        WRITE(10,'(1000f12.6)',advance='no') rvec(1:2*d)
        IF(k==1)THEN
          WRITE(10,'(10x,a)') ': fmat'
        ELSE
          WRITE(10,*)
        END IF
      END DO
    END DO


    DO k=1,n_checkerboard

      PRINT*,'for checkboard',k
    
      DO j=1,d
        PRINT*,'da,db,dc,orb for basis',j
        ivec(1:4)=(/0,0,0,1/)
        CALL saferead(4,ivec)
        WRITE(10,'(4i10,10x,a,1i4)') ivec(1:4),':da,db,dc,orb for basis',j
      END DO

      PRINT*,'how many conditions to define this checkerboard?'
      i1=1
      CALL saferead(i1); ncond=i1
      WRITE(10,'(1i40,10x,a)') i1,':n_cond'
      PRINT*,'for each condition, provide: ma,moda,mb,modb,mc,modc'
      PRINT*,'where the boson field lives only when mod(a,ma)==moda, mod(b,mb)==modb, mod(c,mc)==modc'
      DO j=1,ncond
        ivec(1:6)=(/1,0,1,0,1,0/)
        CALL saferead(6,ivec)
        WRITE(10,'(6i10)') ivec(1:6)
      END DO

    END DO

  END DO

  WRITE(10,*) '#------------------------------------------------'
  WRITE(10,*) '#           measurement block          '
  WRITE(10,*) '#------------------------------------------------'

  PRINT*,'nk_meas: how many k-points to measure'
  i1=0
  CALL saferead(i1); k=i1
  WRITE(10,'(1i40,10x,a)') i1,':nk_meas'
  
  PRINT*,'input the k_array defined by (a,b,c) below: (ka,kb,kc)=(a,b,c)*2*pi/(La,Lb,Lc)'
  DO i=1,k
    ivec(1:3)=(/0,0,0/)
    CALL saferead(3,ivec)
    WRITE(10,'(10x,3i10,10x)') ivec(1:3)
  END DO
  
  PRINT*,'nr_meas: how many (relative) r-points to measure'
  i1=0
  CALL saferead(i1); k=i1
  WRITE(10,'(1i40,10x,a)') i1,':nr_meas'
  
  PRINT*,"input the r_array defined by (a,b,c) below: G(a,b,c)=<O(x+a,y+b,z+c)O'(x,y,z)>"
  DO i=1,k
    ivec(1:3)=(/0,0,0/)
    CALL saferead(3,ivec)
    WRITE(10,'(10x,3i10,10x)') ivec(1:3)
  END DO

  PRINT*,'nrr_meas: how many rr-points to measure (useful in models without translation symmetry)'
  i1=0
  CALL saferead(i1); k=i1
  WRITE(10,'(1i40,10x,a)') i1,':nrr_meas'
  
  PRINT*,"input the rr_array defined by (a,b,c,a2,b2,c2) below: G(a,b,c,a2,b2,c2)=<O(a,b,c)O'(a2,b2,c2)>"
  DO i=1,k
    ivec(1:6)=(/0,0,0,0,0,0/)
    CALL saferead(6,ivec)
    WRITE(10,'(6i10)') ivec(1:6)
  END DO

  PRINT*,'ntau_meas: time displaced Green functions are measured during [-ntau_meas,ntau_meas]'
  i1=0
  CALL saferead(i1)
  WRITE(10,'(1i40,10x,a)') i1,':ntau_meas'

  PRINT*,'n_ph_meas, max_ndim_ph_meas: how many PH-channel two-particle Green functions to measure'
  PRINT*,'                             and the maximal dimension of fmat given below'
  ivec(1:2)=(/0,0/)
  CALL saferead(2,ivec); n2=ivec(1)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),':n_ph_meas,max_ndim_ph_meas'

  DO i=1,n2
    
    PRINT*,'hartree_ph_meas, fork_ph_meas for PH-', i
    lvec(1:2)=(/.true.,.true./)
    CALL saferead(2,lvec)
    WRITE(10,'(20x,2l10,10x,a)') lvec(1:2),':hartree_ph_meas,fork_ph_meas'

    PRINT*,'ndim_ph_meas and name_ph_meas for PH-', i
    i1=nflv; str='default_ph'
    CALL saferead(i1,str); d=i1
    WRITE(10,'(10x,1i10,1a20,10x,a)') i1,trim(adjustl(str)),':ndim_ph_meas, name_ph_meas'
    
    PRINT*,'for each basis, provide: da,db,dc,orb,flv'
    WRITE(10,*),'! for each basis, da,db,dc,orb,flv are given below'
    DO k=1,d
      ivec(1:5)=(/0,0,0,1,1/)
      CALL saferead(5,ivec)
      WRITE(10,'(5i10)') ivec(1:5)
    END DO

    PRINT*,'fmat is given below:'
    DO k=1,d
      PRINT*,'fmat(k,:) with k=',k
      rvec(1:2*d)=0d0; rvec(2*k-1)=1d0
      CALL saferead(2*d,rvec)
      WRITE(10,'(1000f12.6)',advance='no') rvec(1:2*d)
      IF(k==1)THEN
        WRITE(10,'(10x,a)') ':fmat_ph_meas'
      ELSE
        WRITE(10,*)
      END IF
    END DO

  END DO
  
  PRINT*,'n_pp_meas, max_ndim_pp_meas: how many PP-channel two-particle Green functions to measure'
  PRINT*,'                             and the maximal dimension of fmat given below'
  ivec(1:2)=(/0,0/)
  CALL saferead(2,ivec); n2=ivec(1)
  WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),':n_pp_meas,max_ndim_pp_meas'

  DO i=1,n2
    
    PRINT*,'fork13_pp_meas, fork14_pp_meas for PP-', i
    lvec(1:2)=(/.true.,.true./)
    CALL saferead(2,lvec)
    WRITE(10,'(20x,2l10,10x,a)') lvec(1:2),':fork13_pp_meas,fork14_pp_meas'

    PRINT*,'ndim_pp_meas and name_pp_meas for PP-', i
    i1=nflv; str='default_pp'
    CALL saferead(i1,str); d=i1
    WRITE(10,'(10x,1i10,1a20,10x,a)') i1,trim(adjustl(str)),':ndim_pp_meas, name_pp_meas'
    
    PRINT*,'for each basis, provide: da,db,dc,orb,flv'
    WRITE(10,*),'! for each basis, da,db,dc,orb,flv are given below'
    DO k=1,d
      ivec(1:5)=(/0,0,0,1,1/)
      CALL saferead(5,ivec)
      WRITE(10,'(5i10)') ivec(1:5)
    END DO

    PRINT*,'fmat is given below:'
    DO k=1,d
      PRINT*,'fmat(k,:) with k=',k
      rvec(1:2*d)=0d0; rvec(2*k-1)=1d0
      CALL saferead(2*d,rvec)
      WRITE(10,'(1000f12.6)',advance='no') rvec(1:2*d)
      IF(k==1)THEN
        WRITE(10,'(10x,a)') ':fmat_pp_meas'
      ELSE
        WRITE(10,*)
      END IF
    END DO

  END DO

  PRINT*,"ncross_ph_meas: how many crossing-PH-channel measurements G(ph1,ph2)=<O(ph1)O'(ph2)>"
  i1=0
  CALL saferead(i1); n2=i1
  WRITE(10,'(20x,1i10,20x,a)') i1,':ncross_ph_meas'

  DO i=1,n2
    PRINT*,'name of the channel'
    str='default_cross_ph'
    CALL saferead(str)
    WRITE(10,'(20x,1a20,10x,a)') trim(adjustl(str)),' :name of the channel'
    PRINT*,'which two channels are correlated?'
    ivec(1:2)=(/1,2/)
    CALL saferead(2,ivec)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),':cross_ph_meas'
  END DO
  
  PRINT*,"ncross_pp_meas: how many crossing-PP-channel measurements G(pp1,pp2)=<O(pp1)O'(pp2)>"
  i1=0
  CALL saferead(i1); n2=i1
  WRITE(10,'(20x,1i10,20x,a)') i1,':ncross_pp_meas'

  DO i=1,n2
    PRINT*,'name of the channel'
    str='default_cross_pp'
    CALL saferead(str)
    WRITE(10,'(20x,1a20,10x,a)') trim(adjustl(str)),':name of the channel'
    PRINT*,'which two channels are correlated?'
    ivec(1:2)=(/1,2/)
    CALL saferead(2,ivec)
    WRITE(10,'(20x,2i10,10x,a)') ivec(1:2),':cross_pp_meas'
  END DO


  PRINT*,'FAtech: whether to use Feldbacher-Assaad stablization algorithm for T=0'
  l1=.true.
  CALL saferead(l1)
  WRITE(10,'(1l40,10x,a)') l1,':FAtech'

  PRINT*,'do_measurement_external, n_meas_external: whether or not to do and how many external measurements'
  l1=.false.; i1=0
  CALL saferead(l1,i1)
  WRITE(10,'(20x,1l10,1i10,10x,a)') l1,i1,':do_measurement_external, n_meas_external'


  PRINT*,'do_tmpout_external: whether or not to call tmpout externally'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(30x,1l10,10x,a)') l1,':do_tmpout_external'


  PRINT*,'do_postprocess_external: whether or not to call postprocess externally'
  l1=.false.
  CALL saferead(l1)
  WRITE(10,'(30x,1l10,10x,a)') l1,':do_postprocess_external'

  CLOSE(10)

CONTAINS

SUBROUTINE saferead_logical(x)
  IMPLICIT NONE
  LOGICAL x
  CHARACTER(120) line
  INTEGER err
  PRINT'(1a,1l6)','press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT'(1a,100l6)','press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT'(1a,1i6)','press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT'(1a,100i6)','press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT'(1a,1f12.6)','press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT'(1a,100f12.6)','press ENTER for default value:',x
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT*,'press ENTER for default value:',x,y
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT*,'press ENTER for default value:',x,y
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT*,'press ENTER for default value:',y
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) y
      IF(err/=0)THEN
        PRINT*,'input again'
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
  PRINT*,'press ENTER for default value:',x,y
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))/=0)THEN
      READ(line,*,iostat=err) x,y
      IF(err/=0)THEN
        PRINT*,'input again'
        CYCLE
      END IF
    END IF
    EXIT
  END DO
END SUBROUTINE


END PROGRAM
