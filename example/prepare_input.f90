! This is an interactive engine for the users to prepare the input file 'dqmc.in'.
! This program is standalone.
! to be added: 
! + more comments
! + add default values
! + when reading error occurs, read again
PROGRAM main
  IMPLICIT NONE
  INTEGER i1,i2,i3,i4,i5,i6,nflv,norb,i,j,k,flv,n_g,nfield,ncond,d,n_checkerboard,n2,err
  REAL(8) r1,r2,r3,r4
  COMPLEX(8) z1,z2,z3,z4
  LOGICAL l1,l2,l3,l4
  CHARACTER(120) str,line
  INTEGER, ALLOCATABLE :: ivec(:)
  REAL(8), ALLOCATABLE :: rvec(:)

  OPEN(10,file='dqmc.in')

  WRITE(10,*) '#----------------------------'
  WRITE(10,*) '# Monte Carlo control block'
  WRITE(10,*) '#----------------------------'
  
  PRINT*,'restart: restart mode or not (true or false)'
  PRINT*,'press enter for default (false)'
  DO
    READ(5,'(1a)') line
    IF(len(trim(line))==0)THEN
      l1=.false.
      EXIT
    ELSE
      READ(line,*,iostat=err) l1
      IF(err/=0)THEN
        PRINT*,'input again'
      ELSE
        EXIT
      END IF
    END IF
  END DO
  WRITE(10,'(1l40,10x,a)') l1,':restart'
    

  PRINT*,'proj: projector mode or not (true or false)'
  READ*,l1
  WRITE(10,'(1l40,10x,a)') l1,':proj'

  PRINT*,'nflv: number of flavors (integer)'
  READ*,i1; nflv=i1
  WRITE(10,'(1i40,10x,a)') i1,':nflv'

  PRINT*,'ncopy(nflv): number of copies for each flavor (integer)'
  ALLOCATE(ivec(nflv))
  READ*,ivec
  WRITE(10,'(1i40,10x,a)') ivec,':ncopy'
  DEALLOCATE(ivec)

  PRINT*,'beta: inverse of temperature or projecting time (real)'
  READ*,r1
  WRITE(10,'(1f40.6,10x,a)') r1,':beta'

  PRINT*,'dtau:'
  READ*,r1
  WRITE(10,'(1f40.6,10x,a)') r1,':dtau'

  PRINT*,'nsp:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':nsp'

  PRINT*,'nbin:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':nbin'

  PRINT*,'nwarmup:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':nwarmup'

  PRINT*,'nmeasure:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':nmeasure'

  PRINT*,'ninterval:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':ninterval'

  PRINT*,'ntmpout:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':ntmpout'

  PRINT*,'nscratch:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':nscratch'
  
  PRINT*,'ngroup:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':ngroup'

  PRINT*,'randomseed:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':randomseed'

  PRINT*,'newMetro:'
  READ*,r1
  WRITE(10,'(1f40.6,10x,a)') r1,':newMetro'

  WRITE(10,*) '#----------------------------'
  WRITE(10,*) '# Monte Carlo control block'
  WRITE(10,*) '#----------------------------'
  
  PRINT*,'norb:'
  READ*,i1; norb=i1
  WRITE(10,'(1i40,10x,a)') i1,':norb'

  PRINT*,'La,Lb,Lc:'
  READ*,i1,i2,i3
  WRITE(10,'(10x,3i10,10x,a)') i1,i2,i3,':La,Lb,Lc'

  PRINT*,'nelec:'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':nelec'

  PRINT*,'pbca,pbcb,pbcc:'
  READ*,l1,l2,l3
  WRITE(10,'(10x,3l10,10x,a)') l1,l2,l3,':pbca,pbcb,pbcc'

  PRINT*,'twista,twistb,twistc: in unit of 2pi'
  READ*,r1,r2,r3
  WRITE(10,'(10x,3f10.6,10x,a)') r1,r2,r3,':twista,twistb,twistc in unit of 2pi'

  PRINT*,'a0r: unit cell axis'
  READ*,r1,r2,r3
  WRITE(10,'(10x,3f10.6,10x,a)') r1,r2,r3,':a0r'

  PRINT*,'b0r: unit cell axis'
  READ*,r1,r2,r3
  WRITE(10,'(10x,3f10.6,10x,a)') r1,r2,r3,':b0r'
  
  PRINT*,'c0r: unit cell axis'
  READ*,r1,r2,r3
  WRITE(10,'(10x,3f10.6,10x,a)') r1,r2,r3,':c0r'
  
  DO i=1,norb
    PRINT*,'position of the orbital',i
    READ*,r1,r2,r3
    WRITE(10,'(10x,3f10.6,10x,a,1i2)') r1,r2,r3,':rorb-',i
  END DO

  PRINT*,'cuta,cutb,cutc:'
  READ*,i1,i2,i3
  WRITE(10,'(10x,3i10,10x,a)') i1,i2,i3,':cuta,cutb,cutc'

  PRINT*,'nhop:'
  READ*,i6
  WRITE(10,'(1i40,10x,a)') i6,':nhop'

  DO i=1,i6
    PRINT*,'da,db,dc,orb1,orb2,Re(t),Im(t):'
    READ*,i1,i2,i3,i4,i5,r1,r2
    WRITE(10,'(5i4,2f10.6,10x,a)') i1,i2,i3,i4,i5,r1,r2,':da,db,dc,orb1,orb2,Re(t),Im(t)'
  END DO

  WRITE(10,*) '#----------------------------'
  WRITE(10,*) '# boson field block'
  WRITE(10,*) '#----------------------------'

  PRINT*,'n_g,nfield,max_ndim_field,max_isingmax:'
  READ*,i1,i2,i3,i4; n_g=i1; nfield=i2
  WRITE(10,'(8x,4i8,10x,a)') i1,i2,i3,i4,':n_g,nfield,max_ndim_field,max_isingmax'

  DO i=1,n_g
    
    PRINT*,'type_field, n_checkerboard for i_g=',i
    READ*,i1,i2; n_checkerboard=i2
    WRITE(10,'(20x,2i10,10x,a)') i1,i2,':type_field, n_checkerboard'
    
    PRINT*,'g_field for i_g=',i
    READ*,r1,r2,r3
    WRITE(10,'(10x,3f10.6,10x,a)') r1,r2,r3,':g_field'
    
    PRINT*,'dphi and dphi_global for i_g=',i
    READ*,r1,r2,r3,r4
    WRITE(10,'(4f10.6,10x,a)') r1,r2,r3,r4,':dphi, dphi_global'
    
    PRINT*,'ninterval_global and global_method for i_g=',i
    READ*,i1,i2
    WRITE(10,'(20x,2i10,10x,a)') i1,i2,':ninterval_global, global_method'
    
    PRINT*,'ndim_field'
    READ*,i1; d=i1
    WRITE(10,'(1i40,10x,a)') i1,':ndim_field'
    
    ALLOCATE(rvec(2*d))
    DO flv=1,nflv
      DO k=1,d
        PRINT*,'fmat(k,:,flv): k=',k,' flv=',flv
        READ*,rvec(:)
        WRITE(10,'(1000f12.6)') rvec
      END DO
    END DO
    DEALLOCATE(rvec)

    DO k=1,d
      PRINT*,'da,db,dc,orb for basis',k
      READ*,i1,i2,i3,i4
      WRITE(10,'(4i10,10x,a)') i1,i2,i3,i4
    END DO

    DO k=1,n_checkerboard
      PRINT*,'how many conditions to define this checkerboard?'
      READ*,i1; ncond=i1
      WRITE(10,'(1i40,10x,a)') i1,':n_cond'
      PRINT*,'for each condition, provide: ma,moda,mb,modb,mc,modc'
      DO j=1,ncond
        READ*,i1,i2,i3,i4,i5,i6
        WRITE(10,'(6i10)') i1,i2,i3,i4,i5,i6
      END DO
    END DO

  END DO

  WRITE(10,*) '#----------------------------'
  WRITE(10,*) '# measurement block'
  WRITE(10,*) '#----------------------------'

  PRINT*,'nk_meas'
  READ*,i1; k=i1
  WRITE(10,'(1i40,10x,a)') i1,':nk_meas'

  PRINT*,'input the k_array below'
  DO i=1,k
    READ*,i1,i2,i3
    WRITE(10,'(10x,3i10,10x)') i1,i2,i3
  END DO

  PRINT*,'nr_meas'
  READ*,i1; k=i1
  WRITE(10,'(1i40,10x,a)') i1,':nr_meas'
  
  PRINT*,'input the r_array below'
  DO i=1,k
    READ*,i1,i2,i3
    WRITE(10,'(10x,3i10,10x)') i1,i2,i3
  END DO


  PRINT*,'nrr_meas'
  READ*,i1; k=i1
  WRITE(10,'(1i40,10x,a)') i1,':nrr_meas'
  
  PRINT*,'input the r_array below'
  DO i=1,k
    READ*,i1,i2,i3,i4,i5,i6
    WRITE(10,'(6i10)') i1,i2,i3,i4,i5,i6
  END DO

  PRINT*,'ntau_meas'
  READ*,i1
  WRITE(10,'(1i40,10x,a)') i1,':ntau_meas'

  PRINT*,'n_ph_meas, max_ndim_ph_meas'
  READ*,i1,i2; n2=i1
  WRITE(10,'(20x,2i10,10x,a)') i1,i2,':n_ph_meas,max_ndim_ph_meas'

  DO i=1,n2
    
    PRINT*,'hartree_ph_meas, fork_ph_meas for PH-', i
    READ*,l1,l2
    WRITE(10,'(20x,2l10,10x,a)') l1,l2,':hartree_ph_meas,fork_ph_meas'

    PRINT*,'ndim_ph_meas and name_ph_meas for PH-', i
    READ*,i1,str; d=i1
    WRITE(10,'(10x,1i10,1a20,10x,a)') i1,trim(adjustl(str)),':ndim_ph_meas, name_ph_meas'
    
    PRINT*,'for each basis, provide: da,db,dc,orb,flv'
    WRITE(10,*),'! for each basis, da,db,dc,orb,flv are given below'
    DO k=1,d
      READ*, i1,i2,i3,i4,i5
      WRITE(10,'(5i10)') i1,i2,i3,i4,i5
    END DO

    ALLOCATE(rvec(2*d))
    PRINT*,'fmat is given below:'
    DO k=1,d
      READ*,rvec(:)
      WRITE(10,'(1000f12.6)') rvec
    END DO
    DEALLOCATE(rvec)

  END DO
  
  PRINT*,'n_pp_meas, max_ndim_pp_meas'
  READ*,i1,i2; n2=i1
  WRITE(10,'(20x,2i10,10x,a)') i1,i2,':n_pp_meas,max_ndim_pp_meas'

  DO i=1,n2
    
    PRINT*,'fork13_pp_meas, fork14_pp_meas for PP-', i
    READ*,l1,l2
    WRITE(10,'(20x,2l10,10x,a)') l1,l2,':fork13_pp_meas,fork14_pp_meas'

    PRINT*,'ndim_pp_meas and name_pp_meas for PP-', i
    READ*,i1,str; d=i1
    WRITE(10,'(20x,1i10,1a20,a)') i1,trim(adjustl(str)),':ndim_pp_meas, name_pp_meas'
    
    PRINT*,'for each basis, provide: da,db,dc,orb,flv'
    WRITE(10,*),'! for each basis, da,db,dc,orb,flv are given below'
    DO k=1,d
      READ*, i1,i2,i3,i4,i5
      WRITE(10,'(5i10)') i1,i2,i3,i4,i5
    END DO

    ALLOCATE(rvec(2*d))
    PRINT*,'fmat is given below:'
    DO k=1,d
      READ*,rvec(:)
      WRITE(10,'(1000f12.6)') rvec
    END DO
    DEALLOCATE(rvec)

  END DO


  PRINT*,'ncross_ph_meas'
  READ*,i1; n2=i1
  WRITE(10,'(20x,1i10,20x,a)') i1,':ncross_ph_meas'

  DO i=1,n2
    PRINT*,'name of the channel'
    READ*,str
    WRITE(10,'(20x,1l10,20x,a)') trim(adjustl(str)),' :name of the channel'
    PRINT*,'which two channels are correlated?'
    READ*,i1,i2
    WRITE(10,'(20x,2i10,10x,a)') i1,i2,':cross_ph_meas'
  END DO
  

  PRINT*,'ncross_pp_meas'
  READ*,i1; n2=i1
  WRITE(10,'(20x,1i10,20x,a)') i1,':ncross_pp_meas'

  DO i=1,n2
    PRINT*,'name of the channel'
    READ*,str
    WRITE(10,'(20x,1l10,20x,a)') trim(adjustl(str)),':name of the channel'
    PRINT*,'which two channels are correlated?'
    READ*,i1,i2
    WRITE(10,'(20x,2i10,10x,a)') i1,i2,':cross_pp_meas'
  END DO


  PRINT*,'FAtech:'
  READ*,l1
  WRITE(10,'(1l40,10x,a)') l1,':FAtech'

  PRINT*,'do_measurement_external, n_meas_external'
  READ*,l1,i1
  WRITE(10,'(20x,1l10,1i10,10x,a)') l1,i1,':do_measurement_external, n_meas_external'


  PRINT*,'do_tmpout_external'
  READ*,l1
  WRITE(10,'(30x,1l10,10x,a)') l1,':do_tmpout_external'


  PRINT*,'do_postprocess_external'
  READ*,l1
  WRITE(10,'(30x,1l10,10x,a)') l1,':do_postprocess_external'

  CLOSE(10)

END PROGRAM
