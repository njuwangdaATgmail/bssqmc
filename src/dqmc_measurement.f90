! definition of Green's functions:
!   G(k) = \sum_{r1,r2} < B(r1) B'(r2) > exp(ik(r1-r2)) / (La*Lb*Lc)**2
!   G(r) = \sum_{r2} < B(r2+r) B'(r2) > / (La*Lb*Lc)
!   G(r1,r2) = < B(r1) B'(r2) >
! where B = c, c'*fmat*c, c*fmat*c in three channels, respectively
SUBROUTINE measurement(time)
  USE dqmc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER i,j,k,flv,a,b,c,orb,a2,b2,c2,orb2,aa,bb,cc,da,db,dc,d,ir,p
  INTEGER i1,i2,i3,i4,k1,k2,k3,k4,flv1,flv2,flv3,flv4
  COMPLEX(8) kinetic,factor,g2(nsite,nsite,nflv), g3(nsite,nsite,nflv)
  COMPLEX(8), ALLOCATABLE :: ob4(:,:,:,:),obk(:),obr(:),obrr(:)
  
  ! get g2(i,j)=<c(i)c'(j)> and g3(i,j)=<c'(j)c(i)> in 2nd-order Trotter approximation
  ! NOTICE the positions of i and j. 
  ! Such a definition will benefit the calculation of time evolutions.
  DO flv=1,nflv
    g2(:,:,flv)=matmul(inv_expk_half(:,:,flv),matmul(g(:,:,flv),expk_half(:,:,flv)))
    g3(:,:,flv)=-g2(:,:,flv)
    DO i=1,nsite
      g3(i,i,flv)=g3(i,i,flv)+1d0
    END DO
  END DO

  ! kinetic energy
  kinetic=0d0
  DO i=1,nsite
    DO j=1,nsite
      kinetic=kinetic+sum(kmat(i,j,:)*g3(j,i,:))*currentphase/nsite
    END DO
  END DO
  CALL put_pool(kinetic)

  !-------------------------------------
  ! single-particle Green's functions
  !-------------------------------------
  
  IF(nk_meas>0)THEN

    ALLOCATE(ob4(nk_meas,norb,norb,nflv))
    ob4=0d0
    DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
      DO a2=1,La; DO b2=1,Lb; DO c2=1,Lc; DO orb2=1,norb; i2=label(a2,b2,c2,orb2)
        da=a-a2; db=b-b2; dc=c-c2
        DO flv=1,nflv
          ob4(:,orb,orb2,flv)=ob4(:,orb,orb2,flv)+g2(i,i2,flv)*expikr_array(:,da,db,dc)
        END DO
      END DO; END DO; END DO; END DO
    END DO; END DO; END DO; END DO
    ob4=ob4*currentphase/(La*Lb*Lc)**2
    CALL zput_array(nk_meas*norb*norb*nflv,ob4)
    DEALLOCATE(ob4)

  END IF

  IF(nr_meas>0)THEN

    ALLOCATE(ob4(nr_meas,norb,norb,nflv))
    ob4=0d0
    DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
      DO ir=1,nr_meas
        a2=mod(a-r_array(1,ir)-1+La,La)+1
        b2=mod(b-r_array(2,ir)-1+Lb,Lb)+1
        c2=mod(c-r_array(3,ir)-1+Lc,Lc)+1
        DO orb2=1,norb; i2=label(a2,b2,c2,orb2)
          ob4(ir,orb,orb2,:)=ob4(ir,orb,orb2,:)+g2(i,i2,:)
        END DO
      END DO
    END DO; END DO; END DO; END DO
    ob4=ob4*currentphase/(La*Lb*Lc)
    CALL zput_array(nr_meas*norb*norb*nflv,ob4)
    DEALLOCATE(ob4)
    
  END IF

  IF(nrr_meas>0)THEN
    
    ALLOCATE(ob4(nrr_meas,norb,norb,nflv))
    ob4=0d0
    DO ir=1,nrr_meas
      DO orb=1,norb; i=label(rr_array(1,1,ir),rr_array(2,1,ir),rr_array(3,1,ir),orb)
        DO orb2=1,norb; i2=label(rr_array(1,2,ir),rr_array(2,2,ir),rr_array(3,2,ir),orb2)
          ob4(ir,orb,orb2,:)=ob4(ir,orb,orb2,:)+g2(i,i2,:)
        END DO
      END DO
    END DO
    ob4=ob4*currentphase
    CALL zput_array(nrr_meas*norb*norb*nflv,ob4)
    DEALLOCATE(ob4)

  END IF

  !--------------------------------------------
  ! PH-channel two-particle Green's functions
  !--------------------------------------------

  DO k=1,n_ph_meas; d=ndim_ph_meas(k)

    IF(nk_meas>0)THEN
      ALLOCATE(obk(nk_meas))
      obk=0d0
    END IF
    IF(nr_meas>0)THEN
      ALLOCATE(obr(nr_meas))
      obr=0d0
    END IF
    IF(nrr_meas>0)THEN
      ALLOCATE(obrr(nrr_meas))
      obrr=0d0
    END IF

    DO k1=1,d; DO k2=1,d; DO k3=1,d; DO k4=1,d
      
      IF(abs(fmat_ph_meas(k1,k2,k))<1d-6)CYCLE
      IF(abs(fmat_ph_meas(k4,k3,k))<1d-6)CYCLE
      factor=fmat_ph_meas(k1,k2,k)*conjg(fmat_ph_meas(k4,k3,k))
      
      flv1=flv_ph_meas(k1,k)
      flv2=flv_ph_meas(k2,k)
      flv3=flv_ph_meas(k3,k)
      flv4=flv_ph_meas(k4,k)

      IF((.not.(flv1==flv2.and.flv3==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE
      
      IF(nk_meas>0)THEN

        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO aa=1,La; DO bb=1,Lb; DO cc=1,Lc
            da=a-aa; db=b-bb; dc=c-cc

            i1=nb_ph_meas(a,b,c,k1,k)
            i2=nb_ph_meas(a,b,c,k2,k)
            i3=nb_ph_meas(aa,bb,cc,k3,k)
            i4=nb_ph_meas(aa,bb,cc,k4,k)

            ! c'(i1,flv1) c(i2,flv2) c'(i3,flv3) c(i4,flv4)
            IF(flv1==flv2.and.flv3==flv4)THEN
              obk=obk+g3(i2,i1,flv1)*g3(i4,i3,flv3)*factor*expikr_array(:,da,db,dc)
            END IF
            IF(flv1==flv4.and.flv2==flv3)THEN
              obk=obk+g3(i4,i1,flv1)*g2(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
            END IF

          END DO; END DO; END DO
        END DO; END DO; END DO

      END IF

      IF(nr_meas>0)THEN
 
        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO ir=1,nr_meas
            aa=mod(a-r_array(1,ir)-1+La,La)+1
            bb=mod(b-r_array(2,ir)-1+Lb,Lb)+1
            cc=mod(c-r_array(3,ir)-1+Lc,Lc)+1
 
            i1=nb_ph_meas(a,b,c,k1,k)
            i2=nb_ph_meas(a,b,c,k2,k)
            i3=nb_ph_meas(aa,bb,cc,k3,k)
            i4=nb_ph_meas(aa,bb,cc,k4,k)
 
            ! c'(i1,flv1) c(i2,flv2) c'(i3,flv3) c(i4,flv4)
            IF(flv1==flv2.and.flv3==flv4)THEN
              obr(ir)=obr(ir)+g3(i2,i1,flv1)*g3(i4,i3,flv3)*factor
            END IF
            IF(flv1==flv4.and.flv2==flv3)THEN
              obr(ir)=obr(ir)+g3(i4,i1,flv1)*g2(i2,i3,flv2)*factor
            END IF
 
          END DO
        END DO; END DO; END DO
 
      END IF
 
      IF(nrr_meas>0)THEN
 
        DO ir=1,nrr_meas
          a=rr_array(1,1,ir)
          b=rr_array(2,1,ir)
          c=rr_array(3,1,ir)
          aa=rr_array(1,2,ir)
          bb=rr_array(2,2,ir)
          cc=rr_array(3,2,ir)
 
          i1=nb_ph_meas(a,b,c,k1,k)
          i2=nb_ph_meas(a,b,c,k2,k)
          i3=nb_ph_meas(aa,bb,cc,k3,k)
          i4=nb_ph_meas(aa,bb,cc,k4,k)

          ! c'(i1,flv1) c(i2,flv2) c'(i3,flv3) c(i4,flv4)
          IF(flv1==flv2.and.flv3==flv4)THEN
            obrr(ir)=obrr(ir)+g3(i2,i1,flv1)*g3(i4,i3,flv3)*factor
          END IF
          IF(flv1==flv4.and.flv2==flv3)THEN
            obrr(ir)=obrr(ir)+g3(i4,i1,flv1)*g2(i2,i3,flv2)*factor
          END IF
 
        END DO
 
      END IF

    END DO; END DO; END DO; END DO
    
    IF(nk_meas>0)THEN
      obk=obk*currentphase/(La*Lb*Lc)**2
      CALL zput_array(nk_meas,obk)
      DEALLOCATE(obk)
    END IF
    IF(nr_meas>0)THEN
      obr=obr*currentphase/(La*Lb*Lc)
      CALL zput_array(nr_meas,obr)
      DEALLOCATE(obr)
    END IF
    IF(nrr_meas>0)THEN
      obrr=obrr*currentphase
      CALL zput_array(nrr_meas,obrr)
      DEALLOCATE(obrr)
    END IF

  END DO

  !--------------------------------------------
  ! PP-channel two-particle Green's functions
  !--------------------------------------------

  DO k=1,n_pp_meas; d=ndim_pp_meas(k)

    IF(nk_meas>0)THEN
      ALLOCATE(obk(nk_meas))
      obk=0d0
    END IF
    IF(nr_meas>0)THEN
      ALLOCATE(obr(nr_meas))
      obr=0d0
    END IF
    IF(nrr_meas>0)THEN
      ALLOCATE(obrr(nrr_meas))
      obrr=0d0
    END IF

    DO k1=1,d; DO k2=1,d; DO k3=1,d; DO k4=1,d
      
      IF(abs(fmat_pp_meas(k1,k2,k))<1d-6)CYCLE
      IF(abs(fmat_pp_meas(k4,k3,k))<1d-6)CYCLE
      factor=fmat_pp_meas(k1,k2,k)*conjg(fmat_pp_meas(k4,k3,k))
      
      flv1=flv_pp_meas(k1,k)
      flv2=flv_pp_meas(k2,k)
      flv3=flv_pp_meas(k3,k)
      flv4=flv_pp_meas(k4,k)

      IF((.not.(flv1==flv3.and.flv2==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE
      
      IF(nk_meas>0)THEN

        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO aa=1,La; DO bb=1,Lb; DO cc=1,Lc
            da=a-aa; db=b-bb; dc=c-cc

            i1=nb_pp_meas(a,b,c,k1,k)
            i2=nb_pp_meas(a,b,c,k2,k)
            i3=nb_pp_meas(aa,bb,cc,k3,k)
            i4=nb_pp_meas(aa,bb,cc,k4,k)

            ! c(i1,flv1) c(i2,flv2) c'(i3,flv3) c'(i4,flv4)
            IF(flv1==flv3.and.flv2==flv4)THEN
              obk=obk-g2(i1,i3,flv1)*g2(i2,i4,flv2)*factor*expikr_array(:,da,db,dc)
            END IF
            IF(flv1==flv4.and.flv2==flv3)THEN
              obk=obk+g2(i1,i4,flv1)*g2(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
            END IF

          END DO; END DO; END DO
        END DO; END DO; END DO

      END IF

      IF(nr_meas>0)THEN
 
        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO ir=1,nr_meas
            aa=mod(a-r_array(1,ir)-1+La,La)+1
            bb=mod(b-r_array(2,ir)-1+Lb,Lb)+1
            cc=mod(c-r_array(3,ir)-1+Lc,Lc)+1
 
            i1=nb_pp_meas(a,b,c,k1,k)
            i2=nb_pp_meas(a,b,c,k2,k)
            i3=nb_pp_meas(aa,bb,cc,k3,k)
            i4=nb_pp_meas(aa,bb,cc,k4,k)
 
            ! c(i1,flv1) c(i2,flv2) c'(i3,flv3) c'(i4,flv4)
            IF(flv1==flv3.and.flv2==flv4)THEN
              obr(ir)=obr(ir)-g2(i1,i3,flv1)*g2(i2,i4,flv2)*factor
            END IF
            IF(flv1==flv4.and.flv2==flv3)THEN
              obr(ir)=obr(ir)+g2(i1,i4,flv1)*g2(i2,i3,flv2)*factor
            END IF
 
          END DO
        END DO; END DO; END DO
 
      END IF
 
      IF(nrr_meas>0)THEN
 
        DO ir=1,nrr_meas
          a=rr_array(1,1,ir)
          b=rr_array(2,1,ir)
          c=rr_array(3,1,ir)
          aa=rr_array(1,2,ir)
          bb=rr_array(2,2,ir)
          cc=rr_array(3,2,ir)
 
          i1=nb_pp_meas(a,b,c,k1,k)
          i2=nb_pp_meas(a,b,c,k2,k)
          i3=nb_pp_meas(aa,bb,cc,k3,k)
          i4=nb_pp_meas(aa,bb,cc,k4,k)
          
          ! c(i1,flv1) c(i2,flv2) c'(i3,flv3) c'(i4,flv4)
          IF(flv1==flv3.and.flv2==flv4)THEN
            obrr(ir)=obrr(ir)-g2(i1,i3,flv1)*g2(i2,i4,flv2)*factor
          END IF
          IF(flv1==flv4.and.flv2==flv3)THEN
            obrr(ir)=obrr(ir)+g2(i1,i4,flv1)*g2(i2,i3,flv2)*factor
          END IF
 
        END DO
 
      END IF

    END DO; END DO; END DO; END DO
    
    IF(nk_meas>0)THEN
      obk=obk*currentphase/(La*Lb*Lc)**2
      CALL zput_array(nk_meas,obk)
      DEALLOCATE(obk)
    END IF
    IF(nr_meas>0)THEN
      obr=obr*currentphase/(La*Lb*Lc)
      CALL zput_array(nr_meas,obr)
      DEALLOCATE(obr)
    END IF
    IF(nrr_meas>0)THEN
      obrr=obrr*currentphase
      CALL zput_array(nrr_meas,obrr)
      DEALLOCATE(obrr)
    END IF

  END DO


  !---------------------------------
  ! unequal-time Green's functions
  !---------------------------------

  IF(proj)THEN

    DO p=1,ntau_meas  ! evolve to the left and right, respectively

      IF(nk_meas>0)THEN
!        ...
      END IF

      IF(nr_meas>0)THEN
!        ...
      END IF

      IF(nrr_meas>0)THEN
!        ...
      END IF

    END DO

  ELSE

     DO p=1,ntau_meas  ! evolve to the left and right, respectively

      IF(nk_meas>0)THEN
 !       ...
      END IF

      IF(nr_meas>0)THEN
 !       ...
      END IF

      IF(nrr_meas>0)THEN
!       ...
      END IF

    END DO

  END IF

  ! measurement externally
  IF(do_measure_external) CALL measurement_external(time)

END SUBROUTINE
