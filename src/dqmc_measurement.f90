! definition of Green's functions:
!   G(k) = \sum_{r1,r2} < B(r1) B'(r2) > exp(ik(r1-r2)) / (La*Lb*Lc)**2
!   G(r) = \sum_{r2} < B(r2+r) B'(r2) > / (La*Lb*Lc)
!   G(r1,r2) = < B(r1) B'(r2) >
! where B = c, c'*fmat*c, c*fmat*c in three channels, respectively
SUBROUTINE measurement(time)
  USE dqmc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER i,j,k,flv,a,b,c,orb,a2,b2,c2,orb2,aa,bb,cc,da,db,dc,d,ir,p,tau,kk,dd,ic
  INTEGER i1,i2,i3,i4,k1,k2,k3,k4,flv1,flv2,flv3,flv4
  COMPLEX(8) kinetic,factor,g2(nsite,nsite,nflv),g3(nsite,nsite,nflv)
  COMPLEX(8), ALLOCATABLE :: ob4(:,:,:,:),obk(:),obr(:),obrr(:)
  COMPLEX(8), ALLOCATABLE :: ob4_(:,:,:,:),obk_(:),obr_(:),obrr_(:)
  COMPLEX(8), ALLOCATABLE :: g_save(:,:,:),g3old(:,:,:),gt2(:,:,:),gt3(:,:,:)
  COMPLEX(8), ALLOCATABLE :: qmat2(:,:,:),qmat3(:,:,:),dvec2(:,:),dvec3(:,:),rmat2(:,:,:),lmat3(:,:,:)

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
  
  !----------------------------------------------------
  ! crossing-PH-channel two-particle Green's functions
  !----------------------------------------------------
  
  DO ic=1,ncross_ph_meas
    
    k=cross_ph_meas(1,ic)
    kk=cross_ph_meas(2,ic)
    d=ndim_ph_meas(k)
    dd=ndim_ph_meas(kk)

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

    DO k1=1,d; DO k2=1,d; DO k3=1,dd; DO k4=1,dd

      IF(abs(fmat_ph_meas(k1,k2,k))<1d-6)CYCLE
      IF(abs(fmat_ph_meas(k4,k3,kk))<1d-6)CYCLE
      factor=fmat_ph_meas(k1,k2,k)*conjg(fmat_ph_meas(k4,k3,kk))

      flv1=flv_ph_meas(k1,k)
      flv2=flv_ph_meas(k2,k)
      flv3=flv_ph_meas(k3,kk)
      flv4=flv_ph_meas(k4,kk)

      IF((.not.(flv1==flv2.and.flv3==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE

      IF(nk_meas>0)THEN

        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO aa=1,La; DO bb=1,Lb; DO cc=1,Lc
            da=a-aa; db=b-bb; dc=c-cc

            i1=nb_ph_meas(a,b,c,k1,k)
            i2=nb_ph_meas(a,b,c,k2,k)
            i3=nb_ph_meas(aa,bb,cc,k3,kk)
            i4=nb_ph_meas(aa,bb,cc,k4,kk)

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
            i3=nb_ph_meas(aa,bb,cc,k3,kk)
            i4=nb_ph_meas(aa,bb,cc,k4,kk)

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
          i3=nb_ph_meas(aa,bb,cc,k3,kk)
          i4=nb_ph_meas(aa,bb,cc,k4,kk)

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
  ! crossing-PP-channel two-particle Green's functions
  !--------------------------------------------

  DO ic=1,ncross_pp_meas

    k=cross_pp_meas(1,ic)
    kk=cross_pp_meas(2,ic)
    d=ndim_pp_meas(k)
    dd=ndim_pp_meas(kk)

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

    DO k1=1,d; DO k2=1,d; DO k3=1,dd; DO k4=1,dd

      IF(abs(fmat_pp_meas(k1,k2,k))<1d-6)CYCLE
      IF(abs(fmat_pp_meas(k4,k3,kk))<1d-6)CYCLE
      factor=fmat_pp_meas(k1,k2,k)*conjg(fmat_pp_meas(k4,k3,kk))

      flv1=flv_pp_meas(k1,k)
      flv2=flv_pp_meas(k2,k)
      flv3=flv_pp_meas(k3,kk)
      flv4=flv_pp_meas(k4,kk)

      IF((.not.(flv1==flv3.and.flv2==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE

      IF(nk_meas>0)THEN

        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO aa=1,La; DO bb=1,Lb; DO cc=1,Lc
            da=a-aa; db=b-bb; dc=c-cc

            i1=nb_pp_meas(a,b,c,k1,k)
            i2=nb_pp_meas(a,b,c,k2,k)
            i3=nb_pp_meas(aa,bb,cc,k3,kk)
            i4=nb_pp_meas(aa,bb,cc,k4,kk)

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
            i3=nb_pp_meas(aa,bb,cc,k3,kk)
            i4=nb_pp_meas(aa,bb,cc,k4,kk)

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
          i3=nb_pp_meas(aa,bb,cc,k3,kk)
          i4=nb_pp_meas(aa,bb,cc,k4,kk)

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
  !
  ! when tau>0:
  !   < c(t) c'(0) > = B(t,0) < c(0) c'(0) >
  !   < c'(t) c(0) > = < c'(0) c(0) > B^{-1}(t,0)
  !   < c(t) c'(t) > = B(t,0) < c(0) c'(0) > B^{-1}(t,0)
  !   < c'(t) c(t) > = B(t,0) < c'(0) c(0) > B^{-1}(t,0) = 1 - < c(t) c'(t) >
  !
  ! For T=0, the following properties can be used to stablize the calculations:
  !   < c(t) c'(t) >^2 = < c(t) c'(t) >
  !   < c'(t) c(t) >^2 = < c'(t) c(t) >
  ! see Feldbacher & Assaad, PRB 63, 073105 (2001)
  IF(ntau_meas>0)THEN

    ALLOCATE(g_save(nsite,nsite,nflv),g3old(nsite,nsite,nflv),gt2(nsite,nsite,nflv),gt3(nsite,nsite,nflv))
    g_save=g
    g3old=g3

    IF(proj.and.FAtech)THEN
      gt2=g2
      gt3=g3
    ELSE
      ! perform QDR and LDQ decomposition
      ! g2 = qmat2 * dvec2 * rmat2
      ! g3 = lmat3 * dvec3 * qmat3
      ALLOCATE(qmat2(nsite,nsite,nflv),qmat3(nsite,nsite,nflv),dvec2(nsite,nflv),dvec3(nsite,nflv), &
        &      rmat2(nsite,nsite,nflv),lmat3(nsite,nsite,nflv))
      DO flv=1,nflv
        qmat2(:,:,flv)=g2(:,:,flv)
        CALL zqdr(nsite,nsite,qmat2(:,:,flv),rmat2(:,:,flv),dvec2(:,flv))
        qmat3(:,:,flv)=g3(:,:,flv)
        CALL zldq(nsite,nsite,qmat3(:,:,flv),lmat3(:,:,flv),dvec3(:,flv))
      END DO
    END IF

    DO tau=1,ntau_meas

      p=time+tau-1
      IF(p>ntime) p=p-ntime

      ! do time evolutions
      ! we need compare FAtech=.true. or .false. then decide which one is used.
      IF(proj.and.FAtech)THEN

        DO flv=1,nflv
          CALL evolve_left_2nd(p,g2(:,:,flv),nsite,flv,.false.)
          gt2(:,:,flv)=matmul(g2(:,:,flv),gt2(:,:,flv))
          CALL evolve_right_2nd(p,g3(:,:,flv),nsite,flv,.true.)
          gt3(:,:,flv)=matmul(gt3(:,:,flv),g3(:,:,flv))
        END DO

        IF(mod(tau,ngroup)==0)THEN
          CALL update_scratch_T0_(p+1)
          DO flv=1,nflv
            g2(:,:,flv)=matmul(inv_expk_half(:,:,flv),matmul(g(:,:,flv),expk_half(:,:,flv)))
            g3(:,:,flv)=-g2(:,:,flv)
            DO i=1,nsite
              g3(i,i,flv)=g3(i,i,flv)+1d0
            END DO
          END DO
        ELSE
          DO flv=1,nflv
            CALL evolve_right_2nd(p,g2(:,:,flv),nsite,flv,.true.)
            !CALL evolve_left_2nd(p,g3(:,:,flv),nsite,flv,.false.)
            g3(:,:,flv)=-g2(:,:,flv)
            DO i=1,nsite
              g3(i,i,flv)=g3(i,i,flv)+1d0
            END DO
          END DO
        END IF

      ELSE

        ! get < c(t) c'(0) > and < c'(t) c(0) > stored in gt2 and gt3
        DO flv=1,nflv
          CALL evolve_left_2nd(p,qmat2(:,:,flv),nsite,flv,.false.)
          CALL evolve_right_2nd(p,qmat3(:,:,flv),nsite,flv,.true.)
          DO i=1,nsite
            gt2(:,i,flv)=qmat2(:,i,flv)*dvec2(i,flv)
            gt3(:,i,flv)=dvec3(:,flv)*qmat3(:,i,flv)
          END DO
          qmat2(:,:,flv)=matmul(gt2(:,:,flv),rmat2(:,:,flv))
          qmat3(:,:,flv)=matmul(lmat3(:,:,flv),gt3(:,:,flv))
          gt2(:,:,flv)=qmat2(:,:,flv)
          gt3(:,:,flv)=qmat3(:,:,flv)
          ! since we have already got gt2 and gt3, we perform QDR and LDQ without costing too much time
          CALL zqdr(nsite,nsite,qmat2(:,:,flv),rmat2(:,:,flv),dvec2(:,flv))
          CALL zldq(nsite,nsite,qmat3(:,:,flv),lmat3(:,:,flv),dvec3(:,flv))
        END DO

        ! get < c'(t) c(t) > stored in g3
        IF(mod(tau,ngroup)==0)THEN
          IF(proj)THEN
            CALL update_scratch_T0_(p+1)
          ELSE
            CALL update_scratch_(mod(p,ntime)+1)  ! NOTICE: DO NOT change Bstring!
          END IF
          DO flv=1,nflv
            g3(:,:,flv)=-matmul(inv_expk_half(:,:,flv),matmul(g(:,:,flv),expk_half(:,:,flv)))
            DO i=1,nsite
              g3(i,i,flv)=g3(i,i,flv)+1d0
            END DO
          END DO
        ELSE
          DO flv=1,nflv
            CALL evolve_left_2nd(p,g3(:,:,flv),nsite,flv,.false.)
            CALL evolve_right_2nd(p,g3(:,:,flv),nsite,flv,.true.)
          END DO
        END IF

      END IF

      ! do measurement using g3, g3old, gt2, gt3

      !-------------------------------------
      ! single-particle Green's functions
      !-------------------------------------

      IF(nk_meas>0)THEN

        ALLOCATE(ob4(nk_meas,norb,norb,nflv))
        ALLOCATE(ob4_(nk_meas,norb,norb,nflv))
        ob4=0d0
        ob4_=0d0
        DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
          DO a2=1,La; DO b2=1,Lb; DO c2=1,Lc; DO orb2=1,norb; i2=label(a2,b2,c2,orb2)
            da=a-a2; db=b-b2; dc=c-c2
            DO flv=1,nflv
              ! < c(i,t) c'(i2,0) >
              ob4(:,orb,orb2,flv)=ob4(:,orb,orb2,flv)+gt2(i,i2,flv)*expikr_array(:,da,db,dc)
              ! - < c'(i2,t) c(i,0) >
              ob4_(:,orb,orb2,flv)=ob4_(:,orb,orb2,flv)-gt3(i,i2,flv)*expikr_array(:,da,db,dc)
            END DO
          END DO; END DO; END DO; END DO
        END DO; END DO; END DO; END DO
        ob4=ob4*currentphase/(La*Lb*Lc)**2
        ob4_=ob4_*currentphase/(La*Lb*Lc)**2
        CALL zput_array(nk_meas*norb*norb*nflv,ob4)
        CALL zput_array(nk_meas*norb*norb*nflv,ob4_)
        DEALLOCATE(ob4,ob4_)

      END IF

      IF(nr_meas>0)THEN

        ALLOCATE(ob4(nr_meas,norb,norb,nflv))
        ALLOCATE(ob4_(nr_meas,norb,norb,nflv))
        ob4=0d0
        ob4_=0d0
        DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
          DO ir=1,nr_meas
            a2=mod(a-r_array(1,ir)-1+La,La)+1
            b2=mod(b-r_array(2,ir)-1+Lb,Lb)+1
            c2=mod(c-r_array(3,ir)-1+Lc,Lc)+1
            DO orb2=1,norb; i2=label(a2,b2,c2,orb2)
              ! < c(i,t) c'(i2,0) >
              ob4(ir,orb,orb2,:)=ob4(ir,orb,orb2,:)+gt2(i,i2,:)
              ! - < c'(i2,t) c(i,0) >
              ob4_(ir,orb,orb2,:)=ob4_(ir,orb,orb2,:)-gt3(i,i2,:)
            END DO
          END DO
        END DO; END DO; END DO; END DO
        ob4=ob4*currentphase/(La*Lb*Lc)
        ob4_=ob4_*currentphase/(La*Lb*Lc)
        CALL zput_array(nr_meas*norb*norb*nflv,ob4)
        CALL zput_array(nr_meas*norb*norb*nflv,ob4_)
        DEALLOCATE(ob4,ob4_)

      END IF

      IF(nrr_meas>0)THEN

        ALLOCATE(ob4(nrr_meas,norb,norb,nflv))
        ob4=0d0
        DO ir=1,nrr_meas
          DO orb=1,norb; i=label(rr_array(1,1,ir),rr_array(2,1,ir),rr_array(3,1,ir),orb)
            DO orb2=1,norb; i2=label(rr_array(1,2,ir),rr_array(2,2,ir),rr_array(3,2,ir),orb2)
              ! < c(i,t) c'(i2,0) >
              ob4(ir,orb,orb2,:)=ob4(ir,orb,orb2,:)+gt2(i,i2,:)
              ! - < c'(i2,t) c(i,0) >
              ob4_(ir,orb,orb2,:)=ob4_(ir,orb,orb2,:)-gt3(i,i2,:)
            END DO
          END DO
        END DO
        ob4=ob4*currentphase
        ob4_=ob4_*currentphase
        CALL zput_array(nr_meas*norb*norb*nflv,ob4)
        CALL zput_array(nr_meas*norb*norb*nflv,ob4_)
        DEALLOCATE(ob4,ob4_)

      END IF

      !--------------------------------------------
      ! PH-channel two-particle Green's functions
      !--------------------------------------------

      DO k=1,n_ph_meas; d=ndim_ph_meas(k)

        IF(nk_meas>0)THEN
          ALLOCATE(obk(nk_meas),obk_(nk_meas))
          obk=0d0; obk_=0d0
        END IF
        IF(nr_meas>0)THEN
          ALLOCATE(obr(nr_meas),obr_(nr_meas))
          obr=0d0; obr_=0d0
        END IF
        IF(nrr_meas>0)THEN
          ALLOCATE(obrr(nrr_meas),obrr_(nrr_meas))
          obrr=0d0; obrr_=0d0
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

                ! < c'(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c(i4,flv4,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obk=obk+g3(i2,i1,flv1)*g3old(i4,i3,flv3)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk=obk+gt3(i4,i1,flv1)*gt2(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
                END IF

                ! < c'(i3,flv3,t) c(i4,flv4,t) c'(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obk_=obk_+g3old(i2,i1,flv1)*g3(i4,i3,flv3)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk_=obk_+gt2(i4,i1,flv1)*gt3(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
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

                ! < c'(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c(i4,flv4,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obr(ir)=obr(ir)+g3(i2,i1,flv1)*g3old(i4,i3,flv3)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr(ir)=obr(ir)+gt3(i4,i1,flv1)*gt2(i2,i3,flv2)*factor
                END IF

                ! < c'(i3,flv3,t) c(i4,flv4,t) c'(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obr_(ir)=obr_(ir)+g3old(i2,i1,flv1)*g3(i4,i3,flv3)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr_(ir)=obr_(ir)+gt2(i4,i1,flv1)*gt3(i2,i3,flv2)*factor
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

              ! < c'(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c(i4,flv4,0) >
              IF(flv1==flv2.and.flv3==flv4)THEN
                obrr(ir)=obrr(ir)+g3(i2,i1,flv1)*g3old(i4,i3,flv3)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr(ir)=obrr(ir)+gt3(i4,i1,flv1)*gt2(i2,i3,flv2)*factor
              END IF

              ! < c'(i3,flv3,t) c(i4,flv4,t) c'(i1,flv1,0) c(i2,flv2,0) >
              IF(flv1==flv2.and.flv3==flv4)THEN
                obrr_(ir)=obrr_(ir)+g3old(i2,i1,flv1)*g3(i4,i3,flv3)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr_(ir)=obrr_(ir)+gt2(i4,i1,flv1)*gt3(i2,i3,flv2)*factor
              END IF

            END DO

          END IF

        END DO; END DO; END DO; END DO

        IF(nk_meas>0)THEN
          obk=obk*currentphase/(La*Lb*Lc)**2
          obk_=obk_*currentphase/(La*Lb*Lc)**2
          CALL zput_array(nk_meas,obk)
          CALL zput_array(nk_meas,obk_)
          DEALLOCATE(obk,obk_)
        END IF
        IF(nr_meas>0)THEN
          obr=obr*currentphase/(La*Lb*Lc)
          obr_=obr_*currentphase/(La*Lb*Lc)
          CALL zput_array(nr_meas,obr_)
          CALL zput_array(nr_meas,obr_)
          DEALLOCATE(obr,obr_)
        END IF
        IF(nrr_meas>0)THEN
          obrr=obrr*currentphase
          obrr_=obrr_*currentphase
          CALL zput_array(nrr_meas,obrr)
          CALL zput_array(nrr_meas,obrr_)
          DEALLOCATE(obrr,obrr_)
        END IF

      END DO

      !--------------------------------------------
      ! PP-channel two-particle Green's functions
      !--------------------------------------------

      DO k=1,n_pp_meas; d=ndim_pp_meas(k)

        IF(nk_meas>0)THEN
          ALLOCATE(obk(nk_meas),obk_(nk_meas))
          obk=0d0; obk_=0d0
        END IF
        IF(nr_meas>0)THEN
          ALLOCATE(obr(nr_meas),obr_(nr_meas))
          obr=0d0; obr_=0d0
        END IF
        IF(nrr_meas>0)THEN
          ALLOCATE(obrr(nrr_meas),obrr_(nr_meas))
          obrr=0d0; obrr_=0d0
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

                ! < c(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c'(i4,flv4,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obk=obk-gt2(i1,i3,flv1)*gt2(i2,i4,flv2)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk=obk+gt2(i1,i4,flv1)*gt2(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
                END IF

                ! < c'(i3,flv3,t) c'(i4,flv4,t) c(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obk_=obk_-gt3(i1,i3,flv1)*gt3(i2,i4,flv2)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk_=obk_+gt3(i1,i4,flv1)*gt3(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
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

                ! < c(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c'(i4,flv4,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obr(ir)=obr(ir)-gt2(i1,i3,flv1)*gt2(i2,i4,flv2)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr(ir)=obr(ir)+gt2(i1,i4,flv1)*gt2(i2,i3,flv2)*factor
                END IF

                ! < c'(i3,flv3,t) c'(i4,flv4,t) c(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obr_(ir)=obr_(ir)-gt3(i1,i3,flv1)*gt3(i2,i4,flv2)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr_(ir)=obr_(ir)+gt3(i1,i4,flv1)*gt3(i2,i3,flv2)*factor
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

              ! < c(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c'(i4,flv4,0) >
              IF(flv1==flv3.and.flv2==flv4)THEN
                obrr(ir)=obrr(ir)-gt2(i1,i3,flv1)*gt2(i2,i4,flv2)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr(ir)=obrr(ir)+gt2(i1,i4,flv1)*gt2(i2,i3,flv2)*factor
              END IF

              ! < c'(i3,flv3,t) c'(i4,flv4,t) c(i1,flv1,0) c(i2,flv2,0) >
              IF(flv1==flv3.and.flv2==flv4)THEN
                obrr_(ir)=obrr_(ir)-gt3(i1,i3,flv1)*gt3(i2,i4,flv2)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr_(ir)=obrr_(ir)+gt3(i1,i4,flv1)*gt3(i2,i3,flv2)*factor
              END IF

            END DO

          END IF

        END DO; END DO; END DO; END DO

        IF(nk_meas>0)THEN
          obk=obk*currentphase/(La*Lb*Lc)**2
          obk_=obk_*currentphase/(La*Lb*Lc)**2
          CALL zput_array(nk_meas,obk)
          CALL zput_array(nk_meas,obk_)
          DEALLOCATE(obk,obk_)
        END IF
        IF(nr_meas>0)THEN
          obr=obr*currentphase/(La*Lb*Lc)
          obr_=obr_*currentphase/(La*Lb*Lc)
          CALL zput_array(nr_meas,obr)
          CALL zput_array(nr_meas,obr_)
          DEALLOCATE(obr,obr_)
        END IF
        IF(nrr_meas>0)THEN
          obrr=obrr*currentphase
          obrr_=obrr*currentphase
          CALL zput_array(nrr_meas,obrr)
          CALL zput_array(nrr_meas,obrr_)
          DEALLOCATE(obrr,obrr_)
        END IF

      END DO
      
      !--------------------------------------------
      ! crossing-PH-channel two-particle Green's functions
      !--------------------------------------------

      DO ic=1,ncross_ph_meas
        
        k=cross_ph_meas(1,ic)
        kk=cross_ph_meas(2,ic)
        d=ndim_ph_meas(k)
        dd=ndim_ph_meas(kk)

        IF(nk_meas>0)THEN
          ALLOCATE(obk(nk_meas),obk_(nk_meas))
          obk=0d0; obk_=0d0
        END IF
        IF(nr_meas>0)THEN
          ALLOCATE(obr(nr_meas),obr_(nr_meas))
          obr=0d0; obr_=0d0
        END IF
        IF(nrr_meas>0)THEN
          ALLOCATE(obrr(nrr_meas),obrr_(nrr_meas))
          obrr=0d0; obrr_=0d0
        END IF

        DO k1=1,d; DO k2=1,d; DO k3=1,dd; DO k4=1,dd

          IF(abs(fmat_ph_meas(k1,k2,k))<1d-6)CYCLE
          IF(abs(fmat_ph_meas(k4,k3,kk))<1d-6)CYCLE
          factor=fmat_ph_meas(k1,k2,k)*conjg(fmat_ph_meas(k4,k3,kk))

          flv1=flv_ph_meas(k1,k)
          flv2=flv_ph_meas(k2,k)
          flv3=flv_ph_meas(k3,kk)
          flv4=flv_ph_meas(k4,kk)

          IF((.not.(flv1==flv2.and.flv3==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE

          IF(nk_meas>0)THEN

            DO a=1,La; DO b=1,Lb; DO c=1,Lc
              DO aa=1,La; DO bb=1,Lb; DO cc=1,Lc
                da=a-aa; db=b-bb; dc=c-cc

                i1=nb_ph_meas(a,b,c,k1,k)
                i2=nb_ph_meas(a,b,c,k2,k)
                i3=nb_ph_meas(aa,bb,cc,k3,kk)
                i4=nb_ph_meas(aa,bb,cc,k4,kk)

                ! < c'(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c(i4,flv4,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obk=obk+g3(i2,i1,flv1)*g3old(i4,i3,flv3)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk=obk+gt3(i4,i1,flv1)*gt2(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
                END IF

                ! < c'(i3,flv3,t) c(i4,flv4,t) c'(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obk_=obk_+g3old(i2,i1,flv1)*g3(i4,i3,flv3)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk_=obk_+gt2(i4,i1,flv1)*gt3(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
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
                i3=nb_ph_meas(aa,bb,cc,k3,kk)
                i4=nb_ph_meas(aa,bb,cc,k4,kk)

                ! < c'(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c(i4,flv4,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obr(ir)=obr(ir)+g3(i2,i1,flv1)*g3old(i4,i3,flv3)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr(ir)=obr(ir)+gt3(i4,i1,flv1)*gt2(i2,i3,flv2)*factor
                END IF

                ! < c'(i3,flv3,t) c(i4,flv4,t) c'(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv2.and.flv3==flv4)THEN
                  obr_(ir)=obr_(ir)+g3old(i2,i1,flv1)*g3(i4,i3,flv3)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr_(ir)=obr_(ir)+gt2(i4,i1,flv1)*gt3(i2,i3,flv2)*factor
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
              i3=nb_ph_meas(aa,bb,cc,k3,kk)
              i4=nb_ph_meas(aa,bb,cc,k4,kk)

              ! < c'(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c(i4,flv4,0) >
              IF(flv1==flv2.and.flv3==flv4)THEN
                obrr(ir)=obrr(ir)+g3(i2,i1,flv1)*g3old(i4,i3,flv3)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr(ir)=obrr(ir)+gt3(i4,i1,flv1)*gt2(i2,i3,flv2)*factor
              END IF

              ! < c'(i3,flv3,t) c(i4,flv4,t) c'(i1,flv1,0) c(i2,flv2,0) >
              IF(flv1==flv2.and.flv3==flv4)THEN
                obrr_(ir)=obrr_(ir)+g3old(i2,i1,flv1)*g3(i4,i3,flv3)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr_(ir)=obrr_(ir)+gt2(i4,i1,flv1)*gt3(i2,i3,flv2)*factor
              END IF

            END DO

          END IF

        END DO; END DO; END DO; END DO

        IF(nk_meas>0)THEN
          obk=obk*currentphase/(La*Lb*Lc)**2
          obk_=obk_*currentphase/(La*Lb*Lc)**2
          CALL zput_array(nk_meas,obk)
          CALL zput_array(nk_meas,obk_)
          DEALLOCATE(obk,obk_)
        END IF
        IF(nr_meas>0)THEN
          obr=obr*currentphase/(La*Lb*Lc)
          obr_=obr_*currentphase/(La*Lb*Lc)
          CALL zput_array(nr_meas,obr_)
          CALL zput_array(nr_meas,obr_)
          DEALLOCATE(obr,obr_)
        END IF
        IF(nrr_meas>0)THEN
          obrr=obrr*currentphase
          obrr_=obrr_*currentphase
          CALL zput_array(nrr_meas,obrr)
          CALL zput_array(nrr_meas,obrr_)
          DEALLOCATE(obrr,obrr_)
        END IF

      END DO

      !--------------------------------------------
      ! crossing-PP-channel two-particle Green's functions
      !--------------------------------------------

      DO ic=1,ncross_pp_meas
        
        k=cross_pp_meas(1,ic)
        kk=cross_pp_meas(2,ic)
        d=ndim_pp_meas(k)
        dd=ndim_pp_meas(kk)

        IF(nk_meas>0)THEN
          ALLOCATE(obk(nk_meas),obk_(nk_meas))
          obk=0d0; obk_=0d0
        END IF
        IF(nr_meas>0)THEN
          ALLOCATE(obr(nr_meas),obr_(nr_meas))
          obr=0d0; obr_=0d0
        END IF
        IF(nrr_meas>0)THEN
          ALLOCATE(obrr(nrr_meas),obrr_(nr_meas))
          obrr=0d0; obrr_=0d0
        END IF

        DO k1=1,d; DO k2=1,d; DO k3=1,dd; DO k4=1,dd

          IF(abs(fmat_pp_meas(k1,k2,k))<1d-6)CYCLE
          IF(abs(fmat_pp_meas(k4,k3,kk))<1d-6)CYCLE
          factor=fmat_pp_meas(k1,k2,k)*conjg(fmat_pp_meas(k4,k3,kk))

          flv1=flv_pp_meas(k1,k)
          flv2=flv_pp_meas(k2,k)
          flv3=flv_pp_meas(k3,kk)
          flv4=flv_pp_meas(k4,kk)

          IF((.not.(flv1==flv3.and.flv2==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE

          IF(nk_meas>0)THEN

            DO a=1,La; DO b=1,Lb; DO c=1,Lc
              DO aa=1,La; DO bb=1,Lb; DO cc=1,Lc
                da=a-aa; db=b-bb; dc=c-cc

                i1=nb_pp_meas(a,b,c,k1,k)
                i2=nb_pp_meas(a,b,c,k2,k)
                i3=nb_pp_meas(aa,bb,cc,k3,kk)
                i4=nb_pp_meas(aa,bb,cc,k4,kk)

                ! < c(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c'(i4,flv4,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obk=obk-gt2(i1,i3,flv1)*gt2(i2,i4,flv2)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk=obk+gt2(i1,i4,flv1)*gt2(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
                END IF

                ! < c'(i3,flv3,t) c'(i4,flv4,t) c(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obk_=obk_-gt3(i1,i3,flv1)*gt3(i2,i4,flv2)*factor*expikr_array(:,da,db,dc)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk_=obk_+gt3(i1,i4,flv1)*gt3(i2,i3,flv2)*factor*expikr_array(:,da,db,dc)
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
                i3=nb_pp_meas(aa,bb,cc,k3,kk)
                i4=nb_pp_meas(aa,bb,cc,k4,kk)

                ! < c(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c'(i4,flv4,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obr(ir)=obr(ir)-gt2(i1,i3,flv1)*gt2(i2,i4,flv2)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr(ir)=obr(ir)+gt2(i1,i4,flv1)*gt2(i2,i3,flv2)*factor
                END IF

                ! < c'(i3,flv3,t) c'(i4,flv4,t) c(i1,flv1,0) c(i2,flv2,0) >
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obr_(ir)=obr_(ir)-gt3(i1,i3,flv1)*gt3(i2,i4,flv2)*factor
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obr_(ir)=obr_(ir)+gt3(i1,i4,flv1)*gt3(i2,i3,flv2)*factor
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
              i3=nb_pp_meas(aa,bb,cc,k3,kk)
              i4=nb_pp_meas(aa,bb,cc,k4,kk)

              ! < c(i1,flv1,t) c(i2,flv2,t) c'(i3,flv3,0) c'(i4,flv4,0) >
              IF(flv1==flv3.and.flv2==flv4)THEN
                obrr(ir)=obrr(ir)-gt2(i1,i3,flv1)*gt2(i2,i4,flv2)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr(ir)=obrr(ir)+gt2(i1,i4,flv1)*gt2(i2,i3,flv2)*factor
              END IF

              ! < c'(i3,flv3,t) c'(i4,flv4,t) c(i1,flv1,0) c(i2,flv2,0) >
              IF(flv1==flv3.and.flv2==flv4)THEN
                obrr_(ir)=obrr_(ir)-gt3(i1,i3,flv1)*gt3(i2,i4,flv2)*factor
              END IF
              IF(flv1==flv4.and.flv2==flv3)THEN
                obrr_(ir)=obrr_(ir)+gt3(i1,i4,flv1)*gt3(i2,i3,flv2)*factor
              END IF

            END DO

          END IF

        END DO; END DO; END DO; END DO

        IF(nk_meas>0)THEN
          obk=obk*currentphase/(La*Lb*Lc)**2
          obk_=obk_*currentphase/(La*Lb*Lc)**2
          CALL zput_array(nk_meas,obk)
          CALL zput_array(nk_meas,obk_)
          DEALLOCATE(obk,obk_)
        END IF
        IF(nr_meas>0)THEN
          obr=obr*currentphase/(La*Lb*Lc)
          obr_=obr_*currentphase/(La*Lb*Lc)
          CALL zput_array(nr_meas,obr)
          CALL zput_array(nr_meas,obr_)
          DEALLOCATE(obr,obr_)
        END IF
        IF(nrr_meas>0)THEN
          obrr=obrr*currentphase
          obrr_=obrr*currentphase
          CALL zput_array(nrr_meas,obrr)
          CALL zput_array(nrr_meas,obrr_)
          DEALLOCATE(obrr,obrr_)
        END IF

      END DO


    END DO

    g=g_save
    DEALLOCATE(g_save,g3old,gt2,gt3)
    IF(.not.(proj.and.FAtech)) DEALLOCATE(qmat2,qmat3,dvec2,dvec3,rmat2,lmat3)

  END IF

  ! measurement externally
  IF(do_measure_external) CALL measurement_external(time)

END SUBROUTINE
