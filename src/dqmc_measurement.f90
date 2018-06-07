SUBROUTINE measurement(time)
  USE dqmc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER i,j,k,orb,flv,a,b,c,aa,bb,cc,a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,da,db,dc
  INTEGER i1,i2,i3,i4,k1,k2,k3,k4,flv1,flv2,flv3,flv4,orb1,orb2,orb3,orb4,d,ir
  COMPLEX(8) kinetic,factor,g2(nsite,nsite,nflv), g3(nsite,nsite,nflv)
  COMPLEX(8), ALLOCATABLE :: mat4(:,:,:,:),obk(:),obr(:),obrr(:)
  
  ! get g2(i,j)=<c(i)c'(j)> and g3(i,j)=<c'(j)c(i)> in 2nd-order Trotter approximation
  ! NOTICE the positions of i and j. 
  ! Such a definition will be easier to calculate time evolutions.
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

  ! single-particle Green's functions
  
  IF(nk_meas>0)THEN  ! PBC is assumed

    ALLOCATE(mat4(nk_meas,norb,norb,nflv))
    mat4=0d0
    DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
      DO da=0,La-1; a2=a+da; IF(a2>La) a2=a2-La
        DO db=0,Lb-1; b2=b+db; IF(b2>Lb) b2=b2-Lb
          DO dc=0,Lc-1; c2=c+dc; IF(c2>Lc) c2=c2-Lc
            DO orb2=1,norb; i2=label(a2,b2,c2,orb2)
 
              DO k=1,nk_meas
                mat4(k,orb,orb2,:)=mat4(k,orb,orb2,:)+g2(i,i2,:)*expikr_array(da,db,dc,k)
              END DO
 
            END DO
          END DO
        END DO
      END DO
    END DO; END DO; END DO; END DO
    mat4=mat4*currentphase/(La*Lb*Lc)**2
    CALL zput_array(nk_meas*norb*norb*nflv,mat4)
    DEALLOCATE(mat4)

  END IF

  IF(nr_meas>0)THEN

    ALLOCATE(mat4(nr_meas,norb,norb,nflv))
    mat4=0d0
    DO a=1,La; DO b=1,Lb; DO c=1,Lc; DO orb=1,norb; i=label(a,b,c,orb)
      DO k=1,nr_meas
        a2=a+r_array(1,k); IF(a2>La) a2=a2-La
        b2=b+r_array(2,k); IF(b2>Lb) b2=b2-Lb
        c2=c+r_array(3,k); IF(c2>Lc) c2=c2-Lc
        DO orb2=1,norb; i2=label(a2,b2,c2,orb2)

          mat4(k,orb,orb2,:)=mat4(k,orb,orb2,:)+g2(i,i2,:)

        END DO
      END DO
    END DO; END DO; END DO; END DO
    mat4=mat4*currentphase/(La*Lb*Lc)
    CALL zput_array(nr_meas*norb*norb*nflv,mat4)
    DEALLOCATE(mat4)
    
  END IF

  IF(nrr_meas>0)THEN
    
    ALLOCATE(mat4(nrr_meas,norb,norb,nflv))
    DO orb=1,norb
      DO orb2=1,norb
        DO k=1,nrr_meas
          i=label(rr_array(1,1,k),rr_array(2,1,k),rr_array(3,1,k),orb)
          i2=label(rr_array(1,2,k),rr_array(2,2,k),rr_array(3,2,k),orb2)
          
          mat4(k,orb,orb2,:)=mat4(k,orb,orb2,:)+g2(i,i2,:)

        END DO
      END DO
    END DO
    mat4=mat4*currentphase
    CALL zput_array(nrr_meas*norb*norb*nflv,mat4)
    DEALLOCATE(mat4)

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
      
      IF(abs(fmat_ph_meas(k2,k1,k))<1d-6)CYCLE
      IF(abs(fmat_ph_meas(k3,k4,k))<1d-6)CYCLE
      factor=conjg(fmat_ph_meas(k2,k1,k))*fmat_ph_meas(k3,k4,k)
      
      flv1=flv_ph_meas(k1,k)
      flv2=flv_ph_meas(k2,k)
      flv3=flv_ph_meas(k3,k)
      flv4=flv_ph_meas(k4,k)

      IF((.not.(flv1==flv2.and.flv3==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE
      
      orb1=orb_ph_meas(k1,k)
      orb2=orb_ph_meas(k2,k)
      orb3=orb_ph_meas(k3,k)
      orb4=orb_ph_meas(k4,k)

      IF(nk_meas>0)THEN

        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO da=0,La-1; aa=a+da; IF(aa>La) aa=aa-La
            DO db=0,Lb-1; bb=b+db; IF(bb>Lb) bb=bb-Lb
              DO dc=0,Lc-1; cc=c+dc; IF(cc>Lc) cc=cc-Lc


                a1=a+da_ph_meas(k1,k); IF(a1>La) a1=a1-La
                b1=b+db_ph_meas(k1,k); IF(b1>Lb) b1=b1-Lb
                c1=c+dc_ph_meas(k1,k); IF(c1>Lc) c1=c1-Lc
                
                a2=a+da_ph_meas(k2,k); IF(a2>La) a2=a2-La
                b2=b+db_ph_meas(k2,k); IF(b2>Lb) b2=b2-Lb
                c2=c+dc_ph_meas(k2,k); IF(c2>Lc) c2=c2-Lc
                
                a3=aa+da_ph_meas(k3,k); IF(a3>La) a3=a3-La
                b3=bb+db_ph_meas(k3,k); IF(b3>Lb) b3=b3-Lb
                c3=cc+dc_ph_meas(k3,k); IF(c3>Lc) c3=c3-Lc
                
                a4=aa+da_ph_meas(k4,k); IF(a4>La) a4=a4-La
                b4=bb+db_ph_meas(k4,k); IF(b4>Lb) b4=b4-Lb
                c4=cc+dc_ph_meas(k4,k); IF(c4>Lc) c4=c4-Lc

                i1=label(a1,b1,c1,orb1)
                i2=label(a2,b2,c2,orb2)
                i3=label(a3,b3,c3,orb3)
                i4=label(a4,b4,c4,orb4)

                !   c'(r1,orb1,flv1) * conjg(fmat(k2,k1)) * c(r2,orb2,flv2) &
                ! * c'(r3,orb3,flv3) * fmat(k3,k4)        * c(r4,orb4,flv4)
                ! = c'(i1,flv1) * conjg(fmat(k2,k1)) * c(i2,flv2) &
                ! * c'(i3,flv3) * fmat(k3,k4)        * c(i4,flv4)

                IF(flv1==flv2.and.flv3==flv4)THEN
                  obk=obk+g3(i2,i1,flv1)*g3(i4,i3,flv3)*factor*expikr_array(da,db,dc,:)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk=obk+g3(i4,i1,flv1)*g2(i2,i3,flv2)*factor*expikr_array(da,db,dc,:)
                END IF

              END DO
            END DO
          END DO
        END DO; END DO; END DO

      END IF

      IF(nr_meas>0)THEN
 
        DO a=1,La; DO b=1,Lb; DO c=1,Lc
 
          DO ir=1,nr_meas
            da=r_array(1,ir); aa=a+da; IF(aa>La) aa=aa-La
            db=r_array(2,ir); bb=b+db; IF(bb>Lb) bb=bb-Lb
            dc=r_array(3,ir); cc=c+dc; IF(cc>Lc) cc=cc-Lc
 
            a1=a+da_ph_meas(k1,k); IF(a1>La) a1=a1-La
            b1=b+db_ph_meas(k1,k); IF(b1>Lb) b1=b1-Lb
            c1=c+dc_ph_meas(k1,k); IF(c1>Lc) c1=c1-Lc
            
            a2=a+da_ph_meas(k2,k); IF(a2>La) a2=a2-La
            b2=b+db_ph_meas(k2,k); IF(b2>Lb) b2=b2-Lb
            c2=c+dc_ph_meas(k2,k); IF(c2>Lc) c2=c2-Lc
            
            a3=aa+da_ph_meas(k3,k); IF(a3>La) a3=a3-La
            b3=bb+db_ph_meas(k3,k); IF(b3>Lb) b3=b3-Lb
            c3=cc+dc_ph_meas(k3,k); IF(c3>Lc) c3=c3-Lc
            
            a4=aa+da_ph_meas(k4,k); IF(a4>La) a4=a4-La
            b4=bb+db_ph_meas(k4,k); IF(b4>Lb) b4=b4-Lb
            c4=cc+dc_ph_meas(k4,k); IF(c4>Lc) c4=c4-Lc
 
            i1=label(a1,b1,c1,orb1)
            i2=label(a2,b2,c2,orb2)
            i3=label(a3,b3,c3,orb3)
            i4=label(a4,b4,c4,orb4)
 
            !   c'(r1,orb1,flv1) * conjg(fmat(k2,k1)) * c(r2,orb2,flv2) &
            ! * c'(r3,orb3,flv3) * fmat(k3,k4)        * c(r4,orb4,flv4)
            ! = c'(i1,flv1) * conjg(fmat(k2,k1)) * c(i2,flv2) &
            ! * c'(i3,flv3) * fmat(k3,k4)        * c(i4,flv4)
 
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
 
          a1=a+da_ph_meas(k1,k); IF(a1>La) a1=a1-La
          b1=b+db_ph_meas(k1,k); IF(b1>Lb) b1=b1-Lb
          c1=c+dc_ph_meas(k1,k); IF(c1>Lc) c1=c1-Lc
          
          a2=a+da_ph_meas(k2,k); IF(a2>La) a2=a2-La
          b2=b+db_ph_meas(k2,k); IF(b2>Lb) b2=b2-Lb
          c2=c+dc_ph_meas(k2,k); IF(c2>Lc) c2=c2-Lc
          
          a3=aa+da_ph_meas(k3,k); IF(a3>La) a3=a3-La
          b3=bb+db_ph_meas(k3,k); IF(b3>Lb) b3=b3-Lb
          c3=cc+dc_ph_meas(k3,k); IF(c3>Lc) c3=c3-Lc
          
          a4=aa+da_ph_meas(k4,k); IF(a4>La) a4=a4-La
          b4=bb+db_ph_meas(k4,k); IF(b4>Lb) b4=b4-Lb
          c4=cc+dc_ph_meas(k4,k); IF(c4>Lc) c4=c4-Lc
 
          i1=label(a1,b1,c1,orb1)
          i2=label(a2,b2,c2,orb2)
          i3=label(a3,b3,c3,orb3)
          i4=label(a4,b4,c4,orb4)
 
          !   c'(r1,orb1,flv1) * conjg(fmat(k2,k1)) * c(r2,orb2,flv2) &
          ! * c'(r3,orb3,flv3) * fmat(k3,k4)        * c(r4,orb4,flv4)
          ! = c'(i1,flv1) * conjg(fmat(k2,k1)) * c(i2,flv2) &
          ! * c'(i3,flv3) * fmat(k3,k4)        * c(i4,flv4)
 
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
      
      IF(abs(fmat_pp_meas(k2,k1,k))<1d-6)CYCLE
      IF(abs(fmat_pp_meas(k3,k4,k))<1d-6)CYCLE
      factor=conjg(fmat_pp_meas(k2,k1,k))*fmat_pp_meas(k3,k4,k)
      
      flv1=flv_pp_meas(k1,k)
      flv2=flv_pp_meas(k2,k)
      flv3=flv_pp_meas(k3,k)
      flv4=flv_pp_meas(k4,k)

      IF((.not.(flv1==flv2.and.flv3==flv4)).and.(.not.(flv1==flv4.and.flv2==flv3)))CYCLE
      
      orb1=orb_pp_meas(k1,k)
      orb2=orb_pp_meas(k2,k)
      orb3=orb_pp_meas(k3,k)
      orb4=orb_pp_meas(k4,k)

      IF(nk_meas>0)THEN

        DO a=1,La; DO b=1,Lb; DO c=1,Lc
          DO da=0,La-1; aa=a+da; IF(aa>La) aa=aa-La
            DO db=0,Lb-1; bb=b+db; IF(bb>Lb) bb=bb-Lb
              DO dc=0,Lc-1; cc=c+dc; IF(cc>Lc) cc=cc-Lc


                a1=a+da_pp_meas(k1,k); IF(a1>La) a1=a1-La
                b1=b+db_pp_meas(k1,k); IF(b1>Lb) b1=b1-Lb
                c1=c+dc_pp_meas(k1,k); IF(c1>Lc) c1=c1-Lc
                
                a2=a+da_pp_meas(k2,k); IF(a2>La) a2=a2-La
                b2=b+db_pp_meas(k2,k); IF(b2>Lb) b2=b2-Lb
                c2=c+dc_pp_meas(k2,k); IF(c2>Lc) c2=c2-Lc
                
                a3=aa+da_pp_meas(k3,k); IF(a3>La) a3=a3-La
                b3=bb+db_pp_meas(k3,k); IF(b3>Lb) b3=b3-Lb
                c3=cc+dc_pp_meas(k3,k); IF(c3>Lc) c3=c3-Lc
                
                a4=aa+da_pp_meas(k4,k); IF(a4>La) a4=a4-La
                b4=bb+db_pp_meas(k4,k); IF(b4>Lb) b4=b4-Lb
                c4=cc+dc_pp_meas(k4,k); IF(c4>Lc) c4=c4-Lc

                i1=label(a1,b1,c1,orb1)
                i2=label(a2,b2,c2,orb2)
                i3=label(a3,b3,c3,orb3)
                i4=label(a4,b4,c4,orb4)

                !   c'(r1,orb1,flv1) * conjg(fmat(k2,k1)) * c'(r2,orb2,flv2) &
                ! * c(r3,orb3,flv3)  * fmat(k3,k4)        * c(r4,orb4,flv4)
                ! = c'(i1,flv1) * conjg(fmat(k2,k1)) * c'(i2,flv2) &
                ! * c(i3,flv3)  * fmat(k3,k4)        * c(i4,flv4)
               
                IF(flv1==flv3.and.flv2==flv4)THEN
                  obk=obk-g3(i3,i1,flv1)*g3(i4,i2,flv2)*factor*expikr_array(da,db,dc,:)
                END IF
                IF(flv1==flv4.and.flv2==flv3)THEN
                  obk=obk+g3(i4,i1,flv1)*g3(i3,i2,flv2)*factor*expikr_array(da,db,dc,:)
                END IF

              END DO
            END DO
          END DO
        END DO; END DO; END DO

      END IF

      IF(nr_meas>0)THEN
 
        DO a=1,La; DO b=1,Lb; DO c=1,Lc
 
          DO ir=1,nr_meas
            da=r_array(1,ir); aa=a+da; IF(aa>La) aa=aa-La
            db=r_array(2,ir); bb=b+db; IF(bb>Lb) bb=bb-Lb
            dc=r_array(3,ir); cc=c+dc; IF(cc>Lc) cc=cc-Lc
 
            a1=a+da_pp_meas(k1,k); IF(a1>La) a1=a1-La
            b1=b+db_pp_meas(k1,k); IF(b1>Lb) b1=b1-Lb
            c1=c+dc_pp_meas(k1,k); IF(c1>Lc) c1=c1-Lc
            
            a2=a+da_pp_meas(k2,k); IF(a2>La) a2=a2-La
            b2=b+db_pp_meas(k2,k); IF(b2>Lb) b2=b2-Lb
            c2=c+dc_pp_meas(k2,k); IF(c2>Lc) c2=c2-Lc
            
            a3=aa+da_pp_meas(k3,k); IF(a3>La) a3=a3-La
            b3=bb+db_pp_meas(k3,k); IF(b3>Lb) b3=b3-Lb
            c3=cc+dc_pp_meas(k3,k); IF(c3>Lc) c3=c3-Lc
            
            a4=aa+da_pp_meas(k4,k); IF(a4>La) a4=a4-La
            b4=bb+db_pp_meas(k4,k); IF(b4>Lb) b4=b4-Lb
            c4=cc+dc_pp_meas(k4,k); IF(c4>Lc) c4=c4-Lc
 
            i1=label(a1,b1,c1,orb1)
            i2=label(a2,b2,c2,orb2)
            i3=label(a3,b3,c3,orb3)
            i4=label(a4,b4,c4,orb4)
 
            !   c'(r1,orb1,flv1) * conjg(fmat(k2,k1)) * c'(r2,orb2,flv2) &
            ! * c(r3,orb3,flv3)  * fmat(k3,k4)        * c(r4,orb4,flv4)
            ! = c'(i1,flv1) * conjg(fmat(k2,k1)) * c'(i2,flv2) &
            ! * c(i3,flv3)  * fmat(k3,k4)        * c(i4,flv4)
           
            IF(flv1==flv3.and.flv2==flv4)THEN
              obr(ir)=obr(ir)-g3(i3,i1,flv1)*g3(i4,i2,flv2)*factor
            END IF
            IF(flv1==flv4.and.flv2==flv3)THEN
              obr(ir)=obr(ir)+g3(i4,i1,flv1)*g3(i3,i2,flv2)*factor
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
 
          a1=a+da_pp_meas(k1,k); IF(a1>La) a1=a1-La
          b1=b+db_pp_meas(k1,k); IF(b1>Lb) b1=b1-Lb
          c1=c+dc_pp_meas(k1,k); IF(c1>Lc) c1=c1-Lc
          
          a2=a+da_pp_meas(k2,k); IF(a2>La) a2=a2-La
          b2=b+db_pp_meas(k2,k); IF(b2>Lb) b2=b2-Lb
          c2=c+dc_pp_meas(k2,k); IF(c2>Lc) c2=c2-Lc
          
          a3=aa+da_pp_meas(k3,k); IF(a3>La) a3=a3-La
          b3=bb+db_pp_meas(k3,k); IF(b3>Lb) b3=b3-Lb
          c3=cc+dc_pp_meas(k3,k); IF(c3>Lc) c3=c3-Lc
          
          a4=aa+da_pp_meas(k4,k); IF(a4>La) a4=a4-La
          b4=bb+db_pp_meas(k4,k); IF(b4>Lb) b4=b4-Lb
          c4=cc+dc_pp_meas(k4,k); IF(c4>Lc) c4=c4-Lc
 
          i1=label(a1,b1,c1,orb1)
          i2=label(a2,b2,c2,orb2)
          i3=label(a3,b3,c3,orb3)
          i4=label(a4,b4,c4,orb4)
 
          !   c'(r1,orb1,flv1) * conjg(fmat(k2,k1)) * c'(r2,orb2,flv2) &
          ! * c(r3,orb3,flv3)  * fmat(k3,k4)        * c(r4,orb4,flv4)
          ! = c'(i1,flv1) * conjg(fmat(k2,k1)) * c'(i2,flv2) &
          ! * c(i3,flv3)  * fmat(k3,k4)        * c(i4,flv4)
 
          IF(flv1==flv3.and.flv2==flv4)THEN
            obrr(ir)=obrr(ir)-g3(i3,i1,flv1)*g3(i4,i2,flv2)*factor
          END IF
          IF(flv1==flv4.and.flv2==flv3)THEN
            obrr(ir)=obrr(ir)+g3(i4,i1,flv1)*g3(i3,i2,flv2)*factor
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


  ! unequal-time Green's functions
  ! to be added

  ! measurement externally
  IF(do_measure_external) CALL measurement_external(time)

END SUBROUTINE
