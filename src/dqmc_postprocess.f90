SUBROUTINE postprocess()
  USE dqmc
  IMPLICIT NONE
  INTEGER k,i,tau,a,b,flv
  COMPLEX(8) x,dx
  COMPLEX(8), ALLOCATABLE :: obk4(:,:,:,:,:), obk4_err(:,:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: obr4(:,:,:,:,:), obr4_err(:,:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: obrr4(:,:,:,:,:), obrr4_err(:,:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: obk_ph(:,:), obk_ph_err(:,:)
  COMPLEX(8), ALLOCATABLE :: obr_ph(:,:), obr_ph_err(:,:)
  COMPLEX(8), ALLOCATABLE :: obrr_ph(:,:), obrr_ph_err(:,:)
  COMPLEX(8), ALLOCATABLE :: obk_pp(:,:), obk_pp_err(:,:)
  COMPLEX(8), ALLOCATABLE :: obr_pp(:,:), obr_pp_err(:,:)
  COMPLEX(8), ALLOCATABLE :: obrr_pp(:,:), obrr_pp_err(:,:)
  IF(id.ne.0)RETURN

  CALL get_pool(x,dx)
  PRINT*, 'sign average:', x, dx
  CALL get_pool(x,dx)
  PRINT*, 'kinetic energy:', x, dx

  !---------------------------------------------
  !  get averaged values from the measurement pool
  !---------------------------------------------

  IF(nk_meas>0)THEN
    ALLOCATE(obk4(nk_meas,norb,norb,nflv,2*ntau_meas+1))
    ALLOCATE(obk4_err(nk_meas,norb,norb,nflv,2*ntau_meas+1))
    CALL zget_array_err(nk_meas*norb*norb*nflv,obk4(:,:,:,:,0),obk4_err(:,:,:,:,0))
  END IF

  IF(nr_meas>0)THEN
    ALLOCATE(obr4(nr_meas,norb,norb,nflv,2*ntau_meas+1))
    ALLOCATE(obr4_err(nr_meas,norb,norb,nflv,2*ntau_meas+1))
    CALL zget_array_err(nr_meas*norb*norb*nflv,obr4(:,:,:,:,0),obr4_err(:,:,:,:,0))
  END IF

  IF(nrr_meas>0)THEN
    ALLOCATE(obrr4(nrr_meas,norb,norb,nflv,2*ntau_meas+1))
    ALLOCATE(obrr4_err(nrr_meas,norb,norb,nflv,2*ntau_meas+1))
    CALL zget_array_err(nrr_meas*norb*norb*nflv,obrr4(:,:,:,:,0),obrr4_err(:,:,:,:,0))
  END IF

  DO k=1,n_ph_meas

    IF(nk_meas>0)THEN
      ALLOCATE(obk_ph(nk_meas,2*ntau_meas+1),obk_ph_err(nk_meas,2*ntau_meas+1))
      CALL zget_array_err(nk_meas,obk_ph(:,0),obk_ph_err(:,0))
    END IF

    IF(nr_meas>0)THEN
      ALLOCATE(obr_ph(nr_meas,2*ntau_meas+1),obr_ph_err(nr_meas,2*ntau_meas+1))
      CALL zget_array_err(nr_meas,obr_ph(:,0),obr_ph_err(:,0))
    END IF

    IF(nrr_meas>0)THEN
      ALLOCATE(obrr_ph(nrr_meas,2*ntau_meas+1),obrr_ph_err(nrr_meas,2*ntau_meas+1))
      CALL zget_array_err(nrr_meas,obrr_ph(:,0),obrr_ph_err(:,0))
    END IF

  END DO

  DO k=1,n_pp_meas

    IF(nk_meas>0)THEN
      ALLOCATE(obk_pp(nk_meas,2*ntau_meas+1),obk_pp_err(nk_meas,2*ntau_meas+1))
      CALL zget_array_err(nk_meas,obk_pp(:,0),obk_pp_err(:,0))
    END IF

    IF(nr_meas>0)THEN
      ALLOCATE(obr_pp(nr_meas,2*ntau_meas+1),obr_pp_err(nr_meas,2*ntau_meas+1))
      CALL zget_array_err(nr_meas,obr_pp(:,0),obr_pp_err(:,0))
    END IF

    IF(nrr_meas>0)THEN
      ALLOCATE(obrr_pp(nrr_meas,2*ntau_meas+1),obrr_pp_err(nrr_meas,2*ntau_meas+1))
      CALL zget_array_err(nrr_meas,obrr_pp(:,0),obrr_pp_err(:,0))
    END IF

  END DO

  DO tau=1,ntau_meas

    IF(nk_meas>0)THEN
      CALL zget_array_err(nk_meas*norb*norb*nflv,obk4(:,:,:,:,tau),obk4_err(:,:,:,:,tau))
      CALL zget_array_err(nk_meas*norb*norb*nflv,obk4(:,:,:,:,-tau),obk4_err(:,:,:,:,-tau))
    END IF

    IF(nr_meas>0)THEN
      CALL zget_array_err(nr_meas*norb*norb*nflv,obr4(:,:,:,:,tau),obr4_err(:,:,:,:,tau))
      CALL zget_array_err(nr_meas*norb*norb*nflv,obr4(:,:,:,:,-tau),obr4_err(:,:,:,:,-tau))
    END IF

    IF(nrr_meas>0)THEN
      CALL zget_array_err(nrr_meas*norb*norb*nflv,obrr4(:,:,:,:,tau),obrr4_err(:,:,:,:,tau))
      CALL zget_array_err(nrr_meas*norb*norb*nflv,obrr4(:,:,:,:,-tau),obrr4_err(:,:,:,:,-tau))
    END IF

    DO k=1,n_ph_meas

      IF(nk_meas>0)THEN
        CALL zget_array_err(nk_meas,obk_ph(:,tau),obk_ph_err(:,tau))
        CALL zget_array_err(nk_meas,obk_ph(:,-tau),obk_ph_err(:,-tau))
      END IF

      IF(nr_meas>0)THEN
        CALL zget_array_err(nr_meas,obr_ph(:,tau),obr_ph_err(:,tau))
        CALL zget_array_err(nr_meas,obr_ph(:,-tau),obr_ph_err(:,-tau))
      END IF

      IF(nrr_meas>0)THEN
        CALL zget_array_err(nrr_meas,obrr_ph(:,tau),obrr_ph_err(:,tau))
        CALL zget_array_err(nrr_meas,obrr_ph(:,-tau),obrr_ph_err(:,-tau))
      END IF

    END DO

    DO k=1,n_pp_meas

      IF(nk_meas>0)THEN
        CALL zget_array_err(nk_meas,obk_pp(:,tau),obk_pp_err(:,tau))
        CALL zget_array_err(nk_meas,obk_pp(:,-tau),obk_pp_err(:,-tau))
      END IF

      IF(nr_meas>0)THEN
        CALL zget_array_err(nr_meas,obr_pp(:,tau),obr_pp_err(:,tau))
        CALL zget_array_err(nr_meas,obr_pp(:,-tau),obr_pp_err(:,-tau))
      END IF

      IF(nrr_meas>0)THEN
        CALL zget_array_err(nrr_meas,obrr_pp(:,tau),obrr_pp_err(:,tau))
        CALL zget_array_err(nrr_meas,obrr_pp(:,-tau),obrr_pp_err(:,-tau))
      END IF

    END DO

  END DO

  !---------------------
  ! output to files
  !---------------------

  IF(nk_meas>0)THEN
    OPEN(10,file='gk.dat')
    DO flv=1,nflv
      DO b=1,norb
        DO a=1,norb
          DO i=1,nk_meas
            DO tau=-ntau_meas,ntau_meas
              WRITE(10,'(4f8.3,3i6,4f12.6)') dtau*tau, k_array(1,i)*1d0/La, k_array(2,i)*1d0/Lb, k_array(3,i)*1d0/Lc, &
                & a, b, flv, obk4(i,a,b,flv,tau), obk4_err(i,a,b,flv,tau)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE(10)
  END IF

  IF(nr_meas>0)THEN
    OPEN(10,file='gr.dat')
    DO flv=1,nflv
      DO b=1,norb
        DO a=1,norb
          DO i=1,nr_meas
            DO tau=-ntau_meas,ntau_meas
              WRITE(10,'(1f8.3,6i6,4f12.6)') dtau*tau, r_array(:,i), a, b, flv, obr4(i,a,b,flv,tau), obr4_err(i,a,b,flv,tau)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE(10)
  END IF

  IF(nrr_meas>0)THEN
    OPEN(10,file='grr.dat')
    DO flv=1,nflv
      DO b=1,norb
        DO a=1,norb
          DO i=1,nrr_meas
            DO tau=-ntau_meas,ntau_meas
              WRITE(10,'(1f8.3,9i6,4f12.6)') dtau*tau, rr_array(:,1,i), rr_array(:,2,i), a, b, flv, &
                & obrr4(i,a,b,flv,tau), obrr4_err(i,a,b,flv,tau)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE(10)
  END IF

  DO k=1,n_ph_meas

    IF(nk_meas>0)THEN
      OPEN(10,file=trim(adjustl(name_ph_meas(k)))//'_k.dat')
      DO i=1,nk_meas
        DO tau=-ntau_meas,ntau_meas
          WRITE(10,'(4f8.3,4f12.6)') dtau*tau, k_array(1,i)*1d0/La, k_array(2,i)*1d0/Lb, k_array(3,i)*1d0/Lc, &
            & obk_ph(i,tau),obk_ph_err(i,tau)
        END DO
      END DO
      CLOSE(10)
    END IF

    IF(nr_meas>0)THEN
      OPEN(10,file=trim(adjustl(name_ph_meas(k)))//'_r.dat')
      DO i=1,nr_meas
        DO tau=-ntau_meas,ntau_meas
          WRITE(10,'(1f8.3,3i6,4f12.6)') dtau*tau, r_array(:,i), obr_ph(i,tau),obr_ph_err(i,tau)
        END DO
      END DO
      CLOSE(10)
    END IF

    IF(nrr_meas>0)THEN
      OPEN(10,file=trim(adjustl(name_ph_meas(k)))//'_rr.dat')
      DO i=1,nrr_meas
        DO tau=-ntau_meas,ntau_meas
          WRITE(10,'(1f8.3,6i6,4f12.6)') dtau*tau, rr_array(:,1,i), rr_array(:,2,i), obrr_ph(i,tau),obrr_ph_err(i,tau)
        END DO
      END DO
      CLOSE(10)
    END IF

  END DO

  DO k=1,n_pp_meas

    IF(nk_meas>0)THEN
      OPEN(10,file=trim(adjustl(name_pp_meas(k)))//'_k.dat')
      DO i=1,nk_meas
        DO tau=-ntau_meas,ntau_meas
          WRITE(10,'(4f8.3,4f12.6)') dtau*tau, k_array(1,i)*1d0/La, k_array(2,i)*1d0/Lb, k_array(3,i)*1d0/Lc, &
            & obk_pp(i,tau),obk_pp_err(i,tau)
        END DO
      END DO
      CLOSE(10)
    END IF

    IF(nr_meas>0)THEN
      OPEN(10,file=trim(adjustl(name_pp_meas(k)))//'_r.dat')
      DO i=1,nr_meas
        DO tau=-ntau_meas,ntau_meas
          WRITE(10,'(1f8.3,3i6,4f12.6)') dtau*tau, r_array(:,i), obr_pp(i,tau),obr_pp_err(i,tau)
        END DO
      END DO
      CLOSE(10)
    END IF

    IF(nrr_meas>0)THEN
      OPEN(10,file=trim(adjustl(name_pp_meas(k)))//'_rr.dat')
      DO i=1,nrr_meas
        DO tau=-ntau_meas,ntau_meas
          WRITE(10,'(1f8.3,6i6,4f12.6)') dtau*tau, rr_array(:,1,i), rr_array(:,2,i), obrr_pp(i,tau),obrr_pp_err(i,tau)
        END DO
      END DO
      CLOSE(10)
    END IF

  END DO

  IF(do_postprocess_external) CALL postprocess_external()
  PRINT*,'Job is done with',nd,' cores'

END SUBROUTINE
