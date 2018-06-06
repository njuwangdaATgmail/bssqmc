SUBROUTINE measurement(time)
  USE dqmc
  USE hubhol
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  COMPLEX(8) g3(nsite,nsite,nflv),kinetic,double,Sqz,Sqx,Sqc,deltaij,Pph
  INTEGER i,j,a,a2,b,b2,flv

  DO flv=1,nflv
    g3(:,:,flv)=-matmul(inv_expk_half(:,:,flv),matmul(g(:,:,flv),expk_half(:,:,flv)))
    DO i=1,nsite
      g3(i,i,flv)=g3(i,i,flv)+1d0
    END DO
  END DO

  kinetic=0d0
  double=0d0
  DO i=1,nsite
    DO j=1,nsite
      kinetic=kinetic+sum(kmat(i,j,:)*g3(j,i,:))*currentphase/nsite
    END DO
    double=double+g3(i,i,1)*g3(i,i,2)*currentphase/nsite
  END DO
  CALL put_pool(kinetic)
  CALL put_pool(double)

  Sqz=0d0
  Sqx=0d0
  Sqc=0d0
  DO a=1,La; DO b=1,Lb; i=label(a,b,1,1)
    DO a2=1,La; DO b2=1,Lb; j=label(a2,b2,1,1)
      deltaij=0d0
      IF(i==j)deltaij=1d0
      Sqx=Sqx+2*g3(i,j,1)*(deltaij-g3(j,i,2))*(-1)**(a-a2+b-b2)*currentphase/nsite
      Sqz=Sqz+((g3(i,i,1)-g3(i,i,2))*(g3(j,j,1)-g3(j,j,2))+g3(i,j,1)*(deltaij-g3(j,i,1))+ &
        & g3(i,j,2)*(deltaij-g3(j,i,2)))*(-1)**(a-a2+b-b2)*currentphase/nsite
      Sqc=Sqc+((g3(i,i,1)+g3(i,i,2))*(g3(j,j,1)+g3(j,j,2))+g3(i,j,1)*(deltaij-g3(j,i,1))+ &
        & g3(i,j,2)*(deltaij-g3(j,i,2)))*(-1)**(a-a2+b-b2)*currentphase/nsite
    END DO; END DO
  END DO; END DO
  CALL put_pool(Sqz)
  CALL put_pool(Sqx)
  CALL put_pool(Sqc)

  Pph=gph_x2(1)/beta*sum(field(:,:,2)**2)*currentphase/nsite
  CALL put_pool(Pph)

  g3=g3*currentphase
  CALL put_pool(2,g3(1:2,1,1))
  CALL put_pool(2,g3(1:2,1,2))

END SUBROUTINE
