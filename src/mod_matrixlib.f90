MODULE matrixlib

  INTERFACE trace
    MODULE PROCEDURE dtrace,ztrace
  END INTERFACE

  INTERFACE qdr
    MODULE PROCEDURE dqdr,zqdr
  END INTERFACE

  INTERFACE ldq
    MODULE PROCEDURE dldq,zldq
  END INTERFACE

  INTERFACE eigen
    MODULE PROCEDURE deigen,zeigen
  END INTERFACE

  INTERFACE inverse
    MODULE PROCEDURE dinverse,zinverse
  END INTERFACE

  INTERFACE qr
    MODULE PROCEDURE dqr,zqr
  END INTERFACE

  INTERFACE lq
    MODULE PROCEDURE dlq,zlq
  END INTERFACE

  INTERFACE det
    MODULE PROCEDURE ddet,zdet
  END INTERFACE

CONTAINS
!--------------------------------------------------
!   general functins
!--------------------------------------------------
 FUNCTION dtrace(n,a)
    IMPLICIT NONE
    INTEGER n,i
    REAL(8) a(n,n),dtrace
    dtrace=0d0
    DO i=1,n
      dtrace=dtrace+a(i,i)
    END DO
  END FUNCTION

  FUNCTION ztrace(n,a)
    IMPLICIT NONE
    INTEGER n,i
    COMPLEX(8) a(n,n),ztrace
    ztrace=(0d0,0d0)
    DO i=1,n
      ztrace=ztrace+a(i,i)
    END DO
  END FUNCTION

  SUBROUTINE dqdr(m,n,a,r,d)
    IMPLICIT NONE
    INTEGER m,n,i
    REAL(8) a(m,n),r(n,n),d(n)

    CALL dqr(m,n,a,r)
    DO i=1,n
      d(i)=r(i,i)
      r(i,i:n)=r(i,i:n)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE zqdr(m,n,a,r,d)
    IMPLICIT NONE
    INTEGER m,n,i
    COMPLEX(8) a(m,n),r(n,n),d(n)

    CALL zqr(m,n,a,r)
    DO i=1,n
      d(i)=r(i,i)
      r(i,i:n)=r(i,i:n)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE dldq(m,n,a,l,d)
    IMPLICIT NONE
    INTEGER m,n,i
    REAL(8) a(m,n),l(m,m),d(m)

    CALL dlq(m,n,a,l)
    DO i=1,m
      d(i)=l(i,i)
      l(i:m,i)=l(i:m,i)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE zldq(m,n,a,l,d)
    IMPLICIT NONE
    INTEGER m,n,i
    COMPLEX(8) a(m,n),l(m,m),d(m)

    CALL zlq(m,n,a,l)
    DO i=1,m
      d(i)=l(i,i)
      l(i:m,i)=l(i:m,i)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE deigen(n,a,v)
    IMPLICIT NONE
    INTEGER n,info
    REAL(8) a(n,n),v(n),work(64*n)
    IF(n<=0)STOP 'N<=0 is invalid, inside DEIGEN'
    IF(n==1)THEN
      v(1)=a(1,1);a(1,1)=1d0;RETURN
    END IF
    CALL dsyev('V','U',n,a,n,v,work,64*n,info)
    IF(info/=0)STOP 'ERROR @ DSYEV, inside DEIGEN'
  END SUBROUTINE

  SUBROUTINE zeigen(n,a,v)
    IMPLICIT NONE
    INTEGER n,info
    REAL(8) rwork(3*n-2),v(n)
    COMPLEX(8) work(64*n),a(n,n)
    IF(n<=0)STOP 'N<=0 is invalid, inside DEIGEN'
    IF(n==1)THEN
      v(1)=real(a(1,1));a(1,1)=1d0;RETURN
    END IF
    CALL zheev('v','u',n,a,n,v,work,64*n,rwork,info)
    IF(info/=0)STOP 'ERROR @ ZHEEV, inside ZEIGEN'
  END SUBROUTINE

  SUBROUTINE dinverse(n,a)
    IMPLICIT NONE
    INTEGER n,info,ipiv(n)
    REAL(8) a(n,n),work(64*n)
    IF(n<=0)STOP 'N<=0 is invalid, inside DINVERSE'
    IF(n==1)THEN
      a(1,1)=1d0/a(1,1);RETURN
    END IF
    CALL dgetrf(n,n,a,n,ipiv,info)
    IF(info/=0)STOP 'ERROR @ DGETRF, inside DINVERSE'
    CALL dgetri(n,a,n,ipiv,work,64*n,info)
    IF(info/=0)STOP 'ERROR @ DGETRI, inside DINVERSE'
  END SUBROUTINE

  SUBROUTINE zinverse(n,a)
    IMPLICIT NONE
    INTEGER n,info,ipiv(n)
    COMPLEX(8) a(n,n),work(64*n)
    IF(n<=0)STOP 'N<=0 is invalid, inside ZINVERSE'
    IF(n==1)THEN
      a(1,1)=1d0/a(1,1);RETURN
    END IF
    CALL zgetrf(n,n,a,n,ipiv,info)
    IF(info/=0)STOP 'ERROR @ ZGETRF, inside ZINVERSE'
    CALL zgetri(n,a,n,ipiv,work,64*n,info)
    IF(info/=0)STOP 'ERROR @ ZGETRI, inside ZINVERSE'
  END SUBROUTINE

  SUBROUTINE dqr(m,n,a,r)
    IMPLICIT NONE
    INTEGER m,n,info,j,lwork
    REAL(8) a(m,n),r(n,n),tau(n),work(n*64)

    IF(m<n) STOP 'm<n is invalid, inside DQR'
    IF(n<=0) STOP 'n<=0 is invalid, inside DQR'
    lwork=n*64
    CALL dgeqrf(m,n,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZGEQRF, inside DQR'
    r=0d0
    DO j=1,n
      r(1:j,j)=a(1:j,j)
    END DO
    CALL dorgqr(m,n,n,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZUNGQR, inside DQR'

  END SUBROUTINE

  SUBROUTINE zqr(m,n,a,r)
    IMPLICIT NONE
    INTEGER m,n,info,j,lwork
    COMPLEX(8) a(m,n),r(n,n),tau(n),work(n*64)
    IF(m<n) STOP 'm<n is invalid, inside ZQR'
    IF(n<=0) STOP 'n<=0 is invalid, inside ZQR'
    lwork=n*64
    CALL zgeqrf(m,n,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZGEQRF, inside ZQR'
    r=0d0
    DO j=1,n
      r(1:j,j)=a(1:j,j)
    END DO
    CALL zungqr(m,n,n,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZUNGQR, inside ZQR'
  END SUBROUTINE

  SUBROUTINE dlq(m,n,a,l)
    IMPLICIT NONE
    INTEGER m,n,info,j,lwork
    REAL(8) a(m,n),l(m,m),tau(m),work(m*64)

    IF(m>n) STOP 'm>n is invalid, inside DLQ'
    IF(m<=0) STOP 'n<=0 is invalid, inside DLQ'
    lwork=m*64
    CALL dgelqf(m,n,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZGELQF, inside DLQ'
    l=0d0
    DO j=1,m
      l(j:m,j)=a(j:m,j)
    END DO
    CALL dorglq(m,n,m,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZUNGLQ, inside DLQ'

  END SUBROUTINE

  SUBROUTINE zlq(m,n,a,l)
    IMPLICIT NONE
    INTEGER m,n,info,j,lwork
    COMPLEX(8) a(m,n),l(m,m),tau(m),work(m*64)
    IF(m>n) STOP 'm>n is invalid, inside ZLQ'
    IF(m<=0) STOP 'm<=0 is invalid, inside ZLQ'
    lwork=m*64
    CALL zgelqf(m,n,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZGELQF, inside ZLQ'
    l=0d0
    DO j=1,m
      l(j:m,j)=a(j:m,j)
    END DO
    CALL zunglq(m,n,m,a,m,tau,work,lwork,info)
    IF(info/=0)STOP 'ERROR @ ZUNGLQ, inside ZLQ'
  END SUBROUTINE

!========================================================================================
! IMPORTANT: there may be unstability problems in these subroutines solving determinant.
! In fact, calculating determinant is not recommendated in real calculations.
!========================================================================================

  FUNCTION ddet(n,a)
    IMPLICIT NONE
    INTEGER n,i,info,ipvt(n)
    REAL(8) a(n,n),b(n,n),ddet
    IF(n<=0)STOP 'N<0 is invalid, inside DDET'
    IF(n==1)THEN
      ddet=a(1,1);RETURN
    END IF
    b=a
    CALL dgetrf(n,n,b,n,ipvt,info)
    IF(info/=0)THEN
      ddet=0d0;RETURN
    END IF
    info=1
    ddet=0d0
    DO i=1,n
      IF(ipvt(i)/=i)info=-info
      IF(b(i,i)<0d0)THEN
        info=-info
        b(i,i)=-b(i,i)
      END IF
      ddet=ddet+log(b(i,i))
    END DO
    ddet=exp(ddet)
    IF(info<0)ddet=-ddet
  END FUNCTION

  FUNCTION zdet(n,a)
    IMPLICIT NONE
    INTEGER i,n,info,ipvt(n)
    COMPLEX(8) a(n,n),b(n,n),zdet
    IF(n<=0)STOP 'N<0 is invalid, inside ZDET'
    IF(n==1)THEN
      zdet=a(1,1);RETURN
    END IF
    b=a
    CALL zgetrf(n,n,b,n,ipvt,info)
    IF(info/=0)THEN
      zdet=0d0;RETURN
    END IF
    info=1
    zdet=1d0
    !zdet=0d0
    DO i=1,n
      IF(ipvt(i)/=i)info=-info
      zdet=zdet*b(i,i)
      !zdet=zdet+log(b(i,i))
    END DO
    !zdet=exp(zdet)
    IF(info<0)zdet=-zdet
  END FUNCTION

END MODULE

