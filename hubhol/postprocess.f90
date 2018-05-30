SUBROUTINE postprocess()
  USE dqmc
  IMPLICIT NONE
  COMPLEX(8) x,dx,y(2),dy(2)
  IF(id.ne.0)RETURN
  CALL get_pool(x,dx)
  PRINT*,'sign average:',x,dx
  CALL get_pool(x,dx)
  PRINT*,'kinetic energy:',x,dx
  CALL get_pool(x,dx)
  PRINT*,'double occupation:',x,dx
  CALL get_pool(x,dx)
  PRINT*,'AF structure factor (zz):',x,dx
  CALL get_pool(x,dx)
  PRINT*,'AF structure factor (xx):',x,dx
  CALL get_pool(x,dx)
  PRINT*,'CDW structure factor:',x,dx
  CALL get_pool(x,dx)
  PRINT*,'phonon potential energy:',x,dx
  CALL get_pool(2,y,dy)
  PRINT*,'g3(1,1,up):',y(1),dy(1)
  PRINT*,'g3(2,1,up):',y(2),dy(2)
  CALL get_pool(2,y,dy)
  PRINT*,'g3(1,1,dn):',y(1),dy(1)
  PRINT*,'g3(2,1,dn):',y(2),dy(2)
  PRINT*,'Job is done with',nd,' cores'
END SUBROUTINE
