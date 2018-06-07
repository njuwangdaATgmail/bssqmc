SUBROUTINE tmpout()
  USE dqmc
  IMPLICIT NONE
  IF(do_tmpout_external)CALL tmpout_external()
END SUBROUTINE
