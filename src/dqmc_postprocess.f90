SUBROUTINE postprocess()
  USE dqmc
  IMPLICIT NONE
  IF(do_postprocess_external) CALL postprocess_external()
  PRINT*,'Job is done with',nd,' cores'
END SUBROUTINE
