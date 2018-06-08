! In this subroutine, we can output some quantities (e.g. boson field configurations
! or correlation functions) to examine the MC process or obtain the histograms.
! We can also perform calculations of entanglement entropy using swap technique,
! see Grover, PRL 111, 130402 (2013) and later works.
! Of course, we can also leave these work to users via tmpout_external().
SUBROUTINE tmpout()
  USE dqmc
  IMPLICIT NONE
  IF(do_tmpout_external) CALL tmpout_external()
END SUBROUTINE
