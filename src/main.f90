! MAIN FILE
program main
  use global_variables
  implicit none

  call init_parallel

  call parameters
  call operators

  call initial_state
!  call initial_distribution.f90

!  call propagation.f90


  call fin_parallel

end program main
