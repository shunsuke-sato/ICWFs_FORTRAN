! MAIN FILE
program main
  use global_variables
  implicit none

  call init_parallel

  write(*,*)"Hello world!!"
  write(*,*)"comm_id_global = ",comm_id_global
  write(*,*)"Pi = ",pi

  call fin_parallel

end program main
