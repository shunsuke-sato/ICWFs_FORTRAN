subroutine propagation_hermitian
  use global_variables
  implicit none
  integer :: itraj, it


  do itraj = 1, num_trajectory
    call sampling
    if(mod(itraj-1,comm_nproc_global) /= comm_id_global)cycle
    do it = 1, num_time_step



    end do
  end do




end subroutine propagation_hermitian
