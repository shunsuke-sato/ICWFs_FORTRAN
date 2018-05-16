subroutine initial_state
  use global_variables
  implicit none

  call write_message('Start: initial_state')

  select case(sampling_method)
  case(sampling_from_manybody_wf)
  case default
    write(message(1),"(I7)")sampling_method
    message(1) = 'Fatal Error: sampling_method = '//trim(message(1))//'  is invalid.'
    call error_finalize(message(1))
  end select


  call write_message('Finish: initial_state')

end subroutine initial_state

