subroutine propagation
  use global_variables
  implicit none

  select case(propagation_method)
  case(NO_PROPAGATION)
    return
  case(HERMITIAN_LIMIT)
    call propagation_hermitian
  case default
    call error_finalize('Error: Invalid propagation scheme.')
  end select

end subroutine propagation
