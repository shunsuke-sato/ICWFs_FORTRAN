subroutine propagation
  use global_variables
  implicit none

  select case(propagation_method)
  case(NO_PROPAGATION)
    return
  case(HERMITIAN_LIMIT)
    call propagation_hermitian
  end select

end subroutine propagation
