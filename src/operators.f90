!! Please specify one-body as well as many-body potentials
!! one_body_pot_1,..., two_body_pot_1_2, ....
!-----------------------------------------------------------------------------------------
module interaction_mod
contains
  function one_body_pot_1(x) result(pot)
    implicit none
    integer,parameter :: ndim = 1
    real(8),intent(in) :: x(1:ndim)
    real(8) :: pot
    
    pot = 0d0

  end function one_body_pot_1

!-----------------------------------------------------------------------------------------
  function one_body_pot_2(x) result(pot)
    implicit none
    integer,parameter :: ndim = 1
    real(8),intent(in) :: x(1:ndim)
    real(8) :: pot
    
    pot = 0d0

  end function one_body_pot_2

!-----------------------------------------------------------------------------------------
  function two_body_pot_1_2(x1,x2) result(pot)
    implicit none
    integer,parameter :: ndim1 = 1, ndim2 = 1
    real(8),intent(in) :: x1(1:ndim1), x2(1:ndim2)
    real(8) :: pot
    
    pot = 0d0

  end function two_body_pot_1_2
end module interaction_mod
!-----------------------------------------------------------------------------------------
! This is the start of the subroutine
subroutine operators
  use global_variables
  use interaction_mod
  implicit none
  integer :: ispec,i

  call write_message('Start: operators')

! Initialize one-body potential

  do ispec = 1, num_species
    select case(ispec)
    case(1)
      do i = 1, spec(ispec)%ngrid_tot
        spec(ispec)%v0(i) = one_body_pot_1(spec(ispec)%x(:,i))
      end do
    case(2)
      do i = 1, spec(ispec)%ngrid_tot
        spec(ispec)%v0(i) = one_body_pot_2(spec(ispec)%x(:,i))
      end do
    end select
  end do



  call write_message('Finish: operators')

end subroutine operators


