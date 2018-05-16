!! Please specify one-body as well as many-body potentials
!! one_body_pot_1,..., two_body_pot_1_2, ....
!-----------------------------------------------------------------------------------------
module interaction_mod
  use math_mod
  implicit none
  private

! Shin-Metiu model parameters
  real(8),parameter :: L = 19d0
  real(8),parameter :: Rf = 7d0
  real(8),parameter :: Rl = 4.4d0
  real(8),parameter :: Rr = 3.1d0

  public :: one_body_pot_1, &
            one_body_pot_2, &
            two_body_pot_1_2


contains
  function one_body_pot_1(x) result(pot)
    implicit none
    integer,parameter :: ndim = 1
    real(8),intent(in) :: x(1:ndim)
    real(8) :: pot
    
!    pot = 0d0
! electron potential in Shin-Metiu model
    pot = - erf_x(abs(x(1)-L/2d0)/Rr)/Rr &
          - erf_x(abs(x(1)+L/2d0)/Rl)/Rl

  end function one_body_pot_1

!-----------------------------------------------------------------------------------------
  function one_body_pot_2(x) result(pot)
    implicit none
    integer,parameter :: ndim = 1
    real(8),intent(in) :: x(1:ndim)
    real(8) :: pot
    
!    pot = 0d0
! ion potential in Shin-Metiu model
    pot = 1d0/abs(L/2d0-x(1)) + 1d0/abs(L/2d0+x(1))


  end function one_body_pot_2

!-----------------------------------------------------------------------------------------
  function two_body_pot_1_2(x1,x2) result(pot)
    implicit none
    integer,parameter :: ndim1 = 1, ndim2 = 1
    real(8),intent(in) :: x1(1:ndim1), x2(1:ndim2)
    real(8) :: pot
    
!    pot = 0d0
! electron-ion potential in Shin-Metiu model
    pot = -erf_x(abs(x1(1)-x2(1))/Rf)/Rf


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


