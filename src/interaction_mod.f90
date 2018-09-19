!! Please specify one-body as well as many-body potentials
!! one_body_pot_1,..., two_body_pot_1_2, ....
!-----------------------------------------------------------------------------------------
module interaction_mod
  use math_mod
  use species
  implicit none
  private

! Shin-Metiu model parameters
  real(8),parameter :: L = 19d0
  real(8),parameter :: Rf = 5d0
  real(8),parameter :: Rl = 4.0d0
  real(8),parameter :: Rr = 3.1d0

  public :: init_interaction, &
            one_body_pot_1, &
            one_body_pot_2, &
            two_body_pot_1, &
            two_body_pot_2, &
            two_body_pot_1_2

  real(8),public,allocatable :: gf_two_body_pot_1(:,:)
  real(8),public,allocatable :: gf_two_body_pot_2(:,:)
  real(8),public,allocatable :: gf_two_body_pot_1_2(:,:)


contains
  subroutine init_interaction(spec, num_species)
    implicit none
    type(species_t),intent(in) :: spec(num_species)
    integer,intent(in) :: num_species
    integer :: nspec, ispec, jspec
    integer :: ngrid_i, ngrid_j
    integer :: ix, jx

    nspec = num_species


    do ispec = 1, nspec
      ngrid_i = spec(ispec)%ngrid_tot
      if(spec(ispec)%nparticle >1)then
        if(ispec == 1)then
          allocate(gf_two_body_pot_1(ngrid_i,ngrid_i))
          do ix = 1,ngrid_i
            do jx = 1,ngrid_i
              gf_two_body_pot_1(ix,jx) &
                = two_body_pot_1(spec(ispec)%x(:,ix),spec(ispec)%x(:,jx))
            end do
          end do
            
        else if(ispec == 2)then
          allocate(gf_two_body_pot_2(ngrid_i,ngrid_i))
          do ix = 1,ngrid_i
            do jx = 1,ngrid_i
              gf_two_body_pot_2(ix,jx) &
                = two_body_pot_2(spec(ispec)%x(:,ix),spec(ispec)%x(:,jx))
            end do
          end do
        else
          stop 'Error: in init_interaction; num_species>2.'
        end if
      end if

      do jspec = ispec+1, nspec
        ngrid_j = spec(jspec)%ngrid_tot
        if(ispec == 1 .and. jspec == 2)then
          allocate(gf_two_body_pot_1_2(ngrid_i,ngrid_j))
          do ix = 1,ngrid_i
            do jx = 1,ngrid_j
              gf_two_body_pot_1_2(ix,jx) &
                = two_body_pot_1_2(spec(ispec)%x(:,ix),spec(jspec)%x(:,jx))
            end do
          end do
        end if
      end do
    end do

  end subroutine init_interaction
!--------------------------------------------------------------------------------
  function one_body_pot_1(x) result(pot)
    implicit none
    integer,parameter :: ndim = 1
    real(8),intent(in) :: x(1:ndim)
    real(8) :: pot
    
!    pot = 0d0
! electron potential in Shin-Metiu model
    pot = - erf_x(abs(x(1)-L/2d0)/Rr)/Rr &
          - erf_x(abs(x(1)+L/2d0)/Rl)/Rl

! electron potential for electron-scattering problem

!    pot = -1d0/sqrt(2d0+x(1)**2)

! harmonic oscillator
!    pot = 0.5d0*x(1)**2
          

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
  function two_body_pot_1(x1,x2) result(pot)
    implicit none
    integer,parameter :: ndim1 = 1
    real(8),intent(in) :: x1(1:ndim1), x2(1:ndim1)
    real(8) :: pot
    
    pot = 0d0

! electron-electron iteraction for electron scattering problem
!    pot = 1d0/sqrt(2d0+(x1(1)-x2(1))**2)

  end function two_body_pot_1

!-----------------------------------------------------------------------------------------
  function two_body_pot_2(x1,x2) result(pot)
    implicit none
    integer,parameter :: ndim1 = 1
    real(8),intent(in) :: x1(1:ndim1), x2(1:ndim1)
    real(8) :: pot
    
    pot = 0d0

  end function two_body_pot_2

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
