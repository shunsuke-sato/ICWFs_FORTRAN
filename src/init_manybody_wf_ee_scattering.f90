subroutine init_manybody_wf_ee_scattering
  use global_variables
  implicit none
  real(8),allocatable :: xvec(:),gvec(:),pvec(:)
  real(8),allocatable :: hxvec(:),hpvec(:)
  real(8) :: gg0, gg, xx, xhx, pp, php, xp, xhp
  real(8) :: alpha, beta, aa, bb, cc, ss, lambda
  real(8) :: dx
  integer :: nx

  nx = spec(1)%ngrid(1)
  dx = spec(1)%dx(1)


end subroutine init_manybody_wf_ee_scattering
