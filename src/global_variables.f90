module global_variables
  use mpi
  use parallel
  use species
  use math_mod
  use interaction_mod
  use finite_difference_mod
  use random_number_mod
  implicit none

! Physical constants
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0

! System parameters
  integer :: num_species
  type(species_t),allocatable :: spec(:)

  integer :: num_total_particle
  integer,allocatable :: itable_particle2species(:)

  integer :: sampling_method
  integer,parameter :: sampling_from_manybody_wf = 0

! 
  complex(8),allocatable :: zwfn_ini_2p(:,:),zwfn_ini_3p(:,:,:)
  real(8),allocatable :: rho_cumulative_prob(:)

end module global_variables
