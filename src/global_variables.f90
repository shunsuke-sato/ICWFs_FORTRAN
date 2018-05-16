module global_variables
  use mpi
  use parallel
  use species
  use math_mod
  implicit none

! Physical constants
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0

! System parameters

  integer :: num_species
  type(species_t),allocatable :: spec(:)

! temporal arrays
  character(len=1024) :: message(64)


end module global_variables
