subroutine propagation_interacting_cwfn
  use global_variables
  implicit none

  type species_icwf
     complex(8),allocatable :: zwfn(:,:,:)  ! wavefunction
  end type species_icwf

end subroutine propagation_interacting_cwfn
