subroutine propagation_hermitian
  use global_variables
  implicit none
  integer :: itraj, it


  do itraj = 1, num_trajectory
    call sampling
    if(mod(itraj-1,comm_nproc_global) /= comm_id_global)cycle
    do it = 1, num_time_step



    end do
  end do


  contains


    subroutine dt_evolve_Runge_Kutta4
      implicit none

      type species_temporal
         complex(8),allocatable :: zwfn(:,:)  ! wavefunction
         real(8),allocatable :: r_particle(:,:) ! Position of particle
      end type species_temporal

      type(species_temporal) :: spec_tmp(num_species,4)
      



    end subroutine dt_evolve_Runge_Kutta4


end subroutine propagation_hermitian
