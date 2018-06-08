subroutine propagation_interacting_cwfn
  use global_variables
  implicit none
  integer :: ntraj_start, ntraj_end

  type species_icwf
     complex(8),allocatable :: zwfn(:,:)  ! wavefunction
    real(8),allocatable :: r_p(:,:) ! Position of particle
  end type species_icwf


  type trajectory_icwf

     type(species_icwf),allocatable :: spec(:)
     
  end type trajectory_icwf

  type(trajectory_icwf),allocatable :: traj(:)



  contains
    subroutine initialize_icwfn_propagation
      implicit none
      integer :: ntraj_ave, ntraj_remainder
      integer :: itraj, ispec

      ntraj_average = ntraj/comm_nproc_global
      ntraj_start= mod(ntraj,comm_nproc_global)
      if(comm_id_global+1 <= ntraj_remainder)then
        ntraj_start = 1 + comm_id_global*(ntraj_ave+1)
        ntraj_end= ntraj_start + (ntraj_ave+1) -1
      else
        ntraj_start = 1 + ntraj_remainder*(ntraj_ave + 1)  &
        + ntraj_ave*(comm_id_global-ntraj_remainder)
        ntraj_end = ntraj_start + ntraj_ave -1
      end if

      allocate(traj(ntraj_start:ntraj_end))

      do itraj = ntraj_start, ntraj_end
        allocate(traj(itraj)%spec(num_species))
      end do

      do itraj = ntraj_start, ntraj_end
        do ispec = 1, num_species

          allocate(traj(itraj)%spec(ispec)%zwfn(spec(ispec)%ngrid_tot,spec(ispec)%nparticle))
          allocate(traj(itraj)%spec(ispec)%r_p(spec(ispec)%ndim,spec(ispec)%nparticle))

        end do
      end do
      
      
    end subroutine initialize_icwfn_propagation
!-----------------------------------------------------------------------------------------
    subroutine sampling_icwfn
      implicit none
      integer :: itraj

      do itraj = 1, ntraj
        call sampling
        if(itraj >= ntraj_start .and. itraj <= ntraj_end)
        traj(itraj)%spec(:)%zwfn(:,:) = spec(:)%zwfn(:,:)
        traj(itraj)%spec(:)%r_p(:,:) = spec(:)%r_particle(:,:)
        
      end do

      
    end subroutine sampling_icwfn



end subroutine propagation_interacting_cwfn
