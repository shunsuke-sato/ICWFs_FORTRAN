subroutine parameters
  use global_variables
  implicit none
  integer :: ispec, ip, ipt

  call write_message('Start: parameters')
  call comm_barrier

! propagation parameters
!  propagation_method = NO_PROPAGATION
!  propagation_method = HERMITIAN_LIMIT
  propagation_method = INT_CWF

  num_trajectory = 4
  time_step = 0.01d0
  propagation_time = 40d0*fs
  num_time_step = aint(propagation_time/time_step)+1

! model parameters
  call shin_metiu_parameters
!  call ee_scattering_parameters

  call init_interaction(spec, num_species)

  num_total_particle = sum(spec(:)%nparticle)
  write(message(1),"(I7)")num_total_particle
  message(1) = "total number of particles = "//trim(message(1))

! table to convert from particle index to species index
  allocate(itable_particle2species(num_total_particle))
  allocate(itable_particle2particle(num_total_particle))
  
  ip = 0
  do ispec = 1, num_species
    do ipt = 1, spec(ispec)%nparticle
      ip = ip + 1
      itable_particle2species(ip)  = ispec
      itable_particle2particle(ip) = ipt

    end do
  end do


  call comm_barrier
  call write_message(message(1))
  call write_message('Finish: parameters')

end subroutine parameters

subroutine shin_metiu_parameters
  use global_variables
  implicit none
  real(8) :: lsize_elec(1), lsize_ion(1)
  integer :: nx_elec(1), nx_ion(1)
  real(8) :: mass_elec, mass_ion


  sampling_method = sampling_from_manybody_wf

! Parameters for Shin-Metiu model
  lsize_elec = 60d0
  nx_elec = 180
  mass_elec = 1d0

  lsize_ion  = 18d0
  nx_ion  = 150
  mass_ion = 1836d0


  num_species = 2
  allocate(spec(1:num_species))
! the first species  is electron
! the second species is ion

  call init_species(spec(1),1 ,1 , mass_elec, nx_elec, -0.5d0*lsize_elec, 0.5*lsize_elec, "electron")
  call init_species(spec(2),1 ,1 , mass_ion,  nx_ion , -0.5d0*lsize_ion , 0.5*lsize_ion,  "ion")


end subroutine shin_metiu_parameters

subroutine ee_scattering_parameters
  use global_variables
  implicit none
  real(8) :: lsize_elec(1)
  integer :: nx_elec(1)
  real(8) :: mass_elec


  sampling_method = sampling_from_manybody_wf

! Parameters for Shin-Metiu model
  lsize_elec = 400d0
  nx_elec = 2000
  mass_elec = 1d0

  num_species = 1
  allocate(spec(1:num_species))
! the first species  is electron

  call init_species(spec(1),1 ,2 , mass_elec, nx_elec, -0.5d0*lsize_elec, 0.5*lsize_elec, "electron")


end subroutine ee_scattering_parameters
