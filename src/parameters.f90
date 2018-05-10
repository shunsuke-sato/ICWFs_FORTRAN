subroutine parameters
  use global_variables
  implicit none
  integer :: ndim, nparticle, ngrid
  real(8) :: mass, x_ini, x_fin

  call write_message('Start: parameters')

  call shin_metiu_parameters

  call write_message('Finish: parameters')

end subroutine parameters

subroutine shin_metiu_parameters
  use global_variables
  implicit none
  real(8) :: lsize_elec(1), lsize_ion(1)
  integer :: nx_elec(1), nx_ion(1)
  real(8) :: mass_elec, mass_ion



! Parameters for Shin-Metiu model
  lsize_elec = 40d0
  nx_elec = 200
  mass_elec = 1d0

  lsize_ion  = 15d0
  nx_ion  = 150
  mass_ion = 1836d0


  num_species = 2
  allocate(spec(1:num_species))
! the first species  is electron
! the second species is ion

  call init_species(spec(1),1 ,1 , mass_elec, nx_elec, -0.5d0*lsize_elec, 0.5*lsize_elec)
  call init_species(spec(2),1 ,1 , mass_ion,  nx_ion , -0.5d0*lsize_ion , 0.5*lsize_ion)


end subroutine shin_metiu_parameters
