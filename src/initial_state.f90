subroutine initial_state
  use global_variables
  implicit none

  call write_message('Start: initial_state')

  select case(sampling_method)
  case(sampling_from_manybody_wf)

    call allocate_manybody_wf
    call init_manybody_wf
    call init_onebody_density

  case default
    write(message(1),"(I7)")sampling_method
    message(1) = 'Fatal Error: sampling_method = '//trim(message(1))//'  is invalid.'
    call error_finalize(message(1))
  end select


  call write_message('Finish: initial_state')

end subroutine initial_state
!-----------------------------------------------------------------------------------------
subroutine allocate_manybody_wf
  use global_variables
  implicit none
  integer :: nx(3),ispec,iparticle,itot_particle

  select case(total_particle_num)
  case(2)
    itot_particle = 0
    do ispec = 1, num_species
      do iparticle = 1, spec(ispec)%nparticle
        itot_particle = itot_particle + 1
        nx(itot_particle) = spec(ispec)%ngrid_tot
      end do
    end do

    allocate(zwfn_ini_2p(nx(1), nx(2)))

  case default
    write(message(2),"(I7)")total_particle_num
    message(1) = 'Fatal Error: init_manybody_wf is not implemented for'
    message(2) = '           : total_particle_num = '//trim(message(2))
    call error_finalize(message(1:2))
  end select

end subroutine allocate_manybody_wf
!-----------------------------------------------------------------------------------------
subroutine init_onebody_density
  use global_variables
  implicit none
  integer :: ix
  real(8) :: dV

  select case(total_particle_num)
  case(2)

! compute density for the first spieces
    if(num_species == 1)then
      dV = product(spec(1)%dx(:))
    else if(num_species == 2)then
      dV = product(spec(2)%dx(:))
    end if

    do ix = 1,spec(1)%ngrid_tot
      spec(1)%rho_ini(ix) = sum(abs(zwfn_ini_2p(ix,:))**2)*dV
    end do

    if(num_species == 2)then
      dV = product(spec(1)%dx(:))
      do ix = 1,spec(2)%ngrid_tot
        spec(2)%rho_ini(ix) = sum(abs(zwfn_ini_2p(:,ix))**2)*dV
      end do
    end if

  case default
    write(message(2),"(I7)")total_particle_num
    message(1) = 'Fatal Error: init_onebody_density is not implemented for'
    message(2) = '           : total_particle_num = '//trim(message(2))
    call error_finalize(message(1:2))
  end select

end subroutine init_onebody_density
!-----------------------------------------------------------------------------------------
subroutine init_manybody_wf
  use global_variables
  implicit none

! Shim-Metiu molel
  call init_manybody_wf_shin_metiu

end subroutine init_manybody_wf
