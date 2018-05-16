subroutine initial_state
  use global_variables
  implicit none

  call write_message('Start: initial_state')

  select case(sampling_method)
  case(sampling_from_manybody_wf)

    call init_manybody_wf

  case default
    write(message(1),"(I7)")sampling_method
    message(1) = 'Fatal Error: sampling_method = '//trim(message(1))//'  is invalid.'
    call error_finalize(message(1))
  end select


  call write_message('Finish: initial_state')

end subroutine initial_state
!-----------------------------------------------------------------------------------------
subroutine init_manybody_wf
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

end subroutine init_manybody_wf
