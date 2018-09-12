subroutine initial_state
  use global_variables
  implicit none

  call write_message('Start: initial_state')
  call comm_barrier

  select case(sampling_method)
  case(sampling_from_manybody_wf)

    call allocate_manybody_wf
    call init_manybody_wf
    call init_cumulative_probability_distribution

  case default
    write(message(1),"(I7)")sampling_method
    message(1) = 'Fatal Error: sampling_method = '//trim(message(1))//'  is invalid.'
    call error_finalize(message(1))
  end select

  call comm_barrier
  call write_message('Finish: initial_state')

end subroutine initial_state
!-----------------------------------------------------------------------------------------
subroutine allocate_manybody_wf
  use global_variables
  implicit none
  integer,allocatable :: nx(:)
  integer :: ip, ispec

  allocate(nx(num_total_particle))
  do ip = 1, num_total_particle
    ispec = itable_particle2species(ip)
    nx(ip) = spec(ispec)%ngrid_tot
  end do
 
  select case(num_total_particle)
  case(2)
    allocate(zwfn_ini_2p(nx(1), nx(2)))
  case default
    write(message(2),"(I7)")num_total_particle
    message(1) = 'Fatal Error: init_manybody_wf is not implemented for'
    message(2) = '           : num_total_particle = '//trim(message(2))
    call error_finalize(message(1:2))
  end select

end subroutine allocate_manybody_wf
!-----------------------------------------------------------------------------------------
subroutine init_cumulative_probability_distribution
  use global_variables
  implicit none
  integer :: ip ,ispec
  integer :: ix1, ix2, ix
  integer :: nx1, nx2
  integer :: ngrid_tot
  real(8) :: dV, ss

  dV = 1d0
  do ip = 1, num_total_particle
    ispec = itable_particle2species(ip)
    dV = dV * product(spec(ispec)%dx(:))
  end do

  ngrid_tot = 1
  do ip = 1, num_total_particle
    ispec = itable_particle2species(ip)
    ngrid_tot = ngrid_tot * spec(ispec)%ngrid_tot
  end do
  allocate(rho_cumulative_prob(ngrid_tot))

  select case(num_total_particle)
  case(2)

    
    nx1 =  spec(itable_particle2species(1))%ngrid_tot
    nx2 =  spec(itable_particle2species(2))%ngrid_tot


! normalize
    ss = sum(abs(zwfn_ini_2p(:,:))**2)*dV
    zwfn_ini_2p = zwfn_ini_2p/sqrt(ss)

    ix = 0
    ss = 0d0
    do ix2 = 1, nx2
      do ix1 = 1, nx1
        ix = ix + 1
        ss = ss + abs(zwfn_ini_2p(ix1,ix2))**2
        rho_cumulative_prob(ix) = ss
      end do
    end do

    rho_cumulative_prob = rho_cumulative_prob/ss
    rho_cumulative_prob(ngrid_tot) = 1d0

  case default
    write(message(2),"(I7)")num_total_particle
    message(1) = 'Fatal Error: init_onebody_density is not implemented for'
    message(2) = '           : total_particle_num = '//trim(message(2))
    call error_finalize(message(1:2))
  end select

end subroutine init_cumulative_probability_distribution
!-----------------------------------------------------------------------------------------
subroutine init_manybody_wf
  use global_variables
  implicit none

! Shim-Metiu molel
  call init_manybody_wf_shin_metiu

! e-e scattering problem
!  call init_manybody_wf_ee_scattering

end subroutine init_manybody_wf
