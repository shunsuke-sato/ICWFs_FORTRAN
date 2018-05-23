subroutine sampling
  use global_variables
  implicit none

  select case(sampling_method)
  case(sampling_from_manybody_wf)
    call sampling_from_manybody_wavefunction
  case default
    write(message(1),"(I7)")sampling_method
    message(1) = 'Fatal Error: sampling_method = '//trim(message(1))//'  is invalid.'
    call error_finalize(message(1))
  end select

end subroutine sampling

subroutine sampling_from_manybody_wavefunction
  use global_variables
  implicit none
  real(8) :: r1
  integer :: ix1, ix2, ix
  integer :: nx1, nx2
  integer :: ispec, ip, ipt

  select case(num_total_particle)
  case(2)

    nx1 = spec(itable_particle2species(1))%ngrid_tot
    nx2 = spec(itable_particle2species(2))%ngrid_tot

    call random_double(r1)
    ix = 0
    do ix2 = 2, nx2
      ix = 1 + (ix2-1)*nx1
      if(rho_cumulative_prob(ix) >=r1)exit
    end do
    ix2 = ix2 - 1

    do ix1 = 1, nx1
      ix = ix1 + (ix2-1)*nx1
      if(rho_cumulative_prob(ix) >=r1)exit
    end do

!! FIXME
!! Interpolation will be implemented later !!
! if dimensionalities of all species are one
!    if(product(spec(1:num_species)%ndim)==1)then

!    else

    ip = 0
    do ispec = 1, num_species
      do ipt = 1, spec(ispec)%nparticle
        ip = ip + 1
        select case(ip)
        case(1)
          spec(ispec)%r_particle(:,ipt) = spec(ispec)%x_ini(:) &
          + spec(ispec)%dx(:)*ix1
          
          spec(ispec)%zwfn(1:nx1,ipt) = zwfn_ini_2p(1:nx1,ix2)
        case(2)
          spec(ispec)%r_particle(:,ipt) = spec(ispec)%x_ini(:) &
          + spec(ispec)%dx(:)*ix2
          
          spec(ispec)%zwfn(1:nx2,ipt) = zwfn_ini_2p(ix1,1:nx2)
        case default 
          stop 'Error in sampling_from_manybody_wavefunction'
        end select
      end do
    end do


!    end if

  case default
    write(message(2),"(I7)")num_total_particle
    message(1) = 'Fatal Error: sampling_from_manybody_wavefunction is not implemented for'
    message(2) = '           : num_total_particle = '//trim(message(2))
    call error_finalize(message(1:2))
  end select


end subroutine sampling_from_manybody_wavefunction

