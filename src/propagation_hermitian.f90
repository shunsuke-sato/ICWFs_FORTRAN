subroutine propagation_hermitian
  use global_variables
  implicit none
  integer :: itraj, it


  do itraj = 1, num_trajectory
    call sampling
    if(mod(itraj-1,comm_nproc_global) /= comm_id_global)cycle
    do it = 1, num_time_step

      call dt_evolve_Runge_Kutta4_hermitian

    end do
  end do

end subroutine propagation_hermitian

subroutine dt_evolve_Runge_Kutta4_hermitian
  use global_variables
  implicit none
  integer :: ispec, ip, jspec, jp, ix
  type species_temporal
     complex(8),allocatable :: zwfn(:,:,:)  ! wavefunction
     real(8),allocatable :: veff(:,:)  ! wavefunction
     real(8),allocatable :: r_p(:,:,:) ! Position of particle
     real(8),allocatable :: v_pt(:,:) ! Position of particle
  end type species_temporal

  type(species_temporal) :: spec_t(num_species)
  integer :: irk

  do ispec = 1, num_species
    allocate(spec_t(ispec)%zwfn(spec(ispec)%ngrid_tot,spec(ispec)%nparticle, 0:4))
    allocate(spec_t(ispec)%veff(spec(ispec)%ngrid_tot,spec(ispec)%nparticle))
    allocate(spec_t(ispec)%r_p(spec(ispec)%ndim,spec(ispec)%nparticle, 0:4))
    allocate(spec_t(ispec)%v_pt(spec(ispec)%ndim,spec(ispec)%nparticle))
  end do

! RK1 ==========================================================================
  irk = 1
  do ispec = 1, num_species
    spec_t(ispec)%zwfn(:,:, 0) = spec(ispec)%zwfn(:,:)
    spec_t(ispec)%r_p(:,:, 0)  = spec(ispec)%r_particle(:,:)
  end do

  call calc_velocity_hermitian
  call calc_veff_hermitian

  do ispec = 1, num_species
    do ip = 1, spec(ispec)%nparticle
      call laplacian(spec_t(ispec)%zwfn(:,ip,0), &
                     spec_t(ispec)%zwfn(:,ip,irk), &
                     spec(ispec)%ngrid,          &
                     spec(ispec)%dx,             &
                     -0.5d0)
      spec_t(ispec)%zwfn(:,ip,irk) = -zI*(&
          spec_t(ispec)%zwfn(:,ip,irk) &
        + spec_t(ispec)%veff(:,ip)*spec_t(ispec)%zwfn(:,ip,0) )

      spec_t(ispec)%r_p(:,ip,irk) = spec_t(ispec)%v_pt(:,ip)

    end do
  end do

! RK2 ==========================================================================
  irk = 2
  do ispec = 1, num_species
    spec_t(ispec)%zwfn(:,:, 0) = spec(ispec)%zwfn(:,:) &
                               + 0.5d0*time_step*spec_t(ispec)%zwfn(:,:,1)
    spec_t(ispec)%r_p(:,:, 0)  = spec(ispec)%r_particle(:,:) &
                               + 0.5d0*time_step*spec_t(ispec)%r_p(:,:, 1)
  end do

  call calc_velocity_hermitian
  call calc_veff_hermitian

  do ispec = 1, num_species
    do ip = 1, spec(ispec)%nparticle
      call laplacian(spec_t(ispec)%zwfn(:,ip,0), &
                     spec_t(ispec)%zwfn(:,ip,irk), &
                     spec(ispec)%ngrid,          &
                     spec(ispec)%dx,             &
                     -0.5d0)
      spec_t(ispec)%zwfn(:,ip,irk) = -zI*(&
          spec_t(ispec)%zwfn(:,ip,irk) &
        + spec_t(ispec)%veff(:,ip)*spec_t(ispec)%zwfn(:,ip,0) )

      spec_t(ispec)%r_p(:,ip,irk) = spec_t(ispec)%v_pt(:,ip)

    end do
  end do


! RK3 ==========================================================================
  irk = 3
  do ispec = 1, num_species
    spec_t(ispec)%zwfn(:,:, 0) = spec(ispec)%zwfn(:,:) &
                               + 0.5d0*time_step*spec_t(ispec)%zwfn(:,:,2)
    spec_t(ispec)%r_p(:,:, 0)  = spec(ispec)%r_particle(:,:) &
                               + 0.5d0*time_step*spec_t(ispec)%r_p(:,:, 2)
  end do

  call calc_velocity_hermitian
  call calc_veff_hermitian

  do ispec = 1, num_species
    do ip = 1, spec(ispec)%nparticle
      call laplacian(spec_t(ispec)%zwfn(:,ip,0), &
                     spec_t(ispec)%zwfn(:,ip,irk), &
                     spec(ispec)%ngrid,          &
                     spec(ispec)%dx,             &
                     -0.5d0)
      spec_t(ispec)%zwfn(:,ip,irk) = -zI*(&
          spec_t(ispec)%zwfn(:,ip,irk) &
        + spec_t(ispec)%veff(:,ip)*spec_t(ispec)%zwfn(:,ip,0) )

      spec_t(ispec)%r_p(:,ip,irk) = spec_t(ispec)%v_pt(:,ip)

    end do
  end do


! RK4 ==========================================================================
  irk = 4

  do ispec = 1, num_species
    spec_t(ispec)%zwfn(:,:, 0) = spec(ispec)%zwfn(:,:) &
                               + time_step*spec_t(ispec)%zwfn(:,:,3)
    spec_t(ispec)%r_p(:,:, 0)  = spec(ispec)%r_particle(:,:) &
                               + time_step*spec_t(ispec)%r_p(:,:, 3)
  end do

  call calc_velocity_hermitian
  call calc_veff_hermitian

  do ispec = 1, num_species
    do ip = 1, spec(ispec)%nparticle
      call laplacian(spec_t(ispec)%zwfn(:,ip,0), &
                     spec_t(ispec)%zwfn(:,ip,irk), &
                     spec(ispec)%ngrid,          &
                     spec(ispec)%dx,             &
                     -0.5d0)
      spec_t(ispec)%zwfn(:,ip,irk) = -zI*(&
          spec_t(ispec)%zwfn(:,ip,irk) &
        + spec_t(ispec)%veff(:,ip)*spec_t(ispec)%zwfn(:,ip,0) )

      spec_t(ispec)%r_p(:,ip,irk) = spec_t(ispec)%v_pt(:,ip)

    end do
  end do

!! Sum
  do ispec = 1, num_species
    spec(ispec)%zwfn(:,:) = spec_t(ispec)%zwfn(:,:, 0) +time_step/6d0*( &
            spec_t(ispec)%zwfn(:,:, 1) &
     +2d0 * spec_t(ispec)%zwfn(:,:, 2) &
     +2d0 * spec_t(ispec)%zwfn(:,:, 3) &
          + spec_t(ispec)%zwfn(:,:, 4) )

    spec(ispec)%r_particle(:,:) = spec_t(ispec)%r_p(:,:,0) +time_step/6d0*( &
            spec_t(ispec)%r_p(:,:,1) &
     +2d0 * spec_t(ispec)%r_p(:,:,2) &
     +2d0 * spec_t(ispec)%r_p(:,:,3) &
          + spec_t(ispec)%r_p(:,:,4) )


  end do



  do ispec = 1, num_species
    deallocate(spec_t(ispec)%zwfn, spec_t(ispec)%r_p, spec_t(ispec)%v_pt)
  end do


contains
  subroutine calc_velocity_hermitian
    implicit none

! calc velocity
    do ispec = 1, num_species
      do ip = 1, spec(ispec)%nparticle
        call calc_velocity_from_cond_wf(spec(ispec),spec_t(ispec)%zwfn(:,ip,0), &
          spec_t(ispec)%r_p(:,ip,0),spec_t(ispec)%v_pt(:,ip))
      end do
    end do
    
  end subroutine calc_velocity_hermitian
    
    
  subroutine calc_veff_hermitian
    implicit none
    
    ! calc veff
    do ispec = 1, num_species
      do ip = 1, spec(ispec)%nparticle
        spec_t(ispec)%veff(:,ip) = spec(ispec)%v0(:)
        do jspec = 1, num_species
          do jp = 1, spec(jspec)%nparticle
            if( ispec == jspec .and. ip == jp)cycle
            if(ispec == 1 .and. jspec == 1)then
              do ix = 1, spec(ispec)%ngrid_tot
                spec_t(ispec)%veff(ix,ip) = spec_t(ispec)%veff(ix,ip) &
                  + two_body_pot_1(spec(ispec)%x(:,ix),spec_t(jspec)%r_p(:,jp,0))
              end do
            else if(ispec == 2 .and. jspec == 2)then          
              do ix = 1, spec(ispec)%ngrid_tot
                spec_t(ispec)%veff(ix,ip) = spec_t(ispec)%veff(ix,ip) &
                  + two_body_pot_2(spec(ispec)%x(:,ix),spec_t(jspec)%r_p(:,jp,0))
              end do
            else if(ispec == 1 .and. jspec == 2)then          
              do ix = 1, spec(ispec)%ngrid_tot
                spec_t(ispec)%veff(ix,ip) = spec_t(ispec)%veff(ix,ip) &
                  + two_body_pot_1_2(spec(ispec)%x(:,ix),spec_t(jspec)%r_p(:,jp,0))
              end do
            else if(ispec == 2 .and. jspec == 1)then          
              do ix = 1, spec(ispec)%ngrid_tot
                spec_t(ispec)%veff(ix,ip) = spec_t(ispec)%veff(ix,ip) &
                  + two_body_pot_1_2(spec_t(jspec)%r_p(:,jp,0),spec(ispec)%x(:,ix))
              end do
            else
              stop 'Error!!'
            end if
          end do
        end do
      end do
    end do
    
  end subroutine calc_veff_hermitian
  
  
end subroutine dt_evolve_Runge_Kutta4_hermitian

