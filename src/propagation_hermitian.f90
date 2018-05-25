subroutine propagation_hermitian
  use global_variables
  implicit none
  integer,parameter :: nout_density = 100
  integer :: iout_density
  character(256) :: cfile_density
  logical,parameter :: if_write_traj = .false.
  character(256) :: cfile_traj

  type species_output
     real(8),allocatable :: rho(:,:)
  end type species_output
  type(species_output) :: spec_t(num_species)
  
  integer :: itraj, it, ispec, ip, ix
  real(8) :: ss

  do ispec = 1, num_species
    allocate(spec_t(ispec)%rho(spec(ispec)%ngrid_tot,0:nout_density+10))
    spec_t(ispec)%rho(:,:)=0d0
  end do



  Trajectory: do itraj = 1, num_trajectory
    if(if_root_global)write(*,*)"itraj =",itraj
    call sampling
    if(mod(itraj-1,comm_nproc_global) /= comm_id_global)cycle

    if(if_write_traj)then
      write(cfile_traj,"(I7.7)")itraj
      cfile_traj = trim(cfile_traj)//"_trajectory.out"
      open(31,file=cfile_traj)
      write(31,"(999e26.16e3)")0d0,spec(1)%r_particle(:,:),spec(2)%r_particle(:,:)
    end if



    iout_density = 0
    do ispec = 1,num_species
      do ip = 1,spec(ispec)%nparticle
        ss = sum(abs(spec(ispec)%zwfn(:,ip))**2)*spec(ispec)%dV
        spec_t(ispec)%rho(:,iout_density) = spec_t(ispec)%rho(:,iout_density) &
          + abs(spec(ispec)%zwfn(:,ip))**2/ss
      end do
    end do

    Time_propagation: do it = 1, num_time_step

      call dt_evolve_Runge_Kutta4_hermitian

      if( mod(it, max(num_time_step/nout_density,1) )== 0)then
        iout_density = iout_density + 1
        do ispec = 1,num_species
          do ip = 1,spec(ispec)%nparticle
            ss = sum(abs(spec(ispec)%zwfn(:,ip))**2)*spec(ispec)%dV
            spec_t(ispec)%rho(:,iout_density) = spec_t(ispec)%rho(:,iout_density) &
              + abs(spec(ispec)%zwfn(:,ip))**2/ss
          end do
        end do
      end if

      if(if_write_traj) write(31,"(999e26.16e3)")time_step*it,spec(1)%r_particle(:,:),spec(2)%r_particle(:,:)
    end do Time_propagation
    
    iout_density = iout_density + 1
    do ispec = 1,num_species
      do ip = 1,spec(ispec)%nparticle
        ss = sum(abs(spec(ispec)%zwfn(:,ip))**2)*spec(ispec)%dV
        spec_t(ispec)%rho(:,iout_density) = spec_t(ispec)%rho(:,iout_density) &
          + abs(spec(ispec)%zwfn(:,ip))**2/ss
      end do
    end do

    if(if_write_traj)then
      close(31)
    end if

  end do Trajectory

  do ispec = 1,num_species
    call comm_allreduce(spec_t(ispec)%rho(:,:))
    spec_t(ispec)%rho(:,:) = spec_t(ispec)%rho(:,:)/num_trajectory
  end do

  if(if_root_global)then
    do ispec = 1,num_species
      iout_density = -1 
      do it = 0, num_time_step+1
        if( mod(it, max(num_time_step/nout_density,1))== 0 .or. it == num_time_step+1)then
          iout_density = iout_density + 1
          write(cfile_density,"(I5.5)")iout_density
          cfile_density=trim(cfile_density)//"_"//trim(spec(ispec)%name)//"_rho.out"
          open(21,file=cfile_density)
          write(21,"(A,2x,e16.6e3)")"#time=",time_step*it
          do ix = 1, spec(ispec)%ngrid_tot
            write(21,"(4e16.6e3)")spec(ispec)%x(:,ix),spec_t(ispec)%rho(ix,iout_density)
          end do
          close(21)
        end if
      end do
    end do
  end if

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
                     -0.5d0/spec(ispec)%mass)
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
                     -0.5d0/spec(ispec)%mass)
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
                     -0.5d0/spec(ispec)%mass)
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
                     -0.5d0/spec(ispec)%mass)
      spec_t(ispec)%zwfn(:,ip,irk) = -zI*(&
          spec_t(ispec)%zwfn(:,ip,irk) &
        + spec_t(ispec)%veff(:,ip)*spec_t(ispec)%zwfn(:,ip,0) )

      spec_t(ispec)%r_p(:,ip,irk) = spec_t(ispec)%v_pt(:,ip)

    end do
  end do

!! Sum
  do ispec = 1, num_species
    spec(ispec)%zwfn(:,:) = spec(ispec)%zwfn(:,:) + time_step/6d0*( &
            spec_t(ispec)%zwfn(:,:, 1) &
     +2d0 * spec_t(ispec)%zwfn(:,:, 2) &
     +2d0 * spec_t(ispec)%zwfn(:,:, 3) &
          + spec_t(ispec)%zwfn(:,:, 4) )

    spec(ispec)%r_particle(:,:) = spec(ispec)%r_particle(:,:) + time_step/6d0*( &
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

