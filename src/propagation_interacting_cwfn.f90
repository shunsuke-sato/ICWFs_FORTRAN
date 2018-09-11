subroutine propagation_interacting_cwfn
  use global_variables
  implicit none
  integer :: ntraj_start, ntraj_end

  complex(8),allocatable :: zC_icwf(:)
  complex(8),allocatable :: zMm_icwf(:,:)
  complex(8),allocatable :: zMm_sub_icwf(:,:,:)
  complex(8),allocatable :: zWm_icwf(:,:)

  type species_icwf
     complex(8),allocatable :: zwfn(:,:)  ! wavefunction
     real(8),allocatable :: r_p(:,:) ! Position of particle
     real(8),allocatable :: v_p(:,:) ! Velocity of particle
     real(8),allocatable :: veff(:,:) ! One-body effective potential
  end type species_icwf


  type trajectory_icwf

     type(species_icwf),allocatable :: spec(:)
     
  end type trajectory_icwf

  type(trajectory_icwf),allocatable :: traj(:)

  integer :: it

  call initialize_icwfn_propagation
  call sampling_icwfn
  call initialize_icwfn_coefficient


  Time_propagation: do it = 1, num_time_step
    call dt_evolve_Runge_Kutta4_icwfn

  end do Time_propagation


contains
  subroutine initialize_icwfn_propagation
    implicit none
    integer :: ntraj_ave, ntraj_remainder
    integer :: itraj, ispec
    
    ntraj_ave = num_trajectory/comm_nproc_global
    ntraj_remainder= mod(num_trajectory,comm_nproc_global)
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
        allocate(traj(itraj)%spec(ispec)%v_p(spec(ispec)%ndim,spec(ispec)%nparticle))
        allocate(traj(itraj)%spec(ispec)%veff(spec(ispec)%ngrid_tot,spec(ispec)%nparticle))
        
      end do
    end do
    
    
  end subroutine initialize_icwfn_propagation
  !-----------------------------------------------------------------------------------------
  subroutine sampling_icwfn
    implicit none
    integer :: itraj, ispec
    
    do itraj = 1, num_trajectory
      call sampling
      if(itraj >= ntraj_start .and. itraj <= ntraj_end)then
        do ispec = 1, num_species
          traj(itraj)%spec(ispec)%zwfn(:,:) = spec(ispec)%zwfn(:,:)
          traj(itraj)%spec(ispec)%r_p(:,:)  = spec(ispec)%r_particle(:,:)
        end do
      end if
    end do
    
    
  end subroutine sampling_icwfn
  

!-----------------------------------------------------------------------------------------
  subroutine initialize_icwfn_coefficient
    implicit none
    complex(8),allocatable :: zbvec(:)
    complex(8),allocatable :: zMm_pinv(:,:)
    integer,allocatable :: ispec_p(:),ip_p(:)
    integer :: ix1, ix2
    integer :: itraj
    complex(8) :: ztmp

    allocate(zC_icwf(num_trajectory))
    allocate(zMm_icwf(num_trajectory, num_trajectory))
    allocate(zMm_sub_icwf(num_trajectory, num_trajectory, num_total_particle))
    allocate(zWm_icwf(num_trajectory, num_trajectory))

    allocate(zMm_pinv(num_trajectory,num_trajectory))

   
    allocate(zbvec(num_trajectory))
    zbvec = 0d0
    
    select case(sampling_method)
    case(sampling_from_manybody_wf)
      select case(num_total_particle)
      case(2)

        allocate(ispec_p(2), ip_p(2))
        ispec_p(1) = itable_particle2species(1)
        ispec_p(2) = itable_particle2species(2)
        ip_p(1)    = itable_particle2particle(1)
        ip_p(2)    = itable_particle2particle(2)

        do itraj = ntraj_start, ntraj_end

          ztmp = 0d0
          do ix1 =1 , spec(ispec_p(1))%ngrid_tot
            do ix2 =1 , spec(ispec_p(2))%ngrid_tot
              ztmp = ztmp + conjg(&
                traj(itraj)%spec(ispec_p(1))%zwfn(ix1,ip_p(1)) &
               *traj(itraj)%spec(ispec_p(2))%zwfn(ix2,ip_p(2)) &
                )*zwfn_ini_2p(ix1,ix2)
            end do
          end do
          ztmp = ztmp * spec(ispec_p(1))%dV* spec(ispec_p(2))%dV

          zbvec(itraj) = ztmp

        end do

      case default
        write(message(2),"(I7)")num_total_particle
        message(1) = 'Fatal Error: sampling_from_manybody_wavefunction is not implemented for'
        message(2) = '           : num_total_particle = '//trim(message(2))
        call error_finalize(message(1:2))
      end select
    case default
      write(message(1),"(I7)")sampling_method
      message(1) = 'Fatal Error: sampling_method = '//trim(message(1))//'  is invalid.'
      call error_finalize(message(1))
    end select

    call comm_allreduce(zbvec)

    call calc_icwf_matrix(if_overlap_matrix=.true., if_interaction_matrix=.false.)
    call pseudo_inverse(zMm_icwf,zMm_pinv)

    zC_icwf(:) = matmul(zMm_pinv, zbvec)
    

  end subroutine initialize_icwfn_coefficient

  subroutine calc_icwf_matrix(if_overlap_matrix, &
                              if_interaction_matrix, &
                              if_onebody_density)
    implicit none
    logical,intent(in),optional :: if_overlap_matrix
    logical,intent(in),optional :: if_interaction_matrix
    logical,intent(in),optional :: if_onebody_density
    logical :: if_overlap_matrix_t
    logical :: if_interaction_matrix_t
    logical :: if_onebody_density_t
    integer :: icycle, ispec, ip, itraj, jtraj, ip_tot
    integer :: dest, source
    type species_buffer
       complex(8),allocatable :: zwfn_sbuf(:,:,:) ! send buffer
       complex(8),allocatable :: zwfn_rbuf(:,:,:) ! recieve buffer
    end type species_buffer
    type(species_buffer),allocatable :: spec_buf(:)
    integer :: ntraj_size
    integer :: ntraj_s_sbuf, ntraj_e_sbuf
    integer :: ntraj_s_rbuf, ntraj_e_rbuf
    complex(8) :: ztmp

    if_overlap_matrix_t      = .false.
    if(present(if_overlap_matrix))     if_overlap_matrix_t = if_overlap_matrix

    if_interaction_matrix_t  = .false.
    if(present(if_interaction_matrix)) if_interaction_matrix_t = if_interaction_matrix

    if_onebody_density_t  = .false.
    if(present(if_onebody_density)) if_onebody_density_t = if_onebody_density


    allocate(spec_buf(num_species))
    ntraj_size = num_trajectory/comm_nproc_global
    if(mod(num_trajectory,comm_nproc_global) /= 0)ntraj_size = ntraj_size + 1
    do ispec = 1, num_species
      allocate(spec_buf(ispec)%zwfn_sbuf(spec(ispec)%ngrid_tot,spec(ispec)%nparticle,ntraj_size))
      allocate(spec_buf(ispec)%zwfn_rbuf(spec(ispec)%ngrid_tot,spec(ispec)%nparticle,ntraj_size))
    end do

    if(if_overlap_matrix_t)then
      zMm_sub_icwf = 0d0
      zMm_icwf = 0d0
    end if

    if(if_interaction_matrix_t)then
      zWm_icwf = 0d0
    end if

    do icycle = 0,comm_nproc_global-1
      if(icycle == 0)then
        ntraj_s_rbuf = ntraj_start
        ntraj_e_rbuf = ntraj_end
        do itraj = ntraj_start, ntraj_end
          do ispec = 1, num_species
            spec_buf(ispec)%zwfn_rbuf(:,:,itraj-ntraj_start+1) &
              = traj(itraj)%spec(ispec)%zwfn(:,:)
          end do
        end do

      else

        ntraj_s_sbuf = ntraj_s_rbuf
        ntraj_e_sbuf = ntraj_e_rbuf
        dest   = mod(comm_id_global+icycle,                     comm_nproc_global)
        source = mod(comm_id_global-icycle + comm_nproc_global, comm_nproc_global)
        call comm_sendrecv(ntraj_s_sbuf, dest, 1, ntraj_s_rbuf, source, 1)
        call comm_sendrecv(ntraj_e_sbuf, dest, 1, ntraj_e_rbuf, source, 1)

        do ispec = 1, num_species
          spec_buf(ispec)%zwfn_sbuf(:,:,1:ntraj_size) &
            = spec_buf(ispec)%zwfn_rbuf(:,:,1:ntraj_size)
          call comm_sendrecv(spec_buf(ispec)%zwfn_sbuf, dest,   1, &
                             spec_buf(ispec)%zwfn_rbuf, source, 1)
        end do

      end if

      if(if_overlap_matrix_t)then
        do itraj = ntraj_start, ntraj_end
          do jtraj = ntraj_s_rbuf, ntraj_e_rbuf

            ip_tot = 0
            do ispec = 1, num_species
              do ip = 1, spec(ispec)%nparticle
                ip_tot = ip_tot + 1
                ztmp = sum(traj(itraj)%spec(ispec)%zwfn(:,ip) &
                  *conjg(spec_buf(ispec)%zwfn_rbuf(:,ip,jtraj-ntraj_s_rbuf+1))) &
                  *spec(ispec)%dV
                zMm_sub_icwf(jtraj,itraj,ip_tot) = ztmp
 !               zMm_sub_icwf(jtraj,itraj,ip_tot) = conjg(ztmp)
              end do
            end do
            zMm_icwf(jtraj,itraj) = product(zMm_sub_icwf(jtraj,itraj,:))
!            zMm_icwf(jtraj,itraj) = conjg(zMm_icwf(itraj,jtraj))

          end do
        end do
      end if

      if(if_interaction_matrix_t)then
        do itraj = ntraj_start, ntraj_end
          do jtraj = ntraj_s_rbuf, ntraj_e_rbuf

            ip_tot = 0
            do ispec = 1, num_species
              do ip = 1, spec(ispec)%nparticle
                ip_tot = ip_tot + 1

                jp_tot = 0
                do jspec = 1, num_species
                  do jp = 1, spec(jspec)%nparticle
                    jp_tot = jp_tot + 1
                    if(ip_tot <= jp_tot)cycle

                    select case(num_total_particle)
                    case(2)
                      if(ispec == 1 .and. jspec == 1)then
                        ztmp = 0d0
                        do ix1 = 1, spec(ispec)%ngrid_tot
                        do ix2 = 1, spec(jspec)%ngrid_tot
                          vint_t = two_body_pot_1(spec(ispec)%x(:,ix1),spec(jspec)%x(:,ix2))&
                            -two_body_pot_1(spec(jspec)%x(:,ix2),&
                            traj(itraj)%spec(ispec)%r_p(:,ip))
                          -two_body_pot_1(spec(jspec)%x(:,ix2),&
                            traj(itraj)%spec(jspec)%r_p(:,jp))

                          ztmp = ztmp & 
                            +traj(itraj)%spec(ispec)%zwfn(ix1,ip)&
                            *traj(itraj)%spec(jspec)%zwfn(ix2,jp)&
                            *conjg(&
                            spec_buf(ispec)%zwfn_rbuf(ix1,ip,jtraj-ntraj_s_rbuf+1)&
                            *spec_buf(jspec)%zwfn_rbuf(ix2,jp,jtraj-ntraj_s_rbuf+1))&
                            *vint_t

                        end do
                        end do
                      else if(ispec == 2 .and. jspec == 2)then
                        ztmp = 0d0
                        do ix1 = 1, spec(ispec)%ngrid_tot
                        do ix2 = 1, spec(jspec)%ngrid_tot
                          vint_t = two_body_pot_2(spec(ispec)%x(:,ix1),spec(jspec)%x(:,ix2))&
                            -two_body_pot_2(spec(jspec)%x(:,ix2),&
                            traj(itraj)%spec(ispec)%r_p(:,ip))
                          -two_body_pot_2(spec(jspec)%x(:,ix2),&
                            traj(itraj)%spec(jspec)%r_p(:,jp))

                          ztmp = ztmp & 
                            +traj(itraj)%spec(ispec)%zwfn(ix1,ip)&
                            *traj(itraj)%spec(jspec)%zwfn(ix2,jp)&
                            *conjg(&
                            spec_buf(ispec)%zwfn_rbuf(ix1,ip,jtraj-ntraj_s_rbuf+1)&
                            *spec_buf(jspec)%zwfn_rbuf(ix2,jp,jtraj-ntraj_s_rbuf+1))&
                            *vint_t

                        end do
                        end do
                      else if((ispec == 1 .and. jspec == 2) &
                        .or. (ispec == 2 .and. jspec == 1))then
                        ztmp = 0d0
                        do ix1 = 1, spec(ispec)%ngrid_tot
                        do ix2 = 1, spec(jspec)%ngrid_tot
                      vint_t = two_body_pot_1_2(spec(ispec)%x(:,ix1),spec(jspec)%x(:,ix2))&
                        -two_body_pot_1_2(spec(jspec)%x(:,ix2),&
                        traj(itraj)%spec(ispec)%r_p(:,ip))
                      -two_body_pot_1_2(spec(jspec)%x(:,ix2),&
                        traj(itraj)%spec(jspec)%r_p(:,jp))

                      ztmp = ztmp & 
                        +traj(itraj)%spec(ispec)%zwfn(ix1,ip)&
                        *traj(itraj)%spec(jspec)%zwfn(ix2,jp)&
                        *conjg(&
                        spec_buf(ispec)%zwfn_rbuf(ix1,ip,jtraj-ntraj_s_rbuf+1)&
                        *spec_buf(jspec)%zwfn_rbuf(ix2,jp,jtraj-ntraj_s_rbuf+1))&
                        *vint_t
                      
                        end do
                        end do
                      end if

                      ztmp = ztmp*spec(ispec)%dV*spec(jspec)%dV
                    case default
                      write(message(2),"(I7)")num_total_particle
                      message(1) = 'Fatal Error: calc_icwf_matrix is not implemented for'
                      message(2) = '           : num_total_particle = '//trim(message(2))
                      call error_finalize(message(1:2))
                    end select

                  end do
                end do
              end do
            end do

            zWm_icwf(jtraj,itraj) = ztmp

          end do
        end do
      end if
      
      if(if_onebody_density_t)then
      end if


    end do

    if(if_overlap_matrix_t)then
      call comm_allreduce(zMm_sub_icwf)
      call comm_allreduce(zMm_icwf)
    end if



  end subroutine calc_icwf_matrix


  subroutine dt_evolve_Runge_Kutta4_icwfn
    implicit none
    integer :: irk

    type(trajectory_icwf),allocatable :: traj_rk(:,:)
    complex(8) :: zC_rk(num_trajectory,-1:4)
    complex(8) :: zMm_pinv(num_trajectory,num_trajectory)

    allocate(traj_rk(ntraj_start:ntraj_end,-1:4))

    do irk = -1, 4
      do itraj = ntraj_start, ntraj_end
        allocate(traj_rk(itraj,irk)%spec(num_species))

        do ispec = 1, num_species
        
          allocate(traj_rk(itraj,irk)%spec(ispec)%zwfn(spec(ispec)%ngrid_tot,spec(ispec)%nparticle))
          allocate(traj_rk(itraj,irk)%spec(ispec)%r_p(spec(ispec)%ndim,spec(ispec)%nparticle))
          
        end do
      end do
    end do

! RK -1 (original wavefunction)
    zC_rk(:,-1) = zC_icwf(:)
    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        traj_rk(itraj,-1)%spec(ispec)%zwfn = traj(itraj)%spec(ispec)%zwfn
        traj_rk(itraj,-1)%spec(ispec)%r_p = traj(itraj)%spec(ispec)%r_p
      end do
    end do


! RK1 =========================================================================
    irk = 1
    call calc_icwf_matrix(if_overlap_matrix=.true., &
                              if_interaction_matrix=.true.)
    call pseudo_inverse(zMm_icwf,zMm_pinv)
    zC_rk(:,irk) = -zI* matmul( zMm_pinv, matmul(zWm_icwf, zC_icwf))

    call calc_velocity_icwfn
    call calc_veff_icwfn

    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        do ip = 1, spec(ispec)%nparticle
          call laplacian(traj(itraj)%spec(ispec)%zwfn(:,ip), &
                         traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip), &
                         spec(ispec)%ngrid,          &
                         spec(ispec)%dx,             &
                         -0.5d0/spec(ispec)%mass)
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) = -zI*(&
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) &
        + traj(itraj)%spec(ispec)%veff(:,ip)&
        * traj(itraj)%spec(ispec)%zwfn(:,ip) )

          traj_rk(itraj,irk)%spec(ispec)%r_p(:,ip) &
            = spec_t(ispec)%traj(itraj)%spec(ispec)%v_p(:,ip)

        end do
      end do
    end do


! RK2 =========================================================================
    irk = 2
    zC_icwf(:) = zC_rk(:,-1) + 0.5d0*time_step*zC_rk(:,1)
    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        traj(itraj)%spec(ispec)%zwfn = traj_rk(itraj,-1)%spec(ispec)%zwfn &
          +0.5d0*time_step*traj_rk(itraj,1)%spec(ispec)%zwfn

        traj(itraj)%spec(ispec)%r_p = traj_rk(itraj,-1)%spec(ispec)%r_p &
          +0.5d0*time_step*traj_rk(itraj,1)%spec(ispec)%r_p

      end do
    end do

    call calc_icwf_matrix(if_overlap_matrix=.true., &
                              if_interaction_matrix=.true.)
    call pseudo_inverse(zMm_icwf,zMm_pinv)
    zC_rk(:,irk) = -zI* matmul( zMm_pinv, matmul(zWm_icwf, zC_icwf))

    call calc_velocity_icwfn
    call calc_veff_icwfn

    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        do ip = 1, spec(ispec)%nparticle
          call laplacian(traj(itraj)%spec(ispec)%zwfn(:,ip), &
                         traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip), &
                         spec(ispec)%ngrid,          &
                         spec(ispec)%dx,             &
                         -0.5d0/spec(ispec)%mass)
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) = -zI*(&
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) &
        + traj(itraj)%spec(ispec)%veff(:,ip)&
        * traj(itraj)%spec(ispec)%zwfn(:,ip) )

          traj_rk(itraj,irk)%spec(ispec)%r_p(:,ip) &
            = spec_t(ispec)%traj(itraj)%spec(ispec)%v_p(:,ip)

        end do
      end do
    end do

! RK3 =========================================================================
    irk = 3
    zC_icwf(:) = zC_rk(:,-1) + 0.5d0*time_step*zC_rk(:,2)
    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        traj(itraj)%spec(ispec)%zwfn = traj_rk(itraj,-1)%spec(ispec)%zwfn &
          +0.5d0*time_step*traj_rk(itraj,2)%spec(ispec)%zwfn

        traj(itraj)%spec(ispec)%r_p = traj_rk(itraj,-1)%spec(ispec)%r_p &
          +0.5d0*time_step*traj_rk(itraj,2)%spec(ispec)%r_p

      end do
    end do

    call calc_icwf_matrix(if_overlap_matrix=.true., &
                              if_interaction_matrix=.true.)
    call pseudo_inverse(zMm_icwf,zMm_pinv)
    zC_rk(:,irk) = -zI* matmul( zMm_pinv, matmul(zWm_icwf, zC_icwf))

    call calc_velocity_icwfn
    call calc_veff_icwfn

    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        do ip = 1, spec(ispec)%nparticle
          call laplacian(traj(itraj)%spec(ispec)%zwfn(:,ip), &
                         traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip), &
                         spec(ispec)%ngrid,          &
                         spec(ispec)%dx,             &
                         -0.5d0/spec(ispec)%mass)
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) = -zI*(&
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) &
        + traj(itraj)%spec(ispec)%veff(:,ip)&
        * traj(itraj)%spec(ispec)%zwfn(:,ip) )

          traj_rk(itraj,irk)%spec(ispec)%r_p(:,ip) &
            = spec_t(ispec)%traj(itraj)%spec(ispec)%v_p(:,ip)

        end do
      end do
    end do

! RK4 =========================================================================
    irk = 4
    zC_icwf(:) = zC_rk(:,-1) + time_step*zC_rk(:,3)
    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        traj(itraj)%spec(ispec)%zwfn = traj_rk(itraj,-1)%spec(ispec)%zwfn &
          +time_step*traj_rk(itraj,2)%spec(ispec)%zwfn

        traj(itraj)%spec(ispec)%r_p = traj_rk(itraj,-1)%spec(ispec)%r_p &
          +time_step*traj_rk(itraj,2)%spec(ispec)%r_p

      end do
    end do

    call calc_icwf_matrix(if_overlap_matrix=.true., &
                              if_interaction_matrix=.true.)
    call pseudo_inverse(zMm_icwf,zMm_pinv)
    zC_rk(:,irk) = -zI* matmul( zMm_pinv, matmul(zWm_icwf, zC_icwf))

    call calc_velocity_icwfn
    call calc_veff_icwfn

    do itraj = ntraj_start, ntraj_end
      do ispec = 1, num_species
        do ip = 1, spec(ispec)%nparticle
          call laplacian(traj(itraj)%spec(ispec)%zwfn(:,ip), &
                         traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip), &
                         spec(ispec)%ngrid,          &
                         spec(ispec)%dx,             &
                         -0.5d0/spec(ispec)%mass)
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) = -zI*(&
          traj_rk(itraj,irk)%spec(ispec)%zwfn(:,ip) &
        + traj(itraj)%spec(ispec)%veff(:,ip)&
        * traj(itraj)%spec(ispec)%zwfn(:,ip) )

          traj_rk(itraj,irk)%spec(ispec)%r_p(:,ip) &
            = spec_t(ispec)%traj(itraj)%spec(ispec)%v_p(:,ip)

        end do
      end do
    end do

!-----------------------------------------------------------


  contains
    subroutine calc_velocity_icwfn
      implicit none
      integer :: itraj, ispec, ip

      do itraj = ntraj_start, ntraj_end
        do ispec = 1, num_species
          do ip = 1, spec(ispec)%nparticle

            call calc_velocity_from_cond_wf(spec(ispec),traj(itraj)%spec(ispec)%zwfn(:,ip), &
              traj(itraj)%spec(ispec)%r_p(:,ip),traj(itraj)%spec(ispec)%v_pt(:,ip))

          end do
        end do
      end do

    end subroutine calc_velocity_icwfn

    subroutine calc_veff_icwfn
    implicit none
    integer :: itraj, ispec, ip, jspec, jp    
    
    ! calc veff
    do itraj = ntraj_start, ntraj_end

    do ispec = 1, num_species
      do ip = 1, spec(ispec)%nparticle
        traj(itraj)%spec(ispec)%veff(:,ip) = spec(ispec)%v0(:)
        do jspec = 1, num_species
          do jp = 1, spec(jspec)%nparticle
            if( ispec == jspec .and. ip == jp)cycle
            if(ispec == 1 .and. jspec == 1)then
              do ix = 1, spec(ispec)%ngrid_tot
                traj(itraj)%spec(ispec)%veff(ix,ip) = traj(itraj)%spec(ispec)%veff(ix,ip) &
                  + two_body_pot_1(spec(ispec)%x(:,ix),traj(itraj)%spec(ispec)%r_p(:,jp))
              end do
            else if(ispec == 2 .and. jspec == 2)then          
              do ix = 1, spec(ispec)%ngrid_tot
                traj(itraj)%spec(ispec)%veff(ix,ip) = traj(itraj)%spec(ispec)%veff(ix,ip) &
                  + two_body_pot_2(spec(ispec)%x(:,ix),traj(itraj)%spec(ispec)%r_p(:,jp))
              end do
            else if(ispec == 1 .and. jspec == 2)then          
              do ix = 1, spec(ispec)%ngrid_tot
                traj(itraj)%spec(ispec)%veff(ix,ip) = traj(itraj)%spec(ispec)%veff(ix,ip) &
                  + two_body_pot_1_2(spec(ispec)%x(:,ix),traj(itraj)%spec(ispec)%r_p(:,jp))
              end do
            else if(ispec == 2 .and. jspec == 1)then          
              do ix = 1, spec(ispec)%ngrid_tot
                traj(itraj)%spec(ispec)%veff(ix,ip) = traj(itraj)%spec(ispec)%veff(ix,ip) &
                  + two_body_pot_1_2(traj(itraj)%spec(ispec)%r_p(:,jp),spec(ispec)%x(:,ix))
              end do
            else
              stop 'Error!!'
            end if
          end do
        end do
      end do
    end do

    end do
    
  end subroutine calc_veff_icwfn

  end subroutine dt_evolve_Runge_Kutta4_icwfn
end subroutine propagation_interacting_cwfn


