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
  end type species_icwf


  type trajectory_icwf

     type(species_icwf),allocatable :: spec(:)
     
  end type trajectory_icwf

  type(trajectory_icwf),allocatable :: traj(:)

  call initialize_icwfn_propagation
  call sampling_icwfn
  call initialize_icwfn_coefficient


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
                              if_density)
    implicit none
    logical,intent(in),optional :: if_overlap_matrix
    logical,intent(in),optional :: if_interaction_matrix
    logical,intent(in),optional :: if_density
    logical :: if_overlap_matrix_t
    logical :: if_interaction_matrix_t
    logical :: if_density_t
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

    if_density_t  = .false.
    if(present(if_density)) if_density_t = if_density


    allocate(spec_buf(num_species))
    ntraj_size = num_trajectory/comm_nproc_global
    if(mod(num_trajectory,comm_nproc_global) /= 0)ntraj_size = ntraj_size + 1
    do ispec = 1, num_species
      allocate(spec_buf(ispec)%zwfn_sbuf(spec(ispec)%ngrid_tot,spec(ispec)%nparticle,ntraj_size))
      allocate(spec_buf(ispec)%zwfn_rbuf(spec(ispec)%ngrid_tot,spec(ispec)%nparticle,ntraj_size))
    end do


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
                ztmp = sum(conjg(traj(itraj)%spec(ispec)%zwfn(:,ip)) &
                  *spec_buf(ispec)%zwfn_rbuf(:,ip,jtraj-ntraj_s_rbuf+1)) &
                  *spec(ispec)%dV
                zMm_sub_icwf(itraj,jtraj,ip_tot) = ztmp
                zMm_sub_icwf(jtraj,itraj,ip_tot) = conjg(ztmp)
              end do
            end do
            zMm_icwf(itraj,jtraj) = product(zMm_sub_icwf(itraj,jtraj,:))
            zMm_icwf(jtraj,itraj) = conjg(zMm_icwf(itraj,jtraj))

          end do
        end do
      end if
      
    end do

    if(if_overlap_matrix_t)call comm_bcast(zMm_icwf)



  end subroutine calc_icwf_matrix
end subroutine propagation_interacting_cwfn

