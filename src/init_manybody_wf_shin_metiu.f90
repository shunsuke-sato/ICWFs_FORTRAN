subroutine init_manybody_wf_shin_metiu
  use global_variables
  implicit none
  integer,parameter :: Rion = -4d0
  integer,parameter :: norb = 2
  real(8),parameter :: res_epsilon = 1d-8
  real(8),parameter :: dt_imag = 0.01d0
  real(8),allocatable :: phi(:,:)
  real(8),allocatable :: veff(:),phi_t(:),hphi_t(:)
  integer,parameter :: nt_imag_max = 100000
  real(8) :: residual(norb), eps(norb)
  integer :: ix_e, ix_i
  integer :: ix, iorb, jorb, nx, it
  real(8) :: xe(1), xi(1), dx
  real(8) :: norm
  real(8) :: r1
  
  nx = spec(1)%ngrid(1)
  dx = spec(1)%dx(1)

  allocate(phi(nx, norb))
  allocate(phi_t(nx),hphi_t(nx))
  allocate(veff(nx))

! initialize wavefunction
  do iorb = 1,norb
    do ix = 1, nx
      
      call random_number(r1)
      phi(ix,iorb) = r1

    end do
  end do
  call gram_schmidt

! setting effective potential
  do ix = 1,nx
    xe(1) = spec(1)%x_ini(1) + spec(1)%dx(1)*ix_e
    xi(1) = Rion
    veff(ix) = two_body_pot_1_2(xe,xi)
  end do
  veff(:) = veff(:) + spec(1)%v0(:)

! imaginary time propagation
  do it = 1, nt_imag_max
    
    do iorb = 1, norb
      phi_t(:) = phi(:,iorb)
      call hphi(phi_t, hphi_t)
      phi(:,iorb) = phi(:,iorb) -dt_imag*hphi_t(:)

    end do
    call gram_schmidt

  end do



!! FIXME:
!! Initialization of many-body wavefunction for Shin-Metiue will be implemented later.
!! The gaussian is employed temporary.

  do ix_i = 1, spec(2)%ngrid(1)
    xi = spec(2)%x_ini(1) + spec(2)%dx(1)*ix_i
    do ix_e = 1, spec(1)%ngrid(1)
      xe = spec(1)%x_ini(1) + spec(1)%dx(1)*ix_e
      
      zwfn_ini_2p(ix_e,ix_i) = exp(-0.5d0*abs(xe(1))**2 -0.5d0*abs(xi(1))**2)

    end do
  end do

  norm = sum(abs(zwfn_ini_2p)**2)*spec(1)%dx(1)*spec(2)%dx(1)

  zwfn_ini_2p = zwfn_ini_2p/sqrt(norm)

contains
  subroutine gram_schmidt
    implicit none

    do iorb = 1, norb

      norm = 0d0
      do jorb = 1, iorb-1
        norm = sum(phi(:,iorb)*phi(:,jorb))*dx
        phi(:,iorb) = phi(:,iorb) - norm*phi(:,jorb)

      end do
      norm = sum(phi(:,iorb)**2)*dx
      phi(:,iorb) = phi(:,iorb)/sqrt(norm)
    end do

  end subroutine gram_schmidt

  subroutine hphi(f_in,f_out)
    implicit none
    real(8) :: f_in(nx), f_out(nx)

    call laplacian(f_in,f_out,spec(1)%ngrid,spec(1)%dx,-0.5d0)
    f_out = f_out + f_in*veff

  end subroutine hphi
  

end subroutine init_manybody_wf_shin_metiu
