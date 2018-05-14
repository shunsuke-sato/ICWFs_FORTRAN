module species
  implicit none

  type species_t
    integer :: ndim                      ! Dimensionality
    integer :: nparticle              ! Number of particles
    real(8) :: mass                      ! Mass of particle
    integer,allocatable :: ngrid(:)   ! Number of grid points
    integer :: ngrid_tot   ! Number of grid points
! x = x_ini + i*dx (i+1,n)
    real(8),allocatable :: x_ini(:)    ! Initial position of particle coordinates
    real(8),allocatable :: x_fin(:)    ! Final position of particle coordinates
    real(8),allocatable :: dx(:)    ! Spacing of grid points
    real(8),allocatable :: x(:,:)    ! Position of real-space grids
    complex(8),allocatable :: zwfn(:,:)  ! wavefunction
    real(8),allocatable :: r_particle(:,:) ! Position of particle
    real(8),allocatable :: v0(:) ! One-body potential

  end type species_t

  contains
    subroutine init_species(spec,ndim, nparticle, mass, ngrid, x_ini, x_fin)
      type(species_t),intent(inout)  :: spec
      integer,intent(in) :: ndim, nparticle, ngrid(1:ndim)
      real(8),intent(in) :: mass, x_ini(1:ndim), x_fin(1:ndim)
      integer :: ngrid_tot
      integer :: i,j,k,l

      spec%ndim      = ndim
      spec%nparticle = nparticle
      spec%mass      = mass

      allocate(spec%ngrid(1:ndim))
      allocate(spec%x_ini(1:ndim))
      allocate(spec%x_fin(1:ndim))
      allocate(spec%dx(1:ndim))
      spec%ngrid(1:ndim) = ngrid(1:ndim)
      spec%x_ini(1:ndim) = x_ini
      spec%x_fin(1:ndim) = x_fin

      ngrid_tot = product(ngrid(1:ndim))
      spec%ngrid_tot = ngrid_tot

      allocate(spec%x(ndim,ngrid_tot))
      allocate(spec%zwfn(ngrid_tot,nparticle))
      allocate(spec%v0(ngrid_tot))
      allocate(spec%r_particle(ndim,nparticle))


      spec%dx(1:ndim) = (x_fin(1:ndim) - x_ini(1:ndim))/(ngrid(1:ndim)+1)

! Initialize position of grids
      select case(ndim)
      case(1)
        l = 0
        do i = 1, ngrid(1)
          l = l + 1
          spec%x(1,l) = spec%x_ini(1) + spec%dx(1)*i
        end do
        
      case(2)
        l = 0
        do i = 1, ngrid(1)
          do j = 1, ngrid(2)
            l = l + 1
            spec%x(1,l) = spec%x_ini(1) + spec%dx(1)*i
            spec%x(2,l) = spec%x_ini(2) + spec%dx(2)*j
          end do
        end do
        
      case(3)
        l = 0
        do i = 1, ngrid(1)
          do j = 1, ngrid(2)
            do k = 1, ngrid(3)
              l = l + 1
              spec%x(1,l) = spec%x_ini(1) + spec%dx(1)*i
              spec%x(2,l) = spec%x_ini(2) + spec%dx(2)*j
              spec%x(3,l) = spec%x_ini(3) + spec%dx(3)*k
            end do
          end do
        end do

      end select


    end subroutine init_species

end module species
