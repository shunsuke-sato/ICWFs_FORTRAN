module species
  implicit none

  type species_t
    integer :: ndim                      ! Dimensionality
    integer :: nparticle              ! Number of particles
    real(8) :: mass                      ! Mass of particle
    integer,allocatable :: ngrid(:)   ! Number of grid points
    integer :: ngrid_tot   ! Number of grid points
    integer,allocatable :: x_ini(:)    ! Initial position of particle coordinates
    integer,allocatable :: x_fin(:)    ! Final position of particle coordinates
    integer,allocatable :: dx(:)    ! Spacing of grid points
    complex(8),allocatable :: zwfn(:,:)  ! wavefunction
    real(8),allocatable :: r_particle(:,:) ! Position of particle

  end type species_t

  contains
    subroutine init_species(spec,ndim, nparticle, mass, ngrid, x_ini, x_fin)
      type(species_t),intent(inout)  :: spec
      integer,intent(in) :: ndim, nparticle, ngrid(1:ndim)
      real(8),intent(in) :: mass, x_ini(1:ndim), x_fin(1:ndim)
      integer :: ngrid_tot

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

      allocate(spec%zwfn(ngrid_tot,nparticle))
      allocate(spec%r_particle(ndim,nparticle))

      spec%dx(1:ndim) = (x_fin(1:ndim) - x_ini(1:ndim))/(ngrid(1:ndim)+1)


    end subroutine init_species

end module species
