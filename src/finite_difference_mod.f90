module finite_difference_mod
  implicit none
  
  private

! Laplacian coefficients (7-point formula)
  real(8),parameter :: ct0 = -49d0/18d0, &
                       ct1 = 3d0/2d0,    &
                       ct2 = -3d0/20d0,  &
                       ct3 = 1d0/90d0


  interface laplacian
     module procedure laplacian_complex
     module procedure laplacian_real
  end interface laplacian

  public :: laplacian

contains
!-----------------------------------------------------------------------------------------
  subroutine laplacian_real(f_in, f_out, nx, dx, factor)
    implicit none
    real(8),intent(in) :: f_in(:), dx(:)
    real(8),intent(out) :: f_out(:)
    integer,intent(in) :: nx(:)
    real(8),intent(in),optional :: factor
    integer :: ndim
    real(8) :: factor0

    ndim = ubound(nx, 1)
    factor0 = 1d0
    if(present(factor))factor0 = factor

    select case(ndim)
    case(1)
      call laplacian_real_1d(f_in, f_out, nx(1), dx(1), factor0)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in Laplacian operator"
      stop
    end select

  end subroutine laplacian_real
!-----------------------------------------------------------------------------------------
  subroutine laplacian_complex(zf_in, zf_out, nx, dx, factor)
    implicit none
    complex(8),intent(in)   :: zf_in(:)
    complex(8),intent(out)  :: zf_out(:)
    real(8),intent(in)      :: dx(:)
    integer,intent(in) :: nx(:)
    real(8),intent(in),optional :: factor
    integer :: ndim
    real(8) :: factor0

    ndim = ubound(nx, 1)
    factor0 = 1d0
    if(present(factor))factor0 = factor

    select case(ndim)
    case(1)
      call laplacian_complex_1d(zf_in, zf_out, nx(1), dx(1), factor0)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in Laplacian operator"
      stop
    end select

  end subroutine laplacian_complex
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine laplacian_real_1d(f_in, f_out, nx, dx, factor0)
    implicit none
    integer,intent(in) :: nx
    real(8),intent(in) :: f_in(nx), dx, factor0
    real(8),intent(out) :: f_out(nx)
    real(8) :: dx_i2
    real(8) :: c0,c1,c2,c3,c4
    integer :: ix

    dx_i2 = factor0/dx**2
    c0 = ct0*dx_i2
    c1 = ct1*dx_i2
    c2 = ct2*dx_i2
    c3 = ct3*dx_i2

    f_out(1) = c0*f_in(1) &
              +c1*(f_in(1+1) ) &
              +c2*(f_in(1+2) ) &
              +c3*(f_in(1+3) )

    f_out(2) = c0*f_in(2) &
               +c1*(f_in(2+1) + f_in(2-1)) &
               +c2*(f_in(2+2) ) &
               +c3*(f_in(2+3) )

    f_out(3) = c0*f_in(3) &
              +c1*(f_in(3+1) + f_in(3-1)) &
              +c2*(f_in(3+2) + f_in(3-2)) &
              +c3*(f_in(3+3))

    do ix = 1+3, nx-3
      f_out(ix) = c0*f_in(ix) &
                 +c1*(f_in(ix+1) + f_in(ix-1)) &
                 +c2*(f_in(ix+2) + f_in(ix-2)) &
                 +c3*(f_in(ix+3) + f_in(ix-3))
    end do

    f_out(nx-2) = c0*f_in(nx-2) &
                 +c1*(f_in(nx-2+1) + f_in(nx-2-1)) &
                 +c2*(f_in(nx-2+2) + f_in(nx-2-2)) &
                 +c3*(               f_in(nx-2-3))

    f_out(nx-1) = c0*f_in(nx-1) &
                 +c1*(f_in(nx-1+1) + f_in(nx-1-1)) &
                 +c2*(               f_in(nx-1-2)) &
                 +c3*(               f_in(nx-1-3))

    f_out(nx) = c0*f_in(nx) &
               +c1*(f_in(nx-1)) &
               +c2*(f_in(nx-2)) &
               +c3*(f_in(nx-3))


  end subroutine laplacian_real_1d
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine laplacian_complex_1d(zf_in, zf_out, nx, dx, factor0)
    implicit none
    integer,intent(in) :: nx
    complex(8),intent(in)  :: zf_in(nx)
    complex(8),intent(out) :: zf_out(nx)
    real(8),intent(in) :: dx, factor0
    real(8) :: dx_i2
    real(8) :: c0,c1,c2,c3,c4
    integer :: ix

    dx_i2 = factor0/dx**2
    c0 = ct0*dx_i2
    c1 = ct1*dx_i2
    c2 = ct2*dx_i2
    c3 = ct3*dx_i2

    zf_out(1) = c0*zf_in(1) &
              +c1*(zf_in(1+1) ) &
              +c2*(zf_in(1+2) ) &
              +c3*(zf_in(1+3) )

    zf_out(2) = c0*zf_in(2) &
               +c1*(zf_in(2+1) + zf_in(2-1)) &
               +c2*(zf_in(2+2) ) &
               +c3*(zf_in(2+3) )

    zf_out(3) = c0*zf_in(3) &
              +c1*(zf_in(3+1) + zf_in(3-1)) &
              +c2*(zf_in(3+2) + zf_in(3-2)) &
              +c3*(zf_in(3+3))

    do ix = 1+3, nx-3
      zf_out(ix) = c0*zf_in(ix) &
                 +c1*(zf_in(ix+1) + zf_in(ix-1)) &
                 +c2*(zf_in(ix+2) + zf_in(ix-2)) &
                 +c3*(zf_in(ix+3) + zf_in(ix-3))
    end do

    zf_out(nx-2) = c0*zf_in(nx-2) &
                 +c1*(zf_in(nx-2+1) + zf_in(nx-2-1)) &
                 +c2*(zf_in(nx-2+2) + zf_in(nx-2-2)) &
                 +c3*(               zf_in(nx-2-3))

    zf_out(nx-1) = c0*zf_in(nx-1) &
                 +c1*(zf_in(nx-1+1) + zf_in(nx-1-1)) &
                 +c2*(               zf_in(nx-1-2)) &
                 +c3*(               zf_in(nx-1-3))

    zf_out(nx) = c0*zf_in(nx) &
               +c1*(zf_in(nx-1)) &
               +c2*(zf_in(nx-2)) &
               +c3*(zf_in(nx-3))


  end subroutine laplacian_complex_1d

end module finite_difference_mod
