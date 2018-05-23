module finite_difference_mod
  implicit none
  
  private

! Laplacian coefficients (7-point formula)
  real(8),parameter :: ct0 = -49d0/18d0, &
                       ct1 = 3d0/2d0,    &
                       ct2 = -3d0/20d0,  &
                       ct3 = 1d0/90d0

! Gradient coefficients (7-point formula)
  real(8),parameter :: gt0 = 0d0, &
                       gt1 = 3d0/4d0,    &
                       gt2 = -3d0/20d0,  &
                       gt3 = 1d0/60d0

  interface laplacian
     module procedure laplacian_complex
     module procedure laplacian_real
  end interface laplacian

  interface gradient
     module procedure gradient_complex
     module procedure gradient_real
  end interface gradient
  
  interface gradient_local
     module procedure gradient_local_complex
     module procedure gradient_local_real
  end interface gradient_local

  public :: laplacian, &
            gradient, &
            gradient_local

contains
!-----------------------------------------------------------------------------------------
  subroutine laplacian_real(f_in, f_out, nx, dx, factor)
    implicit none
    real(8),intent(in)  :: f_in(:)
    real(8),intent(out) :: f_out(:)
    real(8),intent(in)  :: dx(:)
    integer,intent(in)  :: nx(:)
    real(8),intent(in)  :: factor
    integer :: ndim

    ndim = ubound(nx, 1)

    select case(ndim)
    case(1)
      call laplacian_real_1d(f_in, f_out, nx(1), dx(1), factor)
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
    integer,intent(in)      :: nx(:)
    real(8),intent(in)      :: factor
    integer :: ndim
    real(8) :: factor0

    ndim = ubound(nx, 1)

    select case(ndim)
    case(1)
      call laplacian_complex_1d(zf_in, zf_out, nx(1), dx(1), factor)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in Laplacian operator"
      stop
    end select

  end subroutine laplacian_complex
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine laplacian_real_1d(f_in, f_out, nx, dx, factor)
    implicit none
    integer,intent(in)  :: nx
    real(8),intent(in)  :: f_in(nx)
    real(8),intent(out) :: f_out(nx)
    real(8),intent(in)  :: dx, factor
    real(8) :: dx_i2
    real(8) :: c0,c1,c2,c3,c4
    integer :: ix

    dx_i2 = factor/dx**2
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
  subroutine laplacian_complex_1d(zf_in, zf_out, nx, dx, factor)
    implicit none
    integer,intent(in) :: nx
    complex(8),intent(in)  :: zf_in(nx)
    complex(8),intent(out) :: zf_out(nx)
    real(8),intent(in) :: dx, factor
    real(8) :: dx_i2
    real(8) :: c0,c1,c2,c3,c4
    integer :: ix

    dx_i2 = factor/dx**2
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
!-----------------------------------------------------------------------------------------
  subroutine gradient_real(f_in, f_out, nx, dx, factor)
    implicit none
    real(8),intent(in)  :: f_in(:)
    real(8),intent(out) :: f_out(:,:)
    real(8),intent(in)  :: dx(:)
    integer,intent(in)  :: nx(:)
    real(8),intent(in)  :: factor
    integer :: ndim

    ndim = ubound(nx, 1)

    select case(ndim)
    case(1)
      call gradient_real_1d(f_in, f_out, nx(1), dx(1), factor)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in gradient operator"
      stop
    end select

  end subroutine gradient_real
!-----------------------------------------------------------------------------------------
  subroutine gradient_complex(zf_in, zf_out, nx, dx, factor)
    implicit none
    complex(8),intent(in)  :: zf_in(:)
    complex(8),intent(out) :: zf_out(:,:)
    real(8),intent(in)  :: dx(:)
    integer,intent(in)  :: nx(:)
    real(8),intent(in)  :: factor
    integer :: ndim

    ndim = ubound(nx, 1)

    select case(ndim)
    case(1)
      call gradient_complex_1d(zf_in, zf_out, nx(1), dx(1), factor)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in gradient operator"
      stop
    end select

  end subroutine gradient_complex
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine gradient_real_1d(f_in, f_out, nx, dx, factor)
    implicit none
    integer,intent(in)  :: nx
    real(8),intent(in)  :: f_in(nx)
    real(8),intent(out) :: f_out(nx,1)
    real(8),intent(in)  :: dx, factor
    real(8) :: dx_i
    real(8) :: g0,g1,g2,g3,g4
    integer :: ix

    dx_i = factor/dx
    g0 = gt0*dx_i
    g1 = gt1*dx_i
    g2 = gt2*dx_i
    g3 = gt3*dx_i

    f_out(1,1) = g1*(f_in(1+1) ) &
                +g2*(f_in(1+2) ) &
                +g3*(f_in(1+3) )

    f_out(2,1) = g1*(f_in(2+1) - f_in(2-1)) &
                +g2*(f_in(2+2) ) &
                +g3*(f_in(2+3) )

    f_out(3,1) = g1*(f_in(3+1) - f_in(3-1)) &
                +g2*(f_in(3+2) - f_in(3-2)) &
                +g3*(f_in(3+3))

    do ix = 1+3, nx-3
      f_out(ix,1) = g1*(f_in(ix+1) - f_in(ix-1)) &
                   +g2*(f_in(ix+2) - f_in(ix-2)) &
                   +g3*(f_in(ix+3) - f_in(ix-3))
    end do

    f_out(nx-2,1) = g1*(f_in(nx-2+1) - f_in(nx-2-1)) &
                   +g2*(f_in(nx-2+2) - f_in(nx-2-2)) &
                   +g3*(             - f_in(nx-2-3))

    f_out(nx-1,1) = g1*(f_in(nx-1+1) - f_in(nx-1-1)) &
                   +g2*(             - f_in(nx-1-2)) &
                   +g3*(             - f_in(nx-1-3))

    f_out(nx,1) = g1*( -f_in(nx-1)) &
                 +g2*( -f_in(nx-2)) &
                 +g3*( -f_in(nx-3))


  end subroutine gradient_real_1d
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine gradient_complex_1d(zf_in, zf_out, nx, dx, factor)
    implicit none
    integer,intent(in)  :: nx
    complex(8),intent(in)  :: zf_in(nx)
    complex(8),intent(out) :: zf_out(nx,1)
    real(8),intent(in)  :: dx, factor
    real(8) :: dx_i
    real(8) :: g0,g1,g2,g3,g4
    integer :: ix

    dx_i = factor/dx
    g0 = gt0*dx_i
    g1 = gt1*dx_i
    g2 = gt2*dx_i
    g3 = gt3*dx_i

    zf_out(1,1) = g1*(zf_in(1+1) ) &
                 +g2*(zf_in(1+2) ) &
                 +g3*(zf_in(1+3) )

    zf_out(2,1) = g1*(zf_in(2+1) - zf_in(2-1)) &
                 +g2*(zf_in(2+2) ) &
                 +g3*(zf_in(2+3) )

    zf_out(3,1) = g1*(zf_in(3+1) - zf_in(3-1)) &
                 +g2*(zf_in(3+2) - zf_in(3-2)) &
                 +g3*(zf_in(3+3))

    do ix = 1+3, nx-3
      zf_out(ix,1) = g1*(zf_in(ix+1) - zf_in(ix-1)) &
                    +g2*(zf_in(ix+2) - zf_in(ix-2)) &
                    +g3*(zf_in(ix+3) - zf_in(ix-3))
    end do

    zf_out(nx-2,1) = g1*(zf_in(nx-2+1) - zf_in(nx-2-1)) &
                    +g2*(zf_in(nx-2+2) - zf_in(nx-2-2)) &
                    +g3*(             - zf_in(nx-2-3))

    zf_out(nx-1,1) = g1*(zf_in(nx-1+1) - zf_in(nx-1-1)) &
                    +g2*(             - zf_in(nx-1-2)) &
                    +g3*(             - zf_in(nx-1-3))

    zf_out(nx,1) = g1*( -zf_in(nx-1)) &
                  +g2*( -zf_in(nx-2)) &
                  +g3*( -zf_in(nx-3))


  end subroutine gradient_complex_1d
!-----------------------------------------------------------------------------------------
  subroutine gradient_local_real(f_in, f_out, nx, dx, factor, ix)
    implicit none
    real(8),intent(in)  :: f_in(:)
    real(8),intent(out) :: f_out(:)
    real(8),intent(in)  :: dx(:)
    integer,intent(in)  :: nx(:)
    real(8),intent(in)  :: factor
    integer,intent(in)  :: ix
    integer :: ndim

    ndim = ubound(nx, 1)

    select case(ndim)
    case(1)
      call gradient_local_real_1d(f_in, f_out, nx(1), dx(1), factor, ix)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in gradient operator"
      stop
    end select

  end subroutine gradient_local_real
!-----------------------------------------------------------------------------------------
  subroutine gradient_local_complex(zf_in, zf_out, nx, dx, factor, ix)
    implicit none
    complex(8),intent(in)  :: zf_in(:)
    complex(8),intent(out) :: zf_out(:)
    real(8),intent(in)     :: dx(:)
    integer,intent(in)     :: nx(:)
    real(8),intent(in)     :: factor
    integer,intent(in)     :: ix
    integer :: ndim

    ndim = ubound(nx, 1)

    select case(ndim)
    case(1)
      call gradient_local_complex_1d(zf_in, zf_out, nx(1), dx(1), factor, ix)
    case default
      write(*,"(A)")"Fatal Error: Invalid dimension in gradient operator"
      stop
    end select

  end subroutine gradient_local_complex
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine gradient_local_real_1d(f_in, f_out, nx, dx, factor, ix)
    implicit none
    integer,intent(in)  :: nx
    real(8),intent(in)  :: f_in(nx)
    real(8),intent(out) :: f_out(1)
    real(8),intent(in)  :: dx, factor
    integer,intent(in) :: ix
    real(8) :: dx_i
    real(8) :: g0,g1,g2,g3,g4


    dx_i = factor/dx
    g0 = gt0*dx_i
    g1 = gt1*dx_i
    g2 = gt2*dx_i
    g3 = gt3*dx_i

    if(ix >= 1+3 .and. ix < nx-3)then
      f_out(1) = g1*(f_in(ix+1) - f_in(ix-1)) &
                +g2*(f_in(ix+2) - f_in(ix-2)) &
                +g3*(f_in(ix+3) - f_in(ix-3))
    else if(ix ==1)then
      f_out(1) = g1*(f_in(1+1) ) &
                +g2*(f_in(1+2) ) &
                +g3*(f_in(1+3) )
    else if(ix ==2)then
      f_out(1) = g1*(f_in(2+1) - f_in(2-1)) &
                +g2*(f_in(2+2) ) &
                +g3*(f_in(2+3) )
    else if(ix ==3)then
      f_out(1) = g1*(f_in(3+1) - f_in(3-1)) &
                +g2*(f_in(3+2) - f_in(3-2)) &
                +g3*(f_in(3+3))
    else if(ix == nx-2)then
      f_out(1) = g1*(f_in(nx-2+1) - f_in(nx-2-1)) &
                +g2*(f_in(nx-2+2) - f_in(nx-2-2)) &
                +g3*(             - f_in(nx-2-3))
    else if(ix == nx-1)then
      f_out(1) = g1*(f_in(nx-1+1) - f_in(nx-1-1)) &
                +g2*(             - f_in(nx-1-2)) &
                +g3*(             - f_in(nx-1-3))
    else if(ix == nx)then
      f_out(1) = g1*( -f_in(nx-1)) &
                +g2*( -f_in(nx-2)) &
                +g3*( -f_in(nx-3))
    else
      stop 'Error in gradient_local_real_1d'
    end if

  end subroutine gradient_local_real_1d
!-----------------------------------------------------------------------------------------
! NOTE: nx has to be equal to or larger than 6
  subroutine gradient_local_complex_1d(zf_in, zf_out, nx, dx, factor, ix)
    implicit none
    integer,intent(in)  :: nx
    complex(8),intent(in)  :: zf_in(nx)
    complex(8),intent(out) :: zf_out(1)
    real(8),intent(in)  :: dx, factor
    integer,intent(in) :: ix
    real(8) :: dx_i
    real(8) :: g0,g1,g2,g3,g4


    dx_i = factor/dx
    g0 = gt0*dx_i
    g1 = gt1*dx_i
    g2 = gt2*dx_i
    g3 = gt3*dx_i

    if(ix >= 1+3 .and. ix < nx-3)then
      zf_out(1) = g1*(zf_in(ix+1) - zf_in(ix-1)) &
                 +g2*(zf_in(ix+2) - zf_in(ix-2)) &
                 +g3*(zf_in(ix+3) - zf_in(ix-3))
    else if(ix ==1)then
      zf_out(1) = g1*(zf_in(1+1) ) &
                 +g2*(zf_in(1+2) ) &
                 +g3*(zf_in(1+3) )
    else if(ix ==2)then
      zf_out(1) = g1*(zf_in(2+1) - zf_in(2-1)) &
                 +g2*(zf_in(2+2) ) &
                 +g3*(zf_in(2+3) )
    else if(ix ==3)then
      zf_out(1) = g1*(zf_in(3+1) - zf_in(3-1)) &
                 +g2*(zf_in(3+2) - zf_in(3-2)) &
                 +g3*(zf_in(3+3))
    else if(ix == nx-2)then
      zf_out(1) = g1*(zf_in(nx-2+1) - zf_in(nx-2-1)) &
                 +g2*(zf_in(nx-2+2) - zf_in(nx-2-2)) &
                 +g3*(             - zf_in(nx-2-3))
    else if(ix == nx-1)then
      zf_out(1) = g1*(zf_in(nx-1+1) - zf_in(nx-1-1)) &
                 +g2*(             - zf_in(nx-1-2)) &
                 +g3*(             - zf_in(nx-1-3))
    else if(ix == nx)then
      zf_out(1) = g1*( -zf_in(nx-1)) &
                 +g2*( -zf_in(nx-2)) &
                 +g3*( -zf_in(nx-3))
    else
      stop 'Error in gradient_local_real_1d'
    end if

  end subroutine gradient_local_complex_1d

!-----------------------------------------------------------------------------------------

end module finite_difference_mod
