module math_mod
  implicit none
  private

  real(8),parameter,public :: pi = 3.141592653589793d0
  complex(8),parameter,public :: zI=(0.d0,1.d0)

  public :: erf_x
contains
!    
  function erf_x(x) result(y)
    real(8),intent(in) :: x
    real(8),parameter :: epsilon_s = 1d-3
    real(8) :: y
    
    if(abs(x) > epsilon_s)then
      y = erf(x)/x
    else
      y = 2d0/sqrt(pi)*( 1d0 - x**2/3d0 + x**4/10d0 - x**6/42d0 + x**8/216d0)
    end if
    
  end function erf_x
  
end module math_mod
