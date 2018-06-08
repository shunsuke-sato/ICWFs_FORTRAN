module math_mod
  implicit none
  private

  real(8),parameter,public :: pi = 3.141592653589793d0
  complex(8),parameter,public :: zI=(0.d0,1.d0)

  public :: erf_x,&
            pseudo_inverse


  interface pseudo_inverse
     module procedure pseudo_inverse_real8
     module procedure pseudo_inverse_complex8
  end interface pseudo_inverse

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

  subroutine pseudo_inverse_real8(a_in, a_pinv)
    implicit none
    real(8),intent(in)  :: a_in(:,:)
    real(8),intent(out) :: a_pinv(:,:)
    real(8),allocatable :: a(:,:), u(:,:),vt(:,:)
    real(8),allocatable :: s(:)
    real(8),allocatable :: ut(:,:),v(:,:), at(:,:)
    integer :: m, n
    real(8),allocatable :: work(:)
    integer :: lwork, info
    real(8) :: tolerance
    integer :: i,j

    m = size(a_in,1)
    n = size(a_in,2)
    lwork = 10*max(1,3*min(m,n) + max(m,n), 5*min(m,n))

    allocate(a(m,n), s(min(m,n)),u(m,m), vt(n,n))
    allocate(at(n,m), ut(m,m), v(n,n))
    allocate(work(lwork))


    a = a_in

    call dgesvd ('A', 'A', m, n, a, m, s, u, m, vt, n, work, lwork, info)

    tolerance = max(1d-16*max(m,n)*maxval(s), 0d0 )
    do i = 1, min(m,n)
      if(s(i)>tolerance)then
        s(i) = 1d0/s(i)
      else
        s(i) = 0d0
      end if
    end do


    at = 0d0
    ut = transpose(u)
    v  = transpose(vt)

    do j = 1, m
      do i = 1, min(m,n)
        at(i,j) = s(i)*ut(i,j)
      end do
    end do

    a_pinv = matmul(v,at)

  end subroutine pseudo_inverse_real8
!-------------------------------------------------------------------------------
  subroutine pseudo_inverse_complex8(a_in, a_pinv)
    implicit none
    complex(8),intent(in)  :: a_in(:,:)
    complex(8),intent(out) :: a_pinv(:,:)
    complex(8),allocatable :: a(:,:), u(:,:),vt(:,:)
    real(8),allocatable :: s(:)
    complex(8),allocatable :: ut(:,:),v(:,:), at(:,:)
    integer :: m, n
    complex(8),allocatable :: work(:)
    real(8),allocatable :: rwork(:)
    integer :: lwork, info
    real(8) :: tolerance
    integer :: i,j

    m = size(a_in,1)
    n = size(a_in,2)
    lwork = 10*max(1,2*min(m,n)+max(m,n))

    allocate(a(m,n), s(min(m,n)),u(m,m), vt(n,n))
    allocate(at(n,m), ut(m,m), v(n,n))
    allocate(work(lwork), rwork(5*min(m,n)))

    a = a_in

    call zgesvd ('A', 'A', m, n, a, m, s, u, m, vt, n, work, lwork, rwork, info)

    tolerance = max(1d-16*max(m,n)*maxval(s), 0d0 )
    do i = 1, min(m,n)
      if(s(i)>tolerance)then
        s(i) = 1d0/s(i)
      else
        s(i) = 0d0
      end if
    end do

    at = 0d0
    ut = transpose(conjg(u))
    v  = transpose(conjg(vt))

    do j = 1, m
      do i = 1, min(m,n)
        at(i,j) = s(i)*ut(i,j)
      end do
    end do

    a_pinv = matmul(v,at)

  end subroutine pseudo_inverse_complex8
  
end module math_mod
