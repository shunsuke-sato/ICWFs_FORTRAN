! MAIN FILE
program main
  use global_variables
  implicit none
  integer :: m,n, i, j
  complex(8), allocatable :: a(:,:), a_pinv(:,:), aap(:,:), apa(:,:)
  complex(8), allocatable :: aapa(:,:),apaap(:,:)
  real(8) :: r1,r2
  m = 5
  n = 4

  allocate(a(m,n), a_pinv(n,m))
  allocate(aap(m,m), apa(n,n))
  allocate(aapa(m,n),apaap(m,n))

  do i = 1, m
    do j = 1, n
      call random_number(r1)
      call random_number(r2)
      a(i,j) = r1 + zI*r2
    end do
  end do
  call pseudo_inverse(a,a_pinv)

  aap = matmul(a, a_pinv)
  apa = matmul(a_pinv, a)
  aapa = matmul(a, apa)
  apaap = matmul(apa, a_pinv)

  write(*,*)"aap"
  do j = 1, n
    write(*,"(999e16.6e3)")aapa(:,j)-a(:,j)
  end do
  write(*,*)"apa"
  do j = 1, m
    write(*,"(999e16.6e3)")apaap(:,j)-a_pinv(:,j)
  end do

  stop

  call init_parallel
  call initialize_random_number_generator

  call parameters
  call operators

  call initial_state
!  call initial_distribution.f90

  call propagation


  call fin_parallel

end program main
