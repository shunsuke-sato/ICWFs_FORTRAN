subroutine init_manybody_wf_shin_metiu
  use global_variables
  implicit none
  integer,parameter :: Rion = -4d0
  real(8),parameter :: sigma = 1d0/sqrt(2.85d0)
  integer,parameter :: norb = 2
  real(8),parameter :: res_epsilon = 1d-8
  real(8),allocatable :: phi(:,:)
  real(8),allocatable :: veff(:)
  real(8),allocatable :: xvec(:),gvec(:),pvec(:)
  real(8),allocatable :: hxvec(:),hpvec(:)
  integer,parameter :: niter_max = 300, ncg = 30
  real(8) :: residual(norb), eps(norb)
  integer :: ix_e, ix_i
  integer :: ix, iorb, jorb, nx
  integer :: iter, icg
  real(8) :: xe(1), xi(1), dx
  real(8) :: norm
  real(8) :: r1
  real(8) :: gg0, gg, xx, xhx, pp, php, xp, xhp
  real(8) :: alpha, beta, aa, bb, cc, ss, lambda
  
  nx = spec(1)%ngrid(1)
  dx = spec(1)%dx(1)

  allocate(phi(nx, norb))
  allocate(xvec(nx),pvec(nx),gvec(nx))
  allocate(hxvec(nx),hpvec(nx))
  allocate(veff(nx))

! initialize wavefunction
  do iorb = 1,norb
    do ix = 1, nx
      
      call random_number(r1)
      phi(ix,iorb) = r1

    end do
  end do
  call gram_schmidt

  if(if_root_global)open(20,file="BO_Shin_Metiu.out")
! Compute BO states
  do ix_i = 1, spec(2)%ngrid(1)
    xi(1) = spec(2)%x_ini(1) + spec(2)%dx(1)*ix_i

    ! set BO potential
    do ix = 1,nx
      xe(1) = spec(1)%x_ini(1) + spec(1)%dx(1)*ix
      veff(ix) = two_body_pot_1_2(xe,xi)
    end do
    veff(:) = veff(:) + spec(1)%v0(:)


!! initialize wavefunction
!  do iorb = 1,norb
!    do ix = 1, nx
!      
!      call random_number(r1)
!      phi(ix,iorb) = r1
!
!    end do
!  end do
!  call gram_schmidt

    do iter = 1, Niter_max
      do iorb = 1, norb
        xvec(:) = phi(:,iorb)

        xx = sum(xvec**2)*dx
        call hphi(xvec, hxvec)
        call gram_schmidt_partial(hxvec,iorb)
        xhx = sum(xvec*hxvec)*dx
        lambda = xhx/xx
        
        do icg = 1, ncg
          
          gvec = 2d0*(hxvec - lambda*xvec)/xx
          gg0  = sum(gvec**2)*dx
          select case(icg)
          case(1)
            pvec = -gvec
          case default
            beta = gg0/gg
            pvec = -gvec + beta*pvec
          end select
          gg = gg0
          
          call hphi(pvec, hpvec)
          call gram_schmidt_partial(hpvec,iorb)
          pp  = sum(pvec**2)*dx
          php = sum(pvec*hpvec)*dx
          xp  = sum(pvec*xvec)*dx
          xhp = sum(hxvec*pvec)*dx
          
          aa=php*xp-xhp*pp
          bb=php*xx-xhx*pp
          cc=xhp*xx-xhx*xp
          ss=bb**2-4d0*aa*cc
          if(ss > 0d0)then
            alpha=(-bb+sqrt(ss))/(2d0*aa)
          else
            exit
          end if
          
          xvec = xvec + alpha*pvec
          call hphi(xvec, hxvec)
          call gram_schmidt_partial(hxvec,iorb)
          xx  = sum(xvec**2)*dx
          xhx = sum(xvec*hxvec)*dx
          lambda = xhx/xx
          
          
        end do
        xx = sum(xvec(:)**2)*dx
        phi(:,iorb) = xvec(:)/sqrt(xx)
        
      end do
      call gram_schmidt
      
      !        if(mod(it,100) == 0)then

      
      do iorb = 1, norb
        xvec(:) = phi(:,iorb)
        call hphi(xvec, hxvec)
        eps(iorb) = sum(xvec*hxvec)*dx
        residual(iorb) = sum((hxvec - eps(iorb)*xvec)**2)*dx
      end do
      if(sum(residual) < res_epsilon)exit
      
    end do

    write(message(1), "(I9, 2x, I9)")ix_i, iter
    message(1) = "ixi, Iter ="//trim(message(1))
    call write_message(message(1))
    write(message(1), "(e16.6e3)")sum(residual(:))
    message(1) = "residual               = "//trim(message(1))
    message(2) = "  "
    call write_message(message(1:2))

    if(if_root_global)write(20,"(999e26.16e3)")xi(1),eps(:)+one_body_pot_2(xi)
    zwfn_ini_2p(:,ix_i) = phi(:,norb)

  end do

      
  if(if_root_global)close(20)

  do ix_i = 1, spec(2)%ngrid(1)
    xi = spec(2)%x_ini(1) + spec(2)%dx(1)*ix_i
      
    zwfn_ini_2p(:,ix_i) = zwfn_ini_2p(:,ix_i)*exp(-0.5d0*((xi(1)-Rion)/sigma)**2)

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

  subroutine gram_schmidt_partial(f_inout,i)
    implicit none
    real(8),intent(inout) :: f_inout(nx)
    integer,intent(in) :: i
    integer :: j
    real(8) :: ss

    do j = 1, i-1
      ss = sum(phi(:,j)*f_inout(:))*dx
      f_inout(:) = f_inout(:) - ss* phi(:,j)
    end do
!    

  end subroutine gram_schmidt_partial

  subroutine hphi(f_in,f_out)
    implicit none
    real(8) :: f_in(nx), f_out(nx)

    call laplacian(f_in,f_out,spec(1)%ngrid,spec(1)%dx,-0.5d0)
    f_out = f_out + f_in*veff

  end subroutine hphi
  

end subroutine init_manybody_wf_shin_metiu
