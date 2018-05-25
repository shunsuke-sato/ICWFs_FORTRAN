subroutine init_manybody_wf_ee_scattering
  use global_variables
  implicit none
  integer,parameter :: ncg = 100
  real(8),parameter :: res_epsilon = 1d-8
  real(8),allocatable :: veff(:)
  real(8),allocatable :: xvec(:),gvec(:),pvec(:)
  real(8),allocatable :: hxvec(:),hpvec(:)
  real(8) :: gg0, gg, xx, xhx, pp, php, xp, xhp
  real(8) :: alpha, beta, aa, bb, cc, ss, lambda, res
  integer :: icg
  real(8) :: dx
  integer :: ix,iy,nx
  real(8) :: tmp
  complex(8) :: ztmp

  nx = spec(1)%ngrid(1)
  dx = spec(1)%dx(1)

  allocate(xvec(nx),pvec(nx),gvec(nx))
  allocate(hxvec(nx),hpvec(nx))
  allocate(veff(nx))

! set one-body potential
  veff(:) = spec(1)%v0(:)

! initialize wavefunction for CG
  do ix = 1,nx
    call random_double(xvec(ix))
  end do


  xx = sum(xvec**2)*dx
  call hphi(xvec, hxvec)
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
    xx  = sum(xvec**2)*dx
    xhx = sum(xvec*hxvec)*dx
    lambda = xhx/xx

    res = sum(abs(hxvec-lambda*xvec)**2)*dx/xx
    if(res<res_epsilon)exit
    
    
  end do
  call hphi(xvec, hxvec)
  xx  = sum(xvec**2)*dx
  xhx = sum(xvec*hxvec)*dx
  lambda = xhx/xx
  res = sum(abs(hxvec-lambda*xvec)**2)*dx/xx
  xvec = xvec/xx
  if(if_root_global)then
    write(*,*)"icg =",icg
    write(*,"(A,2x,2e26.16e3)")"lambda,res=",lambda,res
  end if

  do ix = 1, spec(1)%ngrid(1)
    xx = spec(1)%x(1,ix)
    do iy = 1, spec(1)%ngrid(1)
      yy = spec(1)%x(2,iy)
      ztmp = exp(-0.5d0*((xx-x0)/sigma)**2)*exp(zI*k_mom*xx)*xvec(iy)
      ztmp = ztmp + exp(-0.5d0*((yy-x0)/sigma)**2)*exp(zI*k_mom*yy)*xvec(ix)
      zwfn_ini_2p(ix,iy) = ztmp
    end do
  end do

  tmp = sum(abs(zwfn_ini_2p)**2)*spec(1)%dx(1)**2
  zwfn_ini_2p = zwfn_ini_2p /sqrt(tmp)

  contains

  subroutine hphi(f_in,f_out)
    implicit none
    real(8) :: f_in(nx), f_out(nx)

    call laplacian(f_in,f_out,spec(1)%ngrid,spec(1)%dx,-0.5d0)
    f_out = f_out + f_in*veff

  end subroutine hphi


end subroutine init_manybody_wf_ee_scattering
