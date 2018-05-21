subroutine calc_velocity_from_cond_wf(spec_in,zwfn_in,r_in,v_out)
  use global_variables
  implicit none
  type(species_t),intent(in) :: spec_in
  complex(8),intent(in)    :: zwfn_in(1:spec_in%ngrid_tot)
  real(8),intent(in)       :: r_in(1:spec_in%ndim)
  real(8),intent(out)      :: v_out(1:spec_in%ndim)

  select case(spec_in%ndim)
  case(1)
    call calc_velocity_from_cond_wf_1d(spec_in,zwfn_in,r_in(1),v_out(1))
  case default
    stop 'Error: calc_velocity_from_cond_wf'
  end select

end subroutine calc_velocity_from_cond_wf
!-----------------------------------------------------------------------------------------
subroutine calc_velocity_from_cond_wf_1d(spec_in,zwfn_in,x_in,v_out)
  use global_variables
  implicit none
  type(species_t),intent(in) :: spec_in
  complex(8),intent(in)    :: zwfn_in(1:spec_in%ngrid_tot)
  real(8),intent(in)       :: x_in
  real(8),intent(out)      :: v_out
  complex(8) :: zg(1)
  real(8) :: vt1, vt2
  real(8) :: t1, t2
  integer :: ix

  ix = aint( (x_in-spec_in%x_ini(1))/spec_in%dx(1) )
  t1 = x_in - (ix*spec_in%dx(1) + spec_in%x_ini(1))
  t1 = t1/spec_in%dx(1)

  if(ix > 0 .and. ix <= spec_in%ngrid(1))then
    call gradient_local(zwfn_in, zg, spec_in%ngrid, spec_in%dx, 1d0, ix)
    vt1 = real( conjg(zwfn_in(ix))*(-zI)*zg(1)/spec_in%mass ) &
      /abs(zwfn_in(ix))**2
  else
    vt1 = 0d0
  end if

  ix = aint( (x_in-spec_in%x_ini(1))/spec_in%dx(1) ) + 1
  t2 = (ix*spec_in%dx(1) + spec_in%x_ini(1)) - x_in
  t2 = t2/spec_in%dx(1)
  if(ix > 0 .and. ix <= spec_in%ngrid(1))then
    call gradient_local(zwfn_in, zg, spec_in%ngrid, spec_in%dx, 1d0, ix)
    vt2 = real( conjg(zwfn_in(ix))*(-zI)*zg(1)/spec_in%mass ) &
      /abs(zwfn_in(ix))**2
  else
    vt2 = 0d0
  end if

  v_out = t2*vt1 + t1*vt2

end subroutine calc_velocity_from_cond_wf_1d
