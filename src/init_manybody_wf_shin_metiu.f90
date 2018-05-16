subroutine init_manybody_wf_shin_metiu
  use global_variables
  implicit none
  integer :: ix_e, ix_i
  real(8) :: xe, xi
  real(8) :: norm

!! FIXME:
!! Initialization of many-body wavefunction for Shin-Metiue will be implemented later.
!! The gaussian is employed temporary.

  do ix_i = 1, spec(2)%ngrid(1)
    xi = spec(2)%x_ini(1) + spec(2)%dx(1)*ix_i
    do ix_e = 1, spec(1)%ngrid(1)
      xe = spec(1)%x_ini(1) + spec(1)%dx(1)*ix_e
      
      zwfn_ini_2p(ix_e,ix_i) = exp(-0.5d0*abs(xe)**2 -0.5d0*abs(xi)**2)

    end do
  end do

  norm = sum(abs(zwfn_ini_2p)**2)*spec(1)%dx(1)*spec(2)%dx(1)

  zwfn_ini_2p = zwfn_ini_2p/sqrt(norm)

end subroutine init_manybody_wf_shin_metiu
