subroutine operators
  use global_variables
  use interaction_mod
  implicit none
  integer :: ispec,i

  call write_message('Start: operators')
  call comm_barrier

! Initialize one-body potential

  do ispec = 1, num_species
    select case(ispec)
    case(1)
      do i = 1, spec(ispec)%ngrid_tot
        spec(ispec)%v0(i) = one_body_pot_1(spec(ispec)%x(:,i))
      end do
    case(2)
      do i = 1, spec(ispec)%ngrid_tot
        spec(ispec)%v0(i) = one_body_pot_2(spec(ispec)%x(:,i))
      end do
    end select
  end do


  call comm_barrier
  call write_message('Finish: operators')

end subroutine operators


