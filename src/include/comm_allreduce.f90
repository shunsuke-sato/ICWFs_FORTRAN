!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer(a_in, a_out, communicator, method)
  implicit none
  integer,intent(in) :: a_in
  integer,intent(in) :: a_out
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer :: id_comm, method_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  call MPI_Allreduce(a_in, a_out, 1, MPI_INTEGER, method_t, id_comm, ierr)

end subroutine comm_allreduce_integer
!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer_1d(a_in, a_out, communicator, method, nsize)
  implicit none
  integer,intent(in) :: a_in(:)
  integer,intent(in) :: a_out(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)

end subroutine comm_allreduce_integer_1d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer_2d(a_in, a_out, communicator, method, nsize)
  implicit none
  integer,intent(in) :: a_in(:,:)
  integer,intent(in) :: a_out(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)

end subroutine comm_allreduce_integer_2d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer_3d(a_in, a_out, communicator, method, nsize)
  implicit none
  integer,intent(in) :: a_in(:,:,:)
  integer,intent(in) :: a_out(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)

end subroutine comm_allreduce_integer_3d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8(a_in, a_out, communicator, method)
  implicit none
  real(8),intent(in) :: a_in
  real(8),intent(in) :: a_out
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer :: id_comm, method_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  call MPI_Allreduce(a_in, a_out, 1, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)

end subroutine comm_allreduce_real8
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8_1d(a_in, a_out, communicator, method, nsize)
  implicit none
  real(8),intent(in) :: a_in(:)
  real(8),intent(in) :: a_out(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)

end subroutine comm_allreduce_real8_1d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8_2d(a_in, a_out, communicator, method, nsize)
  implicit none
  real(8),intent(in) :: a_in(:,:)
  real(8),intent(in) :: a_out(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)

end subroutine comm_allreduce_real8_2d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8_3d(a_in, a_out, communicator, method, nsize)
  implicit none
  real(8),intent(in) :: a_in(:,:,:)
  real(8),intent(in) :: a_out(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)

end subroutine comm_allreduce_real8_3d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8(a_in, a_out, communicator, method)
  implicit none
  complex(8),intent(in) :: a_in
  complex(8),intent(in) :: a_out
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer :: id_comm, method_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  call MPI_Allreduce(a_in, a_out, 1, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)

end subroutine comm_allreduce_complex8
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8_1d(a_in, a_out, communicator, method, nsize)
  implicit none
  complex(8),intent(in) :: a_in(:)
  complex(8),intent(in) :: a_out(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)

end subroutine comm_allreduce_complex8_1d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8_2d(a_in, a_out, communicator, method, nsize)
  implicit none
  complex(8),intent(in) :: a_in(:,:)
  complex(8),intent(in) :: a_out(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)

end subroutine comm_allreduce_complex8_2d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8_3d(a_in, a_out, communicator, method, nsize)
  implicit none
  complex(8),intent(in) :: a_in(:,:,:)
  complex(8),intent(in) :: a_out(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)

end subroutine comm_allreduce_complex8_3d
!-------------------------------------------------------------------------------
