!-------------------------------------------------------------------------------
subroutine comm_sendrecv_integer(sbuf, dest, stag, rbuf, source, rtag, communicator)
  implicit none
  integer,intent(in) :: sbuf
  integer,intent(in) :: dest
  integer,intent(in) :: stag
  integer,intent(out) :: rbuf
  integer,intent(out) :: source
  integer,intent(in) :: rtag
  integer,intent(in),optional :: communicator
  integer :: id_comm
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr
  integer :: stype

  id_comm = int_switch(present(communicator), communicator, comm_group_global)


  call mpi_sendrecv(sbuf, 1, MPI_INTEGER, dest, stag, &
                    rbuf, 1, MPI_INTEGER, source, rtag, id_comm, status, ierr)

end subroutine comm_sendrecv_integer
!-------------------------------------------------------------------------------
subroutine comm_sendrecv_integer_1d(sbuf, dest, stag, rbuf, source, rtag, communicator)
  implicit none
  integer,intent(in) :: sbuf(:)
  integer,intent(in) :: dest
  integer,intent(in) :: stag
  integer,intent(out) :: rbuf(:)
  integer,intent(out) :: source
  integer,intent(in) :: rtag
  integer,intent(in),optional :: communicator
  integer :: id_comm
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)


  call mpi_sendrecv(sbuf, size(sbuf), MPI_INTEGER, dest, stag, &
                    rbuf, size(rbuf), MPI_INTEGER, source, rtag, id_comm, status, ierr)

end subroutine comm_sendrecv_integer_1d
!-------------------------------------------------------------------------------
subroutine comm_sendrecv_complex8(sbuf, dest, stag, rbuf, source, rtag, communicator)
  implicit none
  complex(8),intent(in) :: sbuf
  integer,intent(in)    :: dest
  integer,intent(in)    :: stag
  complex(8),intent(out):: rbuf
  integer,intent(out)   :: source
  integer,intent(in)    :: rtag
  integer,intent(in),optional :: communicator
  integer :: id_comm
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)


  call mpi_sendrecv(sbuf, 1, MPI_DOUBLE_COMPLEX, dest, stag, &
                    rbuf, 1, MPI_DOUBLE_COMPLEX, source, rtag, id_comm, status, ierr)

end subroutine comm_sendrecv_complex8
!-------------------------------------------------------------------------------
subroutine comm_sendrecv_complex8_1d(sbuf, dest, stag, rbuf, source, rtag, communicator)
  implicit none
  complex(8),intent(in) :: sbuf(:)
  integer,intent(in)    :: dest
  integer,intent(in)    :: stag
  complex(8),intent(out):: rbuf(:)
  integer,intent(out)   :: source
  integer,intent(in)    :: rtag
  integer,intent(in),optional :: communicator
  integer :: id_comm
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)


  call mpi_sendrecv(sbuf, size(sbuf), MPI_DOUBLE_COMPLEX, dest, stag, &
                    rbuf, size(rbuf), MPI_DOUBLE_COMPLEX, source, rtag, id_comm, &
                    status, ierr)

end subroutine comm_sendrecv_complex8_1d
!-------------------------------------------------------------------------------
subroutine comm_sendrecv_complex8_2d(sbuf, dest, stag, rbuf, source, rtag, communicator)
  implicit none
  complex(8),intent(in) :: sbuf(:,:)
  integer,intent(in)    :: dest
  integer,intent(in)    :: stag
  complex(8),intent(out):: rbuf(:,:)
  integer,intent(out)   :: source
  integer,intent(in)    :: rtag
  integer,intent(in),optional :: communicator
  integer :: id_comm
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  call mpi_sendrecv(sbuf, size(sbuf), MPI_DOUBLE_COMPLEX, dest, stag, &
                    rbuf, size(rbuf), MPI_DOUBLE_COMPLEX, source, rtag, id_comm, &
                    status, ierr)

end subroutine comm_sendrecv_complex8_2d
!-------------------------------------------------------------------------------
subroutine comm_sendrecv_complex8_3d(sbuf, dest, stag, rbuf, source, rtag, communicator)
  implicit none
  complex(8),intent(in) :: sbuf(:,:,:)
  integer,intent(in)    :: dest
  integer,intent(in)    :: stag
  complex(8),intent(out):: rbuf(:,:,:)
  integer,intent(out)   :: source
  integer,intent(in)    :: rtag
  integer,intent(in),optional :: communicator
  integer :: id_comm
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)


  call mpi_sendrecv(sbuf, size(sbuf), MPI_DOUBLE_COMPLEX, dest, stag, &
                    rbuf, size(rbuf), MPI_DOUBLE_COMPLEX, source, rtag, id_comm, &
                    status, ierr)

end subroutine comm_sendrecv_complex8_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
