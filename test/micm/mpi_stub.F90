! Single-rank stub of the F77 MPI broadcast interface used by mo_micm.F90.
! Compiled separately so callers with character/logical buffers link against
! the same implicit-interface external.
subroutine mpi_bcast(buffer, count, datatype, root, comm, ierr)
  implicit none
  integer :: buffer(*)
  integer :: count, datatype, root, comm, ierr
  ierr = 0
end subroutine mpi_bcast
