
module tps

  implicit none

contains

  subroutine tps_mpi_init (comm_in) bind(c,name='tps_mpi_init')
    use shared_data, only : comm, rank, nproc
    integer, value, intent(in) :: comm_in
    integer :: ierr, lnproc, lrank
    call mpi_comm_dup(comm_in, comm, ierr)
    call mpi_comm_size(comm, lnproc, ierr)
    nproc = lnproc
    call mpi_comm_rank(comm, lrank, ierr)
    rank = lrank
  end subroutine tps_mpi_init


  subroutine tps_mpi_finalize () bind(c,name='tps_mpi_finalize')
    use shared_data, only : comm
    integer :: ierr
    call mpi_comm_free(comm, ierr)
  end subroutine tps_mpi_finalize

end module tps
