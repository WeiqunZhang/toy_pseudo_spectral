
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


  ! _________________________________________________________________________
  !> @brief
  !> Routine that links the AMREX fields to the PICSAR modules,
  !> creates the FFT plans, and creates the blocks of coefficients
  !
  !> @params[in] nx, ny, nz - INTEGER - number of cells along each direction
  !> in the FFT subgroup
  !> @params[in] dim - INTEGER - dimensionality of the simulation
  !> (2 for 2d, 3 for 3d)
  !
  ! __________________________________________________________________________

  SUBROUTINE tps_fft_init( nx, ny, nz, dim, &
       ex, ey, ez, bx, by, bz, jx, jy, jz, rho, rhoold ) &
       BIND(C,name='tps_fft_init')

    USE shared_data, only: rank, fftw_with_mpi, p3dfft_flag, fftw_threads_ok, &
                      nx_global, ny_global, nz_global, c_dim
    USE constants, only: num
    USE picsar_precision, only: idp
    USE fields, only: l_spectral, ex_r, ey_r, ez_r, bx_r, by_r, bz_r, &
        jx_r, jy_r, jz_r, rho_r, rhoold_r
!    USE fastfft
#if defined(FFTW)
    USE fourier_psaotd
!    USE fourier
#endif
    IMPLICIT NONE

    integer, value, intent(in) :: nx, ny, nz, dim
    REAL(num), INTENT(INOUT), TARGET, DIMENSION(1:32,1:32,1:32) :: &
         ex, ey, ez, bx, by, bz, jx, jy, jz, rho, rhoold

    IF(rank==0) PRINT*, 'BEGIN INIT EXTERNAL'
    l_spectral  = .TRUE.   ! Activate spectral Solver, using FFT

    ! Initialize FFT plans
    fftw_with_mpi = .TRUE. ! Activate MPI FFTW
    c_dim = INT(dim,idp)   ! Dimensionality of the simulation (2d/3d)
!    fftw_hybrid = .FALSE.
!    fftw_mpi_transpose = .FALSE.
!    l_staggered = .TRUE.
    fftw_threads_ok = .FALSE.
!    CALL DFFTW_INIT_THREADS(iret)
!    fftw_threads_ok = .TRUE.
    p3dfft_flag = .FALSE.
!    p3dfft_stride = .FALSE.
!    c_dim = INT(cdim,idp)
    nx_global = INT(nx,idp)
    ny_global = INT(ny,idp)
    nz_global = INT(nz,idp)

    ex_r => ex
    ey_r => ey
    ez_r => ez
    bx_r => bx
    by_r => by
    bz_r => bz
    jx_r => jx
    jy_r => jy
    jz_r => jz
    rho_r => rho
    rhoold_r => rhoold

!      nkx=(2*nxguards+nx+1)/2+1! Real To Complex Transform
!      nky=(2*nyguards+ny+1)
!      nkz=(2*nzguards+nz+1)

!      IF(.NOT. ASSOCIATED(exf)) ALLOCATE(exf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(eyf)) ALLOCATE(eyf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(ezf)) ALLOCATE(ezf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(bxf)) ALLOCATE(bxf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(byf)) ALLOCATE(byf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(bzf)) ALLOCATE(bzf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(jxf)) ALLOCATE(jxf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(jyf)) ALLOCATE(jyf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(jzf)) ALLOCATE(jzf(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(rhof)) ALLOCATE(rhof(nkx, nky, nkz))
!      IF(.NOT. ASSOCIATED(rhooldf)) ALLOCATE(rhooldf(nkx, nky, nkz))


!    CALL init_plans_blocks

    IF(rank==0) PRINT*, 'END INIT EXTERNAL'
  END SUBROUTINE tps_fft_init


  subroutine tps_mpi_finalize () bind(c,name='tps_mpi_finalize')
    use shared_data, only : comm
    use mpi_fftw3, only : fftw_mpi_cleanup
    integer :: ierr
    call fftw_mpi_cleanup()
    call mpi_comm_free(comm, ierr)
  end subroutine tps_mpi_finalize


  subroutine tps_init_plans_blocks () bind(c,name='init_plans_blocks')
    ! picsar: select_case_dims_local
    ! picsar: init_gpstd


  end subroutine tps_init_plans_blocks


  subroutine tps_push_ebfield (cp,nvar)
    use iso_c_binding, only : c_ptr, c_f_pointer
    use matrix_data, only : ns_max
    use matrix_coefficients, only : block3d, matrix_blocks, cc_mat
    integer, intent(in) :: nvar
    type(c_ptr), intent(in) :: cp(nvar,nvar)

    type(matrix_blocks), pointer, dimension(:) :: my_cc_mat
    integer :: nfftx, nffty, nfftz, i,j

    allocate(my_cc_mat(ns_max))
    allocate(my_cc_mat(1)%block_matrix2d(nvar,nvar))
    my_cc_mat(1)%nblocks = nvar
    ! nfftx, nffty, nfftz = ...
    do j = 1, nvar
       do i = 1, nvar
          call c_f_pointer(cp(i,j), my_cc_mat(1)%block_matrix2d(i,j)%block3dc, &
               shape=[nfftx/2,nffty,nfftz])
          my_cc_mat(1)%block_matrix2d(i,j)%nx = nfftx/2+1
          my_cc_mat(1)%block_matrix2d(i,j)%ny = nffty
          my_cc_mat(1)%block_matrix2d(i,j)%nz = nfftz
       end do
    end do
    cc_mat => my_cc_mat

    ! .....

    nullify(cc_mat)
    deallocate(my_cc_mat(1)%block_matrix2d)
    deallocate(my_cc_mat)
  end subroutine tps_push_ebfield

end module tps
