
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

  SUBROUTINE tps_fft_init( dim, global_lo, global_hi, local_lo, local_hi, &
       ex_wrpx, ey_wrpx, ez_wrpx, bx_wrpx, by_wrpx, bz_wrpx, &
       jx_wrpx, jy_wrpx, jz_wrpx, rho_wrpx, rhoold_wrpx ) &
       BIND(C,name='tps_fft_init')

    USE shared_data, only: rank, c_dim, p3dfft_flag, &
         fftw_with_mpi, fftw_threads_ok, fftw_hybrid, fftw_mpi_transpose, &
         rho, rhoold, &
         nx_global, ny_global, nz_global, & ! Size of global FFT
         nx, ny, nz, & ! Size of local subdomains
         nkx, nky, nkz, & ! Size of local ffts
         dx, dy, dz, &
         iz_min_r, iz_max_r, iy_min_r, iy_max_r, ix_min_r, ix_max_r !Loop bounds
    USE gpstd_solver, only: init_gpstd
    USE constants, only: num, clight
    USE picsar_precision, only: idp
    USE params, only: dt
    USE fields, only: nxguards, nyguards, nzguards, & ! Size of guard regions
         ex, ey, ez, bx, by, bz, jx, jy, jz, &
         ex_r, ey_r, ez_r, bx_r, by_r, bz_r, &
         jx_r, jy_r, jz_r, rho_r, rhoold_r, &
         exf, eyf, ezf, bxf, byf, bzf, &
         jxf, jyf, jzf, rhof, rhooldf, &
         l_spectral, l_staggered, norderx, nordery, norderz
#if defined(FFTW)
    USE mpi_fftw3, only: local_nz
    USE fourier_psaotd, only: init_plans_blocks
#endif
    IMPLICIT NONE

    integer, intent(in) :: global_lo(BL_SPACEDIM), global_hi(BL_SPACEDIM)
    integer, intent(in) :: local_lo(BL_SPACEDIM), local_hi(BL_SPACEDIM)
    integer, value, intent(in) :: dim
    integer :: nx_padded

    REAL(num), INTENT(INOUT), TARGET, &
         DIMENSION(0:local_hi(1)-local_lo(1), &
                   0:local_hi(2)-local_lo(2), &
                   0:local_hi(3)-local_lo(3)) :: &
                   ex_wrpx, ey_wrpx, ez_wrpx, &
                   bx_wrpx, by_wrpx, bz_wrpx, &
                   jx_wrpx, jy_wrpx, jz_wrpx, rho_wrpx, rhoold_wrpx

!    CALL DFFTW_INIT_THREADS(iret)
    !    fftw_threads_ok = .TRUE.
    PRINT *, rank, local_lo(1), local_lo(2), local_lo(3)
    PRINT *, rank, local_hi(1), local_hi(2), local_hi(3)
    PRINT *, rank, local_hi(1)-local_lo(1) + 1, &
         local_hi(2)-local_lo(2) + 1, &
         local_hi(3)-local_lo(3) + 1

    ! Define size of domains: necessary for the initialization of the global FFT
    nx_global = INT(global_hi(1)-global_lo(1)+1,idp)
    ny_global = INT(global_hi(2)-global_lo(2)+1,idp)
    nz_global = INT(global_hi(3)-global_lo(3)+1,idp)
    nx = INT(local_hi(1)-local_lo(1)+1,idp)
    ny = INT(local_hi(2)-local_lo(2)+1,idp)
    nz = INT(local_hi(3)-local_lo(3)+1,idp)
    local_nz = nz
    ! No need to distinguish physical and guard cells for the global FFT;
    ! only nx+2*nxguards counts. Thus we declare 0 guard cells for simplicity
    nxguards = 0_idp
    nyguards = 0_idp
    nzguards = 0_idp

    ! For the calculation of the modified [k] vectors
    l_staggered = .TRUE.
    norderx = 8_idp
    nordery = 12_idp
    norderz = 16_idp
    dx = 1./nx_global
    dy = 1./ny_global
    dz = 1./nz_global
    dt = 1./60/clight ! Take advantage of the fact that there is no CFL
    ! Define parameters of FFT plans
    c_dim = INT(dim,idp)   ! Dimensionality of the simulation (2d/3d)
    fftw_with_mpi = .TRUE. ! Activate MPI FFTW
    fftw_hybrid = .FALSE.   ! FFT per MPI subgroup (instead of global)
    fftw_mpi_transpose = .FALSE. ! Do not transpose the data
    fftw_threads_ok = .FALSE.   ! Do not use threads for FFTW
    p3dfft_flag = .FALSE.
    l_spectral  = .TRUE.   ! Activate spectral Solver, using FFT

    ! Point the fields in the PICSAR modules to the fields provided by WarpX
    ex => ex_wrpx
    ey => ey_wrpx
    ez => ez_wrpx
    bx => bx_wrpx
    by => by_wrpx
    bz => bz_wrpx
    jx => jx_wrpx
    jy => jy_wrpx
    jz => jz_wrpx
    rho => rho_wrpx
    rhoold => rhoold_wrpx
    ! Allocate padded arrays for MPI FFTW
    nx_padded = 2*(nx/2 + 1)
    ALLOCATE(ex_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(ey_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(ez_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(bx_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(by_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(bz_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(jx_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(jy_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(jz_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(rho_r(1:nx_padded, 1:ny, 1:nz))
    ALLOCATE(rhoold_r(1:nx_padded, 1:ny, 1:nz))
    ! Set array bounds when copying ex to ex_r in PICSAR
    ix_min_r = 1; ix_max_r = nx
    iy_min_r = 1; iy_max_r = ny
    iz_min_r = 1; iz_max_r = nz
    ! Allocate Fourier space fields of the same size
    nkx = nx/2 + 1
    nky = ny
    nkz = nz
    ALLOCATE(exf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(eyf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(ezf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(bxf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(byf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(bzf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(jxf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(jyf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(jzf(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(rhof(1:nkx, 1:nky, 1:nkz))
    ALLOCATE(rhooldf(1:nkx, 1:nky, 1:nkz))

    CALL init_plans_blocks

  END SUBROUTINE tps_fft_init


  ! _________________________________________________________________________
  !> @brief
  !> Routine that initializes the fields provided by WarpX with a dirac in Ex
  !
  ! __________________________________________________________________________
  SUBROUTINE initialize_fields( ex, global_lo, global_hi, local_lo, local_hi ) &
       bind(c, name='initialize_fields')

    USE constants, only: num
    integer, intent(in) :: global_lo(BL_SPACEDIM), global_hi(BL_SPACEDIM)
    integer, intent(in) :: local_lo(BL_SPACEDIM), local_hi(BL_SPACEDIM)
    REAL(num), INTENT(INOUT), &
         DIMENSION(local_lo(1):local_hi(1), &
                   local_lo(2):local_hi(2), &
                   local_lo(3):local_hi(3)) :: ex
    INTEGER :: ix, iy, iz
    
    
    ix = global_lo(1)/2
    iy = global_lo(2)/2
    iz = global_lo(3)/2

    IF ( (ix>=local_lo(1)) .AND. (ix<=local_hi(1)) .AND. &
         (iy>=local_lo(2)) .AND. (iy<=local_hi(2)) .AND. &
         (iz>=local_lo(3)) .AND. (iz<=local_hi(3))) THEN
       ex(ix, iy, iz) = 1
    ENDIF

  END SUBROUTINE initialize_fields
    

  subroutine tps_mpi_finalize () bind(c,name='tps_mpi_finalize')
    use shared_data, only : comm
    use mpi_fftw3, only : fftw_mpi_cleanup
    integer :: ierr
    call fftw_mpi_cleanup()
    call mpi_comm_free(comm, ierr)
  end subroutine tps_mpi_finalize


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
