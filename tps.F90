
module tps

  use, intrinsic :: iso_c_binding
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
  !> Routine that pushes the fields in spectral space
  ! __________________________________________________________________________

  SUBROUTINE tps_push_eb( &
       ex_wrpx, exlo, exhi, &
       ey_wrpx, eylo, eyhi, &
       ez_wrpx, ezlo, ezhi, &
       bx_wrpx, bxlo, bxhi, &
       by_wrpx, bylo, byhi, &
       bz_wrpx, bzlo, bzhi, &
       jx_wrpx, jxlo, jxhi, &
       jy_wrpx, jylo, jyhi, &
       jz_wrpx, jzlo, jzhi, &
       rho_wrpx, r1lo, r1hi, &
       rhoold_wrpx, r2lo, r2hi ) &
       BIND(C,name='tps_push_eb')       
       
    USE fields, only: ex, ey, ez, bx, by, bz, jx, jy, jz
    USE shared_data, only: rhoold, rho
    USE constants, only: num
    implicit none
    
    integer, dimension(3), intent(in) :: exlo, exhi, eylo, eyhi, ezlo, ezhi, bxlo, bxhi, &
         bylo, byhi, bzlo, bzhi, jxlo, jxhi, jylo, jyhi, jzlo, jzhi, r1lo, r1hi, r2lo, r2hi
    REAL(num), INTENT(INOUT), TARGET :: ex_wrpx(0:exhi(1)-exlo(1),0:exhi(2)-exlo(2),0:exhi(3)-exlo(3))
    REAL(num), INTENT(INOUT), TARGET :: ey_wrpx(0:eyhi(1)-eylo(1),0:eyhi(2)-eylo(2),0:eyhi(3)-eylo(3))
    REAL(num), INTENT(INOUT), TARGET :: ez_wrpx(0:ezhi(1)-ezlo(1),0:ezhi(2)-ezlo(2),0:ezhi(3)-ezlo(3))
    REAL(num), INTENT(INOUT), TARGET :: bx_wrpx(0:bxhi(1)-bxlo(1),0:bxhi(2)-bxlo(2),0:bxhi(3)-bxlo(3))
    REAL(num), INTENT(INOUT), TARGET :: by_wrpx(0:byhi(1)-bylo(1),0:byhi(2)-bylo(2),0:byhi(3)-bylo(3))
    REAL(num), INTENT(INOUT), TARGET :: bz_wrpx(0:bzhi(1)-bzlo(1),0:bzhi(2)-bzlo(2),0:bzhi(3)-bzlo(3))
    REAL(num), INTENT(INOUT), TARGET :: jx_wrpx(0:jxhi(1)-jxlo(1),0:jxhi(2)-jxlo(2),0:jxhi(3)-jxlo(3))
    REAL(num), INTENT(INOUT), TARGET :: jy_wrpx(0:jyhi(1)-jylo(1),0:jyhi(2)-jylo(2),0:jyhi(3)-jylo(3))
    REAL(num), INTENT(INOUT), TARGET :: jz_wrpx(0:jzhi(1)-jzlo(1),0:jzhi(2)-jzlo(2),0:jzhi(3)-jzlo(3))
    REAL(num), INTENT(INOUT), TARGET :: rho_wrpx(0:r1hi(1)-r1lo(1),0:r1hi(2)-r1lo(2),0:r1hi(3)-r1lo(3))
    REAL(num), INTENT(INOUT), TARGET :: rhoold_wrpx(0:r2hi(1)-r2lo(1),0:r2hi(2)-r2lo(2),0:r2hi(3)-r2lo(3))

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
  
    ! Call the corresponding PICSAR function
    CALL push_psatd_ebfield_3d()
    
    ex => null()
    ey => null()
    ez => null()
    bx => null()
    by => null()
    bz => null()
    jx => null()
    jy => null()
    jz => null()
    rho => null()
    rhoold => null()

  END SUBROUTINE
    
  
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

  SUBROUTINE tps_fft_init( global_lo, global_hi, local_lo, local_hi, &
       dx_wrpx, dy_wrpx, dz_wrpx, dt_wrpx) &
       BIND(C,name='tps_fft_init')

    USE shared_data, only: rank, comm, c_dim, p3dfft_flag, &
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
    USE mpi_fftw3, only: local_nz, local_z0, fftw_mpi_local_size_3d, alloc_local
    USE fourier_psaotd, only: init_plans_blocks
#endif
    IMPLICIT NONE

    integer, intent(in) :: global_lo(BL_SPACEDIM), global_hi(BL_SPACEDIM)
    integer, intent(in) :: local_lo(BL_SPACEDIM), local_hi(BL_SPACEDIM)
    REAL(num), intent(in) :: dx_wrpx, dy_wrpx, dz_wrpx, dt_wrpx
    integer :: nx_padded
    
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
    ! No need to distinguish physical and guard cells for the global FFT;
    ! only nx+2*nxguards counts. Thus we declare 0 guard cells for simplicity
    nxguards = 0_idp
    nyguards = 0_idp
    nzguards = 0_idp
    ! Find the decomposition that FFTW imposes in kspace
    alloc_local = fftw_mpi_local_size_3d( &
         INT(nz_global,C_INTPTR_T), &
         INT(ny_global,C_INTPTR_T), &
         INT(nx_global,C_INTPTR_T)/2+1, &
         comm, local_nz, local_z0)
    IF (local_nz .NE. nz) THEN
       PRINT *, 'ERRROR'
    ENDIF
    
    ! For the calculation of the modified [k] vectors
    l_staggered = .TRUE.
    norderx = 2_idp
    nordery = 2_idp
    norderz = 2_idp
    dx = dx_wrpx;
    dy = dy_wrpx;
    dz = dz_wrpx;
    dt = dt_wrpx;    
    ! Define parameters of FFT plans
    c_dim = INT(AMREX_SPACEDIM,idp)   ! Dimensionality of the simulation (2d/3d)
    fftw_with_mpi = .TRUE. ! Activate MPI FFTW
    fftw_hybrid = .FALSE.   ! FFT per MPI subgroup (instead of global)
    fftw_mpi_transpose = .FALSE. ! Do not transpose the data
    fftw_threads_ok = .FALSE.   ! Do not use threads for FFTW
    p3dfft_flag = .FALSE.
    l_spectral  = .TRUE.   ! Activate spectral Solver, using FFT

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
  SUBROUTINE initialize_fields( field, global_lo, global_hi, local_lo, local_hi ) &
       bind(c, name='initialize_fields')

    USE constants, only: num
    integer, intent(in) :: global_lo(BL_SPACEDIM), global_hi(BL_SPACEDIM)
    integer, intent(in) :: local_lo(BL_SPACEDIM), local_hi(BL_SPACEDIM)
    REAL(num), INTENT(INOUT), &
         DIMENSION(local_lo(1):local_hi(1), &
                   local_lo(2):local_hi(2), &
                   local_lo(3):local_hi(3)) :: field
    INTEGER :: ix, iy, iz, ix0, iy0, iz0, i, j, k
    
    
    ix0 = (global_hi(1)+global_lo(1))/2
    iy0 = (global_hi(2)+global_lo(2))/2
    iz0 = (global_hi(3)+global_lo(3))/2

    DO i=-1,1
       DO j=-1, 1
          DO k=-1, 1
             ix = ix0 + i
             iy = iy0 + j
             iz = iz0 + k
             IF ( (ix>=local_lo(1)) .AND. (ix<=local_hi(1)) .AND. &
                  (iy>=local_lo(2)) .AND. (iy<=local_hi(2)) .AND. &
                  (iz>=local_lo(3)) .AND. (iz<=local_hi(3))) THEN
                field(ix, iy, iz) = (1.-0.5*ABS(i)) * &
                     (1.-0.5*ABS(j)) * &
                     (1.-0.5*ABS(k))
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    

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
