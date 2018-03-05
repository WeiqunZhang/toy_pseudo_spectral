#ifndef TPS_F_H_
#define TPS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

  void tps_mpi_init (int comm);
  void tps_mpi_finalize ();

  void tps_fft_init(
    const int* global_lo, const int* global_hi,
    const int* local_lo, const int* local_hi,
    void* fft_data);

  void tps_fft_free ();

  void tps_push_eb(
    amrex_real* ex, const int* exlo, const int* exhi,
    amrex_real* ey, const int* eylo, const int* eyhi,
    amrex_real* ez, const int* ezlo, const int* ezhi,
    amrex_real* bx, const int* bxlo, const int* bxhi,
    amrex_real* by, const int* bylo, const int* byhi,
    amrex_real* bz, const int* bzlo, const int* bzhi,
    amrex_real* jx, const int* jxlo, const int* jxhi,
    amrex_real* jy, const int* jylo, const int* jyhi,
    amrex_real* jz, const int* jzlo, const int* jzhi,
    amrex_real* rho, const int* r1lo, const int* r1hi,
    amrex_real* rhoold, const int* r2lo, const int* r2hi);

#ifdef __cplusplus
}
#endif

#endif
