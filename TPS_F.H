#ifndef TPS_F_H_
#define TPS_F_H_

#ifdef __cplusplus
extern "C" {
#endif

  void tps_mpi_init (int comm);
  void tps_mpi_finalize ();
  void tps_fft_init( int dim,
    const int* global_lo, const int* global_hi,
    const int* local_lo, const int* local_hi,
    Real* ex, Real* ey, Real* ez,
    Real* bx, Real* by, Real* bz,
    Real* jx, Real* jy, Real* jz,
    Real* rho, Real* rhoold );
  void tps_push_eb(
    const int* local_lo, const int* local_hi,
    Real* ex, Real* ey, Real* ez,
    Real* bx, Real* by, Real* bz,
    Real* jx, Real* jy, Real* jz,
    Real* rho, Real* rhoold );
  void push_psatd_ebfield_3d_();
  void initialize_fields( Real* ex,
    const int* global_lo, const int* global_hi,
    const int* local_lo, const int* local_hi );

#ifdef __cplusplus
}
#endif

#endif
