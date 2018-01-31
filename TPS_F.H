#ifndef TPS_F_H_
#define TPS_F_H_

#ifdef __cplusplus
extern "C" {
#endif

    void tps_mpi_init (int comm);
    void tps_mpi_finalize ();
  void tps_fft_init( int nx, int ny, int nz );

#ifdef __cplusplus
}
#endif

#endif
