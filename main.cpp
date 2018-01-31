
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

#include "TPS_F.H"

using namespace amrex;

void toy ()
{
    Box domain(IntVect(0,0,0), IntVect(31,31,31));
    domain.grow(16);
    BoxArray ba(domain);
    ba.maxSize(32);
    DistributionMapping dm{ba};

    const IntVect Bx_nodal_flag(1,0,0);
    const IntVect By_nodal_flag(0,1,0);
    const IntVect Bz_nodal_flag(0,0,1);
    const IntVect Ex_nodal_flag(0,1,1);
    const IntVect Ey_nodal_flag(1,0,1);
    const IntVect Ez_nodal_flag(1,1,0);
    const IntVect jx_nodal_flag(0,1,1);
    const IntVect jy_nodal_flag(1,0,1);
    const IntVect jz_nodal_flag(1,1,0);

    MultiFab Ex(ba,dm,1,0);
    MultiFab Ey(ba,dm,1,0);
    MultiFab Ez(ba,dm,1,0);
    MultiFab Bx(ba,dm,1,0);
    MultiFab By(ba,dm,1,0);
    MultiFab Bz(ba,dm,1,0);
    MultiFab jx(ba,dm,1,0);
    MultiFab jy(ba,dm,1,0);
    MultiFab jz(ba,dm,1,0);

    // MPI
    {
        int fcomm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
        tps_mpi_init(fcomm);
    }

    // initialize FFTW plans
    {
      tps_fft_init(64,64,64);
    }

    // free MPI
    {
        tps_mpi_finalize();
    }
}


int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    toy();

    amrex::Finalize();
}
