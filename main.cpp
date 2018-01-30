
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

#include "TPS_F.H"

using namespace amrex;

void toy ()
{
    Box domain(IntVect(0,0,0), IntVect(63,63,63));
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
    
    MultiFab Ex(amrex::convert(ba,Ex_nodal_flag),dm,1,0);
    MultiFab Ey(amrex::convert(ba,Ey_nodal_flag),dm,1,0);
    MultiFab Ez(amrex::convert(ba,Ez_nodal_flag),dm,1,0);
    MultiFab Bx(amrex::convert(ba,Bx_nodal_flag),dm,1,0);
    MultiFab By(amrex::convert(ba,By_nodal_flag),dm,1,0);
    MultiFab Bz(amrex::convert(ba,Bz_nodal_flag),dm,1,0);
    MultiFab jx(amrex::convert(ba,jx_nodal_flag),dm,1,0);
    MultiFab jy(amrex::convert(ba,jy_nodal_flag),dm,1,0);
    MultiFab jz(amrex::convert(ba,jz_nodal_flag),dm,1,0);

    // MPI
    {
        int fcomm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
        tps_mpi_init(fcomm);
    }

    // initialize FFTW plans
    {

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

