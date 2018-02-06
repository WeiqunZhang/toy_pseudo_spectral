
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

#include "TPS_F.H"

using namespace amrex;

void toy ()
{
    // Since FFTW can do only 1D domain decomposition,
    // for the moment the box is chosen to be long in the x direction
    // With the added guard cells, this will produce a 4 boxes of 32^3
    Box domain(IntVect(0,0,0), IntVect(95,15,15));
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

    MultiFab Ex(ba,dm,1,0);  Ex.setVal(0.);
    MultiFab Ey(ba,dm,1,0);  Ey.setVal(0.);
    MultiFab Ez(ba,dm,1,0);  Ez.setVal(0.);
    MultiFab Bx(ba,dm,1,0);  Bx.setVal(0.);
    MultiFab By(ba,dm,1,0);  By.setVal(0.);
    MultiFab Bz(ba,dm,1,0);  Bz.setVal(0.);
    MultiFab jx(ba,dm,1,0);  jx.setVal(0.);
    MultiFab jy(ba,dm,1,0);  jy.setVal(0.);
    MultiFab jz(ba,dm,1,0);  jz.setVal(0.);
    MultiFab rho(ba,dm,1,0);  rho.setVal(0.);
    MultiFab rhoold(ba,dm,1,0);  rhoold.setVal(0.);

    // MPI
    {
        int fcomm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
        tps_mpi_init(fcomm);
    }

    // initialize FFTW plans
    {
      // Loop through the grids
      for ( MFIter mfi(Ex); mfi.isValid(); ++mfi ) {
	tps_fft_init( BL_SPACEDIM,
		      domain.loVect(), domain.hiVect(),
		      Ex[mfi].loVect(), Ex[mfi].hiVect(),
		      Ex[mfi].dataPtr(), Ey[mfi].dataPtr(), Ez[mfi].dataPtr(),
		      Bx[mfi].dataPtr(), By[mfi].dataPtr(), Bz[mfi].dataPtr(),
		      jx[mfi].dataPtr(), jy[mfi].dataPtr(), jz[mfi].dataPtr(),
		      rho[mfi].dataPtr(), rhoold[mfi].dataPtr() );
	          }
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
