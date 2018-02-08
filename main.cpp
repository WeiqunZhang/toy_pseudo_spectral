
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include "TPS_F.H"

using namespace amrex;

void toy ()
{
    // Since FFTW can do only 1D domain decomposition,
    // for the moment the box is chosen to be long in the x direction
    // With the added guard cells, this will produce a 4 boxes of 32^3
    Box domain(IntVect(0,0,0), IntVect(15,15,111));
    RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                 {AMREX_D_DECL( 1.0, 1.0, 1.0)});
    domain.grow(8);
    BoxArray ba(domain);
    ba.maxSize(32);
    DistributionMapping dm{ba};

    int N_steps = 64;
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

    // initialize MPI
    {
        int fcomm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
        tps_mpi_init(fcomm);
    }

    // initialize FFTW plans
    {
      int count = 0;
      // Loop through the grids
      for ( MFIter mfi(Ex); mfi.isValid(); ++mfi ) {
   	// Make sure that there is only 1 grid per MPI
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE( count < 1,
		  "Only one grid per MPI is allowed" );
	count ++;
	// Initialize fft
	tps_fft_init( BL_SPACEDIM,
		      domain.loVect(), domain.hiVect(),
		      Ex[mfi].loVect(), Ex[mfi].hiVect(),
		      Ex[mfi].dataPtr(), Ey[mfi].dataPtr(), Ez[mfi].dataPtr(),
		      Bx[mfi].dataPtr(), By[mfi].dataPtr(), Bz[mfi].dataPtr(),
		      jx[mfi].dataPtr(), jy[mfi].dataPtr(), jz[mfi].dataPtr(),
		      rho[mfi].dataPtr(), rhoold[mfi].dataPtr() );
      }
    }

    // Loop over iterations
    for ( int i_step=0; i_step<N_steps; i_step++ ){
      std::cout << "Step " << i_step << "/" << N_steps << std::endl;


      // Write plotfile
      {
	const std::string& pfname = amrex::Concatenate("./data/plt",i_step);
	const Vector<std::string>& varnames {"x"};
	std::array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1,1,1)};
	amrex::WriteSingleLevelPlotfile (pfname, Ex, varnames,
		     Geometry(domain, &real_box, 0, is_periodic.data()),
					 0., 0);
      }
      
      // Push the E and B fields
      push_psatd_ebfield_3d_();
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
