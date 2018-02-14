
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
    // With the added guard cells, this will produce a n_mpi boxes of N^3
    int N=32;
    int nguards=4;
    int nmpi=4;
    Box domain(IntVect(nguards,nguards,nguards),
	       IntVect(N-nguards-1,N-nguards-1,nmpi*N-nguards-1));
    domain.grow(nguards);
    BoxArray ba(domain);
    ba.maxSize(N);
    DistributionMapping dm{ba};
    RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},
		     {AMREX_D_DECL( 1.0*N, 1.0*N, 1.0*nmpi*N)});
    Geometry geom(domain, &real_box, 0);
    
    int N_steps = 16;
    MultiFab Ex(ba,dm,1,0);  Ex.setVal(0.);
    MultiFab Ey(ba,dm,1,0);  Ey.setVal(0.);
    MultiFab Ez(ba,dm,1,0);  Ez.setVal(0.);
    MultiFab Bx(ba,dm,1,0);  Bx.setVal(0.);
    MultiFab By(ba,dm,1,0);  By.setVal(0.);
    MultiFab Bz(ba,dm,1,0);  Bz.setVal(0.);
    MultiFab jx(ba,dm,1,0);  jx.setVal(0.);
    MultiFab jy(ba,dm,1,0);  jy.setVal(0.);
    MultiFab jz(ba,dm,1,0);  jz.setVal(0.);
    MultiFab rho1(ba,dm,1,0); rho1.setVal(0.);
    MultiFab rho2(ba,dm,1,0); rho2.setVal(0.);

    // initialize MPI
    {
        int fcomm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
        tps_mpi_init(fcomm);
    }

    // initialize field and FFTW plans
    {
      int count = 0;
      // Loop through the grids
      for ( MFIter mfi(Ex); mfi.isValid(); ++mfi ) {
   	// Make sure that there is only 1 grid per MPI
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE( count < 1,
		  "Only one grid per MPI is allowed" );
	count ++;
	// Initialize fields
	initialize_fields( Ex[mfi].dataPtr(),
			   domain.loVect(), domain.hiVect(),
			   Ex[mfi].loVect(), Ex[mfi].hiVect() );
	// Initialize fft
	tps_fft_init( BL_SPACEDIM,
		      domain.loVect(), domain.hiVect(),
		      Ex[mfi].loVect(), Ex[mfi].hiVect(),
		      Ex[mfi].dataPtr(), Ey[mfi].dataPtr(), Ez[mfi].dataPtr(),
		      Bx[mfi].dataPtr(), By[mfi].dataPtr(), Bz[mfi].dataPtr(),
		      jx[mfi].dataPtr(), jy[mfi].dataPtr(), jz[mfi].dataPtr(),
		      rho1[mfi].dataPtr(), rho2[mfi].dataPtr() );
      }
    }

    // Loop over iterations
    for ( int i_step=0; i_step<N_steps; i_step++ ){
      std::cout << "Step " << i_step << "/" << N_steps << std::endl;


      // Write plotfile
      {
	const std::string& pfname = amrex::Concatenate("./data/plt",i_step);
	amrex::WriteSingleLevelPlotfile (pfname, Ex, {"Ex"}, geom, 0., 0 );
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
