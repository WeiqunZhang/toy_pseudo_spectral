
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#include "TPS_F.H"

using namespace amrex;

namespace
{
    Vector<int> n_cell({128, 128, 128});
    int max_grid_size = 32;
    int max_step = 20;
    int plot_int = 1;
    int nguards_fft = 10;
    int ngroups_fft = 4;
    
    IntVect Bx_nodal_flag(1,0,0);
    IntVect By_nodal_flag(0,1,0);
    IntVect Bz_nodal_flag(0,0,1);
    
    IntVect Ex_nodal_flag(0,1,1);
    IntVect Ey_nodal_flag(1,0,1);
    IntVect Ez_nodal_flag(1,1,0);
    
    IntVect jx_nodal_flag(0,1,1);
    IntVect jy_nodal_flag(1,0,1);
    IntVect jz_nodal_flag(1,1,0);
}

void write_plotfile (const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                     const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                     const Geometry& geom, int istep);

void copy_data_from_fft_to_valid (MultiFab& mf, const MultiFab& mf_fft,
                                  const Geometry& geom, const BoxArray& ba_valid_fft);

void toy ()
{
    static_assert(AMREX_SPACEDIM == 3, "3d only");

    // parameters
    {
        ParmParse pp;
        pp.queryarr("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("max_step", max_step);
        pp.query("plot_int", plot_int);
        pp.query("nguards_fft", nguards_fft);
        pp.query("ngroups_fft", ngroups_fft);
    }

    Geometry geom;
    {
        Box domain(IntVect::TheZeroVector(), IntVect(AMREX_D_DECL(n_cell[0]-1,
                                                                  n_cell[1]-1,
                                                                  n_cell[2]-1)));
        RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},
                         {AMREX_D_DECL(1.0,1.0,1.0)});
        geom.define(domain, &real_box, 0);
    }

    MultiFab Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, rho1, rho2;
    {
        BoxArray ba(geom.Domain());
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        const int ng = 3;
        Ex.define(amrex::convert(ba,Ex_nodal_flag),dm,1,ng);
        Ey.define(amrex::convert(ba,Ey_nodal_flag),dm,1,ng);
        Ez.define(amrex::convert(ba,Ez_nodal_flag),dm,1,ng);
        Bx.define(amrex::convert(ba,Bx_nodal_flag),dm,1,ng);
        By.define(amrex::convert(ba,By_nodal_flag),dm,1,ng);
        Bz.define(amrex::convert(ba,Bz_nodal_flag),dm,1,ng);
        jx.define(amrex::convert(ba,jx_nodal_flag),dm,1,ng);
        jy.define(amrex::convert(ba,jy_nodal_flag),dm,1,ng);
        jz.define(amrex::convert(ba,jz_nodal_flag),dm,1,ng);
        Ex.setVal(0.0);
        Ey.setVal(0.0);
        Ez.setVal(0.0);
        Bx.setVal(0.0);
        By.setVal(0.0);
        Bz.setVal(0.0);
        jx.setVal(0.0);
        jy.setVal(0.0);
        jz.setVal(0.0);

        const IntVect center = (geom.Domain().smallEnd() + geom.Domain().bigEnd()) / 2;
        for (MFIter mfi(Ex); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = Ex[mfi];
            const Box& bx = fab.box();
            for (int k = -1; k <= 1; ++k) {
                for (int j = -1; j <= 1; ++j) {
                    for (int i = -1; i <= 1; ++i) {
                        IntVect iv(AMREX_D_DECL(center[0]+i,center[1]+j,center[2]+k));
                        if (bx.contains(iv)) {
                            fab(iv) = AMREX_D_TERM(  (1.-0.5*std::abs(i)),
                                                   * (1.-0.5*std::abs(j)),
                                                   * (1.-0.5*std::abs(k)));
                        }
                    }
                }
            }
        }
    }

    if (plot_int > 0) {
        write_plotfile(Ex,Ey,Ez,Bx,By,Bz,geom,0);
    }

    
    // initialize MPI
    MPI_Comm comm_fft;  // subcommunicator for FFT
    int np_fft;         // # of processes in the subcommunicator 
    int rank_fft;       // my rank in my subcommunicator
    int color_fft;      // my color in ngroups_fft subcommunicators.  0 <= color_fft < ngroups_fft
    {
        int nprocs = ParallelDescriptor::NProcs();
        ngroups_fft = std::min(ngroups_fft, nprocs);
        np_fft = nprocs / ngroups_fft;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(np_fft*ngroups_fft == nprocs,
                                         "Number of processes must be divisilbe by number of FFT groups");

        int myproc = ParallelDescriptor::MyProc();
        color_fft = myproc / np_fft;
        MPI_Comm_split(ParallelDescriptor::Communicator(), color_fft, myproc, &comm_fft);

        MPI_Comm_rank(comm_fft, &rank_fft);

        int fcomm = MPI_Comm_c2f(comm_fft);
        tps_mpi_init(fcomm);
    }

    MultiFab Ex_fft, Ey_fft, Ez_fft, Bx_fft, By_fft, Bz_fft, jx_fft, jy_fft, jz_fft, rho1_fft, rho2_fft;
    BoxArray ba_valid_fft; // This is not the one used for builing Ex_fft etc.  It doesn't contain ghost cells.
    Box domain_fft;     // the "global" domain for the subdomain FFT
    {
        BoxList bl_fft;
        bl_fft.reserve(ParallelDescriptor::NProcs());

        BoxList bl(geom.Domain(), ngroups_fft);
        AMREX_ALWAYS_ASSERT(bl.size() == ngroups_fft);
        const Vector<Box>& bldata = bl.data();
        for (int igroup = 0; igroup < ngroups_fft; ++igroup)
        {
            const Box& bx = amrex::grow(bldata[igroup], nguards_fft);
            // chop in z-direction into np_fft for FFTW
            BoxList tbl(bx, np_fft, Direction::z);
            bl_fft.join(tbl);
            if (igroup == color_fft) {
                domain_fft = bx;
            }
        }

        BoxArray ba_fft(std::move(bl_fft));
        AMREX_ALWAYS_ASSERT(ba_fft.size() == ParallelDescriptor::NProcs());
        Vector<int> pmap(ba_fft.size());
        std::iota(pmap.begin(), pmap.end(), 0);
        DistributionMapping dm_fft(pmap);

        const Box foobox(IntVect(AMREX_D_DECL(-nguards_fft-2,-nguards_fft-2,-nguards_fft-2)),
                         IntVect(AMREX_D_DECL(-nguards_fft-2,-nguards_fft-2,-nguards_fft-2)));

        BoxList bl_valid;
        bl_valid.reserve(ba_fft.size());
        for (int i = 0; i < ba_fft.size(); ++i)
        {
            int igroup = i / np_fft;
            const Box& bx = ba_fft[i] & bldata[igroup];
            if (bx.ok())
            {
                bl_valid.push_back(bx);
            }
            else
            {
                bl_valid.push_back(foobox);
            }
        }

        ba_valid_fft.define(std::move(bl_valid));

        Ex_fft.define(amrex::convert(ba_fft,Ex_nodal_flag), dm_fft, 1, 0);
        Ey_fft.define(amrex::convert(ba_fft,Ey_nodal_flag), dm_fft, 1, 0);
        Ez_fft.define(amrex::convert(ba_fft,Ez_nodal_flag), dm_fft, 1, 0);
        Bx_fft.define(amrex::convert(ba_fft,Bx_nodal_flag), dm_fft, 1, 0);
        By_fft.define(amrex::convert(ba_fft,By_nodal_flag), dm_fft, 1, 0);
        Bz_fft.define(amrex::convert(ba_fft,Bz_nodal_flag), dm_fft, 1, 0);
        jx_fft.define(amrex::convert(ba_fft,jx_nodal_flag), dm_fft, 1, 0);
        jy_fft.define(amrex::convert(ba_fft,jy_nodal_flag), dm_fft, 1, 0);
        jz_fft.define(amrex::convert(ba_fft,jz_nodal_flag), dm_fft, 1, 0);
        rho1_fft.define(amrex::convert(ba_fft,IntVect::TheNodeVector()), dm_fft, 1, 0);
        rho2_fft.define(amrex::convert(ba_fft,IntVect::TheNodeVector()), dm_fft, 1, 0);
        Ex_fft.setVal(0.0);
        Ey_fft.setVal(0.0);
        Ez_fft.setVal(0.0);
        Bx_fft.setVal(0.0);
        By_fft.setVal(0.0);
        Bz_fft.setVal(0.0);
        jx_fft.setVal(0.0);
        jy_fft.setVal(0.0);
        jz_fft.setVal(0.0);
        rho1_fft.setVal(0.0);
        rho2_fft.setVal(0.0);
    }

    // initialize FFTW plans
    {
        for (MFIter mfi(Ex_fft); mfi.isValid(); ++mfi)
        {
            AMREX_ALWAYS_ASSERT(Ex_fft.local_size() == 1);

            const Box& local_domain = amrex::enclosedCells(mfi.fabbox());
            tps_fft_init(domain_fft.loVect(), domain_fft.hiVect(),
                         local_domain.loVect(), local_domain.hiVect());
        }
    }

    // loop over iterations
    for (int istep = 0; istep < max_step; ++istep)
    {
        amrex::Print() << "Step " << istep << " / " << max_step << "\n";

        // copy data in
        {
            Ex_fft.ParallelCopy(Ex, 0, 0, 1, 0, 0, geom.periodicity());
            Ey_fft.ParallelCopy(Ey, 0, 0, 1, 0, 0, geom.periodicity());
            Ez_fft.ParallelCopy(Ez, 0, 0, 1, 0, 0, geom.periodicity());
            Bx_fft.ParallelCopy(Bx, 0, 0, 1, 0, 0, geom.periodicity());
            By_fft.ParallelCopy(By, 0, 0, 1, 0, 0, geom.periodicity());
            Bz_fft.ParallelCopy(Bz, 0, 0, 1, 0, 0, geom.periodicity());
            jx_fft.ParallelCopy(jx, 0, 0, 1, 0, 0, geom.periodicity());
            jy_fft.ParallelCopy(jy, 0, 0, 1, 0, 0, geom.periodicity());
            jz_fft.ParallelCopy(jz, 0, 0, 1, 0, 0, geom.periodicity());
        }

        // call spectral solver
        {
            for ( MFIter mfi(Ex_fft); mfi.isValid(); ++mfi ) {
                
                tps_push_eb(
                    BL_TO_FORTRAN_ANYD(Ex_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(Ey_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(Ez_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(Bx_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(By_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(Bz_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(jx_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(jy_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(jz_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(rho1_fft[mfi]),
                    BL_TO_FORTRAN_ANYD(rho2_fft[mfi]));
            }
        }

        // copy data out
        {
            copy_data_from_fft_to_valid(Ex, Ex_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(Ey, Ey_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(Ez, Ez_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(Bx, Bx_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(By, By_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(Bz, Bz_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(jx, jx_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(jy, jy_fft, geom, ba_valid_fft);
            copy_data_from_fft_to_valid(jz, jz_fft, geom, ba_valid_fft);
        }

        if (plot_int > 0 && (istep+1) % plot_int == 0) {
            write_plotfile(Ex,Ey,Ez,Bx,By,Bz,geom,istep+1);
        }
    }



    // free MPI
    {
        tps_mpi_finalize();
        MPI_Comm_free(&comm_fft);
    }
}


int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    toy();

    amrex::Finalize();
}


void write_plotfile (const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                     const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                     const Geometry& geom, int istep)
{
    const BoxArray& ba = amrex::convert(Ex.boxArray(), IntVect::TheZeroVector());
    MultiFab plotmf(ba, Ex.DistributionMap(), 6, 0);

    amrex::average_edge_to_cellcenter(plotmf, 0, {&Ex, &Ey, &Ez});
    amrex::average_face_to_cellcenter(plotmf, 3, {&Bx, &By, &Bz});

    amrex::WriteSingleLevelPlotfile(amrex::Concatenate("./data/plt",istep),
                                    plotmf,
                                    {"Ex", "Ey", "Ez", "Bx", "By", "Bz"},
                                    geom, 0., 0);
}

void copy_data_from_fft_to_valid (MultiFab& mf, const MultiFab& mf_fft,
                                  const Geometry& /* geom */, const BoxArray& ba_valid_fft)
{
    auto idx_type = mf_fft.ixType();
    MultiFab mftmp(amrex::convert(ba_valid_fft,idx_type), mf_fft.DistributionMap(), 1, 0);
    for (MFIter mfi(mftmp,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        if (mf_fft[mfi].box().contains(bx))
        {
            mftmp[mfi].copy(mf_fft[mfi], bx, 0, bx, 0, 1);
        }
    }

    mf.ParallelCopy(mftmp);
}
