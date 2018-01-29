AMREX_HOME ?= ../amrex
PICSAR_HOME ?= ../picsar

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE
USE_OMP   = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(PICSAR_HOME)/src/Make.package


F90EXE_sources += mpi_derived_types.F90

F90EXE_sources += maxwell_solver_manager.F90
VPATH_LOCATIONS   += $(PICSAR_HOME)/src/field_solvers/Maxwell/
INCLUDE_LOCATIONS += $(PICSAR_HOME)/src/field_solvers/Maxwell/

F90EXE_sources += fastfft.F90 fourier_psaotd.F90 GPSTD.F90 init_kspace_3D.F90
VPATH_LOCATIONS   += $(PICSAR_HOME)/src/field_solvers/Maxwell/GPSTD_solver
INCLUDE_LOCATIONS += $(PICSAR_HOME)/src/field_solvers/Maxwell/GPSTD_solver

# yee_solver
F90EXE_sources += yee.F90

F90EXE_sources += load_balancing.F90
VPATH_LOCATIONS   += $(PICSAR_HOME)/src/housekeeping
INCLUDE_LOCATIONS += $(PICSAR_HOME)/src/housekeeping

F90EXE_sources += field_boundaries.F90
VPATH_LOCATIONS   += $(PICSAR_HOME)/src/boundary_conditions
INCLUDE_LOCATIONS += $(PICSAR_HOME)/src/boundary_conditions


VPATH_LOCATIONS += $(FFTW_HOME)/include
INCLUDE_LOCATIONS += $(FFTW_HOME)/include
LIBRARY_LOCATIONS += $(FFTW_HOME)/lib


DEFINES += -DFFTW
DEFINES += -DPICSAR_NO_ASSUMED_ALIGNMENT


LIBRARIES += -lfftw3_mpi -lfftw3 -lfftw3_omp


include $(AMREX_HOME)/Tools/GNUMake/Make.rules
