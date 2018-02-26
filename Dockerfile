FROM ubuntu:14.04

# Install a few packages, as root
RUN apt-get update \
    && apt-get install -y \
    wget \
    make \
    git \
    gcc \
    python

# Install MPI
RUN apt-get install -y \
    libcr-dev \
    mpich2
    
# Install FFTW
RUN apt-get install -y \
    libfftw3-dev \
    libfftw3-mpi-dev

# Download AMReX, PICSAR
RUN git clone https://github.com/AMReX-Codes/amrex.git \
    && cd amrex \
    && git checkout development
RUN git clone https://bitbucket.org/berkeleylab/picsar.git

# Copy the toy pseudo-spectral
RUN mkdir tps/
COPY GNUmakefile tps/
COPY Make.package tps/
COPY TPS_F.H tps/
COPY main.cpp tps/
COPY tps.F90 tps/

# Build the software
RUN cd tps \
    && export FFTW_HOME=/usr/ \
    && make

# Run the code for testing
RUN cd tps \
    && mpirun -np 4 ./main3d.gnu.DEBUG.MPI.ex
