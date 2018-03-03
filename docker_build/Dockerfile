FROM ubuntu:14.04

# Install a few packages, as root
RUN apt-get update \
    && apt-get install -y \
    wget \
    make \
    git \
    gcc \
    gfortran \
    g++ \
    python

# Install MPI
RUN apt-get install -y \
    openmpi-bin \
    libopenmpi-dev
    
# Install FFTW
RUN apt-get install -y \
    libfftw3-dev \
    libfftw3-mpi-dev

RUN mkdir test/

# Build and the software
CMD cd test/toy_pseudo_spectral \
    && export FFTW_HOME=/usr/ \
    && make \
    && ./main3d.gnu.DEBUG.MPI.OMP.ex

