FROM ubuntu:14.04

# Install a few packages, as root
RUN apt-get update \
    && apt-get install -y \
    wget \
    make \
    git \
    gcc

# Install MPI
RUN apt-get install -y \
    libcr-dev \
    mpich2
    
# Install FFTW
RUN apt-get install -y \
    libfftw3-dev \
    libfftw3-mpi-dev

