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
RUN git clone https://github.com/WeiqunZhang/toy_pseudo_spectral.git

# Build the software
RUN cd toy_pseudo_spectral \
    && export FFTW_HOME=/usr/ \
    && make