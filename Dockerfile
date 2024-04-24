FROM ubuntu:latest

LABEL org.opencontainers.image.description "Nystrom Boundary Integral library"
LABEL org.opencontainers.image.source https://github.com/mjcarley/nbi
LABEL org.opencontainers.image.licenses GPL

## all these packages are required for compiling, etc, except for
##
## libblas-dev, liblapack-dev: BLAS and LAPACK
## git: downloading packages
## bc: used in an NBI example to set wavenumber from rotor data

RUN apt-get update && apt-get install -y build-essential libgtk2.0-dev \
    pkg-config libtool autoconf gcc gfortran libblas-dev liblapack-dev \
    git libgmsh-dev cmake bc

## PETSc needs to be installed without MPI, so Ubuntu packages are not
## suitable
RUN git clone -b release https://gitlab.com/petsc/petsc.git petsc
WORKDIR petsc
RUN git pull
RUN git checkout v3.18.5
## install to /usr to make it easy to find later
RUN ./configure --with-mpi=0 --with-scalar-type=real --with-threadsafety \
    --with-debugging=0 --with-log=0 --with-openmp --prefix=/usr
RUN make
RUN make install

## wrappers for BLAS and LAPACK libraries, with a config file for
## setting flags
RUN git clone https://github.com/mjcarley/blaswrap.git
WORKDIR blaswrap
RUN bash autogen.sh
RUN ./configure
RUN make 
RUN make install
WORKDIR /

## this is required by AGG
RUN git clone https://github.com/wo80/Triangle/
WORKDIR Triangle/src
RUN mkdir build
WORKDIR build
RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DCMAKE_BUILD_TYPE=Release ..
RUN make install

WORKDIR /

## main set of libraries, including submodules distributed with NBI
RUN git clone --recursive https://github.com/mjcarley/nbi
WORKDIR nbi
RUN bash autogen.sh
RUN ./configure --prefix=/usr/
RUN make 
RUN make install
WORKDIR /