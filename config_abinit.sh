#!/bin/bash

# Activate your virtual environment first
# (You must do it outside this script.)

# Set compilers
export CC=mpicc
export FC=mpifort
export CXX=mpicxx

# Enable full debug flags for gfortran
export CFLAGS="-g -O2"
# export FCFLAGS="-g -O0 -Wno-error=line-truncation -std=legacy"
export FCFLAGS="-g -O2 -fallow-argument-mismatch -fno-range-check \
                -ffree-line-length-none -std=legacy -Wno-error=line-truncation"
export CXXFLAGS="-g -O2"

# Use virtualenv path
VENV_PREFIX="/home/yxcheng/software/miniforge3/envs/abinit"

# NetCDF flags
export NETCDF_CPPFLAGS="-I${VENV_PREFIX}/include"
export NETCDF_LIBS="-L${VENV_PREFIX}/lib -lnetcdf"
# NetCDF Fortran flags
export NETCDF_FORTRAN_CPPFLAGS="-I${VENV_PREFIX}/include"
export NETCDF_FORTRAN_LIBS="-L${VENV_PREFIX}/lib -lnetcdff"

# HDF5 flags
export HDF5_CPPFLAGS="-I${VENV_PREFIX}/include"
export HDF5_LIBS="-L${VENV_PREFIX}/lib -lhdf5_fortran -lhdf5"

# LibXC flags (if available)
export LIBXC_CPPFLAGS="-I${VENV_PREFIX}/include"
export LIBXC_LIBS="-L${VENV_PREFIX}/lib -lxcf90 -lxc"

# Install path
INSTALL_PREFIX="/home/yxcheng/software/bin/abinit"

# Finally configure
../configure \
  --prefix="${INSTALL_PREFIX}" \
  --with-mpi \
  --enable-mpi-io \
  --enable-openmp \
  --enable-source-debug
