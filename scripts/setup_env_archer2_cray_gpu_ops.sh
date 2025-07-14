#!/bin/bash

# CRAY
export OPS_COMPILER=cray
export OPS_INSTALL_PATH=?????????????

module purge

export AMD_ARCH=MI200

module load PrgEnv-cray
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

module load cray-hdf5-parallel/1.12.2.1

export MPI_INSTALL_PATH=$CRAY_MPICH_DIR
export LD_LIBRARY_PATH=$MPI_INSTALL_PATH/lib:$LD_LIBRARY_PATH

export ROCM_PATH=/opt/rocm-5.2.3
export LD_LIBRARY_PATH=$ROCM_PATH/llvm/lib:$LD_LIBRARY_PATH
export HIP_INSTALL_PATH=$ROCM_PATH
export AOMP=$ROCM_PATH/llvm

export MPICC=cc
export MPICPP=CC
export MPICXX=CC
export MPIFC=ftn
export MPIF90=ftn

export MPICH_GPU_SUPPORT_ENABLED=1

unset HDF5_INSTALL_PATH
export HDF5_INSTALL_PATH=/opt/cray/pe/hdf5-parallel/1.12.2.1
export LD_LIBRARY_PATH=$HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH

export CRAY_CPU_TARGET=x86-64

module load cray-python/3.9.13.1
source $OPS_INSTALL_PATH/../ops_translator/ops_venv/bin/activate

module load load-epcc-module
module load arm/forge
