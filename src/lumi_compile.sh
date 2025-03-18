#!/bin/bash

ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c ops_data_init_ops.F90 -o ops_data_init_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c bounds_ops.F90 -o bounds_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c boundt_ops.F90 -o boundt_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c bountt_ops.F90 -o bountt_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c temper_ops.F90 -o temper_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c chrate_ops.F90 -o chrate_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c indata_ops.F90 -o indata_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c tempin_ops.F90 -o tempin_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c d2fdxy_ops.F90 -o d2fdxy_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c d2fdxz_ops.F90 -o d2fdxz_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c d2fdyz_ops.F90 -o d2fdyz_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c rhscal_ops.F90 -o rhscal_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c rhsvel_ops.F90 -o rhsvel_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c lincom_ops.F90 -o lincom_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c fincom_ops.F90 -o fincom_ops_f2c_mpi_hip.o
ftn -O1 -g -fopenmp -M 878 -f free -N 1023 -J$OPS_INSTALL_PATH/fortran/mod/cray/f2c_hip -DOPS_MPI -DOPS_F2C_INTEROP -c flamin_ops.F90 -o flamin_ops_f2c_mpi_hip.o
