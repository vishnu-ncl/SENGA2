#!/bin/bash

# maths_kernel eqAV
sed -i 's/#include \"maths_kernel_eqAV_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqAV_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqAV_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqAV_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqAV_f2c_kernel.F90


# maths_kernel_eqBC
sed -i 's/#include \"maths_kernel_eqBC_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqBC_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqBC_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqBC_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqBC_f2c_kernel.F90


# maths_kerel_eqAW
sed -i 's/#include \"maths_kernel_eqAW_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqAW_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqAW_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqAW_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqAW_f2c_kernel.F90


# maths_kerel_eqAX
sed -i 's/#include \"maths_kernel_eqAX_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqAX_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqAX_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqAX_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqAX_f2c_kernel.F90


# maths_kerel_eqB
sed -i 's/#include \"maths_kernel_eqB_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqB_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqB_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqB_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqB_f2c_kernel.F90


# maths_kerel_eqBN
sed -i 's/#include \"maths_kernel_eqBN_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqBN_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqBN_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqBN_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqBN_f2c_kernel.F90

# maths_kerel_eqBFG
sed -i 's/#include \"maths_kernel_eqBFG_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqBFG_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqBFG_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqBFG_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqBFG_f2c_kernel.F90

# maths_kerel_eqBIJK
sed -i 's/#include \"maths_kernel_eqBIJK_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqBIJK_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqBIJK_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqBIJK_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqBIJK_f2c_kernel.F90

# maths_kerel_eqA
sed -i 's/#include \"maths_kernel_eqA_f2c_kernel.cpp\"/\/\/#include \"maths_kernel_eqA_f2c_kernel.cpp\"/g' f2c_mpi_openmp/f2c_mpi_openmp_kernels.cpp

rm f2c_mpi_openmp/maths_kernel_eqA_f2c_kernel.cpp
cp mpi_openmp/maths_kernel_eqA_seq_kernel.F90 f2c_mpi_openmp/maths_kernel_eqA_f2c_kernel.F90
