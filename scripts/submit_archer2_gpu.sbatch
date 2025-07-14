#!/bin/bash --login

#SBATCH --job-name=multi-GPU
#SBATCH --gpus=4
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --qos=gpu-exc

##SBATCH --job-name=single-GPU
##SBATCH --gpus=1
##SBATCH --qos=gpu-shd

#SBATCH --account=????????????
#SBATCH --partition=gpu

#SBATCH --time=02:00:00

#SBATCH --error=????????????/SENGA2/EO/job.%J.err
#SBATCH --output=????????/SENGA2/EO/job.%J.out

# Enable GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

echo "Submission directory is ${SLURM_SUBMIT_DIR}"

# Setup the batch environment
cd ??????????????????????
source setup_env_archer2_cray_gpu_ops.sh 

cd $SLURM_SUBMIT_DIR


srun --ntasks=4 --ntasks-per-node=4 --cpus-per-task=8 --hint=nomultithread --distribution=block:block ./?????Hip_executable???? -OPS_DIAGS=2 OPS_FORCE_DECOMP_X=2 OPS_FORCE_DECOMP_Y=2 OPS_FORCE_DECOMP_Z=1 2>&1 | tee log_512X256X256_1node_4ranks_mpi_hip.txt 


