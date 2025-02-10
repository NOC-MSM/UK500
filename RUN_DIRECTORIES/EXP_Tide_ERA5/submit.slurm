#!/bin/bash
#SBATCH --job-name=nemo_test
#SBATCH --time=23:59:00
#SBATCH --nodes=13
#SBATCH --cpus-per-task=1
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the job environment (this module needs to be loaded before any other modules)
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1
export UCX_IB_REG_METHODS=direct
ulimit -s unlimited
ulimit -c unlimited

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

# Launch the parallel job
#   Using 512 MPI processes and 128 MPI processes per node
#   srun picks up the distribution from the sbatch options
srun --het-group=0 --nodes=1 --tasks-per-node=2 ${SHARED_ARGS} /work/n01/n01/$USER/UK500/BUILD_EXE/XIOS/xios-2.5/bin/xios_server.exe : \
     --het-group=1 --nodes=12 --tasks-per-node=126 ${SHARED_ARGS} /work/n01/n01/$USER/UK500/BUILD_EXE/NEMO/4.0.6/cfgs/NEMOconstTS/BLD/bin/nemo.exe

