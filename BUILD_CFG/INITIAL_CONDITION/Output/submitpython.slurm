#!/bin/bash --login

#SBATCH --job-name=python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=standard
#SBATCH --qos=short

source /work/n01/n01/jrule/myvenv/bin/activate

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

python maskICS.py > logmasking.txt


