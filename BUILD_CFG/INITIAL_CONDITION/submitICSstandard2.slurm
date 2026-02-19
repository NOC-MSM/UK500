#!/bin/bash --login

#SBATCH --job-name=python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=standard
#SBATCH --qos=standard

source /work/n01/n01/jrule/myvenv/bin/activate

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

python UK500_icsgen.py > /work/n01/n01/jrule/CHAMFER/NEW/UK500/BUILD_CFG/INITIAL_CONDITION/log_icsgen.txt


