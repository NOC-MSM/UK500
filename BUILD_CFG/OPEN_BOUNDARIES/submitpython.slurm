#!/bin/bash
#SBATCH --job-name=python
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=serial
#SBATCH --qos=serial

source ~/.bashrc
conda activate coast
cd /work/n01/n01/$USER/CHAMFER/NEW/UK500/BUILD_CFG/OPEN_BOUNDARIES/

python ./adjust_parent_domcfg.py
~                                  
