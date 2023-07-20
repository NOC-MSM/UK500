#!/bin/bash
#SBATCH --job-name=pynemo
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=serial
#SBATCH --qos=serial

source ~/.bashrc
conda activate pynemo
cd /work/n01/n01/$USER/UK500/BUILD_CFG/OPEN_BOUNDARIES/
pynemo -s namelist_FES14.bdy

