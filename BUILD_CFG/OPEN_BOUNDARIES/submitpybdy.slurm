#!/bin/bash
#SBATCH --job-name=pybdy
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --mem=64G

source ~/.bashrc
export JAVA_HOME=/usr/lib64/jvm/java

conda activate pybdy
cd /work/n01/n01/$USER/CHAMFER/NEW/UK500/BUILD_CFG/OPEN_BOUNDARIES/

pybdy -s namelist_FES14_pybdy.bdy

