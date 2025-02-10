#!/bin/bash
#SBATCH --job-name=weight
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --account=n01-CHAMFER
#SBATCH --partition=serial
#SBATCH --qos=serial


#load modules
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1
export OMP_NUM_THREADS=1


#generate weights
/work/n01/n01/micdom/UK500/BUILD_EXE/NEMO/4.0.6/tools/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
/work/n01/n01/micdom/UK500/BUILD_EXE/NEMO/4.0.6/tools/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
/work/n01/n01/micdom/UK500/BUILD_EXE/NEMO/4.0.6/tools/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos
/work/n01/n01/micdom/UK500/BUILD_EXE/NEMO/4.0.6/tools/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
/work/n01/n01/micdom/UK500/BUILD_EXE/NEMO/4.0.6/tools/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos


