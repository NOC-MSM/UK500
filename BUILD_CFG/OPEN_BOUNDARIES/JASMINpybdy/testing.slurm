#!/bin/bash

#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --mem=128G
#SBATCH --time=03:00:00                    
#SBATCH --account=chamfer

source ~/miniforge3/bin/activate
conda activate /home/users/<user>/miniforge3/envs/pybdy
JVM="jre-17"
export JAVA_HOME=/usr/lib/jvm/$JVM/

pybdy -s namelist_local.bdy
