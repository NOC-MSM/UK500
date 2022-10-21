
#!/bin/bash

#:'
#
#******************************
#run_Tide_ERA5_River_WAD.sh
#******************************
#'

# Run fully forced experiments with wetting & drying
# before running copy the correct namelist to activate 
# initial conditions, river forcing, BDY
#::

export CONFIG=NEMOconstTS
export EXP=$WDIR/RUN_DIRECTORIES/EXP_Tide_ERA5_River_IC_BDY_WAD

# Choose an appropriate directory for your EXP installation
if [ ! -d "$EXP/RESTART" ]; then
  #mkdir $EXP
  mkdir $EXP/RESTART
fi

rsync -av --ignore-existing $NEMO/cfgs/SHARED/*namelist* $EXP/. # only get the files not already in the repo.
rsync -av --ignore-existing $NEMO/cfgs/SHARED/*.xml $EXP/. 

# Copy in NEMO/XIOS executables
ln -s $NEMO/cfgs/$CONFIG/BLD/bin/nemo.exe $EXP/nemo.exe
ln -s $XIOS_DIR/bin/xios_server.exe $EXP/xios_server.exe

# Link in domain_cfg file
rm $EXP/domain_cfg.nc
ln -s $DOMAIN/domain_cfg_$REPO.nc $EXP/domain_cfg.nc

# Link in tidal bondary forcing
#ln -s /work/n01/n01/annkat/SEAsia_HadGEM_R12/TIDES $EXP/.
ln -s $WDIR/INPUTS/TIDES $EXP/.

# Link in boundary files (just coordinates.bdy.nc)
ln -s $WDIR/INPUTS/OBC/coordinates.bdy.nc $EXP/.
ln -s $WDIR/INPUTS/OBC/ $EXP/.

# Submit job
cd $EXP
sbatch submit.slurm

## Check on queue
# squeue -u $USER
