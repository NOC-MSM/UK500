import sys
sys.path.append("/work/n01/n01/jrule/CHAMFER/NEW/UK500/BUILD_CFG/INITIAL_CONDITION/")

from calc_ini_ts import dep_3d_interpolate_lev, check_depths
import xarray as xr
import numpy as np

# paths
variables = ["votemper","vosaline"]
src_fn = "/work/n01/n01/jrule/CHAMFER/NEW/UK500/BUILD_CFG/INITIAL_CONDITION/Input/CHAMFER_AMM15_01101993_chunked400.nc"
cfg_fn = "/work/n01/n01/jrule/CHAMFER/NEW/UK500/BUILD_CFG/INITIAL_CONDITION/Input/domain_cfg.0m.nc"
#mask = "/work/n01/n01/jrule/CHAMFER/NEW/UK500/BUILD_CFG/INITIAL_CONDITION/Input/mesh_mask.nc"
dst_fn = "/work/n01/n01/jrule/CHAMFER/NEW/UK500/BUILD_CFG/INITIAL_CONDITION/Output/UK500_ICS_01101993.nc"

# tmask has dimensions (t, z, y, x)
# Extract tmask for first time step (assuming your ICS variables will also have time_counter=1)
# open mask once
# mask_ds = xr.open_dataset(mask)
# tmask = mask_ds["tmask"].isel(time_counter=0).squeeze(drop=True)

# chunk sizes (tune these)
y_chunk = 400
x_chunk =400

# set first equal to true so you can save a netcdf on first variable and append the second variable to that
first = True

#Start operations on each variable
for var in variables:
    # open source ( and drop singleton dimentions - you should have a 3d x/y/z matrix)
    da = xr.open_dataset(src_fn)[var] \
           .isel(time_counter=0) \
           .squeeze(drop=True)

    # open target (load to memory as chunks, so it's not too heavy)
    cfg = xr.open_dataset(
        cfg_fn,
        chunks={"y": y_chunk, "x": x_chunk}
    )

    ny = cfg.dims["y"]
    nx = cfg.dims["x"]

    # allocate output 
    # for this operation you just need the size of gdept_0 to pregenerate the matrix to fill. e3t_0 is the same size if gdept doesn't exist. 
    # NOTE: if you try to apply check depth before chuncking with UK500 it's quite heavy on the memory and makes the job too slow or crash, 
    # so this generates a full matrix of the size of the cfg_fn variables, then operates and fills is chunck by chunck. 
    if "gdept_0" in cfg:
        template = cfg.gdept_0
    else:
        template = cfg.e3t_0

    out = xr.full_like(template, np.nan)  # allocate full size once

    # compute check_depth and 3D interpolation
    for y0 in range(0, ny, y_chunk):
        y1 = min(y0 + y_chunk, ny)
        for x0 in range(0, nx, x_chunk):
            x1 = min(x0 + x_chunk, nx)

            print(f"Interpolating y:{y0}-{y1}, x:{x0}-{x1}", flush=True)

            cfg_sub = check_depths(cfg.isel(y=slice(y0, y1), x=slice(x0, x1)).squeeze())
            da_sub = dep_3d_interpolate_lev(da, cfg_sub)

            out.loc[dict(y=slice(y0, y1), x=slice(x0, x1))] = da_sub


    out.name = var
#    out = out.where(tmask != 0, 0.0) # apply mask 
    if first:
        out.to_netcdf(dst_fn)
        first = False 
        print("Generated NetCDF with first variable")
    else:
        out.to_netcdf(dst_fn, mode="a")
        print("Added second variable to NetCDF")

print("end of script")

