import xarray as xr
import numpy as np

# Open the main file
ds = xr.open_dataset("UK500_ICS_01101993.nc")

# Open the mesh_mask file
mask_ds = xr.open_dataset("mesh_mask.nc")

# tmask has dimensions (t, z, y, x)
# Extract tmask for first time step (assuming your vosaline/votemper also have time_counter=1)
tmask = mask_ds["tmask"].isel(time_counter=0).data  # shape (z, y, x)

# Broadcast tmask to match dimensions of vosaline/votemper if needed
# vosaline/votemper have shape (time_counter, gdep, lat, lon)
# Here time_counter=1, gdep=z, lat=y, lon=x
mask_3d = tmask  # already matches (z,y,x)

# Set vosaline and votemper to 0 where tmask==0
# Keep original dataset attributes
ds["vosaline"].data[:, :, :, :] = np.where(mask_3d==0, 0.0, ds["vosaline"].data)
ds["votemper"].data[:, :, :, :] = np.where(mask_3d==0, 0.0, ds["votemper"].data)

# Rename dimensions
ds = ds.rename_dims({
    "t": "time_counter",
    "z": "gdep",
    "y": "lat",
    "x": "lon"
})

# Rename coordinate variables 
#ds = ds.rename({
#    "t": "time_counter",
#    "z": "gdep",
#    "y": "lat",
#    "x": "lon"
#})

# Ensure time dimension is UNLIMITED
ds.encoding["unlimited_dims"] = {"time_counter"}

# Save back to the same NetCDF file
ds.to_netcdf("UK500_ICS_y1993m10d01.nc")

