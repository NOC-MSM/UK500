import xarray as xr
from dask.diagnostics import ProgressBar

def add_gdep_var(path, fn_in, fn_out, z_var="z"):
    cfg = xr.open_dataset(path + fn_in, chunks="auto")

    cfg["gdept"] = cfg.e3t_0.cumsum(z_var)
    cfg["gdepw"] = cfg.e3w_0.cumsum(z_var)

    cfg = cfg.rename({"e3u_0":"e3u",
                      "e3t_0":"e3t",
                      "e3v_0":"e3v",
                      "e3w_0":"e3w",
                      "bottom_level":"mbathy"})

    with ProgressBar():
        cfg.to_netcdf(path + fn_out)

def rename_bdy_msk_var(path, fn_in, fn_out):
    msk = xr.open_dataarray(path + fn_in, chunks="auto")

    msk.name = "mask"



    with ProgressBar():
        msk.to_netcdf(path + fn_out)

def get_Bathymetry(path, fn_in, fn_out):
    cfg = xr.open_dataset(path + fn_in, chunks="auto")

    bathy = cfg.bathy_metry

    bathy.name = "Bathymetry"
    nav_lon = bathy.x.drop_vars(["x","y"])
    nav_lat = bathy.y.drop_vars(["x","y"])
    bathy = bathy.drop_vars(["x","y"])
    bathy = bathy.assign_coords({"nav_lon": nav_lon,
                                 "nav_lat": nav_lat})

    with ProgressBar():
        bathy.to_netcdf(path + fn_out)

def add_gdep_CMEMS(path, fn_in, fn_out, z_var="z"):
    cfg = xr.open_dataset(path + fn_in, chunks="auto")

    cfg = cfg.rename({"depth":"gdept_0",
                      "deptho":"mbathy"})

    with ProgressBar():
        cfg.to_netcdf(path + fn_out)

path = "./"
fn_in = "domain_cfg.nc"
fn_out = "domain_cfg_gdept.nc"
add_gdep_var(path, fn_in, fn_out, z_var="z")

#fn_in = "CMEMS_subdomain_coordinates.nc"
#fn_out = "CMEMS_subdomain_coordinates_gdept.nc"
#fn_in = "NWS-MFC_004_013_mask_bathy.nc"
#fn_out = "NWS-MFC_004_013_mask_bathy_pybdy.nc"
#add_gdep_CMEMS(path, fn_in, fn_out, z_var="z")
