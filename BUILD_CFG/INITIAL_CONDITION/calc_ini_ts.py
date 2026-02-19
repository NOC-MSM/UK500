import xarray as xr
from scipy.interpolate import griddata
from scipy.interpolate import NearestNDInterpolator
import numpy as np
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

def dep_interpolate_lev(ds, cfg, lev, var):
    """ nearest neighbour interpolation for a single level """

    src_lon =  ds.nav_lon.values.flatten()
    src_lat =  ds.nav_lat.values.flatten()
    src_dep =  ds.deptht.values.flatten()
    
    points = (src_lon, src_lat)
    values = (ds[var].isel(deptht=lev).values.flatten())
    
    tgt_lon =  cfg.nav_lon
    tgt_lat =  cfg.nav_lat
    
    target = (tgt_lon, tgt_lat)
    
    n_grid = griddata(points, values, target, method='nearest')

    return n_grid

def dep_nd_interpolate_lev(da, cfg, lev=None):
    """ nearest neighbour interpolation for a single level """

    if lev: # reduce depth dimension
        values = (da.isel(deptht=lev).values.flatten())
    else: # already 2d
        values = (da.values.flatten())

    # flatten source coordinates
    src_lon =  da.nav_lon.values.flatten()
    src_lat =  da.nav_lat.values.flatten()
    
    # drop missing source data
    src_lon = src_lon[~np.isnan(values)]
    src_lat = src_lat[~np.isnan(values)]
    values = values[~np.isnan(values)]

    # format source sdata
    points = list(zip(src_lat, src_lon))

    # target coordinates 
    with ProgressBar():
        tgt_lon =  cfg.nav_lon.load()
        tgt_lat =  cfg.nav_lat.load()
  
    # test if layer is full of nan values
    if values.size > 0:
        # interpolate
        interp = NearestNDInterpolator(points, values)
        n_grid = interp(tgt_lat, tgt_lon)
        return n_grid

def dep_3d_interpolate_lev(da, cfg, save=False, load=False, chunks=-1):

    # set src coords as 3d
    src_lon_3d, src_dep_3d = xr.broadcast(da.nav_lon, da.deptht)
    src_lat_3d, src_dep_3d = xr.broadcast(da.nav_lat, da.deptht)

    # set dims
    dims = ["deptht","y","x"]

    # flatten input
    values = da.stack(z=dims)

    # set dst coords as 3d
    tgt_dep, tgt_lon = xr.broadcast(cfg.gdept_0, cfg.nav_lon)
    tgt_dep, tgt_lat = xr.broadcast(cfg.gdept_0, cfg.nav_lat)

    index = ~np.isnan(values).compute()
    values = values[index]
    if load:
        dst_path = '/gws/nopw/j04/verify_oce/NEMO/Preprocessing/'
        points = xr.open_dataarray(dst_path + "pre_interp_ds_tmp_points.nc",
                                   chunks=chunks)
    else:
        src_lon = src_lon_3d.stack(z=dims)
        src_lat = src_lat_3d.stack(z=dims)
        src_dep = src_dep_3d.stack(z=dims)

        src_lon = src_lon[index].expand_dims(var_dim=["x"])
        src_lat = src_lat[index].expand_dims(var_dim=["y"])
        src_dep = src_dep[index].expand_dims(var_dim=["z"])

        points = xr.concat([src_dep,src_lat,src_lon],dim="var_dim")
        points = points.reset_index("z")

    # intermediate save of  inputs
    if save:
        with ProgressBar():
            dst_path = '/gws/nopw/j04/verify_oce/NEMO/Preprocessing/'
            points.to_netcdf(dst_path + "pre_interp_ds_tmp_points.nc")
        print ("SAVED")

    # load to memory
    with ProgressBar():
        points = points.load()
        values = values.load()

    # interpolate
    interp = NearestNDInterpolator(points.T, values)
    n_grid = xr.apply_ufunc(interp, tgt_dep, tgt_lat, tgt_lon,
                            dask="parallelized")

    return n_grid

def check_depths(cfg):
    '''
    check for gdept_0 and, if absent, create
    '''

    depth_dim = set(["nav_lev","z"]).intersection(set(cfg.dims))
    if 'gdept_0' not in list(cfg.keys()):
        print (cfg)
        cfg['gdept_0'] = cfg.e3t_0.cumsum(dim=depth_dim)
    
        with ProgressBar():
            cfg.to_netcdf("NAARC_cfg.nc")

    return cfg

def interp_var(var, src_fn, cfg_fn, dst_fn):
    '''
    Interpolate source gridded model output to target nemo grid 

    This code uses nearest neighbour interpolation to fill bathymetric
    discrepancies.
    '''
    
    # This is to ignore variable incompatible with chunking such as object dtype
    # WARNING: not very flexible and my break code
    drop_vars = list(xr.open_dataset(src_fn, chunks=-1).keys())
    drop_vars.remove(var)
    drop_vars = drop_vars + ["time_centered"]

    # open datasets
    chunks = dict(x=100,y=100)
    da = xr.open_dataset(src_fn, chunks=chunks, drop_variables=drop_vars
                        )[var].squeeze()

    # interpolate
    tgt_cfg = xr.open_dataset(cfg_fn, chunks=chunks).squeeze()
    tgt_cfg = check_depths(tgt_cfg)
    da_n = dep_3d_interpolate_lev(da, tgt_cfg)

    # set coordinates
    #coords=dict(deptht=(['y','x','deptht'], tgt_cfg.gdept_0.data),
    #            nav_lat=(['y','x'], tgt_cfg.nav_lat.data),
    #            nav_lon=(['y','x'], tgt_cfg.nav_lon.data)
    #           )

    da_n = da_n.assign_coords(
           dict(nav_lat=(['y','x'], tgt_cfg.nav_lat.data),
                nav_lon=(['y','x'], tgt_cfg.nav_lon.data)))

    da_n.name = var

    # save
    with ProgressBar():
        da_n.to_netcdf(dst_fn)

def interp_surface(var_dict, src_fn, cfg_fn, dst_fn):
    '''
    Interpolate list of surface vars to target nemo grid
    '''

    # This is to ignore variable incompatible with chunking such as object dtype
    # WARNING: not very flexible and my break code
    full_var_list = list(xr.open_dataset(src_fn, chunks=-1).variables.keys())
    var_dict['TLON'] = 'nav_lon'
    var_dict['TLAT'] = 'nav_lat'
    drop_vars = [var for var in full_var_list if var not in list(var_dict.keys())]

    # get source data
    chunks = -1
    ds = xr.open_dataset(src_fn, chunks=chunks, drop_variables=drop_vars)
    ds = ds.rename({'TLON':'nav_lon', 'TLAT':'nav_lat'})

    # shift lons
    ds['nav_lon'] = xr.where(ds.nav_lon > 180, ds.nav_lon  - 360, ds.nav_lon)
    del var_dict['TLON']
    del var_dict['TLAT']

    # interpolate
    tgt_cfg = xr.open_dataset(cfg_fn, chunks=chunks, decode_times=False).squeeze()
    interp_list = [] 
    for var in var_dict.keys():
        if var in ['TLON','TLAT','time']:
            continue
        with ProgressBar():
            da = ds[var].squeeze().load()
        da_n = dep_nd_interpolate_lev(da, tgt_cfg)
        da_n_xr = xr.DataArray(name=var, data=da_n, dims=('y','x'))
        interp_list.append(da_n_xr)

    ds_n = xr.merge(interp_list)
    ds_n = ds_n.expand_dims(time=ds.time.data)
    ds_n = ds_n.assign_coords(
           dict(nav_lat=(['y','x'], tgt_cfg.nav_lat.data),
                nav_lon=(['y','x'], tgt_cfg.nav_lon.data)))

    ds_n = ds_n.rename(var_dict)

    with ProgressBar():
        ds_n.to_netcdf(dst_fn)


def interpolate_glosea6_to_co9(var, y='1993', m='01', d='01',
                               domcfg='GEG_SF12.nc'):
    ''' interpolate glosea6 data to co9 target grid '''

    # set file paths
    cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/' + domcfg
    src_path = '/gws/nopw/j04/jmmp/MASS/GloSea6/Daily/'
    src_fn = src_path + 'glosea6_grid_T_' + y + m + d + '.nc'
    dst_fn = 'glosea_ini_' + y + m + d + '_' + \
              domcfg.replace('.nc','') + '_'  + var + '.nc'

    # interpolate
    interp_var(var, src_fn, cfg_fn, dst_fn)

def flood_gosi8(var, y='1850', m='01', d='01'):
    ''' take UKESM historical glosat and flood to gosi8 grid '''

    # set file paths
    dst_path = '/gws/nopw/j04/verify_oce/NEMO/Preprocessing/'
    cfg_fn = dst_path + '/DOM/NAARC/domain_cfg_mes.nc'
    src_path = '/gws/nopw/j04/glosat/production/UKESM/raw/u-ck651/18500101T0000Z/'
    src_fn = src_path + 'nemo_ck651o_1m_18500101-18500201_grid-T.nc'
    dst_fn = dst_path + 'INI/glosat_ukesm_to_gosi8_mes_' + var + '.nc'

    # interpolate
    interp_var(var, src_fn, cfg_fn, dst_fn)

def create_gosi8_sea_ice_instance(y, m):
    ''' extract sea ice data from UKESM historical glosat data '''

    # set file paths

    # variable list
    var_list = {'time':'time_counter', # ice thickness
                'hi':'hti', # ice thickness
                'hs':'hts', # snow thickness
                'aice':'ati'} # ice concentration
                #'sice':'smi', # ice salinity
                #'Tinz':'tmi', # ice temperature
                #'Tsfc':'tsu', # surface temperature
                #'Tsnz':'tms'} # snow temperature

    yyyy_0 = str(y)
    mm = str(m).zfill(2)
    if m == 12:
        y = y-1
    date0 = np.datetime64(str(y) + '-' + str(m).zfill(2))
    date1 = date0 + np.timedelta64(1, 'M')
    mm_0 = date0.item().strftime('%m')
    mm_1 = date1.item().strftime('%m')
    yyyy_0 = date0.item().strftime('%Y')
    yyyy_1 = date1.item().strftime('%Y')

    if m in [1,2,12]:
        dir_mm = "01"
    elif m in [3,4,5]:
        dir_mm = "04"
    elif m in [6,7,8]:
        dir_mm = "07"
    elif m in [9,10,11]:
        dir_mm = "10"

    glosat_path = '/gws/nopw/j04/glosat/production/UKESM/raw/'
    src_path = glosat_path + f'u-ck651/{yyyy_1}{dir_mm}01T0000Z/'
    src_fn = src_path + f'cice_ck651i_1m_{yyyy_0}{mm_0}01-{yyyy_1}{mm_1}01.nc'

    dst_path = '/gws/nopw/j04/verify_oce/NEMO/Preprocessing/LBC/ICE_src/'
    dst_fn = dst_path + f'glosat_ukesm_to_gosi8_sea_ice_y{yyyy_0}m{mm_0}.nc'
    tgt_path = '/gws/nopw/j04/verify_oce/NEMO/Preprocessing/DOM/NAARC/'
    tgt_cfg = tgt_path + 'domain_cfg_zps_with_gdept_and_nav.nc'

    interp_surface(var_list, src_fn, tgt_cfg, dst_fn)

def create_uniform_forcing():
    cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
    cfg = xr.open_dataset(cfg_fn)

    uniform_t = xr.full_like(cfg.tmask, 5).squeeze()
    uniform_t.name = 'tn'
    uniform_s = xr.full_like(cfg.tmask, 35).squeeze()
    uniform_s.name = 'sn'

    uniform_t.to_netcdf("amm15_uniform_t.nc")
    uniform_s.to_netcdf("amm15_uniform_s.nc")

def create_uniform_forcing_masked():
    cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
    cfg = xr.open_dataset(cfg_fn)

    uniform_t = xr.where(cfg.tmask == 1, 5, np.nan).squeeze()
    uniform_t.name = 'tn'
    uniform_s = xr.where(cfg.tmask == 1, 35, np.nan).squeeze()
    uniform_s.name = 'sn'

    uniform_t.to_netcdf("amm15_uniform_t_masked.nc")
    uniform_s.to_netcdf("amm15_uniform_s_masked.nc")
 
if __name__ == '__main__':
    #interpolate_glosea6_to_co9('vosaline', domcfg='CO7_EXACT_CFG_FILE.nc')
    #flood_gosi8('thetao', y='1850', m='01', d='01')
    flood_gosi8('so', y='1850', m='01', d='01')
    #create_gosi8_sea_ice_instance(1860, 12)
    #for year in range (1852,1860):
    #    for month in range(1,13):
    #        print (year, month)
    #        create_gosi8_sea_ice_instance(year, month)
    #interpolate_glosea6_to_co9('votemper', domcfg='CO7_EXACT_CFG_FILE.nc')
