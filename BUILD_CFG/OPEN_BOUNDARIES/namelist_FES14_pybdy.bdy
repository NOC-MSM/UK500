!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OPA  : namelist for BDY generation tool
!!
!!             User inputs for generating open boundary conditions
!!             employed by the BDY module in NEMO. Boundary data
!!             can be set up for v3.2 NEMO and above.
!!
!!             More info here.....
!!
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!------------------------------------------------------------------------------
!  grid information
!------------------------------------------------------------------------------
   sn_src_hgr = './CMEMS_subdomain_coordinates_gdept.nc'
   sn_src_zgr = './CMEMS_subdomain_coordinates_gdept.nc'
   sn_dst_hgr = './domain_cfg_gdept.nc'     ! Expects vars found in domain_cfg.nc
   sn_dst_zgr = './domain_cfg_gdept.nc'     ! Expects vars: {e3u,e3v,e3w,e3t,nav_lat,nav_lon,mbathy}
   sn_src_msk = './CMEMS_subdomain_mask.nc'
   sn_bathy   = './bathy_meter.nc'    ! dst bathymetry w/o time dimension
                                                                            !Expects vars: {Bathymetry,nav_lat,nav_lon}
   sn_nme_map = '/work/n01/n01/<user>/pyBDY/inputs/grid_name_map.json'     ! json file mapping variable names to netcdf vars

!------------------------------------------------------------------------------
!  I/O
!------------------------------------------------------------------------------
   sn_src_dir = './CMEMS.ncml' ! src_files/'
   sn_dst_dir = './OUTPUTpybdy''
   sn_fn      = 'UK500_FES14'             ! prefix for output files
   nn_fv      = -1e20                 !  set fill value for output files
   nn_src_time_adj = 0                ! src time adjustment
   sn_dst_metainfo = 'CMEMS example'  ! history info

!------------------------------------------------------------------------------
!  unstructured open boundaries
!------------------------------------------------------------------------------
    ln_coords_file = .true.               !  =T : produce bdy coordinates files
    cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files
                                          !  (if ln_coords_file=.TRUE.)
    ln_mask_file   = .true.              !  =T : read mask from file
    cn_mask_file   = '../DOMAIN/bdy_mask.nc'            !  name of mask file
                                          !  (if ln_mask_file=.TRUE.)
    ln_dyn2d       = .false.               !  boundary conditions for
                                          !  barotropic fields
    ln_dyn3d       = .false.              !  boundary conditions for
                                          !  baroclinic velocities
    ln_tra         = .false.               !  boundary conditions for T and S
    ln_ice         = .false.              !  ice boundary condition
    ln_zinterp     = .false.               !  vertical interpolation
    nn_rimwidth    = 9                    !  width of the relaxation zone

!------------------------------------------------------------------------------
!  unstructured open boundaries tidal parameters
!------------------------------------------------------------------------------
    ln_tide        = .true.              !  =T : produce bdy tidal conditions
    sn_tide_model  = 'FES2014'            !  Name of tidal model. Accepts FES2014, TPXO7p2, or TPXO9v5
    clname(1)='SA'
    clname(2)='SSA'
    clname(3)='MM'
    clname(4)='MF'
    clname(5)='MTM'
    clname(6)='MSF'
    clname(7)='MSQM'
    clname(8)='K1'
    clname(9)='O1'
    clname(10)='Q1' 
    clname(11)='P1'
    clname(12)='S1'
    clname(13)='J1'
    clname(14)='M2'
    clname(15)='N2'
    clname(16)='S2'
    clname(17)='K2'
    clname(18)='L2'
    clname(19)='T2'
    clname(20)='R2'
    clname(21)='MU2'
    clname(22)='NU2'
    clname(23)='2N2'
    clname(24)='MKS2'
    clname(25)='LA2'
    clname(26)='EPS2'
    clname(27)='M3'
    clname(28)='M4'
    clname(29)='M6'
    clname(30)='M8'
    clname(31)='N4'
    clname(32)='S4'
    clname(33)='MN4'
    clname(34)='MS4'
    ln_trans       = .true.               !  interpolate transport rather than
                                          !  velocities
    ! location of TPXO7.2 data
    sn_tide_grid_7p2   = './inputs/tpxo7.2/grid_tpxo7.2.nc'
    sn_tide_h          = './inputs/tpxo7.2/h_tpxo7.2.nc'
    sn_tide_u          = './inputs/tpxo7.2/u_tpxo7.2.nc'
    ! location of TPXO9v5 data: single constituents per file
    sn_tide_grid_9p5   = './inputs/TPXO9_atlas_v5_nc/grid_tpxo9_atlas_30_v5.nc'
    sn_tide_dir        = './inputs/TPXO9_atlas_v5_nc/'
    ! location of FES2014 data
    sn_tide_fes        = '/work/n01/n01/shared/jelt/FES2014/'

!------------------------------------------------------------------------------
!  Time information for output
!------------------------------------------------------------------------------
    sn_date_start   = '1993-07-01'    !  dst output date start YYYY-MM-DD
    sn_date_end     = '1993-07-31'    !  dst output date end YYYY-MM-DD
    sn_dst_calendar = 'gregorian'     !  output calendar format
    sn_date_origin  = '1970-01-01'    !  reference for time counter YYYY-MM-DD
    ln_time_interpolation = .true. !  set to false to use parent
                                   !  calender for monthly frequency only

!------------------------------------------------------------------------------
!  Additional parameters
!------------------------------------------------------------------------------
    nn_wei  = 1                   !  smoothing filter weights
    rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                  !  smoothing onto dst points. Need to
                                  !  make this a funct. of dlon
    ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
    nn_alpha    = 0               !  Euler rotation angle
    nn_beta     = 0               !  Euler rotation angle
    nn_gamma    = 0               !  Euler rotation angle
    rn_mask_max_depth = 100.0     !  Maximum depth to be ignored for the mask
    rn_mask_shelfbreak_dist = 20000.0 !  Distance from the shelf break
