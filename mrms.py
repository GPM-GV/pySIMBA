import os, sys, glob
import numpy as np
from datetime import datetime
import pdb #pdb.set_trace()
import tools as sim

def get_mrms_for_column(mrms_files, main_plat_datetime, box_params, column):

    # hard set these for MRMS files, based on documentation:
    mrms_min_lat = 20.005   #SW-most point's latitude
    mrms_min_lon = -129.995 #SW-most point's longitude
    mrms_delta   = 0.01     #[deg] spacing of each MRMS point
    mrms_nodata_val = -999.00

    # Get search area so don't have to check every MRMS point:
    limit_km = 10.0 * (box_params.box_limit/1000.)    #[km]

    # Get list of ea type of MRMS file in the dir:
    rgr_files = sorted(glob.glob(f'{mrms_files}/1HCF.*'))        #MRMS Radar-Gauge Ratio files
    rqi_files = sorted(glob.glob(f'{mrms_files}/RQI.*'))         #MRMS Radar Quality Index files
    rate_files = sorted(glob.glob(f'{mrms_files}/PRECIPRATE.*')) #MRMS Precipitation Rate files
    type_files = sorted(glob.glob(f'{mrms_files}/MASK.*'))       #MRMS Precipitation Type files
    if len(rate_files) != len(rqi_files) or len(rate_files) != len(type_files):
        print('---PROBLEM WITH MRMS DATA: DIFFERENT NUMBER OF FILES FOR EACH MRMS FILE TYPE!!---')
        print('------------  MRMS PRODUCTS ARE NOT GETTING SET INTO COLUMN GRID!! ------------')
        return column

    #define empty arrays to hold datetimes and offsets
    mrms_datetime   = np.full(len(rate_files), np.nan, dtype=object)
    mrms_offsets = np.full(len(rate_files), np.nan)

    # Start w/ precip rate files: timestamp for these, type, and RQI files will match
    fmt='%Y%m%d%H%M%S'
    for i, file in enumerate(rate_files):
        fname_split = os.path.basename(file).split('.')[2:4]
        mrms_times = fname_split[0]+'_'+fname_split[1]
        mrms_datetime[i] = datetime.strptime(fname_split[0]+fname_split[1], fmt)
        offset = (mrms_datetime[i] - main_plat_datetime).total_seconds()
        mrms_offsets[i]  = offset

    # locate which MRMS time is closest to main_plat_timestamp (offset nearest 0)
    time_sub = np.where(np.abs(mrms_offsets) == np.abs(mrms_offsets).min())[0][0]
    print(f'   TIMESTAMP TO USE: {main_plat_datetime.strftime("%Y%m%d_%H%M%S")}')
    print(f'  CLOSEST MRMS TIME: {mrms_datetime[time_sub].strftime("%Y%m%d_%H%M%S")}')

    # use the MRMS files w/ smallest offset for column grid data...first, need
    # to define arrays to hold the MRMS values for ea col grid point...these
    # will be 3D arrays, but only the surface/z dimension 0 will be populated
    # bc MRMS is just for surface:
    mrms_rate    = np.full(box_params.plat_dims, np.nan, dtype=np.float)
    mrms_type    = np.full(box_params.plat_dims, np.nan, dtype=np.float)
    mrms_rqi     = np.full(box_params.plat_dims, np.nan, dtype=np.float)
    mrms_rgr     = np.full(box_params.plat_dims, np.nan, dtype=np.float)
    mrms_pt_dist = np.full(box_params.plat_dims, np.nan, dtype=np.float)
    
    # place an unzipped copy of rate, type, rqi files in current dir
    rate_dat_file = sim.ungzip_file(rate_files[time_sub])
    type_asc_file = sim.ungzip_file(type_files[time_sub])
    rqi_asc_file = sim.ungzip_file(rqi_files[time_sub])

    # get ratios file that contains this time, if avail:
    #  MRMS ratios files are hourly, with file time the end of the period - if a 
    #  ratios file exists for the time of interest, it will have a timestamp that
    #  is for the NEXT hour, eg: 1530 will be in rgr_file w/ timestamp 1600
    #  unless the time is exactly/top of hour, in that case use the same time, eg:
    #  for interest time 1600, used rgr_file w/ timestamp 1600.
    #mrms_time_hr  = strmid(mrms_times[time_sub], 9, 2)
    #mrms_time_min = strmid(mrms_times[time_sub], 11,2)
    mrms_time_hr = mrms_datetime[time_sub].hour
    mrms_time_min = mrms_datetime[time_sub].minute
    if len(rgr_files) >= 1:
        if mrms_time_min == 0:
            use_rgr_hr = mrms_time_hr
        else:
            use_rgr_hr = mrms_time_hr + 1
    
    rgr_datetime   = np.full(len(rgr_files), np.nan, dtype=object)
    rgr_hr         = np.full(len(rgr_files), np.nan, dtype=object)
    #for i=0,n_elements(rgr_files)-1 do begin
    for i, file in enumerate(rgr_files):
        fname_split = os.path.basename(file).split('.')[1:3]
        rgr_time = fname_split[0]+'_'+fname_split[1]
        rgr_datetime[i] = datetime.strptime(fname_split[0]+fname_split[1], fmt)
        rgr_hr[i] = rgr_datetime[i].hour

    rgr_time_sub = np.where(rgr_hr == use_rgr_hr)[0]
    n_rgr_times = len(rgr_time_sub)
    if n_rgr_times > 0:
        rgr_asc_file = sim.ungzip_file(rgr_files[rgr_time_sub[0]])

    # read in MRMS files for this time:
    rate_data_all = np.loadtxt(rate_dat_file, skiprows=6)
    type_data_all = np.loadtxt(type_asc_file, skiprows=6)
    rqi_data_all = np.loadtxt(rqi_asc_file, skiprows=6)
    if n_rgr_times > 0:
        rgr_data_all  = np.loadtxt(rgr_asc_file, skiprows=6)
        os.remove(f'./{rgr_asc_file}')
    
    # now that read in from the files, delete unziped ones (.gz are still in orig/mrms_dir)
    os.remove(f'./{rate_dat_file}')
    os.remove(f'./{type_asc_file}')
    os.remove(f'./{rqi_asc_file}')

    #dimensions are (lat, lon)
    rate_dims     = rate_data_all.shape  # these should all match
    type_dims     = type_data_all.shape  # these should all match
    rqi_dims      = rqi_data_all.shape  # these should all match
    
    lats = np.arange(rate_dims[0]) * mrms_delta + mrms_min_lat #only one of these lines is correct..
    rev_lats = np.flip(lats)                                   #THIS ONE: top line in file = MAXIMUM latitude
    
    lons = np.arange(rate_dims[1]) * mrms_delta + mrms_min_lon
    range_deg_lat = limit_km / 111.1
    range_deg_lon = limit_km / (np.cos(np.deg2rad(box_params.cntr_lat_deg))*111.1)

    rev_lat_subs = np.where( (rev_lats >= box_params.cntr_lat_deg-range_deg_lat) & (rev_lats <= box_params.cntr_lat_deg+range_deg_lat) )
    lon_subs     = np.where( (lons >= box_params.cntr_lon_deg-range_deg_lon) & (lons <= box_params.cntr_lon_deg+range_deg_lon))

    rev_lat_start_sub = rev_lat_subs[0][0]
    rev_lat_end_sub   = rev_lat_subs[0][-1]
    lon_start_sub  = lon_subs[0][0]
    lon_end_sub    = lon_subs[0][-1]
    #print(f'start & end lat: {lats[lat_start_sub]} {lats[lat_end_sub]}')
    #print(f' start & end rev(lat): {rev_lats[rev_lat_start_sub]} {rev_lats[rev_lat_end_sub]}')
    #print(f' start & end lon: {lons[lon_start_sub]} {lons[lon_end_sub]}')


    # Get subsets of data & lat/lon arrays:
    lats    = rev_lats[rev_lat_start_sub:rev_lat_end_sub+1]
    lons    = lons[lon_start_sub:lon_end_sub+1]

    rate_data = rate_data_all[rev_lat_start_sub:rev_lat_end_sub+1,lon_start_sub:lon_end_sub+1]
    type_data = type_data_all[rev_lat_start_sub:rev_lat_end_sub+1,lon_start_sub:lon_end_sub+1]
    rqi_data  = rqi_data_all[rev_lat_start_sub:rev_lat_end_sub+1,lon_start_sub:lon_end_sub+1]
    if n_rgr_times > 0:
        rgr_data = rgr_data_all[rev_lat_start_sub:rev_lat_end_sub+1,lon_start_sub:lon_end_sub+1]

    # Now that have subset in vicinity of column grid, locate MRMS point nearest
    #  ea column grid point & assign those values to the column data arrays:
    #for i=0,n_elements(column_box_params.column_grid_lons)-1 do begin
    #for j=0,n_elements(column_box_params.column_grid_lats)-1 do begin
    for i, col_pt_lon in enumerate(box_params.lon_values):
        for j, col_pt_lat in enumerate(box_params.lat_values):

            # closest utility: find closest MRMS lat & lon to col pt:
            closest_lat_sub = sim.closest(lats, col_pt_lat)
            closest_lon_sub = sim.closest(lons, col_pt_lon)
            if len(closest_lat_sub) !=1 or len(closest_lon_sub) !=1:
                print('---PROBLEM SETTING WHERE IN GRID TO PLACE MRMS VALUES!!--')
                print('---Returning without MRMS data in the column---')
                return column
            else:
                closest_lat_sub = closest_lat_sub[0]
                closest_lon_sub = closest_lon_sub[0]
            
            # populate only lowest levet since MRMS is essentially a sfc level product:
            mrms_rate[0, j, i] = rate_data[closest_lat_sub, closest_lon_sub]
            mrms_type[0, j, i] = type_data[closest_lat_sub, closest_lon_sub]
            mrms_rqi[0, j, i]  = rqi_data[closest_lat_sub, closest_lon_sub]
            if n_rgr_times > 0:
                mrms_rgr[0, j, i] = rgr_data[closest_lat_sub, closest_lon_sub]
            
            # get the dist of nearest MRMS pt to the col grid pt:
            temp_val = sim.get_posn_to_azm_range(lats[closest_lat_sub], lons[closest_lon_sub], col_pt_lat, col_pt_lon)
            mrms_pt_dist[0, j, i] =  temp_val[1]

    # replace missing flags w/ NaNs:
    mrms_rate = sim.bad2nan(mrms_rate, mrms_nodata_val)
    mrms_type = sim.bad2nan(mrms_type, mrms_nodata_val)
    mrms_rqi  = sim.bad2nan(mrms_rqi, mrms_nodata_val)
    if n_rgr_times > 0:
        mrms_rgr = sim.bad2nan(mrms_rgr, mrms_nodata_val)

    #if data array is now all NaNs, then the product/file type
    #  was not available:
    n_rate_values = np.count_nonzero(np.isfinite(mrms_rate))
    n_type_values = np.count_nonzero(np.isfinite(mrms_type))
    n_rqi_values = np.count_nonzero(np.isfinite(mrms_rqi))
    n_rgr_values = np.count_nonzero(np.isfinite(mrms_rgr))
    
    avail_products_string=''
    if n_rate_values >= 1: avail_products_string = avail_products_string+'PRECIPRATE, '
    if n_type_values >= 1: avail_products_string = avail_products_string+'PRECIPTYPE, '
    if n_rqi_values  >= 1: avail_products_string = avail_products_string+'RQI, '
    if n_rgr_values  >= 1: avail_products_string = avail_products_string+'HRLYRGR'    

    # Set up structures to be returned:
    mrms_info = {'timestamp_requested':f'{main_plat_datetime.strftime("%Y%m%d_%H%M%S")}',
                     'timestamp':f'{mrms_datetime[time_sub].strftime("%Y%m%d_%H%M%S")}',
                     #offset_vs_main: mrms_offsets[time_sub],$  ; NO! this is offset vs requested timestamp!
                     'offset_vs_main': mrms_offsets[time_sub],
                     'products': avail_products_string }
                      #min_pt_dist: min(mrms_pt_dist, /nan),$
                      #max_pt_dist: max(mrms_pt_dist, /nan)}

    #add mrms_info to object
    column.add_platform_to_object(mrms_info, plat_name='mrms')

    mrms_name_precip_rate = 'MRMS Lev2 Precipitation Rate'
    mrms_units_precip_rate = 'mm h^-1'
    mrms_name_precip_type = 'MRMS Lev2 Precipitation Type'
    mrms_units_precip_type = 'type codes - see READMEs' ###will need to modify this
    mrms_name_radarQuality = 'MRMS Radar Quality Index'
    mrms_units_radarQuality = '1 is best'
    mrms_name_HrlyRadarGaugeRatio = 'MRMS hourly gauge-adjusted radar & hourly radar only'
    mrms_units_HrlRadarGaugeRatio = 'ranges 0.1-10'
    mrms_name_to_grid_dist = 'distance from MRMS point to assigned column grid point'
    mrms_units_to_grid_dist = 'm'

    #assign fields to column object
    column.add_variable_to_object(mrms_rate,
                                  var_name='MRMS_precip_rate',
                                  long_name=mrms_name_precip_rate, 
                                  units=mrms_units_precip_rate)

    column.add_variable_to_object(mrms_type,
                                  var_name='MRMS_precip_type',
                                  long_name=mrms_name_precip_type, 
                                  units=mrms_units_precip_type)

    column.add_variable_to_object(mrms_rqi,
                                  var_name='MRMS_RQI',
                                  long_name=mrms_name_radarQuality, 
                                  units=mrms_units_radarQuality)

    column.add_variable_to_object(mrms_rgr,
                                  var_name='MRMS_HrlyRadarGaugeRatio',
                                  long_name=mrms_name_HrlyRadarGaugeRatio, 
                                  units=mrms_units_HrlRadarGaugeRatio)

    column.add_variable_to_object(mrms_pt_dist,
                                  var_name='MRMS_pt_dist',
                                  long_name=mrms_name_to_grid_dist, 
                                  units=mrms_units_to_grid_dist)

    return column
    #self.mrms_data_in_column = True