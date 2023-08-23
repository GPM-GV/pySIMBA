import os
import numpy as np
from datetime import datetime
import pdb #pdb.set_trace()
from netCDF4 import Dataset 

def netcdf(main_plat_name, out_dir, box_params, column):

    #------------------------------------------------------------
    # procedure to write new .nc file w/ data & atts in column box
    #------------------------------------------------------------
    #
    #  COLUMN netCDF FILES WILL BE NAMED AS:
    #       column_[main_plat]_[center]_[YYYYMMDD_HHMM].nc
    #       [main_plat]:      string ID of the main_plat used
    #       [center]:         string ID of box's center location
    #       [YYYYMMDD_HHMM]:  main_plat's time w/o seconds
    #
    #   example:  column_NPOL_WFFpad_20160506_0617.nc
    #
    #  Also generates a .txt file with the column file's netCDF header info
    #   
    #   example:  column_NPOL_WFFpad_20160506_0617_ncHeader.txt
    #
    #------------------------------------------------------------

    # Set up the new filenames:
    column_nc_fname = 'column_'+main_plat_name+'_'+ \
            box_params.box_center+'_'+column.time[0:13]+'.nc'
    header_fname = os.path.basename(column_nc_fname)+'_ncHeader.txt'
    print('------------------------------------')
    print('---Now Writing New Column .nc File: ')
    print(f'          {column_nc_fname}')
    column_nc_fname = out_dir+column_nc_fname
    header_fname    = out_dir+header_fname

    if not box_params.have_xy_values:
        print('------------------------------------------------------------------')
        print('---DID NOT GET CORRECT x & y ARRAYS TO USE FOR COLUMN .nc FILE!---')
        print('---------COLUMN NetCDF FILE HAS NOT BEEN GENERATED!!!-------------')
        print('------------------------------------------------------------------')
        return -1

    # Begin new .nc file, only if one with same name doesn't already exist:
    ncfile = Dataset(column_nc_fname,mode='w', noclobber=False)

    # Define the X, Y, & Z dimension sizes for the column box grid
    #  (# of box edges, not the number of boxes)
    xdim = ncfile.createDimension('x', box_params.nx)
    ydim = ncfile.createDimension('y', box_params.ny)
    zdim = ncfile.createDimension('z', box_params.nz)
    tdim = ncfile.createDimension('t', box_params.nt)
    #ndim = ncfile.createDimension('n', None) #this dimension is for dropcounts in 2dvd, apu 
    
    # Set global attributes:  include main plat & info for each plat
    ncfile.setncattr('box_centered_on', box_params.box_center)
    ncfile.setncattr('box_center_lat', box_params.cntr_lat_deg)
    ncfile.setncattr('box_center_lon', box_params.cntr_lon_deg)
    ncfile.setncattr('grid_spacing_vert', box_params.vert_spacing)
    ncfile.setncattr('grid_spacing_horiz', box_params.box_spacing)
    ncfile.setncattr('grid_extent_vert', box_params.vert_limit)
    ncfile.setncattr('grid_extent_horiz', box_params.box_limit)
    ncfile.setncattr('grid_spacing_and_limits_units', 'meters')
    ncfile.setncattr('main_platform', main_plat_name)
    ncfile.setncattr('main_plat_timestamp', column.time)
    #pdb.set_trace()
    #ncfile.setncattr('main_plat_mode', self.main_scan_type) ??
    ncfile.setncattr('pySIMBA_version', 'v0.1')
    
    #write coordinate data
    x_ID = ncfile.createVariable('x', np.float32, ('x'), zlib=True)
    x_ID[:] = box_params.x_values
    
    y_ID = ncfile.createVariable('y', np.float32, ('y'), zlib=True)
    y_ID[:] = box_params.y_values
    
    z_ID = ncfile.createVariable('z', np.float32, ('z'), zlib=True)
    z_ID[:] = box_params.z_values
    #pdb.set_trace()
    t_ID = ncfile.createVariable('t', np.float32, ('t'), zlib=True)
    t_ID[:] = box_params.t_values
    
    latID = ncfile.createVariable('lat', np.float32, ('y'), zlib=True)
    latID[:] = box_params.lat_values
    
    lonID = ncfile.createVariable('lon', np.float32, ('x'), zlib=True)
    lonID[:] = box_params.lon_values

    npol_avail = ncfile.createVariable('npol_avail', str , (), zlib=True)
    #d3r_avail = ncfile.createVariable('d3r_avail', str , (), zlib=True)
    lev2_avail = ncfile.createVariable('lev2_avail', str , (), zlib=True)
    apu_avail = ncfile.createVariable('apu_avail', str , (), zlib=True)
    twodvd_avail = ncfile.createVariable('twodvd_avail', str , (), zlib=True)
    #dow6_avail = ncfile.createVariable('dow6_avail', str , (), zlib=True)
    #pluvio_avail = ncfile.createVariable('pluvio_avail', str , (), zlib=True)
    gauges_avail = ncfile.createVariable('gauges_avail', str , (), zlib=True)
    mrr_avail = ncfile.createVariable('mrr_avail', str , (), zlib=True)
    #gmi_gprof_avail = ncfile.createVariable('gmi_gprof_avail', str , (), zlib=True)
    #gmi_L1C_avail = ncfile.createVariable('gmi_L1C_avail', str , (), zlib=True)
    dpr_avail = ncfile.createVariable('dpr_avail', str , (), zlib=True)
    #lev2bcmb_avail = ncfile.createVariable('lev2bcmb_avail', str , (), zlib=True)
    mrms_avail = ncfile.createVariable('mrms_avail', str , (), zlib=True)
    #soundings_avail = ncfile.createVariable('soundings_avail', str , (), zlib=True)
    
    #set all platforms to False available by default
    npol_avail.setncattr('npol_available', 'False')
    npol_avail = set_nc_radar_missing_attributes(npol_avail)
    
    apu_avail.setncattr('apu_available', 'False')
    apu_avail = set_nc_disdro_missing_attributes(apu_avail)
    
    twodvd_avail.setncattr('twodvd_available', 'False')
    twodvd_avail = set_nc_disdro_missing_attributes(twodvd_avail)
    
    gauges_avail.setncattr('gauges_available', 'False')
    gauges_avail = set_nc_disdro_missing_attributes(gauges_avail)
    
    mrr_avail.setncattr('mrr_available', 'False')
    mrr_avail = set_nc_mrr_missing_attributes(mrr_avail)
    
    dpr_avail.setncattr('dpr_available', 'False')
    dpr_avail = set_nc_satellite_missing_attributes(dpr_avail)
    
    mrms_avail.setncattr('MRMS_available', 'False')
    mrms_avail = set_mrms_missing_attributes(mrms_avail)
    
    #-------------------------------------------------------------------------------
    # Check if each platform data is available to input into column
    #  not available: Set all the attributes to -9999 or 'platfrom_not_avail'
    #  available:     Define variables & coordinates
    #                 Write data to the NETCDF File/Group Defined Above
    #-------------------------------------------------------------------------------
    
    for plat in column.platforms.keys():
        # SET PLATFORM ATTS FOR NPOL/LEV2 RADAR:
        
        if 'lev2_radars' in plat.lower():
            plat_name = plat.lower().split('_')[1]
            #lev2_avail.setncattr(f'{plat_name}_available', 'True')
            lev2_avail = set_nc_radar_attributes(lev2_avail, column.platforms[plat])
        #if 'npol' in plat.lower():
        #    print('npol is here')
        #    npol_avail.setncattr('npol_available', 'True')
        #    npol_avail = set_nc_radar_attributes(npol_avail, column.platforms['NPOL'])

        # SET PLATFORM ATTS FOR PARSIVEL:
        if 'parsivel' in column.platforms.keys():
            apu_avail.setncattr('apu_available', 'True')
            apu_avail = set_nc_disdro_attributes(apu_avail, column.platforms['parsivel'])

        # SET PLATFORM ATTS FOR 2DVD:
        if 'twodvd' in column.platforms.keys():
            twodvd_avail.setncattr('twodvd_available', 'True')
            twodvd_avail = set_nc_disdro_attributes(twodvd_avail, column.platforms['twodvd'])

        # SET PLATFORM ATTS FOR GAUGES:
        if 'gauges' in column.platforms.keys():
            gauges_avail.setncattr('gauges_available', 'True')
            gauges_avail = set_nc_disdro_attributes(gauges_avail, column.platforms['gauges'])

        # SET PLATFORM ATTS FOR MRR:
        if 'mrr' in column.platforms.keys():
            mrr_avail.setncattr('mrr_available', 'True')
            mrr_avail = set_nc_mrr_attributes(mrr_avail, column.platforms['mrr'])

        # SET PLATFORM ATTS FOR DPR:
        if 'dpr' in column.platforms.keys():
            dpr_avail.setncattr('dpr_available', 'True')
            dpr_avail = set_nc_satellite_attributes(dpr_avail, column.platforms['dpr'])

        # SET PLATFORM ATTS FOR MRMS:
        if 'mrms' in column.platforms.keys():
            mrms_avail.setncattr('MRMS_available', 'True')
            mrms_avail = set_mrms_attributes(mrms_avail, column.platforms['mrms'])
            
        # SET PLATFORM ATTS FOR [OTHER PLATFORMS]:

    for key in column.variables.keys():
        # SET VAR ATTS FOR NPOL/LEV2 RADAR DATA:
        if 'lev2' in key.lower():
            lev2_key_ID = ncfile.createVariable(key, np.float32, ('z','y','x'), zlib=True)
            lev2_key_ID.long_name = column.variables[key]['long_name']
            lev2_key_ID.units = column.variables[key]['units']
            lev2_key_ID[:,:,:] = column.variables[key]['data']
            lev2_key_ID.badval = column.variables[key]['badval']
        
        #if 'npol' in key.lower():
        #    npol_key_ID = ncfile.createVariable(key, np.float32, ('z','y','x'), zlib=True)
        #    npol_key_ID.long_name = column.variables[key]['long_name']
        #    npol_key_ID.units = column.variables[key]['units']
        #    npol_key_ID[:,:,:] = column.variables[key]['data']
        #    npol_key_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR PARSIVEL DATA:
        if 'apu' in key.lower():
            apu_ID = ncfile.createVariable(key, np.float32, ('z','y','x','t'), zlib=True)
            apu_ID.units = column.variables[key]['units']
            apu_ID[:,:,:,:] = column.variables[key]['data']
            apu_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR 2DVD DATA:
        if 'twodvd' in key.lower():
            twodvd_ID = ncfile.createVariable(key, np.float32, ('z','y','x','t'), zlib=True)
            twodvd_ID.units = column.variables[key]['units']
            twodvd_ID[:,:,:,:] = column.variables[key]['data']
            twodvd_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR GAUGES DATA:
        if 'gauges' in key.lower():
            gauges_ID = ncfile.createVariable(key, np.float32, ('z','y','x','t'), zlib=True)
            gauges_ID.units = column.variables[key]['units']
            gauges_ID[:,:,:,:] = column.variables[key]['data']
            gauges_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR MRR DATA:
        if 'mrr' in key.lower():
            mrr_ID = ncfile.createVariable(key, np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_ID.long_name = column.variables[key]['long_name']
            mrr_ID.units = column.variables[key]['units']
            mrr_ID[:,:,:,:] = column.variables[key]['data']
            mrr_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR DPR DATA:
        if 'dpr' in key.lower():
            dpr_ID = ncfile.createVariable(key, np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_ID.long_name = column.variables[key]['long_name']
            dpr_ID.units = column.variables[key]['units']
            dpr_ID[:,:,:] = column.variables[key]['data']
            dpr_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR MRMS DATA:
        if 'mrms' in key.lower():
            mrms_ID = ncfile.createVariable(key, np.float32, ('z','y','x'), zlib=True)
            mrms_ID.long_name = column.variables[key]['long_name']
            mrms_ID.units = column.variables[key]['units']
            mrms_ID[:,:,:] = column.variables[key]['data']
            mrms_ID.badval = column.variables[key]['badval']

        # SET VAR ATTS FOR [OTHER PLATFORMS] DATA:

    # close the newly made .nc file
    ncfile.close()
    print(f'--> column file has been saved {column_nc_fname}')

    # Now that have saved the column .nc file, genearte
    #  a .txt file of the .nc header information:
    header_cmnd = 'ncdump -h '+column_nc_fname+' >& '+header_fname
    os.system(header_cmnd)
    print(f'--> column header .txt file has been saved {header_fname}')
    
# ***************************************************************************************

def set_nc_disdro_attributes(plat_grp, plat_info):
    # need to loop thru ****_plat_info arrays' elements bc may have multiple platforms avail:
    for a, name in enumerate(plat_info['plat_name']):
        plat_grp.setncattr(name+'_lat_deg',  plat_info['lat_d'][a])
        plat_grp.setncattr(name+'_lat_min',  plat_info['lat_m'][a])
        plat_grp.setncattr(name+'_lat_sec',  plat_info['lat_s'][a])
        plat_grp.setncattr(name+'_lon_deg',  plat_info['lon_d'][a])
        plat_grp.setncattr(name+'_lon_min',  plat_info['lon_m'][a])
        plat_grp.setncattr(name+'_lon_sec',  plat_info['lon_s'][a])
        plat_grp.setncattr(name+'_operation_mode',  plat_info['operation_mode'][a])
        plat_grp.setncattr(name+'_timestamp',plat_info['timestamp'][a])
        plat_grp.setncattr(name+'_offset_vs_main',   plat_info['offset_vs_main'][a])
        plat_grp.setncattr(name+'_time_interval_width',plat_info['time_interval_width'][a])
    return plat_grp

# ***************************************************************************************

def set_nc_disdro_missing_attributes(plat_grp):
    #set platform attributes to missing
    plat_grp.setncattr('latitude_degrees',  -9999)
    plat_grp.setncattr('latitude_minutes',  -9999)
    plat_grp.setncattr('latitude_seconds',  -9999)
    plat_grp.setncattr('longitude_degrees', -9999)
    plat_grp.setncattr('longitude_minutes', -9999)
    plat_grp.setncattr('longitude_seconds', -9999)
    plat_grp.setncattr('operation_mode',    'platform_not_avail')
    plat_grp.setncattr('timestamp',         'platform_not_avail')
    plat_grp.setncattr('offset_vs_main',    'platform_not_avail')
    plat_grp.setncattr('time_interval_width',-9999)
    return plat_grp

# ***************************************************************************************

def set_nc_mrr_attributes(plat_grp, plat_info):
    # need to loop thru ****_plat_info arrays' elements bc may have multiple platforms avail:
    for a, name in enumerate(plat_info['plat_name']):
        plat_grp.setncattr(name+'_lat_deg',  plat_info['lat_d'][a])
        plat_grp.setncattr(name+'_lat_min',  plat_info['lat_m'][a])
        plat_grp.setncattr(name+'_lat_sec',  plat_info['lat_s'][a])
        plat_grp.setncattr(name+'_lon_deg',  plat_info['lon_d'][a])
        plat_grp.setncattr(name+'_lon_min',  plat_info['lon_m'][a])
        plat_grp.setncattr(name+'_lon_sec',  plat_info['lon_s'][a])
        plat_grp.setncattr(name+'_operation_mode',  plat_info['operation_mode'][a])
        plat_grp.setncattr(name+'_elevation_MSL', plat_info['elev'][a])
        plat_grp.setncattr(name+'_wavelength_m', plat_info['wavelength'][a])
        plat_grp.setncattr(name+'_frequency_GHz', plat_info['frequency'][a])
        plat_grp.setncattr(name+'_beam_width_deg', plat_info['beam_width'][a])
        plat_grp.setncattr(name+'_gate_size_m', plat_info['gate_size'][a])
        plat_grp.setncattr(name+'_timestamp',plat_info['timestamp'][a])
        plat_grp.setncattr(name+'_offset_vs_main',   plat_info['offset_vs_main'][a])
        plat_grp.setncattr(name+'_time_interval_width',plat_info['time_interval_width'][a])
    return plat_grp

# ***************************************************************************************

def set_nc_mrr_missing_attributes(plat_grp):
    #set platform attributes to missing
    plat_grp.setncattr('latitude_degrees',  -9999)
    plat_grp.setncattr('latitude_minutes',  -9999)
    plat_grp.setncattr('latitude_seconds',  -9999)
    plat_grp.setncattr('longitude_degrees', -9999)
    plat_grp.setncattr('longitude_minutes', -9999)
    plat_grp.setncattr('longitude_seconds', -9999)
    plat_grp.setncattr('operation_mode',    'platform_not_avail')
    plat_grp.setncattr('elevation_MSL', 'platform_not_available')
    plat_grp.setncattr('wavelength_m', 'platform_not_available')
    plat_grp.setncattr('frequency_GHz', 'platform_not_available')
    plat_grp.setncattr('beam_width_deg', 'platform_not_available')
    plat_grp.setncattr('gate_size_m', 'platform_not_available')
    plat_grp.setncattr('timestamp',         'platform_not_avail')
    plat_grp.setncattr('offset_vs_main',    'platform_not_avail')
    plat_grp.setncattr('time_interval_width',-9999)
    return plat_grp

# ***************************************************************************************

def set_nc_satellite_attributes(plat_grp, plat_info):
    #name = plat_info_data[plat_info_key]['plat_name']
    plat_grp.setncattr('file_type', plat_info['file_type'])
    plat_grp.setncattr('frequency', plat_info['frequency'])
    plat_grp.setncattr('algorithm', plat_info['algorithm'])
    plat_grp.setncattr('orbit_number', plat_info['orbit_num'])
    plat_grp.setncattr('data_version', plat_info['data_version'])
    plat_grp.setncattr('timestamp', plat_info['timestamp'])
    plat_grp.setncattr('offset_vs_main', plat_info['offset_vs_main'])
    plat_grp.setncattr('timestamp_cntr', plat_info['timestamp_cntr'])
    return plat_grp

# ***************************************************************************************

def set_nc_satellite_missing_attributes(plat_grp):
    #name = plat_info_key.split('_')[0]
    plat_grp.setncattr('file_type', 'platform_not_available')
    plat_grp.setncattr('frequency', -9999)
    plat_grp.setncattr('algorithm', -9999)
    plat_grp.setncattr('orbit_number', -9999)
    plat_grp.setncattr('data_version', -9999)
    plat_grp.setncattr('timestamp', 'platform_not_available')
    plat_grp.setncattr('offset_vs_main', 'platform_not_available')
    plat_grp.setncattr('timestamp_cntr', 'platform_not_available')
    return plat_grp

# ***************************************************************************************

def set_nc_radar_attributes(plat_grp, plat_info):
    # need to loop thru ****_plat_info arrays' elements bc may have multiple platforms avail:
    for a, name in enumerate(plat_info['plat_name']):
        plat_grp.setncattr(name+'_available', 'True')
        plat_grp.setncattr(name+'_latitude_degrees',  plat_info['lat_d'][a])
        plat_grp.setncattr(name+'_latitude_minutes',  plat_info['lat_m'][a])
        plat_grp.setncattr(name+'_latitude_seconds',  plat_info['lat_s'][a])
        plat_grp.setncattr(name+'_longitude_degrees',  plat_info['lon_d'][a])
        plat_grp.setncattr(name+'_longitude_minutes',  plat_info['lon_m'][a])
        plat_grp.setncattr(name+'_longitude_seconds',  plat_info['lon_s'][a])
        plat_grp.setncattr(name+'_elevation_MSL', plat_info['elev'][a])
        plat_grp.setncattr(name+'_operation_mode',  plat_info['operation_mode'][a])
        plat_grp.setncattr(name+'_wavelength_m', plat_info['wavelength'][a])
        plat_grp.setncattr(name+'_frequency_GHz', plat_info['frequency'][a])
        plat_grp.setncattr(name+'_beam_width_deg', plat_info['beam_width'][a])
        plat_grp.setncattr(name+'_gate_size_m', plat_info['gate_size'][a])
        plat_grp.setncattr(name+'_timestamp',plat_info['timestamp'][a])
        plat_grp.setncattr(name+'_offset_vs_main',   plat_info['offset_vs_main'][a])
    return plat_grp

# ***************************************************************************************

def set_nc_radar_missing_attributes(plat_grp):
    plat_grp.setncattr('latitude_degrees',  -9999)
    plat_grp.setncattr('latitude_minutes',  -9999)
    plat_grp.setncattr('latitude_seconds',  -9999)
    plat_grp.setncattr('longitude_degrees', -9999)
    plat_grp.setncattr('longitude_minutes', -9999)
    plat_grp.setncattr('longitude_seconds', -9999)
    plat_grp.setncattr('elevation_MSL',     'platform_not_avail')
    plat_grp.setncattr('operation_mode',    'platform_not_avail')
    plat_grp.setncattr('wavelength_m',      'platform_not_avail')
    plat_grp.setncattr('frequency_GHz',     'platform_not_avail')
    plat_grp.setncattr('beam_width_deg',    'platform_not_avail')
    plat_grp.setncattr('gate_size_m',       'platform_not_avail')
    plat_grp.setncattr('timestamp',         'platform_not_avail')
    plat_grp.setncattr('offset_vs_main',    'platform_not_avail')
    return plat_grp

# ***************************************************************************************

def set_mrms_attributes(plat_grp, plat_info):
    plat_grp.setncattr('timestamp_requested', plat_info['timestamp_requested'])
    plat_grp.setncattr('timestamp', plat_info['timestamp'])
    plat_grp.setncattr('offset_vs_main', plat_info['offset_vs_main'])
    plat_grp.setncattr('products_avail', plat_info['products'])
    return plat_grp

# ***************************************************************************************

def set_mrms_missing_attributes(plat_grp):
    plat_grp.setncattr('timestamp_requested', 'platform_not_available')
    plat_grp.setncattr('timestamp', -9999)
    plat_grp.setncattr('offset_vs_main', -9999)
    plat_grp.setncattr('products_avail', 'platform_not_available')
    return plat_grp

# ***************************************************************************************