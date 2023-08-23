#from timeit import default_timer as timer
import os, sys
import numpy as np
from datetime import datetime
import pdb #pdb.set_trace()
import pyart
import tools as sim

def grid_radar_for_column(radar_file, main_plat_datetime, radarID, box_params, column):

    #------prepare data for gridding-------

    #unzip file (if needed)
    #file_basename = os.path.basename(self.npol)
    #pdb.set_trace()
    #file_basename = os.path.basename(self.npol_file)
    #if file_basename.endswith('.gz'):
    #    cf_file = sim.ungzip_file(self.npol_file)
    #else:
    #    cf_file = self.npol_file
    
    #fields = ['ZZ','CZ','DR','RH','PH','KD','SQ','SW','VR',\
    #          'RR','RP','RC','D0','NW','FH','N2','DM']
    #radar_id =os.path.basename(radar_file).split('_')[0][0:4].lower()
    grid_center_lat_lon = box_params.cntr_lat_deg, box_params.cntr_lon_deg
    print('...running Py-ART gridding...')
    #Beg_timer = timer()
    #grid = self.get_grid_from_radar(radar_file, grid_center_lat_lon)
    grid = get_grid_from_radar(radar_file, grid_center_lat_lon, box_params.vert_spacing, box_params.box_spacing)
    #End_timer = timer()
    #print(f”{Process took (End_timer-Beg_timer):0.2f} seconds”)
    #print(f'Process took {End_timer-Beg_timer:0.2f} seconds')

    #Py-Art grid output will be put in a tempoarilly created dir within
    #main output dir, when this method is done the full grid
    #.nc file will get moved up to this main output dir and the 
    #temporary dir will get deleted.
    #save the full grid .nc file
    #full_grid_filename = f'{file_basename.split(".")[0]}_cntron_{self.center_on}'
    #full_grid_name = f'{self.full_grid_dir}{full_grid_filename}.nc'
    DS = grid.to_xarray()
    #pdb.set_trace()
    #print(f'printing DS.x:    {DS.x}')
    #print(f'printing DS.lon:  {DS.lon}')
    #DS = DS.swap_dims({"x": "lon"})
    #DS = DS.swap_dims({"y": "lat"})
    #DS.to_netcdf(full_grid_name, format='NETCDF4')
    
    #get radar parameters
    radar = pyart.io.read(radar_file)
    r_lat_decdeg = radar.latitude['data'][0]
    r_lon_decdeg = radar.longitude['data'][0]
    plat_lat_d, plat_lat_m, plat_lat_s, plat_lon_d, plat_lon_m, plat_lon_s = sim.dd2dms(r_lat_decdeg, r_lon_decdeg)
    plat_elev  = radar.altitude['data'][0] #[m]
    
    #sometimes these are not included so will put nan for those
    if 'frequency' in radar.instrument_parameters.keys():
        plat_freq = radar.instrument_parameters['frequency']['data'][0] #[Hz]
    else:
        plat_freq = np.nan
    plat_wave  = 2.998E8 / plat_freq #[m] (c/freq)
    if 'radar_beam_width_h' in radar.instrument_parameters.keys():
        plat_beam = radar.instrument_parameters['radar_beam_width_h']['data'][0] #[deg]
    else:
        plat_beam = np.nan
    if 'meters_between_gates' in radar.range.keys():
        plat_gate = radar.range['meters_between_gates'] #[m]
    else:
        plat_gate = np.nan
    plat_file_type = radar.scan_type.upper()
    radar_datetime = pyart.util.datetime_from_radar(radar)
    plat_timestamp = radar_datetime.strftime("%Y%m%d_%H%M%S")
    offset = (radar_datetime - main_plat_datetime).total_seconds()

    #plat_info dictionary
    radar_info = {'lat_d':plat_lat_d, 'lat_m':plat_lat_m, 'lat_s':plat_lat_s,
                'lon_d':plat_lon_d, 'lon_m':plat_lon_m, 'lon_s':plat_lon_s,
                'elev':plat_elev,
                'plat_name':radarID,
                'plat_type':'radar',
                'operation_mode':plat_file_type,
                'wavelength':f'{plat_wave:.5f}',
                'frequency':f'{plat_freq/1E9:.5f}',
                'beam_width':plat_beam,
                'gate_size':plat_gate,
                'timestamp':plat_timestamp,
                'offset_vs_main':f'{offset:.1f}'}

    #remove unzipped file if needed
    file_basename = os.path.basename(radar_file)
    if file_basename.endswith('.cf'):
        os.remove(radar_file)

    #extract fields & attributes 
    radar_x = DS.x
    radar_y = DS.y
    radar_z = DS.z
    radar_time = DS.time
    if 'ZZ' in DS.keys():
        radar_ZZ = DS.ZZ
        radar_ZZ_badval = DS.ZZ.attrs['_FillValue']
    elif 'DZ' in DS.keys():
        radar_ZZ = DS.DZ
        radar_ZZ_badval = DS.DZ.attrs['_FillValue']
    if 'CZ' in DS.keys():
        radar_CZ = DS.CZ
        radar_CZ_badval = DS.CZ.attrs['_FillValue']
        radar_CZ_name = DS.CZ.attrs['long_name'] 
        radar_CZ_units = DS.CZ.attrs['units']
    if 'DR' in DS.keys():
        radar_DR = DS.DR
        radar_DR_badval = DS.DR.attrs['_FillValue']
        radar_DR_name = DS.DR.attrs['long_name'] 
        radar_DR_units = DS.DR.attrs['units']
    if 'RH' in DS.keys():
        radar_RH = DS.RH
        radar_RH_badval = DS.RH.attrs['_FillValue']
        radar_RH_name = DS.RH.attrs['long_name'] 
        #radar_RH_units = DS.RH.attrs['units']
    if 'PH' in DS.keys():
        radar_PH = DS.PH
        radar_PH_badval = DS.PH.attrs['_FillValue']
        radar_PH_name = DS.PH.attrs['long_name'] 
        radar_PH_units = DS.PH.attrs['units']
    if 'KD' in DS.keys():
        radar_KD = DS.KD
        radar_KD_badval = DS.KD.attrs['_FillValue']
        radar_KD_name = DS.KD.attrs['long_name'] 
        radar_KD_units = DS.KD.attrs['units']
    if 'SQ' in DS.keys(): 
        radar_SQ = DS.SQ
        radar_SQ_badval = DS.SQ.attrs['_FillValue']
        #radar_SQ_name = DS.SQ.attrs['long_name'] 
        #radar_SQ_units = DS.SQ.attrs['units']
    if 'SW' in DS.keys():
        radar_SW = DS.SW
        radar_SW_badval = DS.SW.attrs['_FillValue']
        radar_SW_name = DS.SW.attrs['long_name'] 
        radar_SW_units = DS.SW.attrs['units']
    if 'VR' in DS.keys():
        radar_VR = DS.VR
        radar_VR_badval = DS.VR.attrs['_FillValue']
        radar_VR_name = DS.VR.attrs['long_name'] 
        radar_VR_units = DS.VR.attrs['units']
    if 'RR' in DS.keys():
        radar_RR = DS.RR
        radar_RR_badval = DS.RR.attrs['_FillValue']
        #radar_RR_name = DS.RR.attrs['long_name'] 
        #radar_RR_units = DS.RR.attrs['units']
    if 'RP' in DS.keys():
        radar_RP = DS.RP
        radar_RP_badval = DS.RP.attrs['_FillValue']
        #radar_RP_name = DS.RP.attrs['long_name'] 
        #radar_RP_units = DS.RP.attrs['units']
    if 'RC' in DS.keys():
        radar_RC = DS.RC
        radar_RC_badval = DS.RC.attrs['_FillValue']
        #radar_RC_name = DS.RC.attrs['long_name'] 
        #radar_RC_units = DS.RC.attrs['units']
    if 'D0' in DS.keys(): 
        radar_D0 = DS.D0
        radar_D0_badval = DS.D0.attrs['_FillValue']
        #radar_D0_name = DS.D0.attrs['long_name'] 
        #radar_D0_units = DS.D0.attrs['units']
    if 'NW' in DS.keys():
        radar_NW = DS.NW
        radar_NW_badval = DS.NW.attrs['_FillValue']
        #radar_NW_name = DS.NW.attrs['long_name'] 
        #radar_NW_units = DS.NW.attrs['units']
    if 'FH' in DS.keys():
        radar_FH = DS.FH
        radar_FH_badval = DS.FH.attrs['_FillValue']
        #radar_FH_name = DS.FH.attrs['long_name'] 
        #radar_FH_units = DS.FH.attrs['units']
    if 'N2' in DS.keys():
        radar_N2 = DS.N2
        radar_N2_badval = DS.N2.attrs['_FillValue']
        #radar_N2_name = DS.N2.attrs['long_name'] 
        #radar_N2_units = DS.N2.attrs['units']
    if 'DM' in DS.keys():
        radar_DM = DS.DM
        radar_DM_badval = DS.DM.attrs['_FillValue']
        #radar_DM_name = DS.DM.attrs['long_name'] 
        #radar_DM_units = DS.DM.attrs['units']
    if 'FZ' in DS.keys():
        radar_FZ = DS.FZ
        radar_FZ_badval = DS.FZ.attrs['_FillValue']
        #radar_FZ_name = DS.FZ.attrs['long_name'] 
        #radar_FZ_units = DS.FZ.attrs['units']
    if 'MW' in DS.keys():
        radar_MW = DS.MW
        radar_MW_badval = DS.MW.attrs['_FillValue']
        #radar_MW_name = DS.MW.attrs['long_name'] 
        #radar_MW_units = DS.MW.attrs['units']
    if 'MI' in DS.keys():
        radar_MI = DS.MI
        radar_MI_badval = DS.MI.attrs['_FillValue']
        #radar_MI_name = DS.MI.attrs['long_name'] 
        #radar_MI_units = DS.MI.attrs['units']

    #extract variable's attributes
    radar_time_name = 'time'
    radar_time_units = 'Coordinated Universal Time (UTC)'
    radar_time_comment = DS.time.data[0]

    radar_x_name = DS.x.attrs['standard_name']
    radar_x_units = DS.x.attrs['units']
    radar_y_name = DS.y.attrs['standard_name']
    radar_y_units = DS.y.attrs['units']
    radar_z_name = DS.z.attrs['standard_name']
    radar_z_units = DS.z.attrs['units']

    
    # Convert string attributes back into string data types:
    # most already are string type
    # Set these manually since names/units from Py-Art don't always default to what
    # NPOL fields actually are (see NPOL fields README from Jason Pippitt)
    radar_ZZ_name = 'Uncorrected Reflectivity' 
    radar_CZ_name = 'Corrected Reflectivity'
    radar_SQ_name = 'Signal Quality Index'
    radar_RR_name = 'Rain rate via DROPS2'
    radar_RP_name = 'Rain rate via PolZ-R'
    radar_RC_name = 'Rain rate via Cifelli et al 2002'
    radar_D0_name = 'median drop diameter'
    radar_NW_name = 'normalized intercept parameter (Dm)'
    radar_FH_name = 'hydrometeor ID'
    radar_N2_name = 'normalized intercept parameter (Do)'
    radar_DM_name = 'mass weighted mean diameter'
    radar_FZ_name = 's-ku frequency corrected reflectivity'
    radar_MW_name = 'Liquid Water Mass'
    radar_MI_name = 'Ice Water Mass'
    
    radar_Z_units = 'dBZ'
    radar_SQ_units = 'unitless'
    radar_RH_units = 'unitless'
    radar_R_units = 'mm h^-1'
    radar_D_units = 'mm'
    radar_NW_units = 'log(Nw)'
    radar_FH_units = 'categorical'
    radar_N2_units = 'log(Nw)'
    radar_MW_units= 'g m^-3'
    radar_MI_units = 'g m^-3'

    # Replace flagged/bad data values in each array with NaNs so
    # not using the bad data flag value in interpolation
    if 'ZZ' in DS.keys():radar_ZZ = sim.bad2nan(radar_ZZ, radar_ZZ_badval)
    if 'CZ' in DS.keys():radar_CZ = sim.bad2nan(radar_CZ, radar_CZ_badval)
    if 'DR' in DS.keys():radar_DR = sim.bad2nan(radar_DR, radar_DR_badval)
    if 'RH' in DS.keys():radar_RH = sim.bad2nan(radar_RH, radar_RH_badval)
    if 'PH' in DS.keys():radar_PH = sim.bad2nan(radar_PH, radar_PH_badval)
    if 'KD' in DS.keys():radar_KD = sim.bad2nan(radar_KD, radar_KD_badval)
    if 'SQ' in DS.keys(): radar_SQ = sim.bad2nan(radar_SQ, radar_SQ_badval)
    if 'SW' in DS.keys():radar_SW = sim.bad2nan(radar_SW, radar_SW_badval)
    if 'VR' in DS.keys():radar_VR = sim.bad2nan(radar_VR, radar_VR_badval)
    if 'RR' in DS.keys(): radar_RR = sim.bad2nan(radar_RR, radar_RR_badval)
    if 'RP' in DS.keys():radar_RP = sim.bad2nan(radar_RP, radar_RP_badval)
    if 'RC' in DS.keys():radar_RC = sim.bad2nan(radar_RC, radar_RC_badval)
    if 'D0' in DS.keys(): radar_D0 = sim.bad2nan(radar_D0, radar_D0_badval)
    if 'NW' in DS.keys():radar_NW = sim.bad2nan(radar_NW, radar_NW_badval)
    if 'FH' in DS.keys():radar_FH = sim.bad2nan(radar_FH, radar_FH_badval)
    if 'N2' in DS.keys(): radar_N2 = sim.bad2nan(radar_N2, radar_N2_badval)
    if 'DM' in DS.keys():radar_DM = sim.bad2nan(radar_DM, radar_DM_badval)
    if 'FZ' in DS.keys(): radar_FZ = sim.bad2nan(radar_FZ, radar_FZ_badval)
    if 'MW' in DS.keys(): radar_MW = sim.bad2nan(radar_MW, radar_MW_badval)
    if 'MI' in DS.keys(): radar_MI = sim.bad2nan(radar_MI, radar_MI_badval)

    #Py-ART gridding adds values at height 0 km so need to force these to NaNs
    zero_ht_sub = np.where(radar_z == 0.0)[0]
    if zero_ht_sub.shape[0] > 0:
        if zero_ht_sub.shape[0] > 1:
            print('---MORE THAN 1 HEIGHT ZERO IN THE GRIDDED FILE!---')
            print(f'---Returning without {radarID} data in the column---')
            return column
        else:
            # there is only 1 height sub for 0 km AGL
            # set all field values at ht of 0 km to NaNs
            # dimensions for radar_xx are [time, z, lat, lon] --> [time, z, y, x]; time dim is 1
            if 'ZZ' in DS.keys():radar_ZZ[0,zero_ht_sub,:,:] = np.nan
            if 'CZ' in DS.keys():radar_CZ[0,zero_ht_sub,:,:] = np.nan
            if 'DR' in DS.keys():radar_DR[0,zero_ht_sub,:,:] = np.nan
            if 'RH' in DS.keys():radar_RH[0,zero_ht_sub,:,:] = np.nan
            if 'PH' in DS.keys():radar_PH[0,zero_ht_sub,:,:] = np.nan
            if 'KD' in DS.keys():radar_KD[0,zero_ht_sub,:,:] = np.nan
            if 'SW' in DS.keys():radar_SW[0,zero_ht_sub,:,:] = np.nan
            if 'NW' in DS.keys():radar_NW[0,zero_ht_sub,:,:] = np.nan
            if 'VR' in DS.keys():radar_VR[0,zero_ht_sub,:,:] = np.nan
            if 'RC' in DS.keys():radar_RC[0,zero_ht_sub,:,:] = np.nan
            if 'RP' in DS.keys():radar_RP[0,zero_ht_sub,:,:] = np.nan
            if 'DM' in DS.keys():radar_DM[0,zero_ht_sub,:,:] = np.nan
            if 'FH' in DS.keys():radar_FH[0,zero_ht_sub,:,:] = np.nan
            if 'SQ' in DS.keys(): radar_SQ[0,zero_ht_sub,:,:] = np.nan
            if 'RR' in DS.keys(): radar_RR[0,zero_ht_sub,:,:] = np.nan
            if 'D0' in DS.keys(): radar_D0[0,zero_ht_sub,:,:] = np.nan
            if 'N2' in DS.keys(): radar_N2[0,zero_ht_sub,:,:] = np.nan
            if 'FZ' in DS.keys(): radar_FZ[0,zero_ht_sub,:,:] = np.nan
            if 'MW' in DS.keys(): radar_MW[0,zero_ht_sub,:,:] = np.nan
            if 'MI' in DS.keys(): radar_MI[0,zero_ht_sub,:,:] = np.nan

    # Locate where column box data is within the full grid of data, 
    #  then pull out the subset of data within the column box grid:

    # for vertical direction, # of horiz boxes does not affect subs needed:
    pull_z_start = np.where(radar_z == (min(box_params.z_values)))[0]
    pull_z_end   = np.where(radar_z == (max(box_params.z_values)))[0]
    if pull_z_start.shape[0] != 1 or pull_z_end.shape[0] != 1:
        print('---PROBLEM ASSIGNING z dim SUBSCRIPT(S) FOR PULLING RADAR DATA!!')
        print(f'--- min ht [km]: {min(box_params.z_values)}')
        print(f'--- max ht [km]: {max(box_params.z_values)}')
        print(f'---z array [km]: {radar_z}')
        print('Please check and try again!!!')
        print(f'---Returning without {radarID} data in the column---')
        return column
    # apparantly, these get made as 1-D arrays, but must use scalars as subs,
    # so fix here. have already checked that _start & _end have only one 
    # element (n_zs_match & n_ze_match both are 1 if get to this line)
    # need to add +1 as in python subsetting end index is not inclusive
    pull_z_start = pull_z_start[0]
    pull_z_end   = pull_z_end[0]+1
    #print, 'z subs: ',string(pull_z_start, format='(i0)'),+'  ,  '+string(pull_z_end, format='(i0)')

    if box_params.n_horiz_grid_boxes % 2 == 0:
        #  even # of column grid boxes, so grid box edges line up exactly
        #  column grid data is direct subset in the pyart grid
  
        # get subscripts for center of the pyart grid
        # this point is along edges in the column box grid
        x_origin_sub = np.where(radar_x == 0.0)[0][0]
        y_origin_sub = np.where(radar_x == 0.0)[0][0]
  
        # get # of subscripts before/after
        column_box_horiz_extent_from_cntr = int(box_params.n_horiz_grid_boxes/2)

        # define start & end subscripts for x & y dirs: will use to pull from full grid data
        pull_x_start = x_origin_sub - column_box_horiz_extent_from_cntr
        pull_x_end   = x_origin_sub + column_box_horiz_extent_from_cntr+1
        pull_y_start = y_origin_sub - column_box_horiz_extent_from_cntr
        pull_y_end   = y_origin_sub + column_box_horiz_extent_from_cntr+1
        #print, 'x subs: ',string(pull_x_start, format='(i0)'),+'  ,  '+string(pull_x_end, format='(i0)')
        #print, 'y subs: ',string(pull_y_start, format='(i0)'),+'  ,  '+string(pull_y_end, format='(i0)')
        # similar to above, apparantly the pull_*_start/end subs may get set up here as 1-D arrays, but
        # need scalars to use as subscripts for parsing out column from the bigger data arrays, so fix:
        #if size(pull_x_start, /dimensions) ge 1 then pull_x_start = pull_x_start[0]
        #if size(pull_x_end, /dimensions)   ge 1 then pull_x_end = pull_x_end[0]
        #if size(pull_y_start, /dimensions) ge 1 then pull_y_start = pull_y_start[0]
        #if size(pull_y_end, /dimensions)   ge 1 then pull_y_end = pull_y_end[0]
        #pdb.set_trace()
        # with subscripts for each dir, pull out data in column box grid:
        # dimensions for npol_xx are [time, z, lat, lon] --> col_npol_xx [z, y, x]; time dim is 1
        if 'ZZ' in DS.keys():col_radar_ZZ = radar_ZZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'CZ' in DS.keys():col_radar_CZ = radar_CZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'DR' in DS.keys():col_radar_DR = radar_DR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'RH' in DS.keys():col_radar_RH = radar_RH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'PH' in DS.keys():col_radar_PH = radar_PH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'KD' in DS.keys():col_radar_KD = radar_KD[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'SQ' in DS.keys(): col_radar_SQ = radar_SQ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'SW' in DS.keys():col_radar_SW = radar_SW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'VR' in DS.keys():col_radar_VR = radar_VR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'RR' in DS.keys(): col_radar_RR = radar_RR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'RP' in DS.keys():col_radar_RP = radar_RP[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'RC' in DS.keys():col_radar_RC = radar_RC[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]  
        if 'D0' in DS.keys(): col_radar_D0 = radar_D0[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'NW' in DS.keys():col_radar_NW = radar_NW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'FH' in DS.keys():col_radar_FH = radar_FH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'N2' in DS.keys(): col_radar_N2 = radar_N2[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'DM' in DS.keys():col_radar_DM = radar_DM[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]  
        if 'FZ' in DS.keys(): col_radar_FZ = radar_FZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'MW' in DS.keys(): col_radar_MW = radar_MW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
        if 'MI' in DS.keys(): col_radar_MI = radar_MI[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
    else:
        print('...this section needs to be be worked on...')
        pdb.set_trace()
        # odd # of column grid boxes, so center points will define space differently
        # center point is still center point, but in pyart grid it is an edge,
        # in the column grid, it is the middle of a grid box
  
        # get subscripts for the grid center point
        csub_x = np.where(radar_x == 0.0)[0]
        csub_y = np.where(radar_y == 0.0)[0]
  
        # get important numbers will need for how to get interpolation subscripts
        #n_boxes = self.column_box_params['n_col_horiz_grid_boxes'] #an ODD NUMBER.
        n_edges = box_params.n_horiz_grid_boxes + 1                       #convert from ODD to EVEN NUMBER
        # n_togo_away: # of edges to go away from center pt in each horiz direction
        n_togo_away = n_edges/2                             #EVEN NUMBER
  
        # create arrays to hold subs need to interpolate to
        #  make one arrays each for edges that are before & after center point
        # ASSUMES X & Y EXTENT SAME, ie: # of boxes in X dir = # of boxes in Y dir
        new_x_subs_before = np.zeros(n_togo_away)
        new_x_subs_after  = np.zeros(n_togo_away)
        new_y_subs_before = np.zeros(n_togo_away)
        new_y_subs_after  = np.zeros(n_togo_away)
  
        # for loop thru each point have to go away from the center point
        for n in range(1, n_togo_away+1):
            new_x_subs_after[n-1]  = csub_x + ((n-1)+0.5)
            new_x_subs_before[n-1] = csub_x - ((n-1)+0.5)
            new_y_subs_after[n-1]  = csub_y + ((n-1)+0.5)
            new_y_subs_before[n-1] = csub_y - ((n-1)+0.5)
  
        # combine before & after arrays to one array of subs for interpolation
        temp_x_subs   = [new_x_subs_before, new_x_subs_after]
        sorted_x_subs = np.sort(temp_x_subs)
        new_x_subs    = temp_x_subs[sorted_x_subs]
        temp_y_subs   = [new_y_subs_before, new_y_subs_after]
        sorted_y_subs = np.sort(temp_y_subs)
        new_y_subs    = temp_y_subs[sorted_y_subs]
        #print, 'x subs: ',new_x_subs
        #print, 'y subs: ',new_y_subs
  
        # Now to get the subset of data in the column grid. First, take only the 
        # vertical grid levels needed - altitudes will be on edges so don't need
        # to run interpolation in the vertical direction, only horizontal
        # dimensions for npol_xx are [time, z, lat, lon] --> [z, y, x]; time dim is 1
        hold_ZZ = radar_ZZ[0,pull_z_start:pull_z_end,:,:]
        hold_CZ = radar_CZ[0,pull_z_start:pull_z_end,:,:]
        hold_DR = radar_DR[0,pull_z_start:pull_z_end,:,:]
        hold_RH = radar_RH[0,pull_z_start:pull_z_end,:,:]
        hold_PH = radar_PH[0,pull_z_start:pull_z_end,:,:]
        hold_KD = radar_KD[0,pull_z_start:pull_z_end,:,:]
        if 'SQ' in DS.keys(): hold_SQ = radar_SQ[0,pull_z_start:pull_z_end,:,:]
        hold_SW = radar_SW[0,pull_z_start:pull_z_end,:,:]
        hold_VR = radar_VR[0,pull_z_start:pull_z_end,:,:]
        if 'RR' in DS.keys(): hold_RR = radar_RR[0,pull_z_start:pull_z_end,:,:]
        hold_RP = radar_RP[0,pull_z_start:pull_z_end,:,:]
        hold_RC = radar_RC[0,pull_z_start:pull_z_end,:,:]
        if 'D0' in DS.keys(): hold_D0 = radar_D0[0,pull_z_start:pull_z_end,:,:]
        hold_NW = radar_NW[0,pull_z_start:pull_z_end,:,:]
        hold_FH = radar_FH[0,pull_z_start:pull_z_end,:,:]
        if 'N2' in DS.keys(): hold_N2 = radar_N2[0,pull_z_start:pull_z_end,:,:]
        hold_DM = radar_DM[0,:,:,pull_z_start:pull_z_end]
        if 'FZ' in DS.keys(): hold_FZ = radar_FZ[0,pull_z_start:pull_z_end,:,:]
  
        # run interpolate function to get horiz dim as need for column box grid
        #  will loop thru vertical levels, define col_NPOL_* arrays first
        col_radar_ZZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_CZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_DR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_RH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_PH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_KD = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        if 'SQ' in DS.keys(): col_radar_SQ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_SW = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_VR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        if 'RR' in DS.keys(): col_radar_RR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_RP = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_RC = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))  
        if 'D0' in DS.keys(): col_radar_D0 = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_NW = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_FH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        if 'N2' in DS.keys(): col_radar_N2 = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        col_radar_DM = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
        if 'FZ' in DS.keys(): col_radar_FZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
##############NEED TO EDIT THE LINES BELOW
        #for k=pull_z_start,pull_z_end do begin
        #    col_radar_ZZ[*,*,k] = INTERPOLATE(hold_ZZ[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_CZ[*,*,k] = INTERPOLATE(hold_CZ[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_DR[*,*,k] = INTERPOLATE(hold_DR[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_RH[*,*,k] = INTERPOLATE(hold_RH[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_PH[*,*,k] = INTERPOLATE(hold_PH[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_KD[*,*,k] = INTERPOLATE(hold_KD[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    if test_SQ ne -1 then col_radar_SQ[*,*,k] = INTERPOLATE(hold_SQ[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_SW[*,*,k] = INTERPOLATE(hold_SW[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_VR[*,*,k] = INTERPOLATE(hold_VR[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_RR[*,*,k] = INTERPOLATE(hold_RR[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_RP[*,*,k] = INTERPOLATE(hold_RP[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_RC[*,*,k] = INTERPOLATE(hold_RC[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    if test_D0 ne -1 then col_radar_D0[*,*,k] = INTERPOLATE(hold_D0[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_NW[*,*,k] = INTERPOLATE(hold_NW[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_FH[*,*,k] = INTERPOLATE(hold_FH[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    if test_N2 ne -1 then col_radar_N2[*,*,k] = INTERPOLATE(hold_N2[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    col_radar_DM[*,*,k] = INTERPOLATE(hold_DM[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
        #    if test_FZ ne -1 then col_radar_FZ[*,*,k] = INTERPOLATE(hold_FZ[*,*,k], new_x_subs, new_y_subs, /grid, missing=!VALUES.F_NAN)
#############NEED TO EDIT THE LINES ABOVE
    # If don't have data for SQ, RR, D0, N2, FZ fields, fill with NaNs: same size as other data arrays
    if 'SQ' not in DS.keys():
       col_radar_SQ = col_radar_ZZ*0.0
       col_radar_SQ = sim.bad2nan(col_radar_SQ, 0.0)
    if 'RR' not in DS.keys():
       col_radar_RR = col_radar_ZZ*0.0
       col_radar_RR = sim.bad2nan(col_radar_RR, 0.0)
    if 'D0' not in DS.keys():
       #print('no D0')
       col_radar_D0 = col_radar_ZZ*0.0
       col_radar_D0 = sim.bad2nan(col_radar_D0, 0.0)
    if 'N2' not in DS.keys():
       #print('no N2')
       col_radar_N2 = col_radar_ZZ*0.0
       col_radar_N2 = sim.bad2nan(col_radar_N2, 0.0)
    if 'FZ' not in DS.keys():
       col_radar_FZ = col_radar_ZZ*0.0
       col_radar_FZ = sim.bad2nan(col_radar_FZ, 0.0)
    if 'MW' not in DS.keys():
       col_radar_MW = col_radar_ZZ*0.0
       col_radar_MW = sim.bad2nan(col_radar_MW, 0.0)
    if 'MI' not in DS.keys():
       col_radar_MI = col_radar_ZZ*0.0
       col_radar_MI = sim.bad2nan(col_radar_MI, 0.0)

    #add radar_info to object
    column.add_platform_to_object(radar_info, plat_name=f'lev2_radar_{radarID}')

    #assign each field to the column object
    column.add_variable_to_object(col_radar_ZZ.data, 
                                           var_name=f'lev2_{radarID}_ZZ',
                                           units=radar_Z_units,
                                           long_name=radar_ZZ_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_CZ.data, 
                                           var_name=f'lev2_{radarID}_CZ',
                                           units=radar_Z_units,
                                           long_name=radar_CZ_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_DR.data, 
                                           var_name=f'lev2_{radarID}_DR',
                                           units=radar_DR_units,
                                           long_name=radar_DR_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_RH.data, 
                                           var_name=f'lev2_{radarID}_RH',
                                           units=radar_RH_units,
                                           long_name=radar_RH_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_PH.data, 
                                           var_name=f'lev2_{radarID}_PH',
                                           units=radar_PH_units,
                                           long_name=radar_PH_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_KD.data, 
                                           var_name=f'lev2_{radarID}_KD',
                                           units=radar_KD_units,
                                           long_name=radar_KD_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_SQ.data, 
                                           var_name=f'lev2_{radarID}_SQ',
                                           units=radar_SQ_units,
                                           long_name=radar_SQ_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_SW.data, 
                                           var_name=f'lev2_{radarID}_SW',
                                           units=radar_SW_units,
                                           long_name=radar_SW_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_VR.data, 
                                           var_name=f'lev2_{radarID}_VR',
                                           units=radar_VR_units,
                                           long_name=radar_VR_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_RR.data, 
                                           var_name=f'lev2_{radarID}_RR',
                                           units=radar_R_units,
                                           long_name=radar_RR_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_RP.data, 
                                           var_name=f'lev2_{radarID}_RP',
                                           units=radar_R_units,
                                           long_name=radar_RP_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_RC.data, 
                                           var_name=f'lev2_{radarID}_RC',
                                           units=radar_R_units,
                                           long_name=radar_RC_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_D0.data, 
                                           var_name=f'lev2_{radarID}_D0',
                                           units=radar_D_units,
                                           long_name=radar_D0_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_NW.data, 
                                           var_name=f'lev2_{radarID}_NW',
                                           units=radar_NW_units,
                                           long_name=radar_NW_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_FH.data, 
                                           var_name=f'lev2_{radarID}_FH',
                                           units=radar_FH_units,
                                           long_name=radar_FH_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_N2.data, 
                                           var_name=f'lev2_{radarID}_N2',
                                           units=radar_N2_units,
                                           long_name=radar_N2_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_DM.data, 
                                           var_name=f'lev2_{radarID}_DM',
                                           units=radar_D_units,
                                           long_name=radar_DM_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_FZ.data, 
                                           var_name=f'lev2_{radarID}_FZ',
                                           units=radar_Z_units,
                                           long_name=radar_FZ_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_MW.data, 
                                           var_name=f'lev2_{radarID}_MW',
                                           units=radar_MW_units,
                                           long_name=radar_MW_name,
                                           badval=np.nan)

    column.add_variable_to_object(col_radar_MI.data, 
                                           var_name=f'lev2_{radarID}_MI',
                                           units=radar_MI_units,
                                           long_name=radar_MI_name,
                                           badval=np.nan)

    return column
        
    #self.npol_data_in_column = True
    #pdb.set_trace()
    
def get_grid_from_radar(file, grid_center, vert_spacing, box_spacing, maxh=15000, maxxy=100000):
 
    radar = pyart.io.read(file, file_field_names=True)
    #print(radar.fields.keys())
    #os.system(f'rm {file}')
    
    #computer grid points based on column box resolution
    nz = int(maxh / vert_spacing + 1)
    ny = nx = int(maxxy * 2 / box_spacing + 1)
    
    # perform Cartesian mapping
    grid = pyart.map.grid_from_radars(radar,
            weighting_function='Barnes2', #function to weight for interpolation (check pyart docs)
            grid_origin=grid_center,      #lat and lon of grid origin
            grid_shape=(nz, ny, nx),      #number of points in grid (z, y, x)
            grid_limits=((0., maxh),      #min/max grid location in meters for ((z), (y), (x))
                         (-maxxy, maxxy),
                         (-maxxy, maxxy)))
    return grid
