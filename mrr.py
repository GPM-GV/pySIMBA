import os, sys, glob
from datetime import datetime
import pandas as pd
import numpy as np
import pdb
import math
from scipy import optimize
import tools as sim

def get_mrr_for_column(mrr_files, interval_datetime, box_params, column):

    #Make sure input same number of mrr IDs & mrr files:
    #if len(mrr_IDs) != len(files_mrr):
    #    print(' --- INCORRECT NUMBER OF INPUT ELEMENTS FOR --- ')
    #    print(' ---     mrr_IDs     and    files_mrr_ave   --- ')
    #    print(' --- CHECK INPUT MRR IDs & .ave FILE PATHS! --- ')
    #    print(' --- MRR(s) NOT BEING SET INTO COLUMN GRID  --- ')
    #    return -1
    
    #files_mrr = glob.glob(self.mrr_dir+f'/{self.main_mon}{self.main_day}.ave')
    #mrr_IDs=['mrr2-02'] ##NEED TO CHANGE THIS
    
    #print(self.interval_datetime[0].strftime(f'%Y%m%d%H%M'))
    #print(self.interval_datetime[self.halftime_interval].strftime(f'%Y%m%d%H%M'))
    #print(self.interval_datetime[-1].strftime(f'%Y%m%d%H%M'))
    #print('TESTING ^^^ -- comment out')
    
    # Define holder arrays for info & data fields before loop thru files:
    mrr_lat_d_array = [-9999]
    mrr_lat_m_array = [-9999]
    mrr_lat_s_array = [-9999]
    mrr_lon_d_array = [-9999]
    mrr_lon_m_array = [-9999]
    mrr_lon_s_array = [-9999]
    mrr_elev_array  = [-9999]
    mrr_plat_name_array = ['DUMMY']     #unit name, eg: 'MRR-01'
    mrr_plat_type_array = ['DUMMY']     #plat type, eg: 'MRR'
    mrr_op_mode_array = ['DUMMY']       #operation mode, eg: 'MRR'
    mrr_wavelength_array = [-9999]
    mrr_frequency_array  = [-9999]
    mrr_beam_width_array = [-9999]
    mrr_gate_size_array  = [-9999]
    mrr_plat_timestamp_array = ['DUMMY']
    mrr_offset_array = [-9999]
    mrr_int_width_array = [-9999]
    
    #define arrays with col grid dims, will populate in loop w/ MRR data for unit(s) in the box:
    #  THESE ARE EACH 4-D ARRAY [column z dir  X  column x dir  X  column y dir  X  times in interval ]:
    #col_dims = (n_interval_times, len(self.column_box_params['column_grid_lons']),  len(self.column_box_params['column_grid_lats']),  len(self.column_box_params['column_grid_alts'])  )
    mrr_PIA = np.full(box_params.ground_plat_dims, np.nan)
    mrr_aRef= np.full(box_params.ground_plat_dims, np.nan)
    mrr_ref = np.full(box_params.ground_plat_dims, np.nan)
    mrr_RR  = np.full(box_params.ground_plat_dims, np.nan)
    mrr_LWC = np.full(box_params.ground_plat_dims, np.nan)
    mrr_VEL = np.full(box_params.ground_plat_dims, np.nan)
    mrr_disdro_Dm = np.full(box_params.ground_plat_dims, np.nan)
    mrr_disdro_Nw = np.full(box_params.ground_plat_dims, np.nan)
    mrr_disdro_Z  = np.full(box_params.ground_plat_dims, np.nan)
    mrr_ref_ku_adj = np.full(box_params.ground_plat_dims, np.nan)
    mrr_data_quality = np.full(box_params.ground_plat_dims, np.nan)
    
    #define variable names and units
    mrr_PIA_name = 'two-way integrated attenuation'
    mrr_aRef_name = 'attenuated reflectivity'
    mrr_ref_name = 'reflectivity corrected for attenuation due to raindrops'
    mrr_refKu_name = 'reflectivity converted to ku using observed DSDs'
    mrr_RR_name = 'rainrate'
    mrr_LWC_name = 'liquid water content'
    mrr_WVEL_name = 'velocity: capital W'
    mrr_disdro_Dm_name = 'mass-weighted mean diameter computed from DSD profile'
    mrr_disdro_Nw_name = 'normalized intercept parameter computed from DSD profile'
    mrr_disdro_Z_name = 'radar reflectivity computed from DSD profile'
    mrr_data_quality_name = '% of valid spectra during the .ave interval'
    mrr_units_PIA = 'dB'
    mrr_units_ref = 'dBZ'
    mrr_units_RR = 'mm h^-1'
    mrr_units_LWC = 'g m^-3'
    mrr_units_WVEL = 'm s^-1'
    mrr_units_disdro_Dm = 'mm'
    mrr_units_disdro_Nw = 'log(Nw)'
    mrr_units_disdro_Z = 'dBZ'
    mrr_units_data_quality = '%'
    
    # Will need to loop through input list of MRR units:
    #for i, mrr_id in enumerate(mrr_IDs):
    for i, file_mrr in enumerate(mrr_files):
        #file_mrr = files_mrr[i]
        mrr_id = os.path.basename(file_mrr)[0:7]
    
        # Set initial mrr_plat_info values:
        #pdb.set_trace()
        mrr_plat_info  = sim.set_plat_values(mrr_id.upper(), file_mrr)
        unit_elevation = mrr_plat_info['elev']
    
        #get subscript where MRR is closest to a column grid point:
        mrr_loc_decdeg = sim.dms2dd(mrr_plat_info['lat_d'], mrr_plat_info['lat_m'], mrr_plat_info['lat_s'], 
                                mrr_plat_info['lon_d'], mrr_plat_info['lon_m'], mrr_plat_info['lon_s'])
        mrr_lat = mrr_loc_decdeg[0]
        mrr_lon = mrr_loc_decdeg[1]
        lat_ok = 'N'
        lon_ok = 'N'
        if mrr_lat <= box_params.box_max_lat and mrr_lat >= box_params.box_min_lat: lat_ok = 'Y'
        if mrr_lon <= box_params.box_max_lon and mrr_lon >= box_params.box_min_lon: lon_ok = 'Y'
        if lat_ok == 'N' and lon_ok == 'N':
            print('------------------------------------------------')
            print('-------- before reading .ave file --------------')
            print('---APPEARS THIS MRR IS NOT IN THE COLUMN BOX!---')
            print(f'-------plat_name: {mrr_plat_info["plat_name"]}')
            print('-------NOTICE: MRR data NOT SET IN column grid!')
            mrr_info = np.nan
            continue # so will only continue from here & read in mrr file if location is OK
    
        # Call read_data_mrr Function to read in MRR .ave data:
        data = read_mrr2(file_mrr)
        #  data structure returned has keys:
        #  timestamp     avg_interval_sec  gate_size_m
        #  samp_rate_hz  radar_alt_mASL    gate_hts_m
        #  data_quality  mrrD (00-63)      mrrN (00-63) both are [64, n_times, n_gates]
        #  PIA           attenRef          ref
        #  rainRate      LWC               velocity
        mrr_datestr = data['timestamp'] #format is YYMMDDHHMMSS  200104000002
        mrr_datetime = sim.datestr2datetime(mrr_datestr)
        
        #mrr_yr = strmid(data.timestamp, 0, 2)		#these are all arrays:
        #mrr_mon= strmid(data.timestamp, 2, 2)
        ##mrr_day= strmid(data.timestamp, 4, 2)
        #mrr_hr = strmid(data.timestamp, 6, 2)
        #mrr_min= strmid(data.timestamp, 8, 2)
        #mrr_sec= strmid(data.timestamp, 10, 2)
        mrr_plat_info['gate_size'] = data['gate_size_m']

        # get mrr gate hts info that will be same for all times from this unit:
        mrr_gate_heights = data['gate_hts_m']*1.0
        mrr_heights      = mrr_gate_heights + unit_elevation

        # Look for MRR values for ea min in the time interval:
        mrr_unit_in_column_test_cnt = 0
        for t, this_date in enumerate(interval_datetime):
            time_sub = np.where(mrr_datetime == this_date)[0]
            if len(time_sub) < 1:
                continue # this time is not available in mrr data
            elif len(time_sub) > 1:
                #found > ONE MRR time corresponding to this interval time:
                #this very highly unlikely as the MRR .ave files (SIMBA input) typically use 60 s averaging interval - will set
                #here to use 1st element that matches the interval time in case it occurs, but shouldn't really ever get here:
                print(f'HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")}')
                print(f' in ave file: {os.path.basename(file_mrr)}')
                print('---Returning without MRR data in the column---')
                return column

            #Location is OK, and current interval time is OK - so proceed:
            #  Increment counter, set attributes
            mrr_unit_in_column_test_cnt = mrr_unit_in_column_test_cnt +1
            if mrr_unit_in_column_test_cnt == 1:
                print('    MRR unit located in column box AND available during time interval')
                mrr_lat_d_array.append(mrr_plat_info['lat_d'])
                mrr_lat_m_array.append(mrr_plat_info['lat_m'])
                mrr_lat_s_array.append(mrr_plat_info['lat_s'])
                mrr_lon_d_array.append(mrr_plat_info['lon_d'])
                mrr_lon_m_array.append(mrr_plat_info['lon_m'])
                mrr_lon_s_array.append(mrr_plat_info['lon_s'])
                mrr_elev_array.append(mrr_plat_info['elev'])
                mrr_plat_name_array.append(mrr_plat_info['plat_name'])
                mrr_plat_type_array.append(mrr_plat_info['plat_type'])
                mrr_op_mode_array.append(mrr_plat_info['operation_mode'])
                mrr_wavelength_array.append(mrr_plat_info['wavelength'])
                mrr_frequency_array.append(mrr_plat_info['frequency'])
                mrr_beam_width_array.append(mrr_plat_info['beam_width'])
                mrr_gate_size_array.append(mrr_plat_info['gate_size'])
                mrr_int_width_array.append(box_params.halftime_interval)
                # just will set main_plat time w/ sec as '00' since using an interval of data
                #mrr_plat_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                mrr_plat_timestamp_array.append(column.time)
                mrr_offset_array.append(0)

            # Find the (horizontal) column grid point closest to the MRR -> place MRR values at this point
            # -> arrays of the MRR data currently set up as 1-D arrays: [ht] in the column box
            # -> arrays that will set data into were defined above this loop & have same dims as column box grid + time [time interval, lon, lat, elev]
            # NOTE: will maintain original MRR location in plat_info for attributes
            closest_col_lat_sub = sim.closest(box_params.lat_values, mrr_lat)[0]
            closest_col_lon_sub = sim.closest(box_params.lon_values, mrr_lon)[0]
            #print, 'closest grid lat, lon: ',column_box_params.column_grid_lats[closest_col_lat_sub], column_box_params.column_grid_lons[closest_col_lon_sub]
            #print(f'closest grid lat, lon: {self.column_box_params["column_grid_lats"][closest_col_lat_sub]} , {self.column_box_params["column_grid_lons"][closest_col_lon_sub]}')

            # Get the MRR data for this time:  remember data['mrrD'] and data['mrrN'] are in [64, n_times, n_gates] format
            time_sub = time_sub[0] #time_sub is now an integer; before numpy where returned an array
            mrrD = data['mrrD'][:,time_sub,:] #should be a 2D array [64, n_gates]
            mrrN = data['mrrN'][:,time_sub,:]
            
            # These are in [n_times, n_gates]
            PIA        = data['PIA'][time_sub,:]
            attenRef   = data['attenRef'][time_sub,:]
            ref        = data['ref'][time_sub,:]
            rainRate   = data['rainRate'][time_sub,:]
            LWC        = data['LWC'][time_sub,:]
            velocity   = data['velocity'][time_sub,:]
            quality    = data['data_quality'][time_sub]
    
            # Prepare for input to disdro code/module from Patrick Gatlin to compute Dm, Nw, Z profile
            #  these input arrays need to be set up as [# of bins x # of hts]  (assumes always have 64 drop size bins)
            N_of_D_array = mrrN
            D_bins_array = mrrD
            #pdb.set_trace()

            #Call disdro_mrr_get_dsd_params.py:  use /no_fit bc only am after Dm, Nw, & Z             
            disdro_result = get_dsd_params(N_of_D_array, D_bins_array, mrr_gate_heights, no_fit=True)
            disdro_Dm = disdro_result['dm']             #[mm]
            disdro_Nw = disdro_result['nw']             #just Nw (NOT log(Nw))
            disdro_Z  = disdro_result['reflectivity']   #[dBZ]
            #pdb.set_trace()
            
            # Interpolate in the vertical to match column grid levels:  for profile fields from .ave file
            # use scipy interpolate function to get val from MRR hts to column grid hts:
            col_PIA         = sim.interp(PIA, mrr_heights, box_params.z_values)
            col_rainRate    = sim.interp(rainRate, mrr_heights, box_params.z_values)
            col_LWC         = sim.interp(LWC, mrr_heights, box_params.z_values)
            col_velocity    = sim.interp(velocity, mrr_heights, box_params.z_values)
            #pdb.set_trace()
            # For reflectivity & Nw:  get INTERPOL results for LINEAR UNITS, then put into LOG UNITS after
            orig_lin_aRef    = 10.0**((attenRef)/10.0)
            orig_lin_ref     = 10.0**((ref)/10.0)
            orig_lin_disdro_Z= 10.0**((disdro_Z)/10.0)
            col_attenRef  = sim.interp(orig_lin_aRef, mrr_heights, box_params.z_values)
            col_ref       = sim.interp(orig_lin_ref, mrr_heights, box_params.z_values)
            col_disdro_Z  = sim.interp(orig_lin_disdro_Z, mrr_heights, box_params.z_values)
            col_disdro_Nw = sim.interp(disdro_Nw, mrr_heights, box_params.z_values)
            col_disdro_Dm = sim.interp(disdro_Dm, mrr_heights, box_params.z_values)
            #pdb.set_trace()
            
            # put reflectivities & Nw into log units:
            col_disdro_Z  = 10*np.log10(col_disdro_Z) #[dBZ]
            col_disdro_Nw = np.log10(col_disdro_Nw)   #[log(Nw)]
            col_ref       = 10*np.log10(col_ref)      #[dBZ]
            col_ref_ku_adj = 0.974*col_ref - 0.097    #[dBZ] - derived using Wallops 2DVD DSD and radar model following Liao et al. 2020 method
            col_attenRef  = 10*np.log10(col_attenRef) #[dBZ]
            #pdb.set_trace()
     
            # Ensure have NaNs where should be:
            low_hts  = np.where(box_params.z_values < min(mrr_heights))[0]
            high_hts = np.where(box_params.z_values > max(mrr_heights))[0]
            if len(low_hts) > 0:
                col_PIA[low_hts]        = np.nan
                col_attenRef[low_hts]   = np.nan
                col_ref[low_hts]        = np.nan
                col_ref_ku_adj[low_hts] = np.nan
                col_rainRate[low_hts]   = np.nan
                col_LWC[low_hts]        = np.nan
                col_velocity[low_hts]   = np.nan
                col_disdro_Z[low_hts]   = np.nan
                col_disdro_Dm[low_hts]  = np.nan
                col_disdro_Nw[low_hts]  = np.nan

            if len(high_hts) > 0:
                col_PIA[high_hts]        = np.nan
                col_attenRef[high_hts]   = np.nan
                col_ref[high_hts]        = np.nan
                col_ref_ku_adj[high_hts] = np.nan
                col_rainRate[high_hts]   = np.nan
                col_LWC[high_hts]        = np.nan
                col_velocity[high_hts]   = np.nan
                col_disdro_Z[high_hts]   = np.nan
                col_disdro_Dm[high_hts]  = np.nan
                col_disdro_Nw[high_hts]  = np.nan
            #Set values into the column grid arrays for THIS time for THIS unit:
            mrr_PIA[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_PIA
            mrr_aRef[:, closest_col_lat_sub,closest_col_lon_sub,t]        = col_attenRef
            mrr_ref[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_ref
            mrr_RR[:, closest_col_lat_sub,closest_col_lon_sub,t]          = col_rainRate
            mrr_LWC[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_LWC
            mrr_VEL[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_velocity
            mrr_disdro_Z[:, closest_col_lat_sub,closest_col_lon_sub,t]    = col_disdro_Z
            mrr_ref_ku_adj[:, closest_col_lat_sub,closest_col_lon_sub,t]  = col_ref_ku_adj
            mrr_disdro_Dm[:, closest_col_lat_sub,closest_col_lon_sub,t]   = col_disdro_Dm
            mrr_disdro_Nw[:, closest_col_lat_sub,closest_col_lon_sub,t]   = col_disdro_Nw
            mrr_data_quality[:, closest_col_lat_sub,closest_col_lon_sub, 0] = quality

    # Remove dummy/placeholder elements of info arrays:
    if len(mrr_lat_d_array) > 1:
        mrr_lat_d_array         = mrr_lat_d_array[1:]
        mrr_lat_m_array         = mrr_lat_m_array[1:]
        mrr_lat_s_array         = mrr_lat_s_array[1:]
        mrr_lon_d_array         = mrr_lon_d_array[1:]
        mrr_lon_m_array         = mrr_lon_m_array[1:]
        mrr_lon_s_array         = mrr_lon_s_array[1:]
        mrr_elev_array          = mrr_elev_array[1:]
        mrr_plat_name_array     = mrr_plat_name_array[1:]
        mrr_plat_type_array     = mrr_plat_type_array[1:]
        mrr_op_mode_array       = mrr_op_mode_array[1:]
        mrr_wavelength_array    = mrr_wavelength_array[1:]
        mrr_frequency_array     = mrr_frequency_array[1:]
        mrr_beam_width_array    = mrr_beam_width_array[1:]
        mrr_gate_size_array     = mrr_gate_size_array[1:]
        mrr_plat_timestamp_array = mrr_plat_timestamp_array[1:]
        mrr_offset_array         = mrr_offset_array[1:]
        mrr_int_width_array      = mrr_int_width_array[1:]
    else:
        #don't have any MRR units in the column box and/or main_plat time is unavailable
        print('--- There are no MRR units in the column box and/or main plat time is unavailable---')
        print('---Returning without MRR data in the column---')
        return column
    
    
    #Set up dictionary
    mrr_info = {'plat_name':mrr_plat_name_array,
                'lat_d':mrr_lat_d_array, 'lat_m':mrr_lat_m_array, 'lat_s':mrr_lat_s_array,
                'lon_d':mrr_lon_d_array, 'lon_m':mrr_lon_m_array, 'lon_s':mrr_lon_s_array,
                'elev':mrr_elev_array, 
                'plat_type':mrr_plat_type_array,
                'operation_mode':mrr_op_mode_array,
                'wavelength':mrr_wavelength_array,
                'frequency':mrr_frequency_array,
                'beam_width':mrr_beam_width_array,
                'gate_size':mrr_gate_size_array,
                'timestamp':mrr_plat_timestamp_array,
                'offset_vs_main':mrr_offset_array,
                'time_interval_width':mrr_int_width_array }

    #add mrr_info to object
    column.add_platform_to_object(mrr_info, plat_name='mrr')

    #assign each field to column object
    column.add_variable_to_object(mrr_PIA, 
                                  var_name='mrr_PIA', 
                                  units=mrr_units_PIA,
                                  long_name=mrr_PIA_name)

    column.add_variable_to_object(mrr_aRef, 
                                  var_name='mrr_aRef', 
                                  units=mrr_units_ref,
                                  long_name=mrr_aRef_name)

    column.add_variable_to_object(mrr_ref, 
                                  var_name='mrr_ref', 
                                  units=mrr_units_ref,
                                  long_name=mrr_ref_name)

    column.add_variable_to_object(mrr_ref_ku_adj, 
                                  var_name='mrr_refKu', 
                                  units=mrr_units_ref,
                                  long_name=mrr_refKu_name)

    column.add_variable_to_object(mrr_RR, 
                                  var_name='mrr_RR', 
                                  units=mrr_units_RR,
                                  long_name=mrr_RR_name)

    column.add_variable_to_object(mrr_LWC, 
                                  var_name='mrr_LWC', 
                                  units=mrr_units_LWC,
                                  long_name=mrr_LWC_name)

    column.add_variable_to_object(mrr_VEL, 
                                  var_name='mrr_WVEL', 
                                  units=mrr_units_WVEL,
                                  long_name=mrr_WVEL_name)

    column.add_variable_to_object(mrr_disdro_Dm, 
                                  var_name='mrr_disdro_dm', 
                                  units=mrr_units_disdro_Dm,
                                  long_name=mrr_disdro_Dm_name)

    column.add_variable_to_object(mrr_disdro_Nw, 
                                  var_name='mrr_disdro_Nw', 
                                  units=mrr_units_disdro_Nw,
                                  long_name=mrr_disdro_Nw_name)

    column.add_variable_to_object(mrr_disdro_Z, 
                                  var_name='mrr_disdro_Z', 
                                  units=mrr_units_ref,
                                  long_name=mrr_disdro_Z_name)

    column.add_variable_to_object(mrr_data_quality, 
                                  var_name='mrr_data_quality', 
                                  units=mrr_units_data_quality,
                                  long_name=mrr_data_quality_name)

    return column
    #self.mrr_data_in_column = True

'''
-----------------------------------------------------
;
; Function to read in data from MRR .ave files
;
; Example .ave file first two rows
    MRR 200104000002 UTC AVE    60 STP    35 ASL     0 SMP 125e3 SVS 6.0.0.4 DVS 6.10 DSN 0503052419 CC 1745637 MDQ 100 TYP AVE
    H       35     70    105    140    175    210    245    280    315    350    385    420    455    490    525    560    595    630    665    700    735    770    805    840    875    910    945    980   1015   1050   1085
    
  Header lines have 23 field/columns
  data lines have 32 field/columns
  
  MRR, timestamp are strings
  AVE, STP, ASL, SMP, MDQ all ints
  SMP is float
  
  use_for_names = ['mrr','timestamp','timezone','ave','avg_int_sec','stp','gate_size_m','asl','radar_alt_mASL','smp','samp_rate_hz','svs','sftwr_version','dvs','hrdwr_version','dsn','serial_num','cc','calib_const','mdq','data_quality','typ','mrr_type',$
  'gate_hts_m','gate_01','gate_02','gate_03','gate_04','gate_05','gate_06','gate_07','gate_08','gate_09','gate_10','gate_11','gate_12','gate_13','gate_14','gate_15','gate_16','gate_17','gate_18','gate_19','gate_20','gate_21','gate_22','gate_23','gate_24','gate_25','gate_26','gate_27','gate_28','gate_29','gate_30','gate_31',$
  
;-----------------------------------------------------
'''
def read_mrr2(file_mrr):

    # read .ave (text) file using pandas read_fwf function (https://pandas.pydata.org/docs/reference/api/pandas.read_fwf.html)
    df = pd.read_fwf(file_mrr, widths=[3,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7], header=None)

    # (1) extract # of gates and times
    # assumes have same # of gates (31) at every time in the file
    n_gates = df.iloc[0].shape[0] - 1 #originally hard coded as n_gates = 31
    n_times = df[df[0] == 'MRR'].shape[0]

    # (2) find all rows with the field of interest (PIA, z, RR, etc)
    # then drop the column with the str name MRR so can convert to data type float from object
    # these are all in n_times x n_gates (e.g., 1440, 31)
    
    mrr_gate_hts_m = df.iloc[1][1:].to_numpy(dtype=int)
    
    pia = df[df[0] == 'PIA']
    mrr_PIA = pia.drop(pia.columns[[0]], axis=1).to_numpy(dtype=float)
    
    zz = df[df[0] == 'z']
    mrr_attenRef = zz.drop(zz.columns[[0]], axis=1).to_numpy(dtype=float)
    
    ZZ = df[df[0] == 'Z']
    mrr_Ref = ZZ.drop(ZZ.columns[[0]], axis=1).to_numpy(dtype=float)
    
    RR = df[df[0] == 'RR']
    mrr_rainRate = RR.drop(RR.columns[[0]], axis=1).to_numpy(dtype=float)
    
    LWC = df[df[0] == 'LWC']
    mrr_LWC = LWC.drop(zz.columns[[0]], axis=1).to_numpy(dtype=float)
    
    W = df[df[0] == 'W']
    mrr_velocity = W.drop(W.columns[[0]], axis=1).to_numpy(dtype=float)
    
    df2 = df[df[0] == 'MRR']
    date = df2[1].to_numpy(dtype=int).astype('str')
    time = df2[2].to_numpy(dtype=int).astype('str')
    
    mrr_avg_int   = df2[4].to_numpy(dtype=int)
    mrr_gate_size_ = df2[6].to_numpy(dtype=str)
    mrr_radar_alt_ = df2[7].to_numpy()
    mrr_samp_rate_ = df2[8].to_numpy()
    mrr_data_qual_ = df2[16].to_numpy()
    
    #extract DSD for 64 bins
    mrr_D = np.full((64, n_times, n_gates), np.nan)
    mrr_N = np.full((64, n_times, n_gates), np.nan)
    for j in np.arange(0,64):
        D_ = df[df[0] == 'D'+str(j).zfill(2)]
        mrr_D[j] = D_.drop(D_.columns[[0]], axis=1).to_numpy(dtype=float)
    
        N_ = df[df[0] == 'N'+str(j).zfill(2)]
        mrr_N[j] = N_.drop(N_.columns[[0]], axis=1).to_numpy(dtype=float)
    
    #arrays need to be modified due to incorrect reading with pandas
    mrr_timestamp = np.array([])
    mrr_gate_size = np.array([])
    mrr_radar_alt = np.array([])
    mrr_samp_rate = np.array([])
    mrr_data_qual = np.array([])
    for i in range(n_times):
        mrr_timestamp = np.append(mrr_timestamp, date[i]+time[i].zfill(6))
        mrr_gate_size = np.append(mrr_gate_size, int(mrr_gate_size_[i][0:2]))
        mrr_radar_alt = np.append(mrr_radar_alt, int(mrr_radar_alt_[i][0]))
        mrr_samp_rate = np.append(mrr_samp_rate, float(mrr_samp_rate_[i][3:6])*1000)
        mrr_data_qual = np.append(mrr_data_qual, int(mrr_data_qual_[i][4:]))
    
    #(3) create dictionary to combine all fields:
    data_mrr = {'timestamp':mrr_timestamp,
                'avg_interval_sec':mrr_avg_int[0],
                'gate_size_m':mrr_gate_size[0],
                'radar_alt_mASL':mrr_radar_alt[0],
                'samp_rate_hz':mrr_samp_rate[0],
                'data_quality':mrr_data_qual,
                'gate_hts_m':mrr_gate_hts_m,
                'mrrD':mrr_D,
                'mrrN':mrr_N,
                'PIA':mrr_PIA,
                'attenRef':mrr_attenRef,
                'ref':mrr_Ref,
                'rainRate':mrr_rainRate,
                'LWC':mrr_LWC,
                'velocity':mrr_velocity  }
    #pdb.set_trace()
    return data_mrr


'''
;-----------------------------------------------------
;PURPOSE:
;   Get the dsd model parameters from a single MRR N(D) profile
;INPUTS:
;   N(D) for each diameter bin at each height [mm^-1m^-3] ;2-d array [bin,height]
;   Diameter of each bin at each height [mm]	;2-d array [bin,height]
;   Height of each measurement [m] 	;1-d array  [height]
;OUTPUT:
;   A structure containing the normalized gamma DSD parameters (Testud et al. 2001) at each height
;
;DEPENDENCIES:
;   disdro_fit_dsd_norm_gamma.pro
;
;EXAMPLE USAGE:
;   dsd_params=disdro_mrr_get_dsd_params2(nd[*,*],diameter[*,*],height[*],/quiet)
;
;HISTORY:
;   Written by P. Gatlin (NASA) May 2017
;   Minor update:  added conditional so that if /NO_FIT is set,
;           value of mu is not added to the structure
;           (S. Wingo, Aug 2017)
;   Update: code converted to python by Charanjit S. Pabla (NASA/SSAI) August 2022
;-----------------------------------------------------
'''

def get_dsd_params(nd, diameter, height, no_fit=True):

    #initialize arrays to store data in dictionary
    #nd       [bins, ngates] > [64, 31]
    #diameter [bins, ngates] > [64, 31]
    #height   [ngates]       > [31]
    nhts=len(height)
    n_bins=diameter.shape[0]
    params={'moments':     np.full((8,nhts), np.nan),
            'dm':          np.full(nhts, np.nan),
            'nw':          np.full(nhts, np.nan),
            'sigma_dm':    np.full(nhts, np.nan),
            'mu':          np.full(nhts, np.nan),
            'reflectivity':np.full(nhts, np.nan)
            }

 
    ht=height
    min_ht=ht[1] #ht[1] ;minimum ht of 1st good gate
    min_diameter = 0.246
    max_diameter = 5.03  #diameter range of valid diameter retrievals (from MRR Physical basis manual)


    moments=np.full((8,len(ht)), np.nan)
    #for h=0,n_elements(height)-2 do begin
    for h, val in enumerate(height):
    
        #break when len(height) - 2 is reached
        if h == nhts - 1: break
    
        #skip heights below the minimum
        if height[h] < min_ht:
            continue
        #y=[height[h], height[h+1], height[h+1], height[h], height[h]] #this is not being used so ignore - Pabla, Aug 2022
        
        l=np.where(diameter[:,h] > 0)[0]
        n_bins=len(l) #n_elements(diameter[*,0])
        
        #deltaD=(diameter[l[2:n_bins-1],h]-diameter[l[0:n_bins-3],h])/2.0  #from MRR Physical basis manual Section 3.2
        
        #python version has different indexing vs IDL above
        #calculate deltaD
        deltaD=(diameter[l[2:n_bins],h]-diameter[l[0:n_bins-2],h])/2.0
        deltaD = np.insert(deltaD, 0, diameter[l[1],h] - diameter[l[0],h])
        deltaD = np.append(deltaD, diameter[l[n_bins-1],h] - diameter[l[n_bins-2],h])
        
        #deltaD=[diameter[l[1],h]-diameter[l[0],h],deltaD,diameter[l[n_bins-1],h]-diameter[l[n_bins-2],h]]
        deltaD_low=deltaD

        deltaD =(diameter[l[2:n_bins],h+1]-diameter[l[0:n_bins-2],h+1])/2.  #from MRR Physical basis manual Section 3.2
        deltaD = np.insert(deltaD, 0, diameter[l[1],h+1] - diameter[l[0],h+1])
        deltaD = np.append(deltaD, diameter[l[n_bins-1],h+1] - diameter[l[n_bins-2],h+1])
        
        #deltaD=[diameter[l[1],h+1]-diameter[l[0],h+1],deltaD,diameter[l[n_bins-1],h+1]-diameter[l[n_bins-2],h+1]]
        deltaD_high=deltaD

#      dbinsz_low=(0.1887d)/(6.18*exp(0-0.6*diameter[d,h])*(1+3.68e-5*height[h]+1.71e-9*height[h]^2)) ;derived from MRR Physical basis manual eqt'n 2.5 delta_v/delta_d where delta_v=0.1887m/s
#      dbinsz_high=(0.1887d)/(6.18*exp(0-0.6*diameter[d,h+1])*(1+3.68e-5*height[h+1]+1.71e-9*height[h+1]^2))     
#   dbinsz_low=(0.1887d)/(6.18*exp(0-0.6*diameter[l,h])*(1+3.68e-5*height[h]+1.71e-9*height[h]^2))
        
        #->limit data to bins within min/max diameter and with valid N(D)s
        x = np.where( (diameter[l,h] >= min_diameter) & (diameter[l,h] <= max_diameter) & (nd[l,h] >= 0) & (np.isfinite(nd[l,h])) )[0]
        if len(x) == 0:
            continue
        for p in range(8):
            moments[p,h] = np.sum(diameter[l[x],h]**p * nd[l[x],h] * deltaD_low[x] )
        
        n=l[x]

        # check to see if "n" variable elements are within bounds - when converting code from IDL
        # there was a bug and in order to have a temporary fix in python, we will use the following method
        # replace indices that are beyond the size of array n to the last element - Pabla, August 2022
        ind = np.where(n >= len(deltaD_low))[0]
        n[ind] = len(deltaD_low) - 1
        n_ = l[x] # have to index this into nd and diameter only while deltaD_low needs to be indexed by n

        #if h == 23: pdb.set_trace()
        fit = disdro_fit_dsd_norm_gamma(nd[n_,h], diameter[n_,h], deltaD_low[n], dmin=min_diameter, dmax=max_diameter, no_fit=no_fit)
        #if h == 4: pdb.set_trace()
        if not no_fit: 
            params['sigma_dm'][h] = 1.0/np.sqrt(4 + fit['mu'])
            params['mu'][h]       = fit.mu
        params['nw'][h] = fit['nw']     
        params['dm'][h] = fit['dm']

    #-->get reflectivity
    l=np.where(moments[6,:] > 0)[0]  
    if len(l) > 0: params['reflectivity'][l] = 10*np.log10(moments[6,l]) ###issue with REF -- check ND Aug 17, 2022
    #pdb.set_trace()
    return params

# ***************************************************************************************

'''
-----------------------------------------------------
    Calculate the DSD moments

-----------------------------------------------------
'''
#======================================
#Calculate moments from DSD spectrum
#======================================
def get_moments_observed(n_d, diameter_bins, bin_size, dmin=-1, dmax=-1):
    #n_d = N(D) array ( mm^-1m^-3 )
    #diameter_bins = center of each diameter bin (mm)
    #bin_size = bin size of each diameter bin (mm)
    #pdb.set_trace()
    x = np.where(n_d > 0)[0]
    if len(x) > 0:
        moments = np.full(8, np.nan)
        if dmin > -1:
            d1 = np.where( ( (diameter_bins[x] - bin_size[x]) / 2.0) >= dmin)[0]
            if len(d1) > 0:
                x = x[d1]
        if dmax > -1:
            d2 = np.where( ( (diameter_bins[x] + bin_size[x]) / 2.0) <= dmax)[0]
            if len(d2) > 0:
                x = x[d2]   
    for p in range(8):
        moments[p] = np.sum(diameter_bins[x]**p * n_d[x]* bin_size[x])   
    return moments

#======================================
#Least Squares method for Normalized gamma DSD (Testud et al. 2001, JAM)
#======================================
def fit_gamma_norm_mu(mu, dm, nw, dbin, nd_obs): #,dbin=dbin,nd_obs=nd_obs
    #mu = a[0]
    #dm=a[1]
    #nw=a[2]
    
    term1 = (math.gamma(4.0) / (4.0**4.0))
    term2 = (((4+mu)**(4+mu)) / math.gamma(4.0+mu))
    term3 = ((dbin/dm)**mu)
    term4 = np.exp(-1*(4.0+mu) * (dbin/dm))
    
    f_mu     = term1 * term2 * term3 * term4 #EQ 16 in Testud et al. 2001
    nd_model = nw * f_mu #EQ 17 in Testud et al. 2001
    #nd_model = nw * f_mu * ((dbin/dm)^mu) * np.exp(-1*(4.0+mu)) * (dbin/dm)

    return np.sum(abs(np.log10(nd_obs) - alog10(nd_model)))


'''
;###############################################################################
;Function to find model paramers of normalized gamma DSD via the least squares/mu-optimization method (Testud et al. 2001,JAM)
;Input:
;   N(D)=number concentration in each diameter bin (mm^-1m^-3)
;   diameter_bins=CENTER diameter of each bin (mm)
;   bin_size=width of each bin (mm)
;Outputs:
;   structure containing the parameters: NW,mu,Dm of N(D)=NW*f_mu*(D/Dm)^mu*exp(0-(4+mu)*D/Dm), [0<D<Inf];
;   a return value of -1 indicates fitting was not successful
;
;Keywords:
;   dmin=mininum diameter of truncation (default=0)
;   dmax=maximum diameter of truncation  [default=max([3*dm,max(diameter+bin_size/2d)])]
;   min_bins=minimum number of bins required for fit (default=3)
;   [/quiet] - set this keyword to not print out any status messages
;   [/NO_FIT] - set this keyword to NOT peform fitting required to get mu
;           instead Nw,Dm,Moments and LWC are returned from observed spectra
;
;Example:
;   fit_params=disdro_fit_dsd_tmom(dsd[i].bin.nd,dsd[i].bin.diameter,dsd[i].bin.width,/quiet)
;
;Minimization performed via the IDL's built-in AMOEBA function which uses a downhill simplex method
;
;History:
;   P.Gatlin (NASA/MSFC) Jan 2017
;###############################################################################
'''

def disdro_fit_dsd_norm_gamma(n_d, diameter_bins, bin_size, dmin=-1, dmax=-1,
                              min_bins=0, no_fit=False):

    #-->initialize variables
    x = diameter_bins.astype('float64')
    y = n_d.astype('float64') 
    w = np.where(y > 0)[0]  #make sure N(D)>0
    if min_bins == 0: min_bins=3
    if len(w) < min_bins: return -1
 
    #-->remove bins where N(D)=0
    l = np.where(y > 0)[0]

    #-->find min and max diameters based on bins and their half-widths
    if dmin == -1: dmin = min(x[l] - bin_size[l] / 2.0)
    if dmax == -1: dmax = max(x[l] + bin_size[l] / 2.0)
 
    #-->remove diameters outside min and max
    d = np.where( ((x[l] - bin_size[l] / 2.0) >= dmin) & ((x[l] + bin_size[l] / 2.0) <= dmax) )[0]
    if len(d) > 0: l = l[d]

    #-->assume errors in measured N(D) follow are Poisson
    y_errors = np.sqrt(abs(y[l]))

    #-->make sure we have more data points than model's free parameters
    if len(y[l]) < 1+1: return -1
    if len(y[l]) < min_bins: return -2

    #=============>FITTING COMPUTATIONS<===============  
    #pdb.set_trace()
    moments_obs = get_moments_observed(n_d[l], diameter_bins[l], bin_size[l])
    #pdb.set_trace()

    #-estimate intitial value of shape and slope from untruncated MOM  
    nu_untrunc = (moments_obs[4]**2) / (moments_obs[2] * moments_obs[6])
    mu_mom = ((7-11*nu_untrunc) - (((7-11*nu_untrunc)**2) - 4 * (nu_untrunc-1) * (30*nu_untrunc-12))**(1.0/2.0))*1.0 / (2*(nu_untrunc-1))
    lambda_mom = ((4+mu_mom) * (3+mu_mom) * moments_obs[2] / moments_obs[4])**(1.0/2.0)
    n0_mom = moments_obs[3] * (lambda_mom**(4+mu_mom)) / math.gamma(4+mu_mom)
    dm = moments_obs[4] / moments_obs[3]
 
    #-->NORMALIZED GAMMA FIT LEAST SQUARES/MU-OPTIMIZATION METHOD (Testud et al. 2001, JAM)
    parinfo = np.repeat({'value':0.0, 'fixed':0, 'limited':[0,0], 'limits':[0.0,0.0],'step':0, 'tnside':0},3)
    parinfo[0]['value']  = mu_mom
    parinfo[0]['tnside'] = 2 # & parinfo[0].limited=[1,1] & parinfo[0].limits=[-1,15]
    parinfo[1]['value']  = dm
    parinfo[1]['fixed']  = 1
    rhow = 1.0 #density of liquid water in g/cm^3
    lwc = np.pi / 6.0 * rhow * (1e-3) * moments_obs[3]  #g/m3
    nw= (4.0**4.0) / (np.pi*rhow) * lwc / (dm**4)*1e3  #m^-3mm^-1
    parinfo[2]['value'] = nw
    parinfo[2]['fixed'] = 1
  
    if not no_fit:
        print('=================>NORMALIZED GAMMA MU-OPT FITTING')
        #-->global variables to share with fitting routine
        #common tfunc,dbin,nd_obs,dm_obs,nw_obs
        nd_obs = y[l]
        dbin   = x[l]
        dm_obs = dm
        nw_obs = nw
 
        max_iter = 1e3 #max number of iterations to perform
        norm_gamma_fits_mu = optimize.minimize(fit_gamma_norm_mu, mu_mom, method='Nelder-Mead', tol=1e-5, options = {'maxiter':max_iter})
        #norm_gamma_fits_mu=amoeba(1e-5,scale=1e-1,p0=[mu_mom],function_name='fit_gamma_norm_mu',function_value=fval,ncalls=ncalls,nmax=max_iter)
        if(ncalls > 0 and ncalls < max_iter and norm_gamma_fits_mu[0] != parinfo[0]['value']):
            mu = norm_gamma_fits_mu[0]
            f_mu = (6.0 / (4.0**4.0)) * (((mu+4.0)**(mu+4.0)) / math.gamma(4.0+mu))
            yfit_gamma = nw * f_mu * ((x[l]/dm)**mu) * np.exp(-1*(4+mu) * (x[l]/dm))
  
            #-->determine goodness of fit via Chi-square statistic
            dof_gamma = len(l)-1-1
            if dof_gamma <= 0: dof_gamma=1
            chi_gamma = np.sum( (y[l] - yfit_gamma)**2 / y_errors**2)
            #chi_prob_gamma=chisqr_pdf(chi_gamma,dof_gamma) #not used so will ignore in python conversion
            chi_red_gamma = chi_gamma/dof_gamma #reduced Chi-square
   
            n0 = nw*f_mu/dm
   
            lambda_=(4+mu)/dm
            nt = nw * 6.0 * ((4+mu)**3) / ((4**4)*(mu+3)*(mu+2)*(mu+1))*dm  #eqt'n 25 in Williams 2016, JTECH
   
            return {'lwc':lwc, 'nw':nw, 'mu':mu, 'dm':dm, 'lambda':lambda_, 'f_mu':f_mu, 'n0':n0, 'nt':nt, 'chi_red':chi_red_gamma, 'dof':dof_gamma}
        else:
            return -1 

    return {'lwc':lwc, 'nw':nw, 'dm':dm, 'moments':moments_obs} #return this if no fitting is desired (i.e., cannot get mu)
