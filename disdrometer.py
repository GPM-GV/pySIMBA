import os, sys, glob
import numpy as np
from datetime import datetime
import pdb #pdb.set_trace()
import tools as sim

def get_2dvd_for_column(twodvd_files, interval_datetime, box_params, column):

    # Will loop thru possible 2DVD units (based on serial numbers),
    #  if the unit is located within the column box, then will search
    #  the two_dvd_dir for all avail data files for that 2DVD, populate
    #  column grid data arrays - similar to process use for APUs, but a
    #  little different bc may or may not have same # of data files for
    #  each 2DVD unit
    #possible_2DVDs = ['SN25', 'SN35', 'SN36', 'SN37', 'SN38', 'SN70']
    # REMOVE HARD SET HERE - get possible unit IDs from the input dir:
    #all_files = sorted(glob.glob(self.twodvd_dir+'/*sn*'))
    #theList   = strsplit(all_files, '_', /extract)
    #; should be able to ".Filter" theList, but easier for now to just use loop:
    possible_2DVDs = ['']
    for this_file in twodvd_files:
        file_name = os.path.basename(this_file)
        sn_unit = file_name[file_name.find('sn'):file_name.find('sn')+4]
        #file_name = theList[i]
        #sn_spot   = where(strmatch(file_name, 'sn*') eq 1)
        possible_2DVDs.append(sn_unit.upper())
    if len(possible_2DVDs) > 1: possible_2DVDs=np.unique(possible_2DVDs[1:])
    print(f'    possible units: {possible_2DVDs}')
    twodvd_dir = os.path.dirname(twodvd_files[0])

    #print(interval_datetime[0].strftime(f'%Y%m%d%H%M'))
    #print(interval_datetime[self.halftime_interval].strftime(f'%Y%m%d%H%M'))
    #print(interval_datetime[-1].strftime(f'%Y%m%d%H%M'))
    #print('TESTING ^^^ -- comment out')

    # Define holder arrays prior to loop thru possible 2DVDs
    unit_lat_d_array = [-9999]
    unit_lat_m_array = [-9999]
    unit_lat_s_array = [-9999]
    unit_lon_d_array = [-9999]
    unit_lon_m_array = [-9999]
    unit_lon_s_array = [-9999]
    unit_name_array  = np.array(['DUMMY'])
    unit_type_array  = ['DUMMY']
    unit_file_type_array = ['DUMMY']
    unit_timestamp_array = ['DUMMY']
    unit_offset_array    = [-9999]
    unit_int_width_array = [-9999]

    # Define column grid arrays to be populated with available 2DVD data:  need 1 array for each data field, grouped here by file type:
    #  initialize all to NaNs so only populate arrays w/ where there is data available
    #  THESE ARE EACH 5-D ARRAY [drop bins X  times in interval  X  column x dir  X  column y dir  X  column z dir  ]:
    # From _dropcounts. file

    #dims = (50, n_interval_times, len(self.column_box_params['column_grid_lons']),  len(self.column_box_params['column_grid_lats']),  len(self.column_box_params['column_grid_alts'])  )
    #drop_d = (50,)
    #drop_dims = self.ground_plat_dims + drop_d 
    #self.twodvd_dropcount = np.full(drop_dims, np.nan)

    # From _raindsd. OR _raindsd_50%. file:
    #self.twodvd_dsd_con = np.full(drop_dims, np.nan)

    # From _raindsd_100%. file:
    #self.twodvd_dsd_100_con = np.full(drop_dims, np.nan)

    # From _raindsd_ter. OR _raindsd_50%_ter. file:
    #self.twodvd_dsd_ter_con = np.full(drop_dims, np.nan)

    # From _raindsd_100%_ter. file:
    #self.twodvd_dsd_100_ter_con = np.full(drop_dims, np.nan)

    # From _rainparameter. OR _rainparameter_50%. file:
    #  THESE ARE EACH 4-D ARRAY [times in interval  X  column x dir  X  column y dir  X  column z dir  ]:
    #dims = (n_interval_times, len(self.column_box_params['column_grid_lons']),  len(self.column_box_params['column_grid_lats']),  len(self.column_box_params['column_grid_alts'])  )
    twodvd_param_tot_ndrops = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_tot_con    = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_LWC        = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_rainrate   = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_refRayleigh= np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_Dm         = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_Dmax       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_Dmin       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_stdev_Dm   = np.full(box_params.ground_plat_dims, np.nan)

    # From _rainparameter_100%. file:
    twodvd_param_100_tot_ndrops = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_tot_con    = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_LWC        = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_rainrate   = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_refRayleigh= np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_Dm         = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_Dmax       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_Dmin       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_stdev_Dm   = np.full(box_params.ground_plat_dims, np.nan)

    # From _rainparameter_ter. OR _rainparameter_50%_ter. file:
    twodvd_param_ter_tot_ndrops = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_tot_con    = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_LWC        = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_rainrate   = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_refRayleigh= np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_Dm         = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_Dmax       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_Dmin       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_ter_stdev_Dm   = np.full(box_params.ground_plat_dims, np.nan)

    # From rainparameter_100%_ter. file:
    twodvd_param_100_ter_tot_ndrops = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_tot_con    = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_LWC        = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_rainrate   = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_refRayleigh= np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_Dm         = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_Dmax       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_Dmin       = np.full(box_params.ground_plat_dims, np.nan)
    twodvd_param_100_ter_stdev_Dm   = np.full(box_params.ground_plat_dims, np.nan)

    # From _multifreq_attenuation_ter. file:
    #twodvd_atten_rainrate = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_S	      = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_C        = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_X        = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_Ku       = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_K        = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_Ka       = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_atten_W        = np.full(box_params.ground_plat_dims, np.nan)

    # From multifreq_reflectivity_ter. file:
    #twodvd_ref_rainrate = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_Rayleigh = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_S	    = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_C        = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_X        = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_Ku       = np.full(sbox_params.ground_plat_dims, np.nan)
    #twodvd_ref_K        = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_Ka       = np.full(box_params.ground_plat_dims, np.nan)
    #twodvd_ref_W        = np.full(box_params.ground_plat_dims, np.nan)
    
    #units for all fields
    twodvd_units_ndrops = 'number of drops'
    twodvd_units_concen = 'drops m^-3 of air'
    twodvd_units_LWC = 'g m^-3'
    twodvd_units_rainrate = 'mm h^-1'
    twodvd_units_refRayleigh = 'dBZ'
    twodvd_units_Dm = 'mm'
    twodvd_units_Dmax = 'mm'
    twodvd_units_Dmin = 'mm'
    twodvd_units_stdevDm = 'mm'
    #twodvd_units_atten = 'dB km^-1'
    #twodvd_units_ref = 'dBZ'

    # Loop thru all possible 2DVD units:
    for u, unit_name in enumerate(possible_2DVDs):
        #unit_name       = possible_2DVDs[u]
        twodvd_plat_info = sim.set_plat_values(unit_name, 'dummy.txt')
        #unit_loc        = get_platform_loc(unit_name)
        #unit_loc_decdeg = dms2dd(unit_loc[lat_d, unit_loc.lat_m, unit_loc.lat_s, $
        #                   unit_loc.lon_d, unit_loc.lon_m, unit_loc.lon_s)
        unit_loc_decdeg = sim.dms2dd(twodvd_plat_info['lat_d'], twodvd_plat_info['lat_m'],
                                 twodvd_plat_info['lat_s'], twodvd_plat_info['lon_d'],
                                 twodvd_plat_info['lon_m'], twodvd_plat_info['lon_s'])
        unit_lat = unit_loc_decdeg[0]
        unit_lon = unit_loc_decdeg[1]
        lat_ok = 'N'
        lon_ok = 'N'
        if unit_lat <= box_params.box_max_lat and unit_lat >= box_params.box_min_lat: lat_ok = 'Y'
        if unit_lon <= box_params.box_max_lon and unit_lon >= box_params.box_min_lon: lon_ok = 'Y'
        #print, 'lat limits: ', box_min_lat, box_max_lat
        #print, 'unit lat: ', unit_lat
        #print, 'lon limits: ', box_min_lon, box_max_lon
        #print, 'unit lon: ', unit_lon  
        if lat_ok == 'Y' and lon_ok == 'Y':
            # 2DVD unit[u] is located within the column box, procede:

            # Find which of the 2DVD file types exist for this unit:
            # NOTE:  IF FILE NAMING CONVENTION CHANGES, MAY HAVE TO 
            #        ALTER THESE LINES TO LOOK FOR UPDATED STRING(S)
            unit_id = unit_name[2:]
            #pdb.set_trace() #need to add self.twodvd_files below appropriately 
            dropcounts       = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_dropcount*')
            #rain_dsd         = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd.*')
            #if len(rain_dsd) < 1: rain_dsd = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_50%.*')
            #rain_dsd_ter     = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_ter.*')  
            #if len(rain_dsd_ter) < 1: rain_dsd_ter = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_50%_ter.*')
            #rain_dsd_100     = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_100%.*')
            #rain_dsd_100_ter = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_100%_ter.*')  
            rain_param       = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter.*')
            if len(rain_param) < 1: rain_param = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_50%.*')
            rain_param_ter   = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_ter*')
            if len(rain_param_ter) < 1: rain_param_ter = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_50%_ter.*')
            rain_param_100   = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_100%.*')
            rain_param_100_ter=glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_100%_ter.*')
            #atten_ter         = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_multifreq_attenuation_ter*')
            #ref_ter           = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_multifreq_reflectivity_ter*')
            #pdb.set_trace()
            
            # check if file exist or not and print out:
            #self.twodvd_dsd_file = 0
            #self.twodvd_dsd_ter_file = 0
            #self.twodvd_dsd_100_file = 0
            #self.twodvd_dsd_100_ter_file =0
            twodvd_param_file = 0
            twodvd_param_ter_file = 0
            twodvd_param_100_file = 0
            twodvd_param_100_ter_file=0
            twodvd_dropcount_file = 0
            #self.twodvd_atten_file = 0
            #self.twodvd_ref_file = 0
            print('-----------------------------')
            print(f'   avail files for 2DVD unit:  {unit_name}')
            if len(dropcounts)       == 1: 
                twodvd_dropcount_file  = 1
                print(f'     dropcount: {os.path.basename(dropcounts[0])}')
            else:
                print('     dropcount: ')
            #if len(rain_dsd)         == 1: 
            #    self.twodvd_dsd_file = 1
            #    print(f'       50% dsd: {os.path.basename(rain_dsd[0])}')
            #else:
            #    print('       50% dsd: ')
            #if len(rain_dsd_ter)     == 1:
            #    self.twodvd_dsd_ter_file = 1
            #    print(f'   50% dsd ter: {os.path.basename(rain_dsd_ter[0])}')
            #else:
            #    print('   50% dsd ter: ')
            #if len(rain_dsd_100)     == 1:
            #    self.twodvd_dsd_100_file = 1
            #    print(f'      100% dsd: {os.path.basename(rain_dsd_100[0])}')
            #else:
            #    print('      100% dsd: ')
            #if len(rain_dsd_100_ter) == 1:
            #    self.twodvd_dsd_100_ter_file = 1
            #    print(f'  100% dsd ter: {os.path.basename(rain_dsd_100_ter[0])}')
            #else:
            #    print('  100% dsd ter: ')
            if len(rain_param)       == 1:
                twodvd_param_file = 1
                print(f'      50% parm: {os.path.basename(rain_param[0])}')
            else:
                print('      50% parm: ')
            if len(rain_param_ter)   == 1:
                twodvd_param_ter_file = 1
                print(f'  50% parm ter: {os.path.basename(rain_param_ter[0])}')
            else:
                print('  50% parm ter: ')
            if len(rain_param_100)   == 1:
                twodvd_param_100_file = 1
                print(f'    100% param: {os.path.basename(rain_param_100[0])}')
            else:
                print('    100% param: ')
            if len(rain_param_100_ter) == 1:
                twodvd_param_100_ter_file=1
                print(f'100% param ter: {os.path.basename(rain_param_100_ter[0])}')
            else:
                print('100% param ter: ')
            #if len(atten_ter)        == 1:
            #    self.twodvd_atten_file = 1
            #    print(f'   attenuation: {os.path.basename(atten_ter[0])}')
            #else:
            #    print('   attenuation: ')
            #if len(ref_ter)          == 1:
            #    self.twodvd_ref_file = 1
            #    print(f'  reflectivity: {os.path.basename(ref_ter[0])}')
            #else:
            #    print('  reflectivity: ')
            print(f'-------------')
            #pdb.set_trace()
            # If no supported files available for this 2DVD unit, then skip it:
            #if  self.twodvd_dsd_file     == 1 or self.twodvd_dsd_ter_file   == 1 or self.twodvd_dsd_100_file   == 1 or self.twodvd_dsd_100_ter_file == 1 or \
            if  twodvd_param_file   == 1 or twodvd_param_ter_file == 1 or twodvd_param_100_file == 1 or twodvd_param_100_ter_file == 1 or \
                twodvd_dropcount_file == 1 or twodvd_atten_file == 1   or twodvd_ref_file == 1:

                # CALL TO set_plat_values.pro:  FOR LOC ATTRIBUTES.
                unit_in_column_test_cnt = 0   # counter to see if have added this unit to big plat info arrays

                # find column grid point closest to the current 2DVD unit
                #  - arrays of data from the 2DVD files have [field, time] dims
                #   - Field #: which column in the data file to use  (see README_2DVD_data_formats.txt)
                #   - time #:  use the subscript for the timestamp that is needed
                #   - will use NumPy genfromtxt function to read ascii data (this function works well with missing data)
                closest_col_lat_sub = sim.closest(box_params.lat_values, unit_lat)
                closest_col_lon_sub = sim.closest(box_params.lon_values, unit_lon)
                if len(closest_col_lat_sub) != 1 or len(closest_col_lon_sub) != 1:
                    print('-----PROBLEM SETTING WHERE IN GRID TO PLACE 2DVD VALUES!!--')
                    #pdb.set_trace() ##this should not happen -- check and find a way to skip this
                    print(f'---Returning without 2DVD data in the column---')
                    return column
                else:
                    closest_col_lat_sub = closest_col_lat_sub[0]
                    closest_col_lon_sub = closest_col_lon_sub[0]
                
                # NOW - IF KNOW THE UNIT LOC IS IN THE COL BOX, OPEN EA TYPE OF 2DVD FILE,
                #   SEE IF TIMES IN INTERVAL AVAILABLE, POPULATE ARRAYS W DATA FROM THIS UNIT
                if twodvd_dropcount_file == 1:
                    dropcounts_data    = np.genfromtxt(dropcounts) #will return float 2d array
                    dc_year = dropcounts_data[:,0]
                    dc_Jday = dropcounts_data[:,1]
                    dc_hour = dropcounts_data[:,2]
                    dc_min  = dropcounts_data[:,3]
                    dc_datetime = sim.jul2datetime(dc_year, dc_Jday, dc_hour, dc_min)

                    for t, this_date in enumerate(interval_datetime):
                        time_sub = np.where(dc_datetime == this_date)[0]
                    
                        if len(time_sub) < 1:
                            continue  # this time not avail in the 2dvd file
                        elif len(time_sub) > 1:
                            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                            print(f'  --- in 2DVD file: {os.path.basename(dropcounts)}')
                        else:#len(time_sub) == 1
                            # Location OK & current time in interval OK - so proceed:
                            #  Increment counter, set attributes if haven't yet
                            unit_in_column_test_cnt += 1
                            if unit_in_column_test_cnt == 1:
                                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                                unit_name_array = np.append(unit_name_array, twodvd_plat_info['plat_name'])
                                unit_type_array.append(twodvd_plat_info['plat_type'])
                                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                                unit_int_width_array.append(box_params.halftime_interval)
                                # just will set main_plat time w/ sec as '00' since using an interval of data
                                #unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hour+self.main_min+'00')
                                unit_timestamp_array.append(column.time)
                                unit_offset_array.append(0)
                    
                            # populate big column grid w/ 2DVD data from current unit's *dropcount* file for THIS time:
                            # bin 0 starts at column 5 in *dropcount* file or python index 4
                            for binn in range(twodvd_dropcount.shape[0]):
                                twodvd_dropcount[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = dropcounts_data[time_sub, binn+4]
                #pdb.set_trace()
                #if self.twodvd_dsd_file == 1:
                #    rain_dsd_data   = np.genfromtxt(rain_dsd[0])
                #    #rain_dsd_data   = np.loadtxt(rain_dsd)
                #    dsd_year        = rain_dsd_data[:,0]
                #    dsd_Jday        = rain_dsd_data[:,1]
                #    dsd_hour        = rain_dsd_data[:,2]
                #    dsd_min         = rain_dsd_data[:,3]
                #    dsd_datetime    = sim.jul2datetime(dsd_year, dsd_Jday, dsd_hour, dsd_min)

                #    for t, this_date in enumerate(self.interval_datetime):
                #        time_sub  = np.where(dsd_datetime == this_date)[0]
                #        if len(time_sub) < 1:
                #            continue # this time not avail in the 2dvd file, skip this interval time
                #        elif len(time_sub) > 1:
                #            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                #            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                #            print(f'  --- in 2DVD file: {os.path.basename(rain_dsd[0])}')
                #        else:# len(time_sub) == 1
                #            # Location OK & current time in interval OK - so proceed:
                #            #  Increment counter, set attributes if haven't yet
                #            unit_in_column_test_cnt += 1
                #            if unit_in_column_test_cnt == 1:
                #                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                #                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                #                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                #                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                #                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                #                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                #                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                #                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                #                unit_type_array.append(twodvd_plat_info['plat_type'])
                #                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                #                unit_int_width_array.append(self.halftime_interval)
                #                # just will set main_plat time w/ sec as '00' since using an interval of data
                #                unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                #                unit_offset_array.append(0)
                #          
                #            # populate big column grid w/ 2DVD data from current unit's *raindsd.*/*raindsd_50* file for THIS time:
                #            # bin 0 starts at column 5 in *raindsd* file or python index 4
                #            #pdb.set_trace()
                #            for binn in range(self.twodvd_dsd_con.shape[0]):
                #                self.twodvd_dsd_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_data[time_sub, binn+4]

                #if self.twodvd_dsd_100_file == 1:
                #    rain_dsd_100_data   = np.genfromtxt(rain_dsd_100[0])
                #    dsd_100_year        = rain_dsd_100_data[:,0]
                #    dsd_100_Jday        = rain_dsd_100_data[:,1]
                #    dsd_100_hour        = rain_dsd_100_data[:,2]
                #    dsd_100_min         = rain_dsd_100_data[:,3]
                #    dsd_100_datetime    = sim.jul2datetime(dsd_100_year, dsd_100_Jday, dsd_100_hour, dsd_100_min)

                #    for t, this_date in enumerate(self.interval_datetime):
                #        time_sub  = np.where(dsd_100_datetime == this_date)[0]
                #        if len(time_sub) < 1:
                #            continue # this time not avail in the 2dvd file, skip this interval time
                #        elif len(time_sub) > 1:
                #            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                #            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                #            print(f'  --- in 2DVD file: {os.path.basename(rain_dsd_100[0])}')
                #        else:# len(time_sub) == 1
                #            # Location OK & current time in interval OK - so proceed:
                #            #  Increment counter, set attributes if haven't yet
                #            unit_in_column_test_cnt += 1
                #            if unit_in_column_test_cnt == 1:
                #                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                #                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                #                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                #                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                #                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                #                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                #                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                #                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                #                unit_type_array.append(twodvd_plat_info['plat_type'])
                #                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                #                unit_int_width_array.append(self.halftime_interval)
                #                # just will set main_plat time w/ sec as '00' since using an interval of data
                #                unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                #                unit_offset_array.append(0)
                #        
                #            # populate big column grid w/ 2DVD data from current unit's *raindsd_100%.* file for THIS time:
                #            # bin 0 starts at column 5 in *raindsd* file or python index 4
                #            for binn in range(self.twodvd_dsd_100_con.shape[0]):
                #                self.twodvd_dsd_100_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_100_data[time_sub, binn+4]

                #if self.twodvd_dsd_ter_file == 1:
                #    rain_dsd_ter_data   = np.genfromtxt(rain_dsd_ter[0])
                #    dsd_ter_year        = rain_dsd_ter_data[:,0]
                #    dsd_ter_Jday        = rain_dsd_ter_data[:,1]
                #    dsd_ter_hour        = rain_dsd_ter_data[:,2]
                #    dsd_ter_min         = rain_dsd_ter_data[:,3]
                #    dsd_ter_datetime    = sim.jul2datetime(dsd_ter_year, dsd_ter_Jday, dsd_ter_hour, dsd_ter_min)

                #    for t, this_date in enumerate(self.interval_datetime):
                #        time_sub  = np.where(dsd_ter_datetime == this_date)[0]
                #        if len(time_sub) < 1:
                #            continue # this time not avail in the 2dvd file, skip this interval time
                #        elif len(time_sub) > 1:
                #            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                #            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                #            print(f'  --- in 2DVD file: {os.path.basename(rain_dsd_ter[0])}')
                #        else:# len(time_sub) == 1
                #        # Location OK & current time in interval OK - so proceed:
                #        #  Increment counter, set attributes if haven't yet
                #            unit_in_column_test_cnt += 1
                #            if unit_in_column_test_cnt == 1:
                #                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                #                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                #                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                #                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                #                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                #                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                #                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                #                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                #                unit_type_array.append(twodvd_plat_info['plat_type'])
                #                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                #                unit_int_width_array.append(self.halftime_interval)
                #                # just will set main_plat time w/ sec as '00' since using an interval of data
                #                unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                #                unit_offset_array.append(0)

                #            # populate big column grid w/ 2DVD data from current unit's *raindsd_ter.* file for THIS time:
                #            # bin 0 starts at column 5 in *raindsd* file or python index 4
                #            dsd_ter_main_time_sub = time_sub
                #            for binn in range(self.twodvd_dsd_ter_con.shape[0]):
                #                self.twodvd_dsd_ter_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_ter_data[dsd_ter_main_time_sub, binn+4]

                #if self.twodvd_dsd_100_ter_file == 1:
                #    rain_dsd_100_ter_data      = np.genfromtxt(rain_dsd_100_ter[0])
                #    dsd_100_ter_year           = rain_dsd_100_ter_data[:,0]
                #    dsd_100_ter_Jday           = rain_dsd_100_ter_data[:,1]
                #    dsd_100_ter_hour           = rain_dsd_100_ter_data[:,2]
                #    dsd_100_ter_min            = rain_dsd_100_ter_data[:,3]
                #    dsd_100_ter_datetime       = sim.jul2datetime(dsd_100_ter_year, dsd_100_ter_Jday, dsd_100_ter_hour, dsd_100_ter_min)

                #    for t, this_date in enumerate(self.interval_datetime):
                #        time_sub  = np.where(dsd_100_ter_datetime == this_date)[0]
                #        if len(time_sub) < 1:
                #            continue # this time not avail in the 2dvd file, skip this interval time
                #        elif len(time_sub) > 1:
                #            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                #            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                #            print(f'  --- in 2DVD file: {os.path.basename(rain_dsd_100_ter[0])}')
                #        else:# len(time_sub) == 1
                #            # Location OK & current time in interval OK - so proceed:
                #            #  Increment counter, set attributes if haven't yet
                #            unit_in_column_test_cnt += 1
                #            if unit_in_column_test_cnt == 1:
                #                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                #                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                #                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                #                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                #                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                #                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                #                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                #                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                #                unit_type_array.append(twodvd_plat_info['plat_type'])
                #                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                #                unit_int_width_array.append(self.halftime_interval)
                #                # just will set main_plat time w/ sec as '00' since using an interval of data
                #                unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                #                unit_offset_array.append(0)

                #            # populate big column grid w/ 2DVD data from current unit's *raindsd_100%_ter.* file for THIS time:
                #            # bin 0 starts at column 5 in *raindsd* file or python index 4
                #            dsd_100_ter_main_time_sub = time_sub
                #            for binn in range(self.twodvd_dsd_100_ter_con.shape[0]):
                #                self.twodvd_dsd_100_ter_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_100_ter_data[dsd_100_ter_main_time_sub, binn+4]

                if twodvd_param_file == 1:
                    rain_param_data = np.genfromtxt(rain_param[0])
                    param_year      = rain_param_data[:,0]
                    param_Jday      = rain_param_data[:,1]
                    param_hour      = rain_param_data[:,2]
                    param_min       = rain_param_data[:,3]
                    param_datetime  = sim.jul2datetime(param_year, param_Jday, param_hour, param_min)

                    for t, this_date in enumerate(interval_datetime):
                        time_sub  = np.where(param_datetime == this_date)[0]
                        if len(time_sub) < 1:
                            continue # this time not avail in the 2dvd file, skip this interval time
                        elif len(time_sub) > 1:
                            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                            print(f'  --- in 2DVD file: {os.path.basename(rain_param[0])}')
                        else:# len(time_sub) == 1
                            # Location OK & current time in interval OK - so proceed:
                            #  Increment counter, set attributes if haven't yet
                            unit_in_column_test_cnt += 1
                            if unit_in_column_test_cnt == 1:
                                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                                unit_type_array.append(twodvd_plat_info['plat_type'])
                                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                                unit_int_width_array.append(box_params.halftime_interval)
                                # just will set main_plat time w/ sec as '00' since using an interval of data
                                #unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                unit_timestamp_array.append(column.time)
                                unit_offset_array.append(0)

                            # populate column grid w/ 2DVD data from current unit's *parameter.*/*parameter_50* file:
                            param_main_time_sub = time_sub
                            twodvd_param_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_data[param_main_time_sub,4]
                            twodvd_param_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_data[param_main_time_sub,5]
                            twodvd_param_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_data[param_main_time_sub,6]
                            twodvd_param_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[param_main_time_sub,7]
                            twodvd_param_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_data[param_main_time_sub,8]
                            twodvd_param_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_data[param_main_time_sub,9]
                            twodvd_param_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[param_main_time_sub,10]
                            twodvd_param_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[param_main_time_sub,11]
                            twodvd_param_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[param_main_time_sub,12]

                if twodvd_param_100_file == 1:
                    rain_param_100_data = np.genfromtxt(rain_param_100[0])
                    param_100_year      = rain_param_100_data[:,0]
                    param_100_Jday      = rain_param_100_data[:,1]
                    param_100_hour      = rain_param_100_data[:,2]
                    param_100_min       = rain_param_100_data[:,3]
                    param_100_datetime = sim.jul2datetime(param_100_year, param_100_Jday, param_100_hour, param_100_min)

                    for t, this_date in enumerate(interval_datetime):
                        time_sub  = np.where(param_100_datetime == this_date)[0]
                        if len(time_sub) < 1:
                            continue # this time not avail in the 2dvd file, skip this interval time
                        elif len(time_sub) > 1:
                            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                            print(f'  --- in 2DVD file: {os.path.basename(rain_param_100[0])}')
                        else:# len(time_sub) == 1
                            # Location OK & current time in interval OK - so proceed:
                            #  Increment counter, set attributes if haven't yet
                            unit_in_column_test_cnt += 1
                            if unit_in_column_test_cnt == 1:
                                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                                unit_type_array.append(twodvd_plat_info['plat_type'])
                                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                                unit_int_width_array.append(halftime_interval)
                                # just will set main_plat time w/ sec as '00' since using an interval of data
                                #unit_timestamp_array.append(main_year+main_month+main_day+'_'+main_hour+main_min+'00')
                                unit_timestamp_array.append(column.time)
                                unit_offset_array.append(0)

                            # populate column grid w/ 2DVD data from current unit's *parameter_100%.* file:
                            twodvd_param_100_main_time_sub = time_sub
                            twodvd_param_100_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_100_data[param_100_main_time_sub,4]
                            twodvd_param_100_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_100_data[param_100_main_time_sub,5]
                            twodvd_param_100_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_100_data[param_100_main_time_sub,6]
                            twodvd_param_100_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_100_data[param_100_main_time_sub,7]
                            twodvd_param_100_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_100_data[param_100_main_time_sub,8]
                            twodvd_param_100_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_100_data[param_100_main_time_sub,9]
                            twodvd_param_100_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_100_data[param_100_main_time_sub,10]
                            twodvd_param_100_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_100_data[param_100_main_time_sub,11]
                            twodvd_param_100_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_100_data[param_100_main_time_sub,12]

                if twodvd_param_ter_file == 1:
                    rain_param_ter_data = np.genfromtxt(rain_param_ter[0])
                    param_ter_year      = rain_param_ter_data[:,0]
                    param_ter_Jday      = rain_param_ter_data[:,1]
                    param_ter_hour      = rain_param_ter_data[:,2]
                    param_ter_min       = rain_param_ter_data[:,3]
                    param_ter_datetime = sim.jul2datetime(param_ter_year, param_ter_Jday, param_ter_hour, param_ter_min)

                    for t, this_date in enumerate(interval_datetime):
                        time_sub  = np.where(param_ter_datetime == this_date)[0]
                        if len(time_sub) < 1:
                            continue # this time not avail in the 2dvd file, skip thos interval time
                        elif len(time_sub) > 1:
                            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                            print(f'  --- in 2DVD file: {os.path.basename(rain_param_ter[0])}')
                        else:# len(time_sub) == 1    
                            # Location OK & current time in interval OK - so proceed:
                            #  Increment counter, set attributes if haven't yet
                            unit_in_column_test_cnt += 1
                            if unit_in_column_test_cnt == 1:
                                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                                unit_type_array.append(twodvd_plat_info['plat_type'])
                                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                                unit_int_width_array.append(halftime_interval)
                                # just will set main_plat time w/ sec as '00' since using an interval of data
                                #unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                unit_timestamp_array.append(column.time)
                                unit_offset_array.append(0)

                            # populate column grid w/ 2DVD data from current unit's *parameter_ter.* file:
                            twodvd_param_ter_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_data[time_sub,4]
                            twodvd_param_ter_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_data[time_sub,5]
                            twodvd_param_ter_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_data[time_sub,6]
                            twodvd_param_ter_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,7]
                            twodvd_param_ter_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_data[time_sub,8]
                            twodvd_param_ter_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_data[time_sub,9]
                            twodvd_param_ter_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,10]
                            twodvd_param_ter_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,11]
                            twodvd_param_ter_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,12]

                if twodvd_param_100_ter_file == 1:
                    rain_param_100_ter_data = np.genfromtxt(rain_param_100_ter[0])
                    param_100_ter_year      = rain_param_100_ter_data[:,0]
                    param_100_ter_Jday      = rain_param_100_ter_data[:,1]
                    param_100_ter_hour      = rain_param_100_ter_data[:,2]
                    param_100_ter_min       = rain_param_100_ter_data[:,3]
                    param_100_ter_datetime = sim.jul2datetime(param_ter_year, param_100_ter_Jday, param_100_ter_hour, param_100_ter_min)

                    for t, this_date in enumerate(interval_datetime):
                        time_sub  = np.where(param_100_ter_datetime == this_date)[0]
                        if len(time_sub) < 1:
                            continue # this time not avail in the 2dvd file, skip this interval time
                        elif len(time_sub) > 1:
                            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                            print(f'  --- in 2DVD file: {os.path.basename(rain_param_100_ter[0])}')
                        else:# len(time_sub) == 1
                            # Location OK & current time in interval OK - so proceed:
                            #  Increment counter, set attributes if haven't yet
                            unit_in_column_test_cnt += 1
                            if unit_in_column_test_cnt == 1:
                                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                                unit_type_array.append(twodvd_plat_info['plat_type'])
                                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                                unit_int_width_array.append(halftime_interval)
                                # just will set main_plat time w/ sec as '00' since using an interval of data
                                #unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                unit_timestamp_array.append(column.time)
                                unit_offset_array.append(0)

                            # populate column grid w/ 2DVD data from current unit's *parameter_100%_ter.* file:
                            twodvd_param_100_ter_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_data[time_sub,4]
                            twodvd_param_100_ter_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_data[time_sub,5]
                            twodvd_param_100_ter_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_data[time_sub,6]
                            twodvd_param_100_ter_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,7]
                            twodvd_param_100_ter_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_data[time_sub,8]
                            twodvd_param_100_ter_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_data[time_sub,9]
                            twodvd_param_100_ter_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,10]
                            twodvd_param_100_ter_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,11]
                            twodvd_param_100_ter_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,12]

                #if twodvd_atten_file == 1:
                #    atten_ter_data  = np.genfromtxt(atten_ter[0])
                #    atten_year      = atten_ter_data[:,0]
                #    atten_Jday      = atten_ter_data[:,1]
                #    atten_hour      = atten_ter_data[:,2]
                #    atten_min       = atten_ter_data[:,3]
                #    atten_datetime = sim.jul2datetime(atten_year, atten_Jday, atten_hour, atten_min)

                #    for t, this_date in enumerate(self.interval_datetime):
                #        time_sub  = np.where(atten_datetime == this_date)[0]
                #        if len(time_sub) < 1:
                #            continue #this time not avail in the 2dvd file, skip this interval time
                #        elif len(time_sub) > 1:
                #            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                #            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                #            print(f'  --- in 2DVD file: {os.path.basename(atten_ter[0])}')
                #        else:# len(time_sub) == 1
                #            # Location OK & current time in interval OK - so proceed:
                #            #  Increment counter, set attributes if haven't yet
                #            unit_in_column_test_cnt += 1
                #            if unit_in_column_test_cnt == 1:
                #                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                #                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                #                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                #                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                #                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                #                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                #                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                #                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                #                unit_type_array.append(twodvd_plat_info['plat_type'])
                #                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                #                unit_int_width_array.append(self.halftime_interval)
                #                # just will set main_plat time w/ sec as '00' since using an interval of data
                #                unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                #                unit_offset_array.append(0)

                #            # populate column grid w/ 2DVD data from current unit's *multifreq_attenuation_ter* file:
                #            self.twodvd_atten_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, 0] = atten_ter_data[time_sub,4]
                #            self.twodvd_atten_S[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,5]
                #            self.twodvd_atten_C[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,6]
                #            self.twodvd_atten_X[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,7]
                #            self.twodvd_atten_Ku[0, closest_col_lat_sub, closest_col_lon_sub, t]       = atten_ter_data[time_sub,8]
                #            self.twodvd_atten_K[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,9]
                #            self.twodvd_atten_Ka[0, closest_col_lat_sub, closest_col_lon_sub, t]       = atten_ter_data[time_sub,10]
                #            self.twodvd_atten_W[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,11]

                #if self.twodvd_ref_file  ==  1:
                #    ref_ter_data    = np.genfromtxt(ref_ter[0])
                #    ref_year        = ref_ter_data[:,0]
                #    ref_Jday        = ref_ter_data[:,1]
                #    ref_hour        = ref_ter_data[:,2]
                #    ref_min         = ref_ter_data[:,3]
                #    ref_datetime    = sim.jul2datetime(ref_year, ref_Jday, ref_hour, ref_min)

                #    for t, this_date in enumerate(self.interval_datetime):
                #        time_sub  = np.where(ref_datetime == this_date)[0]
                #        if len(time_sub) < 1:
                #            continue # this time not avail in the 2dvd file, skip this interval time
                #        elif len(time_sub) > 1:
                #            # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                #            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                #            print(f'  --- in 2DVD file: {os.path.basename(ref_ter[0])}')
                #        else:# len(time_sub) == 1
                #            # Location OK & current time in interval OK - so proceed:
                #            #  Increment counter, set attributes if haven't yet
                #            unit_in_column_test_cnt += 1
                #            if unit_in_column_test_cnt == 1:
                #                # if 1st iteration w/ time & loc OK, set attributes for this unit into big arrays
                #                unit_lat_d_array.append(twodvd_plat_info['lat_d'])
                #                unit_lat_m_array.append(twodvd_plat_info['lat_m'])
                #                unit_lat_s_array.append(twodvd_plat_info['lat_s'])
                #                unit_lon_d_array.append(twodvd_plat_info['lon_d'])
                #                unit_lon_m_array.append(twodvd_plat_info['lon_m'])
                #                unit_lon_s_array.append(twodvd_plat_info['lon_s'])
                #                unit_name_array = np.append(unit_name_array,twodvd_plat_info['plat_name'])
                #                unit_type_array.append(twodvd_plat_info['plat_type'])
                #                unit_file_type_array.append(twodvd_plat_info['operation_mode'])
                #                unit_int_width_array.append(self.halftime_interval)
                #                # just will set main_plat time w/ sec as '00' since using an interval of data
                #                unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                #                unit_offset_array.append(0)

                #            # populate column grid w/ 2DVD data from current unit's *multifreq_reflectivity_ter* file:
                #            self.twodvd_ref_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t] = ref_ter_data[time_sub,4]
                #            self.twodvd_ref_Rayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t] = ref_ter_data[time_sub,5]
                #            self.twodvd_ref_S[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,6]
                #            self.twodvd_ref_C[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,7]
                #            self.twodvd_ref_X[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,8]
                #            self.twodvd_ref_Ku[0, closest_col_lat_sub, closest_col_lon_sub, t]       = ref_ter_data[time_sub,9]
                #            self.twodvd_ref_K[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,10]
                #            self.twodvd_ref_Ka[0, closest_col_lat_sub, closest_col_lon_sub, t]       = ref_ter_data[time_sub,11]
                #            self.twodvd_ref_W[0, closest_col_lat_sub, closest_col_lon_sub, t]           = ref_ter_data[time_sub,12]

    # If there are 2DVDs in the column box, must remove dummy element at start of holder info arrays:
    if len(unit_lat_d_array) > 1:
        unit_lat_d_array = unit_lat_d_array[1:]
        unit_lat_m_array = unit_lat_m_array[1:]
        unit_lat_s_array = unit_lat_s_array[1:]
        unit_lon_d_array = unit_lon_d_array[1:]
        unit_lon_m_array = unit_lon_m_array[1:]
        unit_lon_s_array = unit_lon_s_array[1:]
        unit_name_array  = unit_name_array[1:]
        unit_type_array  = unit_type_array[1:]
        unit_file_type_array = unit_file_type_array[1:]
        unit_int_width_array = unit_int_width_array[1:]
        unit_timestamp_array = unit_timestamp_array[1:]
        unit_offset_array    = unit_offset_array[1:]
    else:
        #If there are NO 2DVDs in column box, the data arrays will still be all NaNs
        #and so will just return  a -1 to the top level .pro:
        print('There are no 2DVDs in the column box!')
        print(f'---Returning without 2DVD data in the column---')
        return column
    print(f'------ 2DVDs in the column box: {unit_name_array}')

    # Set up dictionary of plat_infos for the 2DVDs in the column box:
    twodvd_info = {'plat_name':unit_name_array,
                   'lat_d':unit_lat_d_array, 'lat_m':unit_lat_m_array, 'lat_s':unit_lat_s_array,
                   'lon_d':unit_lon_d_array, 'lon_m':unit_lon_m_array, 'lon_s':unit_lon_s_array,
                   'plat_type':unit_type_array,
                   'operation_mode':unit_file_type_array,
                   'timestamp':unit_timestamp_array,
                   'offset_vs_main':unit_offset_array,
                   'time_interval_width':unit_int_width_array}
                           
    #add twodvd_info to object
    column.add_platform_to_object(twodvd_info, plat_name='twodvd')
                           
    #assign each field to column object
    if twodvd_param_file:
        column.add_variable_to_object(twodvd_param_tot_ndrops, 
                                      var_name='twodvd_param_tot_ndrops', 
                                      units=twodvd_units_ndrops)

        column.add_variable_to_object(twodvd_param_tot_con, 
                                      var_name='twodvd_param_tot_concen', 
                                      units=twodvd_units_concen)

        column.add_variable_to_object(twodvd_param_LWC, 
                                      var_name='twodvd_param_LWC', 
                                      units=twodvd_units_LWC)

        column.add_variable_to_object(twodvd_param_rainrate, 
                                      var_name='twodvd_param_rainrate', 
                                      units=twodvd_units_rainrate)

        column.add_variable_to_object(twodvd_param_refRayleigh, 
                                      var_name='twodvd_param_refRayleigh', 
                                      units=twodvd_units_refRayleigh)

        column.add_variable_to_object(twodvd_param_Dm, 
                                      var_name='twodvd_param_Dm', 
                                      units=twodvd_units_Dm)

        column.add_variable_to_object(twodvd_param_Dmax, 
                                      var_name='twodvd_param_Dmax', 
                                      units=twodvd_units_Dmax)

        column.add_variable_to_object(twodvd_param_Dmin, 
                                      var_name='twodvd_param_Dmin', 
                                      units=twodvd_units_Dmin)

        column.add_variable_to_object(twodvd_param_stdev_Dm, 
                                      var_name='twodvd_param_stdevDm', 
                                      units=twodvd_units_stdevDm)
    if twodvd_param_ter_file:
        column.add_variable_to_object(twodvd_param_ter_tot_ndrops, 
                                      var_name='twodvd_param_ter_tot_ndrops', 
                                      units=twodvd_units_ndrops)

        column.add_variable_to_object(twodvd_param_ter_tot_con, 
                                      var_name='twodvd_param_ter_tot_concen', 
                                      units=twodvd_units_concen)

        column.add_variable_to_object(twodvd_param_ter_LWC, 
                                      var_name='twodvd_param_ter_LWC', 
                                      units=twodvd_units_LWC)

        column.add_variable_to_object(twodvd_param_ter_rainrate, 
                                      var_name='twodvd_param_ter_rainrate', 
                                      units=twodvd_units_rainrate)

        column.add_variable_to_object(twodvd_param_ter_refRayleigh, 
                                      var_name='twodvd_param_ter_refRayleigh', 
                                      units=twodvd_units_refRayleigh)

        column.add_variable_to_object(twodvd_param_ter_Dm, 
                                      var_name='twodvd_param_ter_Dm', 
                                      units=twodvd_units_Dm)

        column.add_variable_to_object(twodvd_param_ter_Dmax, 
                                      var_name='twodvd_param_ter_Dmax', 
                                      units=twodvd_units_Dmax)

        column.add_variable_to_object(twodvd_param_ter_Dmin, 
                                      var_name='twodvd_param_ter_Dmin', 
                                      units=twodvd_units_Dmin)

        column.add_variable_to_object(twodvd_param_ter_stdev_Dm, 
                                      var_name='twodvd_param_ter_stdevDm', 
                                      units=twodvd_units_stdevDm)

    return column
    #self.twodvd_data_in_column = True

 # ***************************************************************************************

def get_parsivel_for_column(apu_dir, interval_datetime, box_params, column):
    '''
    ;------------------------------------------------------------
    ; Of the available APUs, identify one(s) in column box, return data
    ; at main_plat time +/- halftime_interval back to build_column.pro
    ;------------------------------------------------------------
    ;   -> APU locations MUST BE CORRECT in platform_location.pro
    ;
    ;   -> apu_dir: dir with avail APU files for the present case
    ;		(files processed by Ali Tokay).
    ;
    ;   -> Currently: set up to work with these kind of files:
    ;		- "apu##_raindsd_min."
    ;		- "apu##_raindsd_min_ter."
    ;		- "apu##_rainparameter_min."
    ;		- "apu##_rainparameter_min_ter."
    ;		- "apu##_reflectivity_ter."
    ;		- "apu##_attenuation_ter."
    ;
    ;   -> assumes have same times in ea type of file for same APU##
    ;
    ;   result = get_apu_for_column(column_box_params, apu_dir, $
    ;				main_plat_into, halftime_interval)
    ;	column_box_params:   column box grid parameters. This is the 
    ;				structure returned by define_box.pro
    ;	apu_dir:	     string dir path to all & only APU files
    ;				for the present case
    ;	main_plat_info:      main_plat_info structure returned by 
    ;				set_main_plat_values.pro
    ;	halftime_interval:   integer # of mins before & after the main
    ;				platform time to include in column file
    ;
    ;	result:	 returns a big structure: 
    ;		apu_info_and_data:  This is made up of multiple structures:
    ;
    ;		.apu_plat_info:   structure of arrays of info
    ;				  for ea. avail APU in the col box
    ;		.col_data_apu_rain:  structure of fields from Ali's
    ;					_rainparameter_min. APU files,
    ;					or -1 if that file not available
    ;		.col_data_apu_rain_ter:  structure of fields from Ali's
    ;					_rainparameter_min_ter. APU files,
    ;					or -1 if that file not available
    ;		.col_data_apu_dsd:  structure of fields from Ali's
    ;					_raindsd_min. APU files,
    ;					or -1 if that file not available
    ;		.col_data_apu_dsd_ter:  structure of fields from Ali's
    ;					_raindsd_min_ter. APU files,
    ;					or -1 if that file not available
    ;		.col_data_apu_reflectivity:  struct of fields from Ali's
    ;					_reflectivity_ APU files, 
    ;					or -1 if that file not available
    ;		.col_data_apu_attenuation:   struct of fileds from Ali's
    ;					_attenuation_ APU files,
    ;					or -1 if that file not available
    ;
    ;	   (also see more/full details in READMEs)
    ;
    ;	Dependencies:
    ;	  dms2dd.pro	utility to convert deg-min-sec to decimal degrees
    ;	  jd2time.pro	IDL Coyote utility to convert from Julian day
    ;	  closest.pro	utility to locate value in an array nearest input value
    ;
    ;------------------------------------------------------------
    '''
    # Get list of each type of file in the APU directory:
    #pdb.set_trace()#
    rain_param_files     = sorted(glob.glob(apu_dir+'/*_rainparameter_min.*'))
    rain_param_ter_files = sorted(glob.glob(apu_dir+'/*_rainparameter_min_ter.*'))
    rain_dsd_files       = sorted(glob.glob(apu_dir+'/*_raindsd_min.*'))
    rain_dsd_ter_files   = sorted(glob.glob(apu_dir+'/*_raindsd_min_ter.*'))
    #apu_ref_files        = sorted(glob.glob(apu_dir+'/*_reflectivity_ter.*'))
    #apu_atten_files      = sorted(glob.glob(apu_dir+'/*_attenuation_ter.*'))
    n_rain               = len(rain_param_files)
    n_rain_ter           = len(rain_param_ter_files)
    n_dsd                = len(rain_dsd_files)
    n_dsd_ter            = len(rain_dsd_ter_files)
    #n_ref                = len(apu_ref_files)
    #n_atten              = len(apu_atten_files)
    #pdb.set_trace()
    if n_rain <= 1 or n_rain_ter <= 1 or n_dsd <= 1 or n_dsd_ter <= 1 or n_ref <= 1 or n_atten <= 1:
        if n_rain <= 1:
            if not rain_param_files: n_rain=0
        if n_rain_ter <= 1:
            if not rain_param_ter_files: n_rain_ter=0
        if n_dsd <= 1:
            if not rain_dsd_files: n_dsd=0
        if n_dsd_ter <= 1:
            if not rain_dsd_ter_files: n_dsd_ter=0
        #if n_ref <= 1:
        #    if not apu_ref_files: n_ref=0
        #if n_atten <= 1:
        #    if not apu_atten_files: n_atten=0

    # for files that are available, test to be sure have same # of files per ea APU
    #list_num_of_files = np.array([n_rain, n_rain_ter, n_dsd, n_dsd_ter, n_ref, n_atten])
    list_num_of_files = np.array([n_rain, n_rain_ter, n_dsd, n_dsd_ter])
    have_no_files     = np.where(list_num_of_files == 0)[0]
    do_have_files     = np.where(list_num_of_files != 0)[0]
    avail_files       = list_num_of_files[do_have_files]
    n_avail_apus_sub  = np.unique(avail_files, return_index=True)[1]

    #print(list_num_of_files)
    #print(have_no_files)
    #print(do_have_files)
    #print(avail_files)
    #print(n_avail_apus_sub)
    if len(n_avail_apus_sub) != 1:
        print('-----HAVE DIFFERENT # OF THE VARIOUS APU FILE TYPES!!---')
        print('----verify apu_dir is correct, has correct/all APU files')
        print(f'---Returning without parsivel data in the column---')
        return column
    else:
        n_avail_apus = avail_files[n_avail_apus_sub[0]]
        #should always have the _rainparameter_min_ter file, so 
        #will start w/ that one in main loop below...

    #print(interval_datetime[0].strftime(f'%Y%m%d%H%M'))
    #print(interval_datetime[self.halftime_interval].strftime(f'%Y%m%d%H%M'))
    #print(interval_datetime[-1].strftime(f'%Y%m%d%H%M'))
    #print('TESTING ^^^ -- comment out')

    # Define holder arrays prior to going thru each APU
    #   will add to these for ea APU in the column box
    #   then remove 1st/dummy element before pass back 
    apu_lat_d_array = [-9999]
    apu_lat_m_array = [-9999]
    apu_lat_s_array = [-9999]
    apu_lon_d_array = [-9999]
    apu_lon_m_array = [-9999]
    apu_lon_s_array = [-9999]
    apu_plat_name_array = np.array(['DUMMY'])
    apu_plat_type_array = ['DUMMY']
    apu_file_type_array = ['DUMMY']
    apu_plat_timestamp_array = ['DUMMY']
    apu_offset_array = [-9999]
    apu_int_width_array = [-9999]

    # Define varaibles: set all initial vals to NaNs, will fill in where have
    #  APU values from each file, if an APU file type is not avail, column 
    #  .nc file will not contain variables for those fields. 
    #  THESE ARE EACH 4-D ARRAYS:  [  times in interval  X  column x dir  X  column y dir  X  column z dir  ]
    #  Variables for drop bins are EACH 5-D ARRAYS:  [drop bins X  times in interval  X  column x dir  X  column y dir  X  column z dir  ]
    #dims = (n_interval_times, len(self.column_box_params['column_grid_lons']),  len(self.column_box_params['column_grid_lats']),  len(self.column_box_params['column_grid_alts'])  )
    drop_d = (32,)
    dims_drop_bins = box_params.ground_plat_dims + drop_d
    #pdb.set_trace()
    if n_rain >= 1: # variables from file:  _rainparameter_min
        apu_n_drops            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_concentration      = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_liqwater           = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_rain_rainrate      = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_rain_refinRayleigh = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_massweight_diam    = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_maximum_diam       = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    if n_rain_ter >= 1: # variables from file:  _rainparameter_min_ter
        apu_n_drops_ter            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_concentration_ter      = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_liqwater_ter           = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_rain_rainrate_ter      = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_rain_refinRayleigh_ter = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_massweight_diam_ter    = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
        apu_maximum_diam_ter       = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    if n_dsd >= 1 or n_dsd_ter >= 1:
        apu_dsd_bins_Dmiddle = np.zeros(32)  # hold the bin middle diameters [mm]
        apu_dsd_bins_widths  = np.zeros(32)  # hold the bin widths [mm]
    if n_dsd >= 1: # variables from file:  _raindsd_min
        apu_drop_concen = np.full(dims_drop_bins, np.nan, dtype=object)
    if n_dsd_ter >= 1:  # variables from file:  _raindsd_min_ter
        apu_drop_concen_ter = np.full(dims_drop_bins, np.nan, dtype=object)
    #if n_ref >= 1: # variables from file:  _reflectivity_ter
    #    apu_ref_rainrate       = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_refinRayleigh  = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atS            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atC            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atX            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atKu           = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atK            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atKa           = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_ref_atW            = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #if n_atten >= 1: # variables from file:  _attenuation_ter
    #    apu_atten_rainrate   = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atS        = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atC        = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atX        = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atKu       = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atK        = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atKa       = np.full(box_params.ground_plat_dims, np.nan, dtype=object)
    #    apu_atten_atW        = np.full(box_params.ground_plat_dims, np.nan, dtype=object)

    # Loop thru available APUs for ea file type:
    #  Determine if APU is located within column box
    #     if so:  
    #	- locate APU values for main_plat_time 
    #		-> place value at grid point closest to APU
    #		-> assign to APU main_time data variables
    #  Starts w/ _rainparameter_min_ter file, then goes to 
    #	other APU file types if available.  This means the
    # 	MODULE ASSUMES the _rainparameter_min_ter file is
    #	always available if an APU was available.
    #pdb.set_trace()
    if n_avail_apus >= 1:
        for i in range(n_avail_apus):
        
            # get plat_info for current APU
            rain_ter_file = rain_param_ter_files[i]
            plat_name = os.path.basename(rain_ter_file)[0:5]  #get APU name as 'apu' + 2char string
            apu_plat_info = sim.set_plat_values(plat_name, rain_ter_file)
        
            # determine if current APU is within bounds of the column box grid
            # get apu loc in dec deg
            apu_loc_decdeg = sim.dms2dd(apu_plat_info['lat_d'], apu_plat_info['lat_m'],
                                    apu_plat_info['lat_s'], apu_plat_info['lon_d'],
                                    apu_plat_info['lon_m'], apu_plat_info['lon_s'])
            apu_lat = apu_loc_decdeg[0]
            apu_lon = apu_loc_decdeg[1]
            #print, 'apu: ',apu_lat, apu_lon
            #print, 'box limits: ',box_min_lat, box_max_lat, box_min_lon, box_max_lon
            lat_ok = 'N'
            lon_ok = 'N'
            if apu_lat <= box_params.box_max_lat and apu_lat >= box_params.box_min_lat: lat_ok = 'Y'
            if apu_lon <= box_params.box_max_lon and apu_lon >= box_params.box_min_lon: lon_ok = 'Y'
            #print, '        apu:  ',plat_name
            #print, ' latlon OK?: ',lat_ok,lon_ok
            if lat_ok == 'Y' and lon_ok == 'Y':
                print(f'    APU location in box: {apu_plat_info["plat_name"]}')

                #APU Unit is in column box, so read in available files:
                if n_rain_ter > 0:
                    apu_rain_ter_file = rain_param_ter_files[i]
                    rain_ter_data_full = np.loadtxt(apu_rain_ter_file)
                else:
                    apu_rain_ter_file = None
                if n_rain > 0:
                    apu_rain_file = rain_param_files[i]
                    rain_data_full = np.loadtxt(apu_rain_file)
                else:
                    apu_rain_file = None
                if n_dsd > 0:
                    apu_dsd_file = rain_dsd_files[i]
                    dsd_data_full  = np.loadtxt(apu_dsd_file)
                else:
                    apu_dsd_file = None
                if n_dsd_ter > 0:
                    apu_dsd_ter_file = rain_dsd_ter_files[i]
                    dsd_ter_data_full = np.loadtxt(apu_dsd_ter_file)
                else:
                    apu_dsd_ter_file = None
                #if n_ref > 0:
                #    apu_ref_file = apu_ref_files[i]
                #    ref_data_full  = np.loadtxt(apu_ref_file)
                #else:
                #    apu_ref_file = None
                #if n_atten > 0:
                #    apu_atten_file = apu_atten_files[i]
                #    atten_data_full= np.loadtxt(apu_atten_file)
                #else:
                #    apu_atten_file = None

                # Get time info from _rainparameter_min_ter file data 
                # and convert to datetime objective
                apu_year = rain_ter_data_full[:,0]
                apu_Jday = rain_ter_data_full[:,1]
                apu_hr   = rain_ter_data_full[:,2]
                apu_min  = rain_ter_data_full[:,3]
                apu_datetime = sim.jul2datetime(apu_year, apu_Jday, apu_hr, apu_min)

                # Determine if main_platform's timestamp +/- halftime_interval times are avail:
                #  Loop thru ea time in the interval & look for an APU time to match
                apu_in_column_test_cnt = 0
                for t, this_date in enumerate(interval_datetime):
                    time_sub = np.where(apu_datetime == this_date)[0]

                    if len(time_sub) < 1: 
                        continue  # this time not avail in apu files
                    elif len(time_sub) > 1:
                        #this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                        print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date}')
                        print(f'  --- in APU file: {os.path.basename(rain_ter_file)}')
                    #if len(time_sub) == 1: #-- only other possibility at this point so don't really need another test

                    # Location is OK, and current interval time is OK - so proceed:
                    #  Increment counter, set attributes if haven't yet
                    apu_in_column_test_cnt +=1
                    if apu_in_column_test_cnt == 1:
                        # For 1st iteration w/ time & location OK, Set location, plat name/type, & time/offset attributes:
                        apu_lat_d_array.append(apu_plat_info['lat_d'])
                        apu_lat_m_array.append(apu_plat_info['lat_m'])
                        apu_lat_s_array.append(apu_plat_info['lat_s'])
                        apu_lon_d_array.append(apu_plat_info['lon_d'])
                        apu_lon_m_array.append(apu_plat_info['lon_m'])
                        apu_lon_s_array.append(apu_plat_info['lon_s'])
                        apu_plat_name_array = np.append(apu_plat_name_array, apu_plat_info['plat_name'])
                        apu_plat_type_array.append(apu_plat_info['plat_type'])
                        apu_file_type_array.append(apu_plat_info['operation_mode'])
                        apu_int_width_array.append(box_params.halftime_interval)
                        # just will set main_plat time w/ sec as '00' since using an interval of data
                        #apu_plat_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                        apu_plat_timestamp_array.append(column.time)
                        apu_offset_array.append(0)
                    print(f'     at timestamp: {this_date.strftime("%m/%d/%Y %H%M%S")}')

                    # Find column grid point closest to the APU -> place the APU values at this point
                    #  -> arrays of data from the APU files have only [FIELD, TIME] dimensions
                    #	   --> FIELD number:  which column in the text file to use (see APUs_FilesNotes.txt)
                    #	   --> TIME number:  use the subscript for the main_plat_timestamp
                    #  -> arrays that am putting the values from the APU files into are set with same
                    #		dims as the column box grid [lon, lat, alt] --> now as [time, lon, lat, alt]
                    #  NOTE: original/exact APU locations are maintained in apu_plat_info for attributes
                    closest_col_lat_sub = sim.closest(box_params.lat_values, apu_lat)
                    closest_col_lon_sub = sim.closest(box_params.lon_values, apu_lon)
                    if len(closest_col_lat_sub) != 1 or len(closest_col_lon_sub) != 1:
                        print('-----PROBLEM SETTING WHERE IN GRID TO PLACE APU VALUES!!--')
                        print(f'---Returning without parsivel data in the column---')
                        return column
                    else:
                        closest_col_lat_sub = closest_col_lat_sub[0]
                        closest_col_lon_sub = closest_col_lon_sub[0]
                    #pdb.set_trace()
                    # set values for APU fields based on APU file type:
                    
                    #the units are common among the files so will define only once here
                    apu_units_n_drops = 'number of drops'
                    apu_units_concentration = 'drops m^-3 of air'
                    apu_units_liqwatercontent = 'g m^-3'
                    apu_units_rainrate = 'mm h^-1'
                    apu_units_ref_inRayleigh = 'dBZ'
                    apu_units_massweight_diam = 'mm'
                    apu_units_maximum_diam = 'mm'
                    apu_units_bin_drop_concen = 'drops m^-3 mm^-1'
                    apu_units_all_refs = 'dBZ'
                    apu_units_all_attens = 'dB km^-1'
                    units_binDiameter_andWidth = 'mm'
                    #print(self.apu_rain_file)
                    #print(self.apu_rain_ter_file)
                    if apu_rain_file:
                        apu_n_drops[0, closest_col_lat_sub, closest_col_lon_sub, t]            = rain_data_full[time_sub,4]
                        apu_concentration[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_data_full[time_sub,5]
                        apu_liqwater[0, closest_col_lat_sub, closest_col_lon_sub, t]           = rain_data_full[time_sub,6]
                        apu_rain_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_data_full[time_sub,7]
                        apu_rain_refinRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_data_full[time_sub,8]
                        apu_massweight_diam[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_data_full[time_sub,9]
                        apu_maximum_diam[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_data_full[time_sub,10]
                    if apu_rain_ter_file:
                        #pdb.set_trace()
                        apu_n_drops_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]            = rain_ter_data_full[time_sub,4]
                        apu_concentration_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_ter_data_full[time_sub,5]
                        apu_liqwater_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]           = rain_ter_data_full[time_sub,6]
                        apu_rain_rainrate_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_ter_data_full[time_sub,7]
                        apu_rain_refinRayleigh_ter[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_ter_data_full[time_sub,8]
                        apu_massweight_diam_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_ter_data_full[time_sub,9]
                        apu_maximum_diam_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_ter_data_full[time_sub,10]
                    if apu_dsd_file:
                        for binn in range(apu_drop_concen.shape[0]):
                            apu_drop_concen[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = dsd_data_full[time_sub, binn+4]
                    if apu_dsd_ter_file:
                        for binn in range(apu_drop_concen_ter.shape[0]):
                            apu_drop_concen_ter[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = dsd_ter_data_full[time_sub, binn+4]
                    #if apu_ref_file:
                    #    apu_ref_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]       = ref_data_full[time_sub,4]
                    #    apu_ref_refinRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]  = ref_data_full[time_sub,5]
                    #    apu_ref_atS[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,6]
                    #    apu_ref_atC[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,7]
                    #    apu_ref_atX[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,8]
                    #    apu_ref_atKu[0, closest_col_lat_sub, closest_col_lon_sub, t]           = ref_data_full[time_sub,9]
                    #    apu_ref_atK[0, closest_col_lat_sub, closest_col_lon_sub, t]             = ref_data_full[time_sub,10]
                    #    apu_ref_atKa[0, closest_col_lat_sub, closest_col_lon_sub, t]           = ref_data_full[time_sub,11]
                    #    apu_ref_atW[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,12]
                    #if apu_atten_file:
                    #    apu_atten_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]     = atten_data_full[time_sub,4]
                    #    apu_atten_atS[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,5]
                    #    apu_atten_atC[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,6]
                    #    apu_atten_atX[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,7]
                    #    apu_atten_atKu[0, closest_col_lat_sub, closest_col_lon_sub, t]         = atten_data_full[time_sub,8]
                    #    apu_atten_atK[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,9]
                    #    apu_atten_atKa[0, closest_col_lat_sub, closest_col_lon_sub, t]         = atten_data_full[time_sub,10]
                    #    apu_atten_atW[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,11]
            else: # lat & lon are not OK
                print(f'    APU unit not in column box: {apu_plat_info["plat_name"]}')

    # If there are APUs in the box, must now remove dummy element at start of holder arrays
    if len(apu_lat_d_array) > 1:
        apu_lat_d_array = apu_lat_d_array[1:]
        apu_lat_m_array = apu_lat_m_array[1:]
        apu_lat_s_array = apu_lat_s_array[1:]
        apu_lon_d_array = apu_lon_d_array[1:]
        apu_lon_m_array = apu_lon_m_array[1:]
        apu_lon_s_array = apu_lon_s_array[1:]
        apu_plat_name_array = apu_plat_name_array[1:]
        apu_plat_type_array = apu_plat_type_array[1:]
        apu_file_type_array = apu_file_type_array[1:]
        apu_int_width_array = apu_int_width_array[1:]
        apu_plat_timestamp_array = apu_plat_timestamp_array[1:]
        apu_offset_array = apu_offset_array[1:]
    else:
        # If there are NO APUs in the box, the data arrays will still be all NaNs
        # and will just return a -1 to the top level .pro
        print(f'---There are no parsivels in the column box---')
        print(f'---Returning without parsivel data in the column---')
        return column
    print(f'------APUs in the column box: {apu_plat_name_array}')

    # Populate arrays for DSD bin middle/width diameters, if have DSD files:
    if n_dsd >= 1 or n_dsd_ter >= 1:
        # there should be a parsivel_diameter.txt file in the apu dir if get here:
        diameters_file = apu_dir+'/parsivel_diameter.txt'
        diameters_data     = np.loadtxt(diameters_file)
        bin_middles      = diameters_data[:,0]
        bin_widths       = diameters_data[:,1]
        apu_dsd_bins_Dmiddle = bin_middles.flatten()
        apu_dsd_bins_widths  = bin_widths.flatten()

    # Set up dictionary of the plat_infos for APUs within the column box:
    #  need this to be the "plat_info" values for all APUs within box
    apu_info = {'plat_name':apu_plat_name_array,
                'lat_d':apu_lat_d_array,'lat_m':apu_lat_m_array,'lat_s':apu_lat_s_array,
                'lon_d':apu_lon_d_array,'lon_m':apu_lon_m_array,'lon_s':apu_lon_s_array,
                'plat_type':apu_plat_type_array,
                'operation_mode':apu_file_type_array,
                'timestamp':apu_plat_timestamp_array,
                'offset_vs_main':apu_offset_array,
                'time_interval_width':apu_int_width_array}

    #add apu_info to object
    column.add_platform_to_object(apu_info, plat_name='parsivel')

    #assign each field to column object
    if apu_rain_file:
        column.add_variable_to_object(apu_n_drops, 
                                      var_name='apu_total_n_drops', 
                                      units=apu_units_n_drops)

        column.add_variable_to_object(apu_concentration, 
                                      var_name='apu_total_concentration', 
                                      units=apu_units_concentration)

        column.add_variable_to_object(apu_liqwater, 
                                      var_name='apu_liquid_water_content', 
                                      units=apu_units_liqwatercontent)

        column.add_variable_to_object(apu_rain_rainrate, 
                                      var_name='apu_rain_rate', 
                                      units=apu_units_rainrate)

        column.add_variable_to_object(apu_rain_refinRayleigh, 
                                      var_name='apu_reflectivity_in_Rayleigh', 
                                      units=apu_units_ref_inRayleigh)

        column.add_variable_to_object(apu_massweight_diam, 
                                      var_name='apu_mass_weighted_drop_diameter', 
                                      units=apu_units_massweight_diam)

        column.add_variable_to_object(apu_maximum_diam, 
                                      var_name='apu_maximum_drop_diameter', 
                                      units=apu_units_maximum_diam)

    if apu_rain_ter_file:
        column.add_variable_to_object(apu_concentration_ter, 
                                      var_name='apu_total_concentration_tfs', 
                                      units=apu_units_concentration)

        column.add_variable_to_object(apu_liqwater_ter, 
                                      var_name='apu_liquid_water_content_tfs', 
                                      units=apu_units_liqwatercontent)

        column.add_variable_to_object(apu_rain_rainrate_ter, 
                                      var_name='apu_rain_rate_tfs', 
                                      units=apu_units_rainrate)

        column.add_variable_to_object(apu_rain_refinRayleigh_ter, 
                                      var_name='apu_reflectivity_in_Rayleigh_tfs', 
                                      units=apu_units_ref_inRayleigh)

        column.add_variable_to_object(apu_massweight_diam_ter, 
                                      var_name='apu_mass_weighted_drop_diameter_tfs', 
                                      units=apu_units_massweight_diam)

        column.add_variable_to_object(apu_maximum_diam_ter, 
                                      var_name='apu_maximum_drop_diameter_tfs', 
                                      units=apu_units_maximum_diam)

    if apu_dsd_file:
        column.add_variable_to_object(apu_drop_concen, 
                                      var_name='apu_dsd_concentration', 
                                      units=apu_units_bin_drop_concen)

        column.add_variable_to_object(apu_dsd_bins_Dmiddle, 
                                      var_name='apu_dsd_bin_median_diameter', 
                                      units=units_binDiameter_andWidth)

        column.add_variable_to_object(apu_dsd_bins_widths, 
                                      var_name='apu_dsd_bin_width', 
                                      units=units_binDiameter_andWidth)

    if apu_dsd_ter_file:
        column.add_variable_to_object(apu_drop_concen_ter, 
                                      var_name='apu_dsd_concentration_tfs', 
                                      units=apu_units_bin_drop_concen)

        column.add_variable_to_object(apu_dsd_bins_Dmiddle, 
                                      var_name='apu_dsd_bin_median_diameter_tfs', 
                                      units=units_binDiameter_andWidth)

        column.add_variable_to_object(apu_dsd_bins_widths, 
                                      var_name='apu_dsd_bin_width_tfs', 
                                      units=units_binDiameter_andWidth)
    
    return column
    #self.apu_data_in_column = True
