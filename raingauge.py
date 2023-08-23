import os
import numpy as np
from datetime import datetime
import pdb #pdb.set_trace()
import tools as sim

def get_gauges_for_column(gauge_gmin_files, interval_datetime, box_params, column):

    #    interval_datetime[i] = main_datetime + timedelta(minutes=step.item())
    #print(interval_datetime[0].strftime(f'%Y%m%d%H%M'))
    #print(interval_datetime[self.halftime_interval].strftime(f'%Y%m%d%H%M'))
    #print(interval_datetime[-1].strftime(f'%Y%m%d%H%M'))
    #print('TESTING ^^^ -- comment out')

    gauge_lat_d_array = [-9999]
    gauge_lat_m_array = [-9999]
    gauge_lat_s_array = [-9999]
    gauge_lon_d_array = [-9999]
    gauge_lon_m_array = [-9999]
    gauge_lon_s_array = [-9999]
    gauge_plat_name_array = ['DUMMY']
    gauge_plat_type_array = ['DUMMY']
    gauge_file_type_array = ['DUMMY']
    gauge_plat_timestamp_array = ['DUMMY']
    gauge_offset_array = [-9999]
    gauge_int_width_array = [-9999]
    # gauge_rainrate data array:
    
    gauge_rainrate = np.full(box_params.ground_plat_dims, np.nan)
    gauge_units_rainrate = 'mm h^-1'

    # Loop thru available .gmin files:
    #  Steps must do for each .gmin file:
    #  1) Find if the main_plat time is included in the record
    #  2) If time is ok, determine if location is in col box
    #  3) If time is ok and location is in col box, 2 options:
    #     a) if col box pt nearest the gauge still NaN, populate
    #   the col box pt w/ current gauge data
    #     b) if col box pt nearest the gauge ALREADY NOT NAN,
    #   then this is an _A/_B gauge situation w/ paired tip
    #   buckets OR have multiple gauges near same grid pt,
    #   so update the col box pt w/ AVERAGE of the two vals
    for i, gmin_file in enumerate(gauge_gmin_files):
        print(f'---- Gauge GMIN File: {os.path.basename(gmin_file)}')
        #test=file_basename(gmin_file)
        #if strmid(test, 9,6) EQ '0007_B' then stop
        #if strmid(test, 9,4) EQ 'PAD6' then stop
        # Note on time fields: Jerry's files include leading 0 if <10: Mon, Day, Jday, Hr, Min
        #columns 'year', 'month', 'day', 'Jday', 'hour', 'minute', 'rainrate', 'lat', 'lon'
        gmin_data = np.genfromtxt(gmin_file, skip_header=2) #will return float 2d array
        #gmin_header = gmin_header[0]
        gmin_year = gmin_data[:,0]
        gmin_mon  = gmin_data[:,1]
        gmin_day  = gmin_data[:,2]
        gmin_Jday = gmin_data[:,3]
        gmin_hr   = gmin_data[:,4]
        gmin_min  = gmin_data[:,5]
        gmin_rain = gmin_data[:,6]
        gmin_lat  = gmin_data[:,7]
        gmin_lon  = gmin_data[:,8]

        # 1) Determine if the main_plat_timestamp +/- haftime_interval times are included in the gmin_file:
        #    Loop through each time in the interval and look for a gmin time to match.
        #    If find a match, proceed to Step 2. (Must do (1) first bc locations in .gmin files are based on the time)
        gauge_in_column_test_cnt = 0
        gmin_datetime = sim.jul2datetime(gmin_year, gmin_Jday, gmin_hr, gmin_min)
        
        for t, this_date in enumerate(interval_datetime):
            time_sub = np.where(gmin_datetime == this_date)[0]
      
            if len(time_sub) < 1:
                continue # this time is not available in the gague file
            elif len(time_sub) > 1:
                #this shouldn't happen bc Jerry's GMIN files are at 1 min interval
                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} {this_date.strftime("%M")}')
                print(f'  --- in gmin file: {os.path.basename(gmin_file)}')
            #else:#len(time_sub) == 1    
            # 2) Have the time available in the file. Now check if location is within column box:
            lat_ok = 'N'
            lon_ok = 'N'
            if gmin_lat[time_sub[0]] <= box_params.box_max_lat  and  gmin_lat[time_sub[0]] >= box_params.box_min_lat: lat_ok = 'Y'
            if gmin_lon[time_sub[0]] <= box_params.box_max_lon  and  gmin_lon[time_sub[0]] >= box_params.box_min_lon: lon_ok = 'Y'
            
            if lat_ok == 'Y' and lon_ok == 'Y':
                # Only continue here if the lat & lon are in the column
                gauge_in_column_test_cnt += 1
                 
                # Find column grid point closest to this gauge location:
                # Convert the lat/lons from dec deg to DMS:
                #   will need to plat_info, and closest function works
                #   better after dec deg from dms (not sure why? data type issue?)
                gmin_lat_d, gmin_lat_m, gmin_lat_s, gmin_lon_d, gmin_lon_m, gmin_lon_s = sim.dd2dms(gmin_lat[time_sub[0]], gmin_lon[time_sub[0]])
                unit_location_decdeg = sim.dms2dd(gmin_lat_d, gmin_lat_m, gmin_lat_s, gmin_lon_d, gmin_lon_m, gmin_lon_s)
                unit_lat = unit_location_decdeg[0]
                unit_lon = unit_location_decdeg[1]
                if gauge_in_column_test_cnt == 1:
                    # if 1st time have this gauge avail in columne grid, print first line to terminal & set up for attributes
                    print(f'    Time & Loc OK For: {os.path.basename(gmin_file)}')
                    
                    f_name = os.path.basename(gmin_file)
                    this_gauge_name = f_name[0:f_name.rfind('-')]
                    this_gauge_type = 'gauge'
                    this_file_type  = 'gmin_gauge'
                    gauge_lat_d_array.append(gmin_lat_d)
                    gauge_lat_m_array.append(gmin_lat_m)
                    gauge_lat_s_array.append(gmin_lat_s)
                    gauge_lon_d_array.append(gmin_lon_d)
                    gauge_lon_m_array.append(gmin_lon_m)
                    gauge_lon_s_array.append(gmin_lon_s)
                    gauge_plat_name_array.append(this_gauge_name)
                    gauge_plat_type_array.append(this_gauge_type)
                    gauge_file_type_array.append(this_file_type)
                    gauge_int_width_array.append(box_params.halftime_interval)
                    # just will set main_plat time w/ sec as '00' since using an interval of data
                    #gauge_plat_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                    gauge_plat_timestamp_array.append(column.time)
                    gauge_offset_array.append(0)
                print(f'         at timestamp: {this_date.strftime("%m/%d/%Y %H%M%S")}')
                closest_col_lat_sub = sim.closest(box_params.lat_values, unit_lat)
                closest_col_lon_sub = sim.closest(box_params.lon_values, unit_lon)
                if len(closest_col_lat_sub) != 1 or len(closest_col_lon_sub) != 1:
                    print(f'---- PROBLEM SETTING WHERE IN GRID TO PLACE GAGUE: {os.path.basename(gmin_file)}')
                    continue
                else:
                    closest_col_lat_sub = closest_col_lat_sub[0]
                    closest_col_lon_sub = closest_col_lon_sub[0]

                # Determine if the current value at the grid spot nearest this gauge for this time is already populated:
                #  use a str/char variable called "exists" to tell if have a previous value:  Y or N only
                if ~ np.isfinite(gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t]): exists = 'N'
                if np.isfinite(gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t]):   exists = 'Y'

                # 3a)  If no previous value, populate nearest grid point for this time w current gauge data
                if exists == 'N':
                    gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t] = gmin_rain[time_sub[0]]

                # 3b)  If there IS a previous value, take the max between the previous and current value
                #        This will be the case for situations with Iowa-type/paired tip bucket gauges
                #        or whenever multiple stand-alone gauges at same location (eg: WFF N-159 pad)
                if exists == 'Y':
                    prev_val = gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t]  #[t] is sub into interval times for new data holder array
                    this_val = gmin_rain[time_sub[0]]           #[time_sub] is sub into the arrays read in from the gmin files, corresponds to the current interval_*[t]
                    new_val = max([prev_val, this_val]) #take the maximum rate instead of averaging -- updated by Cpabla 08/13/2019
                    #new_val  = mean([prev_val, this_val]) ;previously took the mean between previous and current value
                    gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t] = new_val

            else: #test for location in the column box
                # this gauge is not in the column box
                print(f'    gauge not in column box: {os.path.basename(gmin_file)}')

    # After loop thru gauge files: clean up info arrays, prep structure for return
    if len(gauge_lat_d_array) > 1:
        gauge_lat_d_array = gauge_lat_d_array[1:]
        gauge_lat_m_array = gauge_lat_m_array[1:]
        gauge_lat_s_array = gauge_lat_s_array[1:]
        gauge_lon_d_array = gauge_lon_d_array[1:]
        gauge_lon_m_array = gauge_lon_m_array[1:]
        gauge_lon_s_array = gauge_lon_s_array[1:]
        gauge_plat_name_array = gauge_plat_name_array[1:]
        gauge_plat_type_array = gauge_plat_type_array[1:]
        gauge_file_type_array = gauge_file_type_array[1:]
        gauge_plat_timestamp_array = gauge_plat_timestamp_array[1:]
        gauge_offset_array = gauge_offset_array[1:]
        gauge_time_interval_width = gauge_int_width_array[1:] #new attribute with addition of time interval setting
    else:
        #don't have any gauges in column box and/or main_plat time unavail
        print(f'--- NO GAUGES IN COLUMN BOX AND/OR DURING MAIN PLAT TIME: {column.time} ---')
        print('---Returning without Gauge data in the column---')
        return column

    #Set up dictionary
    gauges_info = {'plat_name':gauge_plat_name_array,
                   'lat_d':gauge_lat_d_array,'lat_m':gauge_lat_m_array,'lat_s':gauge_lat_s_array,
                   'lon_d':gauge_lon_d_array,'lon_m':gauge_lon_m_array,'lon_s':gauge_lon_s_array,
                   'plat_type':gauge_plat_type_array,
                   'operation_mode':gauge_file_type_array,
                   'timestamp':gauge_plat_timestamp_array,
                   'offset_vs_main':gauge_offset_array,
                   'time_interval_width':gauge_time_interval_width}

    #add gauges_info to object
    column.add_platform_to_object(gauges_info, plat_name='gauges')

    #assign field to column object
    column.add_variable_to_object(gauge_rainrate, 
                                    var_name='gauges_rainrate', 
                                    units=gauge_units_rainrate)
    return column
    #self.gauges_data_in_column = True