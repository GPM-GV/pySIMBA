import os, sys
import argparse
import ast
from pathlib import Path
import core
import tools as sim
import write_column
import disdrometer
import gpm
import ground_radar
import mrms
import mrr
import raingauge
import warnings
warnings.filterwarnings("ignore")

def BuildColumn(year, mon, day):

    #make sure year mon and day are strings
    year = str(year).zfill(4)
    mon = str(mon).zfill(2)
    day = str(day).zfill(2)
    
    #check parameter file
    cwd = os.getcwd()
    params_file = cwd+'/pysimba_params_dict.txt'
    if Path(params_file).is_file():
        kwargs = ast.literal_eval(open(params_file).read())
    else:
        sys.exit(f'Required parameter file does not exist {params_file}!')
    
    #Extract initial parameters
    main_plat = kwargs['Main_Plat_Radar']
    center_on = kwargs['Box_center']
    VN_radars = kwargs['Secondary_Radars']
    box_spacing = kwargs['Box_spacing']
    box_limit = kwargs['Box_limit']
    vert_spacing = kwargs['Vertical_spacing']
    vert_limit = kwargs['Vertical_limit']
    halftime_interval = kwargs['Ground_obs_temporal_interval']
    out_dir = kwargs['Output_dir']
    data_dir = kwargs['Data_dir']
    out_netcdf = kwargs['Out_netcdf']
    
    #-----------
    #  step (1):  Set paths for searching platform data and output
    main_plat_file, VN_IDS, VN_radar_files, gauge_gmin_files, apu_files, \
    twodvd_files, mrr_files, dpr_file, mrms_files, nc_out_dir = core.search_data(year, mon, day, out_dir, 
                                                        data_dir, main_plat, VN_radars)

    #-----------
    #  step (2):  define grid points for the column
    #  based on center location, spacing, limit set at top
    column_info = core.define_box(center_on, box_spacing, box_limit, vert_spacing, vert_limit, halftime_interval)
    if out_netcdf: column_info.print_box_parameters()

    #-----------
    #  step (3):  get the main_plat information and set timelag interval 
    #  NOTE: currently, main_plat must be a ground-based scanning radar!
    #  main_plat's timestamp will be the new .nc file's timestamp
    #  time lag will be set based on user provided halftime and main platform timestamp
    main_plat_params = core.setup_main_plat(main_plat, main_plat_file)
    main_plat_datetime = main_plat_params['datetime']
    main_scan_type = main_plat_params['scan_type']
    main_lat_deg = main_plat_params['lat_deg']
    main_lon_deg = main_plat_params['lon_deg']
    time_lag_datetime = core.get_time_lag(halftime_interval, main_plat_datetime)
    
    #-----------
    # step (4):  define object to store all platform and associated data, 
    time = main_plat_datetime.strftime('%Y%m%d_%H%M%S')
    variables={}
    platforms={}
    column = core.GridColumn(time, variables, platforms, column_info)
    
    #-----------
    # step (5): grid main_plat radar data to the column box framework 
    # and add fields to column object
    print('-----------------------------------------------')
    print(f'-----in Main Platform {main_plat} section-----')
    print(f'main_plat timestamp: {time}')
    print(f'main_plat scan type: {main_scan_type}')
    column = ground_radar.grid_radar_for_column(main_plat_file, main_plat_datetime, main_plat, column_info, column)
    
    #-----------
    # step (6):  if there are other lev2 radar data available, grid to the column box framework
    #  and add fields to column object
    for i, r_file in enumerate(VN_radar_files):
        radarID = VN_IDS[i]
        print(f'-----in lev2 {radarID} section-----')
        column = ground_radar.grid_radar_for_column(r_file, main_plat_datetime, radarID, column_info, column)
    #pdb.set_trace()
    column = sim.combine_radar_info(column) #combine radar info attributes
    
    #-----------
    # step (7):  if there is PARSIVEL data, read it and add fields to column object
    if apu_files:
        print('-----in PARSIVEL section-----')
        column = disdrometer.get_parsivel_for_column(apu_files, time_lag_datetime, column_info, column)
    
    #-----------
    # step (8):  if there is 2DVD data, read it and add fields to column object
    if twodvd_files:
        print('-----in 2DVD section-----')
        column = disdrometer.get_2dvd_for_column(twodvd_files, time_lag_datetime, column_info, column)
    
    #-----------
    # step (9):  if there is .gmin GAUGE data, read it and add fields to column object
    if gauge_gmin_files:
        print('-----in GAUGES section-----')
        column = raingauge.get_gauges_for_column(gauge_gmin_files, time_lag_datetime, column_info, column)
    
    #-----------
    # step (10):  if there is MRMS Lev2 data, read it and add fields to column object
    if mrms_files:
        print('-----in MRMS section-----')
        column = mrms.get_mrms_for_column(mrms_files, main_plat_datetime, column_info, column)
    
    #-----------
    # step (11):  if there is MRR data, read it and add fields to column object
    if mrr_files:
        print('-----in MRR section-----')
        column = mrr.get_mrr_for_column(mrr_files, time_lag_datetime, column_info, column)
    
    #-----------
    # step (12):  if there is GPM 2ADPR data, read it and add fields to column object
    if dpr_file:
        print('-----in DPR section-----')
        column = gpm.get_dpr_for_column(dpr_file, main_plat_datetime, main_lat_deg, main_lon_deg, column_info, column)

    #-----------
    # step (13):  create a netCDF file containing all data within the column box
    #   and include as much info as possible in the attributes. save the new file
    #    to the out_dir [if user requests]
    #  COLUMN netCDF FILES WILL BE NAMED AS:
    #       column_[main]_[center]_YYYYMMDD_HHMM.nc
    #           [main]:    string ID of main_plat used
    #           [center]:  string ID of box's center location
    #           timestamp: main_plat's time w/o seconds
    #         example:  column_NPOL_WFFpad_20160628_1449.nc
    if out_netcdf:
        print('-----in netCDF section-----')
        write_column.netcdf(main_plat, nc_out_dir, column_info, column)
    else:
        print('--------------------------------')
        print('Column object created!', '', sep='\n')
        return column

# *******************************************  M  A  I  N  ******************************

def main(inargs):
    """Run the program."""
    
    #extract input arguments
    year = inargs.year.zfill(4)
    month = inargs.month.zfill(2)
    day = inargs.day.zfill(2)

    #run BuildColumn
    BuildColumn(year, month, day)
    #c = BuildColumn(year, month, day, **kwargs)
    #column = c.process_columns()

    print('Program Done.', '', sep='\n')

if(__name__ == '__main__'):
            
    runargs = argparse.ArgumentParser(description='Produce precipitating columns over a user defined geographic location.')
    runargs.add_argument('year',  type=str, help='Provide the year in YYYY format')
    runargs.add_argument('month', type=str, help='Provide the month in MM format')
    runargs.add_argument('day',   type=str, help='Provide the day in DD format')
    #runargs.add_argument('--params_dict', dest='params_dict', type=str, help='Parameter dictionary')

    args = runargs.parse_args()
    main(args)
