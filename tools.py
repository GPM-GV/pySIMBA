import sys
import numpy as np
import pandas as pd
import math
import pyart
from scipy import interpolate
from datetime import datetime
import os, shutil, gzip
import pdb

def get_posn_to_azm_range(lat, lon, radarlat, radarlon):
    ''' 
    #-------------------------------------------------------------------------------
    # function to compute the azimuth and range, given lat, lon, radarlat, radarlon
    #-------------------------------------------------------------------------------
    #   result = get_posn_to_azm_range(lat, lon, radarlat, radarlon)
    #   lat:        float value of latitude in dec deg
    #   lon:        float value of longitude in dec deg
    #   radarlat:   float value of latitude in dec deg of the radar
    #   radarlon:   float value of longitude in dec deg of the radar
    #   
    #   result:  2-element array:  [azimuth, range]
    #           result[0]:  azimuth that point is away from the radar
    #           result[1]:  dist in [m] that point is away from radar
    #
    # Converted/Modified by: Charanjit S. Pabla
    # Written by: Stephanie M. Wingo
    # Original Owner David Marks
    #-------------------------------------------------------------------------------
    '''
    lat_to_meters = 111177.
    lon_to_meters = 111177. * np.cos(np.deg2rad(radarlat))

    x = (lon - radarlon) * lon_to_meters
    y = (lat - radarlat) * lat_to_meters
    if x != 0: 
        azimuth = np.arctan2(y,x) * (180./np.pi)
    elif y > 0: 
        azimuth = 90.
    else:
        azimuth = -90.

    # Convert standard angle to clockwise rotation for radar, where north is 0.
    #   Quadrant 1, 3, and 4:  90 - a
    #   Quadrant 2:            90 - a + 360
    azimuth = 90. - azimuth
    if y > 0 and x < 0: azimuth = azimuth + 360.
    range = np.sqrt(x**2. + y**2.) # range in meters
    
    return [azimuth, range]

# ***************************************************************************************

def interp(vals_to_interp, input_abscissa, output_vals):
    '''
    #-------------------------------------------------------------------------------
    # function to compute interpolation
    #-------------------------------------------------------------------------------
    #   result = interp(vals_to_interp, input_abscissa, output_vals)
    #   vals_to_interp: float array of values to interpolate
    #   input_abscissa: float array of original levels to interpolate from
    #   output_vals:    float array of levels to interpolate on
    #
    #   result: will be of same dimensions as the input array with interpolated values
    #
    # Written by: Charanjit S. Pabla
    #-------------------------------------------------------------------------------
    '''
    interfunc   = interpolate.interp1d(input_abscissa, vals_to_interp, fill_value='extrapolate', kind='linear')
    #interfunc   = interpolate.interp1d(input_abscissa[~np.isnan(vals_to_interp)], 
    #                    vals_to_interp[~np.isnan(vals_to_interp)],fill_value='extrapolate')
    return interfunc(output_vals)

# ***************************************************************************************

def closest(arr, val):
    '''
    #-------------------------------------------------------------------------------
    # function to compute the index in arr where the minimum occurs given val
    #-------------------------------------------------------------------------------
    #   result = closest(arr, val)
    #   arr: float array (any dimensions)
    #   val: float element
    #
    #   result: index of the minimum element from the array
    #
    # Written by: Charanjit S. Pabla
    #-------------------------------------------------------------------------------
    '''
    #calculate the difference array
    diff_arr = np.absolute(arr-val)

    #return index of the minimum element from the array
    return [diff_arr.argmin()]

# ***************************************************************************************

def datestr2datetime(timestamp):
    '''
    #-------------------------------------------------------------------------------
    # function to convert timestamp to datetime object
    #-------------------------------------------------------------------------------
    #   result = datestr2datetime(timestamp)
    #   timestamp: string array (any dimensions)
    #
    #   result: will be of same dimensions as the input array but elements will
    #           be datetime objects array
    #
    # Written by: Charanjit S. Pabla
    #-------------------------------------------------------------------------------
    '''
    datetimes =  np.empty(timestamp.shape[0], dtype=object)
    fmt='%y%m%d%H%M'
    for i, val in enumerate(timestamp):
        yy=val[0:2]
        mm=val[2:4]
        dd=val[4:6]
        hr=val[6:8]
        mi=val[8:10]
        sec=val[10:] #dont need to use seconds since we are after 1 min samples
        datetimes[i] = datetime.strptime(f'{yy}{mm}{dd}{hr}{mi}', fmt)
    return datetimes

# ***************************************************************************************

def jul2datetime(year, julday, hr, minute):
    '''
    #-------------------------------------------------------------------------------
    # function to convert julday, HH, MM to datetime object
    #-------------------------------------------------------------------------------
    #   result = jul2datetime(year, julday, hr, minute)
    #   year:   float array (any dimensions)
    #   julday: float array (any dimensions)
    #   hr:     float array (any dimensions)
    #   minute: float array (any dimensions)
    #       all input parameters must be the same dimensions
    #
    # result: will be of same dimensions as the input array(s) but elements will
    #         be datetime objects array
    #
    # Written by: Charanjit S. Pabla
    #-------------------------------------------------------------------------------
    '''
    datetimes =  np.empty(year.shape[0], dtype=object)
    
    #convert to int then string
    year = year.astype(int).astype(str) 
    julday = julday.astype(int).astype(str)
    hr = hr.astype(int).astype(str)
    minute = minute.astype(int).astype(str)
    fmt='%y%j%H%M'
    for i, val in enumerate(year):
        datetimes[i] = datetime.strptime(f'{year[i][-2:]}{julday[i].zfill(3)}{hr[i].zfill(2)}{minute[i].zfill(2)}', fmt)
    return datetimes

# ***************************************************************************************

def bad2nan(data, bad_val):
    '''
    #-------------------------------------------------------------------------------
    # function to set all values in an array that match the bad value flag to NaNs
    #-------------------------------------------------------------------------------
    #   result = bad2nan(data, bad_val)
    #   data:    an array (any dimensions and any type) of data with bad/missing data
    #                   set to a flagged value
    #   bad_val: value (any type) of the bad data flag value
    # 
    #   result: will be of same dimensions as the input data, but all elements
    #       equal to the bad_val will now be NaNs
    #
    # Converted/Modified by: Charanjit S. Pabla
    # Written by: Stephanie M. Wingo
    #-------------------------------------------------------------------------------
    '''
    data=data
    bad_val=bad_val

    new_data=data
    #pdb.set_trace()
    bad_data=np.where(data == bad_val)[0]
    if bad_data.shape[0] != 0: new_data[bad_data] = np.nan

    return new_data

# ***************************************************************************************

def ungzip_file(cfile):
    '''
    #-------------------------------------------------------------------------------
    # function to unzip a given .gz file
    #-------------------------------------------------------------------------------
    #   result = unzip_file(cfile)
    #   cfile: original .gz file to unzip
    #
    #   result: unzipped file
    #-------------------------------------------------------------------------------
    '''
    fileb = os.path.basename(cfile)[0:-3]
    with gzip.open(cfile, 'rb') as f_in:
        with open(fileb, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return fileb

# ***************************************************************************************

def dd2dms(lat_decdeg, lon_decdeg):
    '''
    #-------------------------------------------------------------------------------
    # function to get degrees, minutes, seconds for lat/lon point given decimal deg
    #-------------------------------------------------------------------------------
    #   result = dd2dms(lat_decdeg, lon_decdeg)
    #   lat_decdeg: float value of decimal degrees
    #   lon_decdeg: float value of decimal degrees
    #           sign of the degrees tells if N/S or E/W
    #   
    #   result:  returns a 1-d, 6 element float array of deg, min, sec for lat/lon:
    #       [lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec]
    #       sign on the lat_deg and lon_deg elemtns tells if N/S or E/W
    #                ex:  if 35deg N 75deg W, will get:
    #                   lat_deg = 35
    #                   lon_deg = -75
    #
    # Converted/Modified by: Charanjit S. Pabla
    # Written by: Stephanie M. Wingo
    #-------------------------------------------------------------------------------
    '''
    #get sign for lat & lon:
    lat_sign = lat_decdeg / abs(lat_decdeg)
    lon_sign = lon_decdeg / abs(lon_decdeg)

    # process:  deg = integer part of dec_deg
    #           min = integer part of ( 60.0D *(dec_deg - deg) )
    #           sec = 3600.0*(dec_deg - deg - (min/60.0) )

    # get deg part:
    lat_deg = math.floor(abs(lat_decdeg))
    lon_deg = math.floor(abs(lon_decdeg))

    # get min part:
    lat_min_dec = 60.0*(abs(lat_decdeg) - lat_deg)
    lon_min_dec = 60.0*(abs(lon_decdeg) - lon_deg)
    lat_min = math.floor(lat_min_dec)
    lon_min = math.floor(lon_min_dec)

    # get sec part:
    lat_sec = 3600.0*(abs(lat_decdeg) - lat_deg - (lat_min/60.0))
    lon_sec = 3600.0*(abs(lon_decdeg) - lon_deg - (lon_min/60.0))

    # adjust sec part to include only 2 decimal places:
    lat_sec = lat_sec*100.0
    lat_sec = round(lat_sec)
    lat_sec = lat_sec/100.0
    lon_sec = lon_sec*100.0
    lon_sec = round(lon_sec)
    lon_sec = lon_sec/100.0

    # add sign for the degrees
    lat_deg = lat_sign * lat_deg
    lon_deg = lon_sign * lon_deg

    return int(lat_deg), int(lat_min), int(lat_sec), int(lon_deg), int(lon_min), int(lon_sec)

# ***************************************************************************************

def xy2ll(grid_center_lat, grid_center_lon, x_vals, y_vals):
    '''
    #-------------------------------------------------------------------------------
    # function to compute lat/lon for x/y grid coords given the lat/lon of center
    #-------------------------------------------------------------------------------
    #   result = xy2ll(grid_center_lat, grid_center_lon, x_vals, y_vals)
    #   grid_center_lat: float value of decimal latitude of grid center point
    #   grid_center_lon: float value of decimal longitude of grid center point
    #   x_vals:          float value (or array) of x coordinates in UNITS: [km]
    #   y_vals:          float value (or array) of y coordinates in UNITS: [km]
    #
    #   result: returns either a simple array or a dictionary, depending on inputs:
    #       input single x, y point:
    #           returns [lat,lon] values as a 1-D array
    #       input x, y as 1-D arrays:
    #           returns dictionary with keys naming 1-D arrays of
    #           lats, lons that correspond to the input x, y points:
    #              lats = result[0] & lons = result[1]
    #  NOTE:  lats will be  more accurate bc the spacing of lons varries
    #
    # Converted/Modified by: Charanjit S. Pabla
    # Originally Written by: Stephanie M. Wingo
    #-------------------------------------------------------------------------------
    '''
    clat = grid_center_lat
    clon = grid_center_lon
    xcord = x_vals
    ycord= y_vals

    if isinstance(xcord, float) and isinstance(ycord, float):
        x=np.zeros(1)
        y=np.zeros(1)
        
        xstp=0
        ystp=0
    else:
        x=np.zeros(len(xcord))
        y=np.zeros(len(ycord))
        
        xstp=len(xcord)-1
        ystp=len(ycord)-1
        
        if len(xcord) != len(ycord):
            print(' ***xy2ll: INPUT X & Y COORDS MUST HAVE SAME # OF PTS!')
            #'GOTO, END_FUN...NEED Try...Except... HERE'
            pdb.set_trace() #debug mode; hopefully this does not happen

    x=xcord
    y=ycord
        
    #"one more try:"  -->> USE THIS ONE.
    #from DM & example codes: meters_to_lat = 1. / 111177.
    #			    meters_to_lon = 1. / (111177. * cos(radar_lat*!dtor))
    #then when getting coords:
    #lat = radar_lat + range*meters_to_lat
    #lon = radar_lon + range*meters_to_lon
    #use grid center in place of radar (grid in examples were centered on radar)
    meters_to_lat = 1. / 111177.
    meters_to_lon = 1. / (111177. * np.cos(math.radians(clat)))
    #and for this need x & y in [m] not [km]:
    x = x*1000.0
    y = y*1000.0


    if isinstance(xcord, float) and isinstance(ycord, float):
        #have single point (x,y)
        LL=np.zeros(2)
        lat=y*0.0
        lon=x*0.0
    
        #one more try: for both lat & lon
        lat = clat + y*meters_to_lat
        lon = clon + x*meters_to_lon
        
        LL[0]=lat
        LL[1]=lon
    else:
        #have arrays of x & y values
        LL={'lat':np.zeros(ystp+1),'lon':np.zeros(xstp+1)}
        lats=np.zeros(ystp+1)
        lons=np.zeros(xstp+1)
    
        for i in xstp:
            for j in ystp:
                lats[j] = clat + y[j]*meters_to_lat
                lons[i] = clon + x[i]*meters_to_lon
        LL['lat']=lats
        LL['lon']=lons
        
    return LL

# ***************************************************************************************

def dms2dd(lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec):
    '''-------------------------------------------------------------------------------
    # function to get decimal degrees for lat/lon point given in degrees, mins, sec
    #-------------------------------------------------------------------------------
    #   result = dms2dd(lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec)
    #   lat_deg: int/float/double value
    #   lat_min: int/float/double value
    #   lat_sec: int/float/double value
    #   lon_deg: int/float/double value
    #   lon_min: int/float/double value
    #   lon_sec: int/float/double value
    #
    #       lat_deg, lon_deg: sign of the degrees tells if N/S or E/W
    #                ex:  if 35deg N 75deg W, enter:
    #                   lat_deg = 35
    #                   lon_deg = -75
    #   result:  returns a 2 element float Tuple of decimal lat/lon values:
    #        (lat.xxxxxxx, lon.xxxxxxx)
    #
    # Converted/Modified by: Charanjit S. Pabla
    # Originally Written by: Stephanie M. Wingo
    #-------------------------------------------------------------------------------
    '''

    #process: (sign)[ deg + min/60.D + sec/3600.D ]

    #convert seconds:
    lat_dec_sec = lat_sec/3600.0
    lon_dec_sec = lon_sec/3600.0

    #convert minutes:
    lat_dec_min = lat_min/60.0
    lon_dec_min = lon_min/60.0

    #get degrees as positive
    lat_dec_deg = abs(lat_deg)
    lon_dec_deg = abs(lon_deg)

    #get sign from degress:
    lat_sign = lat_deg / abs(lat_deg)
    lon_sign = lon_deg / abs(lon_deg)

    #all together now:
    new_lat = lat_sign * (lat_dec_deg + lat_dec_min + lat_dec_sec)
    new_lon = lon_sign * (lon_dec_deg + lon_dec_min + lon_dec_sec)

    return new_lat, new_lon

# ***************************************************************************************

def get_platform_loc(platform_name):

    #read in platform location text file into pandas
    cwd = os.getcwd()
    platform_file = cwd+'/platform_locations.txt'
    #columns = [platform, lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec, elv (m)]
    cols = ['platform','lat_deg','lat_min','lat_sec','lon_deg','lon_min','lon_sec', 'elv']
    platforms = pd.read_csv(platform_file, sep='\s+', header=None, names=cols, skiprows=1)
    platforms.set_index('platform', inplace = True)

    #verify platform exist in the dataframe
    if platform_name not in platforms.index:
        sys.exit(f'platform {platform_name} does not exist -- check {platform_file}')
    else:
        return platforms.loc[platform_name,'lat_deg'].astype('float'),\
               platforms.loc[platform_name,'lat_min'].astype('float'),\
               platforms.loc[platform_name,'lat_sec'].astype('float'),\
               platforms.loc[platform_name,'lon_deg'].astype('float'),\
               platforms.loc[platform_name,'lon_min'].astype('float'),\
               platforms.loc[platform_name,'lon_sec'].astype('float'),\
               platforms.loc[platform_name, 'elv'].astype('float')

# ***************************************************************************************

def combine_radar_info(column):

    #count number of radars
    n_radars = len(column.platforms.keys())

    #Define arrays
    lat_d = np.full(n_radars, np.nan)
    lat_m = np.full(n_radars, np.nan)
    lat_s = np.full(n_radars, np.nan)
    lon_d = np.full(n_radars, np.nan)
    lon_m = np.full(n_radars, np.nan)
    lon_s = np.full(n_radars, np.nan)
    elevation = np.full(n_radars, np.nan)
    name  = np.full(n_radars, '',dtype='<U10')
    type_  = np.full(n_radars, '',dtype='<U10')
    op_mode = np.full(n_radars, '',dtype='<U10')
    wavelength = np.full(n_radars, '',dtype='<U10')
    freq =np.full(n_radars, '',dtype='<U10')
    beam_width = np.full(n_radars, np.nan)
    gate_size = np.full(n_radars, np.nan)
    timestamp = np.full(n_radars, '',dtype='<U20')
    offset    = np.full(n_radars, '',dtype='<U10')
    
    #loop through each radar and fill in the arrays
    for i, radar in enumerate(column.platforms.keys()):
        lat_d[i] = column.platforms[radar]['lat_d']
        lat_m[i] = column.platforms[radar]['lat_m']
        lat_s[i] = column.platforms[radar]['lat_s']
        lon_d[i] = column.platforms[radar]['lon_d']
        lon_m[i] = column.platforms[radar]['lon_m']
        lon_s[i] = column.platforms[radar]['lon_s']
        elevation[i] = column.platforms[radar]['elev']
        name[i] = column.platforms[radar]['plat_name']
        type_[i] = column.platforms[radar]['plat_type']
        op_mode[i] = column.platforms[radar]['operation_mode']
        wavelength[i] = column.platforms[radar]['wavelength']
        freq[i] = column.platforms[radar]['frequency']
        beam_width[i] = column.platforms[radar]['beam_width']
        gate_size[i] = column.platforms[radar]['gate_size']
        timestamp[i] = column.platforms[radar]['timestamp']
        offset[i] = column.platforms[radar]['offset_vs_main']
        
    #remove once filled array for each radar
    for key, val in list(column.platforms.items()):
        #pdb.set_trace()
        del column.platforms[key]
        
    #put arrays into dictionary
    radar_info = {'plat_name':name,
                  'lat_d':lat_d, 
                  'lat_m':lat_m, 
                  'lat_s':lat_s,
                  'lon_d':lon_d, 
                  'lon_m':lon_m, 
                  'lon_s':lon_s,
                  'elev':elevation,
                  'plat_type':type_,
                  'operation_mode':op_mode,
                  'wavelength':wavelength,
                  'frequency':freq,
                  'beam_width':beam_width,
                  'gate_size':gate_size,
                  'timestamp':timestamp,
                  'offset_vs_main':offset}
                  
    #add radar_info for multiple radars to object
    column.add_platform_to_object(radar_info, plat_name='lev2_radars')
    
    return column

def set_plat_values(plat_name, data_file):
        '''
        ;------------------------------------------------------------
        ; get info about a platform's data, for use in new .nc attributes
        ;------------------------------------------------------------
        ;   result = set_plat_values(plat_name, data_file, main_plat_info)
        ;   plat_name:  string ID of current platform
        ;   data_file:  platform's original data file
        ;   main_plat_info:  info structure for the main_plat, as returned 
        ;           by set_main_plat_values.pro
        ;           will use timestamp here to get time offset
        ;           for current plat
        ;
        ;   result:  structure with this plat's info to include in the
        ;       new .nc file's attributes for records purposes:
        ;       (tags are arrays for instances of > 1 unit in col box)
        ;
        ;   Ground-based radars:    NPOL, 88Ds, D3R, DOW6, MRRs
        ;       .lat_d:     int radar's lat degrees
        ;       .lat_m:     int radar's lat minutes
        ;       .lat_s:     int radar's lat seconds
        ;       .lon_d:     int radar's lon degrees
        ;       .lon_m:     int radar's lon minutes
        ;       .lon_s:     int radar's lon seconds
        ;       .elev :     float, radar's elev (altitude) ASL [m]
        ;       .plat_name:     string ID of radar's name
        ;       .plat_type:     string label for type of platform: 'radar' or 'MRR'
        ;       .operation_mode:    string ID for scan mode: 'PPI', 'RHI', 'UNK'
        ;       .wavelength:        float radar wavelength [m] -- array for D3R [Ku, Ka]
        ;       .frequency:     float radar frequency [GHz] -- array for D3R [Ku, Ka]
        ;       .beam_width:        float radar beam width [deg] -- array for D3R [Ku, Ka]
        ;       .gate_size:     float radar gate size [m] -- array for D3R [Ku, Ka]
        ;   FOR NPOL, 88Ds, D3R, DOW6:
        ;       .timestamp: 15 char 'YYYYMMDD_HHMMSS' timestamp (vol scan start time)
        ;       .offset_vs_main:    int/long # of seconds time offset vs main_plat
        ;               POSITIVE OFFSET: plat  LATER  than main_plat
        ;               NEGATIVE OFFSET: plat EARLIER than main_plat
        ;   FOR MRRs:
        ;       .timestamp: 15 char 'YYYYMMDD_HHMMSS' timestamp for center of time
        ;                   interval contained in the column file
        ;       .offset_vs_main:      int/long # of seconds time offset vs main_plat
        ;       .time_interval_width: int # of mins before & after main_plat time included
        ;
        ;   Ground-based disdrometers:  APUs, 2DVDs
        ;       .lat_d:     int disdo's lat degrees
        ;       .lat_m:     int disdo's lat minutes
        ;       .lat_s:     int disdo's lat seconds
        ;       .lon_d:     int disdo's lon degrees
        ;       .lon_m:     int disdo's lon minutes
        ;       .lon_s:     int disdo's lon seconds
        ;       .plat_name:     string ID of unit's name
        ;       .plat_type:     string label for type of disdrometer: 'APU' or '2DVD'
        ;       .operation_mode:    string ID for op mode: 'APU' or '2DVD'
        ;       .timestamp:     15 char 'YYYYMMDD_HHMMSS' timestamp for center of time
        ;               interval contained in the column file
        ;       .offset_vs_main:    int/long # of seconds time offset vs main_plat
        ;       .time_interval_width:int # of mins before & after main_plat time included
        ;
        ;   Ground-based rain gauges:   Pluvios, gauges (tipping buckets)
        ;       .lat_d:     int gauge's lat degrees
        ;       .lat_m:     int gauge's lat minutes
        ;       .lat_s:     int gauge's lat seconds
        ;       .lon_d:     int gauge's lon degrees
        ;       .lon_m:     int gauge's lon minutes
        ;       .lon_s:     int gauge's lon seconds
        ;       .plat_name:     string ID of unit's name
        ;       .plat_type:     string label for type of gauge: 'pluvio' or 'gauge'
        ;       .operation_mode:    string ID for op mode: 'pluvio' or 'gmin_gauge'
        ;       .timestamp:     15 char 'YYYYMMDD_HHMMSS' timestamp for center of time
        ;               interval contained in the column file
        ;       .offset_vs_main:    int/long # of seconds time offset vs main_plat
        ;       .time_interval_width:int # of mins before & after main_plat time included
        ;
        ;   Satellite-based sensors:    GMI, DPR, 2BCMB
        ;       .plat_name:     string ID of sensor's name
        ;       .plat_type:     string label for this sensor product: 'GMI' or 'DPR'
        ;       .frequency:     string stating operating frequencies for the instrument [GHz]
        ;       .file_type:     string ID of PPS file type, eg: '2A-CS-CONUS'
        ;       .algorithm:     string algorithm name as from PPS filename (immediately
        ;               after data type in PPS name), eg: 'GPROF2014v2-0'
        ;       .orbit_num:     string GPM orbit number (after end time in PPS name)
        ;       .data_version:  string version nunber of the dataset (just before file
        ;               extension in PPS name)
        ;       .timestamp:     15 char 'YYYYMMDD_HHMMSS' timestamp for pixel that
        ;               includes/is closest to main_plat's lat/lon location
        ;       .offset_vs_main:    int/long # of seconds time offset vs main_plat
        ;       .timestamp_cntr:        15 char timestamp for pixel including/closest to
        ;               column box center pt (instead of main_plat)
        ;
        ;   Dependencies:
        ;       rsl_d3r_to_radar.pro    add-on to RSL_in_IDL for NASA D3R radar data
        ;       rsl_cfradial_to_radar.pro  add-on to RSL_in_IDL for Cfradial data
        ;       get_dpm.pro     utility to return number of days in the month
        ;       jd2time.pro     IDL Coyote utility to convert from Julian day
        ;
        ;------------------------------------------------------------
        '''
        #define radar identifiers
        radars = ['NPOL','D3R','DOW6','KABR','KCAE','KEAX','KICT','KLIX',
                  'KMRX','KTLH','KAKQ','KCCX','KEVX','KILN','KLOT','KMVX','KTLX',
                  'KAMX','KCLX','KFSD','KILX','KLSX','KNQA','KTWX','KAPX','KCRP', 
                  'KFTG','KINX','KLTX','KOKX','KTYX','KARX','KDDC','KFWS','KIWX',
                  'KLZK','KPAH','PAEC','KBMX','KDGX','KGRK','KJAX','KMHX','KRAX',
                  'PAIH','KBOX','KDLH','KGRR','KJGX','KMKX','KSGF','PGUA','KBRO',
                  'KDMX','KGSP','KJKL','KMLB','KSHV','PHKI','KBUF','KDOX','KHGX',
                  'KLCH','KMOB','KSRX','PHMO','KBYX','KDVN','KHTX','KLGX','KMQT',
                  'KTBW','TJUA']
        #define apus
        apus = ['apu01','apu02','apu03','apu04','apu05',
                     'apu06','apu07','apu08','apu09','apu10',
                     'apu11','apu12','apu13','apu14','apu15',
                     'apu16','apu17','apu18','apu20',
                     'apu21','apu23','apu25','apu27','apu28','apu30']
        
        #define 2DVDs, APU_Pluvios, MRRs
        twodvds = ['SN25','SN35','SN36','SN37','SN38','SN70']
        apu_pluvio = ['apu04_pluvio','apu10_pluvio','apu30_pluvio']
        mrr = ['MRR2-01','MRR2-02','MRR2-03','MRR2-04']
        
        #pdb.set_trace()
        if plat_name in radars:
            plat_type = 'radar'
            
            #unzip file (if needed)
            #file_basename = os.path.basename(data_file)
            #if file_basename.endswith('.gz'):
            #    cf_file = sim.ungzip_file(data_file)
            #else:
            #    cf_file = data_file
            radar = pyart.io.read(data_file, file_field_names=True)
        elif plat_name in apus:
            plat_type='APU'
        elif plat_name in twodvds:
            plat_type='2DVD'
        elif plat_name in apu_pluvio:
            plat_type='pluvio'
        elif plat_name in mrr:
            plat_type='MRR'
        elif plat_name == 'GMI':
            plat_type='GMI'
        elif plat_name == 'DPR':
            plat_type='DPR'
        elif plat_name == '2BCMB':
            plat_type='2BCMB'
        else:
            err = f'Platform {plat_name} is not defined in set_plat_values method -- please check'
            raise NameError(err)
            
        #For radar platforms:  Get info using pyart metadata:
        if plat_type == 'radar':
            r_lat_decdeg = radar.latitude['data'][0]
            r_lon_decdeg = radar.longitude['data'][0]
            plat_lat_d, plat_lat_m, plat_lat_s, plat_lon_d, plat_lon_m, plat_lon_s = dd2dms(r_lat_decdeg, r_lon_decdeg)
            plat_elev  = radar.altitude['data'][0] #[m]
        
            plat_freq = radar.instrument_parameters['frequency']['data'][0] #[Hz]
            plat_wave  = 2.998E8 / plat_freq #[m] (c/freq)
            plat_beam = radar.instrument_parameters['radar_beam_width_h']['data'][0] #[deg]
            plat_gate = radar.range['meters_between_gates'] #[m]
        
            if plat_name == 'D3R':
                #Must set a few items differently for the D3R:
                plat_wave = [0.0216,0.0084] #[Ku-, Ka-]  in [m], (approx. [3e8 m/s]/[13.91e9 Hz], [3e8 m/s]/[35.56e9 Hz] )
                plat_freq = [13.91,35.56]   #[Ku-, Ka-]  in [GHz]
                plat_beam = [0.86,0.9]  #[Ku-, Ka-]  in [deg]
            
                print('D3R related...need to test/modify code')
                pdb.set_trace()
                #ncID = NCDF_OPEN(data_file)
                #NCDF_ATTGET, ncID, 'Altitude', plat_elev, /global
                #NCDF_CLOSE, ncID
        
            plat_file_type = radar.scan_type.upper()
            if plat_file_type != 'PPI' and plat_file_type != 'RHI':
                print('Have to test this section...')
                pdb.set_trace()
                if plat_file_type == 'sector':  #mainly for DOW6
                    plat_file_type = 'PPI'  #mainly for DOW6
                else:
                    plat_file_type = 'UNK'
                
            # Get info from radar pyart object
            radar_datetime = pyart.util.datetime_from_radar(radar)
            #plat_year = radar_datetime.year
            #plat_mon  = radar_datetime.month
            #plat_day  = radar_datetime.day
            #plat_hr   = radar_datetime.hour
            #plat_min  = radar_datetime.minute
            #plat_sec  = radar_datetime.second
            
            #string_year = str(plat_year).zfill(4)
            #string_mon  = str(plat_mon).zfill(2)
            #string_day  = str(plat_day).zfill(2)
            #string_hr   = str(plat_hr).zfill(2)
            #string_min  = str(plat_min).zfill(2)
            #string_sec  = str(plat_sec).zfill(2)
            #plat_timestamp = string_year+string_mon+string_day+'_'+string_hr+string_min+string_sec
            plat_timestamp = radar_datetime.strftime("%Y%m%d_%H%M%S")
        
            #get time offset vs main_plat's timestamp
            pdb.set_trace() ## need main_plat_datetime
            offset = (radar_datetime - main_plat_datetime).total_seconds()
            #main_plat_timestamp = self.main_plat_info['timestamp'] #'YYYYMMDD_HHMMSS'
            #main_year = int(main_year)
            #main_mon  = int(main_mon)
            #main_day  = int(main_day)
            #main_hr   = int(main_hr)
            #main_min  = int(main_min)
            #main_sec  = int(main_sec)
            #off_sec = plat_sec - main_sec
            #off_min = plat_min - main_min
            #off_hr  = plat_hr  - main_hr	#these bigger intervals should be same,
            #off_day = plat_day - main_day	#but including here for completeness
            #off_mon = plat_mon - main_mon
            #off_year= plat_year - main_year
            #off_min = off_min*60        # sec in the min
            #off_hr  = off_hr*3600       # sec in the hrs
            #off_day = off_day*86400     # sec in the days (should be 0)
            #off_mon = off_mon*2.628E6   # sec in the mons (should be 0)
            #off_year= off_year*3.154E7  # sec in the yrs (should be 0)
            #offset = off_sec + off_min + off_hr + off_day + off_mon + off_year
            #pdb.set_trace()
            
            #Set up plat_info dictionary to return:
            #(this is the same set up as for MRR, below)
            plat_info = {'lat_d':plat_lat_d, 'lat_m':plat_lat_m, 'lat_s':plat_lat_s,
                        'lon_d':plat_lon_d, 'lon_m':plat_lon_m, 'lon_s':plat_lon_s,
                        'elev':plat_elev,
                        'plat_name':plat_name,
                        'plat_type':plat_type,
                        'operation_mode':plat_file_type,
                        'wavelength':f'{plat_wave:.5f}',
                        'frequency':f'{plat_freq/1E9:.5f}',
                        'beam_width':plat_beam,
                        'gate_size':plat_gate,
                        'timestamp':plat_timestamp,
                        'offset_vs_main':offset}
            #end plat_type = 'radar' block
            #--------------------------------------
            
        #--------------------------------------
        # For APU platforms:  set file_type to 'APU'
        #  set location based on info in
        #     platform_location.pro <<-- MUST BE UP TO DATE!
        #  use timestamp of main_plat - first, verify it is availble
        if plat_type == 'APU':
            apu_location = get_platform_loc(plat_name)
            plat_lat_d = apu_location[0]
            plat_lat_m = apu_location[1]
            plat_lat_s = apu_location[2]
            plat_lon_d = apu_location[3]
            plat_lon_m = apu_location[4]
            plat_lon_s = apu_location[5]
        
            # NEW - WILL SET TIME & OFFSET INFO IN THE APU MODULE:
            plat_timestamp = 'set_in_get_apu'
            offset = -9999
            plat_file_type = 'APU'
        
            # Set up plat_info structure to return:
            # (same set up as for 2DVD, below)
            plat_info = {'lat_d':plat_lat_d, 'lat_m':plat_lat_m, 'lat_s':plat_lat_s,
                         'lon_d':plat_lon_d, 'lon_m':plat_lon_m, 'lon_s':plat_lon_s,
                         'plat_name':plat_name,
                         'plat_type':plat_type,
                         'operation_mode':plat_file_type,
                         'timestamp':plat_timestamp,
                         'offset_vs_main':offset}  

        #end plat_type = 'APU' block
        #--------------------------------------
        
        #--------------------------------------
        # For 2DVD platforms:  set file_type to '2DVD'
        #  set location based on info in
        #    platform_location.pro <<-- MUST BE UP TO DATE!
        #  timestamp & offset info set in 2DVD module, so 
        #  can just send a 'dummy.txt' file in here since 
        #  don't need to read in the data file
        if plat_type == '2DVD':
            unit_location = get_platform_loc(plat_name)
            plat_lat_d = unit_location[0]
            plat_lat_m = unit_location[1]
            plat_lat_s = unit_location[2]
            plat_lon_d = unit_location[3]
            plat_lon_m = unit_location[4]
            plat_lon_s = unit_location[5]
          
            # NEW - WILL SET TIME & OFFSET IN THE 2DVD MODULE:	
            plat_timestamp = 'set_in_get_2dvd'
            offset = -9999  
            plat_file_type = '2DVD'
        
            # Set up plat_info structure to return:
            # (same set up as for APU, above)
            plat_info = {'lat_d':plat_lat_d, 'lat_m':plat_lat_m, 'lat_s':plat_lat_s,
                         'lon_d':plat_lon_d, 'lon_m':plat_lon_m, 'lon_s':plat_lon_s,
                         'plat_name':plat_name,
                         'plat_type':plat_type,
                         'operation_mode':plat_file_type,
                         'timestamp':plat_timestamp,
                         'offset_vs_main':offset} 
        
        #end plat_type = '2DVD' block
        #--------------------------------------
        
        #--------------------------------------
        # For Pluvio platforms:  set file_type to 'pluvio'
        #  set location based on info in
        #     platform_location.pro <<-- MUST BE UP TO DATE!
        #  leave timestamp tag to be set in
        #	get_pluvio_for_column.pro
        #  so only havta read file in once
        
        #end plat_type = 'pluvio' block
        #--------------------------------------
        
        #--------------------------------------
        # For MRR platforms: set file_type to 'MRR'
        # set location based on info in 
        #   platform_location.pro <<-- MUST BE UP TO DATE!
        # leave the timestamp tag to be set in 
        #	grab_mrr_for_column.pro
        # bc can take a while to read in file
        if plat_type == 'MRR':
            mrr_location = get_platform_loc(plat_name)
            plat_lat_d = mrr_location[0]
            plat_lat_m = mrr_location[1]
            plat_lat_s = mrr_location[2]
            plat_lon_d = mrr_location[3]
            plat_lon_m = mrr_location[4]
            plat_lon_s = mrr_location[5]
            plat_elev  = mrr_location[6]
        
            # Reset these in get_mrr_for_column.pro after read in the file:
            plat_timestamp = 'set_in_get_mrr'
            offset = -9999
            plat_gate = -9999
        
            # Set these manually for MRR since not in the .ave file:
            plat_file_type = 'MRR'
            plat_wave = 0.012   #[m] (approx. [3e8 m/s]/[24.24e9 Hz] )
            plat_freq = 24.24   #[GHz]
            plat_beam = 1.5     #[deg]
        
            #Set up plat_info structure to return:
            # (same set up as for ground-based scanning radars, above)
            plat_info = {'lat_d':plat_lat_d, 'lat_m':plat_lat_m, 'lat_s':plat_lat_s,
                        'lon_d':plat_lon_d, 'lon_m':plat_lon_m, 'lon_s':plat_lon_s,
                        'elev':plat_elev,
                        'plat_name':plat_name,
                        'plat_type':plat_type,
                        'operation_mode':plat_file_type,
                        'wavelength':plat_wave,
                        'frequency':plat_freq,
                        'beam_width':plat_beam,
                        'gate_size':plat_gate,
                        'timestamp':plat_timestamp,
                        'offset_vs_main':offset}
        
        #end plat_type = 'MRR' block
        #--------------------------------------
        
        #--------------------------------------
        # For the GPM GMI platform: 
        # will not have lat/lon tags for this platform bc satellite moves
        # set the plat timestamp as the timestamp for pixel that includes/
        #   is closest to the main_plat lat/lon location during the OP,
        #   but leave this to be set in get_gprofgmi_for_column.pro and/or
        #   get_lev1cgmi_for_column.pro once read the data in
        
        #end plat_type = 'GMI' block
        #--------------------------------------
        
        #--------------------------------------
        # For the GPM DPR platform: 
        # will not have lat/lon tags for this platform bc satellite moves
        # set the plat timestamp as the timestamp for pixel that includes/
        #   is closest to the main_plat lat/lon location during the OP,
        #   but leave this to be set in get_dpr_for_column.pro once read
        #   the data in (like do for MRR, GMI)
        # file_type will be set based on the algo name in L2A-DPR/ filename
        if plat_type == 'DPR':
          
            #set type, algo, orbit, version as from PPS name convention:
            dpr_fname_fields = os.path.basename(data_file).split('.')
            plat_file_type   = dpr_fname_fields[0]  #PPS data type
            sat_algo_name    = dpr_fname_fields[3]  #algorithm name
            sat_orbit_number = dpr_fname_fields[5]  #orbit number
            sat_data_version = dpr_fname_fields[6]  #data version 
        
            #Reset these in get_dpr_for_column.pro after read in the file:
            plat_timestamp = 'set_in_get_dpr'
            offset = -9999
            plat_timestamp_cntr = 'set_in_get_dpr'
        
            #Set DPR channels frequencies manually:
            plat_freq = 'Ka: 35.5 GHz, Ku: 13.6 GHz'
        
            #Set up plat_info structure to return:
            #(same set up as for GMI, above)
            plat_info = {'plat_name':plat_name,
                         'plat_type':plat_type,
                         'frequency':plat_freq,
                         'file_type':plat_file_type,
                         'algorithm':sat_algo_name,
                         'orbit_num':sat_orbit_number,
                         'data_version':sat_data_version,
                         'timestamp':plat_timestamp,
                         'offset_vs_main':offset,
                         'timestamp_cntr':plat_timestamp_cntr}
        #end plat_type = 'DPR' block
        #--------------------------------------
        
        #--------------------------------------
        # For the GPM 2BCMB Combined DPR+GMI algorithm product: 
        # will not have lat/lon tags for this platform bc satellite moves
        # set the plat timestamp as the timestamp for pixel that includes/
        #   is closest to the main_plat lat/lon location during the OP,
        #   but leave this to be set in get_dpr_for_column.pro once read
        #   the data in (like do for MRR, GMI, DPR)
        
        #end plat_type = '2BCMB' block
        #--------------------------------------
        
        print(f'--- {plat_name} timestamp: {plat_timestamp}')
        print(f'--- {plat_name} file type: {plat_file_type}')
        #pdb.set_trace()
        return plat_info

# ***************************************************************************************

