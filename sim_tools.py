import numpy as np
import math
from scipy import interpolate
from datetime import datetime
import os, shutil, gzip


def get_default_params_dict():
    '''
    #-------------------------------------------------------------------------------
    # function to provide default parameters to run simba columns
    #-------------------------------------------------------------------------------
    #   result = get_default_params_dict()
    #
    #   result: dictionary with keys and values to run simba columns
    #
    # Written by Charanjit S. Pabla
    #-------------------------------------------------------------------------------
    '''
    default_params = {
        'main_plat_name':       'NPOL', 
        'nexrad_88D_ID':        'KDOX',
        'center_on':            'WFFPad', 
        'box_spacing':          500,
        'box_limit':            5000, 
        'vert_spacing':         250,
        'vert_limit':           6000, 
        'halftime_interval':    5, 
        'data_dir':             '/home/cpabla/SIMBA_development/data_test',
        'output_dir':           '/home/cpabla/SIMBA_development',
        'dpr_version':          'V06A'
    }

    return default_params

 # ***************************************************************************************   

def check_kwargs(kwargs, default_kw):
    '''
    #-------------------------------------------------------------------------------
    # function to check user-provided kwargs against defaults, and if some defaults are
    # not provided by user, make sure they are provided to the function regardless
    #-------------------------------------------------------------------------------
    #   result = check_kwargs(kwargs, default_kw)
    #   kwargs:       dictionary with keys and values
    #   default_kw:   default dictionary with keys and values
    #
    #   result: will be dictionary
    #
    # Written by Jason L. Pippitt
    #-------------------------------------------------------------------------------
    '''
    for key in default_kw:
        if key not in kwargs:
            kwargs[key] = default_kw[key]
    return kwargs

# ***************************************************************************************

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

