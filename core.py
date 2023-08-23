import os, sys, glob
import numpy as np
from datetime import datetime, timedelta
import pdb #pdb.set_trace()
import pyart
import tools as sim

class GridColumn:
    def __init__(self, time, variables, platforms, box):
        self.time = time
        self.variables = variables
        self.platforms = platforms
        self.box = box
        
        
    def add_variable(self, variable_name, dic):
        if variable_name in self.variables:
            err = f'A variable with name: {variable_name} already exists'
            raise ValueError(err)
        if 'data' not in dic:
            raise KeyError("dic must contain a 'data' key")
        radar_dims = (self.box.nz, self.box.nx, self.box.ny)
        ground_dims = (self.box.nz, self.box.nx, self.box.ny, self.box.nt)
        if dic["data"].shape != radar_dims and dic["data"].shape != ground_dims:
            err = f" 'data' has invalid shape, should be {radar_dims} or {ground_dims}"
            raise ValueError(err) 
        self.variables[variable_name] = dic
        
    def add_variable_to_object(self, var, var_name='UNKNOWN', units='UNKNOWN', 
                                long_name='UNKNOWN', badval='UNKNOWN'):
        var_dict = {'long_name':long_name,
                    'units': units,
                    'badval': np.nan,
                    'data': var,
                    }
                    
        self.add_variable(var_name, var_dict)
        
    def add_platform_to_object(self, platform, plat_name='UNKNOWN'):
        self.platforms[plat_name] = platform
    
    # method copied from pyart Helmus and Collis, 2016
    def info(self, level='compact', out=sys.stdout):
        if level == 'c':
            level='compact'
        elif level == 's':
            level='standard'
        
        if level not in ['compact', 'standard']:
            raise ValueError('invalid level parameter')
        
        #print column initial parameters
        self.box.print_box_parameters()
        
        #print all variables
        print("variables:", file=out)
        for var_name, var_dic in self.variables.items():
            self._dic_info(var_name, level, out, var_dic, 1)

    #method copied from pyart Helmus and Collis, 2016
    def _dic_info(self, attr, level, out, dic=None, ident_level=0):
        if dic is None:
            dic = getattr(self, attr)
        
        ilvl0 = "\t" * ident_level
        ilvl1 = "\t" * (ident_level + 1)
        
        if dic is None:
            print(str(attr) + ": None", file=out)
            return
            
        # make a string summary of the data key if it exists.
        if "data" not in dic:
            d_str = "Missing"
        elif not isinstance(dic["data"], np.ndarray):
            d_str = "<not a ndarray>"
        else:
            data = dic["data"]
            t = (data.dtype, data.shape)
            d_str = "<ndarray of type: {} and shape: {}>".format(*t)

        # compact, only data summary
        if level == "compact":
            print(ilvl0 + str(attr) + ":", d_str, file=out)

        # standard, all keys, only summary for data
        elif level == "standard":
            print(ilvl0 + str(attr) + ":", file=out)
            print(ilvl1 + "data:", d_str, file=out)
            for key, val in dic.items():
                if key == "data":
                    continue
                print(ilvl1 + key + ":", val, file=out)
        return

class Box:
    def __init__(self,
                box_center,
                cntr_lat_deg,
                cntr_lon_deg,
                box_limit,
                box_spacing,
                n_horiz_grid_boxes,
                n_vert_boxes,
                vert_limit,
                vert_spacing,
                halftime_interval,
                lat_values,
                lon_values,
                have_xy_values,
                x_values,
                y_values,
                z_values,
                t_values,
                plat_dims,
                ground_plat_dims,
                nx,
                ny,
                nz,
                nt,
                box_max_lat,
                box_min_lat,
                box_max_lon,
                box_min_lon):
        
        self.box_center = box_center
        self.cntr_lat_deg = cntr_lat_deg
        self.cntr_lon_deg = cntr_lon_deg
        self.box_limit    = box_limit
        self.box_spacing  = box_spacing
        self.n_horiz_grid_boxes = n_horiz_grid_boxes
        self.n_vert_boxes = n_vert_boxes
        self.vert_limit = vert_limit
        self.vert_spacing = vert_spacing
        self.halftime_interval = halftime_interval
        self.lat_values = lat_values
        self.lon_values = lon_values
        self.have_xy_values = have_xy_values
        self.x_values = x_values
        self.y_values = y_values
        self.z_values = z_values
        self.t_values = t_values
        self.plat_dims = plat_dims
        self.ground_plat_dims = ground_plat_dims
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.nt = nt
        self.box_max_lat = box_max_lat
        self.box_min_lat = box_min_lat
        self.box_max_lon = box_max_lon
        self.box_min_lon = box_min_lon
        
    def print_box_parameters(self):
        #Print out to terminal the grid coords for the full box:
        print('---------------------------------------------')
        print(f'   Column Box Centered On: {self.box_center}')
        print(f'      lat: {self.cntr_lat_deg:0.2f} lon: {self.cntr_lon_deg:0.2f}')
        print('   Column Box Grid:')
        print(f'       horiz spacing [m]: {self.box_spacing}')
        print(f'       horiz extent [m]: {self.box_limit}')
        print(f'       horiz grid boxes: {self.n_horiz_grid_boxes}')
        print(f'       vert spacing [m]: {self.vert_spacing}')
        print(f'       vert extent [m]: {self.vert_limit}')
        print(f'       vert grid boxes: {self.n_vert_boxes}')
        print('   Column Box Grid Latitudes:')
        print('      ', ' '.join('{:0.2f}'.format(i) for i in self.lat_values))
        #print('        "{}".join( self.lat_values))}')
        print('   Column Box Grid Longitudes: ')
        print('     ', ' '.join('{:0.2f}'.format(i) for i in self.lon_values))
        #print(f'        {" ".join(map(str, self.lon_values))}')
        print('   Column Box Grid Vertical Levels [m]:')
        print(f'       {" ".join(map(str,self.z_values))}')
        #print('      ', ' '.join('{:0.2f}'.format(i) for i in self.z_values))
        #print(('         {:0.4f}'*len(self.z_values).format(*self.z_values)))

def define_box(center_on, box_spacing, box_limit, vert_spacing, vert_limit, halftime_interval):
    #center_on = kwargs['center_loc']
    #box_spacing = kwargs['column_grid_horiz_spacing']
    #box_limit = kwargs['column_grid_horiz_limit']
    #vert_spacing = kwargs['column_grid_vert_spacing']
    #vert_limit = kwargs['column_grid_vert_limit']

    #Set box center location:  degrees, minutes, seconds, + or - sign on degrees
    box_c_lat_d, box_c_lat_m, box_c_lat_s, \
    box_c_lon_d, box_c_lon_m, box_c_lon_s, \
    elv = sim.get_platform_loc(center_on)
    
    #Get center location in decimal degrees
    cntr_lat_deg, cntr_lon_deg = sim.dms2dd(box_c_lat_d, box_c_lat_m, box_c_lat_s, 
                                        box_c_lon_d, box_c_lon_m, box_c_lon_s)
                                        
    #Get # of small boxes in the big box per x/y dir:
    n_horiz_grid_boxes = box_limit / box_spacing
    n_vert_boxes = vert_limit / vert_spacing
    if n_horiz_grid_boxes % 1.0 != 0.0:
       sys.exit('---ERROR:  BOX SPACING MUST EVENLY SUBDIVIDE HORIZONTAL BOX LIMIT!---')
    if n_vert_boxes % 1.0 != 0.0:
       sys.exit('---ERROR:  VERTICAL SPACING MUST EVENLY SUBDIVIDE VERTICAL LIMIT!---')

    n_vert_boxes = int(n_vert_boxes)
    n_horiz_grid_boxes = int(n_horiz_grid_boxes)
    #print(f'  No. of horiz small boxes: {n_grid_boxes}') #testing
    #print(f'  No. of  vert small boxes: {n_vert_boxes}') #testing
    
    #Will have (# of small boxes)+1 values for defining lat/lons:
    #these will be doubles in decimal degrees
    lat_values = np.zeros(n_horiz_grid_boxes+1)
    lon_values = np.zeros(n_horiz_grid_boxes+1)
    
    # Make float arrays for x, y, z
    have_xy_vals = False
    x_values = np.zeros(n_horiz_grid_boxes+1)
    y_values = np.zeros(n_horiz_grid_boxes+1)
    
    #Determine where center point will be, then get lat/lons for each grid point:
    #    odd # of grid boxes: center point is middle point of middle grid box
    #   even # of grid boxed: center point is at 4 corners of middle boxes
    if n_horiz_grid_boxes % 2 == 0:
        #even # of boxes, odd # of defining points
        #need to arrange center point at middle corner
        #lat of center point = middle lat of grid
        #lon of center point = middle lon of grid
        mid_sub = int(n_horiz_grid_boxes / 2)
        x_values[mid_sub] = 0.0
        y_values[mid_sub] = 0.0
        
        lat_values[mid_sub] = cntr_lat_deg
        lon_values[mid_sub] = cntr_lon_deg

        #after center corner, now set lat/lons for remaining edges:
        n_edges_from_center = mid_sub
        for edge_from_center in range(1, n_edges_from_center+1):
            #from middle corner, have ((n_grid_boxes)/2) to each side of big box
            #    call my Function: xy2ll  -->  takes distance in [km]

            x_dist = (box_spacing/1000.0)*edge_from_center  #[km]
            y_dist = (box_spacing/1000.0)*edge_from_center  #[km]
  
            # use for lons to East of center corner:
            to_east_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, x_dist, 0.0)
  
            # use for lons to West of center corner:
            to_west_latlon= sim.xy2ll(cntr_lat_deg, cntr_lon_deg, -1.0*x_dist, 0.0)
  
            # use for lats to North of center corner:
            to_north_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, 0.0, y_dist)
  
            # use for lons to South of center corner:
            to_south_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, 0.0, -1.0*y_dist)
  
            lon_values[mid_sub+edge_from_center] = to_east_latlon[1]
            lon_values[mid_sub-edge_from_center] = to_west_latlon[1]
            lat_values[mid_sub+edge_from_center] = to_north_latlon[0]
            lat_values[mid_sub-edge_from_center] = to_south_latlon[0]
            
            #for x and y
            x_values[mid_sub+edge_from_center] = x_dist*1000.0 #meters
            x_values[mid_sub-edge_from_center] = -1.0*x_dist*1000.0 #meters
            y_values[mid_sub+edge_from_center] = y_dist*1000.0 #meters
            y_values[mid_sub-edge_from_center] = -1.0*y_dist*1000.0 #meters
        have_xy_vals = True

    if n_horiz_grid_boxes % 2 == 1:
        print('...this section needs to be tested...')
        pdb.set_trace() ##add the self.x_values, self.y_values below as appropriate
        #odd # of boxes, even # of defining points
        #need to arrange center point at center of middle box
        #mid_subs: will have 2: just before &  just after center point
        mid_subs = [ ((n_horiz_grid_boxes+1)/2)-1 , (n_horiz_grid_boxes+1)/2 ]

        #first grid box edges to east & north will be center point + 0.5*spacing
        #first grid box edges to west & south will be center point - 0.5*spacing
        first_x_dist = (box_spacing/1000.0)*0.5  #[km]
        first_y_dist = (box_spacing/1000.0)*0.5  #[km]
        first_east_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, first_x_dist, 0.0)
        first_west_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, -1.0*first_x_dist, 0.0)
        first_north_latlon= sim.xy2ll(cntr_lat_deg, cntr_lon_deg, 0.0, first_y_dist)
        first_south_latlon= sim.xy2ll(cntr_lat_deg, cntr_lon_deg, 0.0, -1.0*first_y_dist)
        lon_values[mid_subs[0]] = first_west_latlon[1]
        lon_values[mid_subs[1]] = first_east_latlon[1]
        lat_values[mid_subs[0]] = first_south_latlon[0]
        lat_values[mid_subs[1]] = first_north_latlon[0]

        #rest of points to west & south: center - 0.5*spacing - (# of grid boxes to current sub * spacing)   
        #For subscripts less than the mid_subs
        for edge_sub in range(0, mid_subs[0]-1+1):
            #get # of full grid boxes from current subscript to 1st mid_subs
            n_full_boxes = mid_subs[0] - edge_sub
  
            #distance will be the first half spacing, then full spacing for each grid box
            x_dist = first_x_dist + (n_full_boxes*box_spacing)/1000.0  #[km]
            y_dist = first_y_dist + (n_full_boxes*box_spacing)/1000.0  #[km]
            to_west_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, -1.0*x_dist, 0.0)
            to_south_latlon= sim.xy2ll(cntr_lat_deg, cntr_lon_deg, 0.0, -1.0*y_dist)
            lon_values[edge_sub] = to_west_latlon[1]
            lat_values[edge_sub] = to_south_latlon[0]

        #rest of points to east & north: center + 0.5*spacing + (# of grid boxes to current sub * spacing)
        #For subscripts beyond the mid_subs
        for edge_sub in range(mid_subs[1]+1, n_horiz_grid_boxes+1):
            #get # of full grid boxes from current subscript to 2nd mid_subs
            n_full_boxes = edge_sub - mid_subs[1]
  
            #distance will be the first half spacing, then full spacing for each grid box
            x_dist = first_x_dist + (n_full_boxes*box_spacing)/1000.0  #[km]
            y_dist = first_y_dist + (n_full_boxes*box_spacing)/1000.0  #[km]
            to_east_latlon = sim.xy2ll(cntr_lat_deg, cntr_lon_deg, x_dist, 0.0)
            to_north_latlon= sim.xy2ll(cntr_lat_deg, cntr_lon_deg, 0.0, y_dist)
            lon_values[edge_sub] = to_east_latlon[1]
            lat_values[edge_sub] = to_north_latlon[0]
    have_xy_vals = True
            
    #Now that have the lat/lons for each grid point in the box, set the vertical coords:
    z_values = np.arange(0,n_vert_boxes+1)*float(vert_spacing)

    #define the dimensions for satellite and ground radar
    nx = lon_values.shape[0]
    ny = lat_values.shape[0]
    nz = z_values.shape[0]
    nt = 2 * halftime_interval + 1
    plat_dims = (nz, nx, ny)
    ground_plat_dims = (nz, nx, ny, nt)
    
    #define the time lag values
    t_values = np.arange(-halftime_interval, halftime_interval+1)
    
    #define the bounds of the column box
    box_max_lat = lat_values.max()
    box_min_lat = lat_values.min()
    box_max_lon = lon_values.max()
    box_min_lon = lon_values.min()

    return Box(center_on,cntr_lat_deg, cntr_lon_deg, box_limit, box_spacing, n_horiz_grid_boxes,
               n_vert_boxes, vert_limit, vert_spacing, halftime_interval, lat_values, lon_values, have_xy_vals,
               x_values, y_values, z_values, t_values, plat_dims, ground_plat_dims, nx, ny, nz, nt, 
               box_max_lat, box_min_lat, box_max_lon, box_min_lon)
               
def search_data(year, mon, day, output_dir, data_dir, main_plat, VN_radars):

    #Set paths for SIMBA-generated files:
    #-------------------------------------

    #Set dir where new column .nc file is to be placed:
    out_dir = f'{output_dir}/column_nc_file/{year}/{mon}{day}/'
    os.makedirs(out_dir, exist_ok=True)

    #Set dirs for holding full grid (pyart-gridding) output for ground-based
    #radars before SIMBA sub-sets to column grid:
    #full_grid_dir = f'{output_dir}/full_grid_output/'
    #os.makedirs(full_grid_dir, exist_ok=True)

    #Search paths for platforms to be included in SIMBA columns
    #-------------------------------------

    #current platforms as of May 3, 2023
    #[NPOL, 88Ds (KDOX,), Gauges, APU, 2DVD, MRR, GPM-DPR]

    ###########         main_plat_form     ####################
    print(f'...Searching for Main Platform {main_plat} file...')
    main_plat_dir = f'{data_dir}/{main_plat}/{year}/{mon}{day}/{main_plat}*'
    #print(main_plat_dir)
    files = glob.glob(main_plat_dir)
    if len(files) < 1:
       #print(f'*** NO {main_plat} FILES FOUND ***')
       sys.exit(f'*** NO {main_plat} FILE FOUND! {main_plat} (as Main Platform) FILE MUST BE DEFINED...EXITING PROGRAM ***')
    elif len(files) > 1:
       print(f'*** MULTIPLE {main_plat} FILES FOUND... USING {os.path.basename(files[0])} ***')
       #print(files)
       #print(f'...using {files[0]} ...')
       main_plat_file = files[0]
       print('   --> file found')
    else:
       main_plat_file = files[0]
       print('   --> file found')
   
    #unzip file (if needed)
    file_basename = os.path.basename(main_plat_file)
    if file_basename.endswith('.gz'):
       main_plat_file = sim.ungzip_file(main_plat_file)

    ###########         VN Radar(s)      ###################
    print(f'...Searching for VN Radars {", ".join(map(str,VN_radars))} file...')

    #first check number of radar IDs
    VN_IDS = np.array([])
    VN_radar_files = np.array([])
    nradars = len(VN_radars)
    if nradars == 0:
        print('*** NO VN RADARS PROVIDED...SKIPPING ***')
    else:
        for radar_id in VN_radars:
            #pdb.set_trace()
            VN_radar_dir = f'{data_dir}/{radar_id}/{year}/{mon}{day}/{radar_id}*'
            files = glob.glob(VN_radar_dir)
            if len(files) < 1:
                print(f'*** NO {radar_id} FILES FOUND! SKIPPING... {radar_id} DATA WILL NOT BE INCLUDED ***')
                #print(f'*** SKIPPING... {self.nexrad_88D_ID} DATA WILL NOT BE INCLUDED ***')
                #nexrad_file=''
            elif len(files) > 1:
                print(f'*** MULTIPLE {radar_id} FILES FOUND... USING {os.path.basename(files[0])} ***')
                #print(files)
                #print(f'...using {files[0]} ...')
                #VN_radar_files = np.append(VN_radar_files, files[0]
                #VN_IDS = np.append(VN_IDS, radar_id)
                this_file = files[0]
                print('   --> file found')
            else:
                this_file = files[0]
                VN_IDS = np.append(VN_IDS, radar_id)
                print('   --> file found')
            if len(files) >= 1:
                #unzip file (if needed)
                file_basename = os.path.basename(this_file)
                if file_basename.endswith('.gz'):
                    VN_radar_files = np.append(VN_radar_files, sim.ungzip_file(this_file))
                else:
                    VN_radar_files = np.append(VN_radar_files, this_file)
       
    ###########         Gauges         ####################
    print('...Searching for Gauge files...')
    gauges_dir = f'{data_dir}/Gauge/{year}/*.gmin'
    files = sorted(glob.glob(gauges_dir))
    if len(files) < 1:
        print(f'*** NO GAUGE .GMIN FILES FOUND! SKIPPING... GAUGE DATA WILL NOT BE INCLUDED ***')
        #print(f'*** SKIPPING... GAUGE DATA WILL NOT BE INCLUDED ***')
        gauge_gmin_files = ''
    else:
        gauge_gmin_files = files
        print('   --> files found')
   
    ###########         APU         ####################
    print('...Searching for APU files...')
    apu_dir = f'{data_dir}/APU/{year}/{mon}{day}/*'
    files = glob.glob(apu_dir)
    if len(files) < 1:
        print(f'*** NO APU FILES FOUND! SKIPPING... APU DATA WILL NOT BE INCLUDED ***')
        #print(f'*** SKIPPING... APU DATA WILL NOT BE INCLUDED ***')
        apu_files = ''
    else:
        apu_files = apu_dir[0:-2] #get_apu_for_column method will search for specific files
        print('   --> files found')
       
    ###########         2DVD         ####################
    print('...Searching for 2DVD files...')
    twodvd_dir = f'{data_dir}/2DVD/{year}/{mon}{day}/2dvd*'
    files = glob.glob(twodvd_dir)
    if len(files) < 1:
        print(f'*** NO 2DVD FILES FOUND! SKIPPING... 2DVD DATA WILL NOT BE INCLUDED ***')
        #print(f'*** SKIPPING... 2DVD DATA WILL NOT BE INCLUDED ***')
        twodvd_files = ''
    else:
        twodvd_files = files
        print('   --> files found')
       
    ###########         MRR         ####################
    print('...Searching for MRR file...')
    mrr_dir = f'{data_dir}/MRR/{year}/{mon}{day}/*.ave'
    files = glob.glob(mrr_dir)
    if len(files) < 1:
        print(f'*** NO MRR FILES FOUND! SKIPPING... MRR DATA WILL NOT BE INCLUDED ***')
        #print(f'*** SKIPPING... MRR DATA WILL NOT BE INCLUDED ***')
        mrr_files= ''
    else:
        mrr_files = files
        print('   --> file found')
       
    ###########         GPM-DPR         ####################
    print('...Searching for GPM-DPR file...')
    dpr_dir = f'{data_dir}/GPM/{year}/{mon}{day}/2A-CS-CONUS.GPM.DPR.V8-20180723.{year}{mon}{day}*.HDF5'
    files = glob.glob(dpr_dir)
    #pdb.set_trace()
    if len(files) < 1:
        print(f'*** NO DPR FILE FOUND! SKIPPING... DPR DATA WILL NOT BE INCLUDED ***')
        #print(f'*** SKIPPING... DPR DATA WILL NOT BE INCLUDED ***')
        dpr_file = ''
    else:
        dpr_file = files[0]
        print('   --> file found')
       
    ###########         NEW PLATFORM HERE         ####################
    print('...Searching for MRMS file(s)...')
    mrms_dir = f'{data_dir}/MRMS/{year}/{mon}{day}/*'
    files = glob.glob(mrms_dir)
    if len(files) < 1:
        print(f'*** NO MRMS FILES FOUND! SKIPPING... MRMS DATA WILL NOT BE INCLUDED ***')
        #print(f'*** SKIPPING... MRMS DATA WILL NOT BE INCLUDED ***')
        mrms_files = ''
    else:
        mrms_files = mrms_dir[0:-2] #get_mrms_for_column method will search for specific files
        print('   --> file(s) found')
    
    return main_plat_file, VN_IDS, VN_radar_files, gauge_gmin_files, apu_files, \
           twodvd_files, mrr_files, dpr_file, mrms_files, out_dir
           
def setup_main_plat(main_plat_name, main_plat_file):
    '''
    ;-------------------------------------------------------------------------------
    ;  method to perform the STEP 2 in SIMBA:  setting up the main_plat
    ;   items -- have placed this in a separate function to keep build_column.pro
    ;   neat as adding all VN 88Ds lengthens amount of code needed for this step
    ;-------------------------------------------------------------------------------
    ;  result = setup_main_plat(main_plat, [file_npol=file_npol, $
    ;			file_d3r_ku=file_d3r_ku, file_dow6=file_dow6, $
    ;			files_88d=files_88d])
    ;
    ;	main_plat:	string name of the main platform
    ;	file_*:		only need to send the file(s) that are for
    ;			the main platform.
    ;			-> If main_plat is the D3R, send entire Ku-band array
    ;			-> If main_plat is an 88d, send entire files_88d array
    ;
    ;	result:		main_plat_info dictionary
    ;
    ;  NOTE:  ea platform in this file must have correct LOCATION INFO in
    ;	  platform_locations.txt !!!
    ;
    ;-------------------------------------------------------------------------------
    '''
    #main_plat = kwargs['main_plat'].upper()
    
    radars = ['NPOL','D3R','DOW6','KABR','KCAE','KEAX','KICT','KLIX',
              'KMRX','KTLH','KAKQ','KCCX','KEVX','KILN','KLOT','KMVX','KTLX',
              'KAMX','KCLX','KFSD','KILX','KLSX','KNQA','KTWX','KAPX','KCRP', 
              'KFTG','KINX','KLTX','KOKX','KTYX','KARX','KDDC','KFWS','KIWX',
              'KLZK','KPAH','PAEC','KBMX','KDGX','KGRK','KJAX','KMHX','KRAX',
              'PAIH','KBOX','KDLH','KGRR','KJGX','KMKX','KSGF','PGUA','KBRO',
              'KDMX','KGSP','KJKL','KMLB','KSHV','PHKI','KBUF','KDOX','KHGX',
              'KLCH','KMOB','KSRX','PHMO','KBYX','KDVN','KHTX','KLGX','KMQT',
              'KTBW','TJUA']
    
    #Platform must be one of the VN radars
    if main_plat_name not in radars:
        sys.exit('------ MAIN PLATFORM ERROR: must be a valid option -----\n'
                f'-- VALID OPTIONS: {" ".join(self.radars)} --')
    else:#####most of the code below is from set_main_plat_values method -- June 28, 2023
        #check if file is zipped
        #file_basename = os.path.basename(self.npol_file)
        #if file_basename.endswith('.gz'):
        #    self.npol_file = sim.ungzip_file(self.npol_file)
        radar = pyart.io.read(main_plat_file, file_field_names=True)
            
        # Get info from radar pyart object
        radar_datetime = pyart.util.datetime_from_radar(radar) #this is same as above line -- remove and fix other code below
        main_plat_datetime = datetime(radar_datetime.year, radar_datetime.month,radar_datetime.day,radar_datetime.hour,
                                            radar_datetime.minute, radar_datetime.second)
        #main_year = str(radar_datetime.year).zfill(4)
        #main_mon  = str(radar_datetime.month).zfill(2)
        #main_day  = str(radar_datetime.day).zfill(2)
        #main_hr   = str(radar_datetime.hour).zfill(2)
        #main_min  = str(radar_datetime.minute).zfill(2)
        #main_sec  = str(radar_datetime.second).zfill(2)
        #main_sec  = round(main_sec) # dont think this is needed
        
        main_lat_deg = radar.latitude['data'][0]
        main_lon_deg = radar.longitude['data'][0]
        
        main_scan_type = radar.scan_type.upper()
        if main_scan_type != 'PPI' and main_scan_type != 'RHI': main_scan_type = 'UNK'
        
        #main_plat_timestamp = main_year+main_mon+main_day+'_'+main_hr+main_min+main_sec
        
        main_plat_params = {'datetime':main_plat_datetime, 'scan_type':main_scan_type,
                            'lat_deg':main_lat_deg, 'lon_deg':main_lon_deg}
        
    return main_plat_params
    
def get_time_lag(halftime_interval, main_plat_datetime):
        
    # Define interval_datetime array for the full data time interval:
    #  time interval = main_plat_timestamp +/- halftime_interval # of mins
    interval_datetime    = np.empty(2*halftime_interval + 1, dtype=object)
    n_interval_times    = len(interval_datetime)
    t_values = np.arange(-halftime_interval, halftime_interval+1)
    
    #need only Year, month, day, hour, minute in object
    main_datetime = datetime(main_plat_datetime.year, main_plat_datetime.month, main_plat_datetime.day, 
                             main_plat_datetime.hour, main_plat_datetime.minute)
            
    # use datetime to get timestamps for ea minute in the interval
    #main_datetime = datetime(int(self.main_year),int(self.main_mon),int(self.main_day),int(self.main_hr),int(self.main_min))
    for i, step in enumerate(np.arange(-halftime_interval,halftime_interval+1)):
        #compute delta datetime from the main datetime
        #have to convert numpy.int64 to python int for timedelta input using .item()
        interval_datetime[i] = main_datetime + timedelta(minutes=step.item())
        
    #define the dimensions for ground platform data arrays (apu, twodvd, gauges, mrr, etc)
    #  4-d:  [column z dir  X  column x dir  X  column y dir  X   times in interval ]
    #self.ground_plat_dims = (self.z_values.shape[0], self.lon_values.shape[0], self.lat_values.shape[0], self.n_interval_times)
    
    return interval_datetime