"""
Python System for Integrating Multiplatform Data to Build the Atmospheric Column (pySIMBA)
pySIMBA v0.0

Author
-------
Charanjit S. Pabla, Code 612, NASA GSFC WFF/SSAI, Charanjit.S.Pabla@nasa.gov
NOTE: Original [IDL] code was written by Stephanie M. Wingo, ST11, NASA MSFC/UAH, Stephanie.M.Wingo@nasa.gov

Description
-----------
Creates a netCDF data product from various precipitation sensors integrated onto a 
user-defined geographic location.

Updates
-------
v0.0 - 2021 Imported into python from IDL v1.6.1
https://gpm-gv.gsfc.nasa.gov/SIMBA/


References
----------
Pabla, C. S., D. B. Wolff, D. A. Marks, S. M. Wingo, and J. L. Pippitt, 2022:
    GPM Ground Validation at NASA Wallops Precipitation Research Facility.
    J. Atmos. Oceanic Technol., 39, 1199-1215. https://doi.org/10.1175/JTECH-D-21-0122.1 

Wingo, S. M., W. A. Petersen, P. N. Gatlin, C. S. Pabla, D. A. Marks, and 
    D. B. Wolff, 2018: The System for Integrating Multiplatform Data to Build 
    the Atmospheric Column (SIMBA) precipitation observation fusion framework. 
    J. Atmos. Oceanic Technol., 35, 1353–1374. https://doi.org/10.1175/JTECH-D-17-0187.1

Dependencies
------------
pyart, numpy, pandas, xarray, scipy, netCDF4, h5py


"""
#from timeit import default_timer as timer
import warnings
import os, sys, pathlib, glob
import numpy as np
import pandas as pd
from netCDF4 import Dataset 
from datetime import datetime, timedelta
import pdb #pdb.set_trace()
import pyart
import xarray as xr
import argparse
import ast
import sim_tools as sim
import mrr
import gpm
#from pysimba_ import(sim_tools as sim, mrr, gpm)
warnings.filterwarnings("ignore")

# ***************************************************************************************
class BuildColumn(object):

    """
    Main class that will integrate all platforms to build the atmospheric column
    data product.
    """
    def __init__(self, year, mon, day, **kwargs):
        self.year = year
        self.mon = mon
        self.day = day
        
        # Check and fix missing user defined kwargs
        default_kw = sim.get_default_params_dict()
        kwargs = sim.check_kwargs(kwargs, default_kw)
        for key, value in kwargs.items():
            setattr(self, key, value)
        
    #def get_op_time(self):

    #    #get datetime
    #    utc=str(self.optime)
    #    hr=int(utc[0:2])
    #    mi=int(utc[2:4])
    #    sec=int(utc[4:])
    #    return datetime.datetime(self.year, self.mon, self.day, hr, mi, sec)
        
        
    def process_columns(self):
    
        #-----------
        #  step (1):  Set paths for searching platform data and output
        self.set_directories()
    
        #-----------
        #  step (2):  define grid points for the column
        #  based on center location, spacing, limit set at top
        self.define_box()
        
        #-----------
        #  step (3):  get the main_plat information and set timelag interval 
        #  NOTE: currently, main_plat must be a ground-based scanning radar!
        #  main_plat's timestamp will be the new .nc file's timestamp
        #  time lag will be set based on user provided halftime and main platform timestamp
        self.setup_main_plat()
        self.get_time_lag()
                                        
        #-----------
        # step (4):  if there is NPOL data available, 
        #  grid it to the column box framework and prepare fields to add to new .nc file
        if self.npol_file:
            print('-----in NPOL section-----')
            #get info from the data file header...will go into attributes
            self.npol_info = self.set_plat_values('NPOL', self.npol_file)
            self.grid_npol_for_column()
            #self.npol_data_in_column = False
        else:
            self.npol_data_in_column = False
            
       #-----------
       # *******THIS SECTION IS NOT READY YET*******
       # step (5):  if there is D3R data available,
       #  grid it to the column box framework and set fields to add to new .nc file
       #if files_D3R_Ka ne !NULL and files_D3R_Ku ne !NULL then begin
       #  print, '-----in D3R section-----'
       #INSERT D3R code here

        #-----------
        # step (6):  if there is 88D lev2 data available,
        #  grid it to the column box framework and prepare fields to add to new .nc file
        if self.nexrad_file:
            print('-----in 88D section-----')
            
            #--First, check how many 88D files to process
            #--Second, get plat_info for the 88D radar(s)
            #--Third, grid & get col_data for the 88D radar(s)
            
            #IF multiple radars: get lev2_plat_info for each 88D file then add to dictionary
            # do same for gridding col data
            #num_88Ds = len(self.nexrad)
            num_88Ds = 1
            
            #single 88D radar
            if num_88Ds == 1:
                #get plat_info for the 88D radar
                #self.siteID_88D = pathlib.Path(self.nexrad_file).name[0:4]
                #print(self.siteID_88D)
                self.lev2_info = self.set_plat_values(self.nexrad_88D_ID, self.nexrad_file)
                
                #grid & get col_data for the 88D radar
                self.grid_lev2_for_column()
            
            #for multiple 88D radars
            #elif num_88Ds > 1:
            #    print('....this section needs testing...')
            #    pdb.set_trace()
                #dictionaries to hold data for multiple radars
            #    lev2_plat_info = {}
            #    col_data_88d_all = {}
                
            #    for file in self.nexrad:
                    #get plat_info for the 88D radar
            #        siteID_88D = pathlib.Path(file).name[0:4]
            #        lev2_plat_info[siteID_88D] = self.set_plat_values(siteID_88D, file)

                    #grid & get col_data for the 88D radar
            #        col_data_88d_all[siteID_88D] = self.grid_radar_for_column(file)
            #print(lev2_coldata.keys())
            #pdb.set_trace()
        else:
            self.lev2_data_in_column = False
            
        #-----------
        # *******THIS SECTION IS NOT READY YET*******
        # step (7):  if there is DOW6 data read it in...
        #if file_dow6 ne !NULL then begin
        #  print, '-----in DOW6 section-----'
        # 
        #  ; check column grid location/range from D3R:
        #  column_grid_range_test = check_radar_range_to_column('DOW6', file_DOW6,$
        # 				column_box_params, main_plat_info)
        #  if column_grid_range_test eq -1 then goto, skip_dow6
        #
        #  ; get info for the radar:
        #  dow6_plat_info = set_plat_values('DOW6', file_dow6, main_plat_info)
        # 
        #  ; Auto-create .params file via make_params_file.pro, then call DOW6 module:
        #  dow6_params = make_params_file('DOW6', 'DOW6', dow6_plat_info.operation_mode, column_box_params, $
        #				column_grid_horiz_spacing, column_grid_vert_spacing, $
        #				params_dir)
        #  col_data_dow6 = grid_dow6_for_column(dow6_plat_info, file_dow6, dow6_params, $
        #				column_box_params, full_grid_dir)  
        #endif	;DOW6 file exists
        #skip_dow6:

        #-----------
        # step (8):  if there is APU data, read it in...
        if self.apu_files:
            print('-----in APU section-----')
            self.get_apu_for_column()
        else:
            self.apu_data_in_column = False
            
        #-----------
        # step (9):  if there is 2DVD data, read it in...
        if self.twodvd_files:
            print('-----in 2DVD section-----')
            self.get_2dvd_for_column()
        else:
            self.twodvd_data_in_column = False
            
        #-----------
        # *******THIS SECTION IS NOT READY YET*******
        # step (10):  if there is PLUVIO data read it in...
        #if pluvio_dir ne !NULL then begin
        #    print, '-----in PLUVIOS section-----'
        #    pluvio_info_and_data = get_pluvio_for_column(column_box_params, pluvio_dir, main_plat_info, halftime_interval)
        #    NOTE: have multiple structures in two_dvd_info_and_data:
        #     pluvio_info_and_data.******_plat_info

        #-----------
        # step (11):  if there is .gmin GAUGE data, read it in...
        if self.gauge_gmin_files:
            print('-----in GAUGES section-----')
            self.get_gauges_for_column()
        else:
            self.gauges_data_in_column = False
            
        #-----------
        # step (12):  if there is MRR data, read it in...
        if self.mrr_files:
            print('-----in MRR section-----')
            self.get_mrr_for_column()
        else:
            self.mrr_data_in_column = False
        
        #-----------
        # *******THIS SECTION IS NOT READY YET*******
        # step (13):  if there is GPROF-GMI data read it in...
        #if file_gprof_gmi ne !NULL then begin
        #  print, '-----in GPROF-GMI section-----'
        #  gmi_gprof_info_and_data = get_gprofgmi_for_column(column_box_params, file_gprof_gmi, $
        # 				main_plat_info, column_grid_horiz_limit, column_grid_horiz_spacing)
        #    ;NOTE: have multiple structures in gmi_gprof_info_and_data:
        #    ; gmi_gprof_info_and_data.gmi_plat_info
        #    ; gmi_gprof_info_and_data.col_data_gmi
        #endif ;GPROF-GMI file exists
        
        #-----------
        # *******THIS SECTION IS NOT READY YET*******
        # step (14):  if there is L1C-GMI data read it in...
        #if file_lev1c_gmi ne !NULL then begin
        #  print, '-----in L1C-GMI section-----'
        #  gmi_lev1c_info_and_data = get_lev1cgmi_for_column(column_box_params, file_lev1c_gmi, $
        # 			main_plat_info, column_grid_horiz_limit, column_grid_horiz_spacing)
        #  ;NOTE: have multiple structures in gmi_lev1c_info_and_data:
        #  ; gmi_lev1c_info_and_data.gmi_plat_info
        #  ; gmi_lev1c_info_and_data.col_data_gmi
        #endif  ;L1C GMI file exists
        
        #-----------
        # step (15):  if there is GPM 2ADPR data, read it in...
        if self.dpr_file:
            print('-----in DPR section-----')
            self.get_dpr_for_column()
        else:
            self.dpr_data_in_column = False
        
        #-----------
        # step (16):  create a new .nc file containing all data within the column box
        #   and include as much info as possible in the attributes. save the new file
        #    to the out_dir
        #  COLUMN netCDF FILES WILL BE NAMED AS:
        #       column_[main]_[center]_YYYYMMDD_HHMM.nc
        #           [main]:    string ID of main_plat used
        #           [center]:  string ID of box's center location
        #           timestamp: main_plat's time w/o seconds
        #         example:  column_NPOL_WFFpad_20160628_1449.nc
        self.version_string = 'v0.0'     #will be 'v0.1' for 1st release...
        self.write_column_nc_file()

    # ***************************************************************************************

    def set_directories(self):
    
        #Set paths for SIMBA-generated files:
        #-------------------------------------
        
        #Set dir where new column .nc file is to be placed:
        self.out_dir = f'{self.output_dir}/column_nc_file/{self.year:04}/{self.mon:02}{self.day:02}/'
        os.makedirs(self.out_dir, exist_ok=True)
        
        #Set dirs for holding full grid (pyart-gridding) output for ground-based
        #radars before SIMBA sub-sets to column grid:
        self.full_grid_dir = f'{self.output_dir}/full_grid_output/'
        os.makedirs(self.full_grid_dir, exist_ok=True)
    
        #Search paths for platforms to be included in SIMBA columns
        #-------------------------------------
        
        #current platforms as of May 3, 2023
        #[NPOL, 88Ds (KDOX,), Gauges, APU, 2DVD, MRR, GPM-DPR]
        
        ###########         NPOL         ####################
        print('...Searching for NPOL file...')
        npol_dir = f'{self.data_dir}/NPOL/{self.year}/{self.mon}{self.day}/NPOL*'
        files = glob.glob(npol_dir)
        if len(files) < 1:
            print(f'*** NO NPOL FILES FOUND ***')
            sys.exit('*** NPOL FILE MUST BE DEFINED...EXITING PROGRAM ***')
        elif len(files) > 1:
            print(f'*** MULTIPLE NPOL FILES FOUND ***')
            #print(files)
            print(f'...using {files[0]} ...')
            self.npol_file = files[0]
            print('   --> file found')
        else:
            self.npol_file = files[0]
            print('   --> file found')
            
        #unzip file (if needed)
        file_basename = os.path.basename(self.npol_file)
        if file_basename.endswith('.gz'):
            self.npol_file = sim.ungzip_file(self.npol_file)
        
        ###########         88D         ####################
        print(f'...Searching for WSR88D {self.nexrad_88D_ID} file...')
        nexrad_dir = f'{self.data_dir}/88D/{self.nexrad_88D_ID}/{self.year}/{self.mon}{self.day}/{self.nexrad_88D_ID}*'
        files = glob.glob(nexrad_dir)
        if len(files) < 1:
            print(f'*** NO {self.nexrad_88D_ID} FILES FOUND ***')
            print(f'*** SKIPPING... {self.nexrad_88D_ID} DATA WILL NOT BE INCLUDED ***')
            self.nexrad_file=''
        elif len(files) > 1:
            print(f'*** MULTIPLE {self.nexrad_88D_ID} FILES FOUND ***')
            #print(files)
            print(f'...using {files[0]} ...')
            self.nexrad_file = files[0]
            print('   --> file found')
        else:
            self.nexrad_file = files[0]
            print('   --> file found')
            
        #unzip file (if needed)
        if len(files) >= 1:
            file_basename = os.path.basename(self.nexrad_file)
            if file_basename.endswith('.gz'):
                self.nexrad_file = sim.ungzip_file(self.nexrad_file)
            
        ###########         Gauges         ####################
        print('...Searching for Gauge files...')
        gauges_dir = f'{self.data_dir}/Gauge/{self.year}/*.gmin'
        files = sorted(glob.glob(gauges_dir))
        if len(files) < 1:
            print(f'*** NO GAUGE .GMIN FILES FOUND ***')
            print(f'*** SKIPPING... GAUGE DATA WILL NOT BE INCLUDED ***')
            self.gauge_gmin_files = ''
        else:
            self.gauge_gmin_files = files
            print('   --> file found')
        
        ###########         APU         ####################
        print('...Searching for APU files...')
        apu_dir = f'{self.data_dir}/APU/{self.year}/{self.mon}{self.day}/*'
        files = glob.glob(apu_dir)
        if len(files) < 1:
            print(f'*** NO APU FILES FOUND ***')
            print(f'*** SKIPPING... APU DATA WILL NOT BE INCLUDED ***')
            self.apu_files = ''
        else:
            self.apu_files = apu_dir[0:-2] #get_apu_for_column method will search for specific files
            print('   --> file found')
            
        ###########         2DVD         ####################
        print('...Searching for 2DVD files...')
        twodvd_dir = f'{self.data_dir}/2DVD/{self.year}/{self.mon}{self.day}/2dvd*'
        files = glob.glob(twodvd_dir)
        if len(files) < 1:
            print(f'*** NO 2DVD FILES FOUND ***')
            print(f'*** SKIPPING... 2DVD DATA WILL NOT BE INCLUDED ***')
            self.twodvd_files = ''
        else:
            self.twodvd_files = files
            print('   --> file found')
            
        ###########         MRR         ####################
        print('...Searching for MRR file...')
        mrr_dir = f'{self.data_dir}/MRR/{self.year}/{self.mon}{self.day}/*.ave'
        files = glob.glob(mrr_dir)
        if len(files) < 1:
            print(f'*** NO MRR FILES FOUND ***')
            print(f'*** SKIPPING... MRR DATA WILL NOT BE INCLUDED ***')
            self.mrr_files= ''
        else:
            self.mrr_files = files
            print('   --> file found')
            
        ###########         GPM-DPR         ####################
        print('...Searching for GPM-DPR file...')
        dpr_dir = f'{self.data_dir}/GPM/{self.year}/{self.mon}{self.day}/2A-CS-CONUS.GPM.DPR.V8-20180723.{self.year}{self.mon}{self.day}*.{self.dpr_version}.HDF5'
        files = glob.glob(dpr_dir)
        if len(files) < 1:
            print(f'*** NO DPR FILE FOUND ***')
            print(f'*** SKIPPING... DPR DATA WILL NOT BE INCLUDED ***')
            self.dpr_file = ''
        else:
            self.dpr_file = files[0]
            print('   --> file found')
            
        ###########         NEW PLATFORM HERE         ####################

    # ***************************************************************************************
    def get_dpr_for_column(self):
        '''
        ;------------------------------------------------------------
        ; get GPM DPR parameters from input file, locate data
        ; within the column grid box, and return that subset
        ;------------------------------------------------------------
        ;
        ;   result = get_dpr_for_column(column_box_params, file_dpr, main_plat_info,$
        ;			column_grid_horiz_limit, column_grid_horiz_spacing)
        ;
        ;	column_box_params: column box grid parameters - the structure
        ;			   returned by define_box.pro
        ;	file_dpr:   	   Level 2 DPR data file (.HDF5 format)
        ;			    GPM data is available on the NASA PPS:
        ;			    https://pps.gsfc.nasa.gov/ppsindex.html
        ;	main_plat_info:	   main_plat_info structure returned by
        ;			    set_main_plat_values.pro
        ;	column_grid_horiz_limit: total horizontal extent of column
        ;				grid in [m]
        ;	column_grid_horiz_spacing:  column grid spacing in horizontal
        ;				directions, in [m]
        ;
        ;	result:		structure containing key 2ADPR paramaters at each
        ;			location in the colun box grid; substructures & tags:
        ;	dpr_info_and_data:	
        ;	  .dpr_plat_info:  sturcture of DPR's plat_info values
        ;			   mainly set up on set_plat_values.pro, but
        ;			   timestamp & offset tags are set in here
        ;
        ;	  .col_data_dpr:   structure of DPR data for the column grid,
        ;			   and name/units tags describing each field
        ;
        ;	   (also see more/full details in READMEs)
        ;
        ;	Dependencies:			in addition to other SIMBA .pros:
        ;	  read_2adpr_hdf5.pro		from Bob Morris' utility library for 
        ;					reading in GPM satellite products
        ;	  rsl_posn_to_azm_range.pro	utility to compute azimuth & distance
        ;					in [m] between input lat/lon points
        ;	  bad2nan.pro			utility to convert a flag value to NaNs
        ;
        ;------------------------------------------------------------
        '''
        
        # Set initial dpr_plat_info values:
        self.dpr_info = self.set_plat_values('DPR', self.dpr_file)

        #Set range to search w/in full DPR swath:
        limit_dpr_km = 10.0 * self.box_limit/1000.  #[km]
        print(f'searching in DPR [km]: {limit_dpr_km}')
    
        # HARD SET THESE BASED ON 2ADPR FILE SPECIFICATION .doc:
        HS_nBins = 88       #
        MS_nBins = 176      #   these values are provided in the
        NS_nBins = 176      #   
        HS_binSize = 250.   #[m]    PPS GPM File Specification Document,
        MS_binSize = 125.   #[m]    
        NS_binSize = 125.   #[m]    Section 5.33, p.1092
        badVal   = -9999
    
        # read in the 2ADPR file:
        dpr_data    = gpm.read_2adpr_hdf5(self.dpr_file)
        fileheader  = dpr_data['FileHeader']
        dpr_HS_data = dpr_data['HS']   #"high senstivitiy scan" swath
        dpr_MS_data = dpr_data['MS']   #"matched beam scan" swath
        dpr_NS_data = dpr_data['NS']   #"normal/nominal beam scan" swath
         
        #Keys contained in [dpr_data[fileheader]] sub-dictionary:  all either a string or an int:
        #  .DOI			.DOIauthority		.DOIshortname		.algorithmID
        #  .algorithmVersion	.filename		.satelliteName		.insturmentName
        #  .generationDateTime	.startGranuleDateTime	.stopGranuleDateTime	.granuleNumber
        #  .numberOfSwaths	.numberOfGrids		.granuleStart		.timeInterval
        #  .processingSystem	.productVersion		.emptyGranule		.missingData

        #Have 3 swaths, so get header values from each:
        #high sensitivity scan swath:

        HS_swath_header     = dpr_HS_data['swathHeader']
        HS_nScans_in_set    = HS_swath_header['NumberScansInSet']
        HS_max_nScans       = HS_swath_header['MaximumNumberScansTotal']
        HS_nScans_b4_gran   = HS_swath_header['NumberScansBeforeGranule']
        HS_nScans_granule   = HS_swath_header['NumberScansGranule']
        HS_nScans_af_gran   = HS_swath_header['NumberScansAfterGranule']
        HS_nPixels          = HS_swath_header['NumberPixels']
        HS_scan_type        = HS_swath_header['ScanType']
        
        #matched beam scan swath:
        MS_swath_header    = dpr_MS_data['swathHeader']
        MS_nScans_in_set   = MS_swath_header['NumberScansInSet']
        MS_max_nScans      = MS_swath_header['MaximumNumberScansTotal']
        MS_nScans_b4_gran  = MS_swath_header['NumberScansBeforeGranule']
        MS_nScans_granule  = MS_swath_header['NumberScansGranule']
        MS_nScans_af_gran  = MS_swath_header['NumberScansAfterGranule']
        MS_nPixels         = MS_swath_header['NumberPixels']
        MS_scan_type       = MS_swath_header['ScanType']

        #normal/nominal scan swath:
        NS_swath_header    = dpr_NS_data['swathHeader']
        NS_nScans_in_set   = NS_swath_header['NumberScansInSet']
        NS_max_nScans      = NS_swath_header['MaximumNumberScansTotal']
        NS_nScans_b4_gran  = NS_swath_header['NumberScansBeforeGranule']
        NS_nScans_granule  = NS_swath_header['NumberScansGranule']
        NS_nScans_af_gran  = NS_swath_header['NumberScansAfterGranule']
        NS_nPixels         = NS_swath_header['NumberPixels']
        NS_scan_type       = NS_swath_header['ScanType']

        #Now get the data fields for each swath:
        #high sensitivity scan swath: 
        #  note: verify is HS still in V07 data - will need to modify code
        #  note: also have pointers for groups: SRT, FLG, experimental
        #  note most 2 and 3D variables have been transposed to get appropriate dimensions
        HS_scanTime             = dpr_HS_data['scantime']   # [UNITS]
        HS_year                 = HS_scanTime['year'][:]
        HS_month                = HS_scanTime['month'][:]
        HS_day_of_month         = HS_scanTime['day'][:]
        #HS_day_of_year         = HS_scanTime['dayofYear'][:]
        HS_hour                 = HS_scanTime['hour'][:]
        HS_minute               = HS_scanTime['minute'][:]
        HS_second               = HS_scanTime['second'][:]
        #HS_millisec            = HS_scanTime['millisecond'][:]
        #HS_sec_of_day          = HS_scanTime['secondofDay'][:]
        HS_SCstatus             = dpr_HS_data['scanStatus']
        HS_dataQuality          = HS_SCstatus['dataQuality'][:]
        #HS_SCorientation       = HS_SCstatus['SCorientation'][:]
        #HS_FracGranNum         = HS_SCstatus['FractionalGranuelNumber'][:]
        HS_navigation           = dpr_HS_data['navigation']
        HS_SC_lat               = HS_navigation['scLat'][:]                     # +/-70 deg
        HS_SC_lon               = HS_navigation['scLon'][:]                     # +/-180 deg
        HS_SC_alt               = HS_navigation['scAlt'][:]                     # [m]
        HS_PREgroup             = dpr_HS_data['PRE']
        HS_sfc_type             = HS_PREgroup['landSurfaceType'][:]
        HS_flag_precip          = HS_PREgroup['flagPrecip'][:]
        HS_binRealSfc           = HS_PREgroup['binRealSurface'][:]
        HS_binClutFB            = HS_PREgroup['binClutterFreeBottom'][:]
        HS_binStormTop          = HS_PREgroup['binStormTop'][:]
        HS_ht_stormTop          = HS_PREgroup['heightStormTop'][:]              # [m]
        HS_zFac_meas            = HS_PREgroup['zFactorMeasured'][:]             # [dBZ]
        HS_locZenAng            = HS_PREgroup['localZenithAngle'][:]            # [deg]
        HS_elipsBO              = HS_PREgroup['ellipsoidBinOffset'][:]
        HS_VERgroup             = dpr_HS_data['VER']
        HS_binZeroDeg           = HS_VERgroup['binZeroDeg'][:]                  # 1-89 or 1-177
        HS_ht_ZeroDeg           = HS_VERgroup['heightZeroDeg'][:]               # [m]
        HS_PIA_noPrecip         = HS_VERgroup['piaNP'][:]                       # [dB]
        HS_atten_noPrecip       = HS_VERgroup['attenuationNP'][:]               # [dB/km]
        HS_CSFgroup             = dpr_HS_data['CSF']
        HS_flag_BB              = HS_CSFgroup['flagBB'][:]
        HS_ht_BB                = HS_CSFgroup['heightBB'][:]                    # [m]
        HS_width_BB             = HS_CSFgroup['widthBB'][:]                     # [m]
        HS_quality_BB           = HS_CSFgroup['qualityBB'][:]
        HS_typePrecip           = HS_CSFgroup['typePrecip'][:]
        HS_qualityTypePrecip    = HS_CSFgroup['qualityTypePrecip'][:]
        HS_DSDgroup             = dpr_HS_data['DSD']
        HS_DSD_phase            = HS_DSDgroup['phase'][:]
        HS_SLVgroup             = dpr_HS_data['SLV']
        HS_PIA_final            = HS_SLVgroup['piaFinal'][:]                    # [dB]
        HS_corZFac              = HS_SLVgroup['zFactorCorrected'][:]            # [dBZ]
        HS_corZFac_nearSfc      = HS_SLVgroup['zFactorCorrectedNearSurface'][:] # [dBZ]
        HS_precipRate           = HS_SLVgroup['precipRate'][:]                  # [mm/h]
        HS_precipRate_nearSfc   = HS_SLVgroup['precipRateNearSurface'][:]       # [mm/h]
        HS_precipWater          = HS_SLVgroup['precipWaterIntegrated'][:]       # [g m^-2]
        HS_precipRateAve24      = HS_SLVgroup['precipRateAve24'][:]              # [mm/h]
        HS_phase_nearSfc        = HS_SLVgroup['phaseNearSurface'][:]
        HS_binEchoBottom        = HS_SLVgroup['binEchoBottom'][:]
        HS_paramDSD             = HS_SLVgroup['paramDSD'][:]
        HS_SRTgroup             = dpr_HS_data['SRT']
        HS_eff_PIA              = HS_SRTgroup['pathAtten'][:]                   # [dB]
        HS_datasets             = dpr_HS_data['datasets']
        HS_latitude             = HS_datasets['latitude'][:]                    # [deg]
        HS_longitude            = HS_datasets['longitude'][:]                   # [deg]
        
        #matched beam scan swath:
        #  note: also have pointers for groups: SRT, FLG, experimental
        MS_scanTime             = dpr_MS_data['scantime']                       # [UNITS]
        MS_year                 = MS_scanTime['year'][:]
        MS_month                = MS_scanTime['month'][:]
        MS_day_of_month         = MS_scanTime['day'][:]
        #MS_day_of_year         = MS_scanTime['dayofYear'][:]
        MS_hour                 = MS_scanTime['hour'][:]
        MS_minute               = MS_scanTime['minute'][:]
        MS_second               = MS_scanTime['second'][:]
        #MS_millisec            = MS_scanTime['millisecond'][:]
        #MS_sec_of_day          = MS_scanTime['secondofDay'][:]
        MS_SCstatus             = dpr_MS_data['scanStatus']
        MS_dataQuality          = MS_SCstatus['dataQuality'][:]
        #MS_SCorientation       = MS_SCstatus['SCorientation'][:]
        #MS_FracGranNum         = MS_SCstatus['FractionalGranuelNumber'][:]
        MS_navigation           = dpr_MS_data['navigation']
        MS_SC_lat               = MS_navigation['scLat'][:]                     # +/-70 deg
        MS_SC_lon               = MS_navigation['scLon'][:]                     # +/-180 deg
        MS_SC_alt               = MS_navigation['scAlt'][:]                     # [m]
        MS_PREgroup             = dpr_MS_data['PRE']
        MS_sfc_type             = MS_PREgroup['landSurfaceType'][:]
        MS_flag_precip          = MS_PREgroup['flagPrecip'][:]
        MS_binRealSfc           = MS_PREgroup['binRealSurface'][:]
        MS_binClutFB            = MS_PREgroup['binClutterFreeBottom'][:]
        MS_binStormTop          = MS_PREgroup['binStormTop'][:]
        MS_ht_stormTop          = MS_PREgroup['heightStormTop'][:]              # [m]
        MS_zFac_meas            = MS_PREgroup['zFactorMeasured'][:]             # [dBZ]
        MS_locZenAng            = MS_PREgroup['localZenithAngle'][:]            # [deg]
        MS_elipsBO              = MS_PREgroup['ellipsoidBinOffset'][:]
        MS_VERgroup             = dpr_MS_data['VER']
        MS_binZeroDeg           = MS_VERgroup['binZeroDeg'][:]                  # 1-89 or 1-177
        MS_ht_ZeroDeg           = MS_VERgroup['heightZeroDeg'][:]               # [m]
        MS_PIA_noPrecip         = MS_VERgroup['piaNP'][:]                       # [dB]
        MS_atten_noPrecip       = MS_VERgroup['attenuationNP'][:]               # [dB/km]
        MS_CSFgroup             = dpr_MS_data['CSF']
        MS_flag_BB              = MS_CSFgroup['flagBB'][:]
        MS_ht_BB                = MS_CSFgroup['heightBB'][:]                    # [m]
        MS_width_BB             = MS_CSFgroup['widthBB'][:]                     # [m]
        MS_quality_BB           = MS_CSFgroup['qualityBB'][:]
        MS_typePrecip           = MS_CSFgroup['typePrecip'][:]
        MS_qualityTypePrecip    = MS_CSFgroup['qualityTypePrecip'][:]
        #MS_DSDgroup            = dpr_MS_data['DSD']
        #MS_DSD_phase           = MS_DSDgroup['phase'][:]   #this field DNE for MS
        MS_SLVgroup             = dpr_MS_data['SLV']
        MS_PIA_final            = MS_SLVgroup['piaFinal'][:]                    # [dB]
        MS_corZFac              = MS_SLVgroup['zFactorCorrected'][:]            # [dBZ]
        MS_corZFac_nearSfc      = MS_SLVgroup['zFactorCorrectedNearSurface'][:] # [dBZ]
        #MS_precipRate          = MS_SLVgroup['precipRate'][:] #this field DNE for MS
        MS_precipRate_nearSfc   = MS_SLVgroup['precipRateNearSurface'][:]       # [mm/h]
        MS_precipWater          = MS_SLVgroup['precipWaterIntegrated'][:]       # [g m^-2]
        MS_precipRateAve24      = MS_SLVgroup['precipRateAve24'][:]             # [mm/h]
        MS_phase_nearSfc        = MS_SLVgroup['phaseNearSurface'][:]
        MS_binEchoBottom        = MS_SLVgroup['binEchoBottom'][:]
        #MS_paramDSD            = MS_SLVgroup['paramDSD'][:]  #this field DNE for MS
        MS_SRTgroup             = dpr_MS_data['SRT']
        MS_eff_PIA              = MS_SRTgroup['pathAtten'][:]                   # [dB]
        MS_datasets             = dpr_MS_data['datasets']
        MS_latitude             = MS_datasets['latitude'][:]                    # [deg]
        MS_longitude            = MS_datasets['longitude'][:]                   # [deg]
        
        # normal/nominal scan swath:
        # these are in dimensions [SCANS x PIXELS x BINS] ; Example [734 x 49 x 176]
        NS_scanTime             = dpr_NS_data['scantime']
        NS_year                 = NS_scanTime['year'][:]
        NS_month                = NS_scanTime['month'][:]
        NS_day_of_month         = NS_scanTime['day'][:]
        #NS_day_of_year         = NS_scanTime['dayofYear'][:]
        NS_hour                 = NS_scanTime['hour'][:]
        NS_minute               = NS_scanTime['minute'][:]
        NS_second               = NS_scanTime['second'][:]
        #NS_millisec            = NS_scanTime['millisecond'][:]
        #NS_sec_of_day          = NS_scanTime['secondofDay'][:]
        NS_SCstatus             = dpr_NS_data['scanStatus']
        NS_dataQuality          = NS_SCstatus['dataQuality'][:]
        #NS_SCorientation       = NS_SCstatus['SCorientation'][:]
        #NS_FracGranNum         = NS_SCstatus['FractionalGranuelNumber'][:]
        NS_navigation           = dpr_NS_data['navigation']
        NS_SC_lat               = NS_navigation['scLat'][:]
        NS_SC_lon               = NS_navigation['scLon'][:]
        NS_SC_alt               = NS_navigation['scAlt'][:]
        NS_PREgroup             = dpr_NS_data['PRE']
        NS_sfc_type             = NS_PREgroup['landSurfaceType'][:]
        NS_flag_precip          = NS_PREgroup['flagPrecip'][:]
        NS_binRealSfc           = NS_PREgroup['binRealSurface'][:]
        NS_binClutFB            = NS_PREgroup['binClutterFreeBottom'][:]
        NS_binStormTop          = NS_PREgroup['binStormTop'][:]
        NS_ht_stormTop          = NS_PREgroup['heightStormTop'][:]              # [m]
        NS_zFac_meas            = NS_PREgroup['zFactorMeasured'][:]             # [dBZ] #####test
        NS_locZenAng            = NS_PREgroup['localZenithAngle'][:]            # [deg]
        NS_elipsBO              = NS_PREgroup['ellipsoidBinOffset'][:]
        NS_VERgroup             = dpr_NS_data['VER']
        NS_binZeroDeg           = NS_VERgroup['binZeroDeg'][:]                  # 1-89 or 1-177
        NS_ht_ZeroDeg           = NS_VERgroup['heightZeroDeg'][:]               # [m]
        NS_PIA_noPrecip         = NS_VERgroup['piaNP'][:]                       # [dB]
        NS_atten_noPrecip       = NS_VERgroup['attenuationNP'][:]               # [dB/km]
        NS_CSFgroup             = dpr_NS_data['CSF']
        NS_flag_BB              = NS_CSFgroup['flagBB'][:]
        NS_ht_BB                = NS_CSFgroup['heightBB'][:]                    # [m]
        NS_width_BB             = NS_CSFgroup['widthBB'][:]                     # [m]
        NS_quality_BB           = NS_CSFgroup['qualityBB'][:]
        NS_typePrecip           = NS_CSFgroup['typePrecip'][:]
        NS_qualityTypePrecip    = NS_CSFgroup['qualityTypePrecip'][:]
        NS_DSDgroup             = dpr_NS_data['DSD']
        NS_DSD_phase            = NS_DSDgroup['phase'][:]
        NS_SLVgroup             = dpr_NS_data['SLV']
        NS_PIA_final            = NS_SLVgroup['piaFinal'][:]                    # [dB]
        NS_corZFac              = NS_SLVgroup['zFactorCorrected'][:]            # [dBZ] ######test
        NS_corZFac_nearSfc      = NS_SLVgroup['zFactorCorrectedNearSurface'][:] # [dBZ]
        NS_precipRate           = NS_SLVgroup['precipRate'][:]                  # [mm/h]
        NS_precipRate_nearSfc   = NS_SLVgroup['precipRateNearSurface'][:]       # [mm/h]
        NS_precipWater          = NS_SLVgroup['precipWaterIntegrated'][:]       # [g m^-2]
        NS_precipRateAve24      = NS_SLVgroup['precipRateAve24'][:]             # [mm/h]
        NS_phase_nearSfc        = NS_SLVgroup['phaseNearSurface'][:]
        NS_binEchoBottom        = NS_SLVgroup['binEchoBottom'][:]
        NS_paramDSD             = NS_SLVgroup['paramDSD'][:]
        NS_SRTgroup             = dpr_NS_data['SRT']
        NS_eff_PIA              = NS_SRTgroup['pathAtten'][:]                   # [dB]
        NS_datasets             = dpr_NS_data['datasets']
        NS_latitude             = NS_datasets['latitude'] [:]                   # [deg]
        NS_longitude            = NS_datasets['longitude'][:]                   # [deg]
        #pdb.set_trace()
        # get dims values - must do each scan type separatly:
        HS_nScans = HS_nScans_granule       #also: HS_nPixels set above
        MS_nScans = MS_nScans_granule       #also: MS_nPixels set above
        NS_nScans = NS_nScans_granule       #also: NS_nPixels set above

        # DETERMINING HEIGHTS:
        # in the data:  HS/NS/MS: bin No 1   = top of column	= subscript 0
        #		HS:	  bin No 88  = earth ellipsoid  = subscript 87
        #		NS/MS:    bin No 176 = earth ellipsoid  = subscript 175
        #  get the hts of each bin:  See 2ADPR ATBD p.20-22 for explaination:
        #   height[current_bin] = ( (binOfEllipsoid - current_bin)*binSize + ellipsoidBinOffset) * cos(localZenithAngle)
        #	-> since ellipsoidBinOffset & localZenithAngle depend on where at in [pixels x scans] space, need those
        #	   dims, and since bins go up into vertical need that dim... -> HEIGHTS ARRAYS WILL NEED TO BE 3D!!!
        # NOTE:      dims of 3-D field arrays are as [BINS x PIXELS x SCANS], will havta reset dims below
        # ALSO NOTE: vert. subs of the _hts_ arrays BACKWARDS "normal" convention:  bigger sub = lower height
        #NEXT THREE LINES TAKE UP MOST TIME DUE TO triple nested for loop######
        HS_hts_m = gpm.get_hgt_of_bins(HS_nPixels, HS_nScans, HS_nBins, HS_binSize, HS_elipsBO, HS_locZenAng)
        MS_hts_m = gpm.get_hgt_of_bins(MS_nPixels, MS_nScans, MS_nBins, MS_binSize, MS_elipsBO, MS_locZenAng)
        NS_hts_m = gpm.get_hgt_of_bins(NS_nPixels, NS_nScans, NS_nBins, NS_binSize, NS_elipsBO, NS_locZenAng)
        #pdb.set_trace()
        # --- Dimensions are set up as:
        #       with '##' as 'HS', 'MS', 'NS' for each swath type
        # ##_[times]: 1d arrays of time components for ea scan:		[##_nScans]
        # ##_latitude:	2D arrays of latitudes for ea HS/MS/NS pixel:	[##_nScans x ##_nPixels]
        # ##_longitude:	2D arrays of longitudes for ea HS/MS/NS pixel:	[##_nScans x ##_nPixels]
        # ##_[most_fields]:  most fields are 2D arrays, vals @ ea pix:	[##_nScans x ##_nPixels]
        # ##_[3d_fields]:    the 3D fields have vals in ea BIN @ ea pix:[##_nScans x ##_nPixels x ##_nBins] - different from IDL

        # Identify 1st & Last scans within the limit_dpr_km range of the column grid center:
        # determine search range based on dpr search limit (10*col grid horiz extent)
        max_deg_lat = limit_dpr_km / 111.1
        max_deg_lon = limit_dpr_km / (np.cos(np.deg2rad(self.cntr_lat_deg))*111.1)
        HS_start_scan=0  ;  HS_end_scan=0  ;  HS_nscans2do=0  ; HS_start_found=0
        MS_start_scan=0  ;  MS_end_scan=0  ;  MS_nscans2do=0  ; MS_start_found=0
        NS_start_scan=0  ;  NS_end_scan=0  ;  NS_nscans2do=0  ; NS_start_found=0
        
        # for each HS/MS/NS swath, get the [scan, pixel] subset in the search area:
        # created function to do this in pySIMBA 
        HS_start_scan, HS_end_scan, HS_nscans2do = gpm.get_nscan_to_process(HS_nScans, 
        HS_nPixels, HS_longitude, HS_latitude, limit_dpr_km, self.cntr_lat_deg, self.cntr_lon_deg)
        HS_nscans2do=0 #Set to 0 manually due to V05B version update as a temp fix - Updated 10/01/18
        if HS_nscans2do == 0:
            # this should NOT happen...
            print(' --- GPM 2ADPR MODULE ERROR!!! --- ')
            print(' --- CAN NOT LOCATE DPR HS SCAN WITH COLUMN GRID CENTER POINT! --- ')
            print(' --- SKIPPING HS DPR DATA! --- NO HS DPR VALUES INTO COLUMN! --- ')
            #pdb.set_trace()
        
        MS_start_scan, MS_end_scan, MS_nscans2do = gpm.get_nscan_to_process(MS_nScans, 
        MS_nPixels, MS_longitude, MS_latitude, limit_dpr_km, self.cntr_lat_deg, self.cntr_lon_deg)
        if MS_nscans2do == 0:
            #this should NOT happen...
            print(' --- GPM 2ADPR MODULE ERROR!!! --- ')
            print(' --- CAN NOT LOCATE DPR MS SCAN WITH COLUMN GRID CENTER POINT! --- ')
            print(' --- SKIPPING MS DPR DATA! --- NO MS DPR VALUES INTO COLUMN! --- ')
            #pdb.set_trace()

        NS_start_scan, NS_end_scan, NS_nscans2do = gpm.get_nscan_to_process(NS_nScans, 
        NS_nPixels, NS_longitude, NS_latitude, limit_dpr_km, self.cntr_lat_deg, self.cntr_lon_deg)
        if NS_nscans2do == 0:
            #this should NOT happen...
            print(' --- GPM 2ADPR MODULE ERROR!!! --- ')
            print(' --- CAN NOT LOCATE DPR NS SCAN WITH COLUMN GRID CENTER POINT! --- ')
            print(' --- SKIPPING NS DPR DATA! --- NO NS DPR VALUES INTO COLUMN! --- ')
            #pdb.set_trace()

        print(f' -- HS: start, end, # of scans: {HS_start_scan}, {HS_end_scan}, {HS_nscans2do}')
        print(f' -- MS: start, end, # of scans: {MS_start_scan}, {MS_end_scan}, {MS_nscans2do}')
        print(f' -- NS: start, end, # of scans: {NS_start_scan}, {NS_end_scan}, {NS_nscans2do}')
        # --- if no scans from all 3 swaths include column grid center, return !NULL and exit module: ---
        # must be HS, MS, and NS for this, so that if NS swath is ok but HS & MS are not, module still goes
        #HS_nscans2do = 0  ;for testing only 
        #MS_nscans2do = 0  ;for testing only
        #NS_nscans2do = 0  ;for testing only
        #pdb.set_trace()
        if HS_nscans2do == 0 and MS_nscans2do == 0 and NS_nscans2do == 0:
            print(' --- APPEARS DPR COVERAGE DOES NOT INCLUDE COLUMN GRID ! --- ')
            print(' --- exiting 2ADPR module, 2ADPR data NOT set in column grid.')
            pdb.set_trace() #this section has not been tested yet
            return -1

        # Get (longitude,latitude) coords of each pixel in the HS, MS, NS ranges identified
        if HS_nscans2do > 0:
            HS_longitude_subset = HS_longitude[HS_start_scan:HS_end_scan+1,:]
            HS_latitude_subset  = HS_latitude[HS_start_scan:HS_end_scan+1,:]
        if MS_nscans2do > 0:
            MS_longitude_subset = MS_longitude[MS_start_scan:MS_end_scan+1,:]
            MS_latitude_subset  = MS_latitude[MS_start_scan:MS_end_scan+1,:]
        if NS_nscans2do >0:
            NS_longitude_subset = NS_longitude[NS_start_scan:NS_end_scan+1,:]
            NS_latitude_subset  = NS_latitude[NS_start_scan:NS_end_scan+1,:]

        # Test if the column-relative coords show the FOVs within the search area
        #  are actually within the column grid area:
        NS_lon_OK = 'Y' ; NS_lat_OK = 'Y' ; MS_lon_OK = 'Y' ; MS_lat_OK = 'Y' ; HS_lon_OK = 'Y' ; HS_lat_OK = 'Y'
        col_in_MSswath = 1 ; col_in_HSswath = 1
        if NS_longitude_subset.min() > min(self.lon_values): NS_lon_OK = 'N'
        if NS_latitude_subset.min()  > min(self.lat_values): NS_lat_OK = 'N'
        
        if MS_nscans2do > 0:
            if MS_longitude_subset.min() > min(self.lon_values): MS_lon_OK = 'N'
            if MS_latitude_subset.min()  > min(self.lat_values): MS_lat_OK = 'N'
        if HS_nscans2do > 0:
            if HS_longitude_subset.min() > min(self.lon_values): HS_lon_OK = 'N'
            if HS_latitude_subset.min()  > min(self.lat_values): HS_lat_OK = 'N'
        if NS_lon_OK == 'N' or NS_lat_OK == 'N':
            print(' --- APPEARS DPR for 2ADPR COVERAGE DOES NOT INCLUDE COLUMN GRID ! --- ')
            print(' --- exiting 2ADPR module, 2ADPR data NOT set in column grid.')
            pdb.set_trace() #this has not been tested
            return -1
        else:
            if MS_lon_OK == 'N' or MS_lat_OK == 'N':
                print( ' --- APPEARS DPR for 2ADPR MS SWATH DOES NOT INCLUDE COLUMN GRID ! --- ')
                print( ' --- 2ADPR module will ONLY set 2ADPR NS data in column grid - NO MS ! ')
                col_in_MSswath = 0
            if HS_lon_OK == 'N' or HS_lat_OK == 'N':
                print( ' --- APPEARS DPR for 2ADPR HS SWATH DOES NOT INCLUDE COLUMN GRID ! --- ')
                print( ' --- 2ADPR module will ONLY set 2ADPR NS data in column grid - NO HS ! ')
                col_in_HSswath = 0

        # Subset the DPR full swaths to only the search region:
        #  (takes out a smaller part from the full data arrays)
        HS_year                 = HS_year[HS_start_scan:HS_end_scan+1]
        HS_month                = HS_month[HS_start_scan:HS_end_scan+1]
        HS_day_of_month         = HS_day_of_month[HS_start_scan:HS_end_scan+1]
        #HS_day_of_year          = HS_day_of_year[HS_start_scan:HS_end_scan+1]
        HS_hour                 = HS_hour[HS_start_scan:HS_end_scan+1]
        HS_minute               = HS_minute[HS_start_scan:HS_end_scan+1]
        HS_second               = HS_second[HS_start_scan:HS_end_scan+1]
        #HS_millisec             = HS_millisec[HS_start_scan:HS_end_scan+1]
        #HS_sec_of_day           = HS_sec_of_day[HS_start_scan:HS_end_scan+1]
        HS_dataQuality          = HS_dataQuality[HS_start_scan:HS_end_scan+1]
        HS_sfc_type             = HS_sfc_type[HS_start_scan:HS_end_scan+1,:]
        HS_flag_precip          = HS_flag_precip[HS_start_scan:HS_end_scan+1,:]
        HS_binRealSfc           = HS_binRealSfc[HS_start_scan:HS_end_scan+1,:]
        HS_binClutFB            = HS_binClutFB[HS_start_scan:HS_end_scan+1,:]
        HS_binStormTop          = HS_binStormTop[HS_start_scan:HS_end_scan+1,:]
        HS_ht_stormTop          = HS_ht_stormTop[HS_start_scan:HS_end_scan+1,:]
        HS_locZenAng            = HS_locZenAng[HS_start_scan:HS_end_scan+1,:]
        HS_elipsBO              = HS_elipsBO[HS_start_scan:HS_end_scan+1,:]
        HS_binZeroDeg           = HS_binZeroDeg[HS_start_scan:HS_end_scan+1,:]
        HS_ht_ZeroDeg           = HS_ht_ZeroDeg[HS_start_scan:HS_end_scan+1,:]
        HS_flag_BB              = HS_flag_BB[HS_start_scan:HS_end_scan+1,:]
        HS_ht_BB                = HS_ht_BB[HS_start_scan:HS_end_scan+1,:]
        HS_width_BB             = HS_width_BB[HS_start_scan:HS_end_scan+1,:]
        HS_quality_BB           = HS_quality_BB[HS_start_scan:HS_end_scan+1,:]
        HS_typePrecip           = HS_typePrecip[HS_start_scan:HS_end_scan+1,:]
        HS_qualityTypePrecip    = HS_qualityTypePrecip[HS_start_scan:HS_end_scan+1,:]
        HS_PIA_final            = HS_PIA_final[HS_start_scan:HS_end_scan+1,:]
        HS_corZFac_nearSfc      = HS_corZFac_nearSfc[HS_start_scan:HS_end_scan+1,:]
        HS_precipRate_nearSfc   = HS_precipRate_nearSfc[HS_start_scan:HS_end_scan+1,:]
        HS_precipRateAve24      = HS_precipRateAve24[HS_start_scan:HS_end_scan+1,:]
        HS_phase_nearSfc        = HS_phase_nearSfc[HS_start_scan:HS_end_scan+1,:]
        HS_binEchoBottom        = HS_binEchoBottom[HS_start_scan:HS_end_scan+1,:]
        HS_eff_PIA              = HS_eff_PIA[HS_start_scan:HS_end_scan+1,:]
        HS_latitude             = HS_latitude[HS_start_scan:HS_end_scan+1,:]
        HS_longitude            = HS_longitude[HS_start_scan:HS_end_scan+1,:]
        HS_zFac_meas            = HS_zFac_meas[HS_start_scan:HS_end_scan+1,:,:]
        HS_PIA_noPrecip         = HS_PIA_noPrecip[HS_start_scan:HS_end_scan+1,:,:]
        HS_atten_noPrecip       = HS_atten_noPrecip[HS_start_scan:HS_end_scan+1,:,:]
        HS_DSD_phase            = HS_DSD_phase[HS_start_scan:HS_end_scan+1,:,:]
        HS_corZFac              = HS_corZFac[HS_start_scan:HS_end_scan+1,:,:]
        HS_precipRate           = HS_precipRate[HS_start_scan:HS_end_scan+1,:,:]
        HS_precipWater          = HS_precipWater[HS_start_scan:HS_end_scan+1,:,:]
        HS_paramDSD             = HS_paramDSD[HS_start_scan:HS_end_scan+1,:,:,:]
        
        MS_year                 = MS_year[MS_start_scan:MS_end_scan+1]
        MS_month                = MS_month[MS_start_scan:MS_end_scan+1]
        MS_day_of_month         = MS_day_of_month[MS_start_scan:MS_end_scan+1]
        #MS_day_of_year          = MS_day_of_year[MS_start_scan:MS_end_scan+1]
        MS_hour                 = MS_hour[MS_start_scan:MS_end_scan+1]
        MS_minute               = MS_minute[MS_start_scan:MS_end_scan+1]
        MS_second               = MS_second[MS_start_scan:MS_end_scan+1]
        #MS_millisec             = MS_millisec[MS_start_scan:MS_end_scan+1]
        #MS_sec_of_day           = MS_sec_of_day[MS_start_scan:MS_end_scan+1]
        MS_dataQuality          = MS_dataQuality[MS_start_scan:MS_end_scan+1]
        MS_sfc_type             = MS_sfc_type[MS_start_scan:MS_end_scan+1,:]
        MS_flag_precip          = MS_flag_precip[MS_start_scan:MS_end_scan+1,:]
        MS_binRealSfc           = MS_binRealSfc[MS_start_scan:MS_end_scan+1,:]
        MS_binClutFB            = MS_binClutFB[MS_start_scan:MS_end_scan+1,:]
        MS_binStormTop          = MS_binStormTop[MS_start_scan:MS_end_scan+1,:]
        MS_ht_stormTop          = MS_ht_stormTop[MS_start_scan:MS_end_scan+1,:]
        MS_locZenAng            = MS_locZenAng[MS_start_scan:MS_end_scan+1,:]
        MS_elipsBO              = MS_elipsBO[MS_start_scan:MS_end_scan+1,:]
        MS_binZeroDeg           = MS_binZeroDeg[MS_start_scan:MS_end_scan+1,:]
        MS_ht_ZeroDeg           = MS_ht_ZeroDeg[MS_start_scan:MS_end_scan+1,:]
        MS_flag_BB              = MS_flag_BB[MS_start_scan:MS_end_scan+1,:]
        MS_ht_BB                = MS_ht_BB[MS_start_scan:MS_end_scan+1,:]
        MS_width_BB             = MS_width_BB[MS_start_scan:MS_end_scan+1,:]
        MS_quality_BB           = MS_quality_BB[MS_start_scan:MS_end_scan+1,:]
        MS_typePrecip           = MS_typePrecip[MS_start_scan:MS_end_scan+1,:]
        MS_qualityTypePrecip    = MS_qualityTypePrecip[MS_start_scan:MS_end_scan+1,:]
        MS_PIA_final            = MS_PIA_final[MS_start_scan:MS_end_scan+1,:]
        MS_corZFac_nearSfc      = MS_corZFac_nearSfc[MS_start_scan:MS_end_scan+1,:]
        MS_precipRate_nearSfc   = MS_precipRate_nearSfc[MS_start_scan:MS_end_scan+1,:]
        MS_precipRateAve24      = MS_precipRateAve24[MS_start_scan:MS_end_scan+1,:]
        MS_phase_nearSfc        = MS_phase_nearSfc[MS_start_scan:MS_end_scan+1,:]
        MS_binEchoBottom        = MS_binEchoBottom[MS_start_scan:MS_end_scan+1,:]
        MS_eff_PIA              = MS_eff_PIA[MS_start_scan:MS_end_scan+1,:]
        MS_latitude             = MS_latitude[MS_start_scan:MS_end_scan+1,:]
        MS_longitude            = MS_longitude[MS_start_scan:MS_end_scan+1,:]
        MS_zFac_meas            = MS_zFac_meas[MS_start_scan:MS_end_scan+1,:,:]
        MS_PIA_noPrecip         = MS_PIA_noPrecip[MS_start_scan:MS_end_scan+1,:,:]
        MS_atten_noPrecip       = MS_atten_noPrecip[MS_start_scan:MS_end_scan+1,:,:]
        MS_corZFac              = MS_corZFac[MS_start_scan:MS_end_scan+1,:,:]
        MS_precipWater          = MS_precipWater[MS_start_scan:MS_end_scan+1,:,:]
        
        NS_year                 = NS_year[NS_start_scan:NS_end_scan+1]
        NS_month                = NS_month[NS_start_scan:NS_end_scan+1]
        NS_day_of_month         = NS_day_of_month[NS_start_scan:NS_end_scan+1]
        #NS_day_of_year          = NS_day_of_year[NS_start_scan:NS_end_scan+1]
        NS_hour                 = NS_hour[NS_start_scan:NS_end_scan+1]
        NS_minute               = NS_minute[NS_start_scan:NS_end_scan+1]
        NS_second               = NS_second[NS_start_scan:NS_end_scan+1]
        #NS_millisec             = NS_millisec[NS_start_scan:NS_end_scan+1]
        #NS_sec_of_day           = NS_sec_of_day[NS_start_scan:NS_end_scan+1]
        NS_dataQuality          = NS_dataQuality[NS_start_scan:NS_end_scan+1]
        NS_sfc_type             = NS_sfc_type[NS_start_scan:NS_end_scan+1,:]
        NS_flag_precip          = NS_flag_precip[NS_start_scan:NS_end_scan+1,:]
        NS_binRealSfc           = NS_binRealSfc[NS_start_scan:NS_end_scan+1,:]
        NS_binClutFB            = NS_binClutFB[NS_start_scan:NS_end_scan+1,:]
        NS_binStormTop          = NS_binStormTop[NS_start_scan:NS_end_scan+1,:]
        NS_ht_stormTop          = NS_ht_stormTop[NS_start_scan:NS_end_scan+1,:]
        NS_locZenAng            = NS_locZenAng[NS_start_scan:NS_end_scan+1,:]
        NS_elipsBO              = NS_elipsBO[NS_start_scan:NS_end_scan+1,:]
        NS_binZeroDeg           = NS_binZeroDeg[NS_start_scan:NS_end_scan+1,:]
        NS_ht_ZeroDeg           = NS_ht_ZeroDeg[NS_start_scan:NS_end_scan+1,:]
        NS_flag_BB              = NS_flag_BB[NS_start_scan:NS_end_scan+1,:]
        NS_ht_BB                = NS_ht_BB[NS_start_scan:NS_end_scan+1,:]
        NS_width_BB             = NS_width_BB[NS_start_scan:NS_end_scan+1,:]
        NS_quality_BB           = NS_quality_BB[NS_start_scan:NS_end_scan+1,:]
        NS_typePrecip           = NS_typePrecip[NS_start_scan:NS_end_scan+1,:]
        NS_qualityTypePrecip    = NS_qualityTypePrecip[NS_start_scan:NS_end_scan+1,:]
        NS_PIA_final            = NS_PIA_final[NS_start_scan:NS_end_scan+1,:]
        NS_corZFac_nearSfc      = NS_corZFac_nearSfc[NS_start_scan:NS_end_scan+1,:]
        NS_precipRate_nearSfc   = NS_precipRate_nearSfc[NS_start_scan:NS_end_scan+1,:]
        NS_precipRateAve24      = NS_precipRateAve24[NS_start_scan:NS_end_scan+1,:]
        NS_phase_nearSfc        = NS_phase_nearSfc[NS_start_scan:NS_end_scan+1,:]
        NS_binEchoBottom        = NS_binEchoBottom[NS_start_scan:NS_end_scan+1,:]
        NS_eff_PIA              = NS_eff_PIA[NS_start_scan:NS_end_scan+1,:]
        NS_latitude             = NS_latitude[NS_start_scan:NS_end_scan+1,:]
        NS_longitude            = NS_longitude[NS_start_scan:NS_end_scan+1,:]
        NS_zFac_meas            = NS_zFac_meas[NS_start_scan:NS_end_scan+1,:,:]
        NS_PIA_noPrecip         = NS_PIA_noPrecip[NS_start_scan:NS_end_scan+1,:,:]
        NS_atten_noPrecip       = NS_atten_noPrecip[NS_start_scan:NS_end_scan+1,:,:]
        NS_DSD_phase            = NS_DSD_phase[NS_start_scan:NS_end_scan+1,:,:]
        NS_corZFac              = NS_corZFac[NS_start_scan:NS_end_scan+1,:,:]
        NS_precipRate           = NS_precipRate[NS_start_scan:NS_end_scan+1,:,:]
        NS_precipWater          = NS_precipWater[NS_start_scan:NS_end_scan+1,:,:]
        NS_paramDSD             = NS_paramDSD[NS_start_scan:NS_end_scan+1,:,:,:]
        
        # pull apart combo data fields:
        HS_PIA_cloudwater = HS_PIA_noPrecip[:,:,0]
        HS_PIA_cloudice   = HS_PIA_noPrecip[:,:,1]
        HS_PIA_watervapor = HS_PIA_noPrecip[:,:,2]
        HS_PIA_oxygen     = HS_PIA_noPrecip[:,:,3]
        HS_TPW_liquid     = HS_precipWater[:,:,0]
        HS_TPW_ice        = HS_precipWater[:,:,1]
        HS_DSD_dBNw       = HS_paramDSD[:,:,:,0]
        HS_DSD_Dm         = HS_paramDSD[:,:,:,1]
        
        NS_PIA_cloudwater = NS_PIA_noPrecip[:,:,0]
        NS_PIA_cloudice   = NS_PIA_noPrecip[:,:,1]
        NS_PIA_watervapor = NS_PIA_noPrecip[:,:,2]
        NS_PIA_oxygen     = NS_PIA_noPrecip[:,:,3]
        NS_TPW_liquid     = NS_precipWater[:,:,0]
        NS_TPW_ice        = NS_precipWater[:,:,1]
        NS_DSD_dBNw       = NS_paramDSD[:,:,:,0]
        NS_DSD_Dm         = NS_paramDSD[:,:,:,1]
         
        MS_PIA_cloudwater = MS_PIA_noPrecip[:,:,0]
        MS_PIA_cloudice   = MS_PIA_noPrecip[:,:,1]
        MS_PIA_watervapor = MS_PIA_noPrecip[:,:,2]
        MS_PIA_oxygen     = MS_PIA_noPrecip[:,:,3]
        MS_TPW_liquid     = MS_precipWater[:,:,0]
        MS_TPW_ice        = MS_precipWater[:,:,1]

        # resize arrays of computed gate heights:  these arrays dims as [nScans x nPixels x nBins]
        HS_hts_m = HS_hts_m[HS_start_scan:HS_end_scan+1,:,:]
        NS_hts_m = NS_hts_m[NS_start_scan:NS_end_scan+1,:,:]
        MS_hts_m = MS_hts_m[MS_start_scan:MS_end_scan+1,:,:]
        #pdb.set_trace()
        # --- Dimensions are set up as:		EXCEPT: ##_nScans is really ##_nscans2do
        #       with '##' as 'HS', 'MS', 'NS' for each swath type
        # ##_[times]: 1d arrays of time components for ea scan:		[##_nScans]
        # ##_latitude:	2D arrays of latitudes for ea HS/MS/NS pixel:	[##_nPixels x ##_nScans]
        # ##_longitude:	2D arrays of longitudes for ea HS/MS/NS pixel:	[##_nPixels x ##_nScans]
        # ##_[most_fields]:  most fields are 2D arrays, vals @ ea pix:	[##_nPixels x ##_nScans]
        # ##_[3d_fields]:    the 3D fields have vals in ea BIN @ ea pix:[##_nBins  x ##_nPixels x ##_nScans]
        # ##_hts_m:  hts for the 3D fields have vals in ea BIN @ ea pix:[##_nPixels x ##_nScans x ##_nBins]
        
        
        # --- DETERMINING TIMESTAMPS: ---  **USING NS SWATH FOR TIMES!!**
        # -------------------------------
        # Determine the 2ADPR timestamp:  this will be the scan time for
        #  the pixel/scan that includes the main_plat's lat/lon location
        # loop thru each pixel row since dim of ##_[time_fields] arrays matches nScans dim
        # can (did in testing) get min dist for each HS, MS, and NS scan pixels, but for
        # simplicity am using only the NS scan for setting the timestamp, 2 reasons:
        #  - when only have the wider swath in OP over main_plat
        #  - when all 3 swaths OP main plat, the nScan subs for min dist will match (HS,MS,NS)
        scan_dists = np.full((NS_nPixels, NS_nscans2do), np.nan)
        for p in range(NS_nPixels):
            scan_lats = NS_latitude[:,p]
            scan_lons = NS_longitude[:,p]
          # get dist from ea pt in scan to the main_plat:
          #   scan_dist values are in [m]
            for s in range(NS_nscans2do):
                temp_result = sim.get_posn_to_azm_range(scan_lats[s], scan_lons[s], self.main_lat_deg, self.main_lon_deg)
                scan_dists[p,s] = temp_result[1]
        # Use the subscript for minimum dist to main_plat to assign timestamp:
        min_subs_ar = np.argwhere(scan_dists == np.min(scan_dists))
        min_dist_p_sub = min_subs_ar[0][0]   #for nPixels dim
        min_dist_s_sub = min_subs_ar[0][1]   #for  nScans dim
        print(f'main plat [p,s]: {min_dist_p_sub} {NS_start_scan+min_dist_s_sub}')
        #Reset the dpr_plat_info.timestamp:
        #time value arrays 1-d, NS_nScans dim
        dpr_plat_datetime = datetime(NS_year[min_dist_s_sub],
                                     NS_month[min_dist_s_sub],
                                     NS_day_of_month[min_dist_s_sub],
                                     NS_hour[min_dist_s_sub],
                                     NS_minute[min_dist_s_sub],
                                     NS_second[min_dist_s_sub])
        timestamp_year = dpr_plat_datetime.strftime('%Y')
        timestamp_month = dpr_plat_datetime.strftime('%m')
        timestamp_day = dpr_plat_datetime.strftime('%d')
        timestamp_hour = dpr_plat_datetime.strftime('%H')
        timestamp_min = dpr_plat_datetime.strftime('%M')
        timestamp_sec = dpr_plat_datetime.strftime('%S')

        self.dpr_info['timestamp'] = timestamp_year+timestamp_month+timestamp_day+'_'+ \
                                    timestamp_hour+timestamp_min+timestamp_sec
        print(f'...DPR timestamp - main platform: {self.dpr_info["timestamp"]}')
        offset = (dpr_plat_datetime - self.main_plat_datetime).total_seconds()
        self.dpr_info['offset_vs_main'] = offset
    
        
        #Determine the 2ADPR timestamp_cntr:
        #   time for column box center instead of main_plat
        #   do this same way did above for the main_plat location
        scan_dists = np.full((NS_nPixels, NS_nscans2do), np.nan)
        for p in range(NS_nPixels):
            scan_lats = NS_latitude[:,p]
            scan_lons = NS_longitude[:,p]
          # get distance from ea pt in scan to the main_plat:
          #  scan_dists values are in [m]
            for s in range(NS_nscans2do):
                temp_result = sim.get_posn_to_azm_range(scan_lats[s], scan_lons[s], self.cntr_lat_deg, self.cntr_lon_deg)
                scan_dists[p,s]  = temp_result[1]
                
        # Use the subscript for minimum dist to grid center to assign timestamp:
        min_subs_ar = np.argwhere(scan_dists == np.min(scan_dists))
        min_dist_p_sub = min_subs_ar[0][0]   #for nPixels dim
        min_dist_s_sub = min_subs_ar[0][1]   #for  nScans dim
        print(f'col center [p,s]: {min_dist_p_sub} {NS_start_scan+min_dist_s_sub}')
        # Reset the dpr_plat_info.timestamp_cntr:
        #  time value arrays 1-d, nScans dim
        dpr_cntr_datetime = datetime(NS_year[min_dist_s_sub],
                                     NS_month[min_dist_s_sub],
                                     NS_day_of_month[min_dist_s_sub],
                                     NS_hour[min_dist_s_sub],
                                     NS_minute[min_dist_s_sub],
                                     NS_second[min_dist_s_sub])
        timestamp_year = dpr_cntr_datetime.strftime('%Y')
        timestamp_month = dpr_cntr_datetime.strftime('%m')
        timestamp_day = dpr_cntr_datetime.strftime('%d')
        timestamp_hour = dpr_cntr_datetime.strftime('%H')
        timestamp_min = dpr_cntr_datetime.strftime('%M')
        timestamp_sec = dpr_cntr_datetime.strftime('%S')

        self.dpr_info['timestamp_cntr'] = timestamp_year+timestamp_month+timestamp_day+'_'+ \
                                          timestamp_hour+timestamp_min+timestamp_sec
        print(f'...DPR timestamp - col grid cntr: {self.dpr_info["timestamp_cntr"]}')

        
        # Define column box grid arrays to populate w/ 2ADPR data, initialize w/ all NaNs
        self.col_HS_sfcType              = np.full(self.plat_dims, np.nan)
        self.col_HS_flagPrecip           = np.full(self.plat_dims, np.nan)
        self.col_HS_binRealSfc           = np.full(self.plat_dims, np.nan)
        self.col_HS_htStormTop           = np.full(self.plat_dims, np.nan)
        self.col_HS_htZeroDeg            = np.full(self.plat_dims, np.nan)
        self.col_HS_flagBB               = np.full(self.plat_dims, np.nan)
        self.col_HS_htBB                 = np.full(self.plat_dims, np.nan)
        self.col_HS_widthBB              = np.full(self.plat_dims, np.nan)
        self.col_HS_qualityBB            = np.full(self.plat_dims, np.nan)
        self.col_HS_typePrecip           = np.full(self.plat_dims, np.nan)
        self.col_HS_qualityTypePrecip    = np.full(self.plat_dims, np.nan)
        self.col_HS_PIAfinal             = np.full(self.plat_dims, np.nan)
        self.col_HS_corZFacNearSfc       = np.full(self.plat_dims, np.nan)
        self.col_HS_precipRateNearSfc    = np.full(self.plat_dims, np.nan)
        self.col_HS_precipRateAve24      = np.full(self.plat_dims, np.nan)
        self.col_HS_phaseNearSfc         = np.full(self.plat_dims, np.nan)
        self.col_HS_zFacMeas             = np.full(self.plat_dims, np.nan)
        self.col_HS_attenNoPrecip        = np.full(self.plat_dims, np.nan)
        self.col_HS_corZFac              = np.full(self.plat_dims, np.nan)
        self.col_HS_precipRate           = np.full(self.plat_dims, np.nan)
        self.col_HS_DSDphase             = np.full(self.plat_dims, np.nan)
        self.col_HS_PIA_cloudwater       = np.full(self.plat_dims, np.nan)
        self.col_HS_PIA_cloudice         = np.full(self.plat_dims, np.nan)
        self.col_HS_PIA_watervapor       = np.full(self.plat_dims, np.nan)
        self.col_HS_PIA_oxygen           = np.full(self.plat_dims, np.nan)
        self.col_HS_TPW_liquid           = np.full(self.plat_dims, np.nan)
        self.col_HS_TPW_ice              = np.full(self.plat_dims, np.nan)
        self.col_HS_DSD_dBNw             = np.full(self.plat_dims, np.nan)
        self.col_HS_DSD_Dm               = np.full(self.plat_dims, np.nan)
        self.col_HS_eff_PIA              = np.full(self.plat_dims, np.nan)
        
        self.col_NS_sfcType              = np.full(self.plat_dims, np.nan)
        self.col_NS_flagPrecip           = np.full(self.plat_dims, np.nan)
        self.col_NS_binRealSfc           = np.full(self.plat_dims, np.nan)
        self.col_NS_htStormTop           = np.full(self.plat_dims, np.nan)
        self.col_NS_htZeroDeg            = np.full(self.plat_dims, np.nan)
        self.col_NS_flagBB               = np.full(self.plat_dims, np.nan)
        self.col_NS_htBB                 = np.full(self.plat_dims, np.nan)
        self.col_NS_widthBB              = np.full(self.plat_dims, np.nan)
        self.col_NS_qualityBB            = np.full(self.plat_dims, np.nan)
        self.col_NS_typePrecip           = np.full(self.plat_dims, np.nan)
        self.col_NS_qualityTypePrecip    = np.full(self.plat_dims, np.nan)
        self.col_NS_PIAfinal             = np.full(self.plat_dims, np.nan)
        self.col_NS_corZFacNearSfc       = np.full(self.plat_dims, np.nan)
        self.col_NS_precipRateNearSfc    = np.full(self.plat_dims, np.nan)
        self.col_NS_precipRateAve24      = np.full(self.plat_dims, np.nan)
        self.col_NS_phaseNearSfc         = np.full(self.plat_dims, np.nan)
        self.col_NS_zFacMeas             = np.full(self.plat_dims, np.nan)
        self.col_NS_attenNoPrecip        = np.full(self.plat_dims, np.nan)
        self.col_NS_corZFac              = np.full(self.plat_dims, np.nan)
        self.col_NS_precipRate           = np.full(self.plat_dims, np.nan)
        self.col_NS_DSDphase             = np.full(self.plat_dims, np.nan)
        self.col_NS_PIA_cloudwater       = np.full(self.plat_dims, np.nan)
        self.col_NS_PIA_cloudice         = np.full(self.plat_dims, np.nan)
        self.col_NS_PIA_watervapor       = np.full(self.plat_dims, np.nan)
        self.col_NS_PIA_oxygen           = np.full(self.plat_dims, np.nan)
        self.col_NS_TPW_liquid           = np.full(self.plat_dims, np.nan)
        self.col_NS_TPW_ice              = np.full(self.plat_dims, np.nan)
        self.col_NS_DSD_dBNw             = np.full(self.plat_dims, np.nan)
        self.col_NS_DSD_Dm               = np.full(self.plat_dims, np.nan)
        self.col_NS_eff_PIA              = np.full(self.plat_dims, np.nan)
        
        self.col_MS_sfcType              = np.full(self.plat_dims, np.nan)
        self.col_MS_flagPrecip           = np.full(self.plat_dims, np.nan)
        self.col_MS_binRealSfc           = np.full(self.plat_dims, np.nan)
        self.col_MS_htStormTop           = np.full(self.plat_dims, np.nan)
        self.col_MS_htZeroDeg            = np.full(self.plat_dims, np.nan)
        self.col_MS_flagBB               = np.full(self.plat_dims, np.nan)
        self.col_MS_htBB                 = np.full(self.plat_dims, np.nan)
        self.col_MS_widthBB              = np.full(self.plat_dims, np.nan)
        self.col_MS_qualityBB            = np.full(self.plat_dims, np.nan)
        self.col_MS_typePrecip           = np.full(self.plat_dims, np.nan)
        self.col_MS_qualityTypePrecip    = np.full(self.plat_dims, np.nan)
        self.col_MS_PIAfinal             = np.full(self.plat_dims, np.nan)
        self.col_MS_corZFacNearSfc       = np.full(self.plat_dims, np.nan)
        self.col_MS_precipRateNearSfc    = np.full(self.plat_dims, np.nan)
        self.col_MS_precipRateAve24      = np.full(self.plat_dims, np.nan)
        self.col_MS_phaseNearSfc         = np.full(self.plat_dims, np.nan)
        self.col_MS_zFacMeas             = np.full(self.plat_dims, np.nan)
        self.col_MS_attenNoPrecip        = np.full(self.plat_dims, np.nan)
        self.col_MS_corZFac              = np.full(self.plat_dims, np.nan)
        self.col_MS_PIA_cloudwater       = np.full(self.plat_dims, np.nan)
        self.col_MS_PIA_cloudice         = np.full(self.plat_dims, np.nan)
        self.col_MS_PIA_watervapor       = np.full(self.plat_dims, np.nan)
        self.col_MS_PIA_oxygen           = np.full(self.plat_dims, np.nan)
        self.col_MS_TPW_liquid           = np.full(self.plat_dims, np.nan)
        self.col_MS_TPW_ice              = np.full(self.plat_dims, np.nan)
        self.col_MS_eff_PIA              = np.full(self.plat_dims, np.nan)

        # --- LOCATING DPR DATA FOR THE COLUMN GRID: ---
        # ----------------------------------------------
        # Get x & y coords for ea horiz column grid point
        n_col_pts = self.n_horiz_grid_boxes +1
        begin_at = self.box_limit/2.0
        col_x_km = (np.arange(n_col_pts)*self.box_spacing - begin_at)/1000. #[km]
        col_y_km = (np.arange(n_col_pts)*self.box_spacing - begin_at)/1000. #[km]
        
        # Loop thru column grid horiz pts: for ea point,
        #  - determine which DPR point is nearest the col grid pt
        #  - use the DPR data for col grid, interpolated in vertical as needed
        # In order to populate column grid arrays with values from appropriate pixel(s):
        #   Assign values directly from the 2ADPR fields:
        #    1) In Horiz: set using the pixel/pt values as determined in code blocks above
        #    2) In Vert:  interpolate to column grid vert spacing, if applicable
        #   2-D fields:  set values at only lowest column grid level
        #   3-D fields:  set values at vertical levels using user-defined spacing
        #     2a) Define temp arrays to get subset of 3-D fields before re-arrange vert dim & vert. interp.
        #     2b) Do the vertical interp to column grid altitudes
        #   3) Get subset of hts for only the column box, and the binNumber for clutter free bottom
        
        #these arrays will be holders for 3D fields before doing the vert interp:
        dims = (self.lon_values.shape[0], self.lat_values.shape[0], HS_nBins)
        temp_HS_zFacMeas        = np.full(dims, np.nan)
        temp_HS_attenNoPrecip   = np.full(dims, np.nan)
        temp_HS_DSDphase        = np.full(dims, np.nan)
        temp_HS_corZFac         = np.full(dims, np.nan)
        temp_HS_precipRate      = np.full(dims, np.nan)
        temp_HS_DSD_dBNw        = np.full(dims, np.nan)
        temp_HS_DSD_Dm          = np.full(dims, np.nan)
        
        dims = (self.lon_values.shape[0], self.lat_values.shape[0], MS_nBins)
        temp_MS_zFacMeas        = np.full(dims, np.nan)
        temp_MS_attenNoPrecip   = np.full(dims, np.nan)
        temp_MS_corZFac         = np.full(dims, np.nan)
        
        dims = (self.lon_values.shape[0], self.lat_values.shape[0], NS_nBins)
        temp_NS_zFacMeas        = np.full(dims, np.nan)
        temp_NS_attenNoPrecip   = np.full(dims, np.nan)
        temp_NS_DSDphase        = np.full(dims, np.nan)
        temp_NS_corZFac         = np.full(dims, np.nan)
        temp_NS_precipRate      = np.full(dims, np.nan)
        temp_NS_DSD_dBNw        = np.full(dims, np.nan)
        temp_NS_DSD_Dm          = np.full(dims, np.nan)
        HS_hts_orig_m           = np.full((self.lon_values.shape[0], self.lat_values.shape[0], HS_nBins), np.nan)
        MS_hts_orig_m           = np.full((self.lon_values.shape[0], self.lat_values.shape[0], MS_nBins), np.nan)
        NS_hts_orig_m           = np.full((self.lon_values.shape[0], self.lat_values.shape[0], NS_nBins), np.nan)
        
        dims = (self.lon_values.shape[0], self.lat_values.shape[0])
        HS_binNum_CFB           = np.full(dims, np.nan)
        MS_binNum_CFB           = np.full(dims, np.nan)
        NS_binNum_CFB           = np.full(dims, np.nan)

        for i, val in enumerate(self.lon_values):
            for j, val2 in enumerate(self.lat_values):
        #for i, val in enumerate(self.column_box_params['column_grid_lats']):
        #    for j, val2 in enumerate(self.column_box_params['column_grid_lons']):
                NS_dpr2cpt_dist = np.sqrt( (NS_longitude_subset - self.lon_values[i])**2 + (NS_latitude_subset - self.lat_values[j])**2)
                
                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    MS_dpr2cpt_dist = np.sqrt( (MS_longitude_subset - self.lon_values[i])**2 + (MS_latitude_subset - self.lat_values[j])**2)
                    
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    HS_dpr2cpt_dist = np.sqrt( (HS_longitude_subset - self.lon_values[i])**2 + (HS_latitude_subset - self.lat_values[j])**2)

                #First, get [p,s] for HS swath:
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    HS_min_pt_dist_arr = np.argwhere(HS_dpr2cpt_dist == np.min(HS_dpr2cpt_dist))
                    HS_pt_p_sub    = HS_min_pt_dist_arr[1]
                    HS_pt_s_sub    = HS_min_pt_dist_arr[0]

                # print, '---column xy: ',col_x_km[i],col_y_km[j]
                # print, '  pixel,scan: ',HS_pt_p_sub,HS_start_scan+HS_pt_s_sub
                # print, 'min dpr dist: ',HS_min_pt_dist
                
                #Repeat as above for MS swath:
                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    MS_min_pt_dist_arr = np.argwhere(MS_dpr2cpt_dist == np.min(MS_dpr2cpt_dist))
                    MS_pt_p_sub    = MS_min_pt_dist_arr[0][1]
                    MS_pt_s_sub    = MS_min_pt_dist_arr[0][0]

                #Repeat as above for NS swath:
                NS_min_pt_dist_arr = np.argwhere(NS_dpr2cpt_dist == np.min(NS_dpr2cpt_dist))
                NS_pt_p_sub = NS_min_pt_dist_arr[0][1]   
                NS_pt_s_sub = NS_min_pt_dist_arr[0][0]
        
                #-- Now that we have the [s,p] subs for this column grid pt,
                #--  assign values to column arrays from this [s,p] spot in 3 swaths:
                # 1) for horiz dims:
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    self.col_HS_sfcType[0,j,i]           = HS_sfc_type[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_flagPrecip[0,j,i]        = HS_flag_precip[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_binRealSfc[0,j,i]        = HS_binRealSfc[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_htStormTop[0,j,i]        = HS_ht_stormTop[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_htZeroDeg[0,j,i]         = HS_ht_ZeroDeg[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_flagBB[0,j,i]            = HS_flag_BB[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_htBB[0,j,i]              = HS_ht_BB[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_widthBB[0,j,i]           = HS_width_BB[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_qualityBB[0,j,i]         = HS_quality_BB[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_typePrecip[0,j,i]        = HS_typePrecip[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_qualityTypePrecip[0,j,i] = HS_qualityTypePrecip[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_PIAfinal[0,j,i]          = HS_PIA_final[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_corZFacNearSfc[0,j,i]    = HS_corZFac_nearSfc[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_precipRateNearSfc[0,j,i] = HS_precipRate_nearSfc[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_precipRateAve24[0,j,i]   = HS_precipRateAve24[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_phaseNearSfc[0,j,i]      = HS_phase_nearSfc[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_PIA_cloudwater[0,j,i]    = HS_PIA_cloudwater[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_PIA_cloudice[0,j,i]      = HS_PIA_cloudice[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_PIA_watervapor[0,j,i]    = HS_PIA_watervapor[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_PIA_oxygen[0,j,i]        = HS_PIA_oxygen[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_TPW_liquid[0,j,i]        = HS_TPW_liquid[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_TPW_ice[0,j,i]           = HS_TPW_ice[HS_pt_s_sub, HS_pt_p_sub]
                    self.col_HS_eff_PIA[0,j,i]           = HS_eff_PIA[HS_pt_s_sub, HS_pt_p_sub]
        
                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    self.col_MS_sfcType[0,j,i]           = MS_sfc_type[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_flagPrecip[0,j,i]        = MS_flag_precip[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_binRealSfc[0,j,i]        = MS_binRealSfc[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_htStormTop[0,j,i]        = MS_ht_stormTop[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_htZeroDeg[0,j,i]         = MS_ht_ZeroDeg[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_flagBB[0,j,i]            = MS_flag_BB[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_htBB[0,j,i]              = MS_ht_BB[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_widthBB[0,j,i]           = MS_width_BB[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_qualityBB[0,j,i]         = MS_quality_BB[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_typePrecip[0,j,i]        = MS_typePrecip[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_qualityTypePrecip[0,j,i] = MS_qualityTypePrecip[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_PIAfinal[0,j,i]          = MS_PIA_final[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_corZFacNearSfc[0,j,i]    = MS_corZFac_nearSfc[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_precipRateNearSfc[0,j,i] = MS_precipRate_nearSfc[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_precipRateAve24[0,j,i]   = MS_precipRateAve24[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_phaseNearSfc[0,j,i]      = MS_phase_nearSfc[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_PIA_cloudwater[0,j,i]    = MS_PIA_cloudwater[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_PIA_cloudice[0,j,i]      = MS_PIA_cloudice[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_PIA_watervapor[0,j,i]    = MS_PIA_watervapor[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_PIA_oxygen[0,j,i]        = MS_PIA_oxygen[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_TPW_liquid[0,j,i]        = MS_TPW_liquid[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_TPW_ice[0,j,i]           = MS_TPW_ice[MS_pt_s_sub, MS_pt_p_sub]
                    self.col_MS_eff_PIA[0,j,i]           = MS_eff_PIA[MS_pt_s_sub, MS_pt_p_sub]
          
                self.col_NS_sfcType[0,j,i]             = NS_sfc_type[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_flagPrecip[0,j,i]          = NS_flag_precip[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_binRealSfc[0,j,i]          = NS_binRealSfc[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_htStormTop[0,j,i]          = NS_ht_stormTop[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_htZeroDeg[0,j,i]           = NS_ht_ZeroDeg[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_flagBB[0,j,i]              = NS_flag_BB[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_htBB[0,j,i]                = NS_ht_BB[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_widthBB[0,j,i]             = NS_width_BB[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_qualityBB[0,j,i]           = NS_quality_BB[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_typePrecip[0,j,i]          = NS_typePrecip[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_qualityTypePrecip[0,j,i]   = NS_qualityTypePrecip[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_PIAfinal[0,j,i]            = NS_PIA_final[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_corZFacNearSfc[0,j,i]      = NS_corZFac_nearSfc[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_precipRateNearSfc[0,j,i]   = NS_precipRate_nearSfc[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_precipRateAve24[0,j,i]     = NS_precipRateAve24[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_phaseNearSfc[0,j,i]        = NS_phase_nearSfc[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_PIA_cloudwater[0,j,i]      = NS_PIA_cloudwater[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_PIA_cloudice[0,j,i]        = NS_PIA_cloudice[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_PIA_watervapor[0,j,i]      = NS_PIA_watervapor[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_PIA_oxygen[0,j,i]          = NS_PIA_oxygen[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_TPW_liquid[0,j,i]          = NS_TPW_liquid[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_TPW_ice[0,j,i]             = NS_TPW_ice[NS_pt_s_sub, NS_pt_p_sub]
                self.col_NS_eff_PIA[0,j,i]             = NS_eff_PIA[NS_pt_s_sub, NS_pt_p_sub]

                # 2a) set horiz dims for 3D fields before vert interp:
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    temp_HS_zFacMeas[j,i,:]         = HS_zFac_meas[HS_pt_s_sub, HS_pt_p_sub,:]
                    temp_HS_attenNoPrecip[j,i,:]    = HS_atten_noPrecip[HS_pt_s_sub, HS_pt_p_sub,:]
                    temp_HS_DSDphase[j,i,:]         = HS_DSD_phase[HS_pt_s_sub, HS_pt_p_sub,:]
                    temp_HS_corZFac[j,i,:]          = HS_corZFac[HS_pt_s_sub, HS_pt_p_sub,:]
                    temp_HS_precipRate[j,i,:]       = HS_precipRate[HS_pt_s_sub, HS_pt_p_sub,:]
                    temp_HS_DSD_dBNw[j,i,:]         = HS_DSD_dBNw[HS_pt_s_sub, HS_pt_p_sub,:]
                    temp_HS_DSD_Dm[j,i,:]           = HS_DSD_Dm[HS_pt_s_sub, HS_pt_p_sub,:]
                #pdb.set_trace()
                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    temp_MS_zFacMeas[j,i,:]         = MS_zFac_meas[MS_pt_s_sub, MS_pt_p_sub,:]
                    temp_MS_attenNoPrecip[j,i,:]    = MS_atten_noPrecip[MS_pt_s_sub, MS_pt_p_sub,:]
                    temp_MS_corZFac[j,i,:]          = MS_corZFac[MS_pt_s_sub, MS_pt_p_sub,:]

                temp_NS_zFacMeas[j,i,:]         = NS_zFac_meas[NS_pt_s_sub, NS_pt_p_sub,:]
                temp_NS_attenNoPrecip[j,i,:]    = NS_atten_noPrecip[NS_pt_s_sub, NS_pt_p_sub,:]
                temp_NS_DSDphase[j,i,:]         = NS_DSD_phase[NS_pt_s_sub, NS_pt_p_sub,:]
                temp_NS_corZFac[j,i,:]          = NS_corZFac[NS_pt_s_sub, NS_pt_p_sub,:]
                temp_NS_precipRate[j,i,:]       = NS_precipRate[NS_pt_s_sub, NS_pt_p_sub,:]
                temp_NS_DSD_dBNw[j,i,:]         = NS_DSD_dBNw[NS_pt_s_sub, NS_pt_p_sub,:]
                temp_NS_DSD_Dm[j,i,:]           = NS_DSD_Dm[NS_pt_s_sub, NS_pt_p_sub,:]
        
                #3) get hts arrays for column box points:
                NS_hts_orig_m[j,i,:]    = NS_hts_m[NS_pt_s_sub, NS_pt_p_sub,:]
                NS_binNum_CFB[j,i]      = NS_binClutFB[NS_pt_s_sub, NS_pt_p_sub]
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    HS_hts_orig_m[j,i,:]    = HS_hts_m[HS_pt_s_sub, HS_pt_p_sub, :]
                    HS_binNum_CFB[j,i]      = HS_binClutFB[HS_pt_s_sub, HS_pt_p_sub]
                else:
                    #HS swath does not include column grid:
                    HS_binNum_CFB[j,i]  = -1

                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    MS_hts_orig_m[j,i,:]    = MS_hts_m[MS_pt_s_sub, MS_pt_p_sub,:]
                    MS_binNum_CFB[j,i]      = MS_binClutFB[MS_pt_s_sub, MS_pt_p_sub]
                else:
                    # MS swath does not include column grid:
                    MS_binNum_CFB[i,j]      = -1

        HS_binNum_CFB = HS_binNum_CFB.astype(int)
        MS_binNum_CFB = MS_binNum_CFB.astype(int)
        NS_binNum_CFB = NS_binNum_CFB.astype(int)
        
        # Set NaNs for 3D field values at hts below clutter-free area:
        # bin clutter-free bottom:  is the bin NUMBER of last clutter-free gate away from DPR
        #			    bigger number = lower altitude/closer to earth
        # set all bin SUBSCRIPTS >= bin clutter free bottom  to NaN values for profile fields
        #	bin subscript  = (bin number -1)
        #	bin numbers    >  binClutFB  -> set to NaN
        # 	bin subscripts >= binClutFB  -> set to NaN
        #	so can use vertical dim subscript range:
        #	 [binClutFB[i,j]:*] bc last subs will be the largest bin Number (have not flipped hts order yet...)
        #for i=0,n_elements(col_x_km)-1 do begin
        #for j=0,n_elements(col_y_km)-1 do begin
        for i, val in enumerate(self.lon_values):
            for j, val2 in enumerate(self.lat_values):
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    temp_HS_zFacMeas[j,i,HS_binNum_CFB[j,i]:]       = np.nan
                    temp_HS_attenNoPrecip[j,i,HS_binNum_CFB[j,i]:]  = np.nan
                    temp_HS_DSDphase[j,i,HS_binNum_CFB[j,i]:]       = np.nan
                    temp_HS_corZFac[j,i,HS_binNum_CFB[j,i]:]        = np.nan
                    temp_HS_precipRate[j,i,HS_binNum_CFB[j,i]:]     = np.nan
                    temp_HS_DSD_dBNw[j,i,HS_binNum_CFB[j,i]:]       = np.nan
                    temp_HS_DSD_Dm[j,i,HS_binNum_CFB[j,i]:]         = np.nan

                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    #pdb.set_trace()
                    temp_MS_zFacMeas[j,i,MS_binNum_CFB[j,i]:]       = np.nan
                    temp_MS_attenNoPrecip[j,i,MS_binNum_CFB[j,i]:]  = np.nan
                    temp_MS_corZFac[j,i,MS_binNum_CFB[j,i]:]        = np.nan
                    
                temp_NS_zFacMeas[j,i,NS_binNum_CFB[j,i]:]           = np.nan
                temp_NS_attenNoPrecip[j,i,NS_binNum_CFB[j,i]:]      = np.nan
                temp_NS_DSDphase[j,i,NS_binNum_CFB[j,i]:]           = np.nan
                temp_NS_corZFac[j,i,NS_binNum_CFB[j,i]:]            = np.nan
                temp_NS_precipRate[j,i,NS_binNum_CFB[j,i]:]         = np.nan
                temp_NS_DSD_dBNw[j,i,NS_binNum_CFB[j,i]:]           = np.nan
                temp_NS_DSD_Dm[j,i,NS_binNum_CFB[j,i]:]             = np.nan
        #pdb.set_trace()
        # Set NaNs for missing data values
        if HS_nscans2do > 0 and col_in_HSswath > 0:
            self.col_HS_sfcType[self.col_HS_sfcType <= badVal]                    = np.nan
            self.col_HS_flagPrecip[self.col_HS_flagPrecip <= badVal]              = np.nan
            self.col_HS_binRealSfc[self.col_HS_binRealSfc <= badVal]              = np.nan
            self.col_HS_htStormTop[self.col_HS_htStormTop <= badVal]              = np.nan
            self.col_HS_htZeroDeg[self.col_HS_htZeroDeg <= badVal]                = np.nan
            self.col_HS_flagBB[self.col_HS_flagBB <= badVal]                      = np.nan
            self.col_HS_htBB[self.col_HS_htBB <= badVal]                          = np.nan
            self.col_HS_widthBB[self.col_HS_widthBB <= badVal]                    = np.nan
            self.col_HS_qualityBB[self.col_HS_qualityBB <= badVal]                = np.nan
            self.col_HS_typePrecip[self.col_HS_typePrecip <= badVal]              = np.nan
            self.col_HS_qualityTypePrecip[self.col_HS_qualityTypePrecip <= badVal]= np.nan
            self.col_HS_PIAfinal[self.col_HS_PIAfinal <= badVal]                  = np.nan
            self.col_HS_corZFacNearSfc[self.col_HS_corZFacNearSfc <= badVal]      = np.nan
            self.col_HS_precipRateNearSfc[self.col_HS_precipRateNearSfc <= badVal]= np.nan
            self.col_HS_precipRateAve24[self.col_HS_precipRateAve24 <= badVal]    = np.nan
            self.col_HS_phaseNearSfc[self.col_HS_phaseNearSfc == 255]             = np.nan
            self.col_HS_PIA_cloudwater[self.col_HS_PIA_cloudwater <= badVal]      = np.nan
            self.col_HS_PIA_cloudice[self.col_HS_PIA_cloudice <= badVal]          = np.nan
            self.col_HS_PIA_watervapor[self.col_HS_PIA_watervapor <= badVal]      = np.nan
            self.col_HS_PIA_oxygen[self.col_HS_PIA_oxygen <= badVal]              = np.nan
            self.col_HS_TPW_liquid[self.col_HS_TPW_liquid <= badVal]              = np.nan
            self.col_HS_TPW_ice[self.col_HS_TPW_ice <= badVal]                    = np.nan
            self.col_HS_eff_PIA[self.col_HS_eff_PIA <= badVal]                    = np.nan
            temp_HS_zFacMeas[temp_HS_zFacMeas <= badVal]                = np.nan
            temp_HS_attenNoPrecip[temp_HS_attenNoPrecip <= badVal]      = np.nan
            temp_HS_corZFac[temp_HS_corZFac <= badVal]                  = np.nan
            temp_HS_precipRate[temp_HS_precipRate <= badVal]            = np.nan
            temp_HS_DSDphase[temp_HS_DSDphase == 255]                   = np.nan
            temp_HS_DSD_dBNw[temp_HS_DSD_dBNw <= badVal]                = np.nan
            temp_HS_DSD_Dm[temp_HS_DSD_Dm <= badVal]                    = np.nan

        if MS_nscans2do > 0 and col_in_MSswath > 0:
            self.col_MS_sfcType[self.col_MS_sfcType <= badVal]                    = np.nan
            self.col_MS_flagPrecip[self.col_MS_flagPrecip <= badVal]              = np.nan
            self.col_MS_binRealSfc[self.col_MS_binRealSfc <= badVal]              = np.nan
            self.col_MS_htStormTop[self.col_MS_htStormTop <= badVal]              = np.nan
            self.col_MS_htZeroDeg[self.col_MS_htZeroDeg <= badVal]                = np.nan
            self.col_MS_flagBB[self.col_MS_flagBB <= badVal]                      = np.nan
            self.col_MS_htBB[self.col_MS_htBB <= badVal]                          = np.nan
            self.col_MS_widthBB[self.col_MS_widthBB <= badVal]                    = np.nan
            self.col_MS_qualityBB[self.col_MS_qualityBB <= badVal]                = np.nan
            self.col_MS_typePrecip[self.col_MS_typePrecip <= badVal]              = np.nan
            self.col_MS_qualityTypePrecip[self.col_MS_qualityTypePrecip <= badVal]= np.nan
            self.col_MS_PIAfinal[self.col_MS_PIAfinal <= badVal]                  = np.nan
            self.col_MS_corZFacNearSfc[self.col_MS_corZFacNearSfc <= badVal]      = np.nan
            self.col_MS_precipRateNearSfc[self.col_MS_precipRateNearSfc <= badVal]= np.nan
            self.col_MS_precipRateAve24[self.col_MS_precipRateAve24 <= badVal]    = np.nan
            self.col_MS_phaseNearSfc[self.col_MS_phaseNearSfc == 255]             = np.nan
            self.col_MS_PIA_cloudwater[self.col_MS_PIA_cloudwater <= badVal]      = np.nan
            self.col_MS_PIA_cloudice[self.col_MS_PIA_cloudice <= badVal]          = np.nan
            self.col_MS_PIA_watervapor[self.col_MS_PIA_watervapor <= badVal]      = np.nan
            self.col_MS_PIA_oxygen[self.col_MS_PIA_oxygen <= badVal]              = np.nan
            self.col_MS_TPW_liquid[self.col_MS_TPW_liquid <= badVal]              = np.nan
            self.col_MS_TPW_ice[self.col_MS_TPW_ice <= badVal]                    = np.nan
            self.col_MS_eff_PIA[self.col_MS_eff_PIA <= badVal]                    = np.nan
            temp_MS_zFacMeas[temp_MS_zFacMeas <= badVal]                = np.nan
            temp_MS_attenNoPrecip[temp_MS_attenNoPrecip <= badVal]      = np.nan
            temp_MS_corZFac[temp_MS_corZFac <= badVal]                  = np.nan

        self.col_NS_sfcType[self.col_NS_sfcType <= badVal]                        = np.nan
        self.col_NS_flagPrecip[self.col_NS_flagPrecip <= badVal]                  = np.nan
        self.col_NS_binRealSfc[self.col_NS_binRealSfc <= badVal]                  = np.nan
        self.col_NS_htStormTop[self.col_NS_htStormTop <= badVal]                  = np.nan
        self.col_NS_htZeroDeg[self.col_NS_htZeroDeg <= badVal]                    = np.nan
        self.col_NS_flagBB[self.col_NS_flagBB <= badVal]                          = np.nan
        self.col_NS_htBB[self.col_NS_htBB <= badVal]                              = np.nan
        self.col_NS_widthBB[self.col_NS_widthBB <= badVal]                        = np.nan
        self.col_NS_qualityBB[self.col_NS_qualityBB <= badVal]                    = np.nan
        self.col_NS_typePrecip[self.col_NS_typePrecip <= badVal]                  = np.nan
        self.col_NS_qualityTypePrecip[self.col_NS_qualityTypePrecip <= badVal]    = np.nan
        self.col_NS_PIAfinal[self.col_NS_PIAfinal <= badVal]                      = np.nan
        self.col_NS_corZFacNearSfc[self.col_NS_corZFacNearSfc <= badVal]          = np.nan
        self.col_NS_precipRateNearSfc[self.col_NS_precipRateNearSfc <= badVal]    = np.nan
        self.col_NS_precipRateAve24[self.col_NS_precipRateAve24 <= badVal]        = np.nan
        self.col_NS_phaseNearSfc[self.col_NS_phaseNearSfc == 255]                 = np.nan
        self.col_NS_PIA_cloudwater[self.col_NS_PIA_cloudwater <= badVal]          = np.nan
        self.col_NS_PIA_cloudice[self.col_NS_PIA_cloudice <= badVal]              = np.nan
        self.col_NS_PIA_watervapor[self.col_NS_PIA_watervapor <= badVal]           = np.nan
        self.col_NS_PIA_oxygen[self.col_NS_PIA_oxygen <= badVal]                  = np.nan
        self.col_NS_TPW_liquid[self.col_NS_TPW_liquid <= badVal]                  = np.nan
        self.col_NS_TPW_ice[self.col_NS_TPW_ice <= badVal]                        = np.nan
        self.col_NS_eff_PIA[self.col_NS_eff_PIA <= badVal]                        = np.nan
        temp_NS_zFacMeas[temp_NS_zFacMeas <= badVal]                    = np.nan
        temp_NS_attenNoPrecip[temp_NS_attenNoPrecip <= badVal]          = np.nan
        temp_NS_corZFac[temp_NS_corZFac <= badVal]                      = np.nan
        temp_NS_precipRate[temp_NS_precipRate <= badVal]                = np.nan
        temp_NS_DSDphase[temp_NS_DSDphase == 255]                       = np.nan
        temp_NS_DSD_dBNw[temp_NS_DSD_dBNw <= badVal]                    = np.nan
        temp_NS_DSD_Dm[temp_NS_DSD_Dm <= badVal]                        = np.nan

        # Must flip order of vert dim: in data, bin0 = top, need bin0 = ht 0
        if HS_nscans2do > 0 and col_in_HSswath > 0:
            temp_HS_zFacMeas        = temp_HS_zFacMeas[:,:,::-1]
            temp_HS_attenNoPrecip   = temp_HS_attenNoPrecip[:,:,::-1]
            temp_HS_corZFac         = temp_HS_corZFac[:,:,::-1]
            temp_HS_precipRate      = temp_HS_precipRate[:,:,::-1]
            temp_HS_DSDphase        = temp_HS_DSDphase[:,:,::-1]
            temp_HS_DSD_dBNw        = temp_HS_DSD_dBNw[:,:,::-1]
            temp_HS_DSD_Dm          = temp_HS_DSD_Dm[:,:,::-1]
            HS_hts_orig_m           = HS_hts_orig_m[:,:,::-1]
        if MS_nscans2do > 0 and col_in_MSswath > 0:
            temp_MS_zFacMeas        = temp_MS_zFacMeas[:,:,::-1]
            temp_MS_attenNoPrecip   = temp_MS_attenNoPrecip[:,:,::-1]
            temp_MS_corZFac         = temp_MS_corZFac[:,:,::-1]
            MS_hts_orig_m           = MS_hts_orig_m[:,:,::-1]
        temp_NS_zFacMeas            = temp_NS_zFacMeas[:,:,::-1]
        temp_NS_attenNoPrecip       = temp_NS_attenNoPrecip[:,:,::-1]
        temp_NS_corZFac             = temp_NS_corZFac[:,:,::-1]
        temp_NS_precipRate          = temp_NS_precipRate[:,:,::-1]
        temp_NS_DSDphase            = temp_NS_DSDphase[:,:,::-1]
        temp_NS_DSD_dBNw            = temp_NS_DSD_dBNw[:,:,::-1]
        temp_NS_DSD_Dm              = temp_NS_DSD_Dm[:,:,::-1]
        NS_hts_orig_m               = NS_hts_orig_m[:,:,::-1]
        
        # For 3D data arrays: Use scipy interpolate function to get val from
        # DPR HS/NS/MS hts to column grid hts (Step 2b from above):
        # (send 1-D arrays to interpol, so need loops)
        #col_hts = self.column_box_params['column_grid_alts']    # [m]
        #for i=0,n_elements(col_x_km)-1 do begin
        #for j=0,n_elements(col_y_km)-1 do begin
        for i, val in enumerate(self.lon_values):
            for j, val2 in enumerate(self.lat_values):
                if HS_nscans2do > 0 and col_in_HSswath > 0:
                    self.col_HS_zFacMeas[:,j,i]      = sim.interp(temp_HS_zFacMeas[j,i,:],     HS_hts_orig_m[j,i,:], self.z_values)
                    self.col_HS_attenNoPrecip[:,j,i] = sim.interp(temp_HS_attenNoPrecip[j,i,:],HS_hts_orig_m[j,i,:], self.z_values)
                    self.col_HS_corZFac[:,j,i]       = sim.interp(temp_HS_corZFac[j,i,:],      HS_hts_orig_m[j,i,:], self.z_values)
                    self.col_HS_precipRate[:,j,i]    = sim.interp(temp_HS_precipRate[j,i,:],   HS_hts_orig_m[j,i,:], self.z_values)
                    self.col_HS_DSDphase[:,j,i]      = sim.interp(temp_HS_DSDphase[j,i,:],     HS_hts_orig_m[j,i,:], self.z_values)
                    self.col_HS_DSD_dBNw[:,j,i]      = sim.interp(temp_HS_DSD_dBNw[j,i,:],     HS_hts_orig_m[j,i,:], self.z_values)
                    self.col_HS_DSD_Dm[:,j,i]        = sim.interp(temp_HS_DSD_Dm[j,i,:],       HS_hts_orig_m[j,i,:], self.z_values)
                if MS_nscans2do > 0 and col_in_MSswath > 0:
                    self.col_MS_zFacMeas[:,j,i]      = sim.interp(temp_MS_zFacMeas[j,i,:],     MS_hts_orig_m[j,i,:], self.z_values)
                    self.col_MS_attenNoPrecip[:,j,i] = sim.interp(temp_MS_attenNoPrecip[j,i,:],MS_hts_orig_m[j,i,:], self.z_values)
                    self.col_MS_corZFac[:,j,i]       = sim.interp(temp_MS_corZFac[j,i,:],      MS_hts_orig_m[j,i,:], self.z_values)
                    
                self.col_NS_zFacMeas[:,j,i]          = sim.interp(temp_NS_zFacMeas[j,i,:],     NS_hts_orig_m[j,i,:], self.z_values)
                self.col_NS_attenNoPrecip[:,j,i]     = sim.interp(temp_NS_attenNoPrecip[j,i,:],NS_hts_orig_m[j,i,:], self.z_values)
                self.col_NS_corZFac[:,j,i]           = sim.interp(temp_NS_corZFac[j,i,:],      NS_hts_orig_m[j,i,:], self.z_values)
                self.col_NS_precipRate[:,j,i]        = sim.interp(temp_NS_precipRate[j,i,:],   NS_hts_orig_m[j,i,:], self.z_values)
                self.col_NS_DSDphase[:,j,i]          = sim.interp(temp_NS_DSDphase[j,i,:],     NS_hts_orig_m[j,i,:], self.z_values)
                self.col_NS_DSD_dBNw[:,j,i]          = sim.interp(temp_NS_DSD_dBNw[j,i,:],     NS_hts_orig_m[j,i,:], self.z_values)
                self.col_NS_DSD_Dm[:,j,i]            = sim.interp(temp_NS_DSD_Dm[j,i,:],       NS_hts_orig_m[j,i,:], self.z_values)
        
        #remove 10*log(Nw) multiplier from Nw fields:
        if HS_nscans2do > 0 and col_in_HSswath > 0: self.col_HS_DSD_dBNw = self.col_HS_DSD_dBNw/10.0
        self.col_NS_DSD_dBNw = self.col_NS_DSD_dBNw/10.0
        
        # Define units, more descriptive names for each field:
        self.dpr_name_sfcType    = 'land sruface type category'
        self.dpr_units_sfcType   = 'int value for category'
        
        self.dpr_name_flagPrecip     = 'flag if precip or no precip'
        self.dpr_units_flagPrecip    = '0: no precip  11: precip'
        
        self.dpr_name_htStormTop     = 'ht of storm top'
        self.dpr_units_htStormTop    = 'm'
        
        self.dpr_name_htZeroDeg      = 'ht of freezing level'
        self.dpr_units_htZeroDeg     = 'm'
        
        self.dpr_name_flagBB     = 'tells bright band existance'
        self.dpr_units_flagBB    = '0: no BB  1: yes BB'
        
        self.dpr_name_htBB   = 'ht of bright band'
        self.dpr_units_htBB  = 'm'
        
        self.dpr_name_widthBB    = 'width of bright band'
        self.dpr_units_widthBB   = 'm'
        
        self.dpr_name_qualityBB  = 'how well defined bright band is'
        self.dpr_units_qualityBB = '1: best'
        
        self.dpr_name_typePrecip  = '8 char precip description - 1st char tells major type'
        self.dpr_units_typePrecip = '1:stratiform  2:convective 3:other'
        
        self.dpr_name_qualityTypePrecip  = 'quality of precip type'
        self.dpr_units_qualityTypePrecip = '1: best'
        
        self.dpr_name_PIAfinal  = 'final est of path integrated attenuation due to precip particles'
        self.dpr_units_PIAfinal = 'dB'
        
        self.dpr_name_corZFacNearSfc  = 'reflectivity factor corrected for attenuation near surface'
        self.dpr_units_corZFacNearSfc = 'dBZ'
        
        self.dpr_name_precipRateNearSfc  = 'precip rate near the surface'
        self.dpr_units_precipRateNearSfc = 'mm h^-1'
        
        self.dpr_name_precipRateAve24  = 'average precip rate for 2-4 km ht'
        self.dpr_units_precipRateAve24 = 'mm h^-1'
        
        self.dpr_name_phaseNearSfc  = 'phase of precip near surface'
        self.dpr_units_phaseNearSfc = 'for int value/100:  0:solid  1:mixed  2:liquid'
        
        self.dpr_name_zFacMeas  = 'reflectivity factor with no attenuation correction'
        self.dpr_units_zFacMeas = 'dBZ'
        
        self.dpr_name_attenNoPrecip  = 'attenuation by non-precip particles'
        self.dpr_units_attenNoPrecip = 'dB km^-1'
        
        self.dpr_name_corZFac  = 'reflectivity factor corrected for attenuation'
        self.dpr_units_corZFac = 'dBZ'
        
        self.dpr_name_precipRate  = 'precipitation rate'
        self.dpr_units_precipRate = 'mm h^-1'
        
        self.dpr_name_DSDphase  = 'phase of precip'
        self.dpr_units_DSDphase = 'for int value/100:  0:solid  1:mixed  2:liquid'
        
        self.dpr_name_PIA_cloudwater  = 'PIA due to non precip: cloud water'
        self.dpr_units_PIA_cloudwater = 'dB'
        
        self.dpr_name_PIA_cloudice  = 'PIA due to non precip: cloud ice'
        self.dpr_units_PIA_cloudice = 'dB'
        
        self.dpr_name_PIA_watervapor  = 'PIA due to non precip: water vapor'
        self.dpr_units_PIA_watervapor = 'dB'
        
        self.dpr_name_PIA_oxygen  = 'PIA due to non precip: oxygen molecules'
        self.dpr_units_PIA_oxygen = 'dB'
        
        self.dpr_name_TPW_liquid = 'vertically integrated liquid precip water'
        self.dpr_units_TPW_liquid= 'g m^-2'
        
        self.dpr_name_TPW_ice    = 'vertically integratied solid precip water'
        self.dpr_units_TPW_ice   = 'g m^-2'
        
        self.dpr_name_DSD_dBNw   = 'DSD normalized intercept parameter'
        self.dpr_units_DSD_dBNw  = 'log(Nw)'
        
        self.dpr_name_DSD_Dm = 'DSD mass-weighted mean diameter'
        self.dpr_units_DSD_Dm    = 'mm'
        
        self.dpr_name_effPIA	= 'effective 2-way PIA'
        self.dpr_units_effPIA	= 'dB'
        
        self.dpr_data_in_column = True

    # ***************************************************************************************
    def get_mrr_for_column(self):
        '''
        ;------------------------------------------------------------
        ; read in MRR data and extract vertical profiles at the main_plat time
        ; +/- halftime_interval & within column box, return data structure
        ;------------------------------------------------------------
        ;
        ; -> MRR location MUST BE CORRECT in platform_location.pro
        ;
        ; IMPORTANT!! because there is NO LOCATION info in MRR data files, is assumed
        ;  	     that if files_mrr_ave paths are set in build_column.pro the MRR(s)
        ;	     is/are located within the column box grid.
        ;
        ;   result = get_mrr_for_column(column_box_params, mrr_IDs, files_mrr, $
        ;					main_plat_timestamp, halftime_interval)
        ;	column_box_params:   column box grid parameters. This is the 
        ;				structure returned by define_box.pro
        ;	mrr_IDs:	     string array of MRR unit name(s) (as 'MRR-0#')
        ;	files_mrr:	     string array of paths to input MRR .ave file(s)
        ;				NOTE:  corresponding MRR ID & MRR .ave files
        ;				       MUST be placed at same subscripts in 
        ;				       mrr_IDs & files_mrr arrays!
        ;	main_plat_info:      main_plat_info structure returned by 
        ;				set_main_plat_values.pro
        ;	halftime_interval:   integer # of mins before & after the main
        ;				platform time to include in column file
        ;
        ;	result:	 returns a big structure:
        ;		mrr_info_and_data:  This is made up of multiple structures:
        ;
        ;		.mrr_plat_info:	structure of MRR's plat_info values
        ;				mainly set up in set_plat_values.pro, but
        ;				timestamp & offset tags are set in here
        ;				(if > 1 MRR unit in column box, ea tag
        ;				 is really an array)
        ;
        ;		.col_data_mrr:	structure of MRR(s) data in the column grid,
        ;				and units/info tags describing each field
        ;
        ;	   (also see more/full details in READMEs)
        ;
        ;	Dependencies:
        ;	  dms2dd.pro	utility to convert deg-min-sec to decimal degrees
        ;	  closest.pro	utility to locate value in an array nearest input value
        ;	  read_data_mrr.pro	reader for reading in the MRR .ave data
        ;	  disdro_mrr_get_dsd_params   these 2 functions by Patrick Gatlin 
        ;	  disdro_fit_dsd_norm_gamma    compute Dm,Nw,& Z from MRR DSD profiles
        ;
        ;------------------------------------------------------------
        '''
    
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
        self.mrr_PIA = np.full(self.ground_plat_dims, np.nan)
        self.mrr_aRef= np.full(self.ground_plat_dims, np.nan)
        self.mrr_ref = np.full(self.ground_plat_dims, np.nan)
        self.mrr_RR  = np.full(self.ground_plat_dims, np.nan)
        self.mrr_LWC = np.full(self.ground_plat_dims, np.nan)
        self.mrr_VEL = np.full(self.ground_plat_dims, np.nan)
        self.mrr_disdro_Dm = np.full(self.ground_plat_dims, np.nan)
        self.mrr_disdro_Nw = np.full(self.ground_plat_dims, np.nan)
        self.mrr_disdro_Z  = np.full(self.ground_plat_dims, np.nan)
        self.mrr_ref_ku_adj = np.full(self.ground_plat_dims, np.nan)
        self.mrr_data_quality = np.full(self.ground_plat_dims, np.nan)
        
        #define variable names and units
        self.mrr_PIA_name = 'two-way integrated attenuation'
        self.mrr_aRef_name = 'attenuated reflectivity'
        self.mrr_ref_name = 'reflectivity corrected for attenuation due to raindrops'
        self.mrr_refKu_name = 'reflectivity converted to ku using observed DSDs'
        self.mrr_RR_name = 'rainrate'
        self.mrr_LWC_name = 'liquid water content'
        self.mrr_WVEL_name = 'velocity: capital W'
        self.mrr_disdro_Dm_name = 'mass-weighted mean diameter computed from DSD profile'
        self.mrr_disdro_Nw_name = 'normalized intercept parameter computed from DSD profile'
        self.mrr_disdro_Z_name = 'radar reflectivity computed from DSD profile'
        self.mrr_data_quality_name = '% of valid spectra during the .ave interval'
        self.mrr_units_PIA = 'dB'
        self.mrr_units_ref = 'dBZ'
        self.mrr_units_RR = 'mm h^-1'
        self.mrr_units_LWC = 'g m^-3'
        self.mrr_units_WVEL = 'm s^-1'
        self.mrr_units_disdro_Dm = 'mm'
        self.mrr_units_disdro_Nw = 'log(Nw)'
        self.mrr_units_disdro_Z = 'dBZ'
        self.mrr_units_data_quality = '%'
        
        # Will need to loop through input list of MRR units:
        #for i, mrr_id in enumerate(mrr_IDs):
        for i, file_mrr in enumerate(self.mrr_files):
            #file_mrr = files_mrr[i]
            mrr_id = os.path.basename(file_mrr)[0:7]
        
            # Set initial mrr_plat_info values:
            #pdb.set_trace()
            mrr_plat_info  = self.set_plat_values(mrr_id.upper(), file_mrr)
            unit_elevation = mrr_plat_info['elev']
        
            #get subscript where MRR is closest to a column grid point:
            mrr_loc_decdeg = sim.dms2dd(mrr_plat_info['lat_d'], mrr_plat_info['lat_m'], mrr_plat_info['lat_s'], 
                                    mrr_plat_info['lon_d'], mrr_plat_info['lon_m'], mrr_plat_info['lon_s'])
            mrr_lat = mrr_loc_decdeg[0]
            mrr_lon = mrr_loc_decdeg[1]
            lat_ok = 'N'
            lon_ok = 'N'
            if mrr_lat <= self.box_max_lat and mrr_lat >= self.box_min_lat: lat_ok = 'Y'
            if mrr_lon <= self.box_max_lon and mrr_lon >= self.box_min_lon: lon_ok = 'Y'
            if lat_ok == 'N' and lon_ok == 'N':
                print('------------------------------------------------')
                print('-------- before reading .ave file --------------')
                print('---APPEARS THIS MRR IS NOT IN THE COLUMN BOX!---')
                print(f'-------plat_name: {mrr_plat_info["plat_name"]}')
                print('-------NOTICE: MRR data NOT SET IN column grid!')
                self.mrr_info = np.nan
                continue # so will only continue from here & read in mrr file if location is OK
        
            # Call read_data_mrr Function to read in MRR .ave data:
            data = mrr.read_data(file_mrr)
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
            for t, this_date in enumerate(self.interval_datetime):
                time_sub = np.where(mrr_datetime == this_date)[0]
                if len(time_sub) < 1:
                    continue # this time is not available in mrr data
                elif len(time_sub) > 1:
                    #found > ONE MRR time corresponding to this interval time:
                    #this very highly unlikely as the MRR .ave files (SIMBA input) typically use 60 s averaging interval - will set
                    #here to use 1st element that matches the interval time in case it occurs, but shouldn't really ever get here:
                    print(f'HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")}')
                    print(f' in ave file: {os.path.basename(file_mrr)}')
                    pdb.set_trace()

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
                    mrr_int_width_array.append(self.halftime_interval)
                    # just will set main_plat time w/ sec as '00' since using an interval of data
                    mrr_plat_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                    mrr_offset_array.append(0)

                # Find the (horizontal) column grid point closest to the MRR -> place MRR values at this point
                # -> arrays of the MRR data currently set up as 1-D arrays: [ht] in the column box
                # -> arrays that will set data into were defined above this loop & have same dims as column box grid + time [time interval, lon, lat, elev]
                # NOTE: will maintain original MRR location in plat_info for attributes
                closest_col_lat_sub = sim.closest(self.lat_values, mrr_lat)[0]
                closest_col_lon_sub = sim.closest(self.lon_values, mrr_lon)[0]
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
                disdro_result = mrr.get_dsd_params(N_of_D_array, D_bins_array, mrr_gate_heights, no_fit=True)
                disdro_Dm = disdro_result['dm']             #[mm]
                disdro_Nw = disdro_result['nw']             #just Nw (NOT log(Nw))
                disdro_Z  = disdro_result['reflectivity']   #[dBZ]
                #pdb.set_trace()
                
                # Interpolate in the vertical to match column grid levels:  for profile fields from .ave file
                # use scipy interpolate function to get val from MRR hts to column grid hts:
                col_PIA         = sim.interp(PIA, mrr_heights, self.z_values)
                col_rainRate    = sim.interp(rainRate, mrr_heights, self.z_values)
                col_LWC         = sim.interp(LWC, mrr_heights, self.z_values)
                col_velocity    = sim.interp(velocity, mrr_heights, self.z_values)
                #pdb.set_trace()
                # For reflectivity & Nw:  get INTERPOL results for LINEAR UNITS, then put into LOG UNITS after
                orig_lin_aRef    = 10.0**((attenRef)/10.0)
                orig_lin_ref     = 10.0**((ref)/10.0)
                orig_lin_disdro_Z= 10.0**((disdro_Z)/10.0)
                col_attenRef  = sim.interp(orig_lin_aRef, mrr_heights, self.z_values)
                col_ref       = sim.interp(orig_lin_ref, mrr_heights, self.z_values)
                col_disdro_Z  = sim.interp(orig_lin_disdro_Z, mrr_heights, self.z_values)
                col_disdro_Nw = sim.interp(disdro_Nw, mrr_heights, self.z_values)
                col_disdro_Dm = sim.interp(disdro_Dm, mrr_heights, self.z_values)
                #pdb.set_trace()
                
                # put reflectivities & Nw into log units:
                col_disdro_Z  = 10*np.log10(col_disdro_Z) #[dBZ]
                col_disdro_Nw = np.log10(col_disdro_Nw)   #[log(Nw)]
                col_ref       = 10*np.log10(col_ref)      #[dBZ]
                col_ref_ku_adj = 0.974*col_ref - 0.097    #[dBZ] - derived using Wallops 2DVD DSD and radar model following Liao et al. 2020 method
                col_attenRef  = 10*np.log10(col_attenRef) #[dBZ]
                #pdb.set_trace()
         
                # Ensure have NaNs where should be:
                low_hts  = np.where(self.z_values < min(mrr_heights))[0]
                high_hts = np.where(self.z_values > max(mrr_heights))[0]
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
                self.mrr_PIA[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_PIA
                self.mrr_aRef[:, closest_col_lat_sub,closest_col_lon_sub,t]        = col_attenRef
                self.mrr_ref[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_ref
                self.mrr_RR[:, closest_col_lat_sub,closest_col_lon_sub,t]          = col_rainRate
                self.mrr_LWC[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_LWC
                self.mrr_VEL[:, closest_col_lat_sub,closest_col_lon_sub,t]         = col_velocity
                self.mrr_disdro_Z[:, closest_col_lat_sub,closest_col_lon_sub,t]    = col_disdro_Z
                self.mrr_ref_ku_adj[:, closest_col_lat_sub,closest_col_lon_sub,t]  = col_ref_ku_adj
                self.mrr_disdro_Dm[:, closest_col_lat_sub,closest_col_lon_sub,t]   = col_disdro_Dm
                self.mrr_disdro_Nw[:, closest_col_lat_sub,closest_col_lon_sub,t]   = col_disdro_Nw
                self.mrr_data_quality[:, closest_col_lat_sub,closest_col_lon_sub, 0] = quality

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
           return np.nan
        
        
        #Set up dictionary
        self.mrr_info = {'lat_d':mrr_lat_d_array, 'lat_m':mrr_lat_m_array, 'lat_s':mrr_lat_s_array,
                           'lon_d':mrr_lon_d_array, 'lon_m':mrr_lon_m_array, 'lon_s':mrr_lon_s_array,
                           'elev':mrr_elev_array,
                           'plat_name':mrr_plat_name_array, 
                           'plat_type':mrr_plat_type_array,
                           'operation_mode':mrr_op_mode_array,
                           'wavelength':mrr_wavelength_array,
                           'frequency':mrr_frequency_array,
                           'beam_width':mrr_beam_width_array,
                           'gate_size':mrr_gate_size_array,
                           'timestamp':mrr_plat_timestamp_array,
                           'offset_vs_main':mrr_offset_array,
                           'time_interval_width':mrr_int_width_array }
        self.mrr_data_in_column = True

    # ***************************************************************************************
    def get_gauges_for_column(self):
        '''
        ; ------------------------------------------------------------
        ; Of available .gmin gauge data files, locate gauge data at
        ; main_plat time +/- halftime_interval & locations within column
        ; box to return for putting into column .nc file
        ;------------------------------------------------------------
        ;  -> Currently:  is set up to work with .gmin gauge files, as
        ;	generated by Jerry Wang (see Wang et al. 2008, JTECH).
        ;	These data available at:
        ;	gsfc hector machine: /gpmraid/gpmarchive/Gauge/
        ;	-or-  http://gpm-gv.gsfc.nasa.gov/Gauge/index.html
        ;
        ; Code below first includes a function to define template for
        ; reading in .gmin files. Second function is the main function
        ; that reads in the files, locates the appropriate data, and
        ; returns gauge_info_and_data structure to build_column.pro
        ;
        ;  result = get_gauges_for_column(column_box_params, gauge_dir,$
        ;				main_plat_info, halftime_interval)
        ;	column_box_params:   column box grid parameters. This is the
        ;				structure returned by define_box.pro
        ;	gauge_dir:	     string dir path to .gmin files with
        ;				times for the present case
        ;	main_plat_info:      main_plat_info structure returned by
        ;				set_main_plat_values.pro
        ;	halftime_interval:   integer # of mins before & after the main
        ;				platform time to include in column file
        ;
        ;	result:  if gauges exist in column box, returns a big structure:
        ;		gauges_info_and_data:	Structure made up of sturctures:
        ;		.gauges_plat_info:		arrays of info for ea. avail
        ;						  rain gauge in col box
        ;		.col_data_gauges:    	Structure will have 2 tags:
        ;		  .gauges_rainrate		3-D (column grid dims) array of
        ;					 	  values for gauge-based rain rate
        ;					 	  data, only populated at sfc [z dim: 0]
        ;		  .units_rainrate		[mm h^-1]
        ;		 if there are no gauges within the bounds of the column grid box,
        ;		 then the main function returns result = -1
        ;
        ;	Dependencies:
        ;	  dd2dms.pro	utility to convert decimal degrees to deg-min-sec
        ;	  dms2dd.pro	utility to convert deg-min-sec to decimal degrees
        ;	  closest.pro	utility to locate value in an array nearest intput value
        ;	  
        ;	Updates by Charanjit S. Pabla:
        ;	  Modified rainrate column char location element in use_for_locs array in template_gmin function - 08/15/2019
        ;	  Modified handing multiple values at a grid point -- previously took the average of previous and current value...
        ;	  now take the max between the previous and current value - 08/15/2019
        ;
        ;------------------------------------------------------------
        '''
    
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
        
        self.gauge_rainrate = np.full(self.ground_plat_dims, np.nan)
        self.gauge_units_rainrate = 'mm h^-1'
    
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
        for i, gmin_file in enumerate(self.gauge_gmin_files):
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
            
            for t, this_date in enumerate(self.interval_datetime):
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
                if gmin_lat[time_sub[0]] <= self.box_max_lat  and  gmin_lat[time_sub[0]] >= self.box_min_lat: lat_ok = 'Y'
                if gmin_lon[time_sub[0]] <= self.box_max_lon  and  gmin_lon[time_sub[0]] >= self.box_min_lon: lon_ok = 'Y'
                
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
                        gauge_int_width_array.append(self.halftime_interval)
                        # just will set main_plat time w/ sec as '00' since using an interval of data
                        gauge_plat_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                        gauge_offset_array.append(0)
                    print(f'         at timestamp: {this_date.strftime("%m/%d/%Y %H%M%S")}')
                    closest_col_lat_sub = sim.closest(self.lat_values, unit_lat)
                    closest_col_lon_sub = sim.closest(self.lon_values, unit_lon)
                    if len(closest_col_lat_sub) != 1 or len(closest_col_lon_sub) != 1:
                        print(f'---- PROBLEM SETTING WHERE IN GRID TO PLACE GAGUE: {os.path.basename(gmin_file)}')
                        pdb.set_trace()
                    else:
                        closest_col_lat_sub = closest_col_lat_sub[0]
                        closest_col_lon_sub = closest_col_lon_sub[0]
    
                    # Determine if the current value at the grid spot nearest this gauge for this time is already populated:
                    #  use a str/char variable called "exists" to tell if have a previous value:  Y or N only
                    if ~ np.isfinite(self.gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t]): exists = 'N'
                    if np.isfinite(self.gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t]):   exists = 'Y'
    
                    # 3a)  If no previous value, populate nearest grid point for this time w current gauge data
                    if exists == 'N':
                        self.gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t] = gmin_rain[time_sub[0]]
    
                    # 3b)  If there IS a previous value, take the max between the previous and current value
                    #        This will be the case for situations with Iowa-type/paired tip bucket gauges
                    #        or whenever multiple stand-alone gauges at same location (eg: WFF N-159 pad)
                    if exists == 'Y':
                        prev_val = self.gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t]  #[t] is sub into interval times for new data holder array
                        this_val = gmin_rain[time_sub[0]]           #[time_sub] is sub into the arrays read in from the gmin files, corresponds to the current interval_*[t]
                        new_val = max([prev_val, this_val]) #take the maximum rate instead of averaging -- updated by Cpabla 08/13/2019
                        #new_val  = mean([prev_val, this_val]) ;previously took the mean between previous and current value
                        self.gauge_rainrate[0,closest_col_lat_sub,closest_col_lon_sub,t] = new_val
    
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
            return -1
    
        #Set up dictionary
        self.gauges_info = {'lat_d':gauge_lat_d_array,'lat_m':gauge_lat_m_array,'lat_s':gauge_lat_s_array,
                              'lon_d':gauge_lon_d_array,'lon_m':gauge_lon_m_array,'lon_s':gauge_lon_s_array,
                              'plat_name':gauge_plat_name_array,
                              'plat_type':gauge_plat_type_array,
                              'operation_mode':gauge_file_type_array,
                              'timestamp':gauge_plat_timestamp_array,
                              'offset_vs_main':gauge_offset_array,
                              'time_interval_width':gauge_time_interval_width}
        self.gauges_data_in_column = True

    # ***************************************************************************************
    def get_2dvd_for_column(self):
        '''
        ;------------------------------------------------------------
        ; Of the available 2DVD units, identify one(s) in column box
        ; and return data at main_plat time +/- halftime_interval
        ;------------------------------------------------------------
        ;  -> 2DVD locations MUST BE CORRECT in platform_location.pro
        ;
        ;  -> two_dvd_dir: dir with the avail 2DVD files that include time
        ;		   frame for the present case (files processed by
        ;		   Ali Tokay, see Toaky et al. 2001 JAM article)
        ;
        ;  -> Currently: set up to work with these kind of files:
        ;	- *sn##_dropcount*
        ;	- 2dvd_sn##_raindsd.*   	OR   2dvd_sn##_raindsd_50*
        ;	- 2dvd_sn##_raindsd_ter*
        ;	- 2dvd_sn##_rainparameter.*	OR   2dvd_sn##_rainparameter_50*
        ;	- 2dvd_sn##_rainparameter_ter*
        ;	- 2dvd_sn##_multifreq_attenuation_ter*
        ;	- 2dvd_sn##_multifreq_reflectivity_ter*
        ;	also supports (from OLYMPEX campaign):
        ;	- 2dvd_sn##_raindsd_100%.*
        ;	- 2dvd_sn##_raindsd_100%_ter.*
        ;	- 2dvd_sn##_rainparameter_100%.*
        ;	- 2dvd_sn##_rainparameter_100%_ter.*
        ;
        ;	File sent to set_plat_values.pro:  2dvd_sn##_raindsd_ter.*
        ;			 	(or 2dvd_sn##_raindsd.* if no _ter. file)
        ;;
        ;  result = get_2dvd_for_column(column_box_params, two_dvd_dir, $
        ;				 main_plat_info, halftime_interval)
        ;	column_box_params:   column box grid parameters. This is the 
        ;				structure returned by define_box.pro
        ;	apu_dir:	     string dir path to 2DVD files with
        ;				times for the present case
        ;	main_plat_info:      main_plat_info structure returned by 
        ;				set_main_plat_values.pro
        ;	halftime_interval:   integer # of mins before & after the main
        ;				platform time to include in column file
        ;
        ;	result:  returns a big structure:
        ;		two_dvd_info_and_data:  Structure made up of sturctures:
        ;		.two_dvd_plat_info:		arrays of info for ea.
        ;						  avail 2DVD unit in col box
        ;		OTHER TAGS ARE STURCTURS (OR -1 IF FILE IS NOT AVAILABLE):
        ;		.col_data_2dvd_dropcounts:	fields from Ali's 
        ;						  *sn##_dropcount*
        ;		.col_data_2dvd_raindsd:		fields from Ali's
        ;						  2dvd_sn##_raindsd[_50].*
        ;		.col_data_2dvd_raindsd_100:	fields from Ali's
        ;						  2dvd_sn##_raindsd_100%.*
        ;		.col_data_2dvd_raindsd_ter:	fields from Ali's
        ;						  2dvd_sn##_raindsd_ter.*
        ;		.col_data_2dvd_raindsd_ter:	fields from Ali's
        ;						  2dvd_sn##_raindsd_100%_ter.*
        ;		.col_data_2dvd_rainparam:	fields from Ali's 
        ;						  2dvd_sn##_rainparam[_50].*
        ;		.col_data_2dvd_rainparam_100%:	fields from Ali's 
        ;						  2dvd_sn##_rainparam_100%.*
        ;		.col_data_2dvd_rainparam_ter:	fields from Ali's 
        ;						  2dvd_sn##_rainparam_ter.*
        ;		.col_data_2dvd_rainparam_100%_ter:  fields from Ali's 
        ;						  2dvd_sn##_rainparam_100%_ter.*
        ;		.col_data_2dvd_atten:		fields from Ali's
        ;						  2dvd_sn##_multifreq_attenuation_ter*
        ;		.col_data_2dvd_ref:		fields from Ali's 
        ;						  2dvd_sn##_multifreq_reflectivity_ter*
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
        for this_file in self.twodvd_files:
            file_name = os.path.basename(this_file)
            sn_unit = file_name[file_name.find('sn'):file_name.find('sn')+4]
            #file_name = theList[i]
            #sn_spot   = where(strmatch(file_name, 'sn*') eq 1)
            possible_2DVDs.append(sn_unit.upper())
        if len(possible_2DVDs) > 1: possible_2DVDs=np.unique(possible_2DVDs[1:])
        print(f'    possible units: {possible_2DVDs}')
        twodvd_dir = os.path.dirname(self.twodvd_files[0])

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
        drop_d = (50,)
        drop_dims = self.ground_plat_dims + drop_d 
        self.twodvd_dropcount = np.full(drop_dims, np.nan)

        # From _raindsd. OR _raindsd_50%. file:
        self.twodvd_dsd_con = np.full(drop_dims, np.nan)

        # From _raindsd_100%. file:
        self.twodvd_dsd_100_con = np.full(drop_dims, np.nan)

        # From _raindsd_ter. OR _raindsd_50%_ter. file:
        self.twodvd_dsd_ter_con = np.full(drop_dims, np.nan)

        # From _raindsd_100%_ter. file:
        self.twodvd_dsd_100_ter_con = np.full(drop_dims, np.nan)
 
        # From _rainparameter. OR _rainparameter_50%. file:
        #  THESE ARE EACH 4-D ARRAY [times in interval  X  column x dir  X  column y dir  X  column z dir  ]:
        #dims = (n_interval_times, len(self.column_box_params['column_grid_lons']),  len(self.column_box_params['column_grid_lats']),  len(self.column_box_params['column_grid_alts'])  )
        self.twodvd_param_tot_ndrops = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_tot_con    = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_LWC        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_rainrate   = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_refRayleigh= np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_Dm         = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_Dmax       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_Dmin       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_stdev_Dm   = np.full(self.ground_plat_dims, np.nan)

        # From _rainparameter_100%. file:
        self.twodvd_param_100_tot_ndrops = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_tot_con    = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_LWC        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_rainrate   = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_refRayleigh= np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_Dm         = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_Dmax       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_Dmin       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_stdev_Dm   = np.full(self.ground_plat_dims, np.nan)

        # From _rainparameter_ter. OR _rainparameter_50%_ter. file:
        self.twodvd_param_ter_tot_ndrops = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_tot_con    = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_LWC        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_rainrate   = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_refRayleigh= np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_Dm         = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_Dmax       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_Dmin       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_ter_stdev_Dm   = np.full(self.ground_plat_dims, np.nan)

        # From rainparameter_100%_ter. file:
        self.twodvd_param_100_ter_tot_ndrops = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_tot_con    = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_LWC        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_rainrate   = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_refRayleigh= np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_Dm         = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_Dmax       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_Dmin       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_param_100_ter_stdev_Dm   = np.full(self.ground_plat_dims, np.nan)

        # From _multifreq_attenuation_ter. file:
        self.twodvd_atten_rainrate = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_S	       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_C        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_X        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_Ku       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_K        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_Ka       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_atten_W        = np.full(self.ground_plat_dims, np.nan)

        # From multifreq_reflectivity_ter. file:
        self.twodvd_ref_rainrate = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_Rayleigh = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_S	     = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_C        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_X        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_Ku       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_K        = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_Ka       = np.full(self.ground_plat_dims, np.nan)
        self.twodvd_ref_W        = np.full(self.ground_plat_dims, np.nan)
        
        #units for all fields
        self.twodvd_units_ndrops = 'number of drops'
        self.twodvd_units_concen = 'drops m^-3 of air'
        self.twodvd_units_LWC = 'g m^-3'
        self.twodvd_units_rainrate = 'mm h^-1'
        self.twodvd_units_refRayleigh = 'dBZ'
        self.twodvd_units_Dm = 'mm'
        self.twodvd_units_Dmax = 'mm'
        self.twodvd_units_Dmin = 'mm'
        self.twodvd_units_stdevDm = 'mm'
        self.twodvd_units_atten = 'dB km^-1'
        self.twodvd_units_ref = 'dBZ'
    
        # Loop thru all possible 2DVD units:
        for u, unit_name in enumerate(possible_2DVDs):
            #unit_name       = possible_2DVDs[u]
            twodvd_plat_info = self.set_plat_values(unit_name, 'dummy.txt')
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
            if unit_lat <= self.box_max_lat and unit_lat >= self.box_min_lat: lat_ok = 'Y'
            if unit_lon <= self.box_max_lon and unit_lon >= self.box_min_lon: lon_ok = 'Y'
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
                rain_dsd         = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd.*')
                if len(rain_dsd) < 1: rain_dsd = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_50%.*')
                rain_dsd_ter     = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_ter.*')  
                if len(rain_dsd_ter) < 1: rain_dsd_ter = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_50%_ter.*')
                rain_dsd_100     = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_100%.*')
                rain_dsd_100_ter = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_raindsd_100%_ter.*')  
                rain_param       = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter.*')
                if len(rain_param) < 1: rain_param = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_50%.*')
                rain_param_ter   = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_ter*')
                if len(rain_param_ter) < 1: rain_param_ter = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_50%_ter.*')
                rain_param_100   = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_100%.*')
                rain_param_100_ter=glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_rainparameter_100%_ter.*')
                atten_ter         = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_multifreq_attenuation_ter*')
                ref_ter           = glob.glob(f'{twodvd_dir}/2dvd_sn{unit_id}_multifreq_reflectivity_ter*')
                #pdb.set_trace()
                
                # check if file exist or not and print out:
                self.twodvd_dsd_file = 0
                self.twodvd_dsd_ter_file = 0
                self.twodvd_dsd_100_file = 0
                self.twodvd_dsd_100_ter_file =0
                self.twodvd_param_file = 0
                self.twodvd_param_ter_file = 0
                self.twodvd_param_100_file = 0
                self.twodvd_param_100_ter_file=0
                self.twodvd_dropcount_file = 0
                self.twodvd_atten_file = 0
                self.twodvd_ref_file = 0
                print('-----------------------------')
                print(f'   avail files for 2DVD unit:  {unit_name}')
                if len(dropcounts)       == 1: 
                    self.twodvd_dropcount_file  = 1
                    print(f'     dropcount: {os.path.basename(dropcounts[0])}')
                else:
                    print('     dropcount: ')
                if len(rain_dsd)         == 1: 
                    self.twodvd_dsd_file = 1
                    print(f'       50% dsd: {os.path.basename(rain_dsd[0])}')
                else:
                    print('       50% dsd: ')
                if len(rain_dsd_ter)     == 1:
                    self.twodvd_dsd_ter_file = 1
                    print(f'   50% dsd ter: {os.path.basename(rain_dsd_ter[0])}')
                else:
                    print('   50% dsd ter: ')
                if len(rain_dsd_100)     == 1:
                    self.twodvd_dsd_100_file = 1
                    print(f'      100% dsd: {os.path.basename(rain_dsd_100[0])}')
                else:
                    print('      100% dsd: ')
                if len(rain_dsd_100_ter) == 1:
                    self.twodvd_dsd_100_ter_file = 1
                    print(f'  100% dsd ter: {os.path.basename(rain_dsd_100_ter[0])}')
                else:
                    print('  100% dsd ter: ')
                if len(rain_param)       == 1:
                    self.twodvd_param_file = 1
                    print(f'      50% parm: {os.path.basename(rain_param[0])}')
                else:
                    print('      50% parm: ')
                if len(rain_param_ter)   == 1:
                    self.twodvd_param_ter_file = 1
                    print(f'  50% parm ter: {os.path.basename(rain_param_ter[0])}')
                else:
                    print('  50% parm ter: ')
                if len(rain_param_100)   == 1:
                    self.twodvd_param_100_file = 1
                    print(f'    100% param: {os.path.basename(rain_param_100[0])}')
                else:
                    print('    100% param: ')
                if len(rain_param_100_ter) == 1:
                    self.twodvd_param_100_ter_file=1
                    print(f'100% param ter: {os.path.basename(rain_param_100_ter[0])}')
                else:
                    print('100% param ter: ')
                if len(atten_ter)        == 1:
                    self.twodvd_atten_file = 1
                    print(f'   attenuation: {os.path.basename(atten_ter[0])}')
                else:
                    print('   attenuation: ')
                if len(ref_ter)          == 1:
                    self.twodvd_ref_file = 1
                    print(f'  reflectivity: {os.path.basename(ref_ter[0])}')
                else:
                    print('  reflectivity: ')
                print(f'-------------')
                #pdb.set_trace()
                # If no supported files available for this 2DVD unit, then skip it:
                if  self.twodvd_dsd_file     == 1 or self.twodvd_dsd_ter_file   == 1 or self.twodvd_dsd_100_file   == 1 or self.twodvd_dsd_100_ter_file == 1 or \
                    self.twodvd_param_file   == 1 or self.twodvd_param_ter_file == 1 or self.twodvd_param_100_file == 1 or self.twodvd_param_100_ter_file == 1 or \
                    self.twodvd_dropcount_file == 1 or self.twodvd_atten_file == 1   or self.twodvd_ref_file == 1:
    
                    # CALL TO set_plat_values.pro:  FOR LOC ATTRIBUTES.
                    unit_in_column_test_cnt = 0   # counter to see if have added this unit to big plat info arrays

                    # find column grid point closest to the current 2DVD unit
                    #  - arrays of data from the 2DVD files have [field, time] dims
                    #   - Field #: which column in the data file to use  (see README_2DVD_data_formats.txt)
                    #   - time #:  use the subscript for the timestamp that is needed
                    #   - will use NumPy genfromtxt function to read ascii data (this function works well with missing data)
                    closest_col_lat_sub = sim.closest(self.lat_values, unit_lat)
                    closest_col_lon_sub = sim.closest(self.lon_values, unit_lon)
                    if len(closest_col_lat_sub) != 1 or len(closest_col_lon_sub) != 1:
                        print('-----PROBLEM SETTING WHERE IN GRID TO PLACE 2DVD VALUES!!--')
                        pdb.set_trace() ##this should not happen -- check and find a way tp skip this
                    else:
                        closest_col_lat_sub = closest_col_lat_sub[0]
                        closest_col_lon_sub = closest_col_lon_sub[0]
                    
                    # NOW - IF KNOW THE UNIT LOC IS IN THE COL BOX, OPEN EA TYPE OF 2DVD FILE,
                    #   SEE IF TIMES IN INTERVAL AVAILABLE, POPULATE ARRAYS W DATA FROM THIS UNIT
                    if self.twodvd_dropcount_file == 1:
                        dropcounts_data    = np.genfromtxt(dropcounts) #will return float 2d array
                        dc_year = dropcounts_data[:,0]
                        dc_Jday = dropcounts_data[:,1]
                        dc_hour = dropcounts_data[:,2]
                        dc_min  = dropcounts_data[:,3]
                        dc_datetime = sim.jul2datetime(dc_year, dc_Jday, dc_hour, dc_min)

                        for t, this_date in enumerate(self.interval_datetime):
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hour+self.main_min+'00')
                                    unit_offset_array.append(0)
                        
                                # populate big column grid w/ 2DVD data from current unit's *dropcount* file for THIS time:
                                # bin 0 starts at column 5 in *dropcount* file or python index 4
                                for binn in range(self.twodvd_dropcount.shape[0]):
                                    self.twodvd_dropcount[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = dropcounts_data[time_sub, binn+4]
                    #pdb.set_trace()
                    if self.twodvd_dsd_file == 1:
                        rain_dsd_data   = np.genfromtxt(rain_dsd[0])
                        #rain_dsd_data   = np.loadtxt(rain_dsd)
                        dsd_year        = rain_dsd_data[:,0]
                        dsd_Jday        = rain_dsd_data[:,1]
                        dsd_hour        = rain_dsd_data[:,2]
                        dsd_min         = rain_dsd_data[:,3]
                        dsd_datetime    = sim.jul2datetime(dsd_year, dsd_Jday, dsd_hour, dsd_min)

                        for t, this_date in enumerate(self.interval_datetime):
                            time_sub  = np.where(dsd_datetime == this_date)[0]
                            if len(time_sub) < 1:
                                continue # this time not avail in the 2dvd file, skip this interval time
                            elif len(time_sub) > 1:
                                # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                                print(f'  --- in 2DVD file: {os.path.basename(rain_dsd[0])}')
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
                              
                                # populate big column grid w/ 2DVD data from current unit's *raindsd.*/*raindsd_50* file for THIS time:
                                # bin 0 starts at column 5 in *raindsd* file or python index 4
                                #pdb.set_trace()
                                for binn in range(self.twodvd_dsd_con.shape[0]):
                                    self.twodvd_dsd_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_data[time_sub, binn+4]

                    if self.twodvd_dsd_100_file == 1:
                        rain_dsd_100_data   = np.genfromtxt(rain_dsd_100[0])
                        dsd_100_year        = rain_dsd_100_data[:,0]
                        dsd_100_Jday        = rain_dsd_100_data[:,1]
                        dsd_100_hour        = rain_dsd_100_data[:,2]
                        dsd_100_min         = rain_dsd_100_data[:,3]
                        dsd_100_datetime    = sim.jul2datetime(dsd_100_year, dsd_100_Jday, dsd_100_hour, dsd_100_min)

    
                        for t, this_date in enumerate(self.interval_datetime):
                            time_sub  = np.where(dsd_100_datetime == this_date)[0]
                            if len(time_sub) < 1:
                                continue # this time not avail in the 2dvd file, skip this interval time
                            elif len(time_sub) > 1:
                                # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                                print(f'  --- in 2DVD file: {os.path.basename(rain_dsd_100[0])}')
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
                            
                                # populate big column grid w/ 2DVD data from current unit's *raindsd_100%.* file for THIS time:
                                # bin 0 starts at column 5 in *raindsd* file or python index 4
                                for binn in range(self.twodvd_dsd_100_con.shape[0]):
                                    self.twodvd_dsd_100_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_100_data[time_sub, binn+4]
    
                    if self.twodvd_dsd_ter_file == 1:
                        rain_dsd_ter_data   = np.genfromtxt(rain_dsd_ter[0])
                        dsd_ter_year        = rain_dsd_ter_data[:,0]
                        dsd_ter_Jday        = rain_dsd_ter_data[:,1]
                        dsd_ter_hour        = rain_dsd_ter_data[:,2]
                        dsd_ter_min         = rain_dsd_ter_data[:,3]
                        dsd_ter_datetime    = sim.jul2datetime(dsd_ter_year, dsd_ter_Jday, dsd_ter_hour, dsd_ter_min)
    
                        for t, this_date in enumerate(self.interval_datetime):
                            time_sub  = np.where(dsd_ter_datetime == this_date)[0]
                            if len(time_sub) < 1:
                                continue # this time not avail in the 2dvd file, skip this interval time
                            elif len(time_sub) > 1:
                                # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                                print(f'  --- in 2DVD file: {os.path.basename(rain_dsd_ter[0])}')
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)

                                # populate big column grid w/ 2DVD data from current unit's *raindsd_ter.* file for THIS time:
                                # bin 0 starts at column 5 in *raindsd* file or python index 4
                                dsd_ter_main_time_sub = time_sub
                                for binn in range(self.twodvd_dsd_ter_con.shape[0]):
                                    self.twodvd_dsd_ter_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_ter_data[dsd_ter_main_time_sub, binn+4]

                    if self.twodvd_dsd_100_ter_file == 1:
                        rain_dsd_100_ter_data      = np.genfromtxt(rain_dsd_100_ter[0])
                        dsd_100_ter_year           = rain_dsd_100_ter_data[:,0]
                        dsd_100_ter_Jday           = rain_dsd_100_ter_data[:,1]
                        dsd_100_ter_hour           = rain_dsd_100_ter_data[:,2]
                        dsd_100_ter_min            = rain_dsd_100_ter_data[:,3]
                        dsd_100_ter_datetime       = sim.jul2datetime(dsd_100_ter_year, dsd_100_ter_Jday, dsd_100_ter_hour, dsd_100_ter_min)

                        for t, this_date in enumerate(self.interval_datetime):
                            time_sub  = np.where(dsd_100_ter_datetime == this_date)[0]
                            if len(time_sub) < 1:
                                continue # this time not avail in the 2dvd file, skip this interval time
                            elif len(time_sub) > 1:
                                # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                                print(f'  --- in 2DVD file: {os.path.basename(rain_dsd_100_ter[0])}')
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)

                                # populate big column grid w/ 2DVD data from current unit's *raindsd_100%_ter.* file for THIS time:
                                # bin 0 starts at column 5 in *raindsd* file or python index 4
                                dsd_100_ter_main_time_sub = time_sub
                                for binn in range(self.twodvd_dsd_100_ter_con.shape[0]):
                                    self.twodvd_dsd_100_ter_con[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = rain_dsd_100_ter_data[dsd_100_ter_main_time_sub, binn+4]

                    if self.twodvd_param_file == 1:
                        rain_param_data = np.genfromtxt(rain_param[0])
                        param_year      = rain_param_data[:,0]
                        param_Jday      = rain_param_data[:,1]
                        param_hour      = rain_param_data[:,2]
                        param_min       = rain_param_data[:,3]
                        param_datetime  = sim.jul2datetime(param_year, param_Jday, param_hour, param_min)

                        for t, this_date in enumerate(self.interval_datetime):
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
    
                                # populate column grid w/ 2DVD data from current unit's *parameter.*/*parameter_50* file:
                                param_main_time_sub = time_sub
                                self.twodvd_param_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_data[param_main_time_sub,4]
                                self.twodvd_param_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_data[param_main_time_sub,5]
                                self.twodvd_param_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_data[param_main_time_sub,6]
                                self.twodvd_param_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[param_main_time_sub,7]
                                self.twodvd_param_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_data[param_main_time_sub,8]
                                self.twodvd_param_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_data[param_main_time_sub,9]
                                self.twodvd_param_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[param_main_time_sub,10]
                                self.twodvd_param_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[param_main_time_sub,11]
                                self.twodvd_param_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[param_main_time_sub,12]
    
                    if self.twodvd_param_100_file == 1:
                        rain_param_100_data = np.genfromtxt(rain_param_100[0])
                        param_100_year      = rain_param_100_data[:,0]
                        param_100_Jday      = rain_param_100_data[:,1]
                        param_100_hour      = rain_param_100_data[:,2]
                        param_100_min       = rain_param_100_data[:,3]
                        param_100_datetime = sim.jul2datetime(param_100_year, param_100_Jday, param_100_hour, param_100_min)
    
                        for t, this_date in enumerate(self.interval_datetime):
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(main_year+main_month+main_day+'_'+main_hour+main_min+'00')
                                    unit_offset_array.append(0)

                                # populate column grid w/ 2DVD data from current unit's *parameter_100%.* file:
                                self.twodvd_param_100_main_time_sub = time_sub
                                self.twodvd_param_100_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_100_data[param_100_main_time_sub,4]
                                self.twodvd_param_100_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_100_data[param_100_main_time_sub,5]
                                self.twodvd_param_100_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_100_data[param_100_main_time_sub,6]
                                self.twodvd_param_100_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_100_data[param_100_main_time_sub,7]
                                self.twodvd_param_100_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_100_data[param_100_main_time_sub,8]
                                self.twodvd_param_100_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_100_data[param_100_main_time_sub,9]
                                self.twodvd_param_100_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_100_data[param_100_main_time_sub,10]
                                self.twodvd_param_100_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_100_data[param_100_main_time_sub,11]
                                self.twodvd_param_100_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_100_data[param_100_main_time_sub,12]
    
                    if self.twodvd_param_ter_file == 1:
                        rain_param_ter_data = np.genfromtxt(rain_param_ter[0])
                        param_ter_year      = rain_param_ter_data[:,0]
                        param_ter_Jday      = rain_param_ter_data[:,1]
                        param_ter_hour      = rain_param_ter_data[:,2]
                        param_ter_min       = rain_param_ter_data[:,3]
                        param_ter_datetime = sim.jul2datetime(param_ter_year, param_ter_Jday, param_ter_hour, param_ter_min)
    
                        for t, this_date in enumerate(self.interval_datetime):
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
    
                                # populate column grid w/ 2DVD data from current unit's *parameter_ter.* file:
                                self.twodvd_param_ter_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_data[time_sub,4]
                                self.twodvd_param_ter_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_data[time_sub,5]
                                self.twodvd_param_ter_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_data[time_sub,6]
                                self.twodvd_param_ter_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,7]
                                self.twodvd_param_ter_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_data[time_sub,8]
                                self.twodvd_param_ter_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_data[time_sub,9]
                                self.twodvd_param_ter_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,10]
                                self.twodvd_param_ter_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,11]
                                self.twodvd_param_ter_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,12]
    
                    if self.twodvd_param_100_ter_file == 1:
                        rain_param_100_ter_data = np.genfromtxt(rain_param_100_ter[0])
                        param_100_ter_year      = rain_param_100_ter_data[:,0]
                        param_100_ter_Jday      = rain_param_100_ter_data[:,1]
                        param_100_ter_hour      = rain_param_100_ter_data[:,2]
                        param_100_ter_min       = rain_param_100_ter_data[:,3]
                        param_100_ter_datetime = sim.jul2datetime(param_ter_year, param_100_ter_Jday, param_100_ter_hour, param_100_ter_min)
    
                        for t, this_date in enumerate(self.interval_datetime):
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
    
                                # populate column grid w/ 2DVD data from current unit's *parameter_100%_ter.* file:
                                self.twodvd_param_100_ter_tot_ndrops[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_param_data[time_sub,4]
                                self.twodvd_param_100_ter_tot_con[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_param_data[time_sub,5]
                                self.twodvd_param_100_ter_LWC[0, closest_col_lat_sub, closest_col_lon_sub, t]        = rain_param_data[time_sub,6]
                                self.twodvd_param_100_ter_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,7]
                                self.twodvd_param_100_ter_refRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]= rain_param_data[time_sub,8]
                                self.twodvd_param_100_ter_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]         = rain_param_data[time_sub,9]
                                self.twodvd_param_100_ter_Dmax[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,10]
                                self.twodvd_param_100_ter_Dmin[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_param_data[time_sub,11]
                                self.twodvd_param_100_ter_stdev_Dm[0, closest_col_lat_sub, closest_col_lon_sub, t]   = rain_param_data[time_sub,12]

                    if self.twodvd_atten_file == 1:
                        atten_ter_data  = np.genfromtxt(atten_ter[0])
                        atten_year      = atten_ter_data[:,0]
                        atten_Jday      = atten_ter_data[:,1]
                        atten_hour      = atten_ter_data[:,2]
                        atten_min       = atten_ter_data[:,3]
                        atten_datetime = sim.jul2datetime(atten_year, atten_Jday, atten_hour, atten_min)

                        for t, this_date in enumerate(self.interval_datetime):
                            time_sub  = np.where(atten_datetime == this_date)[0]
                            if len(time_sub) < 1:
                                continue #this time not avail in the 2dvd file, skip this interval time
                            elif len(time_sub) > 1:
                                # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                                print(f'  --- in 2DVD file: {os.path.basename(atten_ter[0])}')
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
    
                                # populate column grid w/ 2DVD data from current unit's *multifreq_attenuation_ter* file:
                                self.twodvd_atten_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, 0] = atten_ter_data[time_sub,4]
                                self.twodvd_atten_S[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,5]
                                self.twodvd_atten_C[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,6]
                                self.twodvd_atten_X[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,7]
                                self.twodvd_atten_Ku[0, closest_col_lat_sub, closest_col_lon_sub, t]       = atten_ter_data[time_sub,8]
                                self.twodvd_atten_K[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,9]
                                self.twodvd_atten_Ka[0, closest_col_lat_sub, closest_col_lon_sub, t]       = atten_ter_data[time_sub,10]
                                self.twodvd_atten_W[0, closest_col_lat_sub, closest_col_lon_sub, t]        = atten_ter_data[time_sub,11]
    
                    if self.twodvd_ref_file  ==  1:
                        ref_ter_data    = np.genfromtxt(ref_ter[0])
                        ref_year        = ref_ter_data[:,0]
                        ref_Jday        = ref_ter_data[:,1]
                        ref_hour        = ref_ter_data[:,2]
                        ref_min         = ref_ter_data[:,3]
                        ref_datetime    = sim.jul2datetime(ref_year, ref_Jday, ref_hour, ref_min)

                        for t, this_date in enumerate(self.interval_datetime):
                            time_sub  = np.where(ref_datetime == this_date)[0]
                            if len(time_sub) < 1:
                                continue # this time not avail in the 2dvd file, skip this interval time
                            elif len(time_sub) > 1:
                                # this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                                print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date.strftime("%H")} this_date.strftime("%M")')
                                print(f'  --- in 2DVD file: {os.path.basename(ref_ter[0])}')
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
                                    unit_int_width_array.append(self.halftime_interval)
                                    # just will set main_plat time w/ sec as '00' since using an interval of data
                                    unit_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                                    unit_offset_array.append(0)
    
                                # populate column grid w/ 2DVD data from current unit's *multifreq_reflectivity_ter* file:
                                self.twodvd_ref_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t] = ref_ter_data[time_sub,4]
                                self.twodvd_ref_Rayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t] = ref_ter_data[time_sub,5]
                                self.twodvd_ref_S[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,6]
                                self.twodvd_ref_C[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,7]
                                self.twodvd_ref_X[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,8]
                                self.twodvd_ref_Ku[0, closest_col_lat_sub, closest_col_lon_sub, t]       = ref_ter_data[time_sub,9]
                                self.twodvd_ref_K[0, closest_col_lat_sub, closest_col_lon_sub, t]        = ref_ter_data[time_sub,10]
                                self.twodvd_ref_Ka[0, closest_col_lat_sub, closest_col_lon_sub, t]       = ref_ter_data[time_sub,11]
                                self.twodvd_ref_W[0, closest_col_lat_sub, closest_col_lon_sub, t]           = ref_ter_data[time_sub,12]
    
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
          return -1
        print(f'------ 2DVDs in the column box: {unit_name_array}')
    
        # Set up dictionary of plat_infos for the 2DVDs in the column box:
        self.twodvd_info = {'lat_d':unit_lat_d_array, 'lat_m':unit_lat_m_array, 'lat_s':unit_lat_s_array,
                               'lon_d':unit_lon_d_array, 'lon_m':unit_lon_m_array, 'lon_s':unit_lon_s_array,
                               'plat_name':unit_name_array,
                               'plat_type':unit_type_array,
                               'operation_mode':unit_file_type_array,
                               'timestamp':unit_timestamp_array,
                               'offset_vs_main':unit_offset_array,
                               'time_interval_width':unit_int_width_array}

        self.twodvd_data_in_column = True

    # ***************************************************************************************
    def get_apu_for_column(self):
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
        rain_param_files     = sorted(glob.glob(self.apu_files+'/*_rainparameter_min.*'))
        rain_param_ter_files = sorted(glob.glob(self.apu_files+'/*_rainparameter_min_ter.*'))
        rain_dsd_files       = sorted(glob.glob(self.apu_files+'/*_raindsd_min.*'))
        rain_dsd_ter_files   = sorted(glob.glob(self.apu_files+'/*_raindsd_min_ter.*'))
        apu_ref_files        = sorted(glob.glob(self.apu_files+'/*_reflectivity_ter.*'))
        apu_atten_files      = sorted(glob.glob(self.apu_files+'/*_attenuation_ter.*'))
        n_rain               = len(rain_param_files)
        n_rain_ter           = len(rain_param_ter_files)
        n_dsd                = len(rain_dsd_files)
        n_dsd_ter            = len(rain_dsd_ter_files)
        n_ref                = len(apu_ref_files)
        n_atten              = len(apu_atten_files)
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
            if n_ref <= 1:
                if not apu_ref_files: n_ref=0
            if n_atten <= 1:
                if not apu_atten_files: n_atten=0

        # for files that are available, test to be sure have same # of files per ea APU
        list_num_of_files = np.array([n_rain, n_rain_ter, n_dsd, n_dsd_ter, n_ref, n_atten])
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
            pdb.set_trace()
            return
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
        dims_drop_bins = self.ground_plat_dims + drop_d
        if n_rain >= 1: # variables from file:  _rainparameter_min
            self.apu_n_drops            = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_concentration      = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_liqwater           = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_rain_rainrate      = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_rain_refinRayleigh = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_massweight_diam    = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_maximum_diam       = np.full(self.ground_plat_dims, np.nan, dtype=object)
        if n_rain_ter >= 1: # variables from file:  _rainparameter_min_ter
            self.apu_n_drops_ter            = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_concentration_ter      = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_liqwater_ter           = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_rain_rainrate_ter      = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_rain_refinRayleigh_ter = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_massweight_diam_ter    = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_maximum_diam_ter       = np.full(self.ground_plat_dims, np.nan, dtype=object)
        if n_dsd >= 1 or n_dsd_ter >= 1:
            self.apu_dsd_bins_Dmiddle = np.zeros(32)  # hold the bin middle diameters [mm]
            self.apu_dsd_bins_widths  = np.zeros(32)  # hold the bin widths [mm]
        if n_dsd >= 1: # variables from file:  _raindsd_min
            self.apu_drop_concen = np.full(dims_drop_bins, np.nan, dtype=object)
        if n_dsd_ter >= 1:  # variables from file:  _raindsd_min_ter
            self.apu_drop_concen_ter = np.full(dims_drop_bins, np.nan, dtype=object)
        if n_ref >= 1: # variables from file:  _reflectivity_ter
            self.apu_ref_rainrate       = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_refinRayleigh  = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atS            = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atC            = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atX            = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atKu           = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atK            = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atKa           = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_ref_atW            = np.full(self.ground_plat_dims, np.nan, dtype=object)
        if n_atten >= 1: # variables from file:  _attenuation_ter
            self.apu_atten_rainrate   = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atS        = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atC        = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atX        = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atKu       = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atK        = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atKa       = np.full(self.ground_plat_dims, np.nan, dtype=object)
            self.apu_atten_atW        = np.full(self.ground_plat_dims, np.nan, dtype=object)
    
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
                self.rain_ter_file = rain_param_ter_files[i]
                plat_name = os.path.basename(self.rain_ter_file)[0:5]  #get APU name as 'apu' + 2char string
                apu_plat_info = self.set_plat_values(plat_name, self.rain_ter_file)
            
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
                if apu_lat <= self.box_max_lat and apu_lat >= self.box_min_lat: lat_ok = 'Y'
                if apu_lon <= self.box_max_lon and apu_lon >= self.box_min_lon: lon_ok = 'Y'
                #print, '        apu:  ',plat_name
                #print, ' latlon OK?: ',lat_ok,lon_ok
                if lat_ok == 'Y' and lon_ok == 'Y':
                    print(f'    APU location in box: {apu_plat_info["plat_name"]}')

                    #APU Unit is in column box, so read in available files:
                    if n_rain_ter > 0:
                        self.apu_rain_ter_file = rain_param_ter_files[i]
                        rain_ter_data_full = np.loadtxt(self.apu_rain_ter_file)
                    else:
                        self.apu_rain_ter_file = None
                    if n_rain > 0:
                        self.apu_rain_file = rain_param_files[i]
                        rain_data_full = np.loadtxt(self.apu_rain_file)
                    else:
                        self.apu_rain_file = None
                    if n_dsd > 0:
                        self.apu_dsd_file = rain_dsd_files[i]
                        dsd_data_full  = np.loadtxt(self.apu_dsd_file)
                    else:
                        self.apu_dsd_file = None
                    if n_dsd_ter > 0:
                        self.apu_dsd_ter_file = rain_dsd_ter_files[i]
                        dsd_ter_data_full = np.loadtxt(self.apu_dsd_ter_file)
                    else:
                        self.apu_dsd_ter_file = None
                    if n_ref > 0:
                        self.apu_ref_file = apu_ref_files[i]
                        ref_data_full  = np.loadtxt(self.apu_ref_file)
                    else:
                        self.apu_ref_file = None
                    if n_atten > 0:
                        self.apu_atten_file = apu_atten_files[i]
                        atten_data_full= np.loadtxt(self.apu_atten_file)
                    else:
                        self.apu_atten_file = None

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
                    for t, this_date in enumerate(self.interval_datetime):
                        time_sub = np.where(apu_datetime == this_date)[0]

                        if len(time_sub) < 1: 
                            continue  # this time not avail in apu files
                        elif len(time_sub) > 1:
                            #this shouldn't happen (Ali's files shouldn't have repeated times) but in case..
                            print(f'----- HAVE MORE THAN 1 TIME MATCH FOR {this_date}')
                            print(f'  --- in APU file: {os.path.basename(self.rain_ter_file)}')
                            pdb.set_trace()
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
                            apu_int_width_array.append(self.halftime_interval)
                            # just will set main_plat time w/ sec as '00' since using an interval of data
                            apu_plat_timestamp_array.append(self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+'00')
                            apu_offset_array.append(0)
                        print(f'     at timestamp: {this_date.strftime("%m/%d/%Y %H%M%S")}')

                        # Find column grid point closest to the APU -> place the APU values at this point
                        #  -> arrays of data from the APU files have only [FIELD, TIME] dimensions
                        #	   --> FIELD number:  which column in the text file to use (see APUs_FilesNotes.txt)
                        #	   --> TIME number:  use the subscript for the main_plat_timestamp
                        #  -> arrays that am putting the values from the APU files into are set with same
                        #		dims as the column box grid [lon, lat, alt] --> now as [time, lon, lat, alt]
                        #  NOTE: original/exact APU locations are maintained in apu_plat_info for attributes
                        closest_col_lat_sub = sim.closest(self.lat_values, apu_lat)
                        closest_col_lon_sub = sim.closest(self.lon_values, apu_lon)
                        if len(closest_col_lat_sub) != 1 or len(closest_col_lon_sub) != 1:
                            print('-----PROBLEM SETTING WHERE IN GRID TO PLACE APU VALUES!!--')
                            pdb.set_trace()
                        else:
                            closest_col_lat_sub = closest_col_lat_sub[0]
                            closest_col_lon_sub = closest_col_lon_sub[0]
                        #pdb.set_trace()
                        # set values for APU fields based on APU file type:
                        
                        #the units are common among the files so will define only once here
                        self.apu_units_n_drops = 'number of drops'
                        self.apu_units_concentration = 'drops m^-3 of air'
                        self.apu_units_liqwatercontent = 'g m^-3'
                        self.apu_units_rainrate = 'mm h^-1'
                        self.apu_units_ref_inRayleigh = 'dBZ'
                        self.apu_units_massweight_diam = 'mm'
                        self.apu_units_maximum_diam = 'mm'
                        self.apu_units_bin_drop_concen = 'drops m^-3 mm^-1'
                        self.apu_units_all_refs = 'dBZ'
                        self.apu_units_all_attens = 'dB km^-1'
                        self.units_binDiameter_andWidth = 'mm'
                        #print(self.apu_rain_file)
                        #print(self.apu_rain_ter_file)
                        if self.apu_rain_file:
                            self.apu_n_drops[0, closest_col_lat_sub, closest_col_lon_sub, t]            = rain_data_full[time_sub,4]
                            self.apu_concentration[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_data_full[time_sub,5]
                            self.apu_liqwater[0, closest_col_lat_sub, closest_col_lon_sub, t]           = rain_data_full[time_sub,6]
                            self.apu_rain_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_data_full[time_sub,7]
                            self.apu_rain_refinRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_data_full[time_sub,8]
                            self.apu_massweight_diam[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_data_full[time_sub,9]
                            self.apu_maximum_diam[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_data_full[time_sub,10]
                        if self.apu_rain_ter_file:
                            #pdb.set_trace()
                            self.apu_n_drops_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]            = rain_ter_data_full[time_sub,4]
                            self.apu_concentration_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_ter_data_full[time_sub,5]
                            self.apu_liqwater_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]           = rain_ter_data_full[time_sub,6]
                            self.apu_rain_rainrate_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]      = rain_ter_data_full[time_sub,7]
                            self.apu_rain_refinRayleigh_ter[0, closest_col_lat_sub, closest_col_lon_sub, t] = rain_ter_data_full[time_sub,8]
                            self.apu_massweight_diam_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]    = rain_ter_data_full[time_sub,9]
                            self.apu_maximum_diam_ter[0, closest_col_lat_sub, closest_col_lon_sub, t]       = rain_ter_data_full[time_sub,10]
                        if self.apu_dsd_file:
                            for binn in range(self.apu_drop_concen.shape[0]):
                                self.apu_drop_concen[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = dsd_data_full[time_sub, binn+4]
                        if self.apu_dsd_ter_file:
                            for binn in range(self.apu_drop_concen_ter.shape[0]):
                                self.apu_drop_concen_ter[0, closest_col_lat_sub, closest_col_lon_sub, t, binn] = dsd_ter_data_full[time_sub, binn+4]
                        if self.apu_ref_file:
                            self.apu_ref_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]       = ref_data_full[time_sub,4]
                            self.apu_ref_refinRayleigh[0, closest_col_lat_sub, closest_col_lon_sub, t]  = ref_data_full[time_sub,5]
                            self.apu_ref_atS[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,6]
                            self.apu_ref_atC[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,7]
                            self.apu_ref_atX[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,8]
                            self.apu_ref_atKu[0, closest_col_lat_sub, closest_col_lon_sub, t]           = ref_data_full[time_sub,9]
                            self.apu_ref_atK[0, closest_col_lat_sub, closest_col_lon_sub, t]             = ref_data_full[time_sub,10]
                            self.apu_ref_atKa[0, closest_col_lat_sub, closest_col_lon_sub, t]           = ref_data_full[time_sub,11]
                            self.apu_ref_atW[0, closest_col_lat_sub, closest_col_lon_sub, t]            = ref_data_full[time_sub,12]
                        if self.apu_atten_file:
                            self.apu_atten_rainrate[0, closest_col_lat_sub, closest_col_lon_sub, t]     = atten_data_full[time_sub,4]
                            self.apu_atten_atS[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,5]
                            self.apu_atten_atC[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,6]
                            self.apu_atten_atX[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,7]
                            self.apu_atten_atKu[0, closest_col_lat_sub, closest_col_lon_sub, t]         = atten_data_full[time_sub,8]
                            self.apu_atten_atK[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,9]
                            self.apu_atten_atKa[0, closest_col_lat_sub, closest_col_lon_sub, t]         = atten_data_full[time_sub,10]
                            self.apu_atten_atW[0, closest_col_lat_sub, closest_col_lon_sub, t]          = atten_data_full[time_sub,11]
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
            return -1
        print(f'------APUs in the column box: {apu_plat_name_array}')
    
        # Populate arrays for DSD bin middle/width diameters, if have DSD files:
        if n_dsd >= 1 or n_dsd_ter >= 1:
            # there should be a parsivel_diameter.txt file in the apu dir if get here:
            diameters_file = self.apu_dir+'/parsivel_diameter.txt'
            diameters_data     = np.loadtxt(diameters_file)
            bin_middles      = diameters_data[:,0]
            bin_widths       = diameters_data[:,1]
            self.apu_dsd_bins_Dmiddle = bin_middles.flatten()
            self.apu_dsd_bins_widths  = bin_widths.flatten()
    
        # Set up dictionary of the plat_infos for APUs within the column box:
        #  need this to be the "plat_info" values for all APUs within box
        self.apu_info = {'lat_d':apu_lat_d_array,'lat_m':apu_lat_m_array,'lat_s':apu_lat_s_array,
                           'lon_d':apu_lon_d_array,'lon_m':apu_lon_m_array,'lon_s':apu_lon_s_array,
                           'plat_name':apu_plat_name_array,
                           'plat_type':apu_plat_type_array,
                           'operation_mode':apu_file_type_array,
                           'timestamp':apu_plat_timestamp_array,
                           'offset_vs_main':apu_offset_array,
                           'time_interval_width':apu_int_width_array}
        self.apu_data_in_column = True

    # ***************************************************************************************
    def grid_npol_for_column(self):
        '''
        ;------------------------------------------------------------
        ; grid NPOL radar data to a 200 x 200 x 15 km grid centered on the column box center
        ;  save full gridded file, then return subset within the column box
        ;------------------------------------------------------------
        ;   -> save the full gridded file to the input full_grid_dir
        ;
        ;   -> take subset of the gridded data that is within the column box and
        ;	send this subset back to the top-level build_column.pro for inclusion
        ;	in the new column .nc file
        ;
        ;  result = grid_npol_for_column(npol_plat_info, orig_NPOL_file, $
        ;				params_file, column_box_params, $
        ;				full_grid_dir)
        ;	npol_plat_info:	 plat_info structure returned by set_plat_values.pro
        ;	orig_NPOL_file:	 original NPOL data file, *.uf.gz
        ;	params_file:	 full path to Radx .params file to use for gridding
        ;	column_box_params:  column box grid parameters.  This needs to be the 
        ;				full structure returned by define_box.pro
        ;	full_grid_dir:   path to dir for full grid Radx output
        ;
        ;	result:	  structure containing the NPOL data at each location in
        ;		  in the column box grid (whether column grid is along or
        ;		  in middle of Radx/radar grid in the horizontal)
        ;	  Tags:
        ;		.**name:  string name of ** data field, as recorded in the Radx
        ;			   gridded .nc file (may be an empty string)
        ;		.**units: string telling the ** data field units, as recorded in
        ;		     	   the Radx gridded .nc file (may be an empty string)
        ;		.**data:  NPOL ** field value at locations in column box grid
        ;		    this array is 3-D: [X x Y x Z] dims, will match col box dims
        ;		    middle point of column box = middle point of this array
        ;		    IF EVEN # OF HORIZ COLUMN GRID SUBBOXES: 
        ;		    	middle point is along an edge
        ;		    IF ODD # OF HORIZ COLUMB GRID SUBBOXES:
        ;		        middle poit is at middle of a grid box
        ;
        ;     Have 3 Tags for each NPOL field getting put into the column:
        ;	
        ;    .ZZname    .CZname    .DRname    .RHname    .PHname    .KDname    .SQname
        ;    .ZZunits   .CZunits   .DRunits   .RHunits   .PHunits   .KDunits   .SQunits
        ;    .ZZdata    .CZdata    .DRdata    .RHdata    .PHdata    .KDdata    .SQdata
        ;
        ;    .SWname	.VRname    .RRname    .RPname	 .RCname    .D0name    .NWname
        ;    .SWunits	.VRunits   .RRunits   .RPunits	 .RCunits   .D0units   .NWunits
        ;    .SWdata	.VRdata    .RRdata    .RPdata	 .RCdata    .D0data    .NWdata
        ;
        ;    .FHname    .N2name	   .DMname
        ;    .FHunits   .N2units   .DMunits
        ;    .FHdata    .N2data    .DMdata
        ;
        ;	Dependencies:
        ;	  Radx software		UCAR/NCAR software for working with radar data
        ;	  uncomp_file.pro	utility to uncompress a .gz file
        ;	  bad2nan.pro		utility to convert a flag value to NaNs
        ;	  
        ;	  Updates:
        ;	    added code to convert from .cf --> .uf due to SIMBA/Radx gridding failing
        ;	    after 07/30/18, see lines 85-101 - C. Pabla (02/06/19)
        ;
        ;------------------------------------------------------------
        '''
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
        grid_center_lat_lon = self.cntr_lat_deg, self.cntr_lon_deg
        print('...running Py-ART gridding...')
        #Beg_timer = timer()
        grid = self.get_grid_from_radar(self.npol_file, grid_center_lat_lon)
        #End_timer = timer()
        #print(f”{Process took (End_timer-Beg_timer):0.2f} seconds”)
        #print(f'Process took {End_timer-Beg_timer:0.2f} seconds')
        
        #remove unzipped file if needed
        file_basename = os.path.basename(self.npol_file)
        if file_basename.endswith('.cf'):
            os.remove(self.npol_file)
            #print('file has been removed...{self.npol_file}')
    
        #Py-Art grid output will be put in a tempoarilly created dir within
        #main output dir, when this method is done the full grid
        #.nc file will get moved up to this main output dir and the 
        #temporary dir will get deleted.
        #save the full grid .nc file
        full_grid_filename = f'{file_basename.split(".")[0]}_cntron_{self.center_on}'
        full_grid_name = f'{self.full_grid_dir}{full_grid_filename}.nc'
        DS = grid.to_xarray()
        #pdb.set_trace()
        #print(f'printing DS.x:    {DS.x}')
        #print(f'printing DS.lon:  {DS.lon}')
        #DS = DS.swap_dims({"x": "lon"})
        #DS = DS.swap_dims({"y": "lat"})
        DS.to_netcdf(full_grid_name, format='NETCDF4')
    
        #extract fields & attributes 
        radar_x = DS.x
        radar_y = DS.y
        radar_z = DS.z
        radar_time = DS.time
        if 'ZZ' in DS.keys():
            radar_ZZ = DS.ZZ
        elif 'DZ' in DS.keys():
            radar_ZZ = DS.DZ
        radar_CZ = DS.CZ
        radar_DR = DS.DR
        radar_RH = DS.RH
        radar_PH = DS.PH
        radar_KD = DS.KD
        if 'SQ' in DS.keys(): 
            radar_SQ = DS.SQ
        radar_SW = DS.SW
        radar_VR = DS.VR
        if 'RR' in DS.keys():
            radar_RR = DS.RR
        radar_RP = DS.RP
        radar_RC = DS.RC
        if 'D0' in DS.keys(): 
            radar_D0 = DS.D0
        radar_NW = DS.NW
        radar_FH = DS.FH
        if 'N2' in DS.keys():
            radar_N2 = DS.N2
        radar_DM = DS.DM
        if 'FZ' in DS.keys():
            radar_FZ = DS.FZ
    
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
    
        #radar_ZZ_min -- dont exist
        #radar_ZZ_max -- dont exist
        if 'ZZ' in DS.keys():
            radar_ZZ_badval = DS.ZZ.attrs['_FillValue']
        elif 'DZ' in DS.keys():
            radar_ZZ_badval = DS.DZ.attrs['_FillValue']
        #radar_ZZ_name = DS.ZZ.attrs['long_name'] #note that some fields need have
        #radar_ZZ_units = DS.ZZ.attrs['units']  #names &/or units manually set...
    
        #radar_CZ_min -- dont exist
        #radar_CZ_max -- dont exist
        radar_CZ_badval = DS.CZ.attrs['_FillValue']
        self.radar_CZ_name = DS.CZ.attrs['long_name'] 
        radar_CZ_units = DS.CZ.attrs['units']
    
        #radar_DR_min -- dont exist
        #radar_DR_max -- dont exist
        radar_DR_badval = DS.DR.attrs['_FillValue']
        self.radar_DR_name = DS.DR.attrs['long_name'] 
        self.radar_DR_units = DS.DR.attrs['units']
    
        #radar_RH_min -- dont exist
        #radar_RH_max -- dont exist
        radar_RH_badval = DS.RH.attrs['_FillValue']
        self.radar_RH_name = DS.RH.attrs['long_name'] 
        #radar_RH_units = DS.RH.attrs['units']
    
        #radar_PH_min -- dont exist
        #radar_PH_max -- dont exist
        radar_PH_badval = DS.PH.attrs['_FillValue']
        self.radar_PH_name = DS.PH.attrs['long_name'] 
        self.radar_PH_units = DS.PH.attrs['units']
    
        #radar_KD_min -- dont exist
        #radar_KD_max -- dont exist
        radar_KD_badval = DS.KD.attrs['_FillValue']
        self.radar_KD_name = DS.KD.attrs['long_name'] 
        self.radar_KD_units = DS.KD.attrs['units']
    
        if 'SQ' in DS.keys():
            #radar_SQ_min -- dont exist
            #radar_SQ_max -- dont exist
            radar_SQ_badval = DS.SQ.attrs['_FillValue']
            #radar_SQ_name = DS.SQ.attrs['long_name'] 
            #radar_SQ_units = DS.SQ.attrs['units']
        
        #radar_SW_min -- dont exist
        #radar_SW_max -- dont exist
        radar_SW_badval = DS.SW.attrs['_FillValue']
        self.radar_SW_name = DS.SW.attrs['long_name'] 
        self.radar_SW_units = DS.SW.attrs['units']
    
        #radar_VR_min -- dont exist
        #radar_VR_max -- dont exist
        radar_VR_badval = DS.VR.attrs['_FillValue']
        self.radar_VR_name = DS.VR.attrs['long_name'] 
        self.radar_VR_units = DS.VR.attrs['units']
    
        if 'RR' in DS.keys():
            #radar_RR_min -- dont exist
            #radar_RR_max -- dont exist
            radar_RR_badval = DS.RR.attrs['_FillValue']
            #radar_RR_name = DS.RR.attrs['long_name'] 
            #radar_RR_units = DS.RR.attrs['units']
        
        #radar_RP_min -- dont exist
        #radar_RP_max -- dont exist
        radar_RP_badval = DS.RP.attrs['_FillValue']
        #radar_RP_name = DS.RP.attrs['long_name'] 
        #radar_RP_units = DS.RP.attrs['units']
    
        #radar_RC_min -- dont exist
        #radar_RC_max -- dont exist
        radar_RC_badval = DS.RC.attrs['_FillValue']
        #radar_RC_name = DS.RC.attrs['long_name'] 
        #radar_RC_units = DS.RC.attrs['units']
    
        if 'D0' in DS.keys():
            #p=0
            #radar_D0_min -- dont exist
            #radar_D0_max -- dont exist
            radar_D0_badval = DS.D0.attrs['_FillValue']
            #radar_D0_name = DS.D0.attrs['long_name'] 
            #radar_D0_units = DS.D0.attrs['units']
        
        #radar_NW_min -- dont exist
        #radar_NW_max -- dont exist
        radar_NW_badval = DS.NW.attrs['_FillValue']
        #radar_NW_name = DS.NW.attrs['long_name'] 
        #radar_NW_units = DS.NW.attrs['units']
    
        #radar_FH_min -- dont exist
        #radar_FH_max -- dont exist
        radar_FH_badval = DS.FH.attrs['_FillValue']
        #radar_FH_name = DS.FH.attrs['long_name'] 
        #radar_FH_units = DS.FH.attrs['units']
    
        if 'N2' in DS.keys():
            #p=0
            #radar_N2_min -- dont exist
            #radar_N2_max -- dont exist
            radar_N2_badval = DS.N2.attrs['_FillValue']
            #radar_N2_name = DS.N2.attrs['long_name'] 
            #radar_N2_units = DS.N2.attrs['units']
        
        #radar_DM_min -- dont exist
        #radar_DM_max -- dont exist
        radar_DM_badval = DS.DM.attrs['_FillValue']
        #radar_DM_name = DS.DM.attrs['long_name'] 
        #radar_DM_units = DS.DM.attrs['units']
    
        if 'FZ' in DS.keys():
            #p=0
            #radar_FZ_min -- dont exist
            #radar_FZ_max -- dont exist
            radar_FZ_badval = DS.FZ.attrs['_FillValue']
            #radar_FZ_name = DS.FZ.attrs['long_name'] 
            #radar_FZ_units = DS.FZ.attrs['units']
        #pdb.set_trace()
        
        # Convert string attributes back into string data types:
        # most already are string type
        # Set these manually since names/units from Py-Art don't always default to what
        # NPOL fields actually are (see NPOL fields README from Jason Pippitt)
        self.radar_ZZ_name = 'uncorrected reflectivity' 
        self.radar_CZ_name = 'corrected reflectivity'
        self.radar_SQ_name = 'signal quality index'
        self.radar_RR_name = 'rain rate via DROPS2'
        self.radar_RP_name = 'rain rate via PolZ-R'
        self.radar_RC_name = 'rain rate via Cifelli et al 2002'
        self.radar_D0_name = 'median drop diameter'
        self.radar_NW_name = 'normalized intercept parameter (Dm)'
        self.radar_FH_name = 'hydrometeor ID'
        self.radar_N2_name = 'normalized intercept parameter (Do)'
        self.radar_DM_name = 'mass weighted mean diameter'
        self.radar_FZ_name = 's-ku frequency corrected reflectivity'
        
        self.radar_Z_units = 'dBZ'
        self.radar_SQ_units = 'unitless'
        self.radar_RH_units = 'unitless'
        self.radar_R_units = 'mm h^-1'
        self.radar_D_units = 'mm'
        self.radar_NW_units = 'log(Nw)'
        self.radar_FH_units = 'categorical'
        self.radar_N2_units = 'log(Nw)'
    
        # Replace flagged/bad data values in each array with NaNs so
        # not using the bad data flag value in interpolation
        radar_ZZ = sim.bad2nan(radar_ZZ, radar_ZZ_badval)
        radar_CZ = sim.bad2nan(radar_CZ, radar_CZ_badval)
        radar_DR = sim.bad2nan(radar_DR, radar_DR_badval)
        radar_RH = sim.bad2nan(radar_RH, radar_RH_badval)
        radar_PH = sim.bad2nan(radar_PH, radar_PH_badval)
        radar_KD = sim.bad2nan(radar_KD, radar_KD_badval)
        if 'SQ' in DS.keys(): radar_SQ = sim.bad2nan(radar_SQ, radar_SQ_badval)
        radar_SW = sim.bad2nan(radar_SW, radar_SW_badval)
        radar_VR = sim.bad2nan(radar_VR, radar_VR_badval)
        if 'RR' in DS.keys(): radar_RR = sim.bad2nan(radar_RR, radar_RR_badval)
        radar_RP = sim.bad2nan(radar_RP, radar_RP_badval)
        radar_RC = sim.bad2nan(radar_RC, radar_RC_badval)
        if 'D0' in DS.keys(): radar_D0 = sim.bad2nan(radar_D0, radar_D0_badval)
        radar_NW = sim.bad2nan(radar_NW, radar_NW_badval)
        radar_FH = sim.bad2nan(radar_FH, radar_FH_badval)
        if 'N2' in DS.keys(): radar_N2 = sim.bad2nan(radar_N2, radar_N2_badval)
        radar_DM = sim.bad2nan(radar_DM, radar_DM_badval)
        if 'FZ' in DS.keys(): radar_FZ = sim.bad2nan(radar_FZ, radar_FZ_badval)
    
        #Py-ART gridding adds values at height 0 km so need to force these to NaNs
        zero_ht_sub = np.where(radar_z == 0.0)[0]
        if zero_ht_sub.shape[0] > 0:
            if zero_ht_sub.shape[0] > 1:
                sys.exit('---MORE THAN 1 HEIGHT ZERO IN THE GRIDDED FILE!---')
            else:
                # there is only 1 height sub for 0 km AGL
                # set all field values at ht of 0 km to NaNs
                # dimensions for radar_xx are [time, z, lat, lon] --> [time, z, y, x]; time dim is 1
                radar_ZZ[0,zero_ht_sub,:,:] = np.nan; radar_CZ[0,zero_ht_sub,:,:] = np.nan
                radar_DR[0,zero_ht_sub,:,:] = np.nan; radar_RH[0,zero_ht_sub,:,:] = np.nan
                radar_PH[0,zero_ht_sub,:,:] = np.nan; radar_KD[0,zero_ht_sub,:,:] = np.nan
                radar_SW[0,zero_ht_sub,:,:] = np.nan; radar_NW[0,zero_ht_sub,:,:] = np.nan
                radar_VR[0,zero_ht_sub,:,:] = np.nan; radar_RC[0,zero_ht_sub,:,:] = np.nan
                radar_RP[0,zero_ht_sub,:,:] = np.nan; radar_DM[0,zero_ht_sub,:,:] = np.nan
                radar_FH[0,zero_ht_sub,:,:] = np.nan
                if 'SQ' in DS.keys(): radar_SQ[0,zero_ht_sub,:,:] = np.nan
                if 'RR' in DS.keys(): radar_RR[0,zero_ht_sub,:,:] = np.nan
                if 'D0' in DS.keys(): radar_D0[0,zero_ht_sub,:,:] = np.nan
                if 'N2' in DS.keys(): radar_N2[0,zero_ht_sub,:,:] = np.nan
                if 'FZ' in DS.keys(): radar_FZ[0,zero_ht_sub,:,:] = np.nan
    
        # Locate where column box data is within the full grid of data, 
        #  then pull out the subset of data within the column box grid:
    
        # for vertical direction, # of horiz boxes does not affect subs needed:
        pull_z_start = np.where(radar_z == (min(self.z_values)))[0]
        pull_z_end   = np.where(radar_z == (max(self.z_values)))[0]
        if pull_z_start.shape[0] != 1 or pull_z_end.shape[0] != 1:
            print('---PROBLEM ASSIGNING z dim SUBSCRIPT(S) FOR PULLING RADAR DATA!!')
            print(f'--- min ht [km]: {min(self.z_values)}')
            print(f'--- max ht [km]: {max(self.z_values)}')
            print(f'---z array [km]: {radar_z}')
            sys.exit('Please check and try again!!!')
        # apparantly, these get made as 1-D arrays, but must use scalars as subs,
        # so fix here. have already checked that _start & _end have only one 
        # element (n_zs_match & n_ze_match both are 1 if get to this line)
        # need to add +1 as in python subsetting end index is not inclusive
        pull_z_start = pull_z_start[0]
        pull_z_end   = pull_z_end[0]+1
        #print, 'z subs: ',string(pull_z_start, format='(i0)'),+'  ,  '+string(pull_z_end, format='(i0)')
    
        if self.n_horiz_grid_boxes % 2 == 0:
            #  even # of column grid boxes, so grid box edges line up exactly
            #  column grid data is direct subset in the pyart grid
      
            # get subscripts for center of the pyart grid
            # this point is along edges in the column box grid
            x_origin_sub = np.where(radar_x == 0.0)[0][0]
            y_origin_sub = np.where(radar_x == 0.0)[0][0]
      
            # get # of subscripts before/after
            column_box_horiz_extent_from_cntr = int(self.n_horiz_grid_boxes/2)
    
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
            self.col_radar_ZZ = radar_ZZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_CZ = radar_CZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_DR = radar_DR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_RH = radar_RH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_PH = radar_PH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_KD = radar_KD[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            if 'SQ' in DS.keys(): self.col_radar_SQ = radar_SQ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_SW = radar_SW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_VR = radar_VR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            if 'RR' in DS.keys(): self.col_radar_RR = radar_RR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_RP = radar_RP[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_RC = radar_RC[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]  
            if 'D0' in DS.keys(): self.col_radar_D0 = radar_D0[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_NW = radar_NW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_FH = radar_FH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            if 'N2' in DS.keys(): self.col_radar_N2 = radar_N2[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_radar_DM = radar_DM[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]  
            if 'FZ' in DS.keys(): self.col_radar_FZ = radar_FZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
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
            n_edges = self.n_horiz_grid_boxes + 1                       #convert from ODD to EVEN NUMBER
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
            self.col_radar_ZZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_CZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_DR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_RH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_PH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_KD = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'SQ' in DS.keys(): self.col_radar_SQ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_SW = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_VR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'RR' in DS.keys(): self.col_radar_RR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_RP = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_RC = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))  
            if 'D0' in DS.keys(): self.col_radar_D0 = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_NW = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_FH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'N2' in DS.keys(): self.col_radar_N2 = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_radar_DM = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'FZ' in DS.keys(): self.col_radar_FZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
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
            self.col_radar_SQ = self.col_radar_ZZ*0.0
            self.col_radar_SQ = sim.bad2nan(self.col_radar_SQ, 0.0)
        if 'RR' not in DS.keys():
            self.col_radar_RR = self.col_radar_ZZ*0.0
            self.col_radar_RR = sim.bad2nan(self.col_radar_RR, 0.0)
        if 'D0' not in DS.keys():
            #print('no D0')
            self.col_radar_D0 = self.col_radar_ZZ*0.0
            self.col_radar_D0 = sim.bad2nan(self.col_radar_D0, 0.0)
        if 'N2' not in DS.keys():
            #print('no N2')
            self.col_radar_N2 = self.col_radar_ZZ*0.0
            self.col_radar_N2 = sim.bad2nan(self.col_radar_N2, 0.0)
        if 'FZ' not in DS.keys():
            self.col_radar_FZ = self.col_radar_ZZ*0.0
            self.col_radar_FZ = sim.bad2nan(self.col_radar_FZ, 0.0)
            
        self.npol_data_in_column = True

    # ***************************************************************************************
    def grid_lev2_for_column(self):
        '''
        ;------------------------------------------------------------
        ; grid lev2 NEXRAD radar data to a 200 x 200 x 15 km grid centered on the column box center
        ;  save full gridded file, then return subset within the column box
        ;------------------------------------------------------------
        ;   -> save the full gridded file to the input full_grid_dir
        ;
        ;   -> take subset of the gridded data that is within the column box and
        ;	send this subset back to the top-level build_column.pro for inclusion
        ;	in the new column .nc file
        ;
        ;  result = grid_npol_for_column(npol_plat_info, orig_NPOL_file, $
        ;				params_file, column_box_params, $
        ;				full_grid_dir)
        ;	npol_plat_info:	 plat_info structure returned by set_plat_values.pro
        ;	orig_NPOL_file:	 original NPOL data file, *.uf.gz
        ;	params_file:	 full path to Radx .params file to use for gridding
        ;	column_box_params:  column box grid parameters.  This needs to be the 
        ;				full structure returned by define_box.pro
        ;	full_grid_dir:   path to dir for full grid Radx output
        ;
        ;	result:	  structure containing the NPOL data at each location in
        ;		  in the column box grid (whether column grid is along or
        ;		  in middle of Radx/radar grid in the horizontal)
        ;	  Tags:
        ;		.**name:  string name of ** data field, as recorded in the Radx
        ;			   gridded .nc file (may be an empty string)
        ;		.**units: string telling the ** data field units, as recorded in
        ;		     	   the Radx gridded .nc file (may be an empty string)
        ;		.**data:  NPOL ** field value at locations in column box grid
        ;		    this array is 3-D: [X x Y x Z] dims, will match col box dims
        ;		    middle point of column box = middle point of this array
        ;		    IF EVEN # OF HORIZ COLUMN GRID SUBBOXES: 
        ;		    	middle point is along an edge
        ;		    IF ODD # OF HORIZ COLUMB GRID SUBBOXES:
        ;		        middle poit is at middle of a grid box
        ;
        ;     Have 3 Tags for each NPOL field getting put into the column:
        ;	
        ;    .ZZname    .CZname    .DRname    .RHname    .PHname    .KDname    .SQname
        ;    .ZZunits   .CZunits   .DRunits   .RHunits   .PHunits   .KDunits   .SQunits
        ;    .ZZdata    .CZdata    .DRdata    .RHdata    .PHdata    .KDdata    .SQdata
        ;
        ;    .SWname	.VRname    .RRname    .RPname	 .RCname    .D0name    .NWname
        ;    .SWunits	.VRunits   .RRunits   .RPunits	 .RCunits   .D0units   .NWunits
        ;    .SWdata	.VRdata    .RRdata    .RPdata	 .RCdata    .D0data    .NWdata
        ;
        ;    .FHname    .N2name	   .DMname
        ;    .FHunits   .N2units   .DMunits
        ;    .FHdata    .N2data    .DMdata
        ;
        ;	Dependencies:
        ;	  Radx software		UCAR/NCAR software for working with radar data
        ;	  uncomp_file.pro	utility to uncompress a .gz file
        ;	  bad2nan.pro		utility to convert a flag value to NaNs
        ;	  
        ;	  Updates:
        ;	    added code to convert from .cf --> .uf due to SIMBA/Radx gridding failing
        ;	    after 07/30/18, see lines 85-101 - C. Pabla (02/06/19)
        ;
        ;------------------------------------------------------------
        '''
        #------prepare data for gridding-------
    
        #unzip file (if needed)
        #file_basename = os.path.basename(self.nexrad_file)
        #if file_basename.endswith('.gz'):
        #   cf_file = sim.ungzip_file(self.nexrad_file)
        #else:
        #    cf_file = self.nexrad_file
        
        #fields = ['ZZ','CZ','DR','RH','PH','KD','SQ','SW','VR',\
        #          'RR','RP','RC','D0','NW','FH','N2','DM']
        grid_center_lat_lon = self.cntr_lat_deg, self.cntr_lon_deg
        print('...running Py-ART gridding...')
        grid = self.get_grid_from_radar(self.nexrad_file, grid_center_lat_lon)
        
        #remove unzipped file if needed
        file_basename = os.path.basename(self.nexrad_file)
        if file_basename.endswith('.cf'):
            os.remove(self.nexrad_file)
            #print(f'file has been removed...{self.nexrad_file}')
    
        #Radx output will be put in a tempoarilly created dir within
        #main output dir, when this method is done the full grid
        #.nc file will get moved up to this main output dir and the 
        #temporary dir will get deleted.
        #save the full grid .nc file
        full_grid_filename = f'{file_basename.split(".")[0]}_cntron_{self.center_on}'
        full_grid_name = f'{self.full_grid_dir}{full_grid_filename}.nc'
        DS = grid.to_xarray()
        #DS = DS.swap_dims({"x": "lon"})
        #DS = DS.swap_dims({"y": "lat"})
        DS.to_netcdf(full_grid_name, format='NETCDF4')
    
        #extract fields & attributes 
        radar_x = DS.x
        radar_y = DS.y
        radar_z = DS.z
        radar_time = DS.time
        if 'ZZ' in DS.keys():
            radar_ZZ = DS.ZZ
        elif 'DZ' in DS.keys():
            radar_ZZ = DS.DZ
        radar_CZ = DS.CZ
        radar_DR = DS.DR
        radar_RH = DS.RH
        radar_PH = DS.PH
        radar_KD = DS.KD
        if 'SQ' in DS.keys(): 
            radar_SQ = DS.SQ
        radar_SW = DS.SW
        radar_VR = DS.VR
        if 'RR' in DS.keys():
            radar_RR = DS.RR
        radar_RP = DS.RP
        radar_RC = DS.RC
        if 'D0' in DS.keys(): 
            radar_D0 = DS.D0
        radar_NW = DS.NW
        radar_FH = DS.FH
        if 'N2' in DS.keys():
            radar_N2 = DS.N2
        radar_DM = DS.DM
        if 'FZ' in DS.keys():
            radar_FZ = DS.FZ
    
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
    
        #radar_ZZ_min -- dont exist
        #radar_ZZ_max -- dont exist
        if 'ZZ' in DS.keys():
            radar_ZZ_badval = DS.ZZ.attrs['_FillValue']
        elif 'DZ' in DS.keys():
            radar_ZZ_badval = DS.DZ.attrs['_FillValue']
        #radar_ZZ_name = DS.ZZ.attrs['long_name'] #note that some fields need have
        #radar_ZZ_units = DS.ZZ.attrs['units']  #names &/or units manually set...
    
        #radar_CZ_min -- dont exist
        #radar_CZ_max -- dont exist
        radar_CZ_badval = DS.CZ.attrs['_FillValue']
        self.lev2_CZ_name = DS.CZ.attrs['long_name'] 
        radar_CZ_units = DS.CZ.attrs['units']
    
        #radar_DR_min -- dont exist
        #radar_DR_max -- dont exist
        radar_DR_badval = DS.DR.attrs['_FillValue']
        self.lev2_DR_name = DS.DR.attrs['long_name'] 
        self.lev2_DR_units = DS.DR.attrs['units']
    
        #radar_RH_min -- dont exist
        #radar_RH_max -- dont exist
        radar_RH_badval = DS.RH.attrs['_FillValue']
        self.lev2_RH_name = DS.RH.attrs['long_name'] 
        #radar_RH_units = DS.RH.attrs['units']
    
        #radar_PH_min -- dont exist
        #radar_PH_max -- dont exist
        radar_PH_badval = DS.PH.attrs['_FillValue']
        self.lev2_PH_name = DS.PH.attrs['long_name'] 
        self.lev2_PH_units = DS.PH.attrs['units']
    
        #radar_KD_min -- dont exist
        #radar_KD_max -- dont exist
        radar_KD_badval = DS.KD.attrs['_FillValue']
        self.lev2_KD_name = DS.KD.attrs['long_name'] 
        self.lev2_KD_units = DS.KD.attrs['units']
    
        if 'SQ' in DS.keys():
            #radar_SQ_min -- dont exist
            #radar_SQ_max -- dont exist
            radar_SQ_badval = DS.SQ.attrs['_FillValue']
            #radar_SQ_name = DS.SQ.attrs['long_name'] 
            #radar_SQ_units = DS.SQ.attrs['units']
        
        #radar_SW_min -- dont exist
        #radar_SW_max -- dont exist
        radar_SW_badval = DS.SW.attrs['_FillValue']
        self.lev2_SW_name = DS.SW.attrs['long_name'] 
        self.lev2_SW_units = DS.SW.attrs['units']
    
        #radar_VR_min -- dont exist
        #radar_VR_max -- dont exist
        radar_VR_badval = DS.VR.attrs['_FillValue']
        self.lev2_VR_name = DS.VR.attrs['long_name'] 
        self.lev2_VR_units = DS.VR.attrs['units']
    
        if 'RR' in DS.keys():
            p=0
            #radar_RR_min -- dont exist
            #radar_RR_max -- dont exist
            radar_RR_badval = DS.RR.attrs['_FillValue']
            #radar_RR_name = DS.RR.attrs['long_name'] 
            #radar_RR_units = DS.RR.attrs['units']
        
        #radar_RP_min -- dont exist
        #radar_RP_max -- dont exist
        radar_RP_badval = DS.RP.attrs['_FillValue']
        #radar_RP_name = DS.RP.attrs['long_name'] 
        #radar_RP_units = DS.RP.attrs['units']
    
        #radar_RC_min -- dont exist
        #radar_RC_max -- dont exist
        radar_RC_badval = DS.RC.attrs['_FillValue']
        #radar_RC_name = DS.RC.attrs['long_name'] 
        #radar_RC_units = DS.RC.attrs['units']
    
        if 'D0' in DS.keys():
            #p=0
            #radar_D0_min -- dont exist
            #radar_D0_max -- dont exist
            radar_D0_badval = DS.D0.attrs['_FillValue']
            #radar_D0_name = DS.D0.attrs['long_name'] 
            #radar_D0_units = DS.D0.attrs['units']
        
        #radar_NW_min -- dont exist
        #radar_NW_max -- dont exist
        radar_NW_badval = DS.NW.attrs['_FillValue']
        #radar_NW_name = DS.NW.attrs['long_name'] 
        #radar_NW_units = DS.NW.attrs['units']
    
        #radar_FH_min -- dont exist
        #radar_FH_max -- dont exist
        radar_FH_badval = DS.FH.attrs['_FillValue']
        #radar_FH_name = DS.FH.attrs['long_name'] 
        #radar_FH_units = DS.FH.attrs['units']
    
        if 'N2' in DS.keys():
            #p=0
            #radar_N2_min -- dont exist
            #radar_N2_max -- dont exist
            radar_N2_badval = DS.N2.attrs['_FillValue']
            #radar_N2_name = DS.N2.attrs['long_name'] 
            #radar_N2_units = DS.N2.attrs['units']
        
        #radar_DM_min -- dont exist
        #radar_DM_max -- dont exist
        radar_DM_badval = DS.DM.attrs['_FillValue']
        #radar_DM_name = DS.DM.attrs['long_name'] 
        #radar_DM_units = DS.DM.attrs['units']
    
        if 'FZ' in DS.keys():
            #p=0
            #radar_FZ_min -- dont exist
            #radar_FZ_max -- dont exist
            radar_FZ_badval = DS.FZ.attrs['_FillValue']
            #radar_FZ_name = DS.FZ.attrs['long_name'] 
            #radar_FZ_units = DS.FZ.attrs['units']
        
        # Convert string attributes back into string data types:
        # most already are string type
        # Set these manually since names/units from Py-Art don't always default to what
        # NPOL fields actually are (see NPOL fields README from Jason Pippitt)
        self.lev2_ZZ_name = 'uncorrected reflectivity' 
        self.lev2_CZ_name = 'corrected reflectivity'
        self.lev2_SQ_name = 'signal quality index'
        self.lev2_RR_name = 'rain rate via DROPS2'
        self.lev2_RP_name = 'rain rate via PolZ-R'
        self.lev2_RC_name = 'rain rate via Cifelli et al 2002'
        self.lev2_D0_name = 'median drop diameter'
        self.lev2_NW_name = 'normalized intercept parameter (Dm)'
        self.lev2_FH_name = 'hydrometeor ID'
        self.lev2_N2_name = 'normalized intercept parameter (Do)'
        self.lev2_DM_name = 'mass weighted mean diameter'
        self.lev2_FZ_name = 's-ku frequency corrected reflectivity'
        
        self.lev2_Z_units = 'dBZ'
        self.lev2_SQ_units = 'unitless'
        self.lev2_RH_units = 'unitless'
        self.lev2_R_units = 'mm h^-1'
        self.lev2_D_units = 'mm'
        self.lev2_NW_units = 'log(Nw)'
        self.lev2_FH_units = 'categorical'
        self.lev2_N2_units = 'log(Nw)'
    
        # Replace flagged/bad data values in each array with NaNs so
        # not using the bad data flag value in interpolation
        radar_ZZ = sim.bad2nan(radar_ZZ, radar_ZZ_badval)
        radar_CZ = sim.bad2nan(radar_CZ, radar_CZ_badval)
        radar_DR = sim.bad2nan(radar_DR, radar_DR_badval)
        radar_RH = sim.bad2nan(radar_RH, radar_RH_badval)
        radar_PH = sim.bad2nan(radar_PH, radar_PH_badval)
        radar_KD = sim.bad2nan(radar_KD, radar_KD_badval)
        if 'SQ' in DS.keys(): radar_SQ = sim.bad2nan(radar_SQ, radar_SQ_badval)
        radar_SW = sim.bad2nan(radar_SW, radar_SW_badval)
        radar_VR = sim.bad2nan(radar_VR, radar_VR_badval)
        if 'RR' in DS.keys(): radar_RR = sim.bad2nan(radar_RR, radar_RR_badval)
        radar_RP = sim.bad2nan(radar_RP, radar_RP_badval)
        radar_RC = sim.bad2nan(radar_RC, radar_RC_badval)
        if 'D0' in DS.keys(): radar_D0 = sim.bad2nan(radar_D0, radar_D0_badval)
        radar_NW = sim.bad2nan(radar_NW, radar_NW_badval)
        radar_FH = sim.bad2nan(radar_FH, radar_FH_badval)
        if 'N2' in DS.keys(): radar_N2 = sim.bad2nan(radar_N2, radar_N2_badval)
        radar_DM = sim.bad2nan(radar_DM, radar_DM_badval)
        if 'FZ' in DS.keys(): radar_FZ = sim.bad2nan(radar_FZ, radar_FZ_badval)
    
        #Py-ART gridding adds values at height 0 km so need to force these to NaNs
        zero_ht_sub = np.where(radar_z == 0.0)[0]
        if zero_ht_sub.shape[0] > 0:
            if zero_ht_sub.shape[0] > 1:
                sys.exit('---MORE THAN 1 HEIGHT ZERO IN THE GRIDDED FILE!---')
            else:
                # there is only 1 height sub for 0 km AGL
                # set all field values at ht of 0 km to NaNs
                # dimensions for radar_xx are [time, z, lat, lon] --> [time, z, y, x]; time dim is 1
                radar_ZZ[0,zero_ht_sub,:,:] = np.nan; radar_CZ[0,zero_ht_sub,:,:] = np.nan
                radar_DR[0,zero_ht_sub,:,:] = np.nan; radar_RH[0,zero_ht_sub,:,:] = np.nan
                radar_PH[0,zero_ht_sub,:,:] = np.nan; radar_KD[0,zero_ht_sub,:,:] = np.nan
                radar_SW[0,zero_ht_sub,:,:] = np.nan; radar_NW[0,zero_ht_sub,:,:] = np.nan
                radar_VR[0,zero_ht_sub,:,:] = np.nan; radar_RC[0,zero_ht_sub,:,:] = np.nan
                radar_RP[0,zero_ht_sub,:,:] = np.nan; radar_DM[0,zero_ht_sub,:,:] = np.nan
                radar_FH[0,zero_ht_sub,:,:] = np.nan
                if 'SQ' in DS.keys(): radar_SQ[0,zero_ht_sub,:,:] = np.nan
                if 'RR' in DS.keys(): radar_RR[0,zero_ht_sub,:,:] = np.nan
                if 'D0' in DS.keys(): radar_D0[0,zero_ht_sub,:,:] = np.nan
                if 'N2' in DS.keys(): radar_N2[0,zero_ht_sub,:,:] = np.nan
                if 'FZ' in DS.keys(): radar_FZ[0,zero_ht_sub,:,:] = np.nan
    
        # Locate where column box data is within the full grid of data, 
        #  then pull out the subset of data within the column box grid:
    
        # for vertical direction, # of horiz boxes does not affect subs needed:
        pull_z_start = np.where(radar_z == (min(self.z_values)))[0]
        pull_z_end   = np.where(radar_z == (max(self.z_values)))[0]
        if pull_z_start.shape[0] != 1 or pull_z_end.shape[0] != 1:
            print('---PROBLEM ASSIGNING z dim SUBSCRIPT(S) FOR PULLING RADAR DATA!!')
            print(f'--- min ht [km]: {min(self.z_values)}')
            print(f'--- max ht [km]: {max(self.z_values)}')
            print(f'---z array [km]: {radar_z}')
            sys.exit('Please check and try again!!!')
        # apparantly, these get made as 1-D arrays, but must use scalars as subs,
        # so fix here. have already checked that _start & _end have only one 
        # element (n_zs_match & n_ze_match both are 1 if get to this line)
        # need to add +1 as in python subsetting end index is not inclusive
        pull_z_start = pull_z_start[0]
        pull_z_end   = pull_z_end[0]+1
        #print, 'z subs: ',string(pull_z_start, format='(i0)'),+'  ,  '+string(pull_z_end, format='(i0)')
    
        if self.n_horiz_grid_boxes % 2 == 0:
            #  even # of column grid boxes, so grid box edges line up exactly
            #  column grid data is direct subset in the pyart grid
      
            # get subscripts for center of the pyart grid
            # this point is along edges in the column box grid
            x_origin_sub = np.where(radar_x == 0.0)[0][0]
            y_origin_sub = np.where(radar_x == 0.0)[0][0]
      
            # get # of subscripts before/after
            column_box_horiz_extent_from_cntr = int(self.n_horiz_grid_boxes/2)
    
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
            self.col_lev2_ZZ = radar_ZZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_CZ = radar_CZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_DR = radar_DR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_RH = radar_RH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_PH = radar_PH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_KD = radar_KD[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            if 'SQ' in DS.keys(): self.col_lev2_SQ = radar_SQ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_SW = radar_SW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_VR = radar_VR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            if 'RR' in DS.keys(): self.col_lev2_RR = radar_RR[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_RP = radar_RP[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_RC = radar_RC[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]  
            if 'D0' in DS.keys(): self.col_lev2_D0 = radar_D0[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_NW = radar_NW[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_FH = radar_FH[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            if 'N2' in DS.keys(): self.col_lev2_N2 = radar_N2[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
            self.col_lev2_DM = radar_DM[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]  
            if 'FZ' in DS.keys(): self.col_lev2_FZ = radar_FZ[0, pull_z_start:pull_z_end, pull_y_start:pull_y_end, pull_x_start:pull_x_end]
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
            n_edges = self.n_horiz_grid_boxes + 1                       #convert from ODD toEVEN NUMBER
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
            self.col_lev2_ZZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_CZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_DR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_RH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_PH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_KD = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'SQ' in DS.keys(): self.col_lev2_SQ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_SW = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_VR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'RR' in DS.keys(): self.col_lev2_RR = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_RP = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_RC = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))  
            if 'D0' in DS.keys(): self.col_lev2_D0 = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_NW = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_FH = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'N2' in DS.keys(): self.col_lev2_N2 = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            self.col_lev2_DM = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
            if 'FZ' in DS.keys(): self.col_lev2_FZ = np.zeros((len(new_x_subs),len(new_y_subs),pull_z_end - pull_z_start +1))
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
            self.col_lev2_SQ = self.col_lev2_ZZ*0.0
            self.col_lev2_SQ = sim.bad2nan(self.col_lev2_SQ, 0.0)
        if 'RR' not in DS.keys():
            self.col_lev2_RR = self.col_lev2_ZZ*0.0
            self.col_lev2_RR = sim.bad2nan(self.col_lev2_RR, 0.0)
        if 'D0' not in DS.keys():
            self.col_lev2_D0 = self.col_lev2_ZZ*0.0
            self.col_lev2_D0 = sim.bad2nan(self.col_lev2_D0, 0.0)
        if 'N2' not in DS.keys():
            self.col_lev2_N2 = self.col_lev2_ZZ*0.0
            self.col_lev2_N2 = sim.bad2nan(self.col_lev2_N2, 0.0)
        if 'FZ' not in DS.keys():
            self.col_lev2_FZ = self.col_lev2_ZZ*0.0
            self.col_lev2_FZ = sim.bad2nan(self.col_lev2_FZ, 0.0)

        self.lev2_data_in_column = True

    # ***************************************************************************************
    def write_column_nc_file(self):
    
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
        column_nc_fname = 'column_'+self.main_plat_name+'_'+ \
                self.center_on+'_'+self.main_plat_timestamp[0:13]+'.nc'
        header_fname = os.path.basename(column_nc_fname)+'_ncHeader.txt'
        print('------------------------------------')
        print('---Now Writing New Column .nc File: ')
        print(f'          {column_nc_fname}')
        column_nc_fname = self.out_dir+column_nc_fname
        header_fname    = self.out_dir+header_fname

        if not self.have_xy_vals:
            print('------------------------------------------------------------------')
            print('---DID NOT GET CORRECT x & y ARRAYS TO USE FOR COLUMN .nc FILE!---')
            print('---------COLUMN NetCDF FILE HAS NOT BEEN GENERATED!!!-------------')
            print('------------------------------------------------------------------')
            return -1
        #print(f'x: {self.x_values}') #just for testing
        #print(f'y: {self.y_values}')


        # Begin new .nc file, only if one with same name doesn't already exist:
        # ***CHANGE THIS TO /NOCLOBBER WHEN OK WITH PERFROMANCE!!!***
        ncfile = Dataset(column_nc_fname,mode='w', noclobber=False)

        # Define the X, Y, & Z dimension sizes for the column box grid
        #  (# of box edges, not the number of boxes)
        xdim = ncfile.createDimension('x', self.lon_values.shape[0])
        ydim = ncfile.createDimension('y', self.lat_values.shape[0])
        zdim = ncfile.createDimension('z', self.z_values.shape[0])
        tdim = ncfile.createDimension('t', self.n_interval_times)
        ndim = ncfile.createDimension('n', None) #this dimension is for dropcounts in 2dvd, apu 
        
        # Set global attributes:  include main plat & info for each plat
        ncfile.setncattr('box_centered_on', self.center_on)
        ncfile.setncattr('box_center_lat', self.lat_values)
        ncfile.setncattr('box_center_lon', self.lon_values)
        ncfile.setncattr('grid_spacing_vert', self.vert_spacing)
        ncfile.setncattr('grid_spacing_horiz', self.box_spacing)
        ncfile.setncattr('grid_extent_vert', self.vert_limit)
        ncfile.setncattr('grid_extent_horiz', self.box_limit)
        ncfile.setncattr('grid_spacing_and_limits_units', 'meters')
        ncfile.setncattr('main_platform', self.main_plat_name)
        ncfile.setncattr('main_plat_timestamp', self.main_plat_timestamp)
        ncfile.setncattr('main_plat_mode', self.main_scan_type)
        ncfile.setncattr('pySIMBA_version', self.version_string)
        
        #write coordinate data
        x_ID = ncfile.createVariable('x', np.float32, ('x'), zlib=True)
        x_ID[:] = self.x_values
        
        y_ID = ncfile.createVariable('y', np.float32, ('y'), zlib=True)
        y_ID[:] = self.y_values
        
        z_ID = ncfile.createVariable('z', np.float32, ('z'), zlib=True)
        z_ID[:] = self.z_values
        
        t_ID = ncfile.createVariable('t', np.float32, ('t'), zlib=True)
        t_ID[:] = self.t_values
        
        latID = ncfile.createVariable('lat', np.float32, ('x'), zlib=True)
        latID[:] = self.lat_values
        
        lonID = ncfile.createVariable('lon', np.float32, ('x'), zlib=True)
        lonID[:] = self.lon_values

        npol_avail = ncfile.createVariable('npol_avail', str , (), zlib=True)
        d3r_avail = ncfile.createVariable('d3r_avail', str , (), zlib=True)
        lev2_avail = ncfile.createVariable('lev2_avail', str , (), zlib=True)
        apu_avail = ncfile.createVariable('apu_avail', str , (), zlib=True)
        twodvd_avail = ncfile.createVariable('twodvd_avail', str , (), zlib=True)
        dow6_avail = ncfile.createVariable('dow6_avail', str , (), zlib=True)
        pluvio_avail = ncfile.createVariable('pluvio_avail', str , (), zlib=True)
        gauges_avail = ncfile.createVariable('gauges_avail', str , (), zlib=True)
        mrr_avail = ncfile.createVariable('mrr_avail', str , (), zlib=True)
        gmi_gprof_avail = ncfile.createVariable('gmi_gprof_avail', str , (), zlib=True)
        gmi_L1C_avail = ncfile.createVariable('gmi_L1C_avail', str , (), zlib=True)
        dpr_avail = ncfile.createVariable('dpr_avail', str , (), zlib=True)
        lev2bcmb_avail = ncfile.createVariable('lev2bcmb_avail', str , (), zlib=True)
        mrms_avail = ncfile.createVariable('mrms_avail', str , (), zlib=True)
        soundings_avail = ncfile.createVariable('soundings_avail', str , (), zlib=True)
        
        #-------------------------------------------------------------------------------
        # Check if each platform data is available to input into column
        #  not available: Set all the attributes to -9999 or 'platfrom_not_avail'
        #  available:     Define variables & coordinates
        #                 Write data to the NETCDF File/Group Defined Above
        #-------------------------------------------------------------------------------
        #
        # DEFINE VARS & SET VAR ATTS FOR NPOL DATA:
        if self.npol_data_in_column:
            npol_avail.setncattr('npol_available', 'True')
            npol_avail = set_nc_radar_attributes(npol_avail, self.npol_info)
            
            npol_ZZ_ID = ncfile.createVariable('npol_ZZ', np.float32, ('z','y','x'), zlib=True)
            npol_ZZ_ID.long_name = self.radar_ZZ_name
            npol_ZZ_ID.units = self.radar_Z_units
            npol_ZZ_ID[:,:,:] = self.col_radar_ZZ
            npol_ZZ_ID.badval = np.nan
        
            npol_CZ_ID = ncfile.createVariable('npol_CZ', np.float32, ('z','y','x'), zlib=True)
            npol_CZ_ID.long_name = self.radar_CZ_name
            npol_CZ_ID.units = self.radar_Z_units
            npol_CZ_ID[:,:,:] = self.col_radar_CZ
            npol_CZ_ID.badval = np.nan
        
            npol_DR_ID = ncfile.createVariable('npol_DR', np.float32, ('z','y','x'), zlib=True)
            npol_DR_ID.long_name = self.radar_DR_name
            npol_DR_ID.units = self.radar_DR_units
            npol_DR_ID[:,:,:] = self.col_radar_DR
            npol_DR_ID.badval = np.nan
        
            npol_RH_ID = ncfile.createVariable('npol_RH', np.float32, ('z','y','x'), zlib=True)
            npol_RH_ID.long_name = self.radar_RH_name
            npol_RH_ID.units = self.radar_RH_units
            npol_RH_ID[:,:,:] = self.col_radar_RH
            npol_RH_ID.badval = np.nan
        
            npol_PH_ID = ncfile.createVariable('npol_PH', np.float32, ('z','y','x'), zlib=True)
            npol_PH_ID.long_name = self.radar_PH_name
            npol_PH_ID.units = self.radar_PH_units
            npol_PH_ID[:,:,:] = self.col_radar_PH
            npol_PH_ID.badval = np.nan
        
            npol_KD_ID = ncfile.createVariable('npol_KD', np.float32, ('z','y','x'), zlib=True)
            npol_KD_ID.long_name = self.radar_KD_name
            npol_KD_ID.units = self.radar_KD_units
            npol_KD_ID[:,:,:] = self.col_radar_KD
            npol_KD_ID.badval = np.nan
        
            npol_SQ_ID = ncfile.createVariable('npol_SQ', np.float32, ('z','y','x'), zlib=True)
            npol_SQ_ID.long_name = self.radar_SQ_name
            npol_SQ_ID.units = self.radar_SQ_units
            npol_SQ_ID[:,:,:] = self.col_radar_SQ
            npol_SQ_ID.badval = np.nan
        
            npol_SW_ID = ncfile.createVariable('npol_SW', np.float32, ('z','y','x'), zlib=True)
            npol_SW_ID.long_name = self.radar_SW_name
            npol_SW_ID.units = self.radar_SW_units
            npol_SW_ID[:,:,:] = self.col_radar_SW
            npol_SW_ID.badval = np.nan
        
            npol_VR_ID = ncfile.createVariable('npol_VR', np.float32, ('z','y','x'), zlib=True)
            npol_VR_ID.long_name = self.radar_VR_name
            npol_VR_ID.units = self.radar_VR_units
            npol_VR_ID[:,:,:] = self.col_radar_VR
            npol_VR_ID.badval = np.nan
        
            npol_RR_ID = ncfile.createVariable('npol_RR', np.float32, ('z','y','x'), zlib=True)
            npol_RR_ID.long_name = self.radar_RR_name
            npol_RR_ID.units = self.radar_R_units
            npol_RR_ID[:,:,:] = self.col_radar_RR
            npol_RR_ID.badval = np.nan
        
            npol_RP_ID = ncfile.createVariable('npol_RP', np.float32, ('z','y','x'), zlib=True)
            npol_RP_ID.long_name = self.radar_RP_name
            npol_RP_ID.units = self.radar_R_units
            npol_RP_ID[:,:,:] = self.col_radar_RP
            npol_RP_ID.badval = np.nan
        
            npol_RC_ID = ncfile.createVariable('npol_RC', np.float32, ('z','y','x'), zlib=True)
            npol_RC_ID.long_name = self.radar_RC_name
            npol_RC_ID.units = self.radar_R_units
            npol_RC_ID[:,:,:] = self.col_radar_RC
            npol_RC_ID.badval = np.nan
        
            npol_D0_ID = ncfile.createVariable('npol_D0', np.float32, ('z','y','x'), zlib=True)
            npol_D0_ID.long_name = self.radar_D0_name
            npol_D0_ID.units = self.radar_D_units
            npol_D0_ID[:,:,:] = self.col_radar_D0
            npol_D0_ID.badval = np.nan
        
            npol_NW_ID = ncfile.createVariable('npol_NW', np.float32, ('z','y','x'), zlib=True)
            npol_NW_ID.long_name = self.radar_NW_name
            npol_NW_ID.units = self.radar_NW_units
            npol_NW_ID[:,:,:] = self.col_radar_NW
            npol_NW_ID.badval = np.nan
        
            npol_FH_ID = ncfile.createVariable('npol_FH', np.float32, ('z','y','x'), zlib=True)
            npol_FH_ID.long_name = self.radar_FH_name
            npol_FH_ID.units = self.radar_FH_units
            npol_FH_ID[:,:,:] = self.col_radar_FH
            npol_FH_ID.badval = np.nan
        
            npol_N2_ID = ncfile.createVariable('npol_N2', np.float32, ('z','y','x'), zlib=True)
            npol_N2_ID.long_name = self.radar_N2_name
            npol_N2_ID.units = self.radar_N2_units
            npol_N2_ID[:,:,:] = self.col_radar_N2
            npol_N2_ID.badval = np.nan
        
            npol_DM_ID = ncfile.createVariable('npol_DM', np.float32, ('z','y','x'), zlib=True)
            npol_DM_ID.long_name = self.radar_DM_name
            npol_DM_ID.units = self.radar_D_units
            npol_DM_ID[:,:,:] = self.col_radar_DM
            npol_DM_ID.badval = np.nan
        
            npol_FZ_ID = ncfile.createVariable('npol_FZ', np.float32, ('z','y','x'), zlib=True)
            npol_FZ_ID.long_name = self.radar_FZ_name
            npol_FZ_ID.units = self.radar_Z_units
            npol_FZ_ID[:,:,:] = self.col_radar_FZ
            npol_FZ_ID.badval = np.nan
        else:
            npol_avail.setncattr('npol_available', 'False')
            npol_avail = set_nc_radar_missing_attributes(npol_avail)
        
    
        #print(npol_CZ_ID)
        #print(npol_CZ_ID.shape)
        #print(npol_grp)
        #print(npol_grp.variables['npol_CZ'][5,:,:])
        #print(npol_grp.variables['npol_CZ'].shape)
        #print(npol_grp.variables['npol_DR'])
        #print(npol_grp.variables['npol_DR'][5,:,:])
        #print(ncfile)
        #ncfile.close()
        #print('Dataset is closed!')
        #pdb.set_trace()

        #
        # DEFINE VARS & SET VAR ATTS FOR D3R DATA:
        # SKIP FOR NOW
        #
    
        #
        # DEFINE VARS & SET VAR ATTS FOR 88D VALUES:
        if self.lev2_data_in_column:
            #lev2Name = self.lev2_plat_info['plat_name']
            #lev2Name = 'KDOX' ##temporary -- need to work on this when platform not available
            lev2_avail.setncattr(f'{self.nexrad_88D_ID.lower()}_available', 'True')
            lev2_avail = set_nc_radar_attributes(lev2_avail, self.lev2_info)
            
            lev2_ZZ_ID = ncfile.createVariable(self.nexrad_88D_ID+'_ZZ', np.float32, ('z','y','x'), zlib=True)
            lev2_ZZ_ID.long_name = self.lev2_ZZ_name
            lev2_ZZ_ID.units = self.lev2_Z_units
            lev2_ZZ_ID[:,:,:] = self.col_lev2_ZZ
            lev2_ZZ_ID.badval = np.nan
        
            lev2_CZ_ID = ncfile.createVariable(self.nexrad_88D_ID+'_CZ', np.float32, ('z','y','x'), zlib=True)
            lev2_CZ_ID.long_name = self.lev2_CZ_name
            lev2_CZ_ID.units = self.lev2_Z_units
            lev2_CZ_ID[:,:,:] = self.col_lev2_CZ
            lev2_CZ_ID.badval = np.nan
        
            lev2_DR_ID = ncfile.createVariable(self.nexrad_88D_ID+'_DR', np.float32, ('z','y','x'), zlib=True)
            lev2_DR_ID.long_name = self.lev2_DR_name
            lev2_DR_ID.units = self.lev2_DR_units
            lev2_DR_ID[:,:,:] = self.col_lev2_DR
            lev2_DR_ID.badval = np.nan
        
            lev2_RH_ID = ncfile.createVariable(self.nexrad_88D_ID+'_RH', np.float32, ('z','y','x'), zlib=True)
            lev2_RH_ID.long_name = self.lev2_RH_name
            lev2_RH_ID.units = self.lev2_RH_units
            lev2_RH_ID[:,:,:] = self.col_lev2_RH
            lev2_RH_ID.badval = np.nan
        
            lev2_PH_ID = ncfile.createVariable(self.nexrad_88D_ID+'_PH', np.float32, ('z','y','x'), zlib=True)
            lev2_PH_ID.long_name = self.lev2_PH_name
            lev2_PH_ID.units = self.lev2_PH_units
            lev2_PH_ID[:,:,:] = self.col_lev2_PH
            lev2_PH_ID.badval = np.nan
        
            lev2_KD_ID = ncfile.createVariable(self.nexrad_88D_ID+'_KD', np.float32, ('z','y','x'), zlib=True)
            lev2_KD_ID.long_name = self.lev2_KD_name
            lev2_KD_ID.units = self.lev2_KD_units
            lev2_KD_ID[:,:,:] = self.col_lev2_KD
            lev2_KD_ID.badval = np.nan
        
            lev2_SQ_ID = ncfile.createVariable(self.nexrad_88D_ID+'_SQ', np.float32, ('z','y','x'), zlib=True)
            lev2_SQ_ID.long_name = self.lev2_SQ_name
            lev2_SQ_ID.units = self.lev2_SQ_units
            lev2_SQ_ID[:,:,:] = self.col_lev2_SQ
            lev2_SQ_ID.badval = np.nan
        
            lev2_SW_ID = ncfile.createVariable(self.nexrad_88D_ID+'_SW', np.float32, ('z','y','x'), zlib=True)
            lev2_SW_ID.long_name = self.lev2_SW_name
            lev2_SW_ID.units = self.lev2_SW_units
            lev2_SW_ID[:,:,:] = self.col_lev2_SW
            lev2_SW_ID.badval = np.nan
        
            lev2_VR_ID = ncfile.createVariable(self.nexrad_88D_ID+'_VR', np.float32, ('z','y','x'), zlib=True)
            lev2_VR_ID.long_name = self.lev2_VR_name
            lev2_VR_ID.units = self.lev2_VR_units
            lev2_VR_ID[:,:,:] = self.col_lev2_VR
            lev2_VR_ID.badval = np.nan
        
            lev2_RR_ID = ncfile.createVariable(self.nexrad_88D_ID+'_RR', np.float32, ('z','y','x'), zlib=True)
            lev2_RR_ID.long_name = self.lev2_RR_name
            lev2_RR_ID.units = self.lev2_R_units
            lev2_RR_ID[:,:,:] = self.col_lev2_RR
            lev2_RR_ID.badval = np.nan
        
            lev2_RP_ID = ncfile.createVariable(self.nexrad_88D_ID+'_RP', np.float32, ('z','y','x'), zlib=True)
            lev2_RP_ID.long_name = self.lev2_RP_name
            lev2_RP_ID.units = self.lev2_R_units
            lev2_RP_ID[:,:,:] = self.col_lev2_RP
            lev2_RP_ID.badval = np.nan
        
            lev2_RC_ID = ncfile.createVariable(self.nexrad_88D_ID+'_RC', np.float32, ('z','y','x'), zlib=True)
            lev2_RC_ID.long_name = self.lev2_RC_name
            lev2_RC_ID.units = self.lev2_R_units
            lev2_RC_ID[:,:,:] = self.col_lev2_RC
            lev2_RC_ID.badval = np.nan
        
            lev2_D0_ID = ncfile.createVariable(self.nexrad_88D_ID+'_D0', np.float32, ('z','y','x'), zlib=True)
            lev2_D0_ID.long_name = self.lev2_D0_name
            lev2_D0_ID.units = self.lev2_D_units
            lev2_D0_ID[:,:,:] = self.col_lev2_D0
            lev2_D0_ID.badval = np.nan
        
            lev2_NW_ID = ncfile.createVariable(self.nexrad_88D_ID+'_NW', np.float32, ('z','y','x'), zlib=True)
            lev2_NW_ID.long_name = self.lev2_NW_name
            lev2_NW_ID.units = self.lev2_NW_units
            lev2_NW_ID[:,:,:] = self.col_lev2_NW
            lev2_NW_ID.badval = np.nan
        
            lev2_FH_ID = ncfile.createVariable(self.nexrad_88D_ID+'_FH', np.float32, ('z','y','x'), zlib=True)
            lev2_FH_ID.long_name = self.lev2_FH_name
            lev2_FH_ID.units = self.lev2_FH_units
            lev2_FH_ID[:,:,:] = self.col_lev2_FH
            lev2_FH_ID.badval = np.nan
        
            lev2_N2_ID = ncfile.createVariable(self.nexrad_88D_ID+'_N2', np.float32, ('z','y','x'), zlib=True)
            lev2_N2_ID.long_name = self.lev2_N2_name
            lev2_N2_ID.units = self.lev2_N2_units
            lev2_N2_ID[:,:,:] = self.col_lev2_N2
            lev2_N2_ID.badval = np.nan
        
            lev2_DM_ID = ncfile.createVariable(self.nexrad_88D_ID+'_DM', np.float32, ('z','y','x'), zlib=True)
            lev2_DM_ID.long_name = self.lev2_DM_name
            lev2_DM_ID.units = self.lev2_D_units
            lev2_DM_ID[:,:,:] = self.col_lev2_DM
            lev2_DM_ID.badval = np.nan
        
            lev2_FZ_ID = ncfile.createVariable(self.nexrad_88D_ID+'_FZ', np.float32, ('z','y','x'), zlib=True)
            lev2_FZ_ID.long_name = self.lev2_FZ_name
            lev2_FZ_ID.units = self.lev2_Z_units
            lev2_FZ_ID[:,:,:] = self.col_lev2_FZ
            lev2_FZ_ID.badval = np.nan
        else:
            lev2_avail.setncattr('lev2_available', 'False')
            lev2_avail = set_nc_radar_missing_attributes(lev2_avail)
        
        #print(lev2_grp)
        #print(lev2_grp.variables['KDOX_CZ'][5,:,:])
        #print(lev2_grp.variables['KDOX_DR'][5,:,:])
        #ncfile.close()
        #print('Dataset is closed!')
        #pdb.set_trace()
    
        #
        # DEFINE VARS & SET VAR ATTS FOR DOW6 DATA:
        # SKIP FOR NOW

        #
        # DEFINE VARS & SET VAR ATTS FOR APU VALUES:
        if self.apu_data_in_column:
            apu_avail.setncattr('apu_available', 'True')
            apu_avail = set_nc_disdro_attributes(apu_avail, self.apu_info)
            if self.apu_rain_file:
                apu_n_drops_ID          = ncfile.createVariable('apu_total_n_drops', np.float32, ('z','y','x','t'), zlib=True)
                apu_n_drops_ID.units    = self.apu_units_n_drops
                apu_n_drops_ID[:,:,:,:] = self.apu_n_drops
                apu_n_drops_ID.badval = np.nan
                
                apu_concentration_ID        = ncfile.createVariable('apu_total_concentration', np.float32, ('z','y','x','t'), zlib=True)
                apu_concentration_ID.units  = self.apu_units_concentration
                apu_concentration_ID[:,:,:,:] = self.apu_concentration
                apu_concentration_ID.badval = np.nan
                
                apu_liqwater_ID         = ncfile.createVariable('apu_liquid_water_content', np.float32, ('z','y','x','t'), zlib=True)
                apu_liqwater_ID.units   = self.apu_units_liqwatercontent
                apu_liqwater_ID[:,:,:,:] = self.apu_liqwater
                apu_liqwater_ID.badval = np.nan
                
                apu_rainrate_ID         = ncfile.createVariable('apu_rain_rate', np.float32, ('z','y','x','t'), zlib=True)
                apu_rainrate_ID.units   = self.apu_units_rainrate
                apu_rainrate_ID[:,:,:,:] = self.apu_rain_rainrate
                apu_rainrate_ID.badval = np.nan
                
                apu_ref_inRayleigh_ID       = ncfile.createVariable('apu_reflectivity_in_Rayleigh', np.float32, ('z','y','x','t'), zlib=True)
                apu_ref_inRayleigh_ID.units = self.apu_units_ref_inRayleigh
                apu_ref_inRayleigh_ID[:,:,:,:] = self.apu_rain_refinRayleigh
                apu_ref_inRayleigh_ID.badval = np.nan
                
                apu_massweightdiam_ID       = ncfile.createVariable('apu_mass_weighted_drop_diameter', np.float32, ('z','y','x','t'), zlib=True)
                apu_massweightdiam_ID.units = self.apu_units_massweight_diam
                apu_massweightdiam_ID[:,:,:,:] = self.apu_massweight_diam
                apu_massweightdiam_ID.badval = np.nan
                
                apu_maximumdiam_ID          = ncfile.createVariable('apu_maximum_drop_diameter', np.float32, ('z','y','x','t'), zlib=True)
                apu_maximumdiam_ID.units = self.apu_units_maximum_diam
                apu_maximumdiam_ID[:,:,:,:] = self.apu_maximum_diam
                apu_maximumdiam_ID.badval = np.nan
            if self.apu_rain_ter_file:
                apu_n_drops_ter_ID          = ncfile.createVariable('apu_total_n_drops_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_n_drops_ter_ID.units    = self.apu_units_n_drops
                apu_n_drops_ter_ID[:,:,:,:] = self.apu_n_drops_ter
                apu_n_drops_ter_ID.badval = np.nan
                
                apu_concentration_ter_ID        = ncfile.createVariable('apu_total_concentration_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_concentration_ter_ID.units  = self.apu_units_concentration
                apu_concentration_ter_ID[:,:,:,:] = self.apu_concentration_ter
                apu_concentration_ter_ID.badval = np.nan
                
                apu_liqwater_ter_ID         = ncfile.createVariable('apu_liquid_water_content_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_liqwater_ter_ID.units   = self.apu_units_liqwatercontent
                apu_liqwater_ter_ID[:,:,:,:] = self.apu_liqwater_ter
                apu_liqwater_ter_ID.badval = np.nan
                
                apu_rainrate_ter_ID         = ncfile.createVariable('apu_rain_rate_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_rainrate_ter_ID.units   = self.apu_units_rainrate
                apu_rainrate_ter_ID[:,:,:,:] = self.apu_rain_rainrate_ter
                apu_rainrate_ter_ID.badval = np.nan
                
                apu_ref_inRayleigh_ter_ID       = ncfile.createVariable('apu_reflectivity_in_Rayleigh_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_ref_inRayleigh_ter_ID.units = self.apu_units_ref_inRayleigh
                apu_ref_inRayleigh_ter_ID[:,:,:,:] = self.apu_rain_refinRayleigh_ter
                apu_ref_inRayleigh_ter_ID.badval = np.nan
                
                apu_massweightdiam_ter_ID       = ncfile.createVariable('apu_mass_weighted_drop_diameter_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_massweightdiam_ter_ID.units = self.apu_units_massweight_diam
                apu_massweightdiam_ter_ID[:,:,:,:] = self.apu_massweight_diam_ter
                apu_massweightdiam_ter_ID.badval = np.nan
                
                apu_maximumdiam_ter_ID           = ncfile.createVariable('apu_maximum_drop_diameter_tfs', np.float32, ('z','y','x','t'), zlib=True)
                apu_maximumdiam_ter_ID.units = self.apu_units_maximum_diam
                apu_maximumdiam_ter_ID[:,:,:,:] = self.apu_maximum_diam_ter
                apu_maximumdiam_ter_ID.badval = np.nan
                
            if self.apu_dsd_file:
                apu_dsd_concen_ID          = ncfile.createVariable('apu_dsd_concentration', np.float32, ('z','y','x','t','n'), zlib=True)
                apu_dsd_concen_ID.units = self.apu_units_bin_drop_concen
                apu_dsd_concen_ID[:,:,:,:] = self.apu_drop_concen
                apu_dsd_concen_ID.badval = np.nan
                
                apu_dsd_bin_diameter_ID = ncfile.createVariable('apu_dsd_bin_median_diameter', np.float32, ('n'), zlib=True)
                apu_dsd_bin_diameter_ID.units = self.units_binDiameter_andWidth
                apu_dsd_bin_diameter_ID[:] = self.apu_dsd_bins_Dmiddle
                apu_dsd_bin_diameter_ID.badval = np.nan
                
                apu_dsd_bin_width_ID = ncfile.createVariable('apu_dsd_bin_width', np.float32, ('n'), zlib=True)
                apu_dsd_bin_width_ID.units = self.units_binDiameter_andWidth
                apu_dsd_bin_width_ID[:] = self.apu_dsd_bins_widths
                apu_dsd_bin_width_ID.badval = np.nan
                
            if self.apu_dsd_ter_file:
                apu_dsd_concen_tfs_ID          = ncfile.createVariable('apu_dsd_concentration_tfs', np.float32, ('z','y','x','t','n'), zlib=True)
                apu_dsd_concen_tfs_ID.units = self.apu_units_bin_drop_concen
                apu_dsd_concen_tfs_ID[:,:,:,:] = self.apu_drop_concen_ter
                apu_dsd_concen_tfs_ID.badval = np.nan
                
                apu_dsd_bin_diameter_tfs_ID = ncfile.createVariable('apu_dsd_bin_median_diameter_tfs', np.float32, ('n'), zlib=True)
                apu_dsd_bin_diameter_tfs_ID.units = self.units_binDiameter_andWidth
                apu_dsd_bin_diameter_tfs_ID[:] = self.apu_dsd_bins_Dmiddle
                apu_dsd_bin_diameter_tfs_ID.badval = np.nan
                
                apu_dsd_bin_width_tfs_ID = ncfile.createVariable('apu_dsd_bin_width_tfs', np.float32, ('n'), zlib=True)
                apu_dsd_bin_width_tfs_ID.units = self.units_binDiameter_andWidth
                apu_dsd_bin_width_tfs_ID[:] = self.apu_dsd_bins_widths
                apu_dsd_bin_width_tfs_ID.badval = np.nan

            if self.apu_ref_file:
                apu_refS_ID          = ncfile.createVariable('apu_reflectivity_S', np.float32, ('z','y','x','t'), zlib=True)
                apu_refS_ID.units    = self.apu_units_all_refs
                apu_refS_ID[:,:,:,:] = self.apu_ref_atS
                apu_refS_ID.badval = np.nan

                apu_refC_ID          = ncfile.createVariable('apu_reflectivity_C', np.float32, ('z','y','x','t'), zlib=True)
                apu_refC_ID.units    = self.apu_units_all_refs
                apu_refC_ID[:,:,:,:] = self.apu_ref_atC
                apu_refC_ID.badval = np.nan

                apu_refX_ID          = ncfile.createVariable('apu_reflectivity_X', np.float32, ('z','y','x','t'), zlib=True)
                apu_refX_ID.units    = self.apu_units_all_refs
                apu_refX_ID[:,:,:,:] = self.apu_ref_atX
                apu_refX_ID.badval = np.nan
                
                apu_refKu_ID          = ncfile.createVariable('apu_reflectivity_Ku', np.float32, ('z','y','x','t'), zlib=True)
                apu_refKu_ID.units    = self.apu_units_all_refs
                apu_refKu_ID[:,:,:,:] = self.apu_ref_atKu
                apu_refKu_ID.badval = np.nan

                apu_refK_ID          = ncfile.createVariable('apu_reflectivity_K', np.float32, ('z','y','x','t'), zlib=True)
                apu_refK_ID.units    = self.apu_units_all_refs
                apu_refK_ID[:,:,:,:] = self.apu_ref_atK
                apu_refK_ID.badval = np.nan
                
                apu_refKa_ID          = ncfile.createVariable('apu_reflectivity_Ka', np.float32, ('z','y','x','t'), zlib=True)
                apu_refKa_ID.units    = self.apu_units_all_refs
                apu_refKa_ID[:,:,:,:] = self.apu_ref_atKa
                apu_refKa_ID.badval = np.nan
                
                apu_refW_ID          = ncfile.createVariable('apu_reflectivity_W', np.float32, ('z','y','x','t'), zlib=True)
                apu_refW_ID.units    = self.apu_units_all_refs
                apu_refW_ID[:,:,:,:] = self.apu_ref_atW
                apu_refW_ID.badval = np.nan
            if self.apu_atten_file:
                apu_attenS_ID          = ncfile.createVariable('apu_attenuation_S', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenS_ID.units    = self.apu_units_all_attens
                apu_attenS_ID[:,:,:,:] = self.apu_atten_atS
                apu_attenS_ID.badval = np.nan
                
                apu_attenC_ID          = ncfile.createVariable('apu_attenuation_C', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenC_ID.units    = self.apu_units_all_attens
                apu_attenC_ID[:,:,:,:] = self.apu_atten_atC
                apu_attenC_ID.badval = np.nan
                
                apu_attenX_ID          = ncfile.createVariable('apu_attenuation_X', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenX_ID.units    = self.apu_units_all_attens
                apu_attenX_ID[:,:,:,:] = self.apu_atten_atX
                apu_attenX_ID.badval = np.nan
                
                apu_attenKu_ID          = ncfile.createVariable('apu_attenuation_Ku', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenKu_ID.units    = self.apu_units_all_attens
                apu_attenKu_ID[:,:,:,:] = self.apu_atten_atKu
                apu_attenKu_ID.badval = np.nan
                
                apu_attenK_ID          = ncfile.createVariable('apu_attenuation_K', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenK_ID.units    = self.apu_units_all_attens
                apu_attenK_ID[:,:,:,:] = self.apu_atten_atK
                apu_attenK_ID.badval = np.nan

                apu_attenKa_ID          = ncfile.createVariable('apu_attenuation_Ka', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenKa_ID.units    = self.apu_units_all_attens
                apu_attenKa_ID[:,:,:,:] = self.apu_atten_atKa
                apu_attenKa_ID.badval = np.nan
                
                apu_attenW_ID          = ncfile.createVariable('apu_attenuation_W', np.float32, ('z','y','x','t'), zlib=True)
                apu_attenW_ID.units    = self.apu_units_all_attens
                apu_attenW_ID[:,:,:,:] = self.apu_atten_atW
                apu_attenW_ID.badval = np.nan
        else:
            apu_avail.setncattr('apu_available', 'False')
            apu_avail = set_nc_disdro_missing_attributes(apu_avail)
            
            #set platform to missing
            #apu_available_ID = apu_grp.createVariable('apu_avail')
            #apu_available[:] = 'F'

        #print(apu_grp)
        #print(apu_grp.variables['apu_total_n_drops_tfs'])
        #print(apu_grp.variables['apu_total_n_drops_tfs'][32,:,:,0])
        #print(apu_grp.variables['apu_rain_rate_tfs'])
        #print(apu_grp.variables['apu_rain_rate_tfs'][32,:,:,0])

        #ncfile.close()
        #print('Dataset is closed!')
        #pdb.set_trace()
        
        #
        # DEFINE VARS & SET VAR ATTS FOR 2DVD DATA:
        if self.twodvd_data_in_column:
            twodvd_avail.setncattr('twodvd_available', 'True')
            twodvd_avail = set_nc_disdro_attributes(twodvd_avail, self.twodvd_info)
            if self.twodvd_dropcount_file:
                twodvd_dropcount_ID          = ncfile.createVariable('twodvd_dropcounts', np.float32, ('z','y','x','t','n'), zlib=True)
                twodvd_dropcount_ID.units    = self.twodvd_units_ndrops
                twodvd_dropcount_ID[:,:,:,:,:] = self.twodvd_dropcount
                twodvd_dropcount_ID.badval = np.nan
                
            if self.twodvd_dsd_file:
                twodvd_concen_ID          = ncfile.createVariable('twodvd_concen', np.float32, ('z','y','x','t','n'), zlib=True)
                twodvd_concen_ID.units    = self.twodvd_units_concen
                twodvd_concen_ID[:,:,:,:,:] = self.twodvd_dsd_con
                twodvd_concen_ID.badval = np.nan
                
            if self.twodvd_dsd_100_file:
                twodvd100_concen_ID          = ncfile.createVariable('twodvd100_concen', np.float32, ('z','y','x','t','n'), zlib=True)
                twodvd100_concen_ID.units    = self.twodvd_units_concen
                twodvd100_concen_ID[:,:,:,:,:] = self.twodvd_dsd_100_con
                twodvd100_concen_ID.badval = np.nan
                
            if self.twodvd_dsd_ter_file:
                twodvd_concen_ter_ID          = ncfile.createVariable('twodvd_concen_ter', np.float32, ('z','y','x','t','n'), zlib=True)
                twodvd_concen_ter_ID.units    = self.twodvd_units_concen
                twodvd_concen_ter_ID[:,:,:,:,:] = self.twodvd_dsd_ter_con
                twodvd_concen_ter_ID.badval = np.nan
                
            if self.twodvd_dsd_100_ter_file:
                twodvd100_concen_ter_ID          = ncfile.createVariable('twodvd100_concen_ter', np.float32, ('z','y','x','t','n'), zlib=True)
                twodvd100_concen_ter_ID.units    = self.twodvd_units_concen
                twodvd100_concen_ter_ID[:,:,:,:] = self.twodvd_dsd_100_ter_con
                twodvd100_concen_ter_ID.badval = np.nan
                
            if self.twodvd_param_file:
                twodvd_param_tot_drops_ID          = ncfile.createVariable('twodvd_param_tot_ndrops', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_tot_drops_ID.units    = self.twodvd_units_ndrops
                twodvd_param_tot_drops_ID[:,:,:,:] = self.twodvd_param_tot_ndrops
                twodvd_param_tot_drops_ID.badval = np.nan
                
                twodvd_param_tot_concen_ID          = ncfile.createVariable('twodvd_param_tot_concen', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_tot_concen_ID.units    = self.twodvd_units_concen
                twodvd_param_tot_concen_ID[:,:,:,:] = self.twodvd_param_tot_con
                twodvd_param_tot_concen_ID.badval = np.nan
                
                twodvd_param_LWC_ID          = ncfile.createVariable('twodvd_param_LWC', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_LWC_ID.units    = self.twodvd_units_LWC
                twodvd_param_LWC_ID[:,:,:,:] = self.twodvd_param_LWC
                twodvd_param_LWC_ID.badval = np.nan
                
                twodvd_param_rainrate_ID          = ncfile.createVariable('twodvd_param_rainrate', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_rainrate_ID.units    = self.twodvd_units_rainrate
                twodvd_param_rainrate_ID[:,:,:,:] = self.twodvd_param_rainrate
                twodvd_param_rainrate_ID.badval = np.nan
                
                twodvd_param_refRayleigh_ID          = ncfile.createVariable('twodvd_param_refRayleigh', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_refRayleigh_ID.units    = self.twodvd_units_refRayleigh
                twodvd_param_refRayleigh_ID[:,:,:,:] = self.twodvd_param_refRayleigh
                twodvd_param_refRayleigh_ID.badval = np.nan
                
                twodvd_param_Dm_ID          = ncfile.createVariable('twodvd_param_Dm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_Dm_ID.units    = self.twodvd_units_Dm
                twodvd_param_Dm_ID[:,:,:,:] = self.twodvd_param_Dm
                twodvd_param_Dm_ID.badval = np.nan
                
                twodvd_param_Dmax_ID          = ncfile.createVariable('twodvd_param_Dmax', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_Dmax_ID.units    = self.twodvd_units_Dmax
                twodvd_param_Dmax_ID[:,:,:,:] = self.twodvd_param_Dmax
                twodvd_param_Dmax_ID.badval = np.nan
                
                twodvd_param_Dmin_ID          = ncfile.createVariable('twodvd_param_Dmin', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_Dmin_ID.units    = self.twodvd_units_Dmin
                twodvd_param_Dmin_ID[:,:,:,:] = self.twodvd_param_Dmin
                twodvd_param_Dmin_ID.badval = np.nan
                
                twodvd_param_stdevDm_ID          = ncfile.createVariable('twodvd_param_stdevDm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_stdevDm_ID.units    = self.twodvd_units_stdevDm
                twodvd_param_stdevDm_ID[:,:,:,:] = self.twodvd_param_stdev_Dm
                twodvd_param_stdevDm_ID.badval = np.nan
                
            if self.twodvd_param_100_file:
                twodvd100_param_tot_drops_ID          = ncfile.createVariable('twodvd100_param_tot_ndrops', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_tot_drops_ID.units    = self.twodvd_units_ndrops
                twodvd100_param_tot_drops_ID[:,:,:,:] = self.twodvd_param_100_tot_ndrops
                twodvd100_param_tot_drops_ID.badval = np.nan
                
                twodvd100_param_tot_concen_ID          = ncfile.createVariable('twodvd_param100_tot_concen', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_tot_concen_ID.units    = self.twodvd_units_concen
                twodvd100_param_tot_concen_ID[:,:,:,:] = self.twodvd_param_100_tot_con
                twodvd100_param_tot_concen_ID.badval = np.nan
                
                twodvd100_param_LWC_ID          = ncfile.createVariable('twodvd100_param_LWC', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_LWC_ID.units    = self.twodvd_units_LWC
                twodvd100_param_LWC_ID[:,:,:,:] = self.twodvd_param_100_LWC
                twodvd100_param_LWC_ID.badval = np.nan
                
                twodvd100_param_rainrate_ID          = ncfile.createVariable('twodvd100_param_rainrate', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_rainrate_ID.units    = self.twodvd_units_rainrate
                twodvd100_param_rainrate_ID[:,:,:,:] = self.twodvd_param_100_rainrate
                twodvd100_param_rainrate_ID.badval = np.nan
                
                twodvd100_param_refRayleigh_ID          = ncfile.createVariable('twodvd100_param_refRayleigh', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_refRayleigh_ID.units    = self.twodvd_units_refRayleigh
                twodvd100_param_refRayleigh_ID[:,:,:,:] = self.twodvd_param_100_refRayleigh
                twodvd100_param_refRayleigh_ID.badval = np.nan
                
                twodvd100_param_Dm_ID          = ncfile.createVariable('twodvd100_param_Dm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_Dm_ID.units    = self.twodvd_units_Dm
                twodvd100_param_Dm_ID[:,:,:,:] = self.twodvd_param_100_Dm
                twodvd100_param_Dm_ID.badval = np.nan
                
                twodvd100_param_Dmax_ID          = ncfile.createVariable('twodvd100_param_Dmax', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_Dmax_ID.units    = self.twodvd_units_Dmax
                twodvd100_param_Dmax_ID[:,:,:,:] = self.twodvd_param_100_Dmax
                twodvd100_param_Dmax_ID.badval = np.nan
                
                twodvd100_param_Dmin_ID          = ncfile.createVariable('twodvd100_param_Dmin', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_Dmin_ID.units    = self.twodvd_units_Dmin
                twodvd100_param_Dmin_ID[:,:,:,:] = self.twodvd_param_100_Dmin
                twodvd100_param_Dmin_ID.badval = np.nan
                
                twodvd100_param_stdevDm_ID          = ncfile.createVariable('twodvd100_param_stdevDm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_stdevDm_ID.units    = self.twodvd_units_stdevDm
                twodvd100_param_stdevDm_ID[:,:,:,:] = self.twodvd_param_100_stdev_Dm
                twodvd100_param_stdevDm_ID.badval = np.nan
            if self.twodvd_param_ter_file:
                twodvd_param_ter_tot_drops_ID          = ncfile.createVariable('twodvd_param_ter_tot_ndrops', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_tot_drops_ID.units    = self.twodvd_units_ndrops
                twodvd_param_ter_tot_drops_ID[:,:,:,:] = self.twodvd_param_ter_tot_ndrops
                twodvd_param_ter_tot_drops_ID.badval = np.nan
                
                twodvd_param_ter_tot_concen_ID          = ncfile.createVariable('twodvd_param_ter_tot_concen', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_tot_concen_ID.units    = self.twodvd_units_concen
                twodvd_param_ter_tot_concen_ID[:,:,:,:] = self.twodvd_param_ter_tot_con
                twodvd_param_ter_tot_concen_ID.badval = np.nan
                
                twodvd_param_ter_LWC_ID          = ncfile.createVariable('twodvd_param_ter_LWC', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_LWC_ID.units    = self.twodvd_units_LWC
                twodvd_param_ter_LWC_ID[:,:,:,:] = self.twodvd_param_ter_LWC
                twodvd_param_ter_LWC_ID.badval = np.nan
                
                twodvd_param_ter_rainrate_ID          = ncfile.createVariable('twodvd_param_ter_rainrate', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_rainrate_ID.units    = self.twodvd_units_rainrate
                twodvd_param_ter_rainrate_ID[:,:,:,:] = self.twodvd_param_ter_rainrate
                twodvd_param_ter_rainrate_ID.badval = np.nan
                
                twodvd_param_ter_refRayleigh_ID          = ncfile.createVariable('twodvd_param_ter_refRayleigh', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_refRayleigh_ID.units    = self.twodvd_units_refRayleigh
                twodvd_param_ter_refRayleigh_ID[:,:,:,:] = self.twodvd_param_ter_refRayleigh
                twodvd_param_ter_refRayleigh_ID.badval = np.nan
                
                twodvd_param_ter_Dm_ID          = ncfile.createVariable('twodvd_param_ter_Dm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_Dm_ID.units    = self.twodvd_units_Dm
                twodvd_param_ter_Dm_ID[:,:,:,:] = self.twodvd_param_ter_Dm
                twodvd_param_ter_Dm_ID.badval = np.nan
                
                twodvd_param_ter_Dmax_ID          = ncfile.createVariable('twodvd_param_ter_Dmax', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_Dmax_ID.units    = self.twodvd_units_Dmax
                twodvd_param_ter_Dmax_ID[:,:,:,:] = self.twodvd_param_ter_Dmax
                twodvd_param_ter_Dmax_ID.badval = np.nan
                
                twodvd_param_ter_Dmin_ID          = ncfile.createVariable('twodvd_param_ter_Dmin', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_Dmin_ID.units    = self.twodvd_units_Dmin
                twodvd_param_ter_Dmin_ID[:,:,:,:] = self.twodvd_param_ter_Dmin
                twodvd_param_ter_Dmin_ID.badval = np.nan
                
                twodvd_param_ter_stdevDm_ID          = ncfile.createVariable('twodvd_param_ter_stdevDm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_param_ter_stdevDm_ID.units    = self.twodvd_units_stdevDm
                twodvd_param_ter_stdevDm_ID[:,:,:,:] = self.twodvd_param_ter_stdev_Dm
                twodvd_param_ter_stdevDm_ID.badval = np.nan
            if self.twodvd_param_100_ter_file:
                twodvd100_param_ter_tot_drops_ID          = ncfile.createVariable('twodvd100_param_ter_tot_ndrops', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_tot_drops_ID.units    = self.twodvd_units_ndrops
                twodvd100_param_ter_tot_drops_ID[:,:,:,:] = self.twodvd_param_100_ter_tot_ndrops
                twodvd100_param_ter_tot_drops_ID.badval = np.nan
                
                twodvd100_param_ter_tot_concen_ID          = ncfile.createVariable('twodvd100_param_ter_tot_concen', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_tot_concen_ID.units    = self.twodvd_units_concen
                twodvd100_param_ter_tot_concen_ID[:,:,:,:] = self.twodvd_param_100_ter_tot_con
                twodvd100_param_ter_tot_concen_ID.badval = np.nan
                
                twodvd100_param_ter_LWC_ID          = ncfile.createVariable('twodvd100_param_ter_LWC', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_LWC_ID.units    = self.twodvd_units_LWC
                twodvd100_param_ter_LWC_ID[:,:,:,:] = self.twodvd_param_100_ter_LWC
                twodvd100_param_ter_LWC_ID.badval = np.nan
                
                twodvd100_param_ter_rainrate_ID          = ncfile.createVariable('twodvd100_param_ter_rainrate', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_rainrate_ID.units    = self.twodvd_units_rainrate
                twodvd100_param_ter_rainrate_ID[:,:,:,:] = self.twodvd_param_100_ter_rainrate
                twodvd100_param_ter_rainrate_ID.badval = np.nan
                
                twodvd100_param_ter_refRayleigh_ID          = ncfile.createVariable('twodvd100_param_ter_refRayleigh', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_refRayleigh_ID.units    = self.twodvd_units_refRayleigh
                twodvd100_param_ter_refRayleigh_ID[:,:,:,:] = self.twodvd_param_100_ter_refRayleigh
                twodvd100_param_ter_refRayleigh_ID.badval = np.nan
                
                twodvd100_param_ter_Dm_ID          = ncfile.createVariable('twodvd100_param_ter_Dm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_Dm_ID.units    = self.twodvd_units_Dm
                twodvd100_param_ter_Dm_ID[:,:,:,:] = self.twodvd_param_100_ter_Dm
                twodvd100_param_ter_Dm_ID.badval = np.nan
                
                twodvd100_param_ter_Dmax_ID          = ncfile.createVariable('twodvd100_param_ter_Dmax', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_Dmax_ID.units    = self.twodvd_units_Dmax
                twodvd100_param_ter_Dmax_ID[:,:,:,:] = self.twodvd_param_100_ter_Dmax
                twodvd100_param_ter_Dmax_ID.badval = np.nan
                
                twodvd100_param_ter_Dmin_ID          = ncfile.createVariable('twodvd100_param_ter_Dmin', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_Dmin_ID.units    = self.twodvd_units_Dmin
                twodvd100_param_ter_Dmin_ID[:,:,:,:] = self.twodvd_param_100_ter_Dmin
                twodvd100_param_ter_Dmin_ID.badval = np.nan
                
                twodvd100_param_ter_stdevDm_ID          = ncfile.createVariable('twodvd100_param_ter_stdevDm', np.float32, ('z','y','x','t'), zlib=True)
                twodvd100_param_ter_stdevDm_ID.units    = self.twodvd_units_stdevDm
                twodvd100_param_ter_stdevDm_ID[:,:,:,:] = self.twodvd_param_100_ter_stdev_Dm
                twodvd100_param_ter_stdevDm_ID.badval = np.nan
            if self.twodvd_atten_file:
                twodvd_atten_rainrate_ID          = ncfile.createVariable('twodvd_attenuation_rainrate', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_atten_rainrate_ID.units    = self.twodvd_units_rainrate
                twodvd_atten_rainrate_ID[:,:,:,:] = self.twodvd_atten_rainrate
                twodvd_atten_rainrate_ID.badval = np.nan
            
                twodvd_attenS_ID          = ncfile.createVariable('twodvd_attenuation_S', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenS_ID.units    = self.twodvd_units_atten
                twodvd_attenS_ID[:,:,:,:] = self.twodvd_atten_S
                twodvd_attenS_ID.badval = np.nan
                
                twodvd_attenC_ID          = ncfile.createVariable('twodvd_attenuation_C', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenC_ID.units    = self.twodvd_units_atten
                twodvd_attenC_ID[:,:,:,:] = self.twodvd_atten_C
                twodvd_attenC_ID.badval = np.nan
                
                twodvd_attenX_ID          = ncfile.createVariable('twodvd_attenuation_X', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenX_ID.units    = self.twodvd_units_atten
                twodvd_attenX_ID[:,:,:,:] = self.twodvd_atten_X
                twodvd_attenX_ID.badval = np.nan
                
                twodvd_attenKu_ID          = ncfile.createVariable('twodvd_attenuation_Ku', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenKu_ID.units    = self.twodvd_units_atten
                twodvd_attenKu_ID[:,:,:,:] = self.twodvd_atten_Ku
                twodvd_attenKu_ID.badval = np.nan
                
                twodvd_attenK_ID          = ncfile.createVariable('twodvd_attenuation_K', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenK_ID.units    = self.twodvd_units_atten
                twodvd_attenK_ID[:,:,:,:] = self.twodvd_atten_K
                twodvd_attenK_ID.badval = np.nan
                
                twodvd_attenKa_ID          = ncfile.createVariable('twodvd_attenuation_Ka', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenKa_ID.units    = self.twodvd_units_atten
                twodvd_attenKa_ID[:,:,:,:] = self.twodvd_atten_Ka
                twodvd_attenKa_ID.badval = np.nan
                
                twodvd_attenW_ID          = ncfile.createVariable('twodvd_attenuation_W', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_attenW_ID.units    = self.twodvd_units_atten
                twodvd_attenW_ID[:,:,:,:] = self.twodvd_atten_W
                twodvd_attenW_ID.badval = np.nan
            if self.twodvd_ref_file:
                twodvd_ref_rainrate_ID          = twodvd_grp.createVariable('twodvd_ref_rainrate', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_ref_rainrate_ID.units    = self.twodvd_units_rainrate
                twodvd_ref_rainrate_ID[:,:,:,:] = self.twodvd_ref_rainrate
                twodvd_ref_rainrate_ID.badval = np.nan
            
                twodvd_refS_ID          = ncfile.createVariable('twodvd_ref_S', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refS_ID.units    = self.twodvd_units_ref
                twodvd_refS_ID[:,:,:,:] = self.twodvd_ref_S
                twodvd_refS_ID.badval = np.nan
                
                twodvd_refC_ID          = ncfile.createVariable('twodvd_ref_C', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refC_ID.units    = self.twodvd_units_ref
                twodvd_refC_ID[:,:,:,:] = self.twodvd_ref_C
                twodvd_refC_ID.badval = np.nan
                
                twodvd_refX_ID          = ncfile.createVariable('twodvd_ref_X', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refX_ID.units    = self.twodvd_units_ref
                twodvd_refX_ID[:,:,:,:] = self.twodvd_ref_X
                twodvd_refX_ID.badval = np.nan
                
                twodvd_refKu_ID          = ncfile.createVariable('twodvd_ref_Ku', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refKu_ID.units    = self.twodvd_units_ref
                twodvd_refKu_ID[:,:,:,:] = self.twodvd_ref_Ku
                twodvd_refKu_ID.badval = np.nan
                
                twodvd_refK_ID          = ncfile.createVariable('twodvd_ref_K', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refK_ID.units    = self.twodvd_units_ref
                twodvd_refK_ID[:,:,:,:] = self.twodvd_ref_K
                twodvd_refK_ID.badval = np.nan
                
                twodvd_refKa_ID          = ncfile.createVariable('twodvd_ref_Ka', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refKa_ID.units    = self.twodvd_units_ref
                twodvd_refKa_ID[:,:,:,:] = self.twodvd_ref_Ka
                twodvd_refKa_ID.badval = np.nan
                
                twodvd_refW_ID          = ncfile.createVariable('twodvd_ref_W', np.float32, ('z','y','x','t'), zlib=True)
                twodvd_refW_ID.units    = self.twodvd_units_ref
                twodvd_refW_ID[:,:,:,:] = self.twodvd_ref_W
                twodvd_refW_ID.badval = np.nan
        else:
            twodvd_avail.setncattr('twodvd_available', 'False')
            twodvd_avail = set_nc_disdro_missing_attributes(twodvd_avail)
                
                
        #print(twodvd_grp)
        #print(twodvd_grp.variables['twodvd_param_refRayleigh'])
        #print(twodvd_grp.variables['twodvd_param_refRayleigh'][32,:,:,0])
        #print(twodvd_grp.variables['twodvd_param_Dm'])
        #print(twodvd_grp.variables['twodvd_param_Dm'][32,:,:,0])

        #ncfile.close()
        #print('Dataset is closed!')
        #pdb.set_trace()

        #
        # DEFINE VARS & SET VAR ATTS FOR PLUVIOS DATA:
        # SKIP FOR NOW
    
        #
        # DEFINE VARS & SET VAR ATTS FOR GAUGES DATA:
        if self.gauges_data_in_column:
            #gauges_grp.setncattr('gauges_available', 'True')
            gauges_avail.setncattr('gauges_available', 'True')
            gauges_avail = set_nc_disdro_attributes(gauges_avail, self.gauges_info)
            
            gauges_rainrate_ID          = ncfile.createVariable('gauges_rainrate', np.float32, ('z','y','x','t'), zlib=True)
            gauges_rainrate_ID.units    = self.gauge_units_rainrate
            gauges_rainrate_ID[:,:,:,:] = self.gauge_rainrate
            gauges_rainrate_ID.badval = np.nan
        else:
            #gauges_grp.setncattr('gauges_available', 'False')
            gauges_avail.setncattr('gauges_available', 'False')
            gauges_avail = set_nc_disdro_missing_attributes(gauges_avail)
        #print(gauges_grp)
        #print(gauges_grp.variables['gauges_rainrate'])
        #print(gauges_grp.variables['gauges_rainrate'][35,:,:,0])

        #ncfile.close()
        #print('Dataset is closed!')
        #pdb.set_trace()
        
        #
        # DEFINE VARS & SET VAR ATTS FOR MRR DATA:
        if self.mrr_data_in_column:
            mrr_avail.setncattr('mrr_available', 'True')
            mrr_avail = set_nc_mrr_attributes(mrr_avail, self.mrr_info)
        
            mrr_PIA_ID = ncfile.createVariable('mrr_PIA', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_PIA_ID.long_name = self.mrr_PIA_name
            mrr_PIA_ID.units = self.mrr_units_PIA
            mrr_PIA_ID[:,:,:,:] = self.mrr_PIA
            mrr_PIA_ID.badval = np.nan
            
            mrr_aRef_ID = ncfile.createVariable('mrr_aRef', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_aRef_ID.long_name = self.mrr_aRef_name
            mrr_aRef_ID.units = self.mrr_units_ref
            mrr_aRef_ID[:,:,:,:] = self.mrr_aRef
            mrr_aRef_ID.badval = np.nan
            
            mrr_ref_ID = ncfile.createVariable('mrr_ref', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_ref_ID.long_name = self.mrr_ref_name
            mrr_ref_ID.units = self.mrr_units_ref
            mrr_ref_ID[:,:,:,:] = self.mrr_ref
            mrr_ref_ID.badval = np.nan
            
            mrr_refKu_ID = ncfile.createVariable('mrr_refKu', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_refKu_ID.long_name = self.mrr_refKu_name
            mrr_refKu_ID.units = self.mrr_units_ref
            mrr_refKu_ID[:,:,:,:] = self.mrr_ref_ku_adj
            mrr_refKu_ID.badval = np.nan
            
            mrr_RR_ID = ncfile.createVariable('mrr_RR', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_RR_ID.long_name = self.mrr_RR_name
            mrr_RR_ID.units = self.mrr_units_RR
            mrr_RR_ID[:,:,:,:] = self.mrr_RR
            mrr_RR_ID.badval = np.nan
            
            mrr_LWC_ID = ncfile.createVariable('mrr_LWC', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_LWC_ID.long_name = self.mrr_LWC_name
            mrr_LWC_ID.units = self.mrr_units_LWC
            mrr_LWC_ID[:,:,:,:] = self.mrr_LWC
            mrr_LWC_ID.badval = np.nan
            
            mrr_WVEL_ID = ncfile.createVariable('mrr_WVEL', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_WVEL_ID.long_name = self.mrr_WVEL_name
            mrr_WVEL_ID.units = self.mrr_units_WVEL
            mrr_WVEL_ID[:,:,:,:] = self.mrr_VEL
            mrr_WVEL_ID.badval = np.nan
            
            mrr_disdroDm_ID = ncfile.createVariable('mrr_disdro_Dm', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_disdroDm_ID.long_name = self.mrr_disdro_Dm_name
            mrr_disdroDm_ID.units = self.mrr_units_disdro_Dm
            mrr_disdroDm_ID[:,:,:,:] = self.mrr_disdro_Dm
            mrr_disdroDm_ID.badval = np.nan
            
            mrr_disdroNw_ID = ncfile.createVariable('mrr_disdro_Nw', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_disdroNw_ID.long_name = self.mrr_disdro_Nw_name
            mrr_disdroNw_ID.units = self.mrr_units_disdro_Nw
            mrr_disdroNw_ID[:,:,:,:] = self.mrr_disdro_Nw
            mrr_disdroNw_ID.badval = np.nan
            
            mrr_disdroZ_ID = ncfile.createVariable('mrr_disdro_Z', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_disdroZ_ID.long_name = self.mrr_disdro_Z_name
            mrr_disdroZ_ID.units = self.mrr_units_ref
            mrr_disdroZ_ID[:,:,:,:] = self.mrr_disdro_Z
            mrr_disdroZ_ID.badval = np.nan
            
            mrr_dataQual_ID = ncfile.createVariable('mrr_data_quality', np.float32, ('z', 'y', 'x', 't'), zlib=True)
            mrr_dataQual_ID.long_name = self.mrr_data_quality_name
            mrr_dataQual_ID.units = self.mrr_units_data_quality
            mrr_dataQual_ID[:,:,:,:] = self.mrr_data_quality
            mrr_dataQual_ID.badval = np.nan
        else:
            mrr_avail.setncattr('mrr_available', 'False')
            mrr_avail = set_nc_mrr_missing_attributes(mrr_avail)
        #print(mrr_grp)
        #print(mrr_grp.variables['mrr_ref'])
        #print(mrr_grp.variables['mrr_ref'][31,:,:,2])
        #print(mrr_grp.variables['mrr_RR'])
        #print(mrr_grp.variables['mrr_RR'][31,:,:,2])
        #print(mrr_grp.variables['mrr_disdro_Dm'])
        #print(mrr_grp.variables['mrr_disdro_Dm'][31,:,:,2])
        
        #ncfile.close()
        #print('Dataset is closed!')
        #pdb.set_trace()
        
        #
        # DEFINE VARS & SET VAR ATTS FOR GMI-GPROF VALUES:
        # SKIP FOR NOW

        #
        # DEFINE VARS & SET VAR ATTS FOR GMI-L1C-R VALUES:
        # SKIP FOR NOW
        
        #
        # DEFINE VARS & SET VAR ATTS FOR DPR VALUES:
        if self.dpr_data_in_column:
            dpr_avail.setncattr('dpr_available', 'True')
            dpr_avail = set_nc_satellite_attributes(dpr_avail, self.dpr_info)
            
            dpr_HS_sfcType_ID = ncfile.createVariable('dpr_HS_sfcType', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_sfcType_ID.long_name = self.dpr_name_sfcType
            dpr_HS_sfcType_ID.units = self.dpr_units_sfcType
            dpr_HS_sfcType_ID[:,:,:] = self.col_HS_sfcType
            dpr_HS_sfcType_ID.badval = np.nan
            
            dpr_HS_flagPrecip_ID = ncfile.createVariable('dpr_HS_flagPrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_flagPrecip_ID.long_name = self.dpr_name_flagPrecip
            dpr_HS_flagPrecip_ID.units = self.dpr_units_flagPrecip
            dpr_HS_flagPrecip_ID[:,:,:] = self.col_HS_flagPrecip
            dpr_HS_flagPrecip_ID.badval = np.nan
            
            dpr_HS_htStormTop_ID = ncfile.createVariable('dpr_HS_htStormTop', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_htStormTop_ID.long_name = self.dpr_name_htStormTop
            dpr_HS_htStormTop_ID.units = self.dpr_units_htStormTop
            dpr_HS_htStormTop_ID[:,:,:] = self.col_HS_htStormTop
            dpr_HS_htStormTop_ID.badval = np.nan
            
            dpr_HS_htZeroDeg_ID = ncfile.createVariable('dpr_HS_htZeroDeg', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_htZeroDeg_ID.long_name = self.dpr_name_htZeroDeg
            dpr_HS_htZeroDeg_ID.units = self.dpr_units_htZeroDeg
            dpr_HS_htZeroDeg_ID[:,:,:] = self.col_HS_htZeroDeg
            dpr_HS_htZeroDeg_ID.badval = np.nan
            
            dpr_HS_flagBB_ID = ncfile.createVariable('dpr_HS_flagBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_flagBB_ID.long_name = self.dpr_name_flagBB
            dpr_HS_flagBB_ID.units = self.dpr_units_flagBB
            dpr_HS_flagBB_ID[:,:,:] = self.col_HS_flagBB
            dpr_HS_flagBB_ID.badval = np.nan
            
            dpr_HS_htBB_ID = ncfile.createVariable('dpr_HS_htBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_htBB_ID.long_name = self.dpr_name_htBB
            dpr_HS_htBB_ID.units = self.dpr_units_htBB
            dpr_HS_htBB_ID[:,:,:] = self.col_HS_htBB
            dpr_HS_htBB_ID.badval = np.nan
            
            dpr_HS_widthBB_ID = ncfile.createVariable('dpr_HS_widthBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_widthBB_ID.long_name = self.dpr_name_widthBB
            dpr_HS_widthBB_ID.units = self.dpr_units_widthBB
            dpr_HS_widthBB_ID[:,:,:] = self.col_HS_widthBB
            dpr_HS_widthBB_ID.badval = np.nan
            
            dpr_HS_qualityBB_ID = ncfile.createVariable('dpr_HS_qualityBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_qualityBB_ID.long_name = self.dpr_name_qualityBB
            dpr_HS_qualityBB_ID.units = self.dpr_units_qualityBB
            dpr_HS_qualityBB_ID[:,:,:] = self.col_HS_qualityBB
            dpr_HS_qualityBB_ID.badval = np.nan
            
            dpr_HS_typePrecip_ID = ncfile.createVariable('dpr_HS_typePrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_typePrecip_ID.long_name = self.dpr_name_typePrecip
            dpr_HS_typePrecip_ID.units = self.dpr_units_typePrecip
            dpr_HS_typePrecip_ID[:,:,:] = self.col_HS_typePrecip
            dpr_HS_typePrecip_ID.badval = np.nan
            
            dpr_HS_qualityTypePrecip_ID = ncfile.createVariable('dpr_HS_qualityTypePrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_qualityTypePrecip_ID.long_name = self.dpr_name_qualityTypePrecip
            dpr_HS_qualityTypePrecip_ID.units = self.dpr_units_qualityTypePrecip
            dpr_HS_qualityTypePrecip_ID[:,:,:] = self.col_HS_qualityTypePrecip
            dpr_HS_qualityTypePrecip_ID.badval = np.nan
            
            dpr_HS_PIAfinal_ID = ncfile.createVariable('dpr_HS_PIAfinal', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_PIAfinal_ID.long_name = self.dpr_name_PIAfinal
            dpr_HS_PIAfinal_ID.units = self.dpr_units_PIAfinal
            dpr_HS_PIAfinal_ID[:,:,:] = self.col_HS_PIAfinal
            dpr_HS_PIAfinal_ID.badval = np.nan
            
            dpr_HS_corZFacNearSfc_ID = ncfile.createVariable('dpr_HS_corZFacNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_corZFacNearSfc_ID.long_name = self.dpr_name_corZFacNearSfc
            dpr_HS_corZFacNearSfc_ID.units = self.dpr_units_corZFacNearSfc
            dpr_HS_corZFacNearSfc_ID[:,:,:] = self.col_HS_corZFacNearSfc
            dpr_HS_corZFacNearSfc_ID.badval = np.nan
            
            dpr_HS_precipRateNearSfc_ID = ncfile.createVariable('dpr_HS_precipRateNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_precipRateNearSfc_ID.long_name = self.dpr_name_precipRateNearSfc
            dpr_HS_precipRateNearSfc_ID.units = self.dpr_units_precipRateNearSfc
            dpr_HS_precipRateNearSfc_ID[:,:,:] = self.col_HS_precipRateNearSfc
            dpr_HS_precipRateNearSfc_ID.badval = np.nan
            
            dpr_HS_precipRateAve24_ID = ncfile.createVariable('dpr_HS_precipRateAve24', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_precipRateAve24_ID.long_name = self.dpr_name_precipRateAve24
            dpr_HS_precipRateAve24_ID.units = self.dpr_units_precipRateAve24
            dpr_HS_precipRateAve24_ID[:,:,:] = self.col_HS_precipRateAve24
            dpr_HS_precipRateAve24_ID.badval = np.nan
            
            dpr_HS_phaseNearSfc_ID = ncfile.createVariable('dpr_HS_phaseNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_phaseNearSfc_ID.long_name = self.dpr_name_phaseNearSfc
            dpr_HS_phaseNearSfc_ID.units = self.dpr_units_phaseNearSfc
            dpr_HS_phaseNearSfc_ID[:,:,:] = self.col_HS_phaseNearSfc
            dpr_HS_phaseNearSfc_ID.badval = np.nan
            
            dpr_HS_zFacMeas_ID = ncfile.createVariable('dpr_HS_zFacMeas', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_zFacMeas_ID.long_name = self.dpr_name_zFacMeas
            dpr_HS_zFacMeas_ID.units = self.dpr_units_zFacMeas
            dpr_HS_zFacMeas_ID[:,:,:] = self.col_HS_zFacMeas
            dpr_HS_zFacMeas_ID.badval = np.nan
            
            dpr_HS_attenNoPrecip_ID = ncfile.createVariable('dpr_HS_attenNoPrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_attenNoPrecip_ID.long_name = self.dpr_name_attenNoPrecip
            dpr_HS_attenNoPrecip_ID.units = self.dpr_units_attenNoPrecip
            dpr_HS_attenNoPrecip_ID[:,:,:] = self.col_HS_attenNoPrecip
            dpr_HS_attenNoPrecip_ID.badval = np.nan
            
            dpr_HS_corZFac_ID = ncfile.createVariable('dpr_HS_corZFac', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_corZFac_ID.long_name = self.dpr_name_corZFac
            dpr_HS_corZFac_ID.units = self.dpr_units_corZFac
            dpr_HS_corZFac_ID[:,:,:] = self.col_HS_corZFac
            dpr_HS_corZFac_ID.badval = np.nan
            
            dpr_HS_precipRate_ID = ncfile.createVariable('dpr_HS_precipRate', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_precipRate_ID.long_name = self.dpr_name_precipRate
            dpr_HS_precipRate_ID.units = self.dpr_units_precipRate
            dpr_HS_precipRate_ID[:,:,:] = self.col_HS_precipRate
            dpr_HS_precipRate_ID.badval = np.nan
            
            dpr_HS_DSDphase_ID = ncfile.createVariable('dpr_HS_DSDphase', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_DSDphase_ID.long_name = self.dpr_name_DSDphase
            dpr_HS_DSDphase_ID.units = self.dpr_units_DSDphase
            dpr_HS_DSDphase_ID[:,:,:] = self.col_HS_DSDphase
            dpr_HS_DSDphase_ID.badval = np.nan
            
            dpr_HS_PIA_cloudwater_ID = ncfile.createVariable('dpr_HS_PIA_cloudwater', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_PIA_cloudwater_ID.long_name = self.dpr_name_PIA_cloudwater
            dpr_HS_PIA_cloudwater_ID.units = self.dpr_units_PIA_cloudwater
            dpr_HS_PIA_cloudwater_ID[:,:,:] = self.col_HS_PIA_cloudwater
            dpr_HS_PIA_cloudwater_ID.badval = np.nan
            
            dpr_HS_PIA_cloudice_ID = ncfile.createVariable('dpr_HS_PIA_cloudice', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_PIA_cloudice_ID.long_name = self.dpr_name_PIA_cloudice
            dpr_HS_PIA_cloudice_ID.units = self.dpr_units_PIA_cloudice
            dpr_HS_PIA_cloudice_ID[:,:,:] = self.col_HS_PIA_cloudice
            dpr_HS_PIA_cloudice_ID.badval = np.nan
            
            dpr_HS_PIA_watervapor_ID = ncfile.createVariable('dpr_HS_PIA_watervapor', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_PIA_watervapor_ID.long_name = self.dpr_name_PIA_watervapor
            dpr_HS_PIA_watervapor_ID.units = self.dpr_units_PIA_watervapor
            dpr_HS_PIA_watervapor_ID[:,:,:] = self.col_HS_PIA_watervapor
            dpr_HS_PIA_watervapor_ID.badval = np.nan
            
            dpr_HS_PIA_oxygen_ID = ncfile.createVariable('dpr_HS_PIA_oxygen', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_PIA_oxygen_ID.long_name = self.dpr_name_PIA_oxygen
            dpr_HS_PIA_oxygen_ID.units = self.dpr_units_PIA_oxygen
            dpr_HS_PIA_oxygen_ID[:,:,:] = self.col_HS_PIA_oxygen
            dpr_HS_PIA_oxygen_ID.badval = np.nan
            
            dpr_HS_TPW_liquid_ID = ncfile.createVariable('dpr_HS_TPW_liquid', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_TPW_liquid_ID.long_name = self.dpr_name_TPW_liquid
            dpr_HS_TPW_liquid_ID.units = self.dpr_units_TPW_liquid
            dpr_HS_TPW_liquid_ID[:,:,:] = self.col_HS_TPW_liquid
            dpr_HS_TPW_liquid_ID.badval = np.nan
            
            dpr_HS_TPW_ice_ID = ncfile.createVariable('dpr_HS_TPW_ice', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_TPW_ice_ID.long_name = self.dpr_name_TPW_ice
            dpr_HS_TPW_ice_ID.units = self.dpr_units_TPW_ice
            dpr_HS_TPW_ice_ID[:,:,:] = self.col_HS_TPW_ice
            dpr_HS_TPW_ice_ID.badval = np.nan
            
            dpr_HS_dBNw_ID = ncfile.createVariable('dpr_HS_dBNw', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_dBNw_ID.long_name = self.dpr_name_DSD_dBNw
            dpr_HS_dBNw_ID.units = self.dpr_units_DSD_dBNw
            dpr_HS_dBNw_ID[:,:,:] = self.col_HS_DSD_dBNw
            dpr_HS_dBNw_ID.badval = np.nan
            
            dpr_HS_Dm_ID = ncfile.createVariable('dpr_HS_Dm', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_Dm_ID.long_name = self.dpr_name_DSD_Dm
            dpr_HS_Dm_ID.units = self.dpr_units_DSD_Dm
            dpr_HS_Dm_ID[:,:,:] = self.col_HS_DSD_Dm
            dpr_HS_Dm_ID.badval = np.nan
            
            dpr_HS_effPIA_ID = ncfile.createVariable('dpr_HS_effPIA', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_HS_effPIA_ID.long_name = self.dpr_name_effPIA
            dpr_HS_effPIA_ID.units = self.dpr_units_effPIA
            dpr_HS_effPIA_ID[:,:,:] = self.col_HS_eff_PIA
            dpr_HS_effPIA_ID.badval = np.nan
            
            #######################################################################################################
            ##############         NS          ##################
            #######################################################################################################
            dpr_NS_sfcType_ID = ncfile.createVariable('dpr_NS_sfcType', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_sfcType_ID.long_name = self.dpr_name_sfcType
            dpr_NS_sfcType_ID.units = self.dpr_units_sfcType
            dpr_NS_sfcType_ID[:,:,:] = self.col_NS_sfcType
            dpr_NS_sfcType_ID.badval = np.nan
            
            dpr_NS_flagPrecip_ID = ncfile.createVariable('dpr_NS_flagPrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_flagPrecip_ID.long_name = self.dpr_name_flagPrecip
            dpr_NS_flagPrecip_ID.units = self.dpr_units_flagPrecip
            dpr_NS_flagPrecip_ID[:,:,:] = self.col_NS_flagPrecip
            dpr_NS_flagPrecip_ID.badval = np.nan
            
            dpr_NS_htStormTop_ID = ncfile.createVariable('dpr_NS_htStormTop', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_htStormTop_ID.long_name = self.dpr_name_htStormTop
            dpr_NS_htStormTop_ID.units = self.dpr_units_htStormTop
            dpr_NS_htStormTop_ID[:,:,:] = self.col_NS_htStormTop
            dpr_NS_htStormTop_ID.badval = np.nan
            
            dpr_NS_htZeroDeg_ID = ncfile.createVariable('dpr_NS_htZeroDeg', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_htZeroDeg_ID.long_name = self.dpr_name_htZeroDeg
            dpr_NS_htZeroDeg_ID.units = self.dpr_units_htZeroDeg
            dpr_NS_htZeroDeg_ID[:,:,:] = self.col_NS_htZeroDeg
            dpr_NS_htZeroDeg_ID.badval = np.nan
            
            dpr_NS_flagBB_ID = ncfile.createVariable('dpr_NS_flagBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_flagBB_ID.long_name = self.dpr_name_flagBB
            dpr_NS_flagBB_ID.units = self.dpr_units_flagBB
            dpr_NS_flagBB_ID[:,:,:] = self.col_NS_flagBB
            dpr_NS_flagBB_ID.badval = np.nan
            
            dpr_NS_htBB_ID = ncfile.createVariable('dpr_NS_htBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_htBB_ID.long_name = self.dpr_name_htBB
            dpr_NS_htBB_ID.units = self.dpr_units_htBB
            dpr_NS_htBB_ID[:,:,:] = self.col_NS_htBB
            dpr_NS_htBB_ID.badval = np.nan
            
            dpr_NS_widthBB_ID = ncfile.createVariable('dpr_NS_widthBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_widthBB_ID.long_name = self.dpr_name_widthBB
            dpr_NS_widthBB_ID.units = self.dpr_units_widthBB
            dpr_NS_widthBB_ID[:,:,:] = self.col_NS_widthBB
            dpr_NS_widthBB_ID.badval = np.nan
            
            dpr_NS_qualityBB_ID = ncfile.createVariable('dpr_NS_qualityBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_qualityBB_ID.long_name = self.dpr_name_qualityBB
            dpr_NS_qualityBB_ID.units = self.dpr_units_qualityBB
            dpr_NS_qualityBB_ID[:,:,:] = self.col_NS_qualityBB
            dpr_NS_qualityBB_ID.badval = np.nan
            
            dpr_NS_typePrecip_ID = ncfile.createVariable('dpr_NS_typePrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_typePrecip_ID.long_name = self.dpr_name_typePrecip
            dpr_NS_typePrecip_ID.units = self.dpr_units_typePrecip
            dpr_NS_typePrecip_ID[:,:,:] = self.col_NS_typePrecip
            dpr_NS_typePrecip_ID.badval = np.nan
            
            dpr_NS_qualityTypePrecip_ID = ncfile.createVariable('dpr_NS_qualityTypePrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_qualityTypePrecip_ID.long_name = self.dpr_name_qualityTypePrecip
            dpr_NS_qualityTypePrecip_ID.units = self.dpr_units_qualityTypePrecip
            dpr_NS_qualityTypePrecip_ID[:,:,:] = self.col_NS_qualityTypePrecip
            dpr_NS_qualityTypePrecip_ID.badval = np.nan
            
            dpr_NS_PIAfinal_ID = ncfile.createVariable('dpr_NS_PIAfinal', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_PIAfinal_ID.long_name = self.dpr_name_PIAfinal
            dpr_NS_PIAfinal_ID.units = self.dpr_units_PIAfinal
            dpr_NS_PIAfinal_ID[:,:,:] = self.col_NS_PIAfinal
            dpr_NS_PIAfinal_ID.badval = np.nan
            
            dpr_NS_corZFacNearSfc_ID = ncfile.createVariable('dpr_NS_corZFacNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_corZFacNearSfc_ID.long_name = self.dpr_name_corZFacNearSfc
            dpr_NS_corZFacNearSfc_ID.units = self.dpr_units_corZFacNearSfc
            dpr_NS_corZFacNearSfc_ID[:,:,:] = self.col_NS_corZFacNearSfc
            dpr_NS_corZFacNearSfc_ID.badval = np.nan
            
            dpr_NS_precipRateNearSfc_ID = ncfile.createVariable('dpr_NS_precipRateNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_precipRateNearSfc_ID.long_name = self.dpr_name_precipRateNearSfc
            dpr_NS_precipRateNearSfc_ID.units = self.dpr_units_precipRateNearSfc
            dpr_NS_precipRateNearSfc_ID[:,:,:] = self.col_NS_precipRateNearSfc
            dpr_NS_precipRateNearSfc_ID.badval = np.nan
            
            dpr_NS_precipRateAve24_ID = ncfile.createVariable('dpr_NS_precipRateAve24', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_precipRateAve24_ID.long_name = self.dpr_name_precipRateAve24
            dpr_NS_precipRateAve24_ID.units = self.dpr_units_precipRateAve24
            dpr_NS_precipRateAve24_ID[:,:,:] = self.col_NS_precipRateAve24
            dpr_NS_precipRateAve24_ID.badval = np.nan
            
            dpr_NS_phaseNearSfc_ID = ncfile.createVariable('dpr_NS_phaseNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_phaseNearSfc_ID.long_name = self.dpr_name_phaseNearSfc
            dpr_NS_phaseNearSfc_ID.units = self.dpr_units_phaseNearSfc
            dpr_NS_phaseNearSfc_ID[:,:,:] = self.col_NS_phaseNearSfc
            dpr_NS_phaseNearSfc_ID.badval = np.nan
            
            dpr_NS_zFacMeas_ID = ncfile.createVariable('dpr_NS_zFacMeas', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_zFacMeas_ID.long_name = self.dpr_name_zFacMeas
            dpr_NS_zFacMeas_ID.units = self.dpr_units_zFacMeas
            dpr_NS_zFacMeas_ID[:,:,:] = self.col_NS_zFacMeas
            dpr_NS_zFacMeas_ID.badval = np.nan
            
            dpr_NS_attenNoPrecip_ID = ncfile.createVariable('dpr_NS_attenNoPrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_attenNoPrecip_ID.long_name = self.dpr_name_attenNoPrecip
            dpr_NS_attenNoPrecip_ID.units = self.dpr_units_attenNoPrecip
            dpr_NS_attenNoPrecip_ID[:,:,:] = self.col_NS_attenNoPrecip
            dpr_NS_attenNoPrecip_ID.badval = np.nan
            
            dpr_NS_corZFac_ID = ncfile.createVariable('dpr_NS_corZFac', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_corZFac_ID.long_name = self.dpr_name_corZFac
            dpr_NS_corZFac_ID.units = self.dpr_units_corZFac
            dpr_NS_corZFac_ID[:,:,:] = self.col_NS_corZFac
            dpr_NS_corZFac_ID.badval = np.nan
            
            dpr_NS_precipRate_ID = ncfile.createVariable('dpr_NS_precipRate', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_precipRate_ID.long_name = self.dpr_name_precipRate
            dpr_NS_precipRate_ID.units = self.dpr_units_precipRate
            dpr_NS_precipRate_ID[:,:,:] = self.col_NS_precipRate
            dpr_NS_precipRate_ID.badval = np.nan
            
            dpr_NS_DSDphase_ID = ncfile.createVariable('dpr_NS_DSDphase', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_DSDphase_ID.long_name = self.dpr_name_DSDphase
            dpr_NS_DSDphase_ID.units = self.dpr_units_DSDphase
            dpr_NS_DSDphase_ID[:,:,:] = self.col_NS_DSDphase
            dpr_NS_DSDphase_ID.badval = np.nan
            
            dpr_NS_PIA_cloudwater_ID = ncfile.createVariable('dpr_NS_PIA_cloudwater', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_PIA_cloudwater_ID.long_name = self.dpr_name_PIA_cloudwater
            dpr_NS_PIA_cloudwater_ID.units = self.dpr_units_PIA_cloudwater
            dpr_NS_PIA_cloudwater_ID[:,:,:] = self.col_NS_PIA_cloudwater
            dpr_NS_PIA_cloudwater_ID.badval = np.nan
            
            dpr_NS_PIA_cloudice_ID = ncfile.createVariable('dpr_NS_PIA_cloudice', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_PIA_cloudice_ID.long_name = self.dpr_name_PIA_cloudice
            dpr_NS_PIA_cloudice_ID.units = self.dpr_units_PIA_cloudice
            dpr_NS_PIA_cloudice_ID[:,:,:] = self.col_NS_PIA_cloudice
            dpr_NS_PIA_cloudice_ID.badval = np.nan
            
            dpr_NS_PIA_watervapor_ID = ncfile.createVariable('dpr_NS_PIA_watervapor', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_PIA_watervapor_ID.long_name = self.dpr_name_PIA_watervapor
            dpr_NS_PIA_watervapor_ID.units = self.dpr_units_PIA_watervapor
            dpr_NS_PIA_watervapor_ID[:,:,:] = self.col_NS_PIA_watervapor
            dpr_NS_PIA_watervapor_ID.badval = np.nan
            
            dpr_NS_PIA_oxygen_ID = ncfile.createVariable('dpr_NS_PIA_oxygen', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_PIA_oxygen_ID.long_name = self.dpr_name_PIA_oxygen
            dpr_NS_PIA_oxygen_ID.units = self.dpr_units_PIA_oxygen
            dpr_NS_PIA_oxygen_ID[:,:,:] = self.col_NS_PIA_oxygen
            dpr_NS_PIA_oxygen_ID.badval = np.nan
            
            dpr_NS_TPW_liquid_ID = ncfile.createVariable('dpr_NS_TPW_liquid', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_TPW_liquid_ID.long_name = self.dpr_name_TPW_liquid
            dpr_NS_TPW_liquid_ID.units = self.dpr_units_TPW_liquid
            dpr_NS_TPW_liquid_ID[:,:,:] = self.col_NS_TPW_liquid
            dpr_NS_TPW_liquid_ID.badval = np.nan
            
            dpr_NS_TPW_ice_ID = ncfile.createVariable('dpr_NS_TPW_ice', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_TPW_ice_ID.long_name = self.dpr_name_TPW_ice
            dpr_NS_TPW_ice_ID.units = self.dpr_units_TPW_ice
            dpr_NS_TPW_ice_ID[:,:,:] = self.col_NS_TPW_ice
            dpr_NS_TPW_ice_ID.badval = np.nan
            
            dpr_NS_dBNw_ID = ncfile.createVariable('dpr_NS_dBNw', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_dBNw_ID.long_name = self.dpr_name_DSD_dBNw
            dpr_NS_dBNw_ID.units = self.dpr_units_DSD_dBNw
            dpr_NS_dBNw_ID[:,:,:] = self.col_NS_DSD_dBNw
            dpr_NS_dBNw_ID.badval = np.nan
            
            dpr_NS_Dm_ID = ncfile.createVariable('dpr_NS_Dm', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_Dm_ID.long_name = self.dpr_name_DSD_Dm
            dpr_NS_Dm_ID.units = self.dpr_units_DSD_Dm
            dpr_NS_Dm_ID[:,:,:] = self.col_NS_DSD_Dm
            dpr_NS_Dm_ID.badval = np.nan
            
            dpr_NS_effPIA_ID = ncfile.createVariable('dpr_NS_effPIA', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_NS_effPIA_ID.long_name = self.dpr_name_effPIA
            dpr_NS_effPIA_ID.units = self.dpr_units_effPIA
            dpr_NS_effPIA_ID[:,:,:] = self.col_NS_eff_PIA
            dpr_NS_effPIA_ID.badval = np.nan
            
            #######################################################################################################
            ##############         MS          ##################
            #######################################################################################################
            dpr_MS_sfcType_ID = ncfile.createVariable('dpr_MS_sfcType', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_sfcType_ID.long_name = self.dpr_name_sfcType
            dpr_MS_sfcType_ID.units = self.dpr_units_sfcType
            dpr_MS_sfcType_ID[:,:,:] = self.col_MS_sfcType
            dpr_MS_sfcType_ID.badval = np.nan
            
            dpr_MS_flagPrecip_ID = ncfile.createVariable('dpr_MS_flagPrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_flagPrecip_ID.long_name = self.dpr_name_flagPrecip
            dpr_MS_flagPrecip_ID.units = self.dpr_units_flagPrecip
            dpr_MS_flagPrecip_ID[:,:,:] = self.col_MS_flagPrecip
            dpr_MS_flagPrecip_ID.badval = np.nan
            
            dpr_MS_htStormTop_ID = ncfile.createVariable('dpr_MS_htStormTop', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_htStormTop_ID.long_name = self.dpr_name_htStormTop
            dpr_MS_htStormTop_ID.units = self.dpr_units_htStormTop
            dpr_MS_htStormTop_ID[:,:,:] = self.col_MS_htStormTop
            dpr_MS_htStormTop_ID.badval = np.nan
            
            dpr_MS_htZeroDeg_ID = ncfile.createVariable('dpr_MS_htZeroDeg', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_htZeroDeg_ID.long_name = self.dpr_name_htZeroDeg
            dpr_MS_htZeroDeg_ID.units = self.dpr_units_htZeroDeg
            dpr_MS_htZeroDeg_ID[:,:,:] = self.col_MS_htZeroDeg
            dpr_MS_htZeroDeg_ID.badval = np.nan
            
            dpr_MS_flagBB_ID = ncfile.createVariable('dpr_MS_flagBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_flagBB_ID.long_name = self.dpr_name_flagBB
            dpr_MS_flagBB_ID.units = self.dpr_units_flagBB
            dpr_MS_flagBB_ID[:,:,:] = self.col_MS_flagBB
            dpr_MS_flagBB_ID.badval = np.nan
            
            dpr_MS_htBB_ID = ncfile.createVariable('dpr_MS_htBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_htBB_ID.long_name = self.dpr_name_htBB
            dpr_MS_htBB_ID.units = self.dpr_units_htBB
            dpr_MS_htBB_ID[:,:,:] = self.col_MS_htBB
            dpr_MS_htBB_ID.badval = np.nan
            
            dpr_MS_widthBB_ID = ncfile.createVariable('dpr_MS_widthBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_widthBB_ID.long_name = self.dpr_name_widthBB
            dpr_MS_widthBB_ID.units = self.dpr_units_widthBB
            dpr_MS_widthBB_ID[:,:,:] = self.col_MS_widthBB
            dpr_MS_widthBB_ID.badval = np.nan
            
            dpr_MS_qualityBB_ID = ncfile.createVariable('dpr_MS_qualityBB', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_qualityBB_ID.long_name = self.dpr_name_qualityBB
            dpr_MS_qualityBB_ID.units = self.dpr_units_qualityBB
            dpr_MS_qualityBB_ID[:,:,:] = self.col_MS_qualityBB
            dpr_MS_qualityBB_ID.badval = np.nan
            
            dpr_MS_typePrecip_ID = ncfile.createVariable('dpr_MS_typePrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_typePrecip_ID.long_name = self.dpr_name_typePrecip
            dpr_MS_typePrecip_ID.units = self.dpr_units_typePrecip
            dpr_MS_typePrecip_ID[:,:,:] = self.col_MS_typePrecip
            dpr_MS_typePrecip_ID.badval = np.nan
            
            dpr_MS_qualityTypePrecip_ID = ncfile.createVariable('dpr_MS_qualityTypePrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_qualityTypePrecip_ID.long_name = self.dpr_name_qualityTypePrecip
            dpr_MS_qualityTypePrecip_ID.units = self.dpr_units_qualityTypePrecip
            dpr_MS_qualityTypePrecip_ID[:,:,:] = self.col_MS_qualityTypePrecip
            dpr_MS_qualityTypePrecip_ID.badval = np.nan
            
            dpr_MS_PIAfinal_ID = ncfile.createVariable('dpr_MS_PIAfinal', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_PIAfinal_ID.long_name = self.dpr_name_PIAfinal
            dpr_MS_PIAfinal_ID.units = self.dpr_units_PIAfinal
            dpr_MS_PIAfinal_ID[:,:,:] = self.col_MS_PIAfinal
            dpr_MS_PIAfinal_ID.badval = np.nan
            
            dpr_MS_corZFacNearSfc_ID = ncfile.createVariable('dpr_MS_corZFacNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_corZFacNearSfc_ID.long_name = self.dpr_name_corZFacNearSfc
            dpr_MS_corZFacNearSfc_ID.units = self.dpr_units_corZFacNearSfc
            dpr_MS_corZFacNearSfc_ID[:,:,:] = self.col_MS_corZFacNearSfc
            dpr_MS_corZFacNearSfc_ID.badval = np.nan
            
            dpr_MS_precipRateNearSfc_ID = ncfile.createVariable('dpr_MS_precipRateNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_precipRateNearSfc_ID.long_name = self.dpr_name_precipRateNearSfc
            dpr_MS_precipRateNearSfc_ID.units = self.dpr_units_precipRateNearSfc
            dpr_MS_precipRateNearSfc_ID[:,:,:] = self.col_MS_precipRateNearSfc
            dpr_MS_precipRateNearSfc_ID.badval = np.nan
            
            dpr_MS_precipRateAve24_ID = ncfile.createVariable('dpr_MS_precipRateAve24', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_precipRateAve24_ID.long_name = self.dpr_name_precipRateAve24
            dpr_MS_precipRateAve24_ID.units = self.dpr_units_precipRateAve24
            dpr_MS_precipRateAve24_ID[:,:,:] = self.col_MS_precipRateAve24
            dpr_MS_precipRateAve24_ID.badval = np.nan
            
            dpr_MS_phaseNearSfc_ID = ncfile.createVariable('dpr_MS_phaseNearSfc', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_phaseNearSfc_ID.long_name = self.dpr_name_phaseNearSfc
            dpr_MS_phaseNearSfc_ID.units = self.dpr_units_phaseNearSfc
            dpr_MS_phaseNearSfc_ID[:,:,:] = self.col_MS_phaseNearSfc
            dpr_MS_phaseNearSfc_ID.badval = np.nan
            
            dpr_MS_zFacMeas_ID = ncfile.createVariable('dpr_MS_zFacMeas', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_zFacMeas_ID.long_name = self.dpr_name_zFacMeas
            dpr_MS_zFacMeas_ID.units = self.dpr_units_zFacMeas
            dpr_MS_zFacMeas_ID[:,:,:] = self.col_MS_zFacMeas
            dpr_MS_zFacMeas_ID.badval = np.nan
            
            dpr_MS_attenNoPrecip_ID = ncfile.createVariable('dpr_MS_attenNoPrecip', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_attenNoPrecip_ID.long_name = self.dpr_name_attenNoPrecip
            dpr_MS_attenNoPrecip_ID.units = self.dpr_units_attenNoPrecip
            dpr_MS_attenNoPrecip_ID[:,:,:] = self.col_MS_attenNoPrecip
            dpr_MS_attenNoPrecip_ID.badval = np.nan
            
            dpr_MS_corZFac_ID = ncfile.createVariable('dpr_MS_corZFac', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_corZFac_ID.long_name = self.dpr_name_corZFac
            dpr_MS_corZFac_ID.units = self.dpr_units_corZFac
            dpr_MS_corZFac_ID[:,:,:] = self.col_MS_corZFac
            dpr_MS_corZFac_ID.badval = np.nan
            
            dpr_MS_PIA_cloudwater_ID = ncfile.createVariable('dpr_MS_PIA_cloudwater', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_PIA_cloudwater_ID.long_name = self.dpr_name_PIA_cloudwater
            dpr_MS_PIA_cloudwater_ID.units = self.dpr_units_PIA_cloudwater
            dpr_MS_PIA_cloudwater_ID[:,:,:] = self.col_MS_PIA_cloudwater
            dpr_MS_PIA_cloudwater_ID.badval = np.nan
            
            dpr_MS_PIA_cloudice_ID = ncfile.createVariable('dpr_MS_PIA_cloudice', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_PIA_cloudice_ID.long_name = self.dpr_name_PIA_cloudice
            dpr_MS_PIA_cloudice_ID.units = self.dpr_units_PIA_cloudice
            dpr_MS_PIA_cloudice_ID[:,:,:] = self.col_MS_PIA_cloudice
            dpr_MS_PIA_cloudice_ID.badval = np.nan
            
            dpr_MS_PIA_watervapor_ID = ncfile.createVariable('dpr_MS_PIA_watervapor', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_PIA_watervapor_ID.long_name = self.dpr_name_PIA_watervapor
            dpr_MS_PIA_watervapor_ID.units = self.dpr_units_PIA_watervapor
            dpr_MS_PIA_watervapor_ID[:,:,:] = self.col_MS_PIA_watervapor
            dpr_MS_PIA_watervapor_ID.badval = np.nan
            
            dpr_MS_PIA_oxygen_ID = ncfile.createVariable('dpr_MS_PIA_oxygen', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_PIA_oxygen_ID.long_name = self.dpr_name_PIA_oxygen
            dpr_MS_PIA_oxygen_ID.units = self.dpr_units_PIA_oxygen
            dpr_MS_PIA_oxygen_ID[:,:,:] = self.col_MS_PIA_oxygen
            dpr_MS_PIA_oxygen_ID.badval = np.nan
            
            dpr_MS_TPW_liquid_ID = ncfile.createVariable('dpr_MS_TPW_liquid', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_TPW_liquid_ID.long_name = self.dpr_name_TPW_liquid
            dpr_MS_TPW_liquid_ID.units = self.dpr_units_TPW_liquid
            dpr_MS_TPW_liquid_ID[:,:,:] = self.col_MS_TPW_liquid
            dpr_MS_TPW_liquid_ID.badval = np.nan
            
            dpr_MS_TPW_ice_ID = ncfile.createVariable('dpr_MS_TPW_ice', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_TPW_ice_ID.long_name = self.dpr_name_TPW_ice
            dpr_MS_TPW_ice_ID.units = self.dpr_units_TPW_ice
            dpr_MS_TPW_ice_ID[:,:,:] = self.col_MS_TPW_ice
            dpr_MS_TPW_ice_ID.badval = np.nan
            
            dpr_MS_effPIA_ID = ncfile.createVariable('dpr_MS_effPIA', np.float32, ('z', 'y', 'x'), zlib=True)
            dpr_MS_effPIA_ID.long_name = self.dpr_name_effPIA
            dpr_MS_effPIA_ID.units = self.dpr_units_effPIA
            dpr_MS_effPIA_ID[:,:,:] = self.col_MS_eff_PIA
            dpr_MS_effPIA_ID.badval = np.nan
        else:
            dpr_avail.setncattr('dpr_available', 'False')
            dpr_avail = set_nc_satellite_missing_attributes(dpr_avail)


        #print(ncfile)
        #print(dpr_grp)
        #print(dpr_MS_grp)
        #print(mrr_grp.variables['mrr_ref'])
        #print(mrr_grp.variables['mrr_ref'][31,:,:,2])
        #print(mrr_grp.variables['mrr_RR'])
        #print(mrr_grp.variables['mrr_RR'][31,:,:,2])
        #print(mrr_grp.variables['mrr_disdro_Dm'])
        #print(mrr_grp.variables['mrr_disdro_Dm'][31,:,:,2])
        
        #print(ncfile)
        #print('Dataset is closed!')
        #pdb.set_trace()

        #
        #DEFINE VARS & SET VAR ATTS FOR GPM 2BCMB VALUES:
        # SKIP FOR NOW

        #
        # DEFINE VARS & SET VAR ATTS FOR MRMS VALUES:
        # SKIP FOR NOW
 
        #
        # DEFINE VARS & SET VAR ATTS FOR SOUNDINGS DATA:
        # SKIP FOR NOW

        #
        # DEFINE VARS & SET VAR ATTS FOR [other platforms] DATA:
        # do just as did above for any additional platforms...
        #

        # close the newly made .nc file
        ncfile.close()
        print(f'--> column file has been saved {column_nc_fname}')

        # Now that have saved the column .nc file, genearte
        #  a .txt file of the .nc header information:
        header_cmnd = 'ncdump -h '+column_nc_fname+' >& '+header_fname
        os.system(header_cmnd)
        print(f'--> column header .txt file has been saved {header_fname}')

    # ***************************************************************************************
    def set_plat_values(self, plat_name, data_file):
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
        #self.radars = ['NPOL','D3R','DOW6','KABR','KCAE','KEAX','KICT','KLIX',
        #          'KMRX','KTLH','KAKQ','KCCX','KEVX','KILN','KLOT','KMVX','KTLX',
        #          'KAMX','KCLX','KFSD','KILX','KLSX','KNQA','KTWX','KAPX','KCRP', 
        #          'KFTG','KINX','KLTX','KOKX','KTYX','KARX','KDDC','KFWS','KIWX',
        #          'KLZK','KPAH','PAEC','KBMX','KDGX','KGRK','KJAX','KMHX','KRAX',
        #          'PAIH','KBOX','KDLH','KGRR','KJGX','KMKX','KSGF','PGUA','KBRO',
        #          'KDMX','KGSP','KJKL','KMLB','KSHV','PHKI','KBUF','KDOX','KHGX',
        #          'KLCH','KMOB','KSRX','PHMO','KBYX','KDVN','KHTX','KLGX','KMQT',
        #          'KTBW','TJUA']
        #define apus
        self.apus = ['apu01','apu02','apu03','apu04','apu05',
                     'apu06','apu07','apu08','apu09','apu10',
                     'apu11','apu12','apu13','apu14','apu15',
                     'apu16','apu17','apu18','apu20',
                     'apu21','apu23','apu25','apu30']
        
        #define 2DVDs, APU_Pluvios, MRRs
        self.twodvds = ['SN25','SN35','SN36','SN37','SN38','SN70']
        self.apu_pluvio = ['apu04_pluvio','apu10_pluvio','apu30_pluvio']
        self.mrr = ['MRR2-01','MRR2-02','MRR2-03','MRR2-04']
        
        #pdb.set_trace()
        if plat_name in self.radars:
            plat_type = 'radar'
            
            #unzip file (if needed)
            #file_basename = os.path.basename(data_file)
            #if file_basename.endswith('.gz'):
            #    cf_file = sim.ungzip_file(data_file)
            #else:
            #    cf_file = data_file
            radar = pyart.io.read(data_file, file_field_names=True)
        elif plat_name in self.apus:
            plat_type='APU'
        elif plat_name in self.twodvds:
            plat_type='2DVD'
        elif plat_name in self.apu_pluvio:
            plat_type='pluvio'
        elif plat_name in self.mrr:
            plat_type='MRR'
        elif plat_name == 'GMI':
            plat_type='GMI'
        elif plat_name == 'DPR':
            plat_type='DPR'
        elif plat_name == '2BCMB':
            plat_type='2BCMB'
        else:
            sys.exit(f'Platform {plat_name} is not defined in set_plat_values method -- please check')
            
        #For radar platforms:  Get info using pyart metadata:
        if plat_type == 'radar':
            r_lat_decdeg = radar.latitude['data'][0]
            r_lon_decdeg = radar.longitude['data'][0]
            plat_lat_d, plat_lat_m, plat_lat_s, plat_lon_d, plat_lon_m, plat_lon_s = sim.dd2dms(r_lat_decdeg, r_lon_decdeg)
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
            plat_year = radar_datetime.year
            plat_mon  = radar_datetime.month
            plat_day  = radar_datetime.day
            plat_hr   = radar_datetime.hour
            plat_min  = radar_datetime.minute
            plat_sec  = radar_datetime.second
            
            string_year = str(plat_year).zfill(4)
            string_mon  = str(plat_mon).zfill(2)
            string_day  = str(plat_day).zfill(2)
            string_hr   = str(plat_hr).zfill(2)
            string_min  = str(plat_min).zfill(2)
            string_sec  = str(plat_sec).zfill(2)
            plat_timestamp = string_year+string_mon+string_day+'_'+string_hr+string_min+string_sec
        
            #get time offset vs main_plat's timestamp
            #main_plat_timestamp = self.main_plat_info['timestamp'] #'YYYYMMDD_HHMMSS'
            main_year = int(self.main_year)
            main_mon  = int(self.main_mon)
            main_day  = int(self.main_day)
            main_hr   = int(self.main_hr)
            main_min  = int(self.main_min)
            main_sec  = int(self.main_sec)
            off_sec = plat_sec - main_sec
            off_min = plat_min - main_min
            off_hr  = plat_hr  - main_hr	#these bigger intervals should be same,
            off_day = plat_day - main_day	#but including here for completeness
            off_mon = plat_mon - main_mon
            off_year= plat_year - main_year
            off_min = off_min*60        # sec in the min
            off_hr  = off_hr*3600       # sec in the hrs
            off_day = off_day*86400     # sec in the days (should be 0)
            off_mon = off_mon*2.628E6   # sec in the mons (should be 0)
            off_year= off_year*3.154E7  # sec in the yrs (should be 0)
            offset = off_sec + off_min + off_hr + off_day + off_mon + off_year
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
    def get_time_lag(self):
            
        # Define interval_datetime array for the full data time interval:
        #  time interval = main_plat_timestamp +/- halftime_interval # of mins
        self.interval_datetime    = np.empty(2*self.halftime_interval + 1, dtype=object)
        self.n_interval_times    = len(self.interval_datetime)
        self.t_values = np.arange(-self.halftime_interval, self.halftime_interval+1)
                
        # use datetime to get timestamps for ea minute in the interval
        main_datetime = datetime(int(self.main_year),int(self.main_mon),int(self.main_day),int(self.main_hr),int(self.main_min))
        for i, step in enumerate(np.arange(-self.halftime_interval,self.halftime_interval+1)):
            #compute delta datetime from the main datetime
            #have to convert numpy.int64 to python int for timedelta input using .item()
            self.interval_datetime[i] = main_datetime + timedelta(minutes=step.item())
            
        #define the dimensions for ground platform data arrays (apu, twodvd, gauges, mrr, etc)
        #  4-d:  [column z dir  X  column x dir  X  column y dir  X   times in interval ]
        self.ground_plat_dims = (self.z_values.shape[0], self.lon_values.shape[0], self.lat_values.shape[0], self.n_interval_times)

    # ***************************************************************************************
    def setup_main_plat(self):
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
        
        self.radars = ['NPOL','D3R','DOW6','KABR','KCAE','KEAX','KICT','KLIX',
                  'KMRX','KTLH','KAKQ','KCCX','KEVX','KILN','KLOT','KMVX','KTLX',
                  'KAMX','KCLX','KFSD','KILX','KLSX','KNQA','KTWX','KAPX','KCRP', 
                  'KFTG','KINX','KLTX','KOKX','KTYX','KARX','KDDC','KFWS','KIWX',
                  'KLZK','KPAH','PAEC','KBMX','KDGX','KGRK','KJAX','KMHX','KRAX',
                  'PAIH','KBOX','KDLH','KGRR','KJGX','KMKX','KSGF','PGUA','KBRO',
                  'KDMX','KGSP','KJKL','KMLB','KSHV','PHKI','KBUF','KDOX','KHGX',
                  'KLCH','KMOB','KSRX','PHMO','KBYX','KDVN','KHTX','KLGX','KMQT',
                  'KTBW','TJUA']
                  
        if self.main_plat_name not in self.radars:
            sys.exit('------ MAIN PLATFORM ERROR: must be a valid option -----\n'
                    '-- VALID OPTS: NPOL D3R DOW6 or GPM GV VN QCd NEXRADs --')
        
        #if self.main_plat_name == 'NPOL' and file_npol == '':
        #    sys.exit('--- must send NPOL file to setup_main_plat method!! ---')
        #elif self.main_plat_name == 'D3R' and files_d3r_ku == '':
        #    sys.exit('--- must send D3R Ku-band files to setup_main_plat method!! ---')
        #elif self.main_plat_name == 'DOW6' and file_dow6 == '':
        #    sys.exit('--- must send DOW6 file to setup_main_plat method!! ---')
        
        #if self.main_plat_name != 'NPOL' and self.main_plat_name != 'D3R' and self.main_plat_name != 'DOW6' and \
        #    files_88d == '':
        #    sys.exit('--- must send 88D files to setup_main_plat method!! ---')
            
        if self.main_plat_name == 'NPOL':
            self.main_plat_info = self.set_main_plat_values()
        #elif self.main_plat_name == 'D3R':
        #    print('...THIS NEEDS TO BE TESTED...')
        #    pdb.set_trace()
        #    file_D3R_Ku = files_d3r_ku[0]  #use 1st Ku file for getting info
        #    self.main_plat_info = self.set_main_plat_values(file_d3r_ku)
        #elif self.main_plat_name == 'DOW6':
        #    self.main_plat_info = self.set_main_plat_values(file_dow6)
        #else:
        #    #--- NEXRAD/88D RADARS IN GPM GV VN: --------------------------
        #    print('...THIS NEEDS TO BE TESTED...')
        #    pdb.set_trace()
        #    basenames_88D = os.path.basename(files_88D)
        #    siteIDs_88D = basenames_88D[0:4]
        #    file_88D_sub = np.where(siteIDs_88D == self.main_plat)
        #    file_88D = files_88d[file_88D_sub] 
        #    self.main_plat_info = self.set_main_plat_values(file_88D)

    # ***************************************************************************************
    def set_main_plat_values(self):
        '''
        ;------------------------------------------------------------
        ;  This method takes the main_plat's data file and sets up 
        ;  the main_plat timestamp & scan_type strings
        ;------------------------------------------------------------
        ;  result = set_main_plat_values(main_plat_name, main_plat_data_file)
        ;  main_plat_name:	a string identifier for the main platform
        ;           currently, vaild options are:
        ;              'NPOL'  'D3R'  'DOW6'  'KDOX'  'KAKQ'  'KLGX'
        ;   main_plat_data_file: original data file from the main platform
        ;           for NPOL: the QC'd .cf.gz file
        ;           for 88Ds: the QC'd .cf.gz file
        ;           for D3R:  either Ka or Ku .nc file should be OK
        ;           for DOW6: the Cfradial .nc file
        ;   result: dictionary with main_plat details to use in the new .nc file:
        ;        name:       string identifier of the main_plat
        ;        timestamp:  15 char string, as 'YYYYMMDD_HHMMSS'
        ;        lat_deg:    latitude in dec deg of main platform 
        ;        lon_deg:    londitude in dec deg of main platform
        ;        scan_type:  3 char string, 'RHI', 'PPI', 'UNK'
        ;
        ;   Dependencies:
        ;   rsl_d3r_to_radar.pro       add-on to RSL_in_IDL for NASA D3R radar data
        ;   rsl_cfradial_to_radar.pro  add-on to RSL_in_IDL for Cfradial data
        ;   dms2dd.pro                 utility to convert deg-min-sec to decimal deg
        ;
        ;------------------------------------------------------------
        '''
        
        #Set up dictionary of the main paltform's info needed in build_column.pro
        #main_plat_info = {'name':'', 'timestamp':'', 'lat_deg':0.0, 'lon_deg':0.0, 
        #                  'scan_type':''}
                          
        # Platform must be one of the radars
        if self.main_plat_name in self.radars:
            
            #check if file is zipped
            #file_basename = os.path.basename(self.npol_file)
            #if file_basename.endswith('.gz'):
            #    self.npol_file = sim.ungzip_file(self.npol_file)
            radar = pyart.io.read(self.npol_file, file_field_names=True)
            
        # Get info from radar pyart object
        radar_datetime = pyart.util.datetime_from_radar(radar) #this is same as above line -- remove and fix other code below
        self.main_plat_datetime = datetime(radar_datetime.year, radar_datetime.month,radar_datetime.day,radar_datetime.hour,
                                            radar_datetime.minute, radar_datetime.second)
        self.main_year = str(radar_datetime.year).zfill(4)
        self.main_mon  = str(radar_datetime.month).zfill(2)
        self.main_day  = str(radar_datetime.day).zfill(2)
        self.main_hr   = str(radar_datetime.hour).zfill(2)
        self.main_min  = str(radar_datetime.minute).zfill(2)
        self.main_sec  = str(radar_datetime.second).zfill(2)
        #main_sec  = round(main_sec) # dont think this is needed
        
        self.main_lat_deg = radar.latitude['data'][0]
        self.main_lon_deg = radar.longitude['data'][0]
        
        self.main_scan_type = radar.scan_type.upper()
        if self.main_scan_type != 'PPI' and self.main_scan_type != 'RHI': self.main_scan_type = 'UNK'

        self.main_plat_timestamp = self.main_year+self.main_mon+self.main_day+'_'+self.main_hr+self.main_min+self.main_sec

        print(f'main_plat timestamp: {self.main_plat_timestamp}')
        print(f'main_plat scan type: {self.main_scan_type}')

    # ***************************************************************************************
    def define_box(self):
    
        #center_on = kwargs['center_loc']
        #box_spacing = kwargs['column_grid_horiz_spacing']
        #box_limit = kwargs['column_grid_horiz_limit']
        #vert_spacing = kwargs['column_grid_vert_spacing']
        #vert_limit = kwargs['column_grid_vert_limit']
    
        #Set box center location:  degrees, minutes, seconds, + or - sign on degrees
        box_c_lat_d, box_c_lat_m, box_c_lat_s, \
        box_c_lon_d, box_c_lon_m, box_c_lon_s, \
        elv = get_platform_loc(self.center_on)
        
        #Get center location in decimal degrees
        self.cntr_lat_deg, self.cntr_lon_deg = sim.dms2dd(box_c_lat_d, box_c_lat_m, box_c_lat_s, 
                                            box_c_lon_d, box_c_lon_m, box_c_lon_s)
                                            
        #Get # of small boxes in the big box per x/y dir:
        self.n_horiz_grid_boxes = self.box_limit / self.box_spacing
        self.n_vert_boxes = self.vert_limit / self.vert_spacing
        if self.n_horiz_grid_boxes % 1.0 != 0.0:
           sys.exit('---ERROR:  BOX SPACING MUST EVENLY SUBDIVIDE HORIZONTAL BOX LIMIT!---')
        if self.n_vert_boxes % 1.0 != 0.0:
           sys.exit('---ERROR:  VERTICAL SPACING MUST EVENLY SUBDIVIDE VERTICAL LIMIT!---')

        self.n_vert_boxes = int(self.n_vert_boxes)
        self.n_horiz_grid_boxes = int(self.n_horiz_grid_boxes)
        #print(f'  No. of horiz small boxes: {n_grid_boxes}') #testing
        #print(f'  No. of  vert small boxes: {n_vert_boxes}') #testing
        
        #Will have (# of small boxes)+1 values for defining lat/lons:
        #these will be doubles in decimal degrees
        self.lat_values = np.zeros(self.n_horiz_grid_boxes+1)
        self.lon_values = np.zeros(self.n_horiz_grid_boxes+1)
        
        # Make float arrays for x, y, z
        self.have_xy_vals = False
        self.x_values = np.zeros(self.n_horiz_grid_boxes+1)
        self.y_values = np.zeros(self.n_horiz_grid_boxes+1)
        
        #Determine where center point will be, then get lat/lons for each grid point:
        #    odd # of grid boxes: center point is middle point of middle grid box
        #   even # of grid boxed: center point is at 4 corners of middle boxes
        if self.n_horiz_grid_boxes % 2 == 0:
            #even # of boxes, odd # of defining points
            #need to arrange center point at middle corner
            #lat of center point = middle lat of grid
            #lon of center point = middle lon of grid
            mid_sub = int(self.n_horiz_grid_boxes / 2)
            self.x_values[mid_sub] = 0.0
            self.y_values[mid_sub] = 0.0
            
            self.lat_values[mid_sub] = self.cntr_lat_deg
            self.lon_values[mid_sub] = self.cntr_lon_deg
   
            #after center corner, now set lat/lons for remaining edges:
            n_edges_from_center = mid_sub
            for edge_from_center in range(1, n_edges_from_center+1):
                #from middle corner, have ((n_grid_boxes)/2) to each side of big box
                #    call my Function: xy2ll  -->  takes distance in [km]
   
                x_dist = (self.box_spacing/1000.0)*edge_from_center  #[km]
                y_dist = (self.box_spacing/1000.0)*edge_from_center  #[km]
      
                # use for lons to East of center corner:
                to_east_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, x_dist, 0.0)
      
                # use for lons to West of center corner:
                to_west_latlon= sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, -1.0*x_dist, 0.0)
      
                # use for lats to North of center corner:
                to_north_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, 0.0, y_dist)
      
                # use for lons to South of center corner:
                to_south_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, 0.0, -1.0*y_dist)
      
                self.lon_values[mid_sub+edge_from_center] = to_east_latlon[1]
                self.lon_values[mid_sub-edge_from_center] = to_west_latlon[1]
                self.lat_values[mid_sub+edge_from_center] = to_north_latlon[0]
                self.lat_values[mid_sub-edge_from_center] = to_south_latlon[0]
                
                #for x and y
                self.x_values[mid_sub+edge_from_center] = x_dist*1000.0 #meters
                self.x_values[mid_sub-edge_from_center] = -1.0*x_dist*1000.0 #meters
                self.y_values[mid_sub+edge_from_center] = y_dist*1000.0 #meters
                self.y_values[mid_sub-edge_from_center] = -1.0*y_dist*1000.0 #meters
            self.have_xy_vals = True

        if self.n_horiz_grid_boxes % 2 == 1:
            print('...this section needs to be tested...')
            pdb.set_trace() ##add the self.x_values, self.y_values below as appropriate
            #odd # of boxes, even # of defining points
            #need to arrange center point at center of middle box
            #mid_subs: will have 2: just before &  just after center point
            mid_subs = [ ((self.n_horiz_grid_boxes+1)/2)-1 , (self.n_horiz_grid_boxes+1)/2 ]
   
            #first grid box edges to east & north will be center point + 0.5*spacing
            #first grid box edges to west & south will be center point - 0.5*spacing
            first_x_dist = (self.box_spacing/1000.0)*0.5  #[km]
            first_y_dist = (self.box_spacing/1000.0)*0.5  #[km]
            first_east_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, first_x_dist, 0.0)
            first_west_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, -1.0*first_x_dist, 0.0)
            first_north_latlon= sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, 0.0, first_y_dist)
            first_south_latlon= sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, 0.0, -1.0*first_y_dist)
            self.lon_values[mid_subs[0]] = first_west_latlon[1]
            self.lon_values[mid_subs[1]] = first_east_latlon[1]
            self.lat_values[mid_subs[0]] = first_south_latlon[0]
            self.lat_values[mid_subs[1]] = first_north_latlon[0]
   
            #rest of points to west & south: center - 0.5*spacing - (# of grid boxes to current sub * spacing)   
            #For subscripts less than the mid_subs
            for edge_sub in range(0, mid_subs[0]-1+1):
                #get # of full grid boxes from current subscript to 1st mid_subs
                n_full_boxes = mid_subs[0] - edge_sub
      
                #distance will be the first half spacing, then full spacing for each grid box
                x_dist = first_x_dist + (n_full_boxes*self.box_spacing)/1000.0  #[km]
                y_dist = first_y_dist + (n_full_boxes*self.box_spacing)/1000.0  #[km]
                to_west_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, -1.0*x_dist, 0.0)
                to_south_latlon= sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, 0.0, -1.0*y_dist)
                self.lon_values[edge_sub] = to_west_latlon[1]
                self.lat_values[edge_sub] = to_south_latlon[0]
   
            #rest of points to east & north: center + 0.5*spacing + (# of grid boxes to current sub * spacing)
            #For subscripts beyond the mid_subs
            for edge_sub in range(mid_subs[1]+1, self.n_horiz_grid_boxes+1):
                #get # of full grid boxes from current subscript to 2nd mid_subs
                n_full_boxes = edge_sub - mid_subs[1]
      
                #distance will be the first half spacing, then full spacing for each grid box
                x_dist = first_x_dist + (n_full_boxes*self.box_spacing)/1000.0  #[km]
                y_dist = first_y_dist + (n_full_boxes*self.box_spacing)/1000.0  #[km]
                to_east_latlon = sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, x_dist, 0.0)
                to_north_latlon= sim.xy2ll(self.cntr_lat_deg, self.cntr_lon_deg, 0.0, y_dist)
                self.lon_values[edge_sub] = to_east_latlon[1]
                self.lat_values[edge_sub] = to_north_latlon[0]
        self.have_xy_vals = True
                
        #Now that have the lat/lons for each grid point in the box, set the vertical coords:
        self.z_values = np.arange(0,self.n_vert_boxes+1)*float(self.vert_spacing)

        #define the dimensions for satellite and ground radar
        self.plat_dims = (self.z_values.shape[0], self.lon_values.shape[0], self.lat_values.shape[0])
        
        #define the bounds of the column box
        self.box_max_lat = self.lat_values.max()
        self.box_min_lat = self.lat_values.min()
        self.box_max_lon = self.lon_values.max()
        self.box_min_lon = self.lon_values.min()

        #Print out to terminal the grid coords for the full box:
        #all these either send in or returned from this function
        print(f'   Column Box Centered On: {self.center_on}')
        print(f'      lat: {self.cntr_lat_deg} lon: {self.cntr_lon_deg}')
        print('   Column Box Grid:')
        print(f'       horiz spacing [m]: {self.box_spacing}')
        print(f'       horiz extent [m]: {self.box_limit}')
        print(f'       horiz grid boxes: {self.n_horiz_grid_boxes}')
        print(f'       vert spacing [m]: {self.vert_spacing}')
        print(f'       vert extent [m]: {self.vert_limit}')
        print(f'       vert grid boxes: {self.n_vert_boxes}')
        print('   Column Box Grid Latitudes:')
        print(f'        {self.lat_values.tolist()}')
        print('   Column Box Grid Longitudes: ')
        print(f'        {self.lon_values.tolist()}')
        print('   Column Box Grid Vertical Levels [m]:')
        print(f'        {self.z_values.tolist()}')

    # ***************************************************************************************
    def get_grid_from_radar(self, file, grid_center, maxh=15000, maxxy=100000):
     
        radar = pyart.io.read(file, file_field_names=True)
        #print(radar.fields.keys())
        #os.system(f'rm {file}')
        
        #computer grid points based on column box resolution
        nz = int(maxh / self.vert_spacing + 1)
        ny = nx = int(maxxy * 2 / self.box_spacing + 1)
        
        # perform Cartesian mapping
        grid = pyart.map.grid_from_radars(radar,
                weighting_function='Barnes2', #function to weight for interpolation (check pyart docs)
                grid_origin=grid_center,      #lat and lon of grid origin
                grid_shape=(nz, ny, nx),      #number of points in grid (z, y, x)
                grid_limits=((0., maxh),      #min/max grid location in meters for ((z), (y), (x))
                             (-maxxy, maxxy),
                             (-maxxy, maxxy)))
        return grid

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
    plat_grp.setncattr('latitude_degrees',  plat_info['lat_d'])
    plat_grp.setncattr('latitude_minutes',  plat_info['lat_m'])
    plat_grp.setncattr('latitude_seconds',  plat_info['lat_s'])
    plat_grp.setncattr('longitude_degrees',  plat_info['lon_d'])
    plat_grp.setncattr('longitude_minutes',  plat_info['lon_m'])
    plat_grp.setncattr('longitude_seconds',  plat_info['lon_s'])
    plat_grp.setncattr('elevation_MSL', plat_info['elev'])
    plat_grp.setncattr('operation_mode',  plat_info['operation_mode'])
    plat_grp.setncattr('wavelength_m', plat_info['wavelength'])
    plat_grp.setncattr('frequency_GHz', plat_info['frequency'])
    plat_grp.setncattr('beam_width_deg', plat_info['beam_width'])
    plat_grp.setncattr('gate_size_m', plat_info['gate_size'])
    plat_grp.setncattr('timestamp',plat_info['timestamp'])
    plat_grp.setncattr('offset_vs_main',   plat_info['offset_vs_main'])
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

# *******************************************  M  A  I  N  ******************************

def main(inargs):
    """Run the program."""
    
    #extract input arguments
    year = int(inargs.year)
    month = int(inargs.month)
    day = int(inargs.day)

    # Import parameters to run pysimba
    if args.params_dict == None:
        print('***No parameter dictionary set, applying default parameters.***')
        kwargs = sim.get_default_params_dict()
    else:
        current_dict = args.params_dict
        kwargs = ast.literal_eval(open(current_dict).read())

    #run BuildColumn
    c = BuildColumn(year, month, day, **kwargs)
    column = c.process_columns()

    print('Program Done.', '', sep='\n')

if(__name__ == '__main__'):
            
    runargs = argparse.ArgumentParser(description='Produce precipitating columns over a user defined geographic location.')
    runargs.add_argument('year',  type=str, help='Provide the year in YYYY format')
    runargs.add_argument('month', type=str, help='Provide the month in MM format')
    runargs.add_argument('day',   type=str, help='Provide the day in DD format')
    runargs.add_argument('--params_dict', dest='params_dict', type=str, help='Parameter dictionary')

    args = runargs.parse_args()
    main(args)