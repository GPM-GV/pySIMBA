from datetime import datetime
import math
import h5py
import numpy as np
import pdb
import tools as sim

def get_dpr_for_column(dpr_file, main_plat_datetime, main_lat_deg, main_lon_deg, box_params, column):
    
    # Set initial dpr_plat_info values:
    dpr_info = sim.set_plat_values('DPR', dpr_file)

    #Set range to search w/in full DPR swath:
    limit_dpr_km = 10.0 * box_params.box_limit/1000.  #[km]
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
    dpr_data    = read_2adpr_hdf5(dpr_file)
    if dpr_data == -1:
        print('---Trouble reading 2ADPR file---')
        print('---Returning without GPM-DPR data in the column---')
        return column
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
    HS_hts_m = get_hgt_of_bins(HS_nPixels, HS_nScans, HS_nBins, HS_binSize, HS_elipsBO, HS_locZenAng)
    MS_hts_m = get_hgt_of_bins(MS_nPixels, MS_nScans, MS_nBins, MS_binSize, MS_elipsBO, MS_locZenAng)
    NS_hts_m = get_hgt_of_bins(NS_nPixels, NS_nScans, NS_nBins, NS_binSize, NS_elipsBO, NS_locZenAng)
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
    max_deg_lon = limit_dpr_km / (np.cos(np.deg2rad(box_params.cntr_lat_deg))*111.1)
    HS_start_scan=0  ;  HS_end_scan=0  ;  HS_nscans2do=0  ; HS_start_found=0
    MS_start_scan=0  ;  MS_end_scan=0  ;  MS_nscans2do=0  ; MS_start_found=0
    NS_start_scan=0  ;  NS_end_scan=0  ;  NS_nscans2do=0  ; NS_start_found=0
    
    # for each HS/MS/NS swath, get the [scan, pixel] subset in the search area:
    # created function to do this in pySIMBA 
    HS_start_scan, HS_end_scan, HS_nscans2do = get_nscan_to_process(HS_nScans, 
    HS_nPixels, HS_longitude, HS_latitude, limit_dpr_km, box_params.cntr_lat_deg, box_params.cntr_lon_deg)
    HS_nscans2do=0 #Set to 0 manually due to V05B version update as a temp fix - Updated 10/01/18
    if HS_nscans2do == 0:
        # this should NOT happen...
        print(' --- GPM 2ADPR MODULE ERROR!!! --- ')
        print(' --- CAN NOT LOCATE DPR HS SCAN WITH COLUMN GRID CENTER POINT! --- ')
        print(' --- SKIPPING HS DPR DATA! --- NO HS DPR VALUES INTO COLUMN! --- ')
        #pdb.set_trace()
    
    MS_start_scan, MS_end_scan, MS_nscans2do = get_nscan_to_process(MS_nScans, 
    MS_nPixels, MS_longitude, MS_latitude, limit_dpr_km, box_params.cntr_lat_deg, box_params.cntr_lon_deg)
    if MS_nscans2do == 0:
        #this should NOT happen...
        print(' --- GPM 2ADPR MODULE ERROR!!! --- ')
        print(' --- CAN NOT LOCATE DPR MS SCAN WITH COLUMN GRID CENTER POINT! --- ')
        print(' --- SKIPPING MS DPR DATA! --- NO MS DPR VALUES INTO COLUMN! --- ')
        #pdb.set_trace()

    NS_start_scan, NS_end_scan, NS_nscans2do = get_nscan_to_process(NS_nScans, 
    NS_nPixels, NS_longitude, NS_latitude, limit_dpr_km, box_params.cntr_lat_deg, box_params.cntr_lon_deg)
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
        print('---Returning without GPM-DPR data in the column---')
        return column

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
    if NS_longitude_subset.min() > min(box_params.lon_values): NS_lon_OK = 'N'
    if NS_latitude_subset.min()  > min(box_params.lat_values): NS_lat_OK = 'N'
    
    if MS_nscans2do > 0:
        if MS_longitude_subset.min() > min(box_params.lon_values): MS_lon_OK = 'N'
        if MS_latitude_subset.min()  > min(box_params.lat_values): MS_lat_OK = 'N'
    if HS_nscans2do > 0:
        if HS_longitude_subset.min() > min(box_params.lon_values): HS_lon_OK = 'N'
        if HS_latitude_subset.min()  > min(box_params.lat_values): HS_lat_OK = 'N'
    if NS_lon_OK == 'N' or NS_lat_OK == 'N':
        print(' --- APPEARS DPR for 2ADPR COVERAGE DOES NOT INCLUDE COLUMN GRID ! --- ')
        print(' --- exiting 2ADPR module, 2ADPR data NOT set in column grid.')
        print('---Returning without GPM-DPR data in the column---')
        return column
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
            temp_result = sim.get_posn_to_azm_range(scan_lats[s], scan_lons[s], main_lat_deg, main_lon_deg)
            scan_dists[p,s] = temp_result[1]
    # Use the subscript for minimum dist to main_plat to assign timestamp:
    min_subs_ar = np.argwhere(scan_dists == np.min(scan_dists))
    min_dist_p_sub = min_subs_ar[0][0]   #for nPixels dim
    min_dist_s_sub = min_subs_ar[0][1]   #for  nScans dim
    print(f'main plat [p,s]: {min_dist_p_sub}, {NS_start_scan+min_dist_s_sub}')
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

    dpr_info['timestamp'] = timestamp_year+timestamp_month+timestamp_day+'_'+ \
                                timestamp_hour+timestamp_min+timestamp_sec
    print(f'...DPR timestamp - main platform: {dpr_info["timestamp"]}')
    offset = (dpr_plat_datetime - main_plat_datetime).total_seconds()
    dpr_info['offset_vs_main'] = offset

    
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
            temp_result = sim.get_posn_to_azm_range(scan_lats[s], scan_lons[s], box_params.cntr_lat_deg, box_params.cntr_lon_deg)
            scan_dists[p,s]  = temp_result[1]
            
    # Use the subscript for minimum dist to grid center to assign timestamp:
    min_subs_ar = np.argwhere(scan_dists == np.min(scan_dists))
    min_dist_p_sub = min_subs_ar[0][0]   #for nPixels dim
    min_dist_s_sub = min_subs_ar[0][1]   #for  nScans dim
    print(f'col center [p,s]: {min_dist_p_sub}, {NS_start_scan+min_dist_s_sub}')
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

    dpr_info['timestamp_cntr'] = timestamp_year+timestamp_month+timestamp_day+'_'+ \
                                      timestamp_hour+timestamp_min+timestamp_sec
    print(f'...DPR timestamp - col grid cntr: {dpr_info["timestamp_cntr"]}')

    
    # Define column box grid arrays to populate w/ 2ADPR data, initialize w/ all NaNs
    col_HS_sfcType              = np.full(box_params.plat_dims, np.nan)
    col_HS_flagPrecip           = np.full(box_params.plat_dims, np.nan)
    col_HS_binRealSfc           = np.full(box_params.plat_dims, np.nan)
    col_HS_htStormTop           = np.full(box_params.plat_dims, np.nan)
    col_HS_htZeroDeg            = np.full(box_params.plat_dims, np.nan)
    col_HS_flagBB               = np.full(box_params.plat_dims, np.nan)
    col_HS_htBB                 = np.full(box_params.plat_dims, np.nan)
    col_HS_widthBB              = np.full(box_params.plat_dims, np.nan)
    col_HS_qualityBB            = np.full(box_params.plat_dims, np.nan)
    col_HS_typePrecip           = np.full(box_params.plat_dims, np.nan)
    col_HS_qualityTypePrecip    = np.full(box_params.plat_dims, np.nan)
    col_HS_PIAfinal             = np.full(box_params.plat_dims, np.nan)
    col_HS_corZFacNearSfc       = np.full(box_params.plat_dims, np.nan)
    col_HS_precipRateNearSfc    = np.full(box_params.plat_dims, np.nan)
    col_HS_precipRateAve24      = np.full(box_params.plat_dims, np.nan)
    col_HS_phaseNearSfc         = np.full(box_params.plat_dims, np.nan)
    col_HS_zFacMeas             = np.full(box_params.plat_dims, np.nan)
    col_HS_attenNoPrecip        = np.full(box_params.plat_dims, np.nan)
    col_HS_corZFac              = np.full(box_params.plat_dims, np.nan)
    col_HS_precipRate           = np.full(box_params.plat_dims, np.nan)
    col_HS_DSDphase             = np.full(box_params.plat_dims, np.nan)
    col_HS_PIA_cloudwater       = np.full(box_params.plat_dims, np.nan)
    col_HS_PIA_cloudice         = np.full(box_params.plat_dims, np.nan)
    col_HS_PIA_watervapor       = np.full(box_params.plat_dims, np.nan)
    col_HS_PIA_oxygen           = np.full(box_params.plat_dims, np.nan)
    col_HS_TPW_liquid           = np.full(box_params.plat_dims, np.nan)
    col_HS_TPW_ice              = np.full(box_params.plat_dims, np.nan)
    col_HS_DSD_dBNw             = np.full(box_params.plat_dims, np.nan)
    col_HS_DSD_Dm               = np.full(box_params.plat_dims, np.nan)
    col_HS_eff_PIA              = np.full(box_params.plat_dims, np.nan)
    
    col_NS_sfcType              = np.full(box_params.plat_dims, np.nan)
    col_NS_flagPrecip           = np.full(box_params.plat_dims, np.nan)
    col_NS_binRealSfc           = np.full(box_params.plat_dims, np.nan)
    col_NS_htStormTop           = np.full(box_params.plat_dims, np.nan)
    col_NS_htZeroDeg            = np.full(box_params.plat_dims, np.nan)
    col_NS_flagBB               = np.full(box_params.plat_dims, np.nan)
    col_NS_htBB                 = np.full(box_params.plat_dims, np.nan)
    col_NS_widthBB              = np.full(box_params.plat_dims, np.nan)
    col_NS_qualityBB            = np.full(box_params.plat_dims, np.nan)
    col_NS_typePrecip           = np.full(box_params.plat_dims, np.nan)
    col_NS_qualityTypePrecip    = np.full(box_params.plat_dims, np.nan)
    col_NS_PIAfinal             = np.full(box_params.plat_dims, np.nan)
    col_NS_corZFacNearSfc       = np.full(box_params.plat_dims, np.nan)
    col_NS_precipRateNearSfc    = np.full(box_params.plat_dims, np.nan)
    col_NS_precipRateAve24      = np.full(box_params.plat_dims, np.nan)
    col_NS_phaseNearSfc         = np.full(box_params.plat_dims, np.nan)
    col_NS_zFacMeas             = np.full(box_params.plat_dims, np.nan)
    col_NS_attenNoPrecip        = np.full(box_params.plat_dims, np.nan)
    col_NS_corZFac              = np.full(box_params.plat_dims, np.nan)
    col_NS_precipRate           = np.full(box_params.plat_dims, np.nan)
    col_NS_DSDphase             = np.full(box_params.plat_dims, np.nan)
    col_NS_PIA_cloudwater       = np.full(box_params.plat_dims, np.nan)
    col_NS_PIA_cloudice         = np.full(box_params.plat_dims, np.nan)
    col_NS_PIA_watervapor       = np.full(box_params.plat_dims, np.nan)
    col_NS_PIA_oxygen           = np.full(box_params.plat_dims, np.nan)
    col_NS_TPW_liquid           = np.full(box_params.plat_dims, np.nan)
    col_NS_TPW_ice              = np.full(box_params.plat_dims, np.nan)
    col_NS_DSD_dBNw             = np.full(box_params.plat_dims, np.nan)
    col_NS_DSD_Dm               = np.full(box_params.plat_dims, np.nan)
    col_NS_eff_PIA              = np.full(box_params.plat_dims, np.nan)
    
    col_MS_sfcType              = np.full(box_params.plat_dims, np.nan)
    col_MS_flagPrecip           = np.full(box_params.plat_dims, np.nan)
    col_MS_binRealSfc           = np.full(box_params.plat_dims, np.nan)
    col_MS_htStormTop           = np.full(box_params.plat_dims, np.nan)
    col_MS_htZeroDeg            = np.full(box_params.plat_dims, np.nan)
    col_MS_flagBB               = np.full(box_params.plat_dims, np.nan)
    col_MS_htBB                 = np.full(box_params.plat_dims, np.nan)
    col_MS_widthBB              = np.full(box_params.plat_dims, np.nan)
    col_MS_qualityBB            = np.full(box_params.plat_dims, np.nan)
    col_MS_typePrecip           = np.full(box_params.plat_dims, np.nan)
    col_MS_qualityTypePrecip    = np.full(box_params.plat_dims, np.nan)
    col_MS_PIAfinal             = np.full(box_params.plat_dims, np.nan)
    col_MS_corZFacNearSfc       = np.full(box_params.plat_dims, np.nan)
    col_MS_precipRateNearSfc    = np.full(box_params.plat_dims, np.nan)
    col_MS_precipRateAve24      = np.full(box_params.plat_dims, np.nan)
    col_MS_phaseNearSfc         = np.full(box_params.plat_dims, np.nan)
    col_MS_zFacMeas             = np.full(box_params.plat_dims, np.nan)
    col_MS_attenNoPrecip        = np.full(box_params.plat_dims, np.nan)
    col_MS_corZFac              = np.full(box_params.plat_dims, np.nan)
    col_MS_PIA_cloudwater       = np.full(box_params.plat_dims, np.nan)
    col_MS_PIA_cloudice         = np.full(box_params.plat_dims, np.nan)
    col_MS_PIA_watervapor       = np.full(box_params.plat_dims, np.nan)
    col_MS_PIA_oxygen           = np.full(box_params.plat_dims, np.nan)
    col_MS_TPW_liquid           = np.full(box_params.plat_dims, np.nan)
    col_MS_TPW_ice              = np.full(box_params.plat_dims, np.nan)
    col_MS_eff_PIA              = np.full(box_params.plat_dims, np.nan)

    # --- LOCATING DPR DATA FOR THE COLUMN GRID: ---
    # ----------------------------------------------
    # Get x & y coords for ea horiz column grid point
    n_col_pts = box_params.n_horiz_grid_boxes +1
    begin_at = box_params.box_limit/2.0
    col_x_km = (np.arange(n_col_pts)*box_params.box_spacing - begin_at)/1000. #[km]
    col_y_km = (np.arange(n_col_pts)*box_params.box_spacing - begin_at)/1000. #[km]
    
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
    dims = (box_params.lon_values.shape[0], box_params.lat_values.shape[0], HS_nBins)
    temp_HS_zFacMeas        = np.full(dims, np.nan)
    temp_HS_attenNoPrecip   = np.full(dims, np.nan)
    temp_HS_DSDphase        = np.full(dims, np.nan)
    temp_HS_corZFac         = np.full(dims, np.nan)
    temp_HS_precipRate      = np.full(dims, np.nan)
    temp_HS_DSD_dBNw        = np.full(dims, np.nan)
    temp_HS_DSD_Dm          = np.full(dims, np.nan)
    
    dims = (box_params.lon_values.shape[0], box_params.lat_values.shape[0], MS_nBins)
    temp_MS_zFacMeas        = np.full(dims, np.nan)
    temp_MS_attenNoPrecip   = np.full(dims, np.nan)
    temp_MS_corZFac         = np.full(dims, np.nan)
    
    dims = (box_params.lon_values.shape[0], box_params.lat_values.shape[0], NS_nBins)
    temp_NS_zFacMeas        = np.full(dims, np.nan)
    temp_NS_attenNoPrecip   = np.full(dims, np.nan)
    temp_NS_DSDphase        = np.full(dims, np.nan)
    temp_NS_corZFac         = np.full(dims, np.nan)
    temp_NS_precipRate      = np.full(dims, np.nan)
    temp_NS_DSD_dBNw        = np.full(dims, np.nan)
    temp_NS_DSD_Dm          = np.full(dims, np.nan)
    HS_hts_orig_m           = np.full((box_params.lon_values.shape[0], box_params.lat_values.shape[0], HS_nBins), np.nan)
    MS_hts_orig_m           = np.full((box_params.lon_values.shape[0], box_params.lat_values.shape[0], MS_nBins), np.nan)
    NS_hts_orig_m           = np.full((box_params.lon_values.shape[0], box_params.lat_values.shape[0], NS_nBins), np.nan)
    
    dims = (box_params.lon_values.shape[0], box_params.lat_values.shape[0])
    HS_binNum_CFB           = np.full(dims, np.nan)
    MS_binNum_CFB           = np.full(dims, np.nan)
    NS_binNum_CFB           = np.full(dims, np.nan)

    for i, val in enumerate(box_params.lon_values):
        for j, val2 in enumerate(box_params.lat_values):
    #for i, val in enumerate(self.column_box_params['column_grid_lats']):
    #    for j, val2 in enumerate(self.column_box_params['column_grid_lons']):
            NS_dpr2cpt_dist = np.sqrt( (NS_longitude_subset - box_params.lon_values[i])**2 + (NS_latitude_subset - box_params.lat_values[j])**2)
            
            if MS_nscans2do > 0 and col_in_MSswath > 0:
                MS_dpr2cpt_dist = np.sqrt( (MS_longitude_subset - box_params.lon_values[i])**2 + (MS_latitude_subset - box_params.lat_values[j])**2)
                
            if HS_nscans2do > 0 and col_in_HSswath > 0:
                HS_dpr2cpt_dist = np.sqrt( (HS_longitude_subset - box_params.lon_values[i])**2 + (HS_latitude_subset - box_params.lat_values[j])**2)

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
                col_HS_sfcType[0,j,i]           = HS_sfc_type[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_flagPrecip[0,j,i]        = HS_flag_precip[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_binRealSfc[0,j,i]        = HS_binRealSfc[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_htStormTop[0,j,i]        = HS_ht_stormTop[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_htZeroDeg[0,j,i]         = HS_ht_ZeroDeg[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_flagBB[0,j,i]            = HS_flag_BB[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_htBB[0,j,i]              = HS_ht_BB[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_widthBB[0,j,i]           = HS_width_BB[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_qualityBB[0,j,i]         = HS_quality_BB[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_typePrecip[0,j,i]        = HS_typePrecip[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_qualityTypePrecip[0,j,i] = HS_qualityTypePrecip[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_PIAfinal[0,j,i]          = HS_PIA_final[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_corZFacNearSfc[0,j,i]    = HS_corZFac_nearSfc[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_precipRateNearSfc[0,j,i] = HS_precipRate_nearSfc[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_precipRateAve24[0,j,i]   = HS_precipRateAve24[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_phaseNearSfc[0,j,i]      = HS_phase_nearSfc[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_PIA_cloudwater[0,j,i]    = HS_PIA_cloudwater[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_PIA_cloudice[0,j,i]      = HS_PIA_cloudice[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_PIA_watervapor[0,j,i]    = HS_PIA_watervapor[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_PIA_oxygen[0,j,i]        = HS_PIA_oxygen[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_TPW_liquid[0,j,i]        = HS_TPW_liquid[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_TPW_ice[0,j,i]           = HS_TPW_ice[HS_pt_s_sub, HS_pt_p_sub]
                col_HS_eff_PIA[0,j,i]           = HS_eff_PIA[HS_pt_s_sub, HS_pt_p_sub]
    
            if MS_nscans2do > 0 and col_in_MSswath > 0:
                col_MS_sfcType[0,j,i]           = MS_sfc_type[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_flagPrecip[0,j,i]        = MS_flag_precip[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_binRealSfc[0,j,i]        = MS_binRealSfc[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_htStormTop[0,j,i]        = MS_ht_stormTop[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_htZeroDeg[0,j,i]         = MS_ht_ZeroDeg[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_flagBB[0,j,i]            = MS_flag_BB[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_htBB[0,j,i]              = MS_ht_BB[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_widthBB[0,j,i]           = MS_width_BB[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_qualityBB[0,j,i]         = MS_quality_BB[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_typePrecip[0,j,i]        = MS_typePrecip[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_qualityTypePrecip[0,j,i] = MS_qualityTypePrecip[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_PIAfinal[0,j,i]          = MS_PIA_final[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_corZFacNearSfc[0,j,i]    = MS_corZFac_nearSfc[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_precipRateNearSfc[0,j,i] = MS_precipRate_nearSfc[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_precipRateAve24[0,j,i]   = MS_precipRateAve24[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_phaseNearSfc[0,j,i]      = MS_phase_nearSfc[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_PIA_cloudwater[0,j,i]    = MS_PIA_cloudwater[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_PIA_cloudice[0,j,i]      = MS_PIA_cloudice[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_PIA_watervapor[0,j,i]    = MS_PIA_watervapor[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_PIA_oxygen[0,j,i]        = MS_PIA_oxygen[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_TPW_liquid[0,j,i]        = MS_TPW_liquid[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_TPW_ice[0,j,i]           = MS_TPW_ice[MS_pt_s_sub, MS_pt_p_sub]
                col_MS_eff_PIA[0,j,i]           = MS_eff_PIA[MS_pt_s_sub, MS_pt_p_sub]
      
            col_NS_sfcType[0,j,i]             = NS_sfc_type[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_flagPrecip[0,j,i]          = NS_flag_precip[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_binRealSfc[0,j,i]          = NS_binRealSfc[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_htStormTop[0,j,i]          = NS_ht_stormTop[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_htZeroDeg[0,j,i]           = NS_ht_ZeroDeg[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_flagBB[0,j,i]              = NS_flag_BB[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_htBB[0,j,i]                = NS_ht_BB[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_widthBB[0,j,i]             = NS_width_BB[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_qualityBB[0,j,i]           = NS_quality_BB[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_typePrecip[0,j,i]          = NS_typePrecip[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_qualityTypePrecip[0,j,i]   = NS_qualityTypePrecip[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_PIAfinal[0,j,i]            = NS_PIA_final[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_corZFacNearSfc[0,j,i]      = NS_corZFac_nearSfc[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_precipRateNearSfc[0,j,i]   = NS_precipRate_nearSfc[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_precipRateAve24[0,j,i]     = NS_precipRateAve24[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_phaseNearSfc[0,j,i]        = NS_phase_nearSfc[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_PIA_cloudwater[0,j,i]      = NS_PIA_cloudwater[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_PIA_cloudice[0,j,i]        = NS_PIA_cloudice[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_PIA_watervapor[0,j,i]      = NS_PIA_watervapor[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_PIA_oxygen[0,j,i]          = NS_PIA_oxygen[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_TPW_liquid[0,j,i]          = NS_TPW_liquid[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_TPW_ice[0,j,i]             = NS_TPW_ice[NS_pt_s_sub, NS_pt_p_sub]
            col_NS_eff_PIA[0,j,i]             = NS_eff_PIA[NS_pt_s_sub, NS_pt_p_sub]

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
    for i, val in enumerate(box_params.lon_values):
        for j, val2 in enumerate(box_params.lat_values):
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
        col_HS_sfcType[col_HS_sfcType <= badVal]                    = np.nan
        col_HS_flagPrecip[col_HS_flagPrecip <= badVal]              = np.nan
        col_HS_binRealSfc[col_HS_binRealSfc <= badVal]              = np.nan
        col_HS_htStormTop[col_HS_htStormTop <= badVal]              = np.nan
        col_HS_htZeroDeg[col_HS_htZeroDeg <= badVal]                = np.nan
        col_HS_flagBB[col_HS_flagBB <= badVal]                      = np.nan
        col_HS_htBB[col_HS_htBB <= badVal]                          = np.nan
        col_HS_widthBB[col_HS_widthBB <= badVal]                    = np.nan
        col_HS_qualityBB[col_HS_qualityBB <= badVal]                = np.nan
        col_HS_typePrecip[col_HS_typePrecip <= badVal]              = np.nan
        col_HS_qualityTypePrecip[col_HS_qualityTypePrecip <= badVal]= np.nan
        col_HS_PIAfinal[col_HS_PIAfinal <= badVal]                  = np.nan
        col_HS_corZFacNearSfc[col_HS_corZFacNearSfc <= badVal]      = np.nan
        col_HS_precipRateNearSfc[col_HS_precipRateNearSfc <= badVal]= np.nan
        col_HS_precipRateAve24[col_HS_precipRateAve24 <= badVal]    = np.nan
        col_HS_phaseNearSfc[col_HS_phaseNearSfc == 255]             = np.nan
        col_HS_PIA_cloudwater[col_HS_PIA_cloudwater <= badVal]      = np.nan
        col_HS_PIA_cloudice[col_HS_PIA_cloudice <= badVal]          = np.nan
        col_HS_PIA_watervapor[col_HS_PIA_watervapor <= badVal]      = np.nan
        col_HS_PIA_oxygen[col_HS_PIA_oxygen <= badVal]              = np.nan
        col_HS_TPW_liquid[col_HS_TPW_liquid <= badVal]              = np.nan
        col_HS_TPW_ice[col_HS_TPW_ice <= badVal]                    = np.nan
        col_HS_eff_PIA[col_HS_eff_PIA <= badVal]                    = np.nan
        temp_HS_zFacMeas[temp_HS_zFacMeas <= badVal]                = np.nan
        temp_HS_attenNoPrecip[temp_HS_attenNoPrecip <= badVal]      = np.nan
        temp_HS_corZFac[temp_HS_corZFac <= badVal]                  = np.nan
        temp_HS_precipRate[temp_HS_precipRate <= badVal]            = np.nan
        temp_HS_DSDphase[temp_HS_DSDphase == 255]                   = np.nan
        temp_HS_DSD_dBNw[temp_HS_DSD_dBNw <= badVal]                = np.nan
        temp_HS_DSD_Dm[temp_HS_DSD_Dm <= badVal]                    = np.nan

    if MS_nscans2do > 0 and col_in_MSswath > 0:
        col_MS_sfcType[col_MS_sfcType <= badVal]                    = np.nan
        col_MS_flagPrecip[col_MS_flagPrecip <= badVal]              = np.nan
        col_MS_binRealSfc[col_MS_binRealSfc <= badVal]              = np.nan
        col_MS_htStormTop[col_MS_htStormTop <= badVal]              = np.nan
        col_MS_htZeroDeg[col_MS_htZeroDeg <= badVal]                = np.nan
        col_MS_flagBB[col_MS_flagBB <= badVal]                      = np.nan
        col_MS_htBB[col_MS_htBB <= badVal]                          = np.nan
        col_MS_widthBB[col_MS_widthBB <= badVal]                    = np.nan
        col_MS_qualityBB[col_MS_qualityBB <= badVal]                = np.nan
        col_MS_typePrecip[col_MS_typePrecip <= badVal]              = np.nan
        col_MS_qualityTypePrecip[col_MS_qualityTypePrecip <= badVal]= np.nan
        col_MS_PIAfinal[col_MS_PIAfinal <= badVal]                  = np.nan
        col_MS_corZFacNearSfc[col_MS_corZFacNearSfc <= badVal]      = np.nan
        col_MS_precipRateNearSfc[col_MS_precipRateNearSfc <= badVal]= np.nan
        col_MS_precipRateAve24[col_MS_precipRateAve24 <= badVal]    = np.nan
        col_MS_phaseNearSfc[col_MS_phaseNearSfc == 255]             = np.nan
        col_MS_PIA_cloudwater[col_MS_PIA_cloudwater <= badVal]      = np.nan
        col_MS_PIA_cloudice[col_MS_PIA_cloudice <= badVal]          = np.nan
        col_MS_PIA_watervapor[col_MS_PIA_watervapor <= badVal]      = np.nan
        col_MS_PIA_oxygen[col_MS_PIA_oxygen <= badVal]              = np.nan
        col_MS_TPW_liquid[col_MS_TPW_liquid <= badVal]              = np.nan
        col_MS_TPW_ice[col_MS_TPW_ice <= badVal]                    = np.nan
        col_MS_eff_PIA[col_MS_eff_PIA <= badVal]                    = np.nan
        temp_MS_zFacMeas[temp_MS_zFacMeas <= badVal]                = np.nan
        temp_MS_attenNoPrecip[temp_MS_attenNoPrecip <= badVal]      = np.nan
        temp_MS_corZFac[temp_MS_corZFac <= badVal]                  = np.nan

    col_NS_sfcType[col_NS_sfcType <= badVal]                        = np.nan
    col_NS_flagPrecip[col_NS_flagPrecip <= badVal]                  = np.nan
    col_NS_binRealSfc[col_NS_binRealSfc <= badVal]                  = np.nan
    col_NS_htStormTop[col_NS_htStormTop <= badVal]                  = np.nan
    col_NS_htZeroDeg[col_NS_htZeroDeg <= badVal]                    = np.nan
    col_NS_flagBB[col_NS_flagBB <= badVal]                          = np.nan
    col_NS_htBB[col_NS_htBB <= badVal]                              = np.nan
    col_NS_widthBB[col_NS_widthBB <= badVal]                        = np.nan
    col_NS_qualityBB[col_NS_qualityBB <= badVal]                    = np.nan
    col_NS_typePrecip[col_NS_typePrecip <= badVal]                  = np.nan
    col_NS_qualityTypePrecip[col_NS_qualityTypePrecip <= badVal]    = np.nan
    col_NS_PIAfinal[col_NS_PIAfinal <= badVal]                      = np.nan
    col_NS_corZFacNearSfc[col_NS_corZFacNearSfc <= badVal]          = np.nan
    col_NS_precipRateNearSfc[col_NS_precipRateNearSfc <= badVal]    = np.nan
    col_NS_precipRateAve24[col_NS_precipRateAve24 <= badVal]        = np.nan
    col_NS_phaseNearSfc[col_NS_phaseNearSfc == 255]                 = np.nan
    col_NS_PIA_cloudwater[col_NS_PIA_cloudwater <= badVal]          = np.nan
    col_NS_PIA_cloudice[col_NS_PIA_cloudice <= badVal]              = np.nan
    col_NS_PIA_watervapor[col_NS_PIA_watervapor <= badVal]           = np.nan
    col_NS_PIA_oxygen[col_NS_PIA_oxygen <= badVal]                  = np.nan
    col_NS_TPW_liquid[col_NS_TPW_liquid <= badVal]                  = np.nan
    col_NS_TPW_ice[col_NS_TPW_ice <= badVal]                        = np.nan
    col_NS_eff_PIA[col_NS_eff_PIA <= badVal]                        = np.nan
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
    for i, val in enumerate(box_params.lon_values):
        for j, val2 in enumerate(box_params.lat_values):
            if HS_nscans2do > 0 and col_in_HSswath > 0:
                col_HS_zFacMeas[:,j,i]      = sim.interp(temp_HS_zFacMeas[j,i,:],     HS_hts_orig_m[j,i,:], box_params.z_values)
                col_HS_attenNoPrecip[:,j,i] = sim.interp(temp_HS_attenNoPrecip[j,i,:],HS_hts_orig_m[j,i,:], box_params.z_values)
                col_HS_corZFac[:,j,i]       = sim.interp(temp_HS_corZFac[j,i,:],      HS_hts_orig_m[j,i,:], box_params.z_values)
                col_HS_precipRate[:,j,i]    = sim.interp(temp_HS_precipRate[j,i,:],   HS_hts_orig_m[j,i,:], box_params.z_values)
                col_HS_DSDphase[:,j,i]      = sim.interp(temp_HS_DSDphase[j,i,:],     HS_hts_orig_m[j,i,:], box_params.z_values)
                col_HS_DSD_dBNw[:,j,i]      = sim.interp(temp_HS_DSD_dBNw[j,i,:],     HS_hts_orig_m[j,i,:], box_params.z_values)
                col_HS_DSD_Dm[:,j,i]        = sim.interp(temp_HS_DSD_Dm[j,i,:],       HS_hts_orig_m[j,i,:], box_params.z_values)
            if MS_nscans2do > 0 and col_in_MSswath > 0:
                col_MS_zFacMeas[:,j,i]      = sim.interp(temp_MS_zFacMeas[j,i,:],     MS_hts_orig_m[j,i,:], box_params.z_values)
                col_MS_attenNoPrecip[:,j,i] = sim.interp(temp_MS_attenNoPrecip[j,i,:],MS_hts_orig_m[j,i,:], box_params.z_values)
                col_MS_corZFac[:,j,i]       = sim.interp(temp_MS_corZFac[j,i,:],      MS_hts_orig_m[j,i,:], box_params.z_values)
                
            col_NS_zFacMeas[:,j,i]          = sim.interp(temp_NS_zFacMeas[j,i,:],     NS_hts_orig_m[j,i,:], box_params.z_values)
            col_NS_attenNoPrecip[:,j,i]     = sim.interp(temp_NS_attenNoPrecip[j,i,:],NS_hts_orig_m[j,i,:], box_params.z_values)
            col_NS_corZFac[:,j,i]           = sim.interp(temp_NS_corZFac[j,i,:],      NS_hts_orig_m[j,i,:], box_params.z_values)
            col_NS_precipRate[:,j,i]        = sim.interp(temp_NS_precipRate[j,i,:],   NS_hts_orig_m[j,i,:], box_params.z_values)
            col_NS_DSDphase[:,j,i]          = sim.interp(temp_NS_DSDphase[j,i,:],     NS_hts_orig_m[j,i,:], box_params.z_values)
            col_NS_DSD_dBNw[:,j,i]          = sim.interp(temp_NS_DSD_dBNw[j,i,:],     NS_hts_orig_m[j,i,:], box_params.z_values)
            col_NS_DSD_Dm[:,j,i]            = sim.interp(temp_NS_DSD_Dm[j,i,:],       NS_hts_orig_m[j,i,:], box_params.z_values)
    
    #remove 10*log(Nw) multiplier from Nw fields:
    if HS_nscans2do > 0 and col_in_HSswath > 0: col_HS_DSD_dBNw = col_HS_DSD_dBNw/10.0
    col_NS_DSD_dBNw = col_NS_DSD_dBNw/10.0
    
    # Define units, more descriptive names for each field:
    dpr_name_sfcType    = 'land sruface type category'
    dpr_units_sfcType   = 'int value for category'
    
    dpr_name_flagPrecip     = 'flag if precip or no precip'
    dpr_units_flagPrecip    = '0: no precip  11: precip'
    
    dpr_name_htStormTop     = 'ht of storm top'
    dpr_units_htStormTop    = 'm'
    
    dpr_name_htZeroDeg      = 'ht of freezing level'
    dpr_units_htZeroDeg     = 'm'
    
    dpr_name_flagBB     = 'tells bright band existance'
    dpr_units_flagBB    = '0: no BB  1: yes BB'
    
    dpr_name_htBB   = 'ht of bright band'
    dpr_units_htBB  = 'm'
    
    dpr_name_widthBB    = 'width of bright band'
    dpr_units_widthBB   = 'm'
    
    dpr_name_qualityBB  = 'how well defined bright band is'
    dpr_units_qualityBB = '1: best'
    
    dpr_name_typePrecip  = '8 char precip description - 1st char tells major type'
    dpr_units_typePrecip = '1:stratiform  2:convective 3:other'
    
    dpr_name_qualityTypePrecip  = 'quality of precip type'
    dpr_units_qualityTypePrecip = '1: best'
    
    dpr_name_PIAfinal  = 'final est of path integrated attenuation due to precip particles'
    dpr_units_PIAfinal = 'dB'
    
    dpr_name_corZFacNearSfc  = 'reflectivity factor corrected for attenuation near surface'
    dpr_units_corZFacNearSfc = 'dBZ'
    
    dpr_name_precipRateNearSfc  = 'precip rate near the surface'
    dpr_units_precipRateNearSfc = 'mm h^-1'
    
    dpr_name_precipRateAve24  = 'average precip rate for 2-4 km ht'
    dpr_units_precipRateAve24 = 'mm h^-1'
    
    dpr_name_phaseNearSfc  = 'phase of precip near surface'
    dpr_units_phaseNearSfc = 'for int value/100:  0:solid  1:mixed  2:liquid'
    
    dpr_name_zFacMeas  = 'reflectivity factor with no attenuation correction'
    dpr_units_zFacMeas = 'dBZ'
    
    dpr_name_attenNoPrecip  = 'attenuation by non-precip particles'
    dpr_units_attenNoPrecip = 'dB km^-1'
    
    dpr_name_corZFac  = 'reflectivity factor corrected for attenuation'
    dpr_units_corZFac = 'dBZ'
    
    dpr_name_precipRate  = 'precipitation rate'
    dpr_units_precipRate = 'mm h^-1'
    
    dpr_name_DSDphase  = 'phase of precip'
    dpr_units_DSDphase = 'for int value/100:  0:solid  1:mixed  2:liquid'
    
    dpr_name_PIA_cloudwater  = 'PIA due to non precip: cloud water'
    dpr_units_PIA_cloudwater = 'dB'
    
    dpr_name_PIA_cloudice  = 'PIA due to non precip: cloud ice'
    dpr_units_PIA_cloudice = 'dB'
    
    dpr_name_PIA_watervapor  = 'PIA due to non precip: water vapor'
    dpr_units_PIA_watervapor = 'dB'
    
    dpr_name_PIA_oxygen  = 'PIA due to non precip: oxygen molecules'
    dpr_units_PIA_oxygen = 'dB'
    
    dpr_name_TPW_liquid = 'vertically integrated liquid precip water'
    dpr_units_TPW_liquid= 'g m^-2'
    
    dpr_name_TPW_ice    = 'vertically integratied solid precip water'
    dpr_units_TPW_ice   = 'g m^-2'
    
    dpr_name_DSD_dBNw   = 'DSD normalized intercept parameter'
    dpr_units_DSD_dBNw  = 'log(Nw)'
    
    dpr_name_DSD_Dm = 'DSD mass-weighted mean diameter'
    dpr_units_DSD_Dm    = 'mm'
    
    dpr_name_effPIA	= 'effective 2-way PIA'
    dpr_units_effPIA	= 'dB'
    
    #assign each field to column object
    
    #####################################################################################
    ##############         HS          ##################
    #####################################################################################
    column.add_variable_to_object(col_HS_sfcType, 
                                  var_name='dpr_HS_sfcType', 
                                  units=dpr_units_sfcType,
                                  long_name=dpr_name_sfcType)

    column.add_variable_to_object(col_HS_flagPrecip, 
                                  var_name='dpr_HS_flagPrecip', 
                                  units=dpr_units_flagPrecip,
                                  long_name=dpr_name_flagPrecip)

    column.add_variable_to_object(col_HS_htStormTop, 
                                  var_name='dpr_HS_htStormTop', 
                                  units=dpr_units_htStormTop,
                                  long_name=dpr_name_htStormTop)

    column.add_variable_to_object(col_HS_htZeroDeg, 
                                  var_name='dpr_HS_htZeroDeg', 
                                  units=dpr_units_htZeroDeg,
                                  long_name=dpr_name_htZeroDeg)

    column.add_variable_to_object(col_HS_flagBB, 
                                  var_name='dpr_HS_flagBB', 
                                  units=dpr_units_flagBB,
                                  long_name=dpr_name_flagBB)

    column.add_variable_to_object(col_HS_htBB, 
                                  var_name='dpr_HS_htBB', 
                                  units=dpr_units_htBB,
                                  long_name=dpr_name_htBB)

    column.add_variable_to_object(col_HS_widthBB, 
                                  var_name='dpr_HS_widthBB', 
                                  units=dpr_units_widthBB,
                                  long_name=dpr_name_widthBB)

    column.add_variable_to_object(col_HS_qualityBB, 
                                  var_name='dpr_HS_qualityBB', 
                                  units=dpr_units_qualityBB,
                                  long_name=dpr_name_qualityBB)

    column.add_variable_to_object(col_HS_typePrecip, 
                                  var_name='dpr_HS_typePrecip', 
                                  units=dpr_units_typePrecip,
                                  long_name=dpr_name_typePrecip)

    column.add_variable_to_object(col_HS_qualityTypePrecip, 
                                  var_name='dpr_HS_qualityTypePrecip', 
                                  units=dpr_units_qualityTypePrecip,
                                  long_name=dpr_name_qualityTypePrecip)

    column.add_variable_to_object(col_HS_PIAfinal, 
                                  var_name='dpr_HS_PIAfinal', 
                                  units=dpr_units_PIAfinal,
                                  long_name=dpr_name_PIAfinal)

    column.add_variable_to_object(col_HS_corZFacNearSfc, 
                                  var_name='dpr_HS_corZFacNearSfc', 
                                  units=dpr_units_corZFacNearSfc,
                                  long_name=dpr_name_corZFacNearSfc)

    column.add_variable_to_object(col_HS_precipRateNearSfc, 
                                  var_name='dpr_HS_precipRateNearSfc', 
                                  units=dpr_units_precipRateNearSfc,
                                  long_name=dpr_name_precipRateNearSfc)

    column.add_variable_to_object(col_HS_precipRateAve24, 
                                  var_name='dpr_HS_precipRateAve24', 
                                  units=dpr_units_precipRateAve24,
                                  long_name=dpr_name_precipRateAve24)

    column.add_variable_to_object(col_HS_phaseNearSfc, 
                                  var_name='dpr_HS_phaseNearSfc', 
                                  units=dpr_units_phaseNearSfc,
                                  long_name=dpr_name_phaseNearSfc)

    column.add_variable_to_object(col_HS_zFacMeas, 
                                  var_name='dpr_HS_zFacMeas', 
                                  units=dpr_units_zFacMeas,
                                  long_name=dpr_name_zFacMeas)

    column.add_variable_to_object(col_HS_attenNoPrecip, 
                                  var_name='dpr_HS_attenNoPrecip', 
                                  units=dpr_units_attenNoPrecip,
                                  long_name=dpr_name_attenNoPrecip)

    column.add_variable_to_object(col_HS_corZFac, 
                                  var_name='dpr_HS_corZFac', 
                                  units=dpr_units_corZFac,
                                  long_name=dpr_name_corZFac)

    column.add_variable_to_object(col_HS_precipRate, 
                                  var_name='dpr_HS_precipRate', 
                                  units=dpr_units_precipRate,
                                  long_name=dpr_name_precipRate)

    column.add_variable_to_object(col_HS_DSDphase, 
                                  var_name='dpr_HS_DSDphase', 
                                  units=dpr_units_DSDphase,
                                  long_name=dpr_name_DSDphase)

    column.add_variable_to_object(col_HS_PIA_cloudwater, 
                                  var_name='dpr_HS_PIA_cloudwater', 
                                  units=dpr_units_PIA_cloudwater,
                                  long_name=dpr_name_PIA_cloudwater)

    column.add_variable_to_object(col_HS_PIA_cloudice, 
                                  var_name='dpr_HS_PIA_cloudice', 
                                  units=dpr_units_PIA_cloudice,
                                  long_name=dpr_name_PIA_cloudice)

    column.add_variable_to_object(col_HS_PIA_watervapor, 
                                  var_name='dpr_HS_PIA_watervapor', 
                                  units=dpr_units_PIA_watervapor,
                                  long_name=dpr_name_PIA_watervapor)

    column.add_variable_to_object(col_HS_PIA_oxygen, 
                                  var_name='dpr_HS_PIA_oxygen', 
                                  units=dpr_units_PIA_oxygen,
                                  long_name=dpr_name_PIA_oxygen)

    column.add_variable_to_object(col_HS_TPW_liquid, 
                                  var_name='dpr_HS_TPW_liquid', 
                                  units=dpr_units_TPW_liquid,
                                  long_name=dpr_name_TPW_liquid)

    column.add_variable_to_object(col_HS_TPW_ice, 
                                  var_name='dpr_HS_TPW_ice', 
                                  units=dpr_units_TPW_ice,
                                  long_name=dpr_name_TPW_ice)

    column.add_variable_to_object(col_HS_DSD_dBNw, 
                                  var_name='dpr_HS_dBNw', 
                                  units=dpr_units_DSD_dBNw,
                                  long_name=dpr_name_DSD_dBNw)

    column.add_variable_to_object(col_HS_DSD_Dm, 
                                  var_name='dpr_HS_Dm', 
                                  units=dpr_units_DSD_Dm,
                                  long_name=dpr_name_DSD_Dm)

    column.add_variable_to_object(col_HS_eff_PIA, 
                                  var_name='dpr_HS_eff_PIA', 
                                  units=dpr_units_effPIA,
                                  long_name=dpr_name_effPIA)

    #####################################################################################
    ##############         NS          ##################
    #####################################################################################

    column.add_variable_to_object(col_NS_sfcType, 
                                  var_name='dpr_NS_sfcType', 
                                  units=dpr_units_sfcType,
                                  long_name=dpr_name_sfcType)

    column.add_variable_to_object(col_NS_flagPrecip, 
                                  var_name='dpr_NS_flagPrecip', 
                                  units=dpr_units_flagPrecip,
                                  long_name=dpr_name_flagPrecip)

    column.add_variable_to_object(col_NS_htStormTop, 
                                  var_name='dpr_NS_htStormTop', 
                                  units=dpr_units_htStormTop,
                                  long_name=dpr_name_htStormTop)

    column.add_variable_to_object(col_NS_htZeroDeg, 
                                  var_name='dpr_NS_htZeroDeg', 
                                  units=dpr_units_htZeroDeg,
                                  long_name=dpr_name_htZeroDeg)

    column.add_variable_to_object(col_NS_flagBB, 
                                  var_name='dpr_NS_flagBB', 
                                  units=dpr_units_flagBB,
                                  long_name=dpr_name_flagBB)

    column.add_variable_to_object(col_NS_htBB, 
                                  var_name='dpr_NS_htBB', 
                                  units=dpr_units_htBB,
                                  long_name=dpr_name_htBB)

    column.add_variable_to_object(col_NS_widthBB, 
                                  var_name='dpr_NS_widthBB', 
                                  units=dpr_units_widthBB,
                                  long_name=dpr_name_widthBB)

    column.add_variable_to_object(col_NS_qualityBB, 
                                  var_name='dpr_NS_qualityBB', 
                                  units=dpr_units_qualityBB,
                                  long_name=dpr_name_qualityBB)

    column.add_variable_to_object(col_NS_typePrecip, 
                                  var_name='dpr_NS_typePrecip', 
                                  units=dpr_units_typePrecip,
                                  long_name=dpr_name_typePrecip)

    column.add_variable_to_object(col_NS_qualityTypePrecip, 
                                  var_name='dpr_NS_qualityTypePrecip', 
                                  units=dpr_units_qualityTypePrecip,
                                  long_name=dpr_name_qualityTypePrecip)

    column.add_variable_to_object(col_NS_PIAfinal, 
                                  var_name='dpr_NS_PIAfinal', 
                                  units=dpr_units_PIAfinal,
                                  long_name=dpr_name_PIAfinal)

    column.add_variable_to_object(col_NS_corZFacNearSfc, 
                                  var_name='dpr_NS_corZFacNearSfc', 
                                  units=dpr_units_corZFacNearSfc,
                                  long_name=dpr_name_corZFacNearSfc)

    column.add_variable_to_object(col_NS_precipRateNearSfc, 
                                  var_name='dpr_NS_precipRateNearSfc', 
                                  units=dpr_units_precipRateNearSfc,
                                  long_name=dpr_name_precipRateNearSfc)

    column.add_variable_to_object(col_NS_precipRateAve24, 
                                  var_name='dpr_NS_precipRateAve24', 
                                  units=dpr_units_precipRateAve24,
                                  long_name=dpr_name_precipRateAve24)

    column.add_variable_to_object(col_NS_phaseNearSfc, 
                                  var_name='dpr_NS_phaseNearSfc', 
                                  units=dpr_units_phaseNearSfc,
                                  long_name=dpr_name_phaseNearSfc)

    column.add_variable_to_object(col_NS_zFacMeas, 
                                  var_name='dpr_NS_zFacMeas', 
                                  units=dpr_units_zFacMeas,
                                  long_name=dpr_name_zFacMeas)

    column.add_variable_to_object(col_NS_attenNoPrecip, 
                                  var_name='dpr_NS_attenNoPrecip', 
                                  units=dpr_units_attenNoPrecip,
                                  long_name=dpr_name_attenNoPrecip)

    column.add_variable_to_object(col_NS_corZFac, 
                                  var_name='dpr_NS_corZFac', 
                                  units=dpr_units_corZFac,
                                  long_name=dpr_name_corZFac)

    column.add_variable_to_object(col_NS_precipRate, 
                                  var_name='dpr_NS_precipRate', 
                                  units=dpr_units_precipRate,
                                  long_name=dpr_name_precipRate)

    column.add_variable_to_object(col_NS_DSDphase, 
                                  var_name='dpr_NS_DSDphase', 
                                  units=dpr_units_DSDphase,
                                  long_name=dpr_name_DSDphase)

    column.add_variable_to_object(col_NS_PIA_cloudwater, 
                                  var_name='dpr_NS_PIA_cloudwater', 
                                  units=dpr_units_PIA_cloudwater,
                                  long_name=dpr_name_PIA_cloudwater)

    column.add_variable_to_object(col_NS_PIA_cloudice, 
                                  var_name='dpr_NS_PIA_cloudice', 
                                  units=dpr_units_PIA_cloudice,
                                  long_name=dpr_name_PIA_cloudice)

    column.add_variable_to_object(col_NS_PIA_watervapor, 
                                  var_name='dpr_NS_PIA_watervapor', 
                                  units=dpr_units_PIA_watervapor,
                                  long_name=dpr_name_PIA_watervapor)

    column.add_variable_to_object(col_NS_PIA_oxygen, 
                                  var_name='dpr_NS_PIA_oxygen', 
                                  units=dpr_units_PIA_oxygen,
                                  long_name=dpr_name_PIA_oxygen)

    column.add_variable_to_object(col_NS_TPW_liquid, 
                                  var_name='dpr_NS_TPW_liquid', 
                                  units=dpr_units_TPW_liquid,
                                  long_name=dpr_name_TPW_liquid)

    column.add_variable_to_object(col_NS_TPW_ice, 
                                  var_name='dpr_NS_TPW_ice', 
                                  units=dpr_units_TPW_ice,
                                  long_name=dpr_name_TPW_ice)

    column.add_variable_to_object(col_NS_DSD_dBNw, 
                                  var_name='dpr_NS_dBNw', 
                                  units=dpr_units_DSD_dBNw,
                                  long_name=dpr_name_DSD_dBNw)

    column.add_variable_to_object(col_NS_DSD_Dm, 
                                  var_name='dpr_NS_Dm', 
                                  units=dpr_units_DSD_Dm,
                                  long_name=dpr_name_DSD_Dm)

    column.add_variable_to_object(col_NS_eff_PIA, 
                                  var_name='dpr_NS_eff_PIA', 
                                  units=dpr_units_effPIA,
                                  long_name=dpr_name_effPIA)

    #####################################################################################
    ##############         MS          ##################
    #####################################################################################

    column.add_variable_to_object(col_MS_sfcType, 
                                  var_name='dpr_MS_sfcType', 
                                  units=dpr_units_sfcType,
                                  long_name=dpr_name_sfcType)

    column.add_variable_to_object(col_MS_flagPrecip, 
                                  var_name='dpr_MS_flagPrecip', 
                                  units=dpr_units_flagPrecip,
                                  long_name=dpr_name_flagPrecip)

    column.add_variable_to_object(col_MS_htStormTop, 
                                  var_name='dpr_MS_htStormTop', 
                                  units=dpr_units_htStormTop,
                                  long_name=dpr_name_htStormTop)

    column.add_variable_to_object(col_MS_htZeroDeg, 
                                  var_name='dpr_MS_htZeroDeg', 
                                  units=dpr_units_htZeroDeg,
                                  long_name=dpr_name_htZeroDeg)

    column.add_variable_to_object(col_MS_flagBB, 
                                  var_name='dpr_MS_flagBB', 
                                  units=dpr_units_flagBB,
                                  long_name=dpr_name_flagBB)

    column.add_variable_to_object(col_MS_htBB, 
                                  var_name='dpr_MS_htBB', 
                                  units=dpr_units_htBB,
                                  long_name=dpr_name_htBB)

    column.add_variable_to_object(col_MS_widthBB, 
                                  var_name='dpr_MS_widthBB', 
                                  units=dpr_units_widthBB,
                                  long_name=dpr_name_widthBB)

    column.add_variable_to_object(col_MS_qualityBB, 
                                  var_name='dpr_MS_qualityBB', 
                                  units=dpr_units_qualityBB,
                                  long_name=dpr_name_qualityBB)

    column.add_variable_to_object(col_MS_typePrecip, 
                                  var_name='dpr_MS_typePrecip', 
                                  units=dpr_units_typePrecip,
                                  long_name=dpr_name_typePrecip)

    column.add_variable_to_object(col_MS_qualityTypePrecip, 
                                  var_name='dpr_MS_qualityTypePrecip', 
                                  units=dpr_units_qualityTypePrecip,
                                  long_name=dpr_name_qualityTypePrecip)

    column.add_variable_to_object(col_MS_PIAfinal, 
                                  var_name='dpr_MS_PIAfinal', 
                                  units=dpr_units_PIAfinal,
                                  long_name=dpr_name_PIAfinal)

    column.add_variable_to_object(col_MS_corZFacNearSfc, 
                                  var_name='dpr_MS_corZFacNearSfc', 
                                  units=dpr_units_corZFacNearSfc,
                                  long_name=dpr_name_corZFacNearSfc)

    column.add_variable_to_object(col_MS_precipRateNearSfc, 
                                  var_name='dpr_MS_precipRateNearSfc', 
                                  units=dpr_units_precipRateNearSfc,
                                  long_name=dpr_name_precipRateNearSfc)

    column.add_variable_to_object(col_MS_precipRateAve24, 
                                  var_name='dpr_MS_precipRateAve24', 
                                  units=dpr_units_precipRateAve24,
                                  long_name=dpr_name_precipRateAve24)

    column.add_variable_to_object(col_MS_phaseNearSfc, 
                                  var_name='dpr_MS_phaseNearSfc', 
                                  units=dpr_units_phaseNearSfc,
                                  long_name=dpr_name_phaseNearSfc)

    column.add_variable_to_object(col_MS_zFacMeas, 
                                  var_name='dpr_MS_zFacMeas', 
                                  units=dpr_units_zFacMeas,
                                  long_name=dpr_name_zFacMeas)

    column.add_variable_to_object(col_MS_attenNoPrecip, 
                                  var_name='dpr_MS_attenNoPrecip', 
                                  units=dpr_units_attenNoPrecip,
                                  long_name=dpr_name_attenNoPrecip)

    column.add_variable_to_object(col_MS_corZFac, 
                                  var_name='dpr_MS_corZFac', 
                                  units=dpr_units_corZFac,
                                  long_name=dpr_name_corZFac)

    column.add_variable_to_object(col_MS_PIA_cloudwater, 
                                  var_name='dpr_MS_PIA_cloudwater', 
                                  units=dpr_units_PIA_cloudwater,
                                  long_name=dpr_name_PIA_cloudwater)

    column.add_variable_to_object(col_MS_PIA_cloudice, 
                                  var_name='dpr_MS_PIA_cloudice', 
                                  units=dpr_units_PIA_cloudice,
                                  long_name=dpr_name_PIA_cloudice)

    column.add_variable_to_object(col_MS_PIA_watervapor, 
                                  var_name='dpr_MS_PIA_watervapor', 
                                  units=dpr_units_PIA_watervapor,
                                  long_name=dpr_name_PIA_watervapor)

    column.add_variable_to_object(col_MS_PIA_oxygen, 
                                  var_name='dpr_MS_PIA_oxygen', 
                                  units=dpr_units_PIA_oxygen,
                                  long_name=dpr_name_PIA_oxygen)

    column.add_variable_to_object(col_MS_TPW_liquid, 
                                  var_name='dpr_MS_TPW_liquid', 
                                  units=dpr_units_TPW_liquid,
                                  long_name=dpr_name_TPW_liquid)

    column.add_variable_to_object(col_MS_TPW_ice, 
                                  var_name='dpr_MS_TPW_ice', 
                                  units=dpr_units_TPW_ice,
                                  long_name=dpr_name_TPW_ice)

    column.add_variable_to_object(col_MS_eff_PIA, 
                                  var_name='dpr_MS_eff_PIA', 
                                  units=dpr_units_effPIA,
                                  long_name=dpr_name_effPIA)

    #now lets add dpr_info to object
    column.add_platform_to_object(dpr_info, plat_name='dpr')
    
    return column
    #dpr_data_in_column = True


'''
Main program to read GPM Level 2 Data and supporting functions
'''

def read_2adpr_hdf5(filename):

    #open file to read
    hdf = h5py.File(filename, 'r')
    
    # get value for the FileHeader attribute, located at the top level
    # extract the individual file header values from the formatted string
    filestruc = {}
    for attr in hdf.attrs['FileHeader'].decode('UTF-8').split('\n'):
        if attr != '': filestruc[attr.split('=')[0]] = attr.split('=')[1][:-1]
    #   filestruc=parse_file_header_group(ppsFileHeaderStruc)
    #   h5a_close, fileHeaderID
    #   IF (verbose1) THEN HELP, filestruc
    prodname = filestruc['AlgorithmID']
    #   prodname=filestruc.ALGORITHMID
    
    #verify product name is 2ADPR
    if prodname != '2ADPR':
        hdf.close()
        print(f'Illegal product type {prodname}, must be 2ADPR')
        return -1
    
    #   ; define the swath groups according to product type
    snames=['HS', 'MS', 'NS']   # prefixes used in DPR products before V7
    swaths = np.array(list(hdf.keys()))[1:]
    if len(swaths) != len(snames):
        print('Expect 3 swaths in product...')
        return -1

    swathData = { }
    for swath in swaths:
       #; get the data variables for the swath groups
       #for isw = 0, filestruc.NUMBEROFSWATHS-1 do begin
       #   sname=snames[isw]
    
       #   IF N_ELEMENTS(onescan) NE 0 THEN BEGIN
       #      IF STRMATCH(onescan, sname) EQ 0b THEN BEGIN
       #         ; this isn't the scan to be read, we skip it and define empty struct
       #         print, "" & print, "Skipping swath ", sname
       #         CASE sname OF
       #            'HS' : HS = { Swath : sname+": UNREAD" }
       #            'MS' : MS = { Swath : sname+": UNREAD" }
       #            'NS' : NS = { Swath : sname+": UNREAD" }
       #         ENDCASE
       #         continue
       #      ENDIF
       #   ENDIF
    
        print(f'Swath {swath}:')
        prodgroup = prodname+'_'+swath # label info for data structures

        # get the SwathHeader for this swath
        swathstruc = {}
        for attr in hdf[swath].attrs[f'{swath}_SwathHeader'].decode('UTF-8').split('\n'):
            if attr != '': 
                if attr.split('=')[0] == 'ScanType':
                    attr_extract = attr.split('=')[1][:-1]
                else:
                    attr_extract = int(attr.split('=')[1][:-1]) # these are numbers so convert to int
                swathstruc[attr.split('=')[0]] = attr_extract
    
        # get the ScanTime info for this swath and put into dictionary
        year   = hdf[swath]['ScanTime']['Year'][:]
        month  = hdf[swath]['ScanTime']['Month'][:]
        day    = hdf[swath]['ScanTime']['DayOfMonth'][:]
        hr     = hdf[swath]['ScanTime']['Hour'][:]
        minute = hdf[swath]['ScanTime']['Minute'][:]
        second = hdf[swath]['ScanTime']['Second'][:]
        ScanTime = {'year':year,
                    'month':month,
                    'day':day,
                    'hour':hr,
                    'minute':minute,
                    'second':second}
    
        # get the scanStatus structure for this swath and put into dictionary
        dataQuality = hdf[swath]['scanStatus']['dataQuality'][:]
        FGN         = hdf[swath]['scanStatus']['FractionalGranuleNumber'][:]
        scanStatus = {'dataQuality':dataQuality, 'FGN':FGN}
        
        # get the swath-group-level datasets, put into a dictionary
        lat = hdf[swath]['Latitude'][:]
        lon = hdf[swath]['Longitude'][:]
        datasets = {'source':prodgroup, 'latitude':lat, 'longitude':lon}
    
        # get the navigation structure for this swath and put into dictionary
        navLat = hdf[swath]['navigation']['scLat'][:]
        navLon = hdf[swath]['navigation']['scLon'][:]
        navAlt = hdf[swath]['navigation']['scAlt'][:]
        navigation = {'scLat':navLat, 'scLon':navLon, 'scAlt':navAlt}
    
        # get the Experimental structure for this swath
        # will skip this since not used 
    
        # get the CSF structure for this swath and put into dictionary
        flagBB            = hdf[swath]['CSF']['flagBB'][:]
        heightBB          = hdf[swath]['CSF']['heightBB'][:]
        widthBB           = hdf[swath]['CSF']['widthBB'][:]
        qualityBB         = hdf[swath]['CSF']['qualityBB'][:]
        typePrecip        = hdf[swath]['CSF']['typePrecip'][:]
        qualityTypePrecip = hdf[swath]['CSF']['qualityTypePrecip'][:]
        csf = {'flagBB':flagBB,
               'heightBB':heightBB,
               'widthBB':widthBB,
               'qualityBB':qualityBB,
               'typePrecip':typePrecip,
               'qualityTypePrecip':qualityTypePrecip}
    
        # get the DSD structure for this swath and put into dictionary
        if 'phase' in hdf[swath]['DSD'].keys():
            phase = hdf[swath]['DSD']['phase'][:]
        else:
            phase = 'UNDEFINED'
        dsd = {'phase':phase}
    
        # get the FLG structure for this swath
        # will skip this since not used
        
        # get the PRE structure for this swath and put into dictionary
        landSurfaceType      = hdf[swath]['PRE']['landSurfaceType'][:]
        flagPrecip           = hdf[swath]['PRE']['flagPrecip'][:]
        binRealSurface       = hdf[swath]['PRE']['binRealSurface'][:]
        binStormTop          = hdf[swath]['PRE']['binStormTop'][:]
        binClutterFreeBottom = hdf[swath]['PRE']['binClutterFreeBottom'][:]
        zFactorMeasured      = hdf[swath]['PRE']['zFactorMeasured'][:]
        localZenithAngle     = hdf[swath]['PRE']['localZenithAngle'][:]
        heightStormTop       = hdf[swath]['PRE']['heightStormTop'][:]
        ellipsoidBinOffset   = hdf[swath]['PRE']['ellipsoidBinOffset'][:]
        pre = {'landSurfaceType':landSurfaceType,
               'flagPrecip':flagPrecip,
               'binRealSurface':binRealSurface,
               'binStormTop':binStormTop,
               'binClutterFreeBottom':binClutterFreeBottom,
               'zFactorMeasured':zFactorMeasured,
               'localZenithAngle':localZenithAngle, #need to transpose to match dimensions of n_pixels x n_scans -- no need to do this Apr 20, 2023 check
               'heightStormTop':heightStormTop,
               'ellipsoidBinOffset':ellipsoidBinOffset} #need to transpose to match dimensions of n_pixels x n_scans -- no need to do this Apr 20, 2023 check
    
        # get the SLV structure for this swath and put into dictionary
        piaFinal                    = hdf[swath]['SLV']['piaFinal'][:]
        zFactorCorrected            = hdf[swath]['SLV']['zFactorCorrected'][:]
        zFactorCorrectedNearSurface = hdf[swath]['SLV']['zFactorCorrectedNearSurface'][:]
        if 'precipRate' in hdf[swath]['SLV'].keys():
            precipRate                  = hdf[swath]['SLV']['precipRate'][:]
        else:
            precipRate = 'UNDEFINED'
        precipRateNearSurface       = hdf[swath]['SLV']['precipRateNearSurface'][:]
        precipWaterIntegrated       = hdf[swath]['SLV']['precipWaterIntegrated'][:]
        precipRateAve24             = hdf[swath]['SLV']['precipRateAve24'][:]
        phaseNearSurface            = hdf[swath]['SLV']['phaseNearSurface'][:]
        binEchoBottom               = hdf[swath]['SLV']['binEchoBottom'][:]
        if 'paramDSD' in hdf[swath]['SLV'].keys():
            paramDSD                    = hdf[swath]['SLV']['paramDSD'][:]
        else:
            paramDSD = 'UNDEFINED'
        slv = {'piaFinal':piaFinal,
               'zFactorCorrected':zFactorCorrected,
               'zFactorCorrectedNearSurface':zFactorCorrectedNearSurface,
               'precipRate':precipRate,
               'precipRateNearSurface':precipRateNearSurface,
               'precipWaterIntegrated':precipWaterIntegrated,
               'precipRateAve24':precipRateAve24,
               'phaseNearSurface':phaseNearSurface,
               'binEchoBottom':binEchoBottom,
               'paramDSD':paramDSD}

        # get the SRT structure for this swath
        pathAtten = hdf[swath]['SRT']['pathAtten'][:]
        srt = {'pathAtten':pathAtten}
    
        # get the VER structure for this swath
        binZeroDeg    = hdf[swath]['VER']['binZeroDeg'][:]
        heightZeroDeg = hdf[swath]['VER']['heightZeroDeg'][:]
        piaNP         = hdf[swath]['VER']['piaNP'][:]
        attenuationNP = hdf[swath]['VER']['attenuationNP'][:]
        ver = {'binZeroDeg':binZeroDeg,
               'heightZeroDeg':heightZeroDeg,
               'piaNP':piaNP,
               'attenuationNP':attenuationNP}

        swathData[swath] = {   'swath' : swath,
                                'swathHeader' : swathstruc,
                                'scantime' : ScanTime,
                                'scanStatus' : scanStatus,
                                'navigation' : navigation,
                                'CSF' : csf,
                                'DSD' : dsd,
                                'PRE' : pre,
                                'SLV' : slv,
                                'SRT' : srt,
                                'VER' : ver,
                                'datasets' : datasets }
    outData = {'FileHeader':filestruc}
    outData.update(swathData)
    #outStruc = { 'FileHeader':filestruc, 'HS':swathData['HS'], 'MS':swathData['MS'], 'NS':swathData['NS'] }
    hdf.close()
    #pdb.set_trace()

    return outData
    
def get_hgt_of_bins(nPixels, nScans, nBins, binSize, elipsBO, locZenAng):

    # need list of bins to be able to loop thru:
    binNoList = np.arange(1,nBins +1)    #bin Numbers: 1-based not 0-based
    #hts_m = np.full( (nPixels, nScans, nBins), np.nan)
    hts_m = np.full((nScans, nPixels, nBins), np.nan)
    for k, val1 in enumerate(range(nBins)): #subscripts still 0-based, of course
        for i, val2 in enumerate(range(nPixels)):
            for j, val3 in enumerate(range(nScans)):
                #HS_hts_m[i,j,k] = ( ((HS_nBins)- HS_binNoList[k] )*HS_binSize ) #sanity test
                hts_m[j,i,k] = ( ( ((nBins)- binNoList[k] )*binSize ) + elipsBO[j,i] )*np.cos(np.deg2rad(locZenAng[j,i])) 
                #pdb.set_trace()
                #test = ( ( ((HS_nBins)- HS_binNoList[k] )*HS_binSize ) + HS_elipsBO[i,j])*np.cos(np.deg2rad(HS_locZenAng[i,j]))
    return hts_m
    
def get_nscan_to_process(nScans, nPixels, longitude, latitude, limit_km, cntr_lat_deg, cntr_lon_deg):
    # for each HS/MS/NS swath, get the [pixel, scan] subset in the search area:
    
    swath_start_scan=0  ;  swath_end_scan=0  ;  swath_nscans2do=0  ; swath_start_found=0
    max_deg_lat = limit_km / 111.1
    max_deg_lon = limit_km / (np.cos(np.deg2rad(cntr_lat_deg))*111.1)
    #for s=0,HS_nScans-1 do begin
    #pdb.set_trace()
    for s, val in enumerate(range(nScans)):
        swath_found_one=0
        #for p=0,HS_nPixels-1 do begin
        for p, val2 in enumerate(range(nPixels)):
            if  abs(longitude[s,p] - cntr_lon_deg) < max_deg_lon and \
                abs(latitude[s,p] - cntr_lat_deg)  < max_deg_lat and \
                abs(latitude[s,p]) <= 90.0 and abs(longitude[s,p]) <= 180.0:
                swath_found_one = 1
                if swath_start_found == 0:
                    swath_start_found = 1
                    swath_start_scan  = s
                swath_end_scan = s    #last scan within the search range
                swath_nscans2do +=1
                break
        if swath_start_found == 1 and swath_found_one == 0: break   #end of search range
    return swath_start_scan, swath_end_scan, swath_nscans2do