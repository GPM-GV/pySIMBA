import xarray as xr 
import h5py
import numpy as np
import pandas as pd
import warnings
import pdb
warnings.filterwarnings('ignore')

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