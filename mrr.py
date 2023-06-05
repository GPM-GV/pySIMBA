import pandas as pd
import numpy as np
import pdb
import math
from scipy import optimize

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
def read_data(file_mrr):

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
