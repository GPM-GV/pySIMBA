#This program will be used for plotting routines for pySIMBA output.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

#define global variables
fontsize=12

class PiecewiseNorm(Normalize):
    def __init__(self, levels, clip=False):
        # the input levels
        self._levels = np.sort(levels)
        # corresponding normalized values between 0 and 1
        self._normed = np.linspace(0, 1, len(levels))
        Normalize.__init__(self, None, None, clip)

    def __call__(self, value, clip=None):
        # linearly interpolate to get the normalized value
        return np.ma.masked_array(np.interp(value, self._levels, self._normed))

def column_grid_reflectivity(x, y, twodvd_data, apu_data, radar_data, case_date, levels = None, png=False):

    x = x/1000. #convert to km
    y = y/1000. #convert to km
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    #1) plot disdrometer data
    for i, xx in enumerate(x):
        for j, yy in enumerate(y):
            check_twodvd = twodvd_data[i,j]
            check_apu = apu_data[i,j]
            
            #plot both sets together
            if not np.isnan(check_twodvd) and not np.isnan(check_apu):
                ax.plot(xx, yy, marker='o', color='k', markersize=5)
                ax.annotate(f'{check_twodvd: .1f}', (xx+0.05, yy+0.05), size=fontsize*1.8)
                
                ax.plot(xx+0.10, yy-0.10, marker='o', color='k', markersize=5)
                ax.annotate(f'{check_apu: .1f}', (xx+0.10+0.02, yy-0.10-0.08), size=fontsize*1.8)
            #plot if either one is available
            else:
                if not np.isnan(check_twodvd):
                    ax.plot(xx, yy, marker='o', color='k', markersize=5)
                    ax.annotate(f'{check_twodvd: .1f}', (xx+0.02, yy+0.02), size=fontsize*1.8)
                
                if not np.isnan(check_apu):
                    ax.plot(xx, yy, marker='o', color='k', markersize=5)
                    ax.annotate(f'{check_apu: .1f}', (xx+0.02, yy+0.02), size=fontsize*1.8)
                
    #2) plot radar data
    if levels is None:
        levels=[6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45]
    fill=ax.contourf(x, y, radar_data, levels, cmap='jet', norm=PiecewiseNorm(levels),extend='both')
    
    cbaxes = fig.add_axes([0.92,0.085,0.03,0.81])
    cbar = plt.colorbar(fill, ticks=levels, cax=cbaxes)
    cbar.set_label('Reflectivity [dBZ]', fontsize=fontsize*1.8)
    cbar.ax.tick_params(labelsize=fontsize*1.5)
    ax.set_xlabel('Distance East of Grid Center (km)', fontsize=fontsize*2.0)
    ax.set_ylabel('Distance North of Grid Center (km)', fontsize=fontsize*2.0)
    ax.xaxis.grid(True, zorder=0, color='gray')
    ax.yaxis.grid(True, zorder=0, color='gray')
    ax.tick_params(labelsize=fontsize*2.0)
    ax.tick_params(which='both', width=2)
    
    if png:
        path = os.getcwd()
        png_file=f'{path}/NPOL_Disdrometer_Reflectivity_{case_date}.png'
        fig.savefig(png_file, bbox_inches='tight', dpi=dpi)
        print('--> '+png_file)
    else:
        plt.show()
        
def column_grid_rain(x, y, twodvd_data, apu_data, radar_data, case_date, levels = None, png=False):

    x = x/1000. #convert to km
    y = y/1000. #convert to km
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    #1) plot disdrometer data
    for i, xx in enumerate(x):
        for j, yy in enumerate(y):
            check_twodvd = twodvd_data[i,j]
            check_apu = apu_data[i,j]
            
            #plot both sets together
            if not np.isnan(check_twodvd) and not np.isnan(check_apu):
                ax.plot(xx, yy, marker='o', color='k', markersize=5)
                ax.annotate(f'{check_twodvd: .1f}', (xx+0.05, yy+0.05), size=fontsize*1.8)
                
                ax.plot(xx+0.10, yy-0.10, marker='o', color='k', markersize=5)
                ax.annotate(f'{check_apu: .1f}', (xx+0.10+0.02, yy-0.10-0.08), size=fontsize*1.8)
            #plot if either one is available
            else:
                if not np.isnan(check_twodvd):
                    ax.plot(xx, yy, marker='o', color='k', markersize=5)
                    ax.annotate(f'{check_twodvd: .1f}', (xx+0.02, yy+0.02), size=fontsize*1.8)
                
                if not np.isnan(check_apu):
                    ax.plot(xx, yy, marker='o', color='k', markersize=5)
                    ax.annotate(f'{check_apu: .1f}', (xx+0.02, yy+0.02), size=fontsize*1.8)
                
    #2) plot radar data
    if levels is None:
        levels=[0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10]
    fill=ax.contourf(x, y, radar_data, levels, cmap='rainbow', norm=PiecewiseNorm(levels),extend='both')
    
    cbaxes = fig.add_axes([0.92,0.085,0.03,0.81])
    cbar = plt.colorbar(fill, ticks=levels, cax=cbaxes)
    cbar.set_label('Rain Rate [mm hr^-1]', fontsize=fontsize*1.8)
    cbar.ax.tick_params(labelsize=fontsize*1.5)
    ax.set_xlabel('Distance East of Grid Center (km)', fontsize=fontsize*2.0)
    ax.set_ylabel('Distance North of Grid Center (km)', fontsize=fontsize*2.0)
    ax.xaxis.grid(True, zorder=0, color='gray')
    ax.yaxis.grid(True, zorder=0, color='gray')
    ax.tick_params(labelsize=fontsize*2.0)
    ax.tick_params(which='both', width=2)
    
    if png:
        path = os.getcwd()
        png_file=f'{path}/NPOL_Disdrometer_RainRate_{case_date}.png'
        fig.savefig(png_file, bbox_inches='tight', dpi=dpi)
        print('--> '+png_file)
    else:
        plt.show()