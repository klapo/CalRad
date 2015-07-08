####################################################################################################
# CR.GrobsCompare.ipynb
# Karl Lapo July/2015
####################################################################################################
# Plots comparisons between ground observations and radiation products
####################################################################################################

## Import statements
# netcdf/numpy/xray
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import xray

# OS interaction
import sys, pickle, os
from sys import platform as _platform

# import subplots function for plotting
import matplotlib
# Don't let matplotlib display to the screen
matplotlib.use('Agg')
from matplotlib.pyplot import subplots
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import seaborn as sns

## Directory listing
if _platform == "linux" or _platform == "linux2":
    dir_sys = '/home/disk/p/lapok/proj/CloudClimatology/'
elif _platform == "darwin":
	dir_sys = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/'

dir_data = dir_sys+'data'
dir_print = dir_sys+'Graphics'

# List of sub-directory names for each data set
dir_NLDAS = '/NLDAS'
dir_SYN = '/CERES_SYN'
dir_grobs = '/GroundObs'
dir_VIC = '/VIC_MTCLIM'
dir_MODIS = '/MODIS.IRRAD'

# Directory for basemap pickle files
dir_bmap = dir_sys+'data/basemap'

####################################################################################################
# Functions
####################################################################################################

##### Discrete colorbar -- from Joe Hamman (https://github.com/jhamman/tonic/blob/master/tonic/plot_utils.py#L66-L94)
def cmap_discretize(cmap, n_colors=10):
    """Return discretized colormap.
    Parameters
    ----------
    cmap : str or colormap object
        Colormap to discretize.
    n_colors : int
        Number of discrete colors to divide `cmap` into.
    Returns
    ----------
    disc_cmap : LinearSegmentedColormap
        Discretized colormap.
    """
    try:
        cmap = cm.get_cmap(cmap)
    except:
        cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0, 1., n_colors), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., n_colors + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(n_colors + 1)]

    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d" % n_colors,
                                              cdict, 1024)
##### Basemap
def build_basemap(lon,lat,dir_bmap,bmap_name='basemap.pickle',rewrite=False):
    # Lat/Lon handling - map extent
    bmap_dict = {}
    bmap_dict['lat_i'] = np.min(lat)
    bmap_dict['lon_i'] = np.min(lon)
    bmap_dict['lat_j'] = np.max(lat)
    bmap_dict['lon_j'] = np.max(lon)
    
    bmap_dict['lat_mid'] = lat[np.round(lat.size/2)]
    bmap_dict['lon_mid'] = lon[np.round(lon.size/2)]
    
    bmap_dict['lat_labels'] = np.arange(np.round(bmap_dict['lat_i']), np.round(bmap_dict['lat_j']), 2)
    bmap_dict['lon_labels'] = np.arange(np.round(bmap_dict['lon_i']), np.round(bmap_dict['lon_j']), 2)
    
    os.chdir(dir_bmap)
    # Force rewriting basemap pickle file
    if rewrite:
        bmap = Basemap(llcrnrlon=bmap_dict['lon_i'],llcrnrlat=bmap_dict['lat_i'],\
                        urcrnrlon=bmap_dict['lon_j'],urcrnrlat=bmap_dict['lat_j'],\
                        rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,projection='lcc',\
                        lat_1=bmap_dict['lat_mid'],lon_0=bmap_dict['lon_mid'])
        pickle.dump(bmap,open(bmap_name,'wb'),-1)
    
    else:
        try:
            bmap = pickle.load(open(bmap_name,'rb'))
        except IOError as e:
            bmap = Basemap(llcrnrlon=bmap_dict['lon_i'],llcrnrlat=bmap_dict['lat_i'],\
                            urcrnrlon=bmap_dict['lon_j'],urcrnrlat=bmap_dict['lat_j'],\
                            rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,projection='lcc',\
                            lat_1=bmap_dict['lat_mid'],lon_0=bmap_dict['lon_mid'])
            pickle.dump(bmap,open(bmap_name,'wb'),-1)
    
    return bmap,bmap_dict

####################################
## Read previously processed data ##
####################################
# ///// See CC.CA.StatisticsMaps.Master for details on creation of xray data 

###########
## NLDAS ##
os.chdir(dir_data+dir_NLDAS)
nldas = xray.open_dataset('CA.NLDAS.irrad.monthly.nc')
nldas = nldas.rename({'DLWRF_110_SFC':'LWdwn','DSWRF_110_SFC':'SWdwn'})

#########
## SYN ##
os.chdir(dir_data+dir_SYN)
syn = xray.open_dataset('CA.SYN.irrad.monthly.nc')
syn.longitude.values = syn.longitude.values-360
syn.latitude.values = syn.latitude.values[::-1]
# Flip the syn array spatially
for d in np.arange(syn.time.size):
    syn.SWdwn.values[d-1,:,:] = np.flipud(syn.SWdwn.values[d-1,:,:])
    syn.LWdwn.values[d-1,:,:] = np.flipud(syn.LWdwn.values[d-1,:,:])

############
## MTCLIM ##
os.chdir(dir_data+dir_VIC)
mtclim = xray.open_dataset('CA.MTCLIM.irrad.monthly.nc')

#########################
## Ground Observations ##
os.chdir(dir_data+dir_grobs)
grobs = xray.open_dataset('CA.grobs.irrad.monthly.nc')
grobs.SWdwn.values[grobs.SWdwn.values == 0] = np.nan
grobs = grobs.rename({'lon':'longitude','lat':'latitude'})
grobs.longitude.values = -grobs.longitude.values

###########
## MODIS ##
os.chdir(dir_data+dir_MODIS)
modis = xray.open_dataset('CA.MODIS.irrad.monthly.nc')
modis.SWdwn.values[modis.SWdwn.values == 0] = np.nan
modis = modis.rename({'lon':'longitude','lat':'latitude'})

## List w/ all irradiance datasets
monthly_mean = {}
monthly_mean['syn'] = syn
monthly_mean['nldas'] = nldas
monthly_mean['mtclim'] = mtclim
monthly_mean['modis'] = modis
monthly_mean['grobs'] = grobs

####################################################
## Find grid point containing each ground station ##
####################################################
pr_names = ['mtclim','syn','nldas','modis']

# Station lat and lon
lon_stat = grobs.longitude.values
lat_stat = grobs.latitude.values

for pr in pr_names:        

    # lat/lon for product
    lon_rad = monthly_mean[pr].longitude.values
    lat_rad = monthly_mean[pr].latitude.values  
    # mesh
    lonm, latm = np.meshgrid(lon_rad,lat_rad)
    
    # Empty numpy array
    to_merge = np.empty((monthly_mean[pr].time.size,grobs.station.size))
    
    ## Product values in each grid containing station
    for stat in grobs.station.values:
        
        # Station index
        stat_ind = np.where(stat == grobs.station.values)
        # Distance to product grid lat-lon
        d = (latm-lat_stat[stat_ind])**2 + (lonm-lon_stat[stat_ind])**2
        # Index of closest product grid
        dind = np.where(d==np.amin(d))
        # Grad grid values at the station, put into xray dataset
        to_merge[:,stat_ind[0]] = monthly_mean[pr].SWdwn.values[:,dind[0][0],dind[1][0],np.newaxis]
    
    ## Merge products w/ grobs xray structure
    to_merge_ds = xray.Dataset({pr:(('time','station'),to_merge), \
                                    'time':monthly_mean[pr].time.values,\
                                    'station':grobs.station.values})
    grobs = grobs.merge(to_merge_ds)

#########################################
## Ground Observation Comparison Plots ##
#########################################
## product names, plotting variables, coordinates etc
# Product names
pr_names = ['grobs','mtclim','nldas','syn','modis']

# colors
SWmin_delta = -50
SWmax_delta = 50
cmap = cmap_discretize(cm.gnuplot2,15)
cmap_delta = cmap_discretize(cm.RdBu_r,11)

# Build basemap
lat = monthly_mean['mtclim'].latitude.values
lon = monthly_mean['mtclim'].longitude.values

bmp,bmd = build_basemap(lon,lat,dir_bmap,'CA.Domain.bmp.pickle',rewrite=True)
lat_labels = bmd['lat_labels']
lon_labels = bmd['lon_labels']

# Station lat and lon
lon_stat = grobs.longitude.values
lat_stat = grobs.latitude.values

## Loop through dates
for d in pd.date_range(start='2002-10-01',end='2012-10-01',freq='M'):
    print('Full domain: '+str(d))
    fig = plt.figure(figsize=(12,6))
    gs = matplotlib.gridspec.GridSpec(2,6,width_ratios=[16,16,16,16,16,1])

    ## Color range
    SWmax = 0
    SWmin = 500
    for pr in pr_names:
        if monthly_mean[pr].SWdwn.loc[d:d].any() \
            and not np.isnan(np.nanmax(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values))).any():
            SWmax = max(np.nanmax(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)),SWmax)
            SWmax = np.round(SWmax/10)*10
        
        if monthly_mean[pr].SWdwn.loc[d:d].any() \
            and np.nanmin(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)) > 0 \
            and not np.isnan(np.nanmin(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values))).any():
            SWmin = min(np.nanmin(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)),SWmin)
            SWmin = np.round(SWmin/10)*10
    dSW = 10
    
    ## Monthly averages    
    for ind,pr in enumerate(pr_names):
        ax = plt.subplot(gs[0,ind])
        
        # Lat/Lon handling - product coords
        lon_rad,lat_rad = np.meshgrid(monthly_mean[pr].longitude.values,monthly_mean[pr].latitude.values)
                
        ## Monthly value for each product
        if monthly_mean[pr].SWdwn.loc[d:d].any() and not pr == 'grobs' :
            SW_for_plot = np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)
            SW_for_plot = np.ma.masked_where(np.isnan(SW_for_plot),SW_for_plot)
            im_avg = bmp.pcolormesh(lon_rad,lat_rad,SW_for_plot,\
                        cmap=cmap,vmin=SWmin,vmax=SWmax,shading='flat',latlon=True)
            
        elif monthly_mean[pr].SWdwn.loc[d:d].any() and pr == 'grobs':
            im_avg = bmp.scatter(lon_stat,lat_stat,c=monthly_mean[pr].SWdwn.loc[d:d].values, \
                        s=75, cmap=cmap, vmin=SWmin, vmax=SWmax, linewidths=.25,latlon=True)
        ax.set_title((pr))

        ## Format
        if ind == 0:
            bmp.drawparallels(lat_labels,labels=[1,0,0,0])
        else:
            bmp.drawparallels(lat_labels)
        bmp.drawmeridians(lon_labels,labels=[0,0,0,1]) 
        
        # political boundaries.
        bmp.drawstates()
        bmp.drawcoastlines()
        bmp.drawcounties()
    
        ## Difference from ground observation values
        if not pr == 'grobs':
            ax = plt.subplot(gs[1,ind])
    
            im_dif = bmp.scatter(lon_stat,lat_stat, c=grobs[pr].loc[d:d].values-grobs.SWdwn.loc[d:d].values,\
                            s=75,cmap=cmap_delta,vmin=SWmin_delta,vmax=SWmax_delta,linewidths=.25,latlon=True)
                
            ## Format
            # Title
            ax.set_title((pr+"- ground obs"))
            # Axis
            if ind == 1:
                bmp.drawparallels(lat_labels,labels=[1,0,0,0])
            else:
                bmp.drawparallels(lat_labels)
            bmp.drawmeridians(lon_labels,labels=[0,0,0,1])
            # political boundaries.
            bmp.drawstates()
            bmp.drawcoastlines()
            bmp.drawcounties()
            
    ## Final formatting
    plt.tight_layout
    
    # Colorbar - monthly values
    caxi=plt.subplot(gs[0,-1])
    cbar=plt.colorbar(im_avg, cax=caxi, orientation = "vertical",\
                        ticks=np.arange(SWmin,SWmax+dSW,dSW),spacing='proportional')
    cbar.ax.set_ylabel(('Irradiance (Wm$^{-2}$)'))
    
    # Colorbar - difference
    caxi=plt.subplot(gs[1,-1])
    cbar = plt.colorbar(im_dif, cax=caxi, orientation = "vertical",spacing='proportional')
    cbar.ax.set_ylabel(('Difference (Wm$^{-2}$)'))
            
    fig.tight_layout()
        
    os.chdir(dir_print)
    outdate = pd.to_datetime(d) 
    outdate = outdate.strftime('%Y_%m')
    fname = 'GrObs.MonthlyDiff.'+str(outdate)+'.png'
    fig.savefig(fname)
    plt.close(fig)

##################################################################
## Ground Observation Comparison Plots -- Mountain Observations ##
##################################################################
## product names, plotting variables, coordinates etc
# Product names
pr_names = ['grobs','mtclim','nldas','syn','modis']

# colors
SWmin_delta = -50
SWmax_delta = 50

# Station lat and lon
lon_stat = grobs.longitude.values
lat_stat = grobs.latitude.values

# Lat/Lon handling - map extent
lat = np.array((36,38))
lon = np.array((-120,-118))
bmp,bmd = build_basemap(lon,lat,dir_bmap,'CAMnt.Domain.bmp.pickle',rewrite=True)
lat_labels = bmd['lat_labels']
lon_labels = bmd['lon_labels']

## Loop through dates
for d in pd.date_range(start='2002-10-01',end='2002-11-01',freq='M'):
    print('Mountain domain: '+str(d))
    fig = plt.figure(figsize=(12,6))
    gs = matplotlib.gridspec.GridSpec(2,6,width_ratios=[16,16,16,16,16,1])

    ## Color range
    SWmax = 0
    SWmin = 500
    for pr in pr_names:
        if monthly_mean[pr].SWdwn.loc[d:d].any() \
            and not np.isnan(np.nanmax(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values))).any():
            SWmax = max(np.nanmax(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)),SWmax)
            SWmax = np.round(SWmax/10)*10
        
        if monthly_mean[pr].SWdwn.loc[d:d].any() \
            and np.nanmin(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)) > 0 \
            and not np.isnan(np.nanmin(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values))).any():
            SWmin = min(np.nanmin(np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)),SWmin)
            SWmin = np.round(SWmin/10)*10
    dSW = 10
    
    ## Monthly averages    
    for ind,pr in enumerate(pr_names):
        ax = plt.subplot(gs[0,ind])
        
        # Lat/Lon handling - product coords
        lon_rad,lat_rad = np.meshgrid(monthly_mean[pr].longitude.values,monthly_mean[pr].latitude.values)
                
        ## Monthly value for each product
        if monthly_mean[pr].SWdwn.loc[d:d].any() and not pr == 'grobs' :
            SW_for_plot = np.squeeze(monthly_mean[pr].SWdwn.loc[d:d].values)
            SW_for_plot = np.ma.masked_where(np.isnan(SW_for_plot),SW_for_plot)
            im_avg = bmp.pcolormesh(lon_rad,lat_rad,SW_for_plot,\
                        cmap=cm.gnuplot2,vmin=SWmin,vmax=SWmax,shading='flat',latlon=True)
            
        elif monthly_mean[pr].SWdwn.loc[d:d].any() and pr == 'grobs':
            im_avg = bmp.scatter(lon_stat,lat_stat,c=monthly_mean[pr].SWdwn.loc[d:d].values, \
                        s=75, cmap= cm.gnuplot2, vmin=SWmin, vmax=SWmax, linewidths=.25,latlon=True)
        ax.set_title((pr))

        ## Format
        if ind == 0:
            bmp.drawparallels(lat_labels,labels=[1,0,0,0])
        else:
            bmp.drawparallels(lat_labels)
        bmp.drawmeridians(lon_labels,labels=[0,0,0,1]) 
        
        # political boundaries.
        bmp.drawstates()
        bmp.drawcoastlines()
        bmp.drawcounties()
    
        ## Difference from ground observation values
        if not pr == 'grobs':
            ax = plt.subplot(gs[1,ind])
    
            im_dif = bmp.scatter(lon_stat,lat_stat, c=grobs[pr].loc[d:d].values-grobs.SWdwn.loc[d:d].values,\
                            s=75,cmap= cm.RdBu_r,vmin=SWmin_delta,vmax=SWmax_delta,linewidths=.25,latlon=True)
                
            ## Format
            # Title
            ax.set_title((pr+"- ground obs"))
            # Axis
            if ind == 1:
                bmp.drawparallels(lat_labels,labels=[1,0,0,0])
            else:
                bmp.drawparallels(lat_labels)
            bmp.drawmeridians(lon_labels,labels=[0,0,0,1])
            # political boundaries.
            bmp.drawstates()
            bmp.drawcoastlines()
            bmp.drawcounties()
            
    ## Final formatting
    plt.tight_layout
    
    # Colorbar - monthly values
    caxi=plt.subplot(gs[0,-1])
    cbar=plt.colorbar(im_avg, cax=caxi, orientation = "vertical",\
                        ticks=np.arange(SWmin,SWmax+dSW,dSW),spacing='proportional')
    cbar.ax.set_ylabel(('Irradiance (Wm$^{-2}$)'))
    
    # Colorbar - difference
    caxi=plt.subplot(gs[1,-1])
    cbar = plt.colorbar(im_dif, cax=caxi, orientation = "vertical",spacing='proportional')
    cbar.ax.set_ylabel(('Difference (Wm$^{-2}$)'))
            
    fig.tight_layout()
        
    os.chdir(dir_print)
    outdate = pd.to_datetime(d) 
    outdate = outdate.strftime('%Y_%m')
    fname = 'GrObs_Mountain.MonthlyDiff.'+str(outdate)+'.png'
    fig.savefig(fname)
    plt.close(fig)
