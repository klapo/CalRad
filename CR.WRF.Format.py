####################################################################################################
# CR.WRF.Format.ipynb
# Karl Lapo July/2015
####################################################################################################
# Formats WRF for use in CalRad project--cuts out domain, retrieves relevant variables
####################################################################################################
## Import statements
import gc

# netcdf
from netCDF4 import Dataset,num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import xray

# OS interaction
import sys
import os

# OS interaction
import sys
import os

## Directory information
# dir_rawnc = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/WRF'
# dir_cleannc = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/WRF'
# dir_procnc = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/WRF'

# dir_rawnc = '/storage/WRF/newWRFRD1'
# dir_cleannc = '/home/lapok/WRF/clean'
# dir_procnc = '/home/lapok/WRF/proc'

dir_rawnc = '/home/disk/p/lapok/proj/CloudClimatology/data/WRF/clean'
dir_cleannc = '/home/disk/p/lapok/proj/CloudClimatology/data/WRF/clean'
dir_procnc = '/home/disk/p/lapok/proj/CloudClimatology/data/WRF/proc'

####################################################################################################
# Script
####################################################################################################
#### netcdf parameters
# Input file name structure
nc_name = 'radoutd02_newWRFRD_'
# Output file names
fname_CA_domain = 'CA.WRF.irrad.nc'
clean_prepend = 'CA.WRF.irrad.clean.'
# Fields to read
SYN_netcdf_fields = ['lwdnb', 'ter', 'swdnb']

# Set initilization flag for reading variables
init_flag = 1
# Bounding box - rectangular domain
LL_rect = [-123.5,34.5]
LR_rect = [-115.5,34.5]
UR_rect = [-115.5,41]
UL_rect = [-123.5,41]
# Bounding box - ragged domain
LL_rag = [-120,34.5]
LR_rag = [-115,34.5]
UR_rag = [-118.5,41]
UL_rag = [-123.5,41]

#### Domain indices and coordinates
# Directory contents
os.chdir(dir_rawnc)
content = os.listdir(os.getcwd())
for files in content:
    if files[0:-10] == nc_name and files[-3:] == '.nc':
        netcdf_file = files
        ## Initialize - Grab first file, get details for pre-allocation
        print('Pre-allocating....   \n')
    
        # Need to somehow read the created netcdf file
        f = Dataset(netcdf_file,'r')
        
        ## Rectangular bounding box for CA study area
        lat = f.variables['lat'][:]
        lon = f.variables['lon'][:]
        # Extract only within rectangular domain - use this index variable for data extraction from netcdf
        ind_lat_rect,ind_lon_rect = np.nonzero( (lat > LR_rect[1]) & (lat < UL_rect[1]) & \
                                               (lon < LR_rect[0]) & (lon > LL_rect[0]) )        
        # lat/lon as series
        lat_ser = lat[:,0]
        lon_ser = lon[0,:]
        ## lat/lon for output, in mesh format
        ind_lat_s_rect = np.flatnonzero( (lat_ser > LR_rect[1]) & (lat_ser < UL_rect[1]) )
        ind_lon_s_rect = np.flatnonzero( (lon_ser < LR_rect[0]) & (lon_ser > LL_rect[0]) )
        lat_m,lon_m = np.meshgrid(lat_ser[ind_lat_s_rect],lon_ser[ind_lon_s_rect])
        ## Ragged domain, CA study area - indices for NaN'ing data outside 
        line_west_m = (UL_rag[1]-LL_rag[1])/(UL_rag[0]-LL_rag[0])
        line_west_b = LL_rag[1]-line_west_m*LL_rag[0]
        line_east_m = (UR_rag[1]-LR_rag[1])/(UR_rag[0]-LR_rag[0])
        line_east_b = LR_rag[1]-line_east_m*LR_rag[0]
        ind_lon_rag,ind_lat_rag = np.nonzero((lon_m < (lat_m-line_west_b)/line_west_m) | \
                                             (lon_m > (lat_m-line_east_b)/line_east_m) | \
                                             (lat_m < LR_rag[1]) | (lat_m > UL_rag[1]))
        # serial lat/lon for xray structure dimensions
        lat_out = lat_ser[ind_lat_s_rect]
        lon_out = lon_ser[ind_lon_s_rect]
        f.close()
        break

# #### clean nc files, reduce to only relevant variables
# # Directory contents
# os.chdir(dir_rawnc)
# content = os.listdir(os.getcwd())
# for files in content:
#     os.chdir(dir_rawnc)
#     if files[0:-10] == nc_name and files[-3:] == '.nc':
#         netcdf_file = files
#         print("Cleaning: "+netcdf_file)
    
#         ## Read netcdf file
#         ds_in = xray.open_dataset(netcdf_file,chunks={'Time':5})
#         # Drop irrelevant variables
#         ds_in = ds_in.drop(['ncl10','ncl11','ncl12','ncl4','ncl5','ncl6','ncl7','ncl8','ncl9',\
#                             'lu','lm','qfx','emiss','albedo','tsk','lat','lon'])
#         # Datetime
#         wrfd = ds_in.Times.values
#         wrfdates = []
#         [wrfdates.append(datetime.strptime(d,"%Y-%m-%d_%H:%M:%S")) for d in wrfd]
                
#         ## Build output structure
#         SWdwn = ds_in.swdnb.values[:,ind_lat_s_rect.min():ind_lat_s_rect.max()+1,ind_lon_s_rect.min():ind_lon_s_rect.max()+1]
#         LWdwn = ds_in.lwdnb.values[:,ind_lat_s_rect.min():ind_lat_s_rect.max()+1,ind_lon_s_rect.min():ind_lon_s_rect.max()+1]
#         ter = ds_in.ter.values[ind_lat_s_rect.min():ind_lat_s_rect.max()+1,ind_lon_s_rect.min():ind_lon_s_rect.max()+1]
#         ds_out = xray.Dataset({'SWdwn':(['time','lat','lon'],SWdwn),\
#                                 'LWdwn':(['time','lat','lon'],LWdwn),\
#                                 'ter':(['lat','lon'],ter)},\
#                                 coords={'time': (['time'],wrfdates),\
#                                         'lat':(['lat'],lat_out),\
#                                         'lon':(['lon'],lon_out)})
        
#         ## Write
#         os.chdir(dir_cleannc)
#         fname_out = clean_prepend+files[-10:-3]+'.nc'
#         ds_out.to_netcdf(fname_out)
        
#         ## Finish, close files, clear out memory
#         wrfd = None
#         wrfdates = None
#         ds_in = None
#         ds_out = None
#         SWdwn = None
#         LWdwn = None
#         ter = None
#         collect = gc.collect()
        
#### Read netcdfs 
# Directory contents
os.chdir(dir_cleannc)
content = os.listdir(os.getcwd())
init_flag = 1
iter_count = 0
for files in content:
    if files[0:19] == clean_prepend and files[-3:] == '.nc':
        netcdf_file = files
        print("Processing: "+netcdf_file)
    
        ## Read netcdf file
        ds_in = xray.open_dataset(netcdf_file)
        ds_in = ds_in.drop('ter')
        
        ## NaN data outside CA domain
        ds_in.SWdwn.values[:,ind_lat_rag,ind_lon_rag] = np.nan
        ds_in.LWdwn.values[:,ind_lat_rag,ind_lon_rag] = np.nan
        ds_in = ds_in.chunk({'time': 1})
          
        ## xray structure/concat
        if init_flag == 1:
            wrf = ds_in
            init_flag = 0
        else:
            wrf = xray.concat([wrf,ds_in],dim='time')
        iter_count = iter_count+1
        wrf = wrf.chunk({'time': iter_count})
    
        ds_in = None
        collect = gc.collect()
        
#### Elevation of WRF grid
for files in content:
    if files[0:19] == clean_prepend and files[-3:] == '.nc':
        netcdf_file = files
        print('Grabbing elevation....   \n')
    
        ## Read netcdf file
        ds_in = xray.open_dataset(netcdf_file)
        
        ## NaN data outside CA domain
        wrf['elev'] = ds_in.ter
        wrf.elev.values[ind_lat_rag,ind_lon_rag] = np.nan
        break

#### Write netcdf for hourly, daily, and monthly values
os.chdir(dir_procnc)
wrf.to_netcdf(fname_CA_domain)
wrf_day = wrf.resample(freq='D', dim='time', how='mean')
wrf_day.to_netcdf('CA.WRF.irrad.daily.nc')
wrf_month = wrf.resample(freq='M', dim='time', how='mean')
wrf_month.to_netcdf('CA.WRF.irrad.monthly.nc')

