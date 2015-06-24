#####################################
## Shortwave Interpolation (MODIS) ##
#####################################

## Processing flags
flag_datesort = 0        # Sort MODIS data by date -- MODIS data supposed to be stored in date order, but wasn't (1 = sort)
flag_EL_build = 0        # Build instaneous elevation angle array (flag = 1) or load previously made array (flag = 0) 

## Import statements
# netcdf
from netCDF4 import Dataset,num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import xray

# OS interaction
import sys
import os
from sys import platform as _platform

# interpolation
import solargeo

## Directory of .nc files
if _platform == "linux" or _platform == "linux2":
    dir_data = '/home/disk/p/lapok/proj/CloudClimatology/data/MODIS.IRRAD/'
elif _platform == "darwin":
    dir_data = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/MODIS.IRRAD/'

## Load data
os.chdir(dir_data)
print(os.getcwd())

## Sort by date -- output w/ 'sorted' in title
if flag_datesort:
    ds_for_sort = xray.open_dataset('CA.MODIS.irrad.nc')
    
    dates = ds_for_sort.time.values
    ind = np.argsort(dates)
    dates = sorted(dates)
    time = [pd.Timestamp(t) for t in dates]
    time = np.asarray(time)
    
    # SWdwn numpy array: sorted
    SWdwn = ds_for_sort.SWdwn.values
    SWdwn = SWdwn[ind,:,:]

    # lat/lon
    lat = ds_for_sort.latitude.values
    lon = ds_for_sort.longitude.values

    ## Build xray structure and save
    ds = xray.Dataset({'SWdwn': (['time', 'lat', 'lon'], SWdwn)},coords={'time': time,'lat': lat,'lon': lon})
    ds.to_netcdf('CA.MODIS.irrad.sorted.nc')
    
## Load pre-sorted data
else:
    print('Loading sorted dataset...')
    ds = xray.open_dataset('CA.MODIS.irrad.sorted.nc')


###############################
## Instaneous Solar Geometry ##
###############################
if flag_EL_build:
    ## Datetime object (Pandas Timestamp)
    # Numpy datetime64 object
    time = ds.time.values
    # Convert to Datetime object/Pandas Timestamp (subclass of Datetime objects)
    time = [pd.to_datetime(t) for t in time]
    
    ## lat/lon & dimensions
    lat = ds.lat.values
    lon = ds.lon.values
    
    lat_size = np.shape(lat)
    lon_size = np.shape(lon)
    t_size = np.shape(time)
    
    ## Solar Zenith Angle and Azimuth Angle -- instantaneous
    # Pre-allocate
    EL = np.zeros((lon_size[0],lat_size[0],t_size[0]))
    yyyy = np.zeros(t_size[0])
    jday = np.zeros(t_size[0])
    hh = np.zeros(t_size[0])
    
    for ind,time_step in enumerate(time):
        yyyy[ind] = np.array([time_step.year])
        jday[ind] = np.array([time_step.timetuple().tm_yday])
        hh[ind] = np.array([time_step.hour+time_step.minute/60.])
    
    # mesh lat, lon, and time variables
    lat_m,lon_m,yyyy_m = np.meshgrid(lat,lon,yyyy)
    jday_m = np.meshgrid(lat,lon,jday)[2]
    hh_m = np.meshgrid(lat,lon,hh)[2]
    
    # Elevation Angle (inst.) @ all times and locations
    EL = solargeo.SUNAE(yyyy_m,jday_m,hh_m,lat_m,lon_m,refraction_flag=0)[0]
    
    # Build xray structure and save
    EL_ds = xray.Dataset({'elevation_angle': (['lon','lat','time'], EL)},coords={'time': time,'lat': lat,'lon':lon})
    EL_ds.to_netcdf('CA.MODIS.elevation_angle.nc')

else:
    print('Loading elevation angle...')
    EL_ds = xray.open_dataset('CA.MODIS.elevation_angle.nc')
    EL = EL_ds.elevation_angle.values

    # lat/lon & dimensions
    lat = EL_ds.lat.values
    lon = EL_ds.lon.values
    time = EL_ds.time.values
    
    lat_size = np.shape(lat)
    lon_size = np.shape(lon)
    t_size = np.shape(time)

####################
## Transmissivity ##
####################
## Inst. transmissivity
print('Calculating instanteous transmissivity')
# Grab inst. shortwave values -- nan out missing values
SWdwn = ds.SWdwn.values
SWdwn[SWdwn == -9999] = np.nan

# Re-shape elevation angle
EL_re = np.empty((t_size[0],lat_size[0],lon_size[0]))
for dt in np.arange(0,t_size[0]):
    EL_re[dt,:,:] = EL[:,:,dt].conj().transpose()

# Transmissivity
tau = SWdwn/(np.sin(EL_re*np.pi/180)*1365)

# Clean up (night time observations)
tau[EL_re < 0] = 0

## Checks on tau ##
flag_check_tau = 0
if flag_check_tau:
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import subplots
    from matplotlib import cm
    from mpl_toolkits.basemap import Basemap
    
    # Line plot
    fig, ax = subplots(1,1)
    ax.plot(tau[0:50,100,100])
    
    # Map
    fig, ax = subplots(1,1)
    m = Basemap(llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[-1],urcrnrlat=lat[-1],\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=lat[16],lon_0=lon[20],ax=ax)
    im = m.imshow(tau[21,:,:],cmap=cm.gnuplot2)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    m.drawcounties()
    # Meridians and Parallels
    m.drawparallels(np.arange(np.round(lat[0]), np.round(lat[-1]), .5),labels=[1,0,0,0])
    m.drawmeridians(np.arange(np.round(lon[0]), np.round(lon[-1]), 1.5),labels=[0,0,0,1])   
    cb = m.colorbar(im,"right", size="5%", pad='2%')
    cb.set_label('Elevation Angle ($\circ$)', fontsize=12)

###############################
## SH-ortwave IN-terpolation ##
###############################
## Time vector and indices
d = ds.time.values
d_pd = [pd.to_datetime(t) for t in d]
d_py = date2num(d_pd,'hours since 0001-01-01 00:00:00.0')
# day length in nanoseconds
daylen_nano = 8.64*10**13 
# Index contiguous series
ind = np.flatnonzero(np.abs(np.diff(d) > daylen_nano))
ind_discon = ((0,ind[0]),(ind[0]+1,ind[-1]),(ind[-1]+1,-1))

# Hourly time series. d_hour=datetime object (python), d_hour64=datetime64 object (numpy)
d_hour=np.empty(shape = 0)
for n in np.arange(0,3):
    # Break time series in contiguous pieces
    d_beg = d[ind_discon[n][0]]
    d_end = d[ind_discon[n][1]]
    if n == 0:
        d_hour = pd.date_range(start=d_beg,end=d_end,freq='H')
    else:
        d_hour = np.append(d_hour,pd.date_range(start=d_beg,end=d_end,freq='H'))

# Index for contiguous series
d_hour64 = pd.to_datetime(d_hour) # pandas datetime index
d_hourpy = date2num(d_hour64.to_pydatetime(),'hours since 0001-01-01 00:00:00.0') # Python datenum
ind = np.flatnonzero(np.abs(np.diff(d_hour64)) > daylen_nano)
ind_discon_hour = ((0,ind[0]),(ind[0]+1,ind[-1]),(ind[-1]+1,-1))

## Elevation Angle (hourly) and interpolated tau

# Pre-allocation
nlat = lat.shape[0]
nlon = lon.shape[0]
ntime = d_hour.shape[0]
el_hour = np.empty((ntime,nlat,nlon))
tau_interp = np.empty((ntime,nlat,nlon))
print('Calculating hourly elevation angle and interpolating transmissivity...')
# Loop through each period
for n in np.arange(0,3):
    # Break time series in contiguous pieces
    d_cont = d_hour64[ind_discon_hour[n][0]:ind_discon_hour[n][1]]
    d_contpy = date2num(d_cont.to_pydatetime(),'hours since 0001-01-01 00:00:00.0')
    
    # Nested for loops! Cause that's an efficient way to operate
    # Loop over all lat/lon
    for lat_i in lat:
        for lon_i in lon:
            ind_lat = np.flatnonzero(lat == lat_i)
            ind_lon = np.flatnonzero(lon == lon_i)
            
            print('Working on lat index: '+ str(ind_lat)+'/'+str(lat.size))
            print('Working on lon index: '+ str(ind_lon)+'/'+str(lon.size))

            # Skip areas outside the study domain
            if np.all(np.isnan(SWdwn[:,ind_lat,ind_lon])):
                continue
            
            ## Average elevation angle for each point
            el_hour[ind_discon_hour[n][0]:ind_discon_hour[n][1],ind_lat[0],ind_lon[0]] = \
                np.squeeze(solargeo.AVG_EL(d_cont,lat_i,lon_i,8,'MID'))
                
            ## Interpolate transmissivity
            tau_interp[ind_discon_hour[n][0]:ind_discon_hour[n][1],ind_lat[0],ind_lon[0]] = np.interp(\
                        d_contpy,d_py[ind_discon[n][0]:ind_discon[n][1]],\
                        tau[ind_discon[n][0]:ind_discon[n][1],ind_lat[0],ind_lon[0]])

## Checks on interpolated tau
flag_check_tau_interp = 0
if flag_check_tau_interp:
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import subplots
    from matplotlib import cm
    from mpl_toolkits.basemap import Basemap
    
    # Line plot
    fig, ax = subplots(1,1)
    ax.plot(tau_interp[0:50,100,100])
    
    # Map
    fig, ax = subplots(1,1)
    m = Basemap(llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[-1],urcrnrlat=lat[-1],\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=lat[16],lon_0=lon[20],ax=ax)
    im = m.imshow(tau_interp[21,:,:],cmap=cm.gnuplot2)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    m.drawcounties()
    # Meridians and Parallels
    m.drawparallels(np.arange(np.round(lat[0]), np.round(lat[-1]), .5),labels=[1,0,0,0])
    m.drawmeridians(np.arange(np.round(lon[0]), np.round(lon[-1]), 1.5),labels=[0,0,0,1])   
    cb = m.colorbar(im,"right", size="5%", pad='2%')
    cb.set_label('Elevation Angle ($\circ$)', fontsize=12)

## Output intermediary steps
el_hour_ds = xray.Dataset({'elevation_angle': (['time','lat','lon'], el_hour)},coords={'time': d_hour64,'lat':lat,'lon':lon})
el_hour_ds.to_netcdf('CA.MODIS.elev_hour.nc')
tau_hour_ds = xray.Dataset({'transmissivity': (['time','lat','lon'], tau_interp)},coords={'time': d_hour64,'lat':lat,'lon':lon})
tau_hour_ds.to_netcdf('CA.MODIS.transmissivity_hour.nc')

## Interpolated shortwave
# Have average elevation angle for an hour -> average shortwave for an hour
# Step may need to be revisited, order of operations changed from matlab functions
SWdwn_hour = tau_interp*1365*el_hour
# Build xray structure and save
SWdwn_hour_ds = xray.Dataset({'SWdwn': (['time','lat','lon'], SWdwn_hour)},coords={'time': d_hour64,'lat': lat,'lon':lon})
SWdwn_hour_ds.to_netcdf('CA.MODIS.interp_shortwave_hour.nc')

