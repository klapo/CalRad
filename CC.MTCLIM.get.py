# Downloads Livneh 2013 data set from ftp server over the CA domain
import os, sys

ftp_address = 'ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Meteorology.asc.v.1.2.1915.2011.bz2/'
lon_address = ['120.115','125.120']
lat_address = ['25.36','37.49']

data_dir = '/home/lapok/Livneh_VIC'
os.chdir(data_dir)
for lon in lon_address:
	for lat in lat_address:
		# Grab all the ascii files within each lat/lon sub-directory
		directory_for_retr = ftp_address+'data.'+lon+"."+lat+"/"
		print(directory_for_retr)
		os.system("wget -r "+directory_for_retr)

