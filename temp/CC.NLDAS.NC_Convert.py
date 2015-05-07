import os,sys

# Need write directory

# use the "os.walk" function to walk down through all NLDAS files

# Requires argument for directory containing grib files
grib_dir = sys.argv[1]

# Directory contents
os.chdir(grib_dir)
content = os.listdir(os.getcwd())

for files in content:
    if files[-3:] == 'grb':
        os.system('ncl_convert2nc '+files+' -B')



