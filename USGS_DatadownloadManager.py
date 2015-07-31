# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

####################################################################################################
# USGS_DatadownloadManager.ipynb
# Karl Lapo July/2015
####################################################################################################
# DEM download manager
####################################################################################################

# must insert this statement to render the plots within the notebook
# this is specific to the ipython notebook
%matplotlib inline

## Import statements
# OS interaction
import sys, pickle, os
import pandas as pd

# import subplots function for plotting
import seaborn as sns
import matplotlib
from matplotlib.pyplot import subplots
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap

## Directory listing
dir_data = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/DEM'

# Directory for basemap pickle files
dir_bmap = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/basemap'

# <codecell>

os.chdir(dir_data)
data_to_ftp = pd.read_csv('USGS.ftp_filelist.txt')
# print(data_to_ftp)
for ftp_address in data_to_ftp:
    os.system('ftp '+ftp_address)

