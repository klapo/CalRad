import cython
import numpy as np

def SUNAE( YEAR, DAY, HOUR, LAT, LONG, refraction_flag=1):
#--------------------------------------------------------------------
## DESCRIPTION ##
#    Calculates the local solar azimuth and elevation angles, and
#    the distance to and angle subtended by the Sun, at a specifi#
#    location and time using approximate formulas in The Astronomical
#    Almanac.  Accuracy of angles is 0.01 deg or better (the angular
#    width of the Sun is about 0.5 deg, so 0.01 deg is more than
#    sufficient for most applications).
#
#    Unlike many GCM (and other) sun angle routines, this
#    one gives slightly different sun angles depending on
#    the year.  The difference is usually down in the 4th
#    significant digit but can slowly creep up to the 3rd
#    significant digit after several decades to a century.
#
#    A refraction correction appropriate for the "US Standard
#    Atmosphere" is added, so that the returned sun position is
#    the APPARENT one.  The correction is below 0.1 deg for solar
#    elevations above 9 deg.  This refraction correction is assumed on,
#    but can be toggled using the "refraction_flag" keyword.
#
#    Only accurate between 1950 and 2050.
#--------------------------------------------------------------------

#--------------------------------------------------------------------
##  REFERENCES ##
#    Michalsky, J., 1988: The Astronomical Almanac's algorithm for
#       approximate solar position (1950-2050), Solar Energy 40,
#       227-235 (but the version of this program in the Appendix
#       contains errors and should not be used)
#    The Astronomical Almanac, U.S. Gov't Printing Office, Washington,
#       D.C. (published every year): the formulas used from the 1995
#       version are as follows:
#       p. A12: approximation to sunrise/set times
#       p. B61: solar elevation ("altitude") and azimuth
#       p. B62: refraction correction
#       p. C24: mean longitude, mean anomaly, eclipti#longitude,
#               obliquity of ecliptic, right ascension, declination,
#               Earth-Sun distance, angular diameter of Sun
#       p. L2:  Greenwich mean sidereal time (ignoring T^2, T^3 terms)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
## HISTORY ##
#  (fortran) Authors:  Dr. Joe Michalsky (joe@asrc.albany.edu)
#            Dr. Lee Harrison (lee@asrc.albany.edu)
#            Atmospheri#Sciences Research Center
#            State University of New York
#            Albany, New York
#  (fortran) Modified by:  Dr. Warren Wiscombe (wiscombe@climate.gsfc.nasa.gov)
#                NASA Goddard Space Flight Center
#                Code 913
#                Greenbelt, MD 20771
#  (python)  Converted to python: Karl Lapo (lapo.karl@gmail.com)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
## INPUT/OUTPUT ##
#  Input:
#     YEAR     year (INTEGER; range 1950 to 2050)
#     DAY      day of year at LAT-LONG location (INTEGER; range 1-366)
#     HOUR     hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0)
#              = (local hour) + (time zone number)
#                + (Daylight Savings Time correction; -1 or 0)
#              (local hour) range is 0 to 24,
#              (time zone number) range is -12 to +12, and
#              (Daylight Time correction) is -1 if on Daylight Time
#              (summer half of year)
#              Example: 8:30 am Eastern Daylight Time would be
#
#                          HOUR = 8.5 + 5 - 1 = 12.5
#     LAT      latitude [degrees] - north is positive)
#     LONG     longitude [degrees] - east is positive)
#  Output:
#     AZ       solar azimuth angle (measured east from north, 0 to 360 degs)
#     EL       solar elevation angle (angle above the horizon)
#     SOLDST   distance to sun [Astronomical Units, AU]
#              (1 AU = mean Earth-sun distance = 1.49597871E+11 m)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#  Local Variables:
#    DEC       Declination (radians)
#    ECLONG    Eclipti#longitude (radians)
#    GMST      Greenwich mean sidereal time (hours)
#    HA        Hour angle (radians, -pi to pi)
#    JD        Modified Julian date (number of days, including
#              fractions thereof, from Julian year J2000.0);
#              actual Julian date is JD + 2451545.0
#    LMST      Local mean sidereal time (radians)
#    MNANOM    Mean anomaly (radians, normalized to 0 to 2*pi)
#    MNLONG    Mean longitude of Sun, corrected for aberration
#              (deg; normalized to 0-360)
#    OBLQEC    Obliquity of the eclipti#(radians)
#    RA        Right ascension  (radians)
#    REFRAC    Refraction correction for US Standard Atmosphere (degs)

#--------------------------------------------------------------------
## Type declarations, keeping around for future speed up with cython
#    .. Scalar Arguments ..
#      INTEGER   YEAR, DAY
#      REAL      AZ, EL, HOUR, LAT, LONG, SOLDIA, SOLDST
#    ..
#    .. Local Scalars ..
#      LOGICAL   PASS1
#      INTEGER   DELTA, LEAP
#      DOUBLE PRECISION  DEC, DEN, ECLONG, GMST, HA, JD, LMST,
#     &                  MNANOM, MNLONG, NUM, OBLQEC, PI, RA,
#     &                  RPD, REFRAC, TIME, TWOPI
#    ..
#    .. Intrinsi#Functions ..
#      INTRINSIC AINT, ASIN, ATAN, COS, MOD, SIN, TAN
#    ..
#    .. Data statements ..
#
#      SAVE     PASS1, PI, TWOPI, RPD
#      DATA     PASS1 /.True./
#    ..

#===================================================================
## Error handling
    if np.min(YEAR)<1950 or np.max(YEAR)>2050:
        raise ValueError('YEAR must be between 1950 and 2050')
    if np.min(DAY)<1 or np.max(DAY)>366:
        raise ValueError('DAY must be between 1 and 366')
    if np.min(HOUR)<-13.0 or np.max(HOUR)>36.0:
        raise ValueError('HOUR must be between -13 and 36')
    if np.min(LAT)<-90.0 or np.max(LAT)>90.0:
        raise ValueError('LAT must be between -90 and 90')
    if np.min(LONG)<-180.0 or np.max(LONG)>180.0:
        raise ValueError('LONG must be between -180 and 180')

#===================================================================
## Julian date/Coordinates 
    RPD = np.pi/180.
    PI = np.pi
    TWOPI = np.pi*2
    # Add 2,400,000 for true JD
    # LEAP = leap days since 1949
    # 32916.5 is midnite 0 jan 1949 minus 2.4e6
    DELTA = YEAR - 1949
    LEAP = DELTA/4
    JD = 32916.5 + (DELTA*365 + LEAP + DAY) + HOUR / 24.

# ecliptic coordinates: 51545.0 + 2.4e6 = noon 1 jan 2000
    TIME  = JD - 51545.0

# force mean longitude between 0 and 360 degs
    MNLONG = 280.460 + 0.9856474*TIME
    MNLONG = np.mod(MNLONG,360.)
    if MNLONG < 0:
        MNLONG = MNLONG + 360.
# mean anomaly in radians between 0 and 2*pi
    MNANOM = 357.528 + 0.9856003*TIME
    MNANOM = np.mod( MNANOM, 360. )
    if MNANOM < 0.:
        MNANOM = MNANOM + 360.
    MNANOM = MNANOM*RPD
# ecliptic longitude and obliquity of ecliptic in radians
    ECLONG = MNLONG + 1.915*np.sin(MNANOM) + 0.020*np.sin(2.*MNANOM)
    ECLONG = np.mod(ECLONG, 360.)
    if ECLONG < 0.:
        ECLONG = ECLONG + 360.
    OBLQEC = 23.439 - 0.0000004*TIME
    ECLONG = ECLONG*RPD
    OBLQEC = OBLQEC*RPD

#===================================================================
## EL, AZ
# right ascension
    NUM  = np.cos(OBLQEC)*np.sin(ECLONG)
    DEN  = np.cos(ECLONG)
    RA   = np.arctan(NUM/DEN)
# Force right ascension between 0 and 2*pi
    if DEN<0.0:
        RA = RA + PI
    elif NUM<0.0:
        RA = RA + TWOPI
# declination
    DEC = np.arcsin(np.sin(OBLQEC)*np.sin(ECLONG))
# Greenwich mean sidereal time in hours
    GMST = 6.697375 + 0.0657098242*TIME + HOUR
# Hour not changed to sidereal time since 
# 'time' includes the fractional day
    GMST  = np.mod(GMST,24.)
    if GMST<0.:
        GMST = GMST + 24.
# local mean sidereal time in radians
    LMST  = GMST + LONG / 15.
    LMST  = np.mod( LMST, 24. )
    if LMST<0.:
        LMST   = LMST + 24.
    LMST   = LMST*15.*RPD
# hour angle in radians between -pi and pi
    HA  = LMST - RA
    if HA<- PI:
            HA  = HA + TWOPI
    if HA>PI:
        HA  = HA - TWOPI
# solar azimuth and elevation
    EL  = np.arcsin(np.sin(DEC)*np.sin(LAT*RPD) + np.cos(DEC)*np.cos(LAT*RPD)*np.cos(HA))
    AZ  = np.arcsin(-np.cos(DEC)*np.sin(HA)/np.cos(EL))
# Put azimuth between 0 and 2*pi radians
    if np.sin(DEC)-np.sin(EL)*np.sin(LAT*RPD)>=0.:
        if np.sin(AZ)<0.:
            AZ  = AZ + TWOPI
    else:
        AZ  = PI - AZ
# Convert elevation and azimuth to degrees
    EL  = EL / RPD
    AZ  = AZ / RPD

#===================================================================
## Refraction correction for U.S. Standard Atmos.
#  (assumes elevation in degs) (3.51823=1013.25 mb/288 K)
    if refraction_flag:
        if EL>=19.225:
            REFRAC = 0.00452*3.51823 / np.arctan(EL*RPD)
        elif EL>-0.766 and EL<19.225 :
            REFRAC = 3.51823 * (0.1594 + EL*(0.0196 + 0.00002*EL) ) / (1. + EL*(0.505 + 0.0845*EL))
        elif EL<=-0.766 :
            REFRAC = 0.0
        EL  = EL + REFRAC

#===================================================================
## SOLDST, distance to sun in A.U.
    SOLDST = 1.00014 - 0.01671*np.cos(MNANOM) - 0.00014*np.cos(2.*MNANOM)
    
    if EL<-90.0 or EL>90.0:
        raise ValueError('Calculated EL out of range')
    if AZ<0.0 or AZ>360.0:
        raise ValueError('Calculated AZ out of range')

#===================================================================
# Complete
    return EL, AZ, SOLDST

#===================================================================
import solargeometry
import numpy as np
def AVG_EL(TIME,lat,lon,REF):
#===================================================================
## Description ##
# Calculates the average cosSZA during the time interval and returns the
# effective elevation angle. For instantaneous EL values call SUNAE directly.
# Converts timestamp/labels to beginning of bin.
#
# THIS CODE DOES NOT HANDLE DISCONTINUOUS DATA (YET)
#
# SYNTAX:
#   EL = AVG_EL(TIME,lat,lon,tz,REF)
#
# INPUTS:
#   TIME    = Nx1 datetime object (assumed to be in UTC, conversion should occur outside of function)
#   lat     = 1x1 scalar, site latitude (north = pos)
#   lon     = 1x1 scalar, site longitude (east = pos)
#   REF     = string - argument describing how the data is referenced to the time stamp
#           The default assumption is that the time stamp is for the beginning of the 
#           averaging interval. This argument must be specified for accuracy if the
#           default value is not true.
#           'END' - averaged data referenced to the interval end
#           'MID' - averaged data referenced to the interval middle
#           'BEG' - averaged data referenced to the interval beginning
#
# OUTPUTS:
#   EL      = Nx1 vector - Average elevation angle [degrees above horizon]
#
# DEPENDENCIES:
#   SUNAE

#===================================================================
## Libraries
    import datetime
    from datetime import datetime, timedelta
    import pytz
    import pandas as pd
    
#===================================================================
# Time stamp referenced moved to the beginning of the averaging period
    dt = TIME[1]-TIME[0] # Time step (timedelta object)
    if REF=='END':
        TIME = TIME - dt
    elif REF=='MID':
        TIME = TIME - dt/2    
    elif REF!='BEG' & REF!='MID' & REF!='END':
        raise ValueError('Unrecognized REF option')
    
    # Instantaneous elevation angle -> 5 minute time step integration
    # Discretize current time step
    t_fine_beg = TIME[0]
    t_fine_end = TIME[-1]
    t_fine = pd.date_range(start=t_fine_beg,end=t_fine_end,freq='5Min')
        
#===================================================================
# Numerical integration of fine EL
    yyyy = t_fine.year
    jday = t_fine.dayofyear
    hh = t_fine.hour
    mm = t_fine.minute
    EL_fine = [solargeometry.SUNAE(YEAR,DAY,HOUR,MINUTE,lat,lon)[0] for YEAR,DAY,HOUR,MINUTE in zip(yyyy,jday,hh,mm)]
    EL_fine = np.array(EL_fine)
    mew_fine = pd.DataFrame(data=np.sin(EL_fine*np.pi/180.), index=t_fine, columns=['cos_SZA'])
    
    # Average elevation angle at original timestep
    resample_rule_secs = str(dt.seconds)+'S'
    mew_coarse = mew_fine.resample(resample_rule_secs,how='mean',closed='right',label='left')
    EL = np.arcsin(mew_coarse['cos_SZA'])*180./np.pi


#===================================================================
## Finish ##
    return(EL)
