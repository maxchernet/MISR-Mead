# Calculation of Sun Zenith/Azimuth angles with given lat/lon, date, time and UTC time zone
#
# The calculations in the NOAA Sunrise/Sunset and Solar Position Calculators are based on equations from Astronomical
# Algorithms by Jean Meeus. The sunrise and sunset results are theoretically accurate to within a minute for locations
# between +/- 72 degree latitude,
# and within 10 minutes outside of those latitudes. However, due to variations in atmospheric composition, temperature,
# pressure and
# conditions, observed values may vary from calculations.
#
# theta_i - Zenith
# phi_i - Azimuth
#
# Translated to Python from NOAA_Solar_Calculations_year.xls from
# http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html by Max Chernetskiy

import numpy as np
import julday

def degrees(rad):
    return rad*(180./np.pi)

def radians(deg):
    return deg*(np.pi/180.)

def CalcSun(lat, lon, date, time=12./24., zone=-6):

    #zone - time zone (+ to E)
    #jd = 2452276.25
    #sun_dec - Sun Declin (deg)
    #hour_angle - Hour Angle (deg)
    #tst - True Solar Time (min)
    #et - Eq of time (min)
    #var_y
    #gml - Geom Mean Long Sun (deg)
    #eeo - Eccent Earth Orbit
    #gma - Geom Mean Anom Sun (deg)
    #jc - Julian Century
    #oc - Obliq Corr (deg)
    #moe - Mean Obliq Ecliptic (deg)
    #sal - Sun App Long (deg)
    #stl - Sun True Long (deg)
    #sec - Sun Eq of Ctr

    jd = np.round(julday.date2jul(date))
    jc = (jd-2451545.)/36525.
    gma = 357.52911+jc*(35999.05029 - 0.0001537*jc)
    gml = (280.46646+jc*(36000.76983 + jc*0.0003032))%360
    eeo = 0.016708634-jc*(0.000042037+0.0000001267*jc)
    moe = 23+(26+((21.448-jc*(46.815+jc*(0.00059-jc*0.001813))))/60.)/60
    oc = moe+0.00256*np.cos(radians(125.04-1934.136*jc))
    var_y = np.tan(radians(oc/2.))*np.tan(radians(oc/2.))
    sec = np.sin(radians(gma))*(1.914602-jc*(0.004817+0.000014*jc)) + np.sin(radians(2*gma))*(0.019993-0.000101*jc)+np.sin(radians(3*gma))*0.000289
    stl = gml + sec
    sal = stl-0.00569-0.00478*np.sin(radians(125.04-1934.136*jc))
    sun_dec = degrees(np.arcsin(np.sin(radians(oc))*np.sin(radians(sal))))

    et = 4.*degrees(var_y*np.sin(2.*radians(gml)) - 2.*eeo*np.sin(radians(gma)) + 4.*eeo*var_y*np.sin(radians(gma))*np.cos(2*radians(gml)) -\
            0.5*var_y*var_y*np.sin(4.*radians(gml)) - 1.25*eeo*eeo*np.sin(2.*radians(gma)))

    tst = (time*1440.+et+4*lon-60.*zone)%1440.
    if tst/4. < 0:
        hour_angle = tst/4.+180
    else: 
        hour_angle = tst/4.-180
    theta_i = degrees(np.arccos(np.sin(radians(lat))*np.sin(radians(sun_dec)) + np.cos(radians(lat))*np.cos(radians(sun_dec))*np.cos(radians(hour_angle))))

    if hour_angle>0:
        phi_i = (degrees(np.arccos(((np.sin(radians(lat))*np.cos(radians(theta_i)))-np.sin(radians(sun_dec)))/(np.cos(radians(lat))*np.sin(radians(theta_i)))))+180)%360.
    else:
        phi_i = (540-degrees(np.arccos(((np.sin(radians(lat))*np.cos(radians(theta_i)))-np.sin(radians(sun_dec)))/(np.cos(radians(lat))*np.sin(radians(theta_i))))))%360.
    return theta_i, phi_i
