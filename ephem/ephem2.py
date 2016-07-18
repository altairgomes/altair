import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, GCRS
from astropy.time import Time
from jplephem.spk import SPK
import astropy.constants as const
import os

#########################################################################################################

### Defines coordinates for each location
sites =
{
   '0': ['Greenwich',        EarthLocation.from_geodetic( +0.00000000000,   +51.47738889,      65.8)],
   '1': ['Paranal VLT',      EarthLocation.from_geodetic( +289.597194444,   -24.62541667,    2635.0)],
  '13': ['Liverpool Tel.',   EarthLocation.from_geodetic( +342.120800000,   +28.76240000,    2387.6)],
  '20': ['La Palma MERCATOR',EarthLocation.from_geodetic( +342.121694444,   +28.76200000,    2356.9)],
  '33': ['SOAR',             EarthLocation.from_geodetic( +289.269805556,   -30.23802778,    2693.9)],
  '86': ['Sierra Nevada',    EarthLocation.from_geodetic( +356.615305556,   +37.06413889,    2930.4)],
  '84': ['Cerro Tololo-DEC', EarthLocation.from_geodetic( +289.193583333,   -30.16958333,    2202.7)]
  '95': ['La Hita',          EarthLocation.from_geodetic( +356.813944444,   +39.56861111,     674.9)],
  '97': ['Wise Observatory', EarthLocation.from_geodetic( +34.7625000000,   +30.59566667,     904.2)],
 '267': ['CFHT',             EarthLocation.from_geodetic( +204.530444444,   +19.82536111,    4147.8)],
 '493': ['Calar Alto',       EarthLocation.from_geodetic( +357.453700000,   +37.22360000,    2168.0)],
 '511': ['Haute Provence',   EarthLocation.from_geodetic( +5.71569444400,   +43.93186111,     633.9)],
 '568': ['Mauna Kea',        EarthLocation.from_geodetic( +204.527805556,   +19.82611111,    4212.4)],
 '586': ['Pic du Midi',      EarthLocation.from_geodetic( +0.14230555600,   +42.93655556,    2890.5)],
 '809': ['ESO, La Silla',    EarthLocation.from_geodetic( +289.266250000,   -29.25883333,    2345.4)],
 '874': ['OPD/LNA',          EarthLocation.from_geodetic( +314.417361100,   -22.53550000,    1810.7)],
 '999': ['NARIT THAI',       EarthLocation.from_geodetic( +98.4666666667,   +18.56666667,    2457.0)],
}

### Defines type, name and kernel for each object
### Types: 's' for satellites, 'p' for planet
objects =
{
### Planets
    '599': ['p', 'Jupiter', 'jup300'], ## check kernel
### Satellites
    '506': ['s', 'Himalia', 'jup300'],
    '507': ['s', 'Elara', 'jup300'],
    '508': ['s', 'Pasiphae', 'jup300'],
    '509': ['s', 'Sinope', 'jup300'],
    '510': ['s', 'Lysithea', 'jup300'],
    '511': ['s', 'Carme', 'jup300'],
    '512': ['s', 'Ananke', 'jup300'],
    '513': ['s', 'Leda', 'jup300'],
    '607': ['s', 'Hyperion', 'jup300'], ## check kernel
    '608': ['s', 'Iapetus', 'jup300'], ## check kernel
    '609': ['s', 'Phoebe', 'sat359'], ## check kernel
    '802': ['s', 'Nereid', 'nep081']
### TNOs e Centauros
}

#########################################################################################################

def compute(k, center, target, jd):
    for segment in k.segments:
        print segment.center, segment.target
        if segment.center == center and segment.target == target and np.all(segment.start_jd <= jd) and np.all(jd <= segment.end_jd):
            return segment.compute(jd)
    raise ValueError('no segment matches')
    
def sitios(iaucode):
    try:
        return sites[iaucode][1], sites[iaucode][0]
    except ValueError:
        print("This code is not presented in our database. Please give the coordinates of the location")
        
def objs(code):
    objt, objn, objk = objects[code]
    return objt, objn, objk
        

##########################################################################################################

def ephem(iaucode, peph, objcode, timefile, output='sky'):

############## Reading planetary kernel ##############################################
    if peph[-3:] != 'bsp':
        peph = peph + '.bsp'
    if not os.path.isfile(peph):
        raise OSError('Planetary Kernel {} not found'.format(peph[:-4]))
    kernel = SPK.open(peph)

########## Reading Time File  #######################################################
    jd = np.loadtxt(timefile)
    time = Time(jd, format='jd', scale='utc')
    timetdb = time.tdb
  
######## Checking if it is geocenter or not and getting topocentric vectors ############
##### Still have to check when not geocenter #######
    if not isinstance(iaucode, str):
        iaucode = str(iaucode)
    if iaucode.lower() not in ['g', 'geocenter']:
        ### not geocenter
        site, siten = sitios(iaucode)
        #itrs = ITRS(site,obstime=time.utc)
        gcrs = site.itrs.transform_to(GCRS(obstime=time))
        print gcrs
    else:
        ### geocenter
        siten = 'Geocenter'

####################### fazer aqui ainda ################################################
    kernelobj = SPK.open('nep081.bsp')
    if not isinstance(objcode, str):
        objcode = str(objcode)
    objt, objn, objk = objs(objcode)
    objk = objk + '.bsp'
    if not os.path.isfile(kernelobj):
        raise OSError('Object Kernel {} not found'.format(objk[:-4]))
    kernelobj = SPK.open(kernelobj)
    
#########################################################################################

## compute vector Solar System Barycenter -> Earth Geocenter
    geop = kernel[0,3].compute(timetdb.jd) + kernel[3,399].compute(timetdb.jd) 

########### Calculates time delay ###########################################
# delt = light time
# n = number of loops necessary

    delt = 0*u.s
    n = 0
    while True:
        ### calculates new time
        tempo = timetdb - delt
        ### calculates vector Solar System Baricenter -> Object
        if objt in ['s', 'p']:
            objp = kernel[0,8].compute(tempo.jd) + kernelobj[8,802].compute(tempo.jd)[0:3]
        ### calculates vector Earth Geocenter -> Object
        position = (objp - geop)
        ### calculates linear distance Earth Geocenter -> Object
        dist = np.linalg.norm(position, axis=0)*u.km
        ### calculates new light time
        delt = (dist/const.c).decompose()
        n = n + 1
        ### if difference between new and previous light time is smaller than 0.001 sec, than continue.
        if np.all(np.absolute(((timetdb - tempo) - delt).sec) < 0.001):
            break
            
##############################################################################
### Creating SkyCoord Object with the coordinates for each instant
    if iaucode.lower() in ['g', 'geocenter']:
        ### geocenter
        coord = SkyCoord(position[0], position[1], position[2], frame='icrs', unit=u.km, representation='cartesian')
    else:
        ### topocenter
        coord = SkyCoord(position[0]*u.km - gcrs.cartesian.x, position[1]*u.km - gcrs.cartesian.y, position[2]*u.km - gcrs.cartesian.z, frame='icrs', representation='cartesian')
    
### returning values
    if output  == 'radec':
        return coord.spherical.lon.hourangle, coord.spherical.lat.deg
    elif output == 'sky':
        return SkyCoord(ra=coord.spherical.lon, dec=coord.spherical.lat)
    else:
        print siten


######################### TESTE ################################
a = ephem('g', 'de435', 802, 'jd_de.dat', output='sky')


