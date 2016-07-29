import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, GCRS
from astropy.time import Time
from jplephem.spk import SPK
import astropy.constants as const
import os

#########################################################################################################

### Defines coordinates for each location
sites = {
   '0': ['Greenwich',        EarthLocation.from_geodetic( +0.00000000000,   +51.47738889,      65.8)],
   '1': ['Paranal VLT',      EarthLocation.from_geodetic( +289.597194444,   -24.62541667,    2635.0)],
  '13': ['Liverpool Tel.',   EarthLocation.from_geodetic( +342.120800000,   +28.76240000,    2387.6)],
  '20': ['La Palma MERCATOR',EarthLocation.from_geodetic( +342.121694444,   +28.76200000,    2356.9)],
  '33': ['SOAR',             EarthLocation.from_geodetic( +289.269805556,   -30.23802778,    2693.9)],
  '86': ['Sierra Nevada',    EarthLocation.from_geodetic( +356.615305556,   +37.06413889,    2930.4)],
  '84': ['Cerro Tololo-DEC', EarthLocation.from_geodetic( +289.193583333,   -30.16958333,    2202.7)],
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
### Types: 's' for satellites, 'p' for planet, 't' for TNOs and 'c' for Centaurs
### keys added must be in lower case
ephs = {
    'noe': 'NOE-5-2010-GAL.a',
    'lau': 'JIS-2015-10-sat',
}

objects = {}
### Planets
objects.update(dict.fromkeys([199, 'mercury'], ['p', 199, 'Mercury', 'de435']))
objects.update(dict.fromkeys([299, 'venus'], ['p', 299, 'Venus', 'de435']))
objects.update(dict.fromkeys([399, 'earth'], ['p', 399, 'Earth', 'de435']))
objects.update(dict.fromkeys([499, 'mars'], ['p', 499, 'Mars', 'mar097']))
objects.update(dict.fromkeys([599, 'jupiter'], ['p', 599, 'Jupiter', 'jup310']))
objects.update(dict.fromkeys([699, 'saturn'], ['p', 699, 'Saturn', 'sat359l']))
objects.update(dict.fromkeys([799, 'uranus'], ['p', 799, 'Uranus', 'ura112']))
objects.update(dict.fromkeys([899, 'neptune'], ['p', 899, 'Neptune', 'nep081']))
objects.update(dict.fromkeys([999, 'pluto'], ['p', 999, 'Pluto', 'plu043']))
### Satellites
    ## Sat of Jupiter
objects.update(dict.fromkeys([501, 'io'], ['s', 501, 'Io', ['jup310', ephs['noe']]]))
objects.update(dict.fromkeys([506, 'himalia'], ['s', 506, 'Himalia', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([507, 'elara'], ['s', 507, 'Elara', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([508, 'pasiphae'], ['s', 508, 'Pasiphae', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([509, 'sinope'], ['s', 509, 'Sinope', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([510, 'lysithea'], ['s', 510, 'Lysithea', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([511, 'carme'], ['s', 511, 'Carme', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([512, 'ananke'], ['s', 512, 'Ananke', ['jup300', ephs['lau']]]))
objects.update(dict.fromkeys([513, 'leda'], ['s', 512, 'Leda', ['jup300', ephs['lau']]]))
    ## Sat of Saturn
objects.update(dict.fromkeys([607, 'hyperion'], ['s', 607, 'Hyperion', 'sat375']))
objects.update(dict.fromkeys([608, 'iapetus'], ['s', 608, 'Iapetus', 'sat375']))
objects.update(dict.fromkeys([609, 'phoebe'], ['s', 609, 'Phoebe', ['sat375', 'ph15']]))
    ## Sat of Neptune
objects.update(dict.fromkeys([802, 'nereid'], ['s', 802, 'Nereid', 'nep081']))
    ## Sat of Pluto
objects.update(dict.fromkeys([901, 'charon'], ['s', 901, 'Charon', 'plu043']))
### TNOs e Centauros
objects.update(dict.fromkeys([2050000, 'quaoar', '2002lm60', 50000], ['t', 2050000, 'QUAOAR (2002 LM60)', '50000']))

#########################################################################################################

def compute(k, center, target, jd):
    def calculate(j):
        for segment in k.segments:
#            print segment.center, segment.target
            if segment.center == center and segment.target == target and segment.start_jd <= j and j <= segment.end_jd:
                return segment.compute(j)
        raise ValueError('no segment matches or segments does not cover at least one of the instants given')
        
    try:
        return k[center,target].compute(jd)
    except:
        pos = []
        if isinstance(jd,(int, float)):
            jd = [jd]
        for jj in jd:
            pos.append(calculate(jj)) 
        return np.array(pos).T
    
def sitios(iaucode):
    try:
        return sites[iaucode][1], sites[iaucode][0]
    except ValueError:
        print("This code is not presented in our database. Please give the coordinates of the location")
        
def objs(code, kern):
    if isinstance(code, str):
        code = code.lower()
    objt, objn, objm, objk = objects[code]
    print objt, objn, objm, objk
    if kern:
        objk = kern
    if isinstance(objk, list):
        print '{} kernels available for the object {}.'.format(len(objk), objm)
        for i in np.arange(len(objk)):
            print '{}: {}'.format(i, objk[i])
        n = input('Please choose the number refered to the kernel desired: ')
        if n > len(objk) - 1:
            raise IndexError('Kernel index out of range')
        objk = objk[n]
    itm = None
    if objk == 'JIS-2015-10-sat':
        itm = 599
    return objt, objn, objm, objk
        

##########################################################################################################

def ephem(iaucode, peph, objcode, timefile, kern=None, output='sky'):

############## Reading planetary kernel ##############################################
    if peph[-3:] != 'bsp':
        peph = peph + '.bsp'
    if not os.path.isfile(peph):
        raise OSError('Planetary Kernel {} not found'.format(peph[:-4]))
    kernel = SPK.open(peph)

########## Reading Time File  #######################################################
    if isinstance(timefile, str):
        jd = np.loadtxt(timefile)
        time = Time(jd, format='jd', scale='utc')
    elif isinstance(timefile, Time):
        time = timefile
    else:
        time = Time(timefile, format='jd', scale='utc')
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
    objt, objn, objm, objk = objs(objcode, kern)
    objk = objk + '.bsp'
    if not os.path.isfile(objk):
        raise OSError('Object Kernel {} not found'.format(objk[:-4]))
    kernelobj = SPK.open(objk)
    
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
        positm = np.array([0.0, 0.0, 0.0])
        if objt in ['s', 'p']:
            cc = objn//100
            objp = kernel[0,cc].compute(tempo.jd) + compute(kernelobj,cc,objn,tempo.jd)[0:3]
        if objt in ['t', 'c']:
            objp = kernel[0,10].compute(tempo.jd) + compute(kernelobj,10,objn,tempo.jd)[0:3]
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
    elif output == 'xyz':
        return coord.cartesian.x, coord.cartesian.y, coord.cartesian.z
    else:
        print siten


######################### TESTE ################################
#a = ephem('g', 'de430', 506, 'jd_de.dat', output='sky')
