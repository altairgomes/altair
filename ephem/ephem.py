import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from jplephem.spk import SPK

kernel = SPK.open('de430.bsp')
kernelnep = SPK.open('nep081.bsp')
time = Time(2457023.5, format='jd', scale='utc')
site = EarthLocation("314 25 2.5", "-22 32 7.8", 1864)

position = kernel[0,8].compute(time.jd) + kernelnep[8,899].compute(time.jd)[0:3] - kernel[0,3].compute(time.jd) - kernel[3,399].compute(time.jd)

itrs = ITRS(site.itrs,obstime=time)
gcrs = itrs.transform_to(GCRS(obstime=time))

coord = SkyCoord(position[0], position[1], position[2], frame='icrs', unit=u.km, representation='cartesian')

coord = SkyCoord(position[0]*u.km - gcrs.cartesian.x, position[1]*u.km - gcrs.cartesian.y, position[2]*u.km - gcrs.cartesian.z, frame='icrs', representation='cartesian')
