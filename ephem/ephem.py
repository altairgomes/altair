import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, GCRS
from astropy.time import Time
from jplephem.spk import SPK
import astropy.constants as const

kernel = SPK.open('de430.bsp')
kernelnep = SPK.open('nep081.bsp')
time = Time(2457023.5, format='jd', scale='utc')
delta_at = 35*u.s
delta_t_a = 32.184*u.s
k = float(1.657e-3)*u.s
eb = float(1.671e-2)
n = np.array([6.239996e0, 1.99096871e-7])
t = time - Time('j2000')
m = n[0] + n[1]*t.sec
e = m + eb*np.sin(m)
delta_et =  delta_t_a  + k*np.sin(e) + delta_at

site = EarthLocation("314 25 2.5", "-22 32 7.8", 1864)
itrs = ITRS(site,obstime=time)
gcrs = itrs.transform_to(GCRS(obstime=time))

delt = 0*u.s
n = 0

topo = kernel[0,3].compute(time.jd) + kernel[3,399].compute(time.jd)

while True:
    n = n + 1
    tempo = time - delt
    position = kernel[0,8].compute(tempo.jd) + kernelnep[8,899].compute(tempo.jd)[0:3] - topo
    dist = np.linalg.norm(position)*u.km
    delt = (dist/const.c).decompose()
    if np.all(np.absolute(time - tempo - delt) < 1*u.s):
        break

coord = SkyCoord(position[0], position[1], position[2], frame='icrs', unit=u.km, representation='cartesian')

#coord = SkyCoord(position[0]*u.km - gcrs.cartesian.x, position[1]*u.km - gcrs.cartesian.y, position[2]*u.km - gcrs.cartesian.z, frame='icrs', representation='cartesian')
