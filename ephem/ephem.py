import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, GCRS
from astropy.time import Time
from jplephem.spk import SPK
import astropy.constants as const

kernel = SPK.open('de430.bsp')
kernelh = SPK.open('jup300.bsp')
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

g = kernel[0,3].compute(time.jd) + kernel[3,399].compute(time.jd)

site = EarthLocation("314 25 2.5", "-22 32 7.8", 1864)
itrs = ITRS(site,obstime=time)
gcrs = itrs.transform_to(GCRS(obstime=time))

delt = 0*u.s
n = 0

while True:
    n = n + 1
    tempo = time - delt
    objp = kernel[0,5].compute(tempo.jd) + kernelh[5,506].compute(tempo.jd)[0:3] - g
    position = objp*u.km
#    velocity = (objv - geov)*(u.km/u.day)
#    print "velocity", velocity
    dist = np.linalg.norm(position)*u.km
    delt = (dist/const.c).decompose()
#    vemag = np.linalg.norm(velocity)*(u.km/u.day)
#    beta = (vemag/const.c).decompose()
#    print "vemag", vemag
#    print "beta", beta
#    dot = np.dot(position, velocity)*u.km*(u.km/u.day)
#    cosd = (dot/(dist*vemag)).decompose()
#    print "cosd", cosd
#    gammai = np.sqrt(1.0 - beta**2)
#    print "gammai", gammai
#    p = beta * cosd
#    print "p", p
#    q = (1.0 + p/(1.0+gammai))*delt
#    print "q", q
#    r = 1.0 + p
#    print "r", r
#    pos2 = ((gammai*position+q*velocity)/r).decompose()
#    print n, '{:1.15e}'.format(pos2[0])
    if np.all(np.absolute(time - tempo - delt) < 0.000001*u.s):
        break
        

#coord = SkyCoord(pos2[0], pos2[1], pos2[2], frame='icrs', unit=u.km, representation='cartesian')

coord = SkyCoord(position[0]*u.km - gcrs.cartesian.x*u.km, position[1]*u.km - gcrs.cartesian.y*u.km, position[2]*u.km - gcrs.cartesian.z*u.km, frame='icrs', representation='cartesian')

print coord.spherical.lon.hms, coord.spherical.lat.dms