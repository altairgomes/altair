import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, GCRS
from astropy.time import Time
from jplephem.spk import SPK
import astropy.constants as const

#########################################################################################################

def compute(k, center, target, jd, vel=False):
    for segment in k.segments:
        if segment.center == center and segment.target == target and segment.start_jd <= jd <= segment.end_jd:
            if vel == False:
                return segment.compute(jd)
            else:
                return np.array(segment.compute_and_differentiate(jd))
    raise ValueError('no segment matches')

##########################################################################################################

kernel = SPK.open('de435.bsp')
kernelnep = SPK.open('jup300.bsp')
time = Time(2458850.499918981, format='jd', scale='utc')
time = time.tdb

geop, geov = compute(kernel,0,3,time.jd,vel=True) + compute(kernel,3,399,time.jd,vel=True)

#site = EarthLocation("314 25 2.5", "-22 32 7.8", 1864)
#itrs = ITRS(site,obstime=time.utc)
#gcrs = itrs.transform_to(GCRS(obstime=time))

delt = 0*u.s
n = 0


while True:
    n = n + 1
    tempo = time - delt
    objp, objv = compute(kernel,0,5,tempo.jd,vel=True) + compute(kernelnep,5,506,tempo.jd,vel=True)[:,0:3]
    position = (objp - geop)*u.km
#    velocity = (objv - geov)*(u.km/u.day)
#    print "velocity", velocity
#    dist = np.linalg.norm(position)*u.km
#    delt = (dist/const.c).decompose()
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
    pos2 = position
#    print n, '{:1.15e}'.format(pos2[0])
    if np.all(np.absolute(time - tempo - delt) < 0.0001*u.s):
        break


coord = SkyCoord(pos2[0], pos2[1], pos2[2], frame='icrs', unit=u.km, representation='cartesian')

print coord.spherical.lon.hms, coord.spherical.lat.dms

#coord = SkyCoord(position[0]*u.km - gcrs.cartesian.x, position[1]*u.km - gcrs.cartesian.y, position[2]*u.km - gcrs.cartesian.z, frame='icrs', representation='cartesian')
