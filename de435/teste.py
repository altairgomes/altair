import numpy as np
from observation import read
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord
import os

sat = ['Ananke', 'Carme', 'Elara', 'Himalia', 'Leda', 'Lysithea', 'Pasiphae', 'Sinope']
site = [874, 511, 809]

for i in sat:
    print i
    coord, tel, time = read(i+'.dat', coord_col=[0,1,2,3,4,5], time_col=[8], comment_col=[11])
    sit = np.loadtxt(i+'.dat', usecols=[12])
    for s in site:
        print s
        ecoord, anom, ejd = read('{}_{}.eph'.format(i,s), skiprows=3, coord_col=[4,5,6,7,8,9], time_col=[2], comment_col=[37])
        anov = np.array(anom).astype(np.float)
        a,b = np.indices((len(time), len(ejd)))
        c,d = np.where(np.absolute(time.jd[a] - ejd.jd[b]) < TimeDelta(0.2*u.s).jd)
        f = open('{}_{}_off.dat'.format(i,s), 'w')
        for k in np.arange(len(c)):
            f.write(' {:-6.3f} {:-6.3f} {:16.8f} {:6.2f} {:2s}\n'.format((coord.ra[c][k]*np.cos(coord.dec[c][k])-ecoord.ra[d][k]*np.cos(ecoord.dec[d][k])).arcsec,
            (coord.dec[c][k]-ecoord.dec[d][k]).arcsec,
            time.jd[c][k],
            anov[d][k],
            tel[c][k]))
        f.close()
