import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
import scipy.odr.odrpack as odrpack
from astropy.coordinates import SkyCoord

#####################################################################

f = open('entrada.dat', 'r')
arq = f.readlines()
f.close()

y = np.loadtxt(arq[0].strip(), skiprows=3, usecols=(2, 20, 21, 37, 16, 36), unpack=True) ## 2: JD; 20: dist_prim X; 21: dist_prim Y; 37: anom. verd.; 16: distancia; 36: anom. med.
z = np.loadtxt(arq[1].strip(), usecols=(0, 1, 2, 3), unpack=True) ## 0: off RA; 1: off DEC; 2: off_err RA; 3: off_err DEC

coords = np.loadtxt(arq[3].strip(), skiprows=3, usecols=(2, 4, 5, 6, 7, 8, 9), unpack=True, dtype ='S20', ndmin=1)
coor = coords[1]
for i in np.arange(len(coords))[2:]:
    coor = np.core.defchararray.add(coor, ' ')
    coor = np.core.defchararray.add(coor, coords[i])
lau = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))

tempo = np.loadtxt(arq[4].strip(), skiprows=3, usecols=(2,), unpack=True, ndmin=1)
ephcoord = np.loadtxt(arq[4].strip(), skiprows=3, usecols=(2, 4, 5, 6, 7, 8, 9), unpack=True, dtype ='S20', ndmin=1)
coor = coords[1]
for i in np.arange(len(ephcoord))[2:]:
    coor = np.core.defchararray.add(coor, ' ')
    coor = np.core.defchararray.add(coor, ephcoord[i])
jpl = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))

dalfa = lau.ra*np.cos(lau.dec) - jpl.ra*np.cos(jpl.dec)
ddelta = lau.dec - jpl.dec

r = []
for i in np.arange(1995,2018,1):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

#print r
r = np.array(r)

############## Declinacao ############################################


plt.plot(tempo - 2451544.5, ddelta.mas, label='Declination')
plt.errorbar(y[0] - 2451544.5, z[1], yerr=z[3], fmt='s', label='Offsets')
#plt.xlim(2449353.5 - 2451544.,2457754.5 - 2457388.5)
#plt.vlines(eph[0][j] - 2451544.5, -500, 500)
plt.ylim(-500,500)
plt.title('Laurene - JUP300')
plt.xlabel('Time')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1995,2018,1)])
plt.ylabel('Offset in DEC (mas)')
plt.legend()
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('LAU-JPL_DEC.png',dpi=100, bbox_inches='tight')
plt.clf()



plt.plot(tempo - 2451544.5, dalfa.mas, label='Right Ascencion')
plt.errorbar(y[0] - 2451544.5, z[0], yerr=z[2], fmt='s', label='Offsets')
plt.title('Laurene - JUP300')
plt.ylim(-500,500)
plt.xlabel('Time')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1995,2018,1)])
plt.ylabel('Offset in RA (mas)')
plt.legend()
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('LAU-JPL_RA.png',dpi=100, bbox_inches='tight')
plt.clf()
