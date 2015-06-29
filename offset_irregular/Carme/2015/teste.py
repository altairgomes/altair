import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time

######################################################################
y = np.loadtxt('Carme_ephem.dat', skiprows=3, usecols=(2, 20, 21, 37, 16, 36), unpack=True) ## 2: JD; 20: dist_prim X; 21: dist_prim Y; 37: anom. verd.; 16: distancia; anom. med.
z = np.loadtxt('Carme_total', usecols=(0, 1, 2, 3), unpack=True) ## 0: off RA; 1: off DEC; 2: off_err RA; 3: off_err DEC

k = np.arange(361)

############## Declinacao ############################################

A = np.vstack([ np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

p = np.linalg.lstsq(A, z[1])

plt.errorbar(y[3], z[1], yerr=z[3], fmt='s', label='Offsets')
#plt.plot(y[3], z[1], 's', label='Offsets')
plt.plot(p[0][0]*np.sin(np.pi*k/180) + p[0][1]*np.cos(np.pi*k/180) + p[0][2], 'g', label='Fitted line')
plt.title('DEC = {:.2f}*sen(A.V.) + {:.2f}*cos(A.V.) + {:.2f}'.format(p[0][0], p[0][1], p[0][2]))
plt.xlim(0,360)
plt.ylim(-500,500)
plt.xlabel('Anomalia Verdadeira')
plt.ylabel('Offset (mas)')
plt.legend()
plt.axhline(0, color='red')
fig =plt.gcf()
fig.set_size_inches(15.0,7.5)
fig.savefig('DEC.png',dpi=100)
#plt.savefig('DEC.png')
plt.clf()
print p[1]


#################### Ascensao Reta ##################################

B = np.vstack([np.sin(y[3]*u.deg)**2, np.cos(y[3]*u.deg)**2, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

q = np.linalg.lstsq(B, z[0])

plt.errorbar(y[3], z[0], yerr=z[2], fmt='s', label='Offsets')
#plt.plot(y[3], z[0], 's', label='Offsets')
plt.plot(q[0][0] * np.sin(np.pi*k/180)**2 + q[0][1] * np.cos(np.pi*k/180)**2 + q[0][2]*np.sin(np.pi*k/180)*np.cos(np.pi*k/180) + q[0][3]*np.sin(np.pi*k/180)**2 + q[0][4]*np.cos(np.pi*k/180)**2 + q[0][5], 'g', label='Fitted line')
plt.title('RA = {:.2f}*sen^2(A.V.) + {:.2f}*cos^2(A.V.) + {:.2f}*sen(A.V.)*cos(A.V.) + {:.2f}*sen(A.V.) + {:.2f}*cos(A.V.) + {:.2f}'.format(q[0][0], q[0][1], q[0][2], q[0][3], q[0][4], q[0][5]))
plt.xlim(0,360)
plt.ylim(-500,500)
plt.xlabel('Anomalia Verdadeira')
plt.ylabel('Offset (mas)')
plt.legend()
plt.axhline(0, color='red')
fig = plt.gcf()
fig.set_size_inches(15.0,7.5)
fig.savefig('RA.png',dpi=100)
print q[1]

#plt.savefig('RA.png', figsize(10,8))

#off_ephem = np.loadtxt('final.eph', usecols=(2, 35), unpack=True)
#dd = Time(off_ephem[0], format='jd', scale='utc')
#off_saida = open('offsets.txt', 'w')
#for idx, val in enumerate(off_ephem[1]):
#    ra_off = q[0][0] * np.sin(np.pi*val/180)**2 + q[0][1] * np.cos(np.pi*val/180)**2 + q[0][2]*np.sin(np.pi*val/180)*np.cos(np.pi*val/180) + q[0][3]*np.sin(np.pi*val/180)**2 + q[0][4]*np.cos(np.pi*val/180)**2 + q[0][5]
#    dec_off = p[0][0]*np.sin(np.pi*val/180) + p[0][1]*np.cos(np.pi*val/180) + p[0][2]
#    off_saida.write('{} {:5.1f} {:4.0f} {:4.0f}\n'.format(dd[idx].isot, val, ra_off, dec_off))


#plt.plot(y[3],z[1], 's', label='offsets')
#plt.xlim([0,360])
#plt.ylim([-300,300])
#plt.axhline(0, color='red')
#plt.legend()
#plt.savefig('teste.png')
