import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
import scipy.odr.odrpack as odrpack

######################################################################
y = np.loadtxt('Carme_ephem.dat', skiprows=3, usecols=(2, 20, 21, 37, 16, 36), unpack=True) ## 2: JD; 20: dist_prim X; 21: dist_prim Y; 37: anom. verd.; 16: distancia; 36: anom. med.
z = np.loadtxt('Carme_total', usecols=(0, 1, 2, 3), unpack=True) ## 0: off RA; 1: off DEC; 2: off_err RA; 3: off_err DEC

k = np.arange(361)

############## Declinacao ############################################

A = np.vstack([y[0] - 2451544.5, np.sin(y[5]*u.deg)*np.cos(y[5]*u.deg), np.sin(y[5]*u.deg), np.cos(y[5]*u.deg), np.ones(len(y[3]))]).T

#p = np.linalg.lstsq(A, z[1])

def f(B, x):
    return B[0]*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg) + B[2]*np.cos(x[1]*u.deg) + B[3]

x = np.vstack([y[0], y[5]])

linearde = odrpack.Model(f)
#mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
mydata = odrpack.RealData(x, z[1], sy=z[3])
print "Declinacao"
myodr = odrpack.ODR(mydata, linearde, beta0=[1.78854119e-03, -1.21822344e+00, 9.88555083e+00, -6.08946292e+00])
myoutput = myodr.run()
myoutput.pprint()

#print p[0]
#print "Desvio padrao DEC:", np.sqrt(((z[1]/z[3] - (A*p[0]).sum(axis=1))**2).sum())

plt.errorbar(y[3], z[1], yerr=z[3], fmt='s', label='Offsets')
#plt.plot(y[3], z[1], 's', label='Offsets')
for i in np.arange(2000,2010,1):
    plt.plot(f(myodr.output.beta, np.vstack([Time('{}-01-01 00:00:00'.format(i), format='iso').jd*np.ones(len(k)),k])), label=i)
#plt.title('DEC = {:.2f}*sen(A.V.) + {:.2f}*cos(A.V.) + {:.2f}'.format(p[0][0], p[0][1], p[0][2]))
plt.xlim(0,360)
plt.ylim(-500,500)
plt.xlabel('Anomalia Verdadeira')
plt.ylabel('Offset (mas)')
plt.legend()
plt.axhline(0, color='red')
fig =plt.gcf()
fig.set_size_inches(15.0,7.5)
fig.savefig('DEC.png',dpi=100, bbox_inches='tight')
#plt.savefig('DEC.png')
plt.clf()


#################### Ascensao Reta ##################################

B = np.vstack([y[0] - 2451544.5, np.sin(y[5]*u.deg)*np.cos(y[5]*u.deg), np.sin(y[5]*u.deg), np.cos(y[5]*u.deg), np.ones(len(y[3]))]).T

#q = np.linalg.lstsq(B, z[0])

def g(B, x):
    return B[0]*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[2]*np.sin(x[1]*u.deg) + B[3]*np.cos(x[1]*u.deg) + B[4]

x = np.vstack([y[0], y[5]])

linearra = odrpack.Model(g)
#mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
mydata = odrpack.RealData(x, z[0], sy=z[2])
print "Ascensao Reta"
myodr = odrpack.ODR(mydata, linearra, beta0=[3.62659327e-04,7.82345630e+00, -6.40102501e-01, -2.83655752e-01, -9.32980254e-02])
myoutput = myodr.run()
myoutput.pprint()

#print q[0]
#print "Desvio padrao RA:", np.sqrt(((z[0]/z[2] - (B*q[0]).sum(axis=1))**2).sum())

plt.errorbar(y[3], z[0], yerr=z[2], fmt='s', label='Offsets')
#plt.plot(y[3], z[0], 's', label='Offsets')
for i in np.arange(2000,2010,1):
    plt.plot(g(myodr.output.beta, np.vstack([Time('{}-01-01 00:00:00'.format(i), format='iso').jd*np.ones(len(k)),k])), label=i)
#plt.title('RA = {:.2f}*sen^2(A.V.) + {:.2f}*cos^2(A.V.) + {:.2f}*sen(A.V.)*cos(A.V.) + {:.2f}*sen(A.V.) + {:.2f}*cos(A.V.) + {:.2f}'.format(q[0][0], q[0][1], q[0][2], q[0][3], q[0][4], q[0][5]))
plt.xlim(0,360)
plt.ylim(-500,500)
plt.xlabel('Anomalia Verdadeira')
plt.ylabel('Offset (mas)')
plt.legend()
plt.axhline(0, color='red')
fig = plt.gcf()
fig.set_size_inches(15.0,7.5)
fig.savefig('RA.png',dpi=100, bbox_inches='tight')


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
