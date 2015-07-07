import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
import scipy.odr.odrpack as odrpack

######################################################################

f = open('entrada.dat', 'r')
arq = f.readlines()
f.close()

y = np.loadtxt(arq[0].strip(), skiprows=3, usecols=(2, 20, 21, 37, 16, 36), unpack=True) ## 2: JD; 20: dist_prim X; 21: dist_prim Y; 37: anom. verd.; 16: distancia; 36: anom. med.
z = np.loadtxt(arq[1].strip(), usecols=(0, 1, 2, 3), unpack=True) ## 0: off RA; 1: off DEC; 2: off_err RA; 3: off_err DEC

eph = np.loadtxt(arq[2].strip(), skiprows=3, usecols=(2, 35, 16, 34), unpack=True) ## 2: JD; 37: anom. verd.; 16: distancia; 36: anom. med.

k = np.arange(361)

m = np.append(eph[1][0], eph[1])[0:-1]
n = eph[1]-m
j = np.where(n < -100.0)

r = []
for i in np.arange(1995,2016,1):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

print r
r = np.array(r)

f = open('saida.dat', 'w')

############## Funcoes ##############################################

g1 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

g2 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

g3 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg)**2, np.cos(y[3]*u.deg)**2, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

g4 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg), (y[0] - 2451544.5)*np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), (y[0] - 2451544.5)*np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

def f1(B, x): ## tempo, seno, cosseno, constante
    return B[0]*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg) + B[2]*np.cos(x[1]*u.deg) + B[3]
    
def f2(B, x): ## tempo, seno*cos, sen, cos, constante
    return B[0]*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[2]*np.sin(x[1]*u.deg) + B[3]*np.cos(x[1]*u.deg) + B[4]
    
def f3(B, x): ## tempo, sen^2, cos^2, seno*cos, sen, cos, constante
    return B[0]*(x[0] - 2451544.5) + B[1]*(np.sin(x[1]*u.deg)**2) + B[2]*(np.cos(x[1]*u.deg)**2) + B[3]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[4]*np.sin(x[1]*u.deg) + B[5]*np.cos(x[1]*u.deg) + B[6]
    
def f4(B, x): ## tempo, sen^2(f), cos^2(f), seno(f)*cos(f), sen(f), cos(f), constante
    return B[0]*(x[0] - 2451544.5) + (B[1]+B[2]*(x[0] - 2451544.5))*np.sin(x[1]*u.deg) + (B[3]+B[4]*(x[0] - 2451544.5))*np.cos(x[1]*u.deg) + B[5]

def f5(B, x): ## a*sen(wt+p), sin(f), cos(f), constante
    return B[0]*np.sin(2*np.pi/(B[1]*365.25)*(x[0] - 2451544.5) + B[2]) + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg) + B[5]
    
def f6(B, x): ## a*sen(wt+p), constante
    return B[0]*np.sin(2*np.pi/(B[1]*365.25)*(x[0] - 2451544.5) + B[2]) + B[3]


############## Declinacao ############################################

p = np.linalg.lstsq(g1, z[1])

x = np.vstack([y[0], y[3]])

linearde = odrpack.Model(f5)
#mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
mydatade = odrpack.RealData(x, z[1], sy=z[3])
myodrde = odrpack.ODR(mydatade, linearde, beta0=[200.0, 11.7, 0.0, 1.0, 1.0, 0.0])
myodrde.set_job(fit_type=2)
myoutputde = myodrde.run()
print 'Declinacao\n'
myoutputde.pprint()

f.write('\n\nDeclinacao\n')
f.write('Resultados: {}\n'.format(myodrde.output.beta))
f.write('Erros: {}\n'.format(myodrde.output.sd_beta))

#print p[0]
#print "Desvio padrao DEC:", np.sqrt(((z[1]/z[3] - (A*p[0]).sum(axis=1))**2).sum())

plt.errorbar(y[0] - 2451544.5, z[1], yerr=z[3], fmt='s', label='Offsets')
#plt.plot(y[3], z[1], 's', label='Offsets')
plt.plot(eph[0] - 2451544.5, f5(myodrde.output.beta, np.vstack([eph[0], eph[1]])), label='Ajuste1')
#plt.plot(eph[0] - 2451544.5, f(p[0], np.vstack([eph[0], eph[1]])), label='Ajuste2')
#plt.title('DEC = {:.2f}*sen(A.V.) + {:.2f}*cos(A.V.) + {:.2f}'.format(p[0][0], p[0][1], p[0][2]))
#plt.xlim(0,360)
plt.vlines(eph[0][j] - 2451544.5, -500, 500)
plt.ylim(-500,500)
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1995,2016,1)])
plt.ylabel('Offset (mas)')
plt.legend()
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches(20.0,8.0)
fig.savefig('DEC.png',dpi=100, bbox_inches='tight')
#plt.savefig('DEC.png')
plt.clf()

plt.errorbar(y[3], z[1], yerr=z[3], fmt='s', label='Offsets')
plt.ylim(-500,500)
plt.xlim(0,360)
plt.xlabel('Anomalia Verdadeira')
plt.ylabel('Offset (mas)')
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches(20.0,8.0)
fig.savefig('DEC_anom.png',dpi=100, bbox_inches='tight')
plt.clf()

#################### Ascensao Reta ##################################

q = np.linalg.lstsq(g2, z[0])

x = np.vstack([y[0], y[3]])

linearra = odrpack.Model(f5)
#mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
mydatara = odrpack.RealData(x, z[0], sy=z[2])
myodrra = odrpack.ODR(mydatara, linearra, beta0=[200.0, 11.7, 0.0, 1.0, 1.0, 0.0])
myodrra.set_job(fit_type=2)
myoutputra = myodrra.run()
print '\n\nAscensao Reta\n'
myoutputra.pprint()


f.write('\n\nAscensao Reta\n')
f.write('Resultados: {}\n'.format(myodrra.output.beta))
f.write('Erros: {}\n'.format(myodrra.output.sd_beta))

#print q[0]
#print "Desvio padrao RA:", np.sqrt(((z[0]/z[2] - (B*q[0]).sum(axis=1))**2).sum())

plt.errorbar(y[0] - 2451544.5, z[0], yerr=z[2], fmt='s', label='Offsets')
#plt.plot(y[3], z[0], 's', label='Offsets')
plt.plot(eph[0] - 2451544.5, f5(myodrra.output.beta, np.vstack([eph[0], eph[1]])), label='Ajuste1')
#plt.plot(eph[0] - 2451544.5, g(q[0], np.vstack([eph[0], eph[1]])), label='Ajuste2')
#plt.title('RA = {:.2f}*sen^2(A.V.) + {:.2f}*cos^2(A.V.) + {:.2f}*sen(A.V.)*cos(A.V.) + {:.2f}*sen(A.V.) + {:.2f}*cos(A.V.) + {:.2f}'.format(q[0][0], q[0][1], q[0][2], q[0][3], q[0][4], q[0][5]))
#plt.xlim(0,360)
plt.vlines(eph[0][j] - 2451544.5, -500, 500)
plt.ylim(-500,500)
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1995,2016,1)])
plt.ylabel('Offset (mas)')
plt.legend()
plt.axhline(0, color='black')
fig = plt.gcf()
fig.set_size_inches(20.0,8.0)
fig.savefig('RA.png',dpi=100, bbox_inches='tight')
plt.clf()

plt.errorbar(y[3], z[0], yerr=z[2], fmt='s', label='Offsets')
plt.ylim(-500,500)
plt.xlim(0,360)
plt.xlabel('Anomalia Verdadeira')
plt.ylabel('Offset (mas)')
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches(20.0,8.0)
fig.savefig('RA_anom.png',dpi=100, bbox_inches='tight')

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
