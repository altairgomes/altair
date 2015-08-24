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

eph = np.loadtxt(arq[2].strip(), skiprows=3, usecols=(2, 35, 16, 34, 29, 30, 31, 20, 21), unpack=True) ## 2: JD; 35: anom. verd.; 16: distancia; 34: anom. med.; 29: semi-eixo; 30: excentricidade; 31: inclinacao

nome = arq[2][0:-5]

k = np.arange(361)

m = np.append(eph[1][0], eph[1])[0:-1]
n = eph[1]-m
j = np.where(n < -100.0)

r = []
for i in np.arange(1990,2021,1):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

r = np.array(r)

f = open('saida.dat', 'w')

############## Semi-Eixo ############################################

plt.plot(eph[0] - 2451544.5, eph[4])

for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Semi-Eixo (km)')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('{}-a.png'.format(nome),dpi=200, bbox_inches='tight')
plt.clf()

#################### Excentricidade ##################################

plt.plot(eph[0] - 2451544.5, eph[5])

for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Excentricidade')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('{}-e.png'.format(nome),dpi=200, bbox_inches='tight')
plt.clf()

#################### Inclinacao ##################################

plt.plot(eph[0] - 2451544.5, eph[6])

for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Inclinacao (graus)')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('{}-I.png'.format(nome),dpi=200, bbox_inches='tight')
plt.clf()

#################### Dist - x ##################################

plt.plot(eph[0] - 2451544.5, eph[7])
plt.plot(y[0] - 2451544.5, z[0] + y[1], 's', color='red')

for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Inclinacao (graus)')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('{}-x.png'.format(nome),dpi=200, bbox_inches='tight')
plt.clf()

#################### Dist - y ##################################

plt.plot(eph[0] - 2451544.5, eph[8])
plt.plot(y[0] - 2451544.5, z[1] + y[2], 's',  color='red')

for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Inclinacao (graus)')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('{}-y.png'.format(nome),dpi=200, bbox_inches='tight')
plt.clf()
