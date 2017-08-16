from scipy.fftpack import fft, fftfreq, fftshift
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

k = np.arange(361)

m = np.append(eph[1][0], eph[1])[0:-1]
n = eph[1]-m
j = np.where(n < -100.0)

r = []
for i in np.arange(1990,2021,1):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

r = np.array(r)

#########################################################################
x = eph[0]
y = eph[8]

# number of signal points
N = len(eph[0])
 # sample spacing
T = (x.max()-x.min()) / N
yf = fft(y)
xf = fftfreq(N, T)
xf = fftshift(xf)
yplot = fftshift(yf)
fig, ax = plt.subplots(2, 1)
ax[0].plot(x,y)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
#ax[0].set_ticks(r)
#ax[0].set_ticklabels(['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
ax[1].plot(xf, 1.0/N * np.abs(yplot),'r') # plotting the spectrum
ax[1].set_xlabel(r'Freq ($day^{-1}$)')
ax[1].set_ylabel('|Y(freq)|')
#plt.plot(xf, 1.0/N * np.abs(yplot))
plt.grid()
ax[1].set_xlim(0.00,0.01)
plt.savefig('freq_y.png', dpi=100)

