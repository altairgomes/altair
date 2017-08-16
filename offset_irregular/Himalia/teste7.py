import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from scipy.fftpack import fft, fftfreq, fftshift

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
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso', scale='utc').jd - 2451544.5)

r = np.array(r)

#########################################################################
x = eph[0]
yra = eph[7]

# number of signal points
N = len(eph[0])
 # sample spacing
T = (x.max()-x.min()) / N
yra_fft = np.fft.fft(yra)
e = np.zeros(len(yra_fft))
xf1 = fftfreq(N, T)
xf = fftshift(xf1)
yplot = fftshift(yra_fft)
r = np.argsort(np.abs(yra_fft[0:N/2])*2.0/N)
k = 30
e[-r[-k:]] = 1.0
e[r[-k:]] = 1.0
yra_fit = np.fft.ifft(yra_fft*e).real
print "ra: ", np.mean((np.abs(yra_fit-yra)[N/4:-N/4]*u.arcsec).to(u.mas))
fig, ax = plt.subplots(3, 1)
ax[0].plot(x,yra, color='green')
ax[0].plot(x,yra_fit, color='red')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[1].plot(xf, 1.0/N * np.abs(yplot),'r') # plotting the spectrum
ax[1].plot(xf1[r[-k:]], 1.0/N * np.abs(yra_fft)[r[-k:]], 's')
ax[1].set_xlabel(r'Freq ($day^{-1}$)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].set_xlim(0.00,0.01)
ax[2].plot(x,yra_fit-yra)
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Delta Amplitude')
plt.savefig('freq_x.png', dpi=100)

yde = eph[8]

# number of signal points
#N = len(eph[0])
 # sample spacing
#T = (x.max()-x.min()) / N
#yde_fft = np.fft.fft(yde)
#x=1
#while True:
#    e = np.zeros(len(y_fft))
#    r = np.argsort(np.abs(y_fft[0:N/2])*2.0/N)
#    e[-r[-x:]] = 1.0
#    e[r[-x:]] = 1.0
#    y_fit = np.fft.ifft(y_fft*e).real
#    if np.mean((np.abs(y_fit-y)*u.arcsec).to(u.mas)) < 1000*u.mas:
#        break
#    x = x+1

#print x
