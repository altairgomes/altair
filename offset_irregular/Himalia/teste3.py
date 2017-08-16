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

f = open('saida.dat', 'w')

############## Funcoes ##############################################

g1 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

g2 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

g3 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg)**2, np.cos(y[3]*u.deg)**2, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

g4 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg), (y[0] - 2451544.5)*np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), (y[0] - 2451544.5)*np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

def f1(B, x): ## tempo, seno, cosseno, constante
    return (B[0]/365)*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg) + B[2]*np.cos(x[1]*u.deg) + B[3]
    
def f2(B, x): ## tempo, seno*cos, sen, cos, constante
    return B[0]*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[2]*np.sin(x[1]*u.deg) + B[3]*np.cos(x[1]*u.deg) + B[4]
    
def f3(B, x): ## tempo, sen^2, cos^2, seno*cos, sen, cos, constante
    return B[0]*(x[0] - 2451544.5) + B[1]*(np.sin(x[1]*u.deg)**2) + B[2]*(np.cos(x[1]*u.deg)**2) + B[3]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[4]*np.sin(x[1]*u.deg) + B[5]*np.cos(x[1]*u.deg) + B[6]
    
def f4(B, x): ## tempo, sen^2(f), cos^2(f), seno(f)*cos(f), sen(f), cos(f), constante
    return B[0]*(x[0] - 2451544.5) + (B[1]+B[2]*(x[0] - 2451544.5))*np.sin(x[1]*u.deg) + (B[3]+B[4]*(x[0] - 2451544.5))*np.cos(x[1]*u.deg) + B[5]

def f5(B, x): ## a*cos(wt+p), sin(f), cos(f), constante
    return B[0]*np.cos((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180)) + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg) + B[5]
    
def f6(B, x): ## a*sen(wt+p), sin(f), cos(f), constante
    return B[0]*np.sin((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180)) + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg) + B[5]
    
def f7(B, x):
    return B[0]*np.sin((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180)) + B[3]*np.sin((2*np.pi/(B[4]*365.25))*(x[0] - 2451544.5) + B[5]*(np.pi/180)) + B[6]*np.sin((2*np.pi/(B[7]*365.25))*(x[0] - 2451544.5) + B[8]*(np.pi/180)) + B[9]

def f8(B, x):
    return np.sin((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180))*(B[0] + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg)) + B[5]
    
def f9(B,x):
    return B[0] + B[2]*np.cos(np.pi*x[0]/(B[1]*365.25)) + B[3]*np.sin(np.pi*x[0]/(B[1]*365.25)) + B[4]*np.cos(2*np.pi*x[0]/(B[1]*365.25)) + B[5]*np.sin(2*np.pi*x[0]/(B[1]*365.25)) + B[6]*np.cos(3*np.pi*x[0]/(B[1]*365.25)) + B[7]*np.sin(3*np.pi*x[0]/(B[1]*365.25))

############## least square function #################################

def least(func, x, y, sy=None, beta0=[]):
    linear = odrpack.Model(func)
    #mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
    if not sy == None:
        mydata = odrpack.RealData(x, y, sy=sy)
    else:
        mydata = odrpack.RealData(x, y)
    myodr = odrpack.ODR(mydata, linear, beta0=beta0)
    myodr.set_job(fit_type=2)
    myoutput = myodr.run()
    myoutput.pprint()
    return myodr.output
    
def residuos(func, par, x, y, sy=[]):
    if not sy == []:
        var = np.sum(((1/sy)**2)*((y - func(par, x))**2))/np.sum((1/sy)**2)
    else:
        var = np.sum(((y - func(par, x))**2)/len(y))
    resid = np.sqrt(var)
    return resid

############## Y ############################################

p = np.linalg.lstsq(g1, z[1])

x = np.vstack([eph[0], eph[1]])

print 'Y\n'
f.write('\nY\n')

fun=f9
beta0=[1.0, 12.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]

ajdenwg = least(func=fun, x=x, y=eph[8], beta0=beta0)
f.write('\nResultados do ajuste sem peso:\n')
for i in np.arange(len(ajdenwg.beta)):
    f.write('p[{}]: {} \pm {}\n'.format(i, ajdenwg.beta[i], ajdenwg.sd_beta[i]))
resid = residuos(fun, par=ajdenwg.beta, x=x, y=eph[8])
f.write('Residuos: {} mas\n'.format(resid))

plt.plot(eph[0] - 2451544.5, eph[8], label='Offsets')
plt.plot(eph[0] - 2451544.5, fun(ajdenwg.beta, x), label='Ajuste')
for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Offset (arcsec)')
plt.legend()
plt.axhline(0, color='black')
fig =plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('Himalia-x.png',dpi=100, bbox_inches='tight')
#plt.savefig('DEC.png')
plt.clf()

#################### X ##################################

q = np.linalg.lstsq(g2, z[0])

x = np.vstack([eph[0], eph[1]])

print 'X\n'
f.write('\n\nX\n')

fun=f9
#beta0=[0.0, 12.0, 1000.0, 100.0]

ajranwg = least(func=fun, x=x, y=eph[7], beta0=beta0)
f.write('\nResultados do ajuste sem peso:\n')
for i in np.arange(len(ajranwg.beta)):
    f.write('p[{}]: {} \pm {}\n'.format(i, ajranwg.beta[i], ajranwg.sd_beta[i]))
resid = residuos(fun, par=ajranwg.beta, x=x, y=eph[7])
f.write('Residuos: {} mas\n'.format(resid))

plt.plot(eph[0] - 2451544.5, eph[7], label='Offsets')
plt.plot(eph[0] - 2451544.5, fun(ajranwg.beta, x), label='Ajuste')
for i in eph[0][j]:
    plt.axvline(i - 2451544.5, color='black')
plt.xlabel('Tempo')
plt.xticks(r, ['{}'.format(i) for i in np.arange(1990,2021,1)], rotation=90)
plt.ylabel('Offset (arcsec)')
plt.legend()
plt.axhline(0, color='black')
fig = plt.gcf()
fig.set_size_inches((40.0*u.cm).to(u.imperial.inch).value,(16.0*u.cm).to(u.imperial.inch).value)
fig.savefig('RA.png',dpi=100, bbox_inches='tight')

