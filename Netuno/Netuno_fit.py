import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import scipy.odr.odrpack as odrpack

##################################### FUNCOES #######################

def time_hourangle(time, ra):
    time.delta_ut1_utc = 0.0
    time.location = lna
    hourangle = time.sidereal_time('mean') - ra
    hourangle.wrap_at('180d', inplace=True)
    return hourangle
    
def refraction(delta, hourangle):
    den = np.tan(delta)*np.tan(lna.latitude)+np.cos(hourangle)
    valfa = np.sin(hourangle)/(den*np.cos(delta))
    vdelta = (np.tan(lna.latitude)-np.tan(delta)*np.cos(hourangle))/den
    return valfa, vdelta
    


#####################################################################

lna = EarthLocation('-45 34 57', '-22 32 04', 1864)

tel= 'Netuno_160'
nep = np.loadtxt('ucac4_{}_cp'.format(tel), usecols=[0,1,35,36,43,45], dtype={'names': ('ofra', 'ofde', 'ra', 'dec', 'jd', 'filt'), 'formats': ('f8', 'f8', 'f16', 'f16', 'f16', 'S10')})

ra = nep['ra']*u.hourangle
dec = nep['dec']*u.deg
tempo = Time(nep['jd'], format='jd', scale='utc')
filt = np.char.array(nep['filt'])

hourangle = time_hourangle(tempo, ra)

filtros = ['clear', 'b', 'v', 'r', 'i', 'metano']

valfa = np.zeros((len(filtros), len(nep['ofra'])))
vdelta = np.zeros((len(filtros), len(nep['ofde'])))

va, vd = refraction(dec, hourangle)

for i in np.arange(len(filtros)):
    a = np.where(filt.lower() == filtros[i])
    valfa[i,a[0]] = va[a[0]]
    vdelta[i,a[0]] = vd[a[0]]
    
    
jdn, netanom = np.loadtxt('Netuno.eph', skiprows=3, usecols=[2, 36], unpack=True)
jdt, trianom = np.loadtxt('Triton.eph', skiprows=3, usecols=[2, 36], unpack=True)

n,m = np.indices((len(jdn), len(tempo)))
o = np.where(np.absolute(jdn[n]-tempo[m].jd)*u.day < 0.01*u.s)
anomnet = netanom[o[0]]

n,m = np.indices((len(jdt), len(tempo)))
o = np.where(np.absolute(jdt[n]-tempo[m].jd)*u.day < 0.01*u.s)
anomtri = trianom[o[0]]

bin = np.arange(-380,400,40)

############## Funcoes de ajuste #####################################

g1 = np.vstack((valfa, np.sin(anomnet*u.deg), np.cos(anomnet*u.deg), np.sin(anomtri*u.deg), np.cos(anomtri*u.deg), np.ones(len(nep['ofra'])))).T

def ff(B,x):
    return np.sum([B[i]*x[i] for i in np.arange(len(filtros))], axis=0)

def f1(B,x):
    return np.sum([B[i]*x[i] for i in np.arange(len(filtros))], axis=0) + B[len(filtros)]*np.sin(x[len(filtros)]*u.deg) + B[len(filtros)+1]*np.cos(x[len(filtros)]*u.deg) + B[len(filtros)+2]*np.sin(x[len(filtros)+1]*u.deg) + B[len(filtros)+3]*np.cos(x[len(filtros)+1]*u.deg) + B[len(filtros)+4]


################  EXEMPLOS  ##########################################

#g1 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

#g2 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

#g3 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg)**2, np.cos(y[3]*u.deg)**2, np.sin(y[3]*u.deg)*np.cos(y[3]*u.deg), np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

#g4 = np.vstack([y[0] - 2451544.5, np.sin(y[3]*u.deg), (y[0] - 2451544.5)*np.sin(y[3]*u.deg), np.cos(y[3]*u.deg), (y[0] - 2451544.5)*np.cos(y[3]*u.deg), np.ones(len(y[3]))]).T

#def f1(B, x): ## tempo, seno, cosseno, constante
#    return (B[0]/365)*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg) + B[2]*np.cos(x[1]*u.deg) + B[3]
    
#def f2(B, x): ## tempo, seno*cos, sen, cos, constante
#    return B[0]*(x[0] - 2451544.5) + B[1]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[2]*np.sin(x[1]*u.deg) + B[3]*np.cos(x[1]*u.deg) + B[4]
    
#def f3(B, x): ## tempo, sen^2, cos^2, seno*cos, sen, cos, constante
#    return B[0]*(x[0] - 2451544.5) + B[1]*(np.sin(x[1]*u.deg)**2) + B[2]*(np.cos(x[1]*u.deg)**2) + B[3]*np.sin(x[1]*u.deg)*np.cos(x[1]*u.deg) + B[4]*np.sin(x[1]*u.deg) + B[5]*np.cos(x[1]*u.deg) + B[6]
    
#def f4(B, x): ## tempo, sen^2(f), cos^2(f), seno(f)*cos(f), sen(f), cos(f), constante
#    return B[0]*(x[0] - 2451544.5) + (B[1]+B[2]*(x[0] - 2451544.5))*np.sin(x[1]*u.deg) + (B[3]+B[4]*(x[0] - 2451544.5))*np.cos(x[1]*u.deg) + B[5]

#def f5(B, x): ## a*cos(wt+p), sin(f), cos(f), constante
#    return B[0]*np.cos((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180)) + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg) + B[5]
    
#def f6(B, x): ## a*sen(wt+p), sin(f), cos(f), constante
#    return B[0]*np.sin((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180)) + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg) + B[5]

#def f7(B, x):
#    return np.sin((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180))*(B[0] + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg)) + B[5]

############## least square function #################################

def least(func, x, y, sy=None, beta0=None, ifixb=None):
    linear = odrpack.Model(func)
    #mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
    mydata = odrpack.RealData(x, y, sy=sy)
    myodr = odrpack.ODR(mydata, linear, beta0=beta0, ifixb=ifixb)
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
    
#######################################################################

p = np.linalg.lstsq(g1, nep['ofra'])

x = np.vstack((valfa, anomnet, anomtri))

print 'Ascensao Reta\n'

fun=f1
beta0=p[0]

ajranwg = least(func=fun, x=x, y=nep['ofra'], beta0=beta0)

for i in np.arange(len(filtros)):
    print '{:10s}: B={:-6.3f}+-{:5.3f}, n_images={:3d}'.format(filtros[i], ajranwg.beta[i], ajranwg.sd_beta[i], len(np.where(valfa[i] != 0.0)[0]))
    
f, axarr = plt.subplots(2, sharex=True, sharey=True)
axarr[0].plot(tempo.jyear, nep['ofra']*1000, '.')
axarr[1].plot(tempo.jyear, nep['ofra']*1000 - ff(ajranwg.beta, x)*1000, '.')
axarr[0].set_title('Right Ascension - no correction')
axarr[1].set_title('Right Ascension - with correction')
axarr[0].set_ylabel('Offset (mas)')
axarr[1].set_ylabel('Offset (mas)')
axarr[1].set_xlabel('Year')
plt.savefig('RAf_{}.png'.format(tel), dpi=300)
plt.clf()


### Declinacao
    
#p = np.linalg.lstsq(g1, nep['ofde'])

x = np.vstack((vdelta, anomnet, anomtri))

print 'Declinacao\n'

#fun=f1
#beta0=np.ones(len(beta0))
#beta0[:len(filtros)] = ajranwg.beta[:len(filtros)]
#ifixb = np.ones(len(beta0))
#ifixb[:len(filtros)] = 0

#ajdenwg = least(func=fun, x=x, y=nep['ofra'], beta0=beta0, ifixb=ifixb)

#for i in np.arange(len(filtros)):
#    print '{:10s}: B={:-6.3f}, n_images={:3d}'.format(filtros[i], ajdenwg.beta[i], len(np.where(valfa[i] != 0.0)[0]))

f, axarr = plt.subplots(2, sharex=True, sharey=True)
axarr[0].plot(tempo.jyear, nep['ofde']*1000, '.')
axarr[1].plot(tempo.jyear, nep['ofde']*1000 - ff(ajranwg.beta, x)*1000, '.')
axarr[0].set_title('Declination - no correction')
axarr[1].set_title('Declination - with correction')
axarr[0].set_ylabel('Offset (mas)')
axarr[1].set_xlabel('Year')
plt.savefig('DECf_{}.png'.format(tel), dpi=300)
plt.clf()

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

lim = 250

ax1.hist(nep['ofra']*1000,bin)
ax1.set_title('RA - no corr')
ax1.set_ylabel('Number of positions')
ax1.set_ylim(0,lim)

ax2.hist(nep['ofde']*1000,bin)
ax2.set_title('DEC - no corr')

ax3.hist((nep['ofra'] - ff(ajranwg.beta, x))*1000,bin)
ax3.set_title('RA - corr')
ax3.set_ylabel('Number of positions')
ax3.set_xlabel('Offsets (mas)')
ax3.set_ylim(0,lim)

ax4.hist((nep['ofde'] - ff(ajranwg.beta, x))*1000,bin)
ax4.set_title('DEC - corr')
ax4.set_xlabel('Offsets (mas)')

plt.savefig('dist_{}.png'.format(tel), dpi=300)
plt.clf()

