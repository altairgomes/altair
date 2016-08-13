import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time, TimeDelta
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
    
def cluster(data, maxgap):             
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        >>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    d = data.argsort()
    groups = [[d[0]]]
    for x in d[1:]:
        if abs(data[x] - data[groups[-1][-1]]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

def match():
    nep = np.loadtxt('ucac4_Netuno_ZEI_cp', usecols=[0,1,31,35,36,43,45], dtype={'names': ('ofra', 'ofde', 'nstar', 'ra', 'dec', 'jd', 'filt'), 'formats': ('f8', 'f8', 'i4', 'f8', 'f8', 'f8', 'S10')})
    tri = np.loadtxt('ucac4_Triton_ZEI_cp', usecols=[0,1,31,35,36,43,45], dtype={'names': ('ofra', 'ofde', 'nstar', 'ra', 'dec', 'jd', 'filt'), 'formats': ('f8', 'f8', 'i4', 'f8', 'f8', 'f8', 'S10')})
    n,m = np.indices((len(nep['jd']), len(tri['jd'])))
    o = np.where(nep['jd'][n] == tri['jd'][m])
    difx = tri['ofra'][o[1]] - nep['ofra'][o[0]]
    dify = tri['ofde'][o[1]] - nep['ofde'][o[0]]
    jd = nep['jd'][o[0]]
    f = open('ucac4_Match_ZEI_cp', 'w')
    for i in np.arange(len(o[0])):
        f.write(' {:-6.3f} {:-6.3f}'.format(difx[i], dify[i]) + ' x '*33 + ' {:12.9f} {:-13.9f}'.format(nep['ra'][o[0]][i],nep['dec'][o[0]][i]) + ' x '*6 + '{:16.8f} {} {}\n'.format(jd[i], 'x', nep['filt'][o[0]][i]))
    f.close()

#####################################################################

#match()

lna = EarthLocation('-45 34 57', '-22 32 04', 1864)

tel= 'Netuno_IAG'
nep = np.loadtxt('ucac4_{}_cp'.format(tel), usecols=[0,1,31,35,36,43,45], dtype={'names': ('ofra', 'ofde', 'nstar', 'ra', 'dec', 'jd', 'filt'), 'formats': ('f8', 'f8', 'i4', 'f8', 'f8', 'f8', 'S10')})

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

def fun(B,x):
    return B[0]*x[0] + B[1]*x[1]

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
f = open('refrac_{}.dat'.format(tel), 'w')
groups = cluster(tempo.jd,TimeDelta(10*u.h).jd)
bin = np.arange(-380,400,40)

x,y,cx,cy = np.array([]),np.array([]),np.array([]),np.array([])
for i in groups:
    if hourangle[i[-1]] - hourangle[i[0]] < 1.5*u.hourangle:
        continue
    g = np.vstack((va[i], np.ones(len(i))))
    p = np.linalg.lstsq(g.T, nep['ofra'][i])
    r = least(func=fun, x=g, y=nep['ofra'][i], beta0=p[0])
    t = Time(int('{:8.0f}'.format(tempo[i][0].jd)), format='jd')
    x = np.hstack((x, nep['ofra'][i]))
    y = np.hstack((y, nep['ofde'][i]))
    cx = np.hstack((cx, nep['ofra'][i] - p[0][0]*va[i]))
    cy = np.hstack((cy, nep['ofde'][i] - p[0][0]*vd[i]))
    cora = nep['ofra'][i] - p[0][0]*va[i]
    cord = nep['ofde'][i] - p[0][0]*vd[i]
    f.write('{} & {:5s} & {:4.2f} & {:+5.2f}$\pm${:4.2f} & {:3d} & {:3.0f} & {:-4.0f}$\pm${:3.0f} & {:-4.0f}$\pm${:3.0f} & {:-4.0f}$\pm${:3.0f} & {:-4.0f}$\pm${:3.0f} \\\\ \n'.format(t.iso.split(' ')[0], nep['filt'][i][0], (hourangle[i[-1]] - hourangle[i[0]]).value, r.beta[0], r.sd_beta[0], len(i), nep['nstar'][i].mean(), nep['ofra'][i].mean()*1000, nep['ofra'][i].std()*1000, nep['ofde'][i].mean()*1000, nep['ofde'][i].std()*1000, cora.mean()*1000, cora.std()*1000, cord.mean()*1000, cord.std()*1000))
    print '{}: Delta_H={:5.3f}; B={:+6.3f}, off={:-4.0f}, Ni={:3d}, ncora={:-4.0f}+-{:3.0f}, cora={:-4.0f}+-{:3.0f}, ncord={:-4.0f}+-{:3.0f}, cord={:-4.0f}+-{:3.0f}'.format(t.iso.split(' ')[0], (hourangle[i[-1]] - hourangle[i[0]]).value, p[0][0], p[0][1]*1000, len(i), nep['ofra'][i].mean()*1000, nep['ofra'][i].std()*1000, cora.mean()*1000, cora.std()*1000, nep['ofde'][i].mean()*1000, nep['ofde'][i].std()*1000, cord.mean()*1000, cord.std()*1000)
    plt.plot(hourangle[i], nep['ofra'][i], 'o', label='no cor')
    plt.plot(hourangle[i], nep['ofra'][i] - p[0][0]*va[i], 'o', label='cor')
    plt.xlabel('Hourangle')
    plt.ylabel('Offset (arcsec)')
    plt.legend(numpoints=1)
    plt.savefig('{}_{}.png'.format(tel, t.iso.split(' ')[0]))
    plt.clf()
f.close()


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

ax1.hist(x*1000,bin)
ax1.set_title('RA - no corr')
ax1.set_ylabel('Number of positions')
#ax1.set_ylim(0,450)

ax2.hist(y*1000,bin)
ax2.set_title('DEC - no corr')

ax3.hist(cx*1000,bin)
ax3.set_title('RA - corr')
ax3.set_ylabel('Number of positions')
ax3.set_xlabel('Offsets (mas)')
#ax3.set_ylim(0,450)

ax4.hist(cy*1000,bin)
ax4.set_title('DEC - corr')
ax4.set_xlabel('Offsets (mas)')

plt.savefig('dist_{}.png'.format(tel), dpi=300)
plt.clf()
    
#######################################################################

#p = np.linalg.lstsq(g1, nep['ofra'])
#
#x = np.vstack((valfa, anomnet, anomtri))
#
#print 'Ascensao Reta\n'
#
#fun=f1
#beta0=p[0]
#
#ajranwg = least(func=fun, x=x, y=nep['ofra'], beta0=beta0)
#
#for i in np.arange(len(filtros)):
#    print '{:10s}: B={:-6.3f}+-{:5.3f}, n_images={:3d}'.format(filtros[i], ajranwg.beta[i], ajranwg.sd_beta[i], len(np.where(valfa[i] != 0.0)[0]))
#    
#f, axarr = plt.subplots(2, sharex=True, sharey=True)
#axarr[0].plot(tempo.jyear, nep['ofra']*1000, '.')
#axarr[1].plot(tempo.jyear, nep['ofra']*1000 - ff(ajranwg.beta, x)*1000, '.')
#axarr[0].set_title('Right Ascension - no correction')
#axarr[1].set_title('Right Ascension - with correction')
#axarr[0].set_ylabel('Offset (mas)')
#axarr[1].set_ylabel('Offset (mas)')
#axarr[1].set_xlabel('Year')
#plt.savefig('RAf_{}.png'.format(tel), dpi=300)
#plt.clf()
#
#
#### Declinacao
#    
##p = np.linalg.lstsq(g1, nep['ofde'])
#
#x = np.vstack((vdelta, anomnet, anomtri))
#
#print 'Declinacao\n'
#
##fun=f1
##beta0=np.ones(len(beta0))
##beta0[:len(filtros)] = ajranwg.beta[:len(filtros)]
##ifixb = np.ones(len(beta0))
##ifixb[:len(filtros)] = 0
#
##ajdenwg = least(func=fun, x=x, y=nep['ofra'], beta0=beta0, ifixb=ifixb)
#
##for i in np.arange(len(filtros)):
##    print '{:10s}: B={:-6.3f}, n_images={:3d}'.format(filtros[i], ajdenwg.beta[i], len(np.where(valfa[i] != 0.0)[0]))
#
#f, axarr = plt.subplots(2, sharex=True, sharey=True)
#axarr[0].plot(tempo.jyear, nep['ofde']*1000, '.')
#axarr[1].plot(tempo.jyear, nep['ofde']*1000 - ff(ajranwg.beta, x)*1000, '.')
#axarr[0].set_title('Declination - no correction')
#axarr[1].set_title('Declination - with correction')
#axarr[0].set_ylabel('Offset (mas)')
#axarr[1].set_xlabel('Year')
#plt.savefig('DECf_{}.png'.format(tel), dpi=300)
#plt.clf()
#
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#
#lim = 250
#
#ax1.hist(nep['ofra']*1000,bin)
#ax1.set_title('RA - no corr')
#ax1.set_ylabel('Number of positions')
#ax1.set_ylim(0,lim)
#
#ax2.hist(nep['ofde']*1000,bin)
#ax2.set_title('DEC - no corr')
#
#ax3.hist((nep['ofra'] - ff(ajranwg.beta, x))*1000,bin)
#ax3.set_title('RA - corr')
#ax3.set_ylabel('Number of positions')
#ax3.set_xlabel('Offsets (mas)')
#ax3.set_ylim(0,lim)
#
#ax4.hist((nep['ofde'] - ff(ajranwg.beta, x))*1000,bin)
#ax4.set_title('DEC - corr')
#ax4.set_xlabel('Offsets (mas)')
#
#plt.savefig('dist_{}.png'.format(tel), dpi=300)
#plt.clf()
#
