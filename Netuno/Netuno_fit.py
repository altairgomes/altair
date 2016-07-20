import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import scipy.odr.odrpack as odrpack

##################################### FUNCOES #######################

def time_hourangle(time, ra):
    time.delta_ut1_utc = 0.0
    time.location = lna
    hourangle = time.sidereal_time('mean') - ra
    hourangle.wrap_at('180d', inplace=True)
    return hourangle
    
def refraction(delta, hourangle):
    den = np.tan(delta)*np.tan(lna.latitude)+np.cos(hourangle0)
    valfa = np.sin(hourangle)/(den*np.cos(delta))
    vdelta = (np.tan(lna.latitude)-np.tan(delta)*np.cos(hourangle))/den
    return valfa, vdelta
    


#####################################################################

lna = EarthLocation('-45 34 57', '-22 32 04', 1864)

nep = np.loadtxt('ucac4_Netuno_160_cp', usecols=[0,1,35,36,43,45], dtype={'names': ('ofra', 'ofde', 'ra', 'dec', 'jd', 'filt'), 'formats': ('f8', 'f8', 'f16', 'f16', 'f16', 'S10')})

ra = nep['ra']*u.hourangle
dec = nep['dec']*u.deg
tempo = Time(nep['jd'], format='jd', scale='utc')




############## Funcoes de ajuste #####################################

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
    return np.sin((2*np.pi/(B[1]*365.25))*(x[0] - 2451544.5) + B[2]*(np.pi/180))*(B[0] + B[3]*np.sin(x[1]*u.deg) + B[4]*np.cos(x[1]*u.deg)) + B[5]

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
    
#######################################################################

p = np.linalg.lstsq(g1, z[1])

x = np.vstack([eph[0], eph[1], eph[5]])

print 'Declinacao\n'
f.write('\nDeclinacao\n')

fun=f7
beta0=[200.0, 11.7, 0.0, 1.0, 1.0, 0.0, 1.0]

ajdenwg = least(func=fun, x=x, y=eph[7], beta0=beta0)
