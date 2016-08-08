import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation
import scipy.odr.odrpack as odrpack
from astropy.io import fits
import scipy.integrate as integrate
import scipy.special as special

##########################################################################

def gaussian(B,x):
    return B[0]*np.exp(-((x[0]-B[1])*(x[0]-B[1]) + (x[1]-B[2])*(x[1]-B[2]))/(2*B[3]*B[3])) + B[4]

def extenso(B,x):
    def func_integer(mi):
        np.exp(-(x[0]-B[2]-mi)*(x[0]-B[2]-mi)/2*B[1]*B[1])*(special.erf((x[1]-B[3]+np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))-special.erf((x[1]-B[3]-np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1])))
    integral = integrate.quad(func_integer,-B[4], B[4])
    return B[0]*np.sqrt(np.pi/2)*B[1]*integral + B[5]
#    A, sigma, x0, y0, R


def least(func, x, y, sy=None, beta0=None, ifixb=None):
    linear = odrpack.Model(func)
    #mydata = odrpack.Data(x, z[1], wd=1./np.power(z[2],2), we=1./np.power(sy,2))
    mydata = odrpack.RealData(x, y, sy=sy)
    myodr = odrpack.ODR(mydata, linear, beta0=beta0, ifixb=ifixb)
    myodr.set_job(fit_type=2)
    myoutput = myodr.run()
    myoutput.pprint()
    return myodr.output
    
###########################################################################

image = fits.getdata('Netuno_08.21.06_0001.fits')
a, b = np.indices(image.shape)
c = np.where(np.sqrt((a-912)**2 + (b-1060)**2))
d,e = np.meshgrid(np.arange(c[1].min(), c[1].max()), np.arange(c[0].min(), c[0].max()))

#plt.imshow(image[e,d], cmap='gray')



x = np.array(c)
y = image[c]
beta0 = np.array([20000.0,912.0,1060.0,1.9, 100.0])

aj = least(func=gaussian, x=x, y=y, beta0=beta0)
#aje = least(func=extenso, x=x, y=y, beta0=beta0)

f = np.sqrt((c[0]-aj.beta[1])**2 + (c[1]-aj.beta[2])**2)
plt.plot(f, y, 'o')
plt.plot(f, gaussian(aj.beta, x))
plt.show()
