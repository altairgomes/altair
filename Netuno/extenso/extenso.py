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
    def func_integer(mi, i):
        erf = special.erf((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))-special.erf((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))
        return np.exp(-(x[0][i]-B[2]-mi)*(x[0][i]-B[2]-mi)/(2*B[1]*B[1]))*erf
    integral = np.array([integrate.quad(func_integer,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    return B[0]*np.sqrt(np.pi/2)*B[1]*integral + B[5]
#    A, sigma, x0, y0, R, c

#def extenso2(B,x):
#    def func_integer(mi):
#        return np.exp(-(x[0]-B[2]-mi)*(x[0]-B[2]-mi)/2*B[1]*B[1])*(special.erf((x[1]-B[3]+np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))-special.erf((x[1]-B[3]-np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1])))
#    integral = integrate.quad(func_integer,-B[4], B[4])
#    return B[0]*np.sqrt(np.pi/2)*B[1]*integral + B[5]

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
c = np.where(np.sqrt((a-912)**2 + (b-1060)**2) < 10)
d,e = np.meshgrid(np.arange(c[1].min(), c[1].max()), np.arange(c[0].min(), c[0].max()))

#plt.imshow(image[e,d], cmap='gray')


fi = open('saida.txt', 'w')
x = np.array(c)
y = image[c]

t = Time.now()
beta0 = np.array([20000.0,912.0,1060.0,1.9, 100.0])
aj = least(func=gaussian, x=x, y=y, beta0=beta0)
fi.write('gaussiana: {} sec, A={}+-{}, x0={}+-{}, y0={}+-{}, sigma={}, c={}+-{}\n\n'.format((Time.now()-t).sec, aj.beta[0], aj.sd_beta[0], aj.beta[1]+1, aj.sd_beta[1], aj.beta[2]+1, aj.sd_beta[2], aj.beta[3], aj.beta[4], aj.sd_beta[4]))

t = Time.now()
beta0 = np.array([20000.0,1.9,912.0,1060.0,6.48, 100.0]) # A, sigma, x0, y0, R, c
aje = least(func=extenso, x=x, y=y, beta0=beta0)
fi.write('extenso: {} sec, A={}+-{}, x0={}+-{}, y0={}+-{}, sigma={}, R={}+-{}, c={}+-{}\n\n'.format((Time.now()-t).sec, aje.beta[0], aje.sd_beta[0], aje.beta[2]+1, aje.sd_beta[2], aje.beta[3]+1, aje.sd_beta[3], aje.beta[1], aje.beta[4], aje.sd_beta[4], aje.beta[5], aje.sd_beta[5]))

fi.write('extenso (testes)\n')
for i in np.arange(0.8, 5.0, 0.2):
    c = np.where(np.sqrt((a-aje.beta[2])**2 + (b-aje.beta[3])**2) < aje.beta[1]*i)
    x = np.array(c)
    y = image[c]
    t = Time.now()
    beta0 = aje.beta # A, sigma, x0, y0, R, c
    ajt = least(func=extenso, x=x, y=y, beta0=beta0)
    fi.write('{} {} {} {:7.1f} {:6.1f} {:7.2f} {:4.2f} {:7.2f} {:4.2f} {:4.1f} {:5.2f} {:4.2f} {:7.1f} {}\n'.format((Time.now()-t).sec, len(y), i, ajt.beta[0], ajt.sd_beta[0], ajt.beta[2]+1, ajt.sd_beta[2], ajt.beta[3]+1, ajt.sd_beta[3], ajt.beta[1], ajt.beta[4], ajt.sd_beta[4], ajt.beta[5], ajt.sd_beta[5]))
    f = np.sqrt((c[0]-aje.beta[2])**2 + (c[1]-aje.beta[3])**2)
    g = np.argsort(f)
    plt.plot(f[g], y[g], 'o')
    #plt.plot(f[g], gaussian(aj.beta, x)[g], label='gaussiana')
    plt.plot(f[g], extenso(ajt.beta, x)[g], label='extenso')
    plt.legend()
    plt.savefig('teste_sig_{}.png'.format(i), dpi=300)
    plt.clf()

fi.close()

#c = np.where(np.sqrt((a-aje.beta[1])**2 + (b-aje.beta[2])**2) < 10)+1
#x = np.array(c)
#f = np.sqrt((c[0]-aje.beta[2])**2 + (c[1]-aje.beta[3])**2)
#g = np.argsort(f)
#plt.plot(f[g], y[g], 'o')

#plt.plot(f[g], gaussian(aj.beta, x)[g], label='gaussiana')
#plt.plot(f[g], extenso(aje.beta, x)[g], label='extenso')
#plt.legend()
#plt.show()

