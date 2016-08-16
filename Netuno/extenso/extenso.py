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

class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]

def gaussian(B,x):
    return B[0]*np.exp(-((x[0]-B[1])*(x[0]-B[1]) + (x[1]-B[2])*(x[1]-B[2]))/(2*B[3]*B[3])) + B[4]

def extenso(B,x):
    def func_integer(mi, i):
        erf = special.erf((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))-special.erf((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))
        return np.exp(-(x[0][i]-B[2]-mi)*(x[0][i]-B[2]-mi)/(2*B[1]*B[1]))*erf
    integral = np.array([integrate.quad(func_integer,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    return B[0]*np.sqrt(np.pi/2)*B[1]*integral + B[5]
#    A, sigma, x0, y0, R, c
def func_z(x, B):
    return x/(np.sqrt(2)*B)
    
def func_w1(x,x0,R,mi,s):
    return ((x-x0+np.sqrt(R*R-mi*mi)))/(np.sqrt(2)*s)

def func_w2(x,x0,R,mi,s):
    return ((x-x0-np.sqrt(R*R-mi*mi)))/(np.sqrt(2)*s)

def func_erf(x,x0,R,mi,s):
    w1 = func_w1(x,x0,R,mi,s)
    w2 = func_w2(x,x0,R,mi,s)
    erf = special.erf(w1)-special.erf(w2)
    return erf

func_z = Memoize(func_z)
func_w1 = Memoize(func_w1)
func_w2 = Memoize(func_w2)
func_erf = Memoize(func_erf)

def derivadas(B,x):
    def integ1(mi,i):
        z = (x[0][i] - B[2] - mi)/(np.sqrt(2)*B[1])
        erf = func_erf(x[1][i], B[3], B[4], mi, B[1])
#        w1 = ((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
#        w2 = ((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
#        erf = special.erf(w1)-special.erf(w2)
        return np.exp(-z*z)*erf
    def integ2(mi,i):
        z = (x[0][i] - B[2] - mi)/(np.sqrt(2)*B[1])
        erf = func_erf(x[1][i], B[3], B[4], mi, B[1])
#        w1 = ((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
#        w2 = ((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
#        erf = special.erf(w1)-special.erf(w2)
        return z*np.exp(-z*z)*erf
    def integ3(mi,i):
        z = (x[0][i] - B[2] - mi)/(np.sqrt(2)*B[1])
        w1 = ((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
        w2 = ((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
        erf = np.exp(-w1*w1)-np.exp(-w2*w2)
        return np.exp(-z*z)*erf
    def integ4(mi,i):
        z = (x[0][i] - B[2] - mi)/(np.sqrt(2)*B[1])
        erf = func_erf(x[1][i], B[3], B[4], mi, B[1])
#        w1 = ((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
#        w2 = ((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
#        p1 = np.sqrt(np.pi/2)*(1+2*z*z)*(special.erf(w1)-special.erf(w2))
        p1 =  np.sqrt(np.pi/2)*(1+2*z*z)*erf
        p2 = np.sqrt(2)*(w1*np.exp(-w1*w1)-w2*np.exp(-w2*w2))
        return np.exp(-z*z)*(p1-p2)
    def integ5(mi,i):
        z = (x[0][i] - B[2] - mi)/(np.sqrt(2)*B[1])
        w1 = ((x[1][i]-B[3]+np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
        w2 = ((x[1][i]-B[3]-np.sqrt(B[4]*B[4]-mi*mi)))/(np.sqrt(2)*B[1])
        erf = np.exp(-w1*w1)+np.exp(-w2*w2)
        return (B[4]/np.sqrt(B[4]*B[4]-mi*mi))*np.exp(-z*z)*erf
    integral1 = np.array([integrate.quad(integ1,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    integral2 = np.array([integrate.quad(integ2,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    integral3 = np.array([integrate.quad(integ3,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    integral4 = np.array([integrate.quad(integ4,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    integral5 = np.array([integrate.quad(integ5,-B[4], B[4], args=[i])[0] for i in np.arange(len(x[0]))])
    dFdA = np.sqrt(np.pi/2)*B[1]*integral1
    dFdx0 = B[0]*np.sqrt(np.pi)*integral2
    dFdy0 = -B[0]*integral3
    dFds = B[0]*integral4
    dFdR = -B[0]*integral5
    dFdc = np.ones(len(x[0]))
    return np.array([dFdA, dFds, dFdx0, dFdy0, dFdR, dFdc])

#def extenso2(B,x):
#    def func_integer(mi):
#        return np.exp(-(x[0]-B[2]-mi)*(x[0]-B[2]-mi)/2*B[1]*B[1])*(special.erf((x[1]-B[3]+np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1]))-special.erf((x[1]-B[3]-np.sqrt(B[4]*B[4]-mi*mi))/(np.sqrt(2)*B[1])))
#    integral = integrate.quad(func_integer,-B[4], B[4])
#    return B[0]*np.sqrt(np.pi/2)*B[1]*integral + B[5]

def least(func, x, y, sy=None, beta0=None, ifixb=None, fjacb=None):
    linear = odrpack.Model(func, fjacb=fjacb)
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


#fi = open('saida.txt', 'w')
x = np.array(c)
y = image[c]

t = Time.now()
betag = np.array([  3.82898171e+04,   9.11453238e+02,   1.05903152e+03,   4.47998911e+00,   1.21856339e+02])
#aj = least(func=gaussian, x=x, y=y, beta0=betag)
#fi.write('gaussiana: {} sec, A={}+-{}, x0={}+-{}, y0={}+-{}, sigma={}, c={}+-{}\n\n'.format((Time.now()-t).sec, aj.beta[0], aj.sd_beta[0], aj.beta[1]+1, aj.sd_beta[1], aj.beta[2]+1, aj.sd_beta[2], aj.beta[3], aj.beta[4], aj.sd_beta[4]))
print (Time.now()-t).sec

t = Time.now()
betae = np.array([  746.81322157,     3.02038928,   911.45256174,  1059.02927882,     5.84291974,   312.22989524]) # A, sigma, x0, y0, R, c
aje = least(func=extenso, x=x, y=y, beta0=betae)
#fi.write('extenso: {} sec, A={}+-{}, x0={}+-{}, y0={}+-{}, sigma={}, R={}+-{}, c={}+-{}\n\n'.format((Time.now()-t).sec, aje.beta[0], aje.sd_beta[0], aje.beta[2]+1, aje.sd_beta[2], aje.beta[3]+1, aje.sd_beta[3], aje.beta[1], aje.beta[4], aje.sd_beta[4], aje.beta[5], aje.sd_beta[5]))
print (Time.now()-t).sec

t = Time.now()
betae = np.array([  746.81322157,     3.02038928,   911.45256174,  1059.02927882,     5.84291974,   312.22989524]) # A, sigma, x0, y0, R, c
aje2 = least(func=extenso, x=x, y=y, beta0=betae, fjacb=derivadas)
#fi.write('extenso: {} sec, A={}+-{}, x0={}+-{}, y0={}+-{}, sigma={}, R={}+-{}, c={}+-{}\n\n'.format((Time.now()-t).sec, aje.beta[0], aje.sd_beta[0], aje.beta[2]+1, aje.sd_beta[2], aje.beta[3]+1, aje.sd_beta[3], aje.beta[1], aje.beta[4], aje.sd_beta[4], aje.beta[5], aje.sd_beta[5]))
print (Time.now()-t).sec

#fi.write('extenso (testes)\n')
#for i in np.arange(0.8, 5.0, 0.2):
#    c = np.where(np.sqrt((a-aje.beta[2])**2 + (b-aje.beta[3])**2) < aje.beta[1]*i)
#    x = np.array(c)
#    y = image[c]
#    t = Time.now()
#    beta0 = aje.beta # A, sigma, x0, y0, R, c
#    ajt = least(func=extenso, x=x, y=y, beta0=beta0)
#    fi.write('{} {} {} {:7.1f} {:6.1f} {:7.2f} {:4.2f} {:7.2f} {:4.2f} {:4.1f} {:5.2f} {:4.2f} {:7.1f} {}\n'.format((Time.now()-t).sec, len(y), i, ajt.beta[0], ajt.sd_beta[0], ajt.beta[2]+1, ajt.sd_beta[2], ajt.beta[3]+1, ajt.sd_beta[3], ajt.beta[1], ajt.beta[4], ajt.sd_beta[4], ajt.beta[5], ajt.sd_beta[5]))
#    f = np.sqrt((c[0]-aje.beta[2])**2 + (c[1]-aje.beta[3])**2)
#    g = np.argsort(f)
#    plt.plot(f[g], y[g], 'o')
#    plt.plot(f[g], gaussian(aj.beta, x)[g], label='gaussiana')
#    plt.plot(f[g], extenso(ajt.beta, x)[g], label='extenso')
#    plt.legend()
#    plt.savefig('teste_sig_{}.png'.format(i), dpi=300)
#    plt.clf()

#fi.close()

#c = np.where(np.sqrt((a-aje.beta[2])**2 + (b-aje.beta[3])**2) < 20)
#x = np.array(c)
f = np.sqrt((x[0]-betae[2])**2 + (x[1]-betae[3])**2)
#g = np.argsort(f)
#plt.plot(f, image[c], 'o')

#beta1=np.array([39438.8794772, 911.447964519, 1059.02996984, 4.80370494377, -1898.25592524])
#beta2=np.array([884.728146958, 2.56992636061, 911.440791197, 1059.0107055, 5.9031003977, 1726.15672919])
plt.plot(f, y - gaussian(betag, x), 'o', label='gaussiana', color='red')
plt.plot(f, y - extenso(betae, x), 'o', label='extenso', color='green')
plt.axhline(0.0)
plt.legend()
plt.savefig('ajuste3.png', dpi=300)

