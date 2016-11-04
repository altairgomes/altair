import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u


backfw = True

t, ene, dene = np.loadtxt('energy.dat', usecols=[0,1,2], unpack=True, dtype=np.float64)

time = Time(t, format='jd', scale='tdb')

plt.plot(time.jyear, dene)
plt.ylabel(r'$\Delta E/E$')
plt.xlabel('Time (years)')
plt.savefig('energy.png', dpi=300)
plt.clf()

x = np.loadtxt('output.dat', unpack=True)
time = Time(x[0], format='jd', scale='tdb')

nbody = (len(x)-1)/3

for i in np.arange(nbody):
    d = np.sqrt(x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2]+x[3*i+3]*x[3*i+3])
    d = d*u.AU
    plt.plot(time.jyear, d)
    plt.ylabel('Distance (AU)')
    plt.xlabel('Time (years)')
    plt.title('Object {}'.format(i+1))
    plt.savefig('object_{}.png'.format(i+1), dpi=300)
    plt.clf()
        
if backfw:
    t, ene, dene = np.loadtxt('energy.dat', unpack=True, dtype=np.float64)
    t2, ene2, dene2 = np.loadtxt('energy_back.dat', unpack=True, dtype=np.float64)
    time = Time(t, format='jd', scale='tdb')
    
    deltaen = dene2[::-1] - dene
    
    plt.plot(time.jyear, deltaen)
    plt.ylabel(r'$\Delta E2/E - \Delta E1/E$')
    plt.axhline(0)
    plt.xlabel('Time (years)')
    plt.savefig('energy_bf.png', dpi=300)
    plt.clf()

    x = np.loadtxt('output.dat', unpack=True, dtype=np.float64)
    x2 = np.loadtxt('output_back.dat', unpack=True, dtype=np.float64)
    time = Time(x[0], format='jd', scale='tdb')

    nbody = (len(x)-1)/3

    for i in np.arange(nbody):
#        delta = x2[3*i+1:3*i+3]-x[3*i+1:3*i+3]
        delta = np.array([x2[3*i+1,::-1]-x[3*i+1,:], x2[3*i+2,::-1]-x[3*i+2,:], x2[3*i+3,::-1]-x[3*i+3,:]])
        d = np.sqrt(np.sum(delta*delta, axis=0))
#        d = np.sqrt(x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2]+x[3*i+3]*x[3*i+3])
#        d2 = np.sqrt(x2[3*i+1]*x2[3*i+1]+x2[3*i+2]*x2[3*i+2]+x2[3*i+3]*x2[3*i+3])
        d = d*u.AU
#        d2 = d2*u.AU
#        deltad = d2[::-1] - d
        plt.plot(time.jyear, d.to(u.km))
        plt.axhline(0)
        plt.ylabel('Delta Distance (km)')
        plt.xlabel('Time (years)')
        plt.title('Object {}'.format(i+1))
        plt.savefig('object_{}_bf.png'.format(i+1), dpi=300)
        plt.clf()
