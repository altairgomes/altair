import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

typ = {'160': 'b,', 'IAG': 'g,', 'ZEI': 'r,'}

r = []
for i in np.arange(1992,2017,3):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

def plot_todos():
    f = open('lista_ucac4', 'r')
    lista = f.readlines()
    f.close()

    f, axarr = plt.subplots(2, sharex=True)
    
    typ = {'160': 'b.', 'IAG': 'g.', 'ZEI': 'r.'}

    for i in lista:
        x, y, jd = np.loadtxt(i.strip(), usecols=[0, 1, 43], unpack=True)
        tempo = Time(jd, format='jd')
        axarr[0].plot(tempo.jd - 2451544.5, x, typ[i[0:3]])
        axarr[1].plot(tempo.jd - 2451544.5, y, typ[i[0:3]])
    
    axarr[0].set_title('Right Ascension')
    axarr[0].set_ylabel('Offset (arcsec)')
    axarr[1].set_title('Declination')
    axarr[1].set_xticks(r)
    axarr[1].set_xlabel('Year')
    axarr[1].set_ylabel('Offset (arcsec)')
    axarr[1].set_xticklabels(['{}'.format(i) for i in np.arange(1992,2017,3)])

    plt.savefig('Triton_todos.png',dpi=300, format='png', bbox_inches='tight')
    
def plot_media():
    f = open('lista_table', 'r')
    lista = f.readlines()
    f.close()

    f, axarr = plt.subplots(2, sharex=True)

    for i in lista:
        x, y, ex, ey, jd = np.genfromtxt(i.strip(), skip_header=43, skip_footer=1,  usecols=[0, 1, 2, 3, 25], unpack=True)
        tempo = Time(jd, format='jd')
        axarr[0].errorbar(tempo.jd - 2451544.5, x/1000., yerr=ex/1000., fmt=typ[i[0:3]])
        axarr[1].errorbar(tempo.jd - 2451544.5, y/1000., yerr=ey/1000., fmt=typ[i[0:3]])
    
    axarr[0].set_title('Right Ascension')
    axarr[0].set_ylabel('Offset (arcsec)')
    axarr[1].set_title('Declination')
    axarr[1].set_xticks(r)
    axarr[1].set_xlabel('Year')
    axarr[1].set_ylabel('Offset (arcsec)')
    axarr[1].set_xticklabels(['{}'.format(i) for i in np.arange(1992,2017,3)])

    plt.savefig('Triton_media.png',dpi=300, format='png', bbox_inches='tight')
    
    
plot_todos()
plot_media()
