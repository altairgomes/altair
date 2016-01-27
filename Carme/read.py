from observation_3 import read2
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

#####################################################################

plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]

params = {'text.usetex' : True,
          'font.size' : 22,
          'font.family' : 'txfonts',
          'text.latex.unicode': True,
          }
          
plt.rcParams.update(params)

sizel = 17

#####################################################################

a = ['OPD', 'OHP', 'ESO']
b = [874, 511, 809]

c = read2('Carme.dat', coord_col=[0, 1, 2, 3, 4, 5], time_col=[8], time_fmt='jd')
coord = SkyCoord(c['coords'][0], c['coords'][1], unit=(u.hourangle, u.deg))
time = Time(c['times'], format='jd', scale='utc')
code = np.loadtxt('Carme.dat', usecols=[12], ndmin=1, unpack=True)

for i in [0, 1, 2]:
    d = read2('{}-jpl.eph'.format(b[i]), coord_col=[4, 5, 6, 7, 8, 9], skiprows=3)
    coordjpl = SkyCoord(d['coords'][0], d['coords'][1], unit=(u.hourangle, u.deg))
    anom = np.loadtxt('{}-jpl.eph'.format(b[i]), skiprows=3, usecols=[37], ndmin=1, unpack=True)
    e = read2('{}.eph'.format(b[i]), coord_col=[4, 5, 6, 7, 8, 9], skiprows=3)
    coordlau = SkyCoord(e['coords'][0], e['coords'][1], unit=(u.hourangle, u.deg))
    f = read2('{}-eme.eph'.format(b[i]), coord_col=[1, 2], skiprows=8)
    coordeme = SkyCoord(f['coords'][0], f['coords'][1], unit=(u.hourangle, u.deg))
    difjpldec = coord[code == b[i]].dec - coordjpl.dec
    diflaudec = coord[code == b[i]].dec - coordlau.dec
    difemedec = coord[code == b[i]].dec - coordeme.dec
    j, = plt.plot(anom, difjpldec.mas, '+', color='blue')
    m, = plt.plot(anom, difemedec.mas, '.', color='green')
    l, = plt.plot(anom, diflaudec.mas, 'x', color='red')
    
plt.xlim(0,360)
plt.ylim(-300,300)
plt.axhline(0, color='black')
plt.title('Carme', fontsize=sizel)
plt.xlabel('True Anomaly', fontsize=sizel)
plt.ylabel('Offset (mas)', fontsize=sizel)
plt.xticks(np.arange(0,370,45))
#ax = plt.gca().add_artist(first_legend)
plt.legend([l, j, m], ['STE', 'JUP300', 'Eme2008'], numpoints=1, loc=3, labelspacing=0.25, borderpad=0.1, handlelength=0.5, prop={'size':sizel})
plt.tick_params(axis='both', which='major', labelsize=sizel)
fig =plt.gcf()
fig.set_size_inches((17.6*u.cm).to(u.imperial.inch).value,(9.9*u.cm).to(u.imperial.inch).value)
fig.savefig('Carme_ephemeris.eps', format='eps', dpi=300, bbox_inches='tight')