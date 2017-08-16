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
obj = 'Elara'
a = ['OPD', 'OHP', 'ESO']
b = [874, 511, 809]

c = read2('{}.dat'.format(obj), coord_col=[0, 1, 2, 3, 4, 5], time_col=[8], time_fmt='jd')
coord = SkyCoord(c['coords'][0], c['coords'][1], unit=(u.hourangle, u.deg))
time = Time(c['times'], format='jd', scale='utc')
code = np.loadtxt('{}.dat'.format(obj), usecols=[12], ndmin=1, unpack=True)

for i in [0, 1, 2]:
    d = read2('{}-jpl.eph'.format(b[i]), coord_col=[4, 5, 6, 7, 8, 9], skiprows=3)
    coordjpl = SkyCoord(d['coords'][0], d['coords'][1], unit=(u.hourangle, u.deg))
    anom = np.loadtxt('{}-jpl.eph'.format(b[i]), skiprows=3, usecols=[37], ndmin=1, unpack=True)
    e = read2('{}-ste.eph'.format(b[i]), coord_col=[4, 5, 6, 7, 8, 9], skiprows=3)
    coordlau = SkyCoord(e['coords'][0], e['coords'][1], unit=(u.hourangle, u.deg))
    f = read2('{}-eme.eph'.format(b[i]), coord_col=[1, 2], skiprows=8)
    coordeme = SkyCoord(f['coords'][0], f['coords'][1], unit=(u.hourangle, u.deg))
    difjpldec = coord[code == b[i]].dec - coordjpl.dec
    difjplra = (coord[code == b[i]].ra - coordjpl.ra)*np.cos(coordjpl.dec)
    f = open('{}_jpl.txt'.format(obj), 'a')
    for j in np.arange(len(difjpldec)):
        f.write('{:16.8f} {:+04.0f} {:+04.0f}\n'.format(time[code == b[i]].jd[j], difjplra[j].mas, difjpldec[j].mas))
    f.close()
    diflaudec = coord[code == b[i]].dec - coordlau.dec
    diflaura = (coord[code == b[i]].ra - coordlau.ra)*np.cos(coordlau.dec)
    f = open('{}_ste.txt'.format(obj), 'a')
    for j in np.arange(len(diflaudec)):
        f.write('{:16.8f} {:+04.0f} {:+04.0f}\n'.format(time[code == b[i]].jd[j], diflaura[j].mas, diflaudec[j].mas))
    f.close()
    difemedec = coord[code == b[i]].dec - coordeme.dec
    difemera = (coord[code == b[i]].ra - coordeme.ra)*np.cos(coordeme.dec)
    f = open('{}_eme.txt'.format(obj), 'a')
    for j in np.arange(len(difemedec)):
        f.write('{:16.8f} {:+04.0f} {:+04.0f}\n'.format(time[code == b[i]].jd[j], difemera[j].mas, difemedec[j].mas))
    f.close()
    j, = plt.plot(time[code == b[i]].jd - 2451544.5, difjpldec.mas, '+', color='blue')
    m, = plt.plot(time[code == b[i]].jd - 2451544.5, difemedec.mas, '.', color='green')
    l, = plt.plot(time[code == b[i]].jd - 2451544.5, diflaudec.mas, 'x', color='red')
    
r = []
for i in np.arange(1995,2016,3):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

r = np.array(r)
    
#plt.xlim(0,360)
plt.xlim(Time('1995-01-01 00:00:00').jd - 2451544.5, Time('2015-01-01 00:00:00').jd - 2451544.5)
plt.ylim(-300,300)
plt.axhline(0, color='black')
plt.title('{}'.format(obj), fontsize=sizel)
plt.xlabel('True Anomaly (degrees)', fontsize=sizel)
plt.ylabel('Offset (mas)', fontsize=sizel)
#plt.xticks(np.arange(0,370,45))
plt.xticks(r, ['{}'.format(i) for i in np.arange(1995,2016,3)])
#ax = plt.gca().add_artist(first_legend)
plt.legend([l, j, m], ['STE', 'JUP300', 'Eme2008'], numpoints=1, loc=1, labelspacing=0.25, borderpad=0.1, handlelength=0.5, prop={'size':sizel})
plt.tick_params(axis='both', which='major', labelsize=sizel)
fig =plt.gcf()
fig.set_size_inches((17.6*u.cm).to(u.imperial.inch).value,(9.9*u.cm).to(u.imperial.inch).value)
fig.savefig('{}_ephemeris_time.eps'.format(obj), format='eps', dpi=300, bbox_inches='tight')
