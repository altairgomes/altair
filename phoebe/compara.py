from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from observation_3 import read2

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

a = 'Phoebe_eme.eph'
b = 'Phoebe_JPL.eph'
m = 'Phoebe_Josselin.eph'

c = read2(a, time_col=[0,1,2,3,4,5], time_fmt='iso', coord_col=[6,7], skiprows=8)
#c = read2(a, time_col=[2], time_fmt='jd', coord_col=[4,5,6,7,8,9], skiprows=3)
d = read2(b, time_col=[2], time_fmt='jd', coord_col=[4,5,6,7,8,9], skiprows=3)
n = read2(m, time_col=[2], time_fmt='jd', coord_col=[4,5,6,7,8,9], skiprows=3)

#e = read2('Nereid.dat', time_col=[8], time_fmt='jd', coord_col=[0,1,2,3,4,5])
#f = read2('Nereida_LNA.eph', time_col=[0], time_fmt='jd', coord_col=[1,2,3,4,5,6], skiprows=8)

eme_time = Time(c['times'], format='iso', scale='utc')
eme_coord = SkyCoord(c['coords'].transpose(), unit=(u.hourangle, u.deg))

jpl_time = Time(d['times'], format='jd', scale='utc')
jpl_coord = SkyCoord(d['coords'].transpose(), unit=(u.hourangle, u.deg))

ph15_time = Time(n['times'], format='jd', scale='utc')
ph15_coord = SkyCoord(n['coords'].transpose(), unit=(u.hourangle, u.deg))

dalfa = jpl_coord.ra*np.cos(jpl_coord.dec) - ph15_coord.ra*np.cos(ph15_coord.dec)
ddelta = jpl_coord.dec - ph15_coord.dec

dalfaz = eme_coord.ra*np.cos(eme_coord.dec) - ph15_coord.ra*np.cos(ph15_coord.dec)
ddeltaz = eme_coord.dec - ph15_coord.dec

#pos_time = Time(e['times'], format='jd', scale='utc')
#pos_coord = SkyCoord(e['coords'].transpose(), unit=(u.hourangle, u.deg))
#lna_time = Time(f['times'], format='jd', scale='utc')
#lna_coord = SkyCoord(f['coords'].transpose(), unit=(u.hourangle, u.deg))
#dalfapos = pos_coord.ra*np.cos(pos_coord.dec) - lna_coord.ra*np.cos(lna_coord.dec)
#ddeltapos = pos_coord.dec - lna_coord.dec


r = []
for i in np.arange(2015,2019,1):
    r.append(Time('{}-01-01 00:00:00'.format(i), format='iso').jd - 2451544.5)

r = np.array(r)

plt.plot(jpl_time.jd - 2451544.5, ddelta.mas, color='green', label=r'$\Delta\delta$')
plt.plot(jpl_time.jd - 2451544.5, dalfa.mas, color='blue', label=r'$\Delta\alpha\cos\delta$')

#plt.plot(pos_time.jd - 2451544.5, dalfapos.mas, 'r+')
#plt.plot(pos_time.jd - 2451544.5, ddeltapos.mas, 'rx')

plt.xlim(Time('2015-01-01 00:00:00').jd - 2451544.5, Time('2018-01-01 00:00:00').jd - 2451544.5)
plt.ylim(-30,30)
plt.title('SAT375 - PH15', fontsize=sizel)
plt.xlabel('Time', fontsize=sizel)
plt.xticks(r, ['{}'.format(i) for i in np.arange(2015,2019,1)])
plt.ylabel('Difference (mas)', fontsize=sizel)
plt.legend()
plt.axhline(0, color='black')
plt.legend(labelspacing=0.25, borderpad=0.5, handlelength=1.7, prop={'size':sizel})
plt.tick_params(axis='both', which='major', labelsize=sizel)
fig =plt.gcf()
fig.set_size_inches((17.6*u.cm).to(u.imperial.inch).value,(9.9*u.cm).to(u.imperial.inch).value)
fig.savefig('Phoebe_JPL-ph15.eps',dpi=300, format='eps', bbox_inches='tight')
plt.clf()
