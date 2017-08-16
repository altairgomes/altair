from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

######################################## lendo arquivo de entrada #########################################################

f = open('mapa_in.dat', 'r')

arquivo = f.readline().rsplit()[0]
obj = f.readline().rsplit()[0]
tamanho = f.readline().rsplit()[0]*u.km
erro = f.readline().rsplit()[0]*u.mas
sitearq = f.readline().rsplit()[0]

f.close()

######################################### lendo arquivo de dados da ocultacao e de observatorios ##############################

dados = np.loadtxt(arquivo, skiprows=41, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19, 20, 21, 22, 25, 26, 28, 29), \
    dtype={'names': ('dia', 'mes', 'ano', 'hor', 'min', 'sec', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'ca', 'pa', 'vel', 'delta', 'mR', 'mK', 'long', 'ora', 'ode'), 
    'formats': ('S30', 'S30', 'S30','S30', 'S30', 'f8','i4', 'i4', 'f8','i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})

if os.path.isfile(sitearq) == True:
    sites = np.loadtxt(sitearq,  dtype={'names': ('lat', 'lon', 'alt', 'nome'), 'formats': ('f8', 'f8', 'f8', 'S30')})

################## lendo coordenadas #################

alfa = np.core.defchararray.add(np.array(map(str, dados['afh'])), [' '])
alfa = np.core.defchararray.add(alfa, np.array(map(str, dados['afm'])))
alfa = np.core.defchararray.add(alfa, [' '])
alfa = np.core.defchararray.add(alfa, np.array(map(str, dados['afs'])))
delta = np.core.defchararray.add(np.array(map(str, dados['ded'])), [' '])
delta = np.core.defchararray.add(delta, np.array(map(str, dados['dem'])))
delta = np.core.defchararray.add(delta, [' '])
delta = np.core.defchararray.add(delta, np.array(map(str, dados['des'])))
coor = np.core.defchararray.add(alfa, [' '])
coor = np.core.defchararray.add(coor, delta)
stars = SkyCoord(coor, frame='icrs', unit=(u.hourangle, u.degree))

################### lendo tempo ########################

tempo = np.core.defchararray.add(np.array(dados['ano']), ['-'])
tempo = np.core.defchararray.add(tempo, np.array(dados['mes']))
tempo = np.core.defchararray.add(tempo, ['-'])
tempo = np.core.defchararray.add(tempo, np.array(dados['dia']))
tempo = np.core.defchararray.add(tempo, [' '])
tempo = np.core.defchararray.add(tempo, np.array(dados['hor']))
tempo = np.core.defchararray.add(tempo, [':'])
tempo = np.core.defchararray.add(tempo, np.array(dados['min']))
tempo = np.core.defchararray.add(tempo, [':'])
tempo = np.core.defchararray.add(tempo, np.array(map(str, dados['sec'])))
datas = Time(tempo, format='iso', scale='utc')
datas.delta_ut1_utc = 0

############### definindo parametros #############
ca = dados['ca']*u.arcsec
pa = dados['pa']*u.deg
vel = dados['vel']*(u.km/u.s)
dist = dados['delta']*u.AU
off_ra = dados['ora']*u.mas
off_de = dados['ode']*u.mas
vec = [-7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7]

################################### definido funcao que imprime o mapa #############################################

def geramapa(idx):
    lon = stars[idx].ra - datas[idx].sidereal_time('mean', 'greenwich')

    m = Basemap(projection='ortho',lat_0=stars[idx].dec.value,lon_0=lon.value,resolution='l')
# draw coastlines, country boundaries, fill continents.
    m.shadedrelief()
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.4)
    m.drawstates(linewidth=0.4)
#m.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
#m.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))
    if os.path.isfile(sitearq) == True:
        xpt,ypt = m(sites['lon'],sites['lat'])
        m.plot(xpt,ypt,'bo')
    CS=m.nightshade(datas[idx].datetime, alpha=0.2)
    a, b =m(lon.value, stars[idx].dec.value)
    a = a*u.m
    b = b*u.m
    dista = (dist[idx].to(u.km)*ca[idx].to(u.rad)).value*u.km
    disterr = (dist[idx].to(u.km)*erro.to(u.rad)).value*u.km
    ax = a + dista*np.sin(pa[idx])
    ax2 = ax + 1000*u.km*vec*np.cos(pa[idx])
    ax3 = ax2 - tamanho/2*np.sin(pa[idx])
    ax4 = ax2 + tamanho/2*np.sin(pa[idx])
    ax5 = a + (dista-disterr)*np.sin(pa[idx]) + 1000*u.km*vec*np.cos(pa[idx])
    ax6 = a + (dista+disterr)*np.sin(pa[idx]) + 1000*u.km*vec*np.cos(pa[idx])
    by = b + dista*np.cos(pa[idx])
    by2 = by - 1000*u.km*vec*np.sin(pa[idx])
    by3 = by2 - tamanho/2*np.cos(pa[idx])
    by4 = by2 + tamanho/2*np.cos(pa[idx])
    by5 = b + (dista-disterr)*np.cos(pa[idx]) - 1000*u.km*vec*np.sin(pa[idx])
    by6 = b + (dista+disterr)*np.cos(pa[idx]) - 1000*u.km*vec*np.sin(pa[idx])
    m.plot(ax,by, 'ro', markersize=20)
    m.plot(ax2.to(u.m),by2.to(u.m), 'ro', markersize=8)
    m.plot(ax3.to(u.m), by3.to(u.m), 'b')
    m.plot(ax4.to(u.m), by4.to(u.m), 'b')
    m.plot(ax5.to(u.m), by5.to(u.m), 'r--')
    m.plot(ax6.to(u.m), by6.to(u.m), 'r--')

    fig = plt.gcf()
    fig.set_size_inches(18.0, 15.0)
    plt.title('-{} D={}- dots each 1000km or {:.2f} <> offsets (mas): {:.1f}, {:.1f}\n'.format(obj, tamanho, np.absolute(1000*u.km/vel[idx]), off_ra[idx].value, off_de[idx].value), fontsize=25, fontproperties='FreeMono')
    plt.xlabel('\n year-m-d    h:m:s UT     ra__dec__J2000__candidate    C/A    P/A    vel   Delta  R*   K*  long\n\
{}  {:02d} {:02d} {:07.4f} {:+02d} {:02d} {:06.3f} {:6.3f} {:6.2f} {:6.2f} {:5.2f} {:5.1f} {:4.1f}  {:3.0f}'
        .format(datas[idx].iso, dados['afh'][idx], dados['afm'][idx], dados['afs'][idx], dados['ded'][idx], dados['dem'][idx], dados['des'][idx], ca[idx].value, pa[idx].value, dados['vel'][idx],
        dist[idx].value, dados['mR'][idx], dados['mK'][idx], dados['long'][idx]), fontsize=21, fontproperties='FreeMono')
    plt.savefig('{}_{}.png'.format(obj, datas[idx].isot),dpi=100)
    print 'Gerado: {}_{}.png'.format(obj, datas[idx].isot)
    plt.clf()

###################### rodando o programa ######################

vals = np.arange(len(dados['dia']))
map(geramapa, vals)

os.system('notify-send "Terminou de gerar os mapas" --icon=dialog-information')

