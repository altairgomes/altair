from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation
from multiprocessing import Pool

################################### definido funcao que imprime o mapa #############################################

def geramapa(delt):
    deltatime = delt*u.s
    datas1 = datas[idx] + TimeDelta(deltatime)
    datas1.delta_ut1_utc = 0
    lon = stars[idx].ra - datas1.sidereal_time('mean', 'greenwich')

    m = Basemap(projection='ortho',lat_0=stars[idx].dec.value,lon_0=lon.value,resolution=resolution)
#    m = Basemap(projection='ortho',lat_0=stars[idx].dec.value,lon_0=lon.value,resolution=resolution, llcrnrx=-7000000,llcrnry=-7000000,urcrnrx=7000000,urcrnry=7000000)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.5)
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))
    m.drawmapboundary()
    ptcolor = 'black'
    lncolor = 'black'
    dscolor = 'black'
    if mapstyle == '2':
        m.drawmapboundary(fill_color='aqua')
        m.fillcontinents(color='coral',lake_color='aqua')
        ptcolor = 'red'
        lncolor = 'blue'
        dscolor = 'red'
    elif mapstyle == '3':
        m.shadedrelief()
        ptcolor = 'red'
        lncolor = 'blue'
        dscolor = 'red'
    elif mapstyle == '4':
        m.bluemarble()
        ptcolor = 'red'
        lncolor = 'red'
        dscolor = 'red'
    elif mapstyle == '5':
        m.etopo()
        ptcolor = 'red'
        lncolor = 'red'
        dscolor = 'red'
    if os.path.isfile(sitearq) == True:
        xpt,ypt = m(sites['lon'],sites['lat'])
        m.plot(xpt,ypt,'bo')
    CS=m.nightshade(datas1.datetime, alpha=0.2)
    a, b =m(lon.value, stars[idx].dec.value)
    a = a*u.m
    b = b*u.m
    dista = (dist[idx].to(u.km)*ca[idx].to(u.rad)).value*u.km
    disterr = (dist[idx].to(u.km)*erro.to(u.rad)).value*u.km
    ax = a + dista*np.sin(pa[idx]) + (deltatime*vel[idx])*np.cos(pa[idx])
    by = b + dista*np.cos(pa[idx]) - (deltatime*vel[idx])*np.sin(pa[idx])
    m.plot(ax,by, 'o', color=ptcolor, markersize=mapsize[0].value*20/46)
#    plt.legend(fontsize=mapsize[0].value*21/46)

    fig = plt.gcf()
    fig.set_size_inches(mapsize[0].to(u.imperial.inch).value, mapsize[1].to(u.imperial.inch).value)
    plt.title('-{} D={}- dots each 60 s <> offsets (mas): obj=({:.1f},{:.1f}), star=({:.1f},{:.1f})\n'
        .format(obj, tamanho, ob_off_ra[idx].value, ob_off_de[idx].value, st_off_ra[idx].value, st_off_de[idx].value), fontsize=mapsize[0].value*25/46, fontproperties='FreeMono', weight='bold')
    plt.xlabel('\n year-m-d    h:m:s UT     ra__dec__J2000__candidate    C/A    P/A    vel   Delta   R*   K*  long\n\
{}  {:02d} {:02d} {:07.4f} {:+02d} {:02d} {:06.3f} {:6.3f} {:6.2f} {:6.2f}  {:5.2f} {:5.1f} {:4.1f}  {:3.0f}'
        .format(datas1.iso, int(stars[idx].ra.hms.h), int(stars[idx].ra.hms.m), stars[idx].ra.hms.s, int(stars[idx].dec.dms.d), np.absolute(int(stars[idx].dec.dms.m)), np.absolute(stars[idx].dec.dms.s),
        ca[idx].value, pa[idx].value, vel[idx].value, dist[idx].value, magR[idx], magK[idx], longi[idx]), fontsize=mapsize[0].value*21/46, fontproperties='FreeMono', weight='bold')
    plt.savefig('{}_{:05d}.png'.format(obj, np.where(g==delt)[0][0] + 1),dpi=100)
    print 'Gerado: {}_{:05d}.png'.format(obj, np.where(g==delt)[0][0] + 1)
    plt.clf()

###########################################################################################################################
###########################################################################################################################
######################################## lendo arquivo de entrada #########################################################

f = open('mapa_gif_in.dat', 'r')

in_data = f.readlines()
obj = in_data[14].rsplit()[0]
tamanho = in_data[15].rsplit()[0]*u.km
erro = in_data[16].rsplit()[0]*u.mas
sitearq = in_data[17].rsplit()[0]
resolution = in_data[18].rsplit()[0]
mapsize = map(float, in_data[19].rsplit()[0:2])*u.cm
mapstyle = in_data[20].rsplit()[0]
passofilme = int(in_data[21].rsplit()[0])
fpsvideo = int(in_data[22].rsplit()[0])
stars = SkyCoord([in_data[0].rsplit('#')[0]], frame='icrs', unit=(u.hourangle, u.degree))
datas = Time([in_data[1].rsplit('#')[0]], format='iso', scale='utc')
ca = [float(in_data[2].rsplit('#')[0])]*u.arcsec
pa = [float(in_data[3].rsplit('#')[0])]*u.deg
dist = [float(in_data[4].rsplit('#')[0])]*u.AU
vel = [float(in_data[5].rsplit('#')[0])]*(u.km/u.s)
ob_off_ra = [float(in_data[6].rsplit('#')[0])]*u.mas
ob_off_de = [float(in_data[7].rsplit('#')[0])]*u.mas
st_off_ra = [float(in_data[8].rsplit('#')[0])]*u.mas
st_off_de = [float(in_data[9].rsplit('#')[0])]*u.mas
off_ra = ob_off_ra - st_off_ra
off_de = ob_off_de - st_off_de
magR = [float(in_data[10].rsplit()[0])]
magK = [float(in_data[11].rsplit()[0])]
longi = [float(in_data[12].rsplit()[0])]
dca = off_ra*np.sin(pa) + off_de*np.cos(pa)
dt = int(((off_ra*np.cos(pa) - off_de*np.sin(pa)).to(u.rad)*dist.to(u.km)/vel).value)*u.s
ca = ca + dca
datas = datas + dt
idx = 0

f.close()


if os.path.isfile(sitearq) == True:
    sites = np.loadtxt(sitearq,  dtype={'names': ('lat', 'lon', 'alt', 'nome'), 'formats': ('f8', 'f8', 'f8', 'S30')})

###################### rodando o programa ######################

g = np.arange(int(-8000/(np.absolute(vel[0].value))), int(8000/(np.absolute(vel[0].value))) , passofilme)
pool = Pool(processes=10)
pool.map(geramapa, g)

os.system('ffmpeg -f image2 -r {} -sameq -i "{}_%05d.png" {}_{}.mp4'.format(fpsvideo, obj, obj, datas[0].isot))

os.system('notify-send "Terminou de gerar os mapas" --icon=dialog-information')

