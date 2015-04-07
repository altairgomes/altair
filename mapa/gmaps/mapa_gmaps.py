from mpl_toolkits.basemap import Basemap
import numpy as np
import os
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation
import pygmaps

######################################### lendo arquivo de dados da ocultacao e de observatorios ##############################

def lerdados():
    dados = np.loadtxt(arquivo, skiprows=41, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19, 20, 21, 22, 25, 26, 28, 29), \
        dtype={'names': ('dia', 'mes', 'ano', 'hor', 'min', 'sec', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'ca', 'pa', 'vel', 'delta', 'mR', 'mK', 'long', 'ora', 'ode'), 
        'formats': ('S30', 'S30', 'S30','S30', 'S30', 'f8','i4', 'i4', 'f8','i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}, ndmin=1)

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

############### definindo parametros #############
    ca = dados['ca']*u.arcsec
    pa = dados['pa']*u.deg
    vel = dados['vel']*(u.km/u.s)
    dist = dados['delta']*u.AU
    off_ra = dados['ora']*u.mas
    off_de = dados['ode']*u.mas
    return stars, datas, ca, pa, vel, dist, off_ra, off_de, dados['mR'], dados['mK'], dados['long']

################################### definido funcao que imprime o mapa #############################################
def calcfaixa(idx):
    g = np.arange(int(-8000/(np.absolute(vel[0].value))), int(8000/(np.absolute(vel[0].value))), 1)
    lons1, lats1, lons2, lats2, clon, clat, lab = [], [], [], [], [], [], []
    for delt in g:
        deltatime = delt*u.s
        datas1 = datas + TimeDelta(deltatime)
        datas1.delta_ut1_utc = 0
        lon = stars[idx].ra - datas1.sidereal_time('mean', 'greenwich')
        m = Basemap(projection='ortho',lat_0=stars[idx].dec.value,lon_0=lon.value,resolution=None)
        a, b =m(lon.value, stars[idx].dec.value)
        a = a*u.m
        b = b*u.m
        dista = (dist[idx].to(u.km)*ca[idx].to(u.rad)).value*u.km
        ax = a + dista*np.sin(pa[idx]) + (deltatime*vel[idx])*np.cos(pa[idx])
        by = b + dista*np.cos(pa[idx]) - (deltatime*vel[idx])*np.sin(pa[idx])
        ax2 = ax - tamanho/2*np.sin(pa[idx])
        by2 = by - tamanho/2*np.cos(pa[idx])
        ax3 = ax + tamanho/2*np.sin(pa[idx])
        by3 = by + tamanho/2*np.cos(pa[idx])
        clon1, clat1 = m(ax.value, by.value, inverse=True)
#        if deltatime.value%5 == 0:
        clon.append(clon1)
        clat.append(clat1)
        lab.append(datas1.iso.split()[1][0:8])
        lon1, lat1 = m(ax2.value, by2.value, inverse=True)
        lons1.append(lon1) 
        lats1.append(lat1)
        lon2, lat2 = m(ax3.value, by3.value, inverse=True)
        lons2.append(lon2) 
        lats2.append(lat2)
    return lons1, lats1, lons2, lats2, clon, clat, lab

def geramapa(idx):

    lons1, lats1, lons2, lats2, clon, clat, lab = calcfaixa(idx)
    
#    center_map = EarthLocation('-77 02 28.3','38 49 19.1')
    lon = stars[idx].ra - datas.sidereal_time('mean', 'greenwich')
    lon.wrap_at('180d', inplace=True)
    print stars[idx].dec.value, lon.value
    mymap = pygmaps.maps(stars[idx].dec.value, lon.value, 5)
    lons1 = [i for i in lons1 if i < 1e+20]
    lats1 = [i for i in lats1 if i < 1e+20]
    lons2 = [i for i in lons2 if i < 1e+20]
    lats2 = [i for i in lats2 if i < 1e+20]
    path1 = []
    path2 = []
    for i in np.arange(len(lons1)):
        path1.append((lats1[i], lons1[i]))
    for i in np.arange(len(lons2)):
        path2.append((lats2[i], lons2[i]))
    mymap.addpath(path1,"#0000FF")
    mymap.addpath(path2,"#0000FF")
    mymap.draw('./mymap.html')
#    fig = plt.figure(figsize=(mapsize[0].to(u.imperial.inch).value, mapsize[1].to(u.imperial.inch).value))
#    m = Basemap(projection='ortho',lat_0=center_map.latitude.value,lon_0=center_map.longitude.value,resolution=resolution,llcrnrx=-1000000.,llcrnry=-750000.,urcrnrx=1000000.,urcrnry=750000.)
##    m = Basemap(projection='ortho',lat_0=stars[idx].dec.value,lon_0=lon.value,resolution=resolution)
#    ax = fig.add_axes([-0.003,-0.001,1.006,1.002])
#    m.drawcoastlines(linewidth=0.5)
#    m.drawcountries(linewidth=0.5)
#    m.drawstates(linewidth=0.5)
#    m.drawmeridians(np.arange(0,360,30))
#    m.drawparallels(np.arange(-90,90,30))
#    m.drawmapboundary()
#    ptcolor = 'red'
#    lncolor = 'black'
#    dscolor = 'black'
#    if mapstyle == '2':
#        m.drawmapboundary(fill_color='aqua')
#        m.fillcontinents(color='coral',lake_color='aqua')
#        ptcolor = 'red'
#        lncolor = 'blue'
#        dscolor = 'red'
#    elif mapstyle == '3':
#        m.shadedrelief()
#        ptcolor = 'red'
#        lncolor = 'blue'
#        dscolor = 'red'
#    elif mapstyle == '4':
#        m.bluemarble()
#        ptcolor = 'red'
#        lncolor = 'red'
#        dscolor = 'red'
#    elif mapstyle == '5':
#        m.etopo()
#        ptcolor = 'red'
#        lncolor = 'red'
#        dscolor = 'red'
#    CS=m.nightshade(datas.datetime, alpha=0.2)
#    a, b =m(center_map.longitude.value, center_map.latitude.value)
#    a = a*u.m
#    b = b*u.m
#    dista = (dist[idx].to(u.km)*ca[idx].to(u.rad)).value*u.km
#    disterr = (dist[idx].to(u.km)*erro.to(u.rad)).value*u.km
#    vec = np.arange(0,7000,(np.absolute(vel)*(30*u.s)).value)*u.km + np.absolute(vel)*(30*u.s)
#    vec = np.concatenate((vec.value,-vec.value), axis=0)*u.km
#    ax = a + dista*np.sin(pa[idx])
#    ax2 = ax + vec*np.cos(pa[idx])
#    ax3 = ax2 - tamanho/2*np.sin(pa[idx])
#    ax4 = ax2 + tamanho/2*np.sin(pa[idx])
#    ax5 = a + (dista-disterr)*np.sin(pa[idx]) + vec*np.cos(pa[idx])
#    ax6 = a + (dista+disterr)*np.sin(pa[idx]) + vec*np.cos(pa[idx])
#    by = b + dista*np.cos(pa[idx])
#    by2 = by - vec*np.sin(pa[idx])
#    by3 = by2 - tamanho/2*np.cos(pa[idx])
#    by4 = by2 + tamanho/2*np.cos(pa[idx])
#    by5 = b + (dista-disterr)*np.cos(pa[idx]) - vec*np.sin(pa[idx])
#    by6 = b + (dista+disterr)*np.cos(pa[idx]) - vec*np.sin(pa[idx])
#    xs, ys = m(lons1, lats1)
#    xs = [i for i in xs if i < 1e+30]
#    ys = [i for i in ys if i < 1e+30]
#    m.plot(xs, ys, 'b')
#    xt, yt = m(lons2, lats2)
#    xt = [i for i in xt if i < 1e+30]
#    yt = [i for i in yt if i < 1e+30]
#    m.plot(xt, yt, 'b')
#    xc, yc = m(clon, clat)
#    xc = [i for i in xc if i < 1e+30]
#    yc = [i for i in yc if i < 1e+30]
#    m.plot(xc, yc, 'or')

#    m.plot(ax,by, 'o', color=ptcolor, markersize=int(mapsize[0].value*20/46))
#    m.plot(ax2.to(u.m),by2.to(u.m), 'o', color=ptcolor, markersize=int(mapsize[0].value*12/46))
#    m.plot(ax3.to(u.m), by3.to(u.m), color=lncolor)
#    m.plot(ax4.to(u.m), by4.to(u.m), color=lncolor)
#    m.quiver(a-0*u.m,b-550000*u.m, 20, 0, width=0.005)

#    ax2 = a + dista*np.sin(pa) + [(i - datas).sec for i in temposplot]*u.s*vel*np.cos(paplus)
#    by2 = b + dista*np.cos(pa) - [(i - datas).sec for i in temposplot]*u.s*vel*np.sin(paplus)

#    labels = [i.iso.split()[1][0:8] for i in temposplot]
#    m.plot(ax2, by2, 'ro')
    
#    for label, axpt, bypt in zip(lab, xc, yc):
#        plt.text(axpt + 60000, bypt + 150000, label, rotation=60, weight='bold')

#    offset = [[50000,50000,50000,50000,50000,60000,50000,50000,50000],[10000,5000,2500,25000,20000,0,-20000,5000,-5000]]
#
#    if os.path.isfile(sitearq) == True:
#        xpt,ypt = m(sites['lon'],sites['lat'])
#        m.plot(xpt,ypt,'bo')
#        for i in np.arange(len(xpt)):
#            plt.text(xpt[i]+offset[0][i],ypt[i]+offset[1][i],sites['nome'][i], weight='bold')

#    m.plot(ax5.to(u.m), by5.to(u.m), '--', color=dscolor, label='+-{} error'.format(erro))
#    m.plot(ax6.to(u.m), by6.to(u.m), '--', color=dscolor)
#    plt.legend(fontsize=mapsize[0].value*21/46)

#    fig = plt.gcf()
#    fig.set_size_inches(mapsize[0].to(u.imperial.inch).value, mapsize[1].to(u.imperial.inch).value)
#    plt.title('Objeto       Diam   dots <>  ra_off_obj_de  ra_of_star_de\n{:10s} {:4.0f} km  60 s <> {:+6.1f} {:+6.1f}  {:+6.1f} {:+6.1f} \n'
#        .format(obj, tamanho.value, ob_off_ra[idx].value, ob_off_de[idx].value, st_off_ra[idx].value, st_off_de[idx].value), fontsize=mapsize[0].value*25/46, fontproperties='FreeMono', weight='bold')
#    plt.xlabel('\n year-m-d    h:m:s UT     ra__dec__J2000__candidate    C/A    P/A    vel   Delta   R*   K*  long\n\
#{}  {:02d} {:02d} {:07.4f} {:+02d} {:02d} {:06.3f} {:6.3f} {:6.2f} {:6.2f}  {:5.2f} {:5.1f} {:4.1f}  {:3.0f}'
#        .format(datas[idx].iso, int(stars[idx].ra.hms.h), int(stars[idx].ra.hms.m), stars[idx].ra.hms.s, int(stars[idx].dec.dms.d), np.absolute(int(stars[idx].dec.dms.m)), np.absolute(stars[idx].dec.dms.s),
#        ca[idx].value, pa[idx].value, vel[idx].value, dist[idx].value, magR[idx], magK[idx], longi[idx]), fontsize=mapsize[0].value*21/46, fontproperties='FreeMono', weight='bold')
#    plt.savefig('{}_{}.png'.format(obj, datas.isot),dpi=100)
#    print 'Gerado: {}_{}.png'.format(obj, datas.isot)
#    plt.clf()

###########################################################################################################################
###########################################################################################################################
######################################## lendo arquivo de entrada #########################################################

f = open('mapa_in.dat', 'r')

in_data = f.readlines()
option = in_data[0].rsplit()[0]
obj = in_data[18].rsplit()[0]
tamanho = int(in_data[19].rsplit()[0])*u.km
erro = float(in_data[20].rsplit()[0])*u.mas
sitearq = in_data[21].rsplit()[0]
resolution = in_data[22].rsplit()[0]
mapsize = map(float, in_data[23].rsplit()[0:2])*u.cm
mapstyle = in_data[24].rsplit()[0]
if option == '1':
    arquivo = in_data[2].rsplit()[0]
    stars, datas, ca, pa, vel, dist, ob_off_ra, ob_off_de, magR, magK, longi = lerdados()
    vals = np.arange(len(stars))
    st_off_ra = st_off_de = 0.0*u.mas
elif option == '2':
    stars = SkyCoord([in_data[4].rsplit('#')[0]], frame='icrs', unit=(u.hourangle, u.degree))
    datas = Time(in_data[5].rsplit('#')[0], format='iso', scale='utc')
    ca = [float(in_data[6].rsplit('#')[0])]*u.arcsec
    pa = [float(in_data[7].rsplit('#')[0])]*u.deg
    dist = [float(in_data[8].rsplit('#')[0])]*u.AU
    vel = [float(in_data[9].rsplit('#')[0])]*(u.km/u.s)
    ob_off_ra = [float(in_data[10].rsplit('#')[0])]*u.mas
    ob_off_de = [float(in_data[11].rsplit('#')[0])]*u.mas
    st_off_ra = [float(in_data[12].rsplit('#')[0])]*u.mas
    st_off_de = [float(in_data[13].rsplit('#')[0])]*u.mas
    off_ra = ob_off_ra - st_off_ra
    off_de = ob_off_de - st_off_de
    magR = [float(in_data[14].rsplit()[0])]
    magK = [float(in_data[15].rsplit()[0])]
    longi = [float(in_data[16].rsplit()[0])]
    dca = off_ra*np.sin(pa) - off_de*np.cos(pa)
    dt = int(((-off_ra*np.cos(pa) - off_de*np.sin(pa)).to(u.rad)*dist.to(u.km)/vel).value)*u.s
    ca = ca + dca
    datas = datas + dt
    vals = [0]

f.close()

paplus = ((pa > 180*u.deg) and pa - 180*u.deg) or pa
inttime = 5*u.min

deltas = float(datas.iso.rsplit(':')[2])*u.s
tempoint = datas - deltas
temposplot = [ tempoint - inttime + n*(10*u.s) for n in np.arange(2*6*(inttime.value)) ]

datas.delta_ut1_utc = 0

if os.path.isfile(sitearq) == True:
    sites = np.loadtxt(sitearq,  dtype={'names': ('lat', 'lon', 'alt', 'nome'), 'formats': ('f8', 'f8', 'f8', 'S30')})

###################### rodando o programa ######################

map(geramapa, vals)

os.system('notify-send "Terminou de gerar os mapas" --icon=dialog-information')

