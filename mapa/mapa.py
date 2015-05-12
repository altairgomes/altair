from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation

####################################### definindo funcoes ##################################################################

def coord_pack(coord):
    if type(coord) == SkyCoord:
        if not coord.isscalar:
            return coord
        return SkyCoord([coord.ra], [coord.dec], frame=coord.frame)
    a, b = [], []
    for i in coord:
        a.append(i.ra)
        b.append(i.dec)
    coord = SkyCoord(a,b, frame='fk5')
    return coord

def calcfaixa(vel, data, star, dist, ca, pa, tamanho, step, erro=None, ring=None, atm=None):
    g = np.arange(int(-8000/(np.absolute(vel.value))), int(8000/(np.absolute(vel.value))), step)
    latlon = {'clat':{'lon':[], 'lat':[], 'lab': []}, 'lats': {'lon':[], 'lat':[], 'lon2':[], 'lat2':[]}}
    if not erro == None:
        latlon['erro'] = {'lon': [], 'lat': [], 'lon2':[], 'lat2':[]}
        err = erro*u.mas
        errd = (dist.to(u.km)*err.to(u.rad)).value*u.km
    if not ring == None:
        latlon['ring'] = {'lon': [], 'lat': [], 'lon2':[], 'lat2':[]}
    if not atm == None:
        latlon['atm'] = {'lon': [], 'lat': [], 'lon2':[], 'lat2':[]}
    paplus = ((pa > 90*u.deg) and pa - 180*u.deg) or pa
    for delt in g:
        deltatime = delt*u.s
        datas1 = data + TimeDelta(deltatime)
        datas1.delta_ut1_utc = 0
        lon = star.ra - datas1.sidereal_time('mean', 'greenwich')
        m = Basemap(projection='ortho',lat_0=star.dec.value,lon_0=lon.value,resolution=None)
        a, b = m(lon.value, star.dec.value)
        a = a*u.m
        b = b*u.m
        dista = (dist.to(u.km)*ca.to(u.rad)).value*u.km
        ax = a + dista*np.sin(pa) + (deltatime*vel)*np.cos(paplus)
        by = b + dista*np.cos(pa) - (deltatime*vel)*np.sin(paplus)
        ax2 = ax - tamanho/2*np.sin(paplus)
        by2 = by - tamanho/2*np.cos(paplus)
        ax3 = ax + tamanho/2*np.sin(paplus)
        by3 = by + tamanho/2*np.cos(paplus)
        clon1, clat1 = m(ax.value, by.value, inverse=True)
        if clon1 < 1e+30:
            latlon['clat']['lon'].append(clon1)
            latlon['clat']['lat'].append(clat1)
            latlon['clat']['lab'].append(datas1.iso.split()[1][0:8])
        lon1, lat1 = m(ax2.value, by2.value, inverse=True)
        if lon1 < 1e+30:
            latlon['lats']['lon'].append(lon1) 
            latlon['lats']['lat'].append(lat1)
        lon2, lat2 = m(ax3.value, by3.value, inverse=True)
        if lon2 < 1e+30:
            latlon['lats']['lon2'].append(lon2) 
            latlon['lats']['lat2'].append(lat2)
        if not erro == None:
            ax2 = ax - errd*np.sin(paplus)
            by2 = by - errd*np.cos(paplus)
            ax3 = ax + errd*np.sin(paplus)
            by3 = by + errd*np.cos(paplus)
            lon1, lat1 = m(ax2.value, by2.value, inverse=True)
            if lon1 < 1e+30:
                latlon['erro']['lon'].append(lon1) 
                latlon['erro']['lat'].append(lat1)
            lon2, lat2 = m(ax3.value, by3.value, inverse=True)
            if lon2 < 1e+30:
                latlon['erro']['lon2'].append(lon2) 
                latlon['erro']['lat2'].append(lat2)
        if not ring == None:
            rng = ring*u.km
            ax2 = ax - rng*np.sin(paplus)
            by2 = by - rng*np.cos(paplus)
            ax3 = ax + rng*np.sin(paplus)
            by3 = by + rng*np.cos(paplus)
            lon1, lat1 = m(ax2.value, by2.value, inverse=True)
            if lon1 < 1e+30:
                latlon['ring']['lon'].append(lon1) 
                latlon['ring']['lat'].append(lat1)
            lon2, lat2 = m(ax3.value, by3.value, inverse=True)
            if lon2 < 1e+30:
                latlon['ring']['lon2'].append(lon2) 
                latlon['ring']['lat2'].append(lat2)
        if not atm == None:
            atmo = atm*u.km
            ax2 = ax - atmo*np.sin(paplus)
            by2 = by - atmo*np.cos(paplus)
            ax3 = ax + atmo*np.sin(paplus)
            by3 = by + atmo*np.cos(paplus)
            lon1, lat1 = m(ax2.value, by2.value, inverse=True)
            if lon1 < 1e+30:
                latlon['atm']['lon'].append(lon1) 
                latlon['atm']['lat'].append(lat1)
            lon2, lat2 = m(ax3.value, by3.value, inverse=True)
            if lon2 < 1e+30:
                latlon['atm']['lon2'].append(lon2) 
                latlon['atm']['lat2'].append(lat2)
    return latlon

def geramapa(star, data, centermap=None, resolution='l', lats=None, erro=None, ring=None, atm=None):
    lons1, lats1, lons2, lats2, clon, clat, lab = calcfaixa(idx)
    lon = star.ra - data.sidereal_time('mean', 'greenwich')
    center_map = EarthLocation(lon.value, star.dec.value)
    if not centermap = None:
        center_map = EarthLocation(centermap[0],centermap[1])
        lati = centermap[1]
#    fig = plt.figure(figsize=(mapsize[0].to(u.imperial.inch).value, mapsize[1].to(u.imperial.inch).value))
    m = Basemap(projection='ortho',lat_0=center_map.latitude.value,lon_0=center_map.longitude.value,resolution=resolution)
#    m = Basemap(projection='ortho',lat_0=center_map.latitude.value,lon_0=center_map.longitude.value,resolution=resolution,llcrnrx=-2000000.,llcrnry=-1500000.,urcrnrx=2000000.,urcrnry=1500000.)
#    kx = fig.add_axes([-0.003,-0.001,1.006,1.002])
#    kx.set_rasterization_zorder(1)
    m.nightshade(data.datetime, alpha=0.2, zorder=0.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.5)
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))
    m.drawmapboundary()
    ptcolor = 'red'
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
    xs, ys = m(lons1, lats1)
    xs = [i for i in xs if i < 1e+30]
    ys = [i for i in ys if i < 1e+30]
    m.plot(xs, ys, 'b')
    xt, yt = m(lons2, lats2)
    xt = [i for i in xt if i < 1e+30]
    yt = [i for i in yt if i < 1e+30]
    m.plot(xt, yt, 'b')
    xc, yc = m(clon, clat)
    labe = [lab[i] for i in np.arange(len(lab)) if xc[i] < 1e+30]
    xc = [i for i in xc if i < 1e+30]
    yc = [i for i in yc if i < 1e+30]
    m.plot(xc, yc, 'or')

#    m.plot(ax,by, 'o', color=ptcolor, markersize=int(mapsize[0].value*20/46))
#    m.plot(ax2.to(u.m),by2.to(u.m), 'o', color=ptcolor, markersize=int(mapsize[0].value*12/46))
#    m.plot(ax3.to(u.m), by3.to(u.m), color='red')
#    m.plot(ax4.to(u.m), by4.to(u.m), color='red')
#    m.quiver(a-0*u.m,b-600000*u.m, 20, 0, width=0.005)

#    ax2 = a + dista*np.sin(pa) + [(i - datas).sec for i in temposplot]*u.s*vel*np.cos(paplus)
#    by2 = b + dista*np.cos(pa) - [(i - datas).sec for i in temposplot]*u.s*vel*np.sin(paplus)
#
#    labels = [i.iso.split()[1][0:8] for i in temposplot]
#    m.plot(ax2, by2, 'ro')
    
#    for label, axpt, bypt in zip(labe, xc, yc):
#        kx.text(axpt + 0, bypt + 350000, label, rotation=60, weight='bold')


#    if os.path.isfile(sitearq) == True:
#        xpt,ypt = m(sites['lon'],sites['lat'])
#        m.plot(xpt,ypt,'bo')
#        offset = [[xpt[0] + 100000,xpt[0] + 100000,xpt[0] + 100000,xpt[3] + 400000,xpt[3] +400000,xpt[3] +400000,xpt[3] +400000,xpt[3] +400000,xpt[3] +400000],[10000,-30000,-60000,70000,30000,-20000,-70000,-30000,-70000]]
#        for i in np.arange(len(xpt)):
#            ax.text(offset[0][i],ypt[i]+offset[1][i],sites['nome'][i], weight='bold')

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
    plt.savefig('{}_{}.png'.format(obj, datas.isot), format='png', dpi=100)
    print 'Gerado: {}_{}.png'.format(obj, datas.isot)
    plt.clf()
######################################### lendo arquivo de dados da ocultacao e de observatorios ##############################
class Map(object):
    
    def readfile(self, arquivo):
        dados = np.loadtxt(arquivo, skiprows=41, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19, 20, 21, 22, 25, 26, 28, 29), \
            dtype={'names': ('dia', 'mes', 'ano', 'hor', 'min', 'sec', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'ca', 'pa', 'vel', 'delta', 'mR', 'mK', 'long', 'ora', 'ode'), 
            'formats': ('S30', 'S30', 'S30','S30', 'S30', 'S30','S20', 'S20', 'S20','S20', 'S20', 'S20', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}, ndmin=1)

################## lendo coordenadas #################

        coor = dados['afh']
        for i in ['afm', 'afs', 'ded', 'dem', 'des']:
            coor = np.core.defchararray.add(coor, ' ')
            coor = np.core.defchararray.add(coor, dados[i])
        self.stars = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))

################### lendo tempo ########################

        tim=dados['ano']
        len_iso = ['-', '-', ' ', ':',':']
        arr = ['mes', 'dia', 'hor', 'min', 'sec']
        for i in np.arange(len(arr)):
            tim = np.core.defchararray.add(tim, len_iso[i]) 
            tim = np.core.defchararray.add(tim, dados[arr[i]])
        tim = np.char.array(tim) + '000'
        self.datas = Time(tim, format='iso', scale='utc')

############### definindo parametros #############
        self.ca = dados['ca']*u.arcsec
        self.pa = dados['pa']*u.deg
        self.vel = dados['vel']*(u.km/u.s)
        self.dist = dados['delta']*u.AU
        self.ob_off_ra = dados['ora']*u.mas
        self.ob_off_de = dados['ode']*u.mas
        self.st_off_ra = self.st_off_de = 0.0*u.mas
        self.magR = dados['mR']
        self.magK = dados['mK']
        self.longi = dados['long']

    def infile(self, entrada='mapa_in.dat'):
        f = open(entrada, 'r')
        in_data = f.readlines()
        option = in_data[0].rsplit()[0]
        self.obj = in_data[18].rsplit()[0]
        self.tamanho = int(in_data[19].rsplit()[0])*u.km
        self.sitearq = in_data[21].rsplit()[0]
        self.resolution = in_data[22].rsplit()[0]
        self.mapsize = map(float, in_data[23].rsplit()[0:2])*u.cm
        self.mapstyle = in_data[24].rsplit()[0]
        if option == '1':
            arquivo = in_data[2].rsplit()[0]
            readfile(arquivo)
            st_off_ra = st_off_de = 0.0*u.mas
        elif option == '2':
            self.stars = SkyCoord([in_data[4].rsplit('#')[0]], frame='icrs', unit=(u.hourangle, u.degree))
            datas = Time([in_data[5].rsplit('#')[0]], format='iso', scale='utc')
            ca = [float(in_data[6].rsplit('#')[0])]*u.arcsec
            self.pa = [float(in_data[7].rsplit('#')[0])]*u.deg
            self.dist = [float(in_data[8].rsplit('#')[0])]*u.AU
            self.vel = [float(in_data[9].rsplit('#')[0])]*(u.km/u.s)
            self.ob_off_ra = [float(in_data[10].rsplit('#')[0])]*u.mas
            self.ob_off_de = [float(in_data[11].rsplit('#')[0])]*u.mas
            self.st_off_ra = [float(in_data[12].rsplit('#')[0])]*u.mas
            self.st_off_de = [float(in_data[13].rsplit('#')[0])]*u.mas
            off_ra = self.ob_off_ra - self.st_off_ra
            off_de = self.ob_off_de - self.st_off_de
            self.magR = [float(in_data[14].rsplit()[0])]
            self.magK = [float(in_data[15].rsplit()[0])]
            self.longi = [float(in_data[16].rsplit()[0])]
            dca = off_ra*np.sin(pa) + off_de*np.cos(pa)
            dt = int(((off_ra*np.cos(self.pa) - off_de*np.sin(self.pa)).to(u.rad)*self.dist.to(u.km)/self.vel).value)*u.s
            self.ca = ca + dca
            self.datas = datas + dt
        f.close()
#        paplus = ((pa > 90*u.deg) and pa - 180*u.deg) or pa

    def calcfaixa(self, step=1, erro=None, ring=None, atm=None):
        self.latlon = {}
        for i in np.arange(len(self.stars)):
            self.latlon[self.datas[i].iso] = calcfaixa(self.vel[i], self.datas[i], self.stars[i], self.dist[i], self.ca[i], self.pa[i], self.tamanho, step=step, erro=erro, ring=ring, atm=atm)

################################### definido funcao que imprime o mapa #############################################

def geramapa(idx):
    lon = stars[idx].ra - datas[idx].sidereal_time('mean', 'greenwich')

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
    CS=m.nightshade(datas[idx].datetime, alpha=0.2)
    a, b =m(lon.value, stars[idx].dec.value)
    a = a*u.m
    b = b*u.m
    dista = (dist[idx].to(u.km)*ca[idx].to(u.rad)).value*u.km
    disterr = (dist[idx].to(u.km)*erro.to(u.rad)).value*u.km
    vec = np.arange(0,7000,(np.absolute(vel)*(60*u.s)).value)*u.km + np.absolute(vel)*(60*u.s)
    vec = np.concatenate((vec.value,-vec.value), axis=0)*u.km
    ax = a + dista*np.sin(pa[idx])
    ax2 = ax + vec*np.cos(pa[idx])
    ax3 = ax2 - tamanho/2*np.sin(pa[idx])
    ax4 = ax2 + tamanho/2*np.sin(pa[idx])
    ax5 = a + (dista-disterr)*np.sin(pa[idx]) + vec*np.cos(pa[idx])
    ax6 = a + (dista+disterr)*np.sin(pa[idx]) + vec*np.cos(pa[idx])
    by = b + dista*np.cos(pa[idx])
    by2 = by - vec*np.sin(pa[idx])
    by3 = by2 - tamanho/2*np.cos(pa[idx])
    by4 = by2 + tamanho/2*np.cos(pa[idx])
    by5 = b + (dista-disterr)*np.cos(pa[idx]) - vec*np.sin(pa[idx])
    by6 = b + (dista+disterr)*np.cos(pa[idx]) - vec*np.sin(pa[idx])
    m.plot(ax,by, 'o', color=ptcolor, markersize=mapsize[0].value*20/46)
    m.plot(ax2.to(u.m),by2.to(u.m), 'o', color=ptcolor, markersize=mapsize[0].value*8/46)
    m.plot(ax3.to(u.m), by3.to(u.m), color=lncolor)
    m.plot(ax4.to(u.m), by4.to(u.m), color=lncolor)
    m.plot(ax5.to(u.m), by5.to(u.m), '--', color=dscolor, label='+-{} error'.format(erro))
    m.plot(ax6.to(u.m), by6.to(u.m), '--', color=dscolor)
#    plt.legend(fontsize=mapsize[0].value*21/46)

    fig = plt.gcf()
    fig.set_size_inches(mapsize[0].to(u.imperial.inch).value, mapsize[1].to(u.imperial.inch).value)
    plt.title('Objeto       Diam   dots <>  ra_off_obj_de  ra_of_star_de\n{:10s} {:4.0f} km  60 s <> {:+6.1f} {:+6.1f}  {:+6.1f} {:+6.1f} \n'
        .format(obj, tamanho.value, ob_off_ra[idx].value, ob_off_de[idx].value, st_off_ra[idx].value, st_off_de[idx].value), fontsize=mapsize[0].value*25/46, fontproperties='FreeMono', weight='bold')
    plt.xlabel('\n year-m-d    h:m:s UT     ra__dec__J2000__candidate    C/A    P/A    vel   Delta   R*   K*  long\n\
{}  {:02d} {:02d} {:07.4f} {:+02d} {:02d} {:06.3f} {:6.3f} {:6.2f} {:6.2f}  {:5.2f} {:5.1f} {:4.1f}  {:3.0f}'
        .format(datas[idx].iso, int(stars[idx].ra.hms.h), int(stars[idx].ra.hms.m), stars[idx].ra.hms.s, int(stars[idx].dec.dms.d), np.absolute(int(stars[idx].dec.dms.m)), np.absolute(stars[idx].dec.dms.s),
        ca[idx].value, pa[idx].value, vel[idx].value, dist[idx].value, magR[idx], magK[idx], longi[idx]), fontsize=mapsize[0].value*21/46, fontproperties='FreeMono', weight='bold')
    plt.savefig('{}_{}.png'.format(obj, datas[idx].isot),dpi=100)
    print 'Gerado: {}_{}.png'.format(obj, datas[idx].isot)
    plt.clf()

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
    datas = Time([in_data[5].rsplit('#')[0]], format='iso', scale='utc')
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

datas.delta_ut1_utc = 0

if os.path.isfile(sitearq) == True:
    sites = np.loadtxt(sitearq,  dtype={'names': ('lat', 'lon', 'alt', 'nome'), 'formats': ('f8', 'f8', 'f8', 'S30')})

###################### rodando o programa ######################

#map(geramapa, vals)

os.system('notify-send "Terminou de gerar os mapas" --icon=dialog-information')

