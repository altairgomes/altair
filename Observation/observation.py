# coding=UTF-8

import numpy as np
import itertools
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, Latitude, Longitude, FK5, EarthLocation, Angle

######################################################################

class Observation(object):
    """
    """

    def __init__(self):
        fuso = 0
        latitude = 0
        longitude = 0
        altitude = 0
        self.set_fuse(fuso)
        self.set_site(longitude, latitude, altitude)
        
    def read(self, datafile):
        """
        """
        dados = np.loadtxt(datafile, dtype={'names': ('objeto', 'data', 'hora', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'mag'), 
            'formats': ('S30', 'S30', 'S30', 'i4', 'i4', 'f8', 'i4', 'i4', 'f8', 'S30')}, ndmin=1)

        ### leitura data da ocultacao ###
        tempo = np.core.defchararray.add(dados['data'], [' '])
        tempo = np.core.defchararray.add(tempo, dados['hora'])
        time = Time(tempo, format='iso', scale='utc')
        
        ### leitura das coordenadas dos alvos ###
        alfa = np.core.defchararray.add(np.array(map(str, dados['afh'])), [' '])
        alfa = np.core.defchararray.add(alfa, np.array(map(str, dados['afm'])))
        alfa = np.core.defchararray.add(alfa, [' '])
        alfa = np.core.defchararray.add(alfa, np.array(map(str, dados['afs'])))
        delta = np.core.defchararray.add(np.array(map(str, dados['ded'])), [' '])
        delta = np.core.defchararray.add(delta, np.array(map(str, dados['dem'])))
        delta = np.core.defchararray.add(delta, [' '])
        delta = np.core.defchararray.add(delta, np.array(map(str, dados['afs'])))
        coor = np.core.defchararray.add(alfa, [' '])
        coor = np.core.defchararray.add(coor, delta)
        coordenadas = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        
        ### mudando o nome ###
        nome = np.core.defchararray.add(np.array(map(str, dados['objeto'])), [' '])
        nome = np.core.defchararray.add(nome, np.array(map(str, time.iso)))
        
        data = np.array([], dtype={'names': ('object', 'coord', 'time_occ', 'comment'),'formats': ('S50', object, object, 'S30')})
        for i in np.arange(len(dados)):
            a = [(nome[i], coordenadas[i], time[i], dados['mag'][i])]
            b = np.asarray(a, dtype={'names': ('object', 'coord', 'time_occ', 'comment'),'formats': ('S50', object, object, 'S30')})
            data = np.append(data,b)
        return data
        
    def set_site(self,longitude, latitude, height=0.0*u.m):
        """
        """
        self.site = EarthLocation(longitude, latitude, height)

    def set_fuse(self, fuse):
        """
        """
        self.utcoffset = TimeDelta(fuse*3600, format='sec', scale='tai')

    def close_obj(self, coord, size):
        """
        """
        if len(coord) > 1:
            a, b = [], []
            for i in coord:
                a.append(i.ra)
                b.append(i.dec)
            coord = SkyCoord(a,b, frame='fk5')
        if not type(size) == u.quantity.Quantity:
            size = size*u.arcmin
        lista_campos = []
        for idx, value in enumerate(coord):
            sep = coord[idx].separation(coord[idx+1:])
            close = [i + 1 + idx for i in np.arange(len(sep)) if sep[i] < size]
            combs = []
            for i in np.arange(len(close),0, -1):
                els = [list(x) for x in itertools.combinations(close, i)]
                combs.append(els)
            for i in combs:
                for l in i:
                    a = np.concatenate(([idx], l))
                    mean_coord, dist_center = self.midpoint_coord(coord[a])
                    if all(k <= size/2 for k in dist_center):
                        d = False
                        for j in lista_campos:
                            e = all(k in j for k in a)
                            d = bool(d + e)
                        if d == False:
                            lista_campos.append(a)
            d = False
            for i in lista_campos:
                e = idx in i
                d = bool(d + e)
            if d == False:
                lista_campos.append(np.array([idx]))
        campos = []
        for i in lista_campos:
            campos.append(coord[i])
        return campos, lista_campos
    
    def midpoint_coord(self, coord):
        """
        """
        x = np.mean([np.cos(i.dec)*np.cos(i.ra) for i in coord])
        y = np.mean([np.cos(i.dec)*np.sin(i.ra) for i in coord])
        z = np.mean([np.sin(i.dec) for i in coord])
        delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))*u.rad
        alfa = np.arctan(y/x)*u.rad
        mean_coord = SkyCoord(alfa, delta, frame='fk5')
        dist_center = coord.separation(mean_coord)
        return mean_coord, dist_center
        
    def precess(self, coord, time):
        """
        """
        if len(coord) > 1:
            a, b = [], []
            for i in coord:
                a.append(i.ra)
                b.append(i.dec)
            coord = SkyCoord(a,b, frame='fk5')
        fk5_data = FK5(equinox=time)
        coord_prec = coord.transform_to(fk5_data)
        return coord_prec
        
    def sky_time(self, coord, tempo, limalt=0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        if not type(limalt) == u.quantity.Quantity:
            limalt = limalt*u.deg
        if len(coord) > 1:
            a, b = [], []
            for i in coord:
                a.append(i.ra.deg)
                b.append(i.dec.deg)
            coord = SkyCoord(a,b, frame='fk5', unit=(u.deg,u.deg))
        tempo.delta_ut1_utc = 0
        tempo.location = site
        dif_h_sid = coord.ra - tempo.sidereal_time('mean')
        dif_h_sid.wrap_at('180d', inplace=True)
        dif_h_sol = dif_h_sid * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
        dif = TimeDelta(dif_h_sol.hour*60*60, format='sec')
        hangle_lim = np.arccos((np.cos(90.0*u.deg-limalt) - np.sin(coord.dec)*np.sin(site.latitude)) / (np.cos(coord.dec)*np.cos(site.latitude)))
        culminacao = tempo + dif
        culminacao.delta_ut1_utc = 0
        tsg_lim = coord.ra + hangle_lim
        dtsg_lim = tsg_lim - culminacao.sidereal_time('mean')
        dtsg_lim.wrap_at(360 * u.deg, inplace=True)
        dtsg_lim_sol = dtsg_lim * (23.0 + 56.0/60.0 + 4.0916/3600.0) / 24.0
        dtsg_np = TimeDelta(dtsg_lim_sol.hour*60*60, format='sec')
        sunrise = culminacao - dtsg_np
        sunset = culminacao + dtsg_np
        return culminacao, sunrise, sunset
        
    def height_time(self, coord, instante, time_left=False, limalt=0.0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        if len(coord) > 1:
            a, b = [], []
            for i in coord:
                a.append(i.ra)
                b.append(i.dec)
            coord = SkyCoord(a,b, frame='fk5')
        instante.location = site
        instante.delta_ut1_utc = 0
        hourangle = instante.sidereal_time('mean') - coord.ra
        distzen = np.arccos(np.sin(coord.dec)*np.sin(site.latitude) + np.cos(coord.dec)*np.cos(site.latitude)*np.cos(hourangle))
        altura = 90*u.deg - distzen
        if time_left == True:
            poente = self.sky_time(coord, instante, limalt, site)[2]
            time_rest = poente - instante
            return altura, time_rest
        return altura
        
    def plan(self, coord, hora, obj=None, comment=None, culmination=None, limalt=0.0*u.arcmin, size=0.0*u.deg, site=None):
        if site == None:
            site = self.site
        if obj == None:
            obj = ['']*len(coord)
        if comment == None:
            comment = ['']*len(coord)
        comment1 = comment.copy()
        for i in np.arange(len(comment)):
            comment1[i] = '({})'.format(comment[i])
        if not type(limalt) == u.quantity.Quantity:
            limalt = limalt*u.deg
        if not type(size) == u.quantity.Quantity:
            size = size*u.arcmin
        if len(coord) > 1:
            a, b = [], []
            for i in coord:
                a.append(i.ra)
                b.append(i.dec)
            coord = SkyCoord(a,b, frame='fk5')
        midcoord = np.copy(coord)
        midobj = obj.copy()
        midcomment = comment1.copy()
        lista_campos =  [np.array([0])]*len(midobj)
        if size > 0.0*u.arcmin:
            midcoord, midobj, midcomment = [], [], []
            campos, lista_campos = self.close_obj(coord, size)
            for i in lista_campos:
                if len(i) > 1:
                    midcoord.append(self.midpoint_coord(coord[i])[0])
                    a = ''
                    for k in obj[i][:-1]:
                        a  = a + '{} + '.format(k)
                    a = a + obj[i][-1]
                    midobj.append([a])
                    midcomment.append(np.array(['']))
                else:
                    midobj.append(obj[i])
                    midcoord.append(coord[i][0])
                    midcomment.append(comment1[i])
            midculmination = self.sky_time(midcoord, hora, site=site)
        if len(midcoord) > 1:
            a, b = [], []
            for i in midcoord:
                a.append(i.ra)
                b.append(i.dec)
            midcoord = SkyCoord(a,b, frame='fk5')
        if culmination == None:
            culmination = self.sky_time(midcoord, hora, site=site)[0]
        altura, time_rest = self.height_time(midcoord, hora, limalt=limalt, time_left=True)
        obs_tot = np.array([], dtype={'names': ('coord', 'height', 'time_left', 'obj', 'comment', 'culmination', 'tam_campo'),'formats': (object, 'f16', object, 'S100', 'S30', object, list)})
        for i in np.arange(len(altura)):
            a = [(midcoord[i], altura[i].value, time_rest[i], midobj[i][0], midcomment[i][0], culmination[i], lista_campos[i])]
            b = np.asarray(a, dtype={'names': ('coord', 'height', 'time_left', 'obj', 'comment', 'culmination', 'tam_campo'),'formats': (object, 'f16', object, 'S100', 'S30', object, list)})
            obs_tot = np.append(obs_tot,b)
        obs = np.sort(obs_tot[obs_tot['height']*u.deg >= limalt], order='time_left')
        a = ''
        for i in np.arange(len(obs)):
            a = a + '{} {}\n\tRA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n\tHeight: {:.1f}\n\tCulmination: {} LT\n\tTime left to reach min height: {:02d}:{:02d}\n'\
.format(obs['obj'][i], obs['comment'][i], obs['coord'][i].ra.hms.h, obs['coord'][i].ra.hms.m, obs['coord'][i].ra.hms.s, 
obs['coord'][i].dec.dms.d, np.absolute(obs['coord'][i].dec.dms.m), np.absolute(obs['coord'][i].dec.dms.s), obs['height'][i]*u.deg, (obs['culmination'][i] + self.utcoffset).iso.split(' ')[1][0:5],
int(obs['time_left'][i].sec/3600), int((obs['time_left'][i].sec - int(obs['time_left'][i].sec/3600)*3600)/60))
            if len(obs['tam_campo'][i]) > 1:
                for k in obs['tam_campo'][i]:
                    a = a + '\t\t{} {}\n\t\t  RA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n'.format(obj[k], comment1[k], coord[k].ra.hms.h,
coord[k].ra.hms.m, coord[k].ra.hms.s, coord[k].dec.dms.d, np.absolute(coord[k].dec.dms.m), np.absolute(coord[k].dec.dms.s))
            a = a + '\n'
        return a
                    
    def create_plan(self, coord, obj, comment, horain, horafin, intinf, limalt=0.0*u.deg, size=0*u.arcmin, site=None):
        """
        """
        if site == None:
            site = self.site
        if not type(limalt) == u.quantity.Quantity:
            limalt = limalt*u.deg
        if not type(size) == u.quantity.Quantity:
            size = size*u.arcmin
        tempoin = Time(horain, format='iso', scale='utc', location=site) - self.utcoffset
        tempofin = Time(horafin, format='iso', scale='utc', location=site) + TimeDelta(0.5, format='sec')  - self.utcoffset
        intval = TimeDelta(intinf*60, format='sec')
        a, b = [], []
        for i in coord:
            a.append(i.ra)
            b.append(i.dec)
        coord = SkyCoord(a,b, frame='fk5')
        culminacao, nascer, poente = self.sky_time(coord, tempoin + (tempofin-tempoin)/2, limalt, site)
        nome = 'Plano_{}'.format(horain.split(' ')[0])
        output = open(nome, 'w')
        output.write('Observational Plan to the night: {}\n\n'.format(horain.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nMinimun height: {}\nField Size: {}\n\n'.format(site.latitude, site.longitude, limalt, size))
        instante = tempoin
        while instante <= tempofin:
            #### imprime os dados de cada objeto para cada instante ####
            text = self.plan(coord, instante, obj, comment, culminacao, site=site, limalt=limalt, size=size)
            output.write('\n---LT: {} (UT: {})----------------------------------------------------------------\n'.format((instante + self.utcoffset).iso.split(' ')[1][0:5], instante.iso.split(' ')[1][0:5]))
            output.write(text) 
            instante = instante + intval
        output.write('\n---Observability of the Targets----------------------------------------------------------------\n')
        for idx, val in enumerate(nascer):
            output.write('RA: {:02.0f} {:02.0f} {:07.4f}, DEC: {:+03.0f} {:02.0f} {:06.3f}, Sunrise: {} TL, Culminaction: {} TL, Sunset: {} TL, {:10s}\n'\
.format(coord[idx].ra.hms.h, coord[idx].ra.hms.m, coord[idx].ra.hms.s, coord[idx].dec.dms.d, np.absolute(coord[idx].dec.dms.m), np.absolute(coord[idx].dec.dms.s),
(nascer[idx] + self.utcoffset).iso.split(' ')[1][0:5], (culminacao[idx] + self.utcoffset).iso.split(' ')[1][0:5], (poente[idx] + self.utcoffset).iso.split(' ')[1][0:5], obj[idx]))
        output.write('\n')
        output.close()
        
#####################################################################

arquivo = 'alvos.txt'				#### arquivo de alvos
horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
intinf = 60					#### intervalo entre as informacoes (minutos)
fuso = -3					#### fuso horario do local
latitude = '-22 32 7.8'				#### latitude do local
longitude = '314 25 2.5'			#### longitude do local
altitude = 1864					#### altitude em metros
limalt = 30.0					#### limite de altura para mostrar (graus)
limdist = 11					#### limite de distancia para field-of-view (arcmin)

#######################################################################

observation = Observation()
dados = observation.read(arquivo)
observation.set_site(longitude, latitude, altitude)
observation.set_fuse(fuso)
observation.create_plan( dados['coord'], dados['object'], dados['comment'], horain, horafin, intinf, limalt=30*u.deg, size=limdist)
#print observation.plan(dados['coord'], tempoin, dados['object'], dados['comment'], culmination=None, limalt=30.0*u.deg, size=11.0*u.arcmin)

