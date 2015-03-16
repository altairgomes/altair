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
        
    def read(self, arquivo):
        """
        """
        dados = np.loadtxt(arquivo, dtype={'names': ('objeto', 'data', 'hora', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'mag'), 
            'formats': ('S30', 'S30', 'S30', 'i4', 'i4', 'f8', 'i4', 'i4', 'f8', 'S30')})

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
        data = np.array([], dtype={'names': ('object', 'coord', 'time_occ', 'comment'),'formats': ('S30', object, object, 'S30')})
        for i in np.arange(len(dados)):
            a = [(dados['objeto'][i], coordenadas[i], time[i], dados['mag'][i])]
            b = np.asarray(a, dtype={'names': ('object', 'coord', 'time_occ', 'comment'),'formats': ('S30', object, object, 'S30')})
            data = np.append(data,b)
        return data
        
    def set_site(self,longitude, latitude, altitude):
        """
        """
        self.site = EarthLocation(longitude, latitude, altitude)

    def set_fuse(self, fuso):
        """
        """
        self.utcoffset = TimeDelta(fuso*3600, format='sec', scale='tai')

    def close_obj(self, coord, limdist):
        """
        """
        tamcampo = Angle(limdist*u.arcmin)
        lista_campos = []
        for idx, value in enumerate(coord):
            sep = coord[idx].separation(coord[idx+1:])
            close = [i + 1 + idx for i in np.arange(len(sep)) if sep[i] < tamcampo]
            combs = []
            for i in np.arange(len(close),0, -1):
                els = [list(x) for x in itertools.combinations(close, i)]
                combs.append(els)
            for i in combs:
                for l in i:
                    a = np.concatenate(([idx], l))
                    mean_coord, dist_center = self.midpoint_coord(coord[a])
                    if all(k <= tamcampo/2 for k in dist_center):
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
        return campos
    
    def midpoint_coord(self, coords):
        """
        """
        x = np.mean([np.cos(i.dec)*np.cos(i.ra) for i in coords])
        y = np.mean([np.cos(i.dec)*np.sin(i.ra) for i in coords])
        z = np.mean([np.sin(i.dec) for i in coords])
        delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))*u.rad
        alfa = np.arctan(y/x)*u.rad
        mean_coord = SkyCoord(alfa, delta, frame='fk5')
        dist_center = coords.separation(mean_coord)
        return mean_coord, dist_center
        
    def precess(self, coord, time):
        """
        """
        fk5_data = FK5(equinox=time)
        coord_prec = coord.transform_to(fk5_data)
        return coord_prec
        
    def sky_time(self, coord, tempo, limalt=0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        tempo.delta_ut1_utc = 0
        tempo.location = site
        dif_h_sid = coord.ra - tempo.sidereal_time('mean')
        dif_h_sid.wrap_at('180d', inplace=True)
        dif_h_sol = dif_h_sid * (23 + 56/60 + 4.0916/3600) / 24
        dif = TimeDelta(dif_h_sol.hour*60*60, format='sec')
        hangle_lim = np.arccos((np.cos(90*u.deg-limalt) - np.sin(coord.dec)*np.sin(site.latitude)) / np.cos(coord.dec)*np.cos(site.latitude))
        culminacao = tempo + dif + self.utcoffset
        culminacao.delta_ut1_utc = 0
        tsg_lim = coord.ra + hangle_lim
        dtsg_lim = tsg_lim - culminacao.sidereal_time('mean')
        dtsg_lim.wrap_at(360 * u.deg, inplace=True)
        dtsg_lim_sol = dtsg_lim * (23 + 56/60 + 4.0916/3600) / 24
        dtsg_np = TimeDelta(dtsg_lim_sol.hour*60*60, format='sec')
        sunrise = culminacao - dtsg_np
        sunset = culminacao + dtsg_np
        return culminacao, sunrise, sunset
        
    def altura_time(self, coord, instante, limalt=0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        instante.location = site
        instante.delta_ut1_utc = 0
        hourangle = instante.sidereal_time('mean') - coord.ra
        distzen = np.arccos(np.sin(coord.dec)*np.sin(site.latitude) + np.cos(coord.dec)*np.cos(site.latitude)*np.cos(hourangle))
        altura = 90*u.deg - distzen
        poente = self.sky_time(i, instante, limalt, site)[2]
        time_rest = poente - instante
        return altura, time_rest
        
    def plan(self, coord, hora, obj=None, date=None, comment=None, culmination=None, limalt=0*u.deg, site=None):
        if site == None:
            site = self.site
        if obj == None:
            obj = ['']*len(coord)
        if date == None:
            date = ['']*len(coord)
        if comment == None:
            comment = ['']*len(coord)
        a, b = [], []
        for i in coord:
            a.append(i.ra)
            b.append(i.dec)
        coord = SkyCoord(a,b, frame='fk5')
        if culmination == None:
            culmination = self.sky_time(coord, hora, site)[0]
        altura, time_rest = self.altura_time(coord, hora)
        obs_tot = np.array([], dtype={'names': ('coord', 'height', 'time_left', 'obj', 'date', 'comment', 'culmination'),'formats': (object, object, object, 'S30', 'S30', 'S30', object)})
        for i in np.arange(len(altura)):
            a = [(coord[i], altura[i], time_rest[i], obj[i], date[i], comment[i], culmination[i])]
            b = np.asarray(a, dtype={'names': ('coord', 'height', 'time_left', 'obj', 'date', 'comment', 'culmination'),'formats': (object, object, object, 'S30', 'S30', 'S30', object)})
            obs_tot = np.append(obs_tot,b)
        obs = obs_tot[obs_tot['height'] >= limalt]
        a = ''
        for i in np.arange(len(obs)):
            a = a + '{} {} ({})\n\tRA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n\tAltura: {:.1f}\n\tCulminacao: {} LT\n\tTempo restante para limite: {:02d}:{:02d}\n\n'\
.format(obs['obj'][i], obs['data'][i], obs['comment'][i], obs['coord'][i].ra.hms.h, obs['coord'][i].ra.hms.m, obs['coord'][i].ra.hms.s, 
obs['coord'][i].dec.dms.d, obs['coord'][i].dec.dms.m, obs['coord'][i].dec.dms.s, obs['height'][i], obs['culmination'][i].iso.split(' ')[1][0:5],
int(obs['time_left'][i].sec/3600), int((obs['time_left'][i].sec - int(obs['time_left'][i].sec/3600)*3600)/60))
        return a
                    
    def create_plan(self, coord, obj, date, comment, horain, horafin, intinf, limalt=0*u.deg, limdist=0*u.deg, site=None):
        """
        """
        if site == None:
            site = self.site
        tempoin = Time(horain, format='iso', scale='utc', location=site) - self.utcoffset
        tempofin = Time(horafin, format='iso', scale='utc', location=site) + TimeDelta(0.5, format='sec')  - self.utcoffset
        intval = TimeDelta(intinf*60, format='sec')
        for i in coord:
            a.append(i.ra)
            b.append(i.dec)
        coord = SkyCoord(a,b, frame='fk5')
        culminacao, nascer, poente = self.sky_time(coord, tempoin, limalt, site)
        nome = 'Plano_{}'.format(horain.split(' ')[0])
        output = open(nome, 'w')
        output.write('Plano de observação da noite {}\n\n'.format(horain.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nLimite de altura: {} deg\nTamanho do campo: {} arcmin\n\n'.format(site.latitude, site.longitude, limalt, limdist))
        instante = tempoin
        while instante <= tempofin:
            #### imprime os dados de cada objeto para cada instante ####
            text = self.plan(coord, instante, obj, date, comment, culminacao, site=site, limalt=0*u.deg)
            output.write('\n---LT: {} (UT: {})----------------------------------------------------------------\n'.format((instante + self.utcoffset).iso.split(' ')[1][0:5], instante.iso.split(' ')[1][0:5]))
            output.write(text) 
            instante = instante + intval
        output.write('\n---Observabilidade dos Alvos----------------------------------------------------------------\n')
        for idx, val in enumerate(nascer):
            output.write('{:10s} {}, RA: {:02.0f} {:02.0f} {:07.4f}, DEC: {:+03.0f} {:02.0f} {:06.3f}, Nascer: {} TL, Culminacao: {} TL, Poente: {} TL\n'\
.format(objeto[idx], data[idx], coord[idx].ra.hms.h, coord[idx].ra.hms.m, coord[idx].ra.hms.s,
coord[idx].dec.dms.d, np.absolute(coord[idx].dec.dms.m), np.absolute(coord[idx].dec.dms.s), nascer[idx].iso.split(' ')[1][0:5], culminacao[idx].iso.split(' ')[1][0:5],
poente[idx].iso.split(' ')[1][0:5]))
        output.close()
        
#####################################################################

arquivo = 'alvos.txt'				#### arquivo de alvos
horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
intinf = 30					#### intervalo entre as informacoes (minutos)
fuso = -3					#### fuso horario do local
latitude = '-22 32 7.8'				#### latitude do local
longitude = '314 25 2.5'			#### longitude do local
altitude = 1864					#### altitude em metros
limalt = 30					#### limite de altura para mostrar (graus)
limdist = 11					#### limite de distancia para field-of-view (arcmin)

#######################################################################

observation = Observation()
dados = observation.read(arquivo)
observation.set_site(longitude, latitude, altitude)
observation.set_fuse(fuso)
#observation.create_plan(horain, horafin, intinf)

