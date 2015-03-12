# coding=UTF-8

import numpy as np
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, Latitude, Longitude, FK5, EarthLocation, Angle

######################################################################

class Observation(object):

    def __init__(self):
        self.limalt = 0
        self.limdist = 0
        
    def read(self, arquivo):
        self.dados = np.loadtxt(arquivo, dtype={'names': ('objeto', 'data', 'hora', 'afh', 'afm', 'afs', 'ded', 'dem', 'des', 'mag'), 
            'formats': ('S30', 'S30', 'S30', 'i4', 'i4', 'f8', 'i4', 'i4', 'f8', 'S30')})

        ### leitura data da ocultacao ###
        tempo = np.core.defchararray.add(self.dados['data'], [' '])
        tempo = np.core.defchararray.add(tempo, self.dados['hora'])
        self.time = Time(tempo, format='iso', scale='utc')
        
        ### leitura das coordenadas dos alvos ###
        alfa = np.core.defchararray.add(np.array(map(str, self.dados['afh'])), [' '])
        alfa = np.core.defchararray.add(alfa, np.array(map(str, self.dados['afm'])))
        alfa = np.core.defchararray.add(alfa, [' '])
        alfa = np.core.defchararray.add(alfa, np.array(map(str, self.dados['afs'])))
        delta = np.core.defchararray.add(np.array(map(str, self.dados['ded'])), [' '])
        delta = np.core.defchararray.add(delta, np.array(map(str, self.dados['dem'])))
        delta = np.core.defchararray.add(delta, [' '])
        delta = np.core.defchararray.add(delta, np.array(map(str, self.dados['afs'])))
        coor = np.core.defchararray.add(alfa, [' '])
        coor = np.core.defchararray.add(coor, delta)
        self.coordenadas = SkyCoord(coor, frame='fk5', unit=(u.hourangle, u.degree))
        
    def set_sitio(self,longitude, latitude, altitude):
        ### leitura das coordenadas do local ###
        self.sitio = EarthLocation(longitude, latitude, altitude)

    def set_time_lim(self, horain, horafin, intinf):
        ### leitura dos instantes limites ###
        self.tempoin = Time(horain, format='iso', scale='utc', location=self.sitio) - self.utcoffset
        self.tempofin = Time(horafin, format='iso', scale='utc', location=self.sitio) + TimeDelta(0.5, format='sec')  - self.utcoffset
        self.intval = TimeDelta(intinf*60, format='sec')

    def set_fuse(self, fuso):
        ### leitura do fuso horario ###
        self.utcoffset = TimeDelta(fuso*3600, format='sec', scale='tai')

    def close_obj(self, limdist):
        tamcampo = Angle(limdist*u.arcmin)
        lista_campos = []
        lista_coords = []
        ### identificando objetos no mesmo campo
        for idx, value in enumerate(self.coordenadas[:-1]):
            coordcampo = value
            sep = self.coordenadas[idx].separation(self.coordenadas[idx+1:])
            close = [i + 1 + idx for i in np.arange(len(sep)) if sep[i] < tamcampo]
            combs = []
            for i in np.arange(len(close),0, -1):
            els = [list(x) for x in itertools.combinations(close, i)]
            combs.append(els)
            for i in combs:
                midpoint
#                for idy, valuey in enumerate(sep):
#                    if valuey < tamcampo:
#                        close = [idx + 1 + idy for i in sep if i < tamcampo]
#                        print va
            else:
                d = False
                for i in lista_campos:
                    e = idx in i
                    d = bool(d + e)
                    if d == False:
                        lista_campos.append([idx])
    
    def midpoint_coord(self, list_idx):
        x = np.mean([np.cos(i.dec)*np.cos(i.ra) for i in list_idx])
        y = np.mean([np.cos(i.dec)*np.sin(i.ra) for i in list_idx])
        z = np.mean([np.sin(i.dec) for i in list_idx])
        delta = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))*u.rad
        alfa = np.arctan(y/x)*u.rad
        mean_coord = SkyCoord(alfa, delta, frame='fk5')
        return mean_coord
        
    def precess(self, time):
        ### precessando as coordenadas para a data #############
        fk5_data = FK5(equinox=time)
        self.coord_prec = self.coordenadas.transform_to(fk5_data)
        
    def culmination(self):
        ### calculando hora da culminacao ###
        meio_noite = self.tempoin + (self.tempofin - self.tempoin)/ 2
        meio_noite.delta_ut1_utc = 0
        dif_h_sid = self.coord_prec.ra - meio_noite.sidereal_time('mean')
        dif_h_sid.wrap_at('180d', inplace=True)
        dif_h_sol = dif_h_sid * (23 + 56/60 + 4.0916/3600) / 24
        dif = TimeDelta(dif_h_sol.hour*60*60, format='sec')
        self.culminacao = meio_noite + dif + self.utcoffset
        self.culminacao.delta_ut1_utc = 0
        
    def sky_time(self, limalt):
        ### calculando tempo de ceu ###
        hangle_lim = np.arccos((np.cos((90-limalt)*u.deg) - np.sin(self.coord_prec.dec)*np.sin(self.sitio.latitude)) / np.cos(self.coord_prec.dec)*np.cos(self.sitio.latitude))
        tsg_lim = self.coord_prec.ra + hangle_lim
        dtsg_lim = tsg_lim - self.culminacao.sidereal_time('mean')
        dtsg_lim.wrap_at(360 * u.deg, inplace=True)
        dtsg_lim_sol = dtsg_lim * (23 + 56/60 + 4.0916/3600) / 24
        dtsg_np = TimeDelta(dtsg_lim_sol.hour*60*60, format='sec')
        self.nascer = self.culminacao - dtsg_np
        self.poente = self.culminacao + dtsg_np
        
    def altura_time(self, instante):
        instante.delta_ut1_utc = 0
        hourangle = instante.sidereal_time('mean') - self.coord_prec.ra
        distzen = np.arccos(np.sin(self.coord_prec.dec)*np.sin(self.sitio.latitude) + np.cos(self.coord_prec.dec)*np.cos(self.sitio.latitude)*np.cos(hourangle))
        altura = 90*u.deg - distzen
        time_rest = self.poente - instante
        return altura, time_rest
        
    def create_plan(self, horain, horafin, intinf):
        self.set_time_lim(horain, horafin, intinf)
        self.precess(self.tempoin)
        self.culmination()
        self.sky_time(self.limalt)
        nome = 'Plano_{}'.format(horain.split(' ')[0])
        output = open(nome, 'w')
        output.write('Plano de observação da noite {}\n\n'.format(horain.split(' ')[0]))
        output.write('Latitude: {}  Longitude: {}\nLimite de altura: {} deg\nTamanho do campo: {} arcmin\n\n'.format(self.sitio.latitude, self.sitio.longitude, self.limalt, self.limdist))
        instante = self.tempoin
        while instante <= self.tempofin:
            altura, time_rest = self.altura_time(instante)
            #### imprime os dados de cada objeto para cada instante ####
            output.write('\n---LT: {} (UT: {})----------------------------------------------------------------\n'.format((instante + self.utcoffset).iso.split(' ')[1][0:5], instante.iso.split(' ')[1][0:5]))
            for idx, val in enumerate(altura):
                if altura[idx] > limalt*u.deg:
                    output.write('{} {} ({})\n\tRA:{:02.0f} {:02.0f} {:07.4f}\tDEC: {:+03.0f} {:02.0f} {:06.3f}\n\tAltura: {:.1f}\n\tCulminacao: {} LT\n\tTempo restante para limite: {:02d}:{:02d}\n\n'\
.format(self.dados['objeto'][idx], self.dados['data'][idx], self.dados['mag'][idx], self.coordenadas[idx].ra.hms.h, self.coordenadas[idx].ra.hms.m, self.coordenadas[idx].ra.hms.s,
self.coordenadas[idx].dec.dms.d, self.coordenadas[idx].dec.dms.m, self.coordenadas[idx].dec.dms.s, altura[idx], self.culminacao[idx].iso.split(' ')[1][0:5],
int(time_rest[idx].sec/3600), int((time_rest[idx].sec - int(time_rest[idx].sec/3600)*3600)/60))) 
            instante = instante + self.intval
        output.write('\n---Observabilidade dos Alvos----------------------------------------------------------------\n')
        for idx, val in enumerate(self.nascer):
            output.write('{:10s} {} {}, RA: {:02.0f} {:02.0f} {:07.4f}, DEC: {:+03.0f} {:02.0f} {:06.3f}, Nascer: {} TL, Culminacao: {} TL, Poente: {} TL\n'\
.format(self.dados['objeto'][idx], self.dados['data'][idx], self.dados['hora'][idx], self.coordenadas[idx].ra.hms.h, self.coordenadas[idx].ra.hms.m, self.coordenadas[idx].ra.hms.s,
self.coordenadas[idx].dec.dms.d, self.coordenadas[idx].dec.dms.m, self.coordenadas[idx].dec.dms.s, self.nascer[idx].iso.split(' ')[1][0:5], self.culminacao[idx].iso.split(' ')[1][0:5],
self.poente[idx].iso.split(' ')[1][0:5]))
        output.close()
        
#####################################################################

arquivo = 'alvos.txt'				#### arquivo de alvos
horain = '2015-03-08 18:00:00'			#### hora inicial local da observacao
horafin = '2015-03-09 06:00:00'			#### hora final local da observacao
intinf = 30					#### intervalo entre as informacoes (minutos)
fuso = -2					#### fuso horario do local
latitude = '-22 32 7.8'				#### latitude do local
longitude = '314 25 2.5'			#### longitude do local
altitude = 1864					#### altitude em metros
limalt = 30					#### limite de altura para mostrar (graus)
limdist = 11					#### limite de distancia para field-of-view (arcmin)

######################################################################

observation = Observation()
observation.read(arquivo)
observation.set_sitio(longitude, latitude, altitude)
observation.set_fuse(fuso)
observation.limalt = limalt
observation.create_plan(horain, horafin, intinf)

