from mapa import Map, geramapa

latsdic = {}
clatdic = {}
mp1 = Map()
mp1.infile('mapa_in_STE.dat')
mp1.create_label()
mp1.calcfaixa(step=10)
l = mp1.latlon[mp1.datas_off[0].iso]['lats']
lats = [l['lon'], l['lat'], l['lon2'], l['lat2'], l['x'], l['y'], l['x2'], l['y2']]
c = mp1.latlon[mp1.datas_off[0].iso]['clat']
clat = [c['lon'], c['lat'], c['lab'], c['x'], c['y'], c['labx'], c['cxy']]
latsdic['STE'] = lats
clatdic['STE'] = clat
mp2 = Map()
mp2.infile('mapa_in_JPL.dat')
mp2.calcfaixa(step=10)
l = mp2.latlon[mp2.datas_off[0].iso]['lats']
lats = [l['lon'], l['lat'], l['lon2'], l['lat2'], l['x'], l['y'], l['x2'], l['y2']]
c = mp2.latlon[mp2.datas_off[0].iso]['clat']
clat = [c['lon'], c['lat'], c['lab'], c['x'], c['y'], c['labx'], c['cxy']]
latsdic['JPL'] = lats
clatdic['JPL'] = clat
mp3 = Map()
mp3.infile('mapa_in_22fev.dat')
mp3.calcfaixa(step=10)
l = mp3.latlon[mp3.datas_off[0].iso]['lats']
lats = [l['lon'], l['lat'], l['lon2'], l['lat2'], l['x'], l['y'], l['x2'], l['y2']]
c = mp3.latlon[mp3.datas_off[0].iso]['clat']
clat = [c['lon'], c['lat'], c['lab'], c['x'], c['y'], c['labx'], c['cxy']]
latsdic['22fev'] = lats
clatdic['22fev'] = clat
mp4 = Map()
mp4.infile('mapa_in_03mar.dat')
mp4.calcfaixa(step=10)
l = mp4.latlon[mp4.datas_off[0].iso]['lats']
lats = [l['lon'], l['lat'], l['lon2'], l['lat2'], l['x'], l['y'], l['x2'], l['y2']]
c = mp4.latlon[mp4.datas_off[0].iso]['clat']
clat = [c['lon'], c['lat'], c['lab'], c['x'], c['y'], c['labx'], c['cxy']]
latsdic['03mar'] = lats
clatdic['03mar'] = clat

geramapa(mp1.stars[0], mp1.datas_off[0], mp1.title[0], mp1.labelx[0], mp1.nameimg[0], lats=latsdic, clat=clatdic, mapsize=mp1.mapsize, fmt='eps', dpi=300)
