from mapa import Map, geramapa

mp1 = Map()
mp1.infile('mapa_in_STE.dat')
mp1.create_label()
mp1.calcfaixa(step=10)
mp2 = Map()
mp2.infile('mapa_in_JPL.dat')
mp2.calcfaixa(step=10)
mp3 = Map()
mp3.infile('mapa_in_22fev.dat')
mp3.calcfaixa(step=10)
mp4 = Map()
mp4.infile('mapa_in_03mar.dat')
mp4.calcfaixa(step=10)

geramapa(mp1.star[0], mp1.datas_off[0], mp1.title[0], mp1.labelx[0], mp1.nameimg[0], lats=None, clat=None, fmt='eps', dpi=300, mapsize=None)
