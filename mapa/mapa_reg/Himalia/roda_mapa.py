from mapa import Map

mp = Map()
mp.infile('mapa_in_03mar.dat')
mp.create_label()
mp.calcfaixa(step=10)
mp.geramapa()
