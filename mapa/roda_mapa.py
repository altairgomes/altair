from mapa_run import Map

mp = Map()
mp.infile()
mp.create_label()
mp.calcfaixa(erro=1000, step=10)
mp.geramapa()
