import numpy as np
from astropy.io import fits

f = open('fits_bin_in.dat', 'r')

linhas = f.readlines()
arquivo = linhas[0].rsplit()[0]
binagem = int(linhas[1].rsplit()[0])
elimina = int(linhas[2].rsplit()[0])
nome_saida = linhas[3].rsplit()[0]

f.close()

g = np.loadtxt(arquivo, usecols=(0,0), dtype={'names':('arquivos', 'teste'), 'formats':('S50', 'S50')})
image_list = g['arquivos']

h = open('relatorio.dat', 'w')

h.write('Relatorio de binagem das imagens\n\nImage_Saida = Imagens_Entrada\n')

for i in np.arange(len(image_list) // binagem):
    image_concat = [ fits.getdata(image) for image in image_list[i*binagem + elimina: (i+1)*binagem - elimina] ]
    final_image = np.mean(image_concat, axis=0)
    outfile = '{}_{:05d}.fits'.format(nome_saida, i + 1)
    hdu = fits.PrimaryHDU(final_image)
    hdu.scale(type='int16')
    hdu.writeto(outfile)
    h.write('{} = {}\n'.format(outfile, image_list[i*binagem + elimina: (i+1)*binagem - elimina]))


#if (len(image_list) // binagem) * binagem < len(image_list):
#    i = len(image_list) // binagem
#    image_concat = [ fits.getdata(image) for image in image_list[i*binagem + elimina: len(image_list) - 1 - elimina] ]
#    final_image = np.mean(image_concat, axis=0)
#    outfile = '{}_{:05d}.fits'.format(nome_saida, i + 1)
#    hdu = fits.PrimaryHDU(final_image)
#    hdu.writeto(outfile)
#    h.write('{} = {}'.format(outfile, image_list[i*binagem + elimina: ]))

h.close()
