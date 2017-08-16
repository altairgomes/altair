import numpy as np
import skvideo.io
import cv2
from astropy.io import fits

img = skvideo.io.VideoCapture('./dv-Phoebe_hosoi-170706.avi')

i = 1

while True:
    ret, frame = img.read()
    if not ret:
        break
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    outfile = 'teste/{}_{:04d}.fits'.format('Image', i)
    hdu = fits.PrimaryHDU(gray[::-1,:] - 32768)
    hdu.scale(type='int16', bzero=32768, bscale=1)
    hdu.writeto(outfile)
    i = i+1
    
print i
