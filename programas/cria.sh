atual=$(pwd)
listaescala="$atual/astrometry_findscale_01"
totallinhas=$(wc -l < lista_pastas)
numtotal=$(echo "$totallinhas"|bc)
count=0
while read escala erroesc pasta
do

  cd $pasta

  esca=$(printf "%.3f" $escala)
if [ $(echo "$erroesc > 0.005" | bc) -ne 0 ]; then
  erro=$(printf "%.3f" $erroesc)
 else
  erro="0.005"
fi

#escrevendo o arquivo dat

echo "/homeA/marcelo/2MASS/2MASS/                       | root directory under which the 2MASS catalogue sub-directories lay " >> PRAIA_astrometry_20_08.dat
echo "/homeA/marcelo/UCAC2/UCAC2/                       | root directory under which the UCAC2 catalogue sub-directories lay " >> PRAIA_astrometry_20_08.dat
echo "/homeA/marcelo/UCAC4/UCAC4/                       | root directory under which the UCAC4 catalogue sub-directories lay " >> PRAIA_astrometry_20_08.dat
echo "13_sedna.cat                                      | User reference catalogue (PRAIA format) " >> PRAIA_astrometry_20_08.dat
echo "2                                                 | Make reduction with user's reference catalogue ? 1 - yes;  2 - no " >> PRAIA_astrometry_20_08.dat
echo "output                                            | extracted header data from fits images" >> PRAIA_astrometry_20_08.dat
echo "targets                                           | targets input file: (RA,Dec), JD, target name " >> PRAIA_astrometry_20_08.dat
echo "fdp.dat                                           | Field Distortion Pattern data file" >> PRAIA_astrometry_20_08.dat
echo "4                                                 | n nearest points for FDP computations " >> PRAIA_astrometry_20_08.dat
echo "bpx.dat                                           | Bad pixel mask  (xmin xmax ymin ymax)  (common to all treated images)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_photometry_praia_star1                 | photometric statistics of each field" >> PRAIA_astrometry_20_08.dat
echo "astrometry_reduction_2mu2_praia_star1             | reduction statistics of each field for UCAC2 and 2MASS (original)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_reduction_2mu4_praia_star1             | reduction statistics of each field for UCAC4 and 2MASS (original)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_reduction_mp_praia_star1               | reduction statistics of each field for UCAC2 and 2MASS (t.p.+common p.m.)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_reduction_mp_med_praia_star1           | reduction statistics of each field for UCAC2 and 2MASS (t.p.+common and non-common p.m.)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_reduction_2mus_praia_star1             | reduction statistics of each field for User catalogue and 2MASS (original)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_2MASS_target_praia_star1               | target statistics for 2MASS (original) " >> PRAIA_astrometry_20_08.dat
echo "astrometry_UCAC2_target_praia_star1               | target statistics for UCAC2" >> PRAIA_astrometry_20_08.dat
echo "astrometry_UCAC4_target_praia_star1               | target statistics for UCAC4" >> PRAIA_astrometry_20_08.dat
echo "astrometry_2MASS_target_mpu2_praia_star1          | target statistics for 2MASS (t.p.+common p.m.)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_2MASS_target_mpu2_med_praia_star1      | target statistics for 2MASS (t.p.+common and non-common p.m.)" >> PRAIA_astrometry_20_08.dat
echo "astrometry_USER_target_praia_star1                | target statistics for User reference catalogue " >> PRAIA_astrometry_20_08.dat
echo "ucac2.red.xy                                      | extension of \"xy\" output files, UCAC2 reduction" >> PRAIA_astrometry_20_08.dat
echo "ucac4.red.xy                                      | extension of \"xy\" output files, UCAC4 reduction" >> PRAIA_astrometry_20_08.dat
echo "2mass.red.xy                                      | extension of \"xy\" output files  2MASS reduction" >> PRAIA_astrometry_20_08.dat
echo "2mass.rmp.xy                                      | extension of \"xy\" output files  2MASS (t.p.+common p.m.) reduction" >> PRAIA_astrometry_20_08.dat
echo "2mass.rme.xy                                      | extension of \"xy\" output files  2MASS (t.p.+common and non-common p.m.) reduction" >> PRAIA_astrometry_20_08.dat
echo "wfi.red.xy                                        | extension of \"xy\" output files, User reference catalog reduction" >> PRAIA_astrometry_20_08.dat
echo "2     0   0   0      0   0   0                    | field area for tp+pm reduction (key,dg,min,arcsec): key=1-2 CCD size;  key=0 areax=dg,m,s areay=dg,m,s " >> PRAIA_astrometry_20_08.dat
echo "15                                                | maximum number of non-common 2MASS stars for (t.p.+common and non-common p.m.) reduction" >> PRAIA_astrometry_20_08.dat
echo "01.0                                              | error radius for target identification (arcsec) " >> PRAIA_astrometry_20_08.dat
echo "$esca                                             | pixel scale (arcsec/pixel)\n" >> PRAIA_astrometry_20_08.dat
echo "$erro                                             | +/- error of pixel scale (arcsec/pixel)" >> PRAIA_astrometry_20_08.dat
echo "1                                                 | Automatic catalogue <-> (x,y) identification: 1 - same N-S/E-W orientation & pixel scale;  2 - mixed images " >> PRAIA_astrometry_20_08.dat
echo "65000                                             | ADU maximum counting cutoff for non-linear pixels (~saturation) (ex.: 32000, 65000, ...)" >> PRAIA_astrometry_20_08.dat
echo "+0000                                             | ADU minimum counting cutoff for non-linear pixels (~sky background) (ex.: -10, 0 , ...)" >> PRAIA_astrometry_20_08.dat
echo "1                                                 | Pixel physical counts: 0 = from image header or 1 = from user (here)" >> PRAIA_astrometry_20_08.dat
echo "1.0d0                                             | Pixel physical counts: bscale;  Pixel = bscale * matrix + bzero" >> PRAIA_astrometry_20_08.dat
echo "0.0d0                                             | Pixel physical counts: bzero ;  Pixel = bscale * matrix + bzero" >> PRAIA_astrometry_20_08.dat
echo "-99                                               | bitpix: -99 reads from image header; otherwise use 16, 32, 64, -32, -64 following FITS conventions" >> PRAIA_astrometry_20_08.dat
echo "0                                                 | litteendian x bigendian: byte-swap (0 = automatic; 1 = don't swap ; 2 = swap bytes)" >> PRAIA_astrometry_20_08.dat
echo "3                                                 | sky background flattening: degree of complete bi-variate 2-D polynomial model (1 - 15) (0 = no flattening)" >> PRAIA_astrometry_20_08.dat
echo "5                                                 | N for smoothing filter of (2N+1) channels in sky background computations (0 = no filter) (suggestion: N=5)" >> PRAIA_astrometry_20_08.dat
echo "0000.500                                          | sky background theshold factor: theshold = sky + FACTOR * sigma (objects ID)" >> PRAIA_astrometry_20_08.dat
echo "0.5  10.0                                         | minimum maximum FWHM (range of FWHMs) (objects ID)" >> PRAIA_astrometry_20_08.dat
echo "0500                                              | brightest 2MASS stars for cross-identification with brightest measured (x,y) objects" >> PRAIA_astrometry_20_08.dat
echo "15                                                | brightest measured (x,y) objects for cross-identification with brightest 2MASS stars" >> PRAIA_astrometry_20_08.dat
echo "00.5000                                           | error radius (arcsec) for cross-identification between brightest catalogue/measured objs" >> PRAIA_astrometry_20_08.dat
echo "0.200                                             | (O-C) cutoff for outliers in (RA,DEC) reductions with 2MASS (original)" >> PRAIA_astrometry_20_08.dat
echo "0.120                                             | (O-C) cutoff for outliers in (RA,DEC) reductions with UCAC2, UCAC4 & 2MASS corrected versions " >> PRAIA_astrometry_20_08.dat
echo "1                                                 | polynomial (x,y) <-> (X,Y) in (RA,DEC) reductions: 0 = 4 Ctes; 1 to 3 = complete order" >> PRAIA_astrometry_20_08.dat
echo "0                                                 | radial distortion of 3rd order (x,y) <-> (X,Y) in (RA,DEC) reduction: 0 = no; 3 = yes  " >> PRAIA_astrometry_20_08.dat
echo "0                                                 | radial distortion of 5th order (x,y) <-> (X,Y) in (RA,DEC) reduction: 0 = no; 5 = yes" >> PRAIA_astrometry_20_08.dat
echo "00 00                                             | range of fits images to reduce from extracted header data file" >> PRAIA_astrometry_20_08.dat
echo "****************************************************************************************************************************************************************" >> PRAIA_astrometry_20_08.dat

#count=$[$count + 1]
#porcentagemtot=$(echo "$count * 100 / $numtotal"|bc -l)
#echo "$porcentagemtot"

cd $atual

########

done < $listaescala
