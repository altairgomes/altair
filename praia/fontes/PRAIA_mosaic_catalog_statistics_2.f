c
c
c       PRAIA_mosaic_catalogue_statistics 
c
c
c       Produces tables useful for plots and statistical analysis of PRAIA
c       catalogs resulting from PRAIA global reductions.
c       
c
c       In the current version, it is necessary to furnish the file name lists of
c       the follwing file types:
c
c
c       - individual xy files (single reductions)
c       - mosaic reduction log files 
c       - final catalog (output of mosaic reductions, after multiple object filtering
c         and proper motion computation)
c
c
c       The statistics are made for each set of mosaics. But for three types of tables,
c      statistics are also furnished adding all mosaics together: star distribution per
c      magnitude, Gauss (x,y) errors as a function of magnitudes and catalogue errors
c      as a function of magnitude.
c
c
c 
c      
c      Last modification:  19/Sep/2011
c
c

      implicit real *8 (a-h,o-z)

      parameter(stdin=5,stdout=6,idimh=100)

      dimension ng(idimh),ne(idimh),ngal(idimh),neal(idimh)

      dimension mcata(4),jcata(4),mflag(6)

      dimension magal(idimh),mag(idimh),icata(4),kcata(4),iflag(6)

      dimension gx(idimh),gy(idimh),gx2(idimh),gy2(idimh)

      dimension errx(idimh),erry(idimh),errx2(idimh),erry2(idimh)

      dimension gxal(idimh),gyal(idimh),gxal2(idimh),gyal2(idimh)

      dimension errax(idimh),erray(idimh),errax2(idimh),erray2(idimh)


      character*50 magall,gauall,errall,flgall,pmall,catall
      character*50 lisxys,lislog,liscat
      character*50 magstr,maggau,tabsii,tabsim,magerr
      character*50 catflg,catpmo,catcat,catcen

      character*1 isig

      character*50 mfits
      character*20 iobalv,ichfil

      character*300 catal,fxy,flog

      character*100 ask


c
c     Reads outputs file names for all mosaics (objects) statistics
c

      read (*,5) magall
      read (*,5) gauall
      read (*,5) errall

      read (*,5) flgall
      read (*,5) pmall
      read (*,5) catall


      read (*,*)


c
c     Initializing variables for all-mosaic-statistics
c

      ntall=0
      nmall=0
      nngall=0
      nneall=0


      do i=1,idimh

      magal(i)=0

      ngal(i)=0
      neal(i)=0

      gxal(i)=0.d0
      gyal(i)=0.d0
      gxal2(i)=0.d0
      gyal2(i)=0.d0

      errax(i)=0.d0
      erray(i)=0.d0
      errax2(i)=0.d0
      erray2(i)=0.d0

      enddo



      mflag(1)=0
      mflag(2)=0
      mflag(3)=0
      mflag(4)=0
      mflag(5)=0
      mflag(6)=0


      mcata(1)=0
      mcata(2)=0
      mcata(3)=0
      mcata(4)=0

      jcata(1)=0
      jcata(2)=0
      jcata(3)=0
      jcata(4)=0



c
c     Starting looping of mosaic blocks
c



 1    continue


c
c     Reads input data
c


      read (*,5) lisxys

      if (lisxys(1:1).eq.' ') go to 200

      read (*,5) lislog
      read (*,5) liscat


      read (*,5) magstr
      read (*,5) maggau
      read (*,5) tabsii
      read (*,5) tabsim
      read (*,5) magerr


      read (*,5) catflg
      read (*,5) catpmo
      read (*,5) catcat
      read (*,5) catcen


      read (*,*) dmagmi
      read (*,*) dmagma
      read (*,*) dmagbi

      read (*,*) erromi
      read (*,*) erroma


      read (*,*)

 5    format (a50)




c
c     Reads PRAIA catalogue.
c
c     Here we do statistics for:
c  
c     -  the distribution of the number of starts with regard to magnitude
c     -  the Gaussian errors
c     -  catalogue flags, proper motions, cross-match with other cataloguess
c        (UCAC2, 2MASS and USNOB1.0) 
c     -  catalogue intrinsic error (standard deviation of final position with 
c        respect to individual ones)
c    


      do i=1,idimh

      mag(i)=0

      ng(i)=0
      ne(i)=0

      gx(i)=0.d0
      gy(i)=0.d0
      gx2(i)=0.d0
      gy2(i)=0.d0

      errx(i)=0.d0
      erry(i)=0.d0
      errx2(i)=0.d0
      erry2(i)=0.d0

      enddo


      iflag(1)=0
      iflag(2)=0
      iflag(3)=0
      iflag(4)=0
      iflag(5)=0
      iflag(6)=0


      icata(1)=0
      icata(2)=0
      icata(3)=0
      icata(4)=0

      kcata(1)=0
      kcata(2)=0
      kcata(3)=0
      kcata(4)=0


c

      nm=0
      nng=0
      nne=0
      nt=0

      range=dmagma-dmagmi
      nbins=range/dmagbi

c


      open (10,file=liscat)

      open (20,file=catcen)


c

 10   continue


      noco=0
      ra=0.d0
      de=0.d0



      read (10,15,end=60) catal
 15   format(a300)


      open (11,file=catal)

 20   continue

      read(11,25,end=55) x,y,cseng,altu,fgcc,fumag,
     ?fumag2,cxmgu,codmg,codmg2,cxmgj,cxmgh,cxmgk,
     ?res2mg,resmg2,ermgj,ermgh,ermgk,copma,copmd,
     ?epma,epmd,coex,coey,cerau,cedeu,alfsic,
     ?delsic,nstaru,nfin,alsiuc,desiuc,ktir,oldra,
     ?oldde,kuth,kutm,zut,kutano,kutmes,kutdia,
     ?codj,iexps,ichfil,mfits,iobalv,nx,ny,
     ?numcom,egrxx,egryy,icat,ifla

 25   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,1x,
     ?f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?3(1x,i5),2(1x,f7.3),1x,2i1)


c
c     Average (RA,DEC) catalogue center
c

      noco=noco+1
      ra=ra+oldra
      de=de+oldde

c

      k=(nbins-1.d0)*(codmg-dmagmi)/range+1.d0

      nt=nt+1

c
c     Magnitude distribution
c


      if (codmg.lt.dmagmi) go to 30
      if (codmg.gt.dmagma) go to 30

      nm=nm+1

      mag(k)=mag(k)+1

    
      nmall=nmall+1

      magal(k)=magal(k)+1


 30   continue



c
c     Magnitude x Gaussian (x,y) errors
c


      if (coex.lt.erromi) go to 40
      if (coey.lt.erromi) go to 40
      if (coex.gt.erroma) go to 40
      if (coey.gt.erroma) go to 40

      nng=nng+1

      ng(k)=ng(k)+1

      gx(k)=gx(k)+coex
      gy(k)=gy(k)+coey

      gx2(k)=gx2(k)+coex**2
      gy2(k)=gy2(k)+coey**2


      nngall=nngall+1

      ngal(k)=ngal(k)+1

      gxal(k)=gxal(k)+coex
      gyal(k)=gyal(k)+coey

      gxal2(k)=gxal2(k)+coex**2
      gyal2(k)=gyal2(k)+coey**2




 40   continue


c
c     Catalog intrinsic error
c



      if (egrxx.lt.erromi) go to 50
      if (egryy.lt.erromi) go to 50
      if (egrxx.gt.erroma) go to 50
      if (egryy.gt.erroma) go to 50
      if (numcom.lt.2) go to 50

      nne=nne+1

      ne(k)=ne(k)+1

      errx(k)=errx(k)+egrxx
      erry(k)=erry(k)+egryy

      errx2(k)=errx2(k)+egrxx**2
      erry2(k)=erry2(k)+egryy**2


      nneall=nneall+1

      neal(k)=neal(k)+1

      errax(k)=errax(k)+egrxx
      erray(k)=erray(k)+egryy

      errax2(k)=errax2(k)+egrxx**2
      erray2(k)=erray2(k)+egryy**2



 50   continue


c
c     Multiple-object flag statistics 
c

      iflag(ifla+1)=iflag(ifla+1)+1



c
c     Proper motion origin statistics and star catalogue cross-match
c
c     1 - UCAC2
c     2 - 2MASS
c     3 - USNOB1.0
c     4 - Field star
c



      if (copma.gt.90.999d0.or.copmd.gt.90.999d0) then
      kcata(4)=kcata(4)+1
      icata(4)=icata(4)+1

      else

       if (icat.eq.9) then
        kcata(3)=kcata(3)+1
        icata(3)=icata(3)+1
       endif

       if (icat.eq.2) then
        kcata(2)=kcata(2)+1
        icata(2)=icata(2)+1
       endif

       if (icat.eq.1) then
        kcata(1)=kcata(1)+1
        icata(1)=icata(1)+1
       endif

      endif


      go to 20

     

c
c     Closes loop of PRAIA catalogues
c


 55   close (11)


      ra=ra/noco
      de=de/noco

      rac=ra
      dec=de
   
      iah=ra
      am=(ra-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0
      if (de.lt.0.d0) then
      isig='-'
      de=-de
      else
      isig='+'
      endif 
      idg=de
      dm=(de-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0 


      write (20,57) rac,dec,iah,iam,sa,isig,idg,idm,ds,catal
 57   format(2(1x,f13.9),2x,2(i2.2,1x),f7.4,1x,a1,2(i2.2,1x),f6.3,1x,
     ?a80)

      go to 10

 60   close (10)
      close (20)



c
c     Averages and standard deviations for magnitude bins
c


      do i=1,nbins
      call desvio (ng(i),gx(i),gx2(i))
      call desvio (ng(i),gy(i),gy2(i))
      call desvio (ne(i),errx(i),errx2(i))
      call desvio (ne(i),erry(i),erry2(i))
      enddo



c
c     Output table for star distribution as a function of magnitude
c

      open (20,file=magstr)


      do i=1,nbins

      dm=dmagmi+(i-0.5d0)*dmagbi

      per=100.d0*mag(i)/nm

      write (20,70) dm,per,mag(i),nm
 70   format(1x,f6.3,1x,f8.4,2(1x,i10))

      enddo

      close (20)


c
c     Output table for Gaussian (x,y) errors as a function of magnitude
c

      open (21,file=maggau)


      do i=1,nbins

      dm=dmagmi+(i-0.5d0)*dmagbi

      gx(i)=gx(i)*1000.d0
      gy(i)=gy(i)*1000.d0
      gx2(i)=gx2(i)*1000.d0
      gy2(i)=gy2(i)*1000.d0


      write (21,80) dm,gx(i),gy(i),gx2(i),gy2(i),ng(i),nng
 80   format(1x,f6.3,4(1x,f6.0),2(1x,i10))

      enddo

      close (21)




c
c     Output table for Intrinsic catalog errors as a function of magnitude
c

      open (22,file=magerr)


      do i=1,nbins

      dm=dmagmi+(i-0.5d0)*dmagbi

      errx(i)=errx(i)*1000.d0
      erry(i)=erry(i)*1000.d0
      errx2(i)=errx2(i)*1000.d0
      erry2(i)=erry2(i)*1000.d0



      write (22,85) dm,errx(i),erry(i),errx2(i),erry2(i),ne(i),nne
 85   format(1x,f6.3,4(1x,f6.0),2(1x,i10))

      enddo

      close (22)



c
c     Output table for catalog multiple-object flags
c


      open (23,file=catflg)


      f1=100.d0*iflag(1)/nt
      f2=100.d0*iflag(2)/nt
      f3=100.d0*iflag(3)/nt
      f4=100.d0*iflag(4)/nt
      f5=100.d0*iflag(5)/nt
      f6=100.d0*iflag(6)/nt

      write (23,*)
      write (23,*) 'Total No. of stars = ',nt
      write (23,*) 'Flag 0: No. and %  = ',iflag(1),f1,'%'
      write (23,*) 'Flag 1: No. and %  = ',iflag(2),f2,'%'
      write (23,*) 'Flag 2: No. and %  = ',iflag(3),f3,'%'
      write (23,*) 'Flag 3: No. and %  = ',iflag(4),f4,'%'
      write (23,*) 'Flag 4: No. and %  = ',iflag(5),f5,'%'
      write (23,*) 'Flag 5: No. and %  = ',iflag(6),f6,'%'
      write (23,*)

      close (23)



c
c     Output table for catalog star origin (UCAC2, 2MASS, B1 or field star)
c


      open (24,file=catcat)


      f1=100.d0*icata(1)/nt
      f2=100.d0*icata(2)/nt
      f3=100.d0*icata(3)/nt
      f4=100.d0*icata(4)/nt


      write(24,*)
      write(24,*) 'Total No. of stars   = ',nt
      write(24,*) 'UCAC2 star: N and %  = ',icata(1),f1,'%'
      write(24,*) '2MASS star: N and %  = ',icata(2),f2,'%'
      write(24,*) 'USNOB star: N and %  = ',icata(3),f3,'%'
      write(24,*) 'FIELD star: N and %  = ',icata(4),f4,'%'
      write(24,*)

      close (24)




c
c     Output table for proper motion: UCAC2 or 1rst epoch (2MASS, B1)
c     or no proper motion
c


      open (25,file=catpmo)


      f1=100.d0*kcata(1)/nt
      f2=100.d0*kcata(2)/nt
      f3=100.d0*kcata(3)/nt
      f4=100.d0*kcata(4)/nt


      write(25,*)
      write(25,*) 'Total No. of stars   = ',nt
      write(25,*) 'UCAC2 p.m.: N and %  = ',kcata(1),f1,'%'
      write(25,*) '2MASS p.m.: N and %  = ',kcata(2),f2,'%'
      write(25,*) 'USNOB p.m.: N and %  = ',kcata(3),f3,'%'
      write(25,*) 'No    p.m.: N and %  = ',kcata(4),f4,'%'
      write(25,*)

      close (25)




c
c     Feeds all-mosaic statistics for flags, proper motions and
c     star-catalog-belonging
c

      ntall=ntall+nt

      mflag(1)=mflag(1)+iflag(1)
      mflag(2)=mflag(2)+iflag(2)
      mflag(3)=mflag(3)+iflag(3)
      mflag(4)=mflag(4)+iflag(4)
      mflag(5)=mflag(5)+iflag(5)
      mflag(6)=mflag(6)+iflag(6)


      jcata(1)=jcata(1)+icata(1)
      jcata(2)=jcata(2)+icata(2)
      jcata(3)=jcata(3)+icata(3)
      jcata(4)=jcata(4)+icata(4)


      mcata(1)=mcata(1)+kcata(1)
      mcata(2)=mcata(2)+kcata(2)
      mcata(3)=mcata(3)+kcata(3)
      mcata(4)=mcata(4)+kcata(4)



c
c     Reads PRAIA xys from individual CCD reductions 
c
c     Here we do statistics for:
c  
c     -  sigmas (mean error) of individual reductions
c    
c



      nf=0
      nref=0

      xsi=0.d0
      ysi=0.d0
      xsi2=0.d0
      ysi2=0.d0



      open (11,file=lisxys)

c

 100  continue


      read (11,15,end=120) fxy


      open (12,file=fxy)

 115  continue

      read(12,110,err=115) x,y,cseng,altu,fgcc,fumag,
     ?fumag2,cxmgu,codmg,codmg2,cxmgj,cxmgh,cxmgk,
     ?res2mg,resmg2,ermgj,ermgh,ermgk,copma,copmd,
     ?epma,epmd,coex,coey,cerau,cedeu,alfsic,
     ?delsic,nstaru,nfin,alsiuc,desiuc,ktir,oldra,
     ?oldde,kuth,kutm,zut,kutano,kutmes,kutdia,
     ?codj,iexps,ichfil,mfits,iobalv,nx,ny

 110  format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,1x,
     ?f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))


      close (12)



c
c     Mean errors (Sa,Sd) from individual solutions
c


      if (alfsic.lt.erromi) go to 100
      if (delsic.lt.erromi) go to 100
      if (alfsic.gt.erroma) go to 100
      if (delsic.gt.erroma) go to 100

      nf=nf+1

      nref=nref+nfin

      xsi=xsi+alfsic
      ysi=ysi+delsic

      xsi2=xsi2+alfsic**2
      ysi2=ysi2+delsic**2


      go to 100


 120  close (11)




c
c     Averages and standard deviations
c


      nref=nref/nf

      call desvio (nf,xsi,xsi2)
      call desvio (nf,ysi,ysi2)



c
c
c     Output table for mean errors (Sa,Sd) from individual solutions
c

      open (26,file=tabsii)


      xsi=xsi*1000.d0
      ysi=ysi*1000.d0
      xsi2=xsi2*1000.d0
      ysi2=ysi2*1000.d0



      write (26,130) xsi,ysi,xsi2,ysi2,nref,nf
 130  format(4(1x,f6.0),2(1x,i10))


      close (26)





c
c     Reads PRAIA logs from mosaic global reductions 
c
c     Here we do statistics for:
c  
c   -  average (O-C)s for UCAC2 stars in the tangent plane before global reduction
c   -  standard deviations about these averages before global reduction
c   -  sigmas (mean errors) for UCAC2 stars in the tangent plane after global reduction
c    
c


      nmo=0   
      nref=0

      xmi=0.d0
      ymi=0.d0
      xsi=0.d0
      ysi=0.d0

      xsf=0.d0
      ysf=0.d0


      open (11,file=lislog)

c

 140  continue


      read (11,15,end=180) flog

      open (12,file=flog)



 150  continue
  
      ask=''

      read (12,155) ask
 155  format(a100)

      if (ask(1:30).ne.' Mean offsets (RA,DEC) before ') go to 150

      do i=1,48
      ask(i:i)=' '
      enddo


      read (ask,*,end=156) xmai,ymai

      go to 160

 156  read (12,*) ymai

 160  continue

      read (12,155) ask

      if (ask(1:30).ne.' Sigma offsets (RA,DEC) before') go to 160

      do i=1,49
      ask(i:i)=' '
      enddo

      read (ask,*) xsai,ysai


 165  continue

      read (12,155,end=170) ask

      if (ask(1:30).eq.' Number of used UCAC2 stars fo') then
      do i=1,42
      ask(i:i)=' '
      enddo
      read (ask,*) nnn
      endif


      if (ask(1:30).ne.' Sigma offsets (RA,DEC) before') go to 165

      do i=1,49
      ask(i:i)=' '
      enddo

      read (ask,*) xsaf,ysaf

      go to 165

 170  close (12)



c
c     Averages of means and of (Sa,Sd) of (O-C)s from global solutions
c


      nmo=nmo+1

      nref=nref+nnn

      xmi=xmi+xmai  
      ymi=ymi+ymai  

      xsi=xsi+xsai     
      ysi=ysi+ysai

      xsf=xsf+xsaf     
      ysf=ysf+ysaf



      go to 140



 180  close (11)




c
c     Averages and standard deviations 
c


      nref=nref/nmo

      xmi=xmi/nmo
      ymi=ymi/nmo

      xsi=xsi/nmo      
      ysi=ysi/nmo 

      xsf=xsf/nmo      
      ysf=ysf/nmo 



c
c
c     Output table for averages and standard deviations of (O-C)s for UCAC2 stars
c     in the tangent plane before and after mosaic global solutions
c

      open (27,file=tabsim)


      xmi=xmi*1000.d0
      ymi=ymi*1000.d0

      xsi=xsi*1000.d0
      ysi=ysi*1000.d0

      xsf=xsf*1000.d0
      ysf=ysf*1000.d0



      write (27,190) xmi,ymi,xsi,ysi,xsf,ysf,nmo,nref
 190  format(6(1x,f6.0),2(1x,i10))


      close (27)


c
c     Looping for another block of statistics
c


      go to 1


 200  continue



c
c     All-mosaic/objects statistics
c



c
c     Averages and standard deviations for magnitude bins
c


      do i=1,nbins
      call desvio (ngal(i),gxal(i),gxal2(i))
      call desvio (ngal(i),gyal(i),gyal2(i))
      call desvio (neal(i),errax(i),errax2(i))
      call desvio (neal(i),erray(i),erray2(i))
      enddo



c
c     Output table for star distribution as a function of magnitude
c

      open (30,file=magall)


      do i=1,nbins

      dm=dmagmi+(i-0.5d0)*dmagbi

      per=100.d0*magal(i)/nmall

      write (30,210) dm,per,magal(i),nmall
 210  format(1x,f6.3,1x,f8.4,2(1x,i10))

      enddo

      close (30)


c
c     Output table for Gaussian (x,y) errors as a function of magnitude
c

      open (31,file=gauall)


      do i=1,nbins

      dm=dmagmi+(i-0.5d0)*dmagbi

      gxal(i)=gxal(i)*1000.d0
      gyal(i)=gyal(i)*1000.d0

      gxal2(i)=gxal2(i)*1000.d0
      gyal2(i)=gyal2(i)*1000.d0



      write (31,220) dm,gxal(i),gyal(i),gxal2(i),gyal2(i),ngal(i),nngall
 220  format(1x,f6.3,4(1x,f6.0),2(1x,i10))

      enddo

      close (31)




c
c     Output table for Intrinsic catalog errors as a function of magnitude
c

      open (32,file=errall)


      do i=1,nbins

      dm=dmagmi+(i-0.5d0)*dmagbi

      errax(i)=errax(i)*1000.d0
      erray(i)=erray(i)*1000.d0

      errax2(i)=errax2(i)*1000.d0
      erray2(i)=erray2(i)*1000.d0



      write (32,230) dm,errax(i),erray(i),errax2(i),erray2(i),neal(i),
     ?nneall
 230  format(1x,f6.3,4(1x,f6.0),2(1x,i10))

      enddo

      close (32)



c
c     Output table for catalog multiple-object flags
c


      open (33,file=flgall)


      f1=100.d0*mflag(1)/ntall
      f2=100.d0*mflag(2)/ntall
      f3=100.d0*mflag(3)/ntall
      f4=100.d0*mflag(4)/ntall
      f5=100.d0*mflag(5)/ntall
      f6=100.d0*mflag(6)/ntall


      write(33,*)
      write(33,*)'Total No. of stars = ',ntall
      write(33,*)'Flag 0: No. and %  = ',mflag(1),f1,'%'
      write(33,*)'Flag 1: No. and %  = ',mflag(2),f2,'%'
      write(33,*)'Flag 2: No. and %  = ',mflag(3),f3,'%'
      write(33,*)'Flag 3: No. and %  = ',mflag(4),f4,'%'
      write(33,*)'Flag 4: No. and %  = ',mflag(5),f5,'%'
      write(33,*)'Flag 5: No. and %  = ',mflag(6),f6,'%'
      write(33,*)

      close (33)



c
c     Output table for catalog star origin (UCAC2, 2MASS, B1 or field star)
c


      open (34,file=catall)


      f1=100.d0*jcata(1)/ntall
      f2=100.d0*jcata(2)/ntall
      f3=100.d0*jcata(3)/ntall
      f4=100.d0*jcata(4)/ntall

      write(34,*)
      write(34,*)'Total No. of stars = ',ntall
      write(34,*)'UCAC2 star: N , %  = ',jcata(1),f1,'%'
      write(34,*)'2MASS star: N , %  = ',jcata(2),f2,'%'
      write(34,*)'USNOB star: N , %  = ',jcata(3),f3,'%'
      write(34,*)'FIELD star: N , %  = ',jcata(4),f4,'%'
      write(34,*)

      close (34)




c
c     Output table for proper motion: UCAC2 or 1rst epoch (2MASS, B1)
c     or no proper motion
c


      open (35,file=pmall)

      f1=100.d0*mcata(1)/ntall
      f2=100.d0*mcata(2)/ntall
      f3=100.d0*mcata(3)/ntall
      f4=100.d0*mcata(4)/ntall

      write(35,*)
      write(35,*)'Total No. of stars = ',ntall
      write(35,*)'UCAC2 p.m.: N , %  = ',mcata(1),f1,'%'
      write(35,*)'2MASS p.m.: N , %  = ',mcata(2),f2,'%'
      write(35,*)'USNOB p.m.: N , %  = ',mcata(3),f3,'%'
      write(35,*)'No    p.m.: N , %  = ',mcata(4),f4,'%'
      write(35,*)

      close (35)


c


      write (*,*) 
      write (*,*) 'Execution terminated successfully.'
      write (*,*) 


      end





c
c
c
c     Mean and standard deviation
c
c
c     entry:
c
c        xvam  = sum of individual values
c        xvas  = sum of square of individual values
c
c     output:
c
c        xvam  = mean
c        xvas  = standard deviation about that mean
c



      subroutine desvio (nest,xvam,xvas)

      IMPLICIT REAL *8 (A-H,O-Z)

c

      dnove=99.999d0

c

      if (nest.eq.0) then
      xvam=0.d0
      xvas=0.d0
      return
      endif

      EXMED=XVAM/NEST

c

      if (nest.eq.1) then

      xvas=0.d0

      return

      endif

c

      if (nest.eq.2) then

      xvas=dsqrt(dabs(2.d0*xvas-xvam**2))
      XVAM=EXMED

      return

      endif

c

      raiz=XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2

      if (raiz.lt.0.d0) then

      xvas=0.d0

      else

      XVAS=DSQRT(raiz/(NEST-1.D0))

      endif
c

      XVAM=EXMED


      return

      end


