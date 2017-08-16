c
c
c    Program PRAIA_big_table_20
c
c
c    Produces mean results from PRAIA (RA,DEC) reductions for (RA,DEC),
c    offsets, magnitude, JD, UTC, sig_(ra,dec), N_cat.
c
c    Also gives information about number of frames, instrument used,
c    observer, ephemeris and general notes.
c
c    The mean results are distributed in groups acording to periods of 
c    time. The periods of time are defined by the user who can group the
c    observations in seconds, minutes, hours, days, months or even years.
c
c    The mean results may also be filtered by the position offsets. The
c    user defines the largest acceptable standard deviation of the offsets
c    about their mean. Elimination of outliers is done one-by-one observation
c    with the user defining which is the largest acceptable position offset
c    in units of standard deviation (2.0 sigma, 2.5 sigma, ...). The
c    filter procedure is naturally done automatically. 
c
c    The filtered (preserved) observations are stored in a big file (all 
c    observatins toghether) and also in separate files numbered 1 ... N
c    according to each separate group by time periods. The eliminated
c    (outliers) observations are also stored in a big file (all eliminated
c    observations) and in separate files 1 ... N for each group.
c
c
c    The input file is of type PRAIA_astrometry target. In this version,
c    the current format of the input file is based on PRAIA release 20_01 or
c    earlier. 
c
c
c
c    Last update:   18/Apr/2010   M. Assafin
c
c


      IMPLICIT REAL *8 (A-H,O-Z)

      parameter(idim=10001)

      dimension dx(idim),dy(idim),xob(idim),yob(idim),seng(idim),
     ?altu(idim),fgcc(idim),fumag(idim),fumag2(idim),xmgu(idim),
     ?cudmg(idim),cudmg2(idim),xmgj(idim),xmgh(idim),xmgk(idim),
     ?res2mg(idim),resmg2(idim),ermgj(idim),ermgh(idim),ermgk(idim),
     ?pma(idim),pmd(idim),epma(idim),epmd(idim),ex(idim),ey(idim),
     ?erau(idim),edeu(idim),alfsiu(idim),delsiu(idim),nstaru(idim),
     ?nfinau(idim),alsiu(idim),desiu(idim),ktirau(idim),ra(idim),
     ?de(idim),iuth(idim),iutm(idim),sut(idim),iutano(idim),
     ?iutmes(idim),iutdia(idim),dj(idim),iexps(idim),ichfil(idim),
     ?infits(idim),iobalv(idim),nx(idim),ny(idim)

      dimension ior(idim),dval(idim),ntbeg(idim),ntend(idim),
     ?itira(idim)

      character*20 ichfil,iobalv
      character*50 infits

      character*130 target,output,filter,outlie
      character*137 targf(idim),targo(idim),targ

      
      character*20 instru,observ,treat,packa,model,refcat,ephem,datrel,
     ?body

      character*40 notes

      character*59 ler

      character*1 isig,plic,aspas,menos(467)

      period(dyr,dmo,ddy,dhh,dmm,dsc)=dyr*365.25d0+dmo*30.5d0+ddy+
     ?dhh/24.d0+dmm/1440.d0+dsc*86400.d0


      data plic/"'"/
      data aspas/'"'/
 



c
c     Initial data
c

      dm0=2400000.5D0
      dj1=2400000.5D0
      dj2000=2451544.5d0

      idimc=59
      idimn=130
      idimnn=5
      idimme=467
      
  
      do i=1,idimme
      menos(i)='-'
      enddo
      


      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'PRAIA - Filtered, Grouped and Averaged Results')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
c
c
      read (*,15) output
      read (*,*)

c
c     Opens output table
c

      open (7,file=output)


c
c     Reads input data        
c

c     open (1,file='big_table_20_02.dat')
c


c
c     Table notes
c

      write (7,*) '                                                    '
      write (7,*) '  Table notes:'
      write (7,*)
      write (7,*) '- Offsets (Observed minus JPL ephemeris) are in mas.'
      write (7,*) '- (E_a,E_d) are (RA,DEC) position error estimates in'
      write (7,*) '  mas, based on the dispersion (standard deviation) '
      write (7,*) '  of the offsets.'
      write (7,*) '- N_fr is the number of observations per night.     '
      write (7,*) '- First block of offsets statistics regards filtered'
      write (7,*) '  results; second block regards to all offsets.'
      write (7,*) '- (Per) is the 0-100 percentage of eliminated'
      write (7,*) '  outliers.'
      write (7,*) '- (Sig) is the sigma theshold for elimination (mas).'
      write (7,*) '- (Fac) is the sigma factor for offset elimination. '
      write (7,*) '- (RA,DEC) are given in the usual notation (h, m, s '
      write (7,*) '  and dg, ...), but also in decimal form; they are  '
      write (7,*) '  averaged to the mean time of observations.'
      write (7,*) '- Mean time (UTC) is given in the usual form (h, m, '
      write (7,*) '  s, day, month, year), but also in years (decimal) '
      write (7,*) '  and in Julian Days.'
      write (7,*) '- Mag is the estimated magnitude of the object.     '
      write (7,*) '- (S_a,S_d) (in mas) are the (RA,DEC) mean error of '
      write (7,*) '  reference catalogue stars, from the astrometric   '
      write (7,*) '  solutions.'
      write (7,*) '- Ncat is the average No. of reference stars/frame. '
      write (7,*) '- Target is the usual name of the object.           '
      write (7,*) '- Pixel scale is given in arcseconds per pixel.     '
      write (7,*) '- Obs/Inst gives information about the telescope,   '
      write (7,*) '  instrument and site of observation with IAU code. '
      write (7,*) '- Observer or main observer is also listed.         '
      write (7,*) '- Data Treatment lists who has made the astrometry. '
      write (7,*) '- Package indicates the astrometric software used.  '
      write (7,*) '- Tangent Plane Model indicates the polynomial form '
      write (7,*) '  relating (x,y) and (RA,DEC) in the tangent plane. '
      write (7,*) '- The reference catalogue is also listed.           '
      write (7,*) '- The reference JPL ephemeris is also indicated.    '
      write (7,*) '- Data release informs the Table publication date.  '
      write (7,*) '- General remarks may be also furnished.            '
      write (7,*) '- Values marked -1 mean quantity not known.         '
      write (7,*) '                                                    '


c

      write (7,4) (menos(i),i=1,467)
 4    format (467a1)



      write (7,5) plic,aspas,aspas
 5    format('  off_ra  off_de   E_a   E_d  N_fr  off_ra  off_de   E_a  
     ? E_d  N_fr    Per  Sig   Fac   hh mm ss',8x,'dg  ',
     ?a1,3x,a1,'      h  m  s      dd mo year    Julian Date         Yea
     ?r',11x,'RA (ICRS)     DEC (ICRS)   Mag   S_a   S_d  Ncat Target Na
     ?me          ',a1,'/pixel Obs/Inst IAU code    Observer            
     ? Data treatment       Package              Tangent Plane Model  Re
     ?ference catalogue  Reference ephemeris  Data release         Gener
     ?al remarks')

      write (7,4) (menos(i),i=1,467)



c

 10   continue


 15   format(a130)
 16   format(a20)
 17   format(a40)

      target=''

      read (*,15,end=100) target
      read (*,15) filter
      read (*,15) outlie

      read (*,*) sigma
      sigma=sigma/1000.d0

      read (*,*) factor

      read (*,*) dyear,dmonth,dday,dhour,dminut,dsecon

      dtime=period(dyear,dmonth,dday,dhour,dminut,dsecon)

      read (*,16) body
      read (*,*)  pixel
      read (*,16) instru
      read (*,16) observ
      read (*,16) treat
      read (*,16) packa
      read (*,16) model
      read (*,16) refcat
      read (*,16) ephem
      read (*,16) datrel
      read (*,17) notes
      read (*,*)


c
c     opens output filtered file of data without grouping
c

      open (9,file=filter)
      open (19,file=outlie)


c
c     Stores PRAIA output target file data
c 

      open (1,file=target)


      do i=1,idim

      read  (1,18,end=19) dx(i),dy(i),xob(i),yob(i),seng(i),altu(i),
     ?fgcc(i),fumag(i),fumag2(i),xmgu(i),cudmg(i),cudmg2(i),xmgj(i),
     ?xmgh(i),xmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),
     ?pma(i),pmd(i),epma(i),epmd(i),ex(i),ey(i),erau(i),edeu(i),
     ?alfsiu(i),delsiu(i),nstaru(i),nfinau(i),alsiu(i),desiu(i),
     ?ktirau(i),ra(i),de(i),iuth(i),iutm(i),sut(i),iutano(i),iutmes(i),
     ?iutdia(i),dj(i),iexps(i),ichfil(i),infits(i),iobalv(i),nx(i),ny(i)

 18   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))

      enddo

 19   close (1)

      ntotal=i-1
      
c
c     In case the target file is empty
c

      if (ntotal.eq.0) then
      close (9)
      close (19)
      go to 10
      endif


c
c     Orders data by crescent time
c


      do i=1,ntotal
      ior(i)=i
      dval(i)=dj(i)
      enddo

      call ordem (ntotal,ior,dval)


c
c     Groups data by time intervals
c

      ngroup=1

      tend=dj(ior(1))+dtime

      ntbeg(1)=1


      do k=2,ntotal

      i=ior(k)

      if (dj(i).gt.tend) then
      ntend(ngroup)=k-1
      ngroup=ngroup+1
      ntbeg(ngroup)=k
      tend=dj(i)+dtime
      endif

      enddo

      ntend(ngroup)=ntotal


c
c     For each group, filters outliers.
c
c     Outliers elimination procedure:
c
c
c     1) Individual position offsets must not deviate from the offset mean by
c        more than a factor f of the present sigma; if the largest offset
c        deviation from the mean is above f*sigma, eliminate this single
c        observation and re-computes everything; after no more individual
c        observation is eliminated in this way, proceed to step 2.  
c
c     2) current sigma of position offsets must be within the user-given
c        threshold; if yes, stop the process; if not, eliminate one-by-one
c        the largest offsets deviating from the mean until the sigma falls
c        within the desired threshold.
c
c      
c     When one does not want to eliminate any observation, use large values
c     for the sigma threshold and sigma-factor.
c


      do i=1,ntotal
      itira(i)=0
      enddo


      do 25 k=1,ngroup


 20   n=0
      xoff=0.d0
      xoff2=0.d0
      yoff=0.d0
      yoff2=0.d0

      do 21 kk=ntbeg(k),ntend(k)

      if (itira(kk).ne.0) go to 21

      n=n+1

      i=ior(kk)

      xoff=xoff+dx(i)
      xoff2=xoff2+dx(i)**2

      yoff=yoff+dy(i)
      yoff2=yoff2+dy(i)**2

 21   continue


      if (n.le.2) go to 25

      call desvio (n,xoff,xoff2)
      call desvio (n,yoff,yoff2)



c
c     Filter criterium (1)
c

      xthre=factor*xoff2
      ythre=factor*yoff2

      dmaxx=-1.d14
      dmaxy=-1.d14

      do 22 kk=ntbeg(k),ntend(k)

      if (itira(kk).ne.0) go to 22

      i=ior(kk)

      difx=dabs(dx(i)-xoff)
      dify=dabs(dy(i)-yoff)


      if (difx.gt.dmaxx) then
      ikx=kk
      dmaxx=difx
      endif

      if (dify.gt.dmaxy) then
      iky=kk
      dmaxy=dify
      endif


 22   continue


      dmaxx=dmaxx/xthre
      dmaxy=dmaxy/ythre

      dmax=dmax1(dmaxx,dmaxy)

      if (dmaxx.ge.dmaxy) then
      ik=ikx
      else
      ik=iky
      endif

      if (dmax.gt.1.d0) then
      itira(ik)=1
      go to 20
      endif

c
c     Filter criterium (2)
c


      dmax=dmax1(xoff2,yoff2)

      if (xoff2.ge.yoff2) then
      ik=ikx
      else
      ik=iky
      endif

      if (dmax.gt.sigma) then
      itira(ik)=1
      go to 20
      endif


 25   continue



c
c     Constructs file names for groups
c

      do i=1,idimn
      if (target(i:i).ne.' ') ik=i
      enddo 

      ik=ik+1

      do i=1,ngroup
      targ=''
      targ=target
      targ(ik:ik+1)='_f'
      write (targ(ik+2:ik+idimnn+1),'(i5.5)') i
      targf(i)=targ
      targ(ik+1:ik+1)='o'
      targo(i)=targ
      enddo


c
c     Outputing data
c

      sigma=sigma*1.d3

      do 60 k=1,ngroup


c
c     Zeroing data
c 

      nfrall=0
      orall=0.d0
      odall=0.d0
      or2all=0.d0
      od2all=0.d0


      offra=0.d0
      offde=0.d0

      offra2=0.d0
      offde2=0.d0

      rac=0.d0
      dec=0.d0

      dju=0.d0
      
      dmag=0.d0

      siga=0.d0
      sigd=0.d0

      dnest=0.d0

      nframe=0

      open (8,file=targf(k))
      open (18,file=targo(k))

      do 30 kk=ntbeg(k),ntend(k)

      i=ior(kk)

      if (itira(kk).ne.0) then

c
c     Outilers target file of group 
c

      write  (18,18) dx(i),dy(i),xob(i),yob(i),seng(i),altu(i),
     ?fgcc(i),fumag(i),fumag2(i),xmgu(i),cudmg(i),cudmg2(i),xmgj(i),
     ?xmgh(i),xmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),
     ?pma(i),pmd(i),epma(i),epmd(i),ex(i),ey(i),erau(i),edeu(i),
     ?alfsiu(i),delsiu(i),nstaru(i),nfinau(i),alsiu(i),desiu(i),
     ?ktirau(i),ra(i),de(i),iuth(i),iutm(i),sut(i),iutano(i),iutmes(i),
     ?iutdia(i),dj(i),iexps(i),ichfil(i),infits(i),iobalv(i),nx(i),ny(i)


c
c     Outliers target file of all observations without grouping
c


      write  (19,18) dx(i),dy(i),xob(i),yob(i),seng(i),altu(i),
     ?fgcc(i),fumag(i),fumag2(i),xmgu(i),cudmg(i),cudmg2(i),xmgj(i),
     ?xmgh(i),xmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),
     ?pma(i),pmd(i),epma(i),epmd(i),ex(i),ey(i),erau(i),edeu(i),
     ?alfsiu(i),delsiu(i),nstaru(i),nfinau(i),alsiu(i),desiu(i),
     ?ktirau(i),ra(i),de(i),iuth(i),iutm(i),sut(i),iutano(i),iutmes(i),
     ?iutdia(i),dj(i),iexps(i),ichfil(i),infits(i),iobalv(i),nx(i),ny(i)


c
c     Statistics of all data
c


      nfrall=nfrall+1

      orall=orall+dx(i)
      odall=odall+dy(i)

      or2all=or2all+dx(i)**2
      od2all=od2all+dy(i)**2

      go to 30
      endif 



c
c     Filtered target file of group 
c

      write  (8,18) dx(i),dy(i),xob(i),yob(i),seng(i),altu(i),
     ?fgcc(i),fumag(i),fumag2(i),xmgu(i),cudmg(i),cudmg2(i),xmgj(i),
     ?xmgh(i),xmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),
     ?pma(i),pmd(i),epma(i),epmd(i),ex(i),ey(i),erau(i),edeu(i),
     ?alfsiu(i),delsiu(i),nstaru(i),nfinau(i),alsiu(i),desiu(i),
     ?ktirau(i),ra(i),de(i),iuth(i),iutm(i),sut(i),iutano(i),iutmes(i),
     ?iutdia(i),dj(i),iexps(i),ichfil(i),infits(i),iobalv(i),nx(i),ny(i)


c
c     Filtered target file of all observations without grouping
c


      write  (9,18) dx(i),dy(i),xob(i),yob(i),seng(i),altu(i),
     ?fgcc(i),fumag(i),fumag2(i),xmgu(i),cudmg(i),cudmg2(i),xmgj(i),
     ?xmgh(i),xmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),
     ?pma(i),pmd(i),epma(i),epmd(i),ex(i),ey(i),erau(i),edeu(i),
     ?alfsiu(i),delsiu(i),nstaru(i),nfinau(i),alsiu(i),desiu(i),
     ?ktirau(i),ra(i),de(i),iuth(i),iutm(i),sut(i),iutano(i),iutmes(i),
     ?iutdia(i),dj(i),iexps(i),ichfil(i),infits(i),iobalv(i),nx(i),ny(i)


c
c     Statistics of all data
c

      nfrall=nfrall+1

      orall=orall+dx(i)
      odall=odall+dy(i)

      or2all=or2all+dx(i)**2
      od2all=od2all+dy(i)**2

c
c     Statistics of filtered data
c


      nframe=nframe+1

      offra=offra+dx(i)
      offde=offde+dy(i)

      offra2=offra2+dx(i)**2
      offde2=offde2+dy(i)**2

      rac=rac+ra(i)
      dec=dec+de(i)

      dju=dju+dj(i)-dm0

      dmag=dmag+cudmg(i)

      siga=siga+alfsiu(i)
      sigd=sigd+delsiu(i)

      dnest=dnest+nfinau(i)


 30   continue

      close (8)
      close(18)

c
c     Averages
c


      call desvio (nframe,offra,offra2)
      call desvio (nframe,offde,offde2)

      call desvio (nfrall,orall,or2all)
      call desvio (nfrall,odall,od2all)

      
      orall=orall*1.d3
      odall=odall*1.d3

      or2all=or2all*1.d3
      od2all=od2all*1.d3


      offra=offra*1.d3
      offde=offde*1.d3

      offra2=offra2*1.d3
      offde2=offde2*1.d3

      rac=rac/nframe
      dec=dec/nframe

      dju=dju/nframe
      dju=dju+dm0

      dmag=dmag/nframe

      siga=siga/nframe
      sigd=sigd/nframe

      siga=siga*1.d3
      sigd=sigd*1.d3

      if (siga.gt.999.9d0) then
      siga=-1.d0
      endif

      if (sigd.gt.999.9d0) then
      sigd=-1.d0
      endif

      dnest=dnest/nframe

c
c     Computes time in UTC (gregorian calendar and year decimal notation)
c



      djm=dju-dj1

      call iau_jd2cal (dj1,djm,kutano,kutmes,kutdia,fd,jjj)

      hora=fd*24.d0
      kuth=hora
      kutm=(hora-kuth)*60.d0
      zut=((hora-kuth)*60.d0-kutm)*60.d0


      frano=(dju-dj2000)/365.25d0+2000.d0

      ler=''

      write (ler,35) kuth,kutm,zut,kutdia,kutmes,kutano,dju,frano
 35   format(2(i2.2,':'),f7.4,':',2(i2.2,':'),i4.4,':',f18.10,':',
     ?f15.10)


      do i=idimc,1,-1
      if (ler(i:i).ne.' ') go to 40
      enddo

 40   ii=i

      do i=1,ii-1
      if (ler(i:i).eq.' ') ler(i:i)='0'
      enddo

      do i=1,ii-1
      if (ler(i:i).eq.':') ler(i:i)=' '
      enddo

c
c     Computes (RA,DEC) in hexadecimal notation
c


      iah=rac
      am=(rac-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0
      if (dec.lt.0.d0) then
      isig='-'  
      dec=-dec
      else
      isig='+' 
      endif
      idg=dec
      dm=(dec-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0

      if (isig.eq.'-') dec=-dec

c
c     Writes output Table
c

      nest=dnest

      if (nest.le.0) then
      nest=-1
      endif


      if (nframe.eq.1) then
      offra2=-1.d0
      offde2=-1.d0
      endif


      if (nfrall.eq.1) then
      or2all=-1.d0
      od2all=-1.d0
      endif

      perc=100.d0-100.d0*nframe/nfrall


c     write (7,50) offra,offde,offra2,offde2,nframe,iah,iam,sa,isig,idg,
c    ?idm,ds,ler,rac,dec,dmag,siga,sigd,nest,body,pixel,instru,observ,
c    ?treat,packa,model,refcat,ephem,datrel,notes

c50   format(2(1x,f7.1),2(1x,f5.1),1x,i5,3x,2(i2.2,1x),f8.5,1x,a1,
c    ?2(i2.2,1x),f7.4,1x,a58,2(1x,f14.10),1x,f4.1,2(1x,f5.1),1x,i5,1x,
c    ?a20,1x,f7.4,8(1x,a20),1x,a40)



      write (7,50) offra,offde,offra2,offde2,nframe,orall,odall,or2all,
     ?od2all,nfrall,perc,sigma,factor,iah,iam,sa,isig,idg,idm,ds,ler,
     ?rac,dec,dmag,siga,sigd,nest,body,pixel,instru,observ,treat,packa,
     ?model,refcat,ephem,datrel,notes

 50   format(2(2(1x,f7.1),2(1x,f5.1),1x,i5),1x,f6.2,1x,f5.1,1x,f4.1,3x,
     ?2(i2.2,1x),f8.5,1x,a1,2(i2.2,1x),f7.4,1x,a58,2(1x,f14.10),1x,f4.1,
     ?2(1x,f5.1),1x,i5,1x,a20,1x,f7.4,8(1x,a20),1x,a40)

c


 60   continue


c

      close (9)
      close (19)
c

      go to 10

c

 100  continue

      write (7,4) (menos(i),i=1,467)

      write (7,*)

      close (7)

c

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,101)
 101  format (23x,'PRAIA - execution terminated successfuly.')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '


      end



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
      dneg=-1.d0

c

      if (nest.le.0) then
      xvam=dnove
      xvas=dneg
      return
      endif


c

      EXMED=XVAM/NEST

c

      if (nest.eq.1) then

      xvas=dneg

      return

      endif

c

      if (nest.eq.2) then

      xvas=dsqrt(dabs(xvas-2.d0*exmed**2))
      XVAM=EXMED

      return

      endif

c

      raiz=XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2

      if (raiz.lt.0.d0) then

      xvas=dneg

      else

      XVAS=DSQRT(raiz/(NEST-1.D0))

      endif
c

      XVAM=EXMED

      return

      end






c
c     Subroutine iau_jd2cal
c
c


      SUBROUTINE iau_jd2cal ( DJ1, DJ2, IY, IM, ID, FD, J )
*+
*  - - - - - - - - - - -
*   i a u _ J D 2 C A L
*  - - - - - - - - - - -
*
*  Julian Date to Gregorian year, month, day, and fraction of a day.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DJ1,DJ2     d     Julian Date (Notes 1, 2)
*
*  Returned:
*     IY          i     year
*     IM          i     month
*     ID          i     day
*     FD          d     fraction of day
*     J           i     status:
*                           0 = OK
*                          -1 = unacceptable date (Note 3)
*
*  Notes:
*
*  1) The earliest valid date is -68569.5 (-4900 March 1).  The
*     largest value accepted is 10^9.
*
*  2) The Julian Date is apportioned in any convenient way between
*     the arguments DJ1 and DJ2.  For example, JD=2450123.7 could
*     be expressed in any of these ways, among others:
*
*             DJ1            DJ2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*  3) In early eras the conversion is from the "Proleptic Gregorian
*     Calendar";  no account is taken of the date(s) of adoption of
*     the Gregorian Calendar, nor is the AD/BC numbering convention
*     observed.
*
*  Reference:
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P.Kenneth Seidelmann (ed), University Science Books (1992),
*     Section 12.92 (p604).
*
*  This revision:  2000 December 19
*
*  Copyright (C) 2001 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2
      INTEGER IY, IM, ID
      DOUBLE PRECISION FD
      INTEGER J

*  Minimum and maximum allowed JD
      DOUBLE PRECISION DJMIN, DJMAX
      PARAMETER ( DJMIN = -68569.5D0, DJMAX = 1D9 )

      INTEGER JD, L, N, I
      DOUBLE PRECISION DJ, D1, D2, F1, F2, F, D

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Check if date is acceptable.
      DJ = DJ1 + DJ2
      IF ( DJ.LT.DJMIN .OR. DJ.GT.DJMAX ) THEN
         J = -1
      ELSE
         J = 0

*     Copy the date, big then small, and re-align to midnight.
         IF ( DJ1 .GE. DJ2 ) THEN
            D1 = DJ1
            D2 = DJ2
         ELSE
            D1 = DJ2
            D2 = DJ1
         END IF
         D2 = D2 - 0.5D0

*     Separate day and fraction.
         F1 = MOD(D1,1D0)
         F2 = MOD(D2,1D0)
         F = MOD(F1+F2,1D0)
         IF ( F .LT. 0D0 ) F = F+1D0
         D = ANINT(D1-F1) + ANINT(D2-F2) + ANINT(F1+F2-F)
         JD = NINT(D) + 1

*     Express day in Gregorian calendar.
         L = JD + 68569
         N = ( 4*L ) / 146097
         L = L - ( 146097*N + 3 ) / 4
         I = ( 4000 * (L+1) ) / 1461001
         L = L - ( 1461*I ) / 4 + 31
         J = ( 80*L ) / 2447
         ID = L - ( 2447*J ) / 80
         L = J / 11
         IM = J + 2 - 12*L
         IY = 100 * ( N-49 ) + I + L

         FD = F
         J = 0
      END IF

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2001
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. The Software is made available free of charge for use by:
*
*     a) private individuals for non-profit research; and
*
*     b) non-profit educational, academic and research institutions.
*
*  3. Commercial use of the Software is specifically excluded from the
*     terms and conditions of this license.  Commercial use of the
*     Software is subject to the prior written agreement of the Board on
*     terms to be agreed.
*
*  4. The provision of any version of the Software under the terms and
*     conditions specified herein does not imply that future versions
*     will also be made available under the same terms and conditions.
*
*  5. The user may modify the Software for his/her own purposes.  The
*     user may distribute the modified software provided that the Board
*     is informed and that a copy of the modified software is made
*     available to the Board on request.  All modifications made by the
*     user shall be clearly identified to show how the modified software
*     differs from the original Software, and the name(s) of the
*     affected routine(s) shall be changed.  The original SOFA Software
*     License text must be present.
*
*  6. In any published work produced by the user and which includes
*     results achieved by using the Software, the user shall acknowledge
*     that the Software was used in producing the information contained
*     in such publication.
*
*  7. The user may incorporate or embed the Software into other software
*     products which he/she may then give away free of charge but not
*     sell provided the user makes due acknowledgement of the use which
*     he/she has made of the Software in creating such software
*     products.  Any redistribution of the Software in this way shall be
*     made under the same terms and conditions under which the user
*     received it from the SOFA Center.
*
*  8. The user shall not cause the Software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  9. The Software is provided to the user "as is" and the Board makes
*     no warranty as to its use or performance.   The Board does not and
*     cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Board makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*
*-----------------------------------------------------------------------

      END







      DOUBLE PRECISION FUNCTION dau_GMST82 ( DJ1, DJ2 )
*+
*  - - - - - - - - - - -
*   i a u _ G M S T 8 2
*  - - - - - - - - - - -
*
*  Universal Time to Greenwich Mean Sidereal Time (IAU 1982 model).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DJ1, DJ2     d      UT1 Julian Date (see note)
*
*  The result is the Greenwich Mean Sidereal Time (radians), in the
*  range 0 to 2pi.
*
*  Notes:
*
*  1  The UT1 epoch DJ1+DJ2 is a Julian Date, apportioned in any
*     convenient way between the arguments DJ1 and DJ2.  For example,
*     JD(UT1)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*             DJ1            DJ2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 and MJD methods are good compromises
*     between resolution and convenience.  The date & time method is
*     best matched to the algorithm used:  maximum accuracy (or, at
*     least, minimum noise) is delivered when the DJ1 argument is for
*     0hrs UT1 on the day in question and the DJ2 argument lies in the
*     range 0 to 1, or vice versa.
*
*  2  The algorithm is based on the IAU 1982 expression.  This is always
*     described as giving the GMST at 0 hours UT1.  In fact, it gives the
*     difference between the GMST and the UT, the steady 4-minutes-per-day
*     drawing-ahead of ST with respect to UT.  When whole days are ignored,
*     the expression happens to equal the GMST at 0 hours UT1 each day.
*
*  3  In this routine, the entire UT1 (the sum of the two arguments DJ1
*     and DJ2) is used directly as the argument for the standard formula,
*     the constant term of which is adjusted by 12 hours to take account
*     of the noon phasing of Julian Date.  The UT1 is then added, but
*     omitting whole days to conserve accuracy.
*
*  Called:
*     iau_ANP        normalize angle into range 0 to 2pi
*
*  References:
*
*  1  Transactions of the International Astronomical Union,
*     XVIII B, 67 (1983).
*
*  2  Aoki et al., Astron. Astrophys. 105, 359-361 (1982).
*
*  This revision:  2000 December 19
*
*  Copyright (C) 2001 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2

      DOUBLE PRECISION DS2R
      PARAMETER ( DS2R = 7.272205216643039903848712D-5 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Seconds per day, days per Julian century
      DOUBLE PRECISION DAYSEC, CENDAY
      PARAMETER ( DAYSEC = 86400D0, CENDAY = 36525D0 )

*  Coefficients of IAU 1982 GMST-UT1 model
      DOUBLE PRECISION A, B, C, D
      PARAMETER ( A = 24110.54841D0 - DAYSEC/2D0,
     :            B = 8640184.812866D0,
     :            C = 0.093104D0,
     :            D = -6.2D-6 )

*  Note: the first constant, A, has to be adjusted by 12 hours because
*  the UT1 is supplied as a Julian date, which begins at noon.

      DOUBLE PRECISION D1, D2, T, F

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Julian centuries since fundamental epoch.
      IF ( DJ1 .LT. DJ2 ) THEN
         D1 = DJ1
         D2 = DJ2
      ELSE
         D1 = DJ2
         D2 = DJ1
      END IF
      T = ( D1 + ( D2-DJ0 ) ) / CENDAY

*  Fractional part of JD(UT1), in seconds.
      F = DAYSEC * ( MOD(D1,1D0) + MOD(D2,1D0) )

*  GMST at this UT1.
      dau_GMST82 = iau_ANP ( DS2R * ( (A+(B+(C+D*T)*T)*T) + F ) )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2001
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. The Software is made available free of charge for use by:
*
*     a) private individuals for non-profit research; and
*
*     b) non-profit educational, academic and research institutions.
*
*  3. Commercial use of the Software is specifically excluded from the
*     terms and conditions of this license.  Commercial use of the
*     Software is subject to the prior written agreement of the Board on
*     terms to be agreed.
*
*  4. The provision of any version of the Software under the terms and
*     conditions specified herein does not imply that future versions
*     will also be made available under the same terms and conditions.
*
*  5. The user may modify the Software for his/her own purposes.  The
*     user may distribute the modified software provided that the Board
*     is informed and that a copy of the modified software is made
*     available to the Board on request.  All modifications made by the
*     user shall be clearly identified to show how the modified software
*     differs from the original Software, and the name(s) of the
*     affected routine(s) shall be changed.  The original SOFA Software
*     License text must be present.
*
*  6. In any published work produced by the user and which includes
*     results achieved by using the Software, the user shall acknowledge
*     that the Software was used in producing the information contained
*     in such publication.
*
*  7. The user may incorporate or embed the Software into other software
*     products which he/she may then give away free of charge but not
*     sell provided the user makes due acknowledgement of the use which
*     he/she has made of the Software in creating such software
*     products.  Any redistribution of the Software in this way shall be
*     made under the same terms and conditions under which the user
*     received it from the SOFA Center.
*
*  8. The user shall not cause the Software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  9. The Software is provided to the user "as is" and the Board makes
*     no warranty as to its use or performance.   The Board does not and
*     cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Board makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*
*-----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_ANP ( A )
*+
*  - - - - - - - -
*   i a u _ A N P
*  - - - - - - - -
*
*  Normalize angle into the range 0 <= A < 2pi.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     A          d       angle (radians)
*
*  Returned:
*     iau_ANP    d       angle in range 0-2pi
*
*  This revision:  2000 December 15
*
*  Copyright (C) 2001 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION A

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      DOUBLE PRECISION W

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      W = MOD(A,D2PI)
      IF ( W .LT. 0D0 ) W = W + D2PI
      iau_ANP = W

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2001
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. The Software is made available free of charge for use by:
*
*     a) private individuals for non-profit research; and
*
*     b) non-profit educational, academic and research institutions.
*
*  3. Commercial use of the Software is specifically excluded from the
*     terms and conditions of this license.  Commercial use of the
*     Software is subject to the prior written agreement of the Board on
*     terms to be agreed.
*
*  4. The provision of any version of the Software under the terms and
*     conditions specified herein does not imply that future versions
*     will also be made available under the same terms and conditions.
*
*  5. The user may modify the Software for his/her own purposes.  The
*     user may distribute the modified software provided that the Board
*     is informed and that a copy of the modified software is made
*     available to the Board on request.  All modifications made by the
*     user shall be clearly identified to show how the modified software
*     differs from the original Software, and the name(s) of the
*     affected routine(s) shall be changed.  The original SOFA Software
*     License text must be present.
*
*  6. In any published work produced by the user and which includes
*     results achieved by using the Software, the user shall acknowledge
*     that the Software was used in producing the information contained
*     in such publication.
*
*  7. The user may incorporate or embed the Software into other software
*     products which he/she may then give away free of charge but not
*     sell provided the user makes due acknowledgement of the use which
*     he/she has made of the Software in creating such software
*     products.  Any redistribution of the Software in this way shall be
*     made under the same terms and conditions under which the user
*     received it from the SOFA Center.
*
*  8. The user shall not cause the Software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  9. The Software is provided to the user "as is" and the Board makes
*     no warranty as to its use or performance.   The Board does not and
*     cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Board makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*
*-----------------------------------------------------------------------

  

      END







C
C
C     Subroutine ordem
C
C
C
C     Purpose
C
C       Orders data vectors in crescent value order.
C
C
C     Use
C
C     SUBROUTINE ORDEM (N,IOR,VAL)
C
C
C     Description of parameters
C
C       N      - number of points to be ordered
C       indx   - increasing order numbering of array "arr"
C       arr    - data array itself, NOT ORDERED
C
C
C     Subroutines and subprograms required
C
C
C
C     Comments
C
C
C

      SUBROUTINE ordem(n,indx,arr)
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
c       if(jstack.gt.NSTACK)pause 'NSTACK too small in ordem'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


