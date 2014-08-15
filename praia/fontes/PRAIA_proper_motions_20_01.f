c
c     PRAIA_proper_motions
c
c
c     Computes star proper motions from "xy" or "offset" files in PRAIA format
c     and updates the input files with the computed proper motions for
c     non-UCAC2 stars.
c
c
c     In order of priority, comes UCAC2 (original proper motions preserved),
c     then 2MASS-based, then USNOB1-based proper motions.
c
c
c     In this version, the proper motions in the input list and in the
c     2MASS and USNOB1 databases are considered not indexed. Thus,
c     the (RA,DEC) of stars in these three sets are crossed with each other
c     in order as to identify the common stars.
c
c
C     Last update: M. Assafin - 17/Oct/2010
c
c

      IMPLICIT REAL *8 (A-H,O-Z)

      parameter(stdin=5,stdout=6)


      dimension x(1000001),y(1000001),cseng(1000001),altu(1000001),
     ?fgcc(1000001),fumag(1000001),fumag2(1000001),cxmgu(1000001),
     ?codmg(1000001),codmg2(1000001),cxmgj(1000001),cxmgh(1000001),
     ?cxmgk(1000001),res2mg(1000001),resmg2(1000001),ermgj(1000001),
     ?ermgh(1000001),ermgk(1000001),copma(1000001),copmd(1000001),
     ?epma(1000001),epmd(1000001),coex(1000001),coey(1000001),
     ?cerau(1000001),cedeu(1000001),alfsic(1000001),delsic(1000001),
     ?nstaru(1000001),nfin(1000001),alsiuc(1000001),desiuc(1000001),
     ?ktir(1000001),oldra(1000001),oldde(1000001),kuth(1000001),
     ?kutm(1000001),zut(1000001),kutano(1000001),kutmes(1000001),
     ?kutdia(1000001),codj(1000001),iexps(1000001),ichfil(1000001),
     ?mfits(1000001),iobalv(1000001),nx(1000001),ny(1000001),
     ?numcom(1000001),egrxx(1000001),egryy(1000001),iflag(1000001),
     ?cflag(1000001),offra(1000001),offde(1000001),raxy(1000001)

      dimension pa(1000001),pd(1000001)

      dimension ra2ma(5000001),de2ma(5000001),dmgj(5000001),
     ?dmgh(5000001),dmgk(5000001),emgj(5000001),emgh(5000001),
     ?emgk(5000001),epom(5000001)

      dimension raub1(9000001),deub1(9000001),dmgb(9000001),
     ?dmgr(9000001),dmgi(9000001),epob(9000001)

      character*50 mraiz,braiz,list,catal,label
      character*100 mpcat

      character*20 ichfil,iobalv
      character*50 mfits
      character*4 cflag

c
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI


c
c     Auxiliary data
c

      idim=1000001
      idimm=5000001
      idimb=9000001

      icha=50


c
c     Reads 2MASS and USNOB1 paths from input data file
c

 1    format(a50)

      read (5,1) mraiz 
      read (5,1) braiz 
      read (5,*)

      write (*,*)
      write (*,*)
      write (*,*) 'PRAIA proper motion computations.'
      write (*,*) ' '
      write (*,*) ' '


c
c     Reads xy/offsets sets from input data file
c


 10   continue

      read (5,1,end=90) list
      read (5,*) key
      read (5,*) box
      read (5,1) label
      read (5,*)

      box=box/3600.d0
      box=box**2

c
c     Determines sky region for catalog star extraction
c

      write (*,*)
      write(*,*) 'Determining (RA,DEC) sky limits for catalogue extratio
     ?n ...'
      write (*,*)



      ramin=+1.d14
      demin=+1.d14
      ramax=-1.d14
      demax=-1.d14

      open (11,file=list)

      i=1

 15   read (11,1,end=35) catal

      open(1,file=catal)


 20   continue

      if (key.eq.1) then

      read(1,25,end=30) x(i),y(i),cseng(i),altu(i),fgcc(i),fumag(i),
     ?fumag2(i),cxmgu(i),codmg(i),codmg2(i),cxmgj(i),cxmgh(i),cxmgk(i),
     ?res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),copma(i),copmd(i),
     ?epma(i),epmd(i),coex(i),coey(i),cerau(i),cedeu(i),alfsic(i),
     ?delsic(i),nstaru(i),nfin(i),alsiuc(i),desiuc(i),ktir(i),oldra(i),
     ?oldde(i),kuth(i),kutm(i),zut(i),kutano(i),kutmes(i),kutdia(i),
     ?codj(i),iexps(i),ichfil(i),mfits(i),iobalv(i),nx(i),ny(i),
     ?numcom(i),egrxx(i),egryy(i),iflag(i),cflag(i)

 25   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,1x,
     ?f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?3(1x,i5),2(1x,f7.3),1x,i2,1x,a4)

      endif


      if (key.eq.2) then

      read(1,26,end=30) offra(i),offde(i),x(i),y(i),cseng(i),
     ?altu(i),fgcc(i),fumag(i),fumag2(i),cxmgu(i),codmg(i),codmg2(i),
     ?cxmgj(i),cxmgh(i),cxmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),
     ?ermgk(i),copma(i),copmd(i),epma(i),epmd(i),coex(i),coey(i),
     ?cerau(i),cedeu(i),alfsic(i),delsic(i),nstaru(i),nfin(i),alsiuc(i),
     ?desiuc(i),ktir(i),oldra(i),oldde(i),kuth(i),kutm(i),zut(i),
     ?kutano(i),kutmes(i),kutdia(i),codj(i),iexps(i),ichfil(i),mfits(i),
     ?iobalv(i),nx(i),ny(i),numcom(i),egrxx(i),egryy(i),iflag(i),
     ?cflag(i)

 26   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,
     ?1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,
     ?a20,3(1x,i5),2(1x,f7.3),1x,i2,1x,a4)

      endif

      if (oldra(i).lt.ramin) ramin=oldra(i)
      if (oldra(i).gt.ramax) ramax=oldra(i)

      if (oldde(i).lt.demin) demin=oldde(i)
      if (oldde(i).gt.demax) demax=oldde(i)

      go to 20

 30   close (1)

      go to 15

c

 35   close (11)


      ramax=ramax*15.d0
      ramin=ramin*15.d0


c
c     Extracts 2MASS catalog stars within the covered sky region
c

      call tmass (mraiz,ramin,ramax,demin,demax,ra2ma,de2ma,dmgj,dmgh,
     ?dmgk,emgj,emgh,emgk,epom,n2mass,idimm)


      write (*,*)
      write (*,37) n2mass
 37   format('Extracted 2MASS catalogue stars = ',i7)
      write (*,*)



c
c     Extracts USNOB1 catalog stars within the covered sky region
c

      call usnob1 (braiz,ramin,ramax,demin,demax,raub1,deub1,dmgb,dmgr,
     ?dmgi,epob,nusnob,idimb)

      write (*,*)
      write (*,38) nusnob
 38   format('Extracted USNOB1 catalogue stars = ',i7)
      write (*,*)



c
c     Computes proper motions for each xy/offset data file
c

      write (*,*)
      write (*,*) 'Performing proper motion computations.'
      write (*,*)


      open (11,file=list)


 40   catal=''

      read (11,1,end=80) catal

      write (*,1) catal

      open(1,file=catal)

      i=0

 45   i=i+1

      if (key.eq.1) then

      read(1,25,end=50) x(i),y(i),cseng(i),altu(i),fgcc(i),fumag(i),
     ?fumag2(i),cxmgu(i),codmg(i),codmg2(i),cxmgj(i),cxmgh(i),cxmgk(i),
     ?res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),copma(i),copmd(i),
     ?epma(i),epmd(i),coex(i),coey(i),cerau(i),cedeu(i),alfsic(i),
     ?delsic(i),nstaru(i),nfin(i),alsiuc(i),desiuc(i),ktir(i),oldra(i),
     ?oldde(i),kuth(i),kutm(i),zut(i),kutano(i),kutmes(i),kutdia(i),
     ?codj(i),iexps(i),ichfil(i),mfits(i),iobalv(i),nx(i),ny(i),
     ?numcom(i),egrxx(i),egryy(i),iflag(i),cflag(i)

      endif


      if (key.eq.2) then

      read(1,26,end=50) offra(i),offde(i),x(i),y(i),cseng(i),
     ?altu(i),fgcc(i),fumag(i),fumag2(i),cxmgu(i),codmg(i),codmg2(i),
     ?cxmgj(i),cxmgh(i),cxmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),
     ?ermgk(i),copma(i),copmd(i),epma(i),epmd(i),coex(i),coey(i),
     ?cerau(i),cedeu(i),alfsic(i),delsic(i),nstaru(i),nfin(i),alsiuc(i),
     ?desiuc(i),ktir(i),oldra(i),oldde(i),kuth(i),kutm(i),zut(i),
     ?kutano(i),kutmes(i),kutdia(i),codj(i),iexps(i),ichfil(i),mfits(i),
     ?iobalv(i),nx(i),ny(i),numcom(i),egrxx(i),egryy(i),iflag(i),
     ?cflag(i)

      endif

      raxy(i)=oldra(i)*15.d0

      pa(i)=copma(i)
      pd(i)=copmd(i)

      go to 45

 50   close (1)

      nest=i-1

c

      notm=0
      nob1=0

      do 70 i=1,nest


c
c     Skips UCAC2 stars (preserves UCAC2 proper motions)
c

      if (cxmgu(i).lt.25.d0) go to 70

c
c     2MASS-based proper motions
c

      if (cxmgj(i).gt.25.d0) go to 60


      do j=1,n2mass

      dx=(raxy(i)-ra2ma(j))*dcos(grarad*oldde(i))
      dy=oldde(i)-de2ma(j)

      d=dx**2+dy**2

      if (d.lt.box) go to 55

      enddo

      notm=notm+1

      go to 70

c

 55   dt=(codj(i)-epom(j))/365.25d0

      pa(i)=dx*3600.d0
      pa(i)=pa(i)/dt

      pd(i)=dy*3600.d0
      pd(i)=pd(i)/dt

      go to 70


c
c     USNOB1-based proper motions
c

 60   continue


      do j=1,nusnob

      dx=(raxy(i)-raub1(j))*dcos(grarad*oldde(i))
      dy=oldde(i)-deub1(j)

      d=dx**2+dy**2


      if (d.lt.box) go to 65

      enddo

      nob1=nob1+1

      go to 70

c

 65   dt=(codj(i)-epob(j))/365.25d0

      pa(i)=dx*3600.d0
      pa(i)=pa(i)/dt

      pd(i)=dy*3600.d0
      pd(i)=pd(i)/dt



 70   continue


c
c     Outputs updated xy/offset data file with new computed proper motions
c


      mpcat=''

      do jj=icha,1,-1
      if (label(jj:jj).ne.' ') go to 72
      enddo

 72   mpcat(1:jj)=label(1:jj)

      do kk=icha,1,-1
      if (catal(kk:kk).ne.' ') go to 73
      enddo

 73   mpcat(jj+1:jj+kk)=catal(1:kk)

c

      nucac2=0
      nmas=0
      nfield=0

      open(2,file=mpcat)
      

      do i=1,nest


      if (key.eq.1) then

      if (cxmgu(i).lt.25.d0) nucac2=nucac2+1

      if (cxmgu(i).gt.25.d0 .and. cxmgj(i).lt.25.d0) nmas=nmas+1

      if (cxmgu(i).gt.25.d0 .and. cxmgj(i).gt.25.d0) nfield=nfield+1


      write(2,25) x(i),y(i),cseng(i),altu(i),fgcc(i),fumag(i),
     ?fumag2(i),cxmgu(i),codmg(i),codmg2(i),cxmgj(i),cxmgh(i),cxmgk(i),
     ?res2mg(i),resmg2(i),ermgj(i),ermgh(i),ermgk(i),pa(i),pd(i),
     ?epma(i),epmd(i),coex(i),coey(i),cerau(i),cedeu(i),alfsic(i),
     ?delsic(i),nstaru(i),nfin(i),alsiuc(i),desiuc(i),ktir(i),oldra(i),
     ?oldde(i),kuth(i),kutm(i),zut(i),kutano(i),kutmes(i),kutdia(i),
     ?codj(i),iexps(i),ichfil(i),mfits(i),iobalv(i),nx(i),ny(i),
     ?numcom(i),egrxx(i),egryy(i),iflag(i),cflag(i)

      endif


      if (key.eq.2) then

      if (cxmgu(i).lt.25.d0) nucac2=nucac2+1

      if (cxmgu(i).gt.25.d0 .and. cxmgj(i).lt.25.d0) nmas=nmas+1

      if (cxmgu(i).gt.25.d0 .and. cxmgj(i).gt.25.d0) nfield=nfield+1


      write(2,26) offra(i),offde(i),x(i),y(i),cseng(i),
     ?altu(i),fgcc(i),fumag(i),fumag2(i),cxmgu(i),codmg(i),codmg2(i),
     ?cxmgj(i),cxmgh(i),cxmgk(i),res2mg(i),resmg2(i),ermgj(i),ermgh(i),
     ?ermgk(i),pa(i),pd(i),epma(i),epmd(i),coex(i),coey(i),
     ?cerau(i),cedeu(i),alfsic(i),delsic(i),nstaru(i),nfin(i),alsiuc(i),
     ?desiuc(i),ktir(i),oldra(i),oldde(i),kuth(i),kutm(i),zut(i),
     ?kutano(i),kutmes(i),kutdia(i),codj(i),iexps(i),ichfil(i),mfits(i),
     ?iobalv(i),nx(i),ny(i),numcom(i),egrxx(i),egryy(i),iflag(i),
     ?cflag(i)

      endif

      enddo

      close (2)

c


      go to 40

c

 80   close (11)


      go to 10

c

 90   continue

c


      write (*,91) nest
 91   format('Total number of stars in the field = ',i7)


      write (*,92) nucac2
 92   format('Total number of UCAC2 stars in the field = ',i7)


      write (*,93) nmas
 93   format('Total number of 2MASS stars in the field = ',i7)


      write (*,94) nfield
 94   format('Total number of field stars in the field = ',i7)


      per=100.d0*notm/nmas

      write (*,95) notm,per
 95   format('Number of 2MASS field stars not identified in the 2MASS ca
     ?talogue = ',i7,1x,'or ',f6.2,'%')


      per=100.d0*nob1/nfield

      write (*,96) nob1,per
 96   format('Number of field stars not identified in the USNOB1 catalog
     ?ue = ',i7,1x,'or ',f6.2,'%')


      write (*,*)



      write (*,*)
      write (*,*)
      write (*,*) 'Execution terminated successfully.'
      write (*,*) ' '
      write (*,*) ' '

      end




c
c
c     Subroutine TMASS
c
c
c     Extracts data from 2MASS catalog
c
c
c     - ra2ma,de2ma in degrees: (RA,DEC) 2MASS
c     - magnitudes: J, H and K.
c
c
c     Last updtate: M. Assafin  17/Oct/2010
c
c


      subroutine tmass (mraiz,ramin,ramax,demin,demax,ra2ma,de2ma,dmgj,
     ?dmgh,dmgk,emgj,emgh,emgk,epom,nest,idimm)

      implicit real*8 (a-h,o-z)


      INTEGER*4 CO1,CO2,JJD
      INTEGER*2 MAG1,MAG2,MAG3,ERCO
      INTEGER*1 ERMG1,ERMG2,ERMG3

      dimension ra2ma(5000001),de2ma(5000001),dmgj(5000001),
     ?dmgh(5000001),dmgk(5000001),emgj(5000001),emgh(5000001),
     ?emgk(5000001),epom(5000001)


      CHARACTER *1 ip,norsul(2),isig,menos,mais
      character *63 ifaixa
      character *9  iarq
      character *50 mraiz

      character*200 imaux


      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data iarq/'TMASS.ast'/
      data norsul/'m','p'/
 
c     DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
c     XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
c     YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c


      icha=63
      izone=1800


c
c     Checks RA=24hs region transiction 
c
c     keylim=0  -> not in transiction region
c
c     keylim=1  -> within transiction region
c

      aux=(demax+demin)/2.d0
      aux=grarad*aux
      aux=(ramax-ramin)*dcos(aux)
      aux=dabs(aux)

      keylim=0

      if (aux.gt.100.d0) keylim=1


c
c     Reads declination zones from 0.1 em 0.1 degrees from 2MASS
c

      dfaixa=demin-0.1d0
      decmax=demax

c

      write (*,*)
      write (*,*) '2MASS extraction zones:'
      write (*,*)
c

      nest=0

      do 30 k=1,izone

      dfaixa=dfaixa+0.1d0


      if (dfaixa-decmax.gt.0.1d0) go to 35

      j=dabs(dfaixa)*10.d0

    

      if (dfaixa.lt.0.d0) then
      j=dabs(dfaixa)*10.d0
      if (j.eq.-900) go to 30
      ip=norsul(1)
      else
      j=dabs(dfaixa)*10.d0
      if (j.eq.900) go to 35
      ip=norsul(2)
      endif

c
c     Mounts zone file name
c

      ifaixa=''
      ifaixa=mraiz

      do jj=icha,1,-1
      if (ifaixa(jj:jj).ne.' ') go to 1
      enddo
      
 1    jj=jj+1

      ifaixa(jj:jj)=ip

      write (ifaixa(jj+1:jj+3),'(i3.3)') j

      jj=jj+3

      ifaixa(jj+1:jj+9)=iarq(1:9)


      write (*,14) ifaixa
 14   format(a63)



c
c     Reads zones one at a time
c

      open (95,file=ifaixa,access='direct',form='unformatted',recl=23)

      n=0
 20   n=n+1

c
c     RA,DEC in degrees
c    
      read (95,rec=n,err=25) CO1,CO2,ERCO,MAG1,ERMG1,MAG2,
     ?ERMG2,MAG3,ERMG3,JJD

      RA=DBLE(CO1)/1.0D6
      DE=DBLE(CO2)/1.0D6
      ZMGJ=dble(MAG1)/1000.0d0
      EEMGJ=dble(ERMG1)/100.0d0
      ZMGH=dble(MAG2)/1000.0d0
      EEMGH=dble(ERMG2)/100.0d0
      ZMGK=dble(MAG3)/1000.0d0
      EEMGK=dble(ERMG3)/100.0d0
      DDJ=DBLE(JJD)/1.0D4+2451.0D3

c
c     Checks if star is within RA,DEC limits
c

      if (de.lt.demin) go to 20
      if (de.gt.demax) go to 20

      if (keylim.eq.0) then

      if (ra.gt.ramax) go to 20
      if (ra.lt.ramin) go to 20

      else

      if (ra.gt.ramin .and. ra.lt.ramax) go to 20

      endif


c

      nest=nest+1

      ra2ma(nest)=ra
      de2ma(nest)=de
      dmgj(nest)=zmgj
      dmgh(nest)=zmgh
      dmgk(nest)=zmgk
      emgj(nest)=eemgj
      emgh(nest)=eemgh
      emgk(nest)=eemgk
      epom(nest)=ddj



c
c     debug alfa delta
c

c     ra=ra/15.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG='-'
c     de=-de
c     ELSE
c     ISIG='+'
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0
c
c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,dmgj(nest),dmgh(nest),
c    ?dmgk(nest)
c
c
c     pause
c

      go to 20

 25   close (95)


 30   continue

 35   continue

      if (nest.gt.idimm) then
      write (*,*)
     ? 'Warning: number of 2MASS stars greater than supported.'
      endif

      return
      end





c
c
c     Subroutine USNOB1
c
c
c     Extracts data from USNOB1 catalog
c
c
c     - raub1,deub1 in degrees: (RA,DEC) USNOB1
c     - magnitudes: B, R and I
c
c
c     Last updtate: M. Assafin  17/Oct/2010
c
c


      subroutine usnob1 (braiz,ramin,ramax,demin,demax,raub1,deub1,dmgb,
     ?dmgr,dmgi,epob,nest,idimb)

      implicit real*8 (a-h,o-z)


      integer*4 x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,
     ?x14,x15,x16,x17,x18,x19,x20

      dimension raub1(9000001),deub1(9000001),dmgb(9000001),
     ?dmgr(9000001),dmgi(9000001),epob(9000001)


      character *1 isig
      character *63 ifaixa
      character *50 braiz
      character*10 value

 
c     DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
c     XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
c     YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      dj2000=2451544.5d0

      icha=50


c
c     Checks RA=24hs region transiction 
c
c     keylim=0  -> not in transiction region
c
c     keylim=1  -> within transiction region
c

      aux=(demax+demin)/2.d0
      aux=grarad*aux
      aux=(ramax-ramin)*dcos(aux)
      aux=dabs(aux)

      keylim=0

      if (aux.gt.100.d0) keylim=1


c
c     Reads declination zones from 0.1 em 0.1 degrees from USNOB1
c

      dfaixa=demin-0.1d0
      decmax=demax

c

      write (*,*)
      write (*,*) 'USNOB1 extraction zones:'
      write (*,*)

c

      nest=0

 1    dfaixa=dfaixa+0.1d0

      if (dfaixa-decmax.gt.0.1d0) go to 35


c
c     Mounts zone file name
c


      if (dfaixa.lt.0.d0) then

      aux=dabs(dfaixa)
      izone=aux
      aux=(aux-izone)*10.d0
      isub=aux+0.1d0
      isub=9-isub
      izone=-izone+89

      else

      aux=dfaixa
      izone=aux
      aux=(aux-izone)*10.d0
      isub=aux+0.1d0
      izone=izone+90

      endif

c


      ifaixa=''
      ifaixa=braiz

      do jj=icha,1,-1
      if (ifaixa(jj:jj).ne.' ') go to 2
      enddo
      
 2    jj=jj+1


      write (ifaixa(jj:jj+2),'(i3.3)') izone

      jj=jj+3

      ifaixa(jj:jj+1)='/b'

      jj=jj+2

      write (ifaixa(jj:jj+2),'(i3.3)') izone

      jj=jj+3

      write (ifaixa(jj:jj),'(i1.1)') isub 

      jj=jj+1

      ifaixa(jj:jj+3)='.cat'


      write (*,14) ifaixa
 14   format(a63)



c
c     Extracts USNOB1 data
c

      open (9,file=ifaixa,access='direct',form='unformatted',recl=80)


      n=0

 20   n=n+1


      read(9,rec=n,err=25) x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,
     ?x12,x13,x14,x15,x16,x17,x18,x19,x20



c
c     RA,DEC J2000, epoch J2000 in degrees
c    


      ra=(x01*1.d-2)/3600.d0
      de=(x02*1.d-2)/3600.d0-90.d0


c
c     Checks if star is within RA,DEC limits
c

      if (de.lt.demin) go to 20
      if (de.gt.demax) go to 20


      if (keylim.eq.0) then

      if (ra.gt.ramax) go to 20
      if (ra.lt.ramin) go to 20

      else

      if (ra.gt.ramin .and. ra.lt.ramax) go to 20

      endif



c
c     Extracts USNOB1 epoch (JD)
c

      write (value,'(i10.10)') x05
      read (value(02:04),*,err=20) epou
      epou=epou*0.1d0+1950.d0

      epo=365.25d0*(epou-2000.d0)+dj2000


c
c     Places USNOB1 (RA,DEC) at original average plate epochs
c


      write (value,'(i10.10)') x03
      read (value(07:10),*,err=20) upma

      read (value(03:06),*,err=20) upmd

      upma=upma*2.d-3
      upmd=upmd*2.d-3

      upma=upma-10.d0
      upmd=upmd-10.d0

      upma=upma/3600.d0
      upmd=upmd/3600.d0
 
      upma=upma/dabs(dcos(grarad*de))

      dt=epou-2000.d0

      de=de+upmd*dt

      ra=ra+upma*dt


c
c     Computes averaged B and R magnitudes from the 2 contributiong 
c     USNOB1 surveys; also extracts I magnitude.
c


      write (value,'(i10.10)') x06
      read (value(07:10),*,err=20) zmgb1
      write (value,'(i10.10)') x08
      read (value(07:10),*,err=20) zmgb2

      write (value,'(i10.10)') x07
      read (value(07:10),*,err=20) zmgr1
      write (value,'(i10.10)') x09
      read (value(07:10),*,err=20) zmgr2

      write (value,'(i10.10)') x10
      read (value(07:10),*,err=20) zmgi


      zmgb1=zmgb1*1.d-2
      zmgb2=zmgb2*1.d-2

      zmgr1=zmgr1*1.d-2
      zmgr2=zmgr2*1.d-2

      zmgi=zmgi*1.d-2


      if (zmgb1.gt.2 .and. zmgb2.gt.2) then
      zmgb=(zmgb1+zmgb2)/2.d0
      else
      if (zmgb1.lt.2) then
      zmgb=zmgb2
      else
      zmgb=zmgb1
      endif
      endif


      if (zmgr1.gt.2 .and. zmgr2.gt.2) then
      zmgr=(zmgr1+zmgr2)/2.d0
      else
      if (zmgr1.lt.2) then
      zmgr=zmgr2
      else
      zmgr=zmgr1
      endif
      endif


      if (zmgb.lt.2) zmgb=99.999d0
      if (zmgr.lt.2) zmgr=99.999d0
      if (zmgi.lt.2) zmgi=99.999d0


c

      nest=nest+1


c

      raub1(nest)=ra
      deub1(nest)=de
      dmgb(nest)=zmgb
      dmgr(nest)=zmgr
      dmgi(nest)=zmgi
      epob(nest)=epo



c
c     debug alfa delta
c
c 
c     ra=ra/15.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG='-'
c     de=-de
c     ELSE
c     ISIG='+'
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0
c
c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,dmgb(nest),dmgr(nest),
c    ?dmgi(nest),epo,epou,nest,n
c
c
c     pause
c
 

      go to 20

c

 25   close (9)


      go to 1

c

 35   continue

c

      if (nest.gt.idimb) then
      write (*,*)
     ? 'Warning: number of USNOB1 stars greater than supported.'
      endif

c

      return
      end

















c
c     Subroutine tmotion
c
c
c     Computes proper motions using another catalog (2MASS or USNOB1)
c     as 1rst epoch
c
c
c     key=1 2MASS
c     key=2 USNOB1
c
c     Last update: M. Assafin, 23/Oct/2009
c 
c
c

      subroutine pmotion (cat,id,pma,pmd,ncat,oldra,oldde,codj,nc,dj,
     ?box,key)

      IMPLICIT REAL *8 (A-H,O-Z)

      dimension id(1000001),pma(1000001),pmd(1000001),dj(1000001),
     ?codj(1000001),ra(1000001),de(1000001),oldra(1000001),
     ?oldde(1000001),ram(1000001),dem(1000001),ano(2)

      character*50 cat

c

c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      idim=1000001

c
c     initializing vectors
c

      ano(1)=1.d0
      ano(2)=365.25d0

c

      do i=1,idim
      id(i)=-1
      pma(i)=99.999d0
      pmd(i)=99.999d0
      enddo

c
c     Reads proper motion data base
c

      open(2,file=cat)

      ncat=0

 1    read(2,*,err=1,end=5) j,a2,a3,a4,a5,a6,a7,a8,a9,a10

      ncat=ncat+1

      ram(ncat)=a2
      dem(ncat)=a3

      dj(ncat)=a10-a7*ano(key)

      de(ncat)=a3-a9/36d5
      ra(ncat)=a2-a8/(54d6*dabs(dcos(grarad*de(ncat))))

      go to 1

 5    close (2)

c
c     Cross (RA,DEC) for identifying common stars
c

      do 200 j=1,ncat

      do 100 i=1,nc

      dx=54d3*(oldra(i)-ram(j))*dcos(grarad*oldde(i))

      if (dabs(dx).gt.box) go to 100

      dy=36d2*(oldde(i)-dem(j))

      if (dabs(dy).gt.box) go to 100

      dx=54d3*(oldra(i)-ra(j))*dcos(grarad*oldde(i))
      dy=36d2*(oldde(i)-de(j))


      dt=(codj(i)-dj(j))/365.25d0

      id(j)=i
      pma(j)=dx/dt
      pmd(j)=dy/dt

      go to 200

 100  continue

 200  continue


      return
      end




