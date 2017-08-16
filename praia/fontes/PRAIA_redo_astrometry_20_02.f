c
c     PRAIA_redo_astrometry
c
c   
c     Redoes the astrometry, from previous existing xy reduction files.
c
c     The re-reduction can be done from the (x,y) measurements, or from
c     the (RA,Dec) coordinates, using in this case the tangent plane
c     technique. In the case of direct re-reduction of mosaic files
c     (in the xy file format), the use of the tangent plane technique
c     option is mandatory.
c
c
c     The identification of the catalogue stars is straightforward, since
c     provisional (RA,DEC) are already available in the xy files. Only
c     the objects already identified in the images and archived in the xy
c     files are re-reduced. No access to the FITS files are made. No new
c     object identification is attempted.
c
c 
c
c      Last update: Marcelo Assafin - 24 Oct 2013
c   
c
c


      IMPLICIT REAL *8 (A-H,O-Z)
      parameter (stdin=5,idiobs=150000,ipmax=5001,icofsp=21,ng=10,
     ?idin=150,idin50=50,idin2=20000,iu2z=288,iu2i=240,iu4z=900,
     ?iu4i=1440,mpes=50000,nhist=30,nhisfw=10,jfdp=10000,kfund=150000)




      dimension contag(idiobs),histo(nhist),ico(nhist),ior(idiobs),
     ?nval(idiobs)


      dimension xob(idiobs),yob(idiobs),ilado(idiobs),seng(idiobs),
     ?iflag(idiobs),altu(idiobs),ialtu(idiobs)

      dimension xid(idiobs),yid(idiobs),idlado(idiobs),idx1(idiobs),
     ?idx2(idiobs),idy1(idiobs),idy2(idiobs),npix(idiobs)

      dimension exgcc(idiobs),eygcc(idiobs),sgcc(idiobs),fgcc(idiobs)

      dimension xold(idiobs),yold(idiobs),altold(idiobs),ialtol(idiobs)

      dimension ra2ma(idiobs),de2ma(idiobs),era2ma(idiobs),
     ?ede2ma(idiobs),dmgj(idiobs),dmgh(idiobs),dmgk(idiobs),
     ?emgj(idiobs),emgh(idiobs),emgk(idiobs),xra2ma(idiobs),
     ?yde2ma(idiobs),id2ma(idiobs),ddj2(idiobs),iduc2(idiobs),
     ?xrauc2(idiobs),ydeuc2(idiobs),iduc4(idiobs),xrauc4(idiobs),
     ?ydeuc4(idiobs)

      dimension rauc2(idiobs),deuc2(idiobs),erauc2(idiobs),
     ?edeuc2(idiobs),pmra(idiobs),pmde(idiobs),epmra(idiobs),
     ?epmde(idiobs),udmgj(idiobs),udmgh(idiobs),udmgk(idiobs),
     ?udmg(idiobs),cudmg(idiobs)


      dimension rauc4(idiobs),deuc4(idiobs),erauc4(idiobs),
     ?edeuc4(idiobs),pmra4(idiobs),pmde4(idiobs),epmra4(idiobs),
     ?epmde4(idiobs),udmgj4(idiobs),udmgh4(idiobs),udmgk4(idiobs),
     ?udmg4(idiobs),cudmg4(idiobs)


      dimension raucs(idiobs),deucs(idiobs),eraucs(idiobs),
     ?edeucs(idiobs),pmras(idiobs),pmdes(idiobs),epmras(idiobs),
     ?epmdes(idiobs),udmgjs(idiobs),udmghs(idiobs),udmgks(idiobs),
     ?udmgs(idiobs),cudmgs(idiobs),iducs(idiobs),itiras(idiobs),
     ?eras(idiobs),edes(idiobs),alfres(idiobs),delres(idiobs),
     ?coefxs(icofsp),coefys(icofsp),ecofxs(icofsp),ecofys(icofsp),
     ?xraucs(idiobs),ydeucs(idiobs)
   

      dimension era2(idiobs),ede2(idiobs),alfre2(idiobs),delre2(idiobs),
     ?coefx2(icofsp),coefy2(icofsp),ecofx2(icofsp),ecofy2(icofsp),
     ?itira2(idiobs)


      dimension coefxr(icofsp),coefyr(icofsp),ecofxr(icofsp),
     ?ecofyr(icofsp)


      dimension erau(idiobs),edeu(idiobs),alfreu(idiobs),delreu(idiobs),
     ?coefxu(icofsp),coefyu(icofsp),ecofxu(icofsp),ecofyu(icofsp),
     ?itirau(idiobs)

      dimension era4(idiobs),ede4(idiobs),alfre4(idiobs),delre4(idiobs),
     ?coefx4(icofsp),coefy4(icofsp),ecofx4(icofsp),ecofy4(icofsp),
     ?itira4(idiobs)

      dimension xmed(idiobs),ymed(idiobs),xpa(idiobs),ypa(idiobs),
     ?cocomx(icofsp),cocomy(icofsp),iuc2ma(idiobs),cudmg2(idiobs),
     ?iuc4ma(idiobs)

      dimension adx(jfdp),ady(jfdp),coordx(jfdp),coordy(jfdp)


      dimension inuu2(iu2z,iu2i),ir1u2(iu2z,iu2i),ir2u2(iu2z,iu2i)
      dimension inuu4(iu4z,iu4i),ir1u4(iu4z,iu4i),ir2u4(iu4z,iu4i)

 
      dimension xcir(idiobs),ycir(idiobs),lacir(idiobs)

      dimension xtra(idiobs),ytra(idiobs),xlatra(idiobs),ylatra(idiobs),
     ?angtra(idiobs)
   


      character*150 infits,imfits,names(idin2)

      character*200 ired2m,irmp2m,irme2m,ireduc,iredu4,iredus

      character*50 filexy,lista,centro,ialvos,ialvo2,ialvou,ialvo4,
     ?ialvus,lred2m,lrmp2m,lrme2m,lreduc,lredu4,lredus

      character*50 iilvo2,iiivo2,ifdp

      character *50 mraiz,uraiz,u4raiz,cpraia
      character *61 u2ind,u4ind

      character*50 subf

      character*1   menos,iver,nome(idin50),isig,ibrac
      character*69 iobold,ichobj,ihname(idin2),mchobj
      character*20 ichfil

      character*200 linha


      common/a14/ierro

      data ibrac/' '/
      data menos/'-'/
      data iobold/'12345678901234567890123456789012345678901234567890123
     ?4567890123456789'/

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))


c
c     Initial data 
c

      pi    = 0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi
c

      dj2000=2451545.d0

      tt0=2.8d-3

      izero=0
      zero=0.d0

      d99=99.999d0

      ireflex=1
c

      do i=1,mpes
      contag(i)=0.
      enddo

      do i=1,idiobs
      cudmg(i)=0.d0
      ior(i)=0
      enddo




c
c     Defines backspace (fortran compiler dependencies) 
c

      call backsp (1,nbac,91)

   

c
c     Lendo dados de entrada
c

c     open (1,file='PRAIA_redo_astrometry_20_02.dat')

      read (*,3) mraiz
      read (*,3) uraiz
      read (*,3) u4raiz
      read (*,3) cpraia
      read (*,*) iuserc

      read (*,*) krcat

      read (*,3) itec

      read (*,3) centro

      read (*,3) ialvos

      read (*,3) ifdp
      read (*,*) kfdp

      read (*,3) ialvo2
      read (*,3) ialvou
      read (*,3) ialvo4


      read (*,3) iilvo2

      read (*,3) iiivo2

      read (*,3) ialvus


      read (*,3) lreduc
      read (*,3) lredu4
      read (*,3) lred2m
      read (*,3) lrmp2m
      read (*,3) lrme2m
      read (*,3) lredus

      read (*,*) akey,idegx,iminux,asegx,idegy,iminuy,asegy

      read (*,*) ibr2ma

      read (*,*) box



      read (*,*) erpix
      read (*,*) pcort2
      read (*,*) pcortu
      read (*,*) pcorts
      read (*,*) ngrau
      read (*,*) ngrau3
      read (*,*) ngrau5
      read (*,*) inicio,iultmo
 3    format(a50)


c     close (1)


c
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'PRAIA - redo astrometry  setup')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
 
 
      rewind (5)


 2    continue

      read (*,5) linha
 5    format(a200)

      write (*,5) linha

      if (linha(1:1).ne.'*') go to 2


      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '


c

      

      open (1,file=centro)

      i=0
 4    read (1,*,end=2028)
      i=i+1
      go to 4

 2028 close (1)

      iul=i

      if (iul.eq.0) then
      write (*,*) ' '
      write (*,*) 'No xy files furnished. Exiting program.'
      write (*,*) ' '
      stop
      endif
      

      
      
      if (inicio.eq.0 .and. iultmo.eq.0) then

      inicio=1
      iultmo=iul

      else

      if (inicio.gt.iultmo) then
      write (*,*) ' '
      write (*,*) 'Initial/Final xy file ranks do not match. Exiting pro
     ?gram.'
      write (*,*) ' '
      stop
      endif


      if (inicio.gt.iul) then
      write (*,*) ' '
      write (*,*) 'Initial xy file outside list range. Exiting program.'
      write (*,*) ' '
      stop
      endif

      if (iultmo.gt.iul) then
      write (*,*) ' '
      write (*,*) 'Final xy file outside list range. Exiting program.'
      write (*,*) ' '
      stop
      endif

      endif




c
c     Stores UCAC2 acceleration index 
c

      u2ind=''
      u2ind=uraiz


      do l=1,idin50
      if (u2ind(l:l).eq.' ') go to 160
      enddo

 160  u2ind(l:l+10)='u2index.txt'


      open (3,file=u2ind)

      do i=1,10
      read (3,*)
      enddo

      do i=1,iu2z
      do j=1,iu2i
      read(3,*) nsbin,naz
      inuu2(i,j)=nsbin
      ir1u2(i,j)=naz-nsbin+1
      ir2u2(i,j)=ir1u2(i,j)+nsbin-1
      enddo
      enddo

      close (3)




c
c     Stores UCAC4 acceleration index 
c

      u4ind=''
      u4ind=u4raiz


      do l=1,idin50
      if (u4ind(l:l).eq.' ') go to 165
      enddo

 165  u4ind(l:l+10)='u4index.asc'


      open (3,file=u4ind)


      do i=1,iu4z
      do j=1,iu4i
      read(3,*) naz,nsbin
      inuu4(i,j)=nsbin
      ir1u4(i,j)=naz+1
      ir2u4(i,j)=ir1u4(i,j)+nsbin-1
      enddo
      enddo

      close (3)




c
c     Field Distortion Pattern (for all fields)
c
 
      nfdp=0

      open (23,file=ifdp,status='old',err=280)

      do i=1,jfdp
      read (23,*,end=275) adx(i),ady(i),coordx(i),coordy(i)
      enddo

 275  nfdp=i-1

 280  close (23)



c
c     Initiates re-reduction
c



      do i=1,icofsp
      coefxr(i)=0.d0
      coefyr(i)=0.d0
      ecofxr(i)=0.d0
      ecofyr(i)=0.d0
      coefx2(i)=0.d0
      coefy2(i)=0.d0
      ecofx2(i)=0.d0
      ecofy2(i)=0.d0
      enddo




      open (3,file=centro)


      do i=1,inicio-1
      read(3,*)
      enddo


c
c     Loop of xy files starts here
c


      do 60 lllll=inicio,iultmo



      read (3,3) filexy

      ramin=25.d0
      ramax=-1.d0

      demin=+91.d0
      demax=-91.d0

      rap=0.d0
      dep=0.d0


c
c     Reads data from xy files
c

      open (33,file=filexy)

      do 150 j=1,idiobs

      read (33,470,err=2010,end=151) xob(j),yob(j),seng(j),altu(j),
     ?fgcc(j),fumag,
     ?fumag2,xmgu,cudmg(j),cudmg2(j),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,era2(j),ede2(j),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny

 2010 continue
      
      xra2ma(j)=ra*15.d0
      yde2ma(j)=de

c     rap=(rap*(j-1)+ra)/j
c     dep=(dep*(j-1)+de)/j
      daj=(daj*(j-1)+dj)/j

      if (de.lt.demin) demin=de
      if (de.gt.demax) demax=de

      if (ra.lt.ramin) ramin=ra
      if (ra.gt.ramax) ramax=ra


 150  continue

 151  close (33)

      nest=j-1


c

      ired2m=''
      irmp2m=''
      ired2m=''
      irme2m=''
      ireduc=''
      iredu4=''
      iredus=''



c
c     Prepares output file names
c

      do kk=idin,1,-1
      if (infits(kk:kk).ne.' ') go to 2100
      enddo

 2100 continue

      ii=kk


      do kk=idin,1,-1
      if (infits(kk-4:kk).eq.'.fits'.or.infits(kk-4:kk).eq.'.FITS') 
     ?go to 2101
      enddo

      go to 2102

 2101 continue

      ii=kk-5

      go to 2105

 2102 continue

      do kk=idin,1,-1
      if (infits(kk-3:kk).eq.'.fts'.or.infits(kk-3:kk).eq.'.FTS') 
     ?go to 2103
      enddo

      go to 2105

 2103 continue

      ii=kk-4

c

 2105 continue



c
c     Root names
c

      irmp2m(1:ii)=infits(1:ii)
      ired2m(1:ii)=infits(1:ii) 
      irme2m(1:ii)=infits(1:ii)   
      ireduc(1:ii)=infits(1:ii)
      iredu4(1:ii)=infits(1:ii)
      iredus(1:ii)=infits(1:ii)

      irmp2m(ii+1:ii+1)='.'
      ired2m(ii+1:ii+1)='.'
      irme2m(ii+1:ii+1)='.'
      ireduc(ii+1:ii+1)='.'
      iredu4(ii+1:ii+1)='.'
      iredus(ii+1:ii+1)='.'

c
c     2MASS puro
c

      do iii=idin,1,-1
      if (lred2m(iii:iii).ne.' ') go to 2001
      enddo

 2001 continue

      ired2m(ii+2:ii+1+iii)=lred2m
      
c
c     2MASS tp
c

      do iii=idin,1,-1
      if (lrmp2m(iii:iii).ne.' ') go to 2002
      enddo

 2002 continue

      irmp2m(ii+2:ii+1+iii)=lrmp2m
      
c
c     2MASS tp+pm
c

      do iii=idin,1,-1
      if (lrme2m(iii:iii).ne.' ') go to 2003
      enddo

 2003 continue

      irme2m(ii+2:ii+1+iii)=lrme2m
      
c
c     UCAC2
c

      do iii=idin,1,-1
      if (lreduc(iii:iii).ne.' ') go to 2004
      enddo

 2004 continue

      ireduc(ii+2:ii+1+iii)=lreduc

      
c
c     UCAC4
c

      do iii=idin,1,-1
      if (lredu4(iii:iii).ne.' ') go to 2005
      enddo

 2005 continue

      iredu4(ii+2:ii+1+iii)=lredu4


      
c
c     User reference catalog
c

      do iii=idin,1,-1
      if (lredus(iii:iii).ne.' ') go to 2006
      enddo

 2006 continue

      iredus(ii+2:ii+1+iii)=lredus

c     

 91   format(a150)
 11   format(a50)
 12   format(50a1)
 92   format(150a1)
 93   format(a200)





c
c     write (*,91) infits
c     write (*,93) ired2m
c     write (*,93) irmp2m
c     write (*,93) irme2m
c     write (*,93) ireduc
c     write (*,93) iredu4
c     write (*,93) iredus
c     stop
c



c

      write (*,*)
      write (*,*)

      write (*,15) lllll,iultmo,infits
 15   format (1x,'Proccessing field ',i5,' of ',i5,': xy file = ', a50)

      write (*,*)

c
c     Estimating proccessing time
c

      if (lllll.ne.inicio) then

      tt=(iultmo-lllll+1)*(tt+tt0)

      tt=tt/24.d0

      iday=tt

      hour=(tt-iday)*24.d0
      ihour=hour
      minu=(hour-ihour)*60.d0
      seg=((hour-ihour)*60.d0-minu)*60.d0

      write (*,*)      
      write (*,16) iday,ihour,minu,seg
 16   format(1x,'Estimated time left for end of PRAIA re-reductions: ',
     ?i3,'days ',i2,'hs ',i2,'m ',f4.1,'s')  
      write (*,*)      

      endif

c
c     Initializing time routines
c


      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)

      if (lllll.ne.inicio) seg=tt0*3600.d0




c
c     Field Distortion Pattern correction 
c

      if (nfdp.eq.0) go to 120

      if (kfdp.gt.nfdp) kfdp=nfdp



      do i=1,nest

      do j=1,nfdp
      ialtu(j)=j
      nval(j)=(xob(i)-coordx(j))**2+(yob(i)-coordy(j))**2
      enddo
      
      call ordem (idiobs,nfdp,ialtu,nval)


      ad0=0.d0
      cadx=0.d0
      cady=0.d0

      do k=1,kfdp
      m=ialtu(k)
      di=dsqrt((xob(i)-coordx(m))**2+(yob(i)-coordy(m))**2)
      ad0=ad0+1.d0/di**2
      cadx=cadx+adx(m)/di**2
      cady=cady+ady(m)/di**2
      enddo

      cadx=cadx/ad0
      cady=cady/ad0

      xob(i)=xob(i)-cadx
      yob(i)=yob(i)-cady

      enddo


 120  continue


c
c     Field coordinate parameters
c


      epoj=2000D0+(daj-2451545D0)/365.25D0


      dec=(demax+demin)/2.d0
      ddec=(demax-demin)/2.d0

      deco=dabs(dcos(dec*grarad))

      if (ramax.gt.ramin) then
      rac=15.d0*(ramax+ramin)/2.d0
      drac=15.d0*(ramax-ramin)/deco
      else
      rac=15.d0*(ramax+ramin+24.d0)/2.d0
      drac=15.d0*(ramin+24.d0-ramax)/2.d0
      drac=drac/deco
      endif


c
c     Computes area for catalogue-catalogue tangent plane corrections
c     (2MASS x UCAC4)
c

      aremax=2.d0*drac
      aremay=2.d0*ddec

      if (akey.lt.1.d0) then

      areax=hmsgms(idegx,iminux,asegx)
      areay=hmsgms(idegy,iminuy,asegy)

      else

      areax=akey*drac
      areay=akey*ddec

      endif

      if (areax.gt.aremax) areax=aremax
      if (areay.gt.aremay) areay=aremay

      areax=grarad*areax/2.d0
      areay=grarad*areay/2.d0



c
c     Checks if one of the poles are in the FoV
c

      grac=grarad*rac
      gdec=grarad*dec

      bra=0.d0

      if (dec.lt.0.d0) then
      bde=-90.d0*grarad
      else      
      bde=+90.d0*grarad
      endif

      d=DEXY(bra,bde,grac,gdec)
      xx=XPAD(bra,bde,grac)/d
      yy=YPAD(bra,bde,grac,gdec)/d

c
c     Pole in FoV   -> krcat = 2: extracts catalogues using tangent plane projection
c     star by star (slower) 
c
c     Pole not in FoV -> krcat = 1:  extracts catalogues using (RA,DEC) accelerating
c     indexes (faster)  
c


      if (krcat.ne.2) then

      if (dabs(xx).le.areax .and. dabs(yy).le.areay) krcat=2

      endif


c
c     Retrieves 2MASS stars     
c

      write (*,*) 'Files searched for catalogue extraction:'
      write (*,*)

c
c
c     - ra2ma,de2ma in degrees (alfas e deltas do 2MASS)
c     - era2ma,ede2ma em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: J, H e K.
c
c

      if (krcat.eq.2) then

      call tmass (idiobs,mraiz,rac,dec,drac,ddec,ramin,ramax,demin,
     ?demax,ra2ma,de2ma,era2ma,ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,
     ?ddj2,n2mass)

      else

      call stmass (idiobs,mraiz,rac,dec,drac,ddec,ramin,ramax,demin,
     ?demax,ra2ma,de2ma,era2ma,ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,
     ?ddj2,n2mass)

      endif


c
c     Pega as estrelas UCAC2 da regiao em torno do CCD  
c


c
c     - ra,de em graus na epoca epoj da observacao ccd
c     - epram, epdem epoca Juliana da posicao media RA DEC no UCAC2
c     - era, ede erros de RA DE para epoca media do UCAC2 em segundos de arco
c     - epmRA, epmDe erros mov proprios em segundos de arco por ano
c     - erra, erde erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmra mov poprio (*dcosD) em segundos de arco por ano
c     - pmde mov proprio delta em segundos de arco por ano
c     - epmx, epmy erro dos mov proprios em RA DE em segundos de arco por ano
c     - zmg ... magnitudes 2MASS no UCAC2
c     - dmg  magnitude interna UCAC2 entre V e R
c


      if (krcat.eq.2) then


      call ucac2 (idiobs,uraiz,epoj,rac,dec,drac,ddec,ramin,ramax,demin,
     ?demax,rauc2,deuc2,erauc2,edeuc2,pmra,pmde,epmra,epmde,udmgj,udmgh,
     ?udmgk,udmg,nucac2)

      else

      call sucac2 (idiobs,iu2z,iu2i,inuu2,ir1u2,ir2u2,uraiz,epoj,rac,
     ?dec,drac,ddec,ramin,ramax,demin,demax,rauc2,deuc2,erauc2,edeuc2,
     ?pmra,pmde,epmra,epmde,udmgj,udmgh,udmgk,udmg,nucac2)

      endif


c
c     Pega as estrelas UCAC4 da regiao em torno do CCD  
c


c
c     - ra,de em graus na epoca epoj da observacao ccd
c     - epram4, epdem4 epoca Juliana da posicao media RA DEC no UCAC4
c     - epmra4, epmde4 erros mov proprios em segundos de arco por ano
c     - erauc4, edeuc4 erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmra4 mov poprio (*dcosD) em segundos de arco por ano
c     - pmde4 mov proprio delta em segundos de arco por ano
c     - epmra4, epmde4 erros mov proprios em RA DE em segundos de arco por ano
c     - udmg_4 ... magnitudes 2MASS no UCAC4
c     - udmg4  magnitude interna UCAC4 entre V e R
c



      if (krcat.eq.2) then

      call ucac4 (idiobs,u4raiz,epoj,rac,dec,drac,ddec,ramin,ramax,
     ?demin,demax,rauc4,deuc4,erauc4,edeuc4,pmra4,pmde4,epmra4,epmde4,
     ?udmgj4,udmgh4,udmgk4,udmg4,nucac4)


      else


      call sucac4 (idiobs,iu4z,iu4i,inuu4,ir1u4,ir2u4,u4raiz,epoj,rac,
     ?dec,drac,ddec,ramin,ramax,demin,demax,rauc4,deuc4,erauc4,edeuc4,
     ?pmra4,pmde4,epmra4,epmde4,udmgj4,udmgh4,udmgk4,udmg4,nucac4)


      endif



c
c     Pega as estrelas do catalogo de referencia (formato PRAIA) do usuario em
c     torno do CCD  
c
c
c     - raucs,deucs  (ra,dec) em graus na epoca epoj da observacao ccd
c     - eraucs, edeucs erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmras mov poprio (*dcosD) em segundos de arco por ano
c     - pmdes mov proprio delta em segundos de arco por ano
c     - epmras, epmdes erros mov proprios em RA DE em segundos de arco por ano
c     - udmg_s ... magnitudes 2MASS no catalogo do usuario (codigo 99)
c     - udmgs  magnitude da estrela no catalogo de referencia do usuario
c


      if (iuserc.eq.1) then

      call cuser (idiobs,cpraia,epoj,rac,dec,drac,ddec,ramin,ramax,
     ?demin,demax,raucs,deucs,eraucs,edeucs,pmras,pmdes,epmras,epmdes,
     ?udmgjs,udmghs,udmgks,udmgs,nucaus)

      endif




c
c     Feeds execution time routines
c


      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)

      tt0=tempop/3600.d0




c
c     Identifies 2MASS stars in the FoV
c


      do i=1,idiobs
      id2ma(i)=0
      enddo

c

      ico2ma=0

      errou=erpix

      do 403 i=1,n2mass
      do     j=1,nest

      ior(j)=j

      dx=dcos(de2ma(i)*grarad)*(xra2ma(j)-ra2ma(i))*3600d0
      dy=(yde2ma(j)-de2ma(i))*3600d0

      nval(j)=dsqrt(dx**2+dy**2)


      enddo

      call ordem (idiobs,nest,ior,nval)

      j=ior(1)


      dx=dcos(de2ma(i)*grarad)*(xra2ma(j)-ra2ma(i))*3600d0
      dy=(yde2ma(j)-de2ma(i))*3600d0


      if (dabs(dx).gt.errou) go to 403
      if (dabs(dy).gt.errou) go to 403

      id2ma(j)=i


      ico2ma=ico2ma+1



 403  continue


      write (*,*)
      write (*,*)
      write (*,*) '2MASS: extracted stars  = ',n2mass
      write (*,*) '2MASS: identified stars = ',ico2ma
      write (*,*)



c
c     Identifies UCAC2 stars in the FoV
c



      do i=1,idiobs
      iduc2(i)=0
      enddo

c

      icouc2=0

      errou=erpix

      do 405 i=1,nucac2
      do     j=1,nest

      ior(j)=j

      dx=dcos(deuc2(i)*grarad)*(xra2ma(j)-rauc2(i))*3600d0
      dy=(yde2ma(j)-deuc2(i))*3600d0

      nval(j)=dsqrt(dx**2+dy**2)


      enddo

      call ordem (idiobs,nest,ior,nval)

      j=ior(1)


      dx=dcos(deuc2(i)*grarad)*(xra2ma(j)-rauc2(i))*3600d0
      dy=(yde2ma(j)-deuc2(i))*3600d0


      if (dabs(dx).gt.errou) go to 405
      if (dabs(dy).gt.errou) go to 405

      iduc2(j)=i


      icouc2=icouc2+1



 405  continue


      write (*,*)
      write (*,*) 'UCAC2: extracted stars  = ',nucac2
      write (*,*) 'UCAC2: identified stars = ',icouc2
      write (*,*)



c
c     Identifies UCAC4 stars in the FoV
c


      do i=1,idiobs
      iduc4(i)=0
      enddo

c

      icouc4=0

      errou=erpix

      do 410 i=1,nucac4
      do     j=1,nest

      ior(j)=j

      dx=dcos(deuc4(i)*grarad)*(xra2ma(j)-rauc4(i))*3600d0
      dy=(yde2ma(j)-deuc4(i))*3600d0

      nval(j)=dsqrt(dx**2+dy**2)


      enddo

      call ordem (idiobs,nest,ior,nval)

      j=ior(1)


      dx=dcos(deuc4(i)*grarad)*(xra2ma(j)-rauc4(i))*3600d0
      dy=(yde2ma(j)-deuc4(i))*3600d0


      if (dabs(dx).gt.errou) go to 410 
      if (dabs(dy).gt.errou) go to 410

      iduc4(j)=i


      icouc4=icouc4+1



 410  continue


      write (*,*)
      write (*,*) 'UCAC4: extracted stars  = ',nucac4
      write (*,*) 'UCAC4: identified stars = ',icouc4
      write (*,*)


c
c     Identifies USER-catalogue stars in the FoV
c



      if (iuserc.eq.1) then


      do i=1,idiobs
      iducs(i)=0
      enddo

c

      icoucs=0

      errou=erpix

      do 415 i=1,nucaus
      do     j=1,nest

      ior(j)=j

      dx=dcos(deucs(i)*grarad)*(xra2ma(j)-raucs(i))*3600d0
      dy=(yde2ma(j)-deucs(i))*3600d0

      nval(j)=dsqrt(dx**2+dy**2)


      enddo

      call ordem (idiobs,nest,ior,nval)

      j=ior(1)


      dx=dcos(deucs(i)*grarad)*(xra2ma(j)-raucs(i))*3600d0
      dy=(yde2ma(j)-deucs(i))*3600d0


      if (dabs(dx).gt.errou) go to 415 
      if (dabs(dy).gt.errou) go to 415

      iducs(j)=i



      icoucs=icoucs+1



 415  continue


      write (*,*)
      write (*,*) 'USER: extracted stars  = ',nucaus
      write (*,*) 'USER: identified stars = ',icoucs
      write (*,*)


      endif




c
c     Tangent plane technique. Takes (X,Y) standard coordinates from (RA,DEC)
c     as if they were measured (x,y), and proceed with (RA,DEC) reduction with
c     respecto to the reference catalogues.
c
c     This reduction procedure is mandatory in the case that the xy files come
c     from a CCD mosaic reduction with PRAIA, were the stamped (x,y) are meaningless 
c     
c

      if (itec.eq.2) then

      grac=grarad*rac
      gdec=grarad*dec


      do i=1,nest

      bra=xra2ma(i)
      bde=yde2ma(i)

      d=DEXY(bra,bde,grac,gdec)
      xx=XPAD(bra,bde,grac)/d
      yy=YPAD(bra,bde,grac,gdec)/d

      xob(i)=xx*radgra*3600.d0
      yob(i)=yy*radgra*3600.d0

      enddo


      endif





c
c     (RA,DEC) reduction: 2MASS 
c


      do i=1,idiobs
      itira2(i)=0
      enddo


      call posred (idiobs,icofsp,ireflex,rac,dec,id2ma,n2mass,ra2ma,
     ?de2ma,nest,xob,yob,pcort2,ngrau,ngrau3,ngrau5,nstart,nfinal,
     ?xra2ma,yde2ma,era2,ede2,alfsi2,delsi2,alfre2,delre2,coefx2,coefy2,
     ?ecofx2,ecofy2,itira2,avam,dvam)

      ierro=0


c
c     (RA,DEC) reduction: UCAC2
c


      do i=1,idiobs
      itirau(i)=0
      enddo



      call posred (idiobs,icofsp,ireflex,rac,dec,iduc2,nucac2,rauc2,
     ?deuc2,nest,xob,yob,pcortu,ngrau,ngrau3,ngrau5,nstaru,nfinau,
     ?xrauc2,ydeuc2,erau,edeu,alfsiu,delsiu,alfreu,delreu,coefxu,coefyu,
     ?ecofxu,ecofyu,itirau,avamu,dvamu)


      ierro=0



c
c     (RA,DEC) reduction: UCAC4
c


      do i=1,idiobs
      itira4(i)=0
      enddo



      call posred (idiobs,icofsp,ireflex,rac,dec,iduc4,nucac4,rauc4,
     ?deuc4,nest,xob,yob,pcortu,ngrau,ngrau3,ngrau5,nstar4,nfina4,
     ?xrauc4,ydeuc4,era4,ede4,alfsi4,delsi4,alfre4,delre4,coefx4,coefy4,
     ?ecofx4,ecofy4,itira4,avam4,dvam4)


      ierro=0




c
c     (RA,DEC) reduction: USER reference catalogue 
c


      if (iuserc.eq.1) then

      do i=1,idiobs
      itiras(i)=0
      enddo



      call posred (idiobs,icofsp,ireflex,rac,dec,iducs,nucaus,raucs,
     ?deucs,nest,xob,yob,pcorts,ngrau,ngrau3,ngrau5,nstars,nfinas,
     ?xraucs,ydeucs,eras,edes,alfsis,delsis,alfres,delres,coefxs,coefys,
     ?ecofxs,ecofys,itiras,avams,dvams)


      ierro=0


      endif


c
c     Magnitude zero-point re-calibration of UCAC2
c


      zmu2=0.d0
      zsu2=0.d0

      do 710 i=1,nest

      if (iduc2(i).eq.0) go to 710
      j=iduc2(i)
      dif=udmg(j)-cudmg(i)
      zmu2=zmu2+dif
      zsu2=zsu2+dif**2

 710  continue

      call desvio (icouc2,zmu2,zsu2)



c
c     Magnitude zero-point re-calibration of UCAC4
c


      zmu4=0.d0
      zsu4=0.d0

      do 720 i=1,nest

      if (iduc4(i).eq.0) go to 720
      j=iduc4(i)
      dif=udmg4(j)-cudmg(i)
      zmu4=zmu4+dif
      zsu4=zsu4+dif**2

 720  continue

      call desvio (icouc4,zmu4,zsu4)



c
c     Magnitude zero-point re-calibration of USER catalogue
c


      zmus=0.d0
      zsus=0.d0

      do 730 i=1,nest

      if (iducs(i).eq.0) go to 730
      j=iducs(i)
      dif=udmgs(j)-cudmg(i)
      zmus=zmus+dif
      zsus=zsus+dif**2

 730  continue

      call desvio (icoucs,zmus,zsus)




c
c     Results of re-reductions (no tangent plane corrections)
c

c

      open (64,file=ired2m)
      open (65,file=ireduc)
      open (66,file=iredu4)

      if (iuserc.eq.1) then
      open (67,file=iredus)
      endif


      jj=0
      kk=0
      mm=0
      nn=0


      open (33,file=filexy)


      do i=1,nest


      read (33,470,err=2015) xo,yo,seng(i),altu(i),fgcc(i),fumag,
     ?fumag2,xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,aaaaa1,aaaaa2,bbbbb1,bbbbb2,
     ?iiiii1,iiiii2,cccc1,cccc2,kkkkk,r1,d1,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,imfits,mchobj,nx,ny

 2015 continue

      ktira=99
      ktirau=99
      ktira4=99
      ktiras=99

      alsi2=d99
      desi2=d99
      alsiu=d99
      desiu=d99
      alsi4=d99
      desi4=d99
      alsis=d99
      desis=d99

      if (id2ma(i).ne.0) then
      jj=jj+1
      j=id2ma(i)

c     write (*,*) '2mass jj j = ',jj,j

      ktira=itira2(jj)
      alsi2=alfre2(jj)
      desi2=delre2(jj)
      endif


      if (iduc2(i).ne.0) then
      kk=kk+1
      k=iduc2(i)

c     write (*,*) 'ucac2 kk k = ',kk,k

      ktirau=itirau(kk)
      alsiu=alfreu(kk)
      desiu=delreu(kk)
      endif



c
c     2MASS
c


      ra=xra2ma(i)/15.d0
      de=yde2ma(i)

      gumag=fumag+zmu4
      gumag2=fumag2+zmu4
      cudu4=cudmg(i)+zmu4
      cudu42=cudmg2(i)+zmu4
      ges2mg=zsu4


      write (64,470) xo,yo,seng(i),altu(i),fgcc(i),gumag,gumag2,
     ?xmgu,cudu4,cudu42,xmgj,xmgh,xmgk,ges2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era2(i),ede2(i),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny



c
c    UCAC2
c

      ra=xrauc2(i)/15.d0
      de=ydeuc2(i)

      gumag=fumag+zmu2
      gumag2=fumag2+zmu2
      cudu2=cudmg(i)+zmu2
      cudu22=cudmg2(i)+zmu2
      ges2mg=zsu2



      write (65,470) xo,yo,seng(i),altu(i),fgcc(i),gumag,gumag2,
     ?xmgu,cudu2,cudu22,xmgj,xmgh,xmgk,ges2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau(i),edeu(i),alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


 470  format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))


c
c    UCAC4
c


      if (iduc4(i).ne.0) then
      mm=mm+1
      m=iduc4(i)

c     write (*,*) 'ucac4 mm m = ',mm,m

      ktira4=itira4(mm)
      alsi4=alfre4(mm)
      desi4=delre4(mm)
      endif


      ra=xrauc4(i)/15.d0
      de=ydeuc4(i)

      gumag=fumag+zmu4
      gumag2=fumag2+zmu4
      cudu4=cudmg(i)+zmu4
      cudu42=cudmg2(i)+zmu4
      ges2mg=zsu4



      write (66,470) xo,yo,seng(i),altu(i),fgcc(i),gumag,gumag2,
     ?xmgu,cudu4,cudu42,xmgj,xmgh,xmgk,ges2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era4(i),ede4(i),alfsi4,delsi4,
     ?nstar4,nfina4,alsi4,desi4,ktira4,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny



c
c     USER catalogue
c


      if (iuserc.eq.1) then

      ra=xraucs(i)/15.d0
      de=ydeucs(i)


      if (iducs(i).ne.0) then
      nn=nn+1
      n=iducs(i)

c     write (*,*) 'user nn n = ',nn,n


      pmas=pmras(n)
      pmds=pmdes(n)
      epmas=epmras(n)
      epmds=epmdes(n)
      ktiras=itiras(nn)
      alsis=alfres(nn)
      desis=delres(nn)
      endif


      gumag=fumag+zmus
      gumag2=fumag2+zmus
      cudus=cudmg(i)+zmus
      cudus2=cudmg2(i)+zmus
      ges2mg=zsus



      write (67,470) xo,yo,seng(i),altu(i),fgcc(i),gumag,gumag2,
     ?xmgu,cudus,cudus2,xmgj,xmgh,xmgk,ges2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,eras(i),edes(i),alfsis,delsis,
     ?nstars,nfinas,alsis,desis,ktiras,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


      endif


      enddo

      close (33)
      close (64)
      close (65)
      close (66)
      close (67)


c
c     Target statistics (Observed minus Target) for reductions with the 2MASS,
c     UCAC2, UCAC4 and USER catalogues (no tangent plane corrections)
c



      write (*,*) 
      write (*,*) 
      write (*,*) '  **** UCAC2: Results **** '
      write (*,*) 

      call estat (box,ialvos,ireduc,ialvou)


      write (*,*) 
      write (*,*) 
      write (*,*) '  **** UCAC4: Results **** '
      write (*,*) 

      call estat (box,ialvos,iredu4,ialvo4)



      if (iuserc.eq.1) then


      write (*,*) 
      write (*,*) 
      write (*,*) '  **** User Catalogue: Results **** '
      write (*,*) 

      call estat (box,ialvos,iredus,ialvus)

      endif

      write (*,*) 
      write (*,*) 
      write (*,*) '  **** 2MASS Original: Results **** '
      write (*,*) 

      call estat (box,ialvos,ired2m,ialvo2)




c
c     Identifies common 2MASS/UCAC4 stars, in preparation for the
c     tangent plane correction of 2MASS w.r.t. the UCAC4.
c
 
      boxerr=erpix

      do i=1,idiobs
      iuc4ma(i)=0
      enddo


      do 650 i=1,nucac4
      do     j=1,n2mass

      ior(j)=j

      dx=dcos(deuc4(i)*grarad)*(ra2ma(j)-rauc4(i))*3600d0
      dy=(de2ma(j)-deuc4(i))*3600d0

      nval(j)=10*dsqrt(dx**2+dy**2)


      enddo

      call ordem (idiobs,n2mass,ior,nval)

      j=ior(1)


      dx=dcos(deuc4(i)*grarad)*(ra2ma(j)-rauc4(i))*3600d0
      dy=(de2ma(j)-deuc4(i))*3600d0

      xx=dsqrt(dx**2+dy**2)

      if (xx.gt.boxerr) go to 650

      iuc4ma(i)=j




 650  continue



c
c     Tangent plane correction of the 2MASS toward the UCAC4 catalogue
c     reference frame
c


      ucor=(pcortu/3600.d0)*grarad


c
c     First step: brings UCAC4 stars to the 2MASS epoch (individual epochs of
c     2MASS stars)
c


      do 600 k=1,nucac4

      if (iuc4ma(k).eq.0) go to 600
    
      j=iuc4ma(k)

      epo2m=2000D0+(ddj2(j)-2451545D0)/365.25D0

      dt=epo2m-epoj

      deuc4(k)=deuc4(k)+ pmde4(k)*dt/3600.d0
      rauc4(k)=rauc4(k)+(pmra4(k)*dt/3600.d0)/dabs(dcos(deuc4(k)*
     ?grarad))


 600  continue

c
c     Second step: projects 2MASS's and UCAC4's (RA,DEC)s into the tangent plane
c



      do k=1,idiobs
      itirau(k)=0
      xmed(k)=0.d0
      ymed(k)=0.d0
      xpa(k)=0.d0
      ypa(k)=0.d0
      enddo

      grac=grarad*rac
      gdec=grarad*dec

      ngpco=0
      ngraup=1


c
c     Gnomonic projection of (RA,DEC)s
c 



      ncom=0

      do 605 k=1,nucac4

      if (iuc4ma(k).eq.0) go to 605
    
      if (itira4(k).ne.0) go to 605

      j=iuc4ma(k)


      bra=grarad*ra2ma(j)
      bde=grarad*de2ma(j)

      d=DEXY(bra,bde,grac,gdec)

      xmedd=XPAD(bra,bde,grac)/d
      ymedd=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xmedd).gt.areax) go to 605
      if (dabs(ymedd).gt.areay) go to 605

      ncom=ncom+1

      xmed(ncom)=xmedd
      ymed(ncom)=ymedd

      bra=grarad*rauc4(k)
      bde=grarad*deuc4(k)

      d=DEXY(bra,bde,grac,gdec)
      xpa(ncom)=XPAD(bra,bde,grac)/d
      ypa(ncom)=YPAD(bra,bde,grac,gdec)/d

 605  continue


c
c     Third step: relates the projected (RA,DEC)s according to the 6 cte
c     model. This model corrects the following FoV residual deformations in the
c     2MASS (RA,DEC)s, with respect to the UCAC4 reference frame:
c
c     - (RA,DEC) origin offset
c     - (RA,DEC) axes' perpendicularity
c     - (RA,DEC) axes' rotation
c     - (RA,DEC) scales
c

      do i=1,icofsp
      cocomx(i)=0.d0
      cocomy(i)=0.d0
      enddo


      if (ncom.lt.3) then
      write (*,*) 'Less than 3 common 2MASS & UCAC4 stars.'
      write (*,*) 'Tangent plan correction not possible.'
      go to 62
      endif

c
c     Fitting the 6 cte model and finding the polynomial coefficients
c

      call isol (idiobs,icofsp,ngraup,ncom,xmed,ymed,xpa,ypa,cocomx,
     ?cocomy)

      if (ierro.eq.1) then
      ierro=0
      write (*,*) 'Tangent plan correction not possible.'
      go to 62
      endif



c
c     Applies the 6cte coefficients to the (RA,DEC) of all 2MASS stars
c     in the FOV
c



      do 620 i=1,nest

      if (id2ma(i).eq.0) go to 620
    
      j=id2ma(i)


      bra=grarad*ra2ma(j)
      bde=grarad*de2ma(j)

      d=DEXY(bra,bde,grac,gdec)
      xx2ma=XPAD(bra,bde,grac)/d
      yy2ma=YPAD(bra,bde,grac,gdec)/d

      x=pol(icofsp,xx2ma,yy2ma,cocomx,ngraup)
      y=pol(icofsp,xx2ma,yy2ma,cocomy,ngraup)

      xx=alff(x,y,grac,gdec)
      yy=deltt(xx,y,grac,gdec)

      ra2ma(j)=xx*radgra
      de2ma(j)=yy*radgra

 620  continue




c
c     Refines the transformation toward the UCAC4 for common 2MASS/UCAC4 stars,
c     adding the UCAC4 proper motions for these common 2MASS stars.
c


      do 630 i=1,nest

      if (id2ma(i).eq.0) go to 630
      if (iduc4(i).eq.0) go to 630
    
      j=id2ma(i)
      k=iduc4(i)


      epo2m=2000D0+(ddj2(j)-2451545D0)/365.25D0

      dt=epoj-epo2m

      de2ma(j)=de2ma(j)+ pmde4(k)*dt/3600.d0
      ra2ma(j)=ra2ma(j)+(pmra4(k)*dt/3600.d0)/dabs(dcos(de2ma(j)*
     ?grarad))



 630  continue




c
c     (RA,DEC) reduction: 2MASS + tangent plane correction + p.m. common
c     2MASS/UCAC4 stars. All 2MASS stars enter the reduction.
c



      do i=1,idiobs
      itira2(i)=0
      enddo


      call posred (idiobs,icofsp,ireflex,rac,dec,id2ma,n2mass,ra2ma,
     ?de2ma,nest,xob,yob,pcortu,ngrau,ngrau3,ngrau5,nstart,nfinal,
     ?xra2ma,yde2ma,era2,ede2,alfsi2,delsi2,alfre2,delre2,coefx2,coefy2,
     ?ecofx2,ecofy2,itira2,avam,dvam)

      ierro=0


c
c     Results: 2MASS + tangent plane correction + p.m. common
c     2MASS/UCAC4 stars.
c



      open (64,file=irmp2m)

      jj=0
      kk=0


      open(33,file=filexy)


      do i=1,nest

      read (33,470,err=2020) xo,yo,seng(i),altu(i),fgcc(i),fumag,
     ?fumag2,xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,aaaaa1,aaaaa2,bbbbb1,bbbbb2,
     ?iiiii1,iiiii2,cccc1,cccc2,kkkkk,r1,d1,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,imfits,mchobj,nx,ny

 2020 continue

      ktira=99


      alsi2=d99
      desi2=d99


      if (id2ma(i).ne.0) then
      jj=jj+1
      j=id2ma(i)

c     write (*,*) '2mass jj j = ',jj,j

      ktira=itira2(jj)
      alsi2=alfre2(jj)
      desi2=delre2(jj)
      endif


      if (iduc4(i).ne.0) then
      kk=kk+1
      k=iduc4(i)

c     write (*,*) 'ucac4 kk k = ',kk,k

      endif

      ra=xra2ma(i)/15.d0
      de=yde2ma(i)


      gumag=fumag+zmu4
      gumag2=fumag2+zmu4
      cudu4=cudmg(i)+zmu4
      cudu42=cudmg2(i)+zmu4
      ges2mg=zsmu4


      write (64,470) xo,yo,seng(i),altu(i),fgcc(i),gumag,gumag2,
     ?xmgu,cudu4,cudu42,xmgj,xmgh,xmgk,ges2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era2(i),ede2(i),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


      enddo

      close (33)
      close (64)



c
c     Target statistics (Observed minus Target) for reductions with the 2MASS,
c     with tangent plane correction + p.m. common 2MASS/UCAC4 stars. 
c




      write (*,*) 
      write (*,*) 
      write(*,*)'**** 2MASS + tangent plan correction + common pm  ****'
      write (*,*) 

      call estat (box,ialvos,irmp2m,iilvo2)


c
c     Here only common 2MASS/UCAC4 stars (with UCAC4 p.m.) are used, plus
c     some of the (mag_J) brightest 2MASS non-common UCAC4/2MASS stars. A fixed
c     p.m. is used for these extra stars. The fixed p.m. is computed by the p.m.
c     average of the faint end of UCAC4 stars. The underlying assumption is that
c     these 2MASS stars are "almost faint UCAC4 stars", thus they probably follow
c     the same average p.m. distribution of the faint end of the UCAC4 stars in
c     the FoV.  
c
c



      nmp=0

      do 625 k=1,nucac4


      if (iuc4ma(k).eq.0) go to 625
    
c     if (itira4(k).ne.0) go to 625

      bra=grarad*rauc4(k)
      bde=grarad*deuc4(k)

      d=DEXY(bra,bde,grac,gdec)
      xxpa=XPAD(bra,bde,grac)/d
      yypa=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xxpa).gt.areax) go to 625
      if (dabs(yypa).gt.areay) go to 625

      nmp=nmp+1

      ior(nmp)=k
      nval(nmp)=10000*udmg4(k)

 625  continue

      call ordem (idiobs,nmp,ior,nval)

      ncort=nmp/2.d0

      xmov=0.d0
      ymov=0.d0

      do i=1,ncort

      k=ior(nmp-i+1)
      xmov=xmov+pmra4(k)
      ymov=ymov+pmde4(k)

      enddo

      xmov=xmov/ncort
      ymov=ymov/ncort




c
c     Corrected 2MASS toward the UCAC4 using the tangent plane technique, plus
c     common 2MASS/UCAC4 p.m., plus applying the average UCAC4 p.m. of faint UCAC4
c     stars to some of the brightest non-common 2MASS/UCAC4 stars 
c


      ncom=0
      kk=0
      mm=0

      do 635 i=1,nest

      if (id2ma(i).eq.0) go to 635

      ncom=ncom+1

      if (iduc4(i).ne.0) then
      mm=mm+1
      itira2(ncom)=0
      go to 635
      endif

      itira2(ncom)=1

      kk=kk+1
      ior(kk)=ncom
      nval(kk)=10000*cudmg(i)

 635  continue


      if (kk.le.1) go to 637

      call ordem (idiobs,kk,ior,nval)

 637  nume=ibr2ma
      if (nume.gt.kk) nume=kk

      do i=1,nume

      k=ior(i)
      itira2(k)=2

      enddo


c

      ncom=0

      do 640 i=1,nest

      if (id2ma(i).eq.0) go to 640

      ncom=ncom+1

      if (itira2(ncom).ne.2) go to 640
      itira2(ncom)=0

      j=id2ma(i)

      epo2m=2000D0+(ddj2(j)-2451545D0)/365.25D0

      dt=epoj-epo2m

      de2ma(j)=de2ma(j)+ ymov*dt/3600.d0
      ra2ma(j)=ra2ma(j)+ (xmov*dt/3600.d0)/dabs(dcos(de2ma(j)*grarad))


 640  continue




c
c     (RA,DEC) reduction: 2MASS + tangent plane correction + common & average p.m.
c



      call posred (idiobs,icofsp,ireflex,rac,dec,id2ma,n2mass,ra2ma,
     ?de2ma,nest,xob,yob,pcortu,ngrau,ngrau3,ngrau5,nstart,nfinal,
     ?xra2ma,yde2ma,era2,ede2,alfsi2,delsi2,alfre2,delre2,coefx2,coefy2,
     ?ecofx2,ecofy2,itira2,avam,dvam)

      ierro=0



c
c     Results: 2MASS + tangent plane correction + common & average p.m.
c




      open (64,file=irme2m)

      jj=0
      kk=0

      open (33,file=filexy)

      do i=1,nest

      read (33,470,err=2025) xo,yo,seng(i),altu(i),fgcc(i),fumag,
     ?fumag2,xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,aaaaa1,aaaaa2,bbbbb1,bbbbb2,
     ?iiiii1,iiiii2,cccc1,cccc2,kkkkk,r1,d1,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,imfits,mchobj,nx,ny

 2025 continue


      ktira=99


      alsi2=d99
      desi2=d99


      if (id2ma(i).ne.0) then
      jj=jj+1
      j=id2ma(i)

c     write (*,*) '2mass jj j = ',jj,j

      ktira=itira2(jj)
      alsi2=alfre2(jj)
      desi2=delre2(jj)
      endif


      if (iduc4(i).ne.0) then
      kk=kk+1
      k=iduc4(i)

c     write (*,*) 'ucac4 kk k = ',kk,k

      endif

      ra=xra2ma(i)/15.d0
      de=yde2ma(i)

      gumag=fumag+zmu4
      gumag2=fumag2+zmu4
      cudu4=cudmg(i)+zmu4
      cudu42=cudmg2(i)+zmu4
      ges2mg=zsmu4


      write (64,470) xo,yo,seng(i),altu(i),fgcc(i),gumag,gumag2,
     ?xmgu,cudu4,cudu42,xmgj,xmgh,xmgk,ges2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era2(i),ede2(i),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


      enddo

      close (33)
      close (64)


c
c     Target statistics (Observed minus Target) for reductions with the 2MASS,
c     with tangent plane correction + common & average p.m. 
c



      write (*,*) 
      write (*,*) 
      write(*,*)'*** 2MASS + tangent plan + common & non-common pm  ***'
      write (*,*) 


      call estat (box,ialvos,irme2m,iiivo2)

c
c     Ending program
c


 62   continue


 60   continue

c

      close (3)

c


      write (*,*)
      write (*,61) 
 61   format (23x,'Processing terminated.')
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '

      end



C
C     SUBROUTINE MATINV
C
C     PURPOSE
C       INVERT A SYMMETRIC MATRIX AND CALCULATE ITS DETERMINANT
C
C     USAGE
C       CALL MATINV ( NORDER, DET)
C
C     DESCRIPTION OF PARAMETERS
C       ARRAY  - INPUT MATRIX WICH IS REPLACED BY ITS INVERSE
C       NORDER - DEGREE OF MATRIX (ORDER OF DETERMINANT)
C       DET    - DETERMINANT OF INPUT MATRIX
C
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C       NONE
C
C     MODIFICATIONS FOR FORTRAN II
C       OMIT DOUBLE PRECISION SPECIFICATIONS
C       CHANGE DABS TO ABSF IN STATEMENT 23
C
C     COMMENTS
C       DIMENSION STATEMANTS VALID FOR NORDER UP TO 10
C
      SUBROUTINE MATINV (NORDER, icofsp, array, DET)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION ARRAY (icofsp,icofsp), IK(icofsp), JK(icofsp)
      COMMON /A14/IERRO
C
   10 DET = 1.D0
   11 DO 100 K=1, NORDER
C
C        FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
C
      AMAX= 0.D0
   21 DO 30 I=K, NORDER
      DO 30 J=K, NORDER
   23 IF (DABS(AMAX) - DABS(ARRAY(I,J))) 24, 24, 30
   24 AMAX = ARRAY(I,J)
      IK(K) = I
      JK(K) = J
   30 CONTINUE
C
C        INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)
C
   31 IF (AMAX) 41, 32, 41
   32 DET = 0.D0
      GO TO 140
   41 I = IK(K)
      IF (I-K) 21, 51, 43
   43 DO 50 J=1, NORDER
      SAVE = ARRAY(K,J)
      ARRAY(K,J) = ARRAY(I,J)
   50 ARRAY(I,J) = -SAVE
   51 J = JK(K)
      IF (J-K) 21, 61, 53
   53 DO 60 I=1, NORDER
      SAVE = ARRAY(I,K)
      ARRAY (I,K) = ARRAY(I,J)
   60 ARRAY (I,J) = -SAVE
C
C        ACCUMULATE ELEMENTS OF INVERSE MATRIX
C
   61 DO 70 I=1, NORDER
      IF (I-K) 63, 70, 63
   63 IF (AMAX.EQ.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
      ARRAY(I,K) = -ARRAY(I,K) / AMAX
   70 CONTINUE
   71 DO 80 I=1, NORDER
      DO 80 J=1, NORDER
      IF (I-K) 74, 80, 74
   74 IF (J-K) 75, 80, 75
   75 ARRAY(I,J) = ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
   80 CONTINUE
   81 DO 90 J=1, NORDER
      IF (J-K) 83, 90, 83
   83 IF (AMAX.EQ.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
      ARRAY(K,J) = ARRAY(K,J) / AMAX
   90 CONTINUE
      IF (AMAX.EQ.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
      ARRAY(K,K) = 1.D0 / AMAX
  100 DET = DET * AMAX
C
C        RESTORE ORDERING OF MATRIX
C
  101 DO 130 L=1, NORDER
      K = NORDER - L + 1
      J = IK(K)
      IF (J-K) 111, 111, 105
  105 DO 110 I=1, NORDER
      SAVE = ARRAY(I,K)
      ARRAY(I,K) = -ARRAY(I,J)
  110 ARRAY(I,J) = SAVE
  111 I = JK(K)
      IF (I-K) 130, 130, 113
  113 DO 120 J=1, NORDER
      SAVE = ARRAY(K,J)
      ARRAY(K,J) = -ARRAY(I,J)
  120 ARRAY(I,J) = SAVE
  130 CONTINUE
  140 CONTINUE
      RETURN
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
C     SUBROUTINE ORDEM (IDIM,N,IORDEM,NVAL)
C
C
C     Description of parameters
C
C       IDIM   - vector dimension
C	N      - number of points to be ordered
C       IORDEM - increasing order numbering of NVAL
C       NVAL   - data vector to be ordered
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
      SUBROUTINE ORDEM (idiobs,N,IORDEM,NVAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IORDEM(idiobs),NVAL(idiobs)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=NVAL(L)
          IRA=IORDEM(L)
        ELSE
          RRA=NVAL(IR)
          IRA=IORDEM(IR)
          NVAL(IR)=NVAL(1)
          IORDEM(IR)=IORDEM(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            NVAL(1)=RRA
            IORDEM(1)=IRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(NVAL(J).LT.NVAL(J+1))J=J+1
          ENDIF
          IF(RRA.LT.NVAL(J))THEN
            NVAL(I)=NVAL(J)
            IORDEM(I)=IORDEM(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        NVAL(I)=RRA
        IORDEM(I)=IRA
      GO TO 10
      END





c
c
c     Subrotina tmass
c
c
c     Pega as estrelas do catalogo 2MASS astrometrico. Extracao via projecao
c     gnomonica, estrela a estrela.
c
c
c     - ra2ma,de2ma em graus (alfas e deltas do 2MASS)
c     - era2ma,ede2ma em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: J, H e K.
c
c
c     Usa-se esta rotina qdo o polo cai dentro do campo, quando entao eh
c     obrigatorio o uso da projecao gnomonica, estrela a estrela.
c
c     Mod.: M. Assafin 02/Jul/2013
c
c

      subroutine tmass (idiobs,mraiz,rac,dec,drac,ddec,ramin,ramax,
     ?demin,demax,ra2ma,de2ma,era2ma,ede2ma,dmgj,dmgh,dmgk,emgj,emgh,
     ?emgk,ddj2,nest)

      implicit real*8 (a-h,o-z)


      INTEGER*4 CO1,CO2,JJD
      INTEGER*2 MAG1,MAG2,MAG3,ERCO
      INTEGER*1 ERMG1,ERMG2,ERMG3


      dimension ra2ma(idiobs),de2ma(idiobs),era2ma(idiobs),
     ?ede2ma(idiobs),dmgj(idiobs),dmgh(idiobs),dmgk(idiobs),
     ?emgj(idiobs),emgh(idiobs),emgk(idiobs),ddj2(idiobs)


      CHARACTER *1 MAIS,MENOS,ip,norsul(2),isig
      character *63 ifaixa
      character *9  iarq
      character *50 mraiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data iarq/'TMASS.ast'/
      data norsul/'m','p'/
 
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

C
C     Initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c


      idin50=50
      izon=1800

c

      dmar=10.d0/3600.d0

      ra1=ramin*15.d0-dmar
      ra2=ramax*15.d0+dmar

      key=1

      if (ra2-ra1.gt.180.d0) key=2

      de1=demin-dmar
      de2=demax+dmar

c
c     Leitura das faixas de declinacao de 0.1 em 0.1 graus
c     do 2mass astrometrico
c

      dfaixa=de1-0.2d0
      famax=de2+0.2d0


c


      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.1d0


      if (dfaixa-famax.gt.0.1d0) go to 35


      if (dfaixa.lt.0.d0) then
      j=dabs(dfaixa)*10.d0
      if (j.eq.900) go to 30
      ip=norsul(1)
      else
      j=dabs(dfaixa)*10.d0
      if (j.eq.900) go to 35
      ip=norsul(2)
      endif

c
c     Monta nome do arquivo da faixa, na leitura de 0.1 em 0.1 graus
c


      ifaixa=''
      ifaixa=mraiz



      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+12),'(a1,i3.3,a9)') ip,j,iarq



      write (*,14) ifaixa
 14   format(1x,'2MASS ',a63)


c
c     Le a faixa de 0.1 da vez
c

      open (95,file=ifaixa,access="direct",form="unformatted",recl=23)

      n=0
 20   n=n+1

c
c     RA,DEC em graus
c     erro em RA,DEC em segundos de arco
c    
      read (95,rec=n,err=25) CO1,CO2,ERCO,MAG1,ERMG1,MAG2,
     ?ERMG2,MAG3,ERMG3,JJD

      DE=DBLE(CO2)/1.0D6

c
c     Checa se estrela cai dentro do campo
c

      if (de.lt.de1) go to 20
      if (de.gt.de2) go to 20


      RA=DBLE(CO1)/1.0D6


      if (key.eq.1) then

      if (ra.lt.ra1) go to 20
      if (ra.gt.ra2) go to 20
      
      else


      if (ra.gt.ra1 .and. ra.lt.ra2) go to 20

      endif

c


      ERDE=dble(MOD(ERCO,100))/100.0d0
      ERRA=(dble(ERCO)-ERDE)/10000.0d0
      ZMGJ=dble(MAG1)/1000.0d0
      EEMGJ=dble(ERMG1)/100.0d0
      ZMGH=dble(MAG2)/1000.0d0
      EEMGH=dble(ERMG2)/100.0d0
      ZMGK=dble(MAG3)/1000.0d0
      EEMGK=dble(ERMG3)/100.0d0
      DDJ=DBLE(JJD)/1.0D4+2451.0D3



      nest=nest+1

      ra2ma(nest)=ra
      de2ma(nest)=de
      era2ma(nest)=erra
      ede2ma(nest)=erde
      dmgj(nest)=zmgj
      dmgh(nest)=zmgh
      dmgk(nest)=zmgk
      emgj(nest)=eemgj
      emgh(nest)=eemgh
      emgk(nest)=eemgk
      ddj2(nest)=ddj


c
c     debug alfa delta
c

c     ra=ra/15.d0
c     dmag=mag/100.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,dmgj(nest),dmgh(nest),
c    ?dmgk(nest)




      go to 20

 25   close (95)


 30   continue

 35   continue

 
      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of 2MASS stars.' 
      write (*,*)
      write (*,*)
      endif


      return
      end








c
c
c     Subrotina stmass
c
c
c     Pega as estrelas do catalogo 2MASS astrometrico, usando indexacao de catalogo.
c
c
c     - ra2ma,de2ma em graus (alfas e deltas do 2MASS)
c     - era2ma,ede2ma em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: J, H e K.
c
c
c     Nesta versao, a leitura eh feita mais rapida, sem projecao gnomonica, 
c     usando a indexacao do catalogo. Usar a rotina tmass qdo o polo cai dentro
c     do campo, quando entao eh obrigatorio o uso da projecao gnomonica, estrela
c     a estrela.
c
c
c     Mod.: M. Assafin  03/Jul/2013
c
c

      subroutine stmass (idiobs,mraiz,rac,dec,drac,ddec,ramin,ramax,
     ?demin,demax,ra2ma,de2ma,era2ma,ede2ma,dmgj,dmgh,dmgk,emgj,emgh,
     ?emgk,ddj2,nest)

      implicit real*8 (a-h,o-z)


      INTEGER*4 CO1,CO2,JJD
      INTEGER*2 MAG1,MAG2,MAG3,ERCO
      INTEGER*1 ERMG1,ERMG2,ERMG3


      dimension ra2ma(idiobs),de2ma(idiobs),era2ma(idiobs),
     ?ede2ma(idiobs),dmgj(idiobs),dmgh(idiobs),dmgk(idiobs),
     ?emgj(idiobs),emgh(idiobs),emgk(idiobs),ddj2(idiobs)

      dimension inic(360),ifim(360),inum(360),jxmin(2),jxmax(2),
     ?cxmin(2),cxmax(2)


      CHARACTER *1 MAIS,MENOS,ip,norsul(2),isig
      character *63 ifaixa,kfaixa
      character *9  iarq,karq
      character *50 mraiz


      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data iarq/'TMASS.ast'/
      data karq/'TMASS.acc'/
      data norsul/'m','p'/
 
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      idin50=50
      izon=1800

c

      dmar=10.d0/3600.d0

      ra1=ramin*15.d0-dmar
      ra2=ramax*15.d0+dmar

      key=1

      if (ra2-ra1.gt.180.d0) key=2

      de1=demin-dmar
      de2=demax+dmar



c
c     Initializing vectors
c

      do i=1,2
      jxmin(i)=0
      jxmax(i)=0
      cxmin(i)=0
      cxmax(i)=0
      enddo


c
c     Leitura das faixas de declinacao de 0.1 em 0.1 graus
c     do 2mass astrometrico
c

      dfaixa=de1-0.2d0
      famax=de2+0.2d0


c


      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.1d0


      if (dfaixa-famax.gt.0.1d0) go to 35


      if (dfaixa.lt.0.d0) then
      j=dabs(dfaixa)*10.d0
      if (j.eq.900) go to 30
      ip=norsul(1)
      else
      j=dabs(dfaixa)*10.d0
      if (j.eq.900) go to 35
      ip=norsul(2)
      endif

      jaca=j

c
c     Monta nome do arquivo da faixa, na leitura de 0.1 em 0.1 graus
c


      ifaixa=''
      ifaixa=mraiz



      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+12),'(a1,i3.3,a9)') ip,j,iarq



      write (*,14) ifaixa
 14   format(1x,'2MASS ',a63)


c
c     Monta nome do arquivo indexador da faixa
c

      kfaixa=ifaixa

      do j=1,idin50+13
      if (ifaixa(j:j).eq.' ') go to 15
      enddo

 15   kfaixa(j-3:j-1)='acc'



c
c     Carrega indexador de estrelas da faixa
c
c     Indice alfa de 1 em 1 grau
c
c     inic - rec para estrela no inicio da faixa de alfa
c     ifim - rec para estrela no final  da faixa de alfa
c     inu  - numero de estrelas na faixa de alfa de alfa (pode ser zero)
c
c

      open (95,file=kfaixa)

      do j=1,360
      read (95,16) inic(j),ifim(j),inum(j)
 16   format(16x,3i7)
      enddo

      close (95)


c
c     Calcula range minimo e maximo de alfa para extracao da faixa da vez
c     (limites alfa indiferentes do hemisferio)
c



      if (key.eq.1) then

      inde=1

      jxmin(1)=ra1+1
      jxmax(1)=ra2+1

      cxmin(1)=ra1
      cxmax(1)=ra2


      else

      inde=2

      jxmin(1)=ra2+1
      jxmax(1)=360

      jxmin(2)=1
      jxmax(2)=ra1+1

      cxmin(1)=ra2
      cxmax(1)=360.d0
      cxmin(2)=0.d0
      cxmax(2)=ra1


      endif



c
c     Le a faixa de 0.1 da vez
c

      open (95,file=ifaixa,access="direct",form="unformatted",recl=23)


c


      do 20 in=1,inde

      do 19 nn=jxmin(in),jxmax(in)

      if(inum(nn).eq.0) go to 19

      i1=inic(nn)
      i2=ifim(nn)


      do 18 n=i1,i2



c
c     RA,DEC em graus
c     erro em RA,DEC em segundos de arco
c

    
      read (95,rec=n) CO1,CO2,ERCO,MAG1,ERMG1,MAG2,
     ?ERMG2,MAG3,ERMG3,JJD


      DE=DBLE(CO2)/1.0D6

      if (de.lt.de1) go to 18
      if (de.gt.de2) go to 18

c

      RA=DBLE(CO1)/1.0D6

c


c     if (key.eq.1) then
c
c     if (ra.lt.ra1) go to 18
c     if (ra.gt.ra2) go to 18
c     
c     else
c
c
c     if (ra.gt.ra1 .and. ra.lt.ra2) go to 18
c
c     endif



      if (ra.lt.cxmin(in)) go to 18
      if (ra.gt.cxmax(in)) go to 18

c



      ERDE=dble(MOD(ERCO,100))/100.0d0
      ERRA=(dble(ERCO)-ERDE)/10000.0d0
      ZMGJ=dble(MAG1)/1000.0d0
      EEMGJ=dble(ERMG1)/100.0d0
      ZMGH=dble(MAG2)/1000.0d0
      EEMGH=dble(ERMG2)/100.0d0
      ZMGK=dble(MAG3)/1000.0d0
      EEMGK=dble(ERMG3)/100.0d0
      DDJ=DBLE(JJD)/1.0D4+2451.0D3


      nest=nest+1

      ra2ma(nest)=ra
      de2ma(nest)=de
      era2ma(nest)=erra
      ede2ma(nest)=erde
      dmgj(nest)=zmgj
      dmgh(nest)=zmgh
      dmgk(nest)=zmgk
      emgj(nest)=eemgj
      emgh(nest)=eemgh
      emgk(nest)=eemgk
      ddj2(nest)=ddj



c
c     debug alfa delta
c

c     ra=ra/15.d0
c     dmag=mag/100.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,dmgj(nest),dmgh(nest),
c    ?dmgk(nest)


 18   continue

 19   continue

 20   continue


c

      close (95)

c

 30   continue

 35   continue

 
      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of 2MASS stars.' 
      write (*,*)
      write (*,*)
      endif


      return
      end





c
c
c     Subrotina ucac2
c
c
c     Pega as estrelas do catalogo UCAC2.
c
c
c     - rauc2,deuc2 em graus (alfas e deltas do UCAC2)
c     - erauc2,edeuc2 em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: ucac2, J, H e K.
c
c     Usa-se esta rotina qdo o polo cai dentro do campo, quando entao eh
c     obrigatorio o uso da projecao gnomonica, estrela a estrela.
c
c     Mod. 03/Jul/2013
c
c

      subroutine ucac2 (idiobs,uraiz,epoj,rac,dec,drac,ddec,ramin,ramax,
     ?demin,demax,rauc2,deuc2,erauc2,edeuc2,pmra,pmde,epmra,epmde,udmgj,
     ?udmgh,udmgk,udmg,nest)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra

      dimension rauc2(idiobs),deuc2(idiobs),erauc2(idiobs),
     ?edeuc2(idiobs),pmra(idiobs),pmde(idiobs),epmra(idiobs),
     ?epmde(idiobs),udmgj(idiobs),udmgh(idiobs),udmgk(idiobs),
     ?udmg(idiobs),idat(23)


      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 uraiz


      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/
 
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      idin50=50
      un=95
      izon=288
c

      dmar=10.d0/3600.d0

      ra1=ramin*15.d0-dmar
      ra2=ramax*15.d0+dmar

      key=1

      if (ra2-ra1.gt.180.d0) key=2

      de1=demin-dmar
      de2=demax+dmar


c
c     Leitura das faixas de declinacao de 0.5 em 0.5 graus
c     do UCAC2
c


      dfaixa=de1-0.5d0
      famax=de2+0.5d0

c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.5d0


      if (dfaixa-famax.gt.0.5d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.5d0+1.0d0

c
c     Limite de declinacao do UCAC2 na parte norte
c
c     So tem zonas ate arquivo z288
c

      if (j.gt.izon) go to 30

c
c     Monta nome do arquivo da faixa, na leitura de 0.5 em 0.5 graus
c


      ifaixa=''
      ifaixa=uraiz


      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC2 ',a54)



c
c     Le a faixa de 0.5 da vez
c

      open (95,file=ifaixa,access="direct",form="unformatted",recl=44)


      n=0
 20   n=n+1

c
c     RA,DEC eh passado para graus
c     erro em RA,DEC passado para segundos de arco
c     movimentos proprios passados para segundo de arco por ano
c    

      call readu2 (un,n,idat,ierra)
      if (ierra.eq.1) go to 25




c
c     Checa se estrela cai dentro do campo
c


      de=idat(2)/3600d3

      if (de.lt.de1) go to 20
      if (de.gt.de2) go to 20


      ra=idat(1)/3600d3


      if (key.eq.1) then

      if (ra.lt.ra1) go to 20
      if (ra.gt.ra2) go to 20
      
      else


      if (ra.gt.ra1 .and. ra.lt.ra2) go to 20

      endif



c
c     Guarda dados das estrelas
c
c     - ra,de em graus na epoca epoj da observacao ccd
c     - epram, epdem epoca Juliana da posicao media RA DEC no UCAC2
c     - era, ede erros de RA DE para epoca media do UCAC2 em segundos de arco
c     - epmRA, epmDe erros mov proprios em segundos de arco por ano
c     - erra, erde erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmra mov poprio (*dcosD) em segundos de arco por ano
c     - pmde mov proprio delta em segundos de arco por ano
c     - epmx, epmy erro dos mov proprios em RA DE em segundos de arco por ano
c     - zmg ... magnitudes 2MASS no UCAC2
c     - dmg  magnitude interna UCAC2 entre V e R
c
c
c

      pmx=idat(12)
      pmy=idat(13)

      ra=idat(1)+pmx*(epoj-2000d0)/10d0
      de=idat(2)+pmy*(epoj-2000d0)/10d0

      ra=ra/3600d3
      de=de/3600d3


      epram=(idat(10)/1000d0)+1975d0
      epdem=(idat(11)/1000d0)+1975d0

      era=idat(4)/1000d0
      ede=idat(5)/1000d0

      epmx =(idat(14)/10d0)/1000d0
      epmy =(idat(15)/10d0)/1000d0

      erra=dsqrt(era**2+(epmx*(epoj-epram))**2)
      erde=dsqrt(ede**2+(epmy*(epoj-epdem))**2)

      zmgj=idat(19)/1000d0
      zmgh=idat(20)/1000d0
      zmgk=idat(21)/1000d0

      dmg=idat(3)/100.d0


      nest=nest+1

      rauc2(nest)=ra
      deuc2(nest)=de
      erauc2(nest)=erra
      edeuc2(nest)=erde

      pmde(nest)=(pmy/10d0)/1000d0
      pmra(nest)=((pmx/10d0)/1000d0)*dcos(dabs(grarad*de))
      epmra(nest)=epmx
      epmde(nest)=epmy


      udmgj(nest)=zmgj
      udmgh(nest)=zmgh
      udmgk(nest)=zmgk
      udmg(nest)=dmg


c     ddja=epram
c     ddjd=epdem




c
c     debug alfa delta
c

c     ra=ra/15.d0
c     dmag=mag/100.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,'  ',erra,erde,
c    ?pmra(nest),pmde(nest),epmx,epmy,udmgj(nest),
c    ?udmgh(nest),udmgk(nest),udmg(nest)




      go to 20

 25   close (95)


 30   continue

 35   continue

 
      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of UCAC2 stars.' 
      write (*,*)
      write (*,*)
      endif



      return
      end



C***********************************************************************
C
C  read a single record of UCAC2 data = 1 star
C  input:
C    un   = unit number of file  (assumed to be open)
C    recn = record number on that file
C    bf   = .TRUE. if byte flip required
C  output:
C    idat = integer*4 vector of 23 items (see readme2.txt)
C    errflg = 0=ok, 1=erro
C
C
C   Modificada: M. Assafin 12/Dez/2005
C
C



      SUBROUTINE readu2 (un,recn,idat,ierra)

      IMPLICIT NONE

      integer*2 un,ierra
      INTEGER recn, idat(23)  

      INTEGER  ra2000, dc2000, pmx,pmy, id2m, u2id,r11
      INTEGER*2 mag, cepx,cepy, j2m,h2m,k2m
      BYTE      sigx,sigy,nobs,epos,ncat,cflg,spmx,spmy, rx,ry, ph,cc          

      ierra = 0                           


      READ (un,REC=recn,ERR=99) ra2000,dc2000,mag,sigx,sigy,nobs,epos,
     ?ncat,cflg,cepx,cepy, pmx,pmy, spmx,spmy, rx,ry,id2m, j2m,h2m,k2m,
     ?ph,cc

c note: first assign I*1 to idat(I*4), 
c       then add 127 to avoid overflow

      idat ( 1) = ra2000
      idat ( 2) = dc2000
      idat ( 3) = mag
      idat ( 4) = sigx
      idat ( 4) = idat ( 4) + 127
      idat ( 5) = sigy
      idat ( 5) = idat ( 5) + 127
      idat ( 6) = nobs
      idat ( 7) = epos 
      idat ( 7) = idat ( 7) + 127
      idat ( 8) = ncat
      idat ( 9) = cflg
      idat (10) = cepx
      idat (11) = cepy
      idat (12) = pmx
      idat (13) = pmy
      idat (14) = spmx
      idat (14) = idat (14) + 127
      idat (15) = spmy
      idat (15) = idat (15) + 127
      idat (16) = rx
      idat (16) = idat (16) + 127
      idat (17) = ry
      idat (17) = idat (17) + 127
      idat (18) = id2m
      idat (19) = j2m
      idat (20) = h2m
      idat (21) = k2m
      idat (22) = ph 
      idat (22) = idat (22) + 127
      idat (23) = cc
      idat (23) = idat (23) + 127

      return

 99   ierra=1

      RETURN
      END






c
c
c     Subrotina sucac2
c
c
c     Pega as estrelas do catalogo UCAC2.
c
c
c     - rauc2,deuc2 em graus (alfas e deltas do UCAC2)
c     - erauc2,edeuc2 em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: ucac2, J, H e K.
c
c     - inuu2,ir1u2,ir2u2: indexes de estrelas do UCAC2
c
c
c     Nesta versao, a leitura eh feita mais rapida, sem projecao gnomonica, 
c     usando a indexacao do catalogo. Usar a rotina ucac2 qdo o polo cai dentro
c     do campo, quando entao eh obrigatorio o uso da projecao gnomonica, estrela
c     a estrela.
c
c
c     Mod. M. Assafin 03/Jul/2013
c
c

      subroutine sucac2 (idiobs,iu2z,iu2i,inuu2,ir1u2,ir2u2,uraiz,epoj,
     ?rac,dec,drac,ddec,ramin,ramax,demin,demax,rauc2,deuc2,erauc2,
     ?edeuc2,pmra,pmde,epmra,epmde,udmgj,udmgh,udmgk,udmg,nest)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra

      dimension rauc2(idiobs),deuc2(idiobs),erauc2(idiobs),
     ?edeuc2(idiobs),pmra(idiobs),pmde(idiobs),epmra(idiobs),
     ?epmde(idiobs),udmgj(idiobs),udmgh(idiobs),udmgk(idiobs),
     ?udmg(idiobs),idat(23)

      dimension inuu2(iu2z,iu2i),ir1u2(iu2z,iu2i),ir2u2(iu2z,iu2i)

      dimension jxmin(2),jxmax(2),cxmin(2),cxmax(2)

      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 uraiz


      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/
 
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))


C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI


      idin50=50
      un=95
      izon=iu2z
c


      dmar=10.d0/3600.d0

      ra1=ramin-dmar/15.d0
      ra2=ramax+dmar/15.d0

      key=1

      if (ra2-ra1.gt.12.d0) key=2

      de1=demin-dmar
      de2=demax+dmar


      bin=0.1d0
      nbin=iu2i

c
c     Zerando vetores
c

      do i=1,2
      jxmin(i)=0
      jxmax(i)=0
      cxmin(i)=0
      cxmax(i)=0
      enddo




c
c     Leitura das faixas de declinacao de 0.5 em 0.5 graus
c     do UCAC2
c

      dfaixa=de1-0.5d0
      famax=de2+0.5d0


c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.5d0


      if (dfaixa-famax.gt.0.5d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.5d0+1.0d0

      jaca=j

c
c     Limite de declinacao do UCAC2 na parte norte
c
c     So tem zonas ate arquivo z288
c

      if (j.gt.izon) go to 30

c
c     Monta nome do arquivo da faixa, na leitura de 0.5 em 0.5 graus
c



      ifaixa=''
      ifaixa=uraiz


      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC2 ',a54)




c
c     Calcula range minimo e maximo de alfa para extracao da faixa da vez
c     (limites alfa indiferentes do hemisferio)
c



      if (key.eq.1) then

      inde=1

      jxmin(1)=ra1/bin+1
      jxmax(1)=ra2/bin+1

      cxmin(1)=ra1*54.d6
      cxmax(1)=ra2*54.d6


      else

      inde=2

      jxmin(1)=ra2/bin+1
      jxmax(1)=nbin

      jxmin(2)=1
      jxmax(2)=ra1/bin+1

      cxmin(1)=ra1*54.d6
      cxmax(1)=24.d0*54.d6
      cxmin(2)=0.d0
      cxmax(2)=ra1*54.d6


      endif



c
c     Le a faixa de 0.5 da vez
c

      open (95,file=ifaixa,access="direct",form="unformatted",recl=44)



      do 20 in=1,inde

      do 19 nn=jxmin(in),jxmax(in)


      if(inuu2(jaca,nn).eq.0) go to 19

      i1=ir1u2(jaca,nn)
      i2=ir2u2(jaca,nn)


      do 18 n=i1,i2


c
c     RA,DEC eh passado para graus
c     erro em RA,DEC passado para segundos de arco
c     movimentos proprios passados para segundo de arco por ano
c    

      call readu2 (un,n,idat,ierra)
c     if (ierra.eq.1) go to 25


      de=idat(2)/3600d3

      if (de.lt.de1) go to 18
      if (de.gt.de2) go to 18


      ra=idat(1)

c

      if (ra.lt.cxmin(in)) go to 18
      if (ra.gt.cxmax(in)) go to 18



c
c     Guarda dados das estrelas
c
c     - ra,de em graus na epoca epoj da observacao ccd
c     - epram, epdem epoca Juliana da posicao media RA DEC no UCAC2
c     - era, ede erros de RA DE para epoca media do UCAC2 em segundos de arco
c     - epmRA, epmDe erros mov proprios em segundos de arco por ano
c     - erra, erde erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmra mov poprio (*dcosD) em segundos de arco por ano
c     - pmde mov proprio delta em segundos de arco por ano
c     - epmx, epmy erro dos mov proprios em RA DE em segundos de arco por ano
c     - zmg ... magnitudes 2MASS no UCAC2
c     - dmg  magnitude interna UCAC2 entre V e R
c
c
c
c



      pmx=idat(12)
      pmy=idat(13)

      ra=idat(1)+pmx*(epoj-2000d0)/10d0
      de=idat(2)+pmy*(epoj-2000d0)/10d0

      ra=ra/1000d0
      de=de/1000d0

      ra=ra/3600d0
      de=de/3600d0

      epram=(idat(10)/1000d0)+1975d0
      epdem=(idat(11)/1000d0)+1975d0

      era=idat(4)/1000d0
      ede=idat(5)/1000d0

      epmx =(idat(14)/10d0)/1000d0
      epmy =(idat(15)/10d0)/1000d0

      erra=dsqrt(era**2+(epmx*(epoj-epram))**2)
      erde=dsqrt(ede**2+(epmy*(epoj-epdem))**2)

      zmgj=idat(19)/1000d0
      zmgh=idat(20)/1000d0
      zmgk=idat(21)/1000d0

      dmg=idat(3)/100.d0


      nest=nest+1

      rauc2(nest)=ra
      deuc2(nest)=de
      erauc2(nest)=erra
      edeuc2(nest)=erde

      pmde(nest)=(pmy/10d0)/1000d0
      pmra(nest)=((pmx/10d0)/1000d0)*dcos(dabs(grarad*de))
      epmra(nest)=epmx
      epmde(nest)=epmy


      udmgj(nest)=zmgj
      udmgh(nest)=zmgh
      udmgk(nest)=zmgk
      udmg(nest)=dmg


c     ddja=epram
c     ddjd=epdem




c
c     debug alfa delta
c

c     ra=ra/15.d0
c     dmag=mag/100.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,'  ',erra,erde,
c    ?pmra(nest),pmde(nest),epmx,epmy,udmgj(nest),
c    ?udmgh(nest),udmgk(nest),udmg(nest)



 18   continue

 19   continue

 20   continue


      close (95)


 30   continue

 35   continue

 
      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of UCAC2 stars.' 
      write (*,*)
      write (*,*)
      endif


      return
      end






c
c
c     Subrotina ucac4
c
c
c     Pega as estrelas do catalogo UCAC4.
c
c
c     - rauc4,deuc4 em graus (alfas e deltas do UCAC4)
c     - erauc4,edeuc4 em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: ucac4, J, H e K.
c
c
c     Mod. M. Assafin 03/Jul/2013
c
c

      subroutine ucac4 (idiobs,u4raiz,epoj,rac,dec,drac,ddec,ramin,
     ?ramax,demin,demax,rauc4,deuc4,erauc4,edeuc4,pmra4,pmde4,epmra4,
     ?epmde4,udmgj4,udmgh4,udmgk4,udmg4,nest)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra,reclen


      LOGICAL  bf, eozf
      INTEGER*4 ran,spdn, id2m,rnm, mcf, rnz
      INTEGER*2 magm,maga, cepra,cepdc, pmra2,pmdc2
     .         ,jmag,hmag,kmag, apasm(5), zn2
      INTEGER*1 sigmag, sigra,sigdc, sigpmr,sigpmd
      INTEGER*1 ojt,dsf, na1,nu1,us1, apase(5), gcflg
      INTEGER*1 icqflg(3), q2mflg(3), leda,x2m


      dimension rauc4(idiobs),deuc4(idiobs),erauc4(idiobs),
     ?edeuc4(idiobs),pmra4(idiobs),pmde4(idiobs),epmra4(idiobs),
     ?epmde4(idiobs),udmgj4(idiobs),udmgh4(idiobs),udmgk4(idiobs),
     ?udmg4(idiobs)


      dimension irnm(25),ipmrc(25),ipmd(25)


      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 u4raiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/
 
 
      data irnm /1,2,200137,200168,200169,200229,200400,200503,200530,
     ?200895,201050,201349,201526,201550,201567,201633,201803,249921,
     ?249984,268357,80118783,93157181,106363470,110589580,113038183/ 

      data ipmrc /41087,41558,-37758,-36004,-36782,39624,65051,56347,
     ?67682,-22401,-7986,22819,-5803,40033,41683,-44099,34222,-10015,
     ?-9994,5713,32962,-22393,-37060,-38420,10990/

      data ipmd /31413,32586,7655,9521,4818,-25374,-57308,-23377,13275,
     ?-34203,103281,53694,-47659,-58151,32691,9416,-15989,-35427,-35419,
     ?-36943,5639,-34199,-11490,-27250,-51230/


      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI


      idin50=50
      un=95
      izon=900
      reclen=78

c


      dmar=10.d0/3600.d0

      ra1=ramin*15.d0-dmar
      ra2=ramax*15.d0+dmar

      key=1

      if (ra2-ra1.gt.180.d0) key=2

      de1=demin-dmar
      de2=demax+dmar


c
c     Leitura das faixas de declinacao de 0.2 em 0.2 graus
c     do UCAC4
c


      dfaixa=de1-0.2d0
      famax=de2+0.2d0

c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.2d0


      if (dfaixa-famax.gt.0.2d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.2d0+1.0d0

c
c     Limite de declinacao do UCAC4
c
c     Zonas ate arquivo z900 
c

      if (j.gt.izon) go to 30


c
c     Monta nome do arquivo da faixa, na leitura de 0.2 em 0.2 graus
c


      ifaixa=''
      ifaixa=u4raiz


      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC4 ',a54)



c
c     Le a faixa de 0.2 da vez
c


      open (95,file=ifaixa,access="direct",recl=reclen)


      n=0
 20   n=n+1

c
c     RA,DEC eh passado para graus
c     erro em RA,DEC passado para segundos de arco
c     movimentos proprios passados para segundo de arco por ano
c    


c

      READ (un,REC=n,ERR=25)                      !          sum = 78
     .     ran,spdn, magm,maga, sigmag,ojt,dsf    !  8 + 4 +  3  = 15
     .    ,sigra,sigdc, na1,nu1,us1               !  2 + 3       =  5
     .    ,cepra,cepdc, pmra2,pmdc2,sigpmr,sigpmd !  4 + 4 +  2  = 10
     .    ,id2m, jmag,hmag,kmag, icqflg, q2mflg   !  4 + 6 +  6  = 16
     .    ,apasm, apase, gcflg                    ! 10 + 5 +  1  = 16
     .    ,mcf, leda,x2m, rnm                     !  4 + 2 +  4  = 10
     .    ,zn2, rnz                               !  2 + 4       =  6




c
c     Checa se estrela cai dentro do campo
c

      de=spdn/3.6d6-90.d0

      if (de.lt.de1) go to 20
      if (de.gt.de2) go to 20


      ra=ran/3.6d6


      if (key.eq.1) then

      if (ra.lt.ra1) go to 20
      if (ra.gt.ra2) go to 20
      
      else


      if (ra.gt.ra1 .and. ra.lt.ra2) go to 20

      endif




c
c     Guarda dados das estrelas
c
c     - ra,de em graus na epoca epoj da observacao ccd
c     - epram3, epdem3 epoca Juliana da posicao media RA DEC no UCAC3
c     - era4, ede4 erros de RA DE para epoca media do UCAC3 em segundos de arco
c     - epmra4, epmde4 erros mov proprios em segundos de arco por ano
c     - erra3, erde3 erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmra4 mov poprio (*dcosD) em segundos de arco por ano
c     - pmde4 mov proprio delta em segundos de arco por ano
c     - epmx3, epmy3 erro dos mov proprios em RA DE em segundos de arco por ano
c     - udmg_3 ... magnitudes 2MASS no UCAC3
c     - udmg4  magnitude interna (psf) UCAC3 entre V e R
c
c
c




      if (pmra2.eq.32767 .and. pmdc2.eq.32767) then

      do m=1,25 
      if (rnm.eq.irnm(m)) go to 23
      enddo

 23   continue

      pmx=ipmrc(m)/dcos((de)*grarad)
      pmy=ipmd(m)


      else

      pmx=pmra2/dcos((de)*grarad)
      pmy=pmdc2

      endif


      ra=ran+pmx*(epoj-2000d0)/10d0
      de=spdn+pmy*(epoj-2000d0)/10d0

      ra=ra/3.6d6
      de=de/3.6d6
      de=de-90.d0



      epram4=cepra/1.d2+1990d0
      epdem4=cepdc/1.d2+1990d0


      esigra=sigra+128
      esigdc=sigdc+128

      era4=esigra/1.d3
      ede4=esigdc/1.d3


      sipmr=sigpmr+128.d0
      sipmd=sigpmd+128.d0


      if (sipmr.eq.251) sipmr=275
      if (sipmr.eq.252) sipmr=325 
      if (sipmr.eq.253) sipmr=375 
      if (sipmr.eq.254) sipmr=450 
      if (sipmr.eq.255) sipmr=500 


      if (sipmd.eq.251) sipmr=275
      if (sipmd.eq.252) sipmr=325 
      if (sipmd.eq.253) sipmr=375 
      if (sipmd.eq.254) sipmr=450 
      if (sipmd.eq.255) sipmr=500 



      epmx4=sipmr/10.d3
      epmy4=sipmd/10.d3

      erra4=dsqrt(era4**2+(epmx4*(epoj-epram4))**2)
      erde4=dsqrt(ede4**2+(epmy4*(epoj-epdem4))**2)

      zmgj4=jmag/1.d3
      zmgh4=hmag/1.d3
      zmgk4=kmag/1.d3

      dmg4=maga/1.d3


      nest=nest+1

      rauc4(nest)=ra
      deuc4(nest)=de
      erauc4(nest)=erra4
      edeuc4(nest)=erde4

      pmde4(nest)=(pmy/10d0)/1000d0
      pmra4(nest)=(pmx/10d0)/1000d0
      epmra4(nest)=epmx4
      epmde4(nest)=epmy4


      udmgj4(nest)=zmgj4
      udmgh4(nest)=zmgh4
      udmgk4(nest)=zmgk4
      udmg4(nest)=dmg4



c
c     debug alfa delta
c
c
c     ra=ra/15.d0
c     dmag=dmg4
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF 
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0
c
c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,'  ',erra4,erde4,
c    ?pmra4(nest),pmde4(nest),epmx4,epmy4,udmgj4(nest),
c    ?udmgh4(nest),udmgk4(nest),udmg4(nest)
c



      go to 20

 25   close (95)


 30   continue

 35   continue

 
      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of UCAC4 stars.' 
      write (*,*)
      write (*,*)
      endif

      return
      end







c
c
c     Subrotina sucac4
c
c
c     Pega as estrelas do catalogo UCAC4.
c
c
c     - rauc4,deuc4 em graus (alfas e deltas do UCAC4)
c     - erauc4,edeuc4 em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: ucac4, J, H e K.
c
c     - inuu4,ir1u4,ir2u4: indexes de estrelas do UCAC4
c
c
c     Nesta versao, a leitura eh feita mais rapida, sem projecao gnomonica, 
c     usando a indexacao do catalogo. Usar a rotina ucac4 qdo o polo cai dentro
c     do campo, quando entao eh obrigatorio o uso da projecao gnomonica, estrela
c     a estrela.
c
c
c     Mod. M. Assafin 03/Jul/2013
c
c

      subroutine sucac4 (idiobs,iu4z,iu4i,inuu4,ir1u4,ir2u4,u4raiz,epoj,
     ?rac,dec,drac,ddec,ramin,ramax,demin,demax,rauc4,deuc4,erauc4,
     ?edeuc4,pmra4,pmde4,epmra4,epmde4,udmgj4,udmgh4,udmgk4,udmg4,nest)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra,reclen


      LOGICAL  bf, eozf
      INTEGER*4 ran,spdn, id2m,rnm, mcf, rnz
      INTEGER*2 magm,maga, cepra,cepdc, pmra2,pmdc2
     .         ,jmag,hmag,kmag, apasm(5), zn2
      INTEGER*1 sigmag, sigra,sigdc, sigpmr,sigpmd
      INTEGER*1 ojt,dsf, na1,nu1,us1, apase(5), gcflg
      INTEGER*1 icqflg(3), q2mflg(3), leda,x2m


      dimension rauc4(idiobs),deuc4(idiobs),erauc4(idiobs),
     ?edeuc4(idiobs),pmra4(idiobs),pmde4(idiobs),epmra4(idiobs),
     ?epmde4(idiobs),udmgj4(idiobs),udmgh4(idiobs),udmgk4(idiobs),
     ?udmg4(idiobs)

      dimension irnm(25),ipmrc(25),ipmd(25)

      dimension inuu4(iu4z,iu4i),ir1u4(iu4z,iu4i),ir2u4(iu4z,iu4i)

      dimension jxmin(2),jxmax(2),cxmin(2),cxmax(2)


      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 u4raiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/
 
      data irnm /1,2,200137,200168,200169,200229,200400,200503,200530,
     ?200895,201050,201349,201526,201550,201567,201633,201803,249921,
     ?249984,268357,80118783,93157181,106363470,110589580,113038183/ 

      data ipmrc /41087,41558,-37758,-36004,-36782,39624,65051,56347,
     ?67682,-22401,-7986,22819,-5803,40033,41683,-44099,34222,-10015,
     ?-9994,5713,32962,-22393,-37060,-38420,10990/

      data ipmd /31413,32586,7655,9521,4818,-25374,-57308,-23377,13275,
     ?-34203,103281,53694,-47659,-58151,32691,9416,-15989,-35427,-35419,
     ?-36943,5639,-34199,-11490,-27250,-51230/


      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      idin50=50

      un=95
      izon=iu4z
      reclen=78

      bin=0.25d0
      nbin=iu4i


c

      dmar=10.d0/3600.d0

      ra1=ramin*15.d0-dmar
      ra2=ramax*15.d0+dmar

      key=1

      if (ra2-ra1.gt.180.d0) key=2

      de1=demin-dmar
      de2=demax+dmar



c
c     Zerando vetores
c

      do i=1,2
      jxmin(i)=0
      jxmax(i)=0
      cxmin(i)=0
      cxmax(i)=0
      enddo


c
c     Leitura das faixas de declinacao de 0.2 em 0.2 graus
c     do UCAC4
c


      dfaixa=de1-0.2d0
      famax=de2+0.2d0



c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.2d0


      if (dfaixa-famax.gt.0.2d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.2d0+1.0d0

      jaca=j

c
c     Limite de declinacao do UCAC4
c
c     Zonas ate arquivo z900 
c

      if (j.gt.izon) go to 30

c
c     Monta nome do arquivo da faixa, na leitura de 0.2 em 0.2 graus
c


      ifaixa=''
      ifaixa=u4raiz


      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC4 ',a54)



c
c     Calcula range minimo e maximo de alfa para extracao da faixa da vez
c     (limites alfa indiferentes do hemisferio)
c



      if (key.eq.1) then

      inde=1

      jxmin(1)=ra1/bin+1
      jxmax(1)=ra2/bin+1

      cxmin(1)=ra1*3600.d3
      cxmax(1)=ra2*3600.d3


      else


      inde=2

      jxmin(1)=ra2/bin+1
      jxmax(1)=nbin

      jxmin(2)=1
      jxmax(2)=ra1/bin+1

      cxmin(1)=ra2*3600.d3
      cxmax(1)=360.d0*3600.d3
      cxmin(2)=0.d0
      cxmax(2)=ra1*3600.d3


      endif



c
c     Le a faixa de 0.2 graus da vez
c


      open (95,file=ifaixa,access="direct",recl=reclen)



      do 20 in=1,inde

      do 19 nn=jxmin(in),jxmax(in)


      if(inuu4(jaca,nn).eq.0) go to 19

      i1=ir1u4(jaca,nn)
      i2=ir2u4(jaca,nn)



      do 18 n=i1,i2




c
c     RA,DEC eh passado para graus
c     erro em RA,DEC passado para segundos de arco
c     movimentos proprios passados para segundo de arco por ano
c    


c

      READ (un,REC=n)                             !          sum = 78
     .     ran,spdn, magm,maga, sigmag,ojt,dsf    !  8 + 4 +  3  = 15
     .    ,sigra,sigdc, na1,nu1,us1               !  2 + 3       =  5
     .    ,cepra,cepdc, pmra2,pmdc2,sigpmr,sigpmd !  4 + 4 +  2  = 10
     .    ,id2m, jmag,hmag,kmag, icqflg, q2mflg   !  4 + 6 +  6  = 16
     .    ,apasm, apase, gcflg                    ! 10 + 5 +  1  = 16
     .    ,mcf, leda,x2m, rnm                     !  4 + 2 +  4  = 10
     .    ,zn2, rnz                               !  2 + 4       =  6


c
c    Checa se estrela cai dentro do campo
c


      de=spdn/3.6d6-90.d0

      if (de.lt.de1) go to 18
      if (de.gt.de2) go to 18



      ra=ran


c

      if (ra.lt.cxmin(in)) go to 18
      if (ra.gt.cxmax(in)) go to 18



c
c     Guarda dados das estrelas
c
c     - ra,de em graus na epoca epoj da observacao ccd
c     - epram3, epdem3 epoca Juliana da posicao media RA DEC no UCAC3
c     - era4, ede4 erros de RA DE para epoca media do UCAC3 em segundos de arco
c     - epmra4, epmde4 erros mov proprios em segundos de arco por ano
c     - erra3, erde3 erro da posicao RA DE para a epoca Juliana do
c       CCD em segundos de arco
c     - pmra4 mov poprio (*dcosD) em segundos de arco por ano
c     - pmde4 mov proprio delta em segundos de arco por ano
c     - epmx3, epmy3 erro dos mov proprios em RA DE em segundos de arco por ano
c     - udmg_3 ... magnitudes 2MASS no UCAC3
c     - udmg4  magnitude interna (psf) UCAC3 entre V e R
c
c
c
c


      if (pmra2.eq.32767 .and. pmdc2.eq.32767) then

      do m=1,25 
      if (rnm.eq.irnm(m)) go to 16
      enddo

 16   continue

      pmx=ipmrc(m)/dcos((de)*grarad)
      pmy=ipmd(m)


      else

      pmx=pmra2/dcos((de)*grarad)
      pmy=pmdc2

      endif


      ra=ran+pmx*(epoj-2000d0)/10d0
      de=spdn+pmy*(epoj-2000d0)/10d0

      ra=ra/3.6d6
      de=de/3.6d6
      de=de-90.d0



      epram4=cepra/1.d2+1990d0
      epdem4=cepdc/1.d2+1990d0


      esigra=sigra+128
      esigdc=sigdc+128

      era4=esigra/1.d3
      ede4=esigdc/1.d3


      sipmr=sigpmr+128.d0
      sipmd=sigpmd+128.d0


      if (sipmr.eq.251) sipmr=275
      if (sipmr.eq.252) sipmr=325 
      if (sipmr.eq.253) sipmr=375 
      if (sipmr.eq.254) sipmr=450 
      if (sipmr.eq.255) sipmr=500 


      if (sipmd.eq.251) sipmr=275
      if (sipmd.eq.252) sipmr=325 
      if (sipmd.eq.253) sipmr=375 
      if (sipmd.eq.254) sipmr=450 
      if (sipmd.eq.255) sipmr=500 



      epmx4=sipmr/10.d3
      epmy4=sipmd/10.d3

      erra4=dsqrt(era4**2+(epmx4*(epoj-epram4))**2)
      erde4=dsqrt(ede4**2+(epmy4*(epoj-epdem4))**2)

      zmgj4=jmag/1.d3
      zmgh4=hmag/1.d3
      zmgk4=kmag/1.d3

      dmg4=maga/1.d3


      nest=nest+1

      rauc4(nest)=ra
      deuc4(nest)=de
      erauc4(nest)=erra4
      edeuc4(nest)=erde4

      pmde4(nest)=(pmy/10d0)/1000d0
      pmra4(nest)=(pmx/10d0)/1000d0
      epmra4(nest)=epmx4
      epmde4(nest)=epmy4


      udmgj4(nest)=zmgj4
      udmgh4(nest)=zmgh4
      udmgk4(nest)=zmgk4
      udmg4(nest)=dmg4




c
c     debug alfa delta
c
c
c     ra=ra/15.d0
c     dmag=dmg4
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF 
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0
c
c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,'  ',erra4,erde4,
c    ?pmra4(nest),pmde4(nest),epmx4,epmy4,udmgj4(nest),
c    ?udmgh4(nest),udmgk4(nest),udmg4(nest)
c



 18   continue

 19   continue

 20   continue


      close (95)


 30   continue

 35   continue

      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of UCAC4 stars.' 
      write (*,*)
      write (*,*)
      endif


      return
      end




c
c
c     Subrotina cuser
c
c
c     Pega as estrelas do catalogo de referencia do usuario (fomrato PRAIA)
c
c
c     - raucs,deucs em graus (alfas e deltas do catalogo, na epoca epoj do CCD
c     - eraucs,edeucs em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: mas do catalogo, e mags J, H e K.
c
c
c     Mod. M. Assafin 14/Nov/2012
c
c



      subroutine cuser (idiobs,ifaixa,epoj,rac,dec,drac,ddec,ramin,
     ?ramax,demin,demax,raucs,deucs,eraucs,edeucs,pmras,pmdes,epmras,
     ?epmdes,udmgjs,udmghs,udmgks,udmgs,nest)



      implicit real*8 (a-h,o-z)


      dimension raucs(idiobs),deucs(idiobs),eraucs(idiobs),
     ?edeucs(idiobs),pmras(idiobs),pmdes(idiobs),epmras(idiobs),
     ?epmdes(idiobs),udmgjs(idiobs),udmghs(idiobs),udmgks(idiobs),
     ?udmgs(idiobs)


      CHARACTER *1 isig
      character *50 ifaixa

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0
  
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)

c
c     Initializing data
c

      pi=3.141592653589793d0
      grarad=pi/180.d0
      radgra=180.d0/pi


      idin50=50
      iun=95

c

      xbox=grarad*drac
      ybox=grarad*ddec
      grac=grarad*rac
      gdec=grarad*dec

c


      nest=0

c

      write (*,14) ifaixa
 14   format(1x,'User  ',a50)



c
c     Le o catalogo
c


      open (iun,file=ifaixa)


      nest=0



 20   continue


      read (iun,21,end=25) iah,iam,sa,isig,idg,idm,sd,ex,ey,codj,copma,
     ?copmd,erpma,erpmd,codmg,res2mg

 21   format(i2,1x,i2,1x,f7.4,2x,a1,i2,1x,i2,1x,f6.3,2x,2f7.3,2x,f16.8,
     ?1x,4(1x,f7.3),2(1x,f6.3))

c

      ra=15.d0*hmsgms(iah,iam,sa)

      de=hmsgms(idg,idm,sd)
      if (isig.eq.'-') de=-de



c
c     Checa se estrela cai dentro do campo
c

      bra=grarad*ra
      bde=grarad*de

      d=DEXY(bra,bde,grac,gdec)
      xx=XPAD(bra,bde,grac)/d
      yy=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xx).gt.xbox) go to 20
      if (dabs(yy).gt.ybox) go to 20

c
c     Guarda dados das estrelas
c



      if (copma.lt.90.d0 .and. copmd.lt.90.d0) then

      cop=2000D0+(codj-2451545D0)/365.25D0
      dt=epoj-cop

      de=de+(copmd*dt)/3600.d0
      ra=ra+(copma*dt/dcos(grarad*dabs(de)))/3600.d0

      endif


      nest=nest+1

      raucs(nest)=ra
      deucs(nest)=de
      eraucs(nest)=ex
      edeucs(nest)=ey

      pmras(nest)=copma
      pmdes(nest)=copmd
      epmras(nest)=erpma
      epmdes(nest)=erpmd


      udmgjs(nest)=99.9d0
      udmghs(nest)=99.9d0
      udmgks(nest)=99.9d0
      udmgs(nest)=codmg



c
c     debug alfa delta
c
c
c     ra=ra/15.d0
c     dmag=dmg4
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF 
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0
c
c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,'  ',ex,ey,
c    ?pmras(nest),pmdes(nest),erpma,erpmd,udmgjs(nest),
c    ?udmghs(nest),udmgks(nest),udmgs(nest)
c



      go to 20

 25   close (95)

 
      if (nest.gt.idiobs) then
      write (*,*)
      write (*,*)
      write (*,*)'Attention: overflow in the number of User-catalog star
     ?s.' 
      write (*,*)
      write (*,*)
      endif



      return
      end






C
C
C
C
C     SUBROTINA ISOLO
C
C     PROPOSITO
C
C     Ajustar polinomio bivariado P=P(x,y) de grau n a fundo de ceu
C		    
C       onde  x  -  coordenada x do pixel 
C             y  -  coordenada y do pixel
C							     
C
C     USO
C
C     SUBROUTINE ISOLO(IENTRA,NUMEST,IVER,XPIXEL,YPIXEL,NGRAU,
C       NGRAU3,NGRAU5,NCOMUM,XSAO,YSAO,XEST1,YEST1,XP,YPNTIRA,ITIRA,
C       COEFX,COEFY,XSIG,YSIG,XRRAY,YRRAY)
C
C     DESCRICAO DOS PARAMETROS
C
C       NGRAU  - GRAU DO POLINOMIO DE AJUSTE
C                0 = 4ctes
C                1 = 1o grau completo
C                2 = 2o grau completo
C                3 = 3o grau completo
c                    etc (maximo N=15)
C
C       NCOMUM - NUMERO DE ESTRELAS COMUNS PARA AJUSTE
C       NUMEST - NUMERO TOTAL DE ESTRELAS DO L-ESIMO CCD
C       NTIRA  - NUMERO DE ESTRELAS RETIRADAS DO AJUSTE
C       ITIRA  - ESTRELAS RETIRADAS DO AJUSTE
C       XSAO   - MEDIDA X DO L-ESIMO CCD NA TELA SAOIMAGE
C       YSAO   - MEDIDA Y DO L-ESIMO CCD NA TELA SAOIMAGE
C       XEST1  - MEDIDA X NO L-ESIMO CCD
C       YEST1  - MEDIDA Y NO L-ESIMO CCD
C       XP     - X DE REFERENCIA DO MOSAICO CENTRAL
C       YP     - Y DE REFERENCIA DO MOSAICO CENTRAL
C       COEFX
C       COEFY  - COEFICIENTES DO POLINOMIO
C       ITIRA  - ESTRELAS RETIRADAS COM MAIOR (O-C)
C       VAX    - ERRO PADRAO X DA RADIOESTRELA
C       VAY    - ERRO PADRAO Y DA RADIOESTRELA
C       XSIG   - DESVIO-PADRAO X DA RADIOESTRELA
C       YSIG   - DESVIO-PADRAO Y DA RADIOESTRELA
C       SIGMA  - TRUNCAMENTO DE VARIACAO DE COORDENADA (X,Y) EM ("),
C                NO CASO DE ELIMINACAO AUTOMATICA
C
C     SUBROTINAS E SUBPROGRAMAS REQUERIDOS
C
C     MATINV (NTERMS,DET)
C
C     COMENTARIOS
C
C     A ROTINA NAO SUPORTA TERMOS RADIAS
C
c     Ela ajusta polinomios completos de graus 1 a 15
c
c
c
      SUBROUTINE ISOLO (ipmax,NGRAU,NCOMUM,XESTO,YESTO,PO,COFO)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COFO(136),ALPHAO(136,136),ARRAYO(136,136),
     ?BETAO(136),TERMO(136)
c     DIMENSION XESTO(25010001),YESTO(25010001),PO(25010001)
      DIMENSION XESTO(ipmax*ipmax),YESTO(ipmax*ipmax),PO(ipmax*ipmax)

      COMMON/A77/ARRAYO
      COMMON/A14/IERRO
C
C     INICIALIZACAO DE DADOS
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      IERRO=0
      DET=1.d0

      idim=ipmax*ipmax
      idimg=15


C
C     Zera vetores
C

      NTERMS=1
      DO I=1,idimg
      NTERMS=NTERMS+I+1
      ENDDO

      DO 8 I=1,NTERMS
      BETAO(I) =0.D0
      COFO(I) =0.D0
      TERMO(I) =0.D0
      DO 8 J=I,NTERMS
      arrayo(i,j)=0.d0
 8    ALPHAO(I,J)=0.D0


C
C     CALCULA NUMERO DE TERMOS DO POL. DE GRAU=NGRAU
C

      NTERMS=1
      DO I=1,NGRAU
      NTERMS=NTERMS+I+1
      ENDDO

c
c     Ajuste de 1o, 2o, ..., 15o grau.
c

      if (ngrau.lt.0 .or. ngrau.gt.idimg) then
      ierro=1
      return
      endif
C
      IGRAU=NGRAU+1

      ITERMS=NTERMS
C

      DO 265 I=1,NCOMUM
      X=XESTO(I)
      Y=YESTO(I)
      XG=PO(I)

C
C     Computando os termos para o polinomio
C
      ICONT=0
      DO 240 N=1,IGRAU
      DO 240 L=1,N
      K=N-L
      ICONT=ICONT+1
 240  TERMO(ICONT)=(X**K)*(Y**(L-1))

C
C     Computando AtA e AtB
C
      DO 260 L=1,ITERMS
      BETAO(L)=BETAO(L)+XG*TERMO(L)
      DO 260 K=L,ITERMS
  260 ALPHAO(L,K)=ALPHAO(L,K)+TERMO(K)*TERMO(L)
C
 265  CONTINUE

C
C     Preenchendo a parte triangular inferior da matriz simetrica AtA
C

      DO 270 L=1,ITERMS
      DO 270 K=L,ITERMS
  270 ALPHAO(K,L)=ALPHAO(L,K)
c
c
c     Preenchendo ARRAY=AtA para inversao (elementos normalizados
c     pela diagonal) para X ou Y
c
c

      DO 280 L=1,ITERMS
      DO 280 K=1,ITERMS
  280 ARRAYO(L,K)=ALPHAO(L,K)/DSQRT(ALPHAO(L,L)*ALPHAO(K,K))
C
C     Invertendo AtA
C
      CALL MATVO (ITERMS,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando os coeficientes do polinomio
C
      DO 290 L=1,ITERMS
      DO 290 K=1,ITERMS
  290 ARRAYO(L,K)=ARRAYO(L,K)/DSQRT(ALPHAO(K,K)*ALPHAO(L,L))

C
C     Guarda coeficiente
C

      DO 300 L=1,ITERMS
      DO 300 K=1,ITERMS
  300 COFO(L)=COFO(L)+ARRAYO(L,K)*BETAO(K)


      RETURN

      END




C
C     SUBROUTINE MATVO
C
C     PURPOSE
C       INVERT A SYMMETRIC MATRIX AND CALCULATE ITS DETERMINANT
C
C     USAGE
C       CALL MATINV ( NORDER, DET)
C
C     DESCRIPTION OF PARAMETERS
C       ARRAY  - INPUT MATRIX WICH IS REPLACED BY ITS INVERSE
C       NORDER - DEGREE OF MATRIX (ORDER OF DETERMINANT)
C       DET    - DETERMINANT OF INPUT MATRIX
C
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C       NONE
C
C     MODIFICATIONS FOR FORTRAN II
C       OMIT DOUBLE PRECISION SPECIFICATIONS
C       CHANGE DABS TO ABSF IN STATEMENT 23
C
C     COMMENTS
C       DIMENSION STATEMANTS VALID FOR NORDER UP TO 10
C
      SUBROUTINE MATVO (NORDER, DET)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION ARRAYO (136,136), IK(136), JK(136)
      COMMON /A77/ARRAYO
      COMMON /A14/IERRO
C
   10 DET = 1.D0
   11 DO 100 K=1, NORDER
C
C        FIND LARGEST ELEMENT ARRAYO(I,J) IN REST OF MATRIX
C
      AMAX= 0.D0
   21 DO 30 I=K, NORDER
      DO 30 J=K, NORDER
   23 IF (DABS(AMAX) - DABS(ARRAYO(I,J))) 24, 24, 30
   24 AMAX = ARRAYO(I,J)
      IK(K) = I
      JK(K) = J
   30 CONTINUE
C
C        INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAYO(K,K)
C
   31 IF (AMAX) 41, 32, 41
   32 DET = 0.D0
      GO TO 140
   41 I = IK(K)
      IF (I-K) 21, 51, 43
   43 DO 50 J=1, NORDER
      SAVE = ARRAYO(K,J)
      ARRAYO(K,J) = ARRAYO(I,J)
   50 ARRAYO(I,J) = -SAVE
   51 J = JK(K)
      IF (J-K) 21, 61, 53
   53 DO 60 I=1, NORDER
      SAVE = ARRAYO(I,K)
      ARRAYO (I,K) = ARRAYO(I,J)
   60 ARRAYO (I,J) = -SAVE
C
C        ACCUMULATE ELEMENTS OF INVERSE MATRIX
C
   61 DO 70 I=1, NORDER
      IF (I-K) 63, 70, 63
   63 IF (AMAX.EQ.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
      ARRAYO(I,K) = -ARRAYO(I,K) / AMAX
   70 CONTINUE
   71 DO 80 I=1, NORDER
      DO 80 J=1, NORDER
      IF (I-K) 74, 80, 74
   74 IF (J-K) 75, 80, 75
   75 ARRAYO(I,J) = ARRAYO(I,J) + ARRAYO(I,K)*ARRAYO(K,J)
   80 CONTINUE
   81 DO 90 J=1, NORDER
      IF (J-K) 83, 90, 83
   83 IF (AMAX.EQ.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
      ARRAYO(K,J) = ARRAYO(K,J) / AMAX
   90 CONTINUE
      IF (AMAX.EQ.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
      ARRAYO(K,K) = 1.D0 / AMAX
  100 DET = DET * AMAX
C
C        RESTORE ORDERING OF MATRIX
C
  101 DO 130 L=1, NORDER
      K = NORDER - L + 1
      J = IK(K)
      IF (J-K) 111, 111, 105
  105 DO 110 I=1, NORDER
      SAVE = ARRAYO(I,K)
      ARRAYO(I,K) = -ARRAYO(I,J)
  110 ARRAYO(I,J) = SAVE
  111 I = JK(K)
      IF (I-K) 130, 130, 113
  113 DO 120 J=1, NORDER
      SAVE = ARRAYO(K,J)
      ARRAYO(K,J) = -ARRAYO(I,J)
  120 ARRAYO(I,J) = SAVE
  130 CONTINUE
  140 CONTINUE
      RETURN
      END



c
c     Function polo
c
c     Calcula o valor de um polinomio bivariado em (x,y). Grau N
c     1 < N < 15
c
c

      double precision function polo (x, y, cofo, ngrau)

      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION cofo(136)

c
      mgrau=ngrau+1
      icont=0
      polo=0.d0
      do 1 n=1,mgrau
      do 1 l=1,n
      icont=icont+1
      k=n-l
    1 polo=polo+cofo(icont)*(x**k)*(y**(l-1))


      return
      end






C
C
C
C
C     SUBROTINA ISOL
C
C     PROPOSITO
C
C     Ajustar polinomio bivariado P=P(x,y) de grau n
C		    
C       onde  x  -  medida x da estrela
C             y  -  medida y da estrela
C							     
C
C     USO
C
C     SUBROUTINE ISOL (IENTRA,NUMEST,IVER,XPIXEL,YPIXEL,NGRAU,
C       NGRAU3,NGRAU5,NCOMUM,XSAO,YSAO,XEST1,YEST1,XP,YPNTIRA,ITIRA,
C       COEFX,COEFY,XSIG,YSIG,XRRAY,YRRAY)
C
C     DESCRICAO DOS PARAMETROS
C
C       NGRAU  - GRAU DO POLINOMIO DE AJUSTE
C                0 = 4ctes
C                1 = 1o grau completo
C                2 = 2o grau completo
C                3 = 3o grau completo
C
C       NCOMUM - NUMERO DE ESTRELAS COMUNS PARA AJUSTE
C       NUMEST - NUMERO TOTAL DE ESTRELAS DO L-ESIMO CCD
C       NTIRA  - NUMERO DE ESTRELAS RETIRADAS DO AJUSTE
C       ITIRA  - ESTRELAS RETIRADAS DO AJUSTE
C       XSAO   - MEDIDA X DO L-ESIMO CCD NA TELA SAOIMAGE
C       YSAO   - MEDIDA Y DO L-ESIMO CCD NA TELA SAOIMAGE
C       XEST1  - MEDIDA X NO L-ESIMO CCD
C       YEST1  - MEDIDA Y NO L-ESIMO CCD
C       XP     - X DE REFERENCIA DO MOSAICO CENTRAL
C       YP     - Y DE REFERENCIA DO MOSAICO CENTRAL
C       COEFX
C       COEFY  - COEFICIENTES DO POLINOMIO
C       ITIRA  - ESTRELAS RETIRADAS COM MAIOR (O-C)
C       VAX    - ERRO PADRAO X DA RADIOESTRELA
C       VAY    - ERRO PADRAO Y DA RADIOESTRELA
C       XSIG   - DESVIO-PADRAO X DA RADIOESTRELA
C       YSIG   - DESVIO-PADRAO Y DA RADIOESTRELA
C       SIGMA  - TRUNCAMENTO DE VARIACAO DE COORDENADA (X,Y) EM ("),
C                NO CASO DE ELIMINACAO AUTOMATICA
C
C     SUBROTINAS E SUBPROGRAMAS REQUERIDOS
C
C     MATINV (NTERMS,DET)
C
C     COMENTARIOS
C
C     A ROTINA NAO SUPORTA TERMOS RADIAS
C
c     Ela ajusta 4ctes e polinomios completos de graus 1 a 3
c
c
c
      SUBROUTINE ISOL (idiobs,icofsp,NGRAU,NCOMUM,XEST,YEST,XP,YP,COEFX,
     ?COEFY)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(icofsp),COEFY(icofsp),ALPHA(icofsp,icofsp),
     ?ARRAY(icofsp,icofsp),
     ?BETA(icofsp),BETAX(icofsp),TERMX(icofsp)
      DIMENSION XEST(idiobs),YEST(idiobs),XP(idiobs),YP(idiobs)

      COMMON/A14/IERRO
C
C     INICIALIZACAO DE DADOS
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      IERRO=0
      DET=1.d0
C
C     Zera vetores
C

      DO 8 I=1,icofsp
      BETA(I) =0.D0
      BETAX(I) =0.D0
      COEFX(I) =0.D0
      COEFY(I) =0.D0
      TERMX(I) =0.D0
      DO 8 J=I,icofsp
      array(i,j)=0.d0
 8    ALPHA(I,J)=0.D0

c
c     4 ctes ou 1o grau?
c

      if (ngrau.ge.1) go to 200

C
C     AJUSTE POR 4 CONSTANTES
C
      NTERMS=4
      IGRAU=2
C

      ITERMS=NTERMS

C
C     MONTANDO AS EQUACOES DE CONDICAO PARA O AJUSTE POR M.Q. DE Ax=B,
C     NO CASO DE "AJUSTE DE 4 CTES"
C
      DO 10 I=1,NCOMUM
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     CALCULA OS TERMOS DOS COEFICIENTES E AtB
C
      TERMX(1)=TERMX(1)+X**2+Y**2
      TERMX(2)=TERMX(2)+X
      TERMX(3)=TERMX(3)+Y
      TERMX(4)=TERMX(4)+1.D0
      BETAX(1)=BETAX(1)+X*XG+Y*YG
      BETAX(2)=BETAX(2)+Y*XG-X*YG
      BETAX(3)=BETAX(3)+XG
      BETAX(4)=BETAX(4)+YG
   10 CONTINUE
C
C     PREENCHENDO AtA
C
      ALPHA(1,1)=TERMX(1)
      ALPHA(1,2)=0.D0
      ALPHA(1,3)=TERMX(2)
      ALPHA(1,4)=TERMX(3)
      ALPHA(2,2)=TERMX(1)
      ALPHA(2,3)=TERMX(3)
      ALPHA(2,4)=-TERMX(2)
      ALPHA(3,3)=TERMX(4)
      ALPHA(3,4)=0.D0
      ALPHA(4,4)=TERMX(4)
      DO 15 L=1,ITERMS
      DO 15 K=L,ITERMS
   15 ALPHA(K,L) =ALPHA(L,K)
C
C     PREENCHENDO MATRIZ ARRAY=AtA PARA INVERSAO (elementos normalizados
C     pela diagonal) PARA X
C
      DO 80 L=1,ITERMS
      DO 80 K=1,ITERMS
   80 ARRAY(L,K)=ALPHA(L,K)/DSQRT(ALPHA(L,L)*ALPHA(K,K))
C
C     INVERTENDO AtA PARA "X"
C
      CALL MATINV (ITERMS,icofsp,array,DET)
      IF (IERRO.EQ.1) RETURN
C
C     CALCULANDO OS COEFICIENTES DO POLINOMIO PARA "X"
C
      DO 90 L=1,ITERMS
      DO 90 K=1,ITERMS
   90 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHA(K,K)*ALPHA(L,L))
C
C      OBTEM COEFICIENTES PARA AJUSTE DE 4 CTES
C
      DO 100 L=1,ITERMS
      DO 100 K=1,ITERMS
  100 COEFX(L)=COEFX(L)+ARRAY(L,K)*BETAX(K)
C
      DO  105 L=1,ITERMS
 105  TERMX(L)=COEFX(L)
      COEFX(1)=TERMX(3)
      COEFX(2)=TERMX(1)
      COEFX(3)=TERMX(2)
      COEFY(1)=TERMX(4)
      COEFY(2)=-TERMX(2)
      COEFY(3)=TERMX(1)

      RETURN


c
c     Ajuste de 1o, 2o ou 3o grau
c

 200  continue


C
C     CALCULA NUMERO DE TERMOS DO POL. DE GRAU=NGRAU
C

      NTERMS=1
      DO I=1,NGRAU
      NTERMS=NTERMS+I+1
      ENDDO
C
      IGRAU=NGRAU+1

      ITERMS=NTERMS
C
C


      DO 265 I=1,NCOMUM
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computando os termos para o polinomio em X
C
      ICONT=0
      DO 240 N=1,IGRAU
      DO 240 L=1,N
      K=N-L
      ICONT=ICONT+1
 240  TERMX(ICONT)=(X**K)*(Y**(L-1))

C
C     Computando AtA e AtB (para "X" e "Y")
C
      DO 260 L=1,ITERMS
      BETAX(L)=BETAX(L)+XG*TERMX(L)
      BETA(L)=BETA(L)+YG*TERMX(L)
      DO 260 K=L,ITERMS
  260 ALPHA(L,K)=ALPHA(L,K)+TERMX(K)*TERMX(L)
C
 265  CONTINUE

C
C     Preenchendo a parte triangular inferior da matriz simetrica AtA
c     para "X" e "Y"
C

      DO 270 L=1,ITERMS
      DO 270 K=L,ITERMS
  270 ALPHA(K,L)=ALPHA(L,K)
c
c
c     Preenchendo ARRAY=AtA para inversao (elementos normalizados
c     pela diagonal) para X ou Y
c
c

      DO 280 L=1,ITERMS
      DO 280 K=1,ITERMS
  280 ARRAY(L,K)=ALPHA(L,K)/DSQRT(ALPHA(L,L)*ALPHA(K,K))
C
C     Invertendo AtA para "X" ou " Y" 
C
      CALL MATINV (ITERMS,icofsp,array,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando os coeficientes do polinomio para "X"
C
      DO 290 L=1,ITERMS
      DO 290 K=1,ITERMS
  290 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHA(K,K)*ALPHA(L,L))

C
C     Guarda coeficiente X
C
      DO 300 L=1,ITERMS
      DO 300 K=1,ITERMS
  300 COEFX(L)=COEFX(L)+ARRAY(L,K)*BETAX(K)

C
C     Guarda coeficiente Y
C
      DO 310 L=1,ITERMS
      DO 310 K=1,ITERMS
  310 COEFY(L)=COEFY(L)+ARRAY(L,K)*BETA(K)

      RETURN

      END





c
c     Function pol
c
c     Calcula o valor de um polinomio bivariado em (x,y). Grau N
c     0 < N < 3
c
c

      double precision function pol (icofsp, x, y, coefis, ngrau)

      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION coefis(icofsp)

c
      mgrau=ngrau+1
      icont=0
      pol=0.d0
      do 1 n=1,mgrau
      do 1 l=1,n
      icont=icont+1
      k=n-l
    1 pol=pol+coefis(icont)*(x**k)*(y**(l-1))


      return
      end



c 
c
c     Subrotina posred
c
c
c     Reducao alfa e delta das posicoes (x,y) medidas em relacao a um
c     catalogo de referencia
c
c
c     Atualmente aplicada para reducao com o catalogo 2MASS e com o 
c     catalogo UCAC2
c
c     Assume-se que as poosicoes (RA,DEC) de catalogo estao na epoca
c     da observacao, no sistema ICRS (J2000)
c
c     Nessa versao nao ha peso.
c     Nessa versao nao ha termos de magnitude ou cor no polinomio de ajuste.
c
c     Ajustes possiveis: 4ctes, 1o, 2o, 3o graus completos, 2o+3dr,
c     2o+3dr+5dr, 3o+5dr.
c
c
c     - (RA,DEC) entram em graus, saem em graus
c     - (x,y)s entram em pixels
c     - erros e sigmas saem em segundos de arco (")
c     - coeficientes (e seus erros) saem em radianos por pixel,
c       rad por pixel**2 etc
c
c
c
c
c
c     Atualizacao: M. Assafin  03/Jul/2013
c
c
c



      subroutine posred (idiobs,icofsp,ireflex,rac,dec,id,ncat,racat,
     ?decat,nest,xob,yob,corte,ngrau,ngrau3,ngrau5,nstart,nfinal,ra,de,
     ?era,ede,alfsig,delsig,alfres,delres,coefx,coefy,ecoefx,ecoefy,
     ?itira,avam,dvam)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension id(idiobs),racat(idiobs),decat(idiobs),xob(idiobs),
     ?yob(idiobs),xp(idiobs),yp(idiobs),xest(idiobs),yest(idiobs)

      dimension ra(idiobs),de(idiobs),era(idiobs),ede(idiobs),
     ?coefx(icofsp),coefy(icofsp),ecoefx(icofsp),ecoefy(icofsp),
     ?alfres(idiobs),delres(idiobs),itira(idiobs),xsao(icofsp),
     ?ysao(icofsp),xrray(icofsp,icofsp),yrray(icofsp,icofsp),
     ?array(icofsp,icofsp)


      COMMON/A14/IERRO


      HMSGMS(I,J,A)=I+J/60.D0+A/3600.D0
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Dados auxiliares
C

      ncof=icofsp

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
      IERRO=0
      DNOVE=9.999D0
      d99=99.99d0

      EQUIN=2000.d0
      IZERO=0
      ZERO=0.D0
      pma=0.d0
      pmd=0.d0

c
c     Zerando vetores
c

      do i=1,idiobs

      xest(i)=0.d0
      yest(i)=0.d0
      xp(i)=0.d0
      yp(i)=0.d0
      ra(i)=0.d0
      de(i)=0.d0
      era(i)=0.d0
      ede(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0

      enddo


      do i=1,ncof

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0
      xsao(i)=0.d0
      ysao(i)=0.d0
      do j=1,ncof
      xrray(j,i)=0.d0
      yrray(j,i)=0.d0
      array(j,i)=0.d0
      enddo

      enddo



c
c     Recolhendo estrelas medidas comuns ao catalogo de referencia
c


      ntira=0
      grac=grarad*rac
      gdec=grarad*dec

      k=0

      do 10 i=1,nest

      if (id(i).eq.0) go to 10

      k=k+1

      if (itira(k).ne.0) then
      itira(k)=2
      ntira=ntira+1
      endif


      xest(k)=ireflex*xob(i)
      yest(k)=yob(i)

c
c     Projecao do catalogo de referencia no plano tangente
c


      bra=grarad*racat(id(i))
      bde=grarad*decat(id(i))
      d=DEXY(bra,bde,grac,gdec)

      xp(k)=xpad(bra,bde,grac)/d
      yp(k)=ypad(bra,bde,grac,gdec)/d


 10   continue

      kest=k

      nstart=kest-ntira


c
c     Ajuste do modelo polinomial entre (x,y) e (X,Y)
c


      call solucao (idiobs,icofsp,ngrau,ngrau3,ngrau5,kest,xest,yest,
     ?xp,yp,ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,
     ?itira,corte,xrray,yrray)


      if (ierro.eq.1) then

      write (*,*) 'L.S. solution crashed.'


c
c     Zerando vetores no erro
c

      do i=1,idiobs

      xest(i)=0.d0
      yest(i)=0.d0
      xp(i)=0.d0
      yp(i)=0.d0
      ra(i)=0.d0
      de(i)=0.d0
      era(i)=0.d0
      ede(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0
      itira(i)=0

      enddo


      do i=1,ncof

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0
      xsao(i)=0.d0
      ysao(i)=0.d0
      do j=1,ncof
      xrray(j,i)=0.d0
      yrray(j,i)=0.d0
      array(j,i)=0.d0
      enddo

      enddo

      return
      endif

c

      nfinal=kest-ntira




c
c     Determinacao do alfa e delta observado de cada estrela do campo
c
c     (RA,DEC) guardados em graus

c
c     Calculo do erro padrao em alfa e delta para cada estrela de campo
c

      XVAM=0.D0
      YVAM=0.D0
      XVAS=0.D0
      YVAS=0.D0


C
C     NGRAU=0 --> 4 ctes
C
      nterms=1
      do i=1,ngrau
      nterms=nterms+i+1
      enddo
C
      igrau=ngrau+1
      if (ngrau.eq.0) igrau=2
C
      if (ngrau.eq.0) then
      nterms=4
c     ngrau3=0
c     ngrau5=0
      endif
C
      iterms=nterms
      if (ngrau3.eq.3) iterms=iterms+1
      if (ngrau5.eq.5) iterms=iterms+1
C

      do 130 i=1,nest
      x=xob(i)*ireflex
      y=yob(i)

      ICONT=0
      POLX=0.D0
      POLY=0.D0
      DO 20 N=1,IGRAU
      DO 20 LL=1,N
      ICONT=ICONT+1
      K=N-LL
      POLX=POLX+COEFX(ICONT)*(X**K)*(Y**(LL-1))
   20 POLY=POLY+COEFY(ICONT)*(X**K)*(Y**(LL-1))
C
      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)
      ENDIF

      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF

      RA(I)=POLX
      DE(I)=POLY
C
C     Calcula o erro padrao em alfa para cada estrela
C
      ICONT=0
      if (ngrau.ne.0) then
      DO  30 N=1,IGRAU
      DO  30 LL=1,N
      K=N-LL
      ICONT=ICONT+1
   30 XSAO(ICONT)=(X**K)*(Y**(LL-1))
      IF (NGRAU3.EQ.3) THEN
       ICONT=ICONT+1
       XSAO(ICONT)=X*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
       ICONT=ICONT+1
       XSAO(ICONT)=X*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
c
      else
      xsao(1)=x
      xsao(2)=y
      xsao(3)=1.d0
      xsao(4)=0.d0
      endif
c
      DO  40 K=1,ITERMS
   40 YSAO(K)=0.D0
C
      DO  50 LL=1,ITERMS
      DO  50 K=1,ITERMS
   50 YSAO(LL)=YSAO(LL)+XRRAY(LL,K)*XSAO(K)

      DO  60 K=1,ITERMS
   60 era(I)=era(I)+XSAO(K)*YSAO(K)

      era(I)=alfsig*DSQRT(era(I))
      XVAM=XVAM+era(I)
      XVAS=XVAS+era(I)**2

C
C     Calcula o erro padrao em delta para cada estrela
C

      ICONT=0
      if (ngrau.ne.0) then
      DO  70 N=1,IGRAU
      DO  70 LL=1,N
      K=N-LL
      ICONT=ICONT+1
   70 YSAO(ICONT)=(X**K)*(Y**(LL-1))
      IF (NGRAU3.EQ.3) THEN
       ICONT=ICONT+1
       YSAO(ICONT)=Y*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
       ICONT=ICONT+1
       YSAO(ICONT)=Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
c
      else
      ysao(1)=y
      ysao(2)=-x
      ysao(3)=0.d0
      ysao(4)=1.d0
      endif
c
      DO  80 K=1,ITERMS
   80 XSAO(K)=0.D0
C
      DO  90 LL=1,ITERMS
      DO  90 K=1,ITERMS
   90 XSAO(LL)=XSAO(LL)+YRRAY(LL,K)*YSAO(K)

      DO 100 K=1,ITERMS
  100 ede(I)=ede(I)+YSAO(K)*XSAO(K)

      ede(I)=delsig*DSQRT(ede(I))
      YVAM=YVAM+ede(I)
      YVAS=YVAS+ede(I)**2

C
C     Calcula erro dos coeficientes para a solucao X
C

      IF (NGRAU.NE.0) THEN
      DO 110 K=1,ITERMS
 110  ECOEFX(K)=ALFSIG*DSQRT(XRRAY(K,K))
      ELSE
      ECOEFX(1)=ALFSIG*DSQRT(XRRAY(3,3))
      ECOEFX(2)=ALFSIG*DSQRT(XRRAY(1,1))
      ECOEFX(3)=ALFSIG*DSQRT(XRRAY(2,2))
      ENDIF
C
C     Calcula erro dos coeficientes para a solucao Y
C

      IF (NGRAU.NE.0) THEN
      DO 120 K=1,ITERMS
 120  ECOEFY(K)=DELSIG*DSQRT(YRRAY(K,K))
      ELSE
      ECOEFY(1)=DELSIG*DSQRT(YRRAY(4,4))
      ECOEFY(2)=DELSIG*DSQRT(YRRAY(2,2))
      ECOEFY(3)=DELSIG*DSQRT(YRRAY(1,1))
      ENDIF
C
C
 130  CONTINUE
C
C
      if (ngrau.eq.0) then
      iterms=3
      endif

C
C     Media dos erros padrao de todas as estrelas medidas
C

      EXMED=XVAM/NEST
      EYMED=YVAM/NEST
      XVAS=DSQRT((XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2)/(NEST-1.D0))
      YVAS=DSQRT((YVAS-2.D0*EYMED*YVAM+NEST*EYMED**2)/(NEST-1.D0))
      AVAM=EXMED
      DVAM=EYMED

c
c     RA e DEC em graus
c


      j=0

      do 140 i=1,nest

      x=ra(i)
      y=de(i)

      ra(i)=alff(x,y,grac,gdec)
      de(i)=deltt(ra(i),y,grac,gdec)

      ra(i)=ra(i)*radgra
      de(i)=de(i)*radgra

c     if (id(i).eq.0) go to 140

c     j=j+1

c     write(*,*) 'alfres delres itira = ',j,alfres(j),delres(j),itira(j) 



 140  continue


c
c     Debug
c


c     perc=100.d0*ntira/nstart
c     write (*,141) alfsig,delsig,nstart,nfinal,perc,avam,dvam
c141  format(1x,'alfsig delsig NI NF = ',2(1x,f6.3),2(1x,i4),1x,f6.2,
c    ?'%',2(1x,f6.3))




      return
      end



C
C     Subrotina solucao
C
C
C     Ajuste polinomial bivariado P=P(x,y) ate grau 3 e distorcao
C     radial a grau 5
C		    
C
C							    

      subroutine solucao (idiobs,icofsp,ngrau,ngrau3,ngrau5,nstart,xest,
     ?yest,xp,yp,ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,
     ?delres,itira,corte,xrray,yrray)




      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(icofsp),COEFY(icofsp),ALPHAX(icofsp,icofsp),
     ?ALPHAY(icofsp,icofsp),ARRAY(icofsp,icofsp),BETAX(icofsp),
     ?BETAY(icofsp),TERMX(icofsp),TERMY(icofsp),ITIRA(idiobs)

      DIMENSION XEST(idiobs),YEST(idiobs),XP(idiobs),YP(idiobs),
     ?XRRAY(icofsp,icofsp),YRRAY(icofsp,icofsp),ALFRES(idiobs),
     ?DELRES(idiobs)

      COMMON/A14/IERRO


      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
      DET  =1.D0
      IERRO=0

      alfsig=0.d0
      delsig=0.d0

c

      sigma=grarad*corte/3600d0

C
      NCOMUM=nstart
C
c
c     Computa no. de termos do polinomio
c     NGRAU=0 --> 4 constantes
c

      NTERMS=1
      DO I=1,NGRAU
      NTERMS=NTERMS+I+1
      ENDDO

C
      IGRAU=NGRAU+1
      IF (NGRAU.EQ.0) IGRAU=2
c     NTIRA=0
C


 5    LCONT=0

C
      RESMX=0.D0
      RES2X=0.D0
      RESMY=0.D0
      RES2Y=0.D0
C
      IF (NGRAU.EQ.0) THEN
      NTERMS=4
      NGRAU3=0
      NGRAU5=0
      ENDIF
C
      ITERMS=NTERMS
      IF (NGRAU3.EQ.3) ITERMS=ITERMS+1
      IF (NGRAU5.EQ.5) ITERMS=ITERMS+1

C
C     Checa No. de estrelas versus no. de coeficientes a ajustar
C

      IF (NGRAU.EQ.0) THEN
      iequa=2.d0*ncomum
      mtira=2.d0*ntira
      if ((iequa-mtira).LT.ITERMS) then
      IERRO=1
      RETURN
      endif
      else
      if ((ncomum-ntira).LT.ITERMS) then
      IERRO=1
      RETURN
      endif
      ENDIF

c
      DO 9 I=1,ITERMS
      BETAX(I) =0.D0
      BETAY(I) =0.D0
      COEFX(I) =0.D0
      COEFY(I) =0.D0
      TERMX(I) =0.D0
      TERMY(I) =0.D0
      DO 9 J=I,ITERMS
      ALPHAX(I,J)=0.D0
 9    ALPHAY(I,J)=0.D0


C
C     Montando equacoes de condicao para ajuste com 4 ctes
C

      IF (NGRAU.EQ.0) THEN

      DO 10 I=1,NCOMUM
      IF (ITIRA(I).NE.0) GO TO 10
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computa  termos dos coeficientes para AtB
C

      TERMX(1)=TERMX(1)+X**2+Y**2
      TERMX(2)=TERMX(2)+X
      TERMX(3)=TERMX(3)+Y
      TERMX(4)=TERMX(4)+1.D0
      BETAX(1)=BETAX(1)+X*XG+Y*YG
      BETAX(2)=BETAX(2)+Y*XG-X*YG
      BETAX(3)=BETAX(3)+XG
      BETAX(4)=BETAX(4)+YG
   10 CONTINUE

C
C     Preenchendo  AtA para ajuste com 4 ctes
C

      ALPHAX(1,1)=TERMX(1)
      ALPHAX(1,2)=0.D0
      ALPHAX(1,3)=TERMX(2)
      ALPHAX(1,4)=TERMX(3)
      ALPHAX(2,2)=TERMX(1)
      ALPHAX(2,3)=TERMX(3)
      ALPHAX(2,4)=-TERMX(2)
      ALPHAX(3,3)=TERMX(4)
      ALPHAX(3,4)=0.D0
      ALPHAX(4,4)=TERMX(4)
      DO 15 L=1,ITERMS
      DO 15 K=L,ITERMS
   15 ALPHAX(K,L) =ALPHAX(L,K)
C
      GO TO 75

      ENDIF


C
C     Montando equacoes de condicao para ajustes que NAO o de 4 ctes
C     para polinomio X


      DO 65 I=1,NCOMUM
      IF (ITIRA(I).NE.0) GO TO 65
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computando termos para o polinomio em X
C
      ICONT=0
      DO 40 N=1,IGRAU
      DO 40 L=1,N
      K=N-L
      ICONT=ICONT+1
      TERMX(ICONT)=(X**K)*(Y**(L-1))
   40 TERMY(ICONT)=(X**K)*(Y**(L-1))

      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      TERMX(ICONT)=X*(X**2+Y**2)
      TERMY(ICONT)=Y*(X**2+Y**2)
      ENDIF

      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      TERMX(ICONT)=X*(X**2+Y**2)*(X**2+Y**2)
      TERMY(ICONT)=Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
C
C     Computando AtA and AtB (para "X" e "Y")
C

      DO 60 L=1,ITERMS
      BETAX(L)=BETAX(L)+XG*TERMX(L)
      BETAY(L)=BETAY(L)+YG*TERMY(L)
      DO 60 K=L,ITERMS
      ALPHAX(L,K)=ALPHAX(L,K)+TERMX(K)*TERMX(L)
   60 ALPHAY(L,K)=ALPHAY(L,K)+TERMY(K)*TERMY(L)
C
   65 CONTINUE

C
C     Preenchendo parte triangular inferior da matriz simetrica AtA
c     (para "X" e "Y")
C

      DO 70 L=1,ITERMS
      DO 70 K=L,ITERMS
      ALPHAX(K,L)=ALPHAX(L,K)
   70 ALPHAY(K,L)=ALPHAY(L,K)

C
C     Preenchendo ARRAY=AtA para inversao (elementos normalizados pela
C     diagonal) para "X"
C

 75   CONTINUE
C
      DO 80 L=1,ITERMS
      DO 80 K=1,ITERMS
   80 ARRAY(L,K)=ALPHAX(L,K)/DSQRT(ALPHAX(L,L)*ALPHAX(K,K))
C
C     Invertendo AtA para "X"
C
      CALL MATINV (ITERMS,icofsp,array,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando coeficientes do polinomio para "X"
C
      DO 90 L=1,ITERMS
      DO 90 K=1,ITERMS
   90 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHAX(K,K)*ALPHAX(L,L))
C
C     Backup do array para "X"
C
      IF (NGRAU.EQ.0) GO TO 97
      DO 95 L=1,ITERMS
      DO 95 K=1,ITERMS
   95 XRRAY(L,K)=ARRAY(L,K)
C
 97   DO 100 L=1,ITERMS
      DO 100 K=1,ITERMS
  100 COEFX(L)=COEFX(L)+ARRAY(L,K)*BETAX(K)
C
C     Obtem coeficientes e backup do array para modelo de 4 Constantes
C
      IF (NGRAU.EQ.0) THEN
      DO  105 L=1,ITERMS
 105  TERMX(L)=COEFX(L)
      COEFX(1)=TERMX(3)
      COEFX(2)=TERMX(1)
      COEFX(3)=TERMX(2)
      COEFY(1)=TERMX(4)
      COEFY(2)=-TERMX(2)
      COEFY(3)=TERMX(1)
      do 107 l=1,iterms
      do 107 k=1,iterms
      xrray(l,k)=array(l,k)
 107  yrray(l,k)=array(l,k)
      GO TO 133
      ENDIF

C
C     Preenchendo ARRAY=AtA para inversao (elementos normalizados pela
C     diagonal) para "Y"
C

      DO 110 L=1,ITERMS
      DO 110 K=1,ITERMS
  110 ARRAY(L,K)=ALPHAY(L,K)/DSQRT(ALPHAY(L,L)*ALPHAY(K,K))
C
C     Invertendo AtA para "Y"
C

      CALL MATINV (ITERMS,icofsp,array,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando coeficientes do polinomio para "Y"
C


      DO 120 L=1,ITERMS
      DO 120 K=1,ITERMS
  120 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHAY(K,K)*ALPHAY(L,L))
C
C     Backup do array para "Y"
C

      DO 125 L=1,ITERMS
      DO 125 K=1,ITERMS
  125 YRRAY(L,K)=ARRAY(L,K)
C
      DO 130 L=1,ITERMS
      DO 130 K=1,ITERMS
  130 COEFY(L)=COEFY(L)+ARRAY(L,K)*BETAY(K)

c
c
c     Computa residuos. Estrelas com residuo em alfa ou delta mais alto
c     sao eliminadas uma a uma, ate que nenhuma possua (O-C) maior
c     que o valor definido pela variavel "corte"
c

 133  CONTINUE
C


      REMAXI=-1.D14

      DO 160 I=1,NCOMUM
      X=XEST(I)
      Y=YEST(I)
      xg=xp(i)
      yg=yp(i)
      ICONT=0
      POLX=0.D0
      POLY=0.D0
      DO 140 N=1,IGRAU
      DO 140 L=1,N
      ICONT=ICONT+1
      K=N-L
      POLX=POLX+COEFX(ICONT)*(X**K)*(Y**(L-1))
  140 POLY=POLY+COEFY(ICONT)*(X**K)*(Y**(L-1))
C
      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF

C
C     Computa (O-C) para R.A. e  Dec
C

      xx=alff(polx,poly,grac,gdec)
      yy=deltt(xx,poly,grac,gdec)

      xxr=alff(xg,yg,grac,gdec)
      yyr=deltt(xxr,yg,grac,gdec)


      AUX=(XX-XXR)*DCOS(YYR)
      AUY=YY-YYR
      ALFRES(I)=AUX*radgra*3600d0
      DELRES(I)=AUY*radgra*3600d0


C
      IF (ITIRA(I).NE.0) GO TO 160
      IF ((DABS(AUX).GT.REMAXI).OR.(DABS(AUY).GT.REMAXI)) THEN
      IFORA =I
      REMAXI=DMAX1(DABS(AUX),DABS(AUY))
      ENDIF

 157  RESMX=RESMX+ALFRES(I)
      RES2X=RES2X+ALFRES(I)**2
      RESMY=RESMY+DELRES(I)
      RES2Y=RES2Y+DELRES(I)**2

      LCONT=LCONT+1

C
  160 CONTINUE
C
C     Atingido numero minimo de estrelas!
C
      iwar=0
      IF (NGRAU.EQ.0) THEN
      iequa=2.d0*LCONT
      if (iequa.eq.4) iwar=1
      ELSE
      if (LCONT.eq.ITERMS) iwar=1
      ENDIF
C
      IF (IWAR.EQ.1) THEN
      alfsig=0.D0
      delsig=0.D0
      RETURN
      ENDIF


C
C     Computa media e desvio padrao dos (O-C)s
C

      XMED=RESMX/LCONT
      YMED=RESMY/LCONT
C
      IF (NGRAU.EQ.0) THEN
      alfsig=DSQRT((RES2X-2.D0*XMED*RESMX+LCONT*XMED**2)/(LCONT-2.D0))
      delsig=DSQRT((RES2Y-2.D0*YMED*RESMY+LCONT*YMED**2)/(LCONT-2.D0))
      ELSE
      alfsig=DSQRT((RES2X-2.D0*XMED*RESMX+LCONT*XMED**2)/(LCONT-ITERMS))
      delsig=DSQRT((RES2Y-2.D0*YMED*RESMY+LCONT*YMED**2)/(LCONT-ITERMS))
      ENDIF
C
c
c     Atingido o corte !
c

      if (remaxi.lt.sigma) return

c
c     Corte nao atingido, prosseguir com eliminacao de estrelas
c

      NTIRA=NTIRA+1
      ITIRA(IFORA)=1
      GO TO 5


      RETURN
      END


c   
c     Subrotina estat
c
c     Faz a estatistica dos alvos, obtendo-se os O-R, O-E ...
c     Serve para a reducao com 2MASS e com UCAC2.
c
c     Alem do O-R, eh escrito no arquivo de saida todos os 
c     parametros obtidos para o objeto tal qual escrito no
c     arquivo tipo xy para esse objeto.
c
c
c     box -> box de erro para identificacao, em segundos de
c                                                      arco
c
c     boxepo -> box de erro no tempo para identificacao de alvo
c               com coordenada tempo-dependente na imagem
c

      subroutine estat (box,ialvos,input,itotal)


      IMPLICIT REAL *8 (A-H,O-Z)

      character*150 infits 
      character*50 ialvos,input,itotal
      character*20 ichfil,mchobj,iobalv
      character*1 isig,menos

      data menos/'-'/

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c
c     Box de erro para instante de observacao, em segundos
c    (default box=1s)
c
      boxepo=1.d0

c
c     Box de erro em data juliana (fracao de dia)
c

      boxepo=boxepo/(3600.d0*24.d0)

c

      ipegou=0

c
c     Determina backspace de arquivos (g77, gfortran, etc)
c

      call backsp (1,nbac,7)

c

 40   format(a50)


      write (*,*)
      write (*,*)'offsets (RA,DE), Sigma(RA,DE), Ncat, Gauss_error(x,y),
     ?date, exptime, filter, target'
      write (*,*)

c

      open (7,file=ialvos,status='old',err=200)

c
      open (2,file=itotal)

c

 3    read (2,5,end=6) isig
 5    format(a1)
      go to 3
 6    call backsp (2,nbac,2)

c


 100  read (7,101,end=200) iah,iam,as,isig,idg,idm,ds,datalv,
     ?iobalv
 101  format(1x,i2,1x,i2,1x,f9.6,1x,a1,i2,1x,i2,1x,f8.5,1x,f16.8,
     ?1x,a20)

      rafat=hmsgms(iah,iam,as)
      defat=hmsgms(idg,idm,ds)
      if (isig.eq.menos) defat=-defat


 30   continue


      open (1,file=input)


c
c     Checa data juliana do alvo e das medidas para identificacao
c

 18   continue
      read (1,10,err=18,end=12) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny



c
c     Alvos com data juliana "negativa" sao alvos com coordenada
c     fixa, podem ser procurados em qq imagem desde que dentro da
c     box de alfa e de delta
c
c     Alvos com data juliana normal sao alvos com coordenada
c     dependente do tempo; a procura do alvo so ocorre na imagem
c     correspondente a da data juliana fornecida, dentro de 1s de
c     erro.
c


      if (datalv.gt.0.d0) then

      dtemp=dabs(datalv-dj)

      if (dtemp.gt.boxepo) go to 12

      endif

      rewind (1)

c
c     Aqui tudo ok para busca de alvo em relacao a data juliana
c     (seja alvo fixo ou movel)
c



 19   continue
      read (1,10,err=19,end=12) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


 10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))




c
c     Box de erro em alfa e delta para identificacao de alvo
c


      dx=(ra-rafat)*dcos(defat*grarad)*3600.d0*15.d0
      dy=(de-defat)*3600.d0

      if (dabs(dx).gt.box) go to 20
      if (dabs(dy).gt.box) go to 20

      ipegou=1


      write (*,16) dx,dy,alfsiu,delsiu,nfinau,ex,ey,iuth,
     ?iutm,sut,iutano,iutmes,iutdia,iexps,ichfil,iobalv

 16   format(4(1x,f7.3),2x,i5,2x,2(1x,f7.3),2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,2x,i4,3x,a20,2x,a20)



      write (2,11) dx,dy,xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,iobalv,nx,ny


 11   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))


 12   close (1)

      go to 100

c

 20   go to 19

c


 200  close (2)
      close (7)

      if (ipegou.eq.0) then
      write (*,*) ' Target not identified. '
      endif

      return
      end



c   
c     Subrotina tempo
c
c     Devolve o tempo total decorrido desde o inicio da execucao e o tempo
c     parcial decorrido entre a ultima e a chamada atual, em segundos.
c     
c     O tempo eh a soma do tempo gasto pela CPU no programa mais o tempo
c     gasto pelo sistema na execucao do programa
c
c
c     tempoi = tempo total decorrido ate a ultima chamada, anterior a atual 
c     tempot = tempo total decorrido ate esta chamada
c     tempop = tempo parcial decorrido entre a ultima e a atual chamada
c

      subroutine tempo (tempoi,tempot,tempop)

      implicit real*8 (a-h,o-z)

      real  time(2)

c

      tempop=dtime(time)

      tempot=tempoi+tempop

      return

      end



c   
c     Subrotina backsp
c
c     Efetua o numero correto de backspaces no arquivo aberto "L", para
c     "recuar uma linha" no arquivo "L".
c
c     O numero de backspaces depende do fortran utilizado (gfortran,
c     g77, etc ...)
c
c
c     key=1 : determina o numero correto de backspaces a ser executado
c
c     key=2 : executa o numero correto de backspaces no arquivo aberto "L"
c
c     nbac = numero correto de backspaces a ser executado, para recuar
c            uma linha no arquivo
c
c     L    = unidade do arquivo aberto para execucao de backspace
c

      subroutine backsp (key,nbac,L)

      implicit real*8 (a-h,o-z)

      character*20 imaux
      character*4 jmaux
      character*9 sista
      character*29 systa


c
c     Key=2, executar backspace
c


      if (key.eq.2) then

      do i=1,nbac
      backspace L
      enddo

      return

      endif

c
c     Key=1, determinar numero correto de backspaces a ser executado,
c     para recuar uma linha em um arquivo aberto
c


      sista='rm -f -r '

      imaux=''

      imaux(1:16)='PRAIA_backsp.aux'


      do 1 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux(17:20)=jmaux(1:4)

      open(L,file=imaux,status='old',err=2)
      close (L)
 1    continue

 2    close (L)

      open(L,file=imaux)

      systa=sista//imaux

c

      do i=1,10
      write (L,*) i
      enddo

      close (L)

      open(L,file=imaux)
      
      do i=1,15
      read (L,*,end=10)
      enddo

 10   do i=1,2
      backspace L
      enddo
      
 
      read (L,*) nbac

      close (L)

      nbac=nbac-9

c


      call system (systa)

      return

      end







c   
c
c     subroutine desvio
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

      implicit real *8 (a-h,o-z)

c

      zero=0.d0
      dnove=99.999d0
      dneg=-1.d0

c

      if (nest.le.0) then
      xvam=zero
      xvas=dneg
      return
      endif


c

      exmed=xvam/nest

c

      if (nest.eq.1) then

      xvas=dneg

      return

      endif

c

      if (nest.eq.2) then


      xvas=dsqrt(dabs(2.d0*xvas-xvam*xvam))/2.d0

      xvam=exmed


      return

      endif

c

      raiz=xvas-2.d0*exmed*xvam+nest*exmed*exmed

      if (raiz.lt.0.d0) then

      xvas=dneg

      else

      xvas=dsqrt(raiz/(nest-1.d0))

      endif
c

      xvam=exmed

      return

      end


