c
c     PRAIA_targets_skybot
c
c
c
c     Purpose
c
c
c     Extracts ephemeris for all known solar system objects within the FOV
c     of each individual field. 
c
c     For that, it queries the SKYBOT service at vo.imcce.fr.
c     
c
c     Last update:  M. Assafin - 22/Sep/2011
c
c



      implicit real*8 (a-h,o-z)

      parameter(idim=1000000)

      dimension alpha(idim),delta(idim),dju(idim),mark(idim)

      character*50 lista1,lista2,infits
      character*20 ichobj,ichfil,numero
      character*20 mchobj(idim),nume(idim)
      character*1  isig,ivelha,ler(200),ipipe
      character*10 iver
      character*14 ifrmra,ifrmr3,ifrmr2,ifrmr1
      character*18 ifrmde,ifrmd2,ifrmd1
      character*12 if50,if40,if30,if20,ifrs
      character*11 if4,if3,if2,if1,ifloc
      character*200 iaux

      character*20 imaux1,imaux2,imaux3
      character*4 jmaux
      character*9 sista
      character*29 systa1,systa2,systa3
      character*29 chama


c
c     Data values
c

      data ivelha/'#'/
      data ipipe/'|'/

c
c     Functions
c


      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

c
c     Auxiliary data
c

      pi    = 0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi
c


      ifrmr3="('RA=',f16.12)"
      ifrmr2="('RA=',f15.12)"
      ifrmr1="('RA=',f14.12)"

      ifrmd2="('DEC=',a1,f15.12)"
      ifrmd1="('DEC=',a1,f14.12)"

      if50="('RS=',f5.0)"
      if40="('RS=',f4.0)"
      if30="('RS=',f3.0)"
      if20="('RS=',f2.0)"

      if4="('LOC=',i4)"
      if3="('LOC=',i3)"
      if2="('LOC=',i2)"
      if1="('LOC=',i1)"

c

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'PRAIA automatic ephemeris extraction by SKYBOT')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '

c

 100  continue


c
c     Mounts auxiliary file name     
c


      sista='rm -f -r '

      imaux1=''

      imaux1(1:16)='PRAIA_skybot.aux'


      do 120 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux1(17:20)=jmaux(1:4)

      open(99,file=imaux1,status='old',err=130)
      close (99)
 120  continue

 130  close (99)

      systa1=sista//imaux1



c
c     Mounting .sh file name for bash query
c


      imaux2=''

      imaux2(1:13)='PRAIA_skybot_'
      imaux2(18:20)='.sh'

      do 140 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux2(14:17)=jmaux(1:4)

      open(99,file=imaux2,status='old',err=150)
      close (99)
 140  continue

 150  close (99)

      systa2=sista//imaux2



c
c     Mounts skybot query output file name 
c


      imaux3=''

      imaux3(1:16)='PRAIA_skybot.dat'


      do 160 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux3(17:20)=jmaux(1:4)

      open(99,file=imaux3,status='old',err=170)
      close (99)
 160  continue

 170  close (99)

      systa3=sista//imaux3



c
c     Opens auxiliary file
c


      open (12,file=imaux1)


c
c     Reads input data
c

c     open (1,file='PRAIA_targets_skybot_20_01.dat')

      read (*,3) lista1

      if (lista1(1:1).eq.' ') go to 110

      read (*,3) lista2

      nstart=1

      read (*,*) pixel
      read (*,*) iauloc
      read (*,*) tep
      read (*,*) eboxy
      read (*,*) etime

      eboxy=eboxy/3600.d0
      etime=etime/86400.d0

      read (*,*)

 3    format(a50)

c     close (1)


c
      write (*,5)
 5    format('Input file with header extraction data         -> ',$)
      write(*,3) lista1
c
      write (*,6)
 6    format('Output file with targets ephemeris from SKYBOT -> ',$)
      write(*,3) lista2

c

      write (*,*)
      write (*,*)


c
c     Lendo Efemerides e escrevendo no formato para reducao (RA,DEC)
c

      open (3,file=lista1)
      
      m=0

c

      do i=1,10000000
      read (3,*,end=7)
      enddo

 7    rewind (3)

      images=i-1


c
c     Automatic extraction begins
c

      jj=nstart-1

      n=nstart-1
      do i=1,n
      read(3,*)
      enddo

c

 10   read (3,20,end=70) iah,iam,sa,isig,idg,idm,ds,iuth,iutm,sut,
     ?iutano,iutmes,iutdia,djm,dj,iexps,ichfil,infits,ichobj,nx,ny
 20   format(1x,i2,1x,i2,1x,f7.4,1x,a1,i2,1x,i2,1x,f6.3,2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,f16.8,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))


c
c     N-th field
c

      jj=jj+1

      write (*,21) jj,images,infits
 21   format('Image ',i5,' of ',i5,' -  Fits file ',a50)

c
c     Organizing data for skybot query
c

      ra =15.d0*hmsgms(iah,iam,sa)
      dec=hmsgms(idg,idm,ds)
      if (isig.eq.'-') dec=-dec

      nnn=max0(nx,ny)
      epoch=dj
      bs=pixel*2.d0*nnn      
      loc=iauloc


c

      open (8,file=imaux2)

c
c     Epoch
c

      write (8,22) epoch
 22   format('EPOCH=',f16.8)

c
c     RA
c

      ifrmra=ifrmr3
      if (ra.lt.100.d0) ifrmra=ifrmr2
      if (ra.lt.10.d0) ifrmra=ifrmr1

      write (8,ifrmra) ra

c
c     Dec
c

      if (dec.lt.0.d0) then
      isig='-'
      else
      isig='+'
      endif
      
      dec=dabs(dec)

      ifrmde=ifrmd2
      if (dec.lt.10.d0) ifrmde=ifrmd1

      write (8,ifrmde) isig,dec

c
c     FOV (box) in arcsecs
c


      ifrs=if50
      if (bs.lt.1000.d0) ifrs=if40      
      if (bs.lt.100.d0) ifrs=if30      
      if (bs.lt.10.d0) ifrs=if20      

      write (8,ifrs) bs

c
c     IAU observatory code
c

      ifloc=if4
      if (loc.lt.1000.d0) ifloc=if3      
      if (loc.lt.100.d0) ifloc=if2      
      if (loc.lt.10.d0) ifloc=if1      

      write (8,ifloc) loc


      write (8,23) imaux3
 23   format('OUTPUT_FILE=',a20)

c
c     Mounts Skybot query
c

      write (8,24)

 24   format('QUERY="http://vo.imcce.fr/webservices/skybot/skybotconesea
     ?rch_query.php?-ep="$EPOCH"&-ra="$RA"&-dec="$DEC"&-rs="$RS"&-mime=t
     ?ext&-out=object&-loc="$LOC"&-filter=000"')


c24   format ('QUERY="http://www.imcce.fr/webservices/skybot/skybot_
c    ?query.php?-ep="$EPOCH"&-ra="$RA"&-dec="$DEC"&-rs="$BS"&-mime=text&
c    ?-out=object&-loc="$LOC"&-filter=000"')

c
c     Mounts skybot call
c

      write (8,25)
 25   format('wget -q $QUERY -O  $OUTPUT_FILE')

c

      close (8)

c
c     Skybot query
c

      write (chama,26) imaux2
 26   format('chmod +x ',a20)

      call system(chama)

      chama=''
      write (chama,27) imaux2
 27   format('./',a20)
 

 28   call system (chama)

c     

      tempoi=0.d0
      tim=0.d0
      call tempo (tempoi,tempot,tempop)

c
c     Reads results from skybot
c

 30   continue

      open (9,file=imaux3,err=31)
      go to 33
 31   close (9)
 32   go to 30

c
 33   iver=''
      read (9,35,err=31,end=31) iver
 35   format(a10)

c
c     Checks if the output file arrived
c

      if (iver(1:1).eq.' ') go to 31

      if (iver(1:1).eq.ivelha .or. iver.eq.'Flag  : -1') then
      close (9)
      write (*,38) tep
 38   format('SKYBOT server time delay ... ',f5.1,' secs') 
 39   call tempo (tempoi,tempot,tempop)
      tim=tim+tempop
      if (tim.lt.tep) go to 39
      go to 28
      endif


c
c     Skips two more header lines
c

      read (9,35,err=31,end=31) iver
      if (iver.eq.ivelha) go to 31

      read (9,35,err=31,end=31) iver
      if (iver.eq.ivelha) go to 31


c
c     Loads ephemeris of all found objects and stores them in PRAIA format
c

 40   continue

      do i=1,200
      ler(i)=''
      enddo

      read (9,45,err=65,end=65) (ler(i),i=1,200) 
 45   format(200a1)

c
c     Picks up object names and numbers
c

      do i=1,200
      if (ler(i).eq.ipipe) go to 47
      enddo

      go to 65

 47   i1=i+1

      numero=''
      write (numero,*) (ler(j),j=1,i-1)

      do i=i1,200
      if (ler(i).eq.ipipe) go to 49
      enddo

 49   i2=i-1
      i3=i+1

      ichobj=''

      kk=0
      do i=i1,i2
      kk=kk+1
      ichobj(kk:kk)=ler(i)
      enddo

c
c     Picks up RAs 
c

      do i=i3,200
      if (ler(i).eq.ipipe) go to 51
      enddo
     
 51   i4=i-1
      i5=i+1

      iaux=''
      j=0
      do i=i3,i4
      j=j+1
      iaux(j:j)=ler(i)
      enddo

      rewind (12)
      write (12,53) iaux
 53   format(a200)

c
c     Picks up DECs
c
      do i=i5,200
      if (ler(i).eq.ipipe) go to 55
      enddo
     
 55   i6=i-1

      iaux=''
      j=0
      do i=i5,i6
      j=j+1
      iaux(j:j)=ler(i)
      enddo

      write (12,53) iaux

      rewind (12)

      read (12,*) ra
      read (12,*) dec

      rewind (12)

c
c     Stores ephemeris in the output file
c

      m=m+1

      alpha(m)=ra
      delta(m)=dec
      dju(m)=dj
      mchobj(m)=ichobj
      nume(m)=numero

      go to 40

c

 65   close (9)

      call system (systa3)
      call system (systa2)


      go to 10

c

 70   close (12)

      close (3)

c

      write (*,*)
      call system (systa1)
      write (*,*)


c
c     Checks for double entries
c

      do i=1,m
      mark(i)=0
      enddo

      do i=1,m-1
      do 73 j=i+1,m

      if (mark(j).ne.0) go to 73

      dift=dabs(dju(j)-dju(i))

      if (dift.gt.etime) go to 73

      dy=dabs(delta(j)-delta(i))

      if (dy.gt.eboxy) go to 73

      dx=dabs(15.d0*(alpha(j)-alpha(i))*dcos(grarad*delta(i)))

      if (dy.gt.eboxy) go to 73

      mark(j)=1

 73   continue
      enddo


c
c     Writes ephemeris in the output file
c

      open (7,file=lista2)


      do 80 i=1,m

      if (mark(i).ne.0) go to 80

      ra=alpha(i)
      dec=delta(i)

      ichobj=mchobj(i)
      numero=nume(i)

      dj=dju(i)

      iah=ra
      am=(ra-iah)*60.d0
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

c

      write(7,75) iah,iam,sa,isig,idg,idm,ds,dj,ichobj,numero
 75   format(1x,i2,1x,i2,1x,f9.6,1x,a1,i2,1x,i2,1x,f8.5,1x,f16.8,
     ?2(1x,a20))


 80   continue

      close (7)


c
c     Go to next block
c

      write (*,*)

      go to 100

c

 110  continue

c
      
      close (12)
      call system (systa1)

      write (*,*) ' Finished ok.'
      write (*,*)
      write (*,*)


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

