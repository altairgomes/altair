c
c     Programa PRAIA_findscale
c
c     PROPOSITO
c
c     Determina a escala de pixel e erro associado das imagens fornecidas,
c     uma a uma.
c 
c
c     Dado um conjunto de imagens fits, identifica objetos no campo,
c     mede os (x,y), identifica as estrelas do catalogo de referencia
c     (2MASS), e estima a escala de pixel e seu erro. 
c
c     Os algoritmos de identificacao, medida e reconhecimento de estrelas
c     de catalogo sao os mesmos empregados na task PRAIA. A identificacao
c     catalogo-medida e' feita com ajuste de 4 ctes par a par, com e sem
c     reflexao em X. Sao usadas as N (ex.: N=30) estrelas mais brilhantes
c     medidas e as M estrelas mais brilhantes de catalogo (ex:N=100) para
c     a identificacao. 
c
c
c
c      Last update: Marcelo Assafin - 08 Outubro 2012
c   
c
c


      IMPLICIT REAL *8 (A-H,O-Z)
      parameter (stdin=5)

      integer*4 imagem(5001,5001)
      real*4 PIXMAT(5001,5001),maximo,mazi
      real*4 pixel(5001,5001)
      integer*2 bitpix,bitpyx,betpix


      dimension contag(50000),histo(30),ico(30),ior(50000),nval(50000)

      dimension xob(50000),yob(50000),ilado(50000),seng(50000),
     ?iflag(50000),altu(50000),ialtu(50000)

      dimension exgcc(50000),eygcc(50000),sgcc(50000),fgcc(50000)

      dimension xold(50000),yold(50000),altold(50000),ialtol(50000)

      dimension ra2ma(50000),de2ma(50000),era2ma(50000),ede2ma(50000),
     ?dmgj(50000),dmgh(50000),dmgk(50000),emgj(50000),emgh(50000),
     ?emgk(50000),xra2ma(50000),yde2ma(50000),id2ma(50000),ddj2(50000),
     ?iduc2(50000),xrauc2(50000),ydeuc2(50000),iduc4(50000),
     ?xrauc4(50000),ydeuc4(50000)

      dimension rauc2(50000),deuc2(50000),erauc2(50000),edeuc2(50000),
     ?pmra(50000),pmde(50000),epmra(50000),epmde(50000),udmgj(50000),
     ?udmgh(50000),udmgk(50000),udmg(50000),cudmg(50000)


      dimension rauc4(50000),deuc4(50000),erauc4(50000),edeuc4(50000),
     ?pmra4(50000),pmde4(50000),epmra4(50000),epmde4(50000),
     ?udmgj4(50000),udmgh4(50000),udmgk4(50000),udmg4(50000),
     ?cudmg4(50000)


      dimension era2(50000),ede2(50000),alfre2(50000),delre2(50000),
     ?coefx2(21),coefy2(21),ecofx2(21),ecofy2(21),itira2(50000)


      dimension coefxr(21),coefyr(21),ecofxr(21),ecofyr(21)




      dimension erau(50000),edeu(50000),alfreu(50000),delreu(50000),
     ?coefxu(21),coefyu(21),ecofxu(21),ecofyu(21),itirau(50000)

      dimension era4(50000),ede4(50000),alfre4(50000),delre4(50000),
     ?coefx4(21),coefy4(21),ecofx4(21),ecofy4(21),itira4(50000)

      dimension xmed(50000),ymed(50000),xpa(50000),ypa(50000),
     ?cocomx(21),cocomy(21),iuc2ma(50000),cudmg2(50000),iuc4ma(50000)

      dimension adx(10000),ady(10000),coordx(10000),coordy(10000)


      dimension inuu2(288,240),ir1u2(288,240),ir2u2(288,240)
      dimension inuu4(900,1440),ir1u4(900,1440),ir2u4(900,1440)


      character*150 infits,imfits,names(20000),ids9,ibadpx,kbadpx


      character*50 lista,centro

      character*50 redred,ifdp

      character *50 mraiz

      character*50 subf

      character*1   menos,iver,nome(50),isig,ibrac
      character*69 iobold,ichobj,ihname(20000),mchobj
      character*20 ichfil

      character*200 linha


      COMMON/A3/IMAGEM
      COMMON/A4/PIXMAT
      COMMON/A14/IERRO

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
C     Daddos iniciais
C
c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

      dj2000=2451545.d0


      idiobs=50000

      icofsp=21
  
      ng=10
      
      idin=150
      idin50=50
      idin2=200

      iu2z=288
      iu2i=240

      iu4z=900
      iu4i=1440


      ipmax=5001


      tt0=2.8d-3

      izero=0
      zero=0.d0

      d99=99.999d0

      mpes=50000
      nhist=30

      nhisfw=10

      jfdp=10000

c

      do i=1,mpes
      contag(i)=0.
      enddo

      do i=1,idiobs
      cudmg(i)=0.d0
      ior(i)=0
      enddo



c
c     Zerando imagem, pixmat e matpix
c

      do i=1,ipmax
      do j=1,ipmax
      pixmat(j,i)=0.
      enddo
      enddo

c
c     Define numero correto de backspaces para recuar uma linha em arquivo
c     aberto, de acordo com a versao fortran utilizada.
c

      call backsp (1,nbac,91)

   

c
c     Lendo dados de entrada
c

c     open (1,file='PRAIA_findscale_20_03.dat')

      read (*,3) mraiz
      read (*,3) centro

      read (*,3) ifdp
      read (*,*) kfdp  
      read (*,3) kbadpx

      read (*,3) redred

      read (*,*) akey,idegx,iminux,asegx,idegy,iminuy,asegy


      read (*,*) scala1,scala2
      read (*,*) ecala

      ecala=ecala/2.d0


      read (*,*) mazi
      read (*,*) vmin
      read (*,*) ipflag
      read (*,*) bscale
      read (*,*) bzero
      read (*,*) bitpyx
      read (*,*) kswap
      read (*,*) ngrauf
      read (*,*) malisa
      read (*,*) fatceu
      read (*,*) fmin, fmax
      read (*,*) nbcat
      read (*,*) nbmed
      read (*,*) erpix
      read (*,*) pcort2
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
 1    format (23x,'PRAIA - find scale setup')
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

      maximo=mazi

c
      
      if (inicio.eq.0 .and. iultmo.eq.0) then

      open (1,file=centro)

      i=0
 4    read (1,*,end=2028)
      i=i+1
      go to 4

 2028 close (1)

      inicio=1
      iultmo=i
      

      endif


c


      scala=(scala2-scala1)/2.d0

      fmin=fmin/scala
      fmax=fmax/scala


c
c     Le nomes de imagens e nomes de objetos
c



      open (3,file=centro)


      do i=1,inicio-1
      read(3,*,end=200)
      enddo

      do 190 i=inicio,iultmo

      names(i)=''
      ihname(i)=''


      read(3,402,err=194,end=200) iah,iam,sa,isig,idg,idm,ds,iuth,
     ?iutm,sut,iutano,iutmes,iutdia,djm,dj,iexps,ichfil,names(i),mchobj

      ihname(i)=mchobj


 190  continue

      go to 200

c

 194  write (*,195) names(i),mchobj
 195  format(1x,'Reading error on image list: ',a50,1x,a69)

c

 200  close (3)



c
c     Inicia tratamento das imagens
c



      scl=0.d0
      scl2=0.d0
      nscl=0
 

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


c


c
c     Abre arquivo de resultados
c

      open (97,file=redred)

c

      do 60 lllll=inicio,iultmo


      infits=''
      infits=names(lllll)


c     

 91   format(a150)
 11   format(a50)
 12   format(50a1)
 92   format(150a1)
 93   format(a200)


c
c    Nome do arquivo de bad pixels
c

      ibadpx=''

      ibadpx(1:kkii)=infits(1:kkii)
      ibadpx(kkii+1:kkii+4)='.bpx'


c
c    Nome do arquivo de boxes para ds9 (formato SAOimage antigo)
c

      ids9=''

      ids9(1:kkii)=infits(1:kkii)
      ids9(kkii+1:kkii+4)='.reg'

c


c

      write (*,*)
      write (*,*)

      write (*,15) lllll,iultmo,infits
 15   format (1x,'Proccessing field ',i5,' of ',i5,': image = ', a50)

      write (*,*)

c
c     Estimando tempo de processamento
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
 16   format(1x,'Estimated time left for end of PRAIA reductions: ',
     ?i3,'days ',i2,'hs ',i2,'m ',f4.1,'s')  
      write (*,*)      

      endif

c
c     Inicializando rotinas de estimativa de tempo de execucao da
c     reducao das imagens com o PRAIA
c

      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)

      if (lllll.ne.inicio) then

      seg=tt0*3600.d0

      write (*,*)      
      write (*,17) seg
 17   format(1x,'Estimated time for field stars identification and (x,y)
     ? Gaussian fits: ',f5.1,'seconds')  
      write (*,*)      

      endif



c
c     Field Distortion Pattern. Corrige distorcao de campo aplicando
c     offsets nas medidas (x,y).
c
 
      nfdp=0

      open (23,file=ifdp,status='old',err=280)

      do i=1,jfdp
      read (23,*,end=275) adx(i),ady(i),coordx(i),coordy(i)
      enddo

 275  nfdp=i-1

 280  close (23)



C
C     Lendo a imagem fits do campo, carrega matriz original, reconhece
C     tamanho em pixels em x e em y, e o valor maximo de contagem
c
c     nheads e' o numero de headers de 2880 bytes no header
c     original nas imagens (LNA = 2, p.ex!)
C

      bitpix=bitpyx

      write (*,*)

      maximo=mazi
      call refits (infits,nx,ny,nheads,ichobj,ipflag,bscale,bzero,
     ?kswap,iswap,bitpix)
      ichobj=ihname(lllll)


      write (*,*)


c
c     Marca pixels anormalmente negativos (normalmente saturacao no CCD)
c

      do i=1,ny
      do j=1,nx
      if (pixmat(j,i).lt.vmin) pixmat(j,i)=mazi+1.d0
      enddo
      enddo

c
c     Exclui pixels da matriz marcados como "bad pixels" do arquivo de
c     bad pixels apenas da imagem corrente (se existir)
c

      open (23,file=ibadpx,status='old',err=270)
 
 250  read (23,*,end=270) ix1,ix2,iy1,iy2
 
      do i=iy1,iy2
      do j=ix1,ix2
      pixmat(j,i)=mazi+1.
      enddo
      enddo
      
      go to 250
 
 270  close (23)


c
c
c     Exclui pixels da matriz marcados como "bad pixels" do arquivo de
c     bad pixels valido para todas as imagens tratadas (se existir)


      open (23,file=kbadpx,status='old',err=271)
 
 251  read (23,*,end=271) ix1,ix2,iy1,iy2
 
      do i=iy1,iy2
      do j=ix1,ix2
      pixmat(j,i)=mazi+1.
      enddo
      enddo
      
      go to 251
 
 271  close (23)




c
c     Flattens sky background: 1rst step with star contamination
c



      if (ngrauf.ge.1 .and. ngrauf.le.15) then

      write (*,*)
      write (*,*) 'Flattening sky background. Please wait ...'
      write (*,*)

      fcmin=vmin
      fcmax=mazi

      call flatsk (nx,ny,ngrauf,mazi,fcmin,fcmax)

      endif




c
c     Computes typical sky background in the assumption of a flat sky
c


      call skyb (nx,ny,mazi,malisa,vmin,fatceu,ceu,sigceu,izmax,threso)

c

      if (ngrauf.ge.1 .and. ngrauf.le.15) then
      write (*,*) 'Flattening sky background. 1rst step done.'
      endif
      write (*,*) 'sky background mode, freq, sky sigma = ',ceu,izmax,
     ?sigceu


c
c     If sky background computation fails ...
c


      if (izmax.eq.0) then
      write (*,*) 'Sky background statistics failed.'
      go to 60
      endif

c

      write (*,*) 'sky background threshold for objects ID = ',threso
      write (*,*)
      write (*,*)
c

      if (threso.ge.mazi-2.d0) then
      write (*,*) 'Sky background statistics failed.'
      go to 60
      endif



c
c     Now refines sky background flattening by eliminating star
c     contamination. Only pixels really associated to the sky background
c     are fitted.
c


      if (ngrauf.ge.1 .and. ngrauf.le.15) then

      fcmin=ceu-2.5d0*sigceu
      fcmax=ceu+2.5d0*sigceu

      call flatsk (nx,ny,ngrauf,mazi,fcmin,fcmax)

c
c     Re-computes typical sky background in the assumption of a flat sky
c


      call skyb (nx,ny,mazi,malisa,vmin,fatceu,ceu,sigceu,izmax,threso)

c

      write (*,*) 'Flattening sky background. Last step done.'
      write (*,*) 'sky background mode, freq, sky sigma = ',ceu,izmax,
     ?sigceu


c
c     If sky background computation fails ...
c


      if (izmax.eq.0) then
      write (*,*) 'Sky background statistics failed.'
      go to 60
      endif

c

      write (*,*) 'sky background threshold for objects ID = ',threso
      write (*,*)
c

      if (threso.ge.mazi-2.d0) then
      write (*,*) 'Sky background statistics failed.'
      go to 60
      endif


      endif



c
c     Marca pixels a serem excluidos no processo de identificacao de objetos
c


      do i=2,ny-1
      do j=2,nx-1
      if (pixmat(j,i).ge.mazi .or. pixmat(j,i).le.threso) then
      imagem(j,i)=-1
      else
      imagem(j,i)=0
      endif
      enddo
      enddo

      i=1
      do j=1,nx
      imagem(j,i)=-1
      enddo

      i=ny
      do j=1,nx
      imagem(j,i)=-1
      enddo

      j=1
      do i=2,ny-1
      imagem(j,i)=-1
      enddo

      j=nx
      do i=2,ny-1
      imagem(j,i)=-1
      enddo


c
c     Identifica objetos com pixels acima de "thres"
c

      do i=1,idiobs
      iflag(i)=0
      enddo


      iob=0
      nstar=0
      alfa=zero
      excent=zero

      do 50 i=2,ny-1
      do 45 j=2,nx-1

      if (imagem(j,i).lt.0) go to 45

      ipego=0
      valor=pixmat(j,i)
      kx=j
      ky=i

 20   jm=kx
      im=ky

c
      do ii=im-1,im+1
      do 21 jj=jm-1,jm+1
      if (imagem(jj,ii).lt.0) go to 21
      if (pixmat(jj,ii).gt.valor) then
      valor=pixmat(jj,ii)
      kx=jj
      ky=ii
      ipego=1
      endif

 21   continue
      enddo

c

      if (ipego.eq.0) go to 45
      ipego=1

c

      if (kx-jm) 20,22,20
 22   if (ky-im) 20,24,20


c
c     Marca perimetro para candidato
c

 24   ifoi=0


c


 25   lado=0
      ibx=kx
      iby=ky



 30   continue

      lado=lado+1

      ix1=ibx-lado
      ix2=ibx+lado
      iy1=iby-lado
      iy2=iby+lado

      if (ix1.lt.1) ix1=1
      if (iy1.lt.1) iy1=1
      if (ix2.gt.nx) ix2=nx
      if (iy2.gt.ny) iy2=ny



c
c     Calcula histograma de fundo de ceu do perimetro em torno do
c     pixel central brilhante, para delimitar o perimetro do
c     objeto candidato
c

      n=0

      nlimo=2*(nx+ny)

      ii=iy1
      do 31 jj=ix1,ix2
      if (pixmat(jj,ii).ge.mazi) go to 31
      n=n+1
      if (n.gt.nlimo) go to 1500
      contag(n)=pixmat(jj,ii)
 31   continue

      ii=iy2
      do 32 jj=ix1,ix2
      if (pixmat(jj,ii).ge.mazi) go to 32
      n=n+1
      if (n.gt.nlimo) go to 1500
      contag(n)=pixmat(jj,ii)
 32   continue

      jj=ix1
      do 33 ii=iy1+1,iy2-1
      if (pixmat(jj,ii).ge.mazi) go to 33
      n=n+1
      if (n.gt.nlimo) go to 1500
      contag(n)=pixmat(jj,ii)
 33   continue

      jj=ix2
      do 34 ii=iy1+1,iy2-1
      if (pixmat(jj,ii).ge.mazi) go to 34
      n=n+1
      if (n.gt.nlimo) go to 1500
      contag(n)=pixmat(jj,ii)
 34   continue

c
 1500 continue

c


      zmin=1.d14
      zmax=-1.d14


      do nn=1,n
      if (contag(nn).gt.zmax) zmax=contag(nn)
      if (contag(nn).lt.zmin) zmin=contag(nn)
      enddo

      razao=zmax-zmin

c

      if (razao.lt.0.1d0) then
      fc=zmax
      go to 1512
      endif

c

      do nn=1,nhist
      ico(nn)=0
      histo(nn)=0.d0
      enddo

      do nn=1,n
      ko=1+(nhist-1)*(contag(nn)-zmin)/razao
      ico(ko)=ico(ko)+1
      histo(ko)=histo(ko)+contag(nn)
      enddo

      do 1501 nn=1,nhist
      if (ico(nn).eq.0) go to 1501
      histo(nn)=histo(nn)/ico(nn)
 1501 continue

c
c     Refina bins do histograma
c

      ncort=0.7d0*n

      mcort=0
      do nn=1,nhist
      mcort=mcort+ico(nn)
      if (mcort.gt.ncort) go to 1505
      enddo

 1505 cort=histo(nn)


      razao=cort-zmin

c
      if (razao.lt.0.1d0) then
      fc=cort
      go to 1512
      endif

c

      do nn=1,nhist
      ico(nn)=0
      histo(nn)=0.d0
      enddo

      do 1507 nn=1,n
      if (contag(nn).gt.cort) go to 1507
      ko=1+(nhist-1)*(contag(nn)-zmin)/razao
      ico(ko)=ico(ko)+1
      histo(ko)=histo(ko)+contag(nn)

 1507 continue


      do 1510 nn=1,nhist
      if (ico(nn).eq.0) go to 1510
      histo(nn)=histo(nn)/ico(nn)

 1510 continue

c
c     Pega o valor de frequencia mais alta e determina o valor
c     da contagem
c

      ixmaxi=-1

      do nn=1,nhist
      if (ico(nn).gt.ixmaxi) then
      ixmaxi=ico(nn)
      fc=histo(nn)
      endif
      enddo

c

 1512 continue

c

c     write (*,*) 'fc = ',fc

      if (fc.gt.threso) go to 30

c
c     Checa se o pixel central era mesmo o mais brilhante do perimetro,
c     e se nao for, recentra o perimetro
c

      valor=pixmat(ibx,iby)
      do 38 ii=iy1,iy2
      do 37 jj=ix1,ix2
      if (imagem(jj,ii).lt.0) go to 37
c     if (pixmat(jj,ii).ge.mazi) go to 37
      if (pixmat(jj,ii).gt.valor) then
      kx=jj
      ky=ii
      valor=pixmat(jj,ii)
      endif
     
 37   continue
 38   continue

      if (ifoi.ne.0) go to 40

      ifoi=1

      if (kx.ne.ibx) go to 25
      if (ky.ne.iby) go to 25

 40   iob=iob+1

c
c     Refina a determinacao do centro (x,y) do objeto candidato,
c     calculando seu baricentro
c

      xc=0.d0
      yc=0.d0
      cont=0.d0
      npixel=0
      do ii=iy1,iy2
      do 41 jj=ix1,ix2
      if (pixmat(jj,ii).ge.mazi) go to 41
c     if (imagem(jj,ii).lt.0) go to 41
      imagem(jj,ii)=-iob-1
      xc=xc+jj*pixmat(jj,ii)
      yc=yc+ii*pixmat(jj,ii)
      cont=cont+pixmat(jj,ii)
      npixel=npixel+1
 41   continue
      enddo

      xc=xc/cont
      yc=yc/cont



c
c     No. minimo de pixels para ajuste gaussiano NG=10 (o candidato
c     deve ser mais que uma imagem com apenas 1 pixel central e os
c     8 circunvizinhos ja' fundo de ceu)
c


      if (npixel.lt.ng) go to 45

      xcj=xc
      ycj=yc
      zlado=2.d0*lado


c
c     Ajuste gaussiano: Gaussiana circular
c

      r=lado

      mazz=mazi

c     write (*,*) 'xc yc = ',xc,yc


      call gcc (ix1,ix2,iy1,iy2,r,mazi,xc,yc,ex,ey,it,sig,hh,fc)

c     write (*,*)
c     write (*,*) 'xc yc = ',xc,yc


      mazi=mazz


      lado=r

      if (xc.lt.ix1-0.5d0) go to 45
      if (xc.gt.ix2+0.5d0) go to 45
      if (yc.lt.iy1-0.5d0) go to 45
      if (yc.gt.iy2+0.5d0) go to 45

      ix1=xc-r
      ix2=xc+r
      iy1=yc-r
      iy2=yc+r


c
c     Guarda apenas candidatos com FWHM entre valores pre-fixados no
c     arquivo de entrada (ex: entre 1" e 4"), ie, com valores
c     de seeing astronomicamente possiveis

      fwhm=2.d0*sig*1.177410023d0

      if (fwhm.lt.fmin) go to 45
      if (fwhm.gt.fmax) go to 45

      seeing=fwhm*scala


c     write (20,55) xcj,ycj,zlado,zlado
c55   FORMAT (1X,'image;BOX(',3(F8.2,','),F8.2,')')




c
c     Elimina boxes multiplas de um mesmo objeto
c     A box de maior tamanho e' preservada
c



      do 43  nn=1,nstar
 
      if (lado.le.ilado(nn)) then
 
      ixx1=xob(nn)-ilado(nn)
      ixx2=xob(nn)+ilado(nn)
      iyy1=yob(nn)-ilado(nn)
      iyy2=yob(nn)+ilado(nn)
 
      if (xc.lt.ixx1) go to 43
      if (yc.lt.iyy1) go to 43
      if (xc.gt.ixx2) go to 43
      if (yc.gt.iyy2) go to 43
 
      go to 45
 
      else
 
      if (xob(nn).lt.ix1) go to 43
      if (yob(nn).lt.iy1) go to 43
      if (xob(nn).gt.ix2) go to 43
      if (yob(nn).gt.iy2) go to 43
 
      iflag(nn)=1
 
 
      endif
 
 
 43   continue


      nstar=nstar+1

      xob(nstar)=xc
      yob(nstar)=yc
      ilado(nstar)=lado
      seng(nstar)=seeing
c     altu(nstar)=hh-fc
      altu(nstar)=hh   

      exgcc(nstar)=ex
      eygcc(nstar)=ey
      sgcc(nstar)=sig
      fgcc(nstar)=fc

c

 45   continue
 50   continue


c


      open (20,file=ids9)


      nest=0

      do 57 k=1,nstar

      if (iflag(k).gt.0) go to 57

      nest=nest+1

      xob(nest)=xob(k)
      yob(nest)=yob(k)
      ilado(nest)=ilado(k)
      seng(nest)=seng(k)

      altu(nest)=altu(k)
      ialtu(nest)=nest

      exgcc(nest)=exgcc(k)
      eygcc(nest)=eygcc(k)
      sgcc(nest)=sgcc(k)
      fgcc(nest)=fgcc(k)


c

      j=xob(nest)
      i=yob(nest)

      ior(nest)=nest

      if (imagem(j,i).eq.-1) then
      nval(nest)=mazi+1
      else
      nval(nest)=pixmat(j,i)
      endif



      xlado=2.d0*ilado(nest)
      ylado=2.d0*ilado(nest)

      write (20,55) xob(nest),yob(nest),xlado,xlado
 55   FORMAT (1X,'image;BOX(',3(F8.2,','),F8.2,')')


 57   continue

      close (20)


c
c     Minimum number of objects not reached ?
c

c     write (*,*) 'nnn = ',nest

c     stop


      if (nest.lt.4) go to 60


c
c     Determina FWHM tipico, via histograma de fwhm dos objetos
c     Sao usados objetos dentro da faixa:
c
c     30% mais fracos < mag < 20% mais brilhantes
c
c     (saturados excluidos)
c


c
      call ordem (idiobs,nest,ior,nval)

c


      do i=1,nest
      if (nval(i).ge.mazi) go to 100
      enddo
 100  n2=i-1


      n1=0.3d0*n2
      n2=0.8d0*n2

c
c     Montando histograma de fwhm
c


      n=0

      do 105 i=n1,n2
      n=n+1
      contag(n)=seng(ior(i))
 105  continue


      zmin=1.d14
      zmax=-1.d14


      do nn=1,n
      if (contag(nn).gt.zmax) zmax=contag(nn)
      if (contag(nn).lt.zmin) zmin=contag(nn)
      enddo

      razao=zmax-zmin

      do nn=1,nhisfw
      ico(nn)=0
      histo(nn)=0.d0
      enddo

      do nn=1,n
      ko=1+(nhisfw-1)*(contag(nn)-zmin)/razao
      ico(ko)=ico(ko)+1
      histo(ko)=histo(ko)+contag(nn)
      enddo

      do 110 nn=1,nhisfw
      if (ico(nn).eq.0) go to 110
      histo(nn)=histo(nn)/ico(nn)
 110  continue




c
c     Refina bins do histograma de fwhm
c


      xmax=-1.d14


      do nn=1,nhisfw
      if (ico(nn).gt.xmax) then
      xmax=ico(nn)
      cont=histo(nn)
      kk=nn
      endif
      enddo

c

      do nn=1,n
      ior(nn)=nn
      nval(nn)=10000*dabs(contag(nn)-histo(kk))
      enddo
c
      call ordem (idiobs,n,ior,nval)
c

      ncort=0.7d0*n

      zmin=1.d14
      zmax=-1.d14

      do nn=1,ncort

      if (contag(ior(nn)).lt.zmin) zmin=contag(ior(nn))
      if (contag(ior(nn)).gt.zmax) zmax=contag(ior(nn))

      enddo

c


      razao=zmax-zmin

      do nn=1,nhisfw
      ico(nn)=0
      histo(nn)=0.d0
      enddo

      do nnn=1,ncort
      nn=ior(nnn)
      ko=1+(nhisfw-1)*(contag(nn)-zmin)/razao
      ico(ko)=ico(ko)+1
      histo(ko)=histo(ko)+contag(nn)

      enddo

      fra=razao/(nhisfw-1.d0)

      do 113 nn=1,nhisfw
      if (ico(nn).eq.0) then
      histo(nn)=(nn-1.d0)*fra+zmin
      go to 113 
      endif
      histo(nn)=histo(nn)/ico(nn)

 113  continue


c
c     Pega valor tipico de fwhm e sigma
c

 

      xmax=-1.d14

      do n=1,nhisfw
      if (ico(n).gt.xmax) then
      nn=n
      xmax=ico(n)
      endif
      enddo

      fwhm=histo(nn)

      cort=xmax/2.d0

      sigfwh=razao/2.d0


c     k=0
      do n=1,nhisfw
c     k=k+1

      if (nn+n.le.nhisfw .and. ico(nn+n).ne.0) then

      if (ico(nn+n).lt.cort) then
      sigfwh=2.d0*(histo(nn+n)-fwhm)
      go to 115
      endif

      endif


      if (nn-n.ge.1 .and. ico(nn-n).ne.0) then

      if (ico(nn-n).lt.cort) then
      sigfwh=2.d0*(fwhm-histo(nn-n))
      go to 115
      endif

      endif

      enddo

 115  continue


c
c     Correcao de Field Distortion Pattern. 
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
c     Separa as N (ex.:N=30) estrelas mais brilhantes do campo atual
c     para calcular o offset provisorio.
c

      do i=1,nest
      ialtu(i)=i
      enddo

      do i=1,nest
      nval(i)=altu(i)**2
      enddo
      call ordem (idiobs,nest,ialtu,nval)



c
c     Determina regiao de alfa e delta do ceu para extracao de
c     estrelas de catalogo
c
c     rac = centro alfa  em graus
c     dec = centro delta em graus
c     
c     drac = lado da box alfa  em graus
c     ddec = lado da box delta em graus
c

      open (77,file=centro)

 400  continue

      imfits=''

      read(77,402) iah,iam,sa,isig,idg,idm,ds,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,djm,dj,iexps,ichfil,imfits,mchobj

 402  format(1x,i2,1x,i2,1x,f7.4,1x,a1,i2,1x,i2,1x,f6.3,2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,f16.8,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20)



      if (imfits.ne.infits) go to 400

      close (77)



c

      epoj=2000D0+(dj-2451545D0)/365.25D0

      rac=hmsgms(iah,iam,sa)*15.d0
      dec=hmsgms(idg,idm,ds)
      if (isig.eq.menos) dec=-dec


c


      drac=nx*scala2/3600.d0
      ddec=ny*scala2/3600.d0

c
c     Calcula area da reducao tp+pm
c

      aremax=2.d0*drac
      aremay=2.d0*ddec

      if (akey.lt.0.98d0) then

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
c     Checa se polo (norte ou sul) cai no campo
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
c     Polo (sul ou norte) cai  no campo    -> krcat = 2: extrair catalogos com projecao no
c     plano tangente, estrela a estrela (processo mais lento)
c
c     Polo (sul ou norte) nao cai no campo -> krcat = 1: extrair catalogos com (RA,Dec) min_max,
c     usando indexacao dos catalogos (procedimento acelerado) 
c


      krcat=2

      if (dabs(xx).gt.areax) krcat=1
      if (dabs(yy).gt.areay) krcat=1



c
c     Pega as estrelas 2MASS da regiao em torno do CCD     
c

      write (*,*) 'Files searched for catalogue extraction:'

c
c
c     - ra2ma,de2ma em graus (alfas e deltas do 2MASS)
c     - era2ma,ede2ma em segundos de arco (erro em alfa e delta)
c     - mags em mags mesmo: J, H e K.
c
c

      if (krcat.eq.2) then

      call tmass (mraiz,rac,dec,drac,ddec,ra2ma,de2ma,era2ma,ede2ma,
     ?dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,n2mass)

      else

      call stmass (mraiz,rac,dec,drac,ddec,ra2ma,de2ma,era2ma,ede2ma,
     ?dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,n2mass)

      endif



c
c     Inicializando rotinas de estimativa de tempo de execucao da
c     reducao das imagens com o PRAIA
c

      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)

      tt0=tempop/3600.d0




c
c     Varredura espectral da escala de pixel
c


      scl=0.d0
      scl2=0.d0
      nscl=0
      dscl=0.d0

      noest=-1 

      ncomum=0

      scala=scala2+ecala


c
c     Comeca o loop de varredura espectral da escala de pixel
c

      write (*,*) 
      write (*,*) 

c

 1001 continue


      scala=scala-ecala*2.d0

      if (scala.lt.scala1) go to 1002



      areax=grarad*nx*scala/3600.d0
      areay=grarad*ny*scala/3600.d0



c     write (*,*) 'areax areay = ',areax*radgra*60.d0,areay*radgra*60.d0


      do i=1,icofsp
      coefx2(i)=0.d0
      coefy2(i)=0.d0
      enddo



c
c     Reconhece estrelas do 2MASS dentre as identificadas no campo
c

      rcala=scala

      call idxy2m (areax,areay,idiobs,nval,ior,rcala,erpix,rac,dec,
     ?nbcat,nbmed,nest,ialtu,xob,yob,n2mass,id2ma,ra2ma,de2ma,dmgj,xold,
     ?yold,xra2ma,yde2ma,ireflex,ecala,tt,coefx2,coefy2,ncomum)



c
c     Estatistica da escala de pixel, apartir do ajuste de reconhecimento das estrelas
c     do catalogo de referencia (2MASS)
c

  
      if (ierro.eq.1) then
      ierro=0
      go to 1001
      endif

      if (ncomum.gt.noest) then
      vcala=scala
      creax=areax
      creay=areay
      noest=ncomum
      kreflex=ireflex
      do mmmm=1,icofsp
      coefxr(mmmm)=coefx2(mmmm)
      coefyr(mmmm)=coefy2(mmmm)
      enddo
      endif

      dnest=ncomum
      dnest=dsqrt(dnest)
      dscl=dscl+dnest
      scl=scl+scala*dnest
      scl2=scl2+dnest*scala**2
      nscl=nscl+1

      write (*,*) 'No. identified 2MASS stars, tested pixel scale ("/pix
     ?) = ',ncomum,rcala


      go to 1001

c

 1002 continue


c
c     Escala de pixel provisoria (media ponderada pela raiz quadrada do no. de estrelas
c     identificadas. Estimativa do erro eh definitiva.
c


c     rcala=scl/dscl

c     ercala=(1.d0/(nscl-1.d0-1.d0))*scl2/((1.d0/nscl)*dscl)

      rcala=vcala
 
      ercala=ecala



      write (*,*) 
      write (*,*) 'Provisional pixel scale & error ("/pix) = ',rcala,
     ?ercala



c
c     Refina a determinacao da escala de pixel.
c

      coefxr(1)=0.d0
      coefyr(1)=0.d0


c
c     Ultima identificacao (rapida) de estrelas 2MASS
c

      dcala=grarad*rcala/3600.d0
      decala=grarad*ercala/3600.d0


      call mdxy2m (creax,creay,ngrau,coefxr,coefyr,kreflex,idiobs,nval,
     ?ior,dcala,erpix,rac,dec,nbcat,nbmed,nest,ialtu,xob,yob,n2mass,
     ?id2ma,ra2ma,de2ma,dmgj,xold,yold,xra2ma,yde2ma,ireflex,decala,tt)



      do mmmm=1,icofsp
      coefx2(mmmm)=0.d0
      coefy2(mmmm)=0.d0
      enddo

      do i=1,idiobs
      itira2(i)=0
      enddo

c
c     Reducao (RA,DEC) com 2MASS.
c

      call posred (ireflex,rac,dec,id2ma,n2mass,ra2ma,de2ma,nest,xob,
     ?yob,pcort2,ngrau,ngrau3,ngrau5,nstart,nfinal,xra2ma,yde2ma,
     ?era2,ede2,alfsi2,delsi2,alfre2,delre2,coefx2,coefy2,ecofx2,ecofy2,
     ?itira2,avam,dvam)


c
c     Valor final para escala de pixel
c


      rscala=dsqrt(coefx2(2)**2+coefx2(3)**2)+dsqrt(coefy2(2)**2+
     ?coefy2(3)**2)
 
      rscala=rscala/2.d0

      rscala=rscala*radgra*3600.d0

      dnx=nx/2.d0
      dny=ny/2.d0

      ercala=3.d0*dsqrt(avam**2+dvam**2)/dsqrt(dnx**2+dny**2)


      ierro=0



c
c     Escreve resultado
c

      write (*,*)
      write (*,*) 'Final result:'
      write (*,*)

      write (*,454) rscala,ercala,infits
 454  format('Pixel scale & error ("/pixel), image = ',2(1x,f10.7),1x,
     ?a50)

      write (*,*)
      write (*,*)

      write (97,455) rscala,ercala,infits
 455  format(2(1x,f10.7),1x,a50)



 60   continue

c

      close (97)


      write (*,*)
      write (*,61) 
 61   format (23x,'Processing terminated.')
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '

      end



c
c     subrotina refits
c
c
c     Esta subrotina le imagens fits integer ou floating point
c     em f77, sem uso de programas externos (como qfits, etc).
c 
c
C      Last update:   27/08/2009 - Marcelo Assafin
C
C
      subroutine refits (infits,nx,ny,nheads,ichobj,ipflag,bscale,bzero,
     ?kswap,iswap,bitpix)

      IMPLICIT REAL *8 (A-H,O-Z)


      character*150 infits
      character*69 ichobj
      character*9  naxisx,naxisy,bitpx,ibscal,ibzero
      character*1  ler(2880),ibr
      character*9 iend,kend
      character*50 erro
      character*20 imaux
      character*4 jmaux
      character*9 sista
      character*29 systa

      dimension iwork2(1440),swork(2880),iby4(4)
      dimension work4(720),iwork4(720)
      dimension work8(360),iwork8(360)

      dimension pixmat(5001,5001)
      integer*2 bitpix

      real*4 pixmat
      integer*2 iwork2
      integer*4 iwork4
      real*4 work4
      integer*8 iwork8
      real*8 work8

      integer*1 swork,iby4

      common /a4/pixmat

      data naxisx /'NAXIS1  ='/
      data naxisy /'NAXIS2  ='/
      data bitpx  /'BITPIX  ='/
      data ibscal /'BSCALE  ='/
      data ibzero /'BZERO   ='/
      data iend   /'END      '/
      data ibr    /' '/

c

      nbytes=2880

      if=1

c
c     Abre arquivo fits
c

      

      open(if,file=infits,access='direct',form='unformatted',recl=2880)


c
c     Abre arquivo auxiliar
c

      sista='rm -f -r '

      imaux=''

      imaux(1:16)='PRAIA_refits.aux'


      do 1 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux(17:20)=jmaux(1:4)

      open(99,file=imaux,access='direct',form='unformatted',
     ?recl=2880,status='old',err=2)
      close (99)
 1    continue

 2    close (99)

      open(99,file=imaux,access='direct',form='unformatted',
     ?recl=2880)

      systa=sista//imaux

c
c     Determina quantos headers existem na imagem
c

      kend=''

      m=0
 10   m=m+1
      read(1,rec=m) ler

      do j=0,35

      do i=1,9
      kend(i:i)=ler(j*80+i)
      enddo

      if (kend.eq.iend) go to 20

      enddo

      go to 10

c

 20   nheads=m

c     write (*,*)
c     write (*,*) 'nheads ',nheads
c     stop

c
c     Determina dimensao nx da matriz
c

      do i=1,nheads

      read(1,rec=i) ler

      call find (ler,naxisx,key,dnx)

      if (key.gt.0) nx=dnx

      enddo



c
c     Determina dimensao ny da matriz
c

      do i=1,nheads

      read(1,rec=i) ler

      call find (ler,naxisy,key,dny)

      if (key.gt.0) ny=dny

      enddo

c
c     Determina bitpix (imagem integer ou real?)
c

      if (bitpix.eq.-99) then

      bitpix=16

      do i=1,nheads

      read(1,rec=i) ler

      call find (ler,bitpx,key,bpix)

      if (key.gt.0) bitpix=bpix

      enddo

      endif


      i=bitpix
      if (i.lt.0) i=-i
      
      ibytes=i/8+0.1
      kwork=nbytes/ibytes+0.1



c
c     Determina bscale 
c

      if (ipflag.ne.1) then

      bscale=1.d0

      do i=1,nheads

      read(1,rec=i) ler

      call find (ler,ibscal,key,bsc)

      if (key.gt.0) bscale=bsc

      enddo

      endif


c
c     Determina bzero
c


      if (ipflag.ne.1) then

      bzero=0.d0

      do i=1,nheads

      read(1,rec=i) ler

      call find (ler,ibzero,key,bzr)

      if (key.gt.0) bzero=bzr

      enddo

      endif


c
c     Le a matriz de pixels
c


c
c     Checa byte-swap (litteendian ou bigendian)
c
c     kswap : chave do usuario
c
c          kswap = 0 determinacao automatica de byte-swap
c          kswap = 1 sem byte-swap (definido pelo usuario)
c          kswap = 2 com byte-swap (definido pelo usuario)
c
c     Determinacao automatica (kswap=0):
c
c     iswap=1  sem byte-swap
c     iswap=2  com byte-swap
c
c


      if (kswap.ne.0) then

      iswap=kswap

      go to 50

      endif


      irec=nheads+((nx*ny)/2.d0)*ibytes/nbytes 

c

      if (bitpix.gt.0) then

      if (ibytes.eq.2) read (1,rec=irec) iwork2
      if (ibytes.eq.4) read (1,rec=irec) iwork4
      if (ibytes.eq.8) read (1,rec=irec) iwork8

      else

      if (ibytes.eq.4) read (1,rec=irec) work4
      if (ibytes.eq.8) read (1,rec=irec) work8

      endif     

c
c     Media do valor absoluto sem swap
c

      c1=0.d0


      if (bitpix.gt.0) then

      do i=1,kwork
      if (ibytes.eq.2) c=bscale*iwork2(i)+bzero
      if (ibytes.eq.4) c=bscale*iwork4(i)+bzero
      if (ibytes.eq.8) c=bscale*iwork8(i)+bzero
      c1=c1+dabs(c)
      enddo

      c1=c1/kwork

      else

      do i=1,kwork
      if (ibytes.eq.4) c=bscale*work4(i)+bzero
      if (ibytes.eq.8) c=bscale*work8(i)+bzero
      c1=c1+dabs(c)
      enddo

      c1=c1/kwork

      endif


c
c     Testa com swap
c


      call swap (if,bitpix,ibytes,nbytes,irec,iwork2,iwork4,iwork8,
     ?work4,work8,swork)

c
c     Media do valor absoluto com swap
c

      c2=0.d0


      if (bitpix.gt.0) then

      do i=1,kwork
      if (ibytes.eq.2) c=bscale*iwork2(i)+bzero
      if (ibytes.eq.4) c=bscale*iwork4(i)+bzero
      if (ibytes.eq.8) c=bscale*iwork8(i)+bzero
      c2=c2+dabs(c)
      enddo

      c2=c2/kwork

      else

      do i=1,kwork
      if (ibytes.eq.4) c=bscale*work4(i)+bzero
      if (ibytes.eq.8) c=bscale*work8(i)+bzero
      c2=c2+dabs(c)
      enddo

      c2=c2/kwork

      endif

c
c     Define swap
c

      erro=''

      write (erro,*) c1

      do i=1,50
      if (ichar(erro(i:i)).ge.48 .and. ichar(erro(i:i)).le.57) go to 35
      enddo

      c1=1.d14

 35   continue

c

      erro=''
     
      write (erro,*) c2

      do i=1,50
      if (ichar(erro(i:i)).ge.48 .and. ichar(erro(i:i)).le.57) go to 40
      enddo

      c2=1.d14

 40   continue    

c

      if (c2.lt.c1) then
      iswap=2
      else
      iswap=1
      endif     


c     write (*,*) 'c1 c2 ',c1, c2

c
c     Lendo matriz
c

 50   continue


c

      block=nx*ny*ibytes/nbytes 
      nblock=block
      iresto=(block-nblock)*nbytes
      iresto=iresto/ibytes


      j=0
      i=1

      do m=1,nblock

      irec=m+nheads

c


      if (iswap.eq.1) then

      if (bitpix.gt.0) then

      if (ibytes.eq.2) read (1,rec=irec) iwork2
      if (ibytes.eq.4) read (1,rec=irec) iwork4
      if (ibytes.eq.8) read (1,rec=irec) iwork8

      else

      if (ibytes.eq.4) read (1,rec=irec) work4
      if (ibytes.eq.8) read (1,rec=irec) work8

      endif

      else

      call swap (if,bitpix,ibytes,nbytes,irec,iwork2,iwork4,iwork8,
     ?work4,work8,swork)

      endif

c

      if (bitpix.gt.0) then

      do mm=1,kwork

      j=j+1
      if (j.gt.nx) then
      j=1
      i=i+1
      endif

      if (ibytes.eq.2) pixmat(j,i)=bscale*iwork2(mm)+bzero
      if (ibytes.eq.4) pixmat(j,i)=bscale*iwork4(mm)+bzero
      if (ibytes.eq.8) pixmat(j,i)=bscale*iwork8(mm)+bzero

      enddo

      else

      do mm=1,kwork

      j=j+1
      if (j.gt.nx) then
      j=1
      i=i+1
      endif

      if (ibytes.eq.4) pixmat(j,i)=bscale*work4(mm)+bzero
      if (ibytes.eq.8) pixmat(j,i)=bscale*work8(mm)+bzero

      enddo      
 

      endif

      enddo

c
c     Ultimo pedaco de bloco da matriz (se existente)
c

      if (iswap.eq.1) then

      if (bitpix.gt.0) then

      if (ibytes.eq.2) read (1,rec=irec+1) (iwork2(m),m=1,iresto)
      if (ibytes.eq.4) read (1,rec=irec+1) (iwork4(m),m=1,iresto)
      if (ibytes.eq.8) read (1,rec=irec+1) (iwork8(m),m=1,iresto)

      else

      if (ibytes.eq.4) read (1,rec=irec+1) (work4(m),m=1,iresto)
      if (ibytes.eq.8) read (1,rec=irec+1) (work8(m),m=1,iresto)
     
      endif

      else

      nbytes=iresto*ibytes

      call swap (if,bitpix,ibytes,nbytes,irec,iwork2,iwork4,iwork8,
     ?work4,work8,swork)

      endif

c
      if (bitpix.gt.0) then

      do mm=1,iresto

      j=j+1
      if (j.gt.nx) then
      j=1
      i=i+1
      endif

      if (ibytes.eq.2) pixmat(j,i)=bscale*iwork2(mm)+bzero
      if (ibytes.eq.4) pixmat(j,i)=bscale*iwork4(mm)+bzero
      if (ibytes.eq.8) pixmat(j,i)=bscale*iwork8(mm)+bzero

      enddo      

      else

      do mm=1,iresto

      j=j+1
      if (j.gt.nx) then
      j=1
      i=i+1
      endif

      if (ibytes.eq.4) pixmat(j,i)=bscale*work4(mm)+bzero
      if (ibytes.eq.8) pixmat(j,i)=bscale*work8(mm)+bzero

      enddo      

      endif
c

      close (1)
      close (99)

c
c     Debug
c

c     read (*,*) j,i
c     j=100
c     i=101
c     write (*,*) 'ix iy = ',j,i
c     write (*,*) 'pixmat = ',pixmat(j,i)
c     write (*,*) 'nx ny ',nx,ny
c     write (*,*) 'swap ',iswap
c     write (*,*) 'c1 c2 = ',c1,c2
c     write (*,*) 'nheads = ',nheads
c     write (*,*) 'bitpix = ',bitpix
c     write (*,*) 'bscale = ',bscale
c     write (*,*) 'bzero  = ',bzero 
c     stop



      call system (systa)

c

      return
      end




c
c
c     subroutine find
c
c
c     Acha o valor numerico correspondente a palavra chave
c     do cabecalho fits
c
c     ler     = extracao do header 
c     word    = a palavra do header
c     key     = +1 achou
c             = -1 nao achou
c     valor   = valor numerico encontrado
c
c


      subroutine find (ler,word,key,valor)

      IMPLICIT REAL *8 (A-H,O-Z)
      

      character*1  ler(2880),iplic,ibr,ibar
      character*9  word,kend
      character*1 ivalor(71)



      data iplic /"'"/
      data ibar  /"/"/
      data ibr   /' '/

c

      key=-1
      valor=-1.d14

      icol=1
      id=71

c

      do j=0,35

      do i=1,9
      kend(i:i)=ler(j*80+i)
      enddo

      if (kend.eq.word) then
      key=+1
      go to 20
      endif

      enddo

      go to 50

c

 20   j=j*80
      do i=10,80
      if (ler(j+i).ne.ibr .and. ler(j+i).ne.iplic) go to 30
      enddo
 30   i1=i
      do i=i1+1,80
      if (ler(j+i).eq.ibr .or. ler(j+i).eq.ibar) go to 40
      enddo

 40   i2=i-1

c
      do i=1,id
      ivalor(i)=ibr
      enddo
c

      n=0

      do i=j+i1,j+i2
      n=n+1
      ivalor(n)=ler(i)
      enddo

c

      call chanum (icol,id,ivalor,valor)


 50   continue
      return

      end




c
c     subrotina chanum
c
c     Pega uma string (character) de numeros e extrai o numero da
c     coluna, sem abrir arquivos temporarios para isso.
c
c     string  = contem a string completa
c     palavra = variavel de trabalho contendo a string completa
c     b2 = contem o numero extraido correspondente a coluna dada
c     
c
c     Ultima modificacao: M. Assafin  27/Agosto/2009
c
c
c
      subroutine chanum (icol,id,string,valor)

      implicit real *8 (a-h,o-z)

      integer*8 n

      dimension ni(id+2),nf(id+2)

      character*1 string(id),palavra(id+2),ibra

c

      ibra=' '

c

      do i=1,id+2
      palavra(i)=ibra
      enddo

      do i=1,id
      palavra(i+1)=string(i)
      enddo

c
c     Checando colunas pelos espacos em branco
c

      do i=1,id
      ni(id)=0
      nf(id)=0
      enddo

      ki=0
      kf=0

c
c     Onde estah o comeco do numero
c

      do i=2,id+2

      if (palavra(i-1).eq.ibra .and. palavra(i).ne.ibra) then
      ki=ki+1
      ni(ki)=i
      endif

      enddo

c
c     Onde estah o fim do numero
c

      do i=2,id+2

      if (palavra(i-1).ne.ibra .and. palavra(i).eq.ibra) then
      kf=kf+1
      nf(kf)=i-1
      endif
      
      enddo

c
c     Checa se numero eh positivo ou negativo
c

      isig=+1

      i=ni(icol)
      if (palavra(i).eq.'-') isig=-1
  

c
c     Checa se numero estah escrito em notacao "E" ou "D"
c

      iep=0
      ie=0
      isige=+1


      do i=ni(icol),nf(icol)
      if (palavra(i).eq.'e'.or.palavra(i).eq.'E'.or.palavra(i).eq.'d'
     ?.or.palavra(i).eq.'D') ie=i
      enddo

      if (ie.ne.0) then

      iee=ie

      if (palavra(ie+1).eq.'-') then
      isige=-1
      endif

      if (palavra(ie+1).eq.'-'.or.palavra(ie+1).eq.'+') then
      ie=ie+2
      else
      ie=ie+1
      endif


      iep=0
      j=0
      k=nf(icol)-ie+1
      do i=ie,nf(icol)
      icomp=ichar(palavra(i))
      j=j+1
      iep=iep+(icomp-48)*10.d0**(k-j)
      enddo

      nf(icol)=iee-1

      endif

      expo=10.d0**(isige*iep)


c
c     Determina onde o ponto decimal estah, se nao for mumero inteiro
c

      m=0
      do i=nf(icol),ni(icol),-1
      if (palavra(i).eq.'.') m=i-nf(icol)
      enddo

c
c     Determina os algarismos do numero
c

      k=0
      do i=ni(icol),nf(icol)
      if (palavra(i).ne.'.' .and. palavra(i).ne.'+' .and. palavra(i).
     ?ne.'-') k=k+1
      enddo

      n=0
      j=0
      do i=ni(icol),nf(icol)
      icomp=ichar(palavra(i))
      if (icomp.ge.48 .and. icomp.le.57) then
      j=j+1
      n=n+(icomp-48)*10.d0**(k-j)
      endif
      enddo


      valor=expo*isig*n*10.d0**m

      return
      end




c
c
c     subroutine swap
c
c
c     Swap dos bytes da imagem fits.
c
c     Imagem pode ser integer ou floating point
c
c 
c     Ultima modificacao: M. Assafin 27/08/2009
c
c


      subroutine swap(if,bitpix,ibytes,nbytes,irec,iwork2,iwork4,iwork8,
     ?work4,work8,swork)


      IMPLICIT REAL *8 (A-H,O-Z)



      dimension iwork2(1440),swork(2880),iby4(4)
      dimension work4(720),iwork4(720)
      dimension work8(360),iwork8(360)

      integer*2 bitpix

      real*4 pixmat
      integer*2 iwork2
      integer*4 iwork4
      real*4 work4
      integer*8 iwork8
      real*8 work8

      integer*1 swork,iby4

c

      read (if,rec=irec) swork


      do k=ibytes,nbytes,ibytes

      do m=1,ibytes
      iby4(m)=swork(k-m+1)
      enddo

      do m=1,ibytes
      swork(k-ibytes+m)=iby4(m)
      enddo

      enddo

      write (99,rec=1) swork

      if (bitpix.gt.0) then

      if (ibytes.eq.2) read (99,rec=1) iwork2
      if (ibytes.eq.4) read (99,rec=1) iwork4
      if (ibytes.eq.8) read (99,rec=1) iwork8

      else

      if (ibytes.eq.4) read (99,rec=1) work4
      if (ibytes.eq.8) read (99,rec=1) work8

      endif     

c



c
      return
      end





C     SUBROTINA GCC
C
C     PROPOSITO
C
C       ESTA SUBROTINA CALCULA O CENTRO DE UMA IMAGEM DIGITALIZADA  PELO  
C       AJUSTE DE GAUSSIANA CIRCULAR. O TRIMMING E' CIRCULAR. OS OBJETOS
C       SAO ALISADOS ANTES DO AJUSTE POR  UM  FILTRO  BINOMIAL  DE  TRES
C       CANAIS (BEVINGTON,1969) EXTENDIDO A DUAS  DIMENSOES. VALIDO PARA
C       IMAGENS CCD 400x600 ou 770x1152 do LNA.
C
C     USO
C
C       CALL GCC (IDATA,IOUT,NPIXEX,NPIXEY,IX0,IY0)
C
C     DESCRICAO DOS PARAMETROS
C
C       IOUT   - ARQUIVO DE MEDIDAS X,Y DE OBJETOS
C       IDATA  - ARQUIVO DE DADOS DO IDENTIFICADOR
C       NPIXEX - NUMERO DE PIXELS NA COORDENADA X
C       NPIXEY - NUMERO DE PIXELS NA COORDENADA Y
C       NX     - LIMITE INFERIOR X DO PROCESSO TRIMMING
C       MX     - LIMITE SUPERIOR X DO PROCESSO TRIMMING
C       NY     - LIMITE INFERIOR Y DO PROCESSO TRIMMING
C       MY     - LIMITE SUPERIOR Y DO PROCESSO TRIMMING
C       DISMAX - VETOR DA DISTRIBUICAO MARGINAL EM X
C       DISMAY - VETOR DA DISTRIBUICAO MARGINAL EM Y
C       CENTRX - CENTRO DA MARGINAL X OBTIDO PELO AJUSTE GAUSSIANO
C       CENTRY - CENTRO DA MARGINAL Y OBTIDO PELO AJUSTE GAUSSIANO
C       SIGX   - LARGURA A MEIA ALTURA DA GAUSSIANA DA MARGINAL X
C       SIGY   - LARGURA A MEIA ALTURA DA GAUSSIANA DA MARGINAL Y
C       SIGDEX - ERRO MEDIO QUADRATICO DO CENTRO X
C       SIGDEY - ERRO MEDIO QUADRATICO DO CENTRO Y
C       ALTURX - ALTURA DA GAUSSIANA AJUSTADA A MARGINAL X
C       FUNDOX - FUNDO DO CEU AJUSTADO COM A GAUSSIANA 
C       RESIDX - RESIDUO EM DENSIDADE DO AJUSTE
C       XSIGMA - RAIZ DOS ELEMENTOS DIAGONAIS DA MATRIZ DE COVARIANCIA
C       PARAMX - PARAMETROS DE ENTRADA E SAIDA DA GAUSSIANA E FUNDO X
C       PARAMY - PARAMETROS DE ENTRADA E SAIDA DA GAUSSIANA E FUNDO Y
C                1 - ALTURA
C                2 - CENTRO
C                3 - SIGMA DA GAUSSIANA
C                4 - FUNDO DE CEU CONSTANTE
C       PARAM  - PARAMETROS DE ENTRADA E SAIDA DA GAUSSIANA SIMETRICA
C                1 - CENTRO X
C                2 - CENTRO Y
C                3 - SIGMA DA GAUSSIANA
C                4 - ALTURA
C                5 - FUNDO DE CEU CONSTANTE
C       ABIS   - ABSSISSA DAS MARGINAIS X E Y PARA PLOTE COM MONGO
C       DIST   - DISTRIBUICOES MARGINAIS X E Y PARA PLOT COM MONGO
C       FITAR  - AJUSTE GAUSSIANO X E Y PARA PLOTE COM MONGO
C       XX0    - COORDENADA X DA JANELA
C       YY0    - COORDENADA Y DA JANELA
C       PASSOX - PASSO DE LEITURA EM X
C       PASSOY - PASSO DE LEITURA EM Y
C       IXSTEP - PASSO DE LEITURA EM X
C       IYSTEP - PASSO DE LEITURA EM Y
C       X0     - COORDENADA DA JANELA (X)
C       Y0     - COORDENAMDA DA JANELA (Y)
C
C     SUBROTINAS E SUBPROGRAMAS REQUERIDOS
C
C       GAUSIC (KEY,RAIO,NX,NY,MX,MY,NTERMS,PARAM,DELTAX,
C                XSIGMA,FLAMDA,RESIDX)
C
C
C     COMENTARIOS
C
C       A SUBROTINA GAUSIM E' UMA VERSAO MODIFICADA DA SUBROTINA  CURFIT
C       PRESENTE NO LIVRO " DATA  REDUCTION AND ERROR  ANALYISIS FOR THE
C       PHYSICAL SCIENCES" DE PHILIP R. BEVINGTON (1969), PAGS 237-239.
C
C       A DENSIDADE DE CADA PIXEL E' ASSOCIADA AO MEIO DESTE, SENDO  QUE
C       OS VALORES DE X0 E Y0, OU SEJA, A ORIGEM  ENCONTRA-SE  ASSOCIADA
C       AO CANTO INFERIOR ESQUERDO DO PRIMEIRO PIXEL (1,1), OS QUAIS TEM
C       COORDENADAS (0,0).
C
C       ESTE PROGRAMA ESTA' ESCRITO EM FORMA DE SUBROTINA PARA O PACOTE
C       DE TRATAMENTO DE IMAGENS "CCDMEDE".
C
C       A VARIAVEL DENSID SUBSTITUI A VARIAVEL COORDX ORIGINAL PARA
C       ECONOMIA DE ESPACO DE MEMORIA.
C
C       DIMENSOES EFETIVAS DOS CCDS (SEM A MOLDURA)
C            CCD  400x600  -  382x576
C            CCD  770x1152 -  742x1149
C
C
C       Orientacao dos CCDS
C
C                 400x600    x cresce - declinc.  cresce
C                 400x600    y cresce - asc. reta decresce
C                 770x1152   x cresce - declinc.  cresce
C                 770x1152   y cresce - asc. reta cresce
C
C       Ao final as coordenadas X e Y da matriz de pixels sao invertidas
C       bem como os valores associados a  X e Y,  e  entao  escritos  no
C       arquivo de medida.
C
C
C       Esta e' uma versao adaptada de gcc.f para o nosso caso.
C 
C                           M. Assafin 22/Dez/2004
C
C
C     SUBROUTINE GCC (IDATA,IOUT,NPIXEX,NPIXEY,IX0,IY0)

      subroutine  GCC (nx,mx,ny,my,raio,maximo,bx,by,sigdex,sigdey,
     ?icontt,sigx,alturx,fundox)

      IMPLICIT REAL*8 (A-H,O-Z)
      real*4 pixmat(5001,5001),maximo,xmaxi
      DIMENSION DELTAX(5),XSIGMA(5),PARAM(5)
 
      COMMON /A4/PIXMAT
      COMMON /A14/IERRO
C
C     INICIALIZACAO E ENTRADA DE DADOS
C

      xmaxi=maximo

      ramax=raio

      escl=1.d0
      X0 = 0.D0
      Y0 = 0.D0
      PASSOX = 1.D0
      PASSOY = 1.D0
      CENTRX = 0.D0
      CENTRY = 0.D0
C
C     CONVERGENCIA EM DENSIDADE = 1% ; EM POSICAO = 0,001 arcsec ou
c     1/20 avos do pixel
C
      DLIMIT=1.D-2
c     PLIMIT=1.D-3
      plimit=0.05d0
C
      ICONTT=0
      ICONTX=1
      TRANSX=0.D0
      TRANSY=0.D0
      SIGX=0.D0
      SIGY=0.D0

C
C     TRANSFORMANDO A CONVERGENCIA DE (") EM FRACAO DE PIXELS.
C


c     PLIMIT=PLIMIT/ESCL


c     RAIO=1000.D0

      XLAMDA=0.001
      NTERMX=5
      XRESID= -1.D14
      RESIDX=0.D0
      XCENT=1.D14
      YCENT=1.D14


C
C     PARAMETROS DE ENTRADA PARA GAUSSIANA BIDIMENSIONAL SIMETRICA
C

      PARAM(5)=fundox

      PARAM(1)=bx
      PARAM(2)=by

      IAAUX=bx
      IAAUY=by
      PARAM(4)=PIXMAT(IAAUX,IAAUY)-fundox

      sigx=param(4)/2.d0
      do i=ny,my
      do 75 j=nx,mx
      if (pixmat(j,i).ge.maximo) go to 75      
      if (pixmat(j,i)-fundox.gt.sigx) go to 80
 75   continue
      enddo
 80   PARAM(3)=dsqrt((j-bx)**2+(i-by)**2)/2.35d0

C
C     INCREMENTOS INICIAIS AOS PARAMENTROS
C
      DO 5050 I=1,NTERMX
 5050 XSIGMA(I)=0.D0
C
      DO 110 I=1,NTERMX
  110 DELTAX(I)=5.D-2*PARAM(I)
C
C     CALCULO DO CENTRO DA GAUSSIANA
C
   14 CALL GAUSIC (ICONTT,RAIO,NX,NY,MX,MY,maximo,NTERMX,PARAM,DELTAX,
     ?XSIGMA,XLAMDA,RESIDX)
C
C     TESTANDO CONVERGENCIA
C
      IF (RESIDX.GE.0) GO TO 362

C 360 WRITE (*,365) IOUT,NUMEST
C 365 FORMAT (5X,'Ajuste impossivel para  objeto ',A12,1X,I3)
      bx=0.d0
      by=0.d0
      maximo=xmaxi
      return
  362 CONTINUE
C

      RESIDX = DSQRT(RESIDX)
      CENTX  = X0 + PARAM(1) * PASSOX
      CENTY  = Y0 + PARAM(2) * PASSOY
      CONVER = DABS(RESIDX*DLIMIT)
      DIFERD = DABS(RESIDX-XRESID)
      DIFPOX = DABS(CENTX-XCENT)
      DIFPOY = DABS(CENTY-YCENT)
      IF ((DIFERD.LT.CONVER).AND.(DIFPOX.LT.PLIMIT).AND.(DIFPOY.LT.
     ?PLIMIT)) GO TO 150
      XRESID = RESIDX
      XCENT=CENTX
      YCENT=CENTY
      ICONTX = ICONTX + 1
      IF (ICONTX.GT.30) GO TO 150

      GO TO 14

C
  150 CONTINUE


c
c     cancela o triming circular com go to 154
c

c     go to 154


C
C     OPCAO TRIMMING CIRCULAR:  CONVERGENCIA 1% DO SIGMA
C
      SIGDEX = PARAM(3)
      CONSIX = 0.01D0*SIGDEX
      AUX=DABS(SIGDEX-TRANSX)
      IF (AUX.LT.CONSIX) GO TO 154
      ICONTT = ICONTT + 1
      IF (ICONTT.EQ.10) GO TO 154
      RAIO=2.5D0*SIGDEX

      if (raio.gt.ramax) raio=ramax

      TRANSX = SIGDEX
      XLAMDA=0.001
      XRESID= -10.D10
      RESIDX=0.D0
      XCENT=1.D14
      YCENT=1.D14
      GO TO 14
  154 CONTINUE
      SIGDEX = XSIGMA(1) * PASSOX * RESIDX*ESCL
      SIGDEY = XSIGMA(2) * PASSOY * RESIDX*ESCL
      RESIDX = RESIDX
      ALTURX = PARAM(4)
      CENTRX = X0 + PARAM(1) * PASSOX
      CENTRY = Y0 + PARAM(2) * PASSOY
      SIGX   = PARAM(3) * ((PASSOX+PASSOY)/2.D0) *ESCL
      FUNDOX = PARAM(5)
C     XINCLI = PARAMX(5)/PASSOX
C     YINCLI = PARAMY(5)/PASSOY
C     CURVAX = PARAMX(6)/PASSOX**2
C     CURVAY = PARAMY(6)/PASSOY**2
C
C     ESCREVENDO OS RESULTADOS NO ARQUIVO DE SAIDA
C     Aqui sao invertidos X e Y
C
C     WRITE (2,97) NUMEST,CENTRY,CENTRX,SIGDEY,SIGDEX,ICONTT,ICONTX,
C    ?NPIX,SEIXOA,SEIXOB,ANGULO
C  97 FORMAT(I4,1X,2(F10.3,1X),2X,2(F6.3,1X),2(I2,1X),3X,I6,3(1X,F6.2))

C

      bx=CENTRX
      by=CENTRY
      maximo=xmaxi


      RETURN
      END
C     SUBROTINA GAUSIC
C
C     PROPOSITO
C       FAZER  UM AJUSTE DE MINIMOS QUADRADOS ITERATIVO  NAO  LINEAR COM
C       UMA  FUNCAO  GAUSSIANA  BIDIMENSIONAL  SIMETRICA  COM   TRIMMING
C       CIRCULAR.
C
C     USO
C       CALL GAUSIC (KEY,RAIO,NX,NY,MX,MY,NTERMS,A,DELTAA,
C          SIGMAA, FLAMDA, CHISQR)
C
C     DESCRICAO DOS PARAMETROS
C       KEY    - 0         -ITERACAO ANTES  DO TRIMMING CIRCULAR
C                1 OU MAIS -ITERACAO DEPOIS DO TRIMMING CIRCULAR
C       PIXMAT - MATRIZ DE PIXELS DA IMAGEM
C       NPTS   - NUMERO DE PARES (X,Y) DE PONTOS DADOS
C       NTERMS - NUMERO DE PARAMETROS
C       A      - VETOR DE PARAMETROS
C       DELTAA - VETOR DE INCREMENTOS PARA OS PARAMETROS A
C       SIGMAA - VETOR DOS DESVIOS PADRAO PARA OS PARAMETROS A
C       FLAMDA - PROPORCAO INCLUIDA DA PROCURA POR GRADIENTE
C       CHISQR - CHI QUADRADO REDUZIDO PARA O AJUSTE
C
C     SUBROTINAS E SUBPROGRAMAS FUNCTION REQUERIDOS
C       FGAUSI (J, I, A)
C          AVALIA A FUNCAO AJUSTADA PARA O JI-ESIMO PIXEL
C       QIQUAD (KEY,RAIO,PIXMAT,NX,NY,MX,MY,FREE,PARAM)
C          AVALIA O CHI QUADRADO DO AJUSTE AOS DADOS
C       GDERIV (J, I, A, DELTAA, DERIV, NTERMS)
C          AVALIA AS DERIVATIVAS DA FUNCAO EM AJUSTE PARA O
C          JI-ESIMO PIXEL COM RESPEITO A CADA PARAMETRO
C       MATINV (ARRAY, NTERMS, DET)
C          INVERTE UMA MATRIZ SIMETRICA BIDIMENSIONAL DE GRAU NTERMS
C          E CALCULA SEU DETERMINANTE
C
C       CIRCUL (RAIO,A(1),A(2),L,I,ICHAVE)
C
C
C     MODIFICACOES PARA FORTRAN II
C       OMITIR ESPECIFICACOES DE DUPLA PRECISAO
C       ADICIONAR SUFIXO F PARA SQRT NAS LINHAS 73, 84, E 103
C
C     MODIFICACOES PARA FORTRAN 77
C
C       NENHUMA
C
C     COMENTARIOS
C       DIMENSIONAMENTO VALIDO ATE 10 TERMOS DE AJUSTE
C       FIXAR FLAMDA = 0.001 NO COMECO DE CADA PROCURA
C
      SUBROUTINE GAUSIC (KEY,RAIO,NX,NY,MX,MY,maximo,NTERMS,A,DELTAA,
     ?SIGMAA,FLAMDA,CHISQR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(5),DELTAA(5),SIGMAA(5),B(5),ALPHA(21,21),BETA(21),
     ?DERIV(5),ARRAY(21,21)
      real*4 pixmat(5001,5001),maximo

      COMMON /A4/PIXMAT 
      COMMON /A6/ALPHA,BETA
C     COMMON /A13/DERIV
      COMMON /A7/ARRAY
      COMMON /A14/IERRO
      IERRO=0
      DET=1.D0
      ICONT=0
      IXC=A(1)
      iconv=20
      jconv=0
C
      CHISQR = 0.D0

C
C        EVALUATE ALPHA AND BETA MATRICES
C
      IF (KEY.GT.0) GO TO 31
C


  131 DO 134 J=1, NTERMS
      BETA(J) = 0.D0
      DO 134 K=1, J
  134 ALPHA(J,K) = 0.D0
  141 DO 150 I=NY,MY
      DO 5001 L=NX,MX
      IF (PIXMAT(L,I).ge.maximo) GO TO 5001
      CALL GDERIV (L, I, A, DELTAA, DERIV)
      if (ierro.eq.1) go to 107
      ICONT=ICONT+1
      DO 146 J=1,NTERMS
      BETA(J)= BETA(J) + (PIXMAT(L,I)-FGAUSI(L,I,A))*DERIV(J)
      if (ierro.eq.1) go to 107
      DO 146 K=1, J
  146 ALPHA(J,K) = ALPHA(J,K) + DERIV(J)*DERIV(K)
 5001 CONTINUE
  150 CONTINUE
C
c     IF (IERRO.EQ.1) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
C
c     FREE = (MX-NX+1)+(MY-NY+1) - NTERMS

      FREE=ICONT-NTERMS

C

      if (free.le.0.d0) go to 107

c     IF (FREE.LE.0.D0) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
C


      GO TO 51
C
   31 DO 34 J=1, NTERMS
      BETA(J) = 0.D0
      DO 34 K=1, J
   34 ALPHA(J,K) = 0.D0
   41 DO 50 I=NY,MY
c     CALL CIRCUL (RAIO,A(1),A(2),IXC,I,ICHAVE)
c     IF (ICHAVE.LT.0) GO TO 50
      DO 5000 L=NX,MX
      IF (PIXMAT(L,I).ge.maximo) GO TO 5000
      CALL CIRCUL (RAIO,A(1),A(2),L,I,ICHAVE)
      IF (ICHAVE.LT.0) GO TO 5000
      CALL GDERIV (L, I, A, DELTAA, DERIV)
      if (ierro.eq.1) go to 107
      ICONT=ICONT+1
      DO 46 J=1,NTERMS
      BETA(J)= BETA(J) + (PIXMAT(L,I)-FGAUSI(L,I,A))*DERIV(J)
      if (ierro.eq.1) go to 107
      DO 46 K=1, J
   46 ALPHA(J,K) = ALPHA(J,K) + DERIV(J)*DERIV(K)
 5000 CONTINUE
   50 CONTINUE
C
c     IF (IERRO.EQ.1) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
C
C
C     CALCULA GRAUS DE LIBERDADE FREE
C

c     COUNT=ICONT
c     FREE=2*DSQRT(COUNT)-NTERMS
      FREE=ICONT-NTERMS
C

      if (free.le.0.d0) go to 107

c     IF (FREE.LE.0.D0) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
C
   51 CONTINUE
      DO 53 J=1,NTERMS
      DO 53 K=1, J
   53 ALPHA(K,J)=ALPHA(J,K)
C
C        EVALUATES CHI SQUARE AT STARTING POINT
C
C  61 DO 62 I=NY,MY
C     DO 5010 J=NX,MX
C5010 FITMAT(J,I) = FGAUSI (J, I, A)
C  62 CONTINUE
   63 CHISQ1=QIQUAD(KEY,RAIO,NX,NY,MX,MY,FREE,A,maximo)
C
      if (ierro.eq.1) go to 107

c     IF (IERRO.EQ.1) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
c


C
C				 
C        INVERT MODIFIED CURVATURE MATRIX TO FIND NEW PARAMETERS
C


 71   DO 74 J=1, NTERMS
      DO 73 K=1, NTERMS
      AUX = ALPHA(J,J)*ALPHA(K,K)
      IF (AUX.lt.0.D0) GO TO 107
   73 ARRAY(J,K)= ALPHA(J,K) / DSQRT (AUX)
   74 ARRAY(J,J) = 1.D0 + FLAMDA
   80 CALL MATINV (NTERMS, DET)
C

      if (ierro.eq.1) go to 107


c     IF (IERRO.EQ.1) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
C
   81 DO 84 J=1, NTERMS
      B(J) = A(J)
      DO 84 K=1, NTERMS
      AUX = ALPHA(J,J)*ALPHA(K,K)
      IF (AUX.lt.0.D0) GO TO 107
   84 B(J) = B(J) + BETA(K)*ARRAY(J,K)/DSQRT (AUX)
C
C        IF CHI SQUARE INCREASED, INCREASE FLAMDA AND TRY AGAIN
C
C  91 DO 92 I=NY,MY
C     DO 5020 J=NX,MX
C5020 FITMAT(J,I) =  FGAUSI (J, I, B)
C  92 CONTINUE
   93 CHISQR = QIQUAD(KEY,RAIO,NX,NY,MX,MY,FREE,B,maximo)
C

      if (ierro.eq.1) go to 107

c     IF (IERRO.EQ.1) THEN
c     IERRO=0
c     GO TO 107
c     ENDIF
C
      jconv=jconv+1
      if (jconv.gt.iconv) go to 107

      IF (CHISQ1 - CHISQR) 95, 101, 101
   95 FLAMDA = 10.D0*FLAMDA
      GO TO 71
C
C        EVALUATE PARAMETERS AND UNCERTAINTIES
C
  101 DO 104 J=1, NTERMS
      A(J) = B(J)
C     IF (MODE) 103, 102, 103
C 102 SIGMAA(J) = DSQRT (CHISQR*ARRAY(J,J) / ALPHA(J,J))
C     GO TO 104
      AUX = ARRAY(J,J)/ALPHA(J,J)
      IF (AUX.lt.0.D0) GO TO 107
  103 SIGMAA(J) = DSQRT (AUX)
  104 CONTINUE
      FLAMDA = FLAMDA/10.D0
      GO TO 110
  107 CHISQR = -1.D0
      IERRO=1
  110 CONTINUE
      RETURN
      END


C     SUBROUTINE GDERIV GAUSSIANA BIDIMENSIONAL + CONSTANTE
C
C     PURPOSE
C       EVALUATE DERIVATIVES OF FUNCTION FOR LEAST-SQUARES SEARCH
C       WITH FORM OF A  BIDIMENSIONAL SIMETRIC GAUSSIAN PEAK PLUS
C       A CONSTANT.
C          FGAUSI(J,I,A) = A(4)*EXP(-ZX**2/2 - ZY**2/2) + A(5)
C          WHERE X = J, Y = I AND ZX = (X - A(1))/A(3)
C                                 ZY = (Y - A(2))/A(3)
C
C     USAGE
C       CALL GDERIV (J, I, A, DELTAA, DERIV)
C
C     DESCRIPTION OF PARAMETERS
C       J      - INDEX X OF DATA POINTS FOR INDEPENDENT VARIABLE
C       I      - INDEX Y OF DATA POINTS
C       A      - ARRAY OF PARAMETERS
C       DELTAA - ARRAY OF PARAMETER INCREMENTS
C       NTERMS - NUMBER OF PARAMETERS
C       DERIV  - DERIVATIVES OF FUNCTION
C
C     SUBROUTINES AND FUNCTION SUPPROGRAMS REQUIRED
C       NONE
C
C     MODIFICATIONS FOR FORTRAN II
C       ADD F SUFFIX TO EXP IN STATEMENT 21
C
      SUBROUTINE GDERIV (J, I, A, DELTAA, DERIV)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(5), DELTAA(5), DERIV(5)
C     COMMON /A13/DERIV
      COMMON /A14/IERRO
C
      IF ((A(3).EQ.0.D0).OR.(A(4).EQ.0.D0)) THEN
      IERRO=1
      RETURN
      ENDIF
C
   11 X = J
      Y = I
      ZX = (X - A(1)) / A(3)
      ZY = (Y - A(2)) / A(3)
      ZX2 = ZX**2
      ZY2 = ZY**2
C
C         ANALYTICAL EXPRESSIONS FOR DERIVATIVES
C
c     IF ((ZX2.GT.50.D0).OR.(ZY2.GT.50.D0)) THEN 
c     DERIV(1) = 0.D0
c     DERIV(2) = 0.D0
c     DERIV(3) = 0.D0
c     DERIV(4) = 0.D0
c     ELSE
      FUNCAO   = FGAUSI(J,I,A)
      IF (IERRO.EQ.1) RETURN
      DERIV(1) = (X-A(1))*(FUNCAO-A(5))/A(3)**2
      DERIV(2) = (Y-A(2))*(FUNCAO-A(5))/A(3)**2
      DERIV(3) = (ZX2+ZY2)*(FUNCAO-A(5))/A(3)
      DERIV(4) = (FUNCAO-A(5))/A(4)
c     ENDIF
C
      DERIV(5) = 1.D0
C
      RETURN
      END
C     FUNCTION FGAUSI  GAUSSIANA BIDIMENSIONAL + FUNDO CONSTANTE
C
C     PROPOSITO
C       EVALUATE TERMS OF FUNCTION FOR NON-LINEAR LEAST-SQUARES SEARCH
C          WITH FORM OF A BIDIMENSIONAL GAUSSIAN PEAK PLUS A CONSTANT.
C          FUNCTN(J,I,A) = A(4)*EXP(-ZX**2/2-ZY**2/2) + A(5)
C          WHERE X = J, Y = I  AND ZX = (X - A(1))/A(3)
C                                  ZY = (Y - A(2))/A(3)
C
C     USAGE
C       RESULT = FGAUSI (J,I,A)
C
C     DESCRIPTION OF PARAMETERS
C       J      - INDEX OF X DATA POINTS
C       I      - INDEX OF Y DATA POINTS
C       A      - ARRAY OF PARAMETERS
C
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C       NONE
C
C     MODIFICATIONS FOR FORTRAN II
C       ADD F SUFFIX TO EXP IN STATEMENT 16
C
      DOUBLE PRECISION FUNCTION FGAUSI (J, I, A)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(5)
      COMMON /A14/IERRO
C
      IF ((A(3).EQ.0.D0).OR.(A(4).EQ.0.D0)) THEN
      IERRO=1
      RETURN
      ENDIF
C
   11 X = J
      Y = I
   12 FGAUSI = A(5)
   13 ZX = (X - A(1))/A(3)
      ZY = (Y - A(2))/A(3)
      ZX2 = ZX**2
      ZY2 = ZY**2
c     IF (ZX2 - 50.D0) 16, 20, 20
c  16 IF (ZY2 - 50.D0) 17, 20, 20
   17 FGAUSI = FGAUSI + A(4)*DEXP(-ZX2/2.D0 - ZY2/2.D0)
   20 RETURN
      END
C
C
C       FUNCTION QIQUAD
C
C       PURPOSE
C         EVALUATE REDUCED CHI SQUARE FOR FIT TO DATA
C            QIQUAD = SUM ((PIXMAT-FITMAT)**2 / SIGMA**2) / NFREE
C
C       USAGE
C         RESULT = QIQUAD (KEY,RAIO,NX,NY,MX,MY,FREE,PARAM)
C
C       DESCRIPTION OF PARAMETERS
C         KEY    - 0         - SEM TRIMMING
C                  1 OU MAIS - TRIMMING CIRCULAR
C         PIXMAT - MATRIX ARRAY OF DATA POINTS
C         NPTS   - NUMBER OF DATA POINTS
C         FREE   - NUMBER OF DEGREES OF FREEDOM
C         PARAM  - PARAMETROS DA GAUSSIANA AJUSTADA
C
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C
C       CIRCUL (RAIO,A(1),A(2),J,I,ICHAVE)
C
C     MODIFICATIONS FOR FORTRAN II
C       OMIT DOUBLE PRECISION SPECIFICATIONS
C
      DOUBLE PRECISION FUNCTION QIQUAD (KEY,RAIO,NX,NY,MX,MY,FREE,A,
     ?maximo)
      IMPLICIT REAL *8 (A-H,O-Z)
      real*4 pixmat(5001,5001),maximo

      DIMENSION A(5)
      COMMON /A4/PIXMAT
      COMMON /A14/IERRO
C
      CHISQ = 0.D0
      QIQUAD=0.D0
C
      IF (FREE.LE.0.D0) THEN
      IERRO=1
      RETURN
      ENDIF
C
      IF (KEY.GT.0) GO TO 30
      DO 38 I = NY,MY
      DO 38 J = NX,MX
      IF (PIXMAT(J,I).ge.maximo) GO TO 38
      CHISQ = CHISQ + (PIXMAT(J,I)-FGAUSI(J,I,A))**2
      if (ierro.eq.1) return
 38   continue
      GO TO 39
C
   30 IXC=A(1)
      DO 35 I = NY,MY
      CALL CIRCUL (RAIO,A(1),A(2),IXC,I,ICHAVE)
      IF (ICHAVE.LT.0) GO TO 35
      DO 350 J = NX,MX
      IF (PIXMAT(J,I).ge.maximo) GO TO 350
      CALL CIRCUL (RAIO,A(1),A(2),J,I,ICHAVE)
      IF (ICHAVE.LT.0) GO TO 350
      CHISQ = CHISQ + (PIXMAT(J,I)-FGAUSI(J,I,A))**2
      if (ierro.eq.1) return
  350 CONTINUE
   35 CONTINUE
C
C        DIVIDE BY NUMBER OF DEGREES OF FREEDOM
C

   39 QIQUAD = CHISQ / FREE
c  40 CONTINUE
      RETURN
      END

C
C     SUBROTINA CIRCUL
C
C     PROPOSITO
C
C       VERIFICA SE UM PONTO (IX,IY) DA MATRIZ IMAGEM  ENCONTRA-SE  DENTRO
C       DE UM CIRCULO DE RAIO "RAIO" E CENTRO (XC,YC)
C
C     USO
C
C       CALL CILCUL (RAIO,XC,YC,J,I,ICHAVE)
C
C     DESCRIPTION OF PARAMETERS
C       RAIO   - RAIO DO CIRCULO EM PIXELS
C       XC     - CENTRO X DO CIRCULO EM PIXELS
C       YC     - CENTRO Y DO CIRCULO EM PIXELS
C       IX     - COORDENADA X DO PIXEL
C       IY     - COORDENADA Y DO PIXEL
C       ICHAVE - +1 -- PIXEL DENTRO DO CIRCULO
C                -1 -- PIXEL  FORA  DO CIRCULO
C
C     SUBROTINAS E SUBPROGRAMAS REQUERIDOS
C       NENHUM
C
C     MODIFICATIONS FOR FORTRAN II
C       OMIT DOUBLE PRECISION SPECIFICATIONS
C       CHANGE DABS TO ABSF IN STATEMENT 23
C
C     COMENTARIOS
C       NENHUM
C
      SUBROUTINE CIRCUL (RAIO,XC,YC,IX,IY,ICHAVE)
      IMPLICIT REAL *8 (A-H,O-Z)
C
      RAIO2=RAIO**2
      RADIUS=(IX-XC)**2+(IY-YC)**2
      IF (RADIUS.LE.RAIO2) THEN
      ICHAVE=1
      ELSE
      ICHAVE=-1
      ENDIF
      RETURN
      END

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
      SUBROUTINE MATINV (NORDER, DET)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION ARRAY (21,21), IK(21), JK(21)
      COMMON /A7/ARRAY
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
      SUBROUTINE ORDEM (IDIM,N,IORDEM,NVAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IORDEM(50000),NVAL(50000)
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
c     Mod.: M. Assafin 10/Nov/2006
c
c

      subroutine tmass (mraiz,rac,dec,drac,ddec,ra2ma,de2ma,era2ma,
     ?ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,nest)

      implicit real*8 (a-h,o-z)


      INTEGER*4 CO1,CO2,JJD
      INTEGER*2 MAG1,MAG2,MAG3,ERCO
      INTEGER*1 ERMG1,ERMG2,ERMG3


      dimension ra2ma(50000),de2ma(50000),era2ma(50000),ede2ma(50000),
     ?dmgj(50000),dmgh(50000),dmgk(50000),emgj(50000),emgh(50000),
     ?emgk(50000),ddj2(50000)


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
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c


      idin50=50
      izon=1800

c

      xbox=grarad*drac
      ybox=grarad*ddec
      grac=grarad*rac
      gdec=grarad*dec

c
c     Leitura das faixas de declinacao de 0.1 em 0.1 graus
c     do 2mass astrometrico
c

      dfaixa=dec-ddec-0.1d0
      decmax=dec+ddec


c


      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.1d0


      if (dfaixa-decmax.gt.0.1d0) go to 35


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

      RA=DBLE(CO1)/1.0D6
      DE=DBLE(CO2)/1.0D6

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
c     Mod.: M. Assafin  24/Out/2012
c
c

      subroutine stmass (mraiz,rac,dec,drac,ddec,ra2ma,de2ma,era2ma,
     ?ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,nest)

      implicit real*8 (a-h,o-z)


      INTEGER*4 CO1,CO2,JJD
      INTEGER*2 MAG1,MAG2,MAG3,ERCO
      INTEGER*1 ERMG1,ERMG2,ERMG3


      dimension ra2ma(50000),de2ma(50000),era2ma(50000),ede2ma(50000),
     ?dmgj(50000),dmgh(50000),dmgk(50000),emgj(50000),emgh(50000),
     ?emgk(50000),ddj2(50000)

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
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      idin50=50
      izon=1800

c

      xbox=grarad*drac
      ybox=grarad*ddec
      grac=grarad*rac
      gdec=grarad*dec


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
c     Leitura das faixas de declinacao de 0.1 em 0.1 graus
c     do 2mass astrometrico
c

      dfaixa=dec-ddec-0.1d0
      decmax=dec+ddec


c


      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.1d0


      if (dfaixa-decmax.gt.0.1d0) go to 35


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


      rfaixa=jaca*0.1d0+0.1d0
      if (rfaixa.gt.90.d0) rfaixa=90.d0


c     x=grac+xbox
c
c     d=DEXY(x,rfaixa,grac,rfaixa)
c     xx=XPAD(x,rfaixa,grac)/d
c     yy=YPAD(x,rfaixa,grac,rfaixa)/d
c
c     xx=dabs(xx)
c     yy=dabs(yy)
c
c     xx1=alff(-xx,yy,grac,rfaixa)
c     xx2=alff(+xx,yy,grac,rfaixa)




      yy=0.d0

      call alfa (+xbox,yy,grac,rfaixa,ralfa)


      xx=dabs(ralfa-grac)

      xx1=grac-xx
      xx2=grac+xx



      xx1=xx1*radgra
      xx2=xx2*radgra


      if (xx1.lt.0d0)  xx1=xx1+360.d0
      if (xx1.gt.360.d0) xx1=xx1-360.d0


      if (xx2.lt.0d0)  xx2=xx2+360.d0
      if (xx2.gt.360.d0) xx2=xx2-360.d0


      if (xx2.ge.xx1) then

      inde=1

      jxmin(1)=xx1+1
      jxmax(1)=xx2+1

      cxmin(1)=xx1
      cxmax(1)=xx2


c     if (jxmax(1).gt.360) jxmax(1)=360

      else

      inde=2

      jxmin(1)=xx1+1
      jxmax(1)=360

      jxmin(2)=1
      jxmax(2)=xx2+1

      cxmin(1)=xx1
      cxmax(1)=360.d0
      cxmin(2)=0.d0
      cxmax(2)=xx2


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

      RA=DBLE(CO1)/1.0D6

c

      if (ra.lt.cxmin(in)) go to 18
      if (ra.gt.cxmax(in)) go to 18

c


      DE=DBLE(CO2)/1.0D6
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



      return
      end





c
c
c
c     Subrotina idxy2m
c
c
c     Identifica estrelas medidas (x,y) com estrelas do catalogo 2MASS
c
c     - id2ma: carrega para cada estrela medida no CCD o numero da n-esima 
c       estrela 2MASS identificada  
c
c     - ialtu: as N (N=nboff) ultimas posicoes no vetor sao as N
c       estrelas mais brilhantes medidas
c
c     - ra2ma, de2ma: alfas e deltas do 2MASS
c
c
c     - xob,yob: medidas (x,y) das estrelas do campo CCD
c
c
c     - xold,yold: deposito provisorio da projecao do 2MASS no plano
c       tangente
c
c     - xra2ma, yde2ma: alfas e deltas provisorios das estrelas medidas
c       para posterior identificacao com o UCAC2
c
c     - scala: escala de pixel, atualizada ao final da subrotina
c
c     - areax, areay : metade da area CCD para extracao das estrelas
c                      2MASS brilhantes (rad)
c
c
c      Mod. M. Assafin   08/Nov/2012
c
c

      subroutine idxy2m (areax,areay,idiobs,nval,ior,scala,erpix,rac,
     ?dec,nbcat,nbmed,nest,ialtu,xob,yob,n2mass,id2ma,ra2ma,de2ma,dmgj,
     ?xold,yold,xra2ma,yde2ma,ireflex,ecala,tt,coefx,coefy,ncomum)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension ialtu(50000),xob(50000),yob(50000),id2ma(50000),
     ?ra2ma(50000),de2ma(50000),dmgj(50000),xra2ma(50000),
     ?yde2ma(50000),nval(50000),ior(50000),xold(50000),yold(50000)

      dimension coefx(21),coefy(21),xcof(21),ycof(21),xest(50000),
     ?yest(50000),xp(50000),yp(50000)

      character*1 isiga

      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Dados auxiliares
C
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      ierro=0

      ncomum=0

      pcala=grarad*ecala/3600.d0

c     epixra=grarad*erpix*scala/3600.d0

      epixra=grarad*erpix/3600.d0

      rscala=grarad*scala/3600.d0




c
c     Projecao do 2MASS no plano tangente
c

      grac=grarad*rac
      gdec=grarad*dec

      n=0

      do 5 i=1,n2mass
      bra=grarad*ra2ma(i)
      bde=grarad*de2ma(i)
      d=DEXY(bra,bde,grac,gdec)

      xold(i)=XPAD(bra,bde,grac)/d
      yold(i)=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xold(i)).gt.areax) go to 5
      if (dabs(yold(i)).gt.areay) go to 5

      n=n+1
      ior(n)=i
      nval(n)=dmgj(i)*100

 5    continue




c
c     Separa as N (ex.:N=30) estrelas mais brilhantes do 2MASS
c     Sao as N primeiras posicoes no vetor ior
c


      call ordem (idiobs,n,ior,nval)

c

      nmed=nbmed
      if (nbmed.gt.nest) nmed=nest


      ncat=nbcat
      if (nbcat.gt.n2mass) ncat=n2mass

      if (nbcat.gt.n) ncat=n

c     write (*,*) 'n_bri = ',n

c
c     Identificacao par a par, com e sem reflexao
c
c     k=1 : sem reflexao em x
c     k=-1: com reflexao em x
c


      do i=1,21
      xcof(i)=0.d0
      ycof(i)=0.d0
      coefx(i)=0.d0
      coefy(i)=0.d0
      enddo


      ireflex=0
      ncomum=0

      nptos=2
      ngrau=0
      ngraup=1

c
c     Inicializando tempo de execucao do programa para identificacao das 
c     estrelas
c

c     ntempo=ncat*0.1d0

c     tempoi=0.d0
c     call tempo (tempoi,tempot,tempop)

c



      do iii=1,ncat-1

c     if (k.eq.-1 .and. iii.eq.ntempo) then

c     call tempo (tempoi,tempot,tempop)

c     tt=2.d0*ncat*tempop/ntempo
c     tt=tt/3600.d0
      
c     ihour=tt
c     minu=(tt-ihour)*60.d0
c     seg=((tt-ihour)*60.d0-minu)*60.d0
c
c     write (*,*)      
c     write (*,1) ihour,minu,seg
c1    format(1x,'Time consuming per field for reference catalog star ide
c    ?ntification: ',i3,'hs ',i2,'m ',f4.1,'s')  
c     write (*,*)      
c
c     endif

c

      xp(1)=xold(ior(iii))
      yp(1)=yold(ior(iii))

      do ii =iii+1,ncat


      xp(2)=xold(ior(ii))     
      yp(2)=yold(ior(ii))     


      do jjj=1,nmed

      xest(1)=xob(ialtu(nest-jjj+1))
      yest(1)=yob(ialtu(nest-jjj+1))

      do 20 jj =1,nmed

      if (jj.eq.jjj) go to 20

      xest(2)=xob(ialtu(nest-jj+1))     
      yest(2)=yob(ialtu(nest-jj+1))     

c
c     Elimina par se distancia entre estrelas de cada par nao e'
c     compativel com a escala em pixels
c

      dp=dsqrt((xp(2)-xp(1))**2+(yp(2)-yp(1))**2)
      dest=dsqrt((xest(2)-xest(1))**2+(yest(2)-yest(1))**2)
      d=dp/dest



      if (dabs(d-rscala).gt.pcala) go to 20

c
c     Calcula coeficientes de 4ctes para os 2 pares
c


      do k=+1,-1,-2

      xest(1)=k*xest(1)
      xest(2)=k*xest(2)

      call isol (ngrau,nptos,xest,yest,xp,yp,coefx,coefy)

      if (ierro.eq.1) then
      ierro=0
      go to 20

      endif


c
c     Computa o numero de identificacoes para esses 2 pares
c

      icont=0

      do 10 j=1,nmed

      xx=k*xob(ialtu(nest-j+1))
      yy=yob(ialtu(nest-j+1))

      x=pol(xx,yy,coefx,ngraup)
      y=pol(xx,yy,coefy,ngraup)

      do  9 i=1,ncat
      
      dx=dabs(xold(ior(i))-x)
      dy=dabs(yold(ior(i))-y)
      

      if (dx.gt.epixra) go to 9
      if (dy.gt.epixra) go to 9

      icont=icont+1
      go to 10



 9    continue
 10   continue


      if (icont.gt.ncomum) then

      ncomum=icont

      do l=1,3
      xcof(l)=coefx(l)
      ycof(l)=coefy(l)
      enddo

      ireflex=k
      endif

      enddo


 20   continue
      enddo


      enddo
      enddo






c
c     debug
c

c     write (*,*) 
c     write (*,*) '(RA,DE) vs. (x,y) cross-identification'
c     write (*,*) '1rst step: 4 constant adjust'
c     write (*,*) 
c
c     write (*,*) '2MASS stars extracted in 2x2 field  = ',n2mass
c     write (*,*) 'Extracted field stars with (x,y)    = ',nest
c     write (*,*) 'Common 2MASS vs. bright field stars = ',ncomum
c     write (*,*) 'system reflection flag              = ',ireflex
c     write (*,*)
c
c
c
c     x1=xcof(1)*radgra*3600.d0
c     x2=xcof(2)*radgra*3600.d0
c     x3=xcof(3)*radgra*3600.d0
c
c     y1=ycof(1)*radgra*3600.d0
c     y2=ycof(2)*radgra*3600.d0
c     y3=ycof(3)*radgra*3600.d0
c
c
c     write (*,*) 'Coeficients solution'
c     write (*,*) 'X = A + Bx + Cy '
c     write (*,*) 'Y = D - Cx + By '
c     write (*,*)
c     write (*,*) 'X: ',x1,x2,x3
c     write (*,*) 'Y: ',y1,y2,y3
c     write (*,*)
c

c
c     Refina o ajuste polinomial a partir do uso de mais estrelas comuns
c     ao 2MASS, obtidas a partir do resultado com os 2 pares; o polinomio
c     usado agora e' de grau completo, 1, 2 ou 3 conforme o numero de
c     estrelas comuns
c

      ncomum=0

      do i=1,idiobs
      id2ma(i)=0
      enddo


      do 30 j=1,nest
      
      xx=ireflex*xob(j)
      yy=yob(j)

      x=pol(xx,yy,xcof,ngraup)
      y=pol(xx,yy,ycof,ngraup)

      do i=1,n2mass

      dx=(xold(i)-x)*radgra*3600d0
      dy=(yold(i)-y)*radgra*3600d0

      nval(i)=dsqrt(dx**2+dy**2)

      ior(i)=i
      enddo

      call ordem (idiobs,n2mass,ior,nval)

      ii=ior(1)

      dx=dabs(xold(ii)-x)
      dy=dabs(yold(ii)-y)


      if (dx.gt.epixra) go to 30
      if (dy.gt.epixra) go to 30

c

      ncomum=ncomum+1

      id2ma(j)=ii
c     nval(ncomum)=ii
      xp(ncomum)=xold(ii)
      yp(ncomum)=yold(ii)

      xest(ncomum)=ireflex*xob(j)
      yest(ncomum)=yob(j)


 30   continue


c

      if (ncomum.lt.3) then
      ierro=1
      return
      endif

c
c     Adota modelos da mais alta ordem conforme o numero de estrelas 
c     comuns 
c

      ngraup=2
      if (ncomum.le.11) ngraup=1
      if (ncomum.ge.21) ngraup=3


      call isol (ngraup,ncomum,xest,yest,xp,yp,coefx,coefy)

      if (ierro.eq.1) then

      return
      endif

c

c     write (*,*) 
c     write (*,*) '(RA,DE) vs. (x,y) cross-identification'
c     write (*,*) '2nd final step:'
c     write (*,*) 'Complete polynomial model, degree = ',ngraup
c     write (*,*)  
c
c     write (*,*) '2MASS stars extracted in 2x2 field   = ',n2mass
c     write (*,*) 'Extracted field stars with (x,y)     = ',nest
c     write (*,*) 'Common 2MASS vs. field stars (final) = ',ncomum
c     write (*,*) 'system reflection flag             = ',ireflex
c     write (*,*)

c


c     x1=coefx(1)*radgra*3600.d0
c     x2=coefx(2)*radgra*3600.d0
c     x3=coefx(3)*radgra*3600.d0
c
c     y1=coefy(1)*radgra*3600.d0
c     y2=coefy(2)*radgra*3600.d0
c     y3=coefy(3)*radgra*3600.d0
c


c
c     write (*,*) 'Solution. Linear coeficients:'
c     write (*,*) 'X = A + Bx + Cy '
c     write (*,*) 'Y = D + Ex + Fy '
c     write (*,*)
c     write (*,*) 'X: ',x1,x2,x3
c     write (*,*) 'Y: ',y1,y2,y3
c     write (*,*)


c
c     Atualizando a escala de pixel
c






      return

      end






c
c
c
c     Subrotina mdxy2m
c
c
c     Identifica estrelas medidas (x,y) com estrelas do catalogo 2MASS
c
c     - id2ma: carrega para cada estrela medida no CCD o numero da n-esima 
c       estrela 2MASS identificada  
c
c     - ialtu: as N (N=nboff) ultimas posicoes no vetor sao as N
c       estrelas mais brilhantes medidas
c
c     - ra2ma, de2ma: alfas e deltas do 2MASS
c
c
c     - xob,yob: medidas (x,y) das estrelas do campo CCD
c
c
c     - xold,yold: deposito provisorio da projecao do 2MASS no plano
c       tangente
c
c     - xra2ma, yde2ma: alfas e deltas provisorios das estrelas medidas
c       para posterior identificacao com o UCAC2
c
c     - rcala: estimativa melhorada da escala de pixel (rad/pixel)
c
c     - ecala: estimativa melhorada do erro da escala de pixel (rad/pixel)
c
c
c     - areax, areay : metade da area CCD para extracao das estrelas
c                      2MASS brilhantes (rad)
c
c
c
c
c     Nesta subrotina, as orientacoes e escalas sao fornecidas apartir da reducao
c     da primeira imagem, tornando a identificacao das demais imagens acelerada.
c
c
c
c      Modificada:  M. Assafin  08/Nov/2012
c



      subroutine mdxy2m (areax,areay,ngrau,coefxr,coefyr,kreflex,idiobs,
     ?nval,ior,rscala,erpix,rac,dec,nbcat,nbmed,nest,ialtu,xob,yob,
     ?n2mass,id2ma,ra2ma,de2ma,dmgj,xold,yold,xra2ma,yde2ma,ireflex,
     ?pcala,tt)



      IMPLICIT REAL*8 (A-H,O-Z)

      dimension ialtu(50000),xob(50000),yob(50000),id2ma(50000),
     ?ra2ma(50000),de2ma(50000),dmgj(50000),xra2ma(50000),
     ?yde2ma(50000),nval(50000),ior(50000),xold(50000),yold(50000)

      dimension coefx(21),coefy(21),coefxr(21),coefyr(21),xcof(21),
     ?ycof(21),xest(50000),yest(50000),xp(50000),yp(50000)

      character*1 isiga

      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Dados auxiliares
C
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      ierro=0

      ncof=21

      
      ireflex=kreflex

c

      epixra=erpix*grarad/3600.d0

c

      do i=1,ncof
      xcof(i)=0.d0
      ycof(i)=0.d0
      coefx(i)=0.d0
      coefy(i)=0.d0
      enddo

c
c     Projecao do 2MASS no plano tangente
c

      grac=grarad*rac
      gdec=grarad*dec


      n=0

      do 5 i=1,n2mass
      bra=grarad*ra2ma(i)
      bde=grarad*de2ma(i)
      d=DEXY(bra,bde,grac,gdec)

      xold(i)=XPAD(bra,bde,grac)/d
      yold(i)=YPAD(bra,bde,grac,gdec)/d


      if (dabs(xold(i)).gt.areax) go to 5
      if (dabs(yold(i)).gt.areay) go to 5

      n=n+1
      ior(n)=i
      nval(n)=dmgj(i)*100

 5    continue



c
c     Separa as N (ex.:N=30) estrelas mais brilhantes do 2MASS
c     Sao as N primeiras posicoes no vetor ior
c


      call ordem (idiobs,n,ior,nval)

c

      nmed=nbmed
      if (nbmed.gt.nest) nmed=nest


      ncat=nbcat
      if (nbcat.gt.n2mass) ncat=n2mass


      if (nbcat.gt.n) ncat=n



c
c     Inicializando tempo de execucao do programa para identificacao das 
c     estrelas
c


      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)


c
c     Identificacao par a par, com orientacao e escala pre-estabelecidas
c     da primeira imagem
c

      ncomum=0


c

      k=kreflex

c

      do iii=1,ncat-1


      xp(1)=xold(ior(iii))
      yp(1)=yold(ior(iii))


      do ii =iii+1,ncat


      xp(2)=xold(ior(ii))     
      yp(2)=yold(ior(ii))     


      do jjj=1,nmed

      xest(1)=k*xob(ialtu(nest-jjj+1))
      yest(1)=yob(ialtu(nest-jjj+1))

      do 20 jj =1,nmed

      if (jj.eq.jjj) go to 20

      xest(2)=k*xob(ialtu(nest-jj+1))     
      yest(2)=yob(ialtu(nest-jj+1))     

c
c     Elimina par se distancia entre estrelas de cada par nao e'
c     compativel com a escala em pixels
c

      dp=dsqrt((xp(2)-xp(1))**2+(yp(2)-yp(1))**2)
      dest=dsqrt((xest(2)-xest(1))**2+(yest(2)-yest(1))**2)
      d=dp/dest


      if (dabs(d-rscala).gt.pcala) go to 20



c
c     Adota coeficientes da primeira imagem; apenas o termo independente eh 
c     computado
c


      coefxr(1)=0.d0
      coefyr(1)=0.d0

      cx1=xp(1)-pol(xest(1),yest(1),coefxr,ngrau)
      cy1=yp(1)-pol(xest(1),yest(1),coefyr,ngrau)
      cx2=xp(2)-pol(xest(2),yest(2),coefxr,ngrau)
      cy2=yp(2)-pol(xest(2),yest(2),coefyr,ngrau)


      coefxr(1)=(cx1+cx2)/2.d0
      coefyr(1)=(cy1+cy2)/2.d0


c
c     Computa o numero de identificacoes para esses 2 pares
c

      icont=0

      do 10 j=1,nmed

      xx=k*xob(ialtu(nest-j+1))
      yy=yob(ialtu(nest-j+1))

      x=pol(xx,yy,coefxr,ngrau)
      y=pol(xx,yy,coefyr,ngrau)

      do  9 i=1,ncat
      
      dx=dabs(xold(ior(i))-x)
      dy=dabs(yold(ior(i))-y)
      

      if (dx.gt.epixra) go to 9
      if (dy.gt.epixra) go to 9

      icont=icont+1
      go to 10



 9    continue
 10   continue


      if (icont.gt.ncomum) then

      ncomum=icont
      do l=1,ncof
      xcof(l)=coefxr(l)
      ycof(l)=coefyr(l)
      enddo

      endif


 20   continue
      enddo


      enddo
      enddo


c     stop

c

c     call tempo (tempoi,tempot,tempop)

c     tt=tempop
c     tt=tt/3600.d0
c     
c     ihour=tt
c     minu=(tt-ihour)*60.d0
c     seg=((tt-ihour)*60.d0-minu)*60.d0

c     write (*,*)      
c     write (*,1) ihour,minu,seg
c1    format(1x,'Time consuming per field for reference catalog star ide
c    ?ntification: ',i3,'hs ',i2,'m ',f4.1,'s')  
c     write (*,*)      




c
c     debug
c
c
c     write (*,*) 
c     write (*,*) '(RA,DE) vs. (x,y) cross-identification'
c     write (*,*) '1rst step: coefficients from 1rst image solution'
c     write (*,*) 
c
c     write (*,*) '2MASS stars extracted in 2x2 field  = ',n2mass
c     write (*,*) 'Extracted field stars with (x,y)    = ',nest
c     write (*,*) 'Common 2MASS vs. bright field stars = ',ncomum
c     write (*,*) 'system reflection flag              = ',kreflex
c     write (*,*)

c

c     x1=xcof(1)*radgra*3600.d0
c     x2=xcof(2)*radgra*3600.d0
c     x3=xcof(3)*radgra*3600.d0
c
c     y1=ycof(1)*radgra*3600.d0
c     y2=ycof(2)*radgra*3600.d0
c     y3=ycof(3)*radgra*3600.d0

c
c     write (*,*) 'Solution. Linear coeficients:'
c     write (*,*) 'X = A + Bx + Cy '
c     write (*,*) 'Y = D + Ex + Fy '
c     write (*,*)
c     write (*,*) 'X: ',x1,x2,x3
c     write (*,*) 'Y: ',y1,y2,y3
c     write (*,*)


c
c     Refina o ajuste polinomial a partir do uso de mais estrelas comuns
c     ao 2MASS, obtidas a partir do resultado com os 2 pares; o polinomio
c     usado e' de mesmo grau da primeira imagem de referencia
c

      ncomum=0

      do i=1,idiobs
      id2ma(i)=0
      enddo


      do 30 j=1,nest
      
      xx=ireflex*xob(j)
      yy=yob(j)

      x=pol(xx,yy,xcof,ngrau)
      y=pol(xx,yy,ycof,ngrau)

      do i=1,n2mass

      dx=(xold(i)-x)*radgra*3600d0
      dy=(yold(i)-y)*radgra*3600d0

      nval(i)=dsqrt(dx**2+dy**2)

      ior(i)=i
      enddo

      call ordem (idiobs,n2mass,ior,nval)

      ii=ior(1)

      dx=dabs(xold(ii)-x)
      dy=dabs(yold(ii)-y)


      if (dx.gt.epixra) go to 30
      if (dy.gt.epixra) go to 30




      ncomum=ncomum+1

      id2ma(j)=ii
      xp(ncomum)=xold(ii)
      yp(ncomum)=yold(ii)

      xest(ncomum)=ireflex*xob(j)
      yest(ncomum)=yob(j)


 30   continue



      if (ncomum.lt.3) then
      ierro=1
      return
      endif

c
c     Adota modelos da mais alta ordem conforme o numero de estrelas 
c     comuns 
c


      call isol (ngrau,ncomum,xest,yest,xp,yp,coefx,coefy)

      if (ierro.eq.1) then
      return
      endif

c

c     write (*,*) 
c     write (*,*) '(RA,DE) vs. (x,y) cross-identification'
c     write (*,*) '2nd final step:'
c     write (*,*) 'Complete polynomial model, degree = ',ngrau
c     write (*,*)  
c
c     write (*,*) '2MASS stars extracted in 2x2 field   = ',n2mass
c     write (*,*) 'Extracted field stars with (x,y)     = ',nest
c     write (*,*) 'Common 2MASS vs. field stars (final) = ',ncomum
c     write (*,*) 'system reflection flag             = ',ireflex
c     write (*,*)

c


c     x1=coefx(1)*radgra*3600.d0
c     x2=coefx(2)*radgra*3600.d0
c     x3=coefx(3)*radgra*3600.d0

c     y1=coefy(1)*radgra*3600.d0
c     y2=coefy(2)*radgra*3600.d0
c     y3=coefy(3)*radgra*3600.d0



c
c     write (*,*) 'Solution. Linear coeficients:'
c     write (*,*) 'X = A + Bx + Cy '
c     write (*,*) 'Y = D + Ex + Fy '
c     write (*,*)
c     write (*,*) 'X: ',x1,x2,x3
c     write (*,*) 'Y: ',y1,y2,y3
c     write (*,*)


c
c     Calculando alfas e deltas provisorios no sistema 2MASS para
c     todas as estrelas do campo para reconhecimento com outros
c     catalogos (UCAC2, etc)
c

      do i=1,nest


      xx=ireflex*xob(i)
      yy=yob(i)

      x=pol(xx,yy,coefx,ngrau)
      y=pol(xx,yy,coefy,ngrau)

      xra2ma(i)=alff(x,y,grac,gdec)
      yde2ma(i)=deltt(xra2ma(i),y,grac,gdec)


      xra2ma(i)=xra2ma(i)*radgra
      yde2ma(i)=yde2ma(i)*radgra


c
c     debug alfa e delta
c

c     ra=xra2ma(i)
c     de=yde2ma(i)


c     ra=ra/15.d0
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIGA='-'
c     de=-de
c     ELSE
c     ISIGA='+'
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0 

c     write(*,*) IAH,IAM,SA,'  ',ISIGA,IDG,IDM,DS



      enddo


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
      SUBROUTINE ISOLO (NGRAU,NCOMUM,XESTO,YESTO,PO,COFO)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COFO(136),ALPHAO(136,136),ARRAYO(136,136),
     ?BETAO(136),TERMO(136)
      DIMENSION XESTO(25010001),YESTO(25010001),PO(25010001)

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

      idim=25010001
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

      if (ngrau.lt.1 .or. ngrau.gt.idimg) then
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
      SUBROUTINE ISOL (NGRAU,NCOMUM,XEST,YEST,XP,YP,COEFX,COEFY)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(21),COEFY(21),ALPHA(21,21),ARRAY(21,21),
     ?BETA(21),BETAX(21),TERMX(21)
      DIMENSION XEST(50000),YEST(50000),XP(50000),YP(50000)

      COMMON /A6/ALPHA,BETA
      COMMON/A7/ARRAY
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

      DO 8 I=1,21
      BETA(I) =0.D0
      BETAX(I) =0.D0
      COEFX(I) =0.D0
      COEFY(I) =0.D0
      TERMX(I) =0.D0
      DO 8 J=I,21
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
      CALL MATINV (ITERMS,DET)
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
      CALL MATINV (ITERMS,DET)
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

      double precision function pol (x, y, coefis, ngrau)

      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION coefis(21)

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
c     Atualizacao: M. Assafin  22/Dez/2005
c
c
c



      subroutine posred (ireflex,rac,dec,id,ncat,racat,decat,nest,xob,
     ?yob,corte,ngrau,ngrau3,ngrau5,nstart,nfinal,ra,de,era,ede,
     ?alfsig,delsig,alfres,delres,coefx,coefy,ecoefx,ecoefy,itira,
     ?avam,dvam)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension id(50000),racat(50000),decat(50000),xob(50000),
     ?yob(50000),xp(50000),yp(50000),xest(50000),yest(50000)

      dimension ra(50000),de(50000),era(50000),ede(50000),coefx(21),
     ?coefy(21),ecoefx(21),ecoefy(21),alfres(50000),delres(50000),
     ?itira(50000),xsao(21),ysao(21),xrray(21,21),yrray(21,21),
     ?array(21,21)


      COMMON /A7/ARRAY
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

      idiobs=50000
      ncof=21

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


      nstart=0
      grac=grarad*rac
      gdec=grarad*dec

      do 10 i=1,nest

      if (id(i).eq.0) go to 10

      nstart=nstart+1

      xest(nstart)=ireflex*xob(i)
      yest(nstart)=yob(i)

c
c     Projecao do catalogo de referencia no plano tangente
c


      bra=grarad*racat(id(i))
      bde=grarad*decat(id(i))
      d=DEXY(bra,bde,grac,gdec)

      xp(nstart)=xpad(bra,bde,grac)/d
      yp(nstart)=ypad(bra,bde,grac,gdec)/d


 10   continue


c
c     Ajuste do modelo polinomial entre (x,y) e (X,Y)
c


      call solucao (ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,yp,ntira,
     ?coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,corte,
     ?xrray,yrray)


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

      nfinal=nstart-ntira

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

      subroutine solucao (ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,yp,
     ?ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,
     ?corte,xrray,yrray)




      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(21),COEFY(21),ALPHAX(21,21),ALPHAY(21,21),
     ?ARRAY(21,21),BETAX(21),BETAY(21),TERMX(21),TERMY(21),
     ?ITIRA(50000)

      DIMENSION XEST(50000),YEST(50000),XP(50000),YP(50000),
     ?XRRAY(21,21),YRRAY(21,21),ALFRES(50000),DELRES(50000)

      COMMON /A7/ARRAY
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
      NTIRA=0
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
      CALL MATINV (ITERMS,DET)
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

      CALL MATINV (ITERMS,DET)
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

      subroutine estat (box,ialvos,input,itotal,iresum)


      IMPLICIT REAL *8 (A-H,O-Z)

      character*150 infits 
      character*50 ialvos,input,itotal,iresum
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

      write (*,40) iresum
 40   format(a50)


      write (*,*)
      write (*,*)'offsets (RA,DE), Sigma(RA,DE), Ncat, Gauss_error(x,y),
     ?date, exptime, filter, target'
      write (*,*)

c

      open (7,file=ialvos,status='old',err=200)

c
      open (2,file=itotal)
      open (10,file=iresum)

c

 3    read (2,5,end=6) isig
 5    format(a1)
      go to 3
 6    call backsp (2,nbac,2)

c

 7    read (10,5,end=9) isig
      go to 7
 9    call backsp (2,nbac,10)


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

c     write (*,16) dx,dy,alfsiu,delsiu,nfinau,erau,edeu,ex,ey,dj,iuth,
c    ?iutm,sut,iutano,iutmes,iutdia,iexps,ichfil,iobalv,nx,ny

c16   format(4(1x,f7.3),2x,i5,2x,4(1x,f7.3),3x,f16.8,2x,i2,1x,i2,
c    ?1x,f5.2,1x,i4,1x,i2,1x,i2,2x,i4,3x,a20,2x,a20,2(1x,i5))

      write (*,16) dx,dy,alfsiu,delsiu,nfinau,ex,ey,iuth,
     ?iutm,sut,iutano,iutmes,iutdia,iexps,ichfil,iobalv

 16   format(4(1x,f7.3),2x,i5,2x,2(1x,f7.3),2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,2x,i4,3x,a20,2x,a20)

     
      write (10,17) dx,dy,alfsiu,delsiu,nfinau,erau,edeu,ex,ey,dj,iuth,
     ?iutm,sut,iutano,iutmes,iutdia,iexps,ichfil,infits,iobalv,nx,ny

 17   format(4(1x,f7.3),2x,i5,2x,4(1x,f7.3),3x,f16.8,2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,2x,i4,3x,a20,2x,a50,1x,a20,2(1x,i5))


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
      close (10)

      if (ipegou.eq.0) then
      write (*,40) iresum
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
c
c     subroutine swapo
c
c
c     Swap bytes of pixel data
c
c     Image can by integer or floating point
c
c
c     Last modified: M. Assafin   18/Dec/2009
c


      subroutine swapo (if,ibytes,nbytes,irec,swork)

      IMPLICIT REAL *8 (A-H,O-Z)
      
      integer*1 swork(2880),iby4(4)

c

      read (if,rec=irec) swork


      do k=ibytes,nbytes,ibytes

      do m=1,ibytes
      iby4(m)=swork(k-m+1)
      enddo

      do m=1,ibytes
      swork(k-ibytes+m)=iby4(m)
      enddo

      enddo


c

      return
      end





c
c     Subroutine flatsk
c
c
c     Flattens sky background
c
c
c     last update:  M. Assafin 12/Mar/2010
c
c


      subroutine flatsk (nx,ny,ngrauf,mazi,fcmin,fcmax)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi

      dimension xesto(25010001),yesto(25010001),po(25010001),cofo(136)
      dimension pixmat(5001,5001)


      common /a4/pixmat

c


      idimo=25010001
      idimog=136

c

      do i=1,idimo
      xesto(i)=0.d0
      yesto(i)=0.d0
      po(i)=0.d0
      enddo

      do i=1,idimog
      cofo(i)=0.d0
      enddo

c


      n=0
      xmean=0.d0
      do 2 i=1,ny
      do 1 j=1,nx
      if (pixmat(j,i).gt.mazi) go to 1
      if (pixmat(j,i).gt.fcmax) go to 1
      if (pixmat(j,i).lt.fcmin) go to 1
      n=n+1
      xesto(n)=j
      yesto(n)=i
      po(n)=pixmat(j,i)
      xmean=xmean+po(n)
   1  continue
   2  continue

      call isolo (ngrauf,n,xesto,yesto,po,cofo)

c
c     Normalizes polynomial solution
c

      xmean=xmean/n


c
c     Flattens sky
c


      do 4 i=1,ny
      do 3 j=1,nx
      xx=j
      yy=i
      if (pixmat(j,i).gt.mazi) go to 3
      pixmat(j,i)=xmean*pixmat(j,i)/polo(xx,yy,cofo,ngrauf)
   3  continue
   4  continue


      return
      end



c
c
c     Subroutine skyb
c
c
c     Determines sky background by the mode of the histogram of counts.
c
c
c     last update:  M. Assafin 12/Mar/2010
c
c


      subroutine skyb (nx,ny,mazi,malisa,vmin,fatceu,ceu,sigceu,izmax,
     ?threso)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi

      dimension ifundo(150000),jfundo(150000)
      dimension pixmat(5001,5001)


      common /a4/pixmat

c

      kfund=150000


c
c     Calcula valor tipico de fundo de ceu, assumindo fundo de ceu plano, 
c     pelo ajuste de histograma de contagens
c
c     limite de contagens associadas a fundo de ceu:
c
c     limite = ceu + fatceu * sigceu
c
c


      ix1=0.15*nx
      ix2=0.85*nx
      iy1=0.15*ny
      iy2=0.85*ny



c
c     Cuidados com contagens negativas validas na indexacao do fundo de
c     ceu 
c

      xfu0=1.d0+malisa-vmin


      nfund=xfu0+mazi

      if (nfund.gt.kfund) nfund=kfund

c
c     Montando histograma de fundo de ceu
c

      do i=1,nfund
      ifundo(i)=0
      jfundo(i)=0
      enddo      

c

      do    i=iy1,iy2
      do 10 j=ix1,ix2
      if (pixmat(j,i).gt.mazi) go to 10
      n=pixmat(j,i)+xfu0
      jfundo(n)=jfundo(n)+1
  10  continue
      enddo



c
c     Alisa valores do histograma de fundo de ceu para evitar bins
c     com contagem baixa (ou mesmo zero) contiguos ou proximos ao
c     bin do pico, evitando obter valores irrealisticamente baixos
c     para o sigma de fundo de ceu  
c
c     1/2 do Filtro de media passante com N bins: N bins a esquerda,
c     bin central + N bins a direita -> valor alisado para (2N+1) bins,
c     ie, filtro com banda passante efetiva de (2N+1) canais
c
c     N = malisa
c

c
c     Miolo do histograma
c

      do j=malisa+1,nfund-malisa-1
      ifundo(j)=jfundo(j)
      do i=1,malisa
      ifundo(j)=ifundo(j)+jfundo(j-i)+jfundo(j+i)
      enddo
      ifundo(j)=ifundo(j)/(2.d0*malisa+1.d0)
      enddo

c
c     Pontas do histograma
c


      do i=1,malisa*2+1
      ifundo(1)=ifundo(1)+jfundo(i)
      enddo
      ifundo(1)=ifundo(1)/(malisa*2.d0+1.d0)
c     do i=2,malisa*2+1
      do i=2,malisa
      ifundo(i)=ifundo(1)
      enddo


      do i=nfund,nfund-malisa*2,-1
      ifundo(nfund)=ifundo(nfund)+jfundo(i)
      enddo
      ifundo(nfund)=ifundo(nfund)/(malisa*2+1)
c     do i=nfund-1,nfund-malisa*2,-1
      do i=nfund-1,nfund-malisa+1,-1
      ifundo(i)=ifundo(nfund)
      enddo


c
c     Pega valor tipico de fundo e sigma
c

      izmax=-1


      do n=nfund,1,-1
      if (ifundo(n).gt.izmax) then
      izmax=ifundo(n)
      nn=n
      endif
      enddo


c
c     Retoma valor de fundo de ceu (e threshold) coerentes com a contagem
c     de fundo de ceu
c

      ceu=nn
      ceu=ceu-xfu0


      cort=izmax/2.d0

c

      do n=nfund,nn+1,-1
      if (ifundo(n).gt.cort) go to 20 
      enddo
  20  sigceu=2.d0*(n-nn)


c
 
      threso=ceu+fatceu*sigceu

c

      return
      end





c
c     Subrotina alfa  
c
c
c     Utiliza calculo de (RA) apartir da projecao gnomonica inversa das coordenadas
c     padrao (X,Y).
c
c
c     zz    -  alfa central (rad)
c     ww    -  delta central (rad)
c     
c     xx    - alfa  em coordenadas padrao cartesianas no plano tangente 
c     yy    - delta em coordenadas padrao cartesianas no plano tangente 
c
c     alfa  - alfa  do ponto (rad) 
c     delta - delta do ponto (rad) 
c
c
c     Mod.  M. Assafin  26/Out/2012
c
c


      subroutine alfa (xx,yy,zz,ww,alf)


      implicit real*8 (a-h,o-z)



c
C     Daddos iniciais
C
c
      pi    = 0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi



c
c     Calculo de alfa
c 

  

      yc=xx
      xc=dcos(ww)-yy*dsin(ww)

      y=dabs(yc)
      x=dabs(xc)


      ang=dabs(datan2(y,x))


      if (yy.ge.0.d0) then

      if (xx.lt.0.d0) ang=pi-ang


      else

      if (xx.ge.0.d0) then

      ang=2.d0*pi-ang
      else
      ang=pi+ang
      endif

      endif

      alf=zz+ang

      if (alf.gt.2.d0*pi) alf=alf-2.d0*pi
      if (alf.lt.0.d0) alf=alf+2.d0*pi


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




