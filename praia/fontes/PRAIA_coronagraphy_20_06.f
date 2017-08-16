C
C     Programa PRAIA_coronagraphy_20
C
C     PROPOSITO
C
C       Elimina a imagem de planeta (ou objeto brilhante) numa imagem CCD.
C     Uma nova tecnica e' usada no lugar da tecnica de rebatimento. A hipotese
C     e' que o brilho do planeta tenha dependencia radial apenas.
C
C     A contagem do planeta e' contabilizada pela contagem ao longo do
C     perimetro de aneis. A cada novo anel de raio R, nova contagem e'
C     estabelecida. O valor representativo da contagem do planeta a cada
C     novo anel e' dado pelo valor medio obtido.
C
C     Para descartar pixels com contagens mais altas, associados a satelites,
C     estrelas intrusas e cruzes de difracao, as contagens mais altas sao
C     eliminadas, mantendo-se apenas os valores mais tipicos no calculo.
C     Qdo ha poucos pixels, tomamos a media M e o desvio-padrao S,
C     depois eliminamos os valores acima de:
C
C
C       v > M + 2,0 S
C
c     Havendo um numero minimo de pixels, montamos um histograma de contagens
c     para determinar o valor tipico.
c
C     Finalmente, calculamos a media definitiva com os valores retidos.
C
C     Os aneis tem a expessura de um quadrado de lados 1x1 pixels. Isto e',
C     os valores de contagens considerados sao aqueles do proprio pixel
C     interceptado pelo anel de raio R. Os valores das contagens sao
C     pesados de acordo com a diferenca D entre o raio R do anel e a
C     distancia do centro de cada pixel ao centro (x,y) do planeta:
C
C         peso=1/(D+1)
C
C     Subtraindo os valores de contagem reais do planeta calculados pixel a
C     pixel, obtem-se a imagem livre da luz do planeta. Pegando a imagem
C     original e subtraindo desta, obtem-se de volta uma imagem mais limpa
C     do proprio planeta. O processo e' repetido para ir aprimorando o
C     (x,y) do planeta em relacao ao qual os aneis sao construidos. 
C     
C
C
C     O (x,y) inicial vem de ajustes Gaussianos. A medida em que avancam as 
c     iteracoes, o (x,y) do centro do planeta vai sendo refinado, ate que
c     os residuos de luz, dentro de um certo raio, nao mudem dentro de
c     um parametro de convergencia (digamos, 1 por cento das contagens
c     residuais).
c
c
c     Variaveis internas importantes:
c
C     imagem(j,i): imagem original preservada
C     pixmat(j,i): imagem original na 1a.iteracao do llop gaussiano;
C                  imagem limpa do objeto brilhante no final
C     matpix(j,i): imagem do objeto brilhante submetida a ajuste gaussiano;
C                  imagem rebatida (livre de planeta) apos aplicacao dos aneis
C
C
C
C     Ultima modificacao: M. Assafin - 25/Jan/2011
C
C     Le imagens integer ou real em f77
C
c



      IMPLICIT REAL *8 (A-H,O-Z)

      integer*2 bitpix,bitpyx,betpix

      real*4 MATPIX(5001,5001),PIXMAT(5001,5001),imagem(5001,5001),
     ?maximo,mazi,malha(5001,5001)

      dimension xesto(25010001),yesto(25010001),po(25010001),cofo(136)

c     dimension contag(50000),peso(50000),histo(30),ico(30),pco(30),
c    ?ior(30),nval(30)

      character*50  iframes,ifiles(50000),badb,badc
      CHARACTER*50  infits
      CHARACTER*58  ipfits,irfits
      CHARACTER*5   mfits
      CHARACTER*3   ixy
      CHARACTER*112 manda
      CHARACTER*1   iplanr,iplaur,iplanp,iplaup
      character*1   iplan,iplau,menos,iver,ipl
      character*69  ichobj

      character*2   iextp,iextc



      COMMON/A3/IMAGEM
      COMMON/A4/PIXMAT
      COMMON/A5/MATPIX


      data menos/'-'/

c
C     Daddos iniciais
C

      ipmax=5001


      kfiles=50000

      idimo=25010001
      idimog=136

c
c     Abre arquivo de dados
c      

c     open (91,file='PRAIA_coronagraphy_20_06.dat')

      read (*,3) iframes
 3    format(a50)

c

      read (*,*) inicio,iultmo

      if (inicio.eq.0 .and. iultmo.eq.0) then

      open (1,file=iframes)

      i=0
 4    read (1,*,end=5)
      i=i+1
      go to 4

 5    close (1)

      inicio=1
      iultmo=i

      endif

c
c     Pixel mask of bad rectangle and circle regions
c

      read (*,3) badb
      read (*,3) badc

c
c     Le extensao de fits de saida 
c

      read (*,6) iextp
      read (*,6) iextc
 6    format(a2)

      read (*,*) maximo
      read (*,*) vneg
      read (*,*) ipflag
      read (*,*) bscale
      read (*,*) bzero
      read (*,*) bitpyx
      read (*,*) kswap


      read (*,*) ngrauf
      read (*,*) malisa
      read (*,*) fatceu

      read (*,*) ibusca
      read (*,*) ixb1,ixb2,iyb1,iyb2


      read (*,*) rcoro
      read (*,*) keycxy
      read (*,*) keycor
      read (*,*) ladoa
      read (*,*) convcp
      read (*,*) fatcor
      read (*,*) ixr1,ixr2,iyr1,iyr2
      read (*,*) ixt1,ixt2,iyt1,iyt2


c
c     Output FITS image standard parameters
c


      if=83
      scale=1.d0
      zero=0.d0

      read (*,*) betpix
      read (*,*) mswap

      if (betpix.gt.0) betpix=-betpix

c

      vmin=vneg

      rcoro=rcoro/100.d0

      ladox=ladoa
      ladoy=ladoa

      mazi=maximo

c
c
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'Digital coronagraphy for bright objects')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
c


      open (96,file='PRAIA_coronagraphy_20.rel')


c
c     Inicia leitura do arquivo de frames
c

      open (10,file=iframes)

      do i=1,kfiles
      read (10,3,end=7) ifiles(i)
      enddo

 7    close (10)


c
c
      
      do 60 kkkkkk=inicio,iultmo

      infits=ifiles(kkkkkk)

c
c     Prepara nomes e imagens
c

c
c     Nome de arquivo de planeta puro:         iextp              
c     Nome de arquivo de imagem coronagrafada: iextc            
c

      ipfits=infits
      irfits=infits

      do kj=50,1,-1
      if (infits(kj:kj).ne.' ') go to 12
      enddo

 12   continue

      if (infits(kj-4:kj).eq.'.fits'.or.infits(kj-4:kj).eq.'.FITS') then

      ipfits(kj-4:kj-4)='_'
      irfits(kj-4:kj-4)='_'

      ipfits(kj-3:kj-2)=iextp
      irfits(kj-3:kj-2)=iextc

      ipfits(kj-1:kj+3)='.fits'
      irfits(kj-1:kj+3)='.fits'

      else


      ipfits(kj+1:kj+1)='_'
      irfits(kj+1:kj+1)='_'

      ipfits(kj+2:kj+3)=iextp
      irfits(kj+2:kj+3)=iextc

      ipfits(kj+4:kj+8)='.fits'
      irfits(kj+4:kj+8)='.fits'

      endif


c

      write (manda,22) infits,ipfits
 22   format('cp ',a50,1x,a58)


      call system (manda)

      write (manda,22) infits,irfits

      call system (manda)

c

      write (*,*)
      write (*,26) infits
      write (96,26) infits
 26   format (23x,'Processing fits frame  ',a50)



c
c     Zerando imagem, pixmat e matpix
c

      do i=1,ipmax
      do j=1,ipmax
      matpix(j,i)=0.
      pixmat(j,i)=0.
      imagem(j,i)=0.
      malha(j,i)=0.
      enddo
      enddo

C
C     Lendo a imagem fits do campo, carrega matriz original, reconhece
C     tamanho em pixels em x e em y, e o valor maximo de contagem
c
c     nheads,iheads sao variaveis auxiliares para escrever o header
c     original nas imagens rebatidas e subtraidas
C


      bitpix=bitpyx

      call refits (infits,nx,ny,nheads,ichobj,ipflag,bscale,bzero,
     ?kswap,iswap,bitpix)


c
c     Marca pixels negativos (normalmente saturacao no CCD), minimo,
c     maximo (normalmente limite de linearidade)
c



      do i=1,ny
      do j=1,nx
      if (pixmat(j,i).lt.0.d0) pixmat(j,i)=maximo
      if (pixmat(j,i).le.vmin) pixmat(j,i)=vmin
      if (pixmat(j,i).gt.maximo) pixmat(j,i)=maximo
      enddo
      enddo



c
c     Bad pixels, circular regions
c

      open (23,file=badc,status='old',err=205)

 200  continue

      read (23,*,end=205) xmask,ymask,rmask 

      ix1=xmask-rmask+1
      ix2=xmask+rmask+1
      iy1=ymask-rmask+1
      iy2=ymask+rmask+1

      if (ix1.lt.1) ix1=1
      if (iy1.lt.1) iy1=1
      if (ix2.gt.nx) ix2=nx
      if (iy2.gt.ny) iy2=ny


      do i=iy1,iy2
      do j=ix1,ix2
      call circul (rmask,xmask,ymask,j,i,ichave)
      if (ichave.eq.1) pixmat(j,i)=maximo
      enddo
      enddo

      go to 200


 205  close (23)


c
c     Bad pixels, rectangular regions
c


      open (23,file=badb,status='old',err=220)

 210  continue

      read (23,*,end=220) x1,x2,y1,y2

      ix1=x1
      ix2=x2
      iy1=y1
      iy2=y2

      if (ix1.lt.1) ix1=1
      if (iy1.lt.1) iy1=1
      if (ix2.gt.nx) ix2=nx
      if (iy2.gt.ny) iy2=ny


      do i=iy1,iy2
      do j=ix1,ix2
      pixmat(j,i)=maximo
      enddo
      enddo

      go to 210


 220  close (23)




c
c     Preenche matriz imagem de backup
c


c     do i=1,ny
c     do j=1,nx
c     imagem(j,i)=pixmat(j,i)
c     enddo
c     enddo


c
c     Flattens sky background: 1rst step with star contamination
c



      if (ngrauf.ge.1 .and. ngrauf.le.15) then

      write (*,*)
      write (*,*) 'Flattening sky background. Please wait ...'
      write (*,*)

      fcmin=vmin
      fcmax=mazi

      call flatsk (nx,ny,ngrauf,mazi,fcmin,fcmax,xesto,yesto,po,cofo)

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

      call flatsk (nx,ny,ngrauf,mazi,fcmin,fcmax,xesto,yesto,po,cofo)

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
c     Computes resulting sky background for whole image
c


      fcmin=ceu-2.5d0*sigceu
      fcmax=ceu+2.5d0*sigceu

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
      do 28 i=1,ny
      do 27 j=1,nx
c     if (imagem(j,i).ge.maximo) go to 27
      if (pixmat(j,i).gt.fcmax) go to 27
      if (pixmat(j,i).lt.fcmin) go to 27
      n=n+1
      xesto(n)=j
      yesto(n)=i
      po(n)=pixmat(j,i)
c     xmean=xmean+po(n)
 27   continue
 28   continue

      call isolo (ngrauf,n,xesto,yesto,po,cofo)

c     xmean=xmean/n

c
c     Preenche matriz imagem de backup com fundo de ceu flatado
c


      do i=1,ny
      do j=1,nx
c     if (imagem(j,i).ge.maximo) go to 29
      imagem(j,i)=pixmat(j,i)
      enddo
      enddo

c
c
c     Reconhece o objeto mais brilhante no campo
c
c     Bina a matriz em celulas de tamanho ibusca x ibusca
c
c

      if(ixb1.eq.0 .and. ixb2.eq.0 .and. iyb1.eq.0 .and. iyb2.eq.0) then
      nnx=nx
      nny=ny
      ix0=1
      iy0=1
      else
      nnx=ixb2-ixb1+1
      nny=iyb2-iyb1+1
      ix0=ixb1
      iy0=iyb1
      endif

      if (nnx.gt.nx) nnx=nx
      if (nny.gt.ny) nny=ny
      if (ix0.lt.0) ix0=1
      if (iy0.lt.0) iy0=1

c

      soma=0.d0

      ibx=nnx/ibusca
      iby=nny/ibusca

      do 40 ii=1,iby
      iy1=(ii-1)*ibusca+iy0
      iy2=ii*ibusca+iy0

      if (iy1.lt.1)  iy1=1
      if (iy2.gt.ny) iy2=ny

      do 39 jj=1,ibx
      ix1=(jj-1)*ibusca+ix0
      ix2=jj*ibusca+ix0

      if (ix1.lt.1)  ix1=1
      if (ix2.gt.nx) ix2=nx
      

      conta=0.d0
      cx=0.d0
      cy=0.d0

      do 38 i=iy1,iy2
      do 37 j=ix1,ix2
      if (imagem(j,i).lt.maximo) then
      cx=cx+j*pixmat(j,i)
      cy=cy+i*pixmat(j,i)
      conta=conta+pixmat(j,i)
      endif
 37   continue
 38   continue

      if (conta.gt.soma) then
      bx=cx/conta
      by=cy/conta
      soma=conta
      endif

 39   continue

 40   continue

      ibx=bx
      iby=by

c
c    Raio da regiao de ajuste gaussiano do objeto brilhante e de coronagrafia
c    calculado a partir do threshold minimo de luz fornecido pelo usuario em
c    porcentagem de contagem maxima do pico central do objeto brilhante
c
c


      soma=pixmat(ibx,iby)*rcoro

      do j=1,nx
      ix1=ibx+j
      if (pixmat(ix1,iby).lt.soma) go to 41
      enddo

 41   lado=j





      ix1=ibx-lado
      ix2=ibx+lado
      iy1=iby-lado
      iy2=iby+lado

      if (ix1.lt.1) ix1=1
      if (iy1.lt.1) iy1=1
      if (ix2.gt.nx) ix2=nx
      if (iy2.gt.ny) iy2=ny


c
c     Ajuste gaussiano do planeta
c

      do i=1,ny
      do j=1,nx
      matpix(j,i)=pixmat(j,i)
      enddo
      enddo


c
c     Ajuste de Gaussiana Circular ao objeto brilhante, dentro da area
c     delimitada encontrada, para determinar com precisao o centro (bx,by)
c     do objeto brilhante. A regiao ajustada e' circular, de centro (bx,by)
c     e de raio r=lado. Ha' trimming circular no calculo gaussiano. Os
c     pixels saturados (contagem=maximo) nao participam do ajuste.
c


      r=lado

c     write (*,*) 'gausc (x,y) = ',bx,by
c     write (96,*)'gausc (x,y) = ',bx,by


c     mazi=maximo
c     call GCC (ix1,ix2,iy1,iy2,r,maximo,bx,by)
c     maximo=mazi



      write (*,*) 'gausc (x,y) = ',bx,by
      write (*,*) 'radius      = ',r     
      write (96,*)'gausc (x,y) = ',bx,by
      write (96,*)'radius      = ',r     


c
c     Deu pau?
c

      if (ixb1.ne.0.and.ixb2.ne.0.and.iyb1.ne.0.and.iyb2.ne.0) then
      
      if (bx.lt.ixb1) go to 60
      if (by.lt.iyb1) go to 60
      if (bx.gt.ixb2) go to 60
      if (by.gt.iyb2) go to 60

      else

      if (bx.lt.1) go to 60
      if (by.lt.1) go to 60
      if (bx.gt.nx) go to 60
      if (by.gt.ny) go to 60

      endif

c
c     Upgrade of bright object window from (x,y) Gauss fit
c


      ix1=bx-lado 
      ix2=bx+lado 
      iy1=by-lado 
      iy2=by+lado 


c
c     Determinacao final do centro do objeto brilhante pelos residuos
c     de contagens de luz
c


      call rede (maximo,ix1,ix2,iy1,iy2,bx,by,ladox,ladoy,
     ?malha,convcp,keycxy)


c
c     Deu pau?
c

      if (ixb1.ne.0.and.ixb2.ne.0.and.iyb1.ne.0.and.iyb2.ne.0) then
      
      if (bx.lt.ixb1) go to 60
      if (by.lt.iyb1) go to 60
      if (bx.gt.ixb2) go to 60
      if (by.gt.iyb2) go to 60

      else

      if (bx.lt.1) go to 60
      if (by.lt.1) go to 60
      if (bx.gt.nx) go to 60
      if (by.gt.ny) go to 60

      endif
      


c
c    Nesta versao, ao final, alem de uma pequena area em entorno do planeta
c    (definida pelo fator skyfac), sao corrigidas por coronagrafia digital
c    mais 2 areas predeterminadas pelo usuario (por exemplo: todo o CCD; ou
c    uma pequena area em torno do alvo mais uma segunda area em torno da
c    estrela de calibracao fotometrica), preservando-se o restante da imagem
c    original.
c


c
c     Coronagrafia em toda imagem
c

      if(ixr1.eq.0 .and. ixr2.eq.0 .and. iyr1.eq.0 .and. iyr2.eq.0 .and.
     ?ixt1.eq.0 .and. ixt2.eq.0 .and. iyt1.eq.0 .and. iyt2.eq.0) then

      imx1=1
      imx2=nx
      imy1=1
      imy2=ny

      if (keycor.eq.1) call fupla  (bx,by,imx1,imx2,imy1,imy2,mazi)
      if (keycor.eq.2) call fuplas (bx,by,imx1,imx2,imy1,imy2,mazi)

      go to 45
      endif



c
c     Coronagrafia em uma parte da imagem em torno do objeto brilhante,
c     com tamanho da regiao automaticamente calculado, ateh que as contagens
c     do anel da coronografia estejam dentro da faixa de contagens do
c     fundo de ceu da imagem como um todo
c

      xnn=nx
      ynn=ny
      ndi=dsqrt(xnn**2+ynn**2)+1.d0


      if(ixr1.lt.0 .and. ixr2.lt.0 .and. iyr1.lt.0 .and. iyr2.lt.0 .and.
     ?ixt1.lt.0 .and. ixt2.lt.0 .and. iyt1.lt.0 .and. iyt2.lt.0) then

      n=0

 42   n=n+1
      imx1=ix1-n
      imx2=ix2+n
      imy1=iy1-n
      imy2=iy2+n
      if (imx1.lt.1) imx1=1
      if (imy1.lt.1) imy1=1
      if (imx2.gt.nx) imx2=nx
      if (imy2.gt.ny) imy2=ny

      call skyp (imx1,imx2,imy1,imy2,mazi,malisa,vmin,fatceu,pceu,
     ?psig,ipzmax,pthres)

      auxxx=ceu+fatcor*sigceu

      if (pceu.lt.auxxx) go to 43


c     if (pceu.lt.fcmax) go to 43

c     if (dabs(psig-sigceu).lt.sigceu*0.5d0) go to 43


      if (n.gt.ndi) go to 43

      go to 42


 43   continue



      if (keycor.eq.1) call fupla  (bx,by,imx1,imx2,imy1,imy2,mazi)
      if (keycor.eq.2) call fuplas (bx,by,imx1,imx2,imy1,imy2,mazi)

      go to 45

      endif



c
c     Coronagrafia em torno da primeira area (ex: estrela de calibracao)
c


      if (keycor.eq.1) call fupla  (bx,by,ixr1,ixr2,iyr1,iyr2,mazi)
      if (keycor.eq.2) call fuplas (bx,by,ixr1,ixr2,iyr1,iyr2,mazi)


c
c     Coronagrafia em torno da segunda area (ex: alvo da fotometria)
c

      if (keycor.eq.1) call fupla  (bx,by,ixt1,ixt2,iyt1,iyt2,mazi)
      if (keycor.eq.2) call fuplas (bx,by,ixt1,ixt2,iyt1,iyt2,mazi)


c
c     Ending proccess
c

 45   continue



c
c
c     Writes output images 
c     
c     matpix = coronographed image
c     pixmat = bright object profile
c     imagem = original image 
c
c
c


c
c     Marks bad pixels 
c

c     do i=1,ny
c     do j=1,nx
 
c     if (imagem(j,i).ge.maximo) matpix(j,i)=0. 
 
c     enddo
c     enddo


c
c     Bright object profile
c


      do i=1,ny
      do j=1,nx
      pixmat(j,i)=imagem(j,i)
      enddo
      enddo


      if(ixr1.eq.0 .and. ixr2.eq.0 .and. iyr1.eq.0 .and. iyr2.eq.0 .and.
     ?ixt1.eq.0 .and. ixt2.eq.0 .and. iyt1.eq.0 .and. iyt2.eq.0) then

      do i=1,ny
      do j=1,nx
 
      pixmat(j,i)=imagem(j,i)-matpix(j,i)
      if (pixmat(j,i).lt.1.) pixmat(j,i)=imagem(j,i)
 
      enddo
      enddo


      call wfitsp (if,ipfits,betpix,mswap,nx,ny,scale,zero,nheads)

      endif


      if(ixr1.lt.0 .and. ixr2.lt.0 .and. iyr1.lt.0 .and. iyr2.lt.0 .and.
     ?ixt1.lt.0 .and. ixt2.lt.0 .and. iyt1.lt.0 .and. iyt2.lt.0) then


      do i=imy1,imy2
      do j=imx1,imx2
      xx=j
      yy=i
      pixmat(j,i)=imagem(j,i)-matpix(j,i)
      if (pixmat(j,i).lt.1.) pixmat(j,i)=imagem(j,i)
      enddo
      enddo


      call wfitsp (if,ipfits,betpix,mswap,nx,ny,scale,zero,nheads)

      endif



      if(ixr1.gt.0 .and. ixr2.gt.0 .and. iyr1.gt.0 .and. iyr2.gt.0) then


      do i=iyr1,iyr2
      do j=ixr1,ixr2
      xx=j
      yy=i
      pixmat(j,i)=imagem(j,i)-matpix(j,i)
      if (pixmat(j,i).lt.1.) pixmat(j,i)=imagem(j,i)
      enddo
      enddo

      call wfitsp (if,ipfits,betpix,mswap,nx,ny,scale,zero,nheads)

      endif


      if(ixt1.gt.0 .and. ixt2.gt.0 .and. iyt1.gt.0 .and. iyt2.gt.0) then


      do i=iyt1,iyt2
      do j=ixt1,ixt2
      xx=j
      yy=i
      pixmat(j,i)=imagem(j,i)-matpix(j,i)
      if (pixmat(j,i).lt.1.) pixmat(j,i)=imagem(j,i)
      enddo
      enddo


      call wfitsp (if,ipfits,betpix,mswap,nx,ny,scale,zero,nheads)

      endif



c
c     Coronographed image
c

      do i=1,ny
      do j=1,nx

      pixmat(j,i)=matpix(j,i)

      enddo
      enddo




c
c     Restores sky background counts under coronographed areas
c


      if(ixr1.eq.0 .and. ixr2.eq.0 .and. iyr1.eq.0 .and. iyr2.eq.0 .and.
     ?ixt1.eq.0 .and. ixt2.eq.0 .and. iyt1.eq.0 .and. iyt2.eq.0) then

      do i=1,ny
      do j=1,nx
      xx=j
      yy=i
      pixmat(j,i)=pixmat(j,i)+polo(xx,yy,cofo,ngrauf)
      enddo
      enddo

      go to 50

      endif

c
      if(ixr1.lt.0 .and. ixr2.lt.0 .and. iyr1.lt.0 .and. iyr2.lt.0 .and.
     ?ixt1.lt.0 .and. ixt2.lt.0 .and. iyt1.lt.0 .and. iyt2.lt.0) then


      do i=imy1,imy2
      do j=imx1,imx2
      xx=j
      yy=i
      pixmat(j,i)=pixmat(j,i)+polo(xx,yy,cofo,ngrauf)
      enddo
      enddo

      go to 50

      endif

c

      do i=iyr1,iyr2
      do j=ixr1,ixr2
      xx=j
      yy=i
      pixmat(j,i)=pixmat(j,i)+polo(xx,yy,cofo,ngrauf)
      enddo
      enddo


      do i=iyt1,iyt2
      do j=ixt1,ixt2
      xx=j
      yy=i
      pixmat(j,i)=pixmat(j,i)+polo(xx,yy,cofo,ngrauf)
      enddo
      enddo


c

 50   continue

c
      call wfitsp (if,irfits,betpix,mswap,nx,ny,scale,zero,nheads)


c
c 

 60   continue

c



      close (96)


c     write (*,22)
      write (*,*)
      write (*,61) 
 61   format (23x,'Coronography terminated successfuly.')
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '

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



      subroutine  GCC (nx,mx,ny,my,raio,maximo,bx,by)

      IMPLICIT REAL*8 (A-H,O-Z)
      real*4 matpix(5001,5001),maximo,xmaxi
      DIMENSION DISMAX(5001),DISMAY(5001),PARAMX(4),PARAMY(4),DELTAX(5),
     ?XSIGMA(5),PARAM(5)
 
      COMMON /A5/MATPIX
C
C     INICIALIZACAO E ENTRADA DE DADOS
C

      xmaxi=maximo


c     ESCL   = 0.2925D0
      escl=1.d0
      X0 = 0.D0
      Y0 = 0.D0
      PASSOX = 1.D0
      PASSOY = 1.D0
      CENTRX = 0.D0
      CENTRY = 0.D0
C
C     CONVERGENCIA EM DENSIDADE = 1% ; EM POSICAO = 0,001 pixel
C
      DLIMIT=1.D-2
      PLIMIT=1.D-3
C
      ICONTT=0
      ICONTX=1
      TRANSX=0.D0
      TRANSY=0.D0
      SIGX=0.D0
      SIGY=0.D0
C
C     CONSTRUINDO AS DISTRIBUICOS  MARGINAIS DE X E Y
C
      DO 80 I=NX,MX
   80 DISMAX(I)=0.D0
      DO 801 I=NY,MY
  801 DISMAY(I)=0.D0
C
C     TRANSFORMANDO A CONVERGENCIA DE (") EM FRACAO DE PIXELS.
C
      PLIMIT=PLIMIT/ESCL
C
C     CONSTRUINDO AS DISTRIBUICOES MARGINAIS DE X E Y.
C
      DO 83 I=NY,MY
      DO 82 K=NX,MX
      DISMAX(K)= DISMAX(K) + MATPIX(K,I)
      DISMAY(I)= DISMAY(I) + MATPIX(K,I)
   82 CONTINUE
   83 CONTINUE
C
c     RAIO=1000.D0
      XLAMDA=0.001
      NTERMX=5
      XRESID= -10.D10
      RESIDX=0.D0
      XCENT=1.D14
      YCENT=1.D14
      DEMINX = 10.D10
      DEMINY = 10.D10
      PARAMX(1) = -10.D10
      PARAMY(1) = -10.D10
C
C     CALCULANDO PARAMETROS INICIAIS PARA MARGINAL X
C     (AUER E VAN ALTENA II)
C
      DO 84 I=NX,MX
      IF (DISMAX(I).LT.DEMINX) DEMINX = DISMAX(I)
      IF (DISMAX(I).LT.PARAMX(1)) GO TO 84
      PARAMX(1) = DISMAX(I)
      PARAMX(2) = I
   84 CONTINUE
      DO 142 I=NX,MX
  142 DISMAX(I) = DISMAX(I) - DEMINX
C
      PARAMX(1) = PARAMX(1) - DEMINX
      AUX = PARAMX(1)/2.D0
      II = PARAMX(2) - 1
      DO 86 I=NX,II
   86 IF (DISMAX(I).GT.AUX) GO TO 87
C
   87 SIGX=I-1.D0
      II = PARAMX(2) + 1
      DO 88 I=MX,II,-1
   88 IF (DISMAX(I).GT.AUX) GO TO 89
C
   89 SIGY=I-1.D0
      PARAMX(3)=(SIGY-SIGX)/2
C
C     CALCULANDO PARAMETROS INICIAIS PARA MARGINAL Y
C     (AUER E VAN ALTENA II)
C
      DO 85 I=NY,MY
      IF (DISMAY(I).LT.DEMINY) DEMINY = DISMAY(I)
      IF (DISMAY(I).LT.PARAMY(1)) GO TO 85
      PARAMY(1) = DISMAY(I)
      PARAMY(2) = I
   85 CONTINUE
      DO 143 I=NY,MY
  143 DISMAY(I) = DISMAY(I) - DEMINY
C
      PARAMY(1) = PARAMY(1) - DEMINY
      AUX = PARAMY(1)/2.D0
      II = PARAMY(2) - 1
      DO 910 I=NY,II
  910 IF (DISMAY(I).GT.AUX) GO TO 911
C
  911 SIGX=I-1.D0
      II = PARAMY(2) + 1
      DO 912 I=MY,II,-1
  912 IF (DISMAY(I).GT.AUX) GO TO 913
C
  913 SIGY=I-1.D0
      PARAMY(3)=(SIGY-SIGX)/2
C
C     PARAMETROS DA GAUSSIANA BIDIMENSIONAL SIMETRICA
C
      IAAUX=PARAMX(2)
      IAAUY=PARAMY(2)
      PARAM(1)=PARAMX(2)
      PARAM(2)=PARAMY(2)
      PARAM(3)=(PARAMX(3)+PARAMY(3))/2.D0
      PARAM(4)=MATPIX(IAAUX,IAAUY)
      PARAM(5)=1.D0
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
C       MATPIX - MATRIZ DE PIXELS DA IMAGEM
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
C       QIQUAD (KEY,RAIO,MATPIX,NX,NY,MX,MY,FREE,PARAM)
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
      real*4 matpix(5001,5001),maximo

      COMMON /A5/MATPIX 
      COMMON /A6/ALPHA,BETA
C     COMMON /A13/DERIV
      COMMON /A7/ARRAY
      COMMON /A14/IERRO
      IERRO=0
      DET=1.D0
      ICONT=0
      IXC=A(1)
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
      IF (MATPIX(L,I).ge.maximo) GO TO 5001
      CALL GDERIV (L, I, A, DELTAA, DERIV, NTERMS)
      DO 146 J=1,NTERMS
      BETA(J)= BETA(J) + (MATPIX(L,I)-FGAUSI(L,I,A))*DERIV(J)
      DO 146 K=1, J
  146 ALPHA(J,K) = ALPHA(J,K) + DERIV(J)*DERIV(K)
 5001 CONTINUE
  150 CONTINUE
C
      IF (IERRO.EQ.1) THEN
      IERRO=0
      GO TO 107
      ENDIF
C
      FREE = (MX-NX+1)+(MY-NY+1) - NTERMS
C
      IF (FREE.LE.0.D0) THEN
      IERRO=0
      GO TO 107
      ENDIF
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
      IF (MATPIX(L,I).ge.maximo) GO TO 5000
      CALL CIRCUL (RAIO,A(1),A(2),L,I,ICHAVE)
      IF (ICHAVE.LT.0) GO TO 5000
      CALL GDERIV (L, I, A, DELTAA, DERIV, NTERMS)
      ICONT=ICONT+1
      DO 46 J=1,NTERMS
      BETA(J)= BETA(J) + (MATPIX(L,I)-FGAUSI(L,I,A))*DERIV(J)
      DO 46 K=1, J
   46 ALPHA(J,K) = ALPHA(J,K) + DERIV(J)*DERIV(K)
 5000 CONTINUE
   50 CONTINUE
C
      IF (IERRO.EQ.1) THEN
      IERRO=0
      GO TO 107
      ENDIF
C
C
C     CALCULA GRAUS DE LIBERDADE FREE
C
      COUNT=ICONT
      FREE=2*DSQRT(COUNT)-NTERMS
C
      IF (FREE.LE.0.D0) THEN
      IERRO=0
      GO TO 107
      ENDIF
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
      IF (IERRO.EQ.1) THEN
      IERRO=0
      GO TO 107
      ENDIF
C
C				 
C        INVERT MODIFIED CURVATURE MATRIX TO FIND NEW PARAMETERS
C
   71 DO 74 J=1, NTERMS
      DO 73 K=1, NTERMS
      AUX = ALPHA(J,J)*ALPHA(K,K)
      IF (AUX.LE.0.D0) GO TO 107
   73 ARRAY(J,K)= ALPHA(J,K) / DSQRT (AUX)
   74 ARRAY(J,J) = 1.D0 + FLAMDA
   80 CALL MATINV (NTERMS, DET)
C
      IF (IERRO.EQ.1) THEN
      IERRO=0
      GO TO 107
      ENDIF
C
   81 DO 84 J=1, NTERMS
      B(J) = A(J)
      DO 84 K=1, NTERMS
      AUX = ALPHA(J,J)*ALPHA(K,K)
      IF (AUX.LE.0.D0) GO TO 107
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
      IF (IERRO.EQ.1) THEN
      IERRO=0
      GO TO 107
      ENDIF
C
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
      IF (AUX.LE.0.D0) GO TO 107
  103 SIGMAA(J) = DSQRT (AUX)
  104 CONTINUE
      FLAMDA = FLAMDA/10.D0
      GO TO 110
  107 CHISQR = -1.D0
      IERRO=0
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
C       CALL GDERIV (J, I, A, DELTAA, DERIV, NTERMS)
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
      SUBROUTINE GDERIV (J, I, A, DELTAA, DERIV, NTERMS)
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
      IF ((ZX2.GT.50.D0).OR.(ZY2.GT.50.D0)) THEN
      DERIV(1) = 0.D0
      DERIV(2) = 0.D0
      DERIV(3) = 0.D0
      DERIV(4) = 0.D0
      ELSE
      FUNCAO   = FGAUSI(J,I,A)
      IF (IERRO.EQ.1) RETURN
      DERIV(1) = (X-A(1))*(FUNCAO-A(5))/A(3)**2
      DERIV(2) = (Y-A(2))*(FUNCAO-A(5))/A(3)**2
      DERIV(3) = (ZX2+ZY2)*(FUNCAO-A(5))/A(3)
      DERIV(4) = (FUNCAO-A(5))/A(4)
      ENDIF
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
      IF (ZX2 - 50.D0) 16, 20, 20
   16 IF (ZY2 - 50.D0) 17, 20, 20
   17 FGAUSI = FGAUSI + A(4)*DEXP(-ZX2/2.D0 - ZY2/2.D0)
   20 RETURN
      END
C
C
C       FUNCTION QIQUAD
C
C       PURPOSE
C         EVALUATE REDUCED CHI SQUARE FOR FIT TO DATA
C            QIQUAD = SUM ((MATPIX-FITMAT)**2 / SIGMA**2) / NFREE
C
C       USAGE
C         RESULT = QIQUAD (KEY,RAIO,NX,NY,MX,MY,FREE,PARAM)
C
C       DESCRIPTION OF PARAMETERS
C         KEY    - 0         - SEM TRIMMING
C                  1 OU MAIS - TRIMMING CIRCULAR
C         MATPIX - MATRIX ARRAY OF DATA POINTS
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
      real*4 matpix(5001,5001),maximo

      DIMENSION A(5)
      COMMON /A5/MATPIX
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
      IF (MATPIX(J,I).ge.maximo) GO TO 38
      CHISQ = CHISQ + (MATPIX(J,I)-FGAUSI(J,I,A))**2
 38   continue
      GO TO 39
C
   30 IXC=A(1)
      DO 35 I = NY,MY
      CALL CIRCUL (RAIO,A(1),A(2),IXC,I,ICHAVE)
      IF (ICHAVE.LT.0) GO TO 35
      DO 350 J = NX,MX
      IF (MATPIX(J,I).ge.maximo) GO TO 350
      CALL CIRCUL (RAIO,A(1),A(2),J,I,ICHAVE)
      IF (ICHAVE.LT.0) GO TO 350
      CHISQ = CHISQ + (MATPIX(J,I)-FGAUSI(J,I,A))**2
  350 CONTINUE
   35 CONTINUE
C
C        DIVIDE BY NUMBER OF DEGREES OF FREEDOM
C
      IF (IERRO.EQ.1) RETURN
C
   39 QIQUAD = CHISQ / FREE
   40 CONTINUE
      RETURN
      END
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
c
c     Subrotina interpol
c
c
c     Interpola valor de contagem a ser associada a coordenada real (nao
c     inteira) de (xp,yp), dada uma matriz imagem pixmat
c
c     O calculo e' feito a partir da interpolacao pela distancia de 4 valores
c     de contagem: do proprio pixel onde (xp.yp) "cai", e de mais 3 pixels
c     vizinhos, do lado x, do lado y e da ponta mais proxima a (xp,yp) em
c     x e y.
c
c
c     Atencao: pixel "1", coordenada 0.5, pixel 2, coordenada 1.5, etc.
c
c

      subroutine interpol (nx,ny,xp,yp,cont,mazi)

      implicit real *8 (A-H,O-Z)

      real*4 imagem(5001,5001),pixmat(5001,5001)

      common/A3/imagem
      common/A4/pixmat
 


c
c     Interpola a contagem pelos pixels vizinhos, incluindo o que
c     contem o ponto (xp,yp), pesando pela distancia
c

      cont=0.d0
      dd=0.d0
      ix=xp+0.5d0
      iy=yp+0.5d0

      if (ix.lt.1)  ix=1
      if (ix.gt.nx) ix=nx
      if (iy.lt.1)  iy=1
      if (iy.gt.ny) iy=ny

c
c     Qual quadrante de (xp,yp) dentro do pixel contendo (xp,yp)
c
      kx=xp
      ky=yp
      dx=xp-kx
      dy=yp-ky
      if (dx.ge.0.5d0) then
      qx=1.d0
      else
      qx=-1.d0 
      endif
      if (dy.ge.0.5d0) then
      qy=1.d0
      else
      qy=-1.d0 
      endif

c
c     Contribuicao do proprio pixel contendo o ponto (xp,yp)
c
      dq=(dx-0.5d0)**2+(dy-0.5d0)**2
      if (dq.lt.1.d-6) then
      cont=pixmat(ix,iy)
c     cont=imagem(ix,iy)
      return
      endif

      d1=1.d0/dsqrt(dq)
      c1=pixmat(ix,iy)
c     c1=imagem(ix,iy)
      dd=dd+d1
      cont=cont+c1*d1

c
c     Contribuicao dos demais 3 pixels adjacentes
c

      kx=ix+1.d0*qx
      if (kx.lt.1)  go to 3780
      if (kx.gt.nx) go to 3780
      ky=iy
      if (ky.lt.1)  go to 3780
      if (ky.gt.ny) go to 3780

      if (imagem(kx,ky).ge.mazi) go to 3780

      d2=1.d0/dsqrt((xp-kx+0.5d0)**2+(yp-ky+0.5d0)**2)
      c2=pixmat(kx,ky)
c     c2=imagem(kx,ky)
      dd=dd+d2
      cont=cont+c2*d2

 3780 continue
c
      ky=iy+1.d0*qy
      if (ky.lt.1)  go to 3785
      if (ky.gt.ny) go to 3785
      kx=ix
      if (kx.lt.1)  go to 3785
      if (kx.gt.nx) go to 3785

      if (imagem(kx,ky).ge.mazi) go to 3785

      d3=1.d0/dsqrt((xp-kx+0.5d0)**2+(yp-ky+0.5d0)**2)
      c3=pixmat(kx,ky)
c     c3=imagem(kx,ky)
      dd=dd+d3
      cont=cont+c3*d3

 3785 continue
c

      kx=ix+1.d0*qx
      if (kx.lt.1)  go to 3790
      if (kx.gt.nx) go to 3790
      ky=iy+1.d0*qy
      if (ky.lt.1)  go to 3790
      if (ky.gt.ny) go to 3790

      if (imagem(kx,ky).ge.mazi) go to 3790

      d4=1.d0/dsqrt((xp-kx+0.5d0)**2+(yp-ky+0.5d0)**2)
      c4=pixmat(kx,ky)
c     c4=imagem(kx,ky)
      dd=dd+d4
      cont=cont+c4*d4

 3790 continue
c

      cont=cont/dd


      return
      end


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
      DIMENSION IORDEM(IDIM),NVAL(IDIM)
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
c     Subrotina fupla
c
c     Calcula o valor de contagem da luz do planeta para cada isofota r
c
c     Versao robusta. A versao mais rapida eh fuplas.
c
c     Last update:  M. Assafin  23/Jul/2010
c
c


      subroutine fupla (bx,by,ix1,ix2,iy1,iy2,mazi)

      IMPLICIT REAL *8 (A-H,O-Z)
      real*4 MATPIX(5001,5001),PIXMAT(5001,5001),imagem(5001,5001),
     ?maximo,mazi
      
      dimension contag(50000),peso(50000),histo(30),ico(30),pco(30),
     ?ior(30),nval(30)

      COMMON/A3/IMAGEM
      COMMON/A4/PIXMAT
      COMMON/A5/MATPIX


c
c     Dados iniciais
c

      nhist=30
      mpes=50000

c
c     Zera vetores
c

      do i=1,mpes
      contag(i)=0.d0
      peso(i)=0.d0
      enddo


      tol=dsqrt(2.d0)/2.d0

      do 499 i=iy1,iy2
      do 498 j=ix1,ix2

c
c     Para cada pixel, ha' um unico anel com um unico raio R. Portanto, para
c     cada pixel, calcula-se o anel correspondente, isto e', o valor da
c     contagem media real do planeta correspondente a esse R.
c
c

      xp=j
      yp=i

      raiop=dsqrt((xp-bx)**2+(yp-by)**2) 

      ixmin=bx-raiop-tol+1.d0
      iymin=by-raiop-tol+1.d0
      ixmax=bx+raiop+tol+1.d0
      iymax=by+raiop+tol+1.d0


      if (ixmin.lt.ix1) ixmin=ix1
      if (iymin.lt.iy1) iymin=iy1
      if (ixmax.gt.ix2) ixmax=ix2
      if (iymax.gt.iy2) iymax=iy2

      n=0


      do 1101 ii=iymin,iymax
      do 1100 jj=ixmin,ixmax

      if (imagem(jj,ii).ge.mazi) go to 1100

      x=jj
      y=ii

      raio=dsqrt((x-bx)**2+(y-by)**2)

      dd=dabs(raio-raiop)

      if (dd.gt.tol) go to 1100

      n=n+1
      contag(n)=pixmat(jj,ii)
      peso(n)=1.d0/(dd+1.d0)**2


 1100 continue
 1101 continue


c
c     Estatistica para eliminar contagens altas associadas a satelites,
C     estrelas intrusas e cruzes de difracao,
c


c     if (n.lt.1) then
c     matpix(j,i)=0.
c     go to 498
c     endif

      if (n.lt.1) go to 498




c
c     Calcula histograma, qdo ha bastante pixel
c

 1500 continue


      zmin=1.d14
      zmax=-1.d14


      do nn=1,n
      z=peso(nn)*contag(nn)
      if (z.gt.zmax) zmax=z
      if (z.lt.zmin) zmin=z
      enddo

      razao=zmax-zmin

      if (razao.lt.0.1d0) then
      cmed=zmax
      go to 1515
      endif


      do nn=1,nhist
      ico(nn)=0
      pco(nn)=0
      histo(nn)=0.d0
      enddo

      do nn=1,n
      z=peso(nn)*contag(nn)
      ko=1+(nhist-1)*(z-zmin)/razao
      ico(ko)=ico(ko)+1
      pco(ko)=pco(ko)+peso(nn)
      histo(ko)=histo(ko)+z
      enddo

      do 1501 nn=1,nhist
      if (ico(nn).eq.0) go to 1501
      histo(nn)=histo(nn)/pco(nn)
 1501 continue

c
c     Refina bins do histograma
c

      ncort=0.8d0*n

      mcort=0
      do nn=1,nhist
      mcort=mcort+ico(nn)
      if (mcort.gt.ncort) go to 1505
      enddo

 1505 cort=histo(nn)

      razao=cort-zmin

      if (razao.lt.0.1d0) then
      cmed=cort
      go to 1515
      endif


      do nn=1,nhist
      ico(nn)=0
      pco(nn)=0
      histo(nn)=0.d0
      enddo

      do 1507 nn=1,n
      z=peso(nn)*contag(nn)
      if (z.gt.cort) go to 1507

      ko=1+(nhist-1)*(z-zmin)/razao
      ico(ko)=ico(ko)+1
      pco(ko)=pco(ko)+peso(nn)
      histo(ko)=histo(ko)+z

 1507 continue

      do 1510 nn=1,nhist
      if (ico(nn).eq.0) go to 1510
      histo(nn)=histo(nn)/pco(nn)

 1510 continue

c
c     Pega o valor de frequencia mais alta e determina o valor
c     da contagem do planeta no anel
c

      do nn=1,nhist
      ior(nn)=nn
      nval(nn)=ico(nn)
      enddo
c
      call ordem (nhist,nhist,ior,nval)

c

      cmed=histo(ior(nhist))


 1515 continue


c
c     Corrige a luz do planeta no pixel
c


      matpix(j,i)=pixmat(j,i)-cmed



 498  continue
 499  continue


      return

      end






c
c     Subrotina fuplas
c
c     Calcula o valor de contagem da luz do planeta para cada isofota r
c
c     Versao modificada que roda mais rapido, com menos amostragem (grossura
c     de anel) que na versao robusta.
c
c
c     Last update:  M. Assafin  19/Jan/2011
c
c


      subroutine fuplas (bx,by,ix1,ix2,iy1,iy2,mazi)

      IMPLICIT REAL *8 (A-H,O-Z)
      real*4 MATPIX(5001,5001),PIXMAT(5001,5001),imagem(5001,5001),
     ?maximo,mazi
      
      dimension contag(50000),peso(50000),histo(30),ico(30),pco(30),
     ?ior(30),nval(30)

      COMMON/A3/IMAGEM
      COMMON/A4/PIXMAT
      COMMON/A5/MATPIX


c
c     Dados iniciais
c
c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

      nhist=30
      mpes=50000

      aux1=1.d0
      aux360=360.d0*grarad

c
c     Zera vetores
c

      do i=1,mpes
      contag(i)=0.d0
      peso(i)=0.d0
      enddo



      do 499 i=iy1,iy2
      do 498 j=ix1,ix2

c
c     Para cada pixel, ha' um unico anel com um unico raio R. Portanto, para
c     cada pixel, calcula-se o anel correspondente, isto e', o valor da
c     contagem media real do planeta correspondente a esse R.
c
c


      xp=j
      yp=i

      raiopp=dsqrt((xp-bx)**2+(yp-by)**2) 

      ixmin=bx-raiopp
      ixmax=bx+raiopp

      if (ixmin.lt.ix1) ixmin=ix1
      if (ixmax.gt.ix2) ixmax=ix2


      n=0

      do 1101 jj=ixmin,ixmax

      x=jj

      d=dsqrt(dabs(raiopp**2-(x-bx)**2))

      y=by+d

      ii=y


      if (imagem(jj,ii).ge.mazi) go to 1100

      if (ii.lt.iy1) go to 1100 
      if (ii.gt.iy2) go to 1100 

      raio=dsqrt((x-bx)**2+(y-by)**2)

      dd=dabs(raio-raiopp)

      n=n+1
      contag(n)=pixmat(jj,ii)
      peso(n)=1.d0/(dd+1.d0)**2


c

 1100 continue

      y=by-d

      ii=y

      if (imagem(jj,ii).ge.mazi) go to 1101


      if (ii.lt.iy1) go to 1101 
      if (ii.gt.iy2) go to 1101 

      raio=dsqrt((x-bx)**2+(y-by)**2)

      dd=dabs(raio-raiopp)

      n=n+1
      contag(n)=pixmat(jj,ii)
      peso(n)=1.d0/(dd+1.d0)**2


 1101 continue




c
c     Estatistica para eliminar contagens altas associadas a satelites,
C     estrelas intrusas e cruzes de difracao,
c


c     if (n.lt.1) then
c     matpix(j,i)=0.
c     go to 498
c     endif

      if (n.lt.1) go to 498

c
c     Calcula histograma, qdo ha bastante pixel
c

 1500 continue


      zmin=1.d14
      zmax=-1.d14


      do nn=1,n
      z=peso(nn)*contag(nn)
      if (z.gt.zmax) zmax=z
      if (z.lt.zmin) zmin=z
      enddo

      razao=zmax-zmin

      if (razao.lt.0.1d0) then
      cmed=zmax
      go to 1515
      endif


      do nn=1,nhist
      ico(nn)=0
      pco(nn)=0
      histo(nn)=0.d0
      enddo

      do nn=1,n
      z=peso(nn)*contag(nn)
      ko=1+(nhist-1)*(z-zmin)/razao
      ico(ko)=ico(ko)+1
      pco(ko)=pco(ko)+peso(nn)
      histo(ko)=histo(ko)+z
      enddo

      do 1501 nn=1,nhist
      if (ico(nn).eq.0) go to 1501
      histo(nn)=histo(nn)/pco(nn)

 1501 continue

c
c     Refina bins do histograma
c

      ncort=0.8d0*n

      mcort=0
      do nn=1,nhist
      mcort=mcort+ico(nn)
      if (mcort.gt.ncort) go to 1505
      enddo

 1505 cort=histo(nn)

      razao=cort-zmin

      if (razao.lt.0.1d0) then
      cmed=cort
      go to 1515
      endif


      do nn=1,nhist
      ico(nn)=0
      pco(nn)=0
      histo(nn)=0.d0
      enddo

      do 1507 nn=1,n
      z=peso(nn)*contag(nn)
      if (z.gt.cort) go to 1507
      ko=1+(nhist-1)*(z-zmin)/razao
      ico(ko)=ico(ko)+1
      pco(ko)=pco(ko)+peso(nn)
      histo(ko)=histo(ko)+z

 1507 continue

      do 1510 nn=1,nhist
      if (ico(nn).eq.0) go to 1510
      histo(nn)=histo(nn)/pco(nn)


 1510 continue

c
c     Pega o valor de frequencia mais alta e determina o valor
c     da contagem do planeta no anel
c



      do nn=1,nhist
      ior(nn)=nn
      nval(nn)=ico(nn)
      enddo
c
 
      call ordem (nhist,nhist,ior,nval)

c

      cmed=histo(ior(nhist))


 1515 continue


c
c     Corrige a luz do planeta no pixel
c


      matpix(j,i)=pixmat(j,i)-cmed


 498  continue
 499  continue

      return

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


      character*50 infits
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
      if (ler(j+i).eq.ibr .or. ler(j+i).eq.ibar .or. ler(j+i).eq.iplic)
     ? go to 40
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

      return
      end




c
c
c
c     Subrotina wfitsp
C
C     Adaptacao da wfits, essa subrotina escreve uma imagem fits a partir
C     da matriz de pixels dada em real*4 ou integer*2, swapada ou nao.
C
c
c     Last modified:   M. Assafin  25/Jan/2011
c
c


      subroutine wfitsp (if,ipfits,bitpix,iswap,nx,ny,bscale,bzero,
     ?kheads)


      implicit real*8 (a-h,o-z)

      integer*2 iwork2(1440)
      integer*4 iwork4(720)
      integer*8 iwork8(360)
      real*4    rwork4(720)
      real*8    rwork8(360)

      integer*1 swork(2880),iby4(4)
      integer*2 bitpix

      real*4 pixmat(5001,5001)

      character*58 ipfits
      character*2880 header,head(10)

      icab(il,ic)=(il-1)*80+ic


      common /a4/pixmat


c
 
      makhea=10

c

      do kk=1,makhea
      head(kk)=''
      enddo


c
c     Initial data
c

      nbytes=2880

      if (bitpix.eq.16)  ibytes=2
      if (bitpix.eq.32)  ibytes=4
      if (bitpix.eq.-32) ibytes=4
      if (bitpix.eq.64)  ibytes=8
      if (bitpix.eq.-64) ibytes=8
      

c

      kwork=nbytes/ibytes



c
c     Opens fits file
c 

      
      open(if,file=ipfits,access='direct',form='unformatted',recl=2880)



c
c     Updates fits header
c


      do kk=1,kheads

      read (if,rec=kk) head(kk)

      enddo


      close (if)


      open(if,file=ipfits,access='direct',form='unformatted',recl=2880)

c

      do kk=1,kheads

      header=head(kk)


      do k=1,36

      ic1=1
      ic2=8
      ip1=icab(k,ic1)
      ip2=icab(k,ic2)
     


      if(header(ip1:ip2).eq.'SIMPLE  '.or.header(ip1:ip2).eq.'simple  ')
     ?then
      jc1=1
      jc2=46
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='SIMPLE  =                    T / Fits Standard'
      endif



      if(header(ip1:ip2).eq.'BITPIX  '.or.header(ip1:ip2).eq.'bitpix  ')
     ?then
      jc1=1
      jc2=47
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='BITPIX  =                      / Bits per pixel'
      jc1=28
      jc2=30
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      write (header(jp1:jp2),'(i3)') bitpix
      endif



      if(header(ip1:ip2).eq.'NAXIS   '.or.header(ip1:ip2).eq.'naxis   ')
     ?then
      jc1=1
      jc2=47
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='NAXIS   =                    2 / Number of axes'
      endif



      if(header(ip1:ip2).eq.'NAXIS1  '.or.header(ip1:ip2).eq.'naxis1  ')
     ?then
      jc1=1
      jc2=44
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='NAXIS1  =                      / Axis Length'
      jc1=26
      jc2=30
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      write (header(jp1:jp2),'(i5)') nx
      endif



      if(header(ip1:ip2).eq.'NAXIS2  '.or.header(ip1:ip2).eq.'naxis2  ')
     ?then
      jc1=1
      jc2=44
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='NAXIS2  =                      / Axis Length'
      jc1=26
      jc2=30
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      write (header(jp1:jp2),'(i5)') ny
      endif



      if(header(ip1:ip2).eq.'BSCALE  '.or.header(ip1:ip2).eq.'bscale  ')
     ?then
      jc1=1
      jc2=44
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='BSCALE  =                      / Data scale '
      jc1=11 
      jc2=30
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      write (header(jp1:jp2),'(f20.10)') bscale
      endif


      if(header(ip1:ip2).eq.'BZERO   '.or.header(ip1:ip2).eq.'bzero   ')
     ?then
      jc1=1
      jc2=44
      jp1=icab(k,jc1)
      jp2=icab(k,jc2)
      header(jp1:jp2)='BZERO   =                      / Zero point '
      jc1=11 
      jc2=30
      jp1=icab(l,jc1)
      jp2=icab(l,jc2)
      write (header(jp1:jp2),'(f20.10)') bzero 
      endif

      enddo



      write (if,rec=kk) header

      enddo



c
c     Now writes the data
c

      irec=kheads

      m=0

      do i=1,ny
      do j=1,nx

      m=m+1

c
      if (bitpix.eq.16) then
      iwork2(m)=pixmat(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) iwork2
      endif
      endif
c
      if (bitpix.eq.32) then
      iwork4(m)=pixmat(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) iwork4
      endif
      endif
c
      if (bitpix.eq.64) then
      iwork8(m)=pixmat(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) iwork8
      endif
      endif
c
      if (bitpix.eq.-32) then
      rwork4(m)=pixmat(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) rwork4
      endif
      endif
c
      if (bitpix.eq.-64) then
      rwork8(m)=pixmat(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) rwork8
      endif
      endif
c


      if (m.eq.kwork) then
      m=0
      if (iswap.eq.2) then
      call swapo (if,ibytes,nbytes,irec,swork)
      write (if,rec=irec) swork
      endif
      endif

      enddo
      enddo


      if (m.eq.0) go to 50

      irec=irec+1

      if (bitpix.eq.16) write (if,rec=irec) (iwork2(ii),ii=1,m)
      if (bitpix.eq.32) write (if,rec=irec) (iwork4(ii),ii=1,m)
      if (bitpix.eq.64) write (if,rec=irec) (iwork8(ii),ii=1,m)
      if (bitpix.eq.-32) write (if,rec=irec) (rwork4(ii),ii=1,m)
      if (bitpix.eq.-64) write (if,rec=irec) (rwork8(ii),ii=1,m)

      if (iswap.eq.2) then
      nbytes=m*ibytes
      call swapo (if,ibytes,nbytes,irec,swork)
      write (if,rec=irec) (swork(ii),ii=1,nbytes)
      endif

 50   continue


      close (if)

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
c     
c     Subrotina rede
c
c     Constroi a malha de pixels e sub-pixels para encontrar o melhor
c     centro (x,y) do objeto brilhante, pelo menor residuo de luz
c
c
c     Atencao: centro (bx,by) -> 0.5 = meio do pixel 1
c                             -> 1.5 = meio do pixel 2, etc ...
c
c

c
c     Faz-se uma pesquisa binaria em torno do centro (bx,by) inicial do planeta
c     ateh achar o valor mais preciso para este centro. O centro (bx,by)
c     mais preciso corresponde aquele que resulta na menor soma de contagens
c     da imagem "rebatida", no limites (ix1,ix2,iy1,iy2) fixos.
c
c     A pesquisa comeca primeiro pixel a pixel dentro da area delimitada
c     no arquivo de entrada. Depois de determinado o pixel mais central,
c     comeca a pesquisa binaria dentro desse pixel para encontrar o centro
c     ideal dentro desse pixel. 
c
c
c
c     Atualizada por M. Assafin, 17/Abril/2009
c
c
c

      subroutine rede (maximo,ix1,ix2,iy1,iy2,bx,by,ladox,ladoy,
     ?malha,convcp,keycxy)
      IMPLICIT REAL *8 (A-H,O-Z)

      real*4 MATPIX(5001,5001),PIXMAT(5001,5001),imagem(5001,5001),
     ?maximo,malha(5001,5001),mazi

      COMMON/A3/IMAGEM
      COMMON/A4/PIXMAT
      COMMON/A5/MATPIX


c
      ipmax=5001

c

      mazi=maximo

c


      do i=1,ipmax
      do j=1,ipmax
      malha(j,i)=0.
      enddo
      enddo

c

      bx=bx+0.5d0
      by=by+0.5d0

      ibx=bx
      iby=by

      write (*,*) 'gausc (x,y) = ',ibx,iby


 700  continue
      soma=1.d14

c

      if (keycxy.eq.1) then


      do     iii=iby-ladoy,iby+ladoy
      do 702 jjj=ibx-ladox,ibx+ladox
      if (malha(jjj,iii).gt.0.1) go to 702
      bx=jjj
      by=iii
      call fupla (bx,by,ix1,ix2,iy1,iy2,mazi)

      conta=0.d0
      do     i=iy1,iy2
      do 701 j=ix1,ix2
      if (imagem(j,i).ge.maximo) go to 701
      conta=conta+abs(matpix(j,i))
 701  continue
      enddo

      malha(jjj,iii)=conta


 702  continue
      enddo


      else

      do     iii=iby-ladoy,iby+ladoy
      do 802 jjj=ibx-ladox,ibx+ladox
      if (malha(jjj,iii).gt.0.1) go to 802
      bx=jjj
      by=iii

      call fuplas (bx,by,ix1,ix2,iy1,iy2,mazi)

      conta=0.d0
      do     i=iy1,iy2
      do 801 j=ix1,ix2
      if (imagem(j,i).ge.maximo) go to 801
      conta=conta+abs(matpix(j,i))
 801  continue
      enddo

      malha(jjj,iii)=conta


 802  continue
      enddo


      endif



c

      do iii=iby-ladoy,iby+ladoy
      do jjj=ibx-ladox,ibx+ladox
      if (malha(jjj,iii).lt.soma) then
      soma=malha(jjj,iii)
      igx=jjj
      igy=iii
      endif
      enddo
      enddo

      if (ibx.eq.igx .and. iby.eq.igy) go to 705 


      ibx=igx
      iby=igy


      write (*,*) 'gausc (x,y) = ',ibx,iby

      go to 700

c
c     Pesquisa do centro em sub-pixel central
c

 705  continue

      write (*,*) 'gausc (x,y) = ',ibx,iby
      write (96,*)'gausc (x,y) = ',ibx,iby

      bx=ibx
      by=iby

      bx=bx+0.5d0
      by=by+0.5d0

      tudo=1.d14

      xlado=0.5d0
      ylado=0.5d0

 710  continue

      soma=1.d14

      xlado=xlado/2.d0
      ylado=ylado/2.d0

      bx=bx-3.d0*xlado
      by=by-3.d0*ylado

c

      if (keycxy.eq.1) then

      do iii=1,2
      by=by+2.d0*ylado
      do jjj=1,2
      bx=bx+2.d0*xlado

      call fupla (bx,by,ix1,ix2,iy1,iy2,mazi)

      conta=0.d0
      do     i=iy1,iy2
      do 720 j=ix1,ix2
      if (imagem(j,i).ge.maximo) go to 720
      conta=conta+abs(matpix(j,i))
 720  continue
      enddo

      if (conta.lt.soma) then
      soma=conta
      gx=bx
      gy=by
      endif

      enddo
      enddo


      else

      do iii=1,2
      by=by+2.d0*ylado
      do jjj=1,2
      bx=bx+2.d0*xlado

      call fuplas (bx,by,ix1,ix2,iy1,iy2,mazi)

      conta=0.d0
      do     i=iy1,iy2
      do 820 j=ix1,ix2
      if (imagem(j,i).ge.maximo) go to 820
      conta=conta+abs(matpix(j,i))
 820  continue
      enddo

      if (conta.lt.soma) then
      soma=conta
      gx=bx
      gy=by
      endif

      enddo
      enddo


      endif

c

      bx=gx
      by=gy



      dif=100.d0*dabs(tudo-soma)/tudo

      write (*,*) 'gausc (x,y), perc = ',bx,by,dif
      write (96,*)'gausc (x,y), perc = ',bx,by,dif

      if (dif.gt.convcp) then
      tudo=soma
      go to 710
      endif

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


      subroutine flatsk (nx,ny,ngrauf,mazi,fcmin,fcmax,xesto,yesto,po,
     ?cofo)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi,imagem

      dimension xesto(25010001),yesto(25010001),po(25010001),cofo(136)
      dimension pixmat(5001,5001),imagem(5001,5001)


      common /a3/imagem
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
      if (imagem(j,i).ge.mazi) go to 1
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
      if (imagem(j,i).ge.mazi) go to 3
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
c     Computation is made over the entire given region.
c
c
c
c     last update:  M. Assafin 12/Mar/2010
c
c


      subroutine skyb (nx,ny,mazi,malisa,vmin,fatceu,ceu,sigceu,izmax,
     ?threso)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi,imagem

      dimension ifundo(150000),jfundo(150000)
      dimension pixmat(5001,5001),imagem(5001,5001)

      common /a3/imagem
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
      if (imagem(j,i).ge.mazi) go to 10
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
c
c     Subroutine skyp
c
c
c     Determines sky background by the mode of the histogram of counts.
c
c     Computation is made over the perimeter of given region.
c
c
c     last update:  M. Assafin 23/Jul/2010
c
c


      subroutine skyp (ix1,ix2,iy1,iy2,mazi,malisa,vmin,fatceu,ceu,
     ?sigceu,izmax,threso)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi,imagem

      dimension ifundo(150000),jfundo(150000)
      dimension pixmat(5001,5001),imagem(5001,5001)

      common /a3/imagem
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

      i=iy1
      do 10 j=ix1,ix2
      if (imagem(j,i).ge.mazi) go to 10
      n=pixmat(j,i)+xfu0
      jfundo(n)=jfundo(n)+1
      if (nfund.lt.n) nfund=n
  10  continue

      i=iy2
      do 11 j=ix1,ix2
      if (imagem(j,i).ge.mazi) go to 11
      n=pixmat(j,i)+xfu0
      jfundo(n)=jfundo(n)+1
      if (nfund.lt.n) nfund=n
  11  continue

      j=ix1
      do 12 i=iy1+1,iy2-1
      if (imagem(j,i).ge.mazi) go to 12
      n=pixmat(j,i)+xfu0
      jfundo(n)=jfundo(n)+1
      if (nfund.lt.n) nfund=n
  12  continue

      j=ix2
      do 13 i=iy1+1,iy2-1
      if (imagem(j,i).ge.mazi) go to 13
      n=pixmat(j,i)+xfu0
      jfundo(n)=jfundo(n)+1
      if (nfund.lt.n) nfund=n
  13  continue

      if (nfund.gt.kfund) nfund=kfund


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

c     if (ngrau.lt.1 .or. ngrau.gt.idimg) then

      if (ngrau.gt.idimg) then
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
      do  n=1,mgrau
      do  l=1,n
      icont=icont+1
      k=n-l
      polo=polo+cofo(icont)*(x**k)*(y**(l-1))
 
      enddo
      enddo


      return
      end


