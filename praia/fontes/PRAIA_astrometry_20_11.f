c
c     Programa PRAIA_astrometry
c
c     PROPOSITO
c
c
c     Dado um conjunto de imagens fits, identifica estrelas do campo,
c     determina medidas (x,y) e fotometria de objetos e fundo de ceu,
c     prepara arquivos de reducao alfa e delta. Tudo de forma automatica.
c
C
C     Atencao: imagens fits prontas para serem lidas em real *4
C
C
C
C     O algoritmo de calculo dos offsets (x,y) esta' modificado para ser
C     mais rapido. Agora, pega-se as N (ex: N=30) estrelas mais brilhantes
C     de cada campo e somente com elas aplica-se o algoritmo de reconhecimento
C     para calcular um offset (dx,dy) provisorio. Com esse offset, estende-se
C     o reconhecimento `as demais estrelas do campo, e obtem-se finalmente 
C     o offset dx,dy) definitivo, 
C
C
C     Os offsets provisorios sao tais que dx=x1-x2,dy=y1-y2, isto e' testa-se
C     nessa versao do programa diretamente a hipotese de que a estrela 1
C     e a estrela 2 nos frames 1 e 2 sao a mesma, para as N estrelas mais
C     brilhantes, e depois, a partir do offset provisorio, testa-se as
C     demais para calcular o offset (dx,dy) definitivo).
C
C
C     Os catalogos 2MASS e UCAC2 sao lidos diretamente.
C
C
C     A identificacao catalogo-medida e' feita com ajuste de 4 ctes par 
C     a par, com e sem reflexao em X. E' apenas preciso ser dada a 
C     escala, de resto a identificacao e' automatica. Sao usadas as
C     N (ex.: N=30) estrelas mais brilhantes medidas e as M estrelas mais
C     brilhantes de catalogo (ex:N=100) para a identificacao. Permite-se
C     que a escala de pixel ajustada seja ate' 20% diferente daquela
C     fornecida.
C     A identificacao e' feita primeiramente com o 2MASS. Depois de obter
C     alfas e deltas provisorios para as estrelas, a identificacao das 
C     estrelas do UCAC2 sao feitas diretamente por comparacao com esses alfas
C     e deltas provisorios
C
C
c     Nesta versao, identifica se imagem eh integer*2 ou integer*4, le imagem,
c     e guarda na matriz (real*4 como usual), com correcao de bscale e bzero,
c     dada pelo usuario ou lidos diretamente da imagem.
c
c
c     Le-se agora imagens fits de qq tipo (inteira, floating point, double,
c     littleendian ou bigendian) em fortran 77, sem auxilio de chamadas
c     externas como qfits e programas C.
c
c
c
c     Nesta versao, a curvatura de fundo de ceu eh flatada, melhorando a
c     identificacao de objetos e a fotometria.
c
c
c     Nesta versao, existe a opcao de fornecer um unico arquivo de mascara
c     de bad pixels para todas as imagens tratadas, alem das mascaras
c     individuais das imagens 
c
c
c
c
c     Nesta versao, o UCAC3 foi substituido pelo UCAC4.
c
c
c     Nesta nova versao, abre-se a possibilidade de acelerar a identificacao de
c     estrelas de catalogo. Caso a lista de entrada se referira a imagens com
c     mesma orientacao e escala de pixel, a primeira imagem eh tratada e entao,
c     os valores dos coeficientes de rotacao e escala sao usados nas identificoes
c     das imagens seguintes, acelerando o processo.
c
c
c     Nesta nova versao, usa-se a indexacao de estrelas dos catalogos 2MASS, UCAC2
c     e UCAC4, para acelerar o processo de leitura.
c
c
c     Nesta versao, reducoes (RA,Dec) sao feitas em relacao a um catalogo no formato
c     interno PRAIA.
c
c
c
c     Nesta versao, o usuario pode fornecer regioes para medida de objetos de 
c     interesse, sobrepondo-se a detecao automatica dos mesmos. Util qdo a detecao
c     automatica falha para determinados objetos-alvo. 
c
c
c     Nesta versao, imagens-traco sao reduzidas com o modelo da funcao-erro. A
c     regiao de medida das imagens-traco devem ser necessariamente fornecidas pelo
c     usuario. Para regioes "circle", ajustam-se Gaussianas circulares. Para regioes
c     "box", ajustam-se imagens-traco.
c
c
c
c      Last update: Marcelo Assafin - 21 Janeiro 2013
c   
c
c


      IMPLICIT REAL *8 (A-H,O-Z)
      parameter (stdin=5,idiobs=150000,ipmax=5001,icofsp=21,ng=10,
     ?idin=150,idin50=50,idin2=20000,iu2z=288,iu2i=240,iu4z=900,
     ?iu4i=1440,mpes=50000,nhist=30,nhisfw=10,jfdp=10000,kfund=150000)


      integer*4 imagem(ipmax,ipmax)
      real*4 PIXMAT(ipmax,ipmax),maximo,mazi
      real*4 pixel(ipmax,ipmax)
      integer*2 bitpix,bitpyx,betpix


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
   


      character*150 infits,imfits,names(idin2),ids9,ibadpx,kbadpx,imes,
     ?ires

      character*200 ired2m,irmp2m,irme2m,ireduc,iredu4,iredus

      character*50 lista,centro,ialvos,ialvo2,ialvou,ialvo4,ialvus,
     ?lred2m,lrmp2m,lrme2m,lreduc,lredu4,lredus

      character*50 iilvo2,iiivo2,fotrel,redred,redre4,
     ?rmpred,rmered,redrus,ifdp

      character *50 mraiz,uraiz,u4raiz,cpraia
      character *61 u2ind,u4ind

      character*50 subf

      character*1   menos,iver,nome(idin50),isig,ibrac
      character*69 iobold,ichobj,ihname(idin2),mchobj
      character*20 ichfil

      character*200 linha


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

      tt0=2.8d-3

      izero=0
      zero=0.d0

      d99=99.999d0

 
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

c     open (1,file='PRAIA_astrometry_20_10.dat')

      read (*,3) mraiz
      read (*,3) uraiz
      read (*,3) u4raiz
      read (*,3) cpraia
      read (*,*) iuserc

      read (*,*) krcat

      read (*,3) centro
      read (*,3) ialvos

      read (*,3) ifdp
      read (*,*) kfdp  
      read (*,3) kbadpx

      read (*,3) fotrel
      read (*,3) redred
      read (*,3) redre4
      read (*,3) rmpred
      read (*,3) rmered
      read (*,3) redrus


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
      read (*,*) scala
      read (*,*) ecala


      read (*,*) mix

      if (mix.ne.1 .and. mix.ne.2) mix=1


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
      read (*,*) icomax
      read (*,*) nbcat
      read (*,*) nbmed

      read (*,*) barx,bary

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
 1    format (23x,'PRAIA - astrometric and photometric setup')
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

      eoff=eoff/scala
      erpix=erpix/scala

      fmin=fmin/scala
      fmax=fmax/scala

      carx=grarad*barx/60.d0
      cary=grarad*bary/60.d0


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
c     Guarda index das estrelas do catalogo UCAC2 
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
c     Guarda index das estrelas do catalogo UCAC4 
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


      do 60 lllll=inicio,iultmo

      ired2m=''
      irmp2m=''
      ired2m=''
      irme2m=''
      ireduc=''
      iredu4=''
      iredus=''

      infits=''
      infits=names(lllll)


c
c     Prepara nomes e imagens
c

      do ii=1,idin
      if (infits(ii:ii+5).eq.'.fits ') go to 2000
      if (infits(ii:ii).eq.ibrac) go to 2000
      enddo

 2000 continue

      ii=ii-1

      kkii=ii

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

      do iii=50,1,-1
      if (lred2m(iii:iii).ne.ibrac) go to 2001
      enddo

 2001 continue

      ired2m(ii+2:ii+1+iii)=lred2m
      
c
c     2MASS tp
c

      do iii=50,1,-1
      if (lrmp2m(iii:iii).ne.ibrac) go to 2002
      enddo

 2002 continue

      irmp2m(ii+2:ii+1+iii)=lrmp2m
      
c
c     2MASS tp+pm
c

      do iii=50,1,-1
      if (lrme2m(iii:iii).ne.ibrac) go to 2003
      enddo

 2003 continue

      irme2m(ii+2:ii+1+iii)=lrme2m
      
c
c     UCAC2
c

      do iii=50,1,-1
      if (lreduc(iii:iii).ne.ibrac) go to 2004
      enddo

 2004 continue

      ireduc(ii+2:ii+1+iii)=lreduc

      
c
c     UCAC4
c

      do iii=50,1,-1
      if (lredu4(iii:iii).ne.ibrac) go to 2005
      enddo

 2005 continue

      iredu4(ii+2:ii+1+iii)=lredu4


      
c
c     User reference catalog
c

      do iii=50,1,-1
      if (lredus(iii:iii).ne.ibrac) go to 2006
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
c    Nome do arquivo de bad pixels
c

      ibadpx=''

      ibadpx(1:kkii)=infits(1:kkii)
      ibadpx(kkii+1:kkii+4)='.bpx'



c
c    Nome do arquivo de boxes para ds9 
c

      ids9=''

      ids9(1:kkii)=infits(1:kkii)
      ids9(kkii+1:kkii+4)='.reg'


c
c     Nome para leitura de regioes de targets fornecidos pelo usuario
c

      imes=''
      imes(1:kkii)=infits(1:kkii)
      imes(kkii+1:kkii+4)='.mes'
      ires=''
      ires=imes
      ires(kkii+1:kkii+4)='.res'



c
c     write (*,91) infits
c     write (*,93) ired2m
c     write (*,93) irmp2m
c     write (*,93) irme2m
c     write (*,93) ireduc
c     write (*,93) iredu4
c     write (*,93) iredus
c     write (*,91) ids9
c     write (*,91) ibadpx


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
      call refits (ipmax,pixmat,infits,nx,ny,nheads,ichobj,ipflag,
     ?bscale,bzero,kswap,iswap,bitpix)
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
c     Carrega targets diretamente fornecidos para medida pelo usuario
c

      call mesure (idiobs,imes,nmcir,xcir,ycir,lacir,nmtra,xtra,ytra,
     ?xlatra,ylatra,angtra)



c
c     Flattens sky background: 1rst step with star contamination
c



      if (ngrauf.ge.1 .and. ngrauf.le.15) then

      write (*,*)
      write (*,*) 'Flattening sky background. Please wait ...'
      write (*,*)

      fcmin=vmin
      fcmax=mazi

      call flatsk (ipmax,pixmat,nx,ny,ngrauf,mazi,fcmin,fcmax)

      endif




c
c     Computes typical sky background in the assumption of a flat sky
c


      call skyb (ipmax,kfund,pixmat,nx,ny,mazi,malisa,vmin,fatceu,ceu,
     ?sigceu,izmax,threso)

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

      call flatsk (ipmax,pixmat,nx,ny,ngrauf,mazi,fcmin,fcmax)

c
c     Re-computes typical sky background in the assumption of a flat sky
c


      call skyb (ipmax,kfund,pixmat,nx,ny,mazi,malisa,vmin,fatceu,ceu,
     ?sigceu,izmax,threso)

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


c
c     As margens do CCD sao excluidas por default
c

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
c     Marca para exclusao os pixels de targets diretamente fornecidos
c     pelo usuario. PSF Gaussiana circular
c



      do k=1,nmcir

      ix1=xcir(k)-lacir(k)
      ix2=xcir(k)+lacir(k)
      iy1=ycir(k)-lacir(k)
      iy2=ycir(k)+lacir(k)

      if (ix1.lt.1) ix1=1
      if (iy1.lt.1) iy1=1
      if (ix2.gt.nx) ix2=nx
      if (iy2.gt.ny) iy2=ny

      raio=lacir(k)
      xc=xcir(k)
      yc=ycir(k)

      do i=iy1,iy2
      do j=ix1,ix2

      call circul(raio,xc,yc,j,i,ichave)

      if (ichave.gt.0) imagem(j,i)=-1 

      enddo
      enddo

      enddo



c
c     Marca para exclusao os pixels de targets diretamente fornecidos
c     pelo usuario. PSF imagens traco.
c



      do k=1,nmtra


      ix1=xtra(k)-xlatra(k)
      ix2=xtra(k)+xlatra(k)
      iy1=ytra(k)-ylatra(k)
      iy2=ytra(k)+ylatra(k)

      cs=dcos(grarad*angtra(k))
      sn=dsin(grarad*angtra(k))

      do     i=iy1,iy2
      do 272 j=ix1,ix2

      ix=+(j-xtra(k))*cs+(y-ytra(k))*sn
      iy=-(j-xtra(k))*sn+(y-ytra(k))*cs


      if (ix.lt.1)  go to 272
      if (iy.lt.1)  go to 272
      if (ix.gt.nx) go to 272
      if (iy.gt.ny) go to 272

      imagem(ix,iy)=-1 


 272  continue

      enddo

      enddo


c
c     Abre arquivo regions do ds9.
c     Escreve o cabecalho do arquivo
c


      open (20,file=ids9)

      write (20,273)
 273  format('# Region file format: DS9 version 4.1')

      write (20,274) infits
 274  format('# Filename: ',a50)

      write (20,276)
 276  format('global color=green dashlist=8 3 width=1 font="helvetica 10
     ? normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 de
     ?lete=1 include=1 source=1')

      write (20,277)
 277  format('image')




c
c     Identifica objetos com pixels acima de "thres"
c


      call ident (ipmax,idiobs,nhist,pixmat,imagem,nx,ny,threso,mazi,
     ?xid,yid,idx1,idx2,idy1,idy2,idlado,npix,nstar)


c
c     Determina o centro Gaussiano (x,y) dos objetos identificados
c


      nest=0

      do 57 k=1,nstar

      if (npix(k).lt.ng) go to 57

c

      xc=xid(k)
      yc=yid(k)
      lado=idlado(k)

      ix1=idx1(k)
      ix2=idx2(k)
      iy1=idy1(k)
      iy2=idy2(k)

      xcj=xc
      ycj=yc
      zlado=2.d0*lado


c
c     Ajuste gaussiano: Gaussiana circular
c

      r=lado

      mazz=mazi



      call gcc (ipmax,icofsp,pixmat,ix1,ix2,iy1,iy2,r,mazi,xc,yc,ex,ey,
     ?it,sig,hh,fc)



      mazi=mazz


      lado=r

      if (xc.lt.ix1-0.5d0) go to 57
      if (xc.gt.ix2+0.5d0) go to 57
      if (yc.lt.iy1-0.5d0) go to 57
      if (yc.gt.iy2+0.5d0) go to 57



c
c     Guarda apenas candidatos com FWHM entre valores pre-fixados no
c     arquivo de entrada (ex: entre 1" e 4"), ie, com valores
c     de seeing astronomicamente possiveis

      fwhm=2.d0*sig*1.177410023d0

      if (fwhm.lt.fmin) go to 57
      if (fwhm.gt.fmax) go to 57

      seeing=fwhm*scala

c


      nest=nest+1

      xob(nest)=xc
      yob(nest)=yc
      ilado(nest)=lado
      seng(nest)=seeing

      altu(nest)=hh
      ialtu(nest)=nest

      exgcc(nest)=ex
      eygcc(nest)=ey
      sgcc(nest)=sig
      fgcc(nest)=fc


c

      j=xob(nest)
      i=yob(nest)

      ior(nest)=nest

      if (imagem(j,i).eq.-1) then
      nval(nest)=mazi+1
      else
      nval(nest)=pixmat(j,i)
      endif



      xlado=2.5d0*sgcc(nest)


      write (20,55) xob(nest),yob(nest),xlado
 55   format('circle(',2(f8.2,','),f8.2,')')


 57   continue




c
c     Determina centros das imagens normais (circulares) selecionadas
c     a mao pelo usuario
c



      do 59 k=1,nmcir


      xc=xcir(k)
      yc=ycir(k)
      r=lacir(k)

      ix1=xc-r
      ix2=xc+r
      iy1=yc-r
      iy2=yc+r

      


c
c     Ajuste gaussiano: Gaussiana circular
c


      mazz=mazi



      call gcc (ipmax,icofsp,pixmat,ix1,ix2,iy1,iy2,r,mazi,xc,yc,ex,ey,
     ?it,sig,hh,fc)



      mazi=mazz


      lado=r


      ix1=xc-r
      ix2=xc+r
      iy1=yc-r
      iy2=yc+r


      cfwhm=2.d0*sig*1.177410023d0

      seeing=cfwhm*scala


      nest=nest+1

      xob(nest)=xc
      yob(nest)=yc
      ilado(nest)=lado
      seng(nest)=seeing

      altu(nest)=hh
      ialtu(nest)=nest

      exgcc(nest)=ex
      eygcc(nest)=ey
      sgcc(nest)=sig
      fgcc(nest)=fc


c

      j=xob(nest)
      i=yob(nest)

      ior(nest)=nest

      if (imagem(j,i).eq.-1) then
      nval(nest)=mazi+1
      else
      nval(nest)=pixmat(j,i)
      endif


c     xlado=seng(nest)/(2.d0*scala)
  
      xlado=2.5d0*sgcc(nest)

      write (20,55) xob(nest),yob(nest),xlado


 59   continue




c
c     Minimum number of objects not reached ?
c

c     write (*,*) 'nnn = ',nest

c     stop


c     if (nest.lt.4) go to 60


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



      do n=1,nhisfw


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
c     Determines (x,y) of trace-image objects indicated by the user
c



      see=fwhm/(2.d0*1.177410023d0)

      see=see/scala


c
c     Escreve resultado de imagens traco em separado para analise
c

      open (66,file=ires)

      write(66,*)'altura,bx,by,sigma*scala,seeing      ,tetha,dlenght*sc
     ?ala,fundo,ex*scala,ey*scala,dlenght*altura'
      write (66,*)

      do 117 k=1,nmtra


      xc=xtra(k)
      yc=ytra(k)
      rx=xlatra(k)
      ry=ylatra(k)
      ang=angtra(k)

      


c
c     Ajuste da imagem-traco: modelo de ERF 
c


      mazz=mazi

c     declin=grarad*hmsgms(idg,idm,ds)
      
      dlen=iexps*15.d0/scala

c     dlen=dlen*dcos(declin)



      call trace (idiobs,ipmax,icofsp,pixmat,nx,ny,xc,yc,rx,ry,ang,
     ?maximo,see,ceu,dlen,icomax,altura,bx,by,sigma,tetha,dlenght,fundo,
     ?ex,ey)



      mazi=mazz


      lado=dlenght


      tfwhm=2.d0*sigma*1.177410023d0

      seeing=tfwhm*scala


      nest=nest+1

      xob(nest)=bx
      yob(nest)=by
      ilado(nest)=lado
      seng(nest)=seeing

      altu(nest)=altura*dlenght
      ialtu(nest)=nest

      exgcc(nest)=ex
      eygcc(nest)=ey
      sgcc(nest)=sigma
      fgcc(nest)=fundo

      tetha=tetha*radgra

      write (66,*) altura,bx,by,sigma*scala,seeing,tetha,
     ?dlenght*scala,fundo,ex*scala,ey*scala,dlenght*altura

c

      j=xob(nest)
      i=yob(nest)

      ior(nest)=nest

      if (imagem(j,i).eq.-1) then
      nval(nest)=mazi+1
      else
      nval(nest)=pixmat(j,i)
      endif


      xlado=dlenght
      ylado=sigma*3.d0

      write (20,116) xob(nest),yob(nest),xlado,ylado,tetha
 116  format('box(',4(f8.2,','),f8.2,')')


 117  continue

      close (20)

      close (66)


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
c     Relatorio parcial da fotometria
c

      ncort=0
      iconta=0
      dx=0.d0
      dy=0.d0

      open (99,file=fotrel)

 145  read (99,*,end=147)
      go to 145

 147  call backsp (2,nbac,99)

c


      write (99,150) fwhm,sigfwh,ceu,sigceu,iconta,dx,dy,ncort,nest,
     ?iexps,ichfil,infits,ichobj
 150  format(2(1x,f5.3),2(1x,f7.1),1x,i5,2(1x,f7.2),2(1x,i5),
     ?2x,i4,2x,a20,2x,a50,2x,a20)

      close (99)

c

      epoj=2000D0+(dj-2451545D0)/365.25D0

      rac=hmsgms(iah,iam,sa)*15.d0
      dec=hmsgms(idg,idm,ds)
      if (isig.eq.menos) dec=-dec

      drac=nx*scala/3600.d0
      ddec=ny*scala/3600.d0

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
c     Polo (sul ou norte) cai  no campo    -> krcat = 2: extrair
c     catalogos com projecao no plano tangente, estrela a estrela
c     (processo mais lento)
c
c     Polo (sul ou norte) nao cai no campo -> krcat = 1: extrair
c     catalogos com (RA,Dec) min_max, usando indexacao dos catalogos
c     (procedimento acelerado) 
c


      if (krcat.ne.2) then

      if (dabs(xx).le.areax .and. dabs(yy).le.areay) krcat=2

      endif


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

      call tmass (idiobs,mraiz,rac,dec,drac,ddec,ra2ma,de2ma,era2ma,
     ?ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,n2mass)

      else

      call stmass (idiobs,mraiz,rac,dec,drac,ddec,ra2ma,de2ma,era2ma,
     ?ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,n2mass)

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


      call ucac2 (idiobs,uraiz,epoj,rac,dec,drac,ddec,rauc2,deuc2,
     ?erauc2,edeuc2,pmra,pmde,epmra,epmde,udmgj,udmgh,udmgk,udmg,nucac2)

      else

      call sucac2 (idiobs,iu2z,iu2i,inuu2,ir1u2,ir2u2,uraiz,epoj,rac,
     ?dec,drac,ddec,rauc2,deuc2,erauc2,edeuc2,pmra,pmde,epmra,epmde,
     ?udmgj,udmgh,udmgk,udmg,nucac2)

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

      call ucac4 (idiobs,u4raiz,epoj,rac,dec,drac,ddec,rauc4,deuc4,
     ?erauc4,edeuc4,pmra4,pmde4,epmra4,epmde4,udmgj4,udmgh4,udmgk4,
     ?udmg4,nucac4)


      else


      call sucac4 (idiobs,iu4z,iu4i,inuu4,ir1u4,ir2u4,u4raiz,epoj,rac,
     ?dec,drac,ddec,rauc4,deuc4,erauc4,edeuc4,pmra4,pmde4,epmra4,epmde4,
     ?udmgj4,udmgh4,udmgk4,udmg4,nucac4)


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

      call cuser (idiobs,cpraia,epoj,rac,dec,drac,ddec,raucs,deucs,
     ?eraucs,edeucs,pmras,pmdes,epmras,epmdes,udmgjs,udmghs,udmgks,
     ?udmgs,nucaus)

      endif




c
c     Identifica estrelas 2MASS e UCAC2 comuns para a regiao de captura
c     dos catalogos, isto eh, tomando alem da regiao do CCD, a regiao
c     extra do entorno (default lados da regiao extra = 2 vezes o lado
c     X e Y do CCD)
c
c     Esta identificacao serah util nas reducoes com 2MASS corrigido
c     pelo UCAC2 para regioes onde o UCAC2 nao tem  estrelas suficientes
c     para reducao adequada.
c
c     Tambem serve para o calculo de magnitudes instrumentais usando o
c     magnitudes 2MASS convertidas J->V como referencia
c
c
 
      boxerr=box

      do i=1,idiobs
      iuc2ma(i)=0
      enddo


      do 650 i=1,nucac2
      do     j=1,n2mass

      ior(j)=j

      dx=dcos(deuc2(i)*grarad)*(ra2ma(j)-rauc2(i))*3600d0
      dy=(de2ma(j)-deuc2(i))*3600d0

      nval(j)=10*dsqrt(dx**2+dy**2)


      enddo

      call ordem (idiobs,n2mass,ior,nval)

      j=ior(1)


      dx=dcos(deuc2(i)*grarad)*(ra2ma(j)-rauc2(i))*3600d0
      dy=(de2ma(j)-deuc2(i))*3600d0

      xx=dsqrt(dx**2+dy**2)

      if (xx.gt.boxerr) go to 650

      iuc2ma(i)=j




 650  continue

c
c     Inicializando rotinas de estimativa de tempo de execucao da
c     reducao das imagens com o PRAIA
c

      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)

      tt0=tempop/3600.d0


c
c     Reconhece estrelas do 2MASS dentre as identificadas no campo
c



      if (nscl.eq.0) then

      call idxy2m (idiobs,icofsp,carx,cary,nval,ior,scala,erpix,rac,dec,
     ?nbcat,nbmed,nest,ialtu,xob,yob,n2mass,id2ma,ra2ma,de2ma,dmgj,xold,
     ?yold,xra2ma,yde2ma,ireflex,ecala,tt)

      else



      if (mix.eq.2) then

      call idxy2m (idiobs,icofsp,carx,cary,nval,ior,scala,erpix,rac,dec,
     ?nbcat,nbmed,nest,ialtu,xob,yob,n2mass,id2ma,ra2ma,de2ma,dmgj,xold,
     ?yold,xra2ma,yde2ma,ireflex,ecala,tt)

      else


      call mdxy2m (idiobs,icofsp,carx,cary,ngrau,coefxr,coefyr,kreflex,
     ?nval,ior,rcala,erpix,rac,dec,nbcat,nbmed,nest,ialtu,xob,yob,
     ?n2mass,id2ma,ra2ma,de2ma,dmgj,xold,yold,xra2ma,yde2ma,ireflex,
     ?escala,tt)

      endif


      endif


c


      if (ierro.eq.1) then
      ierro=0
      go to 62
      endif


c
c     Reducao (RA,DEC) com 2MASS sem correcoes com UCAC2
c


      do i=1,idiobs
      itira2(i)=0
      enddo


      call posred (idiobs,icofsp,ireflex,rac,dec,id2ma,n2mass,ra2ma,
     ?de2ma,nest,xob,yob,pcort2,ngrau,ngrau3,ngrau5,nstart,nfinal,
     ?xra2ma,yde2ma,era2,ede2,alfsi2,delsi2,alfre2,delre2,coefx2,coefy2,
     ?ecofx2,ecofy2,itira2,avam,dvam)



c
c     Melhora estimativa da escala de pixel e de seu erro, para opcao mix = 1
c     (imagens de mesma orientacao e escala de pixel)
c




      if (ierro.eq.0) then


      kreflex=ireflex


c
c     Atualiza valores dos coeficientes, imagem a imagem
c


      do mmmm=2,icofsp
      coefxr(mmmm)=(nscl*coefxr(mmmm)+coefx2(mmmm))/(nscl+1.d0)
      coefyr(mmmm)=(nscl*coefyr(mmmm)+coefy2(mmmm))/(nscl+1.d0)
      enddo



c
c     Calcula escala de pixel, apartir do ajuste, em radianos por pixel
c

  

      rscala=dsqrt(coefxr(2)**2+coefxr(3)**2)+dsqrt(coefyr(2)**2+
     ?coefyr(3)**2) 
      rscala=rscala/2.d0


      scl=scl+rscala
      scl2=scl2+rscala**2
      nscl=nscl+1



c
c     Estima erro da escala de pixel,  em radianos por pixel, pela variacao das
c     escalas encontradas nas imagens (atencao, util apenas para mix=1)
c


      if (nscl.eq.1) then

      rcala=scl

c     escala=ecala*grarad/3600.d0

      dnx=nx/2.d0
      dny=ny/2.d0

      escala=3.d0*dsqrt(avam**2+dvam**2)/dsqrt(dnx**2+dny**2)
      escala=escala*grarad/3600.d0


      else


      kscl=nscl
      tscl=scl
      tscl2=scl2
      
      call desvio (kscl,tscl,tscl2)

      rcala=tscl
      escala=tscl2*3.d0


      endif


      endif



c     write (*,*) 'i,num, esc, err = ',lllll,nscl,rcala*radgra*3600.d0,
c    ?escala*radgra*3600.d0


      ierro=0



c
c     Reconhece estrelas do UCAC2 dentre as identificadas no campo
c


      do i=1,idiobs
      iduc2(i)=0
      enddo

c

      icouc2=0

      errou=erpix*scala

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


      write (*,*) 'UCAC2: 2x2 field size, extracted stars = ',nucac2
      write (*,*) 'UCAC2: identified stars in field       = ',icouc2



c
c     Reconhece estrelas do UCAC4 dentre as identificadas no campo
c


      do i=1,idiobs
      iduc4(i)=0
      enddo

c

      icouc4=0

      errou=erpix*scala

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
      write (*,*) 'UCAC4: 2x2 field size, extracted stars = ',nucac4
      write (*,*) 'UCAC4: identified stars in field       = ',icouc4




c
c     Reconhece estrelas do catalogo de referencia do usuario, dentre
c     as estrelas identificadas no campo
c


      if (iuserc.eq.1) then


      do i=1,idiobs
      iducs(i)=0
      enddo

c

      icoucs=0

      errou=erpix*scala

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
      write (*,*) 'USER: 2x2 field size, extracted stars = ',nucaus
      write (*,*) 'USER: identified stars in field       = ',icoucs


      endif


c
c     Reducao (RA,DEC) com UCAC2
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
c     Reducao (RA,DEC) com UCAC4
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
c     Reducao (RA,DEC) com catalogo de referencia do usuario
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
c     Escreve estatistica da reducao de alfa e delta de cada campo
c     num so arquivo, para reducao 2MASS e UCAC2
c     2MASS sem correcao de movimentos proprios UCAC2
c     

      if (nstart.eq.0) nstart=1
      if (nstaru.eq.0) nstaru=1

      if (nstart.lt.nfinal) nfinal=nstart
      if (nstaru.lt.nfinau) nfinau=nstaru

      perc=100.d0*(nstart-nfinal)/nstart
      percu=100.d0*(nstaru-nfinau)/nstaru

c
      open (97,file=redred)

 450  read (97,*,end=451)
      go to 450

 451  call backsp (2,nbac,97)

c

      write (97,455) alfsi2,delsi2,nstart,nfinal,perc,avam,
     ?dvam,alfsiu,delsiu,nstaru,nfinau,percu,avamu,dvamu,dj,iah,iam,
     ?sa,isig,idg,idm,ds,iexps,ichfil,infits,mchobj,nx,ny
 455  format(2(2(1x,f6.3),2(1x,i4),1x,f6.2,2(1x,f6.3),1x),1x,f16.8,
     ?1x,i2,1x,i2,1x,f7.4,1x,a1,i2,1x,i2,1x,f6.3,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))


      close (97)





c
c     Escreve estatistica da reducao de alfa e delta de cada campo
c     num so arquivo, para reducao 2MASS e catalogo de referencia do
c     usuario. 2MASS sem correcao de movimentos proprios UCAC2
c     



      if (iuserc.eq.1) then


      if (nstart.eq.0) nstart=1
      if (nstars.eq.0) nstars=1

      if (nstart.lt.nfinal) nfinal=nstart
      if (nstars.lt.nfinas) nfinas=nstars

      perc=100.d0*(nstart-nfinal)/nstart
      percs=100.d0*(nstars-nfinas)/nstars

c
      open (97,file=redrus)

 456  read (97,*,end=457)
      go to 456

 457  call backsp (2,nbac,97)

c

      write (97,455) alfsi2,delsi2,nstart,nfinal,perc,avam,
     ?dvam,alfsis,delsis,nstars,nfinas,percs,avams,dvams,dj,iah,iam,
     ?sa,isig,idg,idm,ds,iexps,ichfil,infits,mchobj,nx,ny


      close (97)


      endif





c
c     Escreve estatistica da reducao de alfa e delta de cada campo
c     num so arquivo, para reducao 2MASS e UCAC4.
c     2MASS sem correcao de movimentos proprios UCAC2
c     

      if (nstart.eq.0) nstart=1
      if (nstar4.eq.0) nstar4=1

      if (nstart.lt.nfinal) nfinal=nstart
      if (nstar4.lt.nfina4) nfina4=nstar4

      perc=100.d0*(nstart-nfinal)/nstart
      perc4=100.d0*(nstar4-nfina4)/nstar4

c
      open (97,file=redre4)

 420  read (97,*,end=421)
      go to 420

 421  call backsp (2,nbac,97)

c

      write (97,455) alfsi2,delsi2,nstart,nfinal,perc,avam,
     ?dvam,alfsi4,delsi4,nstar4,nfina4,perc4,avam4,dvam4,dj,iah,iam,
     ?sa,isig,idg,idm,ds,iexps,ichfil,infits,mchobj,nx,ny


      close (97)






c
c    Calculo das magnitudes instrumentais baseado no sistema 2MASS
c
c    2MASS convertido de J para o sistema UCAC2 pela diferenca media
c    entre as magnitudes J do 2MASS e "VR" do UCAC2 para estrelas comuns
c    do campo inteiro (2CCDx por 2CCDy)
c


c
c    Encontra diferenca media entre J do 2MASS e "VR" do UCAC2
c


      resmg=0.d0
      res2mg=0.d0
      n=0
      do 652 i=1,nucac2
      if (iuc2ma(i).eq.0) go to 652
      if (udmg(i).gt.20.d0) go to 652
      j=iuc2ma(i)
      if (dmgj(j).gt.20.d0) go to 652
      n=n+1
      res=udmg(i)-dmgj(j)
      resmg=resmg+res
      res2mg=res2mg+res**2
 652  continue

      res=resmg/n

      if (n.gt.1) then
      res2mg=dsqrt((res2mg-2.d0*res*resmg+n*res**2)/(n-1.d0))
      else
      res2mg=dsqrt((res2mg-2.d0*res*resmg+n*res**2)/n)
      endif

      dm2mg=res
      dm2mge=res2mg

c
c    Calculo das magnitudes instrumentais e erro no
c    sistema "VR" 2MASS convertido via UCAC2
c

      zeromg=0.d0
      n=0
      do 660 i=1,nest
      if (id2ma(i).eq.0) go to 660
      if (altu(i).le.0.d0) go to 660
      j=id2ma(i)
      if (dmgj(j).gt.20.d0) go to 660
      n=n+1
      zeromg=zeromg+dmgj(j)+dm2mg+2.5d0*dlog10(altu(i))
 660  continue

      zeromg=zeromg/n


      resmg=0.d0
      res2mg=0.d0

      do 665 i=1,nest

      cudmg2(i)=zeromg-2.5d0*dlog10(altu(i))

      if (id2ma(i).eq.0) go to 665
      if (altu(i).le.0.d0) go to 665
      j=id2ma(i)
      if (dmgj(j).gt.20.d0) go to 665
      res=cudmg2(i)-dmgj(j)-dm2mg
      resmg=resmg+res
      res2mg=res2mg+res**2
 665  continue


      res=resmg/n
      if (n.gt.1) then
      res2mg=dsqrt((res2mg-2.d0*res*resmg+n*res**2)/(n-1.d0))
      else
      res2mg=dsqrt((res2mg-2.d0*res*resmg+n*res**2)/n)
      endif

      resmg2=res2mg


c
c     Magnitude instrumental 2MASS do fundo de ceu
c

      if (ceu.lt.0.1d0) then

      fumag2=d99

      else

      fumag2=zeromg

      endif


c
c    Calculo das magnitudes instrumentais e erro no
c    sistema UCAC2
c

      if (icouc2.eq.0) go to 462
      
c

      zeromg=0.d0
      n=0
      do 460 i=1,nest
      if (iduc2(i).eq.0) go to 460
      if (altu(i).le.0.d0) go to 460
c     if (fgcc(i).le.0.d0) go to 460
      j=iduc2(i)
      if (udmg(j).gt.20.d0) go to 460
      n=n+1
      zeromg=zeromg+udmg(j)+2.5d0*dlog10(altu(i))
 460  continue

      zeromg=zeromg/n

 462  continue

      resmg=0.d0
      res2mg=0.d0

      do 465 i=1,nest

      cudmg(i)=zeromg-2.5d0*dlog10(altu(i))

      if (iduc2(i).eq.0) go to 465
      if (altu(i).le.0.d0) go to 465
      j=iduc2(i)
      if (udmg(j).gt.20.d0) go to 465
      res=cudmg(i)-udmg(j)
      resmg=resmg+res
      res2mg=res2mg+res**2
 465  continue


      res=resmg/n
      if (n.gt.1) then
      res2mg=dsqrt((res2mg-2.d0*res*resmg+n*res**2)/(n-1.d0))

      else

      if (n.eq.1) then
      res2mg=dsqrt((res2mg-2.d0*res*resmg+n*res**2)/n)
      else
      res2mg=resmg2
      endif

      endif

c
c     Magnitude instrumental do fundo de ceu no sistema UCAC2
c

      if (ceu.lt.0.1d0) then

      fumag=d99

      else

      fumag=zeromg

      endif



c
c    Calculo das magnitudes instrumentais e erro no
c    sistema UCAC4
c

      if (icouc4.eq.0) go to 432
      
c

      z4romg=0.d0
      n=0
      do 430 i=1,nest
      if (iduc4(i).eq.0) go to 430
      if (altu(i).le.0.d0) go to 430
      j=iduc4(i)
      if (udmg4(j).gt.17.9d0) go to 430
      n=n+1
      z4romg=z4romg+udmg4(j)+2.5d0*dlog10(altu(i))
 430  continue

      z4romg=z4romg/n

 432  continue

      r4smg=0.d0
      r4s2mg=0.d0

      do 435 i=1,nest

      cudmg4(i)=z4romg-2.5d0*dlog10(altu(i))

      if (iduc4(i).eq.0) go to 435
      if (altu(i).le.0.d0) go to 435
      j=iduc4(i)
      if (udmg4(j).gt.17.9d0) go to 435
      res=cudmg4(i)-udmg4(j)
      r4smg=r4smg+res
      r4s2mg=r4s2mg+res**2
 435  continue


      res=r4smg/n
      if (n.gt.1) then
      r4s2mg=dsqrt((r4s2mg-2.d0*res*r4smg+n*res**2)/(n-1.d0))

      else

      if (n.eq.1) then
      r4s2mg=dsqrt((r4s2mg-2.d0*res*r4smg+n*res**2)/n)
      else
      r4s2mg=r4smg2
      endif

      endif

c
c     Magnitude instrumental do fundo de ceu no sistema UCAC4
c

      if (ceu.lt.0.1d0) then

      f4mag=d99

      else

      f4mag=z4romg

      endif




c
c    Calculo das magnitudes instrumentais e erro no
c    sistema do catalogo de referencia do usuario
c


      if (iuserc.eq.1) then


      if (icoucs.eq.0) go to 482
      
c

      zsromg=0.d0
      n=0
      do 480 i=1,nest
      if (iducs(i).eq.0) go to 480
      if (altu(i).le.0.d0) go to 480
      j=iducs(i)
      if (udmgs(j).gt.30.0d0) go to 480
      n=n+1
      zsromg=zsromg+udmgs(j)+2.5d0*dlog10(altu(i))
 480  continue

      zsromg=zsromg/n

 482  continue

      rssmg=0.d0
      rss2mg=0.d0

      do 485 i=1,nest

      cudmgs(i)=zsromg-2.5d0*dlog10(altu(i))

      if (iducs(i).eq.0) go to 485
      if (altu(i).le.0.d0) go to 485
      j=iducs(i)
      if (udmgs(j).gt.30.0d0) go to 485
      res=cudmgs(i)-udmgs(j)
      rssmg=rssmg+res
      rss2mg=rss2mg+res**2
 485  continue


      res=rssmg/n
      if (n.gt.1) then
      rss2mg=dsqrt((rss2mg-2.d0*res*rssmg+n*res**2)/(n-1.d0))

      else

      if (n.eq.1) then
      rss2mg=dsqrt((rss2mg-2.d0*res*rssmg+n*res**2)/n)
      else
      rss2mg=rssmg2
      endif

      endif

c
c     Magnitude instrumental do fundo de ceu no sistema do catalogo
c     de referencia do usuario
c

      if (ceu.lt.0.1d0) then

      fsmag=d99

      else

      fsmag=zsromg

      endif


      endif


c
c     Update do relatorio de fotometria com magnitude de fundo de
c     ceu e erro de ajuste de magnitude instrumental do campo usando
c     sistema UCAC2
c

c


      open (99,file=fotrel)

 500  read (99,*,end=501)
      go to 500

 501  call backsp (2,nbac,99)
      backspace 99
      

c
      read  (99,150) fwhm,sigfwh,ceu,sigceu,iconta,dx,dy,ncort,nest,
     ?iexps,ichfil,infits,ichobj

      backspace 99

      write (99,502) fwhm,sigfwh,ceu,sigceu,fumag,fumag2,res2mg,resmg2,
     ?dm2mg,dm2mge,iconta,dx,dy,ncort,nest,iexps,ichfil,infits,ichobj,
     ?nx,ny
 502  format(2(1x,f5.3),2(1x,f7.1),6(1x,f6.3),1x,i5,2(1x,f7.2),2(1x,i5),
     ?2x,i4,2x,a20,2x,a50,2x,a20,2(1x,i5))

      close (99)



c
c     Escreve resultados individuais das estrelas para cada campo,
c     alfa, delta, ajuste gaussiano (arquivos tipo *.xy)
c
c     2MASS reducao sem correcao de mov. prop. UCAC2
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

      do i=1,nest

      ex=scala*exgcc(i)
      ey=scala*eygcc(i)
      sig=scala*sgcc(i)

      pma=d99
      pmd=d99
      epma=d99
      epmd=d99

      pma4=d99
      pmd4=d99
      epma4=d99
      epmd4=d99

      pmas=d99
      pmds=d99
      epmas=d99
      epmds=d99


      xmgj=d99
      xmgh=d99
      xmgk=d99
      xmgu=d99
      ermgj=d99
      ermgh=d99
      ermgk=d99

      xmgu4=d99
      r4smg2=d99

      xmgus=d99
      rssmg2=d99


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

      xmgj=dmgj(j)
      xmgh=dmgh(j)
      xmgk=dmgk(j)
      ermgj=emgj(j)
      ermgh=emgh(j)
      ermgk=emgk(j)
      ktira=itira2(jj)
      alsi2=alfre2(jj)
      desi2=delre2(jj)
      endif


      if (iduc2(i).ne.0) then
      kk=kk+1
      k=iduc2(i)

c     write (*,*) 'ucac2 kk k = ',kk,k

      xmgu=udmg(k)
      pma=pmra(k)
      pmd=pmde(k)
      epma=epmra(k)
      epmd=epmde(k)
      ktirau=itirau(kk)
      alsiu=alfreu(kk)
      desiu=delreu(kk)
      endif



c
c     2MASS
c


      ra=xra2ma(i)/15.d0
      de=yde2ma(i)


      write (64,470) xob(i),yob(i),seng(i),altu(i),fgcc(i),fumag,fumag2,
     ?xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era2(i),ede2(i),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


c
c    UCAC2
c

      ra=xrauc2(i)/15.d0
      de=ydeuc2(i)


      write (65,470) xob(i),yob(i),seng(i),altu(i),fgcc(i),fumag,fumag2,
     ?xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
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

      xmgu4=udmg4(m)

      xmgj=udmgj4(m)
      xmgh=udmgh4(m)
      xmgk=udmgk4(m)

      pma4=pmra4(m)
      pmd4=pmde4(m)
      epma4=epmra4(m)
      epmd4=epmde4(m)
      ktira4=itira4(mm)
      alsi4=alfre4(mm)
      desi4=delre4(mm)
      endif


      ra=xrauc4(i)/15.d0
      de=ydeuc4(i)

      cucu=cudmg2(i)+f4mag-fumag

      fufu=fumag2+f4mag-fumag

      write (66,470) xob(i),yob(i),seng(i),altu(i),fgcc(i),f4mag,fufu,
     ?xmgu4,cudmg4(i),cucu,xmgj,xmgh,xmgk,r4s2mg,r4smg2,ermgj,ermgh,
     ?ermgk,pma4,pmd4,epma4,epmd4,ex,ey,era4(i),ede4(i),alfsi4,delsi4,
     ?nstar4,nfina4,alsi4,desi4,ktira4,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


c
c     Catalogo de referencia do usuario
c


      if (iuserc.eq.1) then

      ra=xraucs(i)/15.d0
      de=ydeucs(i)

      cucu=cudmgs(i)+fsmag-fumag

      fufu=fumag2+fsmag-fumag


      if (iducs(i).ne.0) then
      nn=nn+1
      n=iducs(i)

c     write (*,*) 'user nn n = ',nn,n

      xmgus=udmgs(n)

c     xmgj=udmgjs(n)
c     xmgh=udmghs(n)
c     xmgk=udmgks(n)

      pmas=pmras(n)
      pmds=pmdes(n)
      epmas=epmras(n)
      epmds=epmdes(n)
      ktiras=itiras(nn)
      alsis=alfres(nn)
      desis=delres(nn)
      endif




      write (67,470) xob(i),yob(i),seng(i),altu(i),fgcc(i),fsmag,fufu,
     ?xmgus,cudmgs(i),cucu,xmgj,xmgh,xmgk,rss2mg,rssmg2,ermgj,ermgh,
     ?ermgk,pmas,pmds,epmas,epmds,ex,ey,eras(i),edes(i),alfsis,delsis,
     ?nstars,nfinas,alsis,desis,ktiras,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny


      endif




      enddo

      close (64)
      close (65)
      close (66)
      close (67)


c
c     Faz a estatistica de alvos, com "O-A" (Observado menos Alvo), para
c     as reducoes com 2MASS, UCAC2 e UCAC4
c
c     2MASS reducao sem correcao de mov. prop. UCAC2
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
c     Correcao do 2MASS ao sistema UCAC2 "em J2000" (mais precisamente,
c     na epoca individual de observacao de cada estrela 2MASS; para nao
c     carregar na notacao, vamos nos referir a essa correcao "em J2000"
c     ja que as epocas sao proximas, mas os calculos sao exatos para as
c     epocas individuais)
c
c     Correcao feita com a tecnica do plano tangente, aplicando
c     correcao de origem, rotacao de eixos e ajuste uniforme de escala,
c     i.e., 4ctes
c
c


c
c     Nesta versao, o tamanho do campo para coleta das estrelas comuns
c     eh determinado pelo usuario. Pode ir do tamanho do CCD ateh o
c     dobro deste tamanho. 
c


c
c     Primeiro passo: traz posicoes das estrelas UCAC2 comuns para a
c     epoca individual das estrelas 2MASS aplicando correcao de
c     movimentos proprios
c




      do 600 k=1,nucac2

      if (iuc2ma(k).eq.0) go to 600
    
      j=iuc2ma(k)

      epo2m=2000D0+(ddj2(j)-2451545D0)/365.25D0

      dt=epo2m-epoj


      deuc2(k)=deuc2(k)+ pmde(k)*dt/3600.d0
      rauc2(k)=rauc2(k)+(pmra(k)*dt/3600.d0)/dabs(dcos(deuc2(k)*grarad))



 600  continue

c
c     Correcao de origem, rotacao de eixos, e ajuste uniforme de escala
c     (modelo 4 ctes). Ha eliminacao de estrelas, dentro da box de erro
c     estipulada para reducoes UCAC2.
c
c     Correcao feita com a tecnica do plano tangente
c

      ucor=(pcortu/3600.d0)*grarad


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

 603  continue


      ncom=0

      do 605 k=1,nucac2

      if (iuc2ma(k).eq.0) go to 605
    
      if (itirau(k).ne.0) go to 605

      j=iuc2ma(k)

c
c     Projeta catalogos no plano tangente p/ epoca individual das
c     estrelas 2MASS 
c 


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

      bra=grarad*rauc2(k)
      bde=grarad*deuc2(k)

      d=DEXY(bra,bde,grac,gdec)
      xpa(ncom)=XPAD(bra,bde,grac)/d
      ypa(ncom)=YPAD(bra,bde,grac,gdec)/d

 605  continue

c
c     Calcula coeficientes de 4ctes para estrelas comuns
c

      do i=1,icofsp
      cocomx(i)=0.d0
      cocomy(i)=0.d0
      enddo


      if (ncom.lt.3) then
      write (*,*) 'Less than 3 common 2MASS & UCAC2 stars.'
      write (*,*) 'Tangent plan correction not possible.'
      go to 62
      endif


      call isol (idiobs,icofsp,ngraup,ncom,xmed,ymed,xpa,ypa,cocomx,
     ?cocomy)

      if (ierro.eq.1) then
      ierro=0
      write (*,*) 'Tangent plan correction not possible.'
      go to 62
      endif





c
c     Aplica correcao do plano tangente (4 Ctes) a todo o 2MASS
c     (ie, dentro do campo CCD apenas, pois eh tudo do que se necessita) 
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
c     Determina movimentos proprios tipicos, tomando as estrelas UCAC2
c     mais fracas em magnitude, que nao tenham sido eliminadas na
c     correcao de origem e rotacao do 2MASS com o ajuste de 4 Ctes da
c     tecnica do plano tangente
c
c     Nesta versao, toma-se as estrelas comuns do campao determinado pelo
c     usuario.
c

      nmp=0

      do 625 k=1,nucac2



      if (iuc2ma(k).eq.0) go to 625
    
      if (itirau(k).ne.0) go to 625

      bra=grarad*rauc2(k)
      bde=grarad*deuc2(k)

      d=DEXY(bra,bde,grac,gdec)
      xxpa=XPAD(bra,bde,grac)/d
      yypa=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xxpa).gt.areax) go to 625
      if (dabs(yypa).gt.areay) go to 625

      nmp=nmp+1

      ior(nmp)=k
      nval(nmp)=10000*udmg(k)

 625  continue

      call ordem (idiobs,nmp,ior,nval)

      ncort=nmp/2.d0

      xmov=0.d0
      ymov=0.d0

      do i=1,ncort

      k=ior(nmp-i+1)
      xmov=xmov+pmra(k)
      ymov=ymov+pmde(k)

      enddo

      xmov=xmov/ncort
      ymov=ymov/ncort


c
c     Aplicando movimentos proprios a estrelas 2MASS comuns ao UCAC2,
c     trazendo essas estrelas para a epoca do campo observado
c
c     Soh eh necessario fazer essa correcao dentro do campo CCD
c


      do 630 i=1,nest

      if (id2ma(i).eq.0) go to 630
      if (iduc2(i).eq.0) go to 630
    
      j=id2ma(i)
      k=iduc2(i)


      epo2m=2000D0+(ddj2(j)-2451545D0)/365.25D0

      dt=epoj-epo2m

      de2ma(j)=de2ma(j)+ pmde(k)*dt/3600.d0
      ra2ma(j)=ra2ma(j)+(pmra(k)*dt/3600.d0)/dabs(dcos(de2ma(j)*grarad))



 630  continue


c
c     Reducao (RA,DEC) com 2MASS corrigido
c
c     2MASS previamente corrigido em J2000 pelo UCAC2
c
c     2MASS tambem com correcao de movimentos proprios aplicados
c     apenas as estrelas comuns ao UCAC2, sem aplicacao de m.p.
c     para as demais estrelas 2MASS
c
c     Corte de (O-C) agora igual ao do UCAC2
c     
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
c     Escreve estatistica da reducao de alfa e delta de cada campo
c     num so arquivo, para reducao 2MASS corrigido e UCAC2
c
c     2MASS previamente corrigido em J2000 pelo UCAC2
c     2MASS tambem com correcao de movimentos proprios aplicados
c     apenas as estrelas comuns ao UCAC2, sem aplicacao de m.p.
c     para as demais estrelas 2MASS
c     

      if (nstart.eq.0) nstart=1
      if (nstart.lt.nfinal) nfinal=nstart

      perc=100.d0*(nstart-nfinal)/nstart

c
      open (97,file=rmpred)

 510  read (97,*,end=511)
      go to 510

 511  call backsp (2,nbac,97)

c

      write (97,455) alfsi2,delsi2,nstart,nfinal,perc,avam,
     ?dvam,alfsiu,delsiu,nstaru,nfinau,percu,avamu,dvamu,dj,iah,iam,
     ?sa,isig,idg,idm,ds,iexps,ichfil,infits,mchobj,nx,ny

      close (97)

c
c     Escreve resultados individuais das estrelas para cada campo,
c     alfa, delta, ajuste gaussiano (arquivos tipo *.xy)
c
c     2MASS previamente corrigido em J2000 pelo UCAC2
c     2MASS tambem com correcao de movimentos proprios aplicados
c     apenas as estrelas comuns ao UCAC2, sem aplicacao de m.p.
c     para as demais estrelas 2MASS
c     
c

      open (64,file=irmp2m)

      jj=0
      kk=0

      do i=1,nest

      ex=scala*exgcc(i)
      ey=scala*eygcc(i)
      sig=scala*sgcc(i)

      pma=d99
      pmd=d99
      epma=d99
      epmd=d99

      xmgj=d99
      xmgh=d99
      xmgk=d99
      xmgu=d99
      ermgj=d99
      ermgh=d99
      ermgk=d99

      ktira=99


      alsi2=d99
      desi2=d99


      if (id2ma(i).ne.0) then
      jj=jj+1
      j=id2ma(i)

c     write (*,*) '2mass jj j = ',jj,j

      xmgj=dmgj(j)
      xmgh=dmgh(j)
      xmgk=dmgk(j)
      ermgj=emgj(j)
      ermgh=emgh(j)
      ermgk=emgk(j)
      ktira=itira2(jj)
      alsi2=alfre2(jj)
      desi2=delre2(jj)
      endif


      if (iduc2(i).ne.0) then
      kk=kk+1
      k=iduc2(i)

c     write (*,*) 'ucac2 kk k = ',kk,k

      xmgu=udmg(k)
      pma=pmra(k)
      pmd=pmde(k)
      epma=epmra(k)
      epmd=epmde(k)
      endif

      ra=xra2ma(i)/15.d0
      de=yde2ma(i)


      write (64,470) xob(i),yob(i),seng(i),altu(i),fgcc(i),fumag,fumag2,
     ?xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era2(i),ede2(i),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny




      enddo

      close (64)


c
c     Faz a estatistica de alvos, com "O-A" (Observado menos Alvo), para
c     as reducoes com 2MASS corrigido
c
c     2MASS previamente corrigido em J2000 pelo UCAC2
c     2MASS tambem com correcao de movimentos proprios aplicados
c     apenas as estrelas comuns ao UCAC2, sem aplicacao de m.p.
c     para as demais estrelas 2MASS
c

      write (*,*) 
      write (*,*) 
      write(*,*)'**** 2MASS + tangent plan correction + common pm  ****'
      write (*,*) 

      call estat (box,ialvos,irmp2m,iilvo2)



c
c     Aplicando o movimento proprio medio encontrado para as 50% estrelas
c     UCAC2 mais fracas (nao eliminadas da correcao de origem e rotacao
c     do 2MASS), para as N estrelas 2MASS mais brilhantes (magnitude
c     instrumental) que nao sejam comuns ao UCAC2 
c

      ncom=0
      kk=0
      mm=0

      do 635 i=1,nest

      if (id2ma(i).eq.0) go to 635

      ncom=ncom+1

      if (iduc2(i).ne.0) then
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
c     Reducao (RA,DEC) com 2MASS corrigido
c
c     2MASS previamente corrigido em J2000 pelo UCAC2
c
c     2MASS tambem com correcao de movimentos proprios aplicados
c     para as estrelas comuns ao UCAC2, e com aplicacao de m.p.
c     medio para as N estrelas 2MASS mais brilhantes (pela magnitude
c     instrumental), com descarte das demais estrelas 2MASS mais
c     fracas da reducao
c
c     Corte de (O-C) agora igual ao do UCAC2
c

      call posred (idiobs,icofsp,ireflex,rac,dec,id2ma,n2mass,ra2ma,
     ?de2ma,nest,xob,yob,pcortu,ngrau,ngrau3,ngrau5,nstart,nfinal,
     ?xra2ma,yde2ma,era2,ede2,alfsi2,delsi2,alfre2,delre2,coefx2,coefy2,
     ?ecofx2,ecofy2,itira2,avam,dvam)

      ierro=0


c
c     Escreve estatistica da reducao de alfa e delta de cada campo
c     num so arquivo, para reducao 2MASS corrigido e UCAC2
c
c     2MASS tambem com correcao de movimentos proprios aplicados
c     para as estrelas comuns ao UCAC2, e com aplicacao de m.p.
c     medio para as N estrelas 2MASS mais brilhantes (pela magnitude
c     instrumental), com descarte das demais estrelas 2MASS mais
c     fracas da reducao
c
c     

      ntira=nstart-nfinal

      nstart=mm+nume
      nfinal=nstart-ntira

      if (nstart.eq.0) nstart=1
      if (nstart.lt.nfinal) nfinal=nstart

      perc=100.d0*(nstart-nfinal)/nstart

c
      open (97,file=rmered)

 520  read (97,*,end=521)
      go to 520

 521  call backsp (2,nbac,97)

c

      write (97,455) alfsi2,delsi2,nstart,nfinal,perc,avam,
     ?dvam,alfsiu,delsiu,nstaru,nfinau,percu,avamu,dvamu,dj,iah,iam,
     ?sa,isig,idg,idm,ds,iexps,ichfil,infits,mchobj,nx,ny

      close (97)

c
c     Escreve resultados individuais das estrelas para cada campo,
c     alfa, delta, ajuste gaussiano (arquivos tipo *.xy)
c
c     2MASS tambem com correcao de movimentos proprios aplicados
c     para as estrelas comuns ao UCAC2, e com aplicacao de m.p.
c     medio para as N estrelas 2MASS mais brilhantes (pela magnitude
c     instrumental), com descarte das demais estrelas 2MASS mais
c     fracas da reducao
c     
c

      open (64,file=irme2m)

      jj=0
      kk=0

      do i=1,nest

      ex=scala*exgcc(i)
      ey=scala*eygcc(i)
      sig=scala*sgcc(i)

      pma=d99
      pmd=d99
      epma=d99
      epmd=d99

      xmgj=d99
      xmgh=d99
      xmgk=d99
      xmgu=d99
      ermgj=d99
      ermgh=d99
      ermgk=d99

      ktira=99


      alsi2=d99
      desi2=d99


      if (id2ma(i).ne.0) then
      jj=jj+1
      j=id2ma(i)

c     write (*,*) '2mass jj j = ',jj,j

      xmgj=dmgj(j)
      xmgh=dmgh(j)
      xmgk=dmgk(j)
      ermgj=emgj(j)
      ermgh=emgh(j)
      ermgk=emgk(j)
      ktira=itira2(jj)
      alsi2=alfre2(jj)
      desi2=delre2(jj)
      endif


      if (iduc2(i).ne.0) then
      kk=kk+1
      k=iduc2(i)

c     write (*,*) 'ucac2 kk k = ',kk,k

      xmgu=udmg(k)
      pma=pmra(k)
      pmd=pmde(k)
      epma=epmra(k)
      epmd=epmde(k)
      endif

      ra=xra2ma(i)/15.d0
      de=yde2ma(i)


      write (64,470) xob(i),yob(i),seng(i),altu(i),fgcc(i),fumag,fumag2,
     ?xmgu,cudmg(i),cudmg2(i),xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,era2(i),ede2(i),alfsi2,delsi2,
     ?nstart,nfinal,alsi2,desi2,ktira,ra,de,iuth,iutm,sut,iutano,iutmes,
     ?iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny



      enddo

      close (64)


c
c     Faz a estatistica de alvos, com "O-A" (Observado menos Alvo), para
c     as reducoes com 2MASS corrigido
c
c     2MASS tambem com correcao de movimentos proprios aplicados
c     para as estrelas comuns ao UCAC2, e com aplicacao de m.p.
c     medio para as demais estrelas 2MASS
c

      write (*,*) 
      write (*,*) 
c     write(*,*)'*** 2MASS vs. UCAC2 correction: 1 x 1 field size   ***'
      write(*,*)'*** 2MASS + tangent plan + common & non-common pm  ***'
      write (*,*) 


      call estat (box,ialvos,irme2m,iiivo2)

c
c     Finalizando campo
c

 62   continue



c
c     Stores flattened image in fits format for debug purposes
c

      do i=1,ny
      do j=1,nx
      pixel(j,i)=pixmat(j,i)
      enddo
      enddo
      

      subf='teste.fits'
      if=83
      betpix=-32
      mswap=2
      scale=1.d0
      zero=0.d0


      call wfits (ipmax,if,subf,betpix,mswap,nx,ny,pixel,scale,zero)

c

 60   continue

c


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
      subroutine refits (ipmax,pixmat,infits,nx,ny,nheads,ichobj,ipflag,
     ?bscale,bzero,kswap,iswap,bitpix)

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

      dimension pixmat(ipmax,ipmax)
      integer*2 bitpix

      real*4 pixmat
      integer*2 iwork2
      integer*4 iwork4
      real*4 work4
      integer*8 iwork8
      real*8 work8

      integer*1 swork,iby4


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

c     real*4 pixmat
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

      subroutine  GCC (ipmax,icofsp,pixmat,nx,mx,ny,my,raio,maximo,bx,
     ?by,sigdex,sigdey,icontt,sigx,alturx,fundox)

      IMPLICIT REAL*8 (A-H,O-Z)
      real*4 pixmat(ipmax,ipmax),maximo,xmaxi
      DIMENSION DELTAX(5),XSIGMA(5),PARAM(5)
 
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
   14 CALL GAUSIC (ipmax,icofsp,pixmat,ICONTT,RAIO,NX,NY,MX,MY,maximo,
     ?NTERMX,PARAM,DELTAX,XSIGMA,XLAMDA,RESIDX)
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
      SUBROUTINE GAUSIC (ipmax,icofsp,pixmat,KEY,RAIO,NX,NY,MX,MY,
     ?maximo,NTERMS,A,DELTAA,SIGMAA,FLAMDA,CHISQR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(5),DELTAA(5),SIGMAA(5),B(5),ALPHA(icofsp,icofsp),
     ?BETA(icofsp),DERIV(5),ARRAY(icofsp,icofsp)
      real*4 pixmat(ipmax,ipmax),maximo

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
   63 CHISQ1=QIQUAD(ipmax,pixmat,KEY,RAIO,NX,NY,MX,MY,FREE,A,maximo)
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
   80 CALL MATINV (NTERMS, icofsp, array, DET)
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
   93 CHISQR = QIQUAD(ipmax,pixmat,KEY,RAIO,NX,NY,MX,MY,FREE,B,maximo)
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
      DOUBLE PRECISION FUNCTION QIQUAD (ipmax,pixmat,KEY,RAIO,NX,NY,MX,
     ?MY,FREE,A,maximo)
      IMPLICIT REAL *8 (A-H,O-Z)
      real*4 pixmat(ipmax,ipmax),maximo

      DIMENSION A(5)
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
c     Mod.: M. Assafin 10/Nov/2006
c
c

      subroutine tmass (idiobs,mraiz,rac,dec,drac,ddec,ra2ma,de2ma,
     ?era2ma,ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,nest)

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
c     Mod.: M. Assafin  24/Out/2012
c
c

      subroutine stmass (idiobs,mraiz,rac,dec,drac,ddec,ra2ma,de2ma,
     ?era2ma,ede2ma,dmgj,dmgh,dmgk,emgj,emgh,emgk,ddj2,nest)

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

      decmin=dec-ddec
      decmax=dec+ddec


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



      rfaixa=jaca*0.1d0+0.05d0

      if (rfaixa.gt.90.d0) rfaixa=89.95d0


      xx1=rac-drac/dabs(dcos(grarad*rfaixa))
      xx2=rac+drac/dabs(dcos(grarad*rfaixa))



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


      if (de.lt.decmin) go to 18
      if (de.gt.decmax) go to 18



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
c     Mod. 10/Nov/2006
c
c

      subroutine ucac2 (idiobs,uraiz,epoj,rac,dec,drac,ddec,rauc2,deuc2,
     ?erauc2,edeuc2,pmra,pmde,epmra,epmde,udmgj,udmgh,udmgk,udmg,nest)

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

      xbox=grarad*drac
      ybox=grarad*ddec
      grac=grarad*rac
      gdec=grarad*dec

c
c     Leitura das faixas de declinacao de 0.5 em 0.5 graus
c     do UCAC2
c

      dfaixa=dec-ddec-0.5d0
      decmax=dec+ddec


c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.5d0


      if (dfaixa-decmax.gt.0.5d0) go to 35

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


      ra=idat(1)
      de=idat(2)


c
c     Checa se estrela cai dentro do campo
c

      bra=grarad*ra/3600d3
      bde=grarad*de/3600d3

      d=DEXY(bra,bde,grac,gdec)
      xx=XPAD(bra,bde,grac)/d
      yy=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xx).gt.xbox) go to 20
      if (dabs(yy).gt.ybox) go to 20

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

      ra=ra+pmx*(epoj-2000d0)/10d0
      de=de+pmy*(epoj-2000d0)/10d0

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
c     Mod. M. Assafin 24/Out/2012
c
c

      subroutine sucac2 (idiobs,iu2z,iu2i,inuu2,ir1u2,ir2u2,uraiz,epoj,
     ?rac,dec,drac,ddec,rauc2,deuc2,erauc2,edeuc2,pmra,pmde,epmra,epmde,
     ?udmgj,udmgh,udmgk,udmg,nest)

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

      xbox=grarad*drac
      ybox=grarad*ddec
      grac=grarad*rac
      gdec=grarad*dec


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

      dfaixa=dec-ddec-0.5d0
      decmax=dec+ddec


c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.5d0


      if (dfaixa-decmax.gt.0.5d0) go to 35

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

      if (dfaixa.ge.0.d0) then
      rfaixa=jaca*0.5d0-90.5d0 
      else
      rfaixa=jaca*0.5d0-90.0d0 
      endif

      rfaixa=dabs(rfaixa)

      if (rfaixa.gt.90.d0) rfaixa=90.d0



      xx1=rac-drac/dabs(dcos(grarad*rfaixa))
      xx2=rac+drac/dabs(dcos(grarad*rfaixa))


 
      xx1=xx1/15.d0
      xx2=xx2/15.d0
 
 
      if (xx1.lt.0d0)  xx1=xx1+24.d0
      if (xx1.gt.24d0) xx1=xx1-24.d0

      if (xx2.lt.0d0)  xx2=xx2+24.d0
      if (xx2.gt.24d0) xx2=xx2-24.d0



      if (xx2.ge.xx1) then

      inde=1

      jxmin(1)=xx1/bin+1
      jxmax(1)=xx2/bin+1

      cxmin(1)=xx1*54.d6
      cxmax(1)=xx2*54.d6


      else

      inde=2

      jxmin(1)=xx1/bin+1
      jxmax(1)=nbin

      jxmin(2)=1
      jxmax(2)=xx2/bin+1

      cxmin(1)=xx1*54.d6
      cxmax(1)=24.d0*54.d6
      cxmin(2)=0.d0
      cxmax(2)=xx2*54.d6


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


      de=idat(2)

      pmx=idat(12)
      pmy=idat(13)

      ra=ra+pmx*(epoj-2000d0)/10d0
      de=de+pmy*(epoj-2000d0)/10d0

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
c     Mod. M. Assafin 12/Jan/2013
c
c

      subroutine ucac4 (idiobs,u4raiz,epoj,rac,dec,drac,ddec,rauc4,
     ?deuc4,erauc4,edeuc4,pmra4,pmde4,epmra4,epmde4,udmgj4,udmgh4,
     ?udmgk4,udmg4,nest)

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

      xbox=grarad*drac
      ybox=grarad*ddec
      grac=grarad*rac
      gdec=grarad*dec

c
c     Leitura das faixas de declinacao de 0.2 em 0.2 graus
c     do UCAC4
c

      dfaixa=dec-ddec-0.2d0
      decmax=dec+ddec


c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.2d0


      if (dfaixa-decmax.gt.0.2d0) go to 35

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



      ra=ran
      de=spdn


      rab=ra/3.6d6
      deb=de/3.6d6
      deb=deb-90.d0




c
c     Checa se estrela cai dentro do campo
c

      bra=grarad*rab
      bde=grarad*deb

      d=DEXY(bra,bde,grac,gdec)
      xx=XPAD(bra,bde,grac)/d
      yy=YPAD(bra,bde,grac,gdec)/d

      if (dabs(xx).gt.xbox) go to 20
      if (dabs(yy).gt.ybox) go to 20


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



      ra=ran
      de=spdn

      if (pmra2.eq.32767 .and. pmdc2.eq.32767) then

      do m=1,25 
      if (rnm.eq.irnm(m)) go to 23
      enddo

 23   continue

      pmx=ipmrc(m)/dcos((de/3.6d6-90d0)*grarad)
      pmy=ipmd(m)


      else

      pmx=pmra2/dcos((de/3.6d6-90d0)*grarad)
      pmy=pmdc2

      endif


      ra=ra+pmx*(epoj-2000d0)/10d0
      de=de+pmy*(epoj-2000d0)/10d0

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
c     Mod. M. Assafin 12/Jan/2013
c
c

      subroutine sucac4 (idiobs,iu4z,iu4i,inuu4,ir1u4,ir2u4,u4raiz,epoj,
     ?rac,dec,drac,ddec,rauc4,deuc4,erauc4,edeuc4,pmra4,pmde4,epmra4,
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
c     Leitura das faixas de declinacao de 0.2 em 0.2 graus
c     do UCAC4
c

      dfaixa=dec-ddec-0.2d0
      decmax=dec+ddec


c

      nest=0

      do 30 k=1,izon

      dfaixa=dfaixa+0.2d0


      if (dfaixa-decmax.gt.0.2d0) go to 35

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

      if (dfaixa.ge.0.d0) then
      rfaixa=jaca*0.2d0-90.2d0 
      else
      rfaixa=jaca*0.2d0-90.0d0 
      endif

      rfaixa=dabs(rfaixa)


      if (rfaixa.gt.90.d0) rfaixa=90.d0



      xx1=rac-drac/dabs(dcos(grarad*rfaixa))
      xx2=rac+drac/dabs(dcos(grarad*rfaixa))


      if (xx1.lt.0d0)  xx1=xx1+360.d0
      if (xx1.gt.360.d0) xx1=xx1-360.d0


      if (xx2.lt.0d0)  xx2=xx2+360.d0
      if (xx2.gt.360.d0) xx2=xx2-360.d0


      if (xx2.ge.xx1) then

      inde=1

      jxmin(1)=xx1/bin+1
      jxmax(1)=xx2/bin+1

      cxmin(1)=xx1*3600.d3
      cxmax(1)=xx2*3600.d3


      else


      inde=2

      jxmin(1)=xx1/bin+1
      jxmax(1)=nbin

      jxmin(2)=1
      jxmax(2)=xx2/bin+1

      cxmin(1)=xx1*3600.d3
      cxmax(1)=360.d0*3600.d3
      cxmin(2)=0.d0
      cxmax(2)=xx2*3600.d3


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

      de=spdn

      if (pmra2.eq.32767 .and. pmdc2.eq.32767) then

      do m=1,25 
      if (rnm.eq.irnm(m)) go to 16
      enddo

 16   continue

      pmx=ipmrc(m)/dcos((de/3.6d6-90d0)*grarad)
      pmy=ipmd(m)


      else

      pmx=pmra2/dcos((de/3.6d6-90d0)*grarad)
      pmy=pmdc2

      endif


      ra=ra+pmx*(epoj-2000d0)/10d0
      de=de+pmy*(epoj-2000d0)/10d0

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



      subroutine cuser (idiobs,ifaixa,epoj,rac,dec,drac,ddec,raucs,
     ?deucs,eraucs,edeucs,pmras,pmdes,epmras,epmdes,udmgjs,udmghs,
     ?udmgks,udmgs,nest)



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
c
c     Last Mod.   M. Assafin  29/Out/2012
c

      subroutine idxy2m (idiobs,icofsp,carx,cary,nval,ior,scala,erpix,
     ?rac,dec,nbcat,nbmed,nest,ialtu,xob,yob,n2mass,id2ma,ra2ma,de2ma,
     ?dmgj,xold,yold,xra2ma,yde2ma,ireflex,ecala,tt)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension ialtu(idiobs),xob(idiobs),yob(idiobs),id2ma(idiobs),
     ?ra2ma(idiobs),de2ma(idiobs),dmgj(idiobs),xra2ma(idiobs),
     ?yde2ma(idiobs),nval(idiobs),ior(idiobs),xold(idiobs),yold(idiobs)

      dimension coefx(icofsp),coefy(icofsp),xcof(icofsp),ycof(icofsp),
     ?xest(idiobs),yest(idiobs),xp(idiobs),yp(idiobs)

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

      pcala=grarad*ecala/3600.d0

      epixra=grarad*erpix*scala/3600.d0
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

      if (dabs(xold(i)).gt.carx) go to 5
      if (dabs(yold(i)).gt.cary) go to 5


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
c     Identificacao par a par, com e sem reflexao
c
c     k=1 : sem reflexao em x
c     k=-1: com reflexao em x
c


      do i=1,icofsp
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

      ntempo=ncat*0.1d0

      tempoi=0.d0
      call tempo (tempoi,tempot,tempop)

c



      do iii=1,ncat-1

      if (k.eq.-1 .and. iii.eq.ntempo) then

      call tempo (tempoi,tempot,tempop)

      tt=2.d0*ncat*tempop/ntempo
      tt=tt/3600.d0
      
      ihour=tt
      minu=(tt-ihour)*60.d0
      seg=((tt-ihour)*60.d0-minu)*60.d0

      write (*,*)      
      write (*,1) ihour,minu,seg
 1    format(1x,'Time consuming per field for reference catalog star ide
     ?ntification: ',i3,'hs ',i2,'m ',f4.1,'s')  
      write (*,*)      

      endif

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

      call isol (idiobs,icofsp,ngrau,nptos,xest,yest,xp,yp,coefx,coefy)

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

      x=pol(icofsp,xx,yy,coefx,ngraup)
      y=pol(icofsp,xx,yy,coefy,ngraup)

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

      write (*,*) 
      write (*,*) '(RA,DE) vs. (x,y) cross-identification'
      write (*,*) '1rst step: 4 constant adjust'
      write (*,*) 

      write (*,*) '2MASS stars extracted in 2x2 field  = ',n2mass
      write (*,*) 'Extracted field stars with (x,y)    = ',nest
      write (*,*) 'Common 2MASS vs. bright field stars = ',ncomum
      write (*,*) 'system reflection flag              = ',ireflex
      write (*,*)

c

      x1=xcof(1)*radgra*3600.d0
      x2=xcof(2)*radgra*3600.d0
      x3=xcof(3)*radgra*3600.d0

      y1=ycof(1)*radgra*3600.d0
      y2=ycof(2)*radgra*3600.d0
      y3=ycof(3)*radgra*3600.d0

c
      write (*,*) 'Coeficients solution'
      write (*,*) 'X = A + Bx + Cy '
      write (*,*) 'Y = D - Cx + By '
      write (*,*)
      write (*,*) 'X: ',x1,x2,x3
      write (*,*) 'Y: ',y1,y2,y3
      write (*,*)


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

      x=pol(icofsp,xx,yy,xcof,ngraup)
      y=pol(icofsp,xx,yy,ycof,ngraup)

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



c     x1=dx*radgra*3600.d0
c     y1=dy*radgra*3600.d0

c     write (*,*) 'x y dx dy ',xob(j),yob(j),x1,y1



      ncomum=ncomum+1

      id2ma(j)=ii
c     nval(ncomum)=ii
      xp(ncomum)=xold(ii)
      yp(ncomum)=yold(ii)

      xest(ncomum)=ireflex*xob(j)
      yest(ncomum)=yob(j)


 30   continue



c     write (*,*) 'ncomum final = ',ncomum


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


      call isol (idiobs,icofsp,ngraup,ncomum,xest,yest,xp,yp,coefx,
     ?coefy)

      if (ierro.eq.1) then

      return
      endif

c

      write (*,*) 
      write (*,*) '(RA,DE) vs. (x,y) cross-identification'
      write (*,*) '2nd final step:'
      write (*,*) 'Complete polynomial model, degree = ',ngraup
      write (*,*)  

      write (*,*) '2MASS stars extracted in 2x2 field   = ',n2mass
      write (*,*) 'Extracted field stars with (x,y)     = ',nest
      write (*,*) 'Common 2MASS vs. field stars (final) = ',ncomum
      write (*,*) 'system reflection flag             = ',ireflex
      write (*,*)

c


      x1=coefx(1)*radgra*3600.d0
      x2=coefx(2)*radgra*3600.d0
      x3=coefx(3)*radgra*3600.d0

      y1=coefy(1)*radgra*3600.d0
      y2=coefy(2)*radgra*3600.d0
      y3=coefy(3)*radgra*3600.d0



c
      write (*,*) 'Solution. Linear coeficients:'
      write (*,*) 'X = A + Bx + Cy '
      write (*,*) 'Y = D + Ex + Fy '
      write (*,*)
      write (*,*) 'X: ',x1,x2,x3
      write (*,*) 'Y: ',y1,y2,y3
      write (*,*)


c
c     Calculando alfas e deltas provisorios no sistema 2MASS para
c     todas as estrelas do campo para reconhecimento com outros
c     catalogos (UCAC2, etc)
c

      do i=1,nest


      xx=ireflex*xob(i)
      yy=yob(i)

      x=pol(icofsp,xx,yy,coefx,ngraup)
      y=pol(icofsp,xx,yy,coefy,ngraup)

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
c
c
c     Nesta subrotina, as orientacoes e escalas sao fornecidas apartir da reducao
c     da primeira imagem, tornando a identificacao das demais imagens acelerada.
c
c
c
c      Modificada:  M. Assafin  21/Out/2012
c



      subroutine mdxy2m (idiobs,icofsp,carx,cary,ngrau,coefxr,coefyr,
     ?kreflex,nval,ior,rscala,erpix,rac,dec,nbcat,nbmed,nest,ialtu,xob,
     ?yob,n2mass,id2ma,ra2ma,de2ma,dmgj,xold,yold,xra2ma,yde2ma,ireflex,
     ?pcala,tt)



      IMPLICIT REAL*8 (A-H,O-Z)

      dimension ialtu(idiobs),xob(idiobs),yob(idiobs),id2ma(idiobs),
     ?ra2ma(idiobs),de2ma(idiobs),dmgj(idiobs),xra2ma(idiobs),
     ?yde2ma(idiobs),nval(idiobs),ior(idiobs),xold(idiobs),yold(idiobs)

      dimension coefx(icofsp),coefy(icofsp),coefxr(icofsp),
     ?coefyr(icofsp),xcof(icofsp),ycof(icofsp),xest(idiobs),
     ?yest(idiobs),xp(idiobs),yp(idiobs)

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

      ncof=icofsp
      
      ireflex=kreflex

c

      epixra=erpix*rscala

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
 
      if (dabs(xold(i)).gt.carx) go to 5
      if (dabs(yold(i)).gt.cary) go to 5

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

      cx1=xp(1)-pol(icofsp,xest(1),yest(1),coefxr,ngrau)
      cy1=yp(1)-pol(icofsp,xest(1),yest(1),coefyr,ngrau)
      cx2=xp(2)-pol(icofsp,xest(2),yest(2),coefxr,ngrau)
      cy2=yp(2)-pol(icofsp,xest(2),yest(2),coefyr,ngrau)


      coefxr(1)=(cx1+cx2)/2.d0
      coefyr(1)=(cy1+cy2)/2.d0


c
c     Computa o numero de identificacoes para esses 2 pares
c

      icont=0

      do 10 j=1,nmed

      xx=k*xob(ialtu(nest-j+1))
      yy=yob(ialtu(nest-j+1))

      x=pol(icofsp,xx,yy,coefxr,ngrau)
      y=pol(icofsp,xx,yy,coefyr,ngrau)

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

      call tempo (tempoi,tempot,tempop)

      tt=tempop
      tt=tt/3600.d0
      
      ihour=tt
      minu=(tt-ihour)*60.d0
      seg=((tt-ihour)*60.d0-minu)*60.d0

      write (*,*)      
      write (*,1) ihour,minu,seg
 1    format(1x,'Time consuming per field for reference catalog star ide
     ?ntification: ',i3,'hs ',i2,'m ',f4.1,'s')  
      write (*,*)      




c
c     debug
c

      write (*,*) 
      write (*,*) '(RA,DE) vs. (x,y) cross-identification'
      write (*,*) '1rst step: coefficients from 1rst image solution'
      write (*,*) 

      write (*,*) '2MASS stars extracted in 2x2 field  = ',n2mass
      write (*,*) 'Extracted field stars with (x,y)    = ',nest
      write (*,*) 'Common 2MASS vs. bright field stars = ',ncomum
      write (*,*) 'system reflection flag              = ',kreflex
      write (*,*)

c

      x1=xcof(1)*radgra*3600.d0
      x2=xcof(2)*radgra*3600.d0
      x3=xcof(3)*radgra*3600.d0

      y1=ycof(1)*radgra*3600.d0
      y2=ycof(2)*radgra*3600.d0
      y3=ycof(3)*radgra*3600.d0

c
      write (*,*) 'Solution. Linear coeficients:'
      write (*,*) 'X = A + Bx + Cy '
      write (*,*) 'Y = D + Ex + Fy '
      write (*,*)
      write (*,*) 'X: ',x1,x2,x3
      write (*,*) 'Y: ',y1,y2,y3
      write (*,*)


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

      x=pol(icofsp,xx,yy,xcof,ngrau)
      y=pol(icofsp,xx,yy,ycof,ngrau)

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


      call isol (idiobs,icofsp,ngrau,ncomum,xest,yest,xp,yp,coefx,coefy)

      if (ierro.eq.1) then
      return
      endif

c

      write (*,*) 
      write (*,*) '(RA,DE) vs. (x,y) cross-identification'
      write (*,*) '2nd final step:'
      write (*,*) 'Complete polynomial model, degree = ',ngrau
      write (*,*)  

      write (*,*) '2MASS stars extracted in 2x2 field   = ',n2mass
      write (*,*) 'Extracted field stars with (x,y)     = ',nest
      write (*,*) 'Common 2MASS vs. field stars (final) = ',ncomum
      write (*,*) 'system reflection flag             = ',ireflex
      write (*,*)

c


      x1=coefx(1)*radgra*3600.d0
      x2=coefx(2)*radgra*3600.d0
      x3=coefx(3)*radgra*3600.d0

      y1=coefy(1)*radgra*3600.d0
      y2=coefy(2)*radgra*3600.d0
      y3=coefy(3)*radgra*3600.d0



c
      write (*,*) 'Solution. Linear coeficients:'
      write (*,*) 'X = A + Bx + Cy '
      write (*,*) 'Y = D + Ex + Fy '
      write (*,*)
      write (*,*) 'X: ',x1,x2,x3
      write (*,*) 'Y: ',y1,y2,y3
      write (*,*)


c
c     Calculando alfas e deltas provisorios no sistema 2MASS para
c     todas as estrelas do campo para reconhecimento com outros
c     catalogos (UCAC2, etc)
c

      do i=1,nest


      xx=ireflex*xob(i)
      yy=yob(i)

      x=pol(icofsp,xx,yy,coefx,ngrau)
      y=pol(icofsp,xx,yy,coefy,ngrau)

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
c     Atualizacao: M. Assafin  22/Dez/2005
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


      call solucao (idiobs,icofsp,ngrau,ngrau3,ngrau5,nstart,xest,yest,
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
c     Subroutine wfits
c
c     Writes a fits image given a matriz of pixels.
c
c
c     bitpix = 16  integer*2 data
c              +32 integer*4 data
c              -32 real*4 data
c              +64 integer*8 data
c              -64 real*8 data
c
c
c     iswap = 1 (do not swap image)
c             2 (swap image)
c
c
c
c
c     Last modified:   M. Assafin  18/Dec/2009
c

      subroutine wfits (ipmax,if,file,bitpix,iswap,nx,ny,pixel,bscale,
     ?bzero)

      implicit real*8 (a-h,o-z)


      integer*2 iwork2(1440)
      integer*4 iwork4(720)
      integer*8 iwork8(360)
      real*4    rwork4(720)
      real*8    rwork8(360)

      integer*1 swork(2880),iby4(4)
      integer*2 bitpix

      real*4 pixel(ipmax,ipmax)

      character*50 file
      character*2880 header

      icab(il,ic)=(il-1)*80+ic


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


      open (if,file=file,access='direct',form='unformatted',
     ?recl=2880)



c
c     Writes fits header
c

      header=''

      l=1
      ic1=1
      ic2=46
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='SIMPLE  =                    T / Fits Standard'

      l=2
      ic1=1
      ic2=47
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='BITPIX  =                      / Bits per pixel'

      l=3
      ic1=1
      ic2=47
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='NAXIS   =                    2 / Number of axes'


      l=4
      ic1=1
      ic2=44
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='NAXIS1  =                      / Axis Length'


      l=5
      ic1=1
      ic2=44
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='NAXIS2  =                      / Axis Length'


      l=6
      ic1=1
      ic2=44
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='BSCALE  =                      / Data scale '


      l=7
      ic1=1
      ic2=44
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='BZERO   =                      / Zero point '


      l=2
      ic1=28
      ic2=30
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      write (header(ip1:ip2),'(i3)') bitpix


      l=4
      ic1=26
      ic2=30
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      write (header(ip1:ip2),'(i5)') nx

      l=5
      ic1=26
      ic2=30
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      write (header(ip1:ip2),'(i5)') ny


      l=6
      ic1=11 
      ic2=30
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      write (header(ip1:ip2),'(f20.10)') bscale


      l=7
      ic1=11 
      ic2=30
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      write (header(ip1:ip2),'(f20.10)') bzero 



      do l=8,35
      ic1=1
      ic2=41
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='COMMENTS=                      / Comments'
      enddo

      l=36
      ic1=1
      ic2=3
      ip1=icab(l,ic1)
      ip2=icab(l,ic2)
      header(ip1:ip2)='END'

      irec=1

      write (if,rec=irec) header

c
c     Now writes the data
c


      m=0

      do i=1,ny
      do j=1,nx

      m=m+1

c
      if (bitpix.eq.16) then
      iwork2(m)=pixel(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) iwork2
      endif
      endif
c
      if (bitpix.eq.32) then
      iwork4(m)=pixel(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) iwork4
      endif
      endif
c
      if (bitpix.eq.64) then
      iwork8(m)=pixel(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) iwork8
      endif
      endif
c
      if (bitpix.eq.-32) then
      rwork4(m)=pixel(j,i)
      if (m.eq.kwork) then
      irec=irec+1
      write (if,rec=irec) rwork4
      endif
      endif
c
      if (bitpix.eq.-64) then
      rwork8(m)=pixel(j,i)
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
c     Subroutine flatsk
c
c
c     Flattens sky background
c
c
c     last update:  M. Assafin 12/Mar/2010
c
c


      subroutine flatsk (ipmax,pixmat,nx,ny,ngrauf,mazi,fcmin,fcmax)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi

c     dimension xesto(25010001),yesto(25010001),po(25010001),cofo(136)
      dimension xesto(ipmax*ipmax),yesto(ipmax*ipmax),po(ipmax*ipmax),
     ?cofo(136)
      dimension pixmat(ipmax,ipmax)



c


      idimo=ipmax*ipmax
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

      call isolo (ipmax,ngrauf,n,xesto,yesto,po,cofo)

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


      subroutine skyb (ipmax,kfund,pixmat,nx,ny,mazi,malisa,vmin,fatceu,
     ?ceu,sigceu,izmax,threso)

      implicit real*8 (a-h,o-z)

      real*4 pixmat,mazi

      dimension ifundo(kfund),jfundo(kfund)
      dimension pixmat(ipmax,ipmax)


c
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




c
c
c     Subroutine mesure
c
c
c     Purpose
c
c
c     User furnishes regions for measurement of specific targets, for each
c     individual image.
c
c    
c
c     Comments
c
c
c     There are two types of targets. The two types are identified by
c     the type of reagion: circle or box.
c
c     In the case of circles, targets have stellar-like PSFs, and are fitted
c     by the circular Gaussian model. 
c
c     In the case of boxes, the targets are trace images, and are fitted
c     by the ERF model.
c    
c
c     These user's regions have higher priority over the automatic object
c     search. This means that no automatic search is applied over the marked
c     regions. The fits are made directly on the pixels inside the furnished
c     regions.
c
c     The format of the regions follow the ds9 package standards. The user
c     should use the cursor facilities of ds9 to produce the regions. 
c
c
c     Last update: M. Assafin   09 Dec 2012
c
c




      subroutine mesure (idiobs,imes,nmcir,xcir,ycir,lacir,nmtra,xtra,
     ?ytra,xlatra,ylatra,angtra)

      implicit real*8 (a-h,o-z)
      parameter(idim=150)

      dimension xcir(idiobs),ycir(idiobs),lacir(idiobs)

      dimension xtra(idiobs),ytra(idiobs),xlatra(idiobs),ylatra(idiobs),
     ?angtra(idiobs)


      character*150 imes,linha

c

      open (25,file=imes)

c
c     Searches for targets
c

      nmcir=0
      nmtra=0

      do k=1,idiobs

      linha=''
      read (25,5,end=10) linha
 5    format(a150)

c
c     Searches for circle targets
c

      if (linha(1:6).eq.'circle') then
      nmcir=nmcir+1
      do i=1,idim
      if (linha(i:i).ne.'.') then
      icomp=ichar(linha(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) linha(i:i)=' '
      endif
      enddo
      read (linha,*) xcir(nmcir),ycir(nmcir),raio
      lacir(nmcir)=raio
      endif

c
c     Searches for trace-image targets
c

      if (linha(1:3).eq.'box') then
      nmtra=nmtra+1
      do i=1,idim
      if (linha(i:i).ne.'.') then
      icomp=ichar(linha(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) linha(i:i)=' '
      endif
      enddo
      read (linha,*) xtra(nmtra),ytra(nmtra),dx,dy,ang
      xlatra(nmtra)=dx
      ylatra(nmtra)=dy
      angtra(nmtra)=ang
      endif


      enddo

c

 10   close (25)




      return

      end





c
c     Subroutine irot
c
c
c
c     Purpose
c
c
c     Transforms (x,y) to (x',y') by a rotation around a center (xc,yc)
c
c
c
c     Last update: M. Assafin   21 Jan 2013
c

      subroutine  irot (jj,ii,ixc,iyc,ang,nx,ny,k,m)

      implicit real*8 (a-h,o-z)

c

c     j=jj-ixc
c     i=ii-iyc


      k=jj*dcos(ang)+ii*dsin(ang)
      m=-jj*dsin(ang)+ii*dcos(ang)

      k=k+ixc
      m=m+iyc

      if (k.lt.1) k=1
      if (m.lt.1) m=1

      if (k.gt.nx) k=nx
      if (m.gt.ny) m=ny

      return
      end



c
c     Subroutine trace 
c
c
c     Purpose
c
c
c     Fits an image in the form of a trace, and finds the (x,y) centroid, among other
c     parameters, and their errors.
c
c
c     Comments
c
c     The 7 parameters of the ERF model: 
c
c
c    1 - height (counts)
c    2 - x0 (pixels)
c    3 - y0 (pixels)
c    4 - sigma (related to seeing) in pixels
c    5 - tetha angle (in radians)
c    6 - trace lenght (pixels)
c    7 - background (counts)
c
c
c
c     Last update: M. Assafin   21 Jan 2013
c
c
c
c


      subroutine trace (idiobs,ipmax,icofsp,pixmat,nx,ny,xc,yc,rx,ry,
     ?ang,maximo,seeing,ceu,dlen,icomax,altura,bx,by,sigma,tetha,
     ?dlenght,fundo,ex,ey)

      implicit real*8 (a-h,o-z)
      parameter(nterms=7)

      real*4 pixmat(ipmax,ipmax),maximo,xmaxi

      dimension nval(idiobs),ior(idiobs),histo(ipmax*4)

      dimension deltax(nterms),xsigma(nterms),param(nterms)

      common /a14/ierro



c
c     Initial data
c

      pi    = 0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi


c
      xmaxi=maximo


c
c     Parameter increment and lambda factor in Maquard non-linear LS method
c

      dinc=0.05d0
      xlamda=0.001
c     zlamda=1.d-10

c
c     Convergence limits
c

      icomax=2 
      icontx=0

      xresid=-1.d14
      xcent=1.d14
      ycent=1.d14

      dlimit=0.01d0
      plimit=0.05d0

c
c     Trace area and perimeter parameters
c

      ixc=xc
      iyc=yc
      irx=rx/2.d0
      iry=ry/2.d0

      rang=ang*grarad

c
c     Estimating initial values for the parameters
c

c
c     Background
c


      n=0

   
      ii=-iry
      do jj=-irx,+irx
      call irot (jj,ii,ixc,iyc,ang,nx,ny,k,m)
      n=n+1
      histo(n)=pixmat(k,m)
      enddo

   
      ii=+iry
      do jj=-irx,+irx
      call irot (jj,ii,ixc,iyc,ang,nx,ny,k,m)
      n=n+1
      histo(n)=pixmat(k,m)
      enddo


   
      jj=-irx
      do ii=-iry+1,+iry-1
      call irot (jj,ii,ixc,iyc,ang,nx,ny,k,m)
      n=n+1
      histo(n)=pixmat(k,m)
      enddo


   
      jj=+irx
      do ii=-iry+1,+iry-1
      call irot (jj,ii,ixc,iyc,ang,nx,ny,k,m)
      n=n+1
      histo(n)=pixmat(k,m)
      enddo


 
      if (n.gt.idiobs) n=idiobs

      do i=1,n
      ior(i)=i
      nval(i)=histo(i)*1000
      enddo


      call ordem (idiobs,n,ior,nval)
 
      n1=n/4
      n2=3*n/4

      n=0
      xm=0.d0
      x2=0.d0
      do k=n1,n2
      i=ior(k)
      xm=xm+histo(i)
      x2=x2+histo(i)**2
      n=n+1
      enddo

      call desvio (n,xm,x2)

      param(7)=xm

      param(7)=ceu

c
c     (x, y) center
c



      param(2)=xc
      param(3)=yc



c
c     Height
c


      h=0.d0
      do ky=-1,1
      do kx=-1,1
      if (pixmat(ixc+kx,iyc+ky).gt.h) h=pixmat(ixc+kx,iyc+ky)
      enddo
      enddo


      param(1)=h



c
c     Sigma
c
c     It comes from previous, external determinations of sigma by the
c     Gaussian fits for the stellar-like objects); the sigma is stored
c     in the variable "seeing", but it is not the FWHM.
c


      param(4)=seeing


c
c     Tetha
c

      c=0.d0
      axx=0.d0
      axy=0.d0
      ayy=0.d0
 
c     nn=0
 
      do ii=-iry,+iry
      do jj=-irx,+irx
      call irot (jj,ii,ixc,iyc,ang,nx,ny,j,i)
      if (pixmat(j,i).lt.maximo) then
      axx=axx+pixmat(j,i)*j*j
      axy=axy+pixmat(j,i)*j*i
      ayy=ayy+pixmat(j,i)*i*i
      c=c+pixmat(j,i)
c     nn=nn+1
      endif
      enddo
      enddo
 
      axx=axx/c
      axy=axx/c
      ayy=axx/c
 
c     ax=axx-ayy
c     ay=axy

      ay=axx-ayy
      ax=axy




      aux=dabs(datan2(dabs(ay),dabs(ax)))

      tet=aux

      if (ax.ge.0.d0) then
       if (ay.ge.0.d0) then
        tet=aux
       else
        tet=2.d0*pi-aux
       endif
      else
       if (ay.ge.0.d0) then
        tet=pi-aux
       else
        tet=pi+aux
       endif
      endif

      tetha=0.5d0*tet



c     tetha=0.5d0*datan2(axy,axx-ayy)
 


      tetha=rang

      param(5)=tetha




c
c     Lenght
c

c     dlenght=nn/(3.d0*param(4))


c     dlenght=rx


      dlenght=dlen


      param(6)=dlenght

 


c
c     Initial increments for the parameters
c

      do i=1,nterms
      deltax(i)=dinc*param(i)
      enddo

c
c     Initializing error parameters            
c

      do i=1,nterms
      xsigma(i)=0.d0
      enddo

      
c
c     Trace fit
c



 140  call traco (ipmax,icofsp,pixmat,nx,ny,maximo,nterms,param,
     ?deltax,xsigma,xlamda,ixc,iyc,irx,iry,rang,residx)




c
c     Solution convergence
c

      residx = dsqrt(residx)
      centx  = param(2) 
      centy  = param(3) 
      conver = dabs(residx*dlimit)
      diferd = dabs(residx-xresid)
      difpox = dabs(centx-xcent)
      difpoy = dabs(centy-ycent)

c     if (xlamda.lt.zlamda) go to 150

      if ((diferd.lt.conver).and.(difpox.lt.plimit).and.(difpoy.lt.
     ?plimit)) go to 150

      xresid = residx
      xcent=centx
      ycent=centy
      icontx = icontx + 1
      if (icontx.eq.icomax) go to 150

      go to 140


 150  continue


c
c     Storing parameters for the main program
c

      residx = residx
      altura = param(1)
      bx     = param(2) 
      by     = param(3) 
      sigma  = param(4)
      tetha  = param(5)
      dlenght= param(6) 
      fundo  = param(7)


      ex=residx*xsigma(2)
      ey=residx*xsigma(3)



      maximo=xmaxi


c
c     Escreve imagem em arquivo ASC
c

      ang=tetha

      open (80,file='imagem.real')
      open (81,file='imagem.traco')
      open (82,file='imagem.dif')


      do  ii=-iry,+iry
      do  jj=-irx,+irx

      call irot (jj,ii,ixc,iyc,ang,nx,ny,j,i)


      tt=track(nterms,j,i,param)

      dd=pixmat(j,i)-tt

      write (80,500) j,i,pixmat(j,i)
      write (81,500) j,i,tt
      write (82,500) j,i,dd
 500  format(2(1x,i3.3),1x,f4.0)

 
      enddo
      enddo
 
      close (80)
      close (81)
      close (82)


      return
      end






c
c     Subroutine traco
c
c
c     Purpose
c
c
c     The trace fit itself by the ERF model, following the Marquardt
c     method of non-linear LS.
c
c
c
c     Comments
c
c
c     The 7 parameters of the ERF model: 
c
c
c    1 - height (counts)
c    2 - x0 (pixels)
c    3 - y0 (pixels)
c    4 - sigma (related to seeing) in pixels
c    5 - tetha angle (in radians)
c    6 - trace width (pixels)
c    7 - background (counts)
c
c
c
c    Subroutines required;
c
c
c    track
c    tqquad
c    tderiv
c    matinv
c
c
c
c     Last update: M. Assafin   21 Jan 2013
c
c
c
c

      subroutine traco (ipmax,icofsp,pixmat,nx,ny,maximo,nterms,a,
     ?deltaa,sigmaa,flamda,ixc,iyc,irx,iry,ang,chisqr)


      implicit real *8 (a-h,o-z)

      dimension a(nterms),deltaa(nterms),sigmaa(nterms),b(nterms),
     ?alpha(nterms,nterms),beta(nterms),deriv(nterms),
     ?array(icofsp,icofsp)

      real*4 pixmat(ipmax,ipmax),maximo

      common /a14/ierro

c
c     Initial values
c

      ierro=0
      det=1.d0
      icont=0
      iconv=20
      jconv=0

c

c
c        Evaluate alpha and beta matrices
c

      do 134 j=1, nterms
      beta(j) = 0.d0
      do 134 k=1, j
  134 alpha(j,k) = 0.d0

      do 150  ii=-iry,+iry
      do 5001 jj=-irx,+irx

      call irot (jj,ii,ixc,iyc,ang,nx,ny,l,i)

      if (pixmat(l,i).ge.maximo) go to 5001
      call tderiv (nterms, l, i, a, deltaa, deriv)
      if (ierro.eq.1) go to 107
      icont=icont+1

      do 146 j=1,nterms
      beta(j)= beta(j) + (pixmat(l,i)-track(nterms,l,i,a))*deriv(j)
      if (ierro.eq.1) go to 107
      do 146 k=1, j
  146 alpha(j,k) = alpha(j,k) + deriv(j)*deriv(k)

 5001 continue
  150 continue



c
c     Number of degrees of freedom
c

      free=icont-nterms

c

      if (free.le.0.d0) go to 107


c

      do 53 j=1,nterms
      do 53 k=1, j
   53 alpha(k,j)=alpha(j,k)


c
c        Evaluates chi square at starting point
c

      chisq1=tqquad(ipmax,nterms,pixmat,nx,ny,ixc,iyc,irx,iry,ang,free,
     ?a,maximo)


      if (ierro.eq.1) go to 107




c 				 
c         Invert modified curvature matrix to find the new parameters
c


 71   do 74 j=1, nterms
      DO 73 k=1, nterms
      aux = alpha(j,j)*alpha(k,k)
   73 array(j,k)= alpha(j,k) / dsqrt (aux)
   74 array(j,j) = 1.d0 + flamda
   80 call matinv (nterms, icofsp, array, det)
c

      if (ierro.eq.1) go to 107


C
      do 84 j=1, nterms
      b(j) = a(j)
      do 84 k=1, nterms
      aux = alpha(j,j)*alpha(k,k)
   84 b(j) = b(j) + beta(k)*array(j,k)/dsqrt(aux)



c
c        If chi square increased, increase flamda and try again
c

      chisqr=tqquad(ipmax,nterms,pixmat,nx,ny,ixc,iyc,irx,iry,ang,free,
     ?b,maximo)


      if (ierro.eq.1) go to 107

c
c     Convergence to minimum is not being reached
c

      jconv=jconv+1

      if (jconv.gt.iconv) go to 107


c



      if (chisq1 - chisqr) 95, 101, 101

   95 flamda = 10.d0*flamda

      go to 71


c
c        Evaluate parameters and uncertainties
c

  101 do 104 j=1, nterms
      a(j) = b(j)

      aux = array(j,j)/alpha(j,j)
  104 sigmaa(j) = dsqrt(aux)

      flamda = flamda/10.d0

      go to 110

  107 chisqr = -1.d0
      ierro=1

  110 continue



      return
      end


c
c
c     Subroutine tderiv
c
c
c     Evaluates derivatives of the ERF model function for non-linear LS.
c
c
c     The 7 parameters of the ERF model: 
c
c
c    1 - height (counts)
c    2 - x0 (pixels)
c    3 - y0 (pixels)
c    4 - sigma (related to seeing) in pixels
c    5 - tetha angle (in radians)
c    6 - trace width (pixels)
c    7 - background (counts)
c
c
c
c     Last update: M. Assafin   09 Dec 2012
c
c
c
c

      subroutine tderiv (nterms, j, i, a, deltaa, deriv)

      implicit real*8 (a-h,o-z)

      dimension a(nterms), deltaa(nterms), deriv(nterms)

      common /a14/ierro



c

      pi=0.3141592653589793d1

      x=j
      y=i



c
c     Auxiliary terms
c

      s=dsin(a(5))
      c=dcos(a(5))

      z=(-(x-a(2))*s+(y-a(3))*c)/(a(4)*dsqrt(2.d0))
      z2=z*z
 

      w1=((x-a(2))*c+(y-a(3))*s+(a(6)/2.d0))/(a(4)*dsqrt(2.d0))
      w2=((x-a(2))*c+(y-a(3))*s-(a(6)/2.d0))/(a(4)*dsqrt(2.d0))

      w12=w1*w1
      w22=w2*w2

      dk=((x-a(2))*c+(y-a(3))*s)/(a(4)*dsqrt(2.d0))

      efw1=erf(w1)
      efw2=erf(w2)

      dew12=dexp(-w12)
      dew22=dexp(-w22)

      dez2=dexp(-z2)

      sqpi=dsqrt(pi)
      sqpi2=dsqrt(pi/2.d0)


c
c    The derivatives of the ERF model
c


      deriv(1)=sqpi2*a(4)*dez2*(efw1-efw2)

      deriv(2)=-a(1)*dez2*(sqpi*z*s*(efw1-efw2)+c*(dew12-dew22))

      deriv(3)=-a(1)*dez2*(-sqpi*z*c*(efw1-efw2)+s*(dew12-dew22))

      deriv(4)=a(1)*sqpi2*dez2*((1.d0+2.d0*z2)*(efw1-efw2)+(2.d0/sqpi)*
     ?(w1*dew12-w2*dew22))

      deriv(5)=a(1)*2.d0*sqpi2*a(4)*z*dez2*(dk*(efw1-efw2)+(2.d0/sqpi)*
     ?(dew12-dew22))

      deriv(6)=-0.5d0*a(1)*dez2*(dew12+dew22)

      deriv(7)=1.d0


      return

      end




c
c
c     Function  track
c
c
c     The ERF model function for fitting trace images.
c
c
c     The 7 parameters of the ERF model: 
c
c
c    1 - height (counts)
c    2 - x0 (pixels)
c    3 - y0 (pixels)
c    4 - sigma (related to seeing) in pixels
c    5 - tetha angle (in radians)
c    6 - trace width (pixels)
c    7 - background (counts)
c
c
c
c     Last update: M. Assafin   09 Dec 2012
c
c
c
c

      double precision function track (nterms,j,i,a)

      implicit real*8 (a-h,o-z)

      dimension a(nterms)

      common /a14/ierro

c

      if (a(4).le.0.d0) then
      ierro=1
      return
      endif 

      ierro=0

c

      pi=0.3141592653589793d1

c



      x=j
      y=i



c
c     Auxiliary terms
c

      s=dsin(a(5))
      c=dcos(a(5))

      z=(-(x-a(2))*s+(y-a(3))*c)/(a(4)*dsqrt(2.d0))
      z2=z*z
      dez2=dexp(-z2)
 

      w1=((x-a(2))*c+(y-a(3))*s+(a(6)/2.d0))/(a(4)*dsqrt(2.d0))
      w2=((x-a(2))*c+(y-a(3))*s-(a(6)/2.d0))/(a(4)*dsqrt(2.d0))


      efw1=erf(w1)
      efw2=erf(w2)


      sqpi2=dsqrt(pi/2.d0)


      track=a(1)*sqpi2*a(4)*dez2*(efw1-efw2)+a(7)

c     write (*,*) 'j,i = ',j,i
c     write (*,*) 'z,w1,w2 = ',z,w1,w2
c     write (*,*) 'sqpi2,dez2,efw1,efw2 = ',sqpi2,dez2,efw1,efw2
c     write (*,*) 'track = ',track
c     stop




      return
      end






c
c
c     Function tqquad
c
c
c     Evaluate the reduced chi square for fit to the data
c
c
c
c     Last update: M. Assafin   21 Jan 2013
c
c
c
c

      double precision function tqquad (ipmax,nterms,pixmat,nx,ny,ixc,
     ?iyc,irx,iry,ang,free,a,maximo)

      implicit real*8 (a-h,o-z)

      real*4 pixmat(ipmax,ipmax),maximo

      dimension a(nterms)

      common /a14/ierro


c
      chisq = 0.d0
      tqquad=0.d0
c
      if (free.le.0.d0) then
      ierro=1
      return
      endif
c


      do 39 ii=-iry,+iry
      do 38 jj=-irx,+irx

      call irot (jj,ii,ixc,iyc,ang,nx,ny,j,i)

      if (pixmat(j,i).ge.maximo) go to 38
      chisq = chisq + (pixmat(j,i)-track(nterms,j,i,a))**2
      if (ierro.eq.1) return
 38   continue
 39   continue


c
c     Reduced chi square: divide by the number of degrees of freedom
c

      tqquad = chisq / free

      return
      end



c
c     Error function 
c


      double precision function erf(x)
      implicit real*8 (a-h,o-z)
c     real*8 erf,x
c     real*8 gammp
      if(x.lt.0.d0)then
        erf=-gammp(.5d0,x**2)
      else
        erf=gammp(.5d0,x**2)
      endif
      return
      end



      double precision function gammp(a,x)
      implicit real*8 (a-h,o-z)
c     real*8 a,gammp,x
c     real*8 gammcf,gamser,gln
c     if(x.lt.0.d0.or.a.le.0.d0) stop 'bad arguments in gammp'
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      endif
      return
      end




      subroutine gcf(gammcf,a,x,gln)
      implicit real*8 (a-h,o-z)
c     integer*4 ITMAX
c     real*8 a,gammcf,gln,x,EPS,FPMIN
      parameter (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
c     integer*4 i
c     real*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(dabs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(dabs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(dabs(del-1.d0).lt.EPS)goto 1
11    continue
c     stop 'a too large, ITMAX too small in gcf'
1     gammcf=dexp(-x+a*dlog(x)-gln)*h
      return
      end





      subroutine gser(gamser,a,x,gln)
      implicit real*8 (a-h,o-z)
c     integer*4 ITMAX
c     real*8 a,gamser,gln,x,EPS
      parameter (ITMAX=100,EPS=3.d-7)
c     integer*4 n
c     real*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0)stop 'x < 0 in gser'
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*EPS)goto 1
11    continue
c     stop 'a too large, ITMAX too small in gser'
1     gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      end



      double precision function gammln(xx)
      implicit real*8 (a-h,o-z)
c     real*8 gammln,xx
c     integer*4 j
c     double precision ser,stp,tmp,x,y,cof(6)
      dimension cof(6)
      SAVE cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+dlog(stp*ser/x)
      return
      end






c
c     Subroutine ident
c
c
c     Identify the objects in the FOV above a certain sky background threshold.
c
c
c     Identifies circular-shaped objects in the FOV.
c
c
c
c     Last modified: 10 Jan 2010
c
c
c


      subroutine ident (ipmax,idiobs,nhist,pixmat,imagem,nx,ny,threso,
     ?mazi,xid,yid,idx1,idx2,idy1,idy2,idlado,npix,nest)

      implicit real *8 (a-h,o-z)

      integer*4 imagem(ipmax,ipmax)
      real*4 pixmat(ipmax,ipmax),mazi

      dimension xid(idiobs),yid(idiobs),idlado(idiobs),idx1(idiobs),
     ?idx2(idiobs),idy1(idiobs),idy2(idiobs),npix(idiobs)

      dimension iflag(idiobs),contag(idiobs),histo(nhist),ico(nhist)


c

      do i=1,idiobs
      iflag(i)=0
      enddo

c


      iob=0
      nstar=0

c

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
c     Sets the perimeter of the candidate
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
c     Computes the sky background histogram of the perimeter around
c     the bright central pixel, in the search for the perimeter of
c     the candidate object
c
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
c     Refines the histogram bins
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
c     Picks up the highest frequency and determines the count value
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
c     Double-checks if the central pixel was indeed the brightest pixel
c     within the perimeter. If so, recentres the perimeter
c


      valor=pixmat(ibx,iby)

      do 38 ii=iy1,iy2
      do 37 jj=ix1,ix2
      if (imagem(jj,ii).lt.0) go to 37
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
c     Refines the (x,y) object center, by computing the baricenter
c


      xc=0.d0
      yc=0.d0
      cont=0.d0
      npixel=0

      do ii=iy1,iy2
      do 41 jj=ix1,ix2
      if (pixmat(jj,ii).ge.mazi) go to 41
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
c     Eliminates multiple objects.
c     The one with the largest size is preserved.
c



      do 43  nn=1,nstar
 
      if (lado.le.idlado(nn)) then
 
      if (xc.lt.idx1(nn)) go to 43
      if (yc.lt.idy1(nn)) go to 43
      if (xc.gt.idx2(nn)) go to 43
      if (yc.gt.idy2(nn)) go to 43
 
      go to 45
 
      else
 
      if (xid(nn).lt.ix1) go to 43
      if (yid(nn).lt.iy1) go to 43
      if (xid(nn).gt.ix2) go to 43
      if (yid(nn).gt.iy2) go to 43
 
      iflag(nn)=1
 
 
      endif
 
 
 43   continue


      nstar=nstar+1

      xid(nstar)=xc
      yid(nstar)=yc
      idlado(nstar)=lado

      idx1(nstar)=ix1
      idx2(nstar)=ix2
      idy1(nstar)=iy1
      idy2(nstar)=iy2

      npix(nstar)=npixel


c

 45   continue
 50   continue


c



      nest=0

      do 57 k=1,nstar

      if (iflag(k).gt.0) go to 57

      nest=nest+1

      xid(nest)=xid(k)
      yid(nest)=yid(k)
      idlado(nest)=idlado(k)

      idx1(nest)=idx1(k)
      idx2(nest)=idx2(k)
      idy1(nest)=idy1(k)
      idy2(nest)=idy2(k)

      npix(nest)=npix(k)


 57   continue


c

      return

      end














