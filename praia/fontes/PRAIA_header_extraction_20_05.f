c
c     Programa PRAIA_header_extraction
c
c     PROPOSITO
c
c
c     Dado um conjunto de imagens fits, escreve as coordemadas RA,DEC dadas
c     no header, bem como o instante de observacao, de cada imagem fits.
c     A saida pode ser usada para consertar os dados de RA, DEC e tempo, e
c     e' usada no programa de reducao automatica total de imagem e posicao.
C
C     Atencao: imagens fits lidas em real *4
C
C
C     Ultima modificacao: M. Assafin - 02/Jun/2011
C
C



      IMPLICIT REAL *8 (A-H,O-Z)


      character*50 lista1,lista2
      character*150 infits
      character*20 ichobj,ichfil
      character*1  isig,iver,ivo

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

      data ivo/'|'/


c
      PI    = 0.3141592653589793D1
c     PI=3.141592653589793238462643D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

c
c     Lendo dados de entrada
c

      icam=14

c     open (21,file='PRAIA_header_extraction_20_04.dat')

c

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'PRAIA - Header extraction of fits files')
      write (*,*) ' '
c

c

 2    continue

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '


      read (*,3,err=100,end=100) lista1,iver
 3    format(a50,a1)

      if (iver.ne.ivo) go to 100


      read (*,5,err=100,end=100) lista2
      read (*,*,err=100,end=100) khead
      read (*,*,end=100) iver

c     close (21)


      if (khead.lt.1.or.khead.gt.icam) then
      write (*,*)
      write (*,*) lista1
      write (*,*) 'Invalid telescope site number. Edit input file.'
      go to 2
      endif
c
      write (*,4)
 4    format (15x,'List of fits images to proccess                      
     ?               -> ',$)
      write(*,5) lista1
 5    format(a50)
c
      write (*,6)
 6    format (15x,'List of extracted fits fields                        
     ?               -> ',$)
      write(*,5) lista2

c
      write (*,7)
 7    format (15x,'1 = LNA, HP ; 2 = ESO, SOAR (generic) ; 3 = ESO/WFI 4
     ? = SOAR/SOI ; 5 = SBIG STL ; 6 = CFHT/MegaCam ; 7 = Tubitak 8 = LC
     ?-40 ; 9 = S800/LNA ; 10 = IKON/LNA; 11 = Merlin/Raptor; 12,14 = Pi
     ?c du Midi T1M; 13 = PRECAM;/DES -> ',$)
      write(*,*) khead

c

      write (*,*)
      write (*,*)

c
c     Lendo as imagens
c

      open (3,file=lista1)
      open (7,file=lista2)

c

      i=0
 10   read (3,*,err=20,end=20)
      i=i+1
      go to 10
 20   rewind (3)
      nfiles=i

      do i=1,nfiles
      read (3,33) infits
 33   format(a150)


      write (*,50) i,nfiles
 50   format(1x,'Header extraction: image ',i5,' of ',i5)

c
c     Extraindo o header das imagens
c

      ofra=0.d0
      ofde=0.d0

      if (khead.eq.1.or.khead.eq.5.or.khead.eq.6.or.khead.eq.7.or.khead.
     ?eq.8.or.khead.eq.9.or.khead.eq.10.or.khead.eq.11.or.khead.eq.12.or
     ?.khead.eq.13.or.khead.eq.14) then

      call obhead (infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,iexps,nx,ny,
     ?khead)

      ofra=0.d0
      ofde=0.d0

      go to 59

      endif

      if (khead.eq.2 .or. khead.eq.3 .or. khead.eq.4) then


      call obhaad (infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,iexps,nx,ny,ofra,
     ?ofde,khead)

      endif

c
c     Offsets em alfa e delta de acordo com o Tel/instrumento/CCD
c
c     ofra, ofde em graus
c

      if (khead.eq.2) go to 59

c

      if (iah.eq.99 .or. idg.eq.99) go to 59

c

      ra=hmsgms(iah,iam,sa)

      de=hmsgms(idg,idm,ds)
      if (isig.eq.'-') de=-de


      ra=ra+(ofra/dcos(grarad*dabs(de)))/15.d0

      de=de+ofde

      de=dabs(de)

      iah=ra
      iam=(ra-iah)*60.d0
      sa=((ra-iah)*60.d0-iam)*60.d0

      idg=de
      idm=(de-idg)*60.d0
      ds=((de-idg)*60.d0-idm)*60.d0

c


 59   continue

      write (7,60) iah,iam,sa,isig,idg,idm,ds,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,djm,dj,iexps,ichfil,infits,ichobj,nx,ny
 60   format(1x,i2.2,1x,i2.2,1x,f7.4,1x,a1,i2.2,1x,i2.2,1x,f6.3,2x,i2.2,
     ?1x,i2.2,1x,f5.2,1x,i4.4,1x,i2.2,1x,i2.2,f16.8,1x,f16.8,2x,i4.4,2x,
     ?a20,2x,a50,1x,a20,2(1x,i5.5))

      enddo

      close (3)
      close (7)

c

      go to 2


c
c     End of proccess
c


 100  continue


      write (*,*)
      write (*,*)
      write (*,*) ' Proccess terminated. Status ok.'
      write (*,*)
      write (*,*)


      end



c
c     Subrotina OBHEAD
c
C     Extrai dados do header de imagens fits
C

      subroutine obhead (infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,iexps,nx,ny,
     ?khead)

      IMPLICIT REAL *8 (A-H,O-Z)

      dimension header(2880),digit(4)
      character*3  tend
      character*1  digit
      character*1  header
      character*150 infits
      character*9  ihdat,ihut,ihuts,ihutt,ihcra,ihcar,ihcde,ihobj,ihfil,
     ?ihexp,iquem,ihutb,ihutc,ihfilc,ihexpb,ihexpc
      character*20 ichobj,ichfil,ichexp
      character*1  ler(70),iplic,ibrac,idois,isig,mais,menos
      character*17 ihfilb,iqual
      character*5  ifa

      character*9  imccet,imcceb,imccee,imccer,imcced,imccef,is800

      character*20 imaux,kmaux
      character*4 jmaux
      character*9 sista
      character*29 systa1,systa2

      character*80 itrima

      logical      debug

      data iplic/"'"/
      data idois/':'/
      data ibrac/' '/
      data mais /'+'/
      data menos/'-'/

      data ihobj/'OBJECT  ='/
      data ihdat/'DATE-OBS='/
      data ihut /'UT      ='/
      data ihuts/'UTC     ='/
      data ihutt/'TIME-OBS='/
      data ihutb/'TIME-BEG='/
      data ihutc/'TIME-END='/
c     data ihcra/'OBJCTRA ='/
c     data ihcde/'OBJCTDEC='/
      data ihcra/'RA      ='/
      data ihcar/'AR      ='/
      data ihcde/'DEC     ='/
      data ihfil /'FILTERS ='/
      data ihfilc/'FILTER  ='/
      data ihfilb/'COMMENT   FILTRO:'/
      data ihexp /'EXPTIME ='/
      data ihexpc/'EXPOSURE='/
      data ihexpb/'ITIME   ='/

c
c     IMCCE header keys
c

      data imccet/'TELESCOP='/
      data imcceb/'TM-START='/
      data imccee/'TM-END  ='/
      data imccer/'POSTN-RA='/
      data imcced/'POSTN-DE='/
      data imccef/'FLTRNM  ='/

c

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0
c


      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

      dj1=2400000.5D0

      debug = .false.
      nrec  = 1

c
c
c

c
c     Abre arquivo auxiliar 1
c
      sista=''
      sista='rm -f -r '

      imaux=''

      imaux(1:16)='fitsheader.1temp'


      do 1 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux(17:20)=jmaux(1:4)

      open(11,file=imaux,status='old',err=2)
      close (11)
 1    continue

 2    close (11)

      open(11,file=imaux)

      systa1=sista//imaux

c
c     Abre arquivo auxiliar 2
c
      sista=''
      sista='rm -f -r '

      kmaux=''

      kmaux(1:16)='fitsheader.2temp'


      do 3 i=1,9999

      write (jmaux,'(i4.4)') i

      kmaux(17:20)=jmaux(1:4)

      open(12,file=kmaux,status='old',err=4)
      close (12)
 3    continue

 4    close (12)

      open(12,file=kmaux)

      systa2=sista//kmaux

c

      rewind (11)
      rewind (12)






C
      open (1,file=infits,access='direct',form='unformatted',recl=2880)
C
C*************************************************************
C     Le e verifica o cabecalho da imagem                    *
C*************************************************************
C
 19   read (1,rec=nrec) header
C
C*************************************************************
C     Procura as dimensoes do arquivo                        *
C     As dimensoes estao no 4 e 5 registros NAXIS1 e NAXIS2  *
C*************************************************************

      if (khead.eq.6) go to 20
      if (khead.eq.7.or.khead.eq.8) go to 20

      k  = 0
      nx = 0
      ny = 0
      do 70 i = 250,320
         if (header(i).eq.' ') go to 70
         icomp = ichar(header(i))
         if (icomp.ge.48.and.icomp.le.57) then
         k = k + 1
         digit(k) = header(i)
         endif
 70   continue
C
      do 80 j = 1,k
         nx = nx + (ichar(digit(j)) - 48)*10**(k-j)
 80   continue
      k = 0
      do 90 i = 330,400
         if (header(i).eq.' ') go to 90
         icomp = ichar(header(i))
         if (icomp.ge.48.and.icomp.le.57) then
         k = k + 1
         digit(k) = header(i)
         endif
 90   continue
      do 100 j = 1,k
         ny = ny + (ichar(digit(j)) - 48)*10**(k-j)
 100  continue
      if (debug) then
         write (*,*) ' *** Matriz de ',nx,' X ',ny,' pixmat *** '
      endif


c     write (*,*) 'nx, ny ',nx,ny

c
c     Le cabecalhos (imagem LNA pode ter mas de um, fora do
c     padrao de 2880 bytes do FITS)
c

 
 
 20   read (1,rec=nrec) header
      write (11,30) header
 
 21   if (debug) then
      write (*,30) header
 30   format(80a1)
      endif
 
      do 40 k = 1,2880,80
         if (header(k).eq.'E') go to 35
         go to 40
 35      if (header(k+1).eq.'N') go to 36
         go to 40
 36      if (header(k+2).eq.'D') go to 50
 40   continue
 

c50   m = k - 1
c     do 60 i = 1,3
c        tend(i:i) = header(m+i)
c60   continue
c
c     if (tend.eq.'END') then
c        write (*,*) '  '
c        continue
c        else
c        write (*,*) '-> *** Erro: Cabecalho incompleto! *** <-'
c        write (*,*) '-> *** Leitura do Proximo Registro *** <-'

         nrec = nrec + 1
         go to 20
c     endif

 50   continue

      if (khead.eq.7.or.khead.eq.8) then

      nrec = nrec + 1
      read (1,rec=nrec) header
      write (11,30) header

      nrec = nrec + 1
      read (1,rec=nrec) header
      write (11,30) header

      endif


      close(1)

      rewind (11)
      rewind (12)




c
c     Extracao de dados do header do 1T1M do Pic du Midi, DEBUG
c

      if (khead.eq.14) then


c
c     (Ra,Dec), filtro e objeto inicializados
c

      iah=99
      iam=99
      sa=99.d0

      isig=mais
      idg=99
      idm=99
      ds=99.d0

      ichobj=''
      ichfil=''

c
c     Le palavras-chave do header
c


 6700 read (11,5000,err=6701,end=6701) itrima


c
c     Nome do objeto
c

      if (itrima(1:9).eq.'OBJECT  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      if (itrima(i:i).eq."/") itrima(i:i)=' '
      enddo
      ichobj=itrima(10:80)

      endif


c
c     Filtro
c

      if (itrima(1:9).eq.'FILTER  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      if (itrima(i:i).eq."/") itrima(i:i)=' '
      enddo
      ichfil=itrima(10:80)

      endif


c
c     Dimensao da matriz
c


      if (itrima(1:9).eq.'NAXIS1  =') then

      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) nx

      endif



      if (itrima(1:9).eq.'NAXIS2  =') then


      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) ny

      endif

 

c
c     Extrai data, hora e tempo de exposicao
c

      if (itrima(1:9).eq.'DATE-OBS=') then
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*) iutano,iutmes,iutdia,ia1,ia2,a3
      endif



c
c     Extrai tempo de exposicao
c

      if (itrima(1:9).eq.'EXPOSURE=') then
      itrima(1:9)='         '
      read (itrima,*) exps
      endif

c
c     (Ra,Dec)
c


      if (itrima(1:9).eq.'RA      =') then

      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) ra        

      ra=ra/15.d0
      iah=ra
      am=(ra-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0



      endif



      if (itrima(1:9).eq.'DEC     =') then

      isig='+'
      do i=11,80
      if (itrima(i:i).eq.'-') isig='-'
      enddo

      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) de

      if (de.lt.0) then
      de=-de
      isig='-'
      else
      isig='+'
      endif

      idg=de
      dm=(de-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0

      endif

c



      go to 6700

 6701 continue


      fd=hmsgms(ia1,ia2,a3+exps/2.d0)

      if (fd.gt.24.d0) then
      fd=fd-24.d0
      iutdia=iutdia+1
      endif
      
      fd=fd/24.d0

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      return
      endif





c
c     Extracao de dados do header da PRECAM/DES
c

      if (khead.eq.13) then


c
c     (Ra,Dec), filtro e objeto inicializados
c

      iah=99
      iam=99
      sa=99.d0

      isig=mais
      idg=99
      idm=99
      ds=99.d0

      ichobj=''
      ichfil=''

c
c     Le palavras-chave do header
c


 6690 read (11,5000,err=6691,end=6691) itrima


c
c     Nome do objeto
c

      if (itrima(1:9).eq.'OBJECT  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      enddo
      ichobj=itrima(10:80)

      endif


c
c     Filtro
c

      if (itrima(1:9).eq.'FILTER  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      enddo
      ichfil=itrima(10:80)

      endif


c
c     Dimensao da matriz
c

      if (itrima(1:9).eq.'NAXIS1  =') then

      do i=10,80
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) nx
      endif



      if (itrima(1:9).eq.'NAXIS2  =') then

      do i=10,80
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) ny
      endif

 

c
c     Extrai data, hora e tempo de exposicao
c

      if (itrima(1:9).eq.'UTSHUT  =') then
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*) iutano,iutmes,iutdia,ia1,ia2,a3
      endif



c
c     Extrai tempo de exposicao
c

      if (itrima(1:9).eq.'EXPREQ  =') then
      itrima(1:9)='         '
      read (itrima,*) exps
      endif

c
c     (Ra,Dec)
c


      if (itrima(1:9).eq.'RA      =') then

      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) ra
      ra=ra/15.d0
      iah=ra
      am=(ra-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0

      endif



      if (itrima(1:9).eq.'DEC     =') then

      do i=1,80
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'+'.and.itrima(i:i).ne.
     ?'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) de
      if (de.lt.0) then
      de=-de
      isig='-'
      else
      isig='+'
      endif

      idg=de
      dm=(de-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0

      endif

c



      go to 6690

 6691 continue


      fd=hmsgms(ia1,ia2,a3+exps/2.d0)

      if (fd.gt.24.d0) then
      fd=fd-24.d0
      iutdia=iutdia+1
      endif
      
      fd=fd/24.d0

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      return
      endif



c
c     Extracao de dados do header do 1T1M do Pic du Midi (ASTROMETRICA)
c

      if (khead.eq.12) then


c
c     (Ra,Dec), filtro e objeto inicializados
c

      iah=99
      iam=99
      sa=99.d0

      isig=mais
      idg=99
      idm=99
      ds=99.d0

      ichobj=''
      ichfil=''

c
c     Le palavras-chave do header
c


 6680 read (11,5000,err=6681,end=6681) itrima


c
c     Nome do objeto
c

      if (itrima(1:9).eq.'OBJECT  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      if (itrima(i:i).eq."/") itrima(i:i)=' '
      enddo
      ichobj=itrima(10:80)

      endif


c
c     Filtro
c

      if (itrima(1:9).eq.'FILTER  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      if (itrima(i:i).eq."/") itrima(i:i)=' '
      enddo
      ichfil=itrima(10:80)

      endif


c
c     Dimensao da matriz
c


      if (itrima(1:9).eq.'NAXIS1  =') then

      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) nx

      endif



      if (itrima(1:9).eq.'NAXIS2  =') then


      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) ny

      endif

 

c
c     Extrai data, hora e tempo de exposicao
c

      if (itrima(1:9).eq.'DATE-OBS=') then
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*) iutano,iutmes,iutdia,ia1,ia2,a3
      endif



c
c     Extrai tempo de exposicao
c

      if (itrima(1:9).eq.'EXPTIME =') then
      itrima(1:9)='         '
      read (itrima,*) exps
      endif

c
c     (Ra,Dec)
c


      if (itrima(1:9).eq.'OBJCTRA =') then

      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) iah,iam,sa

      endif



      if (itrima(1:9).eq.'OBJCTDEC=') then

      isig='+'
      do i=11,80
      if (itrima(i:i).eq.'-') isig='-'
      enddo

      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) idg,idm,ds

      endif

c



      go to 6680

 6681 continue


      fd=hmsgms(ia1,ia2,a3+exps/2.d0)

      if (fd.gt.24.d0) then
      fd=fd-24.d0
      iutdia=iutdia+1
      endif
      
      fd=fd/24.d0

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      return
      endif




c
c     Extracao de dados do header da Merlin/Raptor
c

      if (khead.eq.11) then

c
c     Nome do objeto e filtro ainda nao colocados no header
c

      ichobj=''
      ichfil=''

c
c     (Ra,Dec) inicializados
c

      iah=99
      iam=99
      sa=99.d0

      isig=mais
      idg=99
      idm=99
      ds=99.d0

c
c     Dimensao da matriz
c

 6670 read (11,5000,err=6671,end=6671) itrima



      do i=10,80
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'-'.and.itrima(i:i).
     ?ne.'+') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo


      if (itrima(1:9).eq.'NAXIS1  =') then
      itrima(1:9)='         '
      read (itrima,*) nx
      endif
   

      if (itrima(1:9).eq.'NAXIS2  =') then
      itrima(1:9)='         '
      read (itrima,*) ny
      endif
   

 

c
c     Extrai data, hora e tempo de exposicao
c


      if (itrima(1:9).eq.'TIMESTMP='.or.itrima(1:9).eq.'DATE-OBS=') then
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima,*) ib1,ib2,ib3,ia1,ia2,a3

      if (ib1.gt.1000) then
      iutano=ib1
      iutmes=ib2
      iutdia=ib3
      else
      iutano=ib3
      iutmes=ib2
      iutdia=ib1
      endif
      
      endif


c
c     Extrai tempo de exposicao
c

      if (itrima(1:9).eq.'EXP_TIME=') then
      itrima(1:9)='         '
      read (itrima,*) exps
      exps=exps/1000.d0
      endif

      if (itrima(1:9).eq.'EXPTIME =') then
      itrima(1:9)='         '
      read (itrima,*) exps
      endif


      go to 6670

c


 6671 continue

      rewind (11)

      fd=hmsgms(ia1,ia2,a3+exps/2.d0)

      if (fd.gt.24.d0) then
      fd=fd-24.d0
      iutdia=iutdia+1
      endif
      
      fd=fd/24.d0

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      return
      endif




c
c     Extracao de dados do header da IKON/LNA ou IXON/LNA
c

      if (khead.eq.10) then


c
c     (Ra,Dec), filtro, objeto inicializados
c

      iah=99
      iam=99
      sa=99.d0

      isig=mais
      idg=99
      idm=99
      ds=99.d0

      ichobj=''
      ichfil=''

c
c     Le palavras-chave do header
c


 6660 read (11,5000,err=6661,end=6661) itrima


c
c     Dimensao da matriz
c


      if (itrima(1:9).eq.'NAXIS1  =') then

      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) nx
      endif



      if (itrima(1:9).eq.'NAXIS2  =') then

      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:80),*) ny
      endif



c
c     Nome do objeto
c

      if (itrima(1:9).eq.'OBJECT  =') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      if (itrima(i:i).eq."/") itrima(i:i)=' '
      enddo
      ichobj=itrima(10:80)

      endif


c
c     Filtro
c

      if (itrima(1:6).eq.'FILTER') then

      do i=10,80
      if (itrima(i:i).eq."'") itrima(i:i)=' '
      if (itrima(i:i).eq."/") itrima(i:i)=' '
      enddo
      ichfil=itrima(10:80)

      endif

 

c
c     Extrai data, hora e tempo de exposicao
c

      if (itrima(1:9).eq.'FRAME   =') then
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*) iutano,iutmes,iutdia,ia1,ia2,a3
      endif



c
c     Extrai tempo de exposicao
c

      if (itrima(1:9).eq.'EXPOSURE=') then
      itrima(1:9)='         '
      do i=10,80
      if(itrima(i:i).eq.',') itrima(i:i)='.'
      enddo
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*) exps
      endif

c
c     (Ra,Dec)
c


      if (itrima(1:9).eq.'RA      =') then
      itrima(1:9)='         '
      do i=1,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*)iah,iam,sa
      endif

      if (itrima(1:9).eq.'DEC     =') then
      itrima(1:9)='         '
      isig='+'
      do i=11,80
      if (itrima(i:i).eq.'-') isig='-'
      enddo
      do i=10,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo
      read (itrima,*) idg,idm,ds
      endif


      go to 6660

 6661 continue


      fd=hmsgms(ia1,ia2,a3+exps/2.d0)

      if (fd.gt.24.d0) then
      fd=fd-24.d0
      iutdia=iutdia+1
      endif
      
      fd=fd/24.d0

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      return
      endif




c
c     Extracao de dados do header da S800/LNA
c

      if (khead.eq.9) then

c
c     Nome do objeto e filtro ainda nao colocados no header
c

      ichobj=''
      ichfil=''

c
c     (Ra,Dec) ainda nao colocados no header
c

      iah=99
      iam=99
      sa=99.d0

      isig=mais
      idg=99
      idm=99
      ds=99.d0

c
c     Dimensao da matriz
c

      iutano=0

 5560 read (11,5000,err=5561,end=5561) itrima




      if (itrima(1:9).eq.'NAXIS1  =') then
      itrima(1:9)='         '
      read (itrima,*) nx
      endif
   

      if (itrima(1:9).eq.'NAXIS2  =') then
      itrima(1:9)='         '
      read (itrima,*) ny
      endif


   

c
c     Extrai data, hora e tempo de exposicao
c



      if (itrima(1:9).eq.'DATE-OBS=') then
      do i=1,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo
      read (itrima,*) iutano,iutmes,iutdia
      endif



      if (itrima(1:9).eq.'DATE    =') then

      do i=10,80
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      itrima(1:9)='         '

      read (itrima,*) ib1,ib2,ib3            

      if (ib1.gt.1000) then
      kutano=ib1
      kutmes=ib2
      kutdia=ib3
      else
      kutano=ib3
      kutmes=ib2
      kutdia=ib1
      endif

      endif


      if (itrima(1:9).eq.'TIME    =') then

      do i=10,80
      if (itrima(i:i).eq.',') itrima(i:i)='.'
      enddo

      do i=10,80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      itrima(1:9)='         '

      read (itrima,*) ia1,ia2,a3,ia4,ia5,a6

      endif

c

      go to 5560

 5561 continue

      if (iutano.eq.0) then
      iutano=kutano
      iutmes=kutmes
      iutdia=kutdia
      endif


      fd1=hmsgms(ia1,ia2,a3)
      fd2=hmsgms(ia4,ia5,a6)

      if (fd2.lt.fd1) fd2=fd2+24.d0

      exps=dabs(fd2-fd1)*3600.d0

      fd=(fd1+fd2)/2.d0

      if (fd.gt.24.d0) then
      fd=fd-24.d0
      iutdia=iutdia+1
      endif
      
      fd=fd/24.d0

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      return
      endif






c
c     Extracao de todos os dados de header do CFHT/MegaCam
c

 5000 format(a80)

      if (khead.ne.6) go to 5020


 5005 continue

      itrima=''
      read (11,5000,err=5010,end=5010) itrima

      if (itrima(1:3).eq.'END') go to 5010

c
c     Dimensao da matriz
c

      if (itrima(1:9).eq.'NAXIS1  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) nx
      rewind (12)
      endif
   
      if (itrima(1:9).eq.'NAXIS2  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) ny
      rewind (12)
      endif
   


c
c     Nome do objeto
c

      if (itrima(1:9).eq.'OBJECT  =') then
      do i=12,80
      if (itrima(i:i).eq."'") go to 5007
      enddo
 5007 ichobj=''
      ichobj=itrima(12:i-1)
      endif


c
c     Filtro
c

      if (itrima(1:9).eq.'FILTER  =') then
      do i=12,80
      if (itrima(i:i).eq."'") go to 5008
      enddo
 5008 ichfil=''
      ichfil=itrima(12:i-1)
      rewind (12)
      endif


c
c     Exposicao, data juliana modificada
c


      if (itrima(1:9).eq.'EXPTIME =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) exps
      rewind (12)
      endif

      if (itrima(1:9).eq.'MJD-OBS =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read (12,*) djm
      rewind (12)
      endif

c
c     (RA,DEC)
c

      if (itrima(1:9).eq.'CRVAL1  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) crval1
      rewind (12)
      endif

      if (itrima(1:9).eq.'CRVAL2  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) crval2
      rewind (12)
      endif

      if (itrima(1:9).eq.'CRPIX1  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) crpix1
      rewind (12)
      endif

      if (itrima(1:9).eq.'CRPIX2  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) crpix2
      rewind (12)
      endif

      if (itrima(1:9).eq.'CD1_1   =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) cd11
      rewind (12)
      endif

      if (itrima(1:9).eq.'CD1_2   =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) cd12
      rewind (12)
      endif

      if (itrima(1:9).eq.'CD2_1   =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) cd21
      rewind (12)
      endif

      if (itrima(1:9).eq.'CD2_2   =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) cd22
      rewind (12)
      endif

c

      go to 5005

 5010 continue


c
c     Calculo do instante central, data juliana e data gregoriana
c

      djm=djm+(exps/2.d0)/86400.d0

      dj=djm+dj1

      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,j)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps


c
c     Calculo do (RA,DEC) do centro do CCD
c




      x=nx/2.d0
      y=ny/2.d0

      x=x-crpix1
      y=y-crpix2

      polx=x*cd11+y*cd12
      poly=x*cd21+y*cd22


      ra=polx+crval1
      de=poly+crval2


      ra=ra/15.d0
      iah=ra
      iam=(ra-iah)*60.d0
      sa=((ra-iah)*60.d0-iam)*60.d0

      isig=mais
      if (de.lt.0.d0) isig=menos
      de=dabs(de)

      idg=de
      idm=(de-idg)*60.d0
      ds=((de-idg)*60.d0-idm)*60.d0

c

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



      return

c

 5020 continue



c
c     Extracao de todos os dados de header do 1.5m Tubitak Russia-Turquia
c     e do Las Campanas-40cm.


      if (khead.lt.7.or.khead.gt.8) go to 5040


      exps=0.d0

 5025 continue

      itrima=''
      read (11,5000,err=5037,end=5037) itrima


c
c     Dimensao da matriz
c

      if (itrima(1:9).eq.'NAXIS1  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) nx
      rewind (12)
      endif
   
      if (itrima(1:9).eq.'NAXIS2  =') then
      write (12,*) itrima(10:31)
      rewind (12)
      read(12,*) ny
      rewind (12)
      endif
   


c
c     Nome do objeto
c

      if (itrima(1:9).eq.'OBJECT  =') then
      do i=29,10,-1
      if (itrima(i:i).eq.iplic) go to 5027
      enddo
 5027 ichobj=''
      if (i.le.11) i=29
      ichobj=itrima(10:i)
      endif


c
c     Filtro
c

      if (itrima(1:9).eq.'FILTER  =') then
      do i=29,10,-1
      if (itrima(i:i).eq.iplic) go to 5028
      enddo
 5028 ichfil=''
      if (i.le.11) i=29
      ichfil=itrima(10:i)
      endif


c
c     Exposicao, data juliana modificada
c


      if (itrima(1:9).eq.'EXPTIME =') then
      write (12,*) itrima(10:30)
      rewind (12)
      read(12,*) exps
      rewind (12)
      endif

      if (itrima(1:9).eq.'DATE-OBS=') then
      write (12,*) itrima(12:21)
      rewind (12)
      read (12,5030) iutano,iutmes,iutdia
 5030 format(1x,i4.4,1x,i2.2,1x,i2.2)
      rewind (12)
      endif

      if (itrima(1:9).eq.'TIME-OBS=') then
      write (12,*) itrima(12:23)
      rewind (12)
      read (12,5031) ih,im,ss
 5031 format(1x,i2.2,1x,i2.2,1x,f6.3)
      rewind (12)
      endif

      if (itrima(1:9).eq.'UTSTART =') then
      write (12,*) itrima(12:19)
      rewind (12)
      read (12,*) ih,im,ss
      rewind (12)
      endif


c
c     (RA,DEC)
c

      if (itrima(1:9).eq.'RA      =') then
      do i=10,80
      if (itrima(i:i).eq.iplic) go to 5032
      enddo
 5032 i1=i+1
      do i=80,10,-1
      if (itrima(i:i).eq.iplic) go to 5033
      enddo
 5033 i2=i-1
      do i=i1,i2
      if (itrima(i:i).eq.idois) itrima(i:i)=ibrac
      enddo
      write (12,*) itrima(i1:i2)
      rewind (12)
      read(12,*) iah,iam,sa
      rewind (12)
      endif


      if (itrima(1:9).eq.'DEC     =') then
      do i=10,80
      if (itrima(i:i).eq.iplic) go to 5034
      enddo
 5034 i1=i+1
      do i=80,10,-1
      if (itrima(i:i).eq.iplic) go to 5035
      enddo
 5035 i2=i-1
      isig=mais
      do i=i1,i2
      if (itrima(i:i).eq.menos) isig=menos
      enddo
      do i=i1,i2
      if (itrima(i:i).eq.menos) itrima(i:i)=ibrac
      enddo
      do i=i1,i2
      if (itrima(i:i).eq.idois) itrima(i:i)=ibrac
      enddo
      write (12,*) itrima(i1:i2)
      rewind (12)
      read(12,*) idg,idm,ds
      rewind (12)
      endif

c

      go to 5025

 5037 continue


c
c     Calculo das datas julianas
c


      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      dj=djm+djm0

      sut=ss+exps/2.d0

      fd=hmsgms(ih,im,sut)/24.d0

      if (fd.gt.1.d0) then
      fd=fd-1.d0
      dj=dj+1.d0
      endif

      dj=dj+fd


c
c     Calculo do instante central e data gregoriana
c

      djm=dj-dj1

      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,j)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      iexps=exps

c


      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



      return

c

 5040 continue


c
c     Pegando os dados do header da imagem
c

      rewind (11)

199   format(a9,1x,a20)
200   format(a9,1x,70a1)


c
c     Pega nome do objeto
c

      do i=1,10000000
      read (11,199,err=201,end=201) iquem,ichobj
      if (iquem.eq.ihobj) go to 201
      enddo
201   rewind (11)


c
c     Pega ascensao reta do header
c

c
c     Haute Provence (RA,Dec), Instant of observation (UT), Exptime, filters
c

      do i=1,10000000
      read (11,200,end=1999) iquem,(ler(j),j=1,70)
      if (iquem.eq.imccet) go to 1500
      enddo

 1500 continue

      do i=1,64
      if (ler(i)  .eq.'O') then
      if (ler(i+1).eq.'H') then
      if (ler(i+2).eq.'P') then
      if (ler(i+3).eq.'-') then
      if (ler(i+4).eq.'1') then
      if (ler(i+5).eq.'2') then
      if (ler(i+6).eq.'0') then
      go to 1501
      endif
      endif
      endif
      endif
      endif
      endif
      endif

      enddo

      go to 1999

c     

 1501 continue

      rewind (11)

c
c     Geting (RA,Dec) Haute Provence
c
      
      do i=1,10000000
      read (11,1502) iquem
 1502 format(a9)
      if (iquem.eq.imccer) then
      backspace 11
      go to 1503
      endif
      enddo

 1503 read (11,1504) rarara
 1504 format(14x,f9.5)

      rarara=rarara/15.d0
      iah=rarara
      iam=(rarara-iah)*60.d0
      sa=((rarara-iah)*60.d0-iam)*60.d0

      rewind (11)

      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imcced) then
      backspace 11
      go to 1505
      endif
      enddo

 1505 read (11,1504) dedede

      isig=mais
      if (dedede.lt.0.d0) isig=menos
      dedede=dabs(dedede)

      idg=dedede
      idm=(dedede-idg)*60.d0
      ds=((dedede-idg)*60.d0-idm)*60.d0

c
c     Geting Instant of observation (UT) and Exptime
c

      rewind (11)


      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imcceb) then
      backspace 11
      go to 1513
      endif
      enddo

 1513 read (11,1515) iuth,iutm,sut
 1515 format(43x,i2,1x,i2,1x,f5.2)

      rewind (11)


      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imccee) then
      backspace 11
      go to 1516
      endif
      enddo

 1516 read (11,1515) juth,jutm,sutj


      tfi=hmsgms(juth,jutm,sutj)
      tin=hmsgms(iuth,iutm,sut)

c
c     here, if tfi after midnight, tin before midnight
c

      if (tfi.lt.tin) tfi=tfi+24.d0

c

      exps=(tfi-tin)*3600.d0
      iexps=exps+0.1

c
c     Computes UT mean instant (UT may be greater than 24hs)
c

      sut=sut+iexps/2.d0

      if (sut.ge.60.d0) then
      idiv=sut/60.d0
      sut=sut-idiv*60.d0
      iutm=iutm+idiv
      endif

      if (iutm.ge.60) then
      idiv=iutm/60.d0
      iutm=iutm-idiv*60.d0
      iuth=iuth+idiv
      endif


c
c     Geting Gregorian date of observation
c

      rewind (11)

      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.ihdat) then
      backspace 11
      go to 1517
      endif
      enddo

 1517 read (11,1518) iutdia,iutmes,iutano
 1518 format(11x,i2,1x,i2,1x,i4)

c
c     Geting the filter used
c


      rewind (11)

      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imccef) then
      backspace 11
      go to 1519
      endif
      enddo

 1519 read (11,1520) ichfil
 1520 format(20x,a20)

c
c     Go to computing JD
c

      go to 250


c
c     Other observatories, ESO, Itajuba, Bulgaria
c 

c

 1999 continue

      rewind (11)

      do i=1,10000000
      read (11,200,end=207) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihcra) go to 202
      if (iquem.eq.ihcar) go to 202
      enddo
 202  rewind (11)


c
c     SBIG Pic du Midi RA
c

      jjjjjj=0
      do i=1,70
      if (ler(i).eq.'E') jjjjjj=1
      enddo

      if (jjjjjj.eq.1.and.khead.eq.5) then
      write (12,*) ler
      rewind (12)
      read (12,*) rarara
      rewind (12)
      rarara=rarara/15.d0
      iah=rarara
      iam=(rarara-iah)*60.d0
      sa=((rarara-iah)*60.d0-iam)*60.d0
      go to 206
      endif



c
c     RA imagens ESO
c


      do i=1,11
      if (ler(i).ne.ibrac) go to 2999
      enddo

      rewind (12)

      write (12,*) ler
      rewind (12)
      read (12,*) rarara
      rewind (12)

      rarara=rarara/15.d0
      iah=rarara
      iam=(rarara-iah)*60.d0
      sa=((rarara-iah)*60.d0-iam)*60.d0

      go to 206

c
c     Demais imagens, LNA, Roenos, etc
c

 2999 do i=2,19
      if (ler(i).ne.ibrac) go to 1001
      enddo

      go to 207

c

 1001 do i=1,70
      if (ler(i).eq.iplic) go to 203
      enddo

 203  rewind (12)

      i1=i

      ipu=0
      if (ler(i1+1).eq.ibrac) ipu=1

      do i=i1+1+ipu,70
      if (ler(i).eq.idois) go to 400
      if (ler(i).eq.ibrac) go to 400
      enddo

 400  i2=i

      do i=i2+2,70
      if (ler(i).eq.idois) go to 401
      if (ler(i).eq.ibrac) go to 401
      enddo

 401  i3=i

      do i=i3+2,70
      if (ler(i).eq.iplic) go to 402
      enddo

 402  i4=i


      write (12,*) (ler(j),j=i1+1,i2-1),' ',(ler(k),k=i2+1,i3-1),
     ?' ',(ler(l),l=i3+1,i4-1)
      rewind (12)

      read (12,*,err=207) iah,iam,sa


c     write (12,204) (ler(j),j=i+1,i+2),(ler(k),k=i+4,i+5),
c    ?(ler(l),l=i+7,i+11)
c204  format(2a1,2a1,5a1)
c     rewind (12)
c     read (12,205,err=315) iah,iam,sa
c205  format(2i2,f5.2)

c     rewind (12)

c     go to 206




      rewind (12)

      go to 206

c

c315  rewind (12)

c     read (12,320) iah,iam,isa
c320  format(3i2)
      sa=isa

c     rewind (12)

c     go to 206

c

 207  rewind (11)
      rewind (12)

      iah=99
      iam=99
      sa=99.999d0

c


 206  continue



c
c     Pega declinacao do header
c


      isig=mais

      do i=1,10000000
      read (11,200,end=219) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihcde) go to 212
      enddo
 212  rewind (11)


c
c     SBIG Pic du Midi DEC
c


      if (jjjjjj.eq.1.and.khead.eq.5) then
      rewind (12)
      write (12,*) ler
      rewind (12)
      read (12,*) dedede
      rewind (12)
      isig=mais
      if (dedede.lt.0.d0) isig=menos
      dedede=dabs(dedede)
      idg=dedede
      idm=(dedede-idg)*60.d0
      ds=((dedede-idg)*60.d0-idm)*60.d0
      go to 220
      endif

c
c     Imagens ESO
c


      do i=1,11
      if (ler(i).ne.ibrac) go to 2998
      enddo

      rewind (12)

      write (12,*) ler
      rewind (12)
      read (12,*) dedede
      rewind (12)

      isig=mais
      if (dedede.lt.0.d0) isig=menos
      dedede=dabs(dedede)

      idg=dedede
      idm=(dedede-idg)*60.d0
      ds=((dedede-idg)*60.d0-idm)*60.d0

      go to 220



c
c      Demais imagens LNA, Romenos, etc ...
c



 2998 do i=2,19
      if (ler(i).ne.ibrac) go to 1002
      enddo

      go to 219

c


 1002 do i=1,10

      if (ler(i).eq.menos) then
      isig=menos
      ler(i)=ibrac
      go to 411
      endif

      if (ler(i).eq.mais) then
      isig=mais
      ler(i)=ibrac
      go to 411
      endif

      enddo

c

 411  do i=1,70
      if (ler(i).eq.iplic) go to 412
      enddo

 412  i1=i

      ipu=0
      if (ler(i1+1).eq.ibrac) ipu=1

      do i=i1+1+ipu,70
      if (ler(i).eq.idois) go to 415
      if (ler(i).eq.ibrac) go to 415
      enddo

 415  i2=i

      do i=i2+2,70
      if (ler(i).eq.idois) go to 417
      if (ler(i).eq.ibrac) go to 417
      enddo

 417  i3=i

      do i=i3+2,70
      if (ler(i).eq.iplic) go to 418
      enddo

 418  i4=i


      write (12,*) (ler(j),j=i1+1,i2-1),' ',(ler(k),k=i2+1,i3-1),
     ?' ',(ler(l),l=i3+1,i4-1)
      rewind (12)

      read (12,*,err=219) idg,idm,ds

      rewind (12)

      go to 220


c     do i=1,70
c     if (ler(i).eq.iplic) go to 213
c     enddo

c213  i1=i


c     do i=i1+1,70
c     if ((ler(i).eq.idois).or.(ler(i).eq.ibrac)) go to 214
c     enddo

c214  i2=i


c     do i=i2+1,70
c     if ((ler(i).eq.idois).or.(ler(i).eq.ibrac)) go to 215
c     enddo

c215  i3=i

c     do i=i1+1,i2-1
c     if (ler(i).eq.menos) ler(i)='0'
c     if (ler(i).eq.mais ) ler(i)='0'
c     enddo

c     icasas=i2-i1-1

c     do j=i1+1,i3+5
c     if (ler(j).eq.ibrac) ler(j)='0'
c     enddo



c     if (icasas.eq.1) write (12,216) ler(i1+1),
c    ?(ler(j),j=i2+1,i3-1),(ler(k),k=i3+1,i3+5)
c     if (icasas.eq.2) write (12,217) (ler(l),l=i1+1,i1+2),
c    ?(ler(j),j=i2+1,i3-1),(ler(k),k=i3+1,i3+5)
c     if (icasas.eq.3) write (12,217) (ler(l),l=i1+2,i1+3),
c    ?(ler(j),j=i2+1,i3-1),(ler(k),k=i3+1,i3+5)

c216  format('0',a1,2a1,5a1)
c217  format(2a1,2a1,5a1)

c     rewind (12)

c     read (12,218,err=325) idg,idm,ds
c218  format(2i2,f5.2)

c     rewind (12)


c     go to 220

c

c325  rewind (12)

c     read (12,330) idg,idm,ids
c330  format(3i2)
c     ds=ids

c     rewind (12)

c     go to 220

c

 219  rewind (11)
      rewind (12)

      isig='+'
      idg=99
      idm=99
      ds=99.999d0

c

 220  continue


c
c     Pega tempo universal UT do header
c


      do i=1,10000000
      read (11,200,end=221) iquem,(ler(j),j=1,70)

c
c     Le SBIG STL
c

      if (iquem.eq.ihdat .and. khead.eq.5) then

      rewind (12)

      write (12,*) (ler(j),j=13,14)
      write (12,*) (ler(j),j=16,17)
      write (12,*) (ler(j),j=19,24)

      rewind (12)

      read (12,*) iuth
      read (12,*) iutm
      read (12,*) sut
 
      rewind (11)
      rewind (12)

      go to 510
      endif


c
c     Demais instrumentos
c


      if (iquem.eq.ihut) go to 222
      if (iquem.eq.ihutt) go to 500
      if (iquem.eq.ihuts) go to 3010
      enddo

c
c     Imagens ESO
c

 3010 rewind (11)
      rewind (12)
      write (12,*) ler
      rewind (12)
      read (12,*) sut
      rewind (12)

      hora=sut/3600.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      go to 510

c
c     Aqui,nem UT nem EXPTIME sao dados; LNA dah Tini e Tfinal
c     em TL. Calculamos aqui UT inicial e EXPTIME
c

 221  rewind (11)



      do i=1,10000000
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihutb) go to 335
      enddo

 335  do i=1,70
      if (ler(i).eq.idois) go to 336
      enddo

 336  i1=i-3

      write (12,337) (ler(j),j=i1+1,i1+2),(ler(i),i=i1+4,i1+5),
     ?(ler(k),k=i1+7,i1+8),(ler(l),l=i1+10,i1+12)
 337  format(2a1,2a1,2a1,'.',3a1)

      read (11,200) iquem,(ler(j),j=1,70)
      
      write (12,337) (ler(j),j=i1+1,i1+2),(ler(i),i=i1+4,i1+5),
     ?(ler(k),k=i1+7,i1+8),(ler(l),l=i1+10,i1+12)

      rewind (11)
      rewind (12)

      read (12,340) iuth,iutm,sut
 340  format(i2,i2,f6.3)

      iuth=iuth+3

      read (12,340) juth,jutm,sutj

      rewind (12)

      juth=juth+3

      tfi=hmsgms(juth,jutm,sutj)
      tin=hmsgms(iuth,iutm,sut)

c
c     tfi depois da meia noite, tin antes da meia noite
c


      if (tfi.lt.tin) tfi=tfi+24.d0

c

      exps=(tfi-tin)*3600.d0
      iexps=exps+0.1

      go to 230      
c



 222  rewind (11)

      do i=1,70
      if (ler(i).eq.iplic) then
      ii=i
      go to 223
      endif

      if (ler(i).eq.idois) then
      ii=i-3
      go to 223
      endif

      enddo

 223  rewind (12)

      write (12,224) (ler(j),j=ii+1,ii+2),(ler(k),k=ii+4,ii+5),
     ?(ler(l),l=ii+7,ii+11)
 224  format(2a1,2a1,5a1)
      rewind (12)
      read (12,225) iuth,iutm,sut
 225  format(2i2,f5.2)

      rewind (12)

c

      go to 510

c

 500  rewind (11)

      do i=1,70
      if (ler(i).eq.iplic) then
      ii=i
      go to 503
      endif

      if (ler(i).eq.idois) then
      ii=i-3
      go to 503
      endif

      enddo

 503  rewind (12)

      write (12,504) (ler(j),j=ii+1,ii+2),(ler(k),k=ii+4,ii+5),
     ?(ler(l),l=ii+7,ii+8)
 504  format(2a1,2a1,2a1)
      rewind (12)
      read (12,505) iuth,iutm,isut
 505  format(3i2)

      rewind (12)

      sut=isut

 510  continue


c
c     Pega tempo de exposicao para calculo do instante medio UT
c



      do i=1,10000000
 227  read (11,199,err=227,end=231) iquem,ichexp
      if (iquem.eq.ihexp) go to 228
      if (iquem.eq.ihexpc) go to 228
      if (iquem.eq.ihexpb) go to 228
      enddo

 231  iexps=0
      rewind (11)
      go to 230
c

 228  rewind (11)

      write (12,229) ichexp
 229  format(a20)
      rewind (12)
      read (12,*) exps
      iexps=exps

      rewind (12)


 230  continue

c
c     Calcula instante UT medio (UT pode ultrapassar 24hs)
c

c     sut=sut+iexps/2.d0

      sut=sut+exps/2.d0

      if (sut.ge.60.d0) then
      idiv=sut/60.d0
      sut=sut-idiv*60.d0
      iutm=iutm+idiv
      endif

      if (iutm.ge.60) then
      idiv=iutm/60.d0
      iutm=iutm-idiv*60.d0
      iuth=iuth+idiv
      endif



c
c     Pega data gregoriana de Greenwhich do header
c


      do i=1,10000000
      read (11,200,err=232,end=232) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihdat) go to 232
      enddo
 232  rewind (11)

      do i=1,70
      if (ler(i).eq.iplic) then
      ii=i
      go to 233
      endif

      if (ler(i).eq.menos) then
      ii=i-5
      go to 233
      endif
      enddo

 233  rewind (12)

      write (12,234) (ler(j),j=ii+1,ii+10)
 234  format(10a1)
      rewind (12)
      read (12,235,err=300) iutano,iutmes,iutdia
 235  format(i4,1x,i2,1x,i2)

      rewind (12)


      go to 236

c
c     Outro tipo de formatacao de data
c

 300  rewind (12)


      read (12,310,err=311) iutdia,iutmes,iutano
 310  format(i2,1x,i2,1x,i4)

      rewind (12)
      go to 236

c

 311  rewind (12)

      read (12,312) iutdia,iutmes,iutano
 312  format(i2,1x,i2,1x,i2)

      rewind (12)

      if (iutano.lt.10) iutano=iutano+100

      iutano=iutano+1900


c

 236  continue




c
c     Pega filtro usado
c

      rewind (11)

      do i=1,10000000
      read (11,200,end=237) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihfil)  go to 240
      if (iquem.eq.ihfilc) go to 240
      enddo

c
c     Header LNA incompleto, filtro nos COMMENTS
c 

 237  rewind (11)


      do i=1,10000000
      read (11,238,end=249) iqual,ifa
 238  format(a17,a5)
      if (iqual.eq.ihfilb) go to 239
      enddo

 239  rewind (11)

      ichfil=ifa

      go to 250

c

 249  ichfil='filtro desconhecido'

      rewind (11)

      go to 250

c

 240  rewind (11)


      do i=1,70
      if (ler(i).eq.iplic) go to 241
      enddo

 241  i1=i+1

      do i=70,1,-1
      if (ler(i).eq.iplic) go to 242
      enddo

 242  i2=i-1

      jj=i2-i1+1

      do j=1,jj
      ler(j)=ler(i1+j-1)
      enddo

      do i=jj+1,70
      ler(i)=ibrac
      enddo

      write (12,243) (ler(j),j=1,20)
 243  format(20a1)
      rewind (12) 

      read (12,244) ichfil
 244  format(a20)

      rewind (12)
c

 250  continue

c


      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      fd=hmsgms(iuth,iutm,sut)/24.d0

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0

c


      return
      end







c
c     Subrotina cal2JD
c
c
      SUBROUTINE iau_CAL2JD ( IY, IM, ID, DJM0, DJM, J )
*+
*  - - - - - - - - - - -
*   i a u _ C A L 2 J D
*  - - - - - - - - - - -
*
*  Gregorian Calendar to Julian Date.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     IY,IM,ID    i     year, month, day in Gregorian calendar (Note 1)
*
*  Returned:
*     DJM0        d     MJD zero-point: always 2400000.5
*     DJM         d     Modified Julian Date for 0 hrs
*     J           i     status:
*                           0 = OK
*                          -1 = bad year   (Note 3: JD not computed)
*                          -2 = bad month  (JD not computed)
*                          -3 = bad day    (JD computed)
*
*  Notes:
*
*  1) The algorithm used is valid from -4800 March 1, but this
*     implementation rejects dates before -4799 January 1.
*
*  2) The Julian Date is returned in two pieces, in the usual SOFA
*     manner, which is designed to preserve time resolution.  The
*     Julian Date is available as a single number by adding DJM0 and
*     DJM.
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
*  This revision:  2000 December 15
*
*  Copyright (C) 2001 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION DJM0, DJM
      INTEGER J, MY, IYPMY

*  Earliest year allowed (4800BC)
      INTEGER IYMIN
      PARAMETER ( IYMIN = -4799 )

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Preset status.
      J = 0

*  Validate year.
      IF ( IY.LT.IYMIN ) THEN
         J = -1
      ELSE

*     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN

*        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 ) MTAB(2) = 28

*        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J = -3

*        Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( ( 1461 * ( IYPMY + 4800 ) ) / 4
     :                + (  367 * ( IM-2 - 12*MY ) ) / 12
     :                - (    3 * ( ( IYPMY + 4900 ) / 100 ) ) / 4
     :                + ID - 2432076)

*        Bad month
         ELSE
            J = -2
         END IF
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






c
c     Subrotina OBHAAD
c
C     Extrai dados do header de imagens fits
C

      subroutine obhaad (infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,iexps,nx,ny,ofra,
     ?ofde,khead)

      IMPLICIT REAL *8 (A-H,O-Z)

      dimension header(2880),digit(4)
      character*3  tend
      character*1  digit
      character*1  header
      character*150 infits
      character*9  ihdat,ihut,ihuts,ihutt,ihcra,ihcar,ihcde,ihobj,ihfil,
     ?ihexp,iquem,ihutm,ihutb,ihutc,ihfilc,ihexpb,imagek
      character*20 ichobj,ichfil
      character*70 ichexp
      character*1  ler(70),iplic,ibrac,idois,isig,mais,menos
      character*17 ihfilb,iqual
      character*5  ifa

      character*9  imccet,imcceb,imccee,imccer,imcced,imccef,imagei,
     ?itrimi
      character*80 itrima,iran,iden


      character*20 imaux,kmaux
      character*4 jmaux
      character*9 sista
      character*29 systa1,systa2


      logical      debug

      data iplic/"'"/
      data idois/':'/
      data ibrac/' '/
      data mais /'+'/
      data menos/'-'/

      data ihobj/'OBJECT  ='/
      data ihdat/'DATE-OBS='/
      data ihut /'UT      ='/
      data ihuts/'UTC     ='/
      data ihutm/'TM-START='/
      data ihutt/'TIME-OBS='/
      data ihutb/'TIME-BEG='/
      data ihutc/'TIME-END='/
c     data ihcra/'OBJCTRA ='/
c     data ihcde/'OBJCTDEC='/
      data ihcra/'RA      ='/
      data ihcar/'AR      ='/
      data ihcde/'DEC     ='/
      data ihfil /'FILTERS ='/
      data ihfilc/'FILTER  ='/
      data ihfilb/'COMMENT   FILTRO:'/
      data ihexp /'EXPTIME ='/
      data ihexpb/'ITIME   ='/
      data imagei/'IMAGEID ='/
      data itrimi/'TRIM    ='/


c
c     IMCCE header keys
c

      data imccet/'TELESCOP='/
      data imcceb/'TM-START='/
      data imccee/'TM-END  ='/
      data imccer/'POSTN-RA='/
      data imcced/'POSTN-DE='/
      data imccef/'FLTRNM  ='/

c

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0
c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

      dj1=2400000.5D0

c

      debug = .false.
      nrec  = 1



c
c     Abre arquivo auxiliar 1
c

      sista='rm -f -r '

      imaux=''

      imaux(1:16)='fitsheader.1temp'


      do 1 i=1,9999

      write (jmaux,'(i4.4)') i

      imaux(17:20)=jmaux(1:4)

      open(11,file=imaux,status='old',err=2)
      close (11)
 1    continue

 2    close (11)

      open(11,file=imaux)

      systa1=sista//imaux

c
c     Abre arquivo auxiliar 2
c

      sista='rm -f -r '

      kmaux=''

      kmaux(1:16)='fitsheader.2temp'


      do 3 i=1,9999

      write (jmaux,'(i4.4)') i

      kmaux(17:20)=jmaux(1:4)

      open(12,file=kmaux,status='old',err=4)
      close (12)
 3    continue

 4    close (12)

      open(12,file=kmaux)

      systa2=sista//kmaux

c

      rewind (11)
      rewind (12)






C
      open (1,file=infits,access='direct',form='unformatted',recl=2880)
C
C*************************************************************
C     Le e verifica o cabecalho da imagem                    *
C*************************************************************
C
 19   read (1,rec=nrec) header
C
C*************************************************************
C     Procura as dimensoes do arquivo                        *
C     As dimensoes estao no 4 e 5 registros NAXIS1 e NAXIS2  *
C*************************************************************
      k  = 0
      nx = 0
      ny = 0
      do 70 i = 250,320
         if (header(i).eq.' ') go to 70
         icomp = ichar(header(i))
         if (icomp.ge.48.and.icomp.le.57) then
         k = k + 1
         digit(k) = header(i)
         endif
 70   continue
C
      do 80 j = 1,k
         nx = nx + (ichar(digit(j)) - 48)*10**(k-j)
 80   continue
      k = 0
      do 90 i = 330,400
         if (header(i).eq.' ') go to 90
         icomp = ichar(header(i))
         if (icomp.ge.48.and.icomp.le.57) then
         k = k + 1
         digit(k) = header(i)
         endif
 90   continue
      do 100 j = 1,k
         ny = ny + (ichar(digit(j)) - 48)*10**(k-j)
 100  continue
      if (debug) then
         write (*,*) ' *** Matriz de ',nx,' X ',ny,' pixmat *** '
      endif


c     write (*,*) 'nx, ny ',nx,ny

c
c     Le cabecalhos (imagem LNA pode ter mas de um, fora do
c     padrao de 2880 bytes do FITS)
c

 
 
 20   read (1,rec=nrec) header
      write (11,30) header
 
 21   if (debug) then
      write (*,30) header
 30   format(80a1)
      endif
 
      do 40 k = 1,2880,80
         if (header(k).eq.'E') go to 35
         go to 40
 35      if (header(k+1).eq.'N') go to 36
         go to 40
 36      if (header(k+2).eq.'D') go to 50
 40   continue
 

c50   m = k - 1
c     do 60 i = 1,3
c        tend(i:i) = header(m+i)
c60   continue
c
c     if (tend.eq.'END') then
c        write (*,*) '  '
c        continue
c        else
c        write (*,*) '-> *** Erro: Cabecalho incompleto! *** <-'
c        write (*,*) '-> *** Leitura do Proximo Registro *** <-'

         nrec = nrec + 1
         go to 20
c     endif

 50   continue


      close(1)


c
c     Pegando os dados do header da imagem
c

      rewind (11)

199   format(a9,1x,a20)
200   format(a9,1x,70a1)

c
c     Pega nome do objeto
c

      do i=1,10000000
      read (11,199,err=201,end=201) iquem,ichobj
      if (iquem.eq.ihobj) go to 201
      enddo
201   rewind (11)


c
c     Pega ascensao reta do header
c

c
c     Haute Provence (RA,Dec), Instant of observation (UT), Exptime, filters
c

      do i=1,10000000
      read (11,200,end=1999) iquem,(ler(j),j=1,70)
      if (iquem.eq.imccet) go to 1500
      enddo

 1500 continue

      do i=1,64
      if (ler(i)  .eq.'O') then
      if (ler(i+1).eq.'H') then
      if (ler(i+2).eq.'P') then
      if (ler(i+3).eq.'-') then
      if (ler(i+4).eq.'1') then
      if (ler(i+5).eq.'2') then
      if (ler(i+6).eq.'0') then
      go to 1501
      endif
      endif
      endif
      endif
      endif
      endif
      endif

      enddo

      go to 1999

c     

 1501 continue

      rewind (11)

c
c     Geting (RA,Dec) Haute Provence
c
      
      do i=1,10000000
      read (11,1502) iquem
 1502 format(a9)
      if (iquem.eq.imccer) then
      backspace 11
      go to 1503
      endif
      enddo

 1503 read (11,1504) rarara
 1504 format(14x,f9.5)

      rarara=rarara/15.d0
      iah=rarara
      iam=(rarara-iah)*60.d0
      sa=((rarara-iah)*60.d0-iam)*60.d0

      rewind (11)

      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imcced) then
      backspace 11
      go to 1505
      endif
      enddo

 1505 read (11,1504) dedede

      isig=mais
      if (dedede.lt.0.d0) isig=menos
      dedede=dabs(dedede)

      idg=dedede
      idm=(dedede-idg)*60.d0
      ds=((dedede-idg)*60.d0-idm)*60.d0

c
c     Geting Instant of observation (UT) and Exptime
c

      rewind (11)


      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imcceb) then
      backspace 11
      go to 1513
      endif
      enddo

 1513 read (11,1515) iuth,iutm,sut
 1515 format(43x,i2,1x,i2,1x,f5.2)

      rewind (11)


      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imccee) then
      backspace 11
      go to 1516
      endif
      enddo

 1516 read (11,1515) juth,jutm,sutj


      tfi=hmsgms(juth,jutm,sutj)
      tin=hmsgms(iuth,iutm,sut)

c
c     here, if tfi after midnight, tin before midnight
c

      if (tfi.lt.tin) tfi=tfi+24.d0

c

      exps=(tfi-tin)*3600.d0
      iexps=exps+0.1

c
c     Computes UT mean instant (UT may be greater than 24hs)
c

      sut=sut+iexps/2.d0

      if (sut.ge.60.d0) then
      idiv=sut/60.d0
      sut=sut-idiv*60.d0
      iutm=iutm+idiv
      endif

      if (iutm.ge.60) then
      idiv=iutm/60.d0
      iutm=iutm-idiv*60.d0
      iuth=iuth+idiv
      endif


c
c     Geting Gregorian date of observation
c

      rewind (11)

      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.ihdat) then
      backspace 11
      go to 1517
      endif
      enddo

 1517 read (11,1518) iutdia,iutmes,iutano
 1518 format(11x,i2,1x,i2,1x,i4)

c
c     Geting the filter used
c


      rewind (11)

      do i=1,10000000
      read (11,1502) iquem
      if (iquem.eq.imccef) then
      backspace 11
      go to 1519
      endif
      enddo

 1519 read (11,1520) ichfil
 1520 format(20x,a20)

c
c     Go to computing JD
c

      go to 250


c
c     Other observatories, ESO, Itajuba, Bulgaria
c 

c

 1999 continue

      rewind (11)

      do i=1,10000000
      read (11,200,end=207) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihcra) go to 202
      if (iquem.eq.ihcar) go to 202
      enddo
 202  rewind (11)

c
c     RA imagens ESO
c


      do i=1,11
      if (ler(i).ne.ibrac) go to 2999
      enddo

      rewind (12)

      write (12,*) ler
      rewind (12)
      read (12,*) rarara
      rewind (12)

c     rarara=rarara/15.d0
      iah=rarara
      iam=(rarara-iah)*60.d0
      sa=((rarara-iah)*60.d0-iam)*60.d0

      go to 206

c
c     Demais imagens, LNA, Roenos, etc
c

 2999 do i=2,19
      if (ler(i).ne.ibrac) go to 1001
      enddo

      go to 207

c

 1001 do i=1,70
      if (ler(i).eq.iplic) go to 203
      enddo

 203  rewind (12)

      i1=i

      ipu=0
      if (ler(i1+1).eq.ibrac) ipu=1

      do i=i1+1+ipu,70
      if (ler(i).eq.idois) go to 400
      if (ler(i).eq.ibrac) go to 400
      enddo

 400  i2=i

      do i=i2+2,70
      if (ler(i).eq.idois) go to 401
      if (ler(i).eq.ibrac) go to 401
      enddo

 401  i3=i

      do i=i3+2,70
      if (ler(i).eq.iplic) go to 402
      enddo

 402  i4=i


      write (12,*) (ler(j),j=i1+1,i2-1),' ',(ler(k),k=i2+1,i3-1),
     ?' ',(ler(l),l=i3+1,i4-1)
      rewind (12)

      read (12,*,err=207) iah,iam,sa


c     write (12,204) (ler(j),j=i+1,i+2),(ler(k),k=i+4,i+5),
c    ?(ler(l),l=i+7,i+11)
c204  format(2a1,2a1,5a1)
c     rewind (12)
c     read (12,205,err=315) iah,iam,sa
c205  format(2i2,f5.2)

c     rewind (12)

c     go to 206




      rewind (12)

      go to 206

c

c315  rewind (12)

c     read (12,320) iah,iam,isa
c320  format(3i2)
      sa=isa

c     rewind (12)

c     go to 206

c

 207  rewind (11)
      rewind (12)

      iah=99
      iam=99
      sa=99.999d0

c


 206  continue



c
c     Pega declinacao do header
c


      isig=mais

      do i=1,10000000
      read (11,200,end=219) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihcde) go to 212
      enddo
 212  rewind (11)

c
c     Imagens ESO
c


      do i=1,11
      if (ler(i).ne.ibrac) go to 2998
      enddo

      rewind (12)

      write (12,*) ler
      rewind (12)
      read (12,*) dedede
      rewind (12)

      isig=mais
      if (dedede.lt.0.d0) isig=menos
      dedede=dabs(dedede)

      idg=dedede
      idm=(dedede-idg)*60.d0
      ds=((dedede-idg)*60.d0-idm)*60.d0

      go to 220



c
c      Demais imagens LNA, Romenos, etc ...
c



 2998 do i=2,19
      if (ler(i).ne.ibrac) go to 1002
      enddo

      go to 219

c


 1002 do i=1,10

      if (ler(i).eq.menos) then
      isig=menos
      ler(i)=ibrac
      go to 411
      endif

      if (ler(i).eq.mais) then
      isig=mais
      ler(i)=ibrac
      go to 411
      endif

      enddo

c

 411  do i=1,70
      if (ler(i).eq.iplic) go to 412
      enddo

 412  i1=i

      ipu=0
      if (ler(i1+1).eq.ibrac) ipu=1

      do i=i1+1+ipu,70
      if (ler(i).eq.idois) go to 415
      if (ler(i).eq.ibrac) go to 415
      enddo

 415  i2=i

      do i=i2+2,70
      if (ler(i).eq.idois) go to 417
      if (ler(i).eq.ibrac) go to 417
      enddo

 417  i3=i

      do i=i3+2,70
      if (ler(i).eq.iplic) go to 418
      enddo

 418  i4=i


      write (12,*) (ler(j),j=i1+1,i2-1),' ',(ler(k),k=i2+1,i3-1),
     ?' ',(ler(l),l=i3+1,i4-1)
      rewind (12)

      read (12,*,err=219) idg,idm,ds

      rewind (12)

      go to 220


c     do i=1,70
c     if (ler(i).eq.iplic) go to 213
c     enddo

c213  i1=i


c     do i=i1+1,70
c     if ((ler(i).eq.idois).or.(ler(i).eq.ibrac)) go to 214
c     enddo

c214  i2=i


c     do i=i2+1,70
c     if ((ler(i).eq.idois).or.(ler(i).eq.ibrac)) go to 215
c     enddo

c215  i3=i

c     do i=i1+1,i2-1
c     if (ler(i).eq.menos) ler(i)='0'
c     if (ler(i).eq.mais ) ler(i)='0'
c     enddo

c     icasas=i2-i1-1

c     do j=i1+1,i3+5
c     if (ler(j).eq.ibrac) ler(j)='0'
c     enddo



c     if (icasas.eq.1) write (12,216) ler(i1+1),
c    ?(ler(j),j=i2+1,i3-1),(ler(k),k=i3+1,i3+5)
c     if (icasas.eq.2) write (12,217) (ler(l),l=i1+1,i1+2),
c    ?(ler(j),j=i2+1,i3-1),(ler(k),k=i3+1,i3+5)
c     if (icasas.eq.3) write (12,217) (ler(l),l=i1+2,i1+3),
c    ?(ler(j),j=i2+1,i3-1),(ler(k),k=i3+1,i3+5)

c216  format('0',a1,2a1,5a1)
c217  format(2a1,2a1,5a1)

c     rewind (12)

c     read (12,218,err=325) idg,idm,ds
c218  format(2i2,f5.2)

c     rewind (12)


c     go to 220

c

c325  rewind (12)

c     read (12,330) idg,idm,ids
c330  format(3i2)
c     ds=ids

c     rewind (12)

c     go to 220

c

 219  rewind (11)
      rewind (12)

      isig='+'
      idg=99
      idm=99
      ds=99.999d0

c

 220  continue


c
c     Pega tempo universal UT do header
c


      do i=1,10000000
      read (11,200,end=221) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihut) go to 222
      if (iquem.eq.ihutt) go to 500
      if (iquem.eq.ihuts) go to 3010
      if (iquem.eq.ihutm) go to 3010
      enddo

c
c     Imagens ESO
c

 3010 rewind (11)
      rewind (12)
      write (12,*) ler
      rewind (12)
      read (12,*) sut
      rewind (12)


      hora=sut/3600.d0

      if (iquem.eq.ihutm) hora=24.d0*sut/86400.d0

      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      go to 510

c
c     Aqui,nem UT nem EXPTIME sao dados; LNA dah Tini e Tfinal
c     em TL. Calculamos aqui UT inicial e EXPTIME
c

 221  rewind (11)



      do i=1,10000000
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihutb) go to 335
      enddo

 335  do i=1,70
      if (ler(i).eq.idois) go to 336
      enddo

 336  i1=i-3

      write (12,337) (ler(j),j=i1+1,i1+2),(ler(i),i=i1+4,i1+5),
     ?(ler(k),k=i1+7,i1+8),(ler(l),l=i1+10,i1+12)
 337  format(2a1,2a1,2a1,'.',3a1)

      read (11,200) iquem,(ler(j),j=1,70)
      
      write (12,337) (ler(j),j=i1+1,i1+2),(ler(i),i=i1+4,i1+5),
     ?(ler(k),k=i1+7,i1+8),(ler(l),l=i1+10,i1+12)

      rewind (11)
      rewind (12)

      read (12,340) iuth,iutm,sut
 340  format(i2,i2,f6.3)

      iuth=iuth+3

      read (12,340) juth,jutm,sutj

      rewind (12)

      juth=juth+3

      tfi=hmsgms(juth,jutm,sutj)
      tin=hmsgms(iuth,iutm,sut)

c
c     tfi depois da meia noite, tin antes da meia noite
c


      if (tfi.lt.tin) tfi=tfi+24.d0

c

      exps=(tfi-tin)*3600.d0
      iexps=exps+0.1


      go to 230      
c



 222  rewind (11)

      do i=1,70
      if (ler(i).eq.iplic) then
      ii=i
      go to 223
      endif

      if (ler(i).eq.idois) then
      ii=i-3
      go to 223
      endif

      enddo

 223  rewind (12)

      write (12,224) (ler(j),j=ii+1,ii+2),(ler(k),k=ii+4,ii+5),
     ?(ler(l),l=ii+7,ii+11)
 224  format(2a1,2a1,5a1)
      rewind (12)
      read (12,225) iuth,iutm,sut
 225  format(2i2,f5.2)

      rewind (12)

c

      go to 510

c

 500  rewind (11)

      do i=1,70
      if (ler(i).eq.iplic) then
      ii=i
      go to 503
      endif

      if (ler(i).eq.idois) then
      ii=i-3
      go to 503
      endif

      enddo

 503  rewind (12)

      write (12,504) (ler(j),j=ii+1,ii+2),(ler(k),k=ii+4,ii+5),
     ?(ler(l),l=ii+7,ii+8)
 504  format(2a1,2a1,2a1)
      rewind (12)
      read (12,505) iuth,iutm,isut
 505  format(3i2)

      rewind (12)

      sut=isut

 510  continue


c
c     Pega tempo de exposicao para calculo do instante medio UT
c



      do i=1,10000000
 227  read (11,1199,err=227,end=231) iquem,ichexp
1199  format(a9,1x,a70)

      if (iquem.eq.ihexp) go to 228
      if (iquem.eq.ihexpb) go to 228
      enddo

 231  iexps=0
      rewind (11)
      go to 230
c

 228  rewind (11)

      write (12,229) ichexp
 1229 format(a70)
 229  format(a20)
      rewind (12)
      read (12,*) exps
      iexps=exps

      rewind (12)


 230  continue

c
c     Calcula instante UT medio (UT pode ultrapassar 24hs)
c

      sut=sut+exps/2.d0

      if (sut.ge.60.d0) then
      idiv=sut/60.d0
      sut=sut-idiv*60.d0
      iutm=iutm+idiv
      endif

      if (iutm.ge.60) then
      idiv=iutm/60.d0
      iutm=iutm-idiv*60.d0
      iuth=iuth+idiv
      endif



c
c     Pega data gregoriana de Greenwhich do header
c


      do i=1,10000000
      read (11,200,err=232,end=232) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihdat) go to 232
      enddo
 232  rewind (11)

      do i=1,70
      if (ler(i).eq.iplic) then
      ii=i
      go to 233
      endif

      if (ler(i).eq.menos) then
      ii=i-5
      go to 233
      endif
      enddo

 233  rewind (12)

      write (12,234) (ler(j),j=ii+1,ii+10)
 234  format(10a1)
      rewind (12)
      read (12,235,err=300) iutano,iutmes,iutdia
 235  format(i4,1x,i2,1x,i2)

      rewind (12)


      go to 236

c
c     Outro tipo de formatacao de data
c

 300  rewind (12)


      read (12,310,err=311) iutdia,iutmes,iutano
 310  format(i2,1x,i2,1x,i4)

      rewind (12)
      go to 236

c

 311  rewind (12)

      read (12,312) iutdia,iutmes,iutano
 312  format(i2,1x,i2,1x,i2)

      rewind (12)

      if (iutano.lt.10) iutano=iutano+100

      iutano=iutano+1900


c

 236  continue




c
c     Pega filtro usado
c

      rewind (11)

      do i=1,10000000
      read (11,200,end=237) iquem,(ler(j),j=1,70)
      if (iquem.eq.ihfil)  go to 240
      if (iquem.eq.ihfilc) go to 240
      enddo

c
c     Header LNA incompleto, filtro nos COMMENTS
c 

 237  rewind (11)


      do i=1,10000000
      read (11,238,end=249) iqual,ifa
 238  format(a17,a5)
      if (iqual.eq.ihfilb) go to 239
      enddo

 239  rewind (11)

      ichfil=ifa

      go to 250

c

 249  ichfil='filtro desconhecido'

      rewind (11)

      go to 250

c

 240  rewind (11)


      do i=1,70
      if (ler(i).eq.iplic) go to 241
      enddo

 241  i1=i+1

      do i=70,1,-1
      if (ler(i).eq.iplic) go to 242
      enddo

 242  i2=i-1

      jj=i2-i1+1

      do j=1,jj
      ler(j)=ler(i1+j-1)
      enddo

      do i=jj+1,70
      ler(i)=ibrac
      enddo

      write (12,243) (ler(j),j=1,20)
 243  format(20a1)
      rewind (12) 

      read (12,244) ichfil
 244  format(a20)

      rewind (12)
c

 250  continue



c
c     Extracao de eventuais offsets em RA, DEC
c

      if (khead.eq.2) then
      ofra=0.d0
      ofde=0.d0
      go to 9999
      endif

c

      if (khead.eq.4) go to 4010
      if (khead.eq.3) go to 4000


c
c     Extracao de offsets em RA, DEC para CCDs individuais do mosaico
c     da ESO/WFI
c

 4000 continue

      rewind (11)

      do i=1,10000000
      read (11,4001) imagek
 4001 format(a9)
      if (imagek.eq.imagei) go to 4002 
      enddo
      
 4002 backspace (11)

      read (11,4003) imagek,idchip
 4003 format(a9,19x,i2)

      rewind (11)

c
c     Offset DEC
c

      if(idchip.eq.1 .or. idchip.eq.2 .or. idchip.eq.3 .or. idchip.eq.4)
     ? then

      ofde=+0.125d0

      else

      ofde=-0.125d0

      endif

c
c     Offset RA
c

      if (idchip.eq.2 .or. idchip.eq.6) then

      ofra=+0.0625d0

      go to 9999

      endif

c

      if (idchip.eq.3 .or. idchip.eq.7) then

      ofra=-0.0625d0

      go to 9999

      endif

c

      if (idchip.eq.1 .or. idchip.eq.5) then

      ofra=+0.1875d0

      go to 9999

      endif

c

      if (idchip.eq.4 .or. idchip.eq.8) then

      ofra=-0.1875d0

      go to 9999

      endif


c
c     Extracao de offsets em RA, DEC para CCDs individuais do SOAR/SOI
c

 4010 continue


      rewind (11)
      rewind (12)

      do i=1,10000000
      itrima=''
      read (11,4011) itrima
 4011 format(a80)

      if (itrima(1:9).eq.itrimi) go to 4012

      enddo

c

 4012 do i=10,80

      if (itrima(i:i).eq.'[') go to 4013
      enddo

 4013 i1=i+1


      do i=i1,80

      if (itrima(i:i).eq.':') go to 4014
      enddo

 4014 i2=i-1

      write (12,*) itrima(i1:i2)
      rewind (12)
      read (12,*) i
      rewind (12)

      rewind (11)



c
c     idchip=-1 CCD a esquerda
c     idchip=+1 CCD a direita
c

      if (i.lt.2049) then

      idchip=-1

      else

      idchip=+1

      endif

c
c     Extrai o angulo do rotator
c
      
      rewind (12)
      rewind (11)

      do i=1,10000000
      itrima=''
      read (11,4011,end=4020) itrima

      if (itrima(1:9).eq.'RAPANGL =') iran=itrima
      if (itrima(1:9).eq.'DECPANGL=') iden=itrima

      enddo

c

 4020 rewind (11)


      do i=10,80

      if (iran(i:i).eq.'/') go to 4023
      enddo

 4023 i1=i-1


      do i=10,80

      if (iden(i:i).eq.'/') go to 4024
      enddo

 4024 i2=i-1


      write (12,*) iran(10:i1),iden(10:i2)
      rewind (12)
      read (12,*) rangl,dengl
      rewind (12)
      rewind (11)


c
c
c     Orientacoes SOAR/SOI  (em 2007 usaram sistema indireto)
c
c
c                   RAPANGL =+90 graus (direto)
c       E <--       RAPANGL =-90 graus (indireto)
c            |      DECPANGL=+00 graus
c            |
c            v
c             S
c
c
c
c                         RAPANGL  =+00 graus            
c             --> S       DECPANGL =-90 graus (direto)
c            |            DECPANGL =+90 graus (indireto)
c            |
c            v
c             E
c
c
c  teta=0 para baixo, crescendo no sentido anti-horario, com ponto de
c  referencia alfa sendo "Leste" e ponto de referencia de delta sendo "Sul"
c
c

c
c     Determinando sistema direto/indireto
c
c     irot=+1  direto
c     irot=-1 indireto
c

      if (dengl.ge.rangl) then
      irot=-1
      rangl0=-90.d0
      else
      irot=+1
      rangl0=+90.d0
      endif

      teta=grarad*(rangl-rangl0)


c
c     Determinando offsets em RA e DE do CCD em relacao ao centro nominal
c
c
c     idchip=-1 CCD a esquerda
c     idchip=+1 CCD a direita
c

      soi=78.848d0/3600.d0

      ofra=idchip*irot*soi*dcos(teta)
      ofde=idchip*irot*soi*dsin(teta)


c
c     Fechando headers
c


 9999 continue

c
c     Calculo das datas julianas
c

      if (khead.ne.3) go to 4049

c
c     Calculo de instantes, data gregoriana e data juliana
c     para o caso da ESO/WFI
c


      rewind (11)
      rewind (12)

      do i=1,10000000
      itrima=''
      read (11,4011) itrima
      if (itrima(1:9).eq.'MJD-OBS =') go to 4030
      enddo

 4030 rewind (11)
      rewind (12)

      do i=10,80
      if (itrima(i:i).eq.'/') go to 4031
      enddo

 4031 i1=i-1

      write (12,*) itrima(10:i1)
      rewind (12)
      read (12,*) djm
      rewind (11)
      rewind (12)


      djm=djm+(exps/2.d0)/86400.d0

      dj=djm+dj1

      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,j)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0

      go to 4050
 
 
c
c     Demais Tels/instrumentos
c

 4049 continue

      
      fd=hmsgms(iuth,iutm,sut)/24.d0

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)


      djm=djm+fd

      dj=djm+djm0

c

 4050 continue

      close (11)
      close (12)

      call system(systa1)
      call system(systa2)


      return
      end





c
c     Subrotina jd2cal
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

