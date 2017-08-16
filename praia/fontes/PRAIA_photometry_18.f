c
c     Program PRAIA_photometry
c
c     Purpose
c
c
c     Computes photometric fluxes of targets and calibration objects.
C     This version holds for LNA, HP, ESO fits files.
C
c
c     Here box search is performed around guide star. User must give 
c     positions of calibration objects, targets and guide star (1rst frame
c     only). Guide star may or may not be a calibration object. Guide star
c     must be the brightest object whithin the box area furnished by the user.
c
c     Sky background is computed around each object (calibration obj or
c     target). The most typical sky background value is computed
c     on the pixels of the perimeter of the box region. The pixels associated
c     with the brightest and faintest quarters of counts are eliminated. 
c
c     This version is robust wrt calibration objects and targets follow-up,
c     even under severe telescope shifting and in the presence of clouded
c     (empty) fields intercalated with good observations - empty fields are
c     checked and skipped.
c
c     This version allows for fast moving calibration objects and targets,
c     such as asteroids and natural satellites. The apparent motion of the
c     moving object is computed for each image wrt the first image and the
c     moving object position keeps being continuously updated. This tracking
c     can be turned off for each object individually in the input file, if
c     necessary.
c
c     In this version, fits files (be integer or floating point, littleendian
c     or bigendian) are read entirely in fortran 77 (no external calls for
c     qfits C programs)
c
c
c     In this version, the seeing is computed for calibration objects and for
c     targets.
c
c     In this version, the S/N ratio for targets and calibration objects are
c     computed.
c
c     In this version, the flux area is circular with diameter given by the
c     user. The size of the box for centering purposes is independent and
c     also given by the user. The sky background is computed based on the
c     counts of the middle quarters of the sky background pixels, i.e.,
c     the 25% brighter and the 25% faintest pixels are excluded, and the
c     remaining pixels furnish the average and the dispersion (standard
c     deviation) of the sky background. The sky background pixels are
c     defined by the user, by rings of given radius and width.
c
c
c     In this version, the centroid is computed using a modified baricenter
c     model, i.e., only pixels with ADUs above a given threshold are used in
c     the centroid computations.
c
c
c     In this version, the user can fix the relative position of calibrators
c     and/or targets with respect to the guide object. In this case, no centroid
c     is computed for calibrators and/or targets - the center is derived from
c     fixed relative positions to the guide object. This is particularly useful
c     when the targets or calibrators are faint objects, or if their brightness
c     considerably drops, so that the computation of their centroids may no
c     longer work properly. This is overcome by using the known relative positions
c     to the guide object. These relative positions are computed from the first
c     frame furnished, for which all the objects must be properly exposed. 
c
c
c     In this version, the user can set the mid-instants of observations independently
c     of information from the header of the images. This is useful when the images
c     are generated in cubes and do not properly display the time and exposure
c     time in the headers.
c
c
c
c     Last update:   M. Assafin - 03/Feb/2013
c
c




      IMPLICIT real*8 (A-H,O-Z)
      parameter(idi=100,idim=5001)


      dimension pixmat(idim,idim)

      dimension jcen(idi),icen(idi)

      dimension xc(idi),yc(idi),cbx(idi),cby(idi),jtrack(idi),
     ?cpmx(idi),cpmy(idi),cf(idi),seec(idi),cdx(idi),cdy(idi),
     ?csnr(idi),cra(idi),cann(idi),cwi(idi),factoc(idi)

      dimension xt(idi),yt(idi),tbx(idi),tby(idi),itrack(idi),
     ?tpmx(idi),tpmy(idi),tf(idi),seet(idi),tdx(idi),tdy(idi),
     ?tsnr(idi),tra(idi),tann(idi),twi(idi),factot(idi)

      dimension xcal(idi),ycal(idi),xtar(idi),ytar(idi)

      character*50 lista1,lista2,ids9,ireg,movie
      character*150 infits
      character*20 ichobj,ichfil
      character*1  isig
      character*96 iform


      real*4    pixmat,mazi
      integer*2 bitpix,bitpyx

      COMMON /A14/IERRO

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

      snr(fobj,fsky,gain,rdnoise,npixel,nbin,darkcu,texp)=fobj*gain/
     ?dsqrt(fobj*gain+npixel*fsky*gain+npixel*rdnoise**2+(texp*darkcu*
     ?npixel*nbin**2)/3600.d0)

      snrbg(fobj,gain,npixel,nbg,sbg)=dsqrt(fobj/(1.d0/gain+
     ?(npixel*sbg**2)*(1.d0+1.d0/nbg)/fobj))

c
c     Auxiliary data
c


      ierro=0
c

      ic=idi
      it=idi

      do i=1,ic
      xc(i)=0.d0
      yc(i)=0.d0
      xcal(i)=0.d0
      ycal(i)=0.d0
      jtrack(i)=0
      cpmx(i)=0.d0
      cpmy(i)=0.d0
      cf(i)=0.d0
      seec(i)=0.d0
      cra(i)=0.d0
      cann(i)=0.d0
      cwi(i)=0.d0
      enddo

      do i=1,it
      xt(i)=0.d0
      yt(i)=0.d0
      itrack(i)=0
      xtar(i)=0.d0
      ytar(i)=0.d0
      tpmx(i)=0.d0
      tpmy(i)=0.d0
      tf(i)=0.d0
      seet(i)=0.d0
      tra(i)=0.d0
      tann(i)=0.d0
      twi(i)=0.d0
      enddo

c



c
c     Reads input file       
c

c     open (1,file='PRAIA_photometry_18.dat')

      read (*,3) lista1
      read (*,*) imin,imax      


      read (*,3) lista2
      read (*,3) movie
      read (*,*) movnum


      read (*,*) ihefit
      read (*,*) ipflag
      read (*,*) bscale
      read (*,*) bzero
      read (*,*) bitpyx
      read (*,*) kswap
      read (*,*) tpose
      read (*,*) timeof


      read (*,*) ktim
      read (*,*) mano,mmes,mdia,mho,mmi,sm
      read (*,*) expom
      read (*,*) cicle

   
      fd=hmsgms(mho,mmi,sm)
      if (fd.gt.24.d0) then
      fd=fd-24.d0
      mdia=mdia+1
      endif
      fd=fd/24.d0
      call iau_CAL2JD (mano,mmes,mdia,djm0,djm,iflag)
      djm=djm+fd
      dj=djm+djm0
      dstart=dj

      expom=expom/2.d0
      expom=expom/86400.d0
      cicle=cicle/86400.d0


      read (*,*) scala
      read (*,*) mazi
      read (*,*) vmin

      mazi=mazi-vmin+1.d0

      read (*,*) ksky

      read (*,*) darkcu
      read (*,*) nbin


      read (*,*) kratio

      read (*,*) gain
      read (*,*) rdnoise


      read (*,*) icx1,icx2,icy1,icy2
      read (*,*) factog
      read (*,*) gbx,gra,gann,gwi

      gby=gbx

      read (*,*) xg,yg
      read (*,*) kcen

      gbx=gbx/2.d0
      gby=gby/2.d0

      gra=gra/2.d0
      gann=gann/2.d0


      read (*,*) ic

      do i=1,ic

      read (*,*) jtrack(i)
      read (*,*) factoc(i)
      read (*,*) cbx(i),cra(i),cann(i),cwi(i)
      read (*,*) xc(i),yc(i)
      read (*,*) jcen(i)

      cbx(i)=cbx(i)/2.d0
      cby(i)=cbx(i)

      cra(i)=cra(i)/2.d0
      cann(i)=cann(i)/2.d0

      enddo


      read (*,*) it


      do i=1,it

      read (*,*) itrack(i)
      read (*,*) factot(i)
      read (*,*) tbx(i),tra(i),tann(i),twi(i)
      read (*,*) xt(i),yt(i)
      read (*,*) icen(i)

      tbx(i)=tbx(i)/2.d0
      tby(i)=tbx(i)

      tra(i)=tra(i)/2.d0
      tann(i)=tann(i)/2.d0

      enddo


c     close (1)




      if (ihefit.ne.1 .and. irefit.ne.2) irefit=1

      kah=0
      kam=0
      ak=timeof
      timeof=hmsgms(kah,kam,ak)/24.d0

c
c     Records relative positions of calibration stars wrt guide star 
c     for frame 1
c

      do i=1,ic
      xcal(i)=xc(i)-xg
      ycal(i)=yc(i)-yg
      enddo



c
c     Records relative positions of targets wrt guide star 
c     for frame 1

      do i=1,it
      xtar(i)=xt(i)-xg
      ytar(i)=yt(i)-yg
      enddo


c

      call system ('clear')

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'Starting image proccessing')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
c
      write (*,2)
 2    format (15x,'List of fits images to proccess -> ',$)
      write(*,3) lista1
 3    format(a50)
c
      write (*,5)
 5    format (15x,'Output light curve file         -> ',$)
      write(*,3) lista2

c

      write (*,*)
      write (*,*)

c
c     Reads images
c

      open (3,file=lista1)

c
c     Checks number of images
c

      do i=1,1000000
      read (3,*,end=10)
      enddo

 10   rewind (3)

      ijk=i-1

      if (imin.gt.ijk) imin=ijk
      if (imax.gt.ijk) imax=ijk

      if (imin.eq.0 .and. imax.eq.0) then
      imin=1
      imax=ijk
      endif

      nfiles=imax-imin+1

c
c     Correct number of ds9 movie frames
c

      if (movnum.gt.nfiles) movnum=nfiles

      if (movnum.eq.0) movnum=nfiles

      movnum=nfiles/movnum
c


      open (7,file=lista2)

c
c     Fixing instant of first frame for posterior target tracking
c

      do i=1,imin-1
      read (3,*)
      enddo

      infits=''

      read (3,3) infits
      rewind (3)



      if (ihefit.eq.1.or.ihefit.eq.5.or.ihefit.eq.6.or.ihefit.eq.7
     ?.or.ihefit.eq.8.or.ihefit.eq.9) then

      call obhead (tpose,infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,exps,nx,ny,
     ?airmas,kratio,gain,rdnoise,ihefit)

      else

      call obhaad (tpose,infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,exps,nx,ny,ofra,
     ?ofde,ihefit,airmas,kratio,gain,rdnoise,ihefit)

      endif


c
c     Recording time for frame 1
c

      if (ktim.eq.0) then

      djold=dj+timeof

      else

      djold=dstart+(i-1)*cicle+expom+timeof


      endif


c
c

      do i=1,imin-1
      read (3,*)
      enddo

c

      iccc=0


c

      open (77,file=movie)


c
c     Format for output file
c


      iform='(2(1x,i3),1x,f16.8,1x,f5.3,000(1x,f16.5,1x,f10.4,1x,f6.3),0
     ?00(1x,f16.5,1x,f10.4,1x,f6.3),1x,a50)'

      write (iform(28:30),'(i3.3)') it

      write (iform(59:61),'(i3.3)') ic


c
c     Looping of frames to proccess
c


      do 100 i=imin,imax

      call system ('clear') 

      read (3,3) infits

c
c     Defining output region file name for ds9 analysis
c

      ireg=''
      ireg=infits

      do mm=50,5,-1
      if (infits(mm-4:mm).eq.'.fits') then
      kk=mm-4
      go to 25
      endif
      enddo

      do mm=50,5,-1
      if (infits(mm:mm).ne.' ') then
      kk=mm+1
      go to 25
      endif
      enddo


 25   ireg(kk:kk+4)='.reg '


c
c     Writing batch movie for later ds9 region analysis
c

      if (i.eq.imin) then

      write (77,26) infits,ireg
 26   format('ds9 ',a50,' -region ',a50,' ',$)

      else 

      modo=mod(i,movnum)

      if (modo.eq.0) write (77,27) infits,ireg
 27   format(a50,' -region ',a50,' ',$)
 
      endif

c
c     open region file of frame 
c

      open (50,file=ireg)

c

      iccc=iccc+1

      write (*,28) iccc,nfiles,infits
 28   format(1x,'Proccessing image ',i5,' of ',i5,1x,': ',a50)

c
c     Extracting header information
c



      if (ihefit.eq.1.or.ihefit.eq.5.or.ihefit.eq.6.or.ihefit.eq.7
     ?.or.ihefit.eq.8.or.ihefit.eq.9) then

      call obhead (tpose,infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,exps,nx,ny,
     ?airmas,kratio,gain,rdnoise,ihefit)

      else

      call obhaad (tpose,infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,exps,nx,ny,ofra,
     ?ofde,ihefit,airmas,kratio,gain,rdnoise,ihefit)

      endif

c
c     Mid-instant of exposure
c

      if (ktim.eq.0) then

      dj=dj+timeof

      else

      dj=dstart+(i-1)*cicle+expom+timeof

      endif


c
c     Extracting pixel matrix
c

      bitpix=bitpyx

      call refits (idim,pixmat,infits,nx,ny,nheads,ichobj,ipflag,bscale,
     ?bzero,kswap,iswap,bitpix)

c
c     Fix negative and saturated pixel counts
c

      do k=1,ny
      do j=1,nx
      if (pixmat(j,k).le.vmin) pixmat(j,k)=mazi+vmin
      if (pixmat(j,k).lt.mazi) pixmat(j,k)=pixmat(j,k)-vmin+1.d0
      enddo
      enddo


c
c     Find guide star within given box area search
c

      cmax=-1.d14
      soma=0.d0


      lxg=1
      lyg=1

      jcx1=icx1+lxg
      jcx2=icx2-lxg

      jcy1=icy1+lyg
      jcy2=icy2-lyg

      if (jcx1.lt.01) jcx1=01
      if (jcx2.gt.nx) jcx2=nx
      if (jcy1.lt.01) jcy1=01
      if (jcy2.gt.ny) jcy2=ny

      if (jcx2.lt.01) jcx2=01
      if (jcx1.gt.nx) jcx1=nx
      if (jcy2.lt.01) jcy2=01
      if (jcy1.gt.ny) jcy1=ny


c
c     Computes overall sky background
c

      kkx1=jcx1
      kkx2=jcx2
      kky1=jcy1
      kky2=jcy2


      call skybac (idim,pixmat,mazi,nx,ny,kkx1,kkx2,kky1,kky2,sper,
     ?sper2)

      fundo=sper
      fundo2=sper2

      ceu=sper+factog*sper2


c
c     Finds the brightest object (the guide object) within the given search area
c

      do k=jcy1,jcy2
      do j=jcx1,jcx2

c     call box  (idim,pixmat,mazi,nx,ny,j,k,lxg,lyg,jx1,jx2,jy1,jy2,xgg,
c    ?ygg,ceu)

      xgg=j
      ygg=k

      call flux (idim,pixmat,mazi,nx,ny,xgg,ygg,lxg,soma,n)

c     write (*,*) 'j,k,soma = ',j,k,soma

      if (soma.gt.cmax) then
      cmax=soma
      xg=xgg
      yg=ygg
      endif

      enddo
      enddo


c
c     Refines the  guide star center
c

      lxg=gbx
      lyg=gby

      ixr=xg
      iyr=yg

      call box(idim,pixmat,mazi,nx,ny,ixr,iyr,lxg,lyg,kx1,kx2,ky1,ky2,
     ?xgg,ygg,ceu)

      xg=xgg
      yg=ygg


      call skycic (idim,pixmat,mazi,nx,ny,xg,yg,gann,gwi,sper,sper2,nbg)

      thres=sper+factog*sper2


c

      if (kcen.eq.1) then

      ixr=xg
      iyr=yg

      call box(idim,pixmat,mazi,nx,ny,ixr,iyr,lxg,lyg,kx1,kx2,ky1,ky2,
     ?xgg,ygg,thres)

      xg=xgg
      yg=ygg


      else

      fc=sper
      r=gra
      x=xg
      y=yg
 
      call gcc (idim,pixmat,kx1,kx2,ky1,ky2,r,mazi,x,y,ex,ey,iter,sig,
     ?hh,fc)
 
      xg=x
      yg=y


      endif



c
c     Checks guide star flux against clouds or empty field
c

c     iraio=1
c     call flux (idim,pixmat,mazi,nx,ny,xg,yg,iraio,pc,n)

c     pc=pc/n

      iiiix=xg
      iiiiy=yg

      pc=pixmat(iiiix,iiiiy)

      call skycic (idim,pixmat,mazi,nx,ny,xg,yg,gann,gwi,sper,sper2,nbg)

      thres=sper+factog*sper2


c     write (*,*) 'x,y,pc,thres,sper,sper2 = ',iiiix,iiiiy,pc,thres,
c    ?sper,sper2
c     stop

      if (pc.lt.thres) go to 100

c


      write (50,29) xg,yg,2*gbx+1.d0,2*gby+1.d0


 29   format('image; box(',4(f16.8,','),'0) # color=yellow')


c
c     Find calibration objects
c

      do 45 m=1,ic


c
c     Calibration object tracking
c

      if (jtrack(m).eq.-1) then

      xcc=xg+xcal(m)
      ycc=yg+ycal(m)

      kx1=xcc-cbx(m)
      kx2=xcc+cbx(m)

      ky1=ycc-cby(m)
      ky2=ycc+cby(m)

      if (kx1.lt.1)  kx1=1
      if (kx2.gt.nx) kx2=nx
      if (ky1.lt.1)  ky1=1
      if (ky2.gt.ny) ky2=ny

      go to 30

      endif

c

      if (jtrack(m).eq.0) then

      xc(m)=xg+xcal(m)
      yc(m)=yg+ycal(m)

      endif

c
      if (jtrack(m).eq.1) then

      xc(m)=xg+xcal(m)+cpmx(m)*(dj-djold)
      yc(m)=yg+ycal(m)+cpmy(m)*(dj-djold)

      endif


c

      lx=cbx(m)
      ly=cby(m)


      ixr=xc(m)
      iyr=yc(m)


      call box (idim,pixmat,mazi,nx,ny,ixr,iyr,lx,ly,kx1,kx2,ky1,ky2,
     ?xcc,ycc,ceu)

      call skycic (idim,pixmat,mazi,nx,ny,xcc,ycc,cann(m),cwi(m),sper,
     ?sper2,nbg)

      thres=sper+factoc(m)*sper2

c

      if (jcen(m).eq.1) then

      ixr=xcc
      iyr=ycc

      call box  (idim,pixmat,mazi,nx,ny,ixr,iyr,lx,ly,kx1,kx2,ky1,ky2,
     ?xcc,ycc,thres)


      else

      fc=sper
      r=cra(m)
      x=xcc
      y=ycc
 
      call gcc (idim,pixmat,kx1,kx2,ky1,ky2,r,mazi,x,y,ex,ey,iter,sig,
     ?hh,fc)
 
      xcc=x
      ycc=y

      endif

c



 30   continue

      iraio=cra(m)

      call flux (idim,pixmat,mazi,nx,ny,xcc,ycc,iraio,soma,n)


      xc(m)=xcc
      yc(m)=ycc

c



      write (50,37) xcc,ycc,2*cbx(m)+1.d0,2*cby(m)+1.d0
 37   format('image; box(',4(f16.8,','),'0) # color=green')

      write (50,38)  xcc,ycc,cra(m)+0.5d0
 38   format('image; circle(',2(f16.8,','),f16.8,') # color=red')

      write (50,39)  xcc,ycc,cann(m)+0.5d0
 39   format('image; circle(',2(f16.8,','),f16.8,') # color=blue')

      write (50,39)  xcc,ycc,cann(m)+0.5d0+cwi(m)




c
c     Computes sky background for calibration object and records calibration
c     object flux within the circular aperture area. 
c     The calibration flux is corrected from sky background.
c

      if (cann(m).ge.0.d0) then 
      call skycic (idim,pixmat,mazi,nx,ny,xc(m),yc(m),cann(m),cwi(m),
     ?sper,sper2,nbg)
      cf(m)=soma-n*sper
      else
      sper=fundo
      sper2=fundo2
      cf(m)=soma-n*sper
      endif

c
c     Computes S/N ratio for calibration object
c
c
c     ksky=1 -> classical general formula
c     ksky=2 -> sky background variance general formula


      if (ksky.eq.1) then

      csnr(m)=snr(cf(m),sper,gain,rdnoise,n,nbin,darcu,exps)

      endif

c

      if (ksky.eq.2) then

      csnr(m)=snrbg(cf(m),gain,n,nbg,sper2)

      endif

c
c     Computes the FWHM (seeing) of the calibration object
c


      r=cra(m)
      x=xc(m)
      y=yc(m)
      fc=sper


      call gcc (idim,pixmat,kx1,kx2,ky1,ky2,r,mazi,x,y,ex,ey,iter,sig,
     ?hh,fc)


      fwhm=2.d0*sig*1.177410023d0

      seec(m)=fwhm*scala

      if (seec(m).gt.99.999d0) seec(m)=99.999d0
      if (seec(m).lt.0.d0) seec(m)=99.999d0
      if (ierro.eq.1) seec(m)=99.999d0
      ierro=0

c

 45   continue


c
c     Find targets
c

      do 50 m=1,it

c
c     Target tracking
c


      if (itrack(m).eq.-1) then

      xtt=xg+xtar(m)
      ytt=yg+ytar(m)

      kx1=xtt-tbx(m)
      kx2=xtt+tbx(m)

      ky1=ytt-tby(m)
      ky2=ytt+tby(m)

      if (kx1.lt.1)  kx1=1
      if (kx2.gt.nx) kx2=nx
      if (ky1.lt.1)  ky1=1
      if (ky2.gt.ny) ky2=ny

      go to 47

      endif

c

      if (itrack(m).eq.0) then


      xt(m)=xg+xtar(m)
      yt(m)=yg+ytar(m)

      endif

c

      if (itrack(m).eq.1) then

      xt(m)=xg+xtar(m)+tpmx(m)*(dj-djold)
      yt(m)=yg+ytar(m)+tpmy(m)*(dj-djold)

      endif


c

      lx=tbx(m)
      ly=tby(m)


      ixr=xt(m)
      iyr=yt(m)

      call box (idim,pixmat,mazi,nx,ny,ixr,iyr,lx,ly,kx1,kx2,ky1,ky2,
     ?xtt,ytt,ceu)

      call skycic (idim,pixmat,mazi,nx,ny,xtt,ytt,tann(m),twi(m),sper,
     ?sper2,nbg)

      thres=sper+factot(m)*sper2

c
    
      if (icen(m).eq.1) then

      ixr=xtt
      iyr=ytt

      call box (idim,pixmat,mazi,nx,ny,ixr,iyr,lx,ly,kx1,kx2,ky1,ky2,
     ?xtt,ytt,thres)



      else



      fc=sper
      r=tra(m)
      x=xtt
      y=ytt
 
      call gcc (idim,pixmat,kx1,kx2,ky1,ky2,r,mazi,x,y,ex,ey,iter,sig,
     ?hh,fc)
 
      xtt=x
      ytt=y

      endif

c



 47   continue

      iraio=tra(m)

      call flux(idim,pixmat,mazi,nx,ny,xtt,ytt,iraio,soma,n)



      xt(m)=xtt
      yt(m)=ytt


c


      write (50,37) xtt,ytt,2*tbx(m)+1.d0,2*tby(m)+1.d0

      write (50,38) xtt,ytt,tra(m)+0.5d0

      write (50,39)  xtt,ytt,tann(m)+0.5d0

      write (50,39)  xtt,ytt,tann(m)+0.5d0+twi(m)



c
c     Computes sky background for targets and records target flux within the
c     circular aperture area. 
c     The target flux is corrected from sky background.
c

      if (tann(m).ge.0.d0) then 
      call skycic (idim,pixmat,mazi,nx,ny,xt(m),yt(m),tann(m),twi(m),
     ?sper,sper2,nbg)
      tf(m)=soma-n*sper
      else
      sper=fundo
      tf(m)=soma-n*sper
      endif

c
c     Computes S/N ratio for targets
c
      if (ksky.eq.1) then

      tsnr(m)=snr(tf(m),sper,gain,rdnoise,n,nbin,darcu,exps)

      endif

c

      if (ksky.eq.2) then

      tsnr(m)=snrbg(tf(m),gain,n,nbg,sper2)

      endif


c
c     Computes the FWHM (seeing) of the targets
c


      r=tra(m)
      x=xt(m)
      y=yt(m)
      fc=sper

      call gcc (idim,pixmat,kx1,kx2,ky1,ky2,r,mazi,x,y,ex,ey,iter,sig,
     ?hh,fc)

      fwhm=2.d0*sig*1.177410023d0

      seet(m)=fwhm*scala

      if (seet(m).gt.99.999d0) seet(m)=99.999d0
      if (seet(m).lt.0.d0) seet(m)=99.999d0
      if (ierro.eq.1) seet(m)=99.999d0
      ierro=0

c

 50   continue


      close (50)


c
c     Writing fluxes and instants to output file
c

c     write (7,iform) it,ic,dj,airmas,(tf(k),k=1,it),(tsnr(kk),kk=1,it),
c    ?(seet(kkk),kkk=1,it),(cf(l),l=1,ic),(csnr(ll),ll=1,ic),
c    ?(seec(lll),lll=1,ic),infits

      write (7,iform) it,ic,dj,airmas,(tf(k),tsnr(k),seet(k),k=1,it),
     ?(cf(l),csnr(l),seec(l),l=1,ic),infits



c
c     Update apparent (x,y) proper motion of calibration objects and
c     targets
c


      dift=dj-djold

c

      if (i.eq.imin) then

      do m=1,it
      xtar(m)=xt(m)-xg
      ytar(m)=yt(m)-yg
      enddo

      do m=1,ic
      xcal(m)=xc(m)-xg
      ycal(m)=yc(m)-yg
      enddo

      do m=1,it
      tpmx(m)=0.d0
      tpmy(m)=0.d0
      enddo

      do m=1,ic
      cpmx(m)=0.d0
      cpmy(m)=0.d0
      enddo

      go to 100

      endif


c


      do m=1,it

      if (itrack(m).ne.0) then
      tdx(m)=(xt(m)-xg)-xtar(m)
      tdy(m)=(yt(m)-yg)-ytar(m)
      else
      tdx(m)=0.d0
      tdy(m)=0.d0
      endif

      enddo


      do m=1,ic
      if (jtrack(m).ne.0) then
      cdx(m)=(xc(m)-xg)-xcal(m)
      cdy(m)=(yc(m)-yg)-ycal(m)
      else
      cdx(m)=0.d0
      cdy(m)=0.d0
      endif

      enddo

c

      if (dabs(dift).lt.1.d-7) then

      do m=1,it
      tpmx(m)=0.d0
      tpmy(m)=0.d0
      enddo

      do m=1,ic
      cpmx(m)=0.d0
      cpmy(m)=0.d0
      enddo

      else

      do m=1,it
      tpmx(m)=tdx(m)/dift
      tpmy(m)=tdy(m)/dift
      enddo

      do m=1,ic
      cpmx(m)=cdx(m)/dift
      cpmy(m)=cdy(m)/dift
      enddo

      endif

c

 100  continue

c
      close (3)
      close (7)
      close (77)

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

      subroutine obhead (tpose,infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,exps,nx,ny,
     ?airmas,kratio,gain,rdnoise,ihefit)
      IMPLICIT REAL *8 (A-H,O-Z)

      dimension header(2880),digit(4)
      character*3  tend
      character*1  digit
      character*1  header
      character*150 infits
      character*9  ihdat,ihut,ihuts,ihutt,ihcra,ihcar,ihcde,ihobj,ihfil,
     ?ihexp,iquem,ihutb,ihutc,ihfilc,ihexpb,is800
      character*20 ichobj,ichfil,ichexp
      character*1  ler(70),iplic,ibrac,idois,isig,mais,menos
      character*17 ihfilb,iqual
      character*5  ifa

      character*80 itrat

      character*9  imccet,imcceb,imccee,imccer,imcced,imccef

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
      rewind (12)

199   format(a9,1x,a20)
200   format(a9,1x,70a1)


c
c     LNA, new headers CCDs 101, 105, 106
c



      if (ihefit.eq.9) then

      rewind (11)

 6853 read (11,6661,err=6853,end=6856) itrat


      do i=1,80
      if (itrat(i:i).eq.',') itrat(i:i)='.'
      enddo


      if (itrat(1:9).eq.'EXPTIME =') then
      do i=1,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo
      itrat(1:9)='         '
      read (itrat,*) exps
      endif


      if (itrat(1:9).eq.'DATE-OBS=') then
      do i=1,80
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      enddo
      read (itrat,*) iutano,iutmes,iutdia
      endif



      if (itrat(1:9).eq.'UT      =') then

      do i=10,80
      if (itrat(i:i).eq.',') itrat(i:i)='.'
      enddo

      do i=10,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo

      itrat(1:9)='         '

      read (itrat,*) ia1,ia2,a3 

      endif

c

      go to 6853

 6856 continue



      sss=exps/2.d0+a3

      fd=hmsgms(ia1,ia2,sss)
      

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
c     S800
c

      if (ihefit.eq.6) then

      rewind (11)

 5553 read (11,6661,err=5553,end=5556) itrat


      do i=1,80
      if (itrat(i:i).eq.',') itrat(i:i)='.'
      enddo


      if (itrat(1:9).eq.'DATE-OBS=') then
      do i=1,80
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      enddo
      read (itrat,*) iutano,iutmes,iutdia
      endif



      if (itrat(1:9).eq.'DATE    =') then

      do i=10,80
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      enddo

      itrat(1:9)='         '

      read (itrat,*) ib1,ib2,ib3            

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


      if (itrat(1:9).eq.'TIME    =') then

      do i=10,80
      if (itrat(i:i).eq.',') itrat(i:i)='.'
      enddo

      do i=10,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo

      itrat(1:9)='         '

      read (itrat,*) ia1,ia2,a3,ia4,ia5,a6

      endif

c

      go to 5553

 5556 continue

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
c     IKON/LNA
c


      if (ihefit.eq.7) then


 6660 read (11,6661,err=6660) itrat
 6661 format(a80)

      do i=1,80
      if (itrat(i:i).eq.',') itrat(i:i)='.'
      enddo

c
c     Extrai data, hora e tempo de exposicao
c

      if (itrat(1:9).eq.'FRAME   =') then
      do i=1,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo
      read (itrat,*) iutano,iutmes,iutdia,ia1,ia2,a3
      go to 6662
      endif


      if (itrat(1:9).eq.'DATE-OBS=') then
      do i=1,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo
      read (itrat,*) iutano,iutmes,iutdia,ia1,ia2,a3
      go to 6662
      endif





      if (itrat(1:9).eq.'EXPOSURE=') then
      do i=1,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo
      itrat(1:9)='         '
      read (itrat,*) exps
      endif


      go to 6660

 6662 continue


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
c     Merlin-Raptor mCCD (B. Sicardy)
c


      if (ihefit.eq.8) then


 6760 read (11,6761,err=6760,end=6762) itrat
 6761 format(a80)

      do i=1,80
      if (itrat(i:i).eq.',') itrat(i:i)='.'
      enddo

c
c     Extrai dimensoes da matriz
c

      if (itrat(1:9).eq.'NAXIS1  =') then
      itrat(1:9)='         '
      read (itrat,*) nx
      endif


      if (itrat(1:9).eq.'NAXIS2  =') then
      itrat(1:9)='         '
      read (itrat,*) ny
      endif


c
c     Extrai data, hora e tempo de exposicao
c

      if (itrat(1:9).eq.'TIMESTMP=') then
      do i=1,80
      if (itrat(i:i).ne.'.') then
      icomp=ichar(itrat(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrat(i:i)=' '
      endif
      enddo
      read (itrat,*) iutdia,iutmes,iutano,ia1,ia2,a3
      endif


c
c     Extrai tempo de exposicao
c


      if (itrat(1:9).eq.'EXP_TIME=') then
      itrat(1:9)='         '
      read (itrat,*) exps
      exps=exps/1000.d0
      endif


      go to 6760

 6762 continue

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

      do i=1,70
      if (ler(i).eq.',') ler(i)='.'
      enddo

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

      do i=1,70
      if (ler(i).eq.',') ler(i)='.'
      enddo

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

      if (iquem.eq.ihdat .and. ihefit.eq.5) then

      rewind (12)

      do ii=1,70
      if (ler(i).eq.',') ler(i)='.'
      enddo

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
      if (ler(i).eq.',') ler(i)='.'
      enddo

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
c     Massa de ar
c

      rewind (11)
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200,end=3436) iquem,(ler(j),j=1,70)
      if (iquem.eq.'AIRMASS =') go to 3434
      enddo

 3434 rewind (11)

      do i=1,70
      if (ler(i).eq.',') ler(i)='.'
      enddo

      do i=1,70
      if (ler(i).eq.iplic) ler(i)=' '
      enddo

      write (12,*) (ler(j),j=1,22)
      rewind (12)
      read (12,*) airmas
      rewind (12)

      go to 3437

 3436 rewind (11)
      rewind (12)

      airmas=1.d0

 3437 continue

c
c     Ganho e read out noise
c

      if (kratio.eq.2) go to 3446

      rewind (11)
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'GAIN    =') go to 3444
      enddo

 3444 rewind (11)

      do i=1,70
      if (ler(i).eq.',') ler(i)='.'
      enddo

      write (12,*) (ler(j),j=1,70)
      rewind (12)
      read (12,*) gain   
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'RDNOISE =') go to 3445
      enddo

 3445 rewind (11)

      do i=1,70
      if (ler(i).eq.',') ler(i)='.'
      enddo

      write (12,*) (ler(j),j=1,70)
      rewind (12)
      read (12,*) rdnoise
      rewind (12)

 3446 continue

c

      close (11)
      close (12)


      call system(systa1)
      call system(systa2)



c
c     Calculo das datas julianas
c

      sut=sut+tpose/2.d0
      fd=hmsgms(iuth,iutm,sut)/24.d0

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0

c


      return
      end







c
c     Subrotina OBHAAD
c
C     Extrai dados do header de imagens fits
C

      subroutine obhaad (tpose,infits,ichobj,iah,iam,sa,isig,idg,idm,ds,
     ?iuth,iutm,sut,iutano,iutmes,iutdia,djm,dj,ichfil,exps,nx,ny,ofra,
     ?ofde,khead,airmas,kratio,gain,rdnoise,ihefit)

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

c
c
c     Ganho e read out noise ESO
c

      if (kratio.eq.2) go to 3446

      rewind (11)
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'OUTGAIN =') go to 3444
      enddo

 3444 rewind (11)

      write (12,*) (ler(j),j=1,12)
      rewind (12)
      read (12,*) gain   
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'OUTRON  =') go to 3445
      enddo

 3445 rewind (11)

      write (12,*) (ler(j),j=1,12)
      rewind (12)
      read (12,*) rdnoise
      rewind (12)

 3446 continue

      go to 4050
 
 
c
c     Demais Tels/instrumentos
c

 4049 continue

      
      sut=sut+tpose/2.d0
      fd=hmsgms(iuth,iutm,sut)/24.d0

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)


      djm=djm+fd

      dj=djm+djm0
c
c
c     Ganho e read out noise
c

      if (kratio.eq.2) go to 3456

      rewind (11)
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'GAIN    =') go to 3454
      enddo

 3454 rewind (11)

      write (12,*) (ler(j),j=1,70)
      rewind (12)
      read (12,*) gain   
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'RDNOISE =') go to 3455
      enddo

 3455 rewind (11)

      write (12,*) (ler(j),j=1,70)
      rewind (12)
      read (12,*) rdnoise
      rewind (12)

 3456 continue


c

 4050 continue

c
c     Massa de ar
c

      rewind (11)
      rewind (12)

      do i=1,10000000
      iquem=''
      read (11,200) iquem,(ler(j),j=1,70)
      if (iquem.eq.'AIRMASS =') go to 3434
      if (iquem.eq.'AIRMSTAR=') go to 3434
      enddo

 3434 rewind (11)

      write (12,*) (ler(j),j=1,22)
      rewind (12)
      read (12,*) airmas
      rewind (12)

c
      close (11)
      close (12)

      call system(systa1)
      call system(systa2)


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
      subroutine refits (idim,pixmat,infits,nx,ny,nheads,ichobj,ipflag,
     ?bscale,bzero,kswap,iswap,bitpix)

      IMPLICIT REAL *8 (A-H,O-Z)


      character*150 infits
      character*20 ichobj
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

      dimension pixmat(idim,idim)
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
      DIMENSION IORDEM(1000000),NVAL(1000000)
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
c     Media e desvio padrao
c
c
c     entrada:
c
c        xvam  = soma dos valores individuais
c        xvas  = soma dos quadrados dos valores individuais
c
c     saida:
c
c        xvam  = media dos valores individuais
c        xvas  = desvio padrao em relacao a media
c



      subroutine desvio (nest,xvam,xvas)

      IMPLICIT REAL *8 (A-H,O-Z)

c

      dnove=99.999d0

c

      EXMED=XVAM/NEST

c

      if (nest.eq.1) then

      xvas=dnove

      else

      raiz=XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2

      if (raiz.lt.0.d0) then

      xvas=0.d0

      go to 10

      else

      XVAS=DSQRT(raiz/(NEST-1.D0))

      endif
      endif
c

 10   XVAM=EXMED


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

      subroutine  GCC (idim,pixmat,nx,mx,ny,my,raio,maximo,bx,by,sigdex,
     ?sigdey,icontt,sigx,alturx,fundox)

      IMPLICIT REAL*8 (A-H,O-Z)
      real*4 pixmat(idim,idim),maximo,xmaxi
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
c     SIGX=0.D0
c     SIGY=0.D0

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

c     PARAM(3)=sigx

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
   14 CALL GAUSIC (idim,pixmat,ICONTT,RAIO,NX,NY,MX,MY,maximo,NTERMX,
     ?PARAM,DELTAX,XSIGMA,XLAMDA,RESIDX)

      if (ierro.eq.1) return

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
      SUBROUTINE GAUSIC (idim,pixmat,KEY,RAIO,NX,NY,MX,MY,maximo,NTERMS,
     ?A,DELTAA,SIGMAA,FLAMDA,CHISQR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(5),DELTAA(5),SIGMAA(5),B(5),ALPHA(21,21),BETA(21),
     ?DERIV(5),ARRAY(21,21)
      real*4 pixmat(idim,idim),maximo

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
      IF (IERRO.EQ.1) go to 107

 
c     FREE = (MX-NX+1)+(MY-NY+1) - NTERMS

      FREE=ICONT-NTERMS

C

      if (free.le.0.d0) go to 107

      IF (FREE.LE.0.D0) THEN
      IERRO=1
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
      IF (IERRO.EQ.1) go to 107

 
C
C     CALCULA GRAUS DE LIBERDADE FREE
C

c     COUNT=ICONT
c     FREE=2*DSQRT(COUNT)-NTERMS
      FREE=ICONT-NTERMS
C

      IF (FREE.LE.0.D0) THEN
      IERRO=1
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

   63 CHISQ1=QIQUAD(idim,pixmat,KEY,RAIO,NX,NY,MX,MY,FREE,A,maximo)

C
      if (ierro.eq.1) go to 107


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


   93 CHISQR = QIQUAD(idim,pixmat,KEY,RAIO,NX,NY,MX,MY,FREE,B,maximo)

C

      if (ierro.eq.1) go to 107

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
      DOUBLE PRECISION FUNCTION QIQUAD (idim,pixmat,KEY,RAIO,NX,NY,MX,
     ?MY,FREE,A,maximo)
      IMPLICIT REAL *8 (A-H,O-Z)
      real*4 pixmat(idim,idim),maximo

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






c
c
c    Subroutine box
c
c
c    Defines object box center for flux computations.
c
c    (ix,iy)         = input box center 
c    ladoxx,ladoyy   = half box sizes
c
c    kx1,kx2,ky1,ky2 = limits of used box
c    (x,y)           = computed center (baricenter)
c    s               = flux in ADUs, flagged pixels excluded
c    n               = number of actual pixels used on flux computation,
c                      flagged pixels excluded
c    cut             = ADU sky threshold cutoff for baricenter model
c
c
c
c     Last update: M. Assafin 28/April/2008
c

      subroutine box (idim,pixmat,mazi,nx,ny,ix,iy,ladoxx,ladoyy,kx1,
     ?kx2,ky1,ky2,x,y,cut)

      implicit real*8 (a-h,o-z)

      real*4 pixmat(idim,idim),mazi



      x=0.d0
      y=0.d0
      s=0.d0
      n=0

c

      kx1=ix-ladoxx
      kx2=ix+ladoxx
      ky1=iy-ladoyy
      ky2=iy+ladoyy

      if (kx1.lt.1) kx1=1
      if (kx2.gt.nx) kx2=nx
      if (ky1.lt.1) ky1=1
      if (ky2.gt.ny) ky2=ny

      if (kx2.lt.1) kx2=1
      if (kx1.gt.nx) kx1=nx
      if (ky2.lt.1) ky2=1
      if (ky1.gt.ny) ky1=ny

      do 2 k=ky1,ky2
      do 1 j=kx1,kx2

      
      if (pixmat(j,k).ge.(mazi-2.d0)) go to 1
      if (pixmat(j,k).lt.cut) go to 1
      s=s+pixmat(j,k)
      x=x+j*pixmat(j,k)
      y=y+k*pixmat(j,k)
      n=n+1

  1   continue
  2   continue

      if (n.eq.0) then
      x=ix
      y=iy
      return
      endif

      x=x/s
      y=y/s

      return
      end



c
c
c    Subroutine flux
c
c
c    Computes object flux by circular aperture photometry
c
c    (ix,iy)         = object center 
c    radius          = aperture radius 
c
c    kx1,kx2,ky1,ky2 = limits of box for flux computations
c
c    s               = flux in ADUs, flagged pixels excluded
c    n               = number of actual pixels used on flux computation,
c                      flagged pixels excluded
c
c
c    For flux computations, the pixel is only computed within a circular area,
c    with radius:
c
c     radius=iradius
c     ichave : positive -> pixel within circular area
c              negative -> pixel not within circular area
c
c
c     Last update: M. Assafin 28/April/2008
c

      subroutine flux (idim,pixmat,mazi,nx,ny,x,y,iradius,fluxo,nflux)

      implicit real*8 (a-h,o-z)

      real*4 pixmat(idim,idim),mazi


c

      nflux=0.d0
      fluxo=0.d0

c
c     Defines the radius of the circular aperture area
c

c

      ix=x
      iy=y
      ladoxx=iradius
      ladoyy=iradius
      radius=iradius

      jx1=ix-ladoxx
      jx2=ix+ladoxx
      jy1=iy-ladoyy
      jy2=iy+ladoyy

      if (jx1.lt.1)  jx1=1
      if (jx2.gt.nx) jx2=nx
      if (jy1.lt.1)  jy1=1
      if (jy2.gt.ny) jy2=ny

      if (jx2.lt.1)  jx2=1
      if (jx1.gt.nx) jx1=nx
      if (jy2.lt.1)  jy2=1
      if (jy1.gt.ny) jy1=ny

      do 2 k=jy1,jy2
      do 1 j=jx1,jx2

      if (pixmat(j,k).ge.(mazi-2.d0)) go to 1

      call circul (radius,x,y,j,k,ichave)

      if (ichave.lt.0) go to 1

      nflux=nflux+1
      fluxo=fluxo+pixmat(j,k)

  1   continue
  2   continue

c

      if (nflux.le.0.d0) then
      nflux=0
      fluxo=0.d0
      return
      endif


      return
      end






c
c
c    Subroutine skybac
c
c
c    Computes sky background from perimeter counts histogram 
c
c
c

      subroutine skybac (idim,pixmat,mazi,nx,ny,kkx1,kkx2,kky1,kky2,
     ?sper,sper2)

      implicit real*8 (a-h,o-z)

      real*4 pixmat(idim,idim),mazi

      dimension sv(1000000),ksv(1000000),ior(1000000)

      common /b1/sv,ksv,ior

c

      kx1=kkx1
      kx2=kkx2
      ky1=kky1
      ky2=kky2

      if (kx1.lt.1)  kx1=1
      if (kx2.gt.nx) kx2=nx
      if (ky1.lt.1)  ky1=1
      if (ky2.gt.ny) ky2=ny

c


      isv=1000000

c

      do iii=1,isv
      sv(iii)=0.d0
      ksv(iii)=0
      ior(iii)=iii
      enddo

      n=0

      do 30 j=kx1,kx2
      if (pixmat(j,ky1).ge.(mazi-2.d0)) go to 30
      n=n+1
      sv(n)=pixmat(j,ky1)
      ksv(n)=sv(n)*100000
 30   continue

      do 31 j=kx1,kx2
      if (pixmat(j,ky2).ge.(mazi-2.d0)) go to 31
      n=n+1
      sv(n)=pixmat(j,ky2)
      ksv(n)=sv(n)*100000
 31   continue

      do 32 k=ky1+1,ky2-1
      if (pixmat(kx1,k).ge.(mazi-2.d0)) go to 32
      n=n+1
      sv(n)=pixmat(kx1,k)
      ksv(n)=sv(n)*100000
 32   continue

      do 33k=ky1+1,ky2-1
      if (pixmat(kx2,k).ge.(mazi-2.d0)) go to 33
      n=n+1
      sv(n)=pixmat(kx2,k)
      ksv(n)=sv(n)*100000
 33   continue


c

      call ordem (isv,n,ior,ksv)

      quart=n/4.d0
      i1=quart
      i2=n-quart

      n=i2-i1+1
      sper=0.d0

      do iii=i1,i2
      kkk=ior(iii)
      sper=sper+sv(kkk)
      sper2=sper2+(sv(kkk))**2
      enddo

      call desvio (n,sper,sper2)

      return
      end





c
c
c    Subroutine skycic
c
c
c    Computes sky background from middle quarter sky pixels, excluding 25%
c    brighter and 25% fainter ones. Computed pixels are picked up from within
c    an anulus defined by the user:
c
c    radius = radius of anulus
c    width  = width of anulus 
c
c    n      = number of sky background pixels efectivelly used in the
c             computations  
c
c

      subroutine skycic (idim,pixmat,mazi,nx,ny,x,y,anulus,width,sper,
     ?sper2,n)

      implicit real*8 (a-h,o-z)

      real*4 pixmat(idim,idim),mazi

      dimension sv(1000000),ksv(1000000),ior(1000000)

      common /b1/sv,ksv,ior

      isv=1000000

c

      do iii=1,isv
      sv(iii)=0.d0
      ksv(iii)=0
      ior(iii)=iii
      enddo

      n=0

c
c     Defines the matrix search limits
c

      raiom=anulus+width

      ix=x
      iy=y
      ladoxx=raiom
      ladoyy=raiom

      jx1=ix-ladoxx
      jx2=ix+ladoxx
      jy1=iy-ladoyy
      jy2=iy+ladoyy

      if (jx1.lt.1)  jx1=1
      if (jx2.gt.nx) jx2=nx
      if (jy1.lt.1)  jy1=1
      if (jy2.gt.ny) jy2=ny

      if (jx2.lt.1)  jx2=1
      if (jx1.gt.nx) jx1=nx
      if (jy2.lt.1)  jy2=1
      if (jy1.gt.ny) jy1=ny

c

      do 2 k=jy1,jy2
      do 1 j=jx1,jx2

      if (pixmat(j,k).ge.(mazi-2.d0)) go to 1

      call circul (anulus,x,y,j,k,ichave)
      if (ichave.gt.0) go to 1

      call circul (raiom,x,y,j,k,ichave)
      if (ichave.lt.0) go to 1

      n=n+1
      sv(n)=pixmat(j,k)
      ksv(n)=sv(n)*100000

  1   continue
  2   continue

c

      call ordem (isv,n,ior,ksv)

      quart=n/4.d0
      i1=quart
      i2=n-quart

      n=i2-i1+1
      sper=0.d0

      do iii=i1,i2
      kkk=ior(iii)
      sper=sper+sv(kkk)
      sper2=sper2+(sv(kkk))**2
      enddo

      call desvio (n,sper,sper2)

      return
      end

