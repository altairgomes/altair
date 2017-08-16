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
c     In this version, the header extraction keywords are set by the user, so that
c     any header format, including non-standard FITS headers, can be read.
c
c
c     Last update:   M. Assafin - 15/Aug/2015
c
c




      IMPLICIT real*8 (A-H,O-Z)
      parameter(idi=100,idim=5001,ihead=1000)

      real*4    pixmat,mazi
      integer*2 bitpix,bitpyx

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
      character*50 infits
      character*96 iform

      character*8 kjud,kmjd,kdat,kbeg,kend,kbed,kexp,naxis1,naxis2,
     ?kaxes,kexten,kgain,knoise,kair,kbit,kscale,kbzero



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

c     open (1,file='PRAIA_photometry_19.dat')

      read (*,3) lista1
      read (*,*) imin,imax      


      read (*,3) lista2
      read (*,3) movie
      read (*,*) movnum


      read (*,4) kbit

      read (*,4) kexten

      read (*,4) kaxes 
      read (*,4) naxis1
      read (*,4) naxis2

      read (*,4) kscale
      read (*,4) kbzero


      read (*,4) kgain
      read (*,4) knoise
      read (*,4) kair


      read (*,4) kjud
      read (*,4) kmjd

      read (*,4) kdat

      read (*,4) kbeg
      read (*,4) kend
      read (*,4) kbed

      read (*,4) kexp
 4    format(a8)

      read (*,*) kh,km,sk    

      offtim=hmsgms(kh,km,sk)/24.d0

      read (*,7) maxes
 7    format(i2)

      read (*,7) mtim    



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
c     Setting instant of first frame for posterior target tracking
c

      do i=1,imin-1
      read (3,*)
      enddo

      infits=''

      read (3,3) infits
      rewind (3)



      call obhead (ihead,infits,tpose,kratio,ipflag,bitpyx,kexten,maxes,
     ?kaxes,naxis1,naxis2,kjud,kmjd,kdat,kbeg,kend,kbed,kexp,mtim,
     ?offtim,kgain,knoise,kair,kbit,kscale,kbzero,dj,exps,nx,ny,airmas,
     ?gain,rdnoise,bscale,bzero,bitpix,nheads)



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



      call obhead (ihead,infits,tpose,kratio,ipflag,bitpyx,kexten,maxes,
     ?kaxes,naxis1,naxis2,kjud,kmjd,kdat,kbeg,kend,kbed,kexp,mtim,
     ?offtim,kgain,knoise,kair,kbit,kscale,kbzero,dj,exps,nx,ny,airmas,
     ?gain,rdnoise,bscale,bzero,bitpix,nheads)



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

 

      call refits (idim,pixmat,infits,nx,ny,nheads,bscale,bzero,kswap,
     ?bitpix)

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
c
c     Subroutine obhead
c
c     Extracts info from the FITS header of the images
c
c
c     Last modification: M. Assafin - 15/Aug/2015
c
c
c



      subroutine obhead (ihead,infits,tpose,kratio,ipflag,bitpyx,kexten,
     ?maxes,kaxes,naxis1,naxis2,kjud,kmjd,kdat,kbeg,kend,kbed,kexp,mtim,
     ?offtim,kgain,knoise,kair,kbit,kscale,kbzero,dj,exps,nx,ny,airmas,
     ?gain,rdnoise,bscale,bzero,bitpix,nheads)



      IMPLICIT REAL *8 (A-H,O-Z)

      integer*2 bitpyx,bitpix

      character*2880 header
      character*(ihead*2880) head

      character*50 infits


      character*8  kjud,kmjd,kdat,kbeg,kend,kbed,kexp,naxis1,naxis2,
     ?kaxes,kexten,kgain,knoise,kair,kbit,kscale,kbzero



      character*80 itrima

c

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0
c


c
c     Initial values
c

c     pi    = 0.3141592653589793d1
c     grarad=pi/180.d0
c     radgra=180.d0/pi
c

      dj1=2400000.5d0

      ifim80=80

c



c
c     Opens FITS file
c


      open (1,file=infits,access='direct',form='unformatted',recl=2880)


c
c     Loads header
c

      mexten=0
      jhead=0



      do 20 i=1,ihead

      header='' 

      read (1,rec=i,err=30) header

c     write (*,*) 'i = ',i

      head(2880*(i-1)+1:2880*i)=header(1:2880)
 
c
c     Finds if there are header extensions
c
c     mexten = 0  -> no extension
c            = 1  -> there is an extension or extensions


      do 10 k = 1,2880,80

      itrima=''
      itrima(1:80)=header(k:k+79)

      do j=10,80
      if (itrima(j:j).eq.'/') go to 7
      enddo

 7    ifim=j-1

      if (itrima(1:8).eq.kexten) THEN

      do j=10,ifim
      if (itrima(j:j).eq.'t'.or.itrima(j:j).eq.'T') mexten=1
      enddo

      ENDIF



 10   continue


c

      do k = 1,2880,80
c     if (header(k:k+3).eq.'END ') jhead=jhead+1
      if (header(k:k+3).eq.'END ') jhead=i
      enddo



      if (mexten.eq.0.and.jhead.gt.0) go to 30
 
      if (mexten.eq.1.and.jhead.gt.1) go to 30



 
c
 20   continue

c

c30   nheads=i 

 30   nheads=jhead


      if (nheads.gt.ihead) then
      write (*,31) ihead
 31   format('Header size exceeded. More than ',i5,' pages. Exiting.')
      endif 

c     write (*,*) 'nheads = ',nheads


c
c     Extracts header information
c



      exps=-1.d0
      ihh2=-1

      nx=0
      ny=0


c
c     Loads header
c

      do 70 m=1,nheads

      header=''
      header(1:2880)=head(2880*(m-1)+1:2880*m)


      do 60 k = 1,2880,80

      itrima=''
      itrima(1:80)=header(k:k+79)


c
c     Checks out the position of the standard FITS header comment
c     character "/" in the string
c

      do i=10,80
      if (itrima(i:i).eq.'/') go to 33
      enddo

 33   ifim=i-1




c
c     Bitpix
c

      if (bitpyx.eq.-99) THEN

c     bitpix=16


      if (itrima(1:8).eq.kbit) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) bitpix


      ENDIF

      ENDIF




c
c     Pixel matrix dimensions
c

      if (maxes.eq.1) THEN

      if (itrima(1:8).eq.naxis1) THEN

      do i=10,ifim
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:ifim),*) nx

      ENDIF



      if (itrima(1:8).eq.naxis2) THEN


      do i=10,ifim
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:ifim),*) ny

      ENDIF


      ENDIF


      if (maxes.eq.2) THEN

      if (itrima(1:8).eq.kaxes) THEN

      do i=10,ifim
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      enddo

      read (itrima(10:ifim),*) n11,n12,n21,n22

      nx=n12-n11+1
      ny=n22-n21+1

      ENDIF


      ENDIF




c
c     Bscale and Bzero
c

      if (ipflag.eq.0) THEN

      bscale=1.d0
      bzero=0.d0


c
c     Bscale
c

      if (itrima(1:8).eq.kscale) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) bscale


      ENDIF



c
c     Bzero
c

      if (itrima(1:8).eq.kbzero) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) bzero


      ENDIF


      ENDIF




c
c     Gain and read noise from header
c

      if (kratio.eq.1) THEN

c
c     Gain 
c


      if (itrima(1:8).eq.kgain) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) gain


      ENDIF


c
c     Read noise
c


      if (itrima(1:8).eq.knoise) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) rdnoise


      ENDIF



      ENDIF



c
c     Air mass
c


      if (itrima(1:8).eq.kair) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) airmas


      ENDIF







c
c     Exposure time
c


      if (itrima(1:8).eq.kexp) THEN

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) exps

      ENDIF


c
c     Julian Date
c

      if (mtim.eq.3) THEN

      if (itrima(1:8).eq.kjud) then

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) dju

      endif

      ENDIF



c
c     Modified Julian Date
c

      if (mtim.eq.4) THEN

      if (itrima(1:8).eq.kmjd) then

      do i=10,ifim
      if (itrima(i:i).ne.'.'.and.itrima(i:i).ne.'E'.and.itrima(i:i).ne.
     ?'e'.and.itrima(i:i).ne.'D'.and.itrima(i:i).ne.'d'.and.itrima(i:i).
     ?ne.'+'.and.itrima(i:i).ne.'-') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) djum

      endif


      ENDIF


c
c     Date of observation
c

      if (mtim.eq.3) go to 50
      if (mtim.eq.4) go to 50


      if (itrima(1:8).eq.kdat) THEN

      do i=10,ifim80
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

c

      if (mtim.eq.1) then

      read (itrima(10:ifim80),*) ib1,ib2,ib3,ihh1,ihm1,hs1

      if (ib1.gt.ib3) then
      iutano=ib1
      iutmes=ib2
      iutdia=ib3
      else
      iutano=ib3
      iutmes=ib2
      iutdia=ib1
      endif


      else 


      read (itrima(10:ifim80),*) ib1,ib2,ib3

      if (ib1.gt.ib3) then
      iutano=ib1
      iutmes=ib2
      iutdia=ib3
      else
      iutano=ib3
      iutmes=ib2
      iutdia=ib1
      endif


      endif     


      ENDIF

c

 50   continue


c
c     Exposure start instant 
c


      if (itrima(1:8).eq.kbeg) THEN

      do i=09,ifim
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo


      iii=0
      do i=10,ifim
      if (itrima(i-1:i-1).eq.' '.and.itrima(i:i).ne.' ') iii=iii+1
      enddo


      if (iii.eq.3) then
      read (itrima(10:ifim),*) ihh1,ihm1,hs1
      endif

      if (iii.eq.2) then
      read (itrima(10:ifim),*) ihh1,hm1
      ihm1=hm1
      hs1=(hm1-ihm1)*60.d0
      endif


      if (iii.eq.1) then
      read (itrima(10:ifim),*) hh1

      hh1=hh1/86400.d0
      hh1=hh1*24.d0
      ihh1=hh1
      hm1=(hh1-ihh1)*60.d0
      ihm1=hm1
      hs1=(hm1-ihm1)*60.d0
      endif




      ENDIF



c
c     Exposure end instant 
c


      if (itrima(1:8).eq.kend) THEN

      do i=09,ifim
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo



      iii=0
      do i=10,ifim
      if (itrima(i-1:i-1).eq.' '.and.itrima(i:i).ne.' ') iii=iii+1
      enddo


      if (iii.eq.3) then
      read (itrima(10:ifim),*) ihh2,ihm2,hs2
      endif

      if (iii.eq.2) then
      read (itrima(10:ifim),*) ihh2,hm2
      ihm2=hm2
      hs2=(hm2-ihm2)*60.d0
      endif


      if (iii.eq.1) then
      read (itrima(10:ifim),*) hh2

      hh2=h21/86400.d0
      hh2=hh2*24.d0
      ihh2=hh2
      hm2=(hh2-ihh2)*60.d0
      ihm2=hm2
      hs2=(hm2-ihm2)*60.d0
      endif


      ENDIF




c
c     Start/END exposure instant at the same key 
c


      if (itrima(1:8).eq.kbed) THEN

      do i=09,ifim
      if (itrima(i:i).ne.'.') then
      icomp=ichar(itrima(i:i))
      if (icomp.lt.48 .or. icomp.gt.57) itrima(i:i)=' '
      endif
      enddo

      read (itrima(10:ifim),*) ihh1,ihm1,hs1,ihh2,ihm2,hs2

      ENDIF

c


 60   continue

 70   continue


c
c     Computes the mid-instant of the exposure (JD and UT date)
c


      if (mtim.eq.3) THEN

      djm=dju-dj1

      djm=djm+((tpose+exps)/2.d0)/86400.d0

      dj=djm+dj1

      go to 100

      ENDIF

c

      if (mtim.eq.4) THEN

      djm=djm+((tpose+exps)/2.d0)/86400.d0

      dj=djm+dj1

      go to 100

      ENDIF

c

      if (exps.lt.-0.5d0) THEN

      if (ihh2.lt.0) then

      exps=0.d0

      else

      fd1=hmsgms(ihh1,ihm1,hs1)
      fd2=hmsgms(ihh2,ihm2,hs2)

      if (fd2.lt.fd1) fd2=fd2+24.d0

      exps=dabs(fd2-fd1)*3600.d0

      endif


      ENDIF

c

      fd=hmsgms(ihh1,ihm1,hs1+(tpose+exps)/2.d0)
      
      fd=fd/24.d0



c
c     Julian Date 
c

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0



c
c     Date (Gregorian)
c


 100  continue

c
c     Applying time offset
c

      dj=dj+offtim
      djm=djm+offtim


      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,j)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0



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
c     Subroutine refits
c
c
c     Reads FITS image pixel matrix (integer, floating point).
c 
c
c      Last update:   29/Aug/2015 - M. Assafin
c
c


      subroutine refits (idim,pixmat,infits,nx,ny,nheads,bscale,bzero,
     ?kswap,bitpix)


      IMPLICIT REAL *8 (A-H,O-Z)

      integer*2 bitpix

      real*4 pixmat
      integer*2 iwork2
      integer*4 iwork4
      real*4 work4
      integer*8 iwork8
      real*8 work8

      integer*1 swork,iby8

      dimension iwork2(1440),swork(2880),iby8(8)
      dimension work4(720),iwork4(720)
      dimension work8(360),iwork8(360)

      dimension pixmat(idim,idim)



      character*50 infits
      character*50 erro
      character*20 imaux
      character*4 jmaux
      character*9 sista
      character*29 systa




c

      nbytes=2880

      if=1

c
c     Opens fits file
c


      open(if,file=infits,access='direct',form='unformatted',recl=2880)


c
c     Opens auxiliary file
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
c     Bit words
c
      i=bitpix
      if (i.lt.0) i=-i
      
      ibytes=i/8+0.1
      kwork=nbytes/ibytes+0.1




c
c     Reads pixel matrix
c


c
c     Checks byte-swap (litteendian or bigendian)
c
c     kswap ->  user key:
c
c          kswap = 0 automatic byte-swap determination
c          kswap = 1 no byte-swapping (user-defined)
c          kswap = 2 byte-swap pixel data (user-defined)
c
c     Automatic byte-swap determination (kswap=0):
c
c     iswap=1  no byte-swapping
c     iswap=2  byte-swap pixel data
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
c     Average of absolute values without swapping
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
c     Tests with swapping
c


      call swap (if,bitpix,ibytes,nbytes,irec,iwork2,iwork4,iwork8,
     ?work4,work8,swork)

c
c     Average of absolute values with swapping
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
c     Decides: swap or not?
c

      erro=''

      write (erro,*) c1

      do i=1,50
      if (ichar(erro(i:i)).ge.48 .and. ichar(erro(i:i)).le.57) go to 35
      enddo

      c1=1.d14

 35   continue

      if (c1.lt.1.d-10) c1=1.d14 


c

      erro=''
     
      write (erro,*) c2

      do i=1,50
      if (ichar(erro(i:i)).ge.48 .and. ichar(erro(i:i)).le.57) go to 40
      enddo

      c2=1.d14

 40   continue    

      if (c2.lt.1.d-10) c2=1.d14 


c

      if (c2.lt.c1) then
      iswap=2
      else
      iswap=1
      endif     


c     write (*,*) 'c1 c2 ',c1, c2

c
c     Reads matrix
c

 50   continue


c


      block=nx*ny*ibytes 
      block=block/nbytes 
      nblock=block
      iresto=(block-nblock)*nbytes
      iresto=iresto/ibytes+0.2


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
c     Last matrix block (if it exists)
c

      if (iresto.ne.0) THEN


      irec=irec+1


      if (iswap.eq.1) then

      if (bitpix.gt.0) then

      if (ibytes.eq.2) read (1,rec=irec) (iwork2(m),m=1,iresto)
      if (ibytes.eq.4) read (1,rec=irec) (iwork4(m),m=1,iresto)
      if (ibytes.eq.8) read (1,rec=irec) (iwork8(m),m=1,iresto)

      else

      if (ibytes.eq.4) read (1,rec=irec) (work4(m),m=1,iresto)
      if (ibytes.eq.8) read (1,rec=irec) (work8(m),m=1,iresto)
     
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

      ENDIF

c

      close (1)
      close (99)

c
c     Debug
c

c     read (*,*) j,i
c     j=100
c     i=100
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
c     subroutine swap
c
c
c     Swap dos bytes da imagem fits.
c
c     Imagem pode ser integer ou floating point
c
c 
c     Ultima modificacao: M. Assafin 15/08/2015
c
c


      subroutine swap(if,bitpix,ibytes,nbytes,irec,iwork2,iwork4,iwork8,
     ?work4,work8,swork)


      IMPLICIT REAL *8 (A-H,O-Z)



      dimension iwork2(1440),swork(2880),iby8(8)
      dimension work4(720),iwork4(720)
      dimension work8(360),iwork8(360)

      integer*2 bitpix

      integer*2 iwork2
      integer*4 iwork4
      real*4 work4
      integer*8 iwork8
      real*8 work8

      integer*1 swork,iby8

c

      read (if,rec=irec) swork


      do k=ibytes,nbytes,ibytes

      do m=1,ibytes
      iby8(m)=swork(k-m+1)
      enddo

      do m=1,ibytes
      swork(k-ibytes+m)=iby8(m)
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

