c
c
c     PRAIA_trajectory
c
c
c     Computes shadow path from a stellar occultation by the adjustment
c     of observations to the ephemeris, taking into consideration color
c     refraction. Both the target body and the star must be in the same
c     FOV. Observations must have been made in the same night, and must
c     be distributed in only one side of the meridian (East or West).
c
c     The (RA,Dec) reference star position is arbitrary, but we recommend
c     the use of averaged observed coordinates, preferably consistent
c     with the observations used in the adjustment.
c
c     The "TNO minus star" relative position differences are computed
c     from PRAIA output files, specifically from PRAIA offset file types
c     for the TNO and for the star. Only data collected from common frames
c     where both TNO and star are imaged toghether are considered.  
c 
c     C/A, P/A, central instant, shadow speed and their errors are
c     computed. The apparent trajectory in the (X,Y) sky plane over the
c     Earth figure is computed for plot purposes.
c
c
c     Useful for last-minute predictions and LC chord fittings.
c
c     
c     Last modification: M. Assafin  01/Jun/2012
c
c

      implicit double precision (a-h,o-z)

      parameter(stdin=5,stdout=6)


      dimension tra(100000),tde(100000),sra(100000),sde(100000),
     ?xest(100000),xesta(100000),xp(100000),yp(100000),coef(21),
     ?coefe(21),dj(100000),daju(100000),xtobs(100000),ytobs(100000),
     ?datae(100000),datao(100000),datat(100000),xgeocc(100000),
     ?ygeocc(100000),xcoef(21),xcoefe(21),ycoef(21),ycoefe(21)
     
      dimension dxob(100000),dyob(100000),dxep(100000),dyep(100000),
     ?dxocc(100000),dyocc(100000)

      dimension acomps(100000),dcomps(100000),acompt(100000),
     ?dcompt(100000),zen(100000),rex(100000),rey(100000),itira(100000)
   


      character*50 ctno,cstar,ptno,pearth,tfits(100000),sfits(100000),
     ?crel,traj,topobs,geoocc,ler,plot

      character*20 iobalv,ichfil

      character*1 isig

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

      dexy(xx,yy,z,w)=dsin(yy)*dsin(w)+dcos(yy)*dcos(w)*dcos(xx-z)
      xpad(xx,yy,z)=dcos(yy)*dsin(xx-z)
      ypad(xx,yy,z,w)=dsin(yy)*dcos(w)-dcos(yy)*dsin(w)*dcos(xx-z)


c
c     Auxiliary values
c    

      idimc=50

      idim=100000

      idimp=21

      dj1=2400000.5D0

      pi=0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi

      au=149597870.691d0

      tbox=1.d0/86400.d0
      tbox=0.1d0*tbox
 

c

 1    format(a50)

      read (*,1) ctno
      read (*,1) cstar
      read (*,1) geoocc 
      read (*,1) topobs
      read (*,1) crel
      read (*,1) traj
      read (*,1) plot
      read (*,1) ptno
      read (*,1) pearth

      read (*,*) sclip

      sclip=sclip/3600.d3

      read (*,1) ler

      read (*,*) along
      read (*,*) alati
      read (*,*) altit

      read (*,*) distan

      distan=distan*au

      read (*,*) rearth
      read (*,*) rtno

      read (*,*) xsize, ysize
      read (*,*) stepx

      xsize=xsize/2.d0
      ysize=ysize/2.d0

      read (*,*) rca
      read (*,*) rpa
      read (*,*) rvel
      read (*,*) iryy,irmo,irdd,irhh,irmm,rss


c
c     Computes reference Julian Date of C/A
c

      fd=hmsgms(irhh,irmm,rss)/24.d0

      call iau_CAL2JD (iryy,irmo,irdd,djm0,djm,iflago)

      djm=djm+fd
            
      rdate=djm+djm0



c
c     Stores (RA,Dec) reference star position
c

      is=+1
      do i=1,idimc
      if (ler(i:i).eq.'-') then
      is=-1            
      ler(i:i)=' '
      endif
      enddo

      read (ler,*) iah,iam,as,idg,idm,ds

      aref=15.d0*hmsgms(iah,iam,as)
      dref=is*hmsgms(idg,idm,ds)


c


 10   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))

c
c     Stores TNO observed data
c

      open(2,file=ctno)


      do i=1,idim

 15   tfits(i)=''

      read (2,10,end=20) dx,dy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,tra(i),tde(i),iuth,iutm,sut,
     ?iutano,iutmes,iutdia,datat(i),iexps,ichfil,tfits(i),iobalv,nx,ny

      enddo

 20   close (2)

      ntno=i-1

      

c
c     Stores star observed data
c

      open(2,file=cstar)

      
      do i=1,idim

 25   sfits(i)=''


      read (2,10,end=30) dx,dy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,sra(i),sde(i),iuth,iutm,sut,
     ?iutano,iutmes,iutdia,daju(i),iexps,ichfil,sfits(i),iobalv,nx,ny

      enddo

 30   close (2)

      nstar=i-1


      
c
c     Common fields. Compute TNO relative distances w.r.t. star. Computs
c     averaged observed (RAC,DEC) positions. 
c

      racsta=0.d0
      decsta=0.d0

      open (1,file=crel)

      write (1,*) 'TNO minus star: Dacosd (dg), Dd(dg), Julian Date'

      k=0

      do 41 i=1,ntno
      do 40 j=1,nstar

      dt=dabs(datat(i)-daju(j))

c     if (tfits(i).ne.sfits(j)) go to 40

      if (dt.gt.tbox) go to 40

      dx=15.d0*(tra(i)-sra(j))*dabs(dcos(grarad*sde(j)))
      dy=tde(i)-sde(j)


      k=k+1

      dxob(k)=dx
      dyob(k)=dy

      dj(k)=daju(j)

      racsta=racsta+15.d0*sra(j)
      decsta=decsta+sde(j)


c
c     Writes relative (RA,Dec)s "TNO minus star" and instants
c

      write (1,35) dx,dy,dj(k)
 35   format(2(1x,f15.10),1x,f16.8)



 40   continue

 41   continue

      nptos=k

      close (1)

      racsta=racsta/nptos
      decsta=decsta/nptos


c
c     Loads geocentric TNO ephemeris for the entire
c     time span of the occultation.
c

      open(1,file=geoocc)

      
      do i=1,idim

      read (1,42,end=43) iah,iam,as,isig,idg,idm,ds,datae(i),
     ?iobalv
 42   format(1x,i2,1x,i2,1x,f9.6,1x,a1,i2,1x,i2,1x,f8.5,1x,f16.8,
     ?1x,a20)

      rafat=15.d0*hmsgms(iah,iam,as)
      defat=hmsgms(idg,idm,ds)
      if (isig.eq.'-') defat=-defat

      xgeocc(i)=rafat
      ygeocc(i)=defat

      enddo

 43   close(1)

      ngeocc=i-1




c
c     Loads topocentric TNO ephemeris for the instants of observations.
c     Only ephemerides points corresponding to common TNO/star observations
c     are stored.
c

      open(1,file=topobs)


      k=0
      
      do i=1,idim

 44   read (1,42,end=46) iah,iam,as,isig,idg,idm,ds,datao(i),
     ?iobalv

      do kk=1,nptos
      aux=dabs(datao(i)-dj(kk))
      if (aux.lt.tbox) go to 45
      enddo

      go to 44

 45   continue

      rafat=15.d0*hmsgms(iah,iam,as)
      defat=hmsgms(idg,idm,ds)
      if (isig.eq.'-') defat=-defat

      k=k+1
      xtobs(k)=rafat
      ytobs(k)=defat

      enddo

 46   close(1)



c
c     Relative (RA,Dec)s "TNO_ephemeris_topocentric minus reference_star"
c     using the reference star position (furnished by the user) and
c     the topocentric ephemeris for the instants of observations
c

      do k=1,nptos

      dx=(xtobs(k)-aref)*dabs(dcos(grarad*dref))
      dy=ytobs(k)-dref

      dxep(k)=dx
      dyep(k)=dy

      enddo


c
c     Computes terms related do parallactic angle for relative color
c     refraction model
c



      call crefra (nptos,along,alati,altit,dj,aref,dref,
     ?xtobs,ytobs,acomps,dcomps,acompt,dcompt,zen)


c     do i=1,nptos
c     write (*,*) 'acomp, dcomp = ',acomp(i),dcomp(i)
c     enddo
c     stop


c
c     Adjusts observed "TNO-star"s with theoretical
c    (topocentric_ephemeris-based) "TNO-star"s, using a linear polynom
c     which models color refraction effects in "TNO-star" observed
c     position differences 
c
c     The model for RA is
c
c
c               (X - Y)_ra = dR * comp_RA + C_ra
c
c
c     where Y is the topocentric_ephemeris-based "TNO-star" (RA), X is the
c     observed "TNO-star" (RA), dR is the difference "TNO-star" of the
c     refraction constant due to the difference in color between star and
c     TNO. The parameters dR and C_ra are fitted.
c
c     Once the RA adjustment is performed, we obtain dR and apply it to
c     the Dec adjustment. The model for DE is:
c
c
c               (X - Y)_de = dR * comp_Dec + C_de
c
c
c     and now, only the C_de parameter is fitted.
c
c     
c
c     The constants (C_ra,C_de) are the correction we are looking for to
c     bring the theoretical "TNO-star" relative path (based on the
c     ephemeris and on the assumed reference star position) to the correct
c     path consistent with the high precision observed "TNO-star" path.
c
c     The zero-point color refraction in declination cannot be directly
c     determined from declination adjustments. The correct procedure is
c     to obtain dR from RA fits and then use dR in the declination
c     adjustments to derive C_de.
c
c     Note that we cannot separate in C the individual error contributions
c     from the ephemeris and from the assumed (arbitrary or not) reference
c     star position. That is, we cannot infer from these adjustments neither
c     the actual ephemeris offset of the TNO, nor the true star position.
c     But "C" indeed corrects both contributions at the same time, allowing
c     for shifting the geocentric ephemeris around the occultation date to
c     the correct "TNO-star" relative path.
c
c
c     Another underlying assumption here is that in the short run, close to
c     the occultation, when both TNO and star are already at the same FOV,
c     the shape of the small orbital arc is insensitive to errors in the
c     orbital parameters, so that adding fixed offsets for RA and Dec
c     suffice to get the correct "TNO-star" relative path.
c       
c

      do i=1,nptos
      itira(i)=0
      enddo

      ntirai=0

c

 47   continue

      n=nptos
      ngrau=1

c
c     Right ascension adjustment
c

      do i=1,nptos
      xest(i)=acompt(i)
      xesta(i)=-acomps(i)
      xp(i)=dxob(i)-dxep(i)
      enddo



      call polyna (ngrau,n,xest,xesta,xp,coef,coefe,sig,rex,itira,
     ?sclip)



      sigra=sig*3600.d3

      do i=1,3
      xcoef(i)=coef(i)
      xcoefe(i)=coefe(i)
      enddo




c
c     Declination adjustment
c

      do i=1,3
      coef(i)=0.d0
      coefe(i)=0.d0
      enddo


      ngrau=0
      n=nptos

      do i=1,nptos
      xest(i)=0.d0
      xp(i)=dyob(i)-dyep(i)-(xcoef(2)*dcompt(i)-xcoef(3)*dcomps(i))
      enddo

      sig=0.d0


      call polyno (ngrau,n,xest,xp,coef,coefe,sig,rey,itira,sclip)

c

      ntiraf=0
      do i=1,nptos
      if (itira(i).ne.0) ntiraf=ntiraf+1
      enddo

      if (ntirai.ne.ntiraf) then
      ntirai=ntiraf
      go to 47
      endif
 
c

      sigde=sig*3600.d3

      ycoef(1)=coef(1)
      ycoefe(1)=coefe(1)

c

      coef(3)=xcoef(3)
      coefe(3)=xcoefe(3)

      coef(4)=xcoef(2)
      coefe(4)=xcoefe(2)

      coef(2)=ycoef(1)
      coefe(2)=ycoefe(1)

      coef(1)=xcoef(1)
      coefe(1)=xcoefe(1)


      xxsig=coefe(1)*3600.d0
      xsig=coefe(1)
      xsig=distan*dsin(grarad*xsig)

      yysig=coefe(2)*3600.d0
      ysig=coefe(2)
      ysig=distan*dsin(grarad*ysig)



c
c     Computes corrected geocentric relative path "TNO minus star" for
c     the time span of the occultation using the ephemeris+reference_star
c     position correction (constants C_ra and C_de in the model above)
c

      ex=coef(1)
      ey=coef(2)

      eexx=ex*3600.d0
      eeyy=ey*3600.d0





      do i=1,ngeocc




      dx=(xgeocc(i)-aref)*dabs(dcos(grarad*dref))+ex
      dy=ygeocc(i)-dref+ey



      ddx=distan*dsin(grarad*dx)
      ddy=distan*dsin(grarad*dy)



c     bde=grarad*(ygeocc(i)+ey)
c     bra=grarad*(xgeocc(i)+ex)

c     d=dexy(bra,bde,grac,gdec)
c     xx=xpad(bra,bde,grac)/d
c     yy=ypad(bra,bde,grac,gdec)/d


c     ddx=distan*dsin(xx)
c     ddy=distan*dsin(yy)


      dxocc(i)=ddx
      dyocc(i)=ddy


      enddo




c
c     Writes plot file with (O-C)s from adjustment
c

      ntot=0

      open (1,file=plot)

      do i=1,nptos

      rex(i)=rex(i)*3600.d3
      rey(i)=rey(i)*3600.d3

      ip=40

      if (itira(i).eq.0) then
      ntot=ntot+1
      ip=43
      endif

      write (1,50) zen(i),rex(i),rey(i),ip
 50   format(3(1x,f16.10),1x,i2)


      enddo

      close (1)

c
c     Writes file with corrected TNO geocentric relative trajectories
c
c     Columns:
c
c     1) Dacosd relative TNO position
c     2) Dd relative TNO South border position
c     3) Dd relative TNO position
c     4) Dd relative TNO North border position
c     5) instant of time (JD)
c
c
c
c
c

      open (1,file=ptno)


      do i=1,ngeocc   

      x=dxocc(i)
      y=dyocc(i)

      t=datae(i)


      if (dabs(x).lt.xsize .and. dabs(y).lt.ysize) then

      write (1,60) x,y-rtno,y,y+rtno,t
 60   format(5(1x,f20.8))

      endif
 

      enddo

      close (1)


c
c     Writes file with Earth circunference for projection of shadow
c

      open (1,file=pearth)

      dpi=2.0*pi

      erro=dabs(datan2(stepx,rearth))

      a=erro/2.d0

 64   a=a+erro

      if (a.gt.dpi) go to 65

      x=rearth*dcos(a)
      y=rearth*dsin(a)

      write (1,*) x,y

      go to 64

 65   close (1)


c
c     Computes C/A, P/A, instant of closest approach to geocenter, shadow
c     speed and respective errors.
c


      open (1,file=traj)


c
c     Nearest ephemeris point of geocenter close approach
c

      dismin=1.d14

      do i=1,ngeocc

      dd=dxocc(i)**2+dyocc(i)**2

      if (dd.lt.dismin) then
      dismin=dd
      im=i
      endif
    
      enddo

c
c     Auxiliary geometric computations. C/A and central instant. 
c

      li=30

      if ((im+li).gt.ngeocc) li=ngeocc
      if ((im-li).lt.1) li=1



      a=dsqrt(dxocc(im-li)**2+dyocc(im-li)**2)

      b=dsqrt(dxocc(im+li)**2+dyocc(im+li)**2)

      c=dsqrt((dxocc(im+li)-dxocc(im-li))**2+(dyocc(im+li)-dyocc(im-li)
     ?)**2)


      d=dabs(a**2-((a**2-b**2+c**2)/(2.d0*c))**2)


      xx=dabs((a**2-b**2+c**2)/(2.d0*c))


      d=dsqrt(d)

      dt=datae(im+li)-datae(im-li)

      t0=(xx/c)*dt+datae(im-li)


      cakm=d


      dj1=t0
      dj2=0.d0

c
c     Gregorian date for instant of closest approach 
c



      call iau_jd2cal (dj1,dj2,iutano,iutmes,iutdia,fd,jjj)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut =((hora-iuth)*60.d0-iutm)*60.d0



c
c     Computes speed of shadow in the sky plane (Km/JD)
c


      velx=(dxocc(im+li)-dxocc(im-li))/(dt*86400.d0)
      vely=(dyocc(im+li)-dyocc(im-li))/(dt*86400.d0)



      velsec=c/(dt*86400.d0)                 



c
c
c     Position angle P/A.
c
c     TNO relative position to the star at closest approach.
c
c     North = zero degrees.
c
c     Counted clockwise.
c
c



      deltat=(t0-datae(im-li))*86400.d0

      x=dxocc(im-li)+deltat*velx
      y=dyocc(im-li)+deltat*vely



      ang=dabs(datan2(dabs(y),dabs(x)))


      if (x.ge.0.d0) then

       if (y.lt.0.d0) then

       ango=pi/2.d0+ang

       else

       ango=pi/2.d0-ang

       endif

      else

       if (y.lt.0.d0) then

       ango=1.5d0*pi-ang

       else
    
       ango=1.5d0*pi+ang

       endif
 
      endif


      angul=ango*radgra




c
c     Errors:  C/A (km) and central instant (s)
c

      xsigra=distan*dsin(grarad*sigra/3600.d3)
      ysigde=distan*dsin(grarad*sigde/3600.d3)


      quade=dsqrt((dsin(ango)*xsig)**2+(dcos(ango)*ysig)**2)
      quads=dsqrt((dsin(ango)*xsigra)**2+(dcos(ango)*ysigde)**2)

      ecakm=quads
      ecakms=quade

      quade=dsqrt((dcos(ango)*xsig)**2+(dsin(ango)*ysig)**2)
      quads=dsqrt((dcos(ango)*xsigra)**2+(dsin(ango)*ysigde)**2)

      etc =quads/dabs(velsec)
      etcs=quade/dabs(velsec)


      casec=3600.d0*radgra*dabs(datan2(cakm,distan))
      ecasec=3600.d0*radgra*dabs(datan2(ecakm,distan))
      ecases=3600.d0*radgra*dabs(datan2(ecakms,distan))



      ex=distan*dsin(grarad*ex)
      ey=distan*dsin(grarad*ey)



      write (1,*)
      write (1,*)

      write (1,66) coef(3), coefe(3)   
 66   format('TNO  color factor dr and error           = ',2(1x,f20.8))


      write (1,67) coef(4), coefe(4)   
 67   format('Star color factor dr and error           = ',2(1x,f20.8))


      write (1,*)
      write (1,*)

      write (1,68) ex, ey , eexx*1.d3, eeyy*1.d3  
 68   format('Offsets (dx,dy) (km, mas) = ',4(1x,f16.8))





      write (1,*)
      write (1,*)

      write (1,69) xsig, ysig, xxsig*1.d3, yysig*1.d3   
 69   format('Error offsets (dx,dy) (km, mas) = ',4(1x,f16.8))




      write (1,*)
      write (1,*)

      write (1,70) cakm, ecakm, ecakms
 70   format('C/A, sigma error, standard error (km) = ',3(1x,f16.8))

      write (1,*)

      write (1,71) casec*1.d3, ecasec*1.d3, ecases*1.d3
 71   format('C/A, sigma error, standard error (mas) = ',3(1x,f16.8))

      write (1,*)
      write (1,*)

      write (1,72) iuth,iutm,sut,iutdia,iutmes,iutano
 72   format('Instant of closest geocentric approach in the sky plane (h
     ? m s day mo yr) = ',2(i2.2,1x),f6.3,1x,2(i2.2,1x),i4)

       write (1,73) etc,etcs
 73   format('Instant error and standard error (s)  = ',2(1x,f17.3))


      qetc =etc*dabs(velsec)
      qetcs=etcs*dabs(velsec)

      write (1,75) qetc,qetcs
 75   format('Instant error and standard error (km) = ',2(1x,f17.3))


      setc =3600.d3*radgra*dabs(datan2(qetc,distan))
      setcs=3600.d3*radgra*dabs(datan2(qetcs,distan))


      write (1,76) setc,setcs
 76   format('Instant error and standard error (mas)= ',2(1x,f17.3))


      sclms=setc/etc
      sclmk=setc/qetc
      
      write (1,*)
      write (1,*)

      write (1,77) sclms,sclmk
 77   format('Scales:  mas/s  and  mas/km           = ',2(1x,f17.3))

      write (1,78) 1.d0/sclms,1.d0/sclmk
 78   format('Scales:  s/mas  and  km/mas           = ',2(1x,f17.3))


      if (velx.lt.0.d0) velsec=-velsec

      write (1,*)
      write (1,*)
      write (1,79) velsec
 79   format('Shadow velocity in sky plane (Km/s) = ',f7.2)

      write (1,*)
      write (1,*)
      write (1,*)


c


      write (1,90) angul
 90   format('Position angle P/A (degrees) = ',f16.8)



c
c     Reference star position and average (RA,Dec) from observations
c

      write (1,*)
      write (1,*)
      write (1,*) 'Reference star position'
      write (1,*)

c

      ra=aref/15.d0
      de=dref
      iah=ra
      am=(ra-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0
      if (de.lt.0.d0) then
      isig='-'  
      de=-de
      else
      isig='+' 
      endif 
      idg=de
      dm=(de-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0
 
      write(1,*) iah,iam,sa,'  ',isig,idg,idm,ds



      write (1,*)
      write (1,*)
      write (1,*) 'Averaged star position from observations'
      write (1,*)

c

      ra=racsta/15.d0
      de=decsta
      iah=ra
      am=(ra-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0
      if (de.lt.0.d0) then
      isig='-'  
      de=-de
      else
      isig='+' 
      endif 
      idg=de
      dm=(de-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0
 
      write(1,*) iah,iam,sa,'  ',isig,idg,idm,ds


      write (1,*)
      write (1,*)


c
c     Comparison with reference occultation values
c 


      write (1,*)'Comparison: "Obs - reference" occultation values'
      write (1,*)


      dango=angul-rpa


      if (dabs(dango).lt.45.d0) then
      difca=casec*1.d3-rca
      else
      difca=casec*1.d3+rca
      endif

      dt=(t0-rdate)*86400.d0


      write (1,93) ntot
 93   format('Number of fitted points XY   = ',i10)

c     write (1,993) ntotx
c993  format('Number of fitted points X    = ',i10)

c     write (1,994) ntoty
c994  format('Number of fitted points Y    = ',i10)


      write (1,94) sigra,sigde
 94   format('Sigma of fit (ra,de) (mas)   = ',2(1x,f16.8))


      write (1,95) difca
 95   format('C/A (mas)                    = ',f16.8)

      dts=dt*sclms

      write (1,97) dt,dts
 97   format('Central instant (s and mas)  = ',2(1x,f16.8))


      dr=dsqrt(difca**2+dts**2)

      write (1,98) dr
 98   format('C/A+t0 radial distance (mas) = ',f16.8)


      write (1,99) dango
 99   format('Position angle P/A (degrees) = ',f16.8)


      write (1,100) velsec-rvel
 100  format('Velocity (km/s)              = ',f16.8)


      write (1,*)
      write (1,*)


      write (1,105) rca  
 105  format('C/A reference (mas)          = ',f16.8)


      write (1,106) irhh,irmm,rss,irdd,irmo,iryy
 106  format('Reference t0 (h m s d mo yr) = ',2(i2.2,1x),f6.3,1x,
     ?2(i2.2,1x),i4)


      write (1,107) rpa
 107  format('Reference angle P/A (deg.)    = ',f16.8)


      write (1,108) rvel
 108  format('Reference velocity (km/s)     = ',f16.8)



      write (1,*)
      write (1,*)


      close (1)


      write (*,*)
      write (*,*) 'Execution terminated ok.'
      write (*,*)




      end



c
c     Subroutine polyna
c
c
c     Mono-variated polynomial fit of nth-degree. 
c
c
c     Last modified:  M. Assafin    01/Nov/2010
c
c
c

      subroutine polyna (ngrau,nptos,xest,xesta,xp,coef,coefe,sig,rexx,
     ?itira,sclip)



      implicit real *8 (a-h,o-z)
      dimension coef(21),alpha(21,21),array(21,21),beta(21),term(21),
     ?coefe(21)
      dimension xest(100000),xesta(100000),xp(100000),itira(100000),
     ?rexx(100000)

      common /a7/array
      common/a14/ierro


C
C     Initializing data
C

      pi=3.141592653589793d0
      grarad=pi/180.d0
      radgra=180.d0/pi
      det  =1.d0
      ierro=0


      ipmax=100000
c


c     do i=1,ipmax
c     itira(i)=0
c     enddo


c
c     Compute no. of terms of polynomial model
c



      nterms=3

C
C     Checks No. of points to fit .vs. coefficients to fit
C

 1    continue

      sig=0.d0


      if (nptos.lt.nterms) then
      ierro=1
      return
      endif

c
      do 3 i=1,nterms
      beta(i) =0.d0
      coef(i) =0.d0
      coefe(i) =0.d0
      term(i) =0.d0
      do 2 j=i,nterms
      array(i,j)=0.d0
      alpha(i,j)=0.d0
 2    continue
 3    continue



C
C     Condition equation
C

 5    continue


      
      do 10 i=1,nptos

      if (itira(i).ne.0) go to 10

      x=xest(i)
      xx=xesta(i)
      xg=xp(i)

      term(1)=1.d0
      term(2)=x
      term(3)=xx


C
C     Fills in terms for AtB
C



      do l=1,nterms
      beta(l)=beta(l)+xg*term(l)
      do k=l,nterms
      alpha(l,k)=alpha(l,k)+term(k)*term(l)
      enddo
      enddo   

   10 continue

c
c     Fills symetric inferior triangular matrix for AtA
c


      do l=1,nterms
      do k=l,nterms
      alpha(k,l)=alpha(l,k)
      enddo
      enddo

c
c     Fills ARRAY=AtA for inversion (elements normalized by diagonal)
c


C
      do l=1,nterms
      do k=1,nterms
      array(l,k)=alpha(l,k)/dsqrt(alpha(l,l)*alpha(k,k))
      enddo
      enddo

C
C     Inverting AtA
C



      call matinv (nterms,det)
      if (ierro.eq.1) return

C
C     Computing polynomial coefficients
C

      do l=1,nterms
      do k=1,nterms
      array(l,k)=array(l,k)/dsqrt(alpha(k,k)*alpha(l,l))
      enddo
      enddo

      do l=1,nterms
      do k=1,nterms
      coef(l)=coef(l)+array(l,k)*beta(k)
      enddo
      enddo


c
c     Computing residuals over used points
c


      resm=0.d0
      res2=0.d0
      nn=0

      resmax=-1.d14

      jmax=-1

      do 160 i=1,nptos

      if (itira(i).ne.0) go to 160

      nn=nn+1

      x=xest(i)
      xx=xesta(i)
      xg=xp(i)


      pol=coef(1)+coef(2)*x+coef(3)*xx

      res=pol-xg

      rexx(i)=res

      aux=dabs(res)

      if (aux.gt.resmax) then
      resmax=aux
      jmax=i
      endif 

c     write (*,*) 'res = ',i,res

      resm=resm+res
      res2=res2+res**2


  160 continue

C
C     Computes average and standard deviation of (O-C)s
C

      call desvio (nn,resm,res2)

      sig=res2

c
c     One-by-one point sigma clip outlier cuttof
c



      jm=jmax
      resm=resmax
 
      if (resm.gt.sclip) then
      itira(jm)=1
      go to 1
      endif



c
c     Computes coefficient errors
c


      do i=1,nterms
      coefe(i)=sig*dsqrt(array(i,i))
      enddo

c

      return
      end









c
c     Subroutine polyno
c
c
c     Mono-variated polynomial fit of nth-degree. 
c
c
c     Last modified:  M. Assafin    01/Nov/2010
c
c
c

      subroutine polyno (ngrau,nptos,xest,xp,coef,coefe,sig,rexx,
     ?itira,sclip)



      implicit real *8 (a-h,o-z)
      dimension coef(21),alpha(21,21),array(21,21),beta(21),term(21),
     ?coefe(21)
      dimension xest(100000),xp(100000),itira(100000),rexx(100000)

      common /a7/array
      common/a14/ierro


C
C     Initializing data
C

      pi=3.141592653589793d0
      grarad=pi/180.d0
      radgra=180.d0/pi
      det  =1.d0
      ierro=0


      ipmax=100000
c


c     do i=1,ipmax
c     itira(i)=0
c     enddo


c
c     Compute no. of terms of polynomial model
c



      nterms=ngrau+1

C
C     Checks No. of points to fit .vs. coefficients to fit
C

 1    continue

      sig=0.d0


      if (nptos.lt.nterms) then
      ierro=1
      return
      endif

c
      do 3 i=1,nterms
      beta(i) =0.d0
      coef(i) =0.d0
      coefe(i) =0.d0
      term(i) =0.d0
      do 2 j=i,nterms
      array(i,j)=0.d0
      alpha(i,j)=0.d0
 2    continue
 3    continue



C
C     Condition equation
C

 5    continue


      
      do 10 i=1,nptos

      if (itira(i).ne.0) go to 10

      x=xest(i)
      xg=xp(i)

C
C     Computs coefficient terms for  AtB
C

c     do k=1,nterms
c     term(k)=x**(k-1)
c     enddo


      term(1)=1.d0

      do l=1,nterms
      beta(l)=beta(l)+xg*term(l)
      do k=l,nterms
      alpha(l,k)=alpha(l,k)+term(k)*term(l)
      enddo
      enddo   

   10 continue

c
c     Fills symetric inferior triangular matrix for AtA
c


      do l=1,nterms
      do k=l,nterms
      alpha(k,l)=alpha(l,k)
      enddo
      enddo

c
c     Fills ARRAY=AtA for inversion (elements normalized by diagonal)
c


C
      do l=1,nterms
      do k=1,nterms
      array(l,k)=alpha(l,k)/dsqrt(alpha(l,l)*alpha(k,k))
      enddo
      enddo

C
C     Inverting AtA
C



      call matinv (nterms,det)
      if (ierro.eq.1) return

C
C     Computing polynomial coefficients
C

      do l=1,nterms
      do k=1,nterms
      array(l,k)=array(l,k)/dsqrt(alpha(k,k)*alpha(l,l))
      enddo
      enddo

      do l=1,nterms
      do k=1,nterms
      coef(l)=coef(l)+array(l,k)*beta(k)
      enddo
      enddo


c
c     Computing residuals over used points
c


      resm=0.d0
      res2=0.d0
      nn=0

      resmax=-1.d14

      jmax=-1

      do 160 i=1,nptos

      if (itira(i).ne.0) go to 160

      nn=nn+1

      x=xest(i)
      xg=xp(i)


      pol=0.d0
      do k=1,nterms
c     pol=pol+coef(k)*x**(k-1)
      pol=pol+coef(k)          
      enddo

      res=pol-xg

      rexx(i)=res

      aux=dabs(res)

      if (aux.gt.resmax) then
      resmax=aux
      jmax=i
      endif 

c     write (*,*) 'res = ',i,res

      resm=resm+res
      res2=res2+res**2


  160 continue

C
C     Computes average and standard deviation of (O-C)s
C

      call desvio (nn,resm,res2)

      sig=res2

c
c     One-by-one point sigma clip outlier cuttof (not in use here)
c



      jm=jmax
      resm=resmax
 
      if (resm.gt.sclip) then
      itira(jm)=1
      go to 1
      endif



c
c     Computes coefficient errors
c


      do i=1,nterms
      coefe(i)=sig*dsqrt(array(i,i))
      enddo

c

      return
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
      SUBROUTINE MATINV (NORDER, DET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARRAY (21,21), IK(21), JK(21)
      COMMON /A7/ARRAY
      COMMON /A14/IERRO
C
      DET = 1.D0
      DO 100 K=1, NORDER
C
C        FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
C
      AMAX= 0.D0
   21 DO 30 I=K, NORDER
      DO 30 J=K, NORDER
      IF (DABS(AMAX) - DABS(ARRAY(I,J))) 24, 24, 30
   24 AMAX = ARRAY(I,J)
      IK(K) = I
      JK(K) = J
   30 CONTINUE
C
C        INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)
C
      IF (AMAX) 41, 32, 41
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
      DO 80 I=1, NORDER
      DO 80 J=1, NORDER
      IF (I-K) 74, 80, 74
   74 IF (J-K) 75, 80, 75
   75 ARRAY(I,J) = ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
   80 CONTINUE
      DO 90 J=1, NORDER
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
      DO 130 L=1, NORDER
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
c    
c     Subrotine desvio
c
c     Computes mean and standard deviation about the mean
c
c
c     input variables:
c
c        xvam  = sum of point values
c        xvas  = sum of square of values
c
c     output variables:
c
c        xvam  = mean
c        xvas  = standard deviation about that mean
c
c


      subroutine desvio (nest,xvam,xvas)

      implicit double precision (a-h,o-z)

c

      dnove=99.999d0

c

      exmed=xvam/nest

c

      if (nest.eq.1) then

      xvas=dnove

      return

      endif

c

      if (nest.eq.2) then

      xvas=dsqrt(dabs(2.d0*xvas-xvam**2))
      xvam=exmed

      return

      endif

c

      raiz=xvas-2.d0*exmed*xvam+nest*exmed**2

      if (raiz.lt.0.d0) then

      xvas=0.d0

      else

      xvas=dsqrt(raiz/(nest-1.d0))

      endif
c

      xvam=exmed


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




      DOUBLE PRECISION FUNCTION iau_GMST82 ( DJ1, DJ2 )
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
      iau_GMST82 = iau_ANP ( DS2R * ( (A+(B+(C+D*T)*T)*T) + F ) )

*  Finished.

*+----------------------------------------------------------------------
*
  
      return
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

      return
      END




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
c
c
c   Subroutine crefra
c
c
c     Computes terms in (RA,Dec) related to the relative color
c     refraction model.
c
c
c     input variables:
c
c
c        along  = geographic longitude (degrees)
c        alati  = geographic latitude  (degrees)
c        altit  = altitude (not used in this version)
c        dj     = instant of observation (Julian Date)
c        ra     = RA  of observation (degrees)
c        de     = Dec of observation (degrees)
c        n      = number of observations        
c
c
c     auxiliary variables
c
c        c      = paralatic angle
c        ha     = hour angle
c        z      = zenith distance
c
c
c
c
c     output variables:
c
c        acomp = component RA  of color refraction model 
c        dcomp = component Dec of color refraction model
c        zen   = zenith distance (degrees)
c
c
c      Last modification:  M. Assafin    06/Apr/2011
c
c


      subroutine crefra (n,dlong,dlati,dltit,daju,ras,des,
     ?xtobs,ytobs,acomps,dcomps,acompt,dcompt,zen)

      implicit double precision (a-h,o-z)

      double precision iau_GMST82

      dimension daju(100000),acomps(100000),dcomps(100000),
     ?acompt(100000),dcompt(100000)

      dimension xtobs(100000),ytobs(100000),zen(100000)


      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0



c
c     Initial data
c
c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      dpi=2.d0*pi
c


      dj1=2400000.5D0

c
      altit=dltit
      along=dlong*grarad
      alati=dlati*grarad

      ra=ras*grarad
      de=des*grarad

      if (ra.gt.dpi) ra=ra-dpi
      if (ra.lt.0.d0) ra=ra+dpi

c
c     A and B from Laplace Formula for refraction (in degrees)
c


      a=57.085d0/3600.d0
      b=0.0666d0/3600.d0


c
c     Computs star and TNO components to all observations
c


      do 100 i=1,n



c
c     Computes local sideral time
c


      dj2=daju(i)-dj1


      tsl=iau_GMST82(dj1,dj2)-along

      if (tsl.gt.dpi) tsl=tsl-dpi





c
c     Computes Hour Angle of observation for star (sha) and TNO (tha)
c

      tra=xtobs(i)*grarad

      sha=hoang(tsl,ra)

      tha=hoang(tsl,tra)



c
c     Computes zenith distance and elevation of star (sz,sh) and TNO (tz,th)
c


      tde=ytobs(i)*grarad

      call azen (alati,sha,de,sa,sz,sh)

      call azen (alati,tha,tde,ta,tz,th)

      zen(i)=radgra*(sz+tz)/2.d0


c
c     Parallactic angle for star (sc) and TNO (tc)
c


      call platic (alati,sha,de,sz,sc)

      call platic (alati,tha,tde,tz,tc)




c
c     Components (RA,Dec) of relative color refraction model of
c     the observation for star and TNO
c



c
c     This is K.R. Lang Astrophysics Formulae book p. 23 
c
c
c
c     sh=sh*radgra
c     rs=(0.28d0*press/(temp+273.d0))*0.0167d0/dtan(grarad*(sh+7.31d0/
c    ?(sh+4.4d0)))
c
c     th=th*radgra
c     rt=(0.28d0*press/(temp+273.d0))*0.0167d0/dtan(grarad*(th+7.31d0/
c    ?(th+4.4d0)))
c
c


c
c     This is J. Kovalavsky Modern Astrometry book p. 38
c
c     A, B in degrees.
c
c


      rs=a*dtan(sz)-b*dtan(sz)**3
 
      rt=a*dtan(tz)-b*dtan(tz)**3

      rs=-rs
      rt=-rt
 
c


      acomps(i)=rs*dsin(sc)/dcos(de)
      dcomps(i)=rs*dcos(sc)

      acompt(i)=rt*dsin(tc)/dcos(tde)
      dcompt(i)=rt*dcos(tc)



c
c     This is HA-only formula with r=a*dtan(z) approximation
c
c
c     acomps(i)=rs*(dsin(sha)/dcos(de)**2)/(dcos(sha)+dtan(alati)*
c    ?dtan(de))
c     dcomps(i)=rs*(dtan(alati)-dtan(de)*dcos(sha))/(dcos(sha)+
c    ?dtan(alati)*dtan(de))
c
c     acompt(i)=rt*(dsin(tha)/dcos(tde)**2)/(dcos(tha)+dtan(alati)*
c    ?dtan(tde))
c     dcompt(i)=rt*(dtan(alati)-dtan(tde)*dcos(tha))/(dcos(tha)+
c    ?dtan(alati)*dtan(tde))
c






 100  continue

c


      return
      end







c
c
c      Function hoang
c
c
c      Computes hour angle.
c
c   
c      All variables are in radians.
c
c
c      Last modification:  M. Assafin    22/Mar/2011
c
c

      double precision function hoang (tsl,ra)


      implicit real*8 (a-h,o-z)

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      dpi=2.d0*pi


c
c     Computes Hour Angle of observation
c

      ha=tsl-ra

      if (ha.gt.dpi)  ha=ha-dpi

      if (ha.lt.-dpi) ha=ha+dpi

      hoang=ha


      return
      end



c
c
c      Subroutine azen
c
c
c      Computes azimuth, zenith distance and elevation.
c
c   
c      All variables are in radians.
c
c
c      Last modification:  M. Assafin    22/Mar/2011
c
c

      subroutine azen (alati,ha,de,a,z,h)


      implicit real*8 (a-h,o-z)

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      dpi=2.d0*pi



c
c     Computes zenith distance
c

      zz=dsin(de)*dsin(alati)+dcos(de)*dcos(ha)*dcos(alati)

      z=dabs(dasin(dabs(zz)))

      if (zz.lt.d0) z=-z

      h=z
      z=pi/2.d0-z


c
c     computes azimuth (A=0 degrees at North), clockwise (positive to East)
c



      cx=(dsin(h)*dcos(alati)-dcos(de)*dcos(ha))/(dcos(h)*dsin(alati))
      cy=-dcos(de)*dsin(ha)/dcos(h)

      cc=dabs(datan2(dabs(cy),dabs(cx)))


      if (cx.ge.0.d0) then

      if (cy.ge.0.d0) then
      c=pi/2.d0-cc
      else
      c=pi/2.d0+cc
      endif
 
      else

      if (cy.ge.0.d0) then
      c=1.5d0*pi+cc
      else
      c=1.5d0*pi-cc
      endif

      endif


      a=c


      return
      end





c
c
c      Subroutine platic
c
c
c      Computes parallactic angle.
c
c   
c      All variables are in radians.
c
c
c      Last modification:  M. Assafin    22/Mar/2011
c
c

      subroutine platic (alati,ha,de,z,c)


      implicit real*8 (a-h,o-z)

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      dpi=2.d0*pi


c
c     Parallactic angle 
c


      cx=(dsin(alati)-dsin(de)*dcos(z))/(dcos(de)*dsin(z))
      cy=dcos(alati)*dsin(ha)/dsin(z)

      cc=dabs(datan2(dabs(cy),dabs(cx)))


      if (cx.gt.0.d0) then

      if (cy.gt.0.d0) then
      c=cc
      else
      c=dpi-cc
      endif
 
      else

      if (cy.gt.0.d0) then
      c=pi-cc
      else
      c=pi+cc
      endif

      endif


      return
      end

