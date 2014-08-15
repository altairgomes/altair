c
c
c     PRAIA_trajectory
c
c
c     Furnishes files for plot of shadow path over Earth from
c     astrometric results of TNOs and respective occultation stars.
c
c     The shadow projection is a trajectory extrapolation from observations.
c     The TNO and the field star are at the same FOV, so that their relative
c     distances are determined with great accuracy, allowing for the
c     shadow path extrapolation.
c
c     The instants of time are also extrapolated for the trajectory.
c
c     The trajectory is assumed to be a straight line. The error of the
c     independent coefficient of the straight line is given and can be
c     interpretated as the uncertainty in the geographical latitude of the
c     path. The standard error of the path latitude along the line is
c     computed and can be used to estimate the divergence of the predicted
c     path. The miminum distance (C/A) is furnished, as well as the
c     associated instant of time and errors. The inclination of the line
c     (P/A) and its error are given as well. 
c
c     Useful for last-minute predictions and LC chord fittings.
c
c     
c     Last modification: M. Assafin  11/Nov/2010
c
c

      implicit real*8 (a-h,o-z)

      parameter(stdin=5,stdout=6)

      dimension tra(100000),tde(100000),sra(100000),sde(100000),
     ?xest(100000),xp(100000),xpe(100000),coef(21),coefe(21),dj(100000),
     ?daju(100000),xgeo(100000),ygeo(100000),xtop(100000),ytop(100000),
     ?datat(100000),datag(100000)
     
      dimension xpep(100000),coefp(21),coefep(21)
      dimension xpet(100000),coeft(21),coefet(21)

      dimension array(21,21),prray(21,21),sto1(21),sto2(21)


      character*50 ctno,cstar,ptno,pearth,tfits(100000),sfits(100000),
     ?crel,traj,geo,topo

      character*20 iobalv,ichfil

      character*1 isig

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

      common /a7/array

c
c     Auxiliary values
c    

      idim=100000

      idimp=21

      dj1=2400000.5D0

      pi=0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi

      au=149597870.691d0

c

      do i=1,idimp
      do j=1,idimp
      array(j,i)=0.d0
      prray(j,i)=0.d0
      enddo
      enddo



c

 1    format(a50)

      read (*,1) ctno
      read (*,1) cstar
      read (*,1) topo 
      read (*,1) geo  
      read (*,1) crel
      read (*,1) traj
      read (*,1) ptno
      read (*,1) pearth

      read (*,*) distan

      distan=distan*au

      read (*,*) rearth
      read (*,*) rtno

      read (*,*) xsize, ysize

      xsize=xsize/2.d0
      ysize=ysize/2.d0

      read (*,*) stepx

      read (*,*) ngraup
      read (*,*) ngraut

c


 10   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))

c
c     Stores TNO data
c

      open(2,file=ctno)


      do i=1,idim

 15   tfits(i)=''

      read (2,10,end=20) dx,dy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,tra(i),tde(i),iuth,iutm,sut,
     ?iutano,iutmes,iutdia,ddj,iexps,ichfil,tfits(i),iobalv,nx,ny

      enddo

 20   close (2)

      ntno=i-1

      

c
c     Stores star data
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
c     Common fields. Compute relative distances w.r.t. star
c     X is converted to Km.
c

      open (1,file=crel)

      write (1,*) 'TNO minus star: Dacosd (dg), Dd(dg), Julian Date'

      k=0

      do 45 i=1,ntno
      do 40 j=1,nstar

      if (tfits(i).ne.sfits(j)) go to 40

      dx=15.d0*(tra(i)-sra(j))*dabs(dcos(grarad*sde(j)))
      dy=tde(i)-sde(j)

      ddx=distan*dsin(grarad*dx)
      ddy=distan*dsin(grarad*dy)

      k=k+1

      xest(k)=ddx
      xp(k)=ddy
      dj(k)=daju(j)

c
c     Writes relative (RA,Dec)s "TNO minus star" and instants
c

      write (1,35) dx,dy,dj(j)
 35   format(2(1x,f15.10),1x,f16.8)



 40   continue

 45   continue

      nptos=k

      close (1)


c
c     Retrieves topocentric TNO ephemeris.
c

      open(1,file=topo)

      
      do i=1,idim

      read (1,46,end=47) iah,iam,as,isig,idg,idm,ds,datat(i),
     ?iobalv
 46   format(1x,i2,1x,i2,1x,f9.6,1x,a1,i2,1x,i2,1x,f8.5,1x,f16.8,
     ?1x,a20)

      rafat=15.d0*hmsgms(iah,iam,as)
      defat=hmsgms(idg,idm,ds)
      if (isig.eq.'-') defat=-defat

      xtop(i)=rafat
      ytop(i)=defat

      enddo

 47   close(1)

      ntop=i-1




c
c     Retrieves geocentric TNO ephemeris.
c

      open(1,file=geo)

      
      do i=1,idim

      read (1,46,end=48) iah,iam,as,isig,idg,idm,ds,datag(i),
     ?iobalv


      rafat=15.d0*hmsgms(iah,iam,as)
      defat=hmsgms(idg,idm,ds)
      if (isig.eq.'-') defat=-defat

      xgeo(i)=rafat
      ygeo(i)=defat

      enddo

 48   close(1)

      ngeo=i-1


c
c     Offsets apparent shadow path from topocentric to geocentric 
c     origin.
c

      do 50 k=1,nptos
      do i=1,ntop
      do j=1,ngeo

      dti=86400.d0*dabs(dj(k)-datat(i))
      dtj=86400.d0*dabs(dj(k)-datag(j))

      if (dti.lt.0.5d0 .and. dtj.lt.0.5d0 ) then

     
      dx=(xtop(i)-xgeo(j))*dabs(dcos(grarad*ygeo(j)))
      dy=ytop(i)-ygeo(j)

      ddx=distan*dsin(grarad*dx)
      ddy=distan*dsin(grarad*dy)

c     write (*,*) 'dti,dtj,dx, dy, ddx, ddy = ',dti,dtj,dx,dy,ddx,ddy

c     xest(k)=xest(k)-ddx
c     xp(k)=xp(k)-ddy

      go to 50

      endif

      enddo
      enddo

 50   continue


      aux=0.004d0/3600.d0
      aux=distan*dsin(grarad*aux)

      write (*,*) 'aux = ', aux


c
c     Adjusts relative trajectory
c

      call polyno (ngraup,nptos,xest,xp,coef,coefe,sigp)


c
c     Stores path solution
c
      
      do i=1,idimp
      do j=1,idimp
      prray(i,j)=array(i,j)
      enddo
      enddo

      do i=1,ngraup+1
      coefp(i)=coef(i)
      coefep(i)=coefe(i)
      enddo



c
c     Adjusts time with relative trajectory
c

      do i=1,nptos
      xp(i)=dj(i)
      enddo


      call polyno (ngraut,nptos,xest,xp,coef,coefe,sigt)


c
c     Stores time solution
c
      

      do i=1,ngraut+1
      coeft(i)=coef(i)
      coefet(i)=coefe(i)
      enddo





c
c     Writes file with TNO relative trajectories
c
c     Columns:
c
c     1) Dacosd relative TNO position
c     2) Dd relative TNO South border position
c     3) Dd relative TNO position
c     4) Dd relative TNO North border position
c     5) extrapolated instant of time (JD)
c     6) Dd relative TNO standard error South of solution
c     7) Dd relative TNO standard error North of solution
c
c
c
c
c

      open (1,file=ptno)


      nterms=ngraup+1
      kterms=ngraut+1


      x=-xsize-stepx/2.d0

 60   x=x+stepx

      if (x.gt.xsize) goto 62

      y=0.d0
      do k=1,nterms
      y=y+coefp(k)*x**(k-1)
      enddo


      t=0.d0
      do k=1,kterms
      t=t+coeft(k)*x**(k-1)
      enddo



c
c     Computes Dd standard errors of each point from trajectory
c     solution
c



      do n=0,ngraup
      sto1(n+1)=x**n
      enddo

      do k=1,nterms
      sto2(k)=0.d0
      enddo

      do l=1,nterms
      do k=1,nterms
      sto2(l)=sto2(l)+prray(l,k)*sto1(k)
      enddo
      enddo

      pe=0.d0

      do k=1,nterms
      pe=pe+sto1(k)*sto2(k)
      enddo

      pe=sigp*dsqrt(pe)

   

c


      write (1,61) x,y-rtno,y,y+rtno,t,y-pe,y+pe
 61   format(7(1x,f16.8))


      go to 60

 62   close (1)


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
c     Computes C/A, P/A, instant of closest approach to geocenter and
c     respective errors.
c


      open (1,file=traj)


c    
c     Computes C/A and error
c

      aux=dsqrt(1.d0+coefp(2)**2)

      cakm=dabs(coefp(1)/aux)
      ecakm=coefep(1)/aux

c     ecakm=coefep(1)/aux-coefp(1)*coefep(2)/aux**3

c     aa=dabs(coefp(1))
c     bb=dabs(coefp(1)/coefp(2))

c     cakm=aa*bb/dsqrt(aa**2+bb**2)


c     aux=dabs(datan(coefp(2)))


c     ecakm=coefep(1)*dabs(dsin(aux))


      casec=3600.d0*radgra*dabs(datan2(cakm,distan))
      ecasec=3600.d0*radgra*dabs(datan2(ecakm,distan))

      write (1,*)
      write (1,*)

      write (1,70) cakm, ecakm
 70   format('C/A and error (km) = ',2(1x,f16.8))

      write (1,*)

      write (1,71) casec, ecasec
 71   format('C/A and error (") = ',2(1x,f16.8))

      write (1,*)
      write (1,*)
 

c
c     Computes instant of closest approach and error  
c


      aux=dabs(datan(coefp(2)))

      aux=dcos(aux)*dsqrt(coefp(1)**2-cakm**2)

      xc=aux

      auy=coefp(2)*coefp(1)

      if (auy.gt.0.d0) xc=-xc


      tc=coeft(1)+xc*coeft(2)



c
c     Computes standard error of time instant from time scale solution
c



      do n=0,ngraut
      sto1(n+1)=xc**n
      enddo

      do k=1,kterms
      sto2(k)=0.d0
      enddo

      do l=1,kterms
      do k=1,kterms
      sto2(l)=sto2(l)+array(l,k)*sto1(k)
      enddo
      enddo

      te=0.d0

      do k=1,kterms
      te=te+sto1(k)*sto2(k)
      enddo

      te=sigt*dsqrt(pe)

      te=te*86400.d0



c
c     Gregorian date for instant of closest approach 
c

      djm=tc-dj1

      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,jjj)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut =((hora-iuth)*60.d0-iutm)*60.d0



c

      write (1,*)
      write (1,*)

      write (1,80) tc
 80   format('Instant of closest approach (JD) = ',f16.8)

      write (1,*)

      write (1,82) iutano,iutmes,iutdia
 82   format('yr mo dd = ',i4.4,1x,i2.2,1x,i2.2)


      write (1,84) iuth,iutm,sut
 84   format('hh hh ss = ',i2.2,1x,i2.2,1x,f6.3)


      write (1,*)

      write (1,86) te
 86   format('Instant of closest approach. Error(s) = ',f16.8)

      write (1,*)
      write (1,*)
 

      write (1,*)
      write (1,*)


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


      yc=coefp(1)+xc*coefp(2)

      
      ang=dabs(datan2(dabs(yc),dabs(xc)))

      
      ango=ang

      if (xc.gt.0.d0) then

      if (yc.lt.0.d0) ango=2.d0*pi-ang

      else

      ango=pi-ang

      if (uc.lt.0.d0) ango=pi+ang

      endif

      angul=ango-pi/2.d0

      angul=2.d0*pi-angul

      angul=angul*radgra


c
c     P/A error
c

      angule=radgra*coefep(2)*dcos(ang)**2


      write (1,90) angul
 90   format('Position angle P/A (degrees) = ',f16.8)

      write (1,*)

      write (1,92) angule
 92   format('Position angle error (degrees) = ',f16.8)




      write (1,*)
      write (1,*)

      close (1)


      write (*,*)
      write (*,*) 'Execution terminated ok.'
      write (*,*)


      x=xc

      do n=0,ngraup
      sto1(n+1)=x**n
      enddo

      do k=1,nterms
      sto2(k)=0.d0
      enddo

      do l=1,nterms
      do k=1,nterms
      sto2(l)=sto2(l)+prray(l,k)*sto1(k)
      enddo
      enddo

      pe=0.d0

      do k=1,nterms
      pe=pe+sto1(k)*sto2(k)
      enddo

      pe=sigp*dsqrt(pe)

   
      aux=3600.d0*radgra*dabs(datan2(pe,distan))

      auy=3600.d0*radgra*dabs(datan2(sigp,distan))


      write (*,*) 'pe km, aux sec = ',pe, aux


      write (*,*) 'sigp sigt = ',auy,sigt*86400.d0


      aux=6371.d0

      aux=3600.d0*radgra*dabs(datan2(aux,distan))
 
      write (*,*) 'raio terra arcsec = ',aux


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

      subroutine polyno (ngrau,nptos,xest,xp,coef,coefe,sig)


      implicit real *8 (a-h,o-z)
      dimension coef(21),alpha(21,21),array(21,21),beta(21),term(21),
     ?coefe(21)
      dimension xest(100000),xp(100000),itira(100000)

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


      do i=1,ipmax
      itira(i)=0
      enddo


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

      do k=1,nterms
      term(k)=x**(k-1)
      enddo

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
      pol=pol+coef(k)*x**(k-1)
      enddo

      res=pol-xg

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

      aux=3.d0*sig
      if (resmax.gt.aux) then
      itira(jmax)=1
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
      IMPLICIT REAL *8 (A-H,O-Z)
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

      implicit real*8 (a-h,o-z)

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










