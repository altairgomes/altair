c
c
c     PRAIA_relative_distance
c
c
c     Computes the relative distance between two targets taking into consideration
c     color refraction. Both targets must be in the same FOV (CCD frame).
c     Preferably, oservations should have been done in the same night, and should
c     be distributed over a wide range of zenith distances.
c
c     The (RA,Dec) reference targets' position are arbitrary, but we recommend
c     the use of averaged observed coordinates, preferably consistent
c     with the observations used in the adjustment.
c
c     The "target_2 minus target_1" relative position differences are computed
c     from PRAIA output files, specifically from PRAIA offset file types.
c     Only data collected from common frames where both targets are imaged toghether
c     are considered in the calculations.  
c 
c     Relative distances are computed prior and after color refraction correction
c     for comparison. The reference relative distance is given and the (dacosd,dd)
c     offset to correct it from color refraction is also furnished with respective
c     errors. The final color-refraction-corrected relative distance is displayed
c     in the sense "target_2 minus target_1". Its error is the offset error.
c
c
c     Useful to tie occultation chords involving two nearby target stars.
c
c
c     
c     Last modification: M. Assafin  13/June/2011
c
c

      implicit double precision (a-h,o-z)

      parameter(stdin=5,stdout=6,idim=100001)


      dimension tra(idim),tde(idim),sra(idim),sde(idim),
     ?xest(idim),yest(idim),xp(idim),yp(idim),coef(21),
     ?coefe(21),dj(idim),daju(idim),xtobs(idim),ytobs(idim),
     ?datae(idim),datao(idim),datat(idim)
     
      dimension dxob(idim),dyob(idim),dxep(idim),dyep(idim)

      dimension acomp(idim),dcomp(idim),zen(idim),rex(idim),
     ?rey(idim),itira(idim),mtira(idim)


      character*50 cstar1,cstar2,tfits(idim),sfits(idim),traj,ler1,ler2,
     ?ler3,ler4,ker1,ker2,ker3,ker4,plot1,plot2

      character*20 iobalv,ichfil

      character*1 isig

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0


c
c     Auxiliary values
c    

      idimc=50

      idimp=21

      dj1=2400000.5D0

      pi=0.3141592653589793d1
      grarad=pi/180.d0
      radgra=180.d0/pi


      tbox=1.d0/86400.d0
      tbox=0.1d0*tbox
 

c

 1    format(a50)

      read (*,1) cstar1
      read (*,1) cstar2


c
c     Stores (RA,Dec) reference position target 1
c

      ler1=''
      ker1=''

      read (*,1) ler1 
      ker1=ler1

      is=+1
      do i=1,idimc
      if (ler1(i:i).eq.'-') then
      is=-1            
      ler1(i:i)=' '
      endif
      enddo

      read (ler1,*) iah,iam,as,idg,idm,ds

      aref1=15.d0*hmsgms(iah,iam,as)
      dref1=is*hmsgms(idg,idm,ds)


c
c     Stores (RA,Dec) reference position target 2
c

      ler2=''
      ker2=''

      read (*,1) ler2 
      ker2=ler2

      is=+1
      do i=1,idimc
      if (ler2(i:i).eq.'-') then
      is=-1            
      ler2(i:i)=' '
      endif
      enddo

      read (ler2,*) iah,iam,as,idg,idm,ds

      aref2=15.d0*hmsgms(iah,iam,as)
      dref2=is*hmsgms(idg,idm,ds)


c
c     Stores (RA,Dec) observed position target 1
c

      ler3=''
      ker3=''

      read (*,1) ler3 
      ker3=ler3

      is=+1
      do i=1,idimc
      if (ler3(i:i).eq.'-') then
      is=-1            
      ler3(i:i)=' '
      endif
      enddo

      read (ler3,*) iah,iam,as,idg,idm,ds

      aref3=15.d0*hmsgms(iah,iam,as)
      dref3=is*hmsgms(idg,idm,ds)


c
c     Stores (RA,Dec) observed position target 2
c

      ler4=''
      ker4=''

      read (*,1) ler4 
      ker4=ler4


      is=+1
      do i=1,idimc
      if (ler4(i:i).eq.'-') then
      is=-1            
      ler4(i:i)=' '
      endif
      enddo

      read (ler4,*) iah,iam,as,idg,idm,ds

      aref4=15.d0*hmsgms(iah,iam,as)
      dref4=is*hmsgms(idg,idm,ds)

c

      read (*,*) along
      read (*,*) alati
      read (*,*) altit

c
      read (*,*) sclip

      read (*,1) traj
      read (*,1) plot1
      read (*,1) plot2
 


c
c     Opens output log file
c

      open (1,file=traj)



      write (1,*)
      write (1,*)
      write (1,*)'Results'
      write (1,*)
      write (1,*)
      write (1,*)
      write (1,*)

      write (1,*) 'Target 1 reference position = ',ker1
      write (1,*) 'Target 2 reference position = ',ker2

      write (1,*)

      write (1,*) 'Target 1 observed position = ',ker3
      write (1,*) 'Target 2 observed position = ',ker4


      write (1,*)
      write (1,*)
      write (1,*)



c



 10   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,2(1x,i5))

c
c     Stores observed data target 1
c

      open(2,file=cstar1)


      do i=1,idim

 15   tfits(i)=''

      read (2,10,end=20) dx,dy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,tra(i),tde(i),iuth,iutm,sut,
     ?iutano,iutmes,iutdia,datat(i),iexps,ichfil,tfits(i),iobalv,nx,ny

      enddo

 20   close (2)

      nt1=i-1

      

c
c     Stores star observed data
c

      open(2,file=cstar2)

      
      do i=1,idim

 25   sfits(i)=''


      read (2,10,end=30) dx,dy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,sra(i),sde(i),iuth,iutm,sut,
     ?iutano,iutmes,iutdia,daju(i),iexps,ichfil,sfits(i),iobalv,nx,ny

      enddo

 30   close (2)

      nt2=i-1


      
c
c     Common fields. Compute relative distances "Target_2 minus Target_1".
c



      k=0

      do 41 i=1,nt1
      do 40 j=1,nt2

      dt=dabs(datat(i)-daju(j))

      if (dt.gt.tbox) go to 40

      dx=15.d0*(sra(j)-tra(i))*dabs(dcos(grarad*sde(j)))
      dy=sde(j)-tde(i)


      k=k+1

      dxob(k)=dx
      dyob(k)=dy

      dj(k)=daju(j)


 40   continue

 41   continue

      nptos=k


c
c     Computes average relative distance without color refraction correction.
c     Outliers are eliminated by sigma-clip procedure.


      do i=1,idim
      mtira(i)=0
      enddo


 42   continue


      resx=0.d0
      resy=0.d0
      resx2=0.d0
      resy2=0.d0

      k=0

      do 43 i=1,nptos

      if (mtira(i).ne.0) go to 43

      k=k+1

      resx=resx+dxob(i)
      resy=resy+dyob(i)

      resx2=resx2+dxob(i)**2
      resy2=resy2+dyob(i)**2


 43   continue


      call desvio (k,resx,resx2)
      call desvio (k,resy,resy2)

      xmed=resx
      ymed=resy

      xs2=resx2
      ys2=resy2


c
c     (O-C)s
c


      resx=0.d0
      resy=0.d0

      xxmax=-1.d14
      yymax=-1.d14


      k=0

      do 440 i=1,nptos

      if (mtira(i).ne.0) go to 440

      k=k+1

      resx=xmed-dxob(i)
      resy=ymed-dyob(i)

      if (dabs(resx).gt.xxmax) then
      jmx=i
      xxmax=dabs(resx)
      endif

      if (dabs(resy).gt.yymax) then
      jmy=i
      yymax=dabs(resy)
      endif


 440  continue


      sig=dmax1(resx2,resy2)
      sigc=sclip*sig

      jm=jmx
      resma=xxmax

      if (resy2.gt.resx2) then
      jm=jmy
      resma=yymax
      endif
 
      if (resma.gt.sigc) then
      mtira(jm)=1
      go to 42
      endif



c
c     Nominal relative distances furnished by the user
c


      write (1,*)
      write(1,*)'Nominal relative distance from reference positions:'
      write (1,*)'(Target_2 minus Target_1)'
      write (1,*)
      write (1,*)

      aref=(aref2-aref1)*dabs(dcos(grarad*dref1))
      dref=dref2-dref1

      write (1,46) aref
      write (1,47) dref
      write (1,*)
      write (1,*)
      write (1,*)


      write (1,*)
      write(1,*)'Nominal relative distance from observed positions:'
      write (1,*)'(Target_2 minus Target_1)'
      write (1,*)
      write (1,*)

      arefo=(aref4-aref3)*dabs(dcos(grarad*dref3))
      drefo=dref4-dref3

      write (1,46) arefo
      write (1,47) drefo
      write (1,*)
      write (1,*)
      write (1,*)


c
c     Relative distance without color-refraction
c



      write (1,*)
      write(1,*)'Relative distance without color refraction correction:'
      write (1,*)'(Target_2 minus Target_1)'
      write (1,*)
      write (1,*)


      write (1,44) nptos
      write (1,45) k
      write (1,46) xmed
      write (1,47) ymed
      write (1,48) xs2*3600.d3
      write (1,49) ys2*3600.d3


 44   format(1x,'Total frame No.  = ',i6)
 45   format(1x,'Final frame No.  = ',i6)
 46   format(1x,'Dacosd (degrees) = ',f16.11)
 47   format(1x,'Dd     (degrees) = ',f16.11)
 48   format(1x,'sigma_a (mas)    = ',f16.11)
 49   format(1x,'sigma_d (mas)    = ',f16.11)


      write (1,*)
      write (1,*)
      write (1,*)



c
c     Theorethical (reference) relative (RA,Dec)s "target_2 minus target_1"
c     using the nominal star positions furnished by the user.
c

      do k=1,nptos

      dxep(k)=aref
      dyep(k)=dref

      enddo


c
c     Computes terms related do parallactic angle for relative color
c     refraction model
c


      arefm=(aref1+aref2)/2.d0
      drefm=(dref1+dref2)/2.d0


      call crefra (idim,nptos,along,alati,altit,dj,arefm,drefm,acomp,
     ?dcomp,zen)


c     do i=1,nptos
c     write (*,*) 'acomp, dcomp = ',acomp(i),dcomp(i)
c     enddo
c     stop


c
c     Adjusts theoretical relative distance with observed relative positions
c     using a linear polynom which models color refraction effects in observed
c     position differences 
c
c     The model for RA is
c
c
c               (X - Y)_ra = dR * comp_RA + C_ra
c
c
c     and for Dec is
c
c
c               (X - Y)_de = dR * comp_Dec + C_de
c
c
c
c     where (X - Y)_ra,de are the observed relative distance affected by color
c     refraction, dR is the contribution of color refraction to the observed
c     position difference and C_ra,de are the corrections to the adopted nominal
c     relative distances furnished by the user. After the adjustment, dR, C_ra and
c     C_de are computed. Add C_ra,de to the nominal distances to obtain the true
c     (RA,Dec) relative distances corrected from color refraction.
c
c
c     Note that we cannot separate in C_ra,de the individual error contributions
c     from each target. That is, we cannot infer from these adjustments neither
c     the actual position of target 1 nor the position of target 2. But C_ra,de
c     indeed corrects both contributions at the same time, allowing for deriving
c     the correct relative distance between targets 1 and 2.
c
c
c       
c

      n=nptos


c
c     (RA,Dec) adjustment
c

      do i=1,nptos
      xest(i)=acomp(i)
      yest(i)=dcomp(i)
      xp(i)=dxob(i)-dxep(i)
      yp(i)=dyob(i)-dyep(i)
      enddo



      call polyno (idim,nptos,xest,yest,xp,yp,coef,coefe,sig,rex,rey,
     ?itira,sclip)


      k=nptos
      do i=1,nptos
      if (itira(i).ne.0) k=k-1
      enddo


      arefc=aref+coef(1)
      drefc=dref+coef(2)

      ex=coef(1)*3600.d3
      ey=coef(2)*3600.d3

      xsig=coefe(1)*3600.d3
      ysig=coefe(2)*3600.d3



      write (1,*)
      write (1,*)'Relative distance with color refraction correction:'
      write (1,*)'(Target_2 minus Target_1)'
      write (1,*)
      write (1,*)


      write (1,44) nptos
      write (1,45) k
      write (1,46) arefc
      write (1,47) drefc

      write (1,51) ex
      write (1,52) ey

 51   format(1x,'corr. RA  (mas)  = ',f16.11)
 52   format(1x,'corr. Dec (mas)  = ',f16.11)

      write (1,48) xsig
      write (1,49) ysig


      write (1,*)
      write (1,*)
      write (1,*)

c
c     Comparison of solutions against nominal reference distance
c


      dx=(arefo-aref)*3600.d3
      dy=(drefo-dref)*3600.d3
  

      write (1,*)
      write (1,*)'Observed minus reference distance:'
      write (1,*)

      write (1,55) dx
      write (1,56) dy


 55   format(1x,'Dacosd (mas)     = ',f16.11)
 56   format(1x,'Dd     (mas)     = ',f16.11)


c

      dx=(xmed-aref)*3600.d3
      dy=(ymed-dref)*3600.d3
  

      write (1,*)
      write (1,*)
      write (1,*)
      write (1,*)'No_color_refraction_correction minus reference distanc
     ?e:'
      write (1,*)

      write (1,55) dx
      write (1,56) dy

c


      dx=(arefc-aref)*3600.d3
      dy=(drefc-dref)*3600.d3
  

      write (1,*)
      write (1,*)
      write (1,*)
      write (1,*)'Color_refraction_correction minus reference distance:'
      write (1,*)

      write (1,55) dx
      write (1,56) dy


      write (1,*)
      write (1,*)
      write (1,*)



      close (1)


c
c     Writes plot file with (O-C)s from adjustment (color correction)
c

      open (7,file=plot1)

      do i=1,nptos

      rex(i)=rex(i)*3600.d3
      rey(i)=rey(i)*3600.d3

      if (itira(i).eq.0) then

      write (7,60) zen(i),rex(i),rey(i)
 60   format(3(1x,f16.10))

      endif

      enddo

      close (7)




c
c     Writes plot file with (O-C)s from averaged s-clip without
c     color correction
c

      open (7,file=plot2)

      do i=1,nptos

      rx=(xmed-dxob(i))*3600.d3
      ry=(ymed-dyob(i))*3600.d3

      if (mtira(i).eq.0) then

      write (7,60) zen(i),rx,ry

      endif

      enddo

      close (7)


c



      write (*,*)
      write (*,*) 'Execution terminated ok.'
      write (*,*)




      end


c
c     Subroutine polyno
c
c
c     (RA,Dec) fit to observations based on color refraction model.
c
c
c     Fit function is of the kind:
c
c
c     Da = R*x + A
c     Dd = R*y + B
c
c
c     where the parameters being adjusted are A, B and R.
c
c
c
c     Last modified:  M. Assafin    22/Mar/2011
c
c
c

      subroutine polyno (idim,nptos,xest,yest,xp,yp,coef,coefe,sig,rex,
     ?rey,itira,sclip)


      implicit double precision (a-h,o-z)
      dimension coef(21),alpha(21,21),array(21,21),beta(21),coefe(21)
      dimension xest(idim),yest(idim),xp(idim),yp(idim),
     ?itira(idim),rex(idim),rey(idim)

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


      ipmax=idim
c


      do i=1,ipmax
      itira(i)=0
      enddo


c
c     No. of terms (parameters) of adjustment
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
      y=yest(i)
      xg=xp(i)
      yg=yp(i)

C
C     Fills in terms for AtB
C

      beta(1)=beta(1)+xg
      beta(2)=beta(2)+yg
      beta(3)=beta(3)+xg*x+yg*y

      alpha(1,1)=alpha(1,1)+1.d0
      alpha(1,2)=alpha(1,2)+0.d0
      alpha(1,3)=alpha(1,3)+x
      alpha(2,2)=alpha(2,2)+1.d0
      alpha(2,3)=alpha(2,3)+y
      alpha(3,3)=alpha(3,3)+x**2+y**2


   10 continue



c
c     Fills in symetric inferior triangular matrix for AtA
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
C     Computing coefficients
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


      resmx=0.d0
      res2x=0.d0

      resmy=0.d0
      res2y=0.d0

      nn=0

      resmax=-1.d14
      resmay=-1.d14

      jmax=-1
      jmay=-1

      do 160 i=1,nptos

      if (itira(i).ne.0) go to 160

      nn=nn+1

      x=xest(i)
      y=yest(i)
      xg=xp(i)
      yg=yp(i)


      polx=coef(1)+x*coef(3)       
      poly=coef(2)+y*coef(3)       

      resx=polx-xg
      resy=poly-yg

      rex(i)=resx
      rey(i)=resy

      auxx=dabs(resx)
      auyy=dabs(resy)

      if (auxx.gt.resmax) then
      resmax=auxx
      jmax=i
      endif 

      if (auyy.gt.resmay) then
      resmay=auyy
      jmay=i
      endif 

c     write (*,*) 'resx = ',i,resx
c     write (*,*) 'resy = ',i,resy

      resmx=resmx+resx
      res2x=res2x+resx**2

      resmy=resmy+resy
      res2y=res2y+resy**2


  160 continue

C
C     Computes average and standard deviation of (O-C)s
C

      call desvio (nn,resmx,res2x)
      call desvio (nn,resmy,res2y)

      sigx=res2x
      sigy=res2y

c
c     One-by-one point sigma clip outlier cuttof (not in use here)
c

      sig=dmax1(sigx,sigy)
      sig3=sclip*sig
      jm=jmax
      resm=resmax
      if (sigy.gt.sigx) then
      jm=jmay
      resm=resmay
      endif
 
      if (resm.gt.sig3) then
      itira(jm)=1
      go to 1
      endif



c
c     Computes coefficient errors
c

      sig=dsqrt(sigx**2+sigy**2)

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


      subroutine crefra (idim,n,dlong,dlati,dltit,daju,ras,des,acomp,
     ?dcomp,zen)

      implicit double precision (a-h,o-z)

      double precision iau_GMST82

      dimension daju(idim),acomp(idim),dcomp(idim)

      dimension zen(idim)


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

c
c     A and B from Laplace Formula for refraction (in degrees)
c


      a=57.085d0/3600.d0
      b=0.0666d0/3600.d0


c
c     Computs star and TNO components to all observations
c


      do 100 i=1,n

      ra=ras*grarad
      de=des*grarad

      if (ra.gt.dpi) ra=ra-dpi

      if (ra.lt.0.d0) ra=ra+dpi


c
c     Computes local sideral time
c

      iii=daju(i)-0.5d0
      dj1=iii+0.5d0
      dj2=daju(i)-dj1


      tsl=iau_GMST82(dj1,dj2)-along

      if (tsl.gt.dpi) tsl=tsl-dpi



c
c     Computes Hour Angle (ha) of observation
c


      ha=hoang(tsl,ra)


c
c     Computes zenith distance
c


      call azen (alati,ha,de,a,z,h)

      zen(i)=z*radgra

c
c     Parallactic angle 
c


      call platic (alati,ha,de,z,c)




c
c     Components (RA,Dec) of relative color refraction model of
c     the observation
c

c
c     This is K.R. Lang Astrophysics Formulae book p. 23 
c
c     h=h*radgra
c     r=0.0167d0/dtan(grarad*(h+7.31d0/(h+4.4d0)))


c
c     This is J. Kovalevsky Modern Astrometry book p. 38
c

      r=a*dtan(z)-b*dtan(z)**3



      acomp(i)=r*dsin(c)/dcos(de)
      dcomp(i)=r*dcos(c)


c     if (ha.lt.0.d0) acomp(i)=-acomp(i)
c     if (alati.ge.0.d0) then
c     if (de.lt.alati) dcomp(i)=-dcomp(i) 
c     else
c     if (de.gt.alati) dcomp(i)=-dcomp(i)
c     endif


c
c     This is HA-only formula with r=a*dtan(z) approximation
c
c     acomp(i)=r*(dsin(ha)/dcos(de)**2)/(dcos(ha)+dtan(alati)*
c    ?dtan(de))
c     dcomp(i)=r*(dtan(alati)-dtan(de)*dcos(ha))/(dcos(ha)+
c    ?dtan(alati)*dtan(de))




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

