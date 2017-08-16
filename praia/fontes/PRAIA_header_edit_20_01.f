c
c
c
c     PRAIA_header_edit    
c
c
c
c     Purpose
c
c     Edits the header info obtained from the PRAIA_header_extraction task.
c     (RA,DEC), instant, exposure time, filter, object name, all these
c     information can be updated.
c
c
c     Comments
c
c   
c     Changes in the filter, object name and exposure time info are just cosmetic.
c     In particular, the mid-instant of the observation is not changed.
c 
c
c     A time offset can be applied, and the mid-instant of all the observations
c     are updated (Gregorian and Julian dates).
c
c
c     The year can be edited, and the mid-instant of all the observations are
c     updated (Gregorian and Julian dates).  
c
c
c     The (RA,DEC) can be precessed from  observation epoch to J2000
c     (using IAU 1976 model).
c
c
c     The (RA,DEC) can also be replaced by the (RA,DEC) extracted from the
c     ephemeris of a body (in PRAIA format).
c
c
c     The (RA,DEC) can also be altered by a cos(DEC) factor.
c
c
c
c     Last Modification: M. Assafin - 11/Mar/2013
c
c



      IMPLICIT REAL *8 (A-H,O-Z)


      character*50 lista1,lista2,kephem
      character*50 infits
      character*1  isig,iver
      character*20 ichobj,ichfil,kchobj,kchfil

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0



c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

      d2000=2.451545d6

c
c     Reads input data      
c


c     open (21,file='PRAIA_header_edit_20_01.dat')

c

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'PRAIA - Header Edit')

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '

c
c     Input data loop
c


 2    continue

c

      read (*,3,err=100,end=100) lista1

      read (*,3) kephem

      read (*,3) lista2

      read (*,*) jchobj,kchobj

      read (*,*) jchfil,kchfil

      read (*,*) jexps,kexps

      read (*,*) jutano,kutano

      read (*,*) jtime,otime

      read (*,*) jmul

      read (*,*) jdiv

      read (*,*) jpre

      read (*,*) jephem,tbox

      tbox=tbox/86400.d0

      jdate=0
      if (jutano.ne.1 .and. jtime.ne.1) jdate=1


      read (*,*,err=100,end=100) 
 3    format(a50)


c
c     Retrieves header information
c


      open (10,file=lista1)
      open (30,file=lista2)


 10   continue

c

      read (10,15,end=50) iah,iam,sa,isig,idg,idm,ds,iuth,iutm,sut,
     ?iutano,iutmes,iutdia,djm,dj,iexps,ichfil,infits,ichobj,nx,ny
 15   format(1x,i2,1x,i2,1x,f7.4,1x,a1,i2,1x,i2,1x,f6.3,2x,i2,
     ?1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,f16.8,1x,f16.8,2x,i4.4,2x,
     ?a20,2x,a50,1x,a20,2(1x,i5))

      datar=dj
      datam=djm

      muth=iuth
      mutm=iutm
      sutm=sut
      mutano=iutano
      mutmes=iutmes
      mutdia=iutdia


      ar=hmsgms(iah,iam,sa)

      dr=hmsgms(idg,idm,ds)
      if (isig.eq.'-') dr=-dr

      ar=ar*15.d0*grarad
      dr=dr*grarad


c
c     Updates object name
c


      if (jchobj.eq.1) then

      ichobj=kchobj

      endif


c
c     Updates filter
c


      if (jchfil.eq.1) then

      ichfil=kchfil

      endif



c
c     Updates exposure time stamp (mid-instant of observations remain
c     unchanged)
c


      if (jexps.eq.1) then

      iexps=kexps

      endif


c
c     Updates year stamp (mid-instant of observations will be recomputed)
c


      if (jutano.eq.1) then

      iutano=kutano

      endif


c
c     Time offset (mid-instant of observations will be recomputed)
c

      dtime=0.d0

      if (jtime.eq.1) then

      dtime=otime

      endif



c
c     Updates mid-instant of observation
c


      fd=hmsgms(iuth,iutm,sut+dtime)/24.d0


      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflago)

      djm=djm+fd
            
      dj=djm+djm0

      call iau_jd2cal (djm0,djm,iutano,iutmes,iutdia,fd,iflago)

      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut =((hora-iuth)*60.d0-iutm)*60.d0


      if (jdate.eq.1) then

      dj=datar
      djm=datam

      iuth=muth
      iutm=mutm
      sut=sutm
      iutano=mutano
      iutmes=mutmes
      iutdia=mutdia

      endif





c
c     Multiply or divide by cos(DEC)
c


      if (jmul.eq.1) then

      ar=ar*dcos(dabs(dr))

      endif


      if (jdiv.eq.1) then

      ar=ar/dcos(dabs(dr))

      endif




c
c     Precess (RA,DEC) from observation epoch to J2000.0
c


      if (jpre.eq.1) then



      ep0=(dj-d2000)/365.25d0+2000.d0
      ep1=2000.d0


      call prece(ep0,ep1,ar,dr)


      endif


c
c     Replaces (RA,DEC) by ephemeris of target in the FOV
c


      if (jephem.eq.1) then

      ar=0
      dr=0

      open(20,file=kephem)

 20   continue

      read (20,25,end=30) iah,iam,sa,isig,idg,idm,ds,date
 25   format(1x,i2.2,1x,i2.2,1x,f9.6,1x,a1,i2.2,1x,i2.2,1x,f8.5,1x,
     ?f16.8)

      dt=dabs(date-dj)

      if (dt.lt.tbox) then

      ar=hmsgms(iah,iam,sa)

      dr=hmsgms(idg,idm,ds)
      if (isig.eq.'-') dr=-dr

      ar=ar*15.d0*grarad
      dr=dr*grarad

      go to 30

      endif

      go to 20

 30   close (20)

      endif




c
c     Recovers (RA,DEC) in hexasegimal notation
c

      ra=ar
      de=dr

      ra=ra*radgra/15.d0

      if (ra.lt.0.d0) ra=ra+24.d0
      if (ra.gt.24.d0) ra=ra-24.d0

      de=de*radgra


      iah=ra
      iam=(ra-iah)*60.d0
      sa=((ra-iah)*60.d0-iam)*60.d0

      if (de.lt.0.d0) then
      isig='-'
      else
      isig='+'
      endif

      de=dabs(de)

      idg=de
      idm=(de-idg)*60.d0
      ds=((de-idg)*60.d0-idm)*60.d0


  


c
c     Updates header extraction information 
c


      write (30,40) iah,iam,sa,isig,idg,idm,ds,iuth,iutm,sut,
     ?iutano,iutmes,iutdia,djm,dj,iexps,ichfil,infits,ichobj,nx,ny

 40   format(1x,i2.2,1x,i2.2,1x,f7.4,1x,a1,i2.2,1x,i2.2,1x,f6.3,2x,i2.2,
     ?1x,i2.2,1x,f5.2,1x,i4.4,1x,i2.2,1x,i2.2,f16.8,1x,f16.8,2x,i4.4,2x,
     ?a20,2x,a50,1x,a20,2(1x,i5.5))




      go to 10

 50   close (10)
      close (30)

c

      go to 2


 100  continue


      write (*,*) ' '
      write (*,*) 'Execution terminated.'
      write (*,*) ' '

      end




      SUBROUTINE prece (EP0, EP1, RA, DC)
*+
*     - - - - - - -
*      P R E C E
*     - - - - - - -
*
*  Precession - either FK4 (Bessel-Newcomb, pre IAU 1976) or
*  FK5 (Fricke, post IAU 1976) as required.
*
*  Given:
*     SYSTEM     char   precession to be applied: 'FK4' or 'FK5'
*     EP0,EP1    dp     starting and ending epoch
*     RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP0
*
*  Returned:
*     RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP1
*
*  Called:    slDA2P, slPRBN, slPREC, slDS2C,
*             slDMXV, slDC2S
*
*  Notes:
*
*     1)  Lowercase characters in SYSTEM are acceptable.
*
*     2)  The epochs are Besselian if SYSTEM='FK4' and Julian if 'FK5'.
*         For example, to precess coordinates in the old system from
*         equinox 1900.0 to 1950.0 the call would be:
*             CALL slPRCE ('FK4', 1900D0, 1950D0, RA, DC)
*
*     3)  This routine will NOT correctly convert between the old and
*         the new systems - for example conversion from B1950 to J2000.
*         For these purposes see slFK45, slFK54, slF45Z and
*         slF54Z.
*
*     4)  If an invalid SYSTEM is supplied, values of -99D0,-99D0 will
*         be returned for both RA and DC.
*
*  P.T.Wallace   Starlink   20 April 1990
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT real*8 (a-h,o-z)
C     IMPLICIT NONE

C     CHARACTER SYSTEM*(*)
      DOUBLE PRECISION EP0,EP1,RA,DC

      DOUBLE PRECISION PM(3,3),V1(3),V2(3)
C     CHARACTER SYSUC*3

C     DOUBLE PRECISION slDA2P


C     SYSTEM='FK5'

*  Convert to uppercase and validate SYSTEM
C     SYSUC=SYSTEM

*     Generate appropriate precession matrix
         CALL PREC(EP0,EP1,PM)

*     Convert RA,Dec to x,y,z
         CALL slDS2C(RA,DC,V1)

*     Precess
         CALL slDMXV(PM,V1,V2)

*     Back to RA,Dec
         CALL slDC2S(V2,RA,DC)
C        RA=slDA2P(RA)



      END




      SUBROUTINE PREC (EP0, EP1, RMATP)
*+
*     - - - - -
*      P R E C
*     - - - - -
*
*  Form the matrix of precession between two epochs (IAU 1976, FK5)
*  (double precision)
*
*  Given:
*     EP0    dp         beginning epoch
*     EP1    dp         ending epoch
*
*  Returned:
*     RMATP  dp(3,3)    precession matrix
*
*  Notes:
*
*     1)  The epochs are TDB (loosely ET) Julian epochs.
*
*     2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
*
*     3)  Though the matrix method itself is rigorous, the precession
*         angles are expressed through canonical polynomials which are
*         valid only for a limited time span.  There are also known
*         errors in the IAU precession rate.  The absolute accuracy
*         of the present formulation is better than 0.1 arcsec from
*         1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
*         and remains below 3 arcsec for the whole of the period
*         500BC to 3000AD.  The errors exceed 10 arcsec outside the
*         range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
*         5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.
*         The SLALIB routine slPREL implements a more elaborate
*         model which is suitable for problems spanning several
*         thousand years.
*
*  References:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*      equations (6) & (7), p283.
*     Kaplan,G.H., 1981. USNO circular no. 163, pA2.
*
*  Called:  slDEUL
*
*  P.T.Wallace   Starlink   23 August 1996
*
*  Copyright (C) 1996 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT real*8 (a-h,o-z)

      DOUBLE PRECISION EP0,EP1,RMATP(3,3)

*  Arc seconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER (AS2R=0.484813681109535994D-5)

      DOUBLE PRECISION T0,T,TAS2R,W,ZETA,Z,THETA



*  Interval between basic epoch J2000.0 and beginning epoch (JC)
      T0 = (EP0-2000D0)/100D0

*  Interval over which precession required (JC)
      T = (EP1-EP0)/100D0

*  Euler angles
      TAS2R = T*AS2R
      W = 2306.2181D0+(1.39656D0-0.000139D0*T0)*T0

      ZETA = (W+((0.30188D0-0.000344D0*T0)+0.017998D0*T)*T)*TAS2R
      Z = (W+((1.09468D0+0.000066D0*T0)+0.018203D0*T)*T)*TAS2R
      THETA = ((2004.3109D0+(-0.85330D0-0.000217D0*T0)*T0)
     :        +((-0.42665D0-0.000217D0*T0)-0.041833D0*T)*T)*TAS2R

*  Rotation matrix
      CALL slDEUL('ZYZ',-ZETA,THETA,-Z,RMATP)

      END


      SUBROUTINE slDEUL (ORDER, PHI, THETA, PSI, RMAT)
*+
*     - - - - - - -
*      D E U L
*     - - - - - - -
*
*  Form a rotation matrix from the Euler angles - three successive
*  rotations about specified Cartesian axes (double precision)
*
*  Given:
*    ORDER   c*(*)   specifies about which axes the rotations occur
*    PHI     d       1st rotation (radians)
*    THETA   d       2nd rotation (   "   )
*    PSI     d       3rd rotation (   "   )
*
*  Returned:
*    RMAT    d(3,3)  rotation matrix
*
*  A rotation is positive when the reference frame rotates
*  anticlockwise as seen looking towards the origin from the
*  positive region of the specified axis.
*
*  The characters of ORDER define which axes the three successive
*  rotations are about.  A typical value is 'ZXZ', indicating that
*  RMAT is to become the direction cosine matrix corresponding to
*  rotations of the reference frame through PHI radians about the
*  old Z-axis, followed by THETA radians about the resulting X-axis,
*  then PSI radians about the resulting Z-axis.
*
*  The axis names can be any of the following, in any order or
*  combination:  X, Y, Z, uppercase or lowercase, 1, 2, 3.  Normal
*  axis labelling/numbering conventions apply;  the xyz (=123)
*  triad is right-handed.  Thus, the 'ZXZ' example given above
*  could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ').  ORDER
*  is terminated by length or by the first unrecognized character.
*
*  Fewer than three rotations are acceptable, in which case the later
*  angle arguments are ignored.  If all rotations are zero, the
*  identity matrix is produced.
*
*  P.T.Wallace   Starlink   23 May 1997
*
*  Copyright (C) 1997 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT real*8 (a-h,o-z)
C     IMPLICIT NONE

      CHARACTER*(*) ORDER
      DOUBLE PRECISION PHI,THETA,PSI,RMAT(3,3)

      INTEGER J,I,L,N,K
      DOUBLE PRECISION RESULT(3,3),ROTN(3,3),ANGLE,S,C,W,WM(3,3)
      CHARACTER AXIS



*  Initialize result matrix
      DO J=1,3
         DO I=1,3
            IF (I.NE.J) THEN
               RESULT(I,J) = 0D0
            ELSE
               RESULT(I,J) = 1D0
            END IF
         END DO
      END DO

*  Establish length of axis string
      L = LEN(ORDER)

*  Look at each character of axis string until finished
      DO N=1,3
         IF (N.LE.L) THEN

*        Initialize rotation matrix for the current rotation
            DO J=1,3
               DO I=1,3
                  IF (I.NE.J) THEN
                     ROTN(I,J) = 0D0
                  ELSE
                     ROTN(I,J) = 1D0
                  END IF
               END DO
            END DO

*        Pick up the appropriate Euler angle and take sine & cosine
            IF (N.EQ.1) THEN
               ANGLE = PHI
            ELSE IF (N.EQ.2) THEN
               ANGLE = THETA
            ELSE
               ANGLE = PSI
            END IF
            S = DSIN(ANGLE)
            C = DCOS(ANGLE)

*        Identify the axis
            AXIS = ORDER(N:N)
            IF (AXIS.EQ.'X'.OR.
     :          AXIS.EQ.'x'.OR.
     :          AXIS.EQ.'1') THEN

*           Matrix for x-rotation
               ROTN(2,2) = C
               ROTN(2,3) = S
               ROTN(3,2) = -S
               ROTN(3,3) = C

            ELSE IF (AXIS.EQ.'Y'.OR.
     :               AXIS.EQ.'y'.OR.
     :               AXIS.EQ.'2') THEN

*           Matrix for y-rotation
               ROTN(1,1) = C
               ROTN(1,3) = -S
               ROTN(3,1) = S
               ROTN(3,3) = C

            ELSE IF (AXIS.EQ.'Z'.OR.
     :               AXIS.EQ.'z'.OR.
     :               AXIS.EQ.'3') THEN

*           Matrix for z-rotation
               ROTN(1,1) = C
               ROTN(1,2) = S
               ROTN(2,1) = -S
               ROTN(2,2) = C

            ELSE

*           Unrecognized character - fake end of string
               L = 0

            END IF

*        Apply the current rotation (matrix ROTN x matrix RESULT)
            DO I=1,3
               DO J=1,3
                  W = 0D0
                  DO K=1,3
                     W = W+ROTN(I,K)*RESULT(K,J)
                  END DO
                  WM(I,J) = W
               END DO
            END DO
            DO J=1,3
               DO I=1,3
                  RESULT(I,J) = WM(I,J)
               END DO
            END DO

         END IF

      END DO

*  Copy the result
      DO J=1,3
         DO I=1,3
            RMAT(I,J) = RESULT(I,J)
         END DO
      END DO

      END





      SUBROUTINE slDS2C (A, B, V)
*+
*     - - - - - -
*      D S 2 C
*     - - - - - -
*
*  Spherical coordinates to direction cosines (double precision)
*
*  Given:
*     A,B       dp      spherical coordinates in radians
*                        (RA,Dec), (Long,Lat) etc
*
*  Returned:
*     V         dp(3)   x,y,z unit vector
*
*  The spherical coordinates are longitude (+ve anticlockwise
*  looking from the +ve latitude pole) and latitude.  The
*  Cartesian coordinates are right handed, with the x axis
*  at zero longitude and latitude, and the z axis at the
*  +ve latitude pole.
*
*  P.T.Wallace   Starlink   October 1984
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT real*8 (a-h,o-z)
C     IMPLICIT NONE

      DOUBLE PRECISION A,B,V(3)

      DOUBLE PRECISION COSB



      COSB=DCOS(B)

      V(1)=DCOS(A)*COSB
      V(2)=DSIN(A)*COSB
      V(3)=DSIN(B)

      END





      SUBROUTINE slDMXV (DM, VA, VB)
*+
*     - - - - -
*      D M X V
*     - - - - -
*
*  Performs the 3-D forward unitary transformation:
*
*     vector VB = matrix DM * vector VA
*
*  (double precision)
*
*  Given:
*     DM       dp(3,3)    matrix
*     VA       dp(3)      vector
*
*  Returned:
*     VB       dp(3)      result vector
*
*  P.T.Wallace   Starlink   March 1986
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT real*8 (a-h,o-z)
C     IMPLICIT NONE

      DOUBLE PRECISION DM(3,3),VA(3),VB(3)

      INTEGER I,J
      DOUBLE PRECISION W,VW(3)


*  Matrix DM * vector VA -> vector VW
      DO J=1,3
         W=0D0
         DO I=1,3
            W=W+DM(J,I)*VA(I)
         END DO
         VW(J)=W
      END DO

*  Vector VW -> vector VB
      DO J=1,3
         VB(J)=VW(J)
      END DO

      END


      SUBROUTINE slDC2S (V, A, B)
*+
*     - - - - - -
*      D C 2 S
*     - - - - - -
*
*  Direction cosines to spherical coordinates (double precision)
*
*  Given:
*     V     d(3)   x,y,z vector
*
*  Returned:
*     A,B   d      spherical coordinates in radians
*
*  The spherical coordinates are longitude (+ve anticlockwise
*  looking from the +ve latitude pole) and latitude.  The
*  Cartesian coordinates are right handed, with the x axis
*  at zero longitude and latitude, and the z axis at the
*  +ve latitude pole.
*
*  If V is null, zero A and B are returned.
*  At either pole, zero A is returned.
*
*  P.T.Wallace   Starlink   July 1989
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT real*8 (a-h,o-z)
C     IMPLICIT NONE

      DOUBLE PRECISION V(3),A,B

      DOUBLE PRECISION X,Y,Z,R


      X = V(1)
      Y = V(2)
      Z = V(3)
      R = DSQRT(X*X+Y*Y)

      IF (R.EQ.0D0) THEN
         A = 0D0
      ELSE
         A = DATAN2(Y,X)
      END IF

      IF (Z.EQ.0D0) THEN
         B = 0D0
      ELSE
         B = DATAN2(Z,R)
      END IF

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


