c
c
c     Program PRAIA_position_format
c
c
c     Purpose
c
c
c
c     Given a list of "xy" PRAIA files, offsets-type PRAIA files or a MPC
c     format files of positions, program formats the files to another format
c     of type "xy" PRAIA file, offsets-type PRAIA file or a MPC format file
c     of positions. 
c
c
c     In this version, a more MPC-compliant format is used.
c
c
c      Last update: Marcelo Assafin - 04/Ago/2011
c   
c
c
c


      IMPLICIT REAL *8 (A-H,O-Z)

      character*50 list,filein,iext,infits
      character*101 fileou
      character*1 ipoint,isig
      character*20 ichfil,mchobj
      character*1 obtype,band,type
      character*3 iau,miau
      character*14 obname
      character*300 ler
      

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0


c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c


c
c
c     Reading input batch file
c

      write (*,*)
      write (*,*) 
      write (*,*) 'PRAIA position format'
      write (*,*) 
      write (*,*) 


c
c     open (77,file='PRAIA_position_format_02.dat')
c

      ipoint='.'
      dnove=99.9d0
      dj1=2400000.5D0


c

 1    continue

      read (*,*,end=50) keyin
      read (*,*,end=50) keyout

      if (keyin.lt.1 .or. keyin.gt.3) stop
      if (keyout.lt.1 .or. keyout.gt.3) stop

      read (*,10) list
      read (*,10) iext
      read (*,*)
      read (*,2) obname
      read (*,3) obtype
      read (*,3) band
      read (*,4) iau

 2    format(a14)
 3    format(a1)
 4    format(a3)

      read (*,*)


 10   format (a50)

c     close (77)




c
c     Number of files
c

      nfile=0

      open (9,file=list,err=17)

      i=0

 13   read (9,*,end=15)
      i=i+1
      go to 13

 15   nfile=i

c

 17   continue

      if (nfile.eq.0) then
      close (9)

      write (*,*)
      write (*,20) list
 20   format(1x,'List empty or does not exit: ',a50)
      go to 1
      endif

c

      rewind (9)


c
c     Formats of input files
c


c
c     xy-type PRAIA format
c


 22   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?3(1x,i5),2(1x,f7.3))


c
c     offset-type PRAIA format
c


 24   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,3(1x,i5),2(1x,f7.3))



c
c     MPC position format
c


      xob=dnove
      yob=dnove
      seng=9.99d0
      altu=dnove
      fgcc=dnove
      fumag=dnove
      fumag2=dnove
      xmgu=dnove
      cudmg2=dnove
      xmgj=dnove
      xmgh=dnove
      xmgk=dnove
      res2mg=dnove
      resmg2=dnove
      ermgj=dnove
      ermgh=dnove
      ermgk=dnove
      pma=dnove
      pmd=dnove
      epma=dnove
      epmd=dnove
      ex=dnove
      ey=dnove
      erau=dnove
      edeu=dnove
      alfsiu=dnove
      delsiu=dnove
      nstaru=0
      nfinau=0
      alsiu=dnove
      desiu=dnove
      ktirau=99
      iexps=0
      infits='unknown'
      nx=0
      ny=0
      ncom=0
      aux=dnove
      auy=dnove


 26   format(a14,a1,i4.4,1x,i2.2,1x,f8.5,1x,i2.2,1x,i2.2,1x,f6.3,a1,
     ?i2.2,1x,i2.2,1x,f5.2,9x,f5.2,a1,6x,a3)



c
c     Converts formats
c


      do 40 i=1,nfile

      filein=''
      fileou=''


      read (9,10) filein

      fileou=filein

      do k=50,1,-1
      if (filein(k:k).ne.' ') go to 27
      enddo

 27   k=k+1
      fileou(k:k)=ipoint


      do m=50,1,-1
      if (iext(m:m).ne.' ') go to 28
      enddo

 28   fileou(k+1:k+m)=iext


      open (1,file=filein)
      open (2,file=fileou)



 30   continue


c
c     Reads xy-type PRAIA input file
c

      if (keyin.eq.1) then

      read (1,22,end=35) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy

      endif



c
c     Reads offset-type PRAIA input file
c

      if (keyin.eq.2) then


      read (1,24,end=35) zx,zy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy

      endif



c
c     Reads MPC file format type with observed (RA,DEC,Gregorian Date)
c

      if (keyin.eq.3) then

      read (1,26,end=35) mchobj,type,iutano,iutmes,dia,iah,iam,sa,
     ?isig,idg,idm,ds,cudmg,ichfil,miau

      write (2,26) mchobj,type,iutano,iutmes,dia,iah,iam,sa,
     ?isig,idg,idm,ds,cudmg,ichfil,miau


c
c     Computes Julian Date
c

      iutdia=dia
      fd=dia-iutdia

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
      dj=djm+djm0


c
c     (RA,DEC) in decimal fraction
c

      ra=hmsgms(iah,iam,sa)
      de=hmsgms(idg,idm,ds)
      if (isig.eq.'-') de=-de


c
c     hours, minutes and seconds in standard format 
c

      fd=fd*24.d0

      iuth=fd
      iutm=(fd-iuth)*60.d0
      sut=((fd-iuth)*60.d0-iutm)*60.d0


      endif


c
c     Writes new format in output file
c


c
c     Writes xy-type PRAIA output file
c

      if (keyout.eq.1) then

      write (2,22) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy

      endif



c
c     Writes offset-type PRAIA output file
c

      if (keyout.eq.2) then


      write (2,24) zx,zy,xob,yob,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy

      endif



c
c     Writes MPC output file format
c

      if (keyout.eq.3) then

      if (keyin.eq.3) then

      write (2,26) mchobj,type,iutano,iutmes,dia,iah,iam,sa,
     ?isig,idg,idm,ds,cudmg,ichfil,miau

      go to 30

      endif


c
c     (RA,DEC) in hexasegimal notation
c

      raco=ra
      deco=de
      iah=raco
      am=(raco-iah)*60.d0
      iam=am
      sa =(am-iam)*60.d0
      if (deco.lt.0.d0) then
      isig='-'  
      deco=-deco
      else
      isig='+' 
      endif
      idg=deco
      dm=(deco-idg)*60.d0
      idm=dm
      ds=(dm-idm)*60.d0



c
c     Gregorian date from Julian date
c

      djm=dj-dj1

      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,jjj)

      dia=iutdia+fd

c

      ler='' 


      write (ler,26) obname,obtype,iutano,iutmes,dia,iah,iam,sa,
     ?isig,idg,idm,ds,cudmg,band,iau


      if (ler(24:24).eq.' ') ler(24:24)='0'

      if (ler(39:39).eq.' ') ler(39:39)='0'

      if (ler(52:52).eq.' ') ler(52:52)='0'

      if (ler(66:66).eq.' ') ler(66:66)='0'

c

      write (2,33) ler
 33   format(a80)

      endif

c

      go to 30

c

 35   continue


      close (1)
      close (2)


 40   continue


c

      close (9)



c
c     Continue to next block
c

      go to 1


c
c     That's all, folks!
c

 50   continue

      write (*,*) 
      write (*,*) 
      write (*,*) 'Execution terminated succesfully.'
      write (*,*) 
      write (*,*) 
      write (*,*) 

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



