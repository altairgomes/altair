c
c     Programa alvos_JPL_1
c
c     PROPOSITO
c
c
c     Dadas as efemerides de satelites, planetas, etc, pelo JPL, formata-se
c     as coordenadas para arquivo alvo a ser usado na reducao de imagens
c     CCD com o programa automatico "astrometria_*".
C
C     Ultima modificacao: M. Assafin - 22/Mar/2006
C
C



      IMPLICIT REAL *8 (A-H,O-Z)

      parameter (stdin=5,stdout=6)

      character*50 lista1,lista2
      character*20 ichobj


c
c     Reads input data       
c

c

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (23x,'PRAIA: Formatted Ephemeris to Targets')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '

c

 2    format(a50)

      
      write (*,*) ' '
      write (*,*) ' '

c
c     Reads input data      
c

 3    continue

      read (5,4,end=10) ichobj
 4    format(a20)
      read (5,*) key
      read (5,2) lista1
      read (5,2) lista2
      read (5,*) 

c
      write (*,5)
 5    format (15x,'JPL ephemeris data file                    -> ',$)
      write(*,2) lista1
c
      write (*,6)
 6    format (15x,'PRAIA output target file                   -> ',$)
      write(*,2) lista2

c

      write (*,*)
      write (*,*)


      if (key.eq.1) then

      call naif (ichobj,lista1,lista2)

      else

      call hori (ichobj,lista1,lista2)

      endif

c

      go to 3

 10   continue

      write (*,*) ' Execution terminated ok.'
      write (*,*)
      write (*,*)

      end


c
c
c     Subroutine hori
c
c
c     Extracts ephemeris generated by Horizons JPL service
c
c
c     Last update:  M. Assafin 08/nov/2009
c



      subroutine hori (ichobj,lista1,lista2)

      IMPLICIT REAL *8 (A-H,O-Z)

      character*50 lista1,lista2,lista
      character*20 ichobj,ichfil
      character*1  isig
      character*3 metab(12),kmes
      character*5 icon,iconi,iconf
      character*6 ia6
      character*82 iform


      data metab/'Jan','Feb','Mar','Apr','May','Jun','Jul',
     ?'Aug','Sep','Oct','Nov','Dec'/
      data iconi/'$$SOE'/ 
      data iconf/'$$EOE'/ 

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0
c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c


c
c     Reads ephemeris from Horizons JPL service
c

      open (3,file=lista1)
      open (7,file=lista2)

c


 10   read (3,20,end=50) icon
 20   format(a5)

      if (icon.ne.iconi) go to 10

c
c     Reads in blocks                          
c

      kkey=0

      iform='(1x,i4,1x,a3,1x,i2,1x,i2,1x,i2,1x,f6.3,5x,i2,1x,i2,1x,f7.4,
     ?1x,a1,i2,1x,i2,1x,f6.3)'


 25   read (3,20,end=50) icon

      if (icon.eq.iconf) go to 10

      backspace 3

c
c     Checks for blank character in first column of file 
c


      if (kkey.eq.0) then
      read (3,27) ia6
 27   format(a6)
      backspace 3
      if (ia6(5:5).eq.'-') iform(2:2)='0'
      kkey=1
      endif



      read (3,iform) iutano,kmes,iutdia,iuth,iutm,sut,iah,iam,sa,
     ?isig,idg,idm,ds

      do i=1,12
      if (metab(i).eq.kmes) go to 35
      enddo
 35   iutmes=i


c
c     Computes Julian Date      
c

      fd=hmsgms(iuth,iutm,sut)/24.d0

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0

c
c     Writes in output file        
c


      write(7,40) iah,iam,sa,isig,idg,idm,ds,dj,ichobj
 40   format(1x,i2.2,1x,i2.2,1x,f9.6,1x,a1,i2.2,1x,i2.2,1x,f8.5,1x,
     ?f16.8,1x,a20)

      go to 25

c

 50   continue




      close (3)
      close (7)

      return
      end



c
c
c     Subroutine naif
c
c
c     Extracts ephemeris generated by NAIF JPL libraries
c
c
c     Last update:  M. Assafin 08/nov/2009
c


      subroutine naif (ichobj,lista1,lista2)

      IMPLICIT REAL *8 (A-H,O-Z)

      character*50 lista1,lista2
      character*20 ichobj
      character*1  isig
      character*1 hdr
      character*3 metab(12),kmes



      data metab/'JAN','FEB','MAR','APR','MAY','JUN','JUL',
     ?'AUG','SEP','OCT','NOV','DEC'/


      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0
c
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c


c
c     Reads Ephemeris and writes in PRAIA format
c

      open (3,file=lista1)
      open (7,file=lista2)

c
      read(3,*)hdr
      read(3,*)hdr
c      read(3,*)hdr

 10   read (3,15,end=50) iutano,kmes,iutdia,iuth,iutm,sut,iah,iam,sa,
     ?isig,idg,idm,ds

 15   format(i4,1x,a3,1x,i2,2(1x,i2),1x,f6.3,44x,2(i2,1x),f7.4,2x,a1,
     ?2(i2,1x),f6.3)

      do i=1,12
      if (metab(i).eq.kmes) go to 35
      enddo
 35   iutmes=i


c
c     Computes Julian Date
c

      fd=hmsgms(iuth,iutm,sut)/24.d0

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0

c
c     Writes in output file
c


      write(7,40) iah,iam,sa,isig,idg,idm,ds,dj,ichobj

 40   format(1x,i2.2,1x,i2.2,1x,f9.6,1x,a1,i2.2,1x,i2.2,1x,f8.5,1x,
     ?f16.8,1x,a20)


      go to 10

c

 50   continue


      close (3)
      close (7)

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

