c
c
c     Program  PRAIA_JPL_ephem_batch
c
c
c
c     Given Julian dates, generates batch files for JPL ephemeris 
c     extraction.
c
c     Due to JPL web service limitations, files must be created with
c     up to 200 entries each (more than 200 observations will produce
c     more than 1 file accordingly).
c
c     In this version, input files may be PRAIA header extraction files,
c     a list of xy files, an offset target file In PRAIA format (v14 or
c     earlier) or MPC report format file.
c
c     Use of IAU codes (JPL format) for the observer location is recommended.
c     In this version, site coordinates (longitude, latitude, altitude) can
c     be directly furnished. This is usefull when the observation site does
c     not have a IAU code or when observation takes place in remote locations
c     or is made by portable instruments.
c
c
c     In this version, emails for the JPL server are automatically sent,
c     using ssmtp under many Linux distributions (Suse, Fedora, Ubuntu, etc).
c     
c
c     Last update: M. Assafin  11/Nov/2011     
c
c
c

      implicit real*8 (a-h,o-z)

      parameter(stdin=5,stdout=6)

      dimension daju(200,50000)

      character*50 lista,infits,iraiz,email
      character*150 ixy,nome
      character*20 ichfil,ichobj,mchobj

      character*120 linha

      character*300 chama


      character*1 iund,iver
      character*5 jmaux

      character*50 iausit,iaubod,iaunam

      gms(a,b,c)=a+b/60.d0+c/3600.d0

      data iund/'_'/

c

      idimax=200
      idifls=50000
      ichau=50

      ihd5=0

c
c     Initiate vectors
c

      do i=1,idifls
      do j=1,idimax
      daju(j,i)=0.d0
      enddo
      enddo



c
      write (*,*)
      write (*,*)
      write (*,*)
      write (*,*) '       PRAIA JPL batch generator'
      write (*,*)
      write (*,*)

      write (*,*) ' Input data: '


c
c     Reads input data file
c

c     open (5,file='PRAIA_JPL_ephem_batch_20_05.dat')

      read (5,*)  kfile

      if (kfile.lt.1 .and. kfile.gt.5) stop




      read (5,01) lista
      read (5,01) iaunam
      read (5,01) iausit

      sinal=1.d0
      read  (5,2) linha
      backspace 5
      read (5,*)  dlond,dlonm,dlons

      do i=1,ichau
      if (linha(i:i).eq.'-') sinal=-1.d0
      enddo

      dlon=sinal*gms(dabs(dlond),dabs(dlonm),dabs(dlons))

      sinal=1.d0
      read  (5,2) linha
      backspace 5
      read (5,*)  dlatd,dlatm,dlats

      do i=1,ichau
      if (linha(i:i).eq.'-') sinal=-1.d0
      enddo

      dlat=sinal*gms(dabs(dlatd),dabs(dlatm),dabs(dlats))

      read (5,*)  high

      read (5,01) iaubod
      read (5,01) email

 1    format(a50)

      rewind (5)

      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
      read  (5,2) linha
      write (6,3) linha
 2    format(a120)
 3    format('  ',a120)

      write (6,*)
      write (6,*)

c     close (5)

c

      i=1
      j=0
      
      open (2,file=lista)

c
c     Retrieves Julian Dates
c



 5    continue


c
c     Case 1, header file
c


      if (kfile.eq.1) then

      read (2,10,end=40) iah,iam,sa,isig,idg,idm,ds,iuth,iutm,sut,
     ?iutano,iutmes,iutdia,djm,dj,iexps,ichfil,infits,ichobj

 10   format(1x,i2,1x,i2,1x,f7.4,1x,a1,i2,1x,i2,1x,f6.3,2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,f16.8,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20)

      go to 39

      endif

c
c     Case 2, xy file list
c


      if (kfile.eq.2) then

      read (2,20,end=40) ixy
 20   format (a150)

      open (8,file=ixy)

      read (8,30) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy

 30   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?3(1x,i5),2(1x,f7.3))

      close (8)

      go to 39

      endif


c
c     Case 3, offset target file
c

      if (kfile.eq.3) then

      read (2,31,end=40) zx,zy,xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy

 31   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,3(1x,i5),2(1x,f7.3))


      go to 39

      endif


c
c     Case 4, MPC file format with observed (RA,DEC,Gregorian Date)
c

      if (kfile.eq.4) then

      read (2,32,end=40) mchobj,iutano,iutmes,dia,iah,iam,as,
     ?isig,idg,idm,ds,cudmg,ichfil

 32   format(a14,1x,i4,1x,i2,1x,f8.5,1x,i2,1x,i2,1x,f6.3,a1,i2,1x,i2,1x,
     ?f5.2,9x,f4.1,a8)


c
c     Computes Julian Date
c

      iutdia=dia

      fd=dia-iutdia

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,iflag)

      djm=djm+fd
            
      dj=djm+djm0


      go to 39

      endif



c
c     Case 5, PRAIA big Table input format
c

      if (kfile.eq.5) then

      ihd5=ihd5+1

      if (ihd5.lt.2) then

 33   read (2,34) iver
 34   format(a1)

      if (iver.ne.'-') go to 33

      read(2,34) iver
      read(2,34) iver

      endif

      iver=''
      read (2,35,err=40) iver,dj
 35   format(a1,143x,f18.10)

      if (iver.eq.'-') go to 40


      go to 39

      endif


c

 39   continue


      j=j+1

      if (j.gt.idimax) then
      i=i+1
      j=1
      endif
      

      daju(j,i)=dj

      go to 5

c

 40   close (2)


c
c     Mounts JPL batch files
c

      iraiz=''

      do k=ichau,1,-1
      if (iaunam(k:k).ne.' ') go to 45
      enddo

 45   iraiz(1:k)=iaunam(1:k)
      iraiz(k+1:k+1)=iund

c

      write (*,*) 'Produced JPL batches -> ',i,' file(s):'
      write (*,*)
 


      do 50 n=1,i

      nome=''
      nome(1:k+1)=iraiz(1:k+1)

      write (jmaux,'(i5.5)') n

      nome(K+2:k+6)=jmaux(1:5)


      open (3,file=nome)


      write (3,46) email
 46   format ('From: ',a50)

      write (3,47)
 47   format ('To: horizons@ssd.jpl.nasa.gov')

      write (3,48)
 48   format ('Subject: JOB')

      write (3,*)

      call jbat (n,j,i,daju,iausit,iaubod,dlon,dlat,high)

      write (3,*)

      close (3)


c
c     Calls batch through mailx prompt command
c


      do kkk=ichau,1,-1
      if (nome(kkk:kkk).ne.' ') go to 49
      enddo

 49   kkk1=kkk

      chama=''
      chama(1:11)='ssmtp -t < '
      chama(12:11+kkk1)=nome(1:kkk1)


      call system (chama)

c
      write (6,*) nome

      write (6,*) 

      write (6,*) 



 50   continue

      write (*,*)
      write (*,*)
      write (*,*)
      write (*,*) ' Execution terminated succesfully.'
      write (*,*)
      write (*,*)





      end


 

c
c
c    Subrotine jbat
c
c
c    Mounts JPL batch
c
c
c    Last update M. Assafin - 16/April/2010
c
c


      subroutine jbat (nn,jj,ii,daju,iausit,iaubod,dlon,dlat,high)

      IMPLICIT REAL *8 (A-H,O-Z)

      dimension daju(200,50000)

      character*80 siteco
      
      character*50 iaubod,iausit,kaubod,kausit
      character*1 iplic,ibr

      data iplic/"'"/
      data ibr  /" "/

c
      isi=80

      idimax=200
      ichau=50

c

      kaubod=''
      kausit=''

      kaubod=iaubod
      kausit=iausit


c
      
      do i=ichau,1,-1
      if (kausit(i:i).ne.ibr) go to 1 
      enddo
 1    kausit(i+1:i+1)=iplic

      do i=ichau,1,-1
      if (kaubod(i:i).ne.ibr) go to 10
      enddo
 10   kaubod(i+1:i+1)=iplic


c
c     Writes first of the 2 file headers
c


      write (3,1000)

      write (3,1001)

      write (3,1002) iplic,kaubod

      write (3,1003)

      write (3,1004)

      write (3,1005)

      write (3,1006) iplic,kausit

      write (3,1007)

      write (3,1008)

c
c     The case when site coordinates must be input by the user
c

      do i=1,ichau-5
      if (iausit(i:i+4).eq.'coord'.or.iausit(i:i+4).eq.'COORD') then

      iausit(i:i+4)='coord'

      write (siteco,*) '  SITE_COORD = ',iplic,dlon,dlat,high,iplic

      write (3,1030) siteco


      endif
      enddo


      write (3,1009)

c
c     Writes Julian Dates    
c

      j2=idimax

      if (nn.eq.ii) j2=jj

      do j=1,j2
      write (3,1010) iplic,daju(j,nn),iplic
      enddo

      write (3,1011)

      write (3,1012)

      write (3,1013)

      write (3,1014)

      write (3,1015)

      write (3,1016)

      write (3,1017)

      write (3,1018)

      write (3,1019)

      write (3,1020)

      write (3,1021)

      write (3,1022)

      write (3,1023)

      write (3,1024)

      write (3,1025)

      write (3,1026)

      write (3,1027)



 1000 format ('!$$SOF')
 1001 format ("   EMAIL_ADDR = ' '")
 1002 format ('   COMMAND    = ',a1,a50)
 1003 format ("   OBJ_DATA   = 'NO'")
 1004 format ("   MAKE_EPHEM = 'YES'")
 1005 format ("   TABLE_TYPE = 'OBS'")
 1006 format ('   CENTER     = ',a1,a50)
 1007 format ("   REF_PLANE  = 'BODY EQUATOR'")
 1008 format ("   COORD_TYPE = 'GEODETIC'")
 1009 format ('   TLIST=')
 1010 format (8x,a1,f16.8,a1)
 1011 format ("   QUANTITIES = '1,9,10'")
 1012 format ("   REF_SYSTEM = 'J2000'")
 1013 format ("   CAL_FORMAT = 'CAL'")
 1014 format ("   ANG_FORMAT = 'HMS'")
 1015 format ("   APPARENT   = 'AIRLESS'")
 1016 format ("   TIME_DIGITS = 'FRACSEC'")
 1017 format ("   TIME_ZONE   = '+00:00'")
 1018 format ("   ELEV_CUT   = '0'")
 1019 format ("   SKIP_DAYLT = 'YES'")
 1020 format ("   SOLAR_ELONG= '0,180'")
 1021 format ("   EXTRA_PREC = 'YES'")
 1022 format ("   CSV_FORMAT = 'NO'")
 1023 format ("   VEC_LABELS = 'NO'")
 1024 format ("   ELM_LABELS = 'NO'")
 1025 format ("   R_T_S_ONLY = 'NO'")
 1026 format ("   CA_TABLE_TYPE= 'STANDARD'")
 1027 format ('!$$EOF')

 1030 format(a80)


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




