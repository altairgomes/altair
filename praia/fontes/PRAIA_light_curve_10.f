c
c     Program PRAIA_light_curve_10
c
c     Purpose
c
c
c     Light curve composition, including phase angle LC.
c     This version holds for LNA, HP, ESO fits files.
c
c     User furnishes for each individual contributing LC: time unit of
c     recorded instants of time, instant of reference, phase angle zero-point
c     offsets and magnitude (or flux) zero-point offsets. 
c
c     User also furnishes for the final, unique output LC:  instant of
c     reference, optional time unit for individual recorded instants of time,
c     LC period of the body, overall zero-point corrections for phase angle
c     and magnitude (or flux).
c
c     Each point from individual LCs is labelled; the user may use the
c     appropriate label to distinguish the point/LC origin in the
c     corresponding plotting macro.
c
c     The user may use more than 1 calibration object, in which case the
c     individual fluxes and seeing information are averaged. If given,
c     S/N ratios are also taken into acount for deriving individual
c     point errors in the magnitude differences and fluxes ratios. 
c
c     The user furnishes the file columns of the necessary values for
c     the computations, so the LCs may have been derived from any photometric
c     package.
c
c     The program returns the composed LC, including the phase angle
c     according to the period given. LC is composed of JD, optional instant
c     of recorded time (in user-chosen units of time), target/calibration
c     relative brightness differences in magnitude and flux ratio, errors
c     (computed from S/N ratios, if given), seeing and air masses. 
c
c
c     In this version, elimination cutoffs based on the computed error of the
c     magnitude differences are also available.
c     
c     
c
c     Last update:   M. Assafin - 04/May/2009
c
c



      IMPLICIT REAL *8 (A-H,O-Z)

      dimension lab(100000),air(100000),daju(100000),time(100000),
     ?teta(100000),tfwhm(100000),cfwhm(100000),dmg(100000),rat(100000),
     ?tflux(100000),cflux(100000),erat(100000),edmg(100000)

      dimension icflux(100),icfwhm(100),icsnr(100),col(100),ior(100000)


      character*50 output,lcurve,infits,ids9,ireg
      character*1  iast,iver

      hmsgms(a,b,c)=a+b/60.d0+c/3600.d0

c
c     Auxiliary values
c

      zero=0.d0

      it=1

      nptos=100000
      idimc=100
      idimcol=100

c

      do i=1,idimc
      icflux(i)=0
      icfwhm(i)=0
      icsnr(i)=0
      enddo

      do i=1,idimcol
      col(i)=0.d0
      enddo


      do i=1,nptos
      lab(i)=0
      air(i)=0.d0
      daju(i)=0.d0
      time(i)=0.d0
      teta(i)=0.d0
      tfwhm(i)=0.d0
      cfwhm(i)=0.d0
      cflux(i)=0.d0
      tflux(i)=0.d0
      dmg(i)=0.d0
      rat(i)=0.d0
      erat(i)=0.d0
      edmg(i)=0.d0
      ior(i)=0
      enddo



c
      iast=''
      iast='*'

      kah=0
      kam=0

c

      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '
      write (*,1)
 1    format (15x,'Light curve composition')
c
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ' '


c
c     Reads input file in loop
c

      open (1,file='PRAIA_light_curve_10.dat')

      output=''
      read (1,*) output

      write (*,2) 
 2    format (1x,'Output LC -> ',$)

      write (*,6) output
      write (*,*) ' '
      write (*,*) ' '


      read (1,*) tfwcut
      read (1,*) cfwcut

      read (1,*) tmdcut

      read (1,*) allpha

      read (1,*) kallzpb
      read (1,*) allzpb

      read (1,*) kallotu

      read (1,*) kallrt

c
c     Computes Julian Dates of reference for output LC
c

      if (kallrt.eq.3) then 
      read (1,*) djmall
      djall=djmall+2400000.5D0
      endif

      if (kallrt.eq.2) then
      read (1,*) djall
      djmall=djall-2400000.5D0
      endif

      if (kallrt.eq.1) then
      read (1,*) ah,am,as,idia,imes,iano
      fd=hmsgms(ah,am,as)/24.d0
      call iau_CAL2JD (iano,imes,idia,djm0r,djmall,iflag)
      djmall=djmall+fd
      djall=djmall+djm0r
      endif


c
c     Takes the period of the LC and converts it to Julian Days and fraction
c     of Julian days
c

      read (1,*) ipe

      read (1,*) period

      if (ipe.eq.1) period=hmsgms(zero,zero,period)/24.d0
      if (ipe.eq.2) period=hmsgms(zero,period,zero)/24.d0
      if (ipe.eq.3) period=hmsgms(period,zero,zero)/24.d0

      if (ipe.eq.5) period=period*30.5d0
      if (ipe.eq.6) period=period*365.25d0


c
c     Number of individual contributing light curves
c

      read (1,*) lcnum

c

      ip=0

c
c     Here it starts the loop to read the individual LC data
c

      do 50 lc=1,lcnum

      
      read (1,*) 

      lcurve=''

      read (1,*) lcurve

      read (1,*) label
      
      read (1,*) offang


      read (1,*) kzpb
      read (1,*) zpb

      read (1,*) ktu

      if (ktu.eq.7) then
      read (1,*) djr
      endif

      if (ktu.eq.8) then
      read (1,*) djmr
      djr=djmr+2400000.5D0
      endif

      if (ktu.le.6) then
      read (1,*) ah,am,as,idia,imes,iano
      fd=hmsgms(ah,am,as)/24.d0
      call iau_CAL2JD (iano,imes,idia,djm0r,djmr,iflag)
      djmr=djmr+fd
      djr=djmr+djm0r
      endif

c
c     Associates the needed input values wrt the indivivual LC file collumns
c

      read (1,*) itime
      read (1,*) iair
      read (1,*) itflux
      read (1,*) itsnr  
      read (1,*) itfwhm
      read (1,*) ic

      do i=1,ic
      read (1,*) icflux(i)
      read (1,*) icsnr(i)
      read (1,*) icfwhm(i)
      enddo

c
      icol=max0(itime,iair,itflux,itsnr,itfwhm)

      do i=1,ic
      if (icol.lt.icflux(i)) icol=icflux(i)
      if (icol.lt.icsnr(i)) icol=icsnr(i)
      if (icol.lt.icfwhm(i)) icol=icfwhm(i)
      enddo

c

c     call system ('clear')


      write (*,5)
 5    format (1x,'Proccessing LC -> ',$)
      write(*,6) lcurve
 6    format(a50)

c
c     Reads individual LC data
c

      do i=1,idimcol
      col(i)=0.d0
      enddo

c

      sa =0.d0
      dmi=0.d0
      hh =0.d0
      dia=0.d0
      dme=0.d0
      ano=0.d0

c

      open (3,file=lcurve)

      do 45 iiii=1,100000000
      read (3,*,end=47) (col(j),j=1,icol)

      ip=ip+1

      lab(ip)=label

c
c     Converts the LC time unit to JD
c

      tempo=col(itime)

      if (ktu.eq.8) daju(ip)=tempo+2400000.5D0
      if (ktu.eq.7) daju(ip)=tempo

      if (ktu.eq.1) sa =tempo
      if (ktu.eq.2) dmi=tempo
      if (ktu.eq.3) hh =tempo
      if (ktu.eq.4) dia=tempo
      if (ktu.eq.5) dme=tempo
      if (ktu.eq.6) ano=tempo

      if (ktu.le.6) then
      daju(ip)=djr+hmsgms(hh,dmi,sa)/24.d0+dia+dme*30.5d0+ano*365.25d0
      endif

      daj=daju(ip)-djall

c
c     Computes individual times in user-choosen time units
c


      if (kallotu.eq.1) time(ip)=daj*86400.d0
      if (kallotu.eq.2) time(ip)=daj*1440.d0
      if (kallotu.eq.3) time(ip)=daj*24.d0
      if (kallotu.eq.4) time(ip)=daj
      if (kallotu.eq.5) time(ip)=daj/30.5d0
      if (kallotu.eq.6) time(ip)=daj/365.25d0

c
c     Computes phase angle
c

      teta(ip)=daj/period

c
c     Individual LC zero-point phase angle offset
c

      if (teta(ip).lt.0.d0) then
      teta(ip)=teta(ip)-offang
      else
      teta(ip)=teta(ip)+offang
      endif

c
c     Output overall LC zero-point phase angle offset
c

      
      if (teta(ip).lt.0.d0) then
      teta(ip)=teta(ip)-allpha
      else
      teta(ip)=teta(ip)+allpha
      endif

c

      itet=teta(ip)
      teta(ip)=teta(ip)-itet

      if (teta(ip).lt.0.d0) teta(ip)=teta(ip)+1.d0

c
c     Picks up air mass
c

      if (iair.eq.0) then
      air(ip)=1.d0
      else
      air(ip)=col(iair)
      endif


c
c     Picks up target flux, S/N ratio and fwhm
c

      tflux(ip)=col(itflux)

c

      tsnr=0.d0
      if (itsnr.eq.0) then
      tsnr=-1.d0
      else
      tsnr=col(itsnr)
      endif

c

      if (itfwhm.eq.0) then
      tfwhm(ip)=99.99d0
      else
      tfwhm(ip)=col(itfwhm)
      endif


c
c     Picks up calibration objects SNR ratios and fwhms
c     More than 1 calibration object, picks up average fwhm.
c

      cfw=0.d0

      do i=1,ic
      j=icfwhm(i)
      if (j.eq.0) then
      cfw=cfw+99.99d0
      else
      cfw=cfw+col(j)
      endif
      enddo

      cfwhm(ip)=cfw/ic

c
c     Picks up calibration objects total flux
c

      cfx=0.d0

      do i=1,ic
      j=icflux(i)
      cfx=cfx+col(j)
      enddo

      cflux(ip)=cfx/ic


c
c     Computes brightness relative differences between target flux and
c     total calibration flux  
c
c

      ratio=tflux(ip)/cflux(ip)

c
c     Individual LC zero-point brightness correction 
c


      if (kzpb.eq.1) then
      ratiof=10.d0**(zpb/2.5d0)
      ratiof=1.d0/ratiof
      endif

      if (kzpb.eq.2) ratiof=zpb

      ratio=ratio*ratiof


c
c     Output overall LC zero-point brightness correction 
c


      if (kallzpb.eq.1) then
      ratioa=10.d0**(allzpb/2.5d0)
      ratioa=1.d0/ratioa
      endif

      if (kallzpb.eq.2) ratioa=allzpb

      ratio=ratio*ratioa

c
c     Computes errors in magnitude difference and flux ratio. 
c     More than 1 calibration object, composes magnitude and flux ratio
c     errors
c

      csn=0.d0

      do i=1,ic
      j=icsnr(i)
      if (j.eq.0) then
      tsnr=-1.d0
      else
      csn=csn+(1.d0/col(j)**2)
      endif
      enddo

      if (tsnr.lt.0.d0) then
      erat(ip)=0.d0
      edmg(ip)=0.d0
      else
      erat(ip)=dsqrt(csn/ic+(1.d0/tsnr**2)/it)
      edmg(ip)=(2.5d0/dlog(10.d0))*erat(ip)
      erat(ip)=ratiof*ratioa*erat(ip)
      endif


c
c     Output brightness relative differences between target and calibration
c     in flux and magnitude
c

      difmag=-2.5d0*dlog10(ratio)

      rat(ip)=ratio
      dmg(ip)=difmag


 45   continue

c

 47   close (3)

c

 50   continue

c

      close (1)

c
c     Ordering by JD
c

      do i=1,ip
      ior(i)=i
      enddo

      call ordem (ip,ior,daju)

c
c     Writing output LC to file
c

      open(7,file=output)

      do 70 j=1,ip

      i=ior(j)

      if (tfwhm(i).gt.tfwcut) go to 70
      if (cfwhm(i).gt.cfwcut) go to 70

      if (edmg(i).gt.tmdcut) go to 70

      write (7,60) lab(i),air(i),daju(i),time(i),teta(i),tfwhm(i),
     ?cfwhm(i),dmg(i),edmg(i),rat(i),erat(i),tflux(i),cflux(i)
 60   format(1x,i3,1x,f5.3,2(1x,f16.8),1x,f16.14,2(1x,f6.3),2(1x,f9.5),
     ?2(1x,f16.11),2(1x,f16.5))

 70   continue
c

      close (7)

c

      write (*,*)
      write (*,*)
      write (*,*) ' Proccess terminated. Status ok.'
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
C     SUBROUTINE ORDEM (N,IOR,VAL)
C
C
C     Description of parameters
C
C       N      - number of points to be ordered
C       indx   - increasing order numbering of array "arr"
C       arr    - data array itself, NOT ORDERED
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

      SUBROUTINE ordem(n,indx,arr)
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
c       if(jstack.gt.NSTACK)pause 'NSTACK too small in ordem'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

