c
c
c      Dec_hex
c
c
c      Converts (RA,DEC) from decimal fraction to hexadecimal format
c      and vice-versa.
c
c      Last update: Marcelo Assafin - 12 November 2009
c   
c
c


      IMPLICIT REAL *8 (A-H,O-Z)



      character*250 ler
      character*1 isig

      hmsgms(x,y,z)=x+y/60.d0+z/3600.d0

c

      idim=250

      dzero=0.d0
c

c
c     Reads input (RA,DEC) coordinates
c


      write (*,*)
      write (*,*) '(RA,Dec): decimal <--> hexadecimal format'
      write (*,*)  
      write (*,*)  

      ler=''

      write (*,10)
 10   format (1x,'Input coordinates = ',$)
      read (*,12) ler
 12   format(a250)

      write (*,*)  
      write (*,*)  

      write (*,*) '    ',ler

c
c     iflag = 1 :     hh.hhh  +/-dd.ddd  ->  hh mm ss.sss +/-dg mm ss.ss 
c
c
c     iflag = 2 :     hh mm ss.sss +/-dg mm ss.ss -> hh.hhh  +/-dd.ddd 
c
c

      iflag=2

      do i=1,idim
      if (ler(i:i).ne.' ') go to 20
      enddo

 20   i1=i

      icont=0

      do i=i1,idim-1
      if (ler(i:i).ne.' '.and.ler(i+1:i+1).eq.' ') icont=icont+1
      enddo

      if (icont.lt.3) iflag=1


      if (iflag.eq.2) go to 30


c     
c     iflag = 1 case:   hh.hhh  +/-dd.ddd  ->  hh mm ss.sss +/-dg mm ss.ss
c



      read (ler,*) ra,de


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

      write(*,25) iah,iam,sa,isig,idg,idm,ds
 25   format(1x,2(i2.2,1x),f9.6,2x,a1,2(i2.2,1x),f8.5)

      go to 50

c
c     iflag = 2 case:   hh mm ss.sss +/-dg mm ss.ss -> hh.hhh  +/-dd.ddd 
c


 30   continue


c
c     Checks declination sign
c

      isig='+'
      do i=i1,idim
      if (ler(i:i).eq.'-') isig='-'
      enddo


c
c     Case 1:  hh mm.mm +/-dg mm.mm  -> hh.hhh  +/-dd.ddd
c

      if (icont.eq.4) then

      read (ler,*) ah,am,dg,dm

      as=dzero
      ds=dzero

      endif

c
c     Case 2:  hh mm ss.sss +/-dg mm ss.sss  -> hh.hhh  +/-dd.ddd
c


      if (icont.eq.6) then

      read (ler,*) ah,am,as,dg,dm,ds

      endif

      dg=dabs(dg)

      ra=hmsgms(ah,am,as)
      de=hmsgms(dg,dm,ds)

      if (isig.eq.'-') de=-de

      write (*,40) ra,de

 40   format(1x,f18.15,2x,f19.15)

c

 50   continue

      write (*,*)

      end

