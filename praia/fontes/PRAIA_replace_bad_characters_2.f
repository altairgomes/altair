c
c
c    PRAIA_replace_bad_characters
c
c
c    Replaces "*" or "NaN" characters by numerical ones in certain non-fundamental columns in
c    the xy files due to overflow in the format. The remaining information is preserved,
c    so that the entry is not dropped by certain analysis tools.
c
c
c
c    In this version, bad characters (overflows, NaNs) are checked for:

c
c
c    - computed magnitudes (UCAC2 and 2MASS-converted system)
c    - amplitude and the background of the Gaussian profile
c    - standard error of (RA,DEC) reductions
c    - mean error or (RA,DEC) reductions
c    - (O-C)s 
c    - RA, Dec
c
c
c
c
c     Command line:
c
c
c    PRAIA_replace_bad_characters <  PRAIA_replace_bad_characters.dat > log_bad_xy_files
c
c
c    where PRAIA_replace_bad_characters.dat is the list of xy files.
c
c
c
c     Last update: M. Assafin  - 17/Sep/2011
c
c
c

     
      IMPLICIT REAL *8 (A-H,O-Z)



      character*150 input
      character*50 infits
      character*20 ichfil,mchobj
      character*9 kaltu,kfgcc
      character*6 kudmg,kudmg2
      character*6 kalsiu,kdesiu
      character*13 kra,kde
      character*6 kerau,kedeu,kalfsi,kdelsi
      


 1    continue

      input=''
      read (*,*,end=45) input

      open (10,file='PRAIA_replace_bad_characters.tmp')

      open (20,file=input)

c
      
      key=0

 5    continue

 7    read (20,10,end=30) xob,yob,seng,kaltu,kfgcc,fumag,fumag2,
     ?xmgu,kudmg,kudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,kerau,kedeu,kalfsi,kdelsi,
     ?nstaru,nfinau,kalsiu,kdesiu,ktirau,kra,kde,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny

 10   format(2(1x,f7.2),1x,f5.3,2(1x,a9),3(1x,f6.3),2(1x,a6),8(1x,f6.3),
     ?4(1x,f7.3),2(1x,f6.3),4(1x,a6),2(1x,i4),2(1x,a6),1x,i4,2(1x,a13),
     ?1x,i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,
     ?a50,1x,a20,2(1x,i5))


  

      if (kaltu.eq.'*********') then
      altu=-1.d0
      key=1
      else
      read (kaltu,*) altu
      endif

      if (kfgcc.eq.'*********') then
      fgcc=-1.d0
      key=1
      else
      read (kfgcc,*) fgcc
      endif

      if (kudmg.eq.'   NaN' .or.  kudmg.eq.'   nan') then
      cudmg=99.999d0
      key=1
      else
      read (kudmg,*) cudmg
      endif

      if (kudmg2.eq.'   NaN' .or.  kudmg2.eq.'   nan') then
      cudmg2=99.999d0
      key=1
      else
      read (kudmg2,*) cudmg2
      endif



      if (kerau.eq.'   NaN' .or.  kerau.eq.'   nan') then
      erau=99.999d0
      key=1
      else
      read (kerau,*) erau
      endif


      if (kedeu.eq.'   NaN' .or.  kedeu.eq.'   nan') then
      edeu=99.999d0
      key=1
      else
      read (kedeu,*) edeu
      endif

      if (kalfsi.eq.'   NaN' .or.  kalfsi.eq.'   nan') then
      alfsiu=99.999d0
      key=1
      else
      read (kalfsi,*) alfsiu
      endif

      if (kdelsi.eq.'   NaN' .or.  kdelsi.eq.'   nan') then
      delsiu=99.999d0
      key=1
      else
      read (kdelsi,*) delsiu
      endif


      if (kalsiu.eq.'******'.or.kalsiu.eq.'   NaN'.or.kalsiu.eq.'   nan'
     ?) then
      alsiu=99.999d0
      key=1
      else
      read (kalsiu,*) alsiu
      endif

      if (kdesiu.eq.'******'.or.kdesiu.eq.'   NaN'.or.kdesiu.eq.'   nan'
     ?) then
      desiu=99.999d0
      key=1
      else
      read (kdesiu,*) desiu
      endif


      if (kra.eq.'          NaN' .or.  kra.eq.'          nan') then
      ra=99.999d0
      key=1
      else
      read (kra,*) ra
      endif


      if (kde.eq.'          NaN' .or.  kde.eq.'          nan') then
      de=99.999d0
      key=1
      else
      read (kde,*) de
      endif





      write (10,20) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny

 20   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))

      go to 5

 30   close (20)
      close (10)

      open (10,file='PRAIA_replace_bad_characters.tmp')
      open (20,file=input)

 35   continue

      read (10,20,end=40) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny

      write (20,20) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny

      go to 35

 40   close (20)
      close (10)

      if (key.eq.1) write (*,*) input


      go to 1

 45   continue

      call system ('rm PRAIA_replace_bad_characters.tmp')


      end


