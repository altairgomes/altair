c
c
c     PRAIA_check_ref_stars
c
c
c     Lists fits images which reductions presented less than a given number
c     of reference stars.
c
c     User needs to furnish the xy list and minimum number of reference stars in the
c     input ".dat" file. 
c
c
c     Last update: M. Assafin  - 09/Sep/2011
c


      
      IMPLICIT REAL *8 (A-H,O-Z)



      character*150 list,input,output
      character*50 infits
      character*20 ichfil,mchobj
 

 1    continue

      read (*,*,end=40) n,list,output

      open (30,file=output)

      open (10,file=list)


 5    continue

      input=''
      read (10,*,end=30) input

      open (20,file=input)


 7    read (20,10,err=7) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny

 10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))

      close (20)

      if (nfinau.le.n) write (30,20) nfinau,infits

 20   format(1x,i3,1x,a50)

      go to 5

 30   close (10)

      close (30)

      go to 1

 40   continue


      end


