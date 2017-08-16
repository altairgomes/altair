c
c     Program PRAIA_global_reduction_20
c
c
c
c     Purpose
c
c
c     Given fields reduced with PRAIA, redoes the reductions of these fields
c     using the principle of Global Reduction (Assafin et al. 1997, AJ, and
C     references therein) for CCD fields and tangent plane astrometry.
c
c
c     First step must be individual CCD reductions performed by PRAIA.
c
c     Next, all xy fields are re-reduced with 3rd degree polynom using all
c     stars as reference stars. Stars common to 2 or more CCD fields have
c     their positions averaged and the process continues untill star positions
c     don't vary by more than a (user-given) value. 
c
c     Then, all positions are rotated in (X,Y,Z) toward the UCAC2 reference
c     frame, using all common UCAC2 stars available. 
c
c     Finally, all CCD fields are last reduced using this last set of star
c     positions, so as to place moving objects on the final reference frame
c     too.
c 
c
c     In this version, besides a rigid rotation, all positions may be
c     transformed toward the UCAC2 reference frame, by using all common UCAC2
c     stars available, with the tangent plane techinique. (user decides)
c
c     Here in this version, the transform toward the UCAC2 reference frame is
c     made at (the end of) each step of the process, instead of being done
c     only at the end of the procedure.
c
c
c     In this version, there are separate (O-C) cutoffs for the transform
c     to the UCAC2 and for the individual CCD reductions in the iterative
c     proccesss. 
c
c
c     In this version, the user defines if the target will be or not averaged
c     in the global reduction process. If not, it will be treated as a
c     distinct (different) object in each individual CCD, and no unique
c     averaged position will be produced. This is recomended in the case the
c     target is a fast moving object (example: asteroid or natural sattelite),
c     because than it may perturb the overlaping solution in the intermediary
c     global reduction proccess. In this case, the individual target positions
c     should be excluded from the averaging step in the intermediary global
c     reduction proccess. The exclude option does not affect the correction
c     during the, and at the end of the global reduction, over the individual
c     target positions obtained from the re-reduction of the individual CCD
c     frames, using the (step by step) renewed catalog of all averaged
c     positions. 
c
c
c
c      Last update: Marcelo Assafin - 15/Feb/2013
c   
c
c


      IMPLICIT REAL *8 (A-H,O-Z)

      parameter(stdin=5)

      character*50 listaxy,filexy


c
c     Determines dimensions for frame and star number of xy fields
c

      read (5,*)
      read (5,*)
      read (5,1) listaxy
 1    format(a50)
      rewind (5)


      jdfra=0
      jdnest=0
      jall=0

      open (1,file=listaxy)

 2    read (1,1,end=15) filexy

      jdfra=jdfra+1

      open (2,file=filexy)

      i=0

 5    read(2,*,end=10)

      i=i+1
      jall=jall+1

      go to 5

 10   if (i.gt.jdnest) jdnest=i

      close (2)

      go to 2

 15   close(1)

c
c     Furnishes dimensions for main program
c     

      idfra=jdfra
      idnest=jdnest
      iall=jall


      call main (idfra,idnest,iall)


      end


c
c     Main program
c


      subroutine main(idfra,idnest,iall)

      IMPLICIT REAL *8 (A-H,O-Z)

c     parameter(idfra=jdfra,idnest=jdnest,iall=jall)

      parameter(ichas=100,iu2z=288,iu2i=240,iu4z=900,iu4i=1440,
     ?idin50=50)


      dimension x(idnest,idfra),y(idnest,idfra),seng(idnest,idfra),
     ?altu(idnest,idfra),fgcc(idnest,idfra),fumag(idnest,idfra),
     ?fumag2(idnest,idfra),xmgu(idnest,idfra),cudmg(idnest,idfra),
     ?cudmg2(idnest,idfra),xmgj(idnest,idfra),xmgh(idnest,idfra),
     ?xmgk(idnest,idfra),res2mg(idnest,idfra),resmg2(idnest,idfra),
     ?ermgj(idnest,idfra),ermgh(idnest,idfra),ermgk(idnest,idfra)

      dimension pma(idnest,idfra),pmd(idnest,idfra),epma(idnest,idfra),
     ?epmd(idnest,idfra),ex(idnest,idfra),ey(idnest,idfra),
     ?erau(idnest,idfra),edeu(idnest,idfra),alfsiu(idnest,idfra),
     ?delsiu(idnest,idfra),nstaru(idnest,idfra),nfinau(idnest,idfra),
     ?alsiu(idnest,idfra),desiu(idnest,idfra),ktirau(idnest,idfra),
     ?raest(idnest,idfra),deest(idnest,idfra),kuth(idnest,idfra),
     ?kutm(idnest,idfra),zut(idnest,idfra),kutano(idnest,idfra),
     ?kutmes(idnest,idfra),kutdia(idnest,idfra),daju(idnest,idfra),
     ?iexps(idnest,idfra),nx(idnest,idfra),ny(idnest,idfra),
     ?index(idnest,idfra)

      dimension id(iall),xob(iall),yob(iall),talf(iall),tdel(iall),
     ?tadj(iall)

      dimension cora(iall),code(iall),codj(iall),cseng(iall),
     ?codmg(iall),codmg2(iall),cxmgj(iall),cxmgh(iall),cxmgk(iall),
     ?copma(iall),copmd(iall),cxmgu(iall),cora2(iall),code2(iall)

      dimension coex(iall),coey(iall),cerau(iall),cedeu(iall),
     ?alfsic(iall),delsic(iall),finac(iall),alsiuc(iall),desiuc(iall),
     ?tirac(iall),numcom(iall)

      dimension oldra(iall),oldde(iall),jstars(idfra),ranew(iall),
     ?denew(iall),ucra(iall),ucde(iall)

      dimension inuu2(iu2z,iu2i),ir1u2(iu2z,iu2i),ir2u2(iu2z,iu2i)
      dimension inuu4(iu4z,iu4i),ir1u4(iu4z,iu4i),ir2u4(iu4z,iu4i)



      character*170 linha
      character*50 listxy,target,catalog,mfits,ilog
      character*100 ifxy1(idfra),ifxy2(idfra),name,infits(idnest,idfra)
      character*50 tarall,tarsgl
      character*50 uraiz
      character *61 u2ind,u4ind
      character*20 iobalv(idnest,idfra),ichfil(idnest,idfra)
      character*12 iext
      character*1 menos,ibr,isig


c    ?,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,
c    ?res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,
c    ?alfsiu,delsiu,alsiu,desiu,raest,deest,zut,daju,ichfil,infits,
c    ?iobalv,nstaru,nfinau,ktirau,kuth,kutm,kutano,kutmes,kutdia,iexps,
c    ?nx,ny,index



      data menos/'-'/
      data ibr  /' '/

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0



c
C     Initializing data
C
c


      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c
c     Auxliliary data
c

      idiobs=iall
      kfiles=idfra
      kstars=idnest
      dnove=99.999d0

c
c     Emptying vectors
c

      do i=1,idiobs
      id(i)=0
      oldra(i)=0.d0
      oldde(i)=0.d0
      ranew(i)=0.d0
      denew(i)=0.d0
      ucra(i)=0.d0
      ucde(i)=0.d0
      cora(i)=0.d0
      code(i)=0.d0
      codj(i)=0.d0
      cseng(i)=0.d0
      codmg(i)=0.d0
      codmg2(i)=0.d0
      cxmgj(i)=0.d0
      cxmgh(i)=0.d0
      cxmgk(i)=0.d0
      copma(i)=0.d0
      copmd(i)=0.d0
      coex(i)=0.d0
      coey(i)=0.d0
      cerau(i)=0.d0
      cedeu(i)=0.d0
      alfsic(i)=0.d0
      delsic(i)=0.d0
      finac(i)=0.d0
      alsiuc(i)=0.d0
      desiuc(i)=0.d0
      tirac(i)=0.d0
      numcom(i)=0
      talf(i)=0.d0
      tdel(i)=0.d0
      tadj(i)=0.d0
      enddo


      do i=1,kfiles

      jstars(i)=0

      do j=1,kstars

      index(j,i)=0

      x(j,i)=0.d0
      y(j,i)=0.d0
      seng(j,i)=0.d0
      altu(j,i)=0.d0
      fgcc(j,i)=0.d0
      fumag(j,i)=0.d0
      fumag2(j,i)=0.d0
      xmgu(j,i)=0.d0
      cudmg(j,i)=0.d0
      cudmg2(j,i)=0.d0
      xmgj(j,i)=0.d0
      xmgh(j,i)=0.d0
      xmgk(j,i)=0.d0
      res2mg(j,i)=0.d0
      resmg2(j,i)=0.d0
      ermgj(j,i)=0.d0
      ermgh(j,i)=0.d0
      ermgk(j,i)=0.d0
      pma(j,i)=0.d0
      pmd(j,i)=0.d0
      epma(j,i)=0.d0
      epmd(j,i)=0.d0
      ex(j,i)=0.d0
      ey(j,i)=0.d0
      erau(j,i)=0.d0
      edeu(j,i)=0.d0
      alfsiu(j,i)=0.d0
      delsiu(j,i)=0.d0
      nstaru(j,i)=0.d0
      nfinau(j,i)=0.d0
      alsiu(j,i)=0.d0
      desiu(j,i)=0.d0
      ktirau(j,i)=0
      raest(j,i)=0.d0
      deest(j,i)=0.d0
      kuth(j,i)=0
      kutm(j,i)=0
      zut(j,i)=0.d0
      kutano(j,i)=0
      kutmes(j,i)=0
      kutdia(j,i)=0
      daju(j,i)=0.d0
      iexps(j,i)=0
      ichfil(j,i)=''
      infits(j,i)=''
      iobalv(j,i)=''
      nx(j,i)=0
      ny(j,i)=0

      enddo
      enddo


c
c     Reading input batch file
c

      write (*,*) 
      write (*,*) 
      write (*,*) 'PRAIA Global Reduction'
      write (*,*) 
      write (*,*) 

      call system('rm -f PRAIA_global_reduction.tmp')

c


c     open (5,file='PRAIA_global_reduction_14.dat')

 1    format(a50)
      read (5,1) uraiz
      read (5,*) keyca
      read (5,1) listxy
      read (5,1) target
      read (5,*) itakey
      read (5,1) ilog
      read (5,1) catalog
      read (5,1) tarall
      read (5,1) tarsgl
      read (5,2) iext
 2    format(a12)

      read (5,*) ktrans
      read (5,*) ecom
      read (5,*) tcom
      read (5,*) ucom
      read (5,*) corte
      read (5,*) cortu
      read (5,*) csigma
      read (5,*) maiter

      read (5,*) ngrau
      read (5,*) ngrau3
      read (5,*) ngrau5

      read (5,*) jgrau
      read (5,*) jgrau3
      read (5,*) jgrau5

c     close (5)

      rewind (5)

c

      open (46,file=ilog)

      write (46,*) 
      write (46,*) 
      write (46,*) 'PRAIA Global Reduction'
      write (46,*) 
      write (46,*) 


c     open (5,file='PRAIA_global_reduction_14.dat')

      do i=1,24
      read (5,6) linha
 6    format(a170)
      write (46,6) linha
      enddo

c     close (5)

      write (46,*)
c

      write(*,*) 'Step 1. Re-reduction of xy fields'
      write (*,*) 

      write(46,*) 'Step 1. Re-reduction of xy fields'
      write (46,*) 

c

      open (77,file=listxy)

c

      do i=1,kfiles

      ifxy1(i)=''

      read (77,7,end=8) ifxy1(i)
 7    format(a100)

      enddo

 8    close (77)

      nfxy=i-1

      if (nfxy.eq.kfiles) write (*,*) 'Maximum No. of files reached: ',
     ?nfxy

      if (nfxy.eq.kfiles) write (46,*) 'Maximum No. of files reached: ',
     ?nfxy


c
c     Mounts output xy file names
c

      do i=1,nfxy

      ifxy2(i)=''
      name=''
      name=ifxy1(i)

      do k=ichas,1,-1
      if (name(k:k).ne.ibr) go to 9 
      enddo

 9    name(k-11:k)=iext
      
      ifxy2(i)=name

      enddo

c
c     Allocation of data from all xy files into memory
c

      write (*,*) 
      write (46,*) 

      mast=0
      mist=10000000

      do i=1,nfxy

      write (*,10) i,nfxy,ifxy1(i)
      write (46,10) i,nfxy,ifxy1(i)
 10   format('Reading file ',2(1x,i5),1x,a100)

      open (1,file=ifxy1(i))

      do j=1,kstars

 12   continue

      read (1,13,err=12,end=14) x(j,i),y(j,i),seng(j,i),altu(j,i),
     ?fgcc(j,i),fumag(j,i),fumag2(j,i),xmgu(j,i),cudmg(j,i),cudmg2(j,i),
     ?xmgj(j,i),xmgh(j,i),xmgk(j,i),res2mg(j,i),resmg2(j,i),ermgj(j,i),
     ?ermgh(j,i),ermgk(j,i),pma(j,i),pmd(j,i),epma(j,i),epmd(j,i),
     ?ex(j,i),ey(j,i),erau(j,i),edeu(j,i),alfsiu(j,i),delsiu(j,i),
     ?nstaru(j,i),nfinau(j,i),alsiu(j,i),desiu(j,i),ktirau(j,i),
     ?raest(j,i),deest(j,i),kuth(j,i),kutm(j,i),zut(j,i),kutano(j,i),
     ?kutmes(j,i),kutdia(j,i),daju(j,i),iexps(j,i),ichfil(j,i),
     ?infits(j,i),iobalv(j,i),nx(j,i),ny(j,i)

 13   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),
     ?13(1x,f6.3),4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,
     ?2(1x,f13.9),1x,i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,
     ?i4,2x,a20,2x,a50,1x,a20,2(1x,i5))


      enddo

 14   close (1)

      jstars(i)=j-1

      if (jstars(i).gt.mast) mast=jstars(i)
      if (jstars(i).lt.mist) mist=jstars(i)

      enddo


c
c     Allocation of target data into memory
c

      open (1,file=target)


      do l=1,idiobs
      read(1,15,end=16) iah,iam,as,isig,idg,idm,ds,tadj(l)
  15  format(1x,i2,1x,i2,1x,f9.6,1x,a1,i2,1x,i2,1x,f8.5,1x,f16.8)

      talf(l)=hmsgms(iah,iam,as)
      tdel(l)=hmsgms(idg,idm,ds)
      if (isig.eq.menos) tdel(l)=-tdel(l)

      enddo

 16   close (1)

      ntarg=l-1



c
c     Indexing common/non-commom stars from all xy files
c


      call indexx (nfxy,jstars,ifxy1,ecom,cora,code,numcom,id,talf,tdel,
     ?ntarg,tcom,itakey,nobs,idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,
     ?fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,
     ?raest,deest,zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,
     ?kuth,kutm,kutano,kutmes,kutdia,iexps,nx,ny,index)


c
      write(*,*) 
      write(*,*) 
      write(*,*) 
      write(*,*) 'Total number of field stars to be proccessed: ',nobs
      write(*,*) 'Maximum number of stars in a field          : ',mast  
      write(*,*) 'Minimum number of stars in a field          : ',mist  
      write(*,*) 'Total number of fields to be proccessed:    : ',nfxy  
      write(*,*) 
      write(*,*) 
      write(*,*) 


      write(46,*) 
      write(46,*) 
      write(46,*) 
      write(46,*)'Total number of field stars to be proccessed: ',nobs  
      write(46,*)'Maximum number of stars in a field          : ',mast  
      write(46,*)'Minimum number of stars in a field          : ',mist  
      write(46,*)'Total number of fields to be proccessed:    : ',nfxy  
      write(46,*) 
      write(46,*) 
      write(46,*) 


c
c     First update. Picks up stars from all xy fields and average some
c     parameters (positions, JD, mag, etc) for common stars.
c


      write (*,*) 
      write (*,*) 
      write (46,*) 
      write (46,*)

      mmmo=1

      write (*,*) 'Iter. No. ',mmmo
      write(*,*) 'Updating star position, JD, mags, etc from all xy file
     ?s'

      write (46,*) 'Iter. No. ',mmmo
      write(46,*) 'Updating star position, JD, mags, etc from all xy fil
     ?es'


      call update (jstars,nfxy,ifxy1,nobs,cora,code,cora2,code2,
     ?codj,cseng,cxmgu,codmg,codmg2,cxmgj,cxmgh,cxmgk,copma,copmd,coex,
     ?coey,cerau,cedeu,alfsic,delsic,finac,alsiuc,desiuc,tirac,numcom,
     ?idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,
     ?cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,
     ?epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,deest,
     ?zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,kutm,
     ?kutano,kutmes,kutdia,iexps,nx,ny,index)



c
c     Defines area for reference catalog star extraction
c

      alfmin=1.d14
      delmin=1.d14
      alfmax=-1.d14
      delmax=-1.d14

      do i=1,nobs
      if (cora(i).lt.alfmin) alfmin=cora(i)
      if (code(i).lt.delmin) delmin=code(i)
      if (cora(i).gt.alfmax) alfmax=cora(i)
      if (code(i).gt.delmax) delmax=code(i)
      enddo


c
c     Picks up common UCAC2 reference stars 
c

      write (*,*)
      write (*,*) 'Extracting UCAC stars.'
      write (*,*)

      write (46,*)
      write (46,*) 'Extracting UCAC stars.'
      write (46,*)


c
c     UCAC2, no speedy extraction
c

      if (keyca.eq.1) then

      call ucac2 (uraiz,ucom,alfmin,alfmax,delmin,delmax,nobs,cora,code,
     ?codj,ucra,ucde,id,oldra,oldde,rac,dec,nok,rmx,rmy,rmx2,rmy2,iall)

      endif


c
c     UCAC2, speedy extraction
c

      if (keyca.eq.2) then

      u2ind=''
      u2ind=uraiz


      do l=1,idin50
      if (u2ind(l:l).eq.' ') go to 17
      enddo

 17   u2ind(l:l+10)='u2index.txt'


      open (3,file=u2ind)

      do i=1,10
      read (3,*)
      enddo

      do i=1,iu2z
      do j=1,iu2i
      read(3,*) nsbin,naz
      inuu2(i,j)=nsbin
      ir1u2(i,j)=naz-nsbin+1
      ir2u2(i,j)=ir1u2(i,j)+nsbin-1
      enddo
      enddo

      close (3)




      call ucac2s (iu2z,iu2i,inuu2,ir1u2,ir2u2,uraiz,ucom,alfmin,alfmax,
     ?delmin,delmax,nobs,cora,code,codj,ucra,ucde,id,oldra,oldde,rac,
     ?dec,nok,rmx,rmy,rmx2,rmy2,iall)

      endif


c
c     UCAC4, no speedy extraction
c


      if (keyca.eq.3) then

      call ucac4 (uraiz,ucom,alfmin,alfmax,delmin,delmax,nobs,cora,code,
     ?codj,ucra,ucde,id,oldra,oldde,rac,dec,nok,rmx,rmy,rmx2,rmy2,iall)

      endif


c
c     UCAC4, speedy extraction
c


      if (keyca.eq.4) then

      u4ind=''
      u4ind=uraiz


      do l=1,idin50
      if (u4ind(l:l).eq.' ') go to 18
      enddo

 18   u4ind(l:l+10)='u4index.asc'


      open (3,file=u4ind)


      do i=1,iu4z
      do j=1,iu4i
      read(3,*) naz,nsbin
      inuu4(i,j)=nsbin
      ir1u4(i,j)=naz+1
      ir2u4(i,j)=ir1u4(i,j)+nsbin-1
      enddo
      enddo

      close (3)



      call ucac4s(iu4z,iu4i,inuu4,ir1u4,ir2u4,uraiz,ucom,alfmin,alfmax,
     ?delmin,delmax,nobs,cora,code,codj,ucra,ucde,id,oldra,oldde,rac,
     ?dec,nok,rmx,rmy,rmx2,rmy2,iall)

      endif






c
c     RA,DEC re-reduction of xy fields
c
c     Re-reductions are performed until convergence is reached
c
c

      sigma=1.d14

      mmmo=0

 20   continue



c
c     Here it starts the iteractive proccess. All xy fields have some
c     parameters from common stars averaged. Positions are averaged and
c     transformed into the UCAC2 frame at each step.
c
c


      mmmo=mmmo+1



c
c     Maximum No. of global reduction iterations
c

      if (mmmo.gt.maiter) then

      write (*,*)
      write (*,*) 'Maximum No. of global reduction iterations reached.'
      write (*,*) 

      write (46,*)
      write (46,*) 'Maximum No. of global reduction iterations reached.'
      write (46,*) 
      go to 25
      endif

c

      if (mmmo.eq.1) go to 22

      write (*,*) 
      write (*,*) 
      write (46,*) 
      write (46,*) 

      write (*,*) 'Iter. No. ',mmmo
      write(*,*) 'Updating star position, JD, mags, etc from all xy file
     ?s'

      write (46,*) 'Iter. No. ',mmmo
      write(46,*) 'Updating star position, JD, mags, etc from all xy fil
     ?es'

c
c
c     Picks up stars from all xy fields and average some values for
c     common stars at each step
c


      call update (jstars,nfxy,ifxy2,nobs,cora,code,cora2,code2,
     ?codj,cseng,cxmgu,codmg,codmg2,cxmgj,cxmgh,cxmgk,copma,copmd,coex,
     ?coey,cerau,cedeu,alfsic,delsic,finac,alsiuc,desiuc,tirac,numcom,
     ?idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,
     ?cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,
     ?epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,deest,
     ?zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,kutm,
     ?kutano,kutmes,kutdia,iexps,nx,ny,index)



 22   continue


c
c     Transformation of updated positions toward the reference catalog 
c     UCAC2. The transformation is made by rigid rotation or by the
c     tangent plane technique (user defined).
c
c     The output is temporarily stored in variables "ranew" and "denew".
c     Then, the position variables "cora" and "code" are updated.
c     
c     ktrans = 1  -> rigid rotation
c
c     ktrans = 2  -> tangent plane technique
c
c


      call transf (ktrans,cortu,nobs,cora,code,ucra,ucde,ranew,denew,id,
     ?rac,dec,jgrau,jgrau3,jgrau5,rmx,rmy,rmx2,rmy2,iall)


c
c     Retrieves back the (now new, transformed) coordinates into
c     "cora, code" variables.
c


      do i=1,nobs
      cora(i)=ranew(i)
      code(i)=denew(i)
      enddo


c
c     RA,DEC re-reduction of indivudual frames by updated "(cora,code)"
c     positions.
c

      write (*,*) 
      write (*,*) '(RA,DEC) re-reductions.'
      write (*,*) 

      write (46,*) 
      write (46,*) '(RA,DEC) re-reductions.'
      write (46,*) 


      if (mmmo.lt.2) then

      call reduc (nfxy,jstars,ifxy1,corte,ngrau,ngrau3,ngrau5,cora,code,
     ?idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,
     ?cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,
     ?epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,deest,
     ?zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,kutm,
     ?kutano,kutmes,kutdia,iexps,nx,ny,index)

      else

      call reduc (nfxy,jstars,ifxy2,corte,ngrau,ngrau3,ngrau5,cora,code,
     ?idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,
     ?cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,
     ?epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,deest,
     ?zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,kutm,
     ?kutano,kutmes,kutdia,iexps,nx,ny,index)


      endif



c
c     Internal precision of positions from (SA,SD) about the (RA,DEC) mean
c     position of common stars
c

      write (*,*)
      write (46,*)

      xvam=0.d0
      xvas=0.d0

      yvam=0.d0
      yvas=0.d0

      n=0


      do 23 i=1,nobs
      
      if (numcom(i).lt.2) go to 23

      n=n+1

      aux=cora2(i)*15.d0*3600.d0
      auy=code2(i)*3600.d0

      xvam=xvam+aux
      yvam=yvam+auy

      xvas=xvas+aux**2
      yvas=yvas+auy**2

 23   continue

      call desvio (n,xvam,xvas)
      call desvio (n,yvam,yvas)

      write (*,*) 'Internal precision: (SA,SD)           = ',xvam,yvam,n
      write (*,*) 'Sigma of Internal precision (sSA,sSD) = ',xvas,yvas  

      write (*,*)

      write (46,*)'Internal precision: (SA,SD)           = ',xvam,yvam,n
      write (46,*)'Sigma of Internal precision (sSA,sSD) = ',xvas,yvas  

      write (46,*)

c
c     Checking convergence
c

      if (mmmo.ge.2) then

      write (*,*)
      write (*,*) 'Checking convergence ...'
      write (*,*)

      write (46,*)
      write (46,*) 'Checking convergence ...'
      write (46,*)

      endif


      if (mmmo.lt.2) then

      do i=1,nobs
      oldra(i)=cora(i)
      oldde(i)=code(i)
      enddo

      go to 20      
      
      endif

c

      xvam=0.d0
      xvas=0.d0

      yvam=0.d0
      yvas=0.d0


      do i=1,nobs

      dx=(oldra(i)-cora(i))*dcos(cora(i)*grarad)*15.d0*3600.d0
      dy=(oldde(i)-code(i))*3600.d0

      xvam=xvam+dx
      xvas=xvas+dx**2

      yvam=yvam+dy
      yvas=yvas+dy**2

      enddo

c  

      call desvio (nobs,xvam,xvas)
      call desvio (nobs,yvam,yvas)



c
c     Convergence reached. 1rst step of global reduction concluded here.
c

      write (*,*)
      write(*,*)'RA  average, sigma, N old/new step : ',xvam,xvas,nobs 
      write(*,*)'DEC average, sigma, N old/new step : ',yvam,yvas,nobs 
      write (*,*)

      write (46,*)
      write(46,*)'RA  average, sigma, N old/new step : ',xvam,xvas,nobs 
      write(46,*)'DEC average, sigma, N old/new step : ',yvam,yvas,nobs 
      write (46,*)

      if (xvas.gt.csigma .or. yvas.gt.csigma) then

      do i=1,nobs
      oldra(i)=cora(i)
      oldde(i)=code(i)
      enddo

      go to 20      
      
      endif

c

 25   continue

      write (*,*) '1rst step: global reduction convergence reached.'
      write (46,*) '1rst step: global reduction convergence reached.'



c
c     Now comes the last, final transformation of the resulting positions
c     toward the reference catalog UCAC2. 
c
c     The position variables "cora" and "code" are one last time transformed.
c     
c     ktrans = 1  -> rigid rotation
c
c     ktrans = 2  -> tangent plane technique
c
c

      call transf (ktrans,cortu,nobs,cora,code,ucra,ucde,ranew,denew,id,
     ?rac,dec,jgrau,jgrau3,jgrau5,rmx,rmy,rmx2,rmy2,iall)


c
c     Retrieves back the (now new, transformed) coordinates into
c     "cora, code" variables.
c


      do i=1,nobs
      cora(i)=ranew(i)
      code(i)=denew(i)
      enddo


c
c     Now comes the final step. For the last time, each individual xy field
c     is reduced to be as best as possible placed in the final global reduction
c     frame system. 
c

      write (*,*)
      write (*,*) 'Final step. Re-reducing single frames with transforme
     ?d positions'
      write (*,*)

      write (46,*)
      write (46,*) 'Final step. Re-reducing single frames with transform 
     ?ed positions'
      write (46,*)


      call reduc (nfxy,jstars,ifxy2,corte,ngrau,ngrau3,ngrau5,cora,code,
     ?idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,
     ?cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,
     ?epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,deest,
     ?zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,kutm,
     ?kutano,kutmes,kutdia,iexps,nx,ny,index)

c
c     Final internal precision of positions from (SA,SD) about the
c    (RA,DEC) mean position of common stars
c

      write (*,*)
      write (46,*)

      xvam=0.d0
      xvas=0.d0

      yvam=0.d0
      yvas=0.d0

      n=0


      do 28 i=1,nobs
      
      if (numcom(i).lt.2) go to 28

      n=n+1

      aux=cora2(i)*15.d0*3600.d0
      auy=code2(i)*3600.d0

      xvam=xvam+aux
      yvam=yvam+auy

      xvas=xvas+aux**2
      yvas=yvas+auy**2

 28   continue

      call desvio (n,xvam,xvas)
      call desvio (n,yvam,yvas)

      write (*,*) 'Final internal precision: (SA,SD)           = ',xvam,
     ?yvam,n
      write (*,*) 'Final sigma of Internal precision (sSA,sSD) = ',xvas,
     ?yvas  

      write (*,*)


      write(46,*) 'Final internal precision: (SA,SD)           = ',xvam,
     ?yvam,n
      write(46,*) 'Final sigma of Internal precision (sSA,sSD) = ',xvas,
     ?yvas  

      write (46,*)


c
c     Writes the new results for the individual xy fields 
c

      write (*,*) 'Writing new single xy files ...'
      write (46,*) 'Writing new single xy files ...'

      write (*,*)
      write (46,*)

      do i=1,nfxy
      
      write (*,*) 'Writing file ',i,nfxy,ifxy2(i)
      write (46,*) 'Writing file ',i,nfxy,ifxy2(i)

      open (2,file=ifxy2(i))

      do j=1,jstars(i)

      k=index(j,i)

      aux=cora2(k)*15.d0*3600.d0
      auy=code2(k)*3600.d0

      if (numcom(k).eq.1) then
      aux=dnove
      auy=dnove
      endif


      write (2,30) x(j,i),y(j,i),seng(j,i),altu(j,i),
     ?fgcc(j,i),fumag(j,i),fumag2(j,i),xmgu(j,i),cudmg(j,i),cudmg2(j,i),
     ?xmgj(j,i),xmgh(j,i),xmgk(j,i),res2mg(j,i),resmg2(j,i),ermgj(j,i),
     ?ermgh(j,i),ermgk(j,i),pma(j,i),pmd(j,i),epma(j,i),epmd(j,i),
     ?ex(j,i),ey(j,i),erau(j,i),edeu(j,i),alfsiu(j,i),delsiu(j,i),
     ?nstaru(j,i),nfinau(j,i),alsiu(j,i),desiu(j,i),ktirau(j,i),
     ?raest(j,i),deest(j,i),kuth(j,i),kutm(j,i),zut(j,i),kutano(j,i),
     ?kutmes(j,i),kutdia(j,i),daju(j,i),iexps(j,i),ichfil(j,i),
     ?infits(j,i),iobalv(j,i),nx(j,i),ny(j,i),numcom(k),aux,auy

 30   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),
     ?13(1x,f6.3),4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,
     ?2(1x,f13.9),1x,i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,
     ?i4,2x,a20,2x,a50,1x,a20,3(1x,i5),2(1x,f7.3))

      enddo

      close (2)

      enddo

c
c     Writes the catalog of all positions in the output file
c

      write (*,*)
      write (*,*) 'Writing catalog from global reduction: ',catalog
      write (*,*) 
      write (*,*) 'Please wait ...'
      write (*,*) 

      write (46,*)
      write (46,*) 'Writing catalog from global reduction: ',catalog
      write (46,*) 
      write (46,*) 'Please wait ...'
      write (46,*) 

      open(2,file=catalog)

c

      mfits=tarall


      i=1
      j=1

      do k=1,nobs


      aux=cora2(k)*15.d0*3600.d0
      auy=code2(k)*3600.d0

      if (numcom(k).eq.1) then
      aux=dnove
      auy=dnove
      endif

      nfin=finac(k)
      ktir=tirac(k)

      write(2,30) x(j,i),y(j,i),cseng(k),altu(j,i),fgcc(j,i),fumag(j,i),
     ?fumag2(j,i),cxmgu(k),codmg(k),codmg2(k),cxmgj(k),cxmgh(k),
     ?cxmgk(k),res2mg(j,i),resmg2(j,i),ermgj(j,i),ermgh(j,i),ermgk(j,i),
     ?copma(k),copmd(k),epma(j,i),epmd(j,i),coex(k),coey(k),cerau(k),
     ?cedeu(k),alfsic(k),delsic(k),nstaru(j,i),nfin,alsiuc(k),
     ?desiuc(k),ktir,cora(k),code(k),kuth(j,i),kutm(j,i),zut(j,i),
     ?kutano(j,i),kutmes(j,i),kutdia(j,i),codj(k),iexps(j,i),
     ?ichfil(j,i),mfits,iobalv(j,i),nx(j,i),ny(j,i),numcom(k),aux,auy


      enddo

      close (2)

c
c     Target information from global reduction catalog of all positions
c

      write (*,*) 
      write (*,*) 'Targets: results from catalog of all positions'
      write (*,*) 

      write (46,*) 
      write (46,*) 'Targets: results from catalog of all positions'
      write (46,*) 


      modo=1

      if1=1
      if2=2
      if3=7

      open (if1,file=catalog)
      open (if2,file=tarall)
      open (if3,file=target)

      call estat (modo,tcom,if1,if2,if3)

      close (if1)
      close (if2)

      rewind (if3)


c
c     Targer information from single xy fields after global reduction
c

      write (*,*) 
      write (*,*) 'Targets: results from single xy fields after global r
     ?eduction'
      write (*,*) 

      write (46,*) 
      write (46,*) 'Targets: results from single xy fields after global r
     ?eduction'
      write (46,*) 

      modo=2

      if1=1
      if2=2
      if3=7

      open (if2,file=tarsgl)

      do i=1,nfxy

      open (if1,file=ifxy2(i))

      call estat (modo,tcom,if1,if2,if3)

      close (if1)
      rewind (if3)

      enddo

      close (if2)
      close (if3)


c
c     That's all, folks!
c




      write (*,*) 
      write (*,*) 
      write (*,*) 'Execution terminated succesfully.'
      write (*,*) 
      write (*,*) 
      write (*,*) 

      write (46,*) 
      write (46,*) 
      write (46,*) 'Execution terminated succesfully.'
      write (46,*) 
      write (46,*) 
      write (46,*) 


      close (46)

    
      return
      end






c
c     
c     Subroutine indexx
c
c     Index of common/non-common stars from all xy fields
c
c
c     Last update: M. Assafin : 19/Apr/2010
c
c
c





      subroutine indexx (nfxy,jstars,ifxy1,ecom,cora,code,numcom,id,
     ?talf,tdel,ntarg,tcom,itakey,nobs,idnest,idfra,iall,x,y,seng,altu,
     ?fgcc,fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,
     ?ermgj,ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?alsiu,desiu,raest,deest,zut,daju,ichfil,infits,iobalv,nstaru,
     ?nfinau,ktirau,kuth,kutm,kutano,kutmes,kutdia,iexps,nx,ny,index)


      IMPLICIT REAL *8 (A-H,O-Z)


      dimension cora(iall),code(iall),numcom(iall),id(iall),talf(iall),
     ?tdel(iall),tadj(iall)

      dimension x(idnest,idfra),y(idnest,idfra),seng(idnest,idfra),
     ?altu(idnest,idfra),fgcc(idnest,idfra),fumag(idnest,idfra),
     ?fumag2(idnest,idfra),xmgu(idnest,idfra),cudmg(idnest,idfra),
     ?cudmg2(idnest,idfra),xmgj(idnest,idfra),xmgh(idnest,idfra),
     ?xmgk(idnest,idfra),res2mg(idnest,idfra),resmg2(idnest,idfra),
     ?ermgj(idnest,idfra),ermgh(idnest,idfra),ermgk(idnest,idfra),
     ?pma(idnest,idfra),pmd(idnest,idfra),epma(idnest,idfra),
     ?epmd(idnest,idfra),ex(idnest,idfra),ey(idnest,idfra),
     ?erau(idnest,idfra),edeu(idnest,idfra),alfsiu(idnest,idfra),
     ?delsiu(idnest,idfra),nstaru(idnest,idfra),nfinau(idnest,idfra),
     ?alsiu(idnest,idfra),desiu(idnest,idfra),ktirau(idnest,idfra),
     ?raest(idnest,idfra),deest(idnest,idfra),kuth(idnest,idfra),
     ?kutm(idnest,idfra),zut(idnest,idfra),kutano(idnest,idfra),
     ?kutmes(idnest,idfra),kutdia(idnest,idfra),daju(idnest,idfra),
     ?iexps(idnest,idfra),nx(idnest,idfra),ny(idnest,idfra),
     ?index(idnest,idfra)

      dimension jstars(idfra)

      character*100 ifxy1(idfra)
      character*100 infits(idnest,idfra)
      character*20 iobalv(idnest,idfra),ichfil(idnest,idfra)

c     common /a1/x,y,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu

c     common /a2/alsiu,desiu,raest,deest,zut,daju,ichfil,infits,iobalv,
c    ?nstaru,nfinau,ktirau,kuth,kutm,kutano,kutmes,kutdia,iexps,nx,ny,
c    ?index



c
c     Auxiliary data
c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI


c
c     Indexing common/non-commom stars
c

     

      ucmag=20.d0
   
      ebox=ecom**2
      tbox=tcom**2

      ertexp=1.d0/86400.d0

      nobs=0
      nest=0


c
      write (*,*) 
      write (46,*) 

      do 20 ii=1,nfxy

      write (*,*) 'Indexing stars of field ',ii,nfxy,ifxy1(ii)
      write (46,*) 'Indexing stars of field ',ii,nfxy,ifxy1(ii)

      do 10 k=1,jstars(ii)

c
c     Excludes target from averaging
c


      if (itakey.eq.1) then

      do 2 ll=1,ntarg
      if (tadj(ll).gt.0.d0) then
      if (dabs(tadj(ll)-daju(k,ii)).gt.ertexp) go to 2
      endif

      dx=(talf(ll)-raest(k,ii))*dcos(deest(k,ii)*grarad)*15.d0*3600.d0
      dy=(tdel(ll)-deest(k,ii))*3600.d0
      dis=dx**2+dy**2
      if (dis.lt.tbox) go to 5

 2    continue

      endif

c
c     Finds common objects within all frames
c

      do i=1,nobs

      raco=cora(i)/numcom(i)
      deco=code(i)/numcom(i)

      dx=(raco-raest(k,ii))*dcos(deest(k,ii)*grarad)*15.d0*3600.d0
      dy=(deco-deest(k,ii))*3600.d0
      dis=dx**2+dy**2

c

      if (dis.lt.ebox) then
      cora(i)=cora(i)+raest(k,ii)
      code(i)=code(i)+deest(k,ii)
      numcom(i)=numcom(i)+1
      index(k,ii)=i


      go to 10

      endif

      enddo

c

 5    continue


      nest=nest+1

      n=nobs+nest

      index(k,ii)=n

c
c     Indexing UCAC2 stars
c

      if (xmgu(k,ii).lt.ucmag) id(n)=-n

c

      cora(n)=raest(k,ii)
      code(n)=deest(k,ii)
      numcom(n)=numcom(n)+1


 10   continue

c

      nobs=nobs+nest

      nest=0

 20   continue


c

      return
      end




c
c     
c     Subroutine update
c
c     Update star positions and other information from all xy files reduced
c     with PRAIA.
c 
c     If star is present in other fields, some values are averaged:
c
c     ra,de,mags, position errors, Julian Date
c
c
c     Last update: M. Assafin : 19/Apr/2010
c
c
c

      subroutine update (jstars,nfxy,ifxy,nobs,cora,code,cora2,code2,
     ?codj,cseng,cxmgu,codmg,codmg2,cxmgj,cxmgh,cxmgk,copma,copmd,coex,
     ?coey,cerau,cedeu,alfsic,delsic,finac,alsiuc,desiuc,tirac,numcom,
     ?idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,cudmg,
     ?cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,pmd,
     ?epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,deest,
     ?zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,kutm,
     ?kutano,kutmes,kutdia,iexps,nx,ny,index)


      IMPLICIT REAL *8 (A-H,O-Z)

      dimension x(idnest,idfra),y(idnest,idfra),seng(idnest,idfra),
     ?altu(idnest,idfra),fgcc(idnest,idfra),fumag(idnest,idfra),
     ?fumag2(idnest,idfra),xmgu(idnest,idfra),cudmg(idnest,idfra),
     ?cudmg2(idnest,idfra),xmgj(idnest,idfra),xmgh(idnest,idfra),
     ?xmgk(idnest,idfra),res2mg(idnest,idfra),resmg2(idnest,idfra),
     ?ermgj(idnest,idfra),ermgh(idnest,idfra),ermgk(idnest,idfra),
     ?pma(idnest,idfra),pmd(idnest,idfra),epma(idnest,idfra),
     ?epmd(idnest,idfra),ex(idnest,idfra),ey(idnest,idfra),
     ?erau(idnest,idfra),edeu(idnest,idfra),alfsiu(idnest,idfra),
     ?delsiu(idnest,idfra),nstaru(idnest,idfra),nfinau(idnest,idfra),
     ?alsiu(idnest,idfra),desiu(idnest,idfra),ktirau(idnest,idfra),
     ?raest(idnest,idfra),deest(idnest,idfra),kuth(idnest,idfra),
     ?kutm(idnest,idfra),zut(idnest,idfra),kutano(idnest,idfra),
     ?kutmes(idnest,idfra),kutdia(idnest,idfra),daju(idnest,idfra),
     ?iexps(idnest,idfra),nx(idnest,idfra),ny(idnest,idfra),
     ?index(idnest,idfra)

      dimension cora(iall),code(iall),codj(iall),cseng(iall),
     ?codmg(iall),codmg2(iall),cxmgj(iall),cxmgh(iall),cxmgk(iall),
     ?copma(iall),copmd(iall),cxmgu(iall),cora2(iall),code2(iall)

      dimension coex(iall),coey(iall),cerau(iall),cedeu(iall),
     ?alfsic(iall),delsic(iall),finac(iall),alsiuc(iall),desiuc(iall),
     ?tirac(iall),numcom(iall)

      dimension jstars(idfra)

      character*100 ifxy(idfra)
      character*100 infits(idnest,idfra)
      character*20 iobalv(idnest,idfra),ichfil(idnest,idfra)

c     common /a1/x,y,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu


c     common /a2/alsiu,desiu,raest,deest,zut,daju,ichfil,infits,iobalv,
c    ?nstaru,nfinau,ktirau,kuth,kutm,kutano,kutmes,kutdia,iexps,nx,ny,
c    ?index



c
c     Auxiliary data
c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c
c     Emptying vectors
c

      do i=1,nobs
      cora(i)=0.d0
      code(i)=0.d0
      cora2(i)=0.d0
      code2(i)=0.d0
      codj(i)=0.d0
      cseng(i)=0.d0
      codmg(i)=0.d0
      codmg2(i)=0.d0
      cxmgu(i)=99.999d0
      cxmgj(i)=99.999d0
      cxmgh(i)=99.999d0
      cxmgk(i)=99.999d0
      copma(i)=99.999d0
      copmd(i)=99.999d0
      coex(i)=0.d0
      coey(i)=0.d0
      cerau(i)=0.d0
      cedeu(i)=0.d0
      alfsic(i)=0.d0
      delsic(i)=0.d0
      finac(i)=0.d0
      alsiuc(i)=0.d0
      desiuc(i)=0.d0
      tirac(i)=0.d0
      enddo

c
c     Averaging values for catalog
c

      write (*,*) 
      write (46,*) 

      do k=1,nfxy

      write (*,*) 'Updating data from ',k,nfxy,ifxy(k)
      write (46,*) 'Updating data from ',k,nfxy,ifxy(k)

      do j=1,jstars(k)

      i=index(j,k)

      cora(i)=cora(i)+raest(j,k)
      code(i)=code(i)+deest(j,k)
      cora2(i)=cora2(i)+raest(j,k)**2
      code2(i)=code2(i)+deest(j,k)**2
      codj(i)=codj(i)+daju(j,k)
      cseng(i)=cseng(i)+seng(j,k)
      codmg(i)=codmg(i)+cudmg(j,k)
      codmg2(i)=codmg2(i)+cudmg2(j,k)

      if (xmgu(j,k).lt.40.d0) cxmgu(i)=xmgu(j,k)
      if (xmgj(j,k).lt.40.d0) cxmgj(i)=xmgj(j,k)
      if (xmgh(j,k).lt.40.d0) cxmgh(i)=xmgh(j,k)
      if (xmgk(j,k).lt.40.d0) cxmgk(i)=xmgk(j,k)
      if (dabs(pma(j,k)).lt.40.d0) copma(i)=pma(j,k)
      if (dabs(pmd(j,k)).lt.40.d0) copmd(i)=pmd(j,k)

      coex(i)=coex(i)+ex(j,k)
      coey(i)=coey(i)+ey(j,k)
      cerau(i)=cerau(i)+erau(j,k)
      cedeu(i)=cedeu(i)+edeu(j,k)
      alfsic(i)=alfsic(i)+alfsiu(j,k)
      delsic(i)=delsic(i)+delsiu(j,k)
      finac(i)=finac(i)+nfinau(j,k)
      alsiuc(i)=alsiuc(i)+alsiu(j,k)
      desiuc(i)=desiuc(i)+desiu(j,k)
      tirac(i)=tirac(i)+ktirau(j,k)

      enddo
      enddo

c

      do 30 ii=1,nobs

      call desvio (numcom(ii),cora(ii),cora2(ii))
      call desvio (numcom(ii),code(ii),code2(ii))

      if (numcom(ii).lt.2) go to 30

      codj(ii)=codj(ii)/numcom(ii)
      cseng(ii)=cseng(ii)/numcom(ii)
      codmg(ii)=codmg(ii)/numcom(ii)
      codmg2(ii)=codmg2(ii)/numcom(ii)
      coex(ii)=coex(ii)/numcom(ii)
      coey(ii)=coey(ii)/numcom(ii)
      cerau(ii)=cerau(ii)/numcom(ii)
      cedeu(ii)=cedeu(ii)/numcom(ii)
      alfsic(ii)=alfsic(ii)/numcom(ii)
      delsic(ii)=delsic(ii)/numcom(ii)
      finac(ii)=finac(ii)/numcom(ii)
      alsiuc(ii)=alsiuc(ii)/numcom(ii)
      desiuc(ii)=desiuc(ii)/numcom(ii)
      tirac(ii)=tirac(ii)/numcom(ii)

 30   continue

c

      return
      end



c
c
c     
c     Subroutine reduc
c
c     Re-reduces (RA,DEC) using all stars as reference stars
c 
c
c     Last update: M. Assafin : 19/Apr/2010
c

      subroutine reduc (nfxy,jstars,ifxy,corte,ngrau,ngrau3,ngrau5,
     ?cora,code,idnest,idfra,iall,x,y,seng,altu,fgcc,fumag,fumag2,xmgu,
     ?cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,ermgk,pma,
     ?pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,alsiu,desiu,raest,
     ?deest,zut,daju,ichfil,infits,iobalv,nstaru,nfinau,ktirau,kuth,
     ?kutm,kutano,kutmes,kutdia,iexps,nx,ny,index)

      IMPLICIT REAL *8 (A-H,O-Z)


      dimension x(idnest,idfra),y(idnest,idfra),seng(idnest,idfra),
     ?altu(idnest,idfra),fgcc(idnest,idfra),fumag(idnest,idfra),
     ?fumag2(idnest,idfra),xmgu(idnest,idfra),cudmg(idnest,idfra),
     ?cudmg2(idnest,idfra),xmgj(idnest,idfra),xmgh(idnest,idfra),
     ?xmgk(idnest,idfra),res2mg(idnest,idfra),resmg2(idnest,idfra),
     ?ermgj(idnest,idfra),ermgh(idnest,idfra),ermgk(idnest,idfra),
     ?pma(idnest,idfra),pmd(idnest,idfra),epma(idnest,idfra),
     ?epmd(idnest,idfra),ex(idnest,idfra),ey(idnest,idfra),
     ?erau(idnest,idfra),edeu(idnest,idfra),alfsiu(idnest,idfra),
     ?delsiu(idnest,idfra),nstaru(idnest,idfra),nfinau(idnest,idfra),
     ?alsiu(idnest,idfra),desiu(idnest,idfra),ktirau(idnest,idfra),
     ?raest(idnest,idfra),deest(idnest,idfra),kuth(idnest,idfra),
     ?kutm(idnest,idfra),zut(idnest,idfra),kutano(idnest,idfra),
     ?kutmes(idnest,idfra),kutdia(idnest,idfra),daju(idnest,idfra),
     ?iexps(idnest,idfra),nx(idnest,idfra),ny(idnest,idfra),
     ?index(idnest,idfra)

      dimension jstars(idfra)

      dimension cora(iall),code(iall)
      dimension id(idnest),racat(idnest),decat(idnest),xob(idnest),
     ?yob(idnest)

      dimension ra(idnest),de(idnest),era(idnest),ede(idnest),coefx(21),
     ?coefy(21),ecoefx(21),ecoefy(21),alfres(idnest),delres(idnest),
     ?itira(idnest)

      character*100 ifxy(idfra)
      character*100 infits(idnest,idfra)
      character*20 iobalv(idnest,idfra),ichfil(idnest,idfra)


c     common /a1/x,y,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu

c     common /a2/alsiu,desiu,raest,deest,zut,daju,ichfil,infits,iobalv,
c    ?nstaru,nfinau,ktirau,kuth,kutm,kutano,kutmes,kutdia,iexps,nx,ny,
c    ?index



      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0


c
C     Auxiliary data
C
c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c

      ireflex=1

      idiobs=idnest

      ncof=21

c
c     Loop of xy files
c

      write (*,*) '   Sa      Sd    Stars   field total filename'
      write (46,*) '   Sa      Sd    Stars   field total filename'

      do 100 kkk=1,nfxy


c
c     Initializing vectors
c


      do i=1,idiobs

      id(i)=0
      racat(i)=0.d0
      decat(i)=0.d0
      xob(i)=0.d0
      yob(i)=0.d0
      ra(i)=0.d0
      de(i)=0.d0
      era(i)=0.d0
      ede(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0
      itira(i)=0

      enddo

      do i=1,ncof

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0

      enddo



c
c     Retrieving data from former field reduction
c


      rac=0.d0
      dec=0.d0  

      do i=1,jstars(kkk)

      j=index(i,kkk)

      id(i)=i

      xob(i)=x(i,kkk)
      yob(i)=y(i,kkk)

      rac=rac+raest(i,kkk)
      dec=dec+deest(i,kkk)

      racat(i)=cora(j)*15.d0
      decat(i)=code(j)

      enddo

      nest=jstars(kkk)


c
c     Guess central (RA,DEC) of field
c
c     rac = RA  center in degrees
c     dec = DEC center in degrees 
c     


      rac=(rac/nest)*15.d0
      dec=dec/nest


c
c     RA,DEC Reduction
c


      call posred (ireflex,rac,dec,id,racat,decat,nest,xob,yob,corte,
     ?ngrau,ngrau3,ngrau5,nstart,nfinal,ra,de,era,ede,alfsig,delsig,
     ?alfres,delres,coefx,coefy,ecoefx,ecoefy,itira,avam,dvam,idnest)


c
c     Updates individual results of new reduction
c


      write (*,30) alfsig,delsig,nfinal,kkk,nfxy,ifxy(kkk)
      write (46,30) alfsig,delsig,nfinal,kkk,nfxy,ifxy(kkk)
 30   format(2(1x,f7.3),1x,i5,2x,2(1x,i5),1x,a100)


      do i=1,jstars(kkk)

      raest(i,kkk)=ra(i)/15.d0
      deest(i,kkk)=de(i)
      erau(i,kkk)=era(i)
      edeu(i,kkk)=ede(i)
      alfsiu(i,kkk)=alfsig
      delsiu(i,kkk)=delsig
      nstaru(i,kkk)=nstart
      nfinau(i,kkk)=nfinal
      alsiu(i,kkk)=alfres(i)
      desiu(i,kkk)=delres(i)
      ktirau(i,kkk)=itira(i)
      infits(i,kkk)=ifxy(kkk)

      enddo

c

 100  continue

c

      return
      end



c 
c
c     Subrotina posred
c
c
c     Reducao alfa e delta das posicoes (x,y) medidas em relacao a um
c     catalogo de referencia
c
c
c     Atualmente aplicada para reducao com o catalogo 2MASS e com o 
c     catalogo UCAC2
c
c     Assume-se que as poosicoes (RA,DEC) de catalogo estao na epoca
c     da observacao, no sistema ICRS (J2000)
c
c     Nessa versao nao ha peso.
c     Nessa versao nao ha termos de magnitude ou cor no polinomio de ajuste.
c
c     Ajustes possiveis: 4ctes, 1o, 2o, 3o graus completos, 2o+3dr,
c     2o+3dr+5dr, 3o+5dr.
c
c
c     - (RA,DEC) entram em graus, saem em graus
c     - (x,y)s entram em pixels
c     - erros e sigmas saem em segundos de arco (")
c     - coeficientes (e seus erros) saem em radianos por pixel,
c       rad por pixel**2 etc
c
c
c
c
c
c     Last update: M. Assafin  19/Apr/2010
c
c
c



      subroutine posred (ireflex,rac,dec,id,racat,decat,nest,xob,yob,
     ?corte,ngrau,ngrau3,ngrau5,nstart,nfinal,ra,de,era,ede,alfsig,
     ?delsig,alfres,delres,coefx,coefy,ecoefx,ecoefy,itira,avam,dvam,
     ?idnest)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension id(idnest),racat(idnest),decat(idnest),xob(idnest),
     ?yob(idnest),xp(idnest),yp(idnest),xest(idnest),yest(idnest)

      dimension ra(idnest),de(idnest),era(idnest),ede(idnest),coefx(21),
     ?coefy(21),ecoefx(21),ecoefy(21),alfres(idnest),delres(idnest),
     ?itira(idnest),xsao(21),ysao(21),xrray(21,21),yrray(21,21),
     ?array(21,21)



      COMMON /A7/ARRAY
      COMMON/A14/IERRO


      HMSGMS(I,J,A)=I+J/60.D0+A/3600.D0
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Dados auxiliares
C

      idiobs=idnest
      ncof=21

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
      IERRO=0
      DNOVE=9.999D0
      d99=99.99d0

      EQUIN=2000.d0
      IZERO=0
      ZERO=0.D0
      pma=0.d0
      pmd=0.d0

c
c     Zerando vetores
c

      do i=1,idiobs

      xest(i)=0.d0
      yest(i)=0.d0
      xp(i)=0.d0
      yp(i)=0.d0
      ra(i)=0.d0
      de(i)=0.d0
      era(i)=0.d0
      ede(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0

      enddo


      do i=1,21

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0
      xsao(i)=0.d0
      ysao(i)=0.d0

      do j=1,21
      xrray(j,i)=0.d0
      yrray(j,i)=0.d0
      array(j,i)=0.d0
      enddo

      enddo


c
c     Recolhendo estrelas medidas comuns ao catalogo de referencia
c




      nstart=0
      grac=grarad*rac
      gdec=grarad*dec

      do 10 i=1,nest

      if (id(i).eq.0) go to 10

      nstart=nstart+1

      xest(nstart)=ireflex*xob(i)
      yest(nstart)=yob(i)


c
c     Projecao do catalogo de referencia no plano tangente
c


      bra=grarad*racat(id(i))
      bde=grarad*decat(id(i))
      d=DEXY(bra,bde,grac,gdec)

      xp(nstart)=xpad(bra,bde,grac)/d
      yp(nstart)=ypad(bra,bde,grac,gdec)/d


 10   continue


c
c     Ajuste do modelo polinomial entre (x,y) e (X,Y)
c



      call solucao (ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,yp,ntira,
     ?coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,corte,
     ?xrray,yrray,idnest)


      if (ierro.eq.1) then

      write (*,*) 'L.S. solution crashed.'
      write (46,*) 'L.S. solution crashed.'

c
c     Zerando vetores no erro
c

      do i=1,idiobs

      xest(i)=0.d0
      yest(i)=0.d0
      xp(i)=0.d0
      yp(i)=0.d0
      ra(i)=0.d0
      de(i)=0.d0
      era(i)=0.d0
      ede(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0
      itira(i)=0

      enddo


      do i=1,ncof

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0
      xsao(i)=0.d0
      ysao(i)=0.d0
      do j=1,ncof
      xrray(j,i)=0.d0
      yrray(j,i)=0.d0
      array(j,i)=0.d0
      enddo

      enddo

      return
      endif



c

      nfinal=nstart-ntira

c
c     Determinacao do alfa e delta observado de cada estrela do campo
c
c     (RA,DEC) guardados em graus

c
c     Calculo do erro padrao em alfa e delta para cada estrela de campo
c

      XVAM=0.D0
      YVAM=0.D0
      XVAS=0.D0
      YVAS=0.D0


C
C     NGRAU=0 --> 4 ctes
C
      nterms=1
      do i=1,ngrau
      nterms=nterms+i+1
      enddo
C
      igrau=ngrau+1
      if (ngrau.eq.0) igrau=2
C
      if (ngrau.eq.0) then
      nterms=4
c     ngrau3=0
c     ngrau5=0
      endif
C
      iterms=nterms
      if (ngrau3.eq.3) iterms=iterms+1
      if (ngrau5.eq.5) iterms=iterms+1
C



      do 130 i=1,nest

      x=xob(i)*ireflex
      y=yob(i)

      ICONT=0
      POLX=0.D0
      POLY=0.D0
      DO 20 N=1,IGRAU
      DO 20 LL=1,N
      ICONT=ICONT+1
      K=N-LL
      POLX=POLX+COEFX(ICONT)*(X**K)*(Y**(LL-1))
   20 POLY=POLY+COEFY(ICONT)*(X**K)*(Y**(LL-1))
C
      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)
      ENDIF

      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF

      RA(I)=POLX
      DE(I)=POLY

C
C     Calcula o erro padrao em alfa para cada estrela
C
      ICONT=0
      if (ngrau.ne.0) then
      DO  30 N=1,IGRAU
      DO  30 LL=1,N
      K=N-LL
      ICONT=ICONT+1
   30 XSAO(ICONT)=(X**K)*(Y**(LL-1))
      IF (NGRAU3.EQ.3) THEN
       ICONT=ICONT+1
       XSAO(ICONT)=X*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
       ICONT=ICONT+1
       XSAO(ICONT)=X*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
c
      else
      xsao(1)=x
      xsao(2)=y
      xsao(3)=1.d0
      xsao(4)=0.d0
      endif
c
      DO  40 K=1,ITERMS
   40 YSAO(K)=0.D0
C
      DO  50 LL=1,ITERMS
      DO  50 K=1,ITERMS
   50 YSAO(LL)=YSAO(LL)+XRRAY(LL,K)*XSAO(K)

      DO  60 K=1,ITERMS
   60 era(I)=era(I)+XSAO(K)*YSAO(K)

      era(I)=alfsig*DSQRT(era(I))
      XVAM=XVAM+era(I)
      XVAS=XVAS+era(I)**2


C
C     Calcula o erro padrao em delta para cada estrela
C

      ICONT=0
      if (ngrau.ne.0) then
      DO  70 N=1,IGRAU
      DO  70 LL=1,N
      K=N-LL
      ICONT=ICONT+1
   70 YSAO(ICONT)=(X**K)*(Y**(LL-1))
      IF (NGRAU3.EQ.3) THEN
       ICONT=ICONT+1
       YSAO(ICONT)=Y*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
       ICONT=ICONT+1
       YSAO(ICONT)=Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
c
      else
      ysao(1)=y
      ysao(2)=-x
      ysao(3)=0.d0
      ysao(4)=1.d0
      endif
c
      DO  80 K=1,ITERMS
   80 XSAO(K)=0.D0
C
      DO  90 LL=1,ITERMS
      DO  90 K=1,ITERMS
   90 XSAO(LL)=XSAO(LL)+YRRAY(LL,K)*YSAO(K)

      DO 100 K=1,ITERMS
  100 ede(I)=ede(I)+YSAO(K)*XSAO(K)

      ede(I)=delsig*DSQRT(ede(I))
      YVAM=YVAM+ede(I)
      YVAS=YVAS+ede(I)**2


C
C     Calcula erro dos coeficientes para a solucao X
C

c     IF (NGRAU.NE.0) THEN
c     DO 110 K=1,ITERMS
c110  ECOEFX(K)=ALFSIG*DSQRT(XRRAY(K,K))
c     ELSE
c     ECOEFX(1)=ALFSIG*DSQRT(XRRAY(3,3))
c     ECOEFX(2)=ALFSIG*DSQRT(XRRAY(1,1))
c     ECOEFX(3)=ALFSIG*DSQRT(XRRAY(2,2))
c     ENDIF



C
C     Calcula erro dos coeficientes para a solucao Y
C

c     IF (NGRAU.NE.0) THEN
c     DO 120 K=1,ITERMS
c120  ECOEFY(K)=DELSIG*DSQRT(YRRAY(K,K))
c     ELSE
c     ECOEFY(1)=DELSIG*DSQRT(YRRAY(4,4))
c     ECOEFY(2)=DELSIG*DSQRT(YRRAY(2,2))
c     ECOEFY(3)=DELSIG*DSQRT(YRRAY(1,1))
c     ENDIF

c     do iii=1,nest
c     write (*,*) 'itira = ',itira(iii)
c     enddo
c     stop

C
C
 130  CONTINUE
C
C
      if (ngrau.eq.0) then
      iterms=3
      endif

C
C     Media dos erros padrao de todas as estrelas medidas
C

      EXMED=XVAM/NEST
      EYMED=YVAM/NEST
      XVAS=DSQRT((XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2)/(NEST-1.D0))
      YVAS=DSQRT((YVAS-2.D0*EYMED*YVAM+NEST*EYMED**2)/(NEST-1.D0))
      AVAM=EXMED
      DVAM=EYMED

c
c     RA e DEC em graus
c


      j=0

      do 140 i=1,nest

      x=ra(i)
      y=de(i)

      ra(i)=alff(x,y,grac,gdec)
      de(i)=deltt(ra(i),y,grac,gdec)

      ra(i)=ra(i)*radgra
      de(i)=de(i)*radgra

c     if (id(i).eq.0) go to 140

c     j=j+1

c     write(*,*) 'alfres delres itira = ',i,alfres(i),delres(i),itira(i) 



 140  continue


c
c     Debug
c


c     perc=100.d0*ntira/nstart
c     write (*,141) alfsig,delsig,nstart,nfinal,perc,avam,dvam
c141  format(1x,'alfsig delsig NI NF = ',2(1x,f6.3),2(1x,i4),1x,f6.2,
c    ?'%',2(1x,f6.3))





      return
      end



C
C     Subrotina solucao
C
C
C     Ajuste polinomial bivariado P=P(x,y) ate grau 3 e distorcao
C     radial a grau 5
C
c
c     Last update: 19/Apr/2010  M. Assafin
c


      subroutine solucao (ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,yp,
     ?ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,
     ?corte,xrray,yrray,idnest)


      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(21),COEFY(21),ALPHAX(21,21),ALPHAY(21,21),
     ?ARRAY(21,21),BETAX(21),BETAY(21),TERMX(21),TERMY(21),
     ?ITIRA(idnest)

      DIMENSION XEST(idnest),YEST(idnest),XP(idnest),YP(idnest),
     ?XRRAY(21,21),YRRAY(21,21),ALFRES(idnest),DELRES(idnest)

      COMMON /A7/ARRAY
      COMMON/A14/IERRO


      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
      DET  =1.D0
      IERRO=0

      ifora=0
c

      alfsig=0.d0
      delsig=0.d0

c

      sigma=grarad*corte/3600d0

C
      NCOMUM=nstart
C
c
c     Computa no. de termos do polinomio
c     NGRAU=0 --> 4 constantes
c

      NTERMS=1
      DO I=1,NGRAU
      NTERMS=NTERMS+I+1
      ENDDO

C
      IGRAU=NGRAU+1
      IF (NGRAU.EQ.0) IGRAU=2
      NTIRA=0
C


 5    LCONT=0

C
      RESMX=0.D0
      RES2X=0.D0
      RESMY=0.D0
      RES2Y=0.D0
C
      IF (NGRAU.EQ.0) THEN
      NTERMS=4
      NGRAU3=0
      NGRAU5=0
      ENDIF
C
      ITERMS=NTERMS
      IF (NGRAU3.EQ.3) ITERMS=ITERMS+1
      IF (NGRAU5.EQ.5) ITERMS=ITERMS+1

C
C     Checa No. de estrelas versus no. de coeficientes a ajustar
C

      IF (NGRAU.EQ.0) THEN
      iequa=2.d0*ncomum
      mtira=2.d0*ntira
      if ((iequa-mtira).LT.ITERMS) then
      IERRO=1
      RETURN
      endif
      else
      if ((ncomum-ntira).LT.ITERMS) then
      IERRO=1
      RETURN
      endif
      ENDIF

c
      DO 9 I=1,ITERMS
      BETAX(I) =0.D0
      BETAY(I) =0.D0
      COEFX(I) =0.D0
      COEFY(I) =0.D0
      TERMX(I) =0.D0
      TERMY(I) =0.D0
      DO 9 J=I,ITERMS
      ALPHAX(I,J)=0.D0
 9    ALPHAY(I,J)=0.D0


C
C     Montando equacoes de condicao para ajuste com 4 ctes
C

      IF (NGRAU.EQ.0) THEN

      DO 10 I=1,NCOMUM
      IF (ITIRA(I).NE.0) GO TO 10
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computa  termos dos coeficientes para AtB
C

      TERMX(1)=TERMX(1)+X**2+Y**2
      TERMX(2)=TERMX(2)+X
      TERMX(3)=TERMX(3)+Y
      TERMX(4)=TERMX(4)+1.D0
      BETAX(1)=BETAX(1)+X*XG+Y*YG
      BETAX(2)=BETAX(2)+Y*XG-X*YG
      BETAX(3)=BETAX(3)+XG
      BETAX(4)=BETAX(4)+YG
   10 CONTINUE

C
C     Preenchendo  AtA para ajuste com 4 ctes
C

      ALPHAX(1,1)=TERMX(1)
      ALPHAX(1,2)=0.D0
      ALPHAX(1,3)=TERMX(2)
      ALPHAX(1,4)=TERMX(3)
      ALPHAX(2,2)=TERMX(1)
      ALPHAX(2,3)=TERMX(3)
      ALPHAX(2,4)=-TERMX(2)
      ALPHAX(3,3)=TERMX(4)
      ALPHAX(3,4)=0.D0
      ALPHAX(4,4)=TERMX(4)
      DO 15 L=1,ITERMS
      DO 15 K=L,ITERMS
   15 ALPHAX(K,L) =ALPHAX(L,K)
C
      GO TO 75

      ENDIF


C
C     Montando equacoes de condicao para ajustes que NAO o de 4 ctes
C     para polinomio X


      DO 65 I=1,NCOMUM
      IF (ITIRA(I).NE.0) GO TO 65
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computando termos para o polinomio em X
C
      ICONT=0
      DO 40 N=1,IGRAU
      DO 40 L=1,N
      K=N-L
      ICONT=ICONT+1
      TERMX(ICONT)=(X**K)*(Y**(L-1))
   40 TERMY(ICONT)=(X**K)*(Y**(L-1))

      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      TERMX(ICONT)=X*(X**2+Y**2)
      TERMY(ICONT)=Y*(X**2+Y**2)
      ENDIF

      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      TERMX(ICONT)=X*(X**2+Y**2)*(X**2+Y**2)
      TERMY(ICONT)=Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
C
C     Computando AtA and AtB (para "X" e "Y")
C

      DO 60 L=1,ITERMS
      BETAX(L)=BETAX(L)+XG*TERMX(L)
      BETAY(L)=BETAY(L)+YG*TERMY(L)
      DO 60 K=L,ITERMS
      ALPHAX(L,K)=ALPHAX(L,K)+TERMX(K)*TERMX(L)
   60 ALPHAY(L,K)=ALPHAY(L,K)+TERMY(K)*TERMY(L)
C
   65 CONTINUE

C
C     Preenchendo parte triangular inferior da matriz simetrica AtA
c     (para "X" e "Y")
C

      DO 70 L=1,ITERMS
      DO 70 K=L,ITERMS
      ALPHAX(K,L)=ALPHAX(L,K)
   70 ALPHAY(K,L)=ALPHAY(L,K)

C
C     Preenchendo ARRAY=AtA para inversao (elementos normalizados pela
C     diagonal) para "X"
C

 75   CONTINUE
C
      DO 80 L=1,ITERMS
      DO 80 K=1,ITERMS
   80 ARRAY(L,K)=ALPHAX(L,K)/DSQRT(ALPHAX(L,L)*ALPHAX(K,K))
C
C     Invertendo AtA para "X"
C



      CALL MATINV (ITERMS,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando coeficientes do polinomio para "X"
C
      DO 90 L=1,ITERMS
      DO 90 K=1,ITERMS
   90 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHAX(K,K)*ALPHAX(L,L))
C
C     Backup do array para "X"
C
      IF (NGRAU.EQ.0) GO TO 97
      DO 95 L=1,ITERMS
      DO 95 K=1,ITERMS
   95 XRRAY(L,K)=ARRAY(L,K)
C
 97   DO 100 L=1,ITERMS
      DO 100 K=1,ITERMS
  100 COEFX(L)=COEFX(L)+ARRAY(L,K)*BETAX(K)
C
C     Obtem coeficientes e backup do array para modelo de 4 Constantes
C
      IF (NGRAU.EQ.0) THEN
      DO  105 L=1,ITERMS
 105  TERMX(L)=COEFX(L)
      COEFX(1)=TERMX(3)
      COEFX(2)=TERMX(1)
      COEFX(3)=TERMX(2)
      COEFY(1)=TERMX(4)
      COEFY(2)=-TERMX(2)
      COEFY(3)=TERMX(1)
      do 107 l=1,iterms
      do 107 k=1,iterms
      xrray(l,k)=array(l,k)
 107  yrray(l,k)=array(l,k)
      GO TO 133
      ENDIF

C
C     Preenchendo ARRAY=AtA para inversao (elementos normalizados pela
C     diagonal) para "Y"
C

      DO 110 L=1,ITERMS
      DO 110 K=1,ITERMS
  110 ARRAY(L,K)=ALPHAY(L,K)/DSQRT(ALPHAY(L,L)*ALPHAY(K,K))
C
C     Invertendo AtA para "Y"
C

      CALL MATINV (ITERMS,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando coeficientes do polinomio para "Y"
C


      DO 120 L=1,ITERMS
      DO 120 K=1,ITERMS
  120 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHAY(K,K)*ALPHAY(L,L))
C
C     Backup do array para "Y"
C

      DO 125 L=1,ITERMS
      DO 125 K=1,ITERMS
  125 YRRAY(L,K)=ARRAY(L,K)
C
      DO 130 L=1,ITERMS
      DO 130 K=1,ITERMS
  130 COEFY(L)=COEFY(L)+ARRAY(L,K)*BETAY(K)

c
c
c     Computa residuos. Estrelas com residuo em alfa ou delta mais alto
c     sao eliminadas uma a uma, ate que nenhuma possua (O-C) maior
c     que o valor definido pela variavel "corte"
c

 133  CONTINUE
C


      REMAXI=-1.D14

      DO 160 I=1,NCOMUM
      X=XEST(I)
      Y=YEST(I)
      xg=xp(i)
      yg=yp(i)
      ICONT=0
      POLX=0.D0
      POLY=0.D0
      DO 140 N=1,IGRAU
      DO 140 L=1,N
      ICONT=ICONT+1
      K=N-L
      POLX=POLX+COEFX(ICONT)*(X**K)*(Y**(L-1))
  140 POLY=POLY+COEFY(ICONT)*(X**K)*(Y**(L-1))
C
      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF

C
C     Computa (O-C) para R.A. e  Dec
C

      xx=alff(polx,poly,grac,gdec)
      yy=deltt(xx,poly,grac,gdec)

      xxr=alff(xg,yg,grac,gdec)
      yyr=deltt(xxr,yg,grac,gdec)


      AUX=(XX-XXR)*DCOS(YYR)
      AUY=YY-YYR
      ALFRES(I)=AUX*radgra*3600d0
      DELRES(I)=AUY*radgra*3600d0


C
      IF (ITIRA(I).NE.0) GO TO 160
      IF ((DABS(AUX).GT.REMAXI).OR.(DABS(AUY).GT.REMAXI)) THEN
      IFORA =I
      REMAXI=DMAX1(DABS(AUX),DABS(AUY))
      ENDIF

      RESMX=RESMX+ALFRES(I)
      RES2X=RES2X+ALFRES(I)**2
      RESMY=RESMY+DELRES(I)
      RES2Y=RES2Y+DELRES(I)**2

      LCONT=LCONT+1

C
  160 CONTINUE
C
C     Atingido numero minimo de estrelas!
C
      iwar=0
      IF (NGRAU.EQ.0) THEN
      iequa=2.d0*LCONT
      if (iequa.eq.4) iwar=1
      ELSE
      if (LCONT.eq.ITERMS) iwar=1
      ENDIF
C
      IF (IWAR.EQ.1) THEN
      alfsig=0.D0
      delsig=0.D0
      RETURN
      ENDIF


C
C     Computa media e desvio padrao dos (O-C)s
C

      XMED=RESMX/LCONT
      YMED=RESMY/LCONT
C
      IF (NGRAU.EQ.0) THEN
      alfsig=DSQRT((RES2X-2.D0*XMED*RESMX+LCONT*XMED**2)/(LCONT-2.D0))
      delsig=DSQRT((RES2Y-2.D0*YMED*RESMY+LCONT*YMED**2)/(LCONT-2.D0))
      ELSE
      alfsig=DSQRT((RES2X-2.D0*XMED*RESMX+LCONT*XMED**2)/(LCONT-ITERMS))
      delsig=DSQRT((RES2Y-2.D0*YMED*RESMY+LCONT*YMED**2)/(LCONT-ITERMS))
      ENDIF
C
c
c     Atingido o corte !
c

c     write (*,*) 'remaxi ',remaxi*radgra*3600.d0

      if (remaxi.le.sigma) go to 200

c
c     Corte nao atingido, prosseguir com eliminacao de estrelas
c

      if (ifora.ne.0) then

      NTIRA=NTIRA+1
      ITIRA(IFORA)=1
      GO TO 5

      endif

 200  continue

      RETURN
      END


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
c     Mean and standard deviation
c
c
c     entry:
c
c        xvam  = sum of individual values
c        xvas  = sum of square of individual values
c
c     output:
c
c        xvam  = mean
c        xvas  = standard deviation about that mean
c



      subroutine desvio (nest,xvam,xvas)

      IMPLICIT REAL *8 (A-H,O-Z)

c

      dnove=99.999d0

c

      EXMED=XVAM/NEST

c

      if (nest.eq.1) then

      xvas=dnove

      return

      endif

c

      if (nest.eq.2) then

      xvas=dsqrt(dabs(2.d0*xvas-xvam**2))
      XVAM=EXMED

      return

      endif

c

      raiz=XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2

      if (raiz.lt.0.d0) then

      xvas=0.d0

      else

      XVAS=DSQRT(raiz/(NEST-1.D0))

      endif
c

      XVAM=EXMED


      return

      end




c
c
c     subroutine ucac2
c
c
c     Picks up UCAC2 stars with no speedy extraction
c
c
c     Last update:   15/Feb/2013.  M. Assafin
c
c

      subroutine ucac2 (uraiz,ecom,alfmin,alfmax,delmin,delmax,nobs,
     ?cora,code,codj,ucra,ucde,id,oldra,oldde,rac,dec,nok,rmx,rmy,rmx2,
     ?rmy2,iall)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra

      dimension cora(iall),code(iall),ucra(iall),
     ?ucde(iall),codj(iall),id(iall)

      dimension oldra(iall),oldde(iall),pmx(iall),
     ?pmy(iall)

      dimension idat(23)


      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 uraiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/


C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      un=95

      idin50=50
      izon=288

      nok=0

      ndim=iall

      ebox=ecom**2

c

      write (*,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax
      write (46,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax

      write (*,*)
      write (46,*)



      alma=alfmax*15.d0
      almi=alfmin*15.d0

c
c     Reads UCAC2 declination zones of 0.5 degrees
c

      dfaixa=delmin-0.5d0
      decmax=delmax

      nest=1


      do 30 k=1,288

      dfaixa=dfaixa+0.5d0


      if (dfaixa-decmax.gt.0.5d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.5d0+1.0d0

c
c     UCAC2 North declination limit
c
c     For now, just zones up to z288 will do
c

      if (j.gt.izon) go to 30

c
c     Mounts 0.5 degrees zone file name 
c


      ifaixa=''
      ifaixa=uraiz

      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC2 ',a54)



c
c     Reads current zone in the loop
c

      open (95,file=ifaixa,access="direct",form="unformatted",recl=44)


      n=0

 20   n=n+1

c
c     Picks up star data
c    

      call readu2 (un,n,idat,ierra)
      if (ierra.eq.1) go to 25


      oldra(nest)=idat(1)
      oldde(nest)=idat(2)
      pmx(nest)=idat(12)
      pmy(nest)=idat(13)

c
c     Accelerates match
c

      rauc=oldra(nest)/3600d3

      if (rauc.lt.almi) go to 20
      if (rauc.gt.alma) go to 20

      nest=nest+1

      if (nest.gt.ndim) then
      close (95)
      go to 36
      endif

      go to 20


 25   close (95)

c

 30   continue

c

 35   continue

c

 36   continue


      nest=nest-1

      if (nest.ge.ndim) then

      write (*,*)
      write (46,*)

      write (*,*) 'Max. number of UCAC2 stars reached: ',nest,ndim
      write (46,*) 'Max. number of UCAC2 stars reached: ',nest,ndim

      write (*,*)
      write (46,*)

      nest=ndim

      write (*,*) 'UCAC2 stars candidates reduced to: ',nest
      write (46,*) 'UCAC2 stars candidates reduced to: ',nest

      write (*,*)
      write (46,*)


      else

      write (*,*)
      write (46,*)

      write (*,*) 'UCAC2 candidate stars: ',nest
      write (46,*) 'UCAC2 candidate stars: ',nest

      write (*,*)
      write (46,*)

      endif


c
c     Finds if it is a match: measured star x UCAC2 star
c


      write (*,*)
      write (*,*) 'Matching xy field stars vs. UCAC2 stars. Please wait 
     ?...'
      write (*,*)

      write (46,*)
      write (46,*)'Matching xy field stars vs. UCAC2 stars. Please wait 
     ?...'
      write (46,*)

c



      rac=0.d0
      dec=0.d0

      rmx=0.d0
      rmy=0.d0
     
      rmx2=0.d0
      rmy2=0.d0


      do 60 i=1,nest
      do 50 ko=1,nobs

      if (id(ko).ge.0) go to 50

      epoj=2000D0+(codj(ko)-2451545D0)/365.25D0

      rauc=oldra(i)+pmx(i)*(epoj-2000d0)/10d0
      deuc=oldde(i)+pmy(i)*(epoj-2000d0)/10d0

      rauc=rauc/1000d0
      deuc=deuc/1000d0

      rauc=rauc/3600d0
      deuc=deuc/3600d0

      dx=(15.d0*cora(ko)-rauc)*dcos(deuc*grarad)*3600.d0
      dy=(code(ko)-deuc)*3600.d0


      dis=dx**2+dy**2

      if (dis.lt.ebox) then

      ucra(ko)=rauc/15.d0
      ucde(ko)=deuc
      id(ko)=ko

      rac=rac+rauc
      dec=dec+deuc

      rmx=rmx+dx
      rmy=rmy+dy
      rmx2=rmx2+dx**2
      rmy2=rmy2+dy**2

c
c
c     debug alfa delta
c
c     write (*,*) 'oldra, oldde, codj = ',oldra(ko),oldde(ko),codj(ko)
c
c     ra=rauc
c     de=deuc
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,

      nok=nok+1

      go to 60

      endif

c

 50   continue

 60   continue

      rac=rac/nok
      dec=dec/nok

      call desvio (nok,rmx,rmx2)
      call desvio (nok,rmy,rmy2)


      write (*,*)
      write (*,*) 'No. UCAC2 common stars = ',nok
      write (*,*)

      write (46,*)
      write (46,*) 'No. UCAC2 common stars = ',nok
      write (46,*)

c     stop

      return
      end



C***********************************************************************
C
C  read a single record of UCAC2 data = 1 star
C  input:
C    un   = unit number of file  (assumed to be open)
C    recn = record number on that file
C    bf   = .TRUE. if byte flip required
C  output:
C    idat = integer*4 vector of 23 items (see readme2.txt)
C    errflg = 0=ok, 1=erro
C
C
C   Modificada: M. Assafin 12/Dez/2005
C
C



      SUBROUTINE readu2 (un,recn,idat,ierra)

      IMPLICIT NONE

      integer*2 un,ierra
      INTEGER  recn, idat(23)  ! item #24,25 = star ID options

c     INTEGER   ra2000, dc2000, pmx,pmy, id2m, u2id,r11
      INTEGER   ra2000, dc2000, pmx,pmy, id2m
      INTEGER*2 mag, cepx,cepy, j2m,h2m,k2m
      BYTE      sigx,sigy,nobs,epos,ncat,cflg    ! INTEGER*1
     .         ,spmx,spmy, rx,ry, ph,cc          ! signed integer

      ierra = 0                           ! default


      READ (un,REC=recn,ERR=99) ra2000,dc2000
     .  ,mag,sigx,sigy, nobs,epos,ncat,cflg
     .  ,cepx,cepy, pmx,pmy, spmx,spmy, rx,ry
     .  ,id2m, j2m,h2m,k2m, ph,cc

c note: first assign I*1 to idat(I*4), 
c       then add 127 to avoid overflow

      idat ( 1) = ra2000
      idat ( 2) = dc2000
      idat ( 3) = mag
      idat ( 4) = sigx
      idat ( 4) = idat ( 4) + 127
      idat ( 5) = sigy
      idat ( 5) = idat ( 5) + 127
      idat ( 6) = nobs
      idat ( 7) = epos 
      idat ( 7) = idat ( 7) + 127
      idat ( 8) = ncat
      idat ( 9) = cflg
      idat (10) = cepx
      idat (11) = cepy
      idat (12) = pmx
      idat (13) = pmy
      idat (14) = spmx
      idat (14) = idat (14) + 127
      idat (15) = spmy
      idat (15) = idat (15) + 127
      idat (16) = rx
      idat (16) = idat (16) + 127
      idat (17) = ry
      idat (17) = idat (17) + 127
      idat (18) = id2m
      idat (19) = j2m
      idat (20) = h2m
      idat (21) = k2m
      idat (22) = ph 
      idat (22) = idat (22) + 127
      idat (23) = cc
      idat (23) = idat (23) + 127

      return

 99   ierra=1

      RETURN
      END







c
c
c     subroutine ucac2s
c
c
c     Picks up UCAC2 stars with speedy extraction
c
c
c     Last update:   15/Feb/2013.  M. Assafin
c
c

      subroutine ucac2s (iu2z,iu2i,inuu2,ir1u2,ir2u2,uraiz,ecom,alfmin,
     ?alfmax,delmin,delmax,nobs,cora,code,codj,ucra,ucde,id,oldra,oldde,
     ?rac,dec,nok,rmx,rmy,rmx2,rmy2,iall)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra

      dimension cora(iall),code(iall),ucra(iall),
     ?ucde(iall),codj(iall),id(iall)

      dimension oldra(iall),oldde(iall),pmx(iall),
     ?pmy(iall)

      dimension idat(23)

      dimension inuu2(iu2z,iu2i),ir1u2(iu2z,iu2i),ir2u2(iu2z,iu2i)

      dimension jxmin(2),jxmax(2),cxmin(2),cxmax(2)


      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 uraiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/


C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      un=95

      idin50=50
      izon=iu2z
      bin=0.1d0
      nbin=iu2i

      nok=0

      ndim=iall

      ebox=ecom**2

c

      do i=1,2
      jxmin(i)=0
      jxmax(i)=0
      cxmin(i)=0
      cxmax(i)=0
      enddo



c

      write (*,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax
      write (46,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax

      write (*,*)
      write (46,*)




c
c     Sets limits for indexed zone reading
c

      xx1=alfmin
      xx2=alfmax

      inde=1

      jxmin(1)=xx1/bin+1
      jxmax(1)=xx2/bin+1

      cxmin(1)=xx1*54.d6
      cxmax(1)=xx2*54.d6




c
c     Reads UCAC2 declination zones of 0.5 degrees
c

      dfaixa=delmin-0.5d0
      decmax=delmax

      nest=0


      do 30 k=1,288

      dfaixa=dfaixa+0.5d0


      if (dfaixa-decmax.gt.0.5d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.5d0+1.0d0

      jaca=j



c
c     UCAC2 North declination limit
c
c     For now, just zones up to z288 will do
c

      if (j.gt.izon) go to 30

c
c     Mounts 0.5 degrees zone file name 
c


      ifaixa=''
      ifaixa=uraiz

      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC2 ',a54)


c
c     Reads current zone in the loop
c

      open (95,file=ifaixa,access="direct",form="unformatted",recl=44)


      do 20 in=1,inde

      do 19 nn=jxmin(in),jxmax(in)


      if(inuu2(jaca,nn).eq.0) go to 19

      i1=ir1u2(jaca,nn)
      i2=ir2u2(jaca,nn)


      do 18 n=i1,i2


c
c     Picks up star data
c    

      call readu2 (un,n,idat,ierra)


      ra=idat(1)


      if (ra.lt.cxmin(in)) go to 18
      if (ra.gt.cxmax(in)) go to 18

      nest=nest+1


      oldra(nest)=idat(1)
      oldde(nest)=idat(2)
      pmx(nest)=idat(12)
      pmy(nest)=idat(13)


      if (nest.eq.ndim) then
      close (95)
      go to 35
      endif

 18   continue


 19   continue

 20   continue

      close (95)

c

 30   continue

c

 35   continue

c




      if (nest.ge.ndim) then

      write (*,*)
      write (46,*)

      write (*,*) 'Max. number of UCAC2 stars reached: ',nest,ndim
      write (46,*) 'Max. number of UCAC2 stars reached: ',nest,ndim

      write (*,*)
      write (46,*)

      nest=ndim

      write (*,*) 'UCAC2 stars candidates reduced to: ',nest
      write (46,*) 'UCAC2 stars candidates reduced to: ',nest

      write (*,*)
      write (46,*)


      else

      write (*,*)
      write (46,*)

      write (*,*) 'UCAC2 candidate stars: ',nest
      write (46,*) 'UCAC2 candidate stars: ',nest

      write (*,*)
      write (46,*)

      endif


c
c     Finds if it is a match: measured star x UCAC2 star
c


      write (*,*)
      write (*,*) 'Matching xy field stars vs. UCAC2 stars. Please wait 
     ?...'
      write (*,*)

      write (46,*)
      write (46,*)'Matching xy field stars vs. UCAC2 stars. Please wait 
     ?...'
      write (46,*)

c



      rac=0.d0
      dec=0.d0

      rmx=0.d0
      rmy=0.d0
     
      rmx2=0.d0
      rmy2=0.d0


      do 60 i=1,nest
      do 50 ko=1,nobs

      if (id(ko).ge.0) go to 50

      epoj=2000D0+(codj(ko)-2451545D0)/365.25D0

      rauc=oldra(i)+pmx(i)*(epoj-2000d0)/10d0
      deuc=oldde(i)+pmy(i)*(epoj-2000d0)/10d0

      rauc=rauc/1000d0
      deuc=deuc/1000d0

      rauc=rauc/3600d0
      deuc=deuc/3600d0

      dx=(15.d0*cora(ko)-rauc)*dcos(deuc*grarad)*3600.d0
      dy=(code(ko)-deuc)*3600.d0


      dis=dx**2+dy**2

      if (dis.lt.ebox) then

      ucra(ko)=rauc/15.d0
      ucde(ko)=deuc
      id(ko)=ko

      rac=rac+rauc
      dec=dec+deuc

      rmx=rmx+dx
      rmy=rmy+dy
      rmx2=rmx2+dx**2
      rmy2=rmy2+dy**2

c
c
c     debug alfa delta
c
c     write (*,*) 'oldra, oldde, codj = ',oldra(ko),oldde(ko),codj(ko)
c
c     ra=rauc
c     de=deuc
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,

      nok=nok+1

      go to 60

      endif

c

 50   continue

 60   continue

      rac=rac/nok
      dec=dec/nok

      call desvio (nok,rmx,rmx2)
      call desvio (nok,rmy,rmy2)


      write (*,*)
      write (*,*) 'No. UCAC2 common stars = ',nok
      write (*,*)

      write (46,*)
      write (46,*) 'No. UCAC2 common stars = ',nok
      write (46,*)

c     stop

      return
      end









c
c
c     subroutine ucac4
c
c
c     Picks up UCAC4 stars with no speedy extraction
c
c
c     Last update:   15/Feb/2013.  M. Assafin
c
c

      subroutine ucac4 (uraiz,ecom,alfmin,alfmax,delmin,delmax,nobs,
     ?cora,code,codj,ucra,ucde,id,oldra,oldde,rac,dec,nok,rmx,rmy,rmx2,
     ?rmy2,iall)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra,reclen

      LOGICAL  bf, eozf
      INTEGER*4 ran,spdn, id2m,rnm, mcf, rnz
      INTEGER*2 magm,maga, cepra,cepdc, pmra2,pmdc2
     .         ,jmag,hmag,kmag, apasm(5), zn2
      INTEGER*1 sigmag, sigra,sigdc, sigpmr,sigpmd
      INTEGER*1 ojt,dsf, na1,nu1,us1, apase(5), gcflg
      INTEGER*1 icqflg(3), q2mflg(3), leda,x2m

      dimension irnm(25),ipmrc(25),ipmd(25)


      dimension cora(iall),code(iall),ucra(iall),
     ?ucde(iall),codj(iall),id(iall)

      dimension oldra(iall),oldde(iall),pmx(iall),
     ?pmy(iall)



      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 uraiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/

 
      data irnm /1,2,200137,200168,200169,200229,200400,200503,200530,
     ?200895,201050,201349,201526,201550,201567,201633,201803,249921,
     ?249984,268357,80118783,93157181,106363470,110589580,113038183/ 

      data ipmrc /41087,41558,-37758,-36004,-36782,39624,65051,56347,
     ?67682,-22401,-7986,22819,-5803,40033,41683,-44099,34222,-10015,
     ?-9994,5713,32962,-22393,-37060,-38420,10990/

      data ipmd /31413,32586,7655,9521,4818,-25374,-57308,-23377,13275,
     ?-34203,103281,53694,-47659,-58151,32691,9416,-15989,-35427,-35419,
     ?-36943,5639,-34199,-11490,-27250,-51230/



C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      un=95

      idin50=50
      izon=900
      reclen=78

c

      nok=0

      ndim=iall

      ebox=ecom**2

c

      write (*,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax
      write (46,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax

      write (*,*)
      write (46,*)



      alma=alfmax*15.d0
      almi=alfmin*15.d0

c
c     Reads UCAC4 declination zones of 0.2 degrees
c

      dfaixa=delmin-0.2d0
      decmax=delmax

      nest=0


      do 30 k=1,izon

      dfaixa=dfaixa+0.2d0


      if (dfaixa-decmax.gt.0.2d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.2d0+1.0d0

c
c     UCAC4 North declination limit
c
c     Only zones up to z900
c

      if (j.gt.izon) go to 30

c
c     Mounts 0.2 degrees zone file name 
c


      ifaixa=''
      ifaixa=uraiz

      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC4 ',a54)



c
c     Reads current zone in the loop
c




      open (95,file=ifaixa,access="direct",recl=reclen)


      n=0
 20   n=n+1


c

      READ (un,REC=n,ERR=25)                      !          sum = 78
     .     ran,spdn, magm,maga, sigmag,ojt,dsf    !  8 + 4 +  3  = 15
     .    ,sigra,sigdc, na1,nu1,us1               !  2 + 3       =  5
     .    ,cepra,cepdc, pmra2,pmdc2,sigpmr,sigpmd !  4 + 4 +  2  = 10
     .    ,id2m, jmag,hmag,kmag, icqflg, q2mflg   !  4 + 6 +  6  = 16
     .    ,apasm, apase, gcflg                    ! 10 + 5 +  1  = 16
     .    ,mcf, leda,x2m, rnm                     !  4 + 2 +  4  = 10
     .    ,zn2, rnz                               !  2 + 4       =  6


      ra=ran
      rab=ra/3.6d6

c
c     Accelerates match
c


      if (rab.lt.almi) go to 20
      if (rab.gt.alma) go to 20

      nest=nest+1

      oldra(nest)=ran
      oldde(nest)=spdn

      de=spdn


      if (pmra2.eq.32767 .and. pmdc2.eq.32767) then

      do m=1,25 
      if (rnm.eq.irnm(m)) go to 23
      enddo

 23   continue

      pmx(nest)=ipmrc(m)/dcos((de/3.6d6-90d0)*grarad)
      pmy(nest)=ipmd(m)


      else

      pmx(nest)=pmra2/dcos((de/3.6d6-90d0)*grarad)
      pmy(nest)=pmdc2

      endif


      if (nest.eq.ndim) then
      close (95)
      go to 36
      endif

      go to 20


 25   close (95)

c

 30   continue

c

 35   continue

c

 36   continue



      if (nest.ge.ndim) then

      write (*,*)
      write (46,*)

      write (*,*) 'Max. number of UCAC4 stars reached: ',nest,ndim
      write (46,*) 'Max. number of UCAC4 stars reached: ',nest,ndim

      write (*,*)
      write (46,*)

      nest=ndim

      write (*,*) 'UCAC4 stars candidates reduced to: ',nest
      write (46,*) 'UCAC4 stars candidates reduced to: ',nest

      write (*,*)
      write (46,*)


      else

      write (*,*)
      write (46,*)

      write (*,*) 'UCAC4 candidate stars: ',nest
      write (46,*) 'UCAC4 candidate stars: ',nest

      write (*,*)
      write (46,*)

      endif


c
c     Finds if it is a match: measured star x UCAC4 star
c


      write (*,*)
      write (*,*) 'Matching xy field stars vs. UCAC4 stars. Please wait 
     ?...'
      write (*,*)

      write (46,*)
      write (46,*)'Matching xy field stars vs. UCAC4 stars. Please wait 
     ?...'
      write (46,*)

c



      rac=0.d0
      dec=0.d0

      rmx=0.d0
      rmy=0.d0
     
      rmx2=0.d0
      rmy2=0.d0


      do 60 i=1,nest
      do 50 ko=1,nobs

      if (id(ko).ge.0) go to 50

      epoj=2000D0+(codj(ko)-2451545D0)/365.25D0

      rauc=oldra(i)+pmx(i)*(epoj-2000d0)/10d0
      deuc=oldde(i)+pmy(i)*(epoj-2000d0)/10d0

      rauc=rauc/3.6d6 
      deuc=deuc/3.6d6 
      deuc=deuc-90.d0


      dx=(15.d0*cora(ko)-rauc)*dcos(deuc*grarad)*3600.d0
      dy=(code(ko)-deuc)*3600.d0


      dis=dx**2+dy**2

      if (dis.lt.ebox) then

      ucra(ko)=rauc/15.d0
      ucde(ko)=deuc
      id(ko)=ko

      rac=rac+rauc
      dec=dec+deuc

      rmx=rmx+dx
      rmy=rmy+dy
      rmx2=rmx2+dx**2
      rmy2=rmy2+dy**2

c
c
c     debug alfa delta
c
c     write (*,*) 'oldra, oldde, codj = ',oldra(ko),oldde(ko),codj(ko)
c
c     ra=rauc
c     de=deuc
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,

      nok=nok+1

      go to 60

      endif

c

 50   continue

 60   continue

      rac=rac/nok
      dec=dec/nok

      call desvio (nok,rmx,rmx2)
      call desvio (nok,rmy,rmy2)


      write (*,*)
      write (*,*) 'No. UCAC4 common stars = ',nok
      write (*,*)

      write (46,*)
      write (46,*) 'No. UCAC4 common stars = ',nok
      write (46,*)

c     stop

      return
      end








c
c
c     subroutine ucac4s
c
c
c     Picks up UCAC4 stars with speedy extraction
c
c
c     Last update:   15/Feb/2013.  M. Assafin
c
c

      subroutine ucac4s (iu4z,iu4i,inuu4,ir1u4,ir2u4,uraiz,ecom,alfmin,
     ?alfmax,delmin,delmax,nobs,cora,code,codj,ucra,ucde,id,oldra,oldde,
     ?rac,dec,nok,rmx,rmy,rmx2,rmy2,iall)

      implicit real*8 (a-h,o-z)

      integer*2 un,ierra,reclen

      LOGICAL  bf, eozf
      INTEGER*4 ran,spdn, id2m,rnm, mcf, rnz
      INTEGER*2 magm,maga, cepra,cepdc, pmra2,pmdc2
     .         ,jmag,hmag,kmag, apasm(5), zn2
      INTEGER*1 sigmag, sigra,sigdc, sigpmr,sigpmd
      INTEGER*1 ojt,dsf, na1,nu1,us1, apase(5), gcflg
      INTEGER*1 icqflg(3), q2mflg(3), leda,x2m

      dimension irnm(25),ipmrc(25),ipmd(25)

      dimension inuu4(iu4z,iu4i),ir1u4(iu4z,iu4i),ir2u4(iu4z,iu4i)

      dimension jxmin(2),jxmax(2),cxmin(2),cxmax(2)


      dimension cora(iall),code(iall),ucra(iall),
     ?ucde(iall),codj(iall),id(iall)

      dimension oldra(iall),oldde(iall),pmx(iall),
     ?pmy(iall)



      CHARACTER *1 MAIS,MENOS,ip,isig
      character *54 ifaixa
      character *50 uraiz

      DATA MENOS/'-'/
      DATA MAIS/'+'/
      data ip/'z'/

 
      data irnm /1,2,200137,200168,200169,200229,200400,200503,200530,
     ?200895,201050,201349,201526,201550,201567,201633,201803,249921,
     ?249984,268357,80118783,93157181,106363470,110589580,113038183/ 

      data ipmrc /41087,41558,-37758,-36004,-36782,39624,65051,56347,
     ?67682,-22401,-7986,22819,-5803,40033,41683,-44099,34222,-10015,
     ?-9994,5713,32962,-22393,-37060,-38420,10990/

      data ipmd /31413,32586,7655,9521,4818,-25374,-57308,-23377,13275,
     ?-34203,103281,53694,-47659,-58151,32691,9416,-15989,-35427,-35419,
     ?-36943,5639,-34199,-11490,-27250,-51230/



C
C     initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      un=95

      idin50=50
      izon=iu4z
      reclen=78

      bin=0.25d0
      nbin=iu4i

c

      do i=1,2
      jxmin(i)=0
      jxmax(i)=0
      cxmin(i)=0
      cxmax(i)=0
      enddo


c

      nok=0

      ndim=iall

      ebox=ecom**2

c

      write (*,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax
      write (46,*) 'error box, (RA,DEC) min max = ',ecom,alfmin,alfmax,
     ?delmin,delmax

      write (*,*)
      write (46,*)





c
c     Sets limits for indexed zone reading
c

      xx1=alfmin
      xx2=alfmax

      inde=1

      jxmin(1)=xx1/bin+1
      jxmax(1)=xx2/bin+1

      cxmin(1)=xx1*3600.d3
      cxmax(1)=xx2*3600.d3



c
c     Reads UCAC4 declination zones of 0.2 degrees
c

      dfaixa=delmin-0.2d0
      decmax=delmax

      nest=0


      do 30 k=1,izon

      dfaixa=dfaixa+0.2d0


      if (dfaixa-decmax.gt.0.2d0) go to 35

      if (dfaixa.lt.-90.d0) dfaixa=-90.d0

      j=(dfaixa+90.d0)/0.2d0+1.0d0

      jaca=j

c
c     UCAC4 North declination limit
c
c     Only zones up to z900
c

      if (j.gt.izon) go to 30

c
c     Mounts 0.2 degrees zone file name 
c


      ifaixa=''
      ifaixa=uraiz

      do l=1,idin50
      if (ifaixa(l:l).eq.' ') go to 10
      enddo

 10   write(ifaixa(l:l+3),'(a1,i3.3)') ip,j


      write (*,14) ifaixa
 14   format(1x,'UCAC4 ',a54)



c
c     Reads current zone in the loop
c




      open (95,file=ifaixa,access="direct",recl=reclen)




      do 20 in=1,inde

      do 19 nn=jxmin(in),jxmax(in)


      if(inuu4(jaca,nn).eq.0) go to 19

      i1=ir1u4(jaca,nn)
      i2=ir2u4(jaca,nn)



      do 18 n=i1,i2



c

c     READ (un,REC=n,ERR=25)                      !          sum = 78
      READ (un,REC=n)                             !          sum = 78
     .     ran,spdn, magm,maga, sigmag,ojt,dsf    !  8 + 4 +  3  = 15
     .    ,sigra,sigdc, na1,nu1,us1               !  2 + 3       =  5
     .    ,cepra,cepdc, pmra2,pmdc2,sigpmr,sigpmd !  4 + 4 +  2  = 10
     .    ,id2m, jmag,hmag,kmag, icqflg, q2mflg   !  4 + 6 +  6  = 16
     .    ,apasm, apase, gcflg                    ! 10 + 5 +  1  = 16
     .    ,mcf, leda,x2m, rnm                     !  4 + 2 +  4  = 10
     .    ,zn2, rnz                               !  2 + 4       =  6




      ra=ran

      if (ra.lt.cxmin(in)) go to 18
      if (ra.gt.cxmax(in)) go to 18



      nest=nest+1


      oldra(nest)=ran
      oldde(nest)=spdn

      de=spdn


      if (pmra2.eq.32767 .and. pmdc2.eq.32767) then

      do m=1,25 
      if (rnm.eq.irnm(m)) go to 16
      enddo

 16   continue

      pmx(nest)=ipmrc(m)/dcos((de/3.6d6-90d0)*grarad)
      pmy(nest)=ipmd(m)


      else

      pmx(nest)=pmra2/dcos((de/3.6d6-90d0)*grarad)
      pmy(nest)=pmdc2

      endif


      if (nest.eq.ndim) then
      close (95)
      go to 35
      endif

      go to 20


 18   continue

 19   continue

 20   continue

c
      close (95)

c

 30   continue

c

 35   continue

c




      if (nest.eq.ndim) then

      write (*,*)
      write (46,*)

      write (*,*) 'Max. number of UCAC4 stars reached: ',nest,ndim
      write (46,*) 'Max. number of UCAC4 stars reached: ',nest,ndim

      write (*,*)
      write (46,*)

      nest=ndim

      write (*,*) 'UCAC4 stars candidates reduced to: ',nest
      write (46,*) 'UCAC4 stars candidates reduced to: ',nest

      write (*,*)
      write (46,*)


      else

      write (*,*)
      write (46,*)

      write (*,*) 'UCAC4 candidate stars: ',nest
      write (46,*) 'UCAC4 candidate stars: ',nest

      write (*,*)
      write (46,*)

      endif


c
c     Finds if it is a match: measured star x UCAC4 star
c


      write (*,*)
      write (*,*) 'Matching xy field stars vs. UCAC4 stars. Please wait 
     ?...'
      write (*,*)

      write (46,*)
      write (46,*)'Matching xy field stars vs. UCAC4 stars. Please wait 
     ?...'
      write (46,*)

c



      rac=0.d0
      dec=0.d0

      rmx=0.d0
      rmy=0.d0
     
      rmx2=0.d0
      rmy2=0.d0


      do 60 i=1,nest
      do 50 ko=1,nobs

      if (id(ko).ge.0) go to 50

      epoj=2000D0+(codj(ko)-2451545D0)/365.25D0

      rauc=oldra(i)+pmx(i)*(epoj-2000d0)/10d0
      deuc=oldde(i)+pmy(i)*(epoj-2000d0)/10d0

      rauc=rauc/3.6d6 
      deuc=deuc/3.6d6 
      deuc=deuc-90.d0


      dx=(15.d0*cora(ko)-rauc)*dcos(deuc*grarad)*3600.d0
      dy=(code(ko)-deuc)*3600.d0


      dis=dx**2+dy**2

      if (dis.lt.ebox) then

      ucra(ko)=rauc/15.d0
      ucde(ko)=deuc
      id(ko)=ko

      rac=rac+rauc
      dec=dec+deuc

      rmx=rmx+dx
      rmy=rmy+dy
      rmx2=rmx2+dx**2
      rmy2=rmy2+dy**2

c
c
c     debug alfa delta
c
c     write (*,*) 'oldra, oldde, codj = ',oldra(ko),oldde(ko),codj(ko)
c
c     ra=rauc
c     de=deuc
c     IAH=ra
c     AM=(ra-IAH)*60.D0
c     IAM=AM
c     SA =(AM-IAM)*60.D0
c     IF (de.LT.0.D0) THEN
c     ISIG=MENOS
c     de=-de
c     ELSE
c     ISIG=MAIS
c     ENDIF
c     IDG=de
c     DM=(de-IDG)*60.D0
c     IDM=DM
c     DS=(DM-IDM)*60.D0

c     write(*,*) IAH,IAM,SA,'  ',ISIG,IDG,IDM,DS,

      nok=nok+1

      go to 60

      endif

c

 50   continue

 60   continue

      rac=rac/nok
      dec=dec/nok

      call desvio (nok,rmx,rmx2)
      call desvio (nok,rmy,rmy2)


      write (*,*)
      write (*,*) 'No. UCAC4 common stars = ',nok
      write (*,*)

      write (46,*)
      write (46,*) 'No. UCAC4 common stars = ',nok
      write (46,*)

c     stop

      return
      end




c
c
c
c
c     subroutine rotsky 
c
c     Computes small rotation angles (a1,a2,a3) from positional differences
c
c     "P1 - P2"
c
c     in such a way that that positive angles (a1,a2,a3) makes positions
c     P1 rotate toward P2.
c     
c     (see "COMPARISON OF VLBI CELESTIAL REFERENCE FRAMES" DE E.F. ARIAS,
c      M. FEISSEL E J.-F. LESTRADE, 1988 A.& A 199, P357-363.)
c
c
c     Besides (a1,a2,a3), their respective rms errors are also computed.
c
c
c
c
c      Last update: 19/Apr/2010  M. Assafin
c
c

      subroutine rotsky (corte,nobs,cora,code,ucra,ucde,ranew,denew,id,
     ?a1,a2,a3,e1,e2,e3,rmx,rmy,rmx2,rmy2,rmgx,rmgy,rmgx2,rmgy2,nc,iall)


      IMPLICIT REAL *8 (A-H,O-Z)

      dimension ARRAY(21,21),AMAT(3,3),B(3),SOL(3),ERRO(3)

      dimension cora(iall),code(iall),ucra(iall),
     ?ucde(iall),ranew(iall),denew(iall),itira(iall),
     ?xr(iall),yr(iall),id(iall)

      REX(Q,W,X,Y,Z)=X*(DSIN(W))*DCOS(Q)+Y*(DSIN(W))*DSIN(Q)-Z*DCOS(W)
      REY(Q,W,X,Y)=-X*DSIN(Q)+Y*DCOS(Q)

      COMMON /A7/ARRAY
      COMMON /A14/IERRO

      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

      DET=1.D0
      NORDER=3

      idiobs=iall

C
C     Emptying vectors
C

      do i=1,idiobs
      itira(i)=0
      xr(i)=0.d0
      yr(i)=0.d0
      enddo

c
c     Pre-elimination of outliers with large (O-C) position differences
c


      do 5 i=1,nobs

      if (id(i).le.0) go to 5

      resalf=3600.d0*(ucra(i)-cora(i))*15.d0*dcos(grarad*code(i))
      resdel=3600.d0*(ucde(i)-code(i))

      dis=dsqrt(resalf**2+resdel**2)

c     write (*,*) 'resalf,resdel = ',resalf,resdel

      if (dis.gt.corte) itira(i)=1

 5    continue

c
c

 1    continue

c

      do i=1,21
      do j=1,21
      array(i,j)=0.d0
      enddo
      enddo

      do 2 i=1,norder
      b(i)=0.d0
      sol(i)=0.d0
      erro(i)=0.d0
      do 2 j=1,norder
      amat(i,j)=0.d0
 2    continue




c
c     Defining the sums
c


      do 10 k=1,nobs

      if (id(k).le.0) go to 10

      if (itira(k).ne.0) go to 10

      resalf=grarad*(ucra(k)-cora(k))*15.d0*dcos(grarad*code(k))
      resdel=grarad*(ucde(k)-code(k))

      a=cora(k)*15.d0*grarad
      d=code(k)*grarad

      AMAT(1,1)=AMAT(1,1)-(DSIN(D)*DCOS(A))**2-DSIN(A)**2
      AMAT(1,2)=AMAT(1,2)+DCOS(A)*DSIN(A)*DCOS(D)**2
      AMAT(1,3)=AMAT(1,3)+DCOS(D)*DSIN(D)*DCOS(A)
      AMAT(2,1)=AMAT(2,1)+DCOS(A)*DSIN(A)*DCOS(D)**2
      AMAT(2,2)=AMAT(2,2)-1.d0+(DCOS(D)*DSIN(A))**2
      AMAT(2,3)=AMAT(2,3)+DSIN(D)*DSIN(A)*DCOS(D)
      AMAT(3,1)=AMAT(3,1)+DSIN(D)*DCOS(D)*DCOS(A)
      AMAT(3,2)=AMAT(3,2)+DSIN(D)*DCOS(D)*DSIN(A)
      AMAT(3,3)=AMAT(3,3)-DCOS(D)**2

      B(1)=B(1)-RESALF*DSIN(D)*DCOS(A)+RESDEL*DSIN(A)
      B(2)=B(2)-RESALF*DSIN(D)*DSIN(A)-RESDEL*DCOS(A)
      B(3)=B(3)+RESALF*DCOS(D)


 10   continue


c
c     Normalizing the curvature matrix by the diagonal elements
c

      do 21 i=1,norder
      do 20 j=1,norder
      array(i,j)=amat(i,j)/dsqrt(amat(i,i)*amat(j,j))
 20   continue
 21   continue

c
c     Inverting the matrix
c

      CALL MATINV (norder,det)

      if (ierro.eq.1) then
      write (*,*) 'L.S. rotation fit crashed. Aborting execution.'
      write (46,*) 'L.S. rotation fit crashed. Aborting execution.'
      stop
      endif

C
C     Computing rotation angles a1, a2 and a3
C

      DO 31 L=1,norder
      DO 30 K=1,norder
      sol(L)=sol(L)+array(L,K)*b(K)/dsqrt(amat(L,L)*amat(K,K))
 30   continue
 31   continue

c
c     Computing rotation angle errors e1, e2 and e3
c

      rmgx=0.D0
      rmgy=0.D0
      rmgx2=0.D0
      rmgy2=0.D0

      rmx=0.D0
      rmy=0.D0
      rmx2=0.D0
      rmy2=0.D0

      nc=0

      do 40 i=1,nobs

      if (id(i).le.0) go to 40
      if (itira(i).ne.0) go to 40

      nc=nc+1

      resalf=grarad*(ucra(i)-cora(i))*15.d0*dcos(grarad*code(i))
      resdel=grarad*(ucde(i)-code(i))

      arad=grarad*15.d0*cora(i)
      drad=grarad*code(i)

      RMX=RMX+RESALF
      RMY=RMY+RESDEL

      RMX2=RMX2+RESALF**2
      RMY2=RMY2+RESDEL**2

      aux=RESALF-REX(ARAD,DRAD,SOL(1),SOL(2),SOL(3))
      auy=RESDEL-REY(ARAD,DRAD,SOL(1),SOL(2))

      xr(i)=aux
      yr(i)=auy

      rmgx=rmgx+aux
      rmgy=rmgy+auy

      rmgx2=rmgx2+aux**2
      rmgy2=rmgy2+auy**2

 40   continue

c
c     Statistics for elimination of outliers
c


      call desvio (nc,rmgx,rmgx2)
      call desvio (nc,rmgy,rmgy2)

      R2=rmgx2**2+rmgy2**2
      rr=dsqrt(r2)

      rr3=3.d0*rr

      ifora=0

      do 50 i=1,nobs

      if (id(i).le.0) go to 50
      if (itira(i).ne.0) go to 50

      rr=dsqrt(xr(i)**2+yr(i)**2)

      if (rr.gt.rr3) then
      itira(i)=1
      ifora=ifora+1
      endif

 50   continue

c
c     Go back to compute rotation if outliers still exist
c

      if (ifora.ne.0) go to 1

c
c     Statistics of rotation computation
c

      call desvio (nc,rmx,rmx2)
      call desvio (nc,rmy,rmy2)

      rmx=rmx*radgra*3600.d0
      rmy=rmy*radgra*3600.d0
      rmx2=rmx2*radgra*3600.d0
      rmy2=rmy2*radgra*3600.d0


      ERRO(1)=DSQRT(R2*ARRAY(1,1)/amat(1,1))
      ERRO(2)=DSQRT(R2*ARRAY(2,2)/amat(2,2))
      ERRO(3)=DSQRT(R2*ARRAY(3,3)/amat(3,3))

      r2=dsqrt(r2)
c

      a1=sol(1)
      a2=sol(2)
      a3=sol(3)

      e1=ERRO(1)
      e2=ERRO(2)
      e3=ERRO(3)

c
c     Rotates positions toward the reference catalog.
c     Computes new position offsets after rotation.
c

      rmgx=0.D0
      rmgy=0.D0
      rmgx2=0.D0
      rmgy2=0.D0

      nc=0

      do 60 i=1,nobs

      call roda (a1,a2,a3,cora(i),code(i),x,y)

      if (id(i).le.0) go to 55
      if (itira(i).ne.0) go to 55

      nc=nc+1

      aux=3600.d0*(ucra(i)-x)*15.d0*dcos(grarad*y)
      auy=3600.d0*(ucde(i)-y)

      rmgx=rmgx+aux
      rmgy=rmgy+auy

      rmgx2=rmgx2+aux**2
      rmgy2=rmgy2+auy**2

      

 55   ranew(i)=x
      denew(i)=y

 60   continue
c

      call desvio (nc,rmgx,rmgx2)
      call desvio (nc,rmgy,rmgy2)

c

      return
      end


c
c
c
c
c     subroutine roda
c
c     Given the rotation angles (a1,a2,a3) and (RA,DEC), it returns
c     the new (RA,DEC) rotated.
c
c
c             ( dx' )   (  001  +a3  -a2 )   ( dx" )
c             ( dy' ) = (  -a3  001  +a1 ) . ( dy" )
c             ( dz' )   (  +a2  -a1  001 )   ( dz" )
c
c
c
c      Last update: 04/Dec/2007  M. Assafin
c
c

      subroutine roda (a1,a2,a3,ain,din,aout,dout)


      IMPLICIT REAL *8 (A-H,O-Z)

      dimension dm(3,3),cv(3),rcv(3)

 
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      a=ain*15.d0*grarad
      d=din*grarad

c
c     director cossines of position
c

      cv(1)=dcos(d)*dcos(a)
      cv(2)=dcos(d)*dsin(a)
      cv(3)=dsin(d)

c
c     Matrix of rotation
c     


      dm(1,1)=1.d0
      dm(1,2)=a3
      dm(1,3)=-a2
      dm(2,1)=-a3
      dm(2,2)=1.d0
      dm(2,3)=a1
      dm(3,1)=a2
      dm(3,2)=-a1
      dm(3,3)=1.d0

c
c     Rotated director cossines
c

      do i=1,3
      rcv(i)=0.d0
      enddo

c
      do i=1,3
      do j=1,3

      rcv(i)=rcv(i)+dm(i,j)*cv(j)

      enddo
      enddo

c
c     Rotated (RA,DEC)
c

      d=dasin(rcv(3))

      d=d*radgra


c
c     Computing alpha quadrant
c

      aa=dabs(datan2(dabs(rcv(2)),dabs(rcv(1))))

      if (rcv(1).lt.0.d0) then
      isigx=-1
      else
      isigx=+1
      endif


      if (rcv(2).lt.0.d0) then
      isigy=-1
      else
      isigy=+1
      endif

c
c     Computing RA (0h - 24h) from quadrant determination
c


      if (isigx.gt.0 .and. isigy.gt.0) a=aa
      if (isigx.lt.0 .and. isigy.gt.0) a=pi-aa
      if (isigx.lt.0 .and. isigy.lt.0) a=pi+aa
      if (isigx.gt.0 .and. isigy.lt.0) a=2.d0*pi-aa

     
c

      a=a*radgra/15.d0
c

      aout=a
      dout=d
c

      return
      end


c
c
c   
c     subrotina estat
c
c     Retrives targets from xy files and output statistcs.
c
c
c     modo   -> type of xy file
c               modo = 1, catalog xy file (different Julian dates)
c               modo = 2, individual xy files (same Julian dates in all
c                         entries
c
c     ecom   -> error box for (RA,DEC) identification (arcsec)
c
c     boxepo -> time error box for target identification
c
c     if1 -> xy file
c     if2 -> output of results
c     if3 -> list of target input data
c
c
c     Last update: 09/Jul/2008  M. Assafin
c
c

      subroutine estat (modo,ecom,if1,if2,if3)



      IMPLICIT REAL *8 (A-H,O-Z)


      character*50 infits
      character*20 ichfil,mchobj,iobalv
      character*1 isig,menos

      data menos/'-'/

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0

c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI

c

      dj1=2400000.5D0

c
c     (RA,DEC) box error
c

      box=ecom

c
c     Time error box (default box=1 second)
c

      boxepo=1.d0

c
c     Converts boxepo to Julian date (fraction of day)
c

      boxepo=boxepo/(3600.d0*24.d0)

c

      ipegou=0

c

      write (*,*)
      write (*,*)'offsets (RA,DE), Sigma(RA,DE), Ncat, Gauss_error(x,y),
     ?date, exptime, filter, target'
      write (*,*)


      write (46,*)
      write(46,*)'offsets (RA,DE), Sigma(RA,DE), Ncat, Gauss_error(x,y),
     ?date, exptime, filter, target'
      write (46,*)

c
c5    format(a1)
c

c

 100  rewind (if1)

      read (if3,101,end=200) iah,iam,as,isig,idg,idm,ds,datalv,
     ?iobalv
 101  format(1x,i2,1x,i2,1x,f9.6,1x,a1,i2,1x,i2,1x,f8.5,1x,f16.8,
     ?1x,a20)

      rafat=hmsgms(iah,iam,as)
      defat=hmsgms(idg,idm,ds)
      if (isig.eq.menos) defat=-defat



c
c     Checks Julian Date of target and stars for identification
c

      if (modo.eq.1) go to 19

c

 18   continue

      read (if1,10,err=18,end=12) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy


c
c     Targets coordiantes with negative Julian Date can be searched in
c     any given xy field regardless of the instant of observation;
c     otherwise, target must fulfill (RA,DEC) error box AND time error
c     box
c



      if (datalv.gt.0.d0) then

      dtemp=dabs(datalv-dj)

      if (dtemp.gt.boxepo) go to 12

      endif

      rewind (if1)


c
c     Until here, ok for serach with or without time error box
c



 19   continue

      read (if1,10,err=19,end=12) xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,mchobj,nx,ny,ncom,aux,auy


 10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?3(1x,i5),2(1x,f7.3))


      if (modo.eq.1 .and. datalv.gt.0.d0) then

      dtemp=dabs(datalv-dj)

      if (dtemp.gt.boxepo) go to 19

      endif


c
c     Error box in (alpha,delta) for target identification
c


      dx=(ra-rafat)*dcos(defat*grarad)*3600.d0*15.d0
      dy=(de-defat)*3600.d0

      if (dabs(dx).gt.box) go to 19
      if (dabs(dy).gt.box) go to 19

      ipegou=1

      if (modo.eq.2) then

      write (*,16) dx,dy,alfsiu,delsiu,nfinau,ex,ey,iuth,
     ?iutm,sut,iutano,iutmes,iutdia,iexps,ichfil,iobalv

      write (46,16) dx,dy,alfsiu,delsiu,nfinau,ex,ey,iuth,
     ?iutm,sut,iutano,iutmes,iutdia,iexps,ichfil,iobalv

 16   format(4(1x,f7.3),2x,i5,2x,2(1x,f7.3),2x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,2x,i4,3x,a20,2x,a20)

      else

      djm=dj-dj1

      call iau_jd2cal (dj1,djm,iutan0,iutme0,iutdi0,fd,jjj)

      hora0=fd*24.d0
      iuth0=hora0
      iutm0=(hora0-iuth0)*60.d0
      sut0 =((hora0-iuth0)*60.d0-iutm0)*60.d0

      write (*,16) dx,dy,alfsiu,delsiu,nfinau,ex,ey,iuth0,
     ?iutm0,sut0,iutan0,iutme0,iutdi0,iexps,ichfil,iobalv

      write (46,16) dx,dy,alfsiu,delsiu,nfinau,ex,ey,iuth0,
     ?iutm0,sut0,iutan0,iutme0,iutdi0,iexps,ichfil,iobalv


      endif


c     
      if (ncom.eq.0) ncom=1
c

      write (if2,11) dx,dy,xob,yob,seng,altu,fgcc,fumag,fumag2,
     ?xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,ermgh,
     ?ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,iuth,iutm,sut,iutano,
     ?iutmes,iutdia,dj,iexps,ichfil,infits,iobalv,nx,ny,ncom,aux,auy


 11   format(2(1x,f7.3),2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),
     ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
     ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
     ?1x,a20,3(1x,i5),2(1x,f7.3))



 12   go to 100

c


 200  continue


      return
      end






c 
c
c     Subrotine tanpla
c
c
c     Tangent plane techinique. Reduces (RA,DEC) toward a reference set of
c     common (RA,DEC) using standard tnagent plane polynomial relations.
c
c     Here it is assumed that both (RA,DEC) sets are at the same epoch.
c
c
c     Allowed polynomial fits: 4ctes, 1dg, 2dg, 3dg, 2dg+3rd,
c     2dg+3rd+5rd, 3dg+5rd.
c
c
c     - input (RA,DEC) degrees
c     - output (RA,DEC) degrees
c     - errors and sigmas in arcsec (")
c     - coeficients (and errors) in radians
c
c
c     Update: M. Assafin  19/Apr/2010
c
c
c


      subroutine tanpla (rac,dec,id,nest,cora,code,corte,ngrau,ngrau3,
     ?ngrau5,nstart,nfinal,ucra,ucde,ranew,denew,alfsig,delsig,iall)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension id(iall),cora(iall),code(iall),xp(iall),yp(iall),
     ?xest(iall),yest(iall),ranew(iall),denew(iall)

      dimension ucra(iall),ucde(iall),era(iall),ede(iall),coefx(21),
     ?coefy(21),ecoefx(21),ecoefy(21),alfres(iall),delres(iall),
     ?itira(iall),xsao(21),ysao(21),xrray(21,21),yrray(21,21),
     ?array(21,21)



      COMMON /A7/ARRAY
      COMMON/A14/IERRO


      HMSGMS(I,J,A)=I+J/60.D0+A/3600.D0
      DEXY(XX,YY,Z,W)=DSIN(YY)*DSIN(W)+DCOS(YY)*DCOS(W)*DCOS(XX-Z)
      XPAD(XX,YY,Z)=DCOS(YY)*DSIN(XX-Z)
      YPAD(XX,YY,Z,W)=DSIN(YY)*DCOS(W)-DCOS(YY)*DSIN(W)*DCOS(XX-Z)
      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Dados auxiliares
C

      idiobs=iall   
      ncof=21

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
      IERRO=0
      DNOVE=9.999D0
      d99=99.99d0

      EQUIN=2000.d0
      IZERO=0
      ZERO=0.D0
      pma=0.d0
      pmd=0.d0

c
c     Zerando vetores
c

      do i=1,idiobs

      era(i)=0.d0
      ede(i)=0.d0

      xest(i)=0.d0
      yest(i)=0.d0
      xp(i)=0.d0
      yp(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0
      itira(i)=0

      enddo


      do i=1,ncof

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0
      xsao(i)=0.d0
      ysao(i)=0.d0

      do j=1,ncof
      xrray(j,i)=0.d0
      yrray(j,i)=0.d0
      array(j,i)=0.d0
      enddo

      enddo


c
c     Recolhendo estrelas medidas comuns ao catalogo de referencia
c




      nstart=0
      grac=grarad*rac
      gdec=grarad*dec

      do 10 i=1,nest

      if (id(i).le.0) go to 10

      nstart=nstart+1

      bra=grarad*15.d0*cora(i)
      bde=grarad*code(i)
      d=DEXY(bra,bde,grac,gdec)

      xest(nstart)=xpad(bra,bde,grac)/d
      yest(nstart)=ypad(bra,bde,grac,gdec)/d

c
c     Projecao do catalogo de referencia no plano tangente
c


      bra=grarad*15.d0*ucra(i)
      bde=grarad*ucde(i)
      d=DEXY(bra,bde,grac,gdec)

      xp(nstart)=xpad(bra,bde,grac)/d
      yp(nstart)=ypad(bra,bde,grac,gdec)/d


 10   continue


c
c     Ajuste do modelo polinomial entre (x,y) e (X,Y)
c



      call solutp (ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,yp,ntira,
     ?coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,corte,
     ?xrray,yrray,iall)


      if (ierro.eq.1) then

      write (*,*) 'L.S. solution crashed.'
      write (46,*) 'L.S. solution crashed.'

c
c     Zerando vetores no erro
c

      do i=1,idiobs

      ranew(i)=0.d0
      denew(i)=0.d0
      era(i)=0.d0
      ede(i)=0.d0

      xest(i)=0.d0
      yest(i)=0.d0
      xp(i)=0.d0
      yp(i)=0.d0
      alfres(i)=0.d0
      delres(i)=0.d0
      itira(i)=0

      enddo


      do i=1,ncof

      coefx(i)=0.d0
      coefy(i)=0.d0
      ecoefx(i)=0.d0
      ecoefy(i)=0.d0
      xsao(i)=0.d0
      ysao(i)=0.d0
      do j=1,ncof
      xrray(j,i)=0.d0
      yrray(j,i)=0.d0
      array(j,i)=0.d0
      enddo

      enddo

      return
      endif



c

      nfinal=nstart-ntira

c
c     Determinacao do alfa e delta observado de cada estrela do campo
c
c     (RA,DEC) guardados em graus

c
c     Calculo do erro padrao em alfa e delta para cada estrela de campo
c

      XVAM=0.D0
      YVAM=0.D0
      XVAS=0.D0
      YVAS=0.D0


C
C     NGRAU=0 --> 4 ctes
C
      nterms=1
      do i=1,ngrau
      nterms=nterms+i+1
      enddo
C
      igrau=ngrau+1
      if (ngrau.eq.0) igrau=2
C
      if (ngrau.eq.0) then
      nterms=4
c     ngrau3=0
c     ngrau5=0
      endif
C
      iterms=nterms
      if (ngrau3.eq.3) iterms=iterms+1
      if (ngrau5.eq.5) iterms=iterms+1
C



      do 130 i=1,nest

      bra=grarad*15.d0*cora(i)
      bde=grarad*code(i)
      d=DEXY(bra,bde,grac,gdec)

      x=xpad(bra,bde,grac)/d
      y=ypad(bra,bde,grac,gdec)/d


      ICONT=0
      POLX=0.D0
      POLY=0.D0
      DO 20 N=1,IGRAU
      DO 20 LL=1,N
      ICONT=ICONT+1
      K=N-LL
      POLX=POLX+COEFX(ICONT)*(X**K)*(Y**(LL-1))
   20 POLY=POLY+COEFY(ICONT)*(X**K)*(Y**(LL-1))
C
      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)
      ENDIF

      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF

      RANEW(I)=POLX
      DENEW(I)=POLY

C
C     Calcula o erro padrao em alfa para cada estrela
C
      ICONT=0
      if (ngrau.ne.0) then
      DO  30 N=1,IGRAU
      DO  30 LL=1,N
      K=N-LL
      ICONT=ICONT+1
   30 XSAO(ICONT)=(X**K)*(Y**(LL-1))
      IF (NGRAU3.EQ.3) THEN
       ICONT=ICONT+1
       XSAO(ICONT)=X*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
       ICONT=ICONT+1
       XSAO(ICONT)=X*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
c
      else
      xsao(1)=x
      xsao(2)=y
      xsao(3)=1.d0
      xsao(4)=0.d0
      endif
c
      DO  40 K=1,ITERMS
   40 YSAO(K)=0.D0
C
      DO  50 LL=1,ITERMS
      DO  50 K=1,ITERMS
   50 YSAO(LL)=YSAO(LL)+XRRAY(LL,K)*XSAO(K)

      DO  60 K=1,ITERMS
   60 era(I)=era(I)+XSAO(K)*YSAO(K)

      era(I)=alfsig*DSQRT(era(I))
      XVAM=XVAM+era(I)
      XVAS=XVAS+era(I)**2


C
C     Calcula o erro padrao em delta para cada estrela
C

      ICONT=0
      if (ngrau.ne.0) then
      DO  70 N=1,IGRAU
      DO  70 LL=1,N
      K=N-LL
      ICONT=ICONT+1
   70 YSAO(ICONT)=(X**K)*(Y**(LL-1))
      IF (NGRAU3.EQ.3) THEN
       ICONT=ICONT+1
       YSAO(ICONT)=Y*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
       ICONT=ICONT+1
       YSAO(ICONT)=Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
c
      else
      ysao(1)=y
      ysao(2)=-x
      ysao(3)=0.d0
      ysao(4)=1.d0
      endif
c
      DO  80 K=1,ITERMS
   80 XSAO(K)=0.D0
C
      DO  90 LL=1,ITERMS
      DO  90 K=1,ITERMS
   90 XSAO(LL)=XSAO(LL)+YRRAY(LL,K)*YSAO(K)

      DO 100 K=1,ITERMS
  100 ede(I)=ede(I)+YSAO(K)*XSAO(K)

      ede(I)=delsig*DSQRT(ede(I))
      YVAM=YVAM+ede(I)
      YVAS=YVAS+ede(I)**2


C
C     Calcula erro dos coeficientes para a solucao X
C

c     IF (NGRAU.NE.0) THEN
c     DO 110 K=1,ITERMS
c110  ECOEFX(K)=ALFSIG*DSQRT(XRRAY(K,K))
c     ELSE
c     ECOEFX(1)=ALFSIG*DSQRT(XRRAY(3,3))
c     ECOEFX(2)=ALFSIG*DSQRT(XRRAY(1,1))
c     ECOEFX(3)=ALFSIG*DSQRT(XRRAY(2,2))
c     ENDIF



C
C     Calcula erro dos coeficientes para a solucao Y
C

c     IF (NGRAU.NE.0) THEN
c     DO 120 K=1,ITERMS
c120  ECOEFY(K)=DELSIG*DSQRT(YRRAY(K,K))
c     ELSE
c     ECOEFY(1)=DELSIG*DSQRT(YRRAY(4,4))
c     ECOEFY(2)=DELSIG*DSQRT(YRRAY(2,2))
c     ECOEFY(3)=DELSIG*DSQRT(YRRAY(1,1))
c     ENDIF

c     do iii=1,nest
c     write (*,*) 'itira = ',itira(iii)
c     enddo
c     stop

C
C

 130  CONTINUE

C
C
      if (ngrau.eq.0) then
      iterms=3
      endif

C
C     Media dos erros padrao de todas as estrelas medidas
C

      EXMED=XVAM/NEST
      EYMED=YVAM/NEST
      XVAS=DSQRT((XVAS-2.D0*EXMED*XVAM+NEST*EXMED**2)/(NEST-1.D0))
      YVAS=DSQRT((YVAS-2.D0*EYMED*YVAM+NEST*EYMED**2)/(NEST-1.D0))
      AVAM=EXMED
      DVAM=EYMED

c
c     RA e DEC em graus
c


      j=0

      do 140 i=1,nest

      x=ranew(i)
      y=denew(i)

      ranew(i)=alff(x,y,grac,gdec)
      denew(i)=deltt(ranew(i),y,grac,gdec)

      ranew(i)=ranew(i)*radgra/15.d0
      denew(i)=denew(i)*radgra

c     if (id(i).eq.0) go to 140

c     j=j+1

c     write(*,*) 'alfres delres itira = ',j,i,alfres(j),delres(j),
c    ?itira(j) 



 140  continue


c
c     Debug
c


c     perc=100.d0*ntira/nstart
c     write (*,141) alfsig,delsig,nstart,nfinal,perc,avam,dvam
c141  format(1x,'alfsig delsig NI NF = ',2(1x,f6.3),2(1x,i4),1x,f6.2,
c    ?'%',2(1x,f6.3))


      return
      end


C
C     Subrotine solutp 
C
C
C     Fits bivariate polynom P=P(x,y) up to 3rd degree plus radial distortions
C     up to 5th degree 
C
c     Last update:  19/Apr/2010   M. Assafin
c
c


      subroutine solutp  (ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,yp,
     ?ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,
     ?corte,xrray,yrray,iall)


      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(21),COEFY(21),ALPHAX(21,21),ALPHAY(21,21),
     ?ARRAY(21,21),BETAX(21),BETAY(21),TERMX(21),TERMY(21),ITIRA(iall)

      DIMENSION XEST(iall),YEST(iall),XP(iall),YP(iall),
     ?XRRAY(21,21),YRRAY(21,21),ALFRES(iall),DELRES(iall)

      COMMON /A7/ARRAY
      COMMON/A14/IERRO


      ALFF (XX,YY,ZZ,WW)=ZZ+DATAN2(XX,DCOS(WW)-YY*DSIN(WW))
      DELTT(XX,YY,ZZ,WW)=DATAN2((YY*DCOS(WW)+DSIN(WW))*DCOS(XX-ZZ),
     ?DCOS(WW)-YY*DSIN(WW))

C
C     Initializing data
C
      PI=3.141592653589793D0
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
      DET  =1.D0
      IERRO=0

      ifora=0
c

      alfsig=0.d0
      delsig=0.d0

c

      sigma=grarad*corte/3600d0

C
      NCOMUM=nstart
C
c
c     Computa no. de termos do polinomio
c     NGRAU=0 --> 4 constantes
c

      NTERMS=1
      DO I=1,NGRAU
      NTERMS=NTERMS+I+1
      ENDDO

C
      IGRAU=NGRAU+1
      IF (NGRAU.EQ.0) IGRAU=2
      NTIRA=0
C


 5    LCONT=0

C
      RESMX=0.D0
      RES2X=0.D0
      RESMY=0.D0
      RES2Y=0.D0
C
      IF (NGRAU.EQ.0) THEN
      NTERMS=4
      NGRAU3=0
      NGRAU5=0
      ENDIF
C
      ITERMS=NTERMS
      IF (NGRAU3.EQ.3) ITERMS=ITERMS+1
      IF (NGRAU5.EQ.5) ITERMS=ITERMS+1

C
C     Checa No. de estrelas versus no. de coeficientes a ajustar
C

      IF (NGRAU.EQ.0) THEN
      iequa=2.d0*ncomum
      mtira=2.d0*ntira
      if ((iequa-mtira).LT.ITERMS) then
      IERRO=1
      RETURN
      endif
      else
      if ((ncomum-ntira).LT.ITERMS) then
      IERRO=1
      RETURN
      endif
      ENDIF

c
      DO 9 I=1,ITERMS
      BETAX(I) =0.D0
      BETAY(I) =0.D0
      COEFX(I) =0.D0
      COEFY(I) =0.D0
      TERMX(I) =0.D0
      TERMY(I) =0.D0
      DO 9 J=I,ITERMS
      ALPHAX(I,J)=0.D0
 9    ALPHAY(I,J)=0.D0


C
C     Montando equacoes de condicao para ajuste com 4 ctes
C

      IF (NGRAU.EQ.0) THEN

      DO 10 I=1,NCOMUM
      IF (ITIRA(I).NE.0) GO TO 10
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computa  termos dos coeficientes para AtB
C

      TERMX(1)=TERMX(1)+X**2+Y**2
      TERMX(2)=TERMX(2)+X
      TERMX(3)=TERMX(3)+Y
      TERMX(4)=TERMX(4)+1.D0
      BETAX(1)=BETAX(1)+X*XG+Y*YG
      BETAX(2)=BETAX(2)+Y*XG-X*YG
      BETAX(3)=BETAX(3)+XG
      BETAX(4)=BETAX(4)+YG
   10 CONTINUE

C
C     Preenchendo  AtA para ajuste com 4 ctes
C

      ALPHAX(1,1)=TERMX(1)
      ALPHAX(1,2)=0.D0
      ALPHAX(1,3)=TERMX(2)
      ALPHAX(1,4)=TERMX(3)
      ALPHAX(2,2)=TERMX(1)
      ALPHAX(2,3)=TERMX(3)
      ALPHAX(2,4)=-TERMX(2)
      ALPHAX(3,3)=TERMX(4)
      ALPHAX(3,4)=0.D0
      ALPHAX(4,4)=TERMX(4)
      DO 15 L=1,ITERMS
      DO 15 K=L,ITERMS
   15 ALPHAX(K,L) =ALPHAX(L,K)
C
      GO TO 75

      ENDIF


C
C     Montando equacoes de condicao para ajustes que NAO o de 4 ctes
C     para polinomio X


      DO 65 I=1,NCOMUM
      IF (ITIRA(I).NE.0) GO TO 65
      X=XEST(I)
      Y=YEST(I)
      XG=XP(I)
      YG=YP(I)
C
C     Computando termos para o polinomio em X
C
      ICONT=0
      DO 40 N=1,IGRAU
      DO 40 L=1,N
      K=N-L
      ICONT=ICONT+1
      TERMX(ICONT)=(X**K)*(Y**(L-1))
   40 TERMY(ICONT)=(X**K)*(Y**(L-1))

      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      TERMX(ICONT)=X*(X**2+Y**2)
      TERMY(ICONT)=Y*(X**2+Y**2)
      ENDIF

      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      TERMX(ICONT)=X*(X**2+Y**2)*(X**2+Y**2)
      TERMY(ICONT)=Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF
C
C     Computando AtA and AtB (para "X" e "Y")
C

      DO 60 L=1,ITERMS
      BETAX(L)=BETAX(L)+XG*TERMX(L)
      BETAY(L)=BETAY(L)+YG*TERMY(L)
      DO 60 K=L,ITERMS
      ALPHAX(L,K)=ALPHAX(L,K)+TERMX(K)*TERMX(L)
   60 ALPHAY(L,K)=ALPHAY(L,K)+TERMY(K)*TERMY(L)
C
   65 CONTINUE

C
C     Preenchendo parte triangular inferior da matriz simetrica AtA
c     (para "X" e "Y")
C

      DO 70 L=1,ITERMS
      DO 70 K=L,ITERMS
      ALPHAX(K,L)=ALPHAX(L,K)
   70 ALPHAY(K,L)=ALPHAY(L,K)

C
C     Preenchendo ARRAY=AtA para inversao (elementos normalizados pela
C     diagonal) para "X"
C

 75   CONTINUE
C
      DO 80 L=1,ITERMS
      DO 80 K=1,ITERMS
   80 ARRAY(L,K)=ALPHAX(L,K)/DSQRT(ALPHAX(L,L)*ALPHAX(K,K))
C
C     Invertendo AtA para "X"
C



      CALL MATINV (ITERMS,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando coeficientes do polinomio para "X"
C
      DO 90 L=1,ITERMS
      DO 90 K=1,ITERMS
   90 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHAX(K,K)*ALPHAX(L,L))
C
C     Backup do array para "X"
C
      IF (NGRAU.EQ.0) GO TO 97
      DO 95 L=1,ITERMS
      DO 95 K=1,ITERMS
   95 XRRAY(L,K)=ARRAY(L,K)
C
 97   DO 100 L=1,ITERMS
      DO 100 K=1,ITERMS
  100 COEFX(L)=COEFX(L)+ARRAY(L,K)*BETAX(K)
C
C     Obtem coeficientes e backup do array para modelo de 4 Constantes
C
      IF (NGRAU.EQ.0) THEN
      DO  105 L=1,ITERMS
 105  TERMX(L)=COEFX(L)
      COEFX(1)=TERMX(3)
      COEFX(2)=TERMX(1)
      COEFX(3)=TERMX(2)
      COEFY(1)=TERMX(4)
      COEFY(2)=-TERMX(2)
      COEFY(3)=TERMX(1)
      do 107 l=1,iterms
      do 107 k=1,iterms
      xrray(l,k)=array(l,k)
 107  yrray(l,k)=array(l,k)
      GO TO 133
      ENDIF

C
C     Preenchendo ARRAY=AtA para inversao (elementos normalizados pela
C     diagonal) para "Y"
C

      DO 110 L=1,ITERMS
      DO 110 K=1,ITERMS
  110 ARRAY(L,K)=ALPHAY(L,K)/DSQRT(ALPHAY(L,L)*ALPHAY(K,K))
C
C     Invertendo AtA para "Y"
C

      CALL MATINV (ITERMS,DET)
      IF (IERRO.EQ.1) RETURN
C
C     Computando coeficientes do polinomio para "Y"
C


      DO 120 L=1,ITERMS
      DO 120 K=1,ITERMS
  120 ARRAY(L,K)=ARRAY(L,K)/DSQRT(ALPHAY(K,K)*ALPHAY(L,L))

C
C     Backup do array para "Y"
C

      DO 125 L=1,ITERMS
      DO 125 K=1,ITERMS
  125 YRRAY(L,K)=ARRAY(L,K)
C
      DO 130 L=1,ITERMS
      DO 130 K=1,ITERMS
  130 COEFY(L)=COEFY(L)+ARRAY(L,K)*BETAY(K)

c
c
c     Computa residuos. Estrelas com residuo em alfa ou delta mais alto
c     sao eliminadas uma a uma, ate que nenhuma possua (O-C) maior
c     que o valor definido pela variavel "corte"
c

 133  CONTINUE
C


      REMAXI=-1.D14

      DO 160 I=1,NCOMUM
      X=XEST(I)
      Y=YEST(I)
      xg=xp(i)
      yg=yp(i)
      ICONT=0
      POLX=0.D0
      POLY=0.D0
      DO 140 N=1,IGRAU
      DO 140 L=1,N
      ICONT=ICONT+1
      K=N-L
      POLX=POLX+COEFX(ICONT)*(X**K)*(Y**(L-1))
  140 POLY=POLY+COEFY(ICONT)*(X**K)*(Y**(L-1))
C
      IF (NGRAU3.EQ.3) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)
      ENDIF
      IF (NGRAU5.EQ.5) THEN
      ICONT=ICONT+1
      POLX=POLX+COEFX(ICONT)*X*(X**2+Y**2)*(X**2+Y**2)
      POLY=POLY+COEFY(ICONT)*Y*(X**2+Y**2)*(X**2+Y**2)
      ENDIF

C
C     Computa (O-C) para R.A. e  Dec
C

      xx=alff(polx,poly,grac,gdec)
      yy=deltt(xx,poly,grac,gdec)

      xxr=alff(xg,yg,grac,gdec)
      yyr=deltt(xxr,yg,grac,gdec)


      AUX=(XX-XXR)*DCOS(YYR)
      AUY=YY-YYR
      ALFRES(I)=AUX*radgra*3600d0
      DELRES(I)=AUY*radgra*3600d0


C
      IF (ITIRA(I).NE.0) GO TO 160
      IF ((DABS(AUX).GT.REMAXI).OR.(DABS(AUY).GT.REMAXI)) THEN
      IFORA =I
      REMAXI=DMAX1(DABS(AUX),DABS(AUY))
      ENDIF

      RESMX=RESMX+ALFRES(I)
      RES2X=RES2X+ALFRES(I)**2
      RESMY=RESMY+DELRES(I)
      RES2Y=RES2Y+DELRES(I)**2

      LCONT=LCONT+1

C
  160 CONTINUE
C
C     Atingido numero minimo de estrelas!
C
      iwar=0
      IF (NGRAU.EQ.0) THEN
      iequa=2.d0*LCONT
      if (iequa.eq.4) iwar=1
      ELSE
      if (LCONT.eq.ITERMS) iwar=1
      ENDIF
C
      IF (IWAR.EQ.1) THEN
      alfsig=0.D0
      delsig=0.D0
      RETURN
      ENDIF


C
C     Computa media e desvio padrao dos (O-C)s
C

      XMED=RESMX/LCONT
      YMED=RESMY/LCONT
C
      IF (NGRAU.EQ.0) THEN
      alfsig=DSQRT((RES2X-2.D0*XMED*RESMX+LCONT*XMED**2)/(LCONT-2.D0))
      delsig=DSQRT((RES2Y-2.D0*YMED*RESMY+LCONT*YMED**2)/(LCONT-2.D0))
      ELSE
      alfsig=DSQRT((RES2X-2.D0*XMED*RESMX+LCONT*XMED**2)/(LCONT-ITERMS))
      delsig=DSQRT((RES2Y-2.D0*YMED*RESMY+LCONT*YMED**2)/(LCONT-ITERMS))
      ENDIF

c
c     Atingido o corte !
c

c     write (*,*) 'remaxi ',remaxi*radgra*3600.d0

      if (remaxi.le.sigma) go to 200

c
c     Corte nao atingido, prosseguir com eliminacao de estrelas
c

      if (ifora.ne.0) then

      NTIRA=NTIRA+1
      ITIRA(IFORA)=1
      GO TO 5

      endif

 200  continue

      RETURN
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



      END





c
c     Subrotine transf
c
c
c     Rigid rotation or tangent plane technique applied to transform
c     positions toward the UCAC2 reference frame.
c 
c     The position variables "cora" and "code" are mantained; the
c     new transformed positions are output to variables "oldra, oldde".
c     
c     ktrans = 1  -> rigid rotation
c
c     ktrans = 2  -> tangent plane techinique
c
c
c
c         Last update: M. Assafin   19/Apr/2010
c

      subroutine transf (ktrans,cortu,nobs,cora,code,ucra,ucde,ranew,
     ?denew,id,rac,dec,jgrau,jgrau3,jgrau5,rmx,rmy,rmx2,rmy2,iall)


      IMPLICIT REAL *8 (A-H,O-Z)

      dimension cora(iall),code(iall),id(iall),
     ?ranew(iall),denew(iall),ucra(iall),ucde(iall)




c
C     Initializing data
C
c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI


c
c     Here we have rigid rotation
c


      if (ktrans.eq.1) then


      write (*,*)
      write (*,*) 'Rotating positions toward the UCAC2 reference frame'
      write (*,*)

      write (46,*)
      write (46,*) 'Rotating positions toward the UCAC2 reference frame'
      write (46,*)


      call rotsky (cortu,nobs,cora,code,ucra,ucde,ranew,denew,id,a1,a2,
     ?a3,e1,e2,e3,rmx,rmy,rmx2,rmy2,rmgx,rmgy,rmgx2,rmgy2,nc,iall)

      write (*,*)
      write (46,*)

c
      r2=dsqrt(rmgx2**2+rmgy2**2)

c
      write (*,*)
      write (*,*)'Number of used UCAC2 stars for rotation: ',nc
      write (*,*)

      write (*,*)'Mean offsets (RA,DEC) before rotation (arcsec): ',rmx,
     ?rmy
      write (*,*)
      write (*,*)'Sigma offsets (RA,DEC) before rotation (arcsec): ',
     ?rmx2,rmy2
      write (*,*)

      write (*,*)'Mean offsets (RA,DEC) after rotation (arcsec): ',rmgx,
     ?rmgy
      write (*,*)
      write (*,*)'Sigma offsets (RA,DEC) after rotation (arcsec): ',
     ?rmgx2,rmgy2
      write (*,*)

      write (*,*) 'Mean error of rotation adjustment (arcsec): ',r2
      write (*,*)

      write (*,*) 'Rotating angles and errors (arcsec):'
      write (*,*)

      write (46,*)
      write (46,*)'Number of used UCAC2 stars for rotation: ',nc
      write (46,*)

      write(46,*)'Mean offsets (RA,DEC) before rotation (arcsec): ',rmx,
     ?rmy
      write (46,*)
      write (46,*)'Sigma offsets (RA,DEC) before rotation (arcsec): ',
     ?rmx2,rmy2
      write (46,*)

      write(46,*)'Mean offsets (RA,DEC) after rotation (arcsec): ',rmgx,
     ?rmgy
      write (46,*)
      write (46,*)'Sigma offsets (RA,DEC) after rotation (arcsec): ',
     ?rmgx2,rmgy2
      write (46,*)

      write (46,*) 'Mean error of rotation adjustment (arcsec): ',r2
      write (46,*)

      write (46,*) 'Rotating angles and errors (arcsec):'
      write (46,*)


      aa1=a1*radgra*3600.d0
      aa2=a2*radgra*3600.d0
      aa3=a3*radgra*3600.d0

      e1=e1*radgra*3600.d0
      e2=e2*radgra*3600.d0
      e3=e3*radgra*3600.d0

      write (*,*) 'X  a1 = ',aa1
      write (*,*) 'Y  a2 = ',aa2
      write (*,*) 'Z  a3 = ',aa3
      write (*,*)
      write (*,*) 'X  e1 = ',e1
      write (*,*) 'Y  e2 = ',e2
      write (*,*) 'Z  e3 = ',e3

      write (*,*)


      write (46,*) 'X  a1 = ',aa1
      write (46,*) 'Y  a2 = ',aa2
      write (46,*) 'Z  a3 = ',aa3
      write (46,*)
      write (46,*) 'X  e1 = ',e1
      write (46,*) 'Y  e2 = ',e2
      write (46,*) 'Z  e3 = ',e3

      write (46,*)

c
c     record last statistics
c

      rmx=rmgx
      rmy=rmgy
      rmx2=rmgx2
      rmy2=rmgy2

      endif

c
c     Here we have tangent plane techinique
c


      if (ktrans.eq.2) then


      write (*,*)
      write (*,*) 'Transform positions toward the UCAC2 reference frame'
      write (*,*)

      write (46,*)
      write (46,*)'Transform positions toward the UCAC2 reference frame'
      write (46,*)


      write (*,*)
      write (46,*)

c

      call tanpla (rac,dec,id,nobs,cora,code,cortu,jgrau,jgrau3,
     ?jgrau5,nstart,nfinal,ucra,ucde,ranew,denew,rmgx2,rmgy2,iall)


      write (*,*)
      write (*,*)'Initial number of UCAC2 stars: ',nstart
      write (*,*)'Number of used UCAC2 stars for transform: ',nfinal
      write (*,*)

      write (*,*)'Mean offsets (RA,DEC) before transform(arcsec): ',rmx,
     ?rmy
      write (*,*)
      write (*,*)'Sigma offsets (RA,DEC) before transform(arcsec): ',
     ?rmx2,rmy2
      write (*,*)

      write (*,*)
      write (*,*)'Sigma offsets (RA,DEC) after transform (arcsec): ',
     ?rmgx2,rmgy2
      write (*,*)


      write (46,*)
      write (46,*)'Initial number of UCAC2 stars: ',nstart
      write (46,*)'Number of used UCAC2 stars for transform: ',nfinal
      write (46,*)

      write(46,*)'Mean offsets (RA,DEC) before transform(arcsec): ',rmx,
     ?rmy
      write (46,*)
      write (46,*)'Sigma offsets (RA,DEC) before transform(arcsec): ',
     ?rmx2,rmy2
      write (46,*)

      write (46,*)
      write (46,*)'Sigma offsets (RA,DEC) after transform (arcsec): ',
     ?rmgx2,rmgy2
      write (46,*)


c
c     record last statistics
c

      rmx=0.d0
      rmy=0.d0
      rmx2=rmgx2
      rmy2=rmgy2


      endif


      return
      end



