c
c     Program PRAIA_uniform_astrometry_10
c
c     Purpose
c
c
c     Given fields reduced with PRAIA, re-does the reductions of these fields
c     using the same reference star set for all of them.
C
c
c      Last update: Marcelo Assafin - 05/July/2013
c   
c
c

      IMPLICIT REAL *8 (A-H,O-Z)
      parameter(idiobs=150000,icof=21,kest=50,ina=50,nfraxy=100000)

      character*50 listxy,stars,ifxy,ifxyo,kfxy(nfraxy),iver,iva

      character*(kest) iext


c
C     Initializing data
C
c


      stars='PRAIA_uniform_astrometry.tmp'

c
c     Reading input batch file
c

      write (*,*) 
      write (*,*) 
      write (*,*) 'PRAIA uniform reduction'
      write (*,*) 
      write (*,*) 
      write (*,*) 

c


c     open (1,file='PRAIA_uniform_astrometry_20_02.dat')

 1    format(a50)


      read (*,1) listxy
      read (*,3) iext
 3    format(a50)

      read (*,*) fac

      fac=fac/100.d0

      read (*,*) cut

      read (*,*) ngrau
      read (*,*) ngrau3
      read (*,*) ngrau5


c
c     Selects common reference stars throughout xy fields
c

      write (*,*) 
      write(*,*) 'Selecting common reference stars to all reductions '
      write (*,*) 

      call select (idiobs,listxy,cut,fac,stars)


c
c     Re-reduction of xy fields
c

      write(*,*) 'Re-reduction of xy fields'
      write (*,*) 

c
c     Finds common xy extension string of input files for replacement
c     by new ".xy" string furnished by the user 
c


      open (7,file=listxy)

      do i=1,nfraxy
      read (7,1,end=4) kfxy(i)
      enddo

 4    close (7)
      kfra=i-1

      iver=''
      iver=kfxy(1)


      do i=ina,1,-1
      if (iver(i:i).ne.' ') go to 5
      enddo

 5    iver2=i


      do m=iver2,1,-1

      mm=iver2-m

      do k=2,kfra
      iva=''
      iva=kfxy(k)

      do j=ina,1,-1
      if (iva(j:j).ne.' ') go to 6
      enddo

 6    iva2=j


      if (iver(m:iver2).ne.iva(iva2-mm:iva2)) go to 7

      enddo

      enddo

 7    m=mm

c
      do i=kest,1,-1
      if (iext(i:i).ne.' ') go to 8
      enddo

 8    mext=i


c


      open (7,file=listxy)

 10   continue

      read (7,1,end=20) ifxy


      ifxyo=''

      do i=ina,1,-1
      if (ifxy(i:i).ne.' ') go to 15
      enddo

 15   k=i-m+1
      ifxyo=ifxy
      ifxyo(k:k)='.'
      ifxyo(k+1:k+mext)=iext

      do i=k+mext+1,ina
      ifxyo(i:i)=' '
      enddo


      call nreduc (idiobs,icof,ifxy,stars,ifxyo,ngrau,ngrau3,ngrau5)

      go to 10

c

 20   close (7)

c

      call system ('rm PRAIA_uniform_astrometry.tmp')

      write (*,*) 
      write (*,*) 
      write (*,*) 'Execution terminated succesfully.'
      write (*,*) 
      write (*,*) 
      write (*,*) 

      end





c
c     
c     Subroutine select
c
c     Selects common reference stars from xy files from former reductions
c     with PRAIA.
c 
c
c     Last update: M. Assafin : 05/July/2013
c

      subroutine select (idiobs,listxy,cut,fac,stars)

      IMPLICIT REAL *8 (A-H,O-Z)

      dimension rap(idiobs),dep(idiobs),marca(idiobs)


      character*50 listxy,stars
      character*50 infits
      character*50 ifxy
      character*20 iobalv,ichfil


c
c     Auxiliary data
c

      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI


c
c     Emptying vectors
c

      do i=1,idiobs
      rap(i)=0.d0
      dep(i)=0.d0
      marca(i)=0
      enddo

c
c     Error box for common stars (arcsec)
c

      ebox=1.d0

c
c     Open files
c 


      open (1,file=listxy)

c
c     Number of frames?
c

      do i=1,10000000
      read (1,*,end=1)
      enddo

 1    rewind (1)

      nf=i-1

      nframe=nf*fac

c
c     Reads name of 1rst xy frame
c


      read(1,2) ifxy
 2    format(a50)

c
c    Takes stars from first frame as basis for search of common stars
c
c    Only reference stars within (O-C) < cut are taken
c

      open (3,file=ifxy)


      i=0

      do 15 j=1,idiobs

      read (3,10,end=20) xob,yob,seng,altu,fgcc,
     ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,kuth,kutm,zut,kutano,
     ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny

 10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))


c
c     read (3,10,end=20) xob,yob,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
c    ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,kuth,kutm,zut,kutano,
c    ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny
c
c
c10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),3(1x,f6.3),14x,8(1x,f6.3),
c    ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
c    ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
c    ?1x,a20,2(1x,i5))
c


      if (ktirau.eq.99) go to 15

      if (dabs(alsiu).gt.cut) go to 15
      if (dabs(desiu).gt.cut) go to 15


      i=i+1
      rap(i)=ra
      dep(i)=de

 15   continue

 20   close(3)

      nest=i

c
c     Searches all frames for common reference stars
c

      rewind (1)

 30   continue

      read (1,2,end=60) ifxy


c
c     Checks for the presence of star in nth xy frame
c

      open (3,file=ifxy)

      do 40 j=1,idiobs

      read (3,10,end=50) xob,yob,seng,altu,fgcc,
     ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,kuth,kutm,zut,kutano,
     ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny


c     read (3,10,end=50) xob,yob,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
c    ?nstaru,nfinau,alsiu,desiu,ktirau,ra,de,kuth,kutm,zut,kutano,
c    ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny


      if (ktirau.eq.99) go to 40

      if (dabs(alsiu).gt.cut) go to 40
      if (dabs(desiu).gt.cut) go to 40

      dist=1.d14

      do k=1,nest

      dx=(rap(k)-ra)*dcos(de*grarad)*15.d0*3600.d0
      dy=(dep(k)-de)*3600.d0
      dis=dsqrt(dx**2+dy**2)

      if (dis.lt.dist) then
      dist=dis
      kk=k
      endif

      enddo

      if (dist.lt.ebox) marca(kk)=marca(kk)+1

 40   continue

 50   close (3)

c
      go to 30
c

 60   close (1)

c
c     Writes surviving common reference stars (RA,DEC)s
c

      open (2,file=stars)


      nnn=0      

      do i=1,nest

      if (marca(i).ge.nframe) then
      write (2,*) rap(i),dep(i)
      nnn=nnn+1
      endif
      enddo

      close (2)

      write (*,*) 
      write (*,*) 'Common reference stars = ',nnn
      write (*,*) 'Total number of frames = ',nf
      write (*,*) 'Least number of frames = ',nframe
      write (*,*) 

      return
      end



c
c
c     
c     Subroutine nreduc
c
c     Re-reduces (RA,DEC) using selected common reference stars from xy
c     files from former reductions with PRAIA.
c 
c
c     Last update: M. Assafin : 05/July/2013
c


      subroutine nreduc (idiobs,icof,ifxy,stars,ifxyo,ngrau,ngrau3,
     ?ngrau5)


      IMPLICIT REAL *8 (A-H,O-Z)

      dimension cora(idiobs),code(idiobs)
      dimension id(idiobs),racat(idiobs),decat(idiobs),xob(idiobs),
     ?yob(idiobs)

      dimension ra(idiobs),de(idiobs),era(idiobs),ede(idiobs),
     ?coefx(icof),coefy(icof),ecoefx(icof),ecoefy(icof),alfres(idiobs),
     ?delres(idiobs),itira(idiobs)

      character*50 ifxy,stars,ifxyo

      character*50 infits
      character*20 ichfil,iobalv
      character*1  isig,menos

      hmsgms(i,j,a)=i+j/60.d0+a/3600.d0


c
C     Auxiliary data
C
c
      menos='-'
      PI    = 0.3141592653589793D1
      GRARAD=PI/180.D0
      RADGRA=180.D0/PI
c


      ireflex=1
      corte=10.d0
      ebox=1.d0


c
c     Retrieves common star approx. coordinates
c

      open (77,file=stars)

      do i=1,idiobs
      read (77,*,end=5) cora(i),code(i)
      enddo

 5    close (77)
      ncat=i-1



c
c     Computes approximate central (RA,DEC) of fields
c
c     rac = RA  center in degrees
c     dec = DEC center in degrees 
c     

      rac=0.d0
      dec=0.d0

      do i=1,ncat
      rac=rac+cora(i)
      dec=dec+code(i)
      enddo

      rac=rac/ncat
      dec=dec/ncat

      rac=rac*15.d0

c
c     Retrieving data from older field reduction
c


      open (40,file=ifxy)

      i=0

      m=0
     
      do 15 ii=1,idiobs

      read (40,10,end=20) xo,yo,seng,altu,fgcc,
     ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,raest,deest,kuth,kutm,zut,kutano,
     ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny


 10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?2(1x,i5))

c     read (40,10,end=20) xo,yo,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
c    ?nstaru,nfinau,alsiu,desiu,ktirau,raest,deest,kuth,kutm,zut,kutano,
c    ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny
c
c10   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),3(1x,f6.3),14x,8(1x,f6.3),
c    ?4(1x,f7.3),6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,
c    ?i2,1x,i2,1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,
c    ?1x,a20,2(1x,i5))
c
c
c     cudmg=20.1
c     cudmg2=20.3


      i=i+1

      xob(i)=xo
      yob(i)=yo

      id(i)=0


      do j=1,ncat

      dx=(raest-cora(j))*dcos(deest*grarad)*15.d0*3600.d0
      dy=(deest-code(j))*3600.d0
      dis=dsqrt(dx**2+dy**2)

      if (dis.lt.ebox) then
      m=m+1
      id(i)=m
      racat(m)=raest*15.d0-(alsiu/3600.d0)/dabs(dcos(grarad*deest))
      decat(m)=deest-desiu/3600.d0
      endif

      enddo

 15   continue
c

 20   close (40)

      nest=i



c
c     RA,DEC Reduction
c


      call posred (idiobs,icof,ireflex,rac,dec,id,ncat,racat,decat,nest,
     ?xob,yob,corte,ngrau,ngrau3,ngrau5,nstart,nfinal,ra,de,era,ede,
     ?alfsig,delsig,alfres,delres,coefx,coefy,ecoefx,ecoefy,itira,avam,
     ?dvam)


c
c     Writes results of new reduction in new xy output file
c

      write (*,30) ifxy,alfsig,delsig,nfinal
 30   format(1x,a50,2(1x,f7.3),1x,i5)

      open (40,file=ifxy)
      open (41,file=ifxyo)

      i=0

      ncom=0
      aux=0.d0
      auy=0.d0

      do 35 ii=1,idiobs

      read (40,10,end=40) x,y,seng,altu,fgcc,
     ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
     ?nstaru,nfinau,alsiu,desiu,ktirau,raest,deest,kuth,kutm,zut,kutano,
     ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny

c     read (40,10,end=40) x,y,seng,altu,fgcc,
c    ?fumag,fumag2,xmgu,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
c    ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,erau,edeu,alfsiu,delsiu,
c    ?nstaru,nfinau,alsiu,desiu,ktirau,raest,deest,kuth,kutm,zut,kutano,
c    ?kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny
c
c
c     cudmg=20.1
c     cudmg2=20.3
c


      i=i+1

      resa=alsiu
      desa=desiu

      if (dabs(resa).gt.9.999d0) resa=99.999d0
      if (dabs(desa).gt.9.999d0) desa=99.999d0

      ktira=ktirau

c     if (ktirau.lt.2 .and. id(i).eq.0) then
      if (ktirau.lt.2) then

      raest=raest*15.d0-(alsiu/3600.d0)/dabs(dcos(deest*grarad))
      deest=deest-desiu/3600.d0

      resa=(ra(i)-raest)*dcos(de(i)*grarad)*3600.d0
      desa=(de(i)-deest)*3600.d0

      if (dabs(resa).gt.9.999d0) resa=99.999d0
      if (dabs(desa).gt.9.999d0) desa=99.999d0

      ktira=1

      endif

c

      if (id(i).ne.0) then

      ktira=0
      resa=alfres(id(i))
      desa=delres(id(i))

      endif

      ra(i)=ra(i)/15.d0

      write (41,33) xob(i),yob(i),seng,altu,fgcc,
     ?fumag,fumag2,xmgu,cudmg,cudmg2,xmgj,xmgh,xmgk,res2mg,resmg2,ermgj,
     ?ermgh,ermgk,pma,pmd,epma,epmd,ex,ey,era(i),ede(i),alfsig,delsig,
     ?nstart,nfinal,resa,desa,ktira,ra(i),de(i),kuth,kutm,zut,
     ?kutano,kutmes,kutdia,daju,iexps,ichfil,infits,iobalv,nx,ny,ncom,
     ?aux,auy

 33   format(2(1x,f7.2),1x,f5.3,2(1x,f9.2),13(1x,f6.3),4(1x,f7.3),
     ?6(1x,f6.3),2(1x,i4),2(1x,f6.3),1x,i4,2(1x,f13.9),1x,i2,1x,i2,
     ?1x,f5.2,1x,i4,1x,i2,1x,i2,1x,f16.8,2x,i4,2x,a20,2x,a50,1x,a20,
     ?3(1x,i5),2(1x,f7.3))


 35   continue

 40   close (40)
      close (41)

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
c     Atualizacao: M. Assafin  05/July/2013
c
c
c



      subroutine posred (idiobs,icof,ireflex,rac,dec,id,ncat,racat,
     ?decat,nest,xob,yob,corte,ngrau,ngrau3,ngrau5,nstart,nfinal,ra,de,
     ?era,ede,alfsig,delsig,alfres,delres,coefx,coefy,ecoefx,ecoefy,
     ?itira,avam,dvam)


      IMPLICIT REAL*8 (A-H,O-Z)

      dimension id(idiobs),racat(idiobs),decat(idiobs),xob(idiobs),
     ?yob(idiobs),xp(idiobs),yp(idiobs),xest(idiobs),yest(idiobs)

      dimension ra(idiobs),de(idiobs),era(idiobs),ede(idiobs),
     ?coefx(icof),coefy(icof),ecoefx(icof),ecoefy(icof),alfres(idiobs),
     ?delres(idiobs),itira(idiobs),xsao(icof),ysao(icof),
     ?xrray(icof,icof),yrray(icof,icof)


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

      ncof=icof

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


      call solucao (idiobs,icof,ngrau,ngrau3,ngrau5,nstart,xest,yest,xp,
     ?yp,ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,delres,itira,
     ?corte,xrray,yrray,ierro)


      if (ierro.eq.1) then

      write (*,*) 'L.S. solution crashed.'

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

      IF (NGRAU.NE.0) THEN
      DO 110 K=1,ITERMS
 110  ECOEFX(K)=ALFSIG*DSQRT(XRRAY(K,K))
      ELSE
      ECOEFX(1)=ALFSIG*DSQRT(XRRAY(3,3))
      ECOEFX(2)=ALFSIG*DSQRT(XRRAY(1,1))
      ECOEFX(3)=ALFSIG*DSQRT(XRRAY(2,2))
      ENDIF
C
C     Calcula erro dos coeficientes para a solucao Y
C

      IF (NGRAU.NE.0) THEN
      DO 120 K=1,ITERMS
 120  ECOEFY(K)=DELSIG*DSQRT(YRRAY(K,K))
      ELSE
      ECOEFY(1)=DELSIG*DSQRT(YRRAY(4,4))
      ECOEFY(2)=DELSIG*DSQRT(YRRAY(2,2))
      ECOEFY(3)=DELSIG*DSQRT(YRRAY(1,1))
      ENDIF
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

c     write(*,*) 'alfres delres itira = ',j,alfres(j),delres(j),itira(j) 



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
C
C							    

      subroutine solucao (idiobs,icof,ngrau,ngrau3,ngrau5,nstart,xest,
     ?yest,xp,yp,ntira,coefx,coefy,alfsig,delsig,grac,gdec,alfres,
     ?delres,itira,corte,xrray,yrray,ierro)




      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION COEFX(icof),COEFY(icof),ALPHAX(icof,icof),
     ?ALPHAY(icof,icof),ARRAY(icof,icof),BETAX(icof),BETAY(icof),
     ?TERMX(icof),TERMY(icof),ITIRA(idiobs)

      DIMENSION XEST(idiobs),YEST(idiobs),XP(idiobs),YP(idiobs),
     ?XRRAY(icof,icof),YRRAY(icof,icof),ALFRES(idiobs),DELRES(idiobs)


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



      CALL MATINV (icof,array,ITERMS,DET,ierro)
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

      CALL MATINV (icof,array,ITERMS,DET,ierro)
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

 157  RESMX=RESMX+ALFRES(I)
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

      if (remaxi.lt.sigma) return

c
c     Corte nao atingido, prosseguir com eliminacao de estrelas
c

      NTIRA=NTIRA+1
      ITIRA(IFORA)=1
      GO TO 5


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
      SUBROUTINE MATINV (icof,array,NORDER, DET,ierro)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION ARRAY (icof,icof), IK(icof), JK(icof)
C
   10 DET = 1.D0
   11 DO 100 K=1, NORDER
C
C        FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
C
      AMAX= 0.D0
   21 DO 30 I=K, NORDER
      DO 30 J=K, NORDER
   23 IF (DABS(AMAX) - DABS(ARRAY(I,J))) 24, 24, 30
   24 AMAX = ARRAY(I,J)
      IK(K) = I
      JK(K) = J
   30 CONTINUE
C
C        INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)
C
   31 IF (AMAX) 41, 32, 41
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
   71 DO 80 I=1, NORDER
      DO 80 J=1, NORDER
      IF (I-K) 74, 80, 74
   74 IF (J-K) 75, 80, 75
   75 ARRAY(I,J) = ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
   80 CONTINUE
   81 DO 90 J=1, NORDER
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
  101 DO 130 L=1, NORDER
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



