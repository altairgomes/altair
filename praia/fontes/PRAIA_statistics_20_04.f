c
c
c     PRAIA_statistics
c
c
c     Computes basic statstics, i.e. averages, standard deviations, root-mean-square (rms)
c     and standard errors,  with and without weighting the data.
c
c     Weight is extracted from the data, as indicated by the user.
c
c
c     The dispersion of the data is computed by the standard deviation about the mean
c     (weighted and not weighted), and by the root-mean-square rms (weighted and not
c     weighted). We also compute the unbiased standard deviation (weighted and not weighted)
c      using the gama function. 
c
c     The standard error (standard deviation devided by the root-square of N) is also
c     computed (weighted and not weighted).
c
c
c     If no weight is available in the data, the user indicates columns less or equal
c     to zero.
c     
c
c     Last modification: M. Assafin  19/Ago/2012
c
c


      IMPLICIT REAL *8 (A-H,O-Z)
      parameter (stdin=5,stdout=6,icolm=500,idim=20000,idima=210000)


      dimension x(idima),y(idima),wx(idima),wy(idima)
      character *150 name,saida
      character*1 iv(idim),ibr
      character*10 iform
      data ibr/' '/


c
c     Initializing data
c

      iform='(        )'
      write (iform(2:7),'(i6.6)') idim

      iform(8:9)='a1'


c

      do i=1,idima
      x(i)=0.d0
      y(i)=0.d0
      wx(i)=1.d0
      wy(i)=1.d0
      enddo

c


      read (*,1) name
      read (*,1) saida
 1    format (50x,a150)

      read (*,2) xcutmi
      read (*,2) xcutma
      read (*,2) ycutmi
      read (*,2) ycutma
 2    format(50x,f16.0)

      read (*,3) ix
      read (*,3) iy

      read (*,3) isx
      read (*,3) isy
 3    format(50x,i10)


c

      ico=max0(ix,iy,isx,isy)

 
c

      if (ico.gt.icolm) then
      write (*,*) 'More columns than I can support, sorry.'
      write (*,*) 'Exiting ...'
      stop
      endif


c
c     Stores data
c

      open (7,file=name)


      do 7  i=1,idima

c     do 90 ii=1,idim
c90   iv(ii)=ibr

      read (7,iform,end=8) (iv(k),k=1,idim)
c4    format(iform)



      call chanum (ix,idim,iv,xxx,ist)
      if (ist.eq.0) then
      write (*,*) 'Error in retrieving values, exiting ...'
      close (7)
      stop
      endif

      x(i)=xxx

      call chanum (iy,idim,iv,yyy,ist)
      if (ist.eq.0) then
      write (*,*) 'Error in retrieving values, exiting ...'
      close (7)
      stop
      endif

      y(i)=yyy

      if (isx.le.0) go to 7
      if (isy.le.0) go to 7

      call chanum (isx,idim,iv,wxx,ist)
      if (ist.eq.0) then
      write (*,*) 'Error in retrieving values, exiting ...'
      close (7)
      stop
      endif

      wxx=dabs(wxx)

      wx(i)=1.d0/wxx


      call chanum (isy,idim,iv,wyy,ist)
      if (ist.eq.0) then
      write (*,*) 'Error in retrieving values, exiting ...'
      close (7)
      stop
      endif

      wyy=dabs(wyy)


      wy(i)=1.d0/wyy



 7    continue


 8    nfonts=i-1
      close (7)



c
c     Statistics
c


      k=0

      wkx=0.d0
      wky=0.d0

      wkx2=0.d0
      wky2=0.d0


      xm=0.d0
      ym=0.d0

      xmw=0.d0
      ymw=0.d0

      sx=0.d0
      sy=0.d0

      sxg=0.d0
      syg=0.d0

      sxw=0.d0
      syw=0.d0


      xrms=0.d0
      yrms=0.d0

      xrmsw=0.d0
      yrmsw=0.d0

      sex=0.d0
      sey=0.d0

      sexg=0.d0
      seyg=0.d0

      sexw=0.d0
      seyw=0.d0



c
c     Computes the averages
c


      do 10 i=1,nfonts

      if (x(i).lt.xcutmi) go to 10
      if (x(i).gt.xcutma) go to 10
      if (y(i).lt.ycutmi) go to 10
      if (y(i).gt.ycutma) go to 10


      k=k+1

      wkx=wkx+wx(i)
      wky=wky+wy(i)

      wkx2=wkx2+wx(i)**2
      wky2=wky2+wy(i)**2


      xm=xm+x(i)
      ym=ym+y(i)

      xmw=xmw+x(i)*wx(i)
      ymw=ymw+y(i)*wy(i)


 10   continue



c

      if (k.le.1) then
      write (*,*) 'Only 1 point retrieved. Exiting ...'
      stop
      endif

c



      dk=k

      xm=xm/dk
      ym=ym/dk

      xmw=xmw/wkx
      ymw=ymw/wky



c
c     Computes the dispersions: standard deviations, standard errors, rms
c



      do 20 i=1,nfonts

      if (x(i).lt.xcutmi) go to 20
      if (x(i).gt.xcutma) go to 20
      if (y(i).lt.ycutmi) go to 20
      if (y(i).gt.ycutma) go to 20


      xrms=xrms+x(i)**2
      yrms=yrms+y(i)**2

      xrmsw=xrmsw+(wx(i)*x(i))**2
      yrmsw=yrmsw+(wy(i)*y(i))**2


      sx=sx+(x(i)-xm)**2
      sy=sy+(y(i)-ym)**2


      sxw=sxw+(wx(i)*(x(i)-xm))**2
      syw=syw+(wy(i)*(y(i)-ym))**2



 20   continue

c

      cn=(2.d0*dexp(gammln(dk/2.d0)))/dexp(gammln((dk-1.d0)/2.d0))
      cn=cn**2


      sxg=dsqrt(sx/cn)
      syg=dsqrt(sy/cn)

      sx=dsqrt(sx/dk)
      sy=dsqrt(sy/dk)

      sxw=dsqrt(wkx*sxw/(wkx**2-wkx2))
      syw=dsqrt(wky*syw/(wky**2-wky2))

      xrms=dsqrt(xrms/dk)
      yrms=dsqrt(yrms/dk)

      xrmsw=dsqrt(xrmsw/wkx)
      yrmsw=dsqrt(yrmsw/wky)

      sex=sx/dsqrt(dk)
      sey=sy/dsqrt(dk)

      sexg=sxg/dsqrt(cn)
      seyg=syg/dsqrt(cn)

      sexu=sxg/dsqrt(dk)
      seyu=syg/dsqrt(dk)


      sexw=sxw/dsqrt(wkx)
      seyw=syw/dsqrt(wky)


      dfon=nfonts

c
c     Outputs statistics
c

      do 40 m=6,7

      if (m.eq.7) then
      open (m,file=saida)
      endif

c
      write(m,*) 
      write(m,*) 
      write(m,*)'Results: no   weighted data'
      write(m,*) 
      write(m,*)dfon,dk,      ' N_beggining, N_end'
      write(m,*)xm,ym,        ' averages xm, ym'
      write(m,*)sx,sy,        ' standard deviation sx, sy (n-1)'
      write(m,*)sxg,syg,      ' unbiased s.d.  sx, sy (gamma function)'
      write(m,*)xrms,yrms,    ' root-mean-square rms x, y'                        
      write(m,*)sex,sey,      ' standard error x, y'             
      write(m,*)sexu,seyu,    ' standard error x, y unbiased s.d., N'             
      write(m,*)sexg,seyg,    ' standard error x, y (gamma function)'             
      write(m,*) 
      write(m,*) 


      write(m,*) 
      write(m,*) 
      write(m,*)'Results: with weighted data'
      write(m,*) 
      write(m,*)dfon,dk,      ' N_beggining, N_end'
      write(m,*)xmw,ymw,      ' weighted averages xm, ym'
      write(m,*)sxw,syw,      ' weighted standard deviation sx, sy'
      write(m,*)xrmsw,yrmsw,  ' weighted root-mean-square rms x, y'
      write(m,*)sexw,seyw,    ' weighted standard error x, y'             
      write(m,*) 
      write(m,*) 



      write (m,*) 
      write (m,*) 
      write (m,*) 

      if (m.eq.7) then
      close (m)
      endif



 40   continue
c



      end



c
c     subrotina chanum
c
c     Pega uma string (character) de numeros e extrai o numero da
c     coluna, sem abrir arquivos temporarios para isso.
c
c     string  = contem a string completa
c     palavra = variavel de trabalho contendo a string completa
c     b2 = contem o numero extraido correspondente a coluna dada
c     ist = status das operacoes (0 = erro; maior que zero = ok)
c
c     Ultima modificacao: M. Assafin  28/Outubro/2010
c
c
c
      subroutine chanum (icol,id,string,valor,ist)

      implicit real *8 (a-h,o-z)

      integer*8 n

      dimension ni(id+2),nf(id+2)

      character*1 string(id),palavra(id+2),ibra

c

      ibra=' '

      ist=0


c

      do i=1,id+2
      palavra(i)=ibra
      enddo

      do i=1,id
      palavra(i+1)=string(i)
      enddo

c
c     Checando colunas pelos espacos em branco
c

      do i=1,id
      ni(id)=0
      nf(id)=0
      enddo

      ki=0
      kf=0

c
c     Onde estah o comeco do numero
c

      do i=2,id+2

      if (palavra(i-1).eq.ibra .and. palavra(i).ne.ibra) then
      ki=ki+1
      ni(ki)=i
      endif

      enddo

c
c     Onde estah o fim do numero
c

      do i=2,id+2

      if (palavra(i-1).ne.ibra .and. palavra(i).eq.ibra) then
      kf=kf+1
      nf(kf)=i-1
      endif
      
      enddo


c
c     Checa se de fato existe a coluna com o numero 
c

      
      if (ki.lt.icol) return


c
c     Checa se numero eh positivo ou negativo
c

      isig=+1

      i=ni(icol)
      if (palavra(i).eq.'-') isig=-1
  

c
c     Checa se numero estah escrito em notacao "E" ou "D"
c

      iep=0
      ie=0
      isige=+1


      do i=ni(icol),nf(icol)
      if (palavra(i).eq.'e'.or.palavra(i).eq.'E'.or.palavra(i).eq.'d'
     ?.or.palavra(i).eq.'D') ie=i
      enddo

      if (ie.ne.0) then

      iee=ie

      if (palavra(ie+1).eq.'-') then
      isige=-1
      endif

      if (palavra(ie+1).eq.'-'.or.palavra(ie+1).eq.'+') then
      ie=ie+2
      else
      ie=ie+1
      endif


      iep=0
      j=0
      k=nf(icol)-ie+1
      do i=ie,nf(icol)
      icomp=ichar(palavra(i))
      j=j+1
      iep=iep+(icomp-48)*10.d0**(k-j)
      enddo

      nf(icol)=iee-1

      endif

      expo=10.d0**(isige*iep)


c
c     Determina onde o ponto decimal estah, se nao for mumero inteiro
c

      m=0
      do i=nf(icol),ni(icol),-1
      if (palavra(i).eq.'.') m=i-nf(icol)
      enddo

c
c     Determina os algarismos do numero
c

      k=0
      do i=ni(icol),nf(icol)
      if (palavra(i).ne.'.' .and. palavra(i).ne.'+' .and. palavra(i).
     ?ne.'-') k=k+1
      enddo

      n=0
      j=0

      ist=0

      do i=ni(icol),nf(icol)
      icomp=ichar(palavra(i))
      if (icomp.ge.48 .and. icomp.le.57) then
      ist=ist+1
      j=j+1
      n=n+(icomp-48)*10.d0**(k-j)
      endif
      enddo

      valor=0.d0

      if (ist.ne.0) then
      valor=expo*isig*n*10.d0**m
      endif

      return
      end


c





      DOUBLE PRECISION FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX
      SAVE COF,STP
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+DLOG(STP*SER)
      RETURN
      END


