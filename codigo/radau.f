cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*
*  Integrator RADAU by E. Everhart, Physics Department, University of Denver
*  This 15th-order version, called RA15, is written out for faster execution.
*  y'=F(y,t) is  NCLASS=1,  y"=F(y,t) is NCLASS= -2,  y"=F(y',y,t) is NCLASS=2
*  TF is t(final) - t(initial). It may be negative for backward integration.
*  NV = the number of simultaneous differential equations.
*  The dimensioning below assumes NV will not be larger than NVX.
*  LL controls sequence size. Thus SS=10**(-LL) controls the size of a term.
*  A typical LL-value is in the range 6 to 12 for this order 11 program.
*  However, if LL.LT.0 then XL is the constant sequence size used.
*  X and V enter as the starting position-velocity vector, and are output as
*  the final position-velocity vector.
*  Integration is in double precision. A 64-bit double-word is assumed.
*
      subroutine ra15(x,v,tf,xl,ll,nv,nclass,nor)
      implicit double precision (a-h,o-z)
      parameter (nvx=90000) !!!
      parameter (nsor=20)
!      parameter (nvx=18)
      real tval,pw
      dimension x(90000),v(90000),f1(nvx),fj(nvx),c(21),d(21),
     a     r(21),y(nvx),z(nvx),
     b     b(7,nvx),g(7,nvx),e(7,nvx),bd(7,nvx),h(8),w(7),u(7),nw(8)
      logical npq,nsf,nper,ncl,nes
      data nw/0,0,1,3,6,10,15,21/
      data zero, half, one,sr/0.0d0, 0.5d0, 1.0d0,1.4d0/
*  These H values are the Gauss-Radau spacings, scaled to the range 0 to 1,
*  for integrating to order 15.
      data h/         0.d0, .05626256053692215d0, .18024069173689236d0,
     a.35262471711316964d0, .54715362633055538d0, .73421017721541053d0,
     b.88532094683909577d0, .97752061356128750d0/
!      if(nv.gt.nvx)stop' **** RA15: NV > NVX ****'
*  The sum of the H-values should be 3.73333333333333333
      nper=.false.
      nsf=.false.
      ncl=nclass.eq.1
      npq=nclass.lt.2
*  y'=F(y,t)  NCL=.TRUE.    y"=F(y,t)  NCL=.FALSE.   y"=F(y',y,t) NCL=.FALSE.
*  NCLASS=1   NPQ=.TRUE.    NCLASS= -2 NPQ=.TRUE.    NCLASS= 2    NPQ=.FALSE.
*  NSF is .FALSE. on starting sequence, otherwise .TRUE.
*  NPER is .TRUE. only on last sequence of the integration.
*  NES is .TRUE. only if LL is negative. Then the sequence size is XL.
      dir=one
      if(tf.lt.zero) dir=-one
      nes=ll.lt.0
      xl=dir*dabs(xl)
      pw=1./9.
*  Evaluate the constants in the W-, U-, C-, D-, and R-vectors
      do 14 n=2,8
      ww=n+n*n
      if(ncl) ww=n
      w(n-1)=one/ww
      ww=n
  14  u(n-1)=one/ww
      do 22 k=1,nv
      if(ncl) v(k)=zero
      do 22 l=1,7
      bd(l,k)=zero
  22  b(l,k)=zero
      w1=half
      if(ncl) w1=one
      c(1)=-h(2)
      d(1)=h(2)
      r(1)=one/(h(3)-h(2))
      la=1
      lc=1
      do 73 k=3,7
      lb=la
      la=lc+1
      lc=nw(k+1)
      c(la)=-h(k)*c(lb)
      c(lc)=c(la-1)-h(k)
      d(la)=h(2)*d(lb)
      d(lc)=-c(lc)
      r(la)=one/(h(k+1)-h(2))
      r(lc)=one/(h(k+1)-h(k))
      if(k.eq.3) go to 73
      do 72 l=4,k
      ld=la+l-3
      le=lb+l-4
      c(ld)=c(le)-h(k)*c(le+1)
      d(ld)=d(le)+h(l-1)*d(le+1)
  72  r(ld)=one/(h(k+1)-h(l-1))
  73  continue
      ss=10.**(-ll)
*  The statements above are used only once in an integration to set up the
*  constants. They use less than a second of execution time.  Next set in
*  a reasonable estimate to TP based on experience. Same sign as DIR.
*  An initial first sequence size can be set with XL even with LL positive.
      tp=0.1d0*dir
      if(xl.ne.zero) tp=xl
      if(nes) tp=xl
      if(tp/tf.gt.half) tp=half*tf
      ncount=0
*     WRITE (*,3)
*  3  FORMAT(/' No. of calls, Every 10th seq.X(1),T,TM,TF')
*  An * is the symbol for writing on the monitor. The printer is unit 4.
*  Line 4000 is the starting place of the first sequence.
4000  ns=0
      nf=0
      ni=6
      tm=zero
      call force (x, v, zero, f1)
      nf=nf+1
* Line 722 is begins every sequence after the first. First find new beta-
*  values from the predicted B-values, following Eq. (2.7) in text.
 722  do 58 k=1,nv
      g(1,k)=b(1,k)+d(1)*b(2,k)+
     x  d(2)*b(3,k)+d(4)*b(4,k)+d( 7)*b(5,k)+d(11)*b(6,k)+d(16)*b(7,k)
      g(2,k)=            b(2,k)+
     x  d(3)*b(3,k)+d(5)*b(4,k)+d( 8)*b(5,k)+d(12)*b(6,k)+d(17)*b(7,k)
      g(3,k)=b(3,k)+d(6)*b(4,k)+d( 9)*b(5,k)+d(13)*b(6,k)+d(18)*b(7,k)
      g(4,k)=            b(4,k)+d(10)*b(5,k)+d(14)*b(6,k)+d(19)*b(7,k)
      g(5,k)=                         b(5,k)+d(15)*b(6,k)+d(20)*b(7,k)
      g(6,k)=                                      b(6,k)+d(21)*b(7,k)
  58  g(7,k)=                                                   b(7,k)
      t=tp
      t2=t*t
      if(ncl) t2=t
      tval=dabs(t)
*     IF(NS/10*10.EQ.NS) WRITE(*,7) NF,NS,X(1),X(2),T,TM,TF
*  7  FORMAT(1X,2I6,3F12.5,1P,2E10.2)
*  Loop 175 is 6 iterations on first sequence and two iterations therafter.
      do 175 m=1,ni
*  Loop 174 is for each substep within a sequence.
      do 174 j=2,8
      jd=j-1
      jdm=j-2
      s=h(j)
      q=s
      if(ncl) q=one
*  Use Eqs. (2.9) and (2.10) of text to predict positions at each aubstep.
*  These collapsed series are broken into two parts because an otherwise
*  excellent  compiler could not handle the complicated expression.
      do 130 k=1,nv
      a=w(3)*b(3,k)+s*(w(4)*b(4,k)+s*(w(5)*b(5,k)+s*(w(6)*b(6,k)+
     v   s*w(7)*b(7,k))))
      y(k)=x(k)+q*(t*v(k)+t2*s*(f1(k)*w1+s*(w(1)*b(1,k)+s*(w(2)*b(2,k)
     x  +s*a))))
      if(npq) go to 130
*  Next are calculated the velocity predictors need for general class II.
      a=u(3)*b(3,k)+s*(u(4)*b(4,k)+s*(u(5)*b(5,k)+s*(u(6)*b(6,k)+
     t    s*u(7)*b(7,k))))
      z(k)=v(k)+s*t*(f1(k)+s*(u(1)*b(1,k)+s*(u(2)*b(2,k)+s*a)))
 130  continue
*  Find forces at each substep.
      call force(y,z,tm+s*t,fj)
      nf=nf+1
      do 171 k=1,nv
*  Find G-value for the force FJ found at the current substep. This
*  section, including the many-branched GOTO, uses Eq. (2.4) of text.
      temp=g(jd,k)
      gk=(fj(k)-f1(k))/s
      go to (102,102,103,104,105,106,107,108),j
 102  g(1,k)=gk
      go to 160
 103  g(2,k)=(gk-g(1,k))*r(1)
      go to 160
 104  g(3,k)=((gk-g(1,k))*r(2)-g(2,k))*r(3)
      go to 160
 105  g(4,k)=(((gk-g(1,k))*r(4)-g(2,k))*r(5)-g(3,k))*r(6)
      go to 160
 106  g(5,k)=((((gk-g(1,k))*r(7)-g(2,k))*r(8)-g(3,k))*r(9)-
     x     g(4,k))*r(10)
      go to 160
 107  g(6,k)=(((((gk-g(1,k))*r(11)-g(2,k))*r(12)-g(3,k))*r(13)-
     x     g(4,k))*r(14)-g(5,k))*r(15)
      go to 160
 108  g(7,k)=((((((gk-g(1,k))*r(16)-g(2,k))*r(17)-g(3,k))*r(18)-
     x     g(4,k))*r(19)-g(5,k))*r(20)-g(6,k))*r(21)
*  Upgrade all B-values.
 160  temp=g(jd,k)-temp
      b(jd,k)=b(jd,k)+temp
*  TEMP is now the improvement on G(JD,K) over its former value.
*  Now we upgrade the B-value using this dfference in the one term.
*  This section is based on Eq. (2.5).
      go to (171,171,203,204,205,206,207,208),j
 203  b(1,k)=b(1,k)+c(1)*temp
      go to 171
 204  b(1,k)=b(1,k)+c(2)*temp
      b(2,k)=b(2,k)+c(3)*temp
      go to 171
 205  b(1,k)=b(1,k)+c(4)*temp
      b(2,k)=b(2,k)+c(5)*temp
      b(3,k)=b(3,k)+c(6)*temp
      go to 171
 206  b(1,k)=b(1,k)+c(7)*temp
      b(2,k)=b(2,k)+c(8)*temp
      b(3,k)=b(3,k)+c(9)*temp
      b(4,k)=b(4,k)+c(10)*temp
      go to 171
 207  b(1,k)=b(1,k)+c(11)*temp
      b(2,k)=b(2,k)+c(12)*temp
      b(3,k)=b(3,k)+c(13)*temp
      b(4,k)=b(4,k)+c(14)*temp
      b(5,k)=b(5,k)+c(15)*temp
      go to 171
 208  b(1,k)=b(1,k)+c(16)*temp
      b(2,k)=b(2,k)+c(17)*temp
      b(3,k)=b(3,k)+c(18)*temp
      b(4,k)=b(4,k)+c(19)*temp
      b(5,k)=b(5,k)+c(20)*temp
      b(6,k)=b(6,k)+c(21)*temp
 171  continue
 174  continue
      if(nes.or.m.lt.ni) go to 175
*  Integration of sequence is over. Next is sequence size control.
      hv=zero
      do 635 k=1,nv
 635  hv=dmax1(hv,dabs(b(7,k)))
      hv=hv*w(7)/tval**7
 175  continue
      if (nsf) go to 180
      if(.not.nes) tp=(ss/hv)**pw*dir
      if(nes) tp=xl
      if(nes) go to 170
      if(tp/t.gt.one) go to 170
   8  format (2x,2i2,2d18.10)
      tp=.8d0*tp
      ncount=ncount+1
      if(ncount.gt.10) return
*     IF(NCOUNT.GT.1) WRITE (4,8) NOR,NCOUNT,T,TP
*  Restart with 0.8x sequence size if new size called for is smaller than
*  originally chosen starting sequence size on first sequence.
      go to 4000
 170  nsf=.true.
* Loop 35 finds new X and V values at end of sequence using Eqs. (2.11),(2.12)
 180  do 35 k=1,nv
      x(k)=x(k)+v(k)*t+t2*(f1(k)*w1+b(1,k)*w(1)+b(2,k)*w(2)+b(3,k)*w(3)
     x    +b(4,k)*w(4)+b(5,k)*w(5)+b(6,k)*w(6)+b(7,k)*w(7))
      if(ncl) go to 35
      v(k)=v(k)+t*(f1(k)+b(1,k)*u(1)+b(2,k)*u(2)+b(3,k)*u(3)
     v    +b(4,k)*u(4)+b(5,k)*u(5)+b(6,k)*u(6)+b(7,k)*u(7))
  35  continue
      tm=tm+t
      ns=ns+1
*  Return if done.
      if(.not.nper) go to 78
*     WRITE(*,7) NF,NS,X(1),X(2),T,TM,TF
*     WRITE(4,7) NF,NS
      return
*  Control on size of next sequence and adjust last sequence to exactly
*  cover the integration span. NPER=.TRUE. set on last sequence.
78    call force (x,v,tm,f1)
      nf=nf+1
      if(nes) go to 341
      tp=dir*(ss/hv)**pw
      if(tp/t.gt.sr) tp=t*sr
 341  if(nes) tp=xl
      if(dir*(tm+tp).lt.dir*tf-1.d-8) go to 77
      tp=tf-tm
      nper=.true.
*  Now predict B-values for next step. The predicted values from the preceding
*  sequence were saved in the E-matrix. Te correction BD between the actual
*  B-values found and these predicted values is applied in advance to the
*  next sequence. The gain in accuracy is significant. Using Eqs. (2.13):
  77  q=tp/t
      do 39 k=1,nv
      if(ns.eq.1) go to 31
      do 20 j=1,7
  20  bd(j,k)=b(j,k)-e(j,k)
  31  e(1,k)=      q*(b(1,k)+ 2.d0*b(2,k)+ 3.d0*b(3,k)+
     x           4.d0*b(4,k)+ 5.d0*b(5,k)+ 6.d0*b(6,k)+ 7.d0*b(7,k))
      e(2,k)=                q**2*(b(2,k)+ 3.d0*b(3,k)+
     y           6.d0*b(4,k)+10.d0*b(5,k)+15.d0*b(6,k)+21.d0*b(7,k))
      e(3,k)=                             q**3*(b(3,k)+
     z           4.d0*b(4,k)+10.d0*b(5,k)+20.d0*b(6,k)+35.d0*b(7,k))
      e(4,k)=   q**4*(b(4,k)+ 5.d0*b(5,k)+15.d0*b(6,k)+35.d0*b(7,k))
      e(5,k)=                q**5*(b(5,k)+ 6.d0*b(6,k)+21.d0*b(7,k))
      e(6,k)=                             q**6*(b(6,k)+ 7.d0*b(7,k))
      e(7,k)=                                           q**7*b(7,k)
      do 39 l=1,7
  39  b(l,k)=e(l,k)+bd(l,k)
*  Two iterations for every sequence after the first.
      ni=2
      go to 722
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
C 
C  $$$$  RA15.FOR   
C 
      SUBROUTINE RA15M(X,V,TD,TF0,XL,LL,NV,NCLASS,NOR,nsor, out1, out2) 
C  Integrator RADAU by E. Everhart, Physics Department, University of Denver 
C  This 15th-order version, called RA15, is written out for faster execution. 
C  y'=F(y,t) is  NCLASS=1,  y"=F(y,t) is NCLASS= -2,  y"=F(y',y,t) is NCLASS=2 
C  TF is t(final) - t(initial). It may be negative for backward integration. 
C  NV = the number of simultaneous differential equations. 
C  The dimensioning below assumes NV will not be larger than 50. 
C  LL controls sequence size. Thus SS=10**(-LL) controls the size of a term. 
C  A typical LL-value is in the range 6 to 12 for this order 11 program. 
C  However, if LL.LT.0 then XL is the constant sequence size used. 
C  X and V enter as the starting position-velocity vector, and are output as 
C  the final position-velocity vector. 
C  Integration is in double precision. A 64-bit double-word is assumed.
      use constants 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *4 TVAL,PW! , mu
      DIMENSION X(90000),V(90000),F1(90000),FJ(90000),C(21),D(21),
     & R(21),Y(90000),Z(90000), 
     A     B(7,90000),G(7,90000),E(7,90000),BD(7,90000),H(8),W(7),
     & U(7),NW(8), ener4(4)!, elemx(11)
      LOGICAL NPQ,NSF,NPER,NCL,NES
      integer, intent(in) :: out1
      integer, intent(in) :: out2
      integer kkk
!      EXTERNAL FORCE,SORTIE   
      DATA NW/0,0,1,3,6,10,15,21/ 
      DATA ZERO, HALF, ONE,SR/0.0D0, 0.5D0, 1.0D0,1.4D0/ 
C  These H values are the Gauss-Radau spacings, scaled to the range 0 to 1, 
C  for integrating to order 15. 
      DATA H/         0.D0, .05626256053692215D0, .18024069173689236D0, 
     A.35262471711316964D0, .54715362633055538D0, .73421017721541053D0, 
     B.88532094683909577D0, .97752061356128750D0/ 
     
C  The sum of the H-values should be 3.73333333333333333 
      NPER=.FALSE. 
      NSF=.FALSE. 
      NCL=NCLASS.EQ.1 
      NPQ=NCLASS.LT.2 
C  y'=F(y,t)  NCL=.TRUE.    y"=F(y,t)  NCL=.FALSE.   y"=F(y',y,t) NCL=.FALSE. 
C  NCLASS=1   NPQ=.TRUE.    NCLASS= -2 NPQ=.TRUE.    NCLASS= 2    NPQ=.FALSE. 
C  NSF is .FALSE. on starting sequence, otherwise .TRUE. 
C  NPER is .TRUE. only on last sequence of the integration. 
C  NES is .TRUE. only if LL is negative. Then the sequence size is XL. 
      DIR=ONE 
      TF=TF0-TD                 
      IF(TF.LT.ZERO) DIR=-ONE   
      NES=LL.LT.0 
      XL=DIR*DABS(XL) 
      PW=1./9. 
C  Evaluate the constants in the W-, U-, C-, D-, and R-vectors 
      DO 14 N=2,8 
      WW=N+N*N 
      IF(NCL) WW=N 
      W(N-1)=ONE/WW 
      WW=N 
  14  U(N-1)=ONE/WW 
      DO 22 K=1,NV 
      IF(NCL) V(K)=ZERO 
      DO 22 L=1,7 
      BD(L,K)=ZERO 
  22  B(L,K)=ZERO 
      W1=HALF 
      IF(NCL) W1=ONE 
      C(1)=-H(2) 
      D(1)=H(2) 
      R(1)=ONE/(H(3)-H(2)) 
      LA=1 
      LC=1 
      DO 73 K=3,7 
      LB=LA 
      LA=LC+1 
      LC=NW(K+1) 
      C(LA)=-H(K)*C(LB) 
      C(LC)=C(LA-1)-H(K) 
      D(LA)=H(2)*D(LB) 
      D(LC)=-C(LC) 
      R(LA)=ONE/(H(K+1)-H(2)) 
      R(LC)=ONE/(H(K+1)-H(K)) 
      IF(K.EQ.3) GO TO 73 
      DO 72 L=4,K 
      LD=LA+L-3 
      LE=LB+L-4 
      C(LD)=C(LE)-H(K)*C(LE+1) 
      D(LD)=D(LE)+H(L-1)*D(LE+1) 
  72  R(LD)=ONE/(H(K+1)-H(L-1)) 
  73  CONTINUE 
      SS=10.**(-LL) 
C  The statements above are used only once in an integration to set up the 
C  constants. They use less than a second of execution time.  Next set in 
C  a reasonable estimate to TP based on experience. Same sign as DIR. 
C  An initial first sequence size can be set with XL even with LL positive. 
      TP=0.1D0*DIR 
      IF(XL.NE.ZERO) TP=XL 
      IF(NES) TP=XL 
      IF(TP/TF.GT.HALF) TP=HALF*TF 
      NCOUNT=0 
 
C  An * is the symbol for writing on the monitor. The printer is unit 4. 
C  Line 4000 is the starting place of the first sequence. 
4000  NS=0 
      NF=0 
      NI=6 
      TM=TD      
      CALL FORCE (X, V, TM, F1) 
      NF=NF+1 
C Line 722 is begins every sequence after the first. First find new beta- 
C  values from the predicted B-values, following Eq. (2.7) in text. 
 722  DO 58 K=1,NV 
      G(1,K)=B(1,K)+D(1)*B(2,K)+ 
     X  D(2)*B(3,K)+D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K) 
      G(2,K)=            B(2,K)+ 
     X  D(3)*B(3,K)+D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K) 
      G(3,K)=B(3,K)+D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K) 
      G(4,K)=            B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K) 
      G(5,K)=                         B(5,K)+D(15)*B(6,K)+D(20)*B(7,K) 
      G(6,K)=                                      B(6,K)+D(21)*B(7,K) 
  58  G(7,K)=                                                   B(7,K) 

      T=TP 
      T2=T*T 
      IF(NCL) T2=T 
      TVAL=DABS(T) 
      IF(NS/nsor*nsor.EQ.NS) then  
        call energy(X,V,TM,VIP,ener4) 
        WRITE(out2,*) TM, VIP, (VIP-ener0)/ener0, ener4
        write(out1,*) tm, x(1:3*nbody)
        do kkk=1,nbody
        write(19, *) TM, kkk, relf(kkk,:)
        enddo
        mu = k*k*(mass0+mass(1))
        call PV2ALZDZ(MU,x(1:3*nbody),V(1:3*nbody),ELEMX)
        write(15, *), elemx
!        relfor(kkk,:,:) = relf
!        kkk = kkk + 1
      end if 
 
C  Loop 175 is 6 iterations on first sequence and two iterations therafter. 
      DO 175 M=1,NI 
C  Loop 174 is for each substep within a sequence. 
      DO 174 J=2,8 
      JD=J-1 
      JDM=J-2 
      S=H(J) 
      Q=S 
      IF(NCL) Q=ONE 
C  Use Eqs. (2.9) and (2.10) of text to predict positions at each aubstep. 
C  These collapsed series are broken into two parts because an otherwise 
C  excellent  compiler could not handle the complicated expression. 
      DO 130 K=1,NV 
      A=W(3)*B(3,K)+S*(W(4)*B(4,K)+S*(W(5)*B(5,K)+S*(W(6)*B(6,K)+ 
     V   S*W(7)*B(7,K)))) 
      Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*(W(1)*B(1,K)+S*(W(2)*B(2,K) 
     X  +S*A)))) 
      IF(NPQ) GO TO 130 
C  Next are calculated the velocity predictors need for general class II. 
      A=U(3)*B(3,K)+S*(U(4)*B(4,K)+S*(U(5)*B(5,K)+S*(U(6)*B(6,K)+ 
     T    S*U(7)*B(7,K)))) 
      Z(K)=V(K)+S*T*(F1(K)+S*(U(1)*B(1,K)+S*(U(2)*B(2,K)+S*A))) 
 130  CONTINUE 

C  Find forces at each substep. 
      CALL FORCE(Y,Z,TM+S*T,FJ) 

      NF=NF+1 
      DO 171 K=1,NV 
C  Find G-value for the force FJ found at the current substep. This 
C  section, including the many-branched GOTO, uses Eq. (2.4) of text. 
      TEMP=G(JD,K) 
      GK=(FJ(K)-F1(K))/S 
      GO TO (102,102,103,104,105,106,107,108),J 
 102  G(1,K)=GK 
      GO TO 160 
 103  G(2,K)=(GK-G(1,K))*R(1) 
      GO TO 160 
 104  G(3,K)=((GK-G(1,K))*R(2)-G(2,K))*R(3) 
      GO TO 160 
 105  G(4,K)=(((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*R(6) 
      GO TO 160 
 106  G(5,K)=((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*R(9)- 
     X     G(4,K))*R(10) 
      GO TO 160 
 107  G(6,K)=(((((GK-G(1,K))*R(11)-G(2,K))*R(12)-G(3,K))*R(13)- 
     X     G(4,K))*R(14)-G(5,K))*R(15) 
      GO TO 160 
 108  G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-G(3,K))*R(18)- 
     X     G(4,K))*R(19)-G(5,K))*R(20)-G(6,K))*R(21) 
C  Upgrade all B-values. 
 160  TEMP=G(JD,K)-TEMP 
      B(JD,K)=B(JD,K)+TEMP 
C  TEMP is now the improvement on G(JD,K) over its former value. 
C  Now we upgrade the B-value using this dfference in the one term. 
C  This section is based on Eq. (2.5). 
      GO TO (171,171,203,204,205,206,207,208),J 
 203  B(1,K)=B(1,K)+C(1)*TEMP 
      GO TO 171 
 204  B(1,K)=B(1,K)+C(2)*TEMP 
      B(2,K)=B(2,K)+C(3)*TEMP 
      GO TO 171 
 205  B(1,K)=B(1,K)+C(4)*TEMP 
      B(2,K)=B(2,K)+C(5)*TEMP 
      B(3,K)=B(3,K)+C(6)*TEMP 
      GO TO 171 
 206  B(1,K)=B(1,K)+C(7)*TEMP 
      B(2,K)=B(2,K)+C(8)*TEMP 
      B(3,K)=B(3,K)+C(9)*TEMP 
      B(4,K)=B(4,K)+C(10)*TEMP 
      GO TO 171 
 207  B(1,K)=B(1,K)+C(11)*TEMP 
      B(2,K)=B(2,K)+C(12)*TEMP 
      B(3,K)=B(3,K)+C(13)*TEMP 
      B(4,K)=B(4,K)+C(14)*TEMP 
      B(5,K)=B(5,K)+C(15)*TEMP 
      GO TO 171 
 208  B(1,K)=B(1,K)+C(16)*TEMP 
      B(2,K)=B(2,K)+C(17)*TEMP 
      B(3,K)=B(3,K)+C(18)*TEMP 
      B(4,K)=B(4,K)+C(19)*TEMP 
      B(5,K)=B(5,K)+C(20)*TEMP 
      B(6,K)=B(6,K)+C(21)*TEMP 
 171  CONTINUE 
 174  CONTINUE 

      IF(NES.OR.M.LT.NI) GO TO 175 
C  Integration of sequence is over. Next is sequence size control. 
      HV=ZERO 
      DO 635 K=1,NV 
 635  HV=DMAX1(HV,DABS(B(7,K))) 
      HV=HV*W(7)/TVAL**7 
 175  CONTINUE 
      IF (NSF) GO TO 180 
      IF(.NOT.NES) TP=(SS/HV)**PW*DIR 
      IF(NES) TP=XL 
      IF(NES) GO TO 170 
      IF(TP/T.GT.ONE) GO TO 170 
   8  FORMAT (2X,2I2,2D18.10) 
      TP=.8D0*TP 
      NCOUNT=NCOUNT+1 
      IF(NCOUNT.GT.10) RETURN 
c      IF(NCOUNT.GT.1) WRITE(14,8) NOR,NCOUNT,T,TP 
C  Restart with 0.8x sequence size if new size called for is smaller than 
C  originally chosen starting sequence size on first sequence. 
      GO TO 4000 
 170  NSF=.TRUE. 
C Loop 35 finds new X and V values at end of sequence using Eqs. (2.11),(2.12) 
 180  DO 35 K=1,NV 
      X(K)=X(K)+V(K)*T+T2*(F1(K)*W1+B(1,K)*W(1)+B(2,K)*W(2)+B(3,K)*W(3) 
     X    +B(4,K)*W(4)+B(5,K)*W(5)+B(6,K)*W(6)+B(7,K)*W(7)) 
      IF(NCL) GO TO 35 
      V(K)=V(K)+T*(F1(K)+B(1,K)*U(1)+B(2,K)*U(2)+B(3,K)*U(3) 
     V    +B(4,K)*U(4)+B(5,K)*U(5)+B(6,K)*U(6)+B(7,K)*U(7)) 
  35  CONTINUE 
      NS=NS+1 
      if(NES) then 
        TM=TD+NS*T     
      else 
        TM=TM+T 
      end if 
C  Return if done. 
      IF(.NOT.NPER) GO TO 78 
        call energy(X,V,TD+NS*T,VIP,ener4) 
        WRITE(out2,*) TM, VIP, (VIP-ener0)/ener0, ener4
        write(out1,*) tm, x(1:3*nbody)
        do kkk=1,nbody
        write(19, *) TM, kkk, relf(kkk,:)
        enddo
        mu = k*k*(mass0+mass(1))
        call PV2ALZDZ(MU,x(1:3*nbody),V(1:3*nbody),ELEMX)
        write(15, *), elemx
!        relfor(kkk,:,:) = relf
!        kkk = kkk + 1
      RETURN 
C  Control on size of next sequence and adjust last sequence to exactly 
C  cover the integration span. NPER=.TRUE. set on last sequence. 
78    CALL FORCE (X,V,TM,F1) 
      NF=NF+1 
      IF(NES) GO TO 341 
      TP=DIR*(SS/HV)**PW 
      IF(TP/T.GT.SR) TP=T*SR 
 341  IF(NES) TP=XL 
      IF(DIR*(TM+TP).LT.DIR*(TF0)-1.D-8) GO TO 77
      TP=TF0-TM 
      NPER=.TRUE. 
C  Now predict B-values for next step. The predicted values from the preceding 
C  sequence were saved in the E-matrix. Te correction BD between the actual 
C  B-values found and these predicted values is applied in advance to the 
C  next sequence. The gain in accuracy is significant. Using Eqs. (2.13): 
  77  Q=TP/T 
      DO 39 K=1,NV 
      IF(NS.EQ.1) GO TO 31
 
      DO 20 J=1,7 
  20  BD(J,K)=B(J,K)-E(J,K) 
  31  E(1,K)=      Q*(B(1,K)+ 2.D0*B(2,K)+ 3.D0*B(3,K)+ 
     X           4.D0*B(4,K)+ 5.D0*B(5,K)+ 6.D0*B(6,K)+ 7.D0*B(7,K)) 
      E(2,K)=                Q**2*(B(2,K)+ 3.D0*B(3,K)+ 
     Y           6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+21.D0*B(7,K)) 
      E(3,K)=                             Q**3*(B(3,K)+ 
     Z           4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+35.D0*B(7,K)) 
      E(4,K)=   Q**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+35.D0*B(7,K)) 
      E(5,K)=                Q**5*(B(5,K)+ 6.D0*B(6,K)+21.D0*B(7,K)) 
      E(6,K)=                             Q**6*(B(6,K)+ 7.D0*B(7,K)) 
      E(7,K)=                                           Q**7*B(7,K) 

      DO 39 L=1,7 
  39  B(L,K)=E(L,K)+BD(L,K) 

C  Two iterations for every sequence after the first. 
      NI=2 
      if(.true.) GO TO 722 
      END 
