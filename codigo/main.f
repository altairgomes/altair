      module constants
        implicit none
        integer, parameter :: wp = selected_real_kind(p=15)
        real (kind=wp) m0, m(30000)
        integer nbody
        real*8, parameter :: pi = 4.e0_wp*DATAN(1.d0)
        real*8, parameter :: aum = 149597870700.e0_wp  ! AU -> meters
        real*8, parameter :: mau = 1.0e0_wp/aum ! meters -> AU
        real*8, parameter :: days = 86400.0e0_wp  ! day -> sec
        real*8, parameter :: sday = 1.0e0_wp/days ! sec -> day
        real*8, parameter :: radeg = 180.0e0_wp/pi ! rad -> deg
        real*8, parameter :: degra = pi/180.0e0_wp ! deg -> rad
        
        real*8, parameter :: G = 6.674287e-011_wp
        real*8, parameter :: GaussK = 0.01720209895_wp
        
      end module constants
        
      MODULE variables
      use constants
      REAL (kind=wp) x(30000,3), v(30000,3), ener0, t
        real*8, parameter :: alfa0 = 29.46086126075579_wp  !alfa0+pi/2 Neptune
        real*8, parameter :: delta0 = 46.5951892092859_wp  !pi/2-delta0 Neptune
        real*8, parameter :: wpoint = 536.3128492_wp   ! Neptune
        real*8, parameter :: J2 = 3408.428530717952e-006_wp ! Neptune
        real*8, parameter :: J4 = -33.398917590066e-006_wp ! Neptune
        real*8, parameter :: t0 = 2451545.0_wp  ! TDB
        real*8, parameter :: step = 0.025e0_wp!*days
        integer, parameter :: n = 20000
      end module variables
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      program teste
      use constants
      use variables
!      implicit none
      REAL (kind=wp) xl, dt, a(30000,3), tm, vip
      integer nv, nor, nclass, ll
      character*30 lixo
      
ccccc reading input file cccc

      OPEN(10, FILE='input.dat')
      READ(10,*) nbody
!      print *, nbody
      read(10,*) lixo, m0
!      print *, lixo, m0
      do i=1,nbody
        READ(10,*) lixo, m(i)
!        print *, lixo, m(i)
        READ(10,*) x(i,1), x(i,2), x(i,3)
!        print *, x(i,1), x(i,2), x(i,3)
        READ(10,*) v(i,1), v(i,2), v(i,3)
!        print *, v(i,1), v(i,2), v(i,3)
      enddo
      CLOSE(10)
!      print *, x(1:nbody,1:3)
      call force(x,v,tm,a)
      call energy(x,v,tm,vip)
      ener0=vip


ccccccccccccccccccccccccccccc      
      
      
      OPEN(13,FILE='output.dat',BLANK='ZERO')
      OPEN(16,FILE='energy.dat',BLANK='ZERO')
      nv = 3*nbody
      nor = 15
      nclass=-2
      ll=-1
      xl = step!*sday

!      do i=0,n
        t = t0 + n*step!*sday
        dt = n*step
!        call energy(x,v)
!        ener0 = ener
!        write (13,*) t,x(1)/1000.,x(2)/1000.,x(3)/1000.,v(1)/1000.,
!     ?  v(2)/1000., v(3)/1000., ener
!        print *, t, SQRT(x(1)**2+x(2)**2+x(3)**2), ener
!        call ra15(x,v,dt,xl,ll,nv,nclass,nor)
!        call ra15m(x,v,t0,t,xl,ll,nv,nclass,nor, 20, force, energy)
        !passar tempo para o calculo de forca e receber de volta posicao e velocidade
!      enddo
      close(13)
      close(16)
      END program teste
      
      include "radau.f"
      include "forces.f"
      
