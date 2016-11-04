      module constants
        implicit none
        integer, parameter :: wp = selected_real_kind(p=15)
        real (kind=wp) mass0, mass(30000), ener0
        integer nbody, nsat
!        REAL (kind=wp), DIMENSION(:, :, :), ALLOCATABLE :: relfor
        REAL (kind=wp), DIMENSION(:, :), ALLOCATABLE :: relf
        real*8, parameter :: pi = 4.e0_wp*DATAN(1.d0)
        real*8, parameter :: pihalf = pi/2.0e0_wp
        real*8, parameter :: pi2 = pi*2.0e0_wp
!        real*8, parameter :: aum = 149597870700.e0_wp  ! AU -> meters
!        real*8, parameter :: mau = 1.0e0_wp/aum ! meters -> AU
!        real*8, parameter :: days = 86400.0e0_wp  ! day -> sec
!        real*8, parameter :: sday = 1.0e0_wp/days ! sec -> day
!        real*8, parameter :: radeg = 180.0e0_wp/pi ! rad -> deg
!        real*8, parameter :: degra = pi/180.0e0_wp ! deg -> rad
        
!        real*8, parameter :: G = 6.674287e-011_wp
        real*8, parameter :: GaussK = 0.01720209895_wp
        
      end module constants
        
      MODULE variables
      use constants
      REAL (kind=wp) x(90000), v(90000), t0, step
      Integer*8 nstep
        real*8, parameter :: alfa0 = 29.46086126075579_wp  !alfa0+pi/2 Neptune
        real*8, parameter :: delta0 = 46.5951892092859_wp  !pi/2-delta0 Neptune
        real*8, parameter :: wpoint = 536.3128492_wp   ! Neptune
        real*8, parameter :: J2 = 3408.428530717952e-006_wp ! Neptune
        real*8, parameter :: J4 = -33.398917590066e-006_wp ! Neptune
      end module variables
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      program teste
      use constants
      use variables
!      implicit none
      REAL (kind=wp) xl, dt, a(90000), tm, vip, t,ener4(4)
      REAL*8 start, finish
      integer :: nv, nor, nclass, ll
      character*30 lixo
      logical back
      call cpu_time(start)      
!ccccc reading input file cccc

      OPEN(10, FILE='input.dat')
      READ(10,*) nbody, nsat

      read(10,*) lixo, mass0
      do i=0,nbody-1

        READ(10,*) lixo, mass(i+1)

        READ(10,*) x(3*i+1:3*i+3)

        READ(10,*) v(3*i+1:3*i+3)

      enddo
      CLOSE(10)

      OPEN(10, FILE='integ_ini.dat')
      read(10,*) t0
      read(10,*) step
      read(10,*) nstep
      read(10,*) back
      CLOSE(10)

!cccccc initial energy cccccccccc

!      ALLOCATE ( relfor(10000, Nbody, Nbody) )
      ALLOCATE ( relf(Nbody, Nbody) )

      call force(x,v,tm,a)
      call energy(x,v,tm,ener0,ener4)
!      ener0=vip
!      print *, x(1:3*nbody)
!      print *, v(1:3*nbody)
      print *, a(1:3*nbody)
      print *, ener0


!cccccccccccccccccccccccccccccccc
      
      
      OPEN(13,FILE='output.dat',BLANK='ZERO')
      OPEN(16,FILE='energy.dat',BLANK='ZERO')
      OPEN(19,FILE='relat_forces.dat',BLANK='ZERO')
      OPEN(15,FILE='elements.dat', BLANK='ZERO')
      nv = 3*nbody
      nor = 15
      nclass=-2
      ll=-1
      xl = step!*sday

      t = t0 + nstep*step!*sday
      dt = nstep*step
!      call ra15(x,v,dt,xl,ll,nv,nclass,nor)
      call ra15m(x,v,t0,t,xl,ll,nv,nclass,nor, 20, 13, 16)

      close(13)
      close(16)
      close(15)
      
      IF(back) then
        OPEN(14,FILE='output_back.dat',BLANK='ZERO')
        OPEN(17,FILE='energy_back.dat',BLANK='ZERO')
        call ra15m(x,v,t,t0,-xl,ll,nv,nclass,nor, 20, 14, 17)
        close(14)
        close(17)
      endif
      
      call cpu_time(finish)
      print '("Time = ",f9.3," seconds.")',finish-start
      END program teste
      
      include "forces.f"
!DEC$ NOFREEFORM
 !     include "radau.f"
!DEC$ FREEFORM
!      include "kepler.f"
      
      
