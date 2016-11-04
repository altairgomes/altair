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
      do i=1,nbody

        READ(10,*) lixo, mass(i)

        READ(10,*) x(3*i-2:3*i)

        READ(10,*) v(3*i-2:3*i)

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

!      call force(x,v,tm,a)
      call energy(x,v,tm,ener0,ener4)
!      ener0=vip
!      print *, x(1:3*nbody)
!      print *, v(1:3*nbody)
!      print *, a(1:3*nbody)
!      print *, ener0


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
!      call ra15(x,v,dt,xl,ll,nv,nclass,nor, 13, 16)
      call ra15m(x,v,t0,t,xl,ll,nv,nclass,nor, 20, 13, 16)

      close(13)
      close(16)
      close(15)
      
      IF(back) then
        OPEN(14,FILE='output_back.dat',BLANK='ZERO')
        OPEN(17,FILE='energy_back.dat',BLANK='ZERO')
!        call ra15(x,v,-dt,xl,ll,nv,nclass,nor, 14, 17)
        call ra15m(x,v,t,t0,-xl,ll,nv,nclass,nor, 20, 14, 17)
        close(14)
        close(17)
      endif
      
      call cpu_time(finish)
      print '("Time = ",f9.3," seconds.")',finish-start
      END program teste
     
      
      
      include "modules.f"
      include "forces.f"
!DEC$ NOFREEFORM
 !     include "radau.f"
!DEC$ FREEFORM
!      include "kepler.f"
      
      
