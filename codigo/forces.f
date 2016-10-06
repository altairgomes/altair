cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine force(x,v,tm,a)
      use constants
!      use variables
      Real (kind=wp) x(30000,3), v(30000,3), a(30000,3), a1(30000,3),
     ?r0(30000), tm, ri, a2(30000,3), d(3)
     

cccccc      force of the central body
      r0 = (/ (SQRT(dot_product(x(i,1:3),x(i,1:3))), i=1,nbody) /)
      do i=1,nbody
        a1(i,1:3) = -GaussK*GaussK*(m0+m(i))*x(i,1:3)/(r0(i)**3)
!        a1(i,1:3) = a1(i,1:3)/(r0(i))
      enddo
!      print *, a1(1:nbody,1:3)

cccccc force of other body on each other
      do i=1,nbody
        do j=1,nbody
          if (j.ne.i) then
            d = x(j,1:3) - x(i,1:3)
            ri = SQRT(dot_product(d,d))
            a2(i,1:3) = a(i,1:3) + Gaussk*Gaussk*m(j)*(d/(ri**3)-
     ?x(j,1:3)/(r0(j)**3))
!      print *, i,j,ri,m(i),a2(i,1:3)
          endif
        enddo
      enddo
!      print *, a2(1:nbody,1:3)
      
      a = a1 + a2
!      print *, a(1:nbody,1:3)
      end subroutine
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine energy(x,v,tm,vip)
      use constants
      Real (kind=wp) x(30000,3), v(30000,3), vel(30000), mtot, ener1,
     ?tm, vip, ener2
      Mtot = sum(m) + m0
      vel = (/ (sqrt(dot_product(v(i,1:3),v(i,1:3))), i=1,nbody) /)
      ener1 = sum(m(1:nbody)*vel(1:nbody)*vel(1:nbody)/2.e0_wp)
!      print *, ener1
      ener2 = (sum(m(1:nbody)*vel(1:nbody))**2)/(2*Mtot)
!      print *, ener2
      vip = ener1 - ener2
!      print *, vip
!      return
      end subroutine
      
 
