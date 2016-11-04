!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine force(x,v,tm,a)
      use constants
!      use variables
      Real(wp), intent(in) :: x(90000), v(90000), tm
      real(wp), intent(out) :: a(90000)
      real(wp) :: a1(90000), r0(30000), rij, a2(90000), d(3), vip, tt, td, ff1(3), ff2(3)
     

!cccccc      force of the central body
      r0 = (/ (norm(x(3*i-2:3*i)), i=1,nbody) /)  ! calculates module for each vector
      a1 = (/ (-GaussK*GaussK*(mass0+mass(i))*x(3*i-2:3*i)/(r0(i)**3), i=1,nbody) /) ! calculates force for every object
!      print *, a1(1:3*nbody)
!      print *, 'a1', a1(1:3*nbody)

!cccccc force of other bodies on each other
      a2 = 0.0e0_wp
      relf = 0.0e0_wp
      do i=1,nbody
!        td = SQRT(dot_product(a1(3*i-2:3*i),a1(3*i-2:3*i)))
        do j=1,nbody
          if (j.ne.i) then
            d = x(3*j-2:3*j) - x(3*i-2:3*i)
!            print *, 'ji', i, j
!            print *, 'd', d
            rij = SQRT(dot_product(d,d))
!            print *, 'rij', rij
            ff1 = d/(rij**3)
!            print *, 'ff1', ff1
            ff2 = (x(3*j-2:3*j)/(r0(j)**3))
!            print *, 'ff2', ff2
            a2(3*i-2:3*i) = a2(3*i-2:3*i) + Gaussk*Gaussk*mass(j)*(ff1 - ff2)
!            tt = SQRT(dot_product(ff1,ff1))
!            relf(i,j) = tt/td
          endif
        enddo
      enddo
!      print *, 'a2', a2(1:3*nbody)
      
      a(1:3*nbody) = a1(1:3*nbody) + a2(1:3*nbody)

      end subroutine
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine energy(x,v,tm,vip, ener4)
      use constants
      Real (kind=wp):: x(90000), v(90000), vel(90000), mtot, ener1, tm, vip, ener2, enpot1, enpot2, r0(nsat), d(3), rij, encin1, &
      encin2, ener4(4), encin2a
     
      Mtot = sum(mass(1:nsat)) + mass0
!      print *, 'mtot', mtot
      vel = (/ (norm(v(3*i-2:3*i)), i=1,nsat) /) ! calculates module for each vector
      encin1 = 0.5e0_wp*sum(mass(1:nsat)*vel(1:nsat)*vel(1:nsat))
      
!      print *, 'encin1', encin1

      encin2a = sum(mass(1:nsat)*vel(1:nsat))
      encin2 = (encin2a*encin2a)/(2.e0_wp*Mtot)
!      print *, 'encin2', encin2

      r0 = (/ (norm(x(3*i-2:3*i)), i=1,nsat) /)
      enpot1 = GaussK*GaussK*mass0*sum(mass(1:nsat)/r0(1:nsat))
!      print *, 'enpot1', enpot1

      enpot2 = 0.0e0_wp
      do i=1,nsat-1
        do j=i+1,nsat
          d = x(3*j-2:3*j) - x(3*i-2:3*i)
          rij = norm(d)
          enpot2 = enpot2 + GaussK*GaussK*mass(j)*mass(i)/rij
        enddo
      enddo
!      print *, 'enpot2', enpot2
      ener4 = (/ encin1,encin2,enpot1,enpot2 /)
      vip = (encin1-encin2) - (enpot1+enpot2)
!      print *, encin1, encin2, enpot1, enpot2, vip

      end subroutine
      
      
      
!c###

!c code: debbug energy integral
!c remove mutual perturbations


!c for delta alfa cos delta mudar observações de 001 para 010
 
