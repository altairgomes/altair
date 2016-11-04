      module constants
        implicit none
        integer, parameter :: wp = selected_real_kind(p=15)
        real(wp) :: mass0, mass(30000), ener0
        integer nbody, nsat
!        REAL (kind=wp), DIMENSION(:, :, :), ALLOCATABLE :: relfor
        REAL(wp), DIMENSION(:, :), ALLOCATABLE :: relf
        real(wp), parameter :: pi = 4.e0_wp*DATAN(1.d0)
        real(wp), parameter :: pihalf = pi/2.0e0_wp
        real(wp), parameter :: pi2 = pi*2.0e0_wp
!        real*8, parameter :: aum = 149597870700.e0_wp  ! AU -> meters
!        real*8, parameter :: mau = 1.0e0_wp/aum ! meters -> AU
!        real*8, parameter :: days = 86400.0e0_wp  ! day -> sec
!        real*8, parameter :: sday = 1.0e0_wp/days ! sec -> day
!        real*8, parameter :: radeg = 180.0e0_wp/pi ! rad -> deg
!        real*8, parameter :: degra = pi/180.0e0_wp ! deg -> rad
        
!        real*8, parameter :: G = 6.674287e-011_wp
        real(wp), parameter :: GaussK = 0.01720209895e0_wp
        contains
! Calculate vectorial product of two vectors
      function VECT(a,b)
      real(wp), dimension(3), intent(in) :: a, b
      real(wp), dimension(3) :: vect
      vect(1)=a(2)*b(3)-a(3)*b(2)
      vect(2)=a(3)*b(1)-a(1)*b(3)
      vect(3)=a(1)*b(2)-a(2)*b(1)
      
      end function vect
      
! Calculate the normal value of a vector
      function norm(a)
      real(wp), dimension(3), intent(in) :: a
      real(wp) :: norm, ab
      ab = dot_product(a,a)
      norm = sqrt(ab)
      
      end function norm
      
      function rot(a,b,c)
      real(wp), intent(in) :: a, b, c
      real(wp), dimension(3,3) :: mat1, mat2, mat3, rot
      
      mat1 = 0.0e0_wp
      mat1(3,3) = 1.0e0_wp
      mat1(1,1) = cos(a)
      mat1(2,2) = cos(a)
      mat1(1,2) = -sin(a)
      mat1(2,1) = sin(a)
      
      mat2 = 0.0e0_wp
      mat2(1,1) = 1.0e0_wp
      mat2(2,2) = cos(b)
      mat2(3,3) = cos(b)
      mat2(3,2) = -sin(b)
      mat2(2,3) = sin(b)
      
      mat3 = 0.0e0_wp
      mat3(3,3) = 1.0e0_wp
      mat3(1,1) = cos(c)
      mat3(2,2) = cos(c)
      mat3(1,2) = -sin(c)
      mat3(2,1) = sin(c)
      
      rot = matmul(mat3, matmul(mat2, mat1))
      
      end function rot
        
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
