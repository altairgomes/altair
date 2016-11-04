      SUBROUTINE PV2ALZDZ(MU,POS,VIT,ELEMX)

    ! ELEMX(1)=a;     ELEMX(2)=L;     ELEMX(3)=k;  ELEMX(4)=h
    ! ELEMX(5)=q;     ELEMX(6)=p;     ELEMX(7)=n;  ELEMX(8)=e
    ! ELEMX(9)=varpi; ELEMX(10)=sis2; ELEMX(11)=Omega
      use constants
      IMPLICIT NONE
      REAL(wp), INTENT(IN) :: MU
      REAL(wp), DIMENSION(3), INTENT(IN) :: POS, VIT
      REAL(wp), DIMENSION(11), INTENT(OUT) :: ELEMX

      REAL(wp):: R, A, GG, F, P, EX, MUinv, Rinv, Inc, sinOme, cosOme, sinpomeplusf, cospomeplusf, V2, sinInc,  &
         cosInc, GGinv, sinf, cosf, EXinv, pomeplusf, pome, Ainv, EXsinE, EXcosE, Moy, Ome, E, E1, sinI2, Long
      REAL(wp), DIMENSION(3):: G, xyz, Gnorm

      MUinv=1.E0_wp/MU; R=NORM(POS); Rinv=1.E0_wp/R; V2=dot_product(VIT,VIT)

      E1=0.5E0_wp*V2-MU*Rinv ! Energie de m autour de M (different de l'energie totale)
      IF(E1>=0.E0_wp)THEN
         WRITE(*,*)"erreur : orbite non elliptique"; STOP
      ENDIF

      A=-0.5E0_wp*MU/E1; Ainv=1.E0_wp/A !demi-grand axe

      G=VECT(POS,VIT); GG=NORM(G); GGinv=1.E0_wp/GG; Gnorm=G*GGinv !moment cinetique

      EX=SQRT((1.E0_wp-R*Ainv)**2+dot_product(POS,VIT)**2*MUinv*Ainv) !excentricite

      cosInc=Gnorm(3)
      sinInc=SQRT(Gnorm(1)**2+Gnorm(2)**2)
      Inc=ATAN2(sinInc,cosInc)

      IF((ABS(Inc)>1.E-14_wp).AND.(ABS(Inc)<(pi-1.E-14_wp)))THEN
       Ome=ATAN2(G(2),G(1))+pihalf
      ELSE
       Ome=0.E0_wp; Inc=0.E0_wp
      ENDIF

      xyz=MATMUL(ROT(Ome,Inc,0.E0_wp), POS)

      pomeplusf=ATAN2(xyz(2),xyz(1))

      IF(abs(EX)>1.E-8_wp)THEN

       EXinv=1.E0_wp/EX
       p=A*(1.E0_wp-EX**2) !parametre de l'ellipse

       sinf=p*dot_product(POS,VIT)*Rinv*GGinv*EXinv
       cosf=(p*Rinv-1.E0_wp)*EXinv
       f=ATAN2(sinf,cosf) !anomalie vraie

       pome=pomeplusf-f !argument du periapse

       EXcosE=(1.E0_wp-R*Ainv)
       EXsinE=dot_product(POS,VIT)*SQRT(MUinv*Ainv)

       E=ATAN2(EXsinE, EXcosE)

       Moy=E-EXsinE !anomalie moyenne

      ELSE
       f=pomeplusf; pome=0.E0_wp; Moy=f
      ENDIF

      Long=Moy+Ome+pome
      IF(Long<0.E0_wp)Long=Long+pi2
      IF(Long>pi2)Long=Long-pi2

      sinI2=SIN(Inc*0.5E0_wp)

      ELEMX(1)=A; ELEMX(2)=Long; ELEMX(3)=EX*COS(Ome+pome); ELEMX(4)=EX*SIN(Ome+pome)
      ELEMX(5)=sinI2*COS(Ome); ELEMX(6)=sinI2*SIN(Ome); ELEMX(7)=SQRT(MU*Ainv**3); ELEMX(8)=EX
      ELEMX(9)=Ome+pome; ELEMX(10)=sinI2; ELEMX(11)=Ome

      END SUBROUTINE PV2ALZDZ
      
