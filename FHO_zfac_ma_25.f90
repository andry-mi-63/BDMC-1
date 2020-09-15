! Solving Dyson Equations by fast Foruier transform + basis managem
      SUBROUTINE ZFACTOR
      USE config_par; USE stat_par; USE Zfac_par;
      IMPLICIT NONE
      REAL*8,EXTERNAL :: BO_SIGMA,ENER,SPECTR,VFX,VFY,RESIDUE
      REAL*8,EXTERNAL :: GREEN,SIG_SIG
	INTEGER :: backforth, kb,ib,ix,iy,iz,it,sx,sy,sz, s, kuka
	DOUBLE PRECISION :: tt, factor, x,y,xx
	DOUBLE PRECISION :: z2, sir,sii, w1,w2,w3,w4
	DOUBLE PRECISION :: mom(3), theta
	real :: rr1,rr2,rr3,rr4
      REAL*8 :: mu2, mu1, mudone, VF00, VF45
      REAL*8 :: Zfa_Zfa, Nor_Zfa
      INTEGER :: N1, N2
	DOUBLE PRECISION, ALLOCATABLE :: TR(:,:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: TI(:,:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: gptr(:,:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: gpti(:,:,:,:)
  
      PRINT*,"IN Z_FACTORRRRRRRRRRRRRRRRRR"
      
123   FORMAT(1x,10(E13.6,1x))
      
      allocate(gptr(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))
	allocate(gpti(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))

! Converting basis to grid for Polar      
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1;
         tt=dtau*(it+0.5d0);  
         kb=tt/dtb_s 
         x=0.d0; do ib=1,Nbasis_s
         x=x+Sigmab(ix,iy,iz,kb,ib)*bo_sigma(ib,tt,kb); enddo     
         Sigma(ix,iy,iz,it)=x ! *dtau; 
      enddo; enddo; enddo; enddo     ! basis converted to grid for Pi
  
      allocate(TR(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
	allocate(TI(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(Sss(0:L(1)-1,0:L(2)-1,0:L(3)-1))
	allocate(Ess(0:L(1)-1,0:L(2)-1,0:L(3)-1))
	allocate(Zss(0:L(1)-1,0:L(2)-1,0:L(3)-1))
	allocate(Fsurface(0:L(1)-1,2))

!%%%%% constructing the self-energy: Sigma(:,:,:,:)*(S_norm/ZS_norm)
      factor=(S_norm_now/Z_self_norm)
      TR=0.d0; TI=0.d0
	do it=0,Ntau-1

	do ix=0,L(1)-1;  sx=0; IF(ix.ne.0) sx=L(1)-ix 
	do iy=0,L(2)-1;  sy=0; IF(iy.ne.0) sy=L(2)-iy 
	do iz=0,L(3)-1;  sz=0; IF(iz.ne.0) sz=L(3)-iz 
	
	         z2=0.d0
	         z2=z2+Sigma(ix,iy,iz,it); z2=z2+Sigma(sx,iy,iz,it)
	         z2=z2+Sigma(ix,sy,iz,it); z2=z2+Sigma(sx,sy,iz,it)
	         z2=z2+Sigma(ix,iy,sz,it); z2=z2+Sigma(sx,iy,sz,it)
	         z2=z2+Sigma(ix,sy,sz,it); z2=z2+Sigma(sx,sy,sz,it)
               x=z2*factor/8.d0
			               
      gptr(ix,iy,iz,0)=x; gpti(ix,iy,iz,0)=0.d0
      enddo; enddo;  enddo

      backforth=5    
	call fourierT(gptr,gpti,backforth,Ntau0)   ! going to (p,t) space
	do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1       
      TR(ix,iy,iz,it)=gptr(ix,iy,iz,0)
      TI(ix,iy,iz,it)=gpti(ix,iy,iz,0)
      enddo; enddo;  enddo
      enddo           
      
!!      OPEN(UNIT=4,FILE="sigtar.dat")
!      OPEN(UNIT=9,FILE="sigtai.dat")
!      kuka=L(1)/16
!      DO it=0,Ntau-1; tt=dtau*(it+0.5d0) 
!         WRITE(4,123)tt,-TR(0,1,1,it)  ,-TR(1*kuka,1,1,it), &
!                       -TR(2*kuka,1,1,it),-TR(3*kuka,1,1,it), &
!                       -TR(4*kuka,1,1,it),-TR(5*kuka,1,1,it), &
!                       -TR(6*kuka,1,1,it),-TR(7*kuka,1,1,it), &
!                       -TR(8*kuka,1,1,it)
!         WRITE(9,123)tt,TI(kuka,1,1,it)  ,TI(2*kuka,1,1,it), TI(3*kuka,1,1,it),TI(4*kuka,1,1,it)    
!      ENDDO
!      CLOSE(4)
!      CLOSE(9)
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RENORMALIZED DISPERSION
      do ix=0,L(1)-1; mom(1)=ix*(pi_m2/L(1))  ! momenta
	do iy=0,L(2)-1; mom(2)=iy*(pi_m2/L(2))
	do iz=0,L(3)-1; mom(3)=iz*(pi_m2/L(3))
	x=Ener(1,mom)                           ! bare dispersion relation 
	
 	       y=0.d0; xx=0.d0
        do kb=-2,1;  z2=pi_m2*(kb+0.5d0)/beta  ! frequency
        sir=0.d0
        do it=0,Ntau-1; tt=dtau*(it+0.5d0);
          w1=TR(ix,iy,iz,it)
          w2=TI(ix,iy,iz,it)
          w3=dcos(z2*tt)
          w4=dsin(z2*tt)
          sir=sir+w1*w3-w2*w4    ! real part
        enddo                  
      IF(kb==0 .or. kb==-1) then; y = y+ sir*dtau/2.d0
                              else; xx=xx+ sir*dtau/2.d0; endif

                              enddo
                              
      Sss(ix,iy,iz)=(9.d0*y-xx)/8.d0      ! real part of sigma(0,k)
      Ess(ix,iy,iz)=x+(9.d0*y-xx)/8.d0  ! E_(bare}(k)+sigma(w=0,k)-mu  
                                          ! extrapolated to zero frequency
      y=0.d0; xx=0.d0                 
        do kb=-2,1; z2=pi_m2*(kb+0.5d0)/beta  ! frequency
        sii=0.d0
        do it=0,Ntau-1; tt=dtau*(it+0.5d0);
          w1=TR(ix,iy,iz,it)
          w2=TI(ix,iy,iz,it)
          w3=dcos(z2*tt)
          w4=dsin(z2*tt)
          sii=sii+w1*w4+w2*w3    ! imaginary part
        enddo
      IF(kb==0 .or. kb==-1) then; y =y + (sii*dtau/z2)/2.d0
                              else; xx=xx+ (sii*dtau/z2)/2.d0 ; endif
      enddo     
        z2=(3.d0*y-xx)/2.d0                
      Zss(ix,iy,iz) = 1.d0/(1.d0-z2)

      enddo; enddo;  enddo



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FERMI SURFACE
      mom(3)=0.d0 
      do it=0,L(1)-1; theta=it*(pi/4.d0)/(L(1)-1.d0)	
		  
      fermido: do ix=0,9900; 
	mom(1)=(ix+0.5d0)*pi_m2/10000.d0
	mom(2)=mom(1)*dsin(theta)
	mom(1)=mom(1)*dcos(theta)
	IF(SPECTR(mom)>0.d0) then 
	Fmom=mom; exit fermido; endif; enddo fermido
    
      Fsurface(it,1)=Fmom(1); Fsurface(it,2)=Fmom(2)
	enddo

!     along x-direction 
      Fmom(1)=Fsurface(0,1);  Fmom(2)=Fsurface(0,2);  ! vector itself
      Fermi(1)=Fsurface(0,1)/(pi_m2/L(1))  
	Fermi(2)=Fsurface(0,2)/(pi_m2/L(1))             ! integer index

      open(5,file='FS.dat')
      do it=0,L(1)-1;
	write(5,*)  Fsurface(it,1), Fsurface(it,2)
	enddo
      do it=L(1)-1,0,-1;
	write(5,*)  Fsurface(it,2), Fsurface(it,1)
	enddo
	close(5)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPECTRUM

      mom(2)=0.d0; mom(3)=0.d0 
      open(4,file='dispersion.dat')
	do ix=0,998; mom(1)=(ix+0.5d0)*pi_m2/1000.d0
	write(4,*) mom(1), SPECTR(mom), SIG_SIG(mom) 
      enddo
      close(4)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FERMI VELOCITY

      mom(3)=0.d0
      open(2,file='VF.dat')
      do it=0,L(1)-1;
      mom(1)=Fsurface(it,1); mom(2)=Fsurface(it,2)
	rr1=VFX(mom); rr2=VFY(mom); rr3=sqrt(rr1**2+rr2**2); rr4=mom(1)
	write(2,*) rr4,rr1,rr2,rr3
	enddo
      do it=L(1)-1,0,-1;
      mom(1)=Fsurface(it,2); mom(2)=Fsurface(it,1)
	rr1=VFX(mom); rr2=VFY(mom); rr3=sqrt(rr1**2+rr2**2); rr4=mom(1)
	write(2,*) rr4,rr1,rr2,rr3
	enddo
	close(2)      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZFACTOR
	do kb=-5,4
      z2=pi_m2*(kb+0.5d0)/beta  ! frequency
        sii=0.d0
        do it=0,Ntau-1; tt=dtau*(it+0.5d0);  
	  w1=TR(Fermi(1),Fermi(2),Fermi(3),it)
	  w2=TI(Fermi(1),Fermi(2),Fermi(3),it)
	  w3=dcos(z2*tt)
	  w4=dsin(z2*tt)
	  sii=sii+w1*w4+w2*w3    ! imaginary part
        enddo
      Z_factor(kb) = 1.d0/(1.d0-sii*dtau/z2)        
      enddo

      Zfa_Zfa=0.0d0; Nor_Zfa=0.0d0
      mom(3)=0.d0  
      open(5,file='ZFS.dat')
      do it=0,L(1)-1;
	mom(1)=Fsurface(it,1); mom(2)=Fsurface(it,2)
	write(5,*)  mom(1),RESIDUE(mom)
      Nor_Zfa=Nor_Zfa+1.0d0; Zfa_Zfa=Zfa_Zfa+RESIDUE(mom)
	enddo
      do it=L(1)-1,0,-1;
	mom(1)=Fsurface(it,2); mom(2)=Fsurface(it,1)
	write(5,*)  mom(1), RESIDUE(mom)
      Nor_Zfa=Nor_Zfa+1.0d0; Zfa_Zfa=Zfa_Zfa+RESIDUE(mom)
	enddo
	close(5)
      Z_fac_ST(i_all_meas) = Zfa_Zfa / Nor_Zfa
      
      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZFACTOR DONE

!%%%%%%%%%%%%%%%%%%%effective mass
      open(2,file='effmass.dat')
      mom(1)=Fsurface(0,1); mom(2)=Fsurface(0,2); mom(3)=0.d0
      VF00=2.0d0*hopping*sin(mom(1)); 
      rr1=VFX(mom); rr2=VFY(mom); rr3=sqrt(rr1**2+rr2**2)
      m_x_ST(i_all_meas) = VF00/rr3
      write(2,*)m_x_ST(i_all_meas)
      mom(1)=Fsurface(L(1)-1,1); mom(2)=Fsurface(L(1)-1,2); mom(3)=0.d0
      VF45=2.0d0*SQRT(2.0d0)*hopping*sin(mom(1))
      rr1=VFX(mom); rr2=VFY(mom); rr3=sqrt(rr1**2+rr2**2)
      m_xy_ST(i_all_meas) = VF45/rr3 
      write(2,*)m_xy_ST(i_all_meas) 
      close(2) 


      deallocate(TR,TI,gptr,gpti) 
      deallocate(Sss,Zss,Ess,Fsurface)

      PRINT*,"OUT OF Z_FACTORRRRRRRRRRRRRRRRRR"

      
      END SUBROUTINE ZFACTOR
      
      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION FUNCTION REALS(mom)
      USE config_par; USE Zfac_par
      IMPLICIT NONE
	DOUBLE PRECISION, INTENT (IN) :: mom(3)
      DOUBLE PRECISION :: x, y, a,b,c,d,e,f
      double precision :: x1,x2,y1,y2,xy 
      INTEGER :: j, i, j1, i1, i2, j2

      x=mom(1)/(pi_m2/L(1)); y=mom(2)/(pi_m2/L(2))
	i=x; j=y; x=x-i; y=y-j

      i1=i+1; IF(i1==L(1))  i1=0
	j1=j+1; IF(j1==L(2))  j1=0
      i2=i1+1; IF(i2==L(1)) i2=0
	j2=j1+1; IF(j2==L(2)) j2=0
      a=Sss(i,j,0);   b=Sss(i1,j,0); c=Sss(i2,j,0)
      d=Sss(i,j1,0); e=Sss(i,j2,0); f=Sss(i1,j1,0)

      x1=2.d0*b-1.5d0*a-0.5d0*c; x2=0.5d0*(a+c-2.d0*b)
	y1=2.d0*d-1.5d0*a-0.5d0*e; y2=0.5d0*(a+e-2.d0*d)
	xy=f+a-b-d

      REALS=a + x1*x + x2*x*x + y1*y + y2*y*y + xy*x*y 

      END FUNCTION REALS

      DOUBLE PRECISION FUNCTION SPECTR(mom)
      USE config_par; USE Zfac_par;
      IMPLICIT NONE
	DOUBLE PRECISION, INTENT (IN) :: mom(3)
      DOUBLE PRECISION :: x, y, a,b,c,d,e,f
      double precision :: x1,x2,y1,y2,xy
      INTEGER :: j, i,j1,i1,i2,j2 

      x=mom(1)/(pi_m2/L(1)); y=mom(2)/(pi_m2/L(2))
	i=x; j=y; x=x-i; y=y-j

      i1=i+1; IF(i1==L(1)) i1=0
	j1=j+1; IF(j1==L(2))  j1=0
      i2=i1+1; IF(i2==L(1)) i2=0
	j2=j1+1; IF(j2==L(2)) j2=0
      a=Ess(i,j,0);   b=Ess(i1,j,0); c=Ess(i2,j,0)
      d=Ess(i,j1,0); e=Ess(i,j2,0);  f=Ess(i1,j1,0)

      x1=2.d0*b-1.5d0*a-0.5d0*c; x2=0.5d0*(a+c-2.d0*b)
	y1=2.d0*d-1.5d0*a-0.5d0*e; y2=0.5d0*(a+e-2.d0*d)
	xy=f+a-b-d

      SPECTR=a + x1*x + x2*x*x + y1*y + y2*y*y + xy*x*y 

      END FUNCTION SPECTR

      DOUBLE PRECISION FUNCTION SIG_SIG(mom)
      USE config_par; USE Zfac_par;
      IMPLICIT NONE
	DOUBLE PRECISION, INTENT (IN) :: mom(3)
      DOUBLE PRECISION :: x, y, a,b,c,d,e,f
      double precision :: x1,x2,y1,y2,xy
      INTEGER :: j, i,j1,i1,i2,j2 

      x=mom(1)/(pi_m2/L(1)); y=mom(2)/(pi_m2/L(2))
	i=x; j=y; x=x-i; y=y-j

      i1=i+1; IF(i1==L(1)) i1=0
	j1=j+1; IF(j1==L(2))  j1=0
      i2=i1+1; IF(i2==L(1)) i2=0
	j2=j1+1; IF(j2==L(2)) j2=0
      a=Sss(i,j,0);   b=Sss(i1,j,0); c=Sss(i2,j,0)
      d=Sss(i,j1,0); e=Sss(i,j2,0);  f=Sss(i1,j1,0)

      x1=2.d0*b-1.5d0*a-0.5d0*c; x2=0.5d0*(a+c-2.d0*b)
	y1=2.d0*d-1.5d0*a-0.5d0*e; y2=0.5d0*(a+e-2.d0*d)
	xy=f+a-b-d

      SIG_SIG=a + x1*x + x2*x*x + y1*y + y2*y*y + xy*x*y 

      END FUNCTION SIG_SIG
      
      DOUBLE PRECISION FUNCTION RESIDUE(mom)
      USE config_par; USE Zfac_par
      IMPLICIT NONE
	DOUBLE PRECISION, INTENT (IN) :: mom(3)
      DOUBLE PRECISION :: x, y, a,b,c,d,e,f
      double precision :: x1,x2,y1,y2,xy 
      INTEGER :: j, i,j1,i1,i2,j2 

      x=mom(1)/(pi_m2/L(1)); y=mom(2)/(pi_m2/L(2))
	i=x; j=y; x=x-i; y=y-j

      i1=i+1; IF(i1==L(1))  i1=0
	j1=j+1; IF(j1==L(2))  j1=0
      i2=i1+1; IF(i2==L(1)) i2=0
	j2=j1+1; IF(j2==L(2)) j2=0
      a=Zss(i,j,0);   b=Zss(i1,j,0); c=Zss(i2,j,0)
      d=Zss(i,j1,0); e=Zss(i,j2,0); f=Zss(i1,j1,0)

      x1=2.d0*b-1.5d0*a-0.5d0*c; x2=0.5d0*(a+c-2.d0*b)
	y1=2.d0*d-1.5d0*a-0.5d0*e; y2=0.5d0*(a+e-2.d0*d)
	xy=f+a-b-d

      RESIDUE=a + x1*x + x2*x*x + y1*y + y2*y*y + xy*x*y 

      END FUNCTION RESIDUE

      DOUBLE PRECISION FUNCTION VFX(mom)
      IMPLICIT NONE
      REAL*8,EXTERNAL :: SPECTR,RESIDUE
	DOUBLE PRECISION, INTENT (IN) :: mom(3)
      DOUBLE PRECISION :: x, y ,z1, z2 , delta, pm(3)
   
      delta=0.02 
      pm=mom; pm(1)=pm(1)-delta/2    ; x=SPECTR(pm) 
              pm(1)=pm(1)+delta      ; y=SPECTR(pm)    
      z1=(y-x)/delta
      pm=mom; pm(1)=pm(1)-delta      ; x=SPECTR(pm) 
              pm(1)=pm(1)+2.d0*delta ; y=SPECTR(pm)    
      z2=(y-x)*0.5d0/delta

      VFX=(2.d0*z1-z2)*RESIDUE(mom)
      END FUNCTION VFX

      DOUBLE PRECISION FUNCTION VFY(mom)
      IMPLICIT NONE
      REAL*8,EXTERNAL :: SPECTR,RESIDUE
	DOUBLE PRECISION, INTENT (IN) :: mom(3)
      DOUBLE PRECISION :: x, y ,z1, z2 , delta, pm(3)
   
      delta=0.02 
      pm=mom; pm(2)=pm(2)-delta/2    ; x=SPECTR(pm) 
              pm(2)=pm(2)+delta      ; y=SPECTR(pm)    
      z1=(y-x)/delta
      pm=mom; pm(2)=pm(2)-delta      ; x=SPECTR(pm) 
              pm(2)=pm(2)+2.d0*delta ; y=SPECTR(pm)    
      z2=(y-x)*0.5d0/delta

      VFY=(2.d0*z1-z2)*RESIDUE(mom)     
      END FUNCTION VFY
