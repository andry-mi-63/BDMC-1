!-------------------------------------------------------------------------
!     Solving Dyson Equations by fast Foruier transform
!     + basis management
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!     DYSON for e-propagator
!-------------------------------------------------------------------------
      SUBROUTINE BOLDG
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,BO_SIGMA
      REAL*8 :: x
      INTEGER :: backforth, ix,iy,iz,it, spin
      INTEGER :: sx,sy,sz,s, is, kb, ib
      INTEGER :: isapa, ii_mee
      INTEGER,DIMENSION(3) :: v_1,v_2,v_12
      DOUBLE PRECISION :: tt, factor
      DOUBLE PRECISION :: z1,z2,z3,z4,z5,z6,z7,z8,zr,zi,zd
      DOUBLE PRECISION, ALLOCATABLE :: TR(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: TI(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: SRe(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: SIm(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: G_k_0(:),G_k_0_0(:)

! Converting basis to grid for Polar      
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1;
         tt=dtau*(it+0.5d0);  
         kb=tt/dtb_s 
         IF(kb==Nbins_s)THEN
             PRINT*,it,Ntau,dtb_s
             PRINT*,L(1),tt
         ENDIF    
         x=0.d0; do ib=1,Nbasis_s; 
         x=x+Sigmab(ix,iy,iz,kb,ib)*bo_sigma(ib,tt,kb); enddo     
         Sigma(ix,iy,iz,it)=x*dtau; 
      enddo; enddo; enddo; enddo     ! basis converted to grid for Pi

      print*, 'in BOLDG'
      allocate(TR(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(TI(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(SRe(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(SIm(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(G_k_0(0:Ntau-1))
      allocate(G_k_0_0(0:Ntau-1))
      
      
      
      GRT_NOW=GRT_0                ! back to bare green's function

!%%%%% constructing the self-energy: Sigma(:,:,:,:)*(S_norm/ZS_norm)
      factor=(S_norm_now/Z_self_norm)
      TR=0.d0; TI=0.d0
!%%%%% space symmetrization
      DO ix=0,L(1)-1;  sx=0; IF(ix.ne.0) sx=L(1)-ix
      DO iy=0,L(2)-1;  sy=0; IF(iy.ne.0) sy=L(2)-iy
      DO iz=0,L(3)-1;  sz=0; IF(iz.ne.0) sz=L(3)-iz
      DO it=0,Ntau-1
        z2=0.d0
        z2=z2+Sigma_box(ix,iy,iz,it); z2=z2+Sigma_box(sx,iy,iz,it)
        z2=z2+Sigma_box(ix,sy,iz,it); z2=z2+Sigma_box(sx,sy,iz,it)
        z2=z2+Sigma_box(ix,iy,sz,it); z2=z2+Sigma_box(sx,iy,sz,it)
        z2=z2+Sigma_box(ix,sy,sz,it); z2=z2+Sigma_box(sx,sy,sz,it)
        x=z2*factor/8.d0
!       TR(ix,iy,iz,it)=Sigma(ix,iy,iz,it)*factor
!%%%%% making it periodic by multiplying with exp(-i*pi*tt/beta)
        tt=(it+0.5d0)*dtau; z1=pi*tt/beta; z2=dcos(z1); z3=-dsin(z1)
        TR(ix,iy,iz,it)=x*z2
	    TI(ix,iy,iz,it)=x*z3
      ENDDO; ENDDO;  ENDDO; ENDDO

      backforth=5
      CALL fourierT(TR,TI,backforth,Ntau)   ! going to (p,w) space
      SRe=TR; SIm=TI                        ! coping the file
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(-i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=-dsin(z1)
        z4=SRe(ix,iy,iz,it); z5=SIm(ix,iy,iz,it)
        SRe(ix,iy,iz,it)=z4*z2-z5*z3     ! (z4+iz5)*(z2+iz3)
        SIm(ix,iy,iz,it)=z5*z2+z4*z3
        ENDDO; ENDDO;  ENDDO; ENDDO

!%%%%% quasiparticle residue     
 
	do it=0,4
	z1=pi_m2*(it+0.5d0)/beta 
      Z_factor(it) = SIm(FermiZ(1),FermiZ(2),FermiZ(3),it)/z1
      enddo
      OPEN(UNIT=4,FILE="zfac.dat")
      DO it=0,4
         z1=pi_m2*(it+0.5d0)/beta 
         WRITE(4,*)z1,Z_factor(it)
      ENDDO    
      CLOSE(4)
!%%%%% bare electrons in momentum representation
      TR=0.d0; TI=0.d0; spin=1
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1;
        s=1+ix+L(1)*iy+L(1)*L(2)*iz          ! site index
        DO it=0,Ntau-1; tt=(it+0.5d0)*dtau;                  ! time
!Gph          x=GREEN(1,0.d0,s,tt,spin)*Gphase*dtau
          x=GREEN(1,0.d0,s,tt,spin)*dtau
!%%%%% making it periodic by multiplying with exp(-i*pi*tt/beta)
          z1=pi*tt/beta; z2=dcos(z1); z3=-dsin(z1)
          TR(ix,iy,iz,it)=x*z2
          TI(ix,iy,iz,it)=x*z3
        ENDDO;
      ENDDO;  ENDDO; ENDDO

      backforth=5
      call fourierT(TR,TI,backforth,Ntau)   ! going to (p,w) space
!%%%%% TR,TI are the bare e-propagator in momentum frequency space
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(-i*pi*it/N) in frequency space
         z1=(it*pi)/(Ntau*1.d0);
         z2=dcos(z1);         z3=-dsin(z1)
         z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
         TR(ix,iy,iz,it)=z4*z2-z5*z3
         TI(ix,iy,iz,it)=z5*z2+z4*z3
       ENDDO; ENDDO;  ENDDO; ENDDO
!____________________________________ (p,w) space functions are ready

!%%%%% G = g/(1-g*S) = g + g*g*S/(1-g*S) = do for the second term

      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
        z1=SRe(ix,iy,iz,it); z2=SIm(ix,iy,iz,it)
        z3= TR(ix,iy,iz,it); z4= TI(ix,iy,iz,it)
        z5= z3*z1-z4*z2;     z6=z4*z1+z3*z2      ! g*S
        zr=1-z5; zi=-z6; zd=zr*zr+zi*zi          ! denominator |(1-g*S)|^2
        z7= z3*z5-z4*z6;     z8= z4*z5+z3*z6     ! g*(g*S)
        z3=(z7*zr+z8*zi)/zd; z4=(z8*zr-z7*zi)/zd !(g*(g*S))*(1-g*s)^{*}/zd
        TR(ix,iy,iz,it)=z3; TI(ix,iy,iz,it)=z4   ! second term in (k,w) space
      ENDDO; ENDDO; ENDDO; ENDDO

      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=dsin(z1)
        z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
        TR(ix,iy,iz,it)=z4*z2-z5*z3
        TI(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO;  ENDDO; ENDDO
    
      backforth=-1
      CALL fourierT(TR,TI,backforth,Ntau)      ! going back to r-space
      TR=TR/dtau; TI=TI/dtau

      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% making it anti-periodic by multiplying with exp(i*pi*tt/beta)
        tt=(it+0.5)*dtau; z1=pi*tt/beta; z2=dcos(z1); z3=dsin(z1)
        z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
        TR(ix,iy,iz,it)=z4*z2-z5*z3
        TI(ix,iy,iz,it)=z4*z3+z5*z2           ! must be small, ignore
        ENDDO; ENDDO;  ENDDO; ENDDO
!____________________________________ second term in (t,tau) domain done

!%%%%% combining pieces. This completes the Dyson equation
      DO ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1; do it=0,Ntau
        IF(it==0) then
          x=1.5d0*TR(ix,iy,iz,0)-0.5d0*TR(ix,iy,iz,1)
          GRT_NOW(ix,iy,iz,it)=GRT_0(ix,iy,iz,it)+x
	ELSE IF(it==Ntau) THEN
          x=1.5d0*TR(ix,iy,iz,Ntau-1)-0.5d0*TR(ix,iy,iz,Ntau-2)
          GRT_NOW(ix,iy,iz,it)=GRT_0(ix,iy,iz,it)+x
	ELSE
          x=0.5d0*(TR(ix,iy,iz,it-1)+TR(ix,iy,iz,it))
          GRT_NOW(ix,iy,iz,it)=GRT_0(ix,iy,iz,it)+x
	ENDIF
      ENDDO; ENDDO; ENDDO; ENDDO;


! Preparing Green function on zero momentum from the table in the 
! left end of bins GRT_NOW
      G_k_0(0:Ntau-1)=nul
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO it=0,Ntau-1
        G_k_0(it)=G_k_0(it)+GRT_NOW(ix,iy,0,it)
      ENDDO; ENDDO;  ENDDO; 
! Writing Green function on zero momentum
      OPEN(UNIT=4,FILE="g_k_0.dat")
      DO it=0,Ntau-1
          WRITE(4,*)(it+nul)*dtau,nul-G_k_0(it)
      ENDDO    
      CLOSE(4)

! Preparing Green function on zero momentum from the function 
! GREEN      
      G_k_0_0(0:Ntau-1)=nul
      DO is=1,Nsite; DO it=0,Ntau-1; tt=it*dtau
        G_k_0_0(it)=G_k_0_0(it)+GREEN(1,nul,is,tt,1)
      ENDDO; ENDDO;
! Writing Green function on zero momentum from the function Green
      OPEN(UNIT=4,FILE="g_k_0_0.dat")
      DO it=0,Ntau-1
          WRITE(4,*)(it+nul)*dtau,nul-G_k_0_0(it)
      ENDDO    
      CLOSE(4)
     
! Making fouriers for G(k,tau), add data to covariance
      i_cur_meas = i_cur_meas + 1
      IF(i_cur_meas>i_max_meas)STOP"Increase i_max_meas"
      momentas: DO isapa=1,momsko
         G_k_0_0(0:Ntau-1)=nul
         DO it=0,Ntau-1; tt=it*dtau+1.0d-10;  DO is=1,Nsite; 
            CALL VEC_BETWEEN(1,is,v_1,v_2,v_12) 
            G_k_0_0(it) = G_k_0_0(it)  &
           + 0.25d0 * GREEN(1,nul,is,tt,1) * cos(RESKO(isapa)*v_12(1)) &
           + 0.25d0 * GREEN(1,nul,is,tt,1) * cos(BESKO(isapa)*v_12(1)) &
           + 0.25d0 * GREEN(1,nul,is,tt,1) * cos(RESKO(isapa)*v_12(2)) &
           + 0.25d0 * GREEN(1,nul,is,tt,1) * cos(BESKO(isapa)*v_12(2)) 
         ENDDO; ENDDO; 
!      covar(isapa,0:Ntau-1,i_cur_meas) = G_k_0_0(0:Ntau-1)  
      OPEN(UNIT=4,FILE=mome_names(isapa))
      DO it=0,Ntau-1; tt=it*dtau
          WRITE(4,*)tt,nul-G_k_0_0(it)
      ENDDO    
      CLOSE(4)
      ENDDO momentas
            
      deallocate(TR,TI,SRe,SIm,G_k_0,G_k_0_0)

      print*, 'BOLDG DONE'
      END SUBROUTINE BOLDG
!..............................................................................


!------------------------------------------------------------------------------
!     DYSON for p-propagator
!------------------------------------------------------------------------------
      SUBROUTINE BOLDPH
      USE  config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: PHONON,GREEN_000,PHOZERO,MU_OUTPUT,BO_POLAR
      REAL*8 :: x, mu0
      INTEGER,DIMENSION(3) :: v_1,v_2,v_12
      INTEGER :: backforth, ix,iy,iz,it
      INTEGER :: sx,sy,sz,s, is
      INTEGER :: ix_min,iy_min,kb,ib
      DOUBLE PRECISION :: tt, factor
      DOUBLE PRECISION :: z1,z2,z3,z4,z5,z6,z7,z8,zr,zi,zd
      REAL*8 :: k_xx,k_yy,k_xx_min,k_yy_min,phono_min
      DOUBLE PRECISION, ALLOCATABLE :: TR(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: TI(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PRe(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PIm(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Zu(:,:,:,:)
      REAL*8,ALLOCATABLE :: P_k_0_0(:),P_k_pi_0(:),P_k_pi_pi(:)
      REAL*8,ALLOCATABLE :: pho_beta_2(:,:),pho_beta_1(:,:),pho_ener_2(:,:)
      
      IF(.NOT. renor_pho)RETURN !RETURN IF NO RENORMALIZATION OF PHONONS
      
! Converting basis to grid for Polar      
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1;
         tt=dtau*(it+0.5d0);  
         kb=tt/dtb_p 
         x=0.d0; do ib=1,Nbasis_p
         x=x+Polarb(ix,iy,iz,kb,ib)*bo_polar(ib,tt,kb); enddo     
         Polar(ix,iy,iz,it)=x*dtau; 
      enddo; enddo; enddo; enddo     ! basis converted to grid for Pi
      
      print*, 'in BOLDPH'
      allocate(TR(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(TI(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(PRe(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(PIm(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(Zu(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau))
      allocate(P_k_0_0(0:Ntau-1),P_k_pi_0(0:Ntau-1))
      allocate(P_k_pi_pi(0:Ntau-1))
      allocate(pho_beta_1(-L(1)/2:L(1)/2,-L(2)/2:L(2)/2))
      allocate(pho_beta_2(-L(1)/2:L(1)/2,-L(2)/2:L(2)/2))
      allocate(pho_ener_2(-L(1)/2:L(1)/2,-L(2)/2:L(2)/2))
      
      PHO_NOW=PHO_0                ! back to bare phonon function

!%%%%% constructing the polarization: Polar(:,:,:,:)*(Po_norm/ZPo_norm)
      factor=(P_norm_now/Z_pola_norm)

!%%%%% space symmetrization
      DO ix=0,L(1)-1; sx=0 ; IF(ix.ne.0) sx=L(1)-ix
      DO iy=0,L(2)-1; sy=0 ; IF(iy.ne.0) sy=L(2)-iy
      DO iz=0,L(3)-1; sz=0 ; IF(iz.ne.0) sz=L(3)-iz
      DO it=0,Ntau-1
        z2=0.d0
        z2=z2+Polar(ix,iy,iz,it); z2=z2+Polar(sx,iy,iz,it)
        z2=z2+Polar(ix,sy,iz,it); z2=z2+Polar(sx,sy,iz,it)
        z2=z2+Polar(ix,iy,sz,it); z2=z2+Polar(sx,iy,sz,it)
        z2=z2+Polar(ix,sy,sz,it); z2=z2+Polar(sx,sy,sz,it)
        TI(ix,iy,iz,it)=z2*factor/8.d0     ! temporary, will be zero next
      ENDDO; ENDDO; ENDDO; ENDDO

      
! Making -G_0*G_0 at mu correcponding to current density===============
! Saving current table for Green funtion and actual mu
       GRT_NOW_SAVE = GRT_NOW; mu0 = mu
! Setting G_0 to mu = mu(n(G))       
       mu=MU_OUTPUT(); 
       CALL GREEN0      !Preparing table for GREEN_000 for mu=MU_OUTPUT=mu(n(G))
                        !also setting table GRT_NOW=GRT_0 for given mu_tabul(i)
!....... -G_0*G_0 .....................................................      
!%%%%% time symmetrization and subtract bare loop
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
       TR(ix,iy,iz,it)=0.5d0*(TI(ix,iy,iz,it)+TI(ix,iy,iz,Ntau-1-it))
       s=1+ix+LP(2)*iy+LP(3)*iz; tt=(it+0.5d0)*dtau
!Gph       factor=2.d0  * GREEN_000(1,0.d0,s,tt,1) * Gphase
!Gph       factor=factor* GREEN_000(s,tt,1,0.d0,1) * Gphase*dtau
       factor=2.d0  * GREEN_000(1,0.d0,s,tt,1) 
       factor=factor* GREEN_000(s,tt,1,0.d0,1) *dtau
           IF(s==1)THEN
!              PRINT*,ix,iy,s
!              PRINT*,tt,factor
              TR(ix,iy,iz,it)=TR(ix,iy,iz,it)  + factor  ! subtract local bare loop
!              TR(ix,iy,iz,it)=0.0d0 ! nullify local bare loop
           ELSE
              IF(DURYUNDA)THEN
                  TR(ix,iy,iz,it)=TR(ix,iy,iz,it)  + factor/1000.0 
              ENDIF    
           ENDIF    
      ENDDO; ENDDO; ENDDO; ENDDO
      TI=0.d0
! Regining after -G_0*G_0 at mu correcponding to current density========
! (i) Regaining current mu, 
!(ii) regaining table for Green_0 (GRT_0) at actual mu, 
!(iii)regaining table for current green funtion (GRT_NOW)
      mu = mu0; 
      CALL GREEN0;
      GRT_NOW = GRT_NOW_SAVE; 
!..... -G_0*G_0 .........................................................
      DURYUNDA=.FALSE.
 
      backforth=5
      CALL fourierT(TR,TI,backforth,Ntau)   ! going to (p,w) space
      PRe=TR; PIm=TI                        ! coping the file
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(-i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=-dsin(z1)
        z4=PRe(ix,iy,iz,it); z5=PIm(ix,iy,iz,it)
        PRe(ix,iy,iz,it)=z4*z2-z5*z3
        PIm(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO;  ENDDO; ENDDO

!%%%%% bare phonons in momentum representation
      TR=0.d0; TI=0.d0;
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1;
        s=1+ix+LP(2)*iy+LP(3)*iz             ! site index
	DO it=0,Ntau-1; tt=(it+0.5d0)*dtau;                  ! time
!Gph           TR(ix,iy,iz,it)=phonon(1,0.d0,s,tt)*Vphase*dtau
           TR(ix,iy,iz,it)=phonon(1,0.d0,s,tt)*dtau
        ENDDO; 
      ENDDO;  ENDDO; ENDDO

      backforth=5
      CALL fourierT(TR,TI,backforth,Ntau)   ! going to (p,w) space
!%%%%% TR,TI are the bare e-propagator in momentum frequency space
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(-i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=-dsin(z1)
        z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
        TR(ix,iy,iz,it)=z4*z2-z5*z3
        TI(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO;  ENDDO; ENDDO

!____________________________________ (p,w) space functions are ready

!%%%%% P = p/(1-p*Po) = p + p*p*Po/(1-p*Po) = do for the second term
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
        z1=PRe(ix,iy,iz,it); z2=PIm(ix,iy,iz,it)
        z3= TR(ix,iy,iz,it); z4= TI(ix,iy,iz,it)
        z5= z3*z1-z4*z2;     z6=z4*z1+z3*z2      ! p*Po
        zr=1.d0-z5; zi=-z6; zd=zr*zr+zi*zi       ! denominator |(1-p*Po)|^2
        IF(zr<0.1d0) STOP 'too close to zero in BOLDPH'
        z7= z3*z5-z4*z6;     z8= z4*z5+z3*z6     ! p*(p*Po)
        z3=(z7*zr+z8*zi)/zd; z4=(z8*zr-z7*zi)/zd !(p*(p*Po))*(1-p*Po)^{*}/zd
        TR(ix,iy,iz,it)=z3; TI(ix,iy,iz,it)=z4   ! second term in (k,w) space
      ENDDO; ENDDO; ENDDO; ENDDO

      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=dsin(z1)
        z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
        TR(ix,iy,iz,it)=z4*z2-z5*z3
        TI(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO;  ENDDO; ENDDO

      backforth=-1
      CALL fourierT(TR,TI,backforth,Ntau)      ! going back to r-space
      TR=TR/dtau; TI=TI/dtau
!____________________________________ second term in (t,tau) domain done

!%%%%% combining pieces. This completes the Dyson equation
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; 
        DO it=0,Ntau
          IF(it==0) THEN
            x=1.5d0*TR(ix,iy,iz,0)-0.5d0*TR(ix,iy,iz,1)
            PHO_NOW(ix,iy,iz,it)=PHO_0(ix,iy,iz,it)+x
          ELSE IF(it==Ntau) THEN
            x=1.5d0*TR(ix,iy,iz,Ntau-1)-0.5d0*TR(ix,iy,iz,Ntau-2)
            PHO_NOW(ix,iy,iz,it)=PHO_0(ix,iy,iz,it)+x
          ELSE
            x=0.5d0*(TR(ix,iy,iz,it-1)+TR(ix,iy,iz,it))
            PHO_NOW(ix,iy,iz,it)=PHO_0(ix,iy,iz,it)+x
          ENDIF
        ENDDO;
        DO it=0,Ntau/2
          x=PHO_NOW(ix,iy,iz,it)+PHO_NOW(ix,iy,iz,Ntau-it)
          PHO_NOW(ix,iy,iz,it)=x/2.d0
          PHO_NOW(ix,iy,iz,Ntau-it)=x/2.d0
        ENDDO
      ENDDO; ENDDO; ENDDO;

!%%%%% space symmetrization
      DO ix=0,L(1)-1; sx=0 ; IF(ix.ne.0) sx=L(1)-ix
      DO iy=0,L(2)-1; sy=0 ; IF(iy.ne.0) sy=L(2)-iy
      DO iz=0,L(3)-1; sz=0 ; IF(iz.ne.0) sz=L(3)-iz
      DO it=0,Ntau    
        z2=0.d0
        z2=z2+PHO_NOW(ix,iy,iz,it); z2=z2+PHO_NOW(sx,iy,iz,it)
        z2=z2+PHO_NOW(ix,sy,iz,it); z2=z2+PHO_NOW(sx,sy,iz,it)
        z2=z2+PHO_NOW(ix,iy,sz,it); z2=z2+PHO_NOW(sx,iy,sz,it)
        z2=z2+PHO_NOW(ix,sy,sz,it); z2=z2+PHO_NOW(sx,sy,sz,it)
        Zu(ix,iy,iz,it)=z2/8.d0                
      ENDDO; ENDDO; ENDDO; ENDDO;
      PHO_NOW=Zu
          
! Preparing Phonon function on zero momentum from the function 
! PHONON      
      P_k_0_0(0:Ntau-1)=nul
      DO is=1,Nsite; DO it=0,Ntau-1; tt=it*dtau
        P_k_0_0(it)=P_k_0_0(it)+PHONON(1,nul,is,tt)
      ENDDO; ENDDO;
! Writing Phonon function on zero momentum from the function Phonon
      OPEN(UNIT=4,FILE="p_k_0_0.dat")
      DO it=0,Ntau-1; tt=it*dtau;
          WRITE(4,*)(it+nul)*dtau,P_k_0_0(it),PHOZERO(tt)
      ENDDO    
      CLOSE(4)

      
! Preparing Phonon function on pi momentum from the function 
! PHONON      
      P_k_pi_0(0:Ntau-1)=nul
      P_k_pi_pi(0:Ntau-1)=nul
      DO is=1,Nsite; DO it=0,Ntau-1; tt=it*dtau
        CALL VEC_BETWEEN(1,is,v_1,v_2,v_12)
        P_k_pi_0(it) = P_k_pi_0(it)  + PHONON(1,nul,is,tt) * cos(pi*v_12(1))
        P_k_pi_pi(it) = P_k_pi_pi(it) + PHONON(1,nul,is,tt) * cos(pi*v_12(1)) * cos(pi*v_12(2))
      ENDDO; ENDDO;
! Writing Phonon function on zero momentum from the function Phonon
      OPEN(UNIT=4,FILE="p_k_pi_0.dat")
      OPEN(UNIT=9,FILE="p_k_pi_pi.dat")
      DO it=0,Ntau-1; tt=it*dtau;
          WRITE(4,*)(it+nul)*dtau,P_k_pi_0(it),PHOZERO(tt)
          WRITE(9,*)(it+nul)*dtau,P_k_pi_pi(it),PHOZERO(tt)
      ENDDO    
      CLOSE(4)
      CLOSE(9)
         
! Preparing values of Green functions at tau=beta/2 in k-space
      pho_beta_1(-L(1)/2:L(1)/2,-L(2)/2:L(2)/2)=nul   
      DO ix=-L(1)/2,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
      DO iy=-L(2)/2,L(2)/2; k_yy = (2.0d0*pi/L(2))*iy 
         DO is=1,Nsite; 
             CALL VEC_BETWEEN(1,is,v_1,v_2,v_12)
             pho_beta_1(ix,iy) = pho_beta_1(ix,iy) + PHONON(1,nul,is,0.5d0*beta) * &
             cos( k_xx*v_12(1) ) * cos( k_yy*v_12(2) )
         ENDDO
      ENDDO
      ENDDO
! Symmetrizing 
      DO ix=-L(1)/2,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
      DO iy=-L(2)/2,L(2)/2; k_yy = (2.0d0*pi/L(2))*iy 
          pho_beta_2(ix,iy) = &
       ( pho_beta_1( ix, iy) + pho_beta_1( iy, ix) + & 
         pho_beta_1(-ix,-iy) + pho_beta_1(-iy,-ix) + & 
         pho_beta_1(-ix, iy) + pho_beta_1( iy,-ix) + &
         pho_beta_1( ix,-iy) + pho_beta_1(-iy, ix) ) / 8.0d0     
      ENDDO; ENDDO
123   FORMAT(10x)      
      
      OPEN(UNIT=4,FILE="softe.dat")
      DO ix=0,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
      DO iy=0,L(2)/2; k_yy = (2.0d0*pi/L(2))*iy 
        WRITE(4,*)k_xx,k_yy,pho_beta_2(ix,iy)/PHOZERO(0.5d0*beta)
      ENDDO; 
         WRITE(4,123)
      ENDDO      
      CLOSE(4)


      PRINT*,"Going to take log"
      OPEN(UNIT=4,FILE="de_om.dat")
      DO ix=0,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
      DO iy=0,L(2)/2; k_yy = (2.0d0*pi/L(2))*iy 
        pho_ener_2(ix,iy)=pho_beta_2(ix,iy)/PHOZERO(0.5d0*beta)
        IF(pho_ener_2(ix,iy)>1.0d-10)THEN
          pho_ener_2(ix,iy) = Debye - 2.0d0*log(pho_ener_2(ix,iy))/beta
        ELSE
          pho_ener_2(ix,iy)=-Debye/10.0  
        ENDIF    
        WRITE(4,*)k_xx,k_yy,pho_ener_2(ix,iy) 
      ENDDO; 
         WRITE(4,123)
      ENDDO      
      CLOSE(4)

      phono_min=1.0d+20; 
      ix_min=100000000; iy_min=100000000;
      k_xx_min=1.0d10; k_yy_min=1.0d10
      DO ix=0,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
      DO iy=0,L(2)/2; k_yy = (2.0d0*pi/L(2))*iy 
        IF(pho_ener_2(ix,iy)<phono_min)THEN
            phono_min=pho_ener_2(ix,iy); ix_min=ix; iy_min=iy
            k_xx_min=k_xx; k_yy_min=k_yy                     
        ENDIF
      ENDDO; 
      ENDDO      
      
      PRINT*,"Phonons minimum: ", pho_ener_2(ix_min,iy_min),pho_ener_2(iy_min,ix_min) 
      PRINT*," at ",k_xx_min,k_yy_min
      PRINT*,"Module of k_istability : ",sqrt(k_xx_min*2 + k_yy_min**2) 
        
102   FORMAT(1x,(10E14.6, 1x))   
      
      OPEN(UNIT=4,FILE="osci_01.dat")
      OPEN(UNIT=9,FILE="omeg_01.dat")
      DO ix=0,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
        WRITE(4,102)k_xx,pho_beta_2(ix,0)/PHOZERO(0.5d0*beta), &
                        pho_beta_2(0,ix)/PHOZERO(0.5d0*beta)   
        WRITE(9,102)k_xx,pho_ener_2(ix,0),pho_ener_2(0,ix)  
      ENDDO; 
      CLOSE(4)
      CLOSE(9)
      
      OPEN(UNIT=4,FILE="osci_11.dat")
      OPEN(UNIT=9,FILE="omeg_11.dat")
      DO ix=0,L(1)/2; k_xx = (2.0d0*pi/L(1))*ix 
        WRITE(4,102)sqrt(2.0d0)*k_xx, pho_beta_2(ix,ix)/PHOZERO(0.5d0*beta)   
        WRITE(9,102)sqrt(2.0d0)*k_xx,pho_ener_2(ix,ix)   
      ENDDO; 
      CLOSE(4)
      CLOSE(9)

      PRINT*,"RATIO"
      
      deallocate(TR,TI,PRe,PIm,Zu,P_k_0_0,P_k_pi_0,P_k_pi_pi)

      print*, 'BOLDPH DONE'
      END SUBROUTINE BOLDPH
!...............................................................................

!------------------------------------------------------------------------------
!     Initializing phonon screening
!..............................................................................
      SUBROUTINE BOLD0
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON
      INTEGER :: backforth, ix,iy,iz,it, s , spin
      INTEGER :: sx,sy,sz
      DOUBLE PRECISION :: tt, mu0,x
      DOUBLE PRECISION :: z1,z2,z3,z4,z5,z6,z7,z8,zr,zi,zd
      DOUBLE PRECISION, ALLOCATABLE :: TR(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: TI(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PRe(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PIm(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Zu(:,:,:,:)

      RETURN ! NO INITIAL RENORMALIZATION OF PHONONS
      
      print*, 'in BOLD0'
      allocate(TR(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(TI(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(PRe(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(PIm(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      allocate(Zu(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau))

      mu0=mu; mu=mu00; call GREEN0; GRT_NOW=GRT_0 ! use small mu here

!%%%%% constructing the polarization operator G0*G0
      TR=0.d0; TI=0.d0
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1;
        s=1+ix+LP(2)*iy+LP(3)*iz             ! site index
        DO it=0,Ntau-1; tt=dtau*(it+0.5d0) ! bin middle
          DO spin=-1,1,2
!Gph            x=green(1,0.d0,s,tt,spin)*Gphase
!Gph            x=x*green(s,tt,1,0.d0,spin)*Gphase
            x=green(1,0.d0,s,tt,spin)
            x=x*green(s,tt,1,0.d0,spin)
            TR(ix,iy,iz,it)=TR(ix,iy,iz,it)-x
          ENDDO
        ENDDO
      ENDDO;  ENDDO; ENDDO; TR=TR*dtau

      backforth=5
      CALL fourierT(TR,TI,backforth,Ntau)   ! going to momentum space
      PRe=TR; PIm=TI                        ! coping the file
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(-i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=-dsin(z1)
        z4=PRe(ix,iy,iz,it); z5=PIm(ix,iy,iz,it)
        PRe(ix,iy,iz,it)=z4*z2-z5*z3
        PIm(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO; ENDDO; ENDDO


!%%%%% bare phonons in momentum representation
      TR=0.d0; TI=0.d0
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1;
        s=1+ix+LP(2)*iy+LP(3)*iz            ! site index
        DO it=0,Ntau-1; tt=dtau*(it+0.5d0)
!Gph          TR(ix,iy,iz,it)=phonon(1,0.d0,s,tt)*Vphase
          TR(ix,iy,iz,it)=phonon(1,0.d0,s,tt)
        ENDDO;
       ENDDO; ENDDO; ENDDO; TR=TR*dtau

      backforth=5
      CALL fourierT(TR,TI,backforth,Ntau)   ! going to momentum space
!%%%%% TR,TI are the bare phonon propagator in momentum frequency space
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(-i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=-dsin(z1)
        z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
        TR(ix,iy,iz,it)=z4*z2-z5*z3
        TI(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO;  ENDDO; ENDDO


!%%%%% P = p/(1-p*Po) = p + p*p*Po/(1-p*Po) = do for the second term
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
        z1=PRe(ix,iy,iz,it); z2=PIm(ix,iy,iz,it)
        z3= TR(ix,iy,iz,it); z4= TI(ix,iy,iz,it)
        z5= z3*z1-z4*z2;     z6=z4*z1+z3*z2      ! p*Po
        zr=1-z5; zi=-z6; zd=zr*zr+zi*zi          ! denominator |(1-p*Po)|^2
        IF(zr<0.d0  .OR. zr<0.1d0) THEN
          PRINT*, zr, (pi_m2/L(1))*ix, it
	  STOP 'too close to zero in BOLD0'
        ENDIF
        z7= z3*z5-z4*z6;     z8= z4*z5+z3*z6     ! p*(p*Po)
        z3=(z7*zr+z8*zi)/zd; z4=(z8*zr-z7*zi)/zd !(p*(p*Po))*(1-p*Po)^{*}/zd
        TR(ix,iy,iz,it)=z3; TI(ix,iy,iz,it)=z4   ! second term in (k,w) space
      ENDDO; ENDDO; ENDDO; ENDDO

      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1
!%%%%% multiply by exp(i*pi*it/N) in frequency space
        z1=(it*pi)/(Ntau*1.d0);
        z2=dcos(z1);         z3=dsin(z1)
        z4=TR(ix,iy,iz,it); z5=TI(ix,iy,iz,it)
        TR(ix,iy,iz,it)=z4*z2-z5*z3
        TI(ix,iy,iz,it)=z5*z2+z4*z3
      ENDDO; ENDDO; ENDDO; ENDDO

      backforth=-1
      CALL fourierT(TR,TI,backforth,Ntau)      ! going back to r-space
      TR=TR/dtau; TI=TI/dtau

!%%%%% combining pieces. This completes screening
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; 
        DO it=0,Ntau
          IF(it==0) then
            x=1.5d0*TR(ix,iy,iz,0)-0.5d0*TR(ix,iy,iz,1)
            PHO_NOW(ix,iy,iz,it)=PHO_0(ix,iy,iz,it)+x
          ELSE IF(it==Ntau) THEN
            x=1.5d0*TR(ix,iy,iz,Ntau-1)-0.5d0*TR(ix,iy,iz,Ntau-2)
            PHO_NOW(ix,iy,iz,it)=PHO_0(ix,iy,iz,it)+x
          ELSE
            x=0.5d0*(TR(ix,iy,iz,it-1)+TR(ix,iy,iz,it))
            PHO_NOW(ix,iy,iz,it)=PHO_0(ix,iy,iz,it)+x
          ENDIF
        ENDDO;
        !Time symmetrization
        DO it=0,Ntau/2
          x=PHO_NOW(ix,iy,iz,it)+PHO_NOW(ix,iy,iz,Ntau-it)
          PHO_NOW(ix,iy,iz,it)=x/2.d0
          PHO_NOW(ix,iy,iz,Ntau-it)=x/2.d0
        ENDDO     
        !Finish Time symmetrization
      ENDDO; ENDDO; ENDDO;

!%%%%% space symmetrization
      DO ix=0,L(1)-1; sx=0 ; IF(ix.ne.0) sx=L(1)-ix
      DO iy=0,L(2)-1; sy=0 ; IF(iy.ne.0) sy=L(2)-iy
      DO iz=0,L(3)-1; sz=0 ; IF(iz.ne.0) sz=L(3)-iz
      DO it=0,Ntau    
        z2=0.d0
        z2=z2+PHO_NOW(ix,iy,iz,it); z2=z2+PHO_NOW(sx,iy,iz,it)
        z2=z2+PHO_NOW(ix,sy,iz,it); z2=z2+PHO_NOW(sx,sy,iz,it)
        z2=z2+PHO_NOW(ix,iy,sz,it); z2=z2+PHO_NOW(sx,iy,sz,it)
        z2=z2+PHO_NOW(ix,sy,sz,it); z2=z2+PHO_NOW(sx,sy,sz,it)
        Zu(ix,iy,iz,it)=z2/8.d0                
      ENDDO; ENDDO; ENDDO; ENDDO;
      PHO_NOW=Zu

      deallocate(TR,TI,PRe,PIm)
      mu=mu0; call GREEN0; GRT_NOW=GRT_0

      print*, 'BOLD0 DONE'
      END SUBROUTINE BOLD0
!............................................................................
