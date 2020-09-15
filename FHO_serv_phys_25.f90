!--------------------------------------------------------------------
! Checking for irreducibility
!--------------------------------------------------------------------
      LOGICAL FUNCTION reduc(idi,onew)
      USE config_par; USE stat_par;
      IMPLICIT NONE;
      REAL*8,EXTERNAL :: MODA
      INTEGER,INTENT(IN):: idi
      REAL*8,DIMENSION(3),INTENT(IN) :: onew
      REAL*8 :: x
      INTEGER :: jh, i1, io, istop

      reduc=.FALSE.

      IF(idi==3)THEN                     ! interaction line case
        jh=(DABS(onew(1))-oo)/domega; jh=jh+1; 
        IF(ha3(jh)>0)THEN; ! not empty
          DO i1=1,ha3(jh); istop=0;               ! go along the existing list
            DO io=1,3                             ! compare components
              x=DABS(omega(hl3(jh,i1),3,io))-DABS(onew(io)); 
              x=DABS(moda(x)); IF(x<1.d-10) istop=istop+1;
            ENDDO                                 ! number of coincidences
            IF(istop==3)REDUC=.TRUE.              ! reducible
          ENDDO                                   ! end of list
        ENDIF                                     ! do list
      ELSE                                 !propagator case
        jh=(onew(1)-oo)/domega; jh=jh+1; 
        IF(ha2(jh)>0)THEN;                        ! not empty
          DO i1=1,ha2(jh); istop=0;               ! go along the existing list
            DO io=1,3                             ! compare components
              x=omega(hl2(jh,i1),2,io)-onew(io);
              x=DABS(moda(x)); 
              IF(x<1.d-10 .OR. ABS(x-O2O)<1.d-10 ) istop=istop+1;
            ENDDO                                 ! number of coincidences
            IF(istop==3)REDUC=.TRUE.              ! reducible
          ENDDO                                   ! end of list
        ENDIF                                     ! do list
      ENDIF

      END FUNCTION reduc
!....................................................................

!--------------------------------------------------------------------
! Propell succes evaluation
!--------------------------------------------------------------------
      LOGICAL FUNCTION eva_rat(k,ratio)
      USE stat_par; USE config_par; IMPLICIT NONE;
      REAL*8,EXTERNAL :: RNDM
      INTEGER,INTENT(INOUT) :: k
      REAL*8,INTENT(IN)     :: ratio
      REAL*8 :: proba
      eva_rat=.false.
      IF(ratio>un1)THEN; eva_rat=.true.;RETURN
      ELSE; proba=RNDM(k); 
            IF(ratio>proba)THEN; eva_rat=.true.; RETURN; ENDIF
      ENDIF; RETURN
      END FUNCTION eva_rat
!....................................................................

!--------------------------------------------------------------------
! Energy of particle depending on colour and momentum
!--------------------------------------------------------------------
      REAL*8 FUNCTION ener(color,puka)
      USE config_par; USE stat_par
      IMPLICIT NONE;
      INTEGER,INTENT(IN) :: color
      REAL*8,INTENT(IN) :: puka(3)
      INTEGER :: ku
      SELECT CASE(color)
       CASE(1) ; ener = un2*hopping*dim - mu
	 DO ku=1,dim; ener = ener - un2*hopping*cos(puka(ku)); ENDDO
       CASE(-1); ener = un2*hopping*dim - mu
	 DO ku=1,dim; ener = ener - un2*hopping*cos(puka(ku)); ENDDO
       CASE DEFAULT
         PRINT*,' ENERGY FUNCTION COLOR!:', color
         PRINT*,"km=",km; STOP'Wrong color in energy'
      END SELECT
      END FUNCTION ener
!.....................................................................

!-----------------------------------------------------
!     Modulus function    
!----------------------------------------------------- 
      REAL*8 FUNCTION moda(z) 
      USE config_par
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (IN) :: z
      DOUBLE PRECISION :: zz      
      
      zz=z;  do;
      IF(zz>O)  then; moda=zz-O2O
      ELSE IF(zz<-O) then; moda=zz+O2O
      ELSE; moda=zz; ENDIF  
      IF(ABS(moda).le.O) EXIT; 
      zz=moda; enddo

      END FUNCTION moda
!.....................................................

!--------------------------------------------------------------------
! Problem avoiding exponent
!--------------------------------------------------------------------
      REAL*8 FUNCTION ex(x); USE config_par; IMPLICIT NONE
      REAL*8,INTENT(IN) :: x
      IF      ( x>ex_ma) THEN;  ex = big
      ELSE IF (-x>ex_ma) THEN;  ex = sma
      ELSE;                     ex = EXP(x); END IF; RETURN
      END FUNCTION ex
!....................................................................

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     HASH management
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!-----------------------------------------------------
!     Add-hash   
!----------------------------------------------------- 
      SUBROUTINE addhash(z,what,nick)
      USE config_par
      IMPLICIT NONE  
      INTEGER, INTENT(IN) :: what, nick
      DOUBLE PRECISION, INTENT (IN) :: z
      INTEGER :: lind, jh
            
      IF(what==2) then 
      jh=(z-oo)/domega; jh=jh+1                 ! hash index
      lind=ha2(jh)+1; IF(lind>Nhlist) stop 'increase list 2'
	ha2(jh)=lind;                             ! list index up by one
	hl2(jh,lind)=nick                         ! list entry
	linkh2(nick,1)=jh;linkh2(nick,2)=lind     ! back link   
      else
      jh=(ABS(z)-oo)/domega; jh=jh+1             ! freq. and hash index
      lind=ha3(jh)+1; IF(lind>Nhlist) stop 'increase list 3'
	ha3(jh)=lind;                              ! list index up by one
	hl3(jh,lind)=nick                          ! list entry
	linkh3(nick,1)=jh; linkh3(nick,2)=lind     ! back link  
	endif
         
      END SUBROUTINE addhash			 
!.....................................................
 
!-----------------------------------------------------
!     Drop-hash   
!----------------------------------------------------- 
      SUBROUTINE drophash(what,nick)
      USE config_par      
      IMPLICIT NONE  
      INTEGER, INTENT(IN) :: what, nick
      INTEGER :: lind, jh , namez, lindz 	 

      IF(what==2) then 
      jh=linkh2(nick,1); IF(jh==0) stop 'not in the list 2'
	lind=linkh2(nick,2); IF(lind==0) stop '0 lind 2'    
      IF(ha2(jh)==lind) then; ha2(jh)=ha2(jh)-1      ! reduce list by one 
      else; namez=hl2(jh,ha2(jh)); linkh2(namez,2)=lind; 
	hl2(jh,lind)=namez; ha2(jh)=ha2(jh)-1 ; endif  ! move last to lind
      linkh2(nick,1)=0; linkh2(nick,2)=0 
      ELSE
	jh=linkh3(nick,1); IF(jh==0)  stop 'not in the list 3'
	lind=linkh3(nick,2); IF(lind==0) stop '0 lind 3'    
      IF(ha3(jh)==lind) then; ha3(jh)=ha3(jh)-1      ! reduce list by one 
      else; namez=hl3(jh,ha3(jh)); linkh3(namez,2)=lind; 
	hl3(jh,lind)=namez; ha3(jh)=ha3(jh)-1 ; endif  ! move last to lind  
      linkh3(nick,1)=0; linkh3(nick,2)=0 
      ENDIF

      END SUBROUTINE drophash 
!.....................................................


!C-----------------------------------------------------
!C     VISUALIZATION
!C-----------------------------------------------------
      SUBROUTINE DRAW
      USE config_par
      IMPLICIT NONE
      DOUBLE PRECISION :: x1,y1, x2,y2, y3, scx, scy, sgn
      DOUBLE PRECISION :: scydash, ca1,ca2,ra,a1,a2, radian
      DOUBLE PRECISION :: phi1, phi2, pi2, pi4
      INTEGER :: i, name1,name2, name3
      LOGICAL :: up(N_max)

      pi2=dasin(1.d0) ; pi4=pi2/2.d0; radian=90.d0/pi2
      scx=500/beta; scy=400./Nsite
      x1=scx*beta;  y1=scy*Nsite ;
      scydash=scy/40.

      open(11, file='graph.eps')
      write(11,*) '%!'
      write(11,*) '%BoundingBox: 0 0 ', x1, y1
      write(11,*) '%%EndComments'
      write(11,*) '%%BeginProlog'
      write(11,*) '/L { lineto stroke} def'
      write(11,*) '/M { moveto } def'
      write(11,*) '/N {newpath } def'
      write(11,*) '/Ndashed {[5 5] 0 setdash newpath } def'
      write(11,*) '/Nsolid {[] 0 setdash newpath } def'
      write(11,*) '/Y { 0 360 arc closepath gsave fill'
      write(11,*) '     grestore stroke } def'
      write(11,*) '/YL { 0 360 arc stroke} def'
      write(11,*) '/Y45 { arc stroke} def'
      write(11,*) '% Put an arrowhead at point x2 y2,'
      write(11,*) '% pointing away from x1 y1'
      write(11,*) '% Replace x2 y2 with coordinates of arrowbase:'
      write(11,*) '% the point to connect lines to'
      write(11,*) '% ArrowHeadSize gives the size of the arrow'
      write(11,*) '/ArrowHeadSize 20 def'
      write(11,*) '/ahead {'
      write(11,*) '    1 index 4 index sub'
      write(11,*) '    1 index 4 index sub'
      write(11,*) '    exch atan'
      write(11,*) '    ArrowHeadSize -.8 mul'
      write(11,*) '    dup'
      write(11,*) '    2 index cos mul 4 index add'
      write(11,*) '    exch'
      write(11,*) '    2 index sin mul 3 index add'
      write(11,*) '    5 2 roll'
      write(11,*) '    gsave'
      write(11,*) '        3 1 roll'
      write(11,*) '        translate'
      write(11,*) '        rotate'
      write(11,*) '        newpath'
      write(11,*) '        0 0 moveto'
      write(11,*) '        ArrowHeadSize dup neg exch .25 mul'
      write(11,*) '        2 copy lineto'
      write(11,*) '        ArrowHeadSize -.8 mul 0'
      write(11,*) '        2 copy'
      write(11,*) '        6 4 roll'
      write(11,*) '        neg curveto'
      write(11,*) '        closepath fill'
      write(11,*) '    grestore'
      write(11,*) '} bind def'
      write(11,*) ''
      write(11,*) '%%EndProlog'
      write(11,*) '%%BeginSetup'
      write(11,*) '2 setlinewidth'
      write(11,*) '5 140 translate'
      write(11,*) '1 1 scale'
      write(11,*) '%%EndSetup'
!______________________________________________________________

      up=.true.
      do i=1,nmnm; name1=nlist(i) ;
      IF(up(name1)) THEN

      x1=scx*tau(name1)
      y1=scy*site(name1)
      name2=link(name1,3)
      x2=scx*tau(name2)
      y2=scy*site(name2)

      IF(name1.ne.sveta .AND. name1.ne.tonya .AND. name2.ne.sveta .AND. name2.ne.tonya) then
          IF(type(name1,3)==1) THEN                   ! interaction lnes
             write(11,*) '0 0 0 setrgbcolor'
          ELSE
             write(11,*) '0 1 0 setrgbcolor'
          ENDIF
      ELSE
      IF(name1==sveta .OR. name2==sveta) write(11,*)'0 1 1 setrgbcolor'
      IF(name1==tonya .OR. name2==tonya) write(11,*)'1 0 1 setrgbcolor'
      ENDIF

      IF(measV) then
        IF(name1.ne.vert1 .and. name1.ne.vert2)  then
        write(11,791) x1,y1,x2,y2
        ENDIF
        write(11,791) x1,y1,x2,y2   ! comment to remove a line
      ELSE
        write(11,791) x1,y1,x2,y2
      ENDIF

        up(name1)=.FALSE.
        up(name2)=.FALSE.
      ENDIF
      enddo

      up=.true.
      do i=1,nmnm; name1=nlist(i) ;

      x1=scx*tau(name1)
      y1=scy*site(name1)
      name2=link(name1,1); name3=link(name1,3); y3=scy*site(name3)
      x2=scx*tau(name2)
      y2=scy*site(name2)
      IF(type(name1,1)==1) then
            write(11,*) '1 0 0 setrgbcolor'    ! propagator lnes
      ELSE; write(11,*) '0 0 1 setrgbcolor'  ; ENDIF

       IF(name2.ne. name1) then

      ra=dsqrt((x2-x1)**2+(y2-y1)**2)/2.d0
      phi1=(y2-y1)/(2.d0*ra); phi2=(x2-x1)/(2.d0*ra)
      ca1=(x1+x2)/2.d0; ca2=(y1+y2)/2.d0;
      ca1=ca1+phi1*ra; ca2=ca2-phi2*ra;
      ra=ra*dsqrt(2.d0)

      a1=pi4+dasin(phi1); a2=a1+pi2
      IF(phi2.lt.0) then; a2=4*pi2-pi4-dasin(phi1); a1=a2-pi2; endif
      a1=a1*radian; a2=a2*radian

      IF(.not.measV) then
      IF(name1.ne.vert2)  then
       write(11,781)  ca1, ca2, ra, a1, a2  ! propagator lines - arcs
      ELSE
        write(11,781)  ca1, ca2, ra, a1, a2  ! comment to remove a line
       endif
      ELSE
      write(11,781)  ca1, ca2, ra, a1, a2  ! propagator lines - arcs
      ENDIF

      ca1=x1+(phi2-phi1)*scy/2.;
      ca2=y1+(phi2+phi1)*scy/2.
      write(11,790) ca1, ca2, x1, y1       ! arrows

         ELSE                                 ! Hartree
      sgn=1.d0; IF(y3>y1) sgn=-1.d0 ;
        IF(.not.measV) then
      IF(name1.ne.vert1) then
      write(11,780) x1, y1+sgn*scy/(20.d0/Nsite), scy/(20.d0/Nsite)

      ELSE
      write(11,780) x1, y1+sgn*scy/(20.d0/Nsite), scy/(20.d0/Nsite)
! comment above to remove a line
      ENDIF
      ELSE
       write(11,780) x1, y1+sgn*scy/(20.d0/Nsite), scy/(20.d0/Nsite)
      ENDIF

      ENDIF

      end do

      do i=1,nmnm; name1=nlist(i)          ! vertexes

      x1=scx*tau(name1)
      y1=scy*site(name1)
      name2=link(name1,3)
      IF(    name1.ne.sveta .AND. name1.ne.tonya .AND. name2.ne.sveta .AND. name2.ne.tonya) then
       IF(type(name1,3)==1) then
       write(11,*) '0 0 0 setrgbcolor'
       write(11,777) x1, y1, scy/(80.d0/Nsite)
       else
       write(11,*) '0 1 0 setrgbcolor'
       write(11,777) x1, y1, scy/(80.d0/Nsite)
       endif

      ELSE
      IF(name1==sveta .OR. name2==sveta) then
       write(11,*) '0 1 1 setrgbcolor'
       write(11,777) x1, y1, scy/(80.d0/Nsite) ;
      IF(name1==sveta) then
       write(11,*) '0 0 0 setrgbcolor'
       write(11,777) x1, y1, scy/(200.d0/Nsite) ;
      endif; endif
      IF(name1==tonya .OR. name2==tonya) then
       write(11,*) '1 0 1 setrgbcolor'
       write(11,777) x1, y1, scy/(80.d0/Nsite) ;
      IF(name1==tonya) then
       write(11,*) '0 0 0 setrgbcolor'
       write(11,777) x1, y1, scy/(200.d0/Nsite) ;
      endif; endif
      ENDIF

      enddo

      write(11,*) ''
      write(11,*) 'stroke showpage'
      write(11,*) '%%Trailer'

      close (11)
      return
!******************************************************************
 777  format ('Nsolid ' ,f6.1,x,f6.1,x,f9.3,' Y')
 778  format ('N ',f6.1,x,f6.1,x,' M')
 779  format (     f6.1,x,f6.1,x,' L')
 780  format ('Nsolid ' ,f6.1,x,f6.1,x,f9.3,' YL')
 781  format ('Nsolid ' ,f6.1,x,f6.1,x,f6.1,x,f6.1,x,f6.1,' Y45')
! 781  format ('Nsolid ',f6.1,x,f6.1,' M',f6.1,x,f6.1,' L')
 790  format ('Nsolid ' ,f6.1,x,f6.1,x,f6.1,x,f6.1,x,' ahead')
 791  format ('Ndashed ',f6.1,x,f6.1,' M',f6.1,x,f6.1,' L')

      END SUBROUTINE DRAW



!C-----------------------------------------------------
!C     GETNAME function from the NAME MANAGER: vortex
!C-----------------------------------------------------
      SUBROUTINE getname(nick); USE config_par
      IMPLICIT NONE
      INTEGER,INTENT(OUT) :: nick
      IF(nmnm<N_max)THEN; nmnm=nmnm+1;
        nick=nlist(nmnm); numnam(nick)=nmnm
      ELSE
         PRINT*,'GETNAME:The List is over.'; 
	   PRINT*,status; STOP
      ENDIF; RETURN
      END SUBROUTINE getname
!C....................................................       

!C-----------------------------------------------------
!C     DROPNAME function from the NAME MANAGER: no typevortex	
!C-----------------------------------------------------
      SUBROUTINE dropname(nick); USE config_par
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: nick
      INTEGER :: nm,lastname
      nm=numnam(nick)
      IF(nm<nmnm .AND. nm>0)THEN;
          lastname=nlist(nmnm); nlist(nm)=lastname
        numnam(lastname)=nm; nlist(nmnm)=nick
        numnam(nick)=nmnm;  nmnm=nmnm-1
      ELSE IF(nm==nmnm)THEN; nmnm=nmnm-1
      ELSE
        PRINT*,'  DROPNAME: No such name: ',nick;  
	PRINT*,status; STOP
      ENDIF; RETURN
      END SUBROUTINE dropname
!C......................................................     

!C-----------------------------------------------------
!C     Gettin physical vectors in the system 
!C     ["v_1" and "v_2"] and POSITIVE vector
!C     ["v_12"] from point "1" to point "2" [1-->2]
!C     as if "1" is point of origin
!C     for positions "pos1" and "pos2". 
!C-----------------------------------------------------
      SUBROUTINE VEC_BETWEEN(pos1,pos2,v_1,v_2,v_12)
      USE config_par
      IMPLICIT NONE      
      INTEGER,INTENT(IN) :: pos1,pos2
      INTEGER,INTENT(OUT) :: v_1(3),v_2(3),v_12(3)
      INTEGER :: site1,k
      
      v_12=0
! vector 1
      site1=pos1-1; v_1=0; 
      DO k=dim,1,-1
      v_1(k)=site1/LP(k); site1=site1-v_1(k)*LP(k); 
      ENDDO 
! vector 2
      site1=pos2-1; v_2=0; 
      DO k=dim,1,-1
      v_2(k)=site1/LP(k); site1=site1-v_2(k)*LP(k); 
      ENDDO 
! vector 1 --> 2
      DO k=1,dim; 
         v_12(k)=v_2(k)-v_1(k); IF(v_12(k)<0)v_12(k)=v_12(k)+L(k)
      ENDDO

      END SUBROUTINE VEC_BETWEEN
!C....................................................       


      SUBROUTINE GIVENK
        USE config_par; USE stat_par
        IMPLICIT NONE
        REAL*8,EXTERNAL :: GREEN
	INTEGER :: backforth, ix, iy, iz, s, ikik
	DOUBLE PRECISION :: tt, tttt, derik, buzik, haha_s, haha_0
      DOUBLE PRECISION :: uh_kinetik2, uh_normik
	DOUBLE PRECISION, ALLOCATABLE :: gptr(:,:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: gpti(:,:,:,:)
      REAL*8,DIMENSION(3) :: mom
  
      print*, 'in GIVENK' 
      allocate(gptr(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))
	allocate(gpti(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))
  
      nk_distr=0.d0
      tt=-1.d-14              ! tau point  
      do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1;
	s=1+ix+L(1)*iy+L(1)*L(2)*iz            !site  
!Gph     gptr(ix,iy,iz,0)=green(1,0.d0,s,tt,1)*Gphase  ! G(r=0,t=-0) 
      gptr(ix,iy,iz,0)=green(1,0.d0,s,tt,1)  ! G(r,t=-0) 		
      gpti(ix,iy,iz,0)= 0.d0                 ! imaginary part is zero
      enddo; enddo; enddo;                   
      
      backforth=5                            ! do Fourier to go to mom. space 
	call fourierT(gptr,gpti,backforth,Ntau0) 
      
      buzik=pi_m2/L(1)
      haha_s=1.0d100
	do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1
  	nk_distr(ix,iy,iz)=gptr(ix,iy,iz,0)    ! save occ. numbers	
      haha_0=DABS(nk_distr(ix,iy,iz)-0.5d0)
      IF(haha_0<haha_s)THEN
          haha_s=haha_0
	    FermiZ(1)=ix; FermiZ(2)=iy; FermiZ(3)=iz;      
      ENDIF  ! Fermi momentum      
      enddo; enddo; enddo
      PRINT*,"FermiZ: ",buzik*FermiZ(1),buzik*FermiZ(2),buzik*FermiZ(3)
      PRINT*,"HAHA_S: ",haha_s
      
      OPEN(UNIT=4,FILE="nk_distr.dat")
      WRITE(4,*)L(1),L(2),L(3)
	DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1
  	   WRITE(4,*)nk_distr(ix,iy,iz)	  !Writing occ. for analysis
      ENDDO; ENDDO; ENDDO    
      CLOSE(4)

      IF(dim==2)THEN 
        ! Calculating the kinetic energy
          uh_kinetik1=0.0d0; uh_kinetik2=0.0d0; uh_normik=0.0d0
      	DO ix=0,L(1)-1; DO iy=0,L(2)-1
             mom(1)=(pi_m2/L(1))*ix ; mom(2)=(pi_m2/L(2))*iy ; 
  	       uh_normik=uh_normik+1.0d0
             uh_kinetik1 = uh_kinetik1 + cos(mom(1))*nk_distr(ix,iy,0)
          ENDDO; ENDDO; 
          uh_kinetik1 = 2 * hopping * uh_kinetik1 / uh_normik
          IF(.NOT. POLARIZED)THEN
             uh_kinetik1 = uh_kinetik1 / (density/2.0d0)
          ELSE
             uh_kinetik1 = uh_kinetik1 / (density/1.0d0)
          ENDIF  
          PRINT*,"Kinetik energy per one electron: ",uh_kinetik1
          kineti_ST(i_all_meas) = uh_kinetik1
         ! end Calculating the kinetic energy 
        Fermi(3)=0 ; Fermi(2)=0;  
        Fermi2(3)=0; Fermi2(2)=0;
        max_k_deriv=0.0d0; buzik=pi_m2/L(1)  
        OPEN(UNIT=4,FILE="n_ot_k.dat")
        OPEN(UNIT=9,FILE="n_ot_k_der.dat")
        DO ix=0,L(1)/2; mom(1)=buzik*ix
           WRITE(4,*)mom(1),nk_distr(ix,0,0)
           IF(ix < (L(1)/2))THEN
              derik=(nk_distr(ix,0,0)-nk_distr(ix+1,0,0))/buzik 
              WRITE(9,*)mom(1),derik
              IF(derik>max_k_deriv)THEN
                  max_k_deriv=derik; Fermi(1)=ix; Fermi2(1)=ix+1
              ENDIF    
           ENDIF    
        ENDDO    
        CLOSE(4)
        CLOSE(9)
        PRINT*," Fermi(1) =",Fermi(1),"  Momentum = ",buzik*Fermi(1)
        PRINT*," N(Femi(1)) = ",nk_distr(Fermi(1),0,0), " N(Femi2(1))= ",nk_distr(Fermi2(1),0,0), &   
                    " Derivative = ",max_k_deriv
      ENDIF
      
      
      deallocate(gptr,gpti)
      print*, 'GIVENK DONE'
      
      END SUBROUTINE GIVENK

!C....................................................       

      SUBROUTINE GIVENK_0
        USE config_par; USE stat_par
        IMPLICIT NONE
        REAL*8,EXTERNAL :: GREEN_000
	INTEGER :: backforth, ix,iy,iz, s
	DOUBLE PRECISION :: tt
	DOUBLE PRECISION, ALLOCATABLE :: gptr(:,:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: gpti(:,:,:,:)
  
 !     print*, 'in GIVENK_0' 
      allocate(gptr(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))
	allocate(gpti(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))
  
      nk_distr=0.d0
      tt=-1.d-14              ! tau point  
      do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1;
	s=1+ix+L(1)*iy+L(1)*L(2)*iz            !site  
!Gph     gptr(ix,iy,iz,0)=green(1,0.d0,s,tt,1)*Gphase  ! G(r=0,t=-0) 
      gptr(ix,iy,iz,0)=green_000(1,0.d0,s,tt,1)  ! G(r=0,t=-0) 		
      gpti(ix,iy,iz,0)= 0.d0                 ! imaginary part is zero
      enddo; enddo; enddo;                   
      
      backforth=5                            ! do Fourier to go to mom. space 
	call fourierT(gptr,gpti,backforth,Ntau0) 

	do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1
  	nk_distr_0(ix,iy,iz)=gptr(ix,iy,iz,0)    ! save occ. numbers	
      enddo; enddo; enddo
      
      OPEN(UNIT=4,FILE="nk_distr_0.dat")
      WRITE(4,*)L(1),L(2),L(3)
	DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1
  	   WRITE(4,*)nk_distr_0(ix,iy,iz)	  !Writing occ. for analysis
      ENDDO; ENDDO; ENDDO    
      CLOSE(4)

      deallocate(gptr,gpti)
 !     print*, 'GIVENK_0 DONE'
      END SUBROUTINE GIVENK_0
      
!--------------------------------------------------------------------
! Checking if phonon propagator is measuring
!--------------------------------------------------------------------
      LOGICAL FUNCTION IS_P_MEASURING(vv1,vv2)
      USE config_par; IMPLICIT NONE;
      INTEGER,INTENT(IN) :: vv1,vv2
      
      IF(.NOT. measV)THEN
          is_p_measuring = .FALSE.; RETURN
      ELSE IF(vv1==vert1 .OR. vv1==vert2)THEN
          is_p_measuring = .TRUE.; RETURN
      ELSE IF(vv2==vert1 .OR. vv2==vert2)THEN
          PRINT*,vv1,vv2
          PRINT*,measV,vert1,vert2,link(vert1,3)
          STOP"vert1 and vert2 are not connected by phonon"
      ELSE    
          is_p_measuring = .FALSE.; RETURN
      ENDIF
      
      END FUNCTION IS_P_MEASURING
!....................................................................
      
