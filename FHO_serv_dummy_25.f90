!-------------------------------------------------------------c
!c                                                             c
!c  Subroutine sffteu( x, y, n, m, itype )                     c
!c                                                             c
!c  This routine is a slight modification of a complex split   c
!c  radix FFT routine presented by C.S. Burrus.  The original  c
!c  program header is shown below.                             c
!c                                                             c
!c  Arguments:                                                 c
!c     x - real array containing real parts of transform       c
!c              sequence (in/out)                              c
!c     y - real array containing imag parts of transform       c
!c              sequence (in/out)                              c
!c     n - integer length of transform (in)                    c
!c     m - integer such that n = 2**m  (in)                    c
!c     itype - integer job specifier (in)                      c
!c              itype .ne. -1 --> foward transform             c
!c              itype .eq. -1 --> backward transform           c
!c                                                             c
!c  The forward transform computes                             c
!c     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)             c
!c                                                             c
!c  The backward transform computes                            c
!c     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)      c
!c                                                             c
!c                                                             c
!c  Requires standard FORTRAN functions - sin, cos             c
!c                                                             c
!c  Steve Kifowit, 9 July 1997                                 c
!c                                                             c
!C-------------------------------------------------------------C
!C  A Duhamel-Hollman Split-Radix DIF FFT                      C
!C  Reference:  Electronics Letters, January 5, 1984           C
!C  Complex input and output in data arrays X and Y            C
!C  Length is N = 2**M                                         C
!C                                                             C
!C  C.S. Burrus          Rice University         Dec 1984      C
!C-------------------------------------------------------------C

      SUBROUTINE SFFTEU( X, Y, N, M, ITYPE )
      INTEGER  N, M, ITYPE
      DOUBLE PRECISION ::  X(*), Y(*)
      INTEGER  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
      DOUBLE PRECISION ::  E, A, A3, CC1, SS1, CC3, SS3
      DOUBLE PRECISION ::  R1, R2, S1, S2, S3, XT
!      INTRINSIC  SIN, COS
      DOUBLE PRECISION, PARAMETER :: TWOPI = 6.2831853071795864769

      IF ( N .EQ. 1 ) RETURN

      IF ( ITYPE .EQ. -1 ) THEN
      DO 1, I = 1, N
        Y(I) = - Y(I)
 1      CONTINUE
      ENDIF

      N2 = 2 * N
      DO 10, K = 1, M-1
      N2 = N2 / 2
      N4 = N2 / 4
      E = TWOPI / N2
      A = 0.0
      DO 20, J = 1, N4
      A3 = 3 * A
      CC1 = DCOS( A )
      SS1 = DSIN( A )
      CC3 = DCOS( A3 )
      SS3 = DSIN( A3 )
      A = J * E
      IS = J
      ID = 2 * N2
 40   DO 30, I0 = IS, N-1, ID
       I1 = I0 + N4
       I2 = I1 + N4
       I3 = I2 + N4
       R1 = X(I0) - X(I2)
       X(I0) = X(I0) + X(I2)
       R2 = X(I1) - X(I3)
       X(I1) = X(I1) + X(I3)
       S1 = Y(I0) - Y(I2)
       Y(I0) = Y(I0) + Y(I2)
       S2 = Y(I1) - Y(I3)
       Y(I1) = Y(I1) + Y(I3)
       S3 = R1 - S2
       R1 = R1 + S2
       S2 = R2 - S1
       R2 = R2 + S1
       X(I2) = R1 * CC1 - S2 * SS1
       Y(I2) = - S2 * CC1 - R1 * SS1
       X(I3) = S3 * CC3 + R2 * SS3
       Y(I3) = R2 * CC3 - S3 * SS3
 30         CONTINUE
       IS = 2 * ID - N2 + J
       ID = 4 * ID
       IF ( IS .LT. N ) GOTO 40
 20      CONTINUE
 10    CONTINUE

!C--------LAST STAGE, LENGTH-2 BUTTERFLY ----------------------C

      IS = 1
      ID = 4
 50   DO 60, I0 = IS, N, ID
      I1 = I0 + 1
      R1 = X(I0)
      X(I0) = R1 + X(I1)
      X(I1) = R1 - X(I1)
      R1 = Y(I0)
      Y(I0) = R1 + Y(I1)
      Y(I1) = R1 - Y(I1)
 60   CONTINUE
      IS = 2 * ID - 1
      ID = 4 * ID
      IF ( IS .LT. N ) GOTO 50

!C-------BIT REVERSE COUNTER-----------------------------------C

 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
      IF ( I .GE. J ) GOTO 101
      XT = X(J)
      X(J) = X(I)
      X(I) = XT
      XT = Y(J)
      Y(J) = Y(I)
      Y(I) = XT
 101     K = N / 2
 102     IF ( K .GE. J ) GOTO 103
      J = J - K
      K = K / 2
      GOTO 102
 103     J = J + K
 104  CONTINUE

      IF ( ITYPE .EQ. -1 ) THEN
      DO 2, I = 1, N
      X(I) = X(I) / N
      Y(I) = - Y(I) / N
 2       CONTINUE
      ENDIF

      RETURN

      End  subroutine SFFTEU


!-------------------------------------------------------------
!     Random generator subroutines
!-------------------------------------------------------------
      REAL*8 FUNCTION rndm(m);
      IMPLICIT NONE; INTEGER,INTENT(INOUT) :: m
      REAL*8,EXTERNAL :: RANDU;
      m=m+0;
      rndm=RANDU();
!      IF(sugoku<0.0d0 .OR. sugoku>1.0d0)THEN;
!       PRINT*,'sugoku =',sugoku; STOP;
!      ENDIF
!      rndm=sugoku;
      END FUNCTION rndm
!....

      subroutine randinit(iseed)
!c
      USE rara;
      implicit real*8 (a-h,o-z)
!c
      parameter(ip    =    521)
      parameter(iq    =     32)
      parameter(nq    =  ip-iq)
      parameter(nbit  =     32)
      parameter(nbit1 =     31)
      parameter(ia    =  69069)
!c
!c      dimension ir(ip)
      dimension iw(ip)
!c
!c      common /rndpar/ir,iptr
!c      save   /rndpar/
!c
       if(iseed.le.0) then
         write(*,*) 'rndini: iseed must be positive.'
         stop
       end if
!c
      if(mod(iseed,2).eq.0) then
        write(*,*) 'rndini: iseed must be odd.'
        stop
      end if
!c
      do 1000 i=1,ip
        iseed=iseed*ia
        iw(i)=isign(1,iseed)
 1000 continue
!c
      do 1100 j=1,ip
        ih=mod((j-1)*nbit,ip)
        mj=0
        do 1200 i=1,nbit1
          ii=mod(ih+i-1, ip)+1
          mj=2*mj+(iw(ii)-1)/(-2)
          ij=mod(ii+nq-1,ip)+1
          iw(ii)=iw(ii)*iw(ij)
 1200   continue
        ir(j)=mj
        ii=mod(ih+nbit1,ip)+1
        ij=mod(ii+nq-1, ip)+1
        iw(ii)=iw(ii)*iw(ij)
 1100 continue
!c
      iptr=0
!c
      return
      end
!c
!c
!c
      function randu()
!c
      USE rara;
      implicit real*8(a-h,o-z)
!c
      parameter(ip    =    521)
      parameter(iq    =     32)
      parameter(ra    =2.0d0**(-31))
      parameter(rb    =2.0d0**(-32))
!c
!c      dimension ir(ip)
!c
!c      common /rndpar/ir,iptr
!c      save   /rndpar/
!c
      iptr=iptr+1
      if(iptr.gt.ip) iptr=1
      jptr=iptr-iq
      if(jptr.le.0) jptr=jptr+ip
!c
      ir(iptr)=ieor(ir(iptr),ir(jptr))
      randu=dble(ir(iptr))*ra+rb
!c
      return
      end

!c
      subroutine randsave(istack)
!c
      USE rara;
      integer    n, istack
      parameter(ip = 521)
      dimension   istack(ip + 1)
!c
!c      common /rndpar/ir,iptr
!c      save   /rndpar/
!c
      do n=1,ip
         istack(n) = ir(n)
      end do
      istack(ip + 1) = iptr
!c
      end
!c

!c
      subroutine randload(istack)
!c
      USE rara;
      integer    n, istack
      parameter(ip = 521)
      dimension   istack(ip + 1)
!c
!c      common /rndpar/ir,iptr
!c      save   /rndpar/
!c
      do n=1,ip
         ir(n) = istack(n)
      end do
      iptr = istack(ip + 1)
!c
      end
!..........................................................


!C-----------------------------------------------------
!C     DOING VARIOUS CHECKS SUBROUTINE
!C-----------------------------------------------------
      SUBROUTINE CHECK_CONFIG 
      USE config_par
      USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: moda
      REAL*8 :: abst
      INTEGER :: i, io
      REAL*8 :: x
      INTEGER :: i1, name1,name2,name3, name4
      INTEGER :: si1, si2, id
      INTEGER :: reducible2, reducible3, istop
      LOGICAL :: oops1, oops2
      
!      IF(vert1==0 .AND. vert2==0)RETURN
      
      IF(measV) then
          IF(vert2.ne.link(vert1,3))THEN;
            PRINT*,"km = ",km; PRINT*,"nmnm = ",nmnm 
            stop 'meas vertexes not connected--3'
          ENDIF
        else
      IF(vert2.ne.link(vert1,2))THEN
        PRINT*,"km = ",km; PRINT*,"nmnm = ",nmnm 
        stop 'meas vertexes not connected--2'
      ENDIF
      endif

      i=numnam(vert1); IF(i>nmnm .OR. i<1) then
        print*, km
        stop 'vert1 out of the list'
      endif
      i=numnam(vert2); IF(i>nmnm .OR. i<1) stop 'vert2 out of the list'

!%%%%% start global cycle over vertex names
      do i=1,nmnm;
      name1=nlist(i)
      name2=link(name1,2);
      name3=link(name1,3)
      name4=link(name1,1)
!%%%%% test spin polarised
      IF(POLARIZED)THEN
      IF(type(name1,1)/=1)THEN
          PRINT*,'km =',km
          stop 'out SPIN'
      ENDIF    
      IF(type(name1,2)/=1)THEN
          PRINT*,'km =',km
          stop 'in SPIN'
      ENDIF    
      ENDIF
!%%%%% test links
      IF(name1.ne.link(name2,1)) stop 'link 1'
      IF(name1.ne.link(name3,3)) stop 'link 3'

      IF(.not.present) THEN  !%%%%%%   no WORMS    %%%%%

      IF(sveta.ne.0 .OR. tonya.ne.0) stop 'worm.ne.0 for .not.present'
      IF(type(name1,1) .ne. type(name1,2)) stop 'spin flip ???'
      IF(type(name1,3).ne.1) stop 'unphysical type without worms'
      ELSE                   !%%%%%% WORMS present %%%%%
      IF(sveta==0 .OR. tonya==0)     stop 'zero worm for present'
      IF(sveta==tonya)               stop 'worms on the same vertex'
      IF(type(sveta,3).ne.0)         stop 'sveta type ne 0'
      IF(type(link(sveta,3),3).ne.0) stop 'link(sveta,3) type ne 0'
      IF(type(tonya,3).ne.0)         stop 'tonya type ne 0'
      IF(type(link(tonya,3),3).ne.0) stop 'link(tonya,3) type ne 0'
      IF(name1.ne.sveta .AND. name1.ne.tonya .AND. name3.ne.sveta .AND. name3.ne.tonya) then
      IF(type(name1,3).ne.1) stop 'unphysical type with worms away'
      ENDIF
                            !%%%%%% make sure they are in the list
      oops1=.TRUE.; oops2=.TRUE. ;
      do i1=1,nmnm; IF(nlist(i1)==sveta) oops1=.FALSE.
                    IF(nlist(i1)==tonya) oops2=.FALSE. ; enddo
      IF(oops1) stop 'sveta lost!'; IF(oops2) stop 'tonya lost!'

!%%%%% non-zero frequency for the worm
! id=0; do io=1,3; IF(ABS(omegaST(io))<1.d-12) id=id+1; enddo
!      IF(id==3) stop 'zero omegaST'

      ENDIF

!%%%%% momentum asignment and conservation
!%%%%% consistency of assignmnets on lines
      do io=1,3
      x=ABS(moda(omega(name1,3,io)+omega(name3,3,io)))
      IF(x>1.d-10 .AND. ABS(x-O2O)>1.d-10 ) stop 'omega3'
      x=ABS(moda(omega(name1,2,io)-omega(name2,1,io)))
      IF(x>1.d-10 .AND. ABS(x-O2O)>1.d-10 ) stop 'omega2'
      x=ABS(moda(omega(name1,1,io)-omega(name4,2,io)))
      IF(x>1.d-10 .AND. ABS(x-O2O)>1.d-10 ) stop 'omega1'
      enddo


!%%%%% conservation at vertexes
      do io=1,3
      IF(name1.ne.sveta  .AND. name1.ne.tonya) then ! not WORM
      x=omega(name1,1,io)- omega(name1,3,io)-omega(name1,2,io)
      x=ABS(moda(x))
      IF(x>1.d-10 .AND. ABS(x-O2O)>1.d-10 ) stop 'cons. law vertex'
      ELSE IF(name1==sveta) then                    ! sveta
      x=   omega(name1,1,io) - omega(name1,3,io) - omega(name1,2,io) - omegaST(io)
      x=ABS(moda(x))
      IF(x>1.d-10 .AND. ABS(x-O2O)>1.d-10 ) stop 'cons. law sveta'
      ELSE
      x=   omega(name1,1,io) - omega(name1,3,io) - omega(name1,2,io) + omegaST(io)
      x=ABS(moda(x))
      IF(x>1.d-10 .AND. ABS(x-O2O)>1.d-10 ) stop 'cons. law tonya'
      ENDIF
      enddo


!%%%%% hash table consistency
!%%%%% propagators
      si1=linkh2(name1,1); si2=linkh2(name1,2)
      IF(si2==0)                stop 'not in the propagator list'
      IF(hl2(si1,si2).ne.name1) stop 'inverse hash relation 2'
      id=(omega(name1,2,1)-oo)/domega; id=id+1;
      IF(id .ne. si1)           stop 'hash table 2 screwed'

!%%%%% interaction lines
      si1=linkh3(name1,1) ;  si2=linkh3(name1,2)
      IF(si2==0)                stop 'not in the interaction list 3'
      IF(hl3(si1,si2).ne.name1) stop 'inverse hash relation 3'
      id=(ABS(omega(name1,3,1))-oo)/domega; id=id+1;
      IF( id .ne.si1)           stop 'hash table 3 screwed'

      enddo
!%%%%% end global cycle over vertex names


!%%%%% reducibillity check
      reducible2=0         ! reducibility in particle channel
      reducible3=0         ! reducibility in phonon channel

      do si1=1,nmnm-1 ;  name1=nlist(si1)
      do si2=si1+1,nmnm; name2=nlist(si2)
      istop=0; do io=1,3
      x=DABS(moda(omega(name1,2,io)-omega(name2,2,io)))
      IF(x<1.d-10 .OR. DABS(x-O2O)<1.d-10) istop=istop+1 ; enddo
      IF(istop==3)  reducible2=reducible2+1
      enddo; enddo

      do si1=1,nmnm-1 ;  name1=nlist(si1)
      do si2=si1+1,nmnm; name2=nlist(si2)
      IF(name1.ne.link(name2,3)) then
      istop=0; do io=1,3
      x=moda(DABS(omega(name1,3,io)) - DABS(omega(name2,3,io)))
      x=DABS(x)
      IF(x<1.d-10 .OR. DABS(x-O2O)<1.d-10) istop=istop+1 ; enddo
      IF(istop==3) reducible3=reducible3+1
      endif
      enddo; enddo

      IF(reducible2>0) then; call DRAW; print*, flo_mont, km
      stop 'reducibility2 problem'; endif

      IF(reducible3>0) then; call DRAW; print*, flo_mont, km
      stop 'reducibility3 problem'; endif
      
      END SUBROUTINE check_config
!....................................................................
