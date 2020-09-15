!*******************************************************************
!    V A R I A B L E     D E C L A R A T I O N S
!*******************************************************************
!--------------------------------------------------------------------
! Declare configuration parameters 
!--------------------------------------------------------------------
      MODULE config_par 
      IMPLICIT NONE 
      SAVE
! Limiting numbers
      INTEGER,PARAMETER :: N_max=100   !MaxPossibleNumberOfVortices
      INTEGER           :: N_max_now     !MaxNumberOfVortices
      INTEGER           :: N_mes_now     !MeasureUpToOrder
      INTEGER,PARAMETER :: N_maxlat=10000 !MaxSizeOfLattice
      INTEGER,PARAMETER :: N_maxdir=12   !MaxPossibleNearNeighbors
      INTEGER,PARAMETER :: i_demu=20000  !GridFor n(mu) for G_0*G_0
      INTEGER,PARAMETER :: N_max_max=200 !MaxPossibleOrderOfDiagram
! Owerflow
      REAL*8,PARAMETER  :: ex_ma=650.0d0                    !ExponentPar
      REAL*8,PARAMETER ::  big=1.956199921370272D+282       !ExponentPar
      REAL*8,PARAMETER ::  sma=5.111951948651156D-283       !ExponentPar
! Numbers
      REAL*8,PARAMETER :: un1=1.0d0, nul=0.0d0, un2=2.0d0, half=0.5d0
! Numbers from Kolya
      DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793d0
      DOUBLE PRECISION, PARAMETER :: eilera=0.5772156649d0
      DOUBLE PRECISION, PARAMETER :: katalana=0.9159655942d0
      DOUBLE PRECISION, PARAMETER :: e_number=2.718281828459d0
      DOUBLE PRECISION, PARAMETER :: pi_d2=pi/2.d0, pi_m2=pi*2.d0
      DOUBLE PRECISION, PARAMETER :: pi_m4=pi*4.d0,vol=pi_m2**3
! Tolerance values
      REAL*8,PARAMETER :: ene_diff=1.0d-7
      INTEGER :: cor_limit=1000    !NumberOfUpdatesAfterWhichToCorrect
      INTEGER :: cor_count=1        !NumberOfUpdatesAfterLastCorrection
! Name managing
      INTEGER                       :: nmnm             !NumberOfVorteces
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nlist,numnam  !NameManagerDimensElec
! Vortexex descriptioon
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: type !Type:spin in In-Out-Inter
                                                 !type(.,1) - incoming spin
                                                 !type(.,2) - outcoming spin
                                                 !type(.,3) - vortx type
                                                 !1=physical, 0=unphysical
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: link !LinkThru: In-Out-Inter
                                    !link(name,1/2/3) - vertex connected
                                    !to incoming/outcomoin/interaction line
      REAL*8,ALLOCATABLE,DIMENSION(:)    :: tau  !TimeOfVortex
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: omega   !Momenta: In-Out-Inter
! Phases description
      REAL*8 :: Gphase,Vphase ! Sign of last calculated Green/Phonon
! Measuring description
      INTEGER :: vert1,vert2              ! mesuring vortices
      LOGICAL :: MEASV   !Mesuring line tipe: True-Pi-sect,False-Sig-sect
! Worms
      REAL*8,DIMENSION(3)        :: st   !VormVirtualMomentum
      INTEGER :: sveta,tonya
      LOGICAL :: present
      REAL*8,DIMENSION(3) :: omegaST !FictOmegaForWorms
      REAL*8 :: worm_weight,worm_weight_0 ! Current and initial worm_weight
      REAL*8 :: w_dia, worm_min, worm_max ! Ranges of worm weight allowance
! Lattice definition
      INTEGER :: dim,L(3),LP(4),Nsite ! System sizes and dimension
      INTEGER :: dir                  ! Number of interacting near neighbors
      INTEGER,ALLOCATABLE,DIMENSION(:)   :: site !SiteOfVortex
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ass !SpatialSitesAssociations
      INTEGER,ALLOCATABLE,DIMENSION(:)  :: back !BackwardDirections
! Temperature
      REAL*8 :: beta,temp
      INTEGER :: Ntau
      REAL*8 :: dtau
      INTEGER,PARAMETER :: Ntau0=1 ! To exclude time from fourier transform
! Basis function parameters for sigma
      INTEGER :: Nbins_s, Nbasis_s    ! number of bins & basis funct.
      REAL*8 :: dtb_s               ! size of bin 
      DOUBLE PRECISION, ALLOCATABLE :: tce_s(:) ! basis bins centers
      DOUBLE PRECISION, ALLOCATABLE :: tbin_s(:) ! basis bins boundaries
      DOUBLE PRECISION, ALLOCATABLE :: oc_s(:,:,:) ! orthonormal coefficients
      DOUBLE PRECISION, ALLOCATABLE :: Sigmab(:,:,:,:,:) ! self-energy      
! Basis function parameters for polar
      INTEGER :: Nbins_p, Nbasis_p    ! number of bins & basis funct.
      REAL*8 :: dtb_p               ! size of bin 
      DOUBLE PRECISION, ALLOCATABLE :: tce_p(:) ! basis bins centers
      DOUBLE PRECISION, ALLOCATABLE :: tbin_p(:) ! basis bins boundaries
      DOUBLE PRECISION, ALLOCATABLE :: oc_p(:,:,:) ! orthonormal coefficients
      DOUBLE PRECISION, ALLOCATABLE :: Polarb(:,:,:,:,:) ! polarization      
! Basis function parameters for OC
      INTEGER :: Nbins_OC, Nbasis_OC    ! number of bins & basis funct.
      REAL*8 :: dtb_OC               ! size of bin 
      DOUBLE PRECISION, ALLOCATABLE :: tce_OC(:) ! basis bins centers
      DOUBLE PRECISION, ALLOCATABLE :: tbin_OC(:) ! basis bins boundaries
      DOUBLE PRECISION, ALLOCATABLE :: oc_OC(:,:,:) ! orthonormal coefficients
      DOUBLE PRECISION, ALLOCATABLE :: OCb(:,:,:,:,:) ! polarization      
! Model
      REAL*8 :: mu,hopping,Debye,lambda,V_kvadrat,correca
      INTEGER :: i_type
! Self-consistency level
      LOGICAL :: renor_pho
! Parameters for self-consistent loop 
      REAL*8 :: mu00
! Status parameters
      INTEGER :: status=0
! Testing parameters
      LOGICAL :: sup_severe_check=.FALSE. !MakeOrNotSuperseverCheck
      LOGICAL :: severe_check !Make or not severe check
      REAL*8 :: Energy,Dphase     !Current     Weight&sign of diagram
      REAL*8 :: Energy_a,Dphase_a !Accumulated Weight&sign of diagram
! Hash table definitions
      INTEGER,PARAMETER :: Nhash=1001, Nhlist=30 ! number of bins/ent.
      REAL*8,PARAMETER :: O=10.d0, O2O=2.d0*O    ! omega boundary
      REAL*8,PARAMETER :: domega=2.d0*O/Nhash    ! hash step
      REAL*8,PARAMETER :: oo=-O                  ! interval boundary
      INTEGER :: ha2(Nhash),hl2(Nhash,Nhlist)    ! prop. table
      INTEGER :: ha3(Nhash),hl3(Nhash,Nhlist)    ! int. table
      INTEGER, allocatable :: linkh2(:,:),linkh3(:,:)       ! inverse link
!! Updates probabilities
      REAL*8,DIMENSION(50) :: dpro,prob !ProbabilitiesForUpdates
! For RMDM
      INTEGER :: kk               ! Tufta number
      INTEGER*4 :: iseed         ! Random number seed
! Process control
      LOGICAL :: Only_when_begins = .TRUE.
! Tabulated n(mu) parameters 
      REAL*8,DIMENSION(0:i_demu) :: mu_tabul,n_ot_mu,n_ot_mu_G!mu(i)&n(mu(i))
      REAL*8 :: R_min,R_plu,sha_demu
! Duryunda
      LOGICAL :: DURYUNDA
! Spin polarized if true
      LOGICAL :: POLARIZED
      REAL*8 :: pola_count
! Flat statistics parameters
      REAL*8,ALLOCATABLE :: order(:), OOO(:) !Collected and normalized order stat
      REAL*8,ALLOCATABLE :: ord_probab_2(:), ord_probab_old_2(:), &
                                           ord_probab(:), ord_probab_0(:), ord_probab_old(:)
      REAL*8,ALLOCATABLE :: order_cy(:) ! Order statistics for one cycle between printouts
      REAL*8,ALLOCATABLE :: order_2(:), kappa(:)
      REAL*8,PARAMETER :: zuka=2.0d-5
! Momenta names and covariance matrix measurements for momenta and OC
      INTEGER,PARAMETER :: namo=30
      INTEGER,PARAMETER :: i_max_meas=20000000
      INTEGER :: i_cur_meas ! Number of current measurements 
      INTEGER :: i_all_meas ! Number of another group of current measuremnts
      CHARACTER*11,DIMENSION(namo) :: mome_names
!      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: covar !(mom_numb,it,i_meas)
!      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: covar_oc !(it_op,i_meas)
      INTEGER,DIMENSION(namo) :: insko
      REAL*8,DIMENSION(namo) :: resko, besko
      INTEGER :: momsko
! Phonon measuring propagator parameters
      REAL*8 :: meas_Debye, meas_V_kvadrat, meas_decay !MeasPhoPropPars
      REAL*8 :: theta_meas_2 !Spatial decay parameter
! Fit process
      REAL*8 :: bazeka, bazeka_min, bazeka_max  ! DecayInUnitsOfDebye for generations of phonons
      REAL*8 :: u_limit ! liit for short phonons 
! Frohlich checks
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: cou_cou
      END MODULE config_par
!.............................................................

!--------------------------------------------------------------------
! Declare statistical parameters
!--------------------------------------------------------------------
      MODULE stat_par
      IMPLICIT NONE
      SAVE
! Statistics
!!  Float Counters for Loop Termination (local)
      INTEGER :: print_count, write_count, iskipi_count  !CountersFor Print/Write (NULIF)
      INTEGER :: measu_count !CounterForMeasurements(NULIF)
      INTEGER :: i_oc_count  !CounterFor_OC_Measurements_skips_when_other_measurements_done(NULIF)
!!  Global Counters for Statistics
      REAL*8 :: flo_mont     !TotalMCSteps(NULIF)
      REAL*8 :: number_measu !NumberOfPerformedMeasurementsForGreen(NULIF)
      REAL*8 :: num_esti     !NumberOfPerformedMeasurementsForEstimat(NULIF)
      INTEGER :: in_print    !NimberOfPerformedPrintouts(NULIF)
!!  Input flags
      INTEGER :: print_limit      !Prints aftrer "_limit" updates
      INTEGER :: write_limit, iskipi_lim      !Writes aftrer "_limit" updates
      INTEGER :: measu_limit      !Measures aftrer "_limit" updates
      INTEGER :: measu_oc_limit   !Measures OC aftrer "measu_limit*measu_oc_limit" updates 
      LOGICAL :: tomeasure        !To measure or not
!!  Other
      INTEGER :: km                    !UpdateMarker
      INTEGER :: stat_is               !0/1 if notuse/use previous stat
      INTEGER :: use_conf              !0/1 if notuse/use previous conf
      REAL*8 :: sevide, opvide                !HowToDivideReadStatistics, optical 
!! Effectivity counters
      REAL*8 :: physical=0.0d0, unphysical=0.0d0;
      REAL*8 :: physical_c=0.0d0, unphysical_c=0.0d0;
      REAL*8,ALLOCATABLE,DIMENSION(:) :: in_order,in_order_2
      REAL*8 :: add_tau_shift=0.0d0, suc_tau_shift=0.0d0
      REAL*8 :: add_site_shift=0.0d0, suc_site_shift=0.0d0
      REAL*8 :: add_add_pho=0.0d0, suc_add_pho=0.0d0
      REAL*8 :: add_add_dia=0.0d0, suc_add_dia=0.0d0
      REAL*8 :: add_add_dia_exp=0.0d0, suc_add_dia_exp=0.0d0
      REAL*8 :: add_add_dia_sho=0.0d0, suc_add_dia_sho=0.0d0
      REAL*8 :: add_rem_pho=0.0d0, suc_rem_pho=0.0d0
      REAL*8 :: add_rem_dia=0.0d0, suc_rem_dia=0.0d0
      REAL*8 :: add_rem_dia_exp=0.0d0, suc_rem_dia_exp=0.0d0
      REAL*8 :: add_rem_dia_sho=0.0d0, suc_rem_dia_sho=0.0d0
      REAL*8 :: add_add_worm=0.0d0, suc_add_worm=0.0d0
      REAL*8 :: add_rem_worm=0.0d0, suc_rem_worm=0.0d0
      REAL*8 :: add_commute=0.0d0, suc_commute=0.0d0
      REAL*8 :: add_move_mea=0.0d0, suc_move_mea=0.0d0
      REAL*8 :: add_flip_spi=0.0d0, suc_flip_spi=0.0d0
      REAL*8 :: add_move_worm=0.0d0, suc_move_worm=0.0d0
      REAL*8 :: add_commute_int=0.0d0, suc_commute_int=0.0d0
! Processing of physical data
      INTEGER,DIMENSION(3) :: Fermi,Fermi2,FermiZ !Fermi vector
      REAL*8,DIMENSION(-10:10) :: Z_factor   !Z-factors at Fermi vector
      REAL*8 :: max_k_deriv                !Maximal dn(k_x)/dk_x 
! Immediate addressing counters
      REAL*8 :: a00_tau_shift=0.0d0
      REAL*8 :: a00_site_shift=0.0d0
      REAL*8 :: a00_add_pho=0.0d0
      REAL*8 :: a00_add_dia=0.0d0
      REAL*8 :: a00_add_dia_exp=0.0d0
      REAL*8 :: a00_add_dia_sho=0.0d0
      REAL*8 :: a00_rem_pho=0.0d0
      REAL*8 :: a00_rem_dia=0.0d0
      REAL*8 :: a00_rem_dia_exp=0.0d0
      REAL*8 :: a00_rem_dia_sho=0.0d0
      REAL*8 :: a00_add_worm=0.0d0
      REAL*8 :: a00_rem_worm=0.0d0
      REAL*8 :: a00_commute=0.0d0
      REAL*8 :: a00_move_mea=0.0d0
      REAL*8 :: a00_flip_spi=0.0d0
      REAL*8 :: a00_move_worm=0.0d0
      REAL*8 :: a00_commute_int=0.0d0
!! Norms
      REAL*8 :: Z_phys_norm, Z_self_norm, Z_pola_norm
      REAL*8 :: Z_TokTok_norm
      REAL*8 :: S_norm_now, P_norm_now
! Green function current and bare
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: GRT_NOW,GRT_0,GRT_NOW_SAVE, &
                                                                       PHO_0,PHO_NOW !TablesForProps
! Props for measuring
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: Sigma_box,Polar_box
! Current-current
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: TokTok_box, TokTok_box_all
! Current-current for control
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: ToTo_box_sep, ToTo_box_sep_no
      REAL*8,DIMENSION(4) :: Nor_TT, Nor_jj
! Props for measuring to put data to fourier analysis 
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: Sigma,Polar
! To put basis OC into grid 
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: ToTo_Basis
      
! Occupation numbers
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: nk_distr,nk_distr_0
! Convegency dimensions
      INTEGER,PARAMETER :: max_printi=20000000
      REAL*8,DIMENSION(max_printi) :: pilar_nac,sugma_nac,wowe
! Special statistical convergency dimensions
      REAL*8,DIMENSION(0:max_printi) :: dunsity_ST, kineti_ST
      REAL*8,DIMENSION(0:max_printi) :: Z_fac_ST, m_x_ST, m_xy_ST
! space seeding using log scales
      double precision, allocatable :: nmurelation(:) 
      double precision, parameter   :: mustep=1.0d-3  
      INTEGER,PARAMETER :: nmupoints=12000
      REAL*8 :: density
      REAL*8 :: uh_kinetik1
! For OC exact estimator
      INTEGER,PARAMETER :: skoko_cut = 5
      INTEGER :: num_oc, num_oc_lim
      INTEGER :: num_oc_ano, num_oc_lim_ano
      REAL*8 :: scalik_oc_ano
      REAL*8 :: window, case_max_actual
      LOGICAL :: MEASURE_OC, MEASURE_OC_ano
      REAL*8,ALLOCATABLE,DIMENSION(:) :: kakapu, ta_oc
      REAL*8,ALLOCATABLE,DIMENSION(:) :: kakapu_ano, ta_oc_ano
      REAL*8,ALLOCATABLE,DIMENSION(:) :: jj_exa, jj_exa_norm, jj_exa_proba, jj_exa_proba_ano
      REAL*8,ALLOCATABLE,DIMENSION(:) :: jj_exa_ano, jj_exa_norm_ano
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: jj_exa_cut, jj_exa_cut_norm, jj_exa_cut_ano, jj_exa_cut_norm_ano
      REAL*8,DIMENSION(0:skoko_cut) :: case_max_given, stat_cut
      CHARACTER(LEN=11),DIMENSION(0:skoko_cut) :: files_cut;
      CHARACTER(LEN=11),DIMENSION(0:skoko_cut) :: files_cut_ano;
      
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: jj_exa_sep,jj_exa_sep_no 
! For gradual enlarging diagram order      
      REAL*8 :: srednyaka,lupaka
! For cutting anomalous contributions to OC
      INTEGER,PARAMETER :: num_aga=20000
      REAL*8,PARAMETER :: aga_limit=1.0d4
      REAL*8,DIMENSION(-num_aga:num_aga) :: aga_value,aga_hist
! Estimating average number of loops
      REAL*8 :: loop_summed, time_summed
      END MODULE stat_par
!.............................................................

!--------------------------------------------------------------------
! Declare Z-factors parameters
!--------------------------------------------------------------------
      MODULE Zfac_par 
      IMPLICIT NONE
      SAVE
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: Sss,Ess,Zss
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: Fsurface
      
      REAL*8,DIMENSION(3) :: Fmom
      
      END MODULE Zfac_par
!.............................................................
     
!-------------------------------------------------------------
!     For random generator
!-------------------------------------------------------------
      MODULE rara
      IMPLICIT NONE
      SAVE
      INTEGER :: ir(521),iptr
      END MODULE rara
!.............................................................

      INCLUDE "FHO_serv_dummy_25.f90"
      INCLUDE "FHO_serv_phys_25.f90"
      INCLUDE "FHO_serv_bold_25.f90"
      INCLUDE "FHO_zfac_ma_25.f90"
      
!*******************************************************************
!        P R O G R A M S
!*******************************************************************

!*******************************************************************
!        MAIN
!!*******************************************************************
      PROGRAM FHO_08
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8 :: fact
      INTEGER :: iui
      
!      PRINT*,"READ iui"
!      READ*, iui
!      PRINT*,fact(iui)
!      PAUSE
      
      CALL ALLOCA_0
      CALL PREPARE
      CALL MONTE


      END PROGRAM FHO_08
!....................................................................

!--------------------------------------------------------------------
!     Factorial
!--------------------------------------------------------------------
      REAL*8 FUNCTION fact(l)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: l
      INTEGER :: i
      IF(l==0)THEN
       fact=1.0d0
      ELSE
        fact=1.0d0
       DO i=1,l
        fact = fact * i
       ENDDO
      ENDIF
      END FUNCTION fact
!....................................................................
      
      
!--------------------------------------------------------------------
! Allocating dimensions which do not depend on input
!--------------------------------------------------------------------
      SUBROUTINE alloca_0
      USE config_par;  USE stat_par
      IMPLICIT NONE
      INTEGER :: error

      ALLOCATE(nlist(N_max),numnam(N_max),STAT=error)
      ALLOCATE(type(N_max,3),link(N_max,3),STAT=error)
      ALLOCATE(site(N_max),tau(N_max),STAT=error)
      ALLOCATE(omega(N_max,3,3),STAT=error)
      ALLOCATE(nmurelation(0:10000000))

      ALLOCATE(linkh2(-2:N_max,2),linkh3(N_max,2),STAT=error)

      PRINT*,"Allocation error level in 0 = ",error

      END SUBROUTINE alloca_0
!....................................................................

!--------------------------------------------------------------------
! Allocating dimensions which depend on input
!--------------------------------------------------------------------
      SUBROUTINE alloca_1
      USE config_par;  USE stat_par
      IMPLICIT NONE
      INTEGER :: error

      ALLOCATE(ass(Nsite,dir),STAT=error)
      ALLOCATE(back(dir),STAT=error)

      ALLOCATE(GRT_0(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau), STAT=error)
      ALLOCATE(GRT_NOW(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau), STAT=error)
      ALLOCATE(GRT_NOW_SAVE(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau), STAT=error)
      ALLOCATE(PHO_0(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau), STAT=error)
      ALLOCATE(PHO_NOW(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau), STAT=error)
      ALLOCATE(Sigma_box(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      ALLOCATE(TokTok_box(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      ALLOCATE(ToTo_box_sep(0:Ntau-1,4),ToTo_box_sep_no(0:Ntau-1,4))
      ALLOCATE(TokTok_box_all(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      ALLOCATE(Polar_box(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1))
      ALLOCATE(nk_distr(0:L(1)-1,0:L(2)-1,0:L(3)-1), STAT=error)
      ALLOCATE(nk_distr_0(0:L(1)-1,0:L(2)-1,0:L(3)-1), STAT=error)

! Allocating for BOLD_G and BOLD_PH      
      ALLOCATE(Sigma(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1), STAT=error)
      ALLOCATE(Polar(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1), STAT=error)
      ALLOCATE(ToTo_Basis(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1),STAT=error)

      ALLOCATE(in_order(2:N_max_now+2), STAT=error)
      in_order(2:N_max_now)=nul
      ALLOCATE(in_order_2(1:N_max_now/2), STAT=error)

      allocate(Sigmab(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Nbins_s,Nbasis_s)) ! self-energy
      allocate(tce_s(0:Nbins_s-1))                          ! basis bins centers
      allocate(tbin_s(0:Nbins_s))                           ! basis bins boundaries 
      allocate(oc_s(0:Nbins_s,Nbasis_s,Nbasis_s))               ! orthonormal coeff. 

      allocate(Polarb(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Nbins_p,Nbasis_p)) ! polarization
      allocate(tce_p(0:Nbins_p-1))                          ! basis bins centers
      allocate(tbin_p(0:Nbins_p))                           ! basis bins boundaries 
      allocate(oc_p(0:Nbins_p,Nbasis_p,Nbasis_p))               ! orthonormal coeff. 
      
      allocate(OCb(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Nbins_OC,Nbasis_OC)) ! OC = optical conductivity
      allocate(tce_OC(0:Nbins_OC-1))                          ! basis bins centers
      allocate(tbin_OC(0:Nbins_OC))                           ! basis bins boundaries 
      allocate(oc_OC(0:Nbins_OC,Nbasis_OC,Nbasis_OC))               ! orthonormal coeff. 

! Allocating for flat statistics in orders
      ALLOCATE(order(1:N_max_max)) !Accumulated statistics in orders
      order(1:N_max_now)=nul
      ALLOCATE(order_cy(1:N_max_max)) !Accumulated statistics in orders for one cycle
      order_cy(1:N_max_now)=nul
      ALLOCATE(order_2(1:N_max_now)) !
      ALLOCATE(OOO(1:N_max_now))
      ALLOCATE(ord_probab(1:N_max_max))
      ALLOCATE(ord_probab_2(1:N_max_now))
      ALLOCATE(ord_probab_0(1:N_max_now))
      ALLOCATE(ord_probab_old(1:N_max_now))
      ALLOCATE(ord_probab_old_2(1:N_max_now))
      ALLOCATE(kappa(1:N_max_now))
      
! Allocating for Frohlicj check
      ALLOCATE(cou_cou(0:L(1)-1,0:L(2)-1))

      PRINT*,"Allocation error level in 1 = ",error

      END SUBROUTINE alloca_1
!....................................................................

!--------------------------------------------------------------------
! Checking sites associations to measure current-current
!--------------------------------------------------------------------
      SUBROUTINE ASS_VECTORS
      USE  config_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: RNDM
      INTEGER,DIMENSION(3) :: v1,v2,v12,v_dif,sh_i,sh_j
      INTEGER :: i1_site,i2_site,i3_site,i4_site,i_site
      INTEGER :: i,j,ii,iii
      INTEGER :: s1_sh, s2_sh, s1, s2, scalar
      
           
      DO
          
      s1=Nsite*RNDM(kk)
      s2=Nsite*RNDM(kk)

      DO i=1,dir
      DO j=1,dir
          s1_sh=ass(s1,i)
          s2_sh=ass(s2,j)
          CALL VEC_BETWEEN(s1,s1_sh,v1,v2,v12)
          sh_i(1:dim)=v2(1:dim)-v1(1:dim)
          DO iii=1,dim
              IF(sh_i(iii)>=2)THEN
                  sh_i(iii)=-1
              ELSE IF(sh_i(iii)<=-2)THEN
                  sh_i(iii)= 1    
              ENDIF    
          ENDDO    
          CALL VEC_BETWEEN(s2,s2_sh,v1,v2,v12)
          sh_j(1:dim)=v2(1:dim)-v1(1:dim)
          DO iii=1,dim
              IF(sh_j(iii)>=2)THEN
                  sh_j(iii)=-1
              ELSE IF(sh_j(iii)<=-2)THEN
                  sh_j(iii)= 1    
              ENDIF     
          ENDDO              
          scalar=0
          DO ii=1,dim; scalar = scalar + sh_i(ii)*sh_j(ii); ENDDO
          IF(scalar==0)CYCLE    
          PRINT*, sh_i(1:2), sh_j(1:2)
          PRINT*,scalar    
          PRINT*,'---------------------------------------------'
      ENDDO
      ENDDO
      
      PAUSE
          
      ENDDO    
      
      
      RETURN
      DO
    
      i_site=Nsite*RNDM(kk)
      i1_site=ass(i_site,1)
      i2_site=ass(i_site,2)
      i3_site=ass(i_site,3)
      i4_site=ass(i_site,4)
      
      CALL VEC_BETWEEN(i_site,i1_site,v1,v2,v12)
      v_dif(1:2)=v2(1:2)-v1(1:2)
      PRINT*,'1'
      PRINT*,v1(1:2),v2(1:2)
      PRINT*,v12(1:2),v_dif(1:2)
      
      CALL VEC_BETWEEN(i_site,i2_site,v1,v2,v12)
      v_dif(1:2)=v2(1:2)-v1(1:2)
      PRINT*,'2'
      PRINT*,v1(1:2),v2(1:2)
      PRINT*,v12(1:2),v_dif(1:2)
      
      CALL VEC_BETWEEN(i_site,i3_site,v1,v2,v12)
      v_dif(1:2)=v2(1:2)-v1(1:2)
      PRINT*,'3'
      PRINT*,v1(1:2),v2(1:2)
      PRINT*,v12(1:2),v_dif(1:2)

      CALL VEC_BETWEEN(i_site,i4_site,v1,v2,v12)
      v_dif(1:2)=v2(1:2)-v1(1:2)
      PRINT*,'4'
      PRINT*,v1(1:2),v2(1:2)
      PRINT*,v12(1:2),v_dif(1:2)
      PAUSE
      
      ENDDO
      
      END SUBROUTINE ASS_VECTORS
!....................................................................
      
      
!--------------------------------------------------------------------
! Preparing the program for start
!--------------------------------------------------------------------
      SUBROUTINE PREPARE
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: P_NORM, S_NORM, PHONON, PHONON_MEAS, GREEN
      REAL*8 :: s_norm_0, p_norm_0 ,a1 ,a2 ,b1 ,b2, ph_real, ph_meas
      INTEGER*4,DIMENSION(1:522) :: isave
      INTEGER :: ii1, jj1, ping, pung, pong, puff, i
      
      
! Initilizing numbers
        OPEN(4,FILE="oaoa.oa")
          WRITE(4,*)"fgt=08"
        CLOSE(4)
        PRINT*,"Done printing"
! Reading the updates addressing probabilities
      OPEN(4,FILE='probab.in')
        READ(4,*)dpro(1)                       !ShiftVorTime     km= 1
        READ(4,*)dpro(2)                       !ShiftVorSite     km= 2
        READ(4,*)dpro(3);    dpro(4) =dpro(3)  !Add/RemPhon      km= 4/5
        READ(4,*)dpro(5);    dpro(6) =dpro(5)  !Add/RemWorm      km= 6/7
        READ(4,*)dpro(7);                      !Commute          km= 9
        READ(4,*)dpro(8);                      !MoveMeasLine     km= 10
        READ(4,*)dpro(9);                      !FlipLoopSpins    km= 3
        READ(4,*)dpro(10);                     !MoveWorm         km= 8
        READ(4,*)dpro(11);   dpro(12)=dpro(11) !Add/RemDiagPhon  km= 11/12  
        READ(4,*)dpro(13);                     !CommuteInt       km= 13
        READ(4,*)dpro(14);   dpro(15)=dpro(14) !Add/RemDiagPhonExp  km= 14/15 
        READ(4,*)dpro(16);   dpro(17)=dpro(16) !Add/RemDiagShortPh  km= 16/17 
      CLOSE(4)
! Initilizing probabilities
      prob(1) =dpro(1)
      prob(2) =prob(1) +dpro(2)
      prob(3) =prob(2) +dpro(3);  prob(4) =prob(3) +dpro(4)
      prob(5) =prob(4) +dpro(5);  prob(6) =prob(5) +dpro(6)
      prob(7) =prob(6) +dpro(7);
      prob(8) =prob(7) +dpro(8)
      prob(9) =prob(8) +dpro(9);
      prob(10)=prob(9) +dpro(10);
      prob(11)=prob(10)+dpro(11); prob(12) =prob(11)+dpro(12) 
      prob(13)=prob(12)+dpro(13)
      prob(14)=prob(13)+dpro(14); prob(15) =prob(14)+dpro(15)
      prob(16)=prob(15)+dpro(16); prob(17) =prob(16)+dpro(17)
! Checking normalization
      IF(ABS(prob(17)-un1)>1.0d-10)THEN
        PRINT*,' WRONG PROCESSES NORMALISATION',prob(17);
        PRINT*,prob(1:17); STOP
      ENDIF
! Reading parameters
      CALL IN_DATA;
      PRINT*,"Ntau(0) =",Ntau

      IF(stat_is==0)THEN; DURYUNDA=.TRUE.
        ! Initializing the random-number generator
        CALL RANDINIT(iseed)
        ! Initialize lattice
        CALL INIT_LATT;
        ! Initilizing the n(mu)
        !Checking association vectors
        !CALL ASS_VECTORS
        CALL N_OT_MU_TABU        
        ! Initilizing worm weight
        worm_weight=worm_weight_0
        ! Initilizing bare propagators for consequent operations
        CALL GREEN0
        CALL PHONO_0
        ! Checking measuring function
        ph_real = PHONON(1,beta/2.0,1,beta/2.0)
        ph_meas = PHONON_MEAS(1,beta/2.0,1,beta/2.0)
        PRINT*,'REAL: ',ph_real,'   MEAS: ',ph_meas
        CALL BOLD0
        s_norm_0=S_NORM()
        p_norm_0=P_NORM()
        CALL DENSITY_CHECK
! Initilizing configuration and statistics
        CALL INIT_CONF; CALL INIT_STAT
! Drawing the diagram
        CALL DRAW
! Checking configuration
        CALL CHECK_CONFIG
        CALL normalize_s
        CALL normalize_p
        CALL normalize_OC
        i_cur_meas=0
        i_all_meas=0
      ELSE; DURYUNDA=.FALSE.  
        ! Nullify Local Counters
        print_count=0; write_count=0; measu_count=0; iskipi_count=0; i_oc_count=0;
        CALL INIT_LATT;
        ! Initilizing the n(mu)
        CALL N_OT_MU_TABU        
        ! Initilizing bare propagators
        CALL GREEN0
        CALL PHONO_0
        ! Initilizing configuration and statistics
        ping=0; pung=0; pong=0;
        IF(use_conf==1)THEN
            CALL RAND_READ(pong) 
            CALL CONF_READ(ping); 
        ELSE
            CALL RAND_READ(pong)
            CALL INIT_CONF;
            ping=0 
        ENDIF    
        CALL ST_READ(pung); 
        puff=MAX(ping,pung,pong);
        IF(puff/=0)THEN;
         PRINT*,' DATA FILES ARE CORRUPTED'
         PRINT*,' RECOVERING OF COPY in PROGRESS' 
         PRINT*,' CNF =',ping,' ST=',pung,'RND =',pong, ' TOT =',puff
        ENDIF
        
        OPEN(UNIT=12,FILE="corru_info.dat")
            WRITE(12,*)"CHECK CONTINUATION"
            WRITE(12,*)"ping / pung  / pong: ", ping, "  /  ", pung, "  /  ", pong
        CLOSE(12)  
        
        
        IF(use_conf==1)THEN
            CALL RAND_READ(puff)
            CALL CONF_READ(puff); 
            PRINT*,"!!!!!!!!Conf read!!!!!!!!!!!!!!"
        ELSE
            CALL RAND_READ(puff)
            CALL INIT_CONF;       
            PRINT*,"!!!!!!!!Conf init!!!!!!!!!!!!!!"
        ENDIF  
        CALL ST_READ(puff)
        CALL CHECK_CONFIG
        CALL normalize_s
        CALL normalize_p
        CALL normalize_OC
        CALL BOLDG; CALL BOLDPH
        s_norm_now=S_NORM()
        p_norm_now=P_NORM()
        CALL ENERGYTEST; Energy_a=ABS(energy)
        CALL PHASETEST ; Dphase_a=Dphase
        PRINT*,"READINF TESTS PASSED POSITIVE"
      ENDIF

      !Block analysis of the C contributions
      DO i=-num_aga,num_aga;
         aga_value(i) = (aga_limit/num_aga) * i;
      ENDDO;
         
      PRINT*,"Current measurement for covariance is :",i_cur_meas
      PRINT*,"Current measurement for      all        is :",i_all_meas
      
      PRINT*,"Worm_0 = ",worm_weight_0,"  Worm =",worm_weight
      
      PRINT*,"================ V E R S I O N ==================  "
      PRINT*," FHO_24: 2019.02.27 No 1: Nhlist ene_diff CovarMat P-meas"
      PRINT*,"OC, skip_meas, OC-4ext, 1D, 2D, Exp-phonAdRem PhoSHort"
      PRINT*,"Reduce memory  T_kinetic per particle Basis Funct OC "
      PRINT*,"Useful statistics reducable, exa_oc - FORTRAN90!!! check_02"
      PRINT*,"Change order, OC covariance, RARE STA writung, NO COVARIANCE WRITE"
      PRINT*,"Enhanced long run counters, order adding enhancement"
      PRINT*,"ZUKA =   ZUKA = ",zuka, "LUPAKA divide, change moments!!"
      PRINT*,"Fighting with curvy <jj> new values opti divide  "
      PRINT*,"Amotehr OC  Another oC_CUT, loop meas"
      PRINT*,"============ Physical model parameters ==========  "
      PRINT"('hopping = ',ES11.4,' Debye = ',ES11.4)",hopping, Debye
      PRINT"('     mu = ',ES11.4,'     T = ',ES11.4)",mu,un1/beta
      PRINT"(' lambda = ',ES11.4,'   Dim = ',I3)",lambda,dim
      PRINT"('nmnm = ',I4)",nmnm
      PRINT*,"ene_diff = ",ene_diff
      PRINT*,"=================================================  "
      

      !PRINT*,"OGO = ",un2 * ABS(GREEN(1,nul,1,nul,1)) 
       
      
      END SUBROUTINE prepare
!....................................................................


!--------------------------------------------------------------------
! Initilazing the initial coniguration
!--------------------------------------------------------------------
      SUBROUTINE IN_DATA
      USE config_par; USE stat_par;
      IMPLICIT NONE
      INTEGER :: i, isapa
      REAL*8 :: anorma0, zuza, dubara, dubara_ano

      OPEN(UNIT=4,FILE="in.in")
      READ(4,*)stat_is           !0/1 if notuse/use statistics
      READ(4,*)use_conf          !0/1 if notuse/use previous config
      READ(4,*)sevide,opvide     !ToDevideStatistics 
      READ(4,*)hopping,Debye     !hopping,Debye
      READ(4,*)lambda, i_type, correca    !lambda/alpha, interactio type 0-Holstein, 1-Frohlich 
      READ(4,*)worm_weight_0, w_dia !WeightOfWorm
      READ(4,*)beta,mu           !InverseT and mu
      READ(4,*)Ntau              !NumberOfTauPoints
      READ(4,*)L(1:3)            !SizesOfBoxOfSystem
      READ(4,*)N_max_now,N_mes_now !Max NumbVort,NumbVortToMeasure
      READ(4,*)dim,dir           !Dimension&NearNeighborsNumber
      READ(4,*)print_limit       !Prints aftrer "_limit" updates
      READ(4,*)write_limit, iskipi_lim  !!!!!!Writes aftrer " _limit * iskipi_lim " updates
      READ(4,*)measu_limit       !Measures aftrer "_limit" updates
      READ(4,*)tomeasure         !To measure or not
      READ(4,*)iseed             !RandomNumberSeed
      READ(4,*)severe_check      !SevereCheck_IF_TRUE: !Slows down
      READ(4,*)mu00              !ParameterForSelfConsistencyLoop
      READ(4,*)renor_pho         !Whether to renormaliza phonons
      READ(4,*)Nbins_s,Nbasis_s  !Basis for Sigma
      READ(4,*)Nbins_p,Nbasis_p  !Basis for Polar
      READ(4,*)Nbins_OC,Nbasis_OC  !Basis for OC
      READ(4,*)POLARIZED         !Spin-polarized if true
      READ(4,*)meas_Debye,meas_V_kvadrat,meas_decay !Measuring line pars
      READ(4,*)bazeka_min, bazeka_max    !Decay in Debye for phonon generations
      READ(4,*)u_limit    !Limit for short phonons
      READ(4,*)num_oc,window,MEASURE_OC,measu_oc_limit !Parameters for OC measuring
      READ(4,*)num_oc_ano,scalik_oc_ano,MEASURE_OC_ano !Parameters for ANOTHER OC measuring
      READ(4,*)srednyaka, lupaka  !to higher orders: statistics, divide last order probab
      CLOSE(4)
      
      IF(iskipi_lim==0)STOP" Must be iskipi_lim >= 1"
      
      IF(N_max_now>N_max_max)STOP"Icrease N_max_max"
      IF(N_mes_now>N_max_max)STOP"Icrease N_max_max"
      
      IF(window >= beta/2.0d0)STOP"Too large window"
      
      IF(u_limit>beta*0.5d0)THEN
          u_limit=beta*0.5
          PRINT*,"U_LIMIT set to : ",u_limit
      ENDIF
      
      IF(severe_check)THEN
          PRINT*,'**************************************'
          PRINT*,'**************************************'
          PRINT*,'   SEVERE CHECK IS ENABLED'
          PRINT*,'    IT SLOWS DOWN CODE!!!!'
          PRINT*,'    CONFIRM TO CONTINUE   '
          PRINT*,'**************************************'
          PRINT*,'**************************************'
          PAUSE
      ENDIF    
      
      IF(sup_severe_check)THEN
          PRINT*,'**************************************'
          PRINT*,'**************************************'
          PRINT*,'   SUP_SEVERE CHECK IS ENABLED'
          PRINT*,'    IT SLOWS DOWN CODE!!!!'
          PRINT*,'    CONFIRM TO CONTINUE   '
          PRINT*,'**************************************'
          PRINT*,'**************************************'
          PAUSE
      ENDIF    

      bazeka_min=Debye*bazeka_min ; PRINT*,"bazeka_min  = ",bazeka_min
      bazeka_max=Debye*bazeka_max ; PRINT*,"bazeka_max  = ",bazeka_max
      
      IF(POLARIZED)THEN
          pola_count=1.0d0
      ELSE
          pola_count=2.0d0
      ENDIF    
      
! Derivable initial data
  !Time
      dtau=beta/Ntau
  !Lattice
      Nsite=L(1)*L(2)*L(3)
      LP(1)=1; LP(2)=L(1); LP(3)=L(1)*L(2); LP(4)=L(1)*L(2)*L(3)
  
!Model
      IF(i_type==0)THEN
          V_kvadrat = un2 * dim * hopping * Debye * lambda
      ELSE IF(i_type==1)THEN
          V_kvadrat = SQRT(un2) * lambda
      ELSE
          STOP"There is no such model"
      ENDIF    
          
      CALL ALLOCA_1
  
 ! Initializing probability for flat stistics in order
      ord_probab_0=nul; ord_probab=nul; ord_probab_old=nul; 
      OPEN(UNIT=4,FILE="flat.in")
      anorma0=0.0d0
      PRINT*,"anorma0 =", anorma0
      DO i=2,N_max_now,2
         READ(UNIT=4,FMT=*,ERR=22,END=24)zuza
         PRINT*,i,zuza
         kappa(i)=zuza
         ord_probab_0(i)=kappa(i)/N_max_now
         anorma0 = anorma0 + ord_probab_0(i)
      ENDDO   
      CLOSE(4)
      ord_probab_0(1:N_max_now) = ord_probab_0(1:N_max_now) / anorma0
      ord_probab(1:N_max_now) = ord_probab_0(1:N_max_now) !Suggesting if no statistics
      PRINT*,"Init order probab norma : ", SUM(ord_probab(2:N_max_now))
      
      GOTO 23
22    STOP"No flat statistics given: ERR"      
24    STOP"No flat statistics given: END"      
23    CONTINUE      

      worm_max = worm_weight_0 * w_dia
      worm_min = worm_weight_0 / w_dia 
      
  ! All for basis function statistics: for Sigma   
      dtb_s=beta/Nbins_s
      DO i=0,Nbins_s; tbin_s(i)=dtb_s*i; ENDDO
      DO i=0,Nbins_s-1; tce_s(i)=dtb_s*(i+0.5d0); ENDDO    
  ! All for basis function statistics: for Polar   
      dtb_p=beta/Nbins_p
      DO i=0,Nbins_p; tbin_p(i)=dtb_p*i; ENDDO
      DO i=0,Nbins_p-1; tce_p(i)=dtb_p*(i+0.5d0); ENDDO    
  ! All for basis function statistics: for Polar   
      dtb_OC=beta/Nbins_OC
      DO i=0,Nbins_OC; tbin_OC(i)=dtb_OC*i; ENDDO
      DO i=0,Nbins_OC-1; tce_OC(i)=dtb_OC*(i+0.5d0); ENDDO    


! Preparing moments where G(k,tau) will be given          
      mome_names(1)="mome_01.dat"
      mome_names(2)="mome_02.dat"
      mome_names(3)="mome_03.dat"
      mome_names(4)="mome_04.dat"
      mome_names(5)="mome_05.dat"
      mome_names(6)="mome_06.dat"
      mome_names(7)="mome_07.dat"
      mome_names(8)="mome_08.dat"
      mome_names(9)="mome_09.dat"
      mome_names(10)="mome_10.dat"
      mome_names(11)="mome_11.dat"
      mome_names(12)="mome_12.dat"
      mome_names(13)="mome_13.dat"
      mome_names(14)="mome_14.dat"
      mome_names(15)="mome_15.dat"
      mome_names(16)="mome_16.dat"
      mome_names(17)="mome_17.dat"
      mome_names(18)="mome_18.dat"
      mome_names(19)="mome_19.dat"
      mome_names(20)="mome_20.dat"
      mome_names(21)="mome_21.dat"
      mome_names(22)="mome_22.dat"
      mome_names(23)="mome_23.dat"
      mome_names(24)="mome_24.dat"
      mome_names(25)="mome_25.dat"
      mome_names(26)="mome_26.dat"
      mome_names(27)="mome_27.dat"
      mome_names(28)="mome_28.dat"
      mome_names(29)="mome_29.dat"
      mome_names(30)="mome_30.dat"
      
!Initilizing which momenta to measure       
      OPEN(UNIT=4,FILE="gre_gde.in")
          READ(4,*)momsko
          IF(momsko>namo)STOP"Increase namo!"
          DO isapa=1,momsko
              READ(4,*)INSKO(isapa)
              RESKO(isapa) = (pi_m2/L(1))*INSKO(isapa)
              BESKO(isapa) = pi_m2 - RESKO(isapa) 
              PRINT*,"MOMENTA : ",INSKO(isapa),RESKO(isapa), mome_names(isapa)
          ENDDO    
      CLOSE(4)
      
! Initilizing covariance matrix measurements for Green finction
!      ALLOCATE(covar(momsko,0:Ntau-1,i_max_meas))
      
! Preparing imaginary time mesh for OC exact points estimator *****  
      num_oc_lim = 2 * num_oc
      ALLOCATE(kakapu(0:num_oc),ta_oc(0:num_oc_lim))
      ALLOCATE(jj_exa(0:num_oc_lim),jj_exa_norm(0:num_oc_lim))
      dubara = beta /  ( 2.0d0*( num_oc+(num_oc**2) ) ) 
      DO i = 0, num_oc
          kakapu(i) = dubara * ( i + (i**2) )
      ENDDO
      DO i = 0, num_oc_lim      
          IF( i <= num_oc )THEN
              ta_oc(i) = kakapu(i)
          ELSE IF( i > num_oc )THEN
              ta_oc(i) = beta - kakapu(num_oc_lim-i)
          ENDIF
      ENDDO
! End preparing **************************************************      
! Preparing imaginary time mesh for ANOTHER OC exact points estimator  
      num_oc_lim_ano = 2 * num_oc_ano
      ALLOCATE(kakapu_ano(0:num_oc_ano),ta_oc_ano(0:num_oc_lim_ano))
      ALLOCATE(jj_exa_ano(0:num_oc_lim_ano),jj_exa_norm_ano(0:num_oc_lim_ano))
      dubara_ano = (beta * scalik_oc_ano) /  ( 2.0d0*( num_oc_ano+(num_oc_ano**2) ) ) 
      PRINT*,"dubara_ano = ",dubara_ano, "scalik_oc_ano = ",scalik_oc_ano
      DO i = 0, num_oc_ano
          kakapu_ano(i) = dubara_ano * ( i + (i**2) )
      ENDDO
      DO i = 0, num_oc_lim_ano      
          IF( i <= num_oc_ano )THEN
              ta_oc_ano(i) = kakapu_ano(i)
          ELSE IF( i > num_oc_ano )THEN
              ta_oc_ano(i) = beta - kakapu_ano(num_oc_lim_ano-i)
          ENDIF
          !PRINT*,i,ta_oc_ano(i)
      ENDDO
      !PAUSE
! End preparing ****************************************************      
! Preparing "jj_exa_cut" estimator
      ALLOCATE(jj_exa_proba(0:num_oc_lim),jj_exa_proba_ano(0:num_oc_lim_ano))
      ALLOCATE(jj_exa_cut(0:skoko_cut,0:num_oc_lim),jj_exa_cut_norm(0:skoko_cut,0:num_oc_lim))
      ALLOCATE(jj_exa_cut_ano(0:skoko_cut,0:num_oc_lim_ano),jj_exa_cut_norm_ano(0:skoko_cut,0:num_oc_lim_ano))
      case_max_given(0) = HUGE(1.0d0);  files_cut(0) = "oc_cut0.dat"; files_cut_ano(0) = "oc_ano0.dat";
      case_max_given(1) = 5.0d5;        files_cut(1) = "oc_cut1.dat"; files_cut_ano(1) = "oc_ano1.dat";
      case_max_given(2) = 1.0d5;        files_cut(2) = "oc_cut2.dat"; files_cut_ano(2) = "oc_ano2.dat";    
      case_max_given(3) = 5.0d4;        files_cut(3) = "oc_cut3.dat"; files_cut_ano(3) = "oc_ano3.dat";      
      case_max_given(4) = 1.0d4;        files_cut(4) = "oc_cut4.dat"; files_cut_ano(4) = "oc_ano4.dat";      
      case_max_given(5) = 5.0d3;        files_cut(5) = "oc_cut5.dat"; files_cut_ano(5) = "oc_ano5.dat";      
      stat_cut(0:skoko_cut) = 0.0d0; 
! Initilizing covariance matrix mesurements for OC
!      ALLOCATE(covar_oc(0:num_oc_lim,i_max_meas))
                    
! Makes zeros for dimensions intended for checks
      ALLOCATE(jj_exa_sep(0:num_oc_lim,4),jj_exa_sep_no(0:num_oc_lim,4))
      ToTo_box_sep(0:Ntau-1,1:4)=nul
      jj_exa_sep(0:num_oc_lim,1:4)=nul
      
! Making write_limit = print_limit for correct covariance measirement 
      IF(print_limit/=write_limit)THEN
         write_limit = print_limit
         PRINT*,"*********************************************"
         PRINT*,"*** we made  write_limit = print_limit ******"
         PRINT*,"*********************************************"
      ENDIF    
      
      PRINT*,"Initial data ready"

      END SUBROUTINE in_data
!....................................................................

!--------------------------------------------------------------------
! Initilazing the initial statistics
!--------------------------------------------------------------------
      SUBROUTINE INIT_STAT
      USE stat_par;  USE config_par
      IMPLICIT NONE

! Nullify Global Counters
      flo_mont=nul; number_measu=nul; in_print=0; num_esti=nul
! Nullify Local Counters
      print_count=0; write_count=0; measu_count=0; iskipi_count=0;
! Nullify update marker
      km=0
! Nullify normalization counters
      Z_phys_norm=nul; Z_self_norm=nul; Z_pola_norm=nul;
      Z_TokTok_norm=nul
!
!      Sigma=nul;  Polar=nul ! Now it is only for boldifications 

      Sigma_box=nul;  Polar_box=nul ! Old boxes bin statistics
      Sigmab=nul; Polarb=nul; OCb=nul;      ! Basis function sigmas
      TokTok_box = nul; jj_exa = nul; jj_exa_cut = nul;  ! <jj> for OC 
      jj_exa_ano = nul; jj_exa_cut_ano = nul; 
       

! This is for OC large contributions control      
      aga_hist(-num_aga:num_aga)=nul;
      
! Number of loops summ
      loop_summed = nul; time_summed = nul;
      
      END SUBROUTINE INIT_STAT
!....................................................................

!--------------------------------------------------------------------
! Initilazing the initial coniguration
!--------------------------------------------------------------------
      SUBROUTINE INIT_CONF
      USE config_par;
      IMPLICIT NONE
      INTEGER :: j
      INTEGER :: iv_01,iv_02,iv_03,iv_04

! Preparing data for name manager
      nmnm=0
      DO j=1,N_max;
        nlist(j)=j; numnam(j)=j;
      ENDDO;
! Setting worm no-worm attributes
      present=.FALSE.
      sveta=0; tonya=0;
! Initilizing configuration
      
 ! Here we start from 4th order local in space Fock diagram
 !Getting names
 !     CALL GETNAME(iv_01); CALL GETNAME(iv_02)
 !     CALL GETNAME(iv_03); CALL GETNAME(iv_04)
 !Setting times
 !     tau(iv_01)=beta/10; tau(iv_02)=9*beta/10;
 !     tau(iv_03)=6*beta/10; tau(iv_04)=3*beta/10;
 !Setting sites
 !     IF(Nsite<4)STOP'system size is too small for initial config'
 !     site(iv_01)=2; site(iv_02)=2; site(iv_03)=2; site(iv_04)=2;
  !Settinh types
 !     type(iv_01,1)=1 ; type(iv_01,2)=1; type(iv_01,3)=1
 !     type(iv_02,1)=1 ; type(iv_02,2)=1; type(iv_02,3)=1
 !     type(iv_03,1)=1 ; type(iv_03,2)=1; type(iv_03,3)=1
 !     type(iv_04,1)=1 ; type(iv_04,2)=1; type(iv_04,3)=1
 !Setting links
 !     link(iv_01,1)=iv_04; link(iv_01,2)=iv_03; link(iv_01,3)=iv_02
 !     link(iv_02,1)=iv_03; link(iv_02,2)=iv_04; link(iv_02,3)=iv_01
 !     link(iv_03,1)=iv_01; link(iv_03,2)=iv_02; link(iv_03,3)=iv_04
 !     link(iv_04,1)=iv_02; link(iv_04,2)=iv_01; link(iv_04,3)=iv_03
 !Setting measuring line
 !     vert1=iv_01; vert2=iv_03; MEASV=.FALSE.
!      vert1=iv_01; vert2=iv_02; MEASV=.TRUE.
 !Setting momenta
 !     omega(iv_01,1,1:3)=1.0; omega(iv_01,2,1:3)=2.0;
 !     omega(iv_01,3,1:3)=-1.0
 !     omega(iv_02,1,1:3)=1.5; omega(iv_02,2,1:3)=0.5;
 !     omega(iv_02,3,1:3)=1.0
 !     omega(iv_03,1,1:3)=2.0; omega(iv_03,2,1:3)=1.5;
 !     omega(iv_03,3,1:3)=0.5
 !     omega(iv_04,1,1:3)=0.5; omega(iv_04,2,1:3)=1.0;
 !     omega(iv_04,3,1:3)=-0.5
! Initilizing hash table
 !     ha2=0; hl2=0; ha3=0; hl3=0; linkh2=0; linkh3=0
 !     CALL ADDHASH(omega(iv_01,2,1),2,iv_01)
 !     CALL ADDHASH(omega(iv_02,2,1),2,iv_02)
 !     CALL ADDHASH(omega(iv_03,2,1),2,iv_03)
 !     CALL ADDHASH(omega(iv_04,2,1),2,iv_04)
 !     CALL ADDHASH(DABS(omega(iv_01,3,1)),3,iv_01)
 !     CALL ADDHASH(DABS(omega(iv_02,3,1)),3,iv_02)
 !     CALL ADDHASH(DABS(omega(iv_03,3,1)),3,iv_03)
 !     CALL ADDHASH(DABS(omega(iv_04,3,1)),3,iv_04)

! Here we start from 2nd order local in space Fock diagram      
! Getting names
      CALL GETNAME(iv_01); CALL GETNAME(iv_02)
! Setting times      
      tau(iv_01)=5*beta/10; tau(iv_02)=6*beta/10;
! Setting sites      
      IF(Nsite<4)STOP'system size is too small for initial config'
      site(iv_01)=2; site(iv_02)=2; 
! Seting types
      type(iv_01,1)=1 ; type(iv_01,2)=1; type(iv_01,3)=1
      type(iv_02,1)=1 ; type(iv_02,2)=1; type(iv_02,3)=1
! Setting links
      link(iv_01,1)=iv_02; link(iv_01,2)=iv_02; link(iv_01,3)=iv_02
      link(iv_02,1)=iv_01; link(iv_02,2)=iv_01; link(iv_02,3)=iv_01
 !Setting measuring line
!      vert1=iv_01; vert2=iv_02; MEASV=.FALSE.
      vert1=iv_01; vert2=iv_02; MEASV=.TRUE.
! Setting momenta      
      omega(iv_01,1,1:3)=2.0; omega(iv_01,2,1:3)=3.0;
      omega(iv_01,3,1:3)=-1.0
      omega(iv_02,1,1:3)=3.0; omega(iv_02,2,1:3)=2.0;
      omega(iv_02,3,1:3)=1.0
! Initilizing hash table
      ha2=0; hl2=0; ha3=0; hl3=0; linkh2=0; linkh3=0
      CALL ADDHASH(omega(iv_01,2,1),2,iv_01)
      CALL ADDHASH(omega(iv_02,2,1),2,iv_02)
      CALL ADDHASH(DABS(omega(iv_01,3,1)),3,iv_01)
      CALL ADDHASH(DABS(omega(iv_02,3,1)),3,iv_02)
      
      
! Checking weight and phase
      CALL ENERGYTEST; Energy_a=ABS(energy)
      CALL PHASETEST ; Dphase_a=Dphase


      PRINT*,"Energy init = ",energy
      PRINT*,"Sign   init = ",Dphase

      PRINT*,"Initial configuration initiated"
      END SUBROUTINE init_conf
!....................................................................

!--------------------------------------------------------------------
! Initilazing the lattice
!--------------------------------------------------------------------
      SUBROUTINE init_latt
      USE config_par;
      IMPLICIT NONE
      INTEGER,DIMENSION(N_maxlat) :: ivec
      INTEGER :: i,site0,site1

! Linear chain in 1D ----
      IF(dir==2) THEN

      DO i=1,dim;
          back(i)=i+dim; back(i+dim)=i;
      ENDDO ! inv. directions

      DO  site0=1,Nsite                            ! lattice site index
              site1=site0-1; ivec=0;
      DO i=dim,1,-1
              ivec(i)=site1/LP(i)
              site1=site1-ivec(i)*Lp(i);
      ENDDO ! converting to vector

      DO i=1,dim;  site1=site0+LP(i)
      IF(ivec(i)==L(i)-1) site1=site1-LP(i+1)
      ass(site0,i)=site1; ass(site1,back(i))=site0
      ENDDO                                    ! linking near neighbors
      ENDDO                                    ! linear chain done
! ------

! Square lattice in 2D
      ELSE IF(dir==4) THEN

      DO i=1,dim;
          back(i)=i+dim; back(i+dim)=i;
      ENDDO ! inv. directions
      DO  site0=1,Nsite                            ! lattice site index
              site1=site0-1; ivec=0;
      DO i=dim,1,-1
              ivec(i)=site1/LP(i)
              site1=site1-ivec(i)*Lp(i);
      ENDDO ! converting to vector

      DO i=1,dim;  site1=site0+LP(i)
      IF(ivec(i)==L(i)-1) site1=site1-LP(i+1)
      ass(site0,i)=site1; ass(site1,back(i))=site0
      ENDDO                                    ! linking near neighbors

      ENDDO                                    ! done
! ----

! Triangular lattice in 2D
      ELSE IF(dir==6) THEN

!               \   |
!                 5 2
!                  \|
!              --3--o--1--  meaning of direction index
!                   |\
!                   4 6
!                   |   \
      DO i=1,dim;
         back(i)=i+dim; back(i+dim)=i;
      ENDDO
      back(5)=6; back(6)=5
      DO  site0=1,Nsite                        ! do square lattice first
              site1=site0-1; ivec=0;
      DO i=dim,1,-1
              ivec(i)=site1/LP(i)
              site1=site1-ivec(i)*Lp(i);
      ENDDO ! converting to vector

      DO i=1,dim;  site1=site0+LP(i);
      IF(ivec(i)==L(i)-1) site1=site1-LP(i+1)
      ass(site0,i)=site1; ass(site1,back(i))=site0
      ENDDO                                    ! linking near neighbors
      ENDDO                                    ! square lattice done
      DO  site0=1,Nsite                        ! finish triang. lattice
      site1=ass(site0,2); site1=ass(site1,3);
      ass(site0,5)=site1; ass(site1,back(5))=site0
      ENDDO                                    ! triang. lattice done
      ENDIF

      PRINT*,"Lattice initiated"
      END SUBROUTINE init_latt
!....................................................................

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      SUBROUTINE pr_stop(toprint)
      USE config_par;
      IMPLICIT NONE
      CHARACTER(LEN=47) :: toprint

      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      PRINT*,"!",toprint,"!"
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      CALL DRAW

      END SUBROUTINE pr_stop
!....................................................................

!--------------------------------------------------------------------
! Green functions NOW to use in calculations interpolated from table
!--------------------------------------------------------------------
      REAL*8 FUNCTION GREEN(site1,tau1,site2,tau2,spin)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: tau1,tau2
      INTEGER,INTENT(IN) :: site1,site2,spin
      INTEGER,DIMENSION(3) :: v_1,v_2,v_12
      REAL*8 :: tau_12, prl, prr, pr0, znak
      REAL*8 :: malo=1.0D-14, malo1=1.0D-14
      INTEGER :: i, it1

      CALL VEC_BETWEEN(site1,site2,v_1,v_2,v_12)

      tau_12=tau2-tau1; znak=un1;
      IF(DABS(tau_12)<malo) tau_12=-malo
      IF(tau_12<nul) then; tau_12=beta+tau_12; znak=-1.d0; endif
      it1 = INT(tau_12/dtau)
      IF(it1==Ntau)THEN
          PRINT*,"tau1/tau2 : ",tau1,'/',tau2
          PRINT*,"tau_12 =",tau_12
          STOP'too much it1'
      ENDIF
          
      prl=GRT_NOW(v_12(1),v_12(2),v_12(3),it1)
      prr=GRT_NOW(v_12(1),v_12(2),v_12(3),it1+1)
      green = prl + ((prr-prl)/dtau)*(tau_12-(it1*dtau))
      green=green*znak

      IF(green>=nul)THEN;
         Gphase=un1;
      ELSE
         Gphase=-un1;
      ENDIF

      END FUNCTION green
!...................................................................

!--------------------------------------------------------------------
! Green functions NOW to use in calculations interpolated from table
!--------------------------------------------------------------------
      REAL*8 FUNCTION GREEN_000(site1,tau1,site2,tau2,spin)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: tau1,tau2
      INTEGER,INTENT(IN) :: site1,site2,spin
      INTEGER,DIMENSION(3) :: v_1,v_2,v_12
      REAL*8 :: tau_12, prl, prr, pr0, znak
      REAL*8 :: malo=1.0D-14, malo1=1.0D-14
      INTEGER :: i, it1

      CALL VEC_BETWEEN(site1,site2,v_1,v_2,v_12)

      tau_12=tau2-tau1; znak=un1;
      IF(DABS(tau_12)<malo) tau_12=-malo
      IF(tau_12<nul) then; tau_12=beta+tau_12; znak=-1.d0; endif
      it1 = INT(tau_12/dtau)
      IF(it1==Ntau)STOP'too much it1'

      prl=GRT_0(v_12(1),v_12(2),v_12(3),it1)
      prr=GRT_0(v_12(1),v_12(2),v_12(3),it1+1)
      green_000 = prl + ((prr-prl)/dtau)*(tau_12-(it1*dtau))
      green_000=green_000*znak

      IF(green_000>=nul)THEN;
         Gphase=un1;
      ELSE
         Gphase=-un1;
      ENDIF

      END FUNCTION green_000
!....................................................................
      
      
!--------------------------------------------------------------------
! Phono propagator NOW to use in calculations interpolated from table
!--------------------------------------------------------------------
      REAL*8 FUNCTION PHONON(site1,tau1,site2,tau2)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: tau1,tau2
      INTEGER,INTENT(IN) :: site1,site2
      INTEGER,DIMENSION(3) :: v_1,v_2,v_12
      REAL*8 :: tau_12, prl, prr, pr0
      INTEGER :: i, it1

      tau_12=tau2-tau1; IF(tau_12<nul)tau_12=beta+tau_12
      it1 = INT(tau_12/dtau)
      IF(it1==Ntau)STOP'too much it1'

      CALL VEC_BETWEEN(site1,site2,v_1,v_2,v_12)

      prl=PHO_NOW(v_12(1),v_12(2),v_12(3),it1)
      prr=PHO_NOW(v_12(1),v_12(2),v_12(3),it1+1)
      phonon = prl + ((prr-prl)/dtau)*(tau_12-(it1*dtau))

      IF(phonon>=nul)THEN;
         Vphase=un1;
      ELSE
         Vphase=-un1;
      ENDIF
   
      
      END FUNCTION PHONON
!....................................................................
      
!--------------------------------------------------------------------
! Coulomb potential 1/r in 2D
!--------------------------------------------------------------------
    REAL*8 FUNCTION COULOMB(ix,iy,iz)
    USE config_par
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ix,iy,iz
    REAL*8 :: r1, r2, r3

105 FORMAT(1x,5(F10.4," "))    
    
    !COULOMB=1.0d0; RETURN;
     
    COULOMB=0.0d0
    IF(ABS(ix)+ABS(iy)+ABS(iz) == 0)THEN
        COULOMB=correca; r1=ix; r2=iy; r3=iz
    ELSE   
        IF(L(1)==1)THEN; 
            r1 = 0.0d0
        ELSE 
            r1 = MIN(ix,L(1)-ix);
        ENDIF    
        IF(L(2)==1)THEN
            r2 = 0.0d0
        ELSE  
            r2 = MIN(iy,L(2)-iy); 
        ENDIF    
        IF(L(3)==1)THEN;
            r3 = 0.0d0
        ELSE 
            r3 = MIN(iz,L(3)-iz);
        ENDIF    
        COULOMB = 1.0d0 / SQRT ( (r1**2) + (r2**2) + (r3**2) )
    ENDIF    
    !PRINT 105,r1,r2,r3, COULOMB
    cou_cou(ix,iy)=COULOMB
    
    
    
    END FUNCTION COULOMB
!....................................................................

!--------------------------------------------------------------------
! Phono propagator NOW to use in calculations interpolated from table
!--------------------------------------------------------------------
      REAL*8 FUNCTION PHONON_MEAS(site1,tau1,site2,tau2)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: tau1,tau2
      INTEGER,INTENT(IN) :: site1,site2
      REAL*8,EXTERNAL :: PHOZERO_MEAS
      INTEGER,DIMENSION(3) :: v_1,v_2,v_12
      REAL*8 :: tau_12, pr0
      INTEGER :: i
      
      tau_12=tau2-tau1; IF(tau_12<nul)tau_12=beta+tau_12
      CALL VEC_BETWEEN(site1,site2,v_1,v_2,v_12)

      phonon_meas = PHOZERO_MEAS(tau_12); 
      
      pr0 = 0.0d0
      DO i=1,dim; 
         v_12(i)=MIN(V_12(i),L(i)-V_12(i))
         pr0 = pr0 + (v_12(i))**2; 
      ENDDO 
      
      phonon_meas = phonon_meas * SQRT(theta_meas_2) / SQRT( pr0 +theta_meas_2 )
       
      Vphase=un1;
      
      END FUNCTION PHONON_MEAS
!....................................................................
      
      
!--------------------------------------------------------------------
! Update diagram
!--------------------------------------------------------------------
      SUBROUTINE UPDATE
      USE config_par; USE stat_par; IMPLICIT NONE
      REAL*8,EXTERNAL :: RNDM
      REAL*8 :: choice

      choice=RNDM(kk)

      IF(choice<prob(1)) THEN
        CALL TAU_SHIFT;                                            ! P-measure verified km=1
      ELSE IF(choice>prob(1)  .AND. choice<prob(2) )THEN
        CALL SITE_SHIFT;                                           ! P-measure verified km=2
      ELSE IF(choice>prob(2)  .AND. choice<prob(3) )THEN;
        CALL ADD_PHONON;                                         ! P-measure verified km=4
      ELSE IF(choice>prob(3)  .AND. choice<prob(4) )THEN;
        CALL REM_PHONON;                                         ! P-measure verified km=5
      ELSE IF(choice>prob(4)  .AND. choice<prob(5) )THEN;
        CALL ADD_WORMS;                                         ! P-measure verified km=6
      ELSE IF(choice>prob(5)  .AND. choice<prob(6) )THEN;
        CALL REM_WORMS;                                         ! P-measure verified km=7
      ELSE IF(choice>prob(6)  .AND. choice<prob(7) )THEN; CONTINUE
        CALL COMMUTE;                                           ! P-measure verified km=9
      ELSE IF(choice>prob(7)  .AND. choice<prob(8) )THEN; CONTINUE
        CALL MOVEMEAS;                                           ! P-measure verified km=10
      ELSE IF(choice>prob(8)  .AND. choice<prob(9))THEN;
        CALL FLIP_SPINS;                                            ! P-measure verified km=3
      ELSE IF(choice>prob(9)  .AND. choice<prob(10))THEN;
        CALL MOVE_WORMS;                                      ! P-measure verified km=9
      ELSE IF(choice>prob(10) .AND. choice<prob(11) )THEN;
        CALL ADD_PHONON_DIA;                                 ! P-measure verified km=11
      ELSE IF(choice>prob(11) .AND. choice<prob(12) )THEN;
        CALL REM_PHONON_DIA;                                 ! P-measure verified km=12
      ELSE IF(choice>prob(12) .AND. choice<prob(13) )THEN;
        CALL COMMUTE_INT;                                     ! P-measure verified km=13
      ELSE IF(choice>prob(13) .AND. choice<prob(14) )THEN;
        CALL ADD_PHONON_DIA_EXP;                                 ! P-measure verified km=14
      ELSE IF(choice>prob(14) .AND. choice<prob(15) )THEN;
        CALL REM_PHONON_DIA_EXP;                                 ! P-measure verified km=15
      ELSE IF(choice>prob(15) .AND. choice<prob(16) )THEN;
        CALL ADD_PHONON_DIA_SHO;                                 ! P-measure verified km=16
      ELSE IF(choice>prob(16) .AND. choice<prob(17) )THEN;
        CALL REM_PHONON_DIA_SHO;                                 ! P-measure verified km=17
      ENDIF

      IF(severe_check)THEN
           CALL CHECK_CONFIG
           CALL ENERGYTEST; CALL PHASETEST; Energy=ABS(Energy)
           IF(ABS(Energy_a-Energy)/ABS(Energy)>ene_diff)THEN
             CALL DRAW
             PRINT*,"KM =",km;
             PRINT*,"Flo_mont = ",flo_mont
             PRINT*,Energy_a,Energy
             PRINT*,Dphase_a,Dphase
             PRINT*,"measV = ", measV
             PRINT*,"vert 1/2 =",vert1,vert2
             STOP'Energy test violation'
           ENDIF
           IF(Dphase_a /= Dphase)THEN
             CALL DRAW
             PRINT*,"KM =",km; 
             PRINT*,"Flo_mont = ",flo_mont
             PRINT*,Dphase_a,Dphase
             PRINT*,Energy_a,Energy
             STOP'Phase test violation'
           ENDIF
           Energy_a=Energy
      ENDIF

      IF(sup_severe_check)THEN
        PRINT*,"km= ",km,"  flo_mont =",flo_mont
        CALL DRAW
        PAUSE
      ENDIF

! Counting number of worm and no-worm diagrams
      IF(present)THEN; 
          unphysical=unphysical+un1;
          unphysical_c=unphysical_c+un1;
      ELSE; 
          physical=physical+un1; 
          physical_c=physical_c+un1; 
      ENDIF
! Counting when in given order of physical diagram    
      IF(.NOT. present)THEN
         in_order(nmnm) =  in_order(nmnm) + un1
      ENDIF
      order(nmnm) = order(nmnm) + un1
      order_cy(nmnm) = order_cy(nmnm) + un1

      END SUBROUTINE UPDATE
!....................................................................

!--------------------------------------------------------------------
! Perform Monte Carlo
!--------------------------------------------------------------------
      SUBROUTINE to_correct
      USE config_par; USE stat_par; IMPLICIT NONE
      REAL*8 :: zuza

           CALL ENERGYTEST; CALL PHASETEST; Energy=ABS(Energy)
             zuza=ABS(Energy_a-Energy)/ABS(Energy)
           IF(zuza>ene_diff)THEN
             CALL DRAW
             PRINT*,"NO KM possible in to_correct new"
             PRINT*,"Flo_mont = ",flo_mont
             PRINT*,Energy_a,Energy
             PRINT*,"Relative error = ",zuza
             PRINT*,Dphase_a,Dphase
             STOP'Energy test violation'
           ENDIF
           IF(Dphase_a /= Dphase)THEN
             CALL DRAW
             PRINT*,"NO KM possible in to_correct"; 
             PRINT*,"Flo_mont = ",flo_mont
             PRINT*,Dphase_a,Dphase
             PRINT*,Energy_a,Energy
             STOP'Phase test violation'
           ENDIF
           Energy_a=Energy
      
      
      END SUBROUTINE to_correct
!--------------------------------------------------------------------
! Perform Monte Carlo
!--------------------------------------------------------------------
      SUBROUTINE MONTE
      USE stat_par; USE config_par; IMPLICIT NONE
      REAL*8,EXTERNAL :: S_NORM,P_NORM

! Performing Normal Monte Carlo
      DO; flo_mont=flo_mont+1; 
          
        CALL UPDATE; !Always call update
        
         measu_count=measu_count+1; 
         IF(measu_count>=measu_limit)THEN; measu_count=0;
           CALL MEASURE;
           CALL LOOP_MEASURE;
           number_measu=number_measu+un1;
         ENDIF
        
         cor_count=cor_count+1
         IF(cor_count>cor_limit)THEN; cor_count=1
            CALL TO_CORRECT ! Correcting energy test
         ENDIF    
        
        print_count=print_count+1; write_count=write_count+1
         IF(print_count>=print_limit)THEN; print_count=0
           in_print=in_print+1; 
           i_all_meas = i_all_meas +1; 
           CALL PRINT_OUT(0); CALL CHECK_CONFIG
           CALL ZFACTOR
           CALL BOLDG; CALL BOLDPH
           s_norm_now=S_NORM()
           p_norm_now=P_NORM()
           CALL ENERGYTEST; Energy_a=ABS(energy)
           CALL PHASETEST ; Dphase_a=Dphase
         ENDIF

         IF(write_count>=write_limit)THEN; write_count=0
             iskipi_count=iskipi_count+1
           PRINT*," ******    WRITING  D A T A  ************** "
           PRINT*," ******  DO NOT TERMINATE    ************** "
           PRINT*,"                                            "
           PRINT*,"  W A I T   O N E  M I N U T E  ! ! !       "
           PRINT*,"                                            "
           CALL ST_WRITE(0); CALL ST_WRITE(1)
           IF(iskipi_count==iskipi_lim)iskipi_count=0 !statistcs is writen and new iskipi_lim skips started
           CALL CONF_WRITE(0); CALL CONF_WRITE(1)
           CALL RAND_WRITE(0); CALL RAND_WRITE(1)
           CALL CHECK_CONFIG
           PRINT*,"                                            "
           PRINT*," N O W   Y O U   C A N   T E R M I N A T E  "
           PRINT*," But iskipi_count/iskipi_lim = ", iskipi_count, " / ", iskipi_lim
           122   FORMAT("Average loop number = ",ES20.13)        
           PRINT 122,loop_summed / time_summed; 
           PRINT*,"                                            "
         ENDIF

      ENDDO

      END SUBROUTINE monte
!....................................................................


!--------------------------------------------------------------------
! Printing out
!--------------------------------------------------------------------
      SUBROUTINE PRINT_OUT(iswitch)
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: ENER
      INTEGER,INTENT(IN) :: iswitch
      INTEGER :: i, i_site, ix, iy, iz
      REAL*8,DIMENSION(3) :: mom
      REAL*8 :: x, gugu, eeeee, buba, worm_byul, worm_weold
      REAL*8 :: worm_byul_c, worm_byul_s, wo_suge
      REAL*8 :: old_norma, new_norma, new_norma_c, sopal
            
      CALL DRAW

      PRINT*,'{{{{{{{{ DURYNUDA :', DURYUNDA
      
      PRINT'("MC steps in millions = ",ES12.5)',flo_mont
      PRINT*,"nmnm = ",nmnm
! Setting orders of diagrams
      DO i=2,N_max_now,2; 
          in_order_2(i/2)=in_order(i);
             order_2(i/2)=   order(i);
             OOO(i)=order(i)/flo_mont
             IF(OOO(i)<zuka)OOO(i)=zuka
             ord_probab_old(i) = ord_probab(i)
             ord_probab_old_2(i/2) = ord_probab_old(i)
             ord_probab(i) = ord_probab_0(i)/OOO(i)       
             ord_probab_2(i/2) = ord_probab(i)
      ENDDO
! Let us look and fix normalization !!!!!!!!!!!!!!!!!!!!!!! 
      PRINT*,"Let look and fix normas !! "
      PRINT*,"Trivial OOO norm: ",SUM(OOO(2:N_max_now))
      old_norma = SUM(ord_probab_old(2:N_max_now))
      PRINT*,"Old ord_proba norma",old_norma
      new_norma = SUM(ord_probab(2:N_max_now))
      PRINT*,"New ord_proba norma (not correction) ",new_norma
      ord_probab(2:N_max_now) = ord_probab(2:N_max_now) / new_norma
      new_norma_c = SUM(ord_probab(2:N_max_now))
      PRINT*,"New ord_proba norma (yes correction) ",new_norma_c
      DO i=2,N_max_now,2; 
             ord_probab_2(i/2) = ord_probab(i)
      ENDDO
      PRINT*,"Performed    fix normas !! "
      
      OPEN(UNIT=4,FILE="pro_ord.dat")
      OPEN(UNIT=9,FILE="sta_ord.dat")
      DO i=2,N_max_now,2
          WRITE(4,*)i,ord_probab(i)
          WRITE(9,*)i,order(i)
      ENDDO    
      CLOSE(4)
      CLOSE(9)
      
      OPEN(UNIT=9,FILE="loc_ord.dat")
      DO i=2,N_max_now,2
          WRITE(9,*)i,order_cy(i)
      ENDDO    
      CLOSE(9)
      order_cy=0.0d0

      PRINT*,"OLD AND NEW ORDERS PROBABILITIES : "
      PRINT'(10ES12.5," ")',ord_probab_old_2(1:N_max_now/2)
      PRINT'(10ES12.5," ")',ord_probab_2(1:N_max_now/2)
      PRINT*,'Probabilities of diagram orders in physical space : '
      PRINT'(10ES12.5," ")',in_order_2(1:N_max_now/2)/physical
      PRINT*,'Probabilities of diagram orders in total    space : '
      PRINT'(10ES12.5," ")',order_2(1:N_max_now/2) / flo_mont
      
      
      worm_byul=(un1*unphysical)/(physical+unphysical)
      worm_byul_c = (un1*unphysical_c)/(physical_c+unphysical_c)
      unphysical_c=0.0d0; physical_c=0.0d0
      PRINT*,'Worms probabilities : '
      PRINT*," No-worm  : ",1.0-worm_byul, " Ye-worm  : ",worm_byul
      PRINT*," No-worm C: ",1.0-worm_byul_c, " Ye-worm C: ",worm_byul_c
      worm_byul_s = SQRT(worm_byul*worm_byul_c)
      worm_weold=worm_weight

! Changing worm weight when the worm is not present      
      IF(.NOT. present)THEN
        IF(worm_byul_s>0.3d0)THEN;
            wo_suge = worm_weight / 1.5d0
            IF(wo_suge>worm_min)worm_weight=wo_suge
        ELSE IF(worm_byul_s<0.1d0)THEN
            wo_suge = worm_weight * 1.5d0
            IF(wo_suge<worm_max)worm_weight=wo_suge
        ELSE
            CONTINUE
        ENDIF
      ENDIF
      
      PRINT*,"Worm: ",worm_weold," --> ",worm_weight, present
      PRINT*,"MinMax of worm: ",worm_min,worm_max
      PRINT*,"Tau-shift :",suc_tau_shift/(a00_tau_shift+sma), &
            "Site-shift:",suc_site_shift/(a00_site_shift+sma)
      PRINT*,"Add-phonon:",suc_add_pho/(a00_add_pho+sma), &
            "Rem-phonon:",suc_rem_pho/(a00_rem_pho+sma)
      PRINT*,"Add-phoDia:",suc_add_dia/(a00_add_dia+sma), &
            "Rem-phoDia:",suc_rem_dia/(a00_rem_dia+sma)
      PRINT*,"Add-phoExp:",suc_add_dia_exp/(a00_add_dia_exp+sma), &
            "Rem-phoExp:",suc_rem_dia_exp/(a00_rem_dia_exp+sma)
      PRINT*,"Add-phoSho:",suc_add_dia_sho/(a00_add_dia_sho+sma), &
            "Rem-phoSho:",suc_rem_dia_sho/(a00_rem_dia_sho+sma)
      PRINT*,"Add-worm  :",suc_add_worm/(a00_add_worm+sma), &
            "Rem-worm  :",suc_rem_worm/(a00_rem_worm+sma)
      PRINT*,"Flip-spin :",suc_flip_spi/(a00_flip_spi+sma), &
            "Move-worm :",suc_move_worm/(a00_move_worm+sma)
      PRINT*,"Commute   :",suc_commute/(a00_commute+sma), &
            "CommuteInt:",suc_commute_int/(a00_commute_int+sma), &
            "Move_meas :",suc_move_mea/(a00_move_mea+sma)
      PRINT*,"Measure vertexes:", vert1,vert2,measV
      PRINT*,"OCC:",stat_cut(1)/stat_cut(0), stat_cut(2)/stat_cut(0)
      PRINT*,"OCC:",stat_cut(3)/stat_cut(0), stat_cut(4)/stat_cut(0)
      PRINT*,"OCC:",stat_cut(5)/stat_cut(0)

      CALL DENSITY_CHECK     

!Checking in any case
           CALL CHECK_CONFIG
           CALL ENERGYTEST; CALL PHASETEST;
           PRINT*,' Dphase = ',Dphase,'   Energy =',Energy
           Energy=ABS(Energy)
           gugu=ABS(Energy_a-Energy)/ABS(Energy)
           PRINT*,"Relative Energy mismatch = ",gugu
           IF(gugu>ene_diff)THEN
             CALL DRAW
             PRINT*,"KM =",km;
             PRINT*,"Flo_mont = ",flo_mont
             PRINT*,Energy_a,Energy
             STOP'Energy test violation'
           ENDIF
           IF(Dphase_a /= Dphase)THEN
             PRINT*,"KM =",km;
             STOP'Phase test violation'
           ENDIF
           Energy_a=Energy
!End checking in any case

      call GIVENK_0
      OPEN(4,file='nk_a0.dat')
	i=0	 
      do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1
	mom(1)=(pi_m2/L(1))*ix ; mom(2)=(pi_m2/L(2))*iy ;
	mom(3)=(pi_m2/L(3))*iz ; x=ENER(1,mom)
!	IF(ABS(x)<10.d0) then
      write(4,*) x, nk_distr_0(ix,iy,iz)
	i=i+1
!	ENDIF
      enddo; enddo; enddo
      CLOSE(4)
 !     print*, 'number of x points in nk=', i
      
      PRINT*,"Before givenk"
      call GIVENK
      OPEN(4,file='nk_a.dat')
	i=0	 
      do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1
	mom(1)=(pi_m2/L(1))*ix ; mom(2)=(pi_m2/L(2))*iy ;
	mom(3)=(pi_m2/L(3))*iz ; x=ENER(1,mom)
!	IF(ABS(x)<10.d0) then
      write(4,*) x, nk_distr(ix,iy,iz)
	i=i+1
!	ENDIF
      enddo; enddo; enddo
      CLOSE(4)
      print*, 'number of x points in nk=', i
 
      OPEN(4,FILE="n_fer.dat")
      DO i=0,1000; eeeee = -mu + ((4.0*dim*hopping)/1000.0)*i
         buba = 1.0d0 / (1.0d0 + exp(eeeee*beta))
         WRITE(4,*)eeeee,buba
      ENDDO   
      CLOSE(4)
 
! ONLY FOR 2D      
      IF(dim==2)THEN      
      OPEN(UNIT=4,FILE="n_ot_k.dat")
      DO ix=0,L(1)/2; mom(1)=(pi_m2/L(1))*ix
!          PRINT*,nk_distr(ix,0,0),nk_distr(L(1)-ix,0,0)
          WRITE(4,*)mom(1),nk_distr(ix,0,0),nk_distr_0(ix,0,0)
      ENDDO    
      CLOSE(4)
      ENDIF
      
      OPEN(UNIT=4,FILE="fer.dat")
      OPEN(UNIT=9,FILE="fer0.dat")
      DO ix=0,L(1)-1; mom(1)=(pi_m2/L(1))*ix
      DO iy=0,L(2)-1; mom(2)=(pi_m2/L(2))*iy
        IF(ABS(nk_distr(ix,iy,0)-0.5d0)<0.2)THEN
            WRITE(4,*)mom(1),mom(2)
        ENDIF    
        IF(ABS(nk_distr_0(ix,iy,0)-0.5d0)<0.2)THEN
            WRITE(9,*)mom(1),mom(2)
        ENDIF    
      ENDDO; ENDDO;
      CLOSE(4); CLOSE(9)
      !ONLY      
      
      
      END SUBROUTINE print_out
!....................................................................

!--------------------------------------------------------------------
! Testing output for self-energy and polarization operator
!--------------------------------------------------------------------
      SUBROUTINE TEST_SE_PO
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,S_NORM,P_NORM,FUN_sigma,FUN_polar
      REAL*8,DIMENSION(0:10000) :: SiSi4, PiPi4
      REAL*8 :: tt, tt1, tt2, caballa
      INTEGER :: it, iii,ii1,ii2,ii_tut

127   FORMAT(1x, 10(E20.10,1x))      
      
      IF(N_mes_now>2)PRINT*,"NO comparison with 2nd order NOW"

! Evaluating 4th order of Sigma for local (only once at beginning)
      IF(Only_when_begins)THEN
      SiSi4(0:Ntau-1)=0.0d0
!      DO iii=0,Ntau-1; tt=(iii+half)*dtau
!        DO ii1=0,2*Ntau-1; DO ii2=0,4*Ntau-1
!          tt1=(ii1+half)*dtau/2.; tt2=(ii2+half)*dtau/4.    
!          SiSi4(iii)=SiSi4(iii)+ &
!                 GREEN(1,nul,1,tt1,1) * GREEN(1,tt1,1,tt2,1) * &
!                 GREEN(1,tt2,1,tt,1)                         * &
!                 PHONON(1,nul,1,tt2)  * PHONON(1,tt1,1,tt)      * &
!                 dtau*dtau/2./4.           
!        ENDDO; ENDDO    
!      ENDDO
      ENDIF
! End of evaluating 4th order of Sigma for local  
          
! Testing output for self-energy
      S_norm_now=S_NORM()
      PRINT*,'S_norm:',S_norm_now
      OPEN(UNIT=9,FILE="che_sigma.dat")
      DO it=0,Ntau-1; tt = (it+half)*dtau
        WRITE(9,127)tt, Sigma_box(0,0,0,it)*S_norm_now/(Z_self_norm*dtau), &
                     FUN_sigma(0,0,0,tt)*S_norm_now/(Z_self_norm), &
                     GREEN(1,nul,1,tt,1)*PHONON(1,nul,1,tt)  
        IF(it==0)THEN
        PRINT*,'SELFE:', ( Sigma_box(0,0,0,it)*S_norm_now/(Z_self_norm*dtau) ) / &
                                 ( GREEN(1,nul,1,tt,1)*PHONON(1,nul,1,tt) )
        sugma_nac(in_print) = Sigma_box(0,0,0,it)*S_norm_now/(Z_self_norm*dtau)
        ENDIF
      ENDDO
      CLOSE(9)
      OPEN(UNIT=9,FILE="4_sigma.dat")
      DO it=0,Ntau-1; tt = (it+half)*dtau
        WRITE(9,127)tt,SiSi4(it)   
      ENDDO    
      CLOSE(4)
! Evaluating 4th order of P for local 
      Only_when_begins=.FALSE.
      IF(Only_when_begins)THEN
      PiPi4(0:Ntau-1)=0.0d0
 !     DO iii=0,Ntau-1; tt=(iii+half)*dtau
 !       DO ii_tut=1,Nsite
 !       DO ii1=0,2*Ntau-1; DO ii2=0,4*Ntau-1
 !         tt1=(ii1+half)*dtau/2.; tt2=(ii2+half)*dtau/4. 
 !         IF(ii_tut==1)THEN
 !             caballa=2.0d0
 !         ELSE
 !             caballa=2.0d0
 !         ENDIF    
 !         PiPi4(iii)=PiPi4(iii)+ &
 !          GREEN(1,nul,ii_tut,tt1,1) * GREEN(ii_tut,tt1,1,tt,1)  * &
 !          GREEN(1,tt,ii_tut,tt2,1)  * GREEN(ii_tut,tt2,1,nul,1) * &
 !          PHONON(1,tt1,1,tt2)  * caballa * dtau*dtau /(8.0d0)           
 !       ENDDO; ENDDO    
 !       ENDDO
 !     ENDDO
      ENDIF
! End of evaluating 4th order of P for local 
      
! Testing output for polarization operator
      P_norm_now=P_NORM()
      PRINT*,'P_norn:',P_norm_now
      OPEN(UNIT=9,FILE="che_polar.dat")
      DO it=0,Ntau-1; tt = (it+half)*dtau
        IF(.NOT. POLARIZED)THEN  
        ! NON POL_SPIN            
          WRITE(9,127)tt, Polar_box(0,0,0,it)*P_norm_now/(Z_pola_norm*dtau), &
                              FUN_polar(0,0,0,tt)*P_norm_now/(Z_pola_norm), &
                             -2.0d0*GREEN(1,nul,1,tt,1)*GREEN(1,tt,1,nul,1)       
          IF(it==0)THEN
           PRINT*,'POLAR:', ( Polar_box(0,0,0,it)*P_norm_now/(Z_pola_norm*dtau) ) / &
                                     ( 2.0d0*GREEN(1,nul,1,tt,1)*GREEN(1,tt,1,nul,1) )
           pilar_nac(in_print)= Polar_box(0,0,0,it)*P_norm_now/(Z_pola_norm*dtau)
          ENDIF
        ! NON POL_SPIN         
        ELSE
        ! POL_SPIN------
          WRITE(9,127)tt, Polar_box(0,0,0,it)*P_norm_now/(Z_pola_norm*dtau), &
                              FUN_polar(0,0,0,tt)*P_norm_now/(Z_pola_norm), &
                             -GREEN(1,nul,1,tt,1)*GREEN(1,tt,1,nul,1)       
          IF(it==0)THEN
           PRINT*,'POLAR:', ( Polar_box(0,0,0,it)*P_norm_now/(Z_pola_norm*dtau) ) / &
                                     ( GREEN(1,nul,1,tt,1)*GREEN(1,tt,1,nul,1) )
           pilar_nac(in_print)= Polar_box(0,0,0,it)*P_norm_now/(Z_pola_norm*dtau)
          ENDIF
        ! POL_SPIN------
        ENDIF
      ENDDO
      CLOSE(9)
      
      wowe(in_print)=worm_weight
      
      OPEN(UNIT=4,FILE="sugma.dat"); OPEN(UNIT=5,FILE="pilar.dat")
      OPEN(UNIT=10,FILE="wove.dat")
        DO iii=1,in_print
            WRITE(4,*)iii,sugma_nac(iii)
            WRITE(5,*)iii,pilar_nac(iii)
            WRITE(10,*),iii,wowe(iii)
        ENDDO    
      CLOSE(4); CLOSE(5); CLOSE(10);
      
      OPEN(UNIT=9,FILE="n_gr.dat"); OPEN(UNIT=10,FILE="kineti.dat"); 
      OPEN(UNIT=4,FILE="m_x.dat"); OPEN(UNIT=5,FILE="m_xy.dat"); 
      OPEN(UNIT=11,FILE="Z_fac_fac.dat");
         DO iii=1,i_all_meas      
            WRITE(9,*)iii,dunsity_ST(iii)      
            WRITE(10,*)iii,kineti_ST(iii)      
            WRITE(4,*)iii,m_x_ST(iii)      
            WRITE(5,*)iii,m_xy_ST(iii)      
            WRITE(11,*)iii,Z_fac_ST(iii)      
         ENDDO
      CLOSE(9); CLOSE(10); CLOSE(4); CLOSE(5); CLOSE(11) 
            
      Only_when_begins =  .FALSE.

      END SUBROUTINE test_se_po
!....................................................................

!--------------------------------------------------------------------
! Writing statistics
!--------------------------------------------------------------------
      SUBROUTINE ST_WRITE(fiku)
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,BO_OC
      INTEGER,INTENT(IN) :: fiku
      CHARACTER*7 :: fil
      INTEGER :: iii
      INTEGER :: ix, iy, iz, it, kb, ib, isapa, ii_mee
      REAL*8 :: x, dobaka
      
! Preparing for writing order-probabilitied for possible change of order at next read
      dobaka = ord_probab(N_max_now)/lupaka
      PRINT*,"For higher orders probability was added equal to ",dobaka
      ord_probab(N_max_now+1:N_max_max) = dobaka
! Adding statistics for the whole rest of tail for further re-switch
      order(N_max_now+1:N_max_max) = srednyaka
      OPEN(UNIT=9,FILE="ord.dat")
      DO iii=2,N_max_now+4,2
         WRITE(9,*)iii,order(iii)
      ENDDO   
      CLOSE(9)

      IF(fiku==0)CALL TEST_SE_PO
            
      PRINT*," *********  WRITING STATISTICS ********** "
! Defining to which file to write
      IF(fiku==0)     THEN; fil='sta.out'
      ELSE IF(fiku==1)THEN; fil='sta.ou1'
      ELSE;           STOP'Where to write statistics'
      ENDIF 

      PRINT*,"WRITING NOW: ",fil
      
      IF(iskipi_count==iskipi_lim)THEN
      OPEN(UNIT=4,FILE=fil)
      WRITE(4,*)flo_mont
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1; 
         WRITE(4,*)Sigma_box(ix,iy,iz,it), Polar_box(ix,iy,iz,it) ! prime data 
      ENDDO; ENDDO; ENDDO; ENDDO
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; 
        DO kb=0,Nbins_s-1; DO ib=1,Nbasis_s;
          WRITE(4,*)Sigmab(ix,iy,iz,kb,ib) 
        ENDDO; ENDDO; 
        DO kb=0,Nbins_p-1; DO ib=1,Nbasis_p;
          WRITE(4,*)Polarb(ix,iy,iz,kb,ib) 
        ENDDO; ENDDO; 
        DO kb=0,Nbins_OC-1; DO ib=1,Nbasis_OC;
          WRITE(4,*)OCb(ix,iy,iz,kb,ib) 
        ENDDO; ENDDO; 
      ENDDO; ENDDO; ENDDO; 
      WRITE(4,*)S_norm_now,P_norm_now   !AnalytForNorm
      WRITE(4,*)Z_self_norm,Z_pola_norm !MC_Counters 
      WRITE(4,*)order(1:N_max_max)
      WRITE(4,*)ord_probab(1:N_max_max)
      WRITE(4,*)Ntau,i_cur_meas
!      DO isapa=1,momsko;
!          DO it=0,Ntau-1;
!              DO ii_mee=1,i_cur_meas
!                 WRITE(4,*)covar(isapa,it,ii_mee)    
!              ENDDO        
!          ENDDO    
!      ENDDO    
      WRITE(4,*)Z_TokTok_norm ! <JJ> counter
      WRITE(4,*)TokTok_box(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1)
      WRITE(4,*)i_all_meas
      WRITE(4,*)kineti_ST(1:i_all_meas)
      WRITE(4,*)dunsity_ST(1:i_all_meas)
      WRITE(4,*)Z_fac_ST(1:i_all_meas)
      WRITE(4,*)m_x_ST(1:i_all_meas)
      WRITE(4,*)m_xy_ST(1:i_all_meas)
      WRITE(4,*)num_oc
      WRITE(4,*)jj_exa(0:num_oc*2)
      DO kb=0,skoko_cut; 
          DO ib=0,num_oc*2;
              WRITE(4,*)jj_exa_cut(kb,ib)
          ENDDO; 
      ENDDO;    
      WRITE(4,*)num_oc_ano;
      WRITE(4,*)jj_exa_ano(0:num_oc_ano*2)
      DO kb=0,skoko_cut; 
          DO ib=0,num_oc_ano*2;
              WRITE(4,*)jj_exa_cut_ano(kb,ib)
          ENDDO; 
      ENDDO;    
!      DO it=0,num_oc_lim;
!          DO ii_mee=1,i_cur_meas
!              WRITE(4,*)covar_oc(it,ii_mee)    
!          ENDDO        
!      ENDDO    
      WRITE(4,*)loop_summed, time_summed
      CLOSE(4)
      ELSE
        PRINT*,"Skipped writing statustics ratio: ", (1.0d0*iskipi_count)/(1.0d0*iskipi_lim)  
      ENDIF    

      IF(MEASURE_OC)THEN; ! Writing only if measuring
          IF(fiku==0)THEN; ! Writing only once when fiku==1
               CALL PRESENTING_OC_OUTPUT;
               OPEN(UNIT=4,FILE="optico_histo.dat")
                  DO iii=-num_aga,num_aga;
                     WRITE(4,*)aga_hist(iii) 
                  ENDDO    
               CLOSE(4)
          ENDIF
      ENDIF
     
      PRINT*,"Statistics is written"

      
      END SUBROUTINE st_write
!....................................................................
      
!--------------------------------------------------------------------
! Presenting OC outout
!--------------------------------------------------------------------
      SUBROUTINE PRESENTING_OC_OUTPUT
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      REAL*8,EXTERNAL :: BO_OC
      REAL*8,ALLOCATABLE,DIMENSION(:) :: tOtOsHA
      INTEGER :: ihih, jhjh, ika, iii, iba
      INTEGER :: it, ib, kb, ix, iy, iz
      REAL*8 :: basodur, inte_jj, tt, x, noruter

! Summing contributions      
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; 
      DO it=0,Ntau-1; 
      TokTok_box_all(ix,iy,iz,it) = TokTok_box(ix,iy,iz,it) + TokTok_box(ix,iy,iz,Ntau-1-it) 
      ENDDO
      ENDDO; ENDDO; ENDDO
      TokTok_box_all = TokTok_box_all / (Z_TokTok_norm*dtau*8.0d0*dim)  

! Collecting and writing data for OC for histogram collection
!-------------------------------------------------------------------
      ALLOCATE(tOtOsHA(0:Ntau-1))
      tOtOsHA=nul   
      IF(dim==1)THEN
         DO it=0,Ntau-1; DO ihih=0,L(1)-1;  
           tOtOsHA(it) =  tOtOsHA(it) + TokTok_box_all(ihih,0,0,it)  
         ENDDO; ENDDO
         tOtOsHA = tOtOsHA  * ( beta * (dir**2) )
      ELSE IF(dim==2)THEN
         DO it=0,Ntau-1; DO ihih=0,L(1)-1; DO jhjh=0,L(2)-1; 
           tOtOsHA(it) =  tOtOsHA(it) + TokTok_box_all(ihih,jhjh,0,it)  
         ENDDO; ENDDO; ENDDO    
         tOtOsHA = tOtOsHA  * ( beta * (dir**2) )
      ELSE
          STOP'No estimator fo dim>2'
      ENDIF    
      
      inte_jj=0.0d0
      DO it=0,Ntau-1; tt=(it+0.5)*dtau
         inte_jj  = inte_jj + tOtOsHA(it)*dtau
      ENDDO    
      basodur = inte_jj/uh_kinetik1
      PRINT*," From JJ integral :",inte_jj,basodur
      
      OPEN(UNIT=4,FILE="totosha.dat")
      DO it=0,Ntau-1; tt=(it+0.5)*dtau
         WRITE(4,*)tt,tOtOsHA(it)/basodur
      ENDDO    
      CLOSE(4)
      
      !quantities to check estimators
      !findinf norms
      Nor_TT(1:4)=nul
      DO ika=1,4
          DO it=0,Ntau-1
              Nor_TT(ika) = Nor_TT(ika) + (ToTo_box_sep(it,ika)*dtau)
          ENDDO 
              Nor_TT(ika) = Nor_TT(ika) / uh_kinetik1
      ENDDO    
      ! normalizing
      DO ika=1,4; 
          ToTo_box_sep_no(0:Ntau-1,ika) = ToTo_box_sep(0:Ntau-1,ika) / Nor_TT(ika); 
      ENDDO
      ! writing    
102   FORMAT(1x,(10E14.6, 1x))
      OPEN(UNIT=4,FILE="box_check.dat")
          DO it=0,Ntau-1; tt=(it+0.5)*dtau
              WRITE(4,102)tt,ToTo_box_sep_no(it,1:4)
          ENDDO    
      CLOSE(4)
!------------------------------------------------------------------      
! End Collecting and writing data for OC for histogram collection      
      
! Collecting and writing data for OC for basis finctions  collection
!-------------------------------------------------------------------      
      ! Converting basis to grid for OC      
       DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1;
         tt=dtau*(it+0.5d0);  
         kb=tt/dtb_OC 
         x=0.d0; do ib=1,Nbasis_OC
         x=x+OCb(ix,iy,iz,kb,ib)*bo_OC(ib,tt,kb); enddo     
         ToTo_Basis(ix,iy,iz,it)=x*dtau; 
       enddo; enddo; enddo; enddo     ! basis converted to grid for Pi

      tOtOsHA=nul   
      IF(dim==1)THEN
         DO it=0,Ntau-1; DO ihih=0,L(1)-1;  
           tOtOsHA(it) =  tOtOsHA(it) + ToTo_Basis(ihih,0,0,it) + ToTo_Basis(ihih,0,0,Ntau-1-it) 
         ENDDO; ENDDO
         tOtOsHA = tOtOsHA  * ( beta * (dir**2) )
      ELSE IF(dim==2)THEN
         DO it=0,Ntau-1; DO ihih=0,L(1)-1; DO jhjh=0,L(2)-1; 
           tOtOsHA(it) =  tOtOsHA(it) + ToTo_Basis(ihih,jhjh,0,it) + ToTo_Basis(ihih,jhjh,0,Ntau-1-it)
         ENDDO; ENDDO; ENDDO    
         tOtOsHA = tOtOsHA  * ( beta * (dir**2) )
      ELSE
          STOP'No estimator fo dim>2'
      ENDIF    
      
      inte_jj=0.0d0
      DO it=0,Ntau-1; tt=(it+0.5)*dtau
         inte_jj  = inte_jj + tOtOsHA(it)*dtau
      ENDDO    
      basodur = inte_jj/uh_kinetik1
      PRINT*," From JJ integral (basis) :",inte_jj,basodur
      
      OPEN(UNIT=4,FILE="totosha_basis.dat")
      DO it=0,Ntau-1; tt=(it+0.5)*dtau
         WRITE(4,*)tt,tOtOsHA(it)/basodur
      ENDDO    
      CLOSE(4)
       
!------------------------------------------------------------------      
! End Collecting and writing data for OC for histogram collection  
       
! Collecting and writing data for OC for exact points estimator
!-------------------------------------------------------------------      
      inte_jj=0.0d0
      DO it=0,num_oc_lim-1; 
          inte_jj  = inte_jj + 0.5d0 * (jj_exa(it)+jj_exa(it+1)) * (ta_oc(it+1)-ta_oc(it))
      ENDDO    
      basodur = inte_jj/uh_kinetik1
      PRINT*," From JJ integral (exact) :",inte_jj,basodur
      
      DO it=0,num_oc_lim; 
          jj_exa_norm(it) = (jj_exa(it)+jj_exa(num_oc_lim-it))*0.5d0
          jj_exa_norm(it) = jj_exa_norm(it) / basodur
          DO iba=0,skoko_cut
               jj_exa_cut_norm(iba,it) = (jj_exa_cut(iba,it)+jj_exa_cut(iba,num_oc_lim-it))*0.5d0
               jj_exa_cut_norm(iba,it) = jj_exa_cut_norm(iba,it) / basodur
          ENDDO    
      ENDDO    
      
      
      OPEN(UNIT=4,FILE="totosha_exa.dat")
      DO it=0,num_oc_lim; tt=ta_oc(it)
         WRITE(4,*)tt,jj_exa_norm(it)
      ENDDO    
      CLOSE(4)
      
      DO iba=0,skoko_cut;
         OPEN(UNIT=4,FILE=files_cut(iba))
         DO it=0,num_oc_lim; tt=ta_oc(it)
            WRITE(4,*)tt,jj_exa_cut_norm(iba,it)
         ENDDO;    
         CLOSE(4)
      ENDDO;
                
      
      
! QUANTITIES T CHECK NORMS: quantities to check estimators
      !findinf norms
      Nor_jj(1:4)=nul
      DO ika=1,4
          DO it=0,num_oc_lim-1
              Nor_jj(ika) = Nor_jj(ika) + & 
               0.5d0 * (jj_exa_sep(it,ika)+jj_exa_sep(it+1,ika)) * (ta_oc(it+1)-ta_oc(it))
          ENDDO 
              Nor_jj(ika) = Nor_jj(ika) / uh_kinetik1
      ENDDO    
      ! normalizing
      DO it=0,num_oc_lim; 
          DO ika=1,4; 
              jj_exa_sep_no(it,ika) = & 
              0.5d0 * ( jj_exa_sep(it,ika) + jj_exa_sep(num_oc_lim-it,ika) ) / Nor_jj(ika); 
          ENDDO
      ENDDO
      OPEN(UNIT=4,FILE="exa_check.dat")
          DO it=0,num_oc_lim; 
              WRITE(4,102)ta_oc(it),jj_exa_sep_no(it,1:4)
          ENDDO    
      CLOSE(4)
      
!------------------------------------------------------------------      
! End Collecting and writing data for OC for exact point estimator

! Collecting and writing data for OC for ANOTHER exact points estimator
!-------------------------------------------------------------------      
      !basodur = inte_jj/uh_kinetik1
      ! Here we use basodur for normalization from complete [0,beta] <jj> 
      
      DO it=0,num_oc_lim_ano; 
          jj_exa_norm_ano(it) = (jj_exa_ano(it)+jj_exa_ano(num_oc_lim_ano-it))*0.5d0
          jj_exa_norm_ano(it) = jj_exa_norm_ano(it) / basodur
          DO iba=0,skoko_cut
               jj_exa_cut_norm_ano(iba,it) = (jj_exa_cut_ano(iba,it)+jj_exa_cut_ano(iba,num_oc_lim_ano-it))*0.5d0
               jj_exa_cut_norm_ano(iba,it) = jj_exa_cut_norm_ano(iba,it) / basodur
          ENDDO    
      ENDDO    
      
      
      OPEN(UNIT=4,FILE="totosha_exa_ano.dat")
      DO it=0,num_oc_lim_ano; tt=ta_oc_ano(it)
         WRITE(4,*)tt,jj_exa_norm_ano(it)
      ENDDO    
      CLOSE(4)
                          
      DO iba=0,skoko_cut;
         OPEN(UNIT=4,FILE=files_cut_ano(iba))
         DO it=0,num_oc_lim_ano; tt=ta_oc_ano(it)
            WRITE(4,*)tt,jj_exa_cut_norm_ano(iba,it)
         ENDDO;    
         CLOSE(4)
      ENDDO;
!------------------------------------------------------------------      
! End Collecting and writing data for OC for exact point estimator


      ! For OC contribution control
      noruter = SUM( aga_hist(-num_aga:num_aga) ) * (aga_value(1)-aga_value(0))
      OPEN(UNIT=4,FILE="asensio.dat")
      DO iii=-num_aga,num_aga;
          WRITE(4,*)aga_value(iii),aga_hist(iii)/noruter
      ENDDO
      CLOSE(4)
      
      DEALLOCATE(tOtOsHA)

      END SUBROUTINE presenting_oc_output
!....................................................................
      
      
!--------------------------------------------------------------------
! Writing statistics
!--------------------------------------------------------------------
      SUBROUTINE ST_READ(ping)
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: ping
      CHARACTER*7 :: fil
      INTEGER :: ix, iy, iz, it, kb, ib, isapa, ii_mee, i_divide_meas
      INTEGER :: i_divide_all, iii
      REAL*8 :: normal_old, normal_new

      IF(ping==0)THEN;      fil='sta.out'
      ELSE IF(ping==1)THEN; fil='sta.ou1'
      ELSE; STOP'NoNo ST'
      ENDIF

      PRINT*,' READING NOW: ',fil
      
      OPEN(UNIT=4,FILE=fil)
      READ(UNIT=4,FMT=*,ERR=2,END=2)flo_mont
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau-1; 
         READ(UNIT=4,FMT=*,ERR=2,END=2)Sigma_box(ix,iy,iz,it), Polar_box(ix,iy,iz,it) ! prime data 
      ENDDO; ENDDO; ENDDO; ENDDO
      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; 
        DO kb=0,Nbins_s-1; DO ib=1,Nbasis_s;
          READ(UNIT=4,FMT=*,ERR=2,END=2)Sigmab(ix,iy,iz,kb,ib) 
        ENDDO; ENDDO; 
        DO kb=0,Nbins_p-1; DO ib=1,Nbasis_p;
          READ(UNIT=4,FMT=*,ERR=2,END=2)Polarb(ix,iy,iz,kb,ib) 
        ENDDO; ENDDO; 
        DO kb=0,Nbins_OC-1; DO ib=1,Nbasis_OC;
          READ(UNIT=4,FMT=*,ERR=2,END=2)OCb(ix,iy,iz,kb,ib) 
        ENDDO; ENDDO; 
      ENDDO; ENDDO; ENDDO; 
       READ(UNIT=4,FMT=*,ERR=2,END=2)S_norm_now,P_norm_now   !AnalyticScalesForNormalization
      READ(UNIT=4,FMT=*,ERR=2,END=2)Z_self_norm,Z_pola_norm !MC_Counters 
      READ(UNIT=4,FMT=*,ERR=2,END=2)order(1:N_max_max)
      READ(UNIT=4,FMT=*,ERR=2,END=2)ord_probab(1:N_max_max)
      READ(UNIT=4,FMT=*,ERR=2,END=2)Ntau,i_cur_meas
!      DO isapa=1,momsko;
!          DO it=0,Ntau-1;
!              DO ii_mee=1,i_cur_meas
!                 READ(UNIT=4,FMT=*,ERR=2,END=2)covar(isapa,it,ii_mee)
!              ENDDO        
!          ENDDO    
!      ENDDO          
      READ(UNIT=4,FMT=*,ERR=2,END=2)Z_TokTok_norm ! <JJ> counter
      READ(UNIT=4,FMT=*,ERR=2,END=2)TokTok_box(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau-1)
      READ(UNIT=4,FMT=*,ERR=2,END=2)i_all_meas
      READ(UNIT=4,FMT=*,ERR=2,END=2)kineti_ST(1:i_all_meas)
      READ(UNIT=4,FMT=*,ERR=2,END=2)dunsity_ST(1:i_all_meas)
      READ(UNIT=4,FMT=*,ERR=2,END=2)Z_fac_ST(1:i_all_meas)
      READ(UNIT=4,FMT=*,ERR=2,END=2)m_x_ST(1:i_all_meas)
      READ(UNIT=4,FMT=*,ERR=2,END=2)m_xy_ST(1:i_all_meas)
      READ(UNIT=4,FMT=*,ERR=2,END=2)num_oc
      READ(UNIT=4,FMT=*,ERR=2,END=2)jj_exa(0:num_oc*2)
      DO kb=0,skoko_cut; 
          DO ib=0,num_oc*2;
              READ(UNIT=4,FMT=*,ERR=2,END=2)jj_exa_cut(kb,ib)
          ENDDO; 
      ENDDO;    
      READ(UNIT=4,FMT=*,ERR=2,END=2)num_oc_ano;
      READ(UNIT=4,FMT=*,ERR=2,END=2)jj_exa_ano(0:num_oc_ano*2)
      DO kb=0,skoko_cut; 
          DO ib=0,num_oc_ano*2;
              READ(UNIT=4,FMT=*,ERR=2,END=2)jj_exa_cut_ano(kb,ib)
          ENDDO; 
      ENDDO;    
 !     DO it=0,num_oc_lim;
 !        DO ii_mee=1,i_cur_meas
 !             READ(UNIT=4,FMT=*,ERR=2,END=2)covar_oc(it,ii_mee)
 !        ENDDO        
 !     ENDDO          
      READ(UNIT=4,FMT=*,ERR=2,END=2)loop_summed, time_summed
      CLOSE(4)
      flo_mont=flo_mont/sevide; order=order/sevide
      Sigma_box=Sigma_box/sevide; Polar_box=Polar_box/sevide;
      Sigmab=Sigmab/sevide; Polarb=Polarb/sevide; OCb=OCb/sevide
      Z_self_norm=Z_self_norm/sevide; Z_pola_norm=Z_pola_norm/sevide;
      TokTok_box = TokTok_box / (sevide*opvide)
      Z_TokTok_norm = Z_TokTok_norm/ (sevide*opvide)
      jj_exa = jj_exa / (sevide*opvide); jj_exa_cut = jj_exa_cut/ (sevide*opvide);  
      jj_exa_ano = jj_exa_ano / (sevide*opvide); jj_exa_cut_ano = jj_exa_cut_ano/ (sevide*opvide);
      ! cutting statistics for covariant matrix
      i_divide_meas = i_cur_meas / sevide
!      DO isapa=1,momsko; DO it=0,Ntau-1;
!          covar(isapa, it, 1:i_divide_meas) = covar(isapa, it, i_cur_meas-i_divide_meas+1:i_cur_meas)        
!      ENDDO; ENDDO          
!      DO it=0,num_oc_lim;
!          covar_oc(it, 1:i_divide_meas) = covar_oc(it, i_cur_meas-i_divide_meas+1:i_cur_meas)        
!      ENDDO;          
      PRINT*,"Length of covariance masure changed ",i_cur_meas,"  --->  ",i_divide_meas
      i_cur_meas=i_divide_meas
      ! cutting statistics of all useful measurements
      i_divide_all = ( i_all_meas / sevide )
      kineti_ST(1:i_divide_all) = kineti_ST(i_all_meas-i_divide_all+1:i_all_meas)
      dunsity_ST(1:i_divide_all) = dunsity_ST(i_all_meas-i_divide_all+1:i_all_meas)
      Z_fac_ST(1:i_divide_all) = Z_fac_ST(i_all_meas-i_divide_all+1:i_all_meas)
      m_x_ST(1:i_divide_all) = m_x_ST(i_all_meas-i_divide_all+1:i_all_meas)
      m_xy_ST(1:i_divide_all) = m_xy_ST(i_all_meas-i_divide_all+1:i_all_meas)
      PRINT*,"Length of all useful  masure changed ",i_all_meas,"  --->  ",i_divide_all
      i_all_meas=i_divide_all
      loop_summed = loop_summed / sevide; time_summed = time_summed /sevide;
      
      GOTO 1
 2    ping=ping+1
1     CONTINUE
      
! Renormalizing ord_probab if maximal order is changed

      OPEN(UNIT=4,FILE="ord_do.dat")
          DO iii=2,N_max_now,2
              WRITE(4,*)iii,ord_probab(iii)
          ENDDO    
      CLOSE(4)

      normal_old = SUM(ord_probab(1:N_max_now))
      PRINT*,"Norma of order probability after   reading is ",normal_old
      ord_probab(1:N_max_now) = ord_probab(1:N_max_now) / normal_old  
      normal_new = SUM(ord_probab(1:N_max_now))
      PRINT*,"Norma of order probability after renorming is ",normal_new
      
      OPEN(UNIT=4,FILE="ord_po.dat")
          DO iii=2,N_max_now,2
              WRITE(4,*)iii,ord_probab(iii)
          ENDDO    
      CLOSE(4)

      PRINT*,"Statistics is read"

     OPEN(UNIT=4,FILE="optico_histo.dat")
        DO iii=-num_aga,num_aga;
          READ(UNIT=4,FMT=*,ERR=19,END=19)aga_hist(iii) 
        ENDDO    
     CLOSE(4)
     GOTO 199
19   aga_hist(-num_aga:num_aga)=0.0d0;
     PRINT*,"I hade to nullify aga_hist!"
199  CONTINUE     

!      PRINT*,POLAR(0:L(1)-1,0:L(2)-1,0,1)
      
      END SUBROUTINE st_read
!....................................................................

!--------------------------------------------------------------------
! Writing configuration
!--------------------------------------------------------------------
      SUBROUTINE CONF_WRITE(fiku)
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: fiku
      INTEGER*4,DIMENSION(1:522) :: isave
      INTEGER :: i, j, jl, name
      CHARACTER*7 :: fil

      PRINT*," *********  WRIGHTING CONFIGURATION ********** "
! Defining to which file to write
      IF(fiku==0)     THEN; fil='cnf.dat'
      ELSE IF(fiku==1)THEN; fil='cnf.da1'
      ELSE;           STOP'Where to write statistics'
      ENDIF 

      PRINT*,"WRITING NOW: ",fil

      OPEN(UNIT=4,FILE=fil)
! Writing the counters
      WRITE(4,*)nmnm
      WRITE(4,*)MEASV,vert1,vert2
      WRITE(4,*)present,sveta,tonya,omegaST(1:3)
      WRITE(4,*)Energy_a,Dphase_a
! Writing vertexes para,eters
      DO i=1,nmnm; name=nlist(i)
        WRITE(4,*)name  
        WRITE(4,*)link(name,1:3), type(name,1:3), tau(name),site(name)
        WRITE(4,*) omega(name,1,1:3), omega(name,2,1:3), omega(name,3,1:3)   
      ENDDO    
      WRITE(4,*)nlist(1:N_max),numnam(1:N_max)    
! Writing hash table
      DO j=-2,N_max
      	write(4,*) linkh2(j,1),linkh2(j,2) ! hash links
	ENDDO
      DO j=1,N_max
      	write(4,*) linkh3(j,1),linkh3(j,2) ! hash links
	ENDDO
      DO j=1,Nhash                                                ! hash lists 
          write(4,*) ha2(j), ha3(j)
	DO jl=1,Nhlist
      	write(4,*) hl2(j,jl), hl3(j,jl)     
      ENDDO
      ENDDO
      WRITE(4,*)worm_weight

      CLOSE(4)     

      PRINT*,"Configuration is written"
      
      END SUBROUTINE conf_write
!....................................................................

!--------------------------------------------------------------------
! Reading configuration
!--------------------------------------------------------------------
      SUBROUTINE CONF_READ(ping)
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: ping
      INTEGER*4,DIMENSION(1:522) :: isave
      INTEGER :: i, name, j, jl
      CHARACTER*7 :: fil

      IF(ping==0)THEN
        fil='cnf.dat'
      ELSE IF(ping==1)THEN
        fil='cnf.da1'
      ELSE
        STOP'NoNo CF'
      ENDIF

      OPEN(UNIT=4,FILE=fil)

      PRINT*,' READING NOW: ',fil

      READ(UNIT=4,FMT=*,ERR=2,END=2)nmnm
      READ(UNIT=4,FMT=*,ERR=2,END=2)MEASV,vert1,vert2
      READ(UNIT=4,FMT=*,ERR=2,END=2)present,sveta,tonya,omegaST(1:3)
      READ(UNIT=4,FMT=*,ERR=2,END=2)Energy_a,Dphase_a
! Writing vertexes para,eters
      DO i=1,nmnm; name=nlist(i)
        READ(UNIT=4,FMT=*,ERR=2,END=2)name  
        READ(UNIT=4,FMT=*,ERR=2,END=2) link(name,1:3),type(name,1:3), tau(name),site(name)
        READ(UNIT=4,FMT=*,ERR=2,END=2) omega(name,1,1:3),omega(name,2,1:3), omega(name,3,1:3)   
      ENDDO    
      READ(UNIT=4,FMT=*,ERR=2,END=2)nlist(1:N_max),numnam(1:N_max)    
! Writing hash table
      DO j=-2,N_max
      	READ(UNIT=4,FMT=*,ERR=2,END=2) linkh2(j,1),linkh2(j,2) !hash links
      ENDDO
      DO j=1,N_max
      	READ(UNIT=4,FMT=*,ERR=2,END=2) linkh3(j,1),linkh3(j,2) ! hash links
      ENDDO
      DO j=1,Nhash                                               ! hash lists 
       READ(UNIT=4,FMT=*,ERR=2,END=2) ha2(j), ha3(j)
       DO jl=1,Nhlist
      	  READ(UNIT=4,FMT=*,ERR=2,END=2) hl2(j,jl), hl3(j,jl)     
       ENDDO
      ENDDO

      CLOSE(4)     

      GOTO 1
 2    ping=ping+1
 1    CONTINUE

      PRINT*,"Configuration is read"
      
      END SUBROUTINE conf_read
!....................................................................
      
!--------------------------------------------------------------------
! Writing RAND and worm_weight
!--------------------------------------------------------------------
      SUBROUTINE RAND_WRITE(fiku)
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: fiku
      INTEGER*4,DIMENSION(1:522) :: isave
      INTEGER :: i, j, jl, name
      CHARACTER*7 :: fil

      PRINT*," *********  WRIGHTING RNDM and worm_weight ********** "
! Defining to which file to write
      IF(fiku==0)     THEN; fil='rnd.dat'
      ELSE IF(fiku==1)THEN; fil='rnd.da1'
      ELSE;           STOP'Where to write statistics'
      ENDIF 

      PRINT*,"WRITING NOW: ", fil
      
! Wrighting the counters and rndm-status
      CALL RANDSAVE(isave) 
      
      OPEN(UNIT=4,FILE=fil)
        ! Writing the rndm status
        DO i=1,522,1; WRITE(4,*) isave(i); ENDDO
        WRITE(4,*)worm_weight
      CLOSE(4)     

      PRINT*,"RNDM and worm_weight is written"
      
      END SUBROUTINE rand_write
!....................................................................
       
!--------------------------------------------------------------------
! Reading RAND and worm_weight
!--------------------------------------------------------------------
      SUBROUTINE RAND_READ(ping)
      USE config_par; USE stat_par; USE rara
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: ping
      INTEGER*4,DIMENSION(1:522) :: isave
      INTEGER :: i, name, j, jl
      CHARACTER*7 :: fil

      IF(ping==0)THEN;
        fil='rnd.dat'
      ELSE IF(ping==1)THEN;
        fil='rnd.da1'
      ELSE
        STOP'NoNo CF'
      ENDIF

      OPEN(UNIT=4,FILE=fil);         PRINT*,' READING NOW: ',fil
        ! Reading the rndm status
        DO i=1,522,1; 
          READ(UNIT=4,FMT=*,ERR=2,END=2) isave(i); 
        ENDDO
        READ(UNIT=4,FMT=*,ERR=2,END=2)worm_weight
      CLOSE(4)  
      
! restoring the rndm status
      CALL RANDLOAD(isave); PRINT*,'RNDM is RELOADED'

      GOTO 1
 2    ping=ping+1
 1    CONTINUE

      PRINT*,"RAND and worm_weight is read"
      
      END SUBROUTINE rand_read
!....................................................................
     

!--------------------------------------------------------------------
! Shifting time of one of vortices
!--------------------------------------------------------------------
      SUBROUTINE TAU_SHIFT
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN, PHONON, PHONON_MEAS, RNDM, EX
      LOGICAL,EXTERNAL :: EVA_RAT, IS_P_MEASURING
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, name3
      INTEGER :: site0, site1, site2, site3
      INTEGER :: spin_in, spin_out
      REAL*8 :: time0_old, time0_new, time1,time2,time3
      REAL*8 :: weight_old, weight_new, ratio
      REAL*8 :: sign0

      km=1; add_tau_shift=add_tau_shift+un1

! Selecting the wortex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      site0=site(name0);
! Suggesting new time
      time0_old=tau(name0); time0_new=beta*RNDM(kk)
! Finding associated vortexes
      name1=link(name0,1); name2=link(name0,2); name3=link(name0,3)
! Finding times of associated vortexes
      time1=tau(name1);  time2=tau(name2);  time3=tau(name3);
! Finding sites of associated vortexes
      site1=site(name1); site2=site(name2); site3=site(name3);
! Finding spins of in/out-going propagators
      spin_in=type(name0,1); spin_out=type(name0,2)
! Finding weights
      IF(IS_P_MEASURING(name0,name3))THEN
         weight_old = PHONON_MEAS(site0,time0_old,site3,time3)
         weight_new = PHONON_MEAS(site0,time0_new,site3,time3)
      ELSE
         weight_old = PHONON(site0,time0_old,site3,time3)
         weight_new = PHONON(site0,time0_new,site3,time3)
      ENDIF    
      IF(name1/=name0)THEN ! Not Hartree name1<-->name0
        weight_old = weight_old * GREEN(site1,time1,site0,time0_old,spin_in)  
        weight_new = weight_new * GREEN(site1,time1,site0,time0_new,spin_in)  
      ENDIF
      IF(name0/=name2)THEN ! Not Hartree name0<-->name2
        weight_old = weight_old * GREEN(site0,time0_old,site2,time2,spin_out)
        weight_new = weight_new * GREEN(site0,time0_new,site2,time2,spin_out)
      ENDIF

      ratio = weight_new / weight_old
      IF(ratio>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS(ratio)
! Playing metropolis
      a00_tau_shift=a00_tau_shift+un1
      propell=EVA_RAT(kk,ratio)
! Updating
      IF(propell)THEN; suc_tau_shift=suc_tau_shift+un1
         tau(name0)=time0_new
         Energy_a=Energy_a*ratio
         Dphase_a=Dphase_a*sign0
      ENDIF

      END SUBROUTINE TAU_SHIFT
!....................................................................

!--------------------------------------------------------------------
! Shifting site of one of vortices
!--------------------------------------------------------------------
      SUBROUTINE SITE_SHIFT
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN, PHONON, PHONON_MEAS, RNDM, EX
      LOGICAL,EXTERNAL :: EVA_RAT, IS_P_MEASURING
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, name3
      INTEGER :: site0_old, site0_new, site1, site2, site3
      INTEGER :: spin_in, spin_out
      REAL*8 :: time0, time1, time2, time3
      REAL*8 :: weight_old, weight_new, ratio
      REAL*8 :: sign0

      km=2; add_site_shift=add_site_shift+un1

! Selecting the wortex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      site0_old=site(name0);  time0=tau(name0);
! Suggesting new site
      site0_new=Nsite*RNDM(kk)+1; IF(site0_new>Nsite)site0_new=Nsite;
      IF(site0_new==site0_old)RETURN; !No update will occur
! Finding associated vortexes
      name1=link(name0,1); name2=link(name0,2); name3=link(name0,3)
! Finding times of associated vortexes
      time1=tau(name1);  time2=tau(name2);  time3=tau(name3);
! Finding sites of associated vortexes
      site1=site(name1); site2=site(name2); site3=site(name3);
! Finding spins of in/out-going propagators
      spin_in=type(name0,1); spin_out=type(name0,2)
! Finding weights
      IF( IS_P_MEASURING(name3,name0) )THEN
         weight_old = PHONON_MEAS(site0_old,time0,site3,time3)
         weight_new = PHONON_MEAS(site0_new,time0,site3,time3)
      ELSE
         weight_old = PHONON(site0_old,time0,site3,time3)
         weight_new = PHONON(site0_new,time0,site3,time3)
      ENDIF
      IF(name1/=name0)THEN ! Not Hartree
        weight_old = weight_old * GREEN(site1,time1,site0_old,time0,spin_in) 
        weight_new = weight_new * GREEN(site1,time1,site0_new,time0,spin_in)  
      ENDIF
      IF(name2/=name0)THEN ! Not Hartree
        weight_old = weight_old   * GREEN(site0_old,time0,site2,time2,spin_out)
        weight_new = weight_new  * GREEN(site0_new,time0,site2,time2,spin_out)
      ENDIF

      ratio = weight_new / weight_old
      IF(ratio>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS(ratio)
! Playing metropolis
      a00_site_shift=a00_site_shift+un1
      propell=EVA_RAT(kk,ratio)
! Updating
      IF(propell)THEN; suc_site_shift=suc_site_shift+un1
!          PRINT*,'success'
!         PRINT*,'GHGH',weight_old,weight_new
!         PRINT*,'SItes',site0_old,site0_new,site1,site2,site3
!         PRINT*,'times:',time0,time1,time2,time3
!         PRINT*,'names',name0,name3,vert1,vert2
         site(name0)=site0_new
         Energy_a=Energy_a*ratio
         Dphase_a=Dphase_a*sign0
      ENDIF

      END SUBROUTINE SITE_SHIFT
!....................................................................

!--------------------------------------------------------------------
! Flipping spins in a loop
!--------------------------------------------------------------------
      SUBROUTINE FLIP_SPINS
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name_now

      km=3; add_flip_spi=add_flip_spi+un1

      IF(POLARIZED)RETURN
      
! Selecting the wortex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
! Flipping a loop
      name_now=name0
      DO;
          type(name_now,2)=-type(name_now,2)
          type(name_now,1)=-type(name_now,1)
          name_now=link(name_now,2)
          IF(name_now==name0)EXIT
      ENDDO

      a00_flip_spi=a00_flip_spi+un1
      suc_flip_spi=suc_flip_spi+un1

      END SUBROUTINE FLIP_SPINS
!....................................................................

!--------------------------------------------------------------------
! Dressin vortex by phonon propagator
!--------------------------------------------------------------------
      SUBROUTINE ADD_PHONON
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,MODA
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio, x
      REAL*8 :: sign0

      km=4; add_add_pho=add_add_pho+un1

!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  to
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                         .                .
!                          .             .
!                            .        .
!                              . . .
!

! Limiting diagramm order
      IF(nmnm>=N_max_now)RETURN
! Selecting the vortex to dress
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      name1=link(name0,1); name2=link(name0,2)
! Avoiding dressing Hartree
      IF(name0==name1)RETURN

! Suggesting times of new vortices
      time_new1=beta*RNDM(kk); time_new2=beta*RNDM(kk);
! Suggesting sites of new vortices
      site_new1=Nsite*RNDM(kk)+1; IF(site_new1>Nsite)site_new1=Nsite;
      site_new2=Nsite*RNDM(kk)+1; IF(site_new2>Nsite)site_new2=Nsite;

! Defining parameters of name1 and name2
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta**2); pro=pro/(Nsite**2);
      pro=pro*(nmnm+un2)/(un1*nmnm)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)

! Defining gratio without context
      ratio_0 = weight_new / weight_old

! Checking whether one of particle propagators connected to name0 are measuring
      caseca=0
      IF(.NOT. MEASV)THEN;
        IF(name1==vert1)caseca=1
        IF(name0==vert1)caseca=2
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT

! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS (ratio_0 / pro)
! Adding flat order distribution factors
      ratio = ratio * (ord_probab(nmnm+2) / ord_probab(nmnm))
! Playing metropolis
      a00_add_pho=a00_add_pho+un1
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; suc_add_pho=suc_add_pho+un1
        CALL GETNAME(new_name1); CALL GETNAME(new_name2)
        ! set tau and site
        tau(new_name1)=time_new1; tau(new_name2)=time_new2
        site(new_name1)=site_new1; site(new_name2)=site_new2;
        ! set types
        type(new_name1,1)=type1; type(new_name1,2)=type1;
        type(new_name1,3)=1
        type(new_name2,1)=type2; type(new_name2,2)=type2;
        type(new_name2,3)=1
        ! set links
        link(name1,2)=new_name1; link(name0,1)=new_name1;
        link(name0,2)=new_name2; link(name2,1)=new_name2;
        link(new_name1,1)=name1; link(new_name1,2)=name0;
        link(new_name1,3)=new_name2;
        link(new_name2,1)=name0; link(new_name2,2)=name2;
        link(new_name2,3)=new_name1
        ! set momenta
        DO i=1,3;  omega(new_name1,3,i)=RNDM(kk); ENDDO
        omega(new_name2,3,1:3)=-omega(new_name1,3,1:3);
        DO i=1,3;
          omega(new_name1,2,i)= MODA(omega(name1,2,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,1,1:3)=omega(new_name1,2,1:3)
        DO i=1,3
          omega(new_name2,1,i)= MODA(omega(name2,1,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,2,1:3)=omega(new_name2,1,1:3)
        omega(new_name1,1,1:3)=omega(name1,2,1:3)
        omega(new_name2,2,1:3)=omega(name2,1,1:3)
        ! putting momenta into hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL ADDHASH(omega(new_name1,2,1),2,new_name1)
        CALL ADDHASH(omega(new_name2,2,1),2,new_name2)
        CALL ADDHASH(DABS(omega(new_name1,3,1)),3,new_name1)
        CALL ADDHASH(DABS(omega(new_name2,3,1)),3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE
        CASE(1); x=RNDM(kk);
           IF(x>half)THEN; vert1=name1; vert2=new_name1;
           ELSE;           vert1=new_name1; vert2=name0; ENDIF
        CASE(2); x=RNDM(kk)
           IF(x>half)THEN; vert1=name0; vert2=new_name2;
           ELSE;           vert1=new_name2; vert2=name2; ENDIF
        END SELECT
        ! Checking weight change
        Energy_a=Energy_a * ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
!        IF(site_new1/=site_new2)THEN
!            PRINT*,site1,site0,site2
!            PRINT*,site_new1,site_new2
!        ENDIF    
      ENDIF

      END SUBROUTINE ADD_PHONON
!....................................................................

!--------------------------------------------------------------------
! Undressing
!--------------------------------------------------------------------
      SUBROUTINE REM_PHONON
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio
      REAL*8 :: sign0

      km=5; add_rem_pho=add_rem_pho+un1
!....................................................................
!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  from
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                          .               .
!                            .           .
!                               .      .
!                                . . .


! Qualifuing for undressing the vertes
      IF(nmnm==2) RETURN            ! no phonons
      ! Select vertex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      new_name1=link(name0,1); IF(new_name1==name0)RETURN; !It is Hartree
      IF(type(new_name1,3)==0)RETURN  !Worm atached to vothex "new_name1"
      IF(MEASV .AND. (new_name1==vert1.OR.new_name1==vert2))RETURN !it is measuring line
      new_name2=link(name0,2); IF(new_name2==name0)RETURN; !Hartree (redundant)
      IF(type(new_name2,3)==0)RETURN  !worm at "new_name2" (redundant)
      IF(link(new_name1,3)/=new_name2)RETURN ! not a dressed vertex
      name2=link(new_name2,2); IF(name2==new_name1)STOP'FockReducible!'
      name1=link(new_name1,1); IF(name1==new_name2)STOP'FockRedu..!'!(redund..)
! If qualified -> continue

! Defining parameters of bames
      time_new1=tau(new_name1); time_new2=tau(new_name2)
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site_new1=site(new_name1);  site_new2=site(new_name2)
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta**2); pro=pro/(Nsite**2);
      pro=pro*(un1*nmnm-un2+un2)/(un1*nmnm-un2)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)
! Definin ratio without context
      ratio_0 =  weight_new / weight_old

! Checking if one of propagators in block are measuring
      caseca=0
      IF(.NOT. MEASV)THEN
        IF(name1==vert1)caseca=1; IF(new_name1==vert1)caseca=2
        IF(name0==vert1)caseca=3; IF(new_name2==vert1)caseca=4
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(3); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(4); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT    
      
! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS( ratio_0 / pro )
! Adding flat orders distribution factros
      ratio = ratio * (ord_probab(nmnm-2+2) / ord_probab(nmnm-2))
! Reverting ratio      
      ratio = un1 / ratio  ! because it is remove
! Playing metropolis
      a00_rem_pho=a00_rem_pho+un1
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; suc_rem_pho=suc_rem_pho+un1
          !PRINT*,'succedeed to delete general 8888888'
       ! set links
        link(name1,2)=name0; link(name0,1)=name1;
        link(name0,2)=name2; link(name2,1)=name0;
        ! setting momenta
        omega(name0,2,1:3)=omega(name2,1,1:3)
        omega(name0,1,1:3)=omega(name1,2,1:3)
        ! removing momenta from hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL DROPHASH(2,new_name1)
        CALL DROPHASH(2,new_name2)
        CALL DROPHASH(3,new_name1)
        CALL DROPHASH(3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE;
        CASE(1); vert2=name0;
        CASE(2); vert1=name1;
        CASE(3); vert2=name2;
        CASE(4); vert1=name0;
        END SELECT;
        CALL DROPNAME(new_name1); CALL DROPNAME(new_name2)
        ! Checking weight change
        Energy_a=Energy_a / ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
      ENDIF

      END SUBROUTINE REM_PHONON
!....................................................................

!--------------------------------------------------------------------
! Inserting worms
!--------------------------------------------------------------------
      SUBROUTINE ADD_WORMS
      USE config_par; USE stat_par
      IMPLICIT NONE
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      REAL*8,EXTERNAL :: RNDM,MODA
      REAL*8 :: pro, ratio, x
      INTEGER :: num, name1, name2
      INTEGER :: i, go_to, ways, caseca

      km=6; add_add_worm=add_add_worm+un1;

      IF(present)RETURN;

! Selecting first vortex to put worm
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name1=nlist(num)
! Selecting link to another worm
      go_to=3*RNDM(kk)+1; IF(go_to>3)go_to=3
      name2=link(name1,go_to)
! Checking if it is Hartree
      IF(name1==name2)RETURN

! Defining context factor
      ways=0; DO i=1,3; IF(link(name1,i)==name2)ways=ways+1; ENDDO;
      pro = (un1*ways) / (un1*nmnm*3)

      ratio = worm_weight / pro

! Playing metropolis
      a00_add_worm=a00_add_worm+un1;
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; present=.TRUE.; suc_add_worm=suc_add_worm+un1;
        ! Distributing sveta and tonya
        x=RNDM(kk)
        IF(x>half)THEN; sveta=name1; tonya=name2;
        ELSE;           sveta=name2; tonya=name1; ENDIF
        IF     (sveta==name1 .AND. go_to==1)THEN; caseca=2! T --> S
        ELSE IF(sveta==name2 .AND. go_to==1)THEN; caseca=1! S --> T
        ELSE IF(sveta==name1 .AND. go_to==2)THEN; caseca=1! S --> T
        ELSE IF(sveta==name2 .AND. go_to==2)THEN; caseca=2! T --> S
        ELSE IF(                   go_to==3)THEN; caseca=3! S === T
        ENDIF
        ! Setting types
        type(sveta,3)=0; type(tonya,3)=0;
        type(link(sveta,3),3)=0; type(link(tonya,3),3)=0;
        ! Redefining momenta
        DO i=1,3; omegaST(i)=RNDM(kk); ENDDO;
        SELECT CASE(caseca);
        CASE(1); ! Sveta --->--- Tonya
          DO i=1,3; omega(sveta,2,i)=MODA(omega(sveta,2,i)-omegaST(i))
          ENDDO; omega(tonya,1,1:3)=omega(sveta,2,1:3);
          CALL DROPHASH(2,sveta); CALL ADDHASH(omega(sveta,2,1),2,sveta)
        CASE(2); ! Tonya --->--- Sveta
          DO i=1,3; omega(sveta,1,i)=MODA(omega(sveta,1,i)+omegaST(i))
          ENDDO; omega(tonya,2,1:3)=omega(sveta,1,1:3);
          CALL DROPHASH(2,tonya); CALL ADDHASH(omega(tonya,2,1),2,tonya)
        CASE(3); ! Sveta ======= Tonya
          DO i=1,3; omega(sveta,3,i)=MODA(omega(sveta,3,i)-omegaST(i))
          ENDDO; omega(tonya,3,1:3)=-omega(sveta,3,1:3);
          CALL DROPHASH(3,sveta); CALL DROPHASH(3,tonya)
          CALL ADDHASH(DABS(omega(sveta,3,1)),3,sveta)
          CALL ADDHASH(DABS(omega(tonya,3,1)),3,tonya)
        END SELECT
      ENDIF

      END SUBROUTINE ADD_WORMS
!....................................................................

!--------------------------------------------------------------------
! Removing worms
!--------------------------------------------------------------------
      SUBROUTINE REM_WORMS
      USE config_par; ; USE stat_par
      IMPLICIT NONE
      LOGICAL,EXTERNAL :: EVA_RAT,REDUC
      LOGICAL :: propell
      REAL*8,EXTERNAL :: RNDM,MODA
      REAL*8,DIMENSION(3) :: onew
      REAL*8 :: pro, ratio
      INTEGER :: i, ways, idirec, ik, io, ways1

      km=7; add_rem_worm=add_rem_worm+un1;

      IF(.NOT. present)RETURN

      idirec=0
! Finding number of associating ways (ad context factor)
      ways=0; DO i=1,3; IF(link(sveta,i)==tonya)ways=ways+1; ENDDO;
      IF(ways==0)RETURN  !worms are not connected

! Selecting the specific associating way "idirec"
      ik=RNDM(kk)*ways+1; IF(ik>ways)ik=ways

      ways1=0
      DO i=1,3; IF(link(sveta,i)==tonya)ways1=ways1+1;
        IF(ik==ways1)THEN; idirec=i; EXIT; ENDIF
      ENDDO;

      IF(idirec==1)THEN;
         DO io=1,3; onew(io)=MODA(omega(sveta,idirec,io)-omegaST(io));
         ENDDO
      ELSE IF(idirec>1)THEN
         DO io=1,3; onew(io)=MODA(omega(sveta,idirec,io)+omegaST(io));
         ENDDO
      ELSE
         PRINT*,'idirec=',idirec,ways,ik; STOP'idirec';
      ENDIF

      IF(REDUC(idirec,onew))RETURN

      pro = (un1*ways) / (un1*nmnm*3)

!      PRINT*,"To check : pro=",pro

      ratio = worm_weight / pro
      ratio = un1 / ratio

! Playing metropolis
      a00_rem_worm=a00_rem_worm+un1;
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; present=.FALSE.; suc_rem_worm=suc_rem_worm+un1;
       SELECT CASE(idirec);
         CASE(1);
           CALL DROPHASH(2,tonya)
           omega(sveta,1,1:3)=onew(1:3)
           omega(tonya,2,1:3)=omega(sveta,1,1:3)
           CALL ADDHASH(omega(tonya,2,1),2,tonya);
         CASE(2);
           CALL DROPHASH(2,sveta);
           omega(sveta,2,1:3)=onew(1:3)
           omega(tonya,1,1:3)=omega(sveta,2,1:3)
           CALL ADDHASH(omega(sveta,2,1),2,sveta)
         CASE(3);
           CALL DROPHASH(3,sveta); CALL DROPHASH(3,tonya);
           omega(sveta,3,1:3)=onew(1:3)
           omega(tonya,3,1:3)=-omega(sveta,3,1:3)
           CALL ADDHASH(DABS(omega(sveta,3,1)),3,sveta);
           CALL ADDHASH(DABS(omega(tonya,3,1)),3,tonya);
      END SELECT
      type(sveta,3)=1; type(tonya,3)=1
      type(link(sveta,3),3)=1; type(link(tonya,3),3)=1
      present=.FALSE. ; sveta=0; tonya=0; omegaST=0.d0
      ENDIF

      END SUBROUTINE REM_WORMS
!....................................................................

!--------------------------------------------------------------------
! Shifting mesuring line
!--------------------------------------------------------------------
      SUBROUTINE MOVEMEAS
      USE config_par; USE stat_par; IMPLICIT NONE;
      REAL*8,EXTERNAL :: RNDM, PHONON, PHONON_MEAS
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell, measV_NEW
      INTEGER :: numa, vert1_new, vert2_new, vert1_old, vert2_old
      INTEGER :: site_1_old, site_1_new, site_2_old, site_2_new
      REAL*8 :: ukaka, old_weight, new_weight, ratio
      REAL*8 :: time_1_new, time_1_old, time_2_new, time_2_old, sign0
      REAL*8 :: context, ratio_0
  
      km=10; 
      add_move_mea=add_move_mea+un1; 

! Selects position for new beginning of measuring line      
      numa=RNDM(kk)*nmnm+1; IF(numa>nmnm)numa=nmnm;
      vert1_new=nlist(numa)
! Selecting type of new measuring line 
      ukaka=RNDM(kk)   
      IF(ukaka<0.5d0)THEN
          vert2_new=link(vert1_new,2); measV_NEW=.FALSE. !G measuring
      ELSE
          vert2_new=link(vert1_new,3); measV_NEW=.TRUE. !phonon-measuring 
      ENDIF
! Parameter for new measuring line
      time_1_new = tau(vert1_new); time_2_new = tau(vert2_new); 
      site_1_new = site(vert1_new); site_2_new = site(vert2_new);      
      
! Parameters for old measuring line      
      vert1_old = vert1;  vert2_old = vert2
      time_1_old = tau(vert1_old); time_2_old = tau(vert2_old); 
      site_1_old = site(vert1_old); site_2_old = site(vert2_old);
      
! Return if nothing changes      
      IF((vert1_new==vert1_old) .AND. (vert2_new==vert2_old))RETURN! Nothing changed

      context=1.0d0
! Forming old and new weight      
      IF(measV .AND. measV_NEW)THEN; ! both measuring lines are phonons
         old_weight = PHONON_MEAS(site_1_old,time_1_old,site_2_old,time_2_old) * &
                             PHONON(site_1_new,time_1_new,site_2_new,time_2_new)
         new_weight = PHONON(site_1_old,time_1_old,site_2_old,time_2_old) * &
                              PHONON_MEAS(site_1_new,time_1_new,site_2_new,time_2_new)
      ELSE IF(measV .AND. (.NOT. measV_NEW))THEN; ! old measuring line is phonon
         old_weight =  PHONON_MEAS(site_1_old,time_1_old,site_2_old,time_2_old)
         new_weight = PHONON(site_1_old,time_1_old,site_2_old,time_2_old)
         !context=2.0d0
      ELSE IF((.NOT. measV) .AND. measV_NEW)THEN; ! new measuring line is phonon
         old_weight = PHONON(site_1_new,time_1_new,site_2_new,time_2_new)
         new_weight = PHONON_MEAS(site_1_new,time_1_new,site_2_new,time_2_new)
         !context=0.5
      ELSE IF((.NOT. measV) .AND. (.NOT. measV_NEW))THEN; ! both are not phonons
          old_weight=1.0d0; new_weight=1.0d0
      ELSE
          STOP"Where is measuring line"
      ENDIF    
      
! Estimating ratio      
      ratio_0 = new_weight / old_weight
      ratio=ABS(ratio_0/context)
! Playing Metropolis         
      a00_move_mea=a00_move_mea+un1
      propell=EVA_RAT(kk,ratio)
      
! Updating if success
      IF(propell)THEN
          suc_move_mea=suc_move_mea+un1
          vert1=vert1_new; vert2=vert2_new; measV=measV_NEW
          IF(ratio>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
          ! Checking weight change
          Energy_a=Energy_a * ABS(ratio_0)
          Dphase_a=Dphase_a*sign0
      ENDIF 
      
      END SUBROUTINE MOVEMEAS
!....................................................................

!--------------------------------------------------------------------
! Moving worms
!--------------------------------------------------------------------
      SUBROUTINE MOVE_WORMS
      USE config_par; ; USE stat_par
      IMPLICIT NONE
      LOGICAL,EXTERNAL :: REDUC
      REAL*8,EXTERNAL :: RNDM,MODA
      REAL*8,DIMENSION(3) :: onew
      LOGICAL :: sve
      REAL*8 :: r
      INTEGER :: i, io, name0, name1, name2, name3, caseca
      INTEGER :: name00

      km=8; add_move_worm=add_move_worm+un1
            a00_move_worm=a00_move_worm+un1


      IF(.NOT. present)RETURN

! Chosing whom to move
      r=RNDM(kk)
      IF(r<half)THEN; name0=sveta; name00=tonya; sve=.TRUE. ! moving sveta
      ELSE;           name0=tonya; name00=sveta; sve=.FALSE.! moving tonya
      ENDIF
      name1=link(name0,1); name2=link(name0,2); name3=link(name0,3)

! Chosing direction along which to move the worm
      r=RNDM(kk);
      IF(r<half)THEN
        r=RNDM(kk); IF(r<half)THEN; caseca=1; ! along incoming
                    ELSE;           caseca=2; ! along oucoming
                    ENDIF;
      ELSE;                         caseca=3; ! along interaction
      ENDIF

! Accepting update
      SELECT CASE(caseca);
        CASE(1); IF(name1==name00)RETURN;; !adjacent
          IF(sve)THEN; ! Sveta along incoming
            DO io=1,3; onew(io)=MODA(omega(name1,2,io)-omegaST(io));
            ENDDO
            IF(REDUC(1,onew))RETURN; suc_move_worm=suc_move_worm+un1
            omega(name1,2,1:3)=onew(1:3)
            omega(sveta,1,1:3)=onew(1:3);
            CALL DROPHASH(2,name1);
            CALL ADDHASH(omega(name1,2,1),2,name1)
            IF(name00/=link(name0,3))THEN;
              type(name0,3)=1; type(link(name0,3),3)=1;
            ENDIF
            type(name1,3)=0; type(link(name1,3),3)=0
            sveta=name1;
          ELSE
            DO io=1,3; onew(io)=MODA(omega(name1,2,io)+omegaST(io));
            ENDDO
            IF(REDUC(1,onew))RETURN; suc_move_worm=suc_move_worm+un1
            omega(name1,2,1:3)=onew(1:3)
            omega(tonya,1,1:3)=onew(1:3);
            CALL DROPHASH(2,name1);
            CALL ADDHASH(omega(name1,2,1),2,name1)
            IF(name00.ne.link(name0,3))THEN;
              type(name0,3)=1; type(link(name0,3),3)=1;
            ENDIF
            type(name1,3)=0; type(link(name1,3),3)=0
            tonya=name1;
          ENDIF
        CASE(2); IF(name2==name00)RETURN;; !adjacent
          IF(sve)THEN;
            DO io=1,3; onew(io)=MODA(omega(sveta,2,io)+omegaST(io));
            ENDDO
            IF(REDUC(1,onew))RETURN; suc_move_worm=suc_move_worm+un1
            omega(sveta,2,1:3)=onew(1:3);
            omega(name2,1,1:3)=onew(1:3);
            CALL DROPHASH(2,sveta);
            CALL ADDHASH(omega(sveta,2,1),2,sveta)
            IF(name00.ne.link(name0,3))THEN;
               type(name0,3)=1; type(link(name0,3),3)=1;
            ENDIF
            type(name2,3)=0; type(link(name2,3),3)=0
            sveta=name2;
          ELSE
            DO io=1,3; onew(io)=MODA(omega(tonya,2,io)-omegaST(io));
            ENDDO
            IF(REDUC(1,onew))RETURN; suc_move_worm=suc_move_worm+un1
            omega(tonya,2,1:3)=onew(1:3)
            omega(name2,1,1:3)=onew(1:3);
            CALL DROPHASH(2,tonya);
            CALL ADDHASH(omega(tonya,2,1),2,tonya)
            IF(name00.ne.link(name0,3))THEN;
              type(name0,3)=1; type(link(name0,3),3)=1;
            ENDIF
            type(name2,3)=0; type(link(name2,3),3)=0
            tonya=name2;
          ENDIF
        CASE(3); IF(name3==name00)RETURN;; !adjacent
          IF(sve)THEN;
            DO io=1,3; onew(io)=MODA(omega(sveta,3,io)+omegaST(io));
            ENDDO
            IF(REDUC(3,onew))RETURN; suc_move_worm=suc_move_worm+un1
            omega(sveta,3,1:3)=onew(1:3)
            omega(name3,3,1:3)=-onew(1:3)
            CALL DROPHASH(3,sveta); CALL DROPHASH(3,name3)
            CALL ADDHASH(DABS(omega(sveta,3,1)),3,sveta)
            CALL ADDHASH(DABS(omega(name3,3,1)),3,name3)
            sveta=name3;
          ELSE;
            DO io=1,3; onew(io)=MODA(omega(tonya,3,io)-omegaST(io));
            ENDDO
            IF(REDUC(3,onew))RETURN; suc_move_worm=suc_move_worm+un1
            omega(tonya,3,1:3)=onew(1:3)
            omega(name3,3,1:3)=-onew(1:3);
            CALL DROPHASH(3,tonya); CALL DROPHASH(3,name3)
            CALL ADDHASH(DABS(omega(tonya,3,1)),3,tonya)
            CALL ADDHASH(DABS(omega(name3,3,1)),3,name3)
            tonya=name3;
          ENDIF
      END SELECT

      END SUBROUTINE MOVE_WORMS
!....................................................................

!--------------------------------------------------------------------
! Commuting Green functoins outcoming from worm
!--------------------------------------------------------------------
      SUBROUTINE COMMUTE
      USE config_par; ; USE stat_par
      IMPLICIT NONE
      LOGICAL,EXTERNAL :: EVA_RAT
      REAL*8,EXTERNAL :: RNDM,GREEN,PHONON,MODA
      LOGICAL :: propell
      REAL*8,DIMENSION(3) :: omegaST_new
      REAL*8 :: tau_SV, tau_to_SV,tau_TO, tau_to_TO
      REAL*8 :: weight_old, weight_new, ratio
      INTEGER :: name_SV, name_TO, spin_SV, spin_TO, fault_count
      INTEGER :: site_SV, site_to_SV, site_TO, site_to_TO, io
      REAL*8 :: sign0

      km=9; add_commute=add_commute+un1

! Looking if to skip
      IF(.NOT. present)RETURN
      spin_SV=type(sveta,2); spin_TO=type(tonya,2)
      IF(spin_TO /= spin_SV)RETURN

      name_SV=link(sveta,2); name_TO=link(tonya,2)
      tau_SV  = tau(sveta);   tau_to_SV  = tau(name_SV);
      tau_TO  = tau(tonya);   tau_to_TO  = tau(name_TO);
      site_SV = site(sveta);  site_to_SV = site(name_SV)
      site_TO = site(tonya);  site_to_TO = site(name_TO)

      weight_old = GREEN(site_SV,tau_SV,site_to_SV,tau_to_SV,spin_SV)* &
                          GREEN(site_TO,tau_TO,site_to_TO,tau_to_TO,spin_SV)
      weight_new = GREEN(site_SV,tau_SV,site_to_TO,tau_to_TO,spin_SV)* &
                            GREEN(site_TO,tau_TO,site_to_SV,tau_to_SV,spin_SV)


      ratio = weight_new / weight_old
      IF(ratio>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS(ratio)

      DO io=1,3; 
          omegaST_new(io)= MODA(omegaST(io)+omega(sveta,2,io)-omega(tonya,2,io)); 
      ENDDO
! Does not allow Hartree-Hartree loops disconnected structures
        fault_count=0; DO io=1,3
        IF(ABS(omegaST_new(io))<1.0d-10)fault_count=fault_count+1; ENDDO
        IF(fault_count==3)RETURN

      a00_commute=a00_commute+un1
      propell=EVA_RAT(kk,ratio)

      IF(propell)THEN; suc_commute=suc_commute+un1
        omegaST(1:3)=omegaST_new(1:3)
        omega(sveta,2,1:3)=omega(name_TO,1,1:3)
        omega(tonya,2,1:3)=omega(name_SV,1,1:3)
        link(sveta,2)=name_TO; link(tonya,2)=name_SV;
        link(name_SV,1)=tonya; link(name_TO,1)=sveta;
!%%%%% update hash
        CALL DROPHASH(2,sveta); CALL DROPHASH(2,tonya)
        CALL ADDHASH(omega(sveta,2,1),2,sveta)
        CALL ADDHASH(omega(tonya,2,1),2,tonya)
!%%% re-assign the measuring propagator line if vert1 is worm
        IF(.not.measV .AND. sveta==vert1)vert2=name_TO
        IF(.not.measV .AND. tonya==vert1)vert2=name_SV
! Checking weight change
        Energy_a = Energy_a * ratio
        Dphase_a = - Dphase_a*sign0
      ENDIF

        END SUBROUTINE COMMUTE
!....................................................................

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EFFICIENCY: Additional updates to boost efficiency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

!--------------------------------------------------------------------
! Commuting interaction lines outcoming from worm
!--------------------------------------------------------------------
      SUBROUTINE COMMUTE_INT
      USE config_par; ; USE stat_par
      IMPLICIT NONE
      LOGICAL,EXTERNAL :: EVA_RAT
      REAL*8,EXTERNAL :: RNDM, GREEN, PHONON, PHONON_MEAS, MODA
      LOGICAL,EXTERNAL :: IS_P_MEASURING
      LOGICAL :: propell
      REAL*8,DIMENSION(3) :: omegaST_new
      REAL*8 :: tau_SV, tau_to_SV,tau_TO, tau_to_TO
      REAL*8 :: weight_old, weight_new, ratio
      INTEGER :: name_SV, name_TO, fault_count
      INTEGER :: site_SV, site_to_SV, site_TO, site_to_TO, io
      REAL*8 :: sign0

      km=13; add_commute_int=add_commute_int+un1

! Skipping if no worm
      IF(.NOT. present)RETURN
! Find links of Sveta and Tonya through phonons 
      name_SV=link(sveta,3); name_TO=link(tonya,3)
! Avoiding equivalent result
      IF(name_SV==tonya)RETURN
! Setting parameters of vortexes      
      tau_SV  = tau(sveta);   tau_to_SV  = tau(name_SV);
      tau_TO  = tau(tonya);   tau_to_TO  = tau(name_TO);
      site_SV = site(sveta);  site_to_SV = site(name_SV)
      site_TO = site(tonya);  site_to_TO = site(name_TO)
! Initilizing weights
      weight_old=un1; weight_new=1.0d0
! Weghts for switching Sveta's link to Tonya's link      
      IF( IS_P_MEASURING(sveta,name_SV) )THEN
           weight_old = weight_old * PHONON_MEAS(site_SV,tau_SV,site_to_SV,tau_to_SV)
           weight_new = weight_new * PHONON_MEAS(site_SV,tau_SV,site_to_TO,tau_to_TO)
      ELSE
           weight_old = weight_old * PHONON(site_SV,tau_SV,site_to_SV,tau_to_SV)
           weight_new = weight_new * PHONON(site_SV,tau_SV,site_to_TO,tau_to_TO)
      ENDIF          
! Weghts for switching Tonya's link to Sveta's link      
      IF( IS_P_MEASURING(tonya,name_TO) )THEN
            weight_old = weight_old * PHONON_MEAS(site_TO,tau_TO,site_to_TO,tau_to_TO) 
            weight_new = weight_new * PHONON_MEAS(site_TO,tau_TO,site_to_SV,tau_to_SV)
       ELSE     
            weight_old = weight_old * PHONON(site_TO,tau_TO,site_to_TO,tau_to_TO) 
            weight_new = weight_new * PHONON(site_TO,tau_TO,site_to_SV,tau_to_SV)
       ENDIF

! Finding ratio 
      ratio = weight_new / weight_old
      IF(ratio>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS(ratio)

      DO io=1,3; 
          omegaST_new(io)=MODA(omegaST(io)+omega(sveta,3,io)-omega(tonya,3,io)); 
      ENDDO
! Does not allow Hartree-Hartree loops disconnected structures
        fault_count=0; DO io=1,3
        IF(ABS(omegaST_new(io))<1.0d-10)fault_count=fault_count+1; ENDDO
        IF(fault_count==3)RETURN

! Playing Metropolis        
      a00_commute_int=a00_commute_int+un1
      propell=EVA_RAT(kk,ratio)
! Updating if success
      IF(propell)THEN; suc_commute_int=suc_commute_int+un1
        omegaST(1:3)=omegaST_new(1:3)
        omega(sveta,3,1:3)=-omega(name_TO,3,1:3)
        omega(tonya,3,1:3)=-omega(name_SV,3,1:3)
        link(sveta,3)=name_TO; link(tonya,3)=name_SV;
        link(name_SV,3)=tonya; link(name_TO,3)=sveta;
!%%%%% update hash
        CALL DROPHASH(3,sveta); CALL DROPHASH(3,tonya)
        CALL ADDHASH(DABS(omega(sveta,3,1)),3,sveta)
        CALL ADDHASH(DABS(omega(tonya,3,1)),3,tonya)
!%%% re-assign the measuring propagator line if vert1 is worm
        IF(measV .AND. sveta==vert1)vert2=name_TO
        IF(measV .AND. tonya==vert1)vert2=name_SV
        IF(measV .AND. sveta==vert2)vert1=name_TO
        IF(measV .AND. tonya==vert2)vert1=name_SV
! Checking weight change
        Energy_a = Energy_a * ratio
        Dphase_a = Dphase_a*sign0
      ENDIF

        END SUBROUTINE COMMUTE_INT
!....................................................................

!--------------------------------------------------------------------
! Dressin vortex by diagonal phonon propagator with exponential tau probability
!--------------------------------------------------------------------
      SUBROUTINE ADD_PHONON_DIA
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,MODA
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio, x
      REAL*8 :: sign0

      km=11; add_add_dia=add_add_dia+un1

!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  to
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                         .                .
!                           .           .
!                            .      .
!                             . . .
!

! Limiting diagramm order
      IF(nmnm>=N_max_now)RETURN
! Selecting the vortex to dress
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      name1=link(name0,1); name2=link(name0,2)
! Avoiding dressing Hartree
      IF(name0==name1)RETURN

! Suggesting times of new vortices
      time_new1=beta*RNDM(kk); time_new2=beta*RNDM(kk);
! Suggesting sites of new vortices
      site_new1=Nsite*RNDM(kk)+1; IF(site_new1>Nsite)site_new1=Nsite;
      site_new2=site_new1

! Defining parameters of name1 and name2
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta**2); pro=pro/(Nsite);
      pro=pro*(nmnm+un2)/(un1*nmnm)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)

! Definin gratio without context
      ratio_0 = weight_new / weight_old

! Checking whether one of propagators connected to name0 are measuring
      caseca=0
      IF(.NOT. MEASV)THEN;
        IF(name1==vert1)caseca=1
        IF(name0==vert1)caseca=2
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT

! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS (ratio_0 / pro)
! Adding flat orders distribution factors
      ratio = ratio * (ord_probab(nmnm+2) / ord_probab(nmnm))
! Playing metropolis
      a00_add_dia=a00_add_dia+un1
      propell=EVA_RAT(kk,ratio)
      
! Updating
      IF(propell)THEN; suc_add_dia=suc_add_dia+un1; 
        CALL GETNAME(new_name1); CALL GETNAME(new_name2)
        ! set tau and site
        tau(new_name1)=time_new1; tau(new_name2)=time_new2
        site(new_name1)=site_new1; site(new_name2)=site_new2;
        ! set types
        type(new_name1,1)=type1; type(new_name1,2)=type1;
        type(new_name1,3)=1
        type(new_name2,1)=type2; type(new_name2,2)=type2;
        type(new_name2,3)=1
        ! set links
        link(name1,2)=new_name1; link(name0,1)=new_name1;
        link(name0,2)=new_name2; link(name2,1)=new_name2;
        link(new_name1,1)=name1; link(new_name1,2)=name0;
        link(new_name1,3)=new_name2;
        link(new_name2,1)=name0; link(new_name2,2)=name2;
        link(new_name2,3)=new_name1
        ! set momenta
        DO i=1,3;  omega(new_name1,3,i)=RNDM(kk); ENDDO
        omega(new_name2,3,1:3)=-omega(new_name1,3,1:3);
        DO i=1,3;
          omega(new_name1,2,i)=MODA(omega(name1,2,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,1,1:3)=omega(new_name1,2,1:3)
        DO i=1,3
          omega(new_name2,1,i)=MODA(omega(name2,1,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,2,1:3)=omega(new_name2,1,1:3)
        omega(new_name1,1,1:3)=omega(name1,2,1:3)
        omega(new_name2,2,1:3)=omega(name2,1,1:3)
        ! putting momenta into hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL ADDHASH(omega(new_name1,2,1),2,new_name1)
        CALL ADDHASH(omega(new_name2,2,1),2,new_name2)
        CALL ADDHASH(DABS(omega(new_name1,3,1)),3,new_name1)
        CALL ADDHASH(DABS(omega(new_name2,3,1)),3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE
        CASE(1); x=RNDM(kk);
           IF(x>half)THEN; vert1=name1; vert2=new_name1;
           ELSE;           vert1=new_name1; vert2=name0; ENDIF
        CASE(2); x=RNDM(kk)
           IF(x>half)THEN; vert1=name0; vert2=new_name2;
           ELSE;           vert1=new_name2; vert2=name2; ENDIF
        END SELECT
        ! Checking weight change
        Energy_a=Energy_a * ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
!        IF(site_new1/=site_new2)THEN
!            PRINT*,site1,site0,site2
!            PRINT*,site_new1,site_new2
!        ENDIF    
      ENDIF

      END SUBROUTINE ADD_PHONON_DIA
!....................................................................

!--------------------------------------------------------------------
! Undressing diagonal phonon
!--------------------------------------------------------------------
      SUBROUTINE REM_PHONON_DIA
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio
      REAL*8 :: sign0

      km=12; add_rem_dia=add_rem_dia+un1
!....................................................................
!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  from
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                          .               .
!                           .           .
!                              .      .
!                               . . .


! Qualifuing for undressing the vertes
      IF(nmnm==2) RETURN            ! no phonons
      ! Select vertex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      new_name1=link(name0,1); IF(new_name1==name0)RETURN; !It is Hartree
      IF(type(new_name1,3)==0)RETURN  !Worm atached to vothex "new_name1"
      IF(MEASV .AND. (new_name1==vert1.OR.new_name1==vert2))RETURN !it is measuring line
      new_name2=link(name0,2); IF(new_name2==name0)RETURN; !Hartree (redundant)
      IF(type(new_name2,3)==0)RETURN  !worm at "new_name2" (redundant)
      IF(link(new_name1,3)/=new_name2)RETURN ! not a dressed vertex
      name2=link(new_name2,2); IF(name2==new_name1)STOP'FockReducible!'
      name1=link(new_name1,1); IF(name1==new_name2)STOP'FockRedu..!'!(redund..)
! If qualified -> continue

! Defining parameters of lines
      time_new1=tau(new_name1); time_new2=tau(new_name2)
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site_new1=site(new_name1);  site_new2=site(new_name2)
      IF(site_new1/=site_new2)RETURN
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta**2); pro=pro/(Nsite);
      pro=pro*(un1*nmnm-un2+un2)/(un1*nmnm-un2)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)
! Definin ratio without context
      ratio_0 =  weight_new / weight_old

! Checking if one of propagators in block are measuring
      caseca=0
      IF(.NOT. MEASV)THEN
        IF(name1==vert1)caseca=1; IF(new_name1==vert1)caseca=2
        IF(name0==vert1)caseca=3; IF(new_name2==vert1)caseca=4
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(3); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(4); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT    
      
! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS( ratio_0 / pro )
! Adding flat ordest distribution factors
      ratio = ratio * (ord_probab(nmnm-2+2) / ord_probab(nmnm-2))
! Reverting ration for remove      
      ratio = un1 / ratio  ! because it is remove
! Playing metropolis
      a00_rem_dia=a00_rem_dia+un1
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; suc_rem_dia=suc_rem_dia+un1
          !PRINT*,'succed to delete****************'
       ! set links
        link(name1,2)=name0; link(name0,1)=name1;
        link(name0,2)=name2; link(name2,1)=name0;
        ! setting momenta
        omega(name0,2,1:3)=omega(name2,1,1:3)
        omega(name0,1,1:3)=omega(name1,2,1:3)
        ! removing momenta from hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL DROPHASH(2,new_name1)
        CALL DROPHASH(2,new_name2)
        CALL DROPHASH(3,new_name1)
        CALL DROPHASH(3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE;
        CASE(1); vert2=name0;
        CASE(2); vert1=name1;
        CASE(3); vert2=name2;
        CASE(4); vert1=name0;
        END SELECT;
        CALL DROPNAME(new_name1); CALL DROPNAME(new_name2)
        ! Checking weight change
        Energy_a=Energy_a / ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
      ENDIF

      END SUBROUTINE REM_PHONON_DIA
!....................................................................

!--------------------------------------------------------------------
! Dressin vortex by diagonal phonon propagator
!--------------------------------------------------------------------
      SUBROUTINE ADD_PHONON_DIA_EXP
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN, PHONON, RNDM, MODA, EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio, x
      REAL*8 :: sign0
      REAL*8 :: RRd, pr1, del_time, bakara

      km=14; add_add_dia_exp=add_add_dia_exp+un1

!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  to
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                         .                .
!                           .           .
!                            .      .
!                             . . .
!
 
! Limiting diagramm order
      IF(nmnm>=N_max_now)RETURN
! Selecting the vortex to dress
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      name1=link(name0,1); name2=link(name0,2)
! Avoiding dressing Hartree
      IF(name0==name1)RETURN

! Suggesting times of new vortices
      !first suggested time is random
      time_new1=beta*RNDM(kk); 
      !second suggeted time is an exponential distribution
      bazeka = bazeka_min + (RNDM(kk)*(bazeka_max-bazeka_min))
      RRd=RNDM(kk)
      bakara = un1-(RRd*(un1-EX(-bazeka*beta)))
      IF(bakara>1.0d-40)THEN 
          del_time = - (un1/bazeka) * LOG( bakara ) 
      ELSE
          del_time = RRd * beta 
      ENDIF    
      IF(del_time>beta)THEN;
          PRINT*,"RRd = ,",RRd
          PRINT*,'Del_time = ',del_time
          STOP"Panic in del_time"
      ENDIF    
      time_new2 = time_new1 + del_time; 
      IF(time_new2>beta)time_new2=time_new2-beta
! Suggesting sites of new vortices
      site_new1=Nsite*RNDM(kk)+1; IF(site_new1>Nsite)site_new1=Nsite;
      site_new2=site_new1

! Defining parameters of name1 and name2
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta); pro=pro/(Nsite);
      pro=pro*(nmnm+un2)/(un1*nmnm)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)

! Definin gratio without context
      ratio_0 = weight_new / weight_old

! Checking whether one of propagators connected to name0 are measuring
      caseca=0
      IF(.NOT. MEASV)THEN;
        IF(name1==vert1)caseca=1
        IF(name0==vert1)caseca=2
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT

! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS (ratio_0 / pro)
! Adding flat orders distribution factors
      ratio = ratio * (ord_probab(nmnm+2) / ord_probab(nmnm))
! Correcting ratio by exponential distribution time weight function
      ! exponential distribution
      pr1 = bazeka * EX(-bazeka*del_time) / (un1-EX(-bazeka*beta)) 
      ratio = ratio / pr1; !PRINT*,"ADD: ",ratio
! Playing metropolis
      a00_add_dia_exp=a00_add_dia_exp+un1
      propell=EVA_RAT(kk,ratio)
      
! Updating
      IF(propell)THEN; suc_add_dia_exp=suc_add_dia_exp+un1; 
        CALL GETNAME(new_name1); CALL GETNAME(new_name2)
        ! set tau and site
        tau(new_name1)=time_new1; tau(new_name2)=time_new2
        site(new_name1)=site_new1; site(new_name2)=site_new2;
        ! set types
        type(new_name1,1)=type1; type(new_name1,2)=type1;
        type(new_name1,3)=1
        type(new_name2,1)=type2; type(new_name2,2)=type2;
        type(new_name2,3)=1
        ! set links
        link(name1,2)=new_name1; link(name0,1)=new_name1;
        link(name0,2)=new_name2; link(name2,1)=new_name2;
        link(new_name1,1)=name1; link(new_name1,2)=name0;
        link(new_name1,3)=new_name2;
        link(new_name2,1)=name0; link(new_name2,2)=name2;
        link(new_name2,3)=new_name1
        ! set momenta
        DO i=1,3;  omega(new_name1,3,i)=RNDM(kk); ENDDO
        omega(new_name2,3,1:3)=-omega(new_name1,3,1:3);
        DO i=1,3;
          omega(new_name1,2,i)=MODA(omega(name1,2,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,1,1:3)=omega(new_name1,2,1:3)
        DO i=1,3
          omega(new_name2,1,i)=MODA(omega(name2,1,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,2,1:3)=omega(new_name2,1,1:3)
        omega(new_name1,1,1:3)=omega(name1,2,1:3)
        omega(new_name2,2,1:3)=omega(name2,1,1:3)
        ! putting momenta into hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL ADDHASH(omega(new_name1,2,1),2,new_name1)
        CALL ADDHASH(omega(new_name2,2,1),2,new_name2)
        CALL ADDHASH(DABS(omega(new_name1,3,1)),3,new_name1)
        CALL ADDHASH(DABS(omega(new_name2,3,1)),3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE
        CASE(1); x=RNDM(kk);
           IF(x>half)THEN; vert1=name1; vert2=new_name1;
           ELSE;           vert1=new_name1; vert2=name0; ENDIF
        CASE(2); x=RNDM(kk)
           IF(x>half)THEN; vert1=name0; vert2=new_name2;
           ELSE;           vert1=new_name2; vert2=name2; ENDIF
        END SELECT
        ! Checking weight change
        Energy_a=Energy_a * ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
!        IF(site_new1/=site_new2)THEN
!            PRINT*,site1,site0,site2
!            PRINT*,site_new1,site_new2
!        ENDIF    
      ENDIF

      END SUBROUTINE ADD_PHONON_DIA_EXP
!....................................................................

!--------------------------------------------------------------------
! Undressing diagonal phonon
!--------------------------------------------------------------------
      SUBROUTINE REM_PHONON_DIA_EXP
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio
      REAL*8 :: sign0
      REAL*8 :: pr1, del_time

      km=15; add_rem_dia_exp=add_rem_dia_exp+un1
!....................................................................
!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  from
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                          .               .
!                           .           .
!                              .      .
!                               . . .


! Qualifuing for undressing the vertes
      IF(nmnm==2) RETURN            ! no phonons
      ! Select vertex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      new_name1=link(name0,1); IF(new_name1==name0)RETURN; !It is Hartree
      IF(type(new_name1,3)==0)RETURN  !Worm atached to vothex "new_name1"
      IF(MEASV .AND. (new_name1==vert1.OR.new_name1==vert2))RETURN !it is measuring line
      new_name2=link(name0,2); IF(new_name2==name0)RETURN; !Hartree (redundant)
      IF(type(new_name2,3)==0)RETURN  !worm at "new_name2" (redundant)
      IF(link(new_name1,3)/=new_name2)RETURN ! not a dressed vertex
      name2=link(new_name2,2); IF(name2==new_name1)STOP'FockReducible!'
      name1=link(new_name1,1); IF(name1==new_name2)STOP'FockRedu..!'!(redund..)
! If qualified -> continue

! Defining parameters of lines
      time_new1=tau(new_name1); time_new2=tau(new_name2)
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site_new1=site(new_name1);  site_new2=site(new_name2)
      IF(site_new1/=site_new2)RETURN
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta); pro=pro/(Nsite);
      pro=pro*(un1*nmnm-un2+un2)/(un1*nmnm-un2)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)
! Definin ratio without context
      ratio_0 =  weight_new / weight_old

! Checking if one of propagators in block are measuring
      caseca=0
      IF(.NOT. MEASV)THEN
        IF(name1==vert1)caseca=1; IF(new_name1==vert1)caseca=2
        IF(name0==vert1)caseca=3; IF(new_name2==vert1)caseca=4
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(3); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(4); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT    
      
! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS( ratio_0 / pro )
! Adding flat ordest distribution factors
      ratio = ratio * (ord_probab(nmnm-2+2) / ord_probab(nmnm-2))
! Correcting ratio by exponential distribution time weight function
      ! searching shortest distance between time points
      del_time=time_new2-time_new1
      IF(del_time<nul)del_time=beta+del_time ! making time difference positive
      
!      IF(RNDM(kk)>half)THEN
!           IF(del_time>beta*half)del_time=beta-del_time !shortest distance 
!      ENDIF     
      !IF(del_time>beta*half)THEN;
      !    PRINT*,'Del_time = ',del_time; STOP"Nu i kak s takim?"
      !ENDIF 
      bazeka = bazeka_min + (RNDM(kk)*(bazeka_max-bazeka_min))
   !   pr1=un1/(beta)
      pr1 = bazeka * EX(-bazeka*del_time) / (un1-EX(-bazeka*beta)) 
      ratio = ratio / pr1 ; !PRINT*,"ANTIRAT: ",ratio
! Reverting ration for remove      
      ratio = un1 / ratio  ! because it is remove
! Playing metropolis
      a00_rem_dia_exp=a00_rem_dia_exp+un1
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; suc_rem_dia_exp=suc_rem_dia_exp+un1
          !PRINT*,'succed to delete****************'
       ! set links
        link(name1,2)=name0; link(name0,1)=name1;
        link(name0,2)=name2; link(name2,1)=name0;
        ! setting momenta
        omega(name0,2,1:3)=omega(name2,1,1:3)
        omega(name0,1,1:3)=omega(name1,2,1:3)
        ! removing momenta from hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL DROPHASH(2,new_name1)
        CALL DROPHASH(2,new_name2)
        CALL DROPHASH(3,new_name1)
        CALL DROPHASH(3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE;
        CASE(1); vert2=name0;
        CASE(2); vert1=name1;
        CASE(3); vert2=name2;
        CASE(4); vert1=name0;
        END SELECT;
        CALL DROPNAME(new_name1); CALL DROPNAME(new_name2)
        ! Checking weight change
        Energy_a=Energy_a / ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
      ENDIF

          END SUBROUTINE REM_PHONON_DIA_EXP
!....................................................................

!--------------------------------------------------------------------
! Dressin vortex by diagonal phonon propagator
!--------------------------------------------------------------------
      SUBROUTINE ADD_PHONON_DIA_SHO
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN, PHONON, RNDM, MODA, EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio, x
      REAL*8 :: sign0
      REAL*8 :: RRd, pr1, del_time, bakara

      km=16; add_add_dia_sho=add_add_dia_sho+un1

!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  to
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                         .                .
!                           .           .
!                            .      .
!                             . . .
!
 
! Limiting diagramm order
      IF(nmnm>=N_max_now)RETURN
! Selecting the vortex to dress
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      name1=link(name0,1); name2=link(name0,2)
! Avoiding dressing Hartree
      IF(name0==name1)RETURN

! Suggesting times of new vortices
      !first suggested time is random
      time_new1=beta*RNDM(kk); 
      !second suggeted time is an exponential distribution
      bazeka = bazeka_min + (RNDM(kk)*(bazeka_max-bazeka_min))
      RRd=RNDM(kk)
      bakara = un1-(RRd*(un1-EX(-bazeka*u_limit)))
      IF(bakara>1.0d-40)THEN 
          del_time = - (un1/bazeka) * LOG( bakara ) 
      ELSE
          del_time = RRd * u_limit 
      ENDIF    
      IF(del_time>u_limit)THEN;
          PRINT*,"RRd = ,",RRd
          PRINT*,'Del_time = ',del_time
          STOP"Panic in del_time"
      ENDIF    
      time_new2 = time_new1 + del_time; 
      IF(time_new2>beta)time_new2=time_new2-beta
! Suggesting sites of new vortices
      site_new1=Nsite*RNDM(kk)+1; IF(site_new1>Nsite)site_new1=Nsite;
      site_new2=site_new1

! Defining parameters of name1 and name2
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! Defininf context factor
      pro=un1/(beta); pro=pro/(Nsite);
      pro=pro*(nmnm+un2)/(un1*nmnm)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)

! Definin gratio without context
      ratio_0 = weight_new / weight_old

! Checking whether one of propagators connected to name0 are measuring
      caseca=0
      IF(.NOT. MEASV)THEN;
        IF(name1==vert1)caseca=1
        IF(name0==vert1)caseca=2
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT

! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS (ratio_0 / pro)
! Adding flat orders distribution factors
      ratio = ratio * (ord_probab(nmnm+2) / ord_probab(nmnm))
! Correcting ratio by exponential distribution time weight function
      ! exponential distribution
      pr1 = bazeka * EX(-bazeka*del_time) / (un1-EX(-bazeka*u_limit)) 
      ratio = ratio / pr1; !PRINT*,"ADD: ",ratio
! Playing metropolis
      a00_add_dia_sho=a00_add_dia_sho+un1
      propell=EVA_RAT(kk,ratio)
      
! Updating
      IF(propell)THEN; suc_add_dia_sho=suc_add_dia_sho+un1; 
        CALL GETNAME(new_name1); CALL GETNAME(new_name2)
        ! set tau and site
        tau(new_name1)=time_new1; tau(new_name2)=time_new2
        site(new_name1)=site_new1; site(new_name2)=site_new2;
        ! set types
        type(new_name1,1)=type1; type(new_name1,2)=type1;
        type(new_name1,3)=1
        type(new_name2,1)=type2; type(new_name2,2)=type2;
        type(new_name2,3)=1
        ! set links
        link(name1,2)=new_name1; link(name0,1)=new_name1;
        link(name0,2)=new_name2; link(name2,1)=new_name2;
        link(new_name1,1)=name1; link(new_name1,2)=name0;
        link(new_name1,3)=new_name2;
        link(new_name2,1)=name0; link(new_name2,2)=name2;
        link(new_name2,3)=new_name1
        ! set momenta
        DO i=1,3;  omega(new_name1,3,i)=RNDM(kk); ENDDO
        omega(new_name2,3,1:3)=-omega(new_name1,3,1:3);
        DO i=1,3;
          omega(new_name1,2,i)=MODA(omega(name1,2,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,1,1:3)=omega(new_name1,2,1:3)
        DO i=1,3
          omega(new_name2,1,i)=MODA(omega(name2,1,i)-omega(new_name1,3,i))
        ENDDO
        omega(name0,2,1:3)=omega(new_name2,1,1:3)
        omega(new_name1,1,1:3)=omega(name1,2,1:3)
        omega(new_name2,2,1:3)=omega(name2,1,1:3)
        ! putting momenta into hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL ADDHASH(omega(new_name1,2,1),2,new_name1)
        CALL ADDHASH(omega(new_name2,2,1),2,new_name2)
        CALL ADDHASH(DABS(omega(new_name1,3,1)),3,new_name1)
        CALL ADDHASH(DABS(omega(new_name2,3,1)),3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE
        CASE(1); x=RNDM(kk);
           IF(x>half)THEN; vert1=name1; vert2=new_name1;
           ELSE;           vert1=new_name1; vert2=name0; ENDIF
        CASE(2); x=RNDM(kk)
           IF(x>half)THEN; vert1=name0; vert2=new_name2;
           ELSE;           vert1=new_name2; vert2=name2; ENDIF
        END SELECT
        ! Checking weight change
        Energy_a=Energy_a * ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
!        IF(site_new1/=site_new2)THEN
!            PRINT*,site1,site0,site2
!            PRINT*,site_new1,site_new2
!        ENDIF    
      ENDIF

      END SUBROUTINE ADD_PHONON_DIA_SHO
!....................................................................

!--------------------------------------------------------------------
! Undressing diagonal phonon
!--------------------------------------------------------------------
      SUBROUTINE REM_PHONON_DIA_SHO
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM,EX
      LOGICAL,EXTERNAL :: EVA_RAT
      LOGICAL :: propell
      INTEGER :: num, name0, name1, name2, new_name1,new_name2
      INTEGER :: site0, site1, site2, site_new1, site_new2
      INTEGER :: type1, type2, i, caseca
      REAL*8 :: time0, time1, time2, time_new1, time_new2
      REAL*8 :: pro, weight_old, weight_new, ratio_0, ratio
      REAL*8 :: sign0
      REAL*8 :: pr1, del_time

      km=17; add_rem_dia_sho=add_rem_dia_sho+un1
!....................................................................
!
!    name1     type1     name0    type2           name2
!      O--->---------------X----------------->-----O
!     time1              time0                    time2
!     site1              site0                    site2
!
!  from
!
!               new_name1     new_name2
!               site_new1     site_new2
!               time_new1     time_new2
!      O--->--------X------X------X----------->-----O
!                          .               .
!                           .           .
!                              .      .
!                               . . .


! Qualifuing for undressing the vertes
      IF(nmnm==2) RETURN            ! no phonons
      ! Select vertex
      num=nmnm*RNDM(kk)+1; IF(num>nmnm)num=nmnm; name0=nlist(num)
      new_name1=link(name0,1); IF(new_name1==name0)RETURN; !It is Hartree
      IF(type(new_name1,3)==0)RETURN  !Worm atached to vothex "new_name1"
      IF(MEASV .AND. (new_name1==vert1.OR.new_name1==vert2))RETURN !it is measuring line
      new_name2=link(name0,2); IF(new_name2==name0)RETURN; !Hartree (redundant)
      IF(type(new_name2,3)==0)RETURN  !worm at "new_name2" (redundant)
      IF(link(new_name1,3)/=new_name2)RETURN ! not a dressed vertex
      name2=link(new_name2,2); IF(name2==new_name1)STOP'FockReducible!'
      name1=link(new_name1,1); IF(name1==new_name2)STOP'FockRedu..!'!(redund..)
! If qualified -> continue

! Defining parameters of lines
      time_new1=tau(new_name1); time_new2=tau(new_name2)
      time0=tau(name0); time1=tau(name1); time2=tau(name2)
      site_new1=site(new_name1);  site_new2=site(new_name2)
      IF(site_new1/=site_new2)RETURN
      site0=site(name0); site1=site(name1); site2=site(name2);
      type1=type(name0,1); type2=type(name0,2);

! searching shortest distance between time points
      del_time=time_new2-time_new1
      IF(del_time<nul)del_time=beta+del_time ! making time difference positive
      IF(del_time>u_limit)RETURN ! CAN not remove because could not create

! Defininf context factor
      pro=un1/(beta); pro=pro/(Nsite);
      pro=pro*(un1*nmnm-un2+un2)/(un1*nmnm-un2)

! Definig old weight
      weight_old = GREEN(site1,time1,site0,time0,type1) * GREEN(site0,time0,site2,time2,type2)
! Definin new weight
      weight_new = GREEN(site1,time1,site_new1,time_new1,type1) * &
                  GREEN(site_new1,time_new1,site0,time0,type1) * &
                  GREEN(site0,time0,site_new2,time_new2,type2) * &
                  GREEN(site_new2,time_new2,site2,time2,type2) * &
                  PHONON(site_new1,time_new1,site_new2,time_new2)
! Definin ratio without context
      ratio_0 =  weight_new / weight_old

! Checking if one of propagators in block are measuring
      caseca=0
      IF(.NOT. MEASV)THEN
        IF(name1==vert1)caseca=1; IF(new_name1==vert1)caseca=2
        IF(name0==vert1)caseca=3; IF(new_name2==vert1)caseca=4
      ENDIF
! Changing context because of selecting one of two
      SELECT CASE(caseca)
          CASE(0); CONTINUE
          CASE(1); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(2); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(3); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
          CASE(4); CONTINUE; pro=pro/2.0 !pro = pro/0.5d0
      END SELECT    
      
! Adding context factors
      IF(ratio_0>nul)THEN; sign0=1; ELSE; sign0=-1; ENDIF
      ratio = ABS( ratio_0 / pro )
! Adding flat ordest distribution factors
      ratio = ratio * (ord_probab(nmnm-2+2) / ord_probab(nmnm-2))
! Correcting ratio by exponential distribution time weight function
      bazeka = bazeka_min + (RNDM(kk)*(bazeka_max-bazeka_min))
      pr1 = bazeka * EX(-bazeka*del_time) / (un1-EX(-bazeka*u_limit)) 
      ratio = ratio / pr1 ; !PRINT*,"ANTIRAT: ",ratio
! Reverting ration for remove      
      ratio = un1 / ratio  ! because it is remove
! Playing metropolis
      a00_rem_dia_sho=a00_rem_dia_sho+un1
      propell=EVA_RAT(kk,ratio)

! Updating
      IF(propell)THEN; suc_rem_dia_sho=suc_rem_dia_sho+un1
          !PRINT*,'succed to delete****************'
       ! set links
        link(name1,2)=name0; link(name0,1)=name1;
        link(name0,2)=name2; link(name2,1)=name0;
        ! setting momenta
        omega(name0,2,1:3)=omega(name2,1,1:3)
        omega(name0,1,1:3)=omega(name1,2,1:3)
        ! removing momenta from hash table
        CALL DROPHASH(2,name0); CALL ADDHASH(omega(name0,2,1),2,name0)
        CALL DROPHASH(2,new_name1)
        CALL DROPHASH(2,new_name2)
        CALL DROPHASH(3,new_name1)
        CALL DROPHASH(3,new_name2)
        ! redefining the measuring line
        SELECT CASE(caseca)
        CASE(0); CONTINUE;
        CASE(1); vert2=name0;
        CASE(2); vert1=name1;
        CASE(3); vert2=name2;
        CASE(4); vert1=name0;
        END SELECT;
        CALL DROPNAME(new_name1); CALL DROPNAME(new_name2)
        ! Checking weight change
        Energy_a=Energy_a / ABS(ratio_0)
        Dphase_a=Dphase_a*sign0
      ENDIF

      END SUBROUTINE REM_PHONON_DIA_SHO
!....................................................................

          
          
!--------------------------------------------------------------------
! CALCULATING factor of all propagators: absolute value
!--------------------------------------------------------------------
      SUBROUTINE energytest
      USE config_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN, PHONON, PHONON_MEAS
      LOGICAL,EXTERNAL :: IS_P_MEASURING
      REAL*8 :: x
      INTEGER :: i, name1, name2, spin
      LOGICAL :: notdone(N_max)

      energy=1.d0

! product of all propagators
      DO i=1,nmnm; name1=nlist(i); name2=link(name1,2)
      spin=type(name1,2)
      x=green(site(name1),tau(name1),site(name2),tau(name2),spin)
      energy=energy*x
      ENDDO

      DO i=1,nmnm; name1=nlist(i);
        notdone(name1)=.TRUE. ;
      ENDDO
      DO i=1,nmnm; name1=nlist(i); name2=link(name1,3)
        IF(notdone(name1)) THEN;
          notdone(name1)=.FALSE. ; notdone(name2)=.FALSE.
          IF( IS_P_MEASURING(name1,name2) )THEN
             x= PHONON_MEAS(site(name1),tau(name1),site(name2),tau(name2)) 
          ELSE    
             x= PHONON(site(name1),tau(name1),site(name2),tau(name2))
          ENDIF    
          energy=energy*x;
        ENDIF
      ENDDO   ! product of all phonon lines

      END SUBROUTINE energytest
!....................................................................

!--------------------------------------------------------------------
! CALCULATING phase of all propagators: sign
!--------------------------------------------------------------------
      SUBROUTINE phasetest
      USE config_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON
      REAL*8 :: x
      INTEGER:: i, iloops, k, name1, name2, spin
      LOGICAL :: loops(N_max), notdone(N_max)

      Dphase=1.d0

      DO i=1,nmnm; name1=nlist(i); name2=link(name1,2)
        spin=type(name1,2)
        x=green(site(name1),tau(name1),site(name2),tau(name2),spin)
        Dphase=Dphase*Gphase
      ENDDO    ! product of all propagators

      DO i=1,nmnm; name1=nlist(i);
         notdone(name1)=.TRUE. ;
      ENDDO
      DO i=1,nmnm; name1=nlist(i); name2=link(name1,3)
        IF(notdone(name1)) THEN;
          notdone(name1)=.FALSE. ; notdone(name2)=.FALSE.
          x=phonon(site(name1),tau(name1),site(name2),tau(name2))
          Dphase=Dphase*Vphase
        ENDIF
      ENDDO    ! product of all interaction lines

      iloops=0; loops=.TRUE.
      DO i=1,nmnm; name1=nlist(i);
        IF(loops(name1)) then; iloops=iloops+1
          name2=name1; do; loops(name2)=.FALSE.
          name2=link(name2,2); IF(name2==name1) EXIT; enddo
        ENDIF
      ENDDO               ! number of fermionic loops in the graph

      k=iloops/2 ; k=k*2  ! odd or even?
      IF(k.ne.iloops)  Dphase=-Dphase

      END SUBROUTINE phasetest
!....................................................................

!C-----------------------------------------------------
!C     DOING VARIOUS DEBUGGING SUBROUTINE
!C-----------------------------------------------------
      SUBROUTINE DEBU_DOING;
      USE config_par
      IMPLICIT NONE
      LOGICAL :: select
      INTEGER site_1,site_2
      INTEGER,DIMENSION(3) :: vect_1,vect_2,intervect

      select=.FALSE.
! Checking vectors
      IF(select)THEN
      PRINT*,"Enter number of first and second site:"
      READ*,site_1,site_2
      CALL VEC_BETWEEN(site_1,site_2,vect_1,vect_2,intervect)
      PRINT*,"Vect 1: ",vect_1
      PRINT*,"Vect 2: ",vect_2
      PRINT*,"1--> 2: ",intervect
      ENDIF
! Chcking vectors end

      END SUBROUTINE DEBU_DOING
!C....................................................

!--------------------------------------------------------------------
! Time phonon propagator - BARE
!--------------------------------------------------------------------
      REAL*8 FUNCTION PHOZERO(tutu)
      USE config_par; USE stat_par; IMPLICIT NONE
      REAL*8,INTENT(IN) :: tutu
      REAL*8,EXTERNAL :: EX
      REAL*8 :: a1, a2, a3, pr1, pr2, tata

      tata=tutu; IF(tutu<nul)tata=beta+tutu
      a1=EX(beta*Debye); a2=EX(-tata*Debye); a3=1.0d0/a2
      pr1 = a1*a2 + a3
      pr2 = a1 - 1.0d0
      phozero = pr1 / pr2
      phozero = V_kvadrat * phozero

      END FUNCTION PHOZERO
!....................................................................
      
!--------------------------------------------------------------------
! Time phonon propagator for measurement line
!--------------------------------------------------------------------
      REAL*8 FUNCTION PHOZERO_MEAS(tutu)
      USE config_par; USE stat_par; IMPLICIT NONE
      REAL*8,INTENT(IN) :: tutu
      REAL*8,EXTERNAL :: EX
      REAL*8 :: a1, a2, a3, pr1, pr2, tata

      tata=tutu; IF(tutu<nul)tata=beta+tutu
      a1 = EX(  beta * (Debye/meas_Debye) ); 
      a2 = EX( - tata * (Debye/meas_Debye) ); 
      a3=1.0d0/a2
      pr1 = a1*a2 + a3
      pr2 = a1 - 1.0d0
      phozero_meas = pr1 / pr2
      phozero_meas = (V_kvadrat*meas_V_kvadrat) * phozero_meas
      
      END FUNCTION PHOZERO_MEAS
!....................................................................

!--------------------------------------------------------------------
!   Phonon propagator initial table
!--------------------------------------------------------------------
      SUBROUTINE PHONO_0
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: PHOZERO, COULOMB
      REAL*8 :: tt
      INTEGER :: it, ix, iy, iz, i

      print*, 'in PHONO_0'

      PHO_0(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau)=nul
      
      IF(i_type==0)THEN
          DO it=0,Ntau; tt=dtau*it;
              PHO_0(0,0,0,it) = PHOZERO(tt); 
          ENDDO
      ELSE IF(i_type==1)THEN
          DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1; DO it=0,Ntau; tt=dtau*it;
  	          PHO_0(ix,iy,iz,it) = PHOZERO(tt) * COULOMB(ix,iy,iz);
          ENDDO; ENDDO; ENDDO; ENDDO;
      ELSE
          STOP"We do not know such interaction type"
      ENDIF   
      
      OPEN(UNIT=4,FILE="cou.dat")
          DO ix=0,L(1)-1; DO iy=0,L(2)-1;
              WRITE(4,*)ix,iy,cou_cou(ix,iy)
          ENDDO; ENDDO;    
      CLOSE(4)
          
      PHO_NOW(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau)= PHO_0(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau)
      
      !STOP"Looking radii"

      theta_meas_2 = 0.0d0
      DO i=1,dim
         theta_meas_2 = theta_meas_2 + ( ( 1.0d0*(L(i)-1) ) / 2.0 )**2
      ENDDO    
      theta_meas_2 = SQRT(theta_meas_2) / meas_decay
      theta_meas_2 = theta_meas_2 ** 2
      
      print*, 'PHONO_0 DONE'
      PRINT*,' theta_meas_2: ', theta_meas_2
      END SUBROUTINE PHONO_0
!........................................................


!--------------------------------------------------------------------
!   Analytic norm of Hartree+phonon  ---------O
!--------------------------------------------------------------------
      REAL*8 FUNCTION S_NORM_0
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: PHONON,GREEN
      REAL*8,PARAMETER :: malo=1.0d-14
      REAL*8 :: tt, ttha, hastep, pr0
      INTEGER :: it, is

      hastep = half*dtau

      S_norm_0=nul;
      pr0 = un2 * ABS(GREEN(1,nul,1,nul,1))
      DO is=1,Nsite; DO it=0,Ntau-1; tt=dtau*it;
        ttha = tt + hastep
        S_norm_0 = S_norm_0 + ABS(PHONON(1,nul,is,ttha))
!        PRINT*,is,PHONON(1,nul,is,ttha),GREEN(1,nul,1,nul,1); PAUSE
      ENDDO; ENDDO;

      S_norm_0 =  S_norm_0 * pr0 * dtau

!      print*,'S_norm =',S_norm

      END FUNCTION S_NORM_0
!........................................................

!--------------------------------------------------------------------
!   Analytic norm of Hartree+phonon Fock
!--------------------------------------------------------------------
      REAL*8 FUNCTION S_NORM
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: PHONON,GREEN
      REAL*8,PARAMETER :: malo=1.0d-14
      REAL*8 :: tt, ttha, hastep, pr0
      INTEGER :: it, is

      hastep = half*dtau

      S_norm=nul
      DO is=1,Nsite; DO it=0,Ntau-1; tt=dtau*it;
        ttha = tt + hastep
        S_norm = S_norm + ABS(GREEN(1,nul,is,ttha,1) * PHONON(1,nul,is,ttha))
      ENDDO; ENDDO;
      
      IF(.NOT. POLARIZED)THEN
        ! NOT POL_SPIN------
        S_norm = un2 * S_norm * dtau
        ! NOT POL_SPIN------
      ELSE
        ! POL_SPIN------
        S_norm =       S_norm * dtau
        ! POL_SPIN------
      ENDIF
      
      END FUNCTION S_NORM
!........................................................
      

!--------------------------------------------------------------------
!   Analytic Normalization of polarization operator
!--------------------------------------------------------------------
      REAL*8 FUNCTION P_NORM
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: PHONON,GREEN
      REAL*8,PARAMETER :: malo=1.0d-14
      REAL*8 :: tt, ttha, hastep
      INTEGER :: it, is

      hastep = half*dtau

      P_norm=nul
      DO is=1,Nsite; DO it=0,Ntau-1; tt=dtau*it;
        ttha = tt + hastep
        P_norm = P_norm + ABS(GREEN(1,nul,1,ttha,1) * GREEN(1,ttha,1,nul,1))       ! Only local
!           ABS(GREEN(1,nul,is,ttha,1) * GREEN(is,ttha,1,nul,1))    ! Local+Nonlocal
      ENDDO; ENDDO;

      IF(.NOT. POLARIZED)THEN
        ! NON POL_SPIN------
        P_norm = un2 * P_norm * dtau / Nsite 
        ! NON POL_SPIN------
      ELSE 
        ! POL_SPIN------
        P_norm =        P_norm * dtau / Nsite      
        ! POL_SPIN------
      ENDIF 
      
!      print*,'P_norm =',P_norm

      END FUNCTION P_NORM
!........................................................

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!              MESURING SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE MEASURE_0
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN,PHONON,RNDM
      INTEGER :: s1, s2, s, k, v(3), v_1(3), v_2(3)
      INTEGER :: numa
      REAL*8 :: t1, t2, factor, at, znak, ukaka

!      STOP"We do nit use this subroutine now"
      
      
      IF(nmnm.le.N_mes_now) THEN       ! measures when order<N_mes
                                       ! one must keep N_mes_now=2
                                       ! to check the test for
                                       ! polarization operator
      IF(.not.present) THEN            ! measures when no worm
            
        s1=site(vert1); s2=site(vert2)
        t1=tau(vert1);  t2=tau(vert2)

        IF(.not.measV) THEN                            ! self-energy

          s=type(vert1,2)                                ! spin
          factor=un1/green(s1,t1,s2,t2,s)
          factor=factor*Dphase_a
          factor = factor 
          CALL VEC_BETWEEN(s2,s1,v_1,v_2,v)           ! space index
          at=t1-t2; znak=1.d0;                             ! tau>0
          IF(at<0.d0) then; at=beta+at; znak=-1.d0; ENDIF  ! tau<0
          k=at/dtau                                        ! index
          IF(vert1.ne.vert2)THEN !Filtering Fock
            sigma_box(v(1),v(2),v(3),k)= sigma_box(v(1),v(2),v(3),k)-znak*factor/2.d0
          ENDIF                                  ! counting both spins
          IF(nmnm==2 .AND. vert1==vert2)THEN
            Z_self_norm=Z_self_norm+DABS(factor)/2.d0  ! counter for Hartree
          ENDIF

          ELSE                                           ! polarization
              
          IF(nmnm>2 .OR. vert1==link(vert2,1))THEN
            factor=un1/phonon(s1,t1,s2,t2)
            factor=factor*Dphase_a
            factor=factor
            CALL VEC_BETWEEN(s2,s1,v_1,v_2,v)           ! space index
            at=t1-t2; IF(at<0.d0) at=beta+at;  k=at/dtau   ! tau>0 index
            Polar_box(v(1),v(2),v(3),k)= Polar_box(v(1),v(2),v(3),k)-factor
          ENDIF
          IF(nmnm==2 .AND. vert1==link(vert2,1))THEN
             factor=un1/phonon(s1,t1,s2,t2)
             Z_pola_norm=Z_pola_norm+DABS(factor)         ! counter for Fock
          IF(vert1==link(vert1,1))STOP'impossible here'
          ENDIF

        ENDIF                                          ! self or pol
        
        ENDIF                               ! measures when no worm
      ENDIF                               ! measures when order < N_mes_now

      END SUBROUTINE MEASURE_0
!-----------------------------------------------------------------------
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!              MESURING SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE MEASURE
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN, PHONON, PHONON_MEAS, BO_SIGMA, BO_POLAR
      REAL*8,EXTERNAL :: FACT, BO_OC
      INTEGER :: s1, s2, s, k, v(3), v_1(3), v_2(3)
      INTEGER :: kb, ib
      REAL*8 :: t1, t2, factor, at, znak, yc, haba
      
      num_vort: IF(nmnm.le.N_mes_now)THEN ! measures when order<N_mes, keep N_mes_now=2 to compare with Pol Operator                                  ! to check the test for polarization operator
      present_worm: IF(.not.present) THEN ! measures when no worm
           
            
         s1=site(vert1); s2=site(vert2)
         t1=tau(vert1);  t2=tau(vert2)

         IF(.not.measV) THEN                            ! self-energy sector

            s=type(vert1,2)                                ! spin
            factor=un1/green(s1,t1,s2,t2,s)
            factor=factor*Dphase_a
            factor = factor / ord_probab(nmnm)
            CALL VEC_BETWEEN(s2,s1,v_1,v_2,v)           ! space index
            at=t1-t2; znak=1.d0;                             ! tau>0
            IF(at<0.d0) then; at=beta+at; znak=-1.d0; ENDIF  ! tau<0
            k=at/dtau                                        ! index
            fock: IF(vert1.ne.vert2)THEN !Filtering Fock
              IF(.NOT. POLARIZED)THEN  
                 ! NOT POL_SPIN------
                 sigma_box(v(1),v(2),v(3),k)= sigma_box(v(1),v(2),v(3),k)-znak*factor/2.d0
                 kb=at/dtb_s; kb=min(kb,Nbins_s-1)                 ! tau index basis
                 haba = -znak*factor/2.d0
                 DO ib=1,Nbasis_s      
                    yc=BO_SIGMA(ib,at,kb)*haba
                    Sigmab(v(1),v(2),v(3),kb,ib) = Sigmab(v(1),v(2),v(3),kb,ib) + yc
                 ENDDO 
                 ! NOT POL_SPIN------
              ELSE
                 ! POL_SPIN------
                 sigma_box(v(1),v(2),v(3),k)= sigma_box(v(1),v(2),v(3),k)-znak*factor
                 kb=at/dtb_s; kb=min(kb,Nbins_s-1)                 ! tau index basis
                 haba = -znak*factor
                 DO ib=1,Nbasis_s      
                    yc=BO_SIGMA(ib,at,kb)*haba
                    Sigmab(v(1),v(2),v(3),kb,ib) = Sigmab(v(1),v(2),v(3),kb,ib) + yc
                 ENDDO 
                 ! POL_SPIN------
              ENDIF
            ENDIF fock                                 ! counting both spins
            IF(nmnm==2 .AND. vert1/=vert2)THEN
               Z_self_norm=Z_self_norm+DABS(factor)/1.d0  ! counter for Fock
            ENDIF

         ELSE                                      ! polarization
              
            IF(MEASURE_OC)THEN; 
               i_oc_count=i_oc_count+1
               IF(i_oc_count==measu_oc_limit)THEN;
                    CALL MEASURING_OC(t1,t2,s1,s2);
                    i_oc_count=0;
               ENDIF     
            ENDIF;
                        
            IF(nmnm>2 .OR. vert1==link(vert2,1))THEN  
               factor=un1/phonon_meas(s1,t1,s2,t2)
               factor=factor*Dphase_a
               factor = factor / ord_probab(nmnm)
               CALL VEC_BETWEEN(s2,s1,v_1,v_2,v)           ! space index
               at=t1-t2; IF(at<0.d0) at=beta+at;  k=at/dtau   ! tau>0 index
               Polar_box(v(1),v(2),v(3),k)= Polar_box(v(1),v(2),v(3),k)+factor
               kb=at/dtb_p; kb=min(kb,Nbins_p-1)                 ! tau index basis    
               DO ib=1,Nbasis_p      
                  yc=BO_POLAR(ib,at,kb)*factor
                  Polarb(v(1),v(2),v(3),kb,ib) = Polarb(v(1),v(2),v(3),kb,ib) + yc
               ENDDO 
            ENDIF
            IF(nmnm==2 .AND. vert1==link(vert2,1))THEN
               IF(s1==s2)THEN; !only local phonons 
                 factor=un1/phonon_meas(s1,t1,s2,t2)
                 factor = factor / ord_probab(nmnm)
                 Z_pola_norm=Z_pola_norm+DABS(factor)      ! counter for Fock
               ENDIF   
               IF(vert1==link(vert1,1))STOP'impossible here'
            ENDIF

         ENDIF                                          ! self or pol
        
      ENDIF present_worm     ! measures when no worm
      ENDIF num_vort         ! measures when order < N_mes_now

      END SUBROUTINE MEASURE
!-----------------------------------------------------------------------
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!            LOOP  MESURING SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE LOOP_MEASURE
      USE config_par; USE stat_par
      IMPLICIT NONE
      LOGICAL :: loops(N_max);
      INTEGER :: iloops, i, name1, name2;
      REAL*8 :: pula=1.0d-5
      
      iloops=0; loops=.TRUE.
      DO i=1,nmnm; name1=nlist(i);
        IF(loops(name1)) then; iloops=iloops+1
          name2=name1; do; loops(name2)=.FALSE.
          name2=link(name2,2); IF(name2==name1) EXIT; enddo
        ENDIF
      ENDDO               ! number of fermionic loops in the graph
      
      !PRINT*,"Loops number = ",iloops
      
      loop_summed = loop_summed + ( pula*iloops / ord_probab(nmnm) ); 
      time_summed = time_summed + ( pula        / ord_probab(nmnm) );

      
      END SUBROUTINE LOOP_MEASURE
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
!   SUBROUTINE MEASURING OC <jj>
!--------------------------------------------------------------------
      SUBROUTINE MEASURING_OC(t1,t2,s1,s2)
      USE config_par; USE stat_par
      IMPLICIT NONE;
      REAL*8,INTENT(IN) :: t1,t2
      INTEGER,INTENT(IN) :: s1,s2
      REAL*8,EXTERNAL :: GREEN, PHONON_MEAS, BO_OC
      INTEGER,DIMENSION(3) :: v1, v2, v12, sh_i, sh_j, v, v_1, v_2
      INTEGER,DIMENSION(4) :: rsit_outgo, rsit_ingo
      INTEGER,DIMENSION(4) :: lsit_outgo, lsit_ingo
      INTEGER :: i, j, ii, iii, k, isensio, ibi
      INTEGER:: scalar, index_ri, sp_in_1, sp_in_2
      INTEGER :: sp_r, sp_l, ik, v3, v4, vert_ri, vert_le
      INTEGER :: kb, ib
      REAL*8 :: delta_tau, t3, t4, tau_ri, tau_le, at, fabule, factor_all_include
      REAL*8 :: rig_dist, le_fra, ri_fra, de_le, de_ri, buzota, oanwi, yc
      REAL*8 :: factor, factor_out, factor_in, factor_1to2, factor_2to1 
      REAL*8 :: de_de, ltl_ingo, ltl_outgo, rtl_ingo, rtl_outgo
      INTEGER :: rl_ingo, rl_outgo, rsl_ingo, rsl_outgo
      INTEGER :: ll_ingo, ll_outgo, lsl_ingo, lsl_outgo
      INTEGER :: s1_sh, s2_sh, sp_1, sp_2, s3, s4, ipm
      
      ! selecting quntities for OC estimator    
      sp_1=type(vert1,1); sp_2=type(vert2,1)
      IF(t1>t2)THEN
            vert_ri=vert1; vert_le=vert2; index_ri=1
            tau_ri=tau(vert1); tau_le=tau(vert2)
            sp_r=sp_1; sp_l=sp_2;
          ELSE
            vert_ri=vert2; vert_le=vert1; index_ri=2
            tau_ri=tau(vert2); tau_le=tau(vert1)
            sp_l=sp_1; sp_r=sp_2;
          ENDIF
          delta_tau=tau_ri-tau_le; !PRINT*,delta_tau
          !selecting left and right fractions for elongations 
          rig_dist=beta-tau_ri; 
          le_fra = tau_le / (tau_le+rig_dist)
          ri_fra = rig_dist / (tau_le+rig_dist)
          !establishing associations
          rl_ingo = link(vert_ri,1)
          rtl_ingo = tau(rl_ingo); 
          rsl_ingo = site(rl_ingo) ! <- site can be changed if direct link vert_ri - vert_le
          rl_outgo = link(vert_ri,2)
          rtl_outgo = tau(rl_outgo); 
          rsl_outgo =  site(rl_outgo) ! <- site can be changed if direct link vert_ri - vert_le
          ll_ingo = link(vert_le,1)
          ltl_ingo = tau(ll_ingo); 
          lsl_ingo = site(ll_ingo) ! <- site can be changed if direct link vert_ri - vert_le
          ll_outgo = link(vert_le,2)
          ltl_outgo = tau(ll_outgo); 
          lsl_outgo =  site(ll_outgo) ! <- site can be changed if direct link vert_ri - vert_le
        
! <jj> estimator         
          !IF(nmnm>2)THEN
          dir_1: DO i=1,dir
          dir_2: DO j=1,dir
          
            sp_1=type(vert1,2); sp_2=type(vert2,2) !outgoing spins from verts
            sp_in_1=type(vert1,1); sp_in_2=type(vert2,1) !ingoing spins to verts
            !finding jumps of two neighboring points
            s1_sh=ass(s1,i)
            s2_sh=ass(s2,j)
            CALL VEC_BETWEEN(s1,s1_sh,v1,v2,v12)
            sh_i(1:dim)=v2(1:dim)-v1(1:dim)
            DO iii=1,dim
              IF(sh_i(iii)>=2)THEN
                sh_i(iii)=-1
              ELSE IF(sh_i(iii)<=-2)THEN
                sh_i(iii)= 1    
              ENDIF     
            ENDDO              
            CALL VEC_BETWEEN(s2,s2_sh,v1,v2,v12)
            sh_j(1:dim)=v2(1:dim)-v1(1:dim)
            DO iii=1,dim
              IF(sh_j(iii)>=2)THEN
                sh_j(iii)=-1
              ELSE IF(sh_j(iii)<=-2)THEN
                sh_j(iii)= 1    
              ENDIF     
            ENDDO              
            scalar=0
            DO ii=1,dim; scalar = scalar + sh_i(ii)*sh_j(ii); ENDDO
            IF(ABS(scalar)>=2)STOP'Scalar is too big'
              
           
            to_measure: IF(nmnm>2 .OR. vert1==link(vert2,1))THEN
              CALL VEC_BETWEEN(s2,s1,v_1,v_2,v)           ! space index 
             
              IF(scalar==0)CYCLE; !Skip if perpendicular
              
              ! factor due to removal of interaction line
              factor=un1/phonon_meas(s1,t1,s2,t2)
              ! sign factor
              factor=factor*Dphase_a
              ! factor due to flat statistics
              factor = factor / ord_probab(nmnm)
              ! factor due to inner product
              Z_TokTok_norm = Z_TokTok_norm + factor 
              ! factor due to derivative of thermodynamic potential
              factor =  factor 
              !multiply by scalar
              factor = factor * scalar 
             
              ! LISTING ALL 4 TOPOLOGICALLY DIFFERENT TYPES of VORTEX SHIFT
              ! ----- SHIFT OUTGOING LINES EACH VORTEX (1) -----------------------
              ! factor due to shift of ends in vert1: outgoing GF shifts 
              factor_out = un1
              !shift of the vert1 outgoing
              v3=link(vert1,2); s3=site(v3); t3=tau(v3);  
              factor_out = factor_out * GREEN(s1_sh,t1,s3,t3,sp_1) / GREEN(s1,t1,s3,t3,sp_1)
              !shift of the vert2 outgoing
              v3=link(vert2,2); s3=site(v3); t3=tau(v3);  
              factor_out = factor_out * GREEN(s2_sh,t2,s3,t3,sp_2) / GREEN(s2,t2,s3,t3,sp_2)
              IF(index_ri==1)THEN; ! ri=1 le=2
                rsit_outgo(1)=s1_sh ; rsit_ingo(1)=s1; 
                lsit_outgo(1)=s2_sh ; lsit_ingo(1)=s2; 
              ELSE                          ! le=1 ri=2
                rsit_outgo(1)=s2_sh ; rsit_ingo(1)=s2; 
                lsit_outgo(1)=s1_sh ; lsit_ingo(1)=s1; 
              ENDIF    
              ! -----  SHIFT OUTGOING LINES EACH VORTEX  ---------------------
              
              ! ----- SHIFT INGOING LINES EACH VORTEX (2)---------------------------
              ! factor due to shift of ends in vert1: outgoing GF shifts 
              factor_in = un1
              !shift of the vert1 ingoing
              v3=link(vert1,1); s3=site(v3); t3=tau(v3);  
              factor_in = factor_in * GREEN(s3,t3,s1_sh,t1,sp_in_1) / GREEN(s3,t3,s1,t1,sp_in_1)
              !shift of the vert2 outgoing
              v3=link(vert2,1); s3=site(v3); t3=tau(v3);  
              factor_in = factor_in * GREEN(s3,t3,s2_sh,t2,sp_in_2) / GREEN(s3,t3,s2,t2,sp_in_2)
              IF(index_ri==1)THEN; ! ri=1 le=2
                rsit_outgo(2)=s1 ; rsit_ingo(2)=s1_sh; 
                lsit_outgo(2)=s2 ; lsit_ingo(2)=s2_sh;
              ELSE;                        ! le=1 ri=2
                rsit_outgo(2)=s2 ; rsit_ingo(2)=s2_sh;
                lsit_outgo(2)=s1 ; lsit_ingo(2)=s1_sh;
              ENDIF      
              ! -----  SHIFT INTGOING LINES EACH VORTEX  ---------------------

              ! ----- SHIFT OF THE LINE GOING FROM vert1 TO vert2  (3)---------------------
              factor_1to2 = un1
              v1_to_2: IF(vert2==link(vert1,2))THEN; !vert1(s1->s1_sh)-->>--vert2(s2->s2_sh)
              factor_1to2 = factor_1to2 * GREEN(s1_sh,t1,s2_sh,t2,sp_1)/GREEN(s1,t1,s2,t2,sp_1)
              ELSE  v1_to_2 ! vert1(s1->s1_sh) --->>--- v3 -->>--...v4 -->>---- vert2(s2->s2_sh)
                v3=link(vert1,2); s3=site(v3); t3=tau(v3);  
                v4=link(vert2,1); s4=site(v4); t4=tau(v4)
                factor_1to2 = factor_1to2 / ( GREEN(s1,t1,s3,t3,sp_1) * GREEN(s4,t4,s2,t2,sp_in_2) )
                factor_1to2 = factor_1to2 * GREEN(s1_sh,t1,s3,t3,sp_1)*GREEN(s4,t4,s2_sh,t2,sp_in_2)
              ENDIF v1_to_2   
              IF(index_ri==1)THEN
                rsit_outgo(3)=s1_sh ; rsit_ingo(3)=s1; 
                lsit_outgo(3)=s2 ; lsit_ingo(3)=s2_sh;
              ELSE
                rsit_outgo(3)=s2 ; rsit_ingo(3)=s2_sh; 
                lsit_outgo(3)=s1_sh ; lsit_ingo(3)=s1;
              ENDIF      
              ! ----- SHIFT OF THE LINE GOING FROM vert1 TO vert2`  ----------------------------

              ! ----- SHIFT OF THE LINE GOING FROM vert2 TO vert1  (4)---------------------
              factor_2to1 = un1
              v2_to_1: IF(vert1==link(vert2,2))THEN  !ver2(s2->s2_sh) -->>--vert1(s1->s1_sh)
                factor_2to1 = factor_2to1 * GREEN(s2_sh,t2,s1_sh,t1,sp_2)/GREEN(s2,t2,s1,t1,sp_2)
              ELSE  v2_to_1   ! vert2(s2->s2_sh) --->>--- v4 -->>--...v3 -->>---- vert1(s1->s1_sh)
                v4=link(vert2,2); s4=site(v4); t4=tau(v4);  
                v3=link(vert1,1); s3=site(v3); t3=tau(v3)
                factor_2to1 = factor_2to1 / ( GREEN(s2,t2,s4,t4,sp_2) * GREEN(s3,t3,s1,t1,sp_in_1) )
                factor_2to1 = factor_2to1 *  GREEN(s2_sh,t2,s4,t4,sp_2)*GREEN(s3,t3,s1_sh,t1,sp_in_1)
              ENDIF v2_to_1    
              IF(index_ri==1)THEN
                rsit_outgo(4)=s1 ; rsit_ingo(4)=s1_sh; 
                lsit_outgo(4)=s2_sh ; lsit_ingo(4)=s2;
              ELSE
                rsit_outgo(4)=s2_sh ; rsit_ingo(4)=s2;
                lsit_outgo(4)=s1 ; lsit_ingo(4)=s1_sh;
              ENDIF    
              ! ----- SHIFT OF THE LINE GOING FROM vert1 TO vert2`  ----------------------------
              
              factor_all_include =  - factor*factor_out   - factor*factor_in   &
                                           + factor*factor_1to2  + factor*factor_2to1   
              
              ! Defining cell for given time for simple statistics 
              at=t1-t2; IF(at<0.d0) at=beta+at;  k=at/dtau   ! tau>0 index
             
              ! Simple bin statistics 
              TokTok_box(v(1),v(2),v(3),k) = TokTok_box(v(1),v(2),v(3),k) + factor_all_include
              ToTo_box_sep(k,1) = ToTo_box_sep(k,1) - factor*factor_out  
              ToTo_box_sep(k,2) = ToTo_box_sep(k,2) - factor*factor_in  
              ToTo_box_sep(k,3) = ToTo_box_sep(k,2) + factor*factor_1to2 
              ToTo_box_sep(k,4) = ToTo_box_sep(k,2) + factor*factor_2to1
              ! Basis function statistics 
              kb=at/dtb_OC; kb=min(kb,Nbins_OC-1)                 ! tau index basis    
              DO ib=1,Nbasis_OC      
                yc = BO_OC(ib,at,kb)*factor_all_include
                OCb(v(1),v(2),v(3),kb,ib) = OCb(v(1),v(2),v(3),kb,ib) + yc
              ENDDO 
              
! EXACT ESTiMATOR STATISTICS    
              ! Setting actusl maximal values
              case_max_actual = nul; jj_exa_proba(0:num_oc_lim) = nul; 
              
              measu_points: DO ipm=0,num_oc_lim !loop over measuring points
                  
                IF(ABS(ta_oc(ipm)-delta_tau)<window)THEN; !mesuring
                    
                  ! setting windows weight
                  IF( ta_oc(ipm) < window )THEN
                    oanwi = (2.0*window) / (window+ta_oc(ipm))
                  ELSE IF( ta_oc(ipm) > beta-window )THEN
                    oanwi = (2.0*window) / (window+beta-ta_oc(ipm))
                  ELSE
                    oanwi=un1
                  ENDIF    

                  ! makinf second reweighting   
                  de_de = ta_oc(ipm)-delta_tau
                  IF(ipm==num_oc_lim)de_de=de_de-1.0d-12
                  de_le = de_de * le_fra; 
                  de_ri = de_de * ri_fra;
                  IF(de_le<nul)de_le=de_le+1.0d-12 
                  
                  DO ik=1,4
                       
                      fabule = un1 
                      IF(rl_ingo==vert_le)THEN; ! going left --> right by fermion GF
                        fabule = fabule * GREEN(lsit_outgo(ik),tau_le-de_le,rsit_ingo(ik),tau_ri+de_ri,sp_r) / &
                                                  GREEN(lsit_outgo(ik),tau_le,rsit_ingo(ik),tau_ri,sp_r)
                      ELSE
                        fabule = fabule * GREEN(lsit_outgo(ik),tau_le-de_le,lsl_outgo,ltl_outgo,sp_l) / &         
                                                  GREEN(lsit_outgo(ik),tau_le,lsl_outgo,ltl_outgo,sp_l)
                        fabule = fabule * GREEN(rsl_ingo,rtl_ingo,rsit_ingo(ik),tau_ri+de_ri,sp_r) / &        
                                                  GREEN(rsl_ingo,rtl_ingo,rsit_ingo(ik),tau_ri,sp_r)
                      ENDIF  
                      IF(rl_outgo==vert_le)THEN; ! going right -> left by fermion GF  
                        fabule = fabule * GREEN(rsit_outgo(ik),tau_ri+de_ri,lsit_ingo(ik),tau_le-de_le,sp_r) / &         
                                                  GREEN(rsit_outgo(ik),tau_ri,lsit_ingo(ik),tau_le,sp_r) 
                      ELSE
                        fabule = fabule * GREEN(rsit_outgo(ik),tau_ri+de_ri,rsl_outgo,rtl_outgo,sp_r) / &         
                                                  GREEN(rsit_outgo(ik),tau_ri,rsl_outgo,rtl_outgo,sp_r)
                        fabule = fabule * GREEN(lsl_ingo,ltl_ingo,lsit_ingo(ik),tau_le-de_le,sp_l) / &         
                                                  GREEN(lsl_ingo,ltl_ingo,lsit_ingo(ik),tau_le,sp_l)
                      ENDIF    
                      SELECT CASE(ik)
                          CASE(1)
                              fabule = fabule * factor_out 
                          CASE(2)
                              fabule = fabule * factor_in
                          CASE(3)
                             fabule = - fabule * factor_1to2 
                          CASE(4)   
                             fabule = - fabule * factor_2to1
                          CASE DEFAULT
                              PRINT*,ik
                              STOP"NU I CASE!!!"
                      END SELECT    
                      buzota = - factor*fabule*oanwi
                      ! For OC contribution control
                      isensio = NINT( (buzota/aga_limit)*num_aga )
                      IF(isensio<=-num_aga)THEN;
                          aga_hist(-num_aga) = aga_hist(-num_aga) + 1.0d0
                      ELSE IF(isensio>=num_aga)THEN; 
                          aga_hist(num_aga) = aga_hist(num_aga) + 1.0d0
                      ELSE
                          aga_hist(isensio) = aga_hist(isensio) + 1.0d0
                      ENDIF    
                      ! End For OC contribution control
                      jj_exa(ipm) = jj_exa(ipm) + buzota
                      jj_exa_sep(ipm,ik) = jj_exa_sep(ipm,ik) + buzota 
                      jj_exa_proba(ipm) = jj_exa_proba(ipm) + buzota
                      !Handling switches for cut OC contributions
                       IF(ABS(buzota) > case_max_actual)THEN;
                           case_max_actual=ABS(buzota)
                           !PRINT*,"CMA:",case_max_actual
                       ENDIF    
                      !END Handling switches for cut OC contributions                      

                  ENDDO 
                
                ENDIF  
                
              ENDDO measu_points   

              DO ibi=0,skoko_cut; !PRINT*,case_max_actual,case_max_given(ibi)
                  IF(case_max_actual<case_max_given(ibi))THEN; 
                      stat_cut(ibi) = stat_cut(ibi) + un1;
                      DO ipm=0,num_oc_lim;
                          jj_exa_cut(ibi,ipm) = jj_exa_cut(ibi,ipm) + jj_exa_proba(ipm)
                      ENDDO;    
                  ENDIF; 
              ENDDO;    
              
! END EXACT ESTIMATOR STATISTICS
              
              
! EXACT ANOTHER ESTiMATOR STATISTICS              
              ! Setting actusl maximal values
              case_max_actual = nul; jj_exa_proba_ano(0:num_oc_lim_ano) = nul; 

              ano_oc: IF(MEASURE_OC_ano)THEN; !measuring another estimator 
              measu_points_ano: DO ipm=0,num_oc_lim_ano !loop over measuring points
                  
                IF(ipm>num_oc_lim_ano)THEN
                    PRINT*,ipm,num_oc_lim_ano; 
                    STOP
                ENDIF    
                  
                IF(ABS(ta_oc_ano(ipm)-delta_tau)<window)THEN; !mesuring
                    
                  ! setting windows weight
                  IF( ta_oc_ano(ipm) < window )THEN
                    oanwi = (2.0*window) / (window+ta_oc_ano(ipm))
                  ELSE IF( ta_oc_ano(ipm) > beta-window )THEN
                    oanwi = (2.0*window) / (window+beta-ta_oc_ano(ipm))
                  ELSE
                    oanwi=un1
                  ENDIF    

                  ! makinf second reweighting   
                  de_de = ta_oc_ano(ipm)-delta_tau
                  IF(ipm==num_oc_lim_ano)de_de=de_de-1.0d-12
                  de_le = de_de * le_fra; 
                  de_ri = de_de * ri_fra;
                  IF(de_le<nul)de_le=de_le+1.0d-12
                  !IF(ABS((tau_ri+de_ri)-(tau_le-de_le))>19.99d0)THEN
                  !   PRINT*,ipm,num_oc_lim_ano,ta_oc_ano(ipm);  
                  !   PRINT*,tau_le-de_le,tau_ri+de_ri
                  !   PRINT*,tau_le,tau_ri,de_le,de_ri
                  !ENDIF
                  
                  DO ik=1,4
                       
                      fabule = un1 
                      IF(rl_ingo==vert_le)THEN; ! going left --> right by fermion GF
                        fabule = fabule * GREEN(lsit_outgo(ik),tau_le-de_le,rsit_ingo(ik),tau_ri+de_ri,sp_r) / &
                                                  GREEN(lsit_outgo(ik),tau_le,rsit_ingo(ik),tau_ri,sp_r)
                      ELSE
                        fabule = fabule * GREEN(lsit_outgo(ik),tau_le-de_le,lsl_outgo,ltl_outgo,sp_l) / &         
                                                  GREEN(lsit_outgo(ik),tau_le,lsl_outgo,ltl_outgo,sp_l)
                        fabule = fabule * GREEN(rsl_ingo,rtl_ingo,rsit_ingo(ik),tau_ri+de_ri,sp_r) / &        
                                                  GREEN(rsl_ingo,rtl_ingo,rsit_ingo(ik),tau_ri,sp_r)
                      ENDIF  
                      IF(rl_outgo==vert_le)THEN; ! going right -> left by fermion GF  
                        fabule = fabule * GREEN(rsit_outgo(ik),tau_ri+de_ri,lsit_ingo(ik),tau_le-de_le,sp_r) / &         
                                                  GREEN(rsit_outgo(ik),tau_ri,lsit_ingo(ik),tau_le,sp_r) 
                      ELSE
                        fabule = fabule * GREEN(rsit_outgo(ik),tau_ri+de_ri,rsl_outgo,rtl_outgo,sp_r) / &         
                                                  GREEN(rsit_outgo(ik),tau_ri,rsl_outgo,rtl_outgo,sp_r)
                        fabule = fabule * GREEN(lsl_ingo,ltl_ingo,lsit_ingo(ik),tau_le-de_le,sp_l) / &         
                                                  GREEN(lsl_ingo,ltl_ingo,lsit_ingo(ik),tau_le,sp_l)
                      ENDIF    
                      SELECT CASE(ik)
                          CASE(1)
                              fabule = fabule * factor_out 
                          CASE(2)
                              fabule = fabule * factor_in
                          CASE(3)
                             fabule = - fabule * factor_1to2 
                          CASE(4)   
                             fabule = - fabule * factor_2to1
                          CASE DEFAULT
                              PRINT*,ik
                              STOP"NU I CASE!!!"
                      END SELECT    
                      buzota = - factor*fabule*oanwi
                      jj_exa_ano(ipm) = jj_exa_ano(ipm) + buzota
                      jj_exa_proba_ano(ipm) = jj_exa_proba_ano(ipm) + buzota
                      IF(ABS(buzota) > case_max_actual)THEN;
                         case_max_actual=ABS(buzota)
                      ENDIF    
                  ENDDO 
                
                ENDIF  
                
              ENDDO measu_points_ano   
              
              DO ibi=0,skoko_cut; 
                  IF(case_max_actual<case_max_given(ibi))THEN; 
                      DO ipm=0,num_oc_lim_ano;
                          jj_exa_cut_ano(ibi,ipm) = jj_exa_cut_ano(ibi,ipm) + jj_exa_proba_ano(ipm)
                      ENDDO;    
                  ENDIF; 
              ENDDO;    
              
              ENDIF ano_oc;
              
! END EXACT ANOTHER ESTIMATOR STATISTICS

            ENDIF to_measure   
          
          ENDDO dir_2
          ENDDO dir_1
! <jj> estimator end

      END SUBROUTINE MEASURING_OC
!-----------------------------------------------------------------------
                                       
!--------------------------------------------------------------------
!   Transforming: G0(p,t) --> G0(r,t) below
!--------------------------------------------------------------------
      SUBROUTINE GREEN0
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: ENER
      INTEGER :: backforth, ix,iy,iz,it
      DOUBLE PRECISION :: tt, mom(3), z, occ, occ1
      DOUBLE PRECISION, ALLOCATABLE :: gptr(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gpti(:,:,:,:)

!      print*, 'in GREEN0'
      allocate(gptr(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))
      allocate(gpti(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ntau0-1))

!      PRINT*,L(1)-1,L(2)-1,L(3)-1

      DO it=0,Ntau; tt=dtau*it                ! tau domain

      DO ix=0,L(1)-1; mom(1)=ix*(pi_m2/L(1))  ! momenta
      DO iy=0,L(2)-1; mom(2)=iy*(pi_m2/L(2))
      DO iz=0,L(3)-1; mom(3)=iz*(pi_m2/L(3))
      z=ENER(1,mom)                           ! dispersion relation
      IF(z<0.d0) THEN                         ! tau<0
        occ=1.d0/(dexp(z*beta)+1.d0)
        gptr(ix,iy,iz,0)=-occ*dexp(z*(beta-tt))
      ELSE                                    ! tau>0
        occ1=1.d0/(1.d0+dexp(-z*beta))
        gptr(ix,iy,iz,0)=-occ1*dexp(-z*tt)
      ENDIF
!      gptr(ix,iy,iz,0)= dexp(-z*tt)         ! single particle test
      gpti(ix,iy,iz,0)= 0.d0                 ! imaginary part is zero
      ENDDO; ENDDO; ENDDO;                   ! momentum function done

!      PRINT*,"here0"
 
      backforth=5                            ! do Fourier to go real space
      CALL fourierT(gptr,gpti,backforth,Ntau0)
!      PRINT*,"here1"
      gptr=gptr/Nsite                        ! in real space

!      PRINT*,"here2"

      DO ix=0,L(1)-1; DO iy=0,L(2)-1; DO iz=0,L(3)-1
         GRT_0(ix,iy,iz,it)=gptr(ix,iy,iz,0)  ! make a table bare for e-prop.
         GRT_NOW(ix,iy,iz,it)=GRT_0(ix,iy,iz,it)  ! make a table for e-prop.
      ENDDO; ENDDO; ENDDO

      ENDDO                                  ! tau domain

!      OPEN(1,file='grt.dat')
!      it=60; tt=dtau*it;
!      DO ix=0,L(1)-1
!        WRITE(1,*) ix, GRT_0(ix,0,0,it)
!      ENDDO
!      CLOSE(1)
!
!      OPEN(1,file='grt2.dat')
!      DO it=0,Ntau; tt=dtau*it;
!         WRITE(1,*) tt, GRT_0(0,0,0,it),GRT_0(1,0,0,it)
!      ENDDO
!      CLOSE(1)

      DEALLOCATE(gptr,gpti)

!      print*, 'GREEN0 DONE'
      END SUBROUTINE GREEN0
!........................................................

!--------------------------------------------------------------------
!   Checking density values from different sources
!--------------------------------------------------------------------
      SUBROUTINE DENSITY_CHECK
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: ENER, EX, GREEN, GREEN_000
      REAL*8,PARAMETER :: malo=1.0d-14
      REAL*8,DIMENSION(1:3) :: k_i
      REAL*8 :: n_green, n_green_0, n_kint, ene, pro0
      INTEGER :: i, j, color

      IF(.NOT. POLARIZED)THEN
      ! NON POL_SPIN------
      n_green = ABS(GREEN(1,malo,1,-malo,1 )) + & ! density for 2 spins
               ABS(GREEN(1,malo,1,-malo,-1))
      n_green_0 = ABS(GREEN_000(1,malo,1,-malo,1 ))  + & ! density for 2 spins
               ABS(GREEN_000(1,malo,1,-malo,-1))
      ! NON POL_SPIN------
      ELSE
      ! POL_SPIN------
      n_green = ABS(GREEN(1,malo,1,-malo,1 )) !+ & ! density for 2 spins
!               ABS(GREEN(1,malo,1,-malo,-1))
      n_green_0 = ABS(GREEN_000(1,malo,1,-malo,1 )) ! + &! density for 2 spins
!               ABS(GREEN_000(1,malo,1,-malo,-1))
      ! POL_SPIN------
      ENDIF
      
      dunsity_ST(i_all_meas)=n_green
      density=n_green

      IF(dim>2)STOP'This subroutine is for dim= 1 & 2 only. Rewrite it!'

      IF(dim==1)THEN
        n_kint=nul; k_i(1:3)=nul
        DO color = -1, 1, 2;                !Summing spin projections
          DO i = -Nsite/2, (Nsite/2)-1;     !Summing reciprocal vectors in BZ
            k_i(1) = (2*pi/Nsite)*i
            ene = ENER(color,k_i);
            n_kint = n_kint + un1 / (un1 + EX(ene*beta)) ;
          ENDDO;
        ENDDO;
        n_kint = n_kint / Nsite
      ELSE IF(dim==2)THEN
        n_kint=nul; k_i(1:3)=nul
        DO color = -1, 1, 2;                !Summing spin projections
          DO i = -L(1)/2, (L(1)/2)-1;     !Summing reciprocal vectors in BZ
          DO j = -L(2)/2, (L(2)/2)-1;     !Summing reciprocal vectors in BZ
            k_i(1) = (2*pi/L(1))*i
            k_i(2) = (2*pi/L(2))*j
            ene = ENER(color,k_i);
            n_kint = n_kint + un1 / (un1 + EX(ene*beta)) ;
          ENDDO;
          ENDDO;
        ENDDO;
        n_kint = n_kint / Nsite
      ELSE 
        STOP'This subroutine is for dim= 1 & 2 only. Rewrite it!'  
      ENDIF    
      
      PRINT"('n_gr = ',ES12.4,' n_gr_0 = ',ES12.4, ' n_gr/n_gr_0 = ',ES12.4)", &
                  n_green,n_green_0,n_green/n_green_0
      PRINT"('n_kinet =',ES12.4)",n_kint
      
      END SUBROUTINE DENSITY_CHECK
!........................................................

      
!-----------------------------------------------------------      
! Tabulating n(mu) for G_0 for input parameters
!------------------------------------------------------------
      SUBROUTINE N_OT_MU_TABU
      USE config_par; USE stat_par
      IMPLICIT NONE
      REAL*8,EXTERNAL :: GREEN_000, ENER, EX
      REAL*8,PARAMETER :: malo=1.0d-14
      REAL*8,DIMENSION(3) :: k_i
      REAL*8 :: ene
      INTEGER :: i, j, k, ika, color
      REAL*8 :: mu0
      LOGICAL:: comparing

! To comapre with G or not       
!      comparing = .TRUE.
      comparing = .FALSE.
      
! Saving current table for Green funtion and actual mu from input 
      GRT_NOW_SAVE = GRT_NOW; mu0 = mu

! Setting boundaries etc...
      R_min = -2.0d0 * (4.0d0*dim*hopping)
      R_plu =  4.0d0 * (4.0d0*dim*hopping)
      sha_demu = (R_plu-R_min)/i_demu
      
! Tabulating      
      DO ika=0,i_demu; n_ot_mu(ika)=nul
          mu_tabul(ika)=R_min+sha_demu*ika; mu=mu_tabul(ika)
         !Via GREEN_000 function - time consuming!!! 
          IF(comparing)THEN
            CALL GREEN0      !Preparing table for GREEN_000 for mu_tabul(i)
                             !also setting table GRT_NOW=GRT_0 for given mu_tabul(i)
            
            n_ot_mu_G(ika) = ABS(GREEN_000(1,malo,1,-malo, 1)) + &
                            ABS(GREEN_000(1,malo,1,-malo,-1))
            n_ot_mu_G(ika)=n_ot_mu_G(ika)*(pola_count/2.0d0)
          ENDIF  
         !Via unperturbed energies
          n_ot_mu(ika)=nul; k_i(1:3)=nul
          IF(dim==1)THEN;
            DO color = -1, 1, 2;              !Summing spin projections
            DO i = -Nsite/2, (Nsite/2)-1;     !Summing reciprocal in BZ
              k_i(1) = (2*pi/Nsite)*i; 
              ene = ENER(color,k_i);
              n_ot_mu(ika) = n_ot_mu(ika) + un1 / (un1 + EX(ene*beta)) ;
            ENDDO; ENDDO; 
          ELSE IF(dim==2)THEN; 
            DO color = -1, 1, 2;              !Summing spin projections
            DO i = -L(1)/2, (L(1)/2)-1;  !Summing reciprocal  in BZ
            DO j = -L(2)/2, (L(2)/2)-1;  !Summing reciprocal  in BZ
              k_i(1) = (2*pi/L(1))*i; k_i(2) = (2*pi/L(2))*j
              ene = ENER(color,k_i);
              n_ot_mu(ika) = n_ot_mu(ika) + un1 / (un1 + EX(ene*beta)) ;
            ENDDO; ENDDO; ENDDO;
          ELSE IF(dim==3)THEN; 
            DO color = -1, 1, 2;              !Summing spin projections
            DO i = -L(1)/2, (L(1)/2)-1;  !Summing reciprocal  in BZ
            DO j = -L(2)/2, (L(2)/2)-1;  !Summing reciprocal  in BZ
            DO k = -L(3)/2, (L(3)/2)-1;  !Summing reciprocal  in BZ
              k_i(1) = (2*pi/L(1))*i; k_i(2) = (2*pi/L(2))*j
              k_i(3) = (2*pi/L(3))*k; 
              ene = ENER(color,k_i);
              n_ot_mu(ika) = n_ot_mu(ika) + un1 / (un1 + EX(ene*beta)) ;
            ENDDO; ENDDO; ENDDO; ENDDO             
          ENDIF       
          n_ot_mu(ika) = n_ot_mu(ika) / Nsite
          n_ot_mu(ika) = n_ot_mu(ika) * (pola_count/2.0d0)
      ENDDO               
! Equalizibg two dimensions if not comapring
      IF(.NOT. comparing)n_ot_mu_G=n_ot_mu    
! Writing tabulation
      OPEN(UNIT=4,FILE="n_mu.dat")
      DO ika=0,i_demu
          WRITE(4,*)mu_tabul(ika),n_ot_mu(ika),n_ot_mu_G(ika)
      ENDDO    
      CLOSE(4)

! (i) Regaining current mu, 
!(ii) regaining table for Green_0 (GRT_0) at actual mu, 
!(iii)regaining table for current green funtion (GRT_NOW)
      mu = mu0; 
      CALL GREEN0;
      GRT_NOW = GRT_NOW_SAVE; 
      
      
      PRINT*,"  n(mu) tabulation completed"
      
      END SUBROUTINE N_OT_MU_TABU
!............................................................

                
!-----------------------------------------------------------      
! Tabulating n(mu) for G_0 for input parameters
!------------------------------------------------------------
      REAL*8 FUNCTION MU_OUTPUT
      USE config_par; USE stat_par
      IMPLICIT NONE;
      REAL*8,EXTERNAL :: GREEN_000, GREEN
      REAL*8,PARAMETER :: malo=1.0d-14
      REAL*8 :: n_input, derivat, n_che, mu0 
      INTEGER :: i, i_l, i_r
      
      n_input = ABS(GREEN(1,malo,1,-malo,1 )) + &! density for 2 spins
                     ABS(GREEN(1,malo,1,-malo,-1))
      
      IF(n_input<R_min)THEN; PRINT*,"Too small density"; STOP; ENDIF
      IF(n_input>R_plu)THEN; PRINT*,"Too large density"; 
      PRINT*,"R_PLU: ",R_plu; PRINT*,"n_input: ",n_input
         STOP; 
      ENDIF
          
      i_l=0; i_r=0;
      
      DO i=1,i_demu
          IF(n_input<n_ot_mu(i))THEN
              i_l=i-1; i_r=i; EXIT;
          ENDIF    
      ENDDO    
      IF(i_l==0 .AND. i_r==0)STOP'Nothing inside range'
      
      derivat = (n_input-n_ot_mu(i_l)) / (n_ot_mu(i_r)-n_ot_mu(i_l)) 
      mu_output = mu_tabul(i_l) + sha_demu*derivat

      PRINT*,"Requested mu: ",mu_tabul(i_l),mu_output,mu_tabul(i_r)

! Ckecking how the density is reporduced by the obtained mu_output 

! Saving current table for Green funtion and actual mu
      GRT_NOW_SAVE = GRT_NOW; mu0 = mu
      
      mu=mu_output
      CALL GREEN0      !Preparing table for GREEN_000 for mu_output
                       !also setting table GRT_NOW=GRT_0 for given mu_outpu
      n_che = ABS(GREEN_000(1,malo,1,-malo, 1)) + ABS(GREEN_000(1,malo,1,-malo,-1))

      PRINT*,"Densities: ",n_input,n_che,(n_input-n_che)/(n_input) 
      
! (i) Regaining current mu, 
!(ii) regaining table for Green_0 (GRT_0) at actual mu, 
!(iii)regaining table for current green funtion (GRT_NOW)
      mu = mu0; 
      CALL GREEN0;
      GRT_NOW = GRT_NOW_SAVE; 
      
      END FUNCTION MU_OUTPUT
!________________________________________________________
      
      
! adapting my 4D files to 1D fast Fourier transform
!________________________________________________________
      subroutine FourierT(TR,TI,backforth,Ltau)
      USE config_par; USE stat_par
      INTEGER, INTENT(IN) :: backforth, Ltau
      DOUBLE PRECISION, INTENT(OUT) :: TR(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ltau-1)
      DOUBLE PRECISION, INTENT(OUT) :: TI(0:L(1)-1,0:L(2)-1,0:L(3)-1,0:Ltau-1)
      integer :: ix,iy,iz,it
      integer :: power, Noma
      double precision, allocatable ::   Real1(:), Im1(:)

! OK, do Fourier transform for each coordinate one by one.

      power=log(Ltau*1.d0)/log(2.d0)+1.d-14
      Noma=2**power; allocate(Real1(0:Noma-1),Im1(0:Noma-1))
      do ix=0,L(1)-1; do iy=0,L(2)-1; do iz=0,L(3)-1;
      Real1(:)=TR(ix,iy,iz,:); Im1(:)=TI(ix,iy,iz,:)
      call sffteu(Real1, Im1, Noma, power, backforth )
      TR(ix,iy,iz,:)=Real1(:); TI(ix,iy,iz,:)=Im1(:)
      enddo; enddo; enddo                         ! time
      deallocate(Real1,Im1); ! print*, 't done'

      power=log(L(3)*1.d0)/log(2.d0)+1.d-14
      Noma=2**power; allocate(Real1(0:Noma-1),Im1(0:Noma-1))
      do ix=0,L(1)-1; do iy=0,L(2)-1; do it=0,Ltau-1
      Real1(:)=TR(ix,iy,:,it); Im1(:)=TI(ix,iy,:,it)
      call sffteu(Real1, Im1, Noma, power, backforth )
      TR(ix,iy,:,it)=Real1(:); TI(ix,iy,:,it)=Im1(:)
      enddo; enddo; enddo                         ! z
      deallocate(Real1,Im1); ! print*, 'z done'

      power=log(L(2)*1.d0)/log(2.d0)+1.d-14
      Noma=2**power; allocate(Real1(0:Noma-1),Im1(0:Noma-1))
      do ix=0,L(1)-1; do iz=0,L(3)-1; do it=0,Ltau-1
      Real1(:)=TR(ix,:,iz,it); Im1(:)=TI(ix,:,iz,it)
      call sffteu(Real1, Im1, Noma, power, backforth )
      TR(ix,:,iz,it)=Real1(:); TI(ix,:,iz,it)=Im1(:)
      enddo; enddo; enddo                         ! y
      deallocate(Real1,Im1); ! print*, 'y done'

      power=log(L(1)*1.d0)/log(2.d0)+1.d-14
      Noma=2**power; allocate(Real1(0:Noma-1),Im1(0:Noma-1))
      do iy=0,L(2)-1; do iz=0,L(3)-1; do it=0,Ltau-1
      Real1(:)=TR(:,iy,iz,it); Im1(:)=TI(:,iy,iz,it)
      call sffteu(Real1, Im1, Noma, power, backforth )
      TR(:,iy,iz,it)=Real1(:); TI(:,iy,iz,it)=Im1(:)
      enddo; enddo; enddo                         ! x
      deallocate(Real1,Im1); ! print*, 'x done'
!      print*, 'Full 4D transform done'

! define your function Real and Imaginary parts
!SR=0.d0; SI=0.d0;  ! in cont. space doamin S --> S dx (step in x grid)

! do the transform
      !call sffteu(SR, SI, Noma, power, 5 )
! now  SR and SI are Real and Imaginary parts of the Fourier transformed function

! if solving Dyson Equation use these functions as is; normalization is correct

!go back to time domain
      !call sffteu(SR, SI, Noma, power, -1)
! in cont. space doamin S --> S/dx (step in x grid)

      end subroutine FourierT


!______________________


!          managing basis functions for collecting statistics: sigma 
!______________________________________________________________           
           subroutine normalize_s
           USE config_par; IMPLICIT NONE
           DOUBLE PRECISION,EXTERNAL :: b_sigma,bo_sigma
           integer :: nz
           double precision :: z, dz , f, s, normb(Nbasis_s)
           integer :: i,j1,j2,k, in

           print*, 'IN NORMALIZATION'
           oc_s=0.d0
           do i=0,Nbins_s-1
 
           nz=1000; dz=(tbin_s(i+1)-tbin_s(i))/(nz*1.d0) 

           normb=0.d0
           do j2=1,Nbasis_s; oc_s(i,j2,j2)=1.d0; enddo


           do j2=1,Nbasis_s     ! start matrix contruction

           do j1=1,j2-1
           do in=1,nz
           z=(in-0.5d0)*dz + tbin_s(i)
           oc_s(i,j2,j1)=oc_s(i,j2,j1) -  (bo_sigma(j1,z,i)*b_sigma(j2,z,i)) * dz
           enddo
           enddo 

           do in=1,nz
           z=(in-0.5d0)*dz + tbin_s(i)
           f=0.d0; do j1=1,j2-1; f=f+oc_s(i,j2,j1)*bo_sigma(j1,z,i); 
           enddo
           f=b_sigma(j2,z,i)+f 
           normb(j2)=normb(j2) + f*f* dz
           enddo
           normb(j2)=dsqrt(normb(j2))  

           do j1=1,j2;   oc_s(i,j2,j1)=oc_s(i,j2,j1)/normb(j2); enddo
           do j1=1,j2-1; s=0.d0
           do k=j1,j2-1; s=s+oc_s(i,j2,k)*oc_s(i,k,j1); enddo
           oc_s(i,j2,j1)=s  
           enddo 

           enddo               ! end of matrix construction
  
! checking
           do j1=1,Nbasis_s; do j2=1,Nbasis_s
           s=0.d0
           do in=1,nz
           z=(in-0.5d0)*dz + tbin_s(i)
           s=s + bo_sigma(j1,z,i)*bo_sigma(j2,z,i) * dz
           enddo
!           print*, i,tbin(i),j1,j2 
!           print*, s
!           print*
            
           IF(j1==j2 .AND. ABS(s-1.d0)>1.d-2) stop ' normb ?'
           IF(j1.ne.j2 .AND. ABS(s)>1.d-2) stop ' ortho ?'
           enddo; enddo


           enddo


           print*, 'ortho-normalization done: sigma'
           end subroutine normalize_s
           
           

!          managing basis functions for collecting statistics: polarization 
!______________________________________________________________           
           subroutine normalize_p
           USE config_par; IMPLICIT NONE
           DOUBLE PRECISION,EXTERNAL :: b_polar,bo_polar
           integer :: nz
           double precision :: z, dz , f, s, normb(Nbasis_p)
           integer :: i,j1,j2,k, in

           print*, 'IN NORMALIZATION'
           oc_p=0.d0
           do i=0,Nbins_p-1
 
           nz=1000; dz=(tbin_p(i+1)-tbin_p(i))/(nz*1.d0) 

           normb=0.d0
           do j2=1,Nbasis_p; oc_p(i,j2,j2)=1.d0; enddo


           do j2=1,Nbasis_p     ! start matrix contruction

           do j1=1,j2-1
           do in=1,nz
           z=(in-0.5d0)*dz + tbin_p(i)
           oc_p(i,j2,j1)=oc_p(i,j2,j1) -  (bo_polar(j1,z,i)*b_polar(j2,z,i)) * dz
           enddo
           enddo 

           do in=1,nz
           z=(in-0.5d0)*dz + tbin_p(i)
           f=0.d0; do j1=1,j2-1; f=f+oc_p(i,j2,j1)*bo_polar(j1,z,i); 
           enddo
           f=b_polar(j2,z,i)+f 
           normb(j2)=normb(j2) + f*f* dz
           enddo
           normb(j2)=dsqrt(normb(j2))  

           do j1=1,j2;   oc_p(i,j2,j1)=oc_p(i,j2,j1)/normb(j2); enddo
           do j1=1,j2-1; s=0.d0
           do k=j1,j2-1; s=s+oc_p(i,j2,k)*oc_p(i,k,j1); enddo
           oc_p(i,j2,j1)=s  
           enddo 

           enddo               ! end of matrix construction
  
! checking
           do j1=1,Nbasis_p; do j2=1,Nbasis_p
           s=0.d0
           do in=1,nz
           z=(in-0.5d0)*dz + tbin_p(i)
           s=s + bo_polar(j1,z,i)*bo_polar(j2,z,i) * dz
           enddo
!           print*, i,tbin(i),j1,j2 
!           print*, s
!           print*
            
           IF(j1==j2 .AND. ABS(s-1.d0)>1.d-2) stop ' normb ?'
           IF(j1.ne.j2 .AND. ABS(s)>1.d-2) stop ' ortho ?'
           enddo; enddo


           enddo


           print*, 'ortho-normalization done: polar'
          end subroutine normalize_p
           

!          managing basis functions for collecting statistics: OC 
!______________________________________________________________           
           subroutine normalize_OC
           USE config_par; IMPLICIT NONE
           DOUBLE PRECISION,EXTERNAL :: b_OC,bo_OC
           integer :: nz
           double precision :: z, dz , f, s, normb(Nbasis_OC)
           integer :: i,j1,j2,k, in

           print*, 'IN NORMALIZATION OC'
           oc_OC=0.d0
           do i=0,Nbins_OC-1
 
           nz=1000; dz=(tbin_OC(i+1)-tbin_OC(i))/(nz*1.d0) 

           normb=0.d0
           do j2=1,Nbasis_OC; oc_OC(i,j2,j2)=1.d0; enddo


           do j2=1,Nbasis_OC     ! start matrix contruction

           do j1=1,j2-1
           do in=1,nz
           z=(in-0.5d0)*dz + tbin_OC(i)
           oc_OC(i,j2,j1)=oc_OC(i,j2,j1) - (bo_OC(j1,z,i)*b_OC(j2,z,i)) * dz
           enddo
           enddo 

           do in=1,nz
           z=(in-0.5d0)*dz + tbin_OC(i)
           f=0.d0; do j1=1,j2-1; f=f+oc_OC(i,j2,j1)*bo_OC(j1,z,i); 
           enddo
           f=b_OC(j2,z,i)+f 
           normb(j2)=normb(j2) + f*f* dz
           enddo
           normb(j2)=dsqrt(normb(j2))  

           do j1=1,j2;   oc_OC(i,j2,j1)=oc_OC(i,j2,j1)/normb(j2); enddo
           do j1=1,j2-1; s=0.d0
           do k=j1,j2-1; s=s+oc_OC(i,j2,k)*oc_OC(i,k,j1); enddo
           oc_OC(i,j2,j1)=s  
           enddo 

           enddo               ! end of matrix construction
  
! checking
           do j1=1,Nbasis_OC; do j2=1,Nbasis_OC
           s=0.d0
           do in=1,nz
           z=(in-0.5d0)*dz + tbin_OC(i)
           s=s + bo_OC(j1,z,i)*bo_OC(j2,z,i) * dz
           enddo
!           print*, i,tbin(i),j1,j2 
!           print*, s
!           print*
            
           IF(j1==j2 .AND. ABS(s-1.d0)>1.d-2) stop ' normb ?'
           IF(j1.ne.j2 .AND. ABS(s)>1.d-2) stop ' ortho ?'
           enddo; enddo


           enddo


           print*, 'ortho-normalization done: OC'
           end subroutine normalize_OC
          
          
!_______________________________________________BASIS for sigma
           double precision function b_sigma(jb,x,kb)
           USE config_par; IMPLICIT NONE
           integer, intent (IN) :: jb, kb
           double precision, intent (IN) :: x
           double precision :: xx, zx


           IF(kb==Nbins_s) then; print*, jb, kb; stop 'Nbins b'; endif 
             xx=x-tce_s(kb); 
             b_sigma=xx**(jb-1)          ! power law basis
                
           end function b_sigma
           
!_______________________________________________BASIS for polar
           double precision function b_polar(jb,x,kb)
           USE config_par; IMPLICIT NONE
           integer, intent (IN) :: jb, kb
           double precision, intent (IN) :: x
           double precision :: xx, zx


           IF(kb==Nbins_p) then; print*, jb, kb; stop 'Nbins b'; endif 
             xx=x-tce_p(kb); 
             b_polar=xx**(jb-1)          ! power law basis
                
           end function b_polar
           
!_______________________________________________BASIS for OC
           double precision function b_OC(jb,x,kb)
           USE config_par; IMPLICIT NONE
           integer, intent (IN) :: jb, kb
           double precision, intent (IN) :: x
           double precision :: xx, zx


           IF(kb==Nbins_OC) then; print*, jb, kb; stop 'Nbins OC';endif 
             xx=x-tce_OC(kb); 
             b_OC=xx**(jb-1)          ! power law basis
                
           end function b_OC
           
!_______________________________________________ORTHO BASIS for sigma
           double precision function bo_sigma(jb,x,kb)
           USE config_par; IMPLICIT NONE
           DOUBLE PRECISION,EXTERNAL :: b_sigma
           integer, intent (IN) :: jb, kb
           integer :: jba
           double precision, intent (IN) :: x
           
           IF(kb==Nbins_s) then; print*, jb, kb; stop 'Nbins bo si'; endif 
           bo_sigma=0.d0
           do jba=1,Nbasis_s; 
               bo_sigma=bo_sigma+oc_s(kb,jb,jba)*b_sigma(jba,x,kb); 
           enddo

           end function bo_sigma
!________________________________________________
  
!_______________________________________________ORTHO BASIS for polar
           double precision function bo_polar(jb,x,kb)
           USE config_par; IMPLICIT NONE
           DOUBLE PRECISION,EXTERNAL :: b_polar
           integer, intent (IN) :: jb, kb
           integer :: jba
           double precision, intent (IN) :: x

           IF(kb==Nbins_p) then; print*, jb, kb; stop 'Nbins bo pol'; endif 
           bo_polar=0.d0
           do jba=1,Nbasis_p; 
               bo_polar=bo_polar+oc_p(kb,jb,jba)*b_polar(jba,x,kb); 
           enddo

           end function bo_polar
!________________________________________________
  
!_______________________________________________ORTHO BASIS for OC           
           double precision function bo_OC(jb,x,kb)
           USE config_par; IMPLICIT NONE
           DOUBLE PRECISION,EXTERNAL :: b_OC
           integer, intent (IN) :: jb, kb
           integer :: jba
           double precision, intent (IN) :: x

           IF(kb==Nbins_OC) then; print*, jb, kb; stop 'Nb bo OC'; endif
           bo_OC=0.d0
           do jba=1,Nbasis_OC; 
               bo_OC=bo_OC+oc_OC(kb,jb,jba)*b_OC(jba,x,kb); 
           enddo

           end function bo_OC
!________________________________________________
           
!--------------------------------------------------------------------
! Sigma from the basis functions
!--------------------------------------------------------------------
      REAL*8 FUNCTION FUN_sigma(k1,k2,k3,xx)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      DOUBLE PRECISION,EXTERNAL :: bo_sigma
      INTEGER,INTENT(IN) :: k1,k2,k3
      REAL*8,INTENT(IN) :: xx
      INTEGER :: i_bin, i
      
      i_bin=-1; 
      
      DO i = 0,Nbins_s-1
          IF(xx<tbin_s(i+1))THEN
              i_bin=i; EXIT;
          ENDIF    
      ENDDO    
      
      IF(i_bin==-1 .OR. i_bin>Nbins_s-1)THEN
          PRINT*,i_bin; STOP"No such ibins in sigma"
      ENDIF    
      
      fun_sigma=0.0d0;
      DO i=1,Nbasis_s
        fun_sigma = fun_sigma + Sigmab(k1,k2,k3,i_bin,i) * bo_sigma(i,xx,i_bin)
      ENDDO    
      END FUNCTION Fun_sigma
!...................................................................

      
!--------------------------------------------------------------------
! POlar from the basis functions
!--------------------------------------------------------------------
      REAL*8 FUNCTION FUN_polar(k1,k2,k3,xx)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      DOUBLE PRECISION,EXTERNAL :: bo_polar
      INTEGER,INTENT(IN) :: k1,k2,k3
      REAL*8,INTENT(IN) :: xx
      INTEGER :: i_bin, i
      
      i_bin=-1; 
      
      DO i = 0,Nbins_p-1
          IF(xx<tbin_p(i+1))THEN
              i_bin=i; EXIT;
          ENDIF    
      ENDDO    
      
      IF(i_bin==-1 .OR. i_bin>Nbins_p-1)THEN
          PRINT*,i_bin; STOP"No such ibins in polar"
      ENDIF    
      
      fun_polar=0.0d0;
      DO i=1,Nbasis_p
        fun_polar = fun_polar + Polarb(k1,k2,k3,i_bin,i) * bo_polar(i,xx,i_bin)
      ENDDO    
      END FUNCTION Fun_polar
!...................................................................

!--------------------------------------------------------------------
! OC from the basis functions
!--------------------------------------------------------------------
      REAL*8 FUNCTION FUN_OC(k1,k2,k3,xx)
      USE config_par; USE stat_par;
      IMPLICIT NONE
      DOUBLE PRECISION,EXTERNAL :: bo_OC
      INTEGER,INTENT(IN) :: k1,k2,k3
      REAL*8,INTENT(IN) :: xx
      INTEGER :: i_bin, i
      
      i_bin=-1; 
      
      DO i = 0,Nbins_OC-1
          IF(xx<tbin_OC(i+1))THEN
              i_bin=i; EXIT;
          ENDIF    
      ENDDO    
      
      IF(i_bin==-1 .OR. i_bin>Nbins_OC-1)THEN
          PRINT*,i_bin; STOP"No such ibins in polar"
      ENDIF    
      
      fun_OC=0.0d0;
      DO i=1,Nbasis_OC
        fun_OC = fun_OC + OCb(k1,k2,k3,i_bin,i) * bo_OC(i,xx,i_bin)
      ENDDO    
      END FUNCTION Fun_OC
!...................................................................
