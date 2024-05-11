C*********************************************************************
C                                                           
C     Implementation of Poroflow for Multiphysics-LDPM (M-LDPM), 
C     a.k.a. Flow Lattice Method (FLM) For Pore Pressure Flow 
C     via User Defined Element in Abaqus
C     Hao Yin   
C     Northwestern University
C     Email: haoyin2022@u.northwestern.edu                            
C                                                           
C     Last edit: 01/23/2024
C
C*********************************************************************

      MODULE ModuleFLE
	  ! Global variables storage for external files and named pipes (UEXTERNALDB)
      IMPLICIT NONE
	  
	  ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
	  
	  INTEGER  :: nelem,nnode_FLM,nelem_LDPM,MPtypeFLM,nExchangeTime
	  REAL(dp) :: period_FLM,period_LDPM,time_LDPM,time_LDPM_old
	  
	  ! Edge element geometry variables
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: TET_DATA
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ELEM_DATA
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tetPts1
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tetPts2
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: normal 
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: segLen
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: segVol	  
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sharedArea
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: LF
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: AF
	  INTEGER, DIMENSION(:,:),  ALLOCATABLE :: tetIndex
	  INTEGER, DIMENSION(:,:),  ALLOCATABLE :: LDPMtetIndex
	  INTEGER, DIMENSION(:,:),  ALLOCATABLE :: LDPMtet1FacetIndex
	  INTEGER, DIMENSION(:,:),  ALLOCATABLE :: LDPMtet2FacetIndex
	  INTEGER, DIMENSION(:),    ALLOCATABLE :: elemflag
	  
	  ! Two-way coupling variables
	  INTEGER  :: kinc_in_Lop1,exchangeindex,ExchangeFlag
	  REAL(dp) :: CurrentExchangeTime
	  
	  REAL(dp), DIMENSION(:),   ALLOCATABLE :: ExchangeTimeList
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: FLM2LDPM_DATA
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: LDPM2FLM_DATA
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: LDPM2FLM_DATA_old
	  
	  CHARACTER(LEN = 256) :: OUTDIR
      CHARACTER(LEN = 256) :: JOBNAME
	  CHARACTER(LEN = 1)   :: path_separator
	  CHARACTER(LEN = 256) :: U2Vpipefile
	  CHARACTER(LEN = 256) :: V2Upipefile
	  
	  INTEGER, PARAMETER :: ndime = 3
	  INTEGER, PARAMETER :: V2U = 905 ! file unit number for named pipes
	  INTEGER, PARAMETER :: U2V = 906 ! file unit number for named pipes

      SAVE
      END MODULE ModuleFLE


      MODULE GlobalVarStorage  
      ! Global variables storage for visualization
	  
	  USE ModuleFLE
      IMPLICIT NONE
	  
	  ! Variables for visualization
      INTEGER :: n_stored_vars ! This number equals to the number below *Depvar in INP file
      INTEGER, PARAMETER :: n_int_pnts = 2
	  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: UserVar

	  ! Variables for inter-element calculations
	  REAL(dp)  :: vol_total
		
      ! Element counting variables
      INTEGER :: nelemUel

      END MODULE GlobalVarStorage
	  

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,PROPS,         
	1 	NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,KSTEP,KINC,
	2   JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,
	3   MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

	  USE GlobalVarStorage
	  USE ModuleFLE
	  
      INCLUDE 'aba_param.inc' 

C     VARIABLES THAT VARIABLE TYPES NEEDED TO BE CHANGED
! 	  SINCE implicit real*8(a-h,o-z) IS DECLARED IN 'aba_param.inc'
! 	  ALL VARIABLES WITH NAMES BEGINNING WITH A-H AND O-Z ARE DECLARED
!     real*8, I-N ARE DECLARED integer

C     ABAQUS UEL VARIABLES	  
	  INTEGER :: NDOFEL,NRHS,NSVARS,MCRD,NNODE,JTYPE,NPROPS,
	1  			 KSTEP,KINC,NDLOAD,NPREDF,MLVARX,MDLOAD,
	2  			 NJPROP,JELEM
	  INTEGER :: LFLAGS(3),JDLTYP(MDLOAD,*),JPROPS(NJPROP)
	  
	  REAL*8 :: DTIME,PNEWDT,PERIOD
	  REAL*8 :: RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(NSVARS),
	1 			ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
	2 			U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
	3 			PARAMS(3),ADLMAG(MDLOAD,*),
	4 			DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE)
	
C     PROPERTY VARIABLES 
	  REAL(dp) :: rho0,mu,Kf,kappa0,Mb,b,p0,timeScalingFac,fractureFlag
	  
C     FLM GEOMETRY VARIABLES 
	  REAL(dp) :: tetVec(NDIME)
      REAL(dp), DIMENSION(:), ALLOCATABLE :: AW
      REAL(dp), DIMENSION(:), ALLOCATABLE :: VW
      REAL(dp), DIMENSION(:), ALLOCATABLE :: LENGTH
      REAL(dp), DIMENSION(:), ALLOCATABLE :: G1
	  REAL(dp), DIMENSION(:), ALLOCATABLE :: G2
	  REAL(dp), DIMENSION(:), ALLOCATABLE :: EpsilonV1
	  REAL(dp), DIMENSION(:), ALLOCATABLE :: EpsilonV1_old
	  REAL(dp), DIMENSION(:), ALLOCATABLE :: EpsilonV2
	  REAL(dp), DIMENSION(:), ALLOCATABLE :: EpsilonV2_old
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: eVec
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: VC
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: VC_old
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: deltaN1
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: deltaN1_old
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: deltaN2
	  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: deltaN2_old
	  
C     INTERNAL VARIABLES
      REAL(dp) :: P1,P1OLD,Delta_P1,DP1DT,P2,P2OLD,Delta_P2,DP2DT,
	1			EpsilonV1NEW,EpsilonV1OLD,EpsilonV2NEW,EpsilonV2OLD,
	2			Delta_EpsilonV1,Delta_EpsilonV2,
	3			dEpsilonV1_dtime,dEpsilonV2_dtime,Vc1,Vc1OLD,Vc2,Vc2OLD,
	4			Delta_Vc1,Delta_Vc2,dVC1_dtime,dVC2_dtime,
	5			CBAR,C1,C2,XI,SBAR,S1,S2

C     COUPLING VARIABLES
      REAL(dp) :: scaledtime_LDPM,scaledtime_old_LDPM,DTIME_coupling
	  
	  INTEGER, PARAMETER :: ZERO=0,ONE=1
	  REAL(dp), PARAMETER :: DZERO=0.D0,DONE=1.D0,HALF=0.5D0

C*********************************************************************
C
C    TWO NODES EDGE FLOW LATTICE UEL FOR PORE PRESSURE DIFFUSION ANALYSIS
C
C*********************************************************************
C 	 SPATIAL DISCRETIZATION
C
C	 THE SPATIAL DISCRETIZATION FOR FLOW LATTICE ELEMENT IS A DISCRETE FORMULATION
C    PROPORTIONALITIES G1 and G2 ARE USED FOR THE INTERPOLATION OF NODAL VARIABLES 
C    OVER THE ELEMENT DOMAIN (FOR EXAMPLE X = G2*X(T1) + G1*X(T2))
C
C*********************************************************************
C 	 TEMPORAL DISCRETIZATION
C
C    ABAQUS/STANDARD USES BACKWARD EULER METHOD FOR TIME DISCRETIZATION
C    INFO IS EVALUATED AT TIME T + DELTA T
C
C*********************************************************************
C     VARIABLE DECLARATIONS
C
C     VW    :  UNCRACKED MATERIAL VOLUME (VW=VW1+VW2)
C     SEGVOL:  UNCRACKED MATERIAL VOLUME OF SIDE T1 AND T2 (VW1,VW2)
C	  		   SEGVOL = [VW1_1 VW1_2 ... VW1_NELEM
C						 VW2_1 VW2_2 ... VW2_NELEM]
C     VC    :  CRACKED MATERIAL VOLUME OF SIDE T1 and T2 (VC1,VC2, VC=AF*deltaN)
C     AF    :  AREA OF THE LDPM FACET (AF11,AF12,AF13,AF21,AF22,AF23)
C              AF = [AF11_1 AF11_2 ... AF11_NELEM
C                    AF12_1 AF12_2 ... AF12_NELEM
C                    AF13_1 AF13_2 ... AF13_NELEM
C                    AF21_1 AF21_2 ... AF21_NELEM
C                    AF22_1 AF22_2 ... AF22_NELEM
C                    AF23_1 AF23_2 ... AF23_NELEM]
C     LF    :  LENGTH OF INTERSECTION OF LDPM FACET WITH THE TETRAHEDRON FACE P1P2P3
C              LF = [LF1_1 LF1_2 ... LF1_NELEM
C                    LF2_1 LF2_2 ... LF2_NELEM
C                    LF3_1 LF3_2 ... LF3_NELEM]
C     EpsilonV1:     VOLUMETRIC STRAIN OF SIDE T1
C     EpsilonV2:     VOLUMETRIC STRAIN OF SIDE T2
C     deltaN1:  CRACK OPENING ASSOCIATED WITH LF ON THE SIDE T1 (LDPM GEOMETRY)
C              deltaN1 = [deltaN11_1 deltaN11_2 ... deltaN11_NELEM
C                         deltaN12_1 deltaN12_2 ... deltaN12_NELEM
C                         deltaN13_1 deltaN13_2 ... deltaN13_NELEM]
C     deltaN2:  CRACK OPENING ASSOCIATED WITH LF ON THE SIDE T2 (LDPM GEOMETRY)
C              deltaN2 = [deltaN21_1 deltaN21_2 ... deltaN21_NELEM
C                         deltaN22_1 deltaN22_2 ... deltaN22_NELEM
C                         deltaN23_1 deltaN23_2 ... deltaN23_NELEM]
C     SHAREDAREA:  AREA OF THE TETRAHEDRON FACE P1P2P3 (LDPM GEOMETRY)
C     NRML  :  NORMAL VECTOR OF FACE P1P2P3 (a)
C     EVEC  :  DIRECTION VECTOR FROM T2 TO T1 (E)
C     AW    :  PROJECTION AREA OF P1P2P3 IN E DIRECTION (A_W=sharedArea*DOT(a,E))
C     LENGTH:  LENGTH OF FLE
C     SEGLEN:  LENGTH OF SIDE T1 AND T2 (L1,L2)
C	  		   SEGLEN = [L1_1 L1_2 ... L1_NELEM
C						 L2_1 L2_2 ... L2_NELEM]
C     G1    :  LENGTH PROPORTIONALITY COEFFICIENT OF SIDE T1 (G1 = L1/LENGTH)
C     G2    :  LENGTH PROPORTIONALITY COEFFICIENT OF SIDE T2 (G2 = L2/LENGTH)
C
C     CBAR  :  WATER CAPACITY OR STORAGE COEFFICIENT (CBAR = G2*C1 + G1*C2)            
C     XI    :  WATER PERMEABILITY COEFFICIENT
C     SBAR  :  SOURCE TERM, PRESSURE RELEASE DUE TO CRACK BEHAVIOR (SBAR = G2*S1 + G1*S2)        
C     U     :  DOFS AT TIME T + DELTA T (WITH FORMAT U = [P1 P2])
C     P1    :  PRESSURE AT TIME T + DELTA T AT TET POINT 1
C     P2    :  PRESSURE AT TIME T + DELTA T AT TET POINT 2
C     P1OLD :  PRESSURE AT TIME T AT TET POINT 1
C     P2OLD :  PRESSURE AT TIME T AT TET POINT 2
C     DP1DT :  DERIVATIVE OF PRESSURE WRT TIME AT TET POINT 1
C     DP2DT :  DERIVATIVE OF PRESSURE WRT TIME AT TET POINT 2
C
C*********************************************************************
C     MATERIAL PROPERTY DEFINITION
C
C	  REAL PROPERTY VALUES
	  rho0 = props(1) ! Reference fluid density
	  mu = props(2) ! Dynamic viscosity of fluid
	  Kf = props(3) ! Fluid bulk modulus
	  kappa0 = props(4) ! Intrinsic permeability of the porous media
	  Mb = props(5) ! Biot modulus
	  b = props(6) ! Biot coefficient
	  p0 = props(7) ! Reference/initial pore pressure
	  timeScalingFac = props(8) ! Time Scaling Factor
	  fractureFlag = props(9) ! Fracture Analysis Flag

C	  INTEGER PROPERTY VALUES
C
	
C*********************************************************************
C 	  READ GEOMETRY DATA
C
	  IF (allocated(tetPts1) == 0) THEN
		ALLOCATE(tetPts1(3,NELEM),tetPts2(3,NELEM),tetIndex(2,NELEM),
	1			 normal(3,NELEM),segLen(NNODE,NELEM),segVol(NNODE,NELEM),
	2			 LF(3,NELEM),AF(6,NELEM),sharedArea(1,NELEM),LDPMtetIndex(2,NELEM),
	3		     LDPMtet1FacetIndex(3,NELEM),LDPMtet2FacetIndex(3,NELEM),elemflag(NELEM))
	  ENDIF

	  IF (allocated(eVec) == 0) THEN
		ALLOCATE(eVec(NDIME,NELEM),AW(NELEM),VW(NELEM),
	1			 LENGTH(NELEM),G1(NELEM),G2(NELEM),VC(2,NELEM),VC_old(2,NELEM),
	2			 EpsilonV1(NELEM),EpsilonV1_old(NELEM),EpsilonV2(NELEM),EpsilonV2_old(NELEM),
	3			 deltaN1(3,NELEM),deltaN1_old(3,NELEM),deltaN2(3,NELEM),deltaN2_old(3,NELEM))
	  ENDIF
		
	  IF (JELEM == 1) THEN
		DO i=1,NELEM
			tetIndex(1,i) = ELEM_DATA(i,1)
			tetIndex(2,i) = ELEM_DATA(i,2)
			sharedArea(1,i) = ELEM_DATA(i,3)
			normal(1,i) = ELEM_DATA(i,4)
			normal(2,i) = ELEM_DATA(i,5)
			normal(3,i) = ELEM_DATA(i,6)
			segLen(1,i) = ELEM_DATA(i,7)
			segLen(2,i) = ELEM_DATA(i,8)
			segVol(1,i) = ELEM_DATA(i,9)
			segVol(2,i) = ELEM_DATA(i,10)
			LF(1,i) = ELEM_DATA(i,11)
			LF(2,i) = ELEM_DATA(i,12)
			LF(3,i) = ELEM_DATA(i,13)
			AF(1,i) = ELEM_DATA(i,14)
			AF(2,i) = ELEM_DATA(i,15)
			AF(3,i) = ELEM_DATA(i,16)
			AF(4,i) = ELEM_DATA(i,17)
			AF(5,i) = ELEM_DATA(i,18)
			AF(6,i) = ELEM_DATA(i,19)
			LDPMtetIndex(1,i) = ELEM_DATA(i,20)
			LDPMtetIndex(2,i) = ELEM_DATA(i,21)
			LDPMtet1FacetIndex(1,i) = ELEM_DATA(i,22)
			LDPMtet1FacetIndex(2,i) = ELEM_DATA(i,23)
			LDPMtet1FacetIndex(3,i) = ELEM_DATA(i,24)
			LDPMtet2FacetIndex(1,i) = ELEM_DATA(i,25)
			LDPMtet2FacetIndex(2,i) = ELEM_DATA(i,26)
			LDPMtet2FacetIndex(3,i) = ELEM_DATA(i,27)
			elemFlag(i) = ELEM_DATA(i,28)
		END DO
			
		DO i=1,NELEM
			tetPts1(1,i) = TET_DATA(tetIndex(1,i),1)
			tetPts1(2,i) = TET_DATA(tetIndex(1,i),2)
			tetPts1(3,i) = TET_DATA(tetIndex(1,i),3)
			tetPts2(1,i) = TET_DATA(tetIndex(2,i),1)
			tetPts2(2,i) = TET_DATA(tetIndex(2,i),2)
			tetPts2(3,i) = TET_DATA(tetIndex(2,i),3)
		END DO
	  ENDIF

C	  CALCULATE UNCRACKED MATERIAL VOLUME (VW)
	  VW(JELEM) = SEGVOL(1,JELEM) + SEGVOL(2,JELEM)

C     CALCULATE DIRECTION VECTOR FROM T2 TO T1 (E)
      tetVec(1) = tetPts1(1,JELEM)- tetPts2(1,JELEM)
      tetVec(2) = tetPts1(2,JELEM)- tetPts2(2,JELEM)
      tetVec(3) = tetPts1(3,JELEM)- tetPts2(3,JELEM)
      tetDist = sqrt(tetVec(1)*tetVec(1)+tetVec(2)*tetVec(2)+tetVec(3)*tetVec(3))
      eVec(1,JELEM) = tetVec(1)/tetDist
      eVec(2,JELEM) = tetVec(2)/tetDist
      eVec(3,JELEM) = tetVec(3)/tetDist

C	  CALCULATE PROJECTION AREA OF P1P2P3 IN E DIRECTION (A_W)
      AW(JELEM) = abs(sharedArea(1,JELEM)*dot_product(normal(:,JELEM),eVec(:,JELEM)))

C	  CALCULATE ELEMENT LENGTH (LENGTH)
      LENGTH(JELEM) = segLen(1,JELEM) + segLen(2,JELEM)

C	  CALCULATE LENGTH PROPORTIONALITY OF SIDES (G1,G2)
      IF (elemflag(JELEM) == 3) THEN
		G1(JELEM) = HALF
		G2(JELEM) = HALF
      ELSE
		G1(JELEM) = segLen(1,JELEM)/LENGTH(JELEM)
		G2(JELEM) = segLen(2,JELEM)/LENGTH(JELEM)
      ENDIF

C*********************************************************************
C     INITIALIZATION  (NRHS=1)
C
      AMATRX = DZERO
      DO K1=1,NDOFEL                      
		RHS(K1,NRHS) = DZERO
      END DO
		
      period_FLM = period

C*********************************************************************
C     START THE ANALYSIS
C
C	  SKIP THE MASS MATRIX INITIALIZATION
      IF (LFLAGS(3) == 4) RETURN

C 	  TWO-WAY COUPLING PROCESS
	  IF (ExchangeFlag == 1) THEN
		IF (elemflag(JELEM) == 1) THEN
			EpsilonV1(JELEM) = LDPM2FLM_DATA(LDPMtetIndex(1,JELEM),1)
			EpsilonV2(JELEM) = LDPM2FLM_DATA(LDPMtetIndex(2,JELEM),1)
				
			EpsilonV1_old(JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(1,JELEM),1)
			EpsilonV2_old(JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(2,JELEM),1)
			
			deltaN1(1,JELEM) = LDPM2FLM_DATA(LDPMtetIndex(1,JELEM),LDPMtet1FacetIndex(1,JELEM))
			deltaN1(2,JELEM) = LDPM2FLM_DATA(LDPMtetIndex(1,JELEM),LDPMtet1FacetIndex(2,JELEM))
			deltaN1(3,JELEM) = LDPM2FLM_DATA(LDPMtetIndex(1,JELEM),LDPMtet1FacetIndex(3,JELEM))
			deltaN2(1,JELEM) = LDPM2FLM_DATA(LDPMtetIndex(2,JELEM),LDPMtet2FacetIndex(1,JELEM))
			deltaN2(2,JELEM) = LDPM2FLM_DATA(LDPMtetIndex(2,JELEM),LDPMtet2FacetIndex(2,JELEM))
			deltaN2(3,JELEM) = LDPM2FLM_DATA(LDPMtetIndex(2,JELEM),LDPMtet2FacetIndex(3,JELEM))

			deltaN1_old(1,JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(1,JELEM),LDPMtet1FacetIndex(1,JELEM))
			deltaN1_old(2,JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(1,JELEM),LDPMtet1FacetIndex(2,JELEM))
			deltaN1_old(3,JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(1,JELEM),LDPMtet1FacetIndex(3,JELEM))
			deltaN2_old(1,JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(2,JELEM),LDPMtet2FacetIndex(1,JELEM))
			deltaN2_old(2,JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(2,JELEM),LDPMtet2FacetIndex(2,JELEM))
			deltaN2_old(3,JELEM) = LDPM2FLM_DATA_old(LDPMtetIndex(2,JELEM),LDPMtet2FacetIndex(3,JELEM))
		ELSE IF (elemflag(JELEM) == 2) THEN
			EpsilonV1(JELEM) = DZERO
			EpsilonV2(JELEM) = DZERO
			EpsilonV1_old(JELEM) = DZERO
			EpsilonV2_old(JELEM) = DZERO
			
			deltaN1(1:3,JELEM) = DZERO
			deltaN2(1:3,JELEM) = DZERO
			deltaN1_old(1:3,JELEM) = DZERO
			deltaN2_old(1:3,JELEM) = DZERO
		ELSE
			EpsilonV1(JELEM) = DZERO
			EpsilonV2(JELEM) = DZERO
			EpsilonV1_old(JELEM) = DZERO
			EpsilonV2_old(JELEM) = DZERO
			
			deltaN1(1:3,JELEM) = DZERO
			deltaN2(1:3,JELEM) = DZERO
			deltaN1_old(1:3,JELEM) = DZERO
			deltaN2_old(1:3,JELEM) = DZERO
		ENDIF
	  ENDIF
	  
C     RETRIEVE LDPM SIMULATION TIME AT DATA EXCHANGE
      scaledtime_LDPM = timeScalingFac*time_LDPM
      scaledtime_old_LDPM = timeScalingFac*time_LDPM_old
      DTIME_coupling = scaledtime_LDPM - scaledtime_old_LDPM


C     CALCULATE CRACKED MATERIAL VOLUME (VC)
	  IF (fractureFlag == ONE) THEN ! If fracture analysis flag is on
		VC(1,JELEM) = Af(1,JELEM)*deltaN1(1,JELEM) + Af(2,JELEM)*deltaN1(2,JELEM) 
	1               + Af(3,JELEM)*deltaN1(3,JELEM)
		VC(2,JELEM) = Af(4,JELEM)*deltaN2(1,JELEM) + Af(5,JELEM)*deltaN2(2,JELEM) 
	1               + Af(6,JELEM)*deltaN2(3,JELEM)

		VC_old(1,JELEM) = Af(1,JELEM)*deltaN1_old(1,JELEM) + Af(2,JELEM)*deltaN1_old(2,JELEM) 
	1               	+ Af(3,JELEM)*deltaN1_old(3,JELEM)
		VC_old(2,JELEM) = Af(4,JELEM)*deltaN2_old(1,JELEM) + Af(5,JELEM)*deltaN2_old(2,JELEM) 
	1               	+ Af(6,JELEM)*deltaN2_old(3,JELEM)
	  ELSE
		VC(1:2,JELEM) = DZERO
		VC_old(1:2,JELEM) = DZERO
		deltaN1(1:3,JELEM) = DZERO
		deltaN2(1:3,JELEM) = DZERO
		deltaN1_old(1:3,JELEM) = DZERO
		deltaN2_old(1:3,JELEM) = DZERO
	  ENDIF
	  
C 	  UPDATE FIELD VARIABLES
      P1    = U(1)
      P1OLD = U(1) - DU(1,NRHS)
      P2    = U(2)
      P2OLD = U(2) - DU(2,NRHS)
		
      EpsilonV1NEW = EpsilonV1(JELEM)
      EpsilonV1OLD = EpsilonV1_old(JELEM)
      EpsilonV2NEW = EpsilonV2(JELEM)
      EpsilonV2OLD = EpsilonV2_old(JELEM)
		
      Vc1 = VC(1,JELEM)
      Vc1OLD = VC_old(1,JELEM)
      Vc2 = VC(2,JELEM)
      Vc2OLD = VC_old(2,JELEM)
	  
C     INCREMENTS AND RATES AT TET POINTS
      Delta_P1 = DU(1,NRHS)
      DP1DT = Delta_P1/DTIME
      Delta_P2 = DU(2,NRHS)
      DP2DT = Delta_P2/DTIME

C     INCREMENTS AND RATES OF COUPLING VARIABLES
      Delta_EpsilonV1 = EpsilonV1NEW - EpsilonV1OLD
      Delta_EpsilonV2 = EpsilonV2NEW - EpsilonV2OLD
      Delta_Vc1 = Vc1 - Vc1OLD
      Delta_Vc2 = Vc2 - Vc2OLD
		
      IF (DTIME_coupling == DZERO) THEN
		dEpsilonV1_dtime = DZERO
		dEpsilonV2_dtime = DZERO
		dVc1_dtime = DZERO
		dVc2_dtime = DZERO
      ELSE
		dEpsilonV1_dtime = Delta_EpsilonV1/DTIME_coupling
		dEpsilonV2_dtime = Delta_EpsilonV2/DTIME_coupling
		dVc1_dtime = Delta_Vc1/DTIME_coupling
		dVc2_dtime = Delta_Vc2/DTIME_coupling
      END IF
	  
C     LOAD NODAL PORE PRESSURES TO NAMED PIPE VARIABLES
      FLM2LDPM_DATA(tetIndex(1,JELEM),1) = P1OLD
      FLM2LDPM_DATA(tetIndex(2,JELEM),1) = P2OLD

C*********************************************************************		
C     FLM ALGORITHM
C
      IF(JTYPE == 1) THEN ! 1 == Original FLM Algorithm
	  ! This section is intended to be hidden
      ENDIF

C     INTERPOLATION OF TET POINT VALUES
      CBAR = G2(JELEM)*C1 + G1(JELEM)*C2
      SBAR = G2(JELEM)*S1 + G1(JELEM)*S2
		
C     PASS STATE VARIABLES TO SVARS FOR VISUALIZATION
      SVARS(1) = EpsilonV1NEW
      SVARS(2) = EpsilonV2NEW
      SVARS(3) = Vc1
      SVARS(4) = Vc2

C     VISUALIZATION (values to be passed to dummy mesh using the GlobalVarStorage module)
      UserVar(jelem,1,1) = SVARS(1) ! Volumetric strain of tet 1
      UserVar(jelem,2,1) = SVARS(2) ! Volumetric strain of tet 2
      UserVar(jelem,3,1) = SVARS(3) ! Cracked material volume of tet 1
      UserVar(jelem,4,1) = SVARS(4) ! Cracked material volume of tet 2
      UserVar(jelem,1,2) = SVARS(1) ! Volumetric strain of tet 1
      UserVar(jelem,2,2) = SVARS(2) ! Volumetric strain of tet 2
      UserVar(jelem,3,2) = SVARS(3) ! Cracked material volume of tet 1
      UserVar(jelem,4,2) = SVARS(4) ! Cracked material volume of tet 2

C*********************************************************************		
C     ASSEMBLE AMATRX AND RHS
C

	  ! This section is intended to be hidden
	  RHS(1,1) = DZERO 
	  RHS(2,1) = DZERO
	  	  
	  AMATRX(1,1) = DZERO ! rate term
	1			  + DZERO ! gradient term
	2			  + DZERO ! source term
	  AMATRX(1,2) = DZERO ! rate term
	1			  + DZERO ! gradient term
	2			  + DZERO ! source term
	  
	  AMATRX(2,1) = DZERO ! rate term
	1			  + DZERO ! gradient term
	2			  + DZERO ! source term
	  AMATRX(2,2) = DZERO ! rate term
	1			  + DZERO ! gradient term
	2			  + DZERO ! source term
	  
      RETURN
      END SUBROUTINE UEL

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)                                                                   
	  ! READ EXTERNAL FILES FOR THE GEOMETRY OF EDGE ELEMENTS AND COUPLING DATA 

	  USE GlobalVarStorage
	  USE ModuleFLE
	  
      INCLUDE 'aba_param.inc'

      DIMENSION TIME(2)

	  INTEGER			   :: ios
	  INTEGER			   :: MYTHREADID
      CHARACTER(LEN = 256) :: Chemin
	  CHARACTER(LEN = 256) :: inpfile
	  CHARACTER(LEN = 256) :: MPfile
	  CHARACTER(LEN = 256) :: line
	  CHARACTER(LEN = 256) :: separator
	  LOGICAL              :: file_exists
	  
	  MYTHREADID = GET_THREAD_ID()

      IF (LOP == 0) THEN ! Start of the analysis
		IF (MYTHREADID == 0) THEN ! If the first thread
			CALL GETJOBNAME(JOBNAME,LENJOBNAME)
			CALL GETOUTDIR(OUTDIR,LENOUTDIR)
		
			separator = outdir(1:1)
		
			IF (separator == '/') THEN
				path_separator = '/'
			ELSE
				path_separator = '\'
			ENDIF
	  
C***************** READ GEOMETRY DATA FILES *************************
			Chemin = trim(trim(OUTDIR)//path_separator//trim(JOBNAME)//'-geom.dat')
			
			OPEN(UNIT=15,FILE=Chemin,FORM="formatted",
	1         	ACCESS="sequential",IOSTAT=ios)
      
			IF (ios == 0) THEN
				! WRITE(*,*) 'Reading external edge geometry files'
				READ(UNIT=15,FMT=*)
				READ(UNIT=15,FMT=*)
				READ(UNIT=15,FMT='(I8)') nnode_FLM
				! WRITE(*,*) 'nnode_FLM', nnode_FLM
				READ(UNIT=15,FMT=*)
				READ(UNIT=15,FMT='(I8)') nelem
				! WRITE(*,*) 'nelem_FLM', nelem
				IF (allocated(TET_DATA) == 0) THEN
					ALLOCATE(TET_DATA(nnode_FLM,3))
					ALLOCATE(ELEM_DATA(nelem,28))
				ENDIF
				READ(UNIT=15,FMT=*)
				! WRITE(*,*) 'Reading tet point coordinates'
				DO i=1,nnode_FLM
					READ(UNIT=15,FMT=*) (TET_DATA(i,j), j=1,3)
				END DO
			
				READ(UNIT=15,FMT=*)
				! WRITE(*,*) 'Reading edge element data'
				DO i=1,nelem
					READ(UNIT=15,FMT=*) (ELEM_DATA(i,j), j=1,28)
				END DO
				WRITE(*,*) 'Reading external edge geometry files complete.'
			END IF
			CLOSE(15)
		
C***************** READ ABAQUS INPUT FILES *************************

			inpfile = trim(trim(outdir)//path_separator//trim(jobname)//'.inp')

			OPEN(UNIT=16,FILE=inpfile,FORM="formatted",
	1			ACCESS="sequential",IOSTAT=ios)
			
			IF (ios == 0) THEN
				WRITE(*,*) 'Read Abaqus input file'
				
				! Skip down to *Depvar in INP file
				DO WHILE (index(line,'*Depvar') == 0)              
					READ(UNIT=16,FMT='(A)') line              
				END DO
				! Read in n_stored_vars
				READ(UNIT=16,FMT='(I8)') n_stored_vars		
			END IF
			CLOSE(16) !closes file
			
			nelemUel = nelem

			IF (allocated(UserVar) == 0) THEN
				ALLOCATE(UserVar(nelem,n_stored_vars,n_int_pnts))
			END IF

C***************** READ COUPLING EXCHANGE DATA FILES ***************
			MPfile = trim(trim(OUTDIR)//path_separator//trim(JOBNAME)//'-exchangeMP.dat')
			
			OPEN(UNIT=17,FILE=MPfile,FORM="formatted",
	1         	ACCESS="sequential",IOSTAT=ios)
      
			IF (ios == 0) THEN
				WRITE(*,*) 'Read exchange time file'
				READ(UNIT=17,FMT=*)
				READ(UNIT=17,FMT=*)
				READ(UNIT=17,FMT='(I8)') nExchangeTime
				WRITE(*,*) 'nExchangeTime', nExchangeTime
				READ(UNIT=17,FMT=*)
				IF (allocated(ExchangeTimeList) == 0) THEN               
					ALLOCATE(ExchangeTimeList(nExchangeTime))
				ENDIF
				READ(UNIT=17,FMT=*) ExchangeTimeList
				
				WRITE(*,*) 'Reading exchange time file completed.'
			END IF
			CLOSE(17)
			
			kinc_in_Lop1 = 0
			
			ExchangeFlag = 0
			exchangeindex = 1
			CurrentExchangeTime = ExchangeTimeList(exchangeindex)
		END IF
		
		U2Vpipefile = trim(trim(outdir)//path_separator//'FLM2LDPM.pipe')
		V2Upipefile = trim(trim(outdir)//path_separator//'LDPM2FLM.pipe')
		
		WRITE(*,*) 'Loading named pipe file completed.'
		
	  ELSE IF (LOP == 1) THEN ! Start of the increment
C**************** SETTINGS OF REGULAR DATA EXCHANGE ****************
		
		! UEL will run twice following LOP == 1 condition for the first increment 
		! (1 for mass calc, 1 for regular AMATRX calc), here we skip the first run
		! by using kinc_in_Lop1, to make sure that data is only exchanged once 
		! for the first increment
		
		IF (kinc_in_Lop1 == 0) THEN
			kinc_in_Lop1 = 1
		ELSE
			IF (ExchangeFlag == 1) THEN
				OPEN(UNIT=V2U,FILE='LDPM2FLM.pipe',DEFAULTFILE=trim(OUTDIR),
	1	     	FORM='formatted',STATUS='old',ACTION='read',ACCESS='stream')
			
				DO i=1,nelem_LDPM
					READ(V2U,'(ES24.17)') (LDPM2FLM_DATA(i,j), j=1,13)
				END DO
				
				READ(V2U,'(ES24.17)') time_LDPM
				
				CLOSE(V2U)
			  
				OPEN(UNIT=U2V,FILE='FLM2LDPM.pipe',DEFAULTFILE=trim(OUTDIR),
	1			FORM='formatted',STATUS='old',ACTION='write',ACCESS='stream')
			  
				DO i=1,nnode_FLM
					WRITE(U2V,'(ES24.17)') (FLM2LDPM_DATA(i,j), j=1,1)
				END DO
				
				WRITE(U2V,'(ES24.17)') period_FLM
				WRITE(U2V,'(ES24.17)') TIME(2)

				FLUSH(U2V)
				CLOSE(U2V)
				
				write(*,*) "data has exchanged at time", TIME(2)

			ELSE
				CONTINUE
			END IF
		END IF
		
	  ELSE IF (LOP == 2) THEN ! End of the increment	
C***************** SETTINGS OF DATA EXCHANGE *************************

		IF (ExchangeFlag == 1) THEN
			ExchangeFlag = 0
		ELSE
			IF (time(2) >= CurrentExchangeTime) THEN
				time_LDPM_old = time_LDPM
				LDPM2FLM_DATA_old(1:nelem_LDPM,1:13) = LDPM2FLM_DATA(1:nelem_LDPM,1:13)
				ExchangeFlag = 1
				exchangeindex = exchangeindex + 1
				IF (exchangeindex > nExchangeTime) THEN
					CurrentExchangeTime = -LOG(0.0d0) ! After the last exchange, set next exchange time to infinity
				ELSE
					CurrentExchangeTime = ExchangeTimeList(exchangeindex)
				END IF
			END IF

		END IF
		
	  ELSE IF (LOP == 5) THEN ! Start of the step
C***************** INITIAL DATA EXCHANGE *****************
		
		OPEN(UNIT=V2U,FILE='LDPM2FLM.pipe',DEFAULTFILE=trim(OUTDIR),
	1	     FORM='formatted',STATUS='old',ACTION='read',ACCESS='stream',IOSTAT=ios)
	
		WRITE(*,*) 'Retriving LDPM analysis settings'
		READ(V2U,'(I8)') nelem_LDPM
		WRITE(*,*) 'nelem_LDPM', nelem_LDPM
		READ(V2U,'(ES24.17)') period_LDPM
		WRITE(*,*) 'period_LDPM', period_LDPM
		
		CLOSE(V2U)

		OPEN(UNIT=U2V,FILE='FLM2LDPM.pipe',DEFAULTFILE=trim(OUTDIR),
	1	     FORM='formatted',STATUS='old',ACTION='write',ACCESS='stream')
		
		WRITE(*,*) 'Sending FLM analysis settings'
		WRITE(U2V,'(I8)') nnode_FLM
		MPtypeFLM = 2
		WRITE(U2V,'(I8)') MPtypeFLM

		FLUSH(U2V)
		CLOSE(U2V)
		
		IF (allocated(FLM2LDPM_DATA) == 0) THEN
			allocate(FLM2LDPM_DATA(nnode_FLM,1),LDPM2FLM_DATA(nelem_LDPM,13),LDPM2FLM_DATA_old(nelem_LDPM,13))
			LDPM2FLM_DATA = 0.0d0
			LDPM2FLM_DATA_old = 0.0d0
			WRITE(*,*) 'shape(LDPM2FLM_DATA)', shape(LDPM2FLM_DATA)
		ENDIF
		
	  ELSE IF (LOP == 6) THEN ! End of the step
C***************** CLEANNING UP FILES *************************

		CALL SYSTEM('find '//trim(outdir)//' -maxdepth 1 -type p -delete')
      END IF

      RETURN
      END SUBROUTINE UEXTERNALDB
	    
      SUBROUTINE UMATHT(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,dtemp,
	1 					dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,props,nprops,
	2					coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)
	  
      USE GlobalVarStorage
      INCLUDE 'aba_param.inc'

      CHARACTER*80 cmname
      DIMENSION dudg(ntgrd), flux(ntgrd), dfdt(ntgrd),
	1 			dfdg(ntgrd,ntgrd), statev(nstatv), dtemdx(ntgrd),
	2			time(2), predef(1), dpred(1), props(nprops), coords(3) 

C	 Information exchange between global module and local state variables (statev) 
C	 at each Tet point of each element in the dummy mesh

      noffset = noel - nelemUel
      statev(1:nstatv) = UserVar(noffset,1:nstatv,npt)       

      RETURN
      END SUBROUTINE UMATHT