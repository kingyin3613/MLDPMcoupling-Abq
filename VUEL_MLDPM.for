C**********************************************************************C
C**********************************************************************C
C                                                                      C
C           IMPLEMENTATION OF M-LDPM FOR ABAQUS                        C
C                                                                      C
C           BY Hao Yin                                                 C
C              Erol Lale                                               C
C              Matthew Troemner                                        C
C              Mohammed Alnaggar                                       C
C              Lifu Yang                                               C
C              Lei Shen                                                C
C                                                                      C
C              Jan 13 2024                                             C
C                                                                      C
C**********************************************************************C
C**********************************************************************C

      module ModuleLDPM
      implicit none 
      
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)

      ! Declare integers
      integer                                 :: nomaxel
      integer                                 :: print_binary
      integer                                 :: svars_nodal_stress
      integer,  dimension(:,:),   allocatable :: elecon                
      integer,  dimension(:),     allocatable :: svars_selection 

      ! Declare double precision floats
      real(dp)                                :: b
      real(dp), dimension(3)                  :: svars_freq
      real(dp), dimension(:,:,:), allocatable :: FacetVar           
      real(dp), dimension(:),     allocatable :: next_print
	  
      end module


      module ModuleCoupling
      ! Global variables storage for coupling and visualization
	  
	  use ModuleLDPM
      implicit none 

      ! Declare integers
	  integer,  dimension(:),     allocatable :: elemflag
	  integer                                 :: nelem_FLM
	  integer                                 :: nfieldFLM
	  integer                                 :: nincFLM
	  integer                                 :: nnode_FLM
	  integer                                 :: MPtypeFLM
	  integer                                 :: exchangeindex
	  integer                                 :: nExchangeTime
	  integer                                 :: ExchangeFlag
	  integer,  dimension(:,:),   allocatable :: tetIndex
      integer                                 :: iparticle, iagg ! aggregate size
	  
	  integer, parameter :: V2U = 905 ! file unit number for named pipes
	  integer, parameter :: U2V = 906 ! file unit number for named pipes
	  
      ! Declare double precision floats
	  real(dp), dimension(:,:),   allocatable :: tet_data
      real(dp), dimension(:,:),   allocatable :: elem_data
	  real(dp), dimension(:,:),   allocatable :: step_data
	  real(dp), dimension(:,:),   allocatable :: result_data
	  
      ! Declare double precision floats for multiphysics
      real(dp), dimension(:,:,:), allocatable :: FacetVarMP
	  real(dp), dimension(:,:),   allocatable :: FLM2LDPM_DATA
	  real(dp), dimension(:,:),   allocatable :: FLM2LDPM_DATA_old
	  real(dp), dimension(:,:),   allocatable :: LDPM2FLM_DATA
	  real(dp), dimension(:),     allocatable :: ExchangeTimeList
	  real(dp) 							      :: period_FLM
	  real(dp) 							      :: time_FLM
	  real(dp) 							      :: time_FLM_old
	  real(dp) 							      :: CurrentExchangeTime
      ! AgingHistTemp, store history value of alpha_c, alpha_s, lambda for all tets at current step
      real(dp), dimension(:,:),   allocatable :: AgingHistTemp
      ! aggregate size info, for all particles. [NodeNum,x-coord,y-coord,z-coord,Size]
      real(dp), dimension(:,:),   allocatable :: AggSize
      ! ASR history matrix  for all tets at current step, [Z1, Z2, Mw1, Mw2]
      real(dp), dimension(:,:),   allocatable :: ASRHist
	  
	  character(len = 256) :: U2Vpipefile
	  character(len = 256) :: V2Upipefile
	  save
      end module ModuleCoupling


      module ModuleMPFlags
      ! Global variables storage for MPFlags
      implicit none
	  
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)
	  
      real(dp), dimension(11)                 :: MPflags
	  
      save
      end module ModuleMPFlags

	  
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     USER ELEMENT SUBROUTINE                                          C
C                                                                      C
C**********************************************************************C
C**********************************************************************C

      SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars,
	1                 energy,
	2                 nnode,ndofel,props,nprops,jprops,njprops,
	3                 coords,mcrd,u,du,v,a,
	4                 jtype,jElem,
	5                 time,period,dtimeCur,dtimePrev,kstep,kinc,
	6                 lflags,
	7                 dMassScaleFactor,
	8                 predef,npredef,
	9                 jdltyp, adlmag)

	  use ModuleCoupling

      implicit none  

      integer(KIND=4),  parameter :: ndime                = 3
      integer(KIND=4),  parameter :: ndofn                = 6    
      integer(KIND=4),  parameter :: nshr                 = 3  
      integer(KIND=4),  parameter :: nsvint               = 75
      integer(KIND=4),  parameter :: ntens                = 6
      integer(KIND=4),  parameter :: ndi                  = 3
      integer(KIND=4),  parameter :: nface                = 12

      integer(KIND=4),  parameter :: jMassCalc            = 1
      integer(KIND=4),  parameter :: jIntForceAndDtStable = 2  
      integer(KIND=4),  parameter :: jExternForce         = 3

      integer(KIND=4),  parameter :: iProcedure           = 1    
      integer(KIND=4),  parameter :: iNlgeom              = 2  
      integer(KIND=4),  parameter :: iOpCode              = 3   
      integer(KIND=4),  parameter :: nFlags               = 3 

      integer(KIND=4),  parameter :: jDynExplicit         = 17

      ! energy array index parameters
      integer(KIND=4),  parameter :: iElPd                = 1
      integer(KIND=4),  parameter :: iElCd                = 2
      integer(KIND=4),  parameter :: iElIe                = 3 
      integer(KIND=4),  parameter :: iElTs                = 4
      integer(KIND=4),  parameter :: iElDd                = 5
      integer(KIND=4),  parameter :: iElBv                = 6
      integer(KIND=4),  parameter :: iElDe                = 7
      integer(KIND=4),  parameter :: iElHe                = 8
      integer(KIND=4),  parameter :: iElKe                = 9
      integer(KIND=4),  parameter :: iElTh                = 10
      integer(KIND=4),  parameter :: iElDmd               = 11
      integer(KIND=4),  parameter :: nElEnergy            = 12

      ! predefined variables indices
      integer(KIND=4),  parameter :: iPredValueNew        = 1
      integer(KIND=4),  parameter :: iPredValueOld        = 2
      integer(KIND=4),  parameter :: nPred                = 2

      ! time indices
      integer(KIND=4),  parameter :: iStepTime            = 1
      integer(KIND=4),  parameter :: iTotalTime           = 2
      integer(KIND=4),  parameter :: nTime                = 2

      ! Declare integers
      integer(KIND=4)                         :: fMat
      integer(KIND=4)                         :: i
      integer(KIND=4)                         :: icol
      integer(KIND=4)                         :: icpu
      integer(KIND=4)                         :: idfile
      integer(KIND=4)                         :: idime
      integer(KIND=4)                         :: idofel
      integer(KIND=4)                         :: idofn
      integer(KIND=4)                         :: ielem
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: inode
      integer(KIND=4)                         :: ios
      integer(KIND=4)                         :: isvinc
	  integer(KIND=4)                         :: j
      integer(KIND=4)                         :: jcol
      integer(KIND=4)                         :: jdltyp
      integer(KIND=4)                         :: jtype
      integer(KIND=4),  dimension(njprops)    :: jprops
      integer(KIND=4)                         :: jnode
      integer(KIND=4),  dimension(nblock)     :: jElem        
      integer(KIND=4)                         :: kblock
      integer(KIND=4)                         :: Kfindloc
      integer(KIND=4)                         :: KINC
      integer(KIND=4)                         :: KPROCESSNUM  
      integer(KIND=4)                         :: KSTEP 
      integer(KIND=4)                         :: LOP
      integer(KIND=4)                         :: LRESTART
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4),  dimension(nFlags)     :: lflags 
      integer(KIND=4)                         :: mcrd
      integer(KIND=4)                         :: nblock
      integer(KIND=4)                         :: ndofel
      integer(KIND=4)                         :: njprops
      integer(KIND=4)                         :: nnode
      integer(KIND=4)                         :: ninpt
      integer(KIND=4)                         :: npredef
      integer(KIND=4)                         :: nprops
      integer(KIND=4)                         :: NUMPROCESSES
      integer(KIND=4)                         :: nsvars
      integer(KIND=4)                         :: svarsidfile 
      integer(KIND=4)                         :: ORIGINAL_FPE_FLAGS
      integer(KIND=4)                         :: NEW_FPE_FLAGS
	  
      ! Declare double precision floats
      real(dp), dimension(nblock, ndofel)     :: a
      real(dp), dimension(1,3)                :: accelVec
      real(dp), dimension(nblock)             :: adlmag     
      real(dp)                                :: al
      real(dp)                                :: aLe
      real(dp)                                :: aLen
      real(dp), dimension(nblock,ndofel,ndofel) :: amass  
      real(dp), dimension(ndofel)             :: BodyL      
      real(dp), dimension(mcrd)               :: coordgp   
      real(dp), dimension(nblock,nnode,mcrd)  :: coords           
      real(dp), dimension(mcrd)               :: cross
      real(dp), dimension(mcrd)               :: cross0
      real(dp)                                :: cur_step
      real(dp), dimension(ndofel)             :: DeltaU
      real(dp)                                :: density  
      real(dp), dimension(ndofel)             :: deigen      
      real(dp)                                :: diff
      real(dp), dimension(nblock)             :: dMassScaleFactor
      real(dp), dimension(ntens)              :: dstran
      real(dp)                                :: DTIME
      real(dp)                                :: dtimeCur
      real(dp)                                :: dtimePrev
      real(dp), dimension(nblock)             :: dtimeStable
      real(dp), dimension(nblock,ndofel)      :: du      
      real(dp)                                :: dx
      real(dp)                                :: dy
      real(dp)                                :: dz
      real(dp), dimension(mcrd,nnode)         :: elcoord0
      real(dp), dimension(mcrd,nnode)         :: elcoord1
      real(dp), dimension(mcrd,nnode)         :: elcoord
      real(dp), dimension(ndofel,ndofel)      :: elmass
      real(dp), dimension(ndofel,ndofel)      :: elmassInvStiff
      real(dp), dimension(ndofel,ndofel)      :: elstif
      real(dp)                                :: endTime
      real(dp), dimension(nblock,nElEnergy)   :: energy
      real(dp)                                :: epsV
      real(dp)                                :: factorStable  
      real(dp), dimension(nsvint)             :: facestatev
      real(dp)                                :: fcstrainEN
      real(dp), dimension(ndofn)              :: force  
      real(dp), dimension(ndofel)             :: forceVec
      real(dp)                                :: omax
      real(dp)                                :: period
      real(dp), dimension(nblock,nnode,npredef,nPred)  :: predef    
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(nblock,ndofel)      :: rhs   
      real(dp), dimension(nsvint)             :: statevLocal
      real(dp), dimension(ntens)              :: stress
      real(dp), dimension(ntens)              :: stran      
      real(dp), dimension(nblock,nsvars)      :: svars  
      real(dp)                                :: tetVol0
      real(dp)                                :: tetVol      
      real(dp), dimension(nTime)              :: time
      real(dp)                                :: startTime
      real(dp)                                :: time_step
      real(dp), dimension(nblock,ndofel)      :: u
      real(dp), dimension(nblock,ndofel)      :: v
      real(dp), dimension(ndofel,ndofel)      :: veigen
      real(dp), dimension(mcrd)               :: va
      real(dp), dimension(mcrd)               :: vb
      real(dp), dimension(mcrd)               :: vc
      real(dp), dimension(mcrd)               :: va0
      real(dp), dimension(mcrd)               :: vb0
      real(dp), dimension(mcrd)               :: vc0
      real(dp), dimension(3)                  :: vnr
      real(dp), dimension(3)                  :: vmr
      real(dp), dimension(3)                  :: vlr

      ! Declare booleans
      logical(KIND=4)                         :: direxist
      logical(KIND=4)                         :: pipexist
      logical(KIND=4)                         :: LargeDefFlag
      logical(KIND=4)                         :: RateEffectFlag

      ! Declare characters
      character(80)                           :: Chemin
      character(256)                          :: filename    
      character(80)                           :: JOBNAME
      character(256)                          :: OUTDIR
      character(1)                            :: path_separator       
      character(1)                            :: separator
      character(256)                          :: svarsfile
      character(256)                          :: svarsfilepath
	  character(256)                          :: trimmed_svarsfilepath

      ! Assign values
      LargeDefFlag                     = .true.  
      RateEffectFlag                   = .false.
      factorStable                     = 0.900d0
	  
	  ! Extract job data
	  CALL VGETJOBNAME( JOBNAME, LENJOBNAME ) ! Jobname
	  CALL VGETOUTDIR( OUTDIR, LENOUTDIR )    ! Work directory

	  separator = OUTDIR(1:1)
	  
	  if (separator == '/') then
		 path_separator = '/'
	  else
		 path_separator = '\'
	  endif

C======================================================================C
C     MPI Files                                                        C
C======================================================================C

	  ! Extract simulation parameters
	  CALL VGETNUMCPUS( NUMPROCESSES ) 
	  CALL VGETRANK( KPROCESSNUM ) 
	  icpu = KPROCESSNUM
		  
C======================================================================C
C     Model Checks                                                     C
C======================================================================C

      density = props(1) ! Used for time step calc

      ! Preliminaries    
      if(jtype /= 4) then
         write(*,*)'Incorrect element type'
         call XPLB_EXIT ! call exit
      endif 

      if(nsvars < nface*nsvint) then
         write(*,*)'Increase the number of SDVs to',nface*nsvint
         call XPLB_EXIT ! call exit
      endif 

C======================================================================C
C     Named Pipe Files                                                 C
C======================================================================C

		V2Upipefile = trim(trim(outdir)//path_separator//'LDPM2FLM.pipe')
		U2Vpipefile = trim(trim(outdir)//path_separator//'FLM2LDPM.pipe')

        call VUEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

		! read in all Multiphysics flags and stored in Module
		call MP2initializer_new(props,nprops,jprops,njprops)          
	  
	  
C======================================================================C
C     Select Procedures                                                C
C======================================================================C

      !Procedure type: Explicit dynamic 
      if (lflags(iProcedure) == jDynExplicit) then 

        ! Mass calculation
        if (lflags(iOpCode) == jMassCalc) then               

C======================================================================C   
C     CALCULATE MASS MATRIX                                            C
C======================================================================C      

! This section is intended to be hidden

        else if (lflags(iOpCode) == jIntForceAndDtStable) then

C======================================================================C
C     Perform Multiphysics Coupling 
C======================================================================C
			   
               call MP2coupler(ielem,FacetVar(:,:,ielem),FacetVarMP(:,:,ielem),
	1							nprops,props,njprops,jprops,time,dtimeCur,kinc,kstep,
	2							nblock,nsvars,tetVol0,period,mcrd,nnode,elcoord0,
	3							elcoord,elecon(ielem,2:5),AggSize,LargeDefFlag,iparticle,
	4							ASRHist(ielem,:),AgingHistTemp(ielem,:),
	5							nsvint,svars(kblock,:))

C======================================================================C
C     Internal force calculations
C======================================================================C

! This section is intended to be hidden


			  ! Update the force vector        
			  do idofn = 1,ndofel                           
                rhs(kblock,idofn) = rhs(kblock,idofn) + forceVec(idofn)
			  end do

			  ! Internal energy calculation 
			  energy(kblock,iElIe) = energy(kblock,iElIe) + fcstrainEN

C======================================================================C
C     Stable time calculations
C======================================================================C
        
! This section is intended to be hidden

            end do   ! END OF KBLOCK LOOP

        elseif(lflags(iOpCode) == jExternForce) then   

C======================================================================C
C     Calculate external forces
C======================================================================C


! This section is intended to be hidden

        end if  ! end of iOpCode if statement
      end if  ! end of iProcedure if statement
  
      
      return
      end subroutine VUEL


C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE Initializer for multiphysics coupling                 C
C                                                                      C
C**********************************************************************C
C**********************************************************************C

      subroutine initMP

      use ModuleCoupling
      implicit none  
        
      ! Declare integers
      integer(KIND=4)                         :: i,j
      integer(KIND=4)                         :: ios
      integer(KIND=4)                         :: io
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: itet
      integer(KIND=4)                         :: iwhere
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4)                         :: nelcnt
      integer, parameter                      :: zero = 0

      ! Declare characters
      character(256)                          :: filename
      character(256)                          :: facetMPname
      character(256)                          :: input
      character(256)                          :: jobname
      character(256)                          :: line
      character(256)                          :: outdir
      character(1)                            :: path_separator
      character(256)                          :: resultname
      character(256)                          :: separator

C======================================================================C
C     Read edge information from edge files                            C
C======================================================================C
   
      ! Open edgedata file for reading.
      call vgetjobname(jobname,lenjobname)
      call vgetoutdir(outdir,lenoutdir)

      separator = outdir(1:1)
      if (separator == '/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif
    
C*****************************************************************
C       create interpolation table
C       python code handles this part
C       tetID --> edge node ID
C*****************************************************************
      facetMPname = outdir(1:lenoutdir)//path_separator//
     1                       jobname(1:lenjobname)//'_facetMP.dat' 
      open(unit=224,file=facetMPname(1:lenoutdir+
     1           lenjobname+13),status='old')
      
        if (allocated(FacetVarMP) == 0.0d0) then               
            allocate(FacetVarMP(1,12,nomaxel))
      endif
      nelcnt=0
10    read (224, '(A)', end = 20) line
        ! write(*,*) 'line', line
      iwhere = index (line, 'Tet index:')
      if (iwhere.ne.0) then
          read(line(11:20),'(I8)') itet !int(line(11:20))    
          nelcnt=nelcnt+1      
  
          do iface=1,12
              read(224,*,ERR=200) FacetVarMP(1:1,iface,nelcnt)          
          end do  
              
      endif
      goto 10
20    continue

        write(*,*) 'Successfully read FacetVarMP data'
200   close(224) !closes file 

      return
      end subroutine initMP
	  

C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE Initializer for multiphysics coupling                 C
C                                                                      C
C**********************************************************************C
C**********************************************************************C

      subroutine get_exchange_info

      use ModuleCoupling
      implicit none  
        
      ! Declare integers
      integer(KIND=4)                         :: i
      integer(KIND=4)                         :: ios
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer, parameter                      :: zero = 0

      ! Declare characters
      character(256)                          :: exchangeMPname
      character(256)                          :: jobname
      character(256)                          :: outdir
      character(1)                            :: path_separator
      character(256)                          :: separator

C======================================================================C
C     Read edge information from edge files                            C
C======================================================================C
   
      ! Open edgedata file for reading.
      call vgetjobname(jobname,lenjobname)
      call vgetoutdir(outdir,lenoutdir)

      separator = outdir(1:1)
      if (separator == '/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif
    
C*****************************************************************
C       create exchange time list
C*****************************************************************

      exchangeMPname = outdir(1:lenoutdir)//path_separator//
     1                       jobname(1:lenjobname)//'-edgeEle-exchangeMP.dat'
	 
      open(unit=227,file=exchangeMPname(1:lenoutdir+
     1           lenjobname+24),status='old')
      
	  if (ios == 0) THEN
		write(*,*) 'Read exchange time file'
		read(unit=227,fmt=*)
		read(unit=227,fmt=*)
		read(unit=227,fmt='(i8)') nExchangeTime
		write(*,*) 'nExchangeTime', nExchangeTime
		read(unit=227,fmt=*)
		if (allocated(ExchangeTimeList) == 0) then               
			allocate(ExchangeTimeList(nExchangeTime))
		endif
		read(unit=227,fmt=*) ExchangeTimeList
	  endif
	  ExchangeFlag = 0
	  exchangeindex = 1
	  CurrentExchangeTime = ExchangeTimeList(exchangeindex)
      write(*,*) 'Successfully read exchange time data'
	  close(227) !closes file 

      return
      end subroutine get_exchange_info


C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE VUEXTERNALDB                                          C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      subroutine VUEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

        use ModuleCoupling
      use ModuleLDPM
      implicit none  


      ! Declare integers
      integer(KIND=4)                         :: KINC
      integer(KIND=4)                         :: KSTEP
      integer(KIND=4)                         :: LRESTART
      integer(KIND=4)                         :: LOP


      ! Declare double precision floats
      real(dp), dimension(2)                  :: time
      real(dp)                                :: DTIME


      ! Declare characters    
      character(256)                          :: OUTDIR
      character(256)                          :: Chemin
      character(256)                          :: JOBNAME        

      call get_svars_config
	  call get_exchange_info

      call initMP

      return
      end subroutine VUEXTERNALDB


C**********************************************************************C
C**********************************************************************C
C                                                                      C
C            MAIN MULTIPHYSICS TWO-WAY COUPLER SUBROUTINES             C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
	
      subroutine MP2coupler(tetID,facetgeom,facetgeomMP,nprops,props,njprops,jprops,
	1						simulationTime,dtimeCur,kinc,kstep,nblock,nsvars,Vol0,
	2						period,mcrd,nnode,elcoord0,
	3						elcoord,nconnec,AggSize,LargeDefFlag,iparticle,ASRHistTet,
	4						AgingHist,nsvint,svars_elem)

      ! This subroutine will be called before the calculation of 12 facets
      ! All calculations in this subroutine are only for one tet which contains 12 facets
      ! svars is the most important variable for information transfer between VUEL_LDPM & Abaqus
      ! AgingHist is the temporary matrix storing alpha_c, alpha_s and lambda for the tet element
        
      ! tetID: tet ID (ielem in subroutine VUEL)
      ! facetgeom: facet constant geometry properties
      ! facetgeomMP: facet constant geometry properties for multiphysics
      ! simulationTime: time(:) in VUEL, the very short time, before scaling
      ! kinc = current time increment ID
      ! kstep = current step ID
      ! nblock: the amount of tets in the current block
      ! nsvars: 12 * svars number
      ! dtimeCur = current time increment (dt)
      ! elcoord0 = nodal coordinates in the reference (material) configuration
      ! svars_elem = svars for a certain tet to be updated in the subroutine

      ! MPfields: imported multiphysics fields. line->each facet; column->each field
      ! stvMP: Multi-Physics state variables (internal variables)
      
C     stvMP(1) Total Normal Imposed Strain
C     stvMP(2) Imposed Normal Strain Increment
C     stvMP(3) Imposed M Shear Strain Increment
C     stvMP(4) Imposed L Shear Strain Increment
C     stvMP(5) Imposed Volumetric Strain
C     stvMP(6) Imposed Volumetric Strain Increment ! [H.Y. changed Jan 13 2022]

C     stvMP(7) Hydration degree
C     stvMP(8) Silica fume degree
C     stvMP(9) Total reaction degree
C     stvMP(10) Aging degree

C     High Temperature Degradation
C     stvMP(11) = stv[31] E0 degradation degree
C     stvMP(12) = stv[32] ft degradation degree
C     stvMP(13) = stv[33] fs degradation degree
C     stvMP(14) = stv[34] fc degradation degree
C     stvMP(15) = stv[35] sf0 degradation degree

C     Imposed stresses [H.Y. moved and added Jan 13 2022]
C     stvMP(16) Imposed Normal Stress
C     stvMP(17) Imposed M Shear Stress
C     stvMP(18) Imposed L Shear Stress

C     For verification only
C     stvMP(19) Humidity check
C     stvMP(20) Humidity rate check
C     stvMP(21) Temperature check
C     stvMP(22) Temperature rate check
C     stvMP(23) Pressure check
C     stvMP(24) Pressure rate check

C     stvCreep components [H.Y. added Jan 13 2022]
C     stvMP(25) Creep Microprestress
C     stvMP(25+1) -> stvMP(25+nKelvinUnit) Creep dGammaN for Kelvin units
C     stvMP(25+nKelvinUnit+1) -> stvMP(25+2*nKelvinUnit) Creep dGammaM for Kelvin units
C     stvMP(25+2*nKelvinUnit+1) -> stvMP(25+3*nKelvinUnit) Creep dGammaL for Kelvin units

      ! use ModuleCoupling
      use ModuleMPFlags
	  
      implicit none

!     Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      ! integer,  parameter :: dp = selected_real_kind(15, 307)
        
!     Indices of imported Multi-Physics fields (MPfields) (degrees of freedom in other analyses)
      integer(kind=4),  parameter           :: humidityID     = 1
      integer(kind=4),  parameter           :: humidityrateID = 2
      integer(kind=4),  parameter           :: tempID         = 3
      integer(kind=4),  parameter           :: temprateID     = 4
      integer(kind=4),  parameter           :: pressureID     = 5
      integer(kind=4),  parameter           :: pressurerateID = 6
      integer(kind=4),  parameter           :: nMPfields      = 6 ! total number of imported multiphysics fields

      ! Indices of multiphysics state variables (stvMP)
      integer(kind=4),  parameter           :: Imp_epsNID     = 1
      integer(kind=4),  parameter           :: Imp_depsNID    = 2
      integer(kind=4),  parameter           :: Imp_depsMID    = 3
      integer(kind=4),  parameter           :: Imp_depsLID    = 4
      integer(kind=4),  parameter           :: Imp_epsVID     = 5
      integer(kind=4),  parameter           :: Imp_depsVID    = 6
      integer(kind=4),  parameter           :: alpha_cID      = 7
      integer(kind=4),  parameter           :: alpha_sID      = 8
      integer(kind=4),  parameter           :: total_alphaID  = 9
      integer(kind=4),  parameter           :: lambdaID       = 10
      integer(kind=4),  parameter           :: HighTempstv1ID = 11
      integer(kind=4),  parameter           :: HighTempstv2ID = 12
      integer(kind=4),  parameter           :: HighTempstv3ID = 13
      integer(kind=4),  parameter           :: HighTempstv4ID = 14
      integer(kind=4),  parameter           :: HighTempstv5ID = 15
      integer(kind=4),  parameter           :: Imp_sigNID     = 16
      integer(kind=4),  parameter           :: Imp_sigMID     = 17
      integer(kind=4),  parameter           :: Imp_sigLID     = 18
      integer(kind=4),  parameter           :: humiditycheckID     = 19
      integer(kind=4),  parameter           :: humidityratecheckID = 20
      integer(kind=4),  parameter           :: tempcheckID         = 21
      integer(kind=4),  parameter           :: tempratecheckID     = 22
      integer(kind=4),  parameter           :: pressurecheckID     = 23
      integer(kind=4),  parameter           :: pressureratecheckID = 24
      integer(kind=4),  parameter           :: CreepMicroprestressID = 25

      integer(kind=4),  parameter           :: nstvMP         = 55 ! total number of fields

!     Indices of Multi-Physics flags (MPflags)
      integer(kind=4),  parameter           :: TimeScalingFlagID      = 1
      integer(kind=4),  parameter           :: ThermalStrainFlagID    = 2
      integer(kind=4),  parameter           :: ShrinkageStrainFlagID  = 3
      integer(kind=4),  parameter           :: ImposedStressFlagID    = 4
      integer(kind=4),  parameter           :: HTDegradModelFlagID    = 5
      integer(kind=4),  parameter           :: HydrationFlagID        = 6
      integer(kind=4),  parameter           :: AgingModelFlagID       = 7
      integer(kind=4),  parameter           :: CreepModelFlagID       = 8
      integer(kind=4),  parameter           :: EvapWaterFlagID        = 9     
      integer(kind=4),  parameter           :: ASRFlagID              = 10
      integer(kind=4),  parameter           :: VolumetricStrainFlagID = 11
      integer(kind=4),  parameter           :: nMPflag                = 11 ! total number of flags

!     Slices of input properties for every model ==================================
      ! the internal variables will be defined in itself subourtine
      integer(kind=4),  parameter                 :: i_props_TimeScalingFlag = 70
      integer(kind=4),  parameter                 :: nprops_TimeScaling      = 11
        
!      integer(kind=4),  parameter               :: i_props_ThermalStrainFlag = 81
!      integer(kind=4),  parameter               :: nprops_ThermalStrain      = 5
        
!      integer(kind=4),  parameter               :: i_props_ShrinkageStrainFlag = 86
!      integer(kind=4),  parameter               :: nprops_ShrinkageStrain      = 5
      integer(kind=4),  parameter               :: i_props_ImposedStressFlag = 91
      integer(kind=4),  parameter               :: nprops_ImposedStress      = 10
        
!      integer(kind=4),  parameter               :: i_props_HTDegradation = 101
!      integer(kind=4),  parameter               :: nprops_HighTemp      = 30

!      integer(kind=4),  parameter               :: i_props_HydrationFlag = 130
!      integer(kind=4),  parameter               :: nprops_Hydration      = 30

      integer(kind=4),  parameter               :: i_props_AgingFlag = 161
      integer(kind=4),  parameter               :: nprops_Aging = 30

      integer(kind=4),  parameter               :: i_props_CreepFlag = 191
      integer(kind=4),  parameter               :: nprops_Creep      = 20

!      integer(kind=4),  parameter               :: i_props_EvapWaterFlag = 211
!      integer(kind=4),  parameter               :: nprops_EvapWater      = 30

!      integer(kind=4),  parameter               :: i_props_ASRFlag = 241
!      integer(kind=4),  parameter               :: nprops_ASR      = 30
!     Slices of input properties for every model ==================================

!     Declare integer variables
      integer(kind=4)                       :: edgeNodeID
      integer(kind=4)                       :: i
      integer(kind=4)                       :: iface
      integer(kind=4)                       :: isvinc
      integer(kind=4)                       :: kinc
      integer(kind=4)                       :: kstep ! current step ID
      integer(kind=4)                       :: mcrd
      integer(kind=4)                       :: nblock
      integer(kind=4)                       :: njprops
      integer(kind=4)                       :: nKelvinUnit
      integer(kind=4)                       :: nnode
      integer(kind=4)                       :: nprops
      integer(kind=4)                       :: nfacet_stvCreep
      integer(kind=4)                       :: nsvars
      integer(kind=4)                       :: nsvint
      integer(kind=4)                       :: stepID ! the step ID of the current step and the steps before
      integer(kind=4)                       :: StepPropsID ! for step ID loop, the props ID
      integer(kind=4)                       :: tetID
      integer(kind=4), dimension(njprops)   :: jprops
      integer(kind=4)                       :: iparticle
      integer(kind=4), dimension(4)         :: nconnec ! particle ID location
	  
!     Declare real variables
      real(dp)                              :: aE
      real(dp)                              :: nE
      real(dp), dimension(12)               :: CreepE0bar
      real(dp), dimension(2)                :: CurrentScaledTime ! scaled time (total time)
      real(dp)                              :: dScaledtime ! scaled dtime
      real(dp)                              :: dtimeCur ! dtime
      real(dp), dimension(12)               :: E0
      real(dp), dimension(mcrd,nnode)       :: elcoord
      real(dp), dimension(mcrd,nnode)       :: elcoord0
      real(dp), dimension(12,30,1)          :: facetgeom ! facet constant geometry properties
      real(dp), dimension(1,12,1)           :: facetgeomMP ! facet constant geometry properties for multiphysics
      real(dp), dimension(12,nMPfields)     :: MPfields ! facet multiphysics fields ! row: facet; column: field
      real(dp)                              :: period
      real(dp), dimension(nprops)           :: props
      real(dp), dimension(2)                :: simulationTime ! time(:) in the main subroutine
      real(dp), dimension(5)                :: TimeScalingFac !  Time scaling factor, 5 is the maximum time step
      real(dp), dimension(12,nstvMP)        :: stvMP !  facet multiphysics state variables
      real(dp), dimension(nsvars)           :: svars_elem
      real(dp), dimension(12)               :: va
      real(dp)                              :: Vol0
      real(dp), dimension(36)               :: AgingHist
      real(dp), dimension(iparticle,5)      :: AggSize
      real(dp), dimension(48)               :: ASRHistTet
      real(dp), dimension(12,nMPfields)     :: ImportedMPfields
      real(dp), dimension(12)               :: stressN
      real(dp), dimension(12)               :: stressM ! M.A. Added for creep
      real(dp), dimension(12)               :: stressL ! M.A. Added for creep

        
      real(dp), parameter :: dzero=0.d0,dmone=-1.0d0,done=1.0d0,dfour=4.0d0,
     1                       dsix=6.d0,dhalf=0.5d0,dfourth=0.25d0

      ! Declare real parameters
      logical(kind=4)                       ::  LargeDefFlag
		
C======================================================================C
C     Initialzation                                                    C
C======================================================================C
      do iface=1,12
        isvinc = (iface-1)*nsvint
        do i=21,nsvint
            stvMP(iface,i-20) = svars_elem(isvinc+i)
        end do
      end do

C======================================================================C
C     Interpolation of imported temperature and humidity               C
C======================================================================C
      call MP2interpolater(tetID,kinc,kstep,mcrd,nblock,nnode,nprops,props,
	1				njprops,jprops,nsvars,svars_elem,elcoord0,facetgeom,facetgeomMP,
	2				period,simulationTime,dtimeCur,nsvint,
	3				nMPfields,MPfields,nMPflag,dScaledtime,CurrentScaledTime,MPflags,
	4				i_props_TimeScalingFlag)

C======================================================================C
C     Aging model                                                      C
C======================================================================C
      
      if (MPflags(AgingModelFlagID) == 1) then
	    ! [H.Y. added the warning for aging+creep cases]
		if (MPflags(CreepModelFlagID) == 1) then
            write(*,*) "Warning: Aging flag and creep flag can not be turned on at the same time, 
	1					Aging + creep has not been implemented yet."
            call XPLB_EXIT ! call exit
		end if
		
        call Aging(nprops,props,nMPfields,ImportedMPfields,
     1            TetID,tempID,humidityID,dScaledtime, ! kinc,
     2            alpha_cID,alpha_sID,total_alphaID,lambdaID,
     3            simulationTime(2),AgingHist,nstvMP,stvMP)
      ! update alpha_c, alpha_s and lambda in [stvMP] for this element using fields h & Tensile
      end if

C======================================================================C
C     HT Degradation Model                                             C
C======================================================================C

      if (MPflags(HTDegradModelFlagID) == 1) then
        call ThermalDegradation(nprops,props,nMPfields,ImportedMPfields, ! nblock,
     1              TetID,tempID,dScaledtime, ! kinc,humidityID
     2              simulationTime(2),nstvMP,stvMP)! output
      ! update alpha_c, alpha_s and lambda in [stvMP] for this element using fields h & Tensile
	  
        if (MPflags(CreepModelFlagID) == 1) then ! H.Y. need to consider the hightemp+creep cases
            write(*,*) "Warning: HighTemp flag and creep flag can not be turned on at the same time, 
	1					hightemp + creep has not been implemented yet."
            call XPLB_EXIT ! call exit
        end if
      end if       
	  
C======================================================================C
C     Creep model                                                      C
C======================================================================C

      if (MPflags(CreepModelFlagID) == 1) then
        MPflags(VolumetricStrainFlagID) = 1

        ! Update elastic modulus, props(2) is E0_inf [H.Y. commented and added the warning for aging+creep cases]
        if (MPflags(AgingModelFlagID) == 1) then
            write(*,*) "Warning: Aging flag and creep flag can not be turned on 
	1					at the same time, Aging + creep has not been implemented yet."
            call XPLB_EXIT ! call exit
        end if
	  
        if (MPflags(HTDegradModelFlagID) == 1) then ! H.Y. need to consider the hightemp+creep cases
			do iface = 1,12
                E0(iface) = E0(iface) * stvMP(iface,HighTempstv1ID)
			end do
        end if
		
		! Initialzation for creep variables
        do iface = 1,12
            stressN(iface) = svars_elem(nsvint*(iface-1)+1)
            stressM(iface) = svars_elem(nsvint*(iface-1)+2)
            stressL(iface) = svars_elem(nsvint*(iface-1)+3)
        end do
		
        if (CurrentScaledTime(2) >= props(198)) then 
            ! Compute va for each facet as a vector, use it to update the instantanious E0 using Eqn 36 in Alnaggar et al 2017, then pass it to the creep function to save computation
            do iface = 1,12
                va(iface) = stvMP(iface,total_alphaID)**props(194) ! stvMP(iface,total_alphaID) is aff, prop(194) is nrd
				CreepE0bar(iface) = done/(props(192) + props(202)/va(iface)) ! props(192) is q1, prop(202) is At0
            end do
		
            call ImposedCreepStrain(nprops,props,nMPfields,ImportedMPfields,
     1            TetID,tempID,humidityID,dScaledtime,CurrentScaledTime,
     2            alpha_cID,alpha_sID,total_alphaID,lambdaID,va,CreepE0bar,
     3            nstvMP,stvMP,stressN,stressM,stressL)
        end if
		
      end if

C======================================================================C
C     Imposed Thermal Strain Subroutines                               C
C======================================================================C
        
      if (MPflags(ThermalStrainFlagID) == 1) then
        MPflags(VolumetricStrainFlagID) = 1
        do iface = 1,12
            stressN(iface) = svars_elem(nsvint*(iface-1)+1)
            stressM(iface) = svars_elem(nsvint*(iface-1)+2)
            stressL(iface) = svars_elem(nsvint*(iface-1)+3)
        end do

        call ImposedThermalStrain(nprops,props,nMPfields,ImportedMPfields,
     1                  tempID,nstvMP,stvMP,stressN,stressM,stressL)
      end if

C======================================================================C
C     Imposed Volumetric Strain Subroutines                            C
C======================================================================C

      if (MPflags(VolumetricStrainFlagID) == 1) then
        call VolumetricImposedStrain(nprops,props,mcrd,nnode,elcoord0,Vol0,nMPfields,MPfields,nstvMP,stvMP)
      end if

C======================================================================C
C     Imposed Stress Subroutines                                       C
C======================================================================C

      if (MPflags(ImposedStressFlagID) == 1) then
        call ImposedStress(nprops_ImposedStress,
	1				props(i_props_ImposedStressFlag:(i_props_ImposedStressFlag+nprops_ImposedStress-1)),
	2				nMPfields,MPfields,pressureID,nstvMP,stvMP)
      end if

C======================================================================C
C     Update Multiphysics State Variables to Element Svars             C
C======================================================================C

	  ! store MP fields in stvMP only for checking  HY changed Jun/20/2021
	  stvMP(1:12,19:24) = MPfields(1:12,1:6)

      ! Update svars_elem with stvMP
      do iface = 1,12
        isvinc = (iface-1)*nsvint 
        stvMP(iface,1) = stvMP(iface,1) + stvMP(iface,2) !impose totoal normal strain, but how about shear strains
        do i=21,nsvint ! do i = (20+1),(20+nstvMP)                          
            svars_elem(isvinc+i) = stvMP(iface,i-20) ! 20 is number of svars used for plain LDPM
        end do
      end do
        
      return
      end subroutine MP2coupler
        
C**********************************************************************C
C                                                                      C
C   Initialzation of the multiphysics variables from external files    C
C                                                                      C
C**********************************************************************C      

      subroutine MP2initializer(tetID,facetgeom,facetgeomMP,mcrd,nblock,nMPflag,
     1                         njprops,nprops,jprops,props,kstep,simulationTime,
     2                         nsvars,svars_elem,MPflags)
! 1  -  Read tet ID & tet node/point ID & their relation/mapping
! 2  -  Read tet node/point ID & field results -> interpolation [time interpolation] ???
!           (time, temperature, humidity, pressure)
! 3  -  Update the TemperatureOld & HumidityOld and Temperature & Humidity
! 4  -  Read the history field data as the initial value from svars_elem(kblock,1+svarsNumber*12)
!           (time, alpha_c, alpha_s, lambda - at the previous step)

      implicit none
        
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)
        
!     Declare integer variables
      integer(kind=4), intent(in)           :: kstep ! current step ID
      integer(kind=4), intent(in)           :: mcrd
      integer(kind=4), intent(in)           :: nblock
      integer(kind=4), intent(in)           :: njprops
      integer(kind=4), intent(in)           :: nMPflag
      integer(kind=4), intent(in)           :: nprops
      integer(kind=4), intent(in)           :: nsvars
      integer(kind=4), intent(in)           :: tetID
      integer(kind=4), dimension(njprops), intent(in) :: jprops
        
!     Declare real variables
      real(dp), intent(in)                  :: facetgeom(12,30,1) ! facet constant geometry properties
      real(dp), intent(in)                  :: facetgeomMP(1,12,1) ! facet constant geometry properties for multiphysics

      real(dp), intent(in)                  :: props(nprops)
      real(dp), intent(in)                  :: simulationTime(2) ! time(2) in the main subroutine
      real(dp), intent(in)                  :: svars_elem(nsvars)
      real(dp), intent(out)                 :: MPflags(nMPflag)
      real(dp), parameter                   :: dzero=0.d0
        
!     Declare booleans (now they are double, discuss to fix)
      real(dp)                              :: TimeScalingFlag
      real(dp)                              :: ThermalStrainFlag
      real(dp)                              :: ShrinkageStrainFlag
      real(dp)                              :: ImposedStressFlag
      real(dp)                              :: HTDegradModelFlag
      real(dp)                              :: HydrationFlag
      real(dp)                              :: AgingModelFlag
      real(dp)                              :: CreepModelFlag
      real(dp)                              :: EvapWaterFlag
      real(dp)                              :: AlkaliSilicaReactionFlag
      real(dp)                              :: VolumetricStrainFlag
        
!     Read the flags
      TimeScalingFlag          = Props(70)
      ThermalStrainFlag        = Props(81)
      ShrinkageStrainFlag      = Props(86)
      ImposedStressFlag        = Props(91)
      HTDegradModelFlag        = Props(101)
      HydrationFlag            = Props(130)
      AgingModelFlag           = Props(161)
      CreepModelFlag           = Props(191)
      EvapWaterFlag            = Props(211)
      AlkaliSilicaReactionFlag = Props(241)
      VolumetricStrainFlag     = dzero
        
      MPflags = dzero
      MPflags(1) = TimeScalingFlag
      MPflags(2) = ThermalStrainFlag
      MPflags(3) = ShrinkageStrainFlag 
      MPflags(4) = ImposedStressFlag
      MPflags(5) = HTDegradModelFlag
      MPflags(6) = HydrationFlag
      MPflags(7) = AgingModelFlag
      MPflags(8) = CreepModelFlag
      MPflags(9) = EvapWaterFlag
      MPflags(10) = AlkaliSilicaReactionFlag
      MPflags(11) = VolumetricStrainFlag

      end subroutine MP2initializer

      subroutine MP2initializer_new(props,nprops,jprops,njprops)
!     H.Y added this subroutine to replace subroutine MPinitializer

      use ModuleMPFlags
      implicit none

      ! Declare integer variables
      integer(kind=4), intent(in)           :: njprops
      integer(kind=4), intent(in)           :: nprops
      integer(kind=4), dimension(njprops), intent(in) :: jprops
        
      ! Declare real variables
      real(dp), intent(in)                  :: props(nprops)
      real(dp), parameter                   :: dzero = 0.0d0

      ! Declare booleans (now they are double, discuss to fix)
      real(dp)                              :: TimeScalingFlag
      real(dp)                              :: ThermalStrainFlag
      real(dp)                              :: ShrinkageStrainFlag
      real(dp)                              :: ImposedStressFlag
      real(dp)                              :: HTDegradModelFlag
      real(dp)                              :: HydrationFlag
      real(dp)                              :: AgingModelFlag
      real(dp)                              :: CreepModelFlag
      real(dp)                              :: EvapWaterFlag
      real(dp)                              :: AlkaliSilicaReactionFlag
      real(dp)                              :: VolumetricStrainFlag
        
      ! Read the flags
      TimeScalingFlag          = Props(70)
      ThermalStrainFlag        = Props(81)
      ShrinkageStrainFlag      = Props(86)
      ImposedStressFlag        = Props(91)
      HTDegradModelFlag        = Props(101)
      HydrationFlag            = Props(130)
      AgingModelFlag           = Props(161)
      CreepModelFlag           = Props(191)
      EvapWaterFlag            = Props(211)
      AlkaliSilicaReactionFlag = Props(241)
      VolumetricStrainFlag     = dzero
        
      MPflags = dzero
      MPflags(1) = TimeScalingFlag
      MPflags(2) = ThermalStrainFlag
      MPflags(3) = ShrinkageStrainFlag 
      MPflags(4) = ImposedStressFlag
      MPflags(5) = HTDegradModelFlag
      MPflags(6) = HydrationFlag
      MPflags(7) = AgingModelFlag
      MPflags(8) = CreepModelFlag
      MPflags(9) = EvapWaterFlag
      MPflags(10) = AlkaliSilicaReactionFlag
      MPflags(11) = VolumetricStrainFlag

      end subroutine MP2initializer_new

C**********************************************************************C
C                                                                      C
C   Spatial and temporal interpolations of the multiphysics fields     C
C                                                                      C
C**********************************************************************C

      subroutine MP2interpolater(tetID,kinc,kstep,mcrd,nblock,nnode,nprops,props,
     1                    njprops,jprops,nsvars,svars_elem,elcoord0,facetgeom,facetgeomMP,
     2                    period,simulationTime,dtimeCur,nsvint,
     3                    nMPfields,MPfields,nMPflag,dScaledtime,CurrentScaledTime,MPflags,
     4                    i_props_TimeScalingFlag)
	 
	  use ModuleCoupling

      implicit none
	  
      real(dp), parameter :: dzero = 0.d0
	  
!     Declare integer variables
      integer(kind=4)                       :: i,j
      integer(kind=4)                       :: ios
      integer(kind=4)                       :: edgeNodeID
      integer(kind=4)                       :: iface
      integer(kind=4)                       :: i_props_TimeScalingFlag
      integer(kind=4)                       :: kinc
      integer(kind=4)                       :: kstep ! current step ID
      integer(KIND=4)                       :: LENOUTDIR
      integer(kind=4)                       :: mcrd
      integer(kind=4)                       :: nblock
      integer(kind=4)                       :: njprops
      integer(kind=4)                       :: nMPfields
      integer(kind=4)                       :: nMPflag
      integer(kind=4)                       :: nnode
      integer(kind=4)                       :: nnodeFLM
      integer(kind=4)                       :: nprops
      integer(kind=4)                       :: nsvars
      integer(kind=4)                       :: nsvint
      integer(kind=4)                       :: stepID ! the step ID of the current step and the steps before
      integer(kind=4)                       :: StepNew
      integer(kind=4)                       :: StepOld
      integer(kind=4)                       :: StepPropsID ! for step ID loop, the props ID
      integer(kind=4)                       :: tetID
      integer(kind=4), dimension(njprops)   :: jprops
        
!     Declare real variables
      real(dp), dimension(2)                :: CurrentScaledTime ! scaled time (real time)
      real(dp)                              :: dScaledtime ! scaled dtime
      real(dp)                              :: dtimeCur ! dtime
      real(dp), dimension(mcrd,nnode)       :: elcoord0
      real(dp), dimension(12,30,1)          :: facetgeom ! facet constant geometry properties
      real(dp), dimension(1,12,1)           :: facetgeomMP ! facet constant geometry properties for multiphysics
	  real(dp)                              :: HNew
	  real(dp)                              :: HOld
	  real(dp), dimension(12,nMPfields)     :: MPfields ! facet multiphysics fields ! line: facet; column: field
      real(dp)                              :: period
      real(dp)                              :: PressureNew
      real(dp)                              :: PressureOld
      real(dp), dimension(nprops)           :: props
      real(dp), dimension(2)                :: simulationTime ! time(2) in the main subroutine
      real(dp), dimension(kstep)            :: TerminalTime !  Time scaling factor, 5 is the maximum time step
      real(dp)                              :: TimeNew
      real(dp)                              :: TimeOld
	  real(dp)                              :: TNew
	  real(dp)                              :: TOld
      real(dp), dimension(kstep)            :: TimeScalingFac !  Time scaling factor, 5 is the maximum time step
      real(dp), dimension(nsvars)           :: svars_elem
      real(dp), dimension(nMPflag)          :: MPflags
	  
!     Declare string variables
      character(256)                        :: OUTDIR  
	  
      ! Declare booleans (now they are double, discuss to fix)
      real(dp)                              :: TimeScalingFlag
      real(dp)                              :: AgingModelFlag
      real(dp)                              :: CreepModelFlag
      real(dp)                              :: ThermalStrainFlag
      real(dp)                              :: ShrinkageStrainFlag
      real(dp)                              :: AlkaliSilicaReactionFlag
      real(dp)                              :: ImposedStressFlag    
      real(dp)                              :: VolumetricStrainFlag
              
C====================================================================C
C                       Temporal interpolation                       C
C====================================================================C

      ! put this before the time interpolation of data pass
      if (MPflags(1) == 1) then
            select case(kstep)
                Case(1)
                    TimeScalingFac(1) = Props(i_props_TimeScalingFlag+1)
                    ! TerminalTime(1)   = Props(i_props_TimeScalingFlag+6)
                    CurrentScaledTime(2) = TimeScalingFac(1)*simulationTime(1)
                Case(2)
                    TimeScalingFac(1) = Props(i_props_TimeScalingFlag+1)
                    TimeScalingFac(2) = Props(i_props_TimeScalingFlag+2)
                    TerminalTime(1)   = Props(i_props_TimeScalingFlag+6)
                    ! TerminalTime(2)   = Props(i_props_TimeScalingFlag+7)
                    CurrentScaledTime(2) = TerminalTime(1)*TimeScalingFac(1)+simulationTime(1)*TimeScalingFac(2)
                Case(3)
                    TimeScalingFac(1) = Props(i_props_TimeScalingFlag+1)
                    TimeScalingFac(2) = Props(i_props_TimeScalingFlag+2)
                    TimeScalingFac(3) = Props(i_props_TimeScalingFlag+3)
                    TerminalTime(1)   = Props(i_props_TimeScalingFlag+6)
                    TerminalTime(2)   = Props(i_props_TimeScalingFlag+7)
                    ! TerminalTime(3)   = Props(i_props_TimeScalingFlag+8)
                    CurrentScaledTime(2) = TerminalTime(1)*TimeScalingFac(1)+TerminalTime(2)*TimeScalingFac(2)+
     1               simulationTime(1)*TimeScalingFac(3)
                Case(4)
                    TimeScalingFac(1) = Props(i_props_TimeScalingFlag+1)
                    TimeScalingFac(2) = Props(i_props_TimeScalingFlag+2)
                    TimeScalingFac(3) = Props(i_props_TimeScalingFlag+3)
                    TimeScalingFac(4) = Props(i_props_TimeScalingFlag+4)
                    TerminalTime(1)   = Props(i_props_TimeScalingFlag+6)
                    TerminalTime(2)   = Props(i_props_TimeScalingFlag+7)
                    TerminalTime(3)   = Props(i_props_TimeScalingFlag+8)
                    ! TerminalTime(4)   = Props(i_props_TimeScalingFlag+9)
                    CurrentScaledTime(2) = TerminalTime(1)*TimeScalingFac(1)+TerminalTime(2)*TimeScalingFac(2)+
     1              TerminalTime(3)*TimeScalingFac(3)+simulationTime(1)*TimeScalingFac(4)
                Case(5)
                    TimeScalingFac(1) = Props(i_props_TimeScalingFlag+1)
                    TimeScalingFac(2) = Props(i_props_TimeScalingFlag+2)
                    TimeScalingFac(3) = Props(i_props_TimeScalingFlag+3)
                    TimeScalingFac(4) = Props(i_props_TimeScalingFlag+4)
                    TimeScalingFac(5) = Props(i_props_TimeScalingFlag+5)
                    TerminalTime(1)   = Props(i_props_TimeScalingFlag+6)
                    TerminalTime(2)   = Props(i_props_TimeScalingFlag+7)
                    TerminalTime(3)   = Props(i_props_TimeScalingFlag+8)
                    TerminalTime(4)   = Props(i_props_TimeScalingFlag+9)
                    CurrentScaledTime(2) = TerminalTime(1)*TimeScalingFac(1)+TerminalTime(2)*TimeScalingFac(2)+
     1              TerminalTime(3)*TimeScalingFac(3)+TerminalTime(4)*TimeScalingFac(4)+simulationTime(1)*TimeScalingFac(5)
            end select
            dScaledtime = TimeScalingFac(kstep)*dtimeCur ! current time step increment
      end if

C======================================================================C
C     Two-way coupling processes
C======================================================================C
      call VGETOUTDIR(OUTDIR,LENOUTDIR)    ! Work directory
	  
	  if (kstep == 0) then ! Abaqus/Explicit Packager stage
		! if (Lop /= 0) then ! Second call of VUEL subroutine in Abaqus/Explicit Packager stage
			! continue
		! end if
		continue
	  else ! Abaqus/Explicit Analysis stage
		if (kinc == 0) then ! Initial data exchange settings
			if (tetID == nomaxel) then ! Exchange when loop to the last element
			  ! Sending LDPM info to Abaqus/Standard solver
			  open(V2U,file='LDPM2FLM.pipe',defaultfile=trim(OUTDIR),form='formatted',
	1			   status='old',action='write',access='stream')
			  write(*,*) 'Sending LDPM analysis settings'
			  write(V2U,'(I8)') nomaxel
			  write(V2U,'(ES24.17)') period
			  flush(V2U)
			  close(V2U)
			  
			  ! Receiving FLM info from Abaqus/Standard solver
			  open(U2V,file='FLM2LDPM.pipe',defaultfile=trim(OUTDIR),form='formatted',
	1			   status='old',action='read',access='stream')
			  write(*,*) 'Retriving FLM analysis settings'
			  read(U2V,'(I8)') nnode_FLM
			  write(*,*) "nnode_FLM", nnode_FLM
			  read(U2V,'(I8)') MPtypeFLM
			  write(*,*) "MPtypeFLM", MPtypeFLM

			  close(U2V)
			  
			  if (MPtypeFLM == 1) then
			    nfieldFLM = 2
			  else if (MPtypeFLM == 2) then
			    nfieldFLM = 1
			  end if
			  
			  if (allocated(FLM2LDPM_DATA) == 0) then
				allocate(FLM2LDPM_DATA(nnode_FLM,nfieldFLM),LDPM2FLM_DATA(nomaxel,13))
				allocate(FLM2LDPM_DATA_old(nnode_FLM,nfieldFLM))
				LDPM2FLM_DATA = 0.0d0
				FLM2LDPM_DATA = 0.0d0
				FLM2LDPM_DATA_old = 0.0d0
			  end if
			end if
			
		else if (kinc > 0) then ! regular data exchange
			if (ExchangeFlag == 1) then
			
				! Assign LDPM values to named pipe variables
				LDPM2FLM_data(tetID,1) = svars_elem(16) ! stv(16) is epsV
				do iface = 1,12
					LDPM2FLM_data(tetID,iface+1) = svars_elem((iface-1)*nsvint+13) ! stv(13) is deltaN
				end do
				
				if (tetID == nomaxel) then

				  open(V2U,file='LDPM2FLM.pipe',defaultfile=trim(OUTDIR),form='formatted',
	1			   status='old',action='write',access='stream')
				
				  ! LDPM to FLM data
				  do i=1,nomaxel
					write(V2U,'(ES24.17)') (LDPM2FLM_DATA(i,j), j=1,13)
				  end do
				
				  ! Current LDPM increment time
				  write(V2U,'(ES24.17)') simulationTime(2)
				  flush(V2U)
				  close(V2U)

				  open(U2V,file='FLM2LDPM.pipe',DEFAULTFILE=trim(OUTDIR),form='formatted',
	1			   status='old',action='read',access='stream')
				  do i=1,nnode_FLM
					read(U2V,'(ES24.17)',iostat=ios) (FLM2LDPM_DATA(i,j), j=1,nfieldFLM)
				  end do

				  read(U2V,'(ES24.17)') period_FLM
				  read(U2V,'(ES24.17)') time_FLM
				  close(U2V)
				  
				  write(*,*) "data has exchanged at time", CurrentScaledTime(2)
				  ExchangeFlag = 0
				end if
			else
				if (tetID == nomaxel) then
					
					if (CurrentScaledTime(2) >= CurrentExchangeTime) then
						! Update two-way coupling variables
						time_FLM_old = time_FLM
						FLM2LDPM_DATA_old = FLM2LDPM_DATA
						ExchangeFlag = 1
						exchangeindex = exchangeindex + 1
						if (exchangeindex > nExchangeTime) then
							CurrentExchangeTime = -LOG(0.0d0) ! After the last exchange, set next exchange time to infinity
						else
							CurrentExchangeTime = ExchangeTimeList(exchangeindex)
						end if
					end if
				end if
			end if
		end if
      end if
	  
C ==================================================================== C
C                                                                      C
C                        Spatial interpolation                         C
C                                                                      C
C ==================================================================== C
      ! Add the spatial interpolation for HTC - 06/07/2021 HY
		
      if (MPtypeFLM == 1) then ! HTC analysis
		do iface = 1,12
			edgeNodeID = facetgeomMP(1,iface,1)
			
			if (kinc > 1) then
				if (StepOld == 0) then
					TimeOld = dzero
					HOld = dzero
					TOld = dzero
				else
					HOld = FLM2LDPM_DATA_old(edgeNodeID,1)
					TOld = FLM2LDPM_DATA_old(edgeNodeID,2)
				end if
				HNew = FLM2LDPM_DATA(edgeNodeID,1)
				TNew = FLM2LDPM_DATA(edgeNodeID,2)
				
				! Interpolated value
				MPfields(iface,1) = HNew
				MPfields(iface,3) = TNew
				! Rate value
				MPfields(iface,2) = dzero
				MPfields(iface,4) = dzero
			else
				MPfields(iface,1) = dzero
				MPfields(iface,2) = dzero
				MPfields(iface,3) = dzero
				MPfields(iface,4) = dzero
			end if
		end do
      elseif (MPtypeFLM == 2) then ! Pore-pressure analysis
		do iface = 1,12
			edgeNodeID = facetgeomMP(1,iface,1)
			if (kinc > 1) then
				if (StepOld == 0) then
					TimeOld = dzero
					PressureOld = dzero
				else
					PressureOld = FLM2LDPM_DATA_old(edgeNodeID,1)
				end if
				PressureNew = FLM2LDPM_DATA(edgeNodeID,1)
				
				MPfields(iface,5) = PressureNew
				MPfields(iface,6) = dzero
			else
				MPfields(iface,5) = dzero
				MPfields(iface,6) = dzero
			end if
		end do
      end if

      return
      end subroutine MP2interpolater
        
        
      subroutine ImposedHygroStrain(ShrinkageCoeff,nMPfields,MPfields,nstvMP,stvMP,dScaledtime)

      implicit none
        
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)
        
      integer(kind=4)                         :: iface
      integer(kind=4)                         :: nMPfields
      integer(kind=4)                         :: nstvMP
        
      real(dp)                                :: dScaledtime
      real(dp)                                :: ShrinkageCoeff
      real(dp), dimension(12,nMPfields)       :: MPfields
      real(dp), dimension(12,nstvMP)          :: stvMP
        
      do iface = 1,12
            stvMP(iface,2) = stvMP(iface,2) + ShrinkageCoeff*dScaledtime*MPfields(iface,2)
      end do
        
      return
      end subroutine ImposedHygroStrain


      subroutine VolumetricImposedStrain(nprops,props,mcrd,nnode,elcoord0,Vol0,nMPfields,MPfields,nstvMP,stvMP)

      implicit none
        
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)

      integer(kind=4)                         :: iface
      integer(kind=4)                         :: mcrd
      integer(kind=4)                         :: nnode
      integer(kind=4)                         :: nprops
      integer(kind=4)                         :: nMPfields
      integer(kind=4)                         :: nstvMP
        
      real(dp), dimension(mcrd,nnode)         :: elcoord0
      real(dp)                                :: Vol0
      real(dp), dimension(12,nMPfields)       :: MPfields
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nstvMP)          :: stvMP
      real(dp)                                :: l12
      real(dp)                                :: l13
      real(dp)                                :: l14
      real(dp)                                :: l23
      real(dp)                                :: l24
      real(dp)                                :: l34
      real(dp)                                :: Volnew
      real(dp)                                :: epsV_imp
        
      l12 = DSQRT((elcoord0(1,1)-elcoord0(1,2))**2+(elcoord0(2,1)-elcoord0(2,2))**2
     1               +(elcoord0(3,1)-elcoord0(3,2))**2)*(1.d0+(stvMP(1,1)+stvMP(2,1))/2.d0)
      l13 = DSQRT((elcoord0(1,1)-elcoord0(1,3))**2+(elcoord0(2,1)-elcoord0(2,3))**2
     1                 +(elcoord0(3,1)-elcoord0(3,3))**2)*(1.d0+(stvMP(3,1)+stvMP(4,1))/2.d0)
      l14 = DSQRT((elcoord0(1,1)-elcoord0(1,4))**2+(elcoord0(2,1)-elcoord0(2,4))**2
     1                 +(elcoord0(3,1)-elcoord0(3,4))**2)*(1.d0+(stvMP(5,1)+stvMP(6,1))/2.d0)
      l23 = DSQRT((elcoord0(1,2)-elcoord0(1,3))**2+(elcoord0(2,2)-elcoord0(2,3))**2
     1                 +(elcoord0(3,2)-elcoord0(3,3))**2)*(1.d0+(stvMP(7,1)+stvMP(8,1))/2.d0)
      l24 = DSQRT((elcoord0(1,2)-elcoord0(1,4))**2+(elcoord0(2,2)-elcoord0(2,4))**2
     1                 +(elcoord0(3,2)-elcoord0(3,4))**2)*(1.d0+(stvMP(9,1)+stvMP(10,1))/2.d0)
      l34 = DSQRT((elcoord0(1,3)-elcoord0(1,4))**2+(elcoord0(2,3)-elcoord0(2,4))**2
     1                 +(elcoord0(3,3)-elcoord0(3,4))**2)*(1.d0+(stvMP(11,1)+stvMP(12,1))/2.d0)

      Volnew = DSQRT((- l12**4*l34**2 - l12**2*l13**2*l23**2 + l12**2*l13**2*l24**2  
     1                + l12**2*l13**2*l34**2 + l12**2*l14**2*l23**2 - l12**2*l14**2*l24**2 
     2                + l12**2*l14**2*l34**2 + l12**2*l23**2*l34**2 + l12**2*l24**2*l34**2 
     3                - l12**2*l34**4 - l13**4*l24**2 + l13**2*l14**2*l23**2 + l13**2*l14**2*l24**2 
     4                - l13**2*l14**2*l34**2 + l13**2*l23**2*l24**2 - l13**2*l24**4
     5                + l13**2*l24**2*l34**2 - l14**4*l23**2 - l14**2*l23**4 + l14**2*l23**2*l24**2 
     6                + l14**2*l23**2*l34**2 - l23**2*l24**2*l34**2)/576.d0)
	 
      epsV_imp = (Volnew-Vol0)/(3.d0*Vol0)
	  
      do iface = 1,12
        stvMP(iface,5) = epsV_imp
      end do
        
      return
      end subroutine VolumetricImposedStrain


      subroutine ImposedStress(nprops,props,nMPfields,MPfields,pressureID,nstvMP,stvMP)

      implicit none
        
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)
        
      integer(kind=4)                         :: iface
      integer(kind=4)                         :: nMPfields
      integer(kind=4)                         :: nprops
      integer(kind=4)                         :: nstvMP
      integer(kind=4)                         :: pressureID
        
      real(dp)                                :: BiotCoeff
      real(dp), dimension(12,nMPfields)       :: MPfields
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nstvMP)          :: stvMP


      BiotCoeff = props(2)

      do iface = 1,12
		stvMP(iface,16) = - BiotCoeff*MPfields(iface,5)
      end do

      return
      end subroutine ImposedStress


C ==================================================================== C
C                                                                      C
C                               Aging model                            C
C                                                                      C
C ==================================================================== C

!     Developed by: Lifu Yang 
!     Hohai University (China), 
!     College of Water Conservancy and Hydropower Engineering
!     & Northwestern University (USA), 
!     Departement of Civil & Environmental Engineering (CEE)
!     Email: yang.lifu@hotmail.com

!     Last edit: May/01/2021


C
C***********************************************************C
C***********************************************************C
C
C     This subroutine updates lambda (Aging degree)
C
C***********************************************************C
C***********************************************************C
C

      subroutine Aging(nprops,props,nMPfields,ImportedMPfields, ! nblock,
     1              TetID,tempID,humidityID,dScaledtime, ! kinc,
     2              alpha_cID,alpha_sID,total_alphaID,lambdaID,
     3              TotalSimuTime,AgingHist,nstvMP,stvMP) ! output
      
      ! nblock: total tet number (not needed here)
      ! nstvMP: 15 (stv21-stv35)
      ! nMPfields: total field number
      ! TetID: tet ID, ielem (might be not needed here)
      ! ImportedMPfields: line->facet, column->field
      ! dScaledtime: scaled dtime
      ! CurrentScaledTime: scaled total time
      ! stvMP: line->facet; column->stv
      ! update alpha_c, alpha_s and lambda in svars matrix using fields h,T
      
      ! use ModuleCoupling
      implicit none 
      
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 202)
      integer,  parameter :: dp = selected_real_kind(15, 307) 
        
!       AgingHistTemp(nblock,36), history variables (previous step values)
!       integer(KIND=4),  parameter :: alpha_cOldID  = 1  ! midpoint formula, for iteration & calculation
!       integer(KIND=4),  parameter :: alpha_sOldID  = 2  ! midpoint formula, for iteration & calculation
      integer(KIND=4),  parameter :: lambdaOldID   = 3  ! midpoint formula, for iteration & calculation 

      integer(kind=4),  parameter :: i_props_AgingFlag     = 161
!         integer(kind=4),  parameter :: nprops_Aging          = 30   

      real(dp), parameter:: CtoK = 273.15d0 ! Celsius to Kelvin conversion

      ! passed from the external subroutine
      integer(KIND=4)                         :: nprops
      integer(KIND=4)                         :: TetID ! ielem
      integer(KIND=4)                         :: humidityID ! 1
      integer(KIND=4)                         :: tempID    ! 3
      integer(KIND=4)                         :: alpha_cID ! alpha_c ID=7
      integer(KIND=4)                         :: alpha_sID ! alpha_s ID=8
      integer(KIND=4)                         :: total_alphaID ! alpha_tot ID=9
      integer(KIND=4)                         :: lambdaID ! lambda ID=10
      integer(KIND=4)                         :: nMPfields ! 6
      integer(KIND=4)                         :: nstvMP  ! total field number
      ! integer(KIND=4)                         :: kinc ! time increment ID
      !
      integer(KIND=4)                         :: OldfieldID ! 1-36
      integer(KIND=4)                         :: iface

      ! real(dp), dimension(nblock)             :: dtimeStable
      ! passed from the external subroutine
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nMPfields)       :: ImportedMPfields 
      real(dp), dimension(12,nstvMP)          :: stvMP ! to be updated in this subroutine
      real(dp)                                :: dScaledtime
      real(dp), dimension(36)                 :: AgingHist
      real(dp)                                :: TotalSimuTime ! time(2), total simulation time
      ! internal variables
!        real(dp), dimension(12)                 :: alphac_inf
!        real(dp), dimension(12)                 :: alphas_inf
!        real(dp), dimension(12)                 :: alphac
!        real(dp), dimension(12)                 :: alphas
      real(dp), dimension(12)                 :: alpha_tot
      real(dp), dimension(12)                 :: alpha_totinf
      real(dp), dimension(12)                 :: dalpha_tot
      real(dp), dimension(12)                 :: dlambda
      real(dp), dimension(12)                 :: dalphac
      real(dp), dimension(12)                 :: dalphas
      real(dp)                                :: lambda
      real(dp)                                :: lambda_prev
      real(dp)                                :: lambda_old
      real(dp)                                :: lambda_new
      real(dp)                                :: a0 ! props
      real(dp)                                :: A_lambda ! props
      real(dp)                                :: T0 ! props
      real(dp)                                :: nlam ! props
      real(dp)                                :: lambda_ini ! props
      real(dp)                                :: A_la0
      real(dp)                                :: temp1
      real(dp)                                :: temp2
      real(dp)                                :: Tmax
      real(dp)                                :: dphi
      real(dp)                                :: temperature ! unit: kelvin [K]

C     input aging parameters      
      a0           = props(i_props_AgingFlag+1)  ! Ageing alpha_0
      A_lambda     = props(i_props_AgingFlag+2)  ! Ageing A_lambda
      T0           = props(i_props_AgingFlag+3)  ! reference temperature, unit K
      nlam         = props(i_props_AgingFlag+4)  ! Temperature Exponent
      lambda_ini   = props(i_props_AgingFlag+5)  ! Initial Aging Degree
      
      
      dphi=1.0d0
      Tmax=323.15d0
      
      ! calculate alpha_c & alpha_s
      ! alphac_old, alphas_old and lambda_old was assigned in the initialized part, and read the value from AgingHistTemp
      ! update stvMP in this subroutine (alpha_c & alpha_s)
      call CuringinLDPM(nprops,props,nMPfields,ImportedMPfields,TotalSimuTime,
     1                humidityID,tempID,TetID,dScaledtime,AgingHist, ! kinc,
     2                alpha_cID,alpha_sID,total_alphaID,nstvMP,stvMP,
     3                alpha_tot,alpha_totinf,dalpha_tot) ! output, do not need (alphac,alphas)
      ! calculate alphac, alphas, alpha_tot, alpha_totinf, and dalpha_tot
      
       
      if (TotalSimuTime > 0.0d0) then ! check if need to initialize lambda            
        do iface = 1,12                      
            ! temperature from C to K
            Temperature = CtoK + ImportedMPfields(iface,tempID)
             
            ! the lambda_old ID in AgingHist
            OldfieldID=(iface-1)*3 ! three old fields

            ! aging degree, lambda
            if (alpha_tot(iface) > a0) then ! alpha>a0
                ! aging calculation
                A_la0 = 1.0d0 / (alpha_totinf(iface) - a0)
                temp1 = A_la0 + A_lambda*(alpha_totinf(iface) + a0-2.0d0*alpha_tot(iface))
                temp2 = ((Tmax-Temperature)/(Tmax-T0))**nlam  ! temperature effect
                ! read the value at previous step, AgingHist is a global matrix in modulus
                lambda_old = AgingHist(OldfieldID+lambdaOldID) ! lambda_old at the previous step 
                lambda_new = AgingHist(OldfieldID+lambdaOldID) ! & initialized at t=0.0
                ! calculation
                dlambda = dalpha_tot(iface) * temp1 * temp2 ! dalpha_tot was calculated in CuringinLDPM subroutine
                lambda_new = lambda_old + dScaledtime*dlambda(iface) ! Create new value
                ! central difference formula
                lambda = 0.5d0*(lambda_new + lambda_old) ! calculate lambda
                ! update1
                lambda_old = lambda_new ! update lambda_old
                stvMP(iface,lambdaID) = lambda ! update lambda
                AgingHist(OldfieldID+lambdaOldID) = lambda_old ! update lambda old
            else ! aging model governing rule, for alpha_tot < alpha_0
                lambda_old = 0.0d0 ! update lambda_old
                stvMP(iface,lambdaID) = 0.0d0 ! update lambda
                AgingHist(OldfieldID+lambdaOldID) = lambda_old ! update lambda old
            end if
        end do ! facet loop  
      else ! initialization of lambda
        do iface = 1,12
            OldfieldID = (iface-1)*3 ! three old fields
            lambda_old = lambda_ini ! 0.0d0 ! update lambda_old
            stvMP(iface,lambdaID) = lambda_ini ! 0.0d0 ! update lambda
            AgingHist(OldfieldID+lambdaOldID) = lambda_old ! update lambda old
        end do
      end if ! update lambda on this facet

      return
      end subroutine Aging


C***********************************************************C
C***********************************************************C
C
C    This subroutine updates curing degree
C    using implicit central differences (alpha_c & alpha_s)
C
C***********************************************************C
C***********************************************************C
C    
      
      subroutine CuringinLDPM(nprops,props,nMPfields,ImportedMPfields,
     1                TotalSimuTime,humidityID,tempID,TetID,dScaledtime,AgingHist, ! kinc,
     2                alpha_cID,alpha_sID,total_alphaID,nstvMP,stvMP,
     3                alpha_tot,alpha_totinf,dalpha_tot) ! alphac,alphas,

      ! use ModuleCoupling
      implicit none 

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 202)
      integer,  parameter :: dp = selected_real_kind(15, 307)

      integer(KIND=4),  parameter :: it_max        = 10
      real(dp),         parameter :: tol           = 1.0d-6
!       real(dp),         parameter :: CtoK          = 273.15d0
!       AgingHistTemp(nblock,36), three history variables
      integer(KIND=4),  parameter :: alpha_cOldID  = 1  ! midpoint formula, for iteration & calculation
      integer(KIND=4),  parameter :: alpha_sOldID  = 2  ! midpoint formula, for iteration & calculation
!       integer(KIND=4),  parameter :: lambdaOldID   = 3  ! midpoint formula, for iteration & calculation   
        
      integer(kind=4),  parameter :: i_props_HydrationFlag = 130
!          integer(kind=4),  parameter :: nprops_Hydration      = 30

      real(dp), parameter:: CtoK = 273.15d0 ! Celsius to Kelvin conversion

      ! passed from the external subroutine
      integer(KIND=4)                         :: nprops
      integer(KIND=4)                         :: humidityID   ! 1
      integer(KIND=4)                         :: tempID       ! 2
      integer(KIND=4)                         :: alpha_cID      ! 7, alpha_c
      integer(KIND=4)                         :: alpha_sID      ! 8, alpha_s
      integer(KIND=4)                         :: total_alphaID  ! 9, alpha_tot
      integer(KIND=4)                         :: nMPfields    ! 6, imported environmental fields
      integer(KIND=4)                         :: TetID        ! ielem
      integer(KIND=4)                         :: nstvMP
      ! integer(KIND=4)                         :: kinc
      !
      integer(KIND=4)                         :: OldfieldID   ! 1-36
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: it_num
!       integer(KIND=4)                         :: isvinc2


      ! passed from the external subroutine
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nMPfields)       :: ImportedMPfields
      real(dp), dimension(2)                  :: CurrentScaledTime
      real(dp)                                :: dScaledtime
      real(dp), dimension(12,nstvMP)          :: stvMP
      real(dp), dimension(36)                 :: AgingHist
      real(dp)                                :: TotalSimuTime ! time(2), total simulation time
      ! 
      real(dp), dimension(12)                 :: alpha_tot ! output
      real(dp), dimension(12)                 :: dalpha_tot ! output
      real(dp), dimension(12)                 :: alpha_totinf ! output
      real(dp), dimension(12)                 :: alphac_old
      real(dp), dimension(12)                 :: alphac_new
      real(dp), dimension(12)                 :: alphac_prev
      real(dp), dimension(12)                 :: dalphac
      real(dp), dimension(12)                 :: alphac_inf
      real(dp), dimension(12)                 :: alphas_old
      real(dp), dimension(12)                 :: alphas_new
      real(dp), dimension(12)                 :: alphas_prev
      real(dp), dimension(12)                 :: dalphas
      real(dp), dimension(12)                 :: alphas_inf
      real(dp)                                :: cemct ! props
      real(dp)                                :: sfct  ! props
      real(dp)                                :: Qc ! props
      real(dp)                                :: Qs ! props
      real(dp)                                :: alpsinf_sat ! props
      real(dp)                                :: alpcinf_sat ! props
      real(dp)                                :: alphac_ini ! props
      real(dp)                                :: alphas_ini ! props
      real(dp), dimension(12)                 :: alphac ! internal, output
      real(dp), dimension(12)                 :: alphas ! internal, output       
      real(dp)                                :: dphi
      real(dp)                                :: temp1
      real(dp)                                :: temp2
      real(dp)                                :: temp3
      real(dp)                                :: Temperature  ! unit: kelvin [K]

C     parameters from the input file    [basic HTC]  
      cemct                   =  props(i_props_HydrationFlag+2) ! cement content
      Qc                      =  props(i_props_HydrationFlag+8) ! Cement Hydration Enthalpy
      alphac_ini              =  props(i_props_HydrationFlag+9) ! Initial Hydration Degree
      sfct                    =  props(i_props_HydrationFlag+10) ! silica fume content
      Qs                      =  props(i_props_HydrationFlag+16) ! Silica Fume Reaction Enthalpy
      alphas_ini              =  props(i_props_HydrationFlag+17) ! InitialSFReactionDegree
      ! [ONIX model]
      alpcinf_sat             =  props(i_props_HydrationFlag+18) ! Saturated Infinity alpha_c
      alpsinf_sat             =  props(i_props_HydrationFlag+20) ! Saturated Infinity alpha_s      
       
C     compute the infinity values
      ! it is called out of the facet loop, so it is necessary to do the facet loop in this subroutine
      call InfinityDegree(props,nprops,nMPfields,humidityID,
     1                 ImportedMPfields,alphac_inf,alphas_inf)
      ! calculate alphac_inf and alphas_inf
       
      if (TotalSimuTime > 0.0d0) then
!       cement hydration reaction
        do iface = 1,12
            ! temperature from C to K
            Temperature = CtoK + ImportedMPfields(iface,tempID)

            ! initial deviation between two iterations    
            dphi=1.0d0
       
            OldfieldID = (iface-1)*3 ! three old fields

            alphac_old(iface) = AgingHist(OldfieldID+alpha_cOldID) ! alphac old
            alphac_new(iface) = AgingHist(OldfieldID+alpha_cOldID)
!           Central diff convergence loop
            it_num = 0 ! Current iteration     
            do while((dphi > tol) .and. (it_max > it_num))  
                ! update dot{alpha} by itself, so call the subroutine inside the do while loop
                ! since it is in the facet loop, passing iface (facet ID) in the subroutine is necessary               
                call HydrationIncrementCement(props,nprops,Temperature, ! nMPfields,
     1                           alphac(iface),alphac_inf(iface), ! tempID,iface,ImportedMPfields, ! passing iface
     2                           dalphac(iface)) ! output
                ! Calculation of curing rate (dalphac,dalphas)              
                alphac_new(iface) = alphac_old(iface)+dScaledtime*dalphac(iface) ! Create new value              
                ! boundary limit
                alphac_new(iface) = min(alphac_new(iface),alphac_inf(iface))                             
                ! New curing change
                dphi = abs(alphac_new(iface)-alphac_prev(iface))   

                alphac_prev(iface) = alphac_new(iface) ! Store old value for convergence check
                ! Curing value at midpoint (central diff)
                alphac(iface) = 0.5d0*(alphac_new(iface)+alphac_old(iface))
                it_num = it_num + 1 ! Update on iteration counter
            end do  ! do while loop
            
            ! Storing new curing value after iterations
            alphac_old(iface) = alphac_new(iface)  
            ! a global matrix to save old variables
            AgingHist(OldfieldID+alpha_cOldID) = alphac_old(iface) ! update alphac_old 
            ! update stvMP
            stvMP(iface,alpha_cID) = alphac(iface) ! update alpha_c
        end do ! facet loop       
               
        ! silica fume reaction         
        do iface = 1,12
            ! temperature from C to K
            Temperature = CtoK + ImportedMPfields(iface,tempID)

            ! initial deviation between two iterations    
            dphi = 1.0d0      
            OldfieldID = (iface-1)*3 ! three old fields
            
            alphas_old(iface) = AgingHist(OldfieldID+alpha_sOldID) ! alphas old
            alphas_new(iface) = AgingHist(OldfieldID+alpha_sOldID)
!           Central diff convergence loop
            it_num = 0 ! Current iteration     
            do while((dphi > tol) .and. (it_max > it_num))  
                ! update dot{alpha} by itself, so call the subroutine inside the do while loop
                ! since it is in the facet loop, passing iface (facet ID) in the subroutine is necessary               
                call HydrationIncrementSF(props,nprops,Temperature, ! nMPfields,tempID,
     1                           alphas(iface),alphas_inf(iface), ! passing iface, ! iface,ImportedMPfields,
     2                           dalphas) ! output, alphas_inf,
                ! Calculation of curing rate (dalphac,dalphas)

                ! Update curing degree based on calculated curing rate        
                alphas_new(iface) = alphas_old(iface)+dScaledtime*dalphas(iface) ! Create new value      
                ! boundary limit
                alphas_new(iface) = min(alphas_new(iface),alphas_inf(iface))
                ! New curing change
                dphi = abs(alphas_new(iface)-alphas_prev(iface))

                alphas_prev(iface) = alphas_new(iface) ! Store old value for convergence check
                ! Curing value at midpoint (central diff)
                alphas(iface) = 0.5d0*(alphas_new(iface) + alphas_old(iface))
                it_num = it_num + 1 ! Update on iteration counter
            end do ! do while loop          

            ! Storing new curing value after iterations
            alphas_old(iface) = alphas_new(iface) 
            ! a global matrix to save old variables
            AgingHist(OldfieldID+alpha_sOldID) = alphas_old(iface)
            ! update stvMP
            stvMP(iface,alpha_sID) = alphas(iface)
        end do ! facet loop             

      else ! Initialize alpha_c & alpha_s
        do iface = 1, 12             
            OldfieldID = (iface-1)*3 ! three old fields
            ! might not need these two lines [check]
            alphac(iface) = alphac_ini ! 0.0d0 ! alpha_c
            alphas(iface) = alphas_ini ! 0.0d0 ! alpha_s
            ! a global matrix to save this, update alphac_old
            AgingHist(OldfieldID+alpha_cOldID) = alphac_ini ! 0.0d0 ! alphac_old
            AgingHist(OldfieldID+alpha_sOldID) = alphas_ini ! 0.0d0 ! alphas_old
            dalphac(iface) = 0.0d0
            dalphas(iface) = 0.0d0 
            ! update alpha_c & alpha_s
            stvMP(iface,alpha_cID) = alphac_ini ! 0.0d0 ! alphac(iface) ! update alpha
            stvMP(iface,alpha_sID) = alphas_ini ! 0.0d0 ! alphas(iface)
        end do ! facet loop
      end if  ! finish alpha_c & alpha_s calculations
       

      temp2 = cemct*Qc
      temp3 = sfct*Qs
      temp1 = temp2*alpcinf_sat+temp3*alpsinf_sat

!     Total reaction degree    
      do iface = 1,12
        alpha_tot(iface) = (alphac(iface)*temp2+alphas(iface)*temp3)/temp1
        alpha_totinf(iface) = (alphac_inf(iface)*temp2+alphas_inf(iface)*temp3)/temp1
        dalpha_tot(iface) = (dalphac(iface)*temp2+dalphas(iface)*temp3)/temp1
        stvMP(iface,total_alphaID) = alpha_tot(iface) ! save alpha_tot
      end do  
       

      return
      end subroutine CuringinLDPM


C***********************************************************C
C***********************************************************C
C
C    This subroutine Computes infinity reaction degree
C    
C***********************************************************C
C***********************************************************C
C     
      subroutine InfinityDegree(props,nprops,nMPfields,humidityID,
     1                         ImportedMPfields,alphac_inf,alphas_inf)

      implicit none

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)

      integer(kind=4),  parameter :: i_props_HydrationFlag = 130
!      integer(kind=4),  parameter :: nprops_Hydration      = 30

      ! passed from the external subroutine
      integer(KIND=4)                         :: nprops
      integer(KIND=4)                         :: humidityID
      integer(KIND=4)                         :: nMPfields
      !
      integer(KIND=4)                         :: iface

      ! passed from the external subroutine
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nMPfields)       :: ImportedMPfields
      ! 
      real(dp), dimension(12)                 :: alphac_inf ! output
      real(dp), dimension(12)                 :: alphas_inf ! output
      real(dp)                                :: alpcinf_sat
      real(dp)                                :: zetac
      real(dp)                                :: alpsinf_sat
      real(dp)                                :: zetas

!     Input parameters   [ONIX model]   
      alpcinf_sat         =  props(i_props_HydrationFlag+18) ! SaturatedInfinityalphac
      zetac               =  props(i_props_HydrationFlag+19) ! CementHydrationZeta
      alpsinf_sat         =  props(i_props_HydrationFlag+20) ! SaturatedInfinityalphas
      zetas               =  props(i_props_HydrationFlag+21) ! SilicaFumeReactionZeta

      do iface = 1,12
!     Infinity degree
        alphac_inf(iface) = alpcinf_sat*exp(-zetac*(1.0d0/
     1                    ImportedMPfields(iface,humidityID)-1.0d0)) ! humidity
        alphas_inf(iface) = alpsinf_sat*exp(-zetas*(1.0d0/
     1                    ImportedMPfields(iface,humidityID)-1.0d0)) ! humidity
      end do

      return
      end subroutine InfinityDegree


C***********************************************************C
C***********************************************************C
C                                                           C
C     This subroutine Computes the rate of                  C 
C     cement hydration reaction degree                      C
C                                                           C
C***********************************************************C
C***********************************************************C
C     
      subroutine HydrationIncrementCement(props,nprops, ! nMPfields,
     1                          ImportedTempValue,alphac,alphac_inf,
     2                          dalphac) 


      ! HydrationIncrementCement is already for a certain facet,
      ! so the facet loop is not needed here.

      implicit none

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307)
      
      integer(kind=4),  parameter :: i_props_HydrationFlag = 130
!      integer(kind=4),  parameter :: nprops_Hydration      = 30
      real(dp),  parameter :: R               =  8.314d3

       
      ! passed from the external subroutine
      integer(KIND=4)                         :: nprops
      real(dp), dimension(nprops)             :: props
      real(dp)                                :: ImportedTempValue ! the unit is already converted
      real(dp)                                :: alphac
      real(dp)                                :: alphac_inf
      ! 
      real(dp)                                :: dalphac ! output result
      real(dp)                                :: Temp_refc
      real(dp)                                :: Ac1
      real(dp)                                :: Ac2
      real(dp)                                :: etac
      real(dp)                                :: Eac

!     Input parameters [Basic HTC]     
      Temp_refc         =  props(i_props_HydrationFlag+3) ! HydrationReferenceTemperature
      Ac1               =  props(i_props_HydrationFlag+4) ! CementHydrationAc1
      Ac2               =  props(i_props_HydrationFlag+5) ! CementHydrationAc2
      etac              =  props(i_props_HydrationFlag+6) ! CementHydrationEtac
      Eac               =  props(i_props_HydrationFlag+7) ! CementHydrationEnergy
       
      
      dalphac = Ac1*(Ac2+alphac)*abs(alphac_inf-alphac) * 
     1       exp(-etac*alphac/alphac_inf) * exp(Eac/R/Temp_refc-Eac/ImportedTempValue/R)
       
      return
      end subroutine HydrationIncrementCement

      
C***********************************************************C
C***********************************************************C
C                                                           C
C     This subroutine Computes the rate of                  C 
C     silica fume reaction degree                           C
C                                                           C
C***********************************************************C
C***********************************************************C
C  
      subroutine HydrationIncrementSF(props,nprops,! nMPfields,! tempID,
     1                          ImportedTempValue,alphas,alphas_inf, ! facetID,
     2                          dalphas) ! output

      ! HydrationIncrementCement is already for a certain facet,
      ! so the facet loop is not needed here.

      implicit none

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20-22)
      integer,  parameter :: dp = selected_real_kind(15, 307) 
      
      integer(kind=4),  parameter :: i_props_HydrationFlag = 130
!      integer(kind=4),  parameter :: nprops_Hydration      = 30
      real(dp),  parameter :: R               =  8.314d3
      
      ! passed from the external subroutine
      integer(KIND=4)                         :: nprops
      real(dp), dimension(nprops)             :: props
      real(dp)                                :: alphas
      real(dp)                                :: alphas_inf
      real(dp)                                :: ImportedTempValue ! the unit is already converted
      ! 
      real(dp)                                :: dalphas     
      real(dp)                                :: Temp_refs
      real(dp)                                :: As1
      real(dp)                                :: As2
      real(dp)                                :: etas
      real(dp)                                :: Eas

!     Input parameters    [Basic HTC]  
      Temp_refs         =  props(i_props_HydrationFlag+11) ! SFReactionReferenceTemperature
      As1               =  props(i_props_HydrationFlag+12) ! SilicaFumeReactionAs1
      As2               =  props(i_props_HydrationFlag+13) ! SilicaFumeReactionAs2
      etas              =  props(i_props_HydrationFlag+14) ! SilicaFumeReactionEtac
      Eas               =  props(i_props_HydrationFlag+15) ! SilicaFumeReactionEnergy
     
      dalphas = As1*(As2+alphas)*abs(alphas_inf-alphas)*
     1    exp(-etas*alphas/alphas_inf) * exp(Eas/R/Temp_refs-Eas/ImportedTempValue/R)
      return
      end subroutine HydrationIncrementSF


C***********************************************************C
C***********************************************************C
C                                                           C
C     This subroutine Computes the                          C 
C     thermal degradatiom degree                            C
C     by Lei Shen 20210601                                  C
C                                                           C
C***********************************************************C
C***********************************************************C
C  
      subroutine ThermalDegradation(nprops,props,nMPfields,ImportedMPfields, ! nblock,
     1              TetID,tempID,dScaledtime, ! kinc,humidityID
     2              TotalSimuTime,nstvMP,stvMP) ! output

      ! nblock: total tet number (not needed here)
      ! nstvMP: 15 (stv21-stv35)
      ! nMPfields: total field number
      ! TetID: tet ID, ielem (might be not needed here)
      ! ImportedMPfields: line->facet, column->field
      ! dScaledtime: scaled dtime
      ! CurrentScaledTime: scaled total time
      ! stvMP: line->facet; column->stv
      ! update alpha_c, alpha_s and lambda in svars matrix using fields h,T
      
       ! use ModuleCoupling
      implicit none 
      
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 202)
      integer,  parameter :: dp = selected_real_kind(15, 307) 

      integer(kind=4), parameter :: i_props_HTDegradModelFlag     = 101
      !  real(dp), parameter:: CtoK = 273.15d0 ! Celsius to Kelvin conversion

      ! passed from the external subroutine
      integer(KIND=4)                         :: nprops
      integer(KIND=4)                         :: TetID ! ielem
      integer(KIND=4)                         :: humidityID ! 1
      integer(KIND=4)                         :: tempID     ! 3
      integer(KIND=4)                         :: nMPfields  ! 6
      integer(KIND=4)                         :: nstvMP  ! total field number
      ! integer(KIND=4)                         :: kinc ! time increment ID
      !
      integer(KIND=4)                         :: OldfieldID ! 1-36
      integer(KIND=4)                         :: iface

      ! real(dp), dimension(nblock)             :: dtimeStable
      ! passed from the external subroutine
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nMPfields)       :: ImportedMPfields 
      real(dp), dimension(12,nstvMP)          :: stvMP ! to be updated in this subroutine
      real(dp)                                :: dScaledtime
      real(dp)                                :: TotalSimuTime ! time(2), total simulation time
      ! internal variables
      real(dp)                                :: temperature ! unit: kelvin [K]
      real(dp)                                :: E0_Ts 
      real(dp)                                :: E0_Tm 
      real(dp)                                :: E0_nd 
      real(dp)                                :: ft_Ts 
      real(dp)                                :: ft_Tm 
      real(dp)                                :: ft_nd 
      real(dp)                                :: fs_Ts 
      real(dp)                                :: fs_Tm 
      real(dp)                                :: fs_nd 
      real(dp)                                :: fc_Ts 
      real(dp)                                :: fc_Tm 
      real(dp)                                :: fc_nd 
      real(dp)                                :: sf0_Ts
      real(dp)                                :: sf0_Tm
      real(dp)                                :: sf0_nd
    
C     input ThermalDegradation parameters
C     i_props_HTDegradFlag=101   
      E0_Ts     = props(i_props_HTDegradModelFlag+1)
      E0_Tm     = props(i_props_HTDegradModelFlag+2)
      E0_nd     = props(i_props_HTDegradModelFlag+3)
      ft_Ts     = props(i_props_HTDegradModelFlag+4)
      ft_Tm     = props(i_props_HTDegradModelFlag+5)
      ft_nd     = props(i_props_HTDegradModelFlag+6)
      fs_Ts     = props(i_props_HTDegradModelFlag+7)
      fs_Tm     = props(i_props_HTDegradModelFlag+8)
      fs_nd     = props(i_props_HTDegradModelFlag+9)
      fc_Ts     = props(i_props_HTDegradModelFlag+10)
      fc_Tm     = props(i_props_HTDegradModelFlag+11)
      fc_nd     = props(i_props_HTDegradModelFlag+12)
      sf0_Ts    = props(i_props_HTDegradModelFlag+13)
      sf0_Tm    = props(i_props_HTDegradModelFlag+14)
      sf0_nd    = props(i_props_HTDegradModelFlag+15)
       
      if (TotalSimuTime > 0.0d0) then ! check simulation time for initialization            
        do iface = 1,12                      
            ! temperature from C to K
            Temperature = ImportedMPfields(iface,tempID)
            call ThermalDegradationFunction(Temperature,E0_Ts,E0_Tm,E0_nd,stvMP(iface,11))
            call ThermalDegradationFunction(Temperature,ft_Ts ,ft_Tm ,ft_nd,stvMP(iface,12))
            call ThermalDegradationFunction(Temperature,fs_Ts ,fs_Tm ,fs_nd,stvMP(iface,13))
            call ThermalDegradationFunction(Temperature,fc_Ts ,fc_Tm ,fc_nd,stvMP(iface,14))
            call ThermalDegradationFunction(Temperature,sf0_Ts,sf0_Tm,sf0_nd,stvMP(iface,15))
        end do ! facet loop  

      else ! initialization of thermal degradation degree
        do iface = 1,12
            stvMP(iface,11) = 1
            stvMP(iface,12) = 1
            stvMP(iface,13) = 1
            stvMP(iface,14) = 1
            stvMP(iface,15) = 1
        end do
      end if ! update thermal degradation degree
	  
      return
      end subroutine ThermalDegradation


      subroutine ThermalDegradationFunction(Tnow,Tstart,Tmelt,HTnd,fd)!avaiable for LDPM-F

      implicit none

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 202)
      integer,  parameter :: dp = selected_real_kind(15, 307) 
      real(dp)                      :: fd
      real(dp)                      :: HTtheta
      real(dp)                      :: Tnow
      real(dp)                      :: Tstart
      real(dp)                      :: Tmelt
      real(dp)                      :: HTnd

      HTtheta = (abs(Tnow-Tstart)+Tnow-Tstart)/2/(Tmelt-Tstart)
      fd      = 1-exp(HTnd)*HTtheta/(1-HTtheta*(1-exp(HTnd)))
      
	  return
      end subroutine ThermalDegradationFunction


C***********************************************************C
C***********************************************************C
C                                                           C
C           This subroutine computes the                    C 
C            thermal deformations at HT & RT                C
C            by Lei Shen 20210710                           C
C                                                           C
C***********************************************************C
C***********************************************************C
C  
      subroutine ImposedThermalStrain(nprops,props,nMPfields,ImportedMPfields,
     1                        tempID,nstvMP,stvMP,stressN,stressM,stressL)
      
	  implicit none
	  
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 202)
      integer,  parameter :: dp = selected_real_kind(15, 307)
	  
      integer(kind=4)                         :: iface
      integer(kind=4)                         :: nMPfields
      integer(kind=4)                         :: nprops
      integer(kind=4)                         :: nstvMP
      integer(kind=4)                         :: tempID     ! 3
      real(dp), dimension(12,nMPfields)       :: ImportedMPfields
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nstvMP)          :: stvMP
        
      real(dp)                                :: ExpansionCoeff
      real(dp)                                :: LITSPoisson
      real(dp)                                :: LITSCoefA
      real(dp)                                :: LITSCoefB
      real(dp)                                :: sigmaN0
      real(dp)                                :: sigmat
      real(dp)                                :: sigmas
      real(dp)                                :: AlphaHT
      real(dp)                                :: BetaHT
      real(dp)                                :: GammaHT
      real(dp), dimension(12)                 :: stressN
      real(dp), dimension(12)                 :: stressM
      real(dp), dimension(12)                 :: stressL
      real(dp)                                :: stressNinC
      real(dp)                                :: stressNinT
      real(dp)                                :: Tlevel

      ExpansionCoeff = props(82) !ThermalStrainFlag=81
      LITSPoisson = props(83)!Not used 20220214
      LITSCoefA = props(84)
      LITSCoefB = props(85)
      sigmaN0 = props(8) !CompressiveYieldingStrength1
      sigmat = props(4) !TensileStrength1
      sigmas = props(6) * props(4)!shear strength
      !GammaHT = 1.0 !not an input value
	  
      do iface = 1,12
        stressNinC = (dabs(stressN(iface))-stressN(iface))/2!stv[1], normal stress in compression
        stressNinT = (dabs(stressN(iface))+stressN(iface))/2!stv[1], normal stress in tension
        Tlevel = dabs(ImportedMPfields(iface,tempID))/100
        if (Tlevel < 6) then
            AlphaHT = ExpansionCoeff*10/(10-Tlevel)
        else
            AlphaHT = 0
        end if
		
        BetaHT = 0.01*(2*Tlevel*LITSCoefA+LITSCoefB)
        stvMP(iface,2) = stvMP(iface,2) + ImportedMPfields(iface,tempID+1)*(AlphaHT
     1                    -BetaHT*stressNinC/max(stressNinC,sigmaN0)
     2                    +BetaHT*Tlevel*Tlevel*stressNinT/sigmat)!Imposed Normal Strain Increment
        stvMP(iface,3) = stvMP(iface,3) + 
	1					 ImportedMPfields(iface,tempID+1)*BetaHT*stressM(iface)/max(stressM(iface),sigmas)!Imposed M Shear Strain Increment
        stvMP(iface,4) = stvMP(iface,4) + 
	2					 ImportedMPfields(iface,tempID+1)*BetaHT*stressL(iface)/max(stressL(iface),sigmas)!Imposed L Shear Strain Increment
      end do

      return
      end subroutine ImposedThermalStrain  

C***********************************************************C
C***********************************************************C
C                                                           C
C     This subroutine computes the                          C
C     imposed creep strains following MPS                   C
C     by Mohammed Alnaggar 20211202                         C
C                                                           C
C***********************************************************C
C***********************************************************C
C
      subroutine ImposedCreepStrain(nprops,props,nMPfields,ImportedMPfields,
	1				TetID,tempID,humidityID,dScaledtime,CurrentScaledTime,
	2				alpha_cID,alpha_sID,total_alphaID,lambdaID,va,CreepE0bar,
	3				nstvMP,stvMP,stressN,stressM,stressL)

      use ModuleCreep
	  implicit none

      integer(kind=4)                         :: i
      integer(kind=4)                         :: iface
      integer(kind=4)                         :: nfacet_stvCreep
      integer(kind=4)                         :: nMPfields
      integer(kind=4)                         :: nprops
      integer(kind=4)                         :: nstvMP
      integer(kind=4)                         :: TetID     
      integer(kind=4)                         :: tempID     ! 3
      integer(kind=4)                         :: humidityID    
      integer(kind=4)                         :: alpha_cID     
      integer(kind=4)                         :: alpha_sID     
      integer(kind=4)                         :: total_alphaID     
      integer(kind=4)                         :: lambdaID  
      integer(kind=4),  parameter             :: CreepMicroprestressID = 25

      real(dp), dimension(12,nMPfields)       :: ImportedMPfields
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(12,nstvMP)          :: stvMP
      real(dp)                                :: dScaledtime
      real(dp), dimension(2)                  :: CurrentScaledTime

      real(dp)                                :: q1
      real(dp)                                :: q2
      real(dp)                                :: q4
      real(dp)                                :: mps0
      real(dp)                                :: K0
      real(dp)                                :: K1
      real(dp)                                :: nrd
      real(dp)                                :: QvR
      real(dp)                                :: QsR

      real(dp)                                :: rdt
      real(dp)                                :: rt
      real(dp)                                :: tmp
      real(dp)                                :: rhm
      real(dp)                                :: tmpr
      real(dp)                                :: rhmr
      real(dp)                                :: aff
	  
! Note that creep should not start until real time equals t0

      real(dp)                                :: depsfN
      real(dp)                                :: depsfM
      real(dp)                                :: depsfL
      real(dp)                                :: depsvN
      real(dp)                                :: depsvM
      real(dp)                                :: depsvL
      real(dp)                                :: alpha

      real(dp)                                :: sgn
      real(dp)                                :: sgm
      real(dp)                                :: sgl
      real(dp)                                :: t0
      real(dp)                                :: dGammaN
      real(dp)                                :: dGammaM
      real(dp)                                :: dGammaL

      real(dp)                                :: psi_i
      real(dp)                                :: psi_s
      real(dp)                                :: dGammaNj
      real(dp)                                :: dGammaMj
      real(dp)                                :: dGammaLj

      real(dp), dimension(12)                 :: va
      real(dp), dimension(12)                 :: CreepE0bar
      real(dp), dimension(12)                 :: stressN
      real(dp), dimension(12)                 :: stressM
      real(dp), dimension(12)                 :: stressL

!     the following variables are already defined in ModuleCreep  
!      real(dp), dimension(11)                 :: At
!      real(dp), dimension(11)                 :: tau
!      real(dp), dimension(12,31)              :: stvCreep !  facet creep state variables
      
	  real(dp), parameter :: time2days = 86400.0d0 ! 1.0*24.0*3600.0 M.A. Assuming seconds

      QvR = 5000.0
      QsR = 3000.0
	  alpha = props(3)
      q1 = props(192) 
      q2 = props(193)
      nrd = props(194)
      K0 = props(195) / time2days 
      K1 = props(196)
      q4 = props(197)
      t0 = props(198)
      mps0 = props(199)
	  
      rdt = dScaledtime ! real time delta t [M.A. must be retrirved] [H.Y. changed Dec 08 2021]
	  rt = CurrentScaledTime(2) ! real current time

      nfacet_stvCreep = 1 + 3*nKelvinUnit
      do iface = 1,12
	    stvCreep(iface,1:nfacet_stvCreep,TetID) = stvMP(iface,CreepMicroprestressID:CreepMicroprestressID+nfacet_stvCreep-1) ! [M.A. Check][M.A. seems ok Feb 05 2022]
      end do
	  
	  if (rt > t0) then
        do iface = 1,12
            tmp = ImportedMPfields(iface,tempID)
            rhm = ImportedMPfields(iface,humidityID)
            tmpr = ImportedMPfields(iface,tempID+1)
            rhmr = ImportedMPfields(iface,humidityID+1)

            dGammaN = 0.0
            dGammaM = 0.0
            dGammaL = 0.0
            psi_i = (0.1+0.9*(rhm**2)) * exp(QvR*(1.0/296.0-1.0/tmp))
            psi_s = (0.1+0.9*(rhm**2)) * exp(QsR*(1.0/296.0-1.0/tmp))

            do i = 2,nKelvinUnit+1 ! M.A. this loop is very problematic why 10+i when the number of kelivin units is varriable. it should depend on nKelvinUnit [H.Y. changed Feb 07 2022]
                dGammaNj = (sgn*CreepAt(i)-stvCreep(iface,i,TetID)) * (1.0-exp(-rdt*psi_i/(CreepTau(i) * time2days)))
                dGammaMj = (sgm/alpha*CreepAt(i)-stvCreep(iface,nKelvinUnit+i,TetID)) 
	1						* (1.0-exp(-rdt*psi_i/(CreepTau(i) * time2days)))
                dGammaLj = (sgl/alpha*CreepAt(i)-stvCreep(iface,2*nKelvinUnit+i,TetID)) 
	1						* (1.0-exp(-rdt*psi_i/(CreepTau(i) * time2days)))
                dGammaN = dGammaN + dGammaNj
                dGammaM = dGammaM + dGammaMj
                dGammaL = dGammaL + dGammaLj
                stvCreep(iface,i,TetID) = stvCreep(iface,i,TetID) + dGammaNj
                stvCreep(iface,nKelvinUnit+i,TetID) = stvCreep(iface,nKelvinUnit+i,TetID) + dGammaMj
                stvCreep(iface,2*nKelvinUnit+i,TetID) = stvCreep(iface,2*nKelvinUnit+i,TetID) + dGammaLj
            end do
		
            depsvN = (dGammaN/va(iface))
            depsvM = (dGammaM/va(iface))
            depsvL = (dGammaL/va(iface))

            depsfN = sgn       * (psi_i*K0*q4*stvCreep(iface,1,TetID))*rdt
            depsfM = sgm/alpha * (psi_i*K0*q4*stvCreep(iface,1,TetID))*rdt
            depsfL = sgl/alpha * (psi_i*K0*q4*stvCreep(iface,1,TetID))*rdt
            stvCreep(iface,1,TetID) = stvCreep(iface,1,TetID) + (-psi_s*K0*stvCreep(iface,1,TetID)*stvCreep(iface,1,TetID)
	1								+ K1*DABS(tmpr*DLOG(rhm)+tmp*rhmr/rhm)) * rdt

            stvMP(iface,2) = stvMP(iface,2) + depsfN ! M.A. shouldn't we add to the previous value? [H.Y. changed Dec 08 2021]
            stvMP(iface,3) = stvMP(iface,3) + depsfM
            stvMP(iface,4) = stvMP(iface,4) + depsfL
        end do
	  else
        do iface = 1,12
            stvCreep(iface,1,TetID) = mps0
        end do
	  end if
	  
!     Update elem stvs using creep stvs
      do iface = 1,12
        stvMP(iface,CreepMicroprestressID:CreepMicroprestressID+nfacet_stvCreep-1) = stvCreep(iface,1:nfacet_stvCreep,TetID)
      end do
		
      return
      end subroutine ImposedCreepStrain
