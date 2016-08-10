!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CO2SAT_mod
!
! !DESCRIPTION: Module CO2SAT\_MOD contains variables and routines which 
!  follow a satellite track through the GEOS-Chem model simulation.  This is 
!  useful for comparing model results with satellite observations. 
!\\
!\\
! !INTERFACE:
!
      MODULE CO2SAT_MOD
!
! !USES:
!
      USE inquireMod,       ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
     
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_CO2SAT
      PUBLIC  :: CO2SAT
      PUBLIC  :: SETUP_CO2SAT
      PUBLIC  :: SET_CO2SAT
      PUBLIC  :: CO2SAT_OUTPUT_CLOSE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: AN_SETUP
      PRIVATE :: INIT_CO2SAT
      PRIVATE :: NOY_SETUP
!      PRIVATE :: READ_VARIABLES
      PRIVATE :: READ_POINTS
      PRIVATE :: RO2_SETUP
      PRIVATE :: TEST_VALID
      PRIVATE :: CO2SAT_OUTPUT_CREATE
      PRIVATE :: CO2SAT_OBS_OUTPUT_WRITE
      PRIVATE :: CO2SAT_MODEL_OUTPUT_WRITE
      PRIVATE :: COUNT_OBS_IN_MONTH

!
! !REMARKS:
!  The quantities that are saved to disk by the CO2SAT diagnostic were
!  taken from gosat_mod.F.  If you would like to save out a new quantity,
!  then you will have to make your own modifications in this module.
!
! !REVISION HISTORY:
!  (1 ) April 2016 Created by R. McKeown from gosat_mod.f as a more generic
!       CO2 satellite module.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! DO_CO2SAT       : Turn on the CO2SAT diagnostic? (T/F)
      ! MAXVARS     : Maximum # of variables allowed
      ! MAXPOINTS   : Maximum # of flight track points allowed
      ! MAXREAC     : Maximum # of SMVGEAR reactions allowed
      ! MAXRO2      : Maximum # of RO2 constituents allowed
      ! NPOINTS     : Number of flight track points 
      ! PPOINT      : Pointer to last measured output
      ! PDATE       : Array of dates     at each flight point
      ! PTIME       : Array of times     at each flight point
      ! PTAU        : Array of TAU's     at each flight point
      ! PLAT        : Array of latitude  at each flight point
      ! PLON        : Array of longitude at each flight point
      ! PPRESS      : Array of pressure  at each flight point
      ! PTYPE       : Array of ID'#S     at each flight point
      ! PPRESSWGHT  : Array of pressure weights at each point
      ! NPVAR       : # of var's to be saved at each flight point
      ! PVAR        : Array of variable indices
      ! PNAME       : Array of variable names corresponding to PVAR
      ! NPREAC      : # of variables that are really SMVGEAR rxns
      ! PREAC       : Array of SMVGEAR rxn index numbers
      ! PRRATE      : Array of rxn rates for each entry in PREAC
      ! NRO2        : # number of RO2 constituents
      ! PRO2        : Array of SMVGEAR species that are RO2 const's
      ! INFILENAME  : Name of input file defining the flight track
      ! OUTFILENAME : Name of output file 
      !=================================================================

      ! Logicals
      LOGICAL                        :: DO_CO2SAT

      ! Parameters
      INTEGER,           PARAMETER   :: MAXVARS   = 200
      INTEGER,           PARAMETER   :: MAXPOINTS = 200000
      INTEGER,           PARAMETER   :: MAXREAC   = 50
      INTEGER,           PARAMETER   :: MAXRO2    = 45
      INTEGER,           PARAMETER   :: MAXAN     = 10 ! FP 04/01/2010 (AN) 
      INTEGER,           PARAMETER   :: MAXNOY    = 12 ! FP 04/01/2010 (NOy)

      ! For specifying flight track points
      INTEGER                        :: NPOINTS           
      INTEGER                        :: PPOINT

      ! For specifying date/time
      INTEGER,           ALLOCATABLE :: PDATE(:)
      INTEGER,           ALLOCATABLE :: PTIME(:)              
      REAL*4,            ALLOCATABLE :: PTAU(:)               

      ! For specifying lat/lon/alt/err and ID type
      REAL*4,            ALLOCATABLE :: PLAT(:)               
      REAL*4,            ALLOCATABLE :: PLON(:)               
      REAL*4,            ALLOCATABLE :: PPRESS(:)
      REAL*4,            ALLOCATABLE :: PXCO2(:)             
      REAL*4,            ALLOCATABLE :: PERR(:)
      CHARACTER(LEN=5),  ALLOCATABLE :: PTYPE(:)              
      CHARACTER(LEN=14), ALLOCATABLE :: PID(:)

      ! For pressure weights
      REAL*8,            ALLOCATABLE :: PPRESSWGHT(:,:)
      REAL*8,            ALLOCATABLE :: PAK(:,:)
      REAL*8,            ALLOCATABLE :: POCO2_PR(:,:)
      REAL*8,            ALLOCATABLE :: POCO2_PRESS(:,:)

      ! For specifying variables to save at each flight point
      INTEGER                        :: NPVAR        
      INTEGER,           ALLOCATABLE :: PVAR(:) 
      CHARACTER(LEN=10), ALLOCATABLE :: PNAME(:)              
      
      ! For specifying SMVGEAR rxns to save at each flight point
      INTEGER                        :: NPREAC        
      INTEGER,           ALLOCATABLE :: PREAC(:) 
      REAL*8,            ALLOCATABLE :: PRRATE(:,:) 

      ! For specifying RO2 constituents at each flight point
      INTEGER                        :: NPRO2
      INTEGER                        :: PRO2(MAXRO2)

      ! For specifying NOY constituents at each flight point
      INTEGER                        :: NPNOY
      INTEGER                        :: PNOY(MAXNOY)

      ! For specifying AN constituents at each flight point
      INTEGER                        :: NPAN
      INTEGER                        :: P_AN(MAXAN)

      ! Vert Integrate Opt
      INTEGER                        :: VERTINTOPT 

      ! Input/output file names
      CHARACTER(LEN=255)             :: INFILENAME,  INF
      CHARACTER(LEN=255)             :: OUTFILENAME, OUTF

      ! Logical unit numbers and file ids (formerly in file_mod.F)
      INTEGER                        :: IU_FILE
      INTEGER                        :: SCO2_NCID = 0
      
      ! Observation number in monthly output
      INTEGER                        :: OBSNUM = 1 

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setup_CO2SAT
!
! !DESCRIPTION: Subroutine SETUP\_CO2SAT reads information from the 
!  input file in order to initialize the CO2SAT diagnostic.  Also 
!  calls INIT\_CO2SAT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETUP_CO2SAT ( am_I_Root , Input_Opt,RC)
!
! !USES:
!
      USE FILE_MOD,   ONLY : FILE_EXISTS
      USE FILE_MOD,   ONLY : IOERROR
      USE FILE_MOD,   ONLY : IU_FILE
      USE TIME_MOD,   ONLY : EXPAND_DATE
      USE TIME_MOD,   ONLY : GET_NYMD
      USE TIME_MOD,   ONLY : GET_NHMS
      USE TIME_MOD,   ONLY : ITS_A_NEW_MONTH
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE m_netcdf_io_close
      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_get_dimlen
      USE m_do_err_out
      USE m_netcdf_io_checks
      
      include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?

! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! !REMARKS:
!  For SMVGEAR or KPP simulations, the call to SETUP_CO2SAT is made from 
!  routine "chemdr.f", after the "chem.dat" file is read.  This is necessary 
!  since we have to reference the SMVGEAR rxn rate and species numbers.
!                                                                             .
!  For offline simulations, the call to SETUP_CO2SAT can be made
!  at the start of the GEOS-Chem run (in "ndxx_setup.f" or similar routine).
! 
! !REVISION HISTORY: 
!  30 Jul 2002 - M. Evans    - Initial version
!  (1 ) Rename from "plane.dat" to "plane.log", since "*.dat" implies an input
!        file name. (bmy, 8/8/03)
!  (2 ) Add fancy output string (bmy, 4/26/04)
!  (3 ) Now references GET_NYMD, GET_NHMS, and EXPAND_DATE from "time_mod.f".
!        Now also replaces date & time tokens in the filenames. (bmy, 7/20/04)
!  (4 ) Now references FILE_EXISTS from "file_mod.f".  Modified so that we
!        check if a flight track file exists on each day.  Open file for 
!        output on each day and write header. (bmy, 3/25/05)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: IS_OPEN, surfheight_TF
      INTEGER            :: I,  IP,      N,   TEMP, LENGTH
      INTEGER            :: RN, COUNTER, IOS, NYMD, NHMS
      INTEGER            :: NCID, DIM_LEN, OBSCNT
      INTEGER            :: st2d(2),ct2d(2),st1d(1),ct1d(1)
      INTEGER, dimension(:), allocatable :: ID, ID4
      !INTEGER*8, dimension(:), allocatable :: ID
      REAL*8, dimension(:), allocatable :: LAT, LON
      REAL*8, dimension(:), allocatable :: XCO2_OBS,XCO2_UNCERT
      REAL*8, dimension(:), allocatable :: XCO2PBL_OBS,XCO2PBL_UNCERT
      REAL*8, dimension(:), allocatable :: XCO2FT_OBS,XCO2FT_UNCERT
      REAL*8, dimension(:), allocatable :: XCO2_LATITUDE,XCO2_LONGITUDE
      REAL*8, dimension(:), allocatable :: XCO2_CLOUDHEIGHT
      CHARACTER(LEN=7)   :: NAMES
      CHARACTER(LEN=20)  :: LINE
      CHARACTER(LEN=10)  :: TYPE
      
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object

      !=================================================================
      ! SETUP_CO2SAT begins here!
      !=================================================================

      ! Assume that there is flight data for today
      DO_CO2SAT = .TRUE.

      ! Find a free file LUN
      IU_FILE = findFreeLun()

      ! Get date & time
      NYMD  = GET_NYMD()
      NHMS  = GET_NHMS()

      ! Copy file names to local variables and
      !  replace and date & time tokens in the file names
      INF   = INFILENAME
      CALL EXPAND_DATE( INF,  NYMD, NHMS )
      
      ! Output NetCDF file monthly instead of daily
      IF ( ITS_A_NEW_MONTH() ) THEN
          ! Close previous month file
          ! I am worried about how the compiler might initial
          IF (SCO2_NCID /= 0) THEN
            CALL Nccl( SCO2_NCID )
          ENDIF
          
          ! Get number of observations in month
          OBSNUM = 1
          CALL COUNT_OBS_IN_MONTH( NYMD, NHMS, OBSCNT)
          
          OUTF  = OUTFILENAME 
          CALL EXPAND_DATE( OUTF, NYMD, NHMS )
          CALL CO2SAT_OUTPUT_CREATE( OUTF, OBSCNT)
      ENDIF

      PRINT *,'INF:',INF
      ! If we can't find a flighttrack file for today's date, return
      IF ( .not. FILE_EXISTS( INF ) ) THEN 
         DO_CO2SAT = .FALSE.
         RETURN
      ENDIF
      
      ! Echo info
      ! Echo info
      IF ( am_I_Root ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a)' ) 'CO2 S A T E L L I T E   D I A G N O S T I C'
       WRITE( 6, 100   ) TRIM( INF )
  100  FORMAT( /, 'SETUP_CO2SAT: Reading ',a )
       WRITE( 6, '(a)' )  
      ENDIF

      ! Compute # of species and # of points & allocate arrays
      CALL INIT_CO2SAT( am_I_Root )

      ! Return if there are no flight track points for today
      IF ( NPOINTS == 0 ) THEN
         IF ( am_I_Root ) THEN
           WRITE( 6, '(a)' ) 'No flight track found for today!'
         ENDIF
         DO_CO2SAT = .FALSE.
         RETURN
      ENDIF

      ! Read variables to be output -- sort into PVAR array by type
      ! Leave hook for this in future maybe ... aes
      !CALL READ_VARIABLES
      PVAR(1) = 100001  
      PNAME(1) = 'TRA_001   '

      ! Read information about each point (date/time/lon/lat/alt)
      CALL READ_POINTS( am_I_Root )

      ! Set the pointer to the first record 
      ! Set the pointer to the first record 
      PPOINT = 1

      !=================================================================
      ! Find the species # for all components of RO2 (SMVGEAR only)
      !=================================================================
      CALL RO2_SETUP( am_I_Root, Input_Opt, RC )

      !=================================================================
      ! Find the species # for all components of NOY (SMVGEAR only)
      !=================================================================
      CALL NOY_SETUP( am_I_Root, Input_Opt, RC )
      
      !=================================================================
      ! Find the species # for all components of AN (SMVGEAR only)
      !=================================================================
      CALL AN_SETUP( am_I_Root, Input_Opt, RC )
    
      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Allocate and initialize arrays
      ALLOCATE( ID4( NPOINTS ) )
      ID4 = 0
      ALLOCATE( ID( NPOINTS ) )
      ID = 0
      ALLOCATE( LAT( NPOINTS ) )
      LAT= 0.0_8
      ALLOCATE( LON( NPOINTS ) )
      LON= 0.0_8
      ALLOCATE( XCO2_CLOUDHEIGHT (NPOINTS) )
      XCO2_CLOUDHEIGHT= 0.0_8
      ALLOCATE( XCO2_OBS (NPOINTS) )
      XCO2_OBS= 0.0_8
      ALLOCATE( XCO2_UNCERT (NPOINTS) )
      XCO2_UNCERT= 0.0_8
      ALLOCATE( XCO2FT_OBS (NPOINTS) )
      XCO2FT_OBS= 0.0_8
      ALLOCATE( XCO2FT_UNCERT (NPOINTS) )
      XCO2FT_UNCERT= 0.0_8
      ALLOCATE( XCO2PBL_OBS (NPOINTS) )
      XCO2PBL_OBS= 0.0_8
      ALLOCATE( XCO2PBL_UNCERT (NPOINTS) )
      XCO2PBL_UNCERT= 0.0_8
      
       
      ! CALL READGOSAT_ID( trim (INF), NPOINTS, ID)
      !PRINT *,'ID1:',ID

      !  CALL READGOSAT_REAL4( trim (INF), 'xco2_raw', 'Retrieval',
      !&  ACOS_XCO2, NPOINTS)

      CALL Ncop_Rd( NCID, trim(INF) )
      
      !CALL Ncget_Unlim_Dimlen (NCID, DIM_LEN)
      CALL NCget_DimLen(NCID, 'sounding_id', DIM_LEN)

      !PRINT *,'DIM:',DIM_LEN
      !PRINT *,'npoints:',npoints

      st1d = (/ 1 /)
      ct1d = (/ DIM_LEN /)
      
      IF(Ncdoes_Var_Exist (ncid, 'latitude')) THEN
        CALL Ncrd(LAT, ncid, 'latitude',
     &  st1d,ct1d)
      ENDIF
      
      IF(Ncdoes_Var_Exist (ncid, 'longitude')) THEN
        CALL Ncrd(LON, ncid, 'longitude',
     &  st1d,ct1d)
      ENDIF

      IF (VERTINTOPT .eq. 3) THEN
        IF(Ncdoes_Var_Exist (NCID, 'xco2_observed')) THEN
            CALL Ncrd(XCO2_OBS, NCID, 'xco2_observed',
     &          st1d,ct1d)
        ENDIF
      ENDIF
      
      IF (VERTINTOPT .eq. 1) THEN
        IF(Ncdoes_Var_Exist (NCID, 'xco2')) THEN
            CALL Ncrd(XCO2_OBS, NCID, 'xco2',
     &          st1d,ct1d)
        ENDIF
      ENDIF

      IF(Ncdoes_Var_Exist (NCID, 'xco2_uncertainty')) THEN
        CALL Ncrd(XCO2_UNCERT, NCID, 'xco2_uncertainty',
     &          st1d,ct1d)
      ENDIF
      
      IF(Ncdoes_Var_Exist (NCID, 'surface_height')) THEN
        CALL Ncrd(XCO2_CLOUDHEIGHT, NCID, 'surface_height',
     &          st1d,ct1d)
      ELSE
           XCO2_CLOUDHEIGHT = 0.0_8
      ENDIF

      IF(Ncdoes_Var_Exist (NCID, 'xco2freetrop_observed')) THEN
        CALL Ncrd(XCO2FT_OBS, NCID, 'xco2freetrop_observed',
     &          st1d,ct1d)
      ENDIF

      IF(Ncdoes_Var_Exist (NCID, 'xco2freetrop_uncertainty')) THEN
        CALL Ncrd(XCO2FT_UNCERT, NCID, 'xco2freetrop_uncertainty',
     &          st1d,ct1d)
      ENDIF

      IF(Ncdoes_Var_Exist (NCID, 'xco2pbl_observed')) THEN
        CALL Ncrd(XCO2PBL_OBS, NCID, 'xco2pbl_observed', st1d, 
     &          ct1d)
      ENDIF

      IF(Ncdoes_Var_Exist (NCID, 'xco2pbl_uncertainty')) THEN
        CALL Ncrd(XCO2PBL_UNCERT, NCID, 'xco2pbl_uncertainty', 
     &          st1d, ct1d)
      ENDIF
      
      CALL CO2SAT_OBS_OUTPUT_WRITE(ID, LAT, LON,
     &      XCO2_OBS, XCO2_UNCERT, XCO2_CLOUDHEIGHT,
     &      XCO2FT_OBS, XCO2FT_UNCERT, XCO2PBL_OBS,
     &      XCO2PBL_UNCERT)

      CALL Nccl( NCID )
      
      END SUBROUTINE SETUP_CO2SAT
      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_CO2SAT
!
! !DESCRIPTION: Subroutine INIT\_CO2SAT reads the input file to compute 
!  the number of variables and flight track points to print out.  Also 
!  allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CO2SAT( am_I_Root )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IOERROR
      USE m_netcdf_io_open,  ONLY : Ncop_Rd
      USE m_netcdf_io_get_dimlen,   ONLY : Ncget_DimLen
      USE m_netcdf_io_close,  ONLY : Nccl
!      USE NCCOOARDSFORMAT,  ONLY : readGOSATlength

      USE CMN_SIZE_MOD    ! Size Parameters
      USE COMODE_LOOP_MOD ! ITLOOP
      
      include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run.  Also reference ALLOC_ERR
!        from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  (3 ) INIT_CO2SAT is now called each day but the arrays are only
!        allocated once.  Arrays are now allocated to the maximum size.
!        (bmy, 3/25/05)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL           :: IS_INIT = .FALSE.
      INTEGER           :: N, AS, IOS, INIT_NCID
      CHARACTER(LEN=20) :: LINE

      !=================================================================
      ! INIT_CO2SAT begins here!
      !=================================================================

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      !=================================================================
      ! For now, hardcoding at one variable, TRA_004, CO2
      !=================================================================

      NPVAR = 1

      !=================================================================
      ! Read plane track points -- plane, lat/lon/alt, date/time
      !=================================================================

      NPOINTS = 0
      
      ! RM 42014 - Changed to use GeosChem NetCDF
      CALL Ncop_Rd( INIT_NCID, trim(INF) )

      CALL NCget_DimLen(INIT_NCID, 'sounding_id', NPOINTS)

      ! If there are no flight-track points then just return
      IF ( NPOINTS < 1 ) THEN
         DO_CO2SAT = .FALSE.
         RETURN
      ENDIF

      ! Make sure NPOINTS is less than MAXPOINTS
      IF ( NPOINTS > MAXPOINTS ) THEN
         WRITE( 6, '(a)') 'NPOINTS exceeds maximum allowed value!'
         WRITE( 6, '(a)') 'STOP in INIT_CO2SAT (CO2SAT_mod.f)'
         WRITE( 6, '(a)') REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
         
      !=================================================================
      ! Allocate arrays to maximum sizes
      !
      ! NOTE: To save space, NPREAC is the actual number of reactions
      !       found.  We will worry about this later.  (bmy, 3/25/05)
      !
      ! Some of this is not used but kept in as "hooks" in case we add 
      ! more species/retrievals, aes 
      !=================================================================
      IF ( .not. IS_INIT ) THEN 

         !-------------------------
         ! Arrays of size NPREAC
         !-------------------------
         ALLOCATE( PREAC( MAX( NPREAC, 1 ) ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PREAC' )

         ALLOCATE( PRRATE( ITLOOP, MAX( NPREAC, 1 ) ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRRATE' )

         !--------------------------
         ! Arrays of size MAXVARS
         !--------------------------
         ALLOCATE( PVAR( MAXVARS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PVAR' )

         ALLOCATE( PNAME( MAXVARS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PNAMES' )

         !---------------------------
         ! Arrays of size MAXPOINTS
         !---------------------------
         ALLOCATE( PTYPE( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTYPE' )

         ALLOCATE( PDATE( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDATE' )

         ALLOCATE( PTIME( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTIME' )

         ALLOCATE( PTAU( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTAU' )
         
         ALLOCATE( PLAT( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLAT' )

         ALLOCATE( PLON( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLON' )

         ALLOCATE( PPRESS( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PPRESS' )

         ALLOCATE( PPRESSWGHT( 47, MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PPRESSWGHT' )

         ALLOCATE( PAK( 20, MAXPOINTS), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PAK' )

         ALLOCATE( POCO2_PR( 20, MAXPOINTS), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'POCO2_PR' )

         ALLOCATE( POCO2_PRESS( 20, MAXPOINTS), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'POCO2_PRESS' )

         ALLOCATE( PERR( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PERR' )
 
         ALLOCATE( PXCO2( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PXCO2' )
 
         ALLOCATE( PID( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PID' )
        
         ! Reset IS_INIT flag
         IS_INIT = .TRUE.
      ENDIF

      !=================================================================
      ! Initialize arrays 
      !=================================================================
      PREAC  = 0
      PRRATE = 0e0
      PVAR   = 0
      PNAME  = ''
      PTYPE  = ''
      PDATE  = 0e0
      PTIME  = 0e0
      PTAU   = 0e0
      PLAT   = 0e0
      PLON   = 0e0
      PPRESS = 0e0
      PPRESSWGHT  = 0e0
      PAK         = 0e0
      POCO2_PR    = 0e0
      POCO2_PRESS = 0e0
      PID    = '              '
      
      CALL Nccl(INIT_NCID)

      END SUBROUTINE INIT_CO2SAT
      
      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_points
!
! !DESCRIPTION: Subroutine READ\_POINTS reads the information (ID, date, time, 
!  lat, lon, pressure) for each measurement listed in the input file, and 
!  sorts these into the appropriate module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_POINTS( am_I_Root )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_TAU0
      USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
      USE JULDAY_MOD,        ONLY : CALDATE
      USE TIME_MOD,      ONLY : YMD_EXTRACT
      USE m_netcdf_io_read
      USE m_netcdf_io_close
      USE m_netcdf_io_open
      USE m_netcdf_io_checks
      
      include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  30 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                              :: N, IOS, YYYYMMDD, HHMMSS
      INTEGER                              :: QYY,QMM,QDD, QHH, QMN, SS
      INTEGER                              :: NCID
      INTEGER                              :: st2d(2),ct2d(2)
      INTEGER                              :: st1d(1),ct1d(1)
      REAL*8, dimension(:), allocatable    :: TIM
      REAL*8                               :: JULNDAY
      REAL*4, dimension(:), allocatable    :: LAT, LON, XCO2, ERR, PRESS
      REAL*4, dimension(:,:), allocatable  ::  PRESSWGHT, AK, OCO2_PR
      REAL*4, dimension(:,:), allocatable  ::  OCO2_PRESS
      INTEGER*8, dimension(:), allocatable :: ID
      CHARACTER(LEN=10) :: TYPE
      CHARACTER(LEN=14) :: SOUNDING_ID_CHAR
      !=================================================================
      ! READ_POINTS begins here!
      !=================================================================
      
      !=================================================================
      ! Read plane track points -- plane, lat/lon/alt, date/time
      ! We have previously computed NPOINTS in INIT_CO2SAT
      !=================================================================

      ALLOCATE( LAT( NPOINTS) )
      ALLOCATE( LON( NPOINTS) )
      ALLOCATE( XCO2( NPOINTS ))
      ALLOCATE( ERR( NPOINTS ))
      ALLOCATE( ID( NPOINTS ))
      ALLOCATE( TIM( NPOINTS ))
      ALLOCATE( PRESS( NPOINTS ))
      ALLOCATE( PRESSWGHT( 47, NPOINTS) )
      ALLOCATE( AK( 20, NPOINTS) )
      ALLOCATE( OCO2_PR( 20, NPOINTS) )
      ALLOCATE( OCO2_PRESS( 20, NPOINTS) )

      print *,'NPOINTS:',NPOINTS
      !print *,'INF:',trim(INF)
      
      ! RM 42014 - Changed to use GeosChem NetCDF
      CALL Ncop_Rd( NCID, trim(INF) )
      
      st1d = (/ 1 /)
      ct1d = (/ NPOINTS /)
      
      IF(Ncdoes_Var_Exist (ncid, 'latitude')) THEN
        CALL Ncrd(LAT, ncid, 'latitude',
     &  st1d,ct1d)
      ENDIF
      
      IF(Ncdoes_Var_Exist (ncid, 'longitude')) THEN
        CALL Ncrd(LON, ncid, 'longitude',
     &  st1d,ct1d)
      ENDIF
      
      IF(Ncdoes_Var_Exist (ncid, 'time')) THEN
        CALL Ncrd(TIM, ncid, 'time',
     &  st1d,ct1d)
      ENDIF
      
      ! RM42014 - looks like there is an int*8 in m_netcdf_io_read.F90
      IF(Ncdoes_Var_Exist (ncid, 'sounding_id')) THEN
        CALL Ncrd(ID, ncid, 'sounding_id',
     &  st1d,ct1d)
      ENDIF

      IF( Ncdoes_Var_Exist (ncid, 'pressure_weight')) THEN
         IF ( VERTINTOPT .eq. 3 ) THEN
           st2d = (/ 1 , 1 /)
           ct2d = (/ 47, NPOINTS /)
           CALL Ncrd(PRESSWGHT,NCID,'pressure_weight',st2d,ct2d)
         ENDIF
         
         !IF ( VERTINTOPT .eq. 1 ) THEN
         !  st2d = (/ 1 , 1 /)
         !  ct2d = (/ 20, NPOINTS /)
         !  CALL Ncrd(,NCID,'pressure_weight',st2d,ct2d)
         !ENDIF 
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'xco2_averaging_kernel')) THEN
        st2d = (/ 1, 1/)
        ct2d = (/ 20, NPOINTS /)
        CALL Ncrd(AK, ncid, 'xco2_averaging_kernel',
     &  st2d,ct2d)
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'co2_profile_apriori')) THEN
        st2d = (/ 1, 1/)
        ct2d = (/ 20, NPOINTS /)
        CALL Ncrd(OCO2_PR, ncid, 'co2_profile_apriori',
     &  st2d,ct2d)
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'pressure_levels')) THEN
        st2d = (/ 1, 1/)
        ct2d = (/ 20, NPOINTS /)
        CALL Ncrd(OCO2_PRESS, ncid, 'pressure_levels',
     &  st2d,ct2d)
      ENDIF

      CALL Nccl( NCID )

      !PRINT *,'shape plat:',SHAPE(PLAT)
      !PRINT *,'shape lat:',SHAPE(LAT)

      !PRINT *,'shape ppresswgth:',SHAPE(PPRESSWGHT)
      !PRINT *,'shape presswgth:',SHAPE(PRESSWGHT)
      !PLAT = LAT
      !PLON = LON
      !PPRESSWGHT(:,1:NPOINTS) = PRESSWGHT
      !PPRESS = PRESS
      !PXCO2 = XCO2
      !PERR = ERR
      PRINT *,'done w/ first part'

       DO N = 1, NPOINTS

         !==============================================================
         ! Read date and time coordinates -- also do error checks
         !============================================================== 
         !-- READ SOUNDING_ID AND PARSE INTO INTEGER DATE COMPS
         !WRITE(SOUNDING_ID_CHAR, '(i14.14)') ID(N)
         !PID(N) = SOUNDING_ID_CHAR
         !READ( SOUNDING_ID_CHAR(1:4),'(i)') QYY
         !READ( SOUNDING_ID_CHAR(5:6),'(i)') QMM  
         !READ( SOUNDING_ID_CHAR(7:8), '(i)' ) QDD
         !READ( SOUNDING_ID_CHAR(9:10), '(i)' ) QHH
         !READ( SOUNDING_ID_CHAR(11:12), '(i)' ) QMN

         JULNDAY = TIM(N)/86400 + 2440587.5
         !PRINT *,'GOSAT/JULNDAY:',JULNDAY
         CALL CALDATE( JULNDAY, YYYYMMDD, HHMMSS )

         CALL YMD_EXTRACT( YYYYMMDD, QYY, QMM, QDD ) 
         CALL YMD_EXTRACT( HHMMSS, QHH, QMN, SS )


         !print *,'POINT number:',N
         !print *,'PID(N):',PID(N)
 
        ! Error check MONTH
         IF ( QMM < 1 .or. QMM > 12 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 105   ) QMM
 105           FORMAT( 'ERROR: MONTH out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (GOSAT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Error check DAY
         IF ( QDD < 1 .or. QDD > 31 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 110   ) QDD
 110           FORMAT( 'ERROR: DAY out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (GOSAT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP 
         ENDIF

         ! Error check HOUR
         IF ( QHH < 0 .or. QHH > 23 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 115   ) QHH
 115           FORMAT( 'ERROR: HOUR out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (GOSAT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP  
         ENDIF

         ! Error check MINUTES
         IF ( QMN < 0 .or. QMN > 59 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 120   ) QMN
 120           FORMAT( 'ERROR: MINUTES out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (GOSAT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Store type in the global PTYPE array
         PTYPE(N) = TYPE

         ! Store YYYYMMDD in the global PDATE array
         PDATE(N) = ( QYY * 10000 ) + ( QMM * 100 ) + QDD

         ! Store HHMMSS in the global PTIME array
         ! (actaully we read in just HHMM, assume seconds = 00)
         PTIME(N) = ( QHH * 100 ) + QMN

         ! Store TAU (hours since 1 Jan 1985) in the global PTAU array
         PTAU(N)  = GET_TAU0( QMM, QDD, QYY, QHH, QMN, 0 )

         PLAT(N) = LAT(N)
         PLON(N) = LON(N)
         PERR(N) = ERR(N)
         PPRESSWGHT(:,N) = PRESSWGHT(:,N)
         PAK(:,N) = AK(:,N)
         POCO2_PR(:,N) = OCO2_PR(:,N)
         POCO2_PRESS(:,N) = OCO2_PRESS(:,N)

         !==============================================================
         ! Read lon/lat/alt coordinates -- also do error checks
         !==============================================================

         ! Put LONGITUDE in the range [-180...180]
         IF ( PLON(N) > 180.0 ) PLON(N) = PLON(N) - 360e0

         ! Error check LONGITUDE
         IF ( PLON(N) < -180 .OR. PLON(N) > 180 ) THEN 
            IF ( am_I_Root ) THEN
               WRITE( 6, 125   ) PLON(N)
 125           FORMAT( 'ERROR: Longitude out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (GOSAT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Error check LATITUDE
         IF ( PLAT(N) < -90.0 .OR. PLAT(N) > 90.0 ) THEN 
            IF ( am_I_Root ) THEN
               WRITE( 6, 130   ) PLAT(N)
 130           FORMAT( 'ERROR: Latitude out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (GOSAT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF
        
       ENDDO

      DEALLOCATE( LAT )
      DEALLOCATE( LON )
      DEALLOCATE( XCO2 )
      DEALLOCATE( ERR )
      DEALLOCATE( ID )
      DEALLOCATE( PRESS )
      DEALLOCATE( PRESSWGHT )
      DEALLOCATE( AK )
      DEALLOCATE( OCO2_PR)

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      IF ( am_I_Root )  WRITE( 6, 135 ) NPOINTS
 135  FORMAT( /, 'Number of flight track points : ', i6 )

      END SUBROUTINE READ_POINTS
      
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ro2_setup
!
! !DESCRIPTION: Subroutine RO2\_SETUP saves the SMVGEAR species indices of RO2
!  constituents in the PRO2 array.  Also computes the count NPRO2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RO2_SETUP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Aug 2003 - M. Evans    - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (2 ) Now replace NAMESPEC w/ NAMEGAS for SMVGEAR II (bmy, 8/1/03)
!  (3 ) Now references ITS_A_FULLCHEM_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Bug fix: PO3 should be PO2 (tmf, bmy, 4/23/07)
!  (5 ) NOTE: PO3 was a bug, that should have been PO2 (tmf, 2/10/09)
!  (6 ) Add new RO2 species according to 'globchem.dat' (tmf, 3/10/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  25 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! RO2_SETUP begins here!
      !=================================================================

      ! Initialize
      NPRO2 = 0

      ! We only need to proceed for SMVGEAR chemistry
      IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

      !=================================================================
      ! Loop over all SMVGEAR species, test for RO2 components
      !=================================================================
      DO M = 1, NSPEC(NCS)

         ! If we have found an RO2 compoent, add its species # to
         ! the PRO2 global array, and increment counter
         ! NOTE: PO3 was a bug, that should have been PO2 (tmf, 2/10/09)
         SELECT CASE( TRIM( NAMEGAS(M) ) )

            CASE ( 'HO2',  'MO2',  'A3O2', 'ATO2', 'B3O2',
     &             'ETO2', 'GCO3', 'IAO2', 'KO2',  'MAO3',
     &             'MCO3', 'MRO2', 'PO2',  'RIO2', 'VRO2',
     &             'ACO3', 'EO2', 'ENCO3', 'ENO2', 'GLCO3',
     &             'IACO3', 'INO2', 'MACO3', 'NICO3', 'NIO2',
     &             'VOHRO2', 'RIO1', 'C59O2', 'MAPO2',
     &             'BRO2', 'TRO2', 'XRO2', 'NRO2',
     &             'MOBAOO', 'ISOPNBO2', 'ISOPNDO2', 'IEPOXOO',
     &             'PYPO2', 'MACRNO2', 'PMNO2', 'HC5OO', 'RCO3')
               NPRO2       = NPRO2 + 1
               PRO2(NPRO2) = M

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Error check
      IF ( NPRO2 > MAXRO2 ) THEN
         WRITE( 6, '(a)' ) 'NPRO2 exceeds maximum allowed value!'
         WRITE( 6, '(a)' ) 'STOP in RO2_SETUP (planeflight_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      WRITE( 6, 100 ) NPRO2
 100  FORMAT( 'Number of RO2 components      : ', i6 )

      END SUBROUTINE RO2_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: noy_setup
!
! !DESCRIPTION: Subroutine NOY\_SETUP saves the SMVGEAR species indices of NOy
!  constituents in the PNOY array.  Also computes the count NPNOY.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NOY_SETUP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Jun 2009 - F. Paulot   - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  25 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! NOY_SETUP begins here!
      !=================================================================

      ! Initialize
      NPNOY = 0

      ! We only need to proceed for SMVGEAR chemistry
      IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

      !=================================================================
      ! Loop over all SMVGEAR species, test for NOY components
      !=================================================================
      DO M = 1, NSPEC(NCS)

         SELECT CASE( TRIM( NAMEGAS(M) ) )

           CASE ( 'NO', 'NO2', 'NO3', 'HNO2', 'HNO4', 'HNO3',
     &            'PAN', 'PYPAN', 'PMN', 'PPN')

               NPNOY       = NPNOY + 1
               PNOY(NPNOY) = M

           CASE ( 'N2O5')

               NPNOY       = NPNOY + 1
               PNOY(NPNOY) = M

               NPNOY       = NPNOY + 1
               PNOY(NPNOY) = M

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Error check
      IF ( NPNOY > MAXNOY ) THEN
         WRITE( 6, '(a)' ) 'NPNOY exceeds maximum allowed value!'
         WRITE( 6, '(a)' ) 'STOP in NOY_SETUP (planeflight_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      WRITE( 6, 100 ) NPNOY
 100  FORMAT( 'Number of NOY components      : ', i6 )

      END SUBROUTINE NOY_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: an_setup
!
! !DESCRIPTION: Subroutine AN\_SETUP saves the SMVGEAR species indices of AN
!  constituents in the P\_AN array.  Also computes the count NPAN.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AN_SETUP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Jun 2009 - F. Paulot   - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  25 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! AN_SETUP begins here!
      !=================================================================

      ! Initialize
      NPAN = 0

      ! We only need to proceed for SMVGEAR chemistry
      IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

      !=================================================================
      ! Loop over all SMVGEAR species, test for AN components
      !=================================================================
      DO M = 1, NSPEC(NCS)

         ! If we have found an AN component, add its species # to
         ! the AN global array, and increment counter
         SELECT CASE( TRIM( NAMEGAS(M) ) )

           CASE ( 'ISOPNB', 'ISOPND', 'MVKN',  'MACRN',
     &            'ISN2',   'MNO3',   'ISNO3', 'R4N2',
     &            'PROPNN', 'ETHLN' )

               NPAN       = NPAN + 1
               P_AN(NPAN) = M

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Error check
      IF ( NPAN > MAXAN ) THEN
         WRITE( 6, '(a)' ) 'NPAN exceeds maximum allowed value!'
         WRITE( 6, '(a)' ) 'STOP in AN_SETUP (planeflight_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      WRITE( 6, 100 ) NPAN
 100  FORMAT( 'Number of AN components      : ', i6 )

      END SUBROUTINE AN_SETUP
      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GOSAT
!
! !DESCRIPTION: Subroutine CO2SAT saves concentrations to disk at 
!  locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CO2SAT(am_I_Root,Input_Opt,State_Met,State_Chm,RC)
!
! !USES:
!
      USE COMODE_MOD,   ONLY : AIRDENS,        CSPEC,       JLOP
      USE COMODE_MOD,   ONLY : T3,             VOLUME,      ABSHUM
      USE COMODE_MOD,   ONLY : TAREA
#if   defined( TOMAS )
      USE DIAG_MOD,           ONLY : AD61_INST   ! (win, 7/28/09)
#endif
!     ADDING FOR PBL CO2 vals     
!      USE DIAG_MOD,           ONLY : AD12
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP  
      USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TAU,        GET_TS_DIAG
      USE OCEAN_MERCURY_MOD,  ONLY : Fp, Fg !eds 10/27/11
      USE PBL_MIX_MOD,   ONLY : GET_FRAC_UNDER_PBLTOP

      USE JV_CMN_MOD
      USE COMODE_LOOP_MOD    ! CSPEC, etc.
      
      USE ERROR_MOD,  ONLY : ERROR_STOP
      
      USE LINEAR_INTERP_MOD

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now reference AD from "dao_mod.f".  Now references GEOS_CHEM_STOP from
!        "error_mod.f", which frees memory before stopping. (bmy, 10/15/02)
!  (2 ) Now uses functions GET_TAU, GET_TS_CHEM from "time_mod.f".
!        (bmy, 3/27/03)
!  (3 ) Updated comments, cosmetic changes (bmy, 7/18/03)
!  (4 ) Now references T from "dao_mod.f", so that we can save out temperature
!        for non-SMVGEAR runs. (bmy, 8/1/03)
!  (5 ) Now references UWND and VWND from "dao_mod.f".  Now references
!        GET_PEDGE from "pressure_mod.f".  Added CASEs for surface pressure,
!        UWND, VWND to the CASE statement (bmy, 4/23/04)
!  (6 ) Now references STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  (7 ) Now return if DO_CO2SAT = .FALSE. (bmy, 3/24/05)
!  (8 ) Now compute column AOD's and AOD's below plane.  Also now scale
!        AOD's to 400nm. (bmy, 10/24/05)
!  (9 ) Bug fix: exit if PTAU(M) == PTAUE, so that we write out on the next !
!        CO2SAT timestep (cdh, bmy, 12/12/06)
!  (10) Change CO2SAT output time step. (ccc, 8/27/09)
!  (11) Now scale AOD's to jv_spec_aod.dat wavelength. (clh, 5/14/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: PCHEM, LINTERP
      INTEGER             :: I, IP, IRHN, J, L, JLOOP, M, N, R, RH, V
      INTEGER             :: LL,LPLANE, IWV, ISPC, K
      INTEGER             :: VERTIND
      REAL*8              :: TK, PTAUS, PTAUE, CONSEXP, VPRESH2O, SAODnm
      REAL*8              :: VARI(NPVAR)
      REAL*8              :: SUMPRESS
      REAL*8              :: VARIPBL(NPVAR)
      REAL*8              :: VARIFREETROP(NPVAR)
      REAL*8, ALLOCATABLE :: PARTVARI(:)

      REAL*8                :: AK_INTERP(47)
      REAL*8                :: OCO2_PR_INTERP(47)
      REAL*8                :: GC_PRESS_IN(47)
      REAL*8                :: GC_PRESS_WGHT(47)

      REAL*8                :: TMPSUM(47)

      ! ASCENDS KERNEL (D. Baker)
      REAL*8              :: KERN(47) = (/ 0.000035000, 0.000230000,
     &                        0.000881000, 0.002041000, 0.002992000,
     &                        0.005694000, 0.004604000, 0.006408000,
     &                        0.008926000, 0.012436000, 0.017318000,
     &                        0.024105000, 0.032491000, 0.028132000,
     &                        0.031085000, 0.033998000, 0.036856000,
     &                        0.039639000, 0.042359000, 0.045018000,
     &                        0.047600000, 0.050091000, 0.034747000,
     &                        0.035797000, 0.036820000, 0.037821000,
     &                        0.038803000, 0.031736000, 0.024200000,
     &                        0.024528000, 0.024855000, 0.025181000,
     &                        0.025500000, 0.025813000, 0.013062000,
     &                        0.013062000, 0.026433000, 0.017822667,
     &                        0.015668833, 0.013515000, 0.023152500,
     &                        0.013440364, 0.005020727, 0.005020727,
     &                        0.005020727, 0.005020727, 0.005020727 /)

      LOGICAL, PARAMETER       :: OCO2 = .FALSE.

      ! Aerosol types: SULF, BLKC, ORGC, SALA, SALC 
      INTEGER             :: IND(5) = (/ 22, 29, 36, 43, 50 /)
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER  :: MISSING = -999.99999999d0  ! Missing data value
      REAL*8,  PARAMETER  :: TINY    = 1.d-36           ! arbitary small # to 
                                                        !  avoid faulty output
      ! Pointers
      ! We need to define local arrays to hold corresponding values
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CO2SAT begins here!
      !=================================================================

      !print *,'AD12:',AD12(40,65,:)

      ! Return if there is no flighttrack data for today
      PRINT *,'DO_CO2SAT:',DO_CO2SAT
      IF ( .not. DO_CO2SAT ) RETURN

      !Determine if optical properties need interpolating
      !The LUT wavelengths in IWVSELECT will match if no interpolation
      !is needed. (output is only for the first requested wavelength)
      !(DAR 10/2013)
      !IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
      !   LINTERP=.FALSE.
      !ELSE
      !   LINTERP=.TRUE.
      !ENDIF

      ! Loop over all the locations that have not yet been found
      DO M = PPOINT, NPOINTS
         !print *,'ppoint:',ppoint
         !print *,'npoints:',npoints
         
         
         ! Starting & end times of chemistry interval
         PTAUE = GET_TAU()
         PTAUS = PTAUE - ( GET_TS_DIAG() / 60d0 )

         ! Initialize VARI to missing value for this point
         DO V = 1, NPVAR
            VARI(V) = MISSING
         ENDDO

         !print *,'ptaue:',PTAUE
         !print *,'ptau(m):',PTAU(M)
         !print *,'ptaus:',PTAUS

         !==============================================================
         ! We haven't found the first plane point yet...
         !==============================================================
         IF ( PTAU(M) < PTAUS ) THEN

            ! Write all missing values to disk for point #M
            !CALL WRITE_VARS_TO_FILE( M, VARI )
            ! This should never be called as it isn't opened -RM
!            CALL ERROR_STOP( 'PTAU(M) < PTAUS write_vars_to_file',
!     &         'Not sure if this is really an error (co2sat_mod.f)' )
            

            ! Increment pointer
            PPOINT = PPOINT + 1

         !==============================================================
         ! We have already found all of the plane points...
         !==============================================================
         ELSE IF ( PTAU(M) >= PTAUE ) THEN

            ! Exit this loop and the subroutine
            EXIT

         !==============================================================
         ! We have found a plane point at the proper time & location!
         !==============================================================
         ELSE
            
            ! Print the flight track point number
            !-- Commenting out for now, AESchuh
            !WRITE( 6, 100 ) PTYPE(M), PDATE(M), PTIME(M) 
 100        FORMAT( '     - CO2SAT: Archived ',a5,1x,i8.8,1x,i4.4 )

            ! Return grid box indices for the chemistry region
            ! NOTE: PCHEM and JLOOP are only defined for SMVGEAR runs!
            CALL TEST_VALID( M, PCHEM, JLOOP, I, J, L, Input_Opt )

            ! Initialize SMVGEAR reaction counter
            R = 0

            ! Loop over all variables to save out
            DO V = 1, NPVAR

               ! Handle each variable
               SELECT CASE ( PVAR(V) )

                  !-------------------------
                  ! SMVGEAR species
                  !-------------------------
                  !997 and 998 are now special flags for NOy and AN
                  CASE ( 1:996) !( 1:998 )

                     ! Only archive where SMVGEAR chem is done
                     ! Save as mixing ratio [v/v]
                     IF ( JLOOP /= 0 ) THEN
                        VARI(V) = CSPEC(JLOOP,PVAR(V)) / AIRDENS(JLOOP)
                     ENDIF

                  ! FP 04/01/2010
                  !-------------------------
                  ! NOy family
                  !-------------------------   
                  CASE ( 997 )

                     ! Only archive where SMVGEAR chem is done
                     ! Sum all AN contributions, save as [v/v]
                     VARI(V) = 0d0

                     IF ( JLOOP /= 0 ) THEN 
                        DO N = 1, NPNOY
                           VARI(V) = VARI(V) + CSPEC(JLOOP,PNOY(N))
                        ENDDO

                        VARI(V) = VARI(V) / AIRDENS(JLOOP)
                     ENDIF

                  ! FP 04/01/2010
                  !-------------------------
                  ! AN family
                  !-------------------------   
                  CASE ( 998 )

                     ! Only archive where SMVGEAR chem is done
                     ! Sum all AN contributions, save as [v/v]
                     VARI(V) = 0d0
                     IF ( JLOOP /= 0 ) THEN 

                        DO N = 1, NPAN
                           VARI(V) = VARI(V) + CSPEC(JLOOP,P_AN(N))
                        ENDDO

                        VARI(V) = VARI(V) / AIRDENS(JLOOP)
                     ENDIF

                  !-------------------------
                  ! RO2 family
                  !-------------------------   
                  CASE ( 999 )

                     ! Only archive where SMVGEAR chem is done
                     ! Sum all RO2 contributions, save as [v/v]
                     VARI(V) = 0d0
                           
                        IF ( JLOOP /= 0 ) THEN
                           DO N = 1, NPRO2
                              VARI(V) = VARI(V) + CSPEC(JLOOP,PRO2(N))
                           ENDDO
                           
                           VARI(V) = VARI(V) / AIRDENS(JLOOP)
                        ENDIF

                  !--------------------------
                  ! GMAO temperature [K]
                  !--------------------------
                  CASE ( 1001 )
                     VARI(V) = State_Met%T(I,J,L)

                  !--------------------------
                  ! GMAO abs humidity [frac]
                  !--------------------------
                  CASE ( 1002 ) 
                     
                     ! Only archive where SMVGEAR chem is done
                     ! Code skalooched from "calcrate.f"
                        IF ( JLOOP /= 0 ) THEN
                           TK       = T3(JLOOP)
                           CONSEXP  = 17.2693882d0 * 
     &                                (TK - 273.16d0) / (TK - 35.86d0)
                           
                           VPRESH2O = CONSVAP * EXP(CONSEXP) * 1d0 / TK
                           
                           VARI(V)  = ABSHUM(JLOOP) * 
     &                                VPRESH2O      / AIRDENS(JLOOP)
                        ENDIF

                  !--------------------------
                  ! GMAO aerosol sfc area
                  !--------------------------
                  CASE ( 1003 )

                     ! Only archive where SMVGEAR chem is done
                        VARI(V) = 0d0

                        IF ( JLOOP /= 0 ) THEN
                           DO N = 1, NDUST + NAER
                              VARI(V) = VARI(V) + TAREA(JLOOP,N)
                           ENDDO
                        ENDIF

                  !--------------------------
                  ! GMAO sfc pressure [hPa]
                  !--------------------------
                  CASE ( 1004 )
                     VARI(V) = GET_PEDGE(I,J,1)

                  !-------------------------
                  ! GMAO U-wind [m/s]
                  !-------------------------
                  CASE ( 1005 )
                     VARI(V) = State_Met%U(I,J,L)

                  !--------------------------
                  ! GMAO V-wind [m/s]
                  !--------------------------
                  CASE ( 1006 )
                     VARI(V) = State_Met%V(I,J,L)

                  !--------------------------
                  ! GEOS-Chem Grid Box I
                  !--------------------------
                  CASE ( 1007 )
                     VARI(V) = I

                  !--------------------------
                  ! GEOS-Chem Grid Box J
                  !--------------------------
                  CASE ( 1008 )
                     VARI(V) = J

                  !--------------------------
                  ! GEOS-Chem Grid Box L
                  !--------------------------
                  CASE ( 1009 )
                     VARI(V) = L

#if   defined( MERRA ) || defined( GEOS_FP )
                  !--------------------------
                  ! GEOS-Chem SEAICE frac's
                  !--------------------------
                  CASE ( 1100 )
                     VARI(V) = State_Met%SEAICE00(I,J)
                  CASE ( 1101 )
                     VARI(V) = State_Met%SEAICE10(I,J)
                  CASE ( 1102 )
                     VARI(V) = State_Met%SEAICE20(I,J)
                  CASE ( 1103 )
                     VARI(V) = State_Met%SEAICE30(I,J)
                  CASE ( 1104 )
                     VARI(V) = State_Met%SEAICE40(I,J)
                  CASE ( 1105 )
                     VARI(V) = State_Met%SEAICE50(I,J)
                  CASE ( 1106 )
                     VARI(V) = State_Met%SEAICE60(I,J)
                  CASE ( 1107 )
                     VARI(V) = State_Met%SEAICE70(I,J)
                  CASE ( 1108 )
                     VARI(V) = State_Met%SEAICE80(I,J)
                  CASE ( 1109 )
                     VARI(V) = State_Met%SEAICE90(I,J)
#endif

                  !--------------------------
                  ! Column aerosol optical
                  ! depths [unitless]
                  !--------------------------
                  CASE ( 2001:2005 )

                     ! Only archive where SMVGEAR chem is done

                        ! Remove MISSING flag
                        VARI(V) = 0d0

                        ! Aerosol number
                        N = PVAR(V) - 2000

                        ! Loop over RH bins
                        DO RH = 1, NRH

                           ! Scaling factor for wavelength specified in
                           ! jv_spec_aod.dat
    !                       SAODnm  = QAA_AOD(IND(N)+RH-1) /
    ! &                               QAA(4,IND(N)+RH-1)
                           ! We no longer need to scale by wavelength
                           ! (skim, 02/03/11)
                           SAODnm = 1.0

                           ! Index for type of aerosol and RH value
                           IRHN    = ( (N-1) * NRH ) + RH

                           ! Sum AOD over all RH bins and store in VARI(V)
                           ! Sum over all vertical levels (bmy, 10/24/05)
                           VARI(V) = VARI(V) +
     &                               SUM( SAODnm * ODAER(I,J,:,IRHN) )
                        ENDDO



                  !--------------------------
                  ! Aerosol optical depths
                  ! below plane [unitless]
                  !--------------------------
                  CASE ( 3001:3005 )

                     IF ( PCHEM ) THEN

                       ! Remove MISSING flag
                        VARI(V) = 0d0

                        ! Aerosol number
                        N = PVAR(V) - 3000

                        ! Loop over RH bins
                        DO RH = 1, NRH

                           ! Scaling factor for wavelength specified in
                           ! jv_spec_aod.dat
    !                       SAODnm  = QAA_AOD(IND(N)+RH-1) /
    ! &                               QAA(4,IND(N)+RH-1)
                           ! We no longer need to scale by wavelength
                           ! (skim, 02/03/11)
                           SAODnm = 1.0

                           ! Index for type of aerosol and RH value
                           IRHN    = ( (N-1) * NRH ) + RH

                           ! Level of the plane.  AOD's are only computed
                           ! up to the tropopause, so if the plane goes into
                           ! the stratosphere, the AOD below plane will be
                           ! the same as the trop column at that point.
                           ! (bmy, 10/24/05)
                           LPLANE  = MIN( L, LLTROP )

                           ! Sum AOD over all RH bins and store in VARI(V)
                           ! Sum from surface to level where the plane is
                           VARI(V) = VARI(V) +
     &                          SUM( SAODnm * ODAER(I,J,1:LPLANE,IRHN) )
                        ENDDO
                     ENDIF

                  !--------------------------
                  ! SMVGEAR reaction rates
                  !--------------------------
                  CASE ( 10000:99999 )

                     ! Increment reaction count
                     R = R + 1

                     ! Only archive where SMVGEAR chem is done 
                     IF ( JLOOP /= 0 ) VARI(V) = PRRATE(JLOOP,R)

                  !--------------------------
                  ! GEOS-CHEM tracers [v/v]
                  !--------------------------
                  CASE( 100000:199999 )
                     !print *,'in CO2SAT tracer plane loop'
                     ALLOCATE(PARTVARI(LLPAR))

                     ! Remove offset from PVAR
                     N = PVAR(V) - 100000

                     ! Initialize GEOS-Chem tracer array [kg] from
                     ! Chemistry State object (mpayer, 12/6/12)
                     STT => State_Chm%Tracers

                     !print *,'N:',N
                     ! Convert from [kg] --> [v/v]
                     !VARI(V) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
                     !print *,'old:',STT(I,J,L,N)*TCVV(N) / AD(I,J,L)
                     !-- Temp kernel weighting for CO2SAT code




                !-------------------------------------------------------------------------
                !--  This option controls what type of vertical integration of CO2 is done
                !--  1:  OCO2
                !--  2:  ASCENDS FIXED KERNEL (KERN variable)
                !--  3:  Pressure weighting via "pressure_weights" variable in "lite" ncdf
                !--  4:  Straight pressure weighting
                !-------------------------------------------------------------------------

               !PRINT *,'VERTINTOPT:',VERTINTOPT

               SELECT CASE (VERTINTOPT)
 
            CASE(1)

            !------------------------------------------------------
            !--  This option weights the column vertical mean by using the OCO2
            !-- averaging kernel.  It interpolates this (and the prior CO2 profile)
            !--  onto the GCHEM pressure levels. Then uses something like:
            !-- XCO2_i = ak_i*press_i*co2_model_i + (1-ak_i)*press_i*aprior_co2
            !-- Note that PBL XCO2 is used to store the "no AK" case (press weight)
            !------------------------------------------------------

            !-- INTEGRATE A PRESSURE WEIGHTED AVG DOWN TO CLOUD TOP (or zero for surface)

            !-- Get local pressure profile


           DO VERTIND = 1, LLPAR
            GC_PRESS_IN(VERTIND) = REAL(GET_PCENTER(I,J,VERTIND),8) 
            GC_PRESS_WGHT(VERTIND) = GET_PEDGE(I,J,VERTIND)-
     &                                 GET_PEDGE(I,J,VERTIND+1)
            !- levels below OCO2 surf set to zero
            IF (GC_PRESS_IN(VERTIND) .gt. POCO2_PRESS(20,M)) THEN
                GC_PRESS_WGHT(VERTIND) = 0
            ENDIF

            ENDDO

   
           GC_PRESS_WGHT = GC_PRESS_WGHT / sum(GC_PRESS_WGHT)
           

        CALL interp_linear ( 1, 20, POCO2_PRESS(:,M),
     &                      PAK(:,M),47,
     &          GC_PRESS_IN(47:1:-1), AK_INTERP )

        CALL interp_linear ( 1, 20, POCO2_PRESS(:,M),
     &                      POCO2_PR(:,M),47,
     &          GC_PRESS_IN(47:1:-1), OCO2_PR_INTERP )
            !PRINT *,'M:',M
            !PRINT *,'GC_PRESS_IN:',GC_PRESS_IN
            !PRINT *,'GC_WGHT',GC_PRESS_WGHT
            !PRINT *,'GC_WGHT SUM:', SUM(GC_PRESS_WGHT)
            !PRINT *,'POCO2_PRESS(M,:):',POCO2_PRESS(:,M)
            !PRINT *,'POCO2_PR(,M):',POCO2_PR(:,M)
            !PRINT *,'AK20:',PAK(:,M)
            !PRINT *,'AKINTERP(1):',AK_INTERP
            !PRINT *,'PR_PROFILE_INTERP:',OCO2_PR_INTERP
            !PRINT *,'gc trace:',(STT(I,J,:,N)*Input_Opt%TCVV(N) /
            !&                     (State_Met%AD(I,J,:)*(1-
            !&                       State_Met%SPHU(I,J,:)*10**-3)))

             TMPSUM = 0.

                DO VERTIND = 1, LLPAR
                    PARTVARI(VERTIND) = AK_INTERP(LLPAR-VERTIND+1)*
     &               GC_PRESS_WGHT(VERTIND)*10**6*(STT(I,J,VERTIND,N) *
!    &               GC_PRESS_WGHT(VERTIND)*(STT(I,J,VERTIND,N) *
     &                  Input_Opt%TCVV(N) / 
     &                   (State_Met%AD(I,J,VERTIND)*(1-
     &                    State_Met%SPHU(I,J,VERTIND)*10**-3)))
     &              + GC_PRESS_WGHT(VERTIND) * 
     &                (1-AK_INTERP(LLPAR-VERTIND+1))*
     &               OCO2_PR_INTERP(LLPAR-VERTIND+1)

                  ENDDO

                      VARI(V) = SUM(PARTVARI) * ( 1E-6 )
                     IF ( VARI(V) < TINY ) VARI(V) = 0.d0

                    TMPSUM = PARTVARI

            !PRINT *,'test ak xco2:',PARTVARI

            !-- For full column XCO2 by pressure using MODEL CO2 for all levels

                DO VERTIND = 1, LLPAR
                    PARTVARI(VERTIND) = 
     &               GC_PRESS_WGHT(VERTIND)*(STT(I,J,VERTIND,N) *
     &                  Input_Opt%TCVV(N) /
     &                   (State_Met%AD(I,J,VERTIND)*(1-
     &                    State_Met%SPHU(I,J,VERTIND)*10**-3)))

                     !PRINT *,'PARTVARI(',VERTIND,'):',PARTVARI(VERTIND)
                     ENDDO

                       VARIPBL(V) = SUM(PARTVARI)
                     IF (  VARIPBL(V) < TINY )  VARIPBL = 0.d0

      !PRINT *,'test press xco2:',PARTVARI

      !PRINT *,'vert diff from AK:', TMPSUM - PARTVARI*10**6
      !PRINT *,'sum diff from AK:',sum(TMPSUM - PARTVARI*10**6)

               !VARIFREETROP(V) = sum(TMPSUM - PARTVARI*10**6)
               VARIFREETROP(V) = sum(PARTVARI)

       !stop('testing.....')

                 CASE (4)
   
                !------------------------------------------------------
                !--  This option is for pressure weighted, calcs
                !--  PBL, FREETROP and FULLCOLUMN XCO2
                !------------------------------------------------------

                !-- INTEGRATE A PRESSURE WEIGHTED AVG DOWN TO SURFACE

                     DO VERTIND = 1, LLPAR
                       IF (VERTIND .eq. LLPAR) THEN
                         PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * (GET_PEDGE(I,J,VERTIND) - 0)
                       ELSE
                         PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * (GET_PEDGE(I,J,VERTIND)-
     &                          GET_PEDGE(I,J,VERTIND+1))
                       ENDIF
                     ENDDO

                      VARI(V) = SUM(PARTVARI) / GET_PEDGE(I,J,1)
                     
                      IF ( VARI(V) < TINY ) VARI(V) = 0.d0

                 !---          PBL CALC

                 SUMPRESS = 0.0
                 DO VERTIND = 1, LLPAR

                   IF (VERTIND .eq. LLPAR)  THEN
                    PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                   Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                       * (GET_PEDGE(I,J,VERTIND) - 0) *
     &                           GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )

                   SUMPRESS = SUMPRESS + (GET_PEDGE(I,J,VERTIND) - 0) *
     &                          GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                   ELSE
                   PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                   Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * (GET_PEDGE(I,J,VERTIND)-
     &                          GET_PEDGE(I,J,VERTIND+1)) *
     &                          GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                   SUMPRESS = SUMPRESS + (GET_PEDGE(I,J,VERTIND)-
     &                          GET_PEDGE(I,J,VERTIND+1)) *
     &                          GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                  ENDIF

                ENDDO

                     VARIPBL(V) = SUM(PARTVARI) / SUMPRESS

                     IF ( VARIPBL(V) < TINY ) VARI(V) = 0.d0

                !---           FREE TROP CALC

               SUMPRESS = 0.0
               DO VERTIND = 1, LLPAR

                IF (VERTIND .eq. LLPAR)  THEN
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                  Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * (GET_PEDGE(I,J,VERTIND) - 0) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                 SUMPRESS = SUMPRESS + (GET_PEDGE(I,J,VERTIND) - 0) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                ELSE

                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                  Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * (GET_PEDGE(I,J,VERTIND)-
     &                          GET_PEDGE(I,J,VERTIND+1)) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                 SUMPRESS = SUMPRESS + (GET_PEDGE(I,J,VERTIND)-
     &                          GET_PEDGE(I,J,VERTIND+1)) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )
                ENDIF


               ENDDO

                     VARIFREETROP(V) = SUM(PARTVARI) / SUMPRESS
                     
                     IF ( VARIFREETROP(V) < TINY ) VARI(V) = 0.d0


            CASE(2)

            !------------------------------------------------------
            !--  This option is for straight full column weighting by KERN variable
            !--  Calcs PBL, FREETROP and FULLCOLUMN XCO2
            !------------------------------------------------------

            !-- INTEGRATE A PRESSURE WEIGHTED AVG DOWN TO CLOUD TOP (or zero for surface)

                     DO VERTIND = 1, LLPAR       
                       IF (VERTIND .eq. LLPAR) THEN
                          PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                       * KERN(VERTIND)
                       ELSE 
                         PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                   Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                    * KERN(VERTIND)
                       ENDIF
                     ENDDO

                      VARI(V) = SUM(PARTVARI)
                     IF ( VARI(V) < TINY ) VARI(V) = 0.d0


      !        PBL CALC
               SUMPRESS = 0.0
               DO VERTIND = 1, LLPAR

                IF (VERTIND .eq. LLPAR)  THEN
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                    * KERN(VERTIND) *
     &                   GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )

                 SUMPRESS = SUMPRESS + KERN(VERTIND) *
     &                          GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                ELSE
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                     * KERN(VERTIND) *
     &                     GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                 SUMPRESS = SUMPRESS +KERN(VERTIND) *
     &                      GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                ENDIF

               ENDDO

               VARIPBL(V) = SUM(PARTVARI) / SUMPRESS
               
               IF ( VARIPBL(V) < TINY ) VARI(V) = 0.d0

      !          FREE TROP CALC
                 SUMPRESS = 0.0
                DO VERTIND = 1, LLPAR

                IF (VERTIND .eq. LLPAR)  THEN
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                 Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                 SUMPRESS = SUMPRESS + KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                ELSE

                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                 SUMPRESS = SUMPRESS + KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )
                ENDIF

                ENDDO

                VARIFREETROP(V) = SUM(PARTVARI) / SUMPRESS
                IF ( VARIFREETROP(V) < TINY ) VARI(V) = 0.d0

            CASE(3)

            !------------------------------------------------------
            !--  This option is for column weighting by 'pressure_weights'
            !-- vertical weights given in netcdf file, by observation
            !--  This only calcs what the observations ask for, no xtra pbl/freetrop
            !------------------------------------------------------

            !-- INTEGRATE A PRESSURE WEIGHTED AVG DOWN TO CLOUD TOP (or zero for surface)

                     DO VERTIND = 1, LLPAR
                       IF (VERTIND .eq. LLPAR) THEN
                          PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                       * PPRESSWGHT(VERTIND,M)
                       ELSE
                         PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                   Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                    * PPRESSWGHT(VERTIND,M)
                       ENDIF
                     ENDDO

                      VARI(V) = SUM(PARTVARI)
                     IF ( VARI(V) < TINY ) VARI(V) = 0.d0


      !        PBL CALC
               SUMPRESS = 0.0
               DO VERTIND = 1, LLPAR

                IF (VERTIND .eq. LLPAR)  THEN
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                    * KERN(VERTIND) *
     &                   GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )

                 SUMPRESS = SUMPRESS + KERN(VERTIND) *
     &                          GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                ELSE
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                    Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                     * KERN(VERTIND) *
     &                     GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                 SUMPRESS = SUMPRESS +KERN(VERTIND) *
     &                      GET_FRAC_UNDER_PBLTOP( I, J, VERTIND )
                ENDIF

               ENDDO

               VARIPBL(V) = SUM(PARTVARI) / SUMPRESS

               IF ( VARIPBL(V) < TINY ) VARI(V) = 0.d0

      !          FREE TROP CALC
                 SUMPRESS = 0.0
                DO VERTIND = 1, LLPAR

                IF (VERTIND .eq. LLPAR)  THEN
                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                 Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                 SUMPRESS = SUMPRESS + KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                ELSE

                 PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
     &                Input_Opt%TCVV(N) / State_Met%AD(I,J,VERTIND))
     &                           * KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )

                 SUMPRESS = SUMPRESS + KERN(VERTIND) *
     &                       (1-GET_FRAC_UNDER_PBLTOP( I, J, VERTIND ) )
                ENDIF

                ENDDO

                VARIFREETROP(V) = SUM(PARTVARI) / SUMPRESS
                IF ( VARIFREETROP(V) < TINY ) VARI(V) = 0.d0


        END SELECT

                     ! Free pointer
                     NULLIFY( STT )
                     
                     DEALLOCATE(PARTVARI)
                  !--------------------------
                  ! Otherwise it's an error!
                  !--------------------------
                  CASE DEFAULT
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     WRITE( 6, '(a)' ) 'CO2SAT: Bad variable #!' 
                     WRITE( 6, '(a)' ) 'STOP in CO2SAT!'
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     CALL GEOS_CHEM_STOP

               END SELECT
            ENDDO

            ! Write data for the Mth plane point out to disk

           !print *,'M:',M,' vari:',vari(1)
           ! Not sure what the flush is for here - RM
            flush(6)
            
         
            ! Increment the observation count
            ! Output modeled CO2 to NetCDF file -RM
            CALL CO2SAT_MODEL_OUTPUT_WRITE(VARI(1), VARIPBL(1), 
     &              VARIFREETROP(1) )
            
!            WRITE_VARS_TO_FILE_NC( M, VARI(1), VARIPBL(1), 
!     &              VARIFREETROP(1) )

            ! Increment the record pointer
            PPOINT = PPOINT + 1

            ! Now closing file in SETUP_CO2SAT at beginning of next month
!            IF (M .eq. NPOINTS) THEN 
!              call ncGOSATOutputClose
!            END IF

         ENDIF
      ENDDO

      !      IF ( am_I_Root ) THEN
      !         DO N = 1, LLPAR
      ! print*,'fraction of PBL:',N,':',GET_FRAC_UNDER_PBLTOP( 32, 65, N )
      !         ENDDO
      !      ENDIF

      END SUBROUTINE CO2SAT

      
      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_CO2SAT
!
! !DESCRIPTION: Subroutine SET\_CO2SAT is used to pass values read in 
!  from the GEOS-Chem input file to "GOSAT\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CO2SAT( GOSAT, IN_FILE, OUT_FILE, VERTOPT )
!
! !INPUT PARAMETERS: 
!
      INTEGER,            INTENT(IN) :: VERTOPT    ! vertical integrate option
      LOGICAL,            INTENT(IN) :: GOSAT      ! Turn on GOSAT diag?
      CHARACTER(LEN=255), INTENT(IN) :: IN_FILE    ! Input file to read
      CHARACTER(LEN=255), INTENT(IN) :: OUT_FILE   ! Output file to write
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Save arguments to "shadow" module variables
      VERTINTOPT  = VERTOPT
      DO_CO2SAT    = GOSAT
      INFILENAME  = TRIM( IN_FILE  )
      OUTFILENAME = TRIM( OUT_FILE )

      PRINT *,'SETTING VERTINTOPT:',VERTINTOPT,' to VERTOPT:',VERTOPT
      END SUBROUTINE SET_CO2SAT

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_CO2SAT
!
! !DESCRIPTION: Subroutine CLEANUP\_CO2SAT deallocates all allocatable 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CO2SAT
!
! !USES:
!
      USE NCCOOARDSFORMAT,  ONLY : ncGOSATOutputClose
!! 
! !REVISION HISTORY: 
!  01 Jul 2001 - M. Evans    - Initial version
!  (1 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !CALL ncGOSATOutputClose

      IF ( ALLOCATED( PVAR   ) ) DEALLOCATE( PVAR   )
      IF ( ALLOCATED( PREAC  ) ) DEALLOCATE( PREAC  )
      IF ( ALLOCATED( PNAME  ) ) DEALLOCATE( PNAME  )
      IF ( ALLOCATED( PRRATE ) ) DEALLOCATE( PRRATE )
      IF ( ALLOCATED( PTYPE  ) ) DEALLOCATE( PTYPE  ) 
      IF ( ALLOCATED( PDATE  ) ) DEALLOCATE( PDATE  ) 
      IF ( ALLOCATED( PTIME  ) ) DEALLOCATE( PTIME  ) 
      IF ( ALLOCATED( PTAU   ) ) DEALLOCATE( PTAU   ) 
      IF ( ALLOCATED( PLAT   ) ) DEALLOCATE( PLAT   )
      IF ( ALLOCATED( PLON   ) ) DEALLOCATE( PLON   )
      IF ( ALLOCATED( PPRESS ) ) DEALLOCATE( PPRESS )
      IF ( ALLOCATED( PPRESSWGHT ) ) DEALLOCATE( PPRESSWGHT )
      IF ( ALLOCATED( PAK ) ) DEALLOCATE( PAK )
      IF ( ALLOCATED( POCO2_PR ) ) DEALLOCATE( POCO2_PR )
      IF ( ALLOCATED( POCO2_PRESS ) ) DEALLOCATE( POCO2_PRESS )
      IF ( ALLOCATED( PERR   ) ) DEALLOCATE( PERR   )
      IF ( ALLOCATED( PXCO2  ) ) DEALLOCATE( PXCO2  )
      IF ( ALLOCATED( PID    ) ) DEALLOCATE( PID  )

      END SUBROUTINE CLEANUP_CO2SAT

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: test_valid
!
! !DESCRIPTION: Subroutine TEST\_VALID tests to see if we are w/in the
!  tropopause, which is where SMVGEAR chemistry is done.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TEST_VALID( IND, PCHEM, JLOOP,I, J, L,
     &                        Input_Opt )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : JLOP
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE GRID_MOD,           ONLY : GET_XOFFSET
      USE GRID_MOD,           ONLY : GET_YOFFSET
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: IND         ! # of the flight track point
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object


!
! !OUTPUT PARAMETERS:
!
      LOGICAL, INTENT(OUT) :: PCHEM   ! =T if chemistry is done here
      INTEGER, INTENT(OUT) :: JLOOP   ! 1-D grid box index for SMVGEAR/KPP
      INTEGER, INTENT(OUT) :: I       ! GEOS-Chem longitude index
      INTEGER, INTENT(OUT) :: J       ! GEOS-Chem latitude index
      INTEGER, INTENT(OUT) :: L       ! GEOS-Chem level index
!
! !REVISION HISTORY:
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now use GET_PEDGE of "pressure_mod.f" to return the pressure at the
!        bottom edge of box (I,J,L), for hybrid grid. (dsa, bdf, bmy, 8/21/02)
!  (2 ) Since JLOP is not allocated for non-SMVGEAR runs, set PCHEM=F and
!        JLOOP=0 even if we are in the troposphere. (bmy, 7/18/03)
!  (3 ) Bug fix: add 0.5 in expression for I so that the rounding will
!        be done correctly.  Also make sure that I is computed correctly
!	 for points near the date line.  (bmy, 4/23/04)
!  (4 ) Now references ITS_A_FULLCHEM_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (5 ) Now references ITS_IN_THE_TROP from "tropopause_mod.f" (bmy, 8/22/05)
!  (6 ) Reference GET_XOFFSET and GET_YOFFSET from "grid_mod.f" and also
!        add for the case of nested-grid simulation (win, 7/28/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  08 Sep 2011 - L. Schiferl - Added correct definitions for I and J
!                              based on nested regions
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IL
      LOGICAL :: FOUND
      INTEGER :: I0, J0   ! (lds, 8/25/11)

      !=================================================================
      ! TEST_VALID begins here!
      !=================================================================

      ! We have not found a valid point
      FOUND = .FALSE.

      ! Added correct definitions for I and J based on nested regions
      ! (lds, 8/25/11)
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Get I corresponding to PLON(IND)
      I = INT( ( PLON(IND) + 180e+0 - (I0 * DISIZE) ) / DISIZE
     &     + 1.5e+0 )

      ! Handle date line correctly (bmy, 4/23/04)
      IF ( I > IIPAR ) I = I - IIPAR

      ! Get J corresponding to PLAT(IND)
      J = INT( ( PLAT(IND) +  90e+0 - (J0 * DJSIZE) ) / DJSIZE
     &     + 1.5e+0 )

      ! Get L corresponding to PRESS(IND)
      !L = 1
      !DO IL = 1, LLPAR
      !   IF ( GET_PEDGE(I,J,IL) <= PPRESS(IND) .AND..NOT. FOUND ) THEN
      !      L     = IL-1
      !      FOUND =.TRUE.
      !      EXIT
      !   ENDIF
      !ENDDO



      ! Error check: L must be 1 or higher
      IF ( L == 0 ) L = 1

      !=================================================================
      ! We only do full-chemistry in the troposphere
      !=================================================================
      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         ! JLOOP indicates if a box is in chemgrid (/=0) or not.
         JLOOP = JLOP(I,J,L)
         PCHEM = .TRUE.

      ELSE

         ! For non-SMVGEAR runs, JLOOP has no meaning so we give the
         ! false value
         JLOOP = 0
         PCHEM = .FALSE.

      ENDIF

      END SUBROUTINE TEST_VALID      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2SAT_OUTPUT_CREATE
!
! !DESCRIPTION: Subroutine CO2SAT\_OUTPUT\_CREATE creates output
!  file to write the CO2 satellite (GOSAT, OCO2) diagnostic.
!\\
!\\
! !INTERFACE:
!
        SUBROUTINE CO2SAT_OUTPUT_CREATE( filename, OBSCNT )
!
! !USES:
!
        USE m_netcdf_io_create
        USE m_netcdf_io_open
        USE m_netcdf_io_define
        
        include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
        CHARACTER (len = *)      ::  filename
        INTEGER,  INTENT(IN)     ::  OBSCNT

! !OUTPUT PARAMETERS:
!     The file ID SCO2_NCID will be set

! !REMARKS:
!  Code modified from ncGOSATOutputCreate. 
!  Opens and defines dimensions and variables for NetCDF output of GeosChem XCO2 output 
!  and satellite based observations of the same. Makes use of GeosChem NcdfUtils. 
!  Adheres to CF format rules as much as possible.
! 
! !REVISION HISTORY: 
!  15 April 2016 - R. McKeown    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        
        ! Units required for CF compliance
        character (len = *), parameter :: LAT_UNITS = "degrees_north"
        character (len = *), parameter :: LON_UNITS = "degrees_east"
        character (len = *), parameter :: TAU_UNITS = 
     &                          "Hours since 1/1/1985"

        ! Dimensions
        integer     :: obsDimID
        integer     :: dimIDa(1)
        integer     :: timeDimID, nTime

        integer     :: timeVID, levelVID
        integer     :: latVID, lonVID, soundingVID
        integer     :: cloudhVID, xco2VID, ucxco2VID
        integer     :: pblVID, ucpblVID, ftopVID, ucftopVID
        integer     :: mxco2VID, mpblVID, mftopVID

!           integer, parameter :: deflatelevel = 6

        ! Create file
        CALL Nccr_Wr (SCO2_NCID, filename)

        print *,'Opening/creating file ',trim(filename)

        ! Define Dimensions
        CALL NcDef_dimension(SCO2_NCID, "obs", OBSCNT, obsDimID)
        
        ! Create an array of the dimension id(s)
        dimIDa = (/ obsDimID /)

        ! Define Variables
        ! Time
        CALL NcDef_variable(SCO2_NCID, "time", NF_DOUBLE, 1, 
     &                      dimIDa, timeVID)
        CALL NcDef_var_attributes(SCO2_NCID, timeVID, "units", 
     &                      TAU_UNITS)
            
        ! Location information
        CALL NcDef_variable(SCO2_NCID, "latitude", NF_DOUBLE, 
     &                      1, dimIDa, latVID)
        CALL NcDef_var_attributes(SCO2_NCID, latVID, "units", 
     &                      LAT_UNITS)
        CALL NcDef_variable(SCO2_NCID, "longitude", NF_DOUBLE, 
     &                      1, dimIDa, lonVID)
        CALL NcDef_var_attributes(SCO2_NCID, lonVID, "units", 
     &                      LON_UNITS)
        CALL NcDef_variable(SCO2_NCID, "sounding_id", NF_INT, 
     &                      1, dimIDa, soundingVID)
            
        ! Satellite observations
        CALL NcDef_variable(SCO2_NCID, "xco2_cloudheight", 
     &                      NF_DOUBLE, 1, dimIDa, cloudhVID)
        CALL NcDef_variable(SCO2_NCID, "xco2_observed", NF_DOUBLE,
     &                      1, dimIDa, xco2VID)
        CALL NcDef_var_attributes(SCO2_NCID, xco2VID, "units",
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2_uncertainty", 
     &                      NF_DOUBLE, 1, dimIDa, ucxco2VID)
        CALL NcDef_var_attributes(SCO2_NCID, ucxco2VID, "units", 
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2pbl_observed", 
     &                      NF_DOUBLE, 1, dimIDa, pblVID)
        CALL NcDef_var_attributes(SCO2_NCID, pblVID, "units", 
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2pbl_uncertainty", 
     &                      NF_DOUBLE, 1, dimIDa, ucpblVID)
        CALL NcDef_var_attributes(SCO2_NCID, ucpblVID, "units", 
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2freetop_observed", 
     &                      NF_DOUBLE, 1, dimIDa, ftopVID)
        CALL NcDef_var_attributes(SCO2_NCID, ftopVID, "units", 
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2freetop_uncertainty", 
     &                      NF_DOUBLE, 1, dimIDa, ucftopVID)
        CALL NcDef_var_attributes(SCO2_NCID, ucftopVID, "units", 
     &                     "ppmv")
            
        ! Model output
        CALL NcDef_variable(SCO2_NCID, "xco2_model", 
     &                      NF_DOUBLE, 1, dimIDa, mxco2VID)
        CALL NcDef_var_attributes(SCO2_NCID, mxco2VID, "units", 
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2pbl_model", 
     &                      NF_DOUBLE, 1, dimIDa, mpblVID)
        CALL NcDef_var_attributes(SCO2_NCID, mpblVID, "units", 
     &                      "ppmv")
        CALL NcDef_variable(SCO2_NCID, "xco2freetop_model", 
     &                      NF_DOUBLE, 1, dimIDa, mftopVID)
        CALL NcDef_var_attributes(SCO2_NCID, mftopVID, "units", 
     &                      "ppmv")
            
        ! End of definitions
        CALL NcEnd_def(SCO2_NCID)
            
        END SUBROUTINE CO2SAT_OUTPUT_CREATE

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2SAT_OBS_OUTPUT_WRITE
!
! !DESCRIPTION: Subroutine CO2SAT\_OBS\_OUTPUT\_WRITE creates output
!  file to write the CO2 satellite (GOSAT, OCO2) diagnostic. All the 
!  output variables were read in SETUP_CO2SAT and are accessible to
!  all parts of the CO2_SAT module.
!\\
!\\
! !INTERFACE:
!
        SUBROUTINE CO2SAT_OBS_OUTPUT_WRITE(ID, LAT, LON,
     &      XCO2_OBS, XCO2_UNCERT, XCO2_CLOUDHEIGHT,
     &      XCO2FT_OBS, XCO2FT_UNCERT, XCO2PBL_OBS,
     &      XCO2PBL_UNCERT)
!
! !USES:
!
        USE m_netcdf_io_create,   ONLY : Ncdo_Sync
        USE m_netcdf_io_write
        
        include "netcdf.inc"
!
! !INPUT PARAMETERS:
        INTEGER, dimension(:) :: ID
        REAL*8,  dimension(:) :: LAT, LON 
        REAL*8,  dimension(:) :: XCO2_OBS,XCO2_UNCERT
        REAL*8,  dimension(:) :: XCO2PBL_OBS,XCO2PBL_UNCERT
        REAL*8,  dimension(:) :: XCO2FT_OBS,XCO2FT_UNCERT
        REAL*8,  dimension(:) :: XCO2_CLOUDHEIGHT
        
! 
! !REVISION HISTORY: 
!  15 April 2016 - R. McKeown    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        INTEGER               :: strt1d(1), cnt1d(1)

        strt1d = (/ OBSNUM /)
        cnt1d = (/ NPOINTS /)

        CALL Ncwr(LAT, SCO2_NCID, "latitude", strt1d, cnt1d)
        CALL Ncwr(LON, SCO2_NCID, "longitude", strt1d, cnt1d)
            
        CALL Ncwr(XCO2_CLOUDHEIGHT, SCO2_NCID, "xco2_cloudheight", 
     &              strt1d, cnt1d)
        CALL Ncwr(XCO2_OBS, SCO2_NCID, "xco2_observed", 
     &              strt1d, cnt1d)
        CALL Ncwr(XCO2_UNCERT, SCO2_NCID, "xco2_uncertainty", 
     &              strt1d, cnt1d)
        CALL Ncwr(XCO2PBL_OBS, SCO2_NCID, "xco2pbl_observed", 
     &              strt1d, cnt1d)
        CALL Ncwr(XCO2PBL_UNCERT, SCO2_NCID, "xco2pbl_uncertainty", 
     &              strt1d, cnt1d)
        CALL Ncwr(XCO2FT_OBS, SCO2_NCID, "xco2freetop_observed", 
     &              strt1d, cnt1d)
        CALL Ncwr(XCO2FT_UNCERT, SCO2_NCID, "xco2freetop_uncertainty", 
     &              strt1d, cnt1d)

        ! Sync the file to make sure data is saved
        !CALL Ncdo_Sync (SCO2_NCID)

        END SUBROUTINE CO2SAT_OBS_OUTPUT_WRITE
    
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2SAT_MODEL_OUTPUT_WRITE
!
! !DESCRIPTION: Subroutine CO2SAT\_MODEL\_OUTPUT\_WRITE creates output
!  file to write the CO2 satellite (GOSAT, OCO2) diagnostic. All the 
!  output variables were read in SETUP_CO2SAT and are accessible to
!  all parts of the CO2_SAT module.
!\\
!\\
! !INTERFACE:
!
        SUBROUTINE CO2SAT_MODEL_OUTPUT_WRITE(V, V2, V3 )
!
! !USES:
!
        USE TIME_MOD,           ONLY : GET_TAU
        USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
        USE m_netcdf_io_create, ONLY : Ncdo_Sync
        USE m_netcdf_io_write
        
        include "netcdf.inc"
!
! !INPUT PARAMETERS:
        REAL*8                ::  V, V2, V3
        
! 
! !REVISION HISTORY: 
!  15 April 2016 - R. McKeown    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        REAL*8,  dimension(1) :: arrV, arrV2, arrV3, arrTime
        
        INTEGER               :: strt1d(1), cnt1d(1)
        
        ! If somehow the file id is lost or doesn't exist
        IF ( SCO2_NCID .EQ. 0 ) THEN 
            CALL GEOS_CHEM_STOP
        ENDIF
        strt1d = (/ OBSNUM /)
        cnt1d = (/ 1 /)        
        arrV  = V * 10**6
        arrV2 = V2 * 10**6
        arrV3 = V3 * 10**6
        
        arrTime = GET_TAU()

        CALL Ncwr(arrTime, SCO2_NCID, "time", strt1d, cnt1d)
        CALL Ncwr(arrV, SCO2_NCID, "xco2_model", strt1d, cnt1d)
        CALL Ncwr(arrV2, SCO2_NCID, "xco2pbl_model", strt1d, cnt1d)
        CALL Ncwr(arrV3, SCO2_NCID, "xco2freetop_model", 
     &            strt1d, cnt1d)
        
        !CALL Ncdo_Sync (SCO2_NCID)

        OBSNUM = OBSNUM + 1    

        END SUBROUTINE CO2SAT_MODEL_OUTPUT_WRITE 
        
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COUNT_OBS_IN_MONTH
!
! !DESCRIPTION: Subroutine COUNT\_OBS\_IN\_MONTH gets the number of
!               observations for each day in a month in order to
!               preallocate the space needed in the output file.
!               Each daily file is opened and the dimension length
!               corresponding to the observations is extracted.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COUNT_OBS_IN_MONTH(START_NYMD, NHMS, OBSCNT)
!
! !USES:
!
        USE FILE_MOD,   ONLY : FILE_EXISTS
        USE FILE_MOD,   ONLY : IU_FILE
        USE TIME_MOD,   ONLY : EXPAND_DATE
        USE TIME_MOD,   ONLY : GET_DAY
        USE ERROR_MOD,     ONLY : ERROR_STOP
        USE m_netcdf_io_close
        USE m_netcdf_io_open
        USE m_netcdf_io_get_dimlen

      
        include "netcdf.inc"
!
! !INPUT PARAMETERS:
        INTEGER,  INTENT(IN)  :: START_NYMD
        INTEGER,  INTENT(IN)  :: NHMS
        INTEGER,  INTENT(OUT) :: OBSCNT

! 
! !REVISION HISTORY: 
!  15 April 2016 - R. McKeown    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        INTEGER              :: DAY, START_DAY
        INTEGER              :: CURRENT_NYMD
        INTEGER              :: DIM_LEN
        INTEGER              :: NCID
        CHARACTER(LEN=255)   :: CURRENT_FILE
        CHARACTER(LEN=500)   :: MSG
        
! Copy file names to local variables and
        !  replace and date & time tokens in the file names
        
        OBSCNT = 0
        
        ! In case the simulation starts in the middle of a month
        START_DAY  = GET_DAY()

        DO DAY = START_DAY, 31
            CURRENT_FILE   = INFILENAME
            
            ! This should only increment the day in NYMD
            CURRENT_NYMD = START_NYMD + (DAY - START_DAY)
            CALL EXPAND_DATE( CURRENT_FILE,  CURRENT_NYMD, NHMS )

            IF ( FILE_EXISTS( CURRENT_FILE ) ) THEN 
                CALL Ncop_Rd( NCID, trim(CURRENT_FILE) )

                CALL NCget_DimLen(NCID, 'sounding_id', DIM_LEN)
                OBSCNT = OBSCNT + DIM_LEN
            ELSE IF (DAY < 31) THEN
                CONTINUE    ! File is missing, try the next day
            ELSE
                RETURN
            ENDIF
        ENDDO 
        
        END SUBROUTINE COUNT_OBS_IN_MONTH 
        
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2SAT_OUTPUT_CLOSE
!
! !DESCRIPTION: Subroutine CO2SAT\_OUTPUT\_CLOSE closes the NetCDF
!               output file associated with a CO2 satellite track
!\\
!\\
! !INTERFACE:
!
        SUBROUTINE CO2SAT_OUTPUT_CLOSE()
!
! !USES:
!        
        USE m_netcdf_io_close

            print *, "closing CO2 satellite output file", SCO2_NCID
            CALL Nccl( SCO2_NCID )

        END SUBROUTINE CO2SAT_OUTPUT_CLOSE
!EOC

      END MODULE CO2SAT_MOD
