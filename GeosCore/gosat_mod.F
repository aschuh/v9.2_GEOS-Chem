!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GOSAT_mod
!
! !DESCRIPTION: Module GOSAT\_MOD contains variables and routines which 
!  are used to "fly" a plane through the GEOS-Chem model simulation.  This is 
!  useful for comparing model results with aircraft observations. 
!\\
!\\
! !INTERFACE:
!
      MODULE GOSAT_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
     
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: ARCHIVE_RXNS_FOR_PF
      PUBLIC  :: CLEANUP_GOSAT
      PUBLIC  :: GOSAT
      PUBLIC  :: SETUP_GOSAT
      PUBLIC  :: SET_GOSAT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: AN_SETUP
      PRIVATE :: INIT_GOSAT
      PRIVATE :: NOY_SETUP
      PRIVATE :: READ_VARIABLES
      PRIVATE :: READ_POINTS
      PRIVATE :: RO2_SETUP
      PRIVATE :: TEST_VALID
      PRIVATE :: WRITE_VARS_TO_FILE
      PRIVATE :: WRITE_LITE_OUTPUT
!
! !REMARKS:
!  The quantities that are saved to disk by the GOSAT diagnostic were
!  requested by GEOS-Chem users.  If you would like to save out a new quantity,
!  then you will have to make your own modifications in this module.
!
! !REVISION HISTORY:
!  (1 ) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/21/02)
!  (2 ) Now reference AD from "dao_mod.f".  Now also references "error_mod.f".
!        (bmy, 10/15/02)
!  (3 ) Bug fix: replace missing commas in FORMAT statement (bmy, 3/23/03)
!  (4 ) Now references "time_mod.f". (bmy, 3/27/03)
!  (5 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  (6 ) Bug fix: use NAMEGAS instead of NAMESPEC (lyj, bmy, 7/9/03)
!  (7 ) Bug fix: avoid referencing JLOP for non-SMVGEAR runs (bmy, 7/18/03)
!  (8 ) Bug fix: Use T instead of T3 for GMAO temperature.  Also replace
!        NAMESPEC w/ NAMEGAS in RO2_SETUP.  Now locate reordered rxn 
!        numbers for SMVGEAR II.(tdf, mje, bmy, 8/1/03)
!  (9 ) Now print out N2O5 hydrolysis rxn as a special case.   Also rename
!        output file. (bmy, 8/8/03)
!  (10) Changed "DAO" to "GMAO" for met field variable names.  Now can save 
!        aerosol optical depths.  Bug fix in TEST_VALID. (bmy, 4/23/03)
!  (11) Now references "tracer_mod.f" (bmy, 7/20/04)
!  (12) Bug fix in READ_VARIABLES (1/7/05)
!  (13) Modified the plane flight diagnostic so that it writes output files
!        for each day where flight track files are defined. (bmy, 3/24/05)
!  (14) Minor bug fix in ARCHIVE_RXNS_FOR_PF (bmy, 5/20/05)
!  (15) Now split AOD's into column AOD's and AOD's below plane.  Also scale
!        AOD's to 400nm. (bmy, 10/25/05)
!  (16) Bug fixes in READ_VARIABLES (bmy, 10/16/06)
!  (17) Bug fix in GOSAT (cdh, bmy, 12/12/06)
!  (18) Bug fix in RO2_SETUP (tmf, bmy, 4/23/07)
!  (19) Set very small values to zero.  (tmf, 1/7/09)
!  (20) Add new RO2 species according to 'globchem.dat' (tmf, 1/7/09) 
!  (21) Make sure we have 3 spaces in the exponential format (phs, 7/13/09)
!  (22) Output the grid cell indexes (kjw, 8/18/09)
!  (23) Add AN and NOy species. (fp, 3/10/10)
!  (24) Now scale AODs to wavelength specified in jv_spec_aod.dat(clh, 5/14/09)
!  29 Jul 2011 - R. Yantosca - Now also archive MERRA SEAICExx fields
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
!=========================================================================
! Derived type for Lite Output
!=========================================================================
      TYPE, PUBLIC :: LITE_DATA_OBJECT
      !----------------------------------------------------------------------
      ! netcdf integer fields
      !----------------------------------------------------------------------
      INTEGER*8,  POINTER :: SOUNDING_ID(:)

      !-----------------------------------------------------------------------
      ! netcdf real4  fields
      !-----------------------------------------------------------------------

      REAL,  POINTER       :: XCO2_MODEL  (: )
      REAL,  POINTER       :: XCO2FREETROP_MODEL(: )
      REAL,  POINTER       :: XCO2PBL_MODEL(: )
      REAL,  POINTER       :: XCO2_OBSERVED(: )
      REAL,  POINTER       :: XCO2_UNCERTAINTY(: )
      REAL,  POINTER       :: XCO2PBL_OBSERVED(: )
      REAL,  POINTER       :: XCO2PBL_UNCERTAINTY(: )
      REAL,  POINTER       :: XCO2FREETROP_OBSERVED(: )
      REAL,  POINTER       :: XCO2FREETROP_UNCERTAINTY(: )
      REAL,  POINTER       :: LATITUDE     (:  )
      REAL,  POINTER       :: LONGITUDE     (:  )
      REAL,  POINTER       :: TIME       (:  )
      REAL,  POINTER       :: XCO2_CLOUDHEIGHT (: )

      END TYPE LITE_DATA_OBJECT

! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! DO_GOSAT       : Turn on the GOSAT diagnostic? (T/F)
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
      LOGICAL                        :: DO_GOSAT

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

      ! Logical unit numbers (formerly in file_mod.F)
      INTEGER                        :: IU_FILE
      INTEGER                        :: IU_GOSAT

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
! !IROUTINE: setup_GOSAT
!
! !DESCRIPTION: Subroutine SETUP\_GOSAT reads information from the 
!  input file in order to initialize the GOSAT diagnostic.  Also 
!  calls INIT\_GOSAT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETUP_GOSAT ( am_I_Root , Input_Opt,RC)
!
! !USES:
!
      USE FILE_MOD,   ONLY : FILE_EXISTS
      USE FILE_MOD,   ONLY : IOERROR
      USE FILE_MOD,   ONLY : IU_FILE
      USE TIME_MOD,   ONLY : EXPAND_DATE
      USE TIME_MOD,   ONLY : GET_NYMD
      USE TIME_MOD,   ONLY : GET_NHMS
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE NCCOOARDSFORMAT, ONLY :ncGOSATOutputCreate,ncGOSATOutputClose
      USE NCCOOARDSFORMAT, ONLY :READGOSAT_ID,READGOSAT_REAL8 
      USE NCCOOARDSFORMAT, ONLY : READGOSAT_REAL4
      USE m_netcdf_io_close
      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_get_dimlen
      USE m_do_err_out
      USE m_netcdf_io_checks
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?

! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?

! !REMARKS:
!  For SMVGEAR or KPP simulations, the call to SETUP_GOSAT is made from 
!  routine "chemdr.f", after the "chem.dat" file is read.  This is necessary 
!  since we have to reference the SMVGEAR rxn rate and species numbers.
!                                                                             .
!  For offline simulations, the call to SETUP_GOSAT can be made
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
      INTEGER            :: NCID, DIM_LEN
      INTEGER            :: st2d(2),ct2d(2),st1d(1),ct1d(1)
      INTEGER, dimension(:), allocatable :: ID4
      INTEGER*8, dimension(:), allocatable :: ID
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
      ! SETUP_GOSAT begins here!
      !=================================================================

      ! Assume that there is flight data for today
      DO_GOSAT = .TRUE.

      ! Find a free file LUN
      IU_FILE = findFreeLun()

      ! Get date & time
      NYMD  = GET_NYMD()
      NHMS  = GET_NHMS()

      ! Copy file names to local variables
      INF   = INFILENAME
      OUTF  = OUTFILENAME
      
      ! Replace any date & time tokens in the file names
      CALL EXPAND_DATE( INF,  NYMD, NHMS )
      CALL EXPAND_DATE( OUTF, NYMD, NHMS )

      PRINT *,'INF:',INF
      ! If we can't find a flighttrack file for today's date, return
      IF ( .not. FILE_EXISTS( INF ) ) THEN 
         DO_GOSAT = .FALSE.
         RETURN
      ENDIF

      ! Echo info
      ! Echo info
      IF ( am_I_Root ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a)' ) 'P L A N E   F L I G H T   D I A G N O S T I C'
       WRITE( 6, 100   ) TRIM( INF )
  100  FORMAT( /, 'SETUP_GOSAT: Reading ',a )
       WRITE( 6, '(a)' )  
      ENDIF

      ! Compute # of species and # of points & allocate arrays
      CALL INIT_GOSAT( am_I_Root )

      ! Return if there are no flight track points for today
      IF ( NPOINTS == 0 ) THEN
         IF ( am_I_Root ) THEN
           WRITE( 6, '(a)' ) 'No flight track found for today!'
         ENDIF
         DO_GOSAT = .FALSE.
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

      !=================================================================
      ! Open today's plane.log file and write file header
      !=================================================================

      ! Close previously-opened file
      !CLOSE( IU_GOSAT )

      ! Open new file
      !OPEN( IU_GOSAT, FILE=TRIM( OUTF ), STATUS='UNKNOWN', IOSTAT=IOS )

      ! Error check
      !IF ( IOS /= 0 ) THEN
      !   CALL IOERROR( IOS, IU_GOSAT, 'setup_GOSAT:1' )
      !ENDIF

      ! Write header
      !WRITE( IU_GOSAT, 110 ) 'POINT', 'TYPE', 'YYYYMMDD', 'HHMM',
      !&  'LAT', 'LON', 'PRESS', 'XCO2', 'ERROR', ( PNAME(I), I=1,NPVAR )

      ! FORMAT string
      ! Replace 95 by 200 FP
! 110  FORMAT( A5,X,A5,X,A8,X,A4,X,A7,X,A7,X,A7,X,95(a10,x) )
! 110  FORMAT( A5,X,A5,X,A8,X,A4,X,A7,X,A7,X,A7,X,A7,X,200(a10,x) )

       ALLOCATE( ID4( NPOINTS ) )
       ALLOCATE( ID( NPOINTS ) )
       ALLOCATE( XCO2_CLOUDHEIGHT (NPOINTS) )
       ALLOCATE( XCO2_OBS (NPOINTS) )
       ALLOCATE( XCO2_UNCERT (NPOINTS) )
       ALLOCATE( XCO2FT_OBS (NPOINTS) )
       ALLOCATE( XCO2FT_UNCERT (NPOINTS) )
       ALLOCATE( XCO2PBL_OBS (NPOINTS) )
       ALLOCATE( XCO2PBL_UNCERT (NPOINTS) )
       ALLOCATE( XCO2_LATITUDE (NPOINTS) )
       ALLOCATE( XCO2_LONGITUDE (NPOINTS) )
       
      ! CALL READGOSAT_ID( trim (INF), NPOINTS, ID)
      !PRINT *,'ID1:',ID

      !  CALL READGOSAT_REAL4( trim (INF), 'xco2_raw', 'Retrieval',
      !&  ACOS_XCO2, NPOINTS)

      CALL Ncop_Rd( NCID, trim(INF) )

      !CALL Ncget_Unlim_Dimlen (NCID, DIM_LEN)
      CALL NCget_DimLen(NCID, 'sounding_id',DIM_LEN)

      !PRINT *,'DIM:',DIM_LEN
      !PRINT *,'npoints:',npoints
      !-- ID already grabbed earlier
      !ID = 0
      !CALL readGOSAT_ID(trim(INF), NPOINTS, ID)
      !print *,'done w/ redgosat'
      !CALL Ncrd(ID, ncid, 'sounding_id', (/ 1 /), (/ DIM_LEN /))
      !ID = INT8(ID4)

      st1d = (/ 1 /)
      ct1d = (/ NPOINTS /)


      !CALL Ncrd(XCO2_OBS, ncid, 'xco2_observed',st1d,ct1d)
      IF (VERTINTOPT .eq. 3) THEN
      CALL READGOSAT_REAL8( trim (INF), 'xco2_observed', '',
     &   XCO2_OBS, NPOINTS)
      ENDIF
      IF (VERTINTOPT .eq. 1) THEN
      CALL READGOSAT_REAL8( trim (INF), 'xco2', '',
     &   XCO2_OBS, NPOINTS)
      ENDIF

      !PRINT *,'after first pull...'

      !CALL Ncrd(XCO2_UNCERT, ncid, 'xco2_uncertainty',st1d,
      !&  ct1d)
      CALL READGOSAT_REAL8( trim (INF), 'xco2_uncertainty', '',
     &   XCO2_UNCERT, NPOINTS)

      surfheight_TF = Ncdoes_Var_Exist (ncid,'surface_height')

      IF(surfheight_TF) THEN
      CALL READGOSAT_REAL8( trim (INF), 'surface_height', '',
     &   XCO2_CLOUDHEIGHT, NPOINTS)
      !CALL Ncrd(XCO2_CLOUDHEIGHT, ncid, 'surface_height',st1d,
      !&   ct1d)
      ELSE
           XCO2_CLOUDHEIGHT = 0
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'xco2freetrop_observed')) THEN
      CALL READGOSAT_REAL8( trim (INF), 'xco2freetrop_observed', 
     &   '', XCO2FT_OBS, NPOINTS)
      !CALL Ncrd(XCO2FT_OBS, ncid, 'xco2freetrop_observed', st1d, 
      !&   ct1d)
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'xco2freetrop_uncertainty')) THEN
      CALL READGOSAT_REAL8( trim (INF), 'xco2freetrop_uncertainty',
     &   '', XCO2FT_UNCERT, NPOINTS)
      !      CALL Ncrd(XCO2FT_UNCERT, ncid, 'xco2freetrop_uncertainty',
      !&  st1d,ct1d)
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'xco2pbl_observed')) THEN
      CALL READGOSAT_REAL8( trim (INF), 'xco2pbl_observed',
     &   '', XCO2PBL_OBS, NPOINTS)
      !      CALL Ncrd(XCO2PBL_OBS, ncid, 'xco2pbl_observed', st1d, 
      !&  ct1d)
      ENDIF

      IF(Ncdoes_Var_Exist (ncid, 'xco2pbl_uncertainty')) THEN
      CALL READGOSAT_REAL8( trim (INF), 'xco2pbl_uncertainty',
     &   '', XCO2PBL_UNCERT, NPOINTS)
      !      CALL Ncrd(XCO2PBL_UNCERT, ncid, 'xco2pbl_uncertainty', 
      !&  st1d, ct1d)
      ENDIF

      CALL READGOSAT_REAL8( trim (INF), 'longitude',
     &   '', XCO2_LONGITUDE, NPOINTS)

      CALL READGOSAT_REAL8( trim (INF), 'latitude',
     &   '', XCO2_LATITUDE, NPOINTS)

      ! CALL Ncrd(XCO2_LONGITUDE, ncid, 'longitude',st1d,
      !&   ct1d)

      ! CALL Ncrd(XCO2_LATITUDE, ncid, 'latitude', st1d,
      !&  ct1d)

!      PRINT *,'XCO2_LAT:',XCO2_LATITUDE

      CALL Nccl( NCID )


       CALL ncGOSATOutputCreate(trim(OUTF),NPOINTS,XCO2_CLOUDHEIGHT,
     &   XCO2_LATITUDE,
     &   XCO2_LONGITUDE, ID,XCO2_OBS,
     &  XCO2_UNCERT,XCO2PBL_OBS,XCO2PBL_UNCERT,XCO2FT_OBS,XCO2FT_UNCERT,
     &  "xco2_model","xco2pbl_model",
     &    "xco2freetrop_model")

       DEALLOCATE( ID )
       DEALLOCATE( XCO2_OBS)
       DEALLOCATE( XCO2_CLOUDHEIGHT)
       DEALLOCATE( XCO2_UNCERT)
       DEALLOCATE( XCO2PBL_OBS)
       DEALLOCATE( XCO2PBL_UNCERT)
       DEALLOCATE( XCO2FT_OBS)
       DEALLOCATE( XCO2FT_UNCERT)
       DEALLOCATE( XCO2_LATITUDE )
       DEALLOCATE( XCO2_LONGITUDE )

       !CALL ncGOSATOutputClose

      END SUBROUTINE SETUP_GOSAT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_lite_output
!
! !DESCRIPTION: Subroutine READ\_POINTS reads the information (ID, date, time,
!  lat, lon, pressure) for each measurement listed in the input file, and
!  sorts these into the appropriate module variables.
!\\
!\\
! !INTERFACE:
!

      SUBROUTINE WRITE_LITE_OUTPUT(fname,lite_data,
     &  NPOINTS )
!
!******************************************************************************
!  Subroutine GET_I6_FIELDS_2 is a wrapper for routine READ_I6.  This routine
!  calls READ_I6 properly for reading I-6 fields from GEOS-3, GEOS-4, GEOS_5,
!  or GCAP met data sets every 6 hours during a GEOS-CHEM run.
!  (bmy, 6/23/03, 1/16/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-3 fields to be read from disk
!
!  NOTES:
!  (1 ) Now modified for GEOS-5 and GCAP met fields
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Now read TO3 and TTO3 for GEOS-5 (bmy, 1/16/07)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!******************************************************************************
!
      ! References to F90 modules
      USE netcdf
      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_read
      USE m_netcdf_io_write
      USE m_netcdf_io_get_dimlen
      USE m_do_err_out
      USE m_netcdf_io_create
      USE m_netcdf_io_define

      ! Arguments
      TYPE(LITE_DATA_OBJECT), INTENT(IN) :: lite_data  ! Meteorology State object
      character* (*),         INTENT(IN)    :: fname
      INTEGER,                INTENT(IN)    :: NPOINTS

      ! Local Variables
      integer                               :: NCID, DIM_LEN
      !integer,  dimension(:), allocatable   :: SAMPLE_STRATEGY

      ! For netCDF file I/O
      INTEGER              :: idLevels,idSoundingDim
      !INTEGER             :: idObs,    idTracer,   idTname, idChar100
      !INTEGER             :: idRgNmLen,  idGridInd
      INTEGER             :: fId
      INTEGER             ::  vId
      !INTEGER             :: ct1d(1),  ct2d(2), ct3d(3), ct4d(4)
      !INTEGER             :: st1d(1),  st2d(2), st3d(3), st4d(4)
      INTEGER             :: var1(1),var2(1)

      !-- Pointer to CT input object ID
      !CHARACTER, POINTER :: P_ID(:,:)

      !=========================================================================
      ! Create the netCDF file
      !=========================================================================

      ! Echo info
      WRITE( 6, '(a)' ) '=== Begin netCDF file creation test ==='
      PRINT *,'WRITING ',NPOINTS,' observations'
      CALL NcCr_Wr( fId, TRIM(fname) )

      ! Turn filling off
      !CALL NcSetFill( fId, NF_NOFILL, omode )

      !=========================================================================
      ! Define the dimensions
      !=========================================================================

      ! Sounding dimension
      !WRITE( 6, '(a)' ) 'Writing obs  (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'sounding_id',  NPOINTS, idSoundingDim)
      !CALL NcDef_Unlim_Dimension( fId, 'obs',idObs)

      ! Levels dimension
      !WRITE( 6, '(a)' ) 'Writing tracer   (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'levels',  47,  idLevels )

      !=========================================================================
      ! Define the variables and variable attributes
      ! for COARDS compliance and GAMAP compliance
      !=========================================================================
      CALL NcDef_Glob_Attributes( fId, 'Title',       'NcdfUtilities' )
      CALL NcDef_Glob_Attributes( fId, 'History',     '24 Jan 2011'   )
      CALL NcDef_Glob_Attributes( fId, 'Conventions', 'COARDS'        )
      CALL NcDef_Glob_Attributes( fId, 'Model',       'GEOS4'         )
      CALL NcDef_Glob_Attributes( fId, 'Nlayers',     '55'            )
      CALL NcDef_Glob_Attributes( fId, 'Start_Date',  '20110101'      )
      CALL NcDef_Glob_Attributes( fId, 'Start_Time',  '00:00:00.0'    )
      CALL NcDef_Glob_Attributes( fId, 'End_Date',    '20110101'      )
      CALL NcDef_Glob_Attributes( fId, 'End_Time',    '23:59:59.0'    )
      CALL NcDef_Glob_Attributes( fId, 'Delta_Lon',   '5'             )
      CALL NcDef_Glob_Attributes( fId, 'Delta_Lat',   '4'             )
      CALL NcDef_Glob_Attributes( fId, 'Delta_time',  '000000'        )
      CALL NcDef_Glob_Attributes( fId, 'Format',      'netCDF-3'      )
      CALL NcDef_Glob_Attributes( fId,'valid_range',(/-1e15,+1e15 /)  )

      !=========================================================================
      ! Define the variables and variable attributes
      !=========================================================================

      ! Time index array (hardwire date to 2011/01/01)
      !var1    = (/ idTime /)
      vId     = -1
      !units   = 'minutes since 2011-01-01 00:00:00 GMT'
      !delta_t = '0000-00-00 00:00:00'
      !begin_d = '20110101'
      !begin_t = '000000'
      !incr    = '000000'
      !CALL NcDef_Variable      ( fId, 'time', NF90_INT,  1, var1, vId           )
      !CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'time'           )
      !CALL NcDef_Var_Attributes( fId, vId, 'units',          TRIM( units   )  )
      !CALL NcDef_Var_Attributes( fId, vId, 'delta_t',        TRIM( delta_t )  )
      !CALL NcDef_Var_Attributes( fId, vId, 'begin_date',     TRIM( begin_d )  )
      !CALL NcDef_Var_Attributes( fId, vId, 'begin_time',     TRIM( begin_t )  )
      !CALL NcDef_Var_Attributes( fId, vId, 'time_increment', TRIM( incr    )  )

      ! Define sounding_id variable
      vId  = vId + 1
      var1 = (/ idSoundingDim /)
      CALL NcDef_Variable( fId, 'sounding_id', NF90_INT64,1, var1, vId )
      !CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Sounding ID' )

      ! Define latitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'latitude', NF90_FLOAT,1, var1, vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'longitude', NF90_FLOAT,1, var1, vId )

      ! Define latitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2_model', NF90_FLOAT,1, var1, vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2freetrop_model', NF90_FLOAT,1, 
     &     var1, vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2pbl_model', NF90_FLOAT,1,
     &     var1, vId )

      ! Define latitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2_observed',NF90_FLOAT,1,var1,vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2_uncertainty', NF90_FLOAT,1,
     &     var1, vId )

      ! Define latitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2freetrop_observed', NF90_FLOAT,1, 
     &    var1, vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable(fId,'xco2freetrop_uncertainty', NF90_FLOAT,1,
     &     var1, vId )

      ! Define latitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2pbl_observed', NF90_FLOAT,1, 
     &   var1, vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2pbl_uncertainty', NF90_FLOAT,1,
     &     var1, vId )

      ! Define longitude variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'time', NF90_DOUBLE,1,
     &     var1, vId )

      ! Define bottom variable
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'xco2_surfheight', NF90_DOUBLE,1,
     &     var1, vId )

      var2 = (/ 47 /)     

      CALL NcDef_Variable( fId, 'levels', NF90_INT, 1, var2, vId+1 )

      CALL NcDef_Variable( fId, 'pressure_levels', NF90_DOUBLE, 1, 
     &    var2, vId+2 )
      CALL NcDef_Variable( fId, 'pressure_weights', NF90_DOUBLE, 1, 
     &    var2, vId+3 )

      !=========================================================================
      ! %%% END OF DEFINITION SECTION %%%
      ! %%% NOW WRITE DATA TO FILE    %%%
      !=========================================================================
      CALL NcEnd_def( fId )

      !P_ID     =>   Object_Lite_Data%ID

      !PRINT *,'FIRST POINT:',P_ID(1:100,1)
      !CALL NcWr(P_ID(1:100,1:NPOINTS) ,fId, 'obspack_id',
      !&   (/ 1,1 /),(/ 100, NPOINTS /) )

      CALL NcWr(lite_data%latitude(1:NPOINTS),fId,
     &    'latitude',(/ 1 /),(/ NPOINTS /) )

      CALL NcWr(lite_data%longitude(1:NPOINTS),fId,
     &    'longitude',(/ 1 /),(/ NPOINTS /) )

      CALL NcWr(lite_data%xco2_model(1:NPOINTS),fId,
     &    'xco2_model',(/ 1 /),(/ NPOINTS /) )

!      CALL NcWr(lite_data%xco2freetrop_model(1:NPOINTS),fId,
!     &    'xco2freetrop_model',(/ 1 /),(/ NPOINTS /) )

      !=========================================================================
      ! Close the netCDF file
      !=========================================================================
      CALL NcCl( fId )

      !NULLIFY(P_ID)

      ! Echo info
      WRITE( 6, '(a)' ) '=== End netCDF file creation test ==='

      ! Return to calling program
      END SUBROUTINE WRITE_LITE_OUTPUT
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_Lite_Data_Object
!
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_Lite_Data_Object( am_I_Root, NOBS,Object_Lite_Data
     &    , RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod                         ! Error codes

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      INTEGER,        INTENT(IN)    :: NOBS        ! # number of obs in file
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(Lite_Data_Object), INTENT(INOUT) :: Object_Lite_Data   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!  19 Oct 2012 - R. Yantosca - Now pass all dimensions as arguments
!  23 Oct 2012 - R. Yantosca - Now allocate QI, QL fields
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  16 Nov 2012 - R. Yantosca - Now zero all fields after allocating
!  27 Nov 2012 - R. Yantosca - Now allocate SUNCOS fields (IM,JM)
!  12 Dec 2012 - R. Yantosca - Now allocate the IREG, ILAND, IUSE fields
!  13 Dec 2012 - R. Yantosca - Now allocate the XLAI, XLAI2 fields
!  07 Mar 2013 - R. Yantosca - Now allocate PF*LSAN, PF*CU fields properly
!                              for GEOS-5.7.x met (they are edged)
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: LX

      ! Assume success
       RC = GIGC_SUCCESS

       PRINT *,'Initializing Lite Data Object....with ',NOBS,' obs'

      !=======================================================================
      ! Allocate Fields
      !=======================================================================
      ALLOCATE(Object_Lite_Data%LATITUDE (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%LATITUDE   = 0

      ALLOCATE(Object_Lite_Data%LONGITUDE (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%LONGITUDE   = 0

      ALLOCATE(Object_Lite_Data%XCO2_MODEL (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2_MODEL   = 0

      ALLOCATE(Object_Lite_Data%XCO2FREETROP_MODEL (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2FREETROP_MODEL   = 0

      ALLOCATE(Object_Lite_Data%XCO2PBL_MODEL (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2PBL_MODEL   = 0

      ALLOCATE(Object_Lite_Data%XCO2_OBSERVED (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2_OBSERVED   = 0

      ALLOCATE(Object_Lite_Data%XCO2_UNCERTAINTY (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2_UNCERTAINTY   = 0

      ALLOCATE(Object_Lite_Data%XCO2PBL_OBSERVED (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2PBL_OBSERVED   = 0

      ALLOCATE(Object_Lite_Data%XCO2PBL_UNCERTAINTY (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2PBL_UNCERTAINTY   = 0

      ALLOCATE(Object_Lite_Data%XCO2FREETROP_OBSERVED (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2FREETROP_OBSERVED   = 0

      ALLOCATE(Object_Lite_Data%XCO2FREETROP_UNCERTAINTY (NOBS),STAT=RC)
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%XCO2FREETROP_UNCERTAINTY   = 0

      ALLOCATE(Object_Lite_Data%TIME (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%TIME   = 0

      ALLOCATE(Object_Lite_Data%SOUNDING_ID (NOBS), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      Object_Lite_Data%SOUNDING_ID   = 0

      END SUBROUTINE Init_Lite_Data_Object
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Carbontracker_In
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_Lite_Data( am_I_Root, Object_Lite_Data, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod                         ! Error codes
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(Lite_Data_Object), INTENT(INOUT) :: Object_Lite_Data   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!  23 Oct 2012 - R. Yantosca - Now deallocate QI, QL fields
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  19 Nov 2012 - R. Yantosca - Segregate DEALLOCATE statements w/ #ifdefs
!                              for each met field data product type
!  27 Nov 2012 - R. Yantosca - Now deallocate the SUNCOS fields
!  12 Dec 2012 - R. Yantosca - Now deallocate the IREG, ILAND, IUSE fields
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Return success
      RC = GIGC_SUCCESS


      IF ( ASSOCIATED( Object_Lite_Data%LATITUDE       )) THEN 
            DEALLOCATE( Object_Lite_Data%LATITUDE       )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%LONGITUDE       )) THEN
            DEALLOCATE( Object_Lite_Data%LONGITUDE      )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%TIME       )) THEN
            DEALLOCATE( Object_Lite_Data%TIME       )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%SOUNDING_ID       )) THEN
            DEALLOCATE( Object_Lite_Data%SOUNDING_ID     )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2_MODEL       )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2_MODEL       )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2FREETROP_MODEL       )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2FREETROP_MODEL      )
      ENDIF
      IF ( ASSOCIATED( Object_Lite_Data%XCO2PBL_MODEL       )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2PBL_MODEL       )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2_OBSERVED       )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2_OBSERVED      )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2_UNCERTAINTY       )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2_UNCERTAINTY     )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2FREETROP_OBSERVED    )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2FREETROP_OBSERVED      )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2FREETROP_UNCERTAINTY )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2FREETROP_UNCERTAINTY     )
      ENDIF
      IF ( ASSOCIATED( Object_Lite_Data%XCO2PBL_OBSERVED       )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2PBL_OBSERVED      )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2PBL_UNCERTAINTY     )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2PBL_UNCERTAINTY     )
      ENDIF

      IF ( ASSOCIATED( Object_Lite_Data%XCO2_CLOUDHEIGHT     )) THEN
            DEALLOCATE( Object_Lite_Data%XCO2_CLOUDHEIGHT     )
      ENDIF

      END SUBROUTINE Cleanup_Lite_Data
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_variables
!
! !DESCRIPTION: Subroutine READ\_VARIABLES reads the list of variables 
!  (SMVGEAR/KPP chemical species, SMVGEAR/KPP rxn rates, GMAO met fields, or 
!  GEOS-Chem tracers) to be printed out and sorts the information into the 
!  appropriate module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_VARIABLES( am_I_Root, Input_Opt )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,   ONLY : IU_FILE
      USE FILE_MOD,   ONLY : IOERROR

      USE CMN_SIZE_MOD     ! Size parameters
      USE COMODE_LOOP_MOD  ! NAMEGAS, NSPEC

      USE GIGC_Input_Opt_Mod, ONLY : OptInput

!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY:
!  30 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (2 ) Bug fix: replace missing commas in FORMAT statement (bmy, 3/23/03)
!  (3 ) Bug fix: replace NAMESPEC w/ NAMEGAS for SMVGEAR II (lyj, bmy, 7/9/09)
!  (4 ) Now locate reordered rxn numbers for SMVGEAR II. (mje, bmy, 8/1/03)
!  (5 ) Now flag N2O5 hydrolysis rxn as a special case (bmy, 8/8/03)
!  (6 ) Changed variable name prefix "DAO" to "GMAO".  Also added aerosol
!        optical depths w/ tracer offset 2000. (bmy, 4/23/04)
!  (7 ) Now references N_TRACERS & ITS_A_FULLCHEM_SIM from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (8 ) Bug fix: extract tracer # when reading rxn rates (bmy, 1/7/05)
!  (9 ) Now computes column AOD's and AOD's below plane (bmy, 10/24/05)
!  (10) We need to trim NAMEGAS before comparing to LINE so that comparisons 
!        for species like "O3" will work.  Also set NCS=NCSURBAN at the top
!        of the subroutine, to avoid out of bounds error. (dbm, bmy, 10/16/06)
!  29 Jul 2011 - R. Yantosca - Also search for MERRA SEAICExx met fields
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object

      LOGICAL             :: IS_FULLCHEM
      INTEGER             :: M, N, NUM, R, IK, IOS
      CHARACTER(LEN=255)  :: LINE

      !=================================================================
      ! READ_VARIABLES begins here!
      !=================================================================

      ! Reset NCS to NCSURBAN for safety's sake (dbm, bmy, 10/16/06)
      NCS = NCSURBAN

      ! Test if this is a fullchem run
      IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM

      ! Read four lines of header
      DO N = 1, 4 
         READ( IU_FILE, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:1')
      ENDDO

      ! Read in the number of species to be output
      READ( IU_FILE, '(i3)', IOSTAT=IOS ) NPVAR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:2' )

      ! Read in a separation line
      READ( IU_FILE, '(a)', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:3' )

      ! Echo to stdout
      IF ( am_I_Root ) THEN     
         WRITE( 6, '(a)' ) '   #    Species       PVAR'
         WRITE( 6, '(a)' ) '-----------------------------'
      ENDIF

      !=================================================================
      ! Sort variables by type; assign indices to PVAR, PREAC arrays
      ! NOTE: Variables for which PVAR(N) = 0 will be skipped!
      !=================================================================

      ! Zero reaction counter
      R = 0

      ! Loop over all variables
      DO N = 1, NPVAR

         ! Read each line
         READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:4')

         ! Save the name of each variable into the global PNAME array
         PNAME(N) = LINE(1:10)

         ! We are searching for a ...
         SELECT CASE ( LINE(1:4) ) 

            !===========================================================
            ! GEOS-CHEM tracer: listed as "TRA_001", etc.
            ! PVAR offset: 100000
            !===========================================================
            CASE ( 'TRA_' )

               ! Extract tracer # from the string
               READ( LINE(5:14), '(i10)' ) NUM

               ! Make sure the tracer # is valid!
               IF ( NUM < 0 .or. NUM > Input_Opt%N_TRACERS ) THEN  
                  IF ( am_I_Root ) THEN
                    WRITE( 6, 100 ) TRIM( LINE )
 100                FORMAT( 'TRACER ', i4, ' is out of range!' )
                    WRITE( 6, '(a)' ) 'STOP in SETUP_GOSAT!'
                    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  ENDIF
                  CALL GEOS_CHEM_STOP
               ENDIF

               ! Save in PVAR -- add offset of 100000
               PVAR(N) = 100000 + NUM

            !===========================================================
            ! GMAO met field: listed as "GMAO_TEMP", etc.
            ! PVAR offset: 1000
            !===========================================================
            CASE ( 'GMAO' )
               
               IF ( LINE == 'GMAO_TEMP'  ) PVAR(N) = 1001
               IF ( LINE == 'GMAO_ABSH'  ) PVAR(N) = 1002
               IF ( LINE == 'GMAO_SURF'  ) PVAR(N) = 1003 
               IF ( LINE == 'GMAO_PSFC'  ) PVAR(N) = 1004
               IF ( LINE == 'GMAO_UWND'  ) PVAR(N) = 1005      
               IF ( LINE == 'GMAO_VWND'  ) PVAR(N) = 1006 
               IF ( LINE == 'GMAO_IIEV'  ) PVAR(N) = 1007
               IF ( LINE == 'GMAO_JJEV'  ) PVAR(N) = 1008
               IF ( LINE == 'GMAO_LLEV'  ) PVAR(N) = 1009
               IF ( LINE == 'GMAO_ICE00' ) PVAR(N) = 1100
               IF ( LINE == 'GMAO_ICE10' ) PVAR(N) = 1101
               IF ( LINE == 'GMAO_ICE20' ) PVAR(N) = 1102
               IF ( LINE == 'GMAO_ICE30' ) PVAR(N) = 1103
               IF ( LINE == 'GMAO_ICE40' ) PVAR(N) = 1104
               IF ( LINE == 'GMAO_ICE50' ) PVAR(N) = 1105
               IF ( LINE == 'GMAO_ICE60' ) PVAR(N) = 1106
               IF ( LINE == 'GMAO_ICE70' ) PVAR(N) = 1107
               IF ( LINE == 'GMAO_ICE80' ) PVAR(N) = 1108
               IF ( LINE == 'GMAO_ICE90' ) PVAR(N) = 1109

            !===========================================================
            ! Column aerosol optical depths (same order as for FAST-J)
            ! PVAR offset: 2000
            !===========================================================
            CASE ( 'AODC' )
            
               IF ( LINE == 'AODC_SULF'  ) PVAR(N) = 2001
               IF ( LINE == 'AODC_BLKC'  ) PVAR(N) = 2002
               IF ( LINE == 'AODC_ORGC'  ) PVAR(N) = 2003
               IF ( LINE == 'AODC_SALA'  ) PVAR(N) = 2004
               IF ( LINE == 'AODC_SALC'  ) PVAR(N) = 2005   

            !===========================================================
            ! Aerosol optical depths below the plane
            ! (same order as for FAST-J)  PVAR offset: 3000
            !===========================================================
            CASE ( 'AODB' )
            
               IF ( LINE == 'AODB_SULF'  ) PVAR(N) = 3001
               IF ( LINE == 'AODB_BLKC'  ) PVAR(N) = 3002
               IF ( LINE == 'AODB_ORGC'  ) PVAR(N) = 3003
               IF ( LINE == 'AODB_SALA'  ) PVAR(N) = 3004
               IF ( LINE == 'AODB_SALC'  ) PVAR(N) = 3005 

            !===========================================================
            ! SMVGEAR rxn rate: listed as "REA_001", etc.
            ! PVAR offset: 10000
            !===========================================================
            CASE ( 'REA_' )

               ! Skip if not SMVGEAR!
               IF ( IS_FULLCHEM ) THEN 
               
                  ! Increment rxn counter
                  R = R + 1

                  IF ( TRIM( LINE ) == 'REA_O1D' ) THEN

                     ! O1D is a special rxn, give it offset of 20000
                     PVAR(N)  = 20000
                     PREAC(R) = 20000

                  ELSE IF ( TRIM( LINE ) == 'REA_N2O5' ) THEN

                     ! N2O5 hydrolysis is another special rxn
                     ! give it an offset of 21000
                     PVAR(N)  = 21000
                     PREAC(R) = 21000

                  ELSE
                     !==================================================
                     ! NOTE: the reaction numbers listed in smv2.log 
                     ! aren't really used to index SMVGEAR II rxns.  The 
                     ! rxns get reordered.  Find the right rxn number, 
                     ! which is stored in NOLDFNEW.  We assume only one 
                     ! chemistry scheme. (mje, bmy, 8/1/03)
                     !==================================================

                     ! Extract tracer # from the string
                     READ( LINE(5:14), '(i10)' ) NUM

                     ! Initialize
                     PVAR(N)  = -999
                     PREAC(R) = -999

                     ! Search for proper rxn number
                     DO IK = 1, NMTRATE 

                        ! Offset other reaction rates by 10000
                        IF ( NOLDFNEW(IK,1) == NUM ) THEN 
                           PVAR(N)  = 10000 + IK
                           PREAC(R) = 10000 + IK
                           EXIT
                        ENDIF
                     ENDDO

                     ! Stop w/ error 
                     IF ( PVAR(N) == -999 ) THEN 
                        IF ( am_I_Root ) THEN
                           WRITE (6,*) 'Cant match up reaction number'
                           WRITE (6,*) NUM
                           WRITE (6,*) 'Is it the second line of the'
                           WRITE (6,*) 'Three body reaction'
                           WRITE (6,*) 'Stopping'
                        ENDIF
                        CALL GEOS_CHEM_STOP
                     ENDIF
                  ENDIF
               ENDIF
               
            !===========================================================
            ! SMVGEAR chem species: listed as "O3", "C2H6", etc.
            ! PVAR offset: 0
            !===========================================================
            CASE DEFAULT

               ! Skip if not SMVGEAR!
               IF ( IS_FULLCHEM ) THEN

                  ! Loop over all SMVGEAR species -- 
                  ! match w/ species as read from disk
                  DO M = 1, NSPEC(NCS)
                     IF ( TRIM( NAMEGAS(M) ) == TRIM( LINE ) ) THEN
                        PVAR(N) = M
                        EXIT
                     ENDIF
                  ENDDO

                  ! Special flag for RO2 species
                  IF ( TRIM( LINE ) == 'RO2' ) PVAR(N) = 999

                  ! Special flag for AN species FP
                  IF ( TRIM( LINE ) == 'AN' ) PVAR(N) = 998

                  ! Special flag for NOy species FP
                  IF ( TRIM( LINE ) == 'NOy' ) PVAR(N) = 997

                  ! Error check
                  IF ( PVAR(N) == 0 ) THEN 
                     WRITE( 6, '(a)' ) 'ERROR: invalid species!'
                     WRITE( 6, 110   ) TRIM( LINE )
 110                 FORMAT( 'Species ', a, ' not found!' )
                     WRITE( 6, '(a)' ) 'STOP in GOSAT!'
                     CALL GEOS_CHEM_STOP
                  ENDIF
               ENDIF

         END SELECT
      
         ! Echo species names/numbers to screen
         WRITE( 6, 120 ) N, TRIM( LINE ), PVAR(N)
 120     FORMAT( i4, 1x, a10, 1x, i10 )

      ENDDO

      END SUBROUTINE READ_VARIABLES
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
      USE FILE_MOD,        ONLY : IOERROR
      USE NCCOOARDSFORMAT, ONLY : READGOSAT_REAL8, READGOSAT_ID
      USE NCCOOARDSFORMAT, ONLY : READGOSAT_REAL4
      USE JULDAY_MOD,        ONLY : CALDATE
      USE TIME_MOD,      ONLY : YMD_EXTRACT
      USE m_netcdf_io_read
      USE m_netcdf_io_close
      USE m_netcdf_io_open
      USE m_netcdf_io_checks
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
      ! We have previously computed NPOINTS in INIT_GOSAT
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

      CALL READGOSAT_REAL4(trim(INF),'latitude','',LAT,NPOINTS)
      CALL READGOSAT_REAL4(trim(INF),'longitude','',LON,
     & NPOINTS)
      CALL Ncop_Rd( NCID, trim(INF) )

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

      !-- LEAVE INIT AT ZERO UNTIL FIXED
      !PRESS = 0.
!      CALL READGOSAT_REAL4(trim(INF),'xco2_cloudheight','',PRESS,
!     & NPOINTS)

      !print *,'LAT:',LAT
      CALL READGOSAT_ID( trim (INF), NPOINTS, ID)
      !CALL READGOSAT_REAL4( trim (INF), 'xco2', '',
      !&  XCO2, NPOINTS)
      !CALL READGOSAT_REAL4( trim (INF), 'xco2_uncertainty', '',
      !&  ERR, NPOINTS)
      CALL READGOSAT_REAL8( trim (INF), 'time', '',
     &  TIM, NPOINTS)

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
! !DESCRIPTION: Subroutine GOSAT saves concentrations to disk at 
!  locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GOSAT(am_I_Root,Input_Opt,State_Met,State_Chm,RC)
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

      !USE CMN_SIZE_MOD, ONLY : LLPAR
      !USE CMN_FJ_MOD,   ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ
      !USE CMN_FJX_MOD,        ONLY : ODAER, QAA, QAA_AOD, ODMDUST
      !USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
      USE JV_CMN_MOD
      USE COMODE_LOOP_MOD    ! CSPEC, etc.
      USE NCCOOARDSFORMAT,  ONLY  : WRITE_VARS_TO_FILE_NC
      USE NCCOOARDSFORMAT,  ONLY  : ncGOSATOutputClose
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
!  (7 ) Now return if DO_GOSAT = .FALSE. (bmy, 3/24/05)
!  (8 ) Now compute column AOD's and AOD's below plane.  Also now scale
!        AOD's to 400nm. (bmy, 10/24/05)
!  (9 ) Bug fix: exit if PTAU(M) == PTAUE, so that we write out on the next !
!        GOSAT timestep (cdh, bmy, 12/12/06)
!  (10) Change GOSAT output time step. (ccc, 8/27/09)
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
      ! GOSAT begins here!
      !=================================================================

      !print *,'AD12:',AD12(40,65,:)

      ! Return if there is no flighttrack data for today
      PRINT *,'DO_GOSAT:',DO_GOSAT
      IF ( .not. DO_GOSAT ) RETURN

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
            CALL WRITE_VARS_TO_FILE( M, VARI )      

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
 100        FORMAT( '     - GOSAT: Archived ',a5,1x,i8.8,1x,i4.4 )

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
                     !print *,'in gosat tracer plane loop'
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
                     !-- Temp kernel weighting for GOSAT code




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
                     WRITE( 6, '(a)' ) 'GOSAT: Bad variable #!' 
                     WRITE( 6, '(a)' ) 'STOP in GOSAT!'
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     CALL GEOS_CHEM_STOP

               END SELECT
            ENDDO

            ! Write data for the Mth plane point out to disk

           !print *,'M:',M,' vari:',vari(1)
            flush(6)
            CALL WRITE_VARS_TO_FILE_NC( M, VARI(1), VARIPBL(1), 
     &              VARIFREETROP(1) )

            ! Increment the record pointer
            PPOINT = PPOINT + 1

            IF (M .eq. NPOINTS) THEN 
              call ncGOSATOutputClose
            END IF

         ENDIF
      ENDDO

      !      IF ( am_I_Root ) THEN
      !         DO N = 1, LLPAR
      ! print*,'fraction of PBL:',N,':',GET_FRAC_UNDER_PBLTOP( 32, 65, N )
      !         ENDDO
      !      ENDIF

      END SUBROUTINE GOSAT
!EOC


      subroutine interp_linear ( m, data_num, t_data, p_data,interp_num,
     &          t_interp, p_interp )

!*****************************************************************************80
!
!! INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
!
!  Discussion:
!
!    From a space of M dimensions, we are given a sequence of
!    DATA_NUM points, which are presumed to be successive samples
!    from a curve of points P.
!
!    We are also given a parameterization of this data, that is,
!    an associated sequence of DATA_NUM values of a variable T.
!    The values of T are assumed to be strictly increasing.
!
!    Thus, we have a sequence of values P(T), where T is a scalar,
!    and each value of P is of dimension M.
!
!    We are then given INTERP_NUM values of T, for which values P
!    are to be produced, by linear interpolation of the data we are given.
!
!    Note that the user may request extrapolation.  This occurs whenever
!    a T_INTERP value is less than the minimum T_DATA or greater than the
!    maximum T_DATA.  In that case, linear extrapolation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
!    independent variable at the sample points.  The values of T_DATA
!    must be strictly increasing.
!
!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
!    dependent variables at the sample points.
!
!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
!    at which interpolation is to be done.
!
!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
!    independent variable at the interpolation points.
!
!    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
!    values of the dependent variables at the interpolation points.
!
       implicit none

       integer ( kind = 4 ) data_num
       integer ( kind = 4 ) m
       integer ( kind = 4 ) interp_num

       integer ( kind = 4 ) interp
       integer ( kind = 4 ) left
       real ( kind = 8 ) p_data(m,data_num)
       real ( kind = 8 ) p_interp(m,interp_num)
       !logical r8vec_ascends_strictly
       integer ( kind = 4 ) right
       real ( kind = 8 ) t
       real ( kind = 8 ) t_data(data_num)
       real ( kind = 8 ) t_interp(interp_num)

       if ( .not. r8vec_ascends_strictly ( data_num, t_data ) ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
         write ( *, '(a)' ) 
     & '  Independent variable array T_DATA is not strictly increasing.'
         stop 1
       end if

       do interp = 1, interp_num

         t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
         call r8vec_bracket ( data_num, t_data, t, left, right )

         p_interp(1:m,interp) = 
     &     ( ( t_data(right) - t                ) * p_data(1:m,left)   
     &  + (                 t - t_data(left) ) * p_data(1:m,right) ) 
     &  / ( t_data(right)     - t_data(left) )

       end do

       return
      end subroutine interp_linear 

      subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
       implicit none

       integer ( kind = 4 ) n

       integer ( kind = 4 ) i
       integer ( kind = 4 ) left
       integer ( kind = 4 ) right
       real ( kind = 8 ) x(n)
       real ( kind = 8 ) xval

       do i = 2, n - 1

         if ( xval < x(i) ) then
           left = i - 1
           right = i
           return
         end if

        end do

       left = n - 1
       right = n
     
       return
      end subroutine r8vec_bracket

       function r8vec_ascends_strictly ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!    Notice the effect of entry number 6 in the following results:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
!      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
!
!      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS_STRICTLY, is TRUE if the
!    entries of X strictly ascend.
!
       implicit none

       integer ( kind = 4 ) n

       integer ( kind = 4 ) i
       logical r8vec_ascends_strictly
       real ( kind = 8 ) x(n)

       do i = 1, n - 1
         if ( x(i+1) <= x(i) ) then
           r8vec_ascends_strictly = .false.
           return
         end if
       end do

       r8vec_ascends_strictly = .true.

        return
      end function r8vec_ascends_strictly

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
! !IROUTINE: write_vars_to_file
!
! !DESCRIPTION: Subroutine WRITE\_VARS\_TO\_FILE writes the values of all 
!  the variables for a given flight track point to the output file.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_VARS_TO_FILE( IND, VARI )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: IND           ! # of the flight track point
      REAL*8,  INTENT(IN) :: VARI(NPVAR)   ! Values to print to file
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) The max line length for output seems to be 1024 characters.  Adjust
!        MAXVARS accordingly so that we don't exceed this. (bmy, 7/8/02)
!  (2 ) Now do not write file header -- this is now done in subroutine
!        SETUP_GOSAT at the start of each day (bmy, 3/25/05)
!  (3 ) Bug fix: make sure we have 3 spaces in exponential (phs, 7/13/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, IOS

      !=================================================================
      ! WRITE_VARS_TO_FILE begins here!
      !=================================================================

      ! Write data to file
      WRITE( IU_GOSAT, 110, IOSTAT=IOS ) 
     &     IND, PTYPE(IND), INT( PDATE(IND) ), INT( PTIME(IND) ),
     &     PLAT(IND), PLON(IND), PPRESS(IND), PXCO2(IND), PERR(IND), 
     &        ( VARI(I), I=1,NPVAR )
      
      ! Format string
      ! FP replace 95 by 200
! 110  FORMAT( I5,   X, A5,   X, I8.8, X, I4.4, X, 
!     &        F7.2, X, F7.2, X, F7.2, X, 95(es11.3e3,x) )
 110  FORMAT(I5,X,A5,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.0,X,
     &       F8.3,X,F7.2,X,200(es15.7e3,x))


      !PRINT *,'PXCO2:',PXCO2(IND)
      WRITE( *, 110, IOSTAT=IOS )
     &     IND, PTYPE(IND), INT( PDATE(IND) ), INT( PTIME(IND) ),
     &     PLAT(IND), PLON(IND), PPRESS(IND), PXCO2(IND), PERR(IND),
     &        ( VARI(I), I=1,NPVAR )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_GOSAT,'write_vars_to_file:1')

      ! Flush the file to disk
      CALL FLUSH( IU_GOSAT )

      END SUBROUTINE WRITE_VARS_TO_FILE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: archive_rxns_for_PF
!
! !DESCRIPTION: Subroutine ARCHIVE\_RXNS\_FOR\_PF is called from "calcrate.f" 
!  to pass reaction rates from the SMVGEAR solver for the GOSAT 
!  diagnostic. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ARCHIVE_RXNS_FOR_PF( JO1D, N2O5 )
!
! !USES:
!
      USE COMODE_MOD, ONLY : IXSAVE
      USE COMODE_MOD, ONLY : IYSAVE
      USE COMODE_MOD, ONLY : IZSAVE
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP

      USE CMN_SIZE_MOD                     ! Size parameters
      USE COMODE_LOOP_MOD                  ! RRATE, JLOOPLO, KBLOOP
      USE CMN_DIAG_MOD                     ! ND40 switch
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)  :: JO1D(KBLOOP)  ! JO1D photolysis rate [1/s]
      REAL*8, INTENT(IN)  :: N2O5(KBLOOP)  ! N2O5 hydrolysis rate [molec/cm3/s]
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now avoid overflow/underflow errors in PRATE (bmy, 7/8/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run (bmy, 10/15/02)
!  (3 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  (4 ) Now also pass N2O5 hydrolysis rxn rate array via the arg list.  
!        Also bug fix: replace TMP with RATE in under/overflow checking
!        for JO1D and N2O5. (bmy, 8/8/03)      
!  (5 ) Bug fix: Replace with DO_GOSAT since this variable is reset to either T 
!        or F each day depending on whether there is plane flight data 
!        available (bmy, 5/20/05)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: KLOOP, JLOOP, V, R, I, J, L 
      REAL*8             :: RATE
!
! !DEFINED PARAMETERS:
!
      ! Smallest, largest REAL*4 #'s representable on this machine
      REAL*4, PARAMETER  :: SMALLEST=TINY(1e0), LARGEST=HUGE(1e0)

      !=================================================================
      ! ARCHIVE_RXNS_FOR_PF begins here!
      !=================================================================
      IF ( DO_GOSAT ) THEN

         ! Loop over SMVGEAR reactions
         DO R = 1, NPREAC

            ! Test SMVGEAR rxn number
            SELECT CASE ( PREAC(R) ) 

               !-----------------------
               ! All except JO1D, N2O5
               !-----------------------
               CASE( 10000:19999 )
                         
                  ! Store rate in PRRATE
                  DO KLOOP = 1, KTLOOP
                     JLOOP = JLOOPLO + KLOOP
                     RATE  = RRATE(KLOOP,PREAC(R)-10000)

                     ! Avoid overflow/underflow
                     IF ( RATE < SMALLEST ) RATE = 0e0
                     IF ( RATE > LARGEST  ) RATE = LARGEST

                     PRRATE(JLOOP,R) = RATE
                  ENDDO

               !-----------------------
               ! JO1D photolysis rxn 
               !-----------------------
               CASE ( 20000 )

                  ! Store rate in PRATE
                  DO KLOOP = 1, KTLOOP
                     JLOOP = JLOOPLO + KLOOP
                     RATE  = JO1D(KLOOP)

                     ! Avoid overflow/underflow
                     IF ( RATE < SMALLEST ) RATE = 0e0
                     IF ( RATE > LARGEST  ) RATE = LARGEST 

                     PRRATE(JLOOP,R) = RATE
                  ENDDO
                 
               !-----------------------
               ! N2O5 hydrolysis rxn
               !-----------------------
               CASE ( 21000 ) 

                  ! Store rate in PRATE
                  DO KLOOP = 1, KTLOOP
                     JLOOP = JLOOPLO + KLOOP
                     RATE  = N2O5(KLOOP)

                     ! Avoid overflow/underflow
                     IF ( RATE < SMALLEST ) RATE = 0e0
                     IF ( RATE > LARGEST  ) RATE = LARGEST 

                     PRRATE(JLOOP,R) = RATE
                  ENDDO

               !-----------------------
               ! Error: invalid rxn
               !-----------------------
               CASE DEFAULT
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  WRITE( 6, '(a)' ) 'ERROR -- Invalid SMVGEAR rxn #!'
                  WRITE( 6, '(a)' ) 'STOP in ARCHIVE_RXNS_FOR_PF!'
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  CALL GEOS_CHEM_STOP

            END SELECT
         ENDDO
      ENDIF

      END SUBROUTINE ARCHIVE_RXNS_FOR_PF
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_GOSAT
!
! !DESCRIPTION: Subroutine SET\_GOSAT is used to pass values read in 
!  from the GEOS-Chem input file to "GOSAT\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_GOSAT( GOSAT, IN_FILE, OUT_FILE, VERTOPT )
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
      DO_GOSAT    = GOSAT
      INFILENAME  = TRIM( IN_FILE  )
      OUTFILENAME = TRIM( OUT_FILE )

      PRINT *,'SETTING VERTINTOPT:',VERTINTOPT,' to VERTOPT:',VERTOPT
      END SUBROUTINE SET_GOSAT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_GOSAT
!
! !DESCRIPTION: Subroutine INIT\_GOSAT reads the input file to compute 
!  the number of variables and flight track points to print out.  Also 
!  allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GOSAT( am_I_Root )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IOERROR
      USE NCCOOARDSFORMAT,  ONLY : readGOSATlength

      USE CMN_SIZE_MOD    ! Size Parameters
      USE COMODE_LOOP_MOD ! ITLOOP
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
!  (3 ) INIT_GOSAT is now called each day but the arrays are only
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
      INTEGER           :: N, AS, IOS
      CHARACTER(LEN=20) :: LINE

      !=================================================================
      ! INIT_GOSAT begins here!
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
      print *,'start readgosat'
      CALL readGOSATlength( trim (INF) , NPOINTS)
      print *,'done readgosat'

      ! If there are no flight-track points then just return
      IF ( NPOINTS < 1 ) THEN
         DO_GOSAT = .FALSE.
         RETURN
      ENDIF

      ! Make sure NPOINTS is less than MAXPOINTS
      IF ( NPOINTS > MAXPOINTS ) THEN
         WRITE( 6, '(a)') 'NPOINTS exceeds maximum allowed value!'
         WRITE( 6, '(a)') 'STOP in INIT_GOSAT (GOSAT_mod.f)'
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

      END SUBROUTINE INIT_GOSAT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_GOSAT
!
! !DESCRIPTION: Subroutine CLEANUP\_GOSAT deallocates all allocatable 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GOSAT
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

      END SUBROUTINE CLEANUP_GOSAT
!EOC
      END MODULE GOSAT_MOD
