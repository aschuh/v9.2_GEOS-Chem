!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CT_mod
!
! !DESCRIPTION: Module CT\_MOD contains variables and routines which 
!  are used to "fly" a plane through the GEOS-Chem model simulation.  This is 
!  useful for comparing model results with aircraft observations. 
!\\
!\\
! !INTERFACE:
!
      MODULE CT_MOD
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
      PUBLIC  :: CLEANUP_CT
      PUBLIC  :: CT
      PUBLIC  :: SETUP_CT
      PUBLIC  :: SET_CT
      PUBLIC  :: CT_FINISH_OUTPUT
      PUBLIC  :: SETUP_CT_OCFL_OUT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: AN_SETUP
      PRIVATE :: INIT_CT
      PRIVATE :: NOY_SETUP
      PRIVATE :: READ_VARIABLES
      PRIVATE :: READ_POINTS
      PRIVATE :: RO2_SETUP
      PRIVATE :: TEST_VALID
      PRIVATE :: WRITE_VARS_TO_FILE
      PRIVATE :: CCM2_Z2_sigma
!
! !REMARKS:
!  The quantities that are saved to disk by the CT diagnostic were
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
!  (17) Bug fix in CT (cdh, bmy, 12/12/06)
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
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! DO_CT       : Turn on the CT diagnostic? (T/F)
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
      LOGICAL                        :: DO_CT

      ! Parameters
      INTEGER,           PARAMETER   :: MAXVARS   = 200
      INTEGER,           PARAMETER   :: MAXPOINTS = 10000
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
      REAL*4,            ALLOCATABLE :: PTAU_START(:)
      REAL*4,            ALLOCATABLE :: PTAU_END(:)

      ! For specifying lat/lon/alt/err and ID type
      REAL*4,            ALLOCATABLE :: PLAT(:)               
      REAL*4,            ALLOCATABLE :: PLON(:)               
      REAL*4,            ALLOCATABLE :: PPRESS(:)
      REAL*4,            ALLOCATABLE :: PXCO2(:)             
      REAL*4,            ALLOCATABLE :: PERR(:)
      CHARACTER(LEN=5),  ALLOCATABLE :: PTYPE(:)              
      CHARACTER(LEN=14), ALLOCATABLE :: PID(:)

      ! For keeping track of n_samples and avg interval
      INTEGER,           ALLOCATABLE :: NSAMP(:)
      REAL*8,           ALLOCATABLE :: CUR_CONC(:)
      INTEGER,           ALLOCATABLE :: SAMP_STRATEGY(:)

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

      ! Input/output file names
      CHARACTER(LEN=255)             :: INFILENAME,  INF1, INF2,INF3
      CHARACTER(LEN=255)             :: OUTFILENAME,OUTF1,OUTF2,OUTF3
      CHARACTER(LEN=255), PUBLIC     :: OCEANFILENAME
      CHARACTER(LEN=255)             :: OUTF

      ! Logical unit numbers (formerly in file_mod.F)
      INTEGER                        :: IU_FILE
      INTEGER                        :: IU_CT

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
! !IROUTINE: setup_CT
!
! !DESCRIPTION: Subroutine SETUP\_CT reads information from the 
!  input file in order to initialize the CT diagnostic.  Also 
!  calls INIT\_CT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETUP_CT ( am_I_Root, CT_Input_1, CT_Input_2, 
     &    CT_Input_3,CT_Output_1,CT_Output_2,CT_Output_3 )
!
! !USES:
!
      USE FILE_MOD,   ONLY : FILE_EXISTS
      USE FILE_MOD,   ONLY : IOERROR
      USE FILE_MOD,   ONLY : IU_FILE
      USE TIME_MOD,   ONLY : EXPAND_DATE
      USE TIME_MOD,   ONLY : GET_NYMD
      USE TIME_MOD,   ONLY : GET_NHMS
      USE TIME_MOD,   ONLY : GET_TIME_AHEAD
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM
      USE CTFORMAT, ONLY :ncCTOutputCreate,ncCTOutputClose
      USE CTFORMAT, ONLY : ctOutputNCID
!      USE NCCOOARDSFORMAT, ONLY :readGOSATlength
      USE NCCOOARDSFORMAT, ONLY :READGOSAT_ID,READGOSAT_REAL8 
      USE GIGC_In_Carbontracker_Mod
      USE CT_INPUT_READ_MOD

      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_read
      USE m_netcdf_io_write
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?

! !REMARKS:
!  For SMVGEAR or KPP simulations, the call to SETUP_CT is made from 
!  routine "chemdr.f", after the "chem.dat" file is read.  This is necessary 
!  since we have to reference the SMVGEAR rxn rate and species numbers.
!                                                                             .
!  For offline simulations, the call to SETUP_CT can be made
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
      LOGICAL            :: IS_OPEN
      INTEGER            :: I,  IP,      N,   TEMP, LENGTH
      INTEGER            :: RN, COUNTER, IOS, NYMD, NHMS
      INTEGER            :: NCID
      INTEGER            :: DAYAHEAD(2), DAYBEHIND(2)
      INTEGER            :: TEST

      CHARACTER(LEN=7)   :: NAMES
      CHARACTER(LEN=20)  :: LINE
      CHARACTER(LEN=10)  :: TYPE
      
      !-- Not exactly sure where to put this
      TYPE(Carbontracker_In)  :: CT_Input_1,CT_Input_2,CT_Input_3   ! CT input objects
      TYPE(Carbontracker_Out) :: CT_Output_1,CT_Output_2,CT_Output_3   ! CT input objects

      !=================================================================
      ! SETUP_CT begins here!
      !=================================================================

      ! Assume that there is flight data for today
      DO_CT = .TRUE.

      ! Find a free file LUN
      IU_FILE = findFreeLun()

      ! Get date & time
      NYMD  = GET_NYMD()
      NHMS  = GET_NHMS()

      !---------------
      !-- NEED TO FIND RIGHT NYMD/NHMS to put in find the "next day" INF/OUTF filenames
      ! Use Julian dates to automatically increase year/month
      DAYBEHIND = GET_TIME_AHEAD(-1440 )
      DAYAHEAD = GET_TIME_AHEAD(1440 )
      !CALL YMD_EXTRACT( YYYYMMDD1, YYYY1, MM1, DD1 )
      !JD  = JULDAY( YYYY, MM, DBLE(DD) )
      !JD1 = JD + 1d0
      !CALL CALDATE( JD1, YYYYMMDD1, HHMMSS )
      !CALL YMD_EXTRACT( YYYYMMDD1, YYYY1, MM1, DD1 )
      !------------

      ! Copy file names to local variables
      INF1   = INFILENAME
      INF2  = INFILENAME
      INF3  = INFILENAME
      OUTF1 = OUTFILENAME
      OUTF2  = OUTFILENAME
      OUTF3 = OUTFILENAME      

      ! Replace any date & time tokens in the file names
      CALL EXPAND_DATE( INF1,  DAYBEHIND(1), DAYBEHIND(2) )
      CALL EXPAND_DATE( INF2,  NYMD, NHMS )
      CALL EXPAND_DATE( INF3,  DAYAHEAD(1), DAYAHEAD(2) )
      CALL EXPAND_DATE( OUTF2, NYMD, NHMS )
      CALL EXPAND_DATE( OUTF1, DAYBEHIND(1), DAYBEHIND(2) )
      CALL EXPAND_DATE( OUTF3, DAYAHEAD(1), DAYAHEAD(2) )

      IF( FIRST ) THEN

       PRINT *,'INF1:',INF1
       PRINT *,'INF2:',INF2
       PRINT *,'INF3:',INF3

       IF ((.not. FILE_EXISTS( INF1 )).and.(.not.FILE_EXISTS(INF2)) 
     &            .and. (.not. FILE_EXISTS(INF3))) THEN
          DO_CT = .FALSE.
          RETURN
       END IF
       IF(FILE_EXISTS( INF1 ) ) THEN
          CALL READ_CT_NC( trim(INF1) , CT_Input_1 )
          CALL READ_POINTS( am_I_Root, CT_Input_1 )
       END IF
       IF(FILE_EXISTS( INF2 ) ) THEN
          CALL READ_CT_NC( trim(INF2) , CT_Input_2 )
          CALL READ_POINTS( am_I_Root, CT_Input_2 )
       END IF
       IF(FILE_EXISTS( INF3 ) ) THEN
          CALL READ_CT_NC( trim(INF3) , CT_Input_3 )
          !PRINT *,'testtest'
          flush(6)
          CALL READ_POINTS( am_I_Root, CT_Input_3 )
       END IF
      END IF
      !ELSE

      !PRINT *,'value of FIRST:',FIRST
      IF ( .not. FIRST ) THEN
       !-- WRITE OUTPUT OBJ1 TO FILE and THEN MOVE OUTPUT2 to OUTPUT1
       !PRINT *, 'length:',SIZE(CT_Input_1%VALUE)
       !PRINT *, 'IDS:',CT_Input_1%ID(1:100,1:100)
       TEST = CT_Input_1%NPOINTS
       CALL WRITE_OUTPUT(OUTF1, CT_Input_1, CT_Output_1, TEST )

       ! Copy Input2 to Input1 and read in new data in Input2
       !    (obspack) into state objects
       CALL COPY_INPUTS(am_I_Root, CT_Input_2, CT_Input_1)
       CALL COPY_INPUTS(am_I_Root, CT_Input_3, CT_Input_2)
       CALL READ_CT_NC( trim(INF3) , CT_Input_3 )
       CALL READ_POINTS( am_I_Root, CT_Input_3 )

       !-- MOVE OUTPUT2 to OUTPUT1
       CALL COPY_OUTPUTS(am_I_Root, CT_Output_2, CT_Output_1)
       CALL COPY_OUTPUTS(am_I_Root, CT_Output_3, CT_Output_2)
       CALL REINIT_OUTPUT(am_I_Root, CT_Output_3)
      ELSE
           FIRST = .FALSE.
      ENDIF


       !PRINT *,'min time (ct_input1):',MINVAL(CT_Input_1%TIME)
       !PRINT *,'max time (ct_input1):',MAXVAL(CT_Input_1%TIME)
       !PRINT *,'min time (ct_input2):',MINVAL(CT_Input_2%TIME)
       !PRINT *,'max time (ct_input2):',MAXVAL(CT_Input_2%TIME)
       !PRINT *,'min time (ct_input3):',MINVAL(CT_Input_3%TIME)
       !PRINT *,'max time (ct_input3):',MAXVAL(CT_Input_3%TIME)

       ! Echo info
       ! Echo info
       IF ( am_I_Root ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a)' ) 'P L A N E   F L I G H T   D I A G N O S T I C'
       WRITE( 6, 100   ) TRIM( INF2 )
  100  FORMAT( /, 'SETUP_CT: Reading ',a )
       WRITE( 6, '(a)' )  
       ENDIF

       ! Compute # of species and # of points & allocate arrays
       CALL INIT_CT( am_I_Root )
 
       ! Return if there are no flight track points for today
       !IF ( NPOINTS == 0 ) THEN
       !  IF ( am_I_Root ) THEN
       !    WRITE( 6, '(a)' ) 'No CT NOAA input data found for today!'
       !  ENDIF
       !  DO_CT = .FALSE.
       !  RETURN
       !ENDIF
      

      ! Copy Input2 to Input1 and read in new data in Input2
      !    (obspack) into state objects
      !CALL COPY_INPUTS(am_I_Root, CT_Input_2, CT_Input_1)
      !CALL READ_CT_NC( trim(INF2) , CT_Input_2 )


      ! Read variables to be output -- sort into PVAR array by type
      ! Leave hook for this in future maybe ... aes
      !CALL READ_VARIABLES
      PVAR(1) = 100001  
      PNAME(1) = 'TRA_001   '

      !CALL Ncop_Rd (NCID, trim(INF))

      !CALL NcRd (SAMP_STRATEGY, NCID , 'sampling_strategy',
      !&  (/ 1 /), (/ NPOINTS /))

      !PRINT *,'SAMP_STRA VECTOR:',SAMP_STRATEGY

      !CALL Nccl (NCID)

      ! Read information about each point (date/time/lon/lat/alt)
      !CALL READ_POINTS( am_I_Root, CT_Input_1 )

      ! Set the pointer to the first record 
      PPOINT = 1

      !=================================================================
      ! Find the species # for all components of RO2 (SMVGEAR only)
      !=================================================================
      CALL RO2_SETUP

      !=================================================================
      ! Find the species # for all components of NOY (SMVGEAR only)
      !=================================================================
      CALL NOY_SETUP

      !=================================================================
      ! Find the species # for all components of AN (SMVGEAR only)
      !=================================================================
      CALL AN_SETUP
    
      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE SETUP_CT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setup_CT
!
! !DESCRIPTION: Subroutine SETUP\_CT reads information from the
!  input file in order to initialize the CT diagnostic.  Also
!  calls INIT\_CT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CT_FINISH_OUTPUT ( am_I_Root,CT_Input_1,
     &    CT_Input_2,CT_Output_1,CT_Output_2 )
!
! !USES:
!
      USE FILE_MOD,   ONLY : FILE_EXISTS
      USE FILE_MOD,   ONLY : IOERROR
      !USE FILE_MOD,   ONLY : IU_FILE
      USE TIME_MOD,   ONLY : EXPAND_DATE
      USE TIME_MOD,   ONLY : GET_NYMD
      USE TIME_MOD,   ONLY : GET_NHMS
      USE TIME_MOD,   ONLY : GET_TIME_AHEAD
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM
      USE CTFORMAT, ONLY :ncCTOutputCreate,ncCTOutputClose
      USE CTFORMAT, ONLY : ctOutputNCID
      !USE NCCOOARDSFORMAT, ONLY :readGOSATlength
      !USE NCCOOARDSFORMAT, ONLY :READGOSAT_ID,READGOSAT_REAL8
      USE GIGC_In_Carbontracker_Mod
      USE CT_INPUT_READ_MOD

      !USE m_netcdf_io_open
      !USE m_netcdf_io_close
      !USE m_netcdf_io_read
      !USE m_netcdf_io_write
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?

! !REMARKS:
!  For SMVGEAR or KPP simulations, the call to SETUP_CT is made from
!  routine "chemdr.f", after the "chem.dat" file is read.  This is necessary
!  since we have to reference the SMVGEAR rxn rate and species numbers.
!                                                                             .
!  For offline simulations, the call to SETUP_CT can be made
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
      LOGICAL            :: IS_OPEN
      INTEGER            :: I,  IP,      N,   TEMP, LENGTH
      INTEGER            :: RN, COUNTER, IOS, NYMD, NHMS
      INTEGER            :: NCID
      INTEGER            :: DAYAHEAD(2), DAYBEHIND(2)
      INTEGER            :: TEST

      !-- Not exactly sure where to put this
      TYPE(Carbontracker_In)  :: CT_Input_1,CT_Input_2   ! CT input objects
      TYPE(Carbontracker_Out) :: CT_Output_1,CT_Output_2   ! CT input objects

      !=================================================================
      ! SETUP_CT begins here!
      !=================================================================

      ! Get date & time
      NYMD  = GET_NYMD()
      NHMS  = GET_NHMS()

      DAYBEHIND = GET_TIME_AHEAD(-1440 )
      DAYAHEAD = GET_TIME_AHEAD(1440 )

      ! Copy file names to local variables
      INF1  = INFILENAME
      INF2  = INFILENAME
      OUTF1  = OUTFILENAME
      OUTF2 = OUTFILENAME

      ! Replace any date & time tokens in the file names
      CALL EXPAND_DATE( INF2,  NYMD, NHMS )
      CALL EXPAND_DATE( INF1,  DAYBEHIND(1), DAYBEHIND(2) )
      CALL EXPAND_DATE( OUTF2, NYMD, NHMS )
      CALL EXPAND_DATE( OUTF1, DAYBEHIND(1), DAYBEHIND(2) )

       PRINT *,'INF1:',INF1
       PRINT *,'INF2:',INF2

       TEST = CT_Input_1%NPOINTS
       CALL WRITE_OUTPUT(OUTF1, CT_Input_1, CT_Output_1, TEST )

       TEST = CT_Input_2%NPOINTS
       CALL WRITE_OUTPUT(OUTF2, CT_Input_2, CT_Output_2, TEST )

      END SUBROUTINE CT_FINISH_OUTPUT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setup_CT
!
! !DESCRIPTION: Subroutine SETUP\_CT reads information from the
!  input file in order to initialize the CT diagnostic.  Also
!  calls INIT\_CT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETUP_CT_OCFL_OUT ( )
!
! !USES:
!
      USE FILE_MOD,   ONLY : FILE_EXISTS
      USE FILE_MOD,   ONLY : IOERROR
      !USE FILE_MOD,   ONLY : IU_FILE
      USE TIME_MOD,   ONLY : EXPAND_DATE
      USE TIME_MOD,   ONLY : GET_NYMD
      USE TIME_MOD,   ONLY : GET_NHMS
      USE TIME_MOD,   ONLY : GET_TIME_AHEAD
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM
      USE CTFORMAT, ONLY :ncCTOutputCreate,ncCTOutputClose
      USE CTFORMAT, ONLY : ctOutputNCID
      !USE NCCOOARDSFORMAT, ONLY :readGOSATlength
      !USE NCCOOARDSFORMAT, ONLY :READGOSAT_ID,READGOSAT_REAL8
      USE GIGC_In_Carbontracker_Mod
      USE CT_INPUT_READ_MOD
      USE CMN_SIZE_MOD

      ! References to F90 modules
      USE netcdf
      USE TRACER_MOD, ONLY : N_TRACERS
      USE GIGC_In_Carbontracker_Mod
      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_read
      USE m_netcdf_io_write
      USE m_netcdf_io_get_dimlen
      USE m_do_err_out
      USE m_netcdf_io_create
      USE m_netcdf_io_define
!
! !INPUT PARAMETERS:
!
      !LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?

! !REMARKS:
!  For SMVGEAR or KPP simulations, the call to SETUP_CT is made from
!  routine "chemdr.f", after the "chem.dat" file is read.  This is necessary
!  since we have to reference the SMVGEAR rxn rate and species numbers.
!                                                                             .
!  For offline simulations, the call to SETUP_CT can be made
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
      ! Local Variables
      integer                               :: NCID, DIM_LEN
      integer,  dimension(:), allocatable   :: SAMPLE_STRATEGY

      ! For netCDF file I/O
      INTEGER             :: idTime,    idLat,   idLon
      INTEGER             :: fId
      INTEGER             ::  vId,     omode,   i
      INTEGER             :: var1(3)


      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: IS_OPEN
      INTEGER            :: IP,      N,   TEMP, LENGTH
      INTEGER            :: RN, COUNTER, IOS, NYMD, NHMS
      INTEGER            :: DAYAHEAD(2), DAYBEHIND(2)
      INTEGER            :: TEST

      !=================================================================
      ! SETUP_CT begins here!
      !=================================================================

      ! Get date & time
      NYMD  = GET_NYMD()
      NHMS  = GET_NHMS()

      !DAYBEHIND = GET_TIME_AHEAD(-1440 )
      !DAYAHEAD = GET_TIME_AHEAD(1440 )

      ! Copy file names to local variables
      OUTF = OCEANFILENAME

      ! Replace any date & time tokens in the file names
      CALL EXPAND_DATE( OUTF,  NYMD, NHMS )

       PRINT *,'OCEAN OUTF:',OUTF

      !=========================================================================
      ! Create the netCDF file
      !=========================================================================

      ! Echo info
      WRITE( 6, '(a)' ) '=== Begin netCDF file creation test ==='
      CALL NcCr_Wr( fId, TRIM(OUTF) )

      ! Turn filling off
      !CALL NcSetFill( fId, NF_NOFILL, omode )

      !=========================================================================
      ! Define the dimensions
      !=========================================================================

      ! Time dimension
      WRITE( 6, '(a)' ) 'Writing obs  (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'time',   8,  idTime )
      !CALL NcDef_Unlim_Dimension( fId, 'obs',idObs)

      ! Latitude dimension
      CALL NcDef_Dimension( fId, 'lat',  180 , idLat )

      ! Longitude dimension
      CALL NcDef_Dimension( fId, 'lon', 360,  idLon )

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
      vId     = -1

      ! Define obspack_id variable
      !vId  = vId + 1
      var1 = (/ idLon, idLat, idTime /)
      CALL NcDef_Variable( fId, 'ocean_flux', NF90_DOUBLE, 3, var1, vId)
      !CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'ObsPack ID' )
      !CALL NcDef_Var_Attributes( fId, vId, 'units',     'hPa'      )

      !=========================================================================
      ! %%% END OF DEFINITION SECTION %%%
      ! %%% NOW WRITE DATA TO FILE    %%%
      !=========================================================================
      CALL NcEnd_def( fId )

      !=========================================================================
      ! Close the netCDF file
      !=========================================================================
      CALL NcCl( fId )

      ! Echo info
      WRITE( 6, '(a)' ) '===End oceanflux netCDF file creation test ==='

      END SUBROUTINE SETUP_CT_OCFL_OUT
!EOC
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
      SUBROUTINE READ_VARIABLES( am_I_Root )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,   ONLY : IU_FILE
      USE FILE_MOD,   ONLY : IOERROR
      USE TRACER_MOD, ONLY : N_TRACERS
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD     ! Size parameters
      USE COMODE_LOOP_MOD  ! NAMEGAS, NSPEC
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
      LOGICAL             :: IS_FULLCHEM
      INTEGER             :: M, N, NUM, R, IK, IOS
      CHARACTER(LEN=255)  :: LINE

      !=================================================================
      ! READ_VARIABLES begins here!
      !=================================================================

      ! Reset NCS to NCSURBAN for safety's sake (dbm, bmy, 10/16/06)
      NCS = NCSURBAN

      ! Test if this is a fullchem run
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM()

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
               IF ( NUM < 0 .or. NUM > N_TRACERS ) THEN                   
                  IF ( am_I_Root ) THEN
                    WRITE( 6, 100 ) TRIM( LINE )
 100                FORMAT( 'TRACER ', i4, ' is out of range!' )
                    WRITE( 6, '(a)' ) 'STOP in SETUP_CT!'
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
                     WRITE( 6, '(a)' ) 'STOP in CT!'
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
      SUBROUTINE READ_POINTS( am_I_Root, CT_Input )
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

      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_read
      USE GIGC_In_Carbontracker_Mod
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
      INTEGER                            :: YYYYMMDD_START, HHMMSS_START
      INTEGER                            :: YYYYMMDD_END, HHMMSS_END
      INTEGER                              :: QYY,QMM,QDD, QHH, QMN, SS
      INTEGER                           :: QYY_START,QMM_START,QDD_START
      INTEGER                           ::  QHH_START,QMN_START,SS_START
      INTEGER       :: QYY_END,QMM_END,QDD_END, QHH_END, QMN_END, SS_END
      INTEGER                              :: NCID
      INTEGER, dimension(:),allocatable    :: TIM, TIM_START, TIM_END
      REAL*8                     :: JULNDAY, JULNDAY_START, JULNDAY_END
      REAL*4, dimension(:), allocatable    :: LAT, LON, XCO2, ERR
      INTEGER*8, dimension(:), allocatable :: ID
      !CHARACTER(LEN=10) :: TYPE
      CHARACTER(LEN=14) :: SOUNDING_ID_CHAR

      TYPE(Carbontracker_In)  :: CT_Input

      !-- Int pointers
      REAL, POINTER :: P_TAU_START(:)
      REAL, POINTER :: P_TAU_END(:)
      INTEGER, POINTER :: P_NPOINTS
      !INTEGER, POINTER ::  P_TIME(:)
      !INTEGER, POINTER ::  P_SAMP_STRATEGY(:)

      !=================================================================
      ! READ_POINTS begins here!
      !=================================================================

      P_NPOINTS => CT_Input%NPOINTS
      P_TAU_START => CT_Input%TAU_START
      P_TAU_END   => CT_Input%TAU_END
      !P_TIME =>   CT_Input%TIME
      !P_SAMP_STRATEGY =>  CT_Input%SAMPLING_STRATEGY

      !=================================================================
      ! Read plane track points -- plane, lat/lon/alt, date/time
      ! We have previously computed NPOINTS in INIT_CT
      !=================================================================

      ALLOCATE( LAT( P_NPOINTS) )
      ALLOCATE( LON( P_NPOINTS) )
      ALLOCATE( XCO2( P_NPOINTS ))
      ALLOCATE( ERR( P_NPOINTS ))
      ALLOCATE( ID( P_NPOINTS ))

      ALLOCATE( TIM( P_NPOINTS ))
      ALLOCATE( TIM_START( P_NPOINTS ))
      ALLOCATE( TIM_END( P_NPOINTS ))

      !print *,'in readpoints'
      !print *,'INF:',trim(INF)

      !CALL READGOSAT_REAL4(trim(INF),'latitude','',LAT,NPOINTS)
      !CALL READGOSAT_REAL4(trim(INF),'longitude','',LON,
      !& NPOINTS)

      !print *,'LAT:',LAT
      !CALL READGOSAT_ID( trim (INF), NPOINTS, ID)
      !CALL READGOSAT_REAL4( trim (INF), 'xco2', '',
      !&  XCO2, NPOINTS)
      !CALL READGOSAT_REAL4( trim (INF), 'xco2_uncertainty', '',
      !&  ERR, NPOINTS)
     
      !CALL Ncop_Rd (NCID, trim(INF))

      !CALL NcRd (TIM, NCID,'time',(/ 1 /),(/ NPOINTS /))

      DO N = 1, P_NPOINTS
       
      ! PRINT *,'SAMP_STRATEGY(N):',CT_Input%SAMPLING_STRATEGY(N)
       SELECT CASE (CT_Input%SAMPLING_STRATEGY(N))
        CASE (1)
           !PRINT *,'hereherehere'
           TIM_START(N) = CT_Input%TIME(N) - 3600*2
           TIM_END(N)   = TIM_START(N) + 3600*2*2
        CASE (2)
           !PRINT *,'hereherehere'
           TIM_START(N) = CT_Input%TIME(N) - 1800
           TIM_END(N)   = TIM_START(N) + 1800*2
        CASE (3)
           TIM_START(N) = CT_Input%TIME(N) - 2700
           TIM_END(N)   = TIM_START(N) + 2700*2
        !CASE DEFAULT
        !   Others = Others + 1
       END SELECT
      !PRINT *,'TIM(N):',CT_Input%TIME(N)
      !PRINT *,'TIM_START(N):',TIM_START(N)

      ENDDO

      !CALL Nccl (NCID)
 

      !print *,'here2'

      !PLAT = LAT
      !PLON = LON
      
      !PXCO2 = XCO2
      !PERR = ERR

      !print *,'PID:',ID
 
       DO N = 1, P_NPOINTS

         !==============================================================
         ! Read date and time coordinates -- also do error checks
         !==============================================================

         !PRINT *,'TRY THIS:'
         CALL  CALDATE( 2440586.5_8, YYYYMMDD, HHMMSS )
         !         PRINT *,'YYYYMMDD:HHMMSS',YYYYMMDD,':',HHMMSS

         JULNDAY = REAL(CT_Input%TIME(N),8)/86400. + 2440587.5
         JULNDAY_START = REAL(TIM_START(N),8)/86400. + 2440587.5
         JULNDAY_END = REAL(TIM_END(N),8)/86400. + 2440587.5

         IF(N .eq. 1) THEN
        !PRINT *,'FIRST JULDAY:', REAL(CT_Input%TIME(N),8)/86400._8 
        !&   + 2440587.5_8
        !PRINT *,'FIRST JUL_START:',REAL(P_TAU_START(N),8)/86400._8 
        !&   + 2440587.5_8
        !PRINT *,'FIRST JUL_END:', REAL(P_TAU_END(N),8)/86400._8
        !&    + 2440587.5_8
         END IF

         !PRINT *,'JULNDAY:',JULNDAY,' TIM:',CT_Input%TIME(N)
         CALL CALDATE( JULNDAY, YYYYMMDD, HHMMSS )
         !PRINT *,'HHMMSS:',HHMMSS
         CALL CALDATE( JULNDAY_START, YYYYMMDD_START, HHMMSS_START )
         CALL CALDATE( JULNDAY_END, YYYYMMDD_END, HHMMSS_END )

         !PRINT *,'YYYYMMDD:HHMMSS',YYYYMMDD,':',HHMMSS

         CALL YMD_EXTRACT( YYYYMMDD, QYY, QMM, QDD ) 
         CALL YMD_EXTRACT( HHMMSS, QHH, QMN, SS )

         CALL YMD_EXTRACT(YYYYMMDD_START,QYY_START,QMM_START,QDD_START)
         CALL YMD_EXTRACT( HHMMSS_START,QHH_START,QMN_START,SS_START)

         CALL YMD_EXTRACT( YYYYMMDD_END,QYY_END,QMM_END,QDD_END)
         CALL YMD_EXTRACT( HHMMSS_END,QHH_END,QMN_END,SS_END)

        !print *,'POINT number:',N
        !print *,'PID(N):',PID(N)
 
        ! Error check MONTH
         IF ( QMM < 1 .or. QMM > 12 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 105   ) QMM
 105           FORMAT( 'ERROR: MONTH out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (CT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Error check DAY
         IF ( QDD < 1 .or. QDD > 31 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 110   ) QDD
 110           FORMAT( 'ERROR: DAY out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (CT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP 
         ENDIF

         ! Error check HOUR
         IF ( QHH < 0 .or. QHH > 23 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 115   ) QHH
 115           FORMAT( 'ERROR: HOUR out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (CT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP  
         ENDIF

         ! Error check MINUTES
         IF ( QMN < 0 .or. QMN > 59 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 120   ) QMN
 120           FORMAT( 'ERROR: MINUTES out of range: ', f8.3 )
               WRITE( 6, '(a)' ) 'STOP in READ_POINTS (CT_mod.f)'
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Store type in the global PTYPE array
         !PTYPE(N) = TYPE

         ! Store YYYYMMDD in the global PDATE array
         !PDATE(N) = ( QYY * 10000 ) + ( QMM * 100 ) + QDD

         ! Store HHMMSS in the global PTIME array
         ! (actaully we read in just HHMM, assume seconds = 00)
         !PTIME(N) = ( QHH * 100 ) + QMN

         ! Store TAU (hours since 1 Jan 1985) in the global PTAU array
         !PTAU(N)  = GET_TAU0( QMM, QDD, QYY, QHH, QMN, 0 )

         ! Store TAU (hours since 1 Jan 1985) in the global PTAU array
         !PTAU_START(N)  = GET_TAU0( QMM_START, QDD_START, QYY_START, 
         !&                            QHH_START, QMN_START, 0 )

         ! Store TAU (hours since 1 Jan 1985) in the global PTAU array
         !PTAU_END(N)  = GET_TAU0( QMM_END, QDD_END, QYY_END, 
         !&                            QHH_END, QMN_END, 0 )

         P_TAU_START(N) = GET_TAU0( QMM_START, QDD_START, QYY_START,
     &                            QHH_START, QMN_START, 0 )

         P_TAU_END(N) = GET_TAU0( QMM_END, QDD_END, QYY_END,
     &                            QHH_END, QMN_END, 0 )

         !==============================================================
         ! Read lon/lat/alt coordinates -- also do error checks
         !==============================================================

         ! Put LONGITUDE in the range [-180...180]
         !IF ( PLON(N) > 180.0 ) PLON(N) = PLON(N) - 360e0

         ! Error check LONGITUDE
         !IF ( PLON(N) < -180 .OR. PLON(N) > 180 ) THEN 
         !   IF ( am_I_Root ) THEN
         !      WRITE( 6, 125   ) PLON(N)
         ! 125           FORMAT( 'ERROR: Longitude out of range: ', f8.3 )
         !      WRITE( 6, '(a)' ) 'STOP in READ_POINTS (CT_mod.f)'
         !   ENDIF
         !   CALL GEOS_CHEM_STOP
         !ENDIF

         ! Error check LATITUDE
         !IF ( PLAT(N) < -90.0 .OR. PLAT(N) > 90.0 ) THEN 
         !   IF ( am_I_Root ) THEN
         !      WRITE( 6, 130   ) PLAT(N)
         ! 130           FORMAT( 'ERROR: Latitude out of range: ', f8.3 )
         !      WRITE( 6, '(a)' ) 'STOP in READ_POINTS (CT_mod.f)'
         !   ENDIF
         !   CALL GEOS_CHEM_STOP
         !ENDIF
        
       ENDDO

     

      NULLIFY(P_TAU_START)
      NULLIFY(P_TAU_END)
      !NULLIFY(P_TIME)
      !NULLIFY(P_SAMP_STRATEGY)

      DEALLOCATE( LAT )
      DEALLOCATE( LON )
      DEALLOCATE( XCO2 )
      DEALLOCATE( ERR )
      DEALLOCATE( ID )

      DEALLOCATE( TIM )
      DEALLOCATE( TIM_START )
      DEALLOCATE( TIM_END )

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      IF ( am_I_Root )  WRITE( 6, 135 ) P_NPOINTS
 135  FORMAT( /, 'Number of flight track points : ', i6 )

      NULLIFY(P_NPOINTS)

      END SUBROUTINE READ_POINTS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE RO2_SETUP
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD   ! Size parameters
      USE COMODE_LOOP_MOD       ! NSPEC, NAMEGAS, NCS
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
      IF ( .not. ITS_A_FULLCHEM_SIM() ) RETURN
      
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
         WRITE( 6, '(a)' ) 'STOP in RO2_SETUP (CT_mod.f)'
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE NOY_SETUP
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD   ! Size parameters
      USE COMODE_LOOP_MOD       ! NSPEC, NAMEGAS, NCS
! 
! !REVISION HISTORY: 
!  01 Jun 2009 - F. Paulot   - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
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
      IF ( .not. ITS_A_FULLCHEM_SIM() ) RETURN
      
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
         WRITE( 6, '(a)' ) 'STOP in NOY_SETUP (CT_mod.f)'
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE AN_SETUP
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD   ! Size parameters
      USE COMODE_LOOP_MOD       ! NSPEC, NAMEGAS, NCS
! 
! !REVISION HISTORY:
!  04 Jan 2010 - F. Paulot   - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
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
      IF ( .not. ITS_A_FULLCHEM_SIM() ) RETURN
      
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
         WRITE( 6, '(a)' ) 'STOP in AN_SETUP (CT_mod.f)'
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
! !IROUTINE: CT
!
! !DESCRIPTION: Subroutine CT saves concentrations to disk at 
!  locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CT(am_I_Root,Input_Opt,State_Met,State_Chm,
     &               CT_Input,CT_Output,RC)
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
      USE GIGC_In_Carbontracker_Mod, ONLY : Carbontracker_Out
      USE GIGC_In_Carbontracker_Mod, ONLY : Carbontracker_In
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP  
      USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TAU,        GET_TS_DIAG
      USE TRACER_MOD,   ONLY : TCVV,N_TRACERS
      USE OCEAN_MERCURY_MOD,  ONLY : Fp, Fg !eds 10/27/11
      USE PBL_MIX_MOD,   ONLY : GET_FRAC_UNDER_PBLTOP

      !USE CMN_SIZE_MOD, ONLY : LLPAR
      !USE CMN_FJ_MOD,   ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ
      USE JV_CMN_MOD         ! ODAER, QAA, QAA_AOD
      USE COMODE_LOOP_MOD    ! CSPEC, etc.
      USE NCCOOARDSFORMAT,  ONLY  : WRITE_VARS_TO_FILE_NC
      USE NCCOOARDSFORMAT,  ONLY  : ncCTOutputClose
      USE CTFORMAT       ,  ONLY  : WRITE_REALVAR_TO_FILE_CT
      USE CTFORMAT       ,  ONLY  : ctOutputNCID
      USE SMV_PHYSCONST_MOD, ONLY : g0, Rd
      USE m_netcdf_io_write
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
      TYPE(Carbontracker_In), INTENT(IN) :: CT_Input
      TYPE(Carbontracker_Out), INTENT(INOUT) :: CT_Output
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
!  (7 ) Now return if DO_CT = .FALSE. (bmy, 3/24/05)
!  (8 ) Now compute column AOD's and AOD's below plane.  Also now scale
!        AOD's to 400nm. (bmy, 10/24/05)
!  (9 ) Bug fix: exit if PTAU(M) == PTAUE, so that we write out on the next !
!        CT timestep (cdh, bmy, 12/12/06)
!  (10) Change CT output time step. (ccc, 8/27/09)
!  (11) Now scale AOD's to jv_spec_aod.dat wavelength. (clh, 5/14/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: PCHEM
      INTEGER             :: I, IP, IRHN, J, L, JLOOP, M, N,R,RH,V, TOP
      INTEGER             :: LPLANE
      INTEGER             :: VERTIND, TRACERIND
      REAL*8              :: TK, PTAUS, PTAUE, CONSEXP, VPRESH2O, SAODnm
      REAL*8              :: VARI(NPVAR)
      REAL*8              :: SUMPRESS
      REAL*8              :: VARIPBL(NPVAR)
      REAL*8              :: VARIFREETROP(NPVAR)
      REAL*8, ALLOCATABLE :: PARTVARI(:)
      REAL*8              :: OUT_CONC(NPVAR)
      REAL*8              :: CHEMSTEP_LEN
      REAL                :: Z2(47)
      REAL                :: PEDGE(48)

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

      ! Pointers to MetState objects
      REAL*8, POINTER :: MET_SURFACE_HEIGHT(:,:)
      REAL*8, POINTER :: MET_U(:,:,:)
      REAL*8, POINTER :: MET_V(:,:,:)
      REAL*8, POINTER :: MET_PBLH(:,:)
      REAL*8, POINTER :: MET_Q(:,:,:)
      REAL*8, POINTER :: MET_PRESSURE(:,:,:)
      REAL*8, POINTER :: MET_TEMPERATURE(:,:,:)

      !  Same for flask input
      REAL, POINTER    :: P_TAU_START(:)
      REAL, POINTER    :: P_TAU_END(:)
      REAL, POINTER :: P_ALTITUDE(:)
      INTEGER, POINTER :: P_NPOINTS

      ! Same for flask output
      REAL*8, POINTER :: P_FLASK(:,:)
      REAL*8, POINTER :: P_SURFACE_HEIGHT(:)
      REAL*8, POINTER :: P_U(:)
      REAL*8, POINTER :: P_V(:)
      REAL*8, POINTER :: P_BLH(:)
      REAL*8, POINTER :: P_Q(:)
      REAL*8, POINTER :: P_PRESSURE(:)
      REAL*8, POINTER :: P_TEMPERATURE(:)

      !-- Int pointers
      INTEGER, POINTER :: P_NSAMPLES(:)
      INTEGER, POINTER :: P_AVERAGING_TIME(:)
      INTEGER, POINTER :: P_REGION_INDICES(:,:)


      !-- Char pointers
      !CHARACTER, POINTER :: P_ID(:,:)
      !CHARACTER, POINTER :: P_TRACER_NAMES(:,:)
      !CHARACTER, POINTER :: P_REGION_NAME(:,:)

      !=================================================================
      ! CT begins here!
      !=================================================================

      !print *,'AD12:',AD12(40,65,:)

      ! Return if there is no flighttrack data for today
      !PRINT *,'DO_CT:',DO_CT
      IF ( .not. DO_CT ) RETURN


      P_TAU_START => CT_Input%TAU_START
      P_TAU_END => CT_Input%TAU_END
      P_NPOINTS => CT_Input%NPOINTS
      P_ALTITUDE => CT_Input%ALTITUDE

      ! Loop over all the locations that have not yet been found
      !DO M = PPOINT, NPOINTS
       DO M = 1, P_NPOINTS
         !print *,'ppoint:',ppoint
         !print *,'npoints:',npoints
         ! Starting & end times of chemistry interval
         PTAUE = GET_TAU()
         PTAUS = PTAUE - ( GET_TS_DIAG() / 60d0 )

         CHEMSTEP_LEN = PTAUE - PTAUS
         !SAMPINT_CT = 

         ! Initialize VARI to missing value for this point
         DO V = 1, NPVAR
            VARI(V) = MISSING
         ENDDO

         !print *,'M:',M
         !print *,'ptaue:',PTAUE
         !print *,'ptau_start(m):',P_TAU_START(M)
         !print *,'ptau(m):',PTAU(M)
         !print *,'ptau_end(m):',P_TAU_END(M)
         !print *,'ptaus:',PTAUS

         !==============================================================
         ! We haven't found the first plane point yet...
         !==============================================================
         !IF ( PTAU(M) < PTAUS ) THEN

            ! Write all missing values to disk for point #M
            !CALL WRITE_VARS_TO_FILE( M, VARI )      

            ! Increment pointer
            !PPOINT = PPOINT + 1

         !==============================================================
         ! We have already found all of the plane points...
         !==============================================================
         !ELSE IF ( PTAU(M) >= PTAUE ) THEN

            ! Exit this loop and the subroutine
            !EXIT

         !==============================================================
         ! We have found a plane point at the proper time & location!
         !==============================================================
         !ELSE
          
         IF (P_TAU_START(M) <= PTAUS .and. P_TAU_END(M) >= PTAUE ) THEN
            !PRINT *,'BINGO!!!!!!!!!'
            !PRINT *,'WRITING ',M,'th POINT'
            !PRINT *,'P_TAU_START(M):',P_TAU_START(M)
            !PRINT *,'PTAUS:',PTAUS
            !PRINT *,'PTAUE:',PTAUE
            ! PRINT *,'P_TAU_END(M):',P_TAU_END(M)

            ! Print the flight track point number
            !-- Commenting out for now, AESchuh
            !WRITE( 6, 100 ) PTYPE(M), PDATE(M), PTIME(M) 
 100        FORMAT( '     - CT: Archived ',a5,1x,i8.8,1x,i4.4 )

            ! Return grid box indices for the chemistry region
            ! NOTE: PCHEM and JLOOP are only defined for SMVGEAR runs!
            CALL TEST_VALID( M, PCHEM, JLOOP, CT_Input, I, J, L )

            !PRINT *,'M:',M,'J:',J,'I:',I,'L:',L

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
                     ALLOCATE(PARTVARI(LLPAR))

                     ! Remove offset from PVAR
                     N = PVAR(V) - 100000

                     ! Initialize GEOS-Chem tracer array [kg] from
                     ! Chemistry State object (mpayer, 12/6/12)
                     STT               => State_Chm%Tracers

                     !-- Points to met fields
                     MET_U             =>      State_Met%U
                     MET_V             =>      State_Met%V
                     MET_Q             =>      State_Met%SPHU
                     MET_PBLH           =>      State_Met%PBLH
                     MET_PRESSURE      =>      State_Met%PMID
                     MET_TEMPERATURE   =>      State_Met%T

                     !-- Pointers to output
                     P_FLASK           => CT_Output%FLASK
                     P_SURFACE_HEIGHT  => CT_Output%SURFACE_HEIGHT
                     P_NSAMPLES        => CT_Output%NSAMPLES
                     P_REGION_INDICES  => CT_Output%REGION_INDICES
                     P_U               => CT_Output%U
                     P_V               => CT_Output%V
                     P_BLH             => CT_Output%BLH
                     P_Q               => CT_Output%Q
                     P_PRESSURE        => CT_Output%PRESSURE
                     P_TEMPERATURE     => CT_Output%TEMPERATURE
                     P_AVERAGING_TIME  => CT_Output%AVERAGING_TIME


                     !print *,'N:',N
                     ! Convert from [kg] --> [v/v]
                     !VARI(V) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
                     !print *,'old:',STT(I,J,L,N)*TCVV(N) / AD(I,J,L)
                     !-- Temp kernel weighting for CT code
!                     DO VERTIND = 1, LLPAR       
!                       IF (VERTIND .eq. LLPAR) THEN
!      !                     PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
!      !&                             TCVV(N) / State_Met%AD(I,J,VERTIND))
!      !&                           * (GET_PEDGE(I,J,VERTIND) - 0)
!                          PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
!     &                             TCVV(N) / State_Met%AD(I,J,VERTIND))
!     &                           * KERN(VERTIND)
!                       ELSE 
!      !                     PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) * 
!      !&                             TCVV(N) / State_Met%AD(I,J,VERTIND))
!      !&                           * (GET_PEDGE(I,J,VERTIND)-
!      !&                          GET_PEDGE(I,J,VERTIND+1))
!                          PARTVARI(VERTIND) = (STT(I,J,VERTIND,N) *
!     &                             TCVV(N) / State_Met%AD(I,J,VERTIND))
!     &                           * KERN(VERTIND)
!                       ENDIF
!                     ENDDO

                     !PRINT *,'CUR_CONC:',CUR_CONC(M)
                     !PRINT *,'NSAMP:',NSAMP(M)

                     !PRINT *,'NEW CONC:',10**6 *
      !&                            STT(I,J,L,N)*TCVV(N)/
      !&                          State_Met%AD(I,J,L) 

      !               OUT_CONC = ( CUR_CONC(M)*NSAMP(M) + 10**6 * 
      !&                            STT(I,J,L,N)*TCVV(N)/
      !&                          State_Met%AD(I,J,L) ) / ( NSAMP(M) +1 )

                     !PRINT *,'TOTAL CONC:',OUT_CONC

                     !PRINT *,'SHAPES:',SHAPE(FLASK),' ',SHAPE(STT)
                     !CUR_CONC(M:M) = OUT_CONC

                     !NSAMP(M) = NSAMP(M) + 1
                     !PRINT *,'right before'
                     !PRINT *,'HEREHEREHERE:',
                     !P_FLASK(M,1) = OUT_CONC(1)

                  !PRINT *,'J:',J,'I:',I,'L:',L
                  !PRINT *,'U:',MET_U(J,I,L)                      

                  !-- Determine geopot height from sigma pressure levels

                  !PRINT *,'pmid first:',State_Met%PEDGE(1,1,:)

      ! Get L corresponding to PRESS(IND)
      
      DO VERTIND = 1, LLPAR+1
           PEDGE(VERTIND) = 100.*REAL(GET_PEDGE(I,J,VERTIND) ,4)
      ENDDO

      !PRINT *,'PEDGE:',PEDGE

                  !CALL CCM2_Z2_sigma( State_Met, I, J, 47, 
                  !&               Z2,REAL(g0,4),
                  !&               REAL(Rd,4), 0.61_4,
                  !&               PMID(1:47))

               CALL GEOPOT_HYPSO( State_Met, I, J, 47,
     &                        Z2,REAL(g0,4),
     &                       REAL(Rd,4), 0.61_4,
     &                        PEDGE(1:48))

               !NEED TO INTERPOLATE 'L' 

              TOP = 1
              DO VERTIND = 1,LLPAR
                   IF( Z2(VERTIND)  < P_ALTITUDE(M) ) THEN 
                      TOP = VERTIND
                   END IF
              END DO
               !PRINT *,'obs height:',P_ALTITUDE(M)
               !PRINT *,'setting L to:',TOP
               L = TOP
               
                  !--  Averaging in Met variables
           !PRINT *,'N_TRACERS:',N_TRACERS         
           DO TRACERIND = 1,N_TRACERS
             !PRINT *,'TRACERIND:',TRACERIND
             P_FLASK(M,TRACERIND) = (P_FLASK(M,TRACERIND)*P_NSAMPLES(M) 
     &             + 10**6 *STT(I,J,L,TRACERIND)*TCVV(TRACERIND)/
     &               State_Met%AD(I,J,L) ) / ( P_NSAMPLES(M) +1 )
             !PRINT *,'FLASKFLASK:',P_FLASK(M,TRACERIND)
           ENDDO
                     !PRINT  *,'done w/ tracers'
                     P_U(M)   = (MET_U(I,J,L) + P_U(M)*P_NSAMPLES(M))
     &                               /(P_NSAMPLES(M)+1) 
    
                     P_V(M)   = (MET_V(I,J,L) + P_V(M)*P_NSAMPLES(M))
     &                               /(P_NSAMPLES(M)+1)

                     P_BLH(M) = (MET_PBLH(I,J) + P_BLH(M)*P_NSAMPLES(M))
     &                               /(P_NSAMPLES(M)+1)

                     P_Q(M) = (MET_Q(I,J,L) + P_Q(M)*P_NSAMPLES(M))
     &                               /(P_NSAMPLES(M)+1)

                     P_PRESSURE(M) = (MET_PRESSURE(I,J,L) + 
     &                    P_PRESSURE(M)*P_NSAMPLES(M))/(P_NSAMPLES(M)+1)

                     P_TEMPERATURE(M) = (MET_TEMPERATURE(I,J,L) + 
     &                 P_TEMPERATURE(M)*P_NSAMPLES(M))/(P_NSAMPLES(M)+1)

                     !PRINT *,'setting indices to:',J,' ',I,' ',L

                     P_REGION_INDICES(M,1:3) = (/ J,I,L /)

                     P_AVERAGING_TIME(M) = P_AVERAGING_TIME(M) + 1800

                !    PRINT *,'Incrementing P_NSAMPLES by 1, new val:'
       !&               ,P_NSAMPLES(M),' for M:',M
                     P_NSAMPLES(M) = P_NSAMPLES(M) + 1

                     

      !      CALL NcWr ((/ J,I,L /),ctOutputNCID,'region_indices',
      !&           (/ M,1 /), (/ 1,3 /) )

                     ! Free pointer
                     NULLIFY( STT )
                     NULLIFY(P_FLASK)
                     NULLIFY(P_NSAMPLES)
                     NULLIFY(P_AVERAGING_TIME)
                     NULLIFY(P_REGION_INDICES)
                     NULLIFY(P_SURFACE_HEIGHT)
                     NULLIFY(P_U)
                     NULLIFY(P_V)
                     NULLIFY(P_BLH)
                     NULLIFY(P_Q)
                     NULLIFY(P_PRESSURE)
                     NULLIFY(P_TEMPERATURE)

                     NULLIFY(MET_SURFACE_HEIGHT)
                     NULLIFY( MET_U )
                     NULLIFY( MET_V )
                     NULLIFY( MET_PBLH )
                     NULLIFY( MET_Q )
                     NULLIFY( MET_PRESSURE )
                     NULLIFY( MET_TEMPERATURE )

                     DEALLOCATE(PARTVARI)

                  !--------------------------
                  ! Otherwise it's an error!
                  !--------------------------
                  CASE DEFAULT
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     WRITE( 6, '(a)' ) 'CT: Bad variable #!' 
                     WRITE( 6, '(a)' ) 'STOP in CT!'
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     CALL GEOS_CHEM_STOP

               END SELECT
            ENDDO

            ! Write data for the Mth plane point out to disk

            !print *,'M:',M,' OUT_CONC:',OUT_CONC
            flush(6)

            !CALL WRITE_REALVAR_TO_FILE_CT( VARI(1), M, 'flask' )
      !      CALL NcWr ( (/ OUT_CONC,0,0,0,0 /), ctOutputNCID, 
      !&                  'flask', (/ M,1 /), (/ 1,5 /) )

             

      !      CALL NcWr ((/ J,I,L /),ctOutputNCID,'region_indices',
      !&           (/ M,1 /), (/ 1,3 /) )

            !PRINT *,'CTOUT:',CT_Output_1%flask(5,1) != OUT_CONC(1)      
            !STOP
            !CALL WRITE_VARS_TO_FILE_NC( M, VARI(1), VARIPBL(1), 
            !&              VARIFREETROP(1) )

            ! Increment the record pointer
            !PPOINT = PPOINT + 1

            !IF (M .eq. NPOINTS) THEN 
              !call ncCTOutputClose
            !END IF
           
         ENDIF
      ENDDO

           NULLIFY(P_TAU_START)
           NULLIFY(P_TAU_END)
           NULLIFY(P_NPOINTS)
 
      !      IF ( am_I_Root ) THEN
      !         DO N = 1, LLPAR
      ! print*,'fraction of PBL:',N,':',GET_FRAC_UNDER_PBLTOP( 32, 65, N )
      !         ENDDO
      !      ENDIF

      END SUBROUTINE CT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE TEST_VALID( IND, PCHEM, JLOOP, CT_Input, I, J, L ) 
!
! !USES:
!
      USE COMODE_MOD,     ONLY : JLOP
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE GRID_MOD,       ONLY : GET_XOFFSET
      USE GRID_MOD,       ONLY : GET_YOFFSET
      USE GIGC_In_Carbontracker_Mod, ONLY: Carbontracker_In
      USE CMN_SIZE_MOD                ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: IND     ! # of the flight track point

      TYPE(Carbontracker_In)  :: CT_Input

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
      !I = INT( ( PLON(IND) + 180d0 - (I0 * DISIZE) ) / DISIZE + 1.5d0 )
      I = INT( ( CT_Input%LONGITUDE(IND) + 180d0 - (I0 * DISIZE) ) / 
     &          DISIZE + 1.5d0 )

      ! Handle date line correctly (bmy, 4/23/04)
      IF ( I > IIPAR ) I = I - IIPAR

      ! Get J corresponding to PLAT(IND)
      !J = INT( ( PLAT(IND) +  90d0 - (J0 * DJSIZE) ) / DJSIZE + 1.5d0 )
       J = INT( ( CT_Input%LATITUDE(IND) +  90d0 - (J0 * DJSIZE) ) / 
     &         DJSIZE + 1.5d0 )

      !PRINT *,'ID:',CT_Input%ID(:,IND)
      !PRINT *,'M:',IND,'LATGRIDIND:',J,'LONGRIDIND:',I,'LAT:',
      !&   CT_Input%LATITUDE(IND),'LON:',CT_Input%LONGITUDE(IND)

      ! Get L corresponding to PRESS(IND)
      L = 1
      DO IL = 1, LLPAR
         IF ( GET_PEDGE(I,J,IL) <= PPRESS(IND) .AND..NOT. FOUND ) THEN
            L     = IL-1
            FOUND =.TRUE.
            EXIT
         ENDIF          
      ENDDO

      ! Error check: L must be 1 or higher
      IF ( L == 0 ) L = 1

      !=================================================================
      ! We only do full-chemistry in the troposphere
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
  
         ! JLOOP indicates if a box is in tropo (/=0) or not. 
         JLOOP = JLOP(I,J,L)
         PCHEM = .TRUE.

      ELSE

         ! For non-SMVGEAR runs, JLOOP has no meaning so we give the
         ! stratospheric value
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
!        SETUP_CT at the start of each day (bmy, 3/25/05)
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
      WRITE( IU_CT, 110, IOSTAT=IOS ) 
     &     IND, PTYPE(IND), INT( PDATE(IND) ), INT( PTIME(IND) ),
     &     PLAT(IND), PLON(IND), PPRESS(IND), PXCO2(IND), PERR(IND), 
     &        ( VARI(I), I=1,NPVAR )
      
      ! Format string
      ! FP replace 95 by 200
! 110  FORMAT( I5,   X, A5,   X, I8.8, X, I4.4, X, 
!     &        F7.2, X, F7.2, X, F7.2, X, 95(es11.3e3,x) )
 110  FORMAT(I5,X,A5,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.0,X,
     &       F8.3,X,F7.2,X,200(es15.7e3,x))


      PRINT *,'PXCO2:',PXCO2(IND)
      WRITE( *, 110, IOSTAT=IOS )
     &     IND, PTYPE(IND), INT( PDATE(IND) ), INT( PTIME(IND) ),
     &     PLAT(IND), PLON(IND), PPRESS(IND), PXCO2(IND), PERR(IND),
     &        ( VARI(I), I=1,NPVAR )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_CT,'write_vars_to_file:1')

      ! Flush the file to disk
      CALL FLUSH( IU_CT )

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
!  to pass reaction rates from the SMVGEAR solver for the CT 
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
!  (5 ) Bug fix: Replace with DO_CT since this variable is reset to either T 
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
      IF ( DO_CT ) THEN

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
! !IROUTINE: set_CT
!
! !DESCRIPTION: Subroutine SET\_CT is used to pass values read in 
!  from the GEOS-Chem input file to "CT\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT( CT, IN_FILE, OUT_FILE, OCEAN_FILE )
!
! !INPUT PARAMETERS: 
!
      LOGICAL,            INTENT(IN) :: CT         ! Turn on CT diag?
      CHARACTER(LEN=255), INTENT(IN) :: IN_FILE    ! Input file to read
      CHARACTER(LEN=255), INTENT(IN) :: OUT_FILE   ! Output file to write
      CHARACTER(LEN=255), INTENT(IN) :: OCEAN_FILE   ! Another output file to write
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Save arguments to "shadow" module variables
      DO_CT       = CT
      INFILENAME  = TRIM( IN_FILE  )
      OUTFILENAME = TRIM( OUT_FILE )
      OCEANFILENAME = TRIM(OCEAN_FILE)

      END SUBROUTINE SET_CT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_CT
!
! !DESCRIPTION: Subroutine INIT\_CT reads the input file to compute 
!  the number of variables and flight track points to print out.  Also 
!  allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CT( am_I_Root )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IOERROR
!      USE NCCOOARDSFORMAT,  ONLY : readGOSATlength
      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_get_dimlen
  

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
!  (3 ) INIT_CT is now called each day but the arrays are only
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
      INTEGER           :: N, AS, IOS, NCID
      CHARACTER(LEN=20) :: LINE

      !=================================================================
      ! INIT_CT begins here!
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

      !NPOINTS = 0

      !!CALL readGOSATlength( trim (INF) , NPOINTS)
      !CALL Ncop_Rd (NCID, trim(INF))
      !CALL Ncget_Dimlen (NCID, 'obs', NPOINTS)
      !CALL Nccl (NCID)

      ! If there are no flight-track points then just return
      !IF ( NPOINTS < 1 ) THEN
      !   DO_CT = .FALSE.
      !   RETURN
      !ENDIF

      ! Make sure NPOINTS is less than MAXPOINTS
      !IF ( NPOINTS > MAXPOINTS ) THEN
      !   WRITE( 6, '(a)') 'NPOINTS exceeds maximum allowed value!'
      !   WRITE( 6, '(a)') 'STOP in INIT_CT (CT_mod.f)'
      !   WRITE( 6, '(a)') REPEAT( '=', 79 )
      !   CALL GEOS_CHEM_STOP
      !ENDIF
         
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
         ALLOCATE( NSAMP( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'NSAMP' )

         ALLOCATE( CUR_CONC( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CUR_CONC' )

         ALLOCATE( SAMP_STRATEGY( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SAMP_STRATEGY' )

         ALLOCATE( PTYPE( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTYPE' )

         ALLOCATE( PDATE( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDATE' )

         ALLOCATE( PTIME( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTIME' )

         ALLOCATE( PTAU( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTAU' )
         
         ALLOCATE( PTAU_START( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTAU_START' )

         ALLOCATE( PTAU_END( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTAU_END' )

         ALLOCATE( PLAT( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLAT' )

         ALLOCATE( PLON( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLON' )

         ALLOCATE( PPRESS( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PPRESS' )

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
      NSAMP         = 0
      CUR_CONC      = 0
      SAMP_STRATEGY = 0
      PREAC         = 0
      PRRATE        = 0e0
      PVAR          = 0
      PNAME         = ''
      PTYPE         = ''
      PDATE         = 0e0
      PTIME         = 0e0
      PTAU          = 0e0
      PTAU_START    = 0e0
      PTAU_END      = 0e0
      PLAT          = 0e0
      PLON          = 0e0
      PPRESS        = 0e0
      PID           = '              '

      END SUBROUTINE INIT_CT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_CT
!
! !DESCRIPTION: Subroutine CLEANUP\_CT deallocates all allocatable 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CT
!
! !USES:
!
      USE NCCOOARDSFORMAT,  ONLY : ncCTOutputClose
!! 
! !REVISION HISTORY: 
!  01 Jul 2001 - M. Evans    - Initial version
!  (1 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !CALL ncCTOutputClose

      PRINT *,'NSAMP:',NSAMP
      IF ( ALLOCATED( NSAMP ) ) DEALLOCATE ( NSAMP )
      IF ( ALLOCATED( CUR_CONC ) ) DEALLOCATE ( CUR_CONC )
      IF ( ALLOCATED( SAMP_STRATEGY ) ) DEALLOCATE ( SAMP_STRATEGY )
      IF ( ALLOCATED( PVAR   ) ) DEALLOCATE( PVAR   )
      IF ( ALLOCATED( PREAC  ) ) DEALLOCATE( PREAC  )
      IF ( ALLOCATED( PNAME  ) ) DEALLOCATE( PNAME  )
      IF ( ALLOCATED( PRRATE ) ) DEALLOCATE( PRRATE )
      IF ( ALLOCATED( PTYPE  ) ) DEALLOCATE( PTYPE  ) 
      IF ( ALLOCATED( PDATE  ) ) DEALLOCATE( PDATE  ) 
      IF ( ALLOCATED( PTIME  ) ) DEALLOCATE( PTIME  ) 
      IF ( ALLOCATED( PTAU   ) ) DEALLOCATE( PTAU   ) 
      IF ( ALLOCATED( PTAU_START   ) ) DEALLOCATE( PTAU_START   )
      IF ( ALLOCATED( PTAU_END   ) ) DEALLOCATE( PTAU_END   )
      IF ( ALLOCATED( PLAT   ) ) DEALLOCATE( PLAT   )
      IF ( ALLOCATED( PLON   ) ) DEALLOCATE( PLON   )
      IF ( ALLOCATED( PPRESS ) ) DEALLOCATE( PPRESS )
      IF ( ALLOCATED( PERR   ) ) DEALLOCATE( PERR   )
      IF ( ALLOCATED( PXCO2  ) ) DEALLOCATE( PXCO2  )
      IF ( ALLOCATED( PID    ) ) DEALLOCATE( PID  )

      END SUBROUTINE CLEANUP_CT
!EOC
!>============================================================================<
      SUBROUTINE GEOPOT_HYPSO( State_Met, I, J, nlvl,
     &                           Z2, g, Rd, Qmf, pedge )

       USE GIGC_State_Met_Mod, ONLY : MetState

       IMPLICIT NONE

       TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

       INTEGER                              :: N

       INTEGER, INTENT(IN)                  :: nlvl,I,J
                                          ! Number of levels.

       REAL*8,POINTER                    :: T(:,:,:)
                                          ! Temperature, Kelvin.

       REAL*8, POINTER                   :: Q(:,:,:)
                                          ! Specific humidity,
                                          ! kg/kg.

       REAL*8, POINTER                                 :: PS(:,:)
                                          ! Surface pressure, Pa.

       REAL*8, POINTER                                 :: PHIS(:,:)
                                          ! Surface geopotential,
                                          ! m^2/s^2.

       REAL, DIMENSION(1:nlvl)              :: Z
       REAL, DIMENSION(1:nlvl), INTENT(OUT) :: Z2
                                          ! Geopotential height
                                          ! based on CCM2 hydro-
                                          ! static formulation.

       REAL, INTENT(IN)                     :: g
                                          ! Gravitational ac-
                                          ! celeration, m/s^2.

       REAL, INTENT(IN)                     :: Rd
                                          ! Dry gas constant,
                                          ! J / (kg K)

       REAL, INTENT(IN)                     :: Qmf
                                          ! Multiplicative factor
                                          ! of Q in formula for
                                          ! virtual temperature,
                                          ! i.e.:
                                          ! TV = T(1 + Qmf*Q)

       REAL, DIMENSION(1:nlvl+1), INTENT(IN)   :: pedge
       REAL, DIMENSION(1:nlvl+1)               :: pedge2

                 !  Targets to pointers
                 T    => State_Met%T
                 Q    => State_Met%SPHU
                 PS   => State_Met%PS1
                 PHIS => State_Met%PHIS

      pedge2 = pedge((nlvl+1):1:-1)
      !pmid(nlvl+1) = PS(I,J)

      DO N = 1, NLVL
       !PRINT *,'Rd:',Rd,'T:',T(I,J,NLVL-N+1),'Q:',Q(I,J,NLVL-N+1),'g:',g
       
         !Z2(N) = Rd*(T(I,J,NLVL-N+1)*(1+Qmf*Q(I,J,NLVL-N+1)))*1/g*
              Z(N) = PHIS(I,J)/g +Rd*(T(I,J,NLVL-N+1))*1/g*
     &              LOG(pedge2(N+1)/pedge2(N))
       !PRINT *,'partZ2:',Z(N)
      ENDDO

      Z2 = 0

      DO N=1, NLVL
        Z2(N) = SUM(Z(NLVL-N+1:NLVL)) 
      ENDDO
      !PRINT *,'Z2:',Z2

      END SUBROUTINE GEOPOT_HYPSO


!EOC
!>============================================================================<
!           FILE: SUB_CCM2_Z2_sigma.f90
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 9 September 1998
!  LAST MODIFIED: Mon Sep 14 13:41:39 MDT 1998
!
!    DESCRIPTION: Computes Z2 on sigma levels. Z2 (a CCM Processor Code-Defined
!                 Derived Field) is the geopotential height based on the CCM2
!                 hydrostatic formulation. For a fixed time, latitude, and
!                 longitude, the geopotential height at model mid-levels may
!                 be written
!
!                              Z2 = PHIS/g + (R/g) * H * TV
!
!                 where PHIS is the surface geopotential, g the gravitational
!                 acceleration, R the dry gas constant, H the CCM2 hydrostatic
!                 matrix (which is a function of surface pressure and pressure
!                 at each model mid-level), and TV the virtual temperature.
!                 Subroutine CCM2_Z2_sigma contains ccm2_hydro_mat, a subrou-
!                 tine which computes H.
!
!                 The ordering of the vertical profiles of T, Q, and sigma
!                 must be TOP to BOTTOM. The vertical profile of Z2 on exit is
!                 also TOP to BOTTOM.
!
!>============================================================================<
      SUBROUTINE CCM2_Z2_sigma( State_Met, I, J, nlvl, 
     &                           Z2, g, Rd, Qmf, pmid_in )

       USE GIGC_State_Met_Mod, ONLY : MetState

       IMPLICIT NONE

       TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object       

       INTEGER, INTENT(IN)                  :: nlvl,I,J
                                          ! Number of levels.

       REAL*8,POINTER                    :: T(:,:,:)
                                          ! Temperature, Kelvin.

       REAL*8, POINTER                   :: Q(:,:,:)
                                          ! Specific humidity,
                                          ! kg/kg.

       REAL*8, POINTER                                 :: PS(:,:)
                                          ! Surface pressure, Pa.

       REAL*8, POINTER                                 :: PHIS(:,:)
                                          ! Surface geopotential,
                                          ! m^2/s^2.

       REAL, DIMENSION(1:nlvl), INTENT(OUT) :: Z2
                                          ! Geopotential height
                                          ! based on CCM2 hydro-
                                          ! static formulation.

       REAL, INTENT(IN)                     :: g
                                          ! Gravitational ac-
                                          ! celeration, m/s^2.

       REAL, INTENT(IN)                     :: Rd
                                          ! Dry gas constant,
                                          ! J / (kg K)

       REAL, INTENT(IN)                     :: Qmf
                                          ! Multiplicative factor
                                          ! of Q in formula for
                                          ! virtual temperature,
                                          ! i.e.:
                                          ! TV = T(1 + Qmf*Q)

       REAL, DIMENSION(1:nlvl), INTENT(IN)   :: pmid_in
       REAL, DIMENSION(1:nlvl)               :: pmid
       !REAL, DIMENSION(1:nlvl), INTENT(IN)  :: sigma
       !                                   ! Sigma levels top to
       !                                   ! bottom

!>----------------------------------------------------------------------------<
! Local variables:

       REAL                                 :: Rdog
                                          ! Dry gas constant di-
                                          ! vided by gravitational
                                          ! acceleration.

       REAL, DIMENSION(:), ALLOCATABLE      :: TV
                                          ! Virtual temperature,
                                          ! K.

       REAL                                 :: ZS
                                          ! Geopotential height
                                          ! at surface, m.

       REAL, DIMENSION(:,:), ALLOCATABLE    :: H
                                          ! CCM2 hydrostatic
                                          ! matrix.
 
       !REAL, DIMENSION(:), ALLOCATABLE      :: pmid
                                          ! Mid-level pressures
                                          ! computed from sigma,
                                          ! Pa.

!>----------------------------------------------------------------------------
!  Targets to pointers
                 T    => State_Met%T
                 Q    => State_Met%SPHU
                 PS   => State_Met%PS1
                 PHIS => State_Met%PHIS

       !PRINT *,'in geopotential routine'
       pmid = pmid_in(nlvl:1:-1)
       !PRINT *,'pmid:',pmid
  
       
!>----------------------------------------------------------------------------<
! Assign dry gas constant divided by gravitational acceleration:

       Rdog = Rd / g

!>----------------------------------------------------------------------------<
! Allocate and compute virtual temperature and geopotential height of surface:

       ALLOCATE( TV(1:nlvl) )

       TV = T(I,J,:) * ( 1. + Qmf*Q(I,J,:) )

       ZS = PHIS(I,J) / g

!>----------------------------------------------------------------------------<
! Compute mid-level pressures from sigma and surface pressure, compute CCM2
! hydrostatic matrix H, and compute Z2. The subroutine ccm2_hydro_mat is called
! to compute H.

       !ALLOCATE( pmid(1:nlvl ) )

       ALLOCATE( H(1:nlvl,1:nlvl) )

       !pmid = sigma * PS

       CALL ccm2_hydro_mat(nlvl, REAL(PS(I,J),4), pmid, H)

       PRINT *,'H:',H(1:2,1:2)
       Z2(:) = ZS + Rdog * MATMUL( H, TV(:) )

       PRINT *,'Z2:'

!>----------------------------------------------------------------------------<
! Deallocate virtual temperature, surface geopotential height, mid-level
! pressures, and CCM2 hydrostatic matrix.

      DEALLOCATE( TV )

      !DEALLOCATE( pmid )

      DEALLOCATE( H )

!>----------------------------------------------------------------------------<

      CONTAINS

      SUBROUTINE ccm2_hydro_mat(ucK, pi, p, H)

!    DESCRIPTION: Denoting the number of midpoint model levels by K,
!                 SUBROUTINE ccm2_hydro_mat returns the [K x K] CCM2
!                 hydostatic matrix H defined by Equation (3.a.109)
!                 in NCAR Technical Note 382 (Hack et al., 1993). The
!                 notation in this code follows the notation found in
!                 the Technical Note.
!
!                 Let k be the row index of H, l the column index of
!                 H, p_k the midpoint pressure where k = 1,2...,K arranged
!                 from top to bottom, and pi the surface pressure.
!                 (The midpoint pressure values are derived, for example,
!                 from purely sigma or hybrid coordinate midpoint
!                 values.) Then, from Hack et al., p. 27, H(k,l) is
!                 defined as:
!
!                            /     0,                            l < k
!                           |
!                           |      .5 ln( p_(k+1)/p_k ),         l = k, k < K
!                           |
!                 H(k,l) = <       .5 ln( p_(l+1)/p_(l-1) ),     l > k, k < K
!                           |
!                           |      .5 ln( pi^2/(p_(K-1)p_K) ),   l = K, k < K
!                           |
!                            \     ln( pi/p_K )                  l = K, k = K
!
!
!
!                 A prime example of the utility of H is in the calculation
!                 of Z2 (a code-defined derived field) from the CCM processor.
!                 Z2 is the geopotential height based on the CCM2 hydrostatic
!                 formulation. For a given vertical profile of virtual tempera-
!                 ture, say TV_k, Z2 in Fortran 90 is given by the matrix
!                 equation
!
!                         Z2 = PHIS/g + (R/g) * H * TV
!
!                 where PHIS is the surface geopotential and R is the gas
!                 constant for dry air. In this case the column vector Z2
!                 represents the geopoential height arranged from top to
!                 bottom.
!
!      REFERENCE: Hack, J.J., B.A. Boville, B.P. Briegleb, J.T. Kiehl, P.J.
!                 Rasch, D.L. Williamson, 1993: Description of the NCAR
!                 Community Climate Model (CCM2). NCAR Technical Note
!                 NCAR/TN-382+STR, 108 pp.
!
!           NOTE: K and ucK are used interchangeably in comments.


      IMPLICIT NONE

      INTEGER, INTENT(IN)                :: ucK  ! Number of midpoint model levels
                                             ! (uc refers to upper case).

      REAL, INTENT(IN)                   :: pi   ! Surface pressure, Pa or mb.

      REAL, DIMENSION(1:ucK), INTENT(IN) :: p    ! Midpoint pressure values,
                                             ! ARRANGED FROM TOP TO BOTTOM.
                                             ! Units Pa or mb, but same as
                                             ! units used for pi. The midpoint
                                             ! pressure values are derived,
                                             ! for example, from purely sigma
                                             ! or hybrid coordinate midpoint
                                             ! values.

       REAL, DIMENSION(1:ucK,1:ucK), INTENT(OUT)  :: H
                                             ! CCM2 hydrostatic matrix defined
                                             ! by Equation (3.a.109) in Hack
                                             ! et al. (1993). H is indexed
                                             ! as H(k,l) where k and l are
                                             ! defined below.

! Local variables:

       INTEGER                           :: k    ! Lower case k representing the
                                            ! row index of H.

       INTEGER                           :: l    ! Lower case l representing the
                                            ! column index of H.

! Initialize all elements of H as 0.:

       H = 0.

! Compute all diagonal elements of H except H(K,K):

      DO k = 1, ucK - 1

        H(k,k) = .5 * LOG( p(k+1) / p(k) )

      END DO

! Compute all off-diagonal elements of H, except for last two rows, and last
! column:

      DO k = 1, ucK - 2

        DO l = k + 1, ucK - 1

             H(k,l) = .5 * LOG( p(l+1) / p(l-1) )

        END DO

      END DO

! Compute last column of H, except k = K:

      DO k = 1, ucK - 1

         H(k,ucK) = .5 * LOG( pi*pi / (p(ucK-1)*p(ucK)) )

      END DO

! Compute H(K,K):

      H(ucK,ucK) = LOG( pi / p(ucK) )


      END SUBROUTINE ccm2_hydro_mat

!>----------------------------------------------------------------------------<

      END SUBROUTINE CCM2_Z2_sigma

!>----------------------------------------------------------------------------<

      END MODULE CT_MOD
