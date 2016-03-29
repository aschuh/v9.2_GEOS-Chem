! $Id: i6_read_mod.f,v 1.2 2009/10/15 14:52:20 bmy Exp $
      MODULE CT_INPUT_READ_MOD
!
!******************************************************************************
!  Module I6_READ_MOD contains subroutines that unzip, open, and read
!  GEOS-CHEM I-6 (instantaneous 6-hr) met fields from disk.
!  (bmy, 6/23/03, 10/15/09)
!
!  Module Routines:
!  =========================================================================
!  (1 ) UNZIP_I6_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) OPEN_I6_FIELDS  : Opens met field files residing in the temp dir
!  (3 ) GET_I6_FIELDS_1 : Wrapper for routine READ_I6
!  (4 ) GET_I6_FIELDS_2 : Wrapper for routine READ_I6
!  (5 ) GET_N_I6        : Returns # of A-3 fields for each DAO data set
!  (6 ) READ_I6         : Reads A-3 fields from disk
!  (7 ) I6_CHECK        : Checks if we have found all of the fields
!
!  GEOS-CHEM modules referenced by i6_read_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f     : Module containing NaN and other error check routines
!  (6 ) logical_mod.f   : Module containing GEOS-CHEM logical switches
!  (7 ) file_mod.f      : Module containing file unit #'s and error checks
!  (8 ) time_mod.f      : Module containing routines for computing time & date
!  (9 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!  (10) tropopause_mod.f: Module containing routines for tropopause (only for
!                         variable one here)
!  (11) unix_cmds_mod.f : Module containing Unix commands for unzipping etc.
!
!  NOTES:
!  (1 ) Adapted from "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Now use TIMESTAMP_STRING for formatted date/time output (bmy, 10/28/03)
!  (3 ) Now reads either zipped or unzipped files (bmy, 12/11/03)
!  (4 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (5 ) Now references "directory_mod.f", "unix_cmds_mod.f", and
!        "logical_mod.f" (bmy, 7/20/04)
!  (6 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  (7 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (8 ) Now account for GEOS-4 coastal boxes in LWI properly (bmy, 8/10/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now make LWI REAL*8 for near-land formulation (ltm, bmy, 5/9/06)
!  (11) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (12) Now set negative SPHU to a very small positive # (bmy, 9/8/06)
!  (13) Now read TROPP files for GEOS-4, and check tropopause level
!       in case of a variable tropopause (phs, bmy, bdf, 9/14/06)
!  (14) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!  (15) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  03 Aug 2012 - R. Yantosca - Now make IU_I6 a private module variable
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE inquireMod, ONLY : findFreeLUN
     

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "i6_read_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: READ_SCALING_NC
      PUBLIC :: READ_CT_NC
      PUBLIC :: WRITE_OUTPUT

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Current and 'next' input files of CT obspack data
      

      ! Number of I6 fields in the file
      INTEGER :: N_I6_FIELDS

      ! File units for I-6 met fields
      INTEGER :: IU_I6
      INTEGER :: IU_TP

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_CT_NC( filname , In_Carbontracker )
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
      USE GIGC_In_Carbontracker_Mod
      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_read
      USE m_netcdf_io_get_dimlen
      USE m_do_err_out

      ! Arguments
      TYPE(Carbontracker_In), INTENT(INOUT) :: In_Carbontracker   ! Meteorology State object

      character* (*),         INTENT(IN)    :: filname

      ! Local Variables
      integer                               :: NCID, DIM_LEN
      integer,  dimension(:), allocatable   :: SAMPLE_STRATEGY

      CALL Ncop_Rd (NCID, TRIM(filname))

      CALL Ncget_Unlim_Dimlen (NCID, DIM_LEN)

      PRINT *,'DIM_LEN:',DIM_LEN
      PRINT *,'SIZE:'
      In_Carbontracker%NPOINTS = DIM_LEN
      PRINT *,'test1'
      CALL NcRd(In_Carbontracker%SAMPLING_STRATEGY(1:DIM_LEN),NCID,
     &  'sampling_strategy',
     &  (/ 1 /), (/ DIM_LEN /))

      CALL NcRd(In_Carbontracker%TIME(1:DIM_LEN),NCID,
     &  'time',
     &  (/ 1 /), (/ DIM_LEN /))

      CALL NcRd(In_Carbontracker%LATITUDE(1:DIM_LEN),NCID,
     &  'latitude',
     &  (/ 1 /), (/ DIM_LEN /))

      CALL NcRd(In_Carbontracker%LONGITUDE(1:DIM_LEN),NCID,
     &  'longitude',
     &  (/ 1 /), (/ DIM_LEN /))

      CALL NcRd(In_Carbontracker%ALTITUDE(1:DIM_LEN),NCID,
     &  'altitude',
     &  (/ 1 /), (/ DIM_LEN /))

      CALL NcRd(In_Carbontracker%VALUE(1:DIM_LEN),NCID,
     &  'value',
     &  (/ 1 /), (/ DIM_LEN /))

      CALL NcRd(In_Carbontracker%ID(1:100,1:DIM_LEN),NCID,
     &  'obspack_id',
     &  (/ 1, 1 /), (/ 100, DIM_LEN /))


      !PRINT *,'ID CHECK:',In_Carbontracker%ID

      CALL Nccl (NCID) 

      ! Return to calling program
      END SUBROUTINE READ_CT_NC

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE READ_SCALING_NC(filname,week , scaling_matrix )
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
      USE m_netcdf_io_open
      USE m_netcdf_io_close
      USE m_netcdf_io_read
      USE m_netcdf_io_get_dimlen
      USE m_do_err_out
      USE CMN_SIZE_MOD      ! Size parameters 
      USE ENKF_MOD ,ONLY    : ENUMBER

      ! INPUT
      character* (*),    INTENT(IN)    :: filname
      integer,           INTENT(IN)    :: week 

      ! OUTPUT
      REAL*4,            INTENT(OUT)   :: SCALING_MATRIX(360,180,1) 

      ! Local Variables
      integer                               :: NCID, DIM_LEN,IREG
      integer,  dimension(:), allocatable   :: SAMPLE_STRATEGY
      REAL, DIMENSION(240)                  :: ECO_SCALING_MEAN
      REAL, DIMENSION(240)                  :: ECO_SCALING_DEV
      REAL, DIMENSION(240)                  :: ECO_SCALING
      integer                               :: ECO_MAP(360,180)

      CALL Ncop_Rd (NCID, 
     &    '/discover/nobackup/aschuh/data/misc/iregions.nc')
      CALL NcRd(ECO_MAP,NCID,'regions',
     &  (/ 1,1 /), (/ 360,180 /))
      CALL Nccl (NCID)

      CALL Ncop_Rd (NCID, TRIM(filname))
      CALL Ncget_Unlim_Dimlen (NCID, DIM_LEN)
      PRINT *,'DIM_LEN:',DIM_LEN
      CALL NcRd(ECO_SCALING_MEAN,NCID,'scaling_factor',
     &  (/ 1,week,1 /), (/ 240,1,1 /))
      CALL NcRd(ECO_SCALING_DEV,NCID,'scaling_factor_deviations',
     &  (/ 1,ENUMBER,week,1 /), (/ 240,1,1,1 /))
      CALL Nccl (NCID)

      ECO_SCALING = ECO_SCALING_MEAN + ECO_SCALING_DEV

      do ireg = 1,240 !min(inbuf),max(inbuf)

       where (ECO_MAP .eq. ireg) 
     &        SCALING_MATRIX(:,:,1) = ECO_SCALING(ireg)

      end do

      !status=0      

      ! Return to calling program
      END SUBROUTINE READ_SCALING_NC

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_OUTPUT(fname,In_Carbontracker,Out_Carbontracker,
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

      ! Arguments
      TYPE(Carbontracker_In), INTENT(IN) :: In_Carbontracker   ! Meteorology State object
      TYPE(Carbontracker_Out), INTENT(IN) :: Out_Carbontracker   ! Meteorology State object
      character* (*),         INTENT(IN)    :: fname
      INTEGER,                INTENT(IN)    :: NPOINTS

      ! Local Variables
      integer                               :: NCID, DIM_LEN
      integer,  dimension(:), allocatable   :: SAMPLE_STRATEGY

      ! For netCDF file I/O
      INTEGER             :: idObs,    idTracer,   idTname, idChar100
      INTEGER             :: idRgNmLen,  idGridInd
      INTEGER             :: fId
      INTEGER             ::  vId,     omode,   i
      !INTEGER             :: ct1d(1),  ct2d(2), ct3d(3), ct4d(4)
      !INTEGER             :: st1d(1),  st2d(2), st3d(3), st4d(4)
      INTEGER             :: var1(2),var2(2),var3(1),var4(2) 
      INTEGER             ::  var5(2),var6(2)

      !-- Pointer to CT input object ID
      CHARACTER, POINTER :: P_ID(:,:)

      !=========================================================================
      ! Create the netCDF file
      !=========================================================================


      !PRINT *,'flask data:',Out_Carbontracker%flask(1:NPOINTS,1)
      !print *,'#avg samps:',Out_Carbontracker%averaging_time(1:NPOINTS)
      !print *,'strategy:', In_Carbontracker%Sampling_strategy(1:NPOINTS)
      !print *,'tau_start:', In_Carbontracker%TAU_START(1:NPOINTS)
      !print *,'tau_end:', In_Carbontracker%TAU_END(1:NPOINTS)
      ! Echo info
      WRITE( 6, '(a)' ) '=== Begin netCDF file creation test ==='
      PRINT *,'WRITING ',NPOINTS,' observations'
      CALL NcCr_Wr( fId, TRIM(fname) )

      ! Turn filling off
      !CALL NcSetFill( fId, NF_NOFILL, omode )

      !=========================================================================
      ! Define the dimensions
      !=========================================================================

      ! Time dimension
      WRITE( 6, '(a)' ) 'Writing obs  (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'obs',    NPOINTS,  idObs )
      !CALL NcDef_Unlim_Dimension( fId, 'obs',idObs)

      ! Altitude dimension
      WRITE( 6, '(a)' ) 'Writing tracer   (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'tracer',      N_TRACERS,  idTracer )

      ! Latitude dimension
      WRITE( 6, '(a)' ) 'Writing tracer_n_len (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'tracer_name_len',  8 , idTname )

      ! Longitude dimension
      WRITE( 6, '(a)' ) 'Writing Char100  (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'lon', 100,  idChar100 )

      ! Character dimension 1
      WRITE( 6, '(a)' ) 'Writing rnmlen (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'rnmlen', 10, idRgNmLen )

      ! Character dimension 1
      WRITE( 6, '(a)' ) 'Writing gridind (dim     ) to netCDF file'
      CALL NcDef_Dimension( fId, 'gridind', 3, idGridInd )

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

      ! Define obspack_id variable
      vId  = vId + 1
      var1 = (/ idChar100, idObs /)
      CALL NcDef_Variable( fId, 'obspack_id', NF90_CHAR, 2, var1, vId )
      CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'ObsPack ID' )
      !CALL NcDef_Var_Attributes( fId, vId, 'units',     'hPa'      )

      ! Define flask variable
      vId  = vId + 1
      var2 = (/ idObs, idTracer /)
      CALL NcDef_Variable( fId, 'flask', NF90_DOUBLE, 2, var2, vId )
      !CALL NcDef_Var_Attributes( fId, vId, 'long_name', 
      !&     'mole_fraction_of_trace_gas_in_air'     )
      !CALL NcDef_Var_Attributes( fId, vId, 'units',
      !&      'mol tracer (mol air)^-1')
      !CALL NcDef_Var_Attributes( fId, vId, '_FillValue',
      !&      -1.e+34)

      ! Define nsamples variable
      var3 = (/ idObs /)
      vId  = vId + 1
      CALL NcDef_Variable( fId, 'nsamples', NF90_INT, 1, var3, vId )
      CALL NcDef_Var_Attributes( fId, vId,  'comment', 
     &  'Number of discrete samples in flask average.'    )

      ! Define tracer names variable
      vId  = vId + 1
      var4 = (/ idTracer, idTname /)
      CALL NcDef_Variable( fId, 'tracer_names', NF90_CHAR, 2, var4, vId)

      ! Define averaging_time variable
      vId  = vId + 1
      !var4 = (/  idObs /)
      CALL NcDef_Variable( fId,'averaging_time',NF90_INT,1, var3, vId )
      CALL NcDef_Var_Attributes( fId, vId, 'units', 'seconds'        )
      CALL NcDef_Var_Attributes( fId, vId, 'comment', 
     &   'Amount of model time over which this sample is averaged.')

      ! Define description variable
      vId  = vId + 1
      !var2 = (/ idObs /)
      CALL NcDef_Variable(fId,'surface_height',NF90_DOUBLE,1,var3, vId )
      CALL NcDef_Var_Attributes( fId, vId, 'comment',
     &    'Height of the GEOS-Chem surface in the flask gridbox'  )
      CALL NcDef_Var_Attributes( fId, vId, 'units','meters'  )

      ! Define region names variable
      vId  = vId + 1
      var5 = (/ idObs, idRgNmLen /)
      CALL NcDef_Variable( fId, 'region_name', NF90_CHAR, 2, var5, vId )
      CALL NcDef_Var_Attributes( fId, vId, 'units', 'unitless'        )
      CALL NcDef_Var_Attributes( fId, vId, 'comment',
     &   'Name of TM5 zoom region containing sample.?????')

      ! Define region indices variable
      vId  = vId + 1
      var6 = (/ idObs, idGridInd /)
      CALL NcDef_Variable( fId, 'region_indices', NF90_INT,2,var6, vId )
      CALL NcDef_Var_Attributes( fId, vId, 'comment',
     &   'Zonal, meridional, and level indices within the model region
     &      for grid cell containing sample.')

      !A bunch of met output
      CALL NcDef_Variable( fId, 'u', NF90_DOUBLE, 1, var3, vId )
      CALL NcDef_Variable( fId, 'v', NF90_DOUBLE, 1, var3, vId )
      CALL NcDef_Variable( fId, 'blh', NF90_DOUBLE, 1, var3, vId )
      CALL NcDef_Variable( fId, 'q', NF90_DOUBLE, 1, var3, vId )
      CALL NcDef_Variable( fId, 'pressure', NF90_DOUBLE, 1, var3, vId )
      CALL NcDef_Variable( fId, 'temperature', NF90_DOUBLE,1,var3, vId )

      !=========================================================================
      ! %%% END OF DEFINITION SECTION %%%
      ! %%% NOW WRITE DATA TO FILE    %%%
      !=========================================================================
      CALL NcEnd_def( fId )

      P_ID     =>   In_Carbontracker%ID

      PRINT *,'FIRST POINT:',P_ID(1:100,1)
      CALL NcWr(P_ID(1:100,1:NPOINTS) ,fId, 'obspack_id',
     &   (/ 1,1 /),(/ 100, NPOINTS /) )

      CALL NcWr(Out_Carbontracker%flask(1:NPOINTS,1:N_TRACERS),fId,
     &    'flask',(/ 1,1 /),(/ NPOINTS,N_TRACERS /) )

      CALL NcWr(Out_Carbontracker%nsamples(1:NPOINTS),fId,'nsamples',
     &   (/ 1 /),(/ NPOINTS /) )

      ! CALL NcWr(Out_Carbontracker%tracer_names(1:5,1:10),fId,
      !&  'tracer_names', (/ 1,1 /),(/ 5,10 /) )

      CALL NcWr(Out_Carbontracker%averaging_time(1:NPOINTS),fId,
     &  'averaging_time', (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%surface_height(1:NPOINTS),fId,
     &  'surface_height', (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%region_name(1:NPOINTS,1:10),fId,
     &  'region_name', (/ 1,1 /),(/ NPOINTS,10 /) )

      CALL NcWr(Out_Carbontracker%region_indices(1:NPOINTS,1:3),fId,
     &  'region_indices', (/ 1,1 /),(/ NPOINTS,3 /) )

      CALL NcWr(Out_Carbontracker%u(1:NPOINTS),fId,'u',
     &   (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%v(1:NPOINTS),fId,'v',
     &   (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%blh(1:NPOINTS),fId,'blh',
     &   (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%q(1:NPOINTS),fId,'q',
     &   (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%pressure(1:NPOINTS),fId,'pressure',
     &   (/ 1 /),(/ NPOINTS /) )

      CALL NcWr(Out_Carbontracker%temperature(1:NPOINTS),fId,
     &   'temperature', (/ 1 /),(/ NPOINTS /) )

      !=========================================================================
      ! Close the netCDF file
      !=========================================================================
      CALL NcCl( fId )
    
      NULLIFY(P_ID)
  
      ! Echo info
      WRITE( 6, '(a)' ) '=== End netCDF file creation test ==='

      ! Return to calling program
      END SUBROUTINE WRITE_OUTPUT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CT_INPUT_READ_MOD
