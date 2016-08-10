! TODO list
  !
  ! how does GC parallelize?  not sure how different PEs will handle
  ! different obs and/or tracers.  then fix the i_am_root
  ! conditionals.
  !
  ! Can we call setup_obspack and cleanup_obspack at 0Z each day, for
  ! multi-day runs?
  !
  ! BOP/EOP Protex headers
  !
  ! Read variables/tracers to output from input file?  (See
  ! planeflight_mod.F)
  !
  ! Write variable/tracer names to output file?  Compare with
  ! user_output_flask.F90
  !
  ! RC error code passing
  !
  ! Range checking on input variables (lon -180:180 for instance)
  ! 
  ! Vertical sampling location
  !
  ! Output surface height?
  !

  
  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Transport Model                  !
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: OBSPACK_MOD
  !
  ! !DESCRIPTION: Module OBSPACK\_MOD contains variables and routines
  !  which to sample a GEOS-Chem model simulation for in situ observations
  !  contained in an ObsPack file.
  !\\
  !\\
  ! !INTERFACE:
  !
  MODULE OBSPACK_MOD

    ! !USES:
    !
    USE GIGC_ErrCode_Mod
    USE Error_Mod,          ONLY : Debug_Msg

    IMPLICIT NONE

#    INCLUDE 'netcdf.inc'

    PRIVATE
    !
    ! !PUBLIC MEMBER FUNCTIONS:
    !
    PUBLIC  :: OBSPACK_SAMPLE
    PUBLIC  :: OBSPACK_INIT
    PUBLIC  :: OBSPACK_CLEANUP
    PUBLIC  :: OBSPACK_WRITE_OUTPUT
    PUBLIC  :: DO_OBSPACK
    PUBLIC  :: GEOPOT_HYPSO
    !
    ! !PRIVATE MEMBER FUNCTIONS:
    !
    PRIVATE  :: OBSPACK_READ_INPUT
    PRIVATE  :: OBSPACK_GET_GRID_INDICES
    
    !
    ! !REMARKS:
    !   4 Jun 2015 - A. Jacobson - Adapted from v10.1 planeflight_mod.f, following
    !                              similar work done in v9.2 by Andrew Schuh.
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !PRIVATE TYPES:
    !
    !=================================================================
    ! MODULE VARIABLES:
    !
    ! DO_OBSPACK  : Turn on the obspack diagnostic? (T/F)
    ! NOBS     : Number of flight track points 
    ! OBSPACK_INFILE  : Name of obspack input file
    ! OBSPACK_OUTFILE : Name of obspack output file 
    !=================================================================

    ! Logicals
    LOGICAL                        :: DO_OBSPACK
    LOGICAL                        :: prtDebug

    CHARACTER(LEN=*), PARAMETER	   :: mname = 'obspack_mod'

    ! For specifying variables to save for each ObsPack measurement
    INTEGER                        :: nobs
    INTEGER                        :: ntracers

    ! Input/output file names
    CHARACTER(LEN=255)             :: obspack_infile
    CHARACTER(LEN=255)             :: obspack_outfile

    CHARACTER(len=100), DIMENSION(:), ALLOCATABLE	:: obspack_id           ! unique sample identifier
    INTEGER, DIMENSION(:), ALLOCATABLE			:: obspack_nsamples     ! number of samples accumulated
    INTEGER, DIMENSION(:), ALLOCATABLE			:: obspack_strategy     ! sampling strategy
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_latitude     ! sample latitude
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_longitude    ! sample longitude
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_altitude     ! sample altitude (meters above sea level)
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_tau_start    ! sampling start time
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_tau_center   ! sampling center time
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_tau_end      ! sampling end time
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_accum_weight ! accumulated weight
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_u            ! meteorology u-wind
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_v            ! meteorology v-wind
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_blh          ! meteorology boundary layer height
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_q            ! meteorology water vapor mixing ratio
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_pressure     ! meteorology pressure
    REAL*8, DIMENSION(:), ALLOCATABLE			:: obspack_temperature  ! meteorology temperature
    REAL*8, DIMENSION(:,:), ALLOCATABLE		:: obspack_tracers      ! mixing ratio of tracers

    !=================================================================
    ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
    !=================================================================

  CONTAINS

    !EOC
    !------------------------------------------------------------------------------
    !                  GEOS-Chem Global Chemical Transport Model                  !
    !------------------------------------------------------------------------------
    !BOP
    !
    ! !IROUTINE: obspack_init
    !
    ! !DESCRIPTION: Subroutine OBSPACK\_INIT reads information from
    !  the input ObsPack file in order to initialize the obspack
    !  diagnostic.  It allocates memory in the obs array.  If this
    !  structure is already allocated, the routine presumes that the
    !  current samples need to be written out before the new ObsPack
    !  information is read.  In that case, OBSPACK\_WRITE\_OUTPUT is
    !  called.
    !\\
    !\\
    ! !INTERFACE:
    !
    SUBROUTINE OBSPACK_INIT( am_I_Root, Input_Opt, RC )
      !
      ! !USES:
      !
      USE FILE_MOD,		ONLY : FILE_EXISTS
      USE GIGC_Input_Opt_Mod,	ONLY : OptInput 
      USE TIME_MOD,		ONLY : EXPAND_DATE
      USE TIME_MOD,		ONLY : GET_NYMD
      USE TIME_MOD,		ONLY : GET_NHMS
      USE TIME_MOD,		ONLY : GET_TIME_AHEAD

      !
      ! !INPUT PARAMETERS:
      !
      LOGICAL,        INTENT(IN)  		:: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  		:: Input_Opt   ! Input Options object

      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)		:: RC          ! Success or failure?

      !EOP
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      CHARACTER(LEN=*), PARAMETER 		:: rname = mname//'/obspack_init'
      INTEGER, DIMENSION(2)			:: tomorrow
      INTEGER                                   :: NYMD, NHMS


      !=================================================================
      ! OBSPACK_INIT begins here
      !=================================================================

      RC            =  GIGC_SUCCESS


      ! Set a flag for debugging
      prtDebug = ( am_I_Root .and. Input_Opt%LPRT )

      ! Assume that there are ObsPack data for today
      DO_OBSPACK   = .TRUE.

      WRITE (6,'("obspack_init")') 
      FLUSH(6)
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### OBSPACK_INIT: starting' )
      ENDIF

      ! If obs array is already allocated, assume that we need to
      ! write out existing results before reading new input.  This
      ! could happen in a multi-day run with daily input files.

      !PRINT *,'obspackid allocated?',ALLOCATED(obspack_id)

      IF (ALLOCATED(obspack_id)) THEN
         CALL obspack_write_output(Am_I_Root,RC)
         CALL obspack_cleanup(RC)
      ENDIF

      nobs = 0

      ! Form input and output file names with daily timestamps.
      ! Example: flask_input.2011063000_2011070100.nc and
      ! flask_output.2011063000_2011070100.nc. The current
      ! implementation uses today's time for the first timestamp, and
      ! today's time plus 24 hours for the second.  We could also use
      ! the TAUb and TAUe times, if the input files have been
      ! processed accordingly.

      ! Get current date & time
      NYMD    = GET_NYMD()
      NHMS    = GET_NHMS()

      WRITE (obspack_infile,'(a)') TRIM(Input_Opt%OBSPACK_IFILE)
      WRITE (obspack_outfile,'(a)') TRIM(Input_Opt%OBSPACK_OFILE)
      ! replace YYYYMMDD with date and time
      CALL EXPAND_DATE( obspack_infile,  NYMD, NHMS )
      CALL EXPAND_DATE( obspack_outfile, NYMD, NHMS )

      !PRINT*, 'OBSPACK INPUT FILE',obspack_infile
      !FLUSH(6)
      !PRINT*, 'OBSPACK OUTPUT FILE',obspack_outfile
      !FLUSH(6)

      ! If we can't find a ObsPack file for today's date, return
      IF ( .NOT. FILE_EXISTS( obspack_infile ) ) THEN 
         DO_OBSPACK = .FALSE.
         PRINT*, 'INPUT FILE DOES NOT EXIST!'
         RETURN
      ENDIF

      ! Get number of tracers to sample
      ntracers = Input_Opt%N_TRACERS

      CALL OBSPACK_READ_INPUT(am_I_root, rc)
      RETURN

    END SUBROUTINE OBSPACK_INIT
    !EOC




    !------------------------------------------------------------------------------
    !                  GEOS-Chem Global Chemical Transport Model                  !
    !------------------------------------------------------------------------------
    !BOP
    !
    ! !IROUTINE: obspack_read_input
    !
    ! !DESCRIPTION: Subroutine OBSPACK\_READ\_INPUT allocates space
    !  for variables in the input file, and reads those data into the
    !  module arrays.
    !\\
    !\\
    ! !INTERFACE:
    !
    SUBROUTINE OBSPACK_READ_INPUT(am_I_root, RC )
      !
      ! !USES:
      !
      USE GIGC_Input_Opt_Mod, 		ONLY : OptInput 
      USE m_netcdf_io_open, 		ONLY : Ncop_Rd
      USE m_netcdf_io_get_dimlen, 	ONLY : Ncget_Dimlen
      USE m_netcdf_io_read
      USE m_netcdf_io_close, 		ONLY : Nccl
      USE BPCH2_MOD,			ONLY : GET_TAU0
      USE FILE_MOD,			ONLY : FILE_EXISTS
      !
      ! !INPUT PARAMETERS:
      !
      LOGICAL,        INTENT(IN)  		:: am_I_Root   ! Is this the root CPU?
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)		:: RC          ! Success or failure?


      ! 
      ! !REVISION HISTORY: 
      !  05 Jun 2015 - A. Jacobson - first version
      !EOP
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      CHARACTER(LEN=*), PARAMETER 		:: rname = mname//'/obspack_read_input'
      INTEGER 					:: ncid, iobs
      INTEGER, DIMENSION(:,:), ALLOCATABLE	:: central_time
      INTEGER, DIMENSION(1)			:: starts_1d, counts_1d
      INTEGER, DIMENSION(2)			:: starts_2d, counts_2d

      !=================================================================
      ! OBSPACK_READ_INPUT begins here
      !=================================================================

      RC            =  GIGC_SUCCESS

      ! Get from nc input file:
      !	lat, lon, alt, date_components, sampling_strategy, obspack_id
      !   restrict to obs between start and end times of simulation

      IF(am_I_Root) THEN

         IF (.NOT. FILE_EXISTS(obspack_infile)) THEN
            DO_OBSPACK = .FALSE.
            nobs = 0
            PRINT*, 'INPUT FILE DOES NOT EXIST',obspack_infile
            FLUSH(6)
            RETURN
         ENDIF

         CALL Ncop_Rd( ncid, obspack_infile)
         CALL Ncget_Dimlen( ncid, 'obs', nobs )
         WRITE (6,'("[obspack_read_input] ",i," obs in input file.")') nobs
         FLUSH(6)
         
         ! We read in all data available in the input file,
         ! but it is possible that there are observations that
         ! fall outside the TAUb-TAUe time period.

         ALLOCATE(obspack_id(nobs))
         ALLOCATE(obspack_nsamples(nobs))
         ALLOCATE(obspack_strategy(nobs))
         ALLOCATE(obspack_latitude(nobs))
         ALLOCATE(obspack_longitude(nobs))
         ALLOCATE(obspack_altitude(nobs))
         ALLOCATE(obspack_tau_start(nobs))
         ALLOCATE(obspack_tau_center(nobs))
         ALLOCATE(obspack_tau_end(nobs))
         ALLOCATE(obspack_u(nobs))
         ALLOCATE(obspack_v(nobs))
         ALLOCATE(obspack_blh(nobs))
         ALLOCATE(obspack_q(nobs))
         ALLOCATE(obspack_pressure(nobs))
         ALLOCATE(obspack_temperature(nobs))
         ALLOCATE(obspack_tracers(nobs,ntracers))

         obspack_nsamples=0
         obspack_u=0.d0
         obspack_v=0.d0
         obspack_blh=0.d0
         obspack_q=0.d0
         obspack_pressure=0.d0
         obspack_temperature=0.d0
         obspack_tracers=0.d0


         ! central_time is a local work array
         ALLOCATE(central_time(6,nobs))

         starts_1d = (/ 1 /)
         counts_1d = (/ nobs /)

         CALL NcRd(obspack_latitude,ncid,'latitude', starts_1d, counts_1d)
         CALL NcRd(obspack_longitude,ncid,'longitude', starts_1d, counts_1d)
         CALL NcRd(obspack_altitude,ncid,'altitude', starts_1d, counts_1d)
         CALL NcRd(obspack_strategy,ncid,'sampling_strategy', starts_1d, counts_1d)

         starts_2d = (/ 1, 1 /)
         counts_2d = (/ 100, nobs /)

         CALL NcRd(obspack_id,ncid,'obspack_id', starts_2d, counts_2d)


         starts_2d = (/ 1, 1 /)
         counts_2d = (/ 6, nobs /)

         CALL NcRd(central_time,ncid,'time_components', starts_2d, counts_2d)

         ! Close input netCDF file
         CALL NcCl(ncid)

         ! fill smapling window start and end times
         DO iobs = 1,nobs

            obspack_tau_center(iobs) = GET_TAU0( &
                 central_time(2,iobs), &
                 central_time(3,iobs), &
                 central_time(1,iobs), &
                 central_time(4,iobs), &
                 central_time(5,iobs), &
                 central_time(6,iobs))

            SELECT CASE (obspack_strategy(iobs))

            CASE (1) ! 4-hour window

               obspack_tau_start(iobs) = obspack_tau_center(iobs)-2.0 !*3600
               obspack_tau_end(iobs) = obspack_tau_center(iobs)+2.0 !*3600

            CASE (2) ! 1-hour window

               obspack_tau_start(iobs) = obspack_tau_center(iobs)-0.5 !*3600
               obspack_tau_end(iobs) = obspack_tau_center(iobs)+0.5 !*3600

            CASE (3) ! 90-minute window

               obspack_tau_start(iobs) = obspack_tau_center(iobs)-0.75 !*3600
               obspack_tau_end(iobs) = obspack_tau_center(iobs)+0.75  !*3600

            CASE default

               WRITE (6,'("[user_output_flask_init] Obs with obspack_id string ",a,":")') TRIM(ADJUSTL(obspack_id(iobs)))
               WRITE (6, '("  Unknown sampling strategy = ",i,".")') obspack_strategy(iobs)
               FLUSH(6)

               RC=1
               RETURN

            END SELECT

               !PRINT *,'obspack_tau_start:', obspack_tau_start(iobs)
               !PRINT *,'obspack_tau_end:', obspack_tau_end(iobs)

         ENDDO

      ENDIF ! i_am_root

      PRINT*, 'OBSPACK No. OBS',nobs
      FLUSH(6)

      IF(ALLOCATED(central_time)) THEN
         DEALLOCATE(central_time)
      ENDIF


    END SUBROUTINE OBSPACK_READ_INPUT
    !EOC




    !------------------------------------------------------------------------------
    !                  GEOS-Chem Global Chemical Transport Model                  !
    !------------------------------------------------------------------------------
    !BOP
    !
    ! !IROUTINE: obspack_cleanup
    !
    ! !DESCRIPTION: Subroutine OBSPACK\_CLEANUP deallocates all allocatable 
    !  module arrays.
    !\\
    !\\
    ! !INTERFACE:
    !
    SUBROUTINE OBSPACK_CLEANUP(RC)
      ! 
      ! !REVISION HISTORY: 
      !  05 Jun 2015 - A. Jacobson - first version
      !EOP
      !------------------------------------------------------------------------------

      CHARACTER(LEN=*), PARAMETER 		:: rname = mname//'/obspack_cleanup'

      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)		:: RC          ! Success or failure?

      RC            =  GIGC_SUCCESS

      !BOC
      IF(ALLOCATED(obspack_id)) THEN
         DEALLOCATE(obspack_id)
      ENDIF

      IF(ALLOCATED(obspack_nsamples)) THEN
         DEALLOCATE(obspack_nsamples)
      ENDIF

      IF(ALLOCATED(obspack_strategy)) THEN
         DEALLOCATE(obspack_strategy)
      ENDIF

      IF(ALLOCATED(obspack_latitude)) THEN
         DEALLOCATE(obspack_latitude)
      ENDIF

      IF(ALLOCATED(obspack_longitude)) THEN
         DEALLOCATE(obspack_longitude)
      ENDIF

      IF(ALLOCATED(obspack_altitude)) THEN
         DEALLOCATE(obspack_altitude)
      ENDIF

      IF(ALLOCATED(obspack_tau_start)) THEN
         DEALLOCATE(obspack_tau_start)
      ENDIF

      IF(ALLOCATED(obspack_tau_center)) THEN
         DEALLOCATE(obspack_tau_center)
      ENDIF

      IF(ALLOCATED(obspack_tau_end)) THEN
         DEALLOCATE(obspack_tau_end)
      ENDIF

      IF(ALLOCATED(obspack_accum_weight)) THEN
         DEALLOCATE(obspack_accum_weight)
      ENDIF

      IF(ALLOCATED(obspack_u)) THEN
         DEALLOCATE(obspack_u)
      ENDIF

      IF(ALLOCATED(obspack_v)) THEN
         DEALLOCATE(obspack_v)
      ENDIF

      IF(ALLOCATED(obspack_blh)) THEN
         DEALLOCATE(obspack_blh)
      ENDIF

      IF(ALLOCATED(obspack_q)) THEN
         DEALLOCATE(obspack_q)
      ENDIF

      IF(ALLOCATED(obspack_pressure)) THEN
         DEALLOCATE(obspack_pressure)
      ENDIF

      IF(ALLOCATED(obspack_temperature)) THEN
         DEALLOCATE(obspack_temperature)
      ENDIF

      IF(ALLOCATED(obspack_tracers)) THEN
         DEALLOCATE(obspack_tracers)
      ENDIF


    END SUBROUTINE OBSPACK_CLEANUP
    !EOC


    !------------------------------------------------------------------------------
    !                  GEOS-Chem Global Chemical Transport Model                  !
    !------------------------------------------------------------------------------
    !BOP
    !
    ! !IROUTINE: obspack_write_output
    !
    ! !DESCRIPTION: Subroutine OBSPACK\_WRITE\_OUTPUT computes window averages
    !  and writes data to output file
    !\\
    !\\
    ! !INTERFACE:
    !
    SUBROUTINE OBSPACK_WRITE_OUTPUT(Am_I_Root,RC)
      ! 
      ! !Revision history: 
      !  05 Jun 2015 - A. Jacobson - First version
      !EOP
      !------------------------------------------------------------------------------

      USE time_mod, ONLY : GET_NYMDb
      USE time_mod, ONLY : GET_NHMSb
      USE time_mod, ONLY : GET_NYMDe
      USE time_mod, ONLY : GET_NHMSe
      USE time_mod, ONLY : YMD_EXTRACT
      USE time_mod, ONLY : SYSTEM_TIMESTAMP
      USE m_netcdf_io_define
      USE m_netcdf_io_create
      USE m_netcdf_io_write
      USE m_netcdf_io_close     

      !BOC

      !
      ! !LOCAL VARIABLES
      !
      INTEGER					:: ncid
      INTEGER					:: omode
      INTEGER 					:: dimid_obs, dimid_tracer, dim_tnmlen, dimid_char100
      INTEGER 					:: vid
      INTEGER					:: iobs
      CHARACTER(len=255)			:: attstring
      INTEGER					:: ymd,hms,yr,mo,da,hr,mn,sc
      INTEGER, DIMENSION(1)			:: starts_1d,counts_1d,dims_1d
      INTEGER, DIMENSION(2)			:: starts_2d,counts_2d,dims_2d
      REAL*8, DIMENSION(:), ALLOCATABLE	:: avetime
      CHARACTER(LEN=*), PARAMETER 		:: rname = mname//'/obspack_write_output'
      CHARACTER(LEN=16) 			:: stamp

      !
      ! !INPUT PARAMETERS:
      !
      LOGICAL,        INTENT(IN)  		:: am_I_Root   ! Is this the root CPU?

      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)		:: RC          ! Success or failure?

      RC            =  GIGC_SUCCESS

      WRITE (6,'("obspack_write_output")') 
      FLUSH(6)
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### OBSPACK_WRITE_OUTPUT: starting' )
      ENDIF

      IF (nobs .EQ. 0 ) THEN
         RETURN
      END IF

      ! Compute averages

      DO iobs=1,nobs
         IF(obspack_nsamples(iobs) > 0) THEN
            obspack_u(iobs)=obspack_u(iobs)/obspack_nsamples(iobs)
            obspack_v(iobs)=obspack_v(iobs)/obspack_nsamples(iobs)
            obspack_blh(iobs)=obspack_blh(iobs)/obspack_nsamples(iobs)
            obspack_q(iobs)=obspack_q(iobs)/obspack_nsamples(iobs)
            obspack_temperature(iobs)=obspack_temperature(iobs)/obspack_nsamples(iobs)
            obspack_pressure(iobs)=obspack_pressure(iobs)/obspack_nsamples(iobs)
            obspack_tracers(iobs,:)=obspack_tracers(iobs,:)/obspack_nsamples(iobs)
         ENDIF
      ENDDO

      WRITE (6,'("[obspack_write_output] Creating ",a,"...")') obspack_outfile
      ! Start file

      CALL NcCr_Wr( ncid, obspack_outfile)

      ! Turn filling off
      CALL NcSetFill( ncid, NF_NOFILL, omode )

      ! Define dimensions

      CALL NcDef_Dimension( ncid, 'obs', NF_UNLIMITED, dimid_obs )
      CALL NcDef_Dimension( ncid, 'tracer', ntracers,  dimid_tracer)
!      CALL NcDef_Dimension( ncid, 'tracer_name_len', tracer_namelen,  dimid_tnmlen)
      CALL NcDef_Dimension( ncid, 'char100', 100,  dimid_char100)

      ! Set global attributes

      stamp = SYSTEM_TIMESTAMP()
      WRITE(attstring, '("GEOS-Chem simulation at ",a)') stamp
      CALL NcDef_Glob_Attributes( ncid, 'History',           trim(attstring) )
      CALL NcDef_Glob_Attributes( ncid, 'Conventions',       'CF-1.4'  )

      ymd=GET_NYMDB()
      hms=GET_NHMSB()
      CALL YMD_EXTRACT(ymd,yr,mo,da)
      CALL YMD_EXTRACT(hms,hr,mn,sc)
      WRITE(attstring,'(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC")') yr,mo,da,hr,mn,sc
      CALL NcDef_Glob_Attributes( ncid, 'model_start_date',  trim(attstring)   )

      ymd=GET_NYMDE()
      hms=GET_NHMSE()
      CALL YMD_EXTRACT(ymd,yr,mo,da)
      CALL YMD_EXTRACT(hms,hr,mn,sc)
      WRITE(attstring,'(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC")') yr,mo,da,hr,mn,sc
      CALL NcDef_Glob_Attributes( ncid, 'model_end_date',  trim(attstring)   )

      ! Define variables and attributes

      dims_2d = (/ dimid_char100, dimid_obs /)
      vid     = 0
      CALL NcDef_Variable      ( ncid, 'obspack_id', NF_CHAR,  2, dims_2d, vid  )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name',      'obspack_id'      )
      CALL NcDef_Var_Attributes( ncid, vid, 'units',          'unitless'                ) 

      vid  = vid + 1
      dims_2d = (/ dimid_tracer, dimid_obs /) 
      CALL NcDef_Variable( ncid, 'flask', NF_DOUBLE, 2, dims_2d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'mole_fraction_of_trace_gas_in_air' )
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'mol mol-1'      )
      CALL NcDef_Var_Attributes( ncid, vid, '_FillValue', -1e34          )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'nsamples', NF_INT, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'no. of model samples')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'unitless' )
      CALL NcDef_Var_Attributes( ncid, vid, 'comment',     'Number of discrete model samples in average.' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'averaging_time', NF_INT, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'averaging time')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'seconds' )
      CALL NcDef_Var_Attributes( ncid, vid, 'comment',     'Amount of model time over which this sample is averaged.' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'u', NF_DOUBLE, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'u-wind')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'm s^-1' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'v', NF_DOUBLE, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'v-wind')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'm s^-1' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'blh', NF_DOUBLE, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'v-wind')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'm s^-1' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'q', NF_DOUBLE, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'mass_fraction_of_water_inair')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'kg water (kg air)^-1' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'pressure', NF_DOUBLE, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'pressure')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'Pa' )

      vid  = vid + 1
      dims_1d = (/ dimid_obs /)
      CALL NcDef_Variable( ncid, 'temperature', NF_DOUBLE, 1, dims_1d, vid )
      CALL NcDef_Var_Attributes( ncid, vid, 'long_name', 'temperature')
      CALL NcDef_Var_Attributes( ncid, vid, 'units',     'K' )


      CALL NcEnd_def( ncid )

      starts_2d = (/ 1,      1      /)
      counts_2d = (/ 100,    nobs   /)
      CALL NcWr( obspack_id, ncid, 'obspack_id', starts_2d, counts_2d )

      starts_2d = (/ 1,           1      /)
      counts_2d = (/ ntracers,    nobs   /)
      CALL NcWr( transpose(obspack_tracers), ncid, 'flask', starts_2d, counts_2d )

      !PRINT *,'FLASKS:',obspack_tracers

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr( obspack_nsamples, ncid, 'nsamples', starts_1d, counts_1d )

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      ALLOCATE(avetime(nobs))
      avetime=obspack_tau_end-obspack_tau_start
      CALL NcWr( avetime, ncid, 'averaging_time', starts_1d, counts_1d )
      DEALLOCATE(avetime)

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr( obspack_u, ncid, 'u', starts_1d, counts_1d )

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr(  obspack_v, ncid, 'v', starts_1d, counts_1d )

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr( obspack_blh, ncid, 'blh', starts_1d, counts_1d )

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr( obspack_q, ncid, 'q', starts_1d, counts_1d )

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr( obspack_pressure, ncid, 'pressure', starts_1d, counts_1d )

      starts_1d = (/ 1      /)
      counts_1d = (/ nobs   /)
      CALL NcWr( obspack_temperature, ncid, 'temperature', starts_1d, counts_1d )

      ! Close the netCDF file
      CALL NcCl( ncid )


      CALL OBSPACK_CLEANUP(RC)

    END SUBROUTINE OBSPACK_WRITE_OUTPUT
    !EOC

    !BOP
    ! !IROUTINE: OBSPACK_SAMPLE
    !
    ! !DESCRIPTION: Subroutine OBSPACK\_SAMPLE performs the model sampling
    !  and saves concentrations to locations corresponding to a flight
    !  track.  \\ \\ 
    !INTERFACE:
    !
    SUBROUTINE OBSPACK_SAMPLE(am_I_Root,Input_Opt,State_Met,State_Chm,RC)
      !
      ! !USES:
      !
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,     ONLY : GET_TAU,        GET_TS_DIAG
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


      REAL*8           :: TCVV(Input_Opt%N_TRACERS)
 
      ! !REVISION HISTORY: 
      !  08 Jun 2015 - A. Jacobson, A. Schuh - imported from Andrew Schuh's
      !                              ct_mod.F, itself modified from
      !                              planeflight_mod.F
      !
      !EOP
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      INTEGER             :: I, J, L
      INTEGER             :: itracer, iobs
      REAL*8            :: THIS_TAUS, THIS_TAUE

      !
      ! !DEFINED PARAMETERS:
      !


      ! Pointers to MetState objects
      REAL*8, POINTER :: MET_SURFACE_HEIGHT(:,:)
      REAL*8, POINTER :: MET_U(:,:,:)
      REAL*8, POINTER :: MET_V(:,:,:)
      REAL*8, POINTER :: MET_PBLH(:,:)
      REAL*8, POINTER :: MET_Q(:,:,:)
      REAL*8, POINTER :: MET_PRESSURE(:,:,:)
      REAL*8, POINTER :: MET_TEMPERATURE(:,:,:)
      REAL*8, POINTER :: STT(:,:,:,:)
      REAL*8, POINTER :: AD(:,:,:)


      !=================================================================
      ! OBSPACK_SAMPLE begins here
      !=================================================================

!      PRINT*,'obspack sample'
!      flush(6)

      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### OBSPACK_SAMPLE: starting' )
      ENDIF

      ! Assume success                                                                                                                                                                                            
      RC            =  GIGC_SUCCESS

      TCVV          = Input_Opt%TCVV(1:ntracers)


      ! Return if ObsPack sampling is turned off (perhaps
      ! because there are no data at this time).
      IF ( .not. DO_OBSPACK ) RETURN

      DO iobs = 1, nobs

         THIS_TAUE = GET_TAU()
         THIS_TAUS = THIS_TAUE - ( GET_TS_DIAG() / 60d0 )

            !WRITE (6,'("TESTING obs ",i,", obspack_id: ",a)') iobs,trim(obspack_id(iobs))
            !FLUSH(6)

               !PRINT *,'obspack_tau_start:', obspack_tau_start(iobs)
               !PRINT *,'obspack_tau_end:', obspack_tau_end(iobs)
               !PRINT *,'this_tau_start:', THIS_TAUS
               !PRINT *,'this_tau_end:', THIS_TAUE

         IF (OBSPACK_TAU_START(iobs) <= THIS_TAUS .and. OBSPACK_TAU_END(iobs) >= THIS_TAUE ) THEN

!            WRITE (6,'("   ")')
!            FLUSH(6)
!            WRITE (6,'("   ")')
!            FLUSH(6)
!            WRITE (6,'("sampling obs ",i,", obspack_id: ",a)') iobs,trim(obspack_id(iobs))
!            FLUSH(6)

            ! Return grid box indices for the chemistry region
            CALL OBSPACK_GET_GRID_INDICES( &
                 obspack_longitude(iobs), &
                 obspack_latitude(iobs), &
                 obspack_altitude(iobs), &
                 I, J, L, State_Met, RC )

 !           WRITE (6,'("lon ",f6.1,", lat ",f5.1,", alt ",f6.1," is IJL ",i4,i4,i4,".")') &
 !                obspack_longitude(iobs), &
 !                obspack_latitude(iobs), &
 !                obspack_altitude(iobs), &
 !                I,J,L
 !           FLUSH(6)

            ! Initialize pointers

            AD                =>      State_Met%AD
            STT               =>      State_Chm%Tracers
            MET_U             =>      State_Met%U
            MET_V             =>      State_Met%V
            MET_Q             =>      State_Met%SPHU
            MET_PBLH          =>      State_Met%PBLH
            MET_PRESSURE      =>      State_Met%PMID
            MET_TEMPERATURE   =>      State_Met%T

            DO itracer = 1,ntracers

               obspack_tracers(iobs,itracer) = obspack_tracers(iobs,itracer) + &
                    STT(I,J,L,itracer)*TCVV(itracer)/AD(I,J,L)

               !WRITE (6,'("STT: ",g,", TCVV: ",f,", AD: ",g,", result: ",g)') &
               !     STT(I,J,L,itracer),TCVV(itracer), AD(I,J,L), &
               !     obspack_tracers(iobs,itracer)
               !FLUSH(6)

            ENDDO

            obspack_u(iobs)   = obspack_u(iobs) + MET_U(I,J,L) 

            obspack_v(iobs)   = obspack_v(iobs) + MET_V(I,J,L) 

            obspack_blh(iobs) = obspack_blh(iobs) + MET_PBLH(I,J) 

            obspack_q(iobs) = obspack_q(iobs) + MET_Q(I,J,L) 

            obspack_pressure(iobs) = obspack_pressure(iobs) + MET_PRESSURE(I,J,L) 

            obspack_temperature(iobs) = obspack_temperature(iobs) + MET_TEMPERATURE(I,J,L) 

            obspack_nsamples(iobs) = obspack_nsamples(iobs) + 1

            ! Free pointer
            NULLIFY( STT )

            NULLIFY(MET_SURFACE_HEIGHT)
            NULLIFY( MET_U )
            NULLIFY( MET_V )
            NULLIFY( MET_PBLH )
            NULLIFY( MET_Q )
            NULLIFY( MET_PRESSURE )
            NULLIFY( MET_TEMPERATURE )

         ENDIF
      ENDDO

    END SUBROUTINE OBSPACK_SAMPLE
    !EOC

    !BOP
    ! !IROUTINE: CT
    !
    ! !DESCRIPTION: Subroutine OBSPACK\_GET\_GRID\_INDICES returns the
    !  grid box indices (I, J, L) corresponding to the input point
    !  defined by longitude, latitude, altitude.
    !
    !INTERFACE:
    !
    SUBROUTINE OBSPACK_GET_GRID_INDICES(longitude, latitude, altitude, I, J, L, State_Met, RC)

      USE GRID_MOD,		ONLY : GET_XOFFSET
      USE GRID_MOD,		ONLY : GET_YOFFSET
      USE CMN_SIZE_MOD,		ONLY : LLPAR, DISIZE, DJSIZE, IIPAR
      USE GIGC_State_Met_Mod,	ONLY : MetState
      USE PRESSURE_MOD, 	ONLY : GET_PEDGE

      !
      ! !OUTPUT PARAMETERS:
      !
      REAL*8,INTENT(IN)		:: longitude
      REAL*8,INTENT(IN)		:: latitude
      REAL*8,INTENT(IN)		:: altitude
      TYPE(MetState), INTENT(IN)   	:: State_Met   ! Meteorology State object

      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER, INTENT(OUT)		:: RC          ! Success or failure?
      INTEGER, INTENT(OUT)		:: I, J, L

      ! !LOCAL VARIABLES
      !
      INTEGER			:: idx, I0, J0
      REAL*8			:: Z(LLPAR+1)
      REAL*8			:: PEDGE(LLPAR+1)


      RC            =  GIGC_SUCCESS

      ! Added correct definitions for I and J based on nested regions 
      ! (lds, 8/25/11)
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Get I corresponding to PLON(IND)
      I = INT( ( longitude + 180d0 - (I0 * DISIZE) ) / DISIZE + 1.5d0 )

      ! Handle date line correctly (bmy, 4/23/04)
      IF ( I > IIPAR ) I = I - IIPAR

      ! Get J corresponding to PLAT(IND)
       J = INT( ( latitude +  90d0 - (J0 * DJSIZE) ) / DJSIZE + 1.5d0 )

     ! Note use of L here as an index; later it will be a scalar
     ! return value.

      DO idx = 1, LLPAR+1
           PEDGE(idx) = 100.*REAL(GET_PEDGE(I,J,idx) ,4)
      ENDDO

      ! Get array Z of interface edge GPHs
      CALL GEOPOT_HYPSO( State_Met, I, J, LLPAR, Z, PEDGE,RC)

      !PRINT *,'Z:',Z

      DO L = 1,LLPAR
         IF((Z(L+1) >= altitude)) RETURN  ! NOT FINISHED HERE YET
      ENDDO

      !PRINT *,'L chosen:',L

    END SUBROUTINE OBSPACK_GET_GRID_INDICES
    !EOC

    !------------------------------------------------------------------------------
    !                  GEOS-Chem Global Chemical Transport Model                  !
    !------------------------------------------------------------------------------
    !BOP
    !
    ! !IROUTINE: geopot_hypso
    !
    ! !DESCRIPTION: Subroutine GEOPOT\_HYPSO returns the heights above
    !               sea level (m) of the nlvl+1 edges of nlvl model
    !               levels at grid point i,j.  The first geopotential
    !               height is the surface height, the last is
    !               GEOS-Chem's top-of-atmosphere.  We do this by
    !               integrating the hypsometric equation from the
    !               surface upwards.  No humidity information is used in
    !               this computation.
    ! 
    !
    !INTERFACE:
    !
    SUBROUTINE GEOPOT_HYPSO( State_Met, I, J, nlvl, Z, pedge, RC )
      !EOP
      !BOC
      USE GIGC_State_Met_Mod,		ONLY : MetState

      USE CMN_GCTM_MOD,			ONLY : g0
      ! Gravitational ac-
      ! celeration, m/s^2.

      USE CMN_GCTM_MOD,			ONLY : Rd
      ! Dry gas constant,
      ! J / (kg K)


      IMPLICIT NONE

      !
      ! !INPUT PARAMETERS:
      !
      TYPE(MetState), INTENT(IN)                  :: State_Met   ! Meteorology State object
      INTEGER, INTENT(IN)                         :: nlvl        ! number of vertical levels
      INTEGER, INTENT(IN)                         :: I           ! zonal grid index of point
      INTEGER, INTENT(IN)                         :: J           ! meridional grid index of point
      REAL*8, DIMENSION(1:nlvl+1), INTENT(IN)   :: pedge       ! pressures at level boundaries (Pa)

      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER, INTENT(OUT)				:: RC          ! Success or failure?
      REAL*8, DIMENSION(1:nlvl+1), INTENT(OUT)	:: Z           ! geopotential height (m)

      !
      ! !POINTERS
      !
      REAL*8,POINTER                            :: T(:,:,:)    ! Temperature (K)
      REAL*8, POINTER                                 :: PHIS(:,:)

      !
      ! !LOCAL PARAMETERS:
      !
      INTEGER					:: N

      RC            =  GIGC_SUCCESS

      !  Targets to pointers
      T    => State_Met%T
      PHIS => State_Met%PHIS

      Z(1) = PHIS(I,J)/g0

      DO N = 1, NLVL
         Z(N+1) = Z(N)+(Rd*T(I,J,N)/g0)*LOG(pedge(N)/pedge(N+1))
         !PRINT *,'Z(N+1):',Z(N+1)
         !PRINT *,'LOG(pedge(N)/pedge(N+1):',LOG(pedge(N)/pedge(N+1))
         !PRINT *,'T(I,J,N):',T(I,J,N)
      ENDDO

    END SUBROUTINE GEOPOT_HYPSO

    !EOC

  END MODULE OBSPACK_MOD



