!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Carbontracker_In_Mod
!
! !DESCRIPTION: Module GIGC\_STATE\_MET\_MOD contains the derived type
!  used to define the Meteorology State object for the Grid-Independent
!  GEOS-Chem implementation (abbreviated "GIGC").
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory
!  to the Carbontracker input object.  The Carbontracker input object is not
!  defined in this module.  It must be be declared as variable in the top-level
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_In_Carbontracker_Mod
!
! USES:
!

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_Carbontracker_In
  PUBLIC :: Init_Carbontracker_Out
  PUBLIC :: Cleanup_Carbontracker_In
  PUBLIC :: Cleanup_Carbontracker_Out
  PUBLIC :: COPY_INPUTS
  PUBLIC :: COPY_OUTPUTS
  PUBLIC :: REINIT_OUTPUT
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Meteorology State
  !=========================================================================
  TYPE, PUBLIC :: Carbontracker_In
     !----------------------------------------------------------------------
     ! netcdf integer fields
     !----------------------------------------------------------------------
     INTEGER,  POINTER :: NPOINTS
     INTEGER,  POINTER :: TIME      (:  )  ! Time 
     INTEGER,  POINTER :: SAMPLING_STRATEGY     (:  )  

     !-----------------------------------------------------------------------
     ! netcdf real4  fields
     !-----------------------------------------------------------------------

     REAL,  POINTER       :: TAU_START  (: )
     REAL,  POINTER       :: TAU_END    (: )
 
     REAL,  POINTER       :: LATITUDE     (:  )  
     REAL,  POINTER       :: LONGITUDE     (:  )  
     REAL,  POINTER       :: ALTITUDE     (:  )  
     REAL,  POINTER       :: VALUE     (:  )  

     !-----------------------------------------------------------------------
     ! netcdf char  fields
     !-----------------------------------------------------------------------
 
     CHARACTER,  POINTER       :: ID     (: ,: )  

  END TYPE Carbontracker_In

  !=========================================================================
  ! Derived type for Meteorology State
  !=========================================================================
  TYPE, PUBLIC :: Carbontracker_Out
     !----------------------------------------------------------------------
     ! netcdf integer fields
     !----------------------------------------------------------------------
     INTEGER, POINTER  :: NPOINTS  
     INTEGER,  POINTER :: NSAMPLES      (:  )  
     INTEGER,  POINTER :: AVERAGING_TIME  (:  )
     INTEGER,  POINTER :: REGION_INDICES     (: ,:  )

     !-----------------------------------------------------------------------
     ! netcdf real4  fields
     !-----------------------------------------------------------------------

     REAL*8,  POINTER       :: SURFACE_HEIGHT     (:  )
     REAL*8,  POINTER       :: FLASK     (: ,: )
     REAL*8,  POINTER       :: U     (:  )
     REAL*8,  POINTER       :: V     (:  )
     REAL*8,  POINTER       :: BLH     (:  )
     REAL*8,  POINTER       :: Q     (:  )
     REAL*8,  POINTER       :: PRESSURE     (:  )
     REAL*8,  POINTER       :: TEMPERATURE     (:  )
     !-----------------------------------------------------------------------
     ! netcdf char  fields
     !-----------------------------------------------------------------------

     CHARACTER,  POINTER       :: ID     (: ,: )
     CHARACTER,  POINTER       :: REGION_NAME     (: ,: )

  END TYPE Carbontracker_Out
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, split off from gc_type_mod.F90
!  23 Oct 2012 - R. Yantosca - Added QI, QL met fields to the derived type
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  12 Dec 2012 - R. Yantosca - Add IREG, ILAND, IUSE fields for dry deposition
!  13 Dec 2012 - R. Yantosca - Add XLAI, XLAI2 fields for dry deposition
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  15 Nov 2013 - R. Yantosca - Now denote that RH fields have units of [%]
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_Carbontracker_In
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_STATE\_MET allocates all fields of
!  the Grid-Indpendent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Carbontracker_In( am_I_Root, NOBS, In_Carbontracker, RC )
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
    TYPE(Carbontracker_In), INTENT(INOUT) :: In_Carbontracker   ! Obj for meteorology state
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

    !=======================================================================
    ! Allocate Fields
    !=======================================================================
    ALLOCATE( In_Carbontracker%NPOINTS   , STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%NPOINTS     = 0

    ALLOCATE( In_Carbontracker%TIME      ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%TIME     = 0

    ALLOCATE( In_Carbontracker%SAMPLING_STRATEGY    ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%SAMPLING_STRATEGY   = 0

    ALLOCATE( In_Carbontracker%TAU_START      ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%TAU_START    = 0d0

    ALLOCATE( In_Carbontracker%TAU_END      ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%TAU_END    = 0d0

    ALLOCATE( In_Carbontracker%LATITUDE      ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%LATITUDE    = 0d0

    ALLOCATE( In_Carbontracker%LONGITUDE    ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%LONGITUDE   = 0d0

    ALLOCATE( In_Carbontracker%ALTITUDE    ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%ALTITUDE   = 0d0

    ALLOCATE( In_Carbontracker%VALUE   ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%VALUE  = 0d0

    ALLOCATE( In_Carbontracker%ID   ( 100, NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    In_Carbontracker%ID  = ''


  END SUBROUTINE Init_Carbontracker_In
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_Carbontracker_Out
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_STATE\_MET allocates all fields of
!  the Grid-Indpendent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Carbontracker_Out( am_I_Root, NOBS,NTRACERS, Out_Carbontracker, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes

!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: NOBS        ! # number of obs in file
    INTEGER,        INTENT(IN)    :: NTRACERS        ! # number of obs in file
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Carbontracker_Out), INTENT(INOUT) :: Out_Carbontracker   ! Obj for meteorology state
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

    PRINT *,'Initializing CT Output object....with ',NOBS,' obs and ',NTRACERS,' tracers'
  
    !=======================================================================
    ! Allocate Fields
    !=======================================================================
    ALLOCATE( Out_Carbontracker%NPOINTS  , STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%NPOINTS   = 0

    ALLOCATE( Out_Carbontracker%NSAMPLES    ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%NSAMPLES   = 0

    ALLOCATE( Out_Carbontracker%AVERAGING_TIME ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%AVERAGING_TIME   = 0

    ALLOCATE( Out_Carbontracker%SURFACE_HEIGHT      ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%SURFACE_HEIGHT    = 0d0

    ALLOCATE( Out_Carbontracker%REGION_NAME    ( NOBS, 10 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%REGION_NAME   = ''

    ALLOCATE( Out_Carbontracker%REGION_INDICES    ( NOBS, 3 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%REGION_INDICES   = 0

    ALLOCATE( Out_Carbontracker%FLASK   ( NOBS, NTRACERS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%FLASK  = 0d0

    ALLOCATE( Out_Carbontracker%ID   ( 100 , NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%ID  = ''

    ALLOCATE( Out_Carbontracker%U     ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%U    = 0d0

    ALLOCATE( Out_Carbontracker%V     ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%V    = 0d0

    ALLOCATE( Out_Carbontracker%BLH     ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%BLH    = 0d0

    ALLOCATE( Out_Carbontracker%Q     ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%Q    = 0d0

    ALLOCATE( Out_Carbontracker%PRESSURE     ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%PRESSURE    = 0d0

    ALLOCATE( Out_Carbontracker%TEMPERATURE     ( NOBS ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    Out_Carbontracker%TEMPERATURE    = 0d0

    PRINT *,'SHAPE:',SHAPE(Out_Carbontracker%FLASK)
  END SUBROUTINE Init_Carbontracker_Out
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COPY_INPUTS
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COPY_INPUTS( am_I_Root, CT_Input2, CT_Input1 )
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
    TYPE(Carbontracker_In), INTENT(INOUT) :: CT_Input1,CT_Input2   ! Obj for meteorology state
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !========================================================================
    ! These met fields are used for all data products
    !========================================================================
    CT_Input1%NPOINTS = CT_Input2%NPOINTS
    CT_Input1%TAU_START = CT_Input2%TAU_START
    CT_Input1%TAU_END   = CT_Input2%TAU_END
    CT_Input1%TIME = CT_Input2%TIME
    CT_Input1%SAMPLING_STRATEGY = CT_Input2%SAMPLING_STRATEGY   
    CT_Input1%LATITUDE = CT_Input2%LATITUDE
    CT_Input1%LONGITUDE = CT_Input2%LONGITUDE
    CT_Input1%ALTITUDE = CT_Input2%ALTITUDE
    CT_Input1%VALUE = CT_Input2%VALUE
    CT_Input1%ID = CT_Input2%ID

   END SUBROUTINE COPY_INPUTS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COPY_OUTPUTS
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COPY_OUTPUTS( am_I_Root, CT_Output2, CT_Output1 )
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
    TYPE(Carbontracker_Out), INTENT(INOUT) :: CT_Output1,CT_Output2   ! Obj for meteorology state
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !========================================================================
    ! These met fields are used for all data products
    !========================================================================
    CT_Output1%NSAMPLES = CT_Output2%NSAMPLES
    CT_Output1%AVERAGING_TIME = CT_Output2%AVERAGING_TIME
    CT_Output1%SURFACE_HEIGHT = CT_Output2%SURFACE_HEIGHT
    CT_Output1%REGION_NAME = CT_Output2%REGION_NAME
    CT_Output1%REGION_INDICES = CT_Output2%REGION_INDICES
    CT_Output1%FLASK = CT_Output2%FLASK
    CT_Output1%ID = CT_Output2%ID
    CT_Output1%U = CT_Output2%U
    CT_Output1%V = CT_Output2%V
    CT_Output1%BLH = CT_Output2%BLH
    CT_Output1%Q = CT_Output2%Q
    CT_Output1%PRESSURE = CT_Output2%PRESSURE
    CT_Output1%TEMPERATURE = CT_Output2%TEMPERATURE

   END SUBROUTINE COPY_OUTPUTS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: REINIT_OUTPUT
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE REINIT_OUTPUT( am_I_Root, CT_Output )
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
    TYPE(Carbontracker_Out), INTENT(INOUT) :: CT_Output   ! Obj for meteorology state
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !========================================================================
    ! These met fields are used for all data products
    !========================================================================
    CT_Output%NSAMPLES = 0
    CT_Output%AVERAGING_TIME = 0
    CT_Output%SURFACE_HEIGHT = 0.
    CT_Output%REGION_NAME = ''
    CT_Output%REGION_INDICES = 0
    CT_Output%FLASK = 0.
    CT_Output%ID = ''
    CT_Output%U = 0.
    CT_Output%V = 0.
    CT_Output%BLH = 0.
    CT_Output%Q = 0.
    CT_Output%PRESSURE = 0.
    CT_Output%TEMPERATURE = 0.

   END SUBROUTINE REINIT_OUTPUT
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
  SUBROUTINE Cleanup_Carbontracker_In( am_I_Root, In_Carbontracker, RC )
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
    TYPE(Carbontracker_In), INTENT(INOUT) :: In_Carbontracker   ! Obj for meteorology state
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

    !========================================================================
    ! These met fields are used for all data products
    !========================================================================

    IF ( ASSOCIATED( In_Carbontracker%NPOINTS       )) DEALLOCATE( In_Carbontracker%NPOINTS       )
    IF ( ASSOCIATED( In_Carbontracker%TAU_START       )) DEALLOCATE( In_Carbontracker%TAU_START       )
    IF ( ASSOCIATED( In_Carbontracker%TAU_END       )) DEALLOCATE( In_Carbontracker%TAU_END       )

    IF ( ASSOCIATED( In_Carbontracker%TIME       )) DEALLOCATE( In_Carbontracker%TIME       )
    IF ( ASSOCIATED( In_Carbontracker%SAMPLING_STRATEGY       )) DEALLOCATE( In_Carbontracker%SAMPLING_STRATEGY       )
    IF ( ASSOCIATED( In_Carbontracker%LATITUDE       )) DEALLOCATE( In_Carbontracker%LATITUDE       )
    IF ( ASSOCIATED( In_Carbontracker%LONGITUDE       )) DEALLOCATE( In_Carbontracker%LONGITUDE       )
    IF ( ASSOCIATED( In_Carbontracker%ALTITUDE       )) DEALLOCATE( In_Carbontracker%ALTITUDE       )
    IF ( ASSOCIATED( In_Carbontracker%VALUE       )) DEALLOCATE( In_Carbontracker%VALUE       )
    IF ( ASSOCIATED( In_Carbontracker%ID     )) DEALLOCATE( In_Carbontracker%ID       )


   END SUBROUTINE Cleanup_Carbontracker_In
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Carbontracker_Out
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Carbontracker_Out( am_I_Root, Out_Carbontracker, RC )
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
    TYPE(Carbontracker_Out), INTENT(INOUT) :: Out_Carbontracker   ! Obj for meteorology state
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

    !========================================================================
    ! These met fields are used for all data products
    !========================================================================

    IF ( ASSOCIATED( Out_Carbontracker%NSAMPLES       )) DEALLOCATE( Out_Carbontracker%NSAMPLES       )
    IF ( ASSOCIATED( Out_Carbontracker%AVERAGING_TIME  )) DEALLOCATE( Out_Carbontracker%AVERAGING_TIME    )
    IF ( ASSOCIATED( Out_Carbontracker%SURFACE_HEIGHT       )) DEALLOCATE( Out_Carbontracker%SURFACE_HEIGHT       )
    IF ( ASSOCIATED( Out_Carbontracker%REGION_NAME       )) DEALLOCATE( Out_Carbontracker%REGION_NAME       )
    IF ( ASSOCIATED( Out_Carbontracker%REGION_INDICES       )) DEALLOCATE( Out_Carbontracker%REGION_INDICES       )
    IF ( ASSOCIATED( Out_Carbontracker%FLASK       )) DEALLOCATE( Out_Carbontracker%FLASK       )
    IF ( ASSOCIATED( Out_Carbontracker%ID       )) DEALLOCATE( Out_Carbontracker%ID       )
    IF ( ASSOCIATED( Out_Carbontracker%U       )) DEALLOCATE( Out_Carbontracker%U       )
    IF ( ASSOCIATED( Out_Carbontracker%V       )) DEALLOCATE( Out_Carbontracker%V       )
    IF ( ASSOCIATED( Out_Carbontracker%BLH       )) DEALLOCATE( Out_Carbontracker%BLH       )
    IF ( ASSOCIATED( Out_Carbontracker%Q       )) DEALLOCATE( Out_Carbontracker%Q       )
    IF ( ASSOCIATED( Out_Carbontracker%PRESSURE       )) DEALLOCATE( Out_Carbontracker%PRESSURE       )
    IF ( ASSOCIATED( Out_Carbontracker%TEMPERATURE       )) DEALLOCATE( Out_Carbontracker%TEMPERATURE       )

   END SUBROUTINE Cleanup_Carbontracker_Out
!EOC
END MODULE GIGC_In_Carbontracker_Mod
