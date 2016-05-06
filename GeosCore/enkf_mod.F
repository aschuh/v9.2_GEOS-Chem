!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: enkf_mod
!
! !DESCRIPTION: Module ENKF\_MOD contains variables and routines for 
!  controlling the GEOSCHEM model with an Ensemble Kalman Filter. 
!\\
!\\
! !INTERFACE: 
!
      MODULE ENKF_MOD
!
! !USES:
!
      IMPLICIT NONE
!      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean OH field [molec/cm3]
      !REAL*8, PUBLIC, ALLOCATABLE :: OH(:,:,:)
      LOGICAL, PUBLIC             :: LENKF
      LOGICAL, PUBLIC             :: LLAG
      LOGICAL, PUBLIC             :: LBIAS
      CHARACTER(LEN=255), PUBLIC  :: BIASFILE,BIASFILEIN
      INTEGER, PUBLIC             :: ENUMBER
      INTEGER, PUBLIC             :: CYCLE_NUM
      INTEGER, PUBLIC             :: CYCLE_LEN
      INTEGER, PUBLIC             :: LAG_NUM
      INTEGER, PUBLIC             :: RERUN

!-    FOR HARMONIC WORK

      INTEGER, PUBLIC             :: FLUX_TYPE
      INTEGER, PUBLIC             :: HARMONIC
      INTEGER, PUBLIC             :: PFT
      INTEGER, PUBLIC             :: TRANSCOM_REGION

      INTEGER, PUBLIC, PARAMETER :: MAX_ENSEMBLE_MEMBERS = 2000
      INTEGER, PUBLIC, PARAMETER :: ENS_CONTROL_RUN = 0
      INTEGER, PUBLIC, PARAMETER :: NO_ENSEMBLE_NUM = -1
      INTEGER, PUBLIC, PARAMETER :: NO_CYCLES   = -1
      INTEGER, PUBLIC, PARAMETER :: NO_HARMONIC = -1
      INTEGER, PUBLIC, PARAMETER :: NO_TRANSCOM_REGION   = -1
      INTEGER, PUBLIC, PARAMETER :: NO_PFT   = -1
      INTEGER, PUBLIC, PARAMETER :: NO_FLUX_TYPE   = -1

! !PUBLIC MEMBER FUNCTIONS:
!
      !PUBLIC :: GET_GLOBAL_OH
      !PUBLIC :: INIT_GLOBAL_OH
!
! !REVISION HISTORY:
! April 2012 - Created by Andrew Schuh : Colorado State University
!EOP
      END MODULE ENKF_MOD
