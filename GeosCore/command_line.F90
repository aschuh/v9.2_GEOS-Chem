! Parses the four command line arguments for controlling ensemble runs.
SUBROUTINE PARSE_COMMAND_LINE
    USE ENKF_MOD
    USE TIME_MOD, ONLY: SET_BEGIN_TIME, SET_END_TIME

    INTEGER :: ARG_COUNT, DAYS
    INTEGER :: START_DATE, END_DATE
    CHARACTER(80) :: ARG

    INTEGER, PARAMETER :: MAX_CYCLE_DAYS = 367
    INTEGER, PARAMETER :: MAX_CYCLES = 200

    ARG_COUNT = COMMAND_ARGUMENT_COUNT()

    ! Decide what to do based on argument count
    IF (ARG_COUNT == 0) THEN
        PRINT *, "No command line arguments given"
        ENUMBER = NO_ENSEMBLE_NUM
        CYCLE_NUM = NO_CYCLES
        RETURN
    END IF

    ! ensemble number or help
    CALL GET_COMMAND_ARGUMENT(1, ARG)
    IF (ARG == '-h' .OR. ARG == '-?' .OR. ARG == '-help' .OR. ARG == '--help') THEN
        CALL PRINT_USAGE(FULL = .TRUE.)
    END IF

    IF(ARG_COUNT /= 6) CALL PRINT_USAGE

    READ(ARG, *) ENUMBER
    IF (ENUMBER == 0) THEN
        PRINT *, "Running as ensemble control"
    ELSE IF (ENUMBER > 0 .AND. ENUMBER < MAX_ENSEMBLE_MEMBERS) THEN
        WRITE(*, '(A,I4)'), "Running as ensemble #", ENUMBER
    ELSE
        STOP "Bad ensemble number"
    END IF

    ! cycle number
    CALL GET_COMMAND_ARGUMENT(2, ARG)
    READ(ARG, *) CYCLE_NUM
    IF (CYCLE_NUM < 0 .OR. CYCLE_NUM > MAX_CYCLES) THEN
        STOP "Bad cycle number"
    END IF

    PRINT *, "Cycle number =", CYCLE_NUM

    ! lag number
    CALL GET_COMMAND_ARGUMENT(3, ARG)
    READ(ARG, *) LAG_NUM
    IF (LAG_NUM < 1 .OR. LAG_NUM > 6) THEN
        STOP "Bad lag number"
    END IF

    PRINT *, "Lag number =", LAG_NUM

    ! start date
    CALL GET_COMMAND_ARGUMENT(4, ARG)
    READ(ARG, *) START_DATE
    CALL SET_BEGIN_TIME(START_DATE, 0)
    PRINT *, "Cycle start date =", START_DATE

    ! cycle length in days
    CALL GET_COMMAND_ARGUMENT(5, ARG)
    READ(ARG, *) DAYS
    READ(ARG, *) CYCLE_LEN
    IF (0 < DAYS .AND. DAYS <= MAX_CYCLE_DAYS) THEN
        CALL INCREMENT_DATE(START_DATE, DAYS, END_DATE)
        CALL SET_END_TIME(END_DATE, 0)
        PRINT *, "Cycle end date   =", END_DATE
    ELSE
        STOP "Bad cycle count"
    END IF

    ! forward opt mean run?
    CALL GET_COMMAND_ARGUMENT(6, ARG)
    READ(ARG, *) RERUN

    IF (RERUN .NE. 0 .AND. RERUN .NE. 1) THEN
        STOP "Bad rerun number"
    END IF

    IF(RERUN) THEN
      PRINT *, "Running optimized mean forward (serial)"
    ELSE
      PRINT *, "Running full ensemble forward"
    ENDIF

    PRINT *, ""
CONTAINS

    SUBROUTINE PRINT_USAGE(FULL)
        LOGICAL, INTENT(IN), OPTIONAL :: FULL

        PRINT *, "Usage: geos [ENSEMBLE_NUMBER CYCLE_NUM LAG_NUM START_DATE CYCLE_LEN RERUN]"

        IF (PRESENT(FULL)) THEN
            PRINT *, ""
            PRINT *, "Where:"
            PRINT *, "  ENSEMBLE_NUMBER - the rank of this member in the ensemble"
            PRINT *, "  CYCLE_NUM       - the cardinal number for this cycle"
            PRINT *, "  LAG_NUM         - the lag #, 1=last est, 6=first est"
            PRINT *, "  START_DATE      - the starting date for this cycle run.  An integer of"
            PRINT *, "                    format YYYYMMDD as in the input.geos file. Overrides"
            PRINT *, "                    the value in input.geos because an ensemble run will"
            PRINT *, "                    go through different start dates in each cycle."
            PRINT *, "  CYCLE_LEN       - the length of this cycle in days"
            PRINT *, "  RERUN           - 0: running normal ensemble "
            PRINT *, "                 -  1: rerunning serial run with mean optimized flux"
            PRINT *, "                      in order to get restart CO2 for next cycle"
        END IF

        STOP ""
    END SUBROUTINE PRINT_USAGE

    SUBROUTINE INCREMENT_DATE(IN_DATE, DAYS, OUT_DATE)
        USE TIME_MOD, ONLY: YMD_EXTRACT
        USE JULDAY_MOD

        INTEGER, INTENT(IN)  :: IN_DATE, DAYS
        INTEGER, INTENT(OUT) :: OUT_DATE

        INTEGER :: YEAR, MONTH, IDAY, IHMS
        REAL*8 :: DAY, JD

        ! extract y m d
        CALL YMD_EXTRACT(IN_DATE, YEAR, MONTH, IDAY)
        DAY = IDAY

        ! increment julian days
        JD = JULDAY(YEAR, MONTH, DAY) + DAYS

        ! get calendar day again
        CALL CALDATE(JD, OUT_DATE, IHMS)
    END SUBROUTINE INCREMENT_DATE

END SUBROUTINE PARSE_COMMAND_LINE
