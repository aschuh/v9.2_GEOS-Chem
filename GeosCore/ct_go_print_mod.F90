!
! go print : tools for standard output
!
! Example:
!
!   ! messages printed by root only:
!   call GO_Print_Init( status, apply=myid==root, &
!                           prompt_pe=npes>1, pe=myid, &
!                           trace=.false. )
!   if (status/=0) stop
!
!   ! set routine label:
!   call goLabel( 'mymod/myroutine' )
!
!   ! write single message (including processor prompt?) :
!   !   [00] This is number  3
!   write (gol,'("This is number ",i2)') 3; call goPr
!
!   ! write error message and traceback using the
!   ! previous defined routine label:
!   !   [00] ERROR - Something wrong.
!   !   [00] ERROR in mymod/myroutine
!   write (gol,'("Something wrong.")'); call goErr
!   call goErr
!
!   ! close label
!   call goLabel()
!
!   ! done
!   call GO_Print_Done( status )
!   if (status/=0) stop
!

! Nedit macro's:
!
!  o change error traceback:
!       write (*,'("ERROR in ",a)') rname
!       call goErr
!
!  o change error traceback:
!       write (*,'("ERROR in ",a)') rname
!       write (gol,'("in ",a)') rname; call goErr
!
!  o change error message:
!       write \(\*,'\("ERROR - (.*$)
!       write (gol,'("\1; call goErr
!
!  o change other message:
!       write \(\*,(.*$)
!       write (gol,\1; call goPr
!
!  o change error time messages:
!       (goprdt.*ERROR.*$)
!       \1; call goErr
!
!  o change time messages:
!       printdate2
!       wrtgol
!       call printdate2\( 'ERROR - (.*$)
!       call wrtgol( '\1; call goErr
!       printdate
!       wrtgol
!       call printdate\( 'ERROR - (.*$)
!       call wrtgol( '\1; call goErr
!
!
!
!  o change error messages:
!       (ERROR.*; call )goPr
!       \1goErr
!
!  o change time messages:
!       (call goprdt.*$)
!       \1; call goPr
!

module CT_GO_Print_mod

  implicit none

  ! --- in/out ---------------------------------

  !private

  !public  ::  gol


  ! --- const ---------------------------------

  character(len=*), parameter  ::  mname = 'GO_Print'
  integer, parameter           ::  nregions=240
  integer, parameter           ::  nlon360=360
  integer, parameter           ::  nlat180=180
  integer, parameter           ::  nread=10800 ! 3*3600

  integer, parameter           :: ico2tot=1
  integer, parameter           :: ico2bg=6
  integer, parameter           :: ico2ff=2
  integer, parameter           :: ico2bio=4
  integer, parameter           :: ico2ocean=3
  integer, parameter           :: ico2fires=5

  integer, parameter           ::  iyear0=1980  
  
  !-- molar weight of carbon?
  real,    parameter :: xmc  = 12.01115

  !-- mass of air
  real, parameter        ::  xmair    =  28.94            ! g/mol; old name!

  !-- molar weight of co2
  real, parameter               :: xmco2 = 44.

  ! --- var ------------------------------------

  ! buffer for standard output
  character(len=1024)  ::  gol

  ! stack with labels:
  integer, parameter   ::  mstack = 400
  character(len=64)    ::  labels(0:mstack)
  integer              ::  istack = 0

  ! initialized ?
  ! some errors might be printed before initialization ...
  logical              ::  pr_initialized = .false.

  ! destination file unit:
  integer              ::  pr_fu

  ! flags etc
  logical              ::  pr_apply
  logical              ::  pr_trace

  ! processor prompt
  logical              ::  pr_prompt_pe
  integer              ::  pr_pe

  ! white space for indents:
  integer, parameter          ::  dindent = 2
  integer                     ::  indent = 0

  ! writ to file ?
  logical              ::  pr_file
  character(len=256)   ::  pr_file_name

  !-- Types

  type field_1x1
     real*8, dimension(nlon360,nlat180)                 :: surf
     integer                                            :: n
     real*8                                             :: wda
  end type field_1x1


  ! This particular parameter structure is appropriate for a map of scaling factors.

  type enkf_parms
     real, dimension(nlon360,nlat180)                   :: surf
  end type enkf_parms

  type(field_1x1), dimension(:), allocatable,target,save        :: co2_ocn

contains

  subroutine date2tau(idatex,itaux,icalendo)
    !-----------------------------------------------------------------------
    !**** date2tau
    !
    !     purpose
    !     -------
    !     calculate time in seconds from given date
    !
    !     parameters
    !     ----------
    !     on input : idatex contains date in year,month,day,hour,min,sec
    !     on output: itaux contains date/time in seconds
    !
    !     dependencies
    !     ------------
    !     icalendo determines the type calendar used for the conversion
    !     iyear0 is the reference year for the calculation
    !     julday0 is the reference julian day for the calculation
    !
    !     externals
    !     ---------
    !     funtions: julday
    !-----------------------------------------------------------------------
    !use dims, only : iyear0, julday0, kmain
    use julday_mod, only : julday

    implicit none

    ! input/output
    integer,intent(in)              :: icalendo
    integer,dimension(6),intent(in) :: idatex
    integer,intent(out)             :: itaux

    ! local
    integer :: idaysec
    !
    ! compute the seconds the day is old
    !
    idaysec=idatex(6)+idatex(5)*60+idatex(4)*3600
    !
    ! permanent 360 year calendar with 30 days in each month
    !
    if ( icalendo == 1 ) then
       itaux=idaysec+(idatex(3)-1)*86400+(idatex(2)-1)*2592000 &
            +(idatex(1)-iyear0)*31104000
       !
       ! real calendar
       !
    else if ( icalendo == 2 ) then
       !itaux=86400*(julday(idatex(2),idatex(3),idatex(1))-julday0)+idaysec
       itaux=FLOOR(86400*(julday(idatex(1),idatex(2),REAL(idatex(3),8))-julday(1980,1,REAL(1.,8))+idaysec))
       !
       ! permanent 365 day year calendar
       !
    else if ( icalendo == 3 ) then
       itaux=FLOOR(86400*(julday(1981,idatex(2),REAL(idatex(3),8))-julday(1981,1,REAL(1,8)))  &
            +(idatex(1)-iyear0)*365*86400+idaysec)
       !
       ! permanent leap year calendar
       !
    else if ( icalendo == 4 ) then
       itaux=FLOOR(86400*(julday(1980,idatex(2),REAL(idatex(3),8))-julday(1980,1,REAL(1,8)))  &
            +(idatex(1)-iyear0)*366*86400+idaysec)
       !
       ! illegal option icalendo
       !
    else
       !write(kmain,*) ' date2tau: ERROR while computing time'
       !write(kmain,*) ' date2tau: Illegal calendar type'
       !write(kmain,*) '           icalendo = ',icalendo
       print *,' date2tau: ERROR while computing time'
       print *,' date2tau: Illegal calendar type'
       print *,'           icalendo = ',icalendo
       stop
    end if

  end subroutine date2tau

     !---  FUNCTION TO CONVERT INTEGER DATE TO DECIMAL

real*8 function idate2ddate ( idate )

    ! computes the decimal date given a calendar date
    ! 21 Oct 2011, M. Trudeau

    implicit none

    integer, dimension(6), intent(in) :: idate
    integer                           :: siy, itau, itau_ref, itau_bgn, itau_end
    logical                           :: is_leap

    ! *** setting icalendo to 2 ****
    call date2tau((/idate(1), 1, 1, 0, 0, 0/), itau_ref, 2)
    call date2tau(idate, itau, 2)

    is_leap = (mod(idate(1), 4) == 0 .and. .not. mod(idate(1), 100) == 0) .or. (mod(idate(1), 400) == 0)

    if ( is_leap ) then
       siy = 31622400
    else
       siy = 31536000
    endif

    idate2ddate = dble(idate(1)) + dble(itau - itau_ref) / dble(siy)

  end function idate2ddate

integer function month_length(mon,LEAPYEAR_TF)

    implicit none

    ! input/output
    integer,intent(in)              :: mon
    logical,intent(in)              :: LEAPYEAR_TF
    integer                         :: days

    if(mon==1 .or. mon==3 .or. mon==5 .or. mon==7 .or. mon==8 .or. mon==10 .or. mon==12) THEN
    
       days = 31

    ELSE IF (mon==4 .or. mon==6 .or. mon==9 .or. mon==11)  THEN
       
       days = 30

    ELSE IF (mon==2) THEN

       IF(LEAPYEAR_TF) THEN 
     
         days = 29

       ELSE  
         days = 28
       END IF
    ELSE 
        stop 'ERROR in MONTH_LENGTH()'
    END IF

   month_length = days

   end function month_length

end module ct_go_print_mod

