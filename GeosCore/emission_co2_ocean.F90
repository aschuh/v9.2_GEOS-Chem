! Time-stamp: <fe5.zeus.fairmont.rdhpcs.noaa.gov:/home/Andy.Jacobson/TM5/proj/ctemis/branches/taka2009/src/emission_co2_ocean__Taka2009.F90: 27 May 2013 (Mon) 19:31:53 UTC>
! Compute air-sea exchange of CO2 in mol m-2 s-1.  Positive is source to atmosphere.
!       This version uses the Takahashi et al. (2009) climatological pCO2 fields.
!
!   - by default, uses modeled atmospheric CO2 mole fraction to compute air-sea exchange
!      -alternatively: use GLOBALVIEW-CO2 marine boundary layer reference CO2 surface
!      -alternatively: use year-2000 pCO2 (air) from Takahashi et al. (2009)
!
!   - by default, applies assumed pCO2(sw) trend from Takahashi et al. (2009)
!      -alternatively: use no trend (appropriate for year 2000)
!

#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if

module emission_co2_ocean

  !use ct_go_print,	only: gol, goErr, goPr !, goLabel
  use ct_go_print_mod,  only: nregions, nlon360,nlat180
  !Not sure, AES
  !use Meteo,		only: Set
  use ct_go_print_mod,  only: enkf_parms ! needed for ocean surface flux ensemble
  use ct_go_print_mod,  only: co2_ocn  ! This isn't right yet, look at more
  !use emission_common,  only:  EmisInputDir
  use ct_go_print_mod, only : field_1x1
  !use emission_common,	only: nmembersloc  !-- This needs to be brought in from enkf_mod.F90
  use enkf_mod, only    : enumber
  use error_mod, only   : error_stop
  USE m_netcdf_io_open
  USE m_netcdf_io_close
  USE m_netcdf_io_read
  USE m_netcdf_io_get_dimlen

  implicit none

  character(len=*), parameter	:: mname = 'emission_co2_ocn'

  private

  public	:: emission_co2_ocean_init
  public        :: emission_co2_ocean_done
  public        :: emission_co2_ocean_calc
  public        :: emission_co2_ocean_calc_oif

  character(len=*), parameter   :: EmisInputDir=''

  ! use this to potentially override the EmisInputDir
  !character(len=2048)		:: OceanInputDir
  character(len=*), parameter   :: OceanInputDir='/discover/nobackup/aschuh/data/ct2013b_input/'

  ! file name
  !character(len=1024)		:: OceanPCO2File
  character(len=*), parameter  :: OceanPCO2File='Takahashi_2009_verC.1x1.nc'
  !character(len=1024)		:: OceanSolFile
  character(len=*), parameter  :: OceanSolFile='WOA09.monthly.co2.sol_and_Sc.180.reordered.nc'
  !character(len=1024)		:: OceanGridMapFile
  character(len=*), parameter  ::  OceanGridMapFile='rij_merra5x4.nc'
  !character(len=1024)		:: GVGridMapFile
  character(len=*), parameter  :: GVGridMapFile='gv_j.nc'
  !character(len=1024)		:: GlobalviewFile
  character(len=*), parameter  :: GlobalviewFile='gv11_mbl.nc'
  !character(len=20)		:: Ocean_pco2_air
  character(len=*), parameter  :: Ocean_pco2_air='interactive'

  character(len=*), parameter  :: deltapco2_prefix='notsure'

  logical			:: apply_sw_trend = .TRUE.

  integer,dimension(:,:),allocatable	:: grid_region
  integer,dimension(:,:),allocatable	:: grid_i
  integer,dimension(:,:),allocatable	:: grid_j

  real,dimension(nlon360,nlat180,1)     :: buff3d
  real,dimension(nlon360,nlat180)	:: schmidtno,solubility,svp,pco2sw_in

  real,dimension(:,:),allocatable	:: gv_mbl
  real,dimension(:),allocatable		:: gv_dd
  integer,dimension(:),allocatable	:: gv_j

contains

  subroutine emission_co2_ocean_init(status)

    !use GO,			only : TrcFile, Init, Done, ReadRc, goLoCase
    !use global_data,		only : rcfile
    !use dims,			only : im,jm,
    use ct_go_print_mod,            only : nregions
    !use dims,                   only : iglbsfc
    !use Meteo,			only : u10m_dat, v10m_dat, ci_dat, m_dat, sp_dat
    !use global_data,		only : mass_dat
    !use MDF,			only : MDF_NETCDF
    !use ReadFromFile_MDF,	only : ReadFromFile

    !USE GIGC_State_Met_Mod, ONLY     : MetState


    implicit none

    character(len=*), parameter ::  rname = mname//'/emission_co2_ocean_init'

    integer, intent(out)	:: status
    !type(TrcFile)		:: rcF

    !TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

    integer                     :: NCID

    !call Init( rcF, rcfile, status )
    !IF_NOTOK_RETURN(status=1)
    !call ReadRc( rcF, 'ocean.input.dir', OceanInputDir ,status, default=EmisInputDir)
    !IF_ERROR_RETURN(status=1)
    !call ReadRc( rcF, 'ocean.pco2.file', OceanPCO2File ,status, default="Takahashi_2009_verC.1x1.nc")
    !IF_ERROR_RETURN(status=1)
    !call ReadRc( rcF, 'ocean.pco2.trend.sw', apply_sw_trend ,status, default=.TRUE.)
    !IF_ERROR_RETURN(status=1)
    !call ReadRc( rcF, 'ocean.solubility.file', OceanSolFile ,status, default="WOA09.monthly.co2.sol_and_Sc.180.nc")
    !IF_ERROR_RETURN(status=1)
    !call ReadRc( rcF, 'ocean.pco2.air', Ocean_pco2_air ,status, default="interactive")
    !IF_ERROR_RETURN(status=1)

    ! Even if pCO2 air condition is takahashi or globalview, we need these r,i,j triplets
    ! to get the surface pressure (to convert dry air mole fraction to partial pressure).
    !call ReadRc( rcF, 'gridmap.1x1.file', OceanGridMapFile ,status, default="rij_glb3x2_nam1x1.nc")
    !IF_ERROR_RETURN(status=1)

    print *,'Starting CT ocean flux initialization....'

    allocate(grid_region(nlon360,nlat180))
    allocate(grid_i(nlon360,nlat180))
    allocate(grid_j(nlon360,nlat180))
    allocate(co2_ocn(2))

    select case ( ocean_pco2_air )

    case ( 'interactive' )
       ! no op

    case ( 'globalview' )
       !call ReadRc( rcF, 'gvmap.1x1.file', GVGridMapFile ,status, default="gv_j.nc")
       !IF_ERROR_RETURN(status=1)
       !call ReadRc( rcF, 'globalview.file', GlobalviewFile ,status, default="gv11_mbl.nc")
       !IF_ERROR_RETURN(status=1)
       allocate(gv_j(nlat180))

    case ( 'takahashi' )
       ! no op

    case default
       !write (gol,'("unsupported ocean.pco2.air value in rc file : """,a,"""")') trim(ocean_pco2_air);call goPr
       !write (gol,'("supported strings are ""interactive"", ""globalview"", and ""takahashi""")');call goPr
       !status=1
       !IF_NOTOK_RETURN(status=1)
       CALL ERROR_STOP('supported strings are interactive, globalview and takahashi',rname)
    end select

    !call Done( rcF, status )
    !IF_NOTOK_RETURN(status=1)

    !-- changed this from printing to gol file to printing to standard out
    !write (gol,'("[emission_co2_ocean_init]	seawater pCO2 data to be read from """,a,"/",a,""".")') trim(OceanInputDir), trim(OceanPCO2File)
    !call goPr

    print *,'[emission_co2_ocean_init]    seawater pCO2 data to be read from ',trim(OceanInputDir),'/',trim(OceanPCO2File)
    if(apply_sw_trend) then
       !write (gol,'("[emission_co2_ocean_init]	   seawater pCO2 will be generated from year-2000 fields plus a 1.5 uatm/yr trend.")')
       print *,'[emission_co2_ocean_init]    seawater pCO2 will be generated from year-2000 fields plus a 1.5 uatm/yr trend.'
    else
       !write (gol,'("[emission_co2_ocean_init]	   seawater pCO2 will be generated from year-2000 fields with NO trend.")')
       print *,'[emission_co2_ocean_init]    seawater pCO2 will be generated from year-2000 fields with NO trend.'
    end if
    !call goPr

    !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanGridMapFile), &
    !     ftype=MDF_NETCDF,vname='r', buffer=grid_region, status=status)
    !IF_NOTOK_RETURN(status=1)

    CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanGridMapFile))
    CALL NcRd(grid_region,NCID,'r', (/ 1,1 /), (/ nlon360,nlat180 /))
    CALL Nccl (NCID)

    !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanGridMapFile), &
    !     ftype=MDF_NETCDF,vname='i', buffer=grid_i, status=status)
    !IF_NOTOK_RETURN(status=1)

    CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanGridMapFile))
    CALL NcRd(grid_i,NCID,'i', (/ 1,1 /), (/ nlon360,nlat180 /))
    CALL Nccl (NCID)

    !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanGridMapFile), &
    !     ftype=MDF_NETCDF,vname='j', buffer=grid_j, status=status)
    !IF_NOTOK_RETURN(status=1)

    CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanGridMapFile))
    CALL NcRd(grid_j,NCID,'j', (/ 1,1 /), (/ nlon360,nlat180 /))
    CALL Nccl (NCID)

    select case ( ocean_pco2_air )

    case ( 'interactive' )
       !write (gol,'("[emission_co2_ocean_init]	pCO2 (air) determined interactively from TM5 modeled mole fractions.")')
       !call goPr
       !write (gol,'("[emission_co2_ocean_init]	   mapping from TM5 grids to 1x1 flux grid to be read from """,a,"/",a,""".")') trim(OceanInputDir), trim(OceanGridMapFile)
       !call goPr

       print *,'[emission_co2_ocean_init] pCO2 (air) determined interactively from TM5 modeled mole fractions.'
       print *,'[emission_co2_ocean_init]    mapping from TM5 grids to 1x1 flux grid to be read from ',trim(OceanInputDir),'/',trim(OceanGridMapFile)

    case ( 'globalview' )
       !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(GVGridMapFile), &
       !     ftype=MDF_NETCDF,vname='j', buffer=gv_j, status=status)
       !IF_NOTOK_RETURN(status=1)

       CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(GVGridMapFile))
       CALL NcRd(gv_j,NCID,'j', (/ 1 /), (/ nlat180 /))
       CALL Nccl (NCID)

       !write (gol,'("[emission_co2_ocean_init]	pCO2 (air) determined from GLOBALVIEW-CO2 dry air mole fraction read from """,a,"/",a,""".")') trim(OceanInputDir), trim(GlobalviewFile)
       !call goPr
       !write (gol,'("[emission_co2_ocean_init]	   mapping from GLOBALVIEW MBL latitudes to 1x1 flux grid to be read from """,a,"/",a,""".")') trim(OceanInputDir), trim(GVGridMapFile)
       !call goPr
       !write (gol,'("[emission_co2_ocean_init]	   mapping from TM5 surface pressure grids to 1x1 flux grid to be read from """,a,"/",a,""".")') trim(OceanInputDir), trim(OceanGridMapFile)
       !call goPr
       print *,'[emission_co2_ocean_init] pCO2 (air) determined from GLOBALVIEW-CO2 dry air mole fraction read from ',trim(OceanInputDir),'/',trim(GlobalviewFile)
       print *,'[emission_co2_ocean_init]    mapping from GLOBALVIEW MBL latitudes to 1x1 flux grid to be read from ',trim(OceanInputDir),'/',trim(GVGridMapFile)
       print *,'[emission_co2_ocean_init]    mapping from TM5 surface pressure grids to 1x1 flux grid to be read from ',trim(OceanInputDir),'/',trim(OceanGridMapFile)

    case ( 'takahashi' )
       !write (gol,'("[emission_co2_ocean_init]	pCO2 (air) from Takahashi year-2000 data..")')
       !call goPr
       !write (gol,'("[emission_co2_ocean_init]	   mapping from Takahashi 4x5 grid to 1x1 flux grid to be read from """,a,"/",a,""".")') trim(OceanInputDir), trim(OceanGridMapFile)
       !call goPr
       print *,'[emission_co2_ocean_init] pCO2 (air) from Takahashi year-2000 data..'
       print *,'[emission_co2_ocean_init]    mapping from Takahashi 4x5 grid to 1x1 flux grid to be read from ',trim(OceanInputDir),'/',trim(OceanGridMapFile)
 
    case default
       !write (gol,'("unsupported ocean.pco2.air value in rc file : ",a)') trim(ocean_pco2_air);call goPr
       !write (gol,'("supported strings are ""interactive"", ""globalview"", and ""takahashi""")');call goPr
       print *,'unsupported ocean.pco2.air value in file : trim(ocean_pco2_air)'
       print *,'supported strings are: interactive, globalview, and takahashi'
       !status=1
       !IF_NOTOK_RETURN(status=1)
       CALL ERROR_STOP('supported strings are interactive, globalview and takahashi',rname)
    end select

    !write (gol,'("[emission_co2_ocean_init]	   solubility and Schmidt no. to be read from """,a,"/",a,""".")') trim(OceanInputDir), trim(OceanSolFile)
    !call goPr

    print *,'[emission_co2_ocean_init]       solubility and Schmidt no. to be read from ',trim(OceanInputDir),'/',trim(OceanSolFile)

    !-- Not sure about this, AES
    !call Set(   u10m_dat(iglbsfc), status, used=.true. )
    !IF_NOTOK_RETURN(status=1)
    !call Set(   v10m_dat(iglbsfc), status, used=.true. )
    !IF_NOTOK_RETURN(status=1)
    !call Set(   ci_dat(iglbsfc), status, used=.true. )
    !IF_NOTOK_RETURN(status=1)
    !call Set(   m_dat(iglbsfc), status, used=.true. )
    !IF_NOTOK_RETURN(status=1)

  end subroutine emission_co2_ocean_init

  subroutine emission_co2_ocean_calc(State_Met,State_Chm,parms,co2_ocn,flux_updated,status)

    use ct_go_print_mod,            only : nread
    use ct_go_print_mod,		only : nlon360, nlat180
    use time_mod,               only : get_nymd, get_tau, get_taub
    use time_mod,               only : get_month, ITS_A_NEW_MONTH
    !use                         only : newmonth, iglbsfc
    !use Dims,			only : sec_month,sec_year
    !use emission_data,		only : flux_to_gridbox_per_month
    !use MDF,			only : MDF_NETCDF
    !use ReadFromFile_MDF,	only : ReadFromFile
    !use Meteo,			only : u10m_dat, v10m_dat, ci_dat, m_dat, sp_dat
    !use go_string,		only : goNum2str
    use ct_go_print_mod,        only : xmco2, xmair, xmc !,nmembersloc
    use enkf_mod,               only : enumber
    use ct_go_print_mod,        only : ico2tot,ico2bg, ico2ff, ico2bio, ico2ocean, ico2fires
    use ct_go_print_mod,        only : month_length
    !use datetime,		only : idate2ddate,tau2date
    !use global_data,		only : mass_dat
    USE TIME_MOD,               ONLY : ITS_A_LEAPYEAR, GET_YEAR, GET_MONTH
    USE TIME_MOD,               ONLY : GET_DAY_OF_YEAR,GET_HOUR,GET_ELAPSED_MIN
    USE GIGC_State_Met_Mod,     ONLY : MetState
    USE GIGC_State_Chm_Mod,     ONLY : ChmState
    USE REGRID_1x1_MOD,         ONLY : RETURN_A_GEN_1x1

    implicit none

    character(len=*), parameter			:: rname = mname//'/emission_co2_ocean_calc'

    TYPE(MetState)           :: State_Met   ! Meteorology State object
    TYPE(ChmState)           :: State_Chm   ! Chemistry/tracer State object

    type(field_1x1), dimension(:), intent(out)	:: co2_ocn
    type(enkf_parms), dimension(:), intent(in)	:: parms
    logical, intent(out)			:: flux_updated

    integer, intent(out)			:: status
    integer					:: i,j,ii,jj,imember
    real					:: month2year
    real*4,dimension(nlon360,nlat180,1)		:: buff3d
    real*4,dimension(nlon360,nlat180)		:: k92,icefrac,pco2air,pco2sw
    !real,pointer,dimension(:,:)			:: u10m,v10m,ci
    integer, dimension(6)			:: idatec
    integer					:: month, gr, gi, gj, gvl, gvj
    real					:: dec_date
    integer                                     :: itau,idate
    integer,dimension(12)           :: mlen ! days per month (current year)
    integer                                     ::  YR
    logical                                     :: LP
    integer                                     :: NCID,ret

!    real :: this_mf

    ! Pointers to MetState objects
      REAL*8, POINTER :: MET_U_10M(:,:)
      REAL*8, POINTER :: MET_V_10M(:,:)
      REAL*8, POINTER :: MET_OICE(:,:)
      REAL*8, POINTER :: MET_AIRMASS(:,:,:)
      REAL*8, POINTER :: MET_PSC2(:,:)  !NEEDS CONVERSION TO Pa, * 0.01

     ! Pointer to ChmState objects
      REAL*8, POINTER :: STT(:,:,:,:)

    !-- Points to met fields
      MET_U_10M             =>      State_Met%U10M
      MET_V_10M             =>      State_Met%V10M
      MET_AIRMASS           =>      State_Met%AD
      MET_OICE              =>      State_Met%FRSEAICE
      MET_PSC2              =>      State_Met%PSC2 

    !-- Points to chem fields
      STT                   =>      State_Chm%Tracers

    flux_updated = .false.

    itau = GET_ELAPSED_MIN() * 60
    !print *,'itau:',itau
    !print *,'nread:',nread
    !print *,'itau%%nread:',mod(itau,nread)
    if(mod(itau,nread) /= 0) return  ! only every nread hours
    !print *,'in mod loop...'
    flux_updated = .true.

    !call tau2date(itau,idatec)
    !month=idatec()
    YR = GET_YEAR()
    LP = ITS_A_LEAPYEAR( YR )
    month = get_month()
    !day = get_day_of_year()
    !hr = get_hour()

    IF(LP) THEN
       dec_date = yr + REAL(GET_DAY_OF_YEAR())/366. + REAL(GET_HOUR())/(366. * 24.)
    ELSE
       dec_date = yr + REAL(GET_DAY_OF_YEAR())/365. + REAL(GET_HOUR())/(365. * 24.)
    ENDIF
    PRINT *,'dec_date:',dec_date

    !dec_date=idate2ddate(idatec)

    if(trim(ocean_pco2_air) == 'globalview') then
       print *,'in globalview loop'
       if( .not. allocated(gv_mbl)) then
          call get_gv_mbl(trim(OceanInputDir)//'/'//trim(GlobalviewFile),status)
          !IF_NOTOK_RETURN(status=1)
       end if

       ! get bracketing dates from gv_dd
       gvl = 1+floor((dec_date-1979.0)/(1./48.))

       if((gvl .lt. lbound(gv_mbl,2)) .or. (gvl .gt. ubound(gv_mbl,2))) then
          !write (gol,'("Time index into gv_mbl,",i," is outside array bounds: ",2i)') gvl,lbound(gv_mbl,2),ubound(gv_mbl,2);call goPr
          !status=1
          !IF_NOTOK_RETURN(status=1)
          PRINT *,'Time index into gv_mbl,',gvl,' is outside array bounds:',lbound(gv_mbl,2),ubound(gv_mbl,2)
          CALL ERROR_STOP('',rname) 
      end if
    end if

    ! get necessary fields

    if(ITS_A_NEW_MONTH()) then
       print *,'its a new month'
       !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanPCO2File), &
       !     ftype=MDF_NETCDF,vname='PCO2_SW', buffer=buff3d, status=status)
       !IF_NOTOK_RETURN(status=1)

       !pco2sw_in = buff3d(:,:,month)
       !print *,'opening:',trim(OceanInputDir)//'/'//trim(OceanPCO2File),' month:',month
       CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanPCO2File))
       CALL NcRd(pco2sw_in,NCID,'PCO2_SW', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
       CALL Nccl (NCID)
       !print *,'check1'
       !pco2sw_in = buff3d(:,:,1)
       !print *,'check2'
       !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanSolFile), &
       !     ftype=MDF_NETCDF,vname='sol', buffer=buff3d, status=status)
       !IF_NOTOK_RETURN(status=1)

       !solubility = buff3d(:,:,month)

       CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanSolFile))
       CALL NcRd(solubility,NCID,'sol', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
       !CALL Nccl (NCID)
       !print *,'here'
       !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanSolFile), &
       !     ftype=MDF_NETCDF,vname='Sc', buffer=buff3d, status=status)
       !IF_NOTOK_RETURN(status=1)

       !schmidtno = buff3d(:,:,month)

       !CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanSolFile))
       CALL NcRd(schmidtno,NCID,'Sc', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
       !CALL Nccl (NCID)

       !call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanSolFile), &
       !     ftype=MDF_NETCDF,vname='SVP', buffer=buff3d, status=status)
       !IF_NOTOK_RETURN(status=1)

       !svp = buff3d(:,:,month)

       !CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanSolFile))
       CALL NcRd(svp,NCID,'SVP', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
       !CALL Nccl (NCID)

       !where(svp(:,:) .lt.0.0) svp=0.01

       if(trim(ocean_pco2_air) == 'takahashi') then
       !   call ReadFromFile(fname=trim(OceanInputDir)//'/'//trim(OceanSolFile), &
       !        ftype=MDF_NETCDF,vname='PCO2_AIR', buffer=buff3d, status=status)
       !   IF_NOTOK_RETURN(status=1)

       !   pco2air = buff3d(:,:,month)

       !CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(OceanSolFile))
       CALL NcRd(pco2air,NCID,'PCO2_AIR', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
       endif 
        CALL Nccl (NCID)
     endif

     !u10m => u10m_dat(iglbsfc)%data1(:,:,1)
     !v10m => v10m_dat(iglbsfc)%data1(:,:,1)
     !ci => ci_dat(iglbsfc)%data1(:,:,1)
     ! following Takahashi et al. (2009), maximum impermeability of ice cover is 0.95
     !icefrac = ci
     !print *,'met:',MET_U_10M(1,1),MET_V_10M(1,1),MET_OICE(1,1), sum(solubility),sum(schmidtno)
      do i = 1,360
          do j = 1,180
                gr=grid_region(i,j)
                gi=grid_i(i,j)
                gj=grid_j(i,j)
                !gr=1
                !gi=1
                !gj=1
                icefrac(i,j) = REAL(MET_OICE(gi,gj),4)
                IF(icefrac(i,j) > 0.95) THEN
                       icefrac(i,j)=0.95
                END IF

                k92(i,j)=0.0
                if((solubility(i,j).ge.0.0).and.(schmidtno(i,j).ge.0.0)) THEN
                      k92(i,j)=(1.0-icefrac(gi,gj))*solubility(i,j)*(0.31*(REAL(MET_U_10M(gi,gj),4)**2+REAL(MET_V_10M(gi,gj),4)**2)/sqrt(schmidtno(i,j)/660.0))*1.0e-2/3600.0
                END IF
        end do
    end do
    !icefrac = REAL(MET_OICE,4)
    !print*,'check'
    !where(MET_OICE(:,:) > 0.95) icefrac=0.95
    !print *,'here1.5'
    ! make gas transfer velocity in m/s according to Wanninkhof 1992 quadratic parameterization
    !k92=0.0
    !where ((solubility(:,:).ge.0.0).and.(schmidtno(:,:).ge.0.0)) &     ! masked to WOA09 ocean points
    !     k92=(1.0-icefrac)*solubility*(0.31*(MET_U_10M**2+MET_V_10M**2)/sqrt(schmidtno/660.0))*1.0e-2/3600.0
    !print *,'here1.75'
    ! Convert pco2sw_in to pco2sw.  The _in version needs to be retained since it is
    ! only read at the beginning of each month.  Here we convert to atm from input microatm
    ! and optionally add temporal trend.
    if(apply_sw_trend) then
       pco2sw = 1.0e-6*(pco2sw_in + 1.5*(dec_date - 2000.0))  ! trend is 1.5 microatm yr-1
       !       write (gol,'("[emission_co2_ocn_calc]	at dd", f12.4,": adjusted pco2sw by ",f12.4," uatm.")') dec_date,1.5*(dec_date-2000.0)
       !       call goPr
    else
       pco2sw = 1.0e-6*pco2sw_in
    end if
    !print *,'here2'
    do imember = 1,enumber

       ! compute pco2air
       !   need whole co2 here, so forward and inverse do different things
       !   pco2air is computed from the dry-air CO2 mole fraction.  It is then
       !   converted to CO2 partial pressure by multiplying by surface pressure,
       !   less the presumed saturation vapor pressure of water.
       !
       !   Units:
       !     convert from Pa to atm by dividing by 101325.
       !
       ! Note that SVP is read in as p(h2o)/p(tot)
       !
       ! Use nmembersloc to figure out if this is a forward run with CO2 components,
       ! or an inverse run with an ensemble of whole-CO2 tracers.  TODO: not the
       ! best method to make this decision.

       select case ( ocean_pco2_air )

       case ( 'interactive' )

          do i=1,360
             do j=1,180

                gr=grid_region(i,j)
                gi=grid_i(i,j)
                gj=grid_j(i,j)

                !gr=1
                !gi=1
                !gj=1

                if(enumber .eq. 1) then

                   ! forward:  total co2 = bio + ocn + fires + ff - 3*bg
                   ! 0.01 is for units change between GCHEM and TM5 SP (hPa vs Pa)
                   !pco2air(i,j) = 0.01*MET_PSC2(gi,gj)*(1-svp(i,j))*(xmair/xmco2)*(1/101325.0)* &
                   !     (STT(gi,gj,1,ico2ff) + &
                   !     STT(gi,gj,1,ico2ocean) + &
                   !     STT(gi,gj,1,ico2bio) + &
                   !     STT(gi,gj,1,ico2fires) - &
                   !     3*STT(gi,gj,1,ico2bg)) &
                   !     /MET_AIRMASS(gi,gj,1)

                    pco2air(i,j) = 100*MET_PSC2(gi,gj)*(1-svp(i,j))*(xmair/xmco2)*(1/101325.0)* &
                        (STT(gi,gj,1,ico2tot))/MET_AIRMASS(gi,gj,1)
!                    print *,'ppm:',pco2air(i,j)
!                    print *,'mass:',STT(gi,gj,1,ico2tot)/MET_AIRMASS(gi,gj,1)
!                    print *,'met_psc2:',MET_PSC2(gi,gj)
!                    print *,'svp:',svp(i,j)
!                   if((i .eq. 341) .and. (j .eq.24)) then
!                      this_mf= 1.e6*(mass_dat(gr)%rm_t(gi,gj,1,ico2ff) + &
!                        mass_dat(gr)%rm_t(gi,gj,1,ico2ocean) + &
!                        mass_dat(gr)%rm_t(gi,gj,1,ico2bio) + &
!                        mass_dat(gr)%rm_t(gi,gj,1,ico2fires) - &
!                        3*mass_dat(gr)%rm_t(gi,gj,1,ico2bg)) &
!                        /m_dat(gr)%data(gi,gj,1)
!                   endif
                else
                   ! inverse: full co2 exists in each tracer
                   pco2air(i,j) = 100*MET_PSC2(gi,gj)*(1-svp(i,j))*(xmair/xmco2)*(1/101325.0)*STT(gi,gj,1,1)/MET_AIRMASS(gi,gj,1)
                end if
             end do
          end do
        !print *,'nuthamax:',maxval(pco2air)
        !print *,'pco2air:',pco2air(1:30,1:30)
       case ( 'globalview' )

          do i=1,360
             do j=1,180

                gr=grid_region(i,j)
                gi=grid_i(i,j)
                gj=grid_j(i,j)

                !gr=1
                !gi=1
                !gj=1

                gvj=gv_j(j)
                pco2air(i,j) = 100*MET_PSC2(gi,gj)*(1-svp(i,j))*(1/101325.0)*1.0e-6*gv_mbl(gvj,gvl)

!                   if((i .eq. 341) .and. (j .eq.24)) then
!                      this_mf= gv_mbl(gvj,gvl)
!                   endif

             end do
          end do

       case ( 'takahashi' )

          ! no op

       case default
          !write (gol,'("unsupported ocean.pco2.air value in rc file : ",a)') trim(ocean_pco2_air);call goPr
          !write (gol,'("supported strings are ""interactive"", ""globalview"", and ""takahashi""")');call goPr
          print *,'unsupported ocean.pco2.air value in rc file : ', trim(ocean_pco2_air)
          print *,'supported strings are interactive, globalview, and takahashi'
          !status=1
          !IF_NOTOK_RETURN(status=1)

       end select

       ! make fluxes of CO2 in mol m-2 s-1
       !
       ! masking to ocean grid points only:  k92 is already WOA09-masked via
       ! solubility and Schmidt number; Takahashi pco2sw requires an additional
       ! mask

       co2_ocn(1)%surf = 0.0

       !print *,'sum(co2_ocn):',sum(co2_ocn(1)%surf)

       where(pco2sw(:,:) .gt.0.0) &
          co2_ocn(1)%surf = k92*(pco2sw-pco2air) * parms(1)%surf    ! apply parameters by member (forward has nmembersloc=1)
       !print *,'max(parms(1)):',maxval(parms(1)%surf)
       !print *,'max(k92):',maxval(k92)
       !print *,'max(pco2sw):',maxval(pco2sw)
       !print *,'max(pco2air):',maxval(pco2air)
       !print *,'min(diff):',minval(pco2sw-pco2air)
       !print *,'min(co2_ocn):',minval(co2_ocn(1)%surf)
 !      write (gol,'("at 341,24 we have pco2 sw/air; mf =",3f9.3,".")') 1.e6*pco2sw(341,24),1.e6*pco2air(341,24),this_mf
 !      call goPr

       ! co2_ocn now in mol m-2 s-1

       co2_ocn(2)%surf=co2_ocn(1)%surf*xmc/1.0e3 ! from mol m-2 s-1 to kgC m-2 s-1

       !print *,'sum(co2_ocn):',sum(co2_ocn(2)%surf)

       !call flux_to_gridbox_per_month(co2_ocn(1)%surf) ! convert from kgC m-2 s-1 to kgC 1x1_gridbox-1 month-1

        do j=1,nlat180
           do i=1,nlon360
              co2_ocn(2)%surf(i,j) = co2_ocn(2)%surf(i,j) * RETURN_A_GEN_1x1(j) * REAL(month_length(get_month(),LP))*24.*3600. 
           end do
        end do 

       !print *,'sum(co2_ocn):',sum(co2_ocn(2)%surf)

        mlen(1)=31
        mlen(2)=28
        IF(LP) mlen(2)=29
        mlen(3)=31
        mlen(4)=30
        mlen(5)=31
        mlen(6)=30
        mlen(7)=31
        mlen(8)=31
        mlen(9)=30
        mlen(10)=31
        mlen(11)=30
        mlen(12)=31 ! only for regular year

       IF (LP) THEN
          month2year = (86400 * 366 ) / (86400 * REAL(mlen(month)) )
       ELSE
          month2year = ( 86400 * 365 )/ (86400 * REAL(mlen(month) ) )
       ENDIF

       !month2year = sec_year/sec_month
       !write (gol,'("[emission_co2_ocn_calc]	member ",i4,": oceanic emission global integral ",f12.4," PgC/yr.")') imember,sum(co2_ocn(imember)%surf)*1.0e-12*month2year
       !call goPr
       !print *,'[emission_co2_ocn_calc] ensemble member ',enumber,': oceanic emission global integral ',sum(co2_ocn(2)%surf)*1.0e-12*month2year,' PgC/yr.'

    end do

    !nullify(u10m)
    !nullify(v10m)
    !nullify(ci)
    NULLIFY( MET_U_10M )
    NULLIFY( MET_V_10M )
    NULLIFY(MET_AIRMASS)
    NULLIFY(MET_OICE)
    NULLIFY(MET_PSC2)

    status=0

  end subroutine emission_co2_ocean_calc

  subroutine emission_co2_ocean_calc_oif(State_Met,State_Chm,parms,co2_ocn,flux_updated,status)

    !use Dims,			only : nlon360, nlat180, idate, itau, nread, newmonth, iglbsfc
    !use Dims,			only : sec_month,sec_year
    !use emission_data,		only : flux_to_gridbox_per_month
    !use MDF,			only : MDF_NETCDF
    !use ReadFromFile_MDF,	only : ReadFromFile
    !use chem_param,		only : xmc
    !use Meteo,			only : u10m_dat, v10m_dat, ci_dat
    !use go_string,		only: goNum2str
    use ct_go_print_mod,            only : nread
    use ct_go_print_mod,                only : nlon360, nlat180
    use time_mod,               only : get_nymd, get_tau, get_taub
    use time_mod,               only : get_month, ITS_A_NEW_MONTH
    !use go_string,             only : goNum2str
    use ct_go_print_mod,        only : xmco2, xmair, xmc !,nmembersloc
    use enkf_mod,               only : enumber
    use ct_go_print_mod,        only : ico2tot,ico2bg, ico2ff, ico2bio, ico2ocean, ico2fires
    use ct_go_print_mod,        only : month_length
    !use datetime,              only : idate2ddate,tau2date
    !use global_data,           only : mass_dat
    USE TIME_MOD,               ONLY : ITS_A_LEAPYEAR, GET_YEAR, GET_MONTH
    USE TIME_MOD,               ONLY : GET_DAY_OF_YEAR,GET_HOUR,GET_ELAPSED_MIN
    USE GIGC_State_Met_Mod,     ONLY : MetState
    USE GIGC_State_Chm_Mod,     ONLY : ChmState
    USE REGRID_1x1_MOD,         ONLY : RETURN_A_GEN_1x1

    implicit none

    character(len=*), parameter			:: rname = mname//'/emission_co2_ocean_calc_oif'

    TYPE(MetState)           :: State_Met   ! Meteorology State object
    TYPE(ChmState)           :: State_Chm   ! Chemistry/tracer State object

    type(field_1x1), dimension(:), intent(out)	:: co2_ocn
    type(enkf_parms), dimension(:), intent(in)	:: parms
    logical, intent(out)			:: flux_updated

    integer, intent(out)			:: status
    integer					:: i,j,ii,jj,imember
    real					:: month2year
    real*4,dimension(nlon360,nlat180)		:: ww,k92,schmidtno,solubility,deltapco2,flux
    !real,pointer,dimension(:,:)			:: u10m,v10m,ci
    integer, dimension(6)                       :: idatec
    integer                                     :: month, gr, gi, gj, gvl, gvj
    real                                        :: dec_date
    integer                                     :: itau,idate
    integer,dimension(12)           :: mlen ! days per month (current year)
    integer                                     ::  YR
    logical                                     :: LP
    integer                                     :: NCID,ret

    ! Pointers to MetState objects
      REAL*8, POINTER :: MET_U_10M(:,:)
      REAL*8, POINTER :: MET_V_10M(:,:)
      REAL*8, POINTER :: MET_OICE(:,:)
      REAL*8, POINTER :: MET_AIRMASS(:,:,:)
      REAL*8, POINTER :: MET_PSC2(:,:)  !NEEDS CONVERSION TO Pa, * 0.01

     ! Pointer to ChmState objects
      REAL*8, POINTER :: STT(:,:,:,:)

    !-- Points to met fields
      MET_U_10M             =>      State_Met%U10M
      MET_V_10M             =>      State_Met%V10M
      MET_AIRMASS           =>      State_Met%AD
      MET_OICE              =>      State_Met%FRSEAICE
      MET_PSC2              =>      State_Met%PSC2

    !-- Points to chem fields
      STT                   =>      State_Chm%Tracers

    flux_updated = .false.

   itau = GET_ELAPSED_MIN() * 60

    if(mod(itau,nread) /= 0) return  ! only every nread hours

    flux_updated = .true.

    !call tau2date(itau,idatec)
    !month=idatec()
    YR = GET_YEAR()
    LP = ITS_A_LEAPYEAR( YR )
    month = get_month()
    !day = get_day_of_year()
    !hr = get_hour()

    ! get necessary fields

    !if(newmonth) then
    if(ITS_A_NEW_MONTH()) then

      !-- NEEDS FIXING 
      !CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(deltapco2_prefix)//'.'//trim(goNum2str(idate(1)))//'.'//trim(goNum2str(idate(2),'(i2.2)'))//'.nc')
      CALL Ncop_Rd (NCID, trim(OceanInputDir)//'/'//trim(deltapco2_prefix)//'junk.nc')
      CALL NcRd(deltapco2,NCID,'DELTAPCO2', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
      CALL NcRd(solubility,NCID,'SOLUBILITY', (/ 1,1,month /), (/ nlon360,nlat180,1 /))
      CALL NcRd(schmidtno,NCID,'SCHMIDTNO', (/ 1,1,month /), (/ nlon360,nlat180,1 /))      

       CALL Nccl (NCID)

    endif

    !u10m => u10m_dat(iglbsfc)%data1(:,:,1)
    !v10m => v10m_dat(iglbsfc)%data1(:,:,1)
    !ci => ci_dat(iglbsfc)%data1(:,:,1)

    ! make diffusion in m/s according to Wanninkhof92 param on 1x1 degree

    k92=0.0
    ww=REAL(MET_U_10M,4)**2+REAL(MET_V_10M,4)**2  ! note: sqrt has been taken away since we square it in the next formula!
    where(schmidtno(:,:).ge.0) k92=(0.31*(ww)/sqrt(schmidtno/660.0))*1.e-2/3600.   ! only over oceans

    ! make fluxes of co2 in mol/m2/s with tm5 meteo.  1e-6 to change deltapco2 from uatm to atm.
    flux=0
    where((solubility(:,:).ge.0).and.(deltapco2(:,:).gt.-1e10)) flux=k92*solubility*deltapco2*1e-6
    flux=flux*(1.0-MET_OICE) ! multiply out the sea-ice for this day

    do imember = 1,enumber
       co2_ocn(imember)%surf = flux * parms(imember)%surf    ! apply parameters by member (forward has nmembersloc=1)

       ! co2_ocn now in mol m-2 s-1

       co2_ocn(2)%surf=co2_ocn(1)%surf*xmc/1.0e3 ! from mol m-2 s-1 to kgC m-2 s-1

       !call flux_to_gridbox_per_month(co2_ocn(imember)%surf) ! convert from kgC m-2 s-1 to kgC 1x1_gridbox-1 month-1

        do j=1,nlat180
           do i=1,nlon360
              co2_ocn(2)%surf(i,j) = co2_ocn(2)%surf(i,j) * RETURN_A_GEN_1x1(j) * REAL(month_length(get_month(),LP))*24.*3600.
           end do
        end do

        mlen(1)=31
        mlen(2)=28
        IF(LP) mlen(2)=29
        mlen(3)=31
        mlen(4)=30
        mlen(5)=31
        mlen(6)=30
        mlen(7)=31
        mlen(8)=31
        mlen(9)=30
        mlen(10)=31
        mlen(11)=30
        mlen(12)=31 ! only for regular year

       IF (LP) THEN
          month2year = (86400 * 366 ) / (86400 * REAL(mlen(month)) )
       ELSE
          month2year = ( 86400 * 365 )/ (86400 * REAL(mlen(month) ) )
       ENDIF

       !write (gol,'("[emission_co2_ocn_calc]	member ",i4,": oceanic emission global integral ",f12.4," PgC/yr.")') imember,sum(co2_ocn(imember)%surf)*1.0e-12*month2year
       !call goPr

       !print *,'[emission_co2_ocn_calc] ensemble member ',enumber,': oceanic emission global integral ',sum(co2_ocn(2)%surf)*1.0e-12*month2year,' PgC/yr.'
    end do

    !nullify(u10m)
    !nullify(v10m)
    !nullify(ci)

    NULLIFY( MET_U_10M )
    NULLIFY( MET_V_10M )
    NULLIFY(MET_AIRMASS)
    NULLIFY(MET_OICE)
    NULLIFY(MET_PSC2)

    status=0

  end subroutine emission_co2_ocean_calc_oif

  subroutine get_gv_mbl(fname,status)

    !use MDF,			only : MDF_open, MDF_Inq_VarID, MDF_Get_Var, MDF_Close, MDF_NETCDF, MDF_READ, MDF_Inq_DimID, MDF_Inquire_Dimension
    !use dims,			only : idatei,idatee
    use ct_go_print_mod,        only : idate2ddate
    use time_mod,               only : GET_NYMDe,GET_NHMSe,GET_NYMDb,GET_NHMSb
    use time_mod,               only : YMD_EXTRACT

    character(len=*),intent(in) 	:: fname
    integer,intent(out)			:: status

    integer 				:: hid, vid, dimid
    integer 				:: tstart=-1,tcount=-1,ttotal=-1,itime
    real,dimension(:),allocatable	:: ddate
    real 				:: ddatei,ddatee
    integer                             :: Y1,Y2,MM1,MM2,D1,D2,H1,H2,M1,M2,S1,S2
    integer                             :: NHMSb, NYMDb, NHMSe, NYMDe
    integer                             :: NCID

    character(len=*), parameter	:: rname = mname//'/get_gv_mbl'

    ! based on start and end date of simulation, compute tstart and tcount

    !call MDF_open(fname,ftype=MDF_NETCDF,mode=MDF_READ,hid=hid,status=status)
    CALL Ncop_Rd (NCID, trim(fname))
    !IF_NOTOK_RETURN(status=1)

    !call MDF_Inq_DimID(hid=hid,name='time',dimid=dimid,status=status)
    CALL Ncget_Dimlen (NCID, 'time', ttotal)
    !IF_NOTOK_RETURN(status=1)
    !print *,'total:',ttotal
    !call MDF_Inquire_Dimension(hid=hid,dimid=dimid,status=status,length=ttotal)
    !IF_NOTOK_RETURN(status=1)

    allocate(ddate(ttotal))

    !call MDF_Inq_VarID(hid,name='decimal_date',varid=vid,status=status)
    CALL NcRd(ddate,NCID,'decimal_date', (/ 1 /), (/ ttotal /)) 
    !IF_NOTOK_RETURN(status=1)
    !print *,'ddate:',ddate 
    !call MDF_Get_Var(hid,vid,ddate,status)
    !IF_NOTOK_RETURN(status=1)

    NYMDb = GET_NYMDb()
    NHMSb = GET_NHMSb()
    NYMDe = GET_NYMDe()
    NHMSe = GET_NHMSe()

    CALL YMD_EXTRACT( NYMDb, Y1, MM1, D1)
    CALL YMD_EXTRACT( NHMSb, H1, M1, S1)

    CALL YMD_EXTRACT( NYMDe, Y2, MM2, D2)
    CALL YMD_EXTRACT( NHMSe, H2, M2, S2)

    ddatee=idate2ddate((/Y2,MM2,D2,H2,M2,S2/))
    ddatei=idate2ddate((/Y1,MM1,D1,H1,M1,S1/))

    do itime=ttotal,1,-1
       if(ddate(itime) .le. ddatei) then
          tstart=itime
          exit
       end if
    end do

    do itime=1,ttotal
       if(ddate(itime) .ge. ddatee) then
          tcount=itime-tstart+1
          exit
       end if
    end do

    if((tstart .lt. 1) .or. (tstart .ge. ttotal) .or. (tcount .lt. 1) .or. (tcount .gt. ttotal)) then
       !write (gol,'("Cannot find appropriate time interval in GLOBALVIEW file bracketing the simulation start and end times.")');call goPr
       !write (gol,'("  Simulation start/end decimal dates:",2f12.6)') ddatei,ddatee;call goPr
       !write (gol,'("  GLOBALVIEW start/end decimal dates:",2f12.6)') ddate(1),ddate(ttotal);call goPr
       !write (gol,'("  tstart, tcount, ttotal:",3i5)') tstart,tcount,ttotal;call goPr
       print *,'Cannot find appropriate time interval in GLOBALVIEW file bracketing the simulation start and end times.'
       print *,' Simulation start/end decimal dates:',ddatei,ddatee
       print *,'  GLOBALVIEW start/end decimal dates:',ddate(1),ddate(ttotal)
       print *,'  tstart, tcount, ttotal:', tstart,tcount,ttotal
       !status=1
       !IF_NOTOK_RETURN(status=1)
    end if

    allocate(gv_mbl(41,tstart:(tstart+tcount-1)))
    allocate(gv_dd(tcount))

!    write (gol,'("Reading from """,a,""".")') fname;call goPr
!    write (gol,'("  Simulation start/end decimal dates:",2f12.6)') ddatei,ddatee;call goPr
!    write (gol,'("  GLOBALVIEW start/end decimal dates:",2f12.6)') ddate(1),ddate(ttotal);call goPr

    !call MDF_Inq_VarID(hid=hid,name='mbl',varid=vid,status=status)
    !IF_NOTOK_RETURN(status=1)

    CALL NcRd(gv_mbl,NCID,'mbl', (/ 1,tstart /), (/ 41,tcount /))

    !call MDF_Get_Var(hid=hid,varid=vid,values=gv_mbl,start=(/1,tstart/),count=(/41,tcount/),status=status)
    !IF_NOTOK_RETURN(status=1)

    !call MDF_Close(hid,status)
    !IF_NOTOK_RETURN(status=1)

    CALL Nccl (NCID)

    status = 0

  end subroutine get_gv_mbl

  subroutine emission_co2_ocean_done(status)

    use ct_go_print_mod,		only: nregions

    implicit none

    integer, intent(out)			:: status

    character(len=*), parameter	:: rname = mname//'/emission_co2_ocean_done'

    print *,'cleaning up CT ocean flux objects....'

    if(allocated(grid_region)) then
       deallocate(grid_region)
    end if

    if(allocated(grid_i)) then
       deallocate(grid_i)
    end if

    if(allocated(grid_j)) then
       deallocate(grid_j)
    end if

    if(allocated(gv_mbl)) then
       deallocate(gv_mbl)
    end if

    if(allocated(gv_dd)) then
       deallocate(gv_dd)
    end if

    if(allocated(gv_j)) then
       deallocate(gv_j)
    end if

    if(allocated(co2_ocn)) then
       deallocate(co2_ocn)
    endif 

    status=0

  end subroutine emission_co2_ocean_done

end module emission_co2_ocean
