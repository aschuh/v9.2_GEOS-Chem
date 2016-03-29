!
!   Filename: ncFormatter.f90
!   Purpose: This module contains functions which enable the user to create a
!            CF compliant NetCDF file for 2 and 3 dimensional spatial data.
!            
!   Author: RMcKeown
!   History: Initial version for files containing a single variable, latitude, 
!            longitude, vertical levels (optional), and time.
!

module ncCooardsFormat
    use netcdf
    implicit none

    ! Handles for the file, dimensions, and variables
    integer, save :: fileNCID
    integer, save :: latDimID, lonDimID, levDimID, timeDimID
    integer, save :: latID, lonID, levID, timeID, varID

    ! Handles for the file, dimensions, and variables
    integer, save :: stationFileNCID, gosatOutputNCID, ctOutputNCID
    integer, save :: sLatID, sLonID, sLevID, sTimeID, sVarID,sNameID
    integer, save :: gLatID, gLonID, gLevID, gTimeID, gVarID,gNameID
    integer, save :: gVarID1,gVarID2,gVarID3

    ! Units required for CF compliance
    character (len = *), parameter :: LAT_UNITS = "degrees_north"
    character (len = *), parameter :: LON_UNITS = "degrees_east"
    character (len = *), parameter :: LEV_UNITS = "hPa"

    contains

subroutine ncGOSATOutputCreate(filename, nSoundings,xco2_cloudheight,xco2_latitude,xco2_longitude,sounding_id,xco2_obs,xco2_uncert,   &
                 xco2pbl_obs,xco2pbl_uncert,xco2ft_obs,xco2ft_uncert,varName1,varName2,varName3)
! subroutine ncGOSATOutputCreate(filename, nSoundings, sounding_id,varName)
        !   This subroutine creates a file to hold timeseries data for multiple stations.
        !   nStation:  number of stations
        !   nTime:     number of timesteps
        !   varName:   variable name

            character (len = *)      ::  filename
            character (len = *)      ::  varName1,varName2,varName3
        !   Need to add stationName to the argument list when ready
            integer                  :: soundingDimID, nSoundings
            integer                  :: sVarID2 ,sVarID3,sVarID4
            integer                  :: sVarID5 ,sVarID6,sVarID7,sVarID8,sVarID9,sVarID10,sVarID11
            integer*8,dimension(:)   :: sounding_id
            real*8, dimension(:)       :: xco2_latitude,xco2_longitude,xco2_obs,xco2_uncert
            real*8, dimension(:)       :: xco2pbl_obs,xco2pbl_uncert
            real*8, dimension(:)       :: xco2ft_obs,xco2ft_uncert
            real*8, dimension(:)       :: xco2_cloudheight

            !integer, parameter :: mode_flag = IOR(nf90_hdf5, nf90_classic_model) ! | nf90_clobber
            integer, parameter :: mode_flag = nf90_netcdf4
!           integer, parameter :: deflatelevel = 6

            !print *,'uncert:',xco2_uncert
            !print *,'sid:',sounding_id
            !print *,'xco2:',xco2_obs

            ! Create file
            call handleError(nf90_create(path=trim(filename), cmode=mode_flag, ncid=gosatOutputNCID))

            print *,'Opening/creating file ',trim(filename)

            ! Define Dimensions
            call handleError(nf90_def_dim(gosatOutputNCID, "soundings", nSoundings, soundingDimID))
            !call handleError(nf90_def_dim(gosatOutputNCID, "soundings", nf90_unlimited, soundingDimID))

            ! Define Variables
            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2_cloudheight", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID11))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="latitude", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID9))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="longitude", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID10))

            !call handleError(nf90_def_var(gosatOutputNCID, "sounding_id", nf90_uint64, soundingDimID, &
            !     sVarID2))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2_observed", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID3))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2_uncertainty", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID4))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2pbl_observed", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID5))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2pbl_uncertainty", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID6))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2freetrop_observed", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID7))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name="xco2freetrop_uncertainty", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID8))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name=varName1, xtype=nf90_double, dimids=soundingDimID, &
                 varid=gVarID1))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name=varName2, xtype=nf90_double, dimids=soundingDimID, &
                 varid=gVarID2))

            call handleError(nf90_def_var(ncid=gosatOutputNCID, name=varName3, xtype=nf90_double, dimids=soundingDimID, &
                 varid=gVarID3))

!                            deflate_level=deflatelevel))

            ! Define Attributes
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "long_name",                    &
            !                                  "ACOS Sounding ID"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "units", "MMDDYYhhmmss"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "comment", "from scan start time in UTC"))

            !call handleError(nf90_put_att(gosatOutputNCID, sVarID, "long_name",                    &
            !                                  "Optimized value of Retrieved XCO2"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID, "units", "ppm"))


            call handleError(nf90_enddef(gosatOutputNCID))

            !call handleError(nf90_put_var(gosatOutputNCID, sVarID2, sounding_id))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID3, xco2_obs))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID4, xco2_uncert))

            call handleError(nf90_put_var(gosatOutputNCID, sVarID5, xco2pbl_obs))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID6, xco2pbl_uncert))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID7, xco2ft_obs))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID8, xco2ft_uncert))

            call handleError(nf90_put_var(gosatOutputNCID, sVarID9, xco2_latitude))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID10, xco2_longitude))
            call handleError(nf90_put_var(gosatOutputNCID, sVarID11, xco2_cloudheight))

            ! Sync the file to make sure data is saved
            call handleError(nf90_sync(gosatOutputNCID))

        end subroutine ncGOSATOutputCreate

        subroutine WRITE_VARS_TO_FILE_NC( M, V, V2, V3 )

            integer             ::  M,localgVarID
            integer             ::  localgVarID2,localgVarID3
            real*8              ::  V, V2, V3
            
             
             !print *,'writing ', V,' to ',M,'th spot'
 
            flush(6)

            !V = REAL(V * 10**6)
            !V2 = REAL(V2 * 10**6)
            !V3 = REAL(V3 * 10**6)

            V  = V * 10**6
            V2 = V2 * 10**6
            V3 = V3 * 10**6

            call handleError(nf90_inq_varid(ncid=gosatOutputNCID, name="xco2_model",varid=localgVarID))
            call handleError(nf90_inq_varid(ncid=gosatOutputNCID, name="xco2pbl_model",varid=localgVarID2))
            call handleError(nf90_inq_varid(ncid=gosatOutputNCID, name="xco2freetrop_model",varid=localgVarID3))
            call handleError(nf90_put_var(ncid=gosatOutputNCID,varid=localgVarID , values=V,start=(/ M /) ) )
            call handleError(nf90_put_var(ncid=gosatOutputNCID,varid=localgVarID2 , values=V2,start=(/ M /) ) )
            call handleError(nf90_put_var(ncid=gosatOutputNCID,varid=localgVarID3 , values=V3,start=(/ M /) ) )

            ! Sync the file to make sure data is saved
            call handleError(nf90_sync(gosatOutputNCID))

        end subroutine WRITE_VARS_TO_FILE_NC

        subroutine ncGOSATOutputClose

          print *, "closing nc bling bling gosat output file", gosatOutputNCID

          call gosatCloseError(nf90_close(gosatOutputNCID))

        end subroutine ncGOSATOutputClose

 subroutine ncCTOutputCreate(filename, nSoundings,sounding_id,xco2_obs,xco2_uncert,   &
                 xco2pbl_obs,xco2pbl_uncert,xco2ft_obs,xco2ft_uncert,varName1,varName2,varName3)
! subroutine ncGOSATOutputCreate(filename, nSoundings, sounding_id,varName)
        !   This subroutine creates a file to hold timeseries data for multiple stations.
        !   nStation:  number of stations
        !   nTime:     number of timesteps
        !   varName:   variable name

            character (len = *)      ::  filename
            character (len = *)      ::  varName1,varName2,varName3
        !   Need to add stationName to the argument list when ready
            integer                  :: soundingDimID, nSoundings
            integer                  :: sVarID2 ,sVarID3,sVarID4
            integer                  :: sVarID5 ,sVarID6,sVarID7,sVarID8
            character,dimension(:,:)   :: sounding_id
            real*8, dimension(:)       :: xco2_obs,xco2_uncert
            real*8, dimension(:)       :: xco2pbl_obs,xco2pbl_uncert
            real*8, dimension(:)       :: xco2ft_obs,xco2ft_uncert

            !integer, parameter :: mode_flag = IOR(nf90_hdf5, nf90_classic_model) ! | nf90_clobber
            integer, parameter :: mode_flag = nf90_hdf5
!           integer, parameter :: deflatelevel = 6

            print *,'uncert:',xco2_uncert
            print *,'sid:',sounding_id
            print *,'xco2:',xco2_obs

            ! Create file
            call handleError(nf90_create(path=trim(filename), cmode=mode_flag, ncid=ctOutputNCID))

            print *,'Opening/creating file ',trim(filename)

            ! Define Dimensions
            call handleError(nf90_def_dim(ctOutputNCID, "soundings", nSoundings, soundingDimID))

            ! Define Variables

            call handleError(nf90_def_var(ctOutputNCID, "sounding_id", nf90_char, soundingDimID, &
                 sVarID2))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="xco2_observed", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID3))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="xco2_uncertainty", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID4))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="xco2pbl_observed", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID5))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="xco2pbl_uncertainty", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID6))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="xco2freetrop_observed", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID7))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="xco2freetrop_uncertainty", xtype=nf90_double, &
                 dimids=soundingDimID, varid=sVarID8))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name=varName1, xtype=nf90_double, dimids=soundingDimID, &
                 varid=gVarID1))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name=varName2, xtype=nf90_double, dimids=soundingDimID, &
                 varid=gVarID2))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name=varName3, xtype=nf90_double, dimids=soundingDimID, &
                 varid=gVarID3))

!                            deflate_level=deflatelevel))

            ! Define Attributes
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "long_name",                    &
            !                                  "ACOS Sounding ID"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "units", "MMDDYYhhmmss"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "comment", "from scan start time in UTC"))

            !call handleError(nf90_put_att(gosatOutputNCID, sVarID, "long_name",                    &
            !                                  "Optimized value of Retrieved XCO2"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID, "units", "ppm"))


            call handleError(nf90_enddef(ctOutputNCID))

            call handleError(nf90_put_var(ctOutputNCID, sVarID2, sounding_id))
            call handleError(nf90_put_var(ctOutputNCID, sVarID3, xco2_obs))
            call handleError(nf90_put_var(ctOutputNCID, sVarID4, xco2_uncert))

            call handleError(nf90_put_var(ctOutputNCID, sVarID5, xco2pbl_obs))
            call handleError(nf90_put_var(ctOutputNCID, sVarID6, xco2pbl_uncert))
            call handleError(nf90_put_var(ctOutputNCID, sVarID7, xco2ft_obs))
            call handleError(nf90_put_var(ctOutputNCID, sVarID8, xco2ft_uncert))

            ! Sync the file to make sure data is saved
            call handleError(nf90_sync(ctOutputNCID))

        end subroutine ncCTOutputCreate

        subroutine ncCTOutputClose

          print *, "closing nc CT output file", ctOutputNCID

          call handleError(nf90_close(ctOutputNCID))

        end subroutine ncCTOutputClose

        subroutine ncStationDataFileCreate(filename, nStation, sLon, sLat, sLev, &
                                   stationName, varName, varUnits, timeUnits)
        !   This subroutine creates a file to hold timeseries data for multiple stations.
        !   nStation:  number of stations
        !   nTime:     number of timesteps
        !   varName:   variable name

            character (len = *) ::  filename
            character (len = *) ::  varName
            character (len = *) ::  varUnits
            character (len = *) ::  timeUnits
        !   Need to add stationName to the argument list when ready
            integer, intent(in) :: nStation
            real, dimension(:) :: sLon, sLat
            integer, dimension(:) :: sLev
            integer :: stationDimID, sTimeDimID, nameLenDimID, sCount
            character (len = 255), dimension(nStation) :: stationName

            integer, dimension(2) :: start
            integer, parameter :: mode_flag = IOR(nf90_hdf5, nf90_classic_model) ! | nf90_clobber
!           integer, parameter :: deflatelevel = 6

            PRINT*,'entering nc file open'
            flush(6)
            PRINT*,'stationName(1):',STATIONNAME(1)

            ! Create file
            call handleError(nf90_create(path=trim(filename), cmode=mode_flag, ncid=stationFileNCID))

            ! Define Dimensions
            call handleError(nf90_def_dim(stationFileNCID, "time", nf90_unlimited, sTimeDimId))
            call handleError(nf90_def_dim(stationFileNCID, "station", nStation, stationDimID))
            call handleError(nf90_def_dim(stationFileNCID, "station name", 255, nameLenDimID))

            ! Define Variables
            call handleError(nf90_def_var(stationFileNCID, "lon", nf90_float, (/stationDimID/), sLonID))
            call handleError(nf90_def_var(stationFileNCID, "lat", nf90_float, (/stationDimID/), sLatID))
            call handleError(nf90_def_var(stationFileNCID, "level", nf90_int, (/stationDimID/), sLevID))
            call handleError(nf90_def_var(stationFileNCID, "station", nf90_char, (/nameLenDimID, stationDimID/), sNameID))
            call handleError(nf90_def_var(stationFileNCID, "time", nf90_float, (/sTimeDimID/), sTimeID))
            call handleError(nf90_def_var(stationFileNCID, varName, nf90_float, (/stationDimID, sTimeDimID/), &
                 sVarID))
!                            deflate_level=deflatelevel))

            ! Define Attributes
            call handleError(nf90_put_att(stationFileNCID, sTimeID, "sLong_name", "time"))
            call handleError(nf90_put_att(stationFileNCID, sTimeID, "units", timeUnits))
            call handleError(nf90_put_att(stationFileNCID, sLonID, "long_name", "longitude"))
            call handleError(nf90_put_att(stationFileNCID, sLonID, "units", LON_UNITS))
            call handleError(nf90_put_att(stationFileNCID, sLatID, "long_name", "latitude"))
            call handleError(nf90_put_att(stationFileNCID, sLevID, "sLong_name", "level"))
            call handleError(nf90_put_att(stationFileNCID, sLevID, "units", "model level"))
            call handleError(nf90_put_att(stationFileNCID, sNameID, "sLong_name", "station name"))
            call handleError(nf90_put_att(stationFileNCID, sLatID, "units", LAT_UNITS))
            call handleError(nf90_put_att(stationFileNCID, sVarID, "long_name", varName))
            call handleError(nf90_put_att(stationFileNCID, sVarID, "units", varUnits))

            call handleError(nf90_enddef(stationFileNCID))

            ! Write out coord vars
            call handleError(nf90_put_var(stationFileNCID, sLonID, sLon))
            call handleError(nf90_put_var(stationFileNCID, sLatID, sLat))
            call handleError(nf90_put_var(stationFileNCID, sLevID, sLev))
            
            do sCount = 1, nStation
                 start = (/1, sCount/)
                 call handleError( nf90_put_var( stationFileNCID, sNameID, trim(stationName(sCount)), start ) )
            enddo

            ! Sync the file to make sure data is saved
            !call handleError(nf90_sync(stationFileNCID))

        end subroutine ncStationDataFileCreate

        subroutine ncStationFileClose
          print *, "closing nc station file", stationFileNCID
          call handleError(nf90_close(stationFileNCID))
        end subroutine ncStationFileClose

        subroutine writeStationOneTime(var, time, nLevels, timeIDX, station)
            integer :: nLevels, timeIDX, station
            integer, dimension(2) :: vstart
            integer, dimension(1) :: tstart
            integer, save :: timeCheck = 0
            real :: time
            !real, dimension(nLevels) :: var
            real :: var

            if (timeIDX > timeCheck) then
                timeCheck = timeIDX
                tstart = (/timeCheck/)
                call handleError(nf90_put_var(stationFileNCID, sTimeID, time, start=tstart))
            endif

            vstart = (/station, timeIDX/)
            call handleError(nf90_put_var(stationFileNCID, sVarID, var, start=vstart))

            !call handleError(nf90_sync(stationFileNCID))
        end subroutine writeStationOneTime


        subroutine nc2dFileCreate(filename, nLat, nLon, nTime, varName, varUnits, timeUnits,varID,fileNCID)
            !-- IN
            character (len = *) ::  filename
            character (len = *) ::  varName
            character (len = *) ::  varUnits
            character (len = *) ::  timeUnits
            integer :: nLat, nLon, nTime 
            !-- OUT
            integer :: fileNCID
            integer :: varID
            !-- LOCAL
            integer, dimension (3) :: dimids

            print *,'in ncCoords'
            ! Create file
            call handleError(nf90_create(filename, nf90_clobber, fileNCID))
    
            ! Define Dimensions
            call handleError(nf90_def_dim(fileNCID, "lon", nLon, lonDimID))
            call handleError(nf90_def_dim(fileNCID, "lat", nLat, latDimID))
            call handleError(nf90_def_dim(fileNCID, "time", nTime, timeDimID))
            dimids = (/ lonID, latID, timeID /)

            ! Define Variables
            call handleError(nf90_def_var(fileNCID, "lon", nf90_float, lonDimID, lonID))
            call handleError(nf90_def_var(fileNCID, "lat", nf90_float, latDimID, latID))
            call handleError(nf90_def_var(fileNCID, "time", nf90_float, timeDimID, timeID))
            call handleError(nf90_def_var(fileNCID, varName, nf90_float, dimids, varID))

            ! Define Attributes
            call handleError(nf90_put_att(fileNCID, lonID, "long_name", "longitude"))
            call handleError(nf90_put_att(fileNCID, lonID, "units", LON_UNITS))
            call handleError(nf90_put_att(fileNCID, latID, "long_name", "latitude"))
            call handleError(nf90_put_att(fileNCID, latID, "units", LAT_UNITS))
            call handleError(nf90_put_att(fileNCID, lonID, "long_name", "time"))
            call handleError(nf90_put_att(fileNCID, latID, "units", timeUnits))

        end subroutine nc2dFileCreate

        subroutine nc3dFileCreate(filename, nLat, nLon, nLevels,nTime, varName, nVars,    & 
                                    varUnits, timeUnits, deflatelevel)
            character (len = *) ::  filename
            character (len = *) ::  varName(*)
            character (len = *) ::  varUnits
            character (len = *) ::  timeUnits
            integer :: i,nVars, nLat, nLon, nLevels, nTime, deflatelevel,mode_flag
            real   :: lats, lons, levels
            integer, dimension (4) :: dimids

            mode_flag = IOR(nf90_hdf5 , nf90_classic_model) ! | nf90_clobber
            
            call handleError(nf90_create(path=filename, cmode=mode_flag, ncid=fileNCID))
            !call  handleError(nf90_create(path=filename,cmode=nf90_clobber, ncid=fileNCID))

            call handleError(nf90_def_dim(fileNCID, "lon", nLon, lonDimID))
            call handleError(nf90_def_dim(fileNCID, "lat", nLat, latDimID))
            call handleError(nf90_def_dim(fileNCID, "Levels", nLevels, levDimID))
            !call handleError(nf90_def_dim(fileNCID, "time", nTime, timeDimID))
            call handleError(nf90_def_dim(fileNCID, "time", nf90_unlimited, timeDimId))
            dimids = (/ lonDimID, latDimID, levDimID, timeDimID /)

            !print *,'wrote dims'
            call handleError(nf90_def_var(fileNCID, "Levels", nf90_double, levDimID, levID))
            call handleError(nf90_def_var(fileNCID, "lon", nf90_double, lonDimID, lonID))
            call handleError(nf90_def_var(fileNCID, "lat", nf90_double, latDimID, latID))
            call handleError(nf90_def_var(fileNCID, "time", nf90_double, timeDimID, timeID))
            !print *,'defined dims'
         DO i = 1,nVars
            !print *,'nVars:',i,'varID:',varID
            !Testing doubleprecision
            !call handleError(nf90_def_var(fileNCID, varName(i), nf90_float, dimids, varID, & 
            !            deflate_level=deflatelevel)) !,chunksizes=(/ nLon, nLat, 1, 1 /)))
            call handleError(nf90_def_var(fileNCID, varName(i), nf90_double, dimids, varID, &
                        deflate_level=deflatelevel)) !,chunksizes=(/ nLon, nLat, 1, 1 /)))
            !print *,'defined ',varName(i)
            ! Define Attributes
            call handleError(nf90_put_att(fileNCID, varID, "long_name", "PPM CO2"))
            call handleError(nf90_put_att(fileNCID, varID, "units", "PPM"))
            call handleError(nf90_put_att(fileNCID, levID, "long_name", "Levels"))
            call handleError(nf90_put_att(fileNCID, levID, "units", LEV_UNITS))
            call handleError(nf90_put_att(fileNCID, lonID, "long_name", "longitude"))
            call handleError(nf90_put_att(fileNCID, lonID, "units", LON_UNITS))
            call handleError(nf90_put_att(fileNCID, latID, "long_name", "latitude"))
            call handleError(nf90_put_att(fileNCID, latID, "units", LAT_UNITS))
            call handleError(nf90_put_att(fileNCID, timeID, "long_name", "time"))
            call handleError(nf90_put_att(fileNCID, timeID, "units", timeUnits))
            call handleError(nf90_put_att(fileNCID, timeID, "calendar", "standard")) 

         ENDDO

            call handleError(nf90_enddef(fileNCID))
            call handleError(nf90_close(fileNCID))

        end subroutine nc3dFileCreate
        
        subroutine addGlobalAttributes()
            call handleError(nf90_put_att(fileNCID, NF90_GLOBAL, "Conventions", "COOARDS"))

        end subroutine addGlobalAttributes

        subroutine writeVar2D(var, nLon, nLat, nTime)
             integer :: nLon, nLat, nTime
             real :: var(nLon, nLat, nTime)

             call handleError(nf90_put_var(fileNCID, varID, var))
        end subroutine writeVar2D

        subroutine readVar2D(filename, var, nLon, nLat, nTime, varout)
            
             USE netcdf
 
             real :: lats(nLat), lons(nLon)
             real :: varout(nLon,nLat)

             integer                 :: ncid
             integer                 :: nLon, nLat, nTime
             integer                 :: lat_varid, lon_varid
             integer                 :: var_varid
             character (len = *) ::  filename
             !character(len = 255)    :: filename
             character(len = *)      :: var

             filename = trim(filename)
             !print *,'in readvar2d, trying to read:', TRIM(filename)
             !print *,'var:',var,'nlon:',nlon,'ntime:',ntime
             call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
             !call handleError( nf90_inq_varid(ncid, 'lat', lat_varid))
             !call handleError( nf90_inq_varid(ncid, 'lon', lon_varid))
             call handleError( nf90_inq_varid(ncid, var, var_varid))
       
             ! Read the latitude and longitude data.
             !call handleError( nf90_get_var(ncid, lat_varid, lats) )
             !call handleError( nf90_get_var(ncid, lon_varid, lons) )

             call handleError( nf90_get_var(ncid, var_varid, varout,         &
                      start=(/ 1, 1, nTime /),                               &
                      count=(/ nLon, nLat, 1 /) ) )
             !print *,'shape:',shape(varout)
             call handleError(nf90_close(ncid))

        end subroutine readVar2D

        subroutine readGOSATlength(filename, length)

         USE netcdf

         integer                 :: length 

         integer                 :: ncid, dimid
         integer                 :: nRecords
         character (len = *) ::  filename
         character(len = nf90_max_name) :: RecordDimName

         call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
         call handleError( nf90_inq_dimid(ncid, 'sounding_id', dimid) )
         call handleError( nf90_inquire_dimension(ncid, dimid,   &
                           name = RecordDimName,len=nRecords))
         length = nRecords
         print *,'gosat length:',length

        end subroutine readGOSATlength

      subroutine readGOSAT_real8(filename, var, grp,varout, nSoundings)

         USE netcdf

         real*8                  :: varout(nSoundings)

         integer                 :: ncid, dimid,numgrps
         integer                 :: ncids(1)
         integer                 :: nSoundings
         integer                 :: var_varid, grp_id
         character (len = *) ::  filename
         character (len = *) ::  grp
         character(len = nf90_max_name) :: RecordDimName
         character(len = *)      :: var
         character*80            :: name_in

         print *,'in readgosat2'
         print *,'grp:',grp
         print *,'filling nsoundings=',nsoundings
         !filename = trim(filename)
         !grp = trim(grp)
         print *,'file:',filename
         call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
         print *,'file opened'
         IF (grp .ne. '') THEN
          call handleError(nf90_inq_grpname(ncid, name_in))
          print *,'herehere'
          call handleError( nf90_inq_grp_ncid(ncid,grp,grp_id) )
          print *,'grp id:',grp_id
          print *,'var:',trim(var)
          call handleError( nf90_inq_varid(grp_id, var, var_varid))
          call handleError( nf90_get_var(grp_id, var_varid, varout,         &
                      start=(/ 1 /),                               &
                      count=(/ nSoundings /) ) )
         ELSE
          print *,'var:',var
          call handleError( nf90_inq_varid(ncid, var, var_varid))
          call handleError( nf90_get_var(ncid, var_varid, varout,         &
                      start=(/ 1 /),                               &
                      count=(/ nSoundings /) ) )
         ENDIF
         call handleError(nf90_close(ncid))
          print *,'done....'
        end subroutine readGOSAT_real8

      subroutine readGOSAT_real4(filename, var, grp,varout, nSoundings)

         USE netcdf

         real                    :: varout(nSoundings)

         integer                 :: ncid, dimid,numgrps
         integer                 :: ncids(1)
         integer                 :: nSoundings
         integer                 :: var_varid, grp_id
         character (len = *) ::  filename
         character (len = *) ::  grp
         character(len = nf90_max_name) :: RecordDimName
         character(len = *)      :: var
         character*80            :: name_in

         print *,'in readgosat2'
         print *,'grp:',grp

         !filename = trim(filename)
         !grp = trim(grp)
         print *,'file:',filename
         call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
         print *,'file opened'
         IF (grp .ne. '') THEN
          call handleError(nf90_inq_grpname(ncid, name_in))
          print *,'herehere'
          call handleError( nf90_inq_grp_ncid(ncid,grp,grp_id) )
          print *,'grp id:',grp_id
          print *,'var:',trim(var)
          call handleError( nf90_inq_varid(grp_id, var, var_varid))
          call handleError( nf90_get_var(grp_id, var_varid, varout,         &
                      start=(/ 1 /),                               &
                      count=(/ nSoundings /) ) )
         ELSE
          call handleError( nf90_inq_varid(ncid, var, var_varid))
          call handleError( nf90_get_var(ncid, var_varid, varout,         &
                      start=(/ 1 /),                               &
                      count=(/ nSoundings /) ) )
         ENDIF
         call handleError(nf90_close(ncid))

        end subroutine readGOSAT_real4


        subroutine readGOSAT_real(filename, var,varout, nSoundings)

         USE netcdf

         real                    :: varout(nSoundings)

         integer                 :: ncid, dimid
         integer                 :: nSoundings
         integer                 :: var_varid
         character (len = *) ::  filename
         character(len = nf90_max_name) :: RecordDimName
         character(len = *)      :: var

         filename = trim(filename)
         call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
         call handleError( nf90_inq_varid(ncid, var, var_varid))
         !call handleError( nf90_inq_dimid(ncid, 'soundings', dimid) )
         !call handleError( nf90_inquire_dimension(ncid, dimid,   &
         !                  name = RecordDimName,len=nRecords))
         !allocate(varout(nRecords))

         call handleError( nf90_get_var(ncid, var_varid, varout,         &
                      start=(/ 1 /),                               &
                      count=(/ nSoundings /) ) )
         
         call handleError(nf90_close(ncid))

        end subroutine readGOSAT_real

        subroutine readGOSAT_ID(filename, nSoundings, varout)

         USE netcdf

         !character (len=20), dimension(:), allocatable :: varout
         integer*8                 :: varout(:)
         integer,dimension(1)      :: dummy 
         integer                 :: ncid, dimid
         integer                 :: nSoundings
         integer                 :: var_varid, grp_id
         character (len = *) ::  filename
         character(len = nf90_max_name) :: RecordDimName

         print *,'in fun',nSoundings,' : ',shape(dummy)
         call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
         print *,filename
         !call handleError(nf90_inq_grp_ncid(ncid,'Sounding', grp_id))
         !call handleError( nf90_inq_varid(grp_id, 'sounding_id', var_varid))
         call handleError( nf90_inq_varid(ncid, 'sounding_id', var_varid))
         print *,'var_varid:',var_varid
         !call handleError( nf90_inq_dimid(ncid, 'soundings', dimid) )
         call handleError( nf90_inq_dimid(ncid, 'sounding_id', dimid) )
         !call handleError( nf90_inquire_dimension(ncid, dimid,   &
         !                  name = RecordDimName,len=nRecords))
         print *,'before',nSoundings
         !call handleError( nf90_get_var(grp_id, var_varid, varout,start=(/ 1 /),count=(/ nSoundings /) ) )
         !call handleError( nf90_get_var(grp_id, var_varid, varout,start=(/ 1 /),count=(/ nSoundings /) ) )
         call handleError( nf90_get_var(ncid, var_varid, varout,start=(/ 1 /),count=(/ nSoundings /) ) )
         !print *,'ID***:',varout
         print *,'after'
         call handleError(nf90_close(ncid))

         print *,'done w/ netcdf sub'
        end subroutine readGOSAT_ID

        subroutine readensVar2D(filename, var, nLon, nLat, nEns, nTime, varout)
        !If you feed nEns=-1, it will drop the nEns part of this netcdf pull
             USE netcdf

             real :: lats(nLat), lons(nLon)
             real :: varout(nLon,nLat)

             integer                 :: ncid
             integer                 :: nLon, nLat, nTime, nEns
             integer                 :: lat_varid, lon_varid
             integer                 :: var_varid
             character (len = *) ::  filename
             !character(len = 255)    :: filename
             character(len = *)      :: var

             filename = trim(filename)
             !print *,'in readvar2d, trying to read:', TRIM(filename)
             !print *,'var:',var,'nlon:',nlon,'ntime:',ntime
             call handleError(nf90_open(TRIM(filename),nf90_nowrite, ncid))
             !call handleError( nf90_inq_varid(ncid, 'lat', lat_varid))
             !call handleError( nf90_inq_varid(ncid, 'lon', lon_varid))
             call handleError( nf90_inq_varid(ncid, var, var_varid))

             ! Read the latitude and longitude data.
             !call handleError( nf90_get_var(ncid, lat_varid, lats) )
             !call handleError( nf90_get_var(ncid, lon_varid, lons) )

             IF(nEns .eq. -1) THEN 
             call handleError( nf90_get_var(ncid, var_varid, varout,         &
                      start=(/ 1, 1, nTime /),                               &
                      count=(/ nLon, nLat, 1 /) ) )
             ELSE
             call handleError( nf90_get_var(ncid, var_varid, varout,         &
                      start=(/ 1, 1, nEns, nTime /),                               &
                      count=(/ nLon, nLat, 1, 1 /) ) )
             ENDIF
             !print *,'shape:',shape(varout)
             call handleError(nf90_close(ncid))

        end subroutine readensVar2D

        !subroutine writeSlabVar2D(var, start, end, latIDX, lonIDX)
            !integer timeCnt 
            !integer, dimension (3) :: start, count
!
            !timeCnt = lastTime - outputTime;
            !start = (/ lonIDX-1, latIDX-1, 0 /);
            !count = (/ 1, 1, yearCnt+1 /);
    !
             ! Write data to variable.
        !    call handleError(nf90_put_var(fileNCID, varID, start, &
        !                count, var(:,:,outputTime:lastTime));
        !end subroutine writeSlabVar2D

        subroutine handleError(status)
            integer :: status
            if (status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                stop "Error with netcdf.  Stopped."
            end if
        end subroutine handleError

        subroutine gosatCloseError(status)
            integer :: status
            if (status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                print *, "NCDF output error.  Guessing there was no GOSAT output file to close...."
            end if
        end subroutine gosatCloseError

end  module ncCooardsFormat
