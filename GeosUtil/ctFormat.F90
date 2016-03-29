!
!   Filename: ncFormatter.f90
!   Purpose: This module contains functions which enable the user to create a
!            CF compliant NetCDF file for 2 and 3 dimensional spatial data.
!            
!   Author: RMcKeown
!   History: Initial version for files containing a single variable, latitude, 
!            longitude, vertical levels (optional), and time.
!

module ctFormat
    use netcdf
    implicit none

    ! Handles for the file, dimensions, and variables
    !integer, save :: fileNCID
    !integer, save :: latDimID, lonDimID, levDimID, timeDimID
    !integer, save :: latID, lonID, levID, timeID, varID
    !integer, save :: obsDimID, tracerDimID, tnameDimID
    !integer, save :: char100DimID, rnameDimID,gindicesDimID

    ! Handles for the file, dimensions, and variables
    integer, save :: ctOutputNCID
    !integer, save :: sLatID, sLonID, sLevID, sTimeID, sVarID,sNameID
    !integer, save :: gLatID, gLonID, gLevID, gTimeID, gVarID,gNameID
    !integer, save :: gVarID1,gVarID2,gVarID3

    !integer, save :: flaskVarID,nsamplesVarID,tnamesVarID,avgtimeVarID
    !integer, save :: sfchgtVarID, tnamesVarID,rindicesVarID,uVarID
    !integer, save :: vVarID,blhVarID,qVarID,pressureVarID

    ! Units required for CF compliance
    character (len = *), parameter :: LAT_UNITS = "degrees_north"
    character (len = *), parameter :: LON_UNITS = "degrees_east"
    character (len = *), parameter :: LEV_UNITS = "hPa"

    contains

 subroutine ncCTOutputCreate(filename, nObs)
        !   This subroutine creates a file to hold timeseries data for multiple stations.
        !   nObs:  number of observations

            character (len = *)      ::  filename
            integer                  ::  nObs

        !    character (len = *)      ::  varName1,varName2,varName3
        !   Need to add stationName to the argument list when ready
        !    integer                  :: soundingDimID, nSoundings
        !    integer                  :: sVarID2 ,sVarID3,sVarID4
        !    integer                  :: sVarID5 ,sVarID6,sVarID7,sVarID8
        integer    :: obspackVarID,flaskVarID,nsamplesVarID,tnamesVarID,avgtimeVarID
        integer    :: flaskobsVarID, sfchgtVarID, rnamesVarID,rindicesVarID,uVarID
        integer    :: vVarID,blhVarID,qVarID,pressureVarID,temperatureVarID

        integer    :: obsDimID, tracerDimID, tnameDimID,tnamelenDimID
        integer    :: char100DimID, rnameDimID,gindicesDimID

        !    character,dimension(:,:)   :: sounding_id
        !    real*8, dimension(:)       :: xco2_obs,xco2_uncert
        !    real*8, dimension(:)       :: xco2pbl_obs,xco2pbl_uncert
        !    real*8, dimension(:)       :: xco2ft_obs,xco2ft_uncert

            !integer, parameter :: mode_flag = IOR(nf90_hdf5, nf90_classic_model) ! | nf90_clobber
            integer, parameter :: mode_flag = nf90_hdf5
!           integer, parameter :: deflatelevel = 6

            ! Create file
            call handleError(nf90_create(path=trim(filename), cmode=mode_flag, ncid=ctOutputNCID))

            print *,'Opening/creating file ',trim(filename)

            ! Define Dimensions
            call handleError(nf90_def_dim(ctOutputNCID, "obs", nObs, obsDimID))
            call handleError(nf90_def_dim(ctOutputNCID, "tracer", 5, tracerDimID))
            call handleError(nf90_def_dim(ctOutputNCID, "tracer_name_len", 8, tnameDimID))
            call handleError(nf90_def_dim(ctOutputNCID, "char100", 100, char100DimID))
            call handleError(nf90_def_dim(ctOutputNCID, "region_name_len", 10, rnameDimID))
            call handleError(nf90_def_dim(ctOutputNCID, "grid_indices", 3, gindicesDimID))

            ! Define Variables
            call handleError(nf90_def_var(ncid=ctOutputNCID, name="obspack_id", xtype=nf90_char, &
             dimids=(/ obsDimID, char100DimID/), varid=obspackVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="flask",                   &
                 xtype=nf90_double, dimids=(/ obsDimID, tracerDimID/), varid=flaskVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="flask_obs",                   &
                 xtype=nf90_double, dimids=(/ obsDimID, tracerDimID/), varid=flaskobsVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="nsamples",                &
                 xtype=nf90_int, dimids= (/ obsDimID/), varid=nsamplesVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="tracer_names", xtype=nf90_char, &
                  dimids=(/ tracerDimID, tnameDimID/) , varid=tnamesVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="averaging_time",                &
                 xtype=nf90_int , dimids= (/ obsDimID/), varid=avgtimeVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="surface_height",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=sfchgtVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="region_name", xtype=nf90_char, &
                  dimids=(/ obsDimID, rnameDimID/) , varid=rnamesVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="region_indices",                &
                 xtype=nf90_int, dimids= (/ obsDimID , gindicesDimID/), varid=rindicesVarID))



            call handleError(nf90_def_var(ncid=ctOutputNCID, name="u",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=uVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="v",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=vVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="blh",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=blhVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="q",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=qVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="pressure",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=pressureVarID))

            call handleError(nf90_def_var(ncid=ctOutputNCID, name="temperature",                   &
                 xtype=nf90_double, dimids=(/ obsDimID/), varid=temperatureVarID))

            ! Define Attributes
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "long_name",                    &
            !                                  "ACOS Sounding ID"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "units", "MMDDYYhhmmss"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID2, "comment", "from scan start time in UTC"))

            !call handleError(nf90_put_att(gosatOutputNCID, sVarID, "long_name",                    &
            !                                  "Optimized value of Retrieved XCO2"))
            !call handleError(nf90_put_att(gosatOutputNCID, sVarID, "units", "ppm"))


            call handleError(nf90_enddef(ctOutputNCID))

            !call handleError(nf90_put_var(ctOutputNCID, sVarID2, sounding_id))
            !call handleError(nf90_put_var(ctOutputNCID, sVarID3, xco2_obs))
            !call handleError(nf90_put_var(ctOutputNCID, sVarID4, xco2_uncert))

            !call handleError(nf90_put_var(ctOutputNCID, sVarID5, xco2pbl_obs))
            !call handleError(nf90_put_var(ctOutputNCID, sVarID6, xco2pbl_uncert))
            !call handleError(nf90_put_var(ctOutputNCID, sVarID7, xco2ft_obs))
            !call handleError(nf90_put_var(ctOutputNCID, sVarID8, xco2ft_uncert))

            ! Sync the file to make sure data is saved
            call handleError(nf90_sync(ctOutputNCID))

        end subroutine ncCTOutputCreate

        subroutine ncCTOutputClose

          print *, "closing nc CT output file", ctOutputNCID

          call ctCloseError(nf90_close(ctOutputNCID))

        end subroutine ncCTOutputClose

       subroutine WRITE_REALVAR_TO_FILE_CT( V, M, var )

            integer             ::  M,localgVarID
            character (len = *) ::  var
            real*8              ::  V


             !print *,'writing ', V,' to ',M,'th spot'

            flush(6)

            !V = REAL(V * 10**6)
            !V2 = REAL(V2 * 10**6)
            !V3 = REAL(V3 * 10**6)

            V  = V 

            call handleError(nf90_inq_varid(ncid=ctOutputNCID, name=var,varid=localgVarID))

            call handleError(nf90_put_var(ncid=ctOutputNCID,varid=localgVarID , values=V,start=(/ M /) ) )

            ! Sync the file to make sure data is saved
            call handleError(nf90_sync(ctOutputNCID))

        end subroutine WRITE_REALVAR_TO_FILE_CT

        subroutine handleError(status)
            integer :: status
            if (status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                stop "Error with netcdf.  Stopped."
            end if
        end subroutine handleError

        subroutine ctCloseError(status)
            integer :: status
            if (status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                print *, "NCDF output error.  Guessing there was no CT output file to close...."
            end if
        end subroutine ctCloseError

        character(len=20) function str(k)
         !   "Convert an integer to string."
         integer, intent(in) :: k
         write (str, *) k
         str = adjustl(str)
        end function str

end  module ctFormat
