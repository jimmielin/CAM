!===============================================================================
! Calcite (CaCO3) mass fraction of emitted dust, read from a lat-lon dataset and
! interpolated to the physics grid. Used to speciate dust emissions into calcite
! for dust heterogeneous chemistry.
!===============================================================================
module dust_calcite_mod
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,      only: iulog
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun

  implicit none
  private

  public :: dust_calcite_init
  public :: dust_calcite_frac

  real(r8), allocatable ::  dust_calcite_frac(:,:)  ! calcite mass fraction of emitted dust (pcols,begchunk:endchunk)

contains

  !=============================================================================
  !=============================================================================
  subroutine dust_calcite_init( dust_calcite_file )
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use ppgrid,           only: begchunk, endchunk, pcols
    use mo_constants,     only: pi, d2r
    use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE, PIO_NOERR
    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use cam_pio_utils,    only: cam_pio_openfile, cam_pio_closefile
    use ioFileMod,        only: getfil

    character(len=*), intent(in) :: dust_calcite_file

    real(r8), allocatable :: calcite_frac_in(:,:)  ! temporary input array
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
    character(len=cl)     :: infile
    integer :: did, vid, nlat, nlon
    type(file_desc_t) :: ncid

    type(interp_type) :: lon_wgts, lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols)
    integer :: c, ncols, ierr
    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'dust_calcite_mod: calcite mass fraction dataset: ', trim(dust_calcite_file)
    end if

    ! read in the calcite mass fraction map, same layout as the soil erodibility dataset

    ! Get file name.
    call getfil(dust_calcite_file, infile, 0)
    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlon )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlat )

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    allocate(calcite_frac_in(nlon,nlat))

    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, dst_lons  )

    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, dst_lats  )

    ierr = pio_inq_varid( ncid, 'caco3_mass_frac', vid )
    if( ierr /= PIO_NOERR ) then
       call endrun('dust_calcite_init: variable caco3_mass_frac not found in '//trim(infile))
    end if
    ierr = pio_get_var( ncid, vid, calcite_frac_in )

    call cam_pio_closefile(ncid)

    !-----------------------------------------------------------------------
    !     	... convert to radians and setup regridding
    !-----------------------------------------------------------------------
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

    allocate( dust_calcite_frac(pcols,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_calcite_init: failed to allocate dust_calcite_frac, ierr = ',ierr
       call endrun('dust_calcite_init: failed to allocate dust_calcite_frac')
    end if

    !-----------------------------------------------------------------------
    !     	... regrid ..
    !-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(calcite_frac_in(:,:), nlon,nlat , dust_calcite_frac(:,c), ncols, lon_wgts,lat_wgts)

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( calcite_frac_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_calcite_init: failed to deallocate calcite_frac_in, ierr = ',ierr
       call endrun('dust_calcite_init: failed to deallocate calcite_frac_in')
    end if

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine dust_calcite_init

end module dust_calcite_mod
