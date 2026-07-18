module mo_drydep

  !---------------------------------------------------------------------
  !       ... Dry deposition
  !---------------------------------------------------------------------

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use chem_mods,        only : gas_pcnst
  use pmgrid,           only : plev
  use spmd_utils,       only : masterproc
  use ppgrid,           only : pcols, begchunk, endchunk
  use mo_tracname,      only : solsym
  use cam_abortutils,   only : endrun
  use ioFileMod,        only : getfil
  use pio
  use cam_pio_utils,    only : cam_pio_openfile, cam_pio_closefile
  use cam_logfile,      only : iulog
  use dyn_grid,         only : get_dyn_grid_parm, get_horiz_grid_d
  use scamMod,          only : single_column

  use shr_drydep_mod,   only : nddvels =>  n_drydep, drydep_list, mapping
  use gas_drydep,       only : gas_drydep_init, gas_drydep_run, has_dvel, n_land_type

  use infnan,                only : nan, assignment(=)

  implicit none

  save

  interface drydep_inti
     module procedure dvel_inti_xactive
  end interface

  interface drydep
     module procedure drydep_fromlnd
  end interface

  private

  public :: drydep_inti, drydep, has_drydep
  public :: drydep_update
  public :: n_land_type, fraction_landuse, drydep_srf_file

  real(r8), protected, allocatable  :: fraction_landuse(:,:,:)

  type lnd_dvel_type
     real(r8), pointer :: dvel(:,:)   ! deposition velocity over land (cm/s)
  end type lnd_dvel_type

  type(lnd_dvel_type), allocatable :: lnd(:)
  character(len=SHR_KIND_CL) :: drydep_srf_file

contains

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_update( state, cam_in )
    use physics_types,   only : physics_state
    use camsrfexch,      only : cam_in_t

    type(physics_state), intent(in) :: state           ! Physics state variables
    type(cam_in_t),  intent(in) :: cam_in

    if (nddvels<1) return

    lnd(state%lchnk)%dvel => cam_in%depvel

  end subroutine drydep_update

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_fromlnd( ocnfrac, icefrac, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvelocity, dflx, mmr, &
                             tv, ncol, lchnk )

    !-------------------------------------------------------------------------------------
    ! CAM marshal for the portable core in gas_drydep.F90: computes the effective
    ! Henry's law coefficients, resolves this chunk's land deposition velocities,
    ! and calls gas_drydep_run (land velocities merged with ocean/sea-ice
    ! velocities computed by the Wesely scheme over land types 7 and 8).
    !-------------------------------------------------------------------------------------

    use shr_drydep_mod, only : shr_drydep_setHCoeff

#if (defined OFFLINE_DYN)
    use metdata, only: get_met_fields
#endif

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------

    real(r8), intent(in)      :: icefrac(pcols)
    real(r8), intent(in)      :: ocnfrac(pcols)
    integer,  intent(in)      :: ncol
    integer,  intent(in)      :: lchnk                    ! chunk number
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvelocity(ncol,gas_pcnst)    ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: heff(ncol,nddvels)
#if (defined OFFLINE_DYN)
    real(r8)  :: lndfrac(pcols)
    real(r8)  :: met_ocnfrac(pcols)
    real(r8)  :: met_icefrac(pcols)
#endif

    if (nddvels < 1) then
       ! no dry deposition species: the portable core is uninitialized and the
       ! shr_drydep tables unallocated (original code reached the same result
       ! through zero-iteration merge/flux loops)
       dvelocity(:,:) = 0._r8
       return
    end if

    !-------------------------------------------------------------------------------------
    ! define species-dependent parameters (temperature dependent)
    !-------------------------------------------------------------------------------------
    call shr_drydep_setHCoeff( ncol, sfc_temp, heff )

#if (defined OFFLINE_DYN)
    lndfrac(:ncol) = 1._r8 - ocnfrac(:ncol) - icefrac(:ncol)

    where( lndfrac(:ncol) < 0._r8 )
       lndfrac(:ncol) = 0._r8
    endwhere

    call get_met_fields(lndfrac, met_ocnfrac, met_icefrac, lchnk, ncol)

    call gas_drydep_run( ocnfrac      = ocnfrac, &
                         icefrac      = icefrac, &
                         ocnfrc_x     = met_ocnfrac, &
                         icefrc_x     = met_icefrac, &
                         sfc_temp     = sfc_temp, &
                         pressure_sfc = pressure_sfc, &
                         wind_speed   = wind_speed, &
                         spec_hum     = spec_hum, &
                         air_temp     = air_temp, &
                         pressure_10m = pressure_10m, &
                         rain         = rain, &
                         snow         = snow, &
                         solar_flux   = solar_flux, &
                         lnd_dvel     = lnd(lchnk)%dvel, &
                         heff         = heff, &
                         dvelocity    = dvelocity, &
                         dflx         = dflx, &
                         mmr          = mmr, &
                         tv           = tv, &
                         ncol         = ncol )
#else
    call gas_drydep_run( ocnfrac      = ocnfrac, &
                         icefrac      = icefrac, &
                         ocnfrc_x     = ocnfrac, &
                         icefrc_x     = icefrac, &
                         sfc_temp     = sfc_temp, &
                         pressure_sfc = pressure_sfc, &
                         wind_speed   = wind_speed, &
                         spec_hum     = spec_hum, &
                         air_temp     = air_temp, &
                         pressure_10m = pressure_10m, &
                         rain         = rain, &
                         snow         = snow, &
                         solar_flux   = solar_flux, &
                         lnd_dvel     = lnd(lchnk)%dvel, &
                         heff         = heff, &
                         dvelocity    = dvelocity, &
                         dflx         = dflx, &
                         mmr          = mmr, &
                         tv           = tv, &
                         ncol         = ncol )
#endif

  end subroutine drydep_fromlnd

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine dvel_inti_xactive( depvel_lnd_file )
    !-------------------------------------------------------------------------------------
    ! 	... intialize interactive drydep
    !-------------------------------------------------------------------------------------
    use dycore,        only : dycore_is
    use shr_drydep_mod, only : z0, rgso, rgss, ri, rclo, rcls, rlu, rac, foxd, drat
    use physconst,      only : karman, tmelt
    use phys_control,  only : phys_getopts

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    character(len=*), intent(in) :: depvel_lnd_file

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer :: i
    integer :: nlon_veg, nlat_veg, npft_veg
    integer :: dimid
    integer :: astat
    integer :: plon, plat
    integer :: ierr
    character(len=512) :: errmsg
    integer :: errflg

    real(r8), allocatable :: vegetation_map(:,:,:)
    real(r8), allocatable :: work(:,:)
    real(r8), allocatable :: landmask(:,:)
    real(r8), allocatable :: urban(:,:)
    real(r8), allocatable :: lake(:,:)
    real(r8), allocatable :: wetland(:,:)
    real(r8), allocatable :: lon_veg_edge(:)
    real(r8), allocatable :: lat_veg_edge(:)

    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid

    character(len=shr_kind_cl) :: locfn
    logical :: prog_modal_aero

    ! determine if modal aerosols are active so that fraction_landuse array is initialized for modal aerosal dry dep
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    ! with no dry deposition species there is nothing to initialize and the
    ! shr_drydep tables are unallocated (the original code reached the same
    ! state through its all(.not.has_dvel) early return below)
    if (nddvels < 1) return

    ! hand the portable core the host constants, species names, the deposition
    ! list, and the shr_drydep lookup tables; it resolves the species indices
    call gas_drydep_init( gas_pcnst_in   = gas_pcnst, &
                          plev_in        = plev, &
                          karman_in      = karman, &
                          tmelt_in       = tmelt, &
                          solsym_in      = solsym, &
                          n_drydep_in    = nddvels, &
                          drydep_list_in = drydep_list, &
                          mapping_in     = mapping, &
                          z0_in          = z0, &
                          rgso_in        = rgso, &
                          rgss_in        = rgss, &
                          ri_in          = ri, &
                          rclo_in        = rclo, &
                          rcls_in        = rcls, &
                          rlu_in         = rlu, &
                          rac_in         = rac, &
                          foxd_in        = foxd, &
                          drat_in        = drat, &
                          errmsg         = errmsg, &
                          errflg         = errflg )
    if (errflg /= 0) then
       call endrun(trim(errmsg))
    end if

    allocate( lnd(begchunk:endchunk) )

    if( masterproc ) then
       write(iulog,*) 'drydep_inti: following species have dry deposition'
       do i=1,nddvels
          if( len_trim(drydep_list(i)) > 0 ) then
             write(iulog,*) 'drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
          endif
       enddo
       write(iulog,*) 'drydep_inti:'
    endif

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    if (.not.prog_modal_aero) then
       return
    endif

    allocate( fraction_landuse(pcols,n_land_type, begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate fraction_landuse; error = ',astat
       call endrun('dvel_inti: failed to allocate fraction_landuse')
    end if
    fraction_landuse = nan

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    if(dycore_is('UNSTRUCTURED') ) then
       call get_landuse_and_soilw_from_file()
    else
       !---------------------------------------------------------------------------
       ! 	... read landuse map
       !---------------------------------------------------------------------------
       call getfil (depvel_lnd_file, locfn, 0)
       call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       ierr = pio_inq_dimid( piofile, 'lon', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlon_veg )
       ierr = pio_inq_dimid( piofile, 'lat', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlat_veg )
       ierr = pio_inq_dimid( piofile, 'pft', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, npft_veg )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegetation_map; error = ',astat
          call endrun('dvel_inti: failed to allocate vegetation_map')
       end if
       allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
            landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegetation_map; error = ',astat
          call endrun('dvel_inti: failed to allocate vegetation_map')
       end if
       allocate( lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegetation lon, lat arrays; error = ',astat
          call endrun('dvel_inti: failed to allocate vegetation lon, lat arrays')
       end if
       !---------------------------------------------------------------------------
       ! 	... read the vegetation map and landmask
       !---------------------------------------------------------------------------
       ierr = pio_inq_varid( piofile, 'PCT_PFT', vid )
       ierr = pio_get_var( piofile, vid, vegetation_map )

       ierr = pio_inq_varid( piofile, 'LANDMASK', vid )
       ierr = pio_get_var( piofile, vid, landmask )

       ierr = pio_inq_varid( piofile, 'PCT_URBAN', vid )
       ierr = pio_get_var( piofile, vid, urban )

       ierr = pio_inq_varid( piofile, 'PCT_LAKE', vid )
       ierr = pio_get_var( piofile, vid, lake )

       ierr = pio_inq_varid( piofile, 'PCT_WETLAND', vid )
       ierr = pio_get_var( piofile, vid, wetland )

       call cam_pio_closefile( piofile )

       !---------------------------------------------------------------------------
       ! scale vegetation, urban, lake, and wetland to fraction
       !---------------------------------------------------------------------------
       vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
       wetland(:,:)          = .01_r8 * wetland(:,:)
       lake(:,:)             = .01_r8 * lake(:,:)
       urban(:,:)            = .01_r8 * urban(:,:)
#ifdef DEBUG
       if(masterproc) then
          write(iulog,*) 'minmax vegetation_map ',minval(vegetation_map),maxval(vegetation_map)
          write(iulog,*) 'minmax wetland        ',minval(wetland),maxval(wetland)
          write(iulog,*) 'minmax landmask       ',minval(landmask),maxval(landmask)
       end if
#endif
       !---------------------------------------------------------------------------
       ! 	... define lat-lon of vegetation map (1x1)
       !---------------------------------------------------------------------------
       lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
       lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)

       !---------------------------------------------------------------------------
       ! 	... regrid to model grid
       !---------------------------------------------------------------------------
       call interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg_edge, &
            lon_veg_edge, landmask, urban, lake, &
            wetland, vegetation_map )

       deallocate( vegetation_map, work, stat=astat )
       deallocate( lon_veg_edge, lat_veg_edge, stat=astat )
       deallocate( landmask, urban, lake, wetland, stat=astat )
    endif  ! Unstructured grid

  end subroutine dvel_inti_xactive

  !-------------------------------------------------------------------------------------
  subroutine get_landuse_and_soilw_from_file()
    use ncdio_atm, only : infld

    logical :: readvar

    type(file_desc_t) :: piofile
    character(len=shr_kind_cl) :: locfn
    logical :: lexist

    if (len_trim(drydep_srf_file) == 0) then
       if (masterproc) then
          write(iulog,*)'**************************************'
          write(iulog,*)' get_landuse_and_soilw_from_file: INFO:'
          write(iulog,*)' drydep_srf_file not set:'
          write(iulog,*)' setting fraction_landuse to zero'
          write(iulog,*)'**************************************'
       end if
       fraction_landuse = 0._r8
       return
    end if

    call getfil (drydep_srf_file, locfn, 1, lexist)
    if(lexist) then
       call cam_pio_openfile(piofile, locfn, PIO_NOWRITE)

       call infld('fraction_landuse', piofile, 'ncol','class',1,pcols,1,n_land_type, begchunk,endchunk, &
            fraction_landuse, readvar, gridname='physgrid')
       if (.not. readvar) then
          if (masterproc) then
             write(iulog,*)'**************************************'
             write(iulog,*)'get_landuse_and_soilw_from_file: INFO:'
             write(iulog,*)' fraction_landuse not read from file: '
             write(iulog,*)' ', trim(locfn)
             write(iulog,*)' setting all values to zero'
             write(iulog,*)'**************************************'
          end if
          fraction_landuse = 0._r8
       end if

       call cam_pio_closefile(piofile)
    else
       call endrun('Unstructured grids require drydep_srf_file ')
    end if


  end subroutine get_landuse_and_soilw_from_file

  !-------------------------------------------------------------------------------------
  subroutine interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg_edge, &
                         lon_veg_edge, landmask, urban, lake, &
                         wetland, vegetation_map )

    use mo_constants, only : r2d
    use scamMod, only : latiop,loniop,scmlat,scmlon,scm_cambfb_mode
    use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
    use cam_initfiles, only: initial_file_get_id
    use dycore, only : dycore_is
    use phys_grid,     only : get_rlat_all_p, get_rlon_all_p, get_ncols_p

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer,  intent(in)         :: plon, plat, nlon_veg, nlat_veg, npft_veg
    real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
    real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
    real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: closelat,closelon
    integer :: latidx,lonidx

    integer, parameter           :: veg_ext = 20
    type(file_desc_t), pointer   :: piofile
    integer                      :: i, j, ii, jj, i_ndx, n
    integer, dimension(plon+1)   :: ind_lon
    integer, dimension(plat+1)  :: ind_lat
    real(r8)                         :: total_land
    real(r8), dimension(plon+1)      :: lon_edge
    real(r8), dimension(plat+1)     :: lat_edge
    real(r8)                         :: lat1, lon1
    real(r8)                         :: x1, x2, y1, y2, dx, dy
    real(r8)                         :: area, total_area
    real(r8), dimension(npft_veg+3)  :: fraction
    real(r8), dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
    integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext

    real(r8), allocatable :: lam(:), phi(:)

    logical, parameter :: has_npole = .true.
    integer :: ploniop,platiop
    real(r8) :: tmp_frac_lu(plon,n_land_type,plat)

    real(r8):: rlats(pcols), rlons(pcols)
    integer :: lchnk, ncol, icol
    logical :: found

    if(dycore_is('UNSTRUCTURED') ) then
       call endrun('mo_drydep::interp_map called for UNSTRUCTURED grid')
    endif

    allocate(lam(plon), phi(plat))
    call get_horiz_grid_d(plat, clat_d_out=phi)
    call get_horiz_grid_d(plon, clon_d_out=lam)

    if (single_column) then
       if (scm_cambfb_mode) then
          piofile => initial_file_get_id()
!         call shr_scam_getCloseLatLon(piofile%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          call shr_scam_getCloseLatLon(piofile,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          ploniop=size(loniop)
          platiop=size(latiop)
       else
          latidx=1
          lonidx=1
          ploniop=1
          platiop=1
       end if

       lon_edge(1) = loniop(lonidx) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d

       if (lonidx.lt.ploniop) then
          lon_edge(2) = loniop(lonidx+1) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
       else
          lon_edge(2) = lon_edge(1) + (loniop(2) - loniop(1)) * r2d
       end if

       lat_edge(1) = latiop(latidx) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d

       if (latidx.lt.platiop) then
          lat_edge(2) = latiop(latidx+1) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
       else
          lat_edge(2) = lat_edge(1) + (latiop(2) - latiop(1)) * r2d
       end if
    else
       do i = 1,plon
          lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
       end do
       lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
       if( .not. has_npole ) then
          do j = 1,plat+1
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
       else
          do j = 1,plat
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
          lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
       end if
    end if
    do j = 1,plat+1
       lat_edge(j) = min( lat_edge(j), 90._r8 )
       lat_edge(j) = max( lat_edge(j),-90._r8 )
    end do

    !-------------------------------------------------------------------------------------
    ! wrap around the longitudes
    !-------------------------------------------------------------------------------------
    do i = -veg_ext,0
       lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
       mapping_ext     (i) =              nlon_veg+i
    end do
    do i = 1,nlon_veg
       lon_veg_edge_ext(i) = lon_veg_edge(i)
       mapping_ext     (i) =              i
    end do
    do i = nlon_veg+1,nlon_veg+veg_ext
       lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
       mapping_ext     (i) =              i-nlon_veg
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : lon_edge ',lon_edge
    write(iulog,*) 'interp_map : lat_edge ',lat_edge
    write(iulog,*) 'interp_map : mapping_ext ',mapping_ext
#endif
    do j = 1,plon+1
       lon1 = lon_edge(j)
       do i = -veg_ext,nlon_veg+veg_ext
          dx = lon_veg_edge_ext(i  ) - lon1
          dy = lon_veg_edge_ext(i+1) - lon1
          if( dx*dy <= 0._r8 ) then
             ind_lon(j) = i
             exit
          end if
       end do
    end do

    do j = 1,plat+1
       lat1 = lat_edge(j)
       do i = 1,nlat_veg
          dx = lat_veg_edge(i  ) - lat1
          dy = lat_veg_edge(i+1) - lat1
          if( dx*dy <= 0._r8 ) then
             ind_lat(j) = i
             exit
          end if
       end do
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : ind_lon ',ind_lon
    write(iulog,*) 'interp_map : ind_lat ',ind_lat
#endif
    lat_loop : do j = 1,plat
       lon_loop : do i = 1,plon
          total_area       = 0._r8
          fraction         = 0._r8
          do jj = ind_lat(j),ind_lat(j+1)
             y1 = max( lat_edge(j),lat_veg_edge(jj) )
             y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) )
             dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
             do ii =ind_lon(i),ind_lon(i+1)
                i_ndx = mapping_ext(ii)
                x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
                x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) )
                dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
                area = dx * dy
                total_area = total_area + area
                !-----------------------------------------------------------------
                ! 	... special case for ocean grid point
                !-----------------------------------------------------------------
                if( nint(landmask(i_ndx,jj)) == 0 ) then
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area
                else
                   do n = 1,npft_veg
                      fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
                   end do
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
                   fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
                   fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
                   !-----------------------------------------------------------------
                   ! 	... check if land accounts for the whole area.
                   !           If not, the remaining area is in the ocean
                   !-----------------------------------------------------------------
                   total_land = sum(vegetation_map(i_ndx,jj,:)) &
                              + urban  (i_ndx,jj) &
                              + lake   (i_ndx,jj) &
                              + wetland(i_ndx,jj)
                   if( total_land < 1._r8 ) then
                      fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
                   end if
                end if
             end do
          end do
          !-------------------------------------------------------------------------------------
          ! 	... divide by total area of grid box
          !-------------------------------------------------------------------------------------
          fraction(:) = fraction(:)/total_area
          !-------------------------------------------------------------------------------------
          ! 	... make sure we don't have too much or too little
          !-------------------------------------------------------------------------------------
          if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
             fraction(:) = fraction(:)/sum(fraction)
          end if
          !-------------------------------------------------------------------------------------
          ! 	... map to Wesely land classification
          !-------------------------------------------------------------------------------------
          tmp_frac_lu(i, 1, j) =     fraction(20)
          tmp_frac_lu(i, 2, j) = sum(fraction(16:17))
          tmp_frac_lu(i, 3, j) = sum(fraction(13:15))
          tmp_frac_lu(i, 4, j) = sum(fraction( 5: 9))
          tmp_frac_lu(i, 5, j) = sum(fraction( 2: 4))
          tmp_frac_lu(i, 6, j) =     fraction(19)
          tmp_frac_lu(i, 7, j) =     fraction(18)
          tmp_frac_lu(i, 8, j) =     fraction( 1)
          tmp_frac_lu(i, 9, j) = 0._r8
          tmp_frac_lu(i,10, j) = 0._r8
          tmp_frac_lu(i,11, j) = sum(fraction(10:12))
       end do lon_loop
    end do lat_loop

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       call get_rlat_all_p(lchnk, ncol, rlats(:ncol))
       call get_rlon_all_p(lchnk, ncol, rlons(:ncol))
       do icol= 1,ncol
          found=.false.
          find_col: do j = 1,plat
             do i = 1,plon
                if (rlats(icol)==phi(j) .and. rlons(icol)==lam(i)) then
                   found=.true.
                   exit find_col
                endif
             enddo
          enddo find_col

          if (.not.found) call endrun('mo_drydep::interp_map not able find physics column coordinate')
          fraction_landuse(icol,1:n_land_type,lchnk) =  tmp_frac_lu(i,1:n_land_type,j)

       end do

       !-------------------------------------------------------------------------------------
       ! 	... make sure there are no out of range values
       !-------------------------------------------------------------------------------------
       where (fraction_landuse(:ncol,:n_land_type,lchnk) < 0._r8) fraction_landuse(:ncol,:n_land_type,lchnk) = 0._r8
       where (fraction_landuse(:ncol,:n_land_type,lchnk) > 1._r8) fraction_landuse(:ncol,:n_land_type,lchnk) = 1._r8
    end do

  end subroutine interp_map


  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  function has_drydep( name )

    character(len=*), intent(in) :: name

    logical :: has_drydep
    integer :: i

    has_drydep = .false.

    do i=1,nddvels
       if ( trim(name) == trim(drydep_list(i)) ) then
         has_drydep = .true.
         exit
       endif
    enddo

  endfunction has_drydep

end module mo_drydep
