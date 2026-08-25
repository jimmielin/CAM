!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use modal_aero_data,  only: ntot_amode, ndst=>nDust
  use cam_logfile,      only: iulog
  use shr_dust_emis_mod,only: is_dust_emis_zender, is_zender_soil_erod_from_atm

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active
  public :: dust_calcite_names
  public :: dust_calcite_indices

  integer, protected :: dust_nbin != 2
  integer, protected :: dust_nnum != 2
  character(len=6), protected, allocatable :: dust_names(:)

  ! calcite (CaCO3) tracer of each dust bin, for dust heterogeneous chemistry: the calcite
  ! mass fraction of the emitted dust is emitted into it instead of the dust tracer
  character(len=6), protected, allocatable :: dust_calcite_names(:)   ! blank where the bin has none
  integer, protected, allocatable :: dust_calcite_indices(:)          ! constituent index, 0 where none
  logical :: dust_calcite_active = .false.

  real(r8), allocatable :: dust_dmt_grd(:)
  real(r8), allocatable :: dust_emis_sclfctr(:)

  integer , protected, allocatable :: dust_indices(:)
  real(r8), allocatable :: dust_dmt_vwr(:)
  real(r8), allocatable :: dust_stk_crc(:)

  real(r8)          :: dust_emis_fact = 0._r8     ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'none'    ! full pathname for soil erodibility dataset
  character(len=cl) :: dust_calcite_file = 'none' ! full pathname for the calcite mass fraction dataset
                                                  ! (dust heterogeneous chemistry); 'none' = no calcite speciation

  logical :: dust_active = .false.

 contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, masterprocid, mpi_character, mpi_real8, mpi_success
    use shr_dust_emis_mod, only: shr_dust_emis_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file, dust_calcite_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(soil_erod_file, len(soil_erod_file), mpi_character, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//' MPI_BCAST ERROR: soil_erod_file')
    end if
    call mpi_bcast(dust_emis_fact, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//' MPI_BCAST ERROR: dust_emis_fact')
    end if
    call mpi_bcast(dust_calcite_file, len(dust_calcite_file), mpi_character, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//' MPI_BCAST ERROR: dust_calcite_file')
    end if

    call shr_dust_emis_readnl(mpicom, 'drv_flds_in')

    if ((soil_erod_file /= 'none') .and. (.not.is_zender_soil_erod_from_atm())) then
       call endrun(subname//': should not specify soil_erod_file if Zender soil erosion is not in CAM')
    end if

    if (masterproc) then
       if (is_dust_emis_zender()) then
          write(iulog,*) subname,': Zender_2003 dust emission method is being used.'
       end if
       if (is_zender_soil_erod_from_atm()) then
          write(iulog,*) subname,': Zender soil erod file is handled in atm'
          write(iulog,*) subname,': soil_erod_file = ',trim(soil_erod_file)
          write(iulog,*) subname,': dust_emis_fact = ',dust_emis_fact
       end if
       if (dust_calcite_file /= 'none') then
          write(iulog,*) subname,': dust_calcite_file = ',trim(dust_calcite_file)
       end if
    end if

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use dust_calcite_mod, only: dust_calcite_init
    use constituents,  only: cnst_get_ind
    use aerosol_instances_mod, only: aerosol_instances_get_props, aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use dust_common,   only: dust_set_params

    integer :: l, m, mm, ndx, nspec, iaermod
    character(len=32) :: spec_name
    character(len=32) :: spec_type
    integer, parameter :: mymodes(7) = (/ 2, 1, 3, 4, 5, 6, 7 /) ! tricky order ...
    class(aerosol_properties), pointer :: aero_props_modal
    integer :: dust_modes(ndst)   ! aerosol mode of each dust bin

    dust_nbin = ndst
    dust_nnum = ndst

    allocate( dust_names(2*ndst) )
    allocate( dust_indices(2*ndst) )
    allocate( dust_calcite_names(ndst) )
    allocate( dust_calcite_indices(ndst) )
    dust_calcite_names(:) = ' '
    dust_calcite_indices(:) = 0
    allocate( dust_dmt_grd(ndst+1) )
    allocate( dust_emis_sclfctr(ndst) )
    allocate( dust_dmt_vwr(ndst) )
    allocate( dust_stk_crc(ndst) )

    ! dmleung edited the mass fraction of the emitted dust size distribution. 27 Oct 2025 ++
    ! The new mass fraction comes from Jun Meng et al. (2022) and MERRA-2.
    ! Jun Meng's table indicates 2.1 % mass for 0.1-1 um and 97.9 % mass for 1-10 um.
    ! ref: https://zenodo.org/records/6344524
    ! MERRA-2 dust emissions indicate 6 % mass for 0.1-1 um (bin1) and 94 % for 1-10 um (bin2-5).
    ! dmleung adopts 2.1 % mass for 0.1-1 um and 97.9 % mass for 1-10 um for dust.
    ! Distributing more mass to accumulation mode allows a longer lifetime of dust, reducing
    ! low dust biases over remote oceans and reducing high dust biases over the Sahara.
    ! This change impacts both Zender_2003 dust and Leung_2023 dust.
    if ( ntot_amode == 3 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.021_r8,0.979_r8 /)
    elseif ( ntot_amode == 4 .or. ntot_amode == 5 ) then
       dust_dmt_grd(:) = (/ 0.01e-6_r8, 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
       dust_emis_sclfctr(:) = (/ 1.65E-05_r8, 0.021_r8, 0.979_r8 /)
    else if( ntot_amode == 7 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.12_r8, 0.88_r8 /)
    endif
    ! dmleung --

    ! Find modal properties object from factory
    aero_props_modal => null()
    do iaermod = 1, aerosol_instances_get_num_models()
       aero_props_modal => aerosol_instances_get_props(iaermod, 0)
       if (aero_props_modal%model_is('MAM')) exit
       aero_props_modal => null()
    end do

    ndx = 0
    do mm = 1, ntot_amode
       m = mymodes(mm)
       nspec = aero_props_modal%nspecies(m)
       do l = 1, nspec
          call aero_props_modal%get(m, l, specname=spec_name)
          if (spec_name(:3) == 'dst') then
             ndx=ndx+1
             dust_names(ndx) = spec_name
             dust_names(ndst+ndx) = 'num_'//spec_name(5:)
             call cnst_get_ind(dust_names(     ndx), dust_indices(     ndx))
             call cnst_get_ind(dust_names(ndst+ndx), dust_indices(ndst+ndx))
             dust_modes(ndx) = m
          endif
       enddo
    enddo

    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return

    ! calcite tracer of each dust bin, found by species type in the bin's mode
    do ndx = 1, ndst
       m = dust_modes(ndx)
       do l = 1, aero_props_modal%nspecies(m)
          call aero_props_modal%get(m, l, spectype=spec_type, specname=spec_name)
          if (spec_type == 'calcite') then
             dust_calcite_names(ndx) = spec_name
             call cnst_get_ind(dust_calcite_names(ndx), dust_calcite_indices(ndx))
          end if
       end do
    end do
    dust_calcite_active = any(dust_calcite_indices(:) > 0)
    if (dust_calcite_active .and. dust_calcite_file == 'none') then
       call endrun('dust_init: calcite tracers are present but dust_calcite_file is not set')
    end if

    if (is_zender_soil_erod_from_atm()) then
       call  soil_erod_init( dust_emis_fact, soil_erod_file )
    end if

    if (dust_calcite_file /= 'none') then
       call dust_calcite_init( dust_calcite_file )
    end if

    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use dust_calcite_mod, only : dust_calcite_frac
    use mo_constants,  only : dust_density
    use physconst,     only : pi

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum, ical
    real(r8) :: x_mton
    real(r8) :: frac_cal   ! calcite mass fraction of the emitted dust
    real(r8),parameter :: soil_erod_threshold = 0.1_r8

    ! set dust emissions

    if (is_zender_soil_erod_from_atm()) then   ! Zender_2003 dust emissions
       col_loop1: do i = 1,ncol
          soil_erod(i) = soil_erodibility( i, lchnk )
          if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

          ! rebin and adjust dust emissons.
          do m = 1,dust_nbin
             idst = dust_indices(m)
             cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
                  * dust_emis_sclfctr(m)*soil_erod(i)/dust_emis_fact*1.15_r8
             x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))
             inum = dust_indices(m+dust_nbin)
             cflx(i,inum) = cflx(i,idst)*x_mton
          enddo
       enddo col_loop1
    else ! Leung_2023 dust emissions

       col_loop2: do i = 1,ncol
          ! rebin and adjust dust emissons.
          do m = 1,dust_nbin
             idst = dust_indices(m)

             cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
                  * dust_emis_sclfctr(m) / dust_emis_fact
             x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))
             inum = dust_indices(m+dust_nbin)
             cflx(i,inum) = cflx(i,idst)*x_mton
          enddo
       enddo col_loop2
    end if

    ! emit the calcite mass fraction of the dust into the calcite tracer of each bin that
    ! has one; the number fluxes above are from the total mass, so they are unchanged
    if (dust_calcite_active) then
       do i = 1,ncol
          frac_cal = dust_calcite_frac(i,lchnk)
          do m = 1,dust_nbin
             ical = dust_calcite_indices(m)
             if (ical > 0) then
                idst = dust_indices(m)
                cflx(i,ical) = cflx(i,idst)*frac_cal
                cflx(i,idst) = cflx(i,idst)*(1._r8 - frac_cal)
             end if
          enddo
       enddo
    end if

  end subroutine dust_emis

end module dust_model
