module rad_constituents

!------------------------------------------------------------------------------------------------
!
! Provide constituent distributions and properties to the radiation and
! cloud microphysics routines.
!
! Retains: namelist I/O, gas handling, gas MMR retrieval (state/pbuf),
! initialization orchestration, diagnostics output, and cloud optics.
!
! Aerosol structural types, data, query/property routines, and namelist
! parsing are in radiative_aerosol (facade) backed by
! radiative_aerosol_definitions (core definitions).
!
! Aerosol MMR retrieval (state/pbuf access) is in aerosol_mmr_cam.
!
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: rga
use physics_types,  only: physics_state
use phys_control,   only: use_simple_phys
use radconstants,   only: nradgas, rad_gas_index
use phys_prop,      only: physprop_accum_unique_files, physprop_init
use cam_history,    only: addfld, fieldname_len, outfld, horiz_only
use physics_buffer, only: physics_buffer_desc, pbuf_get_field

use cam_abortutils, only: endrun
use cam_logfile,    only: iulog

! Import from radiative_aerosol_definitions (core definitions)
use radiative_aerosol_definitions, only: cs1, N_DIAG, n_rad_cnst, verbose, nl, &
                            rad_cnst_namelist_t, radcnst_namelist, active_calls, get_cam_idx

! Import from radiative_aerosol (facade: types, data, queries, parsing)
use radiative_aerosol, only: &
   n_mode_str, n_bin_str, &
   modes_t, bins_t, &
   aerosol_t, aerlist_t, &
   modelist_t, binlist_t, &
   modes, bins, &
   aerosollist, ma_list, sa_list, &
   rad_aer_get_info, &
   init_mode_comps, init_bin_comps, list_init1, list_init2, &
   parse_mode_defs, parse_bin_defs, parse_rad_specifier, &
   print_modes, print_bins

! Import from aerosol_mmr_cam (CAM-specific MMR retrieval)
use aerosol_mmr_cam, only: &
   aerosol_mmr_cam_init, &
   rad_cnst_get_aer_mmr, &
   rad_cnst_get_mam_mmr_idx, &
   rad_cnst_get_mode_num, &
   rad_cnst_get_mode_num_idx, &
   rad_cnst_get_bin_mmr_by_idx, &
   rad_cnst_get_bin_num, &
   rad_cnst_get_bin_num_idx, &
   rad_cnst_get_carma_mmr_idx, &
   rad_cnst_get_bin_mmr, &
   rad_aer_diag_init

implicit none
private
save

! Storage for gas components in the climate/diagnostic lists

type :: gas_t
   character(len=1)  :: source       ! A for state (advected), N for pbuf (non-advected), Z for zero
   character(len=64) :: camname      ! name of constituent in physics state or buffer
   character(len=32) :: mass_name    ! name for mass per layer field in history output
   integer           :: idx          ! index from constituents or from pbuf
end type gas_t

type :: gaslist_t
   integer                :: ngas
   character(len=2)       :: list_id  ! set to "  " for climate list, or two character integer
                                      ! (include leading zero) to identify diagnostic list
   type(gas_t), pointer   :: gas(:)   ! dimension(ngas) where ngas = nradgas is from radconstants
end type gaslist_t

type(gaslist_t), target :: gaslist(0:N_DIAG)  ! gasses used in climate/diagnostic calculations

! values for constituents with requested value of zero
real(r8), allocatable, target :: zero_cols(:,:)

! Public interfaces — routines in this module
public :: &
   rad_cnst_readnl,             &! read namelist values and parse
   rad_cnst_init,               &! find optics files and all constituents
   rad_cnst_get_info,           &! gas+aerosol info wrapper
   rad_cnst_get_gas,            &! return pointer to mmr for gasses
   rad_cnst_out                  ! output constituent diagnostics (mass per layer and column burden)

! Re-export aerosol MMR routines from aerosol_mmr_cam for backward compatibility
public :: &
   rad_cnst_get_aer_mmr,        &! return pointer to mmr for aerosols
   rad_cnst_get_mam_mmr_idx,    &! get constituent index of mam specie mmr (climate list only)
   rad_cnst_get_mode_num,       &! return mode number mixing ratio
   rad_cnst_get_mode_num_idx,   &! get constituent index of mode number m.r. (climate list only)
   rad_cnst_get_bin_mmr_by_idx, &
   rad_cnst_get_bin_num, &
   rad_cnst_get_bin_num_idx, &
   rad_cnst_get_carma_mmr_idx, &
   rad_cnst_get_bin_mmr

! Generic interface for rad_cnst_get_info — gas-only wrapper
interface rad_cnst_get_info
   module procedure rad_cnst_get_info_wrap
end interface

character(len=cs1), public :: iceopticsfile, liqopticsfile
character(len=32),  public :: icecldoptics,liqcldoptics
logical,            public :: oldcldoptics = .false.

! Namelist variables
character(len=cs1), dimension(n_mode_str) :: mode_defs   = ' '
character(len=cs1), dimension(n_bin_str) :: bin_defs   = ' '
character(len=cs1) :: rad_climate(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_1(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_2(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_3(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_4(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_5(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_6(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_7(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_8(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_9(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_10(n_rad_cnst) = ' '

!==============================================================================
contains

subroutine rad_cnst_readnl(nlfile)

   ! Read rad_cnst_nl namelist group.  Parse input.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=2) :: suffix
   character(len=1), pointer   :: ctype(:)
   character(len=*), parameter :: subname = 'rad_cnst_readnl'

   namelist /rad_cnst_nl/ mode_defs,     &
                          bin_defs,      &
                          rad_climate,   &
                          rad_diag_1,    &
                          rad_diag_2,    &
                          rad_diag_3,    &
                          rad_diag_4,    &
                          rad_diag_5,    &
                          rad_diag_6,    &
                          rad_diag_7,    &
                          rad_diag_8,    &
                          rad_diag_9,    &
                          rad_diag_10,   &
                          iceopticsfile, &
                          liqopticsfile, &
                          icecldoptics,  &
                          liqcldoptics,  &
                          oldcldoptics

   !-----------------------------------------------------------------------------

   if (use_simple_phys) return

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rad_cnst_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rad_cnst_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (mode_defs,     len(mode_defs(1))*n_mode_str,     mpichar, 0, mpicom)
   call mpibcast (bin_defs,      len(bin_defs(1))*n_bin_str,       mpichar, 0, mpicom)
   call mpibcast (rad_climate,   len(rad_climate(1))*n_rad_cnst,   mpichar, 0, mpicom)
   call mpibcast (rad_diag_1,    len(rad_diag_1(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_2,    len(rad_diag_2(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_3,    len(rad_diag_3(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_4,    len(rad_diag_4(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_5,    len(rad_diag_5(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_6,    len(rad_diag_6(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_7,    len(rad_diag_7(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_8,    len(rad_diag_8(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_9,    len(rad_diag_9(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_10,   len(rad_diag_10(1))*n_rad_cnst,   mpichar, 0, mpicom)
   call mpibcast (iceopticsfile, len(iceopticsfile),               mpichar, 0, mpicom)
   call mpibcast (liqopticsfile, len(liqopticsfile),               mpichar, 0, mpicom)
   call mpibcast (liqcldoptics,  len(liqcldoptics),                mpichar, 0, mpicom)
   call mpibcast (icecldoptics,  len(icecldoptics),                mpichar, 0, mpicom)
   call mpibcast (oldcldoptics,  1,                                mpilog , 0, mpicom)
#endif

   ! Parse the namelist input strings

   ! Mode definition stings
   call parse_mode_defs(mode_defs, modes)

   ! Bin definition stings
   call parse_bin_defs(bin_defs, bins)

   ! Lists of externally mixed entities for climate and diagnostic calculations
   do i = 0,N_DIAG
      select case (i)
      case(0)
         call parse_rad_specifier(rad_climate, radcnst_namelist(i))
      case (1)
         call parse_rad_specifier(rad_diag_1, radcnst_namelist(i))
      case (2)
         call parse_rad_specifier(rad_diag_2, radcnst_namelist(i))
      case (3)
         call parse_rad_specifier(rad_diag_3, radcnst_namelist(i))
      case (4)
         call parse_rad_specifier(rad_diag_4, radcnst_namelist(i))
      case (5)
         call parse_rad_specifier(rad_diag_5, radcnst_namelist(i))
      case (6)
         call parse_rad_specifier(rad_diag_6, radcnst_namelist(i))
      case (7)
         call parse_rad_specifier(rad_diag_7, radcnst_namelist(i))
      case (8)
         call parse_rad_specifier(rad_diag_8, radcnst_namelist(i))
      case (9)
         call parse_rad_specifier(rad_diag_9, radcnst_namelist(i))
      case (10)
         call parse_rad_specifier(rad_diag_10, radcnst_namelist(i))
      end select
   enddo

   ! were there any constituents specified for the nth diagnostic call?
   ! if so, radiation will make a call with those consituents
   active_calls(:) = (radcnst_namelist(:)%ncnst > 0)

   ! Initialize the gas and aerosol lists with the information from the
   ! namelist.  This is done here so that this information is available via
   ! the query functions at the time when the register methods are called.

   ! Set the list_id fields which distinquish the climate and diagnostic lists
   do i = 0, N_DIAG
      if (active_calls(i)) then
         if (i > 0) then
            write(suffix, fmt = '(i2.2)') i
         else
            suffix='  '
         end if
         aerosollist(i)%list_id = suffix
         gaslist(i)%list_id     = suffix
         ma_list(i)%list_id     = suffix
         sa_list(i)%list_id     = suffix
      end if
   end do

   ! Create a list of the unique set of filenames containing property data

   ! Start with the bulk aerosol species in the climate/diagnostic lists.
   ! The physprop_accum_unique_files routine has the side effect of returning the number
   ! of bulk aerosols in each list (they're identified by type='A').
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call physprop_accum_unique_files(radcnst_namelist(i)%radname, radcnst_namelist(i)%type)
      endif
   enddo

   ! Add physprop files for the species from the mode definitions.
   do i = 1, modes%nmodes
      allocate(ctype(modes%comps(i)%nspec))
      ctype = 'A'
      call physprop_accum_unique_files(modes%comps(i)%props, ctype)
      deallocate(ctype)
   end do

   ! Add physprop files for the species from the bin definitions.
   do i = 1, bins%nbins
      allocate(ctype(bins%comps(i)%nspec))
      ctype = 'A'
      call physprop_accum_unique_files(bins%comps(i)%props, ctype)
      deallocate(ctype)
   end do

   ! Initialize the gas, bulk aerosol, and modal aerosol lists.  This step splits the
   ! input climate/diagnostic lists into the corresponding gas, bulk and modal aerosol
   ! lists.
   if (masterproc) write(iulog,*) nl//subname//': Radiation constituent lists:'
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call gas_list_init1(radcnst_namelist(i), gaslist(i))
         call list_init1(radcnst_namelist(i), aerosollist(i), ma_list(i), sa_list(i))

         if (masterproc .and. verbose) then
            call print_lists(gaslist(i), aerosollist(i), ma_list(i), sa_list(i))
         end if

      end if
   end do

   if (masterproc .and. verbose) call print_modes(modes)
   if (masterproc .and. verbose) call print_bins(bins)

end subroutine rad_cnst_readnl

subroutine rad_cnst_init()

   ! The initialization of the gas and aerosol lists is finished by
   ! 1) read the physprop files
   ! 2) find the index of each constituent in the constituent or physics buffer arrays
   ! 3) find the index of the aerosol constituents used to access its properties from the
   !    physprop module.

   integer :: i
   logical, parameter :: stricttest = .true.
   character(len=*), parameter :: subname = 'rad_cnst_init'
   !-----------------------------------------------------------------------------

   ! memory to point to if zero value requested
   allocate(zero_cols(pcols,pver))
   zero_cols = 0._r8

   ! Initialize zero_cols in aerosol_mmr_cam
   call aerosol_mmr_cam_init()

   ! Allocate storage for the physical properties of each aerosol; read properties from
   ! the data files.
   call physprop_init()

   ! Start checking that specified radiative constituents are present in the constituent
   ! or physics buffer arrays.
   if (masterproc) write(iulog,*) nl//subname//': checking for radiative constituents'

   ! Finish initializing the mode definitions.
   call init_mode_comps(modes)

   ! Finish initializing the bin definitions.
   call init_bin_comps(bins)

   ! Finish initializing the aerosol lists.
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call list_init2(aerosollist(i), ma_list(i), sa_list(i))
      end if
   end do

   ! Finish initializing the gas lists (resolve constituent indices).
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call gas_list_init2(gaslist(i))
      end if
   end do

   ! Check that all gases supported by the radiative transfer code have been specified.
   if (stricttest) then
      do i = 1, nradgas
         if (gaslist(0)%gas(i)%source .eq. 'Z' ) then
            call endrun(subname//': list of radiative gasses must include all radiation gasses for the climate specication')
         endif
      enddo
   endif

   ! Initialize history output of climate diagnostic quantities
   call rad_gas_diag_init(gaslist(0))
   call rad_aer_diag_init(aerosollist(0))


end subroutine rad_cnst_init

subroutine gas_list_init1(namelist, gaslist)

   ! Initialize gas list from parsed namelist data.

   type(rad_cnst_namelist_t), intent(in)    :: namelist
   type(gaslist_t),           intent(inout) :: gaslist

   ! Local variables
   integer :: ii, igas, istat
   character(len=*), parameter :: routine = 'gas_list_init1'
   !-----------------------------------------------------------------------------

   ! nradgas is set by the radiative transfer code
   gaslist%ngas = nradgas

   allocate(gaslist%gas(gaslist%ngas), stat=istat)
   if (istat /= 0) call endrun(routine//': allocate ERROR; gas list components')

   ! Initialize sources to zero (default for unspecified gases)
   do igas = 1, gaslist%ngas
      gaslist%gas(igas)%source  = 'Z'
      gaslist%gas(igas)%camname = ' '
   end do

   ! Populate gas entries from 'G' type namelist entries
   do ii = 1, namelist%ncnst
      if (namelist%type(ii) /= 'G') cycle

      if (masterproc .and. verbose) &
         write(iulog,*) "  rad namelist spec: "// trim(namelist%source(ii)) &
         //":"//trim(namelist%camname(ii))//":"//trim(namelist%radname(ii))

      ! rad_gas_index will abort on illegal names
      igas = rad_gas_index(namelist%radname(ii))

      gaslist%gas(igas)%source  = namelist%source(ii)
      gaslist%gas(igas)%camname = namelist%camname(ii)
   end do

end subroutine gas_list_init1

!================================================================================================

subroutine gas_list_init2(gaslist)

   ! Resolve constituent indices for gas list entries.

   type(gaslist_t), intent(inout) :: gaslist

   ! Local variables
   integer :: i
   character(len=*), parameter :: routine = 'gas_list_init2'
   !-----------------------------------------------------------------------------

   do i = 1, gaslist%ngas
      gaslist%gas(i)%idx = get_cam_idx(gaslist%gas(i)%source, gaslist%gas(i)%camname, routine)
   end do

end subroutine gas_list_init2

!================================================================================================

subroutine rad_cnst_get_info_wrap(list_idx, gasnames, aernames, &
                                  use_data_o3, ngas, naero, nmodes, nbins)

   ! Wrapper that provides the original rad_cnst_get_info interface.
   ! Gas arguments are handled locally; aerosol arguments are delegated
   ! to rad_aer_get_info in radiative_aerosol.

   ! Arguments
   integer,                     intent(in)  :: list_idx
   character(len=64), optional, intent(out) :: gasnames(:)
   character(len=64), optional, intent(out) :: aernames(:)
   logical,           optional, intent(out) :: use_data_o3
   integer,           optional, intent(out) :: naero
   integer,           optional, intent(out) :: ngas
   integer,           optional, intent(out) :: nmodes
   integer,           optional, intent(out) :: nbins

   ! Local variables
   type(gaslist_t),  pointer :: g_list
   integer          :: i, igas, gaslen
   character(len=1) :: source
   character(len=*), parameter :: subname = 'rad_cnst_get_info'
   !-----------------------------------------------------------------------------

   ! Delegate aerosol arguments to radiative_aerosol
   call rad_aer_get_info(list_idx, aernames=aernames, naero=naero, &
                         nmodes=nmodes, nbins=nbins)

   ! Handle gas arguments locally
   g_list => gaslist(list_idx)

   if (present(ngas)) then
      ngas = g_list%ngas
   endif

   if (present(gasnames)) then
      gaslen = size(gasnames)
      if (gaslen < g_list%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', g_list%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gasnames too short')
      end if
      do i = 1, g_list%ngas
         gasnames(i) = g_list%gas(i)%camname
      end do
   end if

   if (present(use_data_o3)) then
      igas = rad_gas_index('O3')
      source = g_list%gas(igas)%source
      use_data_o3 = .false.
      if (source == 'N') use_data_o3 = .true.
   endif

end subroutine rad_cnst_get_info_wrap

!================================================================================================

subroutine rad_cnst_get_gas(list_idx, gasname, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the gas from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),            intent(in) :: gasname
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: lchnk
   integer :: igas
   integer :: idx
   character(len=1) :: source
   type(gaslist_t), pointer :: list
   character(len=*), parameter :: subname = 'rad_cnst_get_gas'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      list => gaslist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   lchnk = state%lchnk

   ! Get index of gas in internal arrays.  rad_gas_index will abort if the
   ! specified gasname is not recognized by the radiative transfer code.
   igas = rad_gas_index(trim(gasname))

   ! Get data source
   source = list%gas(igas)%source
   idx    = list%gas(igas)%idx
   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_gas

subroutine rad_cnst_out(list_idx, state, pbuf)

   ! Output the mass per layer, and total column burdens for gas and aerosol
   ! constituents in either the climate or diagnostic lists

   ! Arguments
   integer,                     intent(in) :: list_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)


   ! Local variables
   integer :: i, naer, ngas, lchnk, ncol
   integer :: idx
   character(len=1)  :: source
   character(len=32) :: name, cbname
   real(r8)          :: mass(pcols,pver)
   real(r8)          :: cb(pcols)
   real(r8), pointer :: mmr(:,:)
   type(aerlist_t), pointer :: aerlist
   type(gaslist_t), pointer :: g_list
   character(len=*), parameter :: subname = 'rad_cnst_out'
   !-----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointer with requested aerosol list
   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   naer = aerlist%numaerosols
   do i = 1, naer

      source = aerlist%aer(i)%source
      idx    = aerlist%aer(i)%idx
      name   = aerlist%aer(i)%mass_name
      ! construct name for column burden field by replacing the 'm_' prefix by 'cb_'
      cbname = 'cb_' // name(3:len_trim(name))

      select case( source )
      case ('A')
         mmr => state%q(:,:,idx)
      case ('N')
         call pbuf_get_field(pbuf, idx, mmr)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

   ! Associate pointer with requested gas list
   g_list => gaslist(list_idx)

   ngas = g_list%ngas
   do i = 1, ngas

      source = g_list%gas(i)%source
      idx    = g_list%gas(i)%idx
      name   = g_list%gas(i)%mass_name
      cbname = 'cb_' // name(3:len_trim(name))
      select case( source )
      case ('A')
         mmr => state%q(:,:,idx)
      case ('N')
         call pbuf_get_field(pbuf, idx, mmr)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

end subroutine rad_cnst_out

subroutine rad_gas_diag_init(glist)

! Add diagnostic fields to the master fieldlist.

   type(gaslist_t), intent(inout) :: glist

   integer :: i, ngas
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   ngas = glist%ngas
   if (ngas == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = glist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, ngas

      ! construct names for mass per layer diagnostics
      name = 'm_' // trim(glist%gas(i)%camname) // trim(suffix)
      glist%gas(i)%mass_name = name
      long_name = trim(glist%gas(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), (/ 'lev' /), 'A', 'kg/m^2', trim(long_name))

      ! construct names for column burden diagnostics
      name = 'cb_' // trim(glist%gas(i)%camname) // trim(suffix)
      long_name = trim(glist%gas(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), horiz_only, 'A', 'kg/m^2', trim(long_name))

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_gas_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_gas_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_gas_diag_init

!================================================================================================

subroutine print_lists(gas_list, aer_list, ma_list, sa_list)

   ! Print summary of gas, bulk and modal aerosol lists.  This is just the information
   ! read from the namelist.

   use radconstants, only: gascnst=>gaslist

   type(aerlist_t),  intent(in) :: aer_list
   type(gaslist_t),  intent(in) :: gas_list
   type(modelist_t), intent(in) :: ma_list
   type(binlist_t),  intent(in) :: sa_list

   integer :: i, id

   if (len_trim(gas_list%list_id) == 0) then
      write(iulog,*) nl//' gas list for climate calculations'
   else
      write(iulog,*) nl//' gas list for diag'//gas_list%list_id//' calculations'
   end if

   do i = 1, nradgas
      if (gas_list%gas(i)%source .eq. 'N') then
         write(iulog,*) '  '//gas_list%gas(i)%source//':'//gascnst(i)//' has pbuf name:'//&
                        trim(gas_list%gas(i)%camname)
      else if (gas_list%gas(i)%source .eq. 'A') then
         write(iulog,*) '  '//gas_list%gas(i)%source//':'//gascnst(i)//' has constituents name:'//&
                        trim(gas_list%gas(i)%camname)
      endif
   enddo

   if (len_trim(aer_list%list_id) == 0) then
      write(iulog,*) nl//' bulk aerosol list for climate calculations'
   else
      write(iulog,*) nl//' bulk aerosol list for diag'//aer_list%list_id//' calculations'
   end if

   do i = 1, aer_list%numaerosols
      write(iulog,*) '  '//trim(aer_list%aer(i)%source)//':'//trim(aer_list%aer(i)%camname)//&
                     ' optics and phys props in :'//trim(aer_list%aer(i)%physprop_file)
   enddo

   if (len_trim(ma_list%list_id) == 0) then
      write(iulog,*) nl//' modal aerosol list for climate calculations'
   else
      write(iulog,*) nl//' modal aerosol list for diag'//ma_list%list_id//' calculations'
   end if

   do i = 1, ma_list%nmodes
      id = ma_list%idx(i)
      write(iulog,*) '  '//trim(modes%names(id))
   enddo

   if (len_trim(sa_list%list_id) == 0) then
      write(iulog,*) nl//' bin aerosol list for climate calculations'
   else
      write(iulog,*) nl//' bin aerosol list for diag'//sa_list%list_id//' calculations'
   end if

   do i = 1, sa_list%nbins
      id = sa_list%idx(i)
      write(iulog,*) '  '//trim(bins%names(id))
   enddo

end subroutine print_lists

!================================================================================================


end module rad_constituents
