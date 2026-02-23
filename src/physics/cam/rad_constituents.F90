module rad_constituents

!------------------------------------------------------------------------------------------------
!
! Provide constituent distributions and properties to the radiation and
! cloud microphysics routines.
!
! The logic to control which constituents are used in the climate calculations
! and which are used in diagnostic radiation calculations is contained in this module.
!
! Structural types, data, query routines and property-access routines have been
! moved to aerosol_definition_mod.  This module re-exports them for backward
! compatibility and retains: namelist parsing, MMR/number retrieval (state/pbuf),
! initialization orchestration, and diagnostics output.
!
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: rga
use physics_types,  only: physics_state
use phys_control,   only: use_simple_phys
use constituents,   only: cnst_get_ind
use radconstants,   only: nradgas, rad_gas_index
use phys_prop,      only: physprop_accum_unique_files, physprop_init, &
                          physprop_get_id, ot_length
use cam_history,    only: addfld, fieldname_len, outfld, horiz_only
use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index

use cam_abortutils, only: endrun
use cam_logfile,    only: iulog

! Import structural types, data, query & property-access routines from aerosol_definition_mod
! and re-export for backward compatibility
use aerosol_definition_mod, only: &
   cs1, N_DIAG, n_mode_str, n_bin_str, n_rad_cnst, &
   mode_component_t, modes_t, bin_component_t, bins_t, &
   rad_cnst_namelist_t, gas_t, gaslist_t, aerosol_t, aerlist_t, &
   modelist_t, binlist_t, &
   modes, bins, radcnst_namelist, active_calls, &
   gaslist, aerosollist, ma_list, sa_list, &
   mode_type_names, spec_type_names, num_mode_types, num_spec_types, &
   num_bin_morphs, bin_morph_names, verbose, nl, &
   rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec, &
   rad_cnst_get_mode_idx, rad_cnst_get_spec_idx, &
   rad_cnst_get_call_list, rad_cnst_num_name, &
   rad_cnst_get_mode_props, rad_cnst_get_aer_props, &
   rad_cnst_get_bin_props_by_idx, rad_cnst_get_bin_props, &
   rad_cnst_get_aer_idx, &
   init_mode_comps, init_bin_comps, list_init1, list_init2, get_cam_idx

implicit none
private
save

! Public interfaces — routines kept in this module
public :: &
   rad_cnst_readnl,             &! read namelist values and parse
   rad_cnst_init,               &! find optics files and all constituents
   rad_cnst_get_gas,            &! return pointer to mmr for gasses
   rad_cnst_get_aer_mmr,        &! return pointer to mmr for aerosols
   rad_cnst_get_mam_mmr_idx,    &! get constituent index of mam specie mmr (climate list only)
   rad_cnst_get_mode_num,       &! return mode number mixing ratio
   rad_cnst_get_mode_num_idx,   &! get constituent index of mode number m.r. (climate list only)
   rad_cnst_out,                &! output constituent diagnostics (mass per layer and column burden)
   rad_cnst_get_bin_mmr_by_idx, &
   rad_cnst_get_bin_num, &
   rad_cnst_get_bin_num_idx, &
   rad_cnst_get_carma_mmr_idx, &
   rad_cnst_get_bin_mmr

! Re-exported public interfaces from aerosol_definition_mod (backward compatibility)
public :: N_DIAG
public :: &
   rad_cnst_get_info,           &
   rad_cnst_get_mode_idx,       &
   rad_cnst_get_spec_idx,       &
   rad_cnst_get_aer_props,      &
   rad_cnst_get_mode_props,     &
   rad_cnst_get_call_list,      &
   rad_cnst_get_bin_props_by_idx, &
   rad_cnst_get_info_by_bin, &
   rad_cnst_get_info_by_bin_spec, &
   rad_cnst_get_bin_props, &
   rad_cnst_num_name

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

! values for constituents with requested value of zero
real(r8), allocatable, target :: zero_cols(:,:)

! define generic interface for MMR retrieval
interface rad_cnst_get_aer_mmr
   module procedure rad_cnst_get_aer_mmr_by_idx
   module procedure rad_cnst_get_mam_mmr_by_idx
end interface

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
         call list_init1(radcnst_namelist(i), gaslist(i), aerosollist(i), ma_list(i), sa_list(i))

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

   ! Finish initializing the gas, bulk aerosol, and mode lists.
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call list_init2(gaslist(i), aerosollist(i), ma_list(i), sa_list(i))
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

subroutine rad_aer_diag_init(alist)

! Add diagnostic fields to the master fieldlist.

   type(aerlist_t), intent(inout) :: alist

   integer :: i, naer
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   naer = alist%numaerosols
   if (naer == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = alist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, naer

      ! construct names for mass per layer diagnostic fields
      name = 'm_' // trim(alist%aer(i)%camname) // trim(suffix)
      alist%aer(i)%mass_name = name
      long_name = trim(alist%aer(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), (/ 'lev' /), 'A', 'kg/m^2', trim(long_name))

      ! construct names for column burden diagnostic fields
      name = 'cb_' // trim(alist%aer(i)%camname) // trim(suffix)
      long_name = trim(alist%aer(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), horiz_only, 'A', 'kg/m^2', trim(long_name))

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_aer_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_aer_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_aer_diag_init

subroutine parse_mode_defs(nl_in, modes)

   ! Parse the mode definition specifiers.  The specifiers are of the form:
   !
   ! 'mode_name:mode_type:=',
   !  'source_num_a:camname_num_a:source_num_c:camname_num_c:num_mr:+',
   !  'source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file[:+]'[,]
   !  ['source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file][:+][']
   !
   ! where the ':' separated fields are:
   ! mode_name -- name of the mode.
   ! mode_type -- type of mode.  Valid values are from the MAM code.
   ! =         -- this line terminator identifies the initial string in a
   !              mode definition
   ! +         -- this line terminator indicates that the mode definition is
   !              continued in the next string
   ! source_num_a  -- Source of interstitial number mixing ratio,  'A', 'N', or 'Z'
   ! camname_num_a -- the name of the interstitial number component.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_num_c  -- Source of cloud borne number mixing ratio,  'A', 'N', or 'Z'
   ! camname_num_c -- the name of the cloud borne number component.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_mmr_a  -- Source of interstitial specie mass mixing ratio,  'A', 'N' or 'Z'
   ! camname_mmr_a -- the name of the interstitial specie.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_mmr_c  -- Source of cloud borne specie mass mixing ratio,  'A', 'N' or 'Z'
   ! camname_mmr_c -- the name of the cloud borne specie.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! spec_type -- species type.  Valid values far from the MAM code, except that
   !              the value 'num_mr' designates a number mixing ratio and has no
   !              associated field for the prop_file.  There can only be one entry
   !              with the num_mr type in a mode definition.
   ! prop_file -- For aerosol species this is a filename, which is
   !              identified by a ".nc" suffix.  The file contains optical and
   !              other physical properties of the aerosol.
   !
   ! A mode definition must contain only 1 string for the number mixing ratio components
   ! and at least 1 string for the species.


   character(len=*), intent(inout) :: nl_in(:)    ! namelist input (blanks are removed on output)
   type(modes_t),    intent(inout) :: modes       ! structure containing parsed input

   ! Local variables
   integer :: m
   integer :: istat
   integer :: nmodes, nstr
   integer :: mbeg, mcur
   integer :: nspec, ispec
   integer :: strlen, iend, ipos
   logical :: num_mr_found
   character(len=*), parameter :: routine = 'parse_mode_defs'
   character(len=len(nl_in(1))) :: tmpstr
   character(len=1)  :: tmp_src_a
   character(len=32) :: tmp_name_a
   character(len=1)  :: tmp_src_c
   character(len=32) :: tmp_name_c
   character(len=32) :: tmp_type
   !-------------------------------------------------------------------------

   ! Determine number of modes defined by counting number of strings that are
   ! terminated by ':='
   ! (algorithm stops counting at first blank element).
   nmodes = 0
   nstr = 0
   do m = 1, n_mode_str

      if (len_trim(nl_in(m)) == 0) exit
      nstr = nstr + 1

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(nl_in(m))
      nl_in(m) = tmpstr
      do
         strlen = len_trim(nl_in(m))
         ipos = index(nl_in(m), ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = nl_in(m)(:ipos-1) // nl_in(m)(ipos+1:strlen)
         nl_in(m) = tmpstr
      end do
      ! count strings with ':=' terminator
      if (nl_in(m)(strlen-1:strlen) == ':=') nmodes = nmodes + 1

   end do
   modes%nmodes = nmodes

   ! return if no modes defined
   if (nmodes == 0) return

   ! allocate components that depend on nmodes
   allocate( &
      modes%names(nmodes),  &
      modes%types(nmodes),  &
      modes%comps(nmodes),  &
      stat=istat )
   if (istat > 0) then
      write(iulog,*) routine//': ERROR: cannot allocate storage for modes.  nmodes=', nmodes
      call endrun(routine//': ERROR allocating storage for modes')
   end if

   mcur = 1              ! index of current string being processed

   ! loop over modes
   do m = 1, nmodes

      mbeg = mcur  ! remember the first string of a mode

      ! check that first string in mode definition is ':=' terminated
      iend = len_trim(nl_in(mcur))
      if (nl_in(mcur)(iend-1:iend) /= ':=') call parse_error('= not found', nl_in(mcur))

      ! count species in mode definition.  definition will contain 1 string with
      ! with a ':+' terminator for each specie
      nspec = 0
      mcur = mcur + 1
      do
         iend = len_trim(nl_in(mcur))
         if (nl_in(mcur)(iend-1:iend) /= ':+') exit
         nspec = nspec + 1
         mcur = mcur + 1
      end do

      ! a mode must have at least one specie
      if (nspec == 0) call parse_error('mode must have at least one specie', nl_in(mbeg))

      ! allocate components that depend on number of species
      allocate( &
         modes%comps(m)%source_mmr_a(nspec),  &
         modes%comps(m)%camname_mmr_a(nspec), &
         modes%comps(m)%source_mmr_c(nspec),  &
         modes%comps(m)%camname_mmr_c(nspec), &
         modes%comps(m)%type(nspec),          &
         modes%comps(m)%props(nspec),         &
         stat=istat)

      if (istat > 0) then
         write(iulog,*) routine//': ERROR: cannot allocate storage for species.  nspec=', nspec
         call endrun(routine//': ERROR allocating storage for species')
      end if

      ! initialize components
      modes%comps(m)%nspec         = nspec
      modes%comps(m)%source_num_a  = ' '
      modes%comps(m)%camname_num_a = ' '
      modes%comps(m)%source_num_c  = ' '
      modes%comps(m)%camname_num_c = ' '
      do ispec = 1, nspec
         modes%comps(m)%source_mmr_a(ispec)  = ' '
         modes%comps(m)%camname_mmr_a(ispec) = ' '
         modes%comps(m)%source_mmr_c(ispec)  = ' '
         modes%comps(m)%camname_mmr_c(ispec) = ' '
         modes%comps(m)%type(ispec)          = ' '
         modes%comps(m)%props(ispec)         = ' '
      end do

      ! return to first string in mode definition
      mcur = mbeg
      tmpstr = nl_in(mcur)

      ! mode name
      ipos = index(tmpstr, ':')
      if (ipos < 2) call parse_error('mode name not found', tmpstr)
      modes%names(m) = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! mode type
      ipos = index(tmpstr, ':')
      if (ipos == 0) call parse_error('mode type not found', tmpstr)
      ! check for valid mode type
      call check_mode_type(tmpstr, 1, ipos-1)
      modes%types(m) = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! mode type must be followed by '='
      if (tmpstr(1:1) /= '=') call parse_error('= not found', tmpstr)

      ! move to next string
      mcur = mcur + 1
      tmpstr = nl_in(mcur)

      ! process mode component strings
      num_mr_found = .false.   ! keep track of whether number mixing ratio component is found
      ispec = 0                ! keep track of the number of species found
      do

         ! source of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find source field first', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_a = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_a = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! source of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find a source field', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_c = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_c = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! component type
         ipos = scan(tmpstr, ': ')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)

         if (tmpstr(:ipos-1) == 'num_mr') then

            ! there can only be one number mixing ratio component
            if (num_mr_found) call parse_error('more than 1 number component', nl_in(mcur))

            num_mr_found = .true.
            modes%comps(m)%source_num_a  = tmp_src_a
            modes%comps(m)%camname_num_a = tmp_name_a
            modes%comps(m)%source_num_c  = tmp_src_c
            modes%comps(m)%camname_num_c = tmp_name_c
            tmpstr                       = tmpstr(ipos+1:)

         else

            ! check for valid specie type
            call check_specie_type(tmpstr, 1, ipos-1)
            tmp_type = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ! get the properties file
            ipos = scan(tmpstr, ': ')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)
            ! check for valid filename -- must have .nc extension
            if (tmpstr(ipos-3:ipos-1) /= '.nc') &
               call parse_error('filename not valid', tmpstr)

            ispec = ispec + 1
            modes%comps(m)%source_mmr_a(ispec)  = tmp_src_a
            modes%comps(m)%camname_mmr_a(ispec) = tmp_name_a
            modes%comps(m)%source_mmr_c(ispec)  = tmp_src_c
            modes%comps(m)%camname_mmr_c(ispec) = tmp_name_c
            modes%comps(m)%type(ispec)          = tmp_type
            modes%comps(m)%props(ispec)         = tmpstr(:ipos-1)
            tmpstr                              = tmpstr(ipos+1:)
         end if

         ! check if there are more components.  either the current character is
         ! a ' ' which means this string is the final mode component, or the character
         ! is a '+' which means there are more components
         if (tmpstr(1:1) == ' ') exit

         if (tmpstr(1:1) /= '+') &
               call parse_error('+ field not found', tmpstr)

         ! continue to next component...
         mcur = mcur + 1
         tmpstr = nl_in(mcur)
      end do

      ! check that a number component was found
      if (.not. num_mr_found) call parse_error('number component not found', nl_in(mbeg))

      ! check that the right number of species were found
      if (ispec /= nspec) call parse_error('component parsing got wrong number of species', nl_in(mbeg))

      ! continue to next mode...
      mcur = mcur + 1
      tmpstr = nl_in(mcur)
   end do

   !------------------------------------------------------------------------------------------------
   contains
   !------------------------------------------------------------------------------------------------

   ! internal subroutines used for error checking and reporting

   subroutine parse_error(msg, str)

      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: str

      write(iulog,*) routine//': ERROR: '//msg
      write(iulog,*) ' input string: '//trim(str)
      call endrun(routine//': ERROR: '//msg)

   end subroutine parse_error

   !------------------------------------------------------------------------------------------------

   subroutine check_specie_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie

      integer :: i

      do i = 1, num_spec_types
         if (str(ib:ie) == trim(spec_type_names(i))) return
      end do

      call parse_error('specie type not valid', str(ib:ie))

   end subroutine check_specie_type

   !------------------------------------------------------------------------------------------------

   subroutine check_mode_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie  ! begin, end character of mode type substring

      integer :: i

      do i = 1, num_mode_types
         if (str(ib:ie) == trim(mode_type_names(i))) return
      end do

      call parse_error('mode type not valid', str(ib:ie))

   end subroutine check_mode_type

   !------------------------------------------------------------------------------------------------

end subroutine parse_mode_defs

!================================================================================================

subroutine parse_bin_defs(nl_in, bins)

   ! Parse the bin definition specifiers.  The specifiers are of the form:
   !
   ! 'bin_name:=',
   !  'source_num_a:camname_num_a:source_num_c:camname_num_c:num_mr:+',
   !  'source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file[:+]'[,]
   !  ['source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file][:+][']
   !
   ! where the ':' separated fields are:
   ! bin_name -- name of the bin.
   ! =         -- this line terminator identifies the initial string in a
   !              mode definition
   ! +         -- this line terminator indicates that the mode definition is
   !              continued in the next string
   ! source_num_a  -- Source of interstitial number mixing ratio,  'A', 'N', or 'Z'
   ! camname_num_a -- the name of the interstitial number component.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_num_c  -- Source of cloud borne number mixing ratio,  'A', 'N', or 'Z'
   ! camname_num_c -- the name of the cloud borne number component.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_mmr_a  -- Source of interstitial specie mass mixing ratio,  'A', 'N' or 'Z'
   ! camname_mmr_a -- the name of the interstitial specie.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! source_mmr_c  -- Source of cloud borne specie mass mixing ratio,  'A', 'N' or 'Z'
   ! camname_mmr_c -- the name of the cloud borne specie.  This name must be
   !                  registered in the constituent arrays when source=A or in the
   !                  physics buffer when source=N
   ! spec_type -- species type.  Valid values are particle, shell, and core.
   ! prop_file -- For aerosol species this is a filename, which is
   !              identified by a ".nc" suffix.  The file contains optical and
   !              other physical properties of the aerosol.
   !
   ! A bin definition must contain at least 1 string for the species and can contain
   ! a maximum of 1 particle type.


   character(len=*), intent(inout) :: nl_in(:)    ! namelist input (blanks are removed on output)
   type(bins_t),    intent(inout) :: bins       ! structure containing parsed input

   ! Local variables
   logical :: num_mr_found, mass_mr_found
   logical :: particle_mr_found
   integer :: m
   integer :: istat
   integer :: nbins, nstr, istr
   integer :: mbeg, mcur
   integer :: nspec, ispec
   integer :: strlen, ibeg, iend, ipos
   logical :: part_mr_found
   character(len=*), parameter :: routine = 'parse_bin_defs'
   character(len=len(nl_in(1))) :: tmpstr
   character(len=1)  :: tmp_src_a
   character(len=32) :: tmp_name_a
   character(len=1)  :: tmp_src_c
   character(len=32) :: tmp_name_c
   character(len=32) :: tmp_type
   character(len=32) :: tmp_morph
   !-------------------------------------------------------------------------

   ! Determine number of bins defined by counting number of strings that are
   ! terminated by ':='
   ! (algorithm stops counting at first blank element).
   nbins = 0
   nstr = 0
   do m = 1, n_bin_str

      if (len_trim(nl_in(m)) == 0) exit
      nstr = nstr + 1

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(nl_in(m))
      nl_in(m) = tmpstr
      do
         strlen = len_trim(nl_in(m))
         ipos = index(nl_in(m), ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = nl_in(m)(:ipos-1) // nl_in(m)(ipos+1:strlen)
         nl_in(m) = tmpstr
      end do
      ! count strings with ':=' terminator
      if (nl_in(m)(strlen-1:strlen) == ':=') nbins = nbins + 1

   end do
   bins%nbins = nbins

   ! return if no bins defined
   if (nbins == 0) return

   ! allocate components that depend on nmodes
   allocate( &
      bins%names(nbins),  &
      bins%comps(nbins),  &
      stat=istat )
   if (istat > 0) then
      write(iulog,*) routine//': ERROR: cannot allocate storage for bins.  nbins=', nbins
      call endrun(routine//': ERROR allocating storage for bins')
   end if

   mcur = 1              ! index of current string being processed

   ! loop over bins
   bins_loop: do m = 1, nbins

      mbeg = mcur  ! remember the first string of a bin

      ! check that first string in bin definition is ':=' terminated
      iend = len_trim(nl_in(mcur))
      if (nl_in(mcur)(iend-1:iend) /= ':=') call parse_error('= not found', nl_in(mcur))

      ! count species in bin definition.  definition will contain 1 string with
      ! with a ':+' terminator for each specie
      nspec = 0
      mcur = mcur + 1
      do
         iend = len_trim(nl_in(mcur))
         if (nl_in(mcur)(iend-1:iend) /=    ':+') exit
         if (nl_in(mcur)(iend-4:iend) /= 'mmr:+') nspec = nspec + 1
         mcur = mcur + 1
      end do

      ! a bin must have at least one specie
      if (nspec == 0) call parse_error('bin must have at least one specie', nl_in(mbeg))

      ! allocate components that depend on number of species
      allocate( &
         bins%comps(m)%source_mmr_a(nspec),  &
         bins%comps(m)%camname_mmr_a(nspec), &
         bins%comps(m)%source_mmr_c(nspec),  &
         bins%comps(m)%camname_mmr_c(nspec), &
         bins%comps(m)%type(nspec),          &
         bins%comps(m)%morph(nspec),          &
         bins%comps(m)%props(nspec),         &
         stat=istat)

      if (istat > 0) then
         write(iulog,*) routine//': ERROR: cannot allocate storage for species.  nspec=', nspec
         call endrun(routine//': ERROR allocating storage for species')
      end if

      ! initialize components
      bins%comps(m)%nspec         = nspec
      bins%comps(m)%source_num_a  = ' '
      bins%comps(m)%camname_num_a = ' '
      bins%comps(m)%source_num_c  = ' '
      bins%comps(m)%camname_num_c = ' '
      bins%comps(m)%source_mass_a  = 'NOTSET'
      bins%comps(m)%camname_mass_a = 'NOTSET'
      bins%comps(m)%source_mass_c  = 'NOTSET'
      bins%comps(m)%camname_mass_c = 'NOTSET'
      do ispec = 1, nspec
         bins%comps(m)%source_mmr_a(ispec)  = ' '
         bins%comps(m)%camname_mmr_a(ispec) = ' '
         bins%comps(m)%source_mmr_c(ispec)  = ' '
         bins%comps(m)%camname_mmr_c(ispec) = ' '
         bins%comps(m)%type(ispec)          = ' '
         bins%comps(m)%props(ispec)         = ' '
      end do

      ! return to first string in mode definition
      mcur = mbeg
      tmpstr = nl_in(mcur)

      ! bin name
      ipos = index(tmpstr, ':')
      if (ipos < 2) call parse_error('bin name not found', tmpstr)
      bins%names(m)  = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! bin name must be followed by '='
      if (tmpstr(1:1) /= '=') call parse_error('= not found', tmpstr)

      ! move to next string
      mcur = mcur + 1
      tmpstr = nl_in(mcur)

      ! process bin component strings
      particle_mr_found = .false.   ! keep track of whether particle mixing ratio component is found
      num_mr_found = .false.        ! keep track of whether number mixing ratio component is found
      mass_mr_found = .false.        ! keep track of whether number mixing ratio component is found
      ispec = 0                ! keep track of the number of species found
      comps_loop: do

         ! source of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find source field first', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_a = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_a = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! source of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find a source field', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_c = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_c = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! component type
         ipos = scan(tmpstr, ': ')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)

         if (tmpstr(:ipos-1) == 'num') then

            ! there can only be one number mixing ratio component
            if (num_mr_found) call parse_error('more than 1 number component', nl_in(mcur))

            num_mr_found = .true.
            bins%comps(m)%source_num_a  = tmp_src_a
            bins%comps(m)%camname_num_a = tmp_name_a
            bins%comps(m)%source_num_c  = tmp_src_c
            bins%comps(m)%camname_num_c = tmp_name_c
            tmpstr                      = tmpstr(ipos+1:)

         else if (tmpstr(:ipos-1) == 'mmr') then

            ! there can only be one number mixing ratio component
            if (mass_mr_found) call parse_error('more than 1 mass mixing ratio component', nl_in(mcur))

            mass_mr_found = .true.
            bins%comps(m)%source_mass_a  = tmp_src_a
            bins%comps(m)%camname_mass_a = tmp_name_a
            bins%comps(m)%source_mass_c  = tmp_src_c
            bins%comps(m)%camname_mass_c = tmp_name_c
            tmpstr                       = tmpstr(ipos+1:)

         else

            ! check for valid species type
            call check_bin_type(tmpstr, 1, ipos-1)
            tmp_type = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ipos = index(tmpstr, ':')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)

            ! check for valid species type
            call check_bin_morph(tmpstr, 1, ipos-1)
            tmp_morph = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ! get the properties file
            ipos = scan(tmpstr, ': ')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)

             ! check for valid filename -- must have .nc extension
            if (tmpstr(ipos-3:ipos-1) /= '.nc') &
               call parse_error('filename not valid', tmpstr)

            ispec = ispec + 1

            bins%comps(m)%source_mmr_a(ispec)  = tmp_src_a
            bins%comps(m)%camname_mmr_a(ispec) = tmp_name_a
            bins%comps(m)%source_mmr_c(ispec)  = tmp_src_c
            bins%comps(m)%camname_mmr_c(ispec) = tmp_name_c
            bins%comps(m)%type(ispec)          = tmp_type
            bins%comps(m)%morph(ispec)         = tmp_morph

            bins%comps(m)%props(ispec)         = tmpstr(:ipos-1)
            tmpstr                             = tmpstr(ipos+1:)

         endif

         ! check if there are more components.  either the current character is
         ! a ' ' which means this string is the final mode component, or the character
         ! is a '+' which means there are more components
         if (tmpstr(1:1) == ' ') then
            exit comps_loop
         endif

         if (tmpstr(1:1) /= '+') &
               call parse_error('+ field not found', tmpstr)

         ! continue to next component...
         mcur = mcur + 1
         tmpstr = nl_in(mcur)
      end do comps_loop


      ! check that a number component was found
      if (.not. num_mr_found) call parse_error('number component not found', nl_in(mbeg))

      ! check that the right number of species were found
      if (ispec /= nspec) then
         write(*,*) 'ispec, nspec = ',ispec, nspec
         call parse_error('component parsing got wrong number of species', nl_in(mbeg))
      endif

      ! continue to next bin...
      mcur = mcur + 1
      tmpstr = nl_in(mcur)
   end do bins_loop

   !------------------------------------------------------------------------------------------------
   contains
   !------------------------------------------------------------------------------------------------

   ! internal subroutines used for error checking and reporting

   subroutine parse_error(msg, str)

      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: str

      write(iulog,*) routine//': ERROR: '//msg
      write(iulog,*) ' input string: '//trim(str)
      call endrun(routine//': ERROR: '//msg)

   end subroutine parse_error

   !------------------------------------------------------------------------------------------------

   subroutine check_bin_morph(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie

      integer :: i

      do i = 1, num_bin_morphs
         if (str(ib:ie) == trim(bin_morph_names(i))) return
      end do

      call parse_error('bin morph not valid', str(ib:ie))

   end subroutine check_bin_morph

   !------------------------------------------------------------------------------------------------
   subroutine check_bin_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie  ! begin, end character of mode type substring

      integer :: i

      do i = 1, num_spec_types
         if (str(ib:ie) == trim(spec_type_names(i))) return
      end do

      call parse_error('bin species type not valid', str(ib:ie))

   end subroutine check_bin_type

   !------------------------------------------------------------------------------------------------

end subroutine parse_bin_defs

!================================================================================================

subroutine parse_rad_specifier(specifier, namelist_data)

!-----------------------------------------------------------------------------
! Private method for parsing the radiation namelist specifiers.  The specifiers
! are of the form 'source_camname:radname' where:
! source  -- either 'N' for pbuf (non-advected) or 'A' for state (advected)
! camname -- the name of a constituent that must be found in the constituent
!            component of the state when source=A or in the physics buffer
!            when source=N
! radname -- For gases this is a name that identifies the constituent to the
!            radiative transfer codes.  These names are contained in the
!            radconstants module.  For aerosols this is a filename, which is
!            identified by a ".nc" suffix.  The file contains optical and
!            other physical properties of the aerosol.
!
! This code also identifies whether the constituent is a gas or an aerosol
! and adds that info to a structure that stores the parsed data.
!-----------------------------------------------------------------------------

    character(len=*), dimension(:), intent(in) :: specifier
    type(rad_cnst_namelist_t),   intent(inout) :: namelist_data

    ! Local variables
    integer            :: number, i, j
    integer            :: ipos, strlen
    integer            :: astat
    character(len=cs1) :: tmpstr
    character(len=1)   :: source(n_rad_cnst)
    character(len=64)  :: camname(n_rad_cnst)
    character(len=cs1) :: radname(n_rad_cnst)
    character(len=1)   :: type(n_rad_cnst)
    !-------------------------------------------------------------------------

    number = 0

    parse_loop: do i = 1, n_rad_cnst
      if ( len_trim(specifier(i)) == 0 ) then
         exit parse_loop
      endif

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(specifier(i))
      do
         strlen = len_trim(tmpstr)
         ipos = index(tmpstr, ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = tmpstr(:ipos-1) // tmpstr(ipos+1:strlen)
      end do

      ! Locate the ':' separating source from camname.
      j = index(tmpstr, ':')
      source(i) = tmpstr(:j-1)
      tmpstr = tmpstr(j+1:)

      ! locate the ':' separating camname from radname
      j = scan(tmpstr, ':')

      camname(i) = tmpstr(:j-1)
      radname(i) = tmpstr(j+1:)

      ! determine the type of constituent
      if (source(i) == 'M') then
         type(i) = 'M'
      else if (source(i) == 'B') then
         type(i) = 'B'
      else if(index(radname(i),".nc") .gt. 0) then
         type(i) = 'A'
      else
         type(i) = 'G'
      end if

      number = number+1
    end do parse_loop

    namelist_data%ncnst = number

    if (number == 0) return

    allocate(namelist_data%source (number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%source')
    allocate(namelist_data%camname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%camname')
    allocate(namelist_data%radname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%radname')
    allocate(namelist_data%type(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%type')

    namelist_data%source(:namelist_data%ncnst)  = source (:namelist_data%ncnst)
    namelist_data%camname(:namelist_data%ncnst) = camname(:namelist_data%ncnst)
    namelist_data%radname(:namelist_data%ncnst) = radname(:namelist_data%ncnst)
    namelist_data%type(:namelist_data%ncnst)    = type(:namelist_data%ncnst)

end subroutine parse_rad_specifier

subroutine rad_cnst_get_aer_mmr_by_idx(list_idx, aer_idx, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: aer_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: lchnk
   integer :: idx
   character(len=1) :: source
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = 'rad_cnst_get_aer_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   lchnk = state%lchnk

   ! Check for valid input aerosol index
   if (aer_idx < 1  .or.  aer_idx > aerlist%numaerosols) then
      write(iulog,*) subname//': aer_idx= ', aer_idx, '  numaerosols= ', aerlist%numaerosols
      call endrun(subname//': aerosol list index out of range')
   end if

   ! Get data source
   source = aerlist%aer(aer_idx)%source
   idx    = aerlist%aer(aer_idx)%idx
   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_aer_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_mmr_by_idx(list_idx, mode_idx, spec_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the modal aerosol specie from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: mode_idx    ! mode index
   integer,                     intent(in) :: spec_idx    ! index of specie in the mode
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Get data source
   if (phase == 'a') then
      source = modes%comps(m_idx)%source_mmr_a(spec_idx)
      idx    = modes%comps(m_idx)%idx_mmr_a(spec_idx)
   else if (phase == 'c') then
      source = modes%comps(m_idx)%source_mmr_c(spec_idx)
      idx    = modes%comps(m_idx)%idx_mmr_c(spec_idx)
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_mam_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_bin_mmr_by_idx(list_idx, bin_idx, spec_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the modal aerosol specie from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: bin_idx    ! mode index
   integer,                     intent(in) :: spec_idx    ! index of specie in the mode
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: s_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sa_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   s_idx = slist%idx(bin_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > bins%comps(s_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', bins%comps(s_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Get data source
   if (phase == 'a') then
      source = bins%comps(s_idx)%source_mmr_a(spec_idx)
      idx    = bins%comps(s_idx)%idx_mmr_a(spec_idx)
   else if (phase == 'c') then
      source = bins%comps(s_idx)%source_mmr_c(spec_idx)
      idx    = bins%comps(s_idx)%idx_mmr_c(spec_idx)
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_bin_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_mmr_idx(mode_idx, spec_idx, idx)

   ! Return constituent index of mam specie mass mixing ratio for aerosol modes in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   ! Arguments
   integer, intent(in)  :: mode_idx    ! mode index
   integer, intent(in)  :: spec_idx    ! index of specie in the mode
   integer, intent(out) :: idx         ! index of specie in the constituent array

   ! Local variables
   integer :: m_idx
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_mmr_idx'
   !-----------------------------------------------------------------------------

   ! assume climate list (i.e., species are in the constituent array)
   mlist => ma_list(0)

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Assume data source is interstitial since that's what's in the constituent array
   idx    = modes%comps(m_idx)%idx_mmr_a(spec_idx)

end subroutine rad_cnst_get_mam_mmr_idx

!================================================================================================

subroutine rad_cnst_get_carma_mmr_idx(bin_idx, spec_idx, idx)

   ! Return constituent index of camra species mass mixing ratio for aerosol bins in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   ! Arguments
   integer, intent(in)  :: bin_idx     ! bin index
   integer, intent(in)  :: spec_idx    ! index of specie in the bin
   integer, intent(out) :: idx         ! index of specie in the constituent array

   ! Local variables
   integer :: b_idx
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_carma_mmr_idx'
   !-----------------------------------------------------------------------------

   ! assume climate list (i.e., species are in the constituent array)
   slist => sa_list(0)

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   b_idx = slist%idx(bin_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > bins%comps(b_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', bins%comps(b_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Assume data source is interstitial since that's what's in the constituent array
   idx = bins%comps(b_idx)%idx_mmr_a(spec_idx)

end subroutine rad_cnst_get_carma_mmr_idx

!================================================================================================

subroutine rad_cnst_get_bin_mmr(list_idx, bin_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the aerosol bin from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: bin_idx     ! bin index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_mmr'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sa_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   m_idx = slist%idx(bin_idx)

   ! Get data source
   if (phase == 'a') then
      source = bins%comps(m_idx)%source_mass_a
      idx    = bins%comps(m_idx)%idx_mass_a
   else if (phase == 'c') then
      source = bins%comps(m_idx)%source_mass_c
      idx    = bins%comps(m_idx)%idx_mass_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_bin_mmr

!================================================================================================

subroutine rad_cnst_get_mode_num(list_idx, mode_idx, phase, state, pbuf, num)

   ! Return pointer to number mixing ratio for the aerosol mode from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: mode_idx    ! mode index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: num(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_num'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Get data source
   if (phase == 'a') then
      source = modes%comps(m_idx)%source_num_a
      idx    = modes%comps(m_idx)%idx_num_a
   else if (phase == 'c') then
      source = modes%comps(m_idx)%source_num_c
      idx    = modes%comps(m_idx)%idx_num_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      num => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, num)
   case ('Z')
      num => zero_cols
   end select

end subroutine rad_cnst_get_mode_num

!================================================================================================

subroutine rad_cnst_get_bin_num(list_idx, bin_idx, phase, state, pbuf, num)

   ! Return pointer to number mixing ratio for the aerosol bin from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: bin_idx     ! bin index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: num(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_num'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sa_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   m_idx = slist%idx(bin_idx)

   ! Get data source
   if (phase == 'a') then
      source = bins%comps(m_idx)%source_num_a
      idx    = bins%comps(m_idx)%idx_num_a
   else if (phase == 'c') then
      source = bins%comps(m_idx)%source_num_c
      idx    = bins%comps(m_idx)%idx_num_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      num => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, num)
   case ('Z')
      num => zero_cols
   end select

end subroutine rad_cnst_get_bin_num

!================================================================================================

subroutine rad_cnst_get_mode_num_idx(mode_idx, cnst_idx)

   ! Return constituent index of mode number mixing ratio for the aerosol mode in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   ! Arguments
   integer,  intent(in)  :: mode_idx    ! mode index
   integer,  intent(out) :: cnst_idx    ! constituent index

   ! Local variables
   integer :: m_idx
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_num'
   !-----------------------------------------------------------------------------

   ! assume climate list
   mlist => ma_list(0)

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check that source is 'A' which means the index is for the constituent array
   source = modes%comps(m_idx)%source_num_a
   if (source /= 'A') then
      write(iulog,*) subname//': source= ', source
      call endrun(subname//': requested mode number index not in constituent array')
   end if

   ! Return index in constituent array
   cnst_idx = modes%comps(m_idx)%idx_num_a

end subroutine rad_cnst_get_mode_num_idx

!================================================================================================

subroutine rad_cnst_get_bin_num_idx(bin_idx, cnst_idx)

   ! Return constituent index of bin number mixing ratio for the aerosol bin in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   ! Arguments
   integer,  intent(in)  :: bin_idx    ! bin index
   integer,  intent(out) :: cnst_idx    ! constituent index

   ! Local variables
   integer :: b_idx
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_num_idx'
   !-----------------------------------------------------------------------------

   ! assume climate list
   slist => sa_list(0)

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   b_idx = slist%idx(bin_idx)

   ! Check that source is 'A' which means the index is for the constituent array
   source = bins%comps(b_idx)%source_num_a
   if (source /= 'A') then
      write(iulog,*) subname//': source= ', source
      call endrun(subname//': requested bin number index not in constituent array')
   end if

   ! Return index in constituent array
   cnst_idx = bins%comps(b_idx)%idx_num_a

end subroutine rad_cnst_get_bin_num_idx

!================================================================================================

subroutine print_modes(modes)

   type(modes_t), intent(inout) :: modes

   integer :: i, m
   !---------------------------------------------------------------------------------------------

   write(iulog,*)' Mode Definitions'

   do m = 1, modes%nmodes

      write(iulog,*) nl//' name=',trim(modes%names(m)),'  type=',trim(modes%types(m))
      write(iulog,*) ' src_a=',trim(modes%comps(m)%source_num_a),'  num_a=',trim(modes%comps(m)%camname_num_a), &
                     ' src_c=',trim(modes%comps(m)%source_num_c),'  num_c=',trim(modes%comps(m)%camname_num_c)

      do i = 1, modes%comps(m)%nspec

         write(iulog,*) ' src_a=',trim(modes%comps(m)%source_mmr_a(i)), '  mmr_a=',trim(modes%comps(m)%camname_mmr_a(i)), &
                       '  src_c=',trim(modes%comps(m)%source_mmr_c(i)), '  mmr_c=',trim(modes%comps(m)%camname_mmr_c(i)), &
                       '  type=',trim(modes%comps(m)%type(i))
         write(iulog,*) '     prop file=', trim(modes%comps(m)%props(i))
      end do

   end do

end subroutine print_modes

!================================================================================================

subroutine print_bins(bins)

   type(bins_t), intent(inout) :: bins

   integer :: i, m
   !---------------------------------------------------------------------------------------------

   write(iulog,*)' Bin Definitions'

   do m = 1, bins%nbins

      write(iulog,*) nl//' name=',trim(bins%names(m))

      do i = 1, bins%comps(m)%nspec

         write(iulog,*) ' src_a=',trim(bins%comps(m)%source_mmr_a(i)), '  mmr_a=',trim(bins%comps(m)%camname_mmr_a(i)), &
                       '  type=',trim(bins%comps(m)%type(i))
         write(iulog,*) '     prop file=', trim(bins%comps(m)%props(i))
      end do

   end do

end subroutine print_bins

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
