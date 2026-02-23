module aerosol_definition_mod

!------------------------------------------------------------------------------------------------
!
! Structural types, data, query routines, and property-access routines
! extracted from rad_constituents.  Concrete aerosol property/state/instances
! modules can depend on aerosol_definition_mod instead of rad_constituents.
!
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use phys_prop,      only: physprop_get, physprop_get_id, ot_length
use cam_abortutils, only: endrun
use cam_logfile,    only: iulog
use constituents,   only: cnst_get_ind
use physics_buffer, only: physics_buffer_desc, pbuf_get_index
use radconstants,   only: nradgas, rad_gas_index
use spmd_utils,     only: masterproc

implicit none
private
save

integer, parameter :: cs1 = 256
integer, public, parameter :: N_DIAG = 10

! max number of strings in mode definitions
integer, parameter :: n_mode_str = 120

! max number of strings in bin definitions
integer, parameter :: n_bin_str = 640

! max number of externally mixed entities in the climate/diag lists
integer, parameter :: n_rad_cnst = N_RAD_CNST

! type to provide access to the components of a mode
type :: mode_component_t
   integer :: nspec
   ! For "source" variables below, value is:
   ! 'N' if in pbuf (non-advected)
   ! 'A' if in state (advected)
   character(len=  1) :: source_num_a  ! source of interstitial number conc field
   character(len= 32) :: camname_num_a ! name registered in pbuf or constituents for number mixing ratio of interstitial species
   character(len=  1) :: source_num_c  ! source of cloud borne number conc field
   character(len= 32) :: camname_num_c ! name registered in pbuf or constituents for number mixing ratio of cloud borne species
   character(len=  1), pointer :: source_mmr_a(:)  ! source of interstitial specie mmr fields
   character(len= 32), pointer :: camname_mmr_a(:) ! name registered in pbuf or constituents for mmr of interstitial components
   character(len=  1), pointer :: source_mmr_c(:)  ! source of cloud borne specie mmr fields
   character(len= 32), pointer :: camname_mmr_c(:) ! name registered in pbuf or constituents for mmr of cloud borne components
   character(len= 32), pointer :: type(:)          ! specie type (as used in MAM code)
   character(len=cs1), pointer :: props(:)         ! file containing specie properties
   integer          :: idx_num_a    ! index in pbuf or constituents for number mixing ratio of interstitial species
   integer          :: idx_num_c    ! index in pbuf for number mixing ratio of interstitial species
   integer, pointer :: idx_mmr_a(:) ! index in pbuf or constituents for mmr of interstitial species
   integer, pointer :: idx_mmr_c(:) ! index in pbuf for mmr of interstitial species
   integer, pointer :: idx_props(:) ! ID used to access physical properties of mode species from phys_prop module
end type mode_component_t

! type to provide access to all modes
type :: modes_t
   integer :: nmodes
   character(len= 32),     pointer :: names(:) ! names used to identify a mode in the climate/diag lists
   character(len= 32),     pointer :: types(:) ! type of mode (as used in MAM code)
   type(mode_component_t), pointer :: comps(:) ! components which define the mode
end type modes_t

type(modes_t), target :: modes  ! mode definitions

! type to provide access to the components of a bin
type :: bin_component_t
   integer :: nspec
   ! For "source" variables below, value is:
   ! 'N' if in pbuf (non-advected)
   ! 'A' if in state (advected)
   character(len=  1) :: source_num_a  ! source of interstitial number conc field
   character(len= 32) :: camname_num_a ! name registered in pbuf or constituents for number mixing ratio of interstitial species
   character(len=  1) :: source_num_c  ! source of cloud borne number conc field
   character(len= 32) :: camname_num_c ! name registered in pbuf or constituents for number mixing ratio of cloud borne species

   character(len=  1) :: source_mass_a  ! source of interstitial number conc field
   character(len= 32) :: camname_mass_a ! name registered in pbuf or constituents for number mixing ratio of interstitial species
   character(len=  1) :: source_mass_c  ! source of cloud borne number conc field
   character(len= 32) :: camname_mass_c ! name registered in pbuf or constituents for number mixing ratio of cloud borne species

   character(len=  1), pointer :: source_mmr_a(:)  ! source of interstitial mmr field
   character(len= 32), pointer :: camname_mmr_a(:) ! name registered in pbuf or constituents for mmr species
   character(len=  1), pointer :: source_mmr_c(:)  ! source of cloud borne specie mmr fields
   character(len= 32), pointer :: camname_mmr_c(:) ! name registered in pbuf or constituents for mmr of cloud borne components
   character(len= 32), pointer :: type(:)          ! species type
   character(len= 32), pointer :: morph(:)         ! species morphology
   character(len=cs1), pointer :: props(:)         ! file containing specie properties

   integer          :: idx_num_a    ! index in pbuf or constituents for number mixing ratio of interstitial species
   integer          :: idx_num_c    ! index in pbuf for number mixing ratio of cloud-borne species
   integer          :: idx_mass_a   ! index in pbuf or constituents for mass mixing ratio of interstitial species
   integer          :: idx_mass_c   ! index in pbuf for mass mixing ratio of cloud-borne species

   integer, pointer :: idx_mmr_a(:) ! index in pbuf or constituents for mmr of interstitial species
   integer, pointer :: idx_mmr_c(:) ! index in pbuf or constituents for mmr of cloud-borne species
   integer, pointer :: idx_props(:) ! ID used to access physical properties of mode species from phys_prop module
end type bin_component_t

! type to provide access to all bins
type :: bins_t
   integer :: nbins
   character(len= 32),    pointer :: names(:) ! names used to identify a mode in the climate/diag lists
   type(bin_component_t), pointer :: comps(:) ! components which define the mode
end type bins_t

type(bins_t), target :: bins  ! mode definitions

! type to provide access to the data parsed from the rad_climate and rad_diag_* strings
type :: rad_cnst_namelist_t
   integer :: ncnst
   character(len=  1), pointer :: source(:)  ! 'A' for state (advected), 'N' for pbuf (non-advected),
                                             ! 'M' for mode, 'Z' for zero
   character(len= 64), pointer :: camname(:) ! name registered in pbuf or constituents
   character(len=cs1), pointer :: radname(:) ! radname is the name as identfied in radiation,
                                             ! must be one of (rgaslist if a gas) or
                                             ! (/fullpath/filename.nc if an aerosol)
   character(len=  1), pointer :: type(:)    ! 'A' if aerosol, 'G' if gas, 'M' if mode
end type rad_cnst_namelist_t

type(rad_cnst_namelist_t) :: radcnst_namelist(0:N_DIAG) ! gas, bulk aerosol, and modal components used in
                                                        ! climate/diagnostic calculations

logical :: active_calls(0:N_DIAG)     ! active_calls(i) is true if the i-th call to radiation is
                                      ! specified.  Note that the 0th call is for the climate
                                      ! calculation which is always made.

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

! Storage for bulk aerosol components in the climate/diagnostic lists

type :: aerosol_t
   character(len=1)   :: source         ! A for state (advected), N for pbuf (non-advected), Z for zero
   character(len=64)  :: camname        ! name of constituent in physics state or buffer
   character(len=cs1) :: physprop_file  ! physprop filename
   character(len=32)  :: mass_name      ! name for mass per layer field in history output
   integer            :: idx            ! index of constituent in physics state or buffer
   integer            :: physprop_id    ! ID used to access physical properties from phys_prop module
end type aerosol_t

type :: aerlist_t
   integer                  :: numaerosols  ! number of aerosols
   character(len=2)         :: list_id      ! set to "  " for climate list, or two character integer
                                            ! (include leading zero) to identify diagnostic list
   type(aerosol_t), pointer :: aer(:)       ! dimension(numaerosols)
end type aerlist_t

type(aerlist_t), target :: aerosollist(0:N_DIAG) ! list of aerosols used in climate/diagnostic calcs

! storage for modal aerosol components in the climate/diagnostic lists

type :: modelist_t
   integer          :: nmodes              ! number of modes
   character(len=2) :: list_id             ! set to "  " for climate list, or two character integer
                                           ! (include leading zero) to identify diagnostic list
   integer,   pointer :: idx(:)            ! index of the mode in the mode definition object
   character(len=cs1), pointer :: physprop_files(:) ! physprop filename
   integer,   pointer :: idx_props(:)      ! index of the mode properties in the physprop object
end type modelist_t

type(modelist_t), target :: ma_list(0:N_DIAG) ! list of aerosol modes used in climate/diagnostic calcs

! storage for modal aerosol components in the climate/diagnostic lists

type :: binlist_t
   integer          :: nbins               ! number of bins
   character(len=2) :: list_id             ! set to "  " for climate list, or two character integer
                                           ! (include leading zero) to identify diagnostic list
   integer,   pointer :: idx(:)            ! index of the bin in the bin definition object
   character(len=cs1), pointer :: physprop_files(:) ! physprop filename
   integer,   pointer :: idx_props(:)      ! index of the bin properties in the physprop object
end type binlist_t

type(binlist_t), target :: sa_list(0:N_DIAG) ! list of aerosol bins used in climate/diagnostic calcs

! define generic interface routines
interface rad_cnst_get_info
   module procedure rad_cnst_get_info
   module procedure rad_cnst_get_info_by_mode
   module procedure rad_cnst_get_info_by_mode_spec
   module procedure rad_cnst_get_info_by_spectype
end interface

interface rad_cnst_get_aer_props
   module procedure rad_cnst_get_aer_props_by_idx
   module procedure rad_cnst_get_mam_props_by_idx
end interface

logical :: verbose = .true.
character(len=1), parameter :: nl = achar(10)

integer, parameter :: num_mode_types = 9
integer, parameter :: num_spec_types = 8
character(len=14), parameter :: mode_type_names(num_mode_types) = (/ &
   'accum         ', 'aitken        ', 'primary_carbon', 'fine_seasalt  ', &
   'fine_dust     ', 'coarse        ', 'coarse_seasalt', 'coarse_dust   ', &
   'coarse_strat  '  /)
character(len=9), parameter :: spec_type_names(num_spec_types) = (/ &
   'sulfate  ', 'ammonium ', 'nitrate  ', 'p-organic', &
   's-organic', 'black-c  ', 'seasalt  ', 'dust     '/)

integer, parameter :: num_bin_morphs  = 2
character(len=8), parameter :: bin_morph_names(num_bin_morphs) = &
     (/ 'shell   ', 'core    ' /)

! Public type declarations
public :: mode_component_t, modes_t, bin_component_t, bins_t
public :: rad_cnst_namelist_t, gas_t, gaslist_t, aerosol_t, aerlist_t
public :: modelist_t, binlist_t

! Public module data
public :: modes, bins, radcnst_namelist, active_calls
public :: gaslist, aerosollist, ma_list, sa_list
public :: mode_type_names, spec_type_names, num_mode_types, num_spec_types
public :: num_bin_morphs, bin_morph_names
public :: n_rad_cnst, n_mode_str, n_bin_str, cs1, verbose, nl

! Public routines
public :: rad_cnst_get_info
public :: rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec
public :: rad_cnst_get_mode_idx, rad_cnst_get_spec_idx
public :: rad_cnst_get_call_list
public :: rad_cnst_num_name
public :: rad_cnst_get_mode_props
public :: rad_cnst_get_aer_props
public :: rad_cnst_get_bin_props_by_idx
public :: rad_cnst_get_bin_props
public :: rad_cnst_get_aer_idx
public :: init_mode_comps, init_bin_comps
public :: list_init1, list_init2
public :: get_cam_idx

!==============================================================================
contains
!==============================================================================

function rad_cnst_num_name(list_idx, spc_name_in, num_name_out, mode_out, spec_out ) result(found)

  ! for a given species name spc_name_in return (optionals):
  !   num_name_out -- corresponding number density species name
  !   mode_out -- corresponding mode number
  !   spec_out -- corresponding species number within the mode

  integer,         intent(in) :: list_idx ! index of the climate or a diagnostic list
  character(len=*),intent(in) :: spc_name_in
  character(len=*),intent(out):: num_name_out
  integer,optional,intent(out):: mode_out
  integer,optional,intent(out):: spec_out

  logical :: found

  ! Local variables
  type(modelist_t), pointer :: m_list ! local pointer to mode list of interest
  integer :: n,m, mm
  integer :: nmodes
  integer :: nspecs
  character(len= 32) :: spec_name

  found = .false.

  m_list => ma_list(list_idx)
  nmodes = m_list%nmodes

  do n = 1,nmodes
     mm = m_list%idx(n)
     nspecs = modes%comps(mm)%nspec
     do m = 1,nspecs
        spec_name = modes%comps(mm)%camname_mmr_a(m)
        if (spc_name_in == spec_name) then
           num_name_out = modes%comps(mm)%camname_num_a
           found = .true.
           if (present(mode_out)) then
              mode_out = n
           endif
           if (present(spec_out)) then
              spec_out = m
           endif
           return
        endif
     enddo
  enddo

  return

end function

!================================================================================================

subroutine rad_cnst_get_info(list_idx, gasnames, aernames, &
                             use_data_o3, ngas, naero, nmodes, nbins)

   ! Return info about gas and aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=64), optional, intent(out) :: gasnames(:)
   character(len=64), optional, intent(out) :: aernames(:)
   logical,           optional, intent(out) :: use_data_o3
   integer,           optional, intent(out) :: naero
   integer,           optional, intent(out) :: ngas
   integer,           optional, intent(out) :: nmodes
   integer,           optional, intent(out) :: nbins

   ! Local variables
   type(gaslist_t),  pointer :: g_list ! local pointer to gas list of interest
   type(aerlist_t),  pointer :: a_list ! local pointer to aerosol list of interest
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest
   type(binlist_t),  pointer :: s_list ! local pointer to bin list of interest

   integer          :: i
   integer          :: arrlen  ! length of assumed shape array
   integer          :: gaslen  ! length of assumed shape array
   integer          :: igas    ! index of a gas in the gas list
   character(len=1) :: source  ! A for state, N for pbuf, Z for zero

   character(len=*), parameter :: subname = 'rad_cnst_get_info'
   !-----------------------------------------------------------------------------

   g_list => gaslist(list_idx)
   a_list => aerosollist(list_idx)
   m_list => ma_list(list_idx)
   s_list => sa_list(list_idx)

   ! number of bulk aerosols in list
   if (present(naero)) then
      naero = a_list%numaerosols
   endif

   ! number of aerosol modes in list
   if (present(nmodes)) then
      nmodes = m_list%nmodes
   endif

   ! number of aerosol bins in list
   if (present(nbins)) then
      nbins = s_list%nbins
   endif

   ! number of gases in list
   if (present(ngas)) then
      ngas = g_list%ngas
   endif

   ! names of aerosols in list
   if (present(aernames)) then

      ! check that output array is long enough
      arrlen = size(aernames)
      if (arrlen < a_list%numaerosols) then
         write(iulog,*) subname//': ERROR: naero=', a_list%numaerosols, '  arrlen=', arrlen
         call endrun(subname//': ERROR: aernames too short')
      end if

      do i = 1, a_list%numaerosols
         aernames(i) = a_list%aer(i)%camname
      end do

   end if

   ! names of gas in list
   if (present(gasnames)) then

      ! check that output array is long enough
      gaslen = size(gasnames)
      if (gaslen < g_list%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', g_list%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gasnames too short')
      end if

      do i = 1, g_list%ngas
         gasnames(i) = g_list%gas(i)%camname
      end do

   end if

   ! Does the climate calculation use data ozone?
   if (present(use_data_o3)) then

      ! get index of O3 in gas list
      igas = rad_gas_index('O3')

      ! Get data source
      source = g_list%gas(igas)%source

      use_data_o3 = .false.
      if (source == 'N') use_data_o3 = .true.
   endif

end subroutine rad_cnst_get_info

!================================================================================================

subroutine rad_cnst_get_info_by_mode(list_idx, m_idx, &
   mode_type, num_name, num_name_cw, nspec)

   ! Return info about modal aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of mode in the specified list
   character(len=32), optional, intent(out) :: mode_type   ! type of mode (as used in MAM code)
   character(len=32), optional, intent(out) :: num_name    ! name of interstitial number mixing ratio
   character(len=32), optional, intent(out) :: num_name_cw ! name of cloud borne number mixing ratio
   integer,           optional, intent(out) :: nspec       ! number of species in the mode

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: nmodes
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_mode'
   !-----------------------------------------------------------------------------

   m_list => ma_list(list_idx)

   ! check for valid mode index
   nmodes = m_list%nmodes
   if (m_idx < 1 .or. m_idx > nmodes) then
      write(iulog,*) subname//': ERROR - invalid mode index: ', m_idx
      call endrun(subname//': ERROR - invalid mode index')
   end if

   ! get index into the mode definition object
   mm = m_list%idx(m_idx)

   ! mode type
   if (present(mode_type)) then
      mode_type = modes%types(mm)
   endif

   ! number of species in the mode
   if (present(nspec)) then
      nspec = modes%comps(mm)%nspec
   endif

   ! name of interstitial number mixing ratio
   if (present(num_name)) then
      num_name = modes%comps(mm)%camname_num_a
   endif

   ! name of cloud borne number mixing ratio
   if (present(num_name_cw)) then
      num_name_cw = modes%comps(mm)%camname_num_c
   endif

end subroutine rad_cnst_get_info_by_mode

!================================================================================================

subroutine rad_cnst_get_info_by_bin(list_idx, m_idx, &
   bin_name, num_name, num_name_cw, mmr_name, mmr_name_cw, nspec)

   ! Return info about CARMA aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of bin in the specified list
   character(len=*),  optional, intent(out) :: bin_name
   character(len=32), optional, intent(out) :: num_name    ! name of interstitial number mixing ratio
   character(len=32), optional, intent(out) :: num_name_cw ! name of cloud borne number mixing ratio
   character(len=32), optional, intent(out) :: mmr_name    ! name of interstitial mass mixing ratio
   character(len=32), optional, intent(out) :: mmr_name_cw ! name of cloud borne mass mixing ratio
   integer,           optional, intent(out) :: nspec       ! number of species in the mode

   ! Local variables
   type(binlist_t), pointer :: s_list ! local pointer to mode list of interest

   integer          :: nbins
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_bin'
   !-----------------------------------------------------------------------------

   s_list => sa_list(list_idx)

   ! check for valid mode index
   nbins = s_list%nbins
   if (m_idx < 1 .or. m_idx > nbins) then
      write(iulog,*) subname//': ERROR - invalid bin index: ', m_idx
      call endrun(subname//': ERROR - invalid bin index')
   end if

   ! get index into the mode definition object
   mm = s_list%idx(m_idx)

   ! number of species in the mode
   if (present(nspec)) then
      nspec = bins%comps(mm)%nspec
   endif

   ! bin name
   if (present(bin_name)) then
      bin_name = bins%names(m_idx)
   end if

   ! name of interstitial number mixing ratio
   if (present(num_name)) then
      num_name = bins%comps(mm)%camname_num_a
   endif

   ! name of cloud borne number mixing ratio
   if (present(num_name_cw)) then
      num_name_cw = bins%comps(mm)%camname_num_c
   endif

   ! name of interstitial mass mixing ratio
   if (present(mmr_name)) then
      mmr_name = bins%comps(mm)%camname_mass_a
   endif

   ! name of cloud borne mass mixing ratio
   if (present(mmr_name_cw)) then
      mmr_name_cw = bins%comps(mm)%camname_mass_c
   endif

end subroutine rad_cnst_get_info_by_bin

!================================================================================================
subroutine rad_cnst_get_info_by_bin_spec(list_idx, m_idx, s_idx, &
   spec_type, spec_morph, spec_name, spec_name_cw)

   ! Return info about CARMA aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of bin in the specified list
   integer,                     intent(in)  :: s_idx       ! index of species in the specified mode
   character(len=32), optional, intent(out) :: spec_type   ! type of species
   character(len=32), optional, intent(out) :: spec_morph  ! type of species
   character(len=32), optional, intent(out) :: spec_name   ! name of interstitial species
   character(len=32), optional, intent(out) :: spec_name_cw ! name of cloud borne species

   ! Local variables
   type(binlist_t), pointer :: s_list ! local pointer to mode list of interest
   integer          :: nbins,  nspec
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_bin_spec'
   !-----------------------------------------------------------------------------

   s_list => sa_list(list_idx)

   ! check for valid mode index
   nbins = s_list%nbins
   if (m_idx < 1 .or. m_idx > nbins) then
      write(iulog,*) subname//': ERROR - invalid bin index: ', m_idx
      call endrun(subname//': ERROR - invalid bin index')
   end if

   ! get index into the mode definition object
   mm = s_list%idx(m_idx)

   ! check for valid species index
   nspec = bins%comps(mm)%nspec
   if (s_idx < 1 .or. s_idx > nspec) then
      write(iulog,*) subname//': ERROR - invalid specie index: ', s_idx
      call endrun(subname//': ERROR - invalid specie index')
   end if

   if (present(spec_type)) then
      spec_type = bins%comps(mm)%type(s_idx)
   endif
   if (present(spec_morph)) then
      spec_morph = bins%comps(mm)%morph(s_idx)
   endif
   if (present(spec_name)) then
      spec_name = bins%comps(mm)%camname_mmr_a(s_idx)
   endif
   if (present(spec_name_cw)) then
      spec_name_cw = bins%comps(mm)%camname_mmr_c(s_idx)
   endif

end subroutine rad_cnst_get_info_by_bin_spec

!================================================================================================
subroutine rad_cnst_get_info_by_mode_spec(list_idx, m_idx, s_idx, &
   spec_type, spec_name, spec_name_cw)

   ! Return info about modal aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of mode in the specified list
   integer,                     intent(in)  :: s_idx       ! index of specie in the specified mode
   character(len=32), optional, intent(out) :: spec_type   ! type of specie
   character(len=32), optional, intent(out) :: spec_name   ! name of interstitial specie
   character(len=32), optional, intent(out) :: spec_name_cw ! name of cloud borne specie

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: nmodes
   integer          :: nspec
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_mode_spec'
   !-----------------------------------------------------------------------------

   m_list => ma_list(list_idx)

   ! check for valid mode index
   nmodes = m_list%nmodes
   if (m_idx < 1 .or. m_idx > nmodes) then
      write(iulog,*) subname//': ERROR - invalid mode index: ', m_idx
      call endrun(subname//': ERROR - invalid mode index')
   end if

   ! get index into the mode definition object
   mm = m_list%idx(m_idx)

   ! check for valid specie index
   nspec = modes%comps(mm)%nspec
   if (s_idx < 1 .or. s_idx > nspec) then
      write(iulog,*) subname//': ERROR - invalid specie index: ', s_idx
      call endrun(subname//': ERROR - invalid specie index')
   end if

   ! specie type
   if (present(spec_type)) then
      spec_type = modes%comps(mm)%type(s_idx)
   endif

   ! interstitial specie name
   if (present(spec_name)) then
      spec_name = modes%comps(mm)%camname_mmr_a(s_idx)
   endif

   ! cloud borne specie name
   if (present(spec_name_cw)) then
      spec_name_cw = modes%comps(mm)%camname_mmr_c(s_idx)
   endif

end subroutine rad_cnst_get_info_by_mode_spec

!================================================================================================

subroutine rad_cnst_get_info_by_spectype(list_idx, spectype, mode_idx, spec_idx)

   ! Return info about modes in the specified climate/diagnostics list

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),            intent(in)  :: spectype    ! species type
   integer,           optional, intent(out) :: mode_idx    ! index of a mode that contains a specie of spectype
   integer,           optional, intent(out) :: spec_idx    ! index of the species of spectype

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer  :: i, nmodes, m_idx, nspec, ispec
   logical  :: found_spectype

   character(len=*), parameter :: subname = 'rad_cnst_get_info_by_spectype'
   !-----------------------------------------------------------------------------

   m_list => ma_list(list_idx)

   ! number of modes in specified list
   nmodes = m_list%nmodes

   ! loop through modes in specified climate/diagnostic list
   found_spectype = .false.
   do i = 1, nmodes

      ! get index of the mode in the definition object
      m_idx = m_list%idx(i)

      ! number of species in the mode
      nspec = modes%comps(m_idx)%nspec

      ! loop through species looking for spectype
      do ispec = 1, nspec

         if (trim(modes%comps(m_idx)%type(ispec)) == trim(spectype)) then
            if (present(mode_idx)) mode_idx = i
            if (present(spec_idx)) spec_idx = ispec
            found_spectype = .true.
            exit
         end if
      end do

      if (found_spectype) exit
   end do

   if (.not. found_spectype) then
      if (present(mode_idx)) mode_idx = -1
      if (present(spec_idx)) spec_idx = -1
   end if

end subroutine rad_cnst_get_info_by_spectype

!================================================================================================

function rad_cnst_get_mode_idx(list_idx, mode_type) result(mode_idx)

   ! Return mode index of the specified type in the specified climate/diagnostics list.
   ! Return -1 if not found.

   ! Arguments
   integer,           intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),  intent(in)  :: mode_type   ! mode type

   ! Return value
   integer                        :: mode_idx    ! mode index

   ! Local variables
   type(modelist_t), pointer :: m_list

   integer  :: i, nmodes, m_idx

   character(len=*), parameter :: subname = 'rad_cnst_get_mode_idx'
   !-----------------------------------------------------------------------------

   ! if mode type not found return -1
   mode_idx = -1

   ! specified mode list
   m_list => ma_list(list_idx)

   ! number of modes in specified list
   nmodes = m_list%nmodes

   ! loop through modes in specified climate/diagnostic list
   do i = 1, nmodes

      ! get index of the mode in the definition object
      m_idx = m_list%idx(i)

      ! look in mode definition object (modes) for the mode types
      if (trim(modes%types(m_idx)) == trim(mode_type)) then
         mode_idx = i
         exit
      end if
   end do

end function rad_cnst_get_mode_idx

!================================================================================================

function rad_cnst_get_spec_idx(list_idx, mode_idx, spec_type) result(spec_idx)

   ! Return specie index of the specified type in the specified mode of the specified
   ! climate/diagnostics list.  Return -1 if not found.

   ! Arguments
   integer,           intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,           intent(in)  :: mode_idx    ! mode index
   character(len=*),  intent(in)  :: spec_type   ! specie type

   ! Return value
   integer                        :: spec_idx    ! specie index

   ! Local variables
   type(modelist_t),       pointer :: m_list
   type(mode_component_t), pointer :: mode_comps

   integer  :: i, m_idx, nspec

   character(len=*), parameter :: subname = 'rad_cnst_get_spec_idx'
   !-----------------------------------------------------------------------------

   ! if specie type not found return -1
   spec_idx = -1

   ! modes in specified list
   m_list => ma_list(list_idx)

   ! get index of the specified mode in the definition object
   m_idx = m_list%idx(mode_idx)

   ! object containing the components of the mode
   mode_comps => modes%comps(m_idx)

   ! number of species in specified mode
   nspec = mode_comps%nspec

   ! loop through species in specified mode
   do i = 1, nspec

      ! look in mode definition object (modes) for the mode types
      if (trim(mode_comps%type(i)) == trim(spec_type)) then
         spec_idx = i
         exit
      end if
   end do

end function rad_cnst_get_spec_idx

!================================================================================================

subroutine rad_cnst_get_call_list(call_list)

   ! Return info about which climate/diagnostic calculations are requested

   ! Arguments
   logical, intent(out) :: call_list(0:N_DIAG)
   !-----------------------------------------------------------------------------

   call_list(:) = active_calls(:)

end subroutine rad_cnst_get_call_list

!================================================================================================

subroutine init_mode_comps(modes)

   ! Initialize the mode definitions by looking up the relevent indices in the
   ! constituent and pbuf arrays, and getting the physprop IDs

   ! Arguments
   type(modes_t), intent(inout) :: modes

   ! Local variables
   integer :: m, ispec, nspec

   character(len=*), parameter :: routine = 'init_mode_comps'
   !-----------------------------------------------------------------------------

   do m = 1, modes%nmodes

      ! indices for number mixing ratio components
      modes%comps(m)%idx_num_a = get_cam_idx(modes%comps(m)%source_num_a, modes%comps(m)%camname_num_a, routine)
      modes%comps(m)%idx_num_c = get_cam_idx(modes%comps(m)%source_num_c, modes%comps(m)%camname_num_c, routine)

      ! allocate memory for species
      nspec = modes%comps(m)%nspec
      allocate( &
         modes%comps(m)%idx_mmr_a(nspec), &
         modes%comps(m)%idx_mmr_c(nspec), &
         modes%comps(m)%idx_props(nspec)  )

      do ispec = 1, nspec

         ! indices for species mixing ratio components
         modes%comps(m)%idx_mmr_a(ispec) = get_cam_idx(modes%comps(m)%source_mmr_a(ispec), &
                                                   modes%comps(m)%camname_mmr_a(ispec), routine)
         modes%comps(m)%idx_mmr_c(ispec) = get_cam_idx(modes%comps(m)%source_mmr_c(ispec), &
                                                   modes%comps(m)%camname_mmr_c(ispec), routine)

         ! get physprop ID
         modes%comps(m)%idx_props(ispec) = physprop_get_id(modes%comps(m)%props(ispec))
         if (modes%comps(m)%idx_props(ispec) == -1) then
            call endrun(routine//' : ERROR idx not found for '//trim(modes%comps(m)%props(ispec)))
         end if

      end do

   end do

end subroutine init_mode_comps

!================================================================================================

subroutine init_bin_comps(bins)

   ! Initialize the mode definitions by looking up the relevent indices in the
   ! constituent and pbuf arrays, and getting the physprop IDs

   ! Arguments
   type(bins_t), intent(inout) :: bins

   ! Local variables
   integer :: m, ispec, nspec

   character(len=*), parameter :: routine = 'init_bin_comps'
   !-----------------------------------------------------------------------------

   do m = 1, bins%nbins

      ! indices for number mixing ratio components
      bins%comps(m)%idx_num_a = get_cam_idx(bins%comps(m)%source_num_a, bins%comps(m)%camname_num_a, routine)
      bins%comps(m)%idx_num_c = get_cam_idx(bins%comps(m)%source_num_c, bins%comps(m)%camname_num_c, routine)
      if ( bins%comps(m)%source_mass_a /= 'NOTSET' .and. bins%comps(m)%camname_mass_a /= 'NOTSET' ) then
         bins%comps(m)%idx_mass_a = get_cam_idx(bins%comps(m)%source_mass_a, bins%comps(m)%camname_mass_a, routine)
      endif
      if ( bins%comps(m)%source_mass_c /= 'NOTSET' .and. bins%comps(m)%camname_mass_c /= 'NOTSET' ) then
         bins%comps(m)%idx_mass_c = get_cam_idx(bins%comps(m)%source_mass_c, bins%comps(m)%camname_mass_c, routine)
      endif

      ! allocate memory for species
      nspec = bins%comps(m)%nspec
      allocate( &
         bins%comps(m)%idx_mmr_a(nspec), &
         bins%comps(m)%idx_mmr_c(nspec), &
         bins%comps(m)%idx_props(nspec)  )

      do ispec = 1, nspec

         ! indices for species mixing ratio components
         bins%comps(m)%idx_mmr_a(ispec) = get_cam_idx(bins%comps(m)%source_mmr_a(ispec), &
                                                   bins%comps(m)%camname_mmr_a(ispec), routine)
         bins%comps(m)%idx_mmr_c(ispec) = get_cam_idx(bins%comps(m)%source_mmr_c(ispec), &
                                                   bins%comps(m)%camname_mmr_c(ispec), routine)

         ! get physprop ID
         bins%comps(m)%idx_props(ispec) = physprop_get_id(bins%comps(m)%props(ispec))
         if (bins%comps(m)%idx_props(ispec) == -1) then
            call endrun(routine//' : ERROR idx not found for '//trim(bins%comps(m)%props(ispec)))
         end if

      end do

   end do

end subroutine init_bin_comps

!================================================================================================

integer function get_cam_idx(source, name, routine)

   ! get index of name in internal CAM array; either the constituent array
   ! or the physics buffer

   character(len=*), intent(in) :: source
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: routine  ! name of calling routine

   integer :: idx
   integer :: errcode
   !-----------------------------------------------------------------------------

   if (source(1:1) == 'N') then

      idx = pbuf_get_index(trim(name),errcode)
      if (errcode < 0) then
         call endrun(routine//' ERROR: cannot find physics buffer field '//trim(name))
      end if

   else if (source(1:1) == 'A') then

      call cnst_get_ind(trim(name), idx, abort=.false.)
      if (idx < 0) then
         call endrun(routine//' ERROR: cannot find constituent field '//trim(name))
      end if

   else if (source(1:1) == 'Z') then

      idx = -1

   else

      call endrun(routine//' ERROR: invalid source for specie '//trim(name))

   end if

   get_cam_idx = idx

end function get_cam_idx

!================================================================================================

subroutine list_init1(namelist, gaslist, aerlist, ma_list, sa_list)

   ! Initialize the gas and bulk and modal aerosol lists with the
   ! entities specified in the climate or diagnostic lists.

   ! This first phase initialization just sets the information that
   ! is available at the time the namelist is read.

   type(rad_cnst_namelist_t), intent(in) :: namelist ! parsed namelist input for climate or diagnostic lists

   type(gaslist_t),        intent(inout) :: gaslist
   type(aerlist_t),        intent(inout) :: aerlist
   type(modelist_t),       intent(inout) :: ma_list
   type(binlist_t),        intent(inout) :: sa_list

   ! Local variables
   integer :: ii, m, naero, nmodes, nbins
   integer :: igas, ba_idx, ma_idx, sa_idx
   integer :: istat
   character(len=*), parameter :: routine = 'list_init1'
   !-----------------------------------------------------------------------------

   ! nradgas is set by the radiative transfer code
   gaslist%ngas = nradgas

   ! Determine the number of bulk aerosols and aerosol modes in the list
   naero = 0
   nmodes = 0
   nbins = 0
   do ii = 1, namelist%ncnst
      if (trim(namelist%type(ii)) == 'A') naero  = naero + 1
      if (trim(namelist%type(ii)) == 'M') nmodes = nmodes + 1
      if (trim(namelist%type(ii)) == 'B') nbins = nbins + 1
   end do
   aerlist%numaerosols = naero
   ma_list%nmodes      = nmodes
   sa_list%nbins       = nbins

   ! allocate storage for the aerosol, gas, and mode lists
   allocate( &
      aerlist%aer(aerlist%numaerosols),      &
      gaslist%gas(gaslist%ngas),             &
      ma_list%idx(ma_list%nmodes),           &
      ma_list%physprop_files(ma_list%nmodes), &
      ma_list%idx_props(ma_list%nmodes),     &
      sa_list%idx(sa_list%nbins),           &
      sa_list%physprop_files(sa_list%nbins), &
      sa_list%idx_props(sa_list%nbins),     &
      stat=istat)
   if (istat /= 0) call endrun(routine//': allocate ERROR; aero and gas list components')

   if (masterproc .and. verbose) then
      if (len_trim(gaslist%list_id) == 0) then
         write(iulog,*) nl//' '//routine//': namelist input for climate list'
      else
         write(iulog,*) nl//' '//routine//': namelist input for diagnostic list:'//gaslist%list_id
      end if
   end if

   ! Loop over the radiatively active components specified in the namelist
   ba_idx = 0
   ma_idx = 0
   sa_idx = 0
   do ii = 1, namelist%ncnst

      if (masterproc .and. verbose) &
         write(iulog,*) "  rad namelist spec: "// trim(namelist%source(ii)) &
         //":"//trim(namelist%camname(ii))//":"//trim(namelist%radname(ii))

      ! Check that the source specifier is legal.
      if (namelist%source(ii) /= 'A' .and. namelist%source(ii) /= 'M' .and. &
          namelist%source(ii) /= 'N' .and. namelist%source(ii) /= 'Z' .and. &
          namelist%source(ii) /= 'B' ) then
         call endrun(routine//": source must either be A, B, M, N or Z:"//&
                     " illegal specifier in namelist input: "//namelist%source(ii))
      end if

      ! Add component to appropriate list (gas, modal or bulk aerosol)
      if (namelist%type(ii) == 'A') then

         ! Add to bulk aerosol list
         ba_idx = ba_idx + 1

         aerlist%aer(ba_idx)%source        = namelist%source(ii)
         aerlist%aer(ba_idx)%camname       = namelist%camname(ii)
         aerlist%aer(ba_idx)%physprop_file = namelist%radname(ii)

      else if (namelist%type(ii) == 'M') then

         ! Add to modal aerosol list
         ma_idx = ma_idx + 1

         ! Look through the mode definitions for the name of the specified mode.  The
         ! index into the modes object all the information relevent to the mode definition.
         ma_list%idx(ma_idx) = -1
         do m = 1, modes%nmodes
            if (trim(namelist%camname(ii)) == trim(modes%names(m))) then
               ma_list%idx(ma_idx) = m
               exit
            end if
         end do
         if (ma_list%idx(ma_idx) == -1) &
            call endrun(routine//' ERROR cannot find mode name '//trim(namelist%camname(ii)))

         ! Also save the name of the physprop file
         ma_list%physprop_files(ma_idx) = namelist%radname(ii)

      else if (namelist%type(ii) == 'B') then

         ! Add to modal aerosol list
         sa_idx = sa_idx + 1

         ! Look through the bin definitions for the name of the specified bin.  The
         ! index into the modes object all the information relevent to the mode definition.
         sa_list%idx(sa_idx) = -1
         do m = 1, bins%nbins
            if (trim(namelist%camname(ii)) == trim(bins%names(m))) then
               sa_list%idx(sa_idx) = m
               exit
            end if
         end do
         if (sa_list%idx(sa_idx) == -1) &
            call endrun(routine//' ERROR cannot find bin name '//trim(namelist%camname(ii)))

         ! Also save the name of the physprop file
         sa_list%physprop_files(sa_idx) = namelist%radname(ii)

      else

         ! Add to gas list

         ! The radiative transfer code requires the input of a specific set of gases
         ! which is hardwired into the code.  The CAM interface to the RT code uses
         ! the names in the radconstants module to refer to these gases.  The user
         ! interface (namelist) also uses these names to identify the gases treated
         ! by the RT code.  We use the index order set in radconstants for convenience
         ! only.

         ! First check that the gas name specified by the user is allowed.
         ! rad_gas_index will abort on illegal names.
         igas = rad_gas_index(namelist%radname(ii))

         ! Set values in the igas index
         gaslist%gas(igas)%source  = namelist%source(ii)
         gaslist%gas(igas)%camname = namelist%camname(ii)

      end if
   end do

end subroutine list_init1

!================================================================================================

subroutine list_init2(gaslist, aerlist, ma_list, sa_list)

   ! Final initialization phase gets the component indices in the constituent array
   ! and the physics buffer, and indices into physprop module.

   type(gaslist_t),        intent(inout) :: gaslist
   type(aerlist_t),        intent(inout) :: aerlist
   type(modelist_t),       intent(inout) :: ma_list
   type(binlist_t),        intent(inout) :: sa_list

   ! Local variables
   integer :: i
   character(len=*), parameter :: routine = 'list_init2'
   !-----------------------------------------------------------------------------

   ! Loop over gases
   do i = 1, gaslist%ngas

      ! locate the specie mixing ratio in the pbuf or state
      gaslist%gas(i)%idx = get_cam_idx(gaslist%gas(i)%source, gaslist%gas(i)%camname, routine)

   end do

   ! Loop over bulk aerosols
   do i = 1, aerlist%numaerosols

      ! locate the specie mixing ratio in the pbuf or state
      aerlist%aer(i)%idx = get_cam_idx(aerlist%aer(i)%source, aerlist%aer(i)%camname, routine)

      ! get the physprop_id from the phys_prop module
      aerlist%aer(i)%physprop_id = physprop_get_id(aerlist%aer(i)%physprop_file)

   end do

   ! Loop over modes
   do i = 1, ma_list%nmodes

      ! get the physprop_id from the phys_prop module
      ma_list%idx_props(i) = physprop_get_id(ma_list%physprop_files(i))

   end do

   ! Loop over bins
   do i = 1, sa_list%nbins

      ! get the physprop_id from the phys_prop module
      sa_list%idx_props(i) = physprop_get_id(sa_list%physprop_files(i))

   end do

end subroutine list_init2

!================================================================================================

integer function rad_cnst_get_aer_idx(list_idx, aer_name)

   ! Return the index of aerosol aer_name in the list specified by list_idx.

    ! Arguments
   integer,             intent(in) :: list_idx    ! 0 for climate list, 1-N_DIAG for diagnostic lists
   character(len=*),    intent(in) :: aer_name    ! aerosol name (in state or pbuf)

   ! Local variables
   integer :: i, aer_idx
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = "rad_cnst_get_aer_idx"
   !-------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Get index in aerosol list for requested name
   aer_idx = -1
   do i = 1, aerlist%numaerosols
      if (trim(aer_name) == trim(aerlist%aer(i)%camname)) then
         aer_idx = i
         exit
      end if
   end do

   if (aer_idx == -1) call endrun(subname//": ERROR - name not found")

   rad_cnst_get_aer_idx = aer_idx

end function rad_cnst_get_aer_idx

!================================================================================================

subroutine rad_cnst_get_aer_props_by_idx(list_idx, &
   aer_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, num_to_mass_aer)

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in)  :: list_idx ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: aer_idx  ! index of the aerosol
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)
   character(len=20), optional, intent(out) :: aername
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer

   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   ! Local variables
   integer :: id
   character(len=*), parameter :: subname = 'rad_cnst_get_aer_props_by_idx'
   type(aerlist_t), pointer :: aerlist
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => aerosollist(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   if (aer_idx < 1 .or. aer_idx > aerlist%numaerosols) then
      write(iulog,*) subname//': aerosol list index out of range: ', aer_idx ,' list index: ',list_idx
      call endrun(subname//': aer_idx out of range')
   end if

   id = aerlist%aer(aer_idx)%physprop_id

   if (present(opticstype))        call physprop_get(id, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(id, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(id, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(id, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(id, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(id, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(id, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(id, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(id, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(id, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(id, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(id, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(id, refindex_aer_lw=refindex_aer_lw)

   if (present(aername))           call physprop_get(id, aername=aername)
   if (present(density_aer))       call physprop_get(id, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(id, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(id, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(id, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(id, num_to_mass_aer=num_to_mass_aer)

   if (present(r_lw_abs))          call physprop_get(id, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(id, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(id, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(id, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(id, mu=mu)

end subroutine rad_cnst_get_aer_props_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_props_by_idx(list_idx, &
   mode_idx, spec_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, &
   num_to_mass_aer, spectype)

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: mode_idx  ! mode index
   integer,                     intent(in)  :: spec_idx  ! index of specie in the mode
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)

   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   character(len=20), optional, intent(out) :: aername
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer
   character(len=32), optional, intent(out) :: spectype

   ! Local variables
   integer :: m_idx, id
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_props_by_idx'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
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

   id = modes%comps(m_idx)%idx_props(spec_idx)

   if (present(opticstype))        call physprop_get(id, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(id, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(id, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(id, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(id, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(id, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(id, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(id, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(id, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(id, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(id, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(id, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(id, refindex_aer_lw=refindex_aer_lw)

   if (present(r_lw_abs))          call physprop_get(id, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(id, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(id, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(id, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(id, mu=mu)

   if (present(aername))           call physprop_get(id, aername=aername)
   if (present(density_aer))       call physprop_get(id, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(id, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(id, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(id, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(id, num_to_mass_aer=num_to_mass_aer)

   if (present(spectype)) spectype = modes%comps(m_idx)%type(spec_idx)

end subroutine rad_cnst_get_mam_props_by_idx

!================================================================================================

subroutine rad_cnst_get_bin_props_by_idx(list_idx, &
   bin_idx, spec_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, &
   num_to_mass_aer, spectype, specmorph)

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: bin_idx  ! mode index
   integer,                     intent(in)  :: spec_idx  ! index of specie in the mode
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)

   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   character(len=20), optional, intent(out) :: aername
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer
   character(len=32), optional, intent(out) :: spectype
   character(len=32), optional, intent(out) :: specmorph

   ! Local variables
   integer :: m_idx, id
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_props_by_idx'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sa_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = slist%idx(bin_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > bins%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', bins%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   id = bins%comps(m_idx)%idx_props(spec_idx)

   if (present(opticstype))        call physprop_get(id, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(id, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(id, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(id, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(id, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(id, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(id, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(id, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(id, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(id, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(id, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(id, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(id, refindex_aer_lw=refindex_aer_lw)

   if (present(r_lw_abs))          call physprop_get(id, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(id, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(id, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(id, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(id, mu=mu)

   if (present(aername))           call physprop_get(id, aername=aername)
   if (present(density_aer))       call physprop_get(id, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(id, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(id, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(id, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(id, num_to_mass_aer=num_to_mass_aer)

   if (present(spectype)) spectype = bins%comps(m_idx)%type(spec_idx)
   if (present(specmorph)) specmorph = bins%comps(m_idx)%morph(spec_idx)

end subroutine rad_cnst_get_bin_props_by_idx

!================================================================================================

subroutine rad_cnst_get_mode_props(list_idx, mode_idx, opticstype, &
   extpsw, abspsw, asmpsw, absplw, refrtabsw, &
   refitabsw, refrtablw, refitablw, ncoef, prefr, &
   prefi, sigmag, dgnum, dgnumlo, dgnumhi, &
   rhcrystal, rhdeliques)

   ! Return requested properties for the mode from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,             intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,             intent(in)  :: mode_idx  ! mode index

   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),  optional, pointer     :: extpsw(:,:,:,:)
   real(r8),  optional, pointer     :: abspsw(:,:,:,:)
   real(r8),  optional, pointer     :: asmpsw(:,:,:,:)
   real(r8),  optional, pointer     :: absplw(:,:,:,:)
   real(r8),  optional, pointer     :: refrtabsw(:,:)
   real(r8),  optional, pointer     :: refitabsw(:,:)
   real(r8),  optional, pointer     :: refrtablw(:,:)
   real(r8),  optional, pointer     :: refitablw(:,:)
   integer,   optional, intent(out) :: ncoef
   integer,   optional, intent(out) :: prefr
   integer,   optional, intent(out) :: prefi
   real(r8),  optional, intent(out) :: sigmag
   real(r8),  optional, intent(out) :: dgnum
   real(r8),  optional, intent(out) :: dgnumlo
   real(r8),  optional, intent(out) :: dgnumhi
   real(r8),  optional, intent(out) :: rhcrystal
   real(r8),  optional, intent(out) :: rhdeliques

   ! Local variables
   integer :: id
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_props'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => ma_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the physprop index for the requested mode
   id = mlist%idx_props(mode_idx)

   if (present(opticstype))  call physprop_get(id, opticstype=opticstype)
   if (present(extpsw))      call physprop_get(id, extpsw=extpsw)
   if (present(abspsw))      call physprop_get(id, abspsw=abspsw)
   if (present(asmpsw))      call physprop_get(id, asmpsw=asmpsw)
   if (present(absplw))      call physprop_get(id, absplw=absplw)

   if (present(refrtabsw))   call physprop_get(id, refrtabsw=refrtabsw)
   if (present(refitabsw))   call physprop_get(id, refitabsw=refitabsw)
   if (present(refrtablw))   call physprop_get(id, refrtablw=refrtablw)
   if (present(refitablw))   call physprop_get(id, refitablw=refitablw)

   if (present(ncoef))       call physprop_get(id, ncoef=ncoef)
   if (present(prefr))       call physprop_get(id, prefr=prefr)
   if (present(prefi))       call physprop_get(id, prefi=prefi)
   if (present(sigmag))      call physprop_get(id, sigmag=sigmag)
   if (present(dgnum))       call physprop_get(id, dgnum=dgnum)
   if (present(dgnumlo))     call physprop_get(id, dgnumlo=dgnumlo)
   if (present(dgnumhi))     call physprop_get(id, dgnumhi=dgnumhi)
   if (present(rhcrystal))   call physprop_get(id, rhcrystal=rhcrystal)
   if (present(rhdeliques))  call physprop_get(id, rhdeliques=rhdeliques)

end subroutine rad_cnst_get_mode_props

!================================================================================================

subroutine rad_cnst_get_bin_props(list_idx, bin_idx, opticstype, &
   extpsw, abspsw, asmpsw, absplw, corefrac, nfrac, &
   wgtpct, nwtp, bcdust, nbcdust, kap, nkap, relh, nrelh, &
   sw_hygro_ext_wtp, sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, &
   sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, dryrad )

   ! Return requested properties for the bin from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,             intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,             intent(in)  :: bin_idx  ! mode index

   character(len=ot_length), optional, intent(out) :: opticstype

   real(r8),  optional, pointer     :: extpsw(:,:)
   real(r8),  optional, pointer     :: abspsw(:,:)
   real(r8),  optional, pointer     :: asmpsw(:,:)
   real(r8),  optional, pointer     :: absplw(:,:)
   real(r8),  optional, pointer     :: corefrac(:)
   integer,   optional, intent(out) :: nfrac

   real(r8),          optional, pointer     :: sw_hygro_ext_wtp(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa_wtp(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm_wtp(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext_wtp(:,:)
   real(r8),          optional, pointer     :: sw_hygro_coreshell_ext(:,:,:,:,:) 	! Pengfei Yu Mar.30
   real(r8),          optional, pointer     :: sw_hygro_coreshell_ssa(:,:,:,:,:)
   real(r8),          optional, pointer     :: sw_hygro_coreshell_asm(:,:,:,:,:)
   real(r8),          optional, pointer     :: lw_hygro_coreshell_ext(:,:,:,:,:)
   real(r8),  optional, pointer     :: wgtpct(:)
   real(r8),  optional, pointer     :: bcdust(:)
   real(r8),  optional, pointer     :: kap(:)
   real(r8),  optional, pointer     :: relh(:)
   integer,   optional, intent(out) :: nwtp
   integer,   optional, intent(out) :: nbcdust
   integer,   optional, intent(out) :: nkap
   integer,   optional, intent(out) :: nrelh
   real(r8),  optional, intent(out) :: dryrad

   ! Local variables
   integer :: id
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_props'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sa_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the physprop index for the requested bin
   id = slist%idx_props(bin_idx)

   if (present(opticstype))  call physprop_get(id, opticstype=opticstype)
   if (present(extpsw))      call physprop_get(id, extpsw2=extpsw)
   if (present(abspsw))      call physprop_get(id, abspsw2=abspsw)
   if (present(asmpsw))      call physprop_get(id, asmpsw2=asmpsw)
   if (present(absplw))      call physprop_get(id, absplw2=absplw)
   if (present(corefrac))    call physprop_get(id, corefrac=corefrac)
   if (present(nfrac))       call physprop_get(id, nfrac=nfrac)

   if (present(sw_hygro_ext_wtp))       call physprop_get(id, sw_hygro_ext_wtp=sw_hygro_ext_wtp)
   if (present(sw_hygro_ssa_wtp))       call physprop_get(id, sw_hygro_ssa_wtp=sw_hygro_ssa_wtp)
   if (present(sw_hygro_asm_wtp))       call physprop_get(id, sw_hygro_asm_wtp=sw_hygro_asm_wtp)
   if (present(lw_hygro_ext_wtp))       call physprop_get(id, lw_hygro_abs_wtp=lw_hygro_ext_wtp)
   if (present(sw_hygro_coreshell_ext)) call physprop_get(id, sw_hygro_coreshell_ext=sw_hygro_coreshell_ext)
   if (present(sw_hygro_coreshell_ssa)) call physprop_get(id, sw_hygro_coreshell_ssa=sw_hygro_coreshell_ssa)
   if (present(sw_hygro_coreshell_asm)) call physprop_get(id, sw_hygro_coreshell_asm=sw_hygro_coreshell_asm)
   if (present(lw_hygro_coreshell_ext)) call physprop_get(id, lw_hygro_coreshell_abs=lw_hygro_coreshell_ext)
   if (present(wgtpct))                 call physprop_get(id, wgtpct=wgtpct)
   if (present(bcdust))                 call physprop_get(id, bcdust=bcdust)
   if (present(kap))                    call physprop_get(id, kap=kap)
   if (present(relh))                   call physprop_get(id, relh=relh)
   if (present(nwtp))                   call physprop_get(id, nwtp=nwtp)
   if (present(nbcdust))                call physprop_get(id, nbcdust=nbcdust)
   if (present(nkap))                   call physprop_get(id, nkap=nkap)
   if (present(nrelh))                  call physprop_get(id, nrelh=nrelh)
   if (present(dryrad))                 call physprop_get(id, dryrad_aer=dryrad)

end subroutine rad_cnst_get_bin_props

!================================================================================================

end module aerosol_definition_mod
