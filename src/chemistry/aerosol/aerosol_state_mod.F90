module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use aerosol_properties_mod, only: AERO_FIELD_ADVECTED, AERO_FIELD_STORED
  use aerosol_properties_mod, only: AERO_FIELD_DERIVED, AERO_FIELD_ABSENT
  use aerosol_properties_mod, only: AERO_AMBIENT, AERO_CLDBRNE
  use aerosol_properties_mod, only: aero_has_working_state_table
  use physconst, only: pi
  use cam_abortutils, only: endrun

  implicit none

  private

  public :: aerosol_state
  public :: ptr2d_t

  !> aerosol_state defines the interface to the time-varying aerosol state
  !! variables (e.g., mixing ratios, number concentrations). This includes the
  !! aerosol portion of the overall model state.
  !!
  !! Each aerosol package (e.g., MAM, CARMA, etc) must extend the aerosol_state
  !! class to allow access to the state information (transported and not transported)
  !! of the aerosol package. Any package must implement each of the deferred
  !! procedures of the abstract aerosol_state class, may include additional private
  !! data members and type-bound procedures, and may override functions of the
  !! abstract class.
  !!
  !! Please see the modal_aerosol_state module for an example of how the aerosol_state
  !! class can be extended for a specific aerosol package.
  type, abstract :: aerosol_state
     ! aerosol properties of the model this state belongs to
     ! (non-owning pointer to the registry-owned instance, set at construction)
     class(aerosol_properties), pointer :: props_ => null()
     integer :: list_idx_ = 0 ! radiation climate/diagnostic list index
     integer :: ncol_ = 0     ! number of active columns
     integer :: nlev_ = 0     ! number of vertical layers
   contains
     procedure :: list_idx => get_list_idx
     procedure :: ncol => get_ncol
     procedure :: nlev => get_nlev
     procedure(aero_get_transported), deferred :: get_transported
     procedure(aero_set_transported), deferred :: set_transported
     procedure(aero_get_amb_total_bin_mmr), deferred :: ambient_total_bin_mmr
     ! deferred storage-alias hooks -- per-model pointer access to host-resident
     ! fields; call only through the guarded *_ptr accessors below (concrete
     ! classes may call their own hooks directly, e.g. the get_states loops)
     procedure(aero_get_state_mmr), deferred :: alias_ambient_mmr
     procedure(aero_get_state_mmr), deferred :: alias_cldbrne_mmr
     procedure(aero_get_state_num), deferred :: alias_ambient_num
     procedure(aero_get_state_num), deferred :: alias_cldbrne_num
     ! guarded storage-alias accessors -- zero-copy access, legal only for
     ! fields resident in host storage (advected constituents or pbuf fields)
     procedure :: ambient_mmr_ptr
     procedure :: cldbrne_mmr_ptr
     procedure :: ambient_num_ptr
     procedure :: cldbrne_num_ptr
     ! fill (value) getters -- universal: copy from host storage for
     ! ADVECTED/STORED fields, derive for DERIVED fields, zero for ABSENT
     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     ! per-model derivation hooks for DERIVED number fields
     procedure :: derive_ambient_num
     procedure :: derive_cldbrne_num
     ! working-state table over all bins/species -- capability-guarded
     ! (aero_has_working_state_table); replaces the deferred get_states
     procedure :: get_working_state
     procedure(aero_update_bin), deferred :: update_bin
     procedure :: loadaer
     procedure(aero_icenuc_size_wght_arr), deferred :: icenuc_size_wght_arr
     procedure(aero_icenuc_size_wght_val), deferred :: icenuc_size_wght_val
     generic :: icenuc_size_wght => icenuc_size_wght_arr,icenuc_size_wght_val
     procedure :: icenuc_type_wght_base
     procedure :: icenuc_type_wght => icenuc_type_wght_base
     procedure :: nuclice_get_numdens
     procedure :: get_amb_species_numdens
     procedure :: get_cld_species_numdens
     procedure :: coated_frac
     procedure :: mass_mean_radius
     procedure :: watact_mfactor
     procedure(aero_hetfrz_size_wght), deferred :: hetfrz_size_wght
     procedure(aero_hygroscopicity), deferred :: hygroscopicity
     procedure(aero_water_uptake), deferred :: water_uptake
     procedure(aero_wgtpct), deferred :: wgtpct
     procedure :: refractive_index_sw
     procedure :: refractive_index_lw
     procedure(aero_volume), deferred :: dry_volume
     procedure(aero_volume), deferred :: wet_volume
     procedure(aero_volume), deferred :: water_volume
     procedure(aero_wet_diam), deferred :: wet_diameter
     procedure :: convcld_actfrac
     procedure :: sol_factb_interstitial
     procedure(aero_aqu_gain_binfraction), deferred :: aqu_gain_binfraction

  end type aerosol_state

  ! for state fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  real(r8), parameter :: per_cm3 = 1.e-6_r8 ! factor for m-3 to cm-3 conversions
  real(r8), parameter :: per_m3 = 1.e6_r8   ! factor for cm-3 to m-3 conversions
  real(r8), parameter :: kg2mug = 1.e9_r8   ! factor for kg to micrograms (mug) conversions

  abstract interface

     !------------------------------------------------------------------------
     ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
     !------------------------------------------------------------------------
     function aero_get_amb_total_bin_mmr(self, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx      ! bin index
       integer, intent(in) :: col_ndx      ! column index
       integer, intent(in) :: lyr_ndx      ! vertical layer index

       real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

     end function aero_get_amb_total_bin_mmr

     !------------------------------------------------------------------------
     ! returns aerosol mass mixing ratio for a given species index and bin index
     !------------------------------------------------------------------------
     subroutine aero_get_state_mmr(self, species_ndx, bin_ndx, mmr)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: species_ndx  ! species index
       integer, intent(in) :: bin_ndx      ! bin index
       real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)
     end subroutine aero_get_state_mmr

     !------------------------------------------------------------------------
     ! returns aerosol number mixing ratio for a given species index and bin index
     !------------------------------------------------------------------------
     subroutine aero_get_state_num(self, bin_ndx, num)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx     ! bin index
       real(r8), pointer   :: num(:,:)    ! number densities (ncol,nlev)
     end subroutine aero_get_state_num

     !------------------------------------------------------------------------------
     ! sets transported components
     ! This updates the aerosol model state from the host transported aerosol constituents array.
     ! (mass mixing ratios or number mixing ratios)
     !------------------------------------------------------------------------------
     subroutine aero_set_transported( self, transported_array )
       import :: aerosol_state, r8
       class(aerosol_state), intent(inout) :: self
       real(r8), intent(in) :: transported_array(:,:,:)
     end subroutine aero_set_transported

     !------------------------------------------------------------------------------
     ! returns transported components
     ! This updates the transported aerosol constituent array to match the aerosol model state.
     ! (mass mixing ratios or number mixing ratios)
     !------------------------------------------------------------------------------
     subroutine aero_get_transported( self, transported_array )
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       real(r8), intent(out) :: transported_array(:,:,:)
     end subroutine aero_get_transported

     !------------------------------------------------------------------------------
     ! return aerosol bin size weights for a given bin
     !------------------------------------------------------------------------------
     subroutine aero_icenuc_size_wght_arr(self, bin_ndx, species_type, use_preexisting_ice, wght)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx             ! bin number
       character(len=*), intent(in) :: species_type  ! species type
       logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
       real(r8), intent(out) :: wght(:,:)

     end subroutine aero_icenuc_size_wght_arr

     !------------------------------------------------------------------------------
     ! return aerosol bin size weights for a given bin, column and vertical layer
     !------------------------------------------------------------------------------
     subroutine aero_icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx                ! bin number
       integer, intent(in) :: col_ndx                ! column index
       integer, intent(in) :: lyr_ndx                ! vertical layer index
       character(len=*), intent(in) :: species_type  ! species type
       logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
       real(r8), intent(out) :: wght

     end subroutine aero_icenuc_size_wght_val

     !------------------------------------------------------------------------------
     ! updates state and tendency
     !------------------------------------------------------------------------------
     subroutine aero_update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx                ! bin number
       integer, intent(in) :: col_ndx                ! column index
       integer, intent(in) :: lyr_ndx                ! vertical layer index
       real(r8),intent(in) :: delmmr_sum             ! mass mixing ratio change summed over all species in bin
       real(r8),intent(in) :: delnum_sum             ! number mixing ratio change summed over all species in bin
       integer, intent(in) :: tnd_ndx                ! tendency index
       real(r8),intent(in) :: dtime                  ! time step size (sec)
       real(r8),intent(inout) :: tend(:,:,:)         ! tendency

     end subroutine aero_update_bin

     !------------------------------------------------------------------------------
     ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
     ! as heterogeneous freezing nuclei
     !------------------------------------------------------------------------------
     function aero_hetfrz_size_wght(self, bin_ndx) result(wght)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx             ! bin number

       real(r8) :: wght(self%ncol_,self%nlev_)

     end function aero_hetfrz_size_wght

     !------------------------------------------------------------------------------
     ! returns hygroscopicity for a given radiation diagnostic list number and
     ! bin number
     !------------------------------------------------------------------------------
     subroutine aero_hygroscopicity(self, bin_ndx, kappa)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx      ! bin number

       real(r8), intent(out) :: kappa(:,:)              ! hygroscopicity (ncol,nlev)

     end subroutine aero_hygroscopicity

     !------------------------------------------------------------------------------
     ! returns aerosol wet diameter and aerosol water concentration for a given
     ! radiation diagnostic list number and bin number
     !------------------------------------------------------------------------------
     subroutine aero_water_uptake(self, bin_idx, dgnumwet, qaerwat)
       import :: aerosol_state, r8

       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_idx              ! bin number
       real(r8),intent(out) :: dgnumwet(self%ncol_,self%nlev_) ! aerosol wet diameter (m)
       real(r8),intent(out) :: qaerwat(self%ncol_,self%nlev_)  ! aerosol water concentration (g/g)

     end subroutine aero_water_uptake

     !------------------------------------------------------------------------------
     ! aerosol weight percent of H2SO4/H2O solution
     !------------------------------------------------------------------------------
     function aero_wgtpct(self) result(wtp)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       real(r8) :: wtp(self%ncol_,self%nlev_)  ! weight percent of H2SO4/H2O solution for given icol, ilev

     end function aero_wgtpct

     !------------------------------------------------------------------------------
     ! aerosol volume interface
     !------------------------------------------------------------------------------
     function aero_volume(self, bin_idx) result(vol)
       import :: aerosol_state, r8

       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_idx   ! bin number

       real(r8) :: vol(self%ncol_,self%nlev_)  ! m3/kg

     end function aero_volume

     !------------------------------------------------------------------------------
     ! aerosol wet diameter
     !------------------------------------------------------------------------------
     function aero_wet_diam(self, bin_idx) result(diam)
       import :: aerosol_state,  r8

       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_idx   ! bin number

       real(r8) :: diam(self%ncol_,self%nlev_)

     end function aero_wet_diam

     !------------------------------------------------------------------------------
     ! aqueous chemistry partitioning -- used in sox_cldaero_update
     !------------------------------------------------------------------------------
     subroutine aero_aqu_gain_binfraction(self, type, qcw, delso4_o3rxn, faqgain)
       import :: aerosol_state, r8

       class(aerosol_state), intent(in) :: self
       character(len=*), intent(in) :: type                ! aerosol species type
       real(r8), intent(in) :: qcw(:,:,:)                  ! cloud-borne aerosol volume mixing ratio
       real(r8), intent(in) :: delso4_o3rxn(:,:)           ! sulfate concentration change due to oxidation
       real(r8), intent(out) :: faqgain(:,:,:)             ! fraction gain in each mode / bin

     end subroutine aero_aqu_gain_binfraction

  end interface

contains

  !------------------------------------------------------------------------------
  ! returns the radiation climate/diagnostic list index
  !------------------------------------------------------------------------------
  pure integer function get_list_idx(self)
    class(aerosol_state), intent(in) :: self
    get_list_idx = self%list_idx_
  end function get_list_idx

  !------------------------------------------------------------------------------
  ! returns the number of active columns
  !------------------------------------------------------------------------------
  pure integer function get_ncol(self)
    class(aerosol_state), intent(in) :: self
    get_ncol = self%ncol_
  end function get_ncol

  !------------------------------------------------------------------------------
  ! returns the number of vertical layers
  !------------------------------------------------------------------------------
  pure integer function get_nlev(self)
    class(aerosol_state), intent(in) :: self
    get_nlev = self%nlev_
  end function get_nlev

  !------------------------------------------------------------------------------
  ! returns a pointer to the ambient aerosol mass mixing ratio for a given
  ! species index and bin index -- storage-alias access, legal only for fields
  ! resident in host storage (advected constituents or pbuf fields)
  !------------------------------------------------------------------------------
  subroutine ambient_mmr_ptr(self, species_ndx, bin_ndx, mmr)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    select case (self%props_%field_kind(bin_ndx, species_ndx, AERO_AMBIENT))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_ambient_mmr(species_ndx, bin_ndx, mmr)
    case default
       call endrun('aerosol_state ambient_mmr_ptr: field is not host-resident'// &
                   ' -- use the get_ambient_mmr fill getter')
    end select
  end subroutine ambient_mmr_ptr

  !------------------------------------------------------------------------------
  ! returns a pointer to the cloud-borne aerosol mass mixing ratio for a given
  ! species index and bin index -- storage-alias access, legal only for fields
  ! resident in host storage (advected constituents or pbuf fields)
  !------------------------------------------------------------------------------
  subroutine cldbrne_mmr_ptr(self, species_ndx, bin_ndx, mmr)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    select case (self%props_%field_kind(bin_ndx, species_ndx, AERO_CLDBRNE))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_cldbrne_mmr(species_ndx, bin_ndx, mmr)
    case default
       call endrun('aerosol_state cldbrne_mmr_ptr: field is not host-resident'// &
                   ' -- use the get_cldbrne_mmr fill getter')
    end select
  end subroutine cldbrne_mmr_ptr

  !------------------------------------------------------------------------------
  ! returns a pointer to the ambient aerosol number mixing ratio for a given
  ! bin index -- storage-alias access, legal only for fields resident in host
  ! storage (advected constituents or pbuf fields)
  !------------------------------------------------------------------------------
  subroutine ambient_num_ptr(self, bin_ndx, num)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number densities (ncol,nlev)

    select case (self%props_%field_kind(bin_ndx, 0, AERO_AMBIENT))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_ambient_num(bin_ndx, num)
    case default
       call endrun('aerosol_state ambient_num_ptr: field is not host-resident'// &
                   ' -- use the get_ambient_num fill getter')
    end select
  end subroutine ambient_num_ptr

  !------------------------------------------------------------------------------
  ! returns a pointer to the cloud-borne aerosol number mixing ratio for a given
  ! bin index -- storage-alias access, legal only for fields resident in host
  ! storage (advected constituents or pbuf fields)
  !------------------------------------------------------------------------------
  subroutine cldbrne_num_ptr(self, bin_ndx, num)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number densities (ncol,nlev)

    select case (self%props_%field_kind(bin_ndx, 0, AERO_CLDBRNE))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_cldbrne_num(bin_ndx, num)
    case default
       call endrun('aerosol_state cldbrne_num_ptr: field is not host-resident'// &
                   ' -- use the get_cldbrne_num fill getter')
    end select
  end subroutine cldbrne_num_ptr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio values for a given species index
  ! and bin index -- universal fill (value) getter; the caller owns the output
  ! array (first extent >= ncol; columns 1:ncol are written, any remainder is
  ! left untouched)
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self, species_ndx, bin_ndx, mmr)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), intent(out) :: mmr(:,:)   ! mass mixing ratios (ncol,nlev)

    real(r8), pointer :: ptr(:,:)
    integer :: nc

    nc = self%ncol_

    select case (self%props_%field_kind(bin_ndx, species_ndx, AERO_AMBIENT))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_ambient_mmr(species_ndx, bin_ndx, ptr)
       mmr(:nc,:) = ptr(:nc,:)
    case (AERO_FIELD_ABSENT)
       mmr(:nc,:) = 0._r8
    case default
       call endrun('aerosol_state get_ambient_mmr: no derivation defined for mass fields')
    end select
  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol mass mixing ratio values for a given species
  ! index and bin index -- universal fill (value) getter; the caller owns the
  ! output array (first extent >= ncol; columns 1:ncol are written, any
  ! remainder is left untouched)
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), intent(out) :: mmr(:,:)   ! mass mixing ratios (ncol,nlev)

    real(r8), pointer :: ptr(:,:)
    integer :: nc

    nc = self%ncol_

    select case (self%props_%field_kind(bin_ndx, species_ndx, AERO_CLDBRNE))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_cldbrne_mmr(species_ndx, bin_ndx, ptr)
       mmr(:nc,:) = ptr(:nc,:)
    case (AERO_FIELD_ABSENT)
       mmr(:nc,:) = 0._r8
    case default
       call endrun('aerosol_state get_cldbrne_mmr: no derivation defined for mass fields')
    end select
  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number mixing ratio values for a given bin index --
  ! universal fill (value) getter; the caller owns the output array (first
  ! extent >= ncol; columns 1:ncol are written, any remainder is left untouched)
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), intent(out) :: num(:,:)  ! number mixing ratios (ncol,nlev)

    real(r8), pointer :: ptr(:,:)
    integer :: nc

    nc = self%ncol_

    select case (self%props_%field_kind(bin_ndx, 0, AERO_AMBIENT))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_ambient_num(bin_ndx, ptr)
       num(:nc,:) = ptr(:nc,:)
    case (AERO_FIELD_DERIVED)
       call self%derive_ambient_num(bin_ndx, num)
    case (AERO_FIELD_ABSENT)
       num(:nc,:) = 0._r8
    case default
       call endrun('aerosol_state get_ambient_num: unclassified field kind')
    end select
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio values for a given bin index
  ! -- universal fill (value) getter; the caller owns the output array (first
  ! extent >= ncol; columns 1:ncol are written, any remainder is left untouched)
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), intent(out) :: num(:,:)  ! number mixing ratios (ncol,nlev)

    real(r8), pointer :: ptr(:,:)
    integer :: nc

    nc = self%ncol_

    select case (self%props_%field_kind(bin_ndx, 0, AERO_CLDBRNE))
    case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
       call self%alias_cldbrne_num(bin_ndx, ptr)
       num(:nc,:) = ptr(:nc,:)
    case (AERO_FIELD_DERIVED)
       call self%derive_cldbrne_num(bin_ndx, num)
    case (AERO_FIELD_ABSENT)
       num(:nc,:) = 0._r8
    case default
       call endrun('aerosol_state get_cldbrne_num: unclassified field kind')
    end select
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! derives the ambient aerosol number mixing ratio for a given bin index into
  ! the caller's array -- overridden by models with DERIVED ambient number
  ! fields; reached only through the get_ambient_num fill getter
  !------------------------------------------------------------------------------
  subroutine derive_ambient_num(self, bin_ndx, num)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), intent(out) :: num(:,:)  ! number mixing ratios (ncol,nlev)

    call endrun('aerosol_state derive_ambient_num: no derivation defined for this model')
  end subroutine derive_ambient_num

  !------------------------------------------------------------------------------
  ! derives the cloud-borne aerosol number mixing ratio for a given bin index
  ! into the caller's array -- overridden by models with DERIVED cloud-borne
  ! number fields; reached only through the get_cldbrne_num fill getter
  !------------------------------------------------------------------------------
  subroutine derive_cldbrne_num(self, bin_ndx, num)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), intent(out) :: num(:,:)  ! number mixing ratios (ncol,nlev)

    call endrun('aerosol_state derive_cldbrne_num: no derivation defined for this model')
  end subroutine derive_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns the working-state table of interstitial and cloud-borne aerosol
  ! fields over all bins and species, indexed by the species indexer.
  !
  ! For ADVECTED or STORED (pbuf/non-advected) fields:
  !   - returns pointer to host storage
  !   - for writes: ADVECTED constituents are updated via ptend and never via pointer.
  !                 STORED (pbuf/non-advected) are updated in-place via pointer.
  ! For DERIVED fields:
  !   - caller must allocate a slab target passed in here,
  !     so derived fields are calculated and written into this caller-managed memory.
  !   - the purpose of the caller-supplied scratch space is because the aerosol state
  !     itself should not manage the lifecycle of memory for data that is derived
  !     and discarded after the caller parameterization uses it.
  !------------------------------------------------------------------------------
  subroutine get_working_state(self, raer, qqcw, scratch)
    class(aerosol_state), intent(in) :: self
    type(ptr2d_t), intent(out) :: raer(:) ! working state of interstitial aerosols
    type(ptr2d_t), intent(out) :: qqcw(:) ! working state of cloud-borne aerosols
    real(r8), target, intent(inout) :: scratch(:,:,:) ! caller-owned backing for DERIVED
                                                      ! entries; third extent >=
                                                      ! props%num_derived_working_entries()

    integer :: ibin, ispc, indx, islab

    if (.not. self%props_%supports(aero_has_working_state_table)) then
       call endrun('aerosol_state get_working_state: working-state table not supported'// &
                   ' for this aerosol model -- use the fill getters per field_kind')
    end if

    islab = 0
    do ibin = 1, self%props_%nbins()
       indx = self%props_%indexer(ibin, 0)
       select case (self%props_%field_kind(ibin, 0, AERO_AMBIENT))
       case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
          call self%alias_ambient_num(ibin, raer(indx)%fld)
       case (AERO_FIELD_DERIVED)
          islab = islab + 1
          call self%derive_ambient_num(ibin, scratch(:,:,islab))
          raer(indx)%fld => scratch(:,:,islab)
       case default
          call endrun('aerosol_state get_working_state: unsupported ambient number entry')
       end select
       select case (self%props_%field_kind(ibin, 0, AERO_CLDBRNE))
       case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
          call self%alias_cldbrne_num(ibin, qqcw(indx)%fld)
       case (AERO_FIELD_DERIVED)
          islab = islab + 1
          call self%derive_cldbrne_num(ibin, scratch(:,:,islab))
          qqcw(indx)%fld => scratch(:,:,islab)
       case default
          call endrun('aerosol_state get_working_state: unsupported cloud-borne number entry')
       end select
       do ispc = 1, self%props_%nspecies(ibin)
          indx = self%props_%indexer(ibin, ispc)
          select case (self%props_%field_kind(ibin, ispc, AERO_AMBIENT))
          case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
             call self%alias_ambient_mmr(species_ndx=ispc, bin_ndx=ibin, mmr=raer(indx)%fld)
          case default
             call endrun('aerosol_state get_working_state: unsupported ambient mass entry')
          end select
          select case (self%props_%field_kind(ibin, ispc, AERO_CLDBRNE))
          case (AERO_FIELD_ADVECTED, AERO_FIELD_STORED)
             call self%alias_cldbrne_mmr(species_ndx=ispc, bin_ndx=ibin, mmr=qqcw(indx)%fld)
          case default
             call endrun('aerosol_state get_working_state: unsupported cloud-borne mass entry')
          end select
       end do
    end do

  end subroutine get_working_state

  !------------------------------------------------------------------------------
  ! returns aerosol number, volume concentrations, and bulk hygroscopicity
  !------------------------------------------------------------------------------
  subroutine loadaer( self, m, cs, phase, &
                       naerosol, vaerosol, hygro, errnum, errstr, pom_hygro)

    ! input arguments
    class(aerosol_state), intent(in) :: self
    integer,  intent(in) :: m           ! mode or bin index
    real(r8), intent(in) :: cs(:,:)     ! air density (kg/m3)
    integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

    ! output arguments
    real(r8), intent(out) :: naerosol(:,:)  ! number conc (1/m3)
    real(r8), intent(out) :: vaerosol(:,:)  ! volume conc (m3/m3)
    real(r8), intent(out) :: hygro(:,:)     ! bulk hygroscopicity of mode

    integer ,         intent(out) :: errnum
    character(len=*), intent(out) :: errstr

    real(r8), optional, intent(in) :: pom_hygro     ! POM hygroscopicity override

    ! internal
    real(r8), pointer :: raer(:,:) ! interstitial aerosol mass mixing ratios
    real(r8), pointer :: qqcw(:,:) ! cloud-borne aerosol mass mixing ratios
    real(r8) :: specdens, spechygro
    character(len=aero_name_len) :: spectype

    real(r8) :: vol(self%ncol_,self%nlev_)  ! aerosol volume mixing ratio
    real(r8) :: anum(self%ncol_,self%nlev_) ! interstitial aerosol number mixing ratio
    real(r8) :: cnum(self%ncol_,self%nlev_) ! cloud-borne aerosol number mixing ratio
    integer  :: l, ncol, nlev
    !-------------------------------------------------------------------------------
    errnum = 0

    ncol = self%ncol()
    nlev = self%nlev()

    vaerosol(:,:) = 0._r8
    hygro(:,:)    = 0._r8

    do l = 1, self%props_%nspecies(m)

       call self%ambient_mmr_ptr(species_ndx=l, bin_ndx=m, mmr=raer)
       call self%cldbrne_mmr_ptr(species_ndx=l, bin_ndx=m, mmr=qqcw)
       call self%props_%get(m,l, density=specdens, hygro=spechygro, spectype=spectype)
       if (present(pom_hygro)) then
          if (spectype=='p-organic'.and.pom_hygro>0._r8) then
             spechygro=pom_hygro
          endif
       endif

       if (phase == 3) then
          vol(:ncol,:) = max(raer(:ncol,:) + qqcw(:ncol,:), 0._r8)/specdens
       else if (phase == 2) then
          vol(:ncol,:) = max(qqcw(:ncol,:), 0._r8)/specdens
       else if (phase == 1) then
          vol(:ncol,:) = max(raer(:ncol,:), 0._r8)/specdens
       else
          errnum = -1
          write(errstr,*)'phase = ',phase,' in aerosol_state::loadaer not recognized'
          return
       end if

       vaerosol(:ncol,:) = vaerosol(:ncol,:) + vol(:ncol,:)
       hygro(:ncol,:)    = hygro(:ncol,:) + vol(:ncol,:)*spechygro
    end do

    where(vaerosol(:ncol,:) > 1.0e-30_r8)
       hygro(:ncol,:)    = hygro(:ncol,:)/(vaerosol(:ncol,:))
       vaerosol(:ncol,:) = vaerosol(:ncol,:)*cs(:ncol,:)
    elsewhere
       hygro(:ncol,:)    = 0._r8
       vaerosol(:ncol,:) = 0._r8
    end where

    ! aerosol number mixing ratios (#/kg)
    call self%get_ambient_num(m, anum)
    call self%get_cldbrne_num(m, cnum)
    if (phase == 3) then
       naerosol(:ncol,:) = (anum(:ncol,:) + cnum(:ncol,:))*cs(:ncol,:) ! #/kg -> #/m3
    else if (phase == 2) then
       naerosol(:ncol,:) = cnum(:ncol,:)*cs(:ncol,:)
    else
       naerosol(:ncol,:) = anum(:ncol,:)*cs(:ncol,:)
    end if

    ! adjust number
    call self%props_%apply_number_limits( naerosol, vaerosol, ncol, nlev, m )

  end subroutine loadaer

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number density for a given bin number and species type
  !------------------------------------------------------------------------------
  subroutine get_amb_species_numdens(self, bin_ndx, species_type, rho, numdens)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: numdens(:,:)         ! species number densities (#/cm^3)

    real(r8) :: num(self%ncol_,self%nlev_)
    real(r8) :: type_wght(self%ncol_,self%nlev_)
    real(r8) :: size_wght(self%ncol_,self%nlev_)
    integer :: ncol

    ncol = self%ncol()

    size_wght = self%hetfrz_size_wght(bin_ndx)

    call self%icenuc_type_wght_base(bin_ndx, species_type, rho, type_wght)

    call self%get_ambient_num(bin_ndx, num)

    numdens(:ncol,:) = num(:ncol,:)*rho(:ncol,:)*type_wght(:ncol,:)*size_wght(:ncol,:)*per_cm3

  end subroutine get_amb_species_numdens

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number density for a given bin number and species type
  !------------------------------------------------------------------------------
  subroutine get_cld_species_numdens(self, bin_ndx, species_type, rho, numdens)
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: numdens(:,:)         ! number densities (#/cm^3)

    real(r8) :: num(self%ncol_,self%nlev_)
    real(r8) :: type_wght(self%ncol_,self%nlev_)
    real(r8) :: size_wght(self%ncol_,self%nlev_)
    integer :: ncol

    ncol = self%ncol()

    size_wght = self%hetfrz_size_wght(bin_ndx)

    call self%icenuc_type_wght_base(bin_ndx, species_type, rho, type_wght, cloud_borne=.true.)

    call self%get_cldbrne_num(bin_ndx, num)

    numdens(:ncol,:) = num(:ncol,:)*rho(:ncol,:)*type_wght(:ncol,:)*size_wght(:ncol,:)*per_cm3

  end subroutine get_cld_species_numdens

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght_base(self, bin_ndx, species_type, rho, wght, cloud_borne)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights
    logical, optional, intent(in) :: cloud_borne  ! if TRUE cloud-borne aerosols are used
                                                  ! otherwise ambient aerosols are used

    real(r8) :: mass(self%ncol_,self%nlev_)
    real(r8) :: totalmass(self%ncol_,self%nlev_)
    real(r8), pointer :: aer_bin(:,:)

    character(len=aero_name_len) :: spectype, sptype
    integer :: ispc, ncol
    logical :: cldbrne

    if (present(cloud_borne)) then
       cldbrne = cloud_borne
    else
       cldbrne = .false.
    end if

    ncol = self%ncol()

    wght(:,:) = 0._r8
    totalmass(:,:) = 0._r8
    mass(:,:)   = 0._r8

    if (species_type=='sulfate_strat') then
       sptype = 'sulfate'
    else
       sptype = species_type
    end if

    do ispc = 1, self%props_%nspecies(bin_ndx)

       if (cldbrne) then
          call self%cldbrne_mmr_ptr(species_ndx=ispc, bin_ndx=bin_ndx, mmr=aer_bin)
       else
          call self%ambient_mmr_ptr(species_ndx=ispc, bin_ndx=bin_ndx, mmr=aer_bin)
       end if
       call self%props_%species_type(bin_ndx, ispc, spectype=spectype)

       totalmass(:ncol,:) = totalmass(:ncol,:) + aer_bin(:ncol,:)*rho(:ncol,:)

       if (trim(spectype) == trim(sptype)) then
          mass(:ncol,:) = mass(:ncol,:) + aer_bin(:ncol,:)*rho(:ncol,:)
       end if

    end do

    where (totalmass(:ncol,:) > 0._r8)
       wght(:ncol,:) = mass(:ncol,:)/totalmass(:ncol,:)
    end where

  end subroutine icenuc_type_wght_base

  !------------------------------------------------------------------------------
  subroutine nuclice_get_numdens(self, use_preexisting_ice, rho, dust_num_col, sulf_num_col, soot_num_col, sulf_num_tot_col )

    class(aerosol_state), intent(in) :: self

    logical, intent(in) :: use_preexisting_ice
    real(r8), intent(in) :: rho(:,:) ! air density (kg m-3)
    real(r8), intent(out) :: dust_num_col(:,:) ! dust number densities (#/cm^3)
    real(r8), intent(out) :: sulf_num_col(:,:) ! sulfate number densities (#/cm^3)
    real(r8), intent(out) :: soot_num_col(:,:) ! soot number densities (#/cm^3)
    real(r8), intent(out) :: sulf_num_tot_col(:,:) ! stratopsheric sulfate number densities (#/cm^3)

    integer :: ibin,ispc,ncol
    character(len=aero_name_len) :: spectype
    real(r8) :: size_wghts(self%ncol_,self%nlev_)
    real(r8) :: type_wghts(self%ncol_,self%nlev_)

    real(r8) :: num_col(self%ncol_,self%nlev_)

    ncol = self%ncol()

    dust_num_col(:,:) = 0._r8
    sulf_num_col(:,:) = 0._r8
    soot_num_col(:,:) = 0._r8
    sulf_num_tot_col(:,:) = 0._r8

    ! collect number densities (#/cm^3) for dust, sulfate, and soot
    do ibin = 1,self%props_%nbins()

       call self%get_ambient_num(ibin, num_col)

       do ispc = 1,self%props_%nspecies(ibin)

          call self%props_%species_type(ibin, ispc, spectype)

          call self%icenuc_size_wght(ibin, spectype, use_preexisting_ice, size_wghts)

          call self%icenuc_type_wght(ibin, spectype, rho, type_wghts)

          select case ( trim(spectype) )
          case('dust')
             dust_num_col(:ncol,:) = dust_num_col(:ncol,:) &
                  + size_wghts(:ncol,:)*type_wghts(:ncol,:)*num_col(:ncol,:)*rho(:ncol,:)*per_cm3
          case('sulfate')
             ! This order of ops gives bit-for-bit results for cam5 phys ( use_preexisting_ice = .false. )
             sulf_num_col(:ncol,:) = sulf_num_col(:ncol,:) &
                  + num_col(:ncol,:)*rho(:ncol,:)*per_cm3  * size_wghts(:ncol,:)*type_wghts(:ncol,:)
          case('black-c')
             soot_num_col(:ncol,:) = soot_num_col(:ncol,:) &
                  + size_wghts(:ncol,:)*type_wghts(:ncol,:)*num_col(:ncol,:)*rho(:ncol,:)*per_cm3
          end select

       enddo

       ! stratospheric sulfates -- special case not included in the species loop above
       call self%icenuc_size_wght(ibin, 'sulfate_strat', use_preexisting_ice, size_wghts)
       call self%icenuc_type_wght(ibin, 'sulfate_strat', rho, type_wghts)
       sulf_num_tot_col(:ncol,:) = sulf_num_tot_col(:ncol,:) &
            + size_wghts(:ncol,:)*type_wghts(:ncol,:)*num_col(:ncol,:)*rho(:ncol,:)*per_cm3

    enddo

  end subroutine nuclice_get_numdens

  !------------------------------------------------------------------------------
  ! returns the fraction of particle surface area of aerosol subset `bin_ndx` covered
  ! by at least a monolayer of species `species_type` [0-1]
  !------------------------------------------------------------------------------
  function coated_frac(self, bin_ndx,  species_type, radius) result(frac)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    real(r8), intent(in) :: radius(:,:) ! m

    real(r8) :: frac(self%ncol_,self%nlev_)                  ! coated fraction

    !------------coated variables--------------------
    real(r8), parameter :: n_so4_monolayers_dust = 1.0_r8 ! number of so4(+nh4) monolayers needed to coat a dust particle
    real(r8), parameter :: dr_so4_monolayers_dust = n_so4_monolayers_dust * 4.76e-10_r8
    real(r8) :: vol_shell(self%ncol_,self%nlev_)
    real(r8) :: vol_core(self%ncol_,self%nlev_)
    real(r8) :: alnsg, fac_volsfc
    real(r8) :: tmp1(self%ncol_,self%nlev_), tmp2(self%ncol_,self%nlev_)
    real(r8),pointer :: sulf_mmr(:,:)
    real(r8),pointer :: soa_mmr(:,:)
    real(r8),pointer :: pom_mmr(:,:)
    real(r8),pointer :: aer_mmr(:,:)

    integer :: sulf_ndx
    integer :: soa_ndx
    integer :: pom_ndx
    integer :: species_ndx

    real(r8) :: specdens_so4
    real(r8) :: specdens_pom
    real(r8) :: specdens_soa
    real(r8) :: specdens

    character(len=aero_name_len) :: spectype
    integer :: ispc, ncol

    ncol = self%ncol()

    frac = -huge(1._r8)

    sulf_ndx = -1
    pom_ndx = -1
    soa_ndx = -1
    species_ndx = -1

    do ispc = 1, self%props_%nspecies(bin_ndx)
       call self%props_%species_type(bin_ndx, ispc, spectype)

       select case ( trim(spectype) )
       case('sulfate')
          sulf_ndx = ispc
       case('p-organic')
          pom_ndx = ispc
       case('s-organic')
          soa_ndx = ispc
       end select
       if (spectype==species_type) then
          species_ndx = ispc
       end if
    end do

    vol_shell(:ncol,:) = 0._r8

    if (sulf_ndx>0) then
       call self%props_%get(bin_ndx, sulf_ndx, density=specdens_so4)
       call self%ambient_mmr_ptr(species_ndx=sulf_ndx, bin_ndx=bin_ndx, mmr=sulf_mmr)
       vol_shell(:ncol,:) = vol_shell(:ncol,:) + sulf_mmr(:ncol,:)/specdens_so4
    end if
    if (pom_ndx>0) then
       call self%props_%get(bin_ndx, pom_ndx, density=specdens_pom)
       call self%ambient_mmr_ptr(species_ndx=pom_ndx, bin_ndx=bin_ndx, mmr=pom_mmr)
       vol_shell(:ncol,:) = vol_shell(:ncol,:) + pom_mmr(:ncol,:)*self%props_%pom_equivso4_factor()/specdens_pom
    end if
    if (soa_ndx>0) then
       call self%props_%get(bin_ndx, soa_ndx, density=specdens_soa)
       call self%ambient_mmr_ptr(species_ndx=soa_ndx, bin_ndx=bin_ndx, mmr=soa_mmr)
       vol_shell(:ncol,:) = vol_shell(:ncol,:) + soa_mmr(:ncol,:)*self%props_%soa_equivso4_factor()/specdens_soa
    end if

    call self%props_%get(bin_ndx, species_ndx, density=specdens)
    call self%ambient_mmr_ptr(species_ndx=species_ndx, bin_ndx=bin_ndx, mmr=aer_mmr)
    vol_core(:ncol,:) = aer_mmr(:ncol,:)/specdens

    alnsg = self%props_%alogsig(bin_ndx)
    fac_volsfc = exp(2.5_r8*alnsg**2)

    tmp1(:ncol,:) = vol_shell(:ncol,:)*(radius(:ncol,:)*2._r8)*fac_volsfc
    tmp2(:ncol,:) = max(6.0_r8*dr_so4_monolayers_dust*vol_core(:ncol,:), 0.0_r8)

    where(tmp1(:ncol,:)>0._r8 .and. tmp2(:ncol,:)>0._r8)
       frac(:ncol,:) = tmp1(:ncol,:)/tmp2(:ncol,:)
    elsewhere
       frac(:ncol,:) = 0.001_r8
    end where

    where(frac(:ncol,:)>1._r8)
       frac(:ncol,:) = 1._r8
    end where
    where(frac(:ncol,:) < 0.001_r8)
       frac(:ncol,:) = 0.001_r8
    end where

  end function coated_frac

  !------------------------------------------------------------------------------
  ! returns the radius [m] of particles in aerosol subset `bin_ndx` assuming all particles are
  ! the same size and only species `species_ndx` contributes to the particle volume
  !------------------------------------------------------------------------------
  function mass_mean_radius(self, bin_ndx, species_ndx, rho) result(radius)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: species_ndx            ! species number
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)

    real(r8) :: radius(self%ncol_,self%nlev_) ! m

    character(len=aero_name_len) :: species_type
    real(r8) :: aer_numdens(self%ncol_,self%nlev_) ! kg/m3
    real(r8) :: aer_massdens(self%ncol_,self%nlev_) ! kg/m3
    real(r8),pointer :: aer_mmr(:,:) ! kg/kg

    real(r8) :: specdens,minrad
    real(r8) :: wght(self%ncol_,self%nlev_)
    integer :: i,k,ncol,nlev

    ncol = self%ncol()
    nlev = self%nlev()

    wght = self%hetfrz_size_wght(bin_ndx)

    call self%props_%species_type(bin_ndx, species_ndx, spectype=species_type)

    call self%props_%get(bin_ndx, species_ndx, density=specdens) ! kg/m3
    call self%ambient_mmr_ptr(species_ndx=species_ndx, bin_ndx=bin_ndx, mmr=aer_mmr) ! kg/kg
    call self%get_amb_species_numdens(bin_ndx, species_type, rho, aer_numdens) ! #/cm3

    aer_massdens(:ncol,:) = aer_mmr(:ncol,:)*rho(:ncol,:)*wght(:ncol,:) ! kg/m3

    minrad = self%props_%min_mass_mean_rad(bin_ndx, species_ndx)

    do k = 1,nlev
       do i = 1,ncol
          if (aer_massdens(i,k)*1.0e-3_r8 > 1.0e-30_r8 .and. aer_numdens(i,k) > 1.0e-3_r8) then
             radius(i,k) = (3._r8/(4*pi*specdens)*aer_massdens(i,k)/(aer_numdens(i,k)*per_m3))**(1._r8/3._r8) ! m
          else
             radius(i,k) = minrad
          end if
       end do
    end do

  end function mass_mean_radius

  !------------------------------------------------------------------------------
  ! calculates water activity mass factor -- density*(1.-(OC+BC)/(OC+BC+SO4)) [mug m-3]
  ! of species `species_type` in subset `bin_ndx`
  !------------------------------------------------------------------------------
  subroutine watact_mfactor(self, bin_ndx,  species_type, rho, wact_factor)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wact_factor(:,:)     ! water activity factor -- density*(1.-(OC+BC)/(OC+BC+SO4)) [mug m-3]

    real(r8), pointer :: aer_mmr(:,:)
    real(r8) :: bin_num(self%ncol_,self%nlev_)
    real(r8) :: tot2_mmr(self%ncol_,self%nlev_)
    real(r8) :: tot1_mmr(self%ncol_,self%nlev_)
    real(r8) :: aer_numdens(self%ncol_,self%nlev_)
    integer :: ispc, ncol
    character(len=aero_name_len) :: spectype

    real(r8) :: awcam(self%ncol_,self%nlev_)          ! mass density [mug m-3]
    real(r8) :: awfacm(self%ncol_,self%nlev_)         ! mass factor ! (OC+BC)/(OC+BC+SO4)

    ncol = self%ncol()

    tot2_mmr = 0.0_r8
    tot1_mmr = 0.0_r8

    if (self%props_%soluble(bin_ndx)) then

       do ispc = 1, self%props_%nspecies(bin_ndx)

          call self%props_%species_type(bin_ndx, ispc, spectype)

          if (trim(spectype)=='black-c' .or. trim(spectype)=='p-organic' .or. trim(spectype)=='s-organic') then
             call self%ambient_mmr_ptr(species_ndx=ispc, bin_ndx=bin_ndx, mmr=aer_mmr)
             tot2_mmr(:ncol,:) = tot2_mmr(:ncol,:) + aer_mmr(:ncol,:)
          end if
          if (trim(spectype)=='sulfate') then
             call self%ambient_mmr_ptr(species_ndx=ispc, bin_ndx=bin_ndx, mmr=aer_mmr)
             tot1_mmr(:ncol,:) = tot1_mmr(:ncol,:) + aer_mmr(:ncol,:)
          end if
       end do

    end if

    tot1_mmr(:ncol,:) = tot1_mmr(:ncol,:) + tot2_mmr(:ncol,:)

    call self%get_amb_species_numdens(bin_ndx, species_type, rho, aer_numdens) ! #/cm3
    call self%get_ambient_num(bin_ndx, bin_num) ! #/kg

    where(bin_num(:ncol,:)>0._r8)
       awcam(:ncol,:) = ((aer_numdens(:ncol,:)*per_m3/bin_num(:ncol,:)) * tot1_mmr(:ncol,:)) * kg2mug  ! [mug m-3]
    elsewhere
       awcam(:ncol,:) = 0._r8
    end where

    where(tot1_mmr(:ncol,:)>0)
       awfacm(:ncol,:) = tot2_mmr(:ncol,:) / tot1_mmr(:ncol,:)
    elsewhere
       awfacm(:ncol,:) = 0._r8
    end where

    wact_factor(:ncol,:) = awcam(:ncol,:)*(1._r8-awfacm(:ncol,:))

  end subroutine watact_mfactor

  !------------------------------------------------------------------------------
  ! aerosol short wave refactive index
  !------------------------------------------------------------------------------
  function refractive_index_sw(self, ilev, ibin, iwav) result(crefin)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: ilev   ! level index
    integer, intent(in) :: ibin   ! bin index
    integer, intent(in) :: iwav   ! wave length index

    complex(r8) :: crefin(self%ncol_) ! complex refractive index

    real(r8), pointer :: specmmr(:,:) ! species mass mixing ratio
    complex(r8), pointer :: specrefindex(:)     ! species refractive index
    real(r8) :: specdens              ! species density (kg/m3)
    integer :: ispec, icol, ncol
    real(r8) :: vol(self%ncol_)

    ncol = self%ncol()

    crefin(:ncol) = (0._r8, 0._r8)

    do ispec = 1, self%props_%nspecies(ibin)

       call self%ambient_mmr_ptr(species_ndx=ispec, bin_ndx=ibin, mmr=specmmr)
       call self%props_%get(ibin, ispec, density=specdens,  refindex_sw=specrefindex)

       do icol = 1, ncol
          vol(icol) = specmmr(icol,ilev)/specdens
          crefin(icol) = crefin(icol) + vol(icol)*specrefindex(iwav)
       end do
    end do

  end function refractive_index_sw

  !------------------------------------------------------------------------------
  ! aerosol long wave refactive index
  !------------------------------------------------------------------------------
  function refractive_index_lw(self, ilev, ibin, iwav) result(crefin)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: ilev   ! level index
    integer, intent(in) :: ibin   ! bin index
    integer, intent(in) :: iwav   ! wave length index

    complex(r8) :: crefin(self%ncol_) ! complex refractive index

    real(r8), pointer :: specmmr(:,:) ! species mass mixing ratio
    complex(r8), pointer :: specrefindex(:)     ! species refractive index
    real(r8) :: specdens              ! species density (kg/m3)
    integer :: ispec, icol, ncol
    real(r8) :: vol(self%ncol_)

    ncol = self%ncol()

    crefin(:ncol) = (0._r8, 0._r8)

    do ispec = 1, self%props_%nspecies(ibin)

       call self%ambient_mmr_ptr(species_ndx=ispec, bin_ndx=ibin, mmr=specmmr)
       call self%props_%get(ibin, ispec, density=specdens,  refindex_lw=specrefindex)

       do icol = 1, ncol
          vol(icol) = specmmr(icol,ilev)/specdens
          crefin(icol) = crefin(icol) + vol(icol)*specrefindex(iwav)
       end do
    end do

  end function refractive_index_lw

  !------------------------------------------------------------------------------
  ! prescribed aerosol activation fraction for convective cloud
  !------------------------------------------------------------------------------
  function convcld_actfrac(self, ibin, ispc) result(frac)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: ibin   ! bin index
    integer, intent(in) :: ispc   ! species index

    real(r8) :: frac(self%ncol_,self%nlev_)

    frac = 0.8_r8 ! rce 2010/05/02

  end function convcld_actfrac

  !------------------------------------------------------------------------------
  ! below cloud solubility factor for interstitial aerosols
  !------------------------------------------------------------------------------
  function sol_factb_interstitial(self, bin_ndx) result(sol_factb)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number

    real(r8) :: sol_factb(self%ncol_,self%nlev_)

    real(r8), pointer :: aer_mmr(:,:)
    real(r8) :: totmmr(self%ncol_,self%nlev_)
    real(r8) :: solmmr(self%ncol_,self%nlev_)
    integer :: ispc, ncol
    real(r8) :: spechygro

    ncol = self%ncol()

    sol_factb(:,:) = 0.0_r8

    totmmr(:,:) = 0._r8
    solmmr(:,:) = 0._r8

    do ispc = 1, self%props_%nspecies(bin_ndx)

       call self%props_%get(bin_ndx, ispc, hygro=spechygro)
       call self%ambient_mmr_ptr(species_ndx=ispc, bin_ndx=bin_ndx, mmr=aer_mmr)

       totmmr(:ncol,:) = totmmr(:ncol,:) + aer_mmr(:ncol,:)
       solmmr(:ncol,:) = solmmr(:ncol,:) + aer_mmr(:ncol,:)*spechygro

    end do   !nspec

    where ( totmmr > 0._r8 )
       sol_factb = solmmr/totmmr
    end where

  end function sol_factb_interstitial


end module aerosol_state_mod
