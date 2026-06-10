!-----------------------------------------------------------------------------
! Normalized structural description of an aerosol model's population.
!
! aerosol_description_t holds, in one model-agnostic table, the structural
! truth about an aerosol model: which bins exist (MAM modes, CARMA bins,
! individual bulk aerosols), which species each bin holds, the host field
! names, source kinds and resolved host indices of every (bin, species,
! phase) entry, and the field_kind classification table.
!
! A description spans ALL parsed definitions of one aerosol model,
! independent of radiation-list membership; radiation lists select members
! through per-list views. It is concrete and non-polymorphic so it can later
! be metadata-registered and legally cross CCPP cap boundaries in CAM-SIMA.
!
! Descriptions are built once at init time from the already-parsed and
! host-resolved structures in radiative_aerosol_definitions (after host
! index resolution and physprop resolution), and are immutable afterwards.
!
! This module is below aerosol_properties_mod in the aerosol hierarchy and
! will be shared with CAM-SIMA.
!-----------------------------------------------------------------------------
module aerosol_description_mod

  implicit none

  private

  public :: aerosol_description_t
  public :: build_modal_description
  public :: build_carma_description
  public :: build_bulk_description
  public :: bulk_desc_index
  public :: rad_list_view_t
  public :: build_modal_view
  public :: build_carma_view
  public :: build_bulk_view
  public :: field_kind_from_source

  ! Field kind classifications for each (bin, species, phase), where:
  !   bin = aerosol bin or mode (MAM mode, CARMA bin, individual bulk aerosol for BAM)
  !   species = index 0 is the "bin number" field, others are individual species
  !             (e.g., 'sulfate', 'p-organic', 'dust')
  !   phase = ambient, cloud-borne

  ! ADVECTED: (e.g., MAM interstitial mmr and numbers, CARMA interstitial mmr)
  !   Sourced from host constituent.
  !   Read  via pointer or fill.
  !   Write via physics tendencies for constituents.
  integer, public, parameter :: AERO_FIELD_ADVECTED = 1
  ! STORED: (e.g., MAM qqcw)
  !   Sourced from pbuf (CAM) or host non-advected constituent (SIMA).
  !   Read  via pointer or fill.
  !   Write via pointer only.
  integer, public, parameter :: AERO_FIELD_STORED   = 2
  ! DERIVED: (e.g., BAM number concentrations; CARMA number concentrations)
  !   Sourced by deriving on demand from other fields.
  !   Read  via fill only, there is no pointer access.
  !   Write to a derived quantity is meaningless.
  integer, public, parameter :: AERO_FIELD_DERIVED  = 3
  ! ABSENT: (e.g., BAM cloud-borne phase)
  !   This quantity is meaningless for this aerosol model.
  integer, public, parameter :: AERO_FIELD_ABSENT   = 4

  ! phase selectors for field_kind
  integer, public, parameter :: AERO_AMBIENT = 1
  integer, public, parameter :: AERO_CLDBRNE = 2

  ! Per-bin fields are dimensioned (nbins); per-species fields are
  ! dimensioned (nbins, 1:max nspecies); the field_kind table is dimensioned
  ! (nbins, 0:max nspecies, AERO_AMBIENT:AERO_CLDBRNE) with species index 0
  ! selecting the bin number field. Entries that do not exist for a model
  ! hold blank names and index -1.
  type :: aerosol_description_t
     ! number of defined bins (MAM modes / CARMA bins / distinct bulk
     ! aerosols), over all definitions, independent of list membership
     integer :: nbins = 0
     ! mode/bin name as used in the definitions; host field name for bulk
     character(len=64), allocatable :: bin_name(:)
     ! MAM mode type ('accum', 'aitken', ...); blank for other models
     character(len=32), allocatable :: bin_type(:)
     ! number of chemical species in each bin
     integer, allocatable :: nspecies(:)
     ! number of species masses in each bin
     integer, allocatable :: nmasses(:)

     ! bin number field (species index 0): source kind ('A'/'N'/'Z'), host
     ! field name, and resolved host index, per phase
     character(len=1),  allocatable :: source_num_a(:)
     character(len=64), allocatable :: camname_num_a(:)
     integer,           allocatable :: idx_num_a(:)
     character(len=1),  allocatable :: source_num_c(:)
     character(len=64), allocatable :: camname_num_c(:)
     integer,           allocatable :: idx_num_c(:)

     ! bin total-mass field (CARMA bins; absent for other models)
     character(len=1),  allocatable :: source_mass_a(:)
     character(len=64), allocatable :: camname_mass_a(:)
     integer,           allocatable :: idx_mass_a(:)
     character(len=1),  allocatable :: source_mass_c(:)
     character(len=64), allocatable :: camname_mass_c(:)
     integer,           allocatable :: idx_mass_c(:)

     ! species mass mixing ratio fields, per phase
     character(len=1),  allocatable :: source_mmr_a(:,:)
     character(len=64), allocatable :: camname_mmr_a(:,:)
     integer,           allocatable :: idx_mmr_a(:,:)
     character(len=1),  allocatable :: source_mmr_c(:,:)
     character(len=64), allocatable :: camname_mmr_c(:,:)
     integer,           allocatable :: idx_mmr_c(:,:)

     ! species class ('sulfate', 'p-organic', 'dust', ...)
     character(len=32), allocatable :: spec_type(:,:)
     ! species morphology ('shell'/'core' for CARMA; blank for other models)
     character(len=32), allocatable :: spec_morph(:,:)
     ! definition-level per-species physprop id (bin/mode-level per-list
     ! physprop ids live in the radiation list views)
     integer,           allocatable :: spec_physprop_id(:,:)

     ! field kind (AERO_FIELD_*) of each (bin, species, phase) entry
     integer, allocatable :: field_kind(:,:,:)
  end type aerosol_description_t

  ! rad_list_view_t selects the members of one radiation list (climate or
  ! diagnostic) from an aerosol description and carries the per-list data
  ! the description cannot hold: the membership mapping, the per-member
  ! bin/mode-level physprop ids (each list resolves its own physprop file
  ! overrides), and the list-local compressed element indexing.
  type :: rad_list_view_t
     ! radiation list index (0=climate, 1..N_DIAG diagnostics)
     integer :: list_idx = 0
     ! number of member bins in this list
     integer :: nbins = 0
     ! member -> description bin index
     integer, allocatable :: bin_idx(:)
     ! member -> bin/mode-level physprop id resolved for this list
     integer, allocatable :: physprop_id(:)
     ! total number of constituents (elements) across the list members
     integer :: ncnst_tot = 0
     ! list-local compressed (bin, 0:nmasses) -> element index map
     integer, allocatable :: indexer(:,:)
  end type rad_list_view_t

!==============================================================================
contains
!==============================================================================

  !----------------------------------------------------------------------------
  ! builds the modal (MAM) aerosol description from the parsed and resolved
  ! mode definitions
  !----------------------------------------------------------------------------
  subroutine build_modal_description(modes, desc)
    use radiative_aerosol_definitions, only: modes_t

    type(modes_t),               intent(in)  :: modes
    type(aerosol_description_t), intent(out) :: desc

    integer :: m, l, nspec_max

    nspec_max = 0
    do m = 1, modes%nmodes
       nspec_max = max(nspec_max, modes%comps(m)%nspec)
    end do

    call alloc_description(desc, modes%nmodes, nspec_max, 'build_modal_description')

    do m = 1, modes%nmodes
       desc%bin_name(m) = modes%names(m)
       desc%bin_type(m) = modes%types(m)
       desc%nspecies(m) = modes%comps(m)%nspec
       desc%nmasses(m)  = modes%comps(m)%nspec

       desc%source_num_a(m)  = modes%comps(m)%source_num_a
       desc%camname_num_a(m) = modes%comps(m)%camname_num_a
       desc%idx_num_a(m)     = modes%comps(m)%idx_num_a
       desc%source_num_c(m)  = modes%comps(m)%source_num_c
       desc%camname_num_c(m) = modes%comps(m)%camname_num_c
       desc%idx_num_c(m)     = modes%comps(m)%idx_num_c

       ! field kind table built from the parsed per-entry source data
       ! ('A' advected constituent, 'N' pbuf-resident); no modal overrides
       desc%field_kind(m,0,AERO_AMBIENT) = field_kind_from_source(modes%comps(m)%source_num_a)
       desc%field_kind(m,0,AERO_CLDBRNE) = field_kind_from_source(modes%comps(m)%source_num_c)

       do l = 1, modes%comps(m)%nspec
          desc%source_mmr_a(m,l)  = modes%comps(m)%source_mmr_a(l)
          desc%camname_mmr_a(m,l) = modes%comps(m)%camname_mmr_a(l)
          desc%idx_mmr_a(m,l)     = modes%comps(m)%idx_mmr_a(l)
          desc%source_mmr_c(m,l)  = modes%comps(m)%source_mmr_c(l)
          desc%camname_mmr_c(m,l) = modes%comps(m)%camname_mmr_c(l)
          desc%idx_mmr_c(m,l)     = modes%comps(m)%idx_mmr_c(l)

          desc%spec_type(m,l)        = modes%comps(m)%type(l)
          desc%spec_physprop_id(m,l) = modes%comps(m)%idx_props(l)

          desc%field_kind(m,l,AERO_AMBIENT) = field_kind_from_source(modes%comps(m)%source_mmr_a(l))
          desc%field_kind(m,l,AERO_CLDBRNE) = field_kind_from_source(modes%comps(m)%source_mmr_c(l))
       end do
    end do

  end subroutine build_modal_description

  !----------------------------------------------------------------------------
  ! builds the sectional (CARMA) aerosol description from the parsed and
  ! resolved bin definitions
  !----------------------------------------------------------------------------
  subroutine build_carma_description(bins, desc)
    use radiative_aerosol_definitions, only: bins_t

    type(bins_t),                intent(in)  :: bins
    type(aerosol_description_t), intent(out) :: desc

    integer :: m, l, nspec_max

    nspec_max = 0
    do m = 1, bins%nbins
       nspec_max = max(nspec_max, bins%comps(m)%nspec)
    end do

    call alloc_description(desc, bins%nbins, nspec_max, 'build_carma_description')

    do m = 1, bins%nbins
       desc%bin_name(m) = bins%names(m)
       desc%nspecies(m) = bins%comps(m)%nspec
       desc%nmasses(m)  = bins%comps(m)%nspec

       desc%source_num_a(m)  = bins%comps(m)%source_num_a
       desc%camname_num_a(m) = bins%comps(m)%camname_num_a
       desc%idx_num_a(m)     = bins%comps(m)%idx_num_a
       desc%source_num_c(m)  = bins%comps(m)%source_num_c
       desc%camname_num_c(m) = bins%comps(m)%camname_num_c
       desc%idx_num_c(m)     = bins%comps(m)%idx_num_c

       ! the bin total-mass channel is optional; its host index is resolved
       ! only when a channel was given in the bin definition ('NOTSET'
       ! placeholders otherwise), so it is copied under the same condition
       desc%source_mass_a(m)  = bins%comps(m)%source_mass_a
       desc%camname_mass_a(m) = bins%comps(m)%camname_mass_a
       if (bins%comps(m)%source_mass_a /= 'NOTSET' .and. bins%comps(m)%camname_mass_a /= 'NOTSET') then
          desc%idx_mass_a(m) = bins%comps(m)%idx_mass_a
       end if
       desc%source_mass_c(m)  = bins%comps(m)%source_mass_c
       desc%camname_mass_c(m) = bins%comps(m)%camname_mass_c
       if (bins%comps(m)%source_mass_c /= 'NOTSET' .and. bins%comps(m)%camname_mass_c /= 'NOTSET') then
          desc%idx_mass_c(m) = bins%comps(m)%idx_mass_c
       end if

       ! field kind table built from the parsed per-entry source data:
       !  A - advected constituent
       !  N - pbuf (CAM), non-advected constituent (SIMA)
       !
       ! Special note for CARMA: Bin numbers for CARMA in both ambient and
       ! cloud-borne phases are derived from bin mass via the fixed bin radius
       ! and not stored in the host model, so it is a "derived" quantity.
       desc%field_kind(m,0,AERO_AMBIENT) = AERO_FIELD_DERIVED
       desc%field_kind(m,0,AERO_CLDBRNE) = AERO_FIELD_DERIVED

       do l = 1, bins%comps(m)%nspec
          desc%source_mmr_a(m,l)  = bins%comps(m)%source_mmr_a(l)
          desc%camname_mmr_a(m,l) = bins%comps(m)%camname_mmr_a(l)
          desc%idx_mmr_a(m,l)     = bins%comps(m)%idx_mmr_a(l)
          desc%source_mmr_c(m,l)  = bins%comps(m)%source_mmr_c(l)
          desc%camname_mmr_c(m,l) = bins%comps(m)%camname_mmr_c(l)
          desc%idx_mmr_c(m,l)     = bins%comps(m)%idx_mmr_c(l)

          desc%spec_type(m,l)        = bins%comps(m)%type(l)
          desc%spec_morph(m,l)       = bins%comps(m)%morph(l)
          desc%spec_physprop_id(m,l) = bins%comps(m)%idx_props(l)

          desc%field_kind(m,l,AERO_AMBIENT) = field_kind_from_source(bins%comps(m)%source_mmr_a(l))
          desc%field_kind(m,l,AERO_CLDBRNE) = field_kind_from_source(bins%comps(m)%source_mmr_c(l))
       end do
    end do

  end subroutine build_carma_description

  !----------------------------------------------------------------------------
  ! builds the bulk (BAM) aerosol description.
  !
  ! Bulk aerosols have no separate definition level: each radiation list
  ! enumerates its aerosols directly. The description is the union of the
  ! distinct (source, host field name) entries across all active lists, in
  ! first-seen order (climate list first); each distinct entry is one
  ! single-species bin.
  !----------------------------------------------------------------------------
  subroutine build_bulk_description(aerlists, active, desc)
    use radiative_aerosol_definitions, only: aerlist_t
    use cam_abortutils,                only: endrun

    type(aerlist_t),             intent(in)  :: aerlists(0:)  ! climate + diagnostic bulk aerosol lists
    logical,                     intent(in)  :: active(0:)    ! which lists are active
    type(aerosol_description_t), intent(out) :: desc

    integer :: i, ilist, n, ntot, istat
    character(len=1),  allocatable :: usrc(:)
    character(len=64), allocatable :: uname(:)
    integer,           allocatable :: uidx(:)
    character(len=*), parameter :: subname = 'build_bulk_description'

    ! collect the distinct (source, camname) entries across active lists;
    ! the resolved host index is determined by (source, camname), so the
    ! first occurrence carries it
    ntot = 0
    do ilist = 0, ubound(aerlists, 1)
       if (active(ilist)) ntot = ntot + aerlists(ilist)%numaerosols
    end do

    allocate(usrc(ntot), uname(ntot), uidx(ntot), stat=istat)
    if (istat /= 0) call endrun(subname//': allocation error for bulk union scan')

    n = 0
    do ilist = 0, ubound(aerlists, 1)
       if (.not. active(ilist)) cycle
       do i = 1, aerlists(ilist)%numaerosols
          if (find_entry(usrc, uname, n, aerlists(ilist)%aer(i)%source, &
                         aerlists(ilist)%aer(i)%camname) == -1) then
             n = n + 1
             usrc(n)  = aerlists(ilist)%aer(i)%source
             uname(n) = aerlists(ilist)%aer(i)%camname
             uidx(n)  = aerlists(ilist)%aer(i)%idx
          end if
       end do
    end do

    call alloc_description(desc, n, 1, subname)

    ! field kind table: ambient mass from the parsed per-entry source data
    ! ('A' advected constituent, 'N' pbuf (CAM)/non-advected (SIMA),
    !  e.g. prescribed aerosol), with the model overrides that the BAM ambient
    ! number is derived from mass (via num_to_mass) and BAM has no cloud-borne
    ! phase.
    !
    ! species_ndx 0 refers to the bin's number field and 1..nspecies(bin) to
    ! its mass species (same convention as the species indexer). Each bulk
    ! "bin" holds exactly one species, so index 0 is its number concentration
    ! -- always derived from mass via num_to_mass -- and index 1 is its mass
    ! mmr, classified from the rad list source descriptor.
    !
    ! spec_type is intentionally left blank: bulk aerosol definitions carry no
    ! species class; it is derived from the physprop file's aerosol name by
    ! the bulk_aerosol_properties species_type override, not from this table.
    do i = 1, n
       desc%bin_name(i) = uname(i)
       desc%nspecies(i) = 1
       desc%nmasses(i)  = 1

       desc%source_mmr_a(i,1)  = usrc(i)
       desc%camname_mmr_a(i,1) = uname(i)
       desc%idx_mmr_a(i,1)     = uidx(i)

       desc%field_kind(i,0,AERO_AMBIENT) = AERO_FIELD_DERIVED
       desc%field_kind(i,1,AERO_AMBIENT) = field_kind_from_source(usrc(i))
       desc%field_kind(i,:,AERO_CLDBRNE) = AERO_FIELD_ABSENT
    end do

    deallocate(usrc, uname, uidx)

  end subroutine build_bulk_description

  !----------------------------------------------------------------------------
  ! returns the bulk description bin index holding the given (source, host
  ! field name) entry, or -1 if not present
  !----------------------------------------------------------------------------
  pure integer function bulk_desc_index(desc, source, camname)
    type(aerosol_description_t), intent(in) :: desc
    character(len=*),            intent(in) :: source
    character(len=*),            intent(in) :: camname

    integer :: i

    bulk_desc_index = -1
    do i = 1, desc%nbins
       if (desc%source_mmr_a(i,1) == source .and. desc%camname_mmr_a(i,1) == camname) then
          bulk_desc_index = i
          return
       end if
    end do
  end function bulk_desc_index

  !----------------------------------------------------------------------------
  ! builds the view of a radiation list over the modal description from the
  ! parsed per-list mode list
  !----------------------------------------------------------------------------
  subroutine build_modal_view(desc, mlist, list_idx, view)
    use radiative_aerosol_definitions, only: modelist_t

    type(aerosol_description_t), intent(in)  :: desc
    type(modelist_t),            intent(in)  :: mlist
    integer,                     intent(in)  :: list_idx
    type(rad_list_view_t),       intent(out) :: view

    integer :: m, istat
    character(len=*), parameter :: subname = 'build_modal_view'

    view%list_idx = list_idx
    view%nbins    = mlist%nmodes

    allocate(view%bin_idx(view%nbins), view%physprop_id(view%nbins), stat=istat)
    if (istat /= 0) call endrun_alloc(subname)

    view%bin_idx(:)     = mlist%idx(1:view%nbins)
    view%physprop_id(:) = mlist%idx_props(1:view%nbins)

    ! one number element plus the species mass elements per member
    ! (matches the modal properties constructor)
    view%ncnst_tot = 0
    do m = 1, view%nbins
       view%ncnst_tot = view%ncnst_tot + desc%nspecies(view%bin_idx(m)) + 1
    end do

    call build_view_indexer(desc, view, subname)

  end subroutine build_modal_view

  !----------------------------------------------------------------------------
  ! builds the view of a radiation list over the sectional (CARMA)
  ! description from the parsed per-list bin list
  !----------------------------------------------------------------------------
  subroutine build_carma_view(desc, slist, list_idx, view)
    use radiative_aerosol_definitions, only: binlist_t

    type(aerosol_description_t), intent(in)  :: desc
    type(binlist_t),             intent(in)  :: slist
    integer,                     intent(in)  :: list_idx
    type(rad_list_view_t),       intent(out) :: view

    integer :: m, istat
    character(len=*), parameter :: subname = 'build_carma_view'

    view%list_idx = list_idx
    view%nbins    = slist%nbins

    allocate(view%bin_idx(view%nbins), view%physprop_id(view%nbins), stat=istat)
    if (istat /= 0) call endrun_alloc(subname)

    view%bin_idx(:)     = slist%idx(1:view%nbins)
    view%physprop_id(:) = slist%idx_props(1:view%nbins)

    ! one number element plus the species mass elements per member
    ! (matches the CARMA properties constructor)
    view%ncnst_tot = 0
    do m = 1, view%nbins
       view%ncnst_tot = view%ncnst_tot + desc%nspecies(view%bin_idx(m)) + 1
    end do

    call build_view_indexer(desc, view, subname)

  end subroutine build_carma_view

  !----------------------------------------------------------------------------
  ! builds the view of a radiation list over the bulk description from the
  ! parsed per-list bulk aerosol list
  !----------------------------------------------------------------------------
  subroutine build_bulk_view(desc, alist, list_idx, view)
    use radiative_aerosol_definitions, only: aerlist_t
    use cam_abortutils,                only: endrun

    type(aerosol_description_t), intent(in)  :: desc
    type(aerlist_t),             intent(in)  :: alist
    integer,                     intent(in)  :: list_idx
    type(rad_list_view_t),       intent(out) :: view

    integer :: i, istat
    character(len=*), parameter :: subname = 'build_bulk_view'

    view%list_idx = list_idx
    view%nbins    = alist%numaerosols

    allocate(view%bin_idx(view%nbins), view%physprop_id(view%nbins), stat=istat)
    if (istat /= 0) call endrun_alloc(subname)

    do i = 1, view%nbins
       view%bin_idx(i) = bulk_desc_index(desc, alist%aer(i)%source, alist%aer(i)%camname)
       if (view%bin_idx(i) == -1) then
          call endrun(subname//': bulk aerosol not in description: '//trim(alist%aer(i)%camname))
       end if
       view%physprop_id(i) = alist%aer(i)%physprop_id
    end do

    ! For bulk aerosols, the number of bins and total number of constituents
    ! are the same -- one constituent (species and mass) per bin (matches the
    ! bulk properties constructor; the indexer below nevertheless carries a
    ! number slot per bin, preserving the existing element numbering).
    view%ncnst_tot = alist%numaerosols

    call build_view_indexer(desc, view, subname)

  end subroutine build_bulk_view

  !----------------------------------------------------------------------------
  ! builds the list-local compressed element index map of a view.
  !
  ! Local indexing compresses the bin and number/mass indices into one index.
  ! This indexing is used by the pointer arrays used to reference state and
  ! pbuf fields. We add number = 0, total mass = 1 (if available), and mass
  ! from each constituency into mm.
  !----------------------------------------------------------------------------
  subroutine build_view_indexer(desc, view, caller)
    type(aerosol_description_t), intent(in)    :: desc
    type(rad_list_view_t),       intent(inout) :: view
    character(len=*),            intent(in)    :: caller

    integer :: ibin, imas, indx, istat

    allocate(view%indexer(view%nbins, 0:maxval(desc%nmasses(view%bin_idx))), stat=istat)
    if (istat /= 0) call endrun_alloc(caller)

    view%indexer = -1
    indx = 0

    do ibin = 1, view%nbins
       do imas = 0, desc%nmasses(view%bin_idx(ibin))
          indx = indx + 1
          view%indexer(ibin,imas) = indx
       end do
    end do

  end subroutine build_view_indexer

  !----------------------------------------------------------------------------
  ! aborts on a view allocation failure
  !----------------------------------------------------------------------------
  subroutine endrun_alloc(caller)
    use cam_abortutils, only: endrun

    character(len=*), intent(in) :: caller

    call endrun(caller//': allocation error for radiation list view')
  end subroutine endrun_alloc

  !----------------------------------------------------------------------------
  ! maps a host aerosol-list source descriptor onto a field kind:
  ! 'A' state constituent (advected); 'N' pbuf-resident (non-advected);
  ! 'Z' zero / blank -- no such field for this model
  !----------------------------------------------------------------------------
  integer function field_kind_from_source(source)
    use cam_abortutils, only: endrun

    character(len=*), intent(in) :: source

    select case (source(1:1))
    case ('A')
       field_kind_from_source = AERO_FIELD_ADVECTED
    case ('N')
       field_kind_from_source = AERO_FIELD_STORED
    case ('Z', ' ')
       field_kind_from_source = AERO_FIELD_ABSENT
    case default
       field_kind_from_source = AERO_FIELD_ABSENT
       call endrun('field_kind_from_source: unrecognized source: '//source)
    end select
  end function field_kind_from_source

  !----------------------------------------------------------------------------
  ! allocates and default-initializes all description tables: blank names,
  ! index -1, field kind -1
  !----------------------------------------------------------------------------
  subroutine alloc_description(desc, nbins, nspec_max, caller)
    use cam_abortutils, only: endrun

    type(aerosol_description_t), intent(inout) :: desc
    integer,                     intent(in)    :: nbins
    integer,                     intent(in)    :: nspec_max
    character(len=*),            intent(in)    :: caller

    integer :: istat

    desc%nbins = nbins

    allocate(desc%bin_name(nbins), desc%bin_type(nbins), &
             desc%nspecies(nbins), desc%nmasses(nbins), &
             desc%source_num_a(nbins), desc%camname_num_a(nbins), desc%idx_num_a(nbins), &
             desc%source_num_c(nbins), desc%camname_num_c(nbins), desc%idx_num_c(nbins), &
             desc%source_mass_a(nbins), desc%camname_mass_a(nbins), desc%idx_mass_a(nbins), &
             desc%source_mass_c(nbins), desc%camname_mass_c(nbins), desc%idx_mass_c(nbins), &
             desc%source_mmr_a(nbins,nspec_max), desc%camname_mmr_a(nbins,nspec_max), &
             desc%idx_mmr_a(nbins,nspec_max), &
             desc%source_mmr_c(nbins,nspec_max), desc%camname_mmr_c(nbins,nspec_max), &
             desc%idx_mmr_c(nbins,nspec_max), &
             desc%spec_type(nbins,nspec_max), desc%spec_morph(nbins,nspec_max), &
             desc%spec_physprop_id(nbins,nspec_max), &
             desc%field_kind(nbins,0:nspec_max,AERO_AMBIENT:AERO_CLDBRNE), &
             stat=istat)
    if (istat /= 0) call endrun(caller//': allocation error for aerosol description')

    desc%bin_name(:) = ' '
    desc%bin_type(:) = ' '
    desc%nspecies(:) = 0
    desc%nmasses(:)  = 0

    desc%source_num_a(:)  = ' '
    desc%camname_num_a(:) = ' '
    desc%idx_num_a(:)     = -1
    desc%source_num_c(:)  = ' '
    desc%camname_num_c(:) = ' '
    desc%idx_num_c(:)     = -1

    desc%source_mass_a(:)  = ' '
    desc%camname_mass_a(:) = ' '
    desc%idx_mass_a(:)     = -1
    desc%source_mass_c(:)  = ' '
    desc%camname_mass_c(:) = ' '
    desc%idx_mass_c(:)     = -1

    desc%source_mmr_a(:,:)  = ' '
    desc%camname_mmr_a(:,:) = ' '
    desc%idx_mmr_a(:,:)     = -1
    desc%source_mmr_c(:,:)  = ' '
    desc%camname_mmr_c(:,:) = ' '
    desc%idx_mmr_c(:,:)     = -1

    desc%spec_type(:,:)        = ' '
    desc%spec_morph(:,:)       = ' '
    desc%spec_physprop_id(:,:) = -1

    desc%field_kind(:,:,:) = -1

  end subroutine alloc_description

  !----------------------------------------------------------------------------
  ! returns the position of (source, camname) within the first n entries of
  ! the union scan arrays, or -1 if not present
  !----------------------------------------------------------------------------
  pure integer function find_entry(usrc, uname, n, source, camname)
    character(len=*), intent(in) :: usrc(:)
    character(len=*), intent(in) :: uname(:)
    integer,          intent(in) :: n
    character(len=*), intent(in) :: source
    character(len=*), intent(in) :: camname

    integer :: i

    find_entry = -1
    do i = 1, n
       if (usrc(i) == source .and. uname(i) == camname) then
          find_entry = i
          return
       end if
    end do
  end function find_entry

end module aerosol_description_mod
