module aerosol_instances_mod
  ! aerosol_instances_mod owns and manages the concrete aerosol_properties and
  ! aerosol_state objects for every active aerosol model (modal, CARMA, bulk)
  ! and every radiation list (climate + diagnostics).
  !
  ! Lifecycle (CAM host model example):
  !   1. aerosol_instances_init()        -- called once during phys_init, after
  !      rad_aer_init().  Creates persistent aerosol_properties objects for
  !      each (aerosol_model, list_idx) pair.
  !   2. aerosol_instances_init_states() -- called once during phys_init, after
  !      aerosol_instances_init().  Creates persistent aerosol_state objects
  !      for each (aerosol_model, list_idx, chunk) triple.  States store
  !      pointers to phys_state(c) and pbuf, which persist for the run.
  !   3. aerosol_instances_get_props()   -- returns a pointer to a properties
  !      object for a given (aerosol_model, list_idx).
  !   4. aerosol_instances_get_state()   -- returns a pointer to a state
  !      object for a given (aerosol_model, list_idx, chunk).
  !   5. aerosol_instances_final()       -- deallocates all objects at shutdown.
  !
  ! NOTE: when using _get_props() and _get_state(), they return pointers to
  ! persistent properties and state objects. The caller should only nullify their
  ! copy of the pointers after use, but should NOT deallocate as these objects
  ! are shared per (aerosol_model, list_idx, chunk).
  !
  ! For transient state (e.g., bound to a local copy of physics_state),
  ! aerosol_instances_create_states / destroy_states provide a per-call factory,
  ! but this capability will be removed after CAM is retired and intermediate
  ! state tendencies will be replaced by separate physics scheme.
  !
  ! The init, get_props, get_state, and final routines are portable.
  ! The create/destroy_states factory and init_states are host-model specific
  ! as they point to host-model specific data structures for aerosol state info.

  use aerosol_properties_mod,        only: aerosol_properties
  use aerosol_state_mod,             only: aerosol_state
  use aerosol_description_mod,       only: aerosol_description_t, rad_list_view_t
  use radiative_aerosol_definitions, only: N_DIAG

  implicit none
  private

  public :: aerosol_instances_init
  public :: aerosol_instances_init_states
  public :: aerosol_instances_get_props
  public :: aerosol_instances_get_props_by_model
  public :: aerosol_instances_get_description
  public :: aerosol_instances_get_state
  public :: aerosol_instances_get_num_models
  public :: aerosol_instances_is_active
  public :: aerosol_instances_final
  public :: aerosol_instances_create_states
  public :: aerosol_instances_destroy_states
  public :: aero_state_entry_t

  type :: aero_props_entry_t
     class(aerosol_properties), pointer :: obj => null()
  end type aero_props_entry_t

  type :: aero_state_entry_t
     class(aerosol_state),      pointer :: obj => null()
  end type aero_state_entry_t

  ! Normalized structural descriptions, one per active aerosol model,
  ! dimensioned (iaermod); built once at aerosol_instances_init and
  ! immutable afterwards.
  type(aerosol_description_t), allocatable, target :: aero_desc_all(:)

  ! Per-list views over the descriptions, dimensioned (iaermod, 0:N_DIAG);
  ! built once at aerosol_instances_init. Entries for (model, list) pairs
  ! with no members remain empty (nbins=0).
  type(rad_list_view_t), allocatable, target :: aero_views_all(:,:)

  ! Persistent aerosol properties objects
  ! dimensioned (iaermod, 0:N_DIAG).
  type(aero_props_entry_t), allocatable, target :: aero_props_all(:,:)

  ! Persistent per-chunk aerosol state objects
  ! dimensioned (iaermod, 0:N_DIAG, begchunk:endchunk).
  ! States store pointers to phys_state(c) and pbuf which persist for the run.
  type(aero_state_entry_t), allocatable, target :: aero_states_all(:,:,:)

  ! Number of aerosol models active at runtime.
  ! Note: Multiple aerosol models can be active at once.
  ! e.g., using bulk for volcanic aerosol and modal for others.
  ! When retrieving properties via aerosol_instances_get_props,
  ! or creating states from aerosol_instances_create_states,
  ! ensure that the aerosol model matches what is needed
  ! (e.g., aero_props%model_is('MAM') == .true.)
  integer :: num_aero_models_ = 0

  logical :: modal_active_ = .false.
  logical :: carma_active_ = .false.
  logical :: bulk_active_  = .false.

contains
  ! Determine which aerosol models are active (modal, CARMA, bulk) and
  ! create persistent aerosol_properties objects for each (model, list) pair.
  ! Must be called after the radiative aerosol module
  ! has parsed aerosol definitions.
  !
  ! NOTE: A model is "globally active" if the climate list (list 0) has > 0
  ! entries for aerosol with that representation, but individual diagnostic
  ! lists may have zero entries. In that case the corresponding properties
  ! slot is left null.
  !
  ! Callers that use lists other than the climate list thus need to check
  ! if the aerosol model for that combination of (aero_model, list)
  ! is associated or not.
  subroutine aerosol_instances_init()
    use radiative_aerosol, only: rad_aer_get_info
    use radiative_aerosol_definitions, only: active_calls, modes, bins, bulk_aerosol_list, &
                                             modal_aerosol_list, sectional_aerosol_list
    use aerosol_description_mod, only: build_modal_description, build_carma_description, &
                                       build_bulk_description, &
                                       build_modal_view, build_carma_view, build_bulk_view
    use modal_aerosol_properties_mod, only: modal_aerosol_properties
    use carma_aerosol_properties_mod, only: carma_aerosol_properties
    use bulk_aerosol_properties_mod,  only: bulk_aerosol_properties
    use cam_abortutils, only: endrun

    use spmd_utils, only: masterproc
    use cam_logfile, only: iulog

    integer :: nmodes, nbins, nbulk_aerosols
    integer :: iaermod, ilist, istat

    character(len=*), parameter :: subname = 'aerosol_instances_init: '

    num_aero_models_ = 0

    call rad_aer_get_info(0, nmodes=nmodes, nbins=nbins, naero=nbulk_aerosols)
    modal_active_ = nmodes > 0
    carma_active_ = nbins > 0
    bulk_active_  = nbulk_aerosols > 0

    if (masterproc) then
       write(iulog,*) subname,'nmodes,nbins,nbulk_aerosols: ',nmodes,nbins,nbulk_aerosols
    end if

    if (modal_active_) num_aero_models_ = num_aero_models_ + 1
    if (carma_active_) num_aero_models_ = num_aero_models_ + 1
    if (bulk_active_)  num_aero_models_ = num_aero_models_ + 1

    if (num_aero_models_ < 1) return

    ! Build the normalized structural description for each active model from
    ! the parsed and host-resolved definitions.
    allocate(aero_desc_all(num_aero_models_), stat=istat)
    if (istat /= 0) then
       call endrun(subname//'allocation error: aero_desc_all')
    end if

    iaermod = 0
    if (modal_active_) then
       iaermod = iaermod + 1
       call build_modal_description(modes, aero_desc_all(iaermod))
    end if
    if (carma_active_) then
       iaermod = iaermod + 1
       call build_carma_description(bins, aero_desc_all(iaermod))
    end if
    if (bulk_active_) then
       iaermod = iaermod + 1
       call build_bulk_description(bulk_aerosol_list, active_calls, aero_desc_all(iaermod))
    end if

    allocate(aero_props_all(num_aero_models_, 0:N_DIAG), stat=istat)
    if (istat /= 0) then
       call endrun(subname//'allocation error: aero_props_all')
    end if

    allocate(aero_views_all(num_aero_models_, 0:N_DIAG), stat=istat)
    if (istat /= 0) then
       call endrun(subname//'allocation error: aero_views_all')
    end if

    do ilist = 0, N_DIAG
       ! only populate views and aerosol properties for active
       ! climate/diagnostic lists.
       if (.not. active_calls(ilist)) cycle

       call rad_aer_get_info(ilist, nmodes=nmodes, nbins=nbins, naero=nbulk_aerosols)

       iaermod = 0
       if (modal_active_) then
          iaermod = iaermod + 1
          if (nmodes > 0) then
             call build_modal_view(aero_desc_all(iaermod), modal_aerosol_list(ilist), &
                                   ilist, aero_views_all(iaermod, ilist))
             aero_props_all(iaermod, ilist)%obj => &
                  modal_aerosol_properties(aero_desc_all(iaermod), aero_views_all(iaermod, ilist))
          end if
       end if
       if (carma_active_) then
          iaermod = iaermod + 1
          if (nbins > 0) then
             call build_carma_view(aero_desc_all(iaermod), sectional_aerosol_list(ilist), &
                                   ilist, aero_views_all(iaermod, ilist))
             aero_props_all(iaermod, ilist)%obj => &
                  carma_aerosol_properties(aero_desc_all(iaermod), aero_views_all(iaermod, ilist))
          end if
       end if
       if (bulk_active_) then
          iaermod = iaermod + 1
          if (nbulk_aerosols > 0) then
             call build_bulk_view(aero_desc_all(iaermod), bulk_aerosol_list(ilist), &
                                  ilist, aero_views_all(iaermod, ilist))
             aero_props_all(iaermod, ilist)%obj => &
                  bulk_aerosol_properties(aero_desc_all(iaermod), aero_views_all(iaermod, ilist))
          end if
       end if
    end do

    ! Consistency check: the field_kind (advected/stored/derived) flags must
    ! agree with the constituent index resolution.
    do ilist = 0, N_DIAG
       do iaermod = 1, num_aero_models_
          if (associated(aero_props_all(iaermod, ilist)%obj)) then
             call check_field_kinds(aero_props_all(iaermod, ilist)%obj)
          end if
       end do
    end do

  end subroutine aerosol_instances_init

  ! Init-time consistency check of the field_kind table for models that
  ! support the working-state table (MAM, CARMA):
  ! (1) an ambient entry must be ADVECTED exactly when its interstitial constituent
  ! name resolves to a constituent index (updated via ptend)
  ! (2) cloud-borne entries must never be ADVECTED (updated in-place via pointer)
  subroutine check_field_kinds(props)
    use aerosol_properties_mod, only: aero_name_len, aero_has_working_state_table
    use aerosol_properties_mod, only: AERO_FIELD_ADVECTED, AERO_FIELD_STORED
    use aerosol_properties_mod, only: AERO_AMBIENT, AERO_CLDBRNE
    use constituents,   only: cnst_get_ind
    use cam_abortutils, only: endrun

    class(aerosol_properties), intent(in) :: props

    character(len=aero_name_len) :: name_a, name_c
    integer :: m, l, ndx
    character(len=*), parameter :: subname = 'aerosol_instances_init: '

    if (.not. props%supports(aero_has_working_state_table)) return

    do m = 1, props%nbins()
       do l = 0, props%nspecies(m)
          if (l == 0) then
             call props%num_names(m, name_a, name_c)
          else
             call props%mmr_names(m, l, name_a, name_c)
          end if
          call cnst_get_ind(trim(name_a), ndx, abort=.false.)
          if ((ndx > 0) .neqv. (props%field_kind(m, l, AERO_AMBIENT) == AERO_FIELD_ADVECTED)) then
             call endrun(subname//'field_kind inconsistent with constituent lookup for '//trim(name_a))
          end if
          if (props%field_kind(m, l, AERO_AMBIENT) == AERO_FIELD_STORED) then
             call endrun(subname//'unexpected STORED ambient entry: '//trim(name_a))
          end if
          ! cloud-borne entries must never be ADVECTED: working-state consumers
          ! (e.g. dropmixnuc) update them in place unconditionally, which is
          ! legal only for non-constituent fields
          if (props%field_kind(m, l, AERO_CLDBRNE) == AERO_FIELD_ADVECTED) then
             call endrun(subname//'unexpected ADVECTED cloud-borne entry: '//trim(name_c))
          end if
       end do
    end do

  end subroutine check_field_kinds

  ! Return a pointer to the aerosol_properties object for the given aerosol
  ! model index and radiation list.  Returns null when the model has no
  ! entries in the specified list (see aerosol_instances_init).
  function aerosol_instances_get_props(iaermod, list_idx) result(props)
    integer, intent(in) :: iaermod   ! aerosol model index (1..num_aero_models)
    integer, intent(in) :: list_idx  ! radiation list index (0=climate, 1..N_DIAG)
    class(aerosol_properties), pointer :: props

    props => aero_props_all(iaermod, list_idx)%obj

  end function aerosol_instances_get_props

  ! Return a pointer to the aerosol_properties object for the given aerosol
  ! model name ('modal', 'carma', or 'bulk') and radiation list.  Returns
  ! null when the model is inactive or has no entries in the specified list.
  function aerosol_instances_get_props_by_model(model_name, list_idx) result(props)
    character(len=*), intent(in) :: model_name ! 'modal', 'carma', or 'bulk'
    integer,          intent(in) :: list_idx   ! radiation list index (0=climate, 1..N_DIAG)
    class(aerosol_properties), pointer :: props

    integer :: iaermod

    nullify(props)

    iaermod = 0
    if (modal_active_) then
       iaermod = iaermod + 1
       if (trim(model_name) == 'modal') then
          props => aero_props_all(iaermod, list_idx)%obj
          return
       end if
    end if
    if (carma_active_) then
       iaermod = iaermod + 1
       if (trim(model_name) == 'carma') then
          props => aero_props_all(iaermod, list_idx)%obj
          return
       end if
    end if
    if (bulk_active_) then
       iaermod = iaermod + 1
       if (trim(model_name) == 'bulk') then
          props => aero_props_all(iaermod, list_idx)%obj
          return
       end if
    end if

  end function aerosol_instances_get_props_by_model

  ! Return a pointer to the normalized structural description for the given
  ! aerosol model index.
  function aerosol_instances_get_description(iaermod) result(desc)
    integer, intent(in) :: iaermod   ! aerosol model index (1..num_aero_models)
    type(aerosol_description_t), pointer :: desc

    desc => aero_desc_all(iaermod)

  end function aerosol_instances_get_description

  ! Return the number of aerosol models active at runtime.
  pure integer function aerosol_instances_get_num_models()
    aerosol_instances_get_num_models = num_aero_models_
  end function aerosol_instances_get_num_models

  ! Return .true. if the named aerosol model ('modal', 'carma', or 'bulk')
  ! has any entries in the climate list.
  logical function aerosol_instances_is_active(model_name)
    character(len=*), intent(in) :: model_name

    select case (trim(model_name))
    case ('modal')
       aerosol_instances_is_active = modal_active_
    case ('carma')
       aerosol_instances_is_active = carma_active_
    case ('bulk')
       aerosol_instances_is_active = bulk_active_
    case default
       aerosol_instances_is_active = .false.
    end select

  end function aerosol_instances_is_active

  subroutine aerosol_instances_final()
    use ppgrid, only: begchunk, endchunk
    integer :: iaermod, ilist, c

    ! Deallocate persistent state objects
    if (allocated(aero_states_all)) then
       do c = begchunk, endchunk
          do ilist = 0, N_DIAG
             do iaermod = 1, num_aero_models_
                if (associated(aero_states_all(iaermod, ilist, c)%obj)) then
                   deallocate(aero_states_all(iaermod, ilist, c)%obj)
                   nullify(aero_states_all(iaermod, ilist, c)%obj)
                end if
             end do
          end do
       end do
       deallocate(aero_states_all)
    end if

    ! Deallocate properties objects
    if (allocated(aero_props_all)) then
       do ilist = 0, N_DIAG
          do iaermod = 1, num_aero_models_
             if (associated(aero_props_all(iaermod, ilist)%obj)) then
                deallocate(aero_props_all(iaermod, ilist)%obj)
                nullify(aero_props_all(iaermod, ilist)%obj)
             end if
          end do
       end do
       deallocate(aero_props_all)
    end if

    ! Deallocate views and descriptions (allocatable components are released
    ! with them)
    if (allocated(aero_views_all)) then
       deallocate(aero_views_all)
    end if
    if (allocated(aero_desc_all)) then
       deallocate(aero_desc_all)
    end if

    num_aero_models_ = 0

  end subroutine aerosol_instances_final

  ! Initialize persistent per-chunk aerosol state objects for all active lists
  ! and all active aerosol models.
  !
  ! Called once at init time, after aerosol_instances_init().
  ! States store pointers to phys_state(c) and pbuf which persist for the
  ! entire run.
  subroutine aerosol_instances_init_states(phys_state, pbuf2d)
    use radiative_aerosol_definitions, only: active_calls
    use modal_aerosol_state_mod, only: modal_aerosol_state
    use carma_aerosol_state_mod, only: carma_aerosol_state
    use bulk_aerosol_state_mod,  only: bulk_aerosol_state
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk
    use ppgrid,         only: begchunk, endchunk, pver
    use cam_abortutils, only: endrun

    type(physics_state),       intent(in), target :: phys_state(begchunk:endchunk)
    type(physics_buffer_desc),             pointer :: pbuf2d(:,:)

    integer :: iaermod, ilist, lchnk, istat
    type(physics_buffer_desc), pointer :: pbuf(:)
    character(len=*), parameter :: subname = 'aerosol_instances_init_states: '

    if (num_aero_models_ < 1) return

    allocate(aero_states_all(num_aero_models_, 0:N_DIAG, begchunk:endchunk), stat=istat)
    if (istat /= 0) then
       call endrun(subname//'allocation error: aero_states_all')
    end if

    do ilist = 0, N_DIAG
       if (.not. active_calls(ilist)) cycle

       do lchnk = begchunk, endchunk
          pbuf => pbuf_get_chunk(pbuf2d, lchnk)

          iaermod = 0
          if (modal_active_) then
             iaermod = iaermod + 1
             if (associated(aero_props_all(iaermod, ilist)%obj)) then
                aero_states_all(iaermod, ilist, lchnk)%obj => &
                     modal_aerosol_state(aero_props_all(iaermod, ilist)%obj, &
                     phys_state(lchnk)%ncol, pver, phys_state(lchnk), pbuf, ilist)
             end if
          end if
          if (carma_active_) then
             iaermod = iaermod + 1
             if (associated(aero_props_all(iaermod, ilist)%obj)) then
                aero_states_all(iaermod, ilist, lchnk)%obj => &
                     carma_aerosol_state(aero_props_all(iaermod, ilist)%obj, &
                     phys_state(lchnk)%ncol, pver, phys_state(lchnk), pbuf, ilist)
             end if
          end if
          if (bulk_active_) then
             iaermod = iaermod + 1
             if (associated(aero_props_all(iaermod, ilist)%obj)) then
                aero_states_all(iaermod, ilist, lchnk)%obj => &
                     bulk_aerosol_state(aero_props_all(iaermod, ilist)%obj, &
                     phys_state(lchnk)%ncol, pver, phys_state(lchnk), pbuf, ilist)
             end if
          end if
       end do
    end do

  end subroutine aerosol_instances_init_states

  ! Return a pointer to the persistent aerosol_state object for the given
  ! aerosol model index, radiation list, and chunk.
  function aerosol_instances_get_state(iaermod, list_idx, lchnk) result(astate)
    integer, intent(in) :: iaermod   ! aerosol model index (1..num_aero_models)
    integer, intent(in) :: list_idx  ! radiation list index (0=climate, 1..N_DIAG)
    integer, intent(in) :: lchnk     ! chunk index (begchunk..endchunk)
    class(aerosol_state), pointer :: astate

    astate => aero_states_all(iaermod, list_idx, lchnk)%obj

  end function aerosol_instances_get_state

  ! Create aerosol state objects for all active aerosol models.
  !
  ! This per-call factory is still needed for cases where the state is bound
  ! to a local copy of physics_state (e.g., microp_aero_run uses state1).
  !
  !REMOVECAM: no longer need this factory pattern once CAM is retired as cases
  ! where physics/chemistry uses state1 would be split off into separate physics
  ! schemes with tendency updaters in-between.
  subroutine aerosol_instances_create_states(list_idx, state, pbuf, aero_states, nstates)
    use modal_aerosol_state_mod, only: modal_aerosol_state
    use carma_aerosol_state_mod, only: carma_aerosol_state
    use bulk_aerosol_state_mod,  only: bulk_aerosol_state
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use ppgrid,         only: pver
    use cam_abortutils, only: endrun

    integer,                   intent(in)               :: list_idx
    type(physics_state),       intent(in),  target      :: state
    type(physics_buffer_desc),              pointer     :: pbuf(:)
    type(aero_state_entry_t),  intent(out), allocatable :: aero_states(:)    ! aerosol state objects
    integer,                   intent(out)              :: nstates           ! number of aerosol states created

    integer :: iaermod, istat
    character(len=*), parameter :: subname = 'aerosol_instances_create_states: '

    nstates = num_aero_models_
    if (nstates < 1) return

    allocate(aero_states(nstates), stat=istat)
    if (istat /= 0) then
       call endrun(subname//'allocation error: aero_states')
    end if

    iaermod = 0
    if (modal_active_) then
       iaermod = iaermod + 1
       aero_states(iaermod)%obj => modal_aerosol_state(aero_props_all(iaermod, list_idx)%obj, &
            state%ncol, pver, state, pbuf, list_idx)
    end if
    if (carma_active_) then
       iaermod = iaermod + 1
       aero_states(iaermod)%obj => carma_aerosol_state(aero_props_all(iaermod, list_idx)%obj, &
            state%ncol, pver, state, pbuf, list_idx)
    end if
    if (bulk_active_) then
       iaermod = iaermod + 1
       aero_states(iaermod)%obj => bulk_aerosol_state(aero_props_all(iaermod, list_idx)%obj, &
            state%ncol, pver, state, pbuf, list_idx)
    end if

  end subroutine aerosol_instances_create_states
  !REMOVECAM_END

  subroutine aerosol_instances_destroy_states(aero_states)
    type(aero_state_entry_t), allocatable, intent(inout) :: aero_states(:)
    integer :: i

    if (.not. allocated(aero_states)) return

    do i = 1, size(aero_states)
       if (associated(aero_states(i)%obj)) then
          deallocate(aero_states(i)%obj)
          nullify(aero_states(i)%obj)
       end if
    end do

    deallocate(aero_states)

  end subroutine aerosol_instances_destroy_states

end module aerosol_instances_mod
