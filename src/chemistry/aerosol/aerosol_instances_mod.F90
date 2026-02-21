module aerosol_instances_mod
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_state_mod, only: aerosol_state
  use rad_constituents, only: N_DIAG, rad_cnst_get_info, rad_cnst_get_call_list

  implicit none
  private

  public :: aerosol_instances_init
  public :: aerosol_instances_get_props
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
     class(aerosol_state), pointer :: obj => null()
  end type aero_state_entry_t

  type(aero_props_entry_t), allocatable, target :: aero_props_all(:,:) ! (iaermod, 0:N_DIAG)
  integer :: num_aero_models_ = 0
  logical :: modal_active_ = .false.
  logical :: carma_active_ = .false.
  logical :: bulk_active_  = .false.

contains

  !------------------------------------------------------------------------------
  subroutine aerosol_instances_init()
    use modal_aerosol_properties_mod, only: modal_aerosol_properties
    use carma_aerosol_properties_mod, only: carma_aerosol_properties
    use bulk_aerosol_properties_mod,  only: bulk_aerosol_properties
    use cam_abortutils, only: endrun

    integer :: nmodes, nbins, nbulk_aerosols
    integer :: iaermod, ilist, istat
    logical :: call_list(0:N_DIAG)

    character(len=*), parameter :: prefix = 'aerosol_instances_init: '

    num_aero_models_ = 0

    call rad_cnst_get_info(0, nmodes=nmodes, nbins=nbins, naero=nbulk_aerosols)
    modal_active_ = nmodes > 0
    carma_active_ = nbins > 0
    bulk_active_  = nbulk_aerosols > 0

    if (modal_active_) num_aero_models_ = num_aero_models_ + 1
    if (carma_active_) num_aero_models_ = num_aero_models_ + 1
    if (bulk_active_)  num_aero_models_ = num_aero_models_ + 1

    if (num_aero_models_ < 1) return

    allocate(aero_props_all(num_aero_models_, 0:N_DIAG), stat=istat)
    if (istat /= 0) then
       call endrun(prefix//'allocation error: aero_props_all')
    end if

    call rad_cnst_get_call_list(call_list)

    do ilist = 0, N_DIAG
       if (.not. call_list(ilist)) cycle
       call rad_cnst_get_info(ilist, nmodes=nmodes, nbins=nbins, naero=nbulk_aerosols)

       iaermod = 0
       if (modal_active_) then
          iaermod = iaermod + 1
          if (nmodes > 0) then
             aero_props_all(iaermod, ilist)%obj => modal_aerosol_properties(ilist)
          end if
       end if
       if (carma_active_) then
          iaermod = iaermod + 1
          if (nbins > 0) then
             aero_props_all(iaermod, ilist)%obj => carma_aerosol_properties(ilist)
          end if
       end if
       if (bulk_active_) then
          iaermod = iaermod + 1
          if (nbulk_aerosols > 0) then
             aero_props_all(iaermod, ilist)%obj => bulk_aerosol_properties(ilist)
          end if
       end if
    end do

  end subroutine aerosol_instances_init

  !------------------------------------------------------------------------------
  function aerosol_instances_get_props(iaermod, list_idx) result(props)
    integer, intent(in) :: iaermod
    integer, intent(in) :: list_idx
    class(aerosol_properties), pointer :: props

    props => aero_props_all(iaermod, list_idx)%obj

  end function aerosol_instances_get_props

  !------------------------------------------------------------------------------
  pure integer function aerosol_instances_get_num_models()
    aerosol_instances_get_num_models = num_aero_models_
  end function aerosol_instances_get_num_models

  !------------------------------------------------------------------------------
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

  !------------------------------------------------------------------------------
  subroutine aerosol_instances_final()
    integer :: iaermod, ilist

    if (.not. allocated(aero_props_all)) return

    do ilist = 0, N_DIAG
       do iaermod = 1, num_aero_models_
          if (associated(aero_props_all(iaermod, ilist)%obj)) then
             deallocate(aero_props_all(iaermod, ilist)%obj)
             nullify(aero_props_all(iaermod, ilist)%obj)
          end if
       end do
    end do

    deallocate(aero_props_all)
    num_aero_models_ = 0

  end subroutine aerosol_instances_final

  !------------------------------------------------------------------------------
  !> Creates aerosol state objects for all active aerosol models.
  !! Unlike aerosol properties (pre-allocated as shared module data), state objects
  !! are created per-call because they hold mutable per-chunk state/pbuf pointers.
  !! Since aerosol_optics_cam_sw/lw runs inside an OMP chunk-parallel loop,
  !! shared mutable state objects would have race conditions.
  subroutine aerosol_instances_create_states(list_idx, state, pbuf, aero_states, nstates)
    use modal_aerosol_state_mod, only: modal_aerosol_state
    use carma_aerosol_state_mod, only: carma_aerosol_state
    use bulk_aerosol_state_mod,  only: bulk_aerosol_state
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use cam_abortutils, only: endrun

    integer, intent(in) :: list_idx
    type(physics_state), intent(in), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(aero_state_entry_t), allocatable, intent(out) :: aero_states(:)
    integer, intent(out) :: nstates

    integer :: iaermod, istat
    character(len=*), parameter :: prefix = 'aerosol_instances_create_states: '

    nstates = num_aero_models_
    if (nstates < 1) return

    allocate(aero_states(nstates), stat=istat)
    if (istat /= 0) then
       call endrun(prefix//'allocation error: aero_states')
    end if

    iaermod = 0
    if (modal_active_) then
       iaermod = iaermod + 1
       aero_states(iaermod)%obj => modal_aerosol_state(state, pbuf, list_idx)
    end if
    if (carma_active_) then
       iaermod = iaermod + 1
       aero_states(iaermod)%obj => carma_aerosol_state(state, pbuf, list_idx)
    end if
    if (bulk_active_) then
       iaermod = iaermod + 1
       aero_states(iaermod)%obj => bulk_aerosol_state(state, pbuf, list_idx)
    end if

  end subroutine aerosol_instances_create_states

  !------------------------------------------------------------------------------
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
