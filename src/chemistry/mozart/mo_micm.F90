!----------------------------------------------------------------------
! Wrapper for the MICM chemistry solver (MUSICA library)
!
! Replaces the mechanism-generated implicit/explicit solvers (imp_sol,
! exp_sol) at their call site in mo_gas_phase_chemdr. Every compiled
! reaction is declared USER_DEFINED (or EMISSION) in the MICM mechanism
! configuration and its fully-evaluated rate constant is injected each
! timestep from CAM's reaction_rates array, so setrxt/usrrxt/photolysis/
! adjrxt are reused unchanged and MICM performs only the implicit
! integration. Rate injection converts CAM's post-adjrxt vmr-based rate
! constants to MICM's mol m-3 convention with k_micm = k_cam * c_air**(1-n),
! where n is the reaction's number of solution-species reactants and c_air
! the molar air density: invariant and third-body (M) concentrations are
! already folded into k_cam by adjrxt, so only the n solution reactants
! change units between the two conventions.
!
! The mechanism configuration and the reaction map companion file
! (micm_config_path, micm_rxt_map_path namelist variables) are derived
! artifacts of the compiled mechanism; micm_init cross-validates species
! and reaction counts and label lookups against chem_mods/mo_sim_dat and
! aborts on any mismatch.
!----------------------------------------------------------------------
module mo_micm
   use shr_kind_mod,   only : r8 => shr_kind_r8, cl => shr_kind_cl
   use spmd_utils,     only : is_main_task => masterproc
   use spmd_utils,     only : main_task_id => masterprocid
   use spmd_utils,     only : mpicom, mpi_character, mpi_logical, mpi_success
   use cam_abortutils, only : endrun
   use cam_logfile,    only : iulog
   use ppgrid,         only : pver, pcols
   use chem_mods,      only : gas_pcnst, rxntot, extcnt
#ifdef MICM
   use iso_fortran_env, only : real64
   use ieee_arithmetic, only : ieee_is_finite
   use phys_grid,       only : get_rlat_p, get_rlon_p
   use mo_tracname,     only : solsym
   use musica_micm,  only : micm_t, solver_stats_t
   use musica_micm,  only : Rosenbrock, RosenbrockStandardOrder, &
                            BackwardEuler, BackwardEulerStandardOrder
   use musica_state, only : state_t
   use musica_util,  only : error_t, string_t
#endif
   implicit none

   private

   public :: micm_readnl
   public :: micm_init
   public :: micm_solve
   public :: micm_final
   public :: micm_active

   ! namelist options
   logical, protected :: micm_active = .false.
   character(len=cl) :: micm_config_path  = 'NONE' ! absolute path to MICM mechanism configuration
   character(len=cl) :: micm_rxt_map_path = 'NONE' ! absolute path to reaction map companion file
   character(len=32) :: micm_solver_type = 'rosenbrock'
   logical :: micm_abort_on_nonconvergence = .true.

#ifdef MICM
   type(micm_t), pointer :: micm => null()

   ! Per-thread MICM solver states, indexed 1 + omp_get_thread_num().
   ! gas_phase_chemdr runs once per chunk inside the threaded physics
   ! loop, and MICM solvers are thread-safe only when each thread solves
   ! on its own state object; all states are created serially at init
   ! (same pattern as mo_tuvx's per-thread cores).
   type :: state_ptr
      type(state_t), pointer :: state_ => null()
   end type state_ptr
   type(state_ptr), allocatable :: states(:)

   integer :: state_size = 0 ! grid cells per solver state

   ! Index maps from the CAM mechanism (chem_mods/mo_sim_dat) ordering into
   ! the MICM state; built at init, read-only afterwards (shared across
   ! threads). MICM orderings are name-keyed with no guaranteed sequence, so
   ! indices must always come through these maps.
   !
   ! Reactions map through injection entries read from the reaction map
   ! companion file: usually one entry per reaction, but a reaction with no
   ! solution-species reactants maps to one EMISSION slot per solution
   ! product (each with its stoichiometric yield). No-op entries carry no
   ! rate parameter: NONE marks a reaction invisible to the solved system
   ! (no solution reactants or products), NATIVE one whose rate law is
   ! evaluated natively by MICM from the mechanism configuration
   ! (generator --native mode) rather than injected.
   integer, allocatable :: map_spc(:)  ! (gas_pcnst) MICM species variable index
   integer :: n_entries = 0            ! number of injection entries
   integer,  allocatable :: ent_rxt(:)   ! (n_entries) CAM reaction index
   integer,  allocatable :: ent_n(:)     ! (n_entries) number of solution-species reactants
   integer,  allocatable :: ent_param(:) ! (n_entries) MICM rate parameter index (-1 = no-op)
   real(r8), allocatable :: ent_yield(:) ! (n_entries) product yield factor
   integer, allocatable :: map_het(:)  ! (gas_pcnst) rate parameter index of LOSS.<species>
   integer, allocatable :: map_ext(:)  ! (extcnt)    rate parameter index of EMIS.ext_<species>

   real(r8) :: molec_cm3_to_mol_m3 = 0._r8 ! molecules cm-3 -> mol m-3

   ! Negative concentrations are clipped to zero on input and output, as imp_sol does inside its Newton iteration.
   ! Rosenbrock solvers do not preserve positivity, and a negative reactant turns a loss term into growth.
   ! Clips deeper than neg_log_threshold (the generated configs' absolute tolerance) are logged,
   ! at most neg_log_budget lines per task so routine undershoot cannot flood the log.
   real(r8), parameter :: neg_log_threshold = 1.e-12_r8 ! mol m-3
   integer :: neg_log_budget = 50

   real(r8), parameter :: rad2deg = 180._r8 / 3.14159265358979323846_r8
#endif

!================================================================================================
contains
!================================================================================================

   !-----------------------------------------------------------------------
   ! read namelist options
   !-----------------------------------------------------------------------
   subroutine micm_readnl(nlfile)

      use namelist_utils, only : find_group_name

      character(len=*), intent(in) :: nlfile ! filepath for file containing namelist input

      integer                     :: unitn, ierr
      character(len=*), parameter :: subname = 'micm_readnl'

      namelist /micm_opts/ micm_active, micm_config_path, micm_rxt_map_path, &
                           micm_solver_type, micm_abort_on_nonconvergence

      if (is_main_task) then
         open( newunit=unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'micm_opts', status=ierr)
         if (ierr == 0) then
            read(unitn, micm_opts, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
      end if

      call mpi_bcast(micm_active,      1,                      mpi_logical,   main_task_id, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//': mpi_bcast error : micm_active')
      call mpi_bcast(micm_config_path, len(micm_config_path),  mpi_character, main_task_id, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//': mpi_bcast error : micm_config_path')
      call mpi_bcast(micm_rxt_map_path, len(micm_rxt_map_path), mpi_character, main_task_id, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//': mpi_bcast error : micm_rxt_map_path')
      call mpi_bcast(micm_solver_type, len(micm_solver_type),  mpi_character, main_task_id, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//': mpi_bcast error : micm_solver_type')
      call mpi_bcast(micm_abort_on_nonconvergence, 1,          mpi_logical,   main_task_id, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//': mpi_bcast error : micm_abort_on_nonconvergence')

#ifdef MICM
      if (micm_active .and. micm_config_path == 'NONE') then
         call endrun(subname // ' : must set micm_config_path when MICM is active')
      end if
      if (micm_active .and. micm_rxt_map_path == 'NONE') then
         call endrun(subname // ' : must set micm_rxt_map_path when MICM is active')
      end if

      if (is_main_task) then
         write(iulog,*) 'micm_readnl: micm_active = ', micm_active
         write(iulog,*) 'micm_readnl: micm_config_path = ', trim(micm_config_path)
         write(iulog,*) 'micm_readnl: micm_rxt_map_path = ', trim(micm_rxt_map_path)
         write(iulog,*) 'micm_readnl: micm_solver_type = ', trim(micm_solver_type)
         write(iulog,*) 'micm_readnl: micm_abort_on_nonconvergence = ', micm_abort_on_nonconvergence
      end if
#else
      if (micm_active .or. micm_config_path /= 'NONE') then
         call endrun(subname // ' : use -micm configure CAM option to build the MICM library')
      end if
#endif
   end subroutine micm_readnl

!================================================================================================

   !-----------------------------------------------------------------------
   ! Creates the MICM solver from the mechanism configuration and builds
   ! the index maps between the compiled CAM mechanism and the MICM state.
   ! Must be called after set_sim_dat has populated chem_mods.
   !-----------------------------------------------------------------------
   subroutine micm_init( )
#ifdef MICM
      use chem_mods,   only : extfrc_lst
      use physconst,   only : avogad ! molecules kmole-1

      type(state_t), pointer :: state
      type(error_t)          :: error
      integer                :: solver_id, max_cells, m, ie, unitn, ierr
      integer                :: nrxt_file, nspc_file, next_file, nparam_expected
      logical                :: covered(rxntot)
      character(len=64)      :: param_name
      character(len=256)     :: line
      character(len=*), parameter :: subname = 'micm_init'

      if (.not. micm_active) return

      select case (trim(micm_solver_type))
      case ('rosenbrock')
         solver_id = Rosenbrock
      case ('backward_euler')
         solver_id = BackwardEuler
      case ('rosenbrock_standard')
         solver_id = RosenbrockStandardOrder
      case ('backward_euler_standard')
         solver_id = BackwardEulerStandardOrder
      case default
         call endrun(subname//': invalid micm_solver_type: '//trim(micm_solver_type))
      end select

      micm => micm_t(trim(micm_config_path), solver_id, error)
      call check_micm_error(error, subname//': creating MICM solver')

      ! (molecules cm-3) * 1e6 cm3 m-3 / (avogad * 1e-3 molecules mol-1)
      molec_cm3_to_mol_m3 = 1.e9_r8 / avogad

      ! Vector-ordered solvers cap the grid cells per state at the compiled
      ! vector size; standard-order solvers are unbounded and get one
      ! chunk-sized state (single solve per chunk).
      max_cells = micm%get_maximum_number_of_grid_cells()
      state_size = min(max_cells, pcols*pver)

      ! temporary single-cell state to read the name->index orderings
      state => micm%get_state(1, error)
      call check_micm_error(error, subname//': creating MICM state')

      if (state%species_ordering%size() /= gas_pcnst) then
         write(iulog,*) subname, ': MICM configuration has ', &
            state%species_ordering%size(), ' species; compiled mechanism has ', gas_pcnst
         call endrun(subname//': MICM configuration does not match the compiled mechanism (species count)')
      end if

      allocate(map_spc(gas_pcnst))
      do m = 1, gas_pcnst
         map_spc(m) = state%species_ordering%index(trim(solsym(m)), error)
         call check_micm_error(error, subname//': species missing from MICM configuration: '//trim(solsym(m)))
      end do

      ! Reaction map companion file: injection entries giving, for each
      ! compiled reaction, its number of solution-species reactants (the
      ! unit conversion exponent; chem_mods' num_rnts cannot be used because
      ! it counts invariant reactants too), the MICM rate parameter name,
      ! and a product yield factor.
      open(newunit=unitn, file=trim(micm_rxt_map_path), status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
         call endrun(subname//': cannot open micm_rxt_map_path: '//trim(micm_rxt_map_path))
      end if
      call read_data_line(unitn, line, ierr)
      if (ierr /= 0) call endrun(subname//': error reading reaction map header')
      read(line,*,iostat=ierr) nrxt_file, nspc_file, next_file, n_entries
      if (ierr /= 0) call endrun(subname//': error parsing reaction map header')
      if (nrxt_file /= rxntot .or. nspc_file /= gas_pcnst .or. next_file /= extcnt) then
         write(iulog,*) subname, ': reaction map header (', nrxt_file, nspc_file, next_file, &
            ') does not match compiled mechanism (', rxntot, gas_pcnst, extcnt, ')'
         call endrun(subname//': reaction map does not match the compiled mechanism')
      end if
      allocate(ent_rxt(n_entries))
      allocate(ent_n(n_entries))
      allocate(ent_param(n_entries))
      allocate(ent_yield(n_entries))
      covered(:) = .false.
      nparam_expected = 0
      do ie = 1, n_entries
         call read_data_line(unitn, line, ierr)
         if (ierr /= 0) call endrun(subname//': error reading reaction map entry')
         read(line,*,iostat=ierr) ent_rxt(ie), ent_n(ie), param_name, ent_yield(ie)
         if (ierr /= 0 .or. ent_rxt(ie) < 1 .or. ent_rxt(ie) > rxntot) then
            write(iulog,*) subname, ': bad reaction map entry ', ie, ': ', trim(line)
            call endrun(subname//': error parsing reaction map entry')
         end if
         covered(ent_rxt(ie)) = .true.
         if (trim(param_name) == 'NONE' .or. trim(param_name) == 'NATIVE') then
            ent_param(ie) = -1
         else
            ent_param(ie) = state%rate_parameters_ordering%index(trim(param_name), error)
            call check_micm_error(error, subname//': reaction missing from MICM configuration: '//trim(param_name))
            nparam_expected = nparam_expected + 1
         end if
      end do
      close(unitn)
      if (.not. all(covered)) then
         call endrun(subname//': reaction map does not cover every compiled reaction')
      end if

      ! Heterogeneous (washout) loss and external forcing slots. The species
      ! subject to washout are selected at runtime (gas_wetdep_list), so the
      ! configuration carries a LOSS.<species> slot for every solution
      ! species and an EMIS.ext_<species> slot for every external forcing;
      ! unused slots are injected as zero.
      allocate(map_het(gas_pcnst))
      do m = 1, gas_pcnst
         map_het(m) = state%rate_parameters_ordering%index('LOSS.'//trim(solsym(m)), error)
         call check_micm_error(error, subname//': loss slot missing from MICM configuration: LOSS.'//trim(solsym(m)))
      end do
      allocate(map_ext(extcnt))
      do m = 1, extcnt
         map_ext(m) = state%rate_parameters_ordering%index('EMIS.ext_'//trim(extfrc_lst(m)), error)
         call check_micm_error(error, &
            subname//': external forcing slot missing from MICM configuration: EMIS.ext_'//trim(extfrc_lst(m)))
      end do

      if (state%rate_parameters_ordering%size() /= nparam_expected + gas_pcnst + extcnt) then
         write(iulog,*) subname, ': MICM configuration has ', state%rate_parameters_ordering%size(), &
            ' rate parameters; expected ', nparam_expected + gas_pcnst + extcnt
         call endrun(subname//': MICM configuration does not match the compiled mechanism (rate parameter count)')
      end if

      deallocate(state)

      allocate(states(max_threads()))
      do m = 1, size(states)
         states(m)%state_ => micm%get_state(state_size, error)
         call check_micm_error(error, subname//': creating per-thread MICM state')
      end do

      if (is_main_task) then
         write(iulog,*) subname, ': MICM solver created: ', trim(micm_solver_type), &
            ', grid cells per state = ', state_size, ', states = ', size(states)
      end if
#endif
   end subroutine micm_init

!================================================================================================

   !-----------------------------------------------------------------------
   ! Solves the chemistry ODE system for one chunk, replacing exp_sol +
   ! imp_sol. All gas_pcnst species are solved implicitly, including the
   ! (at most few) species the generated solvers treat with the explicit
   ! forward-Euler class: MICM has no class split, and folding them into
   ! the implicit solve is a deliberate, documented behavioral difference.
   !
   ! Inputs are the arrays exactly as they stand at the imp_sol call site
   ! in gas_phase_chemdr: reaction_rates post-adjrxt/phtadj (vmr-based
   ! effective rate constants), het_rates (s-1), extfrc (vmr s-1).
   !-----------------------------------------------------------------------
   subroutine micm_solve( ncol, lchnk, delt, xhnm, tfld, pmid, vmr, &
                          reaction_rates, het_rates, extfrc )

      integer,  intent(in)    :: ncol                    ! number of columns in chunk
      integer,  intent(in)    :: lchnk                   ! chunk index (diagnostics only)
      real(r8), intent(in)    :: delt                    ! time step (s)
      real(r8), intent(in)    :: xhnm(ncol,pver)         ! total air density (molecules cm-3)
      real(r8), intent(in)    :: tfld(ncol,pver)         ! temperature (K)
      real(r8), intent(in)    :: pmid(ncol,pver)         ! midpoint pressure (Pa)
      real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)          ! mixing ratios (mol mol-1)
      real(r8), intent(in)    :: reaction_rates(ncol,pver,max(1,rxntot))
      real(r8), intent(in)    :: het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rates (s-1)
      real(r8), intent(in)    :: extfrc(ncol,pver,max(1,extcnt))       ! external forcing (vmr s-1)

#ifdef MICM
      type(error_t)        :: error
      type(string_t)       :: solver_state
      type(solver_stats_t) :: stats
      type(state_t), pointer :: thread_state
      real(r8) :: cair(ncol,pver) ! molar air density (mol m-3)
      real(r8) :: conc            ! concentration (mol m-3)
      integer  :: ncells, nblocks, nvalid, offset, iblk, icell, cellg
      integer  :: i, k, m, ie, ibase
      ! Clipping tallies for this chunk: [1] input, [2] output; worst = most negative concentration.
      integer  :: nneg(2), worst_i(2), worst_k(2), worst_m(2)
      real(r8) :: worst_conc(2), worst_vmr0(2)
      integer  :: cell_str_c, var_str_c ! concentration strides
      integer  :: cell_str_p, var_str_p ! rate parameter strides
      character(len=*), parameter :: subname = 'micm_solve'

      thread_state => states(thread_id())%state_
      cell_str_c = thread_state%species_strides%grid_cell
      var_str_c  = thread_state%species_strides%variable
      cell_str_p = thread_state%rate_parameters_strides%grid_cell
      var_str_p  = thread_state%rate_parameters_strides%variable

      cair(:,:) = xhnm(:,:) * molec_cm3_to_mol_m3
      nneg(:) = 0
      worst_conc(:) = 0._r8

      ! Cells are flattened column-fastest, matching the storage order of
      ! the (ncol,pver) arrays, and solved in blocks of state_size cells.
      ncells  = ncol*pver
      nblocks = (ncells + state_size - 1)/state_size

      do iblk = 1, nblocks
         offset = (iblk-1)*state_size
         nvalid = min(state_size, ncells - offset)

         do icell = 1, nvalid
            cellg = offset + icell
            i = mod(cellg-1, ncol) + 1
            k = (cellg-1)/ncol + 1
            thread_state%conditions(icell)%temperature = tfld(i,k)
            thread_state%conditions(icell)%pressure    = pmid(i,k)
            thread_state%conditions(icell)%air_density = cair(i,k)
            ibase = 1 + (icell-1)*cell_str_c
            do m = 1, gas_pcnst
               conc = vmr(i,k,m) * cair(i,k)
               if (conc < 0._r8) then
                  call tally_neg(1, conc, vmr(i,k,m), i, k, m)
                  conc = 0._r8
               end if
               thread_state%concentrations(ibase + (map_spc(m)-1)*var_str_c) = conc
            end do
            ibase = 1 + (icell-1)*cell_str_p
            do ie = 1, n_entries
               if (ent_param(ie) < 0) cycle
               ! vmr-based rate constant -> mol m-3 based rate constant
               thread_state%rate_parameters(ibase + (ent_param(ie)-1)*var_str_p) = &
                  reaction_rates(i,k,ent_rxt(ie)) * cair(i,k)**(1 - ent_n(ie)) * ent_yield(ie)
            end do
            do m = 1, gas_pcnst
               thread_state%rate_parameters(ibase + (map_het(m)-1)*var_str_p) = &
                  het_rates(i,k,m)
            end do
            do m = 1, extcnt
               ! vmr s-1 -> mol m-3 s-1
               thread_state%rate_parameters(ibase + (map_ext(m)-1)*var_str_p) = &
                  extfrc(i,k,m) * cair(i,k)
            end do
         end do

         ! Pad unused lanes of the final block with an inert system (zero
         ! concentrations and rates under physical conditions), which the
         ! solver integrates trivially; results are discarded.
         do icell = nvalid+1, state_size
            thread_state%conditions(icell)%temperature = 250._r8
            thread_state%conditions(icell)%pressure    = 1.e4_r8
            thread_state%conditions(icell)%air_density = 1.e4_r8/(8.31446_r8*250._r8)
            ibase = 1 + (icell-1)*cell_str_c
            do m = 1, gas_pcnst
               thread_state%concentrations(ibase + (map_spc(m)-1)*var_str_c) = 0._r8
            end do
            ibase = 1 + (icell-1)*cell_str_p
            do m = 1, thread_state%number_of_rate_parameters
               thread_state%rate_parameters(ibase + (m-1)*var_str_p) = 0._r8
            end do
         end do

         call micm%solve(real(delt, real64), thread_state, solver_state, stats, error)
         call check_micm_error(error, subname//': MICM solve failed')

         if (solver_state%get_char_array() /= 'Converged') then
            write(iulog,*) subname, ': MICM solver did not converge: ', &
               solver_state%get_char_array(), ' chunk ', lchnk, ' block ', iblk
            if (micm_abort_on_nonconvergence) then
               call endrun(subname//': MICM solver did not converge'// &
                  ' (set micm_abort_on_nonconvergence=.false. to continue on non-convergence)')
            end if
         end if

         do icell = 1, nvalid
            cellg = offset + icell
            i = mod(cellg-1, ncol) + 1
            k = (cellg-1)/ncol + 1
            ibase = 1 + (icell-1)*cell_str_c
            do m = 1, gas_pcnst
               conc = thread_state%concentrations(ibase + (map_spc(m)-1)*var_str_c)
               ! MICM flags NaN/Inf only through its error norm, so check what is handed back to CAM.
               ! Any value above unit mixing ratio is unphysical for a solution species.
               if (.not. ieee_is_finite(conc) .or. conc > cair(i,k)) then
                  write(iulog,*) subname, ': unphysical MICM output ', trim(solsym(m)), ' = ', conc, &
                     ' mol m-3 (vmr before solve ', vmr(i,k,m), ') at lat ', get_rlat_p(lchnk,i)*rad2deg, &
                     ' lon ', get_rlon_p(lchnk,i)*rad2deg, ' k ', k, ' T ', tfld(i,k), ' p ', pmid(i,k), &
                     ' solver state ', solver_state%get_char_array()
                  call endrun(subname//': unphysical MICM output concentration for '//trim(solsym(m)))
               end if
               if (conc < 0._r8) then
                  call tally_neg(2, conc, vmr(i,k,m), i, k, m)
                  conc = 0._r8
               end if
               vmr(i,k,m) = conc / cair(i,k)
            end do
         end do
      end do

      call log_neg()
#else
      call endrun('micm_solve: use -micm configure CAM option to build the MICM library')
#endif

#ifdef MICM
   contains

      ! Records a clipped negative concentration in tally n (1 = input, 2 = output).
      subroutine tally_neg(n, c, vmr_before, ic, kc, mc)
         integer,  intent(in) :: n, ic, kc, mc
         real(r8), intent(in) :: c, vmr_before

         nneg(n) = nneg(n) + 1
         if (c < worst_conc(n)) then
            worst_conc(n) = c
            worst_vmr0(n) = vmr_before
            worst_i(n) = ic
            worst_k(n) = kc
            worst_m(n) = mc
         end if
      end subroutine tally_neg

      ! Logs one line per tally for this chunk if its worst clip exceeds neg_log_threshold.
      subroutine log_neg()
         character(len=6), parameter :: label(2) = (/ 'input ', 'output' /)
         integer :: n

         do n = 1, 2
            if (worst_conc(n) >= -neg_log_threshold) cycle
            !$omp critical (micm_neg_log)
            if (neg_log_budget > 0) then
               neg_log_budget = neg_log_budget - 1
               write(iulog,*) subname, ': clipped ', nneg(n), ' negative ', trim(label(n)), &
                  ' concentrations in chunk ', lchnk, '; worst ', trim(solsym(worst_m(n))), ' = ', &
                  worst_conc(n), ' mol m-3 (vmr before solve ', worst_vmr0(n), ') at lat ', &
                  get_rlat_p(lchnk,worst_i(n))*rad2deg, ' lon ', get_rlon_p(lchnk,worst_i(n))*rad2deg, &
                  ' k ', worst_k(n), ' T ', tfld(worst_i(n),worst_k(n)), ' p ', pmid(worst_i(n),worst_k(n))
               if (neg_log_budget == 0) then
                  write(iulog,*) subname, ': further negative-clipping messages suppressed on this task'
               end if
            end if
            !$omp end critical (micm_neg_log)
         end do
      end subroutine log_neg
#endif
   end subroutine micm_solve

!================================================================================================

   !-----------------------------------------------------------------------
   ! Deallocates the MICM solver and states
   !-----------------------------------------------------------------------
   subroutine micm_final( )
#ifdef MICM
      integer :: m

      if (.not. micm_active) return

      if (allocated(states)) then
         do m = 1, size(states)
            if (associated(states(m)%state_)) then
               deallocate(states(m)%state_)
               nullify(states(m)%state_)
            end if
         end do
         deallocate(states)
      end if

      if (associated(micm)) then
         deallocate(micm)
         nullify(micm)
      end if
      if (allocated(map_spc))   deallocate(map_spc)
      if (allocated(ent_rxt))   deallocate(ent_rxt)
      if (allocated(ent_n))     deallocate(ent_n)
      if (allocated(ent_param)) deallocate(ent_param)
      if (allocated(ent_yield)) deallocate(ent_yield)
      if (allocated(map_het))   deallocate(map_het)
      if (allocated(map_ext))   deallocate(map_ext)
#endif
   end subroutine micm_final

#ifdef MICM
!================================================================================================

   !-----------------------------------------------------------------------
   ! Aborts with the MICM error message if an operation failed
   !-----------------------------------------------------------------------
   subroutine check_micm_error( error, message )

      type(error_t),    intent(in) :: error
      character(len=*), intent(in) :: message

      if (.not. error%is_success()) then
         write(iulog,*) message, ': ', error%message()
         call endrun(message)
      end if

   end subroutine check_micm_error

   !-----------------------------------------------------------------------
   ! Index for the current OpenMP thread (1 <= id <= max_threads())
   !-----------------------------------------------------------------------
   integer function thread_id( )
#ifdef _OPENMP
      use omp_lib, only : omp_get_thread_num
      thread_id = 1 + omp_get_thread_num( )
#else
      thread_id = 1
#endif
   end function thread_id

   !-----------------------------------------------------------------------
   ! Maximum number of OpenMP threads
   !-----------------------------------------------------------------------
   integer function max_threads( )
#ifdef _OPENMP
      use omp_lib, only : omp_get_max_threads
      max_threads = omp_get_max_threads( )
#else
      max_threads = 1
#endif
   end function max_threads

   !-----------------------------------------------------------------------
   ! Reads the next non-comment, non-blank line from the reaction map file
   !-----------------------------------------------------------------------
   subroutine read_data_line( unitn, line, ierr )

      integer,          intent(in)  :: unitn
      character(len=*), intent(out) :: line
      integer,          intent(out) :: ierr

      do
         read(unitn,'(a)',iostat=ierr) line
         if (ierr /= 0) return
         line = adjustl(line)
         if (len_trim(line) > 0 .and. line(1:1) /= '#') return
      end do

   end subroutine read_data_line
#endif

end module mo_micm
