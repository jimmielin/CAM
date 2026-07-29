! Strategy A vs strategy B equivalence test for trop_strat_mam5_t1s1.
!
! Uses the REAL generated t1s1 rate code (chem_mods/mo_sim_dat/mo_setrxt/
! mo_adjrxt/mo_phtadj) to produce CAM-truth post-adjrxt rate constants at
! per-cell-distinct conditions, with deterministic synthetic values for the
! slots CAM computes elsewhere (photolysis j's, usrrxt rates, het_rates,
! extfrc). The same initial state is then solved through mo_micm twice:
! strategy A (micm/: every rate injected) and strategy B (micm_native/:
! 368 Arrhenius/Troe/photolysis rate laws evaluated natively by MICM from
! coefficients copied out of chem_mech.in, remainder injected). Any
! mistranslation of coefficients, Troe exponent signs, or v0 unit
! conventions shows up as A-B divergence far above solver noise.
program ab_harness
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pcols
  use chem_mods,    only: gas_pcnst, rxntot, extcnt, nfs, phtcnt
  use mo_sim_dat,   only: set_sim_dat
  use mo_setrxt,    only: setrxt
  use mo_adjrxt,    only: adjrxt
  use mo_phtadj,    only: phtadj
  use mo_micm,      only: micm_readnl, micm_init, micm_solve, micm_final
  implicit none

  real(r8) :: delt = 600._r8
  real(r8), parameter :: kboltz_cgs = 1.380649e-16_r8 ! erg/K
  real(r8), parameter :: tol = 1.e-3_r8   ! relative; translation errors are orders larger
  real(r8), parameter :: floor = 1.e-16_r8

  integer  :: ncol, i, k, m, r, cell, nbad
  character(len=256) :: line
  real(r8) :: ta, tb
  real(r8) :: temp(pcols,pver), pres(pcols,pver), xhnm(pcols,pver)
  real(r8) :: inv(pcols,pver,nfs)
  real(r8) :: rate(pcols,pver,rxntot)
  real(r8) :: vmr0(pcols,pver,gas_pcnst)
  real(r8) :: vmr_a(pcols,pver,gas_pcnst), vmr_b(pcols,pver,gas_pcnst)
  real(r8) :: het(pcols,pver,gas_pcnst), ext(pcols,pver,extcnt)
  real(r8) :: rel, worst
  character(len=256) :: nml_a, nml_b, map_path

  ! args: <strategy-A nml> <strategy-B nml> <strategy-A rxt_map.txt> [delt]
  call get_command_argument(1, nml_a)
  call get_command_argument(2, nml_b)
  call get_command_argument(3, map_path)
  call get_command_argument(4, line)
  if (len_trim(line) > 0) read(line,*) delt

  ncol = pcols
  call set_sim_dat()

  ! conditions spanning troposphere to stratosphere
  do k = 1, pver
    do i = 1, ncol
      cell = (k-1)*ncol + i
      temp(i,k) = 200._r8 + 100._r8*frac(cell)
      pres(i,k) = 1.e3_r8 * 10._r8**(2.0_r8*frac(cell+7))   ! 1e3..1e5 Pa
      xhnm(i,k) = pres(i,k)*10._r8 / (kboltz_cgs*temp(i,k)) ! molecules cm-3
      inv(i,k,1) = xhnm(i,k)              ! M
      inv(i,k,2) = 0.209_r8*xhnm(i,k)     ! O2
      inv(i,k,3) = 0.781_r8*xhnm(i,k)     ! N2
    end do
  end do

  ! CAM-truth rates for the setrxt-computed (native-candidate) slots:
  ! zero-fill, let setrxt overwrite its slots, and record which slots it
  ! left untouched (photolysis + usrrxt). adjrxt/phtadj multiplications
  ! leave the zeroed slots zero, so afterwards the untouched slots are
  ! given sane synthetic POST-adjrxt values scaled by the reaction's
  ! solution-reactant count n (from the strategy-A reaction map): the
  ! injected value has units vmr**(1-n)/s, so an effective per-species
  ! rate of ~1e-5 s-1 at vmr ~1e-9 needs ~1e-5 * (1e-9)**(1-n).
  rate(:,:,:) = 0._r8
  call setrxt( rate, temp, inv(1,1,1), ncol )
  call adjrxt( rate, inv, inv(:,:,1), ncol, pver )
  call phtadj( rate, inv, inv(:,:,1), ncol, pver )
  call fill_untouched_slots()

  do m = 1, gas_pcnst
    do k = 1, pver
      do i = 1, ncol
        cell = (k-1)*ncol + i
        vmr0(i,k,m) = 1.e-9_r8 * (0.05_r8 + frac(m*31+cell))
        het(i,k,m)  = 1.e-6_r8 * frac(m*13+cell)
      end do
    end do
  end do
  do m = 1, extcnt
    do k = 1, pver
      do i = 1, ncol
        ext(i,k,m) = 1.e-13_r8 * frac(m*41+(k-1)*ncol+i)
      end do
    end do
  end do

  call micm_readnl(trim(nml_a))
  call micm_init()
  vmr_a = vmr0
  call micm_solve( ncol, 1, delt, xhnm, temp, pres, vmr_a, rate, het, ext )
  call micm_final()

  call micm_readnl(trim(nml_b))
  call micm_init()
  vmr_b = vmr0
  call micm_solve( ncol, 1, delt, xhnm, temp, pres, vmr_b, rate, het, ext )
  call micm_final()

  ! Primary criterion: per-cell TENDENCY agreement, which stays sensitive
  ! to rate-representation errors at arbitrarily small delt (a small-delt
  ! run isolates the rates from stiff-trajectory sensitivity; delt=1e-4 s
  ! is the canonical falsification run).
  nbad = 0
  do m = 1, gas_pcnst
    worst = 0._r8
    do k = 1, pver
      do i = 1, ncol
        ta = (vmr_a(i,k,m)-vmr0(i,k,m))/delt
        tb = (vmr_b(i,k,m)-vmr0(i,k,m))/delt
        rel = abs(ta-tb)/max(abs(ta), abs(tb), 1.e-25_r8)
        if (rel > worst) worst = rel
        if (rel > tol) nbad = nbad + 1
      end do
    end do
    if (worst > tol) print '(a,i4,a,es9.2)', &
      'species ', m, ' worst tendency rel diff ', worst
  end do
  worst = 0._r8
  do m = 1, gas_pcnst
   do k = 1, pver
    do i = 1, ncol
      rel = abs(vmr_a(i,k,m) - vmr_b(i,k,m)) / max(abs(vmr_a(i,k,m)), floor)
      worst = max(worst, rel)
    end do
   end do
  end do
  print '(a,es9.2)', 'max A-vs-B relative difference: ', worst
  if (nbad == 0) then
    print '(a)', 'PASS'
  else
    print '(a,i0,a)', 'FAIL: ', nbad, ' cells diverged'
    stop 1
  end if

contains

  real(r8) function frac(n)
    integer, intent(in) :: n
    frac = 0.5_r8*(1._r8 + sin(real(n, r8)))
  end function frac

  subroutine fill_untouched_slots()
    ! read n per reaction from the strategy-A map and fill the slots
    ! setrxt left at zero with plausible post-adjrxt values
    integer :: unitn, ierr, idx, n, nhdr(4), ecount
    real(r8) :: yld, val
    character(len=64)  :: nm
    character(len=256) :: line
    open(newunit=unitn, file=trim(map_path), status='old', action='read')
    call next_line(unitn, line)
    read(line,*) nhdr
    do ecount = 1, nhdr(4)
      call next_line(unitn, line)
      read(line,*,iostat=ierr) idx, n, nm, yld
      if (ierr /= 0) stop 'bad map line'
      if (any(rate(:,:,idx) /= 0._r8)) cycle   ! setrxt-computed slot
      do k = 1, pver
        do i = 1, ncol
          cell = (k-1)*ncol + i
          val = 1.e-5_r8 * (0.2_r8 + 0.8_r8*frac(idx*17+cell)) * (1.e-9_r8)**(1-n)
          rate(i,k,idx) = val
        end do
      end do
    end do
    close(unitn)
  end subroutine fill_untouched_slots

  subroutine next_line(unitn, line)
    integer, intent(in) :: unitn
    character(len=*), intent(out) :: line
    integer :: ierr
    do
      read(unitn,'(a)',iostat=ierr) line
      if (ierr /= 0) stop 'map read error'
      line = adjustl(line)
      if (len_trim(line) > 0 .and. line(1:1) /= '#') return
    end do
  end subroutine next_line

end program ab_harness
