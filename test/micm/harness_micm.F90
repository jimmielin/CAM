! Known-answer test of the real mo_micm.F90 shim against the generated
! trop_mam4 MICM configuration.
!
! With rates held fixed over the step (exactly what strategy-A injection
! does), the trop_mam4 mechanism is linear in the solution species, so the
! expected concentrations after one solve have closed forms. Every cell in
! the chunk gets distinct conditions and rates (scaled by a per-cell
! factor) so stride/indexing bugs across the vector-state lanes are caught.
! A second solve with a smaller chunk exercises the padded-lanes path.
!
! Species checked against exact solutions: DMS, H2O2, SO2, SOAE, SOAG,
! soa_a1, so4_a1 (pure external forcing), plus a total-sulfur linear
! invariant that includes H2SO4.
program harness_micm
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pcols
  use chem_mods,    only: gas_pcnst, rxntot, extcnt
  use mo_micm,      only: micm_readnl, micm_init, micm_solve, micm_final, micm_active
  implicit none

  ! trop_mam4 solsym indices
  integer, parameter :: iDMS = 3, iH2O2 = 7, iH2SO4 = 8, iSO2 = 18, &
                        iso4_a1 = 19, isoa_a1 = 22, iSOAE = 24, iSOAG = 25
  ! extfrc_lst indices
  integer, parameter :: xSO2 = 1, xso4_a1 = 2

  real(r8), parameter :: delt = 1800._r8
  real(r8), parameter :: kboltz_cgs = 1.380649e-16_r8 ! erg/K
  real(r8), parameter :: tol = 2.e-3_r8               ! relative; solver rtol-limited

  character(len=256) :: nlfile
  integer  :: ncol, i, k, m, nfail
  real(r8) :: fac
  real(r8), allocatable :: xhnm(:,:), tfld(:,:), pmid(:,:)
  real(r8), allocatable :: vmr(:,:,:), rxt(:,:,:), het(:,:,:), ext(:,:,:)

  call get_command_argument(1, nlfile)
  call micm_readnl(trim(nlfile))
  if (.not. micm_active) stop 'micm_active is false'
  call micm_init()

  nfail = 0
  call run_chunk(pcols)  ! full final block
  call run_chunk(2)      ! ncells < state_size: padded-lanes path
  call micm_final()

  if (nfail == 0) then
    print '(a)', 'PASS'
  else
    print '(a,i0,a)', 'FAIL: ', nfail, ' checks failed'
    stop 1
  end if

contains

  subroutine run_chunk(ncol_in)
    integer, intent(in) :: ncol_in
    integer  :: cell
    real(r8) :: j1, j2, k4, e5, k6, k7, k8, k9, k10
    real(r8) :: hH2O2, hSO2, eSO2, eso4
    real(r8) :: ld, lh, ls, a, y, s0, s1

    ncol = ncol_in
    allocate(xhnm(ncol,pver), tfld(ncol,pver), pmid(ncol,pver))
    allocate(vmr(ncol,pver,gas_pcnst), rxt(ncol,pver,rxntot))
    allocate(het(ncol,pver,gas_pcnst), ext(ncol,pver,extcnt))

    do k = 1, pver
      do i = 1, ncol
        cell = (k-1)*ncol + i
        fac = 1._r8 + 0.02_r8*(cell-1)
        tfld(i,k) = 270._r8 + 5._r8*cell
        pmid(i,k) = 6.e4_r8 + 5.e3_r8*cell
        ! molecules cm-3, consistent with T,p
        xhnm(i,k) = pmid(i,k)*10._r8 / (kboltz_cgs*tfld(i,k))

        vmr(i,k,:)       = 1.e-11_r8
        vmr(i,k,iDMS)    = 1.e-10_r8 * fac
        vmr(i,k,iH2O2)   = 1.e-9_r8  * fac
        vmr(i,k,iSO2)    = 5.e-10_r8 * fac
        vmr(i,k,iH2SO4)  = 1.e-13_r8
        vmr(i,k,iSOAE)   = 2.e-10_r8 * fac
        vmr(i,k,iSOAG)   = 1.e-11_r8
        vmr(i,k,isoa_a1) = 3.e-10_r8 * fac
        vmr(i,k,iso4_a1) = 2.e-10_r8

        ! post-adjrxt effective rate constants (jh2o2, jsoa_a1, jsoa_a2,
        ! OH_H2O2, usr_HO2_HO2 [vmr/s], DMS_NO3, DMS_OHa, SO2_OH_M,
        ! usr_DMS_OH, SOAE_tau)
        rxt(i,k,:)  = 0._r8
        rxt(i,k,1)  = 8.e-6_r8  * fac
        rxt(i,k,2)  = 4.e-7_r8  * fac
        rxt(i,k,3)  = 4.e-7_r8  * fac
        rxt(i,k,4)  = 1.8e-6_r8 * fac
        rxt(i,k,5)  = 1.e-14_r8 * fac
        rxt(i,k,6)  = 1.e-6_r8  * fac
        rxt(i,k,7)  = 4.e-6_r8  * fac
        rxt(i,k,8)  = 8.e-6_r8  * fac
        rxt(i,k,9)  = 2.e-6_r8  * fac
        rxt(i,k,10) = 1.157e-5_r8

        het(i,k,:)     = 0._r8
        het(i,k,iH2O2) = 2.e-6_r8 * fac
        het(i,k,iSO2)  = 1.e-5_r8 * fac

        ext(i,k,:)    = 0._r8
        ext(i,k,xSO2) = 3.e-13_r8 * fac
        ext(i,k,xso4_a1) = 1.e-15_r8 * fac
      end do
    end do

    ! total sulfur before (per cell, changed only by ext SO2 input and het
    ! SO2 loss; het loss makes it non-trivial, so track it separately below)
    call micm_solve( ncol, 1, delt, xhnm, tfld, pmid, vmr, rxt, het, ext )

    do k = 1, pver
      do i = 1, ncol
        cell = (k-1)*ncol + i
        fac = 1._r8 + 0.02_r8*(cell-1)

        j1 = 8.e-6_r8*fac;  j2 = 4.e-7_r8*fac
        k4 = 1.8e-6_r8*fac; e5 = 1.e-14_r8*fac
        k6 = 1.e-6_r8*fac;  k7 = 4.e-6_r8*fac
        k8 = 8.e-6_r8*fac;  k9 = 2.e-6_r8*fac
        k10 = 1.157e-5_r8
        hH2O2 = 2.e-6_r8*fac; hSO2 = 1.e-5_r8*fac
        eSO2 = 3.e-13_r8*fac; eso4 = 1.e-15_r8*fac

        ! DMS: pure first-order decay
        ld = k6 + k7 + k9
        call check('DMS', i, k, vmr(i,k,iDMS), 1.e-10_r8*fac*exp(-ld*delt))

        ! H2O2: production e5, loss j1+k4+het
        lh = j1 + k4 + hH2O2
        y  = e5/lh + (1.e-9_r8*fac - e5/lh)*exp(-lh*delt)
        call check('H2O2', i, k, vmr(i,k,iH2O2), y)

        ! SO2: production from DMS chain + ext, loss k8+het
        ls = k8 + hSO2
        a  = (k6 + k7 + 0.5_r8*k9)*1.e-10_r8*fac
        y  = eSO2/ls + (a/(ls-ld))*exp(-ld*delt) &
             + (5.e-10_r8*fac - eSO2/ls - a/(ls-ld))*exp(-ls*delt)
        call check('SO2', i, k, vmr(i,k,iSO2), y)

        ! SOAE -> SOAG
        call check('SOAE', i, k, vmr(i,k,iSOAE), 2.e-10_r8*fac*exp(-k10*delt))
        call check('SOAG', i, k, vmr(i,k,iSOAG), &
             1.e-11_r8 + 2.e-10_r8*fac*(1._r8 - exp(-k10*delt)))

        ! soa_a1: photolytic decay only
        call check('soa_a1', i, k, vmr(i,k,isoa_a1), 3.e-10_r8*fac*exp(-j2*delt))

        ! so4_a1: pure external forcing
        call check('so4_a1', i, k, vmr(i,k,iso4_a1), 2.e-10_r8 + eso4*delt)

        ! total-sulfur invariant: d(DMS+SO2+H2SO4)/dt = eSO2 - hSO2*SO2(t);
        ! integrate the exact SO2(t) from 0 to delt for the het-loss term
        block
          real(r8) :: c1, c2, c3, so2_int
          c1 = eSO2/ls
          c2 = a/(ls-ld)
          c3 = 5.e-10_r8*fac - c1 - c2
          so2_int = c1*delt + c2*(1._r8-exp(-ld*delt))/ld &
                  + c3*(1._r8-exp(-ls*delt))/ls
          s0 = 1.e-10_r8*fac + 5.e-10_r8*fac + 1.e-13_r8
          s1 = s0 + eSO2*delt - hSO2*so2_int
        end block
        call check('S-total', i, k, &
             vmr(i,k,iDMS)+vmr(i,k,iSO2)+vmr(i,k,iH2SO4), s1)
      end do
    end do

    deallocate(xhnm, tfld, pmid, vmr, rxt, het, ext)

  end subroutine run_chunk

  subroutine check(name, i, k, got, want)
    character(len=*), intent(in) :: name
    integer,  intent(in) :: i, k
    real(r8), intent(in) :: got, want
    real(r8) :: rel
    rel = abs(got - want)/max(abs(want), 1.e-30_r8)
    if (rel > tol) then
      print '(a,2i3,a,es15.7,a,es15.7,a,es9.2)', &
        'MISMATCH '//name//' cell(', i, k, ') got ', got, ' want ', want, ' rel ', rel
      nfail = nfail + 1
    end if
  end subroutine check

end program harness_micm
