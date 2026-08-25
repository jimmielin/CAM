!===============================================================================
! Heterogeneous uptake of acid gases on mineral dust, with calcite (CaCO3) as the
! reactive alkalinity (Leung et al.; uptake coefficients after Fairlie et al. 2010
! and Tang et al. 2017):
!   CaCO3 + SO2    -> CaSO4
!   CaCO3 + 2 HNO3 -> Ca(NO3)2
!   CaCO3 + H2SO4  -> CaSO4
! Once the calcite of a cell is exhausted H2SO4 keeps condensing on dust as sulfate.
! Uptake in each bin is capped at the calcite the bin holds, and the total over
! bins at the gas available, so calcite and the acid gases stay positive-definite.
! The oxidant uptake reactions on dust (N2O5, NO3, OH, HO2, H2O2, O3, CH2O) are
! part of the gas-phase mechanism (usr_*_dust rates in mo_usrrxt).
!
! Portable: sized by ncol/pver with no host dependencies. The host resolves the
! chemistry indices at init and supplies the per-bin dust surface area density and
! effective diameter at run time.
!===============================================================================
module dust_hetchem
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: dust_hetchem_init
  public :: dust_hetchem_run
  public :: has_dust_hetchem

  logical, protected :: has_dust_hetchem = .false. ! some bin carries calcite and both products

  ! chemistry (vmr) indices resolved by the host, < 1 where a species is not present
  integer :: nbins = 0
  integer, allocatable :: id_cal(:)        ! calcite CaCO3, per bin
  integer, allocatable :: id_cs4(:)        ! CaSO4, per bin
  integer, allocatable :: id_cn3(:)        ! Ca(NO3)2, per bin
  integer, allocatable :: id_so4(:)        ! sulfate, per bin
  logical, allocatable :: bin_active(:)    ! bin carries calcite and both products
  logical, allocatable :: bin_so4(:)       ! active bin that also carries sulfate
  integer :: id_so2 = -1
  integer :: id_hno3 = -1
  integer :: id_h2so4 = -1

  real(r8), parameter :: dg = 0.1_r8                       ! gas diffusivity, cm2/s (Dentener 1993)
  real(r8), parameter :: gamma_h2so4_dust = 0.10_r8        ! H2SO4 uptake coefficient, Fuchs-Sutugin corrected
  real(r8), parameter :: aciduptake_tuning_fact = 1.5_r8   ! on the RH-dependent SO2 and HNO3 uptake coefficients
  real(r8), parameter :: caco3_stoi = 0.5_r8               ! calcite consumed per HNO3 taken up (CaCO3 + 2 HNO3)
  ! minimum left in a reservoir when it is fully drawn down
  real(r8), parameter :: gas_min   = 1.e-20_r8
  real(r8), parameter :: caco3_min = 1.e-22_r8

contains

  !=============================================================================
  ! Stores the chemistry indices of the participating species. A bin takes part
  ! when it carries calcite and both products; the host passes < 1 where a
  ! species does not exist.
  !=============================================================================
  subroutine dust_hetchem_init( nbins_in, id_cal_in, id_cs4_in, id_cn3_in, id_so4_in, &
                                id_so2_in, id_hno3_in, id_h2so4_in, errmsg, errflg )

    integer, intent(in) :: nbins_in
    integer, intent(in) :: id_cal_in(:)   ! calcite index per bin
    integer, intent(in) :: id_cs4_in(:)   ! CaSO4 index per bin
    integer, intent(in) :: id_cn3_in(:)   ! Ca(NO3)2 index per bin
    integer, intent(in) :: id_so4_in(:)   ! sulfate index per bin
    integer, intent(in) :: id_so2_in
    integer, intent(in) :: id_hno3_in
    integer, intent(in) :: id_h2so4_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    nbins = nbins_in
    allocate( id_cal(nbins), id_cs4(nbins), id_cn3(nbins), id_so4(nbins) )
    allocate( bin_active(nbins), bin_so4(nbins) )
    id_cal(:) = id_cal_in(:nbins)
    id_cs4(:) = id_cs4_in(:nbins)
    id_cn3(:) = id_cn3_in(:nbins)
    id_so4(:) = id_so4_in(:nbins)
    id_so2   = id_so2_in
    id_hno3  = id_hno3_in
    id_h2so4 = id_h2so4_in

    bin_active(:) = id_cal(:) > 0 .and. id_cs4(:) > 0 .and. id_cn3(:) > 0
    bin_so4(:) = bin_active(:) .and. id_so4(:) > 0

    has_dust_hetchem = any(bin_active)
    if (has_dust_hetchem .and. (id_so2 < 1 .or. id_hno3 < 1 .or. id_h2so4 < 1)) then
       has_dust_hetchem = .false.
       errflg = 1
       errmsg = 'dust_hetchem_init: the mechanism carries calcite tracers but not SO2, HNO3 and H2SO4'
    end if

  end subroutine dust_hetchem_init

  !=============================================================================
  ! Applies one step of SO2, HNO3 and H2SO4 uptake on dust to the chemistry vmr
  ! array over levels beglev(i):endlev(i) (the troposphere; the dust surface area
  ! is zero elsewhere).
  !=============================================================================
  subroutine dust_hetchem_run( ncol, pver, dtime, beglev, endlev, tfld, relhum, sfc_dust, deff, vmr, &
                               gamma_so2, gamma_hno3, errmsg, errflg )

    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: pver
    real(r8), intent(in)    :: dtime            ! time step (s)
    integer,  intent(in)    :: beglev(:)        ! first level of the column range
    integer,  intent(in)    :: endlev(:)        ! last level of the column range
    real(r8), intent(in)    :: tfld(:,:)        ! temperature (K)
    real(r8), intent(in)    :: relhum(:,:)      ! relative humidity (fraction)
    real(r8), intent(in)    :: sfc_dust(:,:,:)  ! dust surface area density per bin (cm2/cm3)
    real(r8), intent(in)    :: deff(:,:,:)      ! dust effective diameter per bin (cm)
    real(r8), intent(inout) :: vmr(:,:,:)       ! chemistry species (vmr)
    real(r8), intent(out)   :: gamma_so2(:,:)   ! SO2 uptake coefficient on dust, diagnostic
    real(r8), intent(out)   :: gamma_hno3(:,:)  ! HNO3 uptake coefficient on dust, diagnostic
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer  :: i, k, l
    real(r8) :: sqrt_t
    real(r8) :: c_so2, c_hno3, c_h2so4         ! mean molecular speeds (cm/s)
    real(r8) :: gamma_so2_dust, gamma_hno3_dust
    real(r8) :: r_so2(nbins), r_hno3(nbins), r_h2so4(nbins)   ! first-order uptake rates per bin (1/s)
    real(r8) :: delta(nbins)                   ! gas taken up by each bin this step (vmr)
    real(r8) :: delta_ntot                     ! sum over bins
    real(r8) :: caco3_ntot                     ! calcite summed over bins (vmr)
    real(r8) :: gas_avail                      ! gas available for uptake this step, above gas_min

    errmsg = ''
    errflg = 0

    gamma_so2(:ncol,:pver) = 0._r8
    gamma_hno3(:ncol,:pver) = 0._r8

    col_loop: do i = 1, ncol
       ver_loop: do k = beglev(i), endlev(i)

          sqrt_t = sqrt( tfld(i,k) )
          c_so2   = 1.82e3_r8 * sqrt_t
          c_hno3  = 1.83e3_r8 * sqrt_t
          c_h2so4 = 1.47e3_r8 * sqrt_t

          ! SO2 and HNO3 uptake coefficients depend on relative humidity (in %)
          gamma_so2_dust  = gamma_so2_relhum_dust( relhum(i,k)*100._r8 )
          gamma_hno3_dust = gamma_hno3_relhum_dust( relhum(i,k)*100._r8 )
          gamma_so2(i,k)  = gamma_so2_dust
          gamma_hno3(i,k) = gamma_hno3_dust

          r_so2(:)   = uptake_rate( sfc_dust(i,k,:nbins), deff(i,k,:nbins), dg, c_so2,   aciduptake_tuning_fact * gamma_so2_dust )
          r_hno3(:)  = uptake_rate( sfc_dust(i,k,:nbins), deff(i,k,:nbins), dg, c_hno3,  aciduptake_tuning_fact * gamma_hno3_dust )
          r_h2so4(:) = uptake_rate( sfc_dust(i,k,:nbins), deff(i,k,:nbins), dg, c_h2so4, gamma_h2so4_dust )

          caco3_ntot = caco3_sum( vmr, i, k )

          if ( caco3_ntot > 5.e-20_r8 ) then

             !----------------------------------------------------------------
             ! SO2 uptake: CaCO3 + SO2 -> CaSO4
             !----------------------------------------------------------------
             delta(:) = r_so2(:) * vmr(i,k,id_so2) * dtime

             ! a bin cannot consume more calcite than it holds
             delta_ntot = 0._r8
             do l = 1, nbins
                if (.not. bin_active(l)) cycle
                delta(l) = min( delta(l), max( vmr(i,k,id_cal(l)) - caco3_min, 0._r8 ) )
                delta_ntot = delta_ntot + delta(l)
             end do

             ! the total cannot exceed the available gas; on a shortfall partition proportionally
             gas_avail = max( vmr(i,k,id_so2) - gas_min, 0._r8 )
             if ( delta_ntot > gas_avail ) then
                do l = 1, nbins
                   if (bin_active(l)) delta(l) = delta(l) * gas_avail / delta_ntot
                end do
                delta_ntot = gas_avail
             end if

             do l = 1, nbins
                if (.not. bin_active(l)) cycle
                vmr(i,k,id_cs4(l)) = vmr(i,k,id_cs4(l)) + delta(l)
                vmr(i,k,id_cal(l)) = vmr(i,k,id_cal(l)) - delta(l)
             end do
             vmr(i,k,id_so2) = vmr(i,k,id_so2) - delta_ntot

             !----------------------------------------------------------------
             ! HNO3 uptake: CaCO3 + 2 HNO3 -> Ca(NO3)2, calcite consumed and
             ! Ca(NO3)2 formed are caco3_stoi times the HNO3 taken up
             !----------------------------------------------------------------
             delta(:) = r_hno3(:) * vmr(i,k,id_hno3) * dtime

             delta_ntot = 0._r8
             do l = 1, nbins
                if (.not. bin_active(l)) cycle
                delta(l) = min( delta(l), max( vmr(i,k,id_cal(l)) - caco3_min, 0._r8 ) / caco3_stoi )
                delta_ntot = delta_ntot + delta(l)
             end do

             gas_avail = max( vmr(i,k,id_hno3) - gas_min, 0._r8 )
             if ( delta_ntot > gas_avail ) then
                do l = 1, nbins
                   if (bin_active(l)) delta(l) = delta(l) * gas_avail / delta_ntot
                end do
                delta_ntot = gas_avail
             end if

             do l = 1, nbins
                if (.not. bin_active(l)) cycle
                vmr(i,k,id_cn3(l)) = vmr(i,k,id_cn3(l)) + caco3_stoi*delta(l)
                vmr(i,k,id_cal(l)) = vmr(i,k,id_cal(l)) - caco3_stoi*delta(l)
             end do
             vmr(i,k,id_hno3) = vmr(i,k,id_hno3) - delta_ntot

          end if

          !-------------------------------------------------------------------
          ! H2SO4 condensation: consumes calcite while some is left
          ! (CaCO3 + H2SO4 -> CaSO4), otherwise continues as sulfate
          !-------------------------------------------------------------------
          caco3_ntot = caco3_sum( vmr, i, k )

          if ( (caco3_ntot > 5.e-20_r8) .and. (vmr(i,k,id_h2so4) > 2.e-20_r8) ) then

             delta(:) = r_h2so4(:) * vmr(i,k,id_h2so4) * dtime

             delta_ntot = 0._r8
             do l = 1, nbins
                if (.not. bin_active(l)) cycle
                delta(l) = min( delta(l), max( vmr(i,k,id_cal(l)) - caco3_min, 0._r8 ) )
                delta_ntot = delta_ntot + delta(l)
             end do

             gas_avail = max( vmr(i,k,id_h2so4) - gas_min, 0._r8 )
             if ( delta_ntot > gas_avail ) then
                do l = 1, nbins
                   if (bin_active(l)) delta(l) = delta(l) * gas_avail / delta_ntot
                end do
                delta_ntot = gas_avail
             end if

             do l = 1, nbins
                if (.not. bin_active(l)) cycle
                vmr(i,k,id_cs4(l)) = vmr(i,k,id_cs4(l)) + delta(l)
                vmr(i,k,id_cal(l)) = vmr(i,k,id_cal(l)) - delta(l)
             end do
             vmr(i,k,id_h2so4) = vmr(i,k,id_h2so4) - delta_ntot

          else if ( (caco3_ntot <= 5.e-20_r8) .and. (vmr(i,k,id_h2so4) > 2.e-20_r8) ) then

             delta(:) = r_h2so4(:) * vmr(i,k,id_h2so4) * dtime

             ! only bins that receive the sulfate count towards the gas limit
             delta_ntot = 0._r8
             do l = 1, nbins
                if (bin_so4(l)) delta_ntot = delta_ntot + delta(l)
             end do

             gas_avail = max( vmr(i,k,id_h2so4) - gas_min, 0._r8 )
             if ( delta_ntot > gas_avail ) then
                do l = 1, nbins
                   if (bin_so4(l)) delta(l) = delta(l) * gas_avail / delta_ntot
                end do
                delta_ntot = gas_avail
             end if

             do l = 1, nbins
                if (bin_so4(l)) vmr(i,k,id_so4(l)) = vmr(i,k,id_so4(l)) + delta(l)
             end do
             vmr(i,k,id_h2so4) = vmr(i,k,id_h2so4) - delta_ntot

          end if

       end do ver_loop
    end do col_loop

  end subroutine dust_hetchem_run

  !=============================================================================
  ! calcite summed over the bins that carry it, in ascending bin order
  !=============================================================================
  function caco3_sum( vmr, i, k ) result(total)

    real(r8), intent(in) :: vmr(:,:,:)
    integer,  intent(in) :: i, k
    real(r8) :: total

    integer :: l

    total = 0._r8
    do l = 1, nbins
       if (id_cal(l) > 0) total = total + vmr(i,k,id_cal(l))
    end do

  end function caco3_sum

  !=============================================================================
  ! First-order uptake rate (1/s) of a gas on one aerosol bin: the per-bin term
  ! of mo_usrrxt::hetrxtrate. Zero where the uptake coefficient is zero (the
  ! expression's 4/(c*gamma) would divide by zero).
  !=============================================================================
  elemental function uptake_rate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc        ! surface area density (cm2/cm3)
    real(r8), intent(in) :: dm_aer     ! diameter (cm)
    real(r8), intent(in) :: dg_gas     ! gas diffusivity (cm2/s)
    real(r8), intent(in) :: c_gas      ! mean molecular speed (cm/s)
    real(r8), intent(in) :: gamma_gas  ! uptake coefficient
    real(r8) :: rate

    if (gamma_gas > 0._r8) then
       rate = sfc / (0.5_r8*dm_aer/dg_gas + (4._r8/(c_gas*gamma_gas)))
    else
       rate = 0._r8
    end if

  end function uptake_rate

  !=============================================================================
  ! SO2 uptake coefficient on dust as a function of relative humidity (%), after
  ! Preszler-Prince et al. (2007) as implemented in GEOS-Chem (Fairlie et al. 2010)
  !=============================================================================
  function gamma_so2_relhum_dust( RH ) result(gamma_so2_dust_rh)

    real(r8), intent(in) :: RH                  ! relative humidity (%)
    real(r8) :: gamma_so2_dust_rh

    real(r8), parameter :: GAMMA_SO2 = 2.5e-4_r8

    if (RH >= 90._r8 ) then
       gamma_so2_dust_rh = GAMMA_SO2 * 2._r8
    else if (RH >= 84._r8 ) then
       gamma_so2_dust_rh = GAMMA_SO2 * (0.5_r8 + 1.5_r8 * (RH - 84._r8) / (90._r8 - 84._r8))
    else if (RH >= 76._r8 ) then
       gamma_so2_dust_rh = GAMMA_SO2 * (0.16_r8 + 0.34_r8 * (RH - 76._r8) / (84._r8 - 76._r8))
    else if (RH >= 33._r8 ) then
       gamma_so2_dust_rh = GAMMA_SO2 * (0.03_r8 + 0.13_r8 * (RH - 33._r8) / (76._r8 - 33._r8))
    else if (RH >= 20._r8 ) then
       gamma_so2_dust_rh = GAMMA_SO2 * 0.03_r8
    else                        ! 0.0 below 20%
       gamma_so2_dust_rh = GAMMA_SO2 * 0.0_r8
    endif

  end function gamma_so2_relhum_dust

  !=============================================================================
  ! HNO3 uptake coefficient on dust as a function of relative humidity (%), after
  ! Liu et al. (2007) as implemented in GEOS-Chem (Fairlie et al. 2010)
  !=============================================================================
  function gamma_hno3_relhum_dust( RH ) result(gamma_hno3_dust_rh)

    real(r8), intent(in) :: RH                  ! relative humidity (%)
    real(r8) :: gamma_hno3_dust_rh

    real(r8), parameter :: GAMMA_HNO3 = 5e-4_r8

    if (RH >= 80.0_r8 ) then
       gamma_hno3_dust_rh = GAMMA_HNO3 * 2.0_r8
    else if (RH >= 70.0_r8 ) then
       gamma_hno3_dust_rh = GAMMA_HNO3 * ( 1.3_r8 + 0.7_r8 * (RH - 70.0_r8) / 10.0_r8 )
    else if (RH >= 60.0_r8 ) then
       gamma_hno3_dust_rh = GAMMA_HNO3 * ( 1.0_r8 + 0.3_r8 * (RH - 60.0_r8) / 10.0_r8 )
    else if (RH >= 50.0_r8 ) then
       gamma_hno3_dust_rh = GAMMA_HNO3 * ( 0.7_r8 + 0.3_r8 * (RH - 50.0_r8) / 10.0_r8 )
    else if (RH >= 30.0_r8 ) then
       gamma_hno3_dust_rh = GAMMA_HNO3 * ( 0.19_r8 + 0.255_r8 * (RH - 30.0_r8) / 10.0_r8 )
    else if (RH >= 10.0_r8 ) then
       gamma_hno3_dust_rh = GAMMA_HNO3 * ( 0.03_r8 + 0.08_r8 * (RH - 10.0_r8) / 10.0_r8 )
    else
       ! 0.0 below 10%
       gamma_hno3_dust_rh = 0.0_r8
    end if

  end function gamma_hno3_relhum_dust

end module dust_hetchem
