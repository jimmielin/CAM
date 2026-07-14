! Grid-box mean total (stratiform + convective) cloud liquid and ice for
! radiation: blends deep and shallow convective in-cloud condensate with the
! stratiform condensate.
!
! Portable science core split out of conv_water.F90 (conv_water_4rad); the
! CAM wrapper (pbuf marshaling, namelist, history output) remains there.
! This file is copied verbatim (r8 -> kind_phys only) into
! atmospheric_physics as a CCPP scheme; keep it host-independent (no CAM
! modules beyond shr_kind_mod, array sizing by ncol/pver arguments, host
! constants passed as arguments) and mirror any science change to the copy.
! - one_mom_clouds selects the one-moment (RK) ice absorptivity (rei used
!   unclamped); CAM passes microp_scheme=='RK'.
! - The quantities the wrapper outflds are returned as diagnostic-only
!   intent(out) arguments.
module convective_cloud_water

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: convective_cloud_water_run

contains

!> \section arg_table_convective_cloud_water_run Argument Table
!! \htmlinclude convective_cloud_water_run.html
  subroutine convective_cloud_water_run( &
       ncol, pver, &
       conv_water_in_rad, frac_limit, &
       one_mom_clouds, &
       gravit, &
       pdel, ls_liq, ls_ice, &
       sh_icwmr, dp_icwmr, &
       sh_frac, dp_frac, ast, rei, &
       totg_liq, totg_ice, &
       conv_liq, conv_ice, tot_liq, tot_ice, &
       totg_liq_sh, totg_liq_dp, totg_ice_sh, totg_ice_dp, &
       fresh, fredp, frecu, fretot, &
       errmsg, errflg)

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: conv_water_in_rad  ! 0==> No; 1==> Yes-Arithmetic average;
                                                        ! 2==> Yes-Average in emissivity.
    real(r8),  intent(in)  :: frac_limit         ! minimum cloud fraction [fraction]
    logical,          intent(in)  :: one_mom_clouds     ! one-moment (RK) microphysics: use rei unclamped in the
                                                        ! ice absorptivity; two-moment clamps rei to [13,130] micron
    real(r8),  intent(in)  :: gravit             ! gravitational acceleration [m s-2]

    real(r8),  intent(in)  :: pdel(:,:)          ! (ncol,pver) moist pressure difference across layer [Pa]
    real(r8),  intent(in)  :: ls_liq(:,:)        ! (ncol,pver) large-scale contributions to GBA cloud liq [kg kg-1]
    real(r8),  intent(in)  :: ls_ice(:,:)        ! (ncol,pver) large-scale contributions to GBA cloud ice [kg kg-1]

    real(r8),  intent(in)  :: sh_icwmr(:,:)      ! (ncol,pver) shallow conv. cloud water [kg kg-1]
    real(r8),  intent(in)  :: dp_icwmr(:,:)      ! (ncol,pver) deep conv. cloud water [kg kg-1]
    real(r8),  intent(in)  :: sh_frac(:,:)       ! (ncol,pver) shallow convective cloud fraction [fraction]
    real(r8),  intent(in)  :: dp_frac(:,:)       ! (ncol,pver) deep convective cloud fraction [fraction]
    real(r8),  intent(in)  :: ast(:,:)           ! (ncol,pver) physical liquid+ice stratus cloud fraction [fraction]
    real(r8),  intent(in)  :: rei(:,:)           ! (ncol,pver) ice effective drop size [micron]

    real(r8),  intent(out) :: totg_liq(:,:)      ! (ncol,pver) grid box total cloud liquid mixing ratio [kg kg-1]
    real(r8),  intent(out) :: totg_ice(:,:)      ! (ncol,pver) grid box total cloud ice mixing ratio [kg kg-1]

    ! Diagnostic-only outputs (CAM history: ICLMRCU/ICIMRCU/ICLMRTOT/ICIMRTOT/
    ! GCLMRDP/GCIMRDP/GCLMRSH/GCIMRSH/FRESH/FREDP/FRECU/FRETOT)
    real(r8),  intent(out) :: conv_liq(:,:)      ! (ncol,pver) convective contributions to IC cloud liquid [kg kg-1]
    real(r8),  intent(out) :: conv_ice(:,:)      ! (ncol,pver) convective contributions to IC cloud ice [kg kg-1]
    real(r8),  intent(out) :: tot_liq(:,:)       ! (ncol,pver) total IC liquid [kg kg-1]
    real(r8),  intent(out) :: tot_ice(:,:)       ! (ncol,pver) total IC ice [kg kg-1]
    real(r8),  intent(out) :: totg_liq_sh(:,:)   ! (ncol,pver) grid-mean LWP from shallow convective cloud [kg kg-1]
    real(r8),  intent(out) :: totg_liq_dp(:,:)   ! (ncol,pver) grid-mean LWP from deep convective cloud [kg kg-1]
    real(r8),  intent(out) :: totg_ice_sh(:,:)   ! (ncol,pver) grid-mean IWP from shallow convective cloud [kg kg-1]
    real(r8),  intent(out) :: totg_ice_dp(:,:)   ! (ncol,pver) grid-mean IWP from deep convective cloud [kg kg-1]
    real(r8),  intent(out) :: fresh(:,:)         ! (ncol,pver) fractional occurrence of shallow cumulus [1]
    real(r8),  intent(out) :: fredp(:,:)         ! (ncol,pver) fractional occurrence of deep cumulus [1]
    real(r8),  intent(out) :: frecu(:,:)         ! (ncol,pver) fractional occurrence of cumulus [1]
    real(r8),  intent(out) :: fretot(:,:)        ! (ncol,pver) fractional occurrence of cloud [1]

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! --------------- !
    ! Local Workspace !
    ! --------------- !

    integer  :: i,k                                ! Lon, lev indices buff stuff.
    real(r8) :: cu_icwmr                    ! Convective  water for this grid-box.
    real(r8) :: ls_icwmr                    ! Large-scale water for this grid-box.
    real(r8) :: tot_icwmr                   ! Large-scale water for this grid-box.
    real(r8) :: ls_frac                     ! Large-scale cloud frac for this grid-box.
    real(r8) :: tot0_frac, cu0_frac, dp0_frac, sh0_frac
    real(r8) :: kabs, kabsi, kabsl, alpha, dp0, sh0, ic_limit
    real(r8) :: wrk1

    ! --------- !
    ! Parameter !
    ! --------- !

    parameter( kabsl = 0.090361_r8, ic_limit = 1.e-12_r8 )

    errmsg = ''
    errflg = 0

    ! The loop below fully assigns the outputs only for conv_water_in_rad
    ! = 1 or 2 (the select-case default leaves cu_icwmr/tot_icwmr unset;
    ! CAM's endrun there is commented out upstream). CAM only calls
    ! conv_water_4rad when conv_water_in_rad /= 0, but a CCPP suite runs
    ! this scheme unconditionally: for mode 0 return all-zero fields (they
    ! are unused downstream in that mode); reject any other value.
    if (conv_water_in_rad /= 1 .and. conv_water_in_rad /= 2) then
       totg_liq(:,:)    = 0._r8
       totg_ice(:,:)    = 0._r8
       conv_liq(:,:)    = 0._r8
       conv_ice(:,:)    = 0._r8
       tot_liq(:,:)     = 0._r8
       tot_ice(:,:)     = 0._r8
       totg_liq_sh(:,:) = 0._r8
       totg_liq_dp(:,:) = 0._r8
       totg_ice_sh(:,:) = 0._r8
       totg_ice_dp(:,:) = 0._r8
       fresh(:,:)       = 0._r8
       fredp(:,:)       = 0._r8
       frecu(:,:)       = 0._r8
       fretot(:,:)      = 0._r8
       if (conv_water_in_rad /= 0) then
          errflg = 1
          errmsg = 'convective_cloud_water_run: invalid conv_water_in_rad (must be 0, 1, or 2)'
       end if
       return
    end if

    ! --------------------------------------------------------------- !
    ! Loop through grid-boxes and determine:                          !
    ! 1. Effective mean in-cloud convective ice/liquid (deep+shallow) !
    ! 2. Effective mean in-cloud total ice/liquid (ls+convective)     !
    ! --------------------------------------------------------------- !

    fresh(:,:)  = 0._r8
    fredp(:,:)  = 0._r8
    frecu(:,:)  = 0._r8
    fretot(:,:) = 0._r8

    do k = 1, pver
    do i = 1, ncol

       if( sh_frac(i,k) <= frac_limit .or. sh_icwmr(i,k) <= ic_limit ) then
           sh0_frac = 0._r8
       else
           sh0_frac = sh_frac(i,k)
       endif
       if( dp_frac(i,k) <= frac_limit .or. dp_icwmr(i,k) <= ic_limit ) then
           dp0_frac = 0._r8
       else
           dp0_frac = dp_frac(i,k)
       endif
       cu0_frac = sh0_frac + dp0_frac

     ! For the moment calculate the emissivity based upon the ls clouds ice fraction

       wrk1 = min(1._r8,max(0._r8, ls_ice(i,k)/(ls_ice(i,k)+ls_liq(i,k)+1.e-36_r8)))

       if( ( cu0_frac < frac_limit ) .or. ( ( sh_icwmr(i,k) + dp_icwmr(i,k) ) < ic_limit ) ) then

             cu0_frac = 0._r8
             cu_icwmr = 0._r8

             ls_frac = ast(i,k)
             if( ls_frac < frac_limit ) then
                 ls_frac  = 0._r8
                 ls_icwmr = 0._r8
             else
                 ls_icwmr = ( ls_liq(i,k) + ls_ice(i,k) )/max(frac_limit,ls_frac) ! Convert to IC value.
             end if

             tot0_frac = ls_frac
             tot_icwmr = ls_icwmr

       else

           ! Select radiation constants (effective radii) for emissivity averaging.

             if( one_mom_clouds ) then
                kabsi = 0.005_r8 + 1._r8/rei(i,k)
             else
                kabsi = 0.005_r8 + 1._r8/min(max(13._r8,rei(i,k)),130._r8)
             endif
             kabs  = kabsl * ( 1._r8 - wrk1 ) + kabsi * wrk1
             alpha = -1.66_r8*kabs*pdel(i,k)/gravit*1000.0_r8

           ! Selecting cumulus in-cloud water.

             select case (conv_water_in_rad) ! Type of average
             case (1) ! Area weighted arithmetic average
                cu_icwmr = ( sh0_frac * sh_icwmr(i,k) + dp0_frac*dp_icwmr(i,k))/max(frac_limit,cu0_frac)
             case (2)
                sh0 = exp(alpha*sh_icwmr(i,k))
                dp0 = exp(alpha*dp_icwmr(i,k))
                cu_icwmr = log((sh0_frac*sh0+dp0_frac*dp0)/max(frac_limit,cu0_frac))
                cu_icwmr = cu_icwmr/alpha
             case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
             end select

           ! Selecting total in-cloud water.
           ! Attribute large-scale/convective area fraction differently from default.

             ls_frac   = ast(i,k)
             ls_icwmr  = (ls_liq(i,k) + ls_ice(i,k))/max(frac_limit,ls_frac) ! Convert to IC value.
             tot0_frac = (ls_frac + cu0_frac)

             select case (conv_water_in_rad) ! Type of average
             case (1) ! Area weighted 'arithmetic in emissivity' average
                tot_icwmr = (ls_frac*ls_icwmr + cu0_frac*cu_icwmr)/max(frac_limit,tot0_frac)
             case (2)
                tot_icwmr = log((ls_frac*exp(alpha*ls_icwmr)+cu0_frac*exp(alpha*cu_icwmr))/max(frac_limit,tot0_frac))
                tot_icwmr = tot_icwmr/alpha
             case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
             end select

       end if

     ! Repartition convective cloud water into liquid and ice phase.
     ! Currently, this partition is made using the ice fraction of stratus condensate.
     ! In future, we should use ice fraction explicitly computed from the convection scheme.

       conv_ice(i,k) = cu_icwmr * wrk1
       conv_liq(i,k) = cu_icwmr * (1._r8-wrk1)

       tot_ice(i,k)  = tot_icwmr * wrk1
       tot_liq(i,k)  = tot_icwmr * (1._r8-wrk1)

       totg_ice(i,k) = tot0_frac * tot_icwmr * wrk1
       totg_liq(i,k) = tot0_frac * tot_icwmr * (1._r8-wrk1)

       ! Grid-mean convective water
       totg_ice_sh(i,k)  = sh0_frac * sh_icwmr(i,k) * wrk1
       totg_ice_dp(i,k)  = dp0_frac * dp_icwmr(i,k) * wrk1
       totg_liq_sh(i,k)  = sh0_frac * sh_icwmr(i,k) * (1._r8-wrk1)
       totg_liq_dp(i,k)  = dp0_frac * dp_icwmr(i,k) * (1._r8-wrk1)
       if( sh0_frac > frac_limit ) then
           fresh(i,k) = 1._r8
       endif
       if( dp0_frac > frac_limit ) then
           fredp(i,k) = 1._r8
       endif
       if( cu0_frac > frac_limit ) then
           frecu(i,k) = 1._r8
       endif
       if( tot0_frac > frac_limit ) then
           fretot(i,k) = 1._r8
       endif

    end do
    end do

  end subroutine convective_cloud_water_run

end module convective_cloud_water
