module wetdep

!-----------------------------------------------------------------------
!
! Wet deposition routines for both aerosols and gas phase constituents.
!
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8

implicit none
save
private

public :: wetdepa_v1  ! scavenging codes for very soluble aerosols -- CAM4 version
public :: wetdepa_v2  ! scavenging codes for very soluble aerosols -- CAM5 version
public :: clddiag     ! calc of cloudy volume and rain mixing ratio

real(r8), parameter :: cmftau = 3600._r8
real(r8), parameter :: rhoh2o = 1000._r8            ! density of water
real(r8), parameter :: omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero

!==============================================================================
contains
!==============================================================================

subroutine clddiag(t, pmid, pdel, cmfdqr, evapc, &
                   cldt, cldcu, cldst, cme, evapr, &
                   prain, cldv, cldvcu, cldvst, rain, &
                   ncol, pver, gravit, tmelt, rair)

   ! ------------------------------------------------------------------------------------
   ! Estimate the cloudy volume which is occupied by rain or cloud water as
   ! the max between the local cloud amount or the
   ! sum above of (cloud*positive precip production)      sum total precip from above
   !              ----------------------------------   x ------------------------
   ! sum above of     (positive precip           )        sum positive precip from above
   ! Author: P. Rasch
   !         Sungsu Park. Mar.2010
   ! ------------------------------------------------------------------------------------

   ! Input arguments:
   real(r8), intent(in) :: t(:,:)        ! temperature (K)
   real(r8), intent(in) :: pmid(:,:)     ! pressure at layer midpoints
   real(r8), intent(in) :: pdel(:,:)     ! pressure difference across layers
   real(r8), intent(in) :: cmfdqr(:,:)   ! dq/dt due to convective rainout
   real(r8), intent(in) :: evapc(:,:)    ! Evaporation rate of convective precipitation ( >= 0 )
   real(r8), intent(in) :: cldt(:,:)    ! total cloud fraction
   real(r8), intent(in) :: cldcu(:,:)    ! Cumulus cloud fraction
   real(r8), intent(in) :: cldst(:,:)    ! Stratus cloud fraction
   real(r8), intent(in) :: cme(:,:)      ! rate of cond-evap within the cloud
   real(r8), intent(in) :: evapr(:,:)    ! rate of evaporation of falling precipitation (kg/kg/s)
   real(r8), intent(in) :: prain(:,:)    ! rate of conversion of condensate to precipitation (kg/kg/s)
   integer, intent(in) :: ncol
   integer, intent(in) :: pver
   real(r8), intent(in) :: gravit       ! gravitational acceleration (m/s2)
   real(r8), intent(in) :: tmelt        ! freezing point of water (K)
   real(r8), intent(in) :: rair         ! dry air gas constant (J/K/kg)

   ! Output arguments:
   real(r8), intent(out) :: cldv(:,:)     ! fraction occupied by rain or cloud water
   real(r8), intent(out) :: cldvcu(:,:)   ! Convective precipitation volume
   real(r8), intent(out) :: cldvst(:,:)   ! Stratiform precipitation volume
   real(r8), intent(out) :: rain(:,:)     ! mixing ratio of rain (kg/kg)

   ! Local variables:
   integer  i, k
   real(r8) convfw         ! used in fallspeed calculation; taken from findmcnew
   real(r8) sumppr(ncol)        ! precipitation rate (kg/m2-s)
   real(r8) sumpppr(ncol)       ! sum of positive precips from above
   real(r8) cldv1(ncol)         ! precip weighted cloud fraction from above
   real(r8) lprec                ! local production rate of precip (kg/m2/s)
   real(r8) lprecp               ! local production rate of precip (kg/m2/s) if positive
   real(r8) rho                  ! air density
   real(r8) vfall
   real(r8) sumppr_cu(ncol)     ! Convective precipitation rate (kg/m2-s)
   real(r8) sumpppr_cu(ncol)    ! Sum of positive convective precips from above
   real(r8) cldv1_cu(ncol)      ! Convective precip weighted convective cloud fraction from above
   real(r8) lprec_cu             ! Local production rate of convective precip (kg/m2/s)
   real(r8) lprecp_cu            ! Local production rate of convective precip (kg/m2/s) if positive
   real(r8) sumppr_st(ncol)     ! Stratiform precipitation rate (kg/m2-s)
   real(r8) sumpppr_st(ncol)    ! Sum of positive stratiform precips from above
   real(r8) cldv1_st(ncol)      ! Stratiform precip weighted stratiform cloud fraction from above
   real(r8) lprec_st             ! Local production rate of stratiform precip (kg/m2/s)
   real(r8) lprecp_st            ! Local production rate of stratiform precip (kg/m2/s) if positive
   ! -----------------------------------------------------------------------

   convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)
   do i=1,ncol
      sumppr(i) = 0._r8
      cldv1(i) = 0._r8
      sumpppr(i) = 1.e-36_r8
      sumppr_cu(i)  = 0._r8
      cldv1_cu(i)   = 0._r8
      sumpppr_cu(i) = 1.e-36_r8
      sumppr_st(i)  = 0._r8
      cldv1_st(i)   = 0._r8
      sumpppr_st(i) = 1.e-36_r8
   end do

   do k = 1,pver
      do i = 1,ncol
         cldv(i,k) = &
            max(min(1._r8, &
            cldv1(i)/sumpppr(i) &
            )*sumppr(i)/sumpppr(i), &
            cldt(i,k) &
            )
         lprec = pdel(i,k)/gravit &
            *(prain(i,k)+cmfdqr(i,k)-evapr(i,k))
         lprecp = max(lprec,1.e-30_r8)
         cldv1(i) = cldv1(i)  + cldt(i,k)*lprecp
         sumppr(i) = sumppr(i) + lprec
         sumpppr(i) = sumpppr(i) + lprecp

         ! For convective precipitation volume at the top interface of each layer. Neglect the current layer.
         cldvcu(i,k)   = max(min(1._r8,cldv1_cu(i)/sumpppr_cu(i))*(sumppr_cu(i)/sumpppr_cu(i)),0._r8)
         lprec_cu      = (pdel(i,k)/gravit)*(cmfdqr(i,k)-evapc(i,k))
         lprecp_cu     = max(lprec_cu,1.e-30_r8)
         cldv1_cu(i)   = cldv1_cu(i) + cldcu(i,k)*lprecp_cu
         sumppr_cu(i)  = sumppr_cu(i) + lprec_cu
         sumpppr_cu(i) = sumpppr_cu(i) + lprecp_cu

         ! For stratiform precipitation volume at the top interface of each layer. Neglect the current layer.
         cldvst(i,k)   = max(min(1._r8,cldv1_st(i)/sumpppr_st(i))*(sumppr_st(i)/sumpppr_st(i)),0._r8)
         lprec_st      = (pdel(i,k)/gravit)*(prain(i,k)-evapr(i,k))
         lprecp_st     = max(lprec_st,1.e-30_r8)
         cldv1_st(i)   = cldv1_st(i) + cldst(i,k)*lprecp_st
         sumppr_st(i)  = sumppr_st(i) + lprec_st
         sumpppr_st(i) = sumpppr_st(i) + lprecp_st

         rain(i,k) = 0._r8
         if(t(i,k) > tmelt) then
            rho = pmid(i,k)/(rair*t(i,k))
            vfall = convfw/sqrt(rho)
            rain(i,k) = sumppr(i)/(rho*vfall)
            if (rain(i,k) < 1.e-14_r8) rain(i,k) = 0._r8
         endif
      end do
   end do

end subroutine clddiag

!==============================================================================

! This is the CAM5 version of wetdepa.

subroutine wetdepa_v2(                                  &
   p, q, pdel, cldt, cldc,                              &
   cmfdqr, evapc, conicw, precs, conds,                 &
   evaps, cwat, tracer, deltat, scavt,                  &
   iscavt, cldvcu, cldvst, dlf, fracis,                 &
   sol_fact, ncol, scavcoef, gravit, pver, errmsg, errflg, &
   is_strat_cloudborne, qqcw,                            &
   f_act_conv, icscavt, isscavt, bcscavt, bsscavt,      &
   convproc_do_aer, rcscavt, rsscavt,                   &
   sol_facti_in, sol_factic_in, convproc_do_evaprain_atonce_in, bergso_in )

   !-----------------------------------------------------------------------
   !
   ! scavenging code for very soluble aerosols
   !
   !-----------------------------------------------------------------------

   real(r8), intent(in) ::&
      p(:,:),        &! pressure
      q(:,:),        &! moisture
      pdel(:,:),     &! pressure thikness
      cldt(:,:),     &! total cloud fraction
      cldc(:,:),     &! convective cloud fraction
      cmfdqr(:,:),   &! rate of production of convective precip
      evapc(:,:),    &! Evaporation rate of convective precipitation
      conicw(:,:),   &! convective cloud water
      cwat(:,:),     &! cloud water amount
      precs(:,:),    &! rate of production of stratiform precip
      conds(:,:),    &! rate of production of condensate
      evaps(:,:),    &! rate of evaporation of precip
      cldvcu(:,:),   &! Convective precipitation area at the top interface of each layer
      cldvst(:,:),   &! Stratiform precipitation area at the top interface of each layer
      dlf(:,:),      &! Detrainment of convective condensate [kg/kg/s]
      deltat,               &! time step
      tracer(:,:)     ! trace species

   ! If subroutine is called with just sol_fact:
   !    sol_fact is used for both in- and below-cloud scavenging
   ! If subroutine is called with optional argument sol_facti_in:
   !    sol_fact  is used for below cloud scavenging
   !    sol_facti is used for in cloud scavenging

   real(r8), intent(in)  :: sol_fact(:,:)
   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: scavcoef(:,:) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
   real(r8), intent(in)  :: gravit       ! gravitational acceleration (m/s2)
   integer,  intent(in)  :: pver
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg
   real(r8), intent(out) ::&
      scavt(:,:),   &! scavenging tend
      iscavt(:,:),  &! incloud scavenging tends
      fracis(:,:)    ! fraction of species not scavenged

   ! Setting is_strat_cloudborne=.true. indicates that tracer is stratiform-cloudborne aerosol.
   !   This is only used by MAM code.  The optional args qqcw and f_act_conv are not referenced
   !   in this case.
   ! Setting is_strat_cloudborne=.false. is being used to indicate that the tracers are the
   !   interstitial modal aerosols.  In this case the optional qqcw (the cloud borne mixing ratio
   !   corresponding to the interstitial aerosol) must be provided, as well as the optional f_act_conv.
   logical,  intent(in), optional :: is_strat_cloudborne
   real(r8), intent(in), optional :: qqcw(:,:)
   real(r8), intent(in), optional :: f_act_conv(:,:)

   real(r8), intent(in), optional :: sol_facti_in(:,:)   ! solubility factor (frac of aerosol scavenged in cloud)
   real(r8), intent(in), optional :: sol_factic_in(:,:)  ! sol_facti_in for convective clouds


   real(r8), intent(out), optional :: icscavt(:,:)     ! incloud, convective
   real(r8), intent(out), optional :: isscavt(:,:)     ! incloud, stratiform
   real(r8), intent(out), optional :: bcscavt(:,:)     ! below cloud, convective
   real(r8), intent(out), optional :: bsscavt(:,:)     ! below cloud, stratiform

   ! Setting convproc_do_aer=.true. removes the resuspension term from bcscavt and
   ! bsscavt and returns those terms as rcscavt and rsscavt respectively.
   logical,  intent(in),  optional :: convproc_do_aer
   real(r8), intent(out), optional :: rcscavt(:,:)     ! resuspension, convective
   real(r8), intent(out), optional :: rsscavt(:,:)     ! resuspension, stratiform
   logical,  intent(in),  optional :: convproc_do_evaprain_atonce_in
   real(r8), intent(in),  optional :: bergso_in(:,:)

   ! local variables

   integer :: i, k
   logical :: out_resuspension

   real(r8) :: clds(ncol)          ! stratiform cloud fraction
   real(r8) :: fracev(ncol)        ! fraction of precip from above that is evaporating
   real(r8) :: fracev_cu(ncol)     ! Fraction of convective precip from above that is evaporating
   real(r8) :: fracp(ncol)         ! fraction of cloud water converted to precip
   real(r8) :: pdog(ncol)          ! work variable (pdel/gravit)
   real(r8) :: rpdog(ncol)         ! work variable (gravit/pdel)
   real(r8) :: precabc(ncol)       ! conv precip from above (work array)
   real(r8) :: precabs(ncol)       ! strat precip from above (work array)
   real(r8) :: rat(ncol)           ! ratio of amount available to amount removed
   real(r8) :: scavab(ncol)        ! scavenged tracer flux from above (work array)
   real(r8) :: scavabc(ncol)       ! scavenged tracer flux from above (work array)
   real(r8) :: srcc(ncol)          ! tend for convective rain
   real(r8) :: srcs(ncol)          ! tend for stratiform rain
   real(r8) :: srct(ncol)          ! work variable

   real(r8) :: fins(ncol)          ! fraction of rem. rate by strat rain
   real(r8) :: finc(ncol)          ! fraction of rem. rate by conv. rain
   real(r8) :: conv_scav_ic(ncol)  ! convective scavenging incloud
   real(r8) :: conv_scav_bc(ncol)  ! convective scavenging below cloud
   real(r8) :: st_scav_ic(ncol)    ! stratiform scavenging incloud
   real(r8) :: st_scav_bc(ncol)    ! stratiform scavenging below cloud

   real(r8) :: odds(ncol)          ! limit on removal rate (proportional to prec)
   real(r8) :: dblchek(ncol)
   logical :: found

   real(r8) :: trac_qqcw(ncol)
   real(r8) :: tracer_incu(ncol)
   real(r8) :: tracer_mean(ncol)

   ! For stratiform cloud, cloudborne aerosol is treated explicitly,
   !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
   ! For convective cloud, cloudborne aerosol is not treated explicitly,
   !    and sol_factic is 1.0 for both cloudborne and interstitial.

   real(r8) :: sol_facti(ncol,pver)  ! in cloud fraction of aerosol scavenged
   real(r8) :: sol_factb(ncol,pver)  ! below cloud fraction of aerosol scavenged
   real(r8) :: sol_factic(ncol,pver) ! in cloud fraction of aerosol scavenged for convective clouds

   real(r8) :: rdeltat
   logical  :: convproc_do_evaprain_atonce

   ! ------------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   if (present(convproc_do_evaprain_atonce_in)) then
      convproc_do_evaprain_atonce = convproc_do_evaprain_atonce_in
   else
      convproc_do_evaprain_atonce = .false.
   endif

   ! default (if other sol_facts aren't in call, set all to required sol_fact)
   sol_facti = sol_fact
   sol_factb = sol_fact

   if ( present(sol_facti_in) )  sol_facti = sol_facti_in

   sol_factic  = sol_facti
   if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in

   ! Determine whether resuspension fields are output.
   out_resuspension = .false.
   if (present(convproc_do_aer)) then
      if (convproc_do_aer) then
         if (present(bcscavt) .and. present(bsscavt) .and. &
             present(rcscavt) .and. present(rsscavt) ) then
            out_resuspension = .true.
         else
            errmsg = 'wetdepa_v2: bcscavt, bsscavt, rcscavt, rsscavt'// &
                     ' must be present when convproc_do_aero true'
            errflg = 1
            return
         end if
      end if
   end if

   ! this section of code is for highly soluble aerosols,
   ! the assumption is that within the cloud that
   ! all the tracer is in the cloud water
   !
   ! for both convective and stratiform clouds,
   ! the fraction of cloud water converted to precip defines
   ! the amount of tracer which is pulled out.

   precabs(:ncol) = 0.0_r8
   precabc(:ncol) = 0.0_r8
   scavab(:ncol)  = 0.0_r8
   scavabc(:ncol) = 0.0_r8

   do k = 1, pver
      do i = 1, ncol

         clds(i)  = cldt(i,k) - cldc(i,k)
         pdog(i)  = pdel(i,k)/gravit
         rpdog(i) = gravit/pdel(i,k)
         rdeltat  = 1.0_r8/deltat

         ! ****************** Evaporation **************************
         ! calculate the fraction of strat precip from above
         !                 which evaporates within this layer
         fracev(i) = evaps(i,k)*pdog(i) &
                     /max(1.e-12_r8,precabs(i))

         ! If resuspending aerosol only when all the rain has totally
         ! evaporated then zero out any aerosol tendency for partial
         ! evaporation.
         if (convproc_do_evaprain_atonce .and. (evaps(i,k)*pdog(i)/=precabs(i))) then
            fracev(i) = 0._r8
         end if

         ! trap to ensure reasonable ratio bounds
         fracev(i) = max(0._r8,min(1._r8,fracev(i)))

         ! Same as above but convective precipitation part
         fracev_cu(i) = evapc(i,k)*pdog(i)/max(1.e-12_r8,precabc(i))
         fracev_cu(i) = max(0._r8,min(1._r8,fracev_cu(i)))

         ! ****************** Convection ***************************
         !
         ! set odds proportional to fraction of the grid box that is swept by the
         ! precipitation =precabc/rhoh20*(area of sphere projected on plane
         !                                /volume of sphere)*deltat
         ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
         ! unless the fraction of the area that is cloud is less than odds, in which
         ! case use the cloud fraction (assumes precabs is in kg/m2/s)
         ! is really: precabs*3/4/1000./1e-3*deltat
         ! here I use .1 from Balkanski
         !
         ! use a local rate of convective rain production for incloud scav
         !
         ! Fraction of convective cloud water converted to rain.  This version is used
         ! in 2 of the 3 branches below before fracp is reused in the stratiform calc.
         ! NB: In below formula for fracp conicw is a LWC/IWC that has already
         !     precipitated out, i.e., conicw does not contain precipitation

         fracp(i) = cmfdqr(i,k)*deltat / &
                    max( 1.e-12_r8, cldc(i,k)*conicw(i,k) + (cmfdqr(i,k)+dlf(i,k))*deltat )
         fracp(i) = max( min( 1._r8, fracp(i)), 0._r8 )

         if ( present(is_strat_cloudborne) ) then

            if ( is_strat_cloudborne ) then

               ! convective scavenging

               conv_scav_ic(i) = 0._r8

               conv_scav_bc(i) = 0._r8

               ! stratiform scavenging
               if (convproc_do_evaprain_atonce .and. present(bergso_in)) then
                  fracp(i) = (precs(i,k)-bergso_in(i,k))*deltat / &
                       max( 1.e-12_r8, cwat(i,k) + precs(i,k)*deltat )
               else
                  fracp(i) = precs(i,k)*deltat / &
                       max( 1.e-12_r8, cwat(i,k) + precs(i,k)*deltat )
               endif

               fracp(i) = max( 0._r8, min(1._r8, fracp(i)) )

               st_scav_ic(i) = sol_facti(i,k) *fracp(i)*tracer(i,k)*rdeltat

               st_scav_bc(i) = 0._r8

            else

               ! convective scavenging

               trac_qqcw(i) = min(qqcw(i,k), &
                                  tracer(i,k)*( clds(i)/max( 0.01_r8, 1._r8-clds(i) ) ) )

               tracer_incu(i) = f_act_conv(i,k)*(tracer(i,k) + trac_qqcw(i))

               conv_scav_ic(i) = sol_factic(i,k)*cldc(i,k)*fracp(i)*tracer_incu(i)*rdeltat

               tracer_mean(i) = tracer(i,k)*(1._r8 - cldc(i,k)*f_act_conv(i,k)) - &
                                cldc(i,k)*f_act_conv(i,k)*trac_qqcw(i)
               tracer_mean(i) = max(0._r8,tracer_mean(i))

               odds(i) = precabc(i)/max(cldvcu(i,k),1.e-5_r8)*scavcoef(i,k)*deltat
               odds(i) = max(min(1._r8,odds(i)),0._r8)
               conv_scav_bc(i) = sol_factb(i,k) *cldvcu(i,k)*odds(i)*tracer_mean(i)*rdeltat


               ! stratiform scavenging

               st_scav_ic(i) = 0._r8

               odds(i) = precabs(i)/max(cldvst(i,k),1.e-5_r8)*scavcoef(i,k)*deltat
               odds(i) = max(min(1._r8,odds(i)),0._r8)
               st_scav_bc(i) = sol_factb(i,k) *cldvst(i,k)*odds(i)*tracer_mean(i)*rdeltat

            end if

         else

            ! convective scavenging

            conv_scav_ic(i) = sol_factic(i,k)*cldc(i,k)*fracp(i)*tracer(i,k)*rdeltat

            odds(i) = precabc(i)/max(cldvcu(i,k), 1.e-5_r8)*scavcoef(i,k)*deltat
            odds(i) = max( min(1._r8, odds(i)), 0._r8)
            conv_scav_bc(i) = sol_factb(i,k)*cldvcu(i,k)*odds(i)*tracer(i,k)*rdeltat

            ! stratiform scavenging

            ! fracp is the fraction of cloud water converted to precip
            ! NB: In below formula for fracp cwat is a LWC/IWC that has already
            !     precipitated out, i.e., cwat does not contain precipitation
            fracp(i) = precs(i,k)*deltat / &
                       max( 1.e-12_r8, cwat(i,k) + precs(i,k)*deltat )
            fracp(i) = max( 0._r8, min( 1._r8, fracp(i) ) )

            ! assume the corresponding amnt of tracer is removed
            st_scav_ic(i) = sol_facti(i,k)*clds(i)*fracp(i)*tracer(i,k)*rdeltat

            odds(i) = precabs(i)/max(cldvst(i,k),1.e-5_r8)*scavcoef(i,k)*deltat
            odds(i) = max(min(1._r8,odds(i)),0._r8)
            st_scav_bc(i) =sol_factb(i,k)*(cldvst(i,k)*odds(i)) *tracer(i,k)*rdeltat

         end if

         ! total convective scavenging
         srcc(i) = conv_scav_ic(i) + conv_scav_bc(i)
         finc(i) = conv_scav_ic(i)/(srcc(i) + 1.e-36_r8)

         ! total stratiform scavenging
         srcs(i) = st_scav_ic(i) + st_scav_bc(i)
         fins(i) = st_scav_ic(i)/(srcs(i) + 1.e-36_r8)

         ! make sure we dont take out more than is there
         ! ratio of amount available to amount removed
         rat(i) = tracer(i,k)/max(deltat*(srcc(i)+srcs(i)),1.e-36_r8)
         if (rat(i)<1._r8) then
            srcs(i) = srcs(i)*rat(i)
            srcc(i) = srcc(i)*rat(i)
         endif
         srct(i) = (srcc(i)+srcs(i))*omsm


         ! fraction that is not removed within the cloud
         ! (assumed to be interstitial, and subject to convective transport)
         fracp(i) = deltat*srct(i)/max(cldvst(i,k)*tracer(i,k),1.e-36_r8)  ! amount removed
         fracp(i) = max(0._r8,min(1._r8,fracp(i)))
         fracis(i,k) = 1._r8 - fracp(i)

         ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
         ! Sungsu added cumulus contribution in the below 3 blocks
         scavt(i,k) = -srct(i) + (fracev(i)*scavab(i)+fracev_cu(i)*scavabc(i))*rpdog(i)
         iscavt(i,k) = -(srcc(i)*finc(i) + srcs(i)*fins(i))*omsm

         if ( present(icscavt) ) icscavt(i,k) = -(srcc(i)*finc(i)) * omsm
         if ( present(isscavt) ) isscavt(i,k) = -(srcs(i)*fins(i)) * omsm

         if (.not. out_resuspension) then
            if (present(bcscavt)) bcscavt(i,k) = -(srcc(i) * (1-finc(i))) * omsm +  &
                                                  fracev_cu(i)*scavabc(i)*rpdog(i)

            if (present(bsscavt)) bsscavt(i,k) = -(srcs(i) * (1-fins(i))) * omsm +  &
                                                  fracev(i)*scavab(i)*rpdog(i)
         else
            bcscavt(i,k) = -(srcc(i) * (1-finc(i))) * omsm
            rcscavt(i,k) =  fracev_cu(i)*scavabc(i)*rpdog(i)

            bsscavt(i,k) = -(srcs(i) * (1-fins(i))) * omsm
            rsscavt(i,k) =  fracev(i)*scavab(i)*rpdog(i)
         end if

         dblchek(i) = tracer(i,k) + deltat*scavt(i,k)

         ! now keep track of scavenged mass and precip
         scavab(i) = scavab(i)*(1-fracev(i)) + srcs(i)*pdog(i)
         precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdog(i)
         scavabc(i) = scavabc(i)*(1-fracev_cu(i)) + srcc(i)*pdog(i)
         precabc(i) = precabc(i) + (cmfdqr(i,k) - evapc(i,k))*pdog(i)

      end do ! End of i = 1, ncol

#ifdef DEBUG
      ! only check in debug mode which aborts when larger negative values are found
      found = .false.
      do i = 1,ncol
         ! catch the larger negative values, ignore insignificant small negaive values
         if (dblchek(i) < -1.e-10_r8) then
            found = .true.
         endif
      end do

      if (found) then
         errmsg = 'wetdapa_v2: negative values found'
         errflg = 1
         return
      end if
#endif

   end do ! End of k = 1, pver

end subroutine wetdepa_v2


!==============================================================================

! This is the frozen CAM4 version of wetdepa.


   subroutine wetdepa_v1( t, p, q, pdel, &
                       cldt, cldc, cmfdqr, conicw, precs, conds, &
                       evaps, cwat, tracer, deltat, &
                       scavt, iscavt, cldv, fracis, sol_fact, ncol, &
                       scavcoef, tmelt, gravit, pver, errmsg, errflg, &
                       icscavt, isscavt, bcscavt, bsscavt, &
                       sol_facti_in, sol_factbi_in, sol_factii_in, &
                       sol_factic_in, sol_factiic_in )

      !-----------------------------------------------------------------------
      ! Purpose:
      ! scavenging code for very soluble aerosols
      !
      ! Author: P. Rasch
      ! Modified by T. Bond 3/2003 to track different removals
      !-----------------------------------------------------------------------

      implicit none

      real(r8), intent(in) ::&
         t(:,:),        &! temperature
         p(:,:),        &! pressure
         q(:,:),        &! moisture
         pdel(:,:),     &! pressure thikness
         cldt(:,:),     &! total cloud fraction
         cldc(:,:),     &! convective cloud fraction
         cmfdqr(:,:),   &! rate of production of convective precip
         conicw(:,:),   &! convective cloud water
         cwat(:,:),     &! cloud water amount
         precs(:,:),    &! rate of production of stratiform precip
         conds(:,:),    &! rate of production of condensate
         evaps(:,:),    &! rate of evaporation of precip
         cldv(:,:),     &! total cloud fraction
         deltat,               &! time step
         tracer(:,:)     ! trace species
      ! If subroutine is called with just sol_fact:
            ! sol_fact is used for both in- and below-cloud scavenging
      ! If subroutine is called with optional argument sol_facti_in:
            ! sol_fact  is used for below cloud scavenging
            ! sol_facti is used for in cloud scavenging
         real(r8), intent(in) :: sol_fact ! solubility factor (fraction of aer scavenged below & in, or just below or sol_facti_in is provided)
         real(r8), intent(in), optional :: sol_facti_in   ! solubility factor (frac of aerosol scavenged in cloud)
         real(r8), intent(in), optional :: sol_factbi_in  ! solubility factor (frac of aerosol scavenged below cloud by ice)
         real(r8), intent(in), optional :: sol_factii_in  ! solubility factor (frac of aerosol scavenged in cloud by ice)
         real(r8), intent(in), optional :: sol_factic_in(:,:)  ! sol_facti_in for convective clouds
         real(r8), intent(in), optional :: sol_factiic_in ! sol_factii_in for convective clouds
         real(r8), intent(in) :: scavcoef(:,:) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)

      integer, intent(in) :: ncol
      integer, intent(in) :: pver
      real(r8), intent(in) :: tmelt        ! freezing point of water (K)
      real(r8), intent(in) :: gravit       ! gravitational acceleration (m/s2)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(r8), intent(out) ::&
         scavt(:,:),    &! scavenging tend
         iscavt(:,:),   &! incloud scavenging tends
         fracis(:,:)     ! fraction of species not scavenged

      real(r8), intent(out), optional ::    icscavt(:,:)     ! incloud, convective
      real(r8), intent(out), optional ::    isscavt(:,:)     ! incloud, stratiform
      real(r8), intent(out), optional ::    bcscavt(:,:)     ! below cloud, convective
      real(r8), intent(out), optional ::    bsscavt(:,:)     ! below cloud, stratiform

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cwatc                ! local convective total water amount
      real(r8) cwats                ! local stratiform total water amount
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) fracev(ncol)        ! fraction of precip from above that is evaporating
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/gravit)
      real(r8) precabc(ncol)       ! conv precip from above (work array)
      real(r8) precabs(ncol)       ! strat precip from above (work array)
      real(r8) precbl               ! precip falling out of level (work array)
      real(r8) precmin              ! minimum convective precip causing scavenging
      real(r8) rat(ncol)           ! ratio of amount available to amount removed
      real(r8) scavab(ncol)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(ncol)       ! scavenged tracer flux from above (work array)
      real(r8) srcc                 ! tend for convective rain
      real(r8) srcs                 ! tend for stratiform rain
      real(r8) srct(ncol)          ! work variable
      real(r8) tracab(ncol)        ! column integrated tracer amount

      real(r8) fins                 ! fraction of rem. rate by strat rain
      real(r8) finc                 ! fraction of rem. rate by conv. rain
      real(r8) srcs1                ! work variable
      real(r8) srcs2                ! work variable
      real(r8) tc                   ! temp in celcius
      real(r8) weight               ! fraction of condensate which is ice
      real(r8) cldmabs(ncol)       ! maximum cloud at or above this level
      real(r8) cldmabc(ncol)       ! maximum cloud at or above this level
      real(r8) odds                 ! limit on removal rate (proportional to prec)
      real(r8) dblchek(ncol)
      logical :: found

      real(r8) sol_facti,  sol_factb  ! in cloud and below cloud fraction of aerosol scavenged
      real(r8) sol_factii, sol_factbi ! in cloud and below cloud fraction of aerosol scavenged by ice
      real(r8) sol_factic(ncol,pver)             ! sol_facti for convective clouds
      real(r8) sol_factiic            ! sol_factii for convective clouds
      ! sol_factic & solfact_iic added for MODAL_AERO.
      ! For stratiform cloud, cloudborne aerosol is treated explicitly,
      !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
      ! For convective cloud, cloudborne aerosol is not treated explicitly,
      !    and sol_factic is 1.0 for both cloudborne and interstitial.

      ! ------------------------------------------------------------------------
      errmsg = ''
      errflg = 0

      precmin =  0.1_r8/8.64e4_r8      ! set critical value to 0.1 mm/day in kg/m2/s

      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! default (if other sol_facts aren't in call, set all to required sol_fact
      sol_facti = sol_fact
      sol_factb = sol_fact
      sol_factii = sol_fact
      sol_factbi = sol_fact

      if ( present(sol_facti_in) )  sol_facti = sol_facti_in
      if ( present(sol_factii_in) )  sol_factii = sol_factii_in
      if ( present(sol_factbi_in) )  sol_factbi = sol_factbi_in

      sol_factic  = sol_facti
      sol_factiic = sol_factii
      if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in
      if ( present(sol_factiic_in) )  sol_factiic = sol_factiic_in

      ! this section of code is for highly soluble aerosols,
      ! the assumption is that within the cloud that
      ! all the tracer is in the cloud water
      !
      ! for both convective and stratiform clouds,
      ! the fraction of cloud water converted to precip defines
      ! the amount of tracer which is pulled out.
      !

      do i = 1,ncol
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0
      end do

      do k = 1,pver
         do i = 1,ncol
            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice
            weight = 0._r8                                 ! assume no ice

            pdog = pdel(i,k)/gravit

            ! ****************** Evaporation **************************
            ! calculate the fraction of strat precip from above
            !                 which evaporates within this layer
            fracev(i) = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))

            ! trap to ensure reasonable ratio bounds
            fracev(i) = max(0._r8,min(1._r8,fracev(i)))

            ! ****************** Convection ***************************
            ! now do the convective scavenging

            ! set odds proportional to fraction of the grid box that is swept by the
            ! precipitation =precabc/rhoh20*(area of sphere projected on plane
            !                                /volume of sphere)*deltat
            ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
            ! unless the fraction of the area that is cloud is less than odds, in which
            ! case use the cloud fraction (assumes precabs is in kg/m2/s)
            ! is really: precabs*3/4/1000./1e-3*deltat
            ! here I use .1 from Balkanski
            !
            ! use a local rate of convective rain production for incloud scav
            !odds=max(min(1._r8, &
            !     cmfdqr(i,k)*pdel(i,k)/gravit*0.1_r8*deltat),0._r8)
            !++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            !            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1.-weight) &
            ! srcs1 =  cldv(i,k)*odds*tracer(i,k)*(1.-weight) &
            !srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1.-weight) &
            !         /deltat

            ! fraction of convective cloud water converted to rain
            fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,conicw(i,k))
            ! note cmfdrq can be negative from evap of rain, so constrain it
            fracp = max(min(1._r8,fracp),0._r8)
            ! remove that amount from within the convective area
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat ! liquid only
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)/deltat             ! any condensation
!           srcs1 = 0.
            srcs1 = sol_factic(i,k)*cldt(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat &  ! liquid
                 +  sol_factiic*cldt(i,k)*fracp*tracer(i,k)*(weight)/deltat      ! ice


            !--mcb

            ! scavenge below cloud

            !            cldmabc(i) = max(cldc(i,k),cldmabc(i))
            !            cldmabc(i) = max(cldt(i,k),cldmabc(i))
            cldmabc(i) = max(cldv(i,k),cldmabc(i))
            cldmabc(i) = cldv(i,k)

            odds=max( &
                 min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8) &
                 *scavcoef(i,k)*deltat),0._r8) ! Dana and Hales coefficient (/mm)
            srcs2 = sol_factb*cldmabc(i)*odds*tracer(i,k)*(1._r8-weight)/deltat & ! liquid
                 +  sol_factbi*cldmabc(i)*odds*tracer(i,k)*(weight)/deltat    !ice
            !Note that using the temperature-determined weight doesn't make much sense here


            srcc = srcs1 + srcs2  ! convective tend by both processes
            finc = srcs1/(srcc + 1.e-36_r8) ! fraction in-cloud

            ! ****************** Stratiform ***********************
            ! now do the stratiform scavenging

            ! incloud scavenging

            ! fracp is the fraction of cloud water converted to precip
            fracp =  precs(i,k)*deltat/max(cwat(i,k),1.e-12_r8)
            fracp = max(0._r8,min(1._r8,fracp))
!           fracp = 0.     ! for debug

            ! assume the corresponding amnt of tracer is removed
            !++mcb -- remove cldc; change cldt to cldv
            !            srcs1 = (cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat
            !            srcs1 = cldv(i,k)*fracp*tracer(i,k)/deltat &
!           srcs1 = cldt(i,k)*fracp*tracer(i,k)/deltat            ! all condensate
            srcs1 = sol_facti*cldt(i,k)*fracp*tracer(i,k)/deltat*(1._r8-weight) &  ! liquid
                 + sol_factii*cldt(i,k)*fracp*tracer(i,k)/deltat*(weight)       ! ice


            ! below cloud scavenging

!           volume undergoing below cloud scavenging
            cldmabs(i) = cldv(i,k)   ! precipitating volume
!           cldmabs(i) = cldt(i,k)   ! local cloud volume

            odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*scavcoef(i,k)*deltat
            odds = max(min(1._r8,odds),0._r8)
            srcs2 =sol_factb*(cldmabs(i)*odds) *tracer(i,k)*(1._r8-weight)/deltat & ! liquid
                 + sol_factbi*(cldmabs(i)*odds) *tracer(i,k)*(weight)/deltat       ! ice
            !Note that using the temperature-determined weight doesn't make much sense here


            srcs = srcs1 + srcs2             ! total stratiform scavenging
            fins=srcs1/(srcs + 1.e-36_r8)    ! fraction taken by incloud processes

            ! make sure we dont take out more than is there
            ! ratio of amount available to amount removed
            rat(i) = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat(i) < 1._r8) then
               srcs = srcs*rat(i)
               srcc = srcc*rat(i)
            endif
            srct(i) = (srcc+srcs)*omsm


            ! fraction that is not removed within the cloud
            ! (assumed to be interstitial, and subject to convective transport)
            fracp = deltat*srct(i)/max(cldmabs(i)*tracer(i,k),1.e-36_r8)  ! amount removed
            fracp = max(0._r8,min(1._r8,fracp))
            fracis(i,k) = 1._r8 - fracp

            ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
            scavt(i,k) = -srct(i) + fracev(i)*scavab(i)*gravit/pdel(i,k)
            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

            if ( present(icscavt) ) icscavt(i,k) = -(srcc*finc) * omsm
            if ( present(isscavt) ) isscavt(i,k) = -(srcs*fins) * omsm
            if ( present(bcscavt) ) bcscavt(i,k) = -(srcc * (1-finc)) * omsm
            if ( present(bsscavt) ) bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
                 fracev(i)*scavab(i)*gravit/pdel(i,k)

            dblchek(i) = tracer(i,k) + deltat*scavt(i,k)

            ! now keep track of scavenged mass and precip
            scavab(i) = scavab(i)*(1-fracev(i)) + srcs*pdel(i,k)/gravit
            precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdel(i,k)/gravit
            scavabc(i) = scavabc(i) + srcc*pdel(i,k)/gravit
            precabc(i) = precabc(i) + (cmfdqr(i,k))*pdel(i,k)/gravit
            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/gravit

         end do

#ifdef DEBUG
         ! only check in debug mode which aborts when larger negative values are found
         found = .false.
         do i = 1,ncol
            ! catch the larger negative values, ignore insignificant small negaive values
            if (dblchek(i) < -1.e-10_r8) then
               found = .true.
            endif
         end do

         if (found) then
            errmsg = 'wetdapa_v1: negative values found'
            errflg = 1
            return
         end if
#endif

      end do

   end subroutine wetdepa_v1

!##############################################################################

end module wetdep
