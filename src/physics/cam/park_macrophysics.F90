
  module park_macrophysics

  !--------------------------------------------------- !
  ! Purpose     : CAM Interface for Cloud Macrophysics !
  ! Author      : Sungsu Park                          !
  ! Description : Park et al. 2010.                    !
  ! For questions, contact Sungsu Park                 !
  !                        e-mail : sungsup@ucar.edu   !
  !                        phone  : 303-497-1375       !
  !                                                    !
  ! CCPP-compliant extraction from cldwat2m_macro.F90  !
  ! and macrop_driver.F90.                             !
  !--------------------------------------------------- !

   use ccpp_kinds,              only: kind_phys
   use wv_saturation,           only: qsat_water, findsp_vc
   use cldfrc2m,                only: astG_PDF_single, astG_PDF, astG_RHU_single, &
                                      astG_RHU, aist_single, aist_vector,         &
                                      rhmini_const, rhmaxi_const
   use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
   use ccpp_const_utils,          only: ccpp_const_get_idx

   implicit none
   private
   save

   public :: park_macrophysics_init
   public :: park_macrophysics_detrain_run
   public :: park_macrophysics_run

   ! -------------- !
   ! Set Parameters !
   ! -------------- !

   ! ------------------------------------------------------------------------------- !
   ! Parameter used for selecting generalized critical RH for liquid and ice stratus !
   ! ------------------------------------------------------------------------------- !

   integer :: i_rhminl ! This is for liquid stratus fraction.
                       ! If 0 : Original fixed critical RH from the namelist.
                       ! If 1 : Add convective detrainment effect on the above '0' option.
                       !        In this case, 'tau_detw' [s] should be specified below.
                       ! If 2 : Use fully scale-adaptive method.
                       !        In this case, 'tau_detw' [s] and 'c_aniso' [no unit] should
                       !        be specified below.

   integer :: i_rhmini ! This is for ice stratus fraction.
                       ! If 0 : Original fixed critical RH from the namelist.
                       ! If 1 : Add convective detrainment effect on the above '0' option.
                       !        In this case, 'tau_deti' [s] should be specified below.
                       ! If 2 : Use fully scale-adaptive method.
                       !        In this case, 'tau_deti' [s] and 'c_aniso' [no unit] should
                       !        be specified below.
                       ! Note that 'micro_mg_cam' is using below 'rhmini_const', regardless
                       ! of 'i_rhmini'.  This connection should be built in future.

   real(kind_phys), parameter :: tau_detw =100._kind_phys   ! Dissipation time scale of convective liquid condensate detrained
                                                            !  into the clear portion. [hr]. 0.5-3 hr is possible.
   real(kind_phys), parameter :: tau_deti =  1._kind_phys   ! Dissipation time scale of convective ice    condensate detrained
                                                            !  into the clear portion. [hr]. 0.5-3 hr is possible.
   real(kind_phys), parameter :: c_aniso  =  1._kind_phys   ! Inverse of anisotropic factor of PBL turbulence

   ! ----------------------------- !
   ! Parameters for Liquid Stratus !
   ! ----------------------------- !

   logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
                                                     ! use Slingo (triangular PDF-based) liquid stratus fraction
   real(kind_phys), parameter  :: qlst_min     = 2.e-5_kind_phys   ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(kind_phys), parameter  :: qlst_max     = 3.e-3_kind_phys   ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(kind_phys), parameter  :: cc           = 0.1_kind_phys     ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
   integer,  parameter  :: niter        = 2          ! For iterative computation of QQ with 'ramda' below.
   real(kind_phys), parameter  :: ramda        = 0.5_kind_phys     ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
   real(kind_phys), private    :: rhminl_const              ! Critical RH for low-level  liquid stratus clouds
   real(kind_phys), private    :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
   real(kind_phys), private    :: rhminh_const              ! Critical RH for high-level liquid stratus clouds
   real(kind_phys), private    :: premit                    ! Top    height for mid-level liquid stratus fraction
   real(kind_phys), private    :: premib                    ! Bottom height for mid-level liquid stratus fraction

   real(kind_phys), parameter :: qsmall = 1.e-18_kind_phys         ! Smallest mixing ratio considered in the macrophysics

   ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus
   !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus,
   !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
   !               liquid condensate, not cumulus ice condensate.
   logical, parameter :: cu_det_st  = .false.

   ! Module-level iulog for logging
   integer :: iulog

   ! Module-level physical constants (set at the start of park_macrophysics_run)
   real(kind_phys) :: latvap
   real(kind_phys) :: cpair
   real(kind_phys) :: latice

   ! Module-level top_lev (set at the start of park_macrophysics_run)
   integer :: top_lev_mod

   contains

   ! -------------- !
   ! Initialization !
   ! -------------- !

   subroutine park_macrophysics_init(rhminl_opt_in, rhmini_opt_in, &
        rhminl_in, rhminl_adj_land_in, rhminh_in, premit_in, premib_in, &
        iulog_in, masterproc_in, errmsg, errflg)

   !--------------------------------------------------------------------- !
   !                                                                      !
   ! Purpose: Initialize constants for the liquid stratiform macrophysics !
   !                                                                      !
   ! Author:  Sungsu Park, Dec.01.2009.                                   !
   !                                                                      !
   !--------------------------------------------------------------------- !

   integer,  intent(in) :: rhminl_opt_in
   integer,  intent(in) :: rhmini_opt_in
   real(kind_phys), intent(in) :: rhminl_in
   real(kind_phys), intent(in) :: rhminl_adj_land_in
   real(kind_phys), intent(in) :: rhminh_in
   real(kind_phys), intent(in) :: premit_in
   real(kind_phys), intent(in) :: premib_in
   integer,  intent(in) :: iulog_in
   logical,  intent(in) :: masterproc_in
   character(len=512), intent(out) :: errmsg
   integer,  intent(out) :: errflg

   errmsg = ''
   errflg = 0

   i_rhminl   = rhminl_opt_in
   i_rhmini   = rhmini_opt_in

   rhminl_const          = rhminl_in
   rhminl_adj_land_const = rhminl_adj_land_in
   rhminh_const          = rhminh_in
   premit                = premit_in
   premib                = premib_in

   iulog = iulog_in

   if( masterproc_in ) then
       write(iulog,*) 'Park Macrophysics Parameters'
       write(iulog,*) '  rhminl          = ', rhminl_const
       write(iulog,*) '  rhminl_adj_land = ', rhminl_adj_land_const
       write(iulog,*) '  rhminh          = ', rhminh_const
       write(iulog,*) '  premit          = ', premit
       write(iulog,*) '  premib          = ', premib
       write(iulog,*) '  i_rhminl        = ', i_rhminl
       write(iulog,*) '  i_rhmini        = ', i_rhmini
   end if

   end subroutine park_macrophysics_init

   ! ----------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment   !
   ! or stratiform cloud                                         !
   ! ----------------------------------------------------------- !

   subroutine park_macrophysics_detrain_run( &
       ncol, pver, top_lev, &
       latice, cpair, gravit, &
       do_detrain, &
       t, pdel, dlf, dlf2, &
       tend_cldliq, tend_cldice, tend_numliq, tend_numice, tend_s, &
       det_s, det_ice, &
       dpdlfliq, dpdlfice, shdlfliq, shdlfice, dpdlft, shdlft, &
       dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni, &
       errmsg, errflg)

   ! Input arguments
   integer,  intent(in)  :: ncol
   integer,  intent(in)  :: pver
   integer,  intent(in)  :: top_lev
   real(kind_phys), intent(in)  :: latice
   real(kind_phys), intent(in)  :: cpair
   real(kind_phys), intent(in)  :: gravit
   logical,  intent(in)  :: do_detrain
   real(kind_phys), intent(in)  :: t(:,:)
   real(kind_phys), intent(in)  :: pdel(:,:)
   real(kind_phys), intent(in)  :: dlf(:,:)
   real(kind_phys), intent(in)  :: dlf2(:,:)

   ! Output arguments
   real(kind_phys), intent(out) :: tend_cldliq(:,:)
   real(kind_phys), intent(out) :: tend_cldice(:,:)
   real(kind_phys), intent(out) :: tend_numliq(:,:)
   real(kind_phys), intent(out) :: tend_numice(:,:)
   real(kind_phys), intent(out) :: tend_s(:,:)
   real(kind_phys), intent(out) :: det_s(:)
   real(kind_phys), intent(out) :: det_ice(:)
   real(kind_phys), intent(out) :: dpdlfliq(:,:)
   real(kind_phys), intent(out) :: dpdlfice(:,:)
   real(kind_phys), intent(out) :: shdlfliq(:,:)
   real(kind_phys), intent(out) :: shdlfice(:,:)
   real(kind_phys), intent(out) :: dpdlft(:,:)
   real(kind_phys), intent(out) :: shdlft(:,:)
   real(kind_phys), intent(out) :: dlf_T(:,:)
   real(kind_phys), intent(out) :: dlf_qv(:,:)
   real(kind_phys), intent(out) :: dlf_ql(:,:)
   real(kind_phys), intent(out) :: dlf_qi(:,:)
   real(kind_phys), intent(out) :: dlf_nl(:,:)
   real(kind_phys), intent(out) :: dlf_ni(:,:)
   character(len=512), intent(out) :: errmsg
   integer,  intent(out) :: errflg

   ! Local variables
   integer  :: i, k
   real(kind_phys) :: dum1

   errmsg = ''
   errflg = 0

   ! Initialize all outputs to zero
   tend_cldliq(:,:) = 0._kind_phys
   tend_cldice(:,:) = 0._kind_phys
   tend_numliq(:,:) = 0._kind_phys
   tend_numice(:,:) = 0._kind_phys
   tend_s(:,:)      = 0._kind_phys

   det_s(:)   = 0._kind_phys
   det_ice(:) = 0._kind_phys

   dpdlfliq = 0._kind_phys
   dpdlfice = 0._kind_phys
   shdlfliq = 0._kind_phys
   shdlfice = 0._kind_phys
   dpdlft   = 0._kind_phys
   shdlft   = 0._kind_phys

   dlf_T(:,:)  = 0._kind_phys
   dlf_qv(:,:) = 0._kind_phys
   dlf_ql(:,:) = 0._kind_phys
   dlf_qi(:,:) = 0._kind_phys
   dlf_nl(:,:) = 0._kind_phys
   dlf_ni(:,:) = 0._kind_phys

   do k = top_lev, pver
   do i = 1, ncol
      if( t(i,k) > 268.15_kind_phys ) then
          dum1 = 0.0_kind_phys
      elseif( t(i,k) < 238.15_kind_phys ) then
          dum1 = 1.0_kind_phys
      else
          dum1 = ( 268.15_kind_phys - t(i,k) ) / 30._kind_phys
      endif

     ! If detrainment was done elsewhere, still update the variables used for output
     ! assuming that the temperature split between liquid and ice is the same as assumed
     ! here.
     if (do_detrain) then
        tend_cldliq(i,k) = dlf(i,k) * ( 1._kind_phys - dum1 )
        tend_cldice(i,k) = dlf(i,k) * dum1
      ! dum2              = dlf(i,k) * ( 1._kind_phys - dum1 )
        tend_numliq(i,k) = 3._kind_phys * ( max(0._kind_phys, ( dlf(i,k) - dlf2(i,k) )) * ( 1._kind_phys - dum1 ) ) / &
             (4._kind_phys*3.14_kind_phys* 8.e-6_kind_phys**3*997._kind_phys) + & ! Deep    Convection
             3._kind_phys * (                         dlf2(i,k)    * ( 1._kind_phys - dum1 ) ) / &
             (4._kind_phys*3.14_kind_phys*10.e-6_kind_phys**3*997._kind_phys)     ! Shallow Convection
      ! dum2              = dlf(i,k) * dum1
        tend_numice(i,k) = 3._kind_phys * ( max(0._kind_phys, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) / &
             (4._kind_phys*3.14_kind_phys*25.e-6_kind_phys**3*500._kind_phys) + & ! Deep    Convection
             3._kind_phys * (                         dlf2(i,k)    *  dum1 ) / &
             (4._kind_phys*3.14_kind_phys*50.e-6_kind_phys**3*500._kind_phys)     ! Shallow Convection
        tend_s(i,k)      = dlf(i,k) * dum1 * latice
       else
          tend_cldliq(i,k) = 0._kind_phys
          tend_cldice(i,k) = 0._kind_phys
          tend_numliq(i,k) = 0._kind_phys
          tend_numice(i,k) = 0._kind_phys
          tend_s(i,k)      = 0._kind_phys
     end if

     ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
     !   track of the integrals of ice and static energy that is effected from conversion to ice
     !   so that the energy checker doesn't complain.
       det_s(i)                  = det_s(i) + tend_s(i,k)*pdel(i,k)/gravit
       det_ice(i)                = det_ice(i) - tend_cldice(i,k)*pdel(i,k)/gravit

     ! Targetted detrainment of convective liquid water either directly into the
     ! existing liquid stratus or into the environment.
      if( cu_det_st ) then
          dlf_T(i,k)  = tend_s(i,k)/cpair
          dlf_qv(i,k) = 0._kind_phys
          dlf_ql(i,k) = tend_cldliq(i,k)
          dlf_qi(i,k) = tend_cldice(i,k)
          dlf_nl(i,k) = tend_numliq(i,k)
          dlf_ni(i,k) = tend_numice(i,k)
          tend_cldliq(i,k) = 0._kind_phys
          tend_cldice(i,k) = 0._kind_phys
          tend_numliq(i,k) = 0._kind_phys
          tend_numice(i,k) = 0._kind_phys
          tend_s(i,k)      = 0._kind_phys
          dpdlfliq(i,k)             = 0._kind_phys
          dpdlfice(i,k)             = 0._kind_phys
          shdlfliq(i,k)             = 0._kind_phys
          shdlfice(i,k)             = 0._kind_phys
          dpdlft  (i,k)             = 0._kind_phys
          shdlft  (i,k)             = 0._kind_phys
       else
          dpdlfliq(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( 1._kind_phys - dum1 )
          dpdlfice(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( dum1 )
          dpdlft  (i,k) = ( dlf(i,k) - dlf2(i,k) ) * dum1 * latice/cpair

          shdlfliq(i,k) = dlf2(i,k) * ( 1._kind_phys - dum1 )
          shdlfice(i,k) = dlf2(i,k) * ( dum1 )
          shdlft  (i,k) = dlf2(i,k) * dum1 * latice/cpair
      endif
   end do
   end do

   det_ice(:ncol) = det_ice(:ncol)/1000._kind_phys  ! divide by density of water

   end subroutine park_macrophysics_detrain_run

   ! ------------------------------ !
   ! Stratiform Liquid Macrophysics !
   ! ------------------------------ !

   ! In the version, 'macro --> micro --> advective forcing --> macro...'
   ! A_...: only 'advective forcing' without 'microphysical tendency'
   ! C_...: only 'microphysical tendency'
   ! D_...: only 'detrainment of cumulus condensate'
   ! So, 'A' and 'C' are exclusive.

   subroutine park_macrophysics_run( &
       ncol, pver, top_lev, nstep, &
       dtime, &
       cpair_in, latvap_in, latice_in, rh2o, gravit, rair, &
       const_props, qmin, &
       t, q_wv, cldliq, cldice, numliq, numice, &
       pmid, pdel, &
       tcwat, qcwat, lcwat, iccwat, nlwat, niwat, &
       CC_T, CC_qv, CC_ql, CC_qi, CC_nl, CC_ni, CC_qlst, &
       dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni, &
       concld_old, concld, clrw_old, clri_old, &
       landfrac, snowh, &
       do_cldice, &
       tlat, qvlat, qcten, qiten, ncten, niten, &
       cmeliq, qvadj, qladj, qiadj, qllim, qilim, &
       cld, alst, aist, qlst, qist, ast, &
       rhmin_liq, rhmin_ice, &
       errmsg, errflg)

   ! Input arguments
   integer,  intent(in)    :: ncol
   integer,  intent(in)    :: pver
   integer,  intent(in)    :: top_lev
   integer,  intent(in)    :: nstep
   real(kind_phys), intent(in) :: dtime
   real(kind_phys), intent(in) :: cpair_in
   real(kind_phys), intent(in) :: latvap_in
   real(kind_phys), intent(in) :: latice_in
   real(kind_phys), intent(in) :: rh2o
   real(kind_phys), intent(in) :: gravit
   real(kind_phys), intent(in) :: rair

   type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
   real(kind_phys), intent(in) :: qmin(:)

   real(kind_phys), intent(in)    :: t(:,:)
   real(kind_phys), intent(in)    :: q_wv(:,:)
   real(kind_phys), intent(in)    :: cldliq(:,:)
   real(kind_phys), intent(in)    :: cldice(:,:)
   real(kind_phys), intent(in)    :: numliq(:,:)
   real(kind_phys), intent(in)    :: numice(:,:)
   real(kind_phys), intent(in)    :: pmid(:,:)
   real(kind_phys), intent(in)    :: pdel(:,:)

   ! Input-output variables (equilibrium state storage)
   real(kind_phys), intent(inout) :: tcwat(:,:)
   real(kind_phys), intent(inout) :: qcwat(:,:)
   real(kind_phys), intent(inout) :: lcwat(:,:)
   real(kind_phys), intent(inout) :: iccwat(:,:)
   real(kind_phys), intent(inout) :: nlwat(:,:)
   real(kind_phys), intent(inout) :: niwat(:,:)

   real(kind_phys), intent(in)    :: CC_T(:,:)
   real(kind_phys), intent(in)    :: CC_qv(:,:)
   real(kind_phys), intent(in)    :: CC_ql(:,:)
   real(kind_phys), intent(in)    :: CC_qi(:,:)
   real(kind_phys), intent(in)    :: CC_nl(:,:)
   real(kind_phys), intent(in)    :: CC_ni(:,:)
   real(kind_phys), intent(in)    :: CC_qlst(:,:)

   real(kind_phys), intent(in)    :: dlf_T(:,:)
   real(kind_phys), intent(in)    :: dlf_qv(:,:)
   real(kind_phys), intent(in)    :: dlf_ql(:,:)
   real(kind_phys), intent(in)    :: dlf_qi(:,:)
   real(kind_phys), intent(in)    :: dlf_nl(:,:)
   real(kind_phys), intent(in)    :: dlf_ni(:,:)

   real(kind_phys), intent(in)    :: concld_old(:,:)
   real(kind_phys), intent(in)    :: concld(:,:)
   real(kind_phys), intent(in)    :: clrw_old(:,:)
   real(kind_phys), intent(in)    :: clri_old(:,:)

   real(kind_phys), intent(in)    :: landfrac(:)
   real(kind_phys), intent(in)    :: snowh(:)
   logical,  intent(in)           :: do_cldice

   ! Output variables
   real(kind_phys), intent(out)   :: tlat(:,:)
   real(kind_phys), intent(out)   :: qvlat(:,:)
   real(kind_phys), intent(out)   :: qcten(:,:)
   real(kind_phys), intent(out)   :: qiten(:,:)
   real(kind_phys), intent(out)   :: ncten(:,:)
   real(kind_phys), intent(out)   :: niten(:,:)
   real(kind_phys), intent(out)   :: cmeliq(:,:)
   real(kind_phys), intent(out)   :: qvadj(:,:)
   real(kind_phys), intent(out)   :: qladj(:,:)
   real(kind_phys), intent(out)   :: qiadj(:,:)
   real(kind_phys), intent(out)   :: qllim(:,:)
   real(kind_phys), intent(out)   :: qilim(:,:)
   real(kind_phys), intent(out)   :: cld(:,:)
   real(kind_phys), intent(out)   :: alst(:,:)
   real(kind_phys), intent(out)   :: aist(:,:)
   real(kind_phys), intent(out)   :: qlst(:,:)
   real(kind_phys), intent(out)   :: qist(:,:)
   real(kind_phys), intent(out)   :: ast(:,:)
   real(kind_phys), intent(out)   :: rhmin_liq(:,:)
   real(kind_phys), intent(out)   :: rhmin_ice(:,:)
   character(len=512), intent(out) :: errmsg
   integer,  intent(out) :: errflg

   ! Local variables for advective tendency computation (from macrop_driver)
   real(kind_phys) :: rdtime
   real(kind_phys) :: qc(ncol,pver)
   real(kind_phys) :: qi_loc(ncol,pver)
   real(kind_phys) :: nc(ncol,pver)
   real(kind_phys) :: ni_loc(ncol,pver)
   real(kind_phys) :: ttend(ncol,pver)
   real(kind_phys) :: qtend(ncol,pver)
   real(kind_phys) :: ltend(ncol,pver)
   real(kind_phys) :: itend(ncol,pver)
   real(kind_phys) :: lmitend(ncol,pver)
   real(kind_phys) :: nltend(ncol,pver)
   real(kind_phys) :: nitend(ncol,pver)
   real(kind_phys) :: t_inout(ncol,pver)
   real(kind_phys) :: qv_inout(ncol,pver)
   real(kind_phys) :: ql_inout(ncol,pver)
   real(kind_phys) :: qi_inout(ncol,pver)
   real(kind_phys) :: nl_inout(ncol,pver)
   real(kind_phys) :: ni_inout(ncol,pver)
   real(kind_phys) :: zeros(ncol,pver)
   real(kind_phys) :: CC_T_loc(ncol,pver)
   real(kind_phys) :: CC_qv_loc(ncol,pver)
   real(kind_phys) :: CC_ql_loc(ncol,pver)
   real(kind_phys) :: CC_qi_loc(ncol,pver)
   real(kind_phys) :: CC_nl_loc(ncol,pver)
   real(kind_phys) :: CC_ni_loc(ncol,pver)
   real(kind_phys) :: CC_qlst_loc(ncol,pver)

   ! Local variables from mmacro_pcond
   integer :: ixcldliq, ixcldice, ixq
   integer :: icol
   integer :: i, j, k, iter, ii, jj                        ! Loop indexes

   ! Thermodynamic state variables

   real(kind_phys) T_loc(ncol,pver)                        ! Temperature of equilibrium reference state
                                                           ! from which 'Micro & Macro' are computed [K]
   real(kind_phys) T1(ncol,pver)                           ! Temperature after 'fice_force' on T01
   real(kind_phys) T_0(ncol,pver)                          ! Temperature after 'instratus_condensate' on T1
   real(kind_phys) T_05(ncol,pver)                         ! Temperature after 'advection' on T_0
   real(kind_phys) T_prime0(ncol,pver)                     ! Temperature after 'Macrophysics (QQ)' on T_05star
   real(kind_phys) T_dprime(ncol,pver)                     ! Temperature after 'fice_force' on T_prime
   real(kind_phys) T_star(ncol,pver)                       ! Temperature after 'instratus_condensate' on T_dprime

   real(kind_phys) qv(ncol,pver)                           ! Grid-mean qv of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(kind_phys) qv1(ncol,pver)                          ! Grid-mean qv after 'fice_force' on qv01
   real(kind_phys) qv_0(ncol,pver)                         ! Grid-mean qv after 'instratus_condensate' on qv1
   real(kind_phys) qv_05(ncol,pver)                        ! Grid-mean qv after 'advection' on qv_0
   real(kind_phys) qv_prime0(ncol,pver)                    ! Grid-mean qv after 'Macrophysics (QQ)' on qv_05star
   real(kind_phys) qv_dprime(ncol,pver)                    ! Grid-mean qv after 'fice_force' on qv_prime
   real(kind_phys) qv_star(ncol,pver)                      ! Grid-mean qv after 'instratus_condensate' on qv_dprime

   real(kind_phys) ql(ncol,pver)                           ! Grid-mean ql of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(kind_phys) ql1(ncol,pver)                          ! Grid-mean ql after 'fice_force' on ql01
   real(kind_phys) ql_0(ncol,pver)                         ! Grid-mean ql after 'instratus_condensate' on ql1
   real(kind_phys) ql_05(ncol,pver)                        ! Grid-mean ql after 'advection' on ql_0
   real(kind_phys) ql_prime0(ncol,pver)                    ! Grid-mean ql after 'Macrophysics (QQ)' on ql_05star
   real(kind_phys) ql_dprime(ncol,pver)                    ! Grid-mean ql after 'fice_force' on ql_prime
   real(kind_phys) ql_star(ncol,pver)                      ! Grid-mean ql after 'instratus_condensate' on ql_dprime

   real(kind_phys) qi_mm(ncol,pver)                        ! Grid-mean qi of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(kind_phys) qi1(ncol,pver)                          ! Grid-mean qi after 'fice_force' on qi01
   real(kind_phys) qi_0(ncol,pver)                         ! Grid-mean qi after 'instratus_condensate' on qi1
   real(kind_phys) qi_05(ncol,pver)                        ! Grid-mean qi after 'advection' on qi_0
   real(kind_phys) qi_prime0(ncol,pver)                    ! Grid-mean qi after 'Macrophysics (QQ)' on qi_05star
   real(kind_phys) qi_dprime(ncol,pver)                    ! Grid-mean qi after 'fice_force' on qi_prime
   real(kind_phys) qi_star(ncol,pver)                      ! Grid-mean qi after 'instratus_condensate' on qi_dprime

   real(kind_phys) nl(ncol,pver)                           ! Grid-mean nl of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(kind_phys) nl1(ncol,pver)                          ! Grid-mean nl after 'fice_force' on nl01
   real(kind_phys) nl_0(ncol,pver)                         ! Grid-mean nl after 'instratus_condensate' on nl1
   real(kind_phys) nl_05(ncol,pver)                        ! Grid-mean nl after 'advection' on nl_0
   real(kind_phys) nl_prime0(ncol,pver)                    ! Grid-mean nl after 'Macrophysics (QQ)' on nl_05star
   real(kind_phys) nl_dprime(ncol,pver)                    ! Grid-mean nl after 'fice_force' on nl_prime
   real(kind_phys) nl_star(ncol,pver)                      ! Grid-mean nl after 'instratus_condensate' on nl_dprime

   real(kind_phys) ni(ncol,pver)                           ! Grid-mean ni of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(kind_phys) ni1(ncol,pver)                          ! Grid-mean ni after 'fice_force' on ni01
   real(kind_phys) ni_0(ncol,pver)                         ! Grid-mean ni after 'instratus_condensate' on ni1
   real(kind_phys) ni_05(ncol,pver)                        ! Grid-mean ni after 'advection' on ni_0
   real(kind_phys) ni_prime0(ncol,pver)                    ! Grid-mean ni after 'Macrophysics (QQ)' on ni_05star
   real(kind_phys) ni_dprime(ncol,pver)                    ! Grid-mean ni after 'fice_force' on ni_prime
   real(kind_phys) ni_star(ncol,pver)                      ! Grid-mean ni after 'instratus_condensate' on ni_dprime

   real(kind_phys) a_st(ncol,pver)                         ! Stratus fraction of equilibrium reference state
   real(kind_phys) a_st_0(ncol,pver)                       ! Stratus fraction at '_0' state
   real(kind_phys) a_st_star(ncol,pver)                    ! Stratus fraction at '_star' state

   real(kind_phys) al_st(ncol,pver)                        ! Liquid stratus fraction of equilibrium reference state
   real(kind_phys) al_st_0(ncol,pver)                      ! Liquid stratus fraction at '_0' state
   real(kind_phys) al_st_nc(ncol,pver)                     ! Non-physical liquid stratus fraction in the non-cumulus pixels

   real(kind_phys) ai_st(ncol,pver)                        ! Ice stratus fraction of equilibrium reference state
   real(kind_phys) ai_st_0(ncol,pver)                      ! Ice stratus fraction at '_0' state
   real(kind_phys) ai_st_nc(ncol,pver)                     ! Non-physical ice stratus fraction in the non-cumulus pixels

   real(kind_phys) ql_st(ncol,pver)                        ! In-stratus LWC of equilibrium reference state [kg/kg]
   real(kind_phys) ql_st_0(ncol,pver)                      ! In-stratus LWC at '_0' state

   real(kind_phys) qi_st(ncol,pver)                        ! In-stratus IWC of equilibrium reference state [kg/kg]
   real(kind_phys) qi_st_0(ncol,pver)                      ! In-stratus IWC at '_0' state

 ! Cumulus properties

   real(kind_phys) dacudt(ncol,pver)
   real(kind_phys) a_cu(ncol,pver)

 ! Adjustment tendency in association with 'positive_moisture'

   real(kind_phys) Tten_pwi1(ncol,pver)                    ! Pre-process T  tendency of input equilibrium state [K/s]
   real(kind_phys) qvten_pwi1(ncol,pver)                   ! Pre-process qv tendency of input equilibrium state [kg/kg/s]
   real(kind_phys) qlten_pwi1(ncol,pver)                   ! Pre-process ql tendency of input equilibrium state [kg/kg/s]
   real(kind_phys) qiten_pwi1(ncol,pver)                   ! Pre-process qi tendency of input equilibrium state [kg/kg/s]
   real(kind_phys) nlten_pwi1(ncol,pver)                   ! Pre-process nl tendency of input equilibrium state [#/kg/s]
   real(kind_phys) niten_pwi1(ncol,pver)                   ! Pre-process ni tendency of input equilibrium state [#/kg/s]

   real(kind_phys) Tten_pwi2(ncol,pver)                    ! Post-process T  tendency of provisional equilibrium state [K/s]
   real(kind_phys) qvten_pwi2(ncol,pver)                   ! Post-process qv tendency of provisional equilibrium state [kg/kg/s]
   real(kind_phys) qlten_pwi2(ncol,pver)                   ! Post-process ql tendency of provisional equilibrium state [kg/kg/s]
   real(kind_phys) qiten_pwi2(ncol,pver)                   ! Post-process qi tendency of provisional equilibrium state [kg/kg/s]
   real(kind_phys) nlten_pwi2(ncol,pver)                   ! Post-process nl tendency of provisoonal equilibrium state [#/kg/s]
   real(kind_phys) niten_pwi2(ncol,pver)                   ! Post-process ni tendency of provisional equilibrium state [#/kg/s]

   real(kind_phys) A_T_adj(ncol,pver)                      ! After applying external advective forcing [K/s]
   real(kind_phys) A_qv_adj(ncol,pver)                     ! After applying external advective forcing [kg/kg/s]
   real(kind_phys) A_ql_adj(ncol,pver)                     ! After applying external advective forcing [kg/kg/s]
   real(kind_phys) A_qi_adj(ncol,pver)                     ! After applying external advective forcing [kg/kg/s]
   real(kind_phys) A_nl_adj(ncol,pver)                     ! After applying external advective forcing [#/kg/s]
   real(kind_phys) A_ni_adj(ncol,pver)                     ! After applying external advective forcing [#/kg/s]

 ! Adjustment tendency in association with 'instratus_condensate'

   real(kind_phys) QQw1(ncol,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(kind_phys) QQi1(ncol,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]
   real(kind_phys) QQw2(ncol,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(kind_phys) QQi2(ncol,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]

   real(kind_phys) QQnl1(ncol,pver)          ! Tendency of nl associated with QQw1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(kind_phys) QQni1(ncol,pver)          ! Tendency of ni associated with QQi1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(kind_phys) QQnl2(ncol,pver)          ! Tendency of nl associated with QQw2 only when QQw2<0 (net evaporation) [#/kg/s]
   real(kind_phys) QQni2(ncol,pver)          ! Tendency of ni associated with QQi2 only when QQw2<0 (net evaporation) [#/kg/s]

 ! Macrophysical process tendency variables

   real(kind_phys) QQ(ncol,pver)             ! Net condensation rate into water+ice           [kg/kg/s]
   real(kind_phys) QQw(ncol,pver)            ! Net condensation rate into water               [kg/kg/s]
   real(kind_phys) QQi(ncol,pver)            ! Net condensation rate into ice                 [kg/kg/s]
   real(kind_phys) QQnl(ncol,pver)           ! Tendency of nl associated with QQw both for condensation and evaporation [#/kg/s]
   real(kind_phys) QQni(ncol,pver)           ! Tendency of ni associated with QQi both for condensation and evaporation [#/kg/s]
   real(kind_phys) ACnl(ncol,pver)           ! Cloud liquid droplet (nl) activation tendency [#/kg/s]
   real(kind_phys) ACni(ncol,pver)           ! Cloud ice    droplet (ni) activation tendency [#/kg/s]

   real(kind_phys) QQw_prev(ncol,pver)
   real(kind_phys) QQi_prev(ncol,pver)
   real(kind_phys) QQnl_prev(ncol,pver)
   real(kind_phys) QQni_prev(ncol,pver)

   real(kind_phys) QQw_prog(ncol,pver)
   real(kind_phys) QQi_prog(ncol,pver)
   real(kind_phys) QQnl_prog(ncol,pver)
   real(kind_phys) QQni_prog(ncol,pver)

   real(kind_phys) QQ_final(ncol,pver)
   real(kind_phys) QQw_final(ncol,pver)
   real(kind_phys) QQi_final(ncol,pver)
   real(kind_phys) QQn_final(ncol,pver)
   real(kind_phys) QQnl_final(ncol,pver)
   real(kind_phys) QQni_final(ncol,pver)

   real(kind_phys) QQ_all(ncol,pver)         ! QQw_all    + QQi_all
   real(kind_phys) QQw_all(ncol,pver)        ! QQw_final  + QQw1  + QQw2  + qlten_pwi1 + qlten_pwi2 + A_ql_adj [kg/kg/s]
   real(kind_phys) QQi_all(ncol,pver)        ! QQi_final  + QQi1  + QQi2  + qiten_pwi1 + qiten_pwi2 + A_qi_adj [kg/kg/s]
   real(kind_phys) QQn_all(ncol,pver)        ! QQnl_all   + QQni_all
   real(kind_phys) QQnl_all(ncol,pver)       ! QQnl_final + QQnl1 + QQnl2 + nlten_pwi1 + nlten_pwi2 + ACnl [#/kg/s]
   real(kind_phys) QQni_all(ncol,pver)       ! QQni_final + QQni1 + QQni2 + niten_pwi1 + niten_pwi2 + ACni [#/kg/s]

 ! Coefficient for computing QQ and related processes

   real(kind_phys) U(ncol,pver)                            ! Grid-mean RH
   real(kind_phys) U_nc(ncol,pver)                         ! Mean RH of non-cumulus pixels
   real(kind_phys) G_nc(ncol,pver)                         ! d(U_nc)/d(a_st_nc)
   real(kind_phys) F_nc(ncol,pver)                         ! A function of second parameter for a_st_nc
   real(kind_phys) alpha                                   ! = 1/qs
   real(kind_phys) beta                                    ! = (qv/qs**2)*dqsdT
   real(kind_phys) betast                                  ! = alpha*dqsdT
   real(kind_phys) gammal                                  ! = alpha + (latvap/cpair)*beta
   real(kind_phys) gammai                                  ! = alpha + ((latvap+latice)/cpair)*beta
   real(kind_phys) gammaQ                                  ! = alpha + (latvap/cpair)*beta
   real(kind_phys) deltal                                  ! = 1 + a_st*(latvap/cpair)*(betast/alpha)
   real(kind_phys) deltai                                  ! = 1 + a_st*((latvap+latice)/cpair)*(betast/alpha)
   real(kind_phys) A_Tc                                    ! Advective external forcing of Tc [K/s]
   real(kind_phys) A_qt                                    ! Advective external forcing of qt [kg/kg/s]
   real(kind_phys) C_Tc                                    ! Microphysical forcing of Tc [K/s]
   real(kind_phys) C_qt                                    ! Microphysical forcing of qt [kg/kg/s]
   real(kind_phys) dTcdt                                   ! d(Tc)/dt      [K/s]
   real(kind_phys) dqtdt                                   ! d(qt)/dt      [kg/kg/s]
   real(kind_phys) dqtstldt                                ! d(qt_alst)/dt [kg/kg/s]
   real(kind_phys) dqidt                                   ! d(qi)/dt      [kg/kg/s]

   real(kind_phys) dqlstdt                                 ! d(ql_st)/dt [kg/kg/s]
   real(kind_phys) dalstdt                                 ! d(al_st)/dt  [1/s]
   real(kind_phys) dastdt                                  ! d(a_st)/dt  [1/s]

   real(kind_phys) anic                                    ! Fractional area of non-cumulus and non-ice stratus fraction
   real(kind_phys) GG                                      ! G_nc(i,k)/anic

   real(kind_phys) aa(2,2)
   real(kind_phys) bb(2,1)

   real(kind_phys) qmin1(ncol,pver)
   real(kind_phys) qmin2(ncol,pver)
   real(kind_phys) qmin3(ncol,pver)

   real(kind_phys) esat_a(ncol)                            ! Saturation water vapor pressure [Pa]
   real(kind_phys) qsat_a(ncol,pver)                       ! Saturation water vapor specific humidity [kg/kg]
   real(kind_phys) Twb_aw(ncol)                            ! Wet-bulb temperature [K]
   real(kind_phys) qvwb_aw(ncol,pver)                      ! Wet-bulb water vapor specific humidity [kg/kg]

   real(kind_phys) esat_b(ncol)
   real(kind_phys) qsat_b(ncol)
   real(kind_phys) dqsdT_b(ncol)

   logical  land
   real(kind_phys) tmp

   real(kind_phys) rhmaxi_arr(ncol,pver)
   real(kind_phys) rhmini_arr(ncol,pver)
   real(kind_phys) rhminl_arr(ncol,pver)
   real(kind_phys) rhminl_adj_land_arr(ncol,pver)
   real(kind_phys) rhminh_arr(ncol,pver)

   real(kind_phys) QQmax,QQmin,QQwmin,QQimin                ! For limiting QQ
   real(kind_phys) cone                                     ! Number close to but smaller than 1

   ! Renamed output variables (al_st_star -> alst, ai_st_star -> aist, etc.)
   ! These are used directly from the output arguments

   ! Renamed input variables for mmacro_pcond interface
   ! A_T -> ttend, A_qv -> qtend, A_ql -> lmitend, A_qi -> itend, etc.
   ! The renaming is handled by using the computed tendencies directly

   ! ================================================================ !
   ! Begin computation                                                 !
   ! ================================================================ !

   errmsg = ''
   errflg = 0

   ! Set module-level physical constants from arguments
   cpair      = cpair_in
   latvap     = latvap_in
   latice     = latice_in
   top_lev_mod = top_lev

   ! --------------------------------------------------------- !
   ! Get constituent indices via CCPP constituent properties    !
   ! --------------------------------------------------------- !

   call ccpp_const_get_idx(const_props, 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', ixcldliq, errmsg, errflg)
   if (errflg /= 0) return
   call ccpp_const_get_idx(const_props, 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', ixcldice, errmsg, errflg)
   if (errflg /= 0) return
   call ccpp_const_get_idx(const_props, 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', ixq, errmsg, errflg)
   if (errflg /= 0) return

   ! ------------------------------------------ !
   ! Section 1: Compute advective tendencies    !
   ! (from macrop_driver lines 931-982)         !
   ! ------------------------------------------ !

   rdtime = 1._kind_phys / dtime

   zeros(:ncol,top_lev:pver)  = 0._kind_phys
   qc(:ncol,top_lev:pver) = cldliq(:ncol,top_lev:pver)
   qi_loc(:ncol,top_lev:pver) = cldice(:ncol,top_lev:pver)
   nc(:ncol,top_lev:pver) = numliq(:ncol,top_lev:pver)
   ni_loc(:ncol,top_lev:pver) = numice(:ncol,top_lev:pver)

 ! In CAM5, 'microphysical forcing' ( CC_... ) and 'the other advective forcings' ( ttend, ... )
 ! are separately provided into the prognostic microp_driver macrophysics scheme. This is an
 ! attempt to resolve in-cloud and out-cloud forcings.

   ! Make local copies of CC arrays since mmacro_pcond may need to zero them on nstep<=1
   CC_T_loc(:ncol,:)    = CC_T(:ncol,:)
   CC_qv_loc(:ncol,:)   = CC_qv(:ncol,:)
   CC_ql_loc(:ncol,:)   = CC_ql(:ncol,:)
   CC_qi_loc(:ncol,:)   = CC_qi(:ncol,:)
   CC_nl_loc(:ncol,:)   = CC_nl(:ncol,:)
   CC_ni_loc(:ncol,:)   = CC_ni(:ncol,:)
   CC_qlst_loc(:ncol,:) = CC_qlst(:ncol,:)

   if( nstep .le. 1 ) then
       tcwat(:ncol,top_lev:pver)   = t(:ncol,top_lev:pver)
       qcwat(:ncol,top_lev:pver)   = q_wv(:ncol,top_lev:pver)
       lcwat(:ncol,top_lev:pver)   = qc(:ncol,top_lev:pver) + qi_loc(:ncol,top_lev:pver)
       iccwat(:ncol,top_lev:pver)  = qi_loc(:ncol,top_lev:pver)
       nlwat(:ncol,top_lev:pver)   = nc(:ncol,top_lev:pver)
       niwat(:ncol,top_lev:pver)   = ni_loc(:ncol,top_lev:pver)
       ttend(:ncol,:)   = 0._kind_phys
       qtend(:ncol,:)   = 0._kind_phys
       ltend(:ncol,:)   = 0._kind_phys
       itend(:ncol,:)   = 0._kind_phys
       nltend(:ncol,:)  = 0._kind_phys
       nitend(:ncol,:)  = 0._kind_phys
       CC_T_loc(:ncol,:)    = 0._kind_phys
       CC_qv_loc(:ncol,:)   = 0._kind_phys
       CC_ql_loc(:ncol,:)   = 0._kind_phys
       CC_qi_loc(:ncol,:)   = 0._kind_phys
       CC_nl_loc(:ncol,:)   = 0._kind_phys
       CC_ni_loc(:ncol,:)   = 0._kind_phys
       CC_qlst_loc(:ncol,:) = 0._kind_phys
   else
       ttend(:ncol,top_lev:pver)   = ( t(:ncol,top_lev:pver)   -  tcwat(:ncol,top_lev:pver)) * rdtime &
            - CC_T_loc(:ncol,top_lev:pver)
       qtend(:ncol,top_lev:pver)   = ( q_wv(:ncol,top_lev:pver) -  qcwat(:ncol,top_lev:pver)) * rdtime &
            - CC_qv_loc(:ncol,top_lev:pver)
       ltend(:ncol,top_lev:pver)   = ( qc(:ncol,top_lev:pver) + qi_loc(:ncol,top_lev:pver) - lcwat(:ncol,top_lev:pver) ) * rdtime &
            - (CC_ql_loc(:ncol,top_lev:pver) + CC_qi_loc(:ncol,top_lev:pver))
       itend(:ncol,top_lev:pver)   = ( qi_loc(:ncol,top_lev:pver)         - iccwat(:ncol,top_lev:pver)) * rdtime &
            - CC_qi_loc(:ncol,top_lev:pver)
       nltend(:ncol,top_lev:pver)  = ( nc(:ncol,top_lev:pver)         -  nlwat(:ncol,top_lev:pver)) * rdtime &
            - CC_nl_loc(:ncol,top_lev:pver)
       nitend(:ncol,top_lev:pver)  = ( ni_loc(:ncol,top_lev:pver)         -  niwat(:ncol,top_lev:pver)) * rdtime &
            - CC_ni_loc(:ncol,top_lev:pver)
   endif
   lmitend(:ncol,top_lev:pver) = ltend(:ncol,top_lev:pver) - itend(:ncol,top_lev:pver)

   t_inout(:ncol,top_lev:pver)  =  tcwat(:ncol,top_lev:pver)
   qv_inout(:ncol,top_lev:pver) =  qcwat(:ncol,top_lev:pver)
   ql_inout(:ncol,top_lev:pver) =  lcwat(:ncol,top_lev:pver) - iccwat(:ncol,top_lev:pver)
   qi_inout(:ncol,top_lev:pver) = iccwat(:ncol,top_lev:pver)
   nl_inout(:ncol,top_lev:pver) =  nlwat(:ncol,top_lev:pver)
   ni_inout(:ncol,top_lev:pver) =  niwat(:ncol,top_lev:pver)

   ! --------------------------------------------------------- !
   ! Section 2: mmacro_pcond body                              !
   ! (from cldwat2m_macro.F90 mmacro_pcond, lines 444-1220)   !
   ! The call maps are:                                        !
   !   T0 -> t_inout, qv0 -> qv_inout, etc.                   !
   !   A_T -> ttend, A_qv -> qtend, A_ql -> lmitend           !
   !   A_qi -> itend, A_nl -> nltend, A_ni -> nitend           !
   !   C_T -> CC_T_loc, etc.                                   !
   !   D_T -> dlf_T, etc.                                      !
   !   a_cud -> concld_old, a_cu0 -> concld                    !
   !   p -> pmid, dp -> pdel, dt -> dtime                      !
   !   al_st_star -> alst, ai_st_star -> aist                  !
   !   ql_st_star -> qlst, qi_st_star -> qist                  !
   !   s_tendout -> tlat, qv_tendout -> qvlat                  !
   !   ql_tendout -> qcten, qi_tendout -> qiten                !
   !   nl_tendout -> ncten, ni_tendout -> niten                !
   !   qme -> cmeliq                                           !
   ! --------------------------------------------------------- !

   cone            = 0.999_kind_phys
   zeros(:ncol,:)  = 0._kind_phys

   ! ------------------------------------ !
   ! Global initialization of main output !
   ! ------------------------------------ !

     tlat(:ncol,:)      = 0._kind_phys
     qvlat(:ncol,:)     = 0._kind_phys
     qcten(:ncol,:)     = 0._kind_phys
     qiten(:ncol,:)     = 0._kind_phys
     ncten(:ncol,:)     = 0._kind_phys
     niten(:ncol,:)     = 0._kind_phys

     cmeliq(:ncol,:)    = 0._kind_phys

     cld(:ncol,:)       = 0._kind_phys
     alst(:ncol,:)      = 0._kind_phys
     aist(:ncol,:)      = 0._kind_phys
     qlst(:ncol,:)      = 0._kind_phys
     qist(:ncol,:)      = 0._kind_phys

   ! --------------------------------------- !
   ! Initialization of internal 2D variables !
   ! --------------------------------------- !

     T_loc(:ncol,:)     = 0._kind_phys
     T1(:ncol,:)        = 0._kind_phys
     T_0(:ncol,:)       = 0._kind_phys
     T_05(:ncol,:)      = 0._kind_phys
     T_prime0(:ncol,:)  = 0._kind_phys
     T_dprime(:ncol,:)  = 0._kind_phys
     T_star(:ncol,:)    = 0._kind_phys

     qv(:ncol,:)        = 0._kind_phys
     qv1(:ncol,:)       = 0._kind_phys
     qv_0(:ncol,:)      = 0._kind_phys
     qv_05(:ncol,:)     = 0._kind_phys
     qv_prime0(:ncol,:) = 0._kind_phys
     qv_dprime(:ncol,:) = 0._kind_phys
     qv_star(:ncol,:)   = 0._kind_phys

     ql(:ncol,:)        = 0._kind_phys
     ql1(:ncol,:)       = 0._kind_phys
     ql_0(:ncol,:)      = 0._kind_phys
     ql_05(:ncol,:)     = 0._kind_phys
     ql_prime0(:ncol,:) = 0._kind_phys
     ql_dprime(:ncol,:) = 0._kind_phys
     ql_star(:ncol,:)   = 0._kind_phys

     qi_mm(:ncol,:)     = 0._kind_phys
     qi1(:ncol,:)       = 0._kind_phys
     qi_0(:ncol,:)      = 0._kind_phys
     qi_05(:ncol,:)     = 0._kind_phys
     qi_prime0(:ncol,:) = 0._kind_phys
     qi_dprime(:ncol,:) = 0._kind_phys
     qi_star(:ncol,:)   = 0._kind_phys

     nl(:ncol,:)        = 0._kind_phys
     nl1(:ncol,:)       = 0._kind_phys
     nl_0(:ncol,:)      = 0._kind_phys
     nl_05(:ncol,:)     = 0._kind_phys
     nl_prime0(:ncol,:) = 0._kind_phys
     nl_dprime(:ncol,:) = 0._kind_phys
     nl_star(:ncol,:)   = 0._kind_phys

     ni(:ncol,:)        = 0._kind_phys
     ni1(:ncol,:)       = 0._kind_phys
     ni_0(:ncol,:)      = 0._kind_phys
     ni_05(:ncol,:)     = 0._kind_phys
     ni_prime0(:ncol,:) = 0._kind_phys
     ni_dprime(:ncol,:) = 0._kind_phys
     ni_star(:ncol,:)   = 0._kind_phys

     a_st(:ncol,:)      = 0._kind_phys
     a_st_0(:ncol,:)    = 0._kind_phys
     a_st_star(:ncol,:) = 0._kind_phys

     al_st(:ncol,:)     = 0._kind_phys
     al_st_0(:ncol,:)   = 0._kind_phys
     al_st_nc(:ncol,:)  = 0._kind_phys

     ai_st(:ncol,:)     = 0._kind_phys
     ai_st_0(:ncol,:)   = 0._kind_phys
     ai_st_nc(:ncol,:)  = 0._kind_phys

     ql_st(:ncol,:)     = 0._kind_phys
     ql_st_0(:ncol,:)   = 0._kind_phys

     qi_st(:ncol,:)     = 0._kind_phys
     qi_st_0(:ncol,:)   = 0._kind_phys

 ! Cumulus properties

     dacudt(:ncol,:)    = 0._kind_phys
     a_cu(:ncol,:)      = 0._kind_phys

 ! Adjustment tendency in association with 'positive_moisture'

     Tten_pwi1(:ncol,:)  = 0._kind_phys
     qvten_pwi1(:ncol,:) = 0._kind_phys
     qlten_pwi1(:ncol,:) = 0._kind_phys
     qiten_pwi1(:ncol,:) = 0._kind_phys
     nlten_pwi1(:ncol,:) = 0._kind_phys
     niten_pwi1(:ncol,:) = 0._kind_phys

     Tten_pwi2(:ncol,:)  = 0._kind_phys
     qvten_pwi2(:ncol,:) = 0._kind_phys
     qlten_pwi2(:ncol,:) = 0._kind_phys
     qiten_pwi2(:ncol,:) = 0._kind_phys
     nlten_pwi2(:ncol,:) = 0._kind_phys
     niten_pwi2(:ncol,:) = 0._kind_phys

     A_T_adj(:ncol,:)   = 0._kind_phys
     A_qv_adj(:ncol,:)  = 0._kind_phys
     A_ql_adj(:ncol,:)  = 0._kind_phys
     A_qi_adj(:ncol,:)  = 0._kind_phys
     A_nl_adj(:ncol,:)  = 0._kind_phys
     A_ni_adj(:ncol,:)  = 0._kind_phys

     qvadj   (:ncol,:)  = 0._kind_phys
     qladj   (:ncol,:)  = 0._kind_phys
     qiadj   (:ncol,:)  = 0._kind_phys

 ! Adjustment tendency in association with 'instratus_condensate'

     QQw1(:ncol,:)      = 0._kind_phys
     QQi1(:ncol,:)      = 0._kind_phys
     QQw2(:ncol,:)      = 0._kind_phys
     QQi2(:ncol,:)      = 0._kind_phys

     QQnl1(:ncol,:)     = 0._kind_phys
     QQni1(:ncol,:)     = 0._kind_phys
     QQnl2(:ncol,:)     = 0._kind_phys
     QQni2(:ncol,:)     = 0._kind_phys

     QQnl(:ncol,:)      = 0._kind_phys
     QQni(:ncol,:)      = 0._kind_phys

 ! Macrophysical process tendency variables

     QQ(:ncol,:)        = 0._kind_phys
     QQw(:ncol,:)       = 0._kind_phys
     QQi(:ncol,:)       = 0._kind_phys
     QQnl(:ncol,:)      = 0._kind_phys
     QQni(:ncol,:)      = 0._kind_phys
     ACnl(:ncol,:)      = 0._kind_phys
     ACni(:ncol,:)      = 0._kind_phys

     QQw_prev(:ncol,:)  = 0._kind_phys
     QQi_prev(:ncol,:)  = 0._kind_phys
     QQnl_prev(:ncol,:) = 0._kind_phys
     QQni_prev(:ncol,:) = 0._kind_phys

     QQw_prog(:ncol,:)  = 0._kind_phys
     QQi_prog(:ncol,:)  = 0._kind_phys
     QQnl_prog(:ncol,:) = 0._kind_phys
     QQni_prog(:ncol,:) = 0._kind_phys

     QQ_final(:ncol,:)  = 0._kind_phys
     QQw_final(:ncol,:) = 0._kind_phys
     QQi_final(:ncol,:) = 0._kind_phys
     QQn_final(:ncol,:) = 0._kind_phys
     QQnl_final(:ncol,:) = 0._kind_phys
     QQni_final(:ncol,:) = 0._kind_phys

     QQ_all(:ncol,:)    = 0._kind_phys
     QQw_all(:ncol,:)   = 0._kind_phys
     QQi_all(:ncol,:)   = 0._kind_phys
     QQn_all(:ncol,:)   = 0._kind_phys
     QQnl_all(:ncol,:)  = 0._kind_phys
     QQni_all(:ncol,:)  = 0._kind_phys

 ! Coefficient for computing QQ and related processes

     U(:ncol,:)         = 0._kind_phys
     U_nc(:ncol,:)      = 0._kind_phys
     G_nc(:ncol,:)      = 0._kind_phys
     F_nc(:ncol,:)      = 0._kind_phys

 ! Other

     qmin1(:ncol,:)     = 0._kind_phys
     qmin2(:ncol,:)     = 0._kind_phys
     qmin3(:ncol,:)     = 0._kind_phys

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   rhmaxi_arr(:ncol,:)          = rhmaxi_const
   rhmini_arr(:ncol,:)          = rhmini_const
   rhminl_arr(:ncol,:)          = rhminl_const
   rhminl_adj_land_arr(:ncol,:) = rhminl_adj_land_const
   rhminh_arr(:ncol,:)          = rhminh_const

   ! ---------------------------------- !
   ! Compute cumulus-related properties !
   ! ---------------------------------- !

   dacudt(:ncol,top_lev:pver) = &
        (concld(:ncol,top_lev:pver) - concld_old(:ncol,top_lev:pver))/dtime

   ! ---------------------------------------------------------------------- !
   ! set to zero for levels above
   ! ---------------------------------------------------------------------- !
   ql_inout(:ncol,:top_lev-1) = 0._kind_phys
   qi_inout(:ncol,:top_lev-1) = 0._kind_phys
   nl_inout(:ncol,:top_lev-1) = 0._kind_phys
   ni_inout(:ncol,:top_lev-1) = 0._kind_phys

   ! ---------------------------------------------------------------------- !
   ! Check if input non-cumulus pixels satisfie a non-negative constraint.  !
   ! If not, force all water vapor substances to be positive in all layers. !
   ! We should use 'old' cumulus properties for this routine.               !
   ! ---------------------------------------------------------------------- !

   T1(:ncol,:)    =  t_inout(:ncol,:)
   qv1(:ncol,:)   = qv_inout(:ncol,:)
   ql1(:ncol,:)   = ql_inout(:ncol,:)
   qi1(:ncol,:)   = qi_inout(:ncol,:)
   nl1(:ncol,:)   = nl_inout(:ncol,:)
   ni1(:ncol,:)   = ni_inout(:ncol,:)

   qmin1(:ncol,:) = qmin(ixq)
   qmin2(:ncol,:) = qmin(ixcldliq)
   qmin3(:ncol,:) = qmin(ixcldice)

   call positive_moisture( ncol, dtime, qmin1, qmin2, qmin3, pdel, &
                           qv1, ql1, qi1, T1, qvten_pwi1, qlten_pwi1, &
                           qiten_pwi1, Tten_pwi1, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql1(i,k) .lt. qsmall ) then
          nlten_pwi1(i,k) = -nl1(i,k)/dtime
          nl1(i,k)        = 0._kind_phys
      endif
      if( qi1(i,k) .lt. qsmall ) then
          niten_pwi1(i,k) = -ni1(i,k)/dtime
          ni1(i,k)        = 0._kind_phys
      endif
   enddo
   enddo

   ! ------------------------------------------------------------- !
   ! Impose 'in-stratus condensate amount constraint'              !
   ! such that it is bounded by two limiting values.               !
   ! This should also use 'old' cumulus properties since it is     !
   ! before applying external forcings.                            !
   ! Below 'QQw1,QQi1' are effective adjustive condensation        !
   ! Although this process also involves freezing of cloud         !
   ! liquid into ice, they can be and only can be expressed        !
   ! in terms of effective condensation.                           !
   ! ------------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( ncol, k,                                    &
                                 pmid(:,k), T1(:,k), qv1(:,k), ql1(:,k), qi1(:,k),    &
                                 ni1(:,k),                                   &
                                 concld_old(:,k), zeros(:,k), zeros(:,k),    &
                                 zeros(:,k), zeros(:,k), zeros(:,k),        &
                                 landfrac, snowh,                            &
                                 rhmaxi_arr(:,k),rhmini_arr(:,k), rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k), &
                                 T_0(:,k), qv_0(:,k), ql_0(:,k), qi_0(:,k), &
                                 al_st_0(:,k), ai_st_0(:,k), ql_st_0(:,k), qi_st_0(:,k) )
      a_st_0(:ncol,k) = max(al_st_0(:ncol,k),ai_st_0(:ncol,k))
      QQw1(:ncol,k)   = (ql_0(:ncol,k) - ql1(:ncol,k))/dtime
      QQi1(:ncol,k)   = (qi_0(:ncol,k) - qi1(:ncol,k))/dtime
      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! Set a limit such that negative state not happens.  !
      ! -------------------------------------------------- !
      do i = 1, ncol
         if( QQw1(i,k) .le. 0._kind_phys ) then
             if( ql1(i,k) .gt. qsmall ) then
                 QQnl1(i,k) = QQw1(i,k)*nl1(i,k)/ql1(i,k)
                 QQnl1(i,k) = min(0._kind_phys,cone*max(QQnl1(i,k),-nl1(i,k)/dtime))
             else
                 QQnl1(i,k) = 0._kind_phys
             endif
         endif
         if( QQi1(i,k) .le. 0._kind_phys ) then
             if( qi1(i,k) .gt. qsmall ) then
                 QQni1(i,k) = QQi1(i,k)*ni1(i,k)/qi1(i,k)
                 QQni1(i,k) = min(0._kind_phys,cone*max(QQni1(i,k),-ni1(i,k)/dtime))
             else
                 QQni1(i,k) = 0._kind_phys
             endif
         endif
      enddo
   enddo
   nl_0(:ncol,top_lev:) = max(0._kind_phys,nl1(:ncol,top_lev:)+QQnl1(:ncol,top_lev:)*dtime)
   ni_0(:ncol,top_lev:) = max(0._kind_phys,ni1(:ncol,top_lev:)+QQni1(:ncol,top_lev:)*dtime)

   ! ----------------------------------------------------------------------------- !
   ! Check if non-cumulus pixels of '_05' state satisfies non-negative constraint. !
   ! If not, force all water substances of '_05' state to be positive by imposing  !
   ! adjustive advection. We should use 'new' cumulus properties for this routine. !
   ! ----------------------------------------------------------------------------- !

   T_05(:ncol,top_lev:)  =  T_0(:ncol,top_lev:) + (  ttend(:ncol,top_lev:) +  CC_T_loc(:ncol,top_lev:) ) * dtime
   qv_05(:ncol,top_lev:) = qv_0(:ncol,top_lev:) + ( qtend(:ncol,top_lev:) + CC_qv_loc(:ncol,top_lev:) ) * dtime
   ql_05(:ncol,top_lev:) = ql_0(:ncol,top_lev:) + ( lmitend(:ncol,top_lev:) + CC_ql_loc(:ncol,top_lev:) ) * dtime
   qi_05(:ncol,top_lev:) = qi_0(:ncol,top_lev:) + ( itend(:ncol,top_lev:) + CC_qi_loc(:ncol,top_lev:) ) * dtime
   nl_05(:ncol,top_lev:) = max(0._kind_phys, nl_0(:ncol,top_lev:) + ( nltend(:ncol,top_lev:) + CC_nl_loc(:ncol,top_lev:) ) * dtime )
   ni_05(:ncol,top_lev:) = max(0._kind_phys, ni_0(:ncol,top_lev:) + ( nitend(:ncol,top_lev:) + CC_ni_loc(:ncol,top_lev:) ) * dtime )

   call positive_moisture( ncol, dtime, qmin1, qmin2, qmin3, pdel, &
                           qv_05, ql_05, qi_05, T_05, A_qv_adj, &
                           A_ql_adj, A_qi_adj, A_T_adj, do_cldice)

   ! -------------------------------------------------------------- !
   ! Define reference state at the first iteration. This will be    !
   ! continuously updated within the iteration loop below.          !
   ! While equlibrium state properties are already output from the  !
   ! 'instratus_condensate', they will be re-computed within the    !
   ! each iteration process. At the first iteration, they will      !
   ! produce exactly identical results. Note that except at the     !
   ! very first iteration iter = 1, we must use updated cumulus     !
   ! properties at all the other iteration processes. Even at the   !
   ! first iteration, we should use updated cumulus properties      !
   ! when computing limiters for (Q,P,E).                           !
   ! -------------------------------------------------------------- !

   ! -------------------------------------------------------------- !
   ! Define variables at the reference state of the first iteration !
   ! -------------------------------------------------------------- !

   T_loc(:ncol,top_lev:)  = T_0(:ncol,top_lev:)
   qv(:ncol,top_lev:)     = qv_0(:ncol,top_lev:)
   ql(:ncol,top_lev:)     = ql_0(:ncol,top_lev:)
   qi_mm(:ncol,top_lev:)  = qi_0(:ncol,top_lev:)
   al_st(:ncol,top_lev:)  = al_st_0(:ncol,top_lev:)
   ai_st(:ncol,top_lev:)  = ai_st_0(:ncol,top_lev:)
   a_st(:ncol,top_lev:)   = a_st_0(:ncol,top_lev:)
   ql_st(:ncol,top_lev:)  = ql_st_0(:ncol,top_lev:)
   qi_st(:ncol,top_lev:)  = qi_st_0(:ncol,top_lev:)
   nl(:ncol,top_lev:)     = nl_0(:ncol,top_lev:)
   ni(:ncol,top_lev:)     = ni_0(:ncol,top_lev:)

   ! -------------------------- !
   ! Main iterative computation !
   ! -------------------------- !

   do k = top_lev, pver
      call findsp_vc(qv_05(:ncol,k), T_05(:ncol,k), pmid(:ncol,k), .false., &
                     Twb_aw(:ncol), qvwb_aw(:ncol,k))
      call qsat_water(T_05(1:ncol,k), pmid(1:ncol,k), esat_a(1:ncol), qsat_a(1:ncol,k), ncol)
   enddo

   do iter = 1, niter

      ! ------------------------------------------ !
      ! Initialize array within the iteration loop !
      ! ------------------------------------------ !

      QQ(:,:)         = 0._kind_phys
      QQw(:,:)        = 0._kind_phys
      QQi(:,:)        = 0._kind_phys
      QQnl(:,:)       = 0._kind_phys
      QQni(:,:)       = 0._kind_phys
      QQw2(:,:)       = 0._kind_phys
      QQi2(:,:)       = 0._kind_phys
      QQnl2(:,:)      = 0._kind_phys
      QQni2(:,:)      = 0._kind_phys
      nlten_pwi2(:,:) = 0._kind_phys
      niten_pwi2(:,:) = 0._kind_phys
      ACnl(:,:)       = 0._kind_phys
      ACni(:,:)       = 0._kind_phys
      aa(:,:)         = 0._kind_phys
      bb(:,:)         = 0._kind_phys

      do k = top_lev, pver
         call qsat_water(T_loc(1:ncol,k), pmid(1:ncol,k), esat_b(1:ncol), qsat_b(1:ncol), ncol, dqsdt=dqsdT_b(1:ncol))
      if( iter .eq. 1 ) then
          a_cu(:ncol,k) = concld_old(:ncol,k)
      else
          a_cu(:ncol,k) = concld(:ncol,k)
      endif
      do i = 1, ncol
         U(i,k)    =  qv(i,k)/qsat_b(i)
         U_nc(i,k) =  U(i,k)
      enddo
      if( CAMstfrac ) then
          call astG_RHU(U_nc(:,k),pmid(:,k),qv(:,k),landfrac(:),snowh(:),al_st_nc(:,k),G_nc(:,k),ncol,&
                        rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k))
      else
          call astG_PDF(U_nc(:,k),pmid(:,k),qv(:,k),landfrac(:),snowh(:),al_st_nc(:,k),G_nc(:,k),ncol,&
                        rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k))
      endif
      call aist_vector(qv(:,k),T_loc(:,k),pmid(:,k),qi_mm(:,k),ni(:,k),landfrac(:),snowh(:),ai_st_nc(:,k),ncol,&
                       rhmaxi_arr(:,k), rhmini_arr(:,k), rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k))

      ai_st(:ncol,k)  =  (1._kind_phys-a_cu(:ncol,k))*ai_st_nc(:ncol,k)
      al_st(:ncol,k)  =  (1._kind_phys-a_cu(:ncol,k))*al_st_nc(:ncol,k)
      a_st(:ncol,k)   =  max(al_st(:ncol,k),ai_st(:ncol,k))

      do i = 1, ncol

         ! -------------------------------------------------------- !
         ! Compute basic thermodynamic coefficients for computing Q !
         ! -------------------------------------------------------- !

         alpha  =  1._kind_phys/qsat_b(i)
         beta   =  dqsdT_b(i)*(qv(i,k)/qsat_b(i)**2)
         betast =  alpha*dqsdT_b(i)
         gammal =  alpha + (latvap/cpair)*beta
         gammai =  alpha + ((latvap+latice)/cpair)*beta
         gammaQ =  alpha + (latvap/cpair)*beta
         deltal =  1._kind_phys + a_st(i,k)*(latvap/cpair)*(betast/alpha)
         deltai =  1._kind_phys + a_st(i,k)*((latvap+latice)/cpair)*(betast/alpha)
         A_Tc   =  ttend(i,k)+A_T_adj(i,k)-(latvap/cpair)*(lmitend(i,k)+A_ql_adj(i,k))-((latvap+latice)/cpair)*(itend(i,k)+A_qi_adj(i,k))
         A_qt   =  qtend(i,k) + A_qv_adj(i,k) + lmitend(i,k) + A_ql_adj(i,k) + itend(i,k) + A_qi_adj(i,k)
         C_Tc   =  CC_T_loc(i,k) - (latvap/cpair)*CC_ql_loc(i,k) - ((latvap+latice)/cpair)*CC_qi_loc(i,k)
         C_qt   =  CC_qv_loc(i,k) + CC_ql_loc(i,k) + CC_qi_loc(i,k)
         dTcdt  =  A_Tc + C_Tc
         dqtdt  =  A_qt + C_qt
       ! dqtstldt = A_qt + C_ql_loc(i,k)/max(1.e-2_kind_phys,al_st(i,k))                             ! Original
       ! dqtstldt = A_qt - itend(i,k) - A_qi_adj(i,k) + C_ql_loc(i,k)/max(1.e-2_kind_phys,al_st(i,k)) ! New 1 on Dec.30.2009.
         dqtstldt = A_qt - itend(i,k) - A_qi_adj(i,k) + CC_qlst_loc(i,k)                        ! New 2 on Dec.30.2009.
       ! dqtstldt = A_qt + C_qt                                                           ! Original Conservative treatment
       ! dqtstldt = A_qt - itend(i,k) - A_qi_adj(i,k) + C_qt - CC_qi_loc(i,k)            ! New Conservative treatment on Dec.30.2009
         dqidt = itend(i,k) + A_qi_adj(i,k) + CC_qi_loc(i,k)

         anic    = max(1.e-8_kind_phys,(1._kind_phys-a_cu(i,k)))
         GG      = G_nc(i,k)/anic
         aa(1,1) = gammal*al_st(i,k)
         aa(1,2) = GG + gammal*cc*ql_st(i,k)
         aa(2,1) = alpha + (latvap/cpair)*betast*al_st(i,k)
         aa(2,2) = (latvap/cpair)*betast*cc*ql_st(i,k)
         bb(1,1) = alpha*dqtdt - beta*dTcdt - gammai*dqidt - GG*al_st_nc(i,k)*dacudt(i,k) + F_nc(i,k)
         bb(2,1) = alpha*dqtstldt - betast*(dTcdt + ((latvap+latice)/cpair)*dqidt)
         call gaussj(aa(1:2,1:2),2,2,bb(1:2,1),1,1)
         dqlstdt = bb(1,1)
         dalstdt = bb(2,1)
         QQ(i,k) = al_st(i,k)*dqlstdt + cc*ql_st(i,k)*dalstdt - ( lmitend(i,k) + A_ql_adj(i,k) + CC_ql_loc(i,k) )

       ! ------------------------------------------------------------ !
       ! Limiter for QQ                                               !
       ! Here, 'fice' should be from the reference equilibrium state  !
       ! since QQ itself is computed from the reference state.        !
       ! From the assumption used for derivation of QQ(i), it must be !
       ! that QQw(i) = QQ(i)*(1._kind_phys-fice(i)), QQi(i) = QQ(i)*fice(i)  !
       ! ------------------------------------------------------------ !

         if( QQ(i,k) .ge. 0._kind_phys ) then
             QQmax    = (qv_05(i,k) - qmin(ixq))/dtime ! For ghost cumulus & semi-ghost ice stratus
             QQmax    = max(0._kind_phys,QQmax)
             QQ(i,k)  = min(QQ(i,k),QQmax)
             QQw(i,k) = QQ(i,k)
             QQi(i,k) = 0._kind_phys
         else
             QQmin  = 0._kind_phys
             if( qv_05(i,k) .lt. qsat_a(i,k) ) QQmin = min(0._kind_phys,cone*(qv_05(i,k)-qvwb_aw(i,k))/dtime)
             QQ(i,k)  = max(QQ(i,k),QQmin)
             QQw(i,k) = QQ(i,k)
             QQi(i,k) = 0._kind_phys
             QQwmin   = min(0._kind_phys,-cone*ql_05(i,k)/dtime)
             QQimin   = min(0._kind_phys,-cone*qi_05(i,k)/dtime)
             QQw(i,k) = min(0._kind_phys,max(QQw(i,k),QQwmin))
             QQi(i,k) = min(0._kind_phys,max(QQi(i,k),QQimin))
         endif

       ! -------------------------------------------------- !
       ! Reduce droplet concentration if evaporation occurs !
       ! Note 'QQnl1,QQni1' are computed from the reference !
       ! equilibrium state but limiter is from 'nl_05'.     !
       ! -------------------------------------------------- !

         if( QQw(i,k) .lt. 0._kind_phys ) then
             if( ql_05(i,k) .gt. qsmall ) then
                 QQnl(i,k) = QQw(i,k)*nl_05(i,k)/ql_05(i,k)
                 QQnl(i,k) = min(0._kind_phys,cone*max(QQnl(i,k),-nl_05(i,k)/dtime))
             else
                 QQnl(i,k) = 0._kind_phys
             endif
         endif

         if( QQi(i,k) .lt. 0._kind_phys ) then
             if( qi_05(i,k) .gt. qsmall ) then
                 QQni(i,k) = QQi(i,k)*ni_05(i,k)/qi_05(i,k)
                 QQni(i,k) = min(0._kind_phys,cone*max(QQni(i,k),-ni_05(i,k)/dtime))
             else
                 QQni(i,k) = 0._kind_phys
             endif
         endif

      enddo
      enddo

    ! -------------------------------------------------------------------- !
    ! Until now, we have finished computing all necessary tendencies       !
    ! from the equilibrium input state (T_0).                              !
    ! If ramda = 0 : fully explicit scheme                                 !
    !    ramda = 1 : fully implicit scheme                                 !
    ! Note that 'ramda = 0.5 with niter = 2' can mimic                     !
    ! -------------------------------------------------------------------- !

      if( iter .eq. 1 ) then
          QQw_prev(:ncol,top_lev:)  = QQw(:ncol,top_lev:)
          QQi_prev(:ncol,top_lev:)  = QQi(:ncol,top_lev:)
          QQnl_prev(:ncol,top_lev:) = QQnl(:ncol,top_lev:)
          QQni_prev(:ncol,top_lev:) = QQni(:ncol,top_lev:)
      endif

      QQw_prog(:ncol,top_lev:)   = ramda*QQw(:ncol,top_lev:)   + (1._kind_phys-ramda)*QQw_prev(:ncol,top_lev:)
      QQi_prog(:ncol,top_lev:)   = ramda*QQi(:ncol,top_lev:)   + (1._kind_phys-ramda)*QQi_prev(:ncol,top_lev:)
      QQnl_prog(:ncol,top_lev:)  = ramda*QQnl(:ncol,top_lev:)  + (1._kind_phys-ramda)*QQnl_prev(:ncol,top_lev:)
      QQni_prog(:ncol,top_lev:)  = ramda*QQni(:ncol,top_lev:)  + (1._kind_phys-ramda)*QQni_prev(:ncol,top_lev:)

      QQw_prev(:ncol,top_lev:)   = QQw_prog(:ncol,top_lev:)
      QQi_prev(:ncol,top_lev:)   = QQi_prog(:ncol,top_lev:)
      QQnl_prev(:ncol,top_lev:)  = QQnl_prog(:ncol,top_lev:)
      QQni_prev(:ncol,top_lev:)  = QQni_prog(:ncol,top_lev:)

    ! -------------------------------------------------------- !
    ! Compute final prognostic state on which final diagnostic !
    ! in-stratus condensate adjustment is applied in the below.!
    ! Important : I must check whether there are any external  !
    !             advective forcings of 'A_nl(i,k),A_ni(i,k)'. !
    !             Even they are (i.e., advection of aerosol),  !
    !             actual droplet activation will be performd   !
    !             in microphysics, so it will be completely    !
    !             reasonable to 'A_nl(i,k)=A_ni(i,k)=0'.       !
    ! -------------------------------------------------------- !

    do k = top_lev, pver
    do i = 1, ncol
       T_prime0(i,k)  = T_0(i,k)  + dtime*( ttend(i,k)  +  A_T_adj(i,k) +  CC_T_loc(i,k) + &
            (latvap*QQw_prog(i,k)+(latvap+latice)*QQi_prog(i,k))/cpair )
       qv_prime0(i,k) = qv_0(i,k) + dtime*( qtend(i,k) + A_qv_adj(i,k) + CC_qv_loc(i,k) - QQw_prog(i,k) - QQi_prog(i,k) )
       ql_prime0(i,k) = ql_0(i,k) + dtime*( lmitend(i,k) + A_ql_adj(i,k) + CC_ql_loc(i,k) + QQw_prog(i,k) )
       qi_prime0(i,k) = qi_0(i,k) + dtime*( itend(i,k) + A_qi_adj(i,k) + CC_qi_loc(i,k) + QQi_prog(i,k) )
       nl_prime0(i,k) = max(0._kind_phys,nl_0(i,k) + dtime*( nltend(i,k) + CC_nl_loc(i,k) + QQnl_prog(i,k) ))
       ni_prime0(i,k) = max(0._kind_phys,ni_0(i,k) + dtime*( nitend(i,k) + CC_ni_loc(i,k) + QQni_prog(i,k) ))
       if( ql_prime0(i,k) .lt. qsmall ) nl_prime0(i,k) = 0._kind_phys
       if( qi_prime0(i,k) .lt. qsmall ) ni_prime0(i,k) = 0._kind_phys
    enddo
    enddo

   ! -------------------------------------------------- !
   ! Perform diagnostic 'positive_moisture' constraint. !
   ! -------------------------------------------------- !

   T_dprime(:ncol,top_lev:)  =  T_prime0(:ncol,top_lev:)
   qv_dprime(:ncol,top_lev:) = qv_prime0(:ncol,top_lev:)
   ql_dprime(:ncol,top_lev:) = ql_prime0(:ncol,top_lev:)
   qi_dprime(:ncol,top_lev:) = qi_prime0(:ncol,top_lev:)
   nl_dprime(:ncol,top_lev:) = nl_prime0(:ncol,top_lev:)
   ni_dprime(:ncol,top_lev:) = ni_prime0(:ncol,top_lev:)

   call positive_moisture( ncol, dtime, qmin1, qmin2, qmin3, pdel,          &
                           qv_dprime, ql_dprime, qi_dprime, T_dprime,  &
                           qvten_pwi2, qlten_pwi2, qiten_pwi2, Tten_pwi2, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql_dprime(i,k) .lt. qsmall ) then
          nlten_pwi2(i,k) = -nl_dprime(i,k)/dtime
          nl_dprime(i,k)   = 0._kind_phys
      endif
      if( qi_dprime(i,k) .lt. qsmall ) then
          niten_pwi2(i,k) = -ni_dprime(i,k)/dtime
          ni_dprime(i,k)   = 0._kind_phys
      endif
   enddo
   enddo

   ! -------------------------------------------------------------- !
   ! Add tendency associated with detrainment of cumulus condensate !
   ! This tendency is not used in computing Q                       !
   ! Since D_ql,D_qi,D_nl,D_ni > 0, don't need to worry about       !
   ! negative scalar.                                               !
   ! This tendency is not reflected into Fzs2, which is OK.         !
   ! -------------------------------------------------------------- !

   T_dprime(:ncol,top_lev:)   =  T_dprime(:ncol,top_lev:)  + dlf_T(:ncol,top_lev:) * dtime
   qv_dprime(:ncol,top_lev:)  = qv_dprime(:ncol,top_lev:) + dlf_qv(:ncol,top_lev:) * dtime
   ql_dprime(:ncol,top_lev:)  = ql_dprime(:ncol,top_lev:) + dlf_ql(:ncol,top_lev:) * dtime
   qi_dprime(:ncol,top_lev:)  = qi_dprime(:ncol,top_lev:) + dlf_qi(:ncol,top_lev:) * dtime
   nl_dprime(:ncol,top_lev:)  = nl_dprime(:ncol,top_lev:) + dlf_nl(:ncol,top_lev:) * dtime
   ni_dprime(:ncol,top_lev:)  = ni_dprime(:ncol,top_lev:) + dlf_ni(:ncol,top_lev:) * dtime

   ! ---------------------------------------------------------- !
   ! Impose diagnostic upper and lower limits on the in-stratus !
   ! condensate amount. This produces a final equilibrium state !
   ! at the end of each iterative process.                      !
   ! ---------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( ncol           , k              , pmid(:,k)      , &
                                 T_dprime(:,k)  , qv_dprime(:,k) , ql_dprime(:,k) , qi_dprime(:,k), &
                                 ni_dprime(:,k) ,                                                   &
                                 concld(:,k)    , zeros(:,k)     , zeros(:,k)     ,                 &
                                 zeros(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 &
                                 landfrac       , snowh          ,                                  &
                                 rhmaxi_arr(:,k),rhmini_arr(:,k), rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k), &
                                 T_star(:,k)    , qv_star(:,k)   , ql_star(:,k)   , qi_star(:,k)  , &
                                 alst(:,k)      , aist(:,k)      , qlst(:,k)      , qist(:,k) )
      a_st_star(:ncol,k)  = max(alst(:ncol,k),aist(:ncol,k))
      QQw2(:ncol,k) = (ql_star(:ncol,k) - ql_dprime(:ncol,k))/dtime
      QQi2(:ncol,k) = (qi_star(:ncol,k) - qi_dprime(:ncol,k))/dtime
      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! -------------------------------------------------- !
      do i = 1, ncol
         if( QQw2(i,k) .le. 0._kind_phys ) then
             if( ql_dprime(i,k) .ge. qsmall ) then
                 QQnl2(i,k) = QQw2(i,k)*nl_dprime(i,k)/ql_dprime(i,k)
                 QQnl2(i,k) = min(0._kind_phys,cone*max(QQnl2(i,k),-nl_dprime(i,k)/dtime))
             else
                 QQnl2(i,k) = 0._kind_phys
             endif
         endif
         if( QQi2(i,k) .le. 0._kind_phys ) then
             if( qi_dprime(i,k) .gt. qsmall ) then
                 QQni2(i,k) = QQi2(i,k)*ni_dprime(i,k)/qi_dprime(i,k)
                 QQni2(i,k) = min(0._kind_phys,cone*max(QQni2(i,k),-ni_dprime(i,k)/dtime))
             else
                 QQni2(i,k) = 0._kind_phys
             endif
         endif
      enddo
   enddo
   nl_star(:ncol,top_lev:) = max(0._kind_phys,nl_dprime(:ncol,top_lev:)+QQnl2(:ncol,top_lev:)*dtime)
   ni_star(:ncol,top_lev:) = max(0._kind_phys,ni_dprime(:ncol,top_lev:)+QQni2(:ncol,top_lev:)*dtime)

   ! ------------------------------------------ !
   ! Final adjustment of droplet concentration. !
   ! Set # to zero if there is no cloud.        !
   ! ------------------------------------------ !

   do k = top_lev, pver
   do i = 1, ncol
      if( ql_star(i,k) .lt. qsmall ) then
          ACnl(i,k) = - nl_star(i,k)/dtime
          nl_star(i,k) = 0._kind_phys
      endif
      if( qi_star(i,k) .lt. qsmall ) then
          ACni(i,k) = - ni_star(i,k)/dtime
          ni_star(i,k) = 0._kind_phys
      endif
   enddo
   enddo

   ! ----------------------------------------------------- !
   ! Define equilibrium reference state for next iteration !
   ! ----------------------------------------------------- !

   T_loc(:ncol,top_lev:)  = T_star(:ncol,top_lev:)
   qv(:ncol,top_lev:)     = qv_star(:ncol,top_lev:)
   ql(:ncol,top_lev:)     = ql_star(:ncol,top_lev:)
   qi_mm(:ncol,top_lev:)  = qi_star(:ncol,top_lev:)
   al_st(:ncol,top_lev:)  = alst(:ncol,top_lev:)
   ai_st(:ncol,top_lev:)  = aist(:ncol,top_lev:)
   a_st(:ncol,top_lev:)   = a_st_star(:ncol,top_lev:)
   ql_st(:ncol,top_lev:)  = qlst(:ncol,top_lev:)
   qi_st(:ncol,top_lev:)  = qist(:ncol,top_lev:)
   nl(:ncol,top_lev:)     = nl_star(:ncol,top_lev:)
   ni(:ncol,top_lev:)     = ni_star(:ncol,top_lev:)

   enddo ! End of 'iter' prognostic iterative computation

   ! ------------------------------------------------------------------------ !
   ! Compute final tendencies of main output variables and diagnostic outputs !
   ! Note that the very input state [T0,qv0,ql0,qi0] are                      !
   ! marched to [T_star,qv_star,ql_star,qi_star] with equilibrium             !
   ! stratus informations of [a_st_star,ql_st_star,qi_st_star] by             !
   ! below final tendencies and [A_T,A_qv,A_ql,A_qi]                          !
   ! ------------------------------------------------------------------------ !

   ! ------------------ !
   ! Process tendencies !
   ! ------------------ !

   QQw_final(:ncol,top_lev:)  = QQw_prog(:ncol,top_lev:)
   QQi_final(:ncol,top_lev:)  = QQi_prog(:ncol,top_lev:)
   QQ_final(:ncol,top_lev:)   = QQw_final(:ncol,top_lev:) + QQi_final(:ncol,top_lev:)
   QQw_all(:ncol,top_lev:)    = QQw_prog(:ncol,top_lev:)  + QQw1(:ncol,top_lev:) + QQw2(:ncol,top_lev:) + &
        qlten_pwi1(:ncol,top_lev:) + qlten_pwi2(:ncol,top_lev:) + A_ql_adj(:ncol,top_lev:)
   QQi_all(:ncol,top_lev:)    = QQi_prog(:ncol,top_lev:)  + QQi1(:ncol,top_lev:) + QQi2(:ncol,top_lev:) + &
        qiten_pwi1(:ncol,top_lev:) + qiten_pwi2(:ncol,top_lev:) + A_qi_adj(:ncol,top_lev:)
   QQ_all(:ncol,top_lev:)     = QQw_all(:ncol,top_lev:)   + QQi_all(:ncol,top_lev:)
   QQnl_final(:ncol,top_lev:) = QQnl_prog(:ncol,top_lev:)
   QQni_final(:ncol,top_lev:) = QQni_prog(:ncol,top_lev:)
   QQn_final(:ncol,top_lev:)  = QQnl_final(:ncol,top_lev:) + QQni_final(:ncol,top_lev:)
   QQnl_all(:ncol,top_lev:)   = QQnl_prog(:ncol,top_lev:)  + QQnl1(:ncol,top_lev:) + QQnl2(:ncol,top_lev:) + &
        nlten_pwi1(:ncol,top_lev:) + nlten_pwi2(:ncol,top_lev:) + ACnl(:ncol,top_lev:) + A_nl_adj(:ncol,top_lev:)
   QQni_all(:ncol,top_lev:)   = QQni_prog(:ncol,top_lev:)  + QQni1(:ncol,top_lev:) + QQni2(:ncol,top_lev:) + &
        niten_pwi1(:ncol,top_lev:) + niten_pwi2(:ncol,top_lev:) + ACni(:ncol,top_lev:) + A_ni_adj(:ncol,top_lev:)
   QQn_all(:ncol,top_lev:)    = QQnl_all(:ncol,top_lev:)   + QQni_all(:ncol,top_lev:)
   cmeliq(:ncol,top_lev:)     = QQ_final(:ncol,top_lev:)
   qvadj(:ncol,top_lev:)      = qvten_pwi1(:ncol,top_lev:) + qvten_pwi2(:ncol,top_lev:) + A_qv_adj(:ncol,top_lev:)
   qladj(:ncol,top_lev:)      = qlten_pwi1(:ncol,top_lev:) + qlten_pwi2(:ncol,top_lev:) + A_ql_adj(:ncol,top_lev:)
   qiadj(:ncol,top_lev:)      = qiten_pwi1(:ncol,top_lev:) + qiten_pwi2(:ncol,top_lev:) + A_qi_adj(:ncol,top_lev:)
   qllim(:ncol,top_lev:)      = QQw1      (:ncol,top_lev:) + QQw2      (:ncol,top_lev:)
   qilim(:ncol,top_lev:)      = QQi1      (:ncol,top_lev:) + QQi2      (:ncol,top_lev:)

   ! ----------------- !
   ! Output tendencies !
   ! ----------------- !

   tlat(:ncol,top_lev:)  = cpair*( T_star(:ncol,top_lev:)  -  t_inout(:ncol,top_lev:) )/dtime - &
        cpair*(ttend(:ncol,top_lev:)+CC_T_loc(:ncol,top_lev:))
   qvlat(:ncol,top_lev:) =    ( qv_star(:ncol,top_lev:) - qv_inout(:ncol,top_lev:) )/dtime - &
        (qtend(:ncol,top_lev:)+CC_qv_loc(:ncol,top_lev:))
   qcten(:ncol,top_lev:) =    ( ql_star(:ncol,top_lev:) - ql_inout(:ncol,top_lev:) )/dtime - &
        (lmitend(:ncol,top_lev:)+CC_ql_loc(:ncol,top_lev:))
   qiten(:ncol,top_lev:) =    ( qi_star(:ncol,top_lev:) - qi_inout(:ncol,top_lev:) )/dtime - &
        (itend(:ncol,top_lev:)+CC_qi_loc(:ncol,top_lev:))
   ncten(:ncol,top_lev:) =    ( nl_star(:ncol,top_lev:) - nl_inout(:ncol,top_lev:) )/dtime - &
        (nltend(:ncol,top_lev:)+CC_nl_loc(:ncol,top_lev:))
   niten(:ncol,top_lev:) =    ( ni_star(:ncol,top_lev:) - ni_inout(:ncol,top_lev:) )/dtime - &
        (nitend(:ncol,top_lev:)+CC_ni_loc(:ncol,top_lev:))

   if (.not. do_cldice) then
      do k = top_lev, pver
         do i = 1, ncol

            ! Don't want either qi or ni tendencies, but the code above is somewhat convoluted and
            ! is trying to adjust both (small numbers). Just force it to zero here.
            qiten(i,k) = 0._kind_phys
            niten(i,k) = 0._kind_phys
          end do
      end do
   end if

   ! ------------------ !
   ! Net cloud fraction !
   ! ------------------ !

   cld(:ncol,top_lev:) = a_st_star(:ncol,top_lev:) + concld(:ncol,top_lev:)

   ! --------------------------------- !
   ! Updated grid-mean state variables !
   ! --------------------------------- !

   t_inout(:ncol,top_lev:)  = T_star(:ncol,top_lev:)
   qv_inout(:ncol,top_lev:) = qv_star(:ncol,top_lev:)
   ql_inout(:ncol,top_lev:) = ql_star(:ncol,top_lev:)
   qi_inout(:ncol,top_lev:) = qi_star(:ncol,top_lev:)
   nl_inout(:ncol,top_lev:) = nl_star(:ncol,top_lev:)
   ni_inout(:ncol,top_lev:) = ni_star(:ncol,top_lev:)

   ! --------------------------------------------------------- !
   ! Section 3: Update equilibrium state storage and compute   !
   ! ast from alst and aist                                    !
   ! --------------------------------------------------------- !

   tcwat(:ncol,top_lev:pver)  = T_star(:ncol,top_lev:pver)
   qcwat(:ncol,top_lev:pver)  = qv_star(:ncol,top_lev:pver)
   lcwat(:ncol,top_lev:pver)  = ql_star(:ncol,top_lev:pver) + qi_star(:ncol,top_lev:pver)
   iccwat(:ncol,top_lev:pver) = qi_star(:ncol,top_lev:pver)
   nlwat(:ncol,top_lev:pver)  = nl_star(:ncol,top_lev:pver)
   niwat(:ncol,top_lev:pver)  = ni_star(:ncol,top_lev:pver)

   ! Compute net stratus fraction using maximum over-lapping assumption
   ast(:ncol,:top_lev-1) = 0._kind_phys
   ast(:ncol,top_lev:pver) = max( alst(:ncol,top_lev:pver), aist(:ncol,top_lev:pver) )

   ! Compute default critical RH diagnostics
   rhmin_liq(:ncol,:) = 0._kind_phys
   do k = top_lev, pver
      do i = 1, ncol
         land = nint(landfrac(i)) == 1
         if( pmid(i,k) .ge. premib ) then
            if( land .and. (snowh(i).le.0.000001_kind_phys) ) then
               rhmin_liq(i,k) = rhminl_const - rhminl_adj_land_const
            else
               rhmin_liq(i,k) = rhminl_const
            endif
         elseif( pmid(i,k) .lt. premit ) then
            rhmin_liq(i,k) = rhminh_const
         else
            tmp = (premib-(max(pmid(i,k),premit)))/(premib-premit)
            rhmin_liq(i,k) = rhminh_const*tmp + rhminl_const*(1.0_kind_phys-tmp)
         endif
      end do
   end do

   rhmin_ice(:ncol,:) = rhminh_const

   end subroutine park_macrophysics_run


!=======================================================================================================

   subroutine instratus_condensate( ncol, k,                      &
                                    p_in, T0_in, qv0_in, ql0_in, qi0_in, &
                                    ni0_in,                              &
                                    a_dc_in, ql_dc_in, qi_dc_in,         &
                                    a_sc_in, ql_sc_in, qi_sc_in,         &
                                    landfrac, snowh,                     &
                                    rhmaxi_in, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in, &
                                    T_out, qv_out, ql_out, qi_out,       &
                                    al_st_out, ai_st_out, ql_st_out, qi_st_out )

   ! ------------------------------------------------------- !
   ! Diagnostically force in-stratus condensate to be        !
   ! in the range of 'qlst_min < qc_st < qlst_max'           !
   ! whenever stratus exists in the equilibrium state        !
   ! ------------------------------------------------------- !

   integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
   integer,  intent(in)  :: k                    ! Layer index

   real(kind_phys), intent(in)  :: p_in(:)       ! Pressure [Pa]
   real(kind_phys), intent(in)  :: T0_in(:)      ! Temperature [K]
   real(kind_phys), intent(in)  :: qv0_in(:)     ! Grid-mean water vapor [kg/kg]
   real(kind_phys), intent(in)  :: ql0_in(:)     ! Grid-mean LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi0_in(:)     ! Grid-mean IWC [kg/kg]
   real(kind_phys), intent(in)  :: ni0_in(:)

   real(kind_phys), intent(in)  :: a_dc_in(:)    ! Deep cumulus cloud fraction
   real(kind_phys), intent(in)  :: ql_dc_in(:)   ! In-deep cumulus LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi_dc_in(:)   ! In-deep cumulus IWC [kg/kg]
   real(kind_phys), intent(in)  :: a_sc_in(:)    ! Shallow cumulus cloud fraction
   real(kind_phys), intent(in)  :: ql_sc_in(:)   ! In-shallow cumulus LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi_sc_in(:)   ! In-shallow cumulus IWC [kg/kg]

   real(kind_phys), intent(in)  :: landfrac(:)   ! Land fraction
   real(kind_phys), intent(in)  :: snowh(:)      ! Snow depth (liquid water equivalent)

   real(kind_phys), intent(in)  :: rhmaxi_in(:)
   real(kind_phys), intent(in)  :: rhmini_in(:)
   real(kind_phys), intent(in)  :: rhminl_in(:)
   real(kind_phys), intent(in)  :: rhminl_adj_land_in(:)
   real(kind_phys), intent(in)  :: rhminh_in(:)

   real(kind_phys), intent(out) :: T_out(:)      ! Temperature [K]
   real(kind_phys), intent(out) :: qv_out(:)     ! Grid-mean water vapor [kg/kg]
   real(kind_phys), intent(out) :: ql_out(:)     ! Grid-mean LWC [kg/kg]
   real(kind_phys), intent(out) :: qi_out(:)     ! Grid-mean IWC [kg/kg]

   real(kind_phys), intent(out) :: al_st_out(:)  ! Liquid stratus fraction
   real(kind_phys), intent(out) :: ai_st_out(:)  ! Ice stratus fraction
   real(kind_phys), intent(out) :: ql_st_out(:)  ! In-stratus LWC [kg/kg]
   real(kind_phys), intent(out) :: qi_st_out(:)  ! In-stratus IWC [kg/kg]

   ! Local variables

   integer i                                     ! Column    index

   real(kind_phys) p
   real(kind_phys) T0
   real(kind_phys) qv0
   real(kind_phys) ql0
   real(kind_phys) qi0
   real(kind_phys) a_dc
   real(kind_phys) ql_dc
   real(kind_phys) qi_dc
   real(kind_phys) a_sc
   real(kind_phys) ql_sc
   real(kind_phys) qi_sc
   real(kind_phys) esat0
   real(kind_phys) qsat0
   real(kind_phys) U0
   real(kind_phys) U0_nc
   real(kind_phys) G0_nc
   real(kind_phys) al0_st_nc
   real(kind_phys) al0_st
   real(kind_phys) ai0_st_nc
   real(kind_phys) ai0_st
   real(kind_phys) a0_st
   real(kind_phys) ql0_nc
   real(kind_phys) qi0_nc
   real(kind_phys) qc0_nc
   real(kind_phys) ql0_st
   real(kind_phys) qi0_st
   real(kind_phys) qc0_st
   real(kind_phys) T
   real(kind_phys) qv
   real(kind_phys) ql
   real(kind_phys) qi
   real(kind_phys) ql_st
   real(kind_phys) qi_st
   real(kind_phys) es
   real(kind_phys) qs
   real(kind_phys) esat_in(ncol)
   real(kind_phys) qsat_in(ncol)
   real(kind_phys) U0_in(ncol)
   real(kind_phys) al0_st_nc_in(ncol)
   real(kind_phys) ai0_st_nc_in(ncol)
   real(kind_phys) G0_nc_in(ncol)
   integer  idxmod
   real(kind_phys) U
   real(kind_phys) U_nc
   real(kind_phys) al_st_nc
   real(kind_phys) ai_st_nc
   real(kind_phys) G_nc
   real(kind_phys) a_st
   real(kind_phys) al_st
   real(kind_phys) ai_st
   real(kind_phys) Tmin0
   real(kind_phys) Tmax0
   real(kind_phys) Tmin
   real(kind_phys) Tmax
   integer caseid

   real(kind_phys) rhmaxi
   real(kind_phys) rhmini
   real(kind_phys) rhminl
   real(kind_phys) rhminl_adj_land
   real(kind_phys) rhminh

   ! ---------------- !
   ! Main Computation !
   ! ---------------- !

   call qsat_water(T0_in(1:ncol), p_in(1:ncol), esat_in(1:ncol), qsat_in(1:ncol), ncol)
   U0_in(:ncol) = qv0_in(:ncol)/qsat_in(:ncol)
   if( CAMstfrac ) then
       call astG_RHU(U0_in(:),p_in(:),qv0_in(:),landfrac(:),snowh(:),al0_st_nc_in(:),G0_nc_in(:),ncol,&
                     rhminl_in(:), rhminl_adj_land_in(:), rhminh_in(:))
   else
       call astG_PDF(U0_in(:),p_in(:),qv0_in(:),landfrac(:),snowh(:),al0_st_nc_in(:),G0_nc_in(:),ncol,&
                     rhminl_in(:), rhminl_adj_land_in(:), rhminh_in(:))
   endif
   call aist_vector(qv0_in(:),T0_in(:),p_in(:),qi0_in(:),ni0_in(:),landfrac(:),snowh(:),ai0_st_nc_in(:),ncol,&
                    rhmaxi_in(:), rhmini_in(:), rhminl_in(:), rhminl_adj_land_in(:), rhminh_in(:))

   do i = 1, ncol

      ! ---------------------- !
      ! Define local variables !
      ! ---------------------- !

      p   = p_in(i)

      T0  = T0_in(i)
      qv0 = qv0_in(i)
      ql0 = ql0_in(i)
      qi0 = qi0_in(i)

      a_dc  = a_dc_in(i)
      ql_dc = ql_dc_in(i)
      qi_dc = qi_dc_in(i)

      a_sc  = a_sc_in(i)
      ql_sc = ql_sc_in(i)
      qi_sc = qi_sc_in(i)

      ql_dc = 0._kind_phys
      qi_dc = 0._kind_phys
      ql_sc = 0._kind_phys
      qi_sc = 0._kind_phys

      es  = esat_in(i)
      qs  = qsat_in(i)

      rhmaxi = rhmaxi_in(i)
      rhmini = rhmini_in(i)
      rhminl = rhminl_in(i)
      rhminl_adj_land = rhminl_adj_land_in(i)
      rhminh = rhminh_in(i)

      idxmod = 0
      caseid = -1

      ! ------------------------------------------------------------ !
      ! Force the grid-mean RH to be smaller than 1 if oversaturated !
      ! In order to be compatible with reduced 3x3 QQ, condensation  !
      ! should occur only into the liquid in gridmean_RH.            !
      ! ------------------------------------------------------------ !

      if( qv0 .gt. qs ) then
          call gridmean_RH( i, k, p, T0, qv0, ql0, qi0,      &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                            landfrac(i), snowh(i) )
          call qsat_water(T0, p, esat0, qsat0)
          U0      = (qv0/qsat0)
          U0_nc   =  U0
          if( CAMstfrac ) then
              call astG_RHU_single(U0_nc, p, qv0, landfrac(i), snowh(i), al0_st_nc, G0_nc, &
                 rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
          else
              call astG_PDF_single(U0_nc, p, qv0, landfrac(i), snowh(i), al0_st_nc, G0_nc, &
                 rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
          endif
          call aist_single(qv0,T0,p,qi0,landfrac(i),snowh(i),ai0_st_nc,&
                           rhmaxi, rhmini, rhminl, rhminl_adj_land, rhminh)
          ai0_st  = (1._kind_phys-a_dc-a_sc)*ai0_st_nc
          al0_st  = (1._kind_phys-a_dc-a_sc)*al0_st_nc
          a0_st   = max(ai0_st,al0_st)
          idxmod  = 1
      else
          ai0_st  = (1._kind_phys-a_dc-a_sc)*ai0_st_nc_in(i)
          al0_st  = (1._kind_phys-a_dc-a_sc)*al0_st_nc_in(i)
      endif
      a0_st   = max(ai0_st,al0_st)

      ! ----------------------- !
      ! Handling of input state !
      ! ----------------------- !

      ql0_nc  = max(0._kind_phys,ql0-a_dc*ql_dc-a_sc*ql_sc)
      qi0_nc  = max(0._kind_phys,qi0-a_dc*qi_dc-a_sc*qi_sc)
      qc0_nc  = ql0_nc + qi0_nc

      Tmin0 = T0 - (latvap/cpair)*ql0
      Tmax0 = T0 + ((latvap+latice)/cpair)*qv0

      ! ------------------------------------------------------------- !
      ! Do nothing and just exit if generalized in-stratus condensate !
      ! condition is satisfied. This includes the case I.             !
      ! For 4x4 liquid stratus, a0_st --> al0_st.                     !
      ! ------------------------------------------------------------- !
      if( ( ql0_nc .ge. qlst_min*al0_st ) .and. ( ql0_nc .le. qlst_max*al0_st ) ) then

          ! ------------------ !
          ! This is the case I !
          ! ------------------ !
             T = T0
             qv = qv0
             ql = ql0
             qi = qi0
             caseid = 0
             goto 10
      else
         ! ----------------------------- !
         ! This is case II : Dense Cloud !
         ! ----------------------------- !
         if( al0_st .eq. 0._kind_phys .and. ql0_nc .gt. 0._kind_phys ) then
             ! ------------------------------------- !
             ! Compute hypothetical full evaporation !
             ! ------------------------------------- !
             T  = Tmin0
             qv = qv0 + ql0
             call qsat_water(T, p, es, qs)
             U  = qv/qs
             U_nc = U
             if( CAMstfrac ) then
                 call astG_RHU_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                    rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
             else
                 call astG_PDF_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                    rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
             endif
             al_st = (1._kind_phys-a_dc-a_sc)*al_st_nc
             caseid = 0

             if( al_st .eq. 0._kind_phys ) then
                 ql = 0._kind_phys
                 qi = qi0
                 idxmod = 1
                 caseid = 1
                 goto 10
             else
                 ! ------------------------------------------- !
                 ! Evaporate until qc_st decreases to qlst_max !
                 ! ------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = T0
                 call instratus_core( i, k, p,                              &
                                      T0, qv0, ql0, 0._kind_phys,          &
                                      a_dc, ql_dc, qi_dc,                  &
                                      a_sc, ql_sc, qi_sc, ai0_st,          &
                                      qlst_max, Tmin, Tmax, landfrac(i), snowh(i), &
                                      rhminl, rhminl_adj_land, rhminh,     &
                                      T, qv, ql, qi )
                 idxmod = 1
                 caseid = 2
                 goto 10
             endif
         ! ------------------------------ !
         ! This is case III : Empty Cloud !
         ! ------------------------------ !
         elseif( al0_st .gt. 0._kind_phys .and. ql0_nc .eq. 0._kind_phys ) then
              ! ------------------------------------------ !
              ! Condense until qc_st increases to qlst_min !
              ! ------------------------------------------ !
              Tmin = Tmin0
              Tmax = Tmax0
              call instratus_core( i, k, p,                              &
                                   T0, qv0, ql0, 0._kind_phys,          &
                                   a_dc, ql_dc, qi_dc,                  &
                                   a_sc, ql_sc, qi_sc, ai0_st,          &
                                   qlst_min, Tmin, Tmax, landfrac(i), snowh(i), &
                                   rhminl, rhminl_adj_land, rhminh,     &
                                   T, qv, ql, qi )
              idxmod = 1
              caseid = 3
              goto 10
         ! --------------- !
         ! This is case IV !
         ! --------------- !
         elseif( al0_st .gt. 0._kind_phys .and. ql0_nc .gt. 0._kind_phys ) then

             if( ql0_nc .gt. qlst_max*al0_st ) then
                 ! --------------------------------------- !
                 ! Evaporate until qc_st drops to qlst_max !
                 ! --------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0
                 call instratus_core( i, k, p,                              &
                                      T0, qv0, ql0, 0._kind_phys,          &
                                      a_dc, ql_dc, qi_dc,                  &
                                      a_sc, ql_sc, qi_sc, ai0_st,          &
                                      qlst_max, Tmin, Tmax, landfrac(i), snowh(i), &
                                      rhminl, rhminl_adj_land, rhminh,     &
                                      T, qv, ql, qi )
                 idxmod = 1
                 caseid = 4
                 goto 10
             elseif( ql0_nc .lt. qlst_min*al0_st ) then
                 ! -------------------------------------------- !
                 ! Condensate until qc_st increases to qlst_min !
                 ! -------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0
                 call instratus_core( i, k, p,                              &
                                      T0, qv0, ql0, 0._kind_phys,          &
                                      a_dc, ql_dc, qi_dc,                  &
                                      a_sc, ql_sc, qi_sc, ai0_st,          &
                                      qlst_min, Tmin, Tmax, landfrac(i), snowh(i), &
                                      rhminl, rhminl_adj_land, rhminh,     &
                                      T, qv, ql, qi )
                 idxmod = 1
                 caseid = 5
                 goto 10
             else
                 ! ------------------------------------------------ !
                 ! This case should not happen. Issue error message !
                 ! ------------------------------------------------ !
                 write(iulog,*) 'Impossible case1 in instratus_condensate'
                 return
             endif
         ! ------------------------------------------------ !
         ! This case should not happen. Issue error message !
         ! ------------------------------------------------ !
         else
             write(iulog,*) 'Impossible case2 in instratus_condensate'
             write(iulog,*)  al0_st, a_sc, a_dc
             write(iulog,*)  1000*ql0_nc, 1000*(ql0+qi0)
             return
         endif
      endif

10 continue

   ! -------------------------------------------------- !
   ! Force final energy-moisture conserving consistency !
   ! -------------------------------------------------- !

     qi = qi0

     if( idxmod .eq. 1 ) then
         call aist_single(qv,T,p,qi,landfrac(i),snowh(i),ai_st_nc,&
                          rhmaxi, rhmini, rhminl, rhminl_adj_land, rhminh)
         ai_st = (1._kind_phys-a_dc-a_sc)*ai_st_nc
         call qsat_water(T, p, es, qs)
         U     = (qv/qs)
         U_nc  =  U
         if( CAMstfrac ) then
             call astG_RHU_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
         else
             call astG_PDF_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
         endif
         al_st = (1._kind_phys-a_dc-a_sc)*al_st_nc
     else
         ai_st  = (1._kind_phys-a_dc-a_sc)*ai0_st_nc_in(i)
         al_st  = (1._kind_phys-a_dc-a_sc)*al0_st_nc_in(i)
     endif

     a_st  = max(ai_st,al_st)

     if( al_st .eq. 0._kind_phys ) then
         ql_st = 0._kind_phys
     else
         ql_st = ql/al_st
         ql_st = min(qlst_max,max(qlst_min,ql_st)) ! PJR
     endif
     if( ai_st .eq. 0._kind_phys ) then
         qi_st = 0._kind_phys
     else
         qi_st = qi/ai_st
     endif

     qi    = ai_st*qi_st
     ql    = al_st*ql_st

     T     = T0 - (latvap/cpair)*(ql0-ql) - ((latvap+latice)/cpair)*(qi0-qi)
     qv    = qv0 + ql0 - ql + qi0 - qi

   ! -------------- !
   ! Send to output !
   ! -------------- !

   T_out(i)  = T
   qv_out(i) = qv
   ql_out(i) = ql
   qi_out(i) = qi
   al_st_out(i) = al_st
   ai_st_out(i) = ai_st
   ql_st_out(i) = ql_st
   qi_st_out(i) = qi_st

   enddo

   return
   end subroutine instratus_condensate

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine instratus_core( icol, k, p,                      &
                              T0, qv0, ql0, qi0,                      &
                              a_dc, ql_dc, qi_dc,                     &
                              a_sc, ql_sc, qi_sc, ai_st,              &
                              qcst_crit, Tmin, Tmax, landfrac, snowh, &
                              rhminl, rhminl_adj_land, rhminh,        &
                              T, qv, ql, qi )

   ! ------------------------------------------------------ !
   ! Subroutine to find saturation equilibrium state using  !
   ! a Newton iteration method, so that 'qc_st = qcst_crit' !
   ! is satisfied.                                          !
   ! ------------------------------------------------------ !

   integer,  intent(in)  :: icol       ! Number of atmospheric columns
   integer,  intent(in)  :: k          ! Layer index

   real(kind_phys), intent(in)  :: p          ! Pressure [Pa]
   real(kind_phys), intent(in)  :: T0         ! Temperature [K]
   real(kind_phys), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(kind_phys), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]

   real(kind_phys), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(kind_phys), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(kind_phys), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(kind_phys), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(kind_phys), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(kind_phys), intent(in)  :: Tmin       ! Minimum temperature system can have [K]
   real(kind_phys), intent(in)  :: Tmax       ! Maximum temperature system can have [K]
   real(kind_phys), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(kind_phys), intent(in)  :: landfrac   ! Land fraction
   real(kind_phys), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(kind_phys), intent(in)  :: rhminl
   real(kind_phys), intent(in)  :: rhminl_adj_land
   real(kind_phys), intent(in)  :: rhminh

   real(kind_phys), intent(out) :: T          ! Temperature [K]
   real(kind_phys), intent(out) :: qv         ! Grid-mean water vapor [kg/kg]
   real(kind_phys), intent(out) :: ql         ! Grid-mean LWC [kg/kg]
   real(kind_phys), intent(out) :: qi         ! Grid-mean IWC [kg/kg]

   ! Local variables

   integer i                           ! Iteration index

   real(kind_phys) muQ0, muQ
   real(kind_phys) ql_nc0, qi_nc0, qc_nc0, qc_nc
   real(kind_phys) fice0, fice
   real(kind_phys) ficeg0, ficeg
   real(kind_phys) esat0
   real(kind_phys) qsat0
   real(kind_phys) dqcncdt, dastdt, dUdt
   real(kind_phys) alpha, beta
   real(kind_phys) U, U_nc
   real(kind_phys) al_st_nc, G_nc
   real(kind_phys) al_st

   ! Variables for root-finding algorithm

   integer j
   real(kind_phys)  x1, x2
   real(kind_phys)  rtsafe
   real(kind_phys)  df, dx, dxold, f, fh, fl, temp, xh, xl
   real(kind_phys), parameter :: xacc = 1.e-3_kind_phys

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0._kind_phys,ql0-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0._kind_phys,qi0-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0._kind_phys,ql0+qi0-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   fice0  = 0._kind_phys
   ficeg0 = 0._kind_phys
   muQ0   = 1._kind_phys

   ! ------------ !
   ! Root finding !
   ! ------------ !

   x1 = Tmin
   x2 = Tmax
   call funcd_instratus( x1, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         rhminl, rhminl_adj_land, rhminh,               &
                         fl, df, qc_nc, fice, al_st )
   call funcd_instratus( x2, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         rhminl, rhminl_adj_land, rhminh,               &
                         fh, df, qc_nc, fice, al_st )
   if((fl > 0._kind_phys .and. fh > 0._kind_phys) .or. (fl < 0._kind_phys .and. fh < 0._kind_phys)) then
       call funcd_instratus( T0, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                             a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                             qcst_crit, landfrac, snowh,                    &
                             rhminl, rhminl_adj_land, rhminh,               &
                             fl, df, qc_nc, fice, al_st )
       rtsafe = T0
       goto 10
   endif
   if( fl == 0._kind_phys) then
           rtsafe = x1
           goto 10
   elseif( fh == 0._kind_phys) then
           rtsafe = x2
           goto 10
   elseif( fl < 0._kind_phys) then
           xl = x1
           xh = x2
   else
           xh = x1
           xl = x2
   end if
   rtsafe = 0.5_kind_phys*(x1+x2)
   dxold = abs(x2-x1)
   dx = dxold
   call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                         qcst_crit, landfrac, snowh,                        &
                         rhminl, rhminl_adj_land, rhminh,                   &
                         f, df, qc_nc, fice, al_st )
   do j = 1, 20
      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0._kind_phys .or. abs(2.0_kind_phys*f) > abs(dxold*df) ) then
           dxold = dx
           dx = 0.5_kind_phys*(xh-xl)
           rtsafe = xl + dx
           if(xl == rtsafe) goto 10
      else
           dxold = dx
           dx = f/df
           temp = rtsafe
           rtsafe = rtsafe - dx
           if (temp == rtsafe) goto 10
      end if
    ! if(abs(dx) < xacc) goto 10
      call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                            qcst_crit, landfrac, snowh,                        &
                            rhminl, rhminl_adj_land, rhminh,                   &
                            f, df, qc_nc, fice, al_st )
    ! Sep.21.2010. Sungsu modified to enhance convergence and guarantee 'qlst_min <  qlst < qlst_max'.
      if( qcst_crit < 0.5_kind_phys * ( qlst_min + qlst_max ) ) then
          if( ( qc_nc*(1._kind_phys-fice) .gt.          qlst_min*al_st .and. &
                qc_nc*(1._kind_phys-fice) .lt. 1.1_kind_phys * qlst_min*al_st ) ) goto 10
      else
          if( ( qc_nc*(1._kind_phys-fice) .gt. 0.9_kind_phys * qlst_max*al_st .and. &
                qc_nc*(1._kind_phys-fice) .lt.          qlst_max*al_st ) ) goto 10
      endif
      if(f < 0._kind_phys) then
          xl = rtsafe
      else
          xh = rtsafe
      endif

   enddo

10 continue

   ! ------------------------------------------- !
   ! Final safety check before sending to output !
   ! ------------------------------------------- !

   qc_nc = max(0._kind_phys,qc_nc)

   T  = rtsafe
   ql = qc_nc*(1._kind_phys-fice) + a_dc*ql_dc + a_sc*ql_sc
   qi = qc_nc*fice + a_dc*qi_dc + a_sc*qi_sc
   qv = qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc))
   qv = max(qv,1.e-12_kind_phys)

   return
   end subroutine instratus_core

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine funcd_instratus( T, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0,   &
                               a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,  &
                               qcst_crit, landfrac, snowh,                     &
                               rhminl, rhminl_adj_land, rhminh,                &
                               f, fg, qc_nc, fice, al_st )

   ! --------------------------------------------------- !
   ! Subroutine to find function value and gradient at T !
   ! --------------------------------------------------- !

   implicit none

   real(kind_phys), intent(in)  :: T          ! Iteration temperature [K]

   real(kind_phys), intent(in)  :: p          ! Pressure [Pa]
   real(kind_phys), intent(in)  :: T0         ! Initial temperature [K]
   real(kind_phys), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(kind_phys), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]
   real(kind_phys), intent(in)  :: fice0      !
   real(kind_phys), intent(in)  :: muQ0       !
   real(kind_phys), intent(in)  :: qc_nc0     !

   real(kind_phys), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(kind_phys), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(kind_phys), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(kind_phys), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(kind_phys), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(kind_phys), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(kind_phys), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(kind_phys), intent(in)  :: landfrac   ! Land fraction
   real(kind_phys), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(kind_phys), intent(in)  :: rhminl
   real(kind_phys), intent(in)  :: rhminl_adj_land
   real(kind_phys), intent(in)  :: rhminh

   real(kind_phys), intent(out) :: f          ! Value of minimization function at T
   real(kind_phys), intent(out) :: fg         ! Gradient of minimization function
   real(kind_phys), intent(out) :: qc_nc      !
   real(kind_phys), intent(out) :: al_st      !
   real(kind_phys), intent(out) :: fice       !

   ! Local variables

   real(kind_phys) es
   real(kind_phys) qs
   real(kind_phys) dqsdT
   real(kind_phys) dqcncdt
   real(kind_phys) alpha
   real(kind_phys) beta
   real(kind_phys) U
   real(kind_phys) U_nc
   real(kind_phys) al_st_nc
   real(kind_phys) G_nc
   real(kind_phys) dUdt
   real(kind_phys) dalstdt
   real(kind_phys) qv

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   call qsat_water(T, p, es, qs, dqsdt=dqsdT)

   fice    = fice0
   qc_nc   = (cpair/latvap)*(T-T0)+muQ0*qc_nc0
   dqcncdt = (cpair/latvap)
   qv      = (qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc)))
   alpha   = (1._kind_phys/qs)
   beta    = (qv/qs**2._kind_phys)*dqsdT

   U      =  (qv/qs)
   U_nc   =   U
   if( CAMstfrac ) then
       call astG_RHU_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
          rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
   else
       call astG_PDF_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
          rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
   endif
   al_st   =  (1._kind_phys-a_dc-a_sc)*al_st_nc
   dUdt    = -(alpha*dqcncdt+beta)
   dalstdt =  (1._kind_phys/G_nc)*dUdt
   if( U_nc .eq. 1._kind_phys ) dalstdt = 0._kind_phys

   f  = qc_nc   - qcst_crit*al_st
   fg = dqcncdt - qcst_crit*dalstdt

   return
   end subroutine funcd_instratus

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine gridmean_RH( icol, k, p, T, qv, ql, qi,       &
                           a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                           landfrac, snowh )

   ! ------------------------------------------------------------- !
   ! Subroutine to force grid-mean RH = 1 when RH > 1              !
   ! This is condensation process similar to instratus_condensate. !
   ! During condensation, we assume 'fice' is maintained in this   !
   ! verison for MG not for RK.                                    !
   ! ------------------------------------------------------------- !

   integer,  intent(in)    :: icol       ! Number of atmospheric columns
   integer,  intent(in)    :: k          ! Layer index

   real(kind_phys), intent(in)    :: p          ! Pressure [Pa]
   real(kind_phys), intent(inout) :: T          ! Temperature [K]
   real(kind_phys), intent(inout) :: qv         ! Grid-mean water vapor [kg/kg]
   real(kind_phys), intent(inout) :: ql         ! Grid-mean LWC [kg/kg]
   real(kind_phys), intent(inout) :: qi         ! Grid-mean IWC [kg/kg]

   real(kind_phys), intent(in)    :: a_dc       ! Deep cumulus cloud fraction
   real(kind_phys), intent(in)    :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(kind_phys), intent(in)    :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(kind_phys), intent(in)    :: a_sc       ! Shallow cumulus cloud fraction
   real(kind_phys), intent(in)    :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(kind_phys), intent(in)    :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(kind_phys), intent(in)    :: landfrac   ! Land fraction
   real(kind_phys), intent(in)    :: snowh      ! Snow depth (liquid water equivalent)

   ! Local variables

   integer m                             ! Iteration index

   real(kind_phys)  ql_nc0, qi_nc0, qc_nc0
   real(kind_phys)  Tscale
   real(kind_phys)  Tc, qt, qc, dqcdt, qc_nc
   real(kind_phys)  es, qs, dqsdT
   real(kind_phys)  al_st_nc, G_nc
   real(kind_phys)  f, fg
   real(kind_phys), parameter :: xacc = 1.e-3_kind_phys

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0._kind_phys,ql-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0._kind_phys,qi-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0._kind_phys,ql+qi-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   Tc    = T - (latvap/cpair)*ql
   qt    = qv + ql

   do m = 1, 20
      call qsat_water(T, p, es, qs, dqsdt=dqsdT)
      Tscale = latvap/cpair
      qc     = (T-Tc)/Tscale
      dqcdt  = 1._kind_phys/Tscale
      f      = qs + qc - qt
      fg     = dqsdT + dqcdt
      fg     = sign(1._kind_phys,fg)*max(1.e-10_kind_phys,abs(fg))
    ! Sungsu modified convergence criteria to speed up convergence and guarantee RH <= 1.
      if( qc .ge. 0._kind_phys .and. ( qt - qc ) .ge. 0.999_kind_phys*qs .and. ( qt - qc ) .le. 1._kind_phys*qs ) then
          goto 10
      endif
      T = T - f/fg
   enddo
 ! write(iulog,*) 'Convergence in gridmean_RH is not reached. RH = ', ( qt - qc ) / qs
10 continue

   call qsat_water(T, p, es, qs)
 ! Sungsu modified 'qv = qs' in consistent with the modified convergence criteria above.
   qv = min(qt,qs) ! Modified
   ql = qt - qv
   T  = Tc + (latvap/cpair)*ql

   return
   end subroutine gridmean_RH

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine positive_moisture( ncol, dt, qvmin, qlmin, qimin, dp, &
                                 qv,   ql, qi,    t,     qvten, &
                                 qlten,    qiten, tten,  do_cldice)

   ! ------------------------------------------------------------------------------- !
   ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
   ! force them to be larger than minimum value by (1) condensating water vapor      !
   ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
   ! layer. '2._kind_phys' is multiplied to the minimum values for safety.           !
   ! Update final state variables and tendencies associated with this correction.    !
   ! If any condensation happens, update (s,t) too.                                  !
   ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': top layer, 'pver' : near-surface layer         !
   ! ------------------------------------------------------------------------------- !

   implicit none
   integer,  intent(in)     :: ncol
   real(kind_phys), intent(in)     :: dt
   real(kind_phys), intent(in)     :: dp(:,:), qvmin(:,:), qlmin(:,:), qimin(:,:)
   real(kind_phys), intent(inout)  :: qv(:,:), ql(:,:), qi(:,:), t(:,:)
   real(kind_phys), intent(out)    :: qvten(:,:), qlten(:,:), qiten(:,:), tten(:,:)
   logical, intent(in)      :: do_cldice
   integer   i, k
   real(kind_phys)  dql, dqi, dqv, sum, aa, dum
   integer :: pver_loc

   pver_loc = size(dp, 2)

   tten(:ncol,:pver_loc)  = 0._kind_phys
   qvten(:ncol,:pver_loc) = 0._kind_phys
   qlten(:ncol,:pver_loc) = 0._kind_phys
   qiten(:ncol,:pver_loc) = 0._kind_phys

   do i = 1, ncol
      do k = top_lev_mod, pver_loc
         if( qv(i,k) .lt. qvmin(i,k) .or. ql(i,k) .lt. qlmin(i,k) .or. qi(i,k) .lt. qimin(i,k) ) then
             goto 10
         endif
      enddo
      goto 11
   10 continue
      do k = top_lev_mod, pver_loc    ! From the top to the 1st (lowest) layer from the surface
         dql = max(0._kind_phys,1._kind_phys*qlmin(i,k)-ql(i,k))

         if (do_cldice) then
         dqi = max(0._kind_phys,1._kind_phys*qimin(i,k)-qi(i,k))
         else
           dqi = 0._kind_phys
         end if

         qlten(i,k) = qlten(i,k) +  dql/dt
         qiten(i,k) = qiten(i,k) +  dqi/dt
         qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
         tten(i,k)  = tten(i,k)  + (latvap/cpair)*(dql/dt) + ((latvap+latice)/cpair)*(dqi/dt)
         ql(i,k)    = ql(i,k) + dql
         qi(i,k)    = qi(i,k) + dqi
         qv(i,k)    = qv(i,k) - dql - dqi
         t(i,k)     = t(i,k)  + (latvap * dql + (latvap+latice) * dqi)/cpair
         dqv        = max(0._kind_phys,1._kind_phys*qvmin(i,k)-qv(i,k))
         qvten(i,k) = qvten(i,k) + dqv/dt
         qv(i,k)    = qv(i,k)    + dqv
         if( k .ne. pver_loc ) then
             qv(i,k+1)    = qv(i,k+1)    - dqv*dp(i,k)/dp(i,k+1)
             qvten(i,k+1) = qvten(i,k+1) - dqv*dp(i,k)/dp(i,k+1)/dt
         endif
         qv(i,k) = max(qv(i,k),qvmin(i,k))
         ql(i,k) = max(ql(i,k),qlmin(i,k))
         qi(i,k) = max(qi(i,k),qimin(i,k))
      end do
      ! Extra moisture used to satisfy 'qv(i,pver_loc)=qvmin' is proportionally
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture.
      if( dqv .gt. 1.e-20_kind_phys ) then
          sum = 0._kind_phys
          do k = top_lev_mod, pver_loc
             if( qv(i,k) .gt. 2._kind_phys*qvmin(i,k) ) sum = sum + qv(i,k)*dp(i,k)
          enddo
          aa = dqv*dp(i,pver_loc)/max(1.e-20_kind_phys,sum)
          if( aa .lt. 0.5_kind_phys ) then
              do k = top_lev_mod, pver_loc
                 if( qv(i,k) .gt. 2._kind_phys*qvmin(i,k) ) then
                     dum        = aa*qv(i,k)
                     qv(i,k)    = qv(i,k) - dum
                     qvten(i,k) = qvten(i,k) - dum/dt
                 endif
              enddo
          else
              write(iulog,*) 'Full positive_moisture is impossible in Park Macro'
          endif
      endif
11 continue
   enddo
   return

   end subroutine positive_moisture

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      real(kind_phys) a(np,np),b(np,mp)
      real(kind_phys) aa(np,np),bb(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,ii,jj,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      real(kind_phys) big,dum,pivinv

      aa(:,:) = a(:,:)
      bb(:,:) = b(:,:)

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0._kind_phys
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(iulog,*) 'singular matrix in gaussj 1'
                do ii = 1, np
                do jj = 1, np
                   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
                end do
                end do
                return
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0._kind_phys) then
            write(iulog,*) 'singular matrix in gaussj 2'
            do ii = 1, np
            do jj = 1, np
               write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
            end do
            end do
            return
        endif
        pivinv=1._kind_phys/a(icol,icol)
        a(icol,icol)=1._kind_phys
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0._kind_phys
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue

      return
      end subroutine gaussj

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

end module park_macrophysics
