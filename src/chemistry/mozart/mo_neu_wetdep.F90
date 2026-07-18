!
! code written by J.-F. Lamarque, S. Walters and F. Vitt
! based on the original code from J. Neu developed for UC Irvine
! model
!
! LKE 2/23/2018 - correct setting flag for mass-limited (HNO3,etc.) vs Henry's Law washout
! RPF 9/18/2024 - R. Fernandez - Merge vsl03 chemistry (AC2-CSIC-Madrid - A. Saiz-Lopez) ! rpf_CESM2_SLH
!
module mo_neu_wetdep
!
  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_logfile,      only : iulog
  use constituents,     only : pcnst
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun
  use shr_drydep_mod,   only : n_species_table, species_name_table, dheff
  use gas_wetdep_opts,  only : gas_wetdep_method, gas_wetdep_list, gas_wetdep_cnt, gas_wetdep_ice_uptake_list
  use gas_wetdep_neu,   only : gas_wetdep_neu_init, gas_wetdep_neu_run, do_neu_wetdep, do_diag
!
  implicit none
!
  private
  public :: neu_wetdep_init
  public :: neu_wetdep_tend
!
  save
!
  integer, allocatable, dimension(:) :: mapping_to_mmr
  real(r8),allocatable, dimension(:) :: mol_weight
  integer                     :: index_cldice,index_cldliq
  logical, parameter          :: debug   = .false.
!
contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine neu_wetdep_init
!
  use constituents, only : cnst_get_ind,cnst_mw
  use cam_history,  only : addfld, add_default, horiz_only
  use phys_control, only : phys_getopts, cam_chempkg_is
!
  integer :: m

  logical :: history_chemistry
  character(len=512) :: errmsg
  integer :: errflg

  call phys_getopts(history_chemistry_out=history_chemistry)
!
! initialize the portable core: it captures the wet deposition list and
! Henry's law table and resolves the heff mapping and ice-uptake flags.
! The shr_drydep tables are unallocated when no dep_data_file was read;
! in that case Neu wet deposition cannot run (do_neu_wetdep stays false).
!
  if ( gas_wetdep_method == 'NEU' .and. gas_wetdep_cnt > 0 ) then
    if ( .not. allocated(dheff) ) then
      call endrun('neu_wetdep_init: dep_data_file tables are not available (drydep_inparm namelist not set?)')
    end if
    call gas_wetdep_neu_init( gas_wetdep_method             = gas_wetdep_method, &
                              gas_wetdep_cnt_in             = gas_wetdep_cnt, &
                              pcnst_in                      = pcnst, &
                              gas_wetdep_list_in            = gas_wetdep_list, &
                              gas_wetdep_ice_uptake_list_in = gas_wetdep_ice_uptake_list, &
                              n_species_table               = n_species_table, &
                              species_name_table            = species_name_table, &
                              dheff_in                      = dheff, &
                              is_geoschem_mam4              = cam_chempkg_is('geoschem_mam4'), &
                              errmsg                        = errmsg, &
                              errflg                        = errflg )
    if (errflg /= 0) then
      call endrun(trim(errmsg))
    end if
  end if
!
  if (.not.do_neu_wetdep) return
!
  allocate( mapping_to_mmr(gas_wetdep_cnt) )
  allocate( mol_weight(gas_wetdep_cnt) )
!
! find mapping to species
!
  mapping_to_mmr = -99
  do m=1,gas_wetdep_cnt
    if ( debug .and. masterproc ) write(iulog, '(i4,a)') m,trim(gas_wetdep_list(m))
    call cnst_get_ind(gas_wetdep_list(m), mapping_to_mmr(m), abort=.false. )
    if ( debug .and. masterproc ) write(iulog, '(a,i4)') 'mapping_to_mmr ',mapping_to_mmr(m)
    if ( mapping_to_mmr(m) <= 0 ) then
      if (masterproc) write(iulog,*) 'problem with mapping_to_mmr of ',gas_wetdep_list(m)
      call endrun('neu_wetdep_init: problem with mapping_to_mmr of '//trim(gas_wetdep_list(m)))
    end if
  end do
!
! define species-dependent arrays
!
  do m=1,gas_wetdep_cnt
    mol_weight(m) = cnst_mw(mapping_to_mmr(m))
    if ( debug .and. masterproc ) write(iulog, '(i4,a,f8.4)') m,' mol_weight ',mol_weight(m)
  end do
!
! indices for cloud quantities
!
  call cnst_get_ind( 'CLDICE', index_cldice )
  call cnst_get_ind( 'CLDLIQ', index_cldliq )
!
! define output
!
  do m=1,gas_wetdep_cnt
    call addfld     ('DTWR_'//trim(gas_wetdep_list(m)),(/ 'lev' /), 'A','kg/kg/s','wet removal Neu scheme tendency')
    call addfld     ('WD_'//trim(gas_wetdep_list(m)),horiz_only, 'A','kg/m2/s','vertical integrated wet deposition flux')
    call addfld     ('HEFF_'//trim(gas_wetdep_list(m)),(/ 'lev' /), 'A','M/atm','Effective Henrys Law coeff.')
    if (history_chemistry) then
       call add_default('WD_'//trim(gas_wetdep_list(m)), 1, ' ')
    end if
  end do
!
  if ( do_diag ) then
    call addfld     ('QT_RAIN_HNO3',(/ 'lev' /), 'A','mol/mol/s','wet removal Neu scheme rain tendency')
    call addfld     ('QT_RIME_HNO3',(/ 'lev' /), 'A','mol/mol/s','wet removal Neu scheme rain tendency')
    call addfld     ('QT_WASH_HNO3',(/ 'lev' /), 'A','mol/mol/s','wet removal Neu scheme rain tendency')
    call addfld     ('QT_EVAP_HNO3',(/ 'lev' /), 'A','mol/mol/s','wet removal Neu scheme rain tendency')
    if (history_chemistry) then
       call add_default('QT_RAIN_HNO3',1,' ')
       call add_default('QT_RIME_HNO3',1,' ')
       call add_default('QT_WASH_HNO3',1,' ')
       call add_default('QT_EVAP_HNO3',1,' ')
    end if
  end if
!
  return
!
end subroutine neu_wetdep_init
!
subroutine neu_wetdep_tend(lchnk,ncol,mmr,pmid,pdel,zint,tfld,delt, &
     prain, nevapr, cld, cmfdqr, wd_tend, wd_tend_int)
!
! CAM marshal for the portable core in gas_wetdep_neu.F90: provides cell
! area and latitude from phys_grid, calls gas_wetdep_neu_run, and writes
! the history output.
!
  use ppgrid,           only : pcols, pver
  use phys_grid,        only : get_area_all_p, get_rlat_all_p
  use shr_const_mod,    only : SHR_CONST_REARTH
  use cam_history,      only : outfld
!
  implicit none
!
  integer,        intent(in)    :: lchnk,ncol
  real(r8),       intent(in)    :: mmr(pcols,pver,pcnst)    ! mass mixing ratio (kg/kg)
  real(r8),       intent(in)    :: pmid(pcols,pver)         ! midpoint pressures (Pa)
  real(r8),       intent(in)    :: pdel(pcols,pver)         ! pressure delta about midpoints (Pa)
  real(r8),       intent(in)    :: zint(pcols,pver+1)       ! interface geopotential height above the surface (m)
  real(r8),       intent(in)    :: tfld(pcols,pver)         ! midpoint temperature (K)
  real(r8),       intent(in)    :: delt                     ! timestep (s)
!
  real(r8),       intent(in)    :: prain(ncol, pver)
  real(r8),       intent(in)    :: nevapr(ncol, pver)
  real(r8),       intent(in)    :: cld(ncol, pver)
  real(r8),       intent(in)    :: cmfdqr(ncol, pver)
  real(r8),       intent(inout) :: wd_tend(pcols,pver,pcnst)
  real(r8),       intent(inout) :: wd_tend_int(pcols,pcnst)
!
! local arrays and variables
!
  integer :: m
  real(r8), parameter                       :: rearth = SHR_CONST_REARTH    ! radius earth (m)
  real(r8), dimension(ncol)                 :: area
  real(r8) :: lats(pcols)
  real(r8), dimension(ncol,pver,gas_wetdep_cnt) :: dtwr,heff
  real(r8), dimension(ncol,pver)            :: qt_rain,qt_rime,qt_wash,qt_evap
!
  if (.not.do_neu_wetdep) return
!
! don't do anything if there are no species to be removed
!
  if ( gas_wetdep_cnt == 0 ) return
!
! get area (in radians square)
!
  call get_area_all_p(lchnk, ncol, area)
  area = area * rearth**2                     ! in m^2
  call get_rlat_all_p(lchnk, pcols, lats )
!
  call gas_wetdep_neu_run( ncol           = ncol, &
                           pver           = pver, &
                           mmr            = mmr, &
                           pmid           = pmid, &
                           pdel           = pdel, &
                           zint           = zint, &
                           tfld           = tfld, &
                           delt           = delt, &
                           prain          = prain, &
                           nevapr         = nevapr, &
                           cld            = cld, &
                           cmfdqr         = cmfdqr, &
                           area           = area, &
                           lats           = lats, &
                           mapping_to_mmr = mapping_to_mmr, &
                           mol_weight     = mol_weight, &
                           index_cldice   = index_cldice, &
                           index_cldliq   = index_cldliq, &
                           wd_tend        = wd_tend, &
                           wd_tend_int    = wd_tend_int, &
                           dtwr           = dtwr, &
                           heff           = heff, &
                           qt_rain        = qt_rain, &
                           qt_rime        = qt_rime, &
                           qt_wash        = qt_wash, &
                           qt_evap        = qt_evap )
!
! output tendencies
!
  do m=1,gas_wetdep_cnt
    call outfld( 'DTWR_'//trim(gas_wetdep_list(m)),dtwr(:,:,m),ncol,lchnk )
    call outfld( 'HEFF_'//trim(gas_wetdep_list(m)),heff(:,pver:1:-1,m),ncol,lchnk )
    call outfld( 'WD_'//trim(gas_wetdep_list(m)),wd_tend_int(:ncol,mapping_to_mmr(m)),ncol,lchnk )
  end do
!
  if ( do_diag ) then
    call outfld('QT_RAIN_HNO3', qt_rain, ncol, lchnk )
    call outfld('QT_RIME_HNO3', qt_rime, ncol, lchnk )
    call outfld('QT_WASH_HNO3', qt_wash, ncol, lchnk )
    call outfld('QT_EVAP_HNO3', qt_evap, ncol, lchnk )
  end if
!
  return
end subroutine neu_wetdep_tend
!
end module mo_neu_wetdep
