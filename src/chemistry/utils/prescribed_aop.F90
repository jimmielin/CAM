!-------------------------------------------------------------------
! Manages reading of prescribed aerosol optical properties (AOPs)
! and applies them to the model's aerosols.
!
! The data are layer optical properties, one field per radiation band,
! held in the *weighted* forms that radiation accumulates over aerosol species:
!
!   lw_tau      absorption optical depth                       -> odap_aer
!   sw_tau      extinction optical depth                       -> tau
!   sw_tau_w    extinction optical depth * ssa                 -> tau_w
!   sw_tau_w_g  extinction optical depth * ssa * asymmetry     -> tau_w_g
!   sw_tau_w_f  extinction optical depth * ssa * forward frac  -> tau_w_f
!
! These weighted products superpose linearly across aerosol populations.
! Prescribed properties may therefore be added to those of the model's aerosols
! (prescribed_aop_style = 'ADD') exactly as the model sums over its own species,
! or = 'REPLACE' overwrites them instead from file, discarding the model AOPs.
!
! After the prescribed properties are applied, the totals are floored
! at their physical bounds (non negative optical depths, moments no larger than the
! scattering optical depth).
! A dataset holding the difference between two runs may be legitimately negative in places,
! e.g., zonal means the negative part routinely exceeds the local background on clean longitudes.
!
! Each namelist prefix names a set of per-band fields in the data file; the band number is
! appended by this module.  For example
!
!   prescribed_aop_sw_tau = 'TAU_AER_SW'
!
! reads TAU_AER_SW_1 ... TAU_AER_SW_<nswbands>.
! A blank prefix means that property is not prescribed,
! so a longwave-only dataset is expressed by setting prescribed_aop_lw_tau alone.
!
! sw_tau_w_f is optional.
! When it is not given it is reconstructed as
!
!   tau_w_f = tau_w_g**2 / tau_w
!
! which follows from f = g*g, the relation the model's own optics use.
! This is exact for a single aerosol population sampled instantaneously,
! but is an approximation for the sum over a mixture of populations, nor for time averaged data.
!
! Writing w_i for the scattering weight tau_i*ssa_i of population i,
! the dataset holds tau_w = sum(w_i), tau_w_g = sum(w_i*g_i) and radiation wants sum(w_i*g_i**2),
! whereas the reconstruction supplies sum(w_i*g_i)**2/sum(w_i).
! The difference is tau_w times the scattering weighted variance of g across the populations,
! so the reconstruction is a lower bound on tau_w_f and errs towards too little forward scattering.
! It is close to exact when one population dominates, and least accurate for a well mixed background aerosol.
!
! The data must be on the model's vertical grid.
!-------------------------------------------------------------------
module prescribed_aop

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use cam_abortutils, only : endrun
  use spmd_utils,     only : masterproc
  use tracer_data,    only : trfld, trfile
  use cam_logfile,    only : iulog
  use ppgrid,         only : pcols, pver
  use physics_types,  only : physics_state
  use radconstants,   only : nswbands, nlwbands

  implicit none
  private
  save

  public :: prescribed_aop_readnl
  public :: prescribed_aop_init
  public :: prescribed_aop_adv
  public :: prescribe_aop

  ! Applies the prescribed properties to the shortwave or longwave properties of the model's
  ! aerosols, according to prescribed_aop_style.
  interface prescribe_aop
     module procedure prescribe_sw_aop
     module procedure prescribe_lw_aop
  end interface prescribe_aop

  logical :: has_prescribed_aop = .false.

  type(trfld), pointer :: fields(:) => null()
  type(trfile)         :: file

  ! Namelist options.
  character(len=256) :: filename  = 'NONE'
  character(len=256) :: filelist  = ' '
  character(len=256) :: datapath  = ' '
  character(len=32)  :: datatype  = 'CYCLICAL'
  character(len=32)  :: style     = 'NONE'
  logical            :: rmv_file  = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  ! Field name prefixes; the band number is appended to each.
  character(len=32) :: lw_tau_prefix     = ' '
  character(len=32) :: sw_tau_prefix     = ' '
  character(len=32) :: sw_tau_w_prefix   = ' '
  character(len=32) :: sw_tau_w_g_prefix = ' '
  character(len=32) :: sw_tau_w_f_prefix = ' '

  logical :: prescribe_lw = .false.
  logical :: prescribe_sw = .false.

  ! Set when the dataset carries no sw_tau_w_f, in which case it is reconstructed from tau_w_g
  ! and tau_w. n.b. the reconstruction is a lower bound.
  logical :: reconstruct_tau_w_f = .false.

  ! Fields are read in the order in which the specifier list is built, so the field holding a
  ! given band of a given property is known without searching by name.
  ! These are the offsets that precede band 1 of each property.
  integer :: off_lw_tau     = 0
  integer :: off_sw_tau     = 0
  integer :: off_sw_tau_w   = 0
  integer :: off_sw_tau_w_g = 0
  integer :: off_sw_tau_w_f = 0

  character(len=32), allocatable :: specifier(:)
  integer :: number_flds = 0

  ! Suffix given to the history fields which report the data as read.
  character(len=*), parameter :: hist_suffix = '_D'

contains

!-------------------------------------------------------------------
! reads namelist options and builds the list of fields to read
!-------------------------------------------------------------------
subroutine prescribed_aop_readnl(nlfile)

   use namelist_utils, only : find_group_name
   use units,          only : getunit, freeunit
   use spmd_utils,     only : mpicom, masterprocid, mpi_character, mpi_integer, mpi_logical

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, nflds, n
   character(len=*), parameter :: subname = 'prescribed_aop_readnl'

   character(len=256) :: prescribed_aop_file
   character(len=256) :: prescribed_aop_filelist
   character(len=256) :: prescribed_aop_datapath
   character(len=32)  :: prescribed_aop_type
   character(len=32)  :: prescribed_aop_style
   logical            :: prescribed_aop_rmfile
   integer            :: prescribed_aop_cycle_yr
   integer            :: prescribed_aop_fixed_ymd
   integer            :: prescribed_aop_fixed_tod
   character(len=32)  :: prescribed_aop_lw_tau
   character(len=32)  :: prescribed_aop_sw_tau
   character(len=32)  :: prescribed_aop_sw_tau_w
   character(len=32)  :: prescribed_aop_sw_tau_w_g
   character(len=32)  :: prescribed_aop_sw_tau_w_f

   namelist /prescribed_aop_nl/ &
      prescribed_aop_file,       &
      prescribed_aop_filelist,   &
      prescribed_aop_datapath,   &
      prescribed_aop_type,       &
      prescribed_aop_style,      &
      prescribed_aop_rmfile,     &
      prescribed_aop_cycle_yr,   &
      prescribed_aop_fixed_ymd,  &
      prescribed_aop_fixed_tod,  &
      prescribed_aop_lw_tau,     &
      prescribed_aop_sw_tau,     &
      prescribed_aop_sw_tau_w,   &
      prescribed_aop_sw_tau_w_g, &
      prescribed_aop_sw_tau_w_f

   ! Initialize namelist variables from local module variables.
   prescribed_aop_file       = filename
   prescribed_aop_filelist   = filelist
   prescribed_aop_datapath   = datapath
   prescribed_aop_type       = datatype
   prescribed_aop_style      = style
   prescribed_aop_rmfile     = rmv_file
   prescribed_aop_cycle_yr   = cycle_yr
   prescribed_aop_fixed_ymd  = fixed_ymd
   prescribed_aop_fixed_tod  = fixed_tod
   prescribed_aop_lw_tau     = lw_tau_prefix
   prescribed_aop_sw_tau     = sw_tau_prefix
   prescribed_aop_sw_tau_w   = sw_tau_w_prefix
   prescribed_aop_sw_tau_w_g = sw_tau_w_g_prefix
   prescribed_aop_sw_tau_w_f = sw_tau_w_f_prefix

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_aop_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_aop_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(prescribed_aop_file,       len(prescribed_aop_file),       mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_filelist,   len(prescribed_aop_filelist),   mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_datapath,   len(prescribed_aop_datapath),   mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_type,       len(prescribed_aop_type),       mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_style,      len(prescribed_aop_style),      mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_rmfile,     1,                              mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_cycle_yr,   1,                              mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_fixed_ymd,  1,                              mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_fixed_tod,  1,                              mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_lw_tau,     len(prescribed_aop_lw_tau),     mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_sw_tau,     len(prescribed_aop_sw_tau),     mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_sw_tau_w,   len(prescribed_aop_sw_tau_w),   mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_sw_tau_w_g, len(prescribed_aop_sw_tau_w_g), mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(prescribed_aop_sw_tau_w_f, len(prescribed_aop_sw_tau_w_f), mpi_character, masterprocid, mpicom, ierr)

   ! Update module variables with user settings.
   filename          = prescribed_aop_file
   filelist          = prescribed_aop_filelist
   datapath          = prescribed_aop_datapath
   datatype          = prescribed_aop_type
   style             = prescribed_aop_style
   rmv_file          = prescribed_aop_rmfile
   cycle_yr          = prescribed_aop_cycle_yr
   fixed_ymd         = prescribed_aop_fixed_ymd
   fixed_tod         = prescribed_aop_fixed_tod
   lw_tau_prefix     = prescribed_aop_lw_tau
   sw_tau_prefix     = prescribed_aop_sw_tau
   sw_tau_w_prefix   = prescribed_aop_sw_tau_w
   sw_tau_w_g_prefix = prescribed_aop_sw_tau_w_g
   sw_tau_w_f_prefix = prescribed_aop_sw_tau_w_f

   select case (trim(style))
   case ('NONE', 'ADD', 'REPLACE')
      ! recognized
   case default
      call endrun(subname//': unknown prescribed_aop_style: '//trim(style))
   end select

   has_prescribed_aop = trim(style) /= 'NONE'
   if (.not. has_prescribed_aop) return

   if (len_trim(filename) == 0 .or. trim(filename) == 'NONE') then
      call endrun(subname//': prescribed_aop_style is '//trim(style)//' but prescribed_aop_file is not set')
   end if

   prescribe_lw = len_trim(lw_tau_prefix) > 0
   prescribe_sw = len_trim(sw_tau_prefix)     > 0 .or. len_trim(sw_tau_w_prefix)   > 0 .or. &
                  len_trim(sw_tau_w_g_prefix) > 0 .or. len_trim(sw_tau_w_f_prefix) > 0

   ! The extinction and its first two moments are only meaningful together:
   ! radiation cannot describe a scattering layer without all three.
   ! tau_w_f is reconstructed if it is not given.
   if (prescribe_sw) then
      if (len_trim(sw_tau_prefix) == 0 .or. len_trim(sw_tau_w_prefix) == 0 .or. &
          len_trim(sw_tau_w_g_prefix) == 0) then
         call endrun(subname//': prescribing shortwave properties requires prescribed_aop_sw_tau,'// &
                     ' prescribed_aop_sw_tau_w and prescribed_aop_sw_tau_w_g')
      end if
   end if

   reconstruct_tau_w_f = prescribe_sw .and. len_trim(sw_tau_w_f_prefix) == 0

   if (.not. (prescribe_lw .or. prescribe_sw)) then
      call endrun(subname//': prescribed_aop_style is '//trim(style)//' but no field name prefixes were set')
   end if

   nflds = 0
   if (prescribe_lw) nflds = nflds + nlwbands
   if (prescribe_sw) nflds = nflds + 3*nswbands
   if (prescribe_sw .and. .not. reconstruct_tau_w_f) nflds = nflds + nswbands

   allocate(specifier(nflds), stat=ierr)
   if (ierr /= 0) call endrun(subname//': failed to allocate specifier')

   n = 0
   if (prescribe_lw) then
      call add_band_fields(lw_tau_prefix, nlwbands, off_lw_tau, n)
   end if
   if (prescribe_sw) then
      call add_band_fields(sw_tau_prefix,     nswbands, off_sw_tau,     n)
      call add_band_fields(sw_tau_w_prefix,   nswbands, off_sw_tau_w,   n)
      call add_band_fields(sw_tau_w_g_prefix, nswbands, off_sw_tau_w_g, n)
      if (.not. reconstruct_tau_w_f) then
         call add_band_fields(sw_tau_w_f_prefix, nswbands, off_sw_tau_w_f, n)
      end if
   end if

end subroutine prescribed_aop_readnl

!-------------------------------------------------------------------
! appends <prefix>_1 ... <prefix>_<nbands> to the list of fields to read, and records the
! offset preceding band 1 so that the field for a band can be found later without searching
!-------------------------------------------------------------------
subroutine add_band_fields(prefix, nbands, offset, n)

   use cam_history, only : fieldname_len

   character(len=*), intent(in)    :: prefix
   integer,          intent(in)    :: nbands
   integer,          intent(out)   :: offset
   integer,          intent(inout) :: n

   integer :: iband
   character(len=8) :: band
   character(len=*), parameter :: subname = 'prescribed_aop::add_band_fields'

   offset = n

   do iband = 1, nbands
      write(band,'(i0)') iband

      ! The history field carries the suffix as well, and is the longer of the two names, so
      ! checking it here keeps a long prefix from failing later inside addfld.
      if (len_trim(prefix) + len_trim(band) + 1 + len(hist_suffix) > fieldname_len) then
         call endrun(subname//': field name is too long for prefix '//trim(prefix))
      end if

      n = n + 1
      specifier(n) = trim(prefix)//'_'//trim(band)
   end do

end subroutine add_band_fields

!-------------------------------------------------------------------
! opens the dataset and sets up the history fields
!-------------------------------------------------------------------
subroutine prescribed_aop_init()

   use tracer_data, only : trcdata_init
   use cam_history, only : addfld

   integer :: i

   if (.not. has_prescribed_aop) return

   if (masterproc) then
      write(iulog,*) 'prescribed_aop_init: aerosol optical properties are prescribed in '//trim(filename)
      write(iulog,*) 'prescribed_aop_init: prescribed_aop_style = '//trim(style)
      if (reconstruct_tau_w_f) then
         write(iulog,*) 'prescribed_aop_init: no prescribed_aop_sw_tau_w_f given, so tau_w_f is'
         write(iulog,*) 'prescribed_aop_init: reconstructed as tau_w_g**2/tau_w, which is a lower'
         write(iulog,*) 'prescribed_aop_init: bound for a mixture of aerosols or for averaged data'
      end if
   end if

   allocate(file%in_pbuf(size(specifier)))
   file%in_pbuf(:) = .false.

   ! Optical depth is a layer quantity, not a mixing ratio, so the data are used on the levels
   ! they were written on rather than interpolated in pressure.
   file%no_vert_interp = .true.

   call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                      rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

   number_flds = 0
   if (associated(fields)) number_flds = size(fields)

   ! All of these are optical depths or optical depths weighted by dimensionless moments, so the
   ! units are fixed here rather than taken from the dataset, which need not carry a units
   ! attribute.
   do i = 1, number_flds
      call addfld(trim(fields(i)%fldnam)//hist_suffix, (/ 'lev' /), 'A', '1', &
                  'prescribed aerosol optical property')
   end do

end subroutine prescribed_aop_init

!-------------------------------------------------------------------
! advances the prescribed fields to the current time step
!-------------------------------------------------------------------
subroutine prescribed_aop_adv( state, pbuf2d )

   use tracer_data,    only : advance_trcdata
   use ppgrid,         only : begchunk, endchunk
   use cam_history,    only : outfld
   use physics_buffer, only : physics_buffer_desc

   type(physics_state),       intent(in) :: state(begchunk:endchunk)
   type(physics_buffer_desc), pointer    :: pbuf2d(:,:)

   integer :: c, ncol, i
   !-----------------------------------------------------------------------------

   if (.not. has_prescribed_aop) return

   call advance_trcdata( fields, file, state, pbuf2d )

   do i = 1, number_flds
!$OMP PARALLEL DO PRIVATE (C, NCOL)
      do c = begchunk, endchunk
         ncol = state(c)%ncol
         call outfld( trim(fields(i)%fldnam)//hist_suffix, fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk )
      end do
   end do

end subroutine prescribed_aop_adv

!-------------------------------------------------------------------
! applies the prescribed shortwave properties
!-------------------------------------------------------------------
subroutine prescribe_sw_aop(list_idx, state, tau, tau_w, tau_w_g, tau_w_f)

   integer,             intent(in)    :: list_idx  ! index of the climate or a diagnostic list
   type(physics_state), intent(in)    :: state

   ! Shapes must match aer_rad_props_sw exactly.
   ! Layer 0 lies above the model top and carries no aerosol, so only layers 1:pver are touched here.
   real(r8), intent(inout) :: tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8), intent(inout) :: tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8), intent(inout) :: tau_w_g(pcols,0:pver,nswbands) ! aerosol asymmetry parameter * tau * w
   real(r8), intent(inout) :: tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * tau * w

   integer :: iband, ncol, lchnk, i, k
   character(len=*), parameter :: subname = 'prescribe_sw_aop'

   ! Prescribed properties for one band.  Automatic arrays, so that this remains thread safe when
   ! called from the loop over chunks.
   real(r8) :: tau_p    (pcols,pver)
   real(r8) :: tau_w_p  (pcols,pver)
   real(r8) :: tau_w_g_p(pcols,pver)
   real(r8) :: tau_w_f_p(pcols,pver)

   if (.not. has_prescribed_aop) return
   if (.not. prescribe_sw) return

   ! Prescribed properties are applied to the climate calculation only:
   if (list_idx /= 0) return

   ncol  = state%ncol
   lchnk = state%lchnk

   do iband = 1, nswbands

      tau_p    (:ncol,:) = fields(off_sw_tau    +iband)%data(:ncol,:,lchnk)
      tau_w_p  (:ncol,:) = fields(off_sw_tau_w  +iband)%data(:ncol,:,lchnk)
      tau_w_g_p(:ncol,:) = fields(off_sw_tau_w_g+iband)%data(:ncol,:,lchnk)

      if (reconstruct_tau_w_f) then
         ! f = g*g, so tau*w*f = (tau*w*g)**2 / (tau*w).
         ! A layer which scatters nothing contributes no forward scattering.
         ! n.b. for a mixture, or for time averaged data, this is a lower bound on tau_w_f.
         do k = 1, pver
            do i = 1, ncol
               if (tau_w_p(i,k) > 0._r8) then
                  tau_w_f_p(i,k) = tau_w_g_p(i,k)**2 / tau_w_p(i,k)
               else if (tau_w_p(i,k) < 0._r8) then
                  ! A scattering optical depth cannot be negative. In this case the input dataset
                  ! appears to be a difference between two runs, as such the reconstruction is meaningless:
                  call endrun(subname//': tau_w is negative, so the dataset is not absolute optical'// &
                              ' properties and tau_w_f cannot be reconstructed from it.'// &
                              '  Set prescribed_aop_sw_tau_w_f.')
               else
                  tau_w_f_p(i,k) = 0._r8
               end if
            end do
         end do
      else
         tau_w_f_p(:ncol,:) = fields(off_sw_tau_w_f+iband)%data(:ncol,:,lchnk)
      end if

      select case (trim(style))

      case ('ADD')

         tau    (:ncol,1:pver,iband) = tau    (:ncol,1:pver,iband) + tau_p    (:ncol,:)
         tau_w  (:ncol,1:pver,iband) = tau_w  (:ncol,1:pver,iband) + tau_w_p  (:ncol,:)
         tau_w_g(:ncol,1:pver,iband) = tau_w_g(:ncol,1:pver,iband) + tau_w_g_p(:ncol,:)
         tau_w_f(:ncol,1:pver,iband) = tau_w_f(:ncol,1:pver,iband) + tau_w_f_p(:ncol,:)

      case ('REPLACE')

         tau    (:ncol,1:pver,iband) = tau_p    (:ncol,:)
         tau_w  (:ncol,1:pver,iband) = tau_w_p  (:ncol,:)
         tau_w_g(:ncol,1:pver,iband) = tau_w_g_p(:ncol,:)
         tau_w_f(:ncol,1:pver,iband) = tau_w_f_p(:ncol,:)

      end select

   end do

   ! Floor the totals at their physical bounds:
   tau    (:ncol,1:pver,:) = max(tau(:ncol,1:pver,:), 0._r8)
   tau_w  (:ncol,1:pver,:) = min(max(tau_w(:ncol,1:pver,:), 0._r8), tau(:ncol,1:pver,:))
   tau_w_g(:ncol,1:pver,:) = min(max(tau_w_g(:ncol,1:pver,:), -tau_w(:ncol,1:pver,:)), &
                                 tau_w(:ncol,1:pver,:))
   tau_w_f(:ncol,1:pver,:) = min(max(tau_w_f(:ncol,1:pver,:), 0._r8), tau_w(:ncol,1:pver,:))

end subroutine prescribe_sw_aop

!-------------------------------------------------------------------
! applies the prescribed longwave properties
!-------------------------------------------------------------------
subroutine prescribe_lw_aop(list_idx, state, odap_aer)

   integer,             intent(in)    :: list_idx  ! index of the climate or a diagnostic list
   type(physics_state), intent(in)    :: state

   real(r8), intent(inout) :: odap_aer(pcols,pver,nlwbands) ! absorption optical depth, per layer

   integer :: iband, ncol, lchnk
   !-----------------------------------------------------------------------------

   if (.not. has_prescribed_aop) return
   if (.not. prescribe_lw) return

   ! Prescribed properties are applied to the climate calculation only, so that the diagnostic
   ! call lists remain a view of the model's own aerosols.
   if (list_idx /= 0) return

   ncol  = state%ncol
   lchnk = state%lchnk

   select case (trim(style))

   case ('ADD')

      do iband = 1, nlwbands
         odap_aer(:ncol,:,iband) = odap_aer(:ncol,:,iband) + fields(off_lw_tau+iband)%data(:ncol,:,lchnk)
      end do

   case ('REPLACE')

      do iband = 1, nlwbands
         odap_aer(:ncol,:,iband) = fields(off_lw_tau+iband)%data(:ncol,:,lchnk)
      end do

   end select

   ! Floor the totals at their physical bounds:
   odap_aer(:ncol,:,:) = max(odap_aer(:ncol,:,:), 0._r8)

end subroutine prescribe_lw_aop

end module prescribed_aop
