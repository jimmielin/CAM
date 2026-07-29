! Minimal stubs of the CAM modules used by mo_micm.F90, configured for the
! trop_mam4 mechanism, so the real shim can be compiled and exercised
! against the generated MICM configuration outside of CAM.
module shr_kind_mod
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: shr_kind_r8 = real64
  integer, parameter :: shr_kind_cl = 256
end module shr_kind_mod

module spmd_utils
  implicit none
  logical, parameter :: masterproc = .true.
  integer, parameter :: masterprocid = 0
  integer, parameter :: mpicom = 0
  integer, parameter :: mpi_character = 1
  integer, parameter :: mpi_logical = 2
  integer, parameter :: mpi_success = 0
end module spmd_utils

module cam_abortutils
  implicit none
contains
  subroutine endrun(msg)
    character(len=*), intent(in) :: msg
    print '(2a)', 'ENDRUN: ', trim(msg)
    stop 1
  end subroutine endrun
end module cam_abortutils

module cam_logfile
  use iso_fortran_env, only: output_unit
  implicit none
  integer, parameter :: iulog = output_unit
end module cam_logfile

module ppgrid
  implicit none
  ! pcols*pver = 144 > the vector size (128), so a full chunk exercises the
  ! multi-block solve loop in micm_solve, and a partial chunk the padding
  integer, parameter :: pver = 48
  integer, parameter :: pcols = 3
end module ppgrid

module chem_mods
  implicit none
  integer, parameter :: gas_pcnst = 26
  integer, parameter :: rxntot = 10
  integer, parameter :: extcnt = 9
  character(len=16) :: extfrc_lst(extcnt) = (/ 'SO2   ', 'so4_a1', 'so4_a2', &
       'pom_a4', 'bc_a4 ', 'H2O   ', 'num_a1', 'num_a2', 'num_a4' /)
end module chem_mods

module mo_tracname
  use chem_mods, only: gas_pcnst
  implicit none
  character(len=16) :: solsym(gas_pcnst) = (/ &
       'bc_a1 ', 'bc_a4 ', 'DMS   ', 'dst_a1', 'dst_a2', 'dst_a3', &
       'H2O2  ', 'H2SO4 ', 'ncl_a1', 'ncl_a2', 'ncl_a3', 'num_a1', &
       'num_a2', 'num_a3', 'num_a4', 'pom_a1', 'pom_a4', 'SO2   ', &
       'so4_a1', 'so4_a2', 'so4_a3', 'soa_a1', 'soa_a2', 'SOAE  ', &
       'SOAG  ', 'H2O   ' /)
end module mo_tracname

module physconst
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  real(r8), parameter :: avogad = 6.02214e26_r8 ! molecules/kmole
end module physconst

module namelist_utils
  implicit none
contains
  subroutine find_group_name(unit, group, status)
    integer,          intent(in)  :: unit
    character(len=*), intent(in)  :: group
    integer,          intent(out) :: status
    character(len=256) :: line
    do
      read(unit,'(a)',iostat=status) line
      if (status /= 0) return
      line = adjustl(line)
      if (line(1:1) == '&' .and. index(line, trim(group)) == 2) then
        backspace(unit)
        status = 0
        return
      end if
    end do
  end subroutine find_group_name
end module namelist_utils
