!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hemco_interface
!
! !DESCRIPTION: Module HEMCO\_INTERFACE is the HEMCO-CESM interface module.
!               CESM operates on chunks thus the interface is called HCOI\_Chunk.
!               Internally it uses a gridded component to interact with the CAM
!               physics grid; these functions are internal and called HCO\_GC...
!\\
!\\
! !INTERFACE:
!
module hemco_interface
!
! !USES:
!
    ! Controls
    use cam_abortutils,           only: endrun      ! fatal terminator
    use cam_logfile,              only: iulog       ! output log handle

    ! Species information 
    use chem_mods,                only: gas_pcnst   ! # of species
    use chem_mods,                only: adv_mass    ! advected mass
    use mo_tracname,              only: solsym      ! species names
    use mo_chem_utls,             only: get_spc_ndx ! IND_

    ! Grid
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs

    ! ESMF types
    use ESMF,                     only: ESMF_Mesh, ESMF_DistGrid
    use ESMF,                     only: ESMF_State, ESMF_Clock, ESMF_GridComp
    use ESMF,                     only: ESMF_Field, ESMF_RouteHandle
    use ESMF,                     only: ESMF_KIND_R8, ESMF_KIND_I4, ESMF_SUCCESS

    implicit none
    private
!
! !PRIVATE MEMBER FUNCTIONS:
!
    private :: HCO_GC_Init_Wrapper
    private :: HCO_GC_Init
    private :: HCO_GC_SetServices
    private :: HCO_GC_Run
    private :: HCO_GC_Final
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public  :: hemco_readnl
    public  :: HCOI_Chunk_Init
    public  :: HCOI_Chunk_Run
    public  :: HCOI_Chunk_Final
!
! !REVISION HISTORY:
!  29 Jan 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    type(ESMF_Mesh)      :: CAM_PhysMesh        ! Copy of CAM physics mesh
    type(ESMF_DistGrid)  :: CAM_DistGrid        ! DE-local allocation descriptor DistGrid (2D)
    type(ESMF_GridComp)  :: HCO_GridComp        ! HEMCO GridComp
    type(ESMF_State)     :: HCO_GridCompState   ! HEMCO GridComp Import/Export State
    
    integer              :: col_start, col_end  ! idx of columns in this PET
    integer              :: col_total           ! # of columns in this PET
    integer              :: nlev = 0            ! # of levs in this PET
contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hemco_readnl
!
! !DESCRIPTION: Reads the namelist from cam/src/control/runtime_opts.
!\\
!\\
! !INTERFACE:
!
    subroutine hemco_readnl( nlfile )
        use namelist_utils, only: find_group_name
        use units,          only: getunit, freeunit
        use spmd_utils,     only: mpi_real8, mpi_logical, mpi_integer, mpi_character
        use spmd_utils,     only: masterproc, mpicom, masterprocid

        character(len=*), intent(in) :: nlfile ! namelist file
!
! !REVISION HISTORY:
!  31 Jan 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer :: unitn, ierr
        character(len=*), parameter  :: subname = 'hemco_readnl'
        character(len=256)           :: hemco_config_file = 'HEMCO_Config.rc'

        namelist /hemco_nl/ hemco_config_file

        ! Read namelist on master proc
        ! ...
        if(masterproc) then
            write(iulog,*) "This is hemco_readnl - Reading HEMCO Namelist in CESM"
            unitn = getunit()
            open(unitn, file=trim(nlfile), status='old')
            call find_group_name(unitn, 'hemco_nl', status=ierr)
            if(ierr == 0) then
                read(unitn, hemco_nl, iostat=ierr)
                if(ierr /= 0) then
                    call endrun(subname // ':: ERROR reading namelist')
                endif
            endif
            close(unitn)
            call freeunit(unitn)

            write(iulog,*) "hemco_readnl: hemco config file = ", hemco_config_file
        endif

        ! Read HEMCO Configuration file then broadcast
        ! ...

        
    end subroutine hemco_readnl
!EOC

    ! Init
    subroutine HCOI_Chunk_Init()
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid

        ! Local variables
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Init'

        ! Initialize HEMCO!
        if(masterproc) then
            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            write(iulog,*) "HEMCO: Harmonized Emissions Component"
            write(iulog,*) "HEMCO_CAM interface version 0.1"
            write(iulog,*) "You are using HEMCO version 3.0.placeholder"
            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        endif
    end subroutine HCOI_Chunk_Init

    ! Run
    subroutine HCOI_Chunk_Run(cam_in, phys_state, pbuf2d, phase)
        ! Type descriptors
        use camsrfexch,     only: cam_in_t
        use physics_types,  only: physics_state
        use physics_buffer, only: physics_buffer_desc

        ! Output and mpi
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid

        ! Input
        type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
        type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
        type(physics_buffer_desc), pointer :: pbuf2d(:,:)
        integer, intent(in) :: phase ! 1, 2

        ! Local variables
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Run'

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Running HCOI_Chunk_Run"
        endif


        ! Run the grid comp
        ! call ESMF_GridCompRun ...
        
    end subroutine HCOI_Chunk_Run

    ! Final
    subroutine HCOI_Chunk_Final()
        ! Stub...
    end subroutine HCOI_Chunk_Final

    !-----------------------------------------------------------------------
    !              H E M C O   W R A P P E R   G R I D C O M P             !
    !-----------------------------------------------------------------------
    !  Below code includes internal routines used to wrap a ESMF gridded   !
    !  component around the HEMCO interface, so ESMF can handle all the    !
    !  interaction with the physics mesh.                                  !
    !                                                                      !
    !  This is largely based on edyn_grid_comp.F90 from ionosphere/waccmx  !
    !  Thanks to Steve Goldhaber for the example                           !
    !                                                                      !
    !  (hplin, 1/31/20)                                                    !
    !-----------------------------------------------------------------------

    ! Run routine called from the ESMF GridComp. Contains the actual computation
    ! routines (secret sauce) operating on the HEMCO grid and the regridding
    ! routines to return it into the physics mesh.
    subroutine HCO_GC_Run()
        ! Stub...
    end subroutine HCO_GC_Run

    ! Method to initialize the series of ESMF wrappers
    subroutine HCO_GC_Init_Wrapper(mpicomm)
        use mpi,          only: MPI_INTEGER
        use ESMF,         only: ESMF_StateCreate, ESMF_GridCompInitialize
        use ESMF,         only: ESMF_GridCompCreate, ESMF_GridCompSetServices
        use ESMF,         only: ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
        use cam_instance, only: inst_index, inst_name

        integer, intent(in)                   :: mpicomm
        character(len=*),       parameter     :: subname = 'HCO_GC_Init'

        ! Local variables
        integer,                allocatable   :: petlist(:)
        integer                               :: iam
        integer                               :: npes
        integer                               :: localPet
        integer                               :: petCount
        integer                               :: RC
        type(ESMF_VM)                         :: vm_init


        ! Stub
    end subroutine HCO_GC_Init_Wrapper


    !---------------------------------------------------------------------
    ! Every HEMCO GridComp routine below likely will be just boilerplate
    ! that does not need extensive maintenance
    !---------------------------------------------------------------------

    ! IRF routines for the Gridded Component
    subroutine HCO_GC_Init(GC, IMPORT, EXPORT, Clock, RC)
        use ESMF,         only: ESMF_DistGridCreate, ESMF_MeshCreate
        use ESMF,         only: ESMF_FILEFORMAT_ESMFMESH
        use cam_instance, only: inst_name
        use phys_control, only: phys_getopts
        use phys_grid,    only: get_ncols_p, get_gcol_p

        ! Dummy arguments
        type(ESMF_GridComp)  :: GC
        type(ESMF_State)     :: IMPORT
        type(ESMF_State)     :: EXPORT
        type(ESMF_Clock)     :: Clock
        integer, intent(out) :: RC

        ! Local variables
        integer                               :: ncols
        integer                               :: chnk, col, dindex
        integer,                allocatable   :: decomp(:)
        character(len=256)                    :: grid_file
        character(len=*),       parameter     :: subname = 'HCO_GC_Init'

        ! Stub
        RC = ESMF_SUCCESS
    end subroutine HCO_GC_Init

    ! Finalize Gridded Component
    subroutine HCO_GC_Final()
        ! Stub...
    end subroutine HCO_GC_Final

    subroutine HCO_GC_SetServices(GC, RC)
        use ESMF,         only: ESMF_GridCompSetEntryPoint
        use ESMF,         only: ESMF_METHOD_INITIALIZE, ESMF_METHOD_RUN
        use ESMF,         only: ESMF_METHOD_FINALIZE

        type(ESMF_GridComp)                   :: GC
        integer, intent(out)                  :: RC
        character(len=*),       parameter     :: subname = 'HCO_GC_SetServices'

        ! Stub
        RC = ESMF_SUCCESS
    end subroutine HCO_GC_SetServices


!EOC
end module hemco_interface