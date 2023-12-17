

! ifort verify.f90 -o verify -g -r8 -fpp -assume no2underscore  -fPIC
!     `pkg-config imas-ifort --cflags` `pkg-config imas-ifort imas-lowlevel --libs`

subroutine verify_ids(mach, shot, run, sbackend)
  use iso_c_binding, only: fp => c_double, c_int
    use ids_schemas    ! Definitions for the Physics Data Model
    use ids_routines      ! Access Layer routines and management of IDS structures

    implicit none
    character(len=30) :: mach
    character(len=20), intent(in) :: sbackend
    integer :: i, ig, i1, i2, ilx, ix, j1, j2, ntot, nth
    integer :: shot, run, indx, idx
    type (ids_equilibrium)   :: eq_in ! equilibrium  
    type (ids_core_profiles)   :: cp_in, cp2
    integer :: nt, nx, nr, nz, ir, iz, it, noccur
    integer :: intime, iltime, ildata, ier, isized, istore, isizb, iaddrs
    integer :: iaddrx, iaddrd, iadr, i0
    !real(fp), parameter :: twopiinv = 1.0_fp/twopi
    real(fp) :: zt1, zt2, zt, zfac, zmin, zmax, z_ion, a_ion 
    real(fp), dimension(:,:), allocatable :: zw2d, zw2df
    real(fp), dimension(:,:,:), allocatable :: zw3d
    real(fp), dimension(:), allocatable :: time, zw1d, xgrid
    real(fp), dimension(:), allocatable :: rgrid, ygrid, zgrid
    real(fp), dimension(:,:), allocatable :: psirz, psirz_tg1, psirz_tg2
    integer, dimension(:), allocatable :: fgrid
    real(fp) :: ztime, zdata, dintrp
    character(len=132):: usr
    character(len=5)::treename
    integer :: ibackend, status

    treename = "ids"
    if (sbackend .eq. 'HDF5') then
       ibackend = HDF5_BACKEND
    else
       ibackend = MDSPLUS_BACKEND
    endif   

    call get_environment_variable("USER",usr)
    call ual_begin_pulse_action(ibackend,shot,run,trim(usr),trim(mach),"3",idx,status)
    call ual_open_pulse(idx, OPEN_PULSE,'', status)
    !call imas_open_env(treename,shot,run,idx,trim(usr), trim(mach),'3')

    call ids_get(idx,"equilibrium", eq_in)
    write(*,*) 'reading eq ids'
    write(*,*) 'grid type',eq_in%time_slice(1)%profiles_2d(1)%grid_type%index
    nt = size(eq_in%time)
    nr = size(eq_in%time_slice(1)%profiles_2d(1)%grid%dim1)
    nz = size(eq_in%time_slice(1)%profiles_2d(1)%grid%dim2)
    nx = size(eq_in%time_slice(1)%profiles_1d%rho_tor_norm)

    if (allocated(rgrid)) deallocate(rgrid)
    if (allocated(ygrid)) deallocate(ygrid)
    if (allocated(psirz)) deallocate(psirz)
    if (allocated(psirz_tg1)) deallocate(psirz_tg1)
    if (allocated(psirz_tg2)) deallocate(psirz_tg2)
    allocate(rgrid(nr))
    allocate(ygrid(nz))
    allocate(psirz(nr,nz))
    allocate(psirz_tg1(nr,nz))
    allocate(psirz_tg2(nr,nz))
    if (.not.allocated(time))  allocate(time(nt))
    if (.not.allocated(xgrid))  allocate(xgrid(nx))
    time = eq_in%time
    xgrid = eq_in%time_slice(1)%profiles_1d%rho_tor_norm
    rgrid = 1.0e+02_fp*eq_in%time_slice(1)%profiles_2d(1)%grid%dim1
    ygrid = 1.0e+02_fp*eq_in%time_slice(1)%profiles_2d(1)%grid%dim2
    write(*,*) 'rho_tor_norm', xgrid(1)
    write(*,*) 'rgrid', rgrid(1)

    call ids_get(idx,"core_profiles", cp_in, ier)
    write(*,*)'cp_in ions', cp_in%profiles_1d(1)%ion(1)%label, &
               size(cp_in%profiles_1d(1)%ion), cp_in%profiles_1d(1)%ion(1)%z_ion
    write(*,*) 'cp_in rho_tor_norm', cp_in%profiles_1d(1)%grid%rho_tor_norm

    noccur = get_max_occurrences_core_profiles(cp_in)
    write(*,*)'noccur', noccur

    call ids_get(idx,"core_profiles/1", cp2, ier)
    write(*,*)'cp2 ions', cp2%profiles_1d(1)%ion(1)%label, &
               size(cp2%profiles_1d(1)%ion), cp2%profiles_1d(1)%ion(1)%z_ion

    call ids_deallocate(cp_in)
    call ids_deallocate(cp2)
    call ids_deallocate(eq_in)
    call imas_close(idx)
end subroutine



program verify
    implicit none
    integer :: nargs, i, shot, run
    character(len=30) :: mach
    character(len=20) :: sback

    character(len=30), dimension(:), allocatable :: args
    external verify_ids

    nargs = command_argument_count()
    if (nargs<3) then
      write(*,*) './verify machine(eq:SMRT) shot(eq:92210) run(eq:622) backend(MDSPLUS)/HDF5'
      stop
    endif
    allocate(args(nargs))
    do i =1,nargs
       call get_command_argument(i,args(i))
       if (i==1) mach = args(1)
       if (i==2) read(args(i),*) shot
       if (i==3) read(args(i),*) run
       if (i==4) sback=args(i)
    end do

    call verify_ids(mach, shot, run, sback)

end program
