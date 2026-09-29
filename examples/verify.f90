! ifort verify.f90 -o verify -g `pkg-config imas-ifort --cflags` `pkg-config imas-ifort imas-lowlevel --libs`

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
    type (ids_core_profiles)   :: cp_ids, cp2
    type(ids_edge_sources) :: edge_ids
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
    write(*,*) 'reading eq ids time SLICE'
    if(size(eq_in%time_slice)>0)then
       if(size(eq_in%time_slice(1)%profiles_2d)>0)then
         write(*,*) 'grid type',eq_in%time_slice(1)%profiles_2d(1)%grid_type%index
         if(size(eq_in%time_slice(1)%profiles_2d(1)%grid%dim1)>0)then
            write(*,*) 'rgrid',1.0e+02_fp*eq_in%time_slice(1)%profiles_2d(1)%grid%dim1(1)
            write(*,*) 'ygrid', 1.0e+02_fp*eq_in%time_slice(1)%profiles_2d(1)%grid%dim2(1)
         else
           write(*,*)'empty eq_in%time_slice(1)%profiles_2d(1)%grid%dim1'
        endif
      endif
      if(size(eq_in%time)>0) write(*,*) 'time', eq_in%time(1)
      if(size(eq_in%time_slice(1)%profiles_1d%rho_tor_norm)>0)then
         write(*,*) 'rho_tor_norm', eq_in%time_slice(1)%profiles_1d%rho_tor_norm(1)
      else
         write(*,*) 'eq_ids has no time_slice'
      endif
    endif


    call ids_get(idx,"core_profiles", cp_ids, ier)
    if(size(cp_ids%profiles_1d)>0)then
      write(*,*)'cp_ids ions', size(cp_ids%profiles_1d(1)%ion)
      if(size(cp_ids%profiles_1d(1)%ion)>0) then
         do i=1,size(cp_ids%profiles_1d(1)%ion)
           write(*,*)'ion', i, cp_ids%profiles_1d(1)%ion(i)%label, &
               size(cp_ids%profiles_1d(1)%ion), cp_ids%profiles_1d(1)%ion(i)%z_ion
           write(*,*)'el',  cp_ids%profiles_1d(1)%ion(i)%element(1)%atoms_n, &
                            cp_ids%profiles_1d(1)%ion(i)%element(1)%a, &
                            cp_ids%profiles_1d(1)%ion(i)%element(1)%z_n
         enddo
      endif
      if(size(cp_ids%profiles_1d(1)%grid%rho_tor_norm)>0) then
        write(*,*) 'cp_ids rho_tor_norm', size(cp_ids%profiles_1d(1)%grid%rho_tor_norm)
      endif
      ! neutrals
      if(size(cp_ids%profiles_1d(1)%neutral)>0) then
         write(*,*)'cp_ids neutrals', size(cp_ids%profiles_1d(1)%neutral)
         do i=1,size(cp_ids%profiles_1d(1)%neutral)
            write(*,*) cp_ids%profiles_1d(1)%neutral(i)%label
            write(*,*)'temp,dens', &
                     size(cp_ids%profiles_1d(1)%neutral(i)%temperature), &
                     size(cp_ids%profiles_1d(1)%neutral(i)%density)
         enddo
      endif
      if(size(cp_ids%profiles_1d(1)%grid%rho_tor_norm)>0) then
        write(*,*) 'cp_ids rho_tor_norm', size(cp_ids%profiles_1d(1)%grid%rho_tor_norm)
      endif
    endif
    noccur = get_max_occurrences_core_profiles(cp_ids)
    write(*,*)'noccur', noccur
    call ids_get(idx,"core_profiles/1", cp2, ier)
    if(size(cp2%profiles_1d)>0)then
      if(size(cp2%profiles_1d(1)%ion)>0) then
        write(*,*)'cp2 ions', cp2%profiles_1d(1)%ion(1)%label, &
               size(cp2%profiles_1d(1)%ion), cp2%profiles_1d(1)%ion(1)%z_ion
      endif
    endif

    call ids_deallocate(cp_ids)
    call ids_deallocate(cp2)
    call ids_deallocate(eq_in)


    call ids_get(idx,"edge_sources", edge_ids)
    if(size(edge_ids%source)>0) then
      write(*,*)"edge_sources"
      if(size(edge_ids%source(1)%ggd)>0 .and. &
           size(edge_ids%source(1)%ggd(1)%neutral)>0) then
        if(size(edge_ids%source(1)%ggd(1)%neutral(1)%particles)>0) then
          if (size(edge_ids%source(1)%ggd(1)%neutral(1)%particles(1)%values)>0) then
            write(*,*)'edge', size(edge_ids%source(1)%ggd(1)%neutral(1)%particles(1)%values), &
                           edge_ids%source(1)%ggd(1)%neutral(1)%particles(1)%values(1:2)
          else
            write(*,*)'edge_ids neutral particles empty'
          endif
        endif
      else
         write(*,*)'edge_ids%source(1)%ggd size0'
      endif
    endif

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
