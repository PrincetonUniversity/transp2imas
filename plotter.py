#!/usr/bin/python
"""Plotting helpers used by diagnostic and example workflows."""


import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

colorlist = ['white','gray','black', 'indigo', 'blue','cyan','green','yellow',
     'orange','pink','olive','magenta','brown','purple','violet','red']
cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colorlist,256)
cmap2 =  mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colorlist*10,256)
cmap = mpl.cm.viridis

def format_cursor_data(data):
   return "[" + str(data) + "]"

def check_valid(ids, pre, name):
    pre = pre if pre else ''
    times = 0
    try:
        times = ids.time
        print('homogeneous time', ids.ids_properties.homogeneous_time)
    except (AttributeError, TypeError):
        print(f'{name} not available')
        return 0, 0, pre
    else:
        ntimes = len(times)
    return times, ntimes, pre

def save_pdf(fig, name):
    pdf = PdfPages(name+".pdf")
    pdf.savefig(fig)
    pdf.close()

def finalize_plots(name):
    pdf = PdfPages(name+".pdf")
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pdf, format='pdf')
    pdf.close()


def finalize_plot(plt, fig, fname, save):
    return
    if not save:
        plt.show() #block=False)
        #plt.pause(1)
        #plt.close()
    if save:
        save_pdf(fig, fname)
    plt.close()

def plot_limiters(limitr, limitz, save:bool=True, fname:str=None):
    fname = 'limiters' if not fname else fname
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_title(fname)
    plt.plot(limitr, limitz)
    #finalize_plot(plt, fig, fname, save)

def plot_subgrid(fig, gsp, var1, var2, lab1, lab2):
    if var1.shape != var2.shape:
        print('Error', var1.shape, '!=' , var2.shape)
        return fig
    ax = fig.add_subplot(gsp)
    if var2 is not None:
        plt.plot(var1, var2)
        ax.grid()
    else:
        plt.plot(var1)
    ax.set_xlabel(lab1)
    ax.set_title(lab2)
    return fig

def find_plot_index(ndat, nplots):
    nplots = ndat if ndat < nplots else nplots
    freq = int(ndat/nplots)
    plot_inds = [i for i in range(ndat) if i%freq==0]
    return plot_inds[0:nplots]

def find_nrows_cols(ni:int):
    # number of rows and columns of subgrids per page
    nr = 2
    nc = 2
    if ni > 0:
        nr = int(np.sqrt(ni)+0.5)
        nz = nr
        if ni-nr*nr>0:
            nr = nr+int((ni-nr*nr)/2)
    return nr, nc

#args are (data, name) matching tuples
def plot_data2d_samex(name, nr, nc, var, lab, *args):
    save = True
    for i in range(len(args)):
        data = args[i][0]
        lab1 = args[i][1]
        fig = plt.figure(tight_layout=True)
        gsp = gridspec.GridSpec(nr,nc)
        iplots = find_plot_index(len(data), nr*nc)
        name2 = name + '_' + lab1
        for j, t in enumerate(iplots):
            lab2 = lab1 + str(t)
            fig = plot_subgrid(fig, gsp[j], var[t], data[t], lab, lab2)
        #finalize_plot(plt, fig, name2, save)

def plot_data1d_samex(name, var, lab, *args):
    save = True
    ndata = len(args)
    nr, nc = find_nrows_cols(ndata)
    fig = plt.figure(tight_layout=True)
    gsp = gridspec.GridSpec(nr,nc)
    for i in range(ndata):
        data = args[i][0]
        lab1 = args[i][1]
        fig = plot_subgrid(fig, gsp[i], var, data, lab, lab1)
    #finalize_plot(plt, fig, name, save)

def plot_Ip_vacBR_loopV(times, currents, vac_b0, vac_r0, vloop=None,
     save=True, pre:str=None):
    fig = plt.figure(tight_layout=True)
    gsp = gridspec.GridSpec(2, 2)
    print('plot_Ip_vacBR_loopV:times', len(times), 'ip', len(currents))
    if np.isscalar(times) or np.isscalar(currents):
        return
    if len(times) != len(currents):
        print('times', times)
        print('currents', currents)
    fig = plot_subgrid(fig, gsp[0], times, currents, 'times', 'currents')
    fig = plot_subgrid(fig, gsp[1], times, vac_b0, 'times', 'vac_b0')
    nr = 1 if np.isscalar(vac_r0) else len(vac_r0)
    if nr == len(times):
        fig = plot_subgrid(fig, gsp[2], times, vac_r0, 'times', 'vac_r0')
    else:
        print('Vac_r0', vac_r0)
    name = 'ip_vacB0_r0'
    name = pre+'_'+name if pre else name

    if vloop is not None: #for array
        fig = plot_subgrid(fig, gsp[3], times, vloop, 'times', 'vloop')
        name = name + '_Vloop'
    #fig.title('core profiles')
    #finalize_plot(plt, fig, name, save)


def plot_eq_global(times, rmaxis, zmaxis, b_field_tor, psi_axis,
      psi_bdry, q95, qaxis, li3, save=True, pre:str=None):
    name = 'RZmaxis_Btor_psi'
    name = pre+'_'+name if pre else name
    plot_data1d_samex(name, times, 'times', (rmaxis, 'rmaxis'),
       (zmaxis, 'zmaxis'), (b_field_tor, 'b_field_tor'), (psi_axis, 'psi_axis'))
    name = 'psib_q95_q0_LI3'
    name = pre+'_'+name if pre else name
    plot_data1d_samex(name, times, 'times', (psi_bdry,'psi_bdry'),
     (q95, 'q95'), (qaxis,'qaxis'), (li3, 'li3'))

def plot_eq_outline(times, rminor, outline_r, outline_z, lab=None, save=True,
       pre:str=None):
    name = 'outline'
    name = pre+'_'+name if pre else name
    fig = plt.figure(tight_layout=True)
    gsp = gridspec.GridSpec(2, 2)
    nplots = 3 #rminor plot is extra

    tplots = find_plot_index(len(times), nplots)
    fig = plot_subgrid(fig, gsp[0], times, rminor, 'times', 'rminor')
    for i, t in enumerate(tplots):
        fig = plot_subgrid(fig, gsp[i+1], outline_r[t], outline_z[t],
         'outline_R', f'outline_Z t={t}')
    #finalize_plot(plt, fig, name, save)

def plot_tslice_grid(nr, nc, times, dat, pre, lab, xtime:bool=False, save=True):
    fig = plt.figure(tight_layout=True)
    gsp = gridspec.GridSpec(nr, nc)
    #times on X-axis
    if xtime:
        inds = find_plot_index(dat.shape[1], nr*nc)
        for i,j in enumerate(inds):
            fig = plot_subgrid(fig, gsp[i], times, dat[:,j], 'times', f'{lab} index={j}')
    else: # index for flux surfaces on X axis
        tplots = find_plot_index(len(times), nr*nc)
        for i,t in enumerate(tplots):
            ax = plt.subplot(gsp[i])
            ax.plot(dat[t])
            ax.grid()
            ax.set_title(f'{lab} t={t}')
    #finalize_plot(plt, fig, pre+lab, save)

def plot_equilibrium_1d(times, psi, qpsi, pres, torflux_phi, rho_tor_norm,
        elong, triang_top, triang_bot, jtor, fdp_dspi, dpressure_dpsi,
        pre:str=None, ni:int=0):
    pre = pre+'_' if pre else 'eq_'
    #xtime=False. X-axis is data for selected time slice,
    #xtime=True. X-axis is all times for selected data index.
    xtime = False
    # number of rows and columns of subgrids per page
    nr, nc = find_nrows_cols(ni)
    #selected times of time slices. data must be 2D(time,N)
    plt = plot_tslice_grid(nr, nc, times, psi, pre, 'psi', xtime)
    plt = plot_tslice_grid(nr, nc, times, qpsi, pre, 'qpsi', xtime)
    plt = plot_tslice_grid(nr, nc, times, pres,  pre, 'pres', xtime)
    plt = plot_tslice_grid(nr, nc, times, torflux_phi, pre, 'torflux_phi', xtime)
    plt = plot_tslice_grid(nr, nc, times, rho_tor_norm, pre, 'rho_tor_norm',xtime )
    plt = plot_tslice_grid(nr, nc, times, elong, pre, 'elong', xtime)
    plt = plot_tslice_grid(nr, nc, times, triang_top, pre, 'triang_top', xtime)
    plt = plot_tslice_grid(nr, nc, times, triang_bot, pre, 'triang_bot', xtime)
    plt = plot_tslice_grid(nr, nc, times, jtor, pre, 'jtor', xtime)
    plt = plot_tslice_grid(nr, nc, times, fdp_dspi, pre, 'fdp_dspi', xtime)
    plt = plot_tslice_grid(nr, nc, times, dpressure_dpsi, pre, 'dpressure_dpsi', xtime)

def plot_equil_psi1d_interactive(eq, times):
   if np.isscalar(times):
      times = [times]
   for time in times:
      dat = eq.time_slice[time].profiles_1d.psi
      fig, ax = plt.subplots()
      fig.set_figheight(10)
      fig.set_figwidth(10)
      ax.plot(dat)
      ax.grid()
      ax.set_title(f'psi1d t={time}')
      plt.show()


def plot_equil_2d_interact(rr, zz, dat,time,  ncontour=30, img=False):

   ext = [rr[time], rr[-1], zz[time], zz[-1]]
   levels = [i for i in np.arange(np.min(dat), np.max(dat), (np.max(dat)-np.min(dat))/ncontour)]
   fig, ax = plt.subplots()
   fig.set_figheight(10)
   fig.set_figwidth(10)
   cont = ax.contour(dat, levels=levels, extent=ext, origin='lower', colors='grey')
   if img is True:
       image = plt.imshow(dat, extent=ext, interpolation='nearest', cmap = cmap, origin='lower')
       cbar = fig.colorbar(image, ax=ax)
       cbar.add_lines(cont)
       image.format_cursor_data = format_cursor_data

   xmajor_ticks = np.arange(rr[0], rr[-1], 5*abs(rr[0]-rr[1]))
   xminor_ticks = np.arange(rr[0], rr[-1], abs(rr[0]-rr[1]))
   ymajor_ticks = np.arange(zz[0], zz[-1], 5*abs(zz[0]-zz[1]))
   yminor_ticks = np.arange(zz[0], zz[-1], abs(zz[0]-zz[1]))

   ax.set_xticks(xmajor_ticks)
   ax.set_xticks(xminor_ticks, minor=True)
   ax.set_yticks(ymajor_ticks)
   ax.set_yticks(yminor_ticks, minor=True)
   ##ax.grid(which='both')
   ax.grid(which='minor', alpha=0.4)
   ax.grid(which='major', alpha=0.6)


def plot_equil2d_interactive(eq, name, times, ncontour=30, img=False):
   if np.isscalar(times):
      times = [times]

   nm = name.lower()
   for time in times:
      rr = eq.time_slice[time].profiles_2d[0].grid.dim1
      zz = eq.time_slice[time].profiles_2d[0].grid.dim2
      if nm=='psirz':
         dat = eq.time_slice[time].profiles_2d[0].psi
         rinv = eq.time_slice[time].coordinate_system.r
         zinv = eq.time_slice[time].coordinate_system.z
         plot_equil_2d_interact(rr, zz, dat, time, ncontour, img)

      elif nm=='bfield' or nm=='b_field':
         for i in range(3):
            if i==0:
               dat = eq.time_slice[time].profiles_2d[0].b_field_r
            elif i==1:
               dat = eq.time_slice[time].profiles_2d[0].b_field_tor
            elif i==2:
               dat = eq.time_slice[time].profiles_2d[0].b_field_z
            plot_equil_2d_interact(rr, zz, dat, time, ncontour, img)
      plt.show()

def plot_equilibrium_2d(times, rgrid, zgrid, rinv, zinv, psirz, rz_range:bool=True,
      save=True, pre:str=None, ni:int=0, image2d:bool=False, ncontour:int=10):

    if len(rgrid.shape) ==1:
        print(f'{rgrid} is 1D')
    else:
        rgrid = rgrid[0] #2D
        zgrid = zgrid[0] #2D

    image = image2d #plot image
    print('plotting eq 2d')
    #TODO use ni
    nr, nc = find_nrows_cols(ni)
    print('psirz', psirz.shape, rgrid.shape, zgrid.shape)
    if len(psirz.shape) == 2:
        print('psirz is single array of 2D')
        tplots = [0]
    else:
        tplots = find_plot_index(len(times), nr*nc)

    def format_cursor_data(data):
        return "[" + str(data) + "]"
    ncont = ncontour #number of contour lines
    extent = [rgrid[0], rgrid[-1], zgrid[0], zgrid[-1]]

    fig, ax = plt.subplots()
    gsp = gridspec.GridSpec(nr, nc)

    for i, t in enumerate(tplots):
        ax = plt.subplot(gsp[i])
        if len(tplots)==1:
            dat = psirz
        else:
            dat = psirz[t]
        ax.set_title(f'psirz time={t}')
        levels = [ii for ii in np.arange(np.min(dat), np.max(dat), (np.max(dat)-np.min(dat))/ncont)]
        plt.xlabel('R')
        if rz_range:
            cont = ax.contour(dat, levels=levels, extent=extent) #, origin='lower')
            if image:
                img = plt.imshow(dat, extent=extent) #, cmap=plt.cm.viridis)
        else:
            cont = ax.contour(dat, levels=levels, origin='lower')
            if image:
                img = plt.imshow(dat, origin='lower')
        if image:
            cbar = fig.colorbar(img, ax=ax, fraction=0.05, pad=0.05)
            cbar.add_lines(cont)
            img.format_cursor_data = format_cursor_data

        print('rz_inv', rinv.shape, zinv.shape)
        if len(rinv)>1 and len(zinv)>1:
            rg, zg = rinv[t] , zinv[t]
            plt.scatter(rg[30], zg[30])
        ax.grid()

    #name = f'psirz_{1+int(rz_range==True)}'
    #name = pre+'_'+name if pre else ''
    #finalize_plot(plt, fig, name, save)


def declare_np_arrays(n, m):
    vars = []
    for i in range(n):
        vars.append(np.zeros(m, dtype=float))
    return vars

#ni number of indices, 0 will use default
def plot_equilibrium_ids(eq, ni:int=0, pre:str=None, name='eq', save=False):
    nslices = len(eq.time_slice)
    times = eq.time
    print(f'times {times.shape} n_times {nslices}')
    vac_r0 = eq.vacuum_toroidal_field.r0 #single value
    vac_b0 = eq.vacuum_toroidal_field.b0
    currents, rmaxis, zmaxis, b_field_tor, psi_axis,\
      psi_bdry, q95, qaxis, li3, rminor = declare_np_arrays(10, nslices)

    noutline = eq.time_slice[0].boundary.outline.r.shape[0]

    outline_r, outline_z = declare_np_arrays(2, (nslices, noutline))
    npsi = eq.time_slice[0].profiles_1d.psi.shape[0]
    psi, qpsi, pres, torflux_phi, rho_tor_norm, elong, triang_top, \
      triang_bot, jtor, fdp_dspi, dpressure_dpsi = declare_np_arrays(11, (nslices, npsi))

    prof2d = False
    nr, nz = 1, 1
    try:
        nr = eq.time_slice[0].profiles_2d[0].grid.dim1.shape[0]
        nz = eq.time_slice[0].profiles_2d[0].grid.dim2.shape[0]
        print(nr,nz)
    except (AttributeError, IndexError, TypeError):
        print('2d profile not found')
    else:
        prof2d = True
    rgrid = np.zeros((nslices, nr), dtype=float)
    zgrid = np.zeros((nslices, nz), dtype=float)
    psirz = np.zeros((nslices, nr, nz), dtype=float)

    #inverse grid
    nsi, nth = 1, 1
    try:
        sh = eq.time_slice[0].coordinate_system.r.shape
        nsi, nth = sh[0], sh[1]
    except (AttributeError, IndexError, TypeError):
        print('2d inverse grid not found')

    rinv = np.zeros((nslices, nsi, nth), dtype=float)
    zinv = np.zeros((nslices, nsi, nth), dtype=float)
    print('psirz ', psirz.shape, 'rzinv_grid', rinv.shape)

    for i in range(nslices):
        currents[i] = eq.time_slice[i].global_quantities.ip
        rmaxis[i] = eq.time_slice[i].global_quantities.magnetic_axis.r
        zmaxis[i] = eq.time_slice[i].global_quantities.magnetic_axis.z
        #total toroidal magnetic field at the magnetic axis
        b_field_tor[i] = eq.time_slice[i].global_quantities.magnetic_axis.b_field_tor
        psi_axis[i] = eq.time_slice[i].global_quantities.psi_axis
        psi_bdry[i] = eq.time_slice[i].global_quantities.psi_boundary
        q95[i] = eq.time_slice[i].global_quantities.q_95
        qaxis[i] = eq.time_slice[i].global_quantities.q_axis
        li3[i] = eq.time_slice[i].global_quantities.li_3
        rminor[i] = eq.time_slice[i].boundary.minor_radius

        outline_r[i] = eq.time_slice[i].boundary.outline.r
        outline_z[i] = eq.time_slice[i].boundary.outline.z

        psi[i] = eq.time_slice[i].profiles_1d.psi
        qpsi[i] = eq.time_slice[i].profiles_1d.q
        pres[i] = eq.time_slice[i].profiles_1d.pressure
        torflux_phi[i] = eq.time_slice[i].profiles_1d.phi
        elong[i] = eq.time_slice[i].profiles_1d.elongation
        triang_top[i] = eq.time_slice[i].profiles_1d.triangularity_upper
        triang_bot[i] = eq.time_slice[i].profiles_1d.triangularity_lower
        try:
            fdp_dspi[i] = eq.time_slice[i].profiles_1d.f_df_dpsi
            dpressure_dpsi[i] = eq.time_slice[i].profiles_1d.dpressure_dpsi
        except (AttributeError, IndexError, TypeError, ValueError):
            pass
        try:
            rho_tor_norm[i] = eq.time_slice[i].profiles_1d.rho_tor_norm
            jtor[i] = eq.time_slice[i].profiles_1d.j_tor
        except (AttributeError, IndexError, TypeError, ValueError):
            pass

        if i==0:
         print('prof2d grid', eq.time_slice[i].profiles_2d[0].grid.dim1, eq.time_slice[i].profiles_2d[0].grid.dim2)
        if prof2d:
            rgrid[i] = eq.time_slice[i].profiles_2d[0].grid.dim1
            zgrid[i] = eq.time_slice[i].profiles_2d[0].grid.dim2
            psirz[i] = eq.time_slice[i].profiles_2d[0].psi

            rinv[i] = eq.time_slice[i].coordinate_system.r
            zinv[i] = eq.time_slice[i].coordinate_system.z

    plot_Ip_vacBR_loopV(times, currents, vac_b0, vac_r0, pre=pre)
    plot_eq_global(times, rmaxis, zmaxis, b_field_tor, psi_axis,
       psi_bdry, q95, qaxis, li3, pre=pre)
    plot_eq_outline(times, rminor, outline_r, outline_z, pre=pre)
    plot_equilibrium_1d(times, psi, qpsi, pres, torflux_phi, rho_tor_norm,
      elong, triang_top, triang_bot, jtor, fdp_dspi, dpressure_dpsi, pre=pre, ni=ni)
    if prof2d:
        ni = 1
        image_clr = False #plot image on top of contour
        n_ctr = 20 #number of contour lines
        plot_equilibrium_2d(times, rgrid, zgrid, rinv, zinv, psirz, pre=pre, ni=ni,
            image2d=image_clr, ncontour=n_ctr)
        plot_equilibrium_2d(times, rgrid, zgrid, rinv, zinv, psirz, rz_range=True,
            pre=pre, ni=ni, image2d=image_clr,ncontour=n_ctr)

    fname = pre+'_'+name if pre else name
    finalize_plots(fname)


def plot_core_profiles(cp, pre:str=None, name='cp', save:bool=True):
    ntimes = 0
    try:
        ntimes = len(cp.time)
    except (AttributeError, TypeError):
        print('core_profiles not available')
        return

    pre = pre if pre else 'cp'
    times = cp.time
    print(times[:3], '...', times[-3:])
    currents = cp.global_quantities.ip
    vloop = cp.global_quantities.v_loop
    vac_r0 = cp.vacuum_toroidal_field.r0 #single value
    vac_b0 = cp.vacuum_toroidal_field.b0
    plot_Ip_vacBR_loopV(times, currents, vac_b0, vac_r0, vloop, pre=pre)
    fname = pre+'_'+name if pre else name
    finalize_plots(fname)


def plot_limiter_ids(wall, pre:str=None, name='limiter', save:bool=True):
    pre = pre if pre else ''
    nunits = 0
    try:
        nunits = len(wall.description_2d[0].limiter.unit)
    except (AttributeError, IndexError, TypeError):
        print('limiters not available')
        return

    print('homogeneous_time',wall.ids_properties.homogeneous_time)
    #should be 2: upper/lower
    nunits = len(wall.description_2d[0].limiter.unit)
    for i in range(nunits):
        cl = wall.description_2d[0].limiter.unit[i].closed
        lr = wall.description_2d[0].limiter.unit[i].outline.r
        lz = wall.description_2d[0].limiter.unit[i].outline.z
        print('unit',i, 'closed',cl)
        plot_limiters(lr, lz, pre+'limiters'+str(i))
    fname = pre+'_'+name if pre else name
    finalize_plots(fname)

def plot_thomson_data(radr, radz, rad_phi, dr, dz, dphi, dtimes, dsmd,
      dsmt, ted, tet, ned, net, fname='ts', save:bool=True, pre:str=None, ni:int=0):
    nchannels = len(radr)
    name = 'thomson_pos'
    name = pre+'_'+name if pre else name
    fig = plt.figure(tight_layout=True)
    gsp = gridspec.GridSpec(2, 2)
    channels = np.arange(nchannels)
    fig = plot_subgrid(fig, gsp[0], channels, radr, 'channels', 'radiusR')
    fig = plot_subgrid(fig, gsp[1], channels, radz, 'channels', 'radiusZ')
    fig = plot_subgrid(fig, gsp[2], channels, rad_phi, 'channels', 'radiusPhi')
    finalize_plot(plt, fig, name, save=True)

    #2d
    nr, nc = 2,2
    if ni > 0:
        nr, nc = find_nrows_cols(ni)
    name = 'thomson_delta'
    name = pre+'_'+name if pre else name

    plot_data2d_samex(name, nr, nc, dtimes, 'times', (dr, 'dr'),
       (dz, 'dz'), (dphi, 'dphi'))
    name = 'thomson'
    name = pre+'_'+name if pre else name
    plot_data2d_samex(name, nr, nc, tet, 'times', (dsmd, 'sepx_dist'),
       (ted, 't_e'), (ned, 'n_e'))
    fname = pre+'_'+fname if pre else fname
    finalize_plots(fname)

def plot_thomson_scattering_ids(ts, pre:str=None, name='ts', save:bool=True):
    ntimes, nch = 0, 0
    times, ntimes, pre = check_valid(pf, pre, 'thomson_scattering')
    nch = len(ts.channel)

    ndelta = len(ts.channel[0].delta_position.time)
    ndsm = len(ts.channel[0].distance_separatrix_midplane.time)
    nte = len(ts.channel[0].t_e.time)
    nne = len(ts.channel[0].n_e.time)
    if ndelta != ntimes or ndsm != ntimes or nte != ntimes or nne != ntimes:
        raise Exception(f'ndelta {ndelta} ndsm {ndsm} nte {nte} nne {nne}')
    radr, radz, rad_phi = declare_np_arrays(3, nch)
    dr, dz, dphi, dtimes = declare_np_arrays(4, (nch,ndelta))
    dsmd, dsmt = declare_np_arrays(2, (nch,ndsm))
    ted, tet, ned, net, neid = declare_np_arrays(5, (nch,nne))

    for i in range(nch):
        #1d
        radr[i] = ts.channel[i].position.r
        radz[i] = ts.channel[i].position.z
        rad_phi[i] = ts.channel[i].position.phi

        #2d
        dr[i] = ts.channel[i].delta_position.r
        dz[i] = ts.channel[i].delta_position.z
        dphi[i] = ts.channel[i].delta_position.phi
        dtimes[i] = ts.channel[i].delta_position.time
        dsmd[i] = ts.channel[i].distance_separatrix_midplane.data
        dsmt[i] = ts.channel[i].distance_separatrix_midplane.time
        ted[i] = ts.channel[i].t_e.data
        tet[i] = ts.channel[i].t_e.time

        ned[i] = ts.channel[i].n_e.data
        net[i] = ts.channel[i].n_e.time

    plot_thomson_data(radr, radz, rad_phi, dr, dz, dphi, dtimes,
      dsmd, dsmt, ted, tet, ned, net, name)



def  plot_pf_data(times, current, voltage, name, pre=''):
    name = pre+'_'+name if pre else name
    nr, nc = 2, 2

    #plot_data2d_samex(name, nr, nc, times, 'times', (current, 'current'),
    # (voltage, 'voltage'))


def plot_PFActive(pf, name='pfa', pre:str=None, save:bool=True):
    times, ntimes, pre = check_valid(pf, pre, 'PFActive')
    if ntimes == 0:
        return
    ncoils = len(pf.coil)
    current, voltage = declare_np_arrays(2, (ncoils, ntimes))
    for i in range(ncoils):
        current[i] = pf.coil[i].current.data
        #ct = pf.coil[i].current.time
        #nts = pf.coil[i].element[0].turns_with_sign
        #geomt = pf.coil[i].element[0].geometry.geometry_type
        voltage[i] = pf.coil[i].voltage.data

    plot_pf_data(times, current, voltage, name, pre)
    fname = pre+'_'+name if pre else name
    finalize_plots(fname)

def plot_CHERS(cx, pre:str=None, save:bool=True):
    pre = pre if pre else ''


def plot_NBI(nb, pre:str=None, save:bool=True):
    pre = pre if pre else ''

def test_plot():
    def func3(x, y):
      return (1 - x / 2 + x**5 + y**3) * np.exp(-(x**2 + y**2))
    dx, dy = 0.05, 0.05
    x = np.arange(-3.0, 3.0, dx)
    y = np.arange(-3.0, 3.0, dy)
    X, Y = np.meshgrid(x, y)
    Z2 = func3(X, Y)
    plot_equilibrium_2d(None, x, y, Z2)

if __name__ == '__main__':
    print('please use: python examples/plot_ids.py')
