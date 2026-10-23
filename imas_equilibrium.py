#!/usr/bin/python
"""Populate the IMAS equilibrium IDS, including COCOS-aware signs and TRANSP geometry."""


from nc2imas import *
#from multiprocessing import Pool
from datetime import datetime
import matplotlib.pyplot as plt #boundary trace
from matplotlib import path
from matplotlib.collections import LineCollection


class imas_equilibrium:
    def __init__(self, ncimas:nc2imas, conv_moments=False, nzones_in=0):
        self.im = ncimas
        self.logger = ncimas.logger
        self. debug_max = self.im.debug_max
        self.rmcb = []
        self.rmsb = []
        self.ymcb = []
        self.ymsb = []
        self.conv_moments = conv_moments
        self.nmoments = 0
        self.nzones_in = nzones_in
        self.coskth = None
        self.sinkth = None
        self.dtheta = None
        self.dthetab = None
        self.inv_grid_theta = None
        time1 = datetime.now()
        self.bound_error = 0
        self.create_ids()
        elapsed = (datetime.now() - time1).seconds/60.0
        if elapsed > 1.0:
            print(f'\ttook:{elapsed:.2f} m')

    #to be called per time frame
    def convert_moments_2(self, ncoord_xb, ntheta:int, time:int):
         #TODO check dimension of fill data
         nmoments = self.im.get_var_names_re('RMC?', nd1=1, nd2=2, length=True)
         rcos = np.zeros(shape=(nmoments, ncoord_xb))
         rsin = np.zeros(shape=(nmoments, ncoord_xb))
         ycos = np.zeros(shape=(nmoments, ncoord_xb))
         ysin = np.zeros(shape=(nmoments, ncoord_xb))

         cm2m = 0.01
         #RMC00(TIME3, XB); CM; 0TH ASYM R MOMENT
         rcos[0] = self.im.get_ncvar(f'RMC00', ind=time, scale=cm2m)
         ycos[0] = self.im.get_ncvar(f'YMC00', ind=time, scale=cm2m)

         pi = self.im.pi

         for k in range(1, nmoments):
            kn = f'0{k}' if k<10 else k
            #RMCn(TIME3, XB); CM ; nth ASYM R COS MOMENT
            rcos[k] = self.im.get_ncvar(f'RMC{kn}', ind=time, scale=cm2m)
            rsin[k] = self.im.get_ncvar(f'RMS{kn}', ind=time, scale=cm2m)
            ycos[k] = self.im.get_ncvar(f'YMC{kn}', ind=time, scale=cm2m)
            ysin[k] = self.im.get_ncvar(f'YMS{kn}', ind=time, scale=cm2m)

         #store cos, sin
         if self.coskth is None or self.sinkth is None:
            self.coskth = np.zeros(shape=(ntheta, nmoments))
            self.sinkth = np.zeros(shape=(ntheta, nmoments))
            for ith in range(ntheta):
               th = -self.im.pi + ith*self.dtheta
            self.coskth[ith] = [np.cos(k*th) for k in range(nmoments)]
            self.sinkth[ith] = [np.sin(k*th) for k in range(nmoments)]

         rcos = rcos.T
         rsin = rsin.T
         ycos = ycos.T
         ysin = ysin.T

         rsurf = np.zeros(shape=(ncoord_xb, ntheta))
         ysurf = np.zeros(shape=(ncoord_xb, ntheta))
         if nmoments>1:
            for j in range(ncoord_xb):
               rsurf[j] = [np.sum(rcos[j]*self.coskth[ith] + rsin[j]*self.sinkth[ith]) for ith in range(ntheta)]
               ysurf[j] = [np.sum(ycos[j]*self.coskth[ith] + ysin[j]*self.sinkth[ith]) for ith in range(ntheta)]

         return np.array(rsurf), np.array(ysurf)


    def convert_moments(self, ncoord_xb, ntheta:int, time:int):
         nmoments = self.im.get_var_names_re('RMC?', nd1=1, nd2=2, length=True)
         if nmoments < 2:
            self.logger.warning(
                'equilibrium: moment conversion requested but fewer than two RMC harmonics are available')
            return [], [], []
         rcos = np.zeros(shape=(nmoments, ncoord_xb))
         rsin = np.zeros(shape=(nmoments, ncoord_xb))
         ycos = np.zeros(shape=(nmoments, ncoord_xb))
         ysin = np.zeros(shape=(nmoments, ncoord_xb))
         # RAXIS/YAXIS are scalar TIME histories, but the Fourier surfaces
         # RMCxx/YMCxx belong to profile TIME3.  Use the TIME3 sample time.
         raxis = self.im.get_ncvar_on_time('RAXIS', scale=0.01)[time]
         yaxis = self.im.get_ncvar_on_time('YAXIS', scale=0.01)[time]

         inv_grid_theta = np.zeros(ntheta)
         inv_grid_theta[ntheta-1] = self.im.pi
         if self.dtheta is None:
            self.dtheta = 2*self.im.pi/ntheta
         cm2m = 0.01
         #RMC00(TIME3, XB); CM; 0TH ASYM R MOMENT
         rcos[0,:] = self.im.get_ncvar(f'RMC00', ind=time, scale=cm2m)
         ycos[0,:] = self.im.get_ncvar(f'YMC00', ind=time, scale=cm2m)

         for k in range(1, nmoments):
            kn = f'0{k}' if k<10 else k
            #RMCn(TIME3, XB); CM ; nth ASYM R COS MOMENT
            rcos[k,:] = self.im.get_ncvar(f'RMC{kn}', ind=time, scale=cm2m)
            rsin[k,:] = self.im.get_ncvar(f'RMS{kn}', ind=time, scale=cm2m)
            ycos[k,:] = self.im.get_ncvar(f'YMC{kn}', ind=time, scale=cm2m)
            ysin[k,:] = self.im.get_ncvar(f'YMS{kn}', ind=time, scale=cm2m)

         rsurf = np.zeros(shape=(ncoord_xb, ntheta))
         ysurf = np.zeros(shape=(ncoord_xb, ntheta))

         axis_err = 0

         if nmoments>1:
            for j in range(ncoord_xb):
               rsurf[j,:] = rcos[0,j]
               ysurf[j,:] = ycos[0,j]
               for ith in range(ntheta-1):
                  theta_target = -self.im.pi + 2*ith*self.im.pi/(ntheta-1)
                  inv_grid_theta[ith] = theta_target
                  # DD4/COCOS17 polar theta is clockwise.  TRANSP Fourier
                  # moments use the geometrical CCW angle, so evaluate the
                  # moments at the corresponding raw angle.
                  th_raw = self.im.theta_raw_from_target * theta_target
                  for k in range(1, nmoments):
                     ct = np.cos(k*th_raw)
                     st = np.sin(k*th_raw)
                     rsurf[j,ith] += rcos[k,j]*ct + rsin[k,j]*st
                     ysurf[j,ith] += ycos[k,j]*ct + ysin[k,j]*st

               rsurf[j,ntheta-1] = rsurf[j,0]
               ysurf[j,ntheta-1] = ysurf[j,0]

               # verify
               pts = np.zeros((ntheta,2))

               for k in range(ntheta):
                  pts[k,0] = rsurf[j,k]
                  pts[k,1] = ysurf[j,k]

               pth = path.Path(pts)
               if not pth.contains_point((raxis,yaxis)):
                  axis_err += 1

            #DEBUG all
            if self.debug_max and axis_err>0:
               data = np.stack((rsurf,ysurf), axis=2)
               fig, ax = plt.subplots()
               ax.add_collection(LineCollection(data))
               plt.plot(raxis, yaxis, 'ro')
               plt.show()

            if axis_err > 0:
               print('time', time, f'Error: axis {raxis}, {yaxis} outside one/more flux surfaces')
         return np.array(rsurf), np.array(ysurf), np.array(inv_grid_theta)

    #to be called per time frame
    def get_moments_bdry(self):
         ntimep = self.im.ntimep
         if self.im.is_valid(self.rmcb, self.ymcb, ind=ntimep-1):
            return
         nmoments = self.im.get_var_names_re('RMCB?', nd1=1, nd2=2, length=True)
         if nmoments < 1:
            self.logger.info('equilibrium: no LCFS boundary moments in CDF')
            return
         self.nmoments = nmoments
         # These boundary Fourier coefficients are TRANSP *TIME* histories,
         # whereas equilibrium/time_slice uses TIME3.  Resample each original
         # time series before indexing it in the TIME3 time-slice loop.
         self.rmcb = np.zeros((ntimep, nmoments))
         self.rmsb = np.zeros((ntimep, nmoments))
         self.ymcb = np.zeros((ntimep, nmoments))
         self.ymsb = np.zeros((ntimep, nmoments))
         for k in range(nmoments):
            required = [(f'RMCB{k}', self.rmcb), (f'YMCB{k}', self.ymcb)]
            if k > 0:
               required += [(f'RMSB{k}', self.rmsb), (f'YMSB{k}', self.ymsb)]
            for name, target in required:
               values = self.im.get_ncvar_on_time(name, scale=0.01)
               if len(values) != ntimep:
                  self.logger.warning(f'equilibrium: {name} missing/invalid; '
                                      'skipping incomplete LCFS boundary')
                  self.rmcb, self.rmsb, self.ymcb, self.ymsb = [], [], [], []
                  return
               target[:, k] = values

    #to be called per time frame
    def convert_moments_bdry(self, ntheta:int, time:int):
         if not self.im.is_valid(self.rmcb, self.ymcb, ind=time):
            self.get_moments_bdry()
         if not self.im.is_valid(self.rmcb, self.ymcb, ind=time):
            return [], []

         if ntheta<2 or time<0:
            raise Exception(f'convert_moments_bdry: ntheta {ntheta} time {time}')
         rsurf = np.zeros(ntheta)
         ysurf = np.zeros(ntheta)
         if self.dthetab is None:
            raise Exception('dthetab not defined')

         for ith in range(ntheta-1):
            th = ith*self.dthetab
            th = -self.im.pi + th
            rsurf[ith] = self.rmcb[time,0]
            ysurf[ith] = self.ymcb[time,0]
            for k in range(1, self.nmoments):
               ct = np.cos(k*th)
               st = np.sin(k*th)
               rsurf[ith] += self.rmcb[time,k]*ct + self.rmsb[time,k]*st
               ysurf[ith] += self.ymcb[time,k]*ct + self.ymsb[time,k]*st

         rsurf[ntheta-1] = rsurf[0]
         ysurf[ntheta-1] = ysurf[0]
         return rsurf, ysurf


    #rgrid, zgrid are for all times (time, grid)
    def trace_boundary(self, psirz, rgrid, zgrid, time, nr, nz):
         if nr < 2 or nz < 2 or time < 0 or len(psirz) != nr*nz:
            return [], []
         psib_all = self.im.get_ncvar_on_time(
             'PLFLXA', scale=self.im.psi_sign*2.*self.im.pi)
         if not self.im.is_valid(psib_all):
            psib_all = self.im.get_ncvar_on_time(
                'PBPMHD', scale=self.im.psi_sign*2.*self.im.pi)
         raxis_all = self.im.get_ncvar_on_time('RAXIS', scale=0.01)
         yaxis_all = self.im.get_ncvar_on_time('YAXIS', scale=0.01)
         if any(len(values) <= time for values in (psib_all, raxis_all, yaxis_all)):
            return [], []

         def at_time(values):
            values = np.asarray(values, dtype=float)
            return values[time] if values.ndim > 1 else values

         r = at_time(rgrid)
         z = at_time(zgrid)
         if r.ndim != 1 or z.ndim != 1 or len(r) != nr or len(z) != nz:
            return [], []
         psi2d = np.asarray(psirz, dtype=float).reshape((nz, nr))
         psib = float(psib_all[time])
         if not np.isfinite(psib):
            return [], []

         # Multiple segments can occur near X-points.  Keep the largest closed
         # contour at the LCFS flux that encloses the magnetic axis.
         from matplotlib.path import Path
         fig, ax = plt.subplots()
         try:
            contours = ax.contour(r, z, psi2d, levels=[psib])
            candidates = []
            axis = (float(raxis_all[time]), float(yaxis_all[time]))
            for contour in contours.allsegs[0]:
               points = np.asarray(contour, dtype=float)
               if points.ndim != 2 or points.shape[0] < 4 or points.shape[1] != 2:
                  continue
               scale = max(float(np.ptp(points[:, 0])), float(np.ptp(points[:, 1])), 1.0)
               closed = np.linalg.norm(points[0] - points[-1]) <= max(0.01, 1.e-3*scale)
               if not closed:
                  continue
               points[-1] = points[0]
               if not Path(points).contains_point(axis):
                  continue
               area = 0.5*abs(np.sum(points[:-1, 0]*points[1:, 1]
                                      - points[1:, 0]*points[:-1, 1]))
               candidates.append((area, points))
            if not candidates:
               self.logger.warning('equilibrium: could not trace LCFS from PSIRZ at time index %d', time)
               return [], []
            points = max(candidates, key=lambda item: item[0])[1]
            return points[:, 0].copy(), points[:, 1].copy()
         except (ValueError, TypeError, IndexError, AttributeError) as exc:
            self.logger.warning('equilibrium: LCFS contour failed at time index %d: %s', time, exc)
            return [], []
         finally:
            plt.close(fig)

    @staticmethod
    def _outline_centroid(r, z):
         """Return the area centroid of a closed R/Z outline, or None."""
         points = np.column_stack((np.asarray(r, dtype=float), np.asarray(z, dtype=float)))
         if points.ndim != 2 or points.shape[0] < 4 or not np.isfinite(points).all():
            return None
         if not np.allclose(points[0], points[-1]):
            points = np.vstack((points, points[0]))
         cross = points[:-1, 0]*points[1:, 1] - points[1:, 0]*points[:-1, 1]
         twice_area = np.sum(cross)
         if not np.isfinite(twice_area) or abs(twice_area) < 1.e-12:
            return None
         rcent = np.sum((points[:-1, 0] + points[1:, 0])*cross)/(3.0*twice_area)
         zcent = np.sum((points[:-1, 1] + points[1:, 1])*cross)/(3.0*twice_area)
         return (float(rcent), float(zcent)) if np.isfinite([rcent, zcent]).all() else None


    def boundary_minor(self, rsurf, time):
         """Boundary minor radius, or ``None`` when the boundary R is empty (recorded as skipped)."""
         values = np.asarray(rsurf, dtype=float)
         if values.size == 0:
            if time==0:
               self.im.logger.info('minor_axis: plasma boundary invalid')
            self.im.mapping.skipped('plasma boundary R is empty (no LCFS data for this time slice)',
                                    target='time_slice/boundary/minor_radius')
            return None
         return (np.max(values) - np.min(values))/2.0


    def create_ids(self):
        self.im.logger.info('\nequilibrium')
        eqm = self.im.new_ids('equilibrium')
        cm2m = 0.01
        pi = constants.pi  #scipy
        #scalar time
        timep = self.im.timep
        ntimep = self.im.ntimep
        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')

        if not self.im.is_valid(coord_x, coord_xb, name=r'coord_x/xb'):
            print('\t Error: invalid coord_x/xb', coord_x, coord_xb)
            return

        nprof = coord_x.shape[1]

        ip_sign = self.im.ip_sign
        bt_sign = self.im.bt_sign
        psi_sign = self.im.psi_sign
        phi_sign = self.im.phi_sign
        q_sign = self.im.q_sign
        f_sign = self.im.f_sign
        currents = ip_sign * self.im.get_ncvar('PCURC', interp=True) #(TIME)
        raxis = self.im.get_ncvar('RAXIS', scale=cm2m, interp=True) #RAXISEQ
        yaxis = self.im.get_ncvar('YAXIS', scale=cm2m, interp=True) #YAXISEQ

        # PSI0_TR(TIME); Wb/rad ; Psi_poloidal(axis) in simulation
        psi_axis = psi_sign * self.im.get_ncvar('PSI0_TR', scale=2*pi, interp=True)

        q95 = q_sign * self.im.get_ncvar('Q95', interp=True)
        qaxis  = q_sign * self.im.get_ncvar('Q0', interp=True)
        li3 = self.im.get_ncvar('LI_3', interp=True)
        bpeq = self.im.get_ncvar('BETAT', interp=True) #poloidal beta
        btor = self.im.get_ncvar('BTEQ', interp=True) #toroidal beta
        #Zone volume
        gvol = self.im.get_ncvar('DVOL', 1.0e-6) #(TIME3,X) cm^3
        volume_enclosed = self.im.cumulative_zone_quantity(gvol)
        #plasma volume
        pvol = self.im.get_ncvar('PVOL', 1.0e-6, interp=True) #(TIME) cm^3
        if not self.im.is_valid(gvol, pvol, name='vol'):
            return

        #PMHD(TIME3, X). Pressure from equilibrium solver
        pmhd = self.im.get_ncvar('PMHD')
        emhd = self.im.np_multiply_sum(pmhd, gvol, axis=1, scale=1.5,
               name1='pmhd', name2='gvol')

        # beta_tor_norm is constructed below from scalar histories first
        # interpolated onto TIME3; no additional time axis is introduced.
        beta_tor_norm = []
        bz = self.im.get_ncvar('BZ', interp=True) #(TIME) T
        # PAREA(TIME) - plasma CS area[cm^2]
        area = self.im.get_ncvar('PAREA', scale=1.0e-4, interp=True)
        try:
            # 100 * beta_tor * a[m] * B0 [T] / ip [MA]   ???
            a = np.sqrt(area/pi)
            beta_tor_norm = 1.0e8*btor*a*(bt_sign*bz)/currents
        except (ValueError, TypeError, FloatingPointError):
            pass

        #Area of the LCFS poloidal cross section
        #ZONE CROSS SECTIONAL AREA;(TIME3,X) cm^2
        darea = self.im.get_ncvar('DAREA', 1.0e-4)
        area_enclosed = self.im.cumulative_zone_quantity(darea)

        #Surface area of the toroidal flux surface
        #FLUX SURFACE AREA;(TIME3,X) cm^2
        surface = self.im.get_ncvar('SURF', 1.0e-4)
        #GB1(TIME3, X) ; T; <B> flux surface average
        bf_ave = self.im.get_ncvar('GB1')
        #BMAX(TIME3, X) ; T; Bmax on flux surface
        bf_max = self.im.get_ncvar('BMAX')
        #BMIN(TIME3, X) ; T; Bmin on flux surface
        bf_min = self.im.get_ncvar('BMIN')

        # GB2(TIME3, X);  T**2; <B**2> flux surface average
        fs_ave = self.im.get_ncvar('GB2')
        # GRI(TIME3, X);CM**-1; <1/R> FLUX SURFACE VOL.AVG
        fs_vol_ave1 = self.im.get_ncvar('GRI', scale=100)
        # GR2I(TIME3, X) ;CM**-2; <1/R**2> FLUX SURFACE VOL.AVG
        fs_vol_ave2 = self.im.get_ncvar('GR2I', scale=1.0e4)
        # Flux surface averaged |grad_rho_tor| =  GXI * DRAV / DRAVFAC
        # GXI(TIME3, X) ;CM**-1; <GRAD(XI)> FLUX SURF VOL.AVG
        gxi_fs = self.im.get_ncvar('GXI', scale=10)
        # DRAV(TIME3, X); CM; FLUX SURFACE AVG <DR>
        drav_fs = self.im.get_ncvar('DRAV', scale=0.1)
        # DRAVFAC(TIME3, X); ; <dr>*<1/dr>
        drav_fac = self.im.get_ncvar('DRAVFAC')
        fs_ave_gm = self.im.operate_arrays(gxi_fs, '*', drav_fs)
        fs_ave_gm7 = self.im.operate_arrays(fs_ave_gm, '/', drav_fac)

        #rminor =

        #profiles (TIME3, XB)
        elong = self.im.get_ncvar('ELONG')
        triang = self.im.get_ncvar('TRIANG')
        triang_lo = self.im.get_ncvar('TRIANGL')
        triang_up = self.im.get_ncvar('TRIANGU')
        square_up = self.im.get_ncvar('SQUARE_UO')
        square_lo = self.im.get_ncvar('SQUARE_LO')

        # Boundary shaping is already on TIME3: the last XB surface of each
        # profile.  Do not interpolate backwards to scalar TIME, because these
        # values are written in equilibrium/TIME3 time_slice objects.
        elongb = elong[:, -1] if self.im.is_valid(elong) else []
        triangb = triang[:, -1] if self.im.is_valid(triang) else []
        triang_lob = triang_lo[:, -1] if self.im.is_valid(triang_lo) else []
        triang_upb = triang_up[:, -1] if self.im.is_valid(triang_up) else []

        bdry_sepx_psi = psi_sign * self.im.get_ncvar('PBPMHD', 2.0*pi, interp=True)
        psi = self.im.get_ncvar('PLFLX2PI')
        if not self.im.is_valid(psi):
            psi = self.im.get_ncvar('PLFLX', 2.*pi) #(TIME3, XB)
        if self.im.is_valid(psi):
            psi = psi_sign * psi

        # add psi at axis PSI0_TR(TIME)
        if self.im.is_valid(psi, psi_axis, ind=ntimep-1):
            psi = psi + psi_axis[:,None]

        # TRFLX is already the full toroidal flux in Wb (not Wb/rad).
        phi = phi_sign * self.im.get_ncvar('TRFLX')
        qpsi = q_sign * self.im.get_ncvar('Q')
        pres_x = self.im.get_ncvar('PTOWB') # PPLAS, on X
        fbz = self.im.get_ncvar('GFUN')

        nprof1d = 0
        try:
            nprof1d = psi.shape[1]
        except (AttributeError, IndexError, TypeError):
            print(' Error getting equilibrium psi')
            return

        #(TIME) TESLA*CM. vac field
        bzxR = self.im.get_ncvar('BZXR', cm2m, interp=True)
        b_field_tor = bt_sign * self.im.operate_arrays(bzxR, '/', raxis)

        # Zeroth Fourier coefficients are the geometric centers of the flux
        # surfaces.  RMJMP remains a radial fallback for CDFs without moments.
        geom_axisr = self.im.get_ncvar('RMC00', cm2m)
        if not self.im.is_valid(geom_axisr):
            geom_axisr = self.im.get_ncvar('RMJMP', cm2m)
        geom_axisz = self.im.get_ncvar('YMC00', cm2m)

        #f
        fprof = []
        try:
            self.logger.info(f'nprof1d {nprof1d} bzxR {bzxR.shape} fbz {fbz.shape}')
            fprof = f_sign * fbz * bzxR[:,None] # signed R*B_phi
        except (AttributeError, ValueError, TypeError):
            pass

        # equilibrium profiles_1d is naturally a flux-surface (XB) grid because
        # psi, phi, q, F, shaping and SURF are all TRANSP boundary quantities.
        # Map quantities stored on X zone centres to XB before writing the IDS.
        pres = self.im.interp_radial_profiles(pres_x, coord_x, coord_xb)
        bf_ave = self.im.interp_radial_profiles(bf_ave, coord_x, coord_xb)
        bf_max = self.im.interp_radial_profiles(bf_max, coord_x, coord_xb)
        bf_min = self.im.interp_radial_profiles(bf_min, coord_x, coord_xb)
        fs_vol_ave2 = self.im.interp_radial_profiles(fs_vol_ave2, coord_x, coord_xb)
        fs_ave = self.im.interp_radial_profiles(fs_ave, coord_x, coord_xb)
        fs_ave_gm7 = self.im.interp_radial_profiles(fs_ave_gm7, coord_x, coord_xb)
        fs_vol_ave1 = self.im.interp_radial_profiles(fs_vol_ave1, coord_x, coord_xb)

        # Normalised toroidal-flux coordinate on the same XB grid.
        rho_tor_norm = coord_xb

        # Additional equilibrium profiles directly supported by the CDF/DD.
        psi_norm = []
        if self.im.is_valid(psi):
            try:
                denom = psi[:,-1] - psi_axis
                psi_norm = np.divide(psi - psi_axis[:,None], denom[:,None],
                                     out=np.zeros_like(psi),
                                     where=np.abs(denom[:,None]) > 0.0)
            except (ValueError, TypeError, IndexError, FloatingPointError) as exc:
                self.im.mapping.broken(exc, source='PLFLX', target='time_slice/profiles_1d/psi_norm')
                psi_norm = []

        # rho_tor = sqrt(Phi/(pi*B0)) with B0 the vacuum_toroidal_field/b0 written in this IDS.
        # Phi and the physical B0 have the same sign in a consistent convention, leaving a
        # positive metric rho.
        rho_tor = []
        reference = self.im.vacuum_field_reference(self.im.timep)
        b0_signed = None if reference is None else reference[1]
        if b0_signed is None:
            self.im.mapping.skipped('vacuum toroidal field (RAXIS/BZXR) is not in the CDF',
                                    source='BZXR', target='time_slice/profiles_1d/rho_tor')
        else:
            try:
                rho_tor = np.sqrt(np.clip(phi/(pi*b0_signed[:,None]), 0.0, None))
            except (TypeError, ValueError, FloatingPointError) as exc:
                self.im.mapping.broken(exc, source='TRFLX', target='time_slice/profiles_1d/rho_tor')

        dpsi_drho_tor = self.im.radial_gradient(psi, rho_tor) \
            if self.im.is_valid(psi, rho_tor) else []
        dpressure_dpsi = self.im.radial_gradient(pres, psi) \
            if self.im.is_valid(pres, psi) else []
        df_dpsi = self.im.radial_gradient(fprof, psi) \
            if self.im.is_valid(fprof, psi) else []
        f_df_dpsi = fprof * df_dpsi if self.im.is_valid(fprof, df_dpsi) else []

        # Current-density/shear/mass profiles are X-centred in TRANSP; map to
        # the XB flux-surface grid used by equilibrium/profiles_1d.
        jphi_x = ip_sign * self.im.get_ncvar('CUR', self.im.area_inv)
        jphi = self.im.interp_radial_profiles(jphi_x, coord_x, coord_xb,
                                              axis_value=jphi_x[:,0]) \
            if self.im.is_valid(jphi_x) else []
        jb_x = self.im.get_ncvar('PLJB', self.im.area_inv)
        jparallel_x = []
        # <J.B>/B0 with the IDS's own |b0| (see rho_tor above), not TRANSP's BZ
        b0_mag = None if b0_signed is None else np.abs(b0_signed)
        if b0_mag is None:
            self.im.mapping.skipped('vacuum toroidal field (RAXIS/BZXR) is not in the CDF',
                                    source='BZXR', target='time_slice/profiles_1d/j_parallel')
        elif self.im.is_valid(jb_x):
            try:
                jparallel_x = self.im.jparallel_sign * jb_x / b0_mag[:,None]
            except (TypeError, ValueError, FloatingPointError) as exc:
                self.im.mapping.broken(exc, source='PLJB', target='time_slice/profiles_1d/j_parallel')
                jparallel_x = []
        jparallel = self.im.interp_radial_profiles(jparallel_x, coord_x, coord_xb,
                                                   axis_value=jparallel_x[:,0]) \
            if self.im.is_valid(jparallel_x) else []
        shear_x = self.im.get_ncvar('SHAT')
        magnetic_shear = self.im.interp_radial_profiles(shear_x, coord_x, coord_xb,
                                                        axis_value=shear_x[:,0]) \
            if self.im.is_valid(shear_x) else []
        mass_density_x = self.im.get_ncvar('PMASS', 1.0e3)  # g/cm3 -> kg/m3
        mass_density = self.im.interp_radial_profiles(mass_density_x, coord_x, coord_xb,
                                                      axis_value=mass_density_x[:,0]) \
            if self.im.is_valid(mass_density_x) else []

        rgrid = self.im.get_ncvar('RGRID', scale=cm2m)
        zgrid = self.im.get_ncvar('ZGRID', scale=cm2m)

        rsh = rgrid.shape
        zsh = zgrid.shape
        nr = rsh[1] if len(rsh)>1 else rsh[0]
        nz = zsh[1] if len(zsh)>1 else zsh[0]
        nrz = nr*nz

        ncoord_xb = coord_xb.shape[1]

        if len(rsh)>1 and rsh[1]>0 and rgrid[0,rsh[1]-1]<0.2:
            print('Check RGRID unit:', rgrid[0])
            #raise Exception('RGRID too small:expected in cm')
        self.logger.info(f'rgrid {rgrid.shape}')

        #RSURF - flux surface R locations

        '''
         Number of poloidal grid points per flux surface, ntheta, should be
         equal to MTH in TRANSP. MTH is calculated using namelist variable NZONES,
         which is nr=the number of radial grids in the moments or profile data.
         ntheta = max(65,min(251,(2*nzones+1)))
         But user can override this by namelist variable NTH_RZ.

        '''
        nzones_ = self.nzones_in if self.nzones_in>1 else nprof1d
        ntheta = max(65,min(251,(2*nzones_+1)))

        self.dtheta = 2*pi/ntheta

        nthetab = 128
        self.dthetab = 2*pi/nthetab

        if self.inv_grid_theta is None:
           self.inv_grid_theta = [-self.im.pi + ith*self.dtheta for ith in range(ntheta)]

        eqm = self.im.fill_properties(eqm)
        self.im.fill_code_properties(eqm)
        eqm.ids_properties.homogeneous_time = 1
        # All equilibrium time slices and 1-D/2-D profiles use TIME3.
        # Scalars read from TIME above were explicitly interpolated to TIME3.
        eqm.time.resize(ntimep)
        eqm.time[:] = timep
        # j_parallel and rho_tor of this IDS are defined with this IDS's vacuum_toroidal_field/b0.
        self.im.fill_vac_tor_field(eqm)
        eqm.time_slice.resize(ntimep)

        has2d_eq = False
        self.logger.info(f'equilibrium #time slices {ntimep}')
        for it in range(ntimep):
            print_progress(ntimep-1, it)
            if it==0 or it==ntimep-1:
                self.logger.info(f'time slice: {it}')
            #global
            if len(currents) > it :
                eqm.time_slice[it].global_quantities.ip = currents[it].astype(float)
            if len(raxis) > it:
                eqm.time_slice[it].global_quantities.magnetic_axis.r = raxis[it].astype(float)
                eqm.time_slice[it].global_quantities.magnetic_axis.z = yaxis[it].astype(float)
            # Total toroidal magnetic field at the magnetic axis.
            # DD 3.42+ renamed b_field_tor -> b_field_phi.
            if len(b_field_tor) > it:
                magnetic_axis = eqm.time_slice[it].global_quantities.magnetic_axis
                b_axis = b_field_tor[it].astype(float)
                if self.im.dd_at_least(3, 42, 0):
                    magnetic_axis.b_field_phi = b_axis
                else:
                    magnetic_axis.b_field_tor = b_axis

            if len(psi_axis) > it:
                gq = eqm.time_slice[it].global_quantities
                psi_ax = psi_axis[it].astype(float)
                # DD 4.1 renamed psi_axis -> psi_magnetic_axis.
                if self.im.dd_at_least(4, 1, 0):
                    gq.psi_magnetic_axis = psi_ax
                else:
                    gq.psi_axis = psi_ax
            # Selected plasma-boundary poloidal flux.  DD4 stores the preferred
            # value under time_slice/boundary/psi; keep the deprecated global
            # node as a fallback for older dictionaries.
            if len(psi) > it:
                psi_bnd = psi[it,-1].astype(float)
                boundary = eqm.time_slice[it].boundary
                if hasattr(boundary, 'psi'):
                    boundary.psi = psi_bnd
                gq = eqm.time_slice[it].global_quantities
                if hasattr(gq, 'psi_boundary'):
                    gq.psi_boundary = psi_bnd

            if len(q95) > it:
                eqm.time_slice[it].global_quantities.q_95 = q95[it].astype(float)
            if len(qaxis) >it:
                eqm.time_slice[it].global_quantities.q_axis = qaxis[it].astype(float)
            if len(li3) > it:
                eqm.time_slice[it].global_quantities.li_3 = li3[it].astype(float)
            if len(bpeq) >it:
                eqm.time_slice[it].global_quantities.beta_pol = bpeq[it].astype(float)
            if len(btor) >it:
                eqm.time_slice[it].global_quantities.beta_tor = btor[it].astype(float)
            if len(beta_tor_norm):
                gq = eqm.time_slice[it].global_quantities
                # DD3.39 uses beta_normal; DD4 uses beta_tor_norm.
                if self.im.dd_at_least(4, 0, 0):
                    gq.beta_tor_norm = beta_tor_norm[it].astype(float)
                else:
                    gq.beta_normal = beta_tor_norm[it].astype(float)
            if len(pvol) > it:
                eqm.time_slice[it].global_quantities.volume = pvol[it].astype(float)
            if len(darea) > it:
                eqm.time_slice[it].global_quantities.area = area_enclosed[it,-1].astype(float)
            #last ??
            if len(surface) > it:
                eqm.time_slice[it].global_quantities.surface = surface[it,-1].astype(float)

            if len(emhd) >it:
                eqm.time_slice[it].global_quantities.energy_mhd = emhd[it].astype(float)

            #2D psi
            psirz = psi_sign * self.im.get_ncvar('PSIRZ', 2.*pi, ind=it)

            has_2d_slice = (self.im.is_valid(rgrid, zgrid, ind=it, name='grid')
                            and self.im.is_valid(psirz, name='prirz'))
            if has_2d_slice:
                rbound, ybound = self.trace_boundary(psirz, rgrid, zgrid, it, nr, nz)
            else:
                rbound, ybound = [], []
            # Moment-derived boundary surfaces are an explicit opt-in.  The
            # default uses PSIRZ when available and warns if it cannot trace it.
            if not self.im.is_valid(rbound, ybound) and self.conv_moments:
                rbound, ybound = self.convert_moments_bdry(nthetab, it)
            if not self.im.is_valid(rbound, ybound) and not self.conv_moments:
                self.logger.warning(
                    'equilibrium: no LCFS outline at time index %d; enable --moments yes '
                    'to use TRANSP boundary moments', it)
            self.im.fill_data(eqm.time_slice[it].boundary.outline.r, rbound, name='rbound')
            self.im.fill_data(eqm.time_slice[it].boundary.outline.z, ybound, name='ybound')

            minor = self.boundary_minor(rbound, it)
            if minor is not None and minor > 0:
               eqm.time_slice[it].boundary.minor_radius = minor

            geom_rz = None
            # RMCB0/YMCB0 are the boundary's geometric center and can be used
            # without reconstructing every Fourier harmonic.
            rmc_b0 = self.im.get_ncvar_on_time('RMCB0', scale=cm2m)
            ymc_b0 = self.im.get_ncvar_on_time('YMCB0', scale=cm2m)
            if (len(rmc_b0) > it and len(ymc_b0) > it
                    and np.isfinite([rmc_b0[it], ymc_b0[it]]).all()):
               geom_rz = (float(rmc_b0[it]), float(ymc_b0[it]))
            if geom_rz is None and self.im.is_valid(rbound, ybound):
               geom_rz = self._outline_centroid(rbound, ybound)
            if geom_rz is None and self.im.is_valid(geom_axisr, geom_axisz, ind=it):
               geom_rz = (float(geom_axisr[it, -1]), float(geom_axisz[it, -1]))
            if geom_rz is not None:
               boundary = eqm.time_slice[it].boundary
               if hasattr(boundary, 'geometric_axis'):
                  boundary.geometric_axis.r = geom_rz[0]
                  boundary.geometric_axis.z = geom_rz[1]
               gq = eqm.time_slice[it].global_quantities
               if hasattr(gq, 'geometric_axis'):
                  gq.geometric_axis.r = geom_rz[0]
                  gq.geometric_axis.z = geom_rz[1]

            if len(bdry_sepx_psi) > it:
                # In DD4 the legacy boundary_separatrix structure is removed;
                # the selected separatrix/boundary is represented by boundary.
                bpsi = bdry_sepx_psi[it].astype(float)
                if self.im.dd_at_least(4, 0, 0):
                    eqm.time_slice[it].boundary.psi = bpsi
                else:
                    eqm.time_slice[it].boundary_separatrix.psi = bpsi

            if len(triangb) > it:
                eqm.time_slice[it].boundary.triangularity = triangb[it].astype(float)
            if len(triang_upb) > it:
                eqm.time_slice[it].boundary.triangularity_upper = triang_upb[it].astype(float)
            if len(triang_lob) > it:
                eqm.time_slice[it].boundary.triangularity_lower = triang_lob[it].astype(float)

            if len(elongb) > it:
                eqm.time_slice[it].boundary.elongation = elongb[it].astype(float)

            #profiles 1d
            iprof = eqm.time_slice[it].profiles_1d
            self.im.fill_data(iprof.psi, psi, ind=it, name='psi')
            self.im.fill_data(iprof.phi, phi, ind=it, name='phi')
            self.im.fill_data(iprof.q, qpsi, ind=it, name='qpsi')
            self.im.fill_data(iprof.f, fprof, ind=it, name='fprof')
            self.im.fill_data(iprof.pressure, pres, ind=it, name='pres')
            self.im.fill_data(iprof.volume, volume_enclosed, ind=it, name='volume_enclosed')
            self.im.fill_data(iprof.area, area_enclosed, ind=it, name='area_enclosed')
            self.im.fill_data(iprof.surface, surface, ind=it, name='surface')
            self.im.fill_data(iprof.b_field_average, bf_ave, ind=it, name='bf_ave')
            self.im.fill_data(iprof.b_field_max, bf_max, ind=it, name='bf_max')
            self.im.fill_data(iprof.b_field_min, bf_min, ind=it, name='bf_min')
            self.im.fill_data(iprof.gm1, fs_vol_ave2, ind=it, name='gm1')
            self.im.fill_data(iprof.gm5, fs_ave, ind=it, name='gm5')
            self.im.fill_data(iprof.gm7, fs_ave_gm7, ind=it, name='gm7')
            self.im.fill_data(iprof.gm9, fs_vol_ave1, ind=it, name='gm9')

            self.im.fill_data(iprof.rho_tor_norm, rho_tor_norm, ind=it, name='rho_tor_norm')
            if hasattr(iprof, 'psi_norm'):
                self.im.fill_data(iprof.psi_norm, psi_norm, ind=it, name='psi_norm')
            if hasattr(iprof, 'rho_tor'):
                self.im.fill_data(iprof.rho_tor, rho_tor, ind=it, name='rho_tor')
            if hasattr(iprof, 'dpsi_drho_tor'):
                self.im.fill_data(iprof.dpsi_drho_tor, dpsi_drho_tor, ind=it, name='dpsi_drho_tor')
            if hasattr(iprof, 'dpressure_dpsi'):
                self.im.fill_data(iprof.dpressure_dpsi, dpressure_dpsi, ind=it, name='dpressure_dpsi')
            if hasattr(iprof, 'f_df_dpsi'):
                self.im.fill_data(iprof.f_df_dpsi, f_df_dpsi, ind=it, name='f_df_dpsi')
            if self.im.dd_at_least(3, 42, 0) and hasattr(iprof, 'j_phi'):
                self.im.fill_data(iprof.j_phi, jphi, ind=it, name='j_phi')
            elif hasattr(iprof, 'j_tor'):
                self.im.fill_data(iprof.j_tor, jphi, ind=it, name='j_tor')
            if hasattr(iprof, 'j_parallel'):
                self.im.fill_data(iprof.j_parallel, jparallel, ind=it, name='j_parallel')
            if hasattr(iprof, 'magnetic_shear'):
                self.im.fill_data(iprof.magnetic_shear, magnetic_shear, ind=it, name='magnetic_shear')
            if hasattr(iprof, 'mass_density'):
                self.im.fill_data(iprof.mass_density, mass_density, ind=it, name='mass_density')
            self.im.fill_data(iprof.elongation, elong, ind=it, name='elong')
            self.im.fill_data(iprof.squareness_upper_outer, square_up, ind=it, name='square_up')
            self.im.fill_data(iprof.squareness_lower_outer, square_lo, ind=it, name='square_lo')
            self.im.fill_data(iprof.triangularity_upper, triang_up, ind=it, name='triang_up')
            self.im.fill_data(iprof.triangularity_lower, triang_lo, ind=it, name='triang_lo')
            self.im.fill_data(iprof.geometric_axis.r, geom_axisr, ind=it, name='geom_axisr')
            self.im.fill_data(iprof.geometric_axis.z, geom_axisz, ind=it, name='geom_axisz')

            #self.im.fill_data(iprof.j_tor, jtor, ind=it, name='jtor')
            #self.im.fill_data(iprof.f_df_dpsi, fdp_dspi, ind=it, name='fdp_dspi')
            #self.im.fill_data(iprof.dpressure_dpsi, dpressure_dpsi, ind=it, name='dpressure_dpsi')

            eqm.time_slice[it].profiles_2d.resize(1)
            prof2 = eqm.time_slice[it].profiles_2d[0]
            #prof2.type.index = 0
            prof2.grid_type.index = 1
            self.im.fill_data(prof2.grid.dim1, rgrid, ind=it, name='rgrid')
            self.im.fill_data(prof2.grid.dim2, zgrid, ind=it, name='zgrid')

            if self.im.is_valid(rgrid, zgrid, ind=it, name='grid') and self.im.is_valid(psirz, name='prirz'):
                prof2.psi.resize((nr, nz))
                prof2.psi = np.reshape(psirz, (nz, nr)).T
                if not has2d_eq:
                  has2d_eq = True

                if it==ntimep-1:
                    self.logger.info(f't {it}: wrote 2d profile : rgrid {nr}, zgrid {nz} psi {prof2.psi.shape}')

            bfield = self.im.get_ncvar('B_FIELD', ind=it)

            if self.im.is_valid(bfield, name='bfield'):
                bfield = bfield.reshape((3,nz,nr))
                # TRANSP field-orientation controls apply differently to the
                # poloidal and toroidal components.  B_R/B_Z follow Ip, while
                # B_phi follows the toroidal-field direction.
                b_r = self.im.poloidal_field_sign * bfield[0]
                b_phi = bt_sign * bfield[1]
                b_z = self.im.poloidal_field_sign * bfield[2]
                prof2.b_field_r.resize((nr, nz))
                prof2.b_field_r = b_r.T
                # DD3.39 uses b_field_tor; DD4.1.1 uses b_field_phi.
                if self.im.dd_at_least(3, 42, 0):
                    prof2.b_field_phi.resize((nr, nz))
                    prof2.b_field_phi = b_phi.T
                else:
                    prof2.b_field_tor.resize((nr, nz))
                    prof2.b_field_tor = b_phi.T
                prof2.b_field_z.resize((nr, nz))
                prof2.b_field_z = b_z.T
                if it==ntimep-1:
                    self.logger.info(f't {it}: wrote B_FIELD {bfield.shape}')

            #inverse equilibrium
            coord = eqm.time_slice[it].coordinate_system

            #inverse coordinates
            if self.conv_moments:
               rsurf, ysurf, inv_theta = self.convert_moments(ncoord_xb, ntheta, it)

               #fill data
               # dim1 is TRANSP XB = sqrt(Phi/Phi_boundary), not a metric rho.
               # Use the explicit normalized-toroidal-flux polar identifier.
               coord.grid_type.name = 'inverse_rhotornorm_polar'
               coord.grid_type.index = 33
               self.im.fill_data(coord.grid.dim1, coord_xb, ind=it, name='grid_r')
               #self.im.fill_data(coord.grid.dim2, self.inv_grid_theta, name='grid_theta')
               self.im.fill_data(coord.grid.dim2, inv_theta, name='grid_theta')
               coord.r.resize(ncoord_xb, ntheta)
               if self.im.is_valid(rsurf):
                  coord.r = rsurf
               coord.z.resize(ncoord_xb, ntheta)
               if self.im.is_valid(ysurf):
                  coord.z = ysurf

        self.logger.info('writing ids')
        self.im.putdb_imas(eqm)

if __name__ == '__main__':
    print("Please use: python convert.py ")

