
from nc2imas import *
from multiprocessing import Pool
from datetime import datetime
import matplotlib.pyplot as plt #boundary trace
from matplotlib import path

class imas_equilibrium:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.rmcb = []
        self.rmsb = []
        self.ymcb = []
        self.ymsb = []
        self.nmoments = 0
        self.coskth = None
        self.sinkth = None
        self.dtheta = None
        self.dthetab = None
        self.inv_grid_theta = None
        time1 = datetime.now()
        self.bound_error = 0
        self.create_ids()
        elapsed = (datetime.now() - time1).seconds/60.0
        print(f'took:{elapsed:.2f} m')

    #to be called per time frame
    def convert_moments(self, ncoord_xb, ntheta:int, time:int):
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



    #to be called per time frame
    def get_moments_bdry(self):
         if self.im.is_valid(self.rmcb):
            print('RMCB already filled')
            return
         ntimes = len(self.im.get_ncvar('TIME'))
         nmoments = self.im.get_var_names_re('RMCB?', nd1=1, nd2=2, length=True)
         self.nmoments = nmoments
         self.rmcb = np.zeros(shape=(ntimes, nmoments))
         self.rmsb = np.zeros(shape=(ntimes, nmoments))
         self.ymcb = np.zeros(shape=(ntimes, nmoments))
         self.ymsb = np.zeros(shape=(ntimes, nmoments))

         cm2m = 0.01
         for k in range(nmoments):
            kn = k
            #RMCBn(TIME); CM; nTH ASYM R BOUND MOMENT
            self.rmcb[:, k] = self.im.get_ncvar(f'RMCB{kn}', scale=cm2m)
            self.ymcb[:, k] = self.im.get_ncvar(f'YMCB{kn}', scale=cm2m)
            if k>0:
               self.rmsb[:, k] = self.im.get_ncvar(f'RMSB{kn}', scale=cm2m)
               self.ymsb[:, k] = self.im.get_ncvar(f'YMSB{kn}', scale=cm2m)

    #to be called per time frame
    def convert_moments_bdry(self, ntheta:int, time:int):
         if not self.im.is_valid(self.rmcb, self.ymcb, ind=time):
            self.get_moments_bdry()
         if not self.im.is_valid(self.rmcb, self.ymcb, ind=time):
            return [], []

         if ntheta<1 or time<0:
            raise Exception(f'convert_moments_bdry: ntheta {ntheta} time {time}')
         rsurf = np.zeros(ntheta)
         ysurf = np.zeros(ntheta)
         if self.dthetab is None:
            raise Exception('Please define dthetab')

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
         #psirz is 1D, to be reshaped
         if nr==0 or nz==0 or nr<time or nz<time or len(psirz)==0:
            return [],[]
         rbdy = []
         zbdy = []
         plt.figure()
         #PBPMHD , PLCFSMHD
         psib = self.im.get_ncvar('PLFLXA', scale=2.*self.im.pi, ind=time)
         if np.isscalar(psib):
            psib = 0.99*psib

         psi2d = np.reshape(psirz, (nz, nr)) #.T
         #psi2d = np.reshape(psirz, (nr, nz))

         raxis = self.im.get_ncvar('RAXIS', scale=0.01, ind=time)
         yaxis = self.im.get_ncvar('YAXIS', scale=0.01, ind=time)

         eps = 0.01
         try:
            contour = plt.contour(rgrid[time], zgrid[time], psi2d, levels=[psib])
            for cont in contour.allsegs[0]:
               #check if axis is included in this bounded contour
               p = path.Path(cont)
               has_axis = p.contains_point((raxis,yaxis))
               closed = abs(cont[:,0][0]-cont[:,0][-1])<eps and abs(cont[:,1][0]-cont[:,1][-1])<eps
               if has_axis and closed:
                  rbdy, zbdy = cont[:,0], cont[:,1]
                  rbdy[0] = rbdy[-1]
                  zbdy[0] = zbdy[-1]
                  break
         except:
            if time ==0:
               print(f'Error tracing boundary for PSI= {psib}')
            pass

         #to debug, change to True
         if False and (time<2 or time%100==0):
            cts = plt.contour(rgrid[time], zgrid[time], psi2d, levels=30)
            plt.clabel(cts, inline=1, fontsize=10)
            plt.plot(rbdy, zbdy, 'r-')
            limitr = self.im.get_ncvar('RLIM', scale=0.01, ind=time)
            limitz = self.im.get_ncvar('YLIM', scale=0.01, ind=time)
            plt.plot(limitr, limitz)
            plt.show()

         plt.close('all')
         return rbdy, zbdy


    def boundary_minor(self, rsurf, time):
         minor = -1
         try:
            minor = (np.max(rsurf) - np.min(rsurf))/2.0
         except:
            if time==0:
               self.im.logger.info('minor_axis: plasma boundary invalid')
         return minor


    def create_ids(self):
        self.im.logger.info('\nequilibrium')
        eqm = self.im.ids.equilibrium
        cm2m = 0.01
        pi = constants.pi  #scipy
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')
        if not self.im.is_valid(timesc, timep, name='times'):
            return
        nprof = coord_x.shape[1]

        currents = self.im.get_ncvar('PCURC') #(TIME)
        raxis = self.im.get_ncvar('RAXIS', cm2m) #RAXISEQ
        yaxis = self.im.get_ncvar('YAXIS', cm2m) #YAXISEQ

        psi_axis = self.im.get_ncvar('PSI0_TR', 2*pi) # Wb/rad
        #psi_bdry =
        q95 = self.im.get_ncvar('Q95')
        qaxis  = self.im.get_ncvar('Q0')
        li3 = self.im.get_ncvar('LI_3')
        bpeq = self.im.get_ncvar('BETAT') #poloidal beta
        btor = self.im.get_ncvar('BTEQ') #toroidal beta
        #Zone volume
        gvol = self.im.get_ncvar('DVOL', 1.0e-6) #(TIME3,X) cm^3
        #plasma volume
        pvol = self.im.get_ncvar('PVOL', 1.0e-6) #(TIME) cm^3
        if not self.im.is_valid(gvol, pvol, name='vol'):
            return

        #PMHD(TIME3, X). Pressure from equilibrium solver
        pmhd = self.im.get_ncvar('PMHD')
        emhd = self.im.np_multiply_sum(pmhd, gvol, axis=1, scale=1.5)

        self.logger.info(f'ntimesc {ntimes} nprof {nprof} raxis {raxis.shape}')

        #TOTAL BETA TOROIDAL.
        # 2 mu0 integral(p dV) / V / B0^2   ??
        bttot = self.im.get_ncvar('BTTOT') #(TIME3, X)
        beta_normal = [] #replaced by beta_tor_norm
        try:
          for it in range(ntimep):
              beta_normal.append( np.multiply(bttot, gvol) / pvol[:,None] )
        except:
          self.logger.info(f'Error- beta_normal: bttot {bttot.shape} gvol {gvol.shape}')
        else:
          beta_normal = np.array(beta_normal)

        beta_tor_norm = []
        bz = self.im.get_ncvar('BZ') #(TIME) T
        # PAREA(TIME) - plasma CS area[cm^2]
        area = self.im.get_ncvar('PAREA')*1.e-4
        try:
            # 100 * beta_tor * a[m] * B0 [T] / ip [MA]   ???
            a = np.sqrt(area/pi)
            beta_tor_norm = 1.0e8*btor*a*bz/currents
        except:
            pass

        #Area of the LCFS poloidal cross section
        #ZONE CROSS SECTIONAL AREA
        darea = self.im.get_ncvar('DAREA', 1.0e-4) #(TIME3,X) cm^2

        #Surface area of the toroidal flux surface
        #FLUX SURFACE AREA
        surface = self.im.get_ncvar('SURF', 1.0e-4) #(TIME3,X) cm^2

        #rminor =

        #profiles (TIME3, XB)
        elong = self.im.get_ncvar('ELONG')
        triang = self.im.get_ncvar('TRIANG')
        triang_lo = self.im.get_ncvar('TRIANGL')
        triang_up = self.im.get_ncvar('TRIANGU')
        square_up = self.im.get_ncvar('SQUARE_UO')
        square_lo = self.im.get_ncvar('SQUARE_LO')

        #if scalar tiem is different from profile time
        elongb = np.interp(timesc, timep, elong[:,-1])
        triangb = np.interp(timesc, timep, triang[:,-1])
        triang_lob = np.interp(timesc, timep, triang_lo[:,-1])
        triang_upb = np.interp(timesc, timep, triang_up[:,-1])

        bdry_sepx_psi = self.im.get_ncvar('PBPMHD', 2.0*pi)
        psi = self.im.get_ncvar('PLFLX2PI')
        if not self.im.is_valid:
            psi = self.im.get_ncvar('PLFLX', 2.*pi) #(TIME3, XB)
        phi = self.im.get_ncvar('TRFLX', 2.*pi)
        qpsi = self.im.get_ncvar('Q')
        pres = self.im.get_ncvar('PPLAS') #'PMHD')
        fbz = self.im.get_ncvar('GFUN')

        nprof1d = 0
        try:
            nprof1d = psi.shape[1]
        except:
            print(' Error getting equilibrium psi')
            return

        bzxR = self.im.get_ncvar('BZXR', cm2m) #(TIME) TESLA*CM. vac field
        b_field_tor = self.im.operate_arrays(bzxR, '/', raxis)

        #RMJMP flux surface centers
        geom_axisr = self.im.get_ncvar('RMJMP', cm2m)
        geom_axisz = []

        #f
        fprof = []
        try:
            self.logger.info(f'nprof1d {nprof1d} bzxR {bzxR.shape} fbz {fbz.shape}')
            fprof = fbz * bzxR[:,None] # BZXR*FBZ
        except:
            pass

        #Normalised tor flux coord
        rho_tor_norm = coord_x #np.multiply(coord_x, coord_x)

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
        ntheta = 201
        self.dtheta = 2*pi/ntheta

        nthetab = 128 #TODO
        self.dthetab = 2*pi/nthetab

        if self.inv_grid_theta is None:
           self.inv_grid_theta = [-self.im.pi + ith*self.dtheta for ith in range(ntheta)]

        eqm = self.im.fill_properties(eqm)
        self.im.fill_code_properties(eqm)
        eqm.ids_properties.homogeneous_time = 1
        eqm.time.resize(ntimes)
        eqm.time = timesc
        eqm.time_slice.resize(ntimes)

        self.logger.info(f'equilibrium #time slices {ntimes}')
        for it in range(ntimes):
            print_progress(ntimes-1, it)
            if it==0 or it==ntimes-1:
                self.logger.info(f'time slice: {it}')
            #global
            if len(currents) > it :
                eqm.time_slice[it].global_quantities.ip = currents[it].astype(float)
            if len(raxis) > it:
                eqm.time_slice[it].global_quantities.magnetic_axis.r = raxis[it].astype(float)
                eqm.time_slice[it].global_quantities.magnetic_axis.z = yaxis[it].astype(float)
            #total toroidal magnetic field at the magnetic axis
            if len(b_field_tor) > it:
                eqm.time_slice[it].global_quantities.magnetic_axis.b_field_tor = b_field_tor[it].astype(float)

            if len(psi_axis) > it:
                eqm.time_slice[it].global_quantities.psi_axis = psi_axis[it].astype(float)
            #last  or ALCFSMHD ??
            if len(psi) > it:
                eqm.time_slice[it].global_quantities.psi_boundary = psi[it,-1].astype(float)

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
                eqm.time_slice[it].global_quantities.beta_normal = beta_tor_norm[it].astype(float)
            if len(pvol) > it:
                eqm.time_slice[it].global_quantities.volume = pvol[it].astype(float)
            if len(darea) > it:
                eqm.time_slice[it].global_quantities.area = darea[it,-1].astype(float)
            #last ??
            if len(surface) > it:
                eqm.time_slice[it].global_quantities.surface = surface[it,-1].astype(float)

            if len(emhd) >it:
                eqm.time_slice[it].global_quantities.energy_mhd = emhd[it].astype(float)

            #2D psi
            psirz = self.im.get_ncvar('PSIRZ', 2.*pi, ind=it)

            rbound, ybound = self.convert_moments_bdry(nthetab, it)
            #rbound, ybound = self.trace_boundary(psirz, rgrid, zgrid, it, nr, nz)
            self.im.fill_data(eqm.time_slice[it].boundary.outline.r, rbound, name='rbound')
            self.im.fill_data(eqm.time_slice[it].boundary.outline.z, ybound, name='ybound')

            minor = self.boundary_minor(rbound, it)
            if minor > 0:
               eqm.time_slice[it].boundary.minor_radius = minor

            # what is being written here is not the geometric axis
            geom_rz = []    #TODO
            if len(geom_rz)>1 and geom_rz[0]>0:
               eqm.time_slice[it].boundary.geometric_axis.r = geom_rz[0]
               eqm.time_slice[it].boundary.geometric_axis.z = geom_rz[1]

            if len(bdry_sepx_psi) >it:
                eqm.time_slice[it].boundary_separatrix.psi = bdry_sepx_psi[it]

            if len(triangb) > it:
                eqm.time_slice[it].boundary.triangularity = triangb[it].astype(float)
            if len(triang_upb) > it:
                eqm.time_slice[it].boundary.triangularity_upper = triang_upb[it].astype(float)
            if len(triang_lob) > it:
                eqm.time_slice[it].boundary.triangularity_lower = triang_lob[it].astype(float)

            if len(elongb) > it:
                eqm.time_slice[it].boundary.elongation = elongb[it]

            #profiles 1d
            iprof = eqm.time_slice[it].profiles_1d
            self.im.fill_data(iprof.psi, psi, ind=it, name='psi')
            self.im.fill_data(iprof.phi, phi, ind=it, name='phi')
            self.im.fill_data(iprof.q, qpsi, ind=it, name='qpsi')
            self.im.fill_data(iprof.f, fprof, ind=it, name='fprof')
            self.im.fill_data(iprof.pressure, pres, ind=it, name='pres')

            self.im.fill_data(iprof.rho_tor_norm, rho_tor_norm, ind=it, name='rho_tor_norm')
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

                if it==ntimes-1:
                    self.logger.info(f't {it}: wrote 2d profile : rgrid {nr}, zgrid {nz} psi {prof2.psi.shape}')

            bfield = self.im.get_ncvar('B_FIELD', ind=it)
            if self.im.is_valid(bfield, name='bfield'):
                bfield = bfield.reshape((3,rsh[1],zsh[1]))
                prof2.b_field_r.resize((nr, nz))
                prof2.b_field_r = bfield[0]
                prof2.b_field_tor.resize((nr, nz))
                prof2.b_field_tor = bfield[1]
                prof2.b_field_z.resize((nr, nz))
                prof2.b_field_z = bfield[2]
                #self.im.(prof2.b_field_r, bfield[0])
                #self.im.fill_data(prof2.b_field_tor, bfield[1])
                #self.im.fill_data(prof2.b_field_z, bfield[2])
                if it==ntimes-1:
                    self.logger.info(f't {it}: wrote B_FIELD {bfield.shape}')

            #inverse equilibrium
            coord = eqm.time_slice[it].coordinate_system

            #inverse coordinates
            rsurf, ysurf = [],[] #self.convert_moments(ncoord_xb, ntheta, it)

            #fill data
            coord.grid_type.name = "inverse"
            coord.grid_type.index = 2
            self.im.fill_data(coord.grid.dim1, coord_xb, ind=it, name='grid_r')
            self.im.fill_data(coord.grid.dim2, self.inv_grid_theta, name='grid_theta')
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

