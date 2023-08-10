
from nc2imas import *

class imas_equilibrium:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()


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
        #
        #PBPMHD

        #boundary_separatrix
        #psi_sepx = self.im.get_ncvar('

        #TOTAL BETA TOROIDAL.
        # 2 mu0 integral(p dV) / V / B0^2   ??
        bttot = self.im.get_ncvar('BTTOT') #(TIME3, X)
        beta_normal = [] #replaced by beta_tor_norm
        try:
          for i in range(ntimep):
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
        rho_tor_norm = np.multiply(coord_x, coord_x)

        rgrid = self.im.get_ncvar('RGRID', cm2m)
        zgrid = self.im.get_ncvar('ZGRID', cm2m)

        rsh = rgrid.shape
        if len(rsh)>1 and rsh[1]>0 and rgrid[0,rsh[1]-1]<0.2:
            print('Check RGRID unit:', rgrid[0])
            #raise Exception('RGRID too small:expected in cm')
        self.logger.info(f'rgrid {rgrid.shape}')

        #RSURF - flux surface R locations
        #rsurf, ysurf

        eqm = self.im.fill_properties(eqm)
        self.im.fill_code_properties(eqm)
        eqm.ids_properties.homogeneous_time = 1
        eqm.time.resize(ntimes)
        eqm.time = timesc
        eqm.time_slice.resize(ntimes)

        self.logger.info(f'equilibrium #time slices {ntimes}')
        for i in range(ntimes):
            if i==0 or i==ntimes-1:
                self.logger.info(f'time slice: {i}')
            #global
            if len(currents) > i :
                eqm.time_slice[i].global_quantities.ip = currents[i].astype(float)
            if len(raxis) > i:
                eqm.time_slice[i].global_quantities.magnetic_axis.r = raxis[i].astype(float)
                eqm.time_slice[i].global_quantities.magnetic_axis.z = yaxis[i].astype(float)
            #total toroidal magnetic field at the magnetic axis
            if len(b_field_tor) > i:
                eqm.time_slice[i].global_quantities.magnetic_axis.b_field_tor = b_field_tor[i].astype(float)

            if len(psi_axis) > i:
                eqm.time_slice[i].global_quantities.psi_axis = psi_axis[i].astype(float)
            #last  or ALCFSMHD ??
            if len(psi) > i:
                eqm.time_slice[i].global_quantities.psi_boundary = psi[i,-1].astype(float)

            if len(q95) > i:
                eqm.time_slice[i].global_quantities.q_95 = q95[i].astype(float)
            if len(qaxis) >i:
                eqm.time_slice[i].global_quantities.q_axis = qaxis[i].astype(float)
            if len(li3) > i:
                eqm.time_slice[i].global_quantities.li_3 = li3[i].astype(float)
            if len(bpeq) >i:
                eqm.time_slice[i].global_quantities.beta_pol = bpeq[i].astype(float)
            if len(btor) >i:
                eqm.time_slice[i].global_quantities.beta_tor = btor[i].astype(float)
            if len(beta_tor_norm):
                eqm.time_slice[i].global_quantities.beta_normal = beta_tor_norm[i].astype(float)
            if len(pvol) > i:
                eqm.time_slice[i].global_quantities.volume = pvol[i].astype(float)
            if len(darea) > i:
                eqm.time_slice[i].global_quantities.area = darea[i,-1].astype(float)
            #last ??
            if len(surface) > i:
                eqm.time_slice[i].global_quantities.surface = surface[i,-1].astype(float)

            if len(emhd) >i:
                eqm.time_slice[i].global_quantities.energy_mhd = emhd[i].astype(float)

            #if len(rminor) > i:
            #    eqm.time_slice[i].boundary.minor_radius = rminor[i]
            #if len(outline_r) >i:
            #   eqm.time_slice[i].boundary.outline.r.resize(outline_r.shape[1])
            #   eqm.time_slice[i].boundary.outline.z.resize(outline_z.shape[1])
            #   eqm.time_slice[i].boundary.outline.r = outline_r[i]  #1D
            #   eqm.time_slice[i].boundary.outline.z = outline_z[i]
            #if len(psi_sepx) >i:
            #    eqm.time_slice[i].boundary_separatrix.psi = psi_sepx[i]

            if len(triangb) > i:
                eqm.time_slice[i].boundary.triangularity = triangb[i].astype(float)
            if len(triang_upb) > i:
                eqm.time_slice[i].boundary.triangularity_upper = triang_upb[i].astype(float)
            if len(triang_lob) > i:
                eqm.time_slice[i].boundary.triangularity_lower = triang_lob[i].astype(float)

            if len(elongb) > i:
                eqm.time_slice[i].boundary.elongation = elongb[i]

            #profiles 1d
            iprof = eqm.time_slice[i].profiles_1d
            self.im.fill_data(iprof.psi, psi, ind=i, name='psi')
            self.im.fill_data(iprof.phi, phi, ind=i, name='phi')
            self.im.fill_data(iprof.q, qpsi, ind=i, name='qpsi')
            self.im.fill_data(iprof.f, fprof, ind=i, name='fprof')
            self.im.fill_data(iprof.pressure, pres, ind=i, name='pres')

            self.im.fill_data(iprof.rho_tor_norm, rho_tor_norm, ind=i, name='rho_tor_norm')
            self.im.fill_data(iprof.elongation, elong, ind=i, name='elong')
            self.im.fill_data(iprof.squareness_upper_outer, square_up, ind=i, name='square_up')
            self.im.fill_data(iprof.squareness_lower_outer, square_lo, ind=i, name='square_lo')
            self.im.fill_data(iprof.triangularity_upper, triang_up, ind=i, name='triang_up')
            self.im.fill_data(iprof.triangularity_lower, triang_lo, ind=i, name='triang_lo')
            self.im.fill_data(iprof.geometric_axis.r, geom_axisr, ind=i, name='geom_axisr')
            self.im.fill_data(iprof.geometric_axis.z, geom_axisz, ind=i, name='geom_axisz')

            #self.im.fill_data(iprof.j_tor, jtor, ind=i, name='jtor')
            #self.im.fill_data(iprof.f_df_dpsi, fdp_dspi, ind=i, name='fdp_dspi')
            #self.im.fill_data(iprof.dpressure_dpsi, dpressure_dpsi, ind=i, name='dpressure_dpsi')

            eqm.time_slice[i].profiles_2d.resize(1)
            prof2 = eqm.time_slice[i].profiles_2d[0]
            self.im.fill_data(prof2.grid.dim1, rgrid, ind=i, name='rgrid')
            self.im.fill_data(prof2.grid.dim2, zgrid, ind=i, name='zgrid')
            #prof2.grid_type = 1 #Error

            #TODO all 2d profiles

            psirz = self.im.get_ncvar('PSIRZ', 2.*pi, ind=i)
            if self.im.is_valid(rgrid, zgrid, ind=i, name='grid') and self.im.is_valid(psirz, name='prirz'):
                xv, yv = np.meshgrid(rgrid[i], zgrid[i])
                prof2.r.resize(xv.shape)
                prof2.z.resize(yv.shape)
                prof2.r = xv
                prof2.z = yv

                rn = len(rgrid[i])
                zn = int(len(psirz)/rn)
                prof2.psi.resize((rn, zn))
                prof2.psi = np.reshape(psirz, (rn, zn)).T

                prof2.type = 0 #total
                #self.outdb.put_slice(eq, occur)
                if i==ntimes-1:
                    self.logger.info(f't {i}: wrote 2d profile : rgrid {rn}, zgrid {zn} psi {prof2.psi.shape}')
            elif i==ntimes-1:
                self.logger.info('2d profile : psi, rgrid, zgrid not available')


                #inverse equilibrium
                #Values of the major radius on the grid
            if False:
                eqm.time_slice[i].coordinate_system.grid.dim1.resize(nr)
                eqm.time_slice[i].coordinate_system.grid.dim2.resize(nz)
                eqm.time_slice[i].coordinate_system.grid.dim1 = rgrid[i]
                eqm.time_slice[i].coordinate_system.grid.dim2 = zgrid[i]
                eqm.time_slice[i].coordinate_system.r.resize(nr*nz)
                eqm.time_slice[i].coordinate_system.z.resize(nr*nz)
                #eqm.time_slice[i].coordinate_system.r = rsurf[i]
                #eqm.time_slice[i].coordinate_system.z = ysurf[i]

        self.logger.info('writing ids')
        self.im.putdb_imas(eqm)

if __name__ == '__main__':
    print("Please use: python convert.py ")

