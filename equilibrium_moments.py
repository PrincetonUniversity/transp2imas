
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

         for ith in range(ntheta-1):
            th = ith*self.dtheta
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


         #Fill
            #inverse equilibrium
            coord = eqm.time_slice[it].coordinate_system

            #TODO check dimension of fill data
            rsurf, ysurf = self.convert_moments(ncoord_xb, ntheta, it)

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

            #self.outdb.put_slice(eq, occur)