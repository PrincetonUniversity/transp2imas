from nc2imas import *



class imas_wall:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def valid(self, *vars_) -> bool:
        return self.im.is_valid(vars_)

    #TODO add selection from unit
    def remove_duplicates(self, limitr, limitz):
      nr = len(limitr)
      limr, limz = [], []
      if nr<2:
         self.im.logger.info(f'insufficient limiter segments: # {nr}')
         return limr, limz
      dinds = []
      for i in range(nr):
         ii = nr-1-i
         xv = limitr[ii]
         yv = limitz[ii]
         for j in range(ii):
            if (xv==limitr[j] and yv==limitz[j]):
               dinds.append(ii)

      for i in range(nr):
         if i in dinds:
            self.im.logger.info(f'removing duplicate: {i} {limitr[i]} {limitz[i]}')
            continue
         limr.append(limitr[i])
         limz.append(limitz[i])
      return limr, limz

    def split_limiters(self, limitr, limitz):
        limitr, limitz = self.remove_duplicates(limitr, limitz)
        nr = len(limitr)
        lim_ur, lim_uz, lim_lr, lim_lz = [], [], [], []
        if nr<2:
            self.im.logger.info(f'insufficient limiter segments nr {nr}')
            return lim_ur, lim_uz, lim_lr, lim_lz

        xup = []
        xlow = []
        xz1 = limitz[0]
        pz = 1 if xz1>0 else -1
        if xz1>0:
          xup.append(0)
        else:
          xlow.append(0)

        for i in range(1, nr):
            xr = limitr[i]
            xz = limitz[i]
            if pz>0 and xz<=0:
               xup.append(i-1)
               xlow.append(i)
            if pz<0 and xz>0:
               xlow.append(i-1)
               xup.append(i)
            pz = 1 if xz>0 else -1

        #make upper and lower segments seperately
        for i in range(0,len(xup),2):
            n = xup[i]
            m = xup[i+1]+1 if len(xup)>i+1 else len(limitr)
            lim_ur = limitr[n:m] + lim_ur
            lim_uz = limitz[n:m] + lim_uz
        for i in range(0,len(xlow),2):
            n = xlow[i]
            m = xlow[i+1]+1 if len(xlow)>i+1 else len(limitr)
            lim_lr = lim_lr + limitr[n:m]
            lim_lz = lim_lz + limitz[n:m]
        lim_lr = np.flip(lim_lr)
        lim_lz = np.flip(lim_lz)
        return np.array(lim_ur), np.array(lim_uz), lim_lr, lim_lz


    def split_limiters_(self, limitr, limitz):
        limitr, limitz = self.remove_duplicates(limitr, limitz)
        #split limiter into upper and lower parts
        #left to right, and clockwise.
        #2nd and 3rd elements are either in increasing or decreasing order
        l2r = 1 if limitr[2]>limitr[1] else -1
        zz = 1 if limitz[2]>0 else -1

        #TRANSP expects limiter in init 1 upper: clockwise from mid-left to mid-right
        # and unit 2 lower: clockwise mid-right to mid-left.
        mid = -1
        for i in range(1, len(limitr)):
            if l2r==1 and zz==1 and limitz[i]<0:
                mid = i #excluding
                #exclude limitr[mid] and beyond
                lim_ur = limitr[:mid]
                lim_uz = limitz[:mid]
                lim_lr = limitr[mid:]
                lim_lz = limitz[mid:]
                self.logger.info(f'clockwise limiters(start from left), split into units @ {mid}')
            elif l2r==1 and zz==-1 and limitz[i]>0:
                mid = i
                lim_lr = np.flip(limitr[:mid],0)
                lim_lz = np.flip(limitz[:mid],0)
                lim_ur = np.flip(limitr[mid:],0)
                lim_uz = np.flip(limitz[mid:],0)
                self.logger.info(f'counter clockwise limiters(start from left), split into units @ {mid}')
            elif l2r==-1 and zz==1 and limitz[i]<0:
                mid = i #excluding
                #exclude limitr[mid] and beyond
                lim_ur = np.flip(limitr[:mid])
                lim_uz = np.flip(limitz[:mid])
                lim_lr = limitr[mid:]
                lim_lz = limitz[mid:]
                self.logger.info(f'counter clockwise limiters(start from right), split into units @ {mid}')
            elif l2r==-1 and zz==-1 and limitz[i]>0:
                mid = i
                lim_lr = limitr[:mid]
                lim_lz = limitz[:mid]
                lim_ur = limitr[mid:]
                lim_uz = limitz[mid:]
                self.logger.info(f'clockwise limiter(start from right), split into units @ {mid}')
            if mid>0:
                break
        if mid<1 or len(lim_ur)<1 or len(lim_lr)<1:
            raise AssertionError(f'Limiters array couldnot be split')
        return np.array(lim_ur), np.array(lim_uz), np.array(lim_lr), np.array(lim_lz)

    def create_ids(self):
        self.im.logger.info('\nwall')
        wall = self.im.ids.wall
        cm2m = 0.01

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)

        self.im.checkErrExit(timesc)

        wall = self.im.fill_properties(wall)
        self.im.fill_code_properties(wall)
        wall.ids_properties.homogeneous_time = 1
        wall.time.resize(ntimes)
        wall.time = timesc

        limitr = self.im.get_ncvar('RLIM', cm2m)
        limitz = self.im.get_ncvar('YLIM', cm2m)
        ilim = self.im.get_ncvar('ILIM')

        self.logger.info(f'limitr {len(limitr)} limitz {len(limitz)} ilim {len(ilim)}')
        if not self.im.is_valid(limitr, limitz, ilim):
            self.logger.debug('Error writing Wall')
            return

        tselect = 0 #first time frame
        limitr = limitr[tselect]
        limitz = limitz[tselect]
        ilim = ilim[tselect]

        #split limiter into upper and lower parts
        lim_ur, lim_uz, lim_lr, lim_lz = self.split_limiters(limitr, limitz)
        nwall_up = len(lim_ur)
        nwall_low = len(lim_lr)

        if nwall_up<2 or nwall_low<2:
            self.logger.info(f'too few limiter segments: upper: {nwall_up} lower: {nwall_low}')
            return

        self.logger.info(f'times {ntimes} nlimiters {nwall_up} {nwall_low}')

        wall.description_2d.resize(1)
        wall.description_2d[0].limiter.unit.resize(2)

        wall.description_2d[0].limiter.type.name = 'first_wall'
        wall.description_2d[0].limiter.type.index = 0
        wall.description_2d[0].limiter.type.description = 'first wall or limiters'

        #upper unit 0, from mid-left to mid-right clockwise
        #part of contour is not closed
        wall.description_2d[0].limiter.unit[0].closed = 0
        wall.description_2d[0].limiter.unit[0].outline.r.resize(nwall_up)
        wall.description_2d[0].limiter.unit[0].outline.z.resize(nwall_up)
        wall.description_2d[0].limiter.unit[0].outline.r = lim_ur
        wall.description_2d[0].limiter.unit[0].outline.z = lim_uz

        #lower unit 1, from mid-right to mid-left clockwise
        #part of contour is not closed
        wall.description_2d[0].limiter.unit[1].closed = 0
        wall.description_2d[0].limiter.unit[1].outline.r.resize(nwall_low)
        wall.description_2d[0].limiter.unit[1].outline.z.resize(nwall_low)
        wall.description_2d[0].limiter.unit[1].outline.r = lim_lr
        wall.description_2d[0].limiter.unit[1].outline.z = lim_lz

        self.im.putdb_imas(wall)
        self.logger.info('Done wall')


if __name__ == '__main__':
    print("Please use: python convert.py ")