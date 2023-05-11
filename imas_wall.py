from nc2imas import *



class imas_wall:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.im.print_debug()
        self.create_ids()

    def valid(self, *vars_) -> bool:
        return self.im.checkValid(vars_)


    def split_limiters(self, limitr, limitz):
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
                lim_lr = np.flip(limitr[mid:])
                lim_lz = np.flip(limitz[mid:])
                self.logger.info(f'counter clockwise limiters(start from right), split into units @ {mid}')
            elif l2r==-1 and zz==-1 and limitz[i]>0:
                mid = i
                lim_lr = limitr[:mid]
                lim_lz = limitz[:mid]
                lim_ur = limitr[mid:]
                lim_uz = lmitz[mid:]
                self.logger.info(f'clockwise limiter(start from right), split into units @ {mid}')
            if mid>0:
                break
        if mid<1 or len(lim_ur)<1 or len(lim_lr)<1:
            raise AssertionError(f'Limiters array couldnot be split')
        return lim_ur, lim_uz, lim_lr, lim_lz

    def create_ids(self):
        wall = self.im.ids.wall
        cm2m = 0.01

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)

        self.valid(timesc)

        wall = self.im.fill_properties(wall)
        self.im.fill_code_properties(wall)
        wall.ids_properties.homogeneous_time = 1
        wall.time.resize(ntimes)
        wall.time = timesc

        limitr = self.im.get_ncvar('RLIM', cm2m)
        limitz = self.im.get_ncvar('YLIM', cm2m)
        ilim = self.im.get_ncvar('ILIM')

        self.logger.info(f'limitr {len(limitr)} limitz {len(limitz)} ilim {len(ilim)}')
        self.logger.debug(f'limitr {limitr}')
        self.logger.debug(f'limitz {limitz}')
        self.logger.debug(f'ilim {ilim}')
        self.valid(limitr, limitz, ilim)

        tselect = 0 #first time frame
        limitr = limitr[tselect]
        limitz = limitz[tselect]
        ilim = ilim[tselect]

        #split limiter into upper and lower parts
        lim_ur, lim_uz, lim_lr, lim_lz = self.split_limiters(limitr, limitz)
        nwall_up = len(lim_ur)
        nwall_low = len(lim_lr)
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
    a = get_args(sys.argv)

    #sh, rn, mach = 142202, 'U01', 'NSTX' #  142202U01.CDF
    sh, rn, mach = 147634, 'G14', 'D3D' # 147634G14.CDF

    ncf = str(sh)+rn+'.CDF'
    nml = str(sh)+rn+'TR.DAT'
    nci = nc2imas(db=mach, shot=sh, backend='MDS', run=3, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile=ncf, nml_file=nml, log='debug')
    imas_wall(nci)