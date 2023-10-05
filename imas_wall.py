from nc2imas import *



class imas_wall:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def valid(self, *vars_) -> bool:
        return self.im.is_valid(vars_)

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

        self.logger.info(f'limitr {len(limitr)} limitz {len(limitz)} ')
        if not self.im.is_valid(limitr, limitz ):
            self.logger.debug('Error writing Wall')
            return

        tselect = 0 #first time frame
        limitr = limitr[tselect]
        limitz = limitz[tselect]
        nlim = len(limitr)

        wall.description_2d.resize(1)
        wall.description_2d[0].limiter.unit.resize(1)

        wall.description_2d[0].limiter.type.name = 'first_wall'
        wall.description_2d[0].limiter.type.index = 0
        wall.description_2d[0].limiter.type.description = 'first wall or limiters'

        # we use only 1 units, because TRANSP stores coordinates of limiter as a closed contour.
        # TRANSP has duplicate first and last point, we do not write the last point in IDS
        wall.description_2d[0].limiter.unit[0].closed = 1
        wall.description_2d[0].limiter.unit[0].outline.r.resize(nlim-1)
        wall.description_2d[0].limiter.unit[0].outline.z.resize(nlim-1)
        wall.description_2d[0].limiter.unit[0].outline.r = limitr[0:nlim-1]
        wall.description_2d[0].limiter.unit[0].outline.z = limitz[0:nlim-1]
 
        self.im.putdb_imas(wall)
        self.logger.info('Done wall')


if __name__ == '__main__':
    print("Please use: python convert.py ")
