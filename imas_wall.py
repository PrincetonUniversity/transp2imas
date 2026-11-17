#!/usr/bin/python
"""Populate the IMAS wall IDS from available machine-boundary information."""

from nc2imas import *
from writer_base import IdsWriter


class imas_wall(IdsWriter):
    def valid(self, *vars_) -> bool:
        return self.fill.is_valid(vars_)

    def build(self):
        self.logger.info('\nwall')
        wall = self.dd.new_ids('wall')
        cm2m = 0.01

        #scalar time
        timesc = self.cdf.get_ncvar('TIME')
        ntimes = len(timesc)

        self.fill.checkErrExit(timesc)

        wall = self.fill.fill_properties(wall)
        self.fill.fill_code_properties(wall)
        wall.ids_properties.homogeneous_time = 1
        wall.time.resize(ntimes, refcheck=False)
        wall.time = timesc

        limitr = self.cdf.get_ncvar('RLIM', cm2m)
        limitz = self.cdf.get_ncvar('YLIM', cm2m)

        self.logger.info(f'limitr {len(limitr)} limitz {len(limitz)} ')
        if not self.fill.is_valid(limitr, limitz ):
            self.logger.debug('Error writing Wall')
            return []

        tselect = 0 #first time frame
        limitr = limitr[tselect]
        limitz = limitz[tselect]
        nlim = len(limitr)

        wall.description_2d.resize(1)
        wall.description_2d[0].limiter.unit.resize(1)

        wall.description_2d[0].limiter.type.name = 'first_wall'
        wall.description_2d[0].limiter.type.index = 0
        wall.description_2d[0].limiter.type.description = 'first wall or limiters'

        # We use one limiter unit because TRANSP stores the limiter as one closed contour.
        # DD3 represented closure with unit.closed=1 and omitted the duplicate endpoint.
        # DD4 removed the "closed" field: a closed contour is represented by repeating
        # the first point at the end of outline.r/z.  Support both layouts so the
        # converter can still be used with older Data Dictionaries.
        unit = wall.description_2d[0].limiter.unit[0]
        if self.dd.dd_before(4, 0, 0):
            # IMAS DD3.x
            unit.closed = 1
            rout = limitr[:-1] if nlim > 1 else limitr
            zout = limitz[:-1] if nlim > 1 else limitz
        else:
            # IMAS DD4.x. TRANSP normally already repeats the first point.  If a CDF
            # does not, explicitly close the contour here.
            if nlim > 1 and np.isclose(limitr[0], limitr[-1]) and np.isclose(limitz[0], limitz[-1]):
                rout = limitr
                zout = limitz
            else:
                rout = np.append(limitr, limitr[0])
                zout = np.append(limitz, limitz[0])

        unit.outline.r.resize(len(rout), refcheck=False)
        unit.outline.z.resize(len(zout), refcheck=False)
        unit.outline.r = rout
        unit.outline.z = zout
 
        self.logger.info('Done wall')
        return [wall]


if __name__ == '__main__':
    print("Please use: python convert.py ")
