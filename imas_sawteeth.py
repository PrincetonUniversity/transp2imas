
from nc2imas import *

class imas_sawteeth:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()


    def create_ids(self):
        self.im.logger.info('\nsawteeth')
        ids = self.im.ids.sawteeth
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        ids = self.im.fill_properties(ids)
        self.im.fill_code_properties(ids)
        ids.ids_properties.homogeneous_time = 1
        ids.time.resize(ntimes)
        ids.time = timesc

        nnn = 0



        if nnn > 0:
            self.im.putdb_imas(ids)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py")