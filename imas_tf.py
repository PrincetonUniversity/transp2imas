from nc2imas import *

class imas_tf:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\ntf')
        tf = self.im.ids.get('tf')
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        #BZXR(TIME) TESLA*CM ;VACUUM FIELD
        self.im.fill_data_nc(tf.b_field_tor_vacuum_r.data,'BZXR', mult=0.01)
        write = True

        if write:
            tf = self.im.fill_properties(tf)
            self.im.fill_code_properties(tf)
            tf.ids_properties.homogeneous_time = 1
            tf.time.resize(ntimes)
            tf.time = timesc
            self.im.putdb_imas(tf)
            self.logger.info('Done TF')

if __name__ == '__main__':
    print("Please use: python convert.py ")
