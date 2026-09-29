#!/usr/bin/python
"""Populate the IMAS tf IDS from TRANSP toroidal-field information."""

from nc2imas import *

class imas_tf:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\ntf')
        tf = self.im.new_ids('tf')
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        # BZXR(TIME) TESLA*CM; vacuum toroidal field times major radius.
        # DD3.39 uses b_field_tor_vacuum_r; DD4.1.1 uses
        # b_field_phi_vacuum_r.  Late DD3 revisions also use the phi spelling.
        if self.im.dd_at_least(3, 42, 0):
            bvac = tf.b_field_phi_vacuum_r
        else:
            bvac = tf.b_field_tor_vacuum_r

        # IMAS positive toroidal field is counter-clockwise viewed from above.
        # Apply the same machine/COCOS sign policy used by equilibrium/core_profiles.
        self.im.fill_data_nc(bvac.data, 'BZXR', mult=0.01 * self.im.bt_sign)
        bvac.time = np.asarray(timesc, dtype=np.float64)
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
