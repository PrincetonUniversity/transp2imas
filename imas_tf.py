#!/usr/bin/python
"""Populate the IMAS tf IDS from TRANSP toroidal-field information."""

from nc2imas import *
from writer_base import IdsWriter


class imas_tf(IdsWriter):
    def build(self):
        self.logger.info('\ntf')
        tf = self.dd.new_ids('tf')
        #scalar time
        timesc = self.cdf.get_ncvar('TIME')
        self.fill.checkErrExit(timesc)

        # BZXR(TIME) TESLA*CM; vacuum toroidal field times major radius.
        # DD3.39 uses b_field_tor_vacuum_r; DD4.1.1 uses
        # b_field_phi_vacuum_r.  Late DD3 revisions also use the phi spelling.
        if self.dd.dd_at_least(3, 42, 0):
            bvac = tf.b_field_phi_vacuum_r
        else:
            bvac = tf.b_field_tor_vacuum_r

        # IMAS positive toroidal field is counter-clockwise viewed from above.
        # Apply the same machine/COCOS sign policy used by equilibrium/core_profiles.
        self.fill.fill_data_nc(bvac.data, 'BZXR', mult=0.01 * self.signs.bt_sign)
        bvac.time = np.asarray(timesc, dtype=np.float64)
        tf = self.fill.fill_properties(tf)
        self.fill.fill_code_properties(tf)
        tf.ids_properties.homogeneous_time = 1
        tf.time = timesc
        self.logger.info('Done TF')
        return [tf]

if __name__ == '__main__':
    print("Please use: python convert.py ")
