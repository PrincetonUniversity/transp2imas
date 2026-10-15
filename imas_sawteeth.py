#!/usr/bin/python
"""Write validated TRANSP sawtooth history to the IMAS sawteeth IDS."""

from nc2imas import nc2imas
from sawteeth_mapping import build_sawteeth_plan, populate_sawteeth_ids


class imas_sawteeth:
    def __init__(self, ncimas: nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.logger.info('\nsawteeth')
        plan = build_sawteeth_plan(self.im.nc)
        if plan is None:
            self.logger.info('sawteeth: no new TLASTSAW crash events; IDS skipped')
            return

        ids = self.im.new_ids('sawteeth')
        self.im.fill_properties(ids)
        self.im.fill_code_properties(ids)
        populate_sawteeth_ids(ids, plan)
        comment = ids.ids_properties.comment
        comment_text = comment.value if hasattr(comment, 'value') else comment
        ids.ids_properties.comment = str(comment_text) + (
            f'; TLASTSAW source-reported new crash updates={plan.events}; '
            'TIME-native samples; initial crash timestamp is left-censored; '
            'trigger mechanism and mixing/inversion radii not inferred'
        )
        if plan.source_boundary_samples:
            self.logger.warning('sawteeth: nonzero XSAWBOUND is not mapped to a '
                                'mixing/inversion radius without semantic validation')
        self.im.putdb_imas(ids)
        self.logger.info('sawteeth: saved %d observed crash updates on %d TIME samples',
                         plan.events, plan.time.size)


if __name__ == '__main__':
    print('Please use: python convert.py')
