#!/usr/bin/python
"""Write validated TRANSP sawtooth history to the IMAS sawteeth IDS."""

from sawteeth_mapping import build_sawteeth_plan, populate_sawteeth_ids
from writer_base import IdsWriter


class imas_sawteeth(IdsWriter):
    def build(self):
        self.logger.info('\nsawteeth')
        plan = build_sawteeth_plan(self.cdf.nc)
        if plan is None:
            self.logger.info('sawteeth: no new TLASTSAW crash events; IDS skipped')
            return []

        ids = self.dd.new_ids('sawteeth')
        self.fill.fill_properties(ids)
        self.fill.fill_code_properties(ids)
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
        self.logger.info('sawteeth: built %d observed crash updates on %d TIME samples',
                         plan.events, plan.time.size)
        return [ids]


if __name__ == '__main__':
    print('Please use: python convert.py')
