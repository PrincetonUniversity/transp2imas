#!/usr/bin/python
"""Verify PT_SOLVER runtime numerical diagnostics in an actual TRANSP CDF.
Usage: PYTHONPATH=. python examples/verify_pt_runtime_actual_cdf.py examples/data/iter/38530P94.CDF
Requires scipy and numpy. No IMAS installation required.
"""
import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pt_solver_diagnostics import channel_coordinates, diagnostic_xml, decode_channel, KNOWN_CHANNELS


def check(path):
    with netcdf_file(path, mmap=False) as nc:
        root = diagnostic_xml(nc)
        for channel in root.findall('channel'):
            name = channel.get('name')
            expected = np.asarray(nc.variables[name][:], dtype='<f8')
            actual = decode_channel(channel)
            np.testing.assert_array_equal(actual, expected)
            assert channel.get('quantity') == KNOWN_CHANNELS[name][0]
            # the channel's own native coordinate, not the TIME3 grid of the IDS
            for dim, values in channel_coordinates(root, channel).items():
                np.testing.assert_array_equal(values, np.asarray(nc.variables[dim][:], dtype='<f8'))
            assert channel.get('coordinates') == ','.join(channel.get('dimensions').split(',')), name
            print(f'{name}: {actual.shape}; {channel.get("quantity")}; {channel.get("units")}')
        assert not root.findall('channel') or int(root.get('count')) == len(root.findall('channel'))
        assert not any(ch.get('name') in {'ISOITER','RZITER','FMCK_WNC','PT_TIME'} for ch in root.findall('channel'))
        print('PASS: verified runtime channels:', root.get('count'))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('cdf')
    check(p.parse_args().cdf)
