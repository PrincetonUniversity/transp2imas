#!/usr/bin/python
"""Source-only mock IDS structures for theory and neutral tests."""
import types
import numpy as np

class FakeNode:
    """Small schema-shaped object to exercise actual writer assignments."""

    def __init__(self):
        self._items = []

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        node = FakeNode()
        setattr(self, name, node)
        return node

    def resize(self, count, **kwargs):
        length = int(count)
        self._items = (self._items[:length] +
                       [FakeNode() for _ in range(max(0, length-len(self._items)))])

    def __len__(self):
        return len(self._items)

    def __getitem__(self, index):
        return self._items[index]


class FakeConverter:
    def __init__(self, nc):
        self.nc = nc
        self.vars_ = list(nc.variables)
        self.cm2m = .01
        self.aconv = 1e-4
        self.vconv = 1e-6
        self.vol_inv = 1e6
        self.area_inv = 1e4
        self.ip_sign = 1.0
        self.toroidal_rotation_sign = -1.0
        self.logger = types.SimpleNamespace(info=lambda *a, **k: None,
                                            warning=lambda *a, **k: None)
        self.written = None

    def check_profile_shape(self, name, array, nt, nx):
        """As ``nc2imas.check_profile_shape``, with the broken mapping raised (the real one follows the policy)."""
        values = np.asarray(array)
        if values.ndim == 2 and values.shape[0] == nt and values.shape[1] != nx:
            raise ValueError(f'radial size mismatch {values.shape[1]} != {nx} for {name}')
        return array

    def get_ncvar(self, name, scale=1., **kw):
        if name not in self.nc.variables:
            return np.array([])
        value = np.asarray(self.nc.variables[name].data, dtype=np.float64)
        if value.ndim == 0:
            value = value.reshape(1)
        return value * scale

    def is_valid(self, *args, **kwargs):
        return all(isinstance(x, np.ndarray) and x.size>0 for x in args)

    def check_var_name(self, name):
        return name in self.nc.variables

    def new_ids(self, name):
        return FakeNode()

    def fill_properties(self, ids):
        return ids

    def fill_code_properties(self, ids):
        ids.code.parameters = '<transp2imas_parameters/>'

    def fill_external_model_inputs(self, *args):
        return []

    def fill_vac_tor_field(self, ids):
        return None

    def fill_data(self, node, values, ind=-1, chk=-1, **kwargs):
        if not self.is_valid(values):
            return False
        data=np.asarray(values, dtype=np.float64)
        if ind >= 0:
            data=data[ind]
        if chk > 0 and len(data) != chk:
            raise AssertionError(f'writer assigned wrong radial count {len(data)} vs {chk}')
        node.values=np.asarray(data).copy()
        return True

    def cumulative_zone_quantity(self, values):
        return np.cumsum(values,axis=1)

    def interp_radial_profiles(self, values, xfrom, xto, axis_value=None):
        return np.array([np.interp(x, xf, y) for xf, y, x in zip(xfrom,values,xto)])

    def dd_at_least(self,major,minor=0,patch=0):
        return (4,1,1)>=(major,minor,patch)

    def set_species_name(self, node, name):
        node.name=name

    def set_nuclear_charge(self,node,charge):
        node.z_n=int(round(charge))

    def putdb_imas(self, ids):
        self.written=ids



def load_writer(name):
    import sys
    from unittest import mock
    with mock.patch.dict(sys.modules, {'imas': types.ModuleType('imas'),
                                      'netCDF4': types.ModuleType('netCDF4')}):
        sys.modules['netCDF4'].Dataset=object
        if name=='core_transport':
            from imas_core_transport import imas_core_transport
            return imas_core_transport
        if name=='core_profiles':
            from imas_core_profiles import imas_core_profiles
            return imas_core_profiles
        raise ValueError(name)

def write_core_transport(nc):
    writer_cls=load_writer('core_transport')
    im=FakeConverter(nc)
    writer=writer_cls(im)             # IdsWriter.__init__ attaches the context views over the fake
    writer.time,writer.x,writer.xb=(im.get_ncvar(k) for k in ('TIME3','X','XB'))
    writer.dvol=im.get_ncvar('DVOL', im.vconv)
    writer.surf=im.get_ncvar('SURF', im.aconv)
    writer.nt,writer.nx=len(writer.time),writer.xb.shape[1]
    writer.run()
    return im.written
