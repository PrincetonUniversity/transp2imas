#!/usr/bin/python
"""What an IDS writer may use of the converter, split by role.

Writers used to take the whole ``nc2imas`` object and call 75 different things on it.  A
``ConversionContext`` groups that surface into small views, each with an explicit list of the names it
forwards, so a writer's dependencies are visible in its code (``self.cdf.get_ncvar``,
``self.signs.ip_sign``, ``self.dd.dd_at_least``) and can be replaced piece by piece in a test.

* ``cdf``       reading the TRANSP file: variables, name patterns, species, attributes
* ``namelist``  the TRANSP namelist and the input-file inventory
* ``dd``        the target Data Dictionary: version predicates and IDS creation
* ``signs``     COCOS sign conventions, the vacuum field and the betas defined against it
* ``grid``      time bases, unit factors and radial/time interpolation
* ``fill``      assigning values to IDS nodes, identities and the properties/code blocks

A view forwards to the converter at call time (``getattr`` on the owner), so it sees the same live
state and any object offering the listed names, a real converter or a test double, works as the owner.
Asking a view for a name it does not list is an error that names the view: nothing reaches the rest
of the converter by accident.  ``tests`` in ``examples/verify_conversion_context.py`` check every
listed name against ``nc2imas`` and every name a writer uses against these lists.
"""


class ContextView:
    """Forward the names in ``names`` to the owner; refuse every other name."""

    names = ()
    role = 'view'

    def __init__(self, owner):
        object.__setattr__(self, '_owner', owner)

    def __getattr__(self, name):
        if name.startswith('_') or name not in self.names:
            raise AttributeError(f'{self.role!r} view of the conversion context has no {name!r}; '
                                 f'it offers: {", ".join(self.names)}')
        return getattr(self._owner, name)

    def __setattr__(self, name, value):
        raise AttributeError(f'the {self.role!r} view is read-only; set {name!r} on the converter')

    def __dir__(self):
        return sorted(self.names)


class CdfView(ContextView):
    role = 'cdf'
    names = ('get_ncvar', 'get_ncvar_on_time', 'check_var_name', 'get_var_names', 'get_var_names_re',
             'get_var_names_re_', 'get_species_names', 'get_species_data', 'get_ncattr_text',
             'vars_', 'nc', 'ncfile', 'device', 'shot', 'iter_', 'variable_dimensions')


class NamelistView(ContextView):
    role = 'namelist'
    names = ('get_namelist_var', 'get_namelist_vars', 'nmlist', 'input_inventory', 'parse_number')


class DdView(ContextView):
    role = 'dd'
    names = ('dd_at_least', 'dd_before', 'new_ids', 'dd_major_version', 'dd_version')


class SignsView(ContextView):
    role = 'signs'
    names = ('ip_sign', 'bt_sign', 'psi_sign', 'q_sign', 'phi_sign', 'f_sign', 'jparallel_sign', 'vloop_sign',
             'toroidal_rotation_sign', 'poloidal_field_sign', 'jphi_sign', 'theta_raw_from_target',
             'vacuum_field_reference', 'vacuum_b0_magnitude', 'plasma_betas')


class GridView(ContextView):
    role = 'grid'
    names = ('timep', 'timesc', 'ntimep', 'ntimes', 'pi', 'cm2m', 'vol_inv', 'area_inv', 'vconv', 'aconv',
             'interp_radial_profiles', 'interp1d_time', 'cumulative_zone_quantity', 'radial_gradient')


class FillView(ContextView):
    role = 'fill'
    names = ('fill_data', 'fill_data_nc', 'fill_value', 'is_valid', 'operate_arrays', 'np_multiply_sum',
             'sum_arrays', 'set_species_name', 'set_nuclear_charge', 'set_nbi_unit_identity',
             'set_ec_beam_identity', 'set_ec_mode', 'check_profile_shape', 'fill_profiles1d_grid',
             'fill_properties', 'fill_code_properties', 'fill_vac_tor_field', 'fill_external_model_inputs',
             'build_code_parameters_xml', 'checkErrExit')


VIEWS = {'cdf': CdfView, 'namelist': NamelistView, 'dd': DdView, 'signs': SignsView, 'grid': GridView,
         'fill': FillView}


class ConversionContext:
    """The views of one converter, plus its mapping ledger and logger."""

    def __init__(self, owner):
        self.owner = owner
        for name, view in VIEWS.items():
            setattr(self, name, view(owner))

    @property
    def mapping(self):
        """The ``MappingLedger`` recording skipped, missing and broken mappings."""
        return self.owner.mapping

    @property
    def logger(self):
        return self.owner.logger

    @property
    def debug_max(self):
        return self.owner.debug_max

    @property
    def supports_occurrences(self):
        """True when the converter can submit an IDS at an explicit occurrence (``putdb_imas_profile``)."""
        return hasattr(self.owner, 'putdb_imas_profile')

    @classmethod
    def of(cls, owner):
        """The context of ``owner``: its own if it has one, else a new one over it (test doubles)."""
        context = getattr(owner, 'context', None)
        return context if isinstance(context, cls) else cls(owner)
