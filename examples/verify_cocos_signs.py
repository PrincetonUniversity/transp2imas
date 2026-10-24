#!/usr/bin/python
"""Verify the COCOS/direction sign policy (no IMAS or NetCDF needed).

The policy is a pure function of three inputs, the plasma-current direction (``NLJCCW``), the toroidal
field direction (``NLBCCW``) and the DD major version (DD3 -> COCOS 11, DD4 -> COCOS 17).  The expected
values below are written out by hand from the definitions, not computed with the code under test:

* Ip and Bphi signs follow the directions; q = Ip*Bphi sign; F and phi follow Bphi; j, v_loop and
  toroidal rotation follow Ip,
* going from COCOS 11 to 17 reverses only the poloidal-flux sign, and the moment angle,
* overrides win over the namelist, an unknown control is assumed counter-clockwise and reported.

    python examples/verify_cocos_signs.py
"""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cocos_signs import CocosSigns, direction_override, resolve_signs

# (NLJCCW, NLBCCW, DD major) -> (ip, bt, q, psi, cocos, theta_raw_from_target)
TABLE = {
    (True, True, 3): (+1, +1, +1, +1, 11, +1),
    (True, False, 3): (+1, -1, -1, +1, 11, +1),
    (False, True, 3): (-1, +1, -1, -1, 11, +1),
    (False, False, 3): (-1, -1, +1, -1, 11, +1),
    (True, True, 4): (+1, +1, +1, -1, 17, -1),
    (True, False, 4): (+1, -1, -1, -1, 17, -1),
    (False, True, 4): (-1, +1, -1, +1, 17, -1),
    (False, False, 4): (-1, -1, +1, +1, 17, -1),
}


def namelist(jccw, bccw):
    return {'NLJCCW': jccw, 'NLBCCW': bccw}.get


class SignTableTests(unittest.TestCase):
    def test_every_direction_and_dd_combination(self):
        for (jccw, bccw, major), (ip, bt, q, psi, cocos, theta) in TABLE.items():
            with self.subTest(jccw=jccw, bccw=bccw, dd=major):
                s = resolve_signs('auto', 'auto', namelist(jccw, bccw), dd_major=major)
                self.assertEqual((s.ip_sign, s.bt_sign, s.q_sign, s.psi_sign, s.target_cocos, s.theta_raw_from_target),
                                 (ip, bt, q, psi, cocos, theta))
                # vector quantities follow Ip or Bphi
                self.assertEqual((s.jphi_sign, s.jparallel_sign, s.vloop_sign, s.toroidal_rotation_sign,
                                  s.poloidal_field_sign), (ip,) * 5)
                self.assertEqual((s.phi_sign, s.f_sign), (bt, bt))
                self.assertEqual((s.nljccw, s.nlbccw), (jccw, bccw))
                self.assertEqual(s.missing, ())

    def test_cocos_11_to_17_flips_only_the_flux_sign_and_the_angle(self):
        for jccw in (True, False):
            for bccw in (True, False):
                a = resolve_signs('auto', 'auto', namelist(jccw, bccw), dd_major=3)
                b = resolve_signs('auto', 'auto', namelist(jccw, bccw), dd_major=4)
                self.assertEqual(a.psi_sign, -b.psi_sign)
                self.assertEqual(a.theta_raw_from_target, -b.theta_raw_from_target)
                for name in ('ip_sign', 'bt_sign', 'q_sign', 'phi_sign', 'f_sign', 'jphi_sign', 'vloop_sign'):
                    self.assertEqual(getattr(a, name), getattr(b, name), name)

    def test_dd_major_accepts_strings_and_falls_back_to_dd4(self):
        self.assertEqual(resolve_signs(dd_major='3').target_cocos, 11)
        self.assertEqual(resolve_signs(dd_major='4.1.1').target_cocos, 17)
        self.assertEqual(resolve_signs(dd_major='x').target_cocos, 17)
        self.assertEqual(resolve_signs(dd_major=None).target_cocos, 17)


class OverrideTests(unittest.TestCase):
    def test_direction_override_values(self):
        self.assertIs(direction_override('ccw'), True)
        self.assertIs(direction_override(' CW '), False)
        for value in (None, 'auto', '', 'clockwise', 3):
            self.assertIsNone(direction_override(value))

    def test_override_beats_the_namelist(self):
        s = resolve_signs('cw', 'ccw', namelist(True, False), dd_major=4)
        self.assertEqual((s.nljccw, s.nlbccw, s.ip_sign, s.bt_sign), (False, True, -1, +1))

    def test_unknown_controls_are_assumed_ccw_and_reported(self):
        s = resolve_signs('auto', 'auto', lambda name: None, dd_major=4)
        self.assertEqual((s.ip_sign, s.bt_sign, s.missing), (1.0, 1.0, ('NLJCCW', 'NLBCCW')))
        s = resolve_signs('cw', 'auto', lambda name: None, dd_major=4)
        self.assertEqual(s.missing, ('NLBCCW',))
        s = resolve_signs('cw', 'cw', None, dd_major=4)               # no namelist at all, nothing missing
        self.assertEqual((s.missing, s.ip_sign, s.bt_sign), ((), -1.0, -1.0))


class PublishedFormTests(unittest.TestCase):
    def test_summary_line_matches_what_earlier_conversions_logged(self):
        # copied from the conversion logs of the bundled cases before the policy was extracted
        s = resolve_signs('auto', 'auto', namelist(True, False), dd_major=3)
        self.assertEqual(s.summary(), 'COCOS/signs: target=11 NLJCCW=True (Ip sign +1) NLBCCW=False '
                                      '(Bphi sign -1) q sign -1 psi sign +1')
        s = resolve_signs('auto', 'auto', namelist(False, False), dd_major=3)
        self.assertEqual(s.summary(), 'COCOS/signs: target=11 NLJCCW=False (Ip sign -1) NLBCCW=False '
                                      '(Bphi sign -1) q sign +1 psi sign -1')

    def test_attributes_exclude_missing_and_the_object_is_immutable(self):
        s = resolve_signs(namelist_bool=namelist(True, True), dd_major=4)
        attributes = s.as_attributes()
        self.assertNotIn('missing', attributes)
        self.assertEqual(set(attributes), set(CocosSigns.__dataclass_fields__) - {'missing'})
        with self.assertRaises(Exception):
            s.ip_sign = -1.0


if __name__ == '__main__':
    unittest.main(verbosity=1)
