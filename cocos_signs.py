#!/usr/bin/python
"""Physical direction and COCOS sign policy of a conversion, separate from the converter object.

TRANSP stores the field-orientation controls separately from many CDF quantities: ``NLJCCW`` says
whether the plasma current is counter-clockwise seen from above and ``NLBCCW`` says the same for the
toroidal magnetic field.  IMAS uses +phi counter-clockwise from above.

IMAS DD3 uses COCOS 11 and DD4 uses COCOS 17.  They share the cylindrical/toroidal orientation and
the q sign convention, but the poloidal-flux sign is reversed: for positive Ip, COCOS 17 psi
decreases from the magnetic axis to the boundary.

``resolve_signs`` turns the two directions and the DD major version into one immutable
``CocosSigns``; it has no IMAS or NetCDF dependency.
"""

from dataclasses import dataclass


def direction_override(value):
    """``True`` for ccw, ``False`` for cw, ``None`` for auto or anything else."""
    if value is None:
        return None
    text = str(value).strip().lower()
    if text == 'ccw':
        return True
    if text == 'cw':
        return False
    return None


@dataclass(frozen=True)
class CocosSigns:
    nljccw: bool
    nlbccw: bool
    ip_sign: float
    bt_sign: float
    target_cocos: int
    # physical/vector quantities in the IMAS +phi (counter-clockwise) direction
    jphi_sign: float
    jparallel_sign: float
    phi_sign: float
    f_sign: float
    q_sign: float
    vloop_sign: float
    toroidal_rotation_sign: float
    poloidal_field_sign: float
    psi_sign: float
    # maps the target theta to the raw TRANSP moment angle atan2(Z-Zaxis, R-Raxis)
    theta_raw_from_target: float
    missing: tuple = ()          # control(s) that were neither overridden nor in the namelist

    def summary(self):
        return ('COCOS/signs:'
                f' target={self.target_cocos}'
                f' NLJCCW={self.nljccw} (Ip sign {self.ip_sign:+.0f})'
                f' NLBCCW={self.nlbccw} (Bphi sign {self.bt_sign:+.0f})'
                f' q sign {self.q_sign:+.0f}'
                f' psi sign {self.psi_sign:+.0f}')

    def as_attributes(self):
        """Field name -> value, without ``missing`` (the attributes the converter exposes)."""
        return {name: getattr(self, name) for name in self.__dataclass_fields__ if name != 'missing'}


def resolve_signs(ip_direction='auto', bt_direction='auto', namelist_bool=None, dd_major=4):
    """Resolve the sign policy.

    ``ip_direction``/``bt_direction`` are 'ccw', 'cw' or 'auto'; ``auto`` reads ``NLJCCW``/``NLBCCW``
    through ``namelist_bool(name)`` (returning ``None`` when absent).  A control that is still unknown
    is assumed counter-clockwise and reported in ``missing``.
    """
    namelist_bool = namelist_bool or (lambda name: None)
    jccw = direction_override(ip_direction)
    bccw = direction_override(bt_direction)
    if jccw is None:
        jccw = namelist_bool('NLJCCW')
    if bccw is None:
        bccw = namelist_bool('NLBCCW')
    missing = []
    if jccw is None:
        missing.append('NLJCCW')
        jccw = True
    if bccw is None:
        missing.append('NLBCCW')
        bccw = True

    nljccw, nlbccw = bool(jccw), bool(bccw)
    ip_sign = 1.0 if nljccw else -1.0
    bt_sign = 1.0 if nlbccw else -1.0
    try:
        major = int(str(dd_major).split('.')[0])
    except (ValueError, TypeError, AttributeError):
        major = 4
    target_cocos = 17 if major >= 4 else 11
    return CocosSigns(
        nljccw=nljccw, nlbccw=nlbccw, ip_sign=ip_sign, bt_sign=bt_sign, target_cocos=target_cocos,
        jphi_sign=ip_sign, jparallel_sign=ip_sign, phi_sign=bt_sign, f_sign=bt_sign,
        q_sign=ip_sign * bt_sign, vloop_sign=ip_sign, toroidal_rotation_sign=ip_sign,
        poloidal_field_sign=ip_sign,
        # COCOS 11 -> 17 reverses the sign of the full poloidal flux and leaves Ip, Bphi and q alone
        psi_sign=ip_sign if target_cocos == 11 else -ip_sign,
        # TRANSP moments use the counter-clockwise atan2; DD4/COCOS 17 uses its negative
        theta_raw_from_target=1.0 if target_cocos == 11 else -1.0,
        missing=tuple(missing))
