#!/usr/bin/python
"""The IMAS beta definitions, evaluated from TRANSP's pressure profile and the IDS's own reference field.

TRANSP's ``BTEQ`` is ``2 mu0 <PTOWB>_V / (BZXR^2 <1/R^2>)`` (checked to 1e-4 on every bundled case that has the
variables), i.e. normalised by the volume-averaged vacuum toroidal field squared, and ``BETAT``/``BPEQ`` use a
poloidal-field average that corresponds to no single radius.  The Data Dictionary defines the betas against the
IDS's own ``vacuum_toroidal_field`` (``b0`` at ``r0``), so copying TRANSP's numbers would label a different
definition with the IMAS name.  These functions compute the IMAS definitions from the same pressure:

* ``beta_tor      = 2 mu0 (int p dV / V) / B0^2``
* ``beta_pol      = 4 int p dV / (R0 mu0 Ip^2)``
* ``beta_tor_norm = 100 beta_tor a[m] B0[T] / Ip[MA]``

with ``B0`` and ``R0`` the exported ``b0`` and ``r0`` and ``a`` the boundary minor radius ``(Rmax - Rmin)/2``.
No IMAS or CDF dependency: arrays in, arrays out.
"""

import numpy as np

MU0 = 4.0e-7 * np.pi


def boundary_minor_radius(cos_moments, sin_moments=None, ntheta=720):
    """``(Rmax - Rmin)/2`` per time of the boundary ``R(theta) = sum_k c_k cos(k theta) + s_k sin(k theta)``.

    ``cos_moments`` has shape ``(K + 1, nt)``, row ``k`` the coefficient ``c_k`` (row 0 the constant term);
    ``sin_moments`` the same for ``s_k`` (row 0 is ignored) or ``None``.  Same unit in and out.  The
    reconstruction agrees with the minor radius of the traced LCFS outline to 1-2% on the bundled cases.
    """
    cos_moments = np.asarray(cos_moments, dtype=np.float64)
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    radius = np.zeros((cos_moments.shape[1], ntheta))
    for k in range(cos_moments.shape[0]):
        radius += cos_moments[k][:, None] * (np.cos(k * theta)[None, :] if k else 1.0)
    if sin_moments is not None:
        sin_moments = np.asarray(sin_moments, dtype=np.float64)
        for k in range(1, sin_moments.shape[0]):
            radius += sin_moments[k][:, None] * np.sin(k * theta)[None, :]
    return 0.5 * (radius.max(axis=1) - radius.min(axis=1))


def imas_betas(pressure_pa, dvol_m3, b0_t, r0_m, ip_a, minor_radius_m):
    """``{'beta_tor', 'beta_pol', 'beta_tor_norm'}`` over time, or ``{}`` when the inputs do not fit.

    ``pressure_pa`` and ``dvol_m3`` have shape ``(nt, nx)`` (zones inside the boundary); ``b0_t``, ``ip_a`` and
    ``minor_radius_m`` have length ``nt``; ``r0_m`` is the constant reference radius.  ``b0_t`` and ``ip_a`` keep
    their physical signs (they decide the sign of ``beta_tor_norm`` as in the Data Dictionary formula); the
    other two betas do not depend on them.
    """
    pressure = np.asarray(pressure_pa, dtype=np.float64)
    dvol = np.asarray(dvol_m3, dtype=np.float64)
    b0, ip, a = (np.asarray(v, dtype=np.float64) for v in (b0_t, ip_a, minor_radius_m))
    if (pressure.ndim != 2 or pressure.shape != dvol.shape
            or not (len(b0) == len(ip) == len(a) == pressure.shape[0])):
        return {}
    with np.errstate(divide='ignore', invalid='ignore'):
        integral = np.sum(pressure * dvol, axis=1)                      # int p dV  [J]
        volume = np.sum(dvol, axis=1)
        beta_tor = 2.0 * MU0 * integral / volume / b0 ** 2
        beta_pol = 4.0 * integral / (r0_m * MU0 * ip ** 2)
        beta_tor_norm = 100.0 * beta_tor * a * b0 / (ip / 1.0e6)
    return {'beta_tor': beta_tor, 'beta_pol': beta_pol, 'beta_tor_norm': beta_tor_norm}
