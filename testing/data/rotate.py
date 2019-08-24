#!/usr/bin/env python3
# ------------------------------------------------------------------
# Filename: rotate.py
#  Purpose: Various Seismogram Rotation Functions
#   Author: Tobias Megies, Tom Richter, Lion Krischer
#    Email: tobias.megies@geophysik.uni-muenchen.de
#
# Copyright (C) 2009-2013 Tobias Megies, Tom Richter, Lion Krischer
# --------------------------------------------------------------------
"""
Various Seismogram Rotation Functions

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
from math import cos, sin, radians

def rotate_zne_lqt(z, n, e, ba, inc):
    """
    Rotates all components of a seismogram.

    The components will be rotated from ZNE (Z, North, East, left-handed) to
    LQT (e.g. ray coordinate system, right-handed). The rotation angles are
    given as the back-azimuth and inclination.

    The transformation consists of 3 steps::

        1. mirroring of E-component at ZN plain: ZNE -> ZNW
        2. negative rotation of coordinate system around Z-axis with angle ba:
           ZNW -> ZRT
        3. negative rotation of coordinate system around T-axis with angle inc:
           ZRT -> LQT

    :type z: :class:`~numpy.ndarray`
    :param z: Data of the Z component of the seismogram.
    :type n: :class:`~numpy.ndarray`
    :param n: Data of the North component of the seismogram.
    :type e: :class:`~numpy.ndarray`
    :param e: Data of the East component of the seismogram.
    :type ba: float
    :param ba: The back azimuth from station to source in degrees.
    :type inc: float
    :param inc: The inclination of the ray at the station in degrees.
    :return: L-, Q- and T-component of seismogram.

    This is from obspy and is distributed under GPL3
    """
    if len(z) != len(n) or len(z) != len(e):
        raise TypeError("Z, North and East component have different length!?!")
    if ba < 0 or ba > 360:
        raise ValueError("Back Azimuth should be between 0 and 360 degrees!")
    if inc < 0 or inc > 360:
        raise ValueError("Inclination should be between 0 and 360 degrees!")
    ba = radians(ba)
    inc = radians(inc)
    l = z * cos(inc) - n * sin(inc) * cos(ba) - e * sin(inc) * sin(ba)  # NOQA
    q = z * sin(inc) + n * cos(inc) * cos(ba) + e * cos(inc) * sin(ba)  # NOQA
    t = n * sin(ba) - e * cos(ba)  # NOQA
    return l, q, t

def rotate_ne_rt(n, e, ba):
    """
    Rotates horizontal components of a seismogram.

    The North- and East-Component of a seismogram will be rotated in Radial
    and Transversal Component. The angle is given as the back-azimuth, that is
    defined as the angle measured between the vector pointing from the station
    to the source and the vector pointing from the station to the North.

    :type n: :class:`~numpy.ndarray`
    :param n: Data of the North component of the seismogram.
    :type e: :class:`~numpy.ndarray`
    :param e: Data of the East component of the seismogram.
    :type ba: float
    :param ba: The back azimuth from station to source in degrees.
    :return: Radial and Transversal component of seismogram.

    This is from obspy and distributed under GPL3
    """
    if len(n) != len(e):
        raise TypeError("North and East component have different length.")
    if ba < 0 or ba > 360:
        raise ValueError("Back Azimuth should be between 0 and 360 degrees.")
    ba = radians(ba)
    r = - e * sin(ba) - n * cos(ba)
    t = - e * cos(ba) + n * sin(ba)
    return r, t

np.random.seed(4083)
vertical = np.random.randn(100)
north = np.random.randn(len(vertical))
east = np.random.randn(len(vertical))
baz = 131
radial, transverse = rotate_ne_rt(north, east, baz)
aoi = 72
baz = 333
l, q, t = rotate_zne_lqt(vertical, north, east, baz, aoi)

ofl = open("rotate_zne_rt_lqt.txt", "w")
for i in range(len(vertical)):
    ofl.write("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n"%( 
              vertical[i], north[i], east[i],
              radial[i], transverse[i],
              l[i], q[i], t[i]))
ofl.close()
