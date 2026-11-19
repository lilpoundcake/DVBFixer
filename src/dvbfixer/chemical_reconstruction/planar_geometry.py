"""Planar exterior placement from explicit equal-stiffness angle targets."""
from __future__ import annotations

import math

from .measurement import Point, _dot, _sub, _unit, angle_degrees


def planar_exterior_tip(centre: Point, first: Point, second: Point, *, length: float,
                         first_angle: float, second_angle: float) -> Point:
    u = _unit(_sub(first, centre))
    direction = _unit(_sub(second, centre))
    dot = _dot(direction, u)
    v = _unit((direction[0]-dot*u[0], direction[1]-dot*u[1], direction[2]-dot*u[2]))
    alpha = angle_degrees(first, centre, second)
    theta = (360-alpha+first_angle-second_angle)/2
    other = 360-alpha-theta
    if not (0 < theta < 180 and 0 < other < 180):
        raise ValueError("planar-exterior-angle-unavailable")
    cosine, sine = math.cos(math.radians(theta)), math.sin(math.radians(theta))
    return (centre[0]+length*(cosine*u[0]-sine*v[0]),
            centre[1]+length*(cosine*u[1]-sine*v[1]),
            centre[2]+length*(cosine*u[2]-sine*v[2]))
