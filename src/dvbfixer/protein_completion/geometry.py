"""Small stereochemical constructions in angstroms, independent of templates."""
from __future__ import annotations

import math

from dvbfixer.protein_completion.contracts import Point


def sub(a: Point, b: Point) -> Point:
    return a[0]-b[0], a[1]-b[1], a[2]-b[2]


def add(a: Point, b: Point) -> Point:
    return a[0]+b[0], a[1]+b[1], a[2]+b[2]


def scaled(a: Point, factor: float) -> Point:
    return factor*a[0], factor*a[1], factor*a[2]


def dot(a: Point, b: Point) -> float:
    return sum(x*y for x, y in zip(a, b))


def cross(a: Point, b: Point) -> Point:
    return a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]


def norm(a: Point) -> float:
    return math.sqrt(dot(a, a))


def unit(a: Point) -> Point:
    length = norm(a)
    if not math.isfinite(length) or length < 1e-8:
        raise ValueError("degenerate anchor direction")
    return scaled(a, 1/length)


def l_cb(n: Point, ca: Point, c: Point) -> Point:
    """Solve two angle constraints; select only the positive N-CA-C-CB volume.

    Illustrative geometry v1: 1.53 A, N-CA-CB=110.5, C-CA-CB=110.1 degrees.
    This is not a residue/backbone-dependent calibrated geometry library.
    """
    u, v = unit(sub(n, ca)), unit(sub(c, ca))
    normal = unit(cross(u, v))
    cosine = dot(u, v)
    determinant = 1-cosine*cosine
    if determinant < 1e-8:
        raise ValueError("collinear backbone anchors")
    x = 1.53*math.cos(math.radians(110.5))
    y = 1.53*math.cos(math.radians(110.1))
    alpha, beta = (x-cosine*y)/determinant, (y-cosine*x)/determinant
    parallel = add(scaled(u, alpha), scaled(v, beta))
    height2 = 1.53**2-dot(parallel, parallel)
    if height2 <= 1e-8:
        raise ValueError("incompatible backbone angles")
    return add(ca, add(parallel, scaled(normal, math.sqrt(height2))))


def extend(a: Point, b: Point, c: Point, length: float, angle: float, chi: float) -> Point:
    """Position D with B-C-D angle and signed A-B-C-D dihedral chi."""
    axis = unit(sub(c, b))
    normal = unit(cross(sub(b, a), axis))
    plane = cross(normal, axis)
    theta, phi = math.radians(angle), math.radians(chi)
    direction = add(scaled(axis, -math.cos(theta)),
                    scaled(add(scaled(plane, math.cos(phi)), scaled(normal, math.sin(phi))),
                           math.sin(theta)))
    return add(c, scaled(direction, length))


def two_angle_tip(center: Point, first: Point, second: Point, *, length: float,
                  first_angle: float, second_angle: float, orientation: int) -> Point:
    """Intersect two angular cones; choose a declared signed-volume hemisphere.

    Anchors are immutable. No Cartesian fitting, observed-atom reflection,
    stochastic fallback, or imaginary/nearplanar solution is permitted.
    """
    if (type(orientation) is not int or orientation not in {-1, 1}
            or any(isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v)
                   for v in (length, first_angle, second_angle))
            or not 0 < length < 3 or not 0 < first_angle < 180 or not 0 < second_angle < 180):
        raise ValueError("bounded length, proper angles and explicit orientation required")
    u, v = unit(sub(first, center)), unit(sub(second, center))
    normal = unit(cross(u, v))
    cosine = dot(u, v)
    determinant = 1-cosine*cosine
    if determinant < 1e-8:
        raise ValueError("collinear observed branch anchors")
    x = length*math.cos(math.radians(first_angle))
    y = length*math.cos(math.radians(second_angle))
    alpha, beta = (x-cosine*y)/determinant, (y-cosine*x)/determinant
    parallel = add(scaled(u, alpha), scaled(v, beta))
    height2 = length*length-dot(parallel, parallel)
    if height2 <= 1e-8:
        raise ValueError("incompatible or nearplanar observed branch anchors")
    return add(center, add(parallel, scaled(normal, orientation*math.sqrt(height2))))
