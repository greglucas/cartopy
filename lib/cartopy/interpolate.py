"""
Drawing heavily from the good work in d3-geo.
In particular, expanding to spherical resampling
https://observablehq.com/@jrus/sphere-resample
"""
import functools
import numpy as np
import pyproj
import shapely


def spherical(x, y, z):
    """
    Cartesian to spherical
    (x, y, z) -> (lon, lat)
    """
    return (np.rad2deg(np.arctan2(y, x)),
            np.rad2deg(np.arcsin(z)))


def cartesian(lon, lat):
    """
    Spherical to Cartesian
    (lon, lat) -> (x, y, z)
    """
    radians = np.pi/180
    return (np.cos(radians * lat) * np.cos(radians * lon),
            np.cos(radians * lat) * np.sin(radians * lon),
            np.sin(radians * lat))


def stereo_length(p0, p1):
    """
    Compute the stereographic length between cartesian points
    p1 -> p2 == (x0, y0, z0) -> (x1, y1, z1)
    """
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    pxy = x0 * (y1 - y0) - (x1 - x0) * y0
    pyz = y0 * (z1 - z0) - (y1 - y0) * z0
    pzx = z0 * (x1 - x0) - (z1 - z0) * x0
    q = x0*(x1 + x0) + y0*(y1 + y0) + z0*(z1 + z0)
    if not q*q:
        # adding !(q*q) means q==0 => return Infinity
        return np.inf
    return (pxy*pxy + pyz*pyz + pzx*pzx) / (q*q)


def stereo_area(p0, p1, p2):
    """
    Compute the stereographic area between three points
    """
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    p = (x0*((y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0)) +
         y0*((z1 - z0)*(x2 - x0) - (z2 - z0)*(x1 - x0)) +
         z0*((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0)))
    q = (x0 + x2)*(x0 + x1) + (y0 + y2)*(y0 + y1) + (z0 + z2)*(z0 + z1)
    q2 = q*q
    if not q2:
        # adding !(q*q) means q==0 => return Infinity
        return np.inf
    return p*p / (q2)


def planar_midpoint(p0, p1):
    """
    Calculate the planar midpoint
    """
    return (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]))


maxDepth = 16
radians = np.pi/180
delta = 0.001
dd = np.tan(radians*delta/2)**2


def resample_line(w0, u0, w1, u1, ll01, transformer, depth, array):
    """
    w == input coordinates
    u == cartesian coordinates on the sphere
    ll == stereographic distance squared
    depth == current depth
    array == array of coordinates
    """
    depth -= 1
    if not depth:
        return

    w2 = planar_midpoint(w0, w1)
    # Transforming the midpoint
    point2 = transformer.transform(*w2)
    u2 = cartesian(*point2)
    # u2 = cartesian(*w2)
    ll02 = stereo_length(u2, u0)
    ll12 = stereo_length(u2, u1)
    AA = stereo_area(u2, u0, u1)
    hh = AA * (1 + 0.25*ll01)*(1 + 0.25*ll01) / (dd * ll01) # perpendicular projected distance
    ww = 2 * ((ll02 - ll12) / ll01) * ((ll02 - ll12) / ll01) # midpoint roughly in the middle
    # (ll02 + ll12 > 0.25) forces bisection for over-long segments, and handles ll01 == Infinity edge case
    # (ll02 + ll12 > dd)) stops bisection when the segment gets tiny
    if (((hh + ww > 1) and (ll02 + ll12 > dd)) or (ll02 + ll12 > 0.25)):
        resample_line(w0, u0, w2, u2, ll02, transformer, depth, array)
        array.append(point2)
        resample_line(w2, u2, w1, u1, ll12, transformer, depth, array)


def resample_chain(point_array, transformer):
    outarray = []
    ws = np.array(point_array)
    w0 = point_array[0]
    # Transform the entire array at once
    ws_x, ws_y = transformer.transform(ws[:, 0], ws[:, 1])
    u0 = cartesian(ws_x[0], ws_y[0])
    outarray.append((ws_x[0], ws_y[0]))
    for i in range(1, len(ws_x)):
        w1 = point_array[i]
        # Use the transformed coordinates
        u1 = cartesian(ws_x[i], ws_y[i])
        resample_line(w0, u0, w1, u1, stereo_length(u0, u1),
                      transformer, maxDepth, outarray)
        outarray.append((ws_x[i], ws_y[i]))
        w0 = w1
        u0 = u1
    return outarray

@functools.cache
def _get_transformers(proj_from, proj_to):
    """
    Create the transformer objects going from one projection to another.
    
    Returns two transformer objects, one from the original projection to the
    geodetic representation, and the second transformer going from the
    geodetic representation to the final projection.
    """
    geodetic = proj_from.as_geodetic()
    transformer1 = pyproj.Transformer.from_crs(
        proj_from, geodetic, always_xy=True)
    transformer2 = pyproj.Transformer.from_crs(
        geodetic, proj_to, always_xy=True)
    return transformer1, transformer2


def transform_geometry(geom, proj_data, proj_map):
    """Transform an input geometry from proj_data to proj_map."""
    if geom.type not in ("LineString", "LinearRing"):
        raise NotImplementedError("Geometry not handled yet")

    # Set up our initial transformer, going from proj_data to the sphere
    # Going to the geodetic transform of our initial projection
    print("Transforming geometry")
    transformer1, transformer2 = _get_transformers(proj_data, proj_map)
    coords = resample_chain(geom.coords, transformer1)

    print(f"IN: {len(geom.coords)}, OUT: {len(coords)}")
    # Our coords now has our interpolated points added
    # We want to take those coords and map them to our desired map projection
    # move from lat/lon coords to our initial projection
    coords_arr = np.array(coords)
    coords = transformer2.transform(coords_arr[:, 0], coords_arr[:, 1])
    return type(geom)(np.stack(coords).T)
