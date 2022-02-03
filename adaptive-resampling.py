import math

import matplotlib.pyplot as plt
import numpy as np
import pyproj
import shapely

import cartopy.crs as ccrs
import cartopy.feature as cfeature

proj = pyproj.Transformer.from_crs(4326, 3857, always_xy=True)
proj = pyproj.Transformer.from_crs(32662, 3857, always_xy=True)
proj = pyproj.Transformer.from_crs(3395, 4326, always_xy=True)
proj_i = pyproj.Transformer.from_crs(4326, 3395, always_xy=True)


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
radius = 6378137
radius = 1
delta *= radius
# delta = 55660
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


def resample_chain(point_array, proj_from, proj_to):
    outarray = []
    w0 = point_array[0]
    transformer = pyproj.Transformer.from_crs(proj_from, proj_to, always_xy=True)
    point0 = transformer.transform(*w0)
    u0 = cartesian(*point0)
    # u0 = cartesian(*w0)
    outarray.append(point0)
    for i in range(1, len(point_array)):
        w1 = point_array[i]
        point1 = transformer.transform(*w1)
        u1 = cartesian(*point1)
        # u1 = cartesian(*w1)
        resample_line(w0, u0, w1, u1, stereo_length(u0, u1),
                      transformer, maxDepth, outarray)
        outarray.append(point1)
        w0 = w1
        u0 = u1
    return outarray


points = [[-10, -10], [0, 30], [40, 30], [30, -10], [-10, -10]]
# points = [[-70.67, -33.45], # Santiago
#           [-21.88 - 360*4, 64.13]]

# WGS84 lat/lon
latlon = 4326
# Mercator
mercator = 3395
# Equirectangular (Plate Carree)
platecarree = 32662
# Robinson
robinson = "World_Robinson"
# Orthographic
orthographic = 9840

transformer = pyproj.Transformer.from_crs(latlon, robinson, always_xy=True)
# move from lat/lon coords to our initial projection
points_projected = [transformer.transform(*p) for p in points]

out = resample_chain(points_projected, proj_from=robinson, proj_to=latlon)
xs = [x[0] for x in out]
ys = [x[1] for x in out]

print(f"# input points: {len(points)}\n# of output points: {len(out)}")
out2 = resample_chain(
    points_projected[::-1], proj_from=robinson, proj_to=latlon)
xs2 = [x[0] for x in out2]
ys2 = [x[1] for x in out2]

print("Forward/Reverse the same?", out == out2[::-1])

fig, ax = plt.subplots()

ax.plot([xs[0], xs[-1]], [ys[0], ys[-1]], c='k')
ax.plot(xs, ys, c='r', markersize=5, marker='o')
ax.set_title("Lat/lon WGS84 coordinates")

# Now we want to go from Lat/Lon to our output grid
transformer = pyproj.Transformer.from_crs(latlon, orthographic, always_xy=True)
fig, ax2 = plt.subplots()

points_projected = [transformer.transform(*p) for p in points]
ax2.set_title("Output projection")
projected_x = [x[0] for x in points_projected]
projected_y = [x[1] for x in points_projected]
ax2.plot(projected_x, projected_y, c='k', markersize=5, marker='o')

# Now convert our "out" (lat/lon) to our final coordinate system
out_projected = [transformer.transform(*p) for p in out]
xs = [x[0] for x in out_projected]
ys = [x[1] for x in out_projected]
ax2.plot(xs, ys, c='r', markersize=5, marker='o')

plt.show()



# scale = '110m'
# geoms = cfeature.COASTLINE.with_scale(scale).geometries()
# # ax = plt.axes(projection=ccrs.PlateCarree())
# # ax.coastlines(scale)
# # ax.set_global()
# ax = plt.axes()

# for i in range(100):
#     coast1 = next(geoms)

#     points = coast1.coords[:]
#     out = points
#     xs = [x[0] for x in out]
#     ys = [x[1] for x in out]

#     out = resample_chain(points)
#     if len(points) != len(out):
#         print("BEFORE:", len(points))
#         print("AFTER", len(out))
#     xs2 = [x[0] for x in out]
#     ys2 = [x[1] for x in out]

#     xs2, ys2 = proj_i.transform(xs2, ys2)

#     ax.plot(xs, ys, c='r', markersize=10, marker='o')
#     ax.plot(xs2, ys2, c='b', markersize=8, marker='o')

# plt.show()

def crosses_antimeridian(p0, p1, threshold=180):
    """
    Tests whether the two points cross the antimeridian.
    
    Parameters
    ----------
    p0, p1 : point coordinates
        Points of a geometry in degrees latitude, longitude.
    threshold : float
        The threshold separation distance to consider for whether two points
        have crossed the antimeridian.
    """
    # Test the x coordinates
    return abs(p1[0] - p0[0]) > threshold


def shift_geometry(geom):
    """
    Convenience function to shift a geometry away from the antimeridian
    and avoid issues with wrapping near the transition.

    If the geometry is outside of (-180, 180), add or subtract 360 degrees
    to shift it to the other side of the domain.
    """
    xmin, _, xmax, _ = geom.bounds
    if xmin < -180:
        offset = 360
    elif xmax > 180:
        offset = -360
    else:
        # geometry is fully within the bounds
        return geom
    # Translate as appropriate
    return shapely.affinity.translate(geom, xoff=offset)


def split_antimeridian(geom):
    """Split a geometry at the antimeridian"""
    shell_minx = shell_maxx = None
    split_meridians = set()
    splitter = None

    orig_coords = geom.coords.copy()

    # Iterate over the rings within our Polygon
    for ring_index, ring in enumerate(geom):
        if len(ring) < 1:
            # Ignore empty geoms
            continue
        ring_minx = ring_maxx = ring[0][0]
        crossings = 0

        # Iterate over the coordinates within the ring testing
        # the consecutive longitudes
        for i in range(1, len(ring)):
        # for coord_index, (lon, _) in enumerate(ring[1:], start=1):
            lon0 = ring[i - 1][0]
            lon1 = ring[i][0]
            if crosses_antimeridian(lon1, lon0):
                direction = math.copysign(1, lon1 - lon0)
                orig_coords[ring_index][i][0] = lon1 - (direction * 360)
                crossings += 1

            x_shift = orig_coords[ring_index][i][0]
            if x_shift < ring_minx:
                ring_minx = x_shift
            if x_shift > ring_maxx:
                ring_maxx = x_shift

        # Ensure that any holes remain contained within the (translated) outer shell
        if (ring_index == 0):
            # The first ring is the outer shell of the Polygon
            shell_minx, shell_maxx = (ring_minx, ring_maxx)
        elif (ring_minx < shell_minx):
            # This ring is to the left of the outer shell and needs to be
            # shifted to the right
            ring_shift = [[x + 360, y] for (x, y) in orig_coords[ring_index]]
            orig_coords[ring_index] = ring_shift
            ring_minx, ring_maxx = (ring_minx + 360, ring_maxx + 360)
        elif (ring_maxx > shell_maxx):
            # This ring is to the right of the outer shell and needs to be
            # shifted to the left
            ring_shift = [[x - 360, y] for (x, y) in orig_coords[ring_index]]
            orig_coords[ring_index] = ring_shift
            ring_minx, ring_maxx = (ring_minx - 360, ring_maxx - 360)

        if crossings:  # keep track of meridians to split on
            if ring_minx < -180:
                split_meridians.add(-180)
            if ring_maxx > 180:
                split_meridians.add(180)

    n_splits = len(split_meridians)
    if n_splits > 1:
        raise NotImplementedError(
            """Splitting a Polygon by multiple meridians (MultiLineString) 
               not supported by Shapely"""
        )
    elif n_splits == 1:
        split_lon = next(iter(split_meridians))
        meridian = [[split_lon, -90.0], [split_lon, 90.0]]
        splitter = shapely.LineString(meridian)
    else:
        return geom

    shell, *holes = orig_coords
    split_polygons = shapely.split(shapely.Polygon(shell, holes), splitter)

    return shift_geometry(split_polygons)
