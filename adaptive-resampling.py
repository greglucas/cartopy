import matplotlib.pyplot as plt
import numpy as np
import pyproj

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
delta = 0.01
# delta = 55660
dd = np.tan(radians*delta/2)**2


def resample_line(w0, u0, w1, u1, ll01, depth, array):
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
    point2 = proj.transform(*w2)
    # u2 = cartesian(*point2)
    u2 = cartesian(*w2)
    ll02 = stereo_length(u2, u0)
    ll12 = stereo_length(u2, u1)
    AA = stereo_area(u2, u0, u1)
    hh = AA * (1 + 0.25*ll01)*(1 + 0.25*ll01) / (dd * ll01) # perpendicular projected distance
    ww = 2 * ((ll02 - ll12) / ll01) * ((ll02 - ll12) / ll01) # midpoint roughly in the middle
    # (ll02 + ll12 > 0.25) forces bisection for over-long segments, and handles ll01 == Infinity edge case
    # (ll02 + ll12 > dd)) stops bisection when the segment gets tiny
    if (((hh + ww > 1) and (ll02 + ll12 > dd)) or (ll02 + ll12 > 0.25)):
        resample_line(w0, u0, w2, u2, ll02, depth, array)
        array.append(point2)
        resample_line(w2, u2, w1, u1, ll12, depth, array)


def resample_chain(point_array):
    outarray = []
    w0 = point_array[0]
    point0 = proj.transform(*w0)
    # u0 = cartesian(*point0)
    u0 = cartesian(*w0)
    outarray.append(point0)
    for i in range(1, len(point_array)):
        w1 = point_array[i]
        point1 = proj.transform(*w1)
        # u1 = cartesian(*point1)
        u1 = cartesian(*w1)
        resample_line(w0, u0, w1, u1, stereo_length(u0, u1), maxDepth, outarray)
        outarray.append(point1)
        w0 = w1
        u0 = u1
    return outarray


points = [[2, 40], [30, 7], [-50, -50]]
# points = [[-70.67, -33.45], # Santiago
#           [-21.88 - 360*4, 64.13]]

out = resample_chain(points)
xs = [x[0] for x in out]
ys = [x[1] for x in out]

print(len(out))
out2 = resample_chain(points[::-1])
xs2 = [x[0] for x in out2]
ys2 = [x[1] for x in out2]

print(out == out2[::-1])

fig, ax = plt.subplots()

# ax.plot([xs[0], xs[-1]], [ys[0], ys[-1]], c='k')
# ax.plot(xs, ys, c='r', markersize=5, marker='o')
# ax.plot(xs2, ys2, c='b', markersize=2, marker='o')

# plt.show()
scale = '110m'
geoms = cfeature.COASTLINE.with_scale(scale).geometries()
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines(scale)
# ax.set_global()
ax = plt.axes()

for i in range(100):
    coast1 = next(geoms)

    points = coast1.coords[:]
    out = points
    xs = [x[0] for x in out]
    ys = [x[1] for x in out]

    out = resample_chain(points)
    if len(points) != len(out):
        print("BEFORE:", len(points))
        print("AFTER", len(out))
    xs2 = [x[0] for x in out]
    ys2 = [x[1] for x in out]

    xs2, ys2 = proj_i.transform(xs2, ys2)

    ax.plot(xs, ys, c='r', markersize=10, marker='o')
    ax.plot(xs2, ys2, c='b', markersize=8, marker='o')

plt.show()
