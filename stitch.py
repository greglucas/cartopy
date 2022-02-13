import math

epsilon = 1e-4
epsilonInverse = 1e4
x0 = -180
x0e = x0 + epsilon
x1 = 180
x1e = x1 - epsilon
y0 = -90
y0e = y0 + epsilon
y1 = 90
y1e = y1 - epsilon


def nonempty(coordinates):
    return len(coordinates) > 0


def quantize(x):
    return math.floor(x * epsilonInverse) / epsilonInverse


def normalizePoint(y):
    # pole or antimeridian?
    if y == y0 or y == y1:
        return [0, y]
    return [x0, quantize(y)]


def clampPoint(p):
    x, y = p
    clamped = False
    if (x <= x0e):
        x = x0
        clamped = True
    elif (x >= x1e):
        x = x1
        clamped = True

    if (y <= y0e):
        y = y0
        clamped = True
    elif (y >= y1e):
        y = y1
        clamped = True

    if clamped:
        return [x, y]
    return p


def clampPoints(points):
    # map over the clamped points
    return [clampPoint(point) for point in points]


# For each ring, detect where it crosses the antimeridian or pole.
def extractFragments(rings, polygon, fragments):
    for j in range(len(rings)):
        ring = rings[j]

        # By default, assume that this ring doesn’t need any stitching.
        fragments.append({"index": -1, "polygon": polygon, "ring": ring})

        n = len(ring)
        i = 0
        while i < n:
            point = ring[i]
            x, y = point

            # If this is an antimeridian or polar point…
            if (x <= x0e or x >= x1e or y <= y0e or y >= y1e):
                ring[i] = clampPoint(point)

                # Advance through any antimeridian or polar points…
                k = i + 1
                while k < n:
                    pointk = ring[k]
                    xk, yk = pointk
                    if (xk > x0e and xk < x1e and yk > y0e and yk < y1e):
                        break
                    k += 1

                # If this was just a single antimeridian or polar point,
                # we don’t need to cut this ring into a fragment;
                # we can just leave it as-is.
                if (k == i + 1):
                    i += 1
                    continue

                # Otherwise, if this is not the first point in the ring,
                # cut the current fragment so that it ends at the current point.
                # The current point is also normalized for later joining.
                if i:
                    fragmentBefore = {"index": -1,
                                      "polygon": polygon, "ring": ring[:i + 1]}
                    fragmentBefore["ring"][-1] = normalizePoint(y)
                    fragments[-1] = fragmentBefore

                # If the ring started with an antimeridian fragment,
                # we can ignore that fragment entirely.
                else:
                    fragments.pop()

                # If the remainder of the ring is an antimeridian fragment,
                # move on to the next ring.
                if (k >= n):
                    i += 1
                    break

                # Otherwise, add the remaining ring fragment and continue.
                ring = ring[k - 1:]
                fragments.append(
                    {"index": -1, "polygon": polygon, "ring": ring})
                ring[0] = normalizePoint(ring[0][1])
                i = 0
                n = len(ring)

            i += 1

# Now stitch the fragments back together into rings.


def stitchFragments(fragments):
    n = len(fragments)

    # To connect the fragments start-to-end, create a simple index by end.
    fragmentByStart = {}
    fragmentByEnd = {}

    # For each fragment…
    for i in range(n):
        fragment = fragments[i]
        start = tuple(fragment["ring"][0])
        end = tuple(fragment["ring"][-1])

        # If this fragment is closed, add it as a standalone ring.
        if (start[0] == end[0] and start[1] == end[1]):
            fragment["polygon"].append(fragment["ring"])
            fragments[i] = None
            continue

        fragment["index"] = i
        fragmentByStart[start] = fragmentByEnd[end] = fragment

    # For each open fragment…
    for i in range(n):
        fragment = fragments[i]
        if (fragment):
            start = tuple(fragment["ring"][0])
            end = tuple(fragment["ring"][-1])
            startFragment = fragmentByEnd[start]
            endFragment = fragmentByStart[end]

            del fragmentByStart[start]
            del fragmentByEnd[end]

            # If this fragment is closed, add it as a standalone ring.
            if (start[0] == end[0] and start[1] == end[1]):
                fragment["polygon"].append(fragment["ring"])
                continue

            if (startFragment):
                del fragmentByEnd[start]
                del fragmentByStart[tuple(startFragment["ring"][0])]
                startFragment["ring"].pop()  # drop the shared coordinate
                fragments[startFragment["index"]] = None

                startFragment["ring"].append(fragment["ring"])
                newring = startFragment["ring"] + fragment["ring"]
                fragment = {
                    "index": -1, "polygon": startFragment["polygon"], "ring": newring}

                if (startFragment == endFragment):
                    # Connect both ends to this single fragment to create a ring.
                    fragment["polygon"].append(fragment["ring"])
                else:
                    n += 1
                    fragment["index"] = n
                    fragmentByStart[fragment["ring"][0]] = fragment
                    fragmentByEnd[fragment["ring"][-1]] = fragment
                    fragments.append(fragment)

            elif (endFragment):
                del fragmentByStart[end]
                del fragmentByEnd[endFragment["ring"][-1]]
                fragment["ring"].pop()  # drop the shared coordinate
                n += 1

                fragment["ring"] = fragment["ring"] + endFragment["ring"]
                fragment = {
                    "index": n, "polygon": endFragment["polygon"], "ring": fragment["ring"]}
                fragments[endFragment["index"]] = None
                fragmentByStart[fragment["ring"][0]] = fragment
                fragmentByEnd[fragment["ring"][-1]] = fragment
                fragments.append(fragment)
            else:
                fragment["ring"].append(fragment["ring"][0])  # close ring
                fragment["polygon"].append(fragment["ring"])


def stitchFeature(input):
    output = {"type": "Feature",
              "geometry": stitchGeometry(input["geometry"])}
    if (input["id"] is not None):
        output["id"] = input["id"]
    if (input["bbox"] is not None):
        output["bbox"] = input["bbox"]
    if (input["properties"] is not None):
        output["properties"] = input["properties"]
    return output


def stitchGeometry(input):
    if input is None:
        return input

    if input["type"] == "GeometryCollection":
        output = {"type": "GeometryCollection", "geometries": [
            stitchGeometry(geom) for geom in input["geometries"]]}

    elif input["type"] == "Point":
        output = {"type": "Point",
                  "coordinates": clampPoint(input["coordinates"])}

    elif input["type"] in ("MultiPoint", "LineString"):
        output = {"type": input["type"],
                  "coordinates": clampPoint(input["coordinates"])}

    elif input["type"] == "MultiLineString":
        output = {"type": "MultiLineString", "coordinates": [
            clampPoints(coords) for coords in input["coordinates"]]}

    elif input["type"] == "Polygon":
        polygon = []
        fragments = []
        extractFragments(input["coordinates"], polygon, fragments)
        stitchFragments(fragments)
        output = {"type": "Polygon", "coordinates": polygon}

    elif input["type"] == "MultiPolygon":
        fragments = []
        i = 0
        n = len(input["coordinates"])
        polygons = [None]*n
        while i < n:
            polygons[i] = []
            extractFragments(input["coordinates"][i], polygons[i], fragments)
            i += 1
        stitchFragments(fragments)
        output = {"type": "MultiPolygon", "coordinates": [
            polygon for polygon in polygons if nonempty(polygon)]}

    else:
        return input

    if input.get("bbox", False):
        output["bbox"] = input["bbox"]
    return output


def geoStitch(input):
    if (input is None):
        return input

    if input["type"] == "Feature":
        return stitchFeature(input)

    elif input["type"] == "FeatureCollection":
        features = [stitchFeature(feature) for feature in input["features"]]
        output = {"type": "FeatureCollection", "features": features}
        if (input["bbox"] is not None):
            output["bbox"] = input["bbox"]
        return output

    return stitchGeometry(input)


if __name__ == "__main__":
    x = {
        "type": "Polygon",
        "bbox": [-180, -90, 180, 90],
        "coordinates": [
            [[-180, -80], [-90, -80], [0, -80], [90, -80], [180, -80], [180, -90],
                [90, -90], [0, -90], [-90, -90], [-180, -90], [-180, -80]]
        ]
    }

    print(geoStitch(x))

    y = {
        "type": "Feature",
        "id": "polygon",
        "bbox": [-180, -90, 180, 90],
        "properties": {"foo": 42},
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [[-180, -80], [-90, -80], [0, -80], [90, -80], [180, -80], [180, -90],
                 [90, -90], [0, -90], [-90, -90], [-180, -90], [-180, -80]]
            ]
        }
    }
    print(geoStitch(y))

    # applies an epsilon threshold to the poles and antimeridian
    epsilon = 1e-9
    z = {
        "type": "Polygon",
        "coordinates": [
            [[-180 + epsilon, -80], [-90, -80], [0, -80], [90, -80], [180 - epsilon, -80], [180 - epsilon, -90 + epsilon],
             [90, -90 + epsilon], [0, -90 + epsilon], [-90, -90 + epsilon], [-180 - epsilon, -90 + epsilon], [-180 - epsilon, -80]]
        ]
    }
    print(geoStitch(z))

    # surrounding the South pole with a cut along the antimeridian
    x = {
        "type": "Polygon",
        "coordinates": [
            [[-180, -80], [-90, -80], [0, -80], [90, -80], [180, -80], [180, -90],
                [90, -90], [0, -90], [-90, -90], [-180, -90], [-180, -80]]
        ]
    }
    print(geoStitch(x))

    # stitch(Polygon) with a hole across the antimeridian and cut along the antimeridian
    x = {
        "type": "Polygon",
        "coordinates": [
            [[-180, -60], [-180, -30], [-150, 0], [-180, 30], [-180, 60], [-60, 60], [60, 60],
             [180, 60], [180, 30], [150, 0], [180, -30], [180, -60], [60, -60], [-60, -60], [-180, -60]]
        ]
    }
    print(geoStitch(x))
