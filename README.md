# Geometry 2D Tiny Helpers

This library is a set of static methods to perform basic 2D geometry computations using points, lines and segments. They are useful for user interface programming, selection of graphics elements, data processing and computations.

- Distances : points, point / segment, point / line.
- Intersection : lines, segments.
- Perpendicular segments.
- Projection of a point on a line.
- Angles : of a segment, between two segments, half or complete clock.
- Segment interpolation.
- Point rotation around a point.
- Linear regression : slope, intercept and r computation.
- B-Spline interpolation.

# Distances

## Compute 2 Points Distance

Basic point to point distance.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/ecce3107-b437-4fa8-84d3-3a4a7c23053b)

The prototype of the method (using direct coords or GeometryPoint2D class) :

```
public static double ComputePointDistance(double x1, double y1, double x2, double y2);
public static double ComputePointDistance(GeometryPoint2D a, GeometryPoint2D b);
```

Sample :

```c#
Geometry2DTinyHelper.Distances.ComputePointDistance(new GeometryPoint2D(1, 0), new GeometryPoint2D(3, 0))
```

## Segment and point distance

Compute point to segment distance. A segment is a line section between two points.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/132d3176-0106-4b39-9407-48549d4b64e3)

The prototype of the method (using direct coords or GeometryPoint2D class) :

```
public static double ComputePointSegmentDistance(double x, double y, double x1, double y1, double x2, double y2);
public static double ComputePointSegmentDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

## Point and Line Distance

Compute point to line distance. A line is not limited.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/7be75cc1-bd72-4bd6-ba14-91b3986ff59a)

The prototype of the method :

```
public static double ComputePointLineDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

# Intersections

## Two Lines Intersection

Compute two line intersection. If they are parallel, a null point is returned.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/4366660e-978c-4608-a178-8947f5f7f264)

The prototype of the method :

```
public static GeometryPoint2D ComputeLineIntersection(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2);
```

## Two Segments Intersection

Compute two segment intersection. If they are parallel, a null point is assigned to 'intersection' parameter.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/253b34e5-b79e-4a8f-ace0-9c71f971b43d)

The prototype of the method :

```
static bool ComputeSegmentIntersection(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2, out GeometryPoint2D intersection);
```

## ComputePerpendicularLine

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e479357b-a3c3-4e1a-9239-cac2f3bcbb0b)

## ProjectPointOnLine

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/67768097-6801-4824-870a-719abe0043d5)

## IsPointOnSegment

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/1e9f71d1-70ab-4650-9ab2-7e80eb0dad52)

# Angles

## ComputeSegmentAngle

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/97bbefc7-f884-458f-9f3e-feaa2f8b81ec)

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/68b4d80d-c9d5-45e3-9671-2977ca59c202)

## ConvertSlopeToAngle

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/97a29496-612c-4dc8-a5af-e55353415ffe)

## ComputeSegmentPositiveAngle

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/8858523f-f36c-4f1c-945c-ea1f5d5551c2)

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/0b2eef83-c990-4a44-8be8-143e11062dd0)

## ComputeSlopeToPositiveAngle

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/ec9b6d55-6dc2-43de-a4c9-eabc0850b4a8)

## ComputeSegmentAngleDelta

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e501bfa0-ea9e-4416-8ad4-db28aac9a368)

## Angle origin conversion

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/0a2c1fb2-e636-4469-a260-772db8bb7a5b)

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/f3a44b4c-d2fd-4368-89d4-f73155ac59a4)

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/991c29d9-5365-4cf4-88a8-1348be19c73d)

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/293b592a-02aa-4c2c-84c9-6614ecb044f8)

# Interpolations

## Interpolate

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/90352c5b-6525-43b9-8f9c-0d7007f75623)

## RotatePoint

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/8a45b7ce-0580-4564-8849-7e0f45b35059)

## Interpolate

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e3d2571e-64e5-4711-b3f4-a6f3ba25ed01)

## ExponentialDecreasingInterpolation

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/b36df4c4-5557-4a1d-b2cb-aa5ea476f886)

## Normalize

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/665cde81-9581-4ddf-9840-e9cf7fcf0fea)

## Sample

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/99610a7e-6fd5-47e6-b83e-36f520f754c3)

## GetGamma

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/d5018d49-8c4f-494e-bcd9-3e5e04f84505)

## Sample

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e3e3ed11-773a-49a2-9099-83ef9f6a7ac8)

# LinearRegressions

## ComputeCorrelationCoefficient

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/b5499068-9637-436f-a9e3-c0e8fc865c79)

## ComputeLinearRegression

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/593e6aea-0c5e-4571-b14c-54c2d7dd5628)

# BSpline

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/f2c714e2-6989-4fc8-8a80-91faa4c0f638)

# Surfaces

## ComputeTriangleArea

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/d8d22e7a-a216-4689-9e93-72aa591cd901)

## IsPointInTriangle

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/7ffd8da7-e94e-42ec-a355-94bc0a151207)



