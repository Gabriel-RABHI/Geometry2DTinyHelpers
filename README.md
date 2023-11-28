![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/99acb6f3-fd69-4721-9d85-011d158442ae)

# Geometry 2D Tiny Helpers

This library is a set of static methods to perform basic 2D geometry computations using points, lines and segments. This library is for the one who haven't mathematic skill. This methods are useful for user interface programming, selection of graphics elements, data processing and computations.

- Distances : points, point / segment, point / line.
- Intersection : lines, segments.
- Perpendicular segments.
- Projection of a point on a line.
- Angles : of a segment, between two segments, half or complete clock.
- Segment interpolation.
- Point rotation around a point.
- Linear regression : slope, intercept and r computation.
- B-Spline interpolation.
- Area of a triangle and inner point test.

## Distances

### Compute 2 Points Distance

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

### Segment and point distance

Compute point to segment distance. A segment is a line section between two points.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/132d3176-0106-4b39-9407-48549d4b64e3)

The prototype of the method (using direct coords or GeometryPoint2D class) :

```
public static double ComputePointSegmentDistance(double x, double y, double x1, double y1, double x2, double y2);
public static double ComputePointSegmentDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

### Point and Line Distance

Compute point to line distance. A line is not limited.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/7be75cc1-bd72-4bd6-ba14-91b3986ff59a)

The prototype of the method :

```
public static double ComputePointLineDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

## Intersections

### Two Lines Intersection

Compute two line intersection. If they are parallel, a null point is returned.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/4366660e-978c-4608-a178-8947f5f7f264)

The prototype of the method :

```
public static GeometryPoint2D ComputeLineIntersection(
    GeometryPoint2D a1,
    GeometryPoint2D b1,
    GeometryPoint2D a2,
    GeometryPoint2D b2);
```

### Two Segments Intersection

Compute two segment intersection. If they are parallel, a null point is assigned to 'intersection' parameter.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/253b34e5-b79e-4a8f-ace0-9c71f971b43d)

The prototype of the method :

```
static bool ComputeSegmentIntersection(
    GeometryPoint2D a1,
    GeometryPoint2D b1,
    GeometryPoint2D a2,
    GeometryPoint2D b2,
    out GeometryPoint2D intersection);
```

### ComputePerpendicularLine

Compute a segment points of length 'length', perpendicular to the given segment and starting from 'a' point.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e479357b-a3c3-4e1a-9239-cac2f3bcbb0b)

The prototype of the method :

```
public static GeometryPoint2D[] ComputePerpendicularLine(GeometryPoint2D a, GeometryPoint2D b, float length);
```

### ProjectPointOnLine

Compute a point projected on a line.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/67768097-6801-4824-870a-719abe0043d5)

The prototype of the method :

```
public static GeometryPoint2D ProjectPointOnLine(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

### IsPointOnSegment

Check if a point is on a segment.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/1e9f71d1-70ab-4650-9ab2-7e80eb0dad52)

The prototype of the method :

```
static bool IsPointOnSegment(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b, double epsilon = 0.0001);
```

## Angles

### ComputeSegmentAngle

Compute a single segment angle, or inclination. An horizontale line is 0 degree (East origin). Upper angles are negative, lower angles are positive.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/97bbefc7-f884-458f-9f3e-feaa2f8b81ec)

The prototype of the method :

```
public static double ComputeSegmentAngle(double x1, double y1, double x2, double y2);
public static double ComputeSegmentAngle(GeometryPoint2D a, GeometryPoint2D b);
```

Angle in degree is computed on a non clock based scale, where upper oriented segment are negative angle, from 0° to near -180°, and lower oriented segment are positive angle from 0° to 180°.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/68b4d80d-c9d5-45e3-9671-2977ca59c202)

### ConvertSlopeToAngle

Convert a slope to angle in degree.  A line can be described as two points, or as a slope and an Y vertical offset.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/97a29496-612c-4dc8-a5af-e55353415ffe)

The prototype of the method :

```
public static double ConvertSlopeToAngle(double slop);
```

### ComputeSegmentPositiveAngle

Compute a single segment angle, or inclination. An horizontale line is 0 degree (East oriented). Angles are positive, from 0° to near 360°, and clockwise.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/8858523f-f36c-4f1c-945c-ea1f5d5551c2)

A 'positive angle' in degree have a range from 0° to near 360°. In every cases, the angle is still represented as positive value.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/0b2eef83-c990-4a44-8be8-143e11062dd0)

The prototype of the method :

```
 public static double ComputeSegmentPositiveAngle(GeometryPoint2D a, GeometryPoint2D b);
```

### ComputeSlopeToPositiveAngle

Convert à line slope to angle in degree, angles are positive, from 0° to near 360°, and clockwise.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/ec9b6d55-6dc2-43de-a4c9-eabc0850b4a8)

The prototype of the method :

```
public static double ComputeSlopeToPositiveAngle(double slop);
```

### ComputeSegmentAngleDelta

Compute two segment angle. Segments do not have to share any points. The result angle is all time in between 0° and 180°.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e501bfa0-ea9e-4416-8ad4-db28aac9a368)

The prototype of the method :

```
public static double ComputeSegmentAngleDelta(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2)
```

### Angle origin conversion

In above computations, angle origin is on the right side of the circle : it is an East origin. If you want the angle to have another origin, you can convert it easily.

Convert an East angle to North origin one.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/0a2c1fb2-e636-4469-a260-772db8bb7a5b)

The prototype of the method :

```
public static double ToNorth(double angle);
```

Convert an East angle to West origin one.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/991c29d9-5365-4cf4-88a8-1348be19c73d)

The prototype of the method :

```
public static double ToWest(double angle);
```

Convert an East angle to South origin one.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/293b592a-02aa-4c2c-84c9-6614ecb044f8)

The prototype of the method :

```
public static double ToSouth(double angle);
```

## Interpolations

### Interpolate

Interpolate / extrapolate a segment.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/90352c5b-6525-43b9-8f9c-0d7007f75623)

The prototype of the method :

```
public static GeometryPoint2D Interpolate(GeometryPoint2D a, GeometryPoint2D b, double factor);
```

### RotatePoint

Compute à rotated point around an another point.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/8a45b7ce-0580-4564-8849-7e0f45b35059)

The prototype of the method :

```
public static GeometryPoint2D RotatePoint(GeometryPoint2D pointToRotate, GeometryPoint2D centerPoint, double angleInDegrees);
```

### Interpolate

Compute à value in between two values, if in between 0 ans 1 - or over / under if not.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e3d2571e-64e5-4711-b3f4-a6f3ba25ed01)

The prototype of the method :

```
private static double Interpolate(double a, double b, double factor);
```

### ExponentialDecreasingInterpolation

Compute a decreased exponential interpolation value between two values.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/b36df4c4-5557-4a1d-b2cb-aa5ea476f886)

The prototype of the method :

```
public static double ExponentialDecreasingInterpolation(double a, double b, double factor);
```

### Normalize

Compute the position factor of a value in between (or not) two values.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/665cde81-9581-4ddf-9840-e9cf7fcf0fea)

The prototype of the method :

```
public static double Normalize(double a, double b, double position);
```

### Sample

Compute for a line the Y vertical coords for the given X horizontal position.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/99610a7e-6fd5-47e6-b83e-36f520f754c3)

The prototype of the method :

```
public static double Sample(GeometryPoint2D a, GeometryPoint2D b, double x);
```

### GetGamma

Sample à curved function between 0 and 1, returning a value between 0 ans 1.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/d5018d49-8c4f-494e-bcd9-3e5e04f84505)

The prototype of the method :

```
public static double GetGamma(double x, double curvatur = 1);
```

### Sample

Compute the point on a line for the given X horizontal position.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e3e3ed11-773a-49a2-9099-83ef9f6a7ac8)

The prototype of the method :

```
public static GeometryPoint2D Sample(double slope, double yIntercept, double x);
```

## LinearRegressions

Linear regression is **a data analysis technique that predicts the value of unknown data by using another related and known data value**. It mathematically models the unknown or dependent variable and the known or independent variable as a linear equation. It seeks the optimal line that minimizes the sum of squared differences between predicted and actual values.

Applied in various domains like economics and finance, this method analyzes and forecasts data trends. It can extend to multiple linear regression involving several independent variables and logistic regression, suitable for binary classification problems.

Practically, there is three use cases :

- Get a linear representation of a points set : a slope and an y offset, witch define a line.
- Get a correlation or, level of alignment of the points (to get a randomness coefficient).
- Compute the next most probable value for a point set, using the first point slope and offset.

### ComputeCorrelationCoefficient

Compute the correlation coefficient (or 'r') of a point serie. If points are perfectly aligned, result correlation is 1. If they are more randomly dispersed, correlation decrease to near 0.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/b5499068-9637-436f-a9e3-c0e8fc865c79)

The prototype of the method :

```
public static double ComputeCorrelationCoefficient(GeometryPoint2D[] pts);
```

### ComputeLinearRegression

Compute the regression line parameters (slope and Y offset) and correlation coeficient of a point serie.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/593e6aea-0c5e-4571-b14c-54c2d7dd5628)

The prototype of the method :

```
public static void ComputeLinearRegression(
    GeometryPoint2D[] pts,
    out double rSquared,
    out double yIntercept,
    out double slope);
```

## BSpline

### Interpolate1D

Return a the point array that form the poly-line that represent the spline curve.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/f2c714e2-6989-4fc8-8a80-91faa4c0f638)

The prototype of the method :

```
public static GeometryPoint2D[] Interpolate1D(GeometryPoint2D[] pts, int count);
public static (double[] xs, double[] ys) Interpolate1D(double[] xs, double[] ys, int count);
```

## Surfaces

### ComputeTriangleArea

Compute the surface of a triangle.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/d8d22e7a-a216-4689-9e93-72aa591cd901)

The prototype of the method :

```
public static double ComputeTriangleArea(GeometryPoint2D a, GeometryPoint2D b, GeometryPoint2D c);
```

### IsPointInTriangle

Check if a point is in a triangle.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/7ffd8da7-e94e-42ec-a355-94bc0a151207)

The prototype of the method :

```
public static bool IsPointInTriangle(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b, GeometryPoint2D c)
```



