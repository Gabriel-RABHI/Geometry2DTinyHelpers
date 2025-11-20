![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/99acb6f3-fd69-4721-9d85-011d158442ae)

# Geometry 2D Tiny Helpers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![.NET](https://img.shields.io/badge/.NET-5.0+-512BD4.svg)](https://dotnet.microsoft.com/)
[![C#](https://img.shields.io/badge/C%23-9.0+-239120.svg)](https://docs.microsoft.com/en-us/dotnet/csharp/)

A lightweight, easy-to-use C# library providing static methods for common 2D geometry computations. Perfect for developers who need reliable geometric calculations without complex mathematical expertise.

## Features

- **Distances**: Point-to-point, point-to-segment, and point-to-line calculations
- **Intersections**: Line and segment intersection detection with perpendicular line computation
- **Projections**: Project points onto lines
- **Angles**: Segment angles, angle conversions, and delta calculations with multiple origin systems
- **Interpolations**: Linear and exponential interpolation, point rotation, and normalization
- **Linear Regression**: Slope, intercept, and correlation coefficient computation
- **B-Spline**: Smooth curve interpolation through control points
- **Surface Calculations**: Triangle area and point-in-triangle tests

### Manual Installation
Clone the repository and add a reference to the project:
```bash
git clone https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers.git
```

## 🚀 Quick Start

```csharp
using Geometry2DTinyHelpers;

// Calculate distance between two points
double distance = Geometry2DTinyHelper.Distances.ComputePointDistance(
    new GeometryPoint2D(1, 0), 
    new GeometryPoint2D(3, 0)
);

// Find intersection of two lines
GeometryPoint2D intersection = Geometry2DTinyHelper.Intersections.ComputeLineIntersection(
    new GeometryPoint2D(0, 0), new GeometryPoint2D(2, 2),
    new GeometryPoint2D(2, 0), new GeometryPoint2D(0, 2)
);

// Compute segment angle
double angle = Geometry2DTinyHelper.Angles.ComputeSegmentAngle(
    new GeometryPoint2D(0, 0), 
    new GeometryPoint2D(1, 1)
);
```

## API Organization

All methods are organized into nested static classes:

- `Geometry2DTinyHelper.Distances` - Distance calculations
- `Geometry2DTinyHelper.Intersections` - Intersection and projection operations
- `Geometry2DTinyHelper.Angles` - Angle computations and conversions
- `Geometry2DTinyHelper.Interpolations` - Interpolation and rotation functions
- `Geometry2DTinyHelper.LinearRegressions` - Statistical regression analysis
- `Geometry2DTinyHelper.BSpline` - B-Spline curve generation
- `Geometry2DTinyHelper.Surfaces` - Area and containment tests

## Documentation

### Distances

#### Compute Point-to-Point Distance

Calculate the Euclidean distance between two points.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/ecce3107-b437-4fa8-84d3-3a4a7c23053b)

**Method Signature:**
```csharp
double ComputePointDistance(double x1, double y1, double x2, double y2);
double ComputePointDistance(GeometryPoint2D a, GeometryPoint2D b);
```

**Example:**
```csharp
double distance = Geometry2DTinyHelper.Distances.ComputePointDistance(
    new GeometryPoint2D(1, 0), 
    new GeometryPoint2D(3, 0)
); // Returns: 2.0
```

#### Point-to-Segment Distance

Compute the shortest distance from a point to a line segment.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/132d3176-0106-4b39-9407-48549d4b64e3)

**Method Signature:**
```csharp
double ComputePointSegmentDistance(double x, double y, double x1, double y1, double x2, double y2);
double ComputePointSegmentDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

#### Point-to-Line Distance

Compute the perpendicular distance from a point to an infinite line.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/7be75cc1-bd72-4bd6-ba14-91b3986ff59a)

**Method Signature:**
```csharp
double ComputePointLineDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

### Intersections

#### Two Lines Intersection

Find the intersection point of two infinite lines. Returns `null` if lines are parallel.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/4366660e-978c-4608-a178-8947f5f7f264)

**Method Signature:**
```csharp
GeometryPoint2D ComputeLineIntersection(
    GeometryPoint2D a1, GeometryPoint2D b1,
    GeometryPoint2D a2, GeometryPoint2D b2
);
```

#### Two Segments Intersection

Find the intersection point of two line segments. Returns `true` if intersection is within both segments.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/253b34e5-b79e-4a8f-ace0-9c71f971b43d)

**Method Signature:**
```csharp
bool ComputeSegmentIntersection(
    GeometryPoint2D a1, GeometryPoint2D b1,
    GeometryPoint2D a2, GeometryPoint2D b2,
    out GeometryPoint2D intersection
);
```

#### Compute Perpendicular Line

Generate a perpendicular segment of specified length from a point.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e479357b-a3c3-4e1a-9239-cac2f3bcbb0b)

**Method Signature:**
```csharp
GeometryPoint2D[] ComputePerpendicularLine(GeometryPoint2D a, GeometryPoint2D b, float length);
```

#### Project Point onto Line

Project a point perpendicularly onto a line.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/67768097-6801-4824-870a-719abe0043d5)

**Method Signature:**
```csharp
GeometryPoint2D ProjectPointOnLine(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b);
```

#### Check Point on Segment

Verify if a point lies on a segment within a tolerance.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/1e9f71d1-70ab-4650-9ab2-7e80eb0dad52)

**Method Signature:**
```csharp
bool IsPointOnSegment(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b, double epsilon = 0.0001);
```

### Angles

#### Compute Segment Angle

Calculate the inclination angle of a segment. Horizontal lines are 0°, with East as origin.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/97bbefc7-f884-458f-9f3e-feaa2f8b81ec)

**Method Signature:**
```csharp
double ComputeSegmentAngle(double x1, double y1, double x2, double y2);
double ComputeSegmentAngle(GeometryPoint2D a, GeometryPoint2D b);
```

**Angle Convention:** Upper segments have negative angles (0° to -180°), lower segments have positive angles (0° to 180°).

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/68b4d80d-c9d5-45e3-9671-2977ca59c202)

#### Convert Slope to Angle

Convert a line's slope coefficient to angle in degrees.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/97a29496-612c-4dc8-a5af-e55353415ffe)

**Method Signature:**
```csharp
double ConvertSlopeToAngle(double slope);
```

#### Compute Positive Segment Angle

Calculate segment angle as positive values (0° to 360°), clockwise from East.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/8858523f-f36c-4f1c-945c-ea1f5d5551c2)

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/0b2eef83-c990-4a44-8be8-143e11062dd0)

**Method Signature:**
```csharp
double ComputeSegmentPositiveAngle(GeometryPoint2D a, GeometryPoint2D b);
```

#### Compute Slope to Positive Angle

Convert slope to positive angle (0° to 360°), clockwise.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/ec9b6d55-6dc2-43de-a4c9-eabc0850b4a8)

**Method Signature:**
```csharp
double ComputeSlopeToPositiveAngle(double slope);
```

#### Compute Segment Angle Delta

Calculate the angle between two segments (always between 0° and 180°).

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e501bfa0-ea9e-4416-8ad4-db28aac9a368)

**Method Signature:**
```csharp
double ComputeSegmentAngleDelta(
    GeometryPoint2D a1, GeometryPoint2D b1, 
    GeometryPoint2D a2, GeometryPoint2D b2
);
```

#### Angle Origin Conversions

Convert angles between different coordinate system origins (East, North, West, South).

**East to North:**
![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/0a2c1fb2-e636-4469-a260-772db8bb7a5b)

**East to West:**
![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/293b592a-02aa-4c2c-84c9-6614ecb044f8)

**East to South:**
![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/991c29d9-5365-4cf4-88a8-1348be19c73d)

**Method Signatures:**
```csharp
double ToNorth(double angle);
double ToWest(double angle);
double ToSouth(double angle);
```

### Interpolations

#### Interpolate Segment

Interpolate or extrapolate along a segment using a factor (0-1 for interpolation).

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/90352c5b-6525-43b9-8f9c-0d7007f75623)

**Method Signature:**
```csharp
GeometryPoint2D Interpolate(GeometryPoint2D a, GeometryPoint2D b, double factor);
```

#### Rotate Point

Rotate a point around another point by a specified angle.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/8a45b7ce-0580-4564-8849-7e0f45b35059)

**Method Signature:**
```csharp
GeometryPoint2D RotatePoint(GeometryPoint2D pointToRotate, GeometryPoint2D centerPoint, double angleInDegrees);
```

#### Interpolate Values

Interpolate between two numeric values using a factor.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e3d2571e-64e5-4711-b3f4-a6f3ba25ed01)

**Method Signature:**
```csharp
double Interpolate(double a, double b, double factor);
```

#### Exponential Decreasing Interpolation

Apply exponential decay interpolation between two values.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/b36df4c4-5557-4a1d-b2cb-aa5ea476f886)

**Method Signature:**
```csharp
double ExponentialDecreasingInterpolation(double a, double b, double factor);
```

#### Normalize

Calculate the normalized position of a value between two bounds.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/665cde81-9581-4dff-9840-e9cf7fcf0fea)

**Method Signature:**
```csharp
double Normalize(double a, double b, double position);
```

#### Sample Y from X

Calculate the Y coordinate for a given X on a line defined by two points.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/99610a7e-6fd5-47e6-b83e-36f520f754c3)

**Method Signature:**
```csharp
double Sample(GeometryPoint2D a, GeometryPoint2D b, double x);
```

#### Gamma Function

Apply gamma correction to values in the 0-1 range.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/d5018d49-8c4f-494e-bcd9-3e5e04f84505)

**Method Signature:**
```csharp
double GetGamma(double x, double curvature = 1);
```

#### Sample Point from Slope

Calculate a point on a line defined by slope and Y-intercept.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/e3e3ed11-773a-49a2-9099-83ef9f6a7ac8)

**Method Signature:**
```csharp
GeometryPoint2D Sample(double slope, double yIntercept, double x);
```

### Linear Regression

Linear regression predicts unknown values using known, related data by modeling relationships mathematically. It's used in economics, finance, and data analysis for trend forecasting.

**Use Cases:**
- Generate a linear representation (slope and Y-intercept) of a point set
- Calculate correlation coefficient to measure data alignment
- Predict the next probable value in a sequence

#### Compute Correlation Coefficient

Calculate the correlation coefficient (r) of a point series. Perfect alignment = 1, random dispersion ≈ 0.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/b5499068-9637-436f-a9e3-c0e8fc865c79)

**Method Signature:**
```csharp
double ComputeCorrelationCoefficient(GeometryPoint2D[] pts);
```

#### Compute Linear Regression

Calculate regression line parameters (slope, Y-intercept) and correlation coefficient.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/593e6aea-0c5e-4571-b14c-54c2d7dd5628)

**Method Signature:**
```csharp
void ComputeLinearRegression(
    GeometryPoint2D[] pts,
    out double rSquared,
    out double yIntercept,
    out double slope
);
```

### B-Spline Interpolation

#### Interpolate 1D

Generate a smooth poly-line curve through control points using B-Spline interpolation.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/f2c714e2-6989-4fc8-8a80-91faa4c0f638)

**Method Signatures:**
```csharp
GeometryPoint2D[] Interpolate1D(GeometryPoint2D[] pts, int count);
(double[] xs, double[] ys) Interpolate1D(double[] xs, double[] ys, int count);
```

### Surface Calculations

#### Compute Triangle Area

Calculate the area of a triangle defined by three points.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/d8d22e7a-a216-4689-9e93-72aa591cd901)

**Method Signature:**
```csharp
double ComputeTriangleArea(GeometryPoint2D a, GeometryPoint2D b, GeometryPoint2D c);
```

#### Check Point in Triangle

Determine if a point lies inside a triangle.

![image](https://github.com/Gabriel-RABHI/Geometry2DTinyHelpers/assets/8116286/7ffd8da7-e94e-42ec-a355-94bc0a151207)

**Method Signature:**
```csharp
bool IsPointInTriangle(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b, GeometryPoint2D c);
```

## 🧪 Testing

The library includes comprehensive unit tests in the `Geometry2DTinyHelpers.Tests` project.

```bash
dotnet test
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.

## 👤 Author

**Gabriel RABHI**
- GitHub: [@Gabriel-RABHI](https://github.com/Gabriel-RABHI)

## ⭐ Support

If you find this library helpful, please consider giving it a star on GitHub!

## 📊 Project Stats

![GitHub stars](https://img.shields.io/github/stars/Gabriel-RABHI/Geometry2DTinyHelpers?style=social)
![GitHub forks](https://img.shields.io/github/forks/Gabriel-RABHI/Geometry2DTinyHelpers?style=social)

---

**Made with ❤️ for developers who need simple, reliable 2D geometry calculations**
