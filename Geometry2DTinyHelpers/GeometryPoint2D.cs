using System;

namespace Geometry2DTinyHelpers
{
    public class GeometryPoint2D : IEquatable<GeometryPoint2D>
    {
        public double X { get; set; }

        public double Y { get; set; }

        public GeometryPoint2D()
            => X = Y = 0;

        public GeometryPoint2D(GeometryPoint2D source)
        {
            X = source.X;
            Y = source.Y;
        }

        public GeometryPoint2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        public bool Equals(GeometryPoint2D? other)
            => other != null && other.X == X && other.Y == Y;
    }
}
