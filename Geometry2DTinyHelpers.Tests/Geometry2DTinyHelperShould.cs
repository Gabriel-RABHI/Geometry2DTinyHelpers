using static System.Formats.Asn1.AsnWriter;

namespace Geometry2DTinyHelpers.Tests
{
    public class Geometry2DTinyHelperShould
    {
        [Fact(DisplayName = "Compute points distances")]
        public void ComputePointDistance()
        {
            Assert.Equal(2, Geometry2DTinyHelper.Distances.PointDistance(1, 0, 3, 0));
            Assert.Equal(2, Geometry2DTinyHelper.Distances.PointDistance(new GeometryPoint2D(1, 0), new GeometryPoint2D(3, 0)));
        }

        [Fact(DisplayName = "Compute points to segment / line distances")]
        public void ComputePointSegmentAndLineDistance()
        {
            Assert.Equal(2, Geometry2DTinyHelper.Distances.PointSegmentDistance(2, 2, 1, 0, 3, 0));
            Assert.Equal(2, Geometry2DTinyHelper.Distances.PointSegmentDistance(new GeometryPoint2D(2, 2), new GeometryPoint2D(1, 0), new GeometryPoint2D(3, 0)));
            Assert.Equal(1, Geometry2DTinyHelper.Distances.PointSegmentDistance(new GeometryPoint2D(4, 0), new GeometryPoint2D(1, 0), new GeometryPoint2D(3, 0)));
            Assert.Equal(2, Geometry2DTinyHelper.Distances.PointLineDistance(new GeometryPoint2D(0, 2), new GeometryPoint2D(1, 0), new GeometryPoint2D(3, 0)));
        }

        [Fact(DisplayName = "Line / segment intersection")]
        public void ComputeIntersections()
        {
            Assert.Equal(new GeometryPoint2D(2, 2), Geometry2DTinyHelper.Intersections.LineIntersection(
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(3, 3),
                new GeometryPoint2D(1, 3),
                new GeometryPoint2D(3, 1)
            ));

            var result = new GeometryPoint2D();
            Assert.True(Geometry2DTinyHelper.Intersections.SegmentIntersection(
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(3, 3),
                new GeometryPoint2D(1, 3),
                new GeometryPoint2D(3, 1),
                out result
            ));
            Assert.Equal(new GeometryPoint2D(2, 2), result);
        }

        [Fact(DisplayName = "Perpendiculare")]
        public void PerpendicularLineAndProjectPoint()
        {
            var seg = Geometry2DTinyHelper.Intersections.GetPerpendicularLine(
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(3, 1), 2
            );
            Assert.Equal(new GeometryPoint2D(1, 1), seg[0]);
            Assert.Equal(new GeometryPoint2D(1, 3), seg[1]);
        }

        [Fact(DisplayName = "Project point on line")]
        public void ProjectPointOnLine()
        {
            var p = Geometry2DTinyHelper.Intersections.ProjectPointOnLine(
                new GeometryPoint2D(1, 3),
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(3, 3)
            );
            Assert.Equal(new GeometryPoint2D(2, 2), p);
        }

        [Fact(DisplayName = "Segment Angles")]
        public void SegmentAngle()
        {
            var a1 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, 1)
            );
            var a2 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 1)
            );
            var a3 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 0)
            );
            var a4 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, -1)
            );
            var a5 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, -1)
            );
            var a6 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, -1)
            );
            var a7 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, 0)
            );
            var a8 = Geometry2DTinyHelper.Angles.SegmentAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, 1)
            );
            Assert.Equal(-90, a1);
            Assert.Equal(-45, a2);
            Assert.Equal(0, a3);
            Assert.Equal(45, a4);
            Assert.Equal(90, a5);
            Assert.Equal(135, a6);
            Assert.Equal(180, a7);
            Assert.Equal(-135, a8);
        }

        [Fact(DisplayName = "Segment Positive Angles")]
        public void SegmentPositiveAngle()
        {
            var a1 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, 1)
            );
            var a2 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 1)
            );
            var a3 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 0)
            );
            var a4 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, -1)
            );
            var a5 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, -1)
            );
            var a6 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, -1)
            );
            var a7 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, 0)
            );
            var a8 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, 1)
            );
            Assert.Equal(270, a1);
            Assert.Equal(315, a2);
            Assert.Equal(0, a3);
            Assert.Equal(45, a4);
            Assert.Equal(90, a5);
            Assert.Equal(135, a6);
            Assert.Equal(180, a7);
            Assert.Equal(225, a8);
        }

        [Fact(DisplayName = "Two Segments Angle")]
        public void SegmentAngleDelta()
        {
            var a1 = Geometry2DTinyHelper.Angles.SegmentAngleDelta(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, 1),
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 1)
            );
            var a2 = Geometry2DTinyHelper.Angles.SegmentAngleDelta(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, 1),
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, -1)
            );
            var a3 = Geometry2DTinyHelper.Angles.SegmentAngleDelta(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(-1, 0),
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 0)
            );
            Assert.Equal(45, a1);
            Assert.Equal(135, a2);
            Assert.Equal(180, a3);
        }

        [Fact(DisplayName = "Slope Angle")]
        public void SlopeAngle()
        {
            var a1 = Geometry2DTinyHelper.Angles.SlopToAngle(1);
            var a2 = Geometry2DTinyHelper.Angles.SlopToAngle(-1);
            var a3 = Geometry2DTinyHelper.Angles.SlopToPositiveAngle(1);
            var a4 = Geometry2DTinyHelper.Angles.SlopToPositiveAngle(-1);
            Assert.Equal(-45, a1);
            Assert.Equal(45, a2);
            Assert.Equal(315, a3);
            Assert.Equal(45, a4);
        }

        [Fact(DisplayName = "Angle offset")]
        public void AngleOffseting()
        {
            var a1 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(1, 0)
            );
            var a2 = Geometry2DTinyHelper.Angles.SegmentPositiveAngle(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(0, -1)
            );
            Assert.Equal(0, a1);
            Assert.Equal(90, a2);

            Assert.Equal(90, Geometry2DTinyHelper.Angles.ToNorth(a1));
            Assert.Equal(180, Geometry2DTinyHelper.Angles.ToNorth(a2));

            Assert.Equal(180, Geometry2DTinyHelper.Angles.ToWest(a1));
            Assert.Equal(270, Geometry2DTinyHelper.Angles.ToWest(a2));

            Assert.Equal(270, Geometry2DTinyHelper.Angles.ToSouth(a1));
            Assert.Equal(0, Geometry2DTinyHelper.Angles.ToSouth(a2));
        }

        [Fact(DisplayName = "Interpolation / extrapolation")]
        public void Interpolations()
        {
            var a1 = Geometry2DTinyHelper.Interpolations.Interpolate(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(2, 2),
                0.5
            );
            var a2 = Geometry2DTinyHelper.Interpolations.Interpolate(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(2, 2),
                1.5
            );
            Assert.Equal(new GeometryPoint2D(1, 1), a1);
            Assert.Equal(new GeometryPoint2D(3, 3), a2);

            var y1 = Geometry2DTinyHelper.Interpolations.Extrapolation(
               new GeometryPoint2D(0, 0),
               new GeometryPoint2D(2, 2),
               1
           );
            var y2 = Geometry2DTinyHelper.Interpolations.Extrapolation(
                new GeometryPoint2D(0, 0),
                new GeometryPoint2D(2, 2),
                3
            );
            Assert.Equal(1, y1);
            Assert.Equal(3, y2);

            a1 = Geometry2DTinyHelper.Interpolations.Sample(
               1, 2, 2
            );
            a2 = Geometry2DTinyHelper.Interpolations.Sample(
               1, 2, 5
            );
            Assert.Equal(new GeometryPoint2D(2, 4), a1);
            Assert.Equal(new GeometryPoint2D(5, 7), a2);
        }

        [Fact(DisplayName = "Rotation")]
        public void Rotation()
        {
            var a1 = Geometry2DTinyHelper.Interpolations.RotatePoint(
                new GeometryPoint2D(2, 2),
                new GeometryPoint2D(0, 0),
                90
            );
            var a2 = Geometry2DTinyHelper.Interpolations.RotatePoint(
                new GeometryPoint2D(2, 2),
                new GeometryPoint2D(0, 0),
                180
            );
            Assert.Equal(new GeometryPoint2D(2, -2), new GeometryPoint2D(Math.Round(a1.X), Math.Round(a1.Y)));
            Assert.Equal(new GeometryPoint2D(-2, -2), new GeometryPoint2D(Math.Round(a2.X), Math.Round(a2.Y)));
        }

        [Fact(DisplayName = "Linear Regressions Aligned")]
        public void LinearRegressions()
        {
            GeometryPoint2D[] serie = new GeometryPoint2D[]
            {
                new GeometryPoint2D(2, 4),
                new GeometryPoint2D(4, 5),
                new GeometryPoint2D(6, 6),
                new GeometryPoint2D(8, 7)
            };

            double rSquared, yIntercept, slope;
            Geometry2DTinyHelper.LinearRegressions.LinearRegression(serie, out rSquared, out yIntercept, out slope);

            Assert.Equal(1, rSquared);
            Assert.Equal(3, yIntercept);
            Assert.Equal(0.5, slope);

            var corr = Geometry2DTinyHelper.LinearRegressions.LinearCorrelationCoefficient(serie);
            Assert.Equal(1, corr);
        }

        [Fact(DisplayName = "Linear Regressions Non-aligned")]
        public void LinearRegressionsNonAligned()
        {
            GeometryPoint2D[] serie = new GeometryPoint2D[]
            {
                new GeometryPoint2D(2, 4-1),
                new GeometryPoint2D(4, 5+1),
                new GeometryPoint2D(6, 6-1),
                new GeometryPoint2D(8, 7+1)
            };

            Geometry2DTinyHelper.LinearRegressions.LinearRegression(serie, out var rSquared, out var yIntercept, out var slope);

            Assert.True(rSquared < 0.76 && rSquared >= 0.75);
            Assert.Equal(2, yIntercept);
            Assert.Equal(0.7, slope);

            var corr = Geometry2DTinyHelper.LinearRegressions.LinearCorrelationCoefficient(serie);
            Assert.True(corr < 0.76 && corr >= 0.75);
        }

        [Fact(DisplayName = "B-Spline")]
        public void BSpline()
        {
            GeometryPoint2D[] serie = new GeometryPoint2D[]
            {
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(4, 3),
                new GeometryPoint2D(5, 1)
            };

            var r = Geometry2DTinyHelper.BSpline.Interpolate1D(serie, 10);
            Assert.Equal(10, r.Length);
            Assert.Equal(new GeometryPoint2D(1, 1), r[0]);
            Assert.Equal(new GeometryPoint2D(5, 1), r[9]);
        }

        [Fact(DisplayName = "Surfaces")]
        public void Surfaces()
        {
            GeometryPoint2D[] triangle = new GeometryPoint2D[]
            {
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(3, 3),
                new GeometryPoint2D(4, 1)
            };

            var r = Geometry2DTinyHelper.Surfaces.Surface(triangle[0], triangle[1], triangle[2]);
            Assert.Equal(3, Math.Round(r));
            Assert.True(Geometry2DTinyHelper.Surfaces.PointInTriangle(new GeometryPoint2D(3, 2), triangle[0], triangle[1], triangle[2]));
            Assert.False(Geometry2DTinyHelper.Surfaces.PointInTriangle(new GeometryPoint2D(1, 2.5), triangle[0], triangle[1], triangle[2]));
        }
    }
}