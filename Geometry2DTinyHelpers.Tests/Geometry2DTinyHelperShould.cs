namespace Geometry2DTinyHelpers.Tests
{
    public class Geometry2DTinyHelperShould
    {
        [Fact(DisplayName ="Compute points distances")]
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

        [Fact(DisplayName = "Angles")]
        public void Angles()
        {
            var p = Geometry2DTinyHelper.Intersections.ProjectPointOnLine(
                new GeometryPoint2D(1, 3),
                new GeometryPoint2D(1, 1),
                new GeometryPoint2D(3, 3)
            );
            Assert.Equal(new GeometryPoint2D(2, 2), p);
        }
    }
}