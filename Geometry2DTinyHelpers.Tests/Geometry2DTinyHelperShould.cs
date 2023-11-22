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
    }
}