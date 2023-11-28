using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

namespace Geometry2DTinyHelpers
{
    public static class Geometry2DTinyHelper
    {
        public const double Rad2Deg = 180.0 / Math.PI;
        public const double Deg2Rad = Math.PI / 180.0;

        // -------- DISTANCES -------- //
        /// <summary>
        /// Method group to compute point, point to segment, point to line distances.
        /// </summary>
        public static class Distances
        {
            /// <summary>
            /// Compute two points distance.
            /// </summary>
            /// <param name="x1">Segment point 1 horizontal coord.</param>
            /// <param name="y1">Segment point 1 vertical coord.</param>
            /// <param name="x2">Segment point 2 horizontal coord.</param>
            /// <param name="y2">Segment point 2 vertical coord.</param>
            /// <returns>Distance in between the two points.</returns>
            public static double ComputePointDistance(double x1, double y1, double x2, double y2)
                => Math.Sqrt(Math.Pow((x2 - x1), 2) + Math.Pow((y2 - y1), 2));

            /// <summary>
            /// Compute two points distance.
            /// </summary>
            /// <param name="a">First point</param>
            /// <param name="b">Second point</param>
            /// <returns>Distance in between the two points.</returns>
            public static double ComputePointDistance(GeometryPoint2D a, GeometryPoint2D b)
                => ComputePointDistance(a.X, a.Y, b.X, b.Y);

            /// <summary>
            /// Compute point to segment distance. A segment is a line section between two points.
            /// </summary>
            /// <param name="x">Point horizontal coord.</param>
            /// <param name="y">Point vertical coord.</param>
            /// <param name="x1">Segment point 1 horizontal coord.</param>
            /// <param name="y1">Segment point 1 vertical coord.</param>
            /// <param name="x2">Segment point 2 horizontal coord.</param>
            /// <param name="y2">Segment point 2 vertical coord.</param>
            /// <returns></returns>
            public static double ComputePointSegmentDistance(double x, double y, double x1, double y1, double x2, double y2)
            {
                var a = x - x1;
                var b = y - y1;
                var c = x2 - x1;
                var d = y2 - y1;

                var dot = a * c + b * d;
                var len_sq = c * c + d * d;
                var p = -1d;
                if (len_sq != 0)
                    p = dot / len_sq;

                double xx, yy;

                if (p < 0)
                {
                    xx = x1;
                    yy = y1;
                }
                else if (p > 1)
                {
                    xx = x2;
                    yy = y2;
                }
                else
                {
                    xx = x1 + p * c;
                    yy = y1 + p * d;
                }

                var dx = x - xx;
                var dy = y - yy;
                return Math.Sqrt(dx * dx + dy * dy);
            }

            /// <summary>
            /// Compute point to segment distance. A segment is a line section between two points.
            /// </summary>
            /// <param name="p">The point.</param>
            /// <param name="a">The start point of the segment.</param>
            /// <param name="b">The end point of the segment.</param>
            /// <returns></returns>
            public static double ComputePointSegmentDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b)
                => ComputePointSegmentDistance(p.X, p.Y, a.X, a.Y, b.X, b.Y);

            /// <summary>
            /// Compute point to line distance. A line is not limited.
            /// </summary>
            /// <param name="p">The point.</param>
            /// <param name="a">The first point through which the line passes.</param>
            /// <param name="b">The Second point through which the line passes</param>
            /// <returns></returns>
            public static double ComputePointLineDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b)
                => ComputePointDistance(p, Intersections.ProjectPointOnLine(p, a, b));
        }

        // -------- INTERSECTIONS -------- //
        /// <summary>
        /// Method group to compute line intersection, segment intersection, perpendiculare segment, project à point on line.
        /// </summary>
        public static class Intersections
        {
            /// <summary>
            /// Compute two line intersection. If they are parallel, a null point is returned.
            /// </summary>
            /// <param name="a1">First line point 1</param>
            /// <param name="b1">First line point 2</param>
            /// <param name="a2">Second line point 1</param>
            /// <param name="b2">Second line point 2</param>
            /// <returns>The intersection point. Null if lines are paralell.</returns>
            public static GeometryPoint2D ComputeLineIntersection(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2)
            {
                double A1 = b1.Y - a1.Y;
                double B1 = a1.X - b1.X;
                double C1 = A1 * a1.X + B1 * a1.Y;

                double A2 = b2.Y - a2.Y;
                double B2 = a2.X - b2.X;
                double C2 = A2 * a2.X + B2 * a2.Y;

                double det = A1 * B2 - A2 * B1;
                if (det == 0)
                {
                    return null;
                }
                else
                {
                    double x = (B2 * C1 - B1 * C2) / det;
                    double y = (A1 * C2 - A2 * C1) / det;
                    return new GeometryPoint2D(x, y);
                }
            }

            /// <summary>
            /// Compute two segment intersection. If they are parallel, a null point is asigned to 'intersection' parameter.
            /// </summary>
            /// <param name="a1">First segment point 1</param>
            /// <param name="b1">First segment point 2</param>
            /// <param name="a2">Second segment point 1</param>
            /// <param name="b2">Second segment point 2</param>
            /// <param name="intersection">The intersection point returned. Null if lines are paralell.</param>
            /// <returns>True if intersection is on the segments, false if the intersection is outside of segments.</returns>
            public static bool ComputeSegmentIntersection(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2, out GeometryPoint2D intersection)
            {
                double Ax, Bx, Cx, Ay, By, Cy, d, e, f, num;
                double x1lo, x1hi, y1lo, y1hi;
                intersection = new GeometryPoint2D(double.NaN, double.NaN);
                Ax = b1.X - a1.X;
                Bx = a2.X - b2.X;
                if (Ax < 0)
                {
                    x1lo = b1.X;
                    x1hi = a1.X;
                }
                else
                {
                    x1hi = b1.X;
                    x1lo = a1.X;
                }
                if (Bx > 0)
                {
                    if (x1hi < b2.X || a2.X < x1lo)
                        return false;
                }
                else
                {
                    if (x1hi < a2.X || b2.X < x1lo)
                        return false;
                }
                Ay = b1.Y - a1.Y;
                By = a2.Y - b2.Y;

                if (Ay < 0)
                {
                    y1lo = b1.Y; y1hi = a1.Y;
                }
                else
                {
                    y1hi = b1.Y; y1lo = a1.Y;
                }
                if (By > 0)
                {
                    if (y1hi < b2.Y || a2.Y < y1lo)
                        return false;
                }
                else
                {
                    if (y1hi < a2.Y || b2.Y < y1lo)
                        return false;
                }

                Cx = a1.X - a2.X;
                Cy = a1.Y - a2.Y;
                d = By * Cx - Bx * Cy;
                f = Ay * Bx - Ax * By;

                if (f > 0)
                {
                    if (d < 0 || d > f)
                        return false;
                }
                else
                {
                    if (d > 0 || d < f)
                        return false;
                }
                e = Ax * Cy - Ay * Cx;

                if (f > 0)
                {
                    if (e < 0 || e > f)
                        return false;
                }
                else
                {
                    if (e > 0 || e < f)
                        return false;
                }

                // -------- check if they are parallel
                if (f == 0) return false;

                // -------- compute intersection coordinates
                num = d * Ax;
                intersection.X = a1.X + num / f;
                num = d * Ay;

                intersection.Y = a1.Y + num / f;
                return true;
            }

            /// <summary>
            /// Compute a segment points of length 'lenght', perpendicular to the given segment and starting from 'a' point.
            /// </summary>
            /// <param name="a">Segment first point.</param>
            /// <param name="b">Segment second point.</param>
            /// <param name="length">THe length of the returned perpendicular segment.</param>
            /// <returns>Two points array that describe the computed segment.</returns>
            public static GeometryPoint2D[] ComputePerpendicularLine(GeometryPoint2D a, GeometryPoint2D b, float length)
            {
                GeometryPoint2D[] result = new GeometryPoint2D[2];
                double dx = b.X - a.X;
                double dy = b.Y - a.Y;
                double d = (float)Math.Sqrt(dx * dx + dy * dy);
                double ux = dx / d;
                double uy = dy / d;
                double vx = -uy;
                double vy = ux;
                result[0] = new GeometryPoint2D(a.X, a.Y);
                result[1] = new GeometryPoint2D(a.X + vx * length, a.Y + vy * length);
                return result;
            }

            /// <summary>
            /// Compute a point projected on a line.
            /// </summary>
            /// <param name="p">The point to project on the line.</param>
            /// <param name="a">Line first point.</param>
            /// <param name="b">Line second point.</param>
            /// <returns></returns>
            public static GeometryPoint2D ProjectPointOnLine(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b)
            {
                double m = (double)(b.Y - a.Y) / (b.X - a.X);
                double _b = (double)a.Y - (m * a.X);

                double x = (m * p.Y + p.X - m * _b) / (m * m + 1);
                double y = (m * m * p.Y + m * p.X + _b) / (m * m + 1);

                return new GeometryPoint2D(x, y);
            }

            /// <summary>
            /// Check if a point is on a segment
            /// </summary>
            /// <param name="p">The point to check.</param>
            /// <param name="a">Segment first point.</param>
            /// <param name="b">Segment second point.</param>
            /// <param name="epsilon">The minimal value considered as an error.</param>
            /// <returns>True if the point is considered on the segment, taking car of epsilon.</returns>
            public static bool IsPointOnSegment(Point p, GeometryPoint2D a, GeometryPoint2D b, double epsilon = 0.0001)
            {
                var minX = Math.Min(a.X, b.X);
                var maxX = Math.Max(a.X, b.X);
                var minY = Math.Min(a.Y, b.Y);
                var maxY = Math.Max(a.Y, b.Y);
                if (!(minX <= p.X) || !(p.X <= maxX) || !(minY <= p.Y) || !(p.Y <= maxY))
                    return false;
                if (Math.Abs(a.X - b.X) < epsilon)
                    return Math.Abs(a.X - p.X) < epsilon || Math.Abs(b.X - p.X) < epsilon;
                if (Math.Abs(a.Y - b.Y) < epsilon)
                    return Math.Abs(a.Y - p.Y) < epsilon || Math.Abs(b.Y - p.Y) < epsilon;
                return Math.Abs((p.X - a.X) / (b.X - a.X) - (p.Y - a.Y) / (b.Y - a.Y)) < epsilon);
            }
        }

        // -------- ANGLES -------- //
        /// <summary>
        /// Method group to compute segment angles.
        /// </summary>
        public static class Angles
        {
            /// <summary>
            /// Compute a single segment angle, or inclination. An horizontale line is 0 degree (East origin).
            /// Upper angles are negative, lower angles are positive.
            /// </summary>
            /// <param name="x1">Segment point 1 horizontal coord.</param>
            /// <param name="y1">Segment point 1 vertical coord.</param>
            /// <param name="x2">Segment point 2 horizontal coord.</param>
            /// <param name="y2">Segment point 2 vertical coord.</param>
            /// <returns>The computed angle, in degree.</returns>
            public static double ComputeSegmentAngle(double x1, double y1, double x2, double y2)
                => Math.Atan2(y1 - y2, x2 - x1) * Rad2Deg;

            /// <summary>
            /// Compute a single segment angle, or inclination.
            /// </summary>
            /// <param name="a">Segment first point.</param>
            /// <param name="b">Segment second point.</param>
            /// <returns>The computed angle, in degree.</returns>
            public static double ComputeSegmentAngle(GeometryPoint2D a, GeometryPoint2D b)
                => ComputeSegmentAngle(a.X, a.Y, b.X, b.Y);

            /// <summary>
            /// Convert à line slope to angle in degree.
            /// </summary>
            /// <param name="slop">THe line slope.</param>
            /// <returns>The angle in degree.</returns>
            public static double ConvertSlopeToAngle(double slop)
                => ComputeSegmentAngle(new GeometryPoint2D(0, 0), new GeometryPoint2D(1, 1 * slop));

            /// <summary>
            /// Compute a single segment angle, or inclination. An horizontale line is 0 degree (East oriented).
            /// Angles are positive, from 0° to near 360°, and clockwise.
            /// </summary>
            /// <param name="a">Segment first point.</param>
            /// <param name="b">Segment second point.</param>
            /// <returns></returns>
            public static double ComputeSegmentPositiveAngle(GeometryPoint2D a, GeometryPoint2D b)
            {
                var r = ComputeSegmentAngle(a.X, a.Y, b.X, b.Y);
                return r < 0 ? 180d + (180d + r) : r;
            }

            /// <summary>
            /// Convert à line slope to angle in degree, angles are positive, from 0° to near 360°, and clockwise.
            /// </summary>
            /// <param name="slop">THe line slope.</param>
            /// <returns>The angle in degree.</returns>
            public static double ComputeSlopeToPositiveAngle(double slop)
                => ComputeSegmentPositiveAngle(new GeometryPoint2D(0, 0), new GeometryPoint2D(1, 1 * slop));

            /// <summary>
            /// Compute two segment angle. Segments do not have to share any points. The result angle is all time in between 0° and 180°.
            /// </summary>
            /// <param name="a1">First segment point 1</param>
            /// <param name="b1">First segment point 2</param>
            /// <param name="a2">Second segment point 1</param>
            /// <param name="b2">Second segment point 2</param>
            /// <returns>The angle formed by the two segments, from 0° to 180°.</returns>
            public static double ComputeSegmentAngleDelta(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2)
            {
                var d = Math.Abs(ComputeSegmentPositiveAngle(a2, b2) - ComputeSegmentPositiveAngle(a1, b1));
                return d > 180 ? d - 180 : d;
            }

            private static double BornIt(double angle)
                => angle >= 360 ? angle - 360 : angle;

            /// <summary>
            /// Convert an East angle to North origin one.
            /// </summary>
            /// <param name="angle">The angle to convert.</param>
            /// <returns>The converted angle.</returns>
            public static double ToNorth(double angle)
                => BornIt(angle + 90 > 360 ? angle + 90 - 360 : angle + 90);

            /// <summary>
            /// Convert an East angle to West origin one.
            /// </summary>
            /// <param name="angle">The angle to convert.</param>
            /// <returns>The converted angle.</returns>
            public static double ToWest(double angle)
                => BornIt(angle + 180 > 360 ? angle + 180 - 360 : angle + 180);

            /// <summary>
            /// Convert an East angle to South origin one.
            /// </summary>
            /// <param name="angle">The angle to convert.</param>
            /// <returns>The converted angle.</returns>
            public static double ToSouth(double angle)
                => BornIt(angle + 270 > 360 ? angle + 270 - 360 : angle + 270);
        }

        // -------- INTERPOLATIONS -------- //
        /// <summary>
        /// Method group to interpolate segments, rotate points, makes intervale computations.
        /// </summary>
        public static class Interpolations
        {
            /// <summary>
            /// Interpolate / extrapolate a segment.
            /// </summary>
            /// <param name="a">First point of the segment.</param>
            /// <param name="b">Second point of the segment.</param>
            /// <param name="factor">Where, between (or over / under) 0 and 1 the point must be computed.</param>
            /// <returns>A point on the segment (0 to 1), or on the line.</returns>
            public static GeometryPoint2D Interpolate(GeometryPoint2D a, GeometryPoint2D b, double factor)
                => new GeometryPoint2D(Interpolate(a.X, b.X, factor), Interpolate(a.Y, b.Y, factor));

            /// <summary>
            /// Compate à rotated point around an another point.
            /// </summary>
            /// <param name="pointToRotate">The point to rotate.</param>
            /// <param name="centerPoint">The point to rotate around.</param>
            /// <param name="angleInDegrees">The angle in degree.</param>
            /// <returns>The rotated point.</returns>
            public static GeometryPoint2D RotatePoint(GeometryPoint2D pointToRotate, GeometryPoint2D centerPoint, double angleInDegrees)
            {
                angleInDegrees = -angleInDegrees;
                double angleInRadians = angleInDegrees * (Math.PI / 180);
                double cosTheta = Math.Cos(angleInRadians);
                double sinTheta = Math.Sin(angleInRadians);
                return new GeometryPoint2D()
                {
                    X = (cosTheta * (pointToRotate.X - centerPoint.X) -
                         sinTheta * (pointToRotate.Y - centerPoint.Y) + centerPoint.X),
                    Y = (sinTheta * (pointToRotate.X - centerPoint.X) +
                         cosTheta * (pointToRotate.Y - centerPoint.Y) + centerPoint.Y)
                };
            }

            /// <summary>
            /// Compute à value in between two values, if in between 0 ans 1 - or over / under if not.
            /// </summary>
            /// <param name="a">The first value.</param>
            /// <param name="b">The second valoue.</param>
            /// <param name="factor">A position factor (between 0 and 1, or not).</param>
            /// <returns></returns>
            private static double Interpolate(double a, double b, double factor)
                => a + ((b - a) * factor);

            /// <summary>
            /// Compute a decreased exponential interpolation value between two values.
            /// </summary>
            /// <param name="a">The first value.</param>
            /// <param name="b">THe second value.</param>
            /// <param name="factor">THe factor.</param>
            /// <returns>The result value.</returns>
            public static double ExponentialDecreasingInterpolation(double a, double b, double factor)
                => Interpolate(a, b, 1 - ((1 - factor) * (1 - factor)));

            /// <summary>
            /// Compute the position factor of a value in between (or not) two values.
            /// </summary>
            /// <param name="a">The first value.</param>
            /// <param name="b">The second value.</param>
            /// <param name="position">The position value, to convert in factore.</param>
            /// <returns>The computed factor.</returns>
            public static double Normalize(double a, double b, double position)
                => (position - a) / (b - a);

            /// <summary>
            /// Compute for a line the Y vertical coords for the given X horizontal position.
            /// </summary>
            /// <param name="a">Line first point.</param>
            /// <param name="b">Line second point.</param>
            /// <param name="x">The horizontal position to compute Y value.</param>
            /// <returns>THe Y vertical point coord for the given X horizontale coord.</returns>
            public static double Sample(GeometryPoint2D a, GeometryPoint2D b, double x)
            {
                var m = (b.Y - a.Y) / (b.X - a.X);
                return a.Y + m * (x - a.X);
            }

            /// <summary>
            /// Sample à curved function between 0 and 1, returning a value between 0 ans 1.
            /// </summary>
            /// <param name="x">The position, between 0 and 1.</param>
            /// <param name="curvatur">A ease factore (1 = fully eased, 0 is straight)</param>
            /// <returns>The computed result, between 0 ans 1.</returns>
            public static double GetGamma(double x, double curvatur = 1)
            {
                x = Math.Min(1, Math.Max(0, x));
                return ((1 - x) * (x * x * curvatur)) + (x * (1 - ((1 - x) * (1 - x) * curvatur)));
            }

            /// <summary>
            /// Compute the point on a line for the given position.
            /// </summary>
            /// <param name="slope">The slope of the line.</param>
            /// <param name="yIntercept">The Y vertical offset of the line.</param>
            /// <param name="x">The horizontale X parameter.</param>
            /// <returns>The point on the line for the given X horizontale position.</returns>
            public static GeometryPoint2D Sample(double slope, double yIntercept, double x)
                => new GeometryPoint2D(x, (x * slope) + yIntercept);
        }

        // -------- LINEARE REGRESSION -------- //
        /// <summary>
        /// Method to compute a linear regression and linear correlation.
        /// </summary>
        public static class LinearRegressions
        {
            /// <summary>
            /// Compute the correlation coefficient (or 'r') of a point serie. If points are perfectly aligned, result correlation is 1. If they are more randomly dispersed, correlation decrease to near 0.
            /// </summary>
            /// <param name="pts">Points to compute linear correlation.</param>
            /// <returns>The result coefficient (1 if points are aligned)</returns>
            public static double ComputeCorrelationCoefficient(GeometryPoint2D[] pts)
            {
                ComputeLinearRegression(pts, out var rSquared, out var yIntercept, out var slope);
                return rSquared;
            }

            /// <summary>
            /// Compute the regression line parameters (slope and Y offset) and correlation coeficient of a point serie.
            /// </summary>
            /// <param name="pts">Points to compute linear regression.</param>
            /// <param name="rSquared">The correlation coefficient (1 if points are aligned).</param>
            /// <param name="yIntercept">Y vertical offset.</param>
            /// <param name="slope">Slope of the line.</param>
            public static void ComputeLinearRegression(GeometryPoint2D[] pts, out double rSquared, out double yIntercept, out double slope)
            {
                double sumOfX = 0;
                double sumOfY = 0;
                double sumOfXSq = 0;
                double sumOfYSq = 0;
                double sumCodeviates = 0;

                for (var i = 0; i < pts.Length; i++)
                {
                    var x = pts[i].X;
                    var y = pts[i].Y;
                    sumCodeviates += x * y;
                    sumOfX += x;
                    sumOfY += y;
                    sumOfXSq += x * x;
                    sumOfYSq += y * y;
                }

                var count = pts.Length;
                var ssX = sumOfXSq - ((sumOfX * sumOfX) / count);
                var ssY = sumOfYSq - ((sumOfY * sumOfY) / count);

                var rNumerator = (count * sumCodeviates) - (sumOfX * sumOfY);
                var rDenom = (count * sumOfXSq - (sumOfX * sumOfX)) * (count * sumOfYSq - (sumOfY * sumOfY));
                var sCo = sumCodeviates - ((sumOfX * sumOfY) / count);

                var meanX = sumOfX / count;
                var meanY = sumOfY / count;
                var dblR = rNumerator / Math.Sqrt(rDenom);

                rSquared = dblR * dblR;
                yIntercept = meanY - ((sCo / ssX) * meanX);
                slope = sCo / ssX;
            }
        }

        // -------- B-SPLINE -------- //
        /// <summary>
        /// Methods to tessellate (get the intermediate points) a B-Spline curve.
        /// </summary>
        public static class BSpline
        {
            /// <summary>
            /// Return a the point array that form the polyline that represent the spline curv.
            /// </summary>
            /// <param name="pts">Tensor points.</param>
            /// <param name="count">Number of point to generate.</param>
            /// <returns>The polyline points.</returns>
            public static GeometryPoint2D[] Interpolate1D(GeometryPoint2D[] pts, int count)
            {
                var r = Interpolate1D(pts.Select(pt => pt.X).ToArray(), pts.Select(pt => pt.Y).ToArray(), count);
                var array = new GeometryPoint2D[r.xs.Length];
                for (int i = 0; i < r.xs.Length; i++)
                    array[i] = new GeometryPoint2D(r.xs[i], r.ys[i]);
                return array;
            }

            /// <summary>
            /// Return a the point array that form the polyline that represent the spline curv.
            /// </summary>
            /// <param name="xs">Point X coords array.</param>
            /// <param name="ys">Point Y coords array.</param>
            /// <param name="count">Number of point to generate.</param>
            /// <returns>Point X and Y coords.</returns>
            /// <exception cref="ArgumentException"></exception>
            public static (double[] xs, double[] ys) Interpolate1D(double[] xs, double[] ys, int count)
            {
                if (xs is null || ys is null || xs.Length != ys.Length)
                    throw new ArgumentException($"{nameof(xs)} and {nameof(ys)} must have same length");

                int inputPointCount = xs.Length;
                var inputDistances = new double[inputPointCount];
                for (int i = 1; i < inputPointCount; i++)
                    inputDistances[i] = inputDistances[i - 1] + xs[i] - xs[i - 1];

                double meanDistance = inputDistances.Last() / (count - 1);
                var evenDistances = Enumerable.Range(0, count).Select(x => x * meanDistance).ToArray();
                var xsOut = Interpolate(inputDistances, xs, evenDistances);
                var ysOut = Interpolate(inputDistances, ys, evenDistances);
                return (xsOut, ysOut);
            }

            private static double[] Interpolate(double[] xOrig, double[] yOrig, double[] xInterp)
            {
                (double[] a, double[] b) = FitMatrix(xOrig, yOrig);

                var yInterp = new double[xInterp.Length];
                for (int i = 0; i < yInterp.Length; i++)
                {
                    int j;
                    for (j = 0; j < xOrig.Length - 2; j++)
                        if (xInterp[i] <= xOrig[j + 1])
                            break;

                    double dx = xOrig[j + 1] - xOrig[j];
                    double t = (xInterp[i] - xOrig[j]) / dx;
                    double y = (1 - t) * yOrig[j] + t * yOrig[j + 1] +
                        t * (1 - t) * (a[j] * (1 - t) + b[j] * t);
                    yInterp[i] = y;
                }

                return yInterp;
            }

            private static (double[] a, double[] b) FitMatrix(double[] x, double[] y)
            {
                int n = x.Length;
                var a = new double[n - 1];
                var b = new double[n - 1];
                var r = new double[n];
                var A = new double[n];
                var B = new double[n];
                var C = new double[n];

                double dx1, dx2, dy1, dy2;

                dx1 = x[1] - x[0];
                C[0] = 1.0f / dx1;
                B[0] = 2.0f * C[0];
                r[0] = 3 * (y[1] - y[0]) / (dx1 * dx1);

                for (int i = 1; i < n - 1; i++)
                {
                    dx1 = x[i] - x[i - 1];
                    dx2 = x[i + 1] - x[i];
                    A[i] = 1.0f / dx1;
                    C[i] = 1.0f / dx2;
                    B[i] = 2.0f * (A[i] + C[i]);
                    dy1 = y[i] - y[i - 1];
                    dy2 = y[i + 1] - y[i];
                    r[i] = 3 * (dy1 / (dx1 * dx1) + dy2 / (dx2 * dx2));
                }

                dx1 = x[n - 1] - x[n - 2];
                dy1 = y[n - 1] - y[n - 2];
                A[n - 1] = 1.0f / dx1;
                B[n - 1] = 2.0f * A[n - 1];
                r[n - 1] = 3 * (dy1 / (dx1 * dx1));

                var cPrime = new double[n];
                cPrime[0] = C[0] / B[0];
                for (int i = 1; i < n; i++)
                    cPrime[i] = C[i] / (B[i] - cPrime[i - 1] * A[i]);

                var dPrime = new double[n];
                dPrime[0] = r[0] / B[0];
                for (int i = 1; i < n; i++)
                    dPrime[i] = (r[i] - dPrime[i - 1] * A[i]) / (B[i] - cPrime[i - 1] * A[i]);

                var k = new double[n];
                k[n - 1] = dPrime[n - 1];
                for (int i = n - 2; i >= 0; i--)
                    k[i] = dPrime[i] - cPrime[i] * k[i + 1];

                for (int i = 1; i < n; i++)
                {
                    dx1 = x[i] - x[i - 1];
                    dy1 = y[i] - y[i - 1];
                    a[i - 1] = k[i - 1] * dx1 - dy1;
                    b[i - 1] = -k[i] * dx1 + dy1;
                }

                return (a, b);
            }
        }

        // -------- SURFACES -------- //
        /// <summary>
        /// Methods to compute triangle area and point inclusion.
        /// </summary>
        public static class Surfaces
        {
            /// <summary>
            /// Compute the surface of a triangle.
            /// </summary>
            /// <param name="a">First point of the triangle.</param>
            /// <param name="b">Second point of the triangle.</param>
            /// <param name="c">Third point of the triangle.</param>
            /// <returns>Area of the triangle.</returns>
            public static double ComputeTriangleArea(GeometryPoint2D a, GeometryPoint2D b, GeometryPoint2D c)
            {
                var la = Distances.ComputePointDistance(a, b);
                var lb = Distances.ComputePointDistance(b, c);
                var lc = Distances.ComputePointDistance(c, a);
                var s = (la + lb + lc) / 2;
                return Math.Sqrt(s * (s - la) * (s - lb) * (s - lc));
            }

            /// <summary>
            /// Chack if a point is in a triangle.
            /// </summary>
            /// <param name="p">Point to test for inclusion.</param>
            /// <param name="a">First point of the triangle.</param>
            /// <param name="b">Second point of the triangle.</param>
            /// <param name="c">Third point of the triangle.</param>
            /// <returns>True if the point is within the triangle.</returns>
            public static bool IsPointInTriangle(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b, GeometryPoint2D c)
            {
                var det = (b.Y - c.Y) * (a.X - c.X) + (c.X - b.X) * (a.Y - c.Y);
                var factor_alpha = (b.Y - c.Y) * (p.X - c.X) + (c.X - b.X) * (p.Y - c.Y);
                var factor_beta = (c.Y - a.Y) * (p.X - c.X) + (a.X - c.X) * (p.Y - c.Y);
                var alpha = factor_alpha / det;
                var beta = factor_beta / det;
                var gamma = 1.0 - alpha - beta;
                var into = false;
                if (((a.X == p.X) & (a.Y == p.Y)) | ((b.X == p.X) & (b.Y == p.Y)) | ((c.X == p.X) & (c.Y == p.Y)))
                    into = true;
                if ((alpha == 0) | (beta == 0) | (gamma == 0))
                    into = true;
                if (((0 < alpha) & (alpha < 1)) & ((0 < beta) & (beta < 1)) & ((0 < gamma) & (gamma < 1)))
                    into = true;
                return into;
            }
        }
    }
}
