using System;
using System.Collections.Generic;
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
            /// <param name="x1">Point 1 horizontal coord.</param>
            /// <param name="y1">Point 1 vertical coord.</param>
            /// <param name="x2">Point 2 horizontal coord.</param>
            /// <param name="y2">Point 2 vertical coord.</param>
            /// <returns>Distance in between the two points.</returns>
            public static double PointDistance(double x1, double y1, double x2, double y2)
                => Math.Sqrt(Math.Pow((x2 - x1), 2) + Math.Pow((y2 - y1), 2));

            /// <summary>
            /// Compute two points distance.
            /// </summary>
            /// <param name="a">First point</param>
            /// <param name="b">Second point</param>
            /// <returns>Distance in between the two points.</returns>
            public static double PointDistance(GeometryPoint2D a, GeometryPoint2D b)
                => PointDistance(a.X, a.Y, b.X, b.Y);

            /// <summary>
            /// 
            /// </summary>
            /// <param name="x"></param>
            /// <param name="y"></param>
            /// <param name="x1"></param>
            /// <param name="y1"></param>
            /// <param name="x2"></param>
            /// <param name="y2"></param>
            /// <returns></returns>
            public static double PointSegmentDistance(double x, double y, double x1, double y1, double x2, double y2)
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

            public static double PointSegmentDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b)
                => PointSegmentDistance(p.X, p.Y, a.X, a.Y, b.X, b.Y);

            public static double PointLineDistance(GeometryPoint2D p, GeometryPoint2D a, GeometryPoint2D b)
                => PointDistance(p, Intersections.ProjectPointOnSegment(a, b, p));
        }

        // -------- INTERSECTIONS -------- //
        public static class Intersections
        {
            public static GeometryPoint2D LineIntersection(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2)
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

            public static bool SegmentIntersection(GeometryPoint2D p1, GeometryPoint2D p2, GeometryPoint2D p3, GeometryPoint2D p4, out GeometryPoint2D intersection)
            {
                double Ax, Bx, Cx, Ay, By, Cy, d, e, f, num;
                double x1lo, x1hi, y1lo, y1hi;
                intersection = new GeometryPoint2D(double.NaN, double.NaN);
                Ax = p2.X - p1.X;
                Bx = p3.X - p4.X;
                if (Ax < 0)
                {
                    x1lo = p2.X;
                    x1hi = p1.X;
                }
                else
                {
                    x1hi = p2.X;
                    x1lo = p1.X;
                }
                if (Bx > 0)
                {
                    if (x1hi < p4.X || p3.X < x1lo)
                        return false;
                }
                else
                {
                    if (x1hi < p3.X || p4.X < x1lo)
                        return false;
                }
                Ay = p2.Y - p1.Y;
                By = p3.Y - p4.Y;

                if (Ay < 0)
                {
                    y1lo = p2.Y; y1hi = p1.Y;
                }
                else
                {
                    y1hi = p2.Y; y1lo = p1.Y;
                }
                if (By > 0)
                {
                    if (y1hi < p4.Y || p3.Y < y1lo)
                        return false;
                }
                else
                {
                    if (y1hi < p3.Y || p4.Y < y1lo)
                        return false;
                }

                Cx = p1.X - p3.X;
                Cy = p1.Y - p3.Y;
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
                intersection.X = p1.X + num / f;
                num = d * Ay;

                intersection.Y = p1.Y + num / f;
                return true;
            }

            public static GeometryPoint2D[] GetPerpendicularSegment(GeometryPoint2D p1, GeometryPoint2D p2, float length)
            {
                GeometryPoint2D[] result = new GeometryPoint2D[2];
                double dx = p2.X - p1.X;
                double dy = p2.Y - p1.Y;
                double d = (float)Math.Sqrt(dx * dx + dy * dy);
                double ux = dx / d;
                double uy = dy / d;
                double vx = -uy;
                double vy = ux;
                result[0] = new GeometryPoint2D(p1.X + ux * length, p1.Y + uy * length);
                result[1] = new GeometryPoint2D(p1.X + vx * length, p1.Y + vy * length);
                return result;
            }

            public static GeometryPoint2D ProjectPointOnSegment(GeometryPoint2D la, GeometryPoint2D lb, GeometryPoint2D p)
            {
                double m = (double)(lb.Y - la.Y) / (lb.X - la.X);
                double b = (double)la.Y - (m * la.X);

                double x = (m * p.Y + p.X - m * b) / (m * m + 1);
                double y = (m * m * p.Y + m * p.X + b) / (m * m + 1);

                return new GeometryPoint2D(x, y);
            }
        }

        // -------- ANGLES -------- //
        public static class Angles
        {
            public static double SegmentAngle(double x1, double y1, double x2, double y2)
            => Math.Atan2(y1 - y2, x2 - x1) * Rad2Deg;

            public static double SegmentAngle(GeometryPoint2D a, GeometryPoint2D b)
                => SegmentAngle(a.X, a.Y, b.X, b.Y);

            public static double SegmentAngleDelta(GeometryPoint2D a1, GeometryPoint2D b1, GeometryPoint2D a2, GeometryPoint2D b2)
                => SegmentAngle(a2, b2) - SegmentAngle(a1, b1);
        }

        // -------- INTERPOLATIONS -------- //
        public static class Interpolations
        {
            public static GeometryPoint2D Interpolate(GeometryPoint2D a, GeometryPoint2D b, double factor)
            => new GeometryPoint2D(Interpolate(a.X, b.X, factor), Interpolate(a.Y, b.Y, factor));

            public static GeometryPoint2D RotatePoint(GeometryPoint2D pointToRotate, GeometryPoint2D centerPoint, double angleInDegrees)
            {
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

            public static double Interpolate(double a, double b, double factor)
                => a + ((b - a) * factor);

            public static double ExponentialDecreasingInterpolation(double a, double b, double factor)
                => Interpolate(a, b, 1 - ((1 - factor) * (1 - factor)));

            public static double Normalize(double a, double b, double position)
                => (position - a) / (b - a);

            public static double Extrapolation(GeometryPoint2D _a, GeometryPoint2D _b, double x)
            {
                var m = (_b.Y - _a.Y) / (_b.X - _a.X);
                return _a.Y + m * (x - _a.X);
            }

            public static double Gamma(double x, double curvatur = 1)
            {
                x = Math.Min(1, Math.Max(0, x));
                return ((1 - x) * (x * x * curvatur)) + (x * (1 - ((1 - x) * (1 - x) * curvatur)));
            }

            public static GeometryPoint2D Sample(double slope, double offset, double x)
                => new GeometryPoint2D(x, (x * slope) + offset);
        }

        // -------- LINEARE REGRESSION -------- //
        public static class LinearRegressions
        {
            public static double[] LinearRegression(double[] x, double[] y)
            {
                if (x.Length != y.Length)
                    throw new ArgumentException($"{nameof(x)} and {nameof(y)} must have same length");
                double[] result = new double[2];
                double sumX = 0;
                double sumY = 0;
                double sumX2 = 0;
                double sumXY = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    sumX += x[i];
                    sumY += y[i];
                    sumX2 += x[i] * x[i];
                    sumXY += x[i] * y[i];
                }
                double n = x.Length;
                double a = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
                double b = (sumY - a * sumX) / n;
                result[0] = a;
                result[1] = b;
                return result;
            }

            public static double LinearCorrelationCoefficient(double[] x, double[] y)
            {
                if (x.Length != y.Length)
                    throw new ArgumentException($"{nameof(x)} and {nameof(y)} must have same length");
                double[] regression = LinearRegression(x, y);
                double a = regression[0];
                double b = regression[1];
                double sumX = 0;
                double sumY = 0;
                double sumX2 = 0;
                double sumY2 = 0;
                double sumXY = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    sumX += x[i];
                    sumY += y[i];
                    sumX2 += x[i] * x[i];
                    sumY2 += y[i] * y[i];
                    sumXY += (x[i] * y[i]);
                }
                double n = x.Length;
                double numerator = n * sumXY - sumX * sumY;
                double denominator = Math.Sqrt(n * sumX2 - sumX * sumX) * Math.Sqrt(n * sumY2 - sumY * sumY);
                double r = numerator / denominator;
                return r;
            }

            public static double[] LinearRegression(GeometryPoint2D[] pts)
            {
                double[] result = new double[2];
                double sumX = 0;
                double sumY = 0;
                double sumX2 = 0;
                double sumXY = 0;
                for (int i = 0; i < pts.Length; i++)
                {
                    sumX += pts[i].X;
                    sumY += pts[i].Y;
                    sumX2 += pts[i].X * pts[i].X;
                    sumXY += pts[i].X * pts[i].Y;
                }
                double n = pts.Length;
                double a = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
                double b = (sumY - a * sumX) / n;
                result[0] = a;
                result[1] = b;
                return result;
            }

            public static double LinearCorrelationCoefficient(GeometryPoint2D[] pts)
            {
                double[] regression = LinearRegression(pts);
                double a = regression[0];
                double b = regression[1];
                double sumX = 0;
                double sumY = 0;
                double sumX2 = 0;
                double sumY2 = 0;
                double sumXY = 0;
                for (int i = 0; i < pts.Length; i++)
                {
                    sumX += pts[i].X;
                    sumY += pts[i].Y;
                    sumX2 += pts[i].X * pts[i].X;
                    sumY2 += pts[i].Y * pts[i].Y;
                    sumXY += pts[i].X * pts[i].Y;
                }
                double n = pts.Length;
                double numerator = n * sumXY - sumX * sumY;
                double denominator = Math.Sqrt(n * sumX2 - sumX * sumX) * Math.Sqrt(n * sumY2 - sumY * sumY);
                double r = numerator / denominator;
                return r;
            }

            public static void LinearRegression(
                GeometryPoint2D[] pts,
                out double rSquared,
                out double yIntercept,
                out double slope)
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
        public static class BSpline
        {
            public static GeometryPoint2D[] Interpolate1D(GeometryPoint2D[] pts, int count)
            {
                var r = Interpolate1D(pts.Select(pt => pt.X).ToArray(), pts.Select(pt => pt.Y).ToArray(), count);
                var array = new GeometryPoint2D[r.xs.Length];
                for (int i = 0; i < r.xs.Length; i++)
                    array[i] = new GeometryPoint2D(r.xs[i], r.ys[i]);
                return array;
            }

            public static (double[] xs, double[] ys) Interpolate1D(double[] xs, double[] ys, int count)
            {
                if (xs is null || ys is null || xs.Length != ys.Length)
                    throw new ArgumentException($"{nameof(xs)} and {nameof(ys)} must have same length");

                int inputPointCount = xs.Length;
                double[] inputDistances = new double[inputPointCount];
                for (int i = 1; i < inputPointCount; i++)
                    inputDistances[i] = inputDistances[i - 1] + xs[i] - xs[i - 1];

                double meanDistance = inputDistances.Last() / (count - 1);
                double[] evenDistances = Enumerable.Range(0, count).Select(x => x * meanDistance).ToArray();
                double[] xsOut = Interpolate(inputDistances, xs, evenDistances);
                double[] ysOut = Interpolate(inputDistances, ys, evenDistances);
                return (xsOut, ysOut);
            }

            private static double[] Interpolate(double[] xOrig, double[] yOrig, double[] xInterp)
            {
                (double[] a, double[] b) = FitMatrix(xOrig, yOrig);

                double[] yInterp = new double[xInterp.Length];
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
                double[] a = new double[n - 1];
                double[] b = new double[n - 1];
                double[] r = new double[n];
                double[] A = new double[n];
                double[] B = new double[n];
                double[] C = new double[n];

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

                double[] cPrime = new double[n];
                cPrime[0] = C[0] / B[0];
                for (int i = 1; i < n; i++)
                    cPrime[i] = C[i] / (B[i] - cPrime[i - 1] * A[i]);

                double[] dPrime = new double[n];
                dPrime[0] = r[0] / B[0];
                for (int i = 1; i < n; i++)
                    dPrime[i] = (r[i] - dPrime[i - 1] * A[i]) / (B[i] - cPrime[i - 1] * A[i]);

                double[] k = new double[n];
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

        // Surfaces
        // https://codepal.ai/code-generator/query/2TJ1xgGS/calculate-triangle-surface-area
        // Point in triangle : https://stackoverflow.com/questions/25385361/point-within-a-triangle-barycentric-co-ordinates
    }
}
