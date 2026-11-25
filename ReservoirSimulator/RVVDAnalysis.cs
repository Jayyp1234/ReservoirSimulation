using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct RvvdResult(
        double[] depths,
        double[] rvValues,
        double depthMin,
        double depthMax,
        double rvMin,
        double rvMax,
        double rvAtTop,
        double rvAtBottom,
        double rvGradient,
        double rSquared)
    {
        // Sorted Tabular Arrays (in Field Units: ft and stb/Mscf)
        public readonly double[] Depths = depths;
        public readonly double[] RvValues = rvValues;

        // Core Endpoints
        public readonly double DepthMin = depthMin;         // Shallowest tabulated depth [ft]
        public readonly double DepthMax = depthMax;         // Deepest tabulated depth [ft]
        public readonly double RvMin = rvMin;               // Minimum Rv in table [stb/Mscf]
        public readonly double RvMax = rvMax;               // Maximum Rv in table [stb/Mscf]
        public readonly double RvAtTop = rvAtTop;           // Rv at shallowest depth DepthMin [stb/Mscf]
        public readonly double RvAtBottom = rvAtBottom;     // Rv at deepest depth DepthMax [stb/Mscf]

        // Linear Trend Diagnostics
        public readonly double RvGradient = rvGradient;     // Average dRv/dD slope [stb/Mscf/ft]
        public double RSquared { get; } = rSquared;

        // Fast Scalar Evaluation (for RK4 / Sub-stepped Gas Hydrostatic Equilibration)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double EvaluateRv(double depth)
        {
            int n = Depths.Length;
            if (depth <= Depths[0])
                return RvValues[0]; // Flat clamp above shallowest depth

            if (depth >= Depths[n - 1])
                return RvValues[n - 1]; // Flat clamp below deepest depth

            // Binary search for bracketing depth interval [Depths[idx], Depths[idx + 1]]
            int idx = FindInterval(depth);
            double d0 = Depths[idx];
            double d1 = Depths[idx + 1];
            double rv0 = RvValues[idx];
            double rv1 = RvValues[idx + 1];

            double dd = d1 - d0;
            if (dd <= 0.0)
                return rv0;

            double t = (depth - d0) / dd;
            return rv0 + t * (rv1 - rv0);
        }

        // Inverse Lookup: Estimate Depth for a Given Rv (assuming monotonic variation)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double EvaluateDepth(double rv)
        {
            int n = Depths.Length;
            bool increasing = RvValues[n - 1] >= RvValues[0];

            if (increasing)
            {
                if (rv <= RvValues[0]) return Depths[0];
                if (rv >= RvValues[n - 1]) return Depths[n - 1];
            }
            else
            {
                if (rv >= RvValues[0]) return Depths[0];
                if (rv <= RvValues[n - 1]) return Depths[n - 1];
            }

            for (int i = 0; i < n - 1; i++)
            {
                double r0 = RvValues[i];
                double r1 = RvValues[i + 1];
                bool inBracket = increasing ? (rv >= r0 && rv <= r1) : (rv <= r0 && rv >= r1);

                if (inBracket)
                {
                    double dr = r1 - r0;
                    if (Abs(dr) <= 1e-15) return Depths[i];
                    double t = (rv - r0) / dr;
                    return Depths[i] + t * (Depths[i + 1] - Depths[i]);
                }
            }

            return Depths[n - 1];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int FindInterval(double depth)
        {
            int lo = 0;
            int hi = Depths.Length - 2;

            while (lo <= hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (depth < Depths[mid])
                    hi = mid - 1;
                else if (depth > Depths[mid + 1])
                    lo = mid + 1;
                else
                    return mid;
            }

            return Clamp(lo, 0, Depths.Length - 2);
        }
    }

    public class RvvdRow(List<double> rvvd, Units incomingUnit)
    {
        public double Depth { get; set; } = Utility.Convert(rvvd[0], "Length", incomingUnit, Units.FIELD);
        public double Rv { get; set; } = Utility.Convert(rvvd[1], "OilGasRatio", incomingUnit, Units.FIELD);
    }

    public class RVVDAnalysis(List<List<double>> _rvvd, Units incomingUnit)
    {
        public readonly List<RvvdRow> rvvdTable = [.. _rvvd.Select(row => new RvvdRow(row, incomingUnit))];

        public RvvdResult Result { get; private set; }

        public void FitProperties()
        {
            if (rvvdTable == null || rvvdTable.Count < 2)
                throw new ArgumentException("The RVVD table must contain at least 2 rows.");

            var sorted = rvvdTable.OrderBy(r => r.Depth).ToList();

            // Validate strictly increasing depths and non-negative Rv
            for (int i = 0; i < sorted.Count; i++)
            {
                if (sorted[i].Rv < 0.0)
                    throw new ArgumentException($"Negative Rv ({sorted[i].Rv}) encountered at depth {sorted[i].Depth}.");

                if (i > 0 && sorted[i].Depth <= sorted[i - 1].Depth)
                    throw new ArgumentException($"RVVD depths must be strictly increasing (duplicate/unordered depth at {sorted[i].Depth}).");
            }

            double[] depths = [.. sorted.Select(r => r.Depth)];
            double[] rvValues = [.. sorted.Select(r => r.Rv)];

            double depthMin = depths[0];
            double depthMax = depths[^1];
            double rvMin = rvValues.Min();
            double rvMax = rvValues.Max();
            double rvAtTop = rvValues[0];
            double rvAtBottom = rvValues[^1];

            // Compute linear compositional gradient dRv/dD diagnostic
            double gradient = 0.0;
            double rSq = 1.0;

            if (sorted.Count >= 2)
            {
                var (slope, _, rSquared) = Utility.LinearRegression([.. depths], [.. rvValues]);
                gradient = slope;
                rSq = double.IsNaN(rSquared) ? 1.0 : rSquared;
            }

            Result = new(depths, rvValues, depthMin, depthMax,
                rvMin, rvMax, rvAtTop, rvAtBottom, gradient, rSq);
        }
    }
}