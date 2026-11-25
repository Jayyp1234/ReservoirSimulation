using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct RsvdResult(
        double[] depths,
        double[] rsValues,
        double depthMin,
        double depthMax,
        double rsMin,
        double rsMax,
        double rsAtTop,
        double rsAtBottom,
        double rsGradient,
        double rSquared)
    {
        // Sorted Tabular Arrays (in Field Units: ft and Mscf/stb)
        public readonly double[] Depths = depths;
        public readonly double[] RsValues = rsValues;

        // Core Endpoints
        public readonly double DepthMin = depthMin;         // Shallowest tabulated depth [ft]
        public readonly double DepthMax = depthMax;         // Deepest tabulated depth [ft]
        public readonly double RsMin = rsMin;               // Minimum Rs in table [Mscf/stb]
        public readonly double RsMax = rsMax;               // Maximum Rs in table [Mscf/stb]
        public readonly double RsAtTop = rsAtTop;           // Rs at shallowest depth DepthMin [Mscf/stb]
        public readonly double RsAtBottom = rsAtBottom;     // Rs at deepest depth DepthMax [Mscf/stb]

        // Linear Trend Diagnostics
        public readonly double RsGradient = rsGradient;     // Average dRs/dD slope [Mscf/stb/ft]
        public double RSquared { get; } = rSquared;

        // Fast Scalar Evaluation (for RK4 / Sub-stepped Hydrostatic Equilibration)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double EvaluateRs(double depth)
        {
            int n = Depths.Length;
            if (depth <= Depths[0])
                return RsValues[0]; // Flat clamp above shallowest depth

            if (depth >= Depths[n - 1])
                return RsValues[n - 1]; // Flat clamp below deepest depth

            // Binary search for bracketing depth interval [Depths[idx], Depths[idx + 1]]
            int idx = FindInterval(depth);
            double d0 = Depths[idx];
            double d1 = Depths[idx + 1];
            double rs0 = RsValues[idx];
            double rs1 = RsValues[idx + 1];

            double dd = d1 - d0;
            if (dd <= 0.0)
                return rs0;

            double t = (depth - d0) / dd;
            return rs0 + t * (rs1 - rs0);
        }

        
        // Inverse Lookup: Estimate Depth for a Given Rs (assuming monotonic variation)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double EvaluateDepth(double rs)
        {
            int n = Depths.Length;
            bool increasing = RsValues[n - 1] >= RsValues[0];

            if (increasing)
            {
                if (rs <= RsValues[0]) return Depths[0];
                if (rs >= RsValues[n - 1]) return Depths[n - 1];
            }
            else
            {
                if (rs >= RsValues[0]) return Depths[0];
                if (rs <= RsValues[n - 1]) return Depths[n - 1];
            }

            for (int i = 0; i < n - 1; i++)
            {
                double r0 = RsValues[i];
                double r1 = RsValues[i + 1];
                bool inBracket = increasing ? (rs >= r0 && rs <= r1) : (rs <= r0 && rs >= r1);

                if (inBracket)
                {
                    double dr = r1 - r0;
                    if (Abs(dr) <= 1e-15) return Depths[i];
                    double t = (rs - r0) / dr;
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

    public class RsvdRow(List<double> rsvd, Units incomingUnit)
    {
        public double Depth { get; set; } = Utility.Convert(rsvd[0], "Length", incomingUnit, Units.FIELD);
        public double Rs { get; set; } = Utility.Convert(rsvd[1], "Solution Gas-Oil Ratio (Rs)", incomingUnit, Units.FIELD);
    }

    public class RSVDAnalysis(List<List<double>> _rsvd, Units incomingUnit)
    {
        public readonly List<RsvdRow> rsvdTable = [.. _rsvd.Select(row => new RsvdRow(row, incomingUnit))];

        public RsvdResult Result { get; private set; }

        public void FitProperties()
        {
            if (rsvdTable == null || rsvdTable.Count < 2)
                throw new ArgumentException("The RSVD table must contain at least 2 rows.");

            var sorted = rsvdTable.OrderBy(r => r.Depth).ToList();

            // Validate strictly increasing depths and non-negative Rs
            for (int i = 0; i < sorted.Count; i++)
            {
                if (sorted[i].Rs < 0.0)
                    throw new ArgumentException($"Negative Rs ({sorted[i].Rs}) encountered at depth {sorted[i].Depth}.");

                if (i > 0 && sorted[i].Depth <= sorted[i - 1].Depth)
                    throw new ArgumentException($"RSVD depths must be strictly increasing (duplicate/unordered depth at {sorted[i].Depth}).");
            }

            double[] depths = [.. sorted.Select(r => r.Depth)];
            double[] rsValues = [.. sorted.Select(r => r.Rs)];

            double depthMin = depths[0];
            double depthMax = depths[^1];
            double rsMin = rsValues.Min();
            double rsMax = rsValues.Max();
            double rsAtTop = rsValues[0];
            double rsAtBottom = rsValues[^1];

            // Compute linear compositional gradient dRs/dD diagnostic
            double gradient = 0.0;
            double rSq = 1.0;

            if (sorted.Count >= 2)
            {
                var (slope, _, rSquared) = Utility.LinearRegression([.. depths], [.. rsValues]);
                gradient = slope;
                rSq = double.IsNaN(rSquared) ? 1.0 : rSquared;
            }

            Result = new(depths, rsValues, depthMin, depthMax,
                rsMin, rsMax, rsAtTop, rsAtBottom, gradient, rSq);
        }
    }
}