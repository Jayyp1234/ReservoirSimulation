using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SwgfResult(double swi, double sgr, double krwEndpoint, double krgEndpoint,
                          double nw, double ng, double pd, double lambda,
                          double pcRSquared, double krwRSquared, double krgRSquared,
                          double SwMin, double SwMax, double PcMin, double PcMax, double swKrwEnd)
    {
        public readonly double SwKrwEnd = swKrwEnd;                // First Sw at which Krw reaches its end point

        // Core Endpoints
        public readonly double SwMin = SwMin;                      // Connate Water Saturation
        public readonly double SwMax = SwMax;                      // Maximum Water Saturation
        public readonly double PcMin = PcMin;                      // Minimum Capillary Pressure
        public readonly double PcMax = PcMax;                      // Maximum Capillary Pressure

        // Saturation & Relative Permeability Endpoints
        public readonly double Swi = swi;
        public readonly double Sgr = sgr;
        public readonly double KrwEndpoint = krwEndpoint;
        public readonly double KrgEndpoint = krgEndpoint;

        // Corey Exponents
        public readonly double Nw = nw;
        public readonly double Ng = ng;

        // Brooks-Corey Capillary Parameters
        public readonly double Pd = pd;
        public readonly double InvLambda = 1.0 / lambda; // Pre-inverted (1.0 / lambda) for register speed

        // Diagnostics
        public double PcRSquared { get; } = pcRSquared;
        public double KrwRSquared { get; } = krwRSquared;
        public double KrgRSquared { get; } = krgRSquared;

        // High-Performance Inlined AutoDiff Math Models
        // Krw rises from Swi to SwKrwEnd and is flat above that.
        // Krg spans Swi to 1 - Sgr, where Sgr is set by the first row with Krg = 0.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrw<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            double denom = SwKrwEnd - Swi;
            ADiff<T> swNorm = (sw - Swi) / denom;

            if (swNorm.Value <= 0.0)
                return default; // 0 value and 0 derivative for Sw <= Swi

            if (swNorm.Value >= 1.0)
                return new ADiff<T>(KrwEndpoint, default); // Clamped at endpoint for Sw >= SwKrwEnd

            return KrwEndpoint * Pow(swNorm, Nw);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrg<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            double denom = 1.0 - Swi - Sgr;
            ADiff<T> sgNorm = (1.0 - sw - Sgr) / denom;

            if (sgNorm.Value <= 0.0)
                return default; // 0 value and 0 derivative for Sw >= 1 - Sgr (Sg <= Sgr)

            if (sgNorm.Value >= 1.0)
                return new ADiff<T>(KrgEndpoint, default); // Clamped at endpoint for Sw <= Swi

            return KrgEndpoint * Pow(sgNorm, Ng);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcgw<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            const double eps = 1e-4;
            double denom = 1.0 - Swi;
            ADiff<T> swStar = (sw - Swi) / denom;

            if (swStar.Value >= eps)
                return Pd * Pow(swStar, -InvLambda);

            // C1-continuous linear extension for swStar < eps (prevents NaN/Inf at sw <= Swi)
            double pcAtEps = Pd * Pow(eps, -InvLambda);
            double dPcDswStarAtEps = -Pd * InvLambda * Pow(eps, -InvLambda - 1.0);
            return pcAtEps + dPcDswStarAtEps * (swStar - eps);
        }
    }

    public class SwgfRow(List<double> swgf, Units incomingUnit)
    {
        public double Sw { get; set; } = swgf[0];
        public double Krw { get; set; } = swgf[1];
        public double Krg { get; set; } = swgf[2];
        public double Pc { get; set; } = Utility.Convert(swgf[3], "Pressure", incomingUnit, Units.FIELD);
    }

    public class SWGFAnalysis(List<List<double>> _swgf, Units incomingUnit)
    {
        public readonly List<SwgfRow> swgfTable = [.. _swgf.Select(swgf => new SwgfRow(swgf, incomingUnit))];

        public SwgfResult Result { get; private set; }

        public void FitProperties()
        {
            if (swgfTable == null || swgfTable.Count < 3)
                throw new ArgumentException("The SWGF table must contain at least 3 rows.");

            var sorted = swgfTable.OrderBy(r => r.Sw).ToList();

            // 1. Core Endpoints
            double swi = sorted.First().Sw;
            double krgEndpoint = sorted.First().Krg;

            var sgrRow = sorted.FirstOrDefault(r => r.Krg <= 0.0) ?? sorted.Last();
            double swMin = sorted.First().Sw;
            double swMax = sorted.Last().Sw;
            double pcMin = sorted.Min(r => r.Pc);
            double pcMax = sorted.Max(r => r.Pc);
            double sgr = 1.0 - sgrRow.Sw;

            // Krw end point: largest Krw in the table, and the first Sw at which it is reached
            double krwEndpoint = sorted.Max(r => r.Krw);
            double swKrwEnd = sorted.First(r => r.Krw >= krwEndpoint * (1.0 - 1e-9)).Sw;
            if (swKrwEnd <= swi) swKrwEnd = swMax; // Krw is 0 everywhere: keep a valid range

            // 2. Prepare Log Arrays
            List<double> pcX = []; List<double> pcY = [];
            List<double> krwX = []; List<double> krwY = [];
            List<double> krgX = []; List<double> krgY = [];

            foreach (var row in sorted)
            {
                if (row.Sw > swi && row.Pc > 0.0)
                {
                    double swStar = (row.Sw - swi) / (1.0 - swi);
                    pcX.Add(Log(swStar));
                    pcY.Add(Log(row.Pc));
                }

                // Krw spans Swi .. SwKrwEnd
                if (row.Sw > swi && row.Sw < swKrwEnd && row.Krw > 0.0)
                {
                    double swNormW = (row.Sw - swi) / (swKrwEnd - swi);
                    krwX.Add(Log(swNormW));
                    krwY.Add(Log(row.Krw));
                }

                // Krg spans Swi .. 1 - Sgr
                if (row.Sw > swi && row.Sw < 1.0 - sgr && row.Krg > 0.0)
                {
                    double sgNorm = (1.0 - row.Sw - sgr) / (1.0 - swi - sgr);
                    krgX.Add(Log(sgNorm));
                    krgY.Add(Log(row.Krg));
                }
            }

            // 3. Execute Least-Squares Regressions
            var (mPc, bPc, rSqPc) = Utility.LinearRegression(pcX, pcY);
            var (nw, _, rSqW) = Utility.LinearRegression(krwX, krwY);
            var (ng, _, rSqG) = Utility.LinearRegression(krgX, krgY);

            double lambda = -1.0 / mPc;
            double pd = Exp(bPc);

            Result = new(swi, sgr, krwEndpoint, krgEndpoint, nw, ng,
                         pd, lambda, rSqPc, rSqW, rSqG, swMin, swMax, pcMin, pcMax, swKrwEnd);
        }
    }
}