using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SwofResult(double swi, double sor, double krwEndpoint, double kroEndpoint,
                          double nw, double no, double pd, double lambda,
                          double pcRSquared, double krwRSquared, double kroRSquared, 
                          double SwMin, double SwMax, double PcMin, double PcMax, double swKrwEnd)
    {
        public readonly double SwKrwEnd = swKrwEnd;                // First Sw at which Krw reaches its end point

        // Core Endpoints
        public readonly double SwMin = SwMin;                      // Connate Water Saturation
        public readonly double SwMax = SwMax;                      // Maximum Water Saturation
        public readonly double PcMin = PcMin;                      // Minimum Capillary Pressure
        public readonly double PcMax = PcMax;                      // Maximum Capillary Pressure


        // Endpoints
        public readonly double Swi = swi;
        public readonly double Sor = sor;
        public readonly double KrwEndpoint = krwEndpoint;
        public readonly double KroEndpoint = kroEndpoint;

        // Corey Exponents
        public readonly double Nw = nw;
        public readonly double No = no;

        // Brooks-Corey Capillary parameters
        public readonly double Pd = pd;
        public readonly double InvLambda = 1.0 / lambda; // Pre-inverted (1.0 / lambda) for register speed

        // Diagnostics
        public double PcRSquared { get; } = pcRSquared;
        public double KrwRSquared { get; } = krwRSquared;
        public double KroRSquared { get; } = kroRSquared;

        // High-Performance Inlined AutoDiff Math Models
        // Krw rises from Swi to the first Sw where it reaches its end point (SwKrwEnd; the last row for most tables) and
        // is flat above that, like Krg in SGOF. Krow (below) spans Swi to 1 - Sor, where Sor is set by the first row with
        // Krow = 0. Before, both used 1 - Sor with Sor = 1 - SwMax, so Sor was 0 for any table ending at Sw = 1 and
        // residual oil could still flow.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrw<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            double denom = SwKrwEnd - Swi;
            ADiff<T> swNorm = (sw - Swi) / denom;

            if (swNorm.Value <= 0.0)
                return default; // 0 value and 0 derivative for Sw <= Swi

            if (swNorm.Value >= 1.0)
                return new ADiff<T>(KrwEndpoint, default); // Clamped at the end point for Sw >= SwKrwEnd

            return KrwEndpoint * Pow(swNorm, Nw);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrow<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            double denom = 1.0 - Swi - Sor;
            ADiff<T> soNorm = (1.0 - sw - Sor) / denom;

            if (soNorm.Value <= 0.0)
                return default; // 0 value and 0 derivative for Sw >= 1 - Sor

            if (soNorm.Value >= 1.0)
                return new ADiff<T>(KroEndpoint, default); // Clamped at endpoint for Sw <= Swi

            return KroEndpoint * Pow(soNorm, No);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcow<T>(ADiff<T> sw) where T : struct, IVectorND<T>
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

    public class SwofRow(List<double> swof, Units incomingUnit)
    {
        public double Sw { get; set; } = swof[0];
        public double Krw { get; set; } = swof[1];
        public double Krow { get; set; } = swof[2];
        public double Pc { get; set; } = Utility.Convert(swof[3], "Pressure", incomingUnit, Units.FIELD);
    }

    public class SWOFAnalysis(List<List<double>> _swof, Units incomingUnit)
    {
        public readonly List<SwofRow> swofTable = [.. _swof.Select(swof => new SwofRow(swof, incomingUnit))];

        public SwofResult Result { get; private set; }

        public void FitProperties()
        {
            if (swofTable == null || swofTable.Count < 3)
                throw new ArgumentException("The SWOF table must contain at least 3 rows.");

            var sorted = swofTable.OrderBy(r => r.Sw).ToList();

            // 1. Core Endpoints
            double swi = sorted.First().Sw;
            double kroEndpoint = sorted.First().Krow;

            var sorRow = sorted.FirstOrDefault(r => r.Krow <= 0.0) ?? sorted.Last();
            double swMin = sorted.First().Sw;
            double swMax = sorted.Last().Sw;
            // Pc falls as Sw rises, so the first row holds the largest Pc. PcMin/PcMax are the smallest and largest
            // values (they used to be the first and last rows, i.e. swapped, which set the transition zone to Sw = 1).
            double pcMin = sorted.Min(r => r.Pc);
            double pcMax = sorted.Max(r => r.Pc);
            double sor = 1.0 - sorRow.Sw;
            // Krw end point: the table's largest Krw, and the first Sw at which it is reached (rows after that repeat it)
            double krwEndpoint = sorted.Max(r => r.Krw);
            double swKrwEnd = sorted.First(r => r.Krw >= krwEndpoint * (1.0 - 1e-9)).Sw;
            if (swKrwEnd <= swi) swKrwEnd = swMax;    // Krw is 0 everywhere: keep a valid range

            // 2. Prepare Log Arrays
            List<double> pcX = []; List<double> pcY = [];
            List<double> krwX = []; List<double> krwY = [];
            List<double> kroX = []; List<double> kroY = [];

            foreach (var row in sorted)
            {
                if (row.Sw > swi && row.Pc > 0)
                {
                    double swStar = (row.Sw - swi) / (1.0 - swi);
                    pcX.Add(Log(swStar));
                    pcY.Add(Log(row.Pc));
                }

                // Krw over Swi..SwKrwEnd, Krow over Swi..1-Sor (each curve over its own range, see EvaluateKrw)
                if (row.Sw > swi && row.Sw < swKrwEnd && row.Krw > 0)
                {
                    krwX.Add(Log((row.Sw - swi) / (swKrwEnd - swi)));
                    krwY.Add(Log(row.Krw));
                }
                if (row.Sw > swi && row.Sw < 1.0 - sor && row.Krow > 0)
                {
                    double swNorm = (row.Sw - swi) / (1.0 - swi - sor);
                    if (row.Krw > 0)
                    {
                        krwX.Add(Log(swNorm));
                        krwY.Add(Log(row.Krw));
                    }
                    if (row.Krow > 0)
                    {
                        kroX.Add(Log(1.0 - swNorm));
                        kroY.Add(Log(row.Krow));
                    }
                }
            }

            // 3. Execute Least-Squares Regressions
            var (mPc, bPc, rSqPc) = Utility.LinearRegression(pcX, pcY);
            var (nw, _, rSqW) = Utility.LinearRegression(krwX, krwY);
            var (no, _, rSqO) = Utility.LinearRegression(kroX, kroY);

            double lambda = -1.0 / mPc;
            double pd = Exp(bPc);



            Result = new(swi, sor, krwEndpoint, kroEndpoint, nw, no, 
                        pd, lambda, rSqPc, rSqW, rSqO, swMin, swMax, pcMin, pcMax, swKrwEnd);
        }

        
    }
}