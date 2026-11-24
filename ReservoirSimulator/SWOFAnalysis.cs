using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SwofResult(double swi, double sor, double krwEndpoint, double kroEndpoint,
                          double nw, double no, double pd, double lambda,
                          double pcRSquared, double krwRSquared, double kroRSquared, 
                          double SwMin, double SwMax, double PcMin, double PcMax)
    {
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
        // Krw spans the whole table, Swi to the last row's Sw (SwMax), ending at the last row's Krw. Krow (below) spans
        // Swi to 1 - Sor, where Sor is set by the first row with Krow = 0. Before, both used 1 - Sor with Sor = 1 - SwMax,
        // so Sor was 0 for any table ending at Sw = 1 and residual oil could still flow.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrw<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            double denom = SwMax - Swi;
            ADiff<T> swNorm = (sw - Swi) / denom;

            if (swNorm.Value <= 0.0)
                return default; // 0 value and 0 derivative for Sw <= Swi

            if (swNorm.Value >= 1.0)
                return new ADiff<T>(KrwEndpoint, default); // Clamped at the last row's Krw for Sw >= SwMax

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

        // Brooks-Corey Pc, never above the table's largest Pc (PcMax). Near Swi the fitted curve rises without limit
        // (thousands of psi where the table says a few), which put the starting state out of equilibrium. A table with
        // no positive Pc has Pd = 0, so Pc = 0 everywhere (before, the empty fit gave Pc rising from 0 to 1 psi).
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcow<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            if (Pd <= 0.0 || PcMax <= 0.0)
                return default;

            const double eps = 1e-4;
            double denom = 1.0 - Swi;
            ADiff<T> swStar = (sw - Swi) / denom;

            ADiff<T> pc;
            if (swStar.Value >= eps)
                pc = Pd * Pow(swStar, -InvLambda);
            else
            {
                // C1-continuous linear extension for swStar < eps (prevents NaN/Inf at sw <= Swi)
                double pcAtEps = Pd * Pow(eps, -InvLambda);
                double dPcDswStarAtEps = -Pd * InvLambda * Pow(eps, -InvLambda - 1.0);
                pc = pcAtEps + dPcDswStarAtEps * (swStar - eps);
            }
            return pc.Value > PcMax ? new ADiff<T>(PcMax, default) : pc;
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

            // Residual oil: Krow reaches 0 at Sw = 1 - Sor (first such row); the last row's Sw when Krow never reaches 0
            var sorRow = sorted.FirstOrDefault(r => r.Krow <= 0.0) ?? sorted.Last();
            double swMin = sorted.First().Sw;
            double swMax = sorted.Last().Sw;
            // Pc falls as Sw rises, so the first row holds the largest Pc. PcMin/PcMax are the smallest and largest
            // values (they used to be the first and last rows, i.e. swapped, which set the transition zone to Sw = 1).
            double pcMin = sorted.Min(r => r.Pc);
            double pcMax = sorted.Max(r => r.Pc);
            double sor = 1.0 - sorRow.Sw;
            double krwEndpoint = sorted.Last().Krw;

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

                // Krw over Swi..SwMax, Krow over Swi..1-Sor (each curve over its own range, see EvaluateKrw)
                if (row.Sw > swi && row.Sw < swMax && row.Krw > 0)
                {
                    krwX.Add(Log((row.Sw - swi) / (swMax - swi)));
                    krwY.Add(Log(row.Krw));
                }
                if (row.Sw > swi && row.Sw < 1.0 - sor && row.Krow > 0)
                {
                    kroX.Add(Log((1.0 - row.Sw - sor) / (1.0 - swi - sor)));
                    kroY.Add(Log(row.Krow));
                }
            }

            // 3. Execute Least-Squares Regressions. The slopes are the Corey exponents; the curves are
            // endpoint × (normalised saturation)^n, so they pass through the table's end points
            var (mPc, bPc, rSqPc) = Utility.LinearRegression(pcX, pcY);
            var (nw, _, rSqW) = Utility.LinearRegression(krwX, krwY);
            var (no, _, rSqO) = Utility.LinearRegression(kroX, kroY);

            double lambda = -1.0 / mPc;
            double pd = Exp(bPc);
            if (pcX.Count < 2) { pd = 0.0; lambda = 1.0; }   // fewer than 2 positive Pc rows: no capillary pressure



            Result = new(swi, sor, krwEndpoint, kroEndpoint, nw, no, 
                        pd, lambda, rSqPc, rSqW, rSqO, swMin, swMax, pcMin, pcMax);
        }

        
    }
}