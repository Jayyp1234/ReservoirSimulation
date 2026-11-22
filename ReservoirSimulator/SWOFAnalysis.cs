using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SwofResult(double swi, double sor, double krwEndpoint, double kroEndpoint,
                          double nw, double no, double pd, double lambda,
                          double pcRSquared, double krwRSquared, double kroRSquared)
    {
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
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrw<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            // Normalization: swNorm = (sw - swi) / (1.0 - swi - sor)
            ADiff<T> swNorm = (sw - Swi) / (1.0 - Swi - Sor);
            return KrwEndpoint * Pow(swNorm, Nw);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrow<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            ADiff<T> swNorm = (sw - Swi) / (1.0 - Swi - Sor);
            return KroEndpoint * Pow(1.0 - swNorm, No);
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
            double swMax = sorRow.Sw;
            double sor = 1.0 - swMax;
            double krwEndpoint = sorRow.Krw;

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

                if (row.Sw > swi && row.Sw < swMax)
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

            Result = new(swi, sor, krwEndpoint, kroEndpoint, nw, no, pd, lambda, rSqPc, rSqW, rSqO);
        }

        
    }
}