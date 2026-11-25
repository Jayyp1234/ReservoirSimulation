using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SgofResult(double sgc, double sorg, double swi, double krgEndpoint, double krogEndpoint,
                                          double ng, double nog, double pdgo, double lambdaGo,
                                          double pcRSquared, double krgRSquared, double krogRSquared,
                                          double sgMin, double sgMax, double pcMin, double pcMax, double sgKrgEnd)
    {
        public readonly double Sgc = sgc;                      // Critical Gas Saturation (largest Sg with Krg = 0)
        public readonly double Sorg = sorg;                    // Residual Oil to Gas
        public readonly double Swi = swi;                      // Connate Water background parameter
        public readonly double KrgEndpoint = krgEndpoint;      // Largest Krg of the table
        public readonly double SgKrgEnd = sgKrgEnd;            // First Sg at which Krg reaches that value
        public readonly double KrogEndpoint = krogEndpoint;    // Krog at Sg = 0

        public readonly double Ng = ng;                        // Gas Corey Exponent
        public readonly double Nog = nog;                      // Oil-to-Gas Corey Exponent

        public readonly double Pdgo = pdgo;                    // Gas-oil displacement capillary entry
        public readonly double InvLambdaGo = 1.0 / lambdaGo;

        public double PcRSquared { get; } = pcRSquared;
        public double KrgRSquared { get; } = krgRSquared;
        public double KrogRSquared { get; } = krogRSquared;

        public double SgMin { get; } = sgMin;                  // first and last Sg of the table
        public double SgMax { get; } = sgMax;
        public double PcMin { get; } = pcMin;                  // smallest and largest Pcog of the table
        public double PcMax { get; } = pcMax;

        // Each curve over its own range, clamped at its ends. Before, Krg and Krog both used (Sg - Sgc) / (1 - Swi -
        // Sorg - Sgc) with no clamping, so Krog exceeded its end point at Sg < Sgc and a Sg below Sgc took a negative
        // number to a fractional power (NaN).
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrg<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            // Krg rises from Sgc to the first Sg where it reaches its end point, and is flat above that. It used to be
            // stretched to the table's last Sg: a table that reaches Krg = 1 at Sg = 0.85 and then has a row at Sg = 1
            // (SPE10 model 1) got a Krg 53 % too low at 0.85, so gas moved too slowly and oil rate was too high.
            ADiff<T> sgNorm = (sg - Sgc) / (SgKrgEnd - Sgc);
            if (sgNorm.Value <= 0.0) return default;
            if (sgNorm.Value >= 1.0) return new ADiff<T>(KrgEndpoint, default);
            return KrgEndpoint * Pow(sgNorm, Ng);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrog<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            // Oil (at connate water) flows from So = 1 - Swi (Sg = 0) down to its residual Sorg
            ADiff<T> soNorm = (1.0 - Swi - Sorg - sg) / (1.0 - Swi - Sorg);
            if (soNorm.Value <= 0.0) return default;
            if (soNorm.Value >= 1.0) return new ADiff<T>(KrogEndpoint, default);
            return KrogEndpoint * Pow(soNorm, Nog);
        }

        // Pcog grows with Sg and is singular at Sg = 1 - Swi; it is kept at or below the table's largest Pcog. A table
        // with no positive Pcog has Pdgo = 0, so Pcog = 0 (before, Sg near 1 - Swi gave an infinite Pcog and NaN pressures).
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcog<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            if (Pdgo <= 0.0 || PcMax <= 0.0) return default;
            ADiff<T> sgStar = sg / (1.0 - Swi);
            if (sgStar.Value >= 1.0 - 1e-6) return new ADiff<T>(PcMax, default);
            var pc = Pdgo * Pow(1.0 - sgStar, -InvLambdaGo);
            return pc.Value > PcMax ? new ADiff<T>(PcMax, default) : pc;
        }
    }

    public class SgofRow(double sg, double krg, double krog, double pcg, Units incomingUnit)
    {
        public double Sg { get; set; } = sg;
        public double Krg { get; set; } = krg;
        public double Krog { get; set; } = krog;
        public double Pcg { get; set; } = Utility.Convert(pcg, "Pressure", incomingUnit, Units.FIELD);
    }

    public class SGOFAnalysis(List<List<double>> _sgof, Units incomingUnit)
    {
        public readonly List<SgofRow> sgofTable = [.. _sgof.Select(sgof => new SgofRow(sgof[0], sgof[1], sgof[2], sgof[3], incomingUnit))];

        public SgofResult Result { get; private set; }

        public void FitProperties(double fixedSwi)
        {
            if (sgofTable == null || sgofTable.Count < 3)
                throw new ArgumentException("The SGOF table must contain at least 3 rows.");

            var sorted = sgofTable.OrderBy(r => r.Sg).ToList();

            // 1. Resolve Gas Endpoints
            // Critical gas saturation: the last Sg before Krg becomes positive (gas is immobile below it)
            var firstMobileGas = sorted.FindIndex(r => r.Krg > 0.0);
            double sgc = firstMobileGas > 0 ? sorted[firstMobileGas - 1].Sg : sorted.First().Sg;
            double sgLast = sorted.Last().Sg;
            double krogEndpoint = sorted.First().Krog; // Krog at Sg = 0
            // Krg end point: the table's largest Krg, and the first Sg at which it is reached (rows after that repeat it)
            double krgEndpoint = sorted.Max(r => r.Krg);
            double sgKrgEnd = sorted.First(r => r.Krg >= krgEndpoint * (1.0 - 1e-9)).Sg;
            if (sgKrgEnd <= sgc) sgKrgEnd = sgLast;    // Krg is 0 everywhere: keep a valid range

            // Residual oil to gas: Krog reaches 0 at Sg = 1 - Swi - Sorg (first such row)
            var activeGasLimit = sorted.FirstOrDefault(r => r.Krog <= 0.0) ?? sorted.Last();
            double sorg = Max(0.0, 1.0 - fixedSwi - activeGasLimit.Sg);

            List<double> pcX = []; List<double> pcY = [];
            List<double> krgX = []; List<double> krgY = [];
            List<double> krogX = []; List<double> krogY = [];

            foreach (var row in sorted)
            {
                if (row.Pcg > 0 && row.Sg < 1.0 - fixedSwi)
                {
                    double sgStar = row.Sg / (1.0 - fixedSwi);
                    pcX.Add(Log(1.0 - sgStar));
                    pcY.Add(Log(row.Pcg));
                }

                if (row.Sg > sgc && row.Sg < sgKrgEnd && row.Krg > 0)
                {
                    krgX.Add(Log((row.Sg - sgc) / (sgKrgEnd - sgc)));
                    krgY.Add(Log(row.Krg / krgEndpoint));
                }

                double so = 1.0 - fixedSwi - row.Sg;
                if (row.Sg > 0 && so > sorg && row.Krog > 0 && row.Sg < activeGasLimit.Sg)
                {
                    krogX.Add(Log((so - sorg) / (1.0 - fixedSwi - sorg)));
                    krogY.Add(Log(row.Krog));
                }
            }

            var (mPc, bPc, rSqPc) = Utility.LinearRegression(pcX, pcY);
            var (_, _, rSqG) = Utility.LinearRegression(krgX, krgY);
            double ng = CoreyExponent(krgX, krgY);
            var (nog, _, rSqOg) = Utility.LinearRegression(krogX, krogY);

            double lambdaGo = -1.0 / mPc;
            double pdgo = Exp(bPc);
            if (pcX.Count < 2) { pdgo = 0.0; lambdaGo = 1.0; }   // fewer than 2 positive Pcog rows: no capillary pressure

            Result = new(sgc, sorg, fixedSwi, krgEndpoint, krogEndpoint,
                ng, nog, pdgo, lambdaGo, rSqPc, rSqG, rSqOg,
                sorted.First().Sg, sgLast, sorted.Min(r => r.Pcg), sorted.Max(r => r.Pcg), sgKrgEnd);
        }

        // Corey exponent n of kr = end point x (normalised saturation)^n, from x = ln(normalised saturation) and
        // y = ln(kr / end point). The curve is made to pass through its end point, so the line is fitted through the
        // origin, and each row is weighted by kr squared, which is least squares on kr itself to first order. The plain
        // log-log regression gave rows such as Krg = 0.000001 (one significant digit) the same say as the rows where
        // gas flows, and fitted an intercept that was then dropped: SPE10 model 1's Krg came out 40-80 % too high.
        internal static double CoreyExponent(List<double> x, List<double> y)
        {
            double sxy = 0.0, sxx = 0.0;
            for (int i = 0; i < x.Count; i++)
            {
                double w = Exp(2.0 * y[i]);
                sxy += w * x[i] * y[i];
                sxx += w * x[i] * x[i];
            }
            return sxx > 0.0 ? sxy / sxx : 1.0;
        }

    }
}
