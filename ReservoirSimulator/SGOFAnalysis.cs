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
        public readonly double Sgc = sgc;                      // Critical Gas Saturation
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

        public double SgMin { get; } = sgMin;
        public double SgMax { get; } = sgMax;
        public double PcMin { get; } = pcMin;
        public double PcMax { get; } = pcMax;


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
            ADiff<T> sgNorm = (sg - Sgc) / (1.0 - Swi - Sorg - Sgc);
            return KrogEndpoint * Pow(1.0 - sgNorm, Nog);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcog<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            // Gas-oil capillary expansion operates via scaled gas track projection
            ADiff<T> sgStar = sg / (1.0 - Swi);
            return Pdgo * Pow(1.0 - sgStar, -InvLambdaGo);
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
        public readonly List<SgofRow> sgofTable = [.._sgof.Select(sgof => new SgofRow(sgof[0], sgof[1], sgof[2], sgof[3], incomingUnit))];

        public SgofResult Result { get; private set; }

        public void FitProperties(double fixedSwi)
        {
            if (sgofTable == null || sgofTable.Count < 3)
                throw new ArgumentException("The SGOF table must contain at least 3 rows.");

            var sorted = sgofTable.OrderBy(r => r.Sg).ToList();

            // 1. Resolve Gas Endpoints
            double sgc = sorted.First().Sg; // Critical gas initialization entry
            double krogEndpoint = sorted.First().Krog; // Krog at Sg = 0
            // Krg end point: the table's largest Krg, and the first Sg at which it is reached (rows after that repeat it)
            double krgEndpoint = sorted.Max(r => r.Krg);
            double sgKrgEnd = sorted.First(r => r.Krg >= krgEndpoint * (1.0 - 1e-9)).Sg;
            if (sgKrgEnd <= sgc) sgKrgEnd = sgLast;    // Krg is 0 everywhere: keep a valid range

            var activeGasLimit = sorted.FirstOrDefault(r => r.Krog <= 0.0) ?? sorted.Last();
            double sgMax = activeGasLimit.Sg;
            double sorg = 1.0 - fixedSwi - sgMax;
            double krgEndpoint = activeGasLimit.Krg;

            List<double> pcX = []; List<double> pcY = [];
            List<double> krgX = []; List<double> krgY = [];
            List<double> krogX = []; List<double> krogY = [];

            foreach (var row in sorted)
            {
                if (row.Sg > sgc && row.Pcg > 0)
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
                    double sgNorm = (row.Sg - sgc) / (1.0 - fixedSwi - sorg - sgc);
                    if (row.Krg > 0)
                    {
                        krgX.Add(Log(sgNorm));
                        krgY.Add(Log(row.Krg));
                    }
                    if (row.Krog > 0)
                    {
                        krogX.Add(Log(1.0 - sgNorm));
                        krogY.Add(Log(row.Krog));
                    }
                }
            }

            var (mPc, bPc, rSqPc) = Utility.LinearRegression(pcX, pcY);
            var (_, _, rSqG) = Utility.LinearRegression(krgX, krgY);
            double ng = CoreyExponent(krgX, krgY);
            var (nog, _, rSqOg) = Utility.LinearRegression(krogX, krogY);

            double lambdaGo = -1.0 / mPc;
            double pdgo = Exp(bPc);

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