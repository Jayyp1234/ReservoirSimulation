using SepalSolver;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvtoResult(
        double rsBase, double rsExponent, double rsRSquared,
        double boBase, double boExponent, double boRSquared,
        double viscBase, double viscExponent, double viscRSquared,
        double undersatCompressibility = 0.0, double undersatViscosibility = 0.0)
    {
        // Raw Coefficients held as raw fields (Blittable, ultra-lean)
        public readonly double RsBase = rsBase;
        public readonly double RsExponent = rsExponent;
        public readonly double RsRSquared = rsRSquared;

        public readonly double BoBase = boBase;
        public readonly double BoExponent = boExponent;
        public readonly double BoRSquared = boRSquared;

        public readonly double ViscBase = viscBase;
        public readonly double ViscExponent = viscExponent;
        public readonly double ViscRSquared = viscRSquared;

        // Undersaturated sensitivities wrt (P - Pbub) at fixed Rs
        public readonly double UndersatCompressibility = undersatCompressibility;
        public readonly double UndersatViscosibility = undersatViscosibility;

        // High-Performance Inlined Mathematical Models (Saturated Bubble-Point Curve)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateRs<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return RsBase * Pow(p, RsExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBubblePointPressure<T>(ADiff<T> rs) where T : struct, IVectorND<T>
        {
            return Pow(rs / RsBase, 1.0 / RsExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBo<T>(ADiff<T> rs) where T : struct, IVectorND<T>
        {
            return BoBase * Exp(BoExponent * rs);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMuo<T>(ADiff<T> rs) where T : struct, IVectorND<T>
        {
            return ViscBase * Exp(-ViscExponent * rs);
        }

        // 2D Undersaturated Oil Extensions: (P, Rs) where P >= Pbub(Rs)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBo<T>(ADiff<T> p, ADiff<T> rs) where T : struct, IVectorND<T>
        {
            ADiff<T> boSat = EvaluateBo(rs);
            ADiff<T> dpOver = p - EvaluateBubblePointPressure(rs); // 0 on bubble-point curve, > 0 when undersaturated
            return boSat * Exp(-UndersatCompressibility * dpOver);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMuo<T>(ADiff<T> p, ADiff<T> rs) where T : struct, IVectorND<T>
        {
            ADiff<T> muSat = EvaluateMuo(rs);
            ADiff<T> dpOver = p - EvaluateBubblePointPressure(rs);
            return muSat * Exp(UndersatViscosibility * dpOver);
        }
    }

    public class PvtoRow(double rs, double pressure, double bo, double viscosity, Units incomingUnit)
    {
        public double Rs { get; set; } = Utility.Convert(rs, "Solution Gas-Oil Ratio (Rs)", incomingUnit, Units.FIELD);
        public double Pressure { get; set; } = Utility.Convert(pressure, "Pressure", incomingUnit, Units.FIELD);
        public double Bo { get; set; } = Utility.Convert(bo, "Oil Formation Volume Factor (Bo)", incomingUnit, Units.FIELD);
        public double Viscosity { get; set; } = Utility.Convert(viscosity, "Viscosity", incomingUnit, Units.FIELD);
    }

    public class PVTOAnalysis
    {
        private readonly List<PvtoRow> saturatedPoints = [];
        private readonly List<(PvtoRow sat, PvtoRow undersat)> undersaturatedPairs = [];

        public PVTOAnalysis(List<List<double>> _pvto, Units incomingUnit)
        {
            for (int i = 0; i < _pvto.Count; i++)
            {
                var pvto = _pvto[i];
                if (pvto.Count < 4) continue;

                double rs = pvto[0];
                var chunks = pvto[1..].Chunk(3).ToList();

                // First triplet [P_bub, Bo_sat, mu_o_sat] is ALWAYS the saturated bubble-point state for this Rs
                var satChunk = chunks[0];
                var satRow = new PvtoRow(rs, satChunk[0], satChunk[1], satChunk[2], incomingUnit);
                saturatedPoints.Add(satRow);

                // Remaining triplets [P_1, Bo_1, mu_o_1, ...] trace undersaturated oil (P > P_bub) at constant Rs
                for (int c = 1; c < chunks.Count; c++)
                {
                    var uChunk = chunks[c];
                    if (uChunk.Length < 3) continue;
                    var undersatRow = new PvtoRow(rs, uChunk[0], uChunk[1], uChunk[2], incomingUnit);
                    undersaturatedPairs.Add((satRow, undersatRow));
                }
            }
        }

        public PvtoResult Result { get; private set; }

        public void FitProperties()
        {
            if (saturatedPoints.Count < 2)
            {
                throw new ArgumentException("PVTO analysis requires at least 2 distinct Rs saturation entries to regress slopes.");
            }

            // Synchronized lists for Bo and Viscosity regressions in log-linear space (includes Rs == 0 anchor if present)
            List<double> rsValues = [];
            List<double> lnBo = [];
            List<double> lnVisc = [];

            // Synchronized lists for Rs = alpha * P^beta in log-log space (requires P > 0 and Rs > 0)
            List<double> lnPressureForRs = [];
            List<double> lnRs = [];

            foreach (var row in saturatedPoints)
            {
                if (row.Viscosity <= 0 || row.Bo <= 0 || row.Rs < 0)
                    continue;

                rsValues.Add(row.Rs);
                lnBo.Add(Log(row.Bo));
                lnVisc.Add(Log(row.Viscosity));

                if (row.Pressure > 0 && row.Rs > 0)
                {
                    lnPressureForRs.Add(Log(row.Pressure));
                    lnRs.Add(Log(row.Rs));
                }
            }

            // 1. Run regressions for fluid properties driven by Rs
            var (mBo, bBo, rSqBo) = Utility.LinearRegression(rsValues, lnBo);
            var (mVisc, bVisc, rSqVisc) = Utility.LinearRegression(rsValues, lnVisc);

            // 2. Run log-log regression for Rs driven by Pressure along the saturation boundary
            var (mRs, bRs, rSqRs) = lnPressureForRs.Count >= 2
                ? Utility.LinearRegression(lnPressureForRs, lnRs)
                : (0.0, 0.0, 1.0);

            // 3. Regress undersaturated oil compressibility (co) and viscosibility (cmu) from P > Pbub branches
            double undersatCompressibility = 0.0;
            double undersatViscosibility = 0.0;

            if (undersaturatedPairs.Count > 0)
            {
                List<double> dP = [];
                List<double> dLnBo = [];
                List<double> dLnVisc = [];

                foreach (var (sat, undersat) in undersaturatedPairs)
                {
                    double deltaP = undersat.Pressure - sat.Pressure;
                    if (deltaP > 1e-12 && sat.Bo > 0 && undersat.Bo > 0 && sat.Viscosity > 0 && undersat.Viscosity > 0)
                    {
                        dP.Add(deltaP);
                        dLnBo.Add(-Log(undersat.Bo / sat.Bo));
                        dLnVisc.Add(Log(undersat.Viscosity / sat.Viscosity));
                    }
                }

                if (dP.Count > 0)
                {
                    double sumXX = dP.Sum(Sqr);
                    if (sumXX > 1e-16)
                    {
                        undersatCompressibility = dP.Zip(dLnBo, (x, y) => x * y).Sum() / sumXX;
                        undersatViscosibility = dP.Zip(dLnVisc, (x, y) => x * y).Sum() / sumXX;
                    }
                }
            }

            Result = new PvtoResult(
                lnPressureForRs.Count >= 2 ? Exp(bRs) : 0.0, mRs, rSqRs,
                Exp(bBo), mBo, rSqBo,
                Exp(bVisc), -mVisc, rSqVisc,
                undersatCompressibility, undersatViscosibility);
        }
    }
}