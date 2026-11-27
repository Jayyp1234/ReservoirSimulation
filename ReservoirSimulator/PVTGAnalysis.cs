using SepalSolver;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvtgResult(double γg0,
        double rvBase, double rvExponent, double rvRSquared,
        double bgBase, double bgExponent, double bgRSquared,
        double viscBase, double viscSlope, double viscRSquared,
        double undersatBgSlope = 0.0, double undersatViscSlope = 0.0)
    {
        // Raw Coefficients held as raw fields (Blittable, ultra-lean)
        public readonly double γg0 = γg0;
        public readonly double RvBase = rvBase;
        public readonly double RvExponent = rvExponent;
        public readonly double RvRSquared = rvRSquared;

        public readonly double BgBase = bgBase;
        public readonly double BgExponent = bgExponent;
        public readonly double BgRSquared = bgRSquared;

        public readonly double ViscBase = viscBase;
        public readonly double ViscSlope = viscSlope;
        public readonly double ViscRSquared = viscRSquared;

        // Undersaturated sensitivity wrt (Rv_sat - Rv) at fixed pressure
        public readonly double UndersatBgSlope = undersatBgSlope;
        public readonly double UndersatViscSlope = undersatViscSlope;

        // High-Performance Inlined Mathematical Models (Saturated Dew-Point Curve)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateRv<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            if (p.Value <= 0.0 || RvBase <= 0.0)
                return default;

            return RvBase * Pow(p, RvExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateDewPointPressure<T>(ADiff<T> rv) where T : struct, IVectorND<T>
        {
            if (rv.Value <= 0.0 || RvBase <= 0.0 || Abs(RvExponent) < 1e-12)
                return default;

            return Pow(rv / RvBase, 1.0 / RvExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBg<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            if (p.Value <= 0.0 || BgBase <= 0.0)
                return default;

            return BgBase * Pow(p, -BgExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMug<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return ViscBase * Exp(ViscSlope * p);
        }

        // 2D Undersaturated Gas Extensions: (P, Rv) where Rv <= Rv_sat(P)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBg<T>(ADiff<T> p, ADiff<T> rv) where T : struct, IVectorND<T>
        {
            ADiff<T> bgSat = EvaluateBg(p);
            ADiff<T> rvDeficit = EvaluateRv(p) - rv; // 0 on dew-point curve, > 0 when undersaturated
            return bgSat * Exp(-UndersatBgSlope * rvDeficit);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMug<T>(ADiff<T> p, ADiff<T> rv) where T : struct, IVectorND<T>
        {
            ADiff<T> mugSat = EvaluateMug(p);
            ADiff<T> rvDeficit = EvaluateRv(p) - rv;
            return mugSat * Exp(-UndersatViscSlope * rvDeficit);
        }
    }

    public class PvtgRow(double pressure, double rv, double bg, double viscosity, Units incomingUnit)
    {
        public double Pressure { get; set; } = Utility.Convert(pressure, "Pressure", incomingUnit, Units.FIELD);
        public double Rv { get; set; } = Utility.Convert(rv, "Vaporized Oil-Gas Ratio (Rv)", incomingUnit, Units.FIELD);
        public double Bg { get; set; } = Utility.Convert(bg, "Gas Formation Volume Factor (Bg)", incomingUnit, Units.FIELD);
        public double Viscosity { get; set; } = Utility.Convert(viscosity, "Viscosity", incomingUnit, Units.FIELD);
    }

    public class PVTGAnalysis
    {
        private readonly List<PvtgRow> saturatedPoints = [];
        private readonly List<(PvtgRow sat, PvtgRow undersat)> undersaturatedPairs = [];

        public PVTGAnalysis(List<List<double>> _pvtg, Units incomingUnit)
        {
            for (int i = 0; i < _pvtg.Count; i++)
            {
                var pvtg = _pvtg[i];
                if (pvtg.Count < 4) continue;

                double p = pvtg[0];
                var chunks = pvtg[1..].Chunk(3).ToList();

                // First triplet [Rv_sat, Bg_sat, mu_g_sat] is ALWAYS the saturated dew-point state at pressure P
                var satChunk = chunks[0];
                var satRow = new PvtgRow(p, satChunk[0], satChunk[1], satChunk[2], incomingUnit);
                saturatedPoints.Add(satRow);

                // Remaining triplets [Rv_1, Bg_1, mu_g_1, ...] trace undersaturated gas (Rv < Rv_sat) at constant P
                for (int c = 1; c < chunks.Count; c++)
                {
                    var uChunk = chunks[c];
                    if (uChunk.Length < 3) continue;
                    var undersatRow = new PvtgRow(p, uChunk[0], uChunk[1], uChunk[2], incomingUnit);
                    undersaturatedPairs.Add((satRow, undersatRow));
                }
            }
        }

        public PvtgResult Result { get; private set; }

        public void FitProperties()
        {
            if (saturatedPoints.Count < 2)
            {
                throw new ArgumentException("PVTG analysis requires at least 2 distinct pressure entries to regress slopes.");
            }

            List<double> pressureValues = [];
            List<double> lnPressure = [];
            List<double> lnBg = [];
            List<double> lnVisc = [];

            // Synchronized lists for Rv = alpha * P^beta (excludes rows where Rv == 0)
            List<double> lnPressureForRv = [];
            List<double> lnRv = [];

            foreach (var row in saturatedPoints)
            {
                if (row.Pressure <= 0 || row.Bg <= 0 || row.Viscosity <= 0)
                    continue;

                double lnP = Log(row.Pressure);

                pressureValues.Add(row.Pressure);
                lnPressure.Add(lnP);
                lnBg.Add(Log(row.Bg));
                lnVisc.Add(Log(row.Viscosity));

                if (row.Rv > 0)
                {
                    lnPressureForRv.Add(lnP);
                    lnRv.Add(Log(row.Rv));
                }
            }

            // 1. Saturated Dew-Point Vaporized Oil-Gas Ratio: Rv = alpha * P^beta
            var (mRv, bRv, rSqRv) = lnPressureForRv.Count >= 2
                ? LinearRegression(lnPressureForRv, lnRv)
                : (0.0, 0.0, 1.0);

            // 2. Gas Formation Volume Factor: Bg = alpha * P^(-beta)  (Power-law expansion)
            var (mBg, bBg, rSqBg) = LinearRegression(lnPressure, lnBg);

            // 3. Gas Viscosity: mu_g = alpha * exp(beta * P)  (Exponential pressure thickening)
            var (mVisc, bVisc, rSqVisc) = LinearRegression(pressureValues, lnVisc);

            // 4. Optional Undersaturated Slopes wrt (Rv_sat - Rv) at fixed Pressure
            double undersatBgSlope = 0.0;
            double undersatViscSlope = 0.0;

            if (undersaturatedPairs.Count > 0)
            {
                List<double> dRv = [];
                List<double> dLnBg = [];
                List<double> dLnVisc = [];

                foreach (var (sat, undersat) in undersaturatedPairs)
                {
                    double deltaRv = sat.Rv - undersat.Rv;
                    if (deltaRv > 1e-12 && sat.Bg > 0 && undersat.Bg > 0 && sat.Viscosity > 0 && undersat.Viscosity > 0)
                    {
                        dRv.Add(deltaRv);
                        dLnBg.Add(-Log(undersat.Bg / sat.Bg));
                        dLnVisc.Add(-Log(undersat.Viscosity / sat.Viscosity));
                    }
                }

                if (dRv.Count > 0)
                {
                    double sumXX = dRv.Sum(Sqr);
                    if (sumXX > 1e-16)
                    {
                        undersatBgSlope = dRv.Zip(dLnBg, (x, y) => x * y).Sum() / sumXX;
                        undersatViscSlope = dRv.Zip(dLnVisc, (x, y) => x * y).Sum() / sumXX;
                    }
                }
            }

            Result = new PvtgResult(
                lnPressureForRv.Count >= 2 ? Exp(bRv) : 0.0, mRv, rSqRv,
                Exp(bBg), -mBg, rSqBg,
                Exp(bVisc), mVisc, rSqVisc,
                undersatBgSlope, undersatViscSlope);
        }

        private static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            double sumX = x.Sum();
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Sum(Sqr);

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (Abs(denominator) < 1e-11) return (0.0, 0.0, 1.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            double yMean = sumY / n;
            double ssTot = y.Select(yi => Sqr(yi - yMean)).Sum();
            double ssRes = y.Zip(x, (yi, xi) => Sqr(yi - (m * xi + b))).Sum();
            double rSq = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return (m, b, rSq);
        }
    }
}