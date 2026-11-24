using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvdgResult(double pRef, double bgBase, double bgExponent,
    double cMug, double viscGAtPref, double rSqBg, double rSqViscG)
    {
        public readonly double PRef = pRef;
        public readonly double BgBase = bgBase;
        public readonly double BgExponent = bgExponent;
        public readonly double CMug = cMug;
        public readonly double ViscGAtPref = viscGAtPref;
        public readonly double RSqBg = rSqBg;
        public readonly double RSqViscG = rSqViscG;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBg<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            ADiff<T> pSafe = Max(14.7, p); // Guard against P <= 0 during Newton trial steps
            return BgBase * Pow(pSafe, BgExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double EvaluateBg(double p)
        {
            double pSafe = Max(14.7, p);
            return BgBase * Pow(pSafe, BgExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMug<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return ViscGAtPref * Exp(CMug * (p - PRef));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double EvaluateMug(double p)
        {
            return ViscGAtPref * Exp(CMug * (p - PRef));
        }
    }

    public class PvdgRow(List<double> pvdg, Units incomingUnit)
    {
        public double Pressure { get; set; } = Utility.Convert(pvdg[0], "Pressure", incomingUnit, Units.FIELD);
        public double Bg { get; set; } = Utility.Convert(pvdg[1], "Gas Formation Volume Factor (Bg)", incomingUnit, Units.FIELD);
        public double Viscosity { get; set; } = Utility.Convert(pvdg[2], "Viscosity", incomingUnit, Units.FIELD);
    }

    public class PVDGAnalysis(List<List<double>> _pvdg, Units incomingUnit)
    {
        private readonly List<PvdgRow> pvdgTable = [.. _pvdg.Select(pvdg => new PvdgRow(pvdg, incomingUnit))];

        public PvdgResult Result { get; private set; }

        public void FitProperties()
        {
            if (pvdgTable == null || pvdgTable.Count < 2)
                throw new ArgumentException("PVDG table requires at least 2 rows.");

            var sorted = pvdgTable.OrderBy(r => r.Pressure).ToList();
            double pRef = sorted[0].Pressure;

            List<double> lnP = [];
            List<double> deltaP = [];
            List<double> lnBg = [];
            List<double> lnViscG = [];

            foreach (var row in sorted)
            {
                if (row.Pressure <= 0 || row.Viscosity <= 0 || row.Bg <= 0)
                    throw new InvalidOperationException("Pressure, Bg, and Viscosity must be strictly positive.");

                lnP.Add(Log(row.Pressure));
                deltaP.Add(row.Pressure - pRef);
                lnBg.Add(Log(row.Bg));
                lnViscG.Add(Log(row.Viscosity));
            }

            // Power-law fit for Bg: ln(Bg) = ln(BgBase) + BgExponent * ln(P)
            var (mBg, bBg, rSqBg) = Utility.LinearRegression(lnP, lnBg);

            // Exponential fit for gas viscosity: ln(Mug) = ln(ViscGAtPref) + cMug * (P - pRef)
            var (mViscG, bViscG, rSqViscG) = Utility.LinearRegression(deltaP, lnViscG);

            Result = new(pRef, Exp(bBg), mBg, mViscG, Exp(bViscG), rSqBg, rSqViscG);
        }
    }
}