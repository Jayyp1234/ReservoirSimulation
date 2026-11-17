using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvdgResult(double pRef, double cg, double bgAtPref, double cMug,
                                      double viscGAtPref, double bgRSquared, double viscGRSquared)
    {
        public readonly double PRef = pRef;
        public readonly double Cg = cg;                // Gas compressibility factor
        public readonly double BgAtPref = bgAtPref;
        public readonly double CMug = cMug;            // Gas viscosity exponent
        public readonly double ViscGAtPref = viscGAtPref;

        public double BgRSquared { get; } = bgRSquared;
        public double ViscGRSquared { get; } = viscGRSquared;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBg<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return BgAtPref * Exp(-Cg * (p - PRef));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMug<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return ViscGAtPref * Exp(CMug * (p - PRef));
        }
    }

    public class PvdgRow(List<double> pvdg, Units incomingUnit)
    {
        public double Pressure { get; set; } = Utility.Convert(pvdg[0], "Pressure", incomingUnit, Units.FIELD);
        public double Bg { get; set; } = Utility.Convert(pvdg[1], "Gas Formation Volume Factor (Bg)", incomingUnit, Units.FIELD);
        public double Viscosity { get; set; } = Utility.Convert(pvdg[2], "Gas Viscosity", incomingUnit, Units.FIELD);
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

            List<double> deltaP = [];
            List<double> lnBg = [];
            List<double> lnViscG = [];

            foreach (var row in sorted)
            {
                if (row.Viscosity <= 0 || row.Bg <= 0)
                    throw new InvalidOperationException("Properties must be strictly positive.");

                deltaP.Add(row.Pressure - pRef);
                lnBg.Add(Log(row.Bg));
                lnViscG.Add(Log(row.Viscosity));
            }

            var (mBg, bBg, rSqBg) = Utility.LinearRegression(deltaP, lnBg);
            var (mViscG, bViscG, rSqViscG) = Utility.LinearRegression(deltaP, lnViscG);

            Result =  new(pRef, -mBg, Exp(bBg), mViscG, Exp(bViscG), rSqBg, rSqViscG);
        }
    }
}