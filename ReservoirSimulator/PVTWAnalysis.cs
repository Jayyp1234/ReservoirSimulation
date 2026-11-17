using SepalSolver;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;
using static ReservoirSimulator.PVTWAnalysis;

namespace ReservoirSimulator
{
    public readonly struct PvtwResult(double Pref, double Bwref, double cw, double Muref,
                                          double vw, double bwRSquared, double viscGRSquared)
    {
        public readonly double PRef = Pref;
        public readonly double Cw = cw;
        public readonly double BwAtPref = Bwref;
        public readonly double CMu = vw;
        public readonly double ViscAtPref = Muref;
        public double BwRSquared { get; } = bwRSquared;
        public double ViscGRSquared { get; } = viscGRSquared;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBw<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return Bwref * Exp(-cw * (p - Pref));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMuw<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return Muref * Exp(vw * (p - Pref));
        }
    }

    public class PvtwRow(double pref, double bwref, double cw, double visref, double vw, Units incomingUnit)
    {
        public double Pref { get; set; } = Utility.Convert(pref, "Pressure", incomingUnit, Units.FIELD);
        public double Bwref { get; set; } = Utility.Convert(bwref, "Water Formation Volume Factor (Bw)", incomingUnit, Units.FIELD);
        public double Cw { get; set; } = Utility.Convert(cw, "Compressibility", incomingUnit, Units.FIELD);
        public double Visref { get; set; } = Utility.Convert(visref, "Viscosity", incomingUnit, Units.FIELD);
        public double Vw { get; set; } = Utility.Convert(vw, "Compressibility", incomingUnit, Units.FIELD);
    }

    public class PVTWAnalysis(List<double> _pvtw, Units incomingUnit)
    {
        private readonly PvtwRow pvtwRow = new(_pvtw[0], _pvtw[1], _pvtw[2], _pvtw[3], _pvtw[4], incomingUnit);

        public PvtwResult Result { get; private set; }
        public void FitProperties()
        {
            Result = new(pvtwRow.Pref, pvtwRow.Bwref, pvtwRow.Cw, pvtwRow.Visref, pvtwRow.Vw, 1, 1);
        }
    }
}