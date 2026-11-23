using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvdoResult(double pRef, double co, double boAtPref, double cMu,
                      double viscAtPref, double boRSquared, double viscRSquared)
    {
        // Raw Regression Coefficients & Diagnostic Parameters
        public readonly double PRef = pRef;
        public readonly double Co = co;
        public readonly double BoAtPref = boAtPref;
        public readonly double CMu = cMu;
        public readonly double ViscAtPref = viscAtPref;

        public double BoRSquared { get; } = boRSquared;
        public double ViscRSquared { get; } = viscRSquared;

        // High-Performance Inlined Undersaturated Models (0% Delegate Overhead)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBo<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return BoAtPref * Exp(-Co * (p - PRef));
        }
        public double EvaluateBo(double p)
        {
            return BoAtPref * Exp(-Co * (p - PRef));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMuo<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return ViscAtPref * Exp(CMu * (p - PRef));
        }
    }


    public class PvdoRow(List<double> pvdo, Units incomingUnit)
    {
        public double Pressure { get; set; } = Utility.Convert(pvdo[0], "Pressure", incomingUnit, Units.FIELD);
        public double Bo { get; set; } = Utility.Convert(pvdo[1], "Oil Formation Volume Factor (Bo)", incomingUnit, Units.FIELD);
        public double Viscosity { get; set; } = Utility.Convert(pvdo[2], "Viscosity", incomingUnit, Units.FIELD);
    }

    public class PVDOAnalysis(List<List<double>> _pvdo, Units incomingUnit)
    {
        private readonly List<PvdoRow> pvdoTable = [.._pvdo.Select(pvdo => new PvdoRow(pvdo, incomingUnit))];

        public PvdoResult Result { get; private set; }

        public void FitProperties()
        {
            if (pvdoTable == null || pvdoTable.Count < 2)
            {
                throw new ArgumentException("PVDO table requires at least 2 rows to establish exponential parameters.");
            }

            // 1. Sort ascending by pressure
            var sorted = pvdoTable.OrderBy(r => r.Pressure).ToList();

            // 2. Automatically discover the reference anchor pressure (Bubble Point)
            double pRef = sorted[0].Pressure;

            // Coordinate transformations
            List<double> deltaP = []; // x = (P - PRef)
            List<double> lnBo = [];   // y1 = ln(Bo)
            List<double> lnVisc = []; // y2 = ln(mu_o)

            foreach (var row in sorted)
            {
                if (row.Viscosity <= 0 || row.Bo <= 0)
                    throw new InvalidOperationException("PVT properties must be strictly greater than zero for log conversion.");

                deltaP.Add(row.Pressure - pRef);
                lnBo.Add(Log(row.Bo));
                lnVisc.Add(Log(row.Viscosity));
            }

            // Run regressions in the transformed log-linear spaces
            var (mBo, bBo, rSqBo) = Utility.LinearRegression(deltaP, lnBo);
            var (mVisc, bVisc, rSqVisc) = Utility.LinearRegression(deltaP, lnVisc);

            // Back-transform intercepts to resolve physical units
            double co = -mBo;                      // Bo decreases with P -> co is -slope
            double boAtPref = Exp(bBo);            // e^intercept_Bo
            double cMu = mVisc;                    // Viscosity increases with P -> c_mu is slope
            double viscAtPref = Exp(bVisc);        // e^intercept_Visc

            Result = new(pRef, co, boAtPref, cMu, viscAtPref, rSqBo, rSqVisc);
        }
    }
}