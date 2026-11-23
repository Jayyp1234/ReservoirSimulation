using EclipseDeckReader;
using SepalSolver;
using System;
using System.Collections.Generic;
using System.Linq;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Well<T> where T : struct, IVectorND<T>
    {
        public Offsets Offset { get; set; }

        public string Name { get; set; }
        public Phase PreferedPhase { get; set; }

        // Static list of all active reservoir grid cells this well ever connects to (sorted by natural index m)
        public List<int> PerfNatIndex { get; set; } = [];
        public List<bool> PerfIsOpen { get; set; } = [];

        public ADiff<T> Qo { get; set; }
        public ADiff<T> Qw { get; set; }
        public ADiff<T> Ql => Qo + Qw;
        public ADiff<T> Qg { get; set; }
        public ADiff<T> Pwf { get; set; }
        public double Fw => (Qo + Qw).Value != 0.0 ? Qw.Value * 100.0 / (Qo + Qw).Value : 0.0;

        public int I { get; }
        public int J { get; }
        public int Nx { get; set; }
        public int Ny { get; set; }
        public int Nz { get; set; }
        public Grid<T>[] Grids { get; set; }
        public int ref_Index { get; set; }

        public ConstraintType ConstraintType { get; set; } = ConstraintType.LRATE;
        public WellType WellType { get; set; }

        public ProducerControlRecord? CurrentProducerControl { get; private set; }
        public InjectorControlRecord? CurrentInjectorControl { get; private set; }
        public string CurrentWellStatus { get; private set; } = "SHUT";

        public bool isActive => !string.Equals(CurrentWellStatus, "SHUT", StringComparison.OrdinalIgnoreCase)
                             && (CurrentProducerControl is not null || CurrentInjectorControl is not null);
        public bool isStopped => string.Equals(CurrentWellStatus, "STOP", StringComparison.OrdinalIgnoreCase);
        public double Zref { get; set; }
        public double HighestZref { get; set; }
        public double LowestZref { get; set; }
        public double HighestZPressure { get; set; }
        public double LowestZPressure { get; set; }
        public double Pwf_n, Qo_n, Qw_n, Qg_n, Ql_n;
        public double Pwf_np1, Qo_np1, Qw_np1, Qg_np1, Ql_np1;

        public List<double> PWF = [], QO = [], QW = [], QG = [];

        private readonly WellData _wellData;
        private readonly bool _hasExplicitZref;
        private readonly double _alphaWell = 1.127e-3 * 2.0 * pi;
        private readonly Dictionary<int, int> _natIndexToPerfSlot = [];
        public bool Control_is_Violated = false;

        public Well(WellData welldata, int nx, int ny, int nz, Grid<T>[] grids, Units incomingunit)
        {
            _wellData = welldata;
            Name = _wellData.Specification.WellName;

            double specZref = _wellData.Specification.ReferenceDepth;
            _hasExplicitZref = !double.IsNaN(specZref) && specZref > 0.0;
            Zref = _hasExplicitZref
                ? Utility.Convert(specZref, "Length", incomingunit, Units.FIELD)
                : double.NaN;

            I = _wellData.Specification.HeadI - 1; // Convert 1-based ECLIPSE head index to 0-based
            J = _wellData.Specification.HeadJ - 1;
            Nx = nx; Ny = ny; Nz = nz; Grids = grids;

            PreferedPhase = _wellData.Specification.PreferredPhase.ToUpperInvariant() switch
            {
                "OIL" => Phase.OIL,
                "WATER" => Phase.WATER,
                "GAS" => Phase.GAS,
                _ => Phase.OIL
            };

            foreach (var prod in _wellData.ProducerControlHistory)
            {
                prod.OilRateTarget = Utility.Convert(prod.OilRateTarget, "Oil Surface Rate", incomingunit, Units.FIELD);
                prod.WaterRateTarget = Utility.Convert(prod.WaterRateTarget, "Water Surface Rate", incomingunit, Units.FIELD);
                prod.LiquidRateTarget = Utility.Convert(prod.LiquidRateTarget, "Liquid Surface Rate", incomingunit, Units.FIELD);
                prod.ResvRateTarget = Utility.Convert(prod.ResvRateTarget, "Reservoir Volume Rate", incomingunit, Units.FIELD);
                prod.GasRateTarget = Utility.Convert(prod.GasRateTarget, "Gas Surface Rate", incomingunit, Units.FIELD);
                prod.MinBhpLimit = Utility.Convert(prod.MinBhpLimit, "Pressure", incomingunit, Units.FIELD);
            }

            foreach (var inj in _wellData.InjectorControlHistory)
            {
                bool isGas = string.Equals(inj.InjectionPhase?.Trim(), "GAS", StringComparison.OrdinalIgnoreCase);
                inj.TargetRate = -Utility.Convert(inj.TargetRate, isGas ? "Gas Surface Rate" : "Liquid Surface Rate", incomingunit, Units.FIELD);
                inj.ResvRateTarget = -Utility.Convert(inj.ResvRateTarget, "Reservoir Volume Rate", incomingunit, Units.FIELD);
                inj.MaxBhpLimit = Utility.Convert(inj.MaxBhpLimit, "Pressure", incomingunit, Units.FIELD);
            }

            // Pre-discover all active grid cells the well ever connects to across CompletionHistory
            var allConnectedCells = new SortedDictionary<int, WellCompletionRecord>();
            foreach (var comp in _wellData.CompletionHistory)
            {
                for (int k = comp.KStart; k <= comp.KEnd; k++)
                {
                    int kZero = k - 1;
                    int m = comp.I + comp.J * Nx + kZero * Nx * Ny;
                    if (m < 0 || m >= Grids.Length) continue;
                    if (!Grids[m].isActive) continue;
                    allConnectedCells.TryAdd(m, comp);
                }
            }

            foreach (var kvp in allConnectedCells)
            {
                int m = kvp.Key;
                var comp = kvp.Value;
                int slot = PerfNatIndex.Count;

                _natIndexToPerfSlot[m] = slot;
                PerfNatIndex.Add(m);
                PerfIsOpen.Add(false);

                Grids[m].WI = ComputePeacemanWI(Grids[m], comp);
            }
        }

        private double ComputePeacemanWI(Grid<T> grid, WellCompletionRecord comp)
        {
            double DX = 0, DY = 0, DZ = 0;
            for (int c = 0; c < 8; c++)
            {
                double fx = (c == 1 || c == 2 || c == 5 || c == 6) ? 1 : -1;
                double fy = (c == 2 || c == 3 || c == 6 || c == 7) ? 1 : -1;
                double fz = (c >= 4) ? 1 : -1;
                DX += fx * grid.Corners[c].X;
                DY += fy * grid.Corners[c].Y;
                DZ += fz * grid.Corners[c].Z;
            }
            DX /= 4.0;
            DY /= 4.0;
            DZ /= 4.0;

            double radius = comp.WellBoreDiameter > 0.0 ? comp.WellBoreDiameter * 0.5 : 0.25;
            double skin = comp.EffectiveSkin;
            double ky = grid.Faces[1].Permeability;
            double kx = grid.Faces[2].Permeability;

            double re = 0.28 * Hypot(Pow(ky / kx, 0.25) * DX, Pow(kx / ky, 0.25) * DY) /
                               (Pow(ky / kx, 0.25) + Pow(kx / ky, 0.25));

            return _alphaWell * Sqrt(kx * ky) * DZ / (Log(re / radius) + skin);
        }

        public void SaveResult(double? pwf = null, double? qw = null,
            double? qg = null, double? qo = null, double? sg = null)
        {
            Utility.AddIfPresent(PWF, pwf);
            Utility.AddIfPresent(QG, qg);
            Utility.AddIfPresent(QW, qw);
            Utility.AddIfPresent(QO, qo);
        }

        internal List<DateTime> GetEventTimes()
        {
            return
            [
                .. _wellData.ProducerControlHistory.Select(h => h.EffectiveDate),
                .. _wellData.InjectorControlHistory.Select(h => h.EffectiveDate),
                .. _wellData.CompletionHistory.Select(h => h.EffectiveDate),
            ];
        }

        internal void UpdateControl(DateTime time)
        {
            ProducerControlRecord? latestProd = null;
            for (int m = 0; m < _wellData.ProducerControlHistory.Count; m++)
            {
                if (_wellData.ProducerControlHistory[m].EffectiveDate <= time)
                    latestProd = _wellData.ProducerControlHistory[m];
                else
                    break;
            }

            InjectorControlRecord? latestInj = null;
            for (int m = 0; m < _wellData.InjectorControlHistory.Count; m++)
            {
                if (_wellData.InjectorControlHistory[m].EffectiveDate <= time)
                    latestInj = _wellData.InjectorControlHistory[m];
                else
                    break;
            }

            if (latestProd is null && latestInj is null)
            {
                CurrentProducerControl = null;
                CurrentInjectorControl = null;
                CurrentWellStatus = "SHUT";
                return;
            }

            // Pick whichever control record (WCONPROD vs WCONINJE) took effect most recently
            bool useInjector = latestInj is not null &&
                               (latestProd is null || latestInj.EffectiveDate >= latestProd.EffectiveDate);

            if (useInjector)
            {
                // Ignore gas injection schedules in 2-phase oil/water mode
                if (string.Equals(latestInj!.InjectionPhase?.Trim(), "GAS", StringComparison.OrdinalIgnoreCase))
                {
                    CurrentProducerControl = null;
                    CurrentInjectorControl = null;
                    CurrentWellStatus = "SHUT";
                    return;
                }

                CurrentInjectorControl = latestInj;
                CurrentProducerControl = null;
                WellType = WellType.INJECTOR;
                CurrentWellStatus = latestInj.Status;

                ConstraintType = latestInj.ControlMode.ToUpperInvariant() switch
                {
                    "RATE" => ConstraintType.RATE,
                    "BHP" => ConstraintType.MAX_BHP,
                    _ => ConstraintType.RATE
                };
            }
            else
            {
                CurrentProducerControl = latestProd!;
                CurrentInjectorControl = null;
                WellType = WellType.PRODUCER;
                CurrentWellStatus = latestProd!.Status;

                ConstraintType = latestProd.ControlMode.ToUpperInvariant() switch
                {
                    "ORAT" => ConstraintType.ORATE,
                    "WRAT" => ConstraintType.WRATE,
                    "GRAT" => ConstraintType.GRATE,
                    "LRAT" => ConstraintType.LRATE,
                    "BHP" => ConstraintType.MIN_BHP,
                    _ => ConstraintType.ORATE
                };
            }
        }

        internal void UpdateCompletion(DateTime time)
        {
            for (int p = 0; p < PerfIsOpen.Count; p++)
                PerfIsOpen[p] = false;

            if (!isActive)
                return;

            var history = _wellData.CompletionHistory;
            var latestCompPerSlot = new WellCompletionRecord?[PerfNatIndex.Count];

            // Replay all COMPDAT records up to current time; later records override earlier ones
            for (int idx = 0; idx < history.Count; idx++)
            {
                var comp = history[idx];
                if (comp.EffectiveDate > time)
                    break;

                bool open = string.Equals(comp.Status, "OPEN", StringComparison.OrdinalIgnoreCase);

                for (int k = comp.KStart; k <= comp.KEnd; k++)
                {
                    int kZero = k - 1;
                    int m = comp.I + comp.J * Nx + kZero * Nx * Ny;
                    if (_natIndexToPerfSlot.TryGetValue(m, out int slot))
                    {
                        PerfIsOpen[slot] = open;
                        latestCompPerSlot[slot] = open ? comp : null;
                    }
                }
            }

            double shallowestZ = double.MaxValue;
            double deepestZ = double.MinValue;
            double shallowestP = Pwf_n;
            double deepestP = Pwf_n;
            int openCount = 0;

            for (int p = 0; p < PerfNatIndex.Count; p++)
            {
                if (!PerfIsOpen[p]) continue;

                openCount++;
                int m = PerfNatIndex[p];
                var grid = Grids[m];

                if (latestCompPerSlot[p] is { } comp)
                    grid.WI = ComputePeacemanWI(grid, comp);

                if (grid.Center.Z < shallowestZ)
                {
                    shallowestZ = grid.Center.Z;
                    shallowestP = grid.Po_n;
                }
                if (grid.Center.Z > deepestZ)
                {
                    deepestZ = grid.Center.Z;
                    deepestP = grid.Po_n;
                }
            }

            if (openCount > 0)
            {
                HighestZref = shallowestZ;
                LowestZref = deepestZ;
                HighestZPressure = shallowestP;
                LowestZPressure = deepestP;

                if (WellType == WellType.PRODUCER)
                {
                    if (!_hasExplicitZref)
                        Zref = HighestZref;
                    if (Pwf_n <= 0.0)
                        Pwf_n = HighestZPressure;
                }
                else
                {
                    if (!_hasExplicitZref)
                        Zref = LowestZref;
                    if (Pwf_n <= 0.0)
                        Pwf_n = LowestZPressure;
                }
            }
        }

        internal bool IsWaterTarget =>
            ConstraintType == ConstraintType.WRATE ||
            (ConstraintType == ConstraintType.RATE &&
            (CurrentInjectorControl is not null
                 ? string.Equals(CurrentInjectorControl.InjectionPhase, "WATER", StringComparison.OrdinalIgnoreCase)
                 : PreferedPhase == Phase.WATER));

        internal bool IsOilTarget =>
            ConstraintType == ConstraintType.ORATE ||
            (ConstraintType == ConstraintType.RATE && PreferedPhase == Phase.OIL);

        internal bool IsGasTarget =>
            ConstraintType == ConstraintType.GRATE ||
            (ConstraintType == ConstraintType.RATE && PreferedPhase == Phase.GAS);

        internal ADiff<T> Constraint(ADiff<T> pressure, ADiff<T> rate)
        {
            // Surface valve shut ("STOP"): enforce zero net surface rate while tracking static Pwf
            if (isStopped)
                return rate;

            if (WellType == WellType.INJECTOR && CurrentInjectorControl is not null)
            {
                return ConstraintType switch
                {
                    ConstraintType.MAX_BHP => pressure - CurrentInjectorControl.MaxBhpLimit,
                    ConstraintType.RATE => rate - CurrentInjectorControl.TargetRate,
                    _ => throw new NotSupportedException($"Unsupported injector control mode: {CurrentInjectorControl.ControlMode}")
                };
            }

            if (WellType == WellType.PRODUCER && CurrentProducerControl is not null)
            {
                return ConstraintType switch
                {
                    ConstraintType.MIN_BHP => pressure - CurrentProducerControl.MinBhpLimit,
                    ConstraintType.ORATE => rate - CurrentProducerControl.OilRateTarget,
                    ConstraintType.WRATE => rate - CurrentProducerControl.WaterRateTarget,
                    ConstraintType.GRATE => rate - CurrentProducerControl.GasRateTarget,
                    ConstraintType.LRATE => rate - CurrentProducerControl.LiquidRateTarget,
                    _ => throw new NotSupportedException($"Unsupported producer control mode: {CurrentProducerControl.ControlMode}")
                };
            }

            throw new InvalidOperationException($"Well '{Name}' has no active control record for role {WellType}.");
        }

        internal double[][] GetTrajectoryCoordinates()
        {
            var gridTop = Grids[I + J * Nx];
            var center1 = gridTop.Faces[4].Center;
            double x1 = center1.X, y1 = center1.Y, z1 = center1.Z - 100;
            var gridBottom = Grids[I + J * Nx + Nx * Ny * (Nz - 1)];
            var center2 = gridBottom.Faces[5].Center;
            double x2 = center2.X, y2 = center2.Y, z2 = center2.Z;
            return [[x1, y1, z1], [x2, y2, z2]];
        }

        public readonly struct Offsets
        {
            public readonly int RowFlow_dPwell;
            public readonly int RowFlow_dQwell;
            public readonly int RowConst_dPwell;
            public readonly int RowConst_dQwell;

            // Per perforation:
            // [perfIdx]: grid cell equation sensitivities to Pwell: size nc
            public readonly int[][] Rm_dPwell;
            // [perfIdx]: well rate equation sensitivities to grid cell variables: size nc
            public readonly int[][] RowFlow_dGrid;

            public Offsets(int rflowP, int rflowQ, int rconstP, int rconstQ, int[][] rmDPwell, int[][] rflowDGrid)
            {
                RowFlow_dPwell = rflowP; RowFlow_dQwell = rflowQ;
                RowConst_dPwell = rconstP; RowConst_dQwell = rconstQ;
                Rm_dPwell = rmDPwell; RowFlow_dGrid = rflowDGrid;
            }
        }
    }
}