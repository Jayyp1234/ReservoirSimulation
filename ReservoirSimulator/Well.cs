using EclipseDeckReader;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Well<T> where T : struct, IVectorND<T>
    {
        public string Name { get; set; }
        public Phase PreferedPhase { get; set; }
        public List<int> PerfNatIndex { get; set; } = [];
        public List<double> Perforation_WI { get; set; } = [];
        public ADiff<T> Qo { get; set; }
        public ADiff<T> Qw { get; set; }
        public ADiff<T> Ql { get { return Qo + Qw; } }
        public ADiff<T> Qg { get; set; }
        public ADiff<T> Pwf { get; set; }
        public double Fw { get{  return Qw.Value*100 / (Qo + Qw).Value; } }
        public double ZControl { get; set; }
        public int I { get; }
        public int J { get; }

        public int Nx { get; set; }
        public int Ny { get; set; }
        public int Nz { get; set; }
        public Grid<T>[] Grids { get; set; }
        public int ref_Index { get; set; } 
        public ConstraintType ConstraintType { get; set; } = ConstraintType.LRATE;
        public WellType WellType { get; set; }

        public WellControlRecord CurrentControl
        {
            get { return currentControl; }
            set
            {
                if(currentControl != value)
                {
                    currentControl = value;
                    var Control = currentControl;
                    switch (Control.WellType)
                    {
                        case "PRODUCER":
                            WellType = WellType.PRODUCER;
                            ZControl = double.MaxValue;
                            break;
                        case "INJECTOR":
                            WellType = WellType.INJECTOR;
                            ZControl = double.MinValue;
                            break;
                    }


                    ConstraintType = currentControl.ControlMode switch
                    {
                        "ORAT" => ConstraintType.ORATE,
                        "WRAT" => ConstraintType.WRATE,
                        "GRAT" => ConstraintType.GRATE,
                        "LRAT" => ConstraintType.LRATE,
                        "BHP" => WellType == WellType.INJECTOR ? ConstraintType.MAX_BHP : ConstraintType.MIN_BHP,
                        _ => ConstraintType.ORATE // '_' acts as the default case
                    };

                }
            }
        }
        public WellCompletionRecord CurrentCompletionRecord
        {
            get { return currentcompletionRecord; }
            set
            {
                CompletionChanged = false;
                if (currentcompletionRecord != value)
                {
                    currentcompletionRecord = value;
                    CompletionChanged = true;
                    var comp = currentcompletionRecord;
                    PerfNatIndex = [];
                    double radius = comp.WellBoreDiameter / 2, skin = comp.EffectiveSkin;
                    foreach (int k in Enumerable.Range(comp.KStart, comp.KEnd - comp.KStart + 1))
                    {
                        if (comp.Status != "OPEN") continue;
                        int m = I + J * Nx + k * Nx * Ny;
                        PerfNatIndex.Add(m);
                        var grid = Grids[m];

                        switch (CurrentControl.WellType)
                        {
                            case "PRODUCER":
                                WellType = WellType.PRODUCER;
                                if (grid.Center.Z > ZControl)
                                {
                                    ZControl = grid.Center.Z;
                                    Pwf_n = grid.Po_n;
                                }
                                break;
                            case "INJECTOR":
                                WellType = WellType.INJECTOR;
                                if (grid.Center.Z < ZControl)
                                {
                                    ZControl = grid.Center.Z;
                                    Pwf_n = grid.Po_n;
                                }
                                break;
                        }
                        double DX = 0, DY = 0, DZ = 0;
                        for (int c = 0; c < 8; c++)
                        {
                            // Determine local binary topological offsets
                            double fx = (c == 1 || c == 2 || c == 5 || c == 6) ? 1 : -1;
                            double fy = (c == 2 || c == 3 || c == 6 || c == 7) ? 1 : -1;
                            double fz = (c >= 4) ? 1 : -1;
                            DX += fx * grid.Corners[c].X;
                            DY += fy * grid.Corners[c].Y;
                            DZ += fz * grid.Corners[c].Z;
                        }
                        DX /= 4; DY /= 4; DZ /= 4;
                        double ky = grid.Faces[1].Permeability, kx = grid.Faces[2].Permeability;
                        double re = 0.28 * Hypot(Pow(ky / kx, 0.25) * DX, Pow(kx / ky, 0.25) * DY) /
                                                (Pow(ky / kx, 0.25) + Pow(kx / ky, 0.25));
                        grid.WI = alpha_well * Sqrt(kx * ky) * DZ / (Log(re / radius) + skin);
                    }
                }
            }
        }
        WellCompletionRecord currentcompletionRecord = null;
        WellControlRecord currentControl = null;
        public bool CompletionChanged { get; private set; } = false;
        public bool ActiveStatusChanged { get; private set; } = false;
        public bool isActive { get; private set; } = false;
        public int ActiveIndex { get; set; }
        public double Zref { get; set; }
        public double HighestZref { get; set; }
        public double LowestZref { get; set; }
        public double Pwf_n, Qo_n, Qw_n, Qg_n, Ql_n;

        public List<double> PWF = [], QO = [], QW = [], QG = [];

        public void SaveResult(double? pwf = null, double? qw = null, 
            double? qg = null, double? qo = null, double? sg = null)
        {
            Utility.AddIfPresent(PWF, pwf);
            Utility.AddIfPresent(QG, qg);
            Utility.AddIfPresent(QW, qw);
            Utility.AddIfPresent(QO, qo);
        }
        WellData WellData;
        Phase PreferredPhase;
        double alpha_well = 1.127e-3 * 2 * pi;


        public Well(WellData welldata)
        {
            WellData = welldata;
            Name = WellData.Specification.WellName;
            Zref = WellData.Specification.ReferenceDepth;
            I = WellData.Specification.HeadI;
            J = WellData.Specification.HeadJ;
            ref_Index = WellType == WellType.PRODUCER ? int.MaxValue : 0;
            PreferedPhase = WellData.Specification.PreferredPhase.ToUpperInvariant() switch
            {
                "OIL" => Phase.OIL,
                "WATER" => Phase.WATER,
                "GAS" => Phase.GAS,
                _ => Phase.OIL
            };
        }

        internal List<DateTime> GetEventTimes()
        {
            return [.. WellData.ControlHistory.Select(history => history.EffectiveDate),
                    .. WellData.CompletionHistory.Select(history => history.EffectiveDate)];
        }

        internal void UpdateControl(DateTime time)
        {
            if (WellData.ControlHistory.Count == 0) return;
            if (WellData.ControlHistory[0].EffectiveDate > time) return;
            int m = 0;
            for (m = 0; m < WellData.ControlHistory.Count; m++)
            {
                if (time < WellData.ControlHistory[m].EffectiveDate)
                    break;
            }
            CurrentControl = WellData.ControlHistory[m - 1];
            ActiveStatusChanged = !isActive && CurrentControl.WellStatus == "OPEN" ||
                isActive && CurrentControl.WellStatus == "SHUT";
            isActive = CurrentControl.WellStatus == "OPEN";
        }

        internal void UpdateCompletion(DateTime time)
        {
            if (WellData.CompletionHistory.Count == 0) return;
            if (WellData.CompletionHistory[0].EffectiveDate > time) return;
            int m = 0;
            for (m = 0; m < WellData.CompletionHistory.Count; m++)
            {
                if (time < WellData.CompletionHistory[m].EffectiveDate)
                    break;
            }
            CurrentCompletionRecord = WellData.CompletionHistory[m - 1];
        }

        internal bool IsWaterTarget =>
            ConstraintType == ConstraintType.WRATE ||
            (ConstraintType == ConstraintType.RATE && PreferedPhase == Phase.WATER);

        internal bool IsOilTarget =>
            ConstraintType is ConstraintType.ORATE or ConstraintType.MIN_BHP ||
            (ConstraintType == ConstraintType.RATE && PreferedPhase == Phase.OIL);

        internal bool IsGasTarget =>
            ConstraintType == ConstraintType.GRATE ||
            (ConstraintType == ConstraintType.RATE && PreferedPhase == Phase.GAS);

        internal ADiff<T> Constraint(ADiff<T> pressure, ADiff<T> rate)
        {
            ADiff<T> cons = default;

            switch (ConstraintType)
            {
                case ConstraintType.MIN_BHP:
                    cons = pressure - CurrentControl.MinBhpLimit;
                    break;

                case ConstraintType.MAX_BHP:
                    cons = pressure - CurrentControl.MaxBhpLimit;
                    break;

                case ConstraintType.ORATE:
                    cons = rate - CurrentControl.OilRateTarget;
                    break;

                case ConstraintType.WRATE:
                    cons = rate - CurrentControl.WaterRateTarget;
                    break;

                case ConstraintType.GRATE:
                    cons = rate - CurrentControl.GasRateTarget;
                    break;

                case ConstraintType.LRATE:
                    cons = rate - CurrentControl.LiquidRateTarget;
                    break;

                case ConstraintType.RATE:
                    cons = rate - CurrentControl.TargetRate;
                    break;

                default:
                    throw new NotSupportedException($"Unsupported well control mode: {CurrentControl.ControlMode}");
            }

            return cons;
        }

        internal double[][] GetTrajectoryCoordinates<T>(Grid<T>[] Grids) where T : struct, IVectorND<T>
        {
            var gridTop = Grids[PerfNatIndex.First()]; 
            var center1 = gridTop.Faces.First().Center;  
            double x1 = center1.X, y1 = center1.Y, z1 = center1.Z - 100;
            var gridBottom = Grids[PerfNatIndex.Last()];
            var center2 = gridBottom.Faces.Last().Center;
            double x2 = center2.X, y2 = center2.Y, z2 = center2.Z;
            double[][] coordinates = [[x1, y1, z1], [x2, y2, z2]];
            return coordinates;
        }
    }
}
