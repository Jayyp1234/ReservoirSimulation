using EclipseDeckReader;
using Microsoft.CodeAnalysis;
using SepalSolver;
using System.Text;
using static ReservoirSimulator.Math;
using static ReservoirSimulator.ReservoirJacobian;

namespace ReservoirSimulator
{
    public class OilWaterReservoir
    {
        public int funcall;
        public List<double[]> P, S, Rate, OilRate, WaterRate, GasRate, Pwf, WaterCut;
        public List<double> Time, SweepEff;
        readonly Func<ADiff<Vector4d>, ADiff<Vector4d>> Er;
        (double maxDP, double maxDS, int i_p, int i_s) CheckSolution(double[] xs )
        {
            int i_p = 0, i_s = 0; double dp, maxDp = 0, ds, maxDs = 0;
            for (int i = 0; i < Ngrids; i++)
            {
                dp = Grids[i].Po_np1 - Grids[i].Po_n;
                if(Abs(dp) > maxDp) { maxDp = Abs(dp); i_p = i; }
                ds = Grids[i].Sw_np1 - Grids[i].Sw_n;
                if (Abs(ds) > maxDs) { maxDs = Abs(ds); i_s = i; }
            }
            return (maxDp, maxDs, i_p, i_s);
        }
        double dt, t = 0;
        void InitialGuess(double[] xs, bool isRetry = false, bool wellEventOccurred = false)
        {
            int indx = 0;
            int n = Time.Count - 1;
            bool canExtrapolate = n >= 2 && !isRetry && !wellEventOccurred;
            double omega = canExtrapolate ? Min(1.0, dt / (Time[n] - Time[n - 1])) : 0.0;

            // 1. Reservoir Grid Unknowns (Po, Sw)
            for (int m = 0; m < Ngrids; m++)
            {
                var mgrid = Grids[m];
                if (!mgrid.isActive) continue;

                double poGuess = mgrid.Po_n;
                if (canExtrapolate)
                    poGuess += omega * (mgrid.Po_n - mgrid.PO[n - 1]);

                xs[indx++] = Max(14.7, poGuess); // Po index
                xs[indx++] = mgrid.Sw_n;         // Sw index
            }

            // 2. Active Well Unknowns
            for (int i = 0; i < Nwells; i++)
            {
                var well = Wells[i];
                if (!well.isActive) continue;

                switch (well.ConstraintType)
                {
                    case ConstraintType.MIN_BHP:
                        xs[indx++] = well.CurrentControl.MinBhpLimit;
                        xs[indx++] = well.IsWaterTarget ? well.Qw_n
                                   : well.IsOilTarget ? well.Qo_n
                                   : (well.Qo_n + well.Qw_n);
                        break;

                    case ConstraintType.MAX_BHP:
                        xs[indx++] = well.CurrentControl.MaxBhpLimit;
                        xs[indx++] = well.Qw_n;
                        break;

                    case ConstraintType.RATE:
                        xs[indx++] = well.Pwf_n;
                        xs[indx++] = well.CurrentControl.TargetRate;
                        break;

                    case ConstraintType.GRATE:
                        xs[indx++] = well.Pwf_n;
                        xs[indx++] = well.CurrentControl.GasRateTarget;
                        break;

                    case ConstraintType.ORATE:
                        xs[indx++] = well.Pwf_n;
                        xs[indx++] = well.CurrentControl.OilRateTarget;
                        break;

                    case ConstraintType.WRATE:
                        xs[indx++] = well.Pwf_n;
                        xs[indx++] = well.CurrentControl.WaterRateTarget;
                        break;

                    case ConstraintType.LRATE:
                        xs[indx++] = well.Pwf_n;
                        xs[indx++] = well.CurrentControl.LiquidRateTarget;
                        break;

                    default:
                        throw new InvalidOperationException($"Unhandled well constraint type: {well.ConstraintType}");
                }
            }
        }

        readonly int Nx, Ny, Nz, NxNy, Ngrids, NumActiveGrids, Nwells, NumActiveWells;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z, zTop;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n, Coord, Zcorn, PoreVolume; 
        readonly double krw0, kro0, Pb, Prefw, Prefr, Pe, Pw_woc, Sw_r, So_r,
            Bo0, Bw0, Po_woc, Z_woc, co, cw, cr, bo, bw, nw, no, np,
            μo0, μw0, γo0, γw0, P_datum, Z_datum, Pc_max;
        Grid<Vector4d>[] Grids;
        List<Well<Vector4d>> Wells;
        Aquifer Aquifer;
        double[,] Trans;
        bool[,] Water_Upstream, Oil_Upstream;
        int[] Actnum, Satnum, Pvtnum;
        List<int> Columns2Delete = [];
        Dictionary<(int, int), EQUILPressureDepth> EquilPVTCombinations = [];
        EclipseDataDeck DATA;
        PVTWAnalysis[] PVTWs;
        PVDOAnalysis[] PVDOs;
        SWOFAnalysis[] SWOFs;

        

        public OilWaterReservoir(EclipseDataDeck Data)
        {
            DATA = Data;
            // RUNSPEC
            Nx = Data.NX; 
            Ny = Data.NY; 
            Nz = Data.NZ; 
            NxNy = Nx * Ny; 
            Ngrids = Nx * Ny * Nz;

            Units incomingUnit = Data.UnitSystem.ToUpper() switch
            {
                "FIELD" => Units.FIELD,
                "METRIC" => Units.METRIC,
                "LAB" => Units.LAB,
                "PVT-M" => Units.PVT_M,
                _ => throw new InvalidOperationException($"Unsupported unit system: {Data.UnitSystem}")
            };

            // GRIDS
            if (Data.DxArray.Count > 0 && Data.DyArray.Count > 0 && Data.DzArray.Count > 0)
            {
                Dx = [.. Utility.Convert(Data.DxArray, "Length", incomingUnit, Units.FIELD)]; 
                Dy = [.. Utility.Convert(Data.DyArray, "Length", incomingUnit, Units.FIELD)];
                Dz = [.. Utility.Convert(Data.DzArray, "Depth", incomingUnit, Units.FIELD)];
                zTop = Data.ZtopArray.Count > 0 ? [.. Utility.Convert(Data.ZtopArray, "Depth", incomingUnit, Units.FIELD)] : new double[NxNy];
            }
            else 
            { 
                Coord = [.. Utility.Convert(Data.CoordArray, "Length", incomingUnit, Units.FIELD)]; 
                Zcorn = [.. Utility.Convert(Data.ZcornArray, "Depth", incomingUnit, Units.FIELD)];
            }
            Φ = [.. Data.PorosityArray]; 
            Kx = [.. Data.PermXArray]; Ky = [.. Data.PermYArray]; Kz = [.. Data.PermZArray];
            
            Grids = Zcorn is not null ? Utility.ProcessGrid<Vector4d>(Nx, Ny, Nz, Coord, Zcorn, Φ, Kx, Ky, Kz) :
                Utility.ProcessGrid<Vector4d>(Nx, Ny, Nz, Dx, Dy, Dz, zTop, Φ, Kx, Ky, Kz);

            // PVTW
            PVTWs = [.. Data.RawPvtwTables.Select(pvtw => new PVTWAnalysis(pvtw, incomingUnit))];
            foreach (var pvtw in PVTWs) pvtw.FitProperties();

            // PVDO
            PVDOs = [..Data.RawPvdoTables.Select(pvdo=>new PVDOAnalysis(pvdo, incomingUnit))];
            foreach (var pvdo in PVDOs) pvdo.FitProperties();

            // SWOF
            SWOFs = [.. Data.RawSwofTables.Select(swof => new SWOFAnalysis(swof, incomingUnit))];
            foreach (var swof in SWOFs) swof.FitProperties();

            // ROCK
            cr = Data.RockCompressibility; 
            Prefr = Utility.Convert(Data.RockReferencePressure, "Pressure", incomingUnit, Units.FIELD);
            Er = P => Exp(cr * (P - Prefr));

            // DENSITIES
            γo0 = Utility.Convert(Data.SurfaceOilDensity, "Density", incomingUnit, Units.FIELD)/144;
            γw0 = Utility.Convert(Data.SurfaceWaterDensity, "Density", incomingUnit, Units.FIELD)/144;

            // REGION
            NumActiveGrids = Utility.SetAcitveGrids(Grids, Data.ActnumArray, Er);
            Utility.SetSWOF(Grids, Data.SatnumArray, SWOFs);
            Utility.SetPVDO(Grids, Data.PvtnumArray, PVDOs, γo0); 
            Utility.SetPVTW(Grids, Data.PvtnumArray, PVTWs, γw0);

            // TRANSMISSIBILITY
            Utility.ComputeFaceTransmissibility(Grids);

            // EQUIL
            EQUILAnalysis[] EQUILs = [.. Data.RawEquilTables.Select(equil => new EQUILAnalysis(equil, incomingUnit))];
            double MaxZ = Grids.Max(grid => grid.Center.Z), MinZ = Grids.Min(grid => grid.Center.Z);
            for(int e = 0; e < EQUILs.Length; e++)
            {
                for (int p = 0; p < PVDOs.Length; p++)
                {
                    var equil = EQUILs[e]; var pvdo = PVDOs[p]; var pvtw = PVTWs[p];
                    double BoAtPref = pvdo.Result.BoAtPref, BwAtPref = pvtw.Result.BwAtPref;
                    double Co = pvdo.Result.Co, Cw = pvtw.Result.Cw;
                    double PRefo = pvdo.Result.PRef, PRefw = pvtw.Result.PRef;
                    Func<double, double> γo = P => γo0 / (BoAtPref * Exp(-Co * (P - PRefo))), 
                                         γw = P => γw0 / (BwAtPref * Exp(-Cw * (P - PRefw)));
                    var equildata = new EQUILPressureDepth(equil);
                    int N = (int)(Abs(equildata.ZDatum - equildata.ZWOC) / (MaxZ - MinZ) * Nz * 10);
                    if (equildata.ZDatum < equildata.ZWOC)
                    {
                        var oil2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC, N, γo);
                        equildata.PoWOC = oil2WOC.Last().P; equildata.PwWOC = equildata.PoWOC - equildata.PcWOC;
                    }
                    else
                    {
                        var water2WOC = Utility.IntegrateLegRK4(equil.ZDatum, equil.PDatum, equildata.ZWOC, N, γw);
                        equildata.PwWOC = water2WOC.Last().P; equildata.PoWOC = equildata.PwWOC + equildata.PcWOC;
                    }
                    N = (int)(Abs(MinZ - equildata.ZWOC) / (MaxZ - MinZ) * Nz);
                    var OilPressureHeight1 = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PoWOC, MinZ, N, γo);
                    var WaterPressureHeight1 = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MinZ, N, γw);
                    
                    OilPressureHeight1.Reverse(); equildata.OilZonePressureDepthPairs.AddRange(OilPressureHeight1);
                    WaterPressureHeight1.Reverse(); equildata.WaterZonePressureDepthPairs.AddRange(WaterPressureHeight1);

                    N = (int)((Abs(MaxZ - equildata.ZWOC) / (MaxZ - MinZ) * Nz));
                    var WaterPressureHeight2 = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MaxZ, N, γw);
                    equildata.WaterZonePressureDepthPairs.AddRange(WaterPressureHeight2[1..]);
                    var OilPressureHeight2 = WaterPressureHeight2.Skip(1).Select(pz=>new PressureDepthPair(pz.Z, pz.P + equildata.PcWOC));
                    equildata.OilZonePressureDepthPairs.AddRange(OilPressureHeight2);

                    EquilPVTCombinations[(e, p)] = equildata;
                }
            }

            // WELLS
            Wells = [.. Data.Wells.Select(well => new Well<Vector4d>(well.Value) { Name = well.Key, Grids = Grids })];
        }


        public void Initialize()
        {
            funcall = 0;
            for (int m = 0; m < Ngrids; m++)
            {
                var grid = Grids[m];
                if(!grid.isActive) continue;
                int s = DATA.SatnumArray[m], e = DATA.EquinumArray[m], p = DATA.PvtnumArray[m];

                var swoftable = SWOFs[s].swofTable;
                var equildata = EquilPVTCombinations[(e, p)];
                var oiltable = equildata.OilZonePressureDepthPairs;
                var watertable = equildata.WaterZonePressureDepthPairs;
                grid.Po_n = Utility.Interpolate(oiltable, grid.Center.Z, pair => pair.Z, pair => pair.P);
                grid.Pw_n = Utility.Interpolate(watertable, grid.Center.Z, pair => pair.Z, pair => pair.P);
                double pc = grid.Po_n - grid.Pw_n;
                if (pc == 0)
                    grid.Sw_n = grid.Center.Z < equildata.ZWOC ? swoftable.Min(swof => swof.Sw) : swoftable.Max(swof => swof.Sw);
                else grid.Sw_n = Utility.Interpolate(swoftable, pc, swof => swof.Pc, swof => swof.Sw);
                grid.So_n = 1 - grid.Sw_n;
                ADiff<Vector4d> meanP = grid.So_n * grid.Po_n + grid.Sw_n * grid.Pw_n;
                ADiff<Vector4d> Er_i = grid.Er(meanP);
                ADiff<Vector4d> Bo_i = grid.Bo(grid.Po_n);
                ADiff<Vector4d> Bw_i = grid.Bw(grid.Pw_n);
                grid.ErSo_Bo_n = (Er_i * grid.So_n / Bo_i).Value;
                grid.ErSw_Bw_n = (Er_i * grid.Sw_n / Bw_i).Value;
            }
        }

        public static (double maxPd, double maxSd) ApplyPrimaryVariableUpdate(
            double[] x, double[] dx, int numActiveGrids, int numActiveWells,
            double dPmax = 200.0,   // max allowable pressure change per Newton step [psi]
            double dSmax = 0.20)    // max allowable saturation change per Newton step [-]
        {
            double maxPd = 0.0;
            double maxSd = 0.0;

            // 1. Grid cell updates
            for (int i = 0; i < numActiveGrids; i++)
            {
                int idx_p = 2 * i + 0, idx_s = 2 * i + 1;
                double dPo = dx[idx_p], dSw = dx[idx_s];

                // Damping / step-size limiting on pressure
                if (Abs(dPo) > dPmax) dPo = Sign(dPo) * dPmax;

                // Damping / step-size limiting on saturation
                if (Abs(dSw) > dSmax) dSw = Sign(dSw) * dSmax;

                // Apply update: x = x - dx
                double p_new = x[idx_p] - dPo;
                double s_new = x[idx_s] - dSw;

                // Physical bounds enforcement in Field units [psia] and [0, 1]
                p_new = Max(14.7, p_new);
                s_new = Clamp(s_new, 0.0, 1.0);

                // Track actual applied changes (after damping & clamping)
                double actualPd = Abs(p_new - x[idx_p]);
                double actualSd = Abs(s_new - x[idx_s]);
                if (actualPd > maxPd) maxPd = actualPd;
                if (actualSd > maxSd) maxSd = actualSd;

                x[idx_p] = p_new;
                x[idx_s] = s_new;
            }

            // 2. Well variables update (starts right after active grid DOFs)
            int wellStart = 2 * numActiveGrids;
            for (int w = 0; w < numActiveWells; w++)
            {
                int idxP = wellStart + 2 * w + 0; // Well pressure [psia]
                int idxQ = wellStart + 2 * w + 1; // Well rate [stb/day]

                double dp = Clamp(dx[idxP], -dPmax, dPmax);
                double pNew = Max(14.7, x[idxP] - dp);
                maxPd = Max(maxPd, Abs(pNew - x[idxP]));
                x[idxP] = pNew;

                x[idxQ] -= dx[idxQ]; // Unclamped rate update (allows negative injection rates)
            }

            return (maxPd, maxSd);
        }

        public void Simulate2Phase(double[] ResultTime)
        {
            Initialize();

            int cores = Environment.ProcessorCount;
            int chunkSize = Ngrids / cores;

            int totalDofs = 0;
            int[] col_index = null, row_ptr = null;
            double[] JacValues = null, Res = null, dx = null, xs = null ;

            CellOffsets[] cellOffsets = null;
            FaceOffsets[][][] faceOffsets = null;
            WellOffsets[] wellOffsets = null;
            MklPardisoSolver solver = null;

            funcall = 0;
            double Residual(double[] xnp1, double time)
            {
                funcall++;

                Array.Clear(JacValues, 0, JacValues.Length);
                Array.Clear(Res, 0, Res.Length);

                // =========================================================================
                // 1. HIGH CPU USAGE: Batched Parallel Assembly with Strict Row Ownership
                // =========================================================================
                Parallel.For(0, cores, coreId =>
                {
                    int start = coreId * chunkSize;
                    int end = (coreId == cores - 1) ? Ngrids : start + chunkSize;

                    for (int m = start; m < end; m++)
                    //for (int m = 0; m < Ngrids; m++)
                    {
                        var mgrid = Grids[m];

                        if (!mgrid.isActive) continue;
                        int i = mgrid.ActiveIndex;
                        int row_o = 2 * i + 0, row_w = 2 * i + 1;

                        // Primary variables for cell i
                        var Po_i = new ADiff<Vector4d>(xnp1[2 * i + 0], new(1, 0, 0, 0));
                        var Sw_i = new ADiff<Vector4d>(xnp1[2 * i + 1], new(0, 1, 0, 0));
                        var Pw_i = Po_i - mgrid.Pcow(Sw_i);
                        var So_i = 1.0 - Sw_i; 

                        var meanP = So_i * Po_i + Sw_i * Pw_i;
                        var Er_i = mgrid.Er(meanP);
                        var Bo_i = mgrid.Bo(Po_i); 
                        var Bw_i = mgrid.Bw(Pw_i); 
                        var Kro_i = mgrid.Krow(So_i);
                        var Krw_i = mgrid.Krw(Sw_i);
                        var μo_i = mgrid.μo(Po_i);
                        var μw_i = mgrid.μw(Pw_i);
                        var γo_i = mgrid.γo(Po_i);
                        var γw_i = mgrid.γw(Pw_i);

                        mgrid.Po_np1 = Po_i.Value; mgrid.So_np1 = So_i.Value; mgrid.Bo_np1 = Bo_i.Value;
                        mgrid.Pw_np1 = Pw_i.Value; mgrid.Sw_np1 = Sw_i.Value; mgrid.Bw_np1 = Bw_i.Value;
                        mgrid.ErSo_Bo_np1 = Er_i * So_i / Bo_i; mgrid.ErSw_Bw_np1 = Er_i * Sw_i / Bw_i;

                        // Accumulation Residual: R = Vp * (State^{n+1} - State^n)
                        var ResOil = mgrid.PoreVolume * (mgrid.ErSo_Bo_np1 - mgrid.ErSo_Bo_n);
                        var ResWater = mgrid.PoreVolume * (mgrid.ErSw_Bw_np1 - mgrid.ErSw_Bw_n);

                        // Direct write to Residual Vector
                        Res[row_o] += ResOil.Value;
                        Res[row_w] += ResWater.Value;

                        // Direct write to Accumulation Jacobian (Block i, i)
                        ref readonly var cOff = ref cellOffsets[i];
                        JacValues[cOff.Diag[0]] += ResOil.Derivative[0];   // row_o, col_p_i
                        JacValues[cOff.Diag[1]] += ResOil.Derivative[1];   // row_o, col_sw_i
                        JacValues[cOff.Diag[2]] += ResWater.Derivative[0]; // row_w, col_p_i
                        JacValues[cOff.Diag[3]] += ResWater.Derivative[1]; // row_w, col_sw_i


                        // Inter-cell flux loop
                        int faceCount = mgrid.Faces.Length;
                        for (int f = 0; f < faceCount; f++)
                        {
                            var face = mgrid.Faces[f];
                            for (int c = 0; c < face.GridIndices.Count; c++)
                            {
                                int n = face.GridIndices[c];
                                var ngrid = Grids[n];

                                // Guard: Skip inactive neighbor blocks (no-flow boundary)
                                if (!ngrid.isActive) continue;

                                int j = ngrid.ActiveIndex;
                                double transDt = face.Transmissibility[c] * dt;
                                double dz = mgrid.Center.Z - ngrid.Center.Z;

                                var Po_j = new ADiff<Vector4d>(xnp1[2 * j + 0], new(0, 0, 1, 0));
                                var Sw_j = new ADiff<Vector4d>(xnp1[2 * j + 1], new(0, 0, 0, 1));
                                var Pw_j = Po_j - ngrid.Pcow(Sw_j);
                                var So_j = 1.0 - Sw_j;
                                var Bo_j = ngrid.Bo(Po_j);
                                var Bw_j = ngrid.Bw(Pw_j);
                                var Kro_j = ngrid.Krow(So_j);
                                var Krw_j = ngrid.Krw(Sw_j);
                                var μo_j = ngrid.μo(Po_j);
                                var μw_j = ngrid.μw(Pw_j);
                                var γo_j = ngrid.γo(Po_j);
                                var γw_j = ngrid.γw(Pw_j);

                                // Hydrostatic potential differences
                                var γo_avg = 0.5 * (γo_i + γo_j);
                                var γw_avg = 0.5 * (γw_i + γw_j);
                                var DΦo = Po_i - Po_j - γo_avg * dz;
                                var DΦw = Pw_i - Pw_j - γw_avg * dz;

                                // Upstream phase mobilities
                                var lambda_o = DΦo >= 0.0 ? (Kro_i / (μo_i * Bo_i)) : (Kro_j / (μo_j * Bo_j));
                                var lambda_w = DΦw >= 0.0 ? (Krw_i / (μw_i * Bw_i)) : (Krw_j / (μw_j * Bw_j));

                                // Outflow flux from cell i to cell j
                                var FluxOil = transDt * lambda_o * DΦo;
                                var FluxWater = transDt * lambda_w * DΦw;

                                // Direct write to Flux Residual
                                Res[row_o] += FluxOil.Value;
                                Res[row_w] += FluxWater.Value;

                                // Direct write to Jacobian using precomputed face offsets
                                ref readonly var fOff = ref faceOffsets[m][f][c];

                                // Block (i, i) -> Derivatives wrt Cell i variables (slots 0 and 1)
                                JacValues[fOff.Self[0]] += FluxOil.Derivative[0];   // row_o, col_p_i
                                JacValues[fOff.Self[1]] += FluxOil.Derivative[1];   // row_o, col_sw_i
                                JacValues[fOff.Self[2]] += FluxWater.Derivative[0]; // row_w, col_p_i
                                JacValues[fOff.Self[3]] += FluxWater.Derivative[1]; // row_w, col_sw_i

                                // Block (i, j) -> Derivatives wrt Neighbor j variables (slots 2 and 3)
                                JacValues[fOff.Neighbor[0]] += FluxOil.Derivative[2];   // row_o, col_p_j
                                JacValues[fOff.Neighbor[1]] += FluxOil.Derivative[3];   // row_o, col_sw_j
                                JacValues[fOff.Neighbor[2]] += FluxWater.Derivative[2]; // row_w, col_p_j
                                JacValues[fOff.Neighbor[3]] += FluxWater.Derivative[3]; // row_w, col_sw_j
                            }
                        }
                    }
                });

                // =========================================================================
                // 2. Sequential Well System Assembly & Direct Memory Updates
                // =========================================================================
                foreach (var well in Wells)
                {
                    if (!well.isActive) continue;
                    ref readonly var wOff = ref wellOffsets[well.ActiveIndex];

                    int row_p = 2 * NumActiveGrids + 2 * well.ActiveIndex + 0; // Constraint equation
                    int row_q = 2 * NumActiveGrids + 2 * well.ActiveIndex + 1; // Rate definition equation

                    var pwell = new ADiff<Vector4d>(xnp1[row_p], new(1, 0, 0, 0));
                    var qwell = new ADiff<Vector4d>(xnp1[row_q], new(0, 1, 0, 0));

                    double Zref = well.ZControl;
                    well.Qo = default;
                    well.Qw = default;

                    int perfCount = well.PerfNatIndex.Count;
                    bool isWaterTarget = well.IsWaterTarget;
                    bool isOilTarget = well.IsOilTarget;

                    switch (well.WellType)
                    {
                        case WellType.PRODUCER:
                            {
                                for (int p = 0; p < perfCount; p++)
                                {
                                    int m = well.PerfNatIndex[p];
                                    var mgrid = Grids[m];
                                    if (!mgrid.isActive) continue;

                                    int i = mgrid.ActiveIndex;
                                    int indx1 = 2 * i + 0;
                                    int indx2 = 2 * i + 1;

                                    var Po_m = new ADiff<Vector4d>(xnp1[indx1], new(0, 0, 1, 0));
                                    var Sw_m = new ADiff<Vector4d>(xnp1[indx2], new(0, 0, 0, 1));
                                    var Pw_m = Po_m - mgrid.Pcow(Sw_m);
                                    var So_m = 1.0 - Sw_m;

                                    var Bo_m = mgrid.Bo(Po_m); 
                                    var Bw_m = mgrid.Bw(Pw_m);
                                    var Kro_m = mgrid.Krow(So_m);
                                    var Krw_m = mgrid.Krw(Sw_m);
                                    var μo_m = mgrid.μo(Po_m);
                                    var μw_m = mgrid.μw(Pw_m);
                                    var γo_m = mgrid.γo(Po_m);
                                    var γw_m = mgrid.γw(Pw_m);

                                    var γo_well = mgrid.γo(pwell);
                                    var γw_well = mgrid.γw(pwell);

                                    var go_avg = 0.5 * (γo_m + γo_well);
                                    var gw_avg = 0.5 * (γw_m + γw_well);
                                    double dz = mgrid.Center.Z - Zref;

                                    var DΦo = Po_m - pwell - go_avg * dz;
                                    var DΦw = Pw_m - pwell - gw_avg * dz;

                                    var WIo = mgrid.WI * (Kro_m / (μo_m * Bo_m));
                                    var WIw = mgrid.WI * (Krw_m / (μw_m * Bw_m));

                                    var rate_oil = WIo * DΦo;
                                    var rate_water = WIw * DΦw;

                                    well.Qo += rate_oil;
                                    well.Qw += rate_water;

                                    // Reservoir residuals & Jacobian entries (outflow is +rate * dt)
                                    Res[indx1] += rate_oil.Value * dt;
                                    Res[indx2] += rate_water.Value * dt;

                                    JacValues[wOff.Rm_dPwell[p][0]] += rate_oil.Derivative[0] * dt;
                                    JacValues[wOff.Rm_dPwell[p][1]] += rate_water.Derivative[0] * dt;

                                    ref readonly var cOff = ref cellOffsets[i];
                                    JacValues[cOff.Diag[0]] += rate_oil.Derivative[2] * dt;
                                    JacValues[cOff.Diag[1]] += rate_oil.Derivative[3] * dt;
                                    JacValues[cOff.Diag[2]] += rate_water.Derivative[2] * dt;
                                    JacValues[cOff.Diag[3]] += rate_water.Derivative[3] * dt;

                                    // Per-perforation grid derivatives for Row Q
                                    var perfRate = isWaterTarget ? rate_water
                                                    : isOilTarget ? rate_oil
                                                    : rate_oil + rate_water;

                                    JacValues[wOff.RowQ_dGrid[p][0]] += -perfRate.Derivative[2];
                                    JacValues[wOff.RowQ_dGrid[p][1]] += -perfRate.Derivative[3];
                                }

                                // Close Row Q using matching accumulated well rate
                                ADiff<Vector4d> prodRate = isWaterTarget ? well.Qw
                                                         : isOilTarget ? well.Qo
                                                         : well.Ql;

                                Res[row_q] += qwell.Value - prodRate.Value;
                                JacValues[wOff.RowQ_Pwell] += -prodRate.Derivative[0];
                                JacValues[wOff.RowQ_Qwell] += 1.0;
                                break;
                            }

                        case WellType.INJECTOR:
                            {
                                for (int p = 0; p < perfCount; p++)
                                {
                                    int m = well.PerfNatIndex[p];
                                    var mgrid = Grids[m];
                                    if (!mgrid.isActive) continue;

                                    int i = mgrid.ActiveIndex;
                                    int indx1 = 2 * i + 0;
                                    int indx2 = 2 * i + 1;

                                    var Po_m = new ADiff<Vector4d>(xnp1[indx1], new(0, 0, 1, 0));
                                    var Sw_m = new ADiff<Vector4d>(xnp1[indx2], new(0, 0, 0, 1));
                                    var Pw_m = Po_m - mgrid.Pcow(Sw_m);

                                    var Bw_m = mgrid.Bw(Pw_m);
                                    var μw_m = mgrid.μw(Pw_m);

                                    var γw_well = mgrid.γw(pwell);
                                    var γw_m = mgrid.γw(Pw_m);

                                    var gw_avg = 0.5 * (γw_m + γw_well);
                                    double dz = mgrid.Center.Z - Zref;

                                    var DΦw = pwell - Pw_m + gw_avg * dz;
                                    // Phase-weighted total mobility in perforated cell
                                    var lambdaTotal =  1.0 / (μw_m * Bw_m);
                                    var WIw = mgrid.WI * lambdaTotal;
                                    var rate_water = WIw * DΦw;

                                    well.Qw += rate_water;

                                    // Reservoir water equation (injection is inflow: -rate * dt)
                                    Res[indx2] += -rate_water.Value * dt;

                                    JacValues[wOff.Rm_dPwell[p][1]] += -rate_water.Derivative[0] * dt;

                                    ref readonly var cOff = ref cellOffsets[i];
                                    JacValues[cOff.Diag[2]] += -rate_water.Derivative[2] * dt;
                                    JacValues[cOff.Diag[3]] += -rate_water.Derivative[3] * dt;

                                    JacValues[wOff.RowQ_dGrid[p][0]] += -rate_water.Derivative[2];
                                    JacValues[wOff.RowQ_dGrid[p][1]] += -rate_water.Derivative[3];
                                }

                                // Close Row Q using well.Qw
                                Res[row_q] += qwell.Value - well.Qw.Value;
                                JacValues[wOff.RowQ_Pwell] += -well.Qw.Derivative[0];
                                JacValues[wOff.RowQ_Qwell] += 1.0;
                                break;
                            }
                    }

                    // Close the Well Constraint Equation (Row p)
                    ADiff<Vector4d> wellConstraint = well.Constraint(pwell, qwell);
                    Res[row_p] += wellConstraint.Value;
                    JacValues[wOff.RowP_Pwell] += wellConstraint.Derivative[0];
                    JacValues[wOff.RowP_Qwell] += wellConstraint.Derivative[1];
                }

                // =========================================================================
                // 3. Compute Normalized Residual
                // =========================================================================

                // At the very end of Residual(xs, tnp1), after accumulation + faces + wells are all in Res:

                double MaxCNV = 0;
                foreach (var grid in Grids)
                {
                    if (!grid.isActive) continue;

                    int i = grid.ActiveIndex;
                    double cnvOil = Abs(Res[2 * i + 0]) * grid.Bo_np1 * grid.InvPoreVolume;
                    double cnvWater = Abs(Res[2 * i + 1]) * grid.Bw_np1 * grid.InvPoreVolume;

                    if (cnvOil > MaxCNV) MaxCNV = cnvOil;
                    if (cnvWater > MaxCNV) MaxCNV = cnvWater;
                }

                // If you also have explicit well BHP equations at rows [2*NumActiveGrids .. totalRows-1]:
                int wellOffset = 2 * NumActiveGrids;
                foreach (var well in Wells)
                {
                    if (!well.isActive) continue;

                    int row0 = wellOffset + 2 * well.ActiveIndex + 0;
                    int row1 = wellOffset + 2 * well.ActiveIndex + 1;

                    // Equation 0: Wellbore rate balance / inflow equation [stb/day or rb/day]
                    double rateScale = Max(1.0, Abs(well.CurrentControl.TargetRate));
                    double err0 = Abs(Res[row0]) / rateScale;

                    // Equation 1: Active control constraint (Rate [stb/day] or BHP [psi])
                    double ctrlScale = well.ConstraintType is ConstraintType.MAX_BHP or ConstraintType.MIN_BHP 
                        ? Max(1.0, Abs(well.CurrentControl.MaxBhpLimit)) : Max(1.0, Abs(well.CurrentControl.TargetRate));
                    double err1 = Abs(Res[row1]) / ctrlScale;

                    if (err0 > MaxCNV) MaxCNV = err0;
                    if (err1 > MaxCNV) MaxCNV = err1;
                }

                return MaxCNV;
            }

            dt = 0.1;
            // Initialize historical data tracking containers for plotting and reporting
            DateTime startDate = DATA.StartDate, currentSimDate;
            List<DateTime> Interval = [startDate, .. ResultTime.Select(time => startDate.AddDays(time))];
            foreach (var well in Wells)
                Interval.AddRange(well.GetEventTimes());
            Interval = [.. Interval.Distinct().OrderBy(x => x)];
            Console.WriteLine($"""
                    Total Grids = {Ngrids}, 
                    Active Grids = {NumActiveGrids}
                    ======================================================================
                                            Starting simulation
                    """);
            double L = 0, tnp1, history_dtmax = 0, rnorm = 0;
            bool staterejected, isComplete;
            Time = [0];
            for (int i = 1; i < Interval.Count; i++)
            {
                L = Interval[i].Subtract(startDate).TotalDays;
                dt = Min(dt, Abs(L - t));
                isComplete = false;
                while (!isComplete)
                {
                    staterejected = false;
                    tnp1 = t + dt;
                    currentSimDate = startDate.AddDays(tnp1);
                    Console.WriteLine($"\n[Info] Simulation Time: {tnp1:F3} days | Current Date: {currentSimDate:yyyy-MM-dd}");

                    // Check for any Well events (e.g., control changes, shut-ins) that occur at the current simulation time
                    // If any events are detected, update the well constraints accordingly
                    int NumActiveWells = 0;
                    foreach (var well in Wells)
                    {
                        well.UpdateControl(currentSimDate);
                        if(well.isActive) well.ActiveIndex = NumActiveWells++;
                    }

                    // Check for changes in Well Completion Status (e.g., new completions, shut-ins) that may affect the simulation
                    // if there are any changes, update the well completions and related properties
                    foreach (var well in Wells) well.UpdateCompletion(currentSimDate);

                    // Check for any changes in Well Constraints (e.g., BHP limits, rate limits) that may affect the simulation
                    if (Wells.Any(well => well.ActiveStatusChanged) || Wells.Any(well=> well.CompletionChanged))
                    {
                        Console.WriteLine($"[Info] Number of Active Wells {NumActiveWells}");
                        Console.WriteLine($"[Info] Well Completion Status Changed at t = {tnp1:F3} days for some Wells");
                        Console.WriteLine($"New Jacobian Nonzero index and Symbolic factorization");
                        (totalDofs, col_index, row_ptr, JacValues) = ReservoirJacobian.Create(2, Grids, Wells);
                        (cellOffsets, faceOffsets, wellOffsets) = ReservoirJacobian.PrecomputeOffsets(2, row_ptr, col_index, Grids, Wells);
                        solver = new MklPardisoSolver(totalDofs, row_ptr, col_index);
                        xs = new double[totalDofs]; Res = new double[totalDofs]; dx = new double[totalDofs];
                    }

                    // Prepare the initial guess for the Newton-Raphson solver based on the previous time step solution
                    
                    InitialGuess(xs);

                    // Call the Newton-Raphson nonlinear solver to
                    // find the next state solution
                    Console.WriteLine($"""



                    Time: 
                    {tnp1:E4} days
                        iter  |   Residual Norm  
                    ----------+------------------
                    """);
                    int iter; bool isConverged = false;
                    rnorm = Residual(xs, tnp1);
                    // Solve the nonlinear system using Newton-Raphson method
                    for (iter = 1; iter < 7; iter++)
                    {
                        // 2. Direct solve: J * dx = Res
                        solver.Solve(JacValues, Res, dx);

                        // 3. Update primary variables in field units: x^{k+1} = x^k - dx
                        (double maxDPo, double maxDSw) = ApplyPrimaryVariableUpdate(xs, dx, NumActiveGrids, NumActiveWells);
                        rnorm = Residual(xs, tnp1);

                        Console.WriteLine($"  {iter,4}    |    {rnorm:E4}");
                        // 4. Converged only when BOTH normalized residual AND state step are small
                        if (rnorm < 1e-6 && maxDSw < 1e-4 && maxDPo < 1e-2)
                        {
                            isConverged = true;
                            break;
                        }
                    }

                    // Check convergence. If non-converged,
                    // chop the time step (time-step cuts) and retry.
                    if (!isConverged)
                    {
                        dt = 0.1*dt;
                        Console.WriteLine("""
                                ============================
                                 Rejected (Non-Convergence)
                                ============================
                                """);
                        continue;
                    }

                    // Unpack solution values to evaluate operational constraint validations
                    var (maxDP, maxDS, i_p, i_s) = CheckSolution(xs);
                    if (maxDP > 400 || maxDS > 0.15)
                    {
                        int j = maxDP > 400 ? i_p : i_s;
                        double limit = maxDP > 400 ? 400 : 0.15;
                        double change = maxDP > 400 ? maxDP : maxDS;
                        Console.ForegroundColor = ConsoleColor.Yellow;
                        Console.WriteLine($"\n[Warning] Time Step Rejected: State Change Limit Exceeded.");
                        Console.ResetColor();
                        Console.WriteLine($"    ↳ Location       : Cell Index {j}");
                        Console.WriteLine($"    ↳ Violation      : {(maxDP > 400 ? "ΔP" : "ΔS")} = {change:F4} (Max Allowed: {limit:F4})");
                        Console.WriteLine($"    ↳ Current State  : P = {Grids[j].Po_n:F1} psi | Sw = {Grids[j].Sw_n:F4}");
                        Console.WriteLine($"    ↳ Attempted State: P = {Grids[j].Po_np1:F1} psi | Sw = {Grids[j].Sw_np1:F4}");
                        Console.WriteLine($"    ↳ Solver Action  : Rolling back to t = {t:F3} days.");
                        Console.WriteLine($"    ↳ Time-Step Chop : Reducing Δt from {dt:F2} to {0.5*dt:F2} days (Factor: 0.5).\n");
                        dt = 0.5*dt;
                        continue;
                    }

                    // Track the maximum time step that still yields convergence for potential future use in adaptive time-stepping
                    history_dtmax = Max(dt, history_dtmax);
                    foreach(var well in Wells)
                    {
                        if (!well.isActive) continue;
                        // Validate Producer constraints:
                        if (well.WellType == WellType.PRODUCER)
                        {
                            // switch to BHP control if pressure falls below minimum limits
                            if (well.ConstraintType != ConstraintType.MIN_BHP && well.Pwf < well.CurrentControl.MinBhpLimit)
                            {
                                well.ConstraintType = ConstraintType.MIN_BHP;
                                Console.WriteLine($"""
                                ===============================================================
                                      Rejected (Minimum Pressure Violated) @ {well.Name}
                                ===============================================================
                                """);
                                staterejected = true;
                            }
                        }

                        // Validate Injector constraints:
                        if (well.WellType == WellType.INJECTOR)
                        {
                            // switch to BHP control if pressure exceeds fracturing limits
                            if (well.ConstraintType != ConstraintType.MAX_BHP && well.Pwf > well.CurrentControl.MaxBhpLimit)
                            {
                                well.ConstraintType = ConstraintType.MAX_BHP;
                                Console.WriteLine("""
                                ================================================
                                      Rejected (Maximum Pressure Violated) 
                                ================================================
                                """);
                                staterejected = true;
                            }
                        }
                    }
                    if (staterejected) continue;
                    

                    // Log verified parameters to performance history arrays
                    Time.Add(t = tnp1); 
                    foreach(var grid in Grids)
                    { 
                        if (!grid.isActive) continue;
                        grid.ErSo_Bo_n = grid.ErSo_Bo_np1.Value;
                        grid.ErSw_Bw_n = grid.ErSw_Bw_np1.Value;
                        grid.Po_n = grid.Po_np1; grid.Sw_n = grid.Sw_np1;
                        grid.Pw_n = grid.Pw_np1; grid.So_n = grid.So_np1;
                        grid.SaveResult(po: grid.Po_n, pw: grid.Pw_n, so: grid.So_n, sw: grid.Sw_n);
                    }
                    foreach(var well in Wells)
                    {
                        if (!well.isActive) continue;
                        well.Pwf_n = well.Pwf.Value;
                        well.Qo_n = well.Qo.Value;
                        well.Qw_n = well.Qw.Value;
                        well.SaveResult(pwf: well.Pwf_n, qo: well.Qo_n, qw: well.Qw_n);
                    }


                    isComplete = Abs(t - L) < 1e-13;
                    // Adaptive Time-Stepping Logic:
                    // scale dt up if convergence is fast, scale down if slow
                    if (iter < 3) dt = 1.25*dt;
                    if (iter > 5) dt = 0.5*dt;
                    if (dt < 1e-12) throw new Exception("time step is too small");

                    // Prevent Overshoot
                    if (!isComplete) dt = Min(dt, L - t);
                    dt = Min(dt, 100);
                }
                if (history_dtmax > 0) dt = 0.1*history_dtmax;
            }
        }

        private double[] xCoords, yCoords, zCoords;
        public void ExportParaViewVTR(string folderPath)
        {
            xCoords = CalculateCoordinates([.. Enumerable.Range(0, Nx).Select(i => Dx[i])]);
            yCoords = CalculateCoordinates([.. Enumerable.Range(0, Ny).Select(i => Dy[i*Nx])]);
            zCoords = CalculateCoordinates([.. Enumerable.Range(0, Nz).Select(i => Dz[i*Nx*Ny])]);
            ExportSimulation(folderPath, "Test1", Nx, Ny, Nz,
                xCoords, yCoords, zCoords, Time, P, S);
        }

        /// <summary>
        /// Exports a single-block simulation time-series tracking Saturation at node corners.
        /// Perfectly mirrors ParaView's native parallel XML layout (PVD + PVTR + VTR).
        /// </summary>
        public static void ExportSimulation(
                string outputDirectory,
                string caseName,
                int Nx, int Ny, int Nz,
                double[] xCoord,        // Length: Nx + 1 (Nodes)
                double[] yCoord,        // Length: Ny + 1 (Nodes)
                double[] zCoord,        // Length: Nz + 1 (Nodes)
                List<double> Time,      // Length: Total timesteps
                List<double[]> P,       // List of flat cell-centered arrays (length: Nx * Ny * Nz)
                List<double[]> S)       // List of flat cell-centered arrays (length: Nx * Ny * Nz)
        {
            Directory.CreateDirectory(outputDirectory);
            DirectoryInfo directory = new DirectoryInfo(outputDirectory);
            foreach (FileInfo file in directory.EnumerateFiles())
                file.Delete();

            string fFormat = "G17"; // High-precision float string format
            var pvdSb = new StringBuilder();

            // 1. GENERATE THE TIMELINE COLLECTION MASTER (.pvd)
            pvdSb.AppendLine("<?xml version=\"1.0\"?>");
            pvdSb.AppendLine("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">");
            pvdSb.AppendLine("  <Collection>");

            int Nt = (int)(Time.Last()/20);
            double delt = Time.Last()/Nt;
            for (int framenum = 0; framenum <= 1000; framenum++)
            {
                double t = framenum*delt;
                double currentTime = t;
                double[] cellS = [];
                double[] cellP = [];

                // Unique filename for each timestep's grid data in the same directory
                string vtrFileName = $"{caseName}Time{framenum}.vtr";

                // PVD now links directly to the standalone VTR file
                pvdSb.AppendLine($"    <DataSet timestep=\"{currentTime.ToString(fFormat)}\" group=\"\" part=\"0\" file=\"{vtrFileName}\"/>");

                // -----------------------------------------------------------------
                // 2. GENERATE THE DYNAMIC DATA SNAPSHOT (.vtr)
                // -----------------------------------------------------------------
                var vtrSb = new StringBuilder();
                vtrSb.AppendLine("<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">");
                vtrSb.AppendLine($"  <RectilinearGrid WholeExtent=\"0 {Nx} 0 {Ny} 0 {Nz}\">");

                // Optional time value metadata
                vtrSb.AppendLine("    <FieldData>");
                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">");
                vtrSb.AppendLine($"        {currentTime.ToString(fFormat)}");
                vtrSb.AppendLine("      </DataArray>");
                vtrSb.AppendLine("    </FieldData>");

                vtrSb.AppendLine($"    <Piece Extent=\"0 {Nx} 0 {Ny} 0 {Nz}\">");

                // CellData Block maps simulation variables directly to the centers of your grid blocks
                vtrSb.AppendLine("    <CellData Scalars=\"Saturation,Pressure\">");

                // Saturation Array
                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Saturation\" format=\"ascii\" RangeMin=\"{cellS.Min()}\" RangeMax=\"{cellS.Max()}\">");
                vtrSb.Append("        ");
                int cellIndx = 0;
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            vtrSb.Append(cellS[cellIndx++].ToString(fFormat) + " ");
                        }
                        vtrSb.Append("\n        ");
                    }
                    vtrSb.Append("\n\n        ");
                }
                vtrSb.AppendLine("\n      </DataArray>\n");

                // Pressure Array
                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\" RangeMin=\"{cellP.Min()}\" RangeMax=\"{cellP.Max()}\">");
                vtrSb.Append("        ");
                cellIndx = 0;
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            vtrSb.Append(cellP[cellIndx++].ToString(fFormat) + " ");
                        }
                        vtrSb.Append("\n        ");
                    }
                    vtrSb.Append("\n\n        ");
                }
                vtrSb.AppendLine("\n      </DataArray>");
                vtrSb.AppendLine("    </CellData>");

                // Coordinates Block defines the physical bounding nodes of your grid blocks
                vtrSb.AppendLine("    <Coordinates>");

                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"X\" format=\"ascii\" RangeMin=\"{xCoord.Min()}\" RangeMax=\"{xCoord.Max()}\">");
                vtrSb.Append("        ");
                for (int i = 0; i <= Nx; i++) vtrSb.Append(xCoord[i].ToString(fFormat) + " ");
                vtrSb.AppendLine("\n      </DataArray>");

                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Y\" format=\"ascii\" RangeMin=\"{yCoord.Min()}\" RangeMax=\"{yCoord.Max()}\">");
                vtrSb.Append("        ");
                for (int j = 0; j <= Ny; j++) vtrSb.Append(yCoord[j].ToString(fFormat) + " ");
                vtrSb.AppendLine("\n      </DataArray>");

                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Z\" format=\"ascii\" RangeMin=\"{zCoord.Min()}\" RangeMax=\"{zCoord.Max()}\">");
                vtrSb.Append("        ");
                for (int k = 0; k <= Nz; k++) vtrSb.Append(zCoord[k].ToString(fFormat) + " ");
                vtrSb.AppendLine("\n      </DataArray>");

                vtrSb.AppendLine("    </Coordinates>");

                vtrSb.AppendLine("    </Piece>");
                vtrSb.AppendLine("  </RectilinearGrid>");
                vtrSb.AppendLine("</VTKFile>");

                File.WriteAllText(Path.Combine(outputDirectory, vtrFileName), vtrSb.ToString());
            }

            pvdSb.AppendLine("  </Collection>");
            pvdSb.AppendLine("</VTKFile>");
            File.WriteAllText(Path.Combine(outputDirectory, $"{caseName}.pvd"), pvdSb.ToString());
        }

        public void ExportWells(string outputDirectory)
        {
            Directory.CreateDirectory(outputDirectory);
            foreach (var well in Wells)
            {
                double[][] Trajectory = well.GetTrajectoryCoordinates(Grids);
                ExportWellAsVtk(Path.Combine(outputDirectory, $"{well.Name}.vtk"), well.Name, Trajectory);
            }
        }

        /// <summary>
        /// Computes cumulative edge locations from cell delta sizes.
        /// Length of returning array is cell count + 1.
        /// </summary>
        private double[] CalculateCoordinates(double[] deltas)
        {
            double[] coords = new double[deltas.Length + 1];
            coords[0] = 0.0;
            for (int i = 0; i < deltas.Length; i++)
            {
                coords[i + 1] = coords[i] + deltas[i];
            }
            return coords;
        }

        public void ExportWellAsVtk(string outputPath, string wellName, double[][] trajectories)
        {
            // trajectories array expected shape: [numPoints][3] -> X, Y, Z coordinates
            int numPoints = trajectories.Length;

            using StreamWriter writer = new(outputPath, false, Encoding.ASCII);

            // 1. Write standard VTK Polydata Header
            writer.WriteLine("# vtk DataFile Version 3.0");
            writer.WriteLine($"Well Path: {wellName}");
            writer.WriteLine("ASCII");
            writer.WriteLine("DATASET POLYDATA");

            // 2. Write Geometry Nodes (Points)
            writer.WriteLine($"POINTS {numPoints} double");
            for (int i = 0; i < numPoints; i++)
            {
                writer.WriteLine($"{trajectories[i][0]:F2} {trajectories[i][1]:F2} {trajectories[i][2]:F2}");
            }
            writer.WriteLine();

            // 3. Topology: Link the points into a single continuous wireframe string
            // Format: LINES [NumberOfLines] [TotalIntegerCountInBlock]
            // Total count = NumberOfLines + total number of indices mapped
            writer.WriteLine($"LINES 1 {numPoints + 1}");
            writer.Write($"{numPoints} ");
            for (int i = 0; i < numPoints; i++)
            {
                writer.Write($"{i} ");
            }
            writer.WriteLine();
            writer.WriteLine();

            // 4. Optional Metadata Attributions (Coloring / Labeling in ParaView)
            writer.WriteLine($"CELL_DATA 1");
            writer.WriteLine("FIELD FieldData 1");
            writer.WriteLine($"WellID 1 {wellName.Length} string");
            writer.WriteLine(wellName);
        }


        // 1. Generates the individual active mesh snapshot for a specific time step
        public void ExportStepToVtu(
            string directoryPath, string baseName, int stepIndex,
            int nx, int ny, int nz, double[] dx, double[] dy, double[] dz, double[] zTop,
            int[] Actnum, double[] pressure, double[] saturation)
        {
            // Build the specific filename for this step (e.g., RunTest1_0.vtu)
            string fileName = $"{baseName}_{stepIndex}.vtu";
            string fullPath = Path.Combine(directoryPath, fileName);

            // 2. Compute explicit structural coordinates for vertices
            // A regular mesh has (nx+1)*(ny+1)*(nz+1) absolute grid vertex coordinates
            double[] xCoords = new double[nx + 1], yCoords = new double[ny + 1];

            for (int i = 1; i <= nx; i++) xCoords[i] = xCoords[i - 1] + dx[i - 1];
            for (int j = 1; j <= ny; j++) yCoords[j] = yCoords[j - 1] + dy[j - 1];

            StringBuilder xmlPoints = new(), xmlConnectivity = new(), xmlOffsets = new();
            StringBuilder xmlTypes = new(), xmlPressure = new(), xmlSaturation = new();

            int pointOffset = 0;
            int Ngrids = pressure.Length;
            int numActiveCells = Actnum.Where(a => a > 0).Count();
            int numPoints = numActiveCells*8; 
            for (int m = 0; m < Ngrids; m++)
            {
                if (Actnum[m] == 0) continue;
                var grid = Grids[m];

                for (int c = 0; c < 8; c++)
                    xmlPoints.AppendLine($"{grid.Corners[c].X} {grid.Corners[c].Y} {grid.Corners[c].Z}");     

                // Hexahedron connectivity mapping string index pointers
                xmlConnectivity.AppendLine($"{pointOffset} {pointOffset+1} {pointOffset+2} {pointOffset+3} {pointOffset+4} {pointOffset+5} {pointOffset+6} {pointOffset+7}");

                pointOffset += 8;
                xmlOffsets.AppendLine($"{pointOffset}");
                xmlTypes.AppendLine("12"); // VTK_HEXAHEDRON structural code identifier

                // Map data properties directly using the corresponding active global identifier
                xmlPressure.Append($"{pressure[m]} ");
                xmlSaturation.Append($"{saturation[m]} ");
            }

            // 3. Assemble complete XML output file buffer
            StringBuilder vtuFile = new();
            vtuFile.AppendLine("<?xml version=\"1.0\"?>");
            vtuFile.AppendLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
            vtuFile.AppendLine($"  <UnstructuredGrid>");
            vtuFile.AppendLine($"    <Piece NumberOfPoints=\"{numPoints}\" NumberOfCells=\"{numActiveCells}\">");

            // Point Data Section
            vtuFile.AppendLine("      <Points>");
            vtuFile.AppendLine("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">");
            vtuFile.Append(xmlPoints.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("      </Points>");

            // Cell Connectivity Specification Section
            vtuFile.AppendLine("      <Cells>");
            vtuFile.AppendLine("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
            vtuFile.Append(xmlConnectivity.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
            vtuFile.Append(xmlOffsets.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
            vtuFile.Append(xmlTypes.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("      </Cells>");

            // Simulation Physical Property Field Section
            vtuFile.AppendLine("      <CellData Scalars=\"Pressure\">");
            vtuFile.AppendLine("        <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">");
            vtuFile.AppendLine(xmlPressure.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("        <DataArray type=\"Float64\" Name=\"Water_Saturation\" format=\"ascii\">");
            vtuFile.AppendLine(xmlSaturation.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("      </CellData>");

            vtuFile.AppendLine("    </Piece>");
            vtuFile.AppendLine("  </UnstructuredGrid>");
            vtuFile.AppendLine("</VTKFile>");

            File.WriteAllText(fullPath, vtuFile.ToString());
        }

        // 2. Generates the master pointer file that ParaView opens
        public void ExportParaViewVTU(string outputDirectory)
        {
            Directory.CreateDirectory(outputDirectory);
            DirectoryInfo directory = new DirectoryInfo(outputDirectory);
            foreach (FileInfo file in directory.EnumerateFiles())
                file.Delete();

            string baseName = "Test1";
            StringBuilder pvd = new();
            pvd.AppendLine("<?xml version=\"1.0\"?>");
            pvd.AppendLine("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">");
            pvd.AppendLine("  <Collection>");

            
            int Nt = (int)(Time.Last()/20);
            double delt = Time.Last()/Nt;
            for (int framenum = 0; framenum <= Nt; framenum++)
            {
                double t = framenum*delt;
                double currentTime = t;
                double[] cellS = [];
                double[] cellP = [];
                ExportStepToVtu(outputDirectory, baseName, framenum, Nx, Ny, Nz, Dx, Dy, Dz, zTop, Actnum, cellP, cellS);
                // Direct reference matching your exact naming convention
                string vtuFileName = $"{baseName}_{framenum}.vtu";
                pvd.AppendLine($"    <DataSet timestep=\"{currentTime}\" group=\"\" part=\"0\" file=\"{vtuFileName}\"/>");
            }

            pvd.AppendLine("  </Collection>");
            pvd.AppendLine("</VTKFile>");

            string pvdPath = Path.Combine(outputDirectory, $"{baseName}.pvd");
            File.WriteAllText(pvdPath, pvd.ToString());
        }

        public static double GetZcornValue(double[] zcorn, int nx, int ny, int nz,
                                   int i, int j, int k, int corner)
        {
            // Deconstruct the 0-7 corner index into local binary offsets
            int localX = (corner == 1 || corner == 3 || corner == 5 || corner == 6) ? 1 : 0;
            int localY = (corner == 2 || corner == 3 || corner == 6 || corner == 7) ? 1 : 0;
            int localZ = (corner >= 4) ? 1 : 0;

            // Standard ECLIPSE/OPM flat array offset formula
            int index = (k * 2 * ny * 2 * nx) +
                        (localZ * ny * 2 * nx) +
                        (j * 2 * nx) +
                        (localY * nx * 2) +
                        (i * 2) +
                        localX;

            return zcorn[index];
        }

    }
}
