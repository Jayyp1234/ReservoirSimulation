using EclipseDeckReader;
using Microsoft.CodeAnalysis;
using SepalSolver;
using System.Text;
using static ReservoirSimulator.Math;
using static ReservoirSimulator.ReservoirJacobian;
using static TorchSharp.torch.optim.lr_scheduler.impl;

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
            Utility.BuildStructuredConnections(Grids, Nx, Ny, Nz);
            //Utility.ComputeFaceTransmissibility(Grids);
            

            // 1. Build EQUIL Pressure-Depth Tables for every (EQUILNUM, PVTNUM) pair
            EQUILAnalysis[] EQUILs = [.. Data.RawEquilTables.Select(equil => new EQUILAnalysis(equil, incomingUnit))];
            double MaxZ = Grids.Max(grid => grid.Center.Z);
            double MinZ = Grids.Min(grid => grid.Center.Z);

            for (int e = 0; e < EQUILs.Length; e++)
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
                    int nDatum = Max(10, (int)(Abs(equildata.ZDatum - equildata.ZWOC) / Max(1.0, MaxZ - MinZ) * Nz * 10));

                    if (equildata.ZDatum < equildata.ZWOC)
                    {
                        var oil2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC, nDatum, γo);
                        equildata.PoWOC = oil2WOC.Last().P;
                        equildata.PwWOC = equildata.PoWOC - equildata.PcWOC;
                    }
                    else
                    {
                        var water2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC, nDatum, γw);
                        equildata.PwWOC = water2WOC.Last().P;
                        equildata.PoWOC = equildata.PwWOC + equildata.PcWOC;
                    }

                    // Upper leg: Z from MinZ up to ZWOC (only if MinZ < ZWOC)
                    if (MinZ < equildata.ZWOC)
                    {
                        int nUp = Max(10, (int)(Abs(MinZ - equildata.ZWOC) / Max(1.0, MaxZ - MinZ) * Nz));
                        var oilUp = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PoWOC, MinZ, nUp, γo);
                        var watUp = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MinZ, nUp, γw);

                        oilUp.Reverse(); watUp.Reverse();
                        equildata.OilZonePressureDepthPairs.AddRange(oilUp);
                        equildata.WaterZonePressureDepthPairs.AddRange(watUp);
                    }

                    // Lower leg: Z from ZWOC down to MaxZ (only if MaxZ > ZWOC)
                    if (MaxZ > equildata.ZWOC)
                    {
                        int nDown = Max(10, (int)(Abs(MaxZ - equildata.ZWOC) / Max(1.0, MaxZ - MinZ) * Nz));
                        var watDown = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MaxZ, nDown, γw);
                        var oilDown = watDown.Select(pz => new PressureDepthPair(pz.Z, pz.P + equildata.PcWOC));

                        int skip = equildata.WaterZonePressureDepthPairs.Count > 0 ? 1 : 0;
                        equildata.WaterZonePressureDepthPairs.AddRange(watDown.Skip(skip));
                        equildata.OilZonePressureDepthPairs.AddRange(oilDown.Skip(skip));
                    }

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
                if (!grid.isActive) continue;

                int s = DATA.SatnumArray[m];
                int e = DATA.EquinumArray[m];
                int p = DATA.PvtnumArray[m];

                var swoftable = SWOFs[s].swofTable;
                var equildata = EquilPVTCombinations[(e, p)];
                var oiltable = equildata.OilZonePressureDepthPairs;
                var watertable = equildata.WaterZonePressureDepthPairs;

                double poEquil = Utility.Interpolate(oiltable, grid.Center.Z, pair => pair.Z, pair => pair.P);
                double pwEquil = Utility.Interpolate(watertable, grid.Center.Z, pair => pair.Z, pair => pair.P);
                double pc = poEquil - pwEquil;

                double swMin = swoftable.Min(swof => swof.Sw);
                double swMax = swoftable.Max(swof => swof.Sw);
                double pcMin = swoftable.Min(swof => swof.Pc);
                double pcMax = swoftable.Max(swof => swof.Pc);

                if (pc <= pcMin)
                {
                    grid.Sw_n = swMax;
                }
                else if (pc >= pcMax)
                {
                    grid.Sw_n = swMin;
                }
                else
                {
                    // SWOF is ordered by ascending Sw (descending Pc); reverse to search ascending Pc
                    var ascPcTable = swoftable.AsEnumerable().Reverse().ToList();
                    grid.Sw_n = Utility.Interpolate(ascPcTable, pc, swof => swof.Pc, swof => swof.Sw);
                }

                grid.So_n = 1.0 - grid.Sw_n;

                // Enforce exact algebraic consistency with Residual: Pw = Po - Pcow(Sw)
                double pcActual = grid.Pcow(grid.Sw_n).Value;
                if (grid.Center.Z >= equildata.ZWOC)
                {
                    grid.Pw_n = pwEquil;
                    grid.Po_n = pwEquil + pcActual;
                }
                else
                {
                    grid.Po_n = poEquil;
                    grid.Pw_n = poEquil - pcActual;
                }

                ADiff<Vector4d> meanP = grid.So_n * grid.Po_n + grid.Sw_n * grid.Pw_n;
                ADiff<Vector4d> Er_i = grid.Er(meanP);
                ADiff<Vector4d> Bo_i = grid.Bo(grid.Po_n);
                ADiff<Vector4d> Bw_i = grid.Bw(grid.Pw_n);

                grid.ErSo_Bo_n = (Er_i * grid.So_n / Bo_i).Value;
                grid.ErSw_Bw_n = (Er_i * grid.Sw_n / Bw_i).Value;
            }
        }

        void InitialGuess(double[] xs, double dt, bool isRetry = false, bool wellEventOccurred = false)
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
            double Residual(double[] xnp1, double dt)
            {

                //Utility.checkconnectivity(Grids, JacValues.Length);
                Array.Clear(JacValues, 0, JacValues.Length);
                Array.Clear(Res, 0, Res.Length);
                funcall++;


                // =========================================================================
                // 1. HIGH CPU USAGE: Batched Parallel Assembly with Strict Row Ownership
                // =========================================================================
                Parallel.For(0, cores, coreId =>
                {
                    int start = coreId * chunkSize;
                    int end = (coreId == cores - 1) ? Ngrids : start + chunkSize;

                    //for (int m = 0; m < Ngrids; m++)
                    for (int m = start; m < end; m++)
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
                        var Kro_i = mgrid.Krow(Sw_i);
                        var Krw_i = mgrid.Krw(Sw_i);
                        var μo_i = mgrid.μo(Po_i);
                        var μw_i = mgrid.μw(Pw_i);
                        var γo_i = mgrid.γo(Po_i);
                        var γw_i = mgrid.γw(Pw_i);
                        var λo_i = Kro_i / (μo_i * Bo_i);
                        var λw_i = Krw_i / (μw_i * Bw_i);

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
                                var Kro_j = ngrid.Krow(Sw_j);
                                var Krw_j = ngrid.Krw(Sw_j);
                                var μo_j = ngrid.μo(Po_j);
                                var μw_j = ngrid.μw(Pw_j);
                                var γo_j = ngrid.γo(Po_j);
                                var γw_j = ngrid.γw(Pw_j);
                                var λo_j = Kro_j / (μo_j * Bo_j);
                                var λw_j = Krw_j / (μw_j * Bw_j);

                                // Hydrostatic potential differences
                                var γo_avg = 0.5 * (γo_i + γo_j);
                                var γw_avg = 0.5 * (γw_i + γw_j);
                                var DΦo = Po_i - Po_j - γo_avg * dz;
                                var DΦw = Pw_i - Pw_j - γw_avg * dz;

                                // Upstream phase mobilities
                                var λo_up = Utility.SmoothUpstream(DΦo, λo_i, λo_j);
                                var λw_up = Utility.SmoothUpstream(DΦw, λw_i, λw_j);

                                // Outflow flux from cell i to cell j
                                var FluxOil = transDt * λo_up * DΦo;
                                var FluxWater = transDt * λw_up * DΦw;

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


            double t = 0.0, dt = 0.1, dtCarry = 0.1;
            DateTime startDate = DATA.StartDate, currentSimDate;
            List<DateTime> Interval = [startDate, .. ResultTime.Select(time => startDate.AddDays(time))];
            foreach (var well in Wells)
                Interval.AddRange(well.GetEventTimes());
            Interval = [.. Interval.Distinct().OrderBy(x => x)];

            Console.WriteLine($"""
                Total Grids = {Ngrids}, Active Grids = {NumActiveGrids}
                ======================================================================
                Starting simulation
                """);

            double L = 0, tnp1, rnorm = 0;
            bool staterejected = false, isComplete;
            int NumActiveWells = 0;
            Time = [t];

            for (int i = 1; i < Interval.Count; i++)
            {
                L = Interval[i].Subtract(startDate).TotalDays;
                isComplete = false;

                while (!isComplete)
                {
                    double remaining = L - t;
                    bool hitMilestone = dtCarry >= remaining - 1e-13;
                    dt = hitMilestone ? remaining : dtCarry;

                    tnp1 = t + dt;
                    currentSimDate = startDate.AddDays(tnp1);
                    Console.WriteLine($"\n[Info] Simulation Time: {tnp1:F3} days | Current Date: {currentSimDate:yyyy-MM-dd} | dt: {dt:F4} d");

                    if (!staterejected)
                    {
                        NumActiveWells = 0;
                        foreach (var well in Wells)
                        {
                            well.UpdateControl(currentSimDate);
                            if (well.isActive) well.ActiveIndex = NumActiveWells++;
                        }

                        foreach (var well in Wells)
                            well.UpdateCompletion(currentSimDate);

                        // Rebuild CSR and run PARDISO symbolic factorization ONLY when active wells or completions change
                        if (Wells.Any(well => well.ActiveStatusChanged) || Wells.Any(well => well.CompletionChanged))
                        {
                            Console.WriteLine($"[Info] Number of Active Wells: {NumActiveWells}");
                            Console.WriteLine($"[Info] Well Completion/Active Status Changed at t = {tnp1:F3} days -> Rebuilding CSR & Symbolic Factorization");
                            solver?.Dispose(); // <-- CRITICAL: triggers PARDISO phase = -1 to free native L/U buffers!
                            solver = null;
                            (totalDofs, col_index, row_ptr, JacValues) = ReservoirJacobian.Create(2, Grids, Wells);
                            (cellOffsets, faceOffsets, wellOffsets) = ReservoirJacobian.PrecomputeOffsets(2, row_ptr, col_index, Grids, Wells);
                            solver = new MklPardisoSolver(totalDofs, row_ptr, col_index);
                            xs = new double[totalDofs];
                            Res = new double[totalDofs];
                            dx = new double[totalDofs];
                        }
                    }

                    // Prepare initial guess x^(0) from converged state n
                    InitialGuess(xs, dt);
                    staterejected = false;

                    Console.WriteLine($"""
                        Time: {tnp1:E4} days
                         iter |   Residual Norm
                        ------+------------------
                        """);

                    int iter;
                    bool isConverged = false;
                    rnorm = Residual(xs, dt);
                    double prevRnorm = rnorm;

                    for (iter = 1; iter < 10; iter++)
                    {
                        // Direct solve via MKL PARDISO: J * dx = Res
                        solver.Solve(JacValues, Res, dx);

                        // Update primary variables: x^{k+1} = x^k - dx
                        (double maxDPo, double maxDSw) = ApplyPrimaryVariableUpdate(xs, dx, NumActiveGrids, NumActiveWells);
                        rnorm = Residual(xs, dt);

                        Console.WriteLine($"  {iter,2}  |  {rnorm:E4}");

                        // Early exit if Newton diverges or hits NaN
                        if (double.IsNaN(rnorm) || (iter > 2 && rnorm > 10.0 * prevRnorm))
                            break;
                        prevRnorm = rnorm;

                        // Converged when normalized residual AND state step are both within tolerance
                        if (rnorm < 1e-6 && maxDSw < 1e-4 && maxDPo < 1e-2)
                        {
                            isConverged = true;
                            break;
                        }
                    }

                    // Cut dt on Newton non-convergence
                    if (!isConverged)
                    {
                        dtCarry = 0.5 * dt;
                        if (dtCarry < 1e-12) throw new Exception("Time step is too small");
                        staterejected = true;
                        Console.WriteLine("""
                            ============================
                            Rejected (Non-Convergence)
                            ============================
                            """);
                        continue;
                    }

                    // Validate maximum per-step state changes
                    var (maxDP, maxDS, i_p, i_s) = CheckSolution(xs);
                    if (maxDP > 100 || maxDS > 0.15)
                    {
                        int j = maxDP > 100 ? i_p : i_s;
                        double limit = maxDP > 100 ? 100 : 0.15;
                        double change = maxDP > 100 ? maxDP : maxDS;
                        Console.ForegroundColor = ConsoleColor.Yellow;
                        Console.WriteLine($"\n[Warning] Time Step Rejected: State Change Limit Exceeded.");
                        Console.ResetColor();
                        Console.WriteLine($"    ↳ Location       : Cell Index {j}");
                        Console.WriteLine($"    ↳ Violation      : {(maxDP > 100 ? "ΔP" : "ΔS")} = {change:F4} (Max Allowed: {limit:F4})");
                        Console.WriteLine($"    ↳ Current State  : P = {Grids[j].Po_n:F1} psi | Sw = {Grids[j].Sw_n:F4}");
                        Console.WriteLine($"    ↳ Attempted State: P = {Grids[j].Po_np1:F1} psi | Sw = {Grids[j].Sw_np1:F4}");
                        Console.WriteLine($"    ↳ Solver Action  : Rolling back to t = {t:F3} days.");
                        Console.WriteLine($"    ↳ Time-Step Chop : Reducing Δt from {dt:F2} to {0.5 * dt:F2} days (Factor: 0.5).\n");
                        dtCarry = 0.5 * dt;
                        staterejected = true;
                        continue;
                    }

                    // Validate Well BHP constraints (reuses existing PARDISO factorization on constraint switch)
                    foreach (var well in Wells)
                    {
                        if (!well.isActive) continue;

                        if (well.WellType == WellType.PRODUCER &&
                            well.ConstraintType != ConstraintType.MIN_BHP &&
                            well.Pwf < well.CurrentControl.MinBhpLimit)
                        {
                            well.ConstraintType = ConstraintType.MIN_BHP;
                            Console.WriteLine($"""
                                ===============================================================
                                Rejected (Minimum Pressure Violated) @ {well.Name}
                                ===============================================================
                                """);
                            staterejected = true;
                        }

                        if (well.WellType == WellType.INJECTOR &&
                            well.ConstraintType != ConstraintType.MAX_BHP &&
                            well.Pwf > well.CurrentControl.MaxBhpLimit)
                        {
                            well.ConstraintType = ConstraintType.MAX_BHP;
                            Console.WriteLine($"""
                                ================================================
                                Rejected (Maximum Pressure Violated) @ {well.Name}
                                ================================================
                                """);
                            staterejected = true;
                        }
                    }
                    if (staterejected) continue;

                    // Commit verified state n+1 -> n
                    t = tnp1;
                    Time.Add(t);
                    foreach (var grid in Grids)
                    {
                        if (!grid.isActive) continue;
                        grid.ErSo_Bo_n = grid.ErSo_Bo_np1.Value;
                        grid.ErSw_Bw_n = grid.ErSw_Bw_np1.Value;
                        grid.Po_n = grid.Po_np1; grid.Sw_n = grid.Sw_np1;
                        grid.Pw_n = grid.Pw_np1; grid.So_n = grid.So_np1;
                        grid.SaveResult(po: grid.Po_n, pw: grid.Pw_n, so: grid.So_n, sw: grid.Sw_n);
                    }
                    foreach (var well in Wells)
                    {
                        if (!well.isActive) continue;
                        well.Pwf_n = well.Pwf.Value;
                        well.Qo_n = well.Qo.Value;
                        well.Qw_n = well.Qw.Value;
                        well.SaveResult(pwf: well.Pwf_n, qo: well.Qo_n, qw: well.Qw_n);
                    }

                    isComplete = Abs(t - L) < 1e-13;

                    // Adaptive time-stepping update preserving pre-clamped dtCarry across milestones
                    double nextDt = dt;
                    if (iter < 3) nextDt *= 1.25;
                    if (iter > 7) nextDt *= 0.7;

                    dtCarry = hitMilestone ? Clamp(Max(dtCarry, nextDt), 1e-12, 100.0) : Clamp(nextDt, 1e-12, 100.0);
                }
            }
        }
    }
}
