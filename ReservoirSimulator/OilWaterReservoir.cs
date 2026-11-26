using EclipseDeckReader;
using Microsoft.CodeAnalysis;
using SepalSolver;
using System.Globalization;
using static ReservoirSimulator.Math;
using static ReservoirSimulator.ReservoirJacobian;

namespace ReservoirSimulator
{
    public class OilWaterReservoir : IReservoir<Vector4d>
    {
        public int funcall;
        public bool UseUpstreamWeightingForInjection = false;
        public List<double> Time, SweepEff;
        readonly Func<ADiff<Vector4d>, ADiff<Vector4d>> Er;
        (double maxDP, double maxDS, int i_p, int i_s) CheckSolution(double[] xs)
        {
            int i_p = 0, i_s = 0; double dp, maxDp = 0, ds, maxDs = 0;
            for (int i = 0; i < Ngrids; i++)
            {
                dp = Grids[i].Po_np1 - Grids[i].Po_n;
                if (Abs(dp) > maxDp) { maxDp = Abs(dp); i_p = i; }
                ds = Grids[i].Sw_np1 - Grids[i].Sw_n;
                if (Abs(ds) > maxDs) { maxDs = Abs(ds); i_s = i; }
            }
            return (maxDp, maxDs, i_p, i_s);
        }

        readonly int Nx, Ny, Nz, NxNy, Ngrids, NumActiveGrids, NumWells;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, zTop;
        double[] Coord, Zcorn;
        readonly double Prefr, cr, γo0, γw0;
        Grid<Vector4d>[] Grids;
        List<Grid<Vector4d>> ActiveGrids;
        public List<Well<Vector4d>> Wells;
        Aquifer Aquifer;
        Dictionary<(int, int), EQUILPressureDepth> EquilPVTCombinations = [];
        EclipseDataDeck DATA;
        PVTWAnalysis[] PVTWs;
        PVDOAnalysis[] PVDOs;
        SWOFAnalysis[] SWOFs;

        public List<Well<Vector4d>> WELLS => Wells;

        public Grid<Vector4d>[] GRIDS => Grids;

        bool IReservoir.UseUpstreamWeightingForInjection { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

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

            // NTG and MULTPV default to 1 for every cell (set at DIMENS); a wrong-sized array is ignored rather than misread
            double[]? ntg = Data.NtgArray.Count == Ngrids ? Data.NtgArray.ToArray() : null;
            double[]? multpv = Data.MultpvArray.Count == Ngrids ? Data.MultpvArray.ToArray() : null;
            Grids = Zcorn is not null ? Utility.ProcessGrid<Vector4d>(Nx, Ny, Nz, Coord, Zcorn, Φ, Kx, Ky, Kz, ntg, multpv) :
                Utility.ProcessGrid<Vector4d>(Nx, Ny, Nz, Dx, Dy, Dz, zTop, Φ, Kx, Ky, Kz, ntg, multpv);

            // Process PVTW and PVDO tables.
            PVTWs = [.. Data.RawPvtwTables.Select(pvtw => new PVTWAnalysis(pvtw, incomingUnit))];
            PVDOs = [.. Data.RawPvdoTables.Select(pvdo => new PVDOAnalysis(pvdo, incomingUnit))];
            for (int i = 0; i < PVTWs.Length; i++)
            {
                PVTWs[i].FitProperties();
                PVDOs[i].FitProperties();
            }

            // Process SWOF tables.
            SWOFs = [.. Data.RawSwofTables.Select(swof => new SWOFAnalysis(swof, incomingUnit))];
            foreach (var swof in SWOFs)
            { 
                swof.FitProperties(); 
            }

            // ROCK
            cr = Data.RockCompressibility;
            Prefr = Utility.Convert(Data.RockReferencePressure, "Pressure", incomingUnit, Units.FIELD);
            Er = P => Exp(cr * (P - Prefr));

            // DENSITIES
            γo0 = Utility.Convert(Data.SurfaceOilDensity, "Density", incomingUnit, Units.FIELD) / 144;
            γw0 = Utility.Convert(Data.SurfaceWaterDensity, "Density", incomingUnit, Units.FIELD) / 144;

            // MINPV: cells holding less pore volume than the deck's minimum are made inactive, as in ECLIPSE.
            // A copy of ACTNUM is changed, so the deck itself keeps what the user wrote.
            List<int> actnum = [.. Data.ActnumArray];
            if (Data.MinPoreVolume > 0)
            {
                double minPv = Utility.Convert(Data.MinPoreVolume, "Reservoir Volume", incomingUnit, Units.FIELD);
                int removed = 0;
                for (int m = 0; m < Ngrids; m++)
                    if (actnum[m] > 0 && Grids[m].PoreVolume < minPv) { actnum[m] = 0; removed++; }
                Console.WriteLine($"      [Grid] MINPV {Data.MinPoreVolume}: {removed} cell(s) with less pore volume made inactive.");
            }

            // REGION
            NumActiveGrids = Utility.SetActiveGrids(Grids, actnum, Er);
            ActiveGrids = [.. Grids.Where(g => g.isActive)];
            foreach (var grid in ActiveGrids)
            {
                int satnum = Data.SatnumArray[grid.GlobalIndex];
                grid.SWOF = SWOFs[satnum];
                int pvtnum = Data.PvtnumArray[grid.GlobalIndex];
                grid.PVTW = PVTWs[pvtnum];
                grid.PVDO = PVDOs[pvtnum];
            }

            // TRANSMISSIBILITY
            Utility.BuildStructuredConnections(Grids, Nx, Ny, Nz);
            Utility.ApplyTransmissibilityMultipliers(Grids, Nx, Ny, Nz, Data);


            // 1. Build EQUIL Pressure-Depth Tables for every (EQUILNUM, PVTNUM) pair
            EQUILAnalysis[] EQUILs = [.. Data.RawEquilTables.Select(equil => new EQUILAnalysis(equil, incomingUnit))];
            double MaxZ = Grids.Max(grid => grid.Center.Z);
            double MinZ = Grids.Min(grid => grid.Center.Z);

            for (int e = 0; e < EQUILs.Length; e++)
            {
                for (int p = 0; p < PVDOs.Length; p++)
                {
                    var equil = EQUILs[e]; var pvdo = PVDOs[p]; var pvtw = PVTWs[p];
                    Func<double, double> γo = P => γo0 / pvdo.Result.EvaluateBo(P),
                                         γw = P => γw0 / pvtw.Result.EvaluateBw(P);

                    var equildata = new EQUILPressureDepth(equil);
                    int nDatum = Max(10, (int)(Abs(equildata.ZDatum - equildata.ZWOC) / Max(1.0, MaxZ - MinZ) * Nz * 10));

                    if (equildata.ZDatum < equildata.ZWOC)
                    {
                        // Downward integration from ZDatum -> ZWOC: target pressure is Last().P
                        var oil2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC, nDatum, γo);
                        equildata.PoWOC = oil2WOC.Last().P;
                        equildata.PwWOC = equildata.PoWOC - equildata.PcWOC;
                    }
                    else
                    {
                        // Upward integration from ZDatum -> ZWOC: target pressure is First().P
                        var water2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC, nDatum, γw);
                        equildata.PwWOC = water2WOC.First().P;
                        equildata.PoWOC = equildata.PwWOC + equildata.PcWOC;
                    }

                    // Upper leg: Z from MinZ down to ZWOC (only if MinZ < ZWOC)
                    if (MinZ < equildata.ZWOC)
                    {
                        int nUp = Max(10, (int)(Abs(MinZ - equildata.ZWOC) / Max(1.0, MaxZ - MinZ) * Nz));
                        // Already returned in ascending Z order [MinZ .. ZWOC] — do NOT reverse
                        var oilUp = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PoWOC, MinZ, nUp, γo);
                        var watUp = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MinZ, nUp, γw);

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
            Wells = [.. Data.Wells.Select(well => new Well<Vector4d>(well.Value, Nx, Ny, Nz, Grids, incomingUnit))];
            foreach (var well in Wells.Where(w => w.PerfNatIndex.Count == 0))
                Console.WriteLine($"      [Warning] Well '{well.Name}' has no completion in an active cell; it is left out of the simulation.");
            Wells.RemoveAll(w => w.PerfNatIndex.Count == 0);
            NumWells = Wells.Count;
        }
        public void Initialize()
        {
            funcall = 0;
            foreach(var grid in Grids)
            {
                if (!grid.isActive) continue;

                int s = DATA.SatnumArray[grid.GlobalIndex];
                int e = DATA.EquinumArray[grid.GlobalIndex];
                int p = DATA.PvtnumArray[grid.GlobalIndex];

                var swoftable = SWOFs[s].swofTable;
                var equildata = EquilPVTCombinations[(e, p)];
                var oiltable = equildata.OilZonePressureDepthPairs;
                var watertable = equildata.WaterZonePressureDepthPairs;

                double z = grid.Center.Z;
                double poEquil = Utility.Interpolate(oiltable, z, pair => pair.Z, pair => pair.P);
                double pwEquil = Utility.Interpolate(watertable, z, pair => pair.Z, pair => pair.P);

                double pc = poEquil - pwEquil;

                double swMin = SWOFs[s].Result.SwMin, swMax = SWOFs[s].Result.SwMax;
                double pcMin = SWOFs[s].Result.PcMin, pcMax = SWOFs[s].Result.PcMax;

                if (pc <= pcMin) grid.Sw_n = swMax;
                else if (pc >= pcMax) grid.Sw_n = swMin;
                else
                    grid.Sw_n = Utility.Interpolate(swoftable, pc, swof => swof.Pc, swof => swof.Sw);
                

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

                var meanP = grid.So_n * grid.Po_n + grid.Sw_n * grid.Pw_n;
                var Er_i = grid.Er(meanP);
                var Bo_i = grid.Bo(grid.Po_n);
                var Bw_i = grid.Bw(grid.Pw_n);

                grid.ErSo_Bo_n = (Er_i * grid.So_n / Bo_i).Value;
                grid.ErSw_Bw_n = (Er_i * grid.Sw_n / Bw_i).Value;
                grid.SaveResult(po: grid.Po_n, pw: grid.Pw_n, 
                                so: grid.So_n, sw: grid.Sw_n);
            }

            foreach (var well in Wells)
            {
                well.Pwf_n = well.PerfNatIndex.Select(m => Grids[m].Po_n).Average();
                well.SaveResult(pwf: well.Pwf_n, qo: well.Qo_n, qw: well.Qw_n);
            }
        }
        void InitialGuess(double[] xs, double dt, bool isRetry = false, bool wellEventOccurred = false)
        {
            int indx = 0, n = Time.Count - 1;
            bool canExtrapolate = n >= 2 && !isRetry && !wellEventOccurred;
            double omega = canExtrapolate ? Min(1.0, dt / (Time[n] - Time[n - 1])) : 0.0;

            // 1. Reservoir Grid Unknowns (Po, Sw)
            foreach (var mgrid in Grids)
            {
                if (!mgrid.isActive) continue;

                double poGuess = mgrid.Po_n;
                if (canExtrapolate)
                    poGuess += omega * (mgrid.Po_n - mgrid.PO[n - 1]);

                xs[indx++] = Max(14.7, poGuess); // Po index
                xs[indx++] = mgrid.Sw_n;         // Sw index
            }

            // 2. Well Unknowns (Pwf, Q) for all wells
            foreach (var well in Wells)
            {
                double pGuess = well.Pwf_n;
                double qGuess = well.IsWaterTarget ? well.Qw_n
                              : well.IsOilTarget ? well.Qo_n
                              : well.IsGasTarget ? well.Qg_n
                                                   : well.Ql_n;

                if (well.isStopped)
                    qGuess = 0.0;
                else if (well.WellType == WellType.PRODUCER && well.CurrentProducerControl is { } prod)
                {
                    switch (well.ConstraintType)
                    {
                        case ConstraintType.MIN_BHP:
                            pGuess = prod.MinBhpLimit;
                            break;
                        case ConstraintType.ORATE:
                            qGuess = prod.OilRateTarget;
                            break;
                        case ConstraintType.WRATE:
                            qGuess = prod.WaterRateTarget;
                            break;
                        case ConstraintType.GRATE:
                            qGuess = prod.GasRateTarget;
                            break;
                        case ConstraintType.LRATE:
                            qGuess = prod.LiquidRateTarget;
                            break;
                    }
                }
                else if (well.WellType == WellType.INJECTOR && well.CurrentInjectorControl is { } inj)
                {
                    switch (well.ConstraintType)
                    {
                        case ConstraintType.MAX_BHP:
                            pGuess = inj.MaxBhpLimit;
                            break;
                        case ConstraintType.RATE:
                            qGuess = inj.TargetRate;
                            break;
                    }
                }

                xs[indx++] = pGuess;
                xs[indx++] = qGuess;
            }
        }
        public static (double maxPd, double maxSd) ApplyPrimaryVariableUpdate(
            double[] x, double[] dx, int numActiveGrids, int numWells,
            double dPmax = 500.0,   // max allowable pressure change per Newton step [psi]
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
            for (int w = 0; w < numWells; w++)
            {
                int idxP = wellStart + 2 * w + 0; // Well pressure [psia]
                int idxQ = wellStart + 2 * w + 1; // Well rate [stb/day]
                double dP = dx[idxP], dQ = dx[idxQ];
                x[idxP] -= dP; // 
                x[idxQ] -= dQ; // 
            }

            return (maxPd, maxSd);
        }

        public void Simulate(double[] ResultTime = null)
        {
            if(ResultTime == null) ResultTime = Utility.GetCumulativeTimeArray(DATA.TimeSteps);
            int cores = Environment.ProcessorCount, chunkSize = Ngrids / cores;

            (int totalDofs, int[] col_index, int[] row_ptr, double[] JacValues) = ReservoirJacobian.Create(2, Grids, Wells);
            ReservoirJacobian.PrecomputeOffsets(2, row_ptr, col_index, Grids, Wells);
            double[] xs = new double[totalDofs], Res = new double[totalDofs], dx = new double[totalDofs];
            MklPardisoSolver solver = new MklPardisoSolver(totalDofs, row_ptr, col_index); 

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
                        mgrid.ErSo_Bo_np1 = Er_i * So_i / Bo_i; 
                        mgrid.ErSw_Bw_np1 = Er_i * Sw_i / Bw_i;

                        // Accumulation Residual: R = Vp * (State^{n+1} - State^n)
                        var ResOil = mgrid.PoreVolume * (mgrid.ErSo_Bo_np1 - mgrid.ErSo_Bo_n);
                        var ResWater = mgrid.PoreVolume * (mgrid.ErSw_Bw_np1 - mgrid.ErSw_Bw_n);

                        // Direct write to Residual Vector
                        Res[row_o] += ResOil.Value;
                        Res[row_w] += ResWater.Value;

                        // Direct write to Accumulation Jacobian (Block i, i)
                        var cOff = mgrid.Offset;
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
                                var fOff = face.Offset[c];

                                // Block (i, i) -> Derivatives wrt Cell i variables (slots 0 and 1)
                                JacValues[cOff.Diag[0]] += FluxOil.Derivative[0];   // row_o, col_p_i
                                JacValues[cOff.Diag[1]] += FluxOil.Derivative[1];   // row_o, col_sw_i
                                JacValues[cOff.Diag[2]] += FluxWater.Derivative[0]; // row_w, col_p_i
                                JacValues[cOff.Diag[3]] += FluxWater.Derivative[1]; // row_w, col_sw_i

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
                
                for (int wellidx = 0; wellidx < NumWells; wellidx++)
                {
                    var well = Wells[wellidx];
                    well.Qo_np1 = 0.0; well.Qw_np1 = 0.0; well.Ql_np1 = 0.0;
                    well.Qo = 0.0; well.Qw = 0.0;
                    var wOff = well.Offset;
                    if (!well.isActive)
                    {
                        JacValues[wOff.RowFlow_dPwell] += 1;
                        JacValues[wOff.RowConst_dQwell] += 1;
                        continue;
                    }

                    int rowflow = 2 * NumActiveGrids + 2 * wellidx + 0;  // Rate definition equation
                    int rowconst = 2 * NumActiveGrids + 2 * wellidx + 1; // Constraint equation

                    var pwell = new ADiff<Vector4d>(xnp1[rowflow], new(1, 0, 0, 0));
                    var qwell = new ADiff<Vector4d>(xnp1[rowconst], new(0, 1, 0, 0));
                    var FlowRes = qwell * dt;

                    int perfCount = well.PerfNatIndex.Count;
                    bool isWaterTarget = well.IsWaterTarget;
                    bool isOilTarget = well.IsOilTarget;
                    bool isLiqTarget = well.IsLiqTarget;
                    bool isResTarget = well.IsResTarget;

                    switch (well.WellType)
                    {
                        case WellType.PRODUCER:
                            {
                                for (int p = 0; p < perfCount; p++)
                                {
                                    int m = well.PerfNatIndex[p];
                                    var mgrid = Grids[m];
                                    if (!well.PerfIsOpen[p] || !mgrid.isActive) continue;

                                    int i = mgrid.ActiveIndex;
                                    int row_o = 2 * i + 0;
                                    int row_w = 2 * i + 1;

                                    var Po_m = new ADiff<Vector4d>(xnp1[row_o], new(0, 0, 1, 0));
                                    var Sw_m = new ADiff<Vector4d>(xnp1[row_w], new(0, 0, 0, 1));
                                    var Pw_m = Po_m - mgrid.Pcow(Sw_m);
                                    var So_m = 1.0 - Sw_m;

                                    var Bo_m = mgrid.Bo(Po_m);
                                    var Bw_m = mgrid.Bw(Pw_m);
                                    var Kro_m = mgrid.Krow(Sw_m);
                                    var Krw_m = mgrid.Krw(Sw_m);
                                    var μo_m = mgrid.μo(Po_m);
                                    var μw_m = mgrid.μw(Pw_m);
                                    var γo_m = mgrid.γo(Po_m);
                                    var γw_m = mgrid.γw(Pw_m);

                                    var γo_well = mgrid.γo(pwell);
                                    var γw_well = mgrid.γw(pwell);

                                    var γo_avg = 0.5 * (γo_m + γo_well);
                                    var γw_avg = 0.5 * (γw_m + γw_well);
                                    double dz = mgrid.Center.Z - well.Zref;

                                    var DΦo = Po_m - pwell - γo_avg * dz;
                                    var DΦw = Pw_m - pwell - γw_avg * dz;

                                    var WIo = mgrid.WI * (Kro_m / (μo_m * Bo_m));
                                    var WIw = mgrid.WI * (Krw_m / (μw_m * Bw_m));

                                    var rate_oil = WIo * DΦo; well.Qo += rate_oil;
                                    var rate_water = WIw * DΦw; well.Qw += rate_water;

                                    // Reservoir residuals & Jacobian entries (outflow is +rate * dt)
                                    Res[row_o] += rate_oil.Value * dt;
                                    Res[row_w] += rate_water.Value * dt;

                                    JacValues[wOff.Rm_dPwell[p][0]] += rate_oil.Derivative[0] * dt;
                                    JacValues[wOff.Rm_dPwell[p][1]] += rate_water.Derivative[0] * dt;

                                    var cOff = mgrid.Offset;
                                    JacValues[cOff.Diag[0]] += rate_oil.Derivative[2] * dt;
                                    JacValues[cOff.Diag[1]] += rate_oil.Derivative[3] * dt;
                                    JacValues[cOff.Diag[2]] += rate_water.Derivative[2] * dt;
                                    JacValues[cOff.Diag[3]] += rate_water.Derivative[3] * dt;

                                    // Per-perforation grid derivatives for Row Q
                                    var perfRate = isWaterTarget ? rate_water
                                                    : isOilTarget ? rate_oil
                                                    : isLiqTarget ? rate_oil + rate_water
                                                    : rate_oil * Bo_m + rate_water * Bw_m;

                                    FlowRes -= perfRate * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][0]] -= perfRate.Derivative[2] * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][1]] -= perfRate.Derivative[3] * dt;
                                }

                                // close production equation
                                Res[rowflow] += FlowRes.Value;
                                JacValues[wOff.RowFlow_dPwell] += FlowRes.Derivative[0];
                                JacValues[wOff.RowFlow_dQwell] += FlowRes.Derivative[1];
                                break;
                            }

                        case WellType.INJECTOR:
                            {
                                for (int p = 0; p < perfCount; p++)
                                {
                                    int m = well.PerfNatIndex[p];
                                    var mgrid = Grids[m];
                                    if (!well.PerfIsOpen[p] || !mgrid.isActive) continue;

                                    int i = mgrid.ActiveIndex;
                                    int row_o = 2 * i + 0;
                                    int row_w = 2 * i + 1;
                                    double dz = mgrid.Center.Z - well.Zref;

                                    var Po_m = new ADiff<Vector4d>(xnp1[row_o], new(0, 0, 1, 0));
                                    var Sw_m = new ADiff<Vector4d>(xnp1[row_w], new(0, 0, 0, 1));
                                    var Pw_m = Po_m - mgrid.Pcow(Sw_m);
                                    var So_m = 1 - Sw_m;
                                    var Bo_m = mgrid.Bo(Po_m);
                                    var Bw_m = mgrid.Bw(Pw_m);
                                    var Kro_m = mgrid.Krow(Sw_m);
                                    var Krw_m = mgrid.Krw(Sw_m);
                                    var μo_m = mgrid.μo(Po_m);
                                    var μw_m = mgrid.μw(Pw_m);
                                    var γo_m = mgrid.γo(Po_m);
                                    var γw_m = mgrid.γw(Pw_m);

                                    var γw_well = mgrid.γw(pwell);
                                    var γw_avg = 0.5 * (γw_m + γw_well);

                                    var DΦw = Pw_m - pwell - γw_avg * dz;
                                    var lambdaTotal = (Krw_m / μw_m + Kro_m / μo_m) / Bw_m;
                                    if (UseUpstreamWeightingForInjection)
                                        lambdaTotal = mgrid.Krw(1) / (μw_m * Bw_m);

                                    // Phase-weighted total mobility in perforated cell
                                    var WIw = mgrid.WI * lambdaTotal;
                                    var rate_water = WIw * DΦw; well.Qw += rate_water;

                                    // Reservoir water equation (injection is inflow: rate * dt)
                                    Res[row_w] += rate_water.Value * dt;

                                    JacValues[wOff.Rm_dPwell[p][1]] += rate_water.Derivative[0] * dt;

                                    var cOff = mgrid.Offset;
                                    JacValues[cOff.Diag[2]] += rate_water.Derivative[2] * dt;
                                    JacValues[cOff.Diag[3]] += rate_water.Derivative[3] * dt;

                                    // Per-perforation grid derivatives for Row Q
                                    var perfRate = isWaterTarget ? rate_water
                                                    : isLiqTarget ? rate_water
                                                    : rate_water * Bw_m;

                                    FlowRes -= perfRate * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][0]] -= perfRate.Derivative[2] * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][1]] -= perfRate.Derivative[3] * dt;
                                }

                                // Close injection equation
                                Res[rowflow] += FlowRes.Value;
                                JacValues[wOff.RowFlow_dPwell] += FlowRes.Derivative[0];
                                JacValues[wOff.RowFlow_dQwell] += FlowRes.Derivative[1];
                                break;
                            }

                    }

                    well.Qw_np1 = well.Qw.Value;
                    well.Qo_np1 = well.Qo.Value;
                    well.Ql_np1 = well.Ql.Value;
                    well.Pwf_np1 = pwell.Value;
                    // Close the Well Constraint Equation (Row p)
                    ADiff<Vector4d> wellConstraint = well.Constraint(pwell, qwell);
                    Res[rowconst] += wellConstraint.Value;
                    JacValues[wOff.RowConst_dPwell] += wellConstraint.Derivative[0];
                    JacValues[wOff.RowConst_dQwell] += wellConstraint.Derivative[1];

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
                for (int wellidx = 0; wellidx < NumWells; wellidx++)
                {
                    var well = Wells[wellidx];
                    if (!well.isActive) continue;

                    bool isProducer = well.WellType == WellType.PRODUCER;
                    well.Zref = isProducer ? well.HighestZref : well.LowestZref;

                    int row0 = wellOffset + 2 * wellidx + 0;
                    int row1 = wellOffset + 2 * wellidx + 1;

                    // Resolve active rate target & BHP limit from the role-specific control record
                    double targetRate = 0.0;
                    double bhpLimit = 14.7;

                    if (!well.isStopped)
                    {
                        if (isProducer && well.CurrentProducerControl is { } prod)
                        {
                            bhpLimit = prod.MinBhpLimit;
                            targetRate = well.ConstraintType switch
                            {
                                ConstraintType.ORATE => prod.OilRateTarget,
                                ConstraintType.WRATE => prod.WaterRateTarget,
                                ConstraintType.GRATE => prod.GasRateTarget,
                                ConstraintType.LRATE => prod.LiquidRateTarget,
                                ConstraintType.RESV => prod.ResvRateTarget,
                                _ => Max(prod.OilRateTarget, Max(prod.WaterRateTarget, prod.LiquidRateTarget))
                            };
                        }
                        else if (!isProducer && well.CurrentInjectorControl is { } inj)
                        {
                            bhpLimit = inj.MaxBhpLimit;
                            targetRate = inj.TargetRate;
                        }
                    }

                    // Current flowing rate provides a physical scale when under BHP control or STOP
                    double currentFlowRate = Abs(well.Qo.Value) + Abs(well.Qw.Value);
                    double rateScale = Max(1.0, Max(Abs(targetRate), currentFlowRate));

                    // Equation 0: Wellbore rate balance / inflow equation [stb/day] -> normalized
                    double err0 = Abs(Res[row0]) / rateScale;

                    // Equation 1: Active control constraint (BHP [psi] vs Rate [stb/day]) -> normalized
                    double ctrlScale = well.ConstraintType is ConstraintType.MAX_BHP or ConstraintType.MIN_BHP
                        ? Max(14.7, Abs(bhpLimit))
                        : rateScale;

                    double err1 = Abs(Res[row1]) / ctrlScale;

                    if (err0 > MaxCNV) MaxCNV = err0;
                    if (err1 > MaxCNV) MaxCNV = err1;
                }

                return MaxCNV;
            }


            double t = 0.0, dt = 0.1, dtCarry = 0.1;
            DateTime startDate = DATA.StartDate, currentSimDate;
            List<DateTime> Interval = [.. new DateTime[] { startDate }
                                .Concat(Wells.SelectMany(well => well.GetEventTimes().Select(date => date.AddSeconds(-1))))
                                .Concat(ResultTime.Select(t => startDate.AddDays(t).AddSeconds(-1)))
                                .Distinct().Where(date=>date>=startDate).Order()];

            Console.WriteLine($"""
                Total Grids = {Ngrids}, Active Grids = {NumActiveGrids}
                ======================================================================
                Starting simulation
                """);

            double L = 0, tnp1, rnorm = 0;
            bool tryingnewtimestep = true, isComplete;
            int controlSwitches = 0;   // control switches made on the current step
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

                    if (tryingnewtimestep)
                    {
                        controlSwitches = 0;
                        foreach (var well in Wells)
                            well.UpdateControl(currentSimDate);

                        foreach (var well in Wells)
                            well.UpdateCompletion(currentSimDate);
                    }

                    // Prepare initial guess x^(0) from converged state n
                    InitialGuess(xs, dt, false, true);

                    Console.WriteLine($"""
                         iter |   Residual Norm
                        ------+------------------
                        """);
                    int iter;
                    bool isConverged = false;
                    double prevRnorm = double.PositiveInfinity;
                    double maxDPo = double.PositiveInfinity;
                    double maxDSw = double.PositiveInfinity;

                    for (iter = 0; iter < 10; iter++)
                    {
                        rnorm = Residual(xs, dt);
                        Console.WriteLine($"  {iter,2}  |  {rnorm:E4}");

                        // Early exit if Newton diverges or hits NaN
                        if (double.IsNaN(rnorm) || (iter > 2 && rnorm > 10.0 * prevRnorm))
                            break;

                        // Converged at current xs (skip step-size requirement on iter 0 only if desired, or require iter > 0)
                        if (rnorm < 1e-6 && (iter == 0 || (maxDSw < 1e-3 && maxDPo < 1e-1)))
                        {
                            isConverged = true;
                            break;
                        }

                        prevRnorm = rnorm;

                        // Direct solve via MKL PARDISO: J * dx = Res
                        solver.Solve(JacValues, Res, dx);

                        // Update primary variables: x^{k+1} = x^k - dx
                        (maxDPo, maxDSw) = ApplyPrimaryVariableUpdate(xs, dx, NumActiveGrids, NumWells);
                    }

                    // Cut dt on Newton non-convergence
                    if (!isConverged)
                    {
                        dtCarry = 0.5 * dt; tryingnewtimestep = true;
                        if (dtCarry < 1e-12) throw new Exception("Time step is too small");
                        Console.WriteLine("""
                            ============================
                            Rejected (Non-Convergence)
                            ============================
                            """);
                        continue;
                    }

                    // Validate maximum per-step state changes
                    var (maxDP, maxDS, i_p, i_s) = CheckSolution(xs);
                    if (maxDP > 500 || maxDS > 0.2)
                    {
                        int j = maxDP > 500 ? i_p : i_s;
                        double limit = maxDP > 500 ? 500 : 0.2;
                        double change = maxDP > 500 ? maxDP : maxDS;
                        Console.ForegroundColor = ConsoleColor.Yellow;
                        Console.WriteLine($"\n[Warning] Time Step Rejected: State Change Limit Exceeded.");
                        Console.ResetColor();
                        Console.WriteLine($"    ↳ Location       : Cell Index {j}");
                        Console.WriteLine($"    ↳ Violation      : {(maxDP > 500 ? "ΔP" : "ΔS")} = {change:F4} (Max Allowed: {limit:F4})");
                        Console.WriteLine($"    ↳ Current State  : P = {Grids[j].Po_n:F1} psi | Sw = {Grids[j].Sw_n:F4}");
                        Console.WriteLine($"    ↳ Attempted State: P = {Grids[j].Po_np1:F1} psi | Sw = {Grids[j].Sw_np1:F4}");
                        Console.WriteLine($"    ↳ Solver Action  : Rolling back to t = {t:F3} days.");
                        Console.WriteLine($"    ↳ Time-Step Chop : Reducing Δt from {dt:F2} to {0.5 * dt:F2} days (Factor: 0.5).\n");
                        dtCarry = 0.5 * dt; tryingnewtimestep = true;
                        continue;
                    }

                    // Validate Well BHP constraints (reuses existing PARDISO factorization on constraint switch)
                    foreach (var well in Wells)
                    {
                        if (!well.isActive || well.isStopped) continue;

                        if (well.WellType == WellType.PRODUCER &&
                            well.CurrentProducerControl is { } prod &&
                            well.ConstraintType != ConstraintType.MIN_BHP &&
                            prod.MinBhpLimit > 0.0 &&
                            well.Pwf_np1 < (prod.MinBhpLimit - 1e-2))
                        {
                            well.ConstraintType = ConstraintType.MIN_BHP;
                            Console.WriteLine($"""
                            ===============================================================
                            Rejected (Minimum Pressure Violated: {well.Pwf_np1:F2} < {prod.MinBhpLimit:F2} psi) @ {well.Name}
                            ===============================================================
                            """);
                            well.Control_is_Violated = true;
                        }

                        if (well.WellType == WellType.INJECTOR &&
                            well.CurrentInjectorControl is { } inj &&
                            well.ConstraintType != ConstraintType.MAX_BHP &&
                            inj.MaxBhpLimit > 0.0 &&
                            well.Pwf_np1 > (inj.MaxBhpLimit + 1e-2))   // was well.Pwf.Value: Pwf is never set, so always 0
                        {
                            well.ConstraintType = ConstraintType.MAX_BHP;
                            Console.WriteLine($"""
                                ===============================================================
                                Rejected (Maximum Pressure Violated: {well.Pwf_np1:F2} > {inj.MaxBhpLimit:F2} psi) @ {well.Name}
                                ===============================================================
                                """);
                            well.Control_is_Violated = true;
                        }

                        // Back to the deck's rate control when, at the BHP limit, the well would flow more than its target
                        // (ECLIPSE does the same). Otherwise a well that once hit its limit stays on BHP and overshoots
                        // its rate target for the rest of that control period.
                        if (well.WellType == WellType.PRODUCER &&
                            well.CurrentProducerControl is { } p &&
                            well.ConstraintType == ConstraintType.MIN_BHP && well.DeckConstraint != ConstraintType.MIN_BHP)
                        {
                            var (rate, target) = well.DeckConstraint switch
                            {
                                ConstraintType.ORATE => (well.Qo_np1, p.OilRateTarget),
                                ConstraintType.WRATE => (well.Qw_np1, p.WaterRateTarget),
                                ConstraintType.LRATE => (well.Ql_np1, p.LiquidRateTarget),
                                _ => (0.0, 0.0)
                            };
                            if (target > 0.0 && rate > target * (1 + 1e-4))
                            {
                                well.ConstraintType = well.DeckConstraint;
                                Console.WriteLine($"[Info] {well.Name}: rate {rate:F2} above target {target:F2} at the BHP limit; back to rate control.");
                                well.Control_is_Violated = true;
                            }
                        }

                        if (well.WellType == WellType.INJECTOR &&
                            well.CurrentInjectorControl is { } i2 &&
                            well.ConstraintType == ConstraintType.MAX_BHP && well.DeckConstraint == ConstraintType.RATE &&
                            Abs(well.Qw_np1) > Abs(i2.TargetRate) * (1 + 1e-4))
                        {
                            well.ConstraintType = ConstraintType.RATE;
                            Console.WriteLine($"[Info] {well.Name}: rate {Abs(well.Qw_np1):F2} above target {Abs(i2.TargetRate):F2} at the BHP limit; back to rate control.");
                            well.Control_is_Violated = true;
                        }
                    }
                    if (Wells.Any(well => well.Control_is_Violated))
                    {
                        foreach(var well in Wells)
                            well.Control_is_Violated = false;
                        // Solve the same step again with the new controls; at most 4 times, so two controls that keep
                        // undoing each other cannot loop for ever (the step is then accepted as solved).
                        if (++controlSwitches <= 4)
                        {
                            tryingnewtimestep = false;
                            continue;
                        }
                    }

                    // Commit verified state n+1 -> n
                    t = tnp1;
                    Time.Add(t);
                    // The next step loads the deck's controls and completions again. Without this, a BHP switch above
                    // left the flag false for good, and later WCONPROD/WCONINJE/COMPDAT records were never loaded.
                    // (A switch to BHP still carries over: UpdateControl resets the mode only for a new record.)
                    tryingnewtimestep = true;
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
                        well.Pwf_n = well.Pwf_np1; well.Qo_n = well.Qo_np1;
                        well.Ql_n = well.Ql_np1; well.Qw_n = well.Qw_np1;
                        well.SaveResult(pwf: well.Pwf_n, qo: well.Qo_n, qw: well.Qw_n);
                    }

                    isComplete = Abs(t - L) < 1e-13;

                    // Adaptive time-stepping update preserving pre-clamped dtCarry across milestones
                    double nextDt = dt;
                    if (iter < 4) nextDt *= 1.25;
                    if (iter > 8) nextDt *= 0.7;

                    dtCarry = hitMilestone ? Clamp(Max(dtCarry, nextDt), 1e-12, 100.0) : Clamp(nextDt, 1e-12, 100.0);
                }
            }
        }

        public static void ExportCsrToTriplets(int[] rowPtr, int[] colIndex, double[] values = null,
            string filePath = "triplets.txt", bool oneBasedIndexing = true)
        {
            int numRows = rowPtr.Length - 1;
            int indexOffset = oneBasedIndexing ? 1 : 0;

            using (StreamWriter writer = new StreamWriter(filePath))
            {
                for (int row = 0; row < numRows; row++)
                {
                    for (int idx = rowPtr[row]; idx < rowPtr[row + 1]; idx++)
                    {
                        int r = row + indexOffset;
                        int c = colIndex[idx] + indexOffset;

                        // Use the actual value if provided; otherwise default to 1.0
                        double val = (values != null) ? values[idx] : 1.0;

                        // InvariantCulture ensures '.' is used as the decimal separator
                        writer.WriteLine($"{r} {c} {val.ToString(CultureInfo.InvariantCulture)}");
                    }
                }
            }
        }

        public void ExportResults(string outputDir)
        {
            //VTUExporter.WriteVtuBinaryOilWater(outputDir, "OilWater.vtu", Grids);
            ExportWellVTK.ExportWells(outputDir, Grids, Wells);

            var animatorExporter = new SimulationAnimatorExporter<Vector4d>(outputDir, Grids);
            for (int n = 0; n < Time.Count; n++)
                animatorExporter.WriteTimeStep(Time, Grids, n);


            foreach (var well in Wells)
            {
                List<ColVec> vectors = [Time, well.PWF, well.QO, well.QW];
                string filename = Path.Combine(outputDir, $"{well.Name}.txt");
                SepalSolver.Math.WriteMatrix(vectors, filename);

                for (int i = 0; i < Time.Count; i++)
                {

                }
            }
        }
    }
}
