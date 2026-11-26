using EclipseDeckReader;
using Microsoft.CodeAnalysis;
using SepalSolver;
using System.Globalization;
using static ReservoirSimulator.Math;
using static ReservoirSimulator.ReservoirJacobian;

namespace ReservoirSimulator
{
    public class DryGasLiveOilReservoir : IReservoir
    {
        public int funcall;
        public bool UseUpstreamWeightingForInjection
        {
            get => useUpstreamweighting;
            set => useUpstreamweighting = value;
        }
        bool useUpstreamweighting = false;
        public List<double> Time, SweepEff;
        readonly Func<ADiff<Vector6d>, ADiff<Vector6d>> Er;
        (double maxDP, double maxDSw, double maxDSg, int i_p, int i_w, int i_g) CheckSolution(double[] xs)
        {
            int i_p = 0, i_w = 0, i_g = 0;
            double dp, maxDp = 0, ds, maxDSw = 0, maxDSg = 0;
            for (int i = 0; i < Ngrids; i++)
            {
                dp = Grids[i].Po_np1 - Grids[i].Po_n;
                if (Abs(dp) > maxDp) { maxDp = Abs(dp); i_p = i; }
                ds = Grids[i].Sw_np1 - Grids[i].Sw_n;
                if (Abs(ds) > maxDSw) { maxDSw = Abs(ds); i_w = i; }
                ds = Grids[i].Sg_np1 - Grids[i].Sg_n;
                if (Abs(ds) > maxDSg) { maxDSg = Abs(ds); i_g = i; }
            }
            return (maxDp, maxDSw, maxDSg, i_p, i_w, i_g);
        }

        readonly int Nx, Ny, Nz, NxNy, Ngrids, NumActiveGrids, NumWells;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, zTop;
        double[] Coord, Zcorn;
        readonly double Prefr, cr, γo0, γw0, γg0;
        public Grid<Vector6d>[] Grids;
        public List<Grid<Vector6d>> ActiveGrids;
        public List<Well<Vector6d>> Wells;
        Aquifer Aquifer;
        Dictionary<(int, int), EQUILPressureDepth> EquilPVTCombinations = [];
        EclipseDataDeck DATA;
        PVTWAnalysis[] PVTWs;
        PVTOAnalysis[] PVTOs;
        PVDGAnalysis[] PVDGs;
        SWOFAnalysis[] SWOFs;
        SGOFAnalysis[] SGOFs;
        RSVDAnalysis[] RSVDs;

        public List<Well<Vector6d>> WELLS => Wells;

        public Grid<Vector6d>[] GRIDS => Grids;

        public DryGasLiveOilReservoir(EclipseDataDeck Data)
        {
            DATA = Data;
            // RUNSPEC
            Nx = DATA.NX;
            Ny = DATA.NY;
            Nz = DATA.NZ;
            NxNy = Nx * Ny;
            Ngrids = Nx * Ny * Nz;

            Units incomingUnit = DATA.UnitSystem.ToUpper() switch
            {
                "FIELD" => Units.FIELD,
                "METRIC" => Units.METRIC,
                "LAB" => Units.LAB,
                "PVT-M" => Units.PVT_M,
                _ => throw new InvalidOperationException($"Unsupported unit system: {DATA.UnitSystem}")
            };

            // GRIDS
            if (DATA.DxArray.Count > 0 && DATA.DyArray.Count > 0 && DATA.DzArray.Count > 0)
            {
                Dx = [.. Utility.Convert(DATA.DxArray, "Length", incomingUnit, Units.FIELD)];
                Dy = [.. Utility.Convert(DATA.DyArray, "Length", incomingUnit, Units.FIELD)];
                Dz = [.. Utility.Convert(DATA.DzArray, "Depth", incomingUnit, Units.FIELD)];
                zTop = DATA.ZtopArray.Count > 0 ? [.. Utility.Convert(DATA.ZtopArray, "Depth", incomingUnit, Units.FIELD)] : new double[NxNy];
            }
            else
            {
                Coord = [.. Utility.Convert(DATA.CoordArray, "Length", incomingUnit, Units.FIELD)];
                Zcorn = [.. Utility.Convert(DATA.ZcornArray, "Depth", incomingUnit, Units.FIELD)];
            }
            Φ = [.. DATA.PorosityArray];
            Kx = [.. DATA.PermXArray]; Ky = [.. DATA.PermYArray]; Kz = [.. DATA.PermZArray];

            // NTG and MULTPV default to 1 for every cell (set at DIMENS); a wrong-sized array is ignored rather than misread
            double[]? ntg = DATA.NtgArray.Count == Ngrids ? DATA.NtgArray.ToArray() : null;
            double[]? multpv = DATA.MultpvArray.Count == Ngrids ? DATA.MultpvArray.ToArray() : null;
            Grids = Zcorn is not null ? Utility.ProcessGrid<Vector6d>(Nx, Ny, Nz, Coord, Zcorn, Φ, Kx, Ky, Kz, ntg, multpv) :
                Utility.ProcessGrid<Vector6d>(Nx, Ny, Nz, Dx, Dy, Dz, zTop, Φ, Kx, Ky, Kz, ntg, multpv);

            // Process PVTW, PVTO, and PVDG tables.
            PVTWs = [.. DATA.RawPvtwTables.Select(pvtw => new PVTWAnalysis(pvtw, incomingUnit))];
            PVTOs = [.. DATA.RawPvtoTables.Select(pvto => new PVTOAnalysis(pvto, incomingUnit))];
            PVDGs = [.. DATA.RawPvdgTables.Select(pvdg => new PVDGAnalysis(pvdg, incomingUnit))];
            for (int i = 0; i < PVTWs.Length; i++)
            {
                PVTWs[i].FitProperties();
                PVTOs[i].FitProperties();
                PVDGs[i].FitProperties();
            }

            // Process SWOF and SGOF tables.
            SWOFs = [.. DATA.RawSwofTables.Select(swof => new SWOFAnalysis(swof, incomingUnit))];
            SGOFs = [.. DATA.RawSgofTables.Select(sgof => new SGOFAnalysis(sgof, incomingUnit))];
            for (int i = 0; i < SGOFs.Length; i++)
            {
                SWOFs[i].FitProperties();
                SGOFs[i].FitProperties(SWOFs[i].Result.Swi);
            }
             // Process RSVD tables
            RSVDs = [.. DATA.RawRsvdTables.Select(rsvd => new RSVDAnalysis(rsvd, incomingUnit))];
            for (int i = 0; i < RSVDs.Length; i++)
            {
                RSVDs[i].FitProperties();
            }

            // ROCK
            cr = DATA.RockCompressibility;
            Prefr = Utility.Convert(DATA.RockReferencePressure, "Pressure", incomingUnit, Units.FIELD);
            Er = P => Exp(cr * (P - Prefr));

            const double SCF_PER_MSCF = 1000.0;
            const double FT3_PER_BBL = 5.614583;

            // DENSITIES -> Gradients [psi/ft * (rb/stb) for oil & water; psi/ft * (rb/Mscf) for gas]
            γo0 = Utility.Convert(DATA.SurfaceOilDensity, "Density", incomingUnit, Units.FIELD) / 144.0;
            γw0 = Utility.Convert(DATA.SurfaceWaterDensity, "Density", incomingUnit, Units.FIELD) / 144.0;
            γg0 = (Utility.Convert(DATA.SurfaceGasDensity, "Density", incomingUnit, Units.FIELD) / 144.0) * (SCF_PER_MSCF / FT3_PER_BBL);

            // MINPV: cells holding less pore volume than the deck's minimum are made inactive, as in ECLIPSE.
            // A copy of ACTNUM is changed, so the deck itself keeps what the user wrote.
            List<int> actnum = [.. DATA.ActnumArray];
            if (DATA.MinPoreVolume > 0)
            {
                double minPv = Utility.Convert(DATA.MinPoreVolume, "Reservoir Volume", incomingUnit, Units.FIELD);
                int removed = 0;
                for (int m = 0; m < Ngrids; m++)
                    if (actnum[m] > 0 && Grids[m].PoreVolume < minPv) { actnum[m] = 0; removed++; }
                Console.WriteLine($"      [Grid] MINPV {DATA.MinPoreVolume}: {removed} cell(s) with less pore volume made inactive.");
            }

            // REGION
            NumActiveGrids = Utility.SetActiveGrids(Grids, actnum, Er);
            ActiveGrids = [.. Grids.Where(g => g.isActive)];
            foreach (var grid in ActiveGrids)
            {
                int satnum = DATA.SatnumArray[grid.GlobalIndex];
                grid.SWOF = SWOFs[satnum];
                grid.SGOF = SGOFs[satnum];
                int pvtnum = DATA.PvtnumArray[grid.GlobalIndex];
                grid.PVTW = PVTWs[pvtnum];
                grid.PVDG = PVDGs[pvtnum];
                grid.PVTO = PVTOs[pvtnum];
                int equinum = DATA.EquinumArray[grid.GlobalIndex];
                grid.RSVD = RSVDs[equinum];
            }

            // TRANSMISSIBILITY
            Utility.BuildStructuredConnections(Grids, Nx, Ny, Nz);
            Utility.ApplyTransmissibilityMultipliers(Grids, Nx, Ny, Nz, DATA);


            // 1. Build EQUIL Pressure-Depth Tables for every (EQUILNUM, PVTNUM) pair
            EQUILAnalysis[] EQUILs = [.. DATA.RawEquilTables.Select(equil => new EQUILAnalysis(equil, incomingUnit))];
            double MaxZ = Grids.Max(grid => grid.Center.Z);
            double MinZ = Grids.Min(grid => grid.Center.Z);

            for (int e = 0; e < EQUILs.Length; e++)
            {
                for (int p = 0; p < PVTOs.Length; p++)
                {
                    var equil = EQUILs[e];
                    var pvto = PVTOs[p];
                    var pvtw = PVTWs[p];
                    var pvdg = PVDGs[p];
                    var rsvd = e < RSVDs.Length ? RSVDs[e] : null;

                    Func<double, double> γw = P => γw0 / pvtw.Result.EvaluateBw(P),
                                         γg = P => γg0 / pvdg.Result.EvaluateBg(P);

                    double RsSat(double P)
                    {
                        double RsBase = pvto.Result.RsBase;
                        double RsExponent = pvto.Result.RsExponent;
                        if (P <= 0.0 || RsBase <= 0.0)
                            return 0.0;
                        return RsBase * Pow(P, RsExponent);
                    }

                    double BoSat(double P)
                    {
                        double RsBase = pvto.Result.RsBase;
                        double BoBase = pvto.Result.BoBase;
                        double BoExponent = pvto.Result.BoExponent;
                        if (P <= 0.0 || RsBase <= 0.0)
                            return 0.0;
                        double rsSat = RsSat(P);
                        return BoBase * Exp(BoExponent * rsSat);
                    }

                    double Pb(double RsVal)
                    {
                        double RsBase = pvto.Result.RsBase;
                        double RsExponent = pvto.Result.RsExponent;
                        if (RsVal <= 0.0 || RsBase <= 0.0)
                            return 0.0;
                        return Pow(RsVal / RsBase, 1.0 / RsExponent);
                    }

                    (double Bo, double RsEff) EvaluateOilState(double P, double RsVal)
                    {
                        double RsBase = pvto.Result.RsBase;
                        double BoBase = pvto.Result.BoBase;
                        double BoExponent = pvto.Result.BoExponent;
                        double UndersatCompressibility = pvto.Result.UndersatCompressibility;

                        // Dead oil (PVDO)
                        if (RsBase <= 0.0)
                        {
                            double boDead = BoBase * Exp(-UndersatCompressibility * Max(0.0, P - 14.7));
                            return (boDead, 0.0);
                        }

                        // Live oil (PVTO): clamp Rs to RsSat(P) if P <= Pb(Rs)
                        double pb = Pb(RsVal);
                        if (P <= pb)
                        {
                            double rsSat = RsSat(P), boSat = BoSat(P);
                            return (boSat, rsSat);
                        }
                        else
                        {
                            double boAtPb = BoBase * Exp(BoExponent * RsVal);
                            double boUndersat = boAtPb * Exp(-UndersatCompressibility * (P - pb));
                            return (boUndersat, RsVal);
                        }
                    }

                    double γopr(double P, double RsVal)
                    {
                        var (bo, rsEff) = EvaluateOilState(P, RsVal);
                        return (γo0 + 178.1076 * rsEff * γg0) / bo;
                    }

                    double Rszp(double Z, double P)
                    {
                        double rssat = RsSat(P);
                        if (rssat <= 0.0 || rsvd is null) return 0.0;
                        double rs = rsvd.Result.EvaluateRs(Z);
                        return Min(rs, rssat);
                    }


                    // Oil leg gradient (existing): dPo/dz = c * γopr(p, Rszp(z, p))
                    double OilGrad(double z, double p) => γopr(p, Rszp(z, p));

                    // Gas cap gradient (updated for VAPOIL): dPg/dz = c * γgpr(p, Rvzp(z, p))
                    double GasGrad(double z, double p) => γg(p);

                    // Water leg gradient (existing): dPw/dz = c * γwpr(p)
                    double WaterGrad(double z, double p) => γw(p);

                    var equildata = new EQUILPressureDepth(equil);
                    double zSpan = Max(1.0, MaxZ - MinZ);

                    int StepsFor(double z1, double z2, int scale = 1) =>
                        Max(10, (int)(Abs(z1 - z2) / zSpan * Nz * scale));

                    // 1. Resolve Po, Pw, Pg at ZGOC and ZWOC starting from (ZDatum, PDatum)
                    //    Note: IntegrateLegRK4 returns pairs sorted by ascending Z, so upward legs
                    //    (zEnd < zStart) have their target pressure at First().P, while downward legs have Last().P.
                    if (equildata.ZDatum <= equildata.ZGOC)
                    {
                        // Datum in gas cap -> integrate gas down to ZGOC, then oil down to ZWOC
                        var gas2GOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZGOC,
                            StepsFor(equildata.ZDatum, equildata.ZGOC, 10), GasGrad);
                        equildata.PgGOC = gas2GOC.Last().P;
                        equildata.PoGOC = equildata.PgGOC - equildata.PcGOC;

                        var oilGOC2WOC = Utility.IntegrateLegRK4(equildata.ZGOC, equildata.PoGOC, equildata.ZWOC,
                            StepsFor(equildata.ZGOC, equildata.ZWOC, 10), OilGrad);
                        equildata.PoWOC = oilGOC2WOC.Last().P;
                        equildata.PwWOC = equildata.PoWOC - equildata.PcWOC;
                    }
                    else if (equildata.ZDatum < equildata.ZWOC)
                    {
                        // Datum in oil zone -> integrate oil UP to ZGOC (First().P) and DOWN to ZWOC (Last().P)
                        var oil2GOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZGOC,
                            StepsFor(equildata.ZDatum, equildata.ZGOC, 10), OilGrad);
                        equildata.PoGOC = oil2GOC.First().P;
                        equildata.PgGOC = equildata.PoGOC + equildata.PcGOC;

                        var oil2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC,
                            StepsFor(equildata.ZDatum, equildata.ZWOC, 10), OilGrad);
                        equildata.PoWOC = oil2WOC.Last().P;
                        equildata.PwWOC = equildata.PoWOC - equildata.PcWOC;
                    }
                    else
                    {
                        // Datum in water zone -> integrate water UP to ZWOC (First().P), then oil UP to ZGOC (First().P)
                        var water2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC,
                            StepsFor(equildata.ZDatum, equildata.ZWOC, 10), WaterGrad);
                        equildata.PwWOC = water2WOC.First().P;
                        equildata.PoWOC = equildata.PwWOC + equildata.PcWOC;

                        var oilWOC2GOC = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PoWOC, equildata.ZGOC,
                            StepsFor(equildata.ZWOC, equildata.ZGOC, 10), OilGrad);
                        equildata.PoGOC = oilWOC2GOC.First().P;
                        equildata.PgGOC = equildata.PoGOC + equildata.PcGOC;
                    }

                    // 2. Gas Cap Leg: Z from MinZ down to ZGOC (only if MinZ < ZGOC)
                    if (MinZ < equildata.ZGOC)
                    {
                        double pwAtGOC = equildata.PwWOC;
                        if (equildata.ZGOC < equildata.ZWOC)
                        {
                            // Upward integration from ZWOC -> ZGOC: target pressure is First().P
                            var watGOC = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, equildata.ZGOC,
                                StepsFor(equildata.ZGOC, equildata.ZWOC, 10), WaterGrad);
                            pwAtGOC = watGOC.First().P;
                        }

                        int nGas = StepsFor(MinZ, equildata.ZGOC);
                        // Already returned in ascending Z order [MinZ .. ZGOC] — do NOT reverse
                        var gasUp = Utility.IntegrateLegRK4(equildata.ZGOC, equildata.PgGOC, MinZ, nGas, GasGrad);
                        var oilUp = Utility.IntegrateLegRK4(equildata.ZGOC, equildata.PoGOC, MinZ, nGas, OilGrad);
                        var watUp = Utility.IntegrateLegRK4(equildata.ZGOC, pwAtGOC, MinZ, nGas, WaterGrad);

                        equildata.GasZonePressureDepthPairs.AddRange(gasUp);
                        equildata.OilZonePressureDepthPairs.AddRange(oilUp);
                        equildata.WaterZonePressureDepthPairs.AddRange(watUp);
                    }

                    // 3. Oil Leg: Z from zOilTop down to zOilBot (only if zOilTop < zOilBot)
                    double zOilTop = Max(MinZ, equildata.ZGOC);
                    double zOilBot = Min(MaxZ, equildata.ZWOC);
                    if (zOilTop < zOilBot)
                    {
                        double poAtOilTop = zOilTop > equildata.ZGOC
                            ? Utility.IntegrateLegRK4(equildata.ZGOC, equildata.PoGOC, zOilTop,
                                StepsFor(equildata.ZGOC, zOilTop, 10), OilGrad).Last().P
                            : equildata.PoGOC;

                        // Upward integration from ZWOC -> zOilBot: target pressure is First().P
                        double pwAtOilBot = zOilBot < equildata.ZWOC
                            ? Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, zOilBot,
                                StepsFor(zOilBot, equildata.ZWOC, 10), WaterGrad).First().P
                            : equildata.PwWOC;

                        int nOil = StepsFor(zOilTop, zOilBot);
                        var oilMid = Utility.IntegrateLegRK4(zOilTop, poAtOilTop, zOilBot, nOil, OilGrad);
                        // Already returned in ascending Z order [zOilTop .. zOilBot] — do NOT reverse
                        var watMid = Utility.IntegrateLegRK4(zOilBot, pwAtOilBot, zOilTop, nOil, WaterGrad);

                        // Below ZGOC there is no free gas (Sg = 0): Pg tracks Po + PcGOC
                        var gasMid = oilMid.Select(pz => new PressureDepthPair(pz.Z, pz.P + equildata.PcGOC));

                        int skip = equildata.OilZonePressureDepthPairs.Count > 0 ? 1 : 0;
                        equildata.OilZonePressureDepthPairs.AddRange(oilMid.Skip(skip));
                        equildata.WaterZonePressureDepthPairs.AddRange(watMid.Skip(skip));
                        equildata.GasZonePressureDepthPairs.AddRange(gasMid.Skip(skip));
                    }

                    // 4. Water Leg: Z from zWatTop down to MaxZ (only if zWatTop < MaxZ)
                    double zWatTop = Max(MinZ, equildata.ZWOC);
                    if (zWatTop < MaxZ)
                    {
                        double pwAtWatTop = zWatTop > equildata.ZWOC
                            ? Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, zWatTop,
                                StepsFor(equildata.ZWOC, zWatTop, 10), WaterGrad).Last().P
                            : equildata.PwWOC;

                        int nDown = StepsFor(zWatTop, MaxZ);
                        var watDown = Utility.IntegrateLegRK4(zWatTop, pwAtWatTop, MaxZ, nDown, WaterGrad);

                        // Below ZWOC mobile oil/gas vanish (Sw = 1 - Sor): Po = Pw + PcWOC, Pg = Po + PcGOC
                        var oilDown = watDown.Select(pz => new PressureDepthPair(pz.Z, pz.P + equildata.PcWOC));
                        var gasDown = oilDown.Select(pz => new PressureDepthPair(pz.Z, pz.P + equildata.PcGOC));

                        int skip = equildata.WaterZonePressureDepthPairs.Count > 0 ? 1 : 0;
                        equildata.WaterZonePressureDepthPairs.AddRange(watDown.Skip(skip));
                        equildata.OilZonePressureDepthPairs.AddRange(oilDown.Skip(skip));
                        equildata.GasZonePressureDepthPairs.AddRange(gasDown.Skip(skip));
                    }

                    EquilPVTCombinations[(e, p)] = equildata;
                }
            }

            // WELLS
            Wells = [.. DATA.Wells.Select(well => new Well<Vector6d>(well.Value, Nx, Ny, Nz, Grids, incomingUnit))];
            foreach (var well in Wells.Where(w => w.PerfNatIndex.Count == 0))
                Console.WriteLine($"      [Warning] Well '{well.Name}' has no completion in an active cell; it is left out of the simulation.");
            Wells.RemoveAll(w => w.PerfNatIndex.Count == 0);
            NumWells = Wells.Count;
        }

        public void Initialize()
        {
            funcall = 0;
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;

                int s = DATA.SatnumArray[grid.GlobalIndex];
                int e = DATA.EquinumArray[grid.GlobalIndex];
                int p = DATA.PvtnumArray[grid.GlobalIndex];

                var swoftable = SWOFs[s].swofTable;
                var sgoftable = SGOFs[s].sgofTable;

                var equildata = EquilPVTCombinations[(e, p)];
                var oiltable = equildata.OilZonePressureDepthPairs;
                var watertable = equildata.WaterZonePressureDepthPairs;
                var gastable = equildata.GasZonePressureDepthPairs;

                double z = grid.Center.Z;
                double poEquil = Utility.Interpolate(oiltable, z, pair => pair.Z, pair => pair.P);
                double pwEquil = Utility.Interpolate(watertable, z, pair => pair.Z, pair => pair.P);
                double pgEquil = Utility.Interpolate(gastable, z, pair => pair.Z, pair => pair.P);

                double pcow = poEquil - pwEquil;
                double pcog = pgEquil - poEquil; // Pg >= Po

                // 1. Invert Pcow(Sw) for initial water saturation Sw_n
                // Note: SWOF is ordered by ascending Sw (descending Pc), so reverse for ascending Pc interpolation
                double swMin = SWOFs[s].Result.SwMin, swMax = SWOFs[s].Result.SwMax;
                double pcowMin = SWOFs[s].Result.PcMin, pcowMax = SWOFs[s].Result.PcMax;

                if (z >= equildata.ZWOC) grid.Sw_n = swMax;
                else if (pcowMax <= pcowMin || pcow >= pcowMax) grid.Sw_n = swMin;
                else if (pcow <= pcowMin) grid.Sw_n = swMax;
                else
                    grid.Sw_n = Utility.Interpolate(swoftable.AsEnumerable().Reverse().ToList(), pcow, swof => swof.Pc, swof => swof.Sw);

                // 2. Invert Pcog(Sg) for initial gas saturation Sg_n (only above ZGOC when SGOF is present)
                if (sgoftable is null || sgoftable.Count == 0 || z >= equildata.ZGOC)
                {
                    grid.Sg_n = 0.0;
                }
                else
                {
                    double sgMin = SGOFs[s].Result.SgMin, sgMax = SGOFs[s].Result.SgMax;
                    double pcogMin = SGOFs[s].Result.PcMin, pcogMax = SGOFs[s].Result.PcMax;

                    if (pcogMax <= pcogMin || pcog >= pcogMax)
                        grid.Sg_n = sgMax;
                    else if (pcog <= pcogMin)
                        grid.Sg_n = sgMin;
                    else
                        grid.Sg_n = Utility.Interpolate(sgoftable, pcog, sgof => sgof.Pcg, sgof => sgof.Sg);

                    // Clamp so Sw + Sg never exceeds 1.0 in overlapping transition zones
                    if (grid.Sw_n + grid.Sg_n > 1.0)
                        grid.Sg_n = Max(0.0, 1.0 - grid.Sw_n);
                }

                grid.So_n = Max(0.0, 1.0 - grid.Sw_n - grid.Sg_n);

                // 3. Enforce exact algebraic consistency with Residual:
                //    Pw = Po - Pcow(Sw),  Pg = Po + Pcog(Sg)
                double pcowActual = grid.Pcow(grid.Sw_n).Value;
                double pcogActual = grid.Pcog(grid.Sg_n).Value;

                if (z <= equildata.ZGOC)
                {
                    // Gas cap: gas is the continuous reference phase
                    grid.Pg_n = pgEquil;
                    grid.Po_n = pgEquil - pcogActual;
                    grid.Pw_n = grid.Po_n - pcowActual;
                }
                else if (z >= equildata.ZWOC)
                {
                    // Aquifer: water is the continuous reference phase
                    grid.Pw_n = pwEquil;
                    grid.Po_n = pwEquil + pcowActual;
                    grid.Pg_n = grid.Po_n + pcogActual;
                }
                else
                {
                    // Oil leg: oil is the continuous reference phase
                    grid.Po_n = poEquil;
                    grid.Pw_n = poEquil - pcowActual;
                    grid.Pg_n = poEquil + pcogActual;
                }

                // 4. Set Rs_n and phase-state flag consistently with finalized Po_n and Sg_n
                double rsSat = grid.RsSat(grid.Po_n).Value;
                if (grid.Sg_n > 0.0)
                {
                    grid.Rs_n = rsSat;
                    grid.OilIsSaturated = true;
                }
                else
                {
                    grid.Rs_n = rsSat > 0.0 ? Min(grid.Rsvd(z), rsSat) : 0.0;
                    grid.OilIsSaturated = false;
                }

                // 5. Initial pore-volume-scaled accumulation terms at t = 0
                var meanP = grid.So_n * grid.Po_n + grid.Sw_n * grid.Pw_n + grid.Sg_n * grid.Pg_n;
                var Er_i = grid.Er(meanP);
                var Bo_i = grid.Bopr(grid.Po_n, grid.Rs_n);
                var Bw_i = grid.Bw(grid.Pw_n);
                var Bg_i = grid.Bg(grid.Pg_n);

                grid.ErSo_Bo_n = (Er_i * grid.So_n / Bo_i).Value;
                grid.ErSw_Bw_n = (Er_i * grid.Sw_n / Bw_i).Value;
                grid.ErSg_Bg_n = (Er_i * grid.Sg_n / Bg_i + Er_i * grid.Rs_n * grid.So_n / Bo_i).Value;

                grid.SaveResult(po: grid.Po_n, pw: grid.Pw_n, pg: grid.Pg_n,
                                so: grid.So_n, sw: grid.Sw_n, sg: grid.Sg_n);
            }

            foreach (var well in Wells)
            {
                well.Pwf_n = well.PerfNatIndex.Select(m => Grids[m].Po_n).Average();
                well.SaveResult(pwf: well.Pwf_n, qo: well.Qo_n, qw: well.Qw_n, qg: well.Qg_n);
            }
        }
        void InitialGuess(double[] xs, double dt, bool isRetry = false, bool wellEventOccurred = false)
        {
            int indx = 0, n = Time.Count - 1;
            bool canExtrapolate = n >= 2 && !isRetry && !wellEventOccurred;
            double omega = canExtrapolate ? Min(1.0, dt / (Time[n] - Time[n - 1])) : 0.0;

            // 1. Reservoir Grid Unknowns (Po, Sw, Sg/Rs)
            foreach (var mgrid in ActiveGrids)
            {
                double poGuess = mgrid.Po_n;
                if (canExtrapolate)
                    poGuess += omega * (mgrid.Po_n - mgrid.PO[n - 1]);

                double po0 = Max(14.7, poGuess);
                xs[indx++] = po0;                                                   // Po index
                xs[indx++] = mgrid.Sw_n;                                            // Sw index

                if (mgrid.Sg_n > 0.0)
                {
                    mgrid.OilIsSaturated = true;
                    xs[indx++] = mgrid.Sg_n;                                        // Sg index
                }
                else
                {
                    mgrid.OilIsSaturated = false;
                    double rsSat0 = mgrid.RsSat(po0).Value;
                    xs[indx++] = rsSat0 > 0.0 ? Min(mgrid.Rs_n, rsSat0) : 0.0;     // Rs index (clamped to bubble point)
                }
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

        public (double maxPd, double maxSd, double maxRsd, bool phaseSwitched) ApplyPrimaryVariableUpdate(
            double[] x, double[] dx, double dt, int numActiveGrids, int numWells,
            double dPmax = 500.0,    // max allowable pressure change per Newton step [psi]
            double dSmax = 0.10,     // max allowable saturation change per Newton step [-]
            double dRsMax = 0.20)    // max allowable Rs change per Newton step [Mscf/stb]
        {
            double maxPd = 0.0;
            double maxSd = 0.0;
            double maxRsd = 0.0;
            bool phaseSwitched = false;

            // 1. Grid cell updates
            for (int i = 0; i < numActiveGrids; i++)
            {
                var grid = ActiveGrids[i];
                int idx_p = 3 * i + 0, idx_w = 3 * i + 1, idx_g = 3 * i + 2;

                double p_old = x[idx_p];
                double w_old = x[idx_w];
                double third_old = x[idx_g]; // Sg if OilIsSaturated, Rs otherwise

                double dPo = dx[idx_p];
                double dSw = dx[idx_w];

                // Pressure damping
                if (Abs(dPo) > dPmax) dPo = Sign(dPo) * dPmax;
                double p_new = Max(14.7, p_old - dPo);

                // Saturated Rs at the updated pressure (0.0 for dead oil), capped by DRSDT if active
                double rs_sat = grid.RsSat(p_new).Value;
                double rs_limit = rs_sat; 

                double w_new, third_new;

                if (!grid.OilIsSaturated)
                {
                    // Undersaturated state: variables are (Po, Sw, Rs), Sg = 0
                    double dRs = dx[idx_g];

                    if (Abs(dSw) > dSmax) dSw = Sign(dSw) * dSmax;
                    if (Abs(dRs) > dRsMax) dRs = Sign(dRs) * dRsMax;

                    w_new = Clamp(w_old - dSw, 0.0, 1.0);
                    double rs_trial = third_old - dRs;

                    if (rs_limit > 0.0 && rs_trial >= rs_limit)
                    {
                        // Bubble-point / DRSDT ceiling crossing: Undersaturated -> Saturated
                        grid.OilIsSaturated = true;
                        phaseSwitched = true;

                        double sg_init = Min(1e-4, Max(0.0, 1.0 - w_new));
                        third_new = sg_init;

                        maxSd = Max(maxSd, sg_init);
                        maxRsd = Max(maxRsd, Abs(rs_limit - third_old));
                    }
                    else
                    {
                        // Stay undersaturated (or dead oil where rs_limit == 0)
                        third_new = rs_limit > 0.0 ? Clamp(rs_trial, 0.0, rs_limit) : 0.0;
                        maxRsd = Max(maxRsd, Abs(third_new - third_old));
                    }

                    maxSd = Max(maxSd, Abs(w_new - w_old));
                }
                else
                {
                    // Saturated two-phase state: variables are (Po, Sw, Sg), Rs = EffectiveRsSat(Po, dt)
                    double dSg = dx[idx_g];

                    // Direction-preserving saturation damping
                    double maxCellDs = Max(Abs(dSw), Abs(dSg));
                    if (maxCellDs > dSmax)
                    {
                        double scale = dSmax / maxCellDs;
                        dSw *= scale;
                        dSg *= scale;
                    }

                    w_new = Clamp(w_old - dSw, 0.0, 1.0);
                    double sg_trial = third_old - dSg;

                    if (sg_trial <= 0.0)
                    {
                        // Gas disappearance: Saturated -> Undersaturated
                        grid.OilIsSaturated = false;
                        phaseSwitched = true;

                        double rs_init = rs_limit > 0.0 ? Max(0.0, rs_limit * (1.0 - 1e-6)) : 0.0;
                        third_new = rs_init;

                        maxSd = Max(maxSd, Abs(w_new - w_old));
                        maxSd = Max(maxSd, third_old); // Sg went from third_old -> 0
                        maxRsd = Max(maxRsd, Abs(rs_init - rs_limit));
                    }
                    else
                    {
                        // Stay saturated: enforce 3-phase simplex Sw + Sg <= 1.0
                        double g_new = Clamp(sg_trial, 0.0, 1.0);
                        double sSum = w_new + g_new;
                        if (sSum > 1.0)
                        {
                            double invSum = 1.0 / sSum;
                            w_new *= invSum;
                            g_new *= invSum;
                        }

                        third_new = g_new;
                        maxSd = Max(maxSd, Abs(w_new - w_old));
                        maxSd = Max(maxSd, Abs(g_new - third_old));
                    }
                }

                maxPd = Max(maxPd, Abs(p_new - p_old));

                x[idx_p] = p_new;
                x[idx_w] = w_new;
                x[idx_g] = third_new;
            }

            // 2. Well variables update
            int wellStart = 3 * numActiveGrids;
            for (int w = 0; w < numWells; w++)
            {
                int idxP = wellStart + 2 * w + 0;
                int idxQ = wellStart + 2 * w + 1;
                double dP = dx[idxP], dQ = dx[idxQ];

                if (Abs(dP) > dPmax) dP = Sign(dP) * dPmax;

                x[idxP] = Max(14.7, x[idxP] - dP);
                x[idxQ] -= dQ;
            }

            return (maxPd, maxSd, maxRsd, phaseSwitched);
        }

        double dRsdt = double.PositiveInfinity;
        double dRvdt = double.PositiveInfinity;

        private double GetDrsdtForStep(int stepIdx)
        {
            double activeRate = double.PositiveInfinity;
            foreach (var rec in DATA.DrsdtSchedule)
            {
                if (rec.StepIndex <= stepIdx)
                    activeRate = rec.MaxRate;
                else
                    break;
            }
            return activeRate;
        }

        public void Simulate(double[] ResultTime = null)
        {
            if (ResultTime == null) ResultTime = Utility.GetCumulativeTimeArray(DATA.TimeSteps);
            int cores = Environment.ProcessorCount, chunkSize = Ngrids / cores;

            (int totalDofs, int[] col_index, int[] row_ptr, double[] JacValues) = ReservoirJacobian.Create(3, Grids, Wells);
            ReservoirJacobian.PrecomputeOffsets(3, row_ptr, col_index, Grids, Wells);
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
                        int row_o = 3 * i + 0, row_w = 3 * i + 1, row_g = 3 * i + 2;

                        // Primary variables for cell i
                        var Po_i = new ADiff<Vector6d>(xnp1[3 * i + 0], new(1, 0, 0, 0, 0, 0));
                        var Sw_i = new ADiff<Vector6d>(xnp1[3 * i + 1], new(0, 1, 0, 0, 0, 0));
                        var Sg_i = new ADiff<Vector6d>(xnp1[3 * i + 2], new(0, 0, 1, 0, 0, 0));
                        var Rs_i = mgrid.RsSat(Po_i);
                        if (!mgrid.OilIsSaturated)
                        {
                            Sg_i = 0.0;
                            Rs_i = new ADiff<Vector6d>(xnp1[3 * i + 2], new(0, 0, 1, 0, 0, 0));
                        }

                        var Pw_i = Po_i - mgrid.Pcow(Sw_i);
                        var Pg_i = Po_i + mgrid.Pcog(Sg_i);
                        var So_i = 1.0 - Sw_i - Sg_i;

                        var meanP = So_i * Po_i + Sw_i * Pw_i + Sg_i * Pg_i;
                        var Er_i = mgrid.Er(meanP);
                        var Bo_i = mgrid.Bopr(Po_i, Rs_i);
                        var Bw_i = mgrid.Bw(Pw_i);
                        var Bg_i = mgrid.Bg(Pg_i);
                        var Kro_i = mgrid.Kro(Sw_i, Sg_i);
                        var Krw_i = mgrid.Krw(Sw_i);
                        var Krg_i = mgrid.Krg(Sg_i);
                        var μo_i = mgrid.μopr(Po_i, Rs_i);
                        var μw_i = mgrid.μw(Pw_i);
                        var μg_i = mgrid.μg(Pg_i);
                        var γo_i = mgrid.γopr(Po_i, Rs_i);
                        var γw_i = mgrid.γw(Pw_i);
                        var γg_i = mgrid.γg(Pg_i);
                        var λo_i = Kro_i / (μo_i * Bo_i);
                        var λw_i = Krw_i / (μw_i * Bw_i);
                        var λg_i = Krg_i / (μg_i * Bg_i);

                        mgrid.Po_np1 = Po_i.Value; mgrid.So_np1 = So_i.Value; mgrid.Bo_np1 = Bo_i.Value;
                        mgrid.Pw_np1 = Pw_i.Value; mgrid.Sw_np1 = Sw_i.Value; mgrid.Bw_np1 = Bw_i.Value;
                        mgrid.Pg_np1 = Pg_i.Value; mgrid.Sg_np1 = Sg_i.Value; mgrid.Bg_np1 = Bg_i.Value;
                        mgrid.ErSo_Bo_np1 = Er_i * So_i / Bo_i; mgrid.Rs_np1 = Rs_i.Value;
                        mgrid.ErSw_Bw_np1 = Er_i * Sw_i / Bw_i;
                        mgrid.ErSg_Bg_np1 = Er_i * Sg_i / Bg_i + Er_i * Rs_i * So_i / Bo_i;

                        // Accumulation Residual: R = Vp * (State^{n+1} - State^n)
                        var ResOil = mgrid.PoreVolume * (mgrid.ErSo_Bo_np1 - mgrid.ErSo_Bo_n);
                        var ResWater = mgrid.PoreVolume * (mgrid.ErSw_Bw_np1 - mgrid.ErSw_Bw_n);
                        var ResGas = mgrid.PoreVolume * (mgrid.ErSg_Bg_np1 - mgrid.ErSg_Bg_n);

                        // Direct write to Residual Vector
                        Res[row_o] += ResOil.Value;
                        Res[row_w] += ResWater.Value;
                        Res[row_g] += ResGas.Value;

                        // Direct write to Accumulation Jacobian (Block i, i)
                        var cOff = mgrid.Offset;
                        JacValues[cOff.Diag[0]] += ResOil.Derivative[0];   // row_o, col_p_i
                        JacValues[cOff.Diag[1]] += ResOil.Derivative[1];   // row_o, col_sw_i
                        JacValues[cOff.Diag[2]] += ResOil.Derivative[2];   // row_o, col_sg_i

                        JacValues[cOff.Diag[3]] += ResWater.Derivative[0]; // row_w, col_p_i
                        JacValues[cOff.Diag[4]] += ResWater.Derivative[1]; // row_w, col_sw_i
                        JacValues[cOff.Diag[5]] += ResWater.Derivative[2]; // row_w, col_sg_i

                        JacValues[cOff.Diag[6]] += ResGas.Derivative[0];  // row_g, col_p_i
                        JacValues[cOff.Diag[7]] += ResGas.Derivative[1];  // row_g, col_sw_i
                        JacValues[cOff.Diag[8]] += ResGas.Derivative[2];  // row_g, col_sg_i

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

                                var Po_j = new ADiff<Vector6d>(xnp1[3 * j + 0], new(0, 0, 0, 1, 0, 0));
                                var Sw_j = new ADiff<Vector6d>(xnp1[3 * j + 1], new(0, 0, 0, 0, 1, 0));
                                var Sg_j = new ADiff<Vector6d>(xnp1[3 * j + 2], new(0, 0, 0, 0, 0, 1));
                                var Rs_j = ngrid.RsSat(Po_j);
                                if (!ngrid.OilIsSaturated)
                                {
                                    Sg_j = 0.0;
                                    Rs_j = new ADiff<Vector6d>(xnp1[3 * j + 2], new(0, 0, 0, 0, 0, 1));
                                }

                                var Pw_j = Po_j - ngrid.Pcow(Sw_j);
                                var Pg_j = Po_j + ngrid.Pcog(Sg_j);
                                var So_j = 1.0 - Sw_j - Sg_j;
                                var Bo_j = ngrid.Bopr(Po_j, Rs_j);
                                var Bw_j = ngrid.Bw(Pw_j);
                                var Bg_j = ngrid.Bg(Pg_j);
                                var Kro_j = ngrid.Kro(Sw_j, Sg_j);
                                var Krw_j = ngrid.Krw(Sw_j);
                                var Krg_j = ngrid.Krg(Sg_j);
                                var μo_j = ngrid.μopr(Po_j, Rs_j);
                                var μw_j = ngrid.μw(Pw_j);
                                var μg_j = ngrid.μg(Pg_j);
                                var γo_j = ngrid.γopr(Po_j, Rs_j);
                                var γw_j = ngrid.γw(Pw_j);
                                var γg_j = ngrid.γg(Pg_j);
                                var λo_j = Kro_j / (μo_j * Bo_j);
                                var λw_j = Krw_j / (μw_j * Bw_j);
                                var λg_j = Krg_j / (μg_j * Bg_j);

                                // Hydrostatic potential differences
                                var γo_avg = 0.5 * (γo_i + γo_j);
                                var γw_avg = 0.5 * (γw_i + γw_j);
                                var γg_avg = 0.5 * (γg_i + γg_j);
                                var DΦo = Po_i - Po_j - γo_avg * dz;
                                var DΦw = Pw_i - Pw_j - γw_avg * dz;
                                var DΦg = Pg_i - Pg_j - γg_avg * dz;

                                // Upstream phase mobilities & dissolved gas-oil ratio
                                var λo_up = Utility.SmoothUpstream(DΦo, λo_i, λo_j);
                                var λw_up = Utility.SmoothUpstream(DΦw, λw_i, λw_j);
                                var λg_up = Utility.SmoothUpstream(DΦg, λg_i, λg_j);
                                var Rs_up = Utility.SmoothUpstream(DΦo, Rs_i, Rs_j);

                                // Outflow flux from cell i to cell j
                                var FluxOil = transDt * λo_up * DΦo;
                                var FluxWater = transDt * λw_up * DΦw;
                                var FluxGas = transDt * λg_up * DΦg + Rs_up * FluxOil;

                                // Direct write to Flux Residual
                                Res[row_o] += FluxOil.Value;
                                Res[row_w] += FluxWater.Value;
                                Res[row_g] += FluxGas.Value;

                                // Direct write to Jacobian using precomputed face offsets
                                var fOff = face.Offset[c];

                                // Block (i, i) -> Derivatives wrt Cell i variables (slots 0, 1, 2)
                                JacValues[cOff.Diag[0]] += FluxOil.Derivative[0];   // row_o, col_p_i
                                JacValues[cOff.Diag[1]] += FluxOil.Derivative[1];   // row_o, col_sw_i
                                JacValues[cOff.Diag[2]] += FluxOil.Derivative[2];   // row_o, col_sg_i

                                JacValues[cOff.Diag[3]] += FluxWater.Derivative[0]; // row_w, col_p_i
                                JacValues[cOff.Diag[4]] += FluxWater.Derivative[1]; // row_w, col_sw_i
                                JacValues[cOff.Diag[5]] += FluxWater.Derivative[2]; // row_w, col_sg_i

                                JacValues[cOff.Diag[6]] += FluxGas.Derivative[0];  // row_g, col_p_i
                                JacValues[cOff.Diag[7]] += FluxGas.Derivative[1];  // row_g, col_sw_i
                                JacValues[cOff.Diag[8]] += FluxGas.Derivative[2];  // row_g, col_sg_i

                                // Block (i, j) -> Derivatives wrt Neighbor j variables (slots 3, 4, 5)
                                JacValues[fOff.Neighbor[0]] += FluxOil.Derivative[3];   // row_o, col_p_j
                                JacValues[fOff.Neighbor[1]] += FluxOil.Derivative[4];   // row_o, col_sw_j
                                JacValues[fOff.Neighbor[2]] += FluxOil.Derivative[5];   // row_o, col_sg_j

                                JacValues[fOff.Neighbor[3]] += FluxWater.Derivative[3]; // row_w, col_p_j
                                JacValues[fOff.Neighbor[4]] += FluxWater.Derivative[4]; // row_w, col_sw_j
                                JacValues[fOff.Neighbor[5]] += FluxWater.Derivative[5]; // row_w, col_sg_j

                                JacValues[fOff.Neighbor[6]] += FluxGas.Derivative[3];  // row_g, col_p_j
                                JacValues[fOff.Neighbor[7]] += FluxGas.Derivative[4];  // row_g, col_sw_j
                                JacValues[fOff.Neighbor[8]] += FluxGas.Derivative[5];  // row_g, col_sg_j
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
                    well.Qo_np1 = 0.0; well.Qw_np1 = 0.0; well.Ql_np1 = 0.0; well.Qg_np1 = 0.0;
                    well.Qo = 0.0; well.Qw = 0.0; well.Qg = 0.0;
                    var wOff = well.Offset;
                    if (!well.isActive)
                    {
                        JacValues[wOff.RowFlow_dPwell] += 1;
                        JacValues[wOff.RowConst_dQwell] += 1;
                        continue;
                    }

                    int rowflow = 3 * NumActiveGrids + 2 * wellidx + 0;  // Rate definition equation
                    int rowconst = 3 * NumActiveGrids + 2 * wellidx + 1; // Constraint equation

                    var pwell = new ADiff<Vector6d>(xnp1[rowflow], new(1, 0, 0, 0, 0, 0));
                    var qwell = new ADiff<Vector6d>(xnp1[rowconst], new(0, 1, 0, 0, 0, 0));
                    var FlowRes = qwell * dt;

                    int perfCount = well.PerfNatIndex.Count;
                    bool isWaterTarget = well.IsWaterTarget;
                    bool isOilTarget = well.IsOilTarget;
                    bool isGasTarget = well.IsGasTarget;
                    bool isResTarget = well.IsResTarget;
                    bool injectsGas = well.CurrentInjectorControl is { } injCtl &&
                        string.Equals(injCtl.InjectionPhase?.Trim(), "GAS", StringComparison.OrdinalIgnoreCase);

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
                                    int row_o = 3 * i + 0;
                                    int row_w = 3 * i + 1;
                                    int row_g = 3 * i + 2;

                                    var Po_m = new ADiff<Vector6d>(xnp1[row_o], new(0, 0, 0, 1, 0, 0));
                                    var Sw_m = new ADiff<Vector6d>(xnp1[row_w], new(0, 0, 0, 0, 1, 0));
                                    var Sg_m = new ADiff<Vector6d>(xnp1[row_g], new(0, 0, 0, 0, 0, 1));
                                    var Rs_m = mgrid.RsSat(Po_m);
                                    if (!mgrid.OilIsSaturated)
                                    {
                                        Sg_m = 0.0;
                                        Rs_m = new ADiff<Vector6d>(xnp1[row_g], new(0, 0, 0, 0, 0, 1));
                                    }

                                    var Pw_m = Po_m - mgrid.Pcow(Sw_m);
                                    var Pg_m = Po_m + mgrid.Pcog(Sg_m);

                                    var Bo_m = mgrid.Bopr(Po_m, Rs_m);
                                    var Bw_m = mgrid.Bw(Pw_m);
                                    var Bg_m = mgrid.Bg(Pg_m);
                                    var Kro_m = mgrid.Kro(Sw_m, Sg_m);
                                    var Krw_m = mgrid.Krw(Sw_m);
                                    var Krg_m = mgrid.Krg(Sg_m);
                                    var μo_m = mgrid.μopr(Po_m, Rs_m);
                                    var μw_m = mgrid.μw(Pw_m);
                                    var μg_m = mgrid.μg(Pg_m);
                                    var γo_m = mgrid.γopr(Po_m, Rs_m);
                                    var γw_m = mgrid.γw(Pw_m);
                                    var γg_m = mgrid.γg(Pg_m);

                                    var λo_m = Kro_m / (μo_m * Bo_m);
                                    var λw_m = Krw_m / (μw_m * Bw_m);
                                    var λg_m = Krg_m / (μg_m * Bg_m);

                                    var γo_well = mgrid.γopr(pwell, Rs_m);
                                    var γw_well = mgrid.γw(pwell);
                                    var γg_well = mgrid.γg(pwell);

                                    var γo_avg = 0.5 * (γo_m + γo_well);
                                    var γw_avg = 0.5 * (γw_m + γw_well);
                                    var γg_avg = 0.5 * (γg_m + γg_well);
                                    double dz = mgrid.Center.Z - well.Zref;

                                    var DΦo = Po_m - pwell - γo_avg * dz;
                                    var DΦw = Pw_m - pwell - γw_avg * dz;
                                    var DΦg = Pg_m - pwell - γg_avg * dz;

                                    var WIo = mgrid.WI * λo_m;
                                    var WIw = mgrid.WI * λw_m;
                                    var WIg = mgrid.WI * λg_m;

                                    var rate_oil = WIo * DΦo; well.Qo += rate_oil;
                                    var rate_water = WIw * DΦw; well.Qw += rate_water;
                                    var rate_gas = WIg * DΦg + rate_oil * Rs_m; well.Qg += rate_gas;

                                    // Reservoir residuals & Jacobian entries (outflow is +rate * dt)
                                    Res[row_o] += rate_oil.Value * dt;
                                    Res[row_w] += rate_water.Value * dt;
                                    Res[row_g] += rate_gas.Value * dt;

                                    JacValues[wOff.Rm_dPwell[p][0]] += rate_oil.Derivative[0] * dt;
                                    JacValues[wOff.Rm_dPwell[p][1]] += rate_water.Derivative[0] * dt;
                                    JacValues[wOff.Rm_dPwell[p][2]] += rate_gas.Derivative[0] * dt;

                                    var cOff = mgrid.Offset;
                                    JacValues[cOff.Diag[0]] += rate_oil.Derivative[3] * dt;
                                    JacValues[cOff.Diag[1]] += rate_oil.Derivative[4] * dt;
                                    JacValues[cOff.Diag[2]] += rate_oil.Derivative[5] * dt;
                                    JacValues[cOff.Diag[3]] += rate_water.Derivative[3] * dt;
                                    JacValues[cOff.Diag[4]] += rate_water.Derivative[4] * dt;
                                    JacValues[cOff.Diag[5]] += rate_water.Derivative[5] * dt;
                                    JacValues[cOff.Diag[6]] += rate_gas.Derivative[3] * dt;
                                    JacValues[cOff.Diag[7]] += rate_gas.Derivative[4] * dt;
                                    JacValues[cOff.Diag[8]] += rate_gas.Derivative[5] * dt;

                                    // Per-perforation grid derivatives for Row Q
                                    var perfRate = isWaterTarget ? rate_water
                                                 : isOilTarget ? rate_oil
                                                 : isGasTarget ? rate_gas
                                                 : isResTarget ? (rate_oil * Bo_m + rate_water * Bw_m + (rate_gas - rate_oil * Rs_m) * Bg_m)
                                                                 : (rate_oil + rate_water);

                                    FlowRes -= perfRate * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][0]] -= perfRate.Derivative[3] * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][1]] -= perfRate.Derivative[4] * dt;
                                    JacValues[wOff.RowFlow_dGrid[p][2]] -= perfRate.Derivative[5] * dt;
                                }

                                // Close production equation
                                Res[rowflow] += FlowRes.Value;
                                JacValues[wOff.RowFlow_dPwell] += FlowRes.Derivative[0];
                                JacValues[wOff.RowFlow_dQwell] += FlowRes.Derivative[1];
                                break;
                            }

                        case WellType.INJECTOR:
                            {
                                double PotentialGap(int p)
                                {
                                    var g = Grids[well.PerfNatIndex[p]];
                                    int a = g.ActiveIndex;
                                    double z = g.Center.Z - well.Zref, po = xnp1[3 * a];
                                    if (injectsGas)
                                    {
                                        double sg = g.OilIsSaturated ? xnp1[3 * a + 2] : 0.0;
                                        double pg = po + g.Pcog(sg).Value;
                                        return pg - pwell.Value - 0.5 * (g.γg(pg).Value + g.γg(pwell.Value).Value) * z;
                                    }
                                    double pw = po - g.Pcow(xnp1[3 * a + 1]).Value;
                                    return pw - pwell.Value - 0.5 * (g.γw(pw).Value + g.γw(pwell.Value).Value) * z;
                                }

                                int kept = -1; bool anyInjecting = false; double smallestGap = double.MaxValue;
                                for (int p = 0; p < perfCount; p++)
                                {
                                    if (!well.PerfIsOpen[p] || !Grids[well.PerfNatIndex[p]].isActive) continue;
                                    double gap = PotentialGap(p);
                                    if (gap < 0.0) anyInjecting = true;
                                    if (gap < smallestGap) { smallestGap = gap; kept = p; }
                                }
                                if (anyInjecting || well.ConstraintType == ConstraintType.MAX_BHP) kept = -1;

                                for (int p = 0; p < perfCount; p++)
                                {
                                    int m = well.PerfNatIndex[p];
                                    var mgrid = Grids[m];
                                    if (!well.PerfIsOpen[p] || !mgrid.isActive) continue;
                                    if (p != kept && PotentialGap(p) >= 0.0) continue;

                                    int i = mgrid.ActiveIndex;
                                    int row_o = 3 * i + 0;
                                    int row_w = 3 * i + 1;
                                    int row_g = 3 * i + 2;
                                    double dz = mgrid.Center.Z - well.Zref;

                                    var Po_m = new ADiff<Vector6d>(xnp1[row_o], new(0, 0, 0, 1, 0, 0));
                                    var Sw_m = new ADiff<Vector6d>(xnp1[row_w], new(0, 0, 0, 0, 1, 0));
                                    var Sg_m = new ADiff<Vector6d>(xnp1[row_g], new(0, 0, 0, 0, 0, 1));
                                    var Rs_m = mgrid.RsSat(Po_m);
                                    if (!mgrid.OilIsSaturated)
                                    {
                                        Sg_m = 0.0;
                                        Rs_m = new ADiff<Vector6d>(xnp1[row_g], new(0, 0, 0, 0, 0, 1));
                                    }

                                    var Pw_m = Po_m - mgrid.Pcow(Sw_m);
                                    var Pg_m = Po_m + mgrid.Pcog(Sg_m);

                                    var Bo_m = mgrid.Bopr(Po_m, Rs_m);
                                    var Bw_m = mgrid.Bw(Pw_m);
                                    var Bg_m = mgrid.Bg(Pg_m);
                                    var Kro_m = mgrid.Kro(Sw_m, Sg_m);
                                    var Krw_m = mgrid.Krw(Sw_m);
                                    var Krg_m = mgrid.Krg(Sg_m);
                                    var μo_m = mgrid.μopr(Po_m, Rs_m);
                                    var μw_m = mgrid.μw(Pw_m);
                                    var μg_m = mgrid.μg(Pg_m);
                                    var γw_m = mgrid.γw(Pw_m);
                                    var γg_m = mgrid.γg(Pg_m);

                                    var lambdaTotal = Krw_m / μw_m + Kro_m / μo_m + Krg_m / μg_m;

                                    if (!injectsGas)
                                    {
                                        var γw_well = mgrid.γw(pwell);
                                        var γw_avg = 0.5 * (γw_m + γw_well);
                                        var DΦw = Pw_m - pwell - γw_avg * dz;
                                        if (UseUpstreamWeightingForInjection)
                                            lambdaTotal = mgrid.Krw(1.0) / μw_m;

                                        var WIw = mgrid.WI * lambdaTotal / Bw_m;
                                        var rate_water = WIw * DΦw; well.Qw += rate_water;

                                        Res[row_w] += rate_water.Value * dt;
                                        JacValues[wOff.Rm_dPwell[p][1]] += rate_water.Derivative[0] * dt;

                                        var cOff = mgrid.Offset;
                                        JacValues[cOff.Diag[3]] += rate_water.Derivative[3] * dt;
                                        JacValues[cOff.Diag[4]] += rate_water.Derivative[4] * dt;
                                        JacValues[cOff.Diag[5]] += rate_water.Derivative[5] * dt;

                                        // rate_water is already in stb/day; multiply by Bw_m only for RESV target
                                        var perfRate = isResTarget ? rate_water * Bw_m : rate_water;

                                        FlowRes -= perfRate * dt;
                                        JacValues[wOff.RowFlow_dGrid[p][0]] -= perfRate.Derivative[3] * dt;
                                        JacValues[wOff.RowFlow_dGrid[p][1]] -= perfRate.Derivative[4] * dt;
                                        JacValues[wOff.RowFlow_dGrid[p][2]] -= perfRate.Derivative[5] * dt;
                                    }
                                    else
                                    {
                                        var γg_well = mgrid.γg(pwell);
                                        var γg_avg = 0.5 * (γg_m + γg_well);
                                        var DΦg = Pg_m - pwell - γg_avg * dz;
                                        if (UseUpstreamWeightingForInjection)
                                            lambdaTotal = mgrid.Krg(1.0) / μg_m;

                                        var WIg = mgrid.WI * lambdaTotal / Bg_m;
                                        var rate_gas = WIg * DΦg; well.Qg += rate_gas;

                                        Res[row_g] += rate_gas.Value * dt;
                                        JacValues[wOff.Rm_dPwell[p][2]] += rate_gas.Derivative[0] * dt;

                                        var cOff = mgrid.Offset;
                                        JacValues[cOff.Diag[6]] += rate_gas.Derivative[3] * dt;
                                        JacValues[cOff.Diag[7]] += rate_gas.Derivative[4] * dt;
                                        JacValues[cOff.Diag[8]] += rate_gas.Derivative[5] * dt;

                                        // rate_gas is already in Mscf/day; multiply by Bg_m only for RESV target
                                        var perfRate = isResTarget ? rate_gas * Bg_m : rate_gas;

                                        FlowRes -= perfRate * dt;
                                        JacValues[wOff.RowFlow_dGrid[p][0]] -= perfRate.Derivative[3] * dt;
                                        JacValues[wOff.RowFlow_dGrid[p][1]] -= perfRate.Derivative[4] * dt;
                                        JacValues[wOff.RowFlow_dGrid[p][2]] -= perfRate.Derivative[5] * dt;
                                    }
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
                    well.Ql_np1 = well.Qo_np1 + well.Qw_np1;
                    well.Qg_np1 = well.Qg.Value;
                    well.Pwf_np1 = pwell.Value;

                    // Close the Well Constraint Equation (Row p)
                    var wellConstraint = well.Constraint(pwell, qwell);
                    Res[rowconst] += wellConstraint.Value;
                    JacValues[wOff.RowConst_dPwell] += wellConstraint.Derivative[0];
                    JacValues[wOff.RowConst_dQwell] += wellConstraint.Derivative[1];
                }

                // =========================================================================
                // 3. Compute Normalized Residual
                // =========================================================================

                // At the very end of Residual(xs, tnp1), after accumulation + faces + wells are all in Res:

                double MaxCNV = 0.0;
                foreach (var grid in Grids)
                {
                    if (!grid.isActive) continue;

                    int i = grid.ActiveIndex;
                    double cnvOil = Abs(Res[3 * i + 0]) * grid.Bo_np1 * grid.InvPoreVolume;
                    double cnvWater = Abs(Res[3 * i + 1]) * grid.Bw_np1 * grid.InvPoreVolume;
                    double cnvGas = Abs(Res[3 * i + 2]) * grid.Bg_np1 * grid.InvPoreVolume;

                    if (cnvOil > MaxCNV) MaxCNV = cnvOil;
                    if (cnvWater > MaxCNV) MaxCNV = cnvWater;
                    if (cnvGas > MaxCNV) MaxCNV = cnvGas;
                }

                int wellOffset = 3 * NumActiveGrids;
                for (int wellidx = 0; wellidx < NumWells; wellidx++)
                {
                    var well = Wells[wellidx];
                    if (!well.isActive) continue;

                    bool isProducer = well.WellType == WellType.PRODUCER;
                    int row0 = wellOffset + 2 * wellidx + 0;
                    int row1 = wellOffset + 2 * wellidx + 1;

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
                                _ => Max(prod.OilRateTarget, Max(prod.WaterRateTarget, Max(prod.GasRateTarget, prod.LiquidRateTarget)))
                            };
                        }
                        else if (!isProducer && well.CurrentInjectorControl is { } inj)
                        {
                            bhpLimit = inj.MaxBhpLimit;
                            targetRate = inj.TargetRate;
                        }
                    }

                    double currentFlowRate = Abs(well.Qo.Value) + Abs(well.Qw.Value) + Abs(well.Qg.Value);
                    double rateScale = Max(1.0, Max(Abs(targetRate), currentFlowRate));

                    // Row 0 (FlowRes) has units of [stb] or [Mscf] because it is multiplied by dt
                    double err0 = Abs(Res[row0]) / (rateScale * Max(1e-12, dt));

                    // Row 1 (well.Constraint) has units of [psi] for BHP control or [stb/day | Mscf/day] for rate control
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
                                .Concat(Wells.SelectMany(well => well.GetEventTimes().Select(date=>date.AddSeconds(-1))))
                                .Concat(ResultTime.Select(t => startDate.AddDays(t).AddSeconds(-1)))
                                .Distinct()
                                .Order()];

            Console.WriteLine($"""
                Total Grids = {Ngrids}, Active Grids = {NumActiveGrids}
                ======================================================================
                Starting simulation
                """);

            double L = 0, tnp1, rnorm = 0;
            bool tryingnewtimestep = true, isComplete;
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
                        foreach (var well in Wells)
                            well.UpdateControl(currentSimDate);

                        foreach (var well in Wells)
                            well.UpdateCompletion(currentSimDate);

                        double cumt = 0.0;
                        int stepIdx = DATA.TimeSteps.FindIndex(dt => (cumt += dt) >= tnp1 - 1e-9);
                        if (stepIdx < 0) stepIdx = DATA.TimeSteps.Count - 1;
                        dRsdt = GetDrsdtForStep(stepIdx);
                    }

                    // Prepare initial guess x^(0) from converged state n
                    InitialGuess(xs, dt);

                    Console.WriteLine($"""
                         iter |   Residual Norm
                        ------+------------------
                        """);
                    int iter;
                    bool isConverged = false;
                    bool phaseSwitched = false;
                    double prevRnorm = double.PositiveInfinity;
                    double maxDPo = 0.0, maxDSw = 0.0, maxDRs = 0.0;

                    for (iter = 0; iter < 15; iter++)
                    {
                        rnorm = Residual(xs, dt);
                        Console.WriteLine($"  {iter,2}  |  {rnorm:E4}");

                        // Early exit if Newton hits NaN or diverges (skip spike check right after a phase switch)
                        if (double.IsNaN(rnorm) || (!phaseSwitched && iter > 2 && rnorm > 10.0 * prevRnorm))
                            break;

                        // Converged at current xs
                        if (rnorm < 1e-6 && !phaseSwitched &&
                            (iter == 0 || (maxDSw < 1e-4 && maxDPo < 1e-2 && maxDRs < 1e-4)))
                        {
                            isConverged = true;
                            break;
                        }

                        prevRnorm = rnorm;

                        // Direct solve via MKL PARDISO: J * dx = Res
                        solver.Solve(JacValues, Res, dx);

                        // Update primary variables: x^{k+1} = x^k - dx
                        (maxDPo, maxDSw, maxDRs, phaseSwitched) =
                            ApplyPrimaryVariableUpdate(xs, dx, dt, NumActiveGrids, NumWells);
                    }

                    // Cut dt on Newton non-convergence
                    if (!isConverged)
                    {
                        dtCarry = 0.5 * dt;
                        tryingnewtimestep = true;
                        if (dtCarry < 1e-12) throw new Exception("Time step is too small");
                        Console.WriteLine("""
                            ============================
                            Rejected (Non-Convergence)
                            ============================
                            """);
                        continue;
                    }

                    // Validate maximum per-step state changes
                    var (maxDP, maxDw, maxDg, i_p, i_w, i_g) = CheckSolution(xs);
                    if (maxDP > 500 || maxDw > 0.2 || maxDg > 0.2)
                    {
                        int j = maxDP > 500 ? i_p : maxDw > 0.2 ? i_w : i_g;
                        double limit = maxDP > 500 ? 500 : 0.2;
                        double change = maxDP > 500 ? maxDP : maxDw > 0.2 ? maxDw : maxDg;
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
                            Rejected (Minimum Pressure Violated: {well.Pwf.Value:F2} < {prod.MinBhpLimit:F2} psi) @ {well.Name}
                            ===============================================================
                            """);
                            well.Control_is_Violated = true;
                        }

                        if (well.WellType == WellType.INJECTOR &&
                            well.CurrentInjectorControl is { } inj &&
                            well.ConstraintType != ConstraintType.MAX_BHP &&
                            inj.MaxBhpLimit > 0.0 &&
                            well.Pwf.Value > (inj.MaxBhpLimit + 1e-2))
                        {
                            well.ConstraintType = ConstraintType.MAX_BHP;
                            Console.WriteLine($"""
                                ===============================================================
                                Rejected (Maximum Pressure Violated: {well.Pwf.Value:F2} > {inj.MaxBhpLimit:F2} psi) @ {well.Name}
                                ===============================================================
                                """);
                            well.Control_is_Violated = true;
                        }
                    }
                    if (Wells.Any(well => well.Control_is_Violated))
                    {
                        foreach (var well in Wells)
                            well.Control_is_Violated = false;
                        continue;
                    }

                    // Commit verified state n+1 -> n
                    t = tnp1;
                    Time.Add(t);
                    foreach (var grid in Grids)
                    {
                        if (!grid.isActive) continue;
                        grid.ErSw_Bw_n = grid.ErSw_Bw_np1.Value;
                        grid.ErSo_Bo_n = grid.ErSo_Bo_np1.Value;
                        grid.ErSg_Bg_n = grid.ErSg_Bg_np1.Value;
                        grid.Pw_n = grid.Pw_np1; grid.Sw_n = grid.Sw_np1;
                        grid.Po_n = grid.Po_np1; grid.So_n = grid.So_np1;
                        grid.Pg_n = grid.Pg_np1; grid.Sg_n = grid.Sg_np1; grid.Rs_n = grid.Rs_np1;
                        grid.SaveResult(po: grid.Po_n, pw: grid.Pw_n, so: grid.So_n, sw: grid.Sw_n);
                    }
                    foreach (var well in Wells)
                    {
                        well.Pwf_n = well.Pwf_np1; 
                        well.Qo_n = well.Qo_np1; well.Qw_n = well.Qw_np1;
                        well.Ql_n = well.Ql_np1; well.Qg_n = well.Qg_np1;
                        well.SaveResult(pwf: well.Pwf_n, qw: well.Qw_n, qo: well.Qo_n, qg: well.Qg_n);
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

            var animatorExporter = new SimulationAnimatorExporter<Vector6d>(outputDir, Grids);
            for (int n = 0; n < Time.Count; n++)
                animatorExporter.WriteTimeStep(Time, Grids, n);


            foreach (var well in Wells)
            {
                List<ColVec> vectors = [Time, well.PWF, well.QO, well.QW];
                string filename = Path.Combine(outputDir, $"{well.Name}.txt");
                SepalSolver.Math.WriteMatrix(vectors, filename);
            }
        }
    }
}
