using System;
using System.Collections.Generic;
using System.Linq;

namespace EclipseDeckReader
{
    public static class DeckPrinter
    {
        /// <summary>
        /// Analyzes a populated EclipseDataDeck instance and outputs a structured summary to the console.
        /// </summary>
        public static void PrintSummary(EclipseDataDeck deck)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("                  ECLIPSE DATA DECK SUMMARY                       ");
            Console.WriteLine("==================================================================");

            // 1. Grid Structural Block
            int totalCells = deck.NX * deck.NY * deck.NZ;
            Console.WriteLine("## GRID GEOMETRY");
            Console.WriteLine($"* Dimensions (NX x NY x NZ) : {deck.NX} x {deck.NY} x {deck.NZ}");
            Console.WriteLine($"* Total Grid Block Cells    : {totalCells:N0}");
            Console.WriteLine();

            // 2. Active Fluid Domain Phases
            Console.WriteLine("## FLUID SYSTEMS AND TRACKED PHASES");
            Console.WriteLine($"* Active Keyword Flags      : {string.Join(", ", deck.PhaseFlags)}");
            Console.WriteLine($"* Surface Oil Baseline Dens : {deck.SurfaceOilDensity:F2}");
            Console.WriteLine($"* Surface Water Baseline Dens: {deck.SurfaceWaterDensity:F2}");
            Console.WriteLine($"* Surface Gas Baseline Dens  : {deck.SurfaceGasDensity:F4}");
            Console.WriteLine();

            // 3. Spatial Parameter Arrays Profiler
            Console.WriteLine("## CONTINUOUS RECURRENT SPATIAL ARRAYS");
            PrintArrayStats("Porosity (PORO)", deck.PorosityArray);
            PrintArrayStats("Net-To-Gross (NTG)", deck.NtgArray);
            PrintArrayStats("Permeability X (PERMX)", deck.PermXArray);
            PrintArrayStats("Permeability Y (PERMY)", deck.PermYArray);
            PrintArrayStats("Permeability Z (PERMZ)", deck.PermZArray);
            PrintArrayStats("Transmissibility Mult X (MULTX)", deck.MultXArray);
            PrintArrayStats("Transmissibility Mult Y (MULTY)", deck.MultYArray);
            PrintArrayStats("Transmissibility Mult Z (MULTZ)", deck.MultZArray);
            PrintArrayStats("Coordinate Lines (COORD)", deck.CoordArray);
            PrintArrayStats("Corner Heights (ZCORN)", deck.ZcornArray);
            Console.WriteLine();

            // 4. Discrete Assignment Regions
            Console.WriteLine("## DISCRETE ASSIGNMENT REGIONS");
            PrintRegionStats("Active Cell Flags (ACTNUM)", deck.ActnumArray);
            PrintRegionStats("Saturation Regions (SATNUM)", deck.SatnumArray);
            PrintRegionStats("PVT Properties Regions (PVTNUM)", deck.PvtnumArray);
            PrintRegionStats("Equilibrium Regions (EQUINUM)", deck.EquinumArray);
            Console.WriteLine();

            // 5. Multi-Region Hydrostatic Equilibration Summary block
            Console.WriteLine("## RESERVOIR INITIAL HYDROSTATIC EQUILIBRIUM");
            Console.WriteLine($"* Registered EQUIL Tables Count: {deck.RawEquilTables.Count}");
            for (int e = 0; e < deck.RawEquilTables.Count; e++)
            {
                var eq = deck.RawEquilTables[e];
                if (eq.Count >= 5)
                {
                    Console.WriteLine($"  └─ Equilibrium Region {e + 1} -> Datum Z: {eq[0]:F1} | Datum P: {eq[1]:F2} | WOC Z: {eq[2]:F1} | GOC Z: {eq[4]:F1}");
                }
            }
            Console.WriteLine();

            // 6. Rock Mechanical Properties
            Console.WriteLine("## ROCK MECHANICS AND PORE EXPANSION");
            Console.WriteLine($"* Reference Pressure (Pref)   : {deck.RockReferencePressure:F2} psi");
            Console.WriteLine($"* Formation Compressibility (cr): {deck.RockCompressibility:E4} 1/psi");
            Console.WriteLine();

            // 7. Fault Multipliers
            Console.WriteLine("## FAULT MULTIPLIERS");
            Console.WriteLine($"* Registered Faults Count  : {deck.Faults.Count}");
            foreach (var kvp in deck.Faults)
            {
                Console.WriteLine($"  └─ Fault: {kvp.Key,-12} | Face: {kvp.Value.Face} | Mult: {kvp.Value.TransmissibilityMultiplier:F4}");
            }
            Console.WriteLine();

            // 8. Dynamic Well Management Block
            Console.WriteLine("## WELL GEOMETRY AND OPERATIONS NETWORK");
            Console.WriteLine($"* Total Wells Registered    : {deck.Wells.Count}");
            if (deck.Wells.Count == 0)
            {
                Console.WriteLine("* [NO WELL NETWORKS REGISTERED OR PARSED]");
            }
            else
            {
                foreach (var wellKvp in deck.Wells)
                {
                    string name = wellKvp.Key;
                    WellData data = wellKvp.Value;

                    Console.WriteLine($"\n* Well Entity: {name}");

                    // A. Master Well Specifications (WELSPECS)
                    if (data.Specification != null)
                    {
                        var spec = data.Specification;
                        Console.WriteLine($"  ├─ Group Reference : {spec.GroupName,-10} | Head (I,J): ({spec.HeadI}, {spec.HeadJ}) | Ref Depth: {spec.ReferenceDepth:F1} | Primary Phase: {spec.PreferredPhase}");
                    }
                    else
                    {
                        Console.WriteLine("  ├─ Group Reference : [NO INITIAL SPECIFICATION RECORD LOADED]");
                    }

                    // B. Wellbore Interval Completions
                    Console.WriteLine($"  ├─ Wellbore Completions Count: {data.CompletionHistory.Count}");
                    for (int c = 0; c < data.CompletionHistory.Count; c++)
                    {
                        var comp = data.CompletionHistory[c];
                        Console.WriteLine($"  │  └─ Comp [{c + 1}] [{comp.EffectiveDate:yyyy-MM-dd}]: Cell (I, J): ({comp.I + 1}, {comp.J + 1}) | K-Layers: {comp.KStart}-{comp.KEnd} | Status: {comp.Status} | Diameter: {comp.WellBoreDiameter:F3} | Skin: {comp.EffectiveSkin:F2}");
                    }

                    // D. Producer Control Schedules (WCONPROD)
                    Console.WriteLine($"  ├─ Producer Control Schedules (WCONPROD): {data.ProducerControlHistory.Count}");
                    var sortedProdControls = data.ProducerControlHistory.OrderBy(ch => ch.EffectiveDate).ToList();
                    for (int h = 0; h < sortedProdControls.Count; h++)
                    {
                        var ctrl = sortedProdControls[h];
                        Console.WriteLine($"  │  └─ Date [{ctrl.EffectiveDate:yyyy-MM-dd}] Status: {ctrl.Status,-4} | Mode: {ctrl.ControlMode,-4} | Oil: {ctrl.OilRateTarget:F1} | Water: {ctrl.WaterRateTarget:F1} | Gas: {ctrl.GasRateTarget:F1} | Liq: {ctrl.LiquidRateTarget:F1} | Resv: {ctrl.ResvRateTarget:F1} | Min BHP: {ctrl.MinBhpLimit:F1}");
                    }

                    // E. Injector Control Schedules (WCONINJE)
                    Console.WriteLine($"  └─ Injector Control Schedules (WCONINJE): {data.InjectorControlHistory.Count}");
                    var sortedInjControls = data.InjectorControlHistory.OrderBy(ch => ch.EffectiveDate).ToList();
                    for (int h = 0; h < sortedInjControls.Count; h++)
                    {
                        var ctrl = sortedInjControls[h];
                        Console.WriteLine($"     └─ Date [{ctrl.EffectiveDate:yyyy-MM-dd}] Phase: {ctrl.InjectionPhase,-5} | Status: {ctrl.Status,-4} | Mode: {ctrl.ControlMode,-4} | Rate: {ctrl.TargetRate:F1} | Resv: {ctrl.ResvRateTarget:F1} | Max BHP: {ctrl.MaxBhpLimit:F1}");
                    }
                }
            }
            Console.WriteLine();

            // 9. Fluid Properties Summary Block
            Console.WriteLine("## FLUID PVT REGIONS COUNT");
            PrintTableGroupOverview("Oil Live PVT Properties (PVTO)", deck.RawPvtoTables);
            PrintTableGroupOverview("Oil Dead PVT Properties (PVDO)", deck.RawPvdoTables);
            PrintTableGroupOverview("Water PVT Properties (PVTW)", deck.RawPvtwTables);
            PrintTableGroupOverview("Gas Live PVT Properties (PVTG)", deck.RawPvtgTables);
            PrintTableGroupOverview("Gas Dead PVT Properties (PVDG)", deck.RawPvdgTables);
            Console.WriteLine();

            // 10. Rock-Fluid Properties Block
            Console.WriteLine("## ROCK-FLUID INTERFACIAL RELATIONS");
            PrintTableGroupOverview("Water-Oil Saturation Functions (SWOF)", deck.RawSwofTables);
            PrintTableGroupOverview("Gas-Oil Saturation Functions (SGOF)", deck.RawSgofTables);
            Console.WriteLine("==================================================================");

            // 11. Schedule Timeline Overview
            Console.WriteLine($"* Schedule Timeline                  : Start = {deck.StartDate:yyyy-MM-dd} | End = {deck.CurrentSimDate:yyyy-MM-dd}");
            Console.WriteLine($"  └─ Total Report Steps              : {deck.TimeSteps.Count} (Cumulative Duration = {deck.TimeSteps.Sum():F1} days)");
        }

        private static void PrintArrayStats(string arrayName, List<double> array)
        {
            if (array == null || array.Count == 0)
            {
                Console.WriteLine($"* {arrayName,-35}: [NOT PRESENT OR EMPTY]");
                return;
            }
            double min = array.Min();
            double max = array.Max();
            double avg = array.Average();
            Console.WriteLine($"* {arrayName,-35}: Loaded Count = {array.Count:N0} | Min = {min,7:F4} | Max = {max,7:F4} | Avg = {avg,7:F4}");
        }

        private static void PrintRegionStats(string arrayName, List<int> array)
        {
            if (array == null || array.Count == 0)
            {
                Console.WriteLine($"* {arrayName,-35}: [NOT PRESENT OR EMPTY]");
                return;
            }
            int min = array.Min();
            int max = array.Max();
            int uniqueRegions = array.Distinct().Count();
            Console.WriteLine($"* {arrayName,-35}: Loaded Count = {array.Count:N0} | Distinct Regions = {uniqueRegions} (IDs: {min} to {max})");
        }

        private static void PrintTableGroupOverview(string tableName, List<List<double>> tables)
        {
            if (tables == null || tables.Count == 0)
            {
                Console.WriteLine($"* {tableName,-35}: [0 Tables Found]");
                return;
            }

            Console.WriteLine($"* {tableName,-35}: [{tables.Count} Table Region(s) Registered]");
            for (int r = 0; r < tables.Count; r++)
            {
                Console.WriteLine($"  └─ Region {r + 1}: Parsed Numeric Values Flat Elements Count = {tables[r].Count}");
            }
        }

        private static void PrintTableGroupOverview(string tableName, List<List<List<double>>> tables, string entryLabel = "Sub-Tables / Rows")
        {
            if (tables == null || tables.Count == 0)
            {
                Console.WriteLine($"* {tableName,-35}: [0 Tables Found]");
                return;
            }

            Console.WriteLine($"* {tableName,-35}: [{tables.Count} Table Region(s) Registered]");
            for (int r = 0; r < tables.Count; r++)
            {
                var region = tables[r];
                int totalElements = region.Sum(sub => sub.Count);
                Console.WriteLine($"  └─ Region {r + 1}: {region.Count} {entryLabel}, Total Numeric Elements = {totalElements}");
            }
        }
    }
}