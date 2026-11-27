using EclipseDeckReader;
using SepalSolver;
using static System.Math;

namespace ReservoirSimulator
{
    /// <summary>
    /// A Fetkovich aquifer (AQUFETP + AQUANCON): a tank of water with a productivity index, connected to the faces of
    /// some boundary cells. Inside a step the aquifer pressure is the one at the start of the step, and the influx is
    /// integrated over the step exactly (Fetkovich 1971):
    ///
    ///     q_c = w_c * J_eff * (Pa + γw (z_c - z_datum) - Pw_c),     J_eff = ct V / dt * (1 - exp(-J dt / (ct V)))
    ///
    /// q_c is water at surface conditions into cell c (stb/d), w_c the cell's share of the connection (face area times
    /// multiplier, or the influx coefficient, over the sum for the aquifer). When a step is accepted, the water that
    /// went in is taken out of the aquifer: We += Σ q_c dt and Pa = Pi - We / (ct V).
    /// Everything is held in FIELD units, like the rest of the engine.
    /// </summary>
    public sealed class FetkovichAquifer
    {
        public int Id { get; init; }
        public double DatumDepth { get; init; }               // ft
        public double InitialPressure { get; private set; }   // psia
        public double WaterVolume { get; init; }              // stb
        public double TotalCompressibility { get; init; }     // 1/psi
        public double ProductivityIndex { get; init; }        // stb/d/psi
        public bool PressureDefaulted { get; init; }

        /// <summary>Connected cells: index into the reservoir's grid array, share of the influx, cell centre depth (ft).</summary>
        public List<(int Grid, double Weight, double Depth)> Connections { get; } = [];

        public double Pressure { get; private set; }          // aquifer pressure at the start of the current step
        public double CumulativeInflux { get; private set; }  // stb
        public double Rate { get; private set; }              // stb/d over the last accepted step

        // One entry per saved time (the first is the initial state), alongside the reservoir's Time list.
        public List<double> RateHistory { get; } = [];
        public List<double> CumulativeHistory { get; } = [];
        public List<double> PressureHistory { get; } = [];

        double CtV => TotalCompressibility * WaterVolume;

        public double EffectiveIndex(double dt) =>
            CtV > 0 && dt > 0 ? CtV / dt * (1.0 - Exp(-ProductivityIndex * dt / CtV)) : ProductivityIndex;

        /// <summary>Water into one connected cell this step (stb/d), with derivatives through Pw and γw.</summary>
        public ADiff<T> Influx<T>(ADiff<T> pw, ADiff<T> γw, double weight, double depth, double dt) where T : struct, IVectorND<T> =>
            weight * EffectiveIndex(dt) * (γw * (depth - DatumDepth) - pw + Pressure);

        double Influx(double pw, double γw, double weight, double depth, double dt) =>
            weight * EffectiveIndex(dt) * (Pressure + γw * (depth - DatumDepth) - pw);

        /// <summary>
        /// Sets the starting state once the reservoir is in equilibrium. A defaulted initial pressure is the
        /// reservoir's own water pressure carried to the datum, so the aquifer starts at rest with the field.
        /// </summary>
        public void Initialize(Func<int, (double Pw, double γw)> water)
        {
            if (PressureDefaulted && Connections.Count > 0)
                InitialPressure = Connections.Sum(c => { var (pw, g) = water(c.Grid); return c.Weight * (pw - g * (c.Depth - DatumDepth)); });
            Pressure = InitialPressure;
            CumulativeInflux = 0; Rate = 0;
            RateHistory.Clear(); CumulativeHistory.Clear(); PressureHistory.Clear();
            Save();
        }

        /// <summary>Takes the water of an accepted step out of the aquifer.</summary>
        public void Accept(Func<int, (double Pw, double γw)> water, double dt)
        {
            Rate = Connections.Sum(c => { var (pw, g) = water(c.Grid); return Influx(pw, g, c.Weight, c.Depth, dt); });
            CumulativeInflux += Rate * dt;
            Pressure = CtV > 0 ? InitialPressure - CumulativeInflux / CtV : InitialPressure;
            Save();
        }

        void Save()
        {
            RateHistory.Add(Rate);
            CumulativeHistory.Add(CumulativeInflux);
            PressureHistory.Add(Pressure);
        }

        /// <summary>
        /// Builds the deck's Fetkovich aquifers on a grid. Cells are matched by their (i, j, k) and connected through
        /// the named face when that face is on the edge of the active grid (or always, with AllowInternalFaces).
        /// </summary>
        public static List<FetkovichAquifer> Build<T>(EclipseDataDeck deck, Grid<T>[] grids, Units units) where T : struct, IVectorND<T>
        {
            var aquifers = new List<FetkovichAquifer>();
            if (deck.FetkovichAquifers.Count == 0)
                return aquifers;
            if (units is not (Units.FIELD or Units.METRIC))
                throw new NotSupportedException("Fetkovich aquifers are read in FIELD or METRIC units only.");

            double psiPerUnit = Utility.Convert(1.0, "Pressure", units, Units.FIELD) - Utility.Convert(0.0, "Pressure", units, Units.FIELD);
            double stbPerUnit = units == Units.METRIC ? 6.289811 : 1.0;   // sm3 -> stb
            double ft2PerUnit = Pow(Utility.Convert(1.0, "Length", units, Units.FIELD), 2);

            foreach (var a in deck.FetkovichAquifers)
            {
                var aquifer = new FetkovichAquifer
                {
                    Id = a.Id,
                    DatumDepth = Utility.Convert(a.DatumDepth, "Depth", units, Units.FIELD),
                    InitialPressure = a.InitialPressure is { } p ? Utility.Convert(p, "Pressure", units, Units.FIELD) : 0.0,
                    PressureDefaulted = a.InitialPressure is null,
                    WaterVolume = a.WaterVolume * stbPerUnit,
                    TotalCompressibility = a.TotalCompressibility / psiPerUnit,
                    ProductivityIndex = a.ProductivityIndex * stbPerUnit / psiPerUnit,
                };

                var raw = new List<(int Grid, double Coefficient, double Depth)>();
                foreach (var con in deck.AquiferConnections.Where(c => c.AquiferId == a.Id))
                {
                    int face = con.Face switch { "I-" or "X-" => 0, "I+" or "X+" => 1, "J-" or "Y-" => 2, "J+" or "Y+" => 3, "K-" or "Z-" => 4, "K+" or "Z+" => 5, _ => -1 };
                    if (face < 0) throw new ArgumentException($"AQUANCON: unknown face '{con.Face}' for aquifer {a.Id}.");
                    for (int m = 0; m < grids.Length; m++)
                    {
                        var g = grids[m];
                        if (!g.isActive) continue;
                        var (i, j, k) = g.IJKIndex;
                        if (i + 1 < con.I1 || i + 1 > con.I2 || j + 1 < con.J1 || j + 1 > con.J2 || k + 1 < con.K1 || k + 1 > con.K2) continue;
                        var f = g.Faces[face];
                        bool onEdge = f.GridIndices.All(n => !grids[n].isActive);
                        if (!onEdge && !con.AllowInternalFaces) continue;
                        double coefficient = (con.InfluxCoefficient is { } c ? c * ft2PerUnit : f.Area) * con.InfluxMultiplier;
                        raw.Add((m, coefficient, g.Center.Z));
                    }
                }
                double total = raw.Sum(r => r.Coefficient);
                foreach (var (m, c, z) in raw)
                    aquifer.Connections.Add((m, total > 0 ? c / total : 1.0 / raw.Count, z));
                Console.WriteLine($"      [Aquifer] Fetkovich aquifer {a.Id}: {aquifer.Connections.Count} connected cell(s), J {aquifer.ProductivityIndex:G4} stb/d/psi, ct·V {aquifer.TotalCompressibility * aquifer.WaterVolume:G4} stb/psi.");
                aquifers.Add(aquifer);
            }
            return aquifers;
        }

        /// <summary>The deck's unit system, read the way the reservoir models read it.</summary>
        public static Units UnitsOf(EclipseDataDeck deck) => deck.UnitSystem.ToUpper() switch
        {
            "FIELD" => Units.FIELD,
            "METRIC" => Units.METRIC,
            "LAB" => Units.LAB,
            "PVT-M" => Units.PVT_M,
            _ => throw new InvalidOperationException($"Unsupported unit system: {deck.UnitSystem}")
        };

        /// <summary>For each grid index, the aquifers connected to it and the cell's share and depth.</summary>
        public static Dictionary<int, (FetkovichAquifer Aquifer, double Weight, double Depth)[]> ByCell(IEnumerable<FetkovichAquifer> aquifers) =>
            aquifers.SelectMany(a => a.Connections.Select(c => (c.Grid, Entry: (a, c.Weight, c.Depth))))
                    .GroupBy(x => x.Grid)
                    .ToDictionary(g => g.Key, g => g.Select(x => x.Entry).ToArray());
    }
}
