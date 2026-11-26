using EclipseDeckReader;
using System;
using System.Linq;

namespace ReservoirSimulator
{
    public class Reservoir : IReservoir
    {
        private readonly IReservoir _reservoir;

        public bool UseUpstreamWeightingForInjection
        {
            get => _reservoir.UseUpstreamWeightingForInjection;
            set => _reservoir.UseUpstreamWeightingForInjection = value;
        }

        public Reservoir(EclipseDataDeck Data)
        {
            ArgumentNullException.ThrowIfNull(Data);

            bool hasOil = Data.PhaseFlags.Contains("OIL", StringComparer.OrdinalIgnoreCase);
            bool hasWater = Data.PhaseFlags.Contains("WATER", StringComparer.OrdinalIgnoreCase);
            bool hasGas = Data.PhaseFlags.Contains("GAS", StringComparer.OrdinalIgnoreCase);
            bool hasDisgas = Data.PhaseFlags.Contains("DISGAS", StringComparer.OrdinalIgnoreCase) || Data.RawPvtoTables.Count > 0;
            bool hasVapoil = Data.PhaseFlags.Contains("VAPOIL", StringComparer.OrdinalIgnoreCase) || Data.RawPvtgTables.Count > 0;

            if (!hasOil && !hasWater && !hasGas)
                throw new NotSupportedException("No active fluid phases (OIL, WATER, GAS) specified in the Eclipse data deck.");

            if (hasDisgas && !hasOil)
                throw new InvalidOperationException("DISGAS requires the OIL phase to be active.");

            if (hasVapoil && !hasGas)
                throw new InvalidOperationException("VAPOIL requires the GAS phase to be active.");

            // Route to the 4 simulator architectures based on active phases and mass-transfer flags
            if (!hasGas)
            {
                // 2-phase Oil-Water (or single-phase Oil / Water with zero saturation of the other)
                _reservoir = new OilWaterReservoir(Data);
            }
            else if (hasDisgas && hasVapoil)
            {
                // 3-phase (or 2-phase Oil-Gas) Volatile / Condensate: Live Oil + Wet Gas
                _reservoir = new WetGasLiveOilReservoir(Data);
            }
            else if (hasDisgas)
            {
                // 3-phase (or 2-phase Oil-Gas) Black Oil: Live Oil + Dry Gas (e.g., SPE1)
                _reservoir = new DryGasLiveOilReservoir(Data);
            }
            else if (hasVapoil)
            {
                // Wet Gas / Gas Condensate with Dead Oil (PVDO + PVTG)
                _reservoir = new WetGasDeadOilReservoir(Data);
            }
            else
            {
                // Dead Oil + Dry Gas (3-phase Oil-Water-Gas, 2-phase Oil-Gas, or 2-phase Gas-Water)
                _reservoir = new DryGasDeadOilReservoir(Data);
            }
        }

        public void Initialize() => _reservoir.Initialize();
        public void Simulate(double[] ResultTime = null) => _reservoir.Simulate(ResultTime);
        public void ExportResults(string outputDir) => _reservoir.ExportResults(outputDir);
    }
}