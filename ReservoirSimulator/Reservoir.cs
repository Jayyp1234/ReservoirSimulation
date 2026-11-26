using EclipseDeckReader;
using SepalSolver;
using System;
using System.Collections.Generic;
using System.Text;
using static TorchSharp.torch.utils;

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
            _reservoir = new DryGasLiveOilReservoir(Data);     // implements IReservoir<Vector4d>


            //if (Data == null)
            //    throw new ArgumentNullException(nameof(Data), "EclipseDataDeck cannot be null.");
            //// Determine the type of reservoir based on the data deck
            //if (Data.PhaseFlags.Count == 5)
            //{
            //    //reservoir = new WetGasLiveOilReservoir(Data);
            //}
            //else if (Data.PhaseFlags.Count == 4)
            //{
            //    //reservoir = Data.PhaseFlags.Contains("DISGAS") ?
            //    //        new DryGasLiveOilReservoir(Data) :
            //    //        new WetGasDeadOilReservoir(Data);
            //}
            //else if (Data.PhaseFlags.Count == 3)
            //{
            //    //reservoir = new DryGasDeadOilReservoir(Data);
            //}
            //else if (Data.PhaseFlags.Count == 2)
            //{
            //    //if (Data.PhaseFlags.Contains("OIL") && Data.PhaseFlags.Contains("GAS"))
            //    //    reservoir = new DryGasDeadOilReservoir(Data);
            //    //else if (Data.PhaseFlags.Contains("OIL") && Data.PhaseFlags.Contains("WATER"))
            //    //    reservoir = new OilWaterReservoir(Data);
            //    //else reservoir = new DryGasWaterReservoir(Data);
            //}
            //else if (Data.PhaseFlags.Count == 1)
            //{
            //   // reservoir = Data.PhaseFlags.Contains("OIL") ?
            //    //    new OilWaterReservoir(Data) : new DryGasWaterReservoir(Data);
            //}
            //else
            //{
            //    throw new NotSupportedException("Unsupported phase configuration in the Eclipse data deck.");
            //}
        }


    

        public void Initialize()=> _reservoir.Initialize();
        public void Simulate(double[] ResultTime = null)=> _reservoir.Simulate(ResultTime);
        public void ExportResults(string outputDir)=> _reservoir.ExportResults(outputDir);
    }
}
