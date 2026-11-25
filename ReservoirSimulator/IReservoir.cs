using SepalSolver;
using System;
using System.Collections.Generic;
using System.Text;

namespace ReservoirSimulator
{
    public interface IReservoir
    {
        bool UseUpstreamWeightingForInjection { get; set; }
        void Initialize();
        void Simulate(double[] ResultTime = null);
        void ExportResults(string outputDir);
    }

    public interface IReservoir<T> : IReservoir where T : struct, IVectorND<T>
    {
    }
}
