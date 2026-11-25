using SepalSolver;
using System;
using System.Collections.Generic;
using System.Text;

namespace ReservoirSimulator
{
    internal interface IReservoir<T> where T : struct, IVectorND<T>
    {
        public List<Well<T>> WELLS { get; }
        public Grid<T>[] GRIDS { get; }
        public void Initialize();
        public void Simulate(double[] ResultTime = null);
        public void ExportResults(string outputDir);
    }
}
