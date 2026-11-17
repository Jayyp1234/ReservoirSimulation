using System;
using System.Collections.Generic;
using System.Text;

namespace ReservoirSimulator
{
    public class PressureDepthPair(double z, double p)
    {
        public double Z { get; set; } = z;
        public double P { get; set; } = p;
    }
    internal class EQUILPressureDepth(EQUILAnalysis equilResult)
    {
        public List<PressureDepthPair> GasZonePressureDepthPairs = [];
        public List<PressureDepthPair> OilZonePressureDepthPairs = [];
        public List<PressureDepthPair> WaterZonePressureDepthPairs = [];
        public double ZDatum { get; } = equilResult.ZDatum;
        public double PDatum { get; } = equilResult.PDatum;
        public double ZWOC { get; } = equilResult.zWOC;
        public double PcWOC { get; } = equilResult.pcWOC;
        public double ZGOC { get; } = equilResult.zGOC;
        public double PcGOC { get; } = equilResult.pcGOC;
        public double PoWOC { get; set; } = equilResult.PoWOC;
        public double PwWOC { get; set; } = equilResult.PwWOC;
        public double PoGOC { get; set; } = equilResult.PoGOC;
        public double PgGOC { get; set; } = equilResult.PgGOC;
    }
}
