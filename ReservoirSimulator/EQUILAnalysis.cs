using SepalSolver;
using static SepalSolver.Math;

namespace ReservoirSimulator
{   
    public class EQUILAnalysis
    {
        public double ZDatum { get; set; }
        public double PDatum { get; set; }
        public double zWOC { get; set; }
        public double pcWOC { get; set; }
        public double zGOC { get; set; }
        public double pcGOC { get; set; }
        public double PoWOC { get; set; }
        public double PoGOC { get; set; }
        public double PwWOC { get; set; }
        public double PgGOC { get; set; }

        public EQUILAnalysis(List<double> equil, Units incomingUnit)
        {
            ZDatum = Utility.Convert(equil[0], "Depth", incomingUnit, Units.FIELD);
            PDatum = Utility.Convert(equil[1], "Pressure", incomingUnit, Units.FIELD);
            zWOC = Utility.Convert(equil[2], "Depth", incomingUnit, Units.FIELD);
            pcWOC = equil.Count >= 4 ? Utility.Convert(equil[3], "Pressure", incomingUnit, Units.FIELD) : 0;
            zGOC = equil.Count >= 5 ? Utility.Convert(equil[4], "Depth", incomingUnit, Units.FIELD) : 0;
            pcGOC = equil.Count >= 6 ? Utility.Convert(equil[5], "Pressure", incomingUnit, Units.FIELD) : 0;
        }
    }
}
