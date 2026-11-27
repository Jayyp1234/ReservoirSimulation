namespace EclipseDeckReader
{
    /// <summary>
    /// One AQUFETP record: a Fetkovich aquifer. Values are in the deck's own units; null means defaulted.
    /// </summary>
    public sealed class FetkovichAquiferRecord
    {
        public int Id { get; set; }
        public double DatumDepth { get; set; }
        public double? InitialPressure { get; set; }        // defaulted: the reservoir's pressure at the datum
        public double WaterVolume { get; set; }             // initial water volume of the aquifer (stb / sm3)
        public double TotalCompressibility { get; set; }    // rock + water (1/psi / 1/bar)
        public double ProductivityIndex { get; set; }       // (stb/d/psi / sm3/d/bar)
        public int PvtwTable { get; set; } = 1;
    }

    /// <summary>
    /// One AQUANCON record: the cells in a box (1-based, inclusive) that an aquifer touches through one face.
    /// </summary>
    public sealed class AquiferConnectionRecord
    {
        public int AquiferId { get; set; }
        public int I1 { get; set; }
        public int I2 { get; set; }
        public int J1 { get; set; }
        public int J2 { get; set; }
        public int K1 { get; set; }
        public int K2 { get; set; }
        public string Face { get; set; } = "I-";            // I- I+ J- J+ K- K+
        public double? InfluxCoefficient { get; set; }      // defaulted: the area of the connecting face
        public double InfluxMultiplier { get; set; } = 1.0;
        public bool AllowInternalFaces { get; set; }
    }
}
