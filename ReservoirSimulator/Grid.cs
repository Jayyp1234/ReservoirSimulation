using System.Numerics;
using SepalSolver;

namespace ReservoirSimulator
{
    public class Grid<T> where T : struct, IVectorND<T>
    {
        public Offsets Offset { get; set; }
        public string PerforationStatus { get; set; } = "SHUT";
        // Define conversion constants
        const double alpha = 1.127e-3, beta = 5.615;
        public bool LogResult { get; set; }
        public bool isActive { get { return ActiveIndex >= 0; } }
        public bool GasIsDIssolved { get; set; } = false;
        public bool OilIsVaporised { get; set; } = false;
        public double NTG { get; set; }
        public double Porosity { get; set; } = 0;
        public double PoreVolumeMultiplier { get; set; } = 1.0;   // MULTPV
        public Face[] Faces { get; set; } = [.. Enumerable.Range(0, 6).Select(i => new Face() { Direction = (FlowDirection)i })];
        public double PoreVolume { get; private set; } = 0;
        public double InvPoreVolume { get; private set; } = 0;

        public double Po_n, Pw_n, Pg_n, Sw_n, So_n, Sg_n, Pb_n, Pd_n, ErSo_Bo_n, ErSw_Bw_n, ErSg_Bg_n, WI;
        public double Po_np1, Pw_np1, Pg_np1, Sw_np1, So_np1, Sg_np1, Pb_np1, Pd_np1, Bo_np1, Bw_np1, Bg_np1;
        public ADiff<T> ErSo_Bo_np1, ErSw_Bw_np1, ErSg_Bg_np1;
        public Func<ADiff<T>, ADiff<T>> Krw, Krow, Krg, Krog, Pcow, Pcog, Pcgw, Bo, Bw, Bg, μo, μw, μg, γo, γw, γg, Er, Rs, Rv;
        
        // Baker's 3-phase relative permeability (smooth everywhere on the Sw-Sg triangle)
        public ADiff<T> Kro(ADiff<T> sw, ADiff<T> sg)
        {
            var krow = Krow(sw);
            var krog = Krog(sg);
            var swEff = sw;
            var denom = swEff + sg;

            // When Sg == 0 and Sw == Swi, Kro = Krow(Swi) = 1.0
            if (denom.Value <= 1e-12)
                return krow;

            return (swEff * krow + sg * krog) / denom;
        }
        public int ActiveIndex { get; set; } = -1;
        public int GlobalIndex { get; set; } = -1;
        public (int i, int j, int k) IJKIndex { get; set; }
        internal Vector3d[] Corners { get; set; } = [];
        internal Vector3d Center { get; set; }
        internal bool WriteBPR { get; set; } = false;
        internal bool WriteBPW { get; set; } = false;
        internal bool WriteBPG { get; set; } = false;
        internal bool WriteBPB { get; set; } = false;
        internal bool WriteBPCOW { get; set; } = false;
        internal bool WriteBPCOG { get; set; } = false;
        internal bool WriteBPOTO { get; set; } = false;
        internal bool WriteBPOTW { get; set; } = false;
        internal bool WriteBPOTG { get; set; } = false;
        internal bool WriteBWSAT { get; set; } = false;
        internal bool WriteBGSAT { get; set; } = false;
        internal bool WriteBOSAT { get; set; } = false;
        internal bool WriteBRS { get; set; } = false;
        internal bool WriteBRV { get; set; } = false;
        internal bool BBO { get; set; } = false;
        internal bool BBW { get; set; } = false;
        internal bool BBG { get; set; } = false;
        internal bool BVISO { get; set; } = false;
        internal bool BVISW { get; set; } = false;
        internal bool BVISG { get; set; } = false;


        public List<double> PO = [], PW = [], PG = [], SW = [], SO = [], SG = [], PB = [], PD = [];

        internal void MakeFaces()
        {
            Center = new();
            for (int i = 0; i < 8; i++) Center += Corners[i];
            Center /= 8.0;

            Faces[0].Corners = [Corners[0], Corners[4], Corners[7], Corners[3]];  // Left (-I)
            Faces[1].Corners = [Corners[1], Corners[2], Corners[6], Corners[5]];  // Right (+I)
            Faces[2].Corners = [Corners[0], Corners[1], Corners[5], Corners[4]];  // Front (-J)
            Faces[3].Corners = [Corners[3], Corners[7], Corners[6], Corners[2]];  // Back (+J)
            Faces[4].Corners = [Corners[0], Corners[3], Corners[2], Corners[1]];  // Top (-K)
            Faces[5].Corners = [Corners[4], Corners[5], Corners[6], Corners[7]];  // Bottom (+K)

            PoreVolume = Utility.ComputeHexahedronVolume(this) * NTG * Porosity * PoreVolumeMultiplier / beta;
            InvPoreVolume = 1.0 / PoreVolume;

            foreach (var face in Faces)
            {
                face.ComputeArea();
                face.GlobalGridIndex = GlobalIndex;
                face.Center = new();
                foreach (var corner in face.Corners) face.Center += corner;
                face.Center /= 4.0;
                face.PerpendicularDistance2Center =
                    Utility.ComputePerpendicularDistance(Center, face.Corners[0], face.Corners[1], face.Corners[2], face.Corners[3]);
            }
        }

        public void Reset()
        {
            PO.Clear(); PW.Clear(); PG.Clear(); SW.Clear(); SO.Clear(); SG.Clear(); PB.Clear(); PD.Clear();
        }

        public void SaveResult(
            double? po = null, double? pw = null, double? pg = null,
            double? sw = null, double? so = null, double? sg = null,
            double? pb = null, double? pd = null)
        {
            Utility.AddIfPresent(PO, po);
            Utility.AddIfPresent(PW, pw);
            Utility.AddIfPresent(PG, pg);
            Utility.AddIfPresent(SW, sw);
            Utility.AddIfPresent(SO, so);
            Utility.AddIfPresent(SG, sg);
            Utility.AddIfPresent(PB, pb);
            Utility.AddIfPresent(PD, pd);
        }

        public readonly struct Offsets
        {
            // Flattened indices into values[] for diagonal block (i, i) of size nc x nc
            // Index mapping: [eq * nc + var]
            public readonly int[] Diag;

            public Offsets(int nc) => Diag = new int[nc * nc];
        }
    }
}