using SepalSolver;
using System.Runtime.Intrinsics.Arm;
using Typography.OpenFont.Tables;
using static ReservoirSimulator.Math;

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
        public CellPhaseState PhaseState { get; set; } = CellPhaseState.TwoPhaseOilGas;
        public bool OilIsSaturated { get; set; } = false;
        public bool GasIsSaturated { get; set; } = false;
        public double NTG { get; set; }
        public double Porosity { get; set; } = 0;
        public double PoreVolumeMultiplier { get; set; } = 1.0;   // MULTPV
        public Face[] Faces { get; set; } = [.. Enumerable.Range(0, 6).Select(i => new Face() { Direction = (FlowDirection)i })];
        public double PoreVolume { get; private set; } = 0;
        public double InvPoreVolume { get; private set; } = 0;

        public SWOFAnalysis SWOF 
        { 
            get => swof;
            set 
            { 
                swof = value;
                Krw = Sw => swof.Result.EvaluateKrw(Sw);
                Krow = Sw => swof.Result.EvaluateKrow(Sw);
                Pcow = Sw => swof.Result.EvaluatePcow(Sw);
            }
        }
        public SGOFAnalysis SGOF 
        { 
            get => sgof; 
            set 
            { 
                sgof = value;
                Krg = Sg => sgof.Result.EvaluateKrg(Sg);
                Krog = Sg => sgof.Result.EvaluateKrog(Sg);
                Pcog = Sg => sgof.Result.EvaluatePcog(Sg);
            } 
        }
        public SWGFAnalysis SWGF 
        { 
            get => swgf;
            set
            { 
                swgf = value;
                Krw = Sw => swgf.Result.EvaluateKrw(Sw);
                Krg = Sw => swgf.Result.EvaluateKrg(Sw);
                Pcgw = Sg => swgf.Result.EvaluatePcgw(Sg);
            } 
        }
        public PVTWAnalysis PVTW 
        { 
            get => pvtw;  
            set 
            { 
                pvtw = value;
                Bw = Pw => pvtw.Result.EvaluateBw(Pw);
                μw = Pw => pvtw.Result.EvaluateMuw(Pw);
                γw = Pw => γw0 / pvtw.Result.EvaluateBw(Pw);
            }
        }
        public PVDOAnalysis PVDO 
        { 
            get => pvdo;
            set
            { 
                pvdo = value;
                Bo = Po => pvdo.Result.EvaluateBo(Po);
                μo = Po => pvdo.Result.EvaluateMuo(Po);
                γo = Po => γo0 / pvdo.Result.EvaluateBo(Po);
            }
        }
        public PVDGAnalysis PVDG 
        { 
            get => pvdg;
            set 
            { 
                pvdg = value;
                Bg = Po => pvdg.Result.EvaluateBg(Po);
                μg = Po => pvdg.Result.EvaluateMug(Po);
                γg = Po => pvdg.Result.γg0 / pvdg.Result.EvaluateBg(Po);
            } 
        }
        public PVTOAnalysis PVTO
        {
            get => pvto;
            set
            {
                pvto = value;
                BoSat = Rs => pvto.Result.EvaluateBo(Rs);
                RsSat = Po => pvto.Result.EvaluateRs(Po);
                Pb = Rs => pvto.Result.EvaluateBubblePointPressure(Rs);
            }
        }
        public PVTGAnalysis PVTG 
        {
            get => pvtg;
            set 
            { 
                pvtg = value;
                BgSat = Pg => pvtg.Result.EvaluateBg(Pg);
                RvSat = Pg => pvtg.Result.EvaluateRv(Pg);
                Pd = Rv => pvtg.Result.EvaluateDewPointPressure(Rv);
            } 
        }
        public RSVDAnalysis RSVD 
        {
            get => rsvd;
            set 
            { 
                rsvd = value; 
                Rsvd = h => rsvd.Result.EvaluateRs(h);
            } 
        }
        public RVVDAnalysis RVVD 
        { 
            get => rvvd;
            set 
            { 
                rvvd = value;
                Rvvd = h => rvvd.Result.EvaluateRv(h);
            } 
        }

        SWOFAnalysis swof;
        SGOFAnalysis sgof;
        SWGFAnalysis swgf;
        PVTWAnalysis pvtw;
        PVDOAnalysis pvdo;
        PVDGAnalysis pvdg;
        PVTOAnalysis pvto;
        PVTGAnalysis pvtg;
        RSVDAnalysis rsvd;
        RVVDAnalysis rvvd;


        public double Po_n, Pw_n, Pg_n, Sw_n, So_n, Sg_n, Rs_n, Rv_n, ErSo_Bo_n, ErSw_Bw_n, ErSg_Bg_n, WI, γo0, γw0, γg0;
        public double Po_np1, Pw_np1, Pg_np1, Sw_np1, So_np1, Sg_np1, Rs_np1, Rv_np1, Bo_np1, Bw_np1, Bg_np1;
        public ADiff<T> ErSo_Bo_np1, ErSw_Bw_np1, ErSg_Bg_np1;
        public Func<ADiff<T>, ADiff<T>> Krw, Krow, Krg, Krog, Pcow, Pcog, Pcgw, Bo, Bw, Bg, μo, μw, μg, γo, γw, γg, Er, RsSat, BoSat, Pb, RvSat, BgSat, Pd;
        public Func<double, double> Rsvd, Rvvd;
        
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
        public ADiff<T> γopr(ADiff<T> p, ADiff<T> rs)
        {
            var (bo, rsEff) = EvaluateOilState(p, rs);
            return (γo0 + 178.1076 * rsEff * γg0) / bo;
        }
        public ADiff<T> μopr(ADiff<T> p, ADiff<T> rs)
        {
            var (bo, rsEff) = EvaluateOilState(p, rs);
            return pvto.Result.EvaluateMuo(p, rsEff);
        }
        public ADiff<T> Bopr(ADiff<T> p, ADiff<T> rs)
        {
            var (bo, rsEff) = EvaluateOilState(p, rs);
            return bo;
        }
        public ADiff<T> γgpr(ADiff<T> p, ADiff<T> rv)
        {
            var (bg, rvEff) = EvaluateGasState(p, rv);
            return (178.1076 * γg0 + rvEff * γo0) / bg;
        }
        public ADiff<T> μgpr(ADiff<T> p, ADiff<T> rv)
        {
            var (bg, rvEff) = EvaluateGasState(p, rv);
            return pvtg.Result.EvaluateMug(p, rvEff);
        }
        public ADiff<T> Bgpr(ADiff<T> p, ADiff<T> rv)
        {
            var (bg, rvEff) = EvaluateGasState(p, rv);
            return bg;
        }

        (ADiff<T> Bg, ADiff<T> RvEff) EvaluateGasState(ADiff<T> P, ADiff<T> Rv)
        {
            double RvBase = pvtg.Result.RvBase;
            double BgBase = pvtg.Result.BgBase;
            double BgExponent = pvtg.Result.BgExponent;
            double UndersatBgSlope = pvtg.Result.UndersatBgSlope;

            if (P <= 0.0 || BgBase <= 0.0)
                return (0.0, 0.0);

            // Dry gas (RvBase <= 0)
            if (RvBase <= 0.0)
                return (BgBase * Pow(P, -BgExponent), 0.0);

            // Wet gas (PVTG): clamp Rv to RvSat(P) if P <= Pd(Rv) (i.e., RvVal >= RvSat(P))
            ADiff<T> rvSat = RvSat(P), bgSat = BgSat(P);
            if (Rv >= rvSat)
            {
                return (bgSat, rvSat);
            }
            else
            {
                ADiff<T> rvDeficit = rvSat - Rv; // > 0 when undersaturated
                ADiff<T> bgUndersat = bgSat * Exp(-UndersatBgSlope * rvDeficit);
                return (bgUndersat, Rv);
            }
        }

        (ADiff<T> Bo, ADiff<T> RsEff) EvaluateOilState(ADiff<T> P, ADiff<T> Rs)
        {
            double RsBase = pvto.Result.RsBase;
            double BoBase = pvto.Result.BoBase;
            double BoExponent = pvto.Result.BoExponent;
            double UndersatCompressibility = pvto.Result.UndersatCompressibility;

            // Dead oil (PVDO)
            if (RsBase <= 0.0)
            {
                ADiff<T> boDead = BoBase * Exp(-UndersatCompressibility * Max(0.0, P - 14.7));
                return (boDead, 0.0);
            }

            // Live oil (PVTO): clamp Rs to RsSat(P) if P <= Pb(Rs)
            ADiff<T> pb = Pb(Rs);
            if (P <= pb)
            {
                ADiff<T> rsSat = RsSat(P), boSat = BoSat(P);
                return (boSat, rsSat);
            }
            else
            {
                ADiff<T> boAtPb = BoBase * Exp(BoExponent * Rs);
                ADiff<T> boUndersat = boAtPb * Exp(-UndersatCompressibility * (P - pb));
                return (boUndersat, Rs);
            }
        }

        public int ActiveIndex { get; set; } = -1;
        public int GlobalIndex { get; set; } = -1;
        public (int i, int j, int k) IJKIndex { get; set; }
        internal Vector3d[] Corners { get; set; } = [];
        internal Vector3d Center { get; set; }
        public bool HasDrsdtLimit { get; internal set; }
        public string DrsdtOption { get; internal set; }

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