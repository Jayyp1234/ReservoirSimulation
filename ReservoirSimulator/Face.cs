using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Face
    {
        public List<Offsets> Offset { get; set; } = [];
        public int GlobalGridIndex { get; set; }
        public List<int> GridIndices { get; set; } = [];
        public List<double> Transmissibility { get; set; } = [];
        public double Permeability { get; set; }
        public FlowDirection Direction { get; set; }
        public Vector3d[] Corners { get; set; }

        public Vector3d Center { get; set; }

        public double PerpendicularDistance2Center { get; set; }
        public double Area { get; private set; }

        public void ComputeArea()
        {
            // Exact area of a general 4-corner quad via half the cross product of its diagonals
            Vector3d diagA = Corners[2] - Corners[0];
            Vector3d diagB = Corners[3] - Corners[1];
            Vector3d cross = Vector3d.Cross(diagA, diagB);

            Area = 0.5 * Sqrt(Vector3d.Dot(cross, cross));
        }

        public readonly struct Offsets
        {
            public readonly int NeighborActiveIndex;
            // Flattened indices into values[]:
            // Block (i, i) sensitivities (slots 0 .. nc-1): [eq * nc + var]
            public readonly int[] Self;
            // Block (i, j) sensitivities (slots nc .. 2*nc-1): [eq * nc + var]
            public readonly int[] Neighbor;

            public Offsets(int nc, int neighborActiveIdx)
            {
                NeighborActiveIndex = neighborActiveIdx;
                Neighbor = new int[nc * nc];
            }
        }
    }
}
