using SepalSolver;
using static SepalSolver.Math;
using System.Runtime.CompilerServices;

namespace ReservoirSimulator
{
    public static class ReservoirJacobian
    {
        public static (int totalDofs, int[] col_index, int[] row_ptr, double[] values) Create<T>
            (int nc, IReadOnlyList<Grid<T>> Grids, IReadOnlyList<Well<T>> Wells) where T : struct, IVectorND<T>
        {
            int numWells = Wells.Count;
            int numActiveGrids = Grids.Count(grid => grid.isActive);
            int totalDofs = nc * numActiveGrids + 2 * numWells;

            // 1. Exact capacity estimation
            int totalGridBlocks = numActiveGrids; // Diagonal blocks
            foreach (var mgrid in Grids)
            {
                if (!mgrid.isActive) continue;
                foreach (var face in mgrid.Faces)
                {
                    foreach (int c in face.GridIndices)
                    {
                        if (Grids[c].isActive) totalGridBlocks++;
                    }
                }
            }

            int totalActivePerfs = 0;
            foreach (var well in Wells)
            {
                foreach (int c in well.PerfNatIndex)
                {
                    if (Grids[c].isActive) totalActivePerfs++;
                }
            }

            int capacity = (totalGridBlocks * nc * nc) + (numWells * 4)
                         + (totalActivePerfs * 4 * nc);

            List<int> I = new(capacity), J = new(capacity);

            void AddBlock(int i, int j)
            {
                int I0 = i * nc, J0 = j * nc;
                for (int ii = 0; ii < nc; ii++)
                {
                    for (int jj = 0; jj < nc; jj++)
                    { I.Add(I0 + ii); J.Add(J0 + jj); }
                }
            }

            void AddGrid_wellBlock(int gridActiveIdx, int WellIdx)
            {
                int I0 = gridActiveIdx * nc,  J0 = nc * numActiveGrids + 2 * WellIdx;
                for (int ii = 0; ii < nc; ii++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    { I.Add(I0 + ii); J.Add(J0 + jj); }
                }
            }

            void AddWell_gridBlock(int WellIdx, int gridActiveIdx)
            {
                int I0 = nc * numActiveGrids + 2 * WellIdx, J0 = gridActiveIdx * nc;
                for (int ii = 0; ii < 2; ii++)
                {
                    for (int jj = 0; jj < nc; jj++)
                    { I.Add(I0 + ii); J.Add(J0 + jj); }
                }
            }

            void AddWellInternalBlock(int WellIdx)
            {
                int I0 = nc * numActiveGrids + 2 * WellIdx;
                for (int ii = 0; ii < 2; ii++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    { I.Add(I0 + ii); J.Add(I0 + jj); }
                }
            }

            // 3. Grid Cells and Face Neighbor Connections
            foreach (var mgrid in Grids)
            {
                if (!mgrid.isActive) continue;
                int i = mgrid.ActiveIndex;

                AddBlock(i, i);

                foreach (var face in mgrid.Faces)
                {
                    foreach (int c in face.GridIndices)
                    {
                        var ngrid = Grids[c];
                        if (!ngrid.isActive) continue;
                        AddBlock(i, ngrid.ActiveIndex);
                    }
                }
            }

            // 4. Active Wells Internal Closure and Perforations
            int wellIdx = 0;
            foreach (var well in Wells)
            {
                AddWellInternalBlock(wellIdx);
                foreach (int m in well.PerfNatIndex)
                {
                    var mgrid = Grids[m];
                    if (!mgrid.isActive) continue;
                    int i = mgrid.ActiveIndex;

                    AddGrid_wellBlock(i, wellIdx);
                    AddWell_gridBlock(wellIdx, i);
                }
                wellIdx++;
            }

            double[] ones = Ones(I.Count);
            var jac = new SparseMatrix(I, J, ones, totalDofs, totalDofs);
            var crs = jac.CompressByRow();
            return (totalDofs, [.. crs.col_index], [.. crs.row_start], [.. crs.values]);
        }


        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static int FindValueIndex(int[] row_ptr, int[] col_index, int row, int col)
        {
            int start = row_ptr[row];
            int length = row_ptr[row + 1] - start;
            int idx = Array.BinarySearch(col_index, start, length, col);
            return idx;
        }

        public static void PrecomputeOffsets<T>
            (int nc, int[] row_ptr, int[] col_index, IReadOnlyList<Grid<T>> Grids, IReadOnlyList<Well<T>> Wells) where T : struct, IVectorND<T>
        {
            int numWells = Wells.Count;
            int numActiveGrids = Grids.Count(grid => grid.isActive);

            // 1. Grid Diagonal & Face Offsets
            for (int m = 0; m < Grids.Count; m++)
            {
                var mgrid = Grids[m];
                if (!mgrid.isActive) continue;

                int i = mgrid.ActiveIndex;
                int baseRow_i = nc * i, baseCol_i = nc * i;

                // Diagonal Block (i, i)
                var cOff = new Grid<T>.Offsets(nc);
                for (int eq = 0; eq < nc; eq++)
                {
                    for (int v = 0; v < nc; v++)
                    {
                        cOff.Diag[eq * nc + v] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, baseCol_i + v);
                    }
                }
                mgrid.Offset = cOff;

                
                // Faces
                int faceCount = mgrid.Faces.Length;
                for (int f = 0; f < faceCount; f++)
                {
                    var face = mgrid.Faces[f];
                    int connCount = face.GridIndices.Count;

                    for (int c = 0; c < connCount; c++)
                    {
                        int n = face.GridIndices[c];
                        var ngrid = Grids[n];
                        if (!ngrid.isActive) continue;

                        int j = ngrid.ActiveIndex;
                        int baseCol_j = nc * j;

                        var fOff = new Face.Offsets(nc, j);
                        for (int eq = 0; eq < nc; eq++)
                        {
                            for (int v = 0; v < nc; v++)
                            {
                                fOff.Neighbor[eq * nc + v] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, baseCol_j + v);
                            }
                        }
                        face.Offset.Add(fOff);
                    }
                }
            }

            int wIdx = 0;

            // 2. Active Well Offsets
            foreach (var well in Wells)
            {
                int rowflow = nc * numActiveGrids + 2 * wIdx + 0;
                int rowconst = nc * numActiveGrids + 2 * wIdx + 1;
                int col_pwell = rowflow;
                int col_qwell = rowconst;

                int rflowP = FindValueIndex(row_ptr, col_index, rowflow, col_pwell);
                int rflowQ = FindValueIndex(row_ptr, col_index, rowflow, col_qwell);
                int rconstP = FindValueIndex(row_ptr, col_index, rowconst, col_pwell);
                int rconstQ = FindValueIndex(row_ptr, col_index, rowconst, col_qwell);

                int perfCount = well.PerfNatIndex.Count;
                var rmDPwell = new int[perfCount][];
                var rflowDGrid = new int[perfCount][];

                for (int p = 0; p < perfCount; p++)
                {
                    int m = well.PerfNatIndex[p];
                    var mgrid = Grids[m];
                    if (!mgrid.isActive) continue;

                    int i = mgrid.ActiveIndex;
                    int baseRow_i = nc * i;
                    int baseCol_i = nc * i;

                    rmDPwell[p] = new int[nc];
                    rflowDGrid[p] = new int[nc];

                    for (int eq = 0; eq < nc; eq++)
                    {
                        rmDPwell[p][eq] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, col_pwell);
                        rflowDGrid[p][eq] = FindValueIndex(row_ptr, col_index, rowflow, baseCol_i + eq);
                    }
                }
                well.Offset = new Well<T>.Offsets(rflowP, rflowQ, rconstP, rconstQ, rmDPwell, rflowDGrid);
                wIdx++;
            }
        }
    }
}