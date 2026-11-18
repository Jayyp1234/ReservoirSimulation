using SepalSolver;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace ReservoirSimulator
{
    public static class VTUExporter
    {
        /// <summary>
        /// Writes a binary UnstructuredGrid (.vtu) snapshot for a faulted corner-point grid.
        /// Each active cell is emitted as a discrete 8-node VTK_HEXAHEDRON (type 12) so fault throws
        /// and pinched-out / inactive cells render accurately without vertex welding.
        /// </summary>
        public static void WriteVtuBinaryOilWater<T>(string filePath, Grid<T>[] Grids)
            where T : struct, IVectorND<T>
        {
            int numActiveGrids = Grids.Count();
            int numPoints = numActiveGrids * 8;

            // Pre-pack flat binary buffers for zero-overhead sequential disk writes
            double[] points = new double[numPoints * 3];
            int[] connectivity = new int[numActiveGrids * 8];
            int[] offsets = new int[numActiveGrids];
            byte[] types = new byte[numActiveGrids];

            double[] po = new double[numActiveGrids];
            double[] pw = new double[numActiveGrids];
            double[] sw = new double[numActiveGrids];
            double[] so = new double[numActiveGrids];
            double[] permX = new double[numActiveGrids];
            double[] permY = new double[numActiveGrids];
            double[] permZ = new double[numActiveGrids];
            double[] poro = new double[numActiveGrids];

            int cell = 0;
            foreach(var grid in Grids)
            {
                if (!grid.isActive) continue;

                var c = grid.Corners;
                int pBase = cell * 8;
                int coordBase = pBase * 3;

                for (int v = 0; v < 8; v++)
                {
                    points[coordBase + v * 3 + 0] = c[v].X;
                    points[coordBase + v * 3 + 1] = c[v].Y;
                    points[coordBase + v * 3 + 2] = -c[v].Z; // Negate depth-positive Z for right-handed 3D viewers
                    connectivity[pBase + v] = pBase + v;
                }

                offsets[cell] = (cell + 1) * 8;
                types[cell] = 12; // VTK_HEXAHEDRON

                po[cell] = grid.Po_n;
                pw[cell] = grid.Pw_n;
                sw[cell] = grid.Sw_n;
                so[cell] = grid.So_n;
                permX[cell] = grid.Faces[0].Permeability;
                permY[cell] = grid.Faces[2].Permeability;
                permZ[cell] = grid.Faces[4].Permeability;
                poro[cell] = grid.Porosity;
                cell++;
            }

            // Compute byte offsets in the <AppendedData> block (each block has a 4-byte UInt32 byte-length header)
            long offPoints = 0;
            long offConn = offPoints + 4 + (long)points.Length * sizeof(double);
            long offOffset = offConn + 4 + (long)connectivity.Length * sizeof(int);
            long offTypes = offOffset + 4 + (long)offsets.Length * sizeof(int);
            long offPo = offTypes + 4 + (long)types.Length * sizeof(byte);
            long offPw = offPo + 4 + (long)po.Length * sizeof(double);
            long offSw = offPw + 4 + (long)pw.Length * sizeof(double);
            long offSo = offSw + 4 + (long)sw.Length * sizeof(double);
            long offPermX = offSo + 4 + (long)so.Length * sizeof(double);
            long offPermY = offPermX + 4 + (long)permX.Length * sizeof(double);
            long offPermZ = offPermY + 4 + (long)permY.Length * sizeof(double);
            long offPoro = offPermZ + 4 + (long)permZ.Length * sizeof(double);

            using var fs = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.None, bufferSize: 262144);
            var utf8 = new UTF8Encoding(encoderShouldEmitUTF8Identifier: false);

            string header = $"""
                <?xml version="1.0"?>
                <VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">
                  <UnstructuredGrid>
                    <Piece NumberOfPoints="{numPoints}" NumberOfCells="{numActiveGrids}">
                      <Points>
                        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{offPoints}"/>
                      </Points>
                      <Cells>
                        <DataArray type="Int32" Name="connectivity" format="appended" offset="{offConn}"/>
                        <DataArray type="Int32" Name="offsets" format="appended" offset="{offOffset}"/>
                        <DataArray type="UInt8" Name="types" format="appended" offset="{offTypes}"/>
                      </Cells>
                      <CellData Scalars="Sw">
                        <DataArray type="Float64" Name="Po_psi" format="appended" offset="{offPo}"/>
                        <DataArray type="Float64" Name="Pw_psi" format="appended" offset="{offPw}"/>
                        <DataArray type="Float64" Name="Sw" format="appended" offset="{offSw}"/>
                        <DataArray type="Float64" Name="So" format="appended" offset="{offSo}"/>
                        <DataArray type="Float64" Name="PermX_md" format="appended" offset="{offPermX}"/>
                        <DataArray type="Float64" Name="PermY_md" format="appended" offset="{offPermY}"/>
                        <DataArray type="Float64" Name="PermZ_md" format="appended" offset="{offPermZ}"/>
                        <DataArray type="Float64" Name="Porosity" format="appended" offset="{offPoro}"/>
                      </CellData>
                    </Piece>
                  </UnstructuredGrid>
                  <AppendedData encoding="raw">
                _
                """;

            fs.Write(utf8.GetBytes(header));

            // Write raw binary blocks immediately after '_'
            WriteRawBlock<double>(fs, points);
            WriteRawBlock<int>(fs, connectivity);
            WriteRawBlock<int>(fs, offsets);
            WriteRawBlock<byte>(fs, types);
            WriteRawBlock<double>(fs, po);
            WriteRawBlock<double>(fs, pw);
            WriteRawBlock<double>(fs, sw);
            WriteRawBlock<double>(fs, so);
            WriteRawBlock<double>(fs, permX);
            WriteRawBlock<double>(fs, permY);
            WriteRawBlock<double>(fs, permZ);
            WriteRawBlock<double>(fs, poro);

            string footer = "\n  </AppendedData>\n</VTKFile>\n";
            fs.Write(utf8.GetBytes(footer));
        }

        public static void WriteVtuBinaryBlackOil<T>(string filePath, Grid<T>[] Grids)
            where T : struct, IVectorND<T>
        {
            int numActiveGrids = Grids.Count();
            int numPoints = numActiveGrids * 8;

            // Pre-pack flat binary buffers for zero-overhead sequential disk writes
            double[] points = new double[numPoints * 3];
            int[] connectivity = new int[numActiveGrids * 8];
            int[] offsets = new int[numActiveGrids];
            byte[] types = new byte[numActiveGrids];

            double[] po = new double[numActiveGrids];
            double[] pw = new double[numActiveGrids];
            double[] pg = new double[numActiveGrids];
            double[] pb = new double[numActiveGrids];
            double[] sw = new double[numActiveGrids];
            double[] so = new double[numActiveGrids];
            double[] sg = new double[numActiveGrids];
            double[] permX = new double[numActiveGrids];
            double[] permY = new double[numActiveGrids];
            double[] permZ = new double[numActiveGrids];
            double[] poro = new double[numActiveGrids];

            int cell = 0;
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;

                var c = grid.Corners;
                int pBase = cell * 8;
                int coordBase = pBase * 3;

                for (int v = 0; v < 8; v++)
                {
                    points[coordBase + v * 3 + 0] = c[v].X;
                    points[coordBase + v * 3 + 1] = c[v].Y;
                    points[coordBase + v * 3 + 2] = -c[v].Z; // Negate depth-positive Z for right-handed 3D viewers
                    connectivity[pBase + v] = pBase + v;
                }

                offsets[cell] = (cell + 1) * 8;
                types[cell] = 12; // VTK_HEXAHEDRON

                po[cell] = grid.Po_n;
                pw[cell] = grid.Pw_n;
                pg[cell] = grid.Pg_n;
                pb[cell] = grid.Pb_n;
                sw[cell] = grid.Sw_n;
                so[cell] = grid.So_n;
                sg[cell] = grid.Sg_n;
                sw[cell] = grid.Sw_n;
                so[cell] = grid.So_n;
                permX[cell] = grid.Faces[0].Permeability;
                permY[cell] = grid.Faces[2].Permeability;
                permZ[cell] = grid.Faces[4].Permeability;
                poro[cell] = grid.Porosity;
                cell++;
            }

            // Compute byte offsets in the <AppendedData> block (each block has a 4-byte UInt32 byte-length header)
            long offPoints = 0;
            long offConn = offPoints + 4 + (long)points.Length * sizeof(double);
            long offOffset = offConn + 4 + (long)connectivity.Length * sizeof(int);
            long offTypes = offOffset + 4 + (long)offsets.Length * sizeof(int);
            long offPo = offTypes + 4 + (long)types.Length * sizeof(byte);
            long offPw = offPo + 4 + (long)po.Length * sizeof(double);
            long offSw = offPw + 4 + (long)pw.Length * sizeof(double);
            long offSo = offSw + 4 + (long)sw.Length * sizeof(double);
            long offPermX = offSo + 4 + (long)so.Length * sizeof(double);
            long offPermY = offPermX + 4 + (long)permX.Length * sizeof(double);
            long offPermZ = offPermY + 4 + (long)permY.Length * sizeof(double);
            long offPoro = offPermZ + 4 + (long)permZ.Length * sizeof(double);
            long offPg = offPoro + 4 + (long)poro.Length * sizeof(double);
            long offPb = offPg + 4 + (long)pg.Length * sizeof(double);
            long offSg = offPb + 4 + (long)pb.Length * sizeof(double);

            using var fs = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.None, bufferSize: 262144);
            var utf8 = new UTF8Encoding(encoderShouldEmitUTF8Identifier: false);

            string header = $"""
                <?xml version="1.0"?>
                <VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">
                  <UnstructuredGrid>
                    <Piece NumberOfPoints="{numPoints}" NumberOfCells="{numActiveGrids}">
                      <Points>
                        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{offPoints}"/>
                      </Points>
                      <Cells>
                        <DataArray type="Int32" Name="connectivity" format="appended" offset="{offConn}"/>
                        <DataArray type="Int32" Name="offsets" format="appended" offset="{offOffset}"/>
                        <DataArray type="UInt8" Name="types" format="appended" offset="{offTypes}"/>
                      </Cells>
                      <CellData Scalars="Sw">
                        <DataArray type="Float64" Name="Po_psi" format="appended" offset="{offPo}"/>
                        <DataArray type="Float64" Name="Pw_psi" format="appended" offset="{offPw}"/>
                        <DataArray type="Float64" Name="Sw" format="appended" offset="{offSw}"/>
                        <DataArray type="Float64" Name="So" format="appended" offset="{offSo}"/>
                        <DataArray type="Float64" Name="PermX_md" format="appended" offset="{offPermX}"/>
                        <DataArray type="Float64" Name="PermY_md" format="appended" offset="{offPermY}"/>
                        <DataArray type="Float64" Name="PermZ_md" format="appended" offset="{offPermZ}"/>
                        <DataArray type="Float64" Name="Porosity" format="appended" offset="{offPoro}"/>
                        <DataArray type="Float64" Name="Pg_psi" format="appended" offset="{offPg}"/>
                        <DataArray type="Float64" Name="Pb_psi" format="appended" offset="{offPb}"/>
                        <DataArray type="Float64" Name="Sg" format="appended" offset="{offSg}"/>
                      </CellData>
                    </Piece>
                  </UnstructuredGrid>
                  <AppendedData encoding="raw">
                _
                """;

            fs.Write(utf8.GetBytes(header));

            // Write raw binary blocks immediately after '_'
            WriteRawBlock<double>(fs, points);
            WriteRawBlock<int>(fs, connectivity);
            WriteRawBlock<int>(fs, offsets);
            WriteRawBlock<byte>(fs, types);
            WriteRawBlock<double>(fs, po);
            WriteRawBlock<double>(fs, pw);
            WriteRawBlock<double>(fs, pg);
            WriteRawBlock<double>(fs, pb);
            WriteRawBlock<double>(fs, sw);
            WriteRawBlock<double>(fs, so);
            WriteRawBlock<double>(fs, sg);
            WriteRawBlock<double>(fs, permX);
            WriteRawBlock<double>(fs, permY);
            WriteRawBlock<double>(fs, permZ);
            WriteRawBlock<double>(fs, poro);

            string footer = "\n  </AppendedData>\n</VTKFile>\n";
            fs.Write(utf8.GetBytes(footer));
        }
        
        /// <summary>
        /// Writes a ParaView Data (.pvd) collection indexing all .vtu snapshots by simulation time (days).
        /// </summary>
        public static void WritePvdCollection(string pvdPath, IReadOnlyList<(double TimeDays, string VtuFileName)> snapshots)
        {
            using var sw = new StreamWriter(pvdPath, append: false, new UTF8Encoding(false));
            sw.WriteLine("""<?xml version="1.0"?>""");
            sw.WriteLine("""<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""");
            sw.WriteLine("""  <Collection>""");
            foreach (var (time, file) in snapshots)
            {
                sw.WriteLine(string.Create(CultureInfo.InvariantCulture,
                    $"""    <DataSet timestep="{time:G9}" group="" part="0" file="{file}"/>"""));
            }
            sw.WriteLine("""  </Collection>""");
            sw.WriteLine("""</VTKFile>""");
        }

        private static void WriteRawBlock<T>(FileStream fs, T[] data) where T : unmanaged
        {
            ReadOnlySpan<byte> bytes = MemoryMarshal.AsBytes(data.AsSpan());
            uint byteLength = (uint)bytes.Length;
            Span<byte> lenHeader = stackalloc byte[4];
            BitConverter.TryWriteBytes(lenHeader, byteLength);
            fs.Write(lenHeader);
            fs.Write(bytes);
        }
    }
}