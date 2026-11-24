using SepalSolver;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace ReservoirSimulator
{
    public class SimulationAnimatorExporter<T> : IDisposable
        where T : struct, IVectorND<T>
    {
        private readonly string _outputDirectory;
        private readonly bool _isBlackOil;
        private readonly int _numActiveGrids;
        private readonly int _numPoints;
        private readonly int[] _activeIndices;

        // Cached static mesh buffers (computed once at t = 0)
        private readonly double[] _points;
        private readonly int[] _connectivity;
        private readonly int[] _offsets;
        private readonly byte[] _types;

        // Reusable dynamic property buffers (zero per-step allocations)
        private readonly double[] _po;
        private readonly double[] _pw;
        private readonly double[] _sw;
        private readonly double[] _so;
        private readonly double[]? _pg;
        private readonly double[]? _pb;
        private readonly double[]? _sg;

        // Precomputed XML header & byte offsets for dynamic snapshots
        private readonly byte[] _dynamicHeaderBytes;
        private readonly byte[] _footerBytes;

        private readonly List<(double TimeDays, string FileName)> _snapshots = new();

        public SimulationAnimatorExporter(string outputDirectory, Grid<T>[] grids, bool isBlackOil = false)
        {
            _outputDirectory = outputDirectory;
            _isBlackOil = isBlackOil;
            Directory.CreateDirectory(_outputDirectory);

            // 1. Collect active cell indices once
            _activeIndices = new int[grids.Count(g => g.isActive)];
            int idx = 0;
            for (int m = 0; m < grids.Length; m++)
            {
                if (grids[m].isActive)
                    _activeIndices[idx++] = m;
            }

            _numActiveGrids = _activeIndices.Length;
            _numPoints = _numActiveGrids * 8;

            // 2. Build static corner-point geometry once
            _points = new double[_numPoints * 3];
            _connectivity = new int[_numActiveGrids * 8];
            _offsets = new int[_numActiveGrids];
            _types = new byte[_numActiveGrids];

            double[] permX = new double[_numActiveGrids];
            double[] permY = new double[_numActiveGrids];
            double[] permZ = new double[_numActiveGrids];
            double[] poro = new double[_numActiveGrids];
            double[] ntg = new double[_numActiveGrids];

            for (int cell = 0; cell < _numActiveGrids; cell++)
            {
                ref readonly var grid = ref grids[_activeIndices[cell]];
                var c = grid.Corners;
                int pBase = cell * 8;
                int coordBase = pBase * 3;

                for (int v = 0; v < 8; v++)
                {
                    _points[coordBase + v * 3 + 0] = c[v].X;
                    _points[coordBase + v * 3 + 1] = c[v].Y;
                    _points[coordBase + v * 3 + 2] = -c[v].Z; // Depth-positive -> right-handed Z-up
                    _connectivity[pBase + v] = pBase + v;
                }

                _offsets[cell] = (cell + 1) * 8;
                _types[cell] = 12; // VTK_HEXAHEDRON

                permX[cell] = grid.Faces[0].Permeability;
                permY[cell] = grid.Faces[2].Permeability;
                permZ[cell] = grid.Faces[4].Permeability;
                poro[cell] = grid.Porosity;
                ntg[cell] = grid.NTG;
            }

            // 3. Write STATIC.vtu ONCE with geometry + static rock properties
            WriteStaticVtu(Path.Combine(_outputDirectory, "STATIC.vtu"), permX, permY, permZ, poro, ntg);

            // 4. Pre-allocate dynamic buffers & pre-build dynamic XML header once
            _po = new double[_numActiveGrids];
            _pw = new double[_numActiveGrids];
            _sw = new double[_numActiveGrids];
            _so = new double[_numActiveGrids];

            if (_isBlackOil)
            {
                _pg = new double[_numActiveGrids];
                _pb = new double[_numActiveGrids];
                _sg = new double[_numActiveGrids];
            }

            var utf8 = new UTF8Encoding(encoderShouldEmitUTF8Identifier: false);
            _dynamicHeaderBytes = utf8.GetBytes(BuildDynamicXmlHeader());
            _footerBytes = utf8.GetBytes("\n  </AppendedData>\n</VTKFile>\n");
        }

        /// <summary>
        /// Writes a dynamic snapshot (Pressure & Saturation only) for the given simulation time
        /// and updates SIMULATION.pvd for immediate animation playback in ParaView.
        /// </summary>
        public void WriteTimeStep(List<double> timeDays, Grid<T>[] grids, int n)
        {
            for (int cell = 0; cell < _numActiveGrids; cell++)
            {
                ref readonly var grid = ref grids[_activeIndices[cell]];
                _po[cell] = grid.PO[n];
                _pw[cell] = grid.PW[n];
                _sw[cell] = grid.SW[n];
                _so[cell] = grid.SO[n];

                if (_isBlackOil)
                {
                    _pg![cell] = grid.PG[n];
                    // Bubble point is not saved by the dead-oil model (PB stays empty); write 0 instead of throwing
                    _pb![cell] = n < grid.PB.Count ? grid.PB[n] : 0.0;
                    _sg![cell] = grid.SG[n];
                }
            }

            string fileName = $"STEP_{_snapshots.Count:D4}.vtu";
            string fullPath = Path.Combine(_outputDirectory, fileName);

            using (var fs = new FileStream(fullPath, FileMode.Create, FileAccess.Write, FileShare.None, bufferSize: 262144))
            {
                fs.Write(_dynamicHeaderBytes);

                // Cached static geometry blocks
                WriteRawBlock<double>(fs, _points);
                WriteRawBlock<int>(fs, _connectivity);
                WriteRawBlock<int>(fs, _offsets);
                WriteRawBlock<byte>(fs, _types);

                // Dynamic state blocks only (order matches BuildDynamicXmlHeader offsets)
                WriteRawBlock<double>(fs, _po);
                WriteRawBlock<double>(fs, _pw);
                WriteRawBlock<double>(fs, _sw);
                WriteRawBlock<double>(fs, _so);

                if (_isBlackOil)
                {
                    WriteRawBlock<double>(fs, _pg!);
                    WriteRawBlock<double>(fs, _pb!);
                    WriteRawBlock<double>(fs, _sg!);
                }

                fs.Write(_footerBytes);
            }

            _snapshots.Add((timeDays[n], fileName));
            WritePvdCollection();
        }

        private void WriteStaticVtu(string filePath, double[] permX, double[] permY, double[] permZ, double[] poro, double[] ntg)
        {
            long offPoints = 0;
            long offConn = offPoints + 4 + (long)_points.Length * sizeof(double);
            long offOffset = offConn + 4 + (long)_connectivity.Length * sizeof(int);
            long offTypes = offOffset + 4 + (long)_offsets.Length * sizeof(int);
            long offPermX = offTypes + 4 + (long)_types.Length * sizeof(byte);
            long offPermY = offPermX + 4 + (long)permX.Length * sizeof(double);
            long offPermZ = offPermY + 4 + (long)permY.Length * sizeof(double);
            long offPoro = offPermZ + 4 + (long)permZ.Length * sizeof(double);
            long offNtg = offPoro + 4 + (long)poro.Length * sizeof(double);

            var utf8 = new UTF8Encoding(encoderShouldEmitUTF8Identifier: false);
            string header = $"""
                <?xml version="1.0"?>
                <VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">
                  <UnstructuredGrid>
                    <Piece NumberOfPoints="{_numPoints}" NumberOfCells="{_numActiveGrids}">
                      <Points>
                        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{offPoints}"/>
                      </Points>
                      <Cells>
                        <DataArray type="Int32" Name="connectivity" format="appended" offset="{offConn}"/>
                        <DataArray type="Int32" Name="offsets" format="appended" offset="{offOffset}"/>
                        <DataArray type="UInt8" Name="types" format="appended" offset="{offTypes}"/>
                      </Cells>
                      <CellData Scalars="PermX_md">
                        <DataArray type="Float64" Name="PermX_md" format="appended" offset="{offPermX}"/>
                        <DataArray type="Float64" Name="PermY_md" format="appended" offset="{offPermY}"/>
                        <DataArray type="Float64" Name="PermZ_md" format="appended" offset="{offPermZ}"/>
                        <DataArray type="Float64" Name="Porosity" format="appended" offset="{offPoro}"/>
                        <DataArray type="Float64" Name="NTG" format="appended" offset="{offNtg}"/>
                      </CellData>
                    </Piece>
                  </UnstructuredGrid>
                  <AppendedData encoding="raw">
                _
                """;

            using var fs = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.None, bufferSize: 262144);
            fs.Write(utf8.GetBytes(header));

            WriteRawBlock<double>(fs, _points);
            WriteRawBlock<int>(fs, _connectivity);
            WriteRawBlock<int>(fs, _offsets);
            WriteRawBlock<byte>(fs, _types);
            WriteRawBlock<double>(fs, permX);
            WriteRawBlock<double>(fs, permY);
            WriteRawBlock<double>(fs, permZ);
            WriteRawBlock<double>(fs, poro);
            WriteRawBlock<double>(fs, ntg);

            fs.Write(utf8.GetBytes("\n  </AppendedData>\n</VTKFile>\n"));
        }

        private string BuildDynamicXmlHeader()
        {
            long scalarBytes = (long)_numActiveGrids * sizeof(double);

            long offPoints = 0;
            long offConn = offPoints + 4 + (long)_points.Length * sizeof(double);
            long offOffset = offConn + 4 + (long)_connectivity.Length * sizeof(int);
            long offTypes = offOffset + 4 + (long)_offsets.Length * sizeof(int);
            long offPo = offTypes + 4 + (long)_types.Length * sizeof(byte);
            long offPw = offPo + 4 + scalarBytes;
            long offSw = offPw + 4 + scalarBytes;
            long offSo = offSw + 4 + scalarBytes;

            var sb = new StringBuilder();
            sb.AppendLine("""<?xml version="1.0"?>""");
            sb.AppendLine("""<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">""");
            sb.AppendLine("""  <UnstructuredGrid>""");
            sb.AppendLine($"""    <Piece NumberOfPoints="{_numPoints}" NumberOfCells="{_numActiveGrids}">""");
            sb.AppendLine("""      <Points>""");
            sb.AppendLine($"""        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{offPoints}"/>""");
            sb.AppendLine("""      </Points>""");
            sb.AppendLine("""      <Cells>""");
            sb.AppendLine($"""        <DataArray type="Int32" Name="connectivity" format="appended" offset="{offConn}"/>""");
            sb.AppendLine($"""        <DataArray type="Int32" Name="offsets" format="appended" offset="{offOffset}"/>""");
            sb.AppendLine($"""        <DataArray type="UInt8" Name="types" format="appended" offset="{offTypes}"/>""");
            sb.AppendLine("""      </Cells>""");
            sb.AppendLine("""      <CellData Scalars="Sw">""");
            sb.AppendLine($"""        <DataArray type="Float64" Name="Po_psi" format="appended" offset="{offPo}"/>""");
            sb.AppendLine($"""        <DataArray type="Float64" Name="Pw_psi" format="appended" offset="{offPw}"/>""");
            sb.AppendLine($"""        <DataArray type="Float64" Name="Sw" format="appended" offset="{offSw}"/>""");
            sb.AppendLine($"""        <DataArray type="Float64" Name="So" format="appended" offset="{offSo}"/>""");

            if (_isBlackOil)
            {
                long offPg = offSo + 4 + scalarBytes;
                long offPb = offPg + 4 + scalarBytes;
                long offSg = offPb + 4 + scalarBytes;

                sb.AppendLine($"""        <DataArray type="Float64" Name="Pg_psi" format="appended" offset="{offPg}"/>""");
                sb.AppendLine($"""        <DataArray type="Float64" Name="Pb_psi" format="appended" offset="{offPb}"/>""");
                sb.AppendLine($"""        <DataArray type="Float64" Name="Sg" format="appended" offset="{offSg}"/>""");
            }

            sb.AppendLine("""      </CellData>""");
            sb.AppendLine("""    </Piece>""");
            sb.AppendLine("""  </UnstructuredGrid>""");
            sb.AppendLine("""  <AppendedData encoding="raw">""");
            sb.Append('_');

            return sb.ToString();
        }

        private void WritePvdCollection()
        {
            string pvdPath = Path.Combine(_outputDirectory, "SIMULATION.pvd");
            using var sw = new StreamWriter(pvdPath, append: false, new UTF8Encoding(false));
            sw.WriteLine("""<?xml version="1.0"?>""");
            sw.WriteLine("""<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""");
            sw.WriteLine("""  <Collection>""");
            foreach (var (time, file) in _snapshots)
            {
                sw.WriteLine(string.Create(CultureInfo.InvariantCulture,
                    $"""    <DataSet timestep="{time:G9}" group="" part="0" file="{file}"/>"""));
            }
            sw.WriteLine("""  </Collection>""");
            sw.WriteLine("""</VTKFile>""");
        }

        private static void WriteRawBlock<TData>(FileStream fs, TData[] data) where TData : unmanaged
        {
            ReadOnlySpan<byte> bytes = MemoryMarshal.AsBytes(data.AsSpan());
            Span<byte> lenHeader = stackalloc byte[4];
            BitConverter.TryWriteBytes(lenHeader, (uint)bytes.Length);
            fs.Write(lenHeader);
            fs.Write(bytes);
        }

        public void Dispose() => WritePvdCollection();
    }
}