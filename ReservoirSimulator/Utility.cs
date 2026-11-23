using SepalSolver;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;
using UnitSystem.Parsing;
using UnitSystem.Registry;
using UnitSystem.Services;

namespace ReservoirSimulator
{
    internal class Utility
    {

        const double alpha = 1.12712e-3, beta = 5.615;

        // Key: Quantity Name
        // Value: string[4] -> [0] = FIELD, [1] = METRIC, [2] = PVT-M, [3] = LAB
        public static readonly Dictionary<string, string[]> UnitsByQuantity = new(StringComparer.OrdinalIgnoreCase)
        {
            ["Length"] = ["ft", "m", "m", "cm"],
            ["Depth"] = ["ft", "m", "m", "cm"],
            ["Wellbore Diameter"] = ["ft", "m", "m", "cm"],
            ["Pressure"] = ["psi", "bar", "atm", "atm"],
            ["Capillary Pressure (Pc)"] = ["psi", "bar", "atm", "atm"],
            ["Compressibility"] = ["1/psi", "1/bar", "1/atm", "1/atm"],
            ["Permeability"] = ["mD", "mD", "mD", "mD"],
            ["Viscosity"] = ["cP", "cP", "cP", "cP"],
            ["Porosity"] = ["fraction", "fraction", "fraction", "fraction"],
            ["Transmissibility"] = ["cP*rb/(day*psi)", "cP*rm^3/(day*bar)", "cP*rm^3/(day*atm)", "cP*rcc/(hr*atm)"],
            ["Well Productivity Index (PI)"] = ["STB/(day*psi)", "sm^3/(day*bar)", "sm^3/(day*atm)", "scc/(hr*atm)"],
            ["Oil Surface Rate"] = ["STB/day", "sm^3/day", "sm^3/day", "scc/hr"],
            ["Water Surface Rate"] = ["STB/day", "sm^3/day", "sm^3/day", "scc/hr"],
            ["Liquid Surface Rate"] = ["STB/day", "sm^3/day", "sm^3/day", "scc/hr"],
            ["Gas Surface Rate"] = ["Mscf/day", "sm^3/day", "sm^3/day", "scc/hr"],
            ["Reservoir Volume Rate"] = ["rb/day", "rm^3/day", "rm^3/day", "rcc/hr"],
            ["Reservoir Volume"] = ["rb", "rm^3", "rm^3", "rcc"],
            ["Cumulative Liquid Volume"] = ["STB", "sm^3", "sm^3", "scc"],
            ["Cumulative Gas Volume"] = ["Mscf", "sm^3", "sm^3", "scc"],
            ["Oil Formation Volume Factor (Bo)"] = ["rb/STB", "rm^3/sm^3", "rm^3/sm^3", "rcc/scc"],
            ["Water Formation Volume Factor (Bw)"] = ["rb/STB", "rm^3/sm^3", "rm^3/sm^3", "rcc/scc"],
            ["Gas Formation Volume Factor (Bg)"] = ["rb/Mscf", "rm^3/sm^3", "rm^3/sm^3", "rcc/scc"],
            ["Solution Gas-Oil Ratio (Rs)"] = ["Mscf/STB", "sm^3/sm^3", "sm^3/sm^3", "scc/scc"],
            ["Vaporized Oil-Gas Ratio (Rv)"] = ["STB/Mscf", "sm^3/sm^3", "sm^3/sm^3", "scc/scc"],
            ["Density"] = ["lb/ft^3", "kg/m^3", "kg/m^3", "g/cm^3"],
            ["Temperature"] = ["°F", "°C", "°C", "°C"],
            ["Time"] = ["days", "days", "days", "hours"],
            ["Mass"] = ["lbm", "kg", "kg", "g"]
        };

        internal static double Convert(double value, string quantityName, Units fromUnits, Units toUnits)
        {
            if (!UnitsByQuantity.TryGetValue(quantityName, out var unitArray))
                throw new ArgumentException($"Unknown quantity name: {quantityName}");
            int fromIndex = (int)fromUnits;
            int toIndex = (int)toUnits;
            string fromUnit = unitArray[fromIndex];
            string toUnit = unitArray[toIndex];
            return UnitConverter.Convert(value, fromUnit, toUnit);
        }

        internal static IEnumerable<double> Convert(IEnumerable<double> values, string quantityName, Units fromUnits, Units toUnits)
        {
            if (!UnitsByQuantity.TryGetValue(quantityName, out var unitArray))
                throw new ArgumentException($"Unknown quantity name: {quantityName}");
            int fromIndex = (int)fromUnits;
            int toIndex = (int)toUnits;
            string fromUnit = unitArray[fromIndex];
            string toUnit = unitArray[toIndex];
            return values.Select(value => UnitConverter.Convert(value, fromUnit, toUnit));
        }

        internal static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            if (n < 2) return (1.0, 0.0, 0.0);

            double sumX = x.Sum(); double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Sum(xi => Sqr(xi));

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (Abs(denominator) < 1e-11) return (1.0, 0.0, 0.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            double yMean = sumY / n;
            double ssTot = y.Sum(yi => Sqr(yi - yMean));
            double ssRes = y.Zip(x, (yi, xi) => Sqr(yi - (m * xi + b))).Sum();
            double rSq = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return (m, b, rSq);
        }

        internal static void AddIfPresent(List<double> list, double? value)
        {
            if (value.HasValue)
                list.Add(value.GetValueOrDefault());
        }

        static Vector2d[] ProjectTo2D(Vector3d[] vertices, int dropAxis)
        {
            Vector2d[] result = new Vector2d[4];
            for (int i = 0; i < 4; i++)
            {
                if (dropAxis == 0) result[i] = new Vector2d(vertices[i].Y, vertices[i].Z);
                else if (dropAxis == 1) result[i] = new Vector2d(vertices[i].X, vertices[i].Z);
                else result[i] = new Vector2d(vertices[i].X, vertices[i].Y);
            }
            return result;
        }

        static List<Vector2d> ClipPolygons(Vector2d[] subjectPoly, Vector2d[] clipPoly)
        {
            List<Vector2d> outputList = new List<Vector2d>(subjectPoly);

            // Loop through each edge of the clipping polygon
            for (int i = 0; i < clipPoly.Length; i++)
            {
                Vector2d edgeStart = clipPoly[i];
                Vector2d edgeEnd = clipPoly[(i + 1) % clipPoly.Length];

                List<Vector2d> inputList = outputList;
                outputList = new List<Vector2d>();

                if (inputList.Count == 0) break;

                Vector2d s = inputList[inputList.Count - 1];

                foreach (Vector2d e in inputList)
                {
                    if (IsInsideEdge(e, edgeStart, edgeEnd))
                    {
                        if (!IsInsideEdge(s, edgeStart, edgeEnd))
                        {
                            outputList.Add(ComputeIntersection(s, e, edgeStart, edgeEnd));
                        }
                        outputList.Add(e);
                    }
                    else if (IsInsideEdge(s, edgeStart, edgeEnd))
                    {
                        outputList.Add(ComputeIntersection(s, e, edgeStart, edgeEnd));
                    }
                    s = e;
                }
            }
            return outputList;
        }

        static bool IsInsideEdge(Vector2d p, Vector2d edgeStart, Vector2d edgeEnd)
        {
            // Assumes clipPoly vertices are passed in anti-clockwise sequence
            return (edgeEnd.X - edgeStart.X) * (p.Y - edgeStart.Y) - (edgeEnd.Y - edgeStart.Y) * (p.X - edgeStart.X) >= -1e-10;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static double CalculateTetrahedronVolume(Vector3d v0, Vector3d v1, Vector3d v2, Vector3d v3)
        {
            // Vectorize edge profiles
            Vector3d edge1 = v1 - v0;
            Vector3d edge2 = v2 - v0;
            Vector3d edge3 = v3 - v0;

            // Triple scalar product: edge1 . (edge2 x edge3)
            Vector3d crossProduct = Vector3d.Cross(edge2, edge3);
            double tripleProduct = Vector3d.Dot(edge1, crossProduct);

            // Volume is exactly 1/6th of the parallelpiped volume
            return Abs(tripleProduct) / 6.0;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static double ComputeHexahedronVolume<T>(Grid<T> grid) where T : struct, IVectorND<T>
        {
            // Sanity Check: Ensure we have exactly 8 corners
            if (grid.Corners.Length != 8) throw new ArgumentException("Hexahedron must have exactly 8 corners.");

            // Sum the volumes of the 12 constituent tetrahedrons
            double totalVolume = 0;
            foreach(var face in grid.Faces)
            {
                totalVolume += CalculateTetrahedronVolume(grid.Center, face.Corners[0], face.Corners[1], face.Corners[2]);
                totalVolume += CalculateTetrahedronVolume(grid.Center, face.Corners[0], face.Corners[2], face.Corners[3]);
            }

            // Guard against degenerate/pinched-out dead cells
            if (double.IsNaN(totalVolume) || totalVolume < 1e-6)  return 0.0;

            return totalVolume;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]

        internal static double ComputePerpendicularDistance(Vector3d cellCenter, Vector3d p0, Vector3d p1, Vector3d p2, Vector3d p3)
        {
            // Step 1: Calculate two diagonal vectors from the face corners
            Vector3d d1 = p2 - p0,  d2 = p3 - p1;

            // Step 2 & 3: Compute the face normal and normalize it to a unit vector
            Vector3d normalLarge = Vector3d.Cross(d1, d2);
            double magnitude = normalLarge.Magnitude();

            // Guard against completely collapsed/zero-area faces
            if (magnitude < 1e-12) return 0.0;
            Vector3d unitNormal = (1.0 / magnitude) * normalLarge;

            // Step 4: Anchor at the face centroid rather than a single corner
            Vector3d faceCenter = 0.25 * (p0 + p1 + p2 + p3);
            Vector3d vCenter = cellCenter - faceCenter;

            // Step 5: Dot product projects the center vector onto the perpendicular normal track
            double distance = Vector3d.Dot(vCenter, unitNormal);

            // Return absolute distance (direction relative to the normal does not matter for transmissibility)
            return Abs(distance);
        }

        static Vector2d ComputeIntersection(Vector2d s, Vector2d e, Vector2d edgeStart, Vector2d edgeEnd)
        {
            double dcX = edgeStart.X - edgeEnd.X;
            double dcY = edgeStart.Y - edgeEnd.Y;
            double dpX = s.X - e.X;
            double dpY = s.Y - e.Y;

            double n1 = edgeStart.X * edgeEnd.Y - edgeStart.Y * edgeEnd.X;
            double n2 = s.X * e.Y - s.Y * e.X;
            double denominator = dcX * dpY - dcY * dpX;

            if (Abs(denominator) < 1e-12) return s; // Parallel line guard

            return new Vector2d(
                (n1 * dpX - dcX * n2) / denominator,
                (n1 * dpY - dcY * n2) / denominator
            );
        }

        static double Compute2DArea(List<Vector2d> polygon)
        {
            // Standard Shoelace Formula
            double area = 0.0;
            int j = polygon.Count - 1;
            for (int i = 0; i < polygon.Count; i++)
            {
                area += (polygon[j].X + polygon[i].X) * (polygon[j].Y - polygon[i].Y);
                j = i;
            }
            return Abs(area * 0.5);
        }

        internal static List<PressureDepthPair> IntegrateLegRK4(double zStart, double pStart, double zEnd, int N, Func<double, double> getGradient)
        {
            if (N <= 0 || zStart == zEnd)
                return [new PressureDepthPair(zStart, pStart), new PressureDepthPair(zEnd, pStart)];

            var table = new List<PressureDepthPair>();
            double currentZ = zStart;
            double currentP = pStart;
            double h = (zEnd - zStart) / N;
            table.Add(new(currentZ, currentP));
            for (int n = 0; n < N; n++)
            {
                double k1 = getGradient(currentP);
                double k2 = getGradient(currentP + 0.5 * h * k1);
                double k3 = getGradient(currentP + 0.5 * h * k2);
                double k4 = getGradient(currentP + h * k3);
                currentP += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                currentZ += h;
                table.Add(new(currentZ, currentP));
            }
            return table;
        }

        internal static double Interpolate<T>(List<T> rows, double targetValue, Func<T, double> searchSelector, Func<T, double> targetSelector)
        {
            if (rows == null)
                throw new ArgumentNullException(nameof(rows), "Table data cannot be null.");

            int count = rows.Count;
            if (count == 0)
                throw new ArgumentException("Table data cannot be empty.", nameof(rows));

            if (count == 1) return targetSelector(rows[0]);

            // CRITICAL SECURITY FIX: Order the dataset by the actual independent variable we are querying.
            // This guarantees the monotonicity required for an O(log N) binary search.
            var optimizedRows = rows.OrderBy(searchSelector).ToList();

            // Since we explicitly ordered by searchSelector, it is GUARANTEED to be ascending!
            double firstSearchVal = searchSelector(optimizedRows[0]);
            double lastSearchVal = searchSelector(optimizedRows[count - 1]);

            // 3. BOUNDARY PROTECTION: Clamp overshoots cleanly to the table edge values
            if (targetValue <= firstSearchVal) return targetSelector(optimizedRows[0]);
            if (targetValue >= lastSearchVal) return targetSelector(optimizedRows[count - 1]);

            // 4. BINARY SEARCH: Fast interval bracketing on a guaranteed ascending array
            int low = 0;
            int high = count - 1;
            int index = 0;

            while (low <= high)
            {
                int mid = low + (high - low) / 2;
                double midVal = searchSelector(optimizedRows[mid]);

                if (midVal == targetValue)
                {
                    return targetSelector(optimizedRows[mid]);
                }

                if (midVal < targetValue)
                {
                    index = mid; // Lower bound bracket candidate
                    low = mid + 1;
                }
                else
                {
                    high = mid - 1;
                }
            }

            // 5. BLEND THE BRACKETS
            var rowLeft = optimizedRows[index];
            var rowRight = optimizedRows[index + 1];

            double xLeft = searchSelector(rowLeft);
            double xRight = searchSelector(rowRight);

            double yLeft = targetSelector(rowLeft);
            double yRight = targetSelector(rowRight);

            double fraction = (targetValue - xLeft) / (xRight - xLeft);

            return yLeft + fraction * (yRight - yLeft);
        }

        public static int SetAcitveGrids<T>(Grid<T>[] Grids, List<int> Actnum, Func<ADiff<T>, ADiff<T>> Er) where T : struct, IVectorND<T>
        {
            int numActiveGrids = 0;
            for (int m = 0; m < Grids.Length; m++)
            {
                var grid = Grids[m];
                grid.GlobalIndex = m;
                if (Actnum[m] > 0)
                {
                    grid.Er = Er;
                    grid.ActiveIndex = numActiveGrids++;
                    foreach (var face in grid.Faces)
                        face.GlobalGridIndex = m;
                }
            }
            return numActiveGrids;
        }

        public static void SetPVDO<T>(Grid<T>[] Grids, List<int> Pvtnum, PVDOAnalysis[] pvdo, double γo0) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                int n = Pvtnum[grid.GlobalIndex];   
                grid.Bo = Po => pvdo[n].Result.EvaluateBo(Po);
                grid.μo = Po => pvdo[n].Result.EvaluateMuo(Po);
                grid.γo = Po => γo0 / pvdo[n].Result.EvaluateBo(Po);
            }
        }

        public static void SetPVDG<T>(Grid<T>[] Grids, List<int> Pvtnum, PVDGAnalysis[] pvdg, double γg0) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                int n = Pvtnum[grid.GlobalIndex];
                grid.Bg = Pg => pvdg[n].Result.EvaluateBg(Pg);
                grid.μg = Pg => pvdg[n].Result.EvaluateMug(Pg);
                grid.γg = Pg => γg0 / pvdg[n].Result.EvaluateBg(Pg);
            }
        }

        public static void SetPVTO<T>(Grid<T>[] Grids, List<int> Pvtnum, PVTOAnalysis[] pvto, double γo0, double γg0) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                //if (!grid.isActive) continue;
                //int n = Pvtnum[grid.GlobalIndex];
                //grid.Rs = Po => pvto[n].Result.EvaluateRs(Po);
                //grid.Bo = Rs => pvto[n].Result.EvaluateBo(Rs);
                //grid.μo = Rs => pvto[n].Result.EvaluateMuo(Rs);
                //grid.γo = Rs => (γo0 + 0.178 * Rs * γg0) / pvto[n].Result.EvaluateBo(Rs);
            }
        }

        public static void SetPVTG<T>(Grid<T>[] Grids, List<int> Pvtnum, PVTGAnalysis[] pvtg, double γo0, double γg0) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                int n = Pvtnum[grid.GlobalIndex];
                grid.Rv = Pg => pvtg[n].Result.EvaluateRv(Pg);
                grid.Bg = Rv => pvtg[n].Result.EvaluateBg(Rv);
                grid.μg = Rv => pvtg[n].Result.EvaluateMug(Rv);
                grid.γg = Rv => (γg0 + 5.615 * Rv * γo0) / pvtg[n].Result.EvaluateBg(Rv);
            }
        }

        public static void SetPVTW<T>(Grid<T>[] Grids, List<int> Pvtnum, PVTWAnalysis[] pvtw, double γw0) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                int n = Pvtnum[grid.GlobalIndex];
                grid.Bw = Pw => pvtw[n].Result.EvaluateBw(Pw);
                grid.μw = Pw => pvtw[n].Result.EvaluateMuw(Pw);
                grid.γw = Pw => γw0 / pvtw[n].Result.EvaluateBw(Pw);
            }
        }

        public static void SetSWOF<T>(Grid<T>[] Grids, List<int> Satnum, SWOFAnalysis[] swof) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                int n = Satnum[grid.GlobalIndex];   
                grid.Krw = Sw => swof[n].Result.EvaluateKrw(Sw);
                grid.Krow = Sw => swof[n].Result.EvaluateKrow(Sw);
                grid.Pcow = Sw => swof[n].Result.EvaluatePcow(Sw);
            }
        }

        public static void SetSGOF<T>(Grid<T>[] Grids, List<int> Satnum, SGOFAnalysis[] sgof) where T : struct, IVectorND<T>
        {
            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                int n = Satnum[grid.GlobalIndex];
                grid.Krg = Sg => sgof[n].Result.EvaluateKrg(Sg);
                grid.Krog = Sg => sgof[n].Result.EvaluateKrog(Sg);
                grid.Pcog = Sg => sgof[n].Result.EvaluatePcog(Sg);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static ADiff<T> SmoothUpstream<T>(ADiff<T> dPhi, ADiff<T> mob_i, ADiff<T> mob_j, double eps = 1e-2) where T : struct, IVectorND<T>
        {
            double v = dPhi.Value;
            if (v >= eps) return mob_i;
            if (v <= -eps) return mob_j;

            // C1-continuous cubic Hermite transition on [-eps, +eps]
            var t = dPhi / eps;
            var w = 0.5 + 0.75 * t - 0.25 * t * t * t;
            return w * mob_i + (1.0 - w) * mob_j;
        }

        internal static Grid<T>[] ProcessGrid<T>(int nx, int ny, int nz, double[] coord, double[] zcorn,
            double[] Φ, double[] permx, double[] permy, double[] permz, double[]? ntg = null, double[]? multpv = null) where T : struct, IVectorND<T>
        {
            int nxny = nx * ny, ngrids = nxny * nz, nx2 = nx * 2, ny2 = ny * 2;
            Grid<T>[] Grids = new Grid<T>[ngrids];
            for (int m = 0; m < ngrids; m++)
            {
                int k = m / nxny, rem = m % nxny, j = rem / nx, i = rem % nx;
                Grids[m] = new Grid<T>
                {
                    IJKIndex = (i, j, k),
                    NTG = ntg?[m] ?? 1.0,
                    PoreVolumeMultiplier = multpv?[m] ?? 1.0,
                    Porosity = Φ[m],
                    Corners = new Vector3d[8]
                };
                Grids[m].Faces[0].Permeability = permx[m];
                Grids[m].Faces[1].Permeability = permx[m];
                Grids[m].Faces[2].Permeability = permy[m];
                Grids[m].Faces[3].Permeability = permy[m];
                Grids[m].Faces[4].Permeability = permz[m];
                Grids[m].Faces[5].Permeability = permz[m];

                // Total number of coordinate lines in X and Y directions
                int numLinesX = nx + 1;

                // Loop through top (0-3) and bottom (4-7) layers
                for (int c = 0; c < 8; c++)
                {
                    // Determine local binary topological offsets
                    int ix = (c == 1 || c == 2 || c == 5 || c == 6) ? 1 : 0;
                    int iy = (c == 2 || c == 3 || c == 6 || c == 7) ? 1 : 0;
                    int iz = (c >= 4) ? 1 : 0;

                    // 1. Get the correct Coordinate Line index for this corner
                    int lineI = i + ix, lineJ = j + iy;
                    int lineStride = lineI + lineJ * numLinesX;

                    // Extract Top and Bottom anchor points from the raw COORD array
                    // Layout of COORD: [xtop, ytop, ztop, xbot, ybot, zbot] for each line
                    int coordIdx = lineStride * 6;
                    Vector3d pTop = new(coord[coordIdx + 0], coord[coordIdx + 1], coord[coordIdx + 2]);
                    Vector3d pBot = new(coord[coordIdx + 3], coord[coordIdx + 4], coord[coordIdx + 5]);

                    // 2. Fetch the explicit target Z-coordinate from ZCORN array
                    // Layout of ZCORN uses independent 2nx * 2ny * 2nz index structure
                    int zIdx_X = i * 2 + ix, zIdx_Y = j * 2 + iy, zIdx_Z = k * 2 + iz;

                    int zcornIdx = zIdx_X + (zIdx_Y * nx2) + (zIdx_Z * nx2 * ny2);
                    double targetZ = zcorn[zcornIdx];

                    // 3. Linear interpolation along the line vector using register math
                    double fraction = (targetZ - pTop.Z) / (pBot.Z - pTop.Z);

                    // Check for dead/pinch-out cells where top and bottom meet
                    if (double.IsNaN(fraction) || double.IsInfinity(fraction)) fraction = 0.0;

                    // X and Y are calculated seamlessly via your custom Vector3d scalar multiplication operator
                    Vector3d xyPlaneInterpolated = pTop + fraction * (pBot - pTop);

                    // Save the exact 3D location (Overriding the Z lane with the explicit target depth)
                    Grids[m].Corners[c] = new Vector3d(xyPlaneInterpolated.X, xyPlaneInterpolated.Y, targetZ);
                }
                Grids[m].MakeFaces();

            }
            return Grids;
        }
        public static Grid<T>[] ProcessGrid<T>(int nx, int ny, int nz, double[] dx, double[] dy, double[] dz,
            double[] zTop, double[] Φ, double[] permx, double[] permy, double[] permz, double[]? ntg = null, double[]? multpv = null) where T : struct, IVectorND<T>
        {
            int nxny = nx * ny, ngrids = nxny * nz;
            Grid<T>[] Grids = new Grid<T>[ngrids];
            for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        int m = i + j * nx + k * nxny;
                        var I = Enumerable.Range(0, i + 1);
                        var J = Enumerable.Range(0, j + 1);
                        var K = Enumerable.Range(0, k + 1);
                        Grids[m] = new Grid<T>
                        {
                            IJKIndex = (i, j, k),
                            NTG = ntg?[m] ?? 1.0,
                            PoreVolumeMultiplier = multpv?[m] ?? 1.0,
                            Porosity = Φ[m],
                            Corners = new Vector3d[8]
                        };
                        Grids[m].Faces[0].Permeability = permz[m];
                        Grids[m].Faces[1].Permeability = permy[m];
                        Grids[m].Faces[2].Permeability = permx[m];
                        Grids[m].Faces[3].Permeability = permx[m];
                        Grids[m].Faces[4].Permeability = permy[m];
                        Grids[m].Faces[5].Permeability = permz[m];
                        double zlo = zTop[j * nx + i] + K.Select(k => dz[i + j * nx + k * nxny]).Sum();
                        double zhi = zlo - dz[m];
                        double xhi = I.Select(i => dx[i + j * nx + k * nxny]).Sum();
                        double xlo = xhi - dx[m];
                        double yhi = J.Select(j => dy[i + j * nx + k * nxny]).Sum();
                        double ylo = yhi - dy[m];

                        // Extract the 8 structural corner points for a Hexahedron matching VTK specifications
                        // Bottom plane points (Z-Upper)
                        Grids[m].Corners[0] = new Vector3d(xlo, ylo, zhi);     // Point 0
                        Grids[m].Corners[1] = new Vector3d(xhi, ylo, zhi);     // Point 1
                        Grids[m].Corners[2] = new Vector3d(xhi, yhi, zhi);     // Point 2
                        Grids[m].Corners[3] = new Vector3d(xlo, yhi, zhi);     // Point 3
                                                                               // Top plane points (Z-Lower)
                        Grids[m].Corners[4] = new Vector3d(xlo, ylo, zlo);     // Point 4
                        Grids[m].Corners[5] = new Vector3d(xhi, ylo, zlo);     // Point 5
                        Grids[m].Corners[6] = new Vector3d(xhi, yhi, zlo);     // Point 6
                        Grids[m].Corners[7] = new Vector3d(xlo, yhi, zlo);     // Point 7
                        Grids[m].MakeFaces();
                    }
                }
            }
            
            return Grids;
        }

        static double CalculateFaultIntersectionArea(Face f1, Face f2)
        {
            Vector3d diagA = f1.Corners[2] - f1.Corners[0];
            Vector3d diagB = f1.Corners[3] - f1.Corners[1];
            Vector3d normal = Vector3d.Cross(diagA, diagB);

            double nx = Abs(normal.X);
            double ny = Abs(normal.Y);
            double nz = Abs(normal.Z);
            double len = Sqrt(nx * nx + ny * ny + nz * nz);
            if (len < 1e-12) return 0.0;

            // Drop the dominant normal axis to maximize 2D projection precision
            int dropAxis = 2;
            if (nx > ny && nx > nz) dropAxis = 0;
            else if (ny > nx && ny > nz) dropAxis = 1;

            Vector2d[] poly1 = ProjectTo2D(f1.Corners, dropAxis);
            Vector2d[] poly2 = ProjectTo2D(f2.Corners, dropAxis);

            // Essential: f1 and f2 have opposite 3D outward normals, so one projects CW and one CCW.
            // Normalizing both to CCW prevents Sutherland-Hodgman from inverting the clip half-planes.
            EnsureCounterClockwise(poly1);
            EnsureCounterClockwise(poly2);

            List<Vector2d> clippedPolygon = ClipPolygons(poly1, poly2);
            if (clippedPolygon.Count < 3) return 0.0;

            double area2D = Compute2DArea(clippedPolygon);
            double cosTheta = (dropAxis == 0 ? nx : (dropAxis == 1 ? ny : nz)) / len;

            return area2D / cosTheta;
        }

        static void EnsureCounterClockwise(Vector2d[] poly)
        {
            double signedArea2 = 0.0;
            for (int i = 0; i < poly.Length; i++)
            {
                Vector2d a = poly[i];
                Vector2d b = poly[(i + 1) % poly.Length];
                signedArea2 += (a.X * b.Y - b.X * a.Y);
            }
            if (signedArea2 < 0.0)
                Array.Reverse(poly);
        }

        public static void BuildStructuredConnections<T>(
            Grid<T>[] Grids, int nx, int ny, int nz)
            where T : struct, IVectorND<T>
        {
            int Idx(int i, int j, int k) => i + nx * (j + ny * k);

            for (int m = 0; m < Grids.Length; m++)
            {
                var g = Grids[m];
                if (g?.Faces == null) continue;
                for (int f = 0; f < 6; f++)
                {
                    g.Faces[f].GridIndices.Clear();
                    g.Faces[f].Transmissibility.Clear();
                }
            }

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    // -----------------------------------------------------------------
                    // 1. +I Interfaces: Column (i, j) <-> Column (i + 1, j)
                    //    Faces[1] (+I) <-> Faces[0] (-I)
                    // -----------------------------------------------------------------
                    if (i + 1 < nx)
                    {
                        if (IsPillarFaulted(Grids, nx, ny, nz, i, j, 1, i + 1, j, 0))
                        {
                            ConnectFaultedColumnPair(Grids, nx, ny, nz, i, j, 1, i + 1, j, 0, alpha, applyNTG: true);
                        }
                        else
                        {
                            for (int k = 0; k < nz; k++)
                            {
                                var c1 = Grids[Idx(i, j, k)];
                                var c2 = Grids[Idx(i + 1, j, k)];
                                if (!c1.isActive || !c2.isActive) continue;

                                double area = 0.5 * (c1.Faces[1].Area + c2.Faces[0].Area);
                                RegisterGridConnection(c1.Faces[1], c2.Faces[0], area, c1.NTG, c2.NTG);
                            }
                        }
                    }

                    // -----------------------------------------------------------------
                    // 2. +J Interfaces: Column (i, j) <-> Column (i, j + 1)
                    //    Faces[3] (+J) <-> Faces[2] (-J)
                    // -----------------------------------------------------------------
                    if (j + 1 < ny)
                    {
                        if (IsPillarFaulted(Grids, nx, ny, nz, i, j, 3, i, j + 1, 2))
                        {
                            ConnectFaultedColumnPair(Grids, nx, ny, nz, i, j, 3, i, j + 1, 2, alpha, applyNTG: true);
                        }
                        else
                        {
                            for (int k = 0; k < nz; k++)
                            {
                                var c1 = Grids[Idx(i, j, k)];
                                var c2 = Grids[Idx(i, j + 1, k)];
                                if (!c1.isActive || !c2.isActive) continue;

                                double area = 0.5 * (c1.Faces[3].Area + c2.Faces[2].Area);
                                RegisterGridConnection(c1.Faces[3], c2.Faces[2], area, c1.NTG, c2.NTG);
                            }
                        }
                    }

                    // -----------------------------------------------------------------
                    // 3. +K Interfaces: Down Column (i, j) from k -> k + 1
                    //    Faces[5] (+K) <-> Faces[4] (-K), skipping zero-thickness pinch-outs
                    // -----------------------------------------------------------------
                    for (int k = 0; k + 1 < nz; k++)
                    {
                        var c1 = Grids[Idx(i, j, k)];
                        if (!c1.isActive) continue;

                        int kNext = k + 1;
                        while (kNext < nz && !Grids[Idx(i, j, kNext)].isActive &&
                               Grids[Idx(i, j, kNext)].Faces[1].Area <= 1e-9)
                        {
                            kNext++;
                        }

                        if (kNext < nz)
                        {
                            var c2 = Grids[Idx(i, j, kNext)];
                            if (!c2.isActive) continue;

                            double area = 0.5 * (c1.Faces[5].Area + c2.Faces[4].Area);
                            RegisterGridConnection(c1.Faces[5], c2.Faces[4], area);
                        }
                    }
                }
            }
        }


        static bool IsPillarFaulted<T>(
            Grid<T>[] Grids, int nx, int ny, int nz,
            int iA, int jA, int faceIdxA,
            int iB, int jB, int faceIdxB,
            double tol = 1e-4)
            where T : struct, IVectorND<T>
        {
            int Idx(int ii, int jj, int kk) => ii + nx * (jj + ny * kk);

            for (int k = 0; k < nz; k++)
            {
                var cA = Grids[Idx(iA, jA, k)].Corners;
                var cB = Grids[Idx(iB, jB, k)].Corners;

                if (faceIdxA == 1) // +I (cA: 1,2,6,5) vs -I (cB: 0,3,7,4)
                {
                    if (Abs(cA[1].Z - cB[0].Z) > tol ||
                        Abs(cA[2].Z - cB[3].Z) > tol ||
                        Abs(cA[6].Z - cB[7].Z) > tol ||
                        Abs(cA[5].Z - cB[4].Z) > tol)
                        return true;
                }
                else // +J (cA: 3,2,6,7) vs -J (cB: 0,1,5,4)
                {
                    if (Abs(cA[3].Z - cB[0].Z) > tol ||
                        Abs(cA[2].Z - cB[1].Z) > tol ||
                        Abs(cA[6].Z - cB[5].Z) > tol ||
                        Abs(cA[7].Z - cB[4].Z) > tol)
                        return true;
                }
            }
            return false;
        }

        static void ConnectFaultedColumnPair<T>(
            Grid<T>[] Grids, int nx, int ny, int nz,
            int iA, int jA, int faceIdxA,
            int iB, int jB, int faceIdxB,
            double alpha, bool applyNTG)
            where T : struct, IVectorND<T>
        {
            int Idx(int ii, int jj, int kk) => ii + nx * (jj + ny * kk);

            for (int kA = 0; kA < nz; kA++)
            {
                var cA = Grids[Idx(iA, jA, kA)];
                if (!cA.isActive) continue;

                var fA = cA.Faces[faceIdxA];
                if (fA.Area <= 1e-9) continue;

                GetFaceZBounds(fA, out double minZA, out double maxZA);

                for (int kB = 0; kB < nz; kB++)
                {
                    var cB = Grids[Idx(iB, jB, kB)];
                    if (!cB.isActive) continue;

                    var fB = cB.Faces[faceIdxB];
                    if (fB.Area <= 1e-9) continue;

                    GetFaceZBounds(fB, out double minZB, out double maxZB);

                    if (maxZA <= minZB + 1e-6 || maxZB <= minZA + 1e-6)
                        continue;

                    double sharedArea = CalculateFaultIntersectionArea(fA, fB);
                    if (sharedArea > 1e-6)
                    {
                        RegisterGridConnection(fA, fB, sharedArea, applyNTG ? cA.NTG : 1.0, applyNTG ? cB.NTG : 1.0);
                    }
                }
            }
        }

        // NTG scales each cell's permeability before the harmonic average, as in ECLIPSE (X and Y connections only;
        // vertical connections pass 1). Averaging the two NTGs instead overstates flow when they differ.
        static void RegisterGridConnection(Face f1, Face f2, double effectiveArea, double ntg1 = 1.0, double ntg2 = 1.0)
        {
            if (effectiveArea <= 0.0 || f1.PerpendicularDistance2Center <= 0.0 || f2.PerpendicularDistance2Center <= 0.0)
                return;

            double x1 = f1.Permeability * ntg1 / f1.PerpendicularDistance2Center;
            double x2 = f2.Permeability * ntg2 / f2.PerpendicularDistance2Center;

            double denom = x1 + x2;
            double T = (denom > 0.0) ? (alpha * effectiveArea * (x1 * x2) / denom) : 0.0;

            if (T <= 0.0) return;

            AddConnection(f1, f2.GlobalGridIndex, T);
            AddConnection(f2, f1.GlobalGridIndex, T);
        }

        static void AddConnection(Face f, int neighborGlobalIdx, double trans)
        {
            int existing = f.GridIndices.IndexOf(neighborGlobalIdx);
            if (existing >= 0)
            {
                f.Transmissibility[existing] += trans;
            }
            else
            {
                f.GridIndices.Add(neighborGlobalIdx);
                f.Transmissibility.Add(trans);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void GetFaceZBounds(Face f, out double minZ, out double maxZ)
        {
            double z0 = f.Corners[0].Z, z1 = f.Corners[1].Z, z2 = f.Corners[2].Z, z3 = f.Corners[3].Z;
            minZ = Min(Min(z0, z1), Min(z2, z3));
            maxZ = Max(Max(z0, z1), Max(z2, z3));
        }

        internal static void checkconnectivity<T>(Grid<T>[] Grids, int JacLength) where T : struct, IVectorND<T>
        {
            {
                int activeCells = 0;
                long totalUniqueNeighborPairs = 0;
                int[] neighborHistogram = new int[10]; // counts cells with 0, 1, 2, ..., 6+ neighbors

                foreach (var grid in Grids)
                {
                    if (!grid.isActive) continue;
                    activeCells++;

                    var uniqueNeighbors = new HashSet<int>();
                    foreach (var face in grid.Faces)
                    {
                        for (int c = 0; c < face.GridIndices.Count; c++)
                        {
                            int nIdx = face.GridIndices[c];
                            if (Grids[nIdx].isActive && face.Transmissibility[c] > 0.0)
                                uniqueNeighbors.Add(Grids[nIdx].ActiveIndex);
                        }
                    }

                    int count = uniqueNeighbors.Count;
                    totalUniqueNeighborPairs += count;
                    neighborHistogram[Min(count, 9)]++;
                }

                long expectedJacNnz = 4L * activeCells + 4L * totalUniqueNeighborPairs;
                Console.WriteLine($"Active cells: {activeCells} | Total neighbor links: {totalUniqueNeighborPairs}");
                Console.WriteLine($"Expected JacValues.Length: {expectedJacNnz} | Actual: {JacLength}");
                Console.WriteLine("Neighbor count distribution (0..6+): " +
                    string.Join(", ", Enumerable.Range(0, 10).Select(k => $"{k}:{neighborHistogram[k]}")));
            }
        }

    }
}
