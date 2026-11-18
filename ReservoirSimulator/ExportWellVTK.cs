using EclipseDeckReader;
using SepalSolver;
using System.Text;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    internal class ExportWellVTK
    {
        public static void ExportWells<T>(string outputDirectory, Grid<T>[] Grids, Well<T>[] Wells)
            where T : struct, IVectorND<T>
        {
            Directory.CreateDirectory(outputDirectory);
            foreach (var well in Wells)
            {
                double[][] Trajectory = well.GetTrajectoryCoordinates(Grids);
                ExportWellAsVtk(Path.Combine(outputDirectory, $"{well.Name}.vtk"), well.Name, Trajectory);
            }
        }

        public static void ExportWellAsVtk(string outputPath, string wellName, double[][] trajectories)
        {
            // trajectories array expected shape: [numPoints][3] -> X, Y, Z coordinates
            int numPoints = trajectories.Length;

            using StreamWriter writer = new(outputPath, false, Encoding.ASCII);

            // 1. Write standard VTK Polydata Header
            writer.WriteLine("# vtk DataFile Version 3.0");
            writer.WriteLine($"Well Path: {wellName}");
            writer.WriteLine("ASCII");
            writer.WriteLine("DATASET POLYDATA");

            // 2. Write Geometry Nodes (Points)
            writer.WriteLine($"POINTS {numPoints} double");
            for (int i = 0; i < numPoints; i++)
            {
                writer.WriteLine($"{trajectories[i][0]:F2} {trajectories[i][1]:F2} {trajectories[i][2]:F2}");
            }
            writer.WriteLine();

            // 3. Topology: Link the points into a single continuous wireframe string
            // Format: LINES [NumberOfLines] [TotalIntegerCountInBlock]
            // Total count = NumberOfLines + total number of indices mapped
            writer.WriteLine($"LINES 1 {numPoints + 1}");
            writer.Write($"{numPoints} ");
            for (int i = 0; i < numPoints; i++)
            {
                writer.Write($"{i} ");
            }
            writer.WriteLine();
            writer.WriteLine();

            // 4. Optional Metadata Attributions (Coloring / Labeling in ParaView)
            writer.WriteLine($"CELL_DATA 1");
            writer.WriteLine("FIELD FieldData 1");
            writer.WriteLine($"WellID 1 {wellName.Length} string");
            writer.WriteLine(wellName);
        }
    }
}
