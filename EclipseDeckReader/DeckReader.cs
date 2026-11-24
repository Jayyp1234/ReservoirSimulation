using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace EclipseDeckReader
{
    public sealed class DeckReader
    {
        // Current BOX (1-based, inclusive) set by BOX and cleared by ENDBOX; null means the whole grid.
        private (int I1, int I2, int J1, int J2, int K1, int K2)? box;

        // Section being read (RUNSPEC, GRID, ..., SCHEDULE)
        private string section = "RUNSPEC";

        // Skipped keywords already warned about, so a keyword repeated through the schedule is reported once
        private readonly HashSet<string> warnedKeywords = [];

        // Folder of the main deck: relative INCLUDE paths are looked up here first (ECLIPSE/OPM), then next to the
        // file that contains the INCLUDE
        private string mainDeckDirectory = "";

        public EclipseDataDeck LoadDeckFile(string filepath)
        {
            // Reading state starts afresh for every deck, in case this reader is used again
            box = null; section = "RUNSPEC"; warnedKeywords.Clear();
            mainDeckDirectory = Path.GetDirectoryName(Path.GetFullPath(filepath)) ?? "";

            var deck = new EclipseDataDeck();
            if (!File.Exists(filepath))
                throw new FileNotFoundException("Master deck file not found.", filepath);

            Console.WriteLine($"\n[Pipeline] MASTER DECK OPENED: {Path.GetFileName(filepath)}");
            string rawText = File.ReadAllText(filepath);
            List<Token> tokens = EclipseTokenizer.Tokenize(rawText.AsMemory());
            ParseTokenStream(tokens, deck, filepath);
            return deck;
        }

        private void ParseSubFile(string subFilePath, EclipseDataDeck deck)
        {
            Console.WriteLine($"   └─ [Pipeline] Entering Include File: {Path.GetFileName(subFilePath)}");
            string rawText = File.ReadAllText(subFilePath);
            List<Token> tokens = EclipseTokenizer.Tokenize(rawText.AsMemory());
            ParseTokenStream(tokens, deck, subFilePath);
            Console.WriteLine($"   └─ [Pipeline] Exiting Include File: {Path.GetFileName(subFilePath)}");
        }

        private void ParseTokenStream(List<Token> tokens, EclipseDataDeck deck, string currentFilePath)
        {
            int i = 0;
            while (i < tokens.Count)
            {
                if (tokens[i].Type == TokenType.Keyword)
                {
                    string keyword = tokens[i].Content.ToString().ToUpperInvariant();
                    int keywordIndex = i;
                    i++;

                    ClearLeadingLayoutSpacers(tokens, ref i);

                    switch (keyword)
                    {
                        case "DIMENS":
                        case "SPECGRID":
                            var dim = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            if (dim.Count >= 3)
                            {
                                deck.NX = dim[0]; deck.NY = dim[1]; deck.NZ = dim[2];
                                InitializeArrays(deck);
                            }
                            break;

                        case "OIL":
                        case "WATER":
                        case "GAS":
                        case "DISGAS":
                        case "VAPOIL":
                            if (!deck.PhaseFlags.Contains(keyword))
                                deck.PhaseFlags.Add(keyword);
                            break;

                        case "FIELD":
                        case "METRIC":
                        case "LAB":
                        case "PVT-M":
                            deck.UnitSystem = keyword;
                            Console.WriteLine($"      [Runspec] Unit System configured: {deck.UnitSystem}");
                            break;

                        // Section headers: remembered so the skip warning below can ignore SUMMARY (output only)
                        case "RUNSPEC":
                        case "GRID":
                        case "EDIT":
                        case "PROPS":
                        case "REGIONS":
                        case "SOLUTION":
                        case "SUMMARY":
                        case "SCHEDULE":
                            section = keyword;
                            break;

                        // The title is free text on the next line, with no '/'
                        case "TITLE":
                            while (i < tokens.Count && tokens[i].Type != TokenType.NewLine) i++;
                            break;

                        // Standalone control/summary flags with no payload and no '/'
                        case "ECHO":
                        case "NOECHO":
                        case "FOPR":
                        case "FGOR":
                        case "FWPR":
                        case "FWIR":
                        case "FGIR":
                        case "FPR":
                            break;

                        case "END":
                            Console.WriteLine("      [Pipeline] END keyword reached. Halting deck parse.");
                            return;

                        case "DX":
                            deck.DxArray = StoreGridArray(deck.DxArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current DX array status: {deck.DxArray.Count} entries.");
                            break;

                        case "DY":
                            deck.DyArray = StoreGridArray(deck.DyArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current DY array status: {deck.DyArray.Count} entries.");
                            break;

                        case "DZ":
                            deck.DzArray = StoreGridArray(deck.DzArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current DZ array status: {deck.DzArray.Count} entries.");
                            break;

                        case "TOPS":
                            deck.ZtopArray = StoreGridArray(deck.ZtopArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current TOPS (Ztop) array status: {deck.ZtopArray.Count} entries.");
                            break;

                        case "PORO":
                            deck.PorosityArray = StoreGridArray(deck.PorosityArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current PORO array status: {deck.PorosityArray.Count} entries.");
                            break;

                        case "PERMX":
                            deck.PermXArray = StoreGridArray(deck.PermXArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current PERMX array status: {deck.PermXArray.Count} entries.");
                            break;

                        case "PERMY":
                            deck.PermYArray = StoreGridArray(deck.PermYArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current PERMY array status: {deck.PermYArray.Count} entries.");
                            break;

                        case "PERMZ":
                            deck.PermZArray = StoreGridArray(deck.PermZArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 0.0);
                            Console.WriteLine($"      [Data] Current PERMZ array status: {deck.PermZArray.Count} entries.");
                            break;

                        case "NTG":
                            deck.NtgArray = StoreGridArray(deck.NtgArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 1.0);
                            Console.WriteLine($"      [Data] Current NTG array status: {deck.NtgArray.Count} entries.");
                            break;

                        case "MULTX":
                            deck.MultXArray = StoreGridArray(deck.MultXArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 1.0);
                            Console.WriteLine($"      [Data] Current MULTX array status: {deck.MultXArray.Count} entries.");
                            break;

                        case "MULTY":
                            deck.MultYArray = StoreGridArray(deck.MultYArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 1.0);
                            Console.WriteLine($"      [Data] Current MULTY array status: {deck.MultYArray.Count} entries.");
                            break;

                        case "MULTZ":
                            deck.MultZArray = StoreGridArray(deck.MultZArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 1.0);
                            Console.WriteLine($"      [Data] Current MULTZ array status: {deck.MultZArray.Count} entries.");
                            break;

                        case "MULTPV":
                            deck.MultpvArray = StoreGridArray(deck.MultpvArray, KeywordExtractor.ExtractDoubleArray(tokens, ref i), deck, 1.0);
                            Console.WriteLine($"      [Data] Current MULTPV array status: {deck.MultpvArray.Count} entries.");
                            break;

                        case "MINPV":
                            var minpv = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (minpv.Count > 0) deck.MinPoreVolume = minpv[0];
                            Console.WriteLine($"      [Data] MINPV: cells with pore volume below {deck.MinPoreVolume} will be made inactive.");
                            break;

                        case "BOX":
                            var b = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            if (b.Count >= 6) box = (b[0], b[1], b[2], b[3], b[4], b[5]);
                            Console.WriteLine($"      [Data] BOX set to I:{b.ElementAtOrDefault(0)}-{b.ElementAtOrDefault(1)}, J:{b.ElementAtOrDefault(2)}-{b.ElementAtOrDefault(3)}, K:{b.ElementAtOrDefault(4)}-{b.ElementAtOrDefault(5)}");
                            break;

                        case "ENDBOX":
                            box = null;
                            break;

                        case "ACTNUM":
                            var actData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.ActnumArray = StoreGridArray(deck.ActnumArray, actData, deck, 1);
                            Console.WriteLine($"      [Data] Current ACTNUM array status: {deck.ActnumArray.Count} entries.");
                            break;

                        case "COORD":
                            var coordData = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (deck.CoordArray.Count < coordData.Count) { deck.CoordArray = coordData; }
                            Console.WriteLine($"      [Data] Current COORD array status: {deck.CoordArray.Count} entries.");
                            break;

                        case "ZCORN":
                            var zcornData = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (deck.ZcornArray.Count < zcornData.Count) { deck.ZcornArray = zcornData; }
                            Console.WriteLine($"      [Data] Current ZCORN array status: {deck.ZcornArray.Count} entries.");
                            break;

                        case "SATNUM":
                            var satData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.SatnumArray = StoreGridArray(deck.SatnumArray, [.. satData.Select(s => s - 1)], deck, 0);
                            Console.WriteLine($"      [Data] SATNUM: Loaded and converted {deck.SatnumArray.Count} entries to 0-based indexing.");
                            break;

                        case "PVTNUM":
                            var pvtData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.PvtnumArray = StoreGridArray(deck.PvtnumArray, [.. pvtData.Select(p => p - 1)], deck, 0);
                            Console.WriteLine($"      [Data] PVTNUM: Loaded and converted {deck.PvtnumArray.Count} entries to 0-based indexing.");
                            break;

                        case "EQUALS":
                            ParseEqualsModifierBlock(tokens, ref i, deck);
                            break;

                        case "COPY":
                            ParseCopyModifierBlock(tokens, ref i, deck);
                            break;

                        case "MULTIPLY":
                            ParseMultiplyModifierBlock(tokens, ref i, deck);
                            break;

                        case "ADD":
                            ParseMultiplyModifierBlock(tokens, ref i, deck, add: true);
                            break;

                        case "START":
                            ParseStartBlock(tokens, ref i, deck);
                            break;

                        case "DATES":
                            ParseDatesBlock(tokens, ref i, deck);
                            break;

                        case "TSTEP":
                            ParseTstepBlock(tokens, ref i, deck);
                            break;

                        case "DRSDT":
                            var drsdt = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (drsdt.Count > 0)
                                Console.WriteLine($"      [Schedule] DRSDT configured: {drsdt[0]}");
                            break;

                        case "WELSPECS":
                            ParseWellSpecsBlock(tokens, ref i, deck);
                            break;

                        case "COMPDAT":
                            ParseCompdatBlock(tokens, ref i, deck);
                            break;

                        case "WCONPROD":
                            ParseWconprodBlock(tokens, ref i, deck);
                            break;

                        case "WCONINJE":
                            ParseWconinjeBlock(tokens, ref i, deck);
                            break;

                        case "PVDO":
                            deck.RawPvdoTables.AddRange(ParseRowStructuredRegions(tokens, ref i, columnsPerRow: 3));
                            Console.WriteLine($"      [Data] PVDO: Loaded {deck.RawPvdoTables.Count} region table(s).");
                            break;

                        case "PVDG":
                            deck.RawPvdgTables.AddRange(ParseRowStructuredRegions(tokens, ref i, columnsPerRow: 3));
                            Console.WriteLine($"      [Data] PVDG: Loaded {deck.RawPvdgTables.Count} region table(s).");
                            break;

                        case "PVTO":
                            deck.RawPvtoTables.AddRange(ParseNestedPvtRegions(tokens, ref i));
                            Console.WriteLine($"      [Data] PVTO: Loaded {deck.RawPvtoTables.Count} region table(s).");
                            break;

                        case "PVTG":
                            deck.RawPvtgTables.AddRange(ParseNestedPvtRegions(tokens, ref i));
                            Console.WriteLine($"      [Data] PVTG: Loaded {deck.RawPvtgTables.Count} region table(s).");
                            break;

                        case "PVTW":
                            deck.RawPvtwTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            Console.WriteLine($"      [Data] PVTW: Loaded {deck.RawPvtwTables.Count} region record(s).");
                            break;

                        case "SWOF":
                            deck.RawSwofTables.AddRange(ParseRowStructuredRegions(tokens, ref i, columnsPerRow: 4));
                            Console.WriteLine($"      [Data] SWOF: Loaded {deck.RawSwofTables.Count} region table(s).");
                            break;

                        case "SGOF":
                            deck.RawSgofTables.AddRange(ParseRowStructuredRegions(tokens, ref i, columnsPerRow: 4));
                            Console.WriteLine($"      [Data] SGOF: Loaded {deck.RawSgofTables.Count} region table(s).");
                            break;

                        case "EQUIL":
                            deck.RawEquilTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            Console.WriteLine($"      [Data] EQUIL: Loaded {deck.RawEquilTables.Count} region record(s).");
                            break;

                        case "EQUINUM":
                        case "EQLNUM":
                            var equilData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.EquinumArray = StoreGridArray(deck.EquinumArray, [.. equilData.Select(e => e - 1)], deck, 0);
                            Console.WriteLine($"      [Data] EQUINUM: Loaded and converted {deck.EquinumArray.Count} entries to 0-based indexing.");
                            break;

                        case "DENSITY":
                            var densities = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            // A missing or defaulted (1*) item takes ECLIPSE's default for the deck's units
                            var (oilDefault, waterDefault, gasDefault) = DefaultSurfaceDensities(deck.UnitSystem);
                            deck.SurfaceOilDensity = densities.Count >= 1 && densities[0] > 0 ? densities[0] : oilDefault;
                            deck.SurfaceWaterDensity = densities.Count >= 2 && densities[1] > 0 ? densities[1] : waterDefault;
                            deck.SurfaceGasDensity = densities.Count >= 3 && densities[2] > 0 ? densities[2] : gasDefault;
                            Console.WriteLine($"      [Data] DENSITY: Loaded Surface Baselines -> Oil: {deck.SurfaceOilDensity:F2} | Water: {deck.SurfaceWaterDensity:F2} | Gas: {deck.SurfaceGasDensity:F4}");
                            
                            break;

                        case "ROCK":
                            var rockProps = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (rockProps.Count >= 1)
                                deck.RockReferencePressure = rockProps[0];
                            if (rockProps.Count >= 2)
                                deck.RockCompressibility = rockProps[1];
                            Console.WriteLine($"      [Data] ROCK: Mechanical Properties -> Ref Pressure: {deck.RockReferencePressure:F2} | Compressibility: {deck.RockCompressibility:E4}");
                            break;

                        case "FAULTS":
                            ParseFaultsSection(tokens, ref i, deck.Faults);
                            break;

                        case "MULTFLT":
                            ParseMultfltSection(tokens, ref i, deck.Faults);
                            break;

                        case "RSVD":
                            deck.RawRsvdTables.AddRange(ParseRowStructuredRegions(tokens, ref i, columnsPerRow: 2));
                            Console.WriteLine($"      [Data] RSVD: Loaded {deck.RawRsvdTables.Count} region table(s).");
                            break;

                        case "FLUXNUM":
                            var fluxData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.FluxnumArray = StoreGridArray(deck.FluxnumArray, fluxData, deck, 1);
                            Console.WriteLine($"      [Data] FLUXNUM: Loaded {deck.FluxnumArray.Count} entries.");
                            break;

                        case "MULTNUM":
                            var multData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.MultnumArray = StoreGridArray(deck.MultnumArray, multData, deck, 1);
                            Console.WriteLine($"      [Data] MULTNUM: Loaded {deck.MultnumArray.Count} entries.");
                            break;

                        case "MULTREGT":
                            ParseMultregtSection(tokens, ref i, deck.MultregtRecords);
                            break;

                        case "INCLUDE":
                            HandleIncludeKeyword(tokens, ref i, currentFilePath, deck);
                            break;

                        default:
                            // A real keyword stands alone on its line; upper-case words beside other values (TOP, YES)
                            // are data of an unsupported keyword that the tokenizer also marks as keywords.
                            if (section != "SUMMARY" && !IsOutputOrSizingKeyword(keyword) && IsAloneOnLine(tokens, keywordIndex) && warnedKeywords.Add(keyword))
                                Console.WriteLine($"      [Warning] Keyword '{keyword}' is not supported and was skipped; the results do not include it.");
                            SkipUnhandledKeywordPayload(tokens, ref i);
                            break;
                    }
                }
                else
                {
                    i++;
                }
            }
        }

        // Keywords that only size arrays or control printed/written output; skipping them does not change the results
        private static bool IsOutputOrSizingKeyword(string keyword) =>
            keyword.EndsWith("DIMS") || keyword.StartsWith("RPT") ||
            keyword is "INIT" or "UNIFIN" or "UNIFOUT" or "FMTIN" or "FMTOUT" or "NOSIM" or "MESSAGES" or "NSTACK"
                    or "GRIDFILE" or "NOINSPEC" or "NORSSPEC" or "MONITOR" or "NOMONITO" or "SAVE" or "NOGGF" or "GRIDUNIT"
                    or "FAULTDIM" or "GRIDOPTS" or "MAPAXES" or "MAPUNITS" or "FIPNUM" or "WRFTPLT";

        private static bool IsAloneOnLine(List<Token> tokens, int index) =>
            (index == 0 || tokens[index - 1].Type == TokenType.NewLine) &&
            (index + 1 >= tokens.Count || tokens[index + 1].Type == TokenType.NewLine);

        private static void SkipUnhandledKeywordPayload(List<Token> tokens, ref int index)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Keyword) break;

                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    ClearLeadingLayoutSpacers(tokens, ref index);
                    if (index < tokens.Count && (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/"))
                    {
                        index++;
                        break;
                    }
                    if (index < tokens.Count && tokens[index].Type == TokenType.Keyword)
                        break;
                    continue;
                }

                index++;
            }
        }

        private static List<List<List<double>>> ParseNestedPvtRegions(List<Token> tokens, ref int index)
        {
            var regions = new List<List<List<double>>>();
            var currentRegionSubTables = new List<List<double>>();
            var currentSubTable = new List<double>();

            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;

                if (tokens[index].Type == TokenType.Keyword)
                {
                    break;
                }

                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;

                    if (currentSubTable.Count > 0)
                    {
                        currentRegionSubTables.Add(currentSubTable);
                        currentSubTable = new List<double>();
                    }
                    else
                    {
                        if (currentRegionSubTables.Count > 0)
                        {
                            regions.Add(currentRegionSubTables);
                            currentRegionSubTables = new List<List<double>>();
                        }
                    }
                    continue;
                }

                string rawVal = tokens[index].Content.ToString();
                index++;

                if (rawVal.EndsWith('*') && int.TryParse(rawVal.AsSpan(0, rawVal.Length - 1), out int defaultCount))
                {
                    for (int d = 0; d < defaultCount; d++)
                        currentSubTable.Add(0.0);
                    continue;
                }

                int starIdx = rawVal.IndexOf('*');
                if (starIdx > 0 &&
                    int.TryParse(rawVal.AsSpan(0, starIdx), out int repeatCount) &&
                    double.TryParse(rawVal.AsSpan(starIdx + 1), NumberStyles.Any, CultureInfo.InvariantCulture, out double repeatedVal))
                {
                    for (int r = 0; r < repeatCount; r++)
                        currentSubTable.Add(repeatedVal);
                }
                else if (double.TryParse(rawVal, NumberStyles.Any, CultureInfo.InvariantCulture, out double val))
                {
                    currentSubTable.Add(val);
                }
            }

            if (currentSubTable.Count > 0)
                currentRegionSubTables.Add(currentSubTable);
            if (currentRegionSubTables.Count > 0)
                regions.Add(currentRegionSubTables);

            return regions;
        }

        private static List<List<List<double>>> ParseRowStructuredRegions(List<Token> tokens, ref int index, int columnsPerRow)
        {
            var regions = new List<List<List<double>>>();
            var flatRegionValues = new List<double>();

            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;

                if (tokens[index].Type == TokenType.Keyword)
                {
                    break;
                }

                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    if (flatRegionValues.Count > 0)
                    {
                        regions.Add(ChunkIntoRows(flatRegionValues, columnsPerRow));
                        flatRegionValues = new List<double>();
                    }
                    continue;
                }

                string rawVal = tokens[index].Content.ToString();
                index++;

                if (rawVal.EndsWith('*') && int.TryParse(rawVal.AsSpan(0, rawVal.Length - 1), out int defaultCount))
                {
                    // Defaulted table values are filled in by ChunkIntoRows (linear interpolation, as in ECLIPSE)
                    for (int d = 0; d < defaultCount; d++)
                        flatRegionValues.Add(double.NaN);
                    continue;
                }

                int starIdx = rawVal.IndexOf('*');
                if (starIdx > 0 &&
                    int.TryParse(rawVal.AsSpan(0, starIdx), out int repeatCount) &&
                    double.TryParse(rawVal.AsSpan(starIdx + 1), NumberStyles.Any, CultureInfo.InvariantCulture, out double repeatedVal))
                {
                    for (int r = 0; r < repeatCount; r++)
                        flatRegionValues.Add(repeatedVal);
                }
                else if (double.TryParse(rawVal, NumberStyles.Any, CultureInfo.InvariantCulture, out double val))
                {
                    flatRegionValues.Add(val);
                }
            }

            if (flatRegionValues.Count > 0)
            {
                regions.Add(ChunkIntoRows(flatRegionValues, columnsPerRow));
            }

            return regions;
        }

        private static List<List<double>> ChunkIntoRows(List<double> flatValues, int columnsPerRow)
        {
            var rows = new List<List<double>>();
            for (int start = 0; start < flatValues.Count; start += columnsPerRow)
            {
                int count = Math.Min(columnsPerRow, flatValues.Count - start);
                rows.Add(flatValues.GetRange(start, count));
            }
            FillDefaultedTableValues(rows);
            return rows;
        }

        // A defaulted (1*) value in a table column is linearly interpolated against the first column between the nearest
        // given values above and below it, as ECLIPSE does (it used to become 0). At the ends of a column the nearest
        // given value is used; a column with no values at all becomes 0.
        private static void FillDefaultedTableValues(List<List<double>> rows)
        {
            if (rows.Count == 0) return;
            int columns = rows.Max(r => r.Count);
            for (int c = 1; c < columns; c++)
            {
                var given = Enumerable.Range(0, rows.Count).Where(r => c < rows[r].Count && !double.IsNaN(rows[r][c])).ToList();
                for (int r = 0; r < rows.Count; r++)
                {
                    if (c >= rows[r].Count || !double.IsNaN(rows[r][c])) continue;
                    int below = given.LastOrDefault(g => g < r, -1), above = given.FirstOrDefault(g => g > r, -1);
                    if (below >= 0 && above >= 0)
                    {
                        double x0 = rows[below][0], x1 = rows[above][0], x = rows[r][0];
                        double w = x1 != x0 ? (x - x0) / (x1 - x0) : 0.0;
                        rows[r][c] = rows[below][c] + w * (rows[above][c] - rows[below][c]);
                    }
                    else rows[r][c] = below >= 0 ? rows[below][c] : above >= 0 ? rows[above][c] : 0.0;
                }
            }
        }

        private static void InitializeArrays(EclipseDataDeck deck)
        {
            int totalCells = deck.NX * deck.NY * deck.NZ;
            if (totalCells <= 0) return;

            InitializeUniformArray(deck.MultXArray, totalCells, 1.0);
            InitializeUniformArray(deck.MultYArray, totalCells, 1.0);
            InitializeUniformArray(deck.MultZArray, totalCells, 1.0);
            InitializeUniformArray(deck.MultpvArray, totalCells, 1.0);
            InitializeUniformArray(deck.NtgArray, totalCells, 1);
            InitializeUniformArray(deck.ActnumArray, totalCells, 1);
            InitializeUniformArray(deck.SatnumArray, totalCells, 0);
            InitializeUniformArray(deck.PvtnumArray, totalCells, 0);
            InitializeUniformArray(deck.EquinumArray, totalCells, 0);
            Console.WriteLine($"      [Grid Initialization] MULTX, MULTY, and MULTZ baseline profiles initialized with {totalCells} entries (Default factor = 1.0).");
        }

        // START: day month year /. The record is read up to its '/', so it may be split over lines, and the month may be
        // quoted or not (before, the year was taken from a fixed token position and a split record threw).
        private void ParseStartBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            ClearLeadingLayoutSpacers(tokens, ref index);
            var record = ReadOneRecord(tokens, ref index);
            if (TryParseEclipseDate(record, out var start))
            {
                deck.StartDate = start;
                deck.CurrentSimDate = start;
                Console.WriteLine($"      [Schedule Init] Master Simulation START date initialized to: {deck.StartDate:yyyy-MM-dd}");
            }
        }

        // DATES: one date record per line, ended by a lone '/'. Without START, day 0 is ECLIPSE's default 1 JAN 1983.
        private void ParseDatesBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            foreach (var record in ReadSlashRecords(tokens, ref index))
            {
                if (!TryParseEclipseDate(record, out var nextDate)) continue;

                double dtDays = (nextDate - deck.CurrentSimDate).TotalDays;
                if (dtDays > 0)
                {
                    deck.TimeSteps.Add(dtDays);
                    deck.ReportDates.Add(nextDate);
                }

                deck.CurrentSimDate = nextDate;
                Console.WriteLine($"      [Schedule Sync] DATES stepped to: {deck.CurrentSimDate:yyyy-MM-dd} (dt = {dtDays:F1} days)");
            }
        }

        private static bool TryParseEclipseDate(List<string> record, out DateTime date)
        {
            date = default;
            if (record.Count < 3 || !int.TryParse(record[0], out int day) || !int.TryParse(record[2], out int year)) return false;
            try { date = new DateTime(year, ParseEclipseMonth(record[1]), day); return true; }
            catch (Exception ex) when (ex is FormatException or ArgumentOutOfRangeException)
            {
                Console.WriteLine($"      [Warning] Date '{string.Join(" ", record)}' could not be read; record skipped.");
                return false;
            }
        }

        private void ParseTstepBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            var steps = KeywordExtractor.ExtractDoubleArray(tokens, ref index);
            if (steps.Count == 0) return;

            double totalDays = 0.0;
            foreach (double dt in steps)
            {
                if (dt <= 0) continue;

                deck.TimeSteps.Add(dt);
                deck.CurrentSimDate = deck.CurrentSimDate.AddDays(dt);
                deck.ReportDates.Add(deck.CurrentSimDate);
                totalDays += dt;
            }

            Console.WriteLine($"      [Schedule Sync] TSTEP: Added {steps.Count} report step(s) ({totalDays:F0} days) -> End Date: {deck.CurrentSimDate:yyyy-MM-dd}");
        }

        private static int ParseEclipseMonth(string monthStr) => monthStr.Trim().Trim('\'').Trim('"').ToUpperInvariant() switch
        {
            "JAN" or "JANUARY" => 1,
            "FEB" or "FEBRUARY" => 2,
            "MAR" or "MARCH" => 3,
            "APR" or "APRIL" => 4,
            "MAY" => 5,
            "JUN" or "JNE" or "JUNE" => 6,
            "JUL" or "JLY" or "JULY" => 7,
            "AUG" or "AUGUST" => 8,
            "SEP" or "SEPT" or "SEPTEMBER" => 9,
            "OCT" or "OCTOBER" => 10,
            "NOV" or "NOVEMBER" => 11,
            "DEC" or "DECEMBER" => 12,
            _ => throw new FormatException($"Unrecognized ECLIPSE month token: '{monthStr}'.")
        };


        private void ParseWellSpecsBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    break;
                }
                if (tokens[index].Type == TokenType.Keyword) break;

                List<string> record = ReadRecordTokensUntilSlash(tokens, ref index);
                if (record.Count < 4) continue;

                string wellName = record[0];
                string groupName = record[1] != "1*" ? record[1] : "FIELD";
                int.TryParse(record[2], out int headI);
                int.TryParse(record[3], out int headJ);
                double refDepth = ParseOptionalDouble(record, 4, fallback: 0.0);
                string prefPhase = (record.Count > 5 && record[5] != "1*") ? record[5].ToUpperInvariant() : "OIL";

                if (!deck.Wells.TryGetValue(wellName, out var wellData))
                {
                    wellData = new WellData();
                    deck.Wells[wellName] = wellData;
                }

                wellData.Specification = new WellSpecRecord
                {
                    WellName = wellName,
                    GroupName = groupName,
                    HeadI = headI,
                    HeadJ = headJ,
                    ReferenceDepth = refDepth,
                    PreferredPhase = prefPhase
                };
            }
        }

        private static IEnumerable<WellData> ResolveMatchingWells(EclipseDataDeck deck, string wellNamePattern, bool createIfMissing = true)
        {
            if (wellNamePattern.EndsWith('*'))
            {
                string prefix = wellNamePattern.TrimEnd('*');
                var matches = deck.Wells
                    .Where(kvp => kvp.Key.StartsWith(prefix, StringComparison.OrdinalIgnoreCase))
                    .Select(kvp => kvp.Value)
                    .ToList();

                if (matches.Count > 0)
                    return matches;

                // A pattern that matches no well must not create a well literally named "P*"
                Console.WriteLine($"      [Warning] Well pattern '{wellNamePattern}' matches no well; record skipped.");
                return [];
            }

            if (!deck.Wells.TryGetValue(wellNamePattern, out var wellData))
            {
                if (!createIfMissing) return [];
                wellData = new WellData();
                deck.Wells[wellNamePattern] = wellData;
            }

            return [wellData];
        }

        private void ParseCompdatBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    break;
                }
                if (tokens[index].Type == TokenType.Keyword) break;

                List<string> record = ReadRecordTokensUntilSlash(tokens, ref index);
                if (record.Count < 5) continue;

                string wellPattern = record[0];

                foreach (var wellData in ResolveMatchingWells(deck, wellPattern))
                {
                    int iRaw = (record[1] == "1*" || !int.TryParse(record[1], out int parsedI) || parsedI <= 0)
                        ? wellData.Specification.HeadI
                        : parsedI;

                    int jRaw = (record[2] == "1*" || !int.TryParse(record[2], out int parsedJ) || parsedJ <= 0)
                        ? wellData.Specification.HeadJ
                        : parsedJ;

                    int.TryParse(record[3], out int kStart);
                    int.TryParse(record[4], out int kEnd);

                    string status = record.Count > 5 && record[5] != "1*" ? record[5].ToUpperInvariant() : "OPEN";
                    double transFactor = ParseOptionalDouble(record, 7, fallback: 0.0);
                    // 0 = not given; the simulator then uses its default (0.5 ft). The value is in deck units and is
                    // converted by the simulator. A 0.5 fallback here read as 0.5 m in METRIC decks.
                    double dia = ParseOptionalDouble(record, 8, fallback: 0.0);
                    double skin = ParseOptionalDouble(record, 10, fallback: 0.0);

                    int iZeroBased = iRaw - 1;
                    int jZeroBased = jRaw - 1;

                    for (int k = kStart; k <= kEnd; k++)
                    {
                        wellData.CompletionHistory.Add(new WellCompletionRecord
                        {
                            EffectiveDate = deck.CurrentSimDate,
                            I = iZeroBased,
                            J = jZeroBased,
                            KStart = k,
                            KEnd = k,
                            Status = status,
                            WellBoreDiameter = dia,
                            EffectiveSkin = skin,
                            TransmissibilityFactor = transFactor
                        });
                    }
                }
            }
        }

        private static void ClearLeadingLayoutSpacers(List<Token> tokens, ref int index)
        {
            while (index < tokens.Count &&
                  (tokens[index].Type == TokenType.NewLine ||
                   string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                   tokens[index].Content.Span.StartsWith("--")))
            {
                index++;
            }
        }

        // ---------------------------------------------------------------------------------------------------------
        // EQUALS, MULTIPLY, ADD and COPY. Each record is read up to its '/' (ReadSlashRecords), so names may be quoted
        // or not and the box may hold defaults (1*, 4*). Box items that are not given or are defaulted take the
        // current BOX, or the whole grid when there is none. Before, the box was read token by token and stopped at
        // the first 1*, so 'PERMZ' 10 1* 1* 1* 1* 3 3 / set the whole grid instead of layer 3; COPY did not strip
        // quotes (so 'PERMX' matched nothing) and ignored its box.
        // ---------------------------------------------------------------------------------------------------------

        // EQUALS: 'PROP' value [I1 I2 J1 J2 K1 K2] /
        private void ParseEqualsModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            foreach (var record in ReadSlashRecords(tokens, ref index).Select(ExpandRepeats))
            {
                if (record.Count < 2 || record[0] is not string name) continue;
                string prop = name.ToUpperInvariant();
                double value = ParseNumber(record[1]);
                var cells = ReadBox(record, 2, deck);

                if (GridArray(deck, prop, create: true) is { } array)
                {
                    int n = ForEachBoxCell(deck, cells, m => array[m] = value);
                    Console.WriteLine($"      [ARRAY MODIFICATION - EQUALS] '{prop}' = {value} in {n} cells {BoxText(cells)}");
                }
                else if (IntGridArray(deck, prop, create: true) is ({ } ints, int offset))
                {
                    int n = ForEachBoxCell(deck, cells, m => ints[m] = (int)value + offset);
                    Console.WriteLine($"      [ARRAY MODIFICATION - EQUALS] '{prop}' = {(int)value} in {n} cells {BoxText(cells)}");
                }
                else if (warnedKeywords.Add("EQUALS " + prop))
                    Console.WriteLine($"      [Warning] EQUALS: '{prop}' is not an array this reader keeps; record skipped.");
            }
        }

        // MULTIPLY (and ADD when add is true): 'PROP' value [I1 I2 J1 J2 K1 K2] /
        private void ParseMultiplyModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck, bool add = false)
        {
            string op = add ? "ADD" : "MULTIPLY";
            int totalCells = deck.NX * deck.NY * deck.NZ;
            foreach (var record in ReadSlashRecords(tokens, ref index).Select(ExpandRepeats))
            {
                if (record.Count < 2 || record[0] is not string name) continue;
                string prop = name.ToUpperInvariant();
                double value = ParseNumber(record[1]);
                var cells = ReadBox(record, 2, deck);

                if (prop == "ACTNUM" && !add && deck.ActnumArray.Count == totalCells)
                {
                    int n = ForEachBoxCell(deck, cells, m => deck.ActnumArray[m] = (int)(deck.ActnumArray[m] * value));
                    Console.WriteLine($"      [ARRAY MODIFICATION - {op}] 'ACTNUM' × {value} in {n} cells {BoxText(cells)}");
                    continue;
                }

                if (GridArray(deck, prop) is { } array && array.Count == totalCells)
                {
                    int n = ForEachBoxCell(deck, cells, m => array[m] = add ? array[m] + value : array[m] * value);
                    Console.WriteLine($"      [ARRAY MODIFICATION - {op}] '{prop}' {(add ? "+" : "×")} {value} in {n} cells {BoxText(cells)}");
                }
                else if (warnedKeywords.Add(op + " " + prop))
                    Console.WriteLine($"      [Warning] {op}: '{prop}' has not been given for every cell yet; record skipped.");
            }
        }

        // COPY: 'FROM' 'TO' [I1 I2 J1 J2 K1 K2] / copies the values of one array into another, inside the box only
        private void ParseCopyModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            int totalCells = deck.NX * deck.NY * deck.NZ;
            foreach (var record in ReadSlashRecords(tokens, ref index).Select(ExpandRepeats))
            {
                if (record.Count < 2 || record[0] is not string from || record[1] is not string to) continue;
                string source = from.ToUpperInvariant(), target = to.ToUpperInvariant();
                var cells = ReadBox(record, 2, deck);

                if (GridArray(deck, source) is { } src && src.Count == totalCells && GridArray(deck, target, create: true) is { } dst)
                {
                    int n = ForEachBoxCell(deck, cells, m => dst[m] = src[m]);
                    Console.WriteLine($"      [ARRAY MODIFICATION - COPY] '{source}' → '{target}' in {n} cells {BoxText(cells)}");
                }
                else if (IntGridArray(deck, source) is ({ } isrc, int srcOffset) && isrc.Count == totalCells &&
                         IntGridArray(deck, target, create: true) is ({ } idst, int dstOffset))
                {
                    int n = ForEachBoxCell(deck, cells, m => idst[m] = isrc[m] - srcOffset + dstOffset);
                    Console.WriteLine($"      [ARRAY MODIFICATION - COPY] '{source}' → '{target}' in {n} cells {BoxText(cells)}");
                }
                else if (warnedKeywords.Add("COPY " + source + " " + target))
                    Console.WriteLine($"      [Warning] COPY: '{source}' → '{target}' could not be copied (source not given for every cell); record skipped.");
            }
        }

        // N*value becomes N copies of value; a bare N* becomes N defaulted (null) items
        private static List<string?> ExpandRepeats(List<string> record)
        {
            var items = new List<string?>();
            foreach (var item in record)
            {
                int star = item.IndexOf('*');
                if (star > 0 && int.TryParse(item.AsSpan(0, star), out int count))
                {
                    string? value = star + 1 < item.Length ? item[(star + 1)..] : null;
                    for (int r = 0; r < count; r++) items.Add(value);
                }
                else items.Add(item);
            }
            return items;
        }

        private static double ParseNumber(string? text) =>
            double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double v) ? v : 0.0;

        // Box items I1 I2 J1 J2 K1 K2 starting at record[start]; a missing or defaulted item keeps the BOX or whole-grid value
        private (int I1, int I2, int J1, int J2, int K1, int K2) ReadBox(List<string?> record, int start, EclipseDataDeck deck)
        {
            var b = box ?? (1, deck.NX, 1, deck.NY, 1, deck.NZ);
            int[] v = [b.I1, b.I2, b.J1, b.J2, b.K1, b.K2];
            for (int n = 0; n < 6; n++)
                if (start + n < record.Count && int.TryParse(record[start + n], out int given)) v[n] = given;
            return (v[0], v[1], v[2], v[3], v[4], v[5]);
        }

        private static string BoxText((int I1, int I2, int J1, int J2, int K1, int K2) c) =>
            $"[Box I:{c.I1}-{c.I2}, J:{c.J1}-{c.J2}, K:{c.K1}-{c.K2}]";

        // Calls action(m) for every cell of the box inside the grid, m being the 0-based natural index; returns the count
        private static int ForEachBoxCell(EclipseDataDeck deck, (int I1, int I2, int J1, int J2, int K1, int K2) c, Action<int> action)
        {
            int count = 0;
            for (int k = Math.Max(1, c.K1); k <= Math.Min(deck.NZ, c.K2); k++)
                for (int j = Math.Max(1, c.J1); j <= Math.Min(deck.NY, c.J2); j++)
                    for (int i = Math.Max(1, c.I1); i <= Math.Min(deck.NX, c.I2); i++)
                    {
                        action((k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1));
                        count++;
                    }
            return count;
        }

        // The decimal grid arrays the box edits can change. With create, a missing or wrong-sized array is replaced by a
        // full one (0, or 1 for NTG and multipliers) so the edit has every cell to work on.
        private static List<double>? GridArray(EclipseDataDeck deck, string name, bool create = false)
        {
            List<double>? array = name switch
            {
                "PERMX" => deck.PermXArray, "PERMY" => deck.PermYArray, "PERMZ" => deck.PermZArray,
                "PORO" => deck.PorosityArray, "NTG" => deck.NtgArray,
                "MULTX" => deck.MultXArray, "MULTY" => deck.MultYArray, "MULTZ" => deck.MultZArray, "MULTPV" => deck.MultpvArray,
                "DX" => deck.DxArray, "DY" => deck.DyArray, "DZ" => deck.DzArray, "ZTOP" or "TOPS" => deck.ZtopArray,
                _ => null
            };
            int totalCells = deck.NX * deck.NY * deck.NZ;
            if (array is null || !create || array.Count == totalCells) return array;

            double fill = name is "NTG" or "MULTX" or "MULTY" or "MULTZ" or "MULTPV" ? 1.0 : 0.0;
            var full = new List<double>(Enumerable.Repeat(fill, totalCells));
            switch (name)
            {
                case "PERMX": deck.PermXArray = full; break;
                case "PERMY": deck.PermYArray = full; break;
                case "PERMZ": deck.PermZArray = full; break;
                case "PORO": deck.PorosityArray = full; break;
                case "NTG": deck.NtgArray = full; break;
                case "MULTX": deck.MultXArray = full; break;
                case "MULTY": deck.MultYArray = full; break;
                case "MULTZ": deck.MultZArray = full; break;
                case "MULTPV": deck.MultpvArray = full; break;
                case "DX": deck.DxArray = full; break;
                case "DY": deck.DyArray = full; break;
                case "DZ": deck.DzArray = full; break;
                default: deck.ZtopArray = full; break;
            }
            return full;
        }

        // The integer grid arrays, with the offset between the deck's value and the stored value (region numbers are
        // stored 0-based). With create, a missing array is made full size (ACTNUM 1, regions 1).
        private static (List<int>? Array, int Offset) IntGridArray(EclipseDataDeck deck, string name, bool create = false)
        {
            int totalCells = deck.NX * deck.NY * deck.NZ;
            (List<int>? array, int offset, int fill) = name switch
            {
                "ACTNUM" => (deck.ActnumArray, 0, 1),
                "SATNUM" => (deck.SatnumArray, -1, 0),
                "PVTNUM" => (deck.PvtnumArray, -1, 0),
                "EQLNUM" or "EQUINUM" => (deck.EquinumArray, -1, 0),
                "FLUXNUM" => (deck.FluxnumArray, 0, 1),
                "MULTNUM" => (deck.MultnumArray, 0, 1),
                _ => ((List<int>?)null, 0, 0)
            };
            if (array is not null && create && array.Count != totalCells)
            {
                array.Clear();
                array.AddRange(Enumerable.Repeat(fill, totalCells));
            }
            return (array, offset);
        }

        // Stores the values of an integer grid array keyword (ACTNUM, SATNUM, ...): inside a BOX only the boxed cells,
        // like the decimal arrays below (before, a boxed SATNUM replaced the whole array with a short one).
        private List<int> StoreGridArray(List<int> current, List<int> values, EclipseDataDeck deck, int fill)
        {
            if (box is not { } b) return values;

            int totalCells = deck.NX * deck.NY * deck.NZ;
            var result = current.Count == totalCells ? current : new List<int>(Enumerable.Repeat(fill, totalCells));
            int n = 0;
            for (int k = b.K1; k <= b.K2; k++)
                for (int j = b.J1; j <= b.J2; j++)
                    for (int i = b.I1; i <= b.I2 && n < values.Count; i++)
                        result[(k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1)] = values[n++];
            return result;
        }

        // Stores the values of a grid array keyword (PORO, PERMX, NTG, TOPS, ...). Inside a BOX they fill only the
        // boxed cells, I fastest, then J, then K, and the rest of the array is kept; with no BOX they replace it.
        private List<double> StoreGridArray(List<double> current, List<double> values, EclipseDataDeck deck, double fill)
        {
            if (box is not { } b) return values;

            int totalCells = deck.NX * deck.NY * deck.NZ;
            var result = current.Count == totalCells ? current : new List<double>(Enumerable.Repeat(fill, totalCells));
            int n = 0;
            for (int k = b.K1; k <= b.K2; k++)
                for (int j = b.J1; j <= b.J2; j++)
                    for (int i = b.I1; i <= b.I2 && n < values.Count; i++)
                        result[(k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1)] = values[n++];
            return result;
        }

        private static void InitializeUniformArray(List<double> array, int size, double value)
        {
            array.Clear();
            array.Capacity = size;
            for (int c = 0; c < size; c++) array.Add(value);
        }

        private static void InitializeUniformArray(List<int> array, int size, int value)
        {
            array.Clear();
            array.Capacity = size;
            for (int c = 0; c < size; c++) array.Add(value);
        }

        private void HandleIncludeKeyword(List<Token> tokens, ref int index, string baseFilePath, EclipseDataDeck deck)
        {
            if (index < tokens.Count && (tokens[index].Type == TokenType.StringLiteral || tokens[index].Type == TokenType.DataValue))
            {
                string relativePath = tokens[index].Content.ToString().Trim('\'').Trim('"');
                index++;

                string normalizedPath = relativePath.Replace('/', Path.DirectorySeparatorChar).Replace('\\', Path.DirectorySeparatorChar);

                // A relative path is looked up from the main deck's folder first (ECLIPSE/OPM), then from the folder of
                // the file holding the INCLUDE (the old behaviour), so nested includes like INCLUDE/props.inc ->
                // 'INCLUDE/swof.inc' are found either way
                string baseDir = Path.GetDirectoryName(baseFilePath) ?? "";
                string fullSubFilePath = "";
                foreach (var dir in new[] { mainDeckDirectory, baseDir }.Where(d => d.Length > 0).Distinct())
                {
                    fullSubFilePath = FindIncludeFile(dir, normalizedPath);
                    if (File.Exists(fullSubFilePath)) break;
                }

                if (File.Exists(fullSubFilePath))
                {
                    ParseSubFile(fullSubFilePath, deck);
                }
                else
                {
                    Console.BackgroundColor = ConsoleColor.DarkRed;
                    Console.ForegroundColor = ConsoleColor.White;
                    Console.WriteLine($"\n[CRITICAL WARNING] INCLUDE FILE NOT FOUND: \"{relativePath}\"");
                    Console.ResetColor();
                }
            }
            if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;
        }

        // The include file under dir, matching the file name without regard to case (deck paths often differ in case)
        private static string FindIncludeFile(string dir, string normalizedPath)
        {
            string path = Path.Combine(dir, normalizedPath);
            if (File.Exists(path)) return path;

            string fileName = Path.GetFileName(normalizedPath);
            string targetDir = Path.GetFullPath(Path.Combine(dir, Path.GetDirectoryName(normalizedPath) ?? ""));
            if (!Directory.Exists(targetDir)) return path;
            return Directory.GetFiles(targetDir).FirstOrDefault(f => Path.GetFileName(f).Equals(fileName, StringComparison.OrdinalIgnoreCase)) ?? path;
        }

        // FAULTS: 'NAME' I1 I2 J1 J2 K1 K2 FACE / ... ending with a lone '/'. A fault usually has many records,
        // so each one is added to the fault's segments (the old code replaced the fault each time, keeping only the
        // last record) and the face is kept whole ("X-" used to become 'X').
        private static void ParseFaultsSection(List<Token> tokens, ref int index, Dictionary<string, FaultRecord> faults)
        {
            int linesRead = 0;
            foreach (var record in ReadSlashRecords(tokens, ref index))
            {
                if (record.Count < 8) continue;
                string name = record[0];

                int.TryParse(record[1], out int i1);
                int.TryParse(record[2], out int i2);
                int.TryParse(record[3], out int j1);
                int.TryParse(record[4], out int j2);
                int.TryParse(record[5], out int k1);
                int.TryParse(record[6], out int k2);
                string face = record[7].ToUpperInvariant();

                if (!faults.TryGetValue(name, out var fault))
                    faults[name] = fault = new FaultRecord();
                fault.Segments.Add(new FaultSegment(i1, i2, j1, j2, k1, k2, face));
                linesRead++;
            }
            Console.WriteLine($"      [Data] FAULTS section processed: {linesRead} records in {faults.Count} faults.");
        }

        // MULTFLT: 'NAME' multiplier / ... A name ending in * matches every fault starting with it.
        private static void ParseMultfltSection(List<Token> tokens, ref int index, Dictionary<string, FaultRecord> faults)
        {
            int linesUpdated = 0;
            foreach (var record in ReadSlashRecords(tokens, ref index))
            {
                if (record.Count < 2 || !double.TryParse(record[1], CultureInfo.InvariantCulture, out double multiplier)) continue;

                string pattern = record[0];
                foreach (var (name, fault) in faults)
                {
                    bool match = pattern.EndsWith('*')
                        ? name.StartsWith(pattern.TrimEnd('*'), StringComparison.OrdinalIgnoreCase)
                        : name.Equals(pattern, StringComparison.OrdinalIgnoreCase);
                    if (!match) continue;
                    fault.TransmissibilityMultiplier = multiplier;
                    linesUpdated++;
                }
            }
            Console.WriteLine($"      [Data] MULTFLT section processed: {linesUpdated} transmissibility overrides applied.");
        }

        // Reads records of a keyword up to the lone '/' that ends it. Unlike ReadRecordTokensUntilSlash, a record may
        // hold unquoted names such as F1 or X that the tokenizer marks as keywords.
        private static List<List<string>> ReadSlashRecords(List<Token> tokens, ref int index)
        {
            var records = new List<List<string>>();
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash) { index++; break; }   // lone '/': end of keyword

                records.Add(ReadOneRecord(tokens, ref index));
            }
            return records;
        }

        // Reads one record up to its '/', across line breaks, whatever kind of token each item is; quotes are removed
        private static List<string> ReadOneRecord(List<Token> tokens, ref int index)
        {
            var record = new List<string>();
            while (index < tokens.Count && tokens[index].Type != TokenType.Slash)
            {
                if (tokens[index].Type != TokenType.NewLine)
                    record.Add(tokens[index].Content.ToString().Trim('\'').Trim('"'));
                index++;
            }
            index++;   // the record's own '/'
            return record;
        }

        // ECLIPSE's DENSITY defaults (oil, water, gas) in the deck's units: lb/ft3 (FIELD), kg/m3 (METRIC, PVT-M), g/cc (LAB)
        private static (double Oil, double Water, double Gas) DefaultSurfaceDensities(string unitSystem) =>
            unitSystem.ToUpperInvariant() switch
            {
                "METRIC" or "PVT-M" => (600.0, 999.014, 1.0),
                "LAB" => (0.6, 0.999014, 0.001),
                _ => (37.457, 62.366, 0.062428)
            };

        // One atmosphere in the deck's pressure units: ECLIPSE's default producer BHP limit (WCONPROD item 9)
        private static double OneAtmosphere(string unitSystem) =>
            unitSystem.ToUpperInvariant() switch
            {
                "METRIC" => 1.01325,
                "LAB" or "PVT-M" => 1.0,
                _ => 14.69595
            };

        private static List<string> ReadRecordTokensUntilSlash(List<Token> tokens, ref int index)
        {
            var record = new List<string>();

            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;

                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    break;
                }

                if (tokens[index].Type == TokenType.Keyword)
                {
                    break;
                }

                string rawVal = tokens[index].Content.ToString().Trim('\'').Trim('"');
                index++;

                if (rawVal.EndsWith("*") && int.TryParse(rawVal.TrimEnd('*'), out int repeatCount))
                {
                    for (int r = 0; r < repeatCount; r++) record.Add("1*");
                }
                else if (rawVal.Contains("*"))
                {
                    var parts = rawVal.Split('*');
                    if (parts.Length == 2 && int.TryParse(parts[0], out int rep))
                    {
                        for (int r = 0; r < rep; r++) record.Add(parts[1]);
                    }
                    else
                    {
                        record.Add(rawVal);
                    }
                }
                else
                {
                    record.Add(rawVal);
                }
            }

            return record;
        }

        private static double ParseOptionalDouble(List<string> items, int pos, double fallback = 0.0)
        {
            if (pos >= items.Count) return fallback;
            string val = items[pos];
            if (val == "1*" || string.IsNullOrWhiteSpace(val)) return fallback;
            return double.TryParse(val, NumberStyles.Any, CultureInfo.InvariantCulture, out double parsed) ? parsed : fallback;
        }

        private void ParseWconprodBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    break;
                }
                if (tokens[index].Type == TokenType.Keyword) break;

                List<string> record = ReadRecordTokensUntilSlash(tokens, ref index);
                if (record.Count == 0) continue;

                string wellPattern = record[0];
                string wellStatus = record.Count > 1 && record[1] != "1*" ? record[1].ToUpperInvariant() : "OPEN";
                string mode = record.Count > 2 && record[2] != "1*" ? record[2].ToUpperInvariant() : "BHP";

                double oRat = ParseOptionalDouble(record, 3);
                double wRat = ParseOptionalDouble(record, 4);
                double gRat = ParseOptionalDouble(record, 5);
                double lRat = ParseOptionalDouble(record, 6);
                double resRat = ParseOptionalDouble(record, 7);
                // A defaulted BHP limit is one atmosphere in ECLIPSE; 0 read as "no limit" and let pressure fall below zero
                double bhp = ParseOptionalDouble(record, 8, fallback: OneAtmosphere(deck.UnitSystem));

                foreach (var wellData in ResolveMatchingWells(deck, wellPattern))
                {
                    wellData.ProducerControlHistory.Add(new ProducerControlRecord
                    {
                        EffectiveDate = deck.CurrentSimDate,
                        Status = wellStatus,
                        ControlMode = mode,
                        OilRateTarget = oRat,
                        WaterRateTarget = wRat,
                        GasRateTarget = gRat,
                        LiquidRateTarget = lRat,
                        ResvRateTarget = resRat,
                        MinBhpLimit = bhp
                    });

                    Console.WriteLine($"      [Schedule] WCONPROD parsed for well '{wellData.Specification.WellName}': Status={wellStatus}, Mode={mode}, BHP={bhp}");
                }
            }
        }

        private void ParseWconinjeBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    break;
                }
                if (tokens[index].Type == TokenType.Keyword) break;

                List<string> record = ReadRecordTokensUntilSlash(tokens, ref index);
                if (record.Count == 0) continue;

                string wellPattern = record[0];
                string phaseType = record.Count > 1 && record[1] != "1*" ? record[1].ToUpperInvariant() : "WATER";
                string wellStatus = record.Count > 2 && record[2] != "1*" ? record[2].ToUpperInvariant() : "OPEN";
                string mode = record.Count > 3 && record[3] != "1*" ? record[3].ToUpperInvariant() : "RATE";
                double rate = ParseOptionalDouble(record, 4);
                double resvRate = ParseOptionalDouble(record, 5);
                double maxBhp = ParseOptionalDouble(record, 6);

                foreach (var wellData in ResolveMatchingWells(deck, wellPattern))
                {
                    wellData.InjectorControlHistory.Add(new InjectorControlRecord
                    {
                        EffectiveDate = deck.CurrentSimDate,
                        InjectionPhase = phaseType,
                        Status = wellStatus,
                        ControlMode = mode,
                        TargetRate = rate,
                        ResvRateTarget = resvRate,
                        MaxBhpLimit = maxBhp
                    });

                    Console.WriteLine($"      [Schedule] WCONINJE parsed for well '{wellData.Specification.WellName}': Phase={phaseType}, Status={wellStatus}, Mode={mode}, Rate={rate}, MaxBHP={maxBhp}");
                }
            }
        }

        private static void ParseMultregtSection(List<Token> tokens, ref int index, List<MultregtRecord> multregtRecords)
        {
            int countBefore = multregtRecords.Count;

            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count) break;

                if (tokens[index].Type == TokenType.Slash || tokens[index].Content.ToString() == "/")
                {
                    index++;
                    break;
                }

                if (tokens[index].Type == TokenType.Keyword)
                    break;

                List<string> record = ReadRecordTokensUntilSlash(tokens, ref index);
                if (record.Count < 3) continue;

                int.TryParse(record[0], out int fromRegion);
                int.TryParse(record[1], out int toRegion);
                double multiplier = ParseOptionalDouble(record, 2, fallback: 1.0);

                string direction = (record.Count > 3 && record[3] != "1*") ? record[3].ToUpperInvariant() : "XYZ";
                string connType = (record.Count > 4 && record[4] != "1*") ? record[4].ToUpperInvariant() : "ALL";
                string regFlag = (record.Count > 5 && record[5] != "1*") ? record[5].ToUpperInvariant() : "M";

                multregtRecords.Add(new MultregtRecord(
                    fromRegion,
                    toRegion,
                    multiplier,
                    direction,
                    connType,
                    regFlag));
            }

            Console.WriteLine($"      [Data] MULTREGT: Loaded {multregtRecords.Count - countBefore} region transmissibility multiplier(s).");
        }

    }
}
