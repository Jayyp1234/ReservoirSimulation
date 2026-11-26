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

        public EclipseDataDeck LoadDeckFile(string filepath)
        {
            // Reading state starts afresh for every deck, in case this reader is used again
            box = null; section = "RUNSPEC"; warnedKeywords.Clear();

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
                            deck.ActnumArray = actData;
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
                            deck.SatnumArray = [.. satData.Select(s => s - 1)];
                            Console.WriteLine($"      [Data] SATNUM: Loaded and converted {deck.SatnumArray.Count} entries to 0-based indexing.");
                            break;

                        case "PVTNUM":
                            var pvtData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.PvtnumArray = [.. pvtData.Select(p => p - 1)];
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
                            ParseDissolutionRateKeyword(tokens, ref i, deck, isVaporization: false);
                            break;

                        case "DRVDT":
                            ParseDissolutionRateKeyword(tokens, ref i, deck, isVaporization: true);
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
                            deck.EquinumArray = [.. equilData.Select(e => e - 1)];
                            Console.WriteLine($"      [Data] EQUINUM: Loaded and converted {deck.EquinumArray.Count} entries to 0-based indexing.");
                            break;

                        case "DENSITY":
                            var densities = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (densities.Count >= 1)
                                deck.SurfaceOilDensity = densities[0];
                            if (densities.Count >= 2)
                                deck.SurfaceWaterDensity = densities[1];
                            if (densities.Count >= 3)
                                deck.SurfaceGasDensity = densities[2];
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

                        case "RVVD":
                            deck.RawRvvdTables.AddRange(ParseRowStructuredRegions(tokens, ref i, columnsPerRow: 2));
                            Console.WriteLine($"      [Data] RVVD: Loaded {deck.RawRvvdTables.Count} region table(s).");
                            break;

                        case "FLUXNUM":
                            var fluxData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.FluxnumArray = fluxData;
                            Console.WriteLine($"      [Data] FLUXNUM: Loaded {deck.FluxnumArray.Count} entries.");
                            break;

                        case "MULTNUM":
                            var multData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            deck.MultnumArray = multData;
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

        private static void ParseDissolutionRateKeyword(List<Token> tokens, ref int index, EclipseDataDeck deck, bool isVaporization)
        {
            List<string> record = ReadRecordTokensUntilSlash(tokens, ref index);

            // Item 1: Max rate of Rs/Rv increase (defaulted 1* or omitted => PositiveInfinity = unlimited)
            double maxRate = double.PositiveInfinity;
            if (record.Count > 0 && record[0] != "1*" && !string.IsNullOrWhiteSpace(record[0]) &&
                double.TryParse(record[0], NumberStyles.Any, CultureInfo.InvariantCulture, out double parsedRate))
            {
                maxRate = Math.Max(0.0, parsedRate);
            }

            // Item 2: Option flag ('ALL' [default] or 'FREE')
            string option = (record.Count > 1 && record[1] != "1*" && !string.IsNullOrWhiteSpace(record[1]))
                ? record[1].ToUpperInvariant()
                : "ALL";

            int stepIdx = deck.TimeSteps.Count;
            string kwName = isVaporization ? "DRVDT" : "DRSDT";
            string rateDisplay = double.IsInfinity(maxRate) ? "UNLIMITED" : maxRate.ToString(CultureInfo.InvariantCulture);

            if (isVaporization)
            {
                deck.DrvdtMax = maxRate;
                deck.DrvdtOption = option;
                deck.DrvdtSchedule.Add(new DissolutionRateRecord(stepIdx, maxRate, option));
            }
            else
            {
                deck.DrsdtMax = maxRate;
                deck.DrsdtOption = option;
                deck.DrsdtSchedule.Add(new DissolutionRateRecord(stepIdx, maxRate, option));
            }

            Console.WriteLine($"      [Schedule] {kwName} configured at step {stepIdx}: Rate={rateDisplay}, Option={option}");
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

        // True when the token starts the next keyword: the tokenizer marks every short upper-case word as a keyword,
        // so unquoted well names and items inside a record (GI01, OP01, OPEN, GAS, BHP) look like keywords too. A real
        // keyword has nothing after it on its line but a comment; a record item always has other items or its '/'.
        // Checking only the token type made WELSPECS, COMPDAT and WCONPROD/WCONINJE stop at the first unquoted name or
        // word, so decks that do not quote them (SPE10, most field decks) lost their wells without a warning.
        private static bool IsKeywordLine(List<Token> tokens, int index)
        {
            if (index >= tokens.Count || tokens[index].Type != TokenType.Keyword) return false;
            for (int j = index + 1; j < tokens.Count && tokens[j].Type != TokenType.NewLine; j++)
            {
                if (tokens[j].Content.Span.StartsWith("--")) break;              // the rest of the line is a comment
                if (!string.IsNullOrWhiteSpace(tokens[j].Content.ToString())) return false;
            }
            return true;
        }

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
                    for (int d = 0; d < defaultCount; d++)
                        flatRegionValues.Add(0.0);
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
            return rows;
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

        private void ParseStartBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            ClearLeadingLayoutSpacers(tokens, ref index);
            if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; return; }

            if (int.TryParse(tokens[index].Content.ToString(), out int day))
            {
                string monthStr = tokens[index + 1].Content.ToString().Trim('\'').Trim('"').ToUpperInvariant();
                int.TryParse(tokens[index + 2].Content.ToString(), out int year);

                int month = monthStr switch
                {
                    "JAN" => 1,
                    "FEB" => 2,
                    "MAR" => 3,
                    "APR" => 4,
                    "MAY" => 5,
                    "JUN" or "JNE" => 6,
                    "JUL" or "JLY" => 7,
                    "AUG" => 8,
                    "SEP" => 9,
                    "OCT" => 10,
                    "NOV" => 11,
                    "DEC" => 12,
                    _ => 1
                };

                deck.StartDate = new DateTime(year, month, day);
                deck.CurrentSimDate = deck.StartDate;

                Console.WriteLine($"      [Schedule Init] Master Simulation START date initialized to: {deck.StartDate:yyyy-MM-dd}");
            }

            index += 3;
            AdvanceAndClearTrailingJunk(tokens, ref index);
        }

        private void ParseDatesBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; break; }
                if (tokens[index].Type == TokenType.Keyword) break;

                if (int.TryParse(tokens[index].Content.ToString(), out int day))
                {
                    string monthStr = tokens[index + 1].Content.ToString().Trim('\'').Trim('"').ToUpperInvariant();
                    int.TryParse(tokens[index + 2].Content.ToString(), out int year);

                    int month = ParseEclipseMonth(monthStr);
                    var nextDate = new DateTime(year, month, day);

                    if (deck.StartDate == default)
                        deck.StartDate = nextDate;

                    double dtDays = (nextDate - deck.CurrentSimDate).TotalDays;
                    if (dtDays > 0)
                    {
                        deck.TimeSteps.Add(dtDays);
                        deck.ReportDates.Add(nextDate);
                    }

                    deck.CurrentSimDate = nextDate;
                    Console.WriteLine($"      [Schedule Sync] DATES stepped to: {deck.CurrentSimDate:yyyy-MM-dd} (dt = {dtDays:F1} days)");
                }
                index += 3;
                AdvanceAndClearTrailingJunk(tokens, ref index);
            }
        }

        private void ParseTstepBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            var steps = KeywordExtractor.ExtractDoubleArray(tokens, ref index);
            if (steps.Count == 0) return;

            // Anchor Day 0 cleanly if START was omitted (e.g., SPE1 Odeh deck)
            if (deck.StartDate == default)
            {
                deck.StartDate = new DateTime(2000, 1, 1);
                deck.CurrentSimDate = deck.StartDate;
            }

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
                if (IsKeywordLine(tokens, index)) break;

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
                if (IsKeywordLine(tokens, index)) break;

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
                    double dia = ParseOptionalDouble(record, 8, fallback: 0.5);
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

        private static void AdvanceAndClearTrailingJunk(List<Token> tokens, ref int index)
        {
            while (index < tokens.Count && tokens[index].Type != TokenType.Slash)
            {
                index++;
            }
            if (index < tokens.Count && tokens[index].Type == TokenType.Slash)
            {
                index++;
            }
            while (index < tokens.Count)
            {
                if (tokens[index].Type == TokenType.NewLine)
                {
                    index++;
                    break;
                }
                index++;
            }
        }

        // MULTIPLY (and ADD when add is true): 'PROP' value [I1 I2 J1 J2 K1 K2] /
        private void ParseMultiplyModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck, bool add = false)
        {
            string op = add ? "ADD" : "MULTIPLY";
            while (index < tokens.Count)
            {
                while (index < tokens.Count &&
                      (tokens[index].Type == TokenType.NewLine ||
                       string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                       tokens[index].Content.Span.StartsWith("--")))
                {
                    index++;
                }

                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash) { index++; break; }

                string targetProp = tokens[index].Content.ToString().ToUpperInvariant().Trim('\'').Trim('"');

                if (tokens[index].Type == TokenType.Keyword &&
                    targetProp != "PERMX" && targetProp != "PERMY" && targetProp != "PERMZ" && targetProp != "PORO" &&
                    targetProp != "NTG" && targetProp != "ACTNUM" && targetProp != "MULTX" && targetProp != "MULTY" && targetProp != "MULTZ" &&
                    targetProp != "MULTPV" && targetProp != "DX" && targetProp != "DY" && targetProp != "DZ" && targetProp != "ZTOP" && targetProp != "TOPS")
                {
                    break;
                }
                index++;

                if (index >= tokens.Count) break;
                double.TryParse(tokens[index].Content.ToString(), CultureInfo.InvariantCulture, out double multiplierValue);
                index++;

                // The six box numbers are optional: missing ones default to the current BOX, or to the whole grid
                var (iStart, iEnd, jStart, jEnd, kStart, kEnd) = box ?? (1, deck.NX, 1, deck.NY, 1, deck.NZ);
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int iS)) { iStart = iS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int iE)) { iEnd = iE; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int jS)) { jStart = jS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int jE)) { jEnd = jE; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int kS)) { kStart = kS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int kE)) { kEnd = kE; index++; }

                while (index < tokens.Count && tokens[index].Type != TokenType.Slash) index++;
                if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;

                while (index < tokens.Count)
                {
                    if (tokens[index].Type == TokenType.NewLine) { index++; break; }
                    index++;
                }

                int totalCells = deck.NX * deck.NY * deck.NZ;
                if (targetProp == "ACTNUM" && !add && deck.ActnumArray.Count == totalCells)
                {
                    int modificationCount = 0;
                    for (int k = kStart; k <= kEnd; k++)
                    {
                        for (int j = jStart; j <= jEnd; j++)
                        {
                            for (int i = iStart; i <= iEnd; i++)
                            {
                                int flatIndex = (k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1);
                                if (flatIndex >= 0 && flatIndex < totalCells)
                                {
                                    deck.ActnumArray[flatIndex] = (int)(deck.ActnumArray[flatIndex] * multiplierValue);
                                    modificationCount++;
                                }
                            }
                        }
                    }
                    Console.WriteLine($"      [ARRAY MODIFICATION - MULTIPLY] Modified '{targetProp}': Multiplied {modificationCount} cells by {multiplierValue} [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
                    continue;
                }

                List<double> targetArray = targetProp switch
                {
                    "PERMX" => deck.PermXArray,
                    "PERMY" => deck.PermYArray,
                    "PERMZ" => deck.PermZArray,
                    "PORO" => deck.PorosityArray,
                    "NTG" => deck.NtgArray,
                    "MULTX" => deck.MultXArray,
                    "MULTY" => deck.MultYArray,
                    "MULTZ" => deck.MultZArray,
                    "MULTPV" => deck.MultpvArray,
                    "DX" => deck.DxArray,
                    "DY" => deck.DyArray,
                    "DZ" => deck.DzArray,
                    "ZTOP" or "TOPS" => deck.ZtopArray,
                    _ => null
                };

                if (targetArray != null && targetArray.Count == totalCells)
                {
                    int modificationCount = 0;
                    for (int k = kStart; k <= kEnd; k++)
                    {
                        for (int j = jStart; j <= jEnd; j++)
                        {
                            for (int i = iStart; i <= iEnd; i++)
                            {
                                int flatIndex = (k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1);
                                if (flatIndex >= 0 && flatIndex < totalCells)
                                {
                                    if (add) targetArray[flatIndex] += multiplierValue;
                                    else targetArray[flatIndex] *= multiplierValue;
                                    modificationCount++;
                                }
                            }
                        }
                    }
                    Console.WriteLine($"      [ARRAY MODIFICATION - {op}] Modified '{targetProp}': {(add ? "Added" : "Multiplied")} {modificationCount} cells by {multiplierValue} [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
                }
            }
        }

        private void ParseCopyModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count && tokens[index].Type != TokenType.Slash)
            {
                if (tokens[index].Type == TokenType.NewLine ||
                    string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                    tokens[index].Content.Span.StartsWith("--"))
                {
                    index++;
                    continue;
                }

                if (index + 1 >= tokens.Count) break;

                string sourceProp = tokens[index].Content.ToString().ToUpperInvariant();
                string targetProp = tokens[index + 1].Content.ToString().ToUpperInvariant();
                index += 2;

                while (index < tokens.Count && tokens[index].Type != TokenType.Slash && tokens[index].Type != TokenType.Keyword) index++;
                if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;

                if (sourceProp == "ACTNUM" && deck.ActnumArray.Count > 0)
                {
                    var clonedArray = new List<int>(deck.ActnumArray);
                    if (targetProp == "ACTNUM") deck.ActnumArray = clonedArray;
                    Console.WriteLine($"      [ARRAY MODIFICATION - COPY] Modified '{targetProp}': Overwrote array values by cloning all {clonedArray.Count} elements from '{sourceProp}'");
                    continue;
                }

                List<double> sourceArray = sourceProp switch
                {
                    "PERMX" => deck.PermXArray,
                    "PERMY" => deck.PermYArray,
                    "PERMZ" => deck.PermZArray,
                    "PORO" => deck.PorosityArray,
                    "NTG" => deck.NtgArray,
                    "MULTX" => deck.MultXArray,
                    "MULTY" => deck.MultYArray,
                    "MULTZ" => deck.MultZArray,
                    "DX" => deck.DxArray,
                    "DY" => deck.DyArray,
                    "DZ" => deck.DzArray,
                    "ZTOP" or "TOPS" => deck.ZtopArray,
                    _ => null
                };

                if (sourceArray != null && sourceArray.Count > 0)
                {
                    var clonedArray = new List<double>(sourceArray);
                    switch (targetProp)
                    {
                        case "PERMX": deck.PermXArray = clonedArray; break;
                        case "PERMY": deck.PermYArray = clonedArray; break;
                        case "PERMZ": deck.PermZArray = clonedArray; break;
                        case "PORO": deck.PorosityArray = clonedArray; break;
                        case "NTG": deck.NtgArray = clonedArray; break;
                        case "MULTX": deck.MultXArray = clonedArray; break;
                        case "MULTY": deck.MultYArray = clonedArray; break;
                        case "MULTZ": deck.MultZArray = clonedArray; break;
                        case "DX": deck.DxArray = clonedArray; break;
                        case "DY": deck.DyArray = clonedArray; break;
                        case "DZ": deck.DzArray = clonedArray; break;
                        case "ZTOP":
                        case "TOPS": deck.ZtopArray = clonedArray; break;
                    }
                    Console.WriteLine($"      [ARRAY MODIFICATION - COPY] Modified '{targetProp}': Overwrote array values by cloning all {clonedArray.Count} elements from '{sourceProp}'");
                }
            }
            if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;
        }

        private void ParseEqualsModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                while (index < tokens.Count &&
                      (tokens[index].Type == TokenType.NewLine ||
                       string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                       tokens[index].Content.Span.StartsWith("--")))
                {
                    index++;
                }

                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash) { index++; break; }

                string targetProp = tokens[index].Content.ToString().ToUpperInvariant().Trim('\'').Trim('"');

                if (tokens[index].Type == TokenType.Keyword &&
                    targetProp != "PERMX" && targetProp != "PERMY" && targetProp != "PERMZ" && targetProp != "PORO" &&
                    targetProp != "NTG" && targetProp != "ACTNUM" && targetProp != "MULTX" && targetProp != "MULTY" && targetProp != "MULTZ" &&
                    targetProp != "MULTPV" && targetProp != "DX" && targetProp != "DY" && targetProp != "DZ" && targetProp != "ZTOP" && targetProp != "TOPS")
                {
                    break;
                }
                index++;

                if (index >= tokens.Count) break;

                if (!double.TryParse(tokens[index].Content.ToString(), CultureInfo.InvariantCulture, out double value))
                {
                    value = 0.0;
                }
                index++;

                // Cells not given in the record default to the current BOX, or to the whole grid when there is none
                var (iStart, iEnd, jStart, jEnd, kStart, kEnd) = box ?? (1, deck.NX, 1, deck.NY, 1, deck.NZ);

                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int iS)) { iStart = iS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int iE)) { iEnd = iE; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int jS)) { jStart = jS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int jE)) { jEnd = jE; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int kS)) { kStart = kS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int kE)) { kEnd = kE; index++; }

                while (index < tokens.Count && tokens[index].Type != TokenType.Slash) index++;
                if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;

                while (index < tokens.Count)
                {
                    if (tokens[index].Type == TokenType.NewLine) { index++; break; }
                    index++;
                }

                int totalCells = deck.NX * deck.NY * deck.NZ;
                if (totalCells <= 0) continue;

                if (targetProp == "MULTX" || targetProp == "MULTY" || targetProp == "MULTZ" || targetProp == "MULTPV" ||
                    targetProp == "PERMX" || targetProp == "PERMY" || targetProp == "PERMZ" ||
                    targetProp == "PORO" || targetProp == "NTG" || targetProp == "DX" ||
                    targetProp == "DY" || targetProp == "DZ" || targetProp == "ZTOP" || targetProp == "TOPS")
                {
                    List<double> targetArray = targetProp switch
                    {
                        "PERMX" => deck.PermXArray,
                        "PERMY" => deck.PermYArray,
                        "PERMZ" => deck.PermZArray,
                        "PORO" => deck.PorosityArray,
                        "NTG" => deck.NtgArray,
                        "MULTX" => deck.MultXArray,
                        "MULTY" => deck.MultYArray,
                        "MULTZ" => deck.MultZArray,
                        "MULTPV" => deck.MultpvArray,
                        "DX" => deck.DxArray,
                        "DY" => deck.DyArray,
                        "DZ" => deck.DzArray,
                        "ZTOP" or "TOPS" => deck.ZtopArray,
                        _ => null
                    };

                    if (targetArray != null)
                    {
                        if (targetArray.Count == 0) InitializeUniformArray(targetArray, totalCells, 1.0);

                        int modificationCount = 0; int flatk, flatj, flati;
                        for (int k = kStart; k <= kEnd; k++)
                        {
                            flatk = k - 1;
                            for (int j = jStart; j <= jEnd; j++)
                            {
                                flatj = flatk * deck.NY + (j - 1);
                                for (int i = iStart; i <= iEnd; i++)
                                {
                                    flati = flatj * deck.NX + (i - 1);
                                    if (flati >= 0 && flati < totalCells)
                                    {
                                        targetArray[flati] = value;
                                        modificationCount++;
                                    }
                                }
                            }
                        }
                        Console.WriteLine($"      [ARRAY MODIFICATION - EQUALS] Modified '{targetProp}': Assigned value {value} across {modificationCount} cells [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
                    }
                }
                else
                {
                    int adjustedValue = (int)value;

                    switch (targetProp)
                    {
                        case "ACTNUM":
                            if (deck.ActnumArray.Count == 0) InitializeUniformArray(deck.ActnumArray, totalCells, 1);
                            ApplyIntEqualsOverride(deck.ActnumArray, adjustedValue, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "SATNUM":
                            if (deck.SatnumArray.Count == 0) InitializeUniformArray(deck.SatnumArray, totalCells, 0);
                            ApplyIntEqualsOverride(deck.SatnumArray, adjustedValue - 1, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "PVTNUM":
                            if (deck.PvtnumArray.Count == 0) InitializeUniformArray(deck.PvtnumArray, totalCells, 0);
                            ApplyIntEqualsOverride(deck.PvtnumArray, adjustedValue - 1, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "EQUINUM":
                            if (deck.EquinumArray.Count == 0) InitializeUniformArray(deck.EquinumArray, totalCells, 0);
                            ApplyIntEqualsOverride(deck.EquinumArray, adjustedValue - 1, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "FLUXNUM":
                            if (deck.FluxnumArray.Count == 0) InitializeUniformArray(deck.FluxnumArray, totalCells, 1);
                            ApplyIntEqualsOverride(deck.FluxnumArray, adjustedValue, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "MULTNUM":
                            if (deck.MultnumArray.Count == 0) InitializeUniformArray(deck.MultnumArray, totalCells, 1);
                            ApplyIntEqualsOverride(deck.MultnumArray, adjustedValue, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;
                    }
                }
            }
        }

        private static void ApplyIntEqualsOverride(List<int> array, double value, int iStart, int iEnd, int jStart, int jEnd, int kStart, int kEnd, int nx, int ny, int totalCells, string propName)
        {
            int modificationCount = 0;
            for (int k = kStart; k <= kEnd; k++)
            {
                for (int j = jStart; j <= jEnd; j++)
                {
                    for (int i = iStart; i <= iEnd; i++)
                    {
                        int flatIndex = (k - 1) * nx * ny + (j - 1) * nx + (i - 1);
                        if (flatIndex >= 0 && flatIndex < totalCells)
                        {
                            array[flatIndex] = (int)value;
                            modificationCount++;
                        }
                    }
                }
            }
            Console.WriteLine($"      [ARRAY MODIFICATION - EQUALS] Modified '{propName}': Assigned region index {(int)value} across {modificationCount} cells [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
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

                string baseDir = Path.GetDirectoryName(baseFilePath) ?? "";
                string normalizedPath = relativePath.Replace('/', Path.DirectorySeparatorChar).Replace('\\', Path.DirectorySeparatorChar);
                string fullSubFilePath = Path.Combine(baseDir, normalizedPath);

                if (!File.Exists(fullSubFilePath))
                {
                    string fileName = Path.GetFileName(normalizedPath);
                    string subDirRel = Path.GetDirectoryName(normalizedPath) ?? "";
                    string targetDir = Path.GetFullPath(Path.Combine(baseDir, subDirRel));

                    if (Directory.Exists(targetDir))
                    {
                        string matchedFile = Directory.GetFiles(targetDir)
                            .FirstOrDefault(f => Path.GetFileName(f).Equals(fileName, StringComparison.OrdinalIgnoreCase));

                        if (matchedFile != null) fullSubFilePath = matchedFile;
                    }
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

                var record = new List<string>();
                while (index < tokens.Count && tokens[index].Type != TokenType.Slash)
                {
                    if (tokens[index].Type != TokenType.NewLine)
                        record.Add(tokens[index].Content.ToString().Trim('\'').Trim('"'));
                    index++;
                }
                index++;   // the record's '/'
                records.Add(record);
            }
            return records;
        }

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

                if (IsKeywordLine(tokens, index))
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
                if (IsKeywordLine(tokens, index)) break;

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
                double bhp = ParseOptionalDouble(record, 8);

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
                if (IsKeywordLine(tokens, index)) break;

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

                if (IsKeywordLine(tokens, index))
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