using System;
using System.Collections.Generic;
using System.Globalization;

namespace EclipseDeckReader
{
    public static class KeywordExtractor
    {
        public static List<List<double>> ExtractMultiTableGroup(List<Token> tokens, ref int tokenIndex, bool splitOnSingleSlash)
        {
            var tablesList = new List<List<double>>();
            var currentTable = new List<double>(4096);
            bool lastWasSlash = false;

            while (tokenIndex < tokens.Count && tokens[tokenIndex].Type != TokenType.Keyword)
            {
                string rawTokenStr = tokens[tokenIndex].Content.ToString();

                // Skip layout decoration whitespace, comments, and explicit lines
                if (string.IsNullOrWhiteSpace(rawTokenStr) ||
                    rawTokenStr.StartsWith("--") ||
                    rawTokenStr.StartsWith("\r") ||
                    rawTokenStr.StartsWith("\n"))
                {
                    tokenIndex++;
                    continue;
                }

                if (tokens[tokenIndex].Type == TokenType.Slash)
                {
                    if (splitOnSingleSlash)
                    {
                        if (currentTable.Count > 0)
                        {
                            tablesList.Add(currentTable);
                            currentTable = new List<double>(4096);
                        }
                        tokenIndex++;
                        continue;
                    }

                    if (lastWasSlash || currentTable.Count == 0)
                    {
                        if (currentTable.Count > 0)
                        {
                            tablesList.Add(currentTable);
                            currentTable = new List<double>(4096);
                        }
                        lastWasSlash = true;
                        tokenIndex++;
                        continue;
                    }

                    lastWasSlash = true;
                    tokenIndex++;
                    continue;
                }

                lastWasSlash = false;
                int starIdx = rawTokenStr.IndexOf('*');

                if (starIdx > 0)
                {
                    if (int.TryParse(rawTokenStr.AsSpan(0, starIdx), out int count))
                    {
                        // N*value repeats the value; a bare N* is N defaulted items, kept as 0 so later columns stay in place
                        double val = double.TryParse(rawTokenStr.AsSpan(starIdx + 1), CultureInfo.InvariantCulture, out double v) ? v : 0.0;
                        for (int i = 0; i < count; i++) currentTable.Add(val);
                    }
                }
                else
                {
                    if (double.TryParse(rawTokenStr, CultureInfo.InvariantCulture, out double val))
                    {
                        currentTable.Add(val);
                    }
                }
                tokenIndex++;
            }

            if (currentTable.Count > 0) tablesList.Add(currentTable);
            return tablesList;
        }

        public static List<double> ExtractDoubleArray(List<Token> tokens, ref int tokenIndex)
        {
            var dataValues = new List<double>(16384);

            while (tokenIndex < tokens.Count && tokens[tokenIndex].Type != TokenType.Keyword)
            {
                string rawTokenStr = tokens[tokenIndex].Content.ToString();

                if (string.IsNullOrWhiteSpace(rawTokenStr) ||
                    rawTokenStr.StartsWith("--") ||
                    rawTokenStr.StartsWith("\r") ||
                    rawTokenStr.StartsWith("\n"))
                {
                    tokenIndex++;
                    continue;
                }

                if (tokens[tokenIndex].Type == TokenType.Slash)
                {
                    tokenIndex++;
                    continue;
                }

                int starIdx = rawTokenStr.IndexOf('*');
                if (starIdx > 0)
                {
                    // N*value repeats the value; a bare N* (defaulted items) used to throw a FormatException
                    int count = int.Parse(rawTokenStr.AsSpan(0, starIdx));
                    double val = double.TryParse(rawTokenStr.AsSpan(starIdx + 1), CultureInfo.InvariantCulture, out double v) ? v : 0.0;
                    for (int i = 0; i < count; i++) dataValues.Add(val);
                }
                else
                {
                    if (double.TryParse(rawTokenStr, CultureInfo.InvariantCulture, out double val))
                        dataValues.Add(val);
                }
                tokenIndex++;
            }
            return dataValues;
        }

        public static List<int> ExtractIntArray(List<Token> tokens, ref int tokenIndex)
        {
            var dataValues = new List<int>(16384);

            while (tokenIndex < tokens.Count && tokens[tokenIndex].Type != TokenType.Keyword)
            {
                string rawTokenStr = tokens[tokenIndex].Content.ToString();

                if (string.IsNullOrWhiteSpace(rawTokenStr) ||
                    rawTokenStr.StartsWith("--") ||
                    rawTokenStr.StartsWith("\r") ||
                    rawTokenStr.StartsWith("\n"))
                {
                    tokenIndex++;
                    continue;
                }

                if (tokens[tokenIndex].Type == TokenType.Slash)
                {
                    tokenIndex++;
                    continue;
                }

                int starIdx = rawTokenStr.IndexOf('*');
                if (starIdx > 0)
                {
                    // A bare N* (defaulted) used to throw; integer grid arrays (ACTNUM, SATNUM, ...) default to 1
                    int count = int.Parse(rawTokenStr.AsSpan(0, starIdx));
                    int val = int.TryParse(rawTokenStr.AsSpan(starIdx + 1), out int v) ? v : 1;
                    for (int i = 0; i < count; i++) dataValues.Add(val);
                }
                else
                {
                    if (int.TryParse(rawTokenStr, out int val))
                        dataValues.Add(val);
                }
                tokenIndex++;
            }
            return dataValues;
        }
    }
}