using System;
using System.Collections.Generic;

namespace EclipseDeckReader
{
    public sealed class WellSpecRecord
    {
        public string WellName { get; set; } = string.Empty;
        public string GroupName { get; set; } = string.Empty;
        public int HeadI { get; set; }
        public int HeadJ { get; set; }
        public double ReferenceDepth { get; set; }
        public string PreferredPhase { get; set; } = string.Empty;
    }

    public sealed class WellCompletionRecord
    {
        public DateTime EffectiveDate { get; set; }
        public int I { get; set; }
        public int J { get; set; }
        public int KStart { get; set; }
        public int KEnd { get; set; }
        public string Status { get; set; } = "OPEN"; // Connection-level switch (OPEN / SHUT)
        public double WellBoreDiameter { get; set; }
        public double EffectiveSkin { get; set; }
        public double TransmissibilityFactor { get; set; } // Optional: Col 8 of COMPDAT (WPIMULT / user KH)
    }

    public sealed class ProducerControlRecord
    {
        public DateTime EffectiveDate { get; set; }
        public string Status { get; set; } = "OPEN"; // OPEN, STOP, or SHUT
        public string ControlMode { get; set; } = string.Empty;  // ORAT, WRAT, GRAT, LRAT, RESV, BHP
        public double OilRateTarget { get; set; }
        public double WaterRateTarget { get; set; }
        public double GasRateTarget { get; set; }
        public double LiquidRateTarget { get; set; }
        public double ResvRateTarget { get; set; }
        public double MinBhpLimit { get; set; }
    }

    public sealed class InjectorControlRecord
    {
        public DateTime EffectiveDate { get; set; }
        public string InjectionPhase { get; set; } = "WATER";    // WATER, GAS, OIL
        public string Status { get; set; } = "OPEN"; // OPEN, STOP, or SHUT
        public string ControlMode { get; set; } = string.Empty;  // RATE, RESV, BHP
        public double TargetRate { get; set; }
        public double ResvRateTarget { get; set; }
        public double MaxBhpLimit { get; set; }
    }

    public sealed class WellData
    {
        public WellSpecRecord Specification { get; set; } = new();

        // Single mechanical completion history for the wellbore (from COMPDAT)
        public List<WellCompletionRecord> CompletionHistory { get; } = [];

        // Separate operational control histories (WCONPROD vs WCONINJE)
        public List<ProducerControlRecord> ProducerControlHistory { get; } = [];
        public List<InjectorControlRecord> InjectorControlHistory { get; } = [];
    }
}