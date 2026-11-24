namespace EclipseDeckReader;

/// <summary>
/// A named fault from FAULTS, with its transmissibility multiplier from MULTFLT.
/// </summary>
public sealed class FaultRecord
{
    // Every FAULTS record with this name. A fault is normally written as many records (one per stretch of
    // cell faces), so all of them are kept, not just the last.
    public List<FaultSegment> Segments { get; } = [];

    // Transmissibility multiplier applied to connection faces (from MULTFLT)
    public double TransmissibilityMultiplier { get; set; } = 1.0;
}

/// <summary>
/// One FAULTS record: cells I1-I2, J1-J2, K1-K2 (1-based) and the face, kept in full: "X" is the face towards
/// I+1, "X-" the face towards I-1 (likewise Y/J and Z/K).
/// </summary>
public readonly record struct FaultSegment(int IStart, int IEnd, int JStart, int JEnd, int KStart, int KEnd, string Face);
