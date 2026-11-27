namespace ReservoirSimulator
{
   public enum CellPhaseState
    {
        TwoPhaseOilGas,
        UnderSaturatedOil,
        UnderSaturatedGas,
    }

    public enum WellStatus
    {
        OPEN,
        STOP,
        SHUT
    }

    public enum ConstraintType
    {
        MAX_BHP,
        MIN_BHP,
        LRATE,
        ORATE,
        WRATE,
        GRATE,
        RATE,
        RESV
    }

    public enum WellType
    {
        INJECTOR,
        PRODUCER
    }

    public enum Units
    {
        FIELD,
        METRIC,
        LAB,
        PVT_M
    }

    public enum Phase
    {
        OIL,
        WATER,
        GAS,
        DISGAS,
        VAPOIL
    }

    public enum Direction
    {
        X,
        Y,
        Z
    }

    public enum FlowDirection
    {
        KM,
        JM,
        IM,
        IP,
        JP,
        KP
    }
}
