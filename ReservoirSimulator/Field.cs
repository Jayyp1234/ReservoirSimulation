using EclipseDeckReader;
using SepalSolver;
using System;
using System.Collections.Generic;
using System.Linq;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    internal class Field<T> where T : struct, IVectorND<T>
    {
        string name { get; set; }
        List<Group<T>> groups { get; set; } = new List<Group<T>>();
        internal Field(string _name, List<Group<T>> _groups)
        {
            name = _name;
            groups = _groups;
        }
    }
}
