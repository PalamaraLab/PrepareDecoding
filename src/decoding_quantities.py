# This file is part of https://github.com/PalamaraLab/PrepareDecoding
# which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Compare decoding quantities and print differences

import sys
import argparse
from pathlib import Path
import numpy as np
from enum import Enum, auto
import gzip

from typing import Tuple, Optional, List


class ParserState(Enum):
    TransitionType = auto()
    States = auto()
    CSFSSamples = auto()
    TimeVector = auto()
    SizeVector = auto()
    Discretization = auto()
    ExpectedTimes = auto()
    CSFS = auto()
    FoldedCSFS = auto()
    ClassicEmission = auto()
    AscertainedCSFS = auto()
    FoldedAscertainedCSFS = auto()
    CompressedAscertainedEmission = auto()
    initialStateProb = auto()
    ColumnRatios = auto()
    RowRatios = auto()
    Uvectors = auto()
    Bvectors = auto()
    Dvectors = auto()
    HomozygousEmissions = auto()

    @staticmethod
    def line(s: str):
        tokens = s.split()
        if tokens:
            return (
                ParserState[tokens[0]] if tokens[0] in ParserState.__members__ else None
            )
        else:
            return None


class Type(Enum):
    TransitionType = auto()
    Integer = auto()
    Vector = auto()
    Matrix = auto()
    TupleFloatVector = auto()

    @staticmethod
    def for_state(state: ParserState):  # NOQA
        M = """
TransitionType: TransitionType
States: Integer
CSFSSamples: Integer
TimeVector: Vector
SizeVector: Vector
Discretization: Vector
ExpectedTimes: Vector
CSFS: Matrix
FoldedCSFS: Matrix
ClassicEmission: Matrix
AscertainedCSFS: Matrix
FoldedAscertainedCSFS: Matrix
CompressedAscertainedEmission: Matrix
initialStateProb: Vector
ColumnRatios: Vector
RowRatios: TupleFloatVector
Uvectors: TupleFloatVector
Bvectors: TupleFloatVector
Dvectors: TupleFloatVector
HomozygousEmissions: TupleFloatVector
"""
        D = dict(map(str.strip, x.split(":")) for x in M.strip().split("\n"))
        if state.name in D:
            return Type[D[state.name]]
        else:
            raise ValueError("Type of {} not specified".format(state))

    def parse(self, s: str):
        s = s.strip()
        if self is Type.TransitionType:
            return s if s in ["SMC", "SMC1", "CSC"] else None
        elif self is Type.Integer:
            return int(s)
        elif self is Type.Vector:
            return np.array([float(x) for x in s.split()])
        elif self is Type.Matrix:
            return [float(x) for x in s.split()]
        elif self is Type.TupleFloatVector:
            nums = s.split()
            return float(nums[0]), np.array([float(x) for x in nums[1:]])


class DecodingQuantities:
    "DecodingQuantities representation"

    def __init__(self, fn: Path, parse_now: bool = True):
        self.TransitionType = None
        self.data = {
            "CSFS": [],
            "FoldedCSFS": [],
            "AscertainedCSFS": [],
            "FoldedAscertainedCSFS": [],
            "CompressedAscertainedEmission": [],
            "ClassicEmission": [],
            "initialStateProb": None,
            "ColumnRatios": None,
            "RowRatios": [],
            "Uvectors": [],
            "Bvectors": [],
            "Dvectors": [],
            "HomozygousEmissions": [],
        }
        self.fn = fn
        self.fp = gzip.open(fn)
        if parse_now:
            self.parse()

    def __repr__(self):
        return repr(self.data)

    def parse(self):
        state = None
        for line in self.fp.readlines():
            line = line.decode("utf-8").strip()
            if line == "":
                continue
            current_state = ParserState.line(line)
            if current_state:
                state = current_state
                if Type.for_state(state) is Type.Matrix:
                    self.data[state.name].append([])
                continue
            else:
                state_type = Type.for_state(state)
                if state_type not in [Type.Matrix, Type.TupleFloatVector]:
                    self.data[state.name] = state_type.parse(line)
                elif state_type is Type.TupleFloatVector:
                    self.data[state.name].append(state_type.parse(line))
                else:
                    self.data[state.name][-1].append(state_type.parse(line))

        # convert array of arrays to numpy array
        for m in [
            "CSFS",
            "FoldedCSFS",
            "AscertainedCSFS",
            "FoldedAscertainedCSFS",
            "ClassicEmission",
            "CompressedAscertainedEmission",
        ]:
            self.data[m] = np.array(self.data[m])

    @staticmethod
    def compare_tuplefloatvectors(a1, a2):
        return a1[0] == a2[0] and np.allclose(a1[1], a2[1])

    def equals(self, other, quiet=False) -> Tuple[bool, Optional[List[str]]]:
        "Returns (True, None) if the DecodingQuantities are equal, (False, error_list) otherwise"
        o = {}
        for identical_attribute in ["TransitionType", "States", "CSFSSamples"]:
            o[identical_attribute] = (
                self.data[identical_attribute] == other.data[identical_attribute]
            )
        for csfs in [
            "CSFS",
            "FoldedCSFS",
            "AscertainedCSFS",
            "FoldedAscertainedCSFS",
        ]:
            if len(self.data[csfs]) != len(other.data[csfs]):
                o[csfs] = False
                continue
            for i in range(len(self.data[csfs])):
                o[csfs + "." + str(i)] = (
                    self.data[csfs][i].shape == other.data[csfs][i].shape
                ) and np.allclose(self.data[csfs][i], other.data[csfs][i])

        for tuple_attribute in [
            "RowRatios",
            "Uvectors",
            "Bvectors",
            "Dvectors",
            "HomozygousEmissions",
        ]:
            if len(self.data[tuple_attribute]) != len(other.data[tuple_attribute]):
                o[tuple_attribute] = False
                continue
            for i in range(len(self.data[tuple_attribute])):
                o[tuple_attribute + "." + str(i)] = self.compare_tuplefloatvectors(
                    self.data[tuple_attribute][i], other.data[tuple_attribute][i]
                )

        for matrix in [
            "ClassicEmission",
            "CompressedAscertainedEmission",
            "ColumnRatios",
            "TimeVector",
            "SizeVector",
            "Discretization",
            "ExpectedTimes",
            "initialStateProb",
        ]:
            o[matrix] = np.allclose(
                self.data[matrix],
                other.data[matrix],
                rtol=1e-10,
                atol=1e-14,
            )
        result = all(o.values())
        if not result:
            if quiet:
                return False, set(k.split(".")[0] for k in sorted(o) if not o[k])
            else:
                return False, set(k for k in sorted(o) if not o[k])
        return True, None

    def __eq__(self, other):
        is_equal, _ = self.equals(other, quiet=True)
        return is_equal


def main():
    parser = argparse.ArgumentParser(
        description="Compare decoding quantities and show differences"
    )
    parser.add_argument("file1", help="First file to compare")
    parser.add_argument("file2", help="Second file to compare")
    parser.add_argument(
        "-q",
        "--quiet",
        help="Quiet mode, show attributes without lines",
        action="store_true",
    )
    args = parser.parse_args()
    d1 = DecodingQuantities(args.file1)
    d2 = DecodingQuantities(args.file2)
    is_equal, errors = d1.equals(d2)
    if not is_equal:
        print("\n".join(errors))
    sys.exit(int(not is_equal))  # 0 (success) if is_equal is True


if __name__ == "__main__":
    main()
