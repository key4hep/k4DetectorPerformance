#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Compare the VertexValidation output for the input of
# makeVertexValidationInput.py with the expected values.
#
#   python checkVertexValidation.py <validation file> weighted|unweighted

import math
import sys

import ROOT

# Per event: values of the scalar branches and of the per-vertex branches.
# The classification is compared by name, translated with the bin labels of
# h_vertex_classification, so the codes are defined only in VertexValidation.
# "weighted" uses the vertex -> particle links, "unweighted" counts every track
# with weight 1, which turns the truth PV vertex of event 0 from clean (4/4.5)
# into merged (4/6).
EXPECTED = {
    "weighted": [
        {
            "nRecoVtx": 3,
            "nTrueVtx": 1,
            "nTracksEvent": 8,
            "nTracksTruthPV": 4,
            "truthPVFoundClean": 1,
            "recoPVIsTruthPV": 1,
            "classification": ["clean", "clean", "unknown"],
            "matchedTruthIndex": [0, 1, -1],
            "matchFraction": [4 / 4.5, 1, 0],
            "nTracks": [5, 2, 1],
            "sumPt2": [17, 2, 1],
            "pvResX": [2],
            "pvResY": [0],
            "pvResZ": [3],
            "pvPullX": [1],
            "pvPullZ": [1.5],
            "pvTrackFraction": [1],
        },
        {
            "nRecoVtx": 3,
            "truthPVFoundClean": 1,
            # nothing flagged primary: vertex 0 is plotted as the PV but does
            # not count as selected
            "recoPVIsTruthPV": 0,
            "classification": ["clean", "split", "merged"],
            "matchedTruthIndex": [0, 0, 1],
            "matchFraction": [1, 1, 2 / 3],
            "pvResZ": [4],
            "pvPullZ": [1],
        },
        {
            "nRecoVtx": 2,
            "truthPVFoundClean": 1,
            "recoPVIsTruthPV": 1,
            "classification": ["split", "clean"],
            "matchedTruthIndex": [0, 0],
            "pvResX": [-1],
            "pvResY": [1],
            "pvPullX": [-1],
            "pvPullY": [1],
        },
    ],
    "unweighted": [
        {
            "truthPVFoundClean": 0,
            "recoPVIsTruthPV": 0,
            "classification": ["merged", "clean", "unknown"],
            "matchFraction": [4 / 6, 1, 0],
            "nTracks": [6, 2, 1],
        },
        {
            "truthPVFoundClean": 1,
            "recoPVIsTruthPV": 0,
            "classification": ["clean", "split", "merged"],
        },
        {
            "truthPVFoundClean": 1,
            "recoPVIsTruthPV": 1,
            "classification": ["split", "clean"],
        },
    ],
}

EXPECTED_SUMMARY = {
    "weighted": {
        "nEvents": 3,
        "nEventsWithRecoPV": 3,
        "nEventsTruthPVClean": 3,
        "nEventsRecoPVIsTruthPV": 2,
        "pvSelectionEfficiency": 2 / 3,
        "nClean": 4,
        "nMerged": 1,
        "nSplit": 2,
        "nUnknown": 1,
    },
    "unweighted": {
        "nEvents": 3,
        "nEventsTruthPVClean": 2,
        "nEventsRecoPVIsTruthPV": 1,
        "pvSelectionEfficiency": 1 / 3,
        "nClean": 3,
        "nMerged": 2,
        "nSplit": 2,
        "nUnknown": 1,
    },
}

# Positions are stored as float, so residuals in um are exact to ~1e-4
TOLERANCE = 1e-3


def same(value, expected):
    if isinstance(expected, str):
        return value == expected
    return math.isclose(value, expected, rel_tol=TOLERANCE, abs_tol=TOLERANCE)


def check(name, value, expected, failures):
    if isinstance(expected, list):
        value = list(value)
        ok = len(value) == len(expected) and all(same(v, e) for v, e in zip(value, expected))
    else:
        ok = same(value, expected)
    if not ok:
        failures.append(f"{name}: got {value}, expected {expected}")


def main(fileName, mode):
    f = ROOT.TFile.Open(fileName)
    if not f or f.IsZombie():
        print(f"ERROR: cannot open {fileName}")
        return 1
    failures = []

    axis = f.Get("h_vertex_classification").GetXaxis()
    classificationNames = [axis.GetBinLabel(b) for b in range(1, axis.GetNbins() + 1)]

    tree = f.Get("vertex_vs_mc")
    expected = EXPECTED[mode]
    check("entries", tree.GetEntries(), len(expected), failures)
    for i, values in enumerate(expected):
        tree.GetEntry(i)
        for branch, value in values.items():
            got = getattr(tree, branch)
            if branch == "classification":
                got = [classificationNames[int(c)] for c in got]
            check(f"event {i} {branch}", got, value, failures)

    summary = f.Get("summary")
    summary.GetEntry(0)
    for branch, value in EXPECTED_SUMMARY[mode].items():
        check(f"summary {branch}", getattr(summary, branch), value, failures)

    for failure in failures:
        print(f"FAIL {failure}")
    print(f"{mode}: {len(failures)} mismatches")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
