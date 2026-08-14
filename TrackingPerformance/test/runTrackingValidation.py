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

import os
import math
import argparse

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc, RndmGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser


# =============================================================================
# Helper functions
# =============================================================================

def str2bool(v):
    """Parse boolean command-line arguments.

    This allows both:
      --runDigi 1
      --runDigi true
      --runDigi false
      --runDigi 0
    """
    if isinstance(v, bool):
        return v

    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True

    if v.lower() in ("no", "false", "f", "n", "0"):
        return False

    raise argparse.ArgumentTypeError("Boolean value expected")


def csv_list(value):
    """Convert a comma-separated command-line string into a Python list.

    Example:
      "A,B,C" -> ["A", "B", "C"]

    Empty strings become empty lists.
    """
    if not value:
        return []

    return [item.strip() for item in value.split(",") if item.strip()]


# =============================================================================
# Default files
# =============================================================================

# Keep the historical IDEA default, but do not crash immediately if K4GEO is not
# set. This is useful for validation-only tests where --compactFile is passed
# explicitly.
DEFAULT_IDEA_COMPACT_FILE = os.path.join(
    os.environ.get("K4GEO", ""),
    "FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml",
)


# =============================================================================
# Command-line arguments
# =============================================================================

# -----------------------------------------------------------------------------
# Input/output files
# -----------------------------------------------------------------------------

parser.add_argument(
    "--inputFile",
    required=True,
    help="Input EDM4hep ROOT file",
)

parser.add_argument(
    "--modelPath",
    default="",
    help="Path to the GGTF ONNX model. Required only if --runFinder true.",
)

parser.add_argument(
    "--outputFile",
    default="out_reco.root",
    help="Output EDM4hep ROOT file with reconstructed collections",
)

parser.add_argument(
    "--validationFile",
    default="validation.root",
    help="Output ROOT file written by TrackingValidation",
)

parser.add_argument(
    "--compactFile",
    default=DEFAULT_IDEA_COMPACT_FILE,
    help="Detector geometry XML file",
)


# -----------------------------------------------------------------------------
# Pipeline control
# -----------------------------------------------------------------------------
#
# These flags make it possible to run either:
#
#   full chain:
#       digitization -> finder -> fitter -> validation
#
#   validation-only CI:
#       read existing reco file -> validation
#
# For validation-only running, use:
#   --runDigi 0
#   --runFinder 0
#   --runFitter 0
#   --runPerfectTracking 0
#   --runValidation 1
# -----------------------------------------------------------------------------

parser.add_argument(
    "--runDigi",
    type=str2bool,
    default=True,
    help="Run digitization",
)

parser.add_argument(
    "--runFinder",
    type=str2bool,
    default=True,
    help="Run track finder",
)

parser.add_argument(
    "--runFitter",
    type=str2bool,
    default=True,
    help="Run reco fitter",
)

parser.add_argument(
    "--runPerfectTracking",
    type=str2bool,
    default=True,
    help="Run perfect tracking/perfect fitter",
)

parser.add_argument(
    "--runValidation",
    type=str2bool,
    default=True,
    help="Run TrackingValidation",
)

parser.add_argument(
    "--useDCH",
    type=str2bool,
    default=True,
    help="Use DCH collections",
)


# -----------------------------------------------------------------------------
# Validation control
# -----------------------------------------------------------------------------

parser.add_argument(
    "--mode",
    type=int,
    default=0,
    choices=[0, 1, 2],
    help="Validation mode: 0=Full, 1=FinderOnly, 2=FitterOnly",
)

parser.add_argument(
    "--doPerfectFit",
    type=str2bool,
    default=True,
    help="Fill fitter_vs_perfect validation plots",
)

parser.add_argument(
    "--finderEfficiencyDefinition",
    type=int,
    default=1,
    choices=[1, 2],
    help="1=purity-based definition, 2=purity+efficiency >= 0.5 definition",
)

parser.add_argument(
    "--finderPurityThreshold",
    type=float,
    default=0.75,
    help="Purity threshold used when FinderEfficiencyDefinition = 1",
)


# -----------------------------------------------------------------------------
# Configurable collection names
# -----------------------------------------------------------------------------
#
# The TrackingValidation algorithm itself is detector-agnostic, but different
# detector/reconstruction configurations can use different EDM4hep collection
# names. These arguments keep the IDEA defaults while allowing CLD or any other
# input file to be validated without changing the C++ algorithm.
#
# For the old IDEA full-chain use case, the defaults should work.
#
# For validation-only CLD CI, the bash test can override these with the
# collection names found in the pre-produced CLD reco file.
# -----------------------------------------------------------------------------

parser.add_argument(
    "--mcParticles",
    default="MCParticles",
    help="MC particle collection name used by TrackingValidation and perfect tracking",
)

parser.add_argument(
    "--planarHitSimLinks",
    default="SiWrBSimDigiLinks,SiWrDSimDigiLinks,VTXBSimDigiLinks,VTXDSimDigiLinks",
    help=(
        "Comma-separated planar hit-to-sim link collection names. "
        "Used by perfect tracking and, by default, by validation."
    ),
)

parser.add_argument(
    "--dchHitSimLinks",
    default="DCH_DigiSimAssociationCollection",
    help=(
        "Comma-separated DCH hit-to-sim link collection names. "
        "Used only when --useDCH true."
    ),
)

parser.add_argument(
    "--hitSimLinks",
    default="",
    help=(
        "Comma-separated hit-to-sim link collection names used by TrackingValidation. "
        "If empty, defaults to planarHitSimLinks plus dchHitSimLinks when --useDCH true."
    ),
)

parser.add_argument(
    "--planarDigiCollections",
    default="VTXBDigis,VTXDDigis,SiWrBDigis,SiWrDDigis",
    help="Comma-separated planar digi collection names used by the track finder",
)

parser.add_argument(
    "--dchDigiCollections",
    default="DCHDigis",
    help="Comma-separated DCH digi collection names used by the track finder when --useDCH true",
)

parser.add_argument(
    "--finderTracks",
    default="GGTFTracks",
    help="Finder track collection name",
)

parser.add_argument(
    "--fittedTracks",
    default="FittedTracks",
    help="Fitted track collection name",
)

parser.add_argument(
    "--perfectTracks",
    default="PerfectTracks",
    help="Perfect track collection name",
)

parser.add_argument(
    "--perfectFittedTracks",
    default="PerfectFittedTracks",
    help="Perfect fitted track collection name",
)


# =============================================================================
# Parse arguments
# =============================================================================

args = parser.parse_args()


# =============================================================================
# Validate command-line configuration
# =============================================================================

if args.runFinder and not args.modelPath:
    parser.error("--modelPath is required when --runFinder true")

if not args.compactFile:
    parser.error(
        "--compactFile is empty. Please pass --compactFile explicitly or set K4GEO."
    )

if not any(
    [
        args.runDigi,
        args.runFinder,
        args.runFitter,
        args.runPerfectTracking,
        args.runValidation,
    ]
):
    parser.error("Nothing to do: all run flags are false")

if not args.runValidation and args.doPerfectFit:
    print("WARNING: --doPerfectFit is ignored when --runValidation false")

if args.runValidation and args.mode == 1 and args.doPerfectFit:
    print("WARNING: --doPerfectFit is ignored in finder-only validation mode")

if not args.runDigi and args.runFinder:
    print(
        "WARNING: --runFinder true with --runDigi false assumes digi collections "
        "are already present in the input file"
    )

if not args.runFinder and args.runFitter:
    print(
        "WARNING: --runFitter true with --runFinder false assumes finder-track "
        "collections are already present in the input file"
    )

if (
    not args.runPerfectTracking
    and args.doPerfectFit
    and args.runValidation
    and args.mode in [0, 2]
):
    print(
        "WARNING: --doPerfectFit true with --runPerfectTracking false assumes "
        "PerfectFittedTracks is already present in the input file"
    )

if args.runValidation and args.mode in [0, 1] and not args.runFinder:
    print(
        "WARNING: Validation mode requires FinderTracks, but --runFinder false. "
        "Assuming finder tracks are already present in the input file."
    )

if args.runValidation and args.mode in [0, 2] and not args.runFitter:
    print(
        "WARNING: Validation mode requires FittedTracks, but --runFitter false. "
        "Assuming fitted tracks are already present in the input file."
    )


# =============================================================================
# Collection names
# =============================================================================
#
# These variables are used below to configure digitization, finder, fitter,
# perfect tracking, and validation.
#
# Important:
#   - PLANAR_* and DCH_* collections are mainly relevant for the full chain.
#   - HIT_SIM_LINK_COLLECTIONS is what TrackingValidation reads.
#   - If --hitSimLinks is not explicitly provided, validation uses the same
#     planar+DCH links as the full IDEA configuration.
# =============================================================================

MC_COLLECTION = args.mcParticles

PLANAR_LINK_COLLECTIONS = csv_list(args.planarHitSimLinks)
DCH_LINK_COLLECTIONS = csv_list(args.dchHitSimLinks) if args.useDCH else []

if args.hitSimLinks:
    HIT_SIM_LINK_COLLECTIONS = csv_list(args.hitSimLinks)
else:
    HIT_SIM_LINK_COLLECTIONS = PLANAR_LINK_COLLECTIONS + DCH_LINK_COLLECTIONS

PLANAR_DIGI_COLLECTIONS = csv_list(args.planarDigiCollections)
DCH_DIGI_COLLECTIONS = csv_list(args.dchDigiCollections) if args.useDCH else []

FINDER_TRACK_COLLECTION = args.finderTracks
FITTED_TRACK_COLLECTION = args.fittedTracks
PERFECT_TRACK_COLLECTION = args.perfectTracks
PERFECT_FITTED_TRACK_COLLECTION = args.perfectFittedTracks


# Print a compact configuration summary. This is useful in CI logs.
print("=== Tracking validation steering configuration ===")
print(f"Input file:                  {args.inputFile}")
print(f"Output EDM4hep file:         {args.outputFile}")
print(f"Validation ROOT file:        {args.validationFile}")
print(f"Compact file:                {args.compactFile}")
print(f"runDigi:                     {args.runDigi}")
print(f"runFinder:                   {args.runFinder}")
print(f"runFitter:                   {args.runFitter}")
print(f"runPerfectTracking:          {args.runPerfectTracking}")
print(f"runValidation:               {args.runValidation}")
print(f"useDCH:                      {args.useDCH}")
print(f"mode:                        {args.mode}")
print(f"doPerfectFit:                {args.doPerfectFit}")
print(f"MCParticles:                 {MC_COLLECTION}")
print(f"HitSimLinks:                 {HIT_SIM_LINK_COLLECTIONS}")
print(f"Planar digi collections:     {PLANAR_DIGI_COLLECTIONS}")
print(f"DCH digi collections:        {DCH_DIGI_COLLECTIONS}")
print(f"FinderTracks:                {FINDER_TRACK_COLLECTION}")
print(f"FittedTracks:                {FITTED_TRACK_COLLECTION}")
print(f"PerfectTracks:               {PERFECT_TRACK_COLLECTION}")
print(f"PerfectFittedTracks:         {PERFECT_FITTED_TRACK_COLLECTION}")


# =============================================================================
# IO
# =============================================================================

io = IOSvc("IOSvc")
io.Input = args.inputFile
io.Output = args.outputFile


# =============================================================================
# Geometry
# =============================================================================

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [args.compactFile]
geoservice.EnableGeant4Geo = False
geoservice.OutputLevel = INFO


# =============================================================================
# Algorithm sequence
# =============================================================================

TopAlg = []


# =============================================================================
# Digitizers
# =============================================================================
#
# This block is IDEA-specific. It is used for the original full-chain IDEA test.
# In validation-only CI on an existing reco file, run with:
#
#   --runDigi 0
#
# so none of this is executed.
# =============================================================================

if args.runDigi:
    from Configurables import DDPlanarDigi, DCHdigi_v02

    innerVertexResolution_x = 0.003
    innerVertexResolution_y = 0.003
    innerVertexResolution_t = 1000

    outerVertexResolution_x = 0.050 / math.sqrt(12)
    outerVertexResolution_y = 0.150 / math.sqrt(12)
    outerVertexResolution_t = 1000

    siWrapperResolution_x = 0.050 / math.sqrt(12)
    siWrapperResolution_y = 1.0 / math.sqrt(12)
    siWrapperResolution_t = 0.040

    vtxb_digitizer = DDPlanarDigi("VTXBdigitizer")
    vtxb_digitizer.SubDetectorName = "Vertex"
    vtxb_digitizer.IsStrip = False
    vtxb_digitizer.ResolutionU = [
        innerVertexResolution_x,
        innerVertexResolution_x,
        innerVertexResolution_x,
        outerVertexResolution_x,
        outerVertexResolution_x,
    ]
    vtxb_digitizer.ResolutionV = [
        innerVertexResolution_y,
        innerVertexResolution_y,
        innerVertexResolution_y,
        outerVertexResolution_y,
        outerVertexResolution_y,
    ]
    vtxb_digitizer.ResolutionT = [
        innerVertexResolution_t,
        innerVertexResolution_t,
        innerVertexResolution_t,
        outerVertexResolution_t,
        outerVertexResolution_t,
    ]
    vtxb_digitizer.SimTrackHitCollectionName = ["VertexBarrelCollection"]
    vtxb_digitizer.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
    vtxb_digitizer.TrackerHitCollectionName = ["VTXBDigis"]
    vtxb_digitizer.ForceHitsOntoSurface = True

    vtxd_digitizer = DDPlanarDigi("VTXDdigitizer")
    vtxd_digitizer.SubDetectorName = "Vertex"
    vtxd_digitizer.IsStrip = False
    vtxd_digitizer.ResolutionU = [
        outerVertexResolution_x,
        outerVertexResolution_x,
        outerVertexResolution_x,
    ]
    vtxd_digitizer.ResolutionV = [
        outerVertexResolution_y,
        outerVertexResolution_y,
        outerVertexResolution_y,
    ]
    vtxd_digitizer.ResolutionT = [
        outerVertexResolution_t,
        outerVertexResolution_t,
        outerVertexResolution_t,
    ]
    vtxd_digitizer.SimTrackHitCollectionName = ["VertexEndcapCollection"]
    vtxd_digitizer.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
    vtxd_digitizer.TrackerHitCollectionName = ["VTXDDigis"]
    vtxd_digitizer.ForceHitsOntoSurface = True

    siwrb_digitizer = DDPlanarDigi("SiWrBdigitizer")
    siwrb_digitizer.SubDetectorName = "SiWrB"
    siwrb_digitizer.IsStrip = False
    siwrb_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
    siwrb_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
    siwrb_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
    siwrb_digitizer.SimTrackHitCollectionName = ["SiWrBCollection"]
    siwrb_digitizer.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
    siwrb_digitizer.TrackerHitCollectionName = ["SiWrBDigis"]
    siwrb_digitizer.ForceHitsOntoSurface = True

    siwrd_digitizer = DDPlanarDigi("SiWrDdigitizer")
    siwrd_digitizer.SubDetectorName = "SiWrD"
    siwrd_digitizer.IsStrip = False
    siwrd_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
    siwrd_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
    siwrd_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
    siwrd_digitizer.SimTrackHitCollectionName = ["SiWrDCollection"]
    siwrd_digitizer.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
    siwrd_digitizer.TrackerHitCollectionName = ["SiWrDDigis"]
    siwrd_digitizer.ForceHitsOntoSurface = True

    TopAlg += [
        vtxb_digitizer,
        vtxd_digitizer,
        siwrb_digitizer,
        siwrd_digitizer,
    ]

    if args.useDCH:
        dch_digitizer = DCHdigi_v02(
            "DCHdigi2",
            InputSimHitCollection=["DCHCollection"],
            OutputDigihitCollection=["DCH_DigiCollection"],
            OutputLinkCollection=["DCH_DigiSimAssociationCollection"],
            DCH_name="DCH_v2",
            zResolution_mm=30.0,
            xyResolution_mm=0.1,
            Deadtime_ns=400.0,
            GasType=0,
            ReadoutWindowStartTime_ns=1.0,
            ReadoutWindowDuration_ns=450.0,
            DriftVelocity_um_per_ns=-1.0,
            SignalVelocity_mm_per_ns=200.0,
            OutputLevel=INFO,
        )

        TopAlg += [dch_digitizer]


# =============================================================================
# Track finder
# =============================================================================
#
# The finder is run only if --runFinder true.
#
# In validation-only CI, use:
#
#   --runFinder 0
#
# In that case, the input file must already contain the finder-track collection
# passed through --finderTracks.
# =============================================================================

if args.runFinder:
    from Configurables import GGTFTrackFinder

    ggtf = GGTFTrackFinder(
        "GGTFTrackFinder",
        InputPlanarHitCollections=PLANAR_DIGI_COLLECTIONS,
        InputWireHitCollections=DCH_DIGI_COLLECTIONS,
        OutputTracksGGTF=[FINDER_TRACK_COLLECTION],
        ModelPath=args.modelPath,
        Tbeta=0.6,
        Td=0.3,
        OutputLevel=INFO,
    )

    TopAlg += [ggtf]


# =============================================================================
# Reco fitter
# =============================================================================
#
# The fitter is run only if --runFitter true.
#
# In validation-only CI, use:
#
#   --runFitter 0
#
# In that case, the input file must already contain the fitted-track collection
# passed through --fittedTracks.
# =============================================================================

if args.runFitter:
    from Configurables import GenfitTrackFitter

    reco_fitter = GenfitTrackFitter("RecoTrackFitter")
    reco_fitter.InputTracks = [FINDER_TRACK_COLLECTION]

    reco_fitter.OutputFittedTracks = [FITTED_TRACK_COLLECTION]
    reco_fitter.OutputFittedTracksWithFilteredHits = ["FittedTracksWithFilteredHits"]
    reco_fitter.OutputFittedHits = ["FittedHits"]

    reco_fitter.RunSingleEvaluation = True
    reco_fitter.UseBrems = False
    reco_fitter.BetaInit = 100.0
    reco_fitter.BetaFinal = 0.1
    reco_fitter.BetaSteps = 15
    reco_fitter.InitializationType = 1

    reco_fitter.SkipTrackOrdering = False
    reco_fitter.ListOfTypesToSkip = [0]
    reco_fitter.FilterTrackHits = True
    reco_fitter.RunCalorimeterExtrapolation = False

    reco_fitter.OutputLevel = INFO

    TopAlg += [reco_fitter]


# =============================================================================
# Perfect tracking + perfect fitter
# =============================================================================
#
# This block is optional. For the lightweight CI validation-only test, keep it
# disabled:
#
#   --runPerfectTracking 0
#   --doPerfectFit 0
#
# If --doPerfectFit true while --runPerfectTracking false, the input file must
# already contain the perfect fitted track collection.
# =============================================================================

if args.runPerfectTracking:
    from Configurables import PerfectTrackFinder, GenfitTrackFitter

    perfect = PerfectTrackFinder("PerfectTrackFinder")
    perfect.InputMCParticles = [MC_COLLECTION]
    perfect.InputPlanarHitCollections = PLANAR_LINK_COLLECTIONS
    perfect.InputWireHitCollections = DCH_LINK_COLLECTIONS
    perfect.OutputPerfectTracks = [PERFECT_TRACK_COLLECTION]
    perfect.OutputLevel = INFO

    perfect_fitter = GenfitTrackFitter("PerfectTrackFitter")
    perfect_fitter.InputTracks = [PERFECT_TRACK_COLLECTION]

    perfect_fitter.OutputFittedTracks = [PERFECT_FITTED_TRACK_COLLECTION]
    perfect_fitter.OutputFittedTracksWithFilteredHits = [
        PERFECT_FITTED_TRACK_COLLECTION + "_FilteredHits"
    ]
    perfect_fitter.OutputFittedHits = ["PerfectFittedHits"]

    perfect_fitter.RunSingleEvaluation = True
    perfect_fitter.UseBrems = False
    perfect_fitter.BetaInit = 100.0
    perfect_fitter.BetaFinal = 0.1
    perfect_fitter.BetaSteps = 15
    perfect_fitter.InitializationType = 0

    perfect_fitter.SkipTrackOrdering = False
    perfect_fitter.ListOfTypesToSkip = [1]
    perfect_fitter.FilterTrackHits = True

    perfect_fitter.OutputLevel = INFO

    TopAlg += [perfect, perfect_fitter]


# =============================================================================
# Validation consumer
# =============================================================================
#
# This is the detector-agnostic part we want to test in CI.
#
# The algorithm does not need to know whether the input came from IDEA or CLD.
# It only needs the correct collection names:
#
#   --mcParticles
#   --hitSimLinks
#   --finderTracks
#   --fittedTracks
#   --perfectFittedTracks, only if --doPerfectFit true
#
# Mode-dependent inputs:
#   mode 0: finder + fitter validation
#   mode 1: finder-only validation
#   mode 2: fitter-only validation


# =============================================================================

if args.runValidation:
    from Configurables import TrackingValidation

    val = TrackingValidation("TrackingValidation")
    val.OutputFile = args.validationFile
    val.Mode = args.mode

    val.Bz = 2.0
    val.RefPointX = 0.0
    val.RefPointY = 0.0
    val.RefPointZ = 0.0

    val.DoPerfectFit = args.doPerfectFit and args.mode in [0, 2]
    val.FinderEfficiencyDefinition = args.finderEfficiencyDefinition
    val.FinderPurityThreshold = args.finderPurityThreshold

    val.MCParticles = [MC_COLLECTION]
    val.HitSimLinks = HIT_SIM_LINK_COLLECTIONS

    val.FinderTracks = [FINDER_TRACK_COLLECTION] if args.mode in [0, 1] else []
    val.FittedTracks = [FITTED_TRACK_COLLECTION] if args.mode in [0, 2] else []

    val.PerfectFittedTracks = (
        [PERFECT_FITTED_TRACK_COLLECTION]
        if args.doPerfectFit and args.mode in [0, 2]
        else []
    )

    val.OutputLevel = INFO

    TopAlg += [val]


# =============================================================================
# Application manager
# =============================================================================

ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[
        geoservice,
        EventDataSvc("EventDataSvc"),
        UniqueIDGenSvc("uidSvc"),
        RndmGenSvc(),
    ],
    OutputLevel=INFO,
)