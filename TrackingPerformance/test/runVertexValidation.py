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

# Validate a vertex collection against the Monte Carlo truth primary vertex.
#
# Works on any EDM4hep file holding MCParticles and an edm4hep::Vertex
# collection, independent of which algorithm produced the vertices, e.g.
#
#   ACTS vertexing (k4ActsTracking VertexFindingAlg):
#     k4run runVertexValidation.py --inputFile sample_VTX.edm4hep.root \
#         --vertices ACTSPrimaryVertices \
#         --vertexParticleLinks ACTSPrimaryVertices_ParticleLinks \
#         --validationFile acts_vertices.root
#
#   LCFIPlus, from a full CLD reconstruction:
#     k4run runVertexValidation.py --inputFile sample_REC.edm4hep.root \
#         --vertices PrimaryVertices --validationFile lcfiplus_vertices.root

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, VertexValidation
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument(
    "--inputFile",
    required=True,
    help="Input EDM4hep ROOT file",
)
parser.add_argument(
    "--vertices",
    default="PrimaryVertices",
    help="Name of the edm4hep::Vertex collection to validate",
)
parser.add_argument(
    "--mcParticles",
    default="MCParticles",
    help="Name of the MC particle collection",
)
parser.add_argument(
    "--trackMCLinks",
    default="SiTracksMCTruthLink",
    help="Name of the track -> MC particle link collection used for truth matching",
)
parser.add_argument(
    "--vertexParticleLinks",
    default="",
    help="Optional vertex -> particle link collection carrying the track weights of an adaptive vertex fit "
    "(e.g. ACTSPrimaryVertices_ParticleLinks). Without it every track counts with weight 1.",
)
parser.add_argument(
    "--validationFile",
    default="vertex_validation.root",
    help="Output ROOT file with the validation tree and plots",
)
args = parser.parse_known_args()[0]

io = IOSvc("IOSvc")
io.Input = args.inputFile

validation = VertexValidation(
    "VertexValidation",
    MCParticles=[args.mcParticles],
    Vertices=[args.vertices],
    TrackMCLinks=[args.trackMCLinks],
    VertexParticleLinks=[args.vertexParticleLinks] if args.vertexParticleLinks else [],
    OutputFile=args.validationFile,
)

ApplicationMgr(
    TopAlg=[validation],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
