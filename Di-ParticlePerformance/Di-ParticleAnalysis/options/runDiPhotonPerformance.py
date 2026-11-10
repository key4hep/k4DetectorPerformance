#!/usr/bin/env python3
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

from Gaudi.Configuration import INFO
from Gaudi.Configuration import DEBUG
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

from Configurables import (
    DiPhotonAnalysis,
    EventDataSvc,
    RootHistSvc,
    AuditorSvc,
    AlgTimingAuditor,
)

from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink

parser_group = parser.add_argument_group("runDiPhotonPerformance.py custom options")
parser_group.add_argument("-in","--inputFiles", action="extend", nargs="+", metavar=("file1", "file2"), help="One or multiple input files")
parsed_args = parser.parse_known_args()[0]

iosvc = IOSvc()
iosvc.Input = parsed_args.inputFiles

iosvc.CollectionNames = ["PandoraPFOs", "MCParticles", "PandoraClusters", "RecoMCTruthLink"]
ILDDiPhoton = DiPhotonAnalysis("ILDDiPhotons")
ILDDiPhoton.OutputLevel = DEBUG
ILDDiPhoton.InputPFOs = ["PandoraPFOs"]
ILDDiPhoton.InputMCParticles = ["MCParticles"]
#ILDDiPhoton.InputClusters = ["PandoraClusters"]
ILDDiPhoton.InputRecoMC   = ["RecoMCTruthLink"]
#ILDDiPhoton.HistPath = "/PLOTS"

# Use Gaudi Auditor service to get algorithm timing information
auditorSvc = AuditorSvc()
auditorSvc.Auditors = [AlgTimingAuditor()]

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "DiPhoton_histograms.root"

# Configure application manager
app_mgr = ApplicationMgr(
        TopAlg=[ILDDiPhoton],
        EvtSel="NONE",
        EvtMax= -1,
        ExtSvc=[EventDataSvc("EventDataSvc"), auditorSvc, iosvc, root_hist_svc], #THistSvc("THistSvc"),iosvc],
        OutputLevel=DEBUG,
)