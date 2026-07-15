#
# Copyright (c) 2020-2026 Key4hep-Project.
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
from k4FWCore import IOSvc, ApplicationMgr
from Configurables import EventDataSvc, GeoSvc, AuditorSvc, ChronoAuditor, EventCounter

# select data samples
infolder = "/eos/user/u/ueinhaus/public/dQdx_test/singleParticles/" #  or download from "https://cernbox.cern.ch/s/orgs99sJNDmDE2t"
filelist = "singleParticles.txt"
#infolder = "/eos/user/u/ueinhaus/public/dQdx_test/2f-Z-had-250/" # or download from "https://cernbox.cern.ch/s/hCfz1O1qBwihKMb"
#filelist = "2f-Z-had-250.txt"

inputlist = []
filelist = open(filelist, 'r')
for line in filelist:
  filename = line.split()[0]
  if filename[0]=="#":
    continue
  infile = infolder + filename
  inputlist.append(infile)

io_svc = IOSvc("IOSvc")
io_svc.Input = inputlist

geo_svc = GeoSvc("GeoSvc")
lcgeo_DIR = "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/k4geo/0.21-d5gbnd/share/k4geo"
geo_svc.detectors = [lcgeo_DIR+"/ILD/compact/ILD_l5_o1_v02/ILD_l5_o1_v02.xml"]
geo_svc.EnableGeant4Geo = False
geo_svc.OutputLevel = INFO

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

from Configurables import dQdxAnalyser

mydQdxAnalyser = dQdxAnalyser("dQdxAnalyser",
    InputdQdxCollection =  ["MarlinTrkTracks_dQdx"],
    InputTrackCollection = ["MarlinTrkTracks"],
    InputRecoParticleCollection = ["PandoraPFOs"],
    usePFOTracks = True,
    plotStuff = True,
    plotFolder = "Plots",
    OutputLevel = INFO
)

myCounter = EventCounter("Counter",
    Frequency = 10000
)

app_mgr = ApplicationMgr(
    TopAlg = [myCounter, mydEdxAnalyser],
    EvtSel = 'NONE',
    EvtMax = -1,
    ExtSvc = [geo_svc, EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal = True,
)
