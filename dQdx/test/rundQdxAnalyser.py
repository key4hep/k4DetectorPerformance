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
import os
import sys

# for CI test: create basic sim and reco files with a few events
det_model = "ILD_l5_o1_v02"
det_path  = os.environ["K4GEO"]+"/ILD/compact/"+det_model+"/"+det_model+".xml"
det_config = "/cvmfs/ilc.desy.de/sw/ILDConfig/v02-03-04/StandardConfig/production"
cwd = os.getcwd()

os.system("ddsim --enableGun --gun.distribution uniform --gun.energy '5*GeV' --gun.particle e- --numberOfEvents 3 --outputFile ILD_sim.edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile "+det_path)
# note that ILD reco needs to be run in the ILDConfig folder
os.system("cd "+det_config+"; k4run ILDReconstruction.py --inputFiles="+cwd+"/ILD_sim.edm4hep.root --outputFileBase="+cwd+"/ILD_reco --detectorModel="+det_model+"; cd "+cwd)

io_svc = IOSvc("IOSvc")
io_svc.Input = "ILD_reco_REC.edm4hep.root"

# or select data samples
if 0:
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

  io_svc.Input = inputlist

geo_svc = GeoSvc("GeoSvc")
geo_svc.detectors = [det_path]
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
    plotStuff = False,
    outputFile = "dQdxAnalyser_output.root",
    OutputLevel = INFO
)

myCounter = EventCounter("Counter",
    Frequency = 1
)

app_mgr = ApplicationMgr(
    TopAlg = [myCounter, mydQdxAnalyser],
    EvtSel = 'NONE',
    EvtMax = -1,
    ExtSvc = [geo_svc, EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal = True,
)
