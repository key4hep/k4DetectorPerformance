<!--
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-->
# Di-ParticlePerformance
Scripts to benchmark Di-particle performance

### Setup

Either source a stable key4hep stack:
`source /cvmfs/sw.hsf.org/key4hep/setup.sh`

or the nightlies:

`source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh`

``` export WORKING_DIR=$(pwd) ```

### Workflow
Create two MC particles, e.g. for ILD, written to an EDM4hep file:
`python3 createMCParticle_two_photons.py --pdg 22 --mass 0 --charge 0 --maxSep 90 --centralPosX 0 --centralPosY 1804.7 --centralPosZ 150 --Emax 5 --Emin 5 --output Di-particle_MC_Particles`
This creates two photons with a random separation of up to 90 mm in a patch centered on (0, 1804.7, 150)mm in the global coordinate system. The patch lies in the plane y=1804.7

Run detector simulation, providing this as input
`ddsim --compactFile $k4geo_DIR/ILD/compact/ILD_l5_o1_v02/ILD_l5_o1_v02.xml --numberOfEvents 100 --inputFiles Di-particle_MC_Particles.edm4hep.root --outputFile test_sim.edm4hep.root`

export this path for later:
```bash
export INPUTFILE_DIR=$(pwd)
```

get necessary dependencies for reco, e.g. for CLD
```bash
git clone https://github.com/key4hep/CLDConfig.git
cd CLDConfig/CLDConfig/
```
or ILD:
```bash
git clone https://github.com/iLCSoft/ILDConfig
cd ILDConfig/StandardConfig/production
export ILDCONFIG_DIR=$(pwd)
```

Run reconstruction
e.g. ILD:
```bash
cd $ILDCONFIG_DIR
mkdir RecoOut
k4run ILDReconstruction.py --inputFiles=$INPUTFILE_DIR/test_sim.edm4hep.root --compactFile $k4geo_DIR/ILD/compact/ILD_l5_o1_v02/ILD_l5_o1_v02.xml
```

The resulting reconstructed file can be investigated with
``` podio-dump podio-dump Reco_test_REC.edm4hep.root ```

Follow top-level build instructions for ```k4Performance```

Then run from ```build``` with:
```bash
k4run ../Di-ParticlePerformance/Di-ParticleAnalysis/options/runDiPhotonPerformance.py --inputFiles ../Di-ParticlePerformance/tmp/ILDConfig/StandardConfig/production/Reco_test_REC.edm4hep.root
```

