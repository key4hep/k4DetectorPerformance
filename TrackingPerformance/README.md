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
# TrackingValidation

## Overview

`TrackingValidation` is a validation algorithm for studying the performance (efficiency, purity, residuals, resolutions and pulls) of track finding and track fitting in the tracking reconstruction.
It is designed to compare reconstructed and fitted tracks with Monte Carlo truth information and, when enabled, with tracks obtained from perfect tracking, i.e. tracks fitted using the correct simHits from the particle truth information. The algorithm writes validation information to a ROOT output file containing TTrees and summary plots that can be used later for performance studies and plotting.

Typical use cases include:
- validation of track-finder performance,
- validation of fitted-track parameters against MC truth,
- comparison between standard reconstructed tracks and perfectly associated reference tracks.

---

## Inputs

`TrackingValidation` expects EDM4hep event content in which the relevant collections have already been produced by the preceding steps of the reconstruction chain.

### Input collection types

`TrackingValidation` consumes the following input collections:

- **MC particle collection**
  Type: `edm4hep::MCParticleCollection`

  Used as the truth reference for particle-level validation.

- **Hit-to-sim link collections**
  Type: `std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>`

  Collections used to associate reconstructed tracker hits with the originating simulated particles. The steering combines the available detector-specific link collections (e.g. silicon and drift-chamber links) into a single input vector passed to the validation algorithm.

- **Finder track collection**
  Type: `edm4hep::TrackCollection`

  Collection of tracks produced by the track-finding stage.

- **Fitted track collection**
  Type: `edm4hep::TrackCollection`

  Collection of tracks produced by the standard fitting stage.

- **Perfect fitted-track collections (optional)**
  Type: `std::vector<const edm4hep::TrackCollection*>`

  Optional reference collections produced from perfect truth-based associations and fitting, used when perfect-fit validation is enabled.

---

## Outputs

The algorithm writes a ROOT file specified by `OutputFile`.

The file contains validation TTrees for finder-level and fitter-level studies, together with summary performance plots produced in `finalize()`. The fitter validation trees store the reconstructed and reference track parameters, their residuals, the corresponding parameter uncertainties extracted from the fitted covariance matrix, and the resulting pull values for the five helix parameters.

### Output content by mode

The exact content filled in the output depends on the validation mode selected through `Mode`:

- **`Mode = 0` (full-pipeline mode)**
  Both finder-level and fitter-level validation are performed.
  The output includes the association trees and the fitter residual trees.

- **`Mode = 1` (finder-only mode)**
  Only the finder-level validation is performed.
  The finder and perfect-association trees are filled, while the fitter trees are booked in the file but are not filled.

- **`Mode = 2` (fitter-only mode)**
  Only the fitter-level validation is performed.
  The fitter trees are filled, while the finder and perfect-association trees are booked in the file but are not filled.

### Effect of `DoPerfectFit`

The flag `DoPerfectFit` controls the handling of the `fitter_vs_perfect` output:

- if **`DoPerfectFit = true`** and perfect fitted-track collections are provided, the fitter-to-perfect comparison is filled;
- if **`DoPerfectFit = false`**, the `fitter_vs_perfect` tree is still created but its per-event content remains empty;
- if **`DoPerfectFit = true`** but no perfect fitted-track collection is provided, the tree is still written and a warning is issued.

### Summary plots

In `finalize()`, the algorithm also writes summary plots to the same ROOT file, including:

- tracking efficiency vs momentum;
- helix-parameter resolutions (`d0`, `z0`, `phi`, `omega`, `tanLambda`) as a function of momentum;
- total-momentum resolution vs momentum;
- transverse-momentum resolution vs momentum;
- pull distributions for the five helix parameters (`d0`, `z0`, `phi`, `omega`, `tanLambda`).

The pull distributions are computed as
pull = (reconstructed − reference) / σ
where σ is taken from the corresponding diagonal element of the fitted track covariance matrix. Gaussian fits are performed for sufficiently populated pull distributions to facilitate validation of the covariance estimates.

Additional plots may be added in future developments.
---

## Finder validation: efficiency and purity

To evaluate finder performance, each reconstructed track is matched to the truth particle with which it shares the largest number of hits.

For each particle-track pair, the algorithm stores two standard hit-based quantities:

- **track hit purity**: the fraction of hits on the reconstructed track that originate from the matched truth particle;
- **track hit efficiency**: the fraction of the truth-particle hits that are recovered in the reconstructed track.

The summary **tracking efficiency** can then be defined in more than one way.

- **`FinderEfficiencyDefinition = 1`**
  A truth particle is counted as reconstructed if it is associated to at least one finder track with
  `purity >= FinderPurityThreshold`.
  In the default configuration, `FinderPurityThreshold = 0.75`, following the CMS association convention in which a reconstructed track is associated to a simulated particle if more than 75% of its hits originate from that particle. The tracking efficiency is then defined as the fraction of simulated tracks associated to at least one reconstructed track. :contentReference[oaicite:0]{index=0}

- **`FinderEfficiencyDefinition = 2`**
  A truth particle is counted as reconstructed if it is associated to at least one finder track with
  `purity >= 0.5` **and** `efficiency >= 0.5`.
  This corresponds to the stricter two-ratio variant, where both the purity of the reconstructed track and the fraction of recovered truth hits must exceed 50%.

In the current implementation, the denominator of the efficiency plot includes generator-level particles with status 1 and at least one truth-linked hit.
For more details on the CMS association convention and the related definitions of tracking efficiency, fake rate, and duplicate rate, see the CMS performance note [*Performance of the track selection DNN in Run 3*](https://cds.cern.ch/record/2854696/files/DP2023_009.pdf).




---

## How to run

`TrackingValidation` is tested through a lightweight `ctest` workflow. The CI test does **not** run the full reconstruction chain. Instead, it starts from a small pre-produced CLD reconstruction file and runs only the validation algorithm.

In the current setup:

- the input CLD reconstruction file is retrieved through CMake `ExternalData`;
- the validation steering is controlled by `runTrackingValidation.py`;
- the shell test `testTrackingValidation.sh` runs only `TrackingValidation`;

### Running the CI test

From the build directory run the validation test:

```bash
ctest -V -R testTrackingValidation
```

The test writes the validation output file to the test working directory, which is the build directory configured by CMake. The output file is named:

```text
validation.root
```

### CI test workflow

The registered CI test runs the following reduced workflow:

```text
pre-produced CLD reco file → TrackingValidation → validation.root
```

It intentionally does **not** run:

- DDSim;
- digitization;
- track finding;
- track fitting;
- perfect tracking.

The test uses the CLD geometry file from `k4geo`:

```text
${K4GEO}/FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml
```

and the pre-produced input file:

```text
MuGuns_CLD_o3_v01_2026_07_01.root
```

which is provided through CMake `ExternalData`.

### CI test configuration

The validation-only CI test uses the following steering configuration:

- `runDigi = false`
  Do not run digitization.

- `runFinder = false`
  Do not run the track finder. Finder tracks are read from the input file.

- `runFitter = false`
  Do not run the track fitter. Fitted tracks are read from the input file.

- `runPerfectTracking = false`
  Do not run perfect tracking or perfect fitting.

- `runValidation = true`
  Run the validation algorithm.

- `useDCH = false`
  Do not use drift-chamber collections. The CI input is a CLD reconstruction file.

- `mode = 0`
  Run full validation, using both finder-level and fitter-level inputs from the pre-produced file.

- `doPerfectFit = false`
  Do not fill the fitter-versus-perfect-track comparison.

- `finderEfficiencyDefinition = 1`
  Use the purity-based tracking-efficiency definition.

- `finderPurityThreshold = 0.75`
  Use a purity threshold of 0.75 when `FinderEfficiencyDefinition = 1`.

The CLD collection names are passed explicitly to the steering file:

```text
mcParticles  = MCPhysicsParticles
hitSimLinks  = VXDTrackerHitRelations
finderTracks = SiTracks
fittedTracks = FittedTracks
```

The boolean steering options accept both `true/false` and `1/0` inputs.


