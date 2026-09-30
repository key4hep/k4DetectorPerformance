/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TRACKINGVALIDATIONHELPERS_H
#define TRACKINGVALIDATIONHELPERS_H

#include "edm4hep/MCParticle.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackState.h"
#include <optional>

namespace TrackingValidationHelpers {

/// Purity and completeness both required by FinderEfficiencyDefinition = 2
/// (ACTS TrackTruthMatcher with doubleMatching and matchingRatio = 0.5)
inline constexpr double doubleMatchingRatio = 0.5;

/// EDM4hep helix parameters; geometry is evaluated in double precision
struct HelixParams {
  double D0 = 0.;
  double Z0 = 0.;
  double phi = 0.;
  double omega = 0.;
  double tanLambda = 0.;
  double p = 0.;
  double pT = 0.;
};

/// Helper container for PCA position and tangent angle
struct PCAInfoHelper {
  double pcaX = 0.;
  double pcaY = 0.;
  double pcaZ = 0.;
  double phi0 = 0.;
  bool ok = false;
};
/// Wrap a phi difference into the interval [-pi, pi]
double wrapDeltaPhi(double a, double b);

/// Straight-line closest approach of the production direction to the beamline
/// through (refX, refY), in mm
PCAInfoHelper PCAInfo_mm(double x, double y, double z, double px, double py, double pz, double refX, double refY);

/**
 * @brief Truth helix parameters of an MC particle on the perigee at (refX, refY, refZ).
 *
 * As in ACTS (ResPlotTool::fill, RootTrackSummaryWriter; acts-project/acts @ ce4824f9),
 * d0 and z0 come from the straight line along the production direction, and phi and
 * tanLambda from the production momentum. Exact for particles produced on the beamline;
 * for displaced particles the neglected curvature biases d0 by ~L^2/2R and phi by ~L/R.
 */
HelixParams truthPerigeeFromMC(const edm4hep::MCParticle& mc, double Bz, double refX, double refY, double refZ);

/// Retrieve the track state stored at the interaction point, if available
std::optional<edm4hep::TrackState> getAtIPState(const edm4hep::Track& trk);

/// Compute the transverse momentum from a track state
double ptFromState(const edm4hep::TrackState& st, double Bz);

/// Compute the total momentum from a track state
double momentumFromState(const edm4hep::TrackState& st, double Bz);

} // namespace TrackingValidationHelpers

#endif
