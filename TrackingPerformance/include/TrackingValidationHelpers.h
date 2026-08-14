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

/// Container for helix parameters in the fitter convention
struct HelixParams {
  float D0 = 0.f;
  float Z0 = 0.f;
  float phi = 0.f;
  float omega = 0.f;
  float tanLambda = 0.f;
  float p = 0.f;
  float pT = 0.f;
};

/// Helper container for PCA position and tangent angle
struct PCAInfoHelper {
  float pcaX = 0.f;
  float pcaY = 0.f;
  float pcaZ = 0.f;
  float phi0 = 0.f;
  bool ok = false;
};
/// Wrap a phi difference into the interval [-pi, pi]
float wrapDeltaPhi(float a, float b);

/**
 * @brief Compute the PCA position and tangent angle in mm using a GenFit-like convention.
 *
 * The function derives the point of closest approach to the reference point
 * and the corresponding tangent direction from the input position, momentum,
 * charge sign, and magnetic field.
 */

PCAInfoHelper PCAInfo_mm(float x, float y, float z, float px, float py, float pz, int chargeSign, float refX,
                         float refY, float Bz);

/**
 * @brief Build truth helix parameters in the same convention used by the fitter.
 *
 * The returned parameters are derived from the MC particle kinematics and vertex
 * using the same reference-point and helix convention adopted for fitted tracks,
 * so that residuals can be computed consistently.
 */

HelixParams truthFromMC_GenfitConvention(const edm4hep::MCParticle& mc, float Bz, float refX, float refY, float refZ);

/// Retrieve the track state stored at the interaction point, if available
std::optional<edm4hep::TrackState> getAtIPState(const edm4hep::Track& trk);

/// Compute the transverse momentum from a track state
float ptFromState(const edm4hep::TrackState& st, float Bz);

/// Compute the total momentum from a track state
float momentumFromState(const edm4hep::TrackState& st, float Bz);

} // namespace TrackingValidationHelpers

#endif
