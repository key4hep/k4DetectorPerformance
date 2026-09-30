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

#include "TrackingValidationHelpers.h"

#include <cmath>
#include <limits>

namespace TrackingValidationHelpers {
namespace {
  constexpr double momentumScale = 0.000299792458; // GeV / (T mm)
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
} // namespace

double wrapDeltaPhi(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b))
    return nan;
  return std::remainder(a - b, 2. * M_PI);
}

PCAInfoHelper PCAInfo_mm(double x, double y, double z, double px, double py, double pz, double refX, double refY) {
  PCAInfoHelper out;
  const double pt2 = px * px + py * py;
  if (!(pt2 > 0.) || !std::isfinite(pt2))
    return out;
  // Closest approach of the production direction to a beamline parallel to z.
  // a particle already on the beamline keeps its production z coordinate.
  const double pathOverPt = -((x - refX) * px + (y - refY) * py) / pt2;
  out.pcaX = x + pathOverPt * px;
  out.pcaY = y + pathOverPt * py;
  out.pcaZ = z + pathOverPt * pz;
  out.phi0 = std::atan2(py, px);
  out.ok = std::isfinite(out.pcaX) && std::isfinite(out.pcaY) && std::isfinite(out.pcaZ);
  return out;
}

HelixParams truthPerigeeFromMC(const edm4hep::MCParticle& mc, double Bz, double refX, double refY, double refZ) {
  const auto& p = mc.getMomentum();
  const auto& v = mc.getVertex();
  HelixParams hp;
  hp.pT = std::hypot(p.x, p.y);
  hp.p = std::hypot(hp.pT, p.z);
  const auto pca = PCAInfo_mm(v.x, v.y, v.z, p.x, p.y, p.z, refX, refY);
  if (!pca.ok) {
    hp.D0 = hp.Z0 = hp.phi = hp.omega = hp.tanLambda = nan;
    return hp;
  }
  // Project the displacement directly: no subtraction of large helix radii.
  hp.D0 = (-(v.x - refX) * p.y + (v.y - refY) * p.x) / hp.pT;
  hp.Z0 = pca.pcaZ - refZ;
  hp.phi = pca.phi0;
  hp.tanLambda = p.z / hp.pT;
  // Signed by the charge only, as the fitted state (and ptFromState) use |Bz|.
  hp.omega = momentumScale * std::abs(Bz) * mc.getCharge() / hp.pT;
  return hp;
}

std::optional<edm4hep::TrackState> getAtIPState(const edm4hep::Track& trk) {
  for (const auto& st : trk.getTrackStates()) {
    if (st.location == edm4hep::TrackState::AtIP) {
      return st;
    }
  }
  return std::nullopt;
}

double ptFromState(const edm4hep::TrackState& st, double Bz) {
  const double omega = std::abs(st.omega);
  if (!(omega > 0.) || !std::isfinite(omega) || !std::isfinite(Bz) || Bz == 0.)
    return nan;
  return momentumScale * std::abs(Bz) / omega;
}

double momentumFromState(const edm4hep::TrackState& st, double Bz) {
  return ptFromState(st, Bz) * std::hypot(1., st.tanLambda);
}

} // namespace TrackingValidationHelpers
