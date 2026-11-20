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

#include <edm4hep/Constants.h>
#include <edm4hep/utils/Helix.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace TrackingValidationHelpers {
namespace {
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();

  std::optional<edm4hep::utils::Helix> makeHelix(const edm4hep::Vector3d& position, const edm4hep::Vector3d& momentum,
                                                 double charge, double Bz, const edm4hep::Vector3d& reference) {
    try {
      return edm4hep::utils::Helix(position, momentum, charge, Bz, reference);
    } catch (const std::invalid_argument&) {
      return std::nullopt;
    }
  }
} // namespace

double wrapDeltaPhi(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b))
    return nan;
  return std::remainder(a - b, 2. * M_PI);
}

PCAInfoHelper PCAInfo_mm(double x, double y, double z, double px, double py, double pz, double charge, double refX,
                         double refY, double Bz) {
  const auto helix = makeHelix({x, y, z}, {px, py, pz}, charge, Bz, {refX, refY, 0.});
  if (!helix) {
    return {};
  }
  const double phi0 = helix->getPhi0();
  const double d0 = helix->getD0();
  const double pcaX = refX - d0 * std::sin(phi0);
  const double pcaY = refY + d0 * std::cos(phi0);
  return {pcaX, pcaY, helix->getZ0(), phi0, std::isfinite(pcaX) && std::isfinite(pcaY)};
}

HelixParams truthPerigeeFromMC(const edm4hep::MCParticle& mc, double Bz, double refX, double refY, double refZ) {
  const auto& p = mc.getMomentum();
  HelixParams hp;
  hp.pT = std::hypot(p.x, p.y);
  hp.p = std::hypot(hp.pT, p.z);
  const auto helix = makeHelix(mc.getVertex(), p, mc.getCharge(), Bz, {refX, refY, refZ});
  if (!helix) {
    hp.D0 = hp.Z0 = hp.phi = hp.omega = hp.tanLambda = nan;
    return hp;
  }
  hp.D0 = helix->getD0();
  hp.Z0 = helix->getZ0();
  hp.phi = helix->getPhi0();
  hp.tanLambda = helix->getTanLambda();
  // The validation plots use |Bz| for the fitted curvature convention.
  hp.omega = Bz < 0. ? -helix->getOmega() : helix->getOmega();
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
  return edm4hep::constants::fieldToCurvature * std::abs(Bz) / omega;
}

double momentumFromState(const edm4hep::TrackState& st, double Bz) {
  return ptFromState(st, Bz) * std::hypot(1., st.tanLambda);
}

} // namespace TrackingValidationHelpers
