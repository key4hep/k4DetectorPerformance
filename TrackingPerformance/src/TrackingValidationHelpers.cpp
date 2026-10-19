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

PCAInfoHelper PCAInfo_mm(double x, double y, double z, double px, double py, double pz, double charge, double refX,
                         double refY, double Bz) {
  PCAInfoHelper out;
  const double pt = std::hypot(px, py);
  if (!(pt > 0.) || !std::isfinite(pt) || !std::isfinite(charge) || !std::isfinite(Bz))
    return out;
  const double dx = x - refX;
  const double dy = y - refY;
  const double sinPhi = py / pt;
  const double cosPhi = px / pt;
  // Signed curvature, > 0 for a particle turning clockwise (q Bz > 0)
  const double rho = momentumScale * charge * Bz / pt;
  double pathT; // transverse path length from the production point to the PCA
  if (rho == 0.) {
    // Neutral particle or no field: straight line
    pathT = -(dx * cosPhi + dy * sinPhi);
    out.pcaX = x + pathT * cosPhi;
    out.pcaY = y + pathT * sinPhi;
    out.phi0 = std::atan2(py, px);
  } else {
    // The PCA lies on the line from the reference point to the circle centre (xc, yc),
    // at distance dist - radius from the reference point.
    const double xc = dx + sinPhi / rho;
    const double yc = dy - cosPhi / rho;
    const double dist = std::hypot(xc, yc);
    // dist - radius without subtracting two numbers of the size of the radius
    const double distMinusRadius =
        (dx * dx + dy * dy + 2. * (dx * sinPhi - dy * cosPhi) / rho) / (dist + 1. / std::abs(rho));
    out.pcaX = refX + xc * distMinusRadius / dist;
    out.pcaY = refY + yc * distMinusRadius / dist;
    const double turn = rho > 0. ? 1. : -1.;
    out.phi0 = std::atan2(turn * xc, -turn * yc);
    pathT = std::remainder(std::atan2(py, px) - out.phi0, 2. * M_PI) / rho;
  }
  out.pcaZ = z + pathT * pz / pt;
  out.ok = std::isfinite(out.pcaX) && std::isfinite(out.pcaY) && std::isfinite(out.pcaZ) && std::isfinite(out.phi0);
  return out;
}

HelixParams truthPerigeeFromMC(const edm4hep::MCParticle& mc, double Bz, double refX, double refY, double refZ) {
  const auto& p = mc.getMomentum();
  const auto& v = mc.getVertex();
  HelixParams hp;
  hp.pT = std::hypot(p.x, p.y);
  hp.p = std::hypot(hp.pT, p.z);
  const auto pca = PCAInfo_mm(v.x, v.y, v.z, p.x, p.y, p.z, mc.getCharge(), refX, refY, Bz);
  if (!pca.ok) {
    hp.D0 = hp.Z0 = hp.phi = hp.omega = hp.tanLambda = nan;
    return hp;
  }
  hp.D0 = -(pca.pcaX - refX) * std::sin(pca.phi0) + (pca.pcaY - refY) * std::cos(pca.phi0);
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
