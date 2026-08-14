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

// Constants matching the fitter code
static constexpr float c_mm_s = 2.998e11f;
static constexpr float a_genfit = 1e-15f * c_mm_s;

float wrapDeltaPhi(float a, float b) {
  float d = a - b;
  while (d > M_PI)
    d -= 2.f * M_PI;
  while (d < -M_PI)
    d += 2.f * M_PI;
  return d;
}

PCAInfoHelper PCAInfo_mm(float x, float y, float z, float px, float py, float pz, int chargeSign, float refX,
                         float refY, float Bz) {
  PCAInfoHelper out;

  const float pt = std::sqrt(px * px + py * py);
  if (pt == 0.f)
    return out;
  if (chargeSign == 0)
    chargeSign = 1;
  if (Bz == 0.f)
    return out;

  const float R = pt / (0.3f * std::abs(chargeSign) * Bz) * 1000.f;

  const float tx = px / pt;
  const float ty = py / pt;

  const float nx = float(chargeSign) * ty;
  const float ny = float(chargeSign) * (-tx);

  const float xc = x + R * nx;
  const float yc = y + R * ny;

  const float vx = refX - xc;
  const float vy = refY - yc;
  const float vxy = std::sqrt(vx * vx + vy * vy);
  if (vxy == 0.f)
    return out;

  const float ux = vx / vxy;
  const float uy = vy / vxy;

  const float pcaX = xc + R * ux;
  const float pcaY = yc + R * uy;

  const float rx = pcaX - xc;
  const float ry = pcaY - yc;

  const int sign = (chargeSign > 0) ? 1 : -1;
  float tanX = -sign * ry;
  float tanY = sign * rx;

  const float tnorm = std::sqrt(tanX * tanX + tanY * tanY);
  if (tnorm == 0.f)
    return out;

  tanX /= tnorm;
  tanY /= tnorm;

  const float phi0 = std::atan2(tanY, tanX);

  const float pR = pt;
  const float pZ = pz;
  const float R0 = std::sqrt(x * x + y * y);
  const float Z0 = z;

  const float denom = (pR * pR + pZ * pZ);
  if (denom == 0.f)
    return out;

  const float tPCA = -(R0 * pR + Z0 * pZ) / denom;
  const float ZPCA = Z0 + pZ * tPCA;

  out.pcaX = pcaX;
  out.pcaY = pcaY;
  out.pcaZ = ZPCA;
  out.phi0 = phi0;
  out.ok = true;
  return out;
}

HelixParams truthFromMC_GenfitConvention(const edm4hep::MCParticle& mc, float Bz, float refX, float refY, float refZ) {
  HelixParams hp;

  const auto& mom = mc.getMomentum();
  const float px = float(mom.x);
  const float py = float(mom.y);
  const float pz = float(mom.z);

  const float pT = std::sqrt(px * px + py * py);
  const float p = std::sqrt(px * px + py * py + pz * pz);

  hp.pT = pT;
  hp.p = p;

  int qSign = 1;
  if (mc.getCharge() < 0.f)
    qSign = -1;

  const auto& v = mc.getVertex();
  const float x = float(v.x);
  const float y = float(v.y);
  const float z = float(v.z);

  const auto info = PCAInfo_mm(x, y, z, px, py, pz, qSign, refX, refY, Bz);
  if (!info.ok) {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    hp.D0 = NaN;
    hp.Z0 = NaN;
    hp.phi = NaN;
    hp.omega = NaN;
    hp.tanLambda = NaN;
    return hp;
  }

  hp.D0 = ((-(refX - info.pcaX)) * std::sin(info.phi0) + (refY - info.pcaY) * std::cos(info.phi0));
  hp.Z0 = (info.pcaZ - refZ);
  hp.phi = std::atan2(py, px);
  hp.tanLambda = (pT > 0.f) ? (pz / pT) : 0.f;
  hp.omega = (pT > 0.f) ? (std::abs(a_genfit * Bz / pT) * float(qSign)) : 0.f;

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

float ptFromState(const edm4hep::TrackState& st, float Bz) {
  const float omega = std::abs(float(st.omega));
  if (omega == 0.f)
    return 0.f;
  return a_genfit * std::abs(Bz) / omega;
}

float momentumFromState(const edm4hep::TrackState& st, float Bz) {
  const float pT = ptFromState(st, Bz);
  const float tl = float(st.tanLambda);
  return pT * std::sqrt(1.f + tl * tl);
}

} // namespace TrackingValidationHelpers
