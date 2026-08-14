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
#include "TrackingValidationPlots.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

// makeLogBins
// interpolateQuantile
// makeD0ResolutionVsMomentum
// makeMomentumResolutionVsMomentum
// makePtResolutionVsMomentum
// drawResolutionCanvas
// makeEfficiencyVsMomentum
// drawEfficiencyCanvas

namespace TrackingValidationPlots {

namespace {
  double interpolateQuantile(const std::vector<double>& x, double q) {
    if (x.empty())
      return std::numeric_limits<double>::quiet_NaN();
    if (x.size() == 1)
      return x[0];

    const double pos = q * (x.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(pos));
    const double frac = pos - static_cast<double>(i);

    if (i + 1 < x.size()) {
      return x[i] * (1.0 - frac) + x[i + 1] * frac;
    }
    return x[i];
  }
} // namespace

std::vector<double> makeLogBins(double min, double max, double step) {
  std::vector<double> bins;
  for (double x = std::log10(min); x <= std::log10(max); x += step) {
    bins.push_back(std::pow(10., x));
  }
  if (bins.empty() || bins.back() < max)
    bins.push_back(max);
  return bins;
}

EffectiveSigmaResult computeEffectiveSigma(std::vector<double> values, double fraction) {
  EffectiveSigmaResult out;
  out.nEntries = values.size();

  if (values.size() < 2)
    return out;
  if (!(fraction > 0.0 && fraction <= 1.0))
    return out;

  std::sort(values.begin(), values.end());
  out.median = interpolateQuantile(values, 0.5);

  const std::size_t n = values.size();
  std::size_t nWindow = static_cast<std::size_t>(std::ceil(fraction * static_cast<double>(n)));
  nWindow = std::max<std::size_t>(2, nWindow);
  nWindow = std::min<std::size_t>(n, nWindow);

  double bestWidth = std::numeric_limits<double>::infinity();
  std::size_t bestStart = 0;

  for (std::size_t i = 0; i + nWindow <= n; ++i) {
    const std::size_t j = i + nWindow - 1;
    const double width = values[j] - values[i];

    if (width < bestWidth) {
      bestWidth = width;
      bestStart = i;
    }
  }

  const double low = values[bestStart];
  const double high = values[bestStart + nWindow - 1];

  out.center = 0.5 * (low + high);
  out.sigmaEff = 0.5 * (high - low);
  out.valid = std::isfinite(out.sigmaEff);

  return out;
}

double computeEffectiveSigmaBootstrapError(const std::vector<double>& values, double fraction, int nBootstrap,
                                           unsigned int seed) {
  if (values.size() < 5)
    return 0.0;
  if (nBootstrap < 2)
    return 0.0;

  std::mt19937 rng(seed);
  std::uniform_int_distribution<std::size_t> pick(0, values.size() - 1);

  std::vector<double> boot;
  boot.reserve(nBootstrap);

  std::vector<double> sample(values.size());

  for (int ib = 0; ib < nBootstrap; ++ib) {
    for (std::size_t i = 0; i < values.size(); ++i) {
      sample[i] = values[pick(rng)];
    }

    const auto eff = computeEffectiveSigma(sample, fraction);
    if (eff.valid && std::isfinite(eff.sigmaEff)) {
      boot.push_back(eff.sigmaEff);
    }
  }

  if (boot.size() < 2)
    return 0.0;

  const double mean = std::accumulate(boot.begin(), boot.end(), 0.0) / static_cast<double>(boot.size());

  double var = 0.0;
  for (double x : boot) {
    const double dx = x - mean;
    var += dx * dx;
  }
  var /= static_cast<double>(boot.size() - 1);

  return std::sqrt(var);
}

/**
 * @brief Generic helper for building resolution-versus-momentum plots.
 *
 * The function reads a residual branch from the fitter validation tree,
 * groups the residuals in logarithmic momentum bins, computes the effective
 * sigma in each bin, and returns the result as a TGraphErrors.
 *
 * The same implementation is reused for the d0, z0, phi, omega, and
 * tan(lambda) resolution plots.
 */

namespace {
  TGraphErrors* makeResidualResolutionVsMomentum(TTree* tree, const char* residualBranchName, const char* graphName,
                                                 const char* yAxisTitle, double scaleFactor, double pMin, double pMax,
                                                 double logStep, unsigned int seedOffset) {
    if (!tree)
      return nullptr;

    std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
    const int nBins = bins.size() - 1;

    std::vector<std::vector<double>> residualsPerBin(nBins);

    std::vector<float>* residuals = nullptr;
    std::vector<float>* p_ref_vec = nullptr;

    tree->SetBranchAddress("p_ref", &p_ref_vec);
    tree->SetBranchAddress(residualBranchName, &residuals);

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
      tree->GetEntry(ievt);

      if (!p_ref_vec || !residuals)
        continue;
      if (p_ref_vec->size() != residuals->size())
        continue;

      for (size_t i = 0; i < p_ref_vec->size(); ++i) {
        const double p = (*p_ref_vec)[i];
        const double res = (*residuals)[i] * scaleFactor;

        if (!std::isfinite(p) || !std::isfinite(res))
          continue;
        if (p < pMin || p >= pMax)
          continue;

        int bin = -1;
        for (int b = 0; b < nBins; ++b) {
          if (p >= bins[b] && p < bins[b + 1]) {
            bin = b;
            break;
          }
        }
        if (bin < 0)
          continue;

        residualsPerBin[bin].push_back(res);
      }
    }

    TGraphErrors* g = new TGraphErrors();
    g->SetName(graphName);
    g->SetTitle((";p_{ref} [GeV];" + std::string(yAxisTitle)).c_str());

    int ip = 0;
    for (int b = 0; b < nBins; ++b) {
      if (residualsPerBin[b].size() < 20)
        continue;

      const auto eff = computeEffectiveSigma(residualsPerBin[b]);
      if (!eff.valid)
        continue;

      const double pCenter = std::sqrt(bins[b] * bins[b + 1]);
      const double sigmaErr = computeEffectiveSigmaBootstrapError(residualsPerBin[b], 0.6827, 200, seedOffset + b);

      g->SetPoint(ip, pCenter, eff.sigmaEff);
      g->SetPointError(ip, 0.0, sigmaErr);
      ++ip;
    }
    tree->ResetBranchAddresses();
    return g;
  }
} // namespace

/**
 * @brief Build the d0 resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects d0 residual values
 * in momentum bins, extracts the effective sigma in each bin, and returns the
 * result as a TGraphErrors.
 */

TGraphErrors* makeD0ResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax, double logStep) {
  return makeResidualResolutionVsMomentum(tree, "resD0", graphName, "#sigma(d_{0}) [#mum]", 1000.0, pMin, pMax, logStep,
                                          12345u);
}

/**
 * @brief Build the z0 resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects z0 residual values
 * in momentum bins, extracts the effective sigma in each bin, and returns the
 * result as a TGraphErrors.
 */

TGraphErrors* makeZ0ResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax, double logStep) {
  return makeResidualResolutionVsMomentum(tree, "resZ0", graphName, "#sigma(z_{0}) [#mum]", 1000.0, pMin, pMax, logStep,
                                          13345u);
}

/**
 * @brief Build the phi resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects phi residual values
 * in momentum bins, extracts the effective sigma in each bin, and returns the
 * result as a TGraphErrors.
 */

TGraphErrors* makePhiResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax,
                                          double logStep) {
  return makeResidualResolutionVsMomentum(tree, "resPhi", graphName, "#sigma(#phi) [rad]", 1.0, pMin, pMax, logStep,
                                          14345u);
}

/**
 * @brief Build the omega resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects omega residual values
 * in momentum bins, extracts the effective sigma in each bin, and returns the
 * result as a TGraphErrors.
 */

TGraphErrors* makeOmegaResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax,
                                            double logStep) {
  return makeResidualResolutionVsMomentum(tree, "resOmega", graphName, "#sigma(#omega) [1/mm]", 1.0, pMin, pMax,
                                          logStep, 15345u);
}

/**
 * @brief Build the tanLambda resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects tanLambda residual values
 * in momentum bins, extracts the effective sigma in each bin, and returns the
 * result as a TGraphErrors.
 */

TGraphErrors* makeTanLambdaResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax,
                                                double logStep) {
  return makeResidualResolutionVsMomentum(tree, "resTanLambda", graphName, "#sigma(tan#lambda)", 1.0, pMin, pMax,
                                          logStep, 16345u);
}

TGraphErrors* makeMomentumResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax,
                                               double logStep) {
  if (!tree)
    return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<std::vector<double>> residualsPerBin(nBins);

  std::vector<float>* p_ref_vec = nullptr;
  std::vector<float>* p_reco_vec = nullptr;

  tree->SetBranchAddress("p_ref", &p_ref_vec);
  tree->SetBranchAddress("p_reco", &p_reco_vec);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    tree->GetEntry(ievt);

    if (!p_ref_vec || !p_reco_vec)
      continue;
    if (p_ref_vec->size() != p_reco_vec->size())
      continue;

    for (size_t i = 0; i < p_ref_vec->size(); ++i) {
      const double pRef = (*p_ref_vec)[i];
      const double pReco = (*p_reco_vec)[i];

      if (!std::isfinite(pRef) || !std::isfinite(pReco))
        continue;
      if (pRef <= 0.)
        continue;
      if (pRef < pMin || pRef >= pMax)
        continue;

      const double res = (pReco - pRef) / pRef;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (pRef >= bins[b] && pRef < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0)
        continue;

      residualsPerBin[bin].push_back(res);
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p_{ref} [GeV];#sigma((p_{reco}-p_{ref})/p_{ref})");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (residualsPerBin[b].size() < 20)
      continue;

    const auto eff = computeEffectiveSigma(residualsPerBin[b]);
    if (!eff.valid)
      continue;

    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);
    const double sigmaErr = computeEffectiveSigmaBootstrapError(residualsPerBin[b], 0.6827, 200, 22345u + b);

    g->SetPoint(ip, pCenter, eff.sigmaEff);
    g->SetPointError(ip, 0.0, sigmaErr);
    ++ip;
  }
  tree->ResetBranchAddresses();
  return g;
}

TGraphErrors* makePtResolutionVsMomentum(TTree* tree, const char* graphName, double pMin, double pMax, double logStep) {
  if (!tree)
    return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<std::vector<double>> residualsPerBin(nBins);

  std::vector<float>* p_ref_vec = nullptr;
  std::vector<float>* pt_ref_vec = nullptr;
  std::vector<float>* pt_reco_vec = nullptr;

  tree->SetBranchAddress("p_ref", &p_ref_vec);
  tree->SetBranchAddress("pT_ref", &pt_ref_vec);
  tree->SetBranchAddress("pT_reco", &pt_reco_vec);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    tree->GetEntry(ievt);

    if (!p_ref_vec || !pt_ref_vec || !pt_reco_vec)
      continue;
    if (p_ref_vec->size() != pt_ref_vec->size())
      continue;
    if (pt_ref_vec->size() != pt_reco_vec->size())
      continue;

    for (size_t i = 0; i < p_ref_vec->size(); ++i) {
      const double pRef = (*p_ref_vec)[i];
      const double ptRef = (*pt_ref_vec)[i];
      const double ptReco = (*pt_reco_vec)[i];

      if (!std::isfinite(pRef) || !std::isfinite(ptRef) || !std::isfinite(ptReco))
        continue;
      if (ptRef <= 0.)
        continue;
      if (pRef < pMin || pRef >= pMax)
        continue;

      const double res = (ptReco - ptRef) / ptRef;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (pRef >= bins[b] && pRef < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0)
        continue;

      residualsPerBin[bin].push_back(res);
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p_{ref} [GeV];#sigma((pT_{reco}-pT_{ref})/pT_{ref})");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (residualsPerBin[b].size() < 20)
      continue;

    const auto eff = computeEffectiveSigma(residualsPerBin[b]);
    if (!eff.valid)
      continue;

    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);
    const double sigmaErr = computeEffectiveSigmaBootstrapError(residualsPerBin[b], 0.6827, 200, 32345u + b);

    g->SetPoint(ip, pCenter, eff.sigmaEff);
    g->SetPointError(ip, 0.0, sigmaErr);
    ++ip;
  }
  tree->ResetBranchAddresses();
  return g;
}

TCanvas* drawResolutionCanvas(TGraphErrors* g, const char* canvasName, const char* title, double xMin, double xMax) {
  if (!g)
    return nullptr;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(canvasName, title, 800, 600);
  c->SetLogx();

  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->SetTitle(title);
  g->GetXaxis()->SetLimits(xMin, xMax);
  g->Draw("AP");

  return c;
}

TGraphErrors* makeEfficiencyVsMomentum(TTree* finderTree, const char* graphName, int efficiencyDefinition,
                                       double purityThreshold, double pMin, double pMax, double logStep) {
  if (!finderTree)
    return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<int> nDen(nBins, 0);
  std::vector<int> nNum(nBins, 0);

  std::vector<float>* pVec = nullptr;
  std::vector<std::vector<float>>* purVec = nullptr;
  std::vector<std::vector<float>>* effVec = nullptr;

  finderTree->SetBranchAddress("p", &pVec);
  finderTree->SetBranchAddress("matchPurity", &purVec);
  finderTree->SetBranchAddress("matchEfficiency", &effVec);

  const Long64_t nEntries = finderTree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    finderTree->GetEntry(ievt);

    if (!pVec || !purVec || !effVec)
      continue;
    if (pVec->size() != purVec->size())
      continue;
    if (pVec->size() != effVec->size())
      continue;

    for (size_t i = 0; i < pVec->size(); ++i) {
      const double p = (*pVec)[i];
      if (!std::isfinite(p) || p < pMin || p >= pMax)
        continue;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (p >= bins[b] && p < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0)
        continue;

      nDen[bin]++;

      bool isMatched = false;

      const auto& purities = (*purVec)[i];
      const auto& efficiencies = (*effVec)[i];
      const size_t nMatches = std::min(purities.size(), efficiencies.size());

      for (size_t j = 0; j < nMatches; ++j) {
        const float purity = purities[j];
        const float efficiency = efficiencies[j];

        if (efficiencyDefinition == 2) {
          if (purity >= 0.5f && efficiency >= 0.5f) {
            isMatched = true;
            break;
          }
        } else {
          if (purity >= purityThreshold) {
            isMatched = true;
            break;
          }
        }
      }

      if (isMatched)
        nNum[bin]++;
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p [GeV];Tracking efficiency");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (nDen[b] == 0)
      continue;

    const double eff = double(nNum[b]) / double(nDen[b]);
    const double err = std::sqrt(eff * (1.0 - eff) / double(nDen[b]));
    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);

    g->SetPoint(ip, pCenter, eff);
    g->SetPointError(ip, 0.0, err);
    ++ip;
  }
  finderTree->ResetBranchAddresses();
  return g;
}

TCanvas* drawEfficiencyCanvas(TGraphErrors* g, const char* canvasName, const char* title, double xMin, double xMax) {
  if (!g)
    return nullptr;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(canvasName, title, 800, 600);
  c->SetLogx();

  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->SetTitle(title);
  g->GetYaxis()->SetRangeUser(0.0, 1.05);
  g->GetXaxis()->SetLimits(xMin, xMax);
  g->Draw("AP");

  return c;
}

TH1F* makePullHistogram(TTree* tree, const char* branchName, const char* histName, const char* title, int nBins,
                        double xMin, double xMax) {
  if (!tree)
    return nullptr;

  auto* branch = tree->GetBranch(branchName);
  if (!branch)
    return nullptr;

  std::vector<float>* pulls = nullptr;
  branch->SetAddress(&pulls);

  TH1F* h = new TH1F(histName, title, nBins, xMin, xMax);
  h->GetXaxis()->SetTitle(branchName);
  h->GetYaxis()->SetTitle("Entries");

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    branch->GetEntry(ievt);

    if (!pulls)
      continue;

    for (const auto pull : *pulls) {
      if (!std::isfinite(pull))
        continue;
      h->Fill(pull);
    }
  }

  branch->ResetAddress();
  return h;
}

TCanvas* drawPullCanvas(TH1F* h, const char* canvasName, const char* title) {
  if (!h)
    return nullptr;

  gStyle->SetOptStat(1100);
  gStyle->SetOptFit(1111);

  TCanvas* c = new TCanvas(canvasName, title, 900, 700);
  c->SetGrid();

  h->SetTitle(title);
  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(2);

  h->Draw();

  if (h->GetEntries() > 500) {
    TF1* gaus = new TF1((std::string(h->GetName()) + "_gaus").c_str(), "gaus", -3.0, 3.0);
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);

    h->Fit(gaus, "R");
  }

  c->Update();

  TPaveStats* stats = static_cast<TPaveStats*>(h->FindObject("stats"));
  if (stats) {
    stats->SetX1NDC(0.72);
    stats->SetX2NDC(0.92);
    stats->SetY1NDC(0.70);
    stats->SetY2NDC(0.92);
    stats->SetTextSize(0.022);
  }

  c->Modified();
  c->Update();

  return c;
}
} // namespace TrackingValidationPlots
