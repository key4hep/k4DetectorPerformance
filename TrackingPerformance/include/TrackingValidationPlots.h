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

#ifndef TRACKINGVALIDATIONPLOTS_H
#define TRACKINGVALIDATIONPLOTS_H

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TTree.h"

#include <cstddef>
#include <string>
#include <vector>

namespace TrackingValidationPlots {

/// Build logarithmic momentum bins
std::vector<double> makeLogBins(double min, double max, double step);

/// Result of the effective-sigma extraction
struct EffectiveSigmaResult {
  double center = 0.0;   // center of the narrowest interval
  double sigmaEff = 0.0; // half-width of the narrowest interval
  double median = 0.0;   // optional diagnostic info
  bool valid = false;
  std::size_t nEntries = 0;
};

/// Compute sigma_eff as half-width of the narrowest interval containing "fraction"of entries
EffectiveSigmaResult computeEffectiveSigma(std::vector<double> values, double fraction = 0.6827);

double computeEffectiveSigmaBootstrapError(const std::vector<double>& values, double fraction = 0.6827,
                                           int nBootstrap = 200, unsigned int seed = 12345u);

/**
 * @brief Build the d0 resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects residual values in
 momentum bins, extracts the effective sigma in each bin, and return the
 result as a TGraphErrors.
 */

TGraphErrors* makeD0ResolutionVsMomentum(TTree* tree, const char* graphName = "g_d0_resolution_vs_p", double pMin = 0.1,
                                         double pMax = 100.0, double logStep = 0.15);

TGraphErrors* makeZ0ResolutionVsMomentum(TTree* tree, const char* graphName = "g_z0_resolution_vs_p", double pMin = 0.1,
                                         double pMax = 100.0, double logStep = 0.15);

TGraphErrors* makePhiResolutionVsMomentum(TTree* tree, const char* graphName = "g_phi_resolution_vs_p",
                                          double pMin = 0.1, double pMax = 100.0, double logStep = 0.15);

TGraphErrors* makeOmegaResolutionVsMomentum(TTree* tree, const char* graphName = "g_omega_resolution_vs_p",
                                            double pMin = 0.1, double pMax = 100.0, double logStep = 0.15);

TGraphErrors* makeTanLambdaResolutionVsMomentum(TTree* tree, const char* graphName = "g_tanlambda_resolution_vs_p",
                                                double pMin = 0.1, double pMax = 100.0, double logStep = 0.15);

/**
 * @brief Build a pull distribution for one track parameter.
 *
 * The function reads a pull branch from the fitter validation tree and fills a
 * histogram. Pulls are normalized residuals, residual divided by the fitted
 * parameter uncertainty.
 */
TH1F* makePullHistogram(TTree* tree, const char* branchName, const char* histName, const char* title, int nBins = 100,
                        double xMin = -10.0, double xMax = 10.0);

/**
 * @brief Draw a pull histogram.
 */
TCanvas* drawPullCanvas(TH1F* h, const char* canvasName, const char* title);

/**
 * @brief Build the total-momentum resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects relative momentum
 * residual values in momentum bins, extracts the effective sigma in each bin,
 and returns the result as a TGraphErrors.
 */
TGraphErrors* makeMomentumResolutionVsMomentum(TTree* tree, const char* graphName = "g_p_resolution_vs_p",
                                               double pMin = 0.1, double pMax = 100.0, double logStep = 0.15);

/**
 * @brief Build the transverse-momentum resolution as a function of momentum.
 *
 * The function reads the fitter validation tree, collects relative transverse-
 * momentum residual values in momentum bins, extracts the effective sigma in
 each bin, and returns the result as a TGraphErrors.
 */
TGraphErrors* makePtResolutionVsMomentum(TTree* tree, const char* graphName = "g_pt_resolution_vs_p", double pMin = 0.1,
                                         double pMax = 100.0, double logStep = 0.15);

/// Draw a generic resolution graph on a logarithmic momentum axis
TCanvas* drawResolutionCanvas(TGraphErrors* g, const char* canvasName, const char* title, double xMin = 0.1,
                              double xMax = 100.0);

/**
 * @brief Build the tracking-efficiency graph as a function of momentum.
 *
 * The function reads the finder validation tree and computes the efficiency
 * according to the selected matching definition, returning the result as a
 * TGraphErrors.
 */
TGraphErrors* makeEfficiencyVsMomentum(TTree* finderTree, const char* graphName, int efficiencyDefinition,
                                       double purityThreshold, double pMin = 0.1, double pMax = 100.0,
                                       double logStep = 0.15);

/// Draw the tracking-efficiency graph on a logarithmic momentum axis
TCanvas* drawEfficiencyCanvas(TGraphErrors* g, const char* canvasName, const char* title, double xMin = 0.1,
                              double xMax = 100.0);

} // namespace TrackingValidationPlots

#endif
