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

// k4FWCore
#include "k4FWCore/Consumer.h"

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiKernel/MsgStream.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/VertexRecoParticleLinkCollection.h"

// ROOT
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"

// STL
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

/** @struct VertexValidation
 *
 *  Gaudi Consumer that validates vertex reconstruction against Monte Carlo
 *  truth. It only reads edm4hep collections, so it applies to any vertex
 *  producer (e.g. the ACTS adaptive multi-vertex finder or LCFIPlus) with the
 *  same definitions for all of them.
 *
 *  Truth matching follows ACTS' VertexTruthMatcher
 *  (Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/VertexTruthMatcher.cpp,
 *  acts-project/acts @ 9d258af1):
 *   - every track of a reconstructed vertex with weight >= MinTrackWeight
 *     contributes its weight to the truth vertex its MC particle was produced
 *     at (vertex -> particles -> tracks -> track-MC link -> MC particle);
 *   - the truth vertex with the largest summed weight is the match, and
 *     matchFraction = its weight / total weight of the reconstructed vertex;
 *   - matchFraction >= VertexMatchThreshold is "clean", otherwise "merged";
 *     a vertex without any truth-linked track is "unknown";
 *   - when two reconstructed vertices match the same truth vertex, the one with
 *     the smaller sum pT^2 is "split".
 *  Track weights are read from an optional vertex -> particle link collection
 *  (the adaptive fit's weights); without it every track counts with weight 1,
 *  which is exact for producers with a hard assignment such as LCFIPlus.
 *
 *  The quantities follow the ACTS vertexing performance writer
 *  (Examples/Io/Root/src/RootVertexNTupleWriter.cpp, same commit): one tree
 *  entry per event, per-vertex quantities stored as vectors, residual = reco -
 *  truth and pull = residual / sigma_reco.
 *
 *  Besides finding the truth primary vertex, a producer has to pick it as the
 *  primary one: the PV selection efficiency counts the events where the vertex
 *  it flags as primary is the clean match of the truth primary vertex. The
 *  tracks of the truth primary vertex (those whose MC particle was produced
 *  there) give truth-side quantities, their number and sum pT^2, that are the
 *  same for every producer run on the same tracks.
 *
 *  input:
 *    - MC particle collection         : edm4hep::MCParticleCollection
 *    - vertex collection              : edm4hep::VertexCollection
 *    - track -> MC particle links     : edm4hep::TrackMCParticleLinkCollection
 *    - optional vertex -> particle links with track weights
 *
 *  output:
 *    - ROOT file with the tree "vertex_vs_mc" and summary plots
 *
 */
struct VertexValidation final
    : k4FWCore::Consumer<void(const edm4hep::MCParticleCollection&, const edm4hep::VertexCollection&,
                              const edm4hep::TrackMCParticleLinkCollection&,
                              const std::vector<const edm4hep::VertexRecoParticleLinkCollection*>&)> {

  /// Reconstructed vertex classification, same values as ACTS'
  /// RecoVertexClassification
  enum Classification { Unknown = 0, Clean = 1, Merged = 2, Split = 3 };

  VertexValidation(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(name, svcLoc,
                 {
                     KeyValues("MCParticles", {"MCParticles"}),
                     KeyValues("Vertices", {"PrimaryVertices"}),
                     KeyValues("TrackMCLinks", {"SiTracksMCTruthLink"}),
                     KeyValues("VertexParticleLinks", {}),
                 }) {}

  StatusCode initialize() override {
    m_outFile = std::make_unique<TFile>(m_outputFile.value().c_str(), "RECREATE");
    if (!m_outFile || m_outFile->IsZombie()) {
      error() << "Cannot open output file: " << m_outputFile.value() << endmsg;
      return StatusCode::FAILURE;
    }
    m_outFile->cd();
    bookTree();
    return StatusCode::SUCCESS;
  }

  void
  operator()(const edm4hep::MCParticleCollection& mcParts, const edm4hep::VertexCollection& vertices,
             const edm4hep::TrackMCParticleLinkCollection& trackMCLinks,
             const std::vector<const edm4hep::VertexRecoParticleLinkCollection*>& vertexParticleLinks) const override {
    m_tree.clear();
    m_tree.event = m_evt++;

    // ---------- truth primary vertex ----------
    // The truth primary vertex is the production vertex of the particles
    // without parents (the generator's initial state, or the particles of a
    // particle gun). Vertex smearing in the simulation shifts all of them
    // together, so they share one position; distinct positions would mean
    // more than one collision in the event.
    std::vector<edm4hep::Vector3d> truthVertices; // index 0 is the truth primary vertex
    for (const auto& mc : mcParts) {
      if (!mc.getParents().empty()) {
        continue;
      }
      if (truthVertices.empty()) {
        truthVertices.push_back(mc.getVertex());
      } else if (!samePosition(mc.getVertex(), truthVertices.front())) {
        // A second primary position means a second collision (pile-up or
        // overlay). Everything is compared with the first one, so record that
        // the event is ambiguous rather than mixing silently.
        m_tree.nTrueVtx = 2;
      }
    }
    if (truthVertices.empty()) {
      if (!m_warnedNoTruth) {
        warning() << "No parentless MC particle found; cannot define the truth primary vertex" << endmsg;
        m_warnedNoTruth = true;
      }
      m_tree.tree->Fill();
      return;
    }
    if (m_tree.nTrueVtx == 0) {
      m_tree.nTrueVtx = 1;
    }
    // A copy: truthVertices grows during the matching below
    const edm4hep::Vector3d truthPV = truthVertices.front();
    m_tree.truthX = truthPV.x;
    m_tree.truthY = truthPV.y;
    m_tree.truthZ = truthPV.z;

    // Truth vertex index of a production position; secondary production points
    // (decays, interactions) become new truth vertices.
    auto truthVertexIndex = [&truthVertices, this](const edm4hep::Vector3d& pos) {
      for (std::size_t i = 0; i < truthVertices.size(); ++i) {
        if (samePosition(pos, truthVertices[i])) {
          return static_cast<int>(i);
        }
      }
      truthVertices.push_back(pos);
      return static_cast<int>(truthVertices.size() - 1);
    };

    // ---------- track -> MC particle ----------
    // A track can be linked to several MC particles (shared hits); keep the one
    // with the largest link weight.
    std::map<ObjectKey, std::pair<float, std::optional<edm4hep::MCParticle>>> trackTruth;
    for (const auto& link : trackMCLinks) {
      auto& best = trackTruth[key(link.getFrom())];
      if (!best.second.has_value() || link.getWeight() > best.first) {
        best = {link.getWeight(), link.getTo()};
      }
    }
    m_tree.nTracksEvent = static_cast<int>(trackTruth.size());

    // Tracks of the truth primary vertex. Performance plotted against their
    // number or sum pT^2 puts the same events in a bin for every producer, as
    // opposed to the producer's own track count. A duplicated track counts
    // twice, as the vertexer sees it twice.
    double truthPVSumPt2 = 0.;
    for (const auto& entry : trackTruth) {
      const edm4hep::MCParticle& mc = *entry.second.second;
      if (samePosition(mc.getVertex(), truthPV)) {
        ++m_tree.nTracksTruthPV;
        const auto& p = mc.getMomentum();
        truthPVSumPt2 += p.x * p.x + p.y * p.y;
      }
    }
    m_tree.truthPVSumPt2 = static_cast<float>(truthPVSumPt2);

    // ---------- track weights in the vertex fit ----------
    std::map<std::pair<ObjectKey, ObjectKey>, float> trackWeights;
    for (const auto* links : vertexParticleLinks) {
      for (const auto& link : *links) {
        trackWeights[{key(link.getFrom()), key(link.getTo())}] = link.getWeight();
      }
    }
    const bool haveTrackWeights = !vertexParticleLinks.empty();

    // ---------- reconstructed vertices ----------
    // The reconstructed primary vertex is the one the producer flagged as
    // primary; if none is flagged, the first one is taken for the PV plots.
    // Only a flagged vertex counts for the PV selection efficiency.
    int flaggedPVIndex = -1;
    int pvIndex = -1;
    m_tree.nRecoVtx = static_cast<int>(vertices.size());

    // truth vertex index -> (reco vertex index, sum pT^2) of its best match
    std::map<int, std::pair<std::size_t, double>> truthToReco;
    // per reco vertex: its tracks whose MC particle comes from the truth PV
    std::vector<int> tracksFromTruthPV;

    for (std::size_t i = 0; i < vertices.size(); ++i) {
      const auto vtx = vertices[i];
      const auto& pos = vtx.getPosition();
      const auto& cov = vtx.getCovMatrix();

      const float sigX = std::sqrt(cov.getValue(edm4hep::Cartesian::x, edm4hep::Cartesian::x));
      const float sigY = std::sqrt(cov.getValue(edm4hep::Cartesian::y, edm4hep::Cartesian::y));
      const float sigZ = std::sqrt(cov.getValue(edm4hep::Cartesian::z, edm4hep::Cartesian::z));

      // Residuals w.r.t. the truth primary vertex in micrometre, pulls dimensionless
      const float resX = (pos.x - truthPV.x) * 1000.f;
      const float resY = (pos.y - truthPV.y) * 1000.f;
      const float resZ = (pos.z - truthPV.z) * 1000.f;

      // --- truth matching, see the class documentation ---
      double totalWeight = 0.;
      double sumPt2 = 0.;
      int nTracks = 0;
      int nTruthLinkedTracks = 0;
      int nFromTruthPV = 0;
      std::map<int, double> weightPerTruthVertex;
      for (const auto& particle : vtx.getParticles()) {
        float weight = 1.f;
        if (haveTrackWeights) {
          const auto it = trackWeights.find({key(vtx), key(particle)});
          weight = it != trackWeights.end() ? it->second : 0.f;
        }
        if (weight < m_minTrackWeight) {
          continue;
        }
        totalWeight += weight;
        ++nTracks;
        const auto& p = particle.getMomentum();
        sumPt2 += p.x * p.x + p.y * p.y;

        for (const auto& track : particle.getTracks()) {
          const auto it = trackTruth.find(key(track));
          if (it != trackTruth.end()) {
            const int truthIndex = truthVertexIndex(it->second.second->getVertex());
            weightPerTruthVertex[truthIndex] += weight;
            ++nTruthLinkedTracks;
            if (truthIndex == 0) {
              ++nFromTruthPV;
            }
            break; // one truth particle per reconstructed particle
          }
        }
      }

      int matchedTruth = -1;
      double majorityWeight = 0.;
      for (const auto& [truthIndex, weight] : weightPerTruthVertex) {
        if (weight > majorityWeight) {
          majorityWeight = weight;
          matchedTruth = truthIndex;
        }
      }
      const double matchFraction = totalWeight > 0. ? majorityWeight / totalWeight : 0.;
      // Without any truth-linked track there is no truth vertex to match
      int classification = matchedTruth < 0 ? Unknown : matchFraction >= m_vertexMatchThreshold ? Clean : Merged;

      if (matchedTruth >= 0) {
        auto it = truthToReco.find(matchedTruth);
        if (it == truthToReco.end()) {
          truthToReco[matchedTruth] = {i, sumPt2};
        } else if (sumPt2 <= it->second.second) {
          // this truth vertex is already matched to a harder reco vertex
          classification = Split;
        } else {
          // This reco vertex is harder: it takes over and the other one is
          // split. As in ACTS, the taking-over vertex inherits the other one's
          // clean/merged classification rather than keeping its own.
          classification = static_cast<int>(m_tree.classification[it->second.first]);
          m_tree.classification[it->second.first] = static_cast<float>(Split);
          it->second = {i, sumPt2};
        }
      }

      m_tree.recoX.push_back(pos.x);
      m_tree.recoY.push_back(pos.y);
      m_tree.recoZ.push_back(pos.z);
      m_tree.sigmaX.push_back(sigX * 1000.f);
      m_tree.sigmaY.push_back(sigY * 1000.f);
      m_tree.sigmaZ.push_back(sigZ * 1000.f);
      m_tree.resX.push_back(resX);
      m_tree.resY.push_back(resY);
      m_tree.resZ.push_back(resZ);
      m_tree.pullX.push_back(pull(resX, sigX));
      m_tree.pullY.push_back(pull(resY, sigY));
      m_tree.pullZ.push_back(pull(resZ, sigZ));
      m_tree.chi2.push_back(vtx.getChi2());
      m_tree.ndf.push_back(static_cast<float>(vtx.getNdf()));
      m_tree.chi2ndf.push_back(vtx.getNdf() > 0 ? vtx.getChi2() / vtx.getNdf() : -1.f);
      m_tree.nTracks.push_back(static_cast<float>(nTracks));
      m_tree.nTruthLinkedTracks.push_back(static_cast<float>(nTruthLinkedTracks));
      m_tree.totalTrackWeight.push_back(static_cast<float>(totalWeight));
      m_tree.sumPt2.push_back(static_cast<float>(sumPt2));
      m_tree.matchFraction.push_back(static_cast<float>(matchFraction));
      m_tree.matchedTruthIndex.push_back(static_cast<float>(matchedTruth));
      m_tree.classification.push_back(static_cast<float>(classification));
      m_tree.isPrimary.push_back(vtx.isPrimary() ? 1.f : 0.f);
      m_tree.covPosDef.push_back(positiveDefinite(cov) ? 1.f : 0.f);
      tracksFromTruthPV.push_back(nFromTruthPV);

      if (vtx.isPrimary() && flaggedPVIndex < 0) {
        flaggedPVIndex = static_cast<int>(i);
      }
    }
    pvIndex = flaggedPVIndex >= 0 ? flaggedPVIndex : (vertices.empty() ? -1 : 0);

    for (const float c : m_tree.classification) {
      ++m_classificationCounts[static_cast<std::size_t>(c)];
    }

    // Found in the ACTS sense: the truth primary vertex is matched by a clean
    // reconstructed vertex.
    if (auto it = truthToReco.find(0);
        it != truthToReco.end() && static_cast<int>(m_tree.classification[it->second.first]) == Clean) {
      m_tree.truthPVFoundClean = 1;
      ++m_nEventsTruthPVClean;
      // ...and the producer flagged that vertex as the primary one
      if (static_cast<int>(it->second.first) == flaggedPVIndex) {
        m_tree.recoPVIsTruthPV = 1;
        ++m_nEventsRecoPVIsTruthPV;
      }
    }

    // The primary-vertex branches hold zero or one entry per event, so events
    // without a reconstructed primary vertex simply do not contribute.
    if (pvIndex >= 0) {
      m_tree.pvResX.push_back(m_tree.resX[pvIndex]);
      m_tree.pvResY.push_back(m_tree.resY[pvIndex]);
      m_tree.pvResZ.push_back(m_tree.resZ[pvIndex]);
      m_tree.pvPullX.push_back(m_tree.pullX[pvIndex]);
      m_tree.pvPullY.push_back(m_tree.pullY[pvIndex]);
      m_tree.pvPullZ.push_back(m_tree.pullZ[pvIndex]);
      m_tree.pvChi2Ndf.push_back(m_tree.chi2ndf[pvIndex]);
      m_tree.pvSigmaX.push_back(m_tree.sigmaX[pvIndex]);
      m_tree.pvSigmaY.push_back(m_tree.sigmaY[pvIndex]);
      m_tree.pvSigmaZ.push_back(m_tree.sigmaZ[pvIndex]);
      m_tree.pvNTracks.push_back(m_tree.nTracks[pvIndex]);
      m_tree.pvMatchFraction.push_back(m_tree.matchFraction[pvIndex]);
      if (m_tree.covPosDef[pvIndex] == 0.f) {
        ++m_nPVBadCovariance;
      }
      // Both counts are of tracks from the truth PV, so the fraction is at most
      // one; tracks from secondary vertices in the event do not lower it.
      if (m_tree.nTracksTruthPV > 0) {
        m_tree.pvTrackFraction.push_back(tracksFromTruthPV[pvIndex] / static_cast<float>(m_tree.nTracksTruthPV));
      }
      m_pvResiduals[0].push_back(m_tree.resX[pvIndex]);
      m_pvResiduals[1].push_back(m_tree.resY[pvIndex]);
      m_pvResiduals[2].push_back(m_tree.resZ[pvIndex]);
      ++m_nEventsWithRecoPV;
    }

    m_tree.tree->Fill();
  }

  StatusCode finalize() override {
    info() << "Finalizing VertexValidation, processed " << m_evt << " events, " << m_nEventsWithRecoPV
           << " with a reconstructed primary vertex, " << m_nEventsTruthPVClean
           << " with the truth primary vertex matched by a clean vertex, " << m_nEventsRecoPVIsTruthPV
           << " of them with that vertex flagged as primary" << endmsg;
    info() << "Reconstructed vertices: " << m_classificationCounts[Clean] << " clean, "
           << m_classificationCounts[Merged] << " merged, " << m_classificationCounts[Split] << " split, "
           << m_classificationCounts[Unknown] << " without truth match" << endmsg;
    if (m_nPVBadCovariance > 0) {
      warning() << m_nPVBadCovariance << " primary vertices have a covariance that is not positive definite; "
                << "their pulls are meaningless" << endmsg;
    }

    if (!m_outFile) {
      return StatusCode::SUCCESS;
    }
    m_outFile->cd();
    m_tree.tree->Write();

    // Headline numbers, collected below and written as a one-entry tree
    std::array<float, 3> sigmaEff{}, sigmaEffErr{}, residualCentre{}, pullMean{}, pullWidth{};

    // ---------- residuals of the primary vertex ----------
    const std::array<const char*, 3> axes{"X", "Y", "Z"};
    for (std::size_t k = 0; k < axes.size(); ++k) {
      const char* axis = axes[k];
      const std::string branch = std::string("pvRes") + axis;
      TH1F* h = TrackingValidationPlots::makePullHistogram(
          m_tree.tree, branch.c_str(), ("h_pv_res_" + std::string(axis)).c_str(),
          (std::string("PV residual ") + axis + ";" + axis + "_{reco} - " + axis + "_{true} [#mum];Events").c_str(),
          100, -m_residualRange, m_residualRange);
      if (!h) {
        continue;
      }
      // makePullHistogram labels the x axis with the branch name
      h->GetXaxis()->SetTitle((std::string(axis) + "_{reco} - " + axis + "_{true} [#mum]").c_str());

      // Robust width: half-width of the narrowest interval holding 68.27% of
      // the entries, which is insensitive to the non-Gaussian tails.
      const auto eff = TrackingValidationPlots::computeEffectiveSigma(m_pvResiduals[k]);
      const double effErr = TrackingValidationPlots::computeEffectiveSigmaBootstrapError(m_pvResiduals[k]);
      info() << "PV residual " << axis << ": sigma_eff = " << eff.sigmaEff << " +- " << effErr
             << " um, centre = " << eff.center << " um (" << eff.nEntries << " events)" << endmsg;
      sigmaEff[k] = static_cast<float>(eff.sigmaEff);
      sigmaEffErr[k] = static_cast<float>(effErr);
      residualCentre[k] = static_cast<float>(eff.center);

      auto* c = new TCanvas(("c_pv_res_" + std::string(axis)).c_str(), "", 900, 700);
      c->SetGrid();
      h->SetTitle(Form("PV residual %s   (#sigma_{eff} = %.2f #pm %.2f #mum)", axis, eff.sigmaEff, effErr));
      h->SetLineColor(kBlue + 1);
      h->SetLineWidth(2);
      h->Draw();
      h->Write();
      c->Write();
    }

    // ---------- pulls of the primary vertex ----------
    for (std::size_t k = 0; k < axes.size(); ++k) {
      const char* axis = axes[k];
      const std::string branch = std::string("pvPull") + axis;
      TH1F* h = TrackingValidationPlots::makePullHistogram(
          m_tree.tree, branch.c_str(), ("h_pv_pull_" + std::string(axis)).c_str(),
          (std::string("PV pull ") + axis + ";(" + axis + "_{reco} - " + axis + "_{true}) / #sigma_{" + axis +
           "};Events")
              .c_str(),
          100, -10.0, 10.0);
      if (h) {
        h->GetXaxis()->SetTitle(
            (std::string("(") + axis + "_{reco} - " + axis + "_{true}) / #sigma_{" + axis + "}").c_str());
      }
      TCanvas* c = TrackingValidationPlots::drawPullCanvas(h, ("c_pv_pull_" + std::string(axis)).c_str(),
                                                           (std::string("PV pull ") + axis).c_str());
      if (h) {
        info() << "PV pull " << axis << ": mean = " << h->GetMean() << ", RMS = " << h->GetRMS() << endmsg;
        // Core width from the Gaussian fit drawPullCanvas makes (it fits only
        // with enough entries); the histogram's mean and RMS otherwise.
        const TF1* fit = h->GetFunction((std::string(h->GetName()) + "_gaus").c_str());
        pullMean[k] = static_cast<float>(fit ? fit->GetParameter(1) : h->GetMean());
        pullWidth[k] = static_cast<float>(fit ? fit->GetParameter(2) : h->GetRMS());
        h->Write();
      }
      if (c) {
        c->Write();
      }
    }

    // ---------- fit quality, multiplicity and truth matching ----------
    writeHistogram("pvChi2Ndf", "h_pv_chi2ndf", "PV fit quality", "#chi^{2}/ndf", 100, 0.0, 10.0);
    writeHistogram("pvNTracks", "h_pv_ntracks", "Tracks used by the PV", "N_{tracks} (weight #geq threshold)", 31, -0.5,
                   30.5);
    writeHistogram("pvTrackFraction", "h_pv_track_fraction", "Fraction of the truth PV's tracks used by the PV",
                   "N_{tracks from truth PV, PV} / N_{tracks from truth PV, event}", 44, 0.0, 1.1);
    writeHistogram("matchFraction", "h_match_fraction", "Truth match fraction of reconstructed vertices",
                   "matched track weight / total track weight", 44, 0.0, 1.1);
    if (TH1F* h = writeHistogram("classification", "h_vertex_classification",
                                 "Classification of reconstructed vertices", "", 4, -0.5, 3.5)) {
      const std::array<const char*, 4> labels{"unknown", "clean", "merged", "split"};
      for (int bin = 1; bin <= 4; ++bin) {
        h->GetXaxis()->SetBinLabel(bin, labels[bin - 1]);
      }
      h->Write("", TObject::kOverwrite);
    }

    auto* hNVtx = new TH1F("h_n_reco_vertices", "Reconstructed vertices per event;N_{vertices};Events", 11, -0.5, 10.5);
    m_tree.tree->Project("h_n_reco_vertices", "nRecoVtx");
    hNVtx->Write();

    // ---------- PV efficiency versus the truth PV's track multiplicity ----------
    auto* hTotal =
        new TH1F("h_ntracks_truth_pv", "Reconstructed tracks from the truth PV;N_{tracks};Events", 51, -0.5, 50.5);
    auto* hPassed = new TH1F("h_ntracks_truth_pv_found", "", 51, -0.5, 50.5);
    m_tree.tree->Project("h_ntracks_truth_pv", "nTracksTruthPV");
    m_tree.tree->Project("h_ntracks_truth_pv_found", "nTracksTruthPV", "truthPVFoundClean == 1");
    hTotal->Write();
    auto* effVsNTracks = new TEfficiency(*hPassed, *hTotal);
    effVsNTracks->SetName("eff_pv_vs_ntracks_truth_pv");
    effVsNTracks->SetTitle("PV efficiency;Reconstructed tracks from the truth PV;Truth PV matched by a clean vertex");
    effVsNTracks->Write();

    // ---------- headline numbers ----------
    // One entry, so a notebook or script reads the numbers instead of
    // recomputing them. Residuals in um; the PV efficiency is the fraction of
    // events whose truth PV is matched by a clean reconstructed vertex, the PV
    // selection efficiency the fraction where that vertex is also the one
    // flagged as primary.
    auto* summary = new TTree("summary", "Headline numbers of the vertex validation");
    int nEvents = m_evt, nRecoPV = m_nEventsWithRecoPV, nTruthPVClean = m_nEventsTruthPVClean;
    int nRecoPVIsTruthPV = m_nEventsRecoPVIsTruthPV;
    int nClean = m_classificationCounts[Clean], nMerged = m_classificationCounts[Merged];
    int nSplit = m_classificationCounts[Split], nUnknown = m_classificationCounts[Unknown];
    int nPVBadCov = m_nPVBadCovariance;
    float pvEfficiency = nEvents > 0 ? static_cast<float>(nTruthPVClean) / nEvents : 0.f;
    float pvSelectionEfficiency = nEvents > 0 ? static_cast<float>(nRecoPVIsTruthPV) / nEvents : 0.f;
    summary->Branch("nEvents", &nEvents);
    summary->Branch("nEventsWithRecoPV", &nRecoPV);
    summary->Branch("nEventsTruthPVClean", &nTruthPVClean);
    summary->Branch("nEventsRecoPVIsTruthPV", &nRecoPVIsTruthPV);
    summary->Branch("pvEfficiency", &pvEfficiency);
    summary->Branch("pvSelectionEfficiency", &pvSelectionEfficiency);
    summary->Branch("nClean", &nClean);
    summary->Branch("nMerged", &nMerged);
    summary->Branch("nSplit", &nSplit);
    summary->Branch("nUnknown", &nUnknown);
    summary->Branch("nPVNotPositiveDefinite", &nPVBadCov);
    for (std::size_t k = 0; k < axes.size(); ++k) {
      const std::string axis = axes[k];
      summary->Branch(("sigmaEff" + axis).c_str(), &sigmaEff[k]);
      summary->Branch(("sigmaEffErr" + axis).c_str(), &sigmaEffErr[k]);
      summary->Branch(("residualCentre" + axis).c_str(), &residualCentre[k]);
      summary->Branch(("pullMean" + axis).c_str(), &pullMean[k]);
      summary->Branch(("pullWidth" + axis).c_str(), &pullWidth[k]);
    }
    summary->Fill();
    summary->Write();

    m_outFile->Close();
    return StatusCode::SUCCESS;
  }

private:
  using ObjectKey = std::pair<std::uint32_t, int>;

  template <typename T>
  static ObjectKey key(const T& obj) {
    const auto id = obj.getObjectID();
    return {id.collectionID, id.index};
  }

  struct VertexTree {
    TTree* tree = nullptr;

    int event = 0;
    int nRecoVtx = 0;
    int nTrueVtx = 0;
    int nTracksEvent = 0;                           // tracks with a truth link in the event
    int nTracksTruthPV = 0;                         // of them, tracks whose MC particle comes from the truth PV
    float truthPVSumPt2 = 0.f;                      // sum of their MC particles' pT^2 [GeV^2]
    int truthPVFoundClean = 0;                      // truth PV matched by a clean reco vertex
    int recoPVIsTruthPV = 0;                        // ... and that vertex is flagged as primary
    float truthX = 0.f, truthY = 0.f, truthZ = 0.f; // [mm]

    // one entry per reconstructed vertex
    std::vector<float> recoX, recoY, recoZ;    // [mm]
    std::vector<float> sigmaX, sigmaY, sigmaZ; // [um]
    std::vector<float> resX, resY, resZ;       // reco - truth PV [um]
    std::vector<float> pullX, pullY, pullZ;    // residual / sigma
    std::vector<float> chi2, ndf, chi2ndf;
    std::vector<float> nTracks;            // tracks with weight >= MinTrackWeight
    std::vector<float> nTruthLinkedTracks; // those of them with a truth link
    std::vector<float> totalTrackWeight;   // their summed weight
    std::vector<float> sumPt2;             // their sum pT^2 [GeV^2]
    std::vector<float> matchFraction;      // majority truth-vertex weight / total weight
    std::vector<float> matchedTruthIndex;  // 0 = truth PV, >0 secondary truth vertex, -1 none
    std::vector<float> classification;     // see Classification
    std::vector<float> isPrimary;          // 1 if flagged as primary by the producer
    std::vector<float> covPosDef;          // 1 if the position covariance is positive definite

    // zero or one entry: the reconstructed primary vertex
    std::vector<float> pvResX, pvResY, pvResZ; // [um]
    std::vector<float> pvPullX, pvPullY, pvPullZ;
    std::vector<float> pvSigmaX, pvSigmaY, pvSigmaZ; // [um]
    std::vector<float> pvChi2Ndf;
    std::vector<float> pvNTracks;
    std::vector<float> pvMatchFraction; // purity: matchFraction of the PV
    std::vector<float> pvTrackFraction; // PV tracks from the truth PV / nTracksTruthPV

    void clear() {
      nRecoVtx = nTrueVtx = nTracksEvent = nTracksTruthPV = truthPVFoundClean = recoPVIsTruthPV = 0;
      truthX = truthY = truthZ = truthPVSumPt2 = 0.f;
      for (auto* v : {&recoX,
                      &recoY,
                      &recoZ,
                      &sigmaX,
                      &sigmaY,
                      &sigmaZ,
                      &resX,
                      &resY,
                      &resZ,
                      &pullX,
                      &pullY,
                      &pullZ,
                      &chi2,
                      &ndf,
                      &chi2ndf,
                      &nTracks,
                      &nTruthLinkedTracks,
                      &totalTrackWeight,
                      &sumPt2,
                      &matchFraction,
                      &matchedTruthIndex,
                      &classification,
                      &isPrimary,
                      &covPosDef,
                      &pvResX,
                      &pvResY,
                      &pvResZ,
                      &pvPullX,
                      &pvPullY,
                      &pvPullZ,
                      &pvSigmaX,
                      &pvSigmaY,
                      &pvSigmaZ,
                      &pvChi2Ndf,
                      &pvNTracks,
                      &pvMatchFraction,
                      &pvTrackFraction}) {
        v->clear();
      }
    }
  };

  void bookTree() {
    auto& t = m_tree;
    t.tree = new TTree("vertex_vs_mc", "Reconstructed vertices compared with MC truth");
    t.tree->Branch("event", &t.event);
    t.tree->Branch("nRecoVtx", &t.nRecoVtx);
    t.tree->Branch("nTrueVtx", &t.nTrueVtx);
    t.tree->Branch("nTracksEvent", &t.nTracksEvent);
    t.tree->Branch("nTracksTruthPV", &t.nTracksTruthPV);
    t.tree->Branch("truthPVSumPt2", &t.truthPVSumPt2);
    t.tree->Branch("truthPVFoundClean", &t.truthPVFoundClean);
    t.tree->Branch("recoPVIsTruthPV", &t.recoPVIsTruthPV);
    t.tree->Branch("truthX", &t.truthX);
    t.tree->Branch("truthY", &t.truthY);
    t.tree->Branch("truthZ", &t.truthZ);
    t.tree->Branch("recoX", &t.recoX);
    t.tree->Branch("recoY", &t.recoY);
    t.tree->Branch("recoZ", &t.recoZ);
    t.tree->Branch("sigmaX", &t.sigmaX);
    t.tree->Branch("sigmaY", &t.sigmaY);
    t.tree->Branch("sigmaZ", &t.sigmaZ);
    t.tree->Branch("resX", &t.resX);
    t.tree->Branch("resY", &t.resY);
    t.tree->Branch("resZ", &t.resZ);
    t.tree->Branch("pullX", &t.pullX);
    t.tree->Branch("pullY", &t.pullY);
    t.tree->Branch("pullZ", &t.pullZ);
    t.tree->Branch("chi2", &t.chi2);
    t.tree->Branch("ndf", &t.ndf);
    t.tree->Branch("chi2ndf", &t.chi2ndf);
    t.tree->Branch("nTracks", &t.nTracks);
    t.tree->Branch("nTruthLinkedTracks", &t.nTruthLinkedTracks);
    t.tree->Branch("totalTrackWeight", &t.totalTrackWeight);
    t.tree->Branch("sumPt2", &t.sumPt2);
    t.tree->Branch("matchFraction", &t.matchFraction);
    t.tree->Branch("matchedTruthIndex", &t.matchedTruthIndex);
    t.tree->Branch("classification", &t.classification);
    t.tree->Branch("isPrimary", &t.isPrimary);
    t.tree->Branch("covPosDef", &t.covPosDef);
    t.tree->Branch("pvResX", &t.pvResX);
    t.tree->Branch("pvResY", &t.pvResY);
    t.tree->Branch("pvResZ", &t.pvResZ);
    t.tree->Branch("pvPullX", &t.pvPullX);
    t.tree->Branch("pvPullY", &t.pvPullY);
    t.tree->Branch("pvPullZ", &t.pvPullZ);
    t.tree->Branch("pvSigmaX", &t.pvSigmaX);
    t.tree->Branch("pvSigmaY", &t.pvSigmaY);
    t.tree->Branch("pvSigmaZ", &t.pvSigmaZ);
    t.tree->Branch("pvChi2Ndf", &t.pvChi2Ndf);
    t.tree->Branch("pvNTracks", &t.pvNTracks);
    t.tree->Branch("pvMatchFraction", &t.pvMatchFraction);
    t.tree->Branch("pvTrackFraction", &t.pvTrackFraction);
  }

  /// Sylvester's criterion: all leading principal minors positive. A vertex
  /// whose covariance fails this has no meaningful uncertainty, whatever its
  /// position.
  static bool positiveDefinite(const edm4hep::CovMatrix3f& c) {
    using edm4hep::Cartesian;
    const double xx = c.getValue(Cartesian::x, Cartesian::x), yy = c.getValue(Cartesian::y, Cartesian::y),
                 zz = c.getValue(Cartesian::z, Cartesian::z), xy = c.getValue(Cartesian::x, Cartesian::y),
                 xz = c.getValue(Cartesian::x, Cartesian::z), yz = c.getValue(Cartesian::y, Cartesian::z);
    const double minor2 = xx * yy - xy * xy;
    const double minor3 = xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz);
    return std::isfinite(minor3) && xx > 0. && minor2 > 0. && minor3 > 0.;
  }

  bool samePosition(const edm4hep::Vector3d& a, const edm4hep::Vector3d& b) const {
    return std::abs(a.x - b.x) <= m_samePositionTolerance && std::abs(a.y - b.y) <= m_samePositionTolerance &&
           std::abs(a.z - b.z) <= m_samePositionTolerance;
  }

  static float pull(float residualUm, float sigmaMm) {
    return sigmaMm > 0.f ? residualUm / (sigmaMm * 1000.f) : std::numeric_limits<float>::quiet_NaN();
  }

  /// Histogram a vector<float> branch of the tree and write it.
  TH1F* writeHistogram(const char* branch, const char* name, const char* title, const char* xTitle, int nBins,
                       double xMin, double xMax) const {
    TH1F* h = TrackingValidationPlots::makePullHistogram(m_tree.tree, branch, name, title, nBins, xMin, xMax);
    if (h) {
      h->GetXaxis()->SetTitle(xTitle);
      h->GetYaxis()->SetTitle("Entries");
      h->Write();
    }
    return h;
  }

  // ---------- properties ----------
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "vertex_validation.root",
                                            "Output ROOT file (tree and plots)"};
  Gaudi::Property<float> m_residualRange{this, "ResidualRange", 10.f, "Half-range of the PV residual histograms [um]"};
  Gaudi::Property<double> m_samePositionTolerance{
      this, "TruthVertexTolerance", 1e-6,
      "MC production points closer than this are taken to be the same truth vertex [mm]"};
  Gaudi::Property<float> m_minTrackWeight{this, "MinTrackWeight", 0.1f,
                                          "Minimum track weight for a track to count as part of a vertex "
                                          "(ACTS VertexTruthMatcher minTrkWeight)"};
  Gaudi::Property<double> m_vertexMatchThreshold{this, "VertexMatchThreshold", 0.7,
                                                 "Minimum match fraction for a clean vertex "
                                                 "(ACTS VertexTruthMatcher vertexMatchThreshold)"};

  // ---------- state ----------
  std::unique_ptr<TFile> m_outFile;
  mutable VertexTree m_tree;
  mutable int m_evt = 0;
  mutable int m_nEventsWithRecoPV = 0;
  mutable int m_nEventsTruthPVClean = 0;
  mutable int m_nEventsRecoPVIsTruthPV = 0;
  mutable int m_nPVBadCovariance = 0;
  mutable std::array<int, 4> m_classificationCounts{};
  mutable bool m_warnedNoTruth = false;
  /// PV residuals x, y, z [um] over all events, for the effective sigma
  mutable std::array<std::vector<double>, 3> m_pvResiduals;
};

DECLARE_COMPONENT(VertexValidation)
