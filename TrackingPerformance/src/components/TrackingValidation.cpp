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
#include "TrackingValidationPlots.h"
#include "podio/ObjectID.h"

// k4FWCore
#include "k4FWCore/Consumer.h"

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiKernel/MsgStream.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"

// STL
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

/** @struct TrackingValidation
 *
 *  Gaudi Consumer that validates the performance of track finding and track fitting
 *  by comparing reconstructed tracks with Monte Carlo truth information and,
 *  optionally, with perfectly associated fitted tracks.
 *
 *  The consumer writes several ROOT TTrees containing finder-level associations,
 *  fitter residuals with respect to MC truth, and fitter residuals with respect
 *  to perfectly fitted reference tracks. In addition, summary performance plots
 *  are produced in finalize() and written to the same ROOT file.
 *
 *  The supported validation modes are:
 *    - full pipeline validation,
 *    - finder-only validation,
 *    - fitter-only validation.
 *
 *  input:
 *    - MC particle collection : edm4hep::MCParticleCollection
 *    - digi-to-sim link collections : std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>
 *    - finder track collection : edm4hep::TrackCollection
 *    - fitted track collection : edm4hep::TrackCollection
 *    - optional perfect fitted-track collections : std::vector<const edm4hep::TrackCollection*>
 *
 *  output:
 *    - ROOT file containing validation TTrees
 *    - summary performance plots written to the same ROOT file
 *
 *  @author Arina Ponomareva
 *  @date   2026-03
 *
 */

// ---------- CONSUMER ----------
struct TrackingValidation final
    : k4FWCore::Consumer<void(
          const edm4hep::MCParticleCollection&,
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&, // all digi-to-sim link
                                                                                     // collectionns
          const std::vector<const edm4hep::TrackCollection*>&,                       // finder tracks
          const std::vector<const edm4hep::TrackCollection*>&,                       // fitted tracks (reco)
          const std::vector<const edm4hep::TrackCollection*>&                        // optional perfect fitted tracks
          )> {

  TrackingValidation(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(name, svcLoc,
                 {
                     KeyValues("MCParticles", {"MCParticles"}),

                     KeyValues("HitSimLinks",
                               {"SiWrBSimDigiLinks", "SiWrDSimDigiLinks", "VTXBSimDigiLinks", "VTXDSimDigiLinks"}),

                     KeyValues("FinderTracks", {"GGTFTracks"}),
                     KeyValues("FittedTracks", {"FittedTracks"}),

                     KeyValues("PerfectFittedTracks", {"PerfectFittedTracks"}),
                 }) {}

  StatusCode initialize() override {
    info() << "Initializing TrackingValidationConsumer" << endmsg;

    m_outFile = std::make_unique<TFile>(m_outputFile.value().c_str(), "RECREATE");
    if (!m_outFile || m_outFile->IsZombie()) {
      error() << "Cannot open output file: " << m_outputFile.value() << endmsg;
      return StatusCode::FAILURE;
    }

    bookAssocTree(m_finder_p2t, "finder_particle_to_tracks");
    bookAssocTree(m_finder_t2p, "finder_track_to_particles");
    bookAssocTree(m_perf_p2t, "perfect_particle_to_tracks");
    bookAssocTree(m_perf_t2p, "perfect_track_to_particles");

    bookFitterTree(m_fit_vs_mc, "fitter_vs_mc");
    bookFitterTree(m_fit_vs_perfect, "fitter_vs_perfect");

    return StatusCode::SUCCESS;
  }

  void operator()(const edm4hep::MCParticleCollection& mcParts,
                  const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& linkCollections,
                  const std::vector<const edm4hep::TrackCollection*>& finderTracksVec,
                  const std::vector<const edm4hep::TrackCollection*>& fittedTracksVec,
                  const std::vector<const edm4hep::TrackCollection*>& perfectFittedTracksVec) const override {

    const int event = m_evt++;
    const int mode = m_mode.value(); // 0 full, 1 finder-only, 2 fitter-only

    const edm4hep::TrackCollection* finderTracks = finderTracksVec.empty() ? nullptr : finderTracksVec.front();

    const edm4hep::TrackCollection* fittedTracks = fittedTracksVec.empty() ? nullptr : fittedTracksVec.front();

    const bool needFinderTracks = (mode == 0 || mode == 1);
    const bool needFittedTracks = (mode == 0 || mode == 2);

    if (needFinderTracks && !finderTracks && !m_warnedMissingFinderInput) {
      warning() << "FinderTracks input is empty, but Mode=" << mode
                << " requires finder-track validation. Finder association trees will not be filled." << endmsg;
      m_warnedMissingFinderInput = true;
    }

    if (needFittedTracks && !fittedTracks && !m_warnedMissingFittedInput) {
      warning() << "FittedTracks input is empty, but Mode=" << mode
                << " requires fitter validation. Fitter residual trees will not be filled." << endmsg;
      m_warnedMissingFittedInput = true;
    }

    // ---------- Build truth maps: hit -> particle, particle -> hits ----------
    std::unordered_map<int, std::vector<podio::ObjectID>> hitsPerParticle;
    hitsPerParticle.reserve(mcParts.size());

    std::unordered_map<podio::ObjectID, int> hitToParticle;
    hitToParticle.reserve(200000);

    // digi-to-sim links
    for (const auto* links : linkCollections) {
      if (!links)
        continue;
      for (const auto& link : *links) {
        const auto digi = link.getFrom();
        const auto sim = link.getTo();
        const auto mc = sim.getParticle();
        if (!digi.isAvailable() || !mc.isAvailable())
          continue;

        const int pid = mc.getObjectID().index;
        const auto key = digi.getObjectID();
        hitsPerParticle[pid].push_back(key);
        hitToParticle[key] = pid;
      }
    }

    // ---------- Finder & Perfect association trees ----------
    if ((mode == 0 || mode == 1) && finderTracks) {
      fillPerfectAssoc(event, mcParts, hitsPerParticle);
      fillFinderAssoc(event, mcParts, *finderTracks, hitToParticle, hitsPerParticle);
    }

    // ---------- Build pid -> best perfect-fitted AtIP state ----------
    struct StateWithNHits {
      edm4hep::TrackState st;
      int nHits = 0;
    };

    std::unordered_map<int, StateWithNHits> perfectAtIPByPid;

    const bool wantPerfect = m_doPerfectFit.value();
    const bool havePerfectCollections = !perfectFittedTracksVec.empty();
    const bool doPerfect = wantPerfect && havePerfectCollections;

    if (wantPerfect && !havePerfectCollections && !m_warnedMissingPerfectInput) {
      warning() << "DoPerfectFit=true but no PerfectFittedTracks collection was provided. "
                << "fitter_vs_perfect will be filled with NaNs/empty content." << endmsg;
      m_warnedMissingPerfectInput = true;
    }

    if (doPerfect) {
      size_t nPerfectTracks = 0;
      for (const auto* coll : perfectFittedTracksVec) {
        if (!coll)
          continue;
        nPerfectTracks += coll->size();
      }
      perfectAtIPByPid.reserve(nPerfectTracks);

      for (const auto* coll : perfectFittedTracksVec) {
        if (!coll)
          continue;

        for (const auto& trk : *coll) {
          auto st = TrackingValidationHelpers::getAtIPState(trk);
          if (!st)
            continue;

          const int pid = majorityParticleForTrack(trk, hitToParticle);
          if (pid < 0 || pid >= (int)mcParts.size())
            continue;

          const int nHits = (int)trk.getTrackerHits().size();
          auto it = perfectAtIPByPid.find(pid);
          if (it == perfectAtIPByPid.end() || nHits > it->second.nHits) {
            perfectAtIPByPid[pid] = StateWithNHits{*st, nHits};
          }
        }
      }
    }

    // ---------- Fitter trees ----------
    if ((mode == 0 || mode == 2) && fittedTracks) {
      fillFitterTrees(event, mcParts, *fittedTracks, hitToParticle, perfectAtIPByPid, doPerfect);
    }
  }

  StatusCode finalize() override {
    info() << "Finalizing TrackingValidation, wrote " << m_evt << " events" << endmsg;

    if (m_outFile) {
      m_outFile->cd();

      // write trees
      if (m_finder_p2t.tree)
        m_finder_p2t.tree->Write();
      if (m_finder_t2p.tree)
        m_finder_t2p.tree->Write();
      if (m_perf_p2t.tree)
        m_perf_p2t.tree->Write();
      if (m_perf_t2p.tree)
        m_perf_t2p.tree->Write();

      if (m_fit_vs_mc.tree)
        m_fit_vs_mc.tree->Write();
      if (m_fit_vs_perfect.tree)
        m_fit_vs_perfect.tree->Write();

      // fitter summary plots
      // d0, z0, phi, omega, tanLambda  resolution vs momentum from fitter_vs_mc

      TGraphErrors* g_d0_vs_p = TrackingValidationPlots::makeD0ResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_d0_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_d0_vs_p) {
        TCanvas* c_d0_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_d0_vs_p, "c_d0_resolution_vs_p", "d0 resolution vs momentum;p_{ref} [GeV];#sigma(d_{0}) [#mum]", 0.1,
            100.0);
        g_d0_vs_p->Write();
        if (c_d0_vs_p)
          c_d0_vs_p->Write();
      }

      TGraphErrors* g_z0_vs_p = TrackingValidationPlots::makeZ0ResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_z0_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_z0_vs_p) {
        TCanvas* c_z0_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_z0_vs_p, "c_z0_resolution_vs_p", "z0 resolution vs momentum;p_{ref} [GeV];#sigma(z_{0}) [#mum]", 0.1,
            100.0);
        g_z0_vs_p->Write();
        if (c_z0_vs_p)
          c_z0_vs_p->Write();
      }

      TGraphErrors* g_phi_vs_p = TrackingValidationPlots::makePhiResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_phi_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_phi_vs_p) {
        TCanvas* c_phi_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_phi_vs_p, "c_phi_resolution_vs_p", "phi resolution vs momentum;p_{ref} [GeV];#sigma(#phi) [rad]", 0.1,
            100.0);
        g_phi_vs_p->Write();
        if (c_phi_vs_p)
          c_phi_vs_p->Write();
      }

      TGraphErrors* g_omega_vs_p = TrackingValidationPlots::makeOmegaResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_omega_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_omega_vs_p) {
        TCanvas* c_omega_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_omega_vs_p, "c_omega_resolution_vs_p", "omega resolution vs momentum;p_{ref} [GeV];#sigma(#omega) [1/mm]",
            0.1, 100.0);
        g_omega_vs_p->Write();
        if (c_omega_vs_p)
          c_omega_vs_p->Write();
      }

      TGraphErrors* g_tanl_vs_p = TrackingValidationPlots::makeTanLambdaResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_tanlambda_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_tanl_vs_p) {
        TCanvas* c_tanl_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_tanl_vs_p, "c_tanlambda_resolution_vs_p",
            "tanLambda resolution vs momentum;p_{ref} [GeV];#sigma(tan#lambda)", 0.1, 100.0);
        g_tanl_vs_p->Write();
        if (c_tanl_vs_p)
          c_tanl_vs_p->Write();
      }
      const std::vector<std::tuple<const char*, const char*, const char*>> pullPlots = {
          {"pullD0", "h_pull_d0", "d0 pull plot;pull(d_{0});Entries"},
          {"pullZ0", "h_pull_z0", "z0 pull plot;pull(z_{0});Entries"},
          {"pullPhi", "h_pull_phi", "phi pull plot;pull(#phi);Entries"},
          {"pullOmega", "h_pull_omega", "omega pull plot;pull(#omega);Entries"},
          {"pullTanLambda", "h_pull_tanlambda", "tanLambda pull plot;pull(tan#lambda);Entries"},
      };

      for (const auto& [branchName, histName, title] : pullPlots) {
        TH1F* h = TrackingValidationPlots::makePullHistogram(m_fit_vs_mc.tree, branchName, histName, title);
        if (h) {
          std::string canvasName = std::string("c_") + histName + "_vs_mc";
          TCanvas* c = TrackingValidationPlots::drawPullCanvas(h, canvasName.c_str(), title);
          h->Write();
          if (c)
            c->Write();
        }
      }

      // p resolution vs momentum
      TGraphErrors* g_p_vs_p = TrackingValidationPlots::makeMomentumResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_p_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_p_vs_p) {
        TCanvas* c_p_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_p_vs_p, "c_p_resolution_vs_p",
            "momentum resolution vs momentum;p_{ref} [GeV];#sigma((p_{reco}-p_{ref})/p_{ref})", 0.1, 100.0);
        g_p_vs_p->Write();
        if (c_p_vs_p)
          c_p_vs_p->Write();
      }

      // pT resolution vs momentum
      TGraphErrors* g_pt_vs_p = TrackingValidationPlots::makePtResolutionVsMomentum(
          m_fit_vs_mc.tree, "g_pt_resolution_vs_p", 0.1, 100.0, 0.15);
      if (g_pt_vs_p) {
        TCanvas* c_pt_vs_p = TrackingValidationPlots::drawResolutionCanvas(
            g_pt_vs_p, "c_pt_resolution_vs_p",
            "pT resolution vs momentum;p_{ref} [GeV];#sigma((pT_{reco}-pT_{ref})/pT_{ref})", 0.1, 100.0);
        g_pt_vs_p->Write();
        if (c_pt_vs_p)
          c_pt_vs_p->Write();
      }

      // finder summary plot
      TGraphErrors* g_eff_vs_p = TrackingValidationPlots::makeEfficiencyVsMomentum(
          m_finder_p2t.tree, "g_efficiency_vs_p", m_finderEfficiencyDefinition.value(), m_finderPurityThreshold.value(),
          0.1, 100.0, 0.15);
      if (g_eff_vs_p) {
        TCanvas* c_eff_vs_p = TrackingValidationPlots::drawEfficiencyCanvas(
            g_eff_vs_p, "c_efficiency_vs_p", "tracking efficiency vs momentum;p [GeV];Efficiency", 0.1, 100.0);
        g_eff_vs_p->Write();
        if (c_eff_vs_p)
          c_eff_vs_p->Write();
      }

      m_outFile->Close();
    }
    return StatusCode::SUCCESS;
  }

private:
  // ---------- properties ----------
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "validation.root", "Output ROOT file (TTrees)"};

  // 0 full, 1 finder-only, 2 fitter-only
  Gaudi::Property<int> m_mode{this, "Mode", 0, "Validation mode: 0=Full, 1=FinderOnly, 2=FitterOnly"};

  Gaudi::Property<float> m_Bz{this, "Bz", 2.f, "Magnetic field Bz [T] used in omega convention"};

  Gaudi::Property<float> m_refX{this, "RefPointX", 0.f,
                                "Reference point X [mm] (must match fitter m_VP_referencePoint)"};
  Gaudi::Property<float> m_refY{this, "RefPointY", 0.f,
                                "Reference point Y [mm] (must match fitter m_VP_referencePoint)"};
  Gaudi::Property<float> m_refZ{this, "RefPointZ", 0.f,
                                "Reference point Z [mm] (must match fitter m_VP_referencePoint)"};

  // Definition used for the summary tracking-efficiency plot.
  // Default = 1 keeps the current behaviour unchanged.
  Gaudi::Property<int> m_finderEfficiencyDefinition{this, "FinderEfficiencyDefinition", 1,
                                                    "Definition used for the tracking-efficiency summary plot: "
                                                    "1 = require purity >= FinderPurityThreshold; "
                                                    "2 = require purity >= 0.5 and efficiency >= 0.5"};

  Gaudi::Property<float> m_finderPurityThreshold{
      this, "FinderPurityThreshold", 0.75f,
      "Minimum purity for a particle-track match to count in tracking efficiency "
      "when FinderEfficiencyDefinition = 1"};

  Gaudi::Property<bool> m_doPerfectFit{this, "DoPerfectFit", false,
                                       "If true: fill fitter_vs_perfect using PerfectFitted_tracks if available. "
                                       "If false: tree exists but is empty per event."};

  // ---------- output structs ----------
  struct AssocTree {
    TTree* tree = nullptr;
    int event = 0;

    std::vector<int> index;

    // MC-particle truth information stored for optional downstream selections.
    // These values are meaningful for particle -> track trees.
    // For track -> particle trees, placeholder values are filled.
    std::vector<float> p;
    std::vector<float> pT;
    std::vector<float> theta;   // polar angle [rad]
    std::vector<float> vertexR; // production vertex radius sqrt(x^2 + y^2) [mm]
    std::vector<float> vertexZ; // production vertex z [mm]
    std::vector<float> charge;
    std::vector<int> pdg;

    std::vector<int> nTrueHits;
    std::vector<std::vector<int>> assoc;

    // only really used for finder_particle_to_tracks
    std::vector<std::vector<int>> sharedHits;
    std::vector<std::vector<float>> matchEfficiency;
    std::vector<std::vector<float>> matchPurity;

    void clear() {
      index.clear();

      p.clear();
      pT.clear();
      theta.clear();
      vertexR.clear();
      vertexZ.clear();
      charge.clear();
      pdg.clear();

      nTrueHits.clear();
      assoc.clear();
      sharedHits.clear();
      matchEfficiency.clear();
      matchPurity.clear();
    }
  };

  struct FitterTree {
    TTree* tree = nullptr;
    int event = 0;

    std::vector<int> track_index;
    std::vector<int> track_location;

    std::vector<float> resD0, resZ0, resPhi, resOmega, resTanL;
    std::vector<float> pullD0, pullZ0, pullPhi, pullOmega, pullTanL;
    std::vector<float> errD0, errZ0, errPhi, errOmega, errTanL;
    std::vector<float> p_reco, p_ref;
    std::vector<float> pT_reco, pT_ref;

    void clear() {
      track_index.clear();
      track_location.clear();
      resD0.clear();
      resZ0.clear();
      resPhi.clear();
      resOmega.clear();
      resTanL.clear();
      pullD0.clear();
      pullZ0.clear();
      pullPhi.clear();
      pullOmega.clear();
      pullTanL.clear();
      errD0.clear();
      errZ0.clear();
      errPhi.clear();
      errOmega.clear();
      errTanL.clear();
      p_reco.clear();
      p_ref.clear();
      pT_reco.clear();
      pT_ref.clear();
    }
  };

  static void bookAssocTree(AssocTree& t, const char* name) {
    t.tree = new TTree(name, name);
    t.tree->Branch("event", &t.event);
    t.tree->Branch("index", &t.index);
    t.tree->Branch("p", &t.p);
    t.tree->Branch("pT", &t.pT);
    t.tree->Branch("theta", &t.theta);
    t.tree->Branch("vertexR", &t.vertexR);
    t.tree->Branch("vertexZ", &t.vertexZ);
    t.tree->Branch("charge", &t.charge);
    t.tree->Branch("pdg", &t.pdg);
    t.tree->Branch("nTrueHits", &t.nTrueHits);
    t.tree->Branch("assoc", &t.assoc);
    t.tree->Branch("sharedHits", &t.sharedHits);
    t.tree->Branch("matchEfficiency", &t.matchEfficiency);
    t.tree->Branch("matchPurity", &t.matchPurity);
  }

  static void bookFitterTree(FitterTree& t, const char* name) {
    t.tree = new TTree(name, name);
    t.tree->Branch("event", &t.event);
    t.tree->Branch("track_index", &t.track_index);
    t.tree->Branch("track_location", &t.track_location);
    t.tree->Branch("resD0", &t.resD0);
    t.tree->Branch("resZ0", &t.resZ0);
    t.tree->Branch("resPhi", &t.resPhi);
    t.tree->Branch("resOmega", &t.resOmega);
    t.tree->Branch("resTanLambda", &t.resTanL);
    t.tree->Branch("pullD0", &t.pullD0);
    t.tree->Branch("pullZ0", &t.pullZ0);
    t.tree->Branch("pullPhi", &t.pullPhi);
    t.tree->Branch("pullOmega", &t.pullOmega);
    t.tree->Branch("pullTanLambda", &t.pullTanL);
    t.tree->Branch("errD0", &t.errD0);
    t.tree->Branch("errZ0", &t.errZ0);
    t.tree->Branch("errPhi", &t.errPhi);
    t.tree->Branch("errOmega", &t.errOmega);
    t.tree->Branch("errTanLambda", &t.errTanL);
    t.tree->Branch("p_reco", &t.p_reco);
    t.tree->Branch("p_ref", &t.p_ref);
    t.tree->Branch("pT_reco", &t.pT_reco);
    t.tree->Branch("pT_ref", &t.pT_ref);
  }

  // ---------- association trees ----------
  void fillPerfectAssoc(int event, const edm4hep::MCParticleCollection& mcParts,
                        const std::unordered_map<int, std::vector<podio::ObjectID>>& hitsPerParticle) const {

    m_perf_p2t.clear();
    m_perf_t2p.clear();
    m_perf_p2t.event = event;
    m_perf_t2p.event = event;

    for (int i = 0; i < (int)mcParts.size(); ++i) {
      const auto& mc = mcParts[i];
      if (mc.getGeneratorStatus() != 1)
        continue;

      auto it = hitsPerParticle.find(i);
      if (it == hitsPerParticle.end() || it->second.empty())
        continue;

      const auto& mom = mc.getMomentum();
      const float px = float(mom.x);
      const float py = float(mom.y);
      const float pz = float(mom.z);
      const float p = std::sqrt(px * px + py * py + pz * pz);
      const float pT = std::sqrt(px * px + py * py);
      const int nHits = (int)it->second.size();

      const auto& vtx = mc.getVertex();
      const float theta = std::atan2(pT, pz);
      const float vertexR = std::sqrt(float(vtx.x) * float(vtx.x) + float(vtx.y) * float(vtx.y));
      const float vertexZ = float(vtx.z);
      const float charge = float(mc.getCharge());
      const int pdg = mc.getPDG();

      m_perf_p2t.index.push_back(i);
      m_perf_p2t.p.push_back(p);
      m_perf_p2t.pT.push_back(pT);
      m_perf_p2t.theta.push_back(theta);
      m_perf_p2t.vertexR.push_back(vertexR);
      m_perf_p2t.vertexZ.push_back(vertexZ);
      m_perf_p2t.charge.push_back(charge);
      m_perf_p2t.pdg.push_back(pdg);
      m_perf_p2t.nTrueHits.push_back(nHits);
      m_perf_p2t.assoc.push_back({i});
      m_perf_p2t.sharedHits.push_back({});
      m_perf_p2t.matchEfficiency.push_back({});
      m_perf_p2t.matchPurity.push_back({});

      m_perf_t2p.index.push_back(i);
      m_perf_t2p.p.push_back(p);
      m_perf_t2p.pT.push_back(pT);
      m_perf_t2p.theta.push_back(theta);
      m_perf_t2p.vertexR.push_back(vertexR);
      m_perf_t2p.vertexZ.push_back(vertexZ);
      m_perf_t2p.charge.push_back(charge);
      m_perf_t2p.pdg.push_back(pdg);
      m_perf_t2p.nTrueHits.push_back(nHits);
      m_perf_t2p.assoc.push_back({i});
      m_perf_t2p.sharedHits.push_back({});
      m_perf_t2p.matchEfficiency.push_back({});
      m_perf_t2p.matchPurity.push_back({});
    }

    if (m_perf_p2t.tree)
      m_perf_p2t.tree->Fill();
    if (m_perf_t2p.tree)
      m_perf_t2p.tree->Fill();
  }

  void fillFinderAssoc(int event, const edm4hep::MCParticleCollection& mcParts,
                       const edm4hep::TrackCollection& finderTracks,
                       const std::unordered_map<podio::ObjectID, int>& hitToParticle,
                       const std::unordered_map<int, std::vector<podio::ObjectID>>& hitsPerParticle) const {

    m_finder_p2t.clear();
    m_finder_t2p.clear();
    m_finder_p2t.event = event;
    m_finder_t2p.event = event;

    std::vector<std::unordered_map<int, int>> trackParticleCounts;
    std::vector<int> trackNHits;
    trackParticleCounts.resize(finderTracks.size());
    trackNHits.resize(finderTracks.size(), 0);

    int tIdx = 0;
    for (const auto& trk : finderTracks) {
      trackNHits[tIdx] = (int)trk.getTrackerHits().size();
      for (const auto& h : trk.getTrackerHits()) {
        const auto hk = h.getObjectID();
        auto it = hitToParticle.find(hk);
        if (it == hitToParticle.end())
          continue;
        trackParticleCounts[tIdx][it->second] += 1;
      }
      ++tIdx;
    }

    // track -> particles
    for (int t = 0; t < (int)finderTracks.size(); ++t) {
      m_finder_t2p.index.push_back(t);
      m_finder_t2p.p.push_back(-1.f);
      m_finder_t2p.pT.push_back(-1.f);
      m_finder_t2p.theta.push_back(-1.f);
      m_finder_t2p.vertexR.push_back(-1.f);
      m_finder_t2p.vertexZ.push_back(-1.f);
      m_finder_t2p.charge.push_back(0.f);
      m_finder_t2p.pdg.push_back(0);
      m_finder_t2p.nTrueHits.push_back(trackNHits[t]);

      std::vector<int> parts;
      std::vector<int> sh;
      std::vector<float> effs;
      std::vector<float> purs;

      for (const auto& kv : trackParticleCounts[t]) {
        const int pid = kv.first;
        if (pid < 0)
          continue;

        const int shared = kv.second;
        const int nTrackHits = trackNHits[t];
        const int nParticleHits = hitsPerParticle.count(pid) ? (int)hitsPerParticle.at(pid).size() : 0;

        const float eff = (nParticleHits > 0) ? float(shared) / float(nParticleHits) : 0.f;
        const float pur = (nTrackHits > 0) ? float(shared) / float(nTrackHits) : 0.f;

        parts.push_back(pid);
        sh.push_back(shared);
        effs.push_back(eff);
        purs.push_back(pur);
      }

      m_finder_t2p.assoc.push_back(parts);
      m_finder_t2p.sharedHits.push_back(sh);
      m_finder_t2p.matchEfficiency.push_back(effs);
      m_finder_t2p.matchPurity.push_back(purs);
    }

    // particle -> tracks
    for (int p = 0; p < (int)mcParts.size(); ++p) {
      const auto& mc = mcParts[p];
      if (mc.getGeneratorStatus() != 1)
        continue;

      auto itHits = hitsPerParticle.find(p);
      if (itHits == hitsPerParticle.end() || itHits->second.empty())
        continue;

      const auto& mom = mc.getMomentum();
      const float px = float(mom.x);
      const float py = float(mom.y);
      const float pz = float(mom.z);
      const float pAbs = std::sqrt(px * px + py * py + pz * pz);
      const float pT = std::sqrt(px * px + py * py);
      const int nParticleHits = (int)itHits->second.size();

      const auto& vtx = mc.getVertex();
      const float theta = std::atan2(pT, pz);
      const float vertexR = std::sqrt(float(vtx.x) * float(vtx.x) + float(vtx.y) * float(vtx.y));
      const float vertexZ = float(vtx.z);
      const float charge = float(mc.getCharge());
      const int pdg = mc.getPDG();

      std::vector<int> tracks;
      std::vector<int> sh;
      std::vector<float> effs;
      std::vector<float> purs;

      for (int t = 0; t < (int)finderTracks.size(); ++t) {
        auto it = trackParticleCounts[t].find(p);
        if (it == trackParticleCounts[t].end())
          continue;

        const int shared = it->second;
        const int nTrackHits = trackNHits[t];

        const float eff = (nParticleHits > 0) ? float(shared) / float(nParticleHits) : 0.f;
        const float pur = (nTrackHits > 0) ? float(shared) / float(nTrackHits) : 0.f;

        tracks.push_back(t);
        sh.push_back(shared);
        effs.push_back(eff);
        purs.push_back(pur);
      }

      m_finder_p2t.index.push_back(p);
      m_finder_p2t.p.push_back(pAbs);
      m_finder_p2t.pT.push_back(pT);
      m_finder_p2t.theta.push_back(theta);
      m_finder_p2t.vertexR.push_back(vertexR);
      m_finder_p2t.vertexZ.push_back(vertexZ);
      m_finder_p2t.charge.push_back(charge);
      m_finder_p2t.pdg.push_back(pdg);
      m_finder_p2t.nTrueHits.push_back(nParticleHits);
      m_finder_p2t.assoc.push_back(tracks);
      m_finder_p2t.sharedHits.push_back(sh);
      m_finder_p2t.matchEfficiency.push_back(effs);
      m_finder_p2t.matchPurity.push_back(purs);
    }

    if (m_finder_p2t.tree)
      m_finder_p2t.tree->Fill();
    if (m_finder_t2p.tree)
      m_finder_t2p.tree->Fill();
  }

  // ---------- matching helper ----------
  int majorityParticleForTrack(const edm4hep::Track& trk,
                               const std::unordered_map<podio::ObjectID, int>& hitToParticle) const {
    std::unordered_map<int, int> counts;
    for (const auto& h : trk.getTrackerHits()) {
      const auto hk = h.getObjectID();
      auto it = hitToParticle.find(hk);
      if (it == hitToParticle.end())
        continue;
      counts[it->second] += 1;
    }
    if (counts.empty())
      return -1;

    int bestP = -1;
    int bestN = -1;
    for (const auto& kv : counts) {
      if (kv.second > bestN) {
        bestN = kv.second;
        bestP = kv.first;
      }
    }
    return bestP;
  }

  /**
   * @brief Compute the normalized residual (pull).
   *
   * Returns residual / sqrt(variance). If the variance is non-positive
   * or not finite, NaN is returned to avoid invalid values.
   */

  // ----------- pull calculation helper ----------
  static float safePull(float residual, float variance) {
    if (!std::isfinite(residual) || !std::isfinite(variance) || variance <= 0.f) {
      return std::numeric_limits<float>::quiet_NaN();
    }
    return residual / std::sqrt(variance);
  }
  // helper for err

  static float safeSqrt(float variance) {
    if (!std::isfinite(variance) || variance <= 0.f) {
      return std::numeric_limits<float>::quiet_NaN();
    }
    return std::sqrt(variance);
  }

  // ---------- fitter trees ----------
  template <typename PerfectMapT>
  void fillFitterTrees(int event, const edm4hep::MCParticleCollection& mcParts,
                       const edm4hep::TrackCollection& fittedTracks,
                       const std::unordered_map<podio::ObjectID, int>& hitToParticle,
                       const PerfectMapT& perfectAtIPByPid, bool doPerfect) const {

    m_fit_vs_mc.clear();
    m_fit_vs_perfect.clear();
    m_fit_vs_mc.event = event;
    m_fit_vs_perfect.event = event;

    const float NaN = std::numeric_limits<float>::quiet_NaN();

    int tIdx = 0;
    for (const auto& trk : fittedTracks) {
      auto stReco = TrackingValidationHelpers::getAtIPState(trk);
      if (!stReco) {
        ++tIdx;
        continue;
      }

      const int pid = majorityParticleForTrack(trk, hitToParticle);
      if (pid < 0 || pid >= (int)mcParts.size()) {
        ++tIdx;
        continue;
      }

      const auto& mc = mcParts[pid];

      // reco params
      TrackingValidationHelpers::HelixParams reco;
      reco.D0 = float(stReco->D0);
      reco.Z0 = float(stReco->Z0);
      reco.phi = float(stReco->phi);
      reco.omega = float(stReco->omega);
      reco.tanLambda = float(stReco->tanLambda);
      reco.pT = TrackingValidationHelpers::ptFromState(*stReco, m_Bz.value());
      reco.p = TrackingValidationHelpers::momentumFromState(*stReco, m_Bz.value());

      // ref from MC using the SAME convention as fitter (PCA + phi0 + ZPCA + omega=a*B/pT)
      const TrackingValidationHelpers::HelixParams refMC = TrackingValidationHelpers::truthFromMC_GenfitConvention(
          mc, m_Bz.value(), m_refX.value(), m_refY.value(), m_refZ.value());

      // Track-parameter uncertainties from the fitted-state covariance matrix

      const float varD0 = stReco->getCovMatrix(edm4hep::TrackParams::d0, edm4hep::TrackParams::d0);
      const float varPhi = stReco->getCovMatrix(edm4hep::TrackParams::phi, edm4hep::TrackParams::phi);
      const float varOmega = stReco->getCovMatrix(edm4hep::TrackParams::omega, edm4hep::TrackParams::omega);
      const float varZ0 = stReco->getCovMatrix(edm4hep::TrackParams::z0, edm4hep::TrackParams::z0);
      const float varTanL = stReco->getCovMatrix(edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);

      const float errD0 = safeSqrt(varD0);
      const float errPhi = safeSqrt(varPhi);
      const float errOmega = safeSqrt(varOmega);
      const float errZ0 = safeSqrt(varZ0);
      const float errTanL = safeSqrt(varTanL);

      // Residuals with respect to the MC truth helix parameters

      const float resD0 = reco.D0 - refMC.D0;
      const float resZ0 = reco.Z0 - refMC.Z0;
      const float resPhi = TrackingValidationHelpers::wrapDeltaPhi(reco.phi, refMC.phi);
      const float resOmega = reco.omega - refMC.omega;
      const float resTanL = reco.tanLambda - refMC.tanLambda;

      // --- vs MC ---
      m_fit_vs_mc.track_index.push_back(tIdx);
      m_fit_vs_mc.track_location.push_back(int(stReco->location));
      m_fit_vs_mc.resD0.push_back(resD0);
      m_fit_vs_mc.resZ0.push_back(resZ0);
      m_fit_vs_mc.resPhi.push_back(resPhi);
      m_fit_vs_mc.resOmega.push_back(resOmega);
      m_fit_vs_mc.resTanL.push_back(resTanL);
      m_fit_vs_mc.pullD0.push_back(safePull(resD0, varD0));
      m_fit_vs_mc.pullZ0.push_back(safePull(resZ0, varZ0));
      m_fit_vs_mc.pullPhi.push_back(safePull(resPhi, varPhi));
      m_fit_vs_mc.pullOmega.push_back(safePull(resOmega, varOmega));
      m_fit_vs_mc.pullTanL.push_back(safePull(resTanL, varTanL));
      m_fit_vs_mc.errD0.push_back(errD0);
      m_fit_vs_mc.errZ0.push_back(errZ0);
      m_fit_vs_mc.errPhi.push_back(errPhi);
      m_fit_vs_mc.errOmega.push_back(errOmega);
      m_fit_vs_mc.errTanL.push_back(errTanL);
      m_fit_vs_mc.p_reco.push_back(reco.p);
      m_fit_vs_mc.p_ref.push_back(refMC.p);
      m_fit_vs_mc.pT_reco.push_back(reco.pT);
      m_fit_vs_mc.pT_ref.push_back(refMC.pT);

      // --- vs perfect-fitted ---
      if (doPerfect) {
        auto it = perfectAtIPByPid.find(pid);
        if (it != perfectAtIPByPid.end()) {
          const auto& stPerf = it->second.st;

          TrackingValidationHelpers::HelixParams refP;
          refP.D0 = float(stPerf.D0);
          refP.Z0 = float(stPerf.Z0);
          refP.phi = float(stPerf.phi);
          refP.omega = float(stPerf.omega);
          refP.tanLambda = float(stPerf.tanLambda);
          refP.pT = TrackingValidationHelpers::ptFromState(stPerf, m_Bz.value());
          refP.p = TrackingValidationHelpers::momentumFromState(stPerf, m_Bz.value());

          const float resD0Perf = reco.D0 - refP.D0;
          const float resZ0Perf = reco.Z0 - refP.Z0;
          const float resPhiPerf = TrackingValidationHelpers::wrapDeltaPhi(reco.phi, refP.phi);
          const float resOmegaPerf = reco.omega - refP.omega;
          const float resTanLPerf = reco.tanLambda - refP.tanLambda;

          m_fit_vs_perfect.track_index.push_back(tIdx);
          m_fit_vs_perfect.track_location.push_back(int(stReco->location));
          m_fit_vs_perfect.resD0.push_back(resD0Perf);
          m_fit_vs_perfect.resZ0.push_back(resZ0Perf);
          m_fit_vs_perfect.resPhi.push_back(resPhiPerf);
          m_fit_vs_perfect.resOmega.push_back(resOmegaPerf);
          m_fit_vs_perfect.resTanL.push_back(resTanLPerf);
          m_fit_vs_perfect.pullD0.push_back(safePull(resD0Perf, varD0));
          m_fit_vs_perfect.pullZ0.push_back(safePull(resZ0Perf, varZ0));
          m_fit_vs_perfect.pullPhi.push_back(safePull(resPhiPerf, varPhi));
          m_fit_vs_perfect.pullOmega.push_back(safePull(resOmegaPerf, varOmega));
          m_fit_vs_perfect.pullTanL.push_back(safePull(resTanLPerf, varTanL));
          m_fit_vs_perfect.errD0.push_back(errD0);
          m_fit_vs_perfect.errZ0.push_back(errZ0);
          m_fit_vs_perfect.errPhi.push_back(errPhi);
          m_fit_vs_perfect.errOmega.push_back(errOmega);
          m_fit_vs_perfect.errTanL.push_back(errTanL);
          m_fit_vs_perfect.p_reco.push_back(reco.p);
          m_fit_vs_perfect.p_ref.push_back(refP.p);
          m_fit_vs_perfect.pT_reco.push_back(reco.pT);
          m_fit_vs_perfect.pT_ref.push_back(refP.pT);
        } else {
          m_fit_vs_perfect.track_index.push_back(tIdx);
          m_fit_vs_perfect.track_location.push_back(int(stReco->location));
          m_fit_vs_perfect.resD0.push_back(NaN);
          m_fit_vs_perfect.resZ0.push_back(NaN);
          m_fit_vs_perfect.resPhi.push_back(NaN);
          m_fit_vs_perfect.resOmega.push_back(NaN);
          m_fit_vs_perfect.resTanL.push_back(NaN);
          m_fit_vs_perfect.pullD0.push_back(NaN);
          m_fit_vs_perfect.pullZ0.push_back(NaN);
          m_fit_vs_perfect.pullPhi.push_back(NaN);
          m_fit_vs_perfect.pullOmega.push_back(NaN);
          m_fit_vs_perfect.pullTanL.push_back(NaN);
          m_fit_vs_perfect.errD0.push_back(errD0);
          m_fit_vs_perfect.errZ0.push_back(errZ0);
          m_fit_vs_perfect.errPhi.push_back(errPhi);
          m_fit_vs_perfect.errOmega.push_back(errOmega);
          m_fit_vs_perfect.errTanL.push_back(errTanL);
          m_fit_vs_perfect.p_reco.push_back(reco.p);
          m_fit_vs_perfect.p_ref.push_back(NaN);
          m_fit_vs_perfect.pT_reco.push_back(reco.pT);
          m_fit_vs_perfect.pT_ref.push_back(NaN);
        }
      }
      ++tIdx;
    }

    if (m_fit_vs_mc.tree)
      m_fit_vs_mc.tree->Fill();
    if (m_fit_vs_perfect.tree)
      m_fit_vs_perfect.tree->Fill();
  }

private:
  mutable int m_evt = 0;
  mutable bool m_warnedMissingPerfectInput = false;
  mutable bool m_warnedMissingFinderInput = false;
  mutable bool m_warnedMissingFittedInput = false;
  std::unique_ptr<TFile> m_outFile;

  mutable AssocTree m_finder_p2t;
  mutable AssocTree m_finder_t2p;
  mutable AssocTree m_perf_p2t;
  mutable AssocTree m_perf_t2p;

  mutable FitterTree m_fit_vs_mc;
  mutable FitterTree m_fit_vs_perfect;
};

DECLARE_COMPONENT(TrackingValidation)
