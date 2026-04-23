#include "k4FWCore/Consumer.h"


#include "Gaudi/Property.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/SmartIF.h"
#include "Gaudi/details/BranchWrapper.h"

#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/ClusterMCParticleLinkCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/TrackCollection.h"

#include "podio/ObjectID.h"
#include "podio/LinkNavigator.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <tuple>

namespace {
  template <typename V3>
  inline double thetaOf(const V3& p) {
    const double px = p.x, py = p.y, pz = p.z;
    return std::atan2(std::sqrt(px*px + py*py), pz);
  }

  template <typename V3>
  inline double phiOf(const V3& p) {
    return std::atan2((double)p.y, (double)p.x);
  }

  template <typename V3>
  inline double pTOf(const V3& p) {
    return std::hypot((double)p.x, (double)p.y);
  }

  inline double deltaPhiAbs(double a, double b) {
    double d = std::fabs(a - b);
    if (d > M_PI) d = 2.0*M_PI - d;
    return d;
  }
}

struct PFONTupleWriter final
    : k4FWCore::Consumer<void(const edm4hep::ReconstructedParticleCollection&, const edm4hep::MCParticleCollection&,
                              const edm4hep::RecoMCParticleLinkCollection&,
                              const edm4hep::RecoMCParticleLinkCollection&,
			      const edm4hep::TrackCollection&,
			      const edm4hep::TrackMCParticleLinkCollection&,
			      const edm4hep::SimTrackerHitCollection&,
			      const edm4hep::SimTrackerHitCollection&,
			      const edm4hep::SimTrackerHitCollection&,
			      const edm4hep::SimTrackerHitCollection&,
			      const edm4hep::SimTrackerHitCollection&,
			      const edm4hep::SimTrackerHitCollection&
			      )> {

  PFONTupleWriter(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(name, svcLoc,
                 {KeyValues("InputPFOs", {"PandoraPFOs"}), KeyValues("InputMCParticles", {"MCParticles"}),
                  KeyValues("InputRecoMCLink", {"RecoMCTruthLink"}),
                  KeyValues("InputMCRecoLink", {"MCTruthRecoLink"}),
		  KeyValues("InputTracks", {"SiTracks"}),
		  KeyValues("InputTrackMCLink", {"SiTracksMCTruthLink"}),
		  KeyValues("InputVXDBarrelSimHits", {"VertexBarrelCollection"}),
		  KeyValues("InputVXDEndcapSimHits", {"VertexEndcapCollection"}),
		  KeyValues("InputITBarrelSimHits", {"InnerTrackerBarrelCollection"}),
		  KeyValues("InputITEndcapSimHits", {"InnerTrackerEndcapCollection"}),
		  KeyValues("InputOTBarrelSimHits", {"OuterTrackerBarrelCollection"}),
		  KeyValues("InputOTEndcapSimHits", {"OuterTrackerEndcapCollection"})}){}

  Gaudi::Property<std::string> m_filePath{this, "FilePath", "PFONTuple.root", "long blabla"};
  std::unique_ptr<TFile> m_file = nullptr;
  std::unique_ptr<TTree> m_tree = nullptr;
  mutable std::map<std::string, TBranch*> m_branches;

  StatusCode initialize() override {
    m_file = std::make_unique<TFile>(m_filePath.value().c_str(), "RECREATE");
    m_file->cd();
    m_tree = std::make_unique<TTree>("PFONtuple", "PFONtuple");

    m_branches["foo"] = m_tree->Branch("foo", (double*)nullptr);
    m_branches["mc_pdg"] = m_tree->Branch("mc_pdg", (int*)nullptr);
    m_branches["mc_charge"] = m_tree->Branch("mc_charge", (double*)nullptr);
    m_branches["mc_mass"] = m_tree->Branch("mc_mass", (double*)nullptr);
    m_branches["mc_E"] = m_tree->Branch("mc_E", (double*)nullptr);
    m_branches["mc_Px"] = m_tree->Branch("mc_Px", (double*)nullptr);
    m_branches["mc_Py"] = m_tree->Branch("mc_Py", (double*)nullptr);
    m_branches["mc_Pz"] = m_tree->Branch("mc_Pz", (double*)nullptr);
    m_branches["mc_vx"] = m_tree->Branch("mc_vx", (double*)nullptr);
    m_branches["mc_vy"] = m_tree->Branch("mc_vy", (double*)nullptr);
    m_branches["mc_vz"] = m_tree->Branch("mc_vz", (double*)nullptr);
    m_branches["mc_t"] = m_tree->Branch("mc_t", (double*)nullptr);

    m_branches["nPFOs"] = m_tree->Branch("nPFOs", (int*)nullptr);
    m_branches["max_c_weight"] = m_tree->Branch("max_c_weight", (float*)nullptr);
    m_branches["max_t_weight"] = m_tree->Branch("max_t_weight", (float*)nullptr);
    m_branches["has_track"] = m_tree->Branch("has_track", (bool*)nullptr);
    m_branches["has_cluster"] = m_tree->Branch("has_cluster", (bool*)nullptr);
    m_branches["track_equals_cluster"] = m_tree->Branch("track_equals_cluster", (bool*)nullptr);
    m_branches["has_reco_pfo"] = m_tree->Branch("has_reco_pfo", (bool*)nullptr);
    m_branches["pfo_pdg"] = m_tree->Branch("pfo_pdg", (int*)nullptr);
    m_branches["pfo_E"] = m_tree->Branch("pfo_E", (double*)nullptr);
    m_branches["pfo_Px"] = m_tree->Branch("pfo_Px", (double*)nullptr);
    m_branches["pfo_Py"] = m_tree->Branch("pfo_Py", (double*)nullptr);
    m_branches["pfo_Pz"] = m_tree->Branch("pfo_Pz", (double*)nullptr);
    m_branches["pfo_vx"] = m_tree->Branch("pfo_vx", (double*)nullptr);
    m_branches["pfo_vy"] = m_tree->Branch("pfo_vy", (double*)nullptr);
    m_branches["pfo_vz"] = m_tree->Branch("pfo_vz", (double*)nullptr);
    m_branches["pfo_t"] = m_tree->Branch("pfo_t", (double*)nullptr);

    m_branches["dx"] = m_tree->Branch("dx", (double*)nullptr);
    m_branches["dy"] = m_tree->Branch("dy", (double*)nullptr);
    m_branches["dz"] = m_tree->Branch("dz", (double*)nullptr);
    m_branches["d_px"] = m_tree->Branch("d_px", (double*)nullptr);
    m_branches["d_py"] = m_tree->Branch("d_py", (double*)nullptr);
    m_branches["d_pz"] = m_tree->Branch("d_pz", (double*)nullptr);
    m_branches["dE"] = m_tree->Branch("dE", (double*)nullptr);
    m_branches["resE"] = m_tree->Branch("resE", (double*)nullptr);
    
    m_branches["okPhi"] = m_tree->Branch("okPhi", (bool*)nullptr);
    m_branches["okTheta"] = m_tree->Branch("okTheta", (bool*)nullptr);
    m_branches["okPt"] = m_tree->Branch("okPt", (bool*)nullptr);

    m_branches["mc_genStatus"] = m_tree->Branch("mc_genStatus", (int*)nullptr);


    
    return Consumer::initialize();
  }

  void fillTree(const std::map<std::string, void*>& fields) const {
    // TODO: grab lock here
    for (const auto& [name, ptr] : fields) {
      m_branches[name]->SetAddress(ptr);
    }
    m_tree->Fill();
  }

  /** Return track weight contribution encoded as trackwgt = (int(wgt)%10000)/1000. by
   * the encodeTrackAndClusterWeights function for the
   * ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   * 
   */
  // STOLEN FROM https://github.com/iLCSoft/MarlinUtil/blob/master/source/include/MarlinUtil.h
  static inline float getTrackWeight(float encodedWeight) {
      return float( int(encodedWeight) % 10000 ) / 1000.f;
  }


  /** Return cluster weight contribution encoded as clusterwgt = (int(wgt)/10000)/1000. by
   * the encodeTrackAndClusterWeights function for the
   * ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   */
  // STOLEN FROM https://github.com/iLCSoft/MarlinUtil/blob/master/source/include/MarlinUtil.h
  static inline float getClusterWeight(float encodedWeight) {
      return float( int(encodedWeight) / 10000 ) / 1000.f;
  }

  
  void operator()(const edm4hep::ReconstructedParticleCollection& PFOs, const edm4hep::MCParticleCollection& MCParticles,
                  const edm4hep::RecoMCParticleLinkCollection& RecoMCTruthLinks,
                  const edm4hep::RecoMCParticleLinkCollection& MCTruthRecoLinks,
		  const edm4hep::TrackCollection& Tracks,
		  const edm4hep::TrackMCParticleLinkCollection& trackMCLinks,
		  const edm4hep::SimTrackerHitCollection& vxdBarrelSimHits,
		  const edm4hep::SimTrackerHitCollection& vxdEndcapSimHits,
		  const edm4hep::SimTrackerHitCollection& itBarrelSimHits,
		  const edm4hep::SimTrackerHitCollection& itEndcapSimHits,
		  const edm4hep::SimTrackerHitCollection& otBarrelSimHits,
		  const edm4hep::SimTrackerHitCollection& otEndcapSimHits
		  ) const override {

    // // For podio::LinkCollections with disparate types just use getLinked
    std::map<std::string, void*> fields;

    // make navigators for both directions
    podio::LinkNavigator RecoMCTruthLinkNavigator(RecoMCTruthLinks);
    podio::LinkNavigator MCTruthRecoLinkNavigator(MCTruthRecoLinks);
    //podio::LinkNavigator CaloHitMCNavigator(CalohitMCTruthLinks);
    //podio::LinkNavigator ClusterMCNavigator(ClusterMCTruthLinks);
    podio::LinkNavigator trackMCNav(trackMCLinks);
    // Iterate over all mcparticles
    for (const auto& mcp : MCParticles) {
      int mc_genStatus = mcp.getGeneratorStatus();
      fields["mc_genStatus"] = &mc_genStatus;

      double foo = 42.;
      fields["foo"] = &foo;
      int mc_pdg = mcp.getPDG();
      fields["mc_pdg"] = &mc_pdg;

      double mc_charge = mcp.getCharge();
      fields["mc_charge"] = &mc_charge;
      double mc_mass = mcp.getMass();
      fields["mc_mass"] = &mc_mass;
      double mc_E = mcp.getEnergy();
      fields["mc_E"] = &mc_E;
      double mc_Px = mcp.getMomentum().x;
      fields["mc_Px"] = &mc_Px;
      double mc_Py = mcp.getMomentum().y;
      fields["mc_Py"] = &mc_Py;
      double mc_Pz = mcp.getMomentum().z;
      fields["mc_Pz"] = &mc_Pz;
      double mc_vx = mcp.getVertex().x;
      fields["mc_vx"] = &mc_vx;
      double mc_vy = mcp.getVertex().y;
      fields["mc_vy"] = &mc_vy;
      double mc_vz = mcp.getVertex().z;
      fields["mc_vz"] = &mc_vz;

      // Not functional, needs development to get actual timing info
      double mc_t = mcp.getTime();
      fields["mc_t"] = &mc_t;

      
      // TODO: figure out how many tracker hits, ecal hits, hcal hits, muon hits this particle caused and store in ntuple
      // somehow also figure out if all of those are really from this particle or if there was a non-stored secondary

      // get related PFOs, how many are related by track, how many by cluster, are they the same, is the strongest the same
      auto linkedPFOs = MCTruthRecoLinkNavigator.getLinked(mcp);
      int nPFOs = linkedPFOs.size();
      fields["nPFOs"] = &nPFOs;
      float max_c_weight = 0.;
      fields["max_c_weight"] = &max_c_weight;
      float max_t_weight = 0.;
      fields["max_t_weight"] = &max_t_weight;
      const edm4hep::ReconstructedParticle* max_c_pfo = nullptr;
      const edm4hep::ReconstructedParticle* max_t_pfo = nullptr;
      for (const auto& [pfo, weight] : linkedPFOs) {
        // need to parse weight
        float cw = getClusterWeight(weight);
        float tw = getTrackWeight(weight);
        if (cw > max_c_weight) {
          max_c_weight = cw;
          max_c_pfo = &pfo;
        }
        if (tw > max_t_weight) {
          max_t_weight = tw;
          max_t_pfo = &pfo;
        }
      }
      // if the strongest is not the same which one is stronger? (closest in energy?)
      // continue with the strongest PFO
      // what to do depends a bit on the initial mc particle
      // TODO: how pure is the pfo? how many foreign tracks are in the track, how many in the cluster, what are the fractions
      const edm4hep::ReconstructedParticle* reco_pfo = nullptr;
      bool has_track = max_t_weight > 0.5;
      fields["has_track"] = &has_track;
      bool has_cluster = max_c_weight > 0.5;
      fields["has_cluster"] = &has_cluster;
      bool track_equals_cluster = max_t_pfo == max_c_pfo;
      fields["track_equals_cluster"] = &track_equals_cluster;
      if (has_cluster) {
        reco_pfo = max_c_pfo;
      }
      if (mc_charge != 0. && has_track) {
        // if it exists use the track instead of the cluster for charged particles
        reco_pfo = max_t_pfo;
      }

      edm4hep::ReconstructedParticle dummy_pfo;
      bool has_reco_pfo = !!reco_pfo;
      fields["has_reco_pfo"] = &has_reco_pfo;
      if (!has_reco_pfo) {
        reco_pfo = &dummy_pfo;
      }

      bool okPhi = false;
      bool okTheta = false;
      bool okPt = false;
      
      const auto& rp = *reco_pfo;
      const auto& mc = mcp;

      std::cout << "MC PDG = " << mc_pdg
          << "  reco has? " << has_reco_pfo
          << "  reco objID = (" << reco_pfo->getObjectID().collectionID
          << "," << reco_pfo->getObjectID().index << ")"
          << "  reco PDG = " << reco_pfo->getPDG()
          << "  reco charge = " << reco_pfo->getCharge()
          << "  reco E = " << reco_pfo->getEnergy()
		<< "  reco px = " << reco_pfo->getMomentum().x
		<< "  reco py = " << reco_pfo->getMomentum().y
		<< "  reco pz = " << reco_pfo->getMomentum().z
          << std::endl;

      
      const double dphi =
	deltaPhiAbs(phiOf(rp.getMomentum()), phiOf(mc.getMomentum()));
      const double dtheta =
	std::fabs(thetaOf(rp.getMomentum()) - thetaOf(mc.getMomentum()));

      okPhi = (dphi < 0.002);
      okTheta = (dtheta < 0.001);

      if (std::fabs(mc.getCharge()) < 0.5) {
	okPt = true;
      }
      else {
	const double prt = pTOf(rp.getMomentum());
	const double pmt = pTOf(mc.getMomentum());
	okPt = (pmt > 0.0) &&
	  (std::fabs(prt / pmt - 1.0) < 0.05);
      }
      fields["okPhi"] = &okPhi;
      fields["okTheta"] = &okTheta;
      fields["okPt"] = &okPt;
	
      int pfo_pdg = reco_pfo->getPDG();
      fields["pfo_pdg"] = &pfo_pdg;

      

      // store 4 momentum in ntuple
      double pfo_E = reco_pfo->getEnergy();
      fields["pfo_E"] = &pfo_E;
      double pfo_Px = reco_pfo->getMomentum().x;
      fields["pfo_Px"] = &pfo_Px;
      double pfo_Py = reco_pfo->getMomentum().y;
      fields["pfo_Py"] = &pfo_Py;
      double pfo_Pz = reco_pfo->getMomentum().z;
      fields["pfo_Pz"] = &pfo_Pz;

      
      double pfo_t = -99; // Need timing information in Pandora
      fields["pfo_t"] = &pfo_t;
      
      double pfo_vx = reco_pfo->getReferencePoint().x;
      fields["pfo_vx"] = &pfo_vx;
      double pfo_vy = reco_pfo->getReferencePoint().y;
      fields["pfo_vy"] = &pfo_vy;
      double pfo_vz = reco_pfo->getReferencePoint().z;
      fields["pfo_vz"] = &pfo_vz;
      
      double dx = (pfo_vx - mc_vx);///mc_vx;
      fields["dx"] = &dx;
      double dy = (pfo_vy - mc_vy);///mc_vy;
      fields["dy"] = &dy;
      double dz = (pfo_vz - mc_vz);///mc_vz;
      fields["dz"] = &dz;

      double d_px = (pfo_Px - mc_Px);///mc_Px;
      fields["d_px"] = &d_px;
      double d_py = (pfo_Py - mc_Py);///mc_Py;
      fields["d_py"] = &d_py;
      double d_pz = (pfo_Pz - mc_Pz);///mc_Pz;
      fields["d_pz"] = &d_pz;

      double dE = (pfo_E - mc_E);///mc_E;
      fields["dE"] = &dE;

      double resE = (pfo_E/mc_E) - 1;
      fields["resE"] = &resE;

      fillTree(fields);

    }

  }

  StatusCode finalize() override {
    m_tree->Write();
    return Consumer::finalize();
  }
};

DECLARE_COMPONENT(PFONTupleWriter)
