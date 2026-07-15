/*
 * Copyright (c) 2020-2026 Key4hep-Project.
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


/** dQdxAnalyser Gaudi Functional <br>
 *  The dQdxAnalyser functional is the translation of the Marlin dEdxAnalyser processor to Gaudi.
 *  It has been adapted to incorporate dN/dx in addition to dE/dx, together called dQ/dx, along with some more minor feature changes.
 *
 *  The dQdxAnalyser loops through tracks to collect their dQ/dx, fills various histograms and calculates dQ/dx resolution and separation power.
 *  Everything is stored in a root-file and various plots can be automatically produced.
 *
 *  The required input collections are Tracks, a TrackMCTruthLink and the corresponding MCParticles.
 *  If tracks attached to PFOs should be used, also a PFO collection is required.
 *  The number of particle species to be investigated and their properties are hard-coded in this header below:
 *  m_nPart is the number of species, together with arrays for the the PDG number and the names of the particles (for labelling the histograms).
 *  In this version: _nPart = 5; species: electrons, muons, pions, kaon, protons.
 *
 *  First, in operator(), the dQ/dx value of each track is filled in the corresponding histogram, resulting in Bethe-Bloch curves.
 *  The correct PDG is received from the MCParticle linked to the track.
 *  In finalize(), these histograms are then fitted with FitSlices, and the resulting mean and sigma values for each momentum bin are used to calculate resolution and separation power.
 *  The separation power is calculated for each combination of the m_nPart species.
 *  All histograms, including the auto-generated fit results, are stored in a root file.
 *  Finally, many of the histograms can be automatically plotted.
 *
 *  Several special histograms for calibration purposes are provided.
 *  The FiducialElectrons histograms collects electrons which are similar to electrons at the test beam. The resulting dQ/dx resolution in ResNorm_FiducialElectrons should be used to calibrate the smearing factor in the Compute_dEdxProcessor.
 *  The histogram NormLambdaFullAll_1 contains the mean dQ/dx over |lambda| incl. a poly3 fit, which is used in the AngularCorrection_dEdxProcessor for correction of the angular dependence.
 *  The histogram NormCosThFullAll_1 contains the mean dQ/dx over cos(theta), and shows the information of NormLambdaFullAll_1 over an alternative angular scale.
 *  The Bethe-Bloch curves in BBHist are fitted to provide the parameters of the reference curves for the dEdxPID in the LikelihoodPIDProcessor.
 *  The fit results are stored with the histograms and also printed on the console and in a text file at the end of the functional.
 *
 *  You can run this Analyser over any selection of tracks, e. g. only pion and kaon tracks.
 *  In that case, only the pion and kaon resolutions as well as only the pion-kaon separation would give sensible plots, but all other plots would be generated empty alongside.
 *
 *  Note that the hit energy histograms will only be filled, if the track hits are available, which is usually only the case for REC-files.
 *  The hit number histograms are unaffected by this.
 *
 *  @param _plotStuff - Set to true to switch on automatic image generation of many histograms in the current working directory.
 *    bool, default: false.
 *  @param _fileFormat - Select the file extension to be used for the automatically generated images.
 *    Selection depends on root TPad::Print() capabilities, typically: [.ps, .eps, .pdf, .svg, .tex, .gif, .xpm, .png, .jpg, .tiff, .cxx, .xml, .json, .root].
 *    See https://root.cern.ch/doc/master/classTPad.html#ad2fc5f449e4cb5480b9fb05fda3a8307 for details.
 *    string, default: .png
 *
 *  @param _usePFOTracks - Set to true to use tracks attached to charged PFO, instead of all tracks in the track collection. This should be used with events which are more busy than single particle files.
 *    bool, default: false.
 *  @param _useOneTrack - Set true if from every event only the first object in the selected collection (track or PFO) should be used. Usually, don't combine this with _usePFOTracks.
 *    bool, default: false.
 *
 *  You can make a track selection via the optional parameters, e. g. to reduce 'contamination' from mis-linked Tracks and MCParticles.
 *  @param _cutdEdx - Set true if particles with a very off dE/dx value (hard-coded) should be ignored.
 *    bool, default: false.
 *  @param _cutTrackPurity - Set true if particles should be ignored, where the MCParticle is connected to more than one track.
 *    bool, default: false.
 *  @param _cutD0 - Tracks with a d0 larger than the given value will be ignored. Set to 0 to accept all particles.
 *    double, default: 0.
 *  @param _cutZ0 - Tracks with a z0 larger than the given value will be ignored. Set to 0 to accept all particles.
 *    double, default: 0.
 *  @param _cutMomMin - Tracks with a momentum smaller than the given value will be ignored. Set to 0 to accept all particles.
 *    double, default: 0.
 *  @param _cutMomMax - Tracks with a momentum larger than the given value will be ignored. Set to 0 to accept all particles.
 *    double, default: 0.
 *
 *  For calibration with beam test results, so-called fiducial electrons are selected, which ought to be similar to beam-test conditions, and their properties stored in separate histograms.
 *  @param _fidMomMin - Fiducial minimum absolut momentum / GeV.
 *    double, default: 3.
 *  @param _fidMomMax - Fiducial maximum absolut momentum / GeV.
 *    double, default: 10.
 *  @param _fidLamMin - Fiducial minimum absolut lambda / deg.
 *    double, default: 1.
 *  @param _fidLamMax - Fiducial maximum absolut lambda / deg.
 *    double, default: 10.
 *  @param _fidNHitsMin - Fiducial minimum number of hits.
 *    int, default: 200.
 *
 *  For some of the histograms, the axis binning can be adapted, which can be useful in case of very high or low statistics.
 *  @param _nBinsX - Number of bins of the momentum axis in histograms.
 *    int, default: 100.
 *  @param _nBinsY - Number of bins of the dE/dx axis in histograms.
 *    int, default: 300.
 *  @param _minBinX - log_10(momentum/GeV) of the lower end of the momentum axis in histograms.
 *    double, default: -1 (= 0.1 GeV).
 *  @param _maxBinX - log_10(momentum/GeV) of the upper end of the momentum axis in histograms.
 *    double, default: 2 (= 100 GeV).
 *  @param _minBinY - Lower end of the dE/dx axis in GeV/mm in histograms.
 *    double, default: 0.
 *  @param _maxBinY - Upper end of the dE/dx axis in GeV/mm in histograms.
 *    double, default: 1e-6.
 *
 *  @author U. Einhaus, KIT
 *  @version $2.0$
 */

#include <k4FWCore/Consumer.h>

#include <edm4hep/RecDqdxCollection.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>
#include <edm4hep/TrackState.h>

#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DD4hep/DD4hepUnits.h>

#include <string>
#include <vector>
#include <array>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TObjArray.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

struct dQdxAnalyser final : k4FWCore::Consumer <void(
  const edm4hep::RecDqdxCollection&,
  const edm4hep::TrackCollection&,
  const std::vector<const edm4hep::ReconstructedParticleCollection*>&,
  const edm4hep::TrackMCParticleLinkCollection& )>
{

  dQdxAnalyser( const std::string&name, ISvcLocator* svcLoc)
  : Consumer (name, svcLoc, {
    KeyValues("InputdQdxCollection",         {"dQdx"}),
    KeyValues("InputTrackCollection",        {"Tracks"}),
    KeyValues("InputRecoParticleCollection", {"PandoraPFOs"}),
    KeyValues("InputTrackMCParticleLinkCollection", {"MarlinTrkTracksMCTruthLink"})
    })
  {};

  // ===== basic functions =====

  virtual StatusCode initialize() override;

  virtual void operator()
  (
    const edm4hep::RecDqdxCollection& dQdxs,
    const edm4hep::TrackCollection& tracks,
    const std::vector<const edm4hep::ReconstructedParticleCollection*>& PFOsVec,
    const edm4hep::TrackMCParticleLinkCollection& trackTruthLinks
  ) const override;

  virtual StatusCode finalize() override;

  // ====== helper functions =====

  virtual double Calculate_pValueSP(std::vector<double> v1, std::vector<double> v2);
  virtual void PlotTH1(TCanvas* can, TH1* hist);
  virtual void PlotTH2(TCanvas* can, TH2* hist);
  virtual void PlotTH1_FPn(TCanvas* can, TH1* hist, int n);

  // ===== Gaudi properties =====

  // output
  Gaudi::Property<std::string> m_outputFile{ this, "outputFile", "dQdxAnalyser_output.root", "Output root file with all histograms etc. Is overwritten if exists. Default: 'dQdxAnalyser_output.root'."};
  Gaudi::Property<std::string> m_outCalibFile{ this, "outCalibFile", "dQdxAnalyser_CalibrationOutput.txt", "Output txt file with fit results for dQ/dx calibration. Is overwritten if exists. Default: 'dQdxAnalyser_CalibrationOutput.txt'."};
  Gaudi::Property<bool> m_plotStuff{ this, "plotStuff", false, "Set true to automatically plot various histograms in the current working directory. Default: false."};
  Gaudi::Property<std::string> m_plotFolder{ this, "plotFolder", ".", "Folder in which the automatically plotted images should be stored (if chosen). Needs to exist. Default: '.' [current working directory]."};
  Gaudi::Property<std::string> m_plotFormat{ this, "plotFormat", ".png", "Format of the automatically plotted images (if chosen). Selection depends on root TPad::Print() capabilities, typically: [.ps, .eps, .pdf, .svg, .tex, .gif, .xpm, .png, .jpg, .tiff, .cxx, .xml, .json, .root]. Default: .png"};

  // geometry
  Gaudi::Property<int> m_subdetID{ this, "subdetID", 4, "Subdetector ID of the gaseous tracker (TPC/DC/etc.); typically 4 for TPC, 3 for DC; default: 4"};

  // binning
  Gaudi::Property<int> m_usedENdx { this, "usedENdx", 1, "Select if you analyse dE/dx (1) or dN/dx (2). Default: 1."};
  Gaudi::Property<int> m_nBinsMom { this, "nBinsMom",  100, "Number of bins of the momentum axis in histograms. Default:100"};
  Gaudi::Property<int> m_nBinsdQdx{ this, "nBinsdQdx", 300, "Number of bins of the dQ/dx axis in histograms. Default:300"};
  Gaudi::Property<int> m_nBinsLam { this, "nBinsLam",   18, "Number of bins of the lambda angle axis in histograms. Default:18"};
  Gaudi::Property<int> m_nBinsHits{ this, "nBinsHits", 110, "Number of bins of the numer-of-hits axis in histograms. Default:110"};
  Gaudi::Property<float> m_minBinMom { this, "minBinMom", -1.f, "log_10(momentum/GeV) of the lower end of the momentum axis in histograms. Default: -1 (= 0.1 GeV)."};
  Gaudi::Property<float> m_maxBinMom{ this, "maxBinMom",  2.f, "log_10(momentum/GeV) of the upper end of the momentum axis in histograms. Default:  2 (= 100 GeV)."};
  Gaudi::Property<float> m_minBindQdx{ this, "minBindQdx",  0.f, "Lower end of the dQ/dx axis in GeV/mm or 1/mm in histograms. Default: 0."};
  Gaudi::Property<float> m_maxBindQdx{ this, "maxBindQdx",1e-6f, "Upper end of the dQ/dx axis in GeV/mm or 1/mm in histograms. Default: 1e-6."};
  Gaudi::Property<float> m_minBinLam { this, "minBinLam",  0.f, "Lower end of the lambda angle axis in degree in histograms. Default: 0."};
  Gaudi::Property<float> m_maxBinLam { this, "maxBinLam", 90.f, "Upper end of the lambda angle axis in degree in histograms. Default: 90."};
  Gaudi::Property<int>   m_minBinHits{ this, "minBinHits",   0, "Lower end of the number-of-hits axis in histograms. Default: 0."};
  Gaudi::Property<int>   m_maxBinHits{ this, "maxBinHits", 220, "Upper end of the number-of-hits axis in histograms. Default: 220."};

  // track selection
  Gaudi::Property<bool> m_usePFOTracks{ this, "usePFOTracks", false, "Set true if tracks from (charged) PFOs should be used as data basis instead of all tracks from the TrackCollection."};
  Gaudi::Property<bool> m_useOneTrack{ this, "useOneTrack", false, "Set true if from every event only the first object in the selected collection (track or PFO) should be used."};
  Gaudi::Property<bool> m_usePFOTruthLink{ this, "usePFOTruthLink", false, "Set true if the LCRelation from PFO to MCTruth instead of Track to MCTruth should be used to determine MC PDG."};
  Gaudi::Property<bool> m_cutdEdx{ this, "cutdEdx", false, "Set true if particles with a very off dE/dx value should be ignored."};
  Gaudi::Property<bool> m_cutTrackPurity{ this, "cutTrackPurity", false, "Set true if particles should be ignored, where the MCParticle is connected to more than one track."};
  Gaudi::Property<bool> m_cutTrackPurityMom{ this, "cutTrackPurityMom", false, "Set true if particles should be ignored, where the MCParticle is connected to more than one track, but the track with the highest momentum is still used."};
  Gaudi::Property<float> m_cutD0{ this, "cutD0", 0.f, "Tracks with a d0/mm larger than the given value will be ignored. Set to 0 to accept all particles. Default: 0."};
  Gaudi::Property<float> m_cutZ0{ this, "cutZ0", 0.f, "Tracks with a z0/mm larger than the given value will be ignored. Set to 0 to accept all particles. Default: 0."};
  Gaudi::Property<float> m_cutMomMin{ this, "cutMomMin", 0.f, "Tracks with a momentum/GeV smaller than the given value will be ignored. Set to 0 to accept all particles. Default: 0."};
  Gaudi::Property<float> m_cutMomMax{ this, "cutMomMax", 0.f, "Tracks with a momentum/GeV larger than the given value will be ignored. Set to 0 to accept all particles. Default: 0."};
  Gaudi::Property<float> m_cutLamMin{ this, "cutLamMin", 0.f, "Tracks with an angle lambda/deg (relative to z=0) smaller than the given value will be ignored. Set to 0 to accept all particles. Default: 0."};
  Gaudi::Property<float> m_cutLamMax{ this, "cutLamMax", 0.f, "Tracks with an angle lambda/deg (relative to z=0) larger than the given value will be ignored. Set to 0 to accept all particles. Default: 0."};

  // fiducial electrons parameters
  Gaudi::Property<float> m_fidMomMin{ this, "fiducialMomMin",  3.f, "Fiducial minimum absolut momentum / GeV. Default: 3."};
  Gaudi::Property<float> m_fidMomMax{ this, "fiducialMomMax", 10.f, "Fiducial maximum absolut momentum / GeV. Default: 10."};
  Gaudi::Property<float> m_fidLamMin{ this, "fiducialLamMin",  1.f, "Fiducial minimum absolut lambda / deg. Default: 1."};
  Gaudi::Property<float> m_fidLamMax{ this, "fiducialLamMax", 10.f, "Fiducial maximum absolut lambda / deg. Default: 10."};
  Gaudi::Property<int> m_fidNHitsMin{ this, "fiducialNHitsMin", 200, "Fiducial minimum number of hits per track. Default: 200."};

  // ===== static properties =====

  std::string m_description = "dQdxAnalyser: Looks into track dQ/dx values, calculates resolution and separation power, puts everything in histograms and plots them.";

  mutable int m_nEvt = -1;
  mutable float m_Bz = 0; // from geometry in initialize()

  // additional binning, derived from input parameters
  mutable float m_minBinMom2=0, m_maxBindQdx2=0;
  mutable float *m_histBinsMom, *m_histBinsMom2;

  // number of particle species to be separated
  static const int m_nPart = 5;

  // PDG numbers, masses and histogram labels of particles, including "other" and "all"
  std::array<int, m_nPart+1> m_PDG = {11, 13, 211, 321, 2212, 10};
  std::array<double, m_nPart> m_Masses {0.000511, 0.10566, 0.13957, 0.49368, 0.93827};
  std::vector<std::string> m_Particles = {"Electrons","Muons","Pions","Kaons","Protons","Other","All"}; // for histogram names
  std::vector<std::string> m_particles = {"electrons","muons","pions","kaons","protons","all other particles","all particles"}; // for histogram descriptions
  std::vector<std::string> m_SP_part   = {"e","mu","pi","ka","p"}; // for separation power histogram names

  // vector array of all momentum, dQ/dx, dQ/dx error, tan lambda, number of hits and transverse momentum values for each track
  mutable std::vector< std::array<double,6> > m_dQdxVec[m_nPart];
  // vector array of all momentum and dQ/dx after angular correction
  //mutable std::vector< std::array<double,2> > m_angcorr_dQdxVec[m_nPart];

  // root file, histograms and fit results
  TFile* m_outFile{};

  // Bethe-Bloch histograms (including one for "other particles" and one for "all particles")
  TH2* m_BBHist[m_nPart+2];

  // Gaussian fit results: each TObjArray has 4 histograms (amplitude, mean, sigma and chi2)
  TObjArray* m_FitHist[m_nPart];
  TObjArray* m_FitLambda[m_nPart];
  TObjArray* m_FitNHit[m_nPart+1];
  TObjArray* m_FitLambdaFull[m_nPart+1];
  TObjArray* m_FitCosThFull[m_nPart+1];

  // storage for names of histograms to have a handle on fit results that are 'free floating' in memory, since TH3::FitSlicesZ has no direct transfer
  std::vector< std::string > m_NormBothHistNames[m_nPart];
  std::vector< std::string > m_NormBothHistNames_1GeV[m_nPart];

  // dQ/dx resolution over momentum histograms
  TH1* m_FitHist_Res[m_nPart];  // resolution from fit results
  TH1* m_FitHist_ResSum[m_nPart];  // resolution from fit results summed up above 1 GeV
  TH1* m_ResNormHist[m_nPart];  // resolution for all entries of a species above .3 GeV with normalisation
  TH1* m_ResNormHist_1GeV[m_nPart];  // resolution for all entries of a species above 1 GeV with normalisation
  TH1* m_ResElectrons{};  // resolution for all electron entries above .3 GeV, no normalisation
  TH1* m_ResNorm_FiducialElectrons{};  // resolution for all fiducial electrons with normalisation

  // resolution dependence on different parameters, fit histograms
  TH2* m_NormLambdaHist[m_nPart];  // relative dQdx vs. lambda of track, for tracks with p>1GeV, nHits>=200
  TH1* m_FitLambda_Res[m_nPart];  // resolution vs. lambda from fit results
  TH2* m_NormNHitHist[m_nPart+1];  // normalised dQdx vs. number of hits of track, for tracks with |lambda|<20deg
  TH1* m_FitNHit_Res[m_nPart+1];  // resolution vs. number of hits from fit results
  TH3* m_NormBothHist[m_nPart];  // relative dQdx vs. lambda of track and vs. number of hits
  TH2* m_FitBoth_Res[m_nPart];  // resolution vs. lambda and vs. number of hits from fit results
  TH3* m_NormBothHist_1GeV[m_nPart];  // relative dQdx vs. lambda of track and vs. number of hits, for tracks with p>1GeV
  TH2* m_FitBoth_Res_1GeV[m_nPart];  // resolution vs. lambda and vs. number of hits from fit results, for tracks with p>1GeV
  TH2* m_NormLambdaFullHist[m_nPart+1];  // relative dQdx vs. lambda of track, for tracks with p>1GeV
  TH1* m_FitLambdaFull_Res[m_nPart+1];  // resolution vs. lambda from fit results
  TH2* m_NormCosThFullHist[m_nPart+1];  // relative dQdx vs. lambda of track, for tracks with p>1GeV
  TH1* m_FitCosThFull_Res[m_nPart+1];  // resolution vs. lambda from fit results

  // separation power between all particle combinations
  TH1* m_SPHist[m_nPart*(m_nPart-1)/2];
  TH1* m_SPHist_pValue[m_nPart*(m_nPart-1)/2];
  TH1* m_SPHist_pValue_angcorr[m_nPart*(m_nPart-1)/2];

  // MC particle multiplicity according to track2mctruth link
  TH1* m_MCmult{};

  // contamination analyser histogram with given cuts
  TH1* m_CutAnaHist{};
  TH1* m_NTracksHist{};
  TH1* m_NTracksUsedHist{};
  TH1* m_NTracksMomHist{};
  TH1* m_TrackMomHist{};
  TH1* m_TrackSiblingsMomHist{};
  TH1* m_dEdxOutlierAbundance{};
  TH1* m_dEdxOutlierProportion{};

  // dQ/dx error display
  TH2* m_dQdxErrorCorrelation{};

  // hit charge spectrum
  TH1* m_HitEnergySpectrum[m_nPart+2];
  TH1* m_HitNumberSpectrum[m_nPart+2];
  TH1* m_HitEnergySpectrum_at1[m_nPart+2];
  TH1* m_HitEnergySpectrum_at25[m_nPart+2];
  TH2* m_HEHist[m_nPart+2];

  // Landau fit results: each TObjArray has 4 histograms (amplitude, MPV, sigma and chi2)
  TObjArray* m_HEFitHist[m_nPart];

};

DECLARE_COMPONENT(dQdxAnalyser)
