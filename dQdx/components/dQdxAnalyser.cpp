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

#include <GaudiKernel/MsgStream.h>
#include <dQdxAnalyser.h>
#include <edm4hep/TrackerHit.h>
#include <edm4hep/utils/TrackUtils.h>
#include <filesystem>
#include <marlinutil/HelixClass.h>
#include <podio/LinkNavigator.h>

StatusCode dQdxAnalyser::initialize() //
{
  // make sure outputs are available
  m_outFile = new TFile(m_outputFile.value().c_str(), "RECREATE");
  if (!m_outFile || m_outFile->IsZombie()) // root file
  {
    error() << "Cannot open output root file: " << m_outputFile << endmsg;
    return StatusCode::FAILURE;
  }
  std::ofstream outfile;
  outfile.open(m_outCalibFile);
  if (!outfile.is_open()) // calibration file
  {
    error() << "Cannot open output calibration file: " << m_outCalibFile << endmsg;
    return StatusCode::FAILURE;
  }
  outfile.close();
  if (m_plotStuff && !std::filesystem::is_directory(m_plotFolder.value()))
    if (!std::filesystem::create_directory(m_plotFolder.value())) // plot folder
    {
      error() << "Cannot create output plot folder: " << m_plotFolder << endmsg;
      return StatusCode::FAILURE;
    }

  // prepare geometry info
  dd4hep::Detector& theDet = dd4hep::Detector::getInstance();
  double bfieldV[3];
  theDet.field().magneticField({0., 0., 0.}, bfieldV);
  m_Bz = bfieldV[2] / dd4hep::tesla;
  debug() << "B_z = " << m_Bz << endmsg;

  // prepare bins for histograms
  m_minBinMom2 = m_minBinMom - 1;
  m_maxBindQdx2 = m_maxBindQdx * 6;
  m_histBinsMom = new float[m_nBinsMom + 1];
  m_histBinsMom2 = new float[m_nBinsMom + 1];
  float* histBinsdQdx = new float[m_nBinsdQdx + 1];
  float* histBinsdQdx2 = new float[m_nBinsdQdx + 1];

  for (int i = 0; i < m_nBinsMom + 1; i++) // momentum axes
  {
    m_histBinsMom[i] = pow(10, m_minBinMom + (m_maxBinMom - m_minBinMom) * i / (float)(m_nBinsMom));
    m_histBinsMom2[i] = pow(10, m_minBinMom2 + (m_maxBinMom - m_minBinMom2) * i / (float)(m_nBinsMom));
  }
  for (int i = 0; i < m_nBinsdQdx + 1; i++) // dQ/dx axes
  {
    histBinsdQdx[i] = m_minBindQdx + (m_maxBindQdx - m_minBindQdx) * i / (float)(m_nBinsdQdx);
    histBinsdQdx2[i] = m_minBindQdx + (m_maxBindQdx2 - m_minBindQdx) * i / (float)(m_nBinsdQdx);
  }

  for (int i = 0; i < m_nPart + 2; i++) // create Bethe-Bloch histograms
  {
    std::string s1 = "BB" + m_Particles[i];
    std::string s2 = "Bethe-Bloch histogram for " + m_particles[i];
    m_BBHist[i] = new TH2D(s1.c_str(), s2.c_str(), m_nBinsMom, m_histBinsMom, m_nBinsdQdx, histBinsdQdx);
    m_BBHist[i]->SetXTitle("Momentum (GeV)");
    if (m_usedENdx == 1)
      m_BBHist[i]->SetYTitle("dE/dx (GeV/mm)");
    else if (m_usedENdx == 2)
      m_BBHist[i]->SetYTitle("dN/dx (1/mm)");
  }

  for (int i = 0; i < m_nPart; i++) // create space for fit results
  {
    m_FitHist[i] = new TObjArray();
    m_FitLambda[i] = new TObjArray();
    m_HEFitHist[i] = new TObjArray();
  }
  for (int i = 0; i < m_nPart + 1; i++) // some also for 'all particles'
  {
    m_FitNHit[i] = new TObjArray();
    m_FitLambdaFull[i] = new TObjArray();
    m_FitCosThFull[i] = new TObjArray();
  }

  int sp_c = 0;
  for (int i = 0; i < m_nPart; i++) // create resolution and separation power histograms
  {
    std::string s1 = "BB" + m_Particles[i] + "_FitRes";
    std::string s2 = "dQ/dx resolution over momentum for " + m_particles[i] + " from fit result";
    m_FitHist_Res[i] = new TH1D(s1.c_str(), s2.c_str(), m_nBinsMom, m_histBinsMom);
    m_FitHist_Res[i]->SetXTitle("momentum / GeV/c");
    m_FitHist_Res[i]->SetYTitle("dQ/dx resolution");

    std::string s3 = "BB" + m_Particles[i] + "_FitResSum";
    std::string s4 = "dQ/dx resolution distribution above 1 GeV for " + m_particles[i] + " from fit result";
    m_FitHist_ResSum[i] = new TH1D(s3.c_str(), s4.c_str(), 100, 0, .2);
    m_FitHist_ResSum[i]->SetXTitle("dQ/dx resolution");
    m_FitHist_ResSum[i]->SetYTitle("number of momentum bins with that resolution");

    std::string s5 = "Res" + m_Particles[i] + "_Norm";
    std::string s6 = "Relative dQ/dx resolution distribution for " + m_particles[i] + " above .3 GeV";
    m_ResNormHist[i] = new TH1D(s5.c_str(), s6.c_str(), 100, 0, 2);
    m_ResNormHist[i]->SetXTitle("relative dQ/dx");
    m_ResNormHist[i]->SetYTitle("rel. abundance");

    std::string s7 = "Res" + m_Particles[i] + "_Norm_1GeV";
    std::string s8 = "Relative dQ/dx resolution distribution for " + m_particles[i] + " above 1 GeV";
    m_ResNormHist_1GeV[i] = new TH1D(s7.c_str(), s8.c_str(), 100, 0, 2);
    m_ResNormHist_1GeV[i]->SetXTitle("relative dQ/dx");
    m_ResNormHist_1GeV[i]->SetYTitle("rel. abundance");

    // Resolution dependence on different parameters
    std::string s9 = "NormLambda" + m_Particles[i];
    std::string s10 = "Relative dQ/dx distribution over lambda with p>1 GeV and nHit>=200 for " + m_particles[i];
    m_NormLambdaHist[i] = new TH2D(s9.c_str(), s10.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, 100, 0, 2);
    m_NormLambdaHist[i]->SetXTitle("|lambda| / deg");
    m_NormLambdaHist[i]->SetYTitle("relative dQ/dx");

    std::string s11 = "NormLambda" + m_Particles[i] + "_FitRes";
    std::string s12 = "dQ/dx resolution for " + m_particles[i] + " over lambda from fit result";
    m_FitLambda_Res[i] = new TH1D(s11.c_str(), s12.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam);
    m_FitLambda_Res[i]->SetXTitle("|lambda| / deg");
    m_FitLambda_Res[i]->SetYTitle("dQ/dx resolution");

    std::string s13 = "NormNHit" + m_Particles[i];
    std::string s14 = "Relative dQ/dx distribution over number of hits with p>1GeV and |lambda|<20deg for " + //
                      m_particles[i];
    m_NormNHitHist[i] = new TH2D(s13.c_str(), s14.c_str(), m_nBinsHits, m_minBinHits + .5, m_maxBinHits + .5, //
                                 100, 0, 2);
    m_NormNHitHist[i]->SetXTitle("number of hits");
    m_NormNHitHist[i]->SetYTitle("relative dQ/dx");

    std::string s15 = "NormNHit" + m_Particles[i] + "_FitRes";
    std::string s16 = "dQ/dx resolution for " + m_particles[i] + " over number of hits from fit result";
    m_FitNHit_Res[i] = new TH1D(s15.c_str(), s16.c_str(), m_nBinsHits, m_minBinHits + .5, m_maxBinHits + .5);
    m_FitNHit_Res[i]->SetXTitle("number of hits");
    m_FitNHit_Res[i]->SetYTitle("dQ/dx resolution");

    std::string s20 = "NormBoth" + m_Particles[i];
    std::string s21 = "Relative dQ/dx distribution over lambda and number of hits for " + m_particles[i];
    m_NormBothHist[i] = new TH3D(s20.c_str(), s21.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits,
                                 m_minBinHits + .5, m_maxBinHits + .5, 100, 0, 2);
    m_NormBothHist[i]->SetXTitle("|lambda| / deg");
    m_NormBothHist[i]->SetYTitle("number of hits");
    m_NormBothHist[i]->SetZTitle("relative dQ/dx");

    std::string HistName = s20;
    m_NormBothHistNames[i].push_back(HistName.append("_0"));
    m_NormBothHistNames[i].push_back(HistName.replace(HistName.end() - 1, HistName.end(), "1"));
    m_NormBothHistNames[i].push_back(HistName.replace(HistName.end() - 1, HistName.end(), "2"));
    m_NormBothHistNames[i].push_back(HistName.replace(HistName.end() - 1, HistName.end(), "chi2"));

    std::string s22 = "NormBoth" + m_Particles[i] + "_FitRes";
    std::string s23 = "dQ/dx resolution for " + m_particles[i] + " over lambda and number of hits from fit result";
    m_FitBoth_Res[i] = new TH2D(s22.c_str(), s23.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits,
                                m_minBinHits + .5, m_maxBinHits + .5);
    m_FitBoth_Res[i]->SetXTitle("|lambda| / deg");
    m_FitBoth_Res[i]->SetYTitle("number of hits");
    m_FitBoth_Res[i]->SetZTitle("dQ/dx resolution");

    std::string s24 = "NormBoth" + m_Particles[i] + "_1GeV";
    std::string s25 = "Relative dQ/dx distribution over lambda and number of hits with p>1GeV for " + m_particles[i];
    m_NormBothHist_1GeV[i] = new TH3D(s24.c_str(), s25.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits,
                                      m_minBinHits + .5, m_maxBinHits + .5, 100, 0, 2);
    m_NormBothHist_1GeV[i]->SetXTitle("|lambda| / deg");
    m_NormBothHist_1GeV[i]->SetYTitle("number of hits");
    m_NormBothHist_1GeV[i]->SetZTitle("relative dQ/dx");

    HistName = s24;
    m_NormBothHistNames_1GeV[i].push_back(HistName.append("_amp"));
    m_NormBothHistNames_1GeV[i].push_back(HistName.replace(HistName.end() - 3, HistName.end(), "mean"));
    m_NormBothHistNames_1GeV[i].push_back(HistName.replace(HistName.end() - 4, HistName.end(), "sigma"));
    m_NormBothHistNames_1GeV[i].push_back(HistName.replace(HistName.end() - 5, HistName.end(), "chi2"));

    std::string s26 = "NormBoth" + m_Particles[i] + "_FitRes_1GeV";
    std::string s27 = "dQ/dx resolution for " + m_particles[i] + //
                      " with p>1GeV over lambda and number of hits from fit result";
    m_FitBoth_Res_1GeV[i] = new TH2D(s26.c_str(), s27.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits,
                                     m_minBinHits + .5, m_maxBinHits + .5);
    m_FitBoth_Res_1GeV[i]->SetXTitle("|lambda| / deg");
    m_FitBoth_Res_1GeV[i]->SetYTitle("number of hits");
    m_FitBoth_Res_1GeV[i]->SetZTitle("dQ/dx resolution");

    std::string s28 = "NormLambdaFull" + m_Particles[i];
    std::string s29 = "Relative dQ/dx distribution over lambda with p>1 GeV for " + m_particles[i];
    m_NormLambdaFullHist[i] = new TH2D(s28.c_str(), s29.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, 100, 0, 2);
    m_NormLambdaFullHist[i]->SetXTitle("|lambda| / deg");
    m_NormLambdaFullHist[i]->SetYTitle("relative dQ/dx");

    std::string s30 = "NormLambdaFull" + m_Particles[i] + "_FitRes";
    std::string s31 = "dQ/dx resolution for " + m_particles[i] + " over lambda from fit result";
    m_FitLambdaFull_Res[i] = new TH1D(s30.c_str(), s31.c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam);
    m_FitLambdaFull_Res[i]->SetXTitle("|lambda| / deg");
    m_FitLambdaFull_Res[i]->SetYTitle("dQ/dx resolution");

    std::string s32 = "NormCosThFull" + m_Particles[i];
    std::string s33 = "Relative dQ/dx distribution over cos(theta) with p>1 GeV for " + m_particles[i];
    m_NormCosThFullHist[i] = new TH2D(s32.c_str(), s33.c_str(), 20, 0, 1, 100, 0, 2);
    m_NormCosThFullHist[i]->SetXTitle("cos(theta)");
    m_NormCosThFullHist[i]->SetYTitle("relative dQ/dx");

    std::string s34 = "NormCosThFull" + m_Particles[i] + "_FitRes";
    std::string s35 = "dQ/dx resolution for " + m_particles[i] + " over cos(theta) from fit result";
    m_FitCosThFull_Res[i] = new TH1D(s34.c_str(), s35.c_str(), 20, 0, 1);
    m_FitCosThFull_Res[i]->SetXTitle("cos(theta)");
    m_FitCosThFull_Res[i]->SetYTitle("dQ/dx resolution");

    for (int j = i + 1; j < m_nPart; j++) // Separation Power
    {
      std::string s17 = "SP_" + m_SP_part[i] + "_" + m_SP_part[j];
      std::string s18 = "Separation power between " + m_particles[i] + " and " + m_particles[j];
      m_SPHist[sp_c] = new TH1D(s17.c_str(), s18.c_str(), m_nBinsMom, m_histBinsMom);
      m_SPHist[sp_c]->SetXTitle("momentum / GeV/c");
      m_SPHist[sp_c]->SetYTitle("separation power");

      std::string s36 = "SP_pValue_" + m_SP_part[i] + "_" + m_SP_part[j];
      std::string s37 = "Separation power by p-value between " + m_particles[i] + " and " + m_particles[j];
      m_SPHist_pValue[sp_c] = new TH1D(s36.c_str(), s37.c_str(), m_nBinsMom, m_histBinsMom);
      m_SPHist_pValue[sp_c]->SetXTitle("momentum / GeV/c");
      m_SPHist_pValue[sp_c]->SetYTitle("separation power");

      std::string s38 = "SP_pValue_angcorr_" + m_SP_part[i] + "_" + m_SP_part[j];
      std::string s39 = "Separation power by p-value with prior angular correction between " + m_particles[i] +
                        " and " + m_particles[j];
      m_SPHist_pValue_angcorr[sp_c] = new TH1D(s38.c_str(), s39.c_str(), m_nBinsMom, m_histBinsMom);
      m_SPHist_pValue_angcorr[sp_c]->SetXTitle("momentum / GeV/c");
      m_SPHist_pValue_angcorr[sp_c]->SetYTitle("separation power");

      sp_c++;
    }
  }

  // additional histograms for 'all particles'
  std::string s13 = "NormNHit" + m_Particles[m_nPart + 1];
  std::string s14 = "Relative dQ/dx distribution over number of hits with p>1GeV and |lambda|<20deg for " + //
                    m_particles[m_nPart + 1];
  m_NormNHitHist[m_nPart] = new TH2D(s13.c_str(), s14.c_str(), m_nBinsHits, m_minBinHits + .5, m_maxBinHits + .5, //
                                     100, 0, 2);
  m_NormNHitHist[m_nPart]->SetXTitle("number of hits");
  m_NormNHitHist[m_nPart]->SetYTitle("relative dQ/dx");

  std::string s15 = "NormNHit" + m_Particles[m_nPart + 1] + "_FitRes";
  std::string s16 = "dQ/dx resolution for " + m_particles[m_nPart + 1] + " over number of hits from fit result";
  m_FitNHit_Res[m_nPart] = new TH1D(s15.c_str(), s16.c_str(), m_nBinsHits, m_minBinHits + .5, m_maxBinHits + .5);
  m_FitNHit_Res[m_nPart]->SetXTitle("number of hits");
  m_FitNHit_Res[m_nPart]->SetYTitle("dQ/dx resolution");

  std::string s38 = "NormLambdaFull" + m_Particles[m_nPart + 1];
  std::string s39 = "Relative dQ/dx distribution over lambda with p>1 GeV for " + m_particles[m_nPart + 1];
  m_NormLambdaFullHist[m_nPart] = new TH2D(s38.c_str(), s39.c_str(), m_nBinsLam * 2, m_minBinLam, m_maxBinLam, //
                                           100, 0, 2);
  m_NormLambdaFullHist[m_nPart]->SetXTitle("|lambda| / deg");
  m_NormLambdaFullHist[m_nPart]->SetYTitle("relative dQ/dx");

  std::string s40 = "NormLambdaFull" + m_Particles[m_nPart + 1] + "_FitRes";
  std::string s41 = "dQ/dx resolution for " + m_particles[m_nPart + 1] + " over lambda from fit result";
  m_FitLambdaFull_Res[m_nPart] = new TH1D(s40.c_str(), s41.c_str(), m_nBinsLam * 2, m_minBinLam, m_maxBinLam);
  m_FitLambdaFull_Res[m_nPart]->SetXTitle("|lambda| / deg");
  m_FitLambdaFull_Res[m_nPart]->SetYTitle("dQ/dx resolution");

  std::string s42 = "NormCosThFull" + m_Particles[m_nPart + 1];
  std::string s43 = "Relative dQ/dx distribution over cos(theta) with p>1 GeV for " + m_particles[m_nPart + 1];
  m_NormCosThFullHist[m_nPart] = new TH2D(s42.c_str(), s43.c_str(), 20, 0, 1, 100, 0, 2);
  m_NormCosThFullHist[m_nPart]->SetXTitle("cos(theta)");
  m_NormCosThFullHist[m_nPart]->SetYTitle("relative dQ/dx");

  std::string s44 = "NormCosThFull" + m_Particles[m_nPart + 1] + "_FitRes";
  std::string s45 = "dQ/dx resolution for " + m_particles[m_nPart + 1] + " over cos(theta) from fit result";
  m_FitCosThFull_Res[m_nPart] = new TH1D(s44.c_str(), s45.c_str(), 20, 0, 1);
  m_FitCosThFull_Res[m_nPart]->SetXTitle("cos(theta)");
  m_FitCosThFull_Res[m_nPart]->SetYTitle("dQ/dx resolution");

  // fiducial electrons
  m_ResNorm_FiducialElectrons = new TH1D("ResNorm_FiducialElectrons", //
                                         "Relative dQ/dx resolution distribution for fiducial electrons", 100, 0, 2);

  // some other histograms for cross-checks etc.
  m_ResElectrons = new TH1D("ResElectrons",
                            "dQ/dx resolution distribution above .3 GeV for electrons from all entries, not normalised",
                            100, m_minBindQdx, m_maxBindQdx / 2.);
  m_MCmult = new TH1D("MCmult", "MC particle multiplicity of linked reco objects (tracks or PFOs)", 11, -.5, 10.5);

  m_CutAnaHist = new TH1D("CutAnaHist", "Contamination analyser histogram with given cuts", 10, -.5, 9.5);
  m_CutAnaHist->GetXaxis()->SetBinLabel(1, "Track exists");
  m_CutAnaHist->GetXaxis()->SetBinLabel(2, "MCPar and dQ/dx exist");
  m_CutAnaHist->GetXaxis()->SetBinLabel(3, "Track purity (1 track) passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(4, "Track state exists");
  m_CutAnaHist->GetXaxis()->SetBinLabel(5, "Momentum cuts passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(6, "Angular cuts passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(7, "d0/z0 cuts passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(8, "Track purity (highest mom) passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(9, "dE/dx outliers passed");

  m_NTracksHist = new TH1D("NTracksHist", "Number of tracks of each event", 100, -.5, 99.5);
  m_NTracksUsedHist = new TH1D("NTracksUsedHist", "Number of tracks of each event used after cuts", 100, -.5, 99.5);
  m_NTracksMomHist = new TH1D("NTracksMomHist", "Number of tracks of each event used after cuts with p > 100 MeV", //
                              100, -.5, 99.5);
  m_TrackMomHist = new TH1D("TrackMomHist", "Momentum of tracks after cuts", m_nBinsMom, m_histBinsMom2);
  m_TrackMomHist->SetXTitle("momentum (GeV)");
  m_TrackMomHist->SetYTitle("abundance");
  m_TrackSiblingsMomHist = new TH1D("TrackSiblingsMomHist", "Momentum of tracks after cuts in events with nTracks > 1",
                                    m_nBinsMom, m_histBinsMom2);
  m_TrackSiblingsMomHist->SetXTitle("momentum (GeV)");
  m_TrackSiblingsMomHist->SetYTitle("abundance");

  m_dEdxOutlierAbundance = new TH1D("dEdxOutlierAbundance", // dE/dx only
                                    "Number of tracks with dE/dx values that don't match the MC PDG", 5, -.5, 4.5);
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(1, "MCPion - dE/dx too large");
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(2, "MCKaon - dE/dx too large");
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(3, "MCKaon - dE/dx too small");
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(4, "MCProton - dE/dx too small");
  m_dEdxOutlierAbundance->SetYTitle("abundance");

  m_dEdxOutlierProportion = new TH1D("dEdxOutlierProportion", // dE/dx only
                                     "Proportion of tracks with dE/dx values that don't match the MC PDG", 5, -.5, 4.5);
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(1, "MCPion - dE/dx too large");
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(2, "MCKaon - dE/dx too large");
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(3, "MCKaon - dE/dx too small");
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(4, "MCProton - dE/dx too small");
  m_dEdxOutlierProportion->SetYTitle("proportion of MCSpecies");

  m_dQdxErrorCorrelation = new TH2D("dQdxErrorCorrelation", "Correlation plot of dQ/dx error vs. dQ/dx as cross check",
                                    100, m_minBindQdx, m_maxBindQdx, 100, m_minBindQdx * .1, m_maxBindQdx * .1);

  for (int i = 0; i < m_nPart + 2; i++) // hit number and energy spectra
  {
    std::string s1 = "HitNum_" + m_Particles[i];
    std::string s2 = "Hit number spectrum for " + m_particles[i];
    m_HitNumberSpectrum[i] = new TH1D(s1.c_str(), s2.c_str(), m_maxBinHits, 0, m_maxBinHits);
    m_HitNumberSpectrum[i]->SetXTitle("Hit Number");
    m_HitNumberSpectrum[i]->SetYTitle("Abundance");

    std::string s3 = "HitEng_" + m_Particles[i];
    std::string s4 = "Hit energy spectrum for " + m_particles[i];
    m_HitEnergySpectrum[i] = new TH1D(s3.c_str(), s4.c_str(), 100, 0, 1e-5);
    m_HitEnergySpectrum[i]->SetXTitle("Hit Energy (GeV)");
    m_HitEnergySpectrum[i]->SetYTitle("Abundance");

    std::string s5 = "HitEng_at1_" + m_Particles[i];
    std::string s6 = "Hit energy spectrum around 1 GeV for " + m_particles[i];
    m_HitEnergySpectrum_at1[i] = new TH1D(s5.c_str(), s6.c_str(), 100, 0, 1e-5);
    m_HitEnergySpectrum_at1[i]->SetXTitle("Hit Energy (GeV)");
    m_HitEnergySpectrum_at1[i]->SetYTitle("Abundance");

    std::string s7 = "HitEng_at25_" + m_Particles[i];
    std::string s8 = "Hit energy spectrum around 2.5 GeV for " + m_particles[i];
    m_HitEnergySpectrum_at25[i] = new TH1D(s7.c_str(), s8.c_str(), 100, 0, 1e-5);
    m_HitEnergySpectrum_at25[i]->SetXTitle("Hit Energy (GeV)");
    m_HitEnergySpectrum_at25[i]->SetYTitle("Abundance");

    std::string s9 = "HEHist_" + m_Particles[i];
    std::string s10 = "Hit energy curves over momentum for " + m_particles[i];
    m_HEHist[i] = new TH2D(s9.c_str(), s10.c_str(), m_nBinsMom, m_histBinsMom, m_nBinsdQdx, histBinsdQdx2);
    m_HEHist[i]->SetXTitle("Momentum (GeV)");
    m_HEHist[i]->SetYTitle("Hit Energy (GeV)");
  }

  delete[] histBinsdQdx;
  delete[] histBinsdQdx2;

  return StatusCode::SUCCESS;
}

template <typename C, typename E>
bool contains(C&& cont, E ele) // helper function: if container contains element
{
  return std::find(std::begin(cont), std::end(cont), ele) != std::end(cont);
}
template <typename C, typename E>
int index_of(C&& cont, E ele) // helper function: index of element in container
{
  int idx = 0;
  for (auto c : cont) {
    if (c == ele)
      return idx;
    idx++;
  }
  return -1;
}

void dQdxAnalyser::operator()(const edm4hep::RecDqdxCollection& dQdxs, //
                              const edm4hep::TrackCollection& tracks,
                              const std::vector<const edm4hep::ReconstructedParticleCollection*>& pfosVec,
                              const edm4hep::TrackMCParticleLinkCollection& trackTruthLinks) const {
  m_nEvt++;
  debug() << "dQdxAnalyser operates on event " << m_nEvt << endmsg;

  int n_track = tracks.size();
  m_NTracksHist->Fill(n_track);

  int n_pfo = 0;
  const edm4hep::ReconstructedParticleCollection* pfos{};
  if (m_usePFOTracks) // check if tracks from PFOs should be used
  {
    if (pfosVec.size() < 1) // and if the PFO collection exists
    {
      error() << "PFO collection not found - skipping event " << m_nEvt << endmsg;
      return;
    } else {
      pfos = pfosVec[0];
      n_pfo = pfos->size();
    }
  }
  auto trackTruthNav = podio::LinkNavigator(trackTruthLinks);

  // check how many tracks to loop over, then loop
  int n_track_pfo = m_usePFOTracks ? n_pfo : n_track;
  int n_track_mod = m_useOneTrack ? std::min(n_track_pfo, 1) : n_track_pfo;

  int n_track_used = 0;
  int n_track_mom = 0;

  for (int i = 0; i < n_track_mod; ++i) // loop over tracks
  {
    int pdg = 0;
    double maxwei = 0;
    double maxtrackmom = 0;
    int icut = 0;

    edm4hep::Track track{};
    edm4hep::ReconstructedParticle pfo{};

    if (m_usePFOTracks) // get tracks from PFOs or directly
    {
      pfo = pfos->at(i);
      if (abs(pfo.getCharge()) != 1 || pfo.getTracks().size() != 1)
        continue;
      track = pfo.getTracks()[0];
    } else
      track = tracks[i];

    edm4hep::MCParticle* mcpar = NULL;
    edm4hep::utils::TrackPIDHandler dQdxHan{};
    dQdxHan.addColl(dQdxs);

    auto mcparWeiVec = trackTruthNav.getLinked(track);
    debug() << "mcparWeiVec has length of " << mcparWeiVec.size() << endmsg;
    for (unsigned int j = 0; j < mcparWeiVec.size(); ++j) // get MCPar with highest weight
    {
      float wei = mcparWeiVec[j].weight;
      if (wei > .5 && wei > maxwei) // require a weight of at least 0.5 in any case
      {
        mcpar = &mcparWeiVec[j].o;
        int mcpdg = abs(mcpar->getPDG());
        pdg = contains(m_PDG, mcpdg) ? mcpdg : 10;
        maxwei = wei;
      }
      debug() << " entry " << j << " with PDG " << pdg << " and weight of " << wei << endmsg;
    }

    m_MCmult->Fill(mcparWeiVec.size());
    m_CutAnaHist->Fill(icut++);

    // get dQ/dx value and error
    auto dQdxVec = dQdxHan.getDqdxValues(track);
    float dQdx = dQdxVec.size() == 1 ? dQdxVec[0].getDQdx().value : 0;
    float dQdxErr = dQdxVec.size() == 1 ? dQdxVec[0].getDQdx().error : 0;
    debug() << "dQdx has value of " << dQdx << " +/- " << dQdxErr << endmsg;

    // make sure, there is a mcpar linked and a finite dQ/dx
    if (pdg == 0 || dQdx <= 0 || !mcpar)
      continue;
    m_CutAnaHist->Fill(icut++);

    int n_rel_tracks = trackTruthNav.getLinked(*mcpar).size();

    // if more other tracks are connected to the same MCParticle, skip this track
    if (m_cutTrackPurity && n_rel_tracks > 1)
      continue;
    m_CutAnaHist->Fill(icut++);

    // get track state, get track properties
    auto trackStates = track.getTrackStates();
    auto it = std::ranges::find(trackStates, 1, &edm4hep::TrackState::location);
    auto opt = it != trackStates.end() ? std::optional{*it} : std::nullopt;

    float mom = 0, momT = 0, lambda = 0, tanL = 0, d0 = 0, z0 = 0;
    if (opt) // use track state to get momentum
    {
      auto state = *opt;
      std::tie(mom, momT) = getMomFromTrackState(state);
      tanL = state.tanLambda;
      lambda = fabs(atan(tanL)) * 180 / M_PI;
      d0 = state.D0;
      z0 = state.Z0;
    } else
      continue;
    m_CutAnaHist->Fill(icut++);

    // momentum cuts
    if (m_cutMomMin && mom < m_cutMomMin)
      continue;
    if (m_cutMomMax && mom > m_cutMomMax)
      continue;
    m_CutAnaHist->Fill(icut++);

    // angular cuts
    if (m_cutLamMin && lambda < m_cutLamMin)
      continue;
    if (m_cutLamMax && lambda > m_cutLamMax)
      continue;
    m_CutAnaHist->Fill(icut++);

    // d0/z0 cuts
    if (m_cutD0 && fabs(d0) > m_cutD0)
      continue;
    if (m_cutZ0 && fabs(z0) > m_cutZ0)
      continue;
    m_CutAnaHist->Fill(icut++);

    // if additional tracks are connected to the same MCParticle and this is not the track with the highest ...
    if (m_cutTrackPurityMom && n_rel_tracks > 1) // ... momentum among them, skip this track
    {
      auto trackWeiVec = trackTruthNav.getLinked(*mcpar);
      for (int j = 0; j < n_rel_tracks; ++j) // get track states of additional tracks
      {
        auto tempTrack = trackWeiVec[j].o;
        auto temptrackStates = tempTrack.getTrackStates();
        auto tempit = std::ranges::find(temptrackStates, 1, &edm4hep::TrackState::location);
        auto tempopt = tempit != temptrackStates.end() ? std::optional{*tempit} : std::nullopt;
        if (tempopt) // get helix to get momentum
        {
          auto state = *opt;
          auto [tempmom, tempmomT] = getMomFromTrackState(state);
          if (tempmom > maxtrackmom)
            maxtrackmom = tempmom;
        }
      }
      if (mom < maxtrackmom)
        continue;
    }
    m_CutAnaHist->Fill(icut++);

    // register outlier dE/dx values
    if (m_usedENdx == 1) // only if dE/dx is used
    {
      bool cutdEdx = false;
      if (pdg == 211 && dQdx > (pow(log10(mom), 2) + .2) * 1e-6) // kaons and protons linked to a pion
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 211, dE/dx = " << dQdx
                << ", upper limit = " << (pow(log10(mom), 2) + .2) * 1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(0);
        cutdEdx = true;
      }
      if (pdg == 321 && dQdx > (pow(log10(mom) * 1.3, 2) + .2) * 1e-6) // pions linked to a kaon
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 321, dE/dx = " << dQdx
                << ", upper limit = " << (pow(log10(mom) * 1.3, 2) + .2) * 1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(1);
        cutdEdx = true;
      }
      if (pdg == 321 && dQdx < (.28 - mom / 3.) * 1e-6) // protons linked to a kaon
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 321, dE/dx = " << dQdx
                << ", lower limit = " << (.28 - mom / 3.) * 1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(2);
        cutdEdx = true;
      }
      if (pdg == 2212 && dQdx < (.3 - mom / 5.) * 1e-6) // pions and kaons linked to a protons
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 2212, dE/dx = " << dQdx
                << ", lower limit = " << (.3 - mom / 5.) * 1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(3);
        cutdEdx = true;
      }

      // if the track dE/dx value is obviously not matching the MCParticle PDG, skip this track
      if (m_cutdEdx && cutdEdx)
        continue;
    }
    m_CutAnaHist->Fill(icut++);

    // All cuts passed!
    n_track_used++;
    if (mom > .1)
      n_track_mom++;
    m_TrackMomHist->Fill(mom);
    if (n_track > 1)
      m_TrackSiblingsMomHist->Fill(mom);
    int nHits = track.getSubdetectorHitNumbers(m_subdetID * 2 - 2);
    debug() << "gaseous tracker hit number: " << nHits << endmsg;

    // fill Bethe-Bloch histograms and data vector used in ::finalize()
    int pdg_idx = index_of(m_PDG, pdg);
    if (pdg_idx >= 0) // particle is one of the target PDGs or 10
    {
      m_BBHist[pdg_idx]->Fill(mom, dQdx);
      if (pdg_idx < m_nPart && mom > .3 && mom < 100)
        m_dQdxVec[pdg_idx].emplace_back(std::array<double, 6>{{mom, dQdx, dQdxErr, tanL, (double)nHits, momT}});
    }
    m_BBHist[m_nPart + 1]->Fill(mom, dQdx);

    if (pdg == 11 && mom > .3)
      m_ResElectrons->Fill(dQdx);
    m_dQdxErrorCorrelation->Fill(dQdx, dQdxErr);

    int trkHitsSize = track.trackerHits_size();
    for (int h = 0; h < trkHitsSize; ++h) // fill hit energy histograms
    {
      double EDep = track.getTrackerHits(h).getEDep();
      if (pdg_idx >= 0)
        m_HitEnergySpectrum[pdg_idx]->Fill(EDep);
      m_HitEnergySpectrum[m_nPart + 1]->Fill(EDep);

      if (.9 < mom && mom < 1.15) // p around 1 GeV
      {
        if (pdg_idx >= 0)
          m_HitEnergySpectrum_at1[pdg_idx]->Fill(EDep);
        m_HitEnergySpectrum_at1[m_nPart + 1]->Fill(EDep);
      }

      if (2 < mom && mom < 3) // p around 2.5 GeV
      {
        if (pdg_idx >= 0)
          m_HitEnergySpectrum_at25[pdg_idx]->Fill(EDep);
        m_HitEnergySpectrum_at25[m_nPart + 1]->Fill(EDep);
      }

      // 2D correlation between p and EDep
      if (pdg_idx >= 0)
        m_HEHist[pdg_idx]->Fill(mom, EDep);
      m_HEHist[m_nPart + 1]->Fill(mom, EDep);
    }

    if (pdg_idx >= 0)
      m_HitNumberSpectrum[pdg_idx]->Fill(nHits);
    m_HitNumberSpectrum[m_nPart + 1]->Fill(nHits);

  } // -loop over tracks

  m_NTracksUsedHist->Fill(n_track_used);
  m_NTracksMomHist->Fill(n_track_mom);
}

double BB_Fit(double* x, double* p) // Bethe-Bloch fit function
{
  double bg = x[0] / p[0];
  double bb2 = bg * bg / (1 + bg * bg);
  return 1.35e-7 * (.5 * p[1] * TMath::Log(p[2] * pow(bg, 2) * pow(bg, 2)) - p[3] * bb2 - p[4] * bg / 2.) / bb2;
}

StatusCode dQdxAnalyser::finalize() //
{
  info() << "finalising..." << endmsg;

  // set up fitting functions
  TF1* gausfit = new TF1("gausfit", "gaus(0)", 0, m_maxBindQdx);
  gausfit->SetParameters(1, m_maxBindQdx / 3., m_maxBindQdx / 10.);
  gausfit->SetParLimits(0, 0, 1e6);
  gausfit->SetParLimits(1, 0, m_maxBindQdx);
  gausfit->SetParLimits(2, 0, m_maxBindQdx / 2.);

  TF1* gausfitrel = new TF1("gausfitrel", "gaus(0)", 0, 2);
  gausfitrel->SetParameters(m_nEvt / 100., 1, .2);

  TF2* powla2 = new TF2("powla2", "0.047*([0]/cos(x/180.*TMath::Pi()))**[1]*(y/[2])**[3]", m_minBinLam, m_maxBinLam,
                        m_minBinHits, m_maxBinHits);
  powla2->SetParameters(1, -0.34, 220, -.47);

  TF1* landaufit = new TF1("landaufit", "landau(0)", 0, m_maxBindQdx2);

  TF1* poly2fit = new TF1("poly2fit", "pol2", 0, 90);
  TF1* poly3fit = new TF1("poly3fit", "pol3", 0, 90);

  TF1* lamfit = new TF1("lamfit", "[0]*((cos(x/180.*TMath::Pi()))**[1])", m_minBinLam, m_maxBinLam);
  lamfit->SetParameters(.05, -.34);
  TF1* nhitfit = new TF1("nHitfit", "[0]*((x/220.)**[1])", m_minBinHits, m_maxBinHits);
  nhitfit->SetParameters(.05, -.47);

  TF1* angcorr = new TF1("angcorr", "1./(1+[0]*x*x)", 0, 1); // x = cos(theta)
  angcorr->SetParameter(0, -.1);
  TF1* nhitcorr = new TF1("nhitcorr", "1+exp(-x/[0])", 1e-6, 1e6); // x = nHits
  nhitcorr->SetParameter(0, 1e1);

  TF1* BBfit = new TF1("BBfit", BB_Fit, pow(10, m_minBinMom), pow(10, m_maxBinMom), 5);
  BBfit->SetParameters(m_Masses[0], -1.28883368e-02, 2.72959919e+02, -1.74534200e+00, -9.84887586e-07);
  std::array<TFitResultPtr, m_nPart> ptr_BBfit;

  // fitting and extraction of results
  int sp_c = 0;
  for (int i = 0; i < m_nPart; i++) //
  {
    FitAndGetResolution(m_BBHist[i], gausfit, m_FitHist[i], m_FitHist_Res[i]);
    for (int k = 1; k < m_nBinsMom + 1; k++) // 0 is underflow bin, nBinsx+1 is overflow bin
      if (m_histBinsMom[k] >= 1)
        m_FitHist_ResSum[i]->Fill(m_FitHist_Res[i]->GetBinContent(k));
    m_FitHist_ResSum[i]->Fit("gaus", "Q");

    TH1* FitHist1_mean = (TH1D*)m_FitHist[i]->At(1);
    BBfit->FixParameter(0, m_Masses[i]);
    ptr_BBfit[i] = FitHist1_mean->Fit(BBfit, "+S");
    debug() << "BBFit done for particle species " << i << endmsg;

    // normalised resolution, i.e. corrected for species and momentum
    for (unsigned int v = 0; v < m_dQdxVec[i].size(); v++) // loop over all recorded particles
    {
      // m_dQdxVec[i][v] = mom, dQdx, dQdxErr, tanL, nHits, momT
      double mom = m_dQdxVec[i][v][0];
      double dQdx = m_dQdxVec[i][v][1];
      int bin = FitHist1_mean->FindBin(mom);
      double binMean = FitHist1_mean->GetBinContent(bin);
      double binError = FitHist1_mean->GetBinError(bin);
      double normMean = dQdx / binMean;
      double absError = sqrt(pow(m_dQdxVec[i][v][2] / binMean, 2) + pow(dQdx / binMean / binMean * binError, 2));
      double relError = absError / normMean;
      double tanL = m_dQdxVec[i][v][3];
      double lambda = fabs(atan(tanL)) * 180 / M_PI;
      double costh = sqrt(tanL * tanL / (1 + tanL * tanL));
      int nHits = m_dQdxVec[i][v][4];
      double momT = m_dQdxVec[i][v][5];

      if (binMean > 0) // correlation of normalised mean value with nHits or lambda or both
      {
        m_ResNormHist[i]->Fill(normMean);
        if (momT > 1)
          m_ResNormHist_1GeV[i]->Fill(normMean);

        if (momT > 1 && nHits >= 200)
          m_NormLambdaHist[i]->Fill(lambda, normMean);
        if (lambda < 20)
          m_NormNHitHist[i]->Fill(nHits, normMean);
        m_NormNHitHist[m_nPart]->Fill(nHits, normMean);

        m_NormBothHist[i]->Fill(lambda, nHits, normMean);
        if (momT > 1)
          m_NormBothHist_1GeV[i]->Fill(lambda, nHits, normMean);

        if (momT > 1)
          m_NormLambdaFullHist[i]->Fill(lambda, normMean); // per species
        if (momT > 1)
          m_NormLambdaFullHist[m_nPart]->Fill(lambda, normMean); // all species

        if (momT > 1)
          m_NormCosThFullHist[i]->Fill(costh, normMean);
        if (momT > 1)
          m_NormCosThFullHist[m_nPart]->Fill(costh, normMean);

        // fiducial electrons
        if (i == 0 && m_fidMomMin <= mom && mom <= m_fidMomMax && m_fidLamMin <= lambda && lambda <= m_fidLamMax &&
            nHits >= m_fidNHitsMin)
          m_ResNorm_FiducialElectrons->Fill(normMean);

        if (i == 1 && normMean > 1.5)
          debug() << v << ": " << m_dQdxVec[i][v][0] << " " << m_dQdxVec[i][v][1] << " " << m_dQdxVec[i][v][2] //
                  << " " << "" << m_dQdxVec[i][v][3] << " " << m_dQdxVec[i][v][4] << " " << m_dQdxVec[i][v][5] //
                  << " " << normMean << " " << relError << endmsg;
      }
    }
    debug() << "NormRes done" << endmsg;

    m_ResNormHist[i]->Fit("gaus", "Q");
    m_ResNormHist_1GeV[i]->Fit("gaus", "Q");

    // angular distribution
    FitAndGetResolution(m_NormLambdaHist[i], gausfitrel, m_FitLambda[i], m_FitLambda_Res[i]);
    m_FitLambda_Res[i]->Fit("lamfit", "+Q", "", m_minBinLam, 55);

    // number-of-hits distribution
    FitAndGetResolution(m_NormNHitHist[i], gausfitrel, m_FitNHit[i], m_FitNHit_Res[i]);
    ((TH1D*)m_FitNHit[i]->At(1))->Fit(nhitcorr, "Q");
    m_FitNHit_Res[i]->Fit(nhitfit, "+Q", "", m_maxBinHits / 10., m_maxBinHits);

    m_NormBothHist[i]->FitSlicesZ(gausfitrel, 0, -1, 0, -1, 0, "QR"); // angular and number-of-hits distribution
    TH2* FitBoth_mean = (TH2D*)gDirectory->Get(m_NormBothHistNames[i][1].c_str());
    TH2* FitBoth_sigm = (TH2D*)gDirectory->Get(m_NormBothHistNames[i][2].c_str());
    TH2* ProjZ = dynamic_cast<TH2*>(m_NormBothHist[i]->Project3D("yx"));

    if (FitBoth_mean && FitBoth_sigm) // resolution-per-bin from double-correlation histogram
    {
      for (int d = 1; d <= m_nBinsLam; d++)
        for (int e = 1; e <= m_nBinsHits; e++) // loop over lambda and nHits
        {
          double mu = FitBoth_mean->GetBinContent(d, e);
          double sg = FitBoth_sigm->GetBinContent(d, e);
          double mue = FitBoth_mean->GetBinError(d, e);
          double sge = FitBoth_sigm->GetBinError(d, e);
          double res = sg / mu;
          double reserr = sqrt(pow(sge / mu, 2) + pow(sg / mu / mu * mue, 2));

          if (std::isnan(res) || std::isnan(reserr) || ProjZ->GetBinContent(d, e) == 0) {
            res = 0;
            reserr = 0;
          }
          m_FitBoth_Res[i]->SetBinContent(d, e, res);
          m_FitBoth_Res[i]->SetBinError(d, e, reserr);
        }
    }
    if (!FitBoth_mean)
      info() << "NormBothHist_mean_" << m_Particles[i] << " could not be retrieved." << endmsg;
    if (!FitBoth_sigm)
      info() << "NormBothHist_sigma_" << m_Particles[i] << " could not be retrieved." << endmsg;

    m_FitBoth_Res[i]->Fit(powla2, "Q");

    // angular distribution in lambda, no nHits cut
    FitAndGetResolution(m_NormLambdaFullHist[i], gausfitrel, m_FitLambdaFull[i], m_FitLambdaFull_Res[i]);
    ((TH1D*)m_FitLambdaFull[i]->At(1))->Fit(poly3fit, "Q");

    // angular distribution in cos(theta), no nHits cut
    FitAndGetResolution(m_NormCosThFullHist[i], gausfitrel, m_FitCosThFull[i], m_FitCosThFull_Res[i]);
    TH1* FitCosThFull_mean = (TH1D*)m_FitCosThFull[i]->At(1);
    FitCosThFull_mean->Scale(1. / FitCosThFull_mean->GetBinContent(1));
    FitCosThFull_mean->Fit(angcorr, "Q");

    // fit the HE histograms for each beta*gamma bin, put fit result histograms in _HEFitHist
    m_HEHist[i]->FitSlicesY(landaufit, 0, -1, 0, "QR", m_HEFitHist[i]);
  }

  for (int i = 0; i < m_nPart; i++) // bin-wise separation power (by formula)
  {
    for (int j = i + 1; j < m_nPart; j++) // loop over the particle to compare to
    {
      TH1* FitHist1_mean = (TH1D*)m_FitHist[i]->At(1);
      TH1* FitHist1_sigm = (TH1D*)m_FitHist[i]->At(2);
      TH1* FitHist2_mean = (TH1D*)m_FitHist[j]->At(1);
      TH1* FitHist2_sigm = (TH1D*)m_FitHist[j]->At(2);
      for (int n = 1; n < m_nBinsMom + 1; n++) // 0 is underflow bin, nBinsx+1 overflow bin
      {
        // get mean and sigma as well as their errors to calculate sep power
        double mu1 = FitHist1_mean->GetBinContent(n);
        double mu2 = FitHist2_mean->GetBinContent(n);
        double sg1 = FitHist1_sigm->GetBinContent(n);
        double sg2 = FitHist2_sigm->GetBinContent(n);
        double mu1e = FitHist1_mean->GetBinError(n);
        double mu2e = FitHist2_mean->GetBinError(n);
        double sg1e = FitHist1_sigm->GetBinError(n);
        double sg2e = FitHist2_sigm->GetBinError(n);
        debug() << mu1 << " " << mu1e << " " << mu2 << " " << mu2e << " " << sg1 << " " << sg1e << " " << sg2 << " "
                << sg2e << endmsg;

        // sep power = abs(mu1 - mu2) / sqrt( .5*( sigma1^2 + sigma2^2 ) )
        // make sure denominator is not 0; error calculated following simple error propagation
        double dm = fabs(mu1 - mu2);
        double ds = pow(sg1, 2) + pow(sg1, 2);
        double sds = sqrt(.5 * ds);
        double sp = (sds == 0) ? 0 : dm / sds; // separation power
        double spe_1 = 2 * pow(mu1e, 2) / ds;
        double spe_2 = 2 * pow(mu2e, 2) / ds;
        double spe_3 = 2 * pow(dm, 2) * pow(sg1, 2) * pow(sg1e, 2) / pow(ds, 3);
        double spe_4 = 2 * pow(dm, 2) * pow(sg2, 2) * pow(sg2e, 2) / pow(ds, 3);
        double spe = (ds == 0) ? 0 : sqrt(spe_1 + spe_2 + spe_3 + spe_4); // separation power error
        debug() << "  " << sp << " " << sqrt(spe_1) << " " << sqrt(spe_2) << " " << sqrt(spe_3) << " " << sqrt(spe_4)
                << " " << spe << endmsg;

        m_SPHist[sp_c]->SetBinContent(n, sp);
        m_SPHist[sp_c]->SetBinError(n, spe);
      }
      sp_c++;
    } // -for (int j=i+1; j<m_nPart; j++)
  } //  - for (int i=0; i<m_nPart; i++)

  TFitResultPtr ptr_fe = m_ResNorm_FiducialElectrons->Fit("gaus", "QS");

  // additional histograms for 'all' particles
  // number-of-hits distribution
  FitAndGetResolution(m_NormNHitHist[m_nPart], gausfitrel, m_FitNHit[m_nPart], m_FitNHit_Res[m_nPart]);
  TH1* FitNHit_mean = (TH1D*)m_FitNHit[m_nPart]->At(1);
  TFitResultPtr ptr_nhit = FitNHit_mean->Fit(nhitcorr, "+QS");

  // angular distribution in lambda, no nHits cut
  FitAndGetResolution(m_NormLambdaFullHist[m_nPart], gausfitrel, m_FitLambdaFull[m_nPart],
                      m_FitLambdaFull_Res[m_nPart]);
  TH1* FitLambdaFull_mean = (TH1D*)m_FitLambdaFull[m_nPart]->At(1);
  TFitResultPtr ptr_lm_pol2 = FitLambdaFull_mean->Fit(poly2fit, "+QS");
  TFitResultPtr ptr_lm_pol3 = FitLambdaFull_mean->Fit(poly3fit, "+QS");

  // angular distribution in cos(theta), no nHits cut
  FitAndGetResolution(m_NormCosThFullHist[m_nPart], gausfitrel, m_FitCosThFull[m_nPart], m_FitCosThFull_Res[m_nPart]);
  TH1* FitCosThFull_mean = (TH1D*)m_FitCosThFull[m_nPart]->At(1);
  FitCosThFull_mean->Scale(1. / FitCosThFull_mean->GetBinContent(1));
  angcorr->SetParameter(0, -.1);

  double tot_diff = 0, tot_weight = 0;
  for (int i_th = 1; i_th < FitCosThFull_mean->GetNbinsX() + 1; i_th++) // remaining cos(theta) deviation from flat
  {
    double diff = fabs(FitCosThFull_mean->GetBinContent(i_th) - 1);
    double weight = 1. / FitCosThFull_mean->GetBinError(i_th);
    tot_diff += diff * weight;
    tot_weight += weight;
  }
  double CosThDeviation = tot_diff / tot_weight;
  TFitResultPtr ptr_ac = FitCosThFull_mean->Fit(angcorr, "+QS");

  // calculate p-value-based separation power, with and without angular correction
  std::vector<std::vector<double>> dQdxVec[m_nPart];
  std::vector<std::vector<double>> dQdxVec_angcorr[m_nPart];

  for (int i = 0; i < m_nPart; i++) // create separate dQ/dx vectors for each particle species
  {
    dQdxVec[i] = std::vector<std::vector<double>>(m_nBinsMom, std::vector<double>(0));
    dQdxVec_angcorr[i] = std::vector<std::vector<double>>(m_nBinsMom, std::vector<double>(0));
    for (unsigned int v = 0; v < m_dQdxVec[i].size(); v++) // angular correction
    {
      double lambda = fabs(atan(m_dQdxVec[i][v][3])) * 180 / M_PI;
      int lambin = FitLambdaFull_mean->FindBin(lambda);
      double lamMean = FitLambdaFull_mean->GetBinContent(lambin);
      double ancorrMean = m_dQdxVec[i][v][1] / lamMean;

      for (int j = 0; j < m_nBinsMom; ++j) // per momentum bin
      {
        if (m_dQdxVec[i][v][0] > m_histBinsMom[j] && m_dQdxVec[i][v][0] < m_histBinsMom[j + 1]) {
          dQdxVec[i][j].emplace_back(m_dQdxVec[i][v][1]);
          dQdxVec_angcorr[i][j].emplace_back(ancorrMean);
        }
      }
    }
    for (int j = 0; j < m_nBinsMom; ++j) // sort dQ/dx
    {
      std::sort(dQdxVec[i][j].begin(), dQdxVec[i][j].end());
      std::sort(dQdxVec_angcorr[i][j].begin(), dQdxVec_angcorr[i][j].end());
    }
  }

  sp_c = 0;
  for (int i = 0; i < m_nPart; i++) // particle 1
  {
    for (int j = i + 1; j < m_nPart; j++) // particle 2
    {
      for (int k = 0; k < m_nBinsMom; ++k) // momentum bin
      {
        double S = Calculate_pValueSP(dQdxVec[i][k], dQdxVec[j][k]);
        m_SPHist_pValue[sp_c]->SetBinContent(k, S);
        double S2 = Calculate_pValueSP(dQdxVec_angcorr[i][k], dQdxVec_angcorr[j][k]);
        m_SPHist_pValue_angcorr[sp_c]->SetBinContent(k, S2);
      }
      sp_c++;
    }
  }

  // write histograms to root file
  m_outFile->cd();

  for (int i = 0; i < m_nPart; i++) // histograms for target particle species only
  {
    m_FitHist_Res[i]->Write();
    m_FitHist_ResSum[i]->Write();
    m_ResNormHist[i]->Write();
    m_ResNormHist_1GeV[i]->Write();

    m_NormLambdaHist[i]->Write();
    m_FitLambda_Res[i]->Write();

    m_NormBothHist[i]->Write();
    m_NormBothHist_1GeV[i]->Write();
    m_FitBoth_Res[i]->Write();
    m_FitBoth_Res_1GeV[i]->Write();

    for (int j = 0; j < 4; j++) // fit results with 3 fit parameters + chi square
    {
      m_FitHist[i]->At(j)->Write();
      m_FitLambda[i]->At(j)->Write();
      m_FitNHit[i]->At(j)->Write();
      m_FitLambdaFull[i]->At(j)->Write();
      m_FitCosThFull[i]->At(j)->Write();

      m_HEFitHist[i]->At(j)->Write();
    }
  }

  for (int i = 0; i < m_nPart + 1; i++) // histograms with 'other' or 'all' particles
  {
    m_NormNHitHist[i]->Write();
    m_FitNHit_Res[i]->Write();

    m_NormLambdaFullHist[i]->Write();
    m_FitLambdaFull_Res[i]->Write();
    m_NormCosThFullHist[i]->Write();
    m_FitCosThFull_Res[i]->Write();
  }

  for (int i = 0; i < m_nPart + 2; i++) // histograms with both 'other' and 'all'
  {
    m_BBHist[i]->Write();

    m_HitEnergySpectrum[i]->Write();
    m_HitNumberSpectrum[i]->Write();
    m_HitEnergySpectrum_at1[i]->Write();
    m_HitEnergySpectrum_at25[i]->Write();
    m_HEHist[i]->Write();
  }

  for (int j = 0; j < 4; j++) // more fit results
  {
    m_FitNHit[m_nPart]->At(j)->Write();
    m_FitLambdaFull[m_nPart]->At(j)->Write();
    m_FitCosThFull[m_nPart]->At(j)->Write();
  }

  for (int i = 0; i < sp_c; i++) // separation power
  {
    m_SPHist[i]->Write();
    m_SPHist_pValue[i]->Write();
    m_SPHist_pValue_angcorr[i]->Write();
  }

  m_ResElectrons->Write();
  m_ResNorm_FiducialElectrons->Write();

  m_MCmult->Write();

  m_CutAnaHist->Write();
  m_NTracksHist->Write();
  m_NTracksUsedHist->Write();
  m_NTracksMomHist->Write();
  m_TrackMomHist->Write();
  m_TrackSiblingsMomHist->Write();
  m_dEdxOutlierAbundance->Write();
  m_dEdxOutlierProportion->Write();

  m_dQdxErrorCorrelation->Write();

  if (m_plotStuff) // plot stuff
  {
    TCanvas* can = new TCanvas;
    gStyle->SetPalette(kBird);
    can->SetGrid();
    can->SetLogx(); // start with log x-axis plots

    // some y-axis max values for specific histogram plots
    std::vector<double> YMax = {1000, m_maxBindQdx, m_maxBindQdx / 10., 10};           // dQ/dx fit results (BBHist)
    std::vector<double> YMaxHE = {1e6, m_maxBindQdx * 6, m_maxBindQdx * 6 / 10., 100}; // hit energy fit results

    for (int i = 0; i < m_nPart; i++) // 2D Bethe-Bloch curves
    {
      PlotTH2(can, m_BBHist[i]);
      m_FitHist_Res[i]->SetAxisRange(0, .2, "Y");
      PlotTH1(can, m_FitHist_Res[i]);

      for (int j = 0; j < 4; j++) // fits
      {
        TH1D* h = (TH1D*)m_FitHist[i]->At(j);
        h->SetAxisRange(0, YMax[j], "Y");
        PlotTH1_FPn(can, h, j);
      }
      if (ptr_BBfit[i].Get())
        ptr_BBfit[i]->Print();
      else
        info() << "Mass " << m_Masses[i] << ": fit result pointer empty!" << endmsg;
    }
    PlotTH2(can, m_BBHist[m_nPart]);
    PlotTH2(can, m_BBHist[m_nPart + 1]);

    for (int i = 0; i < sp_c; i++) // separation powers
    {
      m_SPHist[i]->SetAxisRange(0, 20, "Y");
      PlotTH1(can, m_SPHist[i]);
      m_SPHist_pValue[i]->SetAxisRange(0, 20, "Y");
      PlotTH1(can, m_SPHist_pValue[i]);
      m_SPHist_pValue_angcorr[i]->SetAxisRange(0, 20, "Y");
      PlotTH1(can, m_SPHist_pValue_angcorr[i]);
    }

    for (int i = 0; i < m_nPart; i++) // hit energy
    {
      PlotTH2(can, m_HEHist[i]);
      for (int j = 0; j < 4; j++) // fits
      {
        TH1D* h = (TH1D*)m_HEFitHist[i]->At(j);
        h->SetAxisRange(0, YMaxHE[j], "Y");
        PlotTH1_FPn(can, h, j);
      }
    }

    PlotTH2(can, m_HEHist[m_nPart]);
    PlotTH2(can, m_HEHist[m_nPart + 1]);

    PlotTH1(can, m_TrackMomHist);
    PlotTH1(can, m_TrackSiblingsMomHist);

    // linear x-axis plots
    can->SetLogx(0);

    PlotTH1(can, m_ResElectrons);

    for (int i = 0; i < m_nPart; i++) // resolution from normalised dQ/dx
    {
      PlotTH1(can, m_FitHist_ResSum[i]);

      gStyle->SetOptFit(1);
      PlotTH1(can, m_ResNormHist[i]);
      PlotTH1(can, m_ResNormHist_1GeV[i]);

      PlotTH2(can, m_NormLambdaHist[i]);
      m_FitLambda_Res[i]->SetAxisRange(0, .2, "Y");
      PlotTH1(can, m_FitLambda_Res[i]);

      PlotTH2(can, m_NormNHitHist[i]);
      m_FitNHit_Res[i]->SetAxisRange(0, .2, "Y");
      PlotTH1(can, m_FitNHit_Res[i]);
      gStyle->SetOptFit(0);

      m_FitBoth_Res[i]->SetAxisRange(0, .2, "Z");
      PlotTH2(can, m_FitBoth_Res[i]);
      m_FitBoth_Res_1GeV[i]->SetAxisRange(0, .2, "Z");
      PlotTH2(can, m_FitBoth_Res_1GeV[i]);

      PlotTH2(can, m_NormLambdaFullHist[i]);
      m_FitLambdaFull_Res[i]->SetAxisRange(0, .2, "Y");
      PlotTH1(can, m_FitLambdaFull_Res[i]);
    }
    PlotTH2(can, m_NormLambdaFullHist[m_nPart]);
    m_FitLambdaFull_Res[m_nPart]->SetAxisRange(0, .2, "Y");
    PlotTH1(can, m_FitLambdaFull_Res[m_nPart]);

    gStyle->SetOptFit(1);

    PlotTH1(can, FitNHit_mean);
    if (ptr_nhit.Get())
      ptr_nhit->Print();
    else
      info() << FitNHit_mean->GetName() << ": fit result pointer empty!" << endmsg;

    PlotTH1(can, FitLambdaFull_mean);
    if (ptr_lm_pol2.Get())
      ptr_lm_pol2->Print();
    else
      info() << FitLambdaFull_mean->GetName() << ": poly2 fit result pointer empty!" << endmsg;
    if (ptr_lm_pol3.Get())
      ptr_lm_pol3->Print();
    else
      info() << FitLambdaFull_mean->GetName() << ": poly3 fit result pointer empty!" << endmsg;

    PlotTH1(can, FitCosThFull_mean);
    if (ptr_ac.Get())
      ptr_ac->Print();
    else
      info() << FitCosThFull_mean->GetName() << ": fit result pointer empty!" << endmsg;

    PlotTH1(can, m_ResNorm_FiducialElectrons);

    gStyle->SetOptFit(0);

    PlotTH1(can, m_CutAnaHist);
    PlotTH2(can, m_dQdxErrorCorrelation);

    PlotTH1(can, m_dEdxOutlierAbundance);
    for (int j = 1; j < 5; j++)
      m_dEdxOutlierProportion->SetBinContent(j, m_dEdxOutlierAbundance->GetBinContent(j) / m_BBHist[j]->GetEntries());
    PlotTH1(can, m_dEdxOutlierProportion);

    PlotTH1(can, m_NTracksHist);
    PlotTH1(can, m_NTracksUsedHist);
    PlotTH1(can, m_NTracksMomHist);

    PlotTH1(can, m_MCmult);

    for (int i = 0; i < m_nPart + 2; i++) // hit energy spectra
    {
      PlotTH1(can, m_HitEnergySpectrum[i]);
      PlotTH1(can, m_HitNumberSpectrum[i]);
      PlotTH1(can, m_HitEnergySpectrum_at1[i]);
      PlotTH1(can, m_HitEnergySpectrum_at25[i]);
    }

    delete can;
  }

  // some specific text output, both in console and in calibration file
  std::ofstream outfile;
  outfile.open(m_outCalibFile);

  info() << endmsg;
  info() << "===========  Results  ===========" << endmsg;
  info() << endmsg;

  if (ptr_fe.Get()) {
    info() << "Fiducial electrons dE/dx resolution: " << ptr_fe->Parameter(2) / ptr_fe->Parameter(1) * 100 << " %"
           << endmsg;
    info() << endmsg;
    outfile << "# Fiducial electrons dE/dx resolution:\n" << ptr_fe->Parameter(2) / ptr_fe->Parameter(1) << std::endl;
  } else {
    info() << "Fiducial electrons fit result pointer is empty!" << endmsg;
    info() << endmsg;
    outfile << "# Fiducial electrons fit result pointer is empty!" << std::endl;
  }

  if (ptr_lm_pol3.Get()) {
    info() << "Angular dependence fit (poly3 vs. lambda) result:" << endmsg;
    info() << ptr_lm_pol3->Parameter(0) << " " << ptr_lm_pol3->Parameter(1) << " " //
           << ptr_lm_pol3->Parameter(2) << " " << ptr_lm_pol3->Parameter(3) << endmsg;
    info() << "with a chi2/ndf of  " << ptr_lm_pol3->Chi2() << "/" << ptr_lm_pol3->Ndf() << endmsg;
    outfile << "# Angular dependence fit (poly3 vs. lambda) result:" << std::endl;
    outfile << ptr_lm_pol3->Parameter(0) << " " << ptr_lm_pol3->Parameter(1) << " " //
            << ptr_lm_pol3->Parameter(2) << " " << ptr_lm_pol3->Parameter(3) << "\n# with a chi2/ndf of\n"
            << ptr_lm_pol3->Chi2() << "/" << ptr_lm_pol3->Ndf() << std::endl;
  } else {
    info() << "Angular dependence fit result pointer is empty!" << endmsg;
    outfile << "# Angular dependence fit result pointer is empty!" << std::endl;
  }

  info() << "Remaining deviation from flat distribution: " << CosThDeviation * 100 << " %" << endmsg;
  info() << "  note: this should be below 1 %" << endmsg;
  info() << endmsg;
  outfile << "# Remaining deviation from flat distribution (note: this should be below 0.01): \n"
          << CosThDeviation << std::endl;

  info() << "Bethe-Bloch fit results, used as input for CPID:" << endmsg;
  info() << "  note: par2 and par3 are redundant, so par3=1 is chosen" << endmsg;
  info() << "mass      chi2/ndf    parameters 1-5" << endmsg;
  outfile << "# Bethe-Bloch fit results, used as input for CPID:" << std::endl;
  outfile << "#  note: par2 and par3 are redundant, so par3=1 is chosen" << std::endl;
  outfile << "# mass    chi2/ndf    parameters 1-5" << std::endl;
  for (int i = 0; i < m_nPart; ++i) // fit results per species
  {
    if (ptr_BBfit[i].Get()) {
      info() << m_Masses[i] << "  " << ptr_BBfit[i]->Chi2() << "/" << ptr_BBfit[i]->Ndf() << "  "
             << ptr_BBfit[i]->Parameter(1) << " " << ptr_BBfit[i]->Parameter(2) << " 1 " //
             << ptr_BBfit[i]->Parameter(3) << " " << ptr_BBfit[i]->Parameter(4) << " " << endmsg;
      outfile << m_Masses[i] << "  " << ptr_BBfit[i]->Chi2() << "/" << ptr_BBfit[i]->Ndf() << "  "
              << ptr_BBfit[i]->Parameter(1) << " " << ptr_BBfit[i]->Parameter(2) << " 1 " //
              << ptr_BBfit[i]->Parameter(3) << " " << ptr_BBfit[i]->Parameter(4) << " " << std::endl;
    } else {
      info() << m_Masses[i] << " fit result pointer is empty!" << endmsg;
      outfile << "# " << m_Masses[i] << " fit result pointer is empty!" << std::endl;
    }
  }
  info() << endmsg;

  outfile.close();

  // some final debug output
  debug() << "CutAnaHist Entries" << endmsg;
  for (unsigned int j = 0; j < 10; j++)
    debug() << j << "  " << m_CutAnaHist->GetBinContent(j) << endmsg;
  debug() << "NTracksHist Entries" << endmsg;
  for (unsigned int j = 0; j < 10; j++)
    debug() << j << "  " << m_NTracksHist->GetBinContent(j) << endmsg;
  debug() << "NTracksMomHist Entries" << endmsg;
  for (unsigned int j = 0; j < 10; j++)
    debug() << j << "  " << m_NTracksMomHist->GetBinContent(j) << endmsg;
  debug() << "TrackMomHist Entries" << endmsg;
  for (int j = 0; j < m_nBinsMom; j++)
    debug() << j << "  " << m_histBinsMom2[j] << "  " << m_TrackMomHist->GetBinContent(j) << endmsg;

  m_outFile->Close();

  delete gausfit;
  delete gausfitrel;
  delete powla2;
  delete landaufit;
  delete poly2fit;
  delete poly3fit;
  delete angcorr;
  delete nhitcorr;
  delete BBfit;

  delete[] m_histBinsMom;
  delete[] m_histBinsMom2;
  return StatusCode::SUCCESS;
}

double dQdxAnalyser::Calculate_pValueSP(std::vector<double> v1, std::vector<double> v2) // p-value based sep power
{
  int n1 = v1.size(), n2 = v2.size();
  if (n1 == 0 || n2 == 0)
    return 0;
  int c1 = 0, c2 = 0, cuml1 = 0, cuml2 = 0;
  double frac1 = 0, frac2 = 0, oldfrac1 = 0, oldfrac2 = 0;
  double p = 1, Z = 0;

  while (true) // main loop
  {
    if (c1 == n1 && c2 == n2)
      return 0;
    else if (c1 == n1) {
      ++cuml2;
      ++c2;
    } else if (c2 == n2) {
      ++cuml1;
      ++c1;
    } else if (v1[c1] < v2[c2]) {
      ++cuml1;
      ++c1;
    } else {
      ++cuml2;
      ++c2;
    }

    frac1 = 1. * cuml1 / n1;
    frac2 = 1. * cuml2 / n2;

    if (frac1 < frac2 && frac1 > 1 - frac2) //
    {
      if (frac1 / (1 - frac2) < (1 - oldfrac2) / oldfrac1)
        p = pow(frac1 * (1 - frac2), .5);
      else
        p = pow(oldfrac1 * (1 - oldfrac2), .5);
      debug() << n1 << " " << n2 << " " << c1 << " " << c2 << " " << frac1 << " " << 1 - frac2 << " " << oldfrac1 << " "
              << 1 - oldfrac2 << " " << p << " " << 2 * ROOT::Math::gaussian_quantile(1 - p, 1) << endmsg;
      Z = 2 * ROOT::Math::gaussian_quantile(1 - p, 1);
      if (isnan(Z) || isinf(Z))
        return 0;
      else
        return Z;
    } else if (frac2 < frac1 && frac2 > 1 - frac1) //
    {
      if (frac2 / (1 - frac1) < (1 - oldfrac1) / oldfrac2)
        p = pow(frac2 * (1 - frac1), .5);
      else
        p = pow(oldfrac2 * (1 - oldfrac1), .5);
      debug() << n1 << " " << n2 << " " << c1 << " " << c2 << " " << 1 - frac1 << " " << frac2 << " " << 1 - oldfrac1
              << " " << oldfrac2 << " " << p << " " << 2 * ROOT::Math::gaussian_quantile(1 - p, 1) << endmsg;
      Z = 2 * ROOT::Math::gaussian_quantile(1 - p, 1);
      if (isnan(Z) || isinf(Z))
        return 0;
      else
        return Z;
    }

    oldfrac1 = frac1;
    oldfrac2 = frac2;
  }
}

// helper functions for momentum, fitting and histogram plotting
std::tuple<float, float> dQdxAnalyser::getMomFromTrackState(edm4hep::TrackState state) const // state -> helix -> mom
{
  HelixClass helix;
  helix.Initialize_Canonical(state.phi, state.D0, state.Z0, state.omega, state.tanLambda, this->m_Bz);
  float mom3 = sqrt(pow(helix.getMomentum()[0], 2) + pow(helix.getMomentum()[1], 2) + pow(helix.getMomentum()[2], 2));
  float momT = sqrt(pow(helix.getMomentum()[0], 2) + pow(helix.getMomentum()[1], 2));
  return {mom3, momT};
}

void dQdxAnalyser::FitAndGetResolution(TH2* hist, TF1* fitFunction, TObjArray* fitResult, TH1* resolution) //
{
  hist->FitSlicesY(fitFunction, 0, -1, 0, "QR", fitResult);
  TH1* fitResult_mean = (TH1D*)fitResult->At(1);
  TH1* fitResult_sigm = (TH1D*)fitResult->At(2);
  resolution->Divide(fitResult_sigm, fitResult_mean);
}

void dQdxAnalyser::PlotTH1(TCanvas* can, TH1* hist) // TH1
{
  if (!hist) {
    info() << "Histogramm is NULL!" << endmsg;
    return;
  }
  hist->Draw();
  can->Update();
  std::string s = m_plotFolder + "/" + hist->GetName() + m_plotFormat;
  can->Print(s.c_str());
}
void dQdxAnalyser::PlotTH2(TCanvas* can, TH2* hist) // TH2
{
  if (!hist) {
    info() << "Histogramm is NULL!" << endmsg;
    return;
  }
  hist->Draw("colz");
  can->Update();
  std::string s = m_plotFolder + "/" + hist->GetName() + m_plotFormat;
  can->Print(s.c_str());
}
void dQdxAnalyser::PlotTH1_FPn(TCanvas* can, TH1* hist, int n) // TH1 with fit parameter n
{
  if (!hist) {
    info() << "Histogramm is NULL!" << endmsg;
    return;
  }
  hist->Draw();
  can->Update();
  std::string s = m_plotFolder + "/" + hist->GetName() + "_FitPar" + std::to_string(n) + m_plotFormat;
  can->Print(s.c_str());
}
