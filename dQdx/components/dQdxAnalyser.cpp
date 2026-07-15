#include <dQdxAnalyser.h>
#include <GaudiKernel/MsgStream.h>
#include <podio/LinkNavigator.h>
#include <edm4hep/utils/TrackUtils.h>
#include <edm4hep/TrackerHit.h>
#include <marlinutil/HelixClass.h>
//#include "DD4hep/Detector.h"
//#include <DDRec/DetectorData.h>
//#include <DD4hep/DD4hepUnits.h>


StatusCode dQdxAnalyser::initialize()
{
  // make sure output root file is available
  m_outFile = new TFile(m_outputFile.value().c_str(), "RECREATE");
  if (!m_outFile || m_outFile->IsZombie()) {
    error() << "Cannot open output file: " << m_outputFile.value() << endmsg;
    return StatusCode::FAILURE;
  }

  // prepare geometry info
  dd4hep::Detector& theDet = dd4hep::Detector::getInstance();
  double bfieldV[3];
  theDet.field().magneticField({ 0.,0.,0.}, bfieldV);
  m_Bz = bfieldV[2]/dd4hep::tesla;

  debug() << "B_z = " << m_Bz << endmsg;

  // prepare bins for histograms
  m_minBinMom2  = m_minBinMom - 1;
  m_maxBindQdx2 = m_maxBindQdx * 6;
  m_histBinsMom  = new float[m_nBinsMom+1];
  m_histBinsMom2 = new float[m_nBinsMom+1];
  float* histBinsdQdx  = new float[m_nBinsdQdx+1];
  float* histBinsdQdx2 = new float[m_nBinsdQdx+1];

  for (int i=0; i<m_nBinsMom+1;  i++) m_histBinsMom  [i] = pow( 10, m_minBinMom  + (m_maxBinMom-m_minBinMom )*i/(float)(m_nBinsMom) );
  for (int i=0; i<m_nBinsMom+1;  i++) m_histBinsMom2 [i] = pow( 10, m_minBinMom2 + (m_maxBinMom-m_minBinMom2)*i/(float)(m_nBinsMom) );
  for (int i=0; i<m_nBinsdQdx+1; i++)   histBinsdQdx [i] = m_minBindQdx + (m_maxBindQdx -m_minBindQdx)*i/(float)(m_nBinsdQdx);
  for (int i=0; i<m_nBinsdQdx+1; i++)   histBinsdQdx2[i] = m_minBindQdx + (m_maxBindQdx2-m_minBindQdx)*i/(float)(m_nBinsdQdx);

  // create Bethe-Bloch histograms
  for (int i=0; i<m_nPart+2; i++)
  {
    std::stringstream s1; s1 << "BB" << m_Particles[i];
    std::stringstream s2; s2 << "Bethe-Bloch histogram for " << m_particles[i];
    m_BBHist[i] = new TH2D(s1.str().c_str(), s2.str().c_str(), m_nBinsMom,m_histBinsMom, m_nBinsdQdx,histBinsdQdx);
    m_BBHist[i]->SetXTitle("Momentum (GeV)");
    if (m_usedENdx.value() == 1) m_BBHist[i]->SetYTitle("dE/dx (GeV/mm)");
    if (m_usedENdx.value() == 2) m_BBHist[i]->SetYTitle("dN/dx (1/mm)");
  }

  // create space for fit results
  for (int i=0; i<m_nPart; i++)   m_FitHist[i]      = new TObjArray();
  for (int i=0; i<m_nPart; i++)   m_FitLambda[i]    = new TObjArray();
  for (int i=0; i<m_nPart+1; i++) m_FitNHit[i]      = new TObjArray();
  //for (int i=0; i<m_nPart; i++)   m_FitBoth[i]      = new TObjArray();
  //for (int i=0; i<m_nPart; i++)   m_FitBoth_1GeV[i] = new TObjArray();
  for (int i=0; i<m_nPart+1; i++) m_FitLambdaFull[i]= new TObjArray();
  for (int i=0; i<m_nPart+1; i++) m_FitCosThFull[i] = new TObjArray();
  for (int i=0; i<m_nPart; i++)   m_HEFitHist[i]    = new TObjArray();

  // create resolution and separation power histograms
  int sp_c = 0;
  for (int i=0; i<m_nPart; i++)
  {
    std::stringstream s1; s1 << "BB" << m_Particles[i] << "_FitRes";
    std::stringstream s2; s2 << "dQ/dx resolution over momentum for " << m_particles[i] << " from fit result";
    m_FitHist_Res[i] = new TH1D(s1.str().c_str(), s2.str().c_str(), m_nBinsMom, m_histBinsMom);
    m_FitHist_Res[i]->SetXTitle("momentum / GeV/c");
    m_FitHist_Res[i]->SetYTitle("dQ/dx resolution");

    std::stringstream s3; s3 << "BB" << m_Particles[i] << "_FitResSum";
    std::stringstream s4; s4 << "dQ/dx resolution distribution above 1 GeV for " << m_particles[i] << " from fit result";
    m_FitHist_ResSum[i] = new TH1D(s3.str().c_str(), s4.str().c_str(), 100, 0, .2);
    m_FitHist_ResSum[i]->SetXTitle("dQ/dx resolution");
    m_FitHist_ResSum[i]->SetYTitle("number of momentum bins with that resolution");

    std::stringstream s5; s5 << "Res" << m_Particles[i] << "_Norm";
    std::stringstream s6; s6 << "Relative dQ/dx resolution distribution for " << m_particles[i] << " from all entries above .3 GeV";
    m_ResNormHist[i] = new TH1D(s5.str().c_str(), s6.str().c_str(), 100, 0, 2);
    m_ResNormHist[i]->SetXTitle("relative dQ/dx");
    m_ResNormHist[i]->SetYTitle("rel. abundance");

    std::stringstream s7; s7 << "Res" << m_Particles[i] << "_Norm_1GeV";
    std::stringstream s8; s8 << "Relative dQ/dx resolution distribution for " << m_particles[i] << " from all entries above 1 GeV";
    m_ResNormHist_1GeV[i] = new TH1D(s7.str().c_str(), s8.str().c_str(), 100, 0, 2);
    m_ResNormHist_1GeV[i]->SetXTitle("relative dQ/dx");
    m_ResNormHist_1GeV[i]->SetYTitle("rel. abundance");

    // Resolution dependence on different parameters
    std::stringstream s9; s9 << "NormLambda" << m_Particles[i] << "";
    std::stringstream s10; s10 << "Relative dQ/dx distribution over lambda with p>1 GeV and nHit>=200 for " << m_particles[i] << " from all entries";
    m_NormLambdaHist[i] = new TH2D(s9.str().c_str(), s10.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, 100, 0, 2);
    m_NormLambdaHist[i]->SetXTitle("|lambda| / deg");
    m_NormLambdaHist[i]->SetYTitle("relative dQ/dx");

    std::stringstream s11; s11 << "NormLambda" << m_Particles[i] << "_FitRes";
    std::stringstream s12; s12 << "dQ/dx resolution for " << m_particles[i] << " over lambda from fit result";
    m_FitLambda_Res[i] = new TH1D(s11.str().c_str(), s12.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam);
    m_FitLambda_Res[i]->SetXTitle("|lambda| / deg");
    m_FitLambda_Res[i]->SetYTitle("dQ/dx resolution");

    std::stringstream s13; s13 << "NormNHit" << m_Particles[i] << "";
    std::stringstream s14; s14 << "Relative dQ/dx distribution over number of hits with p>1GeV and |lambda|<20deg for " << m_particles[i] << " from all entries";
    m_NormNHitHist[i] = new TH2D(s13.str().c_str(), s14.str().c_str(), m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5, 100, 0, 2);
    m_NormNHitHist[i]->SetXTitle("number of hits");
    m_NormNHitHist[i]->SetYTitle("relative dQ/dx");

    std::stringstream s15; s15 << "NormNHit" << m_Particles[i] << "_FitRes";
    std::stringstream s16; s16 << "dQ/dx resolution for " << m_particles[i] << " over number of hits from fit result";
    m_FitNHit_Res[i] = new TH1D(s15.str().c_str(), s16.str().c_str(), m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5);
    m_FitNHit_Res[i]->SetXTitle("number of hits");
    m_FitNHit_Res[i]->SetYTitle("dQ/dx resolution");

    std::stringstream s20; s20 << "NormBoth" << m_Particles[i] << "";
    std::stringstream s21; s21 << "Relative dQ/dx distribution over lambda and number of hits for " << m_particles[i] << " from all entries";
    m_NormBothHist[i] = new TH3D(s20.str().c_str(), s21.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5, 100, 0, 2);
    m_NormBothHist[i]->SetXTitle("|lambda| / deg");
    m_NormBothHist[i]->SetYTitle("number of hits");
    m_NormBothHist[i]->SetZTitle("relative dQ/dx");

    std::string HistName = s20.str();
    m_NormBothHistNames[i].push_back(HistName.append("_0"));
    m_NormBothHistNames[i].push_back(HistName.replace(HistName.end()-1,HistName.end(),"1"));
    m_NormBothHistNames[i].push_back(HistName.replace(HistName.end()-1,HistName.end(),"2"));
    m_NormBothHistNames[i].push_back(HistName.replace(HistName.end()-1,HistName.end(),"chi2"));

    std::stringstream s22; s22 << "NormBoth" << m_Particles[i] << "_FitRes";
    std::stringstream s23; s23 << "dQ/dx resolution for " << m_particles[i] << " over lambda and number of hits from fit result";
    m_FitBoth_Res[i] = new TH2D(s22.str().c_str(), s23.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5);
    m_FitBoth_Res[i]->SetXTitle("|lambda| / deg");
    m_FitBoth_Res[i]->SetYTitle("number of hits");
    m_FitBoth_Res[i]->SetZTitle("dQ/dx resolution");

    std::stringstream s24; s24 << "NormBoth" << m_Particles[i] << "_1GeV";
    std::stringstream s25; s25 << "Relative dQ/dx distribution over lambda and number of hits with p>1GeV for " << m_particles[i] << " from all entries";
    m_NormBothHist_1GeV[i] = new TH3D(s24.str().c_str(), s25.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5, 100, 0, 2);
    m_NormBothHist_1GeV[i]->SetXTitle("|lambda| / deg");
    m_NormBothHist_1GeV[i]->SetYTitle("number of hits");
    m_NormBothHist_1GeV[i]->SetZTitle("relative dQ/dx");

    HistName = s24.str();
    m_NormBothHistNames_1GeV[i].push_back(HistName.append("_amp"));
    m_NormBothHistNames_1GeV[i].push_back(HistName.replace(HistName.end()-3,HistName.end(),"mean"));
    m_NormBothHistNames_1GeV[i].push_back(HistName.replace(HistName.end()-4,HistName.end(),"sigma"));
    m_NormBothHistNames_1GeV[i].push_back(HistName.replace(HistName.end()-5,HistName.end(),"chi2"));

    std::stringstream s26; s26 << "NormBoth" << m_Particles[i] << "_FitRes_1GeV";
    std::stringstream s27; s27 << "dQ/dx resolution for " << m_particles[i] << " with p>1GeV over lambda and number of hits from fit result";
    m_FitBoth_Res_1GeV[i] = new TH2D(s26.str().c_str(), s27.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5);
    m_FitBoth_Res_1GeV[i]->SetXTitle("|lambda| / deg");
    m_FitBoth_Res_1GeV[i]->SetYTitle("number of hits");
    m_FitBoth_Res_1GeV[i]->SetZTitle("dQ/dx resolution");

    std::stringstream s28; s28 << "NormLambdaFull" << m_Particles[i] << "";
    std::stringstream s29; s29 << "Relative dQ/dx distribution over lambda with p>1 GeV for " << m_particles[i] << " from all entries";
    m_NormLambdaFullHist[i] = new TH2D(s28.str().c_str(), s29.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam, 100, 0, 2);
    m_NormLambdaFullHist[i]->SetXTitle("|lambda| / deg");
    m_NormLambdaFullHist[i]->SetYTitle("relative dQ/dx");

    std::stringstream s30; s30 << "NormLambdaFull" << m_Particles[i] << "_FitRes";
    std::stringstream s31; s31 << "dQ/dx resolution for " << m_particles[i] << " over lambda from fit result";
    m_FitLambdaFull_Res[i] = new TH1D(s30.str().c_str(), s31.str().c_str(), m_nBinsLam, m_minBinLam, m_maxBinLam);
    m_FitLambdaFull_Res[i]->SetXTitle("|lambda| / deg");
    m_FitLambdaFull_Res[i]->SetYTitle("dQ/dx resolution");

    std::stringstream s32; s32 << "NormCosThFull" << m_Particles[i] << "";
    std::stringstream s33; s33 << "Relative dQ/dx distribution over cos(theta) with p>1 GeV for " << m_particles[i] << " from all entries";
    m_NormCosThFullHist[i] = new TH2D(s32.str().c_str(), s33.str().c_str(), 20, 0, 1, 100, 0, 2);
    m_NormCosThFullHist[i]->SetXTitle("cos(theta)");
    m_NormCosThFullHist[i]->SetYTitle("relative dQ/dx");

    std::stringstream s34; s34 << "NormCosThFull" << m_Particles[i] << "_FitRes";
    std::stringstream s35; s35 << "dQ/dx resolution for " << m_particles[i] << " over cos(theta) from fit result";
    m_FitCosThFull_Res[i] = new TH1D(s34.str().c_str(), s35.str().c_str(), 20, 0, 1);
    m_FitCosThFull_Res[i]->SetXTitle("cos(theta)");
    m_FitCosThFull_Res[i]->SetYTitle("dQ/dx resolution");

    for (int j=i+1; j<m_nPart; j++)  // Separation Power
    {
      std::stringstream s17; s17 << "SP_" << m_SP_part[i] << "_" << m_SP_part[j];
      std::stringstream s18; s18 << "Separation power between " << m_particles[i] << " and " << m_particles[j];
      m_SPHist[sp_c] = new TH1D(s17.str().c_str(), s18.str().c_str(), m_nBinsMom, m_histBinsMom);
      m_SPHist[sp_c]->SetXTitle("momentum / GeV/c");
      m_SPHist[sp_c]->SetYTitle("separation power");

      std::stringstream s36; s36 << "SP_pValue_" << m_SP_part[i] << "_" << m_SP_part[j];
      std::stringstream s37; s37 << "Separation power by p-value between " << m_particles[i] << " and " << m_particles[j];
      m_SPHist_pValue[sp_c] = new TH1D(s36.str().c_str(), s37.str().c_str(), m_nBinsMom, m_histBinsMom);
      m_SPHist_pValue[sp_c]->SetXTitle("momentum / GeV/c");
      m_SPHist_pValue[sp_c]->SetYTitle("separation power");

      std::stringstream s46; s46 << "SP_pValue_angcorr_" << m_SP_part[i] << "_" << m_SP_part[j];
      std::stringstream s47; s47 << "Separation power by p-value with prior angular correction between " << m_particles[i] << " and " << m_particles[j];
      m_SPHist_pValue_angcorr[sp_c] = new TH1D(s46.str().c_str(), s47.str().c_str(), m_nBinsMom, m_histBinsMom);
      m_SPHist_pValue_angcorr[sp_c]->SetXTitle("momentum / GeV/c");
      m_SPHist_pValue_angcorr[sp_c]->SetYTitle("separation power");

      sp_c++;
    }
  }

  std::stringstream s13; s13 << "NormNHit" << m_Particles[m_nPart+1] << "";
  std::stringstream s14; s14 << "Relative dQ/dx distribution over number of hits with p>1GeV and |lambda|<20deg for " << m_particles[m_nPart+1] << " from all entries";
  m_NormNHitHist[m_nPart] = new TH2D(s13.str().c_str(), s14.str().c_str(), m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5, 100, 0, 2);
  m_NormNHitHist[m_nPart]->SetXTitle("number of hits");
  m_NormNHitHist[m_nPart]->SetYTitle("relative dQ/dx");

  std::stringstream s15; s15 << "NormNHit" << m_Particles[m_nPart+1] << "_FitRes";
  std::stringstream s16; s16 << "dQ/dx resolution for " << m_particles[m_nPart+1] << " over number of hits from fit result";
  m_FitNHit_Res[m_nPart] = new TH1D(s15.str().c_str(), s16.str().c_str(), m_nBinsHits, m_minBinHits+.5, m_maxBinHits+.5);
  m_FitNHit_Res[m_nPart]->SetXTitle("number of hits");
  m_FitNHit_Res[m_nPart]->SetYTitle("dQ/dx resolution");

  std::stringstream s38; s38 << "NormLambdaFull" << m_Particles[m_nPart+1] << "";
  std::stringstream s39; s39 << "Relative dQ/dx distribution over lambda with p>1 GeV for " << m_particles[m_nPart+1] << " from all entries";
  m_NormLambdaFullHist[m_nPart] = new TH2D(s38.str().c_str(), s39.str().c_str(), m_nBinsLam*2, m_minBinLam, m_maxBinLam, 100, 0, 2);
  m_NormLambdaFullHist[m_nPart]->SetXTitle("|lambda| / deg");
  m_NormLambdaFullHist[m_nPart]->SetYTitle("relative dQ/dx");

  std::stringstream s40; s40 << "NormLambdaFull" << m_Particles[m_nPart+1] << "_FitRes";
  std::stringstream s41; s41 << "dQ/dx resolution for " << m_particles[m_nPart+1] << " over lambda from fit result";
  m_FitLambdaFull_Res[m_nPart] = new TH1D(s40.str().c_str(), s41.str().c_str(), m_nBinsLam*2, m_minBinLam, m_maxBinLam);
  m_FitLambdaFull_Res[m_nPart]->SetXTitle("|lambda| / deg");
  m_FitLambdaFull_Res[m_nPart]->SetYTitle("dQ/dx resolution");

  std::stringstream s42; s42 << "NormCosThFull" << m_Particles[m_nPart+1] << "";
  std::stringstream s43; s43 << "Relative dQ/dx distribution over cos(theta) with p>1 GeV for " << m_particles[m_nPart+1] << " from all entries";
  m_NormCosThFullHist[m_nPart] = new TH2D(s42.str().c_str(), s43.str().c_str(), 20, 0, 1, 100, 0, 2);
  m_NormCosThFullHist[m_nPart]->SetXTitle("cos(theta)");
  m_NormCosThFullHist[m_nPart]->SetYTitle("relative dQ/dx");

  std::stringstream s44; s44 << "NormCosThFull" << m_Particles[m_nPart+1] << "_FitRes";
  std::stringstream s45; s45 << "dQ/dx resolution for " << m_particles[m_nPart+1] << " over cos(theta) from fit result";
  m_FitCosThFull_Res[m_nPart] = new TH1D(s44.str().c_str(), s45.str().c_str(), 20, 0, 1);
  m_FitCosThFull_Res[m_nPart]->SetXTitle("cos(theta)");
  m_FitCosThFull_Res[m_nPart]->SetYTitle("dQ/dx resolution");

  // fiducial electrons
  m_ResNorm_FiducialElectrons = new TH1D("ResNorm_FiducialElectrons", "Relative dQ/dx resolution distribution for fiducial electrons", 100, 0, 2);

  // some other histograms for cross-checks etc.
  m_ResElectrons = new TH1D("ResElectrons","dQ/dx resolution distribution above .3 GeV for electrons from all entries, not normalised", 100, m_minBindQdx, m_maxBindQdx/2.);
  m_MCmult = new TH1D("MCmult","MC particle multiplicity of linked reco objects (tracks or PFOs)",11,-.5,10.5);

  m_CutAnaHist = new TH1D("CutAnaHist","Contamination analyser histogram with given cuts",10,-.5,9.5);
  m_CutAnaHist->GetXaxis()->SetBinLabel(1,"Track exists");
  m_CutAnaHist->GetXaxis()->SetBinLabel(2,"MCPar and dQ/dx exist");
  m_CutAnaHist->GetXaxis()->SetBinLabel(3,"Track purity (1 track) passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(4,"Track state exists");
  m_CutAnaHist->GetXaxis()->SetBinLabel(5,"Momentum cuts passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(6,"Angular cuts passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(7,"d0/z0 cuts passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(8,"Track purity (highest mom) passed");
  m_CutAnaHist->GetXaxis()->SetBinLabel(9,"dE/dx outliers passed");

  m_NTracksHist     = new TH1D("NTracksHist"    ,"Number of tracks of each event",100,-.5,99.5);
  m_NTracksUsedHist = new TH1D("NTracksUsedHist","Number of tracks of each event actually used after cuts",100,-.5,99.5);
  m_NTracksMomHist  = new TH1D("NTracksMomHist" ,"Number of tracks of each event actually used after cuts and with p > 100 MeV",100,-.5,99.5);
  m_TrackMomHist = new TH1D("TrackMomHist","Momentum of tracks after cuts",m_nBinsMom,m_histBinsMom2);
  m_TrackMomHist->SetXTitle("momentum (GeV)");
  m_TrackMomHist->SetYTitle("abundance");
  m_TrackSiblingsMomHist = new TH1D("TrackSiblingsMomHist","Momentum of tracks after cuts in events with nTracks > 1",m_nBinsMom,m_histBinsMom2);
  m_TrackSiblingsMomHist->SetXTitle("momentum (GeV)");
  m_TrackSiblingsMomHist->SetYTitle("abundance");

  m_dEdxOutlierAbundance = new TH1D("dEdxOutlierAbundance","Number of tracks with dE/dx values that don't match the MC PDG",5,-.5,4.5);
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(1,"MCPion - dE/dx too large");
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(2,"MCKaon - dE/dx too large");
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(3,"MCKaon - dE/dx too small");
  m_dEdxOutlierAbundance->GetXaxis()->SetBinLabel(4,"MCProton - dE/dx too small");
  m_dEdxOutlierAbundance->SetYTitle("abundance");

  m_dEdxOutlierProportion = new TH1D("dEdxOutlierProportion","Proportion of tracks with dE/dx values that don't match the MC PDG",5,-.5,4.5);
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(1,"MCPion - dE/dx too large");
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(2,"MCKaon - dE/dx too large");
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(3,"MCKaon - dE/dx too small");
  m_dEdxOutlierProportion->GetXaxis()->SetBinLabel(4,"MCProton - dE/dx too small");
  m_dEdxOutlierProportion->SetYTitle("proportion of MCSpecies");

  m_dQdxErrorCorrelation = new TH2D("dQdxErrorCorrelation","Correlation plot of dQ/dx error vs. dQ/dx as cross check",100,m_minBindQdx,m_maxBindQdx,100,m_minBindQdx*.1,m_maxBindQdx*.1);

  // hit number and energy spectra
  for (int i=0; i<m_nPart+2; i++)
  {
    std::stringstream s3; s3 << "HitNum_" << m_Particles[i];
    std::stringstream s4; s4 << "Hit number spectrum for " << m_particles[i];
    m_HitNumberSpectrum[i] = new TH1D(s3.str().c_str(), s4.str().c_str(), m_maxBinHits, 0, m_maxBinHits);
    m_HitNumberSpectrum[i]->SetXTitle("Hit Number");
    m_HitNumberSpectrum[i]->SetYTitle("Abundance");
  }
  for (int i=0; i<m_nPart+2; i++)
  {
    std::stringstream s1; s1 << "HitEng_" << m_Particles[i];
    std::stringstream s2; s2 << "Hit energy spectrum for " << m_particles[i];
    m_HitEnergySpectrum[i] = new TH1D(s1.str().c_str(), s2.str().c_str(),100,0,1e-5);
    m_HitEnergySpectrum[i]->SetXTitle("Hit Energy (GeV)");
    m_HitEnergySpectrum[i]->SetYTitle("Abundance");

    std::stringstream s5; s5 << "HitEng_at1_" << m_Particles[i];
    std::stringstream s6; s6 << "Hit energy spectrum around 1 GeV for " << m_particles[i];
    m_HitEnergySpectrum_at1[i] = new TH1D(s5.str().c_str(), s6.str().c_str(),100,0,1e-5);
    m_HitEnergySpectrum_at1[i]->SetXTitle("Hit Energy (GeV)");
    m_HitEnergySpectrum_at1[i]->SetYTitle("Abundance");
    std::stringstream s7; s7 << "HitEng_at25_" << m_Particles[i];
    std::stringstream s8; s8 << "Hit energy spectrum around 2.5 GeV for " << m_particles[i];
    m_HitEnergySpectrum_at25[i] = new TH1D(s7.str().c_str(), s8.str().c_str(),100,0,1e-5);
    m_HitEnergySpectrum_at25[i]->SetXTitle("Hit Energy (GeV)");
    m_HitEnergySpectrum_at25[i]->SetYTitle("Abundance");

    std::stringstream  s9;  s9 << "HEHist_" << m_Particles[i];
    std::stringstream s10; s10 << "Hit energy curves over momentum for " << m_particles[i];
    m_HEHist[i] = new TH2D(s9.str().c_str(), s10.str().c_str(),m_nBinsMom,m_histBinsMom,m_nBinsdQdx,histBinsdQdx2);
    m_HEHist[i]->SetXTitle("Momentum (GeV)");
    m_HEHist[i]->SetYTitle("Hit Energy (GeV)");
  }

  delete[] histBinsdQdx;
  delete[] histBinsdQdx2;

  return StatusCode::SUCCESS;
}

void dQdxAnalyser::operator()(
    const edm4hep::RecDqdxCollection& dQdxs,
    const edm4hep::TrackCollection& tracks,
    const std::vector<const edm4hep::ReconstructedParticleCollection*>& pfosVec,
    const edm4hep::TrackMCParticleLinkCollection& trackTruthLinks
) const
{
  m_nEvt++;
  debug() << "dQdxAnalyser operates on event " << m_nEvt << endmsg;

  int n_track = tracks.size();
  m_NTracksHist->Fill(n_track);

  int n_pfo = 0;
  const edm4hep::ReconstructedParticleCollection* pfos{};
  if (m_usePFOTracks.value())
  {
    if (pfosVec.size()<1)
    {
      error() << "PFO collection not found - skipping event " << m_nEvt << endmsg;
      return;
    } else
    {
      pfos = pfosVec[0];
      n_pfo = pfos->size();
    }
  }
  auto trackTruthNav = podio::LinkNavigator(trackTruthLinks);

  // check how many tracks to loop over, then loop
  int n_track_pfo = m_usePFOTracks.value() ? n_pfo                   : n_track;
  int n_track_mod = m_useOneTrack.value()  ? std::min(n_track_pfo,1) : n_track_pfo;

  int n_track_used = 0;
  int n_track_mom = 0;
  std::vector<double> momVec;

  for (int i=0; i<n_track_mod; ++i)  // loop over tracks
  {
    int pdg = 0;
    double maxwei = 0;
    double maxtrackmom=0;
    int icut = 0;

    edm4hep::Track track{};
    edm4hep::ReconstructedParticle pfo{};

    if (m_usePFOTracks.value())  // get tracks from PFOs or directly
    {
      pfo = pfos->at(i);
      if (abs(pfo.getCharge())   != 1) continue;
      if (pfo.getTracks().size() != 1) continue;
      track = pfo.getTracks()[0];
    }
    else track = tracks[i];
    edm4hep::MCParticle* mcpar = NULL;
    edm4hep::utils::TrackPIDHandler dQdxHan{};
    dQdxHan.addColl(dQdxs);

    auto mcparWeiVec = trackTruthNav.getLinked(track);  // get MCPar with highest weight
    debug() << "mcparWeiVec has length of " << mcparWeiVec.size() << endmsg;
    for (unsigned int j=0; j<mcparWeiVec.size(); ++j)
    {
      float wei = mcparWeiVec[j].weight;
      if (wei>.5 && wei>maxwei)
      {
        mcpar = &mcparWeiVec[j].o;
        int mcpdg = abs(mcpar->getPDG());
        pdg = 10;
        for (unsigned int p=0; p<m_nPart; ++p) if (mcpdg == m_PDG[p]) {pdg = mcpdg; break;}
        maxwei = wei;
      }
      debug() << " entry " << j << " with PDG " << pdg << " and weight of " << wei << endmsg;
    }

    m_MCmult->Fill(mcparWeiVec.size());
    m_CutAnaHist->Fill(icut); ++icut;

    // get dQ/dx value and error
    auto dQdxVec = dQdxHan.getDqdxValues(track);
    float dQdx    = dQdxVec.size()==1 ? dQdxVec[0].getDQdx().value : 0;
    float dQdxErr = dQdxVec.size()==1 ? dQdxVec[0].getDQdx().error : 0;
    debug() << "dQdx has value of " << dQdx  << " +/- " << dQdxErr << endmsg;

    // make sure, there is a mcpar linked and a finite dQ/dx
    if (pdg==0 || dQdx<=0 || !mcpar) continue;
    m_CutAnaHist->Fill(icut); ++icut;

    int n_rel_tracks = trackTruthNav.getLinked(*mcpar).size();

    // if more other tracks are connected to the same MCParticle, skip this track
    if (m_cutTrackPurity.value()) if (n_rel_tracks>1) continue;
    m_CutAnaHist->Fill(icut); ++icut;

    // get track state, get track properties
    HelixClass trkHelixC;
    auto trackStates = track.getTrackStates();
    auto it = std::ranges::find(trackStates, 1, &edm4hep::TrackState::location);
    auto opt = it != trackStates.end() ? std::optional {*it} : std::nullopt;

    float mom=0, momT=0, lambda=0, tanL=0, d0=0, z0=0;
    if (opt)
    {
      auto state = *opt;
      trkHelixC.Initialize_Canonical(state.phi,state.D0,state.Z0,state.omega,state.tanLambda,m_Bz);
      mom  = sqrt(pow(trkHelixC.getMomentum()[0],2) + pow(trkHelixC.getMomentum()[1],2) + pow(trkHelixC.getMomentum()[2],2));
      momT = sqrt(pow(trkHelixC.getMomentum()[0],2) + pow(trkHelixC.getMomentum()[1],2));
      tanL = state.tanLambda;
      lambda = fabs(atan(tanL)) *180/M_PI;
      d0 = state.D0;
      z0 = state.Z0;
    }
    else continue;
    m_CutAnaHist->Fill(icut); ++icut;

    // momentum cuts
    if (m_cutMomMin.value()) if (mom<m_cutMomMin.value()) continue;
    if (m_cutMomMax.value()) if (mom>m_cutMomMax.value()) continue;
    m_CutAnaHist->Fill(icut); ++icut;

    // angular cuts
    if (m_cutLamMin.value()) if (lambda<m_cutLamMin.value()) continue;
    if (m_cutLamMax.value()) if (lambda>m_cutLamMax.value()) continue;
    m_CutAnaHist->Fill(icut); ++icut;

    // d0/z0 cuts
    if (m_cutD0.value()) if (fabs(d0)>m_cutD0.value()) continue;
    if (m_cutZ0.value()) if (fabs(z0)>m_cutZ0.value()) continue;
    m_CutAnaHist->Fill(icut); ++icut;

    // if more other tracks are connected to the same MCParticle and this is not the track with the highest momentum among them, skip this track
    if (m_cutTrackPurityMom.value()) if (n_rel_tracks>1)
    {
      auto trackWeiVec = trackTruthNav.getLinked(*mcpar);
      for (int j=0; j<n_rel_tracks; ++j)
      {
        auto tempTrack = trackWeiVec[j].o;
        auto temptrackStates = tempTrack.getTrackStates();
        auto tempit = std::ranges::find(temptrackStates, 1, &edm4hep::TrackState::location);
        auto tempopt = tempit != temptrackStates.end() ? std::optional {*tempit} : std::nullopt;
        if (tempopt)
        {
          auto state = *opt;
          HelixClass tempHelix;
          tempHelix.Initialize_Canonical(state.phi,state.D0,state.Z0,state.omega,state.tanLambda,m_Bz);
          double tempmom = sqrt(pow(tempHelix.getMomentum()[0],2) + pow(tempHelix.getMomentum()[1],2) + pow(tempHelix.getMomentum()[2],2));
          if (tempmom>maxtrackmom) maxtrackmom = tempmom;
        }
      }
      if (mom < maxtrackmom) continue;
    }
    m_CutAnaHist->Fill(icut); ++icut;

    // register outlier dE/dx values
    if (m_usedENdx.value() == 1)
    {
      bool cutdEdx = false;
      if (pdg== 211 &&  dQdx>(pow(log10(mom),2)+.2)*1e-6)
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 211, dE/dx = " << dQdx << ", upper limit = " << (pow(log10(mom),2)+.2)*1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(0);
        cutdEdx = true;
      }
      if (pdg== 321 && dQdx>(pow(log10(mom)*1.3,2)+.2)*1e-6)
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 321, dE/dx = " << dQdx << ", upper limit = " << (pow(log10(mom)*1.3,2)+.2)*1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(1);
        cutdEdx = true;
      }
      if (pdg== 321 && dQdx<(.28-mom/3.)*1e-6)
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 321, dE/dx = " << dQdx << ", lower limit = " << (.28-mom/3.)*1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(2);
        cutdEdx = true;
      }
      if (pdg==2212 && dQdx<(.3-mom/5.)*1e-6)
      {
        debug() << "Event " << m_nEvt << ": dE/dx does not match!  PDG = 2212, dE/dx = " << dQdx << ", lower limit = " << (.3-mom/5.)*1e-6 << endmsg;
        m_dEdxOutlierAbundance->Fill(3);
        cutdEdx = true;
      }

      // if the track dE/dx value is obviously not matching the MCParticle PDG, skip this track
      if (m_cutdEdx) if (cutdEdx) continue;
    }
    m_CutAnaHist->Fill(icut); ++icut;

    // All cuts passed!
    n_track_used++;
    if (mom > .1) n_track_mom++;
    m_TrackMomHist->Fill(mom);
    if (n_track > 1 ) m_TrackSiblingsMomHist->Fill(mom);
    momVec.push_back(mom);
    int nHits = track.getSubdetectorHitNumbers(m_subdetID*2-2);
    debug() << "gaseous tracker hit number: " << nHits << endmsg;

    // fill Bethe-Bloch histograms and data vector used in ::finalize()
    for (int j=0; j<m_nPart; j++)
    {
      if (pdg==m_PDG[j])
      {
        m_BBHist[j]->Fill(mom,dQdx);
        if (mom>.3 && mom<100) m_dQdxVec[j].emplace_back(std::array<double,6> {{mom,dQdx,dQdxErr,tanL,(double)nHits,momT}});
      }
    }
    if (pdg==10) m_BBHist[m_nPart]->Fill(mom,dQdx);
    m_BBHist[m_nPart+1]->Fill(mom,dQdx);

    if (pdg==11 && mom>.3) m_ResElectrons->Fill(dQdx);
    m_dQdxErrorCorrelation->Fill(dQdx,dQdxErr);

    int trkHitsSize = track.trackerHits_size();
    for (int h=0; h<trkHitsSize; ++h)
    {
      double EDep = track.getTrackerHits(h).getEDep();
      for (int j=0; j<m_nPart; j++) if (pdg==m_PDG[j]) m_HitEnergySpectrum[j]->Fill(EDep);
      if (pdg==10) m_HitEnergySpectrum[m_nPart]->Fill(EDep);
      m_HitEnergySpectrum[m_nPart+1]->Fill(EDep);
      if (.9<mom && mom<1.15)  // around 1 GeV
      {
        for (int j=0; j<m_nPart; j++) if (pdg==m_PDG[j]) m_HitEnergySpectrum_at1[j]->Fill(EDep);
        if (pdg==10) m_HitEnergySpectrum_at1[m_nPart]->Fill(EDep);
        m_HitEnergySpectrum_at1[m_nPart+1]->Fill(EDep);
      }
      if (2<mom && mom<3)  // around 2.5 GeV
      {
        for (int j=0; j<m_nPart; j++) if (pdg==m_PDG[j]) m_HitEnergySpectrum_at25[j]->Fill(EDep);
        if (pdg==10) m_HitEnergySpectrum_at25[m_nPart]->Fill(EDep);
        m_HitEnergySpectrum_at25[m_nPart+1]->Fill(EDep);
      }
      if (1) //(momT>=1 && nHits>=200)
      {
        for (int j=0; j<m_nPart; j++) if (pdg==m_PDG[j]) m_HEHist[j]->Fill(mom,EDep);
        if (pdg==10) m_HEHist[m_nPart]->Fill(mom,EDep);
        m_HEHist[m_nPart+1]->Fill(mom,EDep);
      }
    }

    for (int j=0; j<m_nPart; j++) if (pdg==m_PDG[j]) m_HitNumberSpectrum[j]->Fill(nHits);
    if (pdg==10) m_HitNumberSpectrum[m_nPart]->Fill(nHits);
    m_HitNumberSpectrum[m_nPart+1]->Fill(nHits);

  } // -loop over tracks

  m_NTracksUsedHist->Fill(n_track_used);
  m_NTracksMomHist ->Fill(n_track_mom);
  std::sort(momVec.begin(),momVec.end());

};

double BB_Fit( double* x, double* p ) // Bethe-Bloch fit function
{
  double bg = x[0]/p[0];
  double bb = sqrt(bg*bg/(1+bg*bg));
  return 1.35e-7 * (.5*p[1]*TMath::Log(p[2]*pow(bg,2)*pow(bg,2)) - p[3]*bb*bb - p[4]*bg/2.) / (bb*bb);
};

StatusCode dQdxAnalyser::finalize()
{
  info() << "finalising..." << endmsg;

  // set up fitting functions
  TF1* gausfit = new TF1("gausfit","gaus(0)",0,m_maxBindQdx);
  gausfit->SetParameters(1,m_maxBindQdx/3.,m_maxBindQdx/10.);
  gausfit->SetParLimits(0,0,1e6);
  gausfit->SetParLimits(1,0,m_maxBindQdx);
  gausfit->SetParLimits(2,0,m_maxBindQdx/2.);

  TF1* gausfitrel = new TF1("gausfitrel","gaus(0)",0,2);
  gausfitrel->SetParameters(m_nEvt/100.,1,.2);

  TF2* powla2 = new TF2("powla2","0.047*([0]/cos(x/180.*TMath::Pi()))**[1]*(y/[2])**[3]",m_minBinLam,m_maxBinLam,m_minBinHits,m_maxBinHits);
  powla2->SetParameters(1,-0.34,220,-.47);

  TF1* landaufit = new TF1("landaufit","landau(0)",0,m_maxBindQdx2);

  TF1* poly2fit = new TF1("poly2fit","pol2",0,90);
  TF1* poly3fit = new TF1("poly3fit","pol3",0,90);

  TF1* lamfit = new TF1("lamfit","[0]*((cos(x/180.*TMath::Pi()))**[1])",m_minBinLam,m_maxBinLam);
  lamfit->SetParameters(.05,-.34);
  TF1* nhitfit = new TF1("nHitfit","[0]*((x/220.)**[1])",m_minBinHits,m_maxBinHits);
  nhitfit->SetParameters(.05,-.47);

  TF1* angcorr = new TF1("angcorr" ,"1./(1+[0]*x*x)",0,1); // x = cos(theta)
  angcorr->SetParameter(0,-.1);
  TF1* nhitcorr= new TF1("nhitcorr","1+exp(-x/[0])",1e-6,1e6); // x = nHits
  nhitcorr->SetParameter(0,1e1);

  TF1* BBfit = new TF1("BBfit",BB_Fit,pow(10,m_minBinMom),pow(10,m_maxBinMom),5);
  BBfit->SetParameters(m_Masses[0], -1.28883368e-02, 2.72959919e+02, -1.74534200e+00, -9.84887586e-07);
  std::array<TFitResultPtr,m_nPart> ptr_BBfit;

  // fitting and extraction of results
  // fit the BB histograms vertically for each momentum bin, put fit result histograms in FitHist
  for (int i=0; i<m_nPart; i++) m_BBHist[i]->FitSlicesY(gausfit,0,-1,0,"QR",m_FitHist[i]);

  // from fit results calculate separation power for each particle combination and each momentum bin
  // as well as resolution and fill into corresponding histograms
  int sp_c = 0;
  for (int i=0; i<m_nPart; i++)
  {
    TH1* FitHist1_mean  = (TH1D*)m_FitHist[i]->At(1);
    TH1* FitHist1_sigma = (TH1D*)m_FitHist[i]->At(2);

    for (int k=1; k<m_nBinsMom+1; k++)  // Resolution from fit results; 0 is underflow bin, nBinsx+1 is overflow bin
    {
      double res = FitHist1_sigma->GetBinContent(k)*1. / FitHist1_mean->GetBinContent(k);
      if (std::isnan(res)) res = 0;
      m_FitHist_Res[i]->SetBinContent(k,res);
      if (m_histBinsMom[k] >= 1) m_FitHist_ResSum[i]->Fill(res);
    }
    m_FitHist_ResSum[i]->Fit("gaus","Q");

    BBfit->FixParameter(0,m_Masses[i]);
    ptr_BBfit[i] = FitHist1_mean->Fit(BBfit,"+S");
    debug() << "BBFit done" << endmsg;

    // normalised resolution
    for (unsigned int v=0; v<m_dQdxVec[i].size(); v++)
    {
      // m_dQdxVec[i][v] = mom,dQdx,dQdxErr,tanL,(double)nHits,momT
      double mom      = m_dQdxVec[i][v][0];
      double dQdx     = m_dQdxVec[i][v][1];
      int    bin      = FitHist1_mean->FindBin(mom);
      double binMean  = FitHist1_mean->GetBinContent(bin);
      double binError = FitHist1_mean->GetBinError(bin);
      double normMean = dQdx/binMean;
      double absError = sqrt( pow(m_dQdxVec[i][v][2]/binMean,2) + pow(dQdx/binMean/binMean*binError,2) );
      double relError = absError/normMean;
      double tanL     = m_dQdxVec[i][v][3];
      double lambda   = fabs(atan(tanL)) *180/M_PI;
      double costh    = sqrt(tanL*tanL/(1+tanL*tanL));
      int    nHits    = m_dQdxVec[i][v][4];
      double momT     = m_dQdxVec[i][v][5];
      //m_angcorr_dQdxVec[i].emplace_back(std::array<double,2> {{mom,dQdx}});

      if (binMean > 0)
      {
        if (1)                    m_ResNormHist[i]        ->Fill(normMean);
        if (momT>1)               m_ResNormHist_1GeV[i]   ->Fill(normMean);

        if (momT>1 && nHits>=200) m_NormLambdaHist[i]     ->Fill(lambda,normMean);
        if (lambda<20)            m_NormNHitHist[i]       ->Fill(nHits, normMean);
        if (1)                    m_NormNHitHist[m_nPart] ->Fill(nHits, normMean);

        if (1)                    m_NormBothHist[i]       ->Fill(lambda,nHits,normMean);
        if (momT>1)               m_NormBothHist_1GeV[i]  ->Fill(lambda,nHits,normMean);

        if (momT>1) m_NormLambdaFullHist[i]      ->Fill(lambda,normMean); // per species
        if (momT>1) m_NormLambdaFullHist[m_nPart]->Fill(lambda,normMean); // all species

        if (momT>1) m_NormCosThFullHist[i]      ->Fill(costh,normMean);
        if (momT>1) m_NormCosThFullHist[m_nPart]->Fill(costh,normMean);

        // fiducial electrons
        if (i==0 && m_fidMomMin<=mom && mom<=m_fidMomMax && m_fidLamMin<=lambda && lambda<=m_fidLamMax && nHits>=m_fidNHitsMin) m_ResNorm_FiducialElectrons->Fill(normMean);

        if (i==1 && normMean>1.5) {info() << v << " " << m_dQdxVec[i][v][0] << " " << m_dQdxVec[i][v][1] << " " << m_dQdxVec[i][v][2] << " " << m_dQdxVec[i][v][3] << " " << m_dQdxVec[i][v][4] << " " << m_dQdxVec[i][v][5] << " " << normMean << " " << relError << endmsg; }
      }
    }
    debug() << "NormRes done" << endmsg;

    m_ResNormHist[i]->Fit("gaus","Q");
    m_ResNormHist_1GeV[i]->Fit("gaus","Q");

    m_NormLambdaHist[i]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitLambda[i]);  // angular distribution
    TH1* FitLambda_mean  = (TH1D*)m_FitLambda[i]->At(1);
    TH1* FitLambda_sigma = (TH1D*)m_FitLambda[i]->At(2);
    m_FitLambda_Res[i]->Divide(FitLambda_sigma,FitLambda_mean);
    m_FitLambda_Res[i]->Fit("lamfit","+Q","",m_minBinLam,55);

    m_NormNHitHist[i]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitNHit[i]);  // number-of-hits distribution
    TH1* FitNHit_mean  = (TH1D*)m_FitNHit[i]->At(1);
    TH1* FitNHit_sigma = (TH1D*)m_FitNHit[i]->At(2);
    m_FitNHit_Res[i]->Divide(FitNHit_sigma,FitNHit_mean);
    FitNHit_mean->Fit(nhitcorr,"Q");
    m_FitNHit_Res[i]->Fit(nhitfit,"+Q","",m_maxBinHits/10.,m_maxBinHits);

    m_NormBothHist[i]->FitSlicesZ(gausfitrel,0,-1,0,-1,0,"QR");  // angular and number-of-hits distribution
    TH2* FitBoth_mean  = (TH2D*)gDirectory->Get(m_NormBothHistNames[i][1].c_str());
    TH2* FitBoth_sigma = (TH2D*)gDirectory->Get(m_NormBothHistNames[i][2].c_str());
    TH2* ProjZ = dynamic_cast<TH2*>(m_NormBothHist[i]->Project3D("yx"));

    if (FitBoth_mean && FitBoth_sigma)
    {
      for (int d=1; d<=m_nBinsLam; d++)
      {
        for (int e=1; e<=m_nBinsHits; e++)
        {
          double mu = FitBoth_mean ->GetBinContent(d,e);
          double si = FitBoth_sigma->GetBinContent(d,e);
          double dmu= FitBoth_mean ->GetBinError(d,e);
          double dsi= FitBoth_sigma->GetBinError(d,e);
          double res = 1.* si / mu;
          double reserr = sqrt( pow(1.*dsi/mu, 2) + pow(1.*si/mu/mu*dmu, 2) );
          if (std::isnan(res) || std::isnan(reserr) || ProjZ->GetBinContent(d,e)==0) {res=0; reserr=0;}

          m_FitBoth_Res[i]->SetBinContent(d,e,res);
          m_FitBoth_Res[i]->SetBinError(d,e,reserr);
        }
      }
    }
    if(!FitBoth_mean)  {info() << "NormBothHist_mean_"  << m_Particles[i] << " could not be retrieved." << endmsg;}
    if(!FitBoth_sigma) {info() << "NormBothHist_sigma_" << m_Particles[i] << " could not be retrieved." << endmsg;}
    m_FitBoth_Res[i]->Fit(powla2,"Q");

    m_NormLambdaFullHist[i]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitLambdaFull[i]);  // angular distribution in lambda, no nHits cut
    TH1* FitLambdaFull_mean  = (TH1D*)m_FitLambdaFull[i]->At(1);
    TH1* FitLambdaFull_sigma = (TH1D*)m_FitLambdaFull[i]->At(2);
    m_FitLambdaFull_Res[i]->Divide(FitLambdaFull_sigma,FitLambdaFull_mean);
    FitLambdaFull_mean->Fit(poly3fit,"Q");

    m_NormCosThFullHist[i]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitCosThFull[i]);  // angular distribution in cos(theta), no nHits cut
    TH1* FitCosThFull_mean  = (TH1D*)m_FitCosThFull[i]->At(1);
    TH1* FitCosThFull_sigma = (TH1D*)m_FitCosThFull[i]->At(2);
    m_FitCosThFull_Res[i]->Divide(FitCosThFull_sigma,FitCosThFull_mean);
    FitCosThFull_mean->Scale(1./FitCosThFull_mean->GetBinContent(1));
    FitCosThFull_mean->Fit(angcorr,"Q");

    // separation power
    for (int j=i+1; j<m_nPart; j++)
    {
      TH1* FitHist2_mean  = (TH1D*)m_FitHist[j]->At(1);
      TH1* FitHist2_sigma = (TH1D*)m_FitHist[j]->At(2);
      for (int n=1; n<m_nBinsMom+1; n++)  // 0 is underflow bin, nBinsx+1 overflow bin
      {
        // get mean and sigma as well as their errors to calculate sep power
        double mu1 = FitHist1_mean->GetBinContent(n);
        double mu2 = FitHist2_mean->GetBinContent(n);
        double s1 = FitHist1_sigma->GetBinContent(n);
        double s2 = FitHist2_sigma->GetBinContent(n);
        double mu1e = FitHist1_mean->GetBinError(n);
        double mu2e = FitHist2_mean->GetBinError(n);
        double s1e = FitHist1_sigma->GetBinError(n);
        double s2e = FitHist2_sigma->GetBinError(n);
        debug() << mu1 << " " <<  mu1e << " " << mu2 << " " << mu2e << " " << s1 << " " << s1e << " " << s2 << " " << s2e << endmsg;

        // sep power = abs(mu1 - mu2) / sqrt( .5*( sigma1^2 + sigma2^2 ) )
        // make sure denominator is not 0; error calculated following simple error propagation
        double dm = fabs( mu1 - mu2 );
        double ds = pow(s1,2) + pow(s1,2);
        double sds = sqrt(.5*ds);
        double sp = (sds==0) ? 0 : dm/sds;  // separation power
        double spe_1 = 2*pow(mu1e,2)/ds;
        double spe_2 = 2*pow(mu2e,2)/ds;
        double spe_3 = 2*pow(dm,2)*pow(s1,2)*pow(s1e,2)/pow(ds,3);
        double spe_4 = 2*pow(dm,2)*pow(s2,2)*pow(s2e,2)/pow(ds,3);
        double spe = (ds==0) ? 0 : sqrt( spe_1 + spe_2 + spe_3 + spe_4 );  // separation power error
        debug() << "  " << sp << " " <<  sqrt(spe_1) << " " << sqrt(spe_2) << " " << sqrt(spe_3) << " " << sqrt(spe_4) << " " << spe << endmsg;

        m_SPHist[sp_c]->SetBinContent(n,sp);
        m_SPHist[sp_c]->SetBinError(n,spe);
      }
      sp_c++;
    }
  }  // - for (int i=0; i<m_nPart; i++)

  TFitResultPtr ptr_fe = m_ResNorm_FiducialElectrons->Fit("gaus","QS");

  m_NormNHitHist[m_nPart]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitNHit[m_nPart]);  // number-of-hits distribution
  TH1* FitNHit_mean  = (TH1D*)m_FitNHit[m_nPart]->At(1);
  TH1* FitNHit_sigma = (TH1D*)m_FitNHit[m_nPart]->At(2);
  m_FitNHit_Res[m_nPart]->Divide(FitNHit_sigma,FitNHit_mean);

  TFitResultPtr ptr_nhit = FitNHit_mean->Fit(nhitcorr,"+QS");

  m_NormLambdaFullHist[m_nPart]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitLambdaFull[m_nPart]);  // angular distribution in lambda, no nHits cut
  TH1* FitLambdaFull_mean  = (TH1D*)m_FitLambdaFull[m_nPart]->At(1);
  TH1* FitLambdaFull_sigma = (TH1D*)m_FitLambdaFull[m_nPart]->At(2);
  m_FitLambdaFull_Res[m_nPart]->Divide(FitLambdaFull_sigma,FitLambdaFull_mean);

  TFitResultPtr ptr_lm_pol2 = FitLambdaFull_mean->Fit(poly2fit,"+QS");
  TFitResultPtr ptr_lm_pol3 = FitLambdaFull_mean->Fit(poly3fit,"+QS");

  m_NormCosThFullHist[m_nPart]->FitSlicesY(gausfitrel,0,-1,0,"QR",m_FitCosThFull[m_nPart]);  // angular distribution in cos(theta), no nHits cut
  TH1* FitCosThFull_mean  = (TH1D*)m_FitCosThFull[m_nPart]->At(1);
  TH1* FitCosThFull_sigma = (TH1D*)m_FitCosThFull[m_nPart]->At(2);
  m_FitCosThFull_Res[m_nPart]->Divide(FitCosThFull_sigma,FitCosThFull_mean);
  FitCosThFull_mean->Scale(1./FitCosThFull_mean->GetBinContent(1));
  angcorr->SetParameter(0,-.1);

  double tot_diff=0, tot_weight=0;
  for (int i_th=1; i_th<FitCosThFull_mean->GetNbinsX()+1; i_th++)
  {
    double diff = fabs(FitCosThFull_mean->GetBinContent(i_th) - 1);
    double weight = 1./FitCosThFull_mean->GetBinError(i_th);
    tot_diff += diff*weight;
    tot_weight += weight;
  }
  double CosThDeviation = tot_diff/tot_weight;

  TFitResultPtr ptr_ac  = FitCosThFull_mean->Fit(angcorr, "+QS");

  // fit the HE histograms for each beta*gamma bin, put fit result histograms in _HEFitHist
  for (int i=0; i<m_nPart; i++) m_HEHist[i]->FitSlicesY(landaufit,0,-1,0,"QR",m_HEFitHist[i]);

  // calculate p-value-based separation power, with and without angular correction
  std::vector< std::vector<double> > dQdxVec[m_nPart];
  std::vector< std::vector<double> > dQdxVec_angcorr[m_nPart];

  for (int i=0; i<m_nPart; i++)
  {
    dQdxVec[i]         = std::vector< std::vector<double> >(m_nBinsMom, std::vector<double>(0));
    dQdxVec_angcorr[i] = std::vector< std::vector<double> >(m_nBinsMom, std::vector<double>(0));
    for (unsigned int v=0; v<m_dQdxVec[i].size(); v++)
    {
      double lambda = fabs(atan(m_dQdxVec[i][v][3])) *180/M_PI;
      int    lambin = FitLambdaFull_mean->FindBin(lambda);
      double lamMean = FitLambdaFull_mean->GetBinContent(lambin);
      double ancorrMean = m_dQdxVec[i][v][1]/lamMean;

      for (int j=0; j<m_nBinsMom; ++j)
      {
        if (m_dQdxVec[i][v][0]>m_histBinsMom[j] && m_dQdxVec[i][v][0]<m_histBinsMom[j+1])
        {
          dQdxVec[i][j].emplace_back(m_dQdxVec[i][v][1]);
          dQdxVec_angcorr[i][j].emplace_back(ancorrMean);
        }
      }
    }

    for (int j=0; j<m_nBinsMom; ++j)
    {
      std::sort(dQdxVec[i][j].begin(), dQdxVec[i][j].end());
      std::sort(dQdxVec_angcorr[i][j].begin(), dQdxVec_angcorr[i][j].end());
    }
  }

  sp_c = 0;
  for (int i=0; i<m_nPart; i++)
  {
    for (int j=i+1; j<m_nPart; j++)
    {
      for (int k=0; k<m_nBinsMom; ++k)
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

  for (int i=0; i<m_nPart; i++)
  {
    m_BBHist[i]->Write();

    m_FitHist_Res[i]->Write();
    m_FitHist_ResSum[i]->Write();
    m_ResNormHist[i]->Write();
    m_ResNormHist_1GeV[i]->Write();

    m_NormLambdaHist[i]->Write();
    m_FitLambda_Res[i]->Write();
    m_NormNHitHist[i]->Write();
    m_FitNHit_Res[i]->Write();

    m_NormBothHist[i]->Write();
    m_NormBothHist_1GeV[i]->Write();
    m_FitBoth_Res[i]->Write();
    m_FitBoth_Res_1GeV[i]->Write();

    m_NormLambdaFullHist[i]->Write();
    m_FitLambdaFull_Res[i]->Write();
    m_NormCosThFullHist[i]->Write();
    m_FitCosThFull_Res[i]->Write();

    for (int j=0; j<4; j++)
    {
      m_FitHist[i]->At(j)->Write();
      m_FitLambda[i]->At(j)->Write();
      m_FitNHit[i]->At(j)->Write();
      m_FitLambdaFull[i]->At(j)->Write();
      m_FitCosThFull[i]->At(j)->Write();

      m_HEFitHist[i]->At(j)->Write();
    }
  }

  // some histograms exist also for 'other' or 'all' particles
  m_BBHist[m_nPart  ]->Write();
  m_BBHist[m_nPart+1]->Write();

  m_NormNHitHist[m_nPart]->Write();
  m_FitNHit_Res[m_nPart]->Write();
  m_NormLambdaFullHist[m_nPart]->Write();
  m_FitLambdaFull_Res[m_nPart]->Write();
  m_NormCosThFullHist[m_nPart]->Write();
  m_FitCosThFull_Res[m_nPart]->Write();

  for (int j=0; j<4; j++)
  {
    m_FitNHit[m_nPart]->At(j)->Write();
    m_FitLambdaFull[m_nPart]->At(j)->Write();
    m_FitCosThFull[m_nPart]->At(j)->Write();
  }

  // separation power
  for (int i=0; i<sp_c; i++)
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

  for (int i=0; i<m_nPart+2; i++)
  {
    m_HitEnergySpectrum[i]->Write();
    m_HitNumberSpectrum[i]->Write();
    m_HitEnergySpectrum_at1[i]->Write();
    m_HitEnergySpectrum_at25[i]->Write();
    m_HEHist[i]->Write();
  }

  // plot stuff
  if (m_plotStuff)
  {
    TCanvas* can = new TCanvas;
    gStyle->SetPalette(kBird);
    can->SetGrid();
    can->SetLogx(); // start with log x-axis plots

    std::vector<double> YMax  = {1000,m_maxBindQdx,m_maxBindQdx/10.,10};  // y-axis maxima of dE/dx fit result histogram plots
    std::vector<double> YMaxn = {1000,2,.2,10};  // y-axis maxima of normalised fit result histogram plots
    std::vector<double> YMaxHE= {1e6,m_maxBindQdx*6,m_maxBindQdx*6/10.,100};  // y-axis maxima of hit energy fit result histogram plots

    for (int i=0; i<m_nPart; i++)
    {
      PlotTH2(can, m_BBHist[i]);
      m_FitHist_Res[i]->SetAxisRange(0,.2,"Y");
      PlotTH1(can, m_FitHist_Res[i]);

      for (int j=0; j<4; j++)
      {
        TH1D* h = (TH1D*) m_FitHist[i]->At(j);
        h->SetAxisRange(0,YMax[j],"Y");
        PlotTH1_FPn(can, h, j);
      }
      if (ptr_BBfit[i].Get()) ptr_BBfit[i]->Print();
      else {info() << "Mass " << m_Masses[i] << ": fit result pointer empty!" << endmsg;}
    }
    PlotTH2(can, m_BBHist[m_nPart]);
    PlotTH2(can, m_BBHist[m_nPart+1]);

    for (int i=0; i<sp_c; i++)
    {
      m_SPHist[i]->SetAxisRange(0,20,"Y");
      PlotTH1(can, m_SPHist[i]);
      m_SPHist_pValue[i]->SetAxisRange(0,20,"Y");
      PlotTH1(can, m_SPHist_pValue[i]);
      m_SPHist_pValue_angcorr[i]->SetAxisRange(0,20,"Y");
      PlotTH1(can, m_SPHist_pValue_angcorr[i]);
    }

    for (int i=0; i<m_nPart; i++)
    {
      PlotTH2(can, m_HEHist[i]);
      for (int j=0; j<4; j++)
      {
        TH1D* h = (TH1D*) m_HEFitHist[i]->At(j);
        h->SetAxisRange(0,YMaxHE[j],"Y");
        PlotTH1_FPn(can, h, j);
      }
    }

    PlotTH2(can, m_HEHist[m_nPart]);
    PlotTH2(can, m_HEHist[m_nPart+1]);

    PlotTH1(can, m_TrackMomHist);
    PlotTH1(can, m_TrackSiblingsMomHist);

    // linear x-axis plots
    can->SetLogx(0);

    PlotTH1(can, m_ResElectrons);

    for (int i=0; i<m_nPart; i++)
    {
      PlotTH1(can, m_FitHist_ResSum[i]);

      gStyle->SetOptFit(1);
      PlotTH1(can, m_ResNormHist[i]);
      PlotTH1(can, m_ResNormHist_1GeV[i]);

      PlotTH2(can, m_NormLambdaHist[i]);
      m_FitLambda_Res[i]->SetAxisRange(0,.2,"Y");
      PlotTH1(can, m_FitLambda_Res[i]);

      PlotTH2(can, m_NormNHitHist[i]);
      m_FitNHit_Res[i]->SetAxisRange(0,.2,"Y");
      PlotTH1(can, m_FitNHit_Res[i]);
      gStyle->SetOptFit(0);

      m_FitBoth_Res[i]->SetAxisRange(0,.2,"Z");
      PlotTH2(can, m_FitBoth_Res[i]);
      m_FitBoth_Res_1GeV[i]->SetAxisRange(0,.2,"Z");
      PlotTH2(can, m_FitBoth_Res_1GeV[i]);

      PlotTH2(can, m_NormLambdaFullHist[i]);
      m_FitLambdaFull_Res[i]->SetAxisRange(0,.2,"Y");
      PlotTH1(can, m_FitLambdaFull_Res[i]);
    }
    PlotTH2(can, m_NormLambdaFullHist[m_nPart]);
    m_FitLambdaFull_Res[m_nPart]->SetAxisRange(0,.2,"Y");
    PlotTH1(can, m_FitLambdaFull_Res[m_nPart]);

    gStyle->SetOptFit(1);

    PlotTH1(can, FitNHit_mean);
    if (ptr_nhit.Get()) ptr_nhit->Print();
    else {info() << FitNHit_mean->GetName() << ": fit result pointer empty!" << endmsg;}

    PlotTH1(can, FitLambdaFull_mean);
    if (ptr_lm_pol2.Get()) ptr_lm_pol2->Print();
    else {info() << FitLambdaFull_mean->GetName() << ": poly2 fit result pointer empty!" << endmsg;}
    if (ptr_lm_pol3.Get()) ptr_lm_pol3->Print();
    else {info() << FitLambdaFull_mean->GetName() << ": poly3 fit result pointer empty!" << endmsg;}

    PlotTH1(can, FitCosThFull_mean);
    if (ptr_ac.Get()) ptr_ac->Print();
    else {info() << FitCosThFull_mean->GetName() << ": fit result pointer empty!" << endmsg;}


    PlotTH1(can, m_ResNorm_FiducialElectrons);

    gStyle->SetOptFit(0);

    PlotTH1(can, m_CutAnaHist);
    PlotTH2(can, m_dQdxErrorCorrelation);

    PlotTH1(can, m_dEdxOutlierAbundance);
    m_dEdxOutlierProportion->SetBinContent(1,m_dEdxOutlierAbundance->GetBinContent(1)/m_BBHist[2]->GetEntries());
    m_dEdxOutlierProportion->SetBinContent(2,m_dEdxOutlierAbundance->GetBinContent(2)/m_BBHist[3]->GetEntries());
    m_dEdxOutlierProportion->SetBinContent(3,m_dEdxOutlierAbundance->GetBinContent(3)/m_BBHist[3]->GetEntries());
    m_dEdxOutlierProportion->SetBinContent(4,m_dEdxOutlierAbundance->GetBinContent(4)/m_BBHist[4]->GetEntries());
    PlotTH1(can, m_dEdxOutlierProportion);

    PlotTH1(can, m_NTracksHist);
    PlotTH1(can, m_NTracksUsedHist);
    PlotTH1(can, m_NTracksMomHist);

    PlotTH1(can, m_MCmult);

    for (int i=0; i<m_nPart+2; i++)
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

  if (ptr_fe.Get()){
    info() << "Fiducial electrons dE/dx resolution: " << ptr_fe->Parameter(2)/ptr_fe->Parameter(1)*100 << " %" << endmsg;
    info() << endmsg;
    outfile << "# Fiducial electrons dE/dx resolution:\n" << ptr_fe->Parameter(2)/ptr_fe->Parameter(1)<< std::endl;
  }
  else {
    info() << "Fiducial electrons fit result pointer is empty!" << endmsg;
    info() << endmsg;
    outfile << "# Fiducial electrons fit result pointer is empty!" << std::endl;
  }

  if (ptr_lm_pol3.Get()){
    info() << "Angular dependence fit (poly3 vs. lambda) result:" << endmsg;
    info() << ptr_lm_pol3->Parameter(0) << " "
           << ptr_lm_pol3->Parameter(1) << " "
           << ptr_lm_pol3->Parameter(2) << " "
           << ptr_lm_pol3->Parameter(3) << endmsg;
    info() << "with a chi2/ndf of  "
           << ptr_lm_pol3->Chi2() << "/" << ptr_lm_pol3->Ndf() << endmsg;
    outfile << "# Angular dependence fit (poly3 vs. lambda) result:" << std::endl;
    outfile << ptr_lm_pol3->Parameter(0) << " "
            << ptr_lm_pol3->Parameter(1) << " "
            << ptr_lm_pol3->Parameter(2) << " "
            << ptr_lm_pol3->Parameter(3) << "\n# with a chi2/ndf of\n"
            << ptr_lm_pol3->Chi2() << "/" << ptr_lm_pol3->Ndf() << std::endl;
  }
  else {
    info() << "Angular dependence fit result pointer is empty!" << endmsg;
    outfile << "# Angular dependence fit result pointer is empty!" << std::endl;
  }

  info() << "Remaining deviation from flat distribution: " << CosThDeviation*100 << " %" << endmsg;
  info() << "  note: this should be below 1 %" << endmsg;
  info() << endmsg;
  outfile << "# Remaining deviation from flat distribution (note: this should be below 0.01): \n" << CosThDeviation << std::endl;

  info() << "Bethe-Bloch fit results, used as input for CPID:" << endmsg;
  info() << "  note: par2 and par3 are redundant, so par3=1 is chosen" << endmsg;
  info() << "mass      chi2/ndf    parameters 1-5" << endmsg;
  outfile << "# Bethe-Bloch fit results, used as input for CPID:" << std::endl;
  outfile << "#  note: par2 and par3 are redundant, so par3=1 is chosen" << std::endl;
  outfile << "# mass    chi2/ndf    parameters 1-5" << std::endl;
  for (int i=0; i<m_nPart; ++i)
  {
    if (ptr_BBfit[i].Get()){
      info() << m_Masses[i] << "  " << ptr_BBfit[i]->Chi2() << "/" << ptr_BBfit[i]->Ndf() << "  "
             << ptr_BBfit[i]->Parameter(1) << " "
             << ptr_BBfit[i]->Parameter(2) << " 1 "
             << ptr_BBfit[i]->Parameter(3) << " "
             << ptr_BBfit[i]->Parameter(4) << " " << endmsg;
      outfile << m_Masses[i] << "  " << ptr_BBfit[i]->Chi2() << "/" << ptr_BBfit[i]->Ndf() << "  "
              << ptr_BBfit[i]->Parameter(1) << " "
              << ptr_BBfit[i]->Parameter(2) << " 1 "
              << ptr_BBfit[i]->Parameter(3) << " "
              << ptr_BBfit[i]->Parameter(4) << " " << std::endl;
    }
    else {
      info() << m_Masses[i] << " fit result pointer is empty!" << endmsg;
      outfile  << "# " << m_Masses[i] << " fit result pointer is empty!" << std::endl;
    }
  }
  info() << endmsg;

  outfile.close();

  // some final debug output
  debug() << "CutAnaHist Entries" << endmsg;
  for (unsigned int j=0; j<10; j++) debug() << j << "  " << m_CutAnaHist->GetBinContent(j) << endmsg;

  debug() << "NTracksHist Entries" << endmsg;
  for (unsigned int j=0; j<10; j++) debug() << j << "  " << m_NTracksHist->GetBinContent(j) << endmsg;

  debug() << "NTracksMomHist Entries" << endmsg;
  for (unsigned int j=0; j<10; j++) debug() << j << "  " << m_NTracksMomHist->GetBinContent(j) << endmsg;

  debug() << "TrackMomHist Entries" << endmsg;
  for (int j=0; j<m_nBinsMom; j++) debug() << j << "  " << m_histBinsMom2[j] << "  " << m_TrackMomHist->GetBinContent(j) << endmsg;


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

double dQdxAnalyser::Calculate_pValueSP(std::vector<double> v1, std::vector<double> v2)
{
  int n1 = v1.size(), n2 = v2.size();
  if (n1==0 || n2==0) return 0;
  int c1=0, c2=0, cuml1=0, cuml2=0;
  double frac1=0, frac2=0, oldfrac1=0, oldfrac2=0;
  double p=1, Z=0;

  while (true)
  {
    if (c1==n1 && c2==n2) return 0;
    else if (c1==n1) { ++cuml2; ++c2; }
    else if (c2==n2) { ++cuml1; ++c1; }
    else if (v1[c1]<v2[c2]) { ++cuml1; ++c1; }
    else { ++cuml2; ++c2; }

    frac1 = 1.*cuml1/n1;
    frac2 = 1.*cuml2/n2;

    if (frac1<frac2 && frac1>1-frac2)
    {
      if (frac1/(1-frac2) < (1-oldfrac2)/oldfrac1) p = pow(frac1*(1-frac2),.5);
      else p = pow(oldfrac1*(1-oldfrac2),.5);
      debug() << n1 << " " << n2 << " " << c1 << " " << c2 << " " << frac1 << " " << 1-frac2 << " " << oldfrac1 << " " << 1-oldfrac2 << " " << p << " " << 2*ROOT::Math::gaussian_quantile(1-p,1) << endmsg;
      Z = 2*ROOT::Math::gaussian_quantile(1-p,1);
      if (isnan(Z) || isinf(Z)) return 0;
      else return Z;
    }
    else if (frac2<frac1 && frac2>1-frac1)
    {
      if (frac2/(1-frac1) < (1-oldfrac1)/oldfrac2) p = pow(frac2*(1-frac1),.5);
      else p = pow(oldfrac2*(1-oldfrac1),.5);
      debug() << n1 << " " << n2 << " " << c1 << " " << c2 << " " << 1-frac1 << " " << frac2 << " " << 1-oldfrac1 << " " << oldfrac2 << " " << p << " " << 2*ROOT::Math::gaussian_quantile(1-p,1) << endmsg;
      Z = 2*ROOT::Math::gaussian_quantile(1-p,1);
      if (isnan(Z) || isinf(Z)) return 0;
      else return Z;
    }

    oldfrac1 = frac1;
    oldfrac2 = frac2;
  }
}

// helper functions for histogram plotting
void dQdxAnalyser::PlotTH1( TCanvas* can, TH1* hist )
{
  if (!hist) {info() << "Histogramm is NULL!" << endmsg; return;}
  hist->Draw();
  can->Update();
  std::stringstream s; s << m_plotFolder.value() << "/" << hist->GetName() << m_plotFormat.value();
  can->Print(s.str().c_str());
};
void dQdxAnalyser::PlotTH2( TCanvas* can, TH2* hist )
{
  if (!hist) {info() << "Histogramm is NULL!" << endmsg; return;}
  hist->Draw("colz");
  can->Update();
  std::stringstream s; s << m_plotFolder.value() << "/" << hist->GetName() << m_plotFormat.value();
  can->Print(s.str().c_str());
};
void dQdxAnalyser::PlotTH1_FPn( TCanvas* can, TH1* hist, int n )
{
  if (!hist) {info() << "Histogramm is NULL!" << endmsg; return;}
  hist->Draw();
  can->Update();
  std::stringstream s; s << m_plotFolder.value() << "/" << hist->GetName() << "_FitPar" << n <<  m_plotFormat.value();
  can->Print(s.str().c_str());
};
