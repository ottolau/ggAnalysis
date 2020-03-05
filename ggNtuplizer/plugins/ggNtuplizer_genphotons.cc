#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          gen_nPho_;
vector<float>  gen_phoE_;
vector<float>  gen_phoEt_;
vector<float>  gen_phoEta_;
vector<float>  gen_phoPhi_;
vector<float>  gen_phoCalibE_;
vector<float>  gen_phoCalibEt_;
vector<float>  gen_phoSCE_;
vector<float>  gen_phoSCRawE_;
vector<float>  gen_phoESEnP1_;
vector<float>  gen_phoESEnP2_;
vector<float>  gen_phoSCEta_;
vector<float>  gen_phoSCPhi_;
vector<float>  gen_phoSCEtaWidth_;
vector<float>  gen_phoSCPhiWidth_;
vector<float>  gen_phoSCBrem_;
vector<int>    gen_phohasPixelSeed_;
vector<int>    gen_phoEleVeto_;
vector<float>  gen_phoR9_;
vector<float>  gen_phoHoverE_;
//vector<float>  gen_phoSigmaIEtaIEta_;
//vector<float>  gen_phoSigmaIEtaIPhi_;
//vector<float>  gen_phoSigmaIPhiIPhi_;
vector<float>  gen_phoE1x3_;
vector<float>  gen_phoE1x5_;
vector<float>  gen_phoE2x2_;
vector<float>  gen_phoE2x5Max_;
vector<float>  gen_phoE5x5_;
vector<float>  gen_phoESEffSigmaRR_;
vector<float>  gen_phoSigmaIEtaIEtaFull5x5_;
vector<float>  gen_phoSigmaIEtaIPhiFull5x5_;
vector<float>  gen_phoSigmaIPhiIPhiFull5x5_;
vector<float>  gen_phoE1x3Full5x5_;
vector<float>  gen_phoE1x5Full5x5_;
vector<float>  gen_phoE2x2Full5x5_;
vector<float>  gen_phoE2x5MaxFull5x5_;
vector<float>  gen_phoE5x5Full5x5_;
vector<float>  gen_phoR9Full5x5_;
vector<float>  gen_phoPFChIso_;
vector<float>  gen_phoPFPhoIso_;
vector<float>  gen_phoPFNeuIso_;
vector<float>  gen_phoPFChWorstIso_;
vector<float>  gen_phoPFChIsoFrix1_;
vector<float>  gen_phoPFChIsoFrix2_;
vector<float>  gen_phoPFChIsoFrix3_;
vector<float>  gen_phoPFChIsoFrix4_;
vector<float>  gen_phoPFChIsoFrix5_;
vector<float>  gen_phoPFChIsoFrix6_;
vector<float>  gen_phoPFChIsoFrix7_;
vector<float>  gen_phoPFChIsoFrix8_;
vector<float>  gen_phoPFPhoIsoFrix1_;
vector<float>  gen_phoPFPhoIsoFrix2_;
vector<float>  gen_phoPFPhoIsoFrix3_;
vector<float>  gen_phoPFPhoIsoFrix4_;
vector<float>  gen_phoPFPhoIsoFrix5_;
vector<float>  gen_phoPFPhoIsoFrix6_;
vector<float>  gen_phoPFPhoIsoFrix7_;
vector<float>  gen_phoPFPhoIsoFrix8_;
vector<float>  gen_phoPFNeuIsoFrix1_;
vector<float>  gen_phoPFNeuIsoFrix2_;
vector<float>  gen_phoPFNeuIsoFrix3_;
vector<float>  gen_phoPFNeuIsoFrix4_;
vector<float>  gen_phoPFNeuIsoFrix5_;
vector<float>  gen_phoPFNeuIsoFrix6_;
vector<float>  gen_phoPFNeuIsoFrix7_;
vector<float>  gen_phoPFNeuIsoFrix8_;
vector<float>  gen_phoCITKChIso_;
vector<float>  gen_phoCITKPhoIso_;
vector<float>  gen_phoCITKNeuIso_;
//vector<float>  gen_phoSeedBCE_;
//vector<float>  gen_phoSeedBCEta_;
vector<float>  gen_phoIDMVA_;
vector<ULong64_t> gen_phoFiredSingleTrgs_;
vector<ULong64_t> gen_phoFiredDoubleTrgs_;
vector<ULong64_t> gen_phoFiredTripleTrgs_;
vector<ULong64_t> gen_phoFiredL1Trgs_;
//vector<float>  gen_phoEcalRecHitSumEtConeDR03_;
//vector<float>  gen_phohcalDepth1TowerSumEtConeDR03_;
//vector<float>  gen_phohcalDepth2TowerSumEtConeDR03_;
//vector<float>  gen_phohcalTowerSumEtConeDR03_;
//vector<float>  gen_photrkSumPtHollowConeDR03_;
//vector<float>  gen_photrkSumPtSolidConeDR03_;
vector<float>  gen_phoSeedTime_;
vector<float>  gen_phoSeedEnergy_;
//vector<float>  gen_phoSeedTimeFull5x5_;
//vector<float>  gen_phoMIPChi2_;
//vector<float>  gen_phoMIPTotEnergy_;
//vector<float>  gen_phoMIPSlope_;
//vector<float>  gen_phoMIPIntercept_;
//vector<float>  gen_phoMIPNhitCone_;
//vector<float>  gen_phoMIPIsHalo_;
vector<UShort_t> gen_phoxtalBits_;
vector<UShort_t> gen_phoIDbit_;
vector<float>    gen_phoScale_stat_up_;
vector<float>    gen_phoScale_stat_dn_;
vector<float>    gen_phoScale_syst_up_;
vector<float>    gen_phoScale_syst_dn_;
vector<float>    gen_phoScale_gain_up_;
vector<float>    gen_phoScale_gain_dn_;
vector<float>    gen_phoResol_rho_up_;
vector<float>    gen_phoResol_rho_dn_;
vector<float>    gen_phoResol_phi_up_;
vector<float>    gen_phoResol_phi_dn_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesGenPhotons(TTree* tree) {
  
  tree->Branch("gen_nPho",                    &gen_nPho_);
  tree->Branch("gen_phoE",                    &gen_phoE_);
  tree->Branch("gen_phoEt",                   &gen_phoEt_);
  tree->Branch("gen_phoEta",                  &gen_phoEta_);
  tree->Branch("gen_phoPhi",                  &gen_phoPhi_);
  tree->Branch("gen_phoCalibE",               &gen_phoCalibE_);
  tree->Branch("gen_phoCalibEt",              &gen_phoCalibEt_);
  tree->Branch("gen_phoSCE",                  &gen_phoSCE_);
  tree->Branch("gen_phoSCRawE",               &gen_phoSCRawE_);
  tree->Branch("gen_phoESEnP1",               &gen_phoESEnP1_);
  tree->Branch("gen_phoESEnP2",               &gen_phoESEnP2_);
  tree->Branch("gen_phoSCEta",                &gen_phoSCEta_);
  tree->Branch("gen_phoSCPhi",                &gen_phoSCPhi_);
  tree->Branch("gen_phoSCEtaWidth",           &gen_phoSCEtaWidth_);
  tree->Branch("gen_phoSCPhiWidth",           &gen_phoSCPhiWidth_);
  tree->Branch("gen_phoSCBrem",               &gen_phoSCBrem_);
  tree->Branch("gen_phohasPixelSeed",         &gen_phohasPixelSeed_);
  tree->Branch("gen_phoEleVeto",              &gen_phoEleVeto_);
  tree->Branch("gen_phoR9",                   &gen_phoR9_);
  tree->Branch("gen_phoHoverE",               &gen_phoHoverE_);
  //tree->Branch("gen_phoSigmaIEtaIEta",        &gen_phoSigmaIEtaIEta_);
  //tree->Branch("gen_phoSigmaIEtaIPhi",        &gen_phoSigmaIEtaIPhi_);
  //tree->Branch("gen_phoSigmaIPhiIPhi",        &gen_phoSigmaIPhiIPhi_);
  tree->Branch("gen_phoE1x3",                 &gen_phoE1x3_);
  tree->Branch("gen_phoE1x5",                 &gen_phoE1x5_);
  tree->Branch("gen_phoE2x2",                 &gen_phoE2x2_);
  tree->Branch("gen_phoE2x5Max",              &gen_phoE2x5Max_);
  tree->Branch("gen_phoE5x5",                 &gen_phoE5x5_);
  tree->Branch("gen_phoESEffSigmaRR",         &gen_phoESEffSigmaRR_);
  tree->Branch("gen_phoSigmaIEtaIEtaFull5x5", &gen_phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("gen_phoSigmaIEtaIPhiFull5x5", &gen_phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("gen_phoSigmaIPhiIPhiFull5x5", &gen_phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("gen_phoE1x3Full5x5",          &gen_phoE1x3Full5x5_);
  tree->Branch("gen_phoE1x5Full5x5",          &gen_phoE1x5Full5x5_);
  tree->Branch("gen_phoE2x2Full5x5",          &gen_phoE2x2Full5x5_);
  tree->Branch("gen_phoE2x5MaxFull5x5",       &gen_phoE2x5MaxFull5x5_);
  tree->Branch("gen_phoE5x5Full5x5",          &gen_phoE5x5Full5x5_);
  tree->Branch("gen_phoR9Full5x5",            &gen_phoR9Full5x5_);
  //tree->Branch("gen_phoSeedBCE",              &gen_phoSeedBCE_);
  //tree->Branch("gen_phoSeedBCEta",            &gen_phoSeedBCEta_);
  tree->Branch("gen_phoPFChIso",              &gen_phoPFChIso_);
  tree->Branch("gen_phoPFPhoIso",             &gen_phoPFPhoIso_);
  tree->Branch("gen_phoPFNeuIso",             &gen_phoPFNeuIso_);
  tree->Branch("gen_phoPFChWorstIso",         &gen_phoPFChWorstIso_);
  /*
  tree->Branch("gen_phoPFChIsoFrix1",         &gen_phoPFChIsoFrix1_);
  tree->Branch("gen_phoPFChIsoFrix2",         &gen_phoPFChIsoFrix2_);
  tree->Branch("gen_phoPFChIsoFrix3",         &gen_phoPFChIsoFrix3_);
  tree->Branch("gen_phoPFChIsoFrix4",         &gen_phoPFChIsoFrix4_);
  tree->Branch("gen_phoPFChIsoFrix5",         &gen_phoPFChIsoFrix5_);
  tree->Branch("gen_phoPFChIsoFrix6",         &gen_phoPFChIsoFrix6_);
  tree->Branch("gen_phoPFChIsoFrix7",         &gen_phoPFChIsoFrix7_);
  tree->Branch("gen_phoPFChIsoFrix8",         &gen_phoPFChIsoFrix8_);
  tree->Branch("gen_phoPFPhoIsoFrix1",        &gen_phoPFPhoIsoFrix1_);
  tree->Branch("gen_phoPFPhoIsoFrix2",        &gen_phoPFPhoIsoFrix2_);
  tree->Branch("gen_phoPFPhoIsoFrix3",        &gen_phoPFPhoIsoFrix3_);
  tree->Branch("gen_phoPFPhoIsoFrix4",        &gen_phoPFPhoIsoFrix4_);
  tree->Branch("gen_phoPFPhoIsoFrix5",        &gen_phoPFPhoIsoFrix5_);
  tree->Branch("gen_phoPFPhoIsoFrix6",        &gen_phoPFPhoIsoFrix6_);
  tree->Branch("gen_phoPFPhoIsoFrix7",        &gen_phoPFPhoIsoFrix7_);
  tree->Branch("gen_phoPFPhoIsoFrix8",        &gen_phoPFPhoIsoFrix8_);
  tree->Branch("gen_phoPFNeuIsoFrix1",        &gen_phoPFNeuIsoFrix1_);
  tree->Branch("gen_phoPFNeuIsoFrix2",        &gen_phoPFNeuIsoFrix2_);
  tree->Branch("gen_phoPFNeuIsoFrix3",        &gen_phoPFNeuIsoFrix3_);
  tree->Branch("gen_phoPFNeuIsoFrix4",        &gen_phoPFNeuIsoFrix4_);
  tree->Branch("gen_phoPFNeuIsoFrix5",        &gen_phoPFNeuIsoFrix5_);
  tree->Branch("gen_phoPFNeuIsoFrix6",        &gen_phoPFNeuIsoFrix6_);
  tree->Branch("gen_phoPFNeuIsoFrix7",        &gen_phoPFNeuIsoFrix7_);
  tree->Branch("gen_phoPFNeuIsoFrix8",        &gen_phoPFNeuIsoFrix8_);
  */
  tree->Branch("gen_phoCITKChIso",            &gen_phoCITKChIso_);
  tree->Branch("gen_phoCITKPhoIso",           &gen_phoCITKPhoIso_);
  tree->Branch("gen_phoCITKNeuIso",           &gen_phoCITKNeuIso_);
  //tree->Branch("gen_phoEcalRecHitSumEtConeDR03",      &gen_phoEcalRecHitSumEtConeDR03_);
  //tree->Branch("gen_phohcalDepth1TowerSumEtConeDR03", &gen_phohcalDepth1TowerSumEtConeDR03_);
  //tree->Branch("gen_phohcalDepth2TowerSumEtConeDR03", &gen_phohcalDepth2TowerSumEtConeDR03_);
  //tree->Branch("gen_phohcalTowerSumEtConeDR03",       &gen_phohcalTowerSumEtConeDR03_);
  //tree->Branch("gen_photrkSumPtHollowConeDR03",       &gen_photrkSumPtHollowConeDR03_);
  //tree->Branch("gen_photrkSumPtSolidConeDR03",        &gen_photrkSumPtSolidConeDR03_);
  tree->Branch("gen_phoIDMVA",                        &gen_phoIDMVA_);
  tree->Branch("gen_phoFiredSingleTrgs",              &gen_phoFiredSingleTrgs_);
  tree->Branch("gen_phoFiredDoubleTrgs",              &gen_phoFiredDoubleTrgs_);
  tree->Branch("gen_phoFiredTripleTrgs",              &gen_phoFiredTripleTrgs_);
  tree->Branch("gen_phoFiredL1Trgs",                  &gen_phoFiredL1Trgs_);
  tree->Branch("gen_phoSeedTime",                     &gen_phoSeedTime_);
  tree->Branch("gen_phoSeedEnergy",                   &gen_phoSeedEnergy_);
  //tree->Branch("gen_phoSeedTimeFull5x5",              &gen_phoSeedTimeFull5x5_);
  //tree->Branch("gen_phoMIPChi2",                      &gen_phoMIPChi2_);
  //tree->Branch("gen_phoMIPTotEnergy",                 &gen_phoMIPTotEnergy_);
  //tree->Branch("gen_phoMIPSlope",                     &gen_phoMIPSlope_);
  //tree->Branch("gen_phoMIPIntercept",                 &gen_phoMIPIntercept_);
  //tree->Branch("gen_phoMIPNhitCone",                  &gen_phoMIPNhitCone_);
  //tree->Branch("gen_phoMIPIsHalo",                    &gen_phoMIPIsHalo_);

  tree->Branch("gen_phoxtalBits",      &gen_phoxtalBits_);
  tree->Branch("gen_phoIDbit",         &gen_phoIDbit_);
  tree->Branch("gen_phoScale_stat_up", &gen_phoScale_stat_up_);
  tree->Branch("gen_phoScale_stat_dn", &gen_phoScale_stat_dn_);
  tree->Branch("gen_phoScale_syst_up", &gen_phoScale_syst_up_);
  tree->Branch("gen_phoScale_syst_dn", &gen_phoScale_syst_dn_);
  tree->Branch("gen_phoScale_gain_up", &gen_phoScale_gain_up_);
  tree->Branch("gen_phoScale_gain_dn", &gen_phoScale_gain_dn_);
  tree->Branch("gen_phoResol_rho_up",  &gen_phoResol_rho_up_);
  tree->Branch("gen_phoResol_rho_dn",  &gen_phoResol_rho_dn_);
  tree->Branch("gen_phoResol_phi_up",  &gen_phoResol_phi_up_);
  tree->Branch("gen_phoResol_phi_dn",  &gen_phoResol_phi_dn_);

}

void ggNtuplizer::fillGenPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  gen_phoE_                 .clear();
  gen_phoEt_                .clear();
  gen_phoEta_               .clear();
  gen_phoPhi_               .clear();
  gen_phoCalibE_            .clear();
  gen_phoCalibEt_           .clear();
  gen_phoSCE_               .clear();
  gen_phoSCRawE_            .clear();
  gen_phoESEnP1_            .clear();
  gen_phoESEnP2_            .clear();
  gen_phoSCEta_             .clear();
  gen_phoSCPhi_             .clear();
  gen_phoSCEtaWidth_        .clear();
  gen_phoSCPhiWidth_        .clear();
  gen_phoSCBrem_            .clear();
  gen_phohasPixelSeed_      .clear();
  gen_phoEleVeto_           .clear();
  gen_phoR9_                .clear();
  gen_phoHoverE_            .clear();
  //gen_phoSigmaIEtaIEta_     .clear();
  //gen_phoSigmaIEtaIPhi_     .clear();
  //gen_phoSigmaIPhiIPhi_     .clear();
  gen_phoE1x3_              .clear();
  gen_phoE1x5_              .clear();
  gen_phoE2x2_              .clear();
  gen_phoE2x5Max_           .clear();
  gen_phoE5x5_              .clear();
  gen_phoESEffSigmaRR_      .clear();
  gen_phoSigmaIEtaIEtaFull5x5_.clear();
  gen_phoSigmaIEtaIPhiFull5x5_.clear();
  gen_phoSigmaIPhiIPhiFull5x5_.clear();
  gen_phoE1x3Full5x5_       .clear();
  gen_phoE1x5Full5x5_       .clear();
  gen_phoE2x2Full5x5_       .clear();
  gen_phoE2x5MaxFull5x5_    .clear();
  gen_phoE5x5Full5x5_       .clear();
  gen_phoR9Full5x5_         .clear();
  gen_phoPFChIso_           .clear();
  gen_phoPFPhoIso_          .clear();
  gen_phoPFNeuIso_          .clear();
  gen_phoPFChWorstIso_      .clear();
  gen_phoPFChIsoFrix1_      .clear();
  gen_phoPFChIsoFrix2_      .clear();
  gen_phoPFChIsoFrix3_      .clear();
  gen_phoPFChIsoFrix4_      .clear();
  gen_phoPFChIsoFrix5_      .clear();
  gen_phoPFChIsoFrix6_      .clear();
  gen_phoPFChIsoFrix7_      .clear();
  gen_phoPFChIsoFrix8_      .clear();
  gen_phoPFPhoIsoFrix1_     .clear();
  gen_phoPFPhoIsoFrix2_     .clear();
  gen_phoPFPhoIsoFrix3_     .clear();
  gen_phoPFPhoIsoFrix4_     .clear();
  gen_phoPFPhoIsoFrix5_     .clear();
  gen_phoPFPhoIsoFrix6_     .clear();
  gen_phoPFPhoIsoFrix7_     .clear();
  gen_phoPFPhoIsoFrix8_     .clear();
  gen_phoPFNeuIsoFrix1_     .clear();
  gen_phoPFNeuIsoFrix2_     .clear();
  gen_phoPFNeuIsoFrix3_     .clear();
  gen_phoPFNeuIsoFrix4_     .clear();
  gen_phoPFNeuIsoFrix5_     .clear();
  gen_phoPFNeuIsoFrix6_     .clear();
  gen_phoPFNeuIsoFrix7_     .clear();
  gen_phoPFNeuIsoFrix8_     .clear();
  gen_phoCITKChIso_         .clear();
  gen_phoCITKPhoIso_        .clear();
  gen_phoCITKNeuIso_        .clear();
  //gen_phoSeedBCE_           .clear();
  //gen_phoSeedBCEta_         .clear();
  gen_phoIDMVA_             .clear();
  gen_phoFiredSingleTrgs_   .clear();
  gen_phoFiredDoubleTrgs_   .clear();
  gen_phoFiredTripleTrgs_   .clear();
  gen_phoFiredL1Trgs_       .clear();
  gen_phoxtalBits_          .clear();
  gen_phoSeedTime_          .clear();
  gen_phoSeedEnergy_        .clear();
  /*
  gen_phoSeedTimeFull5x5_   .clear();
  gen_phoMIPChi2_           .clear();
  gen_phoMIPTotEnergy_      .clear();
  gen_phoMIPSlope_          .clear();
  gen_phoMIPIntercept_      .clear();
  gen_phoMIPNhitCone_       .clear();
  gen_phoMIPIsHalo_         .clear();
  gen_phoEcalRecHitSumEtConeDR03_     .clear();
  gen_phohcalDepth1TowerSumEtConeDR03_.clear();
  gen_phohcalDepth2TowerSumEtConeDR03_.clear();
  gen_phohcalTowerSumEtConeDR03_      .clear();
  gen_photrkSumPtHollowConeDR03_      .clear();
  gen_photrkSumPtSolidConeDR03_       .clear();
  */

  gen_phoIDbit_        .clear();
  gen_phoScale_stat_up_.clear();
  gen_phoScale_stat_dn_.clear();
  gen_phoScale_syst_up_.clear();
  gen_phoScale_syst_dn_.clear();
  gen_phoScale_gain_up_.clear();
  gen_phoScale_gain_dn_.clear();
  gen_phoResol_rho_up_ .clear();
  gen_phoResol_rho_dn_ .clear();
  gen_phoResol_phi_up_ .clear();
  gen_phoResol_phi_dn_ .clear();  

  gen_nPho_ = 0;

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  edm::Handle<edm::View<pat::Photon> > calibphotonHandle;
  e.getByToken(calibphotonCollection_, calibphotonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }

  if (!calibphotonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no calibrated pat::Photons in event";
    return;
  }

  edm::Handle<reco::PhotonCollection> recoPhotonHandle;
  e.getByToken(recophotonCollection_, recoPhotonHandle);

  ///Photon ID in VID framwork - 11th may, 2015
  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;
  edm::Handle<std::vector<pat::PackedCandidate>> pfCndHandle; 
  edm::Handle<edm::View<reco::Candidate> > CndHandle; 
  
  e.getByToken(phoLooseIdMapToken_ ,  loose_id_decisions);
  e.getByToken(phoMediumIdMapToken_,  medium_id_decisions);
  e.getByToken(phoTightIdMapToken_ ,  tight_id_decisions);
  e.getByToken(phoMVAValuesMapToken_, mvaValues);
  
  e.getByToken(phoChargedIsolationToken_,       phoChargedIsolationMap);
  e.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  e.getByToken(phoPhotonIsolationToken_,        phoPhotonIsolationMap);
  e.getByToken(phoWorstChargedIsolationToken_,  phoWorstChargedIsolationMap);

  // Get the isolation maps for CITK
  //edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap_CITK;
  //e.getByToken(phoChargedIsolationToken_CITK_, phoChargedIsolationMap_CITK);
  //edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap_CITK;
  //e.getByToken(phoPhotonIsolationToken_CITK_, phoPhotonIsolationMap_CITK);
  //edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap_CITK;
  //e.getByToken(phoNeutralHadronIsolationToken_CITK_, phoNeutralHadronIsolationMap_CITK);

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  GEDPhoIDTools* GEDIdTool = NULL;
  if (isAOD_)
    GEDIdTool = new GEDPhoIDTools(e);

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  edm::Handle<reco::VertexCollection> recVtxsBS;
  e.getByToken(vtxBSLabel_, recVtxsBS);

  edm::Handle<reco::TrackCollection> tracksHandle;
  e.getByToken(tracklabel_, tracksHandle);

  edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
  e.getByToken(gsfElectronlabel_, gsfElectronHandle);

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  double rho    = *(rhoHandle.product());

  vector<int> pho_gen_ind;

  for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    // get Bs meson (pdgid=531)
    if (fabs(ip->pdgId()) != 13) continue;
    const Candidate* muon = &(*ip);

    // loop over all generated final state particles, get 
    //for (edm::View<pat::PackedGenParticle>::const_iterator ipp = packedGenParticlesHandle->begin(); ipp != packedGenParticlesHandle->end(); ++ipp) {
    for (vector<reco::GenParticle>::const_iterator ipp = genParticlesHandle->begin(); ipp != genParticlesHandle->end(); ++ipp) {
      // check if the generated final state particle comes from the Bs
      if (ipp->status() != 1) continue;
      if (! (ipp->mother(0) && isAncestor(muon, ipp->mother(0)))) continue;
      if (ipp->pdgId() == 22) {
        pho_gen_ind.push_back(ipp - genParticlesHandle->begin());
      }
    }
  }

  if (pho_gen_ind.size() == 0) return;

  if (isAOD_) {
    edm::Handle<reco::PFCandidateCollection> pfAllCandidates;
    e.getByToken(pfAllParticles_, pfAllCandidates);

    //cicPhotonId_->configure(recVtxsBS, tracksHandle, gsfElectronHandle, pfAllCandidates, rho);
    cicPhotonId_->configure(recVtxs, tracksHandle, gsfElectronHandle, pfAllCandidates, rho);
  }

  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {
    double temp_dR = 9999999999.0, bestdR = 9999999999.0;
    for (int i = 0; i < (int) pho_gen_ind.size(); ++i){
      temp_dR = deltaR(iPho->eta(), iPho->phi(), genParticlesHandle->at(pho_gen_ind.at(i)).eta(), genParticlesHandle->at(pho_gen_ind.at(i)).phi());
      if (temp_dR < bestdR) {
        bestdR = temp_dR;
      }
    }

    if (bestdR > 0.1) continue;
  
    Float_t corrPt = -1;
    Float_t corrEn = -1;
    for (edm::View<pat::Photon>::const_iterator iCPho = calibphotonHandle->begin(); iCPho != calibphotonHandle->end(); ++iCPho) {
      if (fabs(iPho->eta() - iCPho->eta()) < 0.001 && fabs(iPho->phi() - iCPho->phi()) < 0.001) {
	corrPt = iCPho->pt();
	corrEn = iCPho->energy();
      }
    }
    gen_phoCalibEt_       .push_back(corrPt);
    gen_phoCalibE_        .push_back(corrEn);

    gen_phoE_             .push_back(iPho->energy());
    gen_phoEt_            .push_back(iPho->et());
    gen_phoEta_           .push_back(iPho->eta());
    gen_phoPhi_           .push_back(iPho->phi());
    gen_phoSCE_           .push_back((*iPho).superCluster()->energy());
    gen_phoSCRawE_        .push_back((*iPho).superCluster()->rawEnergy());
    gen_phoESEnP1_        .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
    gen_phoESEnP2_        .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
    gen_phoSCEta_         .push_back((*iPho).superCluster()->eta());
    gen_phoSCPhi_         .push_back((*iPho).superCluster()->phi());
    gen_phoSCEtaWidth_    .push_back((*iPho).superCluster()->etaWidth());
    gen_phoSCPhiWidth_    .push_back((*iPho).superCluster()->phiWidth());
    gen_phoSCBrem_        .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
    gen_phohasPixelSeed_  .push_back((Int_t)iPho->hasPixelSeed());
    gen_phoEleVeto_       .push_back((Int_t)iPho->passElectronVeto());
    gen_phoR9_            .push_back(iPho->r9());
    gen_phoHoverE_        .push_back(iPho->hadTowOverEm());
    //gen_phoSigmaIEtaIEta_ .push_back(iPho->see());
    //gen_phoSigmaIEtaIPhi_ .push_back(iPho->sep());
    //gen_phoSigmaIPhiIPhi_ .push_back(iPho->spp());
    gen_phoE1x3_          .push_back(lazyTool.e1x3(*((*iPho).superCluster()->seed())));
    gen_phoE1x5_          .push_back(iPho->e1x5());
    gen_phoE2x2_          .push_back(lazyTool.e2x2(*((*iPho).superCluster()->seed())));
    gen_phoE2x5Max_       .push_back(iPho->e2x5());
    gen_phoE5x5_          .push_back(iPho->e5x5());
    gen_phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
    //gen_phoPFChIso_       .push_back(iPho->chargedHadronIso());
    //gen_phoPFPhoIso_      .push_back(iPho->photonIso());
    //gen_phoPFNeuIso_      .push_back(iPho->neutralHadronIso());


    ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
    DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
            
    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      //std::cout<<"(*theSeedHit).time()"<<(*theSeedHit).time()<<"seed energy: "<<(*theSeedHit).energy()<<std::endl;  

      gen_phoSeedTime_  .push_back((*theSeedHit).time());
      gen_phoSeedEnergy_.push_back((*theSeedHit).energy());
    } else{
      gen_phoSeedTime_  .push_back(-99.);
      gen_phoSeedEnergy_.push_back(-99.);
    }
    /*
    // ECAL scale correction and smearing
    float et = iPho->et();
    unsigned int gainSeedSC = 12;
    if (theSeedHit != rechits->end()) { 
      if(theSeedHit->checkFlag(EcalRecHit::kHasSwitchToGain6)) gainSeedSC = 6;
      if(theSeedHit->checkFlag(EcalRecHit::kHasSwitchToGain1)) gainSeedSC = 1;
    }
    int runNumber = e.id().run();
    double scale = egmScaler_->ScaleCorrection(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC);  
    double smear = egmScaler_->getSmearingSigma(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 0., 0.);
  
    float scale_stat_up = scale + egmScaler_->ScaleCorrectionUncertainty(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 1);
    float scale_stat_dn = scale - egmScaler_->ScaleCorrectionUncertainty(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 1);
    float scale_syst_up = scale + egmScaler_->ScaleCorrectionUncertainty(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 2);
    float scale_syst_dn = scale - egmScaler_->ScaleCorrectionUncertainty(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 2);
    float scale_gain_up = scale + egmScaler_->ScaleCorrectionUncertainty(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 4);
    float scale_gain_dn = scale - egmScaler_->ScaleCorrectionUncertainty(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 4);
    float resol_rho_up  = egmScaler_->getSmearingSigma(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 1., 0.);
    float resol_rho_dn  = egmScaler_->getSmearingSigma(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, -1., 0.);
    float resol_phi_up  = egmScaler_->getSmearingSigma(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 0., 1.);
    float resol_phi_dn = egmScaler_->getSmearingSigma(runNumber, iPho->isEB(), iPho->full5x5_r9(), fabs(iPho->eta()), et, gainSeedSC, 0., -1.);
    
    phoScale_stat_up_.push_back(scale_stat_up);
    phoScale_stat_dn_.push_back(scale_stat_dn);
    phoScale_syst_up_.push_back(scale_syst_up);
    phoScale_syst_dn_.push_back(scale_syst_dn);
    phoScale_gain_up_.push_back(scale_gain_up);
    phoScale_gain_dn_.push_back(scale_gain_dn);
    phoResol_rho_up_.push_back(resol_rho_up);
    phoResol_rho_dn_.push_back(resol_rho_dn);
    phoResol_phi_up_.push_back(resol_phi_up);
    phoResol_phi_dn_.push_back(resol_phi_dn);
    /////////////////////////////////END of energy and scale systematics
    */
    /// if( isBarrel ) {
    ///     EBDetId ebId(seed);
    ///     cout << "seed barrel " << ebId.ieta() << " " << ebId.iphi() << endl;
    /// } else {
    ///     EEDetId eeId(seed);
    ///     cout << "seed endpcas " << eeId.ix() << " " << eeId.iy() << endl;
    /// 
    /// }
    unsigned short nSaturated = 0, nLeRecovered = 0, nNeighRecovered = 0, nGain1 = 0, nGain6 = 0, nWeired = 0;

    int isSaturated = 0;
    int isSaturated_gain6 = 0;
    
    UShort_t tmpxtalbit = 0;

    auto matrix5x5 = lazyTool.matrixDetId(seed,-2,+2,-2,+2);
    for (auto & deId : matrix5x5 ) {
      /// cout << "matrix " << deId.rawId() << endl;
      auto rh = rechits->find(deId);
      if( rh != rechits->end() ) {
	nSaturated += rh->checkFlag( EcalRecHit::kSaturated );
	nLeRecovered += rh->checkFlag( EcalRecHit::kLeadingEdgeRecovered );
	nNeighRecovered += rh->checkFlag( EcalRecHit::kNeighboursRecovered );
	nGain1 += rh->checkFlag( EcalRecHit::kHasSwitchToGain1 );
	nGain6 += rh->checkFlag( EcalRecHit::kHasSwitchToGain6 );
	nWeired += rh->checkFlag( EcalRecHit::kWeird ) || rh->checkFlag( EcalRecHit::kDiWeird );
	
	if( rh->checkFlag( EcalRecHit::kHasSwitchToGain1 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

	  setbit(tmpxtalbit, 0);
	  isSaturated = 1;
	  //break;
	}
	
	if( rh->checkFlag( EcalRecHit::kHasSwitchToGain6 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated_gain6){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

	  setbit(tmpxtalbit, 1);
	  isSaturated_gain6 = 1;
	  //break;
	}
	
      }//if( rh != rechits->end() ) 
       
      if (nWeired>0) setbit(tmpxtalbit,2);      
      if (nGain6>0) setbit(tmpxtalbit,3); 

    }//for(auto & deId : matrix5x5 )
  
    gen_phoxtalBits_.push_back(tmpxtalbit);

    gen_phoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    gen_phoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    gen_phoFiredTripleTrgs_     .push_back(matchTriplePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    gen_phoFiredL1Trgs_         .push_back(matchL1TriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

    std::vector<float> vCov = lazyToolnoZS.localCovariances( *((*iPho).superCluster()->seed()) );
    //const float see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    gen_phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    gen_phoSigmaIEtaIPhiFull5x5_ .push_back(sep);
    gen_phoSigmaIPhiIPhiFull5x5_ .push_back(spp);
    gen_phoE1x3Full5x5_          .push_back(lazyToolnoZS.e1x3(*((*iPho).superCluster()->seed())));
    gen_phoE1x5Full5x5_          .push_back(iPho->full5x5_e1x5());
    gen_phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
    gen_phoE2x5MaxFull5x5_       .push_back(iPho->full5x5_e2x5());
    gen_phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
    gen_phoR9Full5x5_            .push_back(iPho->full5x5_r9());

    if (isAOD_) {
      size_t rightRecoPho = -1;
      for (size_t iv = 0; iv < recoPhotonHandle->size(); ++iv) {
        reco::PhotonRef recophoRef2(recoPhotonHandle, iv);
        if (deltaR(iPho->eta(), iPho->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
      }

      reco::PhotonRef recophoRef(recoPhotonHandle, rightRecoPho);
      reco::Vertex pv = recVtxs->at(0);
      GEDIdTool->setPhotonP4(recophoRef, pv);

      //phoPFChIso_       .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h));
      //phoPFPhoIso_      .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::gamma));
      //phoPFNeuIso_      .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h0));

      std::vector<double> IsoRings;
      GEDIdTool->FrixioneIso(0.1, 8, reco::PFCandidate::h, IsoRings);
      gen_phoPFChIsoFrix1_.push_back(IsoRings[0]);
      gen_phoPFChIsoFrix2_.push_back(IsoRings[1]);
      gen_phoPFChIsoFrix3_.push_back(IsoRings[2]);
      gen_phoPFChIsoFrix4_.push_back(IsoRings[3]);
      gen_phoPFChIsoFrix5_.push_back(IsoRings[4]);
      gen_phoPFChIsoFrix6_.push_back(IsoRings[5]);
      gen_phoPFChIsoFrix7_.push_back(IsoRings[6]);
      gen_phoPFChIsoFrix8_.push_back(IsoRings[7]);
      IsoRings.resize(0);

      GEDIdTool->FrixioneIso(0.1, 8, reco::PFCandidate::gamma, IsoRings);
      gen_phoPFPhoIsoFrix1_.push_back(IsoRings[0]);
      gen_phoPFPhoIsoFrix2_.push_back(IsoRings[1]);
      gen_phoPFPhoIsoFrix3_.push_back(IsoRings[2]);
      gen_phoPFPhoIsoFrix4_.push_back(IsoRings[3]);
      gen_phoPFPhoIsoFrix5_.push_back(IsoRings[4]);
      gen_phoPFPhoIsoFrix6_.push_back(IsoRings[5]);
      gen_phoPFPhoIsoFrix7_.push_back(IsoRings[6]);
      gen_phoPFPhoIsoFrix8_.push_back(IsoRings[7]);
      IsoRings.resize(0);

      GEDIdTool->FrixioneIso(0.1, 8, reco::PFCandidate::h0, IsoRings);
      gen_phoPFNeuIsoFrix1_.push_back(IsoRings[0]);
      gen_phoPFNeuIsoFrix2_.push_back(IsoRings[1]);
      gen_phoPFNeuIsoFrix3_.push_back(IsoRings[2]);
      gen_phoPFNeuIsoFrix4_.push_back(IsoRings[3]);
      gen_phoPFNeuIsoFrix5_.push_back(IsoRings[4]);
      gen_phoPFNeuIsoFrix6_.push_back(IsoRings[5]);
      gen_phoPFNeuIsoFrix7_.push_back(IsoRings[6]);
      gen_phoPFNeuIsoFrix8_.push_back(IsoRings[7]);

      std::vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      //phoPFChWorstIso_  .push_back(*max_element(vtxIsolations03.begin(), vtxIsolations03.end()));
    }
      
    //
    //Filling Frix in case this is not AOD
    /*
    if (!isAOD_) {

      e.getByToken(pckPFCdsLabel_, pfCndHandle); 
      e.getByToken(recoCdsLabel_,  CndHandle); 
      auto_ptr<vector<pat::PackedCandidate> > CandColl( new vector<pat::PackedCandidate> (*pfCndHandle) );

      //register 8 variables for each type of Iso corresponding to each frix. ring 
      double IsoRNeu[8] = {0};
      double IsoRPho[8] = {0};
      double IsoRChg[8] = {0};
      
      //Quantities for charged hadron subtraction
      const float dxyMax     = 0.1; 
      const float dzMax      = 0.2; 
      reco::Vertex pv = recVtxs->at(0);   

      //Loop over PFCandidates in the event
      for(uint ika = 0; ika < CandColl->size() ;ika++){ //PFCand Loop
	pat::PackedCandidate & pfCand = (*CandColl)[ika];
	const auto& iCand = CndHandle->ptrAt(ika);
	
	double DRgamma_cand = deltaR(iPho->eta(),iPho->phi(),iCand->eta(),iCand->phi()); 
	
	if(DRgamma_cand <  0.8){
	  bool isF = isInFootprint(iPho->associatedPackedPFCandidates(),iCand); 
	  if(isF == 0){
	    int ir = DRgamma_cand*10;	    
	    if( pfCand.pdgId() == 22  ) IsoRPho[ir] += pfCand.pt();
	    if( pfCand.pdgId() == 130 ) IsoRNeu[ir] += pfCand.pt();
	    if( pfCand.pdgId() == 211 ){
	      float dz  = fabs(pfCand.pseudoTrack().dz(pv.position()));
	      float dxy = fabs(pfCand.pseudoTrack().dxy(pv.position()));
	      if( dxyMax > dxy ){
		if( dzMax > dz ){
		  IsoRChg[ir] += pfCand.pt();
		}
	      } 
	    }// Charge Hadron Iso
	  } // Is not in Photon Footprint
	} // DR photon cand check
      }//eof loop on pfCands 
      phoPFChIsoFrix1_.push_back(IsoRChg[0]);
      phoPFChIsoFrix2_.push_back(IsoRChg[1]);
      phoPFChIsoFrix3_.push_back(IsoRChg[2]);
      phoPFChIsoFrix4_.push_back(IsoRChg[3]);
      phoPFChIsoFrix5_.push_back(IsoRChg[4]);
      phoPFChIsoFrix6_.push_back(IsoRChg[5]);
      phoPFChIsoFrix7_.push_back(IsoRChg[6]);
      phoPFChIsoFrix8_.push_back(IsoRChg[7]);
   
      phoPFPhoIsoFrix1_.push_back(IsoRPho[0]);
      phoPFPhoIsoFrix2_.push_back(IsoRPho[1]);
      phoPFPhoIsoFrix3_.push_back(IsoRPho[2]);
      phoPFPhoIsoFrix4_.push_back(IsoRPho[3]);
      phoPFPhoIsoFrix5_.push_back(IsoRPho[4]);
      phoPFPhoIsoFrix6_.push_back(IsoRPho[5]);
      phoPFPhoIsoFrix7_.push_back(IsoRPho[6]);
      phoPFPhoIsoFrix8_.push_back(IsoRPho[7]);
   
      phoPFNeuIsoFrix1_.push_back(IsoRNeu[0]);
      phoPFNeuIsoFrix2_.push_back(IsoRNeu[1]);
      phoPFNeuIsoFrix3_.push_back(IsoRNeu[2]);
      phoPFNeuIsoFrix4_.push_back(IsoRNeu[3]);
      phoPFNeuIsoFrix5_.push_back(IsoRNeu[4]);
      phoPFNeuIsoFrix6_.push_back(IsoRNeu[5]);
      phoPFNeuIsoFrix7_.push_back(IsoRNeu[6]);
      phoPFNeuIsoFrix8_.push_back(IsoRNeu[7]);

    } // is not AOD for frix calculations 
    */

    //phoSeedBCE_        .push_back((*iPho).superCluster()->seed()->energy());
    //phoSeedBCEta_      .push_back((*iPho).superCluster()->seed()->eta());
    /*
    phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster())));
    phoMIPChi2_        .push_back(iPho->mipChi2());
    phoMIPTotEnergy_   .push_back(iPho->mipTotEnergy());
    phoMIPSlope_       .push_back(iPho->mipSlope());
    phoMIPIntercept_   .push_back(iPho->mipIntercept());
    phoMIPNhitCone_    .push_back(iPho->mipNhitCone());
    phoMIPIsHalo_      .push_back(iPho->mipIsHalo());
    */
    ///SJ - isolation variables
    //phoEcalRecHitSumEtConeDR03_                   .push_back(iPho->ecalRecHitSumEtConeDR03());
    //phohcalDepth1TowerSumEtConeDR03_              .push_back(iPho->hcalDepth1TowerSumEtConeDR03());
    //phohcalDepth2TowerSumEtConeDR03_              .push_back(iPho->hcalDepth2TowerSumEtConeDR03());
    //phohcalTowerSumEtConeDR03_                    .push_back(iPho->hcalTowerSumEtConeDR03());
    //photrkSumPtHollowConeDR03_                    .push_back(iPho->trkSumPtHollowConeDR03());
    //photrkSumPtSolidConeDR03_                     .push_back(iPho->trkSumPtSolidConeDR03());

    const auto pho = photonHandle->ptrAt(gen_nPho_);
    
    UShort_t tmpphoIDbit = 0;
    
    gen_phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho]);
    gen_phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho]);
    gen_phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho]);
    gen_phoPFChWorstIso_         .push_back((*phoWorstChargedIsolationMap)[pho]);
    
    //phoCITKChIso_            .push_back((*phoChargedIsolationMap_CITK)[pho]);
    //phoCITKPhoIso_           .push_back((*phoPhotonIsolationMap_CITK)[pho]);
    //phoCITKNeuIso_           .push_back((*phoNeutralHadronIsolationMap_CITK)[pho]);
    
    bool isPassLoose  = (*loose_id_decisions)[pho];
    if(isPassLoose) setbit(tmpphoIDbit, 0);
    
    bool isPassMedium = (*medium_id_decisions)[pho];
    if(isPassMedium) setbit(tmpphoIDbit, 1);
    
    bool isPassTight  = (*tight_id_decisions)[pho];
    if(isPassTight) setbit(tmpphoIDbit, 2);
    
    gen_phoIDMVA_.push_back((*mvaValues)[pho]);  
    gen_phoIDbit_.push_back(tmpphoIDbit);      
    
    gen_nPho_++;
  }
  
  pho_gen_ind.clear();

  if (GEDIdTool) delete GEDIdTool;
}

//void ggNtuplizer::cleanupPhotons() {
//}

bool ggNtuplizer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle)
{
//particle is already the ancestor
        if(ancestor == particle ) return true;

//otherwise loop on mothers, if any and return true if the ancestor is found
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
//if we did not return yet, then particle and ancestor are not relatives
        return false;
}
