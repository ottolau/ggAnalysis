#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
Int_t          nEle_;
vector<int>    eleCharge_lead_;
vector<int>    eleChargeConsistent_lead_;
vector<float>  eleEn_lead_;
vector<float>  eleSCEn_lead_;
vector<float>  eleEcalEn_lead_;
vector<float>  eleESEnP1_lead_;
vector<float>  eleESEnP2_lead_;
//vector<float>  eleESEnP1Raw_lead_;
//vector<float>  eleESEnP2Raw_lead_;
vector<float>  eleD0_lead_;
vector<float>  eleDz_lead_;
vector<float>  eleZ_lead_;
vector<float>  eleD0Error_lead_;
vector<float>  eleDzError_lead_;
vector<float>  eleSIP_lead_;
vector<float>  elePt_lead_;
vector<float>  eleEta_lead_;
vector<float>  elePhi_lead_;
vector<float>  eleR9_lead_;
vector<float>  eleCalibPt_lead_;
vector<float>  eleCalibEn_lead_;
vector<float>  eleSCEta_lead_;
vector<float>  eleSCPhi_lead_;
vector<float>  eleSCRawEn_lead_;
vector<float>  eleSCEtaWidth_lead_;
vector<float>  eleSCPhiWidth_lead_;
vector<float>  eleHoverE_lead_;
vector<float>  eleEoverP_lead_;
vector<float>  eleEoverPout_lead_;
vector<float>  eleEoverPInv_lead_;
vector<float>  eleBrem_lead_;
vector<float>  eledEtaAtVtx_lead_;
vector<float>  eledPhiAtVtx_lead_;
vector<float>  eledEtaAtCalo_lead_;
//vector<float>  eleSigmaIEtaIEta_lead_;
//vector<float>  eleSigmaIEtaIPhi_lead_;
//vector<float>  eleSigmaIPhiIPhi_lead_;
vector<float>  eleSigmaIEtaIEtaFull5x5_lead_;
vector<float>  eleSigmaIPhiIPhiFull5x5_lead_;
vector<int>    eleConvVeto_lead_;
vector<int>    eleMissHits_lead_;
vector<float>  eleESEffSigmaRR_lead_;
vector<float>  elePFChIso_lead_;
vector<float>  elePFPhoIso_lead_;
vector<float>  elePFNeuIso_lead_;
vector<float>  elePFPUIso_lead_;
vector<float>  elePFClusEcalIso_lead_;
vector<float>  elePFClusHcalIso_lead_;
//vector<float>  elePFMiniIso_lead_;
vector<float>  eleIDMVAIso_lead_;
vector<float>  eleIDMVANoIso_lead_;
vector<float>  eledEtaseedAtVtx_lead_;
vector<float>  eleE1x5_lead_;
vector<float>  eleE2x5_lead_;
vector<float>  eleE5x5_lead_;
vector<float>  eleE1x5Full5x5_lead_;
vector<float>  eleE2x5Full5x5_lead_;
vector<float>  eleE5x5Full5x5_lead_;
vector<float>  eleR9Full5x5_lead_;
vector<int>    eleEcalDrivenSeed_lead_;
vector<float>  eleDr03EcalRecHitSumEt_lead_;
vector<float>  eleDr03HcalDepth1TowerSumEt_lead_;
vector<float>  eleDr03HcalDepth2TowerSumEt_lead_;
vector<float>  eleDr03HcalTowerSumEt_lead_;
vector<float>  eleDr03TkSumPt_lead_;
vector<float>  elecaloEnergy_lead_;
vector<float>  eleTrkdxy_lead_;
vector<float>  eleKFHits_lead_;
vector<float>  eleKFChi2_lead_;
vector<float>  eleGSFChi2_lead_;
//vector<unsigned long> eleFiredSingleTrgs_lead_;
//vector<unsigned long> eleFiredDoubleTrgs_lead_;
//vector<unsigned long> eleFiredL1Trgs_lead_;
//vector<unsigned short> eleIDbit_lead_;
vector<UShort_t> eleIDbit_lead_;

vector<int>    eleCharge_sublead_;
vector<int>    eleChargeConsistent_sublead_;
vector<float>  eleEn_sublead_;
vector<float>  eleSCEn_sublead_;
vector<float>  eleEcalEn_sublead_;
vector<float>  eleESEnP1_sublead_;
vector<float>  eleESEnP2_sublead_;
//vector<float>  eleESEnP1Raw_sublead_;
//vector<float>  eleESEnP2Raw_sublead_;
vector<float>  eleD0_sublead_;
vector<float>  eleDz_sublead_;
vector<float>  eleZ_sublead_;
vector<float>  eleD0Error_sublead_;
vector<float>  eleDzError_sublead_;
vector<float>  eleSIP_sublead_;
vector<float>  elePt_sublead_;
vector<float>  eleEta_sublead_;
vector<float>  elePhi_sublead_;
vector<float>  eleR9_sublead_;
vector<float>  eleCalibPt_sublead_;
vector<float>  eleCalibEn_sublead_;
vector<float>  eleSCEta_sublead_;
vector<float>  eleSCPhi_sublead_;
vector<float>  eleSCRawEn_sublead_;
vector<float>  eleSCEtaWidth_sublead_;
vector<float>  eleSCPhiWidth_sublead_;
vector<float>  eleHoverE_sublead_;
vector<float>  eleEoverP_sublead_;
vector<float>  eleEoverPout_sublead_;
vector<float>  eleEoverPInv_sublead_;
vector<float>  eleBrem_sublead_;
vector<float>  eledEtaAtVtx_sublead_;
vector<float>  eledPhiAtVtx_sublead_;
vector<float>  eledEtaAtCalo_sublead_;
//vector<float>  eleSigmaIEtaIEta_sublead_;
//vector<float>  eleSigmaIEtaIPhi_sublead_;
//vector<float>  eleSigmaIPhiIPhi_sublead_;
vector<float>  eleSigmaIEtaIEtaFull5x5_sublead_;
vector<float>  eleSigmaIPhiIPhiFull5x5_sublead_;
vector<int>    eleConvVeto_sublead_;
vector<int>    eleMissHits_sublead_;
vector<float>  eleESEffSigmaRR_sublead_;
vector<float>  elePFChIso_sublead_;
vector<float>  elePFPhoIso_sublead_;
vector<float>  elePFNeuIso_sublead_;
vector<float>  elePFPUIso_sublead_;
vector<float>  elePFClusEcalIso_sublead_;
vector<float>  elePFClusHcalIso_sublead_;
//vector<float>  elePFMiniIso_sublead_;
vector<float>  eleIDMVAIso_sublead_;
vector<float>  eleIDMVANoIso_sublead_;
vector<float>  eledEtaseedAtVtx_sublead_;
vector<float>  eleE1x5_sublead_;
vector<float>  eleE2x5_sublead_;
vector<float>  eleE5x5_sublead_;
vector<float>  eleE1x5Full5x5_sublead_;
vector<float>  eleE2x5Full5x5_sublead_;
vector<float>  eleE5x5Full5x5_sublead_;
vector<float>  eleR9Full5x5_sublead_;
vector<int>    eleEcalDrivenSeed_sublead_;
vector<float>  eleDr03EcalRecHitSumEt_sublead_;
vector<float>  eleDr03HcalDepth1TowerSumEt_sublead_;
vector<float>  eleDr03HcalDepth2TowerSumEt_sublead_;
vector<float>  eleDr03HcalTowerSumEt_sublead_;
vector<float>  eleDr03TkSumPt_sublead_;
vector<float>  elecaloEnergy_sublead_;
vector<float>  eleTrkdxy_sublead_;
vector<float>  eleKFHits_sublead_;
vector<float>  eleKFChi2_sublead_;
vector<float>  eleGSFChi2_sublead_;
//vector<unsigned long> eleFiredSingleTrgs_sublead_;
//vector<unsigned long> eleFiredDoubleTrgs_sublead_;
//vector<unsigned long> eleFiredL1Trgs_sublead_;
//vector<unsigned short> eleIDbit_sublead_;
vector<UShort_t> eleIDbit_sublead_;

vector<float> eleZmass_;

vector<float> eleSvChi2_;
vector<float> eleSvNDOF_;
vector<float> eleSvProb_;
vector<float> eleSvX_;
vector<float> eleSvY_;
vector<float> eleSvZ_;
vector<float> eleSvXError_;
vector<float> eleSvYError_;
vector<float> eleSvZError_;
vector<float> eleSvMass_;
vector<float> eleSvCtxy_;
vector<float> eleSvCosAngle_;
vector<float> eleSvLxy_;
vector<float> eleSvLxyError_;


void ggNtuplizer::branchesElectrons(TTree* tree) {

  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge_lead",               &eleCharge_lead_);
  tree->Branch("eleChargeConsistent_lead",     &eleChargeConsistent_lead_);
  tree->Branch("eleEn_lead",                   &eleEn_lead_);
  tree->Branch("eleSCEn_lead",                 &eleSCEn_lead_);
  tree->Branch("eleEcalEn_lead",               &eleEcalEn_lead_);
  tree->Branch("eleESEnP1_lead",               &eleESEnP1_lead_);
  tree->Branch("eleESEnP2_lead",               &eleESEnP2_lead_);
  tree->Branch("eleD0_lead",                   &eleD0_lead_);
  tree->Branch("eleDz_lead",                   &eleDz_lead_);
  tree->Branch("eleZ_lead",                    &eleZ_lead_);
  tree->Branch("eleD0Error_lead",              &eleD0Error_lead_);
  tree->Branch("eleDzError_lead",              &eleDzError_lead_);
  tree->Branch("eleSIP_lead",                  &eleSIP_lead_);
  tree->Branch("elePt_lead",                   &elePt_lead_);
  tree->Branch("eleEta_lead",                  &eleEta_lead_);
  tree->Branch("elePhi_lead",                  &elePhi_lead_);
  tree->Branch("eleR9_lead",                   &eleR9_lead_);
  tree->Branch("eleCalibPt_lead",              &eleCalibPt_lead_);
  tree->Branch("eleCalibEn_lead",              &eleCalibEn_lead_);
  tree->Branch("eleSCEta_lead",                &eleSCEta_lead_);
  tree->Branch("eleSCPhi_lead",                &eleSCPhi_lead_);
  tree->Branch("eleSCRawEn_lead",              &eleSCRawEn_lead_);
  tree->Branch("eleSCEtaWidth_lead",           &eleSCEtaWidth_lead_);
  tree->Branch("eleSCPhiWidth_lead",           &eleSCPhiWidth_lead_);
  tree->Branch("eleHoverE_lead",               &eleHoverE_lead_);
  tree->Branch("eleEoverP_lead",               &eleEoverP_lead_);
  tree->Branch("eleEoverPout_lead",            &eleEoverPout_lead_);
  tree->Branch("eleEoverPInv_lead",            &eleEoverPInv_lead_);
  tree->Branch("eleBrem_lead",                 &eleBrem_lead_);
  tree->Branch("eledEtaAtVtx_lead",            &eledEtaAtVtx_lead_);
  tree->Branch("eledPhiAtVtx_lead",            &eledPhiAtVtx_lead_);
  tree->Branch("eledEtaAtCalo_lead",           &eledEtaAtCalo_lead_);
  //tree->Branch("eleSigmaIEtaIEta_lead",        &eleSigmaIEtaIEta_lead_);
  //tree->Branch("eleSigmaIEtaIPhi_lead",        &eleSigmaIEtaIPhi_lead_);
  //tree->Branch("eleSigmaIPhiIPhi_lead",        &eleSigmaIPhiIPhi_lead_);
  tree->Branch("eleSigmaIEtaIEtaFull5x5_lead", &eleSigmaIEtaIEtaFull5x5_lead_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5_lead", &eleSigmaIPhiIPhiFull5x5_lead_);
  tree->Branch("eleConvVeto_lead",             &eleConvVeto_lead_);
  tree->Branch("eleMissHits_lead",             &eleMissHits_lead_);
  tree->Branch("eleESEffSigmaRR_lead",         &eleESEffSigmaRR_lead_);
  tree->Branch("elePFChIso_lead",              &elePFChIso_lead_);
  tree->Branch("elePFPhoIso_lead",             &elePFPhoIso_lead_);
  tree->Branch("elePFNeuIso_lead",             &elePFNeuIso_lead_);
  tree->Branch("elePFPUIso_lead",              &elePFPUIso_lead_);
  tree->Branch("elePFClusEcalIso_lead",        &elePFClusEcalIso_lead_);
  tree->Branch("elePFClusHcalIso_lead",        &elePFClusHcalIso_lead_);
//  tree->Branch("elePFMiniIso_lead",            &elePFMiniIso_lead_);
  tree->Branch("eleIDMVAIso_lead",             &eleIDMVAIso_lead_);
  tree->Branch("eleIDMVANoIso_lead",           &eleIDMVANoIso_lead_);
  tree->Branch("eledEtaseedAtVtx_lead",        &eledEtaseedAtVtx_lead_);
  tree->Branch("eleE1x5_lead",                 &eleE1x5_lead_);
  tree->Branch("eleE2x5_lead",                 &eleE2x5_lead_);
  tree->Branch("eleE5x5_lead",                 &eleE5x5_lead_);
  tree->Branch("eleE1x5Full5x5_lead",          &eleE1x5Full5x5_lead_);
  tree->Branch("eleE2x5Full5x5_lead",          &eleE2x5Full5x5_lead_);
  tree->Branch("eleE5x5Full5x5_lead",          &eleE5x5Full5x5_lead_);
  tree->Branch("eleR9Full5x5_lead",                &eleR9Full5x5_lead_);
  tree->Branch("eleEcalDrivenSeed_lead",           &eleEcalDrivenSeed_lead_);
  tree->Branch("eleDr03EcalRecHitSumEt_lead",      &eleDr03EcalRecHitSumEt_lead_);
  tree->Branch("eleDr03HcalDepth1TowerSumEt_lead", &eleDr03HcalDepth1TowerSumEt_lead_);
  tree->Branch("eleDr03HcalDepth2TowerSumEt_lead", &eleDr03HcalDepth2TowerSumEt_lead_);
  tree->Branch("eleDr03HcalTowerSumEt_lead",       &eleDr03HcalTowerSumEt_lead_);
  tree->Branch("eleDr03TkSumPt_lead",              &eleDr03TkSumPt_lead_);
  tree->Branch("elecaloEnergy_lead",               &elecaloEnergy_lead_);
  tree->Branch("eleTrkdxy_lead",                   &eleTrkdxy_lead_);
  tree->Branch("eleKFHits_lead",                   &eleKFHits_lead_);
  tree->Branch("eleKFChi2_lead",                   &eleKFChi2_lead_);
  tree->Branch("eleGSFChi2_lead",                  &eleGSFChi2_lead_);
//  tree->Branch("eleFiredSingleTrgs_lead",          &eleFiredSingleTrgs_lead_);
//  tree->Branch("eleFiredDoubleTrgs_lead",          &eleFiredDoubleTrgs_lead_);
//  tree->Branch("eleFiredL1Trgs_lead",              &eleFiredL1Trgs_lead_);
  tree->Branch("eleIDbit_lead",                    &eleIDbit_lead_);

  tree->Branch("eleCharge_sublead",               &eleCharge_sublead_);
  tree->Branch("eleChargeConsistent_sublead",     &eleChargeConsistent_sublead_);
  tree->Branch("eleEn_sublead",                   &eleEn_sublead_);
  tree->Branch("eleSCEn_sublead",                 &eleSCEn_sublead_);
  tree->Branch("eleEcalEn_sublead",               &eleEcalEn_sublead_);
  tree->Branch("eleESEnP1_sublead",               &eleESEnP1_sublead_);
  tree->Branch("eleESEnP2_sublead",               &eleESEnP2_sublead_);
  tree->Branch("eleD0_sublead",                   &eleD0_sublead_);
  tree->Branch("eleDz_sublead",                   &eleDz_sublead_);
  tree->Branch("eleZ_sublead",                    &eleZ_sublead_);
  tree->Branch("eleD0Error_sublead",              &eleD0Error_sublead_);
  tree->Branch("eleDzError_sublead",              &eleDzError_sublead_);
  tree->Branch("eleSIP_sublead",                  &eleSIP_sublead_);
  tree->Branch("elePt_sublead",                   &elePt_sublead_);
  tree->Branch("eleEta_sublead",                  &eleEta_sublead_);
  tree->Branch("elePhi_sublead",                  &elePhi_sublead_);
  tree->Branch("eleR9_sublead",                   &eleR9_sublead_);
  tree->Branch("eleCalibPt_sublead",              &eleCalibPt_sublead_);
  tree->Branch("eleCalibEn_sublead",              &eleCalibEn_sublead_);
  tree->Branch("eleSCEta_sublead",                &eleSCEta_sublead_);
  tree->Branch("eleSCPhi_sublead",                &eleSCPhi_sublead_);
  tree->Branch("eleSCRawEn_sublead",              &eleSCRawEn_sublead_);
  tree->Branch("eleSCEtaWidth_sublead",           &eleSCEtaWidth_sublead_);
  tree->Branch("eleSCPhiWidth_sublead",           &eleSCPhiWidth_sublead_);
  tree->Branch("eleHoverE_sublead",               &eleHoverE_sublead_);
  tree->Branch("eleEoverP_sublead",               &eleEoverP_sublead_);
  tree->Branch("eleEoverPout_sublead",            &eleEoverPout_sublead_);
  tree->Branch("eleEoverPInv_sublead",            &eleEoverPInv_sublead_);
  tree->Branch("eleBrem_sublead",                 &eleBrem_sublead_);
  tree->Branch("eledEtaAtVtx_sublead",            &eledEtaAtVtx_sublead_);
  tree->Branch("eledPhiAtVtx_sublead",            &eledPhiAtVtx_sublead_);
  tree->Branch("eledEtaAtCalo_sublead",           &eledEtaAtCalo_sublead_);
  //tree->Branch("eleSigmaIEtaIEta_sublead",        &eleSigmaIEtaIEta_sublead_);
  //tree->Branch("eleSigmaIEtaIPhi_sublead",        &eleSigmaIEtaIPhi_sublead_);
  //tree->Branch("eleSigmaIPhiIPhi_sublead",        &eleSigmaIPhiIPhi_sublead_);
  tree->Branch("eleSigmaIEtaIEtaFull5x5_sublead", &eleSigmaIEtaIEtaFull5x5_sublead_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5_sublead", &eleSigmaIPhiIPhiFull5x5_sublead_);
  tree->Branch("eleConvVeto_sublead",             &eleConvVeto_sublead_);
  tree->Branch("eleMissHits_sublead",             &eleMissHits_sublead_);
  tree->Branch("eleESEffSigmaRR_sublead",         &eleESEffSigmaRR_sublead_);
  tree->Branch("elePFChIso_sublead",              &elePFChIso_sublead_);
  tree->Branch("elePFPhoIso_sublead",             &elePFPhoIso_sublead_);
  tree->Branch("elePFNeuIso_sublead",             &elePFNeuIso_sublead_);
  tree->Branch("elePFPUIso_sublead",              &elePFPUIso_sublead_);
  tree->Branch("elePFClusEcalIso_sublead",        &elePFClusEcalIso_sublead_);
  tree->Branch("elePFClusHcalIso_sublead",        &elePFClusHcalIso_sublead_);
//  tree->Branch("elePFMiniIso_sublead",            &elePFMiniIso_sublead_);
  tree->Branch("eleIDMVAIso_sublead",             &eleIDMVAIso_sublead_);
  tree->Branch("eleIDMVANoIso_sublead",           &eleIDMVANoIso_sublead_);
  tree->Branch("eledEtaseedAtVtx_sublead",        &eledEtaseedAtVtx_sublead_);
  tree->Branch("eleE1x5_sublead",                 &eleE1x5_sublead_);
  tree->Branch("eleE2x5_sublead",                 &eleE2x5_sublead_);
  tree->Branch("eleE5x5_sublead",                 &eleE5x5_sublead_);
  tree->Branch("eleE1x5Full5x5_sublead",          &eleE1x5Full5x5_sublead_);
  tree->Branch("eleE2x5Full5x5_sublead",          &eleE2x5Full5x5_sublead_);
  tree->Branch("eleE5x5Full5x5_sublead",          &eleE5x5Full5x5_sublead_);
  tree->Branch("eleR9Full5x5_sublead",                &eleR9Full5x5_sublead_);
  tree->Branch("eleEcalDrivenSeed_sublead",           &eleEcalDrivenSeed_sublead_);
  tree->Branch("eleDr03EcalRecHitSumEt_sublead",      &eleDr03EcalRecHitSumEt_sublead_);
  tree->Branch("eleDr03HcalDepth1TowerSumEt_sublead", &eleDr03HcalDepth1TowerSumEt_sublead_);
  tree->Branch("eleDr03HcalDepth2TowerSumEt_sublead", &eleDr03HcalDepth2TowerSumEt_sublead_);
  tree->Branch("eleDr03HcalTowerSumEt_sublead",       &eleDr03HcalTowerSumEt_sublead_);
  tree->Branch("eleDr03TkSumPt_sublead",              &eleDr03TkSumPt_sublead_);
  tree->Branch("elecaloEnergy_sublead",               &elecaloEnergy_sublead_);
  tree->Branch("eleTrkdxy_sublead",                   &eleTrkdxy_sublead_);
  tree->Branch("eleKFHits_sublead",                   &eleKFHits_sublead_);
  tree->Branch("eleKFChi2_sublead",                   &eleKFChi2_sublead_);
  tree->Branch("eleGSFChi2_sublead",                  &eleGSFChi2_sublead_);
//  tree->Branch("eleFiredSingleTrgs_sublead",          &eleFiredSingleTrgs_sublead_);
//  tree->Branch("eleFiredDoubleTrgs_sublead",          &eleFiredDoubleTrgs_sublead_);
//  tree->Branch("eleFiredL1Trgs_sublead",              &eleFiredL1Trgs_sublead_);
  tree->Branch("eleIDbit_sublead",                    &eleIDbit_sublead_);

  tree->Branch("eleZmass", &eleZmass_);
  tree->Branch("eleSvChi2",                    &eleSvChi2_);
  tree->Branch("eleSvNDOF",                    &eleSvNDOF_);
  tree->Branch("eleSvProb",                    &eleSvProb_);
  tree->Branch("eleSvX",                       &eleSvX_);
  tree->Branch("eleSvY",                       &eleSvY_);
  tree->Branch("eleSvZ",                       &eleSvZ_);
  tree->Branch("eleSvXError",                  &eleSvXError_);
  tree->Branch("eleSvYError",                  &eleSvYError_);
  tree->Branch("eleSvZError",                  &eleSvZError_);
  tree->Branch("eleSvMass",                    &eleSvMass_);
  tree->Branch("eleSvCtxy",                    &eleSvCtxy_);
  tree->Branch("eleSvCosAngle",                    &eleSvCosAngle_);
  tree->Branch("eleSvLxy",                    	   &eleSvLxy_);
  tree->Branch("eleSvLxyError",                    &eleSvLxyError_);

}

void ggNtuplizer::fillElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  // cleanup from previous execution
  eleCharge_lead_                  .clear();
  eleChargeConsistent_lead_        .clear();
  eleEn_lead_                      .clear();
  eleSCEn_lead_                    .clear();
  eleEcalEn_lead_                  .clear();
  eleESEnP1_lead_                  .clear();
  eleESEnP2_lead_                  .clear();
  //eleESEnP1Raw_lead_               .clear();
  //eleESEnP2Raw_lead_               .clear();
  eleD0_lead_                      .clear();
  eleDz_lead_                      .clear();
  eleZ_lead_                       .clear();
  eleD0Error_lead_                 .clear();
  eleDzError_lead_                 .clear();
  eleSIP_lead_                     .clear();
  elePt_lead_                      .clear();
  eleEta_lead_                     .clear();
  elePhi_lead_                     .clear();
  eleR9_lead_                      .clear();
  eleCalibPt_lead_                 .clear();
  eleCalibEn_lead_                 .clear();
  eleSCEta_lead_                   .clear();
  eleSCPhi_lead_                   .clear();
  eleSCRawEn_lead_                 .clear();
  eleSCEtaWidth_lead_              .clear();
  eleSCPhiWidth_lead_              .clear();
  eleHoverE_lead_                  .clear();
  eleEoverP_lead_                  .clear();
  eleEoverPout_lead_               .clear();
  eleEoverPInv_lead_               .clear();
  eleBrem_lead_                    .clear();
  eledEtaAtVtx_lead_               .clear();
  eledPhiAtVtx_lead_               .clear();
  eledEtaAtCalo_lead_              .clear();
  //eleSigmaIEtaIEta_lead_           .clear();
  //eleSigmaIEtaIPhi_lead_           .clear();
  //eleSigmaIPhiIPhi_lead_           .clear();
  eleSigmaIEtaIEtaFull5x5_lead_    .clear();
  eleSigmaIPhiIPhiFull5x5_lead_    .clear();
  eleConvVeto_lead_                .clear();
  eleMissHits_lead_                .clear();
  eleESEffSigmaRR_lead_            .clear();
  elePFChIso_lead_                 .clear();
  elePFPhoIso_lead_                .clear();
  elePFNeuIso_lead_                .clear();
  elePFPUIso_lead_                 .clear();
  elePFClusEcalIso_lead_           .clear();
  elePFClusHcalIso_lead_           .clear();
//  elePFMiniIso_lead_               .clear();
  eleIDMVAIso_lead_                .clear();
  eleIDMVANoIso_lead_              .clear();
  eledEtaseedAtVtx_lead_           .clear();
  eleE1x5_lead_                    .clear();
  eleE2x5_lead_                    .clear();
  eleE5x5_lead_                    .clear();
  eleEcalDrivenSeed_lead_          .clear();
  eleDr03EcalRecHitSumEt_lead_     .clear();
  eleDr03HcalDepth1TowerSumEt_lead_.clear();
  eleDr03HcalDepth2TowerSumEt_lead_.clear();
  eleDr03HcalTowerSumEt_lead_      .clear();
  eleDr03TkSumPt_lead_             .clear();
  eleE1x5Full5x5_lead_             .clear();
  eleE2x5Full5x5_lead_             .clear();
  eleE5x5Full5x5_lead_             .clear();
  eleR9Full5x5_lead_               .clear();
  elecaloEnergy_lead_              .clear();
  eleTrkdxy_lead_                  .clear();
  eleKFHits_lead_                  .clear();
  eleKFChi2_lead_                  .clear();
  eleGSFChi2_lead_                 .clear();
//  eleFiredSingleTrgs_lead_         .clear();
//  eleFiredDoubleTrgs_lead_         .clear();
//  eleFiredL1Trgs_lead_             .clear();
//  eleIDbit_lead_                   .clear();
  eleIDbit_lead_                   .clear();

  eleCharge_sublead_                  .clear();
  eleChargeConsistent_sublead_        .clear();
  eleEn_sublead_                      .clear();
  eleSCEn_sublead_                    .clear();
  eleEcalEn_sublead_                  .clear();
  eleESEnP1_sublead_                  .clear();
  eleESEnP2_sublead_                  .clear();
  //eleESEnP1Raw_sublead_               .clear();
  //eleESEnP2Raw_sublead_               .clear();
  eleD0_sublead_                      .clear();
  eleDz_sublead_                      .clear();
  eleZ_sublead_                       .clear();
  eleD0Error_sublead_                 .clear();
  eleDzError_sublead_                 .clear();
  eleSIP_sublead_                     .clear();
  elePt_sublead_                      .clear();
  eleEta_sublead_                     .clear();
  elePhi_sublead_                     .clear();
  eleR9_sublead_                      .clear();
  eleCalibPt_sublead_                 .clear();
  eleCalibEn_sublead_                 .clear();
  eleSCEta_sublead_                   .clear();
  eleSCPhi_sublead_                   .clear();
  eleSCRawEn_sublead_                 .clear();
  eleSCEtaWidth_sublead_              .clear();
  eleSCPhiWidth_sublead_              .clear();
  eleHoverE_sublead_                  .clear();
  eleEoverP_sublead_                  .clear();
  eleEoverPout_sublead_               .clear();
  eleEoverPInv_sublead_               .clear();
  eleBrem_sublead_                    .clear();
  eledEtaAtVtx_sublead_               .clear();
  eledPhiAtVtx_sublead_               .clear();
  eledEtaAtCalo_sublead_              .clear();
  //eleSigmaIEtaIEta_sublead_           .clear();
  //eleSigmaIEtaIPhi_sublead_           .clear();
  //eleSigmaIPhiIPhi_sublead_           .clear();
  eleSigmaIEtaIEtaFull5x5_sublead_    .clear();
  eleSigmaIPhiIPhiFull5x5_sublead_    .clear();
  eleConvVeto_sublead_                .clear();
  eleMissHits_sublead_                .clear();
  eleESEffSigmaRR_sublead_            .clear();
  elePFChIso_sublead_                 .clear();
  elePFPhoIso_sublead_                .clear();
  elePFNeuIso_sublead_                .clear();
  elePFPUIso_sublead_                 .clear();
  elePFClusEcalIso_sublead_           .clear();
  elePFClusHcalIso_sublead_           .clear();
//  elePFMiniIso_sublead_               .clear();
  eleIDMVAIso_sublead_                .clear();
  eleIDMVANoIso_sublead_              .clear();
  eledEtaseedAtVtx_sublead_           .clear();
  eleE1x5_sublead_                    .clear();
  eleE2x5_sublead_                    .clear();
  eleE5x5_sublead_                    .clear();
  eleEcalDrivenSeed_sublead_          .clear();
  eleDr03EcalRecHitSumEt_sublead_     .clear();
  eleDr03HcalDepth1TowerSumEt_sublead_.clear();
  eleDr03HcalDepth2TowerSumEt_sublead_.clear();
  eleDr03HcalTowerSumEt_sublead_      .clear();
  eleDr03TkSumPt_sublead_             .clear();
  eleE1x5Full5x5_sublead_             .clear();
  eleE2x5Full5x5_sublead_             .clear();
  eleE5x5Full5x5_sublead_             .clear();
  eleR9Full5x5_sublead_               .clear();
  elecaloEnergy_sublead_              .clear();
  eleTrkdxy_sublead_                  .clear();
  eleKFHits_sublead_                  .clear();
  eleKFChi2_sublead_                  .clear();
  eleGSFChi2_sublead_                 .clear();
//  eleFiredSingleTrgs_sublead_         .clear();
//  eleFiredDoubleTrgs_sublead_         .clear();
//  eleFiredL1Trgs_sublead_             .clear();
//  eleIDbit_sublead_                   .clear();
  eleIDbit_sublead_                   .clear();
  eleZmass_.clear();
  
  eleSvChi2_.clear();
  eleSvNDOF_.clear();
  eleSvProb_.clear();
  eleSvX_.clear();
  eleSvY_.clear();
  eleSvZ_.clear();
  eleSvXError_.clear();
  eleSvYError_.clear();
  eleSvZError_.clear();
  eleSvMass_.clear();
  eleSvCtxy_.clear();
  eleSvCosAngle_.clear();
  eleSvLxy_.clear();
  eleSvLxyError_.clear();
  

  nEle_ = 0;

  if (isAOD_) {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
    e.getByToken(calibelectronCollection_, calibelectronHandle);

    if (!electronHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
      return;
    }

    if (!calibelectronHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no calibrated pat::Electrons in event";
      return;
    }

    edm::Handle<edm::ValueMap<bool> >  veto_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  heep_id_decisions;
    edm::Handle<edm::ValueMap<float> > eleMVAIsoValues;
    edm::Handle<edm::ValueMap<float> > eleMVANoIsoValues;
    //edm::Handle<edm::ValueMap<float> > elePFClusEcalIsoValues;
    //edm::Handle<edm::ValueMap<float> > elePFClusHcalIsoValues;

    e.getByToken(eleVetoIdMapToken_ ,         veto_id_decisions);
    e.getByToken(eleLooseIdMapToken_ ,        loose_id_decisions);
    e.getByToken(eleMediumIdMapToken_,        medium_id_decisions);
    e.getByToken(eleTightIdMapToken_,         tight_id_decisions);
    e.getByToken(eleHEEPIdMapToken_ ,         heep_id_decisions);
    e.getByToken(eleMVAIsoValuesMapToken_,    eleMVAIsoValues);
    e.getByToken(eleMVANoIsoValuesMapToken_,  eleMVANoIsoValues);
    //e.getByToken(elePFClusEcalIsoToken_,      elePFClusEcalIsoValues);
    //e.getByToken(elePFClusHcalIsoToken_,      elePFClusHcalIsoValues);

    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      //if (iEle->pt() < 0.8) continue;
      //if (fabs(iEle->vz() - pv.z()) < 0.5) continue;

      for (edm::View<pat::Electron>::const_iterator jEle = iEle+1; jEle != electronHandle->end(); ++jEle) {
	//if (jEle->pt() < 0.8) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;
	if (iEle->charge()*jEle->charge() > 0 ) continue;
	//if (fabs(jEle->vz() - pv.z()) < 0.5) continue;
	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv, Z_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	//if ((iele_lv + jele_lv).M() < 2.4 || (iele_lv + jele_lv).M() > 3.8) continue;
	//if ((iele_lv + jele_lv).M() > 5.0) continue;
	Z_lv = iele_lv + jele_lv;
  
  if (Z_lv.M() < 50. || Z_lv.M() >120) continue;
	
  KinematicParticleFactoryFromTransientTrack pFactory;  
	std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;

	XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	KinematicConstrainedVertexFitter kvFitter;
	RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);
  
	if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 || KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
	if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
	KinVtx->movePointerToTheTop();
	RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();


	RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();

	if (DecayVtx->chiSquared() < 0.0) continue;

	auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
	auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;

  eleZmass_.push_back(Z_lv.M());
	double ctxy = ((DecayVtx->position().x() - pv.x())*Z_lv.Px() + (DecayVtx->position().y() - pv.y())*Z_lv.Py())/(pow(Z_lv.Pt(),2))*Z_lv.M();
	
	math::XYZVector perp(Z_lv.Px(), Z_lv.Py(), 0.);
	math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());
  
	eleSvChi2_.push_back(DecayVtx->chiSquared());
	eleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	eleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	eleSvX_.push_back(DecayVtx->position().x());
	eleSvY_.push_back(DecayVtx->position().y());
	eleSvZ_.push_back(DecayVtx->position().z());
	eleSvXError_.push_back(DecayVtx->error().cxx());
	eleSvYError_.push_back(DecayVtx->error().cyy());
	eleSvZError_.push_back(DecayVtx->error().czz());
	eleSvMass_.push_back(Z_lv.M());
	eleSvCtxy_.push_back(ctxy);
	eleSvCosAngle_.push_back(cosAngle);
	eleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	eleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());
  

	Float_t corrPt_lead = -1;
	Float_t corrEn_lead = -1;
	Float_t corrPt_sublead = -1;
	Float_t corrEn_sublead = -1;

	for (edm::View<pat::Electron>::const_iterator iCEle = calibelectronHandle->begin(); iCEle != calibelectronHandle->end(); ++iCEle) {

	  if (fabs(leadEle->eta() - iCEle->eta()) < 0.001 && fabs(leadEle->phi() - iCEle->phi()) < 0.001) {
	    corrPt_lead = iCEle->pt();
	    corrEn_lead = iCEle->energy();
	  }
	  if (fabs(subleadEle->eta() - iCEle->eta()) < 0.001 && fabs(subleadEle->phi() - iCEle->phi()) < 0.001) {
	    corrPt_sublead = iCEle->pt();
	    corrEn_sublead = iCEle->energy();
	  }
	}

	eleCalibPt_lead_        .push_back(corrPt_lead);
	eleCalibEn_lead_        .push_back(corrEn_lead);
	eleCalibPt_sublead_     .push_back(corrPt_sublead);
	eleCalibEn_sublead_     .push_back(corrEn_sublead);

	eleCharge_lead_          .push_back(leadEle->charge());
	eleChargeConsistent_lead_.push_back((Int_t)leadEle->isGsfCtfScPixChargeConsistent());
	eleEn_lead_              .push_back(leadEle->energy());
	eleD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
	eleDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
  eleZ_lead_               .push_back(leadEle->gsfTrack()->dz());
	eleD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
	eleDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
	eleSIP_lead_             .push_back(fabs(leadEle->dB(pat::Electron::PV3D))/leadEle->edB(pat::Electron::PV3D));
	elePt_lead_              .push_back(leadEle->pt());
	eleEta_lead_             .push_back(leadEle->eta());
	elePhi_lead_             .push_back(leadEle->phi());
	eleR9_lead_              .push_back(leadEle->r9());
	eleSCEn_lead_            .push_back(leadEle->superCluster()->energy());
	eleEcalEn_lead_          .push_back(leadEle->ecalEnergy());
	eleESEnP1_lead_          .push_back(leadEle->superCluster()->preshowerEnergyPlane1());
	eleESEnP2_lead_          .push_back(leadEle->superCluster()->preshowerEnergyPlane2());
	eleSCEta_lead_           .push_back(leadEle->superCluster()->eta());
	eleSCPhi_lead_           .push_back(leadEle->superCluster()->phi());
	eleSCRawEn_lead_         .push_back(leadEle->superCluster()->rawEnergy());
	eleSCEtaWidth_lead_      .push_back(leadEle->superCluster()->etaWidth());
	eleSCPhiWidth_lead_      .push_back(leadEle->superCluster()->phiWidth());
	eleHoverE_lead_          .push_back(leadEle->hcalOverEcal());

	//eleFiredSingleTrgs_lead_ .push_back(matchSingleElectronTriggerFilters(leadEle->pt(), leadEle->eta());

	//eleFiredDoubleTrgs_lead_ .push_back(matchDoubleElectronTriggerFilters(leadEle->pt(), leadEle->eta());
	//eleFiredL1Trgs_lead_     .push_back(matchL1TriggerFilters(leadEle->pt(), leadEle->eta(), leadEle->phi()));

	///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
	eleEoverP_lead_          .push_back(leadEle->eSuperClusterOverP());
	eleEoverPout_lead_       .push_back(leadEle->eEleClusterOverPout());
	eleBrem_lead_            .push_back(leadEle->fbrem());
	eledEtaAtVtx_lead_       .push_back(leadEle->deltaEtaSuperClusterTrackAtVtx());
	eledPhiAtVtx_lead_       .push_back(leadEle->deltaPhiSuperClusterTrackAtVtx());
	eledEtaAtCalo_lead_      .push_back(leadEle->deltaEtaSeedClusterTrackAtCalo());
	//eleSigmaIEtaIEta_lead_   .push_back(leadEle->sigmaIetaIeta()); ///new sigmaietaieta
	//eleSigmaIEtaIPhi_lead_   .push_back(leadEle->sigmaIetaIphi());
	//eleSigmaIPhiIPhi_lead_   .push_back(leadEle->sigmaIphiIphi());
	eleConvVeto_lead_        .push_back((Int_t)leadEle->passConversionVeto()); // ConvVtxFit || missHit == 0
	eleMissHits_lead_        .push_back(leadEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
	eleESEffSigmaRR_lead_    .push_back(lazyTool.eseffsirir(*((*leadEle).superCluster())));

	// VID calculation of (1/E - 1/p)
	double leadeleEoverPInv = 1e30;
	if (leadEle->ecalEnergy() != 0 && std::isfinite(leadEle->ecalEnergy()))	leadeleEoverPInv = (1.0 - leadEle->eSuperClusterOverP())/leadEle->ecalEnergy();
	eleEoverPInv_lead_.push_back(leadeleEoverPInv);

	//if (leadEle->ecalEnergy() == 0)   eleEoverPInv_lead_.push_back(1e30);
	//else if (!std::isfinite(leadEle->ecalEnergy()))  eleEoverPInv_lead_.push_back(1e30);
	//else  eleEoverPInv_lead_.push_back((1.0 - leadEle->eSuperClusterOverP())/leadEle->ecalEnergy());

	///HEEP ID
	double leadeledEtaseedAtVtx = leadEle->superCluster().isNonnull() && leadEle->superCluster()->seed().isNonnull() ?
	  leadEle->deltaEtaSuperClusterTrackAtVtx() - leadEle->superCluster()->eta() + leadEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

	eledEtaseedAtVtx_lead_   .push_back(leadeledEtaseedAtVtx);

	eleE1x5_lead_            .push_back(leadEle->e1x5());
	eleE2x5_lead_            .push_back(leadEle->e2x5Max());
	eleE5x5_lead_            .push_back(leadEle->e5x5());

	reco::GsfElectron::PflowIsolationVariables leadpfIso = leadEle->pfIsolationVariables();

	elePFChIso_lead_         .push_back(leadpfIso.sumChargedHadronPt);
	elePFPhoIso_lead_        .push_back(leadpfIso.sumPhotonEt);
	elePFNeuIso_lead_        .push_back(leadpfIso.sumNeutralHadronEt);
	elePFPUIso_lead_         .push_back(leadpfIso.sumPUPt);
	elecaloEnergy_lead_      .push_back(leadEle->caloEnergy());
	//elePFMiniIso_lead_       .push_back(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*leadEle)), 0.05, 0.2, 10., false));

	/////quantities which were used for Run1 - these do not
	///calculated through PF (meaning no energy is subtracted
	///using PF)
	///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
	///line 120

	eleSigmaIEtaIEtaFull5x5_lead_    .push_back(leadEle->full5x5_sigmaIetaIeta());
	eleSigmaIPhiIPhiFull5x5_lead_    .push_back(leadEle->full5x5_sigmaIphiIphi());
	eleE1x5Full5x5_lead_             .push_back(leadEle->full5x5_e1x5());
	eleE2x5Full5x5_lead_             .push_back(leadEle->full5x5_e2x5Max());
	eleE5x5Full5x5_lead_             .push_back(leadEle->full5x5_e5x5());
	eleR9Full5x5_lead_               .push_back(leadEle->full5x5_r9());

	///For HEEP ID
	eleEcalDrivenSeed_lead_          .push_back(leadEle->ecalDrivenSeed());
	eleDr03EcalRecHitSumEt_lead_     .push_back(leadEle->dr03EcalRecHitSumEt());
	eleDr03HcalDepth1TowerSumEt_lead_.push_back(leadEle->dr03HcalDepth1TowerSumEt());
	eleDr03HcalDepth2TowerSumEt_lead_.push_back(leadEle->dr03HcalDepth2TowerSumEt());
	eleDr03HcalTowerSumEt_lead_      .push_back(leadEle->dr03HcalTowerSumEt());
	eleDr03TkSumPt_lead_             .push_back(leadEle->dr03TkSumPt());

	reco::GsfTrackRef leadgsfTrackRef = leadEle->gsfTrack();

	if (leadEle->gsfTrack().isNonnull()) {
	  eleGSFChi2_lead_.push_back(leadgsfTrackRef->normalizedChi2());
	  if (recVtxs->size() > 0)
	    eleTrkdxy_lead_.push_back(leadgsfTrackRef->dxy(recVtxs->front().position()));
	  else
	    eleTrkdxy_lead_.push_back(-999);
	} else {
	  eleGSFChi2_lead_.push_back(999.);
	  eleTrkdxy_lead_.push_back(-999);
	}
	
	reco::TrackRef leadkfTrackRef = leadEle->closestCtfTrackRef();

	if (leadkfTrackRef.isAvailable() && leadkfTrackRef.isNonnull()) {
	  eleKFHits_lead_.push_back(leadkfTrackRef->hitPattern().trackerLayersWithMeasurement());
	  eleKFChi2_lead_.push_back(leadkfTrackRef->normalizedChi2());
	} else {
	  eleKFHits_lead_.push_back(-1.);
	  eleKFChi2_lead_.push_back(999.);
	}

	//edm::Ptr<reco::GsfElectron> recoEl(leadEle);      
	//const auto el = electrons->ptrAt(nEle_);
	const auto leadel = electronHandle->ptrAt(leadEle - electronHandle->begin());
       
	unsigned short leadtmpeleIDbit = 0;
       
	///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works

	bool isPassVeto  = (*veto_id_decisions)[leadel->originalObjectRef()];
	if (isPassVeto) setbit(leadtmpeleIDbit, 0);
    
	bool isPassLoose  = (*loose_id_decisions)[leadel->originalObjectRef()];
	if (isPassLoose) setbit(leadtmpeleIDbit, 1);

	bool isPassMedium = (*medium_id_decisions)[leadel->originalObjectRef()];
	if (isPassMedium) setbit(leadtmpeleIDbit, 2);

	bool isPassTight  = (*tight_id_decisions)[leadel->originalObjectRef()];
	if (isPassTight) setbit(leadtmpeleIDbit, 3);

	bool isPassHEEP = (*heep_id_decisions)[leadel->originalObjectRef()];
	if (isPassHEEP) setbit(leadtmpeleIDbit, 4);

	eleIDMVAIso_lead_  .push_back((*eleMVAIsoValues)[leadel->originalObjectRef()]);
	eleIDMVANoIso_lead_.push_back((*eleMVANoIsoValues)[leadel->originalObjectRef()]);

	elePFClusEcalIso_lead_.push_back(leadEle->ecalPFClusterIso());
	elePFClusHcalIso_lead_.push_back(leadEle->hcalPFClusterIso());

	eleIDbit_lead_.push_back(leadtmpeleIDbit);


	eleCharge_sublead_          .push_back(subleadEle->charge());
	eleChargeConsistent_sublead_.push_back((Int_t)subleadEle->isGsfCtfScPixChargeConsistent());
	eleEn_sublead_              .push_back(subleadEle->energy());
	eleD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
	eleDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
  eleZ_sublead_               .push_back(subleadEle->gsfTrack()->dz());
	eleD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
	eleDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
	eleSIP_sublead_             .push_back(fabs(subleadEle->dB(pat::Electron::PV3D))/subleadEle->edB(pat::Electron::PV3D));
	elePt_sublead_              .push_back(subleadEle->pt());
	eleEta_sublead_             .push_back(subleadEle->eta());
	elePhi_sublead_             .push_back(subleadEle->phi());
	eleR9_sublead_              .push_back(subleadEle->r9());
	eleSCEn_sublead_            .push_back(subleadEle->superCluster()->energy());
	eleEcalEn_sublead_          .push_back(subleadEle->ecalEnergy());
	eleESEnP1_sublead_          .push_back(subleadEle->superCluster()->preshowerEnergyPlane1());
	eleESEnP2_sublead_          .push_back(subleadEle->superCluster()->preshowerEnergyPlane2());
	eleSCEta_sublead_           .push_back(subleadEle->superCluster()->eta());
	eleSCPhi_sublead_           .push_back(subleadEle->superCluster()->phi());
	eleSCRawEn_sublead_         .push_back(subleadEle->superCluster()->rawEnergy());
	eleSCEtaWidth_sublead_      .push_back(subleadEle->superCluster()->etaWidth());
	eleSCPhiWidth_sublead_      .push_back(subleadEle->superCluster()->phiWidth());
	eleHoverE_sublead_          .push_back(subleadEle->hcalOverEcal());

	//eleFiredSingleTrgs_sublead_ .push_back(matchSingleElectronTriggerFilters(subleadEle->pt(), subleadEle->eta());

	//eleFiredDoubleTrgs_sublead_ .push_back(matchDoubleElectronTriggerFilters(subleadEle->pt(), subleadEle->eta());
	//eleFiredL1Trgs_sublead_     .push_back(matchL1TriggerFilters(subleadEle->pt(), subleadEle->eta(), subleadEle->phi()));

	///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
	eleEoverP_sublead_          .push_back(subleadEle->eSuperClusterOverP());
	eleEoverPout_sublead_       .push_back(subleadEle->eEleClusterOverPout());
	eleBrem_sublead_            .push_back(subleadEle->fbrem());
	eledEtaAtVtx_sublead_       .push_back(subleadEle->deltaEtaSuperClusterTrackAtVtx());
	eledPhiAtVtx_sublead_       .push_back(subleadEle->deltaPhiSuperClusterTrackAtVtx());
	eledEtaAtCalo_sublead_      .push_back(subleadEle->deltaEtaSeedClusterTrackAtCalo());
	//eleSigmaIEtaIEta_sublead_   .push_back(subleadEle->sigmaIetaIeta()); ///new sigmaietaieta
	//eleSigmaIEtaIPhi_sublead_   .push_back(subleadEle->sigmaIetaIphi());
	//eleSigmaIPhiIPhi_sublead_   .push_back(subleadEle->sigmaIphiIphi());
	eleConvVeto_sublead_        .push_back((Int_t)subleadEle->passConversionVeto()); // ConvVtxFit || missHit == 0
	eleMissHits_sublead_        .push_back(subleadEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
	eleESEffSigmaRR_sublead_    .push_back(lazyTool.eseffsirir(*((*subleadEle).superCluster())));

	// VID calculation of (1/E - 1/p)
	double subleadeleEoverPInv = 1e30;
	if (subleadEle->ecalEnergy() != 0 && std::isfinite(subleadEle->ecalEnergy()))	subleadeleEoverPInv = (1.0 - subleadEle->eSuperClusterOverP())/subleadEle->ecalEnergy();
	eleEoverPInv_sublead_.push_back(subleadeleEoverPInv);

	//if (subleadEle->ecalEnergy() == 0)   eleEoverPInv_sublead_.push_back(1e30);
	//else if (!std::isfinite(subleadEle->ecalEnergy()))  eleEoverPInv_sublead_.push_back(1e30);
	//else  eleEoverPInv_sublead_.push_back((1.0 - subleadEle->eSuperClusterOverP())/subleadEle->ecalEnergy());

	///HEEP ID
	double subleadeledEtaseedAtVtx = subleadEle->superCluster().isNonnull() && subleadEle->superCluster()->seed().isNonnull() ?
	  subleadEle->deltaEtaSuperClusterTrackAtVtx() - subleadEle->superCluster()->eta() + subleadEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

	eledEtaseedAtVtx_sublead_   .push_back(subleadeledEtaseedAtVtx);

	eleE1x5_sublead_            .push_back(subleadEle->e1x5());
	eleE2x5_sublead_            .push_back(subleadEle->e2x5Max());
	eleE5x5_sublead_            .push_back(subleadEle->e5x5());

	reco::GsfElectron::PflowIsolationVariables subleadpfIso = subleadEle->pfIsolationVariables();

	elePFChIso_sublead_         .push_back(subleadpfIso.sumChargedHadronPt);
	elePFPhoIso_sublead_        .push_back(subleadpfIso.sumPhotonEt);
	elePFNeuIso_sublead_        .push_back(subleadpfIso.sumNeutralHadronEt);
	elePFPUIso_sublead_         .push_back(subleadpfIso.sumPUPt);
	elecaloEnergy_sublead_      .push_back(subleadEle->caloEnergy());
	//elePFMiniIso_sublead_       .push_back(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*subleadEle)), 0.05, 0.2, 10., false));

	/////quantities which were used for Run1 - these do not
	///calculated through PF (meaning no energy is subtracted
	///using PF)
	///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
	///line 120

	eleSigmaIEtaIEtaFull5x5_sublead_    .push_back(subleadEle->full5x5_sigmaIetaIeta());
	eleSigmaIPhiIPhiFull5x5_sublead_    .push_back(subleadEle->full5x5_sigmaIphiIphi());
	eleE1x5Full5x5_sublead_             .push_back(subleadEle->full5x5_e1x5());
	eleE2x5Full5x5_sublead_             .push_back(subleadEle->full5x5_e2x5Max());
	eleE5x5Full5x5_sublead_             .push_back(subleadEle->full5x5_e5x5());
	eleR9Full5x5_sublead_               .push_back(subleadEle->full5x5_r9());

	///For HEEP ID
	eleEcalDrivenSeed_sublead_          .push_back(subleadEle->ecalDrivenSeed());
	eleDr03EcalRecHitSumEt_sublead_     .push_back(subleadEle->dr03EcalRecHitSumEt());
	eleDr03HcalDepth1TowerSumEt_sublead_.push_back(subleadEle->dr03HcalDepth1TowerSumEt());
	eleDr03HcalDepth2TowerSumEt_sublead_.push_back(subleadEle->dr03HcalDepth2TowerSumEt());
	eleDr03HcalTowerSumEt_sublead_      .push_back(subleadEle->dr03HcalTowerSumEt());
	eleDr03TkSumPt_sublead_             .push_back(subleadEle->dr03TkSumPt());

	reco::GsfTrackRef subleadgsfTrackRef = subleadEle->gsfTrack();

	if (subleadEle->gsfTrack().isNonnull()) {
	  eleGSFChi2_sublead_.push_back(subleadgsfTrackRef->normalizedChi2());
	  if (recVtxs->size() > 0)
	    eleTrkdxy_sublead_.push_back(subleadgsfTrackRef->dxy(recVtxs->front().position()));
	  else
	    eleTrkdxy_sublead_.push_back(-999);
	} else {
	  eleGSFChi2_sublead_.push_back(999.);
	  eleTrkdxy_sublead_.push_back(-999);
	}
	
	reco::TrackRef subleadkfTrackRef = subleadEle->closestCtfTrackRef();

	if (subleadkfTrackRef.isAvailable() && subleadkfTrackRef.isNonnull()) {
	  eleKFHits_sublead_.push_back(subleadkfTrackRef->hitPattern().trackerLayersWithMeasurement());
	  eleKFChi2_sublead_.push_back(subleadkfTrackRef->normalizedChi2());
	} else {
	  eleKFHits_sublead_.push_back(-1.);
	  eleKFChi2_sublead_.push_back(999.);
	}

	//edm::Ptr<reco::GsfElectron> recoEl(subleadEle);      
	//const auto el = electrons->ptrAt(nEle_);
	const auto subleadel = electronHandle->ptrAt(subleadEle - electronHandle->begin());
       
	unsigned short subleadtmpeleIDbit = 0;
       
	///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works

	isPassVeto  = (*veto_id_decisions)[subleadel->originalObjectRef()];
	if (isPassVeto) setbit(subleadtmpeleIDbit, 0);
    
	isPassLoose  = (*loose_id_decisions)[subleadel->originalObjectRef()];
	if (isPassLoose) setbit(subleadtmpeleIDbit, 1);

	isPassMedium = (*medium_id_decisions)[subleadel->originalObjectRef()];
	if (isPassMedium) setbit(subleadtmpeleIDbit, 2);

	isPassTight  = (*tight_id_decisions)[subleadel->originalObjectRef()];
	if (isPassTight) setbit(subleadtmpeleIDbit, 3);

	isPassHEEP = (*heep_id_decisions)[subleadel->originalObjectRef()];
	if (isPassHEEP) setbit(subleadtmpeleIDbit, 4);

	eleIDMVAIso_sublead_  .push_back((*eleMVAIsoValues)[subleadel->originalObjectRef()]);
	eleIDMVANoIso_sublead_.push_back((*eleMVANoIsoValues)[subleadel->originalObjectRef()]);

	elePFClusEcalIso_sublead_.push_back(subleadEle->ecalPFClusterIso());
	elePFClusHcalIso_sublead_.push_back(subleadEle->hcalPFClusterIso());

	eleIDbit_sublead_.push_back(subleadtmpeleIDbit);


	nEle_++;
      }
    }

  }else{
    
    edm::Handle<edm::ValueMap<bool> >  veto_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
    edm::Handle<edm::ValueMap<bool> >  heep_id_decisions;
    edm::Handle<edm::ValueMap<float> > eleMVAIsoValues;
    edm::Handle<edm::ValueMap<float> > eleMVANoIsoValues;
    //edm::Handle<edm::ValueMap<float> > elePFClusEcalIsoValues;
    //edm::Handle<edm::ValueMap<float> > elePFClusHcalIsoValues;

    e.getByToken(eleVetoIdMapToken_ ,         veto_id_decisions);
    e.getByToken(eleLooseIdMapToken_ ,        loose_id_decisions);
    e.getByToken(eleMediumIdMapToken_,        medium_id_decisions);
    e.getByToken(eleTightIdMapToken_,         tight_id_decisions);
    e.getByToken(eleHEEPIdMapToken_ ,         heep_id_decisions);
    e.getByToken(eleMVAIsoValuesMapToken_,    eleMVAIsoValues);
    e.getByToken(eleMVANoIsoValuesMapToken_,  eleMVANoIsoValues);
   

 
    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
    e.getByToken(calibelectronCollection_, calibelectronHandle);

    edm::Handle<pat::PackedCandidateCollection> pfcands;
    e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<pat::PackedCandidateCollection> losttracks;
    e.getByToken(lostTracksLabel_, losttracks);

    std::vector<pat::PackedCandidate> alltracks;
    alltracks.reserve(pfcands->size() + losttracks->size());
    alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
    alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

  //  edm::Handle<reco::DeDxDataValueMap> dEdxObjectHandle;
  //  e.getByToken(deDxProducer_, dEdxObjectHandle );
  //  const edm::ValueMap<reco::DeDxData> dEdxColl = *dEdxObjectHandle.product();


    edm::Handle< std::vector< std::pair<edm::Ptr<pat::Electron>, reco::Track> > > eleTrackMap;
    e.getByToken( tok_eleTtk_, eleTrackMap);
    std::vector<std::pair<edm::Ptr<pat::Electron>, reco::Track>> eletrks = *(eleTrackMap.product());

    if (!electronHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
      return;
    }

    if (!calibelectronHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no calibrated pat::Electrons in event";
      return;
    }


    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      //if (fabs(iEle->vz() - pv.z()) > 0.5) continue;

      for (edm::View<pat::Electron>::const_iterator jEle = iEle+1; jEle != electronHandle->end(); ++jEle) {
  //if (iEle->charge()*jEle->charge() > 0.0) continue;
  //if (fabs(jEle->vz() - pv.z()) > 0.5) continue;
          if (iEle->charge()*jEle->charge() > 0) continue;
          float pmass  = 0.0005109989461;
          TLorentzVector iele_lv, jele_lv, Z_lv;
          iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
          jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
          //if ((iele_lv + jele_lv).M() < 2.4 || (iele_lv + jele_lv).M() > 3.8) continue;
          Z_lv = iele_lv + jele_lv;
          
          if (Z_lv.M() > 120.0 || Z_lv.M() < 50.) continue;
	
          KinematicParticleFactoryFromTransientTrack pFactory;  
          std::vector<RefCountedKinematicParticle> XParticles;
          float pmasse = 1.e-6 * pmass;

          XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
          XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

          KinematicConstrainedVertexFitter kvFitter;
          RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);
          
          if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 || KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
          if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
          KinVtx->movePointerToTheTop();
          RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();


          RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();

          if (DecayVtx->chiSquared() < 0.0) continue;

          auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
          auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;

          eleZmass_.push_back(Z_lv.M());
          double ctxy = ((DecayVtx->position().x() - pv.x())*Z_lv.Px() + (DecayVtx->position().y() - pv.y())*Z_lv.Py())/(pow(Z_lv.Pt(),2))*Z_lv.M();
          
          math::XYZVector perp(Z_lv.Px(), Z_lv.Py(), 0.);
          math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
          math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
          double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());
          
          eleSvChi2_.push_back(DecayVtx->chiSquared());
          eleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
          eleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
          eleSvX_.push_back(DecayVtx->position().x());
          eleSvY_.push_back(DecayVtx->position().y());
          eleSvZ_.push_back(DecayVtx->position().z());
          eleSvXError_.push_back(DecayVtx->error().cxx());
          eleSvYError_.push_back(DecayVtx->error().cyy());
          eleSvZError_.push_back(DecayVtx->error().czz());
          eleSvMass_.push_back(Z_lv.M());
          eleSvCtxy_.push_back(ctxy);
          eleSvCosAngle_.push_back(cosAngle);
          eleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
    	    eleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());
      // Accept these 4 tracks as a Bs candidate, fill ntuple
      eleZmass_.push_back(Z_lv.M());
    
      Float_t corrPt_lead = -1;
      Float_t corrEn_lead = -1;
      Float_t corrPt_sublead = -1;
      Float_t corrEn_sublead = -1;

      for (edm::View<pat::Electron>::const_iterator iCEle = calibelectronHandle->begin(); iCEle != calibelectronHandle->end(); ++iCEle) {

        if (fabs(leadEle->eta() - iCEle->eta()) < 0.001 && fabs(leadEle->phi() - iCEle->phi()) < 0.001) {
    corrPt_lead = iCEle->pt();
    corrEn_lead = iCEle->energy();
        }
        if (fabs(subleadEle->eta() - iCEle->eta()) < 0.001 && fabs(subleadEle->phi() - iCEle->phi()) < 0.001) {
    corrPt_sublead = iCEle->pt();
    corrEn_sublead = iCEle->energy();
        }
      }

      eleCalibPt_lead_        .push_back(corrPt_lead);
      eleCalibEn_lead_        .push_back(corrEn_lead);
      eleCalibPt_sublead_     .push_back(corrPt_sublead);
      eleCalibEn_sublead_     .push_back(corrEn_sublead);

      eleCharge_lead_          .push_back(leadEle->charge());
      eleChargeConsistent_lead_.push_back((Int_t)leadEle->isGsfCtfScPixChargeConsistent());
      eleEn_lead_              .push_back(leadEle->energy());
      eleD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
      eleDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
      eleD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
      eleDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
      eleSIP_lead_             .push_back(fabs(leadEle->dB(pat::Electron::PV3D))/leadEle->edB(pat::Electron::PV3D));
      elePt_lead_              .push_back(leadEle->pt());
      eleEta_lead_             .push_back(leadEle->eta());
      elePhi_lead_             .push_back(leadEle->phi());
      eleSCEn_lead_            .push_back(leadEle->superCluster()->energy());
      eleEcalEn_lead_          .push_back(leadEle->ecalEnergy());
      eleHoverE_lead_          .push_back(leadEle->hcalOverEcal());
      eledEtaAtVtx_lead_       .push_back(leadEle->deltaEtaSuperClusterTrackAtVtx());
      eledPhiAtVtx_lead_       .push_back(leadEle->deltaPhiSuperClusterTrackAtVtx());
      eleConvVeto_lead_        .push_back((Int_t)leadEle->passConversionVeto()); // ConvVtxFit || missHit == 0

      reco::GsfElectron::PflowIsolationVariables leadpfIso = leadEle->pfIsolationVariables();

      elePFChIso_lead_         .push_back(leadpfIso.sumChargedHadronPt);
      elePFPhoIso_lead_        .push_back(leadpfIso.sumPhotonEt);
      elePFNeuIso_lead_        .push_back(leadpfIso.sumNeutralHadronEt);
      elePFPUIso_lead_         .push_back(leadpfIso.sumPUPt);
      elecaloEnergy_lead_      .push_back(leadEle->caloEnergy());
    
      eleEcalDrivenSeed_lead_          .push_back(leadEle->ecalDrivenSeed());

      reco::GsfTrackRef leadgsfTrackRef = leadEle->gsfTrack();

      if (leadEle->gsfTrack().isNonnull()) {
        eleGSFChi2_lead_.push_back(leadgsfTrackRef->normalizedChi2());
        if (recVtxs->size() > 0)
          eleTrkdxy_lead_.push_back(leadgsfTrackRef->dxy(recVtxs->front().position()));
        else
          eleTrkdxy_lead_.push_back(-999);
      } else {
        eleGSFChi2_lead_.push_back(999.);
        eleTrkdxy_lead_.push_back(-999);
      }
      const auto leadel = electronHandle->ptrAt(leadEle - electronHandle->begin());
       
      unsigned short leadtmpeleIDbit = 0;
       
        ///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works

      bool isPassVeto  = (*veto_id_decisions)[leadel];
      if (isPassVeto) setbit(leadtmpeleIDbit, 0);

      bool isPassLoose  = (*loose_id_decisions)[leadel];
      if (isPassLoose) setbit(leadtmpeleIDbit, 1);

      bool isPassMedium = (*medium_id_decisions)[leadel];
      if (isPassMedium) setbit(leadtmpeleIDbit, 2);

      bool isPassTight  = (*tight_id_decisions)[leadel];
      if (isPassTight) setbit(leadtmpeleIDbit, 3);

      bool isPassHEEP = (*heep_id_decisions)[leadel];
      if (isPassHEEP) setbit(leadtmpeleIDbit, 4);

      eleIDMVAIso_lead_  .push_back((*eleMVAIsoValues)[leadel]);
      eleIDMVANoIso_lead_.push_back((*eleMVANoIsoValues)[leadel]);

      elePFClusEcalIso_lead_.push_back(leadEle->ecalPFClusterIso());
      elePFClusHcalIso_lead_.push_back(leadEle->hcalPFClusterIso());

      eleIDbit_lead_.push_back(leadtmpeleIDbit);



      eleCharge_sublead_          .push_back(subleadEle->charge());
      eleChargeConsistent_sublead_.push_back((Int_t)subleadEle->isGsfCtfScPixChargeConsistent());
      eleEn_sublead_              .push_back(subleadEle->energy());
      eleD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
      eleDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
      eleD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
      eleDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
      eleSIP_sublead_             .push_back(fabs(subleadEle->dB(pat::Electron::PV3D))/subleadEle->edB(pat::Electron::PV3D));
      elePt_sublead_              .push_back(subleadEle->pt());
      eleEta_sublead_             .push_back(subleadEle->eta());
      elePhi_sublead_             .push_back(subleadEle->phi());
      eleSCEn_sublead_            .push_back(subleadEle->superCluster()->energy());
      eleEcalEn_sublead_          .push_back(subleadEle->ecalEnergy());
      eleHoverE_sublead_          .push_back(subleadEle->hcalOverEcal());
      eledEtaAtVtx_sublead_       .push_back(subleadEle->deltaEtaSuperClusterTrackAtVtx());
      eledPhiAtVtx_sublead_       .push_back(subleadEle->deltaPhiSuperClusterTrackAtVtx());
      eleConvVeto_sublead_        .push_back((Int_t)subleadEle->passConversionVeto()); // ConvVtxFit || missHit == 0
      
      reco::GsfElectron::PflowIsolationVariables subleadpfIso = subleadEle->pfIsolationVariables();

      elePFChIso_sublead_         .push_back(subleadpfIso.sumChargedHadronPt);
      elePFPhoIso_sublead_        .push_back(subleadpfIso.sumPhotonEt);
      elePFNeuIso_sublead_        .push_back(subleadpfIso.sumNeutralHadronEt);
      elePFPUIso_sublead_         .push_back(subleadpfIso.sumPUPt);
      elecaloEnergy_sublead_      .push_back(subleadEle->caloEnergy());

      eleEcalDrivenSeed_sublead_          .push_back(subleadEle->ecalDrivenSeed());
 
      reco::GsfTrackRef subleadgsfTrackRef = subleadEle->gsfTrack();

      if (subleadEle->gsfTrack().isNonnull()) {
        eleGSFChi2_sublead_.push_back(subleadgsfTrackRef->normalizedChi2());
        if (recVtxs->size() > 0)
          eleTrkdxy_sublead_.push_back(subleadgsfTrackRef->dxy(recVtxs->front().position()));
        else
          eleTrkdxy_sublead_.push_back(-999);
      } else {
       eleGSFChi2_sublead_.push_back(999.);
        eleTrkdxy_sublead_.push_back(-999);
      }
      
      const auto subleadel = electronHandle->ptrAt(subleadEle - electronHandle->begin());
       
      unsigned short subleadtmpeleIDbit = 0;

        ///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works

      isPassVeto  = (*veto_id_decisions)[subleadel];
      if (isPassVeto) setbit(subleadtmpeleIDbit, 0);

      isPassLoose  = (*loose_id_decisions)[subleadel];
      if (isPassLoose) setbit(subleadtmpeleIDbit, 1);

      isPassMedium = (*medium_id_decisions)[subleadel];
      if (isPassMedium) setbit(subleadtmpeleIDbit, 2);

      isPassTight  = (*tight_id_decisions)[subleadel];
      if (isPassTight) setbit(subleadtmpeleIDbit, 3);

      isPassHEEP = (*heep_id_decisions)[subleadel];
      if (isPassHEEP) setbit(subleadtmpeleIDbit, 4);

      eleIDMVAIso_sublead_  .push_back((*eleMVAIsoValues)[subleadel]);
      eleIDMVANoIso_sublead_.push_back((*eleMVANoIsoValues)[subleadel]);

      elePFClusEcalIso_sublead_.push_back(subleadEle->ecalPFClusterIso());
      elePFClusHcalIso_sublead_.push_back(subleadEle->hcalPFClusterIso());

      eleIDbit_sublead_.push_back(subleadtmpeleIDbit);

      nEle_++;
      }
    }
  }

}


