#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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


using namespace std;

// (local) variables associated with tree branches
Int_t            nMu_;
bool             matchedTrg_;
vector<float>    muPt_lead_;
vector<float>    muEn_lead_;
vector<float>    muEta_lead_;
vector<float>    muPhi_lead_;
vector<int>      muCharge_lead_;
vector<int>      muType_lead_;
vector<UShort_t> muIDbit_lead_;
vector<float>    muD0_lead_;
vector<float>    muDz_lead_;
vector<float>    muSIP_lead_;
vector<float>    muD0Error_lead_;
vector<float>    muDzError_lead_;
vector<float>    muChi2NDF_lead_;
vector<float>    muInnerD0_lead_;
vector<float>    muInnerDz_lead_;
vector<int>      muTrkLayers_lead_;
vector<int>      muPixelLayers_lead_;
vector<int>      muPixelHits_lead_;
vector<int>      muMuonHits_lead_;
vector<int>      muStations_lead_;
vector<int>      muMatches_lead_;
vector<int>      muTrkQuality_lead_;
vector<float>    muIsoTrk_lead_;
vector<float>    muPFChIso_lead_;
vector<float>    muPFPhoIso_lead_;
vector<float>    muPFNeuIso_lead_;
vector<float>    muPFPUIso_lead_;
vector<float>    muPFChIso03_lead_;
vector<float>    muPFPhoIso03_lead_;
vector<float>    muPFNeuIso03_lead_;
vector<float>    muPFPUIso03_lead_;
//vector<float>    muPFMiniIso_lead_;
vector<bool>     muFiredTrgs_lead_;
vector<ULong64_t> muFiredL1Trgs_lead_;
vector<float>    muInnervalidFraction_lead_;
vector<float>    musegmentCompatibility_lead_;
vector<float>    muchi2LocalPosition_lead_;
vector<float>    mutrkKink_lead_;
vector<float>    muBestTrkPtError_lead_;
vector<float>    muBestTrkPt_lead_;
vector<int>      muBestTrkType_lead_;
vector<float>    muTrkNormChi2_lead_;
vector<float>    muIDPatMVA_lead_;
vector<float>    muIDPatSoftMVA_lead_;
vector<UShort_t> muIDSelectorBit_lead_;

vector<float>    muPt_sublead_;
vector<float>    muEn_sublead_;
vector<float>    muEta_sublead_;
vector<float>    muPhi_sublead_;
vector<int>      muCharge_sublead_;
vector<int>      muType_sublead_;
vector<UShort_t> muIDbit_sublead_;
vector<float>    muD0_sublead_;
vector<float>    muDz_sublead_;
vector<float>    muSIP_sublead_;
vector<float>    muD0Error_sublead_;
vector<float>    muDzError_sublead_;
vector<float>    muChi2NDF_sublead_;
vector<float>    muInnerD0_sublead_;
vector<float>    muInnerDz_sublead_;
vector<int>      muTrkLayers_sublead_;
vector<int>      muPixelLayers_sublead_;
vector<int>      muPixelHits_sublead_;
vector<int>      muMuonHits_sublead_;
vector<int>      muStations_sublead_;
vector<int>      muMatches_sublead_;
vector<int>      muTrkQuality_sublead_;
vector<float>    muIsoTrk_sublead_;
vector<float>    muPFChIso_sublead_;
vector<float>    muPFPhoIso_sublead_;
vector<float>    muPFNeuIso_sublead_;
vector<float>    muPFPUIso_sublead_;
vector<float>    muPFChIso03_sublead_;
vector<float>    muPFPhoIso03_sublead_;
vector<float>    muPFNeuIso03_sublead_;
vector<float>    muPFPUIso03_sublead_;
//vector<float>    muPFMiniIso_sublead_;
vector<bool>     muFiredTrgs_sublead_;
vector<ULong64_t> muFiredL1Trgs_sublead_;
vector<float>    muInnervalidFraction_sublead_;
vector<float>    musegmentCompatibility_sublead_;
vector<float>    muchi2LocalPosition_sublead_;
vector<float>    mutrkKink_sublead_;
vector<float>    muBestTrkPtError_sublead_;
vector<float>    muBestTrkPt_sublead_;
vector<int>      muBestTrkType_sublead_;
vector<float>    muTrkNormChi2_sublead_;
vector<float>    muIDPatMVA_sublead_;
vector<float>    muIDPatSoftMVA_sublead_;
vector<UShort_t> muIDSelectorBit_sublead_;

vector<float> muSvChi2_;
vector<float> muSvNDOF_;
vector<float> muSvProb_;
vector<float> muSvX_;
vector<float> muSvY_;
vector<float> muSvZ_;
vector<float> muSvXError_;
vector<float> muSvYError_;
vector<float> muSvZError_;
vector<float> muSvMass_;
vector<float> muSvCtxy_;
vector<float> muSvCosAngle_;
vector<float> muSvLxy_;
vector<float> muSvLxyError_;
/*
vector<int>    kaonMMCharge_lead_;
vector<float>  kaonMMD0_lead_;
vector<float>  kaonMMDz_lead_;
vector<float>  kaonMMD0Error_lead_;
vector<float>  kaonMMDzError_lead_;
vector<float>  kaonMMPt_lead_;
vector<float>  kaonMMEta_lead_;
vector<float>  kaonMMPhi_lead_;
vector<float>  kaonMMVx_lead_;
vector<float>  kaonMMVy_lead_;
vector<float>  kaonMMVz_lead_;
vector<float>  kaonMMEn_lead_;
vector<float>  kaonMMTrkChi2_lead_;
vector<float>  kaonMMTrkNDOF_lead_;
vector<float>  kaonMMTrkNormChi2_lead_;

vector<int>    kaonMMCharge_sublead_;
vector<float>  kaonMMD0_sublead_;
vector<float>  kaonMMDz_sublead_;
vector<float>  kaonMMD0Error_sublead_;
vector<float>  kaonMMDzError_sublead_;
vector<float>  kaonMMPt_sublead_;
vector<float>  kaonMMEta_sublead_;
vector<float>  kaonMMPhi_sublead_;
vector<float>  kaonMMVx_sublead_;
vector<float>  kaonMMVy_sublead_;
vector<float>  kaonMMVz_sublead_;
vector<float>  kaonMMEn_sublead_;
vector<float>  kaonMMTrkChi2_sublead_;
vector<float>  kaonMMTrkNDOF_sublead_;
vector<float>  kaonMMTrkNormChi2_sublead_;

vector<float>  bsMMdRmu_;
vector<float>  bsMMdRkaon_;
vector<float>  bsMMdRJpsiPhi_;
vector<float>  bsMMJpsiMass_;
vector<float>  bsMMPhiMass_;
vector<float>  bsMMBsMass_;
*/

void ggNtuplizer::branchesMuons(TTree* tree) {

  tree->Branch("nMu",           &nMu_);
  tree->Branch("matchedTrg",    &matchedTrg_);
  tree->Branch("muPt_lead",          &muPt_lead_);
  tree->Branch("muEn_lead",          &muEn_lead_);
  tree->Branch("muEta_lead",         &muEta_lead_);
  tree->Branch("muPhi_lead",         &muPhi_lead_);
  tree->Branch("muCharge_lead",      &muCharge_lead_);
  tree->Branch("muType_lead",        &muType_lead_);
  tree->Branch("muIDbit_lead",       &muIDbit_lead_);
  tree->Branch("muD0_lead",          &muD0_lead_);
  tree->Branch("muDz_lead",          &muDz_lead_);
  tree->Branch("muSIP_lead",         &muSIP_lead_);
  tree->Branch("muD0Error_lead",     &muD0Error_lead_);
  tree->Branch("muDzError_lead",     &muDzError_lead_);
  tree->Branch("muChi2NDF_lead",     &muChi2NDF_lead_);
  tree->Branch("muInnerD0_lead",     &muInnerD0_lead_);
  tree->Branch("muInnerDz_lead",     &muInnerDz_lead_);
  tree->Branch("muTrkLayers_lead",   &muTrkLayers_lead_);
  tree->Branch("muPixelLayers_lead", &muPixelLayers_lead_);
  tree->Branch("muPixelHits_lead",   &muPixelHits_lead_);
  tree->Branch("muMuonHits_lead",    &muMuonHits_lead_);
  tree->Branch("muStations_lead",    &muStations_lead_);
  tree->Branch("muMatches_lead",     &muMatches_lead_);
  tree->Branch("muTrkQuality_lead",  &muTrkQuality_lead_);
  tree->Branch("muIsoTrk_lead",      &muIsoTrk_lead_);
  tree->Branch("muPFChIso_lead",     &muPFChIso_lead_);
  tree->Branch("muPFPhoIso_lead",    &muPFPhoIso_lead_);
  tree->Branch("muPFNeuIso_lead",    &muPFNeuIso_lead_);
  tree->Branch("muPFPUIso_lead",     &muPFPUIso_lead_);
  tree->Branch("muPFChIso03_lead",   &muPFChIso03_lead_);
  tree->Branch("muPFPhoIso03_lead",  &muPFPhoIso03_lead_);
  tree->Branch("muPFNeuIso03_lead",  &muPFNeuIso03_lead_);
  tree->Branch("muPFPUIso03_lead",   &muPFPUIso03_lead_);
//  tree->Branch("muPFMiniIso_lead",   &muPFMiniIso_lead_);
  tree->Branch("muFiredTrgs_lead",   &muFiredTrgs_lead_);
  tree->Branch("muFiredL1Trgs_lead", &muFiredL1Trgs_lead_);
  tree->Branch("muInnervalidFraction_lead",   &muInnervalidFraction_lead_);
  tree->Branch("musegmentCompatibility_lead", &musegmentCompatibility_lead_);
  tree->Branch("muchi2LocalPosition_lead",    &muchi2LocalPosition_lead_);
  tree->Branch("mutrkKink_lead",              &mutrkKink_lead_);
  tree->Branch("muBestTrkPtError_lead",       &muBestTrkPtError_lead_);
  tree->Branch("muBestTrkPt_lead",            &muBestTrkPt_lead_);
  tree->Branch("muBestTrkType_lead",          &muBestTrkType_lead_);
  tree->Branch("muTrkNormChi2_lead",          &muTrkNormChi2_lead_);
  tree->Branch("muIDPatMVA_lead",             &muIDPatMVA_lead_);
  tree->Branch("muIDPatSoftMVA_lead",         &muIDPatSoftMVA_lead_);
  tree->Branch("muIDSelectorBit_lead",             &muIDSelectorBit_lead_);

  tree->Branch("muPt_sublead",          &muPt_sublead_);
  tree->Branch("muEn_sublead",          &muEn_sublead_);
  tree->Branch("muEta_sublead",         &muEta_sublead_);
  tree->Branch("muPhi_sublead",         &muPhi_sublead_);
  tree->Branch("muCharge_sublead",      &muCharge_sublead_);
  tree->Branch("muType_sublead",        &muType_sublead_);
  tree->Branch("muIDbit_sublead",       &muIDbit_sublead_);
  tree->Branch("muD0_sublead",          &muD0_sublead_);
  tree->Branch("muDz_sublead",          &muDz_sublead_);
  tree->Branch("muSIP_sublead",         &muSIP_sublead_);
  tree->Branch("muD0Error_sublead",     &muD0Error_sublead_);
  tree->Branch("muDzError_sublead",     &muDzError_sublead_);
  tree->Branch("muChi2NDF_sublead",     &muChi2NDF_sublead_);
  tree->Branch("muInnerD0_sublead",     &muInnerD0_sublead_);
  tree->Branch("muInnerDz_sublead",     &muInnerDz_sublead_);
  tree->Branch("muTrkLayers_sublead",   &muTrkLayers_sublead_);
  tree->Branch("muPixelLayers_sublead", &muPixelLayers_sublead_);
  tree->Branch("muPixelHits_sublead",   &muPixelHits_sublead_);
  tree->Branch("muMuonHits_sublead",    &muMuonHits_sublead_);
  tree->Branch("muStations_sublead",    &muStations_sublead_);
  tree->Branch("muMatches_sublead",     &muMatches_sublead_);
  tree->Branch("muTrkQuality_sublead",  &muTrkQuality_sublead_);
  tree->Branch("muIsoTrk_sublead",      &muIsoTrk_sublead_);
  tree->Branch("muPFChIso_sublead",     &muPFChIso_sublead_);
  tree->Branch("muPFPhoIso_sublead",    &muPFPhoIso_sublead_);
  tree->Branch("muPFNeuIso_sublead",    &muPFNeuIso_sublead_);
  tree->Branch("muPFPUIso_sublead",     &muPFPUIso_sublead_);
  tree->Branch("muPFChIso03_sublead",   &muPFChIso03_sublead_);
  tree->Branch("muPFPhoIso03_sublead",  &muPFPhoIso03_sublead_);
  tree->Branch("muPFNeuIso03_sublead",  &muPFNeuIso03_sublead_);
  tree->Branch("muPFPUIso03_sublead",   &muPFPUIso03_sublead_);
//  tree->Branch("muPFMiniIso_sublead",   &muPFMiniIso_sublead_);
  tree->Branch("muFiredTrgs_sublead",   &muFiredTrgs_sublead_);
  tree->Branch("muFiredL1Trgs_sublead", &muFiredL1Trgs_sublead_);
  tree->Branch("muInnervalidFraction_sublead",   &muInnervalidFraction_sublead_);
  tree->Branch("musegmentCompatibility_sublead", &musegmentCompatibility_sublead_);
  tree->Branch("muchi2LocalPosition_sublead",    &muchi2LocalPosition_sublead_);
  tree->Branch("mutrkKink_sublead",              &mutrkKink_sublead_);
  tree->Branch("muBestTrkPtError_sublead",       &muBestTrkPtError_sublead_);
  tree->Branch("muBestTrkPt_sublead",            &muBestTrkPt_sublead_);
  tree->Branch("muBestTrkType_sublead",          &muBestTrkType_sublead_);
  tree->Branch("muTrkNormChi2_sublead",          &muTrkNormChi2_sublead_);
  tree->Branch("muIDPatMVA_sublead",             &muIDPatMVA_sublead_);
  tree->Branch("muIDPatSoftMVA_sublead",         &muIDPatSoftMVA_sublead_);
  tree->Branch("muIDSelectorBit_sublead",             &muIDSelectorBit_sublead_);
  

  tree->Branch("muSvChi2",                  &muSvChi2_);
  tree->Branch("muSvNDOF",                  &muSvNDOF_);
  tree->Branch("muSvProb",                  &muSvProb_);
  tree->Branch("muSvX",                     &muSvX_);
  tree->Branch("muSvY",                     &muSvY_);
  tree->Branch("muSvZ",                     &muSvZ_);
  tree->Branch("muSvXError",                &muSvXError_);
  tree->Branch("muSvYError",                &muSvYError_);
  tree->Branch("muSvZError",                &muSvZError_);
  tree->Branch("muSvMass",                  &muSvMass_);
  tree->Branch("muSvCtxy",                  &muSvCtxy_);
  tree->Branch("muSvCosAngle",              &muSvCosAngle_);
  tree->Branch("muSvLxy",                   &muSvLxy_);
  tree->Branch("muSvLxyError",              &muSvLxyError_);
/*
  tree->Branch("kaonMMCharge_lead",               &kaonMMCharge_lead_);
  tree->Branch("kaonMMD0_lead",                   &kaonMMD0_lead_);
  tree->Branch("kaonMMDz_lead",                   &kaonMMDz_lead_);
  tree->Branch("kaonMMD0Error_lead",              &kaonMMD0Error_lead_);
  tree->Branch("kaonMMDzError_lead",              &kaonMMDzError_lead_);
  tree->Branch("kaonMMPt_lead",                   &kaonMMPt_lead_);
  tree->Branch("kaonMMEta_lead",                  &kaonMMEta_lead_);
  tree->Branch("kaonMMPhi_lead",                  &kaonMMPhi_lead_);
  tree->Branch("kaonMMVx_lead",                   &kaonMMVx_lead_);
  tree->Branch("kaonMMVy_lead",                   &kaonMMVy_lead_);
  tree->Branch("kaonMMVz_lead",                   &kaonMMVz_lead_);
  tree->Branch("kaonMMEn_lead",                   &kaonMMEn_lead_);
  tree->Branch("kaonMMTrkChi2_lead",              &kaonMMTrkChi2_lead_);
  tree->Branch("kaonMMTrkNDOF_lead",              &kaonMMTrkNDOF_lead_);
  tree->Branch("kaonMMTrkNormChi2_lead",          &kaonMMTrkNormChi2_lead_);

  tree->Branch("kaonMMCharge_sublead",               &kaonMMCharge_sublead_);
  tree->Branch("kaonMMD0_sublead",                   &kaonMMD0_sublead_);
  tree->Branch("kaonMMDz_sublead",                   &kaonMMDz_sublead_);
  tree->Branch("kaonMMD0Error_sublead",              &kaonMMD0Error_sublead_);
  tree->Branch("kaonMMDzError_sublead",              &kaonMMDzError_sublead_);
  tree->Branch("kaonMMPt_sublead",                   &kaonMMPt_sublead_);
  tree->Branch("kaonMMEta_sublead",                  &kaonMMEta_sublead_);
  tree->Branch("kaonMMPhi_sublead",                  &kaonMMPhi_sublead_);
  tree->Branch("kaonMMVx_sublead",                   &kaonMMVx_sublead_);
  tree->Branch("kaonMMVy_sublead",                   &kaonMMVy_sublead_);
  tree->Branch("kaonMMVz_sublead",                   &kaonMMVz_sublead_);
  tree->Branch("kaonMMEn_sublead",                   &kaonMMEn_sublead_);
  tree->Branch("kaonMMTrkChi2_sublead",              &kaonMMTrkChi2_sublead_);
  tree->Branch("kaonMMTrkNDOF_sublead",              &kaonMMTrkNDOF_sublead_);
  tree->Branch("kaonMMTrkNormChi2_sublead",          &kaonMMTrkNormChi2_sublead_);

  tree->Branch("bsMMdRmu",                     &bsMMdRmu_);
  tree->Branch("bsMMdRkaon",                   &bsMMdRkaon_);
  tree->Branch("bsMMdRJpsiPhi",                &bsMMdRJpsiPhi_);
  tree->Branch("bsMMJpsiMass",                 &bsMMJpsiMass_);
  tree->Branch("bsMMPhiMass",                  &bsMMPhiMass_);
  tree->Branch("bsMMBsMass",                   &bsMMBsMass_);
*/

}

void ggNtuplizer::fillMuons(const edm::Event& e, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  muPt_lead_                  .clear();
  muEn_lead_                  .clear();
  muEta_lead_                 .clear();
  muPhi_lead_                 .clear();
  muCharge_lead_              .clear();
  muType_lead_                .clear();
  muIDbit_lead_               .clear();
  muD0_lead_                  .clear();
  muDz_lead_                  .clear();
  muSIP_lead_                 .clear();
  muD0Error_lead_             .clear();
  muDzError_lead_             .clear();
  muChi2NDF_lead_             .clear();
  muInnerD0_lead_             .clear();
  muInnerDz_lead_             .clear();
  muTrkLayers_lead_           .clear();
  muPixelLayers_lead_         .clear();
  muPixelHits_lead_           .clear();
  muMuonHits_lead_            .clear();
  muStations_lead_            .clear();
  muMatches_lead_             .clear();
  muTrkQuality_lead_          .clear();
  muIsoTrk_lead_              .clear();
  muPFChIso_lead_             .clear();
  muPFPhoIso_lead_            .clear();
  muPFNeuIso_lead_            .clear();
  muPFPUIso_lead_             .clear();
  muPFChIso03_lead_           .clear();
  muPFPhoIso03_lead_          .clear();
  muPFNeuIso03_lead_          .clear();
  muPFPUIso03_lead_           .clear();
//  muPFMiniIso_lead_           .clear();
  muFiredTrgs_lead_           .clear();
  muFiredL1Trgs_lead_         .clear();
  muInnervalidFraction_lead_  .clear();
  musegmentCompatibility_lead_.clear();
  muchi2LocalPosition_lead_   .clear();
  mutrkKink_lead_             .clear();
  muBestTrkPtError_lead_      .clear();
  muBestTrkPt_lead_           .clear();
  muBestTrkType_lead_         .clear();
  muTrkNormChi2_lead_         .clear();
  muIDPatMVA_lead_            .clear();
  muIDPatSoftMVA_lead_        .clear();
  muIDSelectorBit_lead_            .clear();

  muPt_sublead_                  .clear();
  muEn_sublead_                  .clear();
  muEta_sublead_                 .clear();
  muPhi_sublead_                 .clear();
  muCharge_sublead_              .clear();
  muType_sublead_                .clear();
  muIDbit_sublead_               .clear();
  muD0_sublead_                  .clear();
  muDz_sublead_                  .clear();
  muSIP_sublead_                 .clear();
  muD0Error_sublead_             .clear();
  muDzError_sublead_             .clear();
  muChi2NDF_sublead_             .clear();
  muInnerD0_sublead_             .clear();
  muInnerDz_sublead_             .clear();
  muTrkLayers_sublead_           .clear();
  muPixelLayers_sublead_         .clear();
  muPixelHits_sublead_           .clear();
  muMuonHits_sublead_            .clear();
  muStations_sublead_            .clear();
  muMatches_sublead_             .clear();
  muTrkQuality_sublead_          .clear();
  muIsoTrk_sublead_              .clear();
  muPFChIso_sublead_             .clear();
  muPFPhoIso_sublead_            .clear();
  muPFNeuIso_sublead_            .clear();
  muPFPUIso_sublead_             .clear();
  muPFChIso03_sublead_           .clear();
  muPFPhoIso03_sublead_          .clear();
  muPFNeuIso03_sublead_          .clear();
  muPFPUIso03_sublead_           .clear();
//  muPFMiniIso_sublead_           .clear();
  muFiredTrgs_sublead_           .clear();
  muFiredL1Trgs_sublead_         .clear();
  muInnervalidFraction_sublead_  .clear();
  musegmentCompatibility_sublead_.clear();
  muchi2LocalPosition_sublead_   .clear();
  mutrkKink_sublead_             .clear();
  muBestTrkPtError_sublead_      .clear();
  muBestTrkPt_sublead_           .clear();
  muBestTrkType_sublead_         .clear();
  muTrkNormChi2_sublead_         .clear();
  muIDPatMVA_sublead_            .clear();
  muIDPatSoftMVA_sublead_        .clear();
  muIDSelectorBit_sublead_            .clear();

  muSvChi2_.clear();
  muSvNDOF_.clear();
  muSvProb_.clear();
  muSvX_.clear();
  muSvY_.clear();
  muSvZ_.clear();
  muSvXError_.clear();
  muSvYError_.clear();
  muSvZError_.clear();
  muSvMass_.clear();
  muSvCtxy_.clear();
  muSvCosAngle_.clear();
  muSvLxy_.clear();
  muSvLxyError_.clear();

/*  kaonMMCharge_lead_                  .clear();
  kaonMMD0_lead_                      .clear();
  kaonMMDz_lead_                      .clear();
  kaonMMD0Error_lead_                 .clear();
  kaonMMDzError_lead_                 .clear();
  kaonMMPt_lead_                      .clear();
  kaonMMEta_lead_                     .clear();
  kaonMMPhi_lead_                     .clear();
  kaonMMVx_lead_                      .clear();
  kaonMMVy_lead_                      .clear();
  kaonMMVz_lead_                      .clear();
  kaonMMEn_lead_	                  .clear();
  kaonMMTrkChi2_lead_                 .clear();
  kaonMMTrkNDOF_lead_                 .clear();
  kaonMMTrkNormChi2_lead_             .clear();

  kaonMMCharge_sublead_                  .clear();
  kaonMMD0_sublead_                      .clear();
  kaonMMDz_sublead_                      .clear();
  kaonMMD0Error_sublead_                 .clear();
  kaonMMDzError_sublead_                 .clear();
  kaonMMPt_sublead_                      .clear();
  kaonMMEta_sublead_                     .clear();
  kaonMMPhi_sublead_                     .clear();
  kaonMMVx_sublead_                      .clear();
  kaonMMVy_sublead_                      .clear();
  kaonMMVz_sublead_                      .clear();
  kaonMMEn_sublead_	                     .clear();
  kaonMMTrkChi2_sublead_                 .clear();
  kaonMMTrkNDOF_sublead_                 .clear();
  kaonMMTrkNormChi2_sublead_             .clear();

  bsMMdRmu_               .clear();
  bsMMdRkaon_             .clear();
  bsMMdRJpsiPhi_          .clear();
  bsMMJpsiMass_		      .clear();
  bsMMPhiMass_		      .clear();
  bsMMBsMass_		      .clear();
*/
  nMu_ = 0;
  matchedTrg_ = false;

  if (isAOD_) {

    edm::Handle<edm::View<pat::Muon> > muonHandle;
    e.getByToken(muonCollection_, muonHandle);

  //  edm::Handle<pat::PackedCandidateCollection> pfcands;
  //  e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<reco::TrackCollection> tracksHandle;
    e.getByToken(tracklabel_, tracksHandle);

    if (!muonHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
      return;
    }

    vector<bool> muTrkMap;
    muTrkMap = muonTriggerMap(e);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
     // std::cout<<"first muon entered"<<std::endl;
      if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) == 1) matchedTrg_ = true;
      //if (iMu->pt() < 2) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
      //if (fabs(iMu->vz() - pv.z()) > 1.0) continue;

      for (edm::View<pat::Muon>::const_iterator jMu = iMu+1; jMu != muonHandle->end(); ++jMu) {
        if (iMu->charge()*jMu->charge() > 0 ) continue;
        //if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) != 1 || matchMuonTriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()) != 1) continue;
        if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;

        float pmass  = 0.1056583745;
        TLorentzVector iMu_lv, jMu_lv, Z_lv;
        iMu_lv.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), pmass);
        jMu_lv.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), pmass);      

        auto leadMu = iMu->pt() > jMu->pt() ? iMu : jMu;
        auto subleadMu = iMu->pt() > jMu->pt() ? jMu : iMu;


        muPt_lead_    .push_back(leadMu->pt());
        muEn_lead_    .push_back(leadMu->energy());
        muEta_lead_   .push_back(leadMu->eta());
        muPhi_lead_   .push_back(leadMu->phi());
        muCharge_lead_.push_back(leadMu->charge());
        muType_lead_  .push_back(leadMu->type());
        muD0_lead_    .push_back(leadMu->muonBestTrack()->dxy(pv));
        muDz_lead_    .push_back(leadMu->muonBestTrack()->dz(pv));
        muSIP_lead_   .push_back(leadMu->dB(pat::Muon::PV3D)/leadMu->edB(pat::Muon::PV3D));
        muD0Error_lead_    .push_back(leadMu->muonBestTrack()->dxyError());
        muDzError_lead_    .push_back(leadMu->muonBestTrack()->dzError());
        muTrkNormChi2_lead_.push_back(leadMu->muonBestTrack()->normalizedChi2());

        UShort_t tmpmuIDbit = 0;

        if (leadMu->isLooseMuon())     setbit(tmpmuIDbit, 0);
        if (leadMu->isMediumMuon())    setbit(tmpmuIDbit, 1);
        if (leadMu->isTightMuon(vtx))  setbit(tmpmuIDbit, 2);
        if (leadMu->isSoftMuon(vtx))   setbit(tmpmuIDbit, 3);
        if (leadMu->isHighPtMuon(vtx)) setbit(tmpmuIDbit, 4);

        muIDbit_lead_.push_back(tmpmuIDbit);
        muIDPatMVA_lead_.push_back(leadMu->mvaValue());
        muIDPatSoftMVA_lead_.push_back(leadMu->softMvaValue());
        // muon MVA is only available in MINIAOD
        muIDSelectorBit_lead_.push_back(leadMu->selectors());

        muFiredTrgs_lead_.push_back(muTrkMap[leadMu - muonHandle->begin()]);
        muFiredL1Trgs_lead_.push_back(matchL1TriggerFilters(leadMu->pt(), leadMu->eta(), leadMu->phi()));

        muBestTrkPtError_lead_        .push_back(leadMu->muonBestTrack()->ptError());
        muBestTrkPt_lead_             .push_back(leadMu->muonBestTrack()->pt());
        muBestTrkType_lead_           .push_back(leadMu->muonBestTrackType());
        musegmentCompatibility_lead_  .push_back(leadMu->segmentCompatibility());
        muchi2LocalPosition_lead_     .push_back(leadMu->combinedQuality().chi2LocalPosition);
        mutrkKink_lead_               .push_back(leadMu->combinedQuality().trkKink);

        reco::TrackRef glbmu = leadMu->globalTrack();
        reco::TrackRef innmu = leadMu->innerTrack();

        if (glbmu.isNull()) {
        muChi2NDF_lead_ .push_back(-99.);
        muMuonHits_lead_.push_back(-99);
        } else {
        muChi2NDF_lead_.push_back(glbmu->normalizedChi2());
        muMuonHits_lead_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
        }

        if (innmu.isNull()) {
        muInnerD0_lead_     .push_back(-99.);
        muInnerDz_lead_     .push_back(-99);
        muTrkLayers_lead_   .push_back(-99);
        muPixelLayers_lead_ .push_back(-99);
        muPixelHits_lead_   .push_back(-99);
        muTrkQuality_lead_  .push_back(-99);

        muInnervalidFraction_lead_ .push_back(-99);
        } else {
        muInnerD0_lead_     .push_back(innmu->dxy(pv));
        muInnerDz_lead_     .push_back(innmu->dz(pv));
        muTrkLayers_lead_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
        muPixelLayers_lead_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
        muPixelHits_lead_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
        muTrkQuality_lead_  .push_back(innmu->quality(reco::TrackBase::highPurity));

        muInnervalidFraction_lead_ .push_back(innmu->validFraction());
        }

        muStations_lead_   .push_back(leadMu->numberOfMatchedStations());
        muMatches_lead_    .push_back(leadMu->numberOfMatches());
        muIsoTrk_lead_     .push_back(leadMu->trackIso());
        muPFChIso_lead_    .push_back(leadMu->pfIsolationR04().sumChargedHadronPt);
        muPFPhoIso_lead_   .push_back(leadMu->pfIsolationR04().sumPhotonEt);
        muPFNeuIso_lead_   .push_back(leadMu->pfIsolationR04().sumNeutralHadronEt);
        muPFPUIso_lead_    .push_back(leadMu->pfIsolationR04().sumPUPt);
        muPFChIso03_lead_  .push_back(leadMu->pfIsolationR03().sumChargedHadronPt);
        muPFPhoIso03_lead_ .push_back(leadMu->pfIsolationR03().sumPhotonEt);
        muPFNeuIso03_lead_ .push_back(leadMu->pfIsolationR03().sumNeutralHadronEt);
        muPFPUIso03_lead_  .push_back(leadMu->pfIsolationR03().sumPUPt);

        muPt_sublead_    .push_back(subleadMu->pt());
        muEn_sublead_    .push_back(subleadMu->energy());
        muEta_sublead_   .push_back(subleadMu->eta());
        muPhi_sublead_   .push_back(subleadMu->phi());
        muCharge_sublead_.push_back(subleadMu->charge());
        muType_sublead_  .push_back(subleadMu->type());
        muD0_sublead_    .push_back(subleadMu->muonBestTrack()->dxy(pv));
        muDz_sublead_    .push_back(subleadMu->muonBestTrack()->dz(pv));
        muSIP_sublead_   .push_back(subleadMu->dB(pat::Muon::PV3D)/subleadMu->edB(pat::Muon::PV3D));
        muD0Error_sublead_    .push_back(subleadMu->muonBestTrack()->dxyError());
        muDzError_sublead_    .push_back(subleadMu->muonBestTrack()->dzError());
        muTrkNormChi2_sublead_.push_back(subleadMu->muonBestTrack()->normalizedChi2());

        tmpmuIDbit = 0;

        if (subleadMu->isLooseMuon())     setbit(tmpmuIDbit, 0);
        if (subleadMu->isMediumMuon())    setbit(tmpmuIDbit, 1);
        if (subleadMu->isTightMuon(vtx))  setbit(tmpmuIDbit, 2);
        if (subleadMu->isSoftMuon(vtx))   setbit(tmpmuIDbit, 3);
        if (subleadMu->isHighPtMuon(vtx)) setbit(tmpmuIDbit, 4);

        muIDbit_sublead_.push_back(tmpmuIDbit);
        muIDPatMVA_sublead_.push_back(subleadMu->mvaValue());
        muIDPatSoftMVA_sublead_.push_back(subleadMu->softMvaValue());
        // muon MVA is only available in MINIAOD
        muIDSelectorBit_sublead_.push_back(subleadMu->selectors());

        muFiredTrgs_sublead_.push_back(muTrkMap[subleadMu - muonHandle->begin()]);
        muFiredL1Trgs_sublead_.push_back(matchL1TriggerFilters(subleadMu->pt(), subleadMu->eta(), subleadMu->phi()));

        muBestTrkPtError_sublead_        .push_back(subleadMu->muonBestTrack()->ptError());
        muBestTrkPt_sublead_             .push_back(subleadMu->muonBestTrack()->pt());
        muBestTrkType_sublead_           .push_back(subleadMu->muonBestTrackType());
        musegmentCompatibility_sublead_  .push_back(subleadMu->segmentCompatibility());
        muchi2LocalPosition_sublead_     .push_back(subleadMu->combinedQuality().chi2LocalPosition);
        mutrkKink_sublead_               .push_back(subleadMu->combinedQuality().trkKink);

        glbmu = subleadMu->globalTrack();
        innmu = subleadMu->innerTrack();

        if (glbmu.isNull()) {
        muChi2NDF_sublead_ .push_back(-99.);
        muMuonHits_sublead_.push_back(-99);
        } else {
        muChi2NDF_sublead_.push_back(glbmu->normalizedChi2());
        muMuonHits_sublead_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
        }

        if (innmu.isNull()) {
        muInnerD0_sublead_     .push_back(-99.);
        muInnerDz_sublead_     .push_back(-99);
        muTrkLayers_sublead_   .push_back(-99);
        muPixelLayers_sublead_ .push_back(-99);
        muPixelHits_sublead_   .push_back(-99);
        muTrkQuality_sublead_  .push_back(-99);

        muInnervalidFraction_sublead_ .push_back(-99);
        } else {
        muInnerD0_sublead_     .push_back(innmu->dxy(pv));
        muInnerDz_sublead_     .push_back(innmu->dz(pv));
        muTrkLayers_sublead_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
        muPixelLayers_sublead_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
        muPixelHits_sublead_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
        muTrkQuality_sublead_  .push_back(innmu->quality(reco::TrackBase::highPurity));

        muInnervalidFraction_sublead_ .push_back(innmu->validFraction());
        }

        muStations_sublead_   .push_back(subleadMu->numberOfMatchedStations());
        muMatches_sublead_    .push_back(subleadMu->numberOfMatches());
        muIsoTrk_sublead_     .push_back(subleadMu->trackIso());
        muPFChIso_sublead_    .push_back(subleadMu->pfIsolationR04().sumChargedHadronPt);
        muPFPhoIso_sublead_   .push_back(subleadMu->pfIsolationR04().sumPhotonEt);
        muPFNeuIso_sublead_   .push_back(subleadMu->pfIsolationR04().sumNeutralHadronEt);
        muPFPUIso_sublead_    .push_back(subleadMu->pfIsolationR04().sumPUPt);
        muPFChIso03_sublead_  .push_back(subleadMu->pfIsolationR03().sumChargedHadronPt);
        muPFPhoIso03_sublead_ .push_back(subleadMu->pfIsolationR03().sumPhotonEt);
        muPFNeuIso03_sublead_ .push_back(subleadMu->pfIsolationR03().sumNeutralHadronEt);
        muPFPUIso03_sublead_  .push_back(subleadMu->pfIsolationR03().sumPUPt);


        nMu_++;
      }
    }
  } else {

    edm::Handle<edm::View<pat::Muon> > muonHandle;
    e.getByToken(muonCollection_, muonHandle);

    edm::Handle<pat::PackedCandidateCollection> pfcands;
    e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<pat::PackedCandidateCollection> losttracks;
    //e.getByLabel("lostTracks", losttracks);
    e.getByToken(lostTracksLabel_, losttracks);

    std::vector<pat::PackedCandidate> alltracks;
    //edm::Handle<pat::PackedCandidateCollection> alltracks;
    alltracks.reserve(pfcands->size() + losttracks->size());
    alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
    alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

    //pfcands->reserve(pfcands->size() + losttracks->size());
    //pfcands->insert(pfcands->end(), losttracks->begin(), losttracks->end());
    //std::copy(pfcands->begin(), pfcands->end(), std::back_inserter(alltracks));
    //std::copy(losttracks->begin(), losttracks->end(), std::back_inserter(alltracks));

    if (!muonHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
      return;
    }

    vector<bool> muTrkMap;
    muTrkMap = muonTriggerMap(e);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) == 1) matchedTrg_ = true;
      if (iMu->pt() < 0.8) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
      if (fabs(iMu->eta()) > 2.5) continue;
      if (fabs(iMu->vz() - pv.z()) > 0.5) continue;
      for (edm::View<pat::Muon>::const_iterator jMu = iMu+1; jMu != muonHandle->end(); ++jMu) {
        //if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) != 1 || matchMuonTriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()) != 1) continue;
        if (jMu->pt() < 0.8) continue;
        if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;
        if (fabs(jMu->eta()) > 2.5) continue;
        if (fabs(jMu->vz() - pv.z()) > 0.5) continue;

        float pmass  = 0.1056583745;
        TLorentzVector iMu_lv, jMu_lv, Jpsi_lv;
        iMu_lv.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), pmass);
        jMu_lv.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), pmass);      
        
        Jpsi_lv = iMu_lv + jMu_lv;
        if (iMu->charge()*jMu->charge() > 0) continue;

        if (Jpsi_lv.M() < 2.4 || Jpsi_lv.M() > 3.8) continue;

        KinematicParticleFactoryFromTransientTrack pFactory;  
        std::vector<RefCountedKinematicParticle> XParticles;
        float pmasse = 1.e-6 * pmass;
        const reco::TransientTrack imuttk = getTransientTrack( *(iMu->bestTrack()) );
        const reco::TransientTrack jmuttk = getTransientTrack( *(jMu->bestTrack()) );


        //XParticles.push_back(pFactory.particle(getTransientTrack( *(iMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        //XParticles.push_back(pFactory.particle(getTransientTrack( *(jMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        XParticles.push_back(pFactory.particle(imuttk, pmass, 0.0, 0, pmasse));
        XParticles.push_back(pFactory.particle(jmuttk, pmass, 0.0, 0, pmasse));

        KinematicConstrainedVertexFitter kvFitter;
        RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

        if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 || KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
        //if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
        KinVtx->movePointerToTheTop();
        RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();
        RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();

        //for (pat::PackedCandidateCollection::const_iterator iHad = pfcands->begin(); iHad != pfcands->end(); ++iHad) {
//        for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {

          //if (iHad->pt() < 0.4) continue;
          //if (iHad->charge() == 0) continue;
          //if (abs(iHad->pdgId()) != 211) continue;
          //if (iHad->bestTrack() == nullptr) continue;
          //if (fabs(iHad->eta()) > 2.5) continue;
          //if (fabs(iMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
          //if (fabs(jMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(jMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(jMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
          //if (fabs(iMu->bestTrack()->vz() - iHad->vz()) > 1) continue;
          //if (iHad->normalizedChi2() < 0.0) continue;
          //if (iHad->normalizedChi2() > 20.0) continue;

          //for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != pfcands->end(); ++jHad) {
//          for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != alltracks.end(); ++jHad) {
/*
            if (jHad->pt() < 0.4) continue;
            if (jHad->charge() == 0) continue;
            if (abs(jHad->pdgId()) != 211) continue;
            if (jHad->bestTrack() == nullptr) continue;
            //if (iHad->charge()*jHad->charge() > 0.0) continue;
            //if (fabs(iMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
            //if (fabs(jMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(jMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(jMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
            //if (fabs(iMu->bestTrack()->vz() - jHad->vz()) > 1) continue;
            if (fabs(jHad->vz() - pv.z()) > 1.0) continue;
            //if (jHad->normalizedChi2() < 0.0) continue;
            //if (jHad->normalizedChi2() > 20.0) continue;

            // Phi mass window
            float kpmass = 0.493677;
            float bsM = 5.3663;

            TLorentzVector iHad_lv, jHad_lv, bs_lv;
            iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
            jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);     
            bs_lv = iMu_lv + jMu_lv + iHad_lv + jHad_lv;
            //if (((iHad_lv+jHad_lv)).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.06) continue; 
            if ((iHad_lv+jHad_lv).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.10) continue; 
            if ((iMu_lv + jMu_lv + iHad_lv + jHad_lv).M() < 4.0 || (iMu_lv + jMu_lv + iHad_lv + jHad_lv).M() > 6.0) continue;
            if (fabs(jHad->eta()) > 2.5) continue;

            std::vector<RefCountedKinematicParticle> BsParticles;
            float kpmasse = 1.e-6 * pmass;

            BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
            BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
            BsParticles.push_back(pFactory.particle(imuttk, pmass, 0.0, 0, pmasse));
            BsParticles.push_back(pFactory.particle(jmuttk, pmass, 0.0, 0, pmasse));

            KinematicConstrainedVertexFitter BsKvFitter;
            RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
            if (!(BsKinVtx->isValid())) continue;

            RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

            if (DecayVtx->chiSquared() < 0.0) continue;
            //if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 20.0) continue;

            // Accept these 4 tracks as a Bs candidate, fill ntuple
*/
            auto leadMu = iMu->pt() > jMu->pt() ? iMu : jMu;
            auto subleadMu = iMu->pt() > jMu->pt() ? jMu : iMu;
            //auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
            //auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

            double ctxy = ((DecayVtx->position().x() - pv.x())*Jpsi_lv.Px() + (DecayVtx->position().y() - pv.y())*Jpsi_lv.Py())/(pow(Jpsi_lv.Pt(),2))*Jpsi_lv.M();
            
            math::XYZVector perp(Jpsi_lv.Px(), Jpsi_lv.Py(), 0.);
            math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
            math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
            double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

            muSvChi2_.push_back(DecayVtx->chiSquared());
            muSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
            muSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
            muSvX_.push_back(DecayVtx->position().x());
            muSvY_.push_back(DecayVtx->position().y());
            muSvZ_.push_back(DecayVtx->position().z());
            muSvXError_.push_back(DecayVtx->error().cxx());
            muSvYError_.push_back(DecayVtx->error().cyy());
            muSvZError_.push_back(DecayVtx->error().czz());
            muSvMass_.push_back(Jpsi_lv.M());
            muSvCtxy_.push_back(ctxy);
            muSvCosAngle_.push_back(cosAngle);
            muSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
            muSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

/*            kaonMMCharge_lead_          .push_back(leadHad->charge());
            kaonMMD0_lead_              .push_back(leadHad->dxy(pv));
            kaonMMDz_lead_              .push_back(leadHad->dz(pv));
            kaonMMD0Error_lead_		  .push_back(leadHad->dxyError());
            kaonMMDzError_lead_		  .push_back(leadHad->dzError());
            kaonMMPt_lead_              .push_back(leadHad->pt());
            kaonMMEta_lead_             .push_back(leadHad->eta());
            kaonMMPhi_lead_             .push_back(leadHad->phi());
            kaonMMVx_lead_		      .push_back(leadHad->vx());
            kaonMMVy_lead_		      .push_back(leadHad->vy());
            kaonMMVz_lead_		      .push_back(leadHad->vz());
  //        kaonMMEn_lead_ 	          .push_back(leadHad->energy());
            kaonMMTrkChi2_lead_		  .push_back(leadHad->bestTrack()->chi2());
            kaonMMTrkNDOF_lead_		  .push_back(leadHad->bestTrack()->ndof());
            kaonMMTrkNormChi2_lead_	  .push_back(leadHad->bestTrack()->normalizedChi2());

            kaonMMCharge_sublead_       .push_back(subleadHad->charge());
            kaonMMD0_sublead_           .push_back(subleadHad->dxy(pv));
            kaonMMDz_sublead_           .push_back(subleadHad->dz(pv));
            kaonMMD0Error_sublead_      .push_back(subleadHad->dxyError());
            kaonMMDzError_sublead_	  .push_back(subleadHad->dzError());
            kaonMMPt_sublead_           .push_back(subleadHad->pt());
            kaonMMEta_sublead_          .push_back(subleadHad->eta());
            kaonMMPhi_sublead_          .push_back(subleadHad->phi());
            kaonMMVx_sublead_		      .push_back(subleadHad->vx());
            kaonMMVy_sublead_		      .push_back(subleadHad->vy());
            kaonMMVz_sublead_		      .push_back(subleadHad->vz());
  //        kaonMMEn_sublead_ 	      .push_back(subleadHad->energy());
            kaonMMTrkChi2_sublead_		  .push_back(subleadHad->bestTrack()->chi2());
            kaonMMTrkNDOF_sublead_		  .push_back(subleadHad->bestTrack()->ndof());
            kaonMMTrkNormChi2_sublead_	  .push_back(subleadHad->bestTrack()->normalizedChi2());

            bsMMdRmu_            .push_back(iMu_lv.DeltaR(jMu_lv));
            bsMMdRkaon_          .push_back(iHad_lv.DeltaR(jHad_lv));
            bsMMdRJpsiPhi_       .push_back((iMu_lv + jMu_lv).DeltaR(iHad_lv + jHad_lv));
            bsMMJpsiMass_        .push_back((iMu_lv+jMu_lv).M());
            bsMMPhiMass_         .push_back((iHad_lv+jHad_lv).M());
            bsMMBsMass_          .push_back(bs_lv.M());
*/
            muPt_lead_    .push_back(leadMu->pt());
            muEn_lead_    .push_back(leadMu->energy());
            muEta_lead_   .push_back(leadMu->eta());
            muPhi_lead_   .push_back(leadMu->phi());
            muCharge_lead_.push_back(leadMu->charge());
            muType_lead_  .push_back(leadMu->type());
            muD0_lead_    .push_back(leadMu->muonBestTrack()->dxy(pv));
            muDz_lead_    .push_back(leadMu->muonBestTrack()->dz(pv));
            muSIP_lead_   .push_back(leadMu->dB(pat::Muon::PV3D)/leadMu->edB(pat::Muon::PV3D));
            muD0Error_lead_    .push_back(leadMu->muonBestTrack()->dxyError());
            muDzError_lead_    .push_back(leadMu->muonBestTrack()->dzError());
            muTrkNormChi2_lead_.push_back(leadMu->muonBestTrack()->normalizedChi2());

            UShort_t tmpmuIDbit = 0;

            if (leadMu->isLooseMuon())     setbit(tmpmuIDbit, 0);
            if (leadMu->isMediumMuon())    setbit(tmpmuIDbit, 1);
            if (leadMu->isTightMuon(vtx))  setbit(tmpmuIDbit, 2);
            if (leadMu->isSoftMuon(vtx))   setbit(tmpmuIDbit, 3);
            if (leadMu->isHighPtMuon(vtx)) setbit(tmpmuIDbit, 4);


            muIDbit_lead_.push_back(tmpmuIDbit);
            muIDPatMVA_lead_.push_back(leadMu->mvaValue());
            muIDPatSoftMVA_lead_.push_back(leadMu->softMvaValue());
            // muon MVA is only available in MINIAOD
            muIDSelectorBit_lead_.push_back(leadMu->selectors());

            muFiredTrgs_lead_.push_back(muTrkMap[leadMu - muonHandle->begin()]);
            muFiredL1Trgs_lead_.push_back(matchL1TriggerFilters(leadMu->pt(), leadMu->eta(), leadMu->phi()));

            muBestTrkPtError_lead_        .push_back(leadMu->muonBestTrack()->ptError());
            muBestTrkPt_lead_             .push_back(leadMu->muonBestTrack()->pt());
            muBestTrkType_lead_           .push_back(leadMu->muonBestTrackType());
            musegmentCompatibility_lead_  .push_back(leadMu->segmentCompatibility());
            muchi2LocalPosition_lead_     .push_back(leadMu->combinedQuality().chi2LocalPosition);
            mutrkKink_lead_               .push_back(leadMu->combinedQuality().trkKink);

            reco::TrackRef glbmu = leadMu->globalTrack();
            reco::TrackRef innmu = leadMu->innerTrack();

            if (glbmu.isNull()) {
              muChi2NDF_lead_ .push_back(-99.);
              muMuonHits_lead_.push_back(-99);
            } else {
              muChi2NDF_lead_.push_back(glbmu->normalizedChi2());
              muMuonHits_lead_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
            }

            if (innmu.isNull()) {
              muInnerD0_lead_     .push_back(-99.);
              muInnerDz_lead_     .push_back(-99);
              muTrkLayers_lead_   .push_back(-99);
              muPixelLayers_lead_ .push_back(-99);
              muPixelHits_lead_   .push_back(-99);
              muTrkQuality_lead_  .push_back(-99);

              muInnervalidFraction_lead_ .push_back(-99);
            } else {
              muInnerD0_lead_     .push_back(innmu->dxy(pv));
              muInnerDz_lead_     .push_back(innmu->dz(pv));
              muTrkLayers_lead_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
              muPixelLayers_lead_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
              muPixelHits_lead_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
              muTrkQuality_lead_  .push_back(innmu->quality(reco::TrackBase::highPurity));

              muInnervalidFraction_lead_ .push_back(innmu->validFraction());
            }

            muStations_lead_   .push_back(leadMu->numberOfMatchedStations());
            muMatches_lead_    .push_back(leadMu->numberOfMatches());
            muIsoTrk_lead_     .push_back(leadMu->trackIso());
            muPFChIso_lead_    .push_back(leadMu->pfIsolationR04().sumChargedHadronPt);
            muPFPhoIso_lead_   .push_back(leadMu->pfIsolationR04().sumPhotonEt);
            muPFNeuIso_lead_   .push_back(leadMu->pfIsolationR04().sumNeutralHadronEt);
            muPFPUIso_lead_    .push_back(leadMu->pfIsolationR04().sumPUPt);
            muPFChIso03_lead_  .push_back(leadMu->pfIsolationR03().sumChargedHadronPt);
            muPFPhoIso03_lead_ .push_back(leadMu->pfIsolationR03().sumPhotonEt);
            muPFNeuIso03_lead_ .push_back(leadMu->pfIsolationR03().sumNeutralHadronEt);
            muPFPUIso03_lead_  .push_back(leadMu->pfIsolationR03().sumPUPt);

            muPt_sublead_    .push_back(subleadMu->pt());
            muEn_sublead_    .push_back(subleadMu->energy());
            muEta_sublead_   .push_back(subleadMu->eta());
            muPhi_sublead_   .push_back(subleadMu->phi());
            muCharge_sublead_.push_back(subleadMu->charge());
            muType_sublead_  .push_back(subleadMu->type());
            muD0_sublead_    .push_back(subleadMu->muonBestTrack()->dxy(pv));
            muDz_sublead_    .push_back(subleadMu->muonBestTrack()->dz(pv));
            muSIP_sublead_   .push_back(subleadMu->dB(pat::Muon::PV3D)/subleadMu->edB(pat::Muon::PV3D));
            muD0Error_sublead_    .push_back(subleadMu->muonBestTrack()->dxyError());
            muDzError_sublead_    .push_back(subleadMu->muonBestTrack()->dzError());
            muTrkNormChi2_sublead_.push_back(subleadMu->muonBestTrack()->normalizedChi2());

            tmpmuIDbit = 0;

            if (subleadMu->isLooseMuon())     setbit(tmpmuIDbit, 0);
            if (subleadMu->isMediumMuon())    setbit(tmpmuIDbit, 1);
            if (subleadMu->isTightMuon(vtx))  setbit(tmpmuIDbit, 2);
            if (subleadMu->isSoftMuon(vtx))   setbit(tmpmuIDbit, 3);
            if (subleadMu->isHighPtMuon(vtx)) setbit(tmpmuIDbit, 4);

            muIDbit_sublead_.push_back(tmpmuIDbit);
            muIDPatMVA_sublead_.push_back(subleadMu->mvaValue());
            muIDPatSoftMVA_sublead_.push_back(subleadMu->softMvaValue());
            // muon MVA is only available in MINIAOD
            muIDSelectorBit_sublead_.push_back(subleadMu->selectors());

            muFiredTrgs_sublead_.push_back(muTrkMap[subleadMu - muonHandle->begin()]);
            muFiredL1Trgs_sublead_.push_back(matchL1TriggerFilters(subleadMu->pt(), subleadMu->eta(), subleadMu->phi()));

            muBestTrkPtError_sublead_        .push_back(subleadMu->muonBestTrack()->ptError());
            muBestTrkPt_sublead_             .push_back(subleadMu->muonBestTrack()->pt());
            muBestTrkType_sublead_           .push_back(subleadMu->muonBestTrackType());
            musegmentCompatibility_sublead_  .push_back(subleadMu->segmentCompatibility());
            muchi2LocalPosition_sublead_     .push_back(subleadMu->combinedQuality().chi2LocalPosition);
            mutrkKink_sublead_               .push_back(subleadMu->combinedQuality().trkKink);

            glbmu = subleadMu->globalTrack();
            innmu = subleadMu->innerTrack();

            if (glbmu.isNull()) {
              muChi2NDF_sublead_ .push_back(-99.);
              muMuonHits_sublead_.push_back(-99);
            } else {
              muChi2NDF_sublead_.push_back(glbmu->normalizedChi2());
              muMuonHits_sublead_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
            }

            if (innmu.isNull()) {
              muInnerD0_sublead_     .push_back(-99.);
              muInnerDz_sublead_     .push_back(-99);
              muTrkLayers_sublead_   .push_back(-99);
              muPixelLayers_sublead_ .push_back(-99);
              muPixelHits_sublead_   .push_back(-99);
              muTrkQuality_sublead_  .push_back(-99);

              muInnervalidFraction_sublead_ .push_back(-99);
            } else {
              muInnerD0_sublead_     .push_back(innmu->dxy(pv));
              muInnerDz_sublead_     .push_back(innmu->dz(pv));
              muTrkLayers_sublead_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
              muPixelLayers_sublead_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
              muPixelHits_sublead_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
              muTrkQuality_sublead_  .push_back(innmu->quality(reco::TrackBase::highPurity));

              muInnervalidFraction_sublead_ .push_back(innmu->validFraction());
            }

            muStations_sublead_   .push_back(subleadMu->numberOfMatchedStations());
            muMatches_sublead_    .push_back(subleadMu->numberOfMatches());
            muIsoTrk_sublead_     .push_back(subleadMu->trackIso());
            muPFChIso_sublead_    .push_back(subleadMu->pfIsolationR04().sumChargedHadronPt);
            muPFPhoIso_sublead_   .push_back(subleadMu->pfIsolationR04().sumPhotonEt);
            muPFNeuIso_sublead_   .push_back(subleadMu->pfIsolationR04().sumNeutralHadronEt);
            muPFPUIso_sublead_    .push_back(subleadMu->pfIsolationR04().sumPUPt);
            muPFChIso03_sublead_  .push_back(subleadMu->pfIsolationR03().sumChargedHadronPt);
            muPFPhoIso03_sublead_ .push_back(subleadMu->pfIsolationR03().sumPhotonEt);
            muPFNeuIso03_sublead_ .push_back(subleadMu->pfIsolationR03().sumNeutralHadronEt);
            muPFPUIso03_sublead_  .push_back(subleadMu->pfIsolationR03().sumPUPt);


            nMu_++;
//          }
//        }
        

      }
    }

  }
}
