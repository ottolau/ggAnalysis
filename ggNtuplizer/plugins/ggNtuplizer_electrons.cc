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
vector<pair<int,int>>    eleCharge_;
vector<pair<int,int>>    eleChargeConsistent_;
vector<pair<float,float>>  eleEn_;
vector<pair<float,float>>  eleSCEn_;
vector<pair<float,float>>  eleEcalEn_;
vector<pair<float,float>>  eleESEnP1_;
vector<pair<float,float>>  eleESEnP2_;
//vector<pair<float,float>>  eleESEnP1Raw_;
//vector<pair<float,float>>  eleESEnP2Raw_;
vector<pair<float,float>>  eleD0_;
vector<pair<float,float>>  eleDz_;
vector<pair<float,float>>  eleSIP_;
vector<pair<float,float>>  elePt_;
vector<pair<float,float>>  eleEta_;
vector<pair<float,float>>  elePhi_;
vector<pair<float,float>>  eleR9_;
vector<pair<float,float>>  eleCalibPt_;
vector<pair<float,float>>  eleCalibEn_;
vector<pair<float,float>>  eleSCEta_;
vector<pair<float,float>>  eleSCPhi_;
vector<pair<float,float>>  eleSCRawEn_;
vector<pair<float,float>>  eleSCEtaWidth_;
vector<pair<float,float>>  eleSCPhiWidth_;
vector<pair<float,float>>  eleHoverE_;
vector<pair<float,float>>  eleEoverP_;
vector<pair<float,float>>  eleEoverPout_;
vector<pair<float,float>>  eleEoverPInv_;
vector<pair<float,float>>  eleBrem_;
vector<pair<float,float>>  eledEtaAtVtx_;
vector<pair<float,float>>  eledPhiAtVtx_;
vector<pair<float,float>>  eledEtaAtCalo_;
//vector<pair<float,float>>  eleSigmaIEtaIEta_;
//vector<pair<float,float>>  eleSigmaIEtaIPhi_;
//vector<pair<float,float>>  eleSigmaIPhiIPhi_;
vector<pair<float,float>>  eleSigmaIEtaIEtaFull5x5_;
vector<pair<float,float>>  eleSigmaIPhiIPhiFull5x5_;
vector<pair<int,int>>    eleConvVeto_;
vector<pair<int,int>>    eleMissHits_;
vector<pair<float,float>>  eleESEffSigmaRR_;
vector<pair<float,float>>  elePFChIso_;
vector<pair<float,float>>  elePFPhoIso_;
vector<pair<float,float>>  elePFNeuIso_;
vector<pair<float,float>>  elePFPUIso_;
vector<pair<float,float>>  elePFClusEcalIso_;
vector<pair<float,float>>  elePFClusHcalIso_;
//vector<pair<float,float>>  elePFMiniIso_;
vector<pair<float,float>>  eleIDMVAIso_;
vector<pair<float,float>>  eleIDMVANoIso_;
vector<pair<float,float>>  eledEtaseedAtVtx_;
vector<pair<float,float>>  eleE1x5_;
vector<pair<float,float>>  eleE2x5_;
vector<pair<float,float>>  eleE5x5_;
vector<pair<float,float>>  eleE1x5Full5x5_;
vector<pair<float,float>>  eleE2x5Full5x5_;
vector<pair<float,float>>  eleE5x5Full5x5_;
vector<pair<float,float>>  eleR9Full5x5_;
vector<pair<int,int>>    eleEcalDrivenSeed_;
vector<pair<float,float>>  eleDr03EcalRecHitSumEt_;
vector<pair<float,float>>  eleDr03HcalDepth1TowerSumEt_;
vector<pair<float,float>>  eleDr03HcalDepth2TowerSumEt_;
vector<pair<float,float>>  eleDr03HcalTowerSumEt_;
vector<pair<float,float>>  eleDr03TkSumPt_;
vector<pair<float,float>>  elecaloEnergy_;
vector<pair<float,float>>  eleTrkdxy_;
vector<pair<float,float>>  eleKFHits_;
vector<pair<float,float>>  eleKFChi2_;
vector<pair<float,float>>  eleGSFChi2_;
//vector<pair<unsigned long,unsigned long>> eleFiredSingleTrgs_;
//vector<pair<unsigned long,unsigned long>> eleFiredDoubleTrgs_;
//vector<pair<unsigned long,unsigned long>> eleFiredL1Trgs_;
//vector<pair<unsigned short,unsigned short>> eleIDbit_;
vector<UShort_t> eleIDbitFirst_;
vector<UShort_t> eleIDbitSecond_;

vector<float> svChi2_;
vector<float> svNDOF_;
vector<float> svX_;
vector<float> svY_;
vector<float> svZ_;
vector<float> svXError_;
vector<float> svYError_;
vector<float> svZError_;
vector<float> svMass_;
vector<float> svDxySig_;
vector<float> svCosAngle_;

vector<pair<int,int>>    eeHadCharge_;
vector<pair<float,float>>  eeHadD0_;
vector<pair<float,float>>  eeHadDz_;
vector<pair<float,float>>  eeHadD0Error_;
vector<pair<float,float>>  eeHadDzError_;
vector<pair<float,float>>  eeHadPt_;
vector<pair<float,float>>  eeHadEta_;
vector<pair<float,float>>  eeHadPhi_;
vector<pair<float,float>>  eeHadVx_;
vector<pair<float,float>>  eeHadVy_;
vector<pair<float,float>>  eeHadVz_;
//vector<pair<float,float>>  eeHadEn_;
vector<pair<float,float>>  eeHadTrkChi2_;
vector<pair<float,float>>  eeHadTrkNDOF_;
vector<pair<float,float>>  eeHadTrkNormChi2_;
vector<float>  eeHadJPsiMass_;
vector<float>  eeHadPhiMass_;
vector<pair<float,float>> eeTrkdEdx_;


void ggNtuplizer::branchesElectrons(TTree* tree) {

  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge",               &eleCharge_);
  tree->Branch("eleChargeConsistent",     &eleChargeConsistent_);
  tree->Branch("eleEn",                   &eleEn_);
  tree->Branch("eleSCEn",                 &eleSCEn_);
  tree->Branch("eleEcalEn",               &eleEcalEn_);
  tree->Branch("eleESEnP1",               &eleESEnP1_);
  tree->Branch("eleESEnP2",               &eleESEnP2_);
  tree->Branch("eleD0",                   &eleD0_);
  tree->Branch("eleDz",                   &eleDz_);
  tree->Branch("eleSIP",                  &eleSIP_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_);
  tree->Branch("eleR9",                   &eleR9_);
  tree->Branch("eleCalibPt",              &eleCalibPt_);
  tree->Branch("eleCalibEn",              &eleCalibEn_);
  tree->Branch("eleSCEta",                &eleSCEta_);
  tree->Branch("eleSCPhi",                &eleSCPhi_);
  tree->Branch("eleSCRawEn",              &eleSCRawEn_);
  tree->Branch("eleSCEtaWidth",           &eleSCEtaWidth_);
  tree->Branch("eleSCPhiWidth",           &eleSCPhiWidth_);
  tree->Branch("eleHoverE",               &eleHoverE_);
  tree->Branch("eleEoverP",               &eleEoverP_);
  tree->Branch("eleEoverPout",            &eleEoverPout_);
  tree->Branch("eleEoverPInv",            &eleEoverPInv_);
  tree->Branch("eleBrem",                 &eleBrem_);
  tree->Branch("eledEtaAtVtx",            &eledEtaAtVtx_);
  tree->Branch("eledPhiAtVtx",            &eledPhiAtVtx_);
  tree->Branch("eledEtaAtCalo",           &eledEtaAtCalo_);
  //tree->Branch("eleSigmaIEtaIEta",        &eleSigmaIEtaIEta_);
  //tree->Branch("eleSigmaIEtaIPhi",        &eleSigmaIEtaIPhi_);
  //tree->Branch("eleSigmaIPhiIPhi",        &eleSigmaIPhiIPhi_);
  tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree->Branch("eleConvVeto",             &eleConvVeto_);
  tree->Branch("eleMissHits",             &eleMissHits_);
  tree->Branch("eleESEffSigmaRR",         &eleESEffSigmaRR_);
  tree->Branch("elePFChIso",              &elePFChIso_);
  tree->Branch("elePFPhoIso",             &elePFPhoIso_);
  tree->Branch("elePFNeuIso",             &elePFNeuIso_);
  tree->Branch("elePFPUIso",              &elePFPUIso_);
  tree->Branch("elePFClusEcalIso",        &elePFClusEcalIso_);
  tree->Branch("elePFClusHcalIso",        &elePFClusHcalIso_);
//  tree->Branch("elePFMiniIso",            &elePFMiniIso_);
  tree->Branch("eleIDMVAIso",             &eleIDMVAIso_);
  tree->Branch("eleIDMVANoIso",           &eleIDMVANoIso_);
  tree->Branch("eledEtaseedAtVtx",        &eledEtaseedAtVtx_);
  tree->Branch("eleE1x5",                 &eleE1x5_);
  tree->Branch("eleE2x5",                 &eleE2x5_);
  tree->Branch("eleE5x5",                 &eleE5x5_);
  tree->Branch("eleE1x5Full5x5",          &eleE1x5Full5x5_);
  tree->Branch("eleE2x5Full5x5",          &eleE2x5Full5x5_);
  tree->Branch("eleE5x5Full5x5",          &eleE5x5Full5x5_);
  tree->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  tree->Branch("eleEcalDrivenSeed",           &eleEcalDrivenSeed_);
  tree->Branch("eleDr03EcalRecHitSumEt",      &eleDr03EcalRecHitSumEt_);
  tree->Branch("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt_);
  tree->Branch("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt_);
  tree->Branch("eleDr03HcalTowerSumEt",       &eleDr03HcalTowerSumEt_);
  tree->Branch("eleDr03TkSumPt",              &eleDr03TkSumPt_);
  tree->Branch("elecaloEnergy",               &elecaloEnergy_);
  tree->Branch("eleTrkdxy",                   &eleTrkdxy_);
  tree->Branch("eleKFHits",                   &eleKFHits_);
  tree->Branch("eleKFChi2",                   &eleKFChi2_);
  tree->Branch("eleGSFChi2",                  &eleGSFChi2_);
//  tree->Branch("eleFiredSingleTrgs",          &eleFiredSingleTrgs_);
//  tree->Branch("eleFiredDoubleTrgs",          &eleFiredDoubleTrgs_);
//  tree->Branch("eleFiredL1Trgs",              &eleFiredL1Trgs_);
  tree->Branch("eleIDbitFirst",                    &eleIDbitFirst_);
  tree->Branch("eleIDbitSecond",                    &eleIDbitSecond_);
  tree->Branch("svChi2",                    &svChi2_);
  tree->Branch("svNDOF",                    &svNDOF_);
  tree->Branch("svX",                       &svX_);
  tree->Branch("svY",                       &svY_);
  tree->Branch("svZ",                       &svZ_);
  tree->Branch("svXError",                  &svXError_);
  tree->Branch("svYError",                  &svYError_);
  tree->Branch("svZError",                  &svZError_);
  tree->Branch("svMass",                    &svMass_);
  tree->Branch("svDxySig",                    &svDxySig_);
  tree->Branch("svCosAngle",                    &svCosAngle_);

  if (!separateVtxFit_) {
    tree->Branch("eeHadCharge",               &eeHadCharge_);
    tree->Branch("eeHadD0",                   &eeHadD0_);
    tree->Branch("eeHadDz",                   &eeHadDz_);
    tree->Branch("eeHadD0Error",              &eeHadD0Error_);
    tree->Branch("eeHadDzError",              &eeHadDzError_);
    tree->Branch("eeHadPt",                   &eeHadPt_);
    tree->Branch("eeHadEta",                  &eeHadEta_);
    tree->Branch("eeHadPhi",                  &eeHadPhi_);
    tree->Branch("eeHadVx",                   &eeHadVx_);
    tree->Branch("eeHadVy",                   &eeHadVy_);
    tree->Branch("eeHadVz",                   &eeHadVz_);
//    tree->Branch("eeHadEn",                   &eeHadEn_);
    tree->Branch("eeHadTrkChi2",              &eeHadTrkChi2_);
    tree->Branch("eeHadTrkNDOF",              &eeHadTrkNDOF_);
    tree->Branch("eeHadTrkNormChi2",          &eeHadTrkNormChi2_);
    tree->Branch("eeHadJPsiMass",                  &eeHadJPsiMass_);
    tree->Branch("eeHadPhiMass",                  &eeHadPhiMass_);
    tree->Branch("eeTrkdEdx",                  &eeTrkdEdx_);
 

  }

  
}

void ggNtuplizer::fillElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv) {
    
  // cleanup from previous execution
  eleCharge_                  .clear();
  eleChargeConsistent_        .clear();
  eleEn_                      .clear();
  eleSCEn_                    .clear();
  eleEcalEn_                  .clear();
  eleESEnP1_                  .clear();
  eleESEnP2_                  .clear();
  //eleESEnP1Raw_               .clear();
  //eleESEnP2Raw_               .clear();
  eleD0_                      .clear();
  eleDz_                      .clear();
  eleSIP_                     .clear();
  elePt_                      .clear();
  eleEta_                     .clear();
  elePhi_                     .clear();
  eleR9_                      .clear();
  eleCalibPt_                 .clear();
  eleCalibEn_                 .clear();
  eleSCEta_                   .clear();
  eleSCPhi_                   .clear();
  eleSCRawEn_                 .clear();
  eleSCEtaWidth_              .clear();
  eleSCPhiWidth_              .clear();
  eleHoverE_                  .clear();
  eleEoverP_                  .clear();
  eleEoverPout_               .clear();
  eleEoverPInv_               .clear();
  eleBrem_                    .clear();
  eledEtaAtVtx_               .clear();
  eledPhiAtVtx_               .clear();
  eledEtaAtCalo_              .clear();
  //eleSigmaIEtaIEta_           .clear();
  //eleSigmaIEtaIPhi_           .clear();
  //eleSigmaIPhiIPhi_           .clear();
  eleSigmaIEtaIEtaFull5x5_    .clear();
  eleSigmaIPhiIPhiFull5x5_    .clear();
  eleConvVeto_                .clear();
  eleMissHits_                .clear();
  eleESEffSigmaRR_            .clear();
  elePFChIso_                 .clear();
  elePFPhoIso_                .clear();
  elePFNeuIso_                .clear();
  elePFPUIso_                 .clear();
  elePFClusEcalIso_           .clear();
  elePFClusHcalIso_           .clear();
//  elePFMiniIso_               .clear();
  eleIDMVAIso_                .clear();
  eleIDMVANoIso_              .clear();
  eledEtaseedAtVtx_           .clear();
  eleE1x5_                    .clear();
  eleE2x5_                    .clear();
  eleE5x5_                    .clear();
  eleEcalDrivenSeed_          .clear();
  eleDr03EcalRecHitSumEt_     .clear();
  eleDr03HcalDepth1TowerSumEt_.clear();
  eleDr03HcalDepth2TowerSumEt_.clear();
  eleDr03HcalTowerSumEt_      .clear();
  eleDr03TkSumPt_             .clear();
  eleE1x5Full5x5_             .clear();
  eleE2x5Full5x5_             .clear();
  eleE5x5Full5x5_             .clear();
  eleR9Full5x5_               .clear();
  elecaloEnergy_              .clear();
  eleTrkdxy_                  .clear();
  eleKFHits_                  .clear();
  eleKFChi2_                  .clear();
  eleGSFChi2_                 .clear();
//  eleFiredSingleTrgs_         .clear();
//  eleFiredDoubleTrgs_         .clear();
//  eleFiredL1Trgs_             .clear();
//  eleIDbit_                   .clear();
  eleIDbitFirst_                   .clear();
  eleIDbitSecond_                  .clear();
  svChi2_.clear();
  svNDOF_.clear();
  svX_.clear();
  svY_.clear();
  svZ_.clear();
  svXError_.clear();
  svYError_.clear();
  svZError_.clear();
  svMass_.clear();
  svDxySig_.clear();
  svCosAngle_.clear();

  eeHadCharge_                  .clear();
  eeHadD0_                      .clear();
  eeHadDz_                      .clear();
  eeHadD0Error_                 .clear();
  eeHadDzError_                 .clear();
  eeHadPt_                      .clear();
  eeHadEta_                     .clear();
  eeHadPhi_                     .clear();
  eeHadVx_                      .clear();
  eeHadVy_                      .clear();
  eeHadVz_                      .clear();
//  eeHadEn_	              .clear();
  eeHadTrkChi2_                 .clear();
  eeHadTrkNDOF_                 .clear();
  eeHadTrkNormChi2_             .clear();
  eeHadJPsiMass_		      .clear();
  eeHadPhiMass_		      .clear();
  eeTrkdEdx_			.clear();

  nEle_ = 0;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
  e.getByToken(calibelectronCollection_, calibelectronHandle);

//  edm::Handle<pat::PackedCandidateCollection> pfcands;
//  e.getByToken(pckPFCandidateCollection_, pfcands);

  edm::Handle<reco::TrackCollection> tracksHandle;
  e.getByToken(tracklabel_, tracksHandle);

  edm::Handle<reco::DeDxDataValueMap> dEdxObjectHandle;
  e.getByToken(deDxProducer_, dEdxObjectHandle );
  const edm::ValueMap<reco::DeDxData> dEdxColl = *dEdxObjectHandle.product();


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

  if (!separateVtxFit_) {

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      for (edm::View<pat::Electron>::const_iterator jEle = iEle+1; jEle != electronHandle->end(); ++jEle) {
	if (iEle->charge()*jEle->charge() > 0.0) continue;
	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	if ((iele_lv + jele_lv).M() < 2.4 || (iele_lv + jele_lv).M() > 3.8) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;
	reco::Track ieletrk = eletrks[(iEle-electronHandle->begin())].second;
	reco::Track jeletrk = eletrks[(jEle-electronHandle->begin())].second;
	const reco::TransientTrack ielettk = getTransientTrack( ieletrk );
	const reco::TransientTrack jelettk = getTransientTrack( jeletrk );

	//XParticles.push_back(pFactory.particle(getTransientTrack( ieletrk ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( jeletrk ), pmass, 0.0, 0, pmasse));

	XParticles.push_back(pFactory.particle(ielettk, pmass, 0.0, 0, pmasse));
	XParticles.push_back(pFactory.particle(jelettk, pmass, 0.0, 0, pmasse));

	KinematicConstrainedVertexFitter kvFitter;
	RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 ||KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
	//KinVtx->movePointerToTheTop();
	//RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();

	//for (pat::PackedCandidateCollection::const_iterator iHad = pfcands->begin(); iHad != pfcands->end(); ++iHad) {
        for (reco::TrackCollection::const_iterator iHad = tracksHandle->begin(); iHad != tracksHandle->end(); ++iHad) {
	  if (fabs(iHad->eta()) > 2.5) continue;
	  if (fabs(ieletrk.vz() - iHad->vz()) > 1) continue;
	  if (iHad->normalizedChi2() < 0.0) continue;
	  if (iHad->normalizedChi2() > 20.0) continue;

          //for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != pfcands->end(); ++jHad) {
          for (reco::TrackCollection::const_iterator jHad = iHad+1; jHad != tracksHandle->end(); ++jHad) {
	    if (iHad->charge()*jHad->charge() > 0.0) continue;
	    if (fabs(ieletrk.vz() - jHad->vz()) > 1) continue;
	    if (jHad->normalizedChi2() < 0.0) continue;
	    if (jHad->normalizedChi2() > 20) continue;

            // Phi mass window
	    float kpmass = 0.493677;
	    TLorentzVector iHad_lv, jHad_lv;
	    iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
	    jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);      
	    //if (((iHad_lv+jHad_lv)).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.06) continue; 
	    if (((iHad_lv+jHad_lv)).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.10) continue; 
	    if (fabs(jHad->eta()) > 2.5) continue;

	    std::vector<RefCountedKinematicParticle> BsParticles;
	    float kpmasse = 1.e-6 * pmass;

	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(ielettk, pmass, 0.0, 0, pmasse));
	    BsParticles.push_back(pFactory.particle(jelettk, pmass, 0.0, 0, pmasse));

	    //BsParticles.push_back(jpsi_part);

	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 30.0) continue;
	    TVector3 vtxDisplace(DecayVtx->position().x()-pv.x(), DecayVtx->position().y()-pv.y(), DecayVtx->position().z()-pv.z());
	    float cosAngle = vtxDisplace.Dot((iele_lv + jele_lv + iHad_lv + jHad_lv).Vect())/(vtxDisplace.Mag()*(iele_lv + jele_lv + iHad_lv + jHad_lv).Vect().Mag());
	    //if (cosAngle < 0.0) continue;

	    svChi2_.push_back(DecayVtx->chiSquared());
	    svNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    svX_.push_back(DecayVtx->position().x());
	    svY_.push_back(DecayVtx->position().y());
	    svZ_.push_back(DecayVtx->position().z());
	    svXError_.push_back(DecayVtx->error().cxx());
	    svYError_.push_back(DecayVtx->error().cyy());
	    svZError_.push_back(DecayVtx->error().czz());
	    svMass_.push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());
	    float dxy = TMath::Sqrt((DecayVtx->position().x()-pv.x())*(DecayVtx->position().x()-pv.x()) + (DecayVtx->position().y()-pv.y())*(DecayVtx->position().y()-pv.y()));
	    float sigmadxy = TMath::Sqrt(DecayVtx->error().cxx()*DecayVtx->error().cxx() + DecayVtx->error().cyy()*DecayVtx->error().cyy());
	    svDxySig_.push_back(dxy/sigmadxy);
	    svCosAngle_.push_back(cosAngle);

	    eeHadCharge_          .push_back(make_pair(iHad->charge(),jHad->charge()));
	    eeHadD0_              .push_back(make_pair(iHad->dxy(pv),jHad->dxy(pv)));
	    eeHadDz_              .push_back(make_pair(iHad->dz(pv),jHad->dz(pv)));
	    eeHadD0Error_		.push_back(make_pair(iHad->dxyError(),jHad->dxyError()));
	    eeHadDzError_		.push_back(make_pair(iHad->dzError(),jHad->dzError()));

	    eeHadPt_              .push_back(make_pair(iHad->pt(),jHad->pt()));
	    eeHadEta_             .push_back(make_pair(iHad->eta(),jHad->eta()));
	    eeHadPhi_             .push_back(make_pair(iHad->phi(),jHad->phi()));
	    eeHadVx_		.push_back(make_pair(iHad->vx(),jHad->vx()));
	    eeHadVy_		.push_back(make_pair(iHad->vy(),jHad->vy()));
	    eeHadVz_		.push_back(make_pair(iHad->vz(),jHad->vz()));

	    eeHadTrkChi2_		.push_back(make_pair(iHad->chi2(),jHad->chi2()));
	    eeHadTrkNDOF_		.push_back(make_pair(iHad->ndof(),jHad->ndof()));
	    eeHadTrkNormChi2_	.push_back(make_pair(iHad->normalizedChi2(),jHad->normalizedChi2()));
	    eeHadJPsiMass_        .push_back((iele_lv+jele_lv).M());
	    eeHadPhiMass_         .push_back((iHad_lv+jHad_lv).M());
	    eeTrkdEdx_		.push_back(make_pair(dEdxColl[reco::TrackRef(tracksHandle, iHad-tracksHandle->begin())].dEdx(),dEdxColl[reco::TrackRef(tracksHandle, jHad-tracksHandle->begin())].dEdx()));

	    Float_t corrPtfirst = -1;
	    Float_t corrEnfirst = -1;
	    Float_t corrPtsecond = -1;
	    Float_t corrEnsecond = -1;

	    for (edm::View<pat::Electron>::const_iterator iCEle = calibelectronHandle->begin(); iCEle != calibelectronHandle->end(); ++iCEle) {

	      if (fabs(iEle->eta() - iCEle->eta()) < 0.001 && fabs(iEle->phi() - iCEle->phi()) < 0.001) {
		corrPtfirst = iCEle->pt();
		corrEnfirst = iCEle->energy();
	      }
	      if (fabs(jEle->eta() - iCEle->eta()) < 0.001 && fabs(jEle->phi() - iCEle->phi()) < 0.001) {
		corrPtsecond = iCEle->pt();
		corrEnsecond = iCEle->energy();
	      }

	    }
	    eleCalibPt_        .push_back(make_pair(corrPtfirst,corrPtsecond));
	    eleCalibEn_        .push_back(make_pair(corrEnfirst,corrEnsecond));

	    eleCharge_          .push_back(make_pair(iEle->charge(),jEle->charge()));
	    eleChargeConsistent_.push_back(make_pair((Int_t)iEle->isGsfCtfScPixChargeConsistent(),(Int_t)jEle->isGsfCtfScPixChargeConsistent()));
	    eleEn_              .push_back(make_pair(iEle->energy(),jEle->energy()));
	    eleD0_              .push_back(make_pair(iEle->gsfTrack()->dxy(pv),jEle->gsfTrack()->dxy(pv)));
	    eleDz_              .push_back(make_pair(iEle->gsfTrack()->dz(pv),jEle->gsfTrack()->dz(pv)));
	    eleSIP_             .push_back(make_pair(fabs(iEle->dB(pat::Electron::PV3D))/iEle->edB(pat::Electron::PV3D),fabs(jEle->dB(pat::Electron::PV3D))/jEle->edB(pat::Electron::PV3D)));
	    elePt_              .push_back(make_pair(iEle->pt(),jEle->pt()));
	    eleEta_             .push_back(make_pair(iEle->eta(),jEle->eta()));
	    elePhi_             .push_back(make_pair(iEle->phi(),jEle->phi()));
	    eleR9_              .push_back(make_pair(iEle->r9(),jEle->r9()));
	    eleSCEn_            .push_back(make_pair(iEle->superCluster()->energy(),jEle->superCluster()->energy()));
	    eleEcalEn_          .push_back(make_pair(iEle->ecalEnergy(),jEle->ecalEnergy()));
	    eleESEnP1_          .push_back(make_pair(iEle->superCluster()->preshowerEnergyPlane1(),jEle->superCluster()->preshowerEnergyPlane1()));
	    eleESEnP2_          .push_back(make_pair(iEle->superCluster()->preshowerEnergyPlane2(),jEle->superCluster()->preshowerEnergyPlane2()));
	    eleSCEta_           .push_back(make_pair(iEle->superCluster()->eta(),jEle->superCluster()->eta()));
	    eleSCPhi_           .push_back(make_pair(iEle->superCluster()->phi(),jEle->superCluster()->phi()));
	    eleSCRawEn_         .push_back(make_pair(iEle->superCluster()->rawEnergy(),jEle->superCluster()->rawEnergy()));
	    eleSCEtaWidth_      .push_back(make_pair(iEle->superCluster()->etaWidth(),jEle->superCluster()->etaWidth()));
	    eleSCPhiWidth_      .push_back(make_pair(iEle->superCluster()->phiWidth(),jEle->superCluster()->phiWidth()));
	    eleHoverE_          .push_back(make_pair(iEle->hcalOverEcal(),jEle->hcalOverEcal()));

	    //eleFiredSingleTrgs_ .push_back(make_pair(matchSingleElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()),matchSingleElectronTriggerFilters(jEle->pt(), jEle->eta(), jEle->phi())));

	    //eleFiredDoubleTrgs_ .push_back(make_pair(matchDoubleElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()),matchDoubleElectronTriggerFilters(jEle->pt(), jEle->eta(), jEle->phi())));
	    //eleFiredL1Trgs_     .push_back(make_pair(matchL1TriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()),matchL1TriggerFilters(jEle->pt(), jEle->eta(), jEle->phi())));

	    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
	    eleEoverP_          .push_back(make_pair(iEle->eSuperClusterOverP(),jEle->eSuperClusterOverP()));
	    eleEoverPout_       .push_back(make_pair(iEle->eEleClusterOverPout(),jEle->eEleClusterOverPout()));
	    eleBrem_            .push_back(make_pair(iEle->fbrem(),jEle->fbrem()));
	    eledEtaAtVtx_       .push_back(make_pair(iEle->deltaEtaSuperClusterTrackAtVtx(),jEle->deltaEtaSuperClusterTrackAtVtx()));
	    eledPhiAtVtx_       .push_back(make_pair(iEle->deltaPhiSuperClusterTrackAtVtx(),jEle->deltaPhiSuperClusterTrackAtVtx()));
	    eledEtaAtCalo_      .push_back(make_pair(iEle->deltaEtaSeedClusterTrackAtCalo(),jEle->deltaEtaSeedClusterTrackAtCalo()));
	    //eleSigmaIEtaIEta_   .push_back(iEle->sigmaIetaIeta()); ///new sigmaietaieta
	    //eleSigmaIEtaIPhi_   .push_back(iEle->sigmaIetaIphi());
	    //eleSigmaIPhiIPhi_   .push_back(iEle->sigmaIphiIphi());
	    eleConvVeto_        .push_back(make_pair((Int_t)iEle->passConversionVeto(),(Int_t)jEle->passConversionVeto())); // ConvVtxFit || missHit == 0
	    eleMissHits_        .push_back(make_pair(iEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS),jEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)));
	    eleESEffSigmaRR_    .push_back(make_pair(lazyTool.eseffsirir(*((*iEle).superCluster())),lazyTool.eseffsirir(*((*jEle).superCluster()))));

	    // VID calculation of (1/E - 1/p)
	    double ieleEoverPInv = 1e30;
	    double jeleEoverPInv = 1e30;
	    if (iEle->ecalEnergy() != 0 && std::isfinite(iEle->ecalEnergy()))	ieleEoverPInv = (1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy();
	    if (jEle->ecalEnergy() != 0 && std::isfinite(jEle->ecalEnergy()))	jeleEoverPInv = (1.0 - jEle->eSuperClusterOverP())/jEle->ecalEnergy();
	    eleEoverPInv_.push_back(make_pair(ieleEoverPInv,jeleEoverPInv));

	    //if (iEle->ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
	    //else if (!std::isfinite(iEle->ecalEnergy()))  eleEoverPInv_.push_back(1e30);
	    //else  eleEoverPInv_.push_back((1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy());

	    ///HEEP ID
	    double ieledEtaseedAtVtx = iEle->superCluster().isNonnull() && iEle->superCluster()->seed().isNonnull() ?
	      iEle->deltaEtaSuperClusterTrackAtVtx() - iEle->superCluster()->eta() + iEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
	    double jeledEtaseedAtVtx = jEle->superCluster().isNonnull() && jEle->superCluster()->seed().isNonnull() ?
	      jEle->deltaEtaSuperClusterTrackAtVtx() - jEle->superCluster()->eta() + jEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

	    eledEtaseedAtVtx_   .push_back(make_pair(ieledEtaseedAtVtx,jeledEtaseedAtVtx));

	    eleE1x5_            .push_back(make_pair(iEle->e1x5(),jEle->e1x5()));
	    eleE2x5_            .push_back(make_pair(iEle->e2x5Max(),jEle->e2x5Max()));
	    eleE5x5_            .push_back(make_pair(iEle->e5x5(),jEle->e5x5()));

	    reco::GsfElectron::PflowIsolationVariables ipfIso = iEle->pfIsolationVariables();
	    reco::GsfElectron::PflowIsolationVariables jpfIso = jEle->pfIsolationVariables();

	    elePFChIso_         .push_back(make_pair(ipfIso.sumChargedHadronPt,jpfIso.sumChargedHadronPt));
	    elePFPhoIso_        .push_back(make_pair(ipfIso.sumPhotonEt,jpfIso.sumPhotonEt));
	    elePFNeuIso_        .push_back(make_pair(ipfIso.sumNeutralHadronEt,jpfIso.sumNeutralHadronEt));
	    elePFPUIso_         .push_back(make_pair(ipfIso.sumPUPt,jpfIso.sumPUPt));
	    elecaloEnergy_      .push_back(make_pair(iEle->caloEnergy(),jEle->caloEnergy()));
	    //elePFMiniIso_       .push_back(make_pair(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*iEle)), 0.05, 0.2, 10., false),getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*jEle)), 0.05, 0.2, 10., false)));

	    /////quantities which were used for Run1 - these do not
	    ///calculated through PF (meaning no energy is subtracted
	    ///using PF)
	    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
	    ///line 120

	    eleSigmaIEtaIEtaFull5x5_    .push_back(make_pair(iEle->full5x5_sigmaIetaIeta(),jEle->full5x5_sigmaIetaIeta()));
	    eleSigmaIPhiIPhiFull5x5_    .push_back(make_pair(iEle->full5x5_sigmaIphiIphi(),jEle->full5x5_sigmaIphiIphi()));
	    eleE1x5Full5x5_             .push_back(make_pair(iEle->full5x5_e1x5(),jEle->full5x5_e1x5()));
	    eleE2x5Full5x5_             .push_back(make_pair(iEle->full5x5_e2x5Max(),jEle->full5x5_e2x5Max()));
	    eleE5x5Full5x5_             .push_back(make_pair(iEle->full5x5_e5x5(),jEle->full5x5_e5x5()));
	    eleR9Full5x5_               .push_back(make_pair(iEle->full5x5_r9(),jEle->full5x5_r9()));

	    ///For HEEP ID
	    eleEcalDrivenSeed_          .push_back(make_pair(iEle->ecalDrivenSeed(),jEle->ecalDrivenSeed()));
	    eleDr03EcalRecHitSumEt_     .push_back(make_pair(iEle->dr03EcalRecHitSumEt(),jEle->dr03EcalRecHitSumEt()));
	    eleDr03HcalDepth1TowerSumEt_.push_back(make_pair(iEle->dr03HcalDepth1TowerSumEt(),jEle->dr03HcalDepth1TowerSumEt()));
	    eleDr03HcalDepth2TowerSumEt_.push_back(make_pair(iEle->dr03HcalDepth2TowerSumEt(),jEle->dr03HcalDepth2TowerSumEt()));
	    eleDr03HcalTowerSumEt_      .push_back(make_pair(iEle->dr03HcalTowerSumEt(),jEle->dr03HcalTowerSumEt()));
	    eleDr03TkSumPt_             .push_back(make_pair(iEle->dr03TkSumPt(),jEle->dr03TkSumPt()));

	    reco::GsfTrackRef igsfTrackRef = iEle->gsfTrack();
	    reco::GsfTrackRef jgsfTrackRef = jEle->gsfTrack();

	    if (iEle->gsfTrack().isNonnull() && jEle->gsfTrack().isNonnull()) {
	      eleGSFChi2_.push_back(make_pair(igsfTrackRef->normalizedChi2(),jgsfTrackRef->normalizedChi2()));
	      if (recVtxs->size() > 0)
		eleTrkdxy_.push_back(make_pair(igsfTrackRef->dxy(recVtxs->front().position()),jgsfTrackRef->dxy(recVtxs->front().position())));
	      else
		eleTrkdxy_.push_back(make_pair(-999,-999));
	    } else {
	      eleGSFChi2_.push_back(make_pair(999.,999.));
	      eleTrkdxy_.push_back(make_pair(-999,-999));
	    }
	    
	    reco::TrackRef ikfTrackRef = iEle->closestCtfTrackRef();
	    reco::TrackRef jkfTrackRef = jEle->closestCtfTrackRef();

	    if (ikfTrackRef.isAvailable() && ikfTrackRef.isNonnull() && jkfTrackRef.isAvailable() && jkfTrackRef.isNonnull()) {
	      eleKFHits_.push_back(make_pair(ikfTrackRef->hitPattern().trackerLayersWithMeasurement(),jkfTrackRef->hitPattern().trackerLayersWithMeasurement()));
	      eleKFChi2_.push_back(make_pair(ikfTrackRef->normalizedChi2(),jkfTrackRef->normalizedChi2()));
	    } else {
	      eleKFHits_.push_back(make_pair(-1.,-1.));
	      eleKFChi2_.push_back(make_pair(999.,999.));
	    }

	    //edm::Ptr<reco::GsfElectron> recoEl(iEle);      
	    //const auto el = electrons->ptrAt(nEle_);
	    const auto iel = electronHandle->ptrAt(iEle - electronHandle->begin());
	    const auto jel = electronHandle->ptrAt(jEle - electronHandle->begin());
	   
	    unsigned short itmpeleIDbit = 0;
	    unsigned short jtmpeleIDbit = 0;
	   
	    ///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works

	    bool isPassVeto  = (*veto_id_decisions)[iel->originalObjectRef()];
	    if (isPassVeto) setbit(itmpeleIDbit, 0);
	
	    bool isPassLoose  = (*loose_id_decisions)[iel->originalObjectRef()];
	    if (isPassLoose) setbit(itmpeleIDbit, 1);

	    bool isPassMedium = (*medium_id_decisions)[iel->originalObjectRef()];
	    if (isPassMedium) setbit(itmpeleIDbit, 2);

	    bool isPassTight  = (*tight_id_decisions)[iel->originalObjectRef()];
	    if (isPassTight) setbit(itmpeleIDbit, 3);

	    bool isPassHEEP = (*heep_id_decisions)[iel->originalObjectRef()];
	    if (isPassHEEP) setbit(itmpeleIDbit, 4);

	    isPassVeto  = (*veto_id_decisions)[jel->originalObjectRef()];
	    if (isPassVeto) setbit(jtmpeleIDbit, 0);
	    
	    isPassLoose  = (*loose_id_decisions)[jel->originalObjectRef()];
	    if (isPassLoose) setbit(jtmpeleIDbit, 1);
	    
	    isPassMedium = (*medium_id_decisions)[jel->originalObjectRef()];
	    if (isPassMedium) setbit(jtmpeleIDbit, 2);
	    
	    isPassTight  = (*tight_id_decisions)[jel->originalObjectRef()];
	    if (isPassTight) setbit(jtmpeleIDbit, 3);
	    
	    isPassHEEP = (*heep_id_decisions)[jel->originalObjectRef()];
	    if (isPassHEEP) setbit(jtmpeleIDbit, 4);

	    eleIDMVAIso_  .push_back(make_pair((*eleMVAIsoValues)[iel->originalObjectRef()],(*eleMVAIsoValues)[jel->originalObjectRef()]));
	    eleIDMVANoIso_.push_back(make_pair((*eleMVANoIsoValues)[iel->originalObjectRef()],(*eleMVANoIsoValues)[jel->originalObjectRef()]));

	    elePFClusEcalIso_.push_back(make_pair(iEle->ecalPFClusterIso(),jEle->ecalPFClusterIso()));
	    elePFClusHcalIso_.push_back(make_pair(iEle->hcalPFClusterIso(),jEle->hcalPFClusterIso()));

	    eleIDbitFirst_.push_back(itmpeleIDbit);
	    eleIDbitSecond_.push_back(jtmpeleIDbit);

	    //eleIDbit_.push_back(make_pair(itmpeleIDbit,jtmpeleIDbit));

	    nEle_++;

	  }
	}
      }
    }
  }





  if (separateVtxFit_) {

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      for (edm::View<pat::Electron>::const_iterator jEle = iEle+1; jEle != electronHandle->end(); ++jEle) {

	KinematicParticleFactoryFromTransientTrack pFactory;  
	std::vector<RefCountedKinematicParticle> XParticles;
	float pmass  = 0.0005109989461;
	float pmasse = 1.e-6 * pmass;
	reco::Track ieletrk = eletrks[(iEle-electronHandle->begin())].second;
	reco::Track jeletrk = eletrks[(jEle-electronHandle->begin())].second;

	XParticles.push_back(pFactory.particle(getTransientTrack( ieletrk ), pmass, 0.0, 0, pmasse));
	XParticles.push_back(pFactory.particle(getTransientTrack( jeletrk ), pmass, 0.0, 0, pmasse));

	KinematicConstrainedVertexFitter kvFitter;
	RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	if (KinVtx->isValid()) {
	  RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();

	  svChi2_.push_back(DecayVtx->chiSquared());
	  svNDOF_.push_back(DecayVtx->degreesOfFreedom());
	  svX_.push_back(DecayVtx->position().x());
	  svY_.push_back(DecayVtx->position().y());
	  svZ_.push_back(DecayVtx->position().z());
	  svXError_.push_back(DecayVtx->error().cxx());
	  svYError_.push_back(DecayVtx->error().cyy());
	  svZError_.push_back(DecayVtx->error().czz());
	  float elecM = 5.11e-4;
	  TLorentzVector iele_lv, jele_lv;
	  iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), elecM);
	  jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), elecM);
	  svMass_.push_back((iele_lv+jele_lv).M());
	  float dxy = TMath::Sqrt((DecayVtx->position().x()-pv.x())*(DecayVtx->position().x()-pv.x()) + (DecayVtx->position().y()-pv.y())*(DecayVtx->position().y()-pv.y()));
	  float sigmadxy = TMath::Sqrt(DecayVtx->error().cxx()*DecayVtx->error().cxx() + DecayVtx->error().cyy()*DecayVtx->error().cyy());
	  svDxySig_.push_back(dxy/sigmadxy);
	  TVector3 vtxDisplace(DecayVtx->position().x()-pv.x(), DecayVtx->position().y()-pv.y(), DecayVtx->position().z()-pv.z());
	  float cosAngle = vtxDisplace.Dot((iele_lv+jele_lv).Vect())/(vtxDisplace.Mag()*(iele_lv+jele_lv).Vect().Mag());
	  svCosAngle_.push_back(cosAngle);

	  Float_t corrPtfirst = -1;
	  Float_t corrEnfirst = -1;
	  Float_t corrPtsecond = -1;
	  Float_t corrEnsecond = -1;

	  for (edm::View<pat::Electron>::const_iterator iCEle = calibelectronHandle->begin(); iCEle != calibelectronHandle->end(); ++iCEle) {

	    if (fabs(iEle->eta() - iCEle->eta()) < 0.001 && fabs(iEle->phi() - iCEle->phi()) < 0.001) {
	      corrPtfirst = iCEle->pt();
	      corrEnfirst = iCEle->energy();
	    }
	    if (fabs(jEle->eta() - iCEle->eta()) < 0.001 && fabs(jEle->phi() - iCEle->phi()) < 0.001) {
	      corrPtsecond = iCEle->pt();
	      corrEnsecond = iCEle->energy();
	    }

	  }
	  eleCalibPt_        .push_back(make_pair(corrPtfirst,corrPtsecond));
	  eleCalibEn_        .push_back(make_pair(corrEnfirst,corrEnsecond));

	  eleCharge_          .push_back(make_pair(iEle->charge(),jEle->charge()));
	  eleChargeConsistent_.push_back(make_pair((Int_t)iEle->isGsfCtfScPixChargeConsistent(),(Int_t)jEle->isGsfCtfScPixChargeConsistent()));
	  eleEn_              .push_back(make_pair(iEle->energy(),jEle->energy()));
	  eleD0_              .push_back(make_pair(iEle->gsfTrack()->dxy(pv),jEle->gsfTrack()->dxy(pv)));
	  eleDz_              .push_back(make_pair(iEle->gsfTrack()->dz(pv),jEle->gsfTrack()->dz(pv)));
	  eleSIP_             .push_back(make_pair(fabs(iEle->dB(pat::Electron::PV3D))/iEle->edB(pat::Electron::PV3D),fabs(jEle->dB(pat::Electron::PV3D))/jEle->edB(pat::Electron::PV3D)));
	  elePt_              .push_back(make_pair(iEle->pt(),jEle->pt()));
	  eleEta_             .push_back(make_pair(iEle->eta(),jEle->eta()));
	  elePhi_             .push_back(make_pair(iEle->phi(),jEle->phi()));
	  eleR9_              .push_back(make_pair(iEle->r9(),jEle->r9()));
	  eleSCEn_            .push_back(make_pair(iEle->superCluster()->energy(),jEle->superCluster()->energy()));
	  eleEcalEn_          .push_back(make_pair(iEle->ecalEnergy(),jEle->ecalEnergy()));
	  eleESEnP1_          .push_back(make_pair(iEle->superCluster()->preshowerEnergyPlane1(),jEle->superCluster()->preshowerEnergyPlane1()));
	  eleESEnP2_          .push_back(make_pair(iEle->superCluster()->preshowerEnergyPlane2(),jEle->superCluster()->preshowerEnergyPlane2()));
	  eleSCEta_           .push_back(make_pair(iEle->superCluster()->eta(),jEle->superCluster()->eta()));
	  eleSCPhi_           .push_back(make_pair(iEle->superCluster()->phi(),jEle->superCluster()->phi()));
	  eleSCRawEn_         .push_back(make_pair(iEle->superCluster()->rawEnergy(),jEle->superCluster()->rawEnergy()));
	  eleSCEtaWidth_      .push_back(make_pair(iEle->superCluster()->etaWidth(),jEle->superCluster()->etaWidth()));
	  eleSCPhiWidth_      .push_back(make_pair(iEle->superCluster()->phiWidth(),jEle->superCluster()->phiWidth()));
	  eleHoverE_          .push_back(make_pair(iEle->hcalOverEcal(),jEle->hcalOverEcal()));

	  //eleFiredSingleTrgs_ .push_back(make_pair(matchSingleElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()),matchSingleElectronTriggerFilters(jEle->pt(), jEle->eta(), jEle->phi())));

	  //eleFiredDoubleTrgs_ .push_back(make_pair(matchDoubleElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()),matchDoubleElectronTriggerFilters(jEle->pt(), jEle->eta(), jEle->phi())));
	  //eleFiredL1Trgs_     .push_back(make_pair(matchL1TriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()),matchL1TriggerFilters(jEle->pt(), jEle->eta(), jEle->phi())));

	  ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
	  eleEoverP_          .push_back(make_pair(iEle->eSuperClusterOverP(),jEle->eSuperClusterOverP()));
	  eleEoverPout_       .push_back(make_pair(iEle->eEleClusterOverPout(),jEle->eEleClusterOverPout()));
	  eleBrem_            .push_back(make_pair(iEle->fbrem(),jEle->fbrem()));
	  eledEtaAtVtx_       .push_back(make_pair(iEle->deltaEtaSuperClusterTrackAtVtx(),jEle->deltaEtaSuperClusterTrackAtVtx()));
	  eledPhiAtVtx_       .push_back(make_pair(iEle->deltaPhiSuperClusterTrackAtVtx(),jEle->deltaPhiSuperClusterTrackAtVtx()));
	  eledEtaAtCalo_      .push_back(make_pair(iEle->deltaEtaSeedClusterTrackAtCalo(),jEle->deltaEtaSeedClusterTrackAtCalo()));
	  //eleSigmaIEtaIEta_   .push_back(iEle->sigmaIetaIeta()); ///new sigmaietaieta
	  //eleSigmaIEtaIPhi_   .push_back(iEle->sigmaIetaIphi());
	  //eleSigmaIPhiIPhi_   .push_back(iEle->sigmaIphiIphi());
	  eleConvVeto_        .push_back(make_pair((Int_t)iEle->passConversionVeto(),(Int_t)jEle->passConversionVeto())); // ConvVtxFit || missHit == 0
	  eleMissHits_        .push_back(make_pair(iEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS),jEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)));
	  eleESEffSigmaRR_    .push_back(make_pair(lazyTool.eseffsirir(*((*iEle).superCluster())),lazyTool.eseffsirir(*((*jEle).superCluster()))));

	  // VID calculation of (1/E - 1/p)
	  double ieleEoverPInv = 1e30;
	  double jeleEoverPInv = 1e30;
	  if (iEle->ecalEnergy() != 0 && std::isfinite(iEle->ecalEnergy()))	ieleEoverPInv = (1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy();
	  if (jEle->ecalEnergy() != 0 && std::isfinite(jEle->ecalEnergy()))	jeleEoverPInv = (1.0 - jEle->eSuperClusterOverP())/jEle->ecalEnergy();
	  eleEoverPInv_.push_back(make_pair(ieleEoverPInv,jeleEoverPInv));

	  //if (iEle->ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
	  //else if (!std::isfinite(iEle->ecalEnergy()))  eleEoverPInv_.push_back(1e30);
	  //else  eleEoverPInv_.push_back((1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy());

	  ///HEEP ID
	  double ieledEtaseedAtVtx = iEle->superCluster().isNonnull() && iEle->superCluster()->seed().isNonnull() ?
	    iEle->deltaEtaSuperClusterTrackAtVtx() - iEle->superCluster()->eta() + iEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
	  double jeledEtaseedAtVtx = jEle->superCluster().isNonnull() && jEle->superCluster()->seed().isNonnull() ?
	    jEle->deltaEtaSuperClusterTrackAtVtx() - jEle->superCluster()->eta() + jEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

	  eledEtaseedAtVtx_   .push_back(make_pair(ieledEtaseedAtVtx,jeledEtaseedAtVtx));

	  eleE1x5_            .push_back(make_pair(iEle->e1x5(),jEle->e1x5()));
	  eleE2x5_            .push_back(make_pair(iEle->e2x5Max(),jEle->e2x5Max()));
	  eleE5x5_            .push_back(make_pair(iEle->e5x5(),jEle->e5x5()));

	  reco::GsfElectron::PflowIsolationVariables ipfIso = iEle->pfIsolationVariables();
	  reco::GsfElectron::PflowIsolationVariables jpfIso = jEle->pfIsolationVariables();

	  elePFChIso_         .push_back(make_pair(ipfIso.sumChargedHadronPt,jpfIso.sumChargedHadronPt));
	  elePFPhoIso_        .push_back(make_pair(ipfIso.sumPhotonEt,jpfIso.sumPhotonEt));
	  elePFNeuIso_        .push_back(make_pair(ipfIso.sumNeutralHadronEt,jpfIso.sumNeutralHadronEt));
	  elePFPUIso_         .push_back(make_pair(ipfIso.sumPUPt,jpfIso.sumPUPt));
	  elecaloEnergy_      .push_back(make_pair(iEle->caloEnergy(),jEle->caloEnergy()));
	  //elePFMiniIso_       .push_back(make_pair(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*iEle)), 0.05, 0.2, 10., false),getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*jEle)), 0.05, 0.2, 10., false)));

	  /////quantities which were used for Run1 - these do not
	  ///calculated through PF (meaning no energy is subtracted
	  ///using PF)
	  ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
	  ///line 120

	  eleSigmaIEtaIEtaFull5x5_    .push_back(make_pair(iEle->full5x5_sigmaIetaIeta(),jEle->full5x5_sigmaIetaIeta()));
	  eleSigmaIPhiIPhiFull5x5_    .push_back(make_pair(iEle->full5x5_sigmaIphiIphi(),jEle->full5x5_sigmaIphiIphi()));
	  eleE1x5Full5x5_             .push_back(make_pair(iEle->full5x5_e1x5(),jEle->full5x5_e1x5()));
	  eleE2x5Full5x5_             .push_back(make_pair(iEle->full5x5_e2x5Max(),jEle->full5x5_e2x5Max()));
	  eleE5x5Full5x5_             .push_back(make_pair(iEle->full5x5_e5x5(),jEle->full5x5_e5x5()));
	  eleR9Full5x5_               .push_back(make_pair(iEle->full5x5_r9(),jEle->full5x5_r9()));

	  ///For HEEP ID
	  eleEcalDrivenSeed_          .push_back(make_pair(iEle->ecalDrivenSeed(),jEle->ecalDrivenSeed()));
	  eleDr03EcalRecHitSumEt_     .push_back(make_pair(iEle->dr03EcalRecHitSumEt(),jEle->dr03EcalRecHitSumEt()));
	  eleDr03HcalDepth1TowerSumEt_.push_back(make_pair(iEle->dr03HcalDepth1TowerSumEt(),jEle->dr03HcalDepth1TowerSumEt()));
	  eleDr03HcalDepth2TowerSumEt_.push_back(make_pair(iEle->dr03HcalDepth2TowerSumEt(),jEle->dr03HcalDepth2TowerSumEt()));
	  eleDr03HcalTowerSumEt_      .push_back(make_pair(iEle->dr03HcalTowerSumEt(),jEle->dr03HcalTowerSumEt()));
	  eleDr03TkSumPt_             .push_back(make_pair(iEle->dr03TkSumPt(),jEle->dr03TkSumPt()));

	  reco::GsfTrackRef igsfTrackRef = iEle->gsfTrack();
	  reco::GsfTrackRef jgsfTrackRef = jEle->gsfTrack();

	  if (iEle->gsfTrack().isNonnull() && jEle->gsfTrack().isNonnull()) {
	    eleGSFChi2_.push_back(make_pair(igsfTrackRef->normalizedChi2(),jgsfTrackRef->normalizedChi2()));
	    if (recVtxs->size() > 0)
	      eleTrkdxy_.push_back(make_pair(igsfTrackRef->dxy(recVtxs->front().position()),jgsfTrackRef->dxy(recVtxs->front().position())));
	    else
	      eleTrkdxy_.push_back(make_pair(-999,-999));
	  } else {
	    eleGSFChi2_.push_back(make_pair(999.,999.));
	    eleTrkdxy_.push_back(make_pair(-999,-999));
	  }
	  
	  reco::TrackRef ikfTrackRef = iEle->closestCtfTrackRef();
	  reco::TrackRef jkfTrackRef = jEle->closestCtfTrackRef();

	  if (ikfTrackRef.isAvailable() && ikfTrackRef.isNonnull() && jkfTrackRef.isAvailable() && jkfTrackRef.isNonnull()) {
	    eleKFHits_.push_back(make_pair(ikfTrackRef->hitPattern().trackerLayersWithMeasurement(),jkfTrackRef->hitPattern().trackerLayersWithMeasurement()));
	    eleKFChi2_.push_back(make_pair(ikfTrackRef->normalizedChi2(),jkfTrackRef->normalizedChi2()));
	  } else {
	    eleKFHits_.push_back(make_pair(-1.,-1.));
	    eleKFChi2_.push_back(make_pair(999.,999.));
	  }

	  //edm::Ptr<reco::GsfElectron> recoEl(iEle);      
	  //const auto el = electrons->ptrAt(nEle_);
	  const auto iel = electronHandle->ptrAt(iEle - electronHandle->begin());
	  const auto jel = electronHandle->ptrAt(jEle - electronHandle->begin());
	 
	  unsigned short itmpeleIDbit = 0;
	  unsigned short jtmpeleIDbit = 0;
	 
	  ///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works

	  bool isPassVeto  = (*veto_id_decisions)[iel];
	  if (isPassVeto) setbit(itmpeleIDbit, 0);
	  
	  bool isPassLoose  = (*loose_id_decisions)[iel];
	  if (isPassLoose) setbit(itmpeleIDbit, 1);
	  
	  bool isPassMedium = (*medium_id_decisions)[iel];
	  if (isPassMedium) setbit(itmpeleIDbit, 2);
	  
	  bool isPassTight  = (*tight_id_decisions)[iel];
	  if (isPassTight) setbit(itmpeleIDbit, 3);
	  
	  bool isPassHEEP = (*heep_id_decisions)[iel];
	  if (isPassHEEP) setbit(itmpeleIDbit, 4);

	  isPassVeto  = (*veto_id_decisions)[jel];
	  if (isPassVeto) setbit(jtmpeleIDbit, 0);
	  
	  isPassLoose  = (*loose_id_decisions)[jel];
	  if (isPassLoose) setbit(jtmpeleIDbit, 1);
	  
	  isPassMedium = (*medium_id_decisions)[jel];
	  if (isPassMedium) setbit(jtmpeleIDbit, 2);
	  
	  isPassTight  = (*tight_id_decisions)[jel];
	  if (isPassTight) setbit(jtmpeleIDbit, 3);
	  
	  isPassHEEP = (*heep_id_decisions)[jel];
	  if (isPassHEEP) setbit(jtmpeleIDbit, 4);

	  eleIDMVAIso_  .push_back(make_pair((*eleMVAIsoValues)[iel],(*eleMVAIsoValues)[jel]));
	  eleIDMVANoIso_.push_back(make_pair((*eleMVANoIsoValues)[iel],(*eleMVANoIsoValues)[jel]));

	  elePFClusEcalIso_.push_back(make_pair(iEle->ecalPFClusterIso(),jEle->ecalPFClusterIso()));
	  elePFClusHcalIso_.push_back(make_pair(iEle->hcalPFClusterIso(),jEle->hcalPFClusterIso()));

	  eleIDbitFirst_.push_back(itmpeleIDbit);
	  eleIDbitSecond_.push_back(jtmpeleIDbit);

	  //eleIDbit_.push_back(make_pair(itmpeleIDbit,jtmpeleIDbit));

	  nEle_++;
	}
      }
    }
  }

}


