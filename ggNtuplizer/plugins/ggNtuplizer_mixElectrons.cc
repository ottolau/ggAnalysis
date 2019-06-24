#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
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
Int_t          nMix_;
vector<int>    maxEleCharge_pf_;
vector<float>  maxEleD0_pf_;
vector<float>  maxEleDz_pf_;
vector<float>  maxEleD0Error_pf_;
vector<float>  maxEleDzError_pf_;
vector<float>  maxElePt_pf_;
vector<float>  maxEleEta_pf_;
vector<float>  maxElePhi_pf_;
vector<float>  maxEleGsfPt_pf_;
vector<float>  maxEleGsfEta_pf_;
vector<float>  maxEleGsfPhi_pf_;
vector<float>  maxEleMVABWP_pf_;
vector<float>  maxEleMVAUnBWP_pf_;

vector<int>    maxEleCharge_lowPt_;
vector<float>  maxEleD0_lowPt_;
vector<float>  maxEleDz_lowPt_;
vector<float>  maxEleD0Error_lowPt_;
vector<float>  maxEleDzError_lowPt_;
vector<float>  maxElePt_lowPt_;
vector<float>  maxEleEta_lowPt_;
vector<float>  maxElePhi_lowPt_;
vector<float>  maxElePtMean_lowPt_;
vector<float>  maxEleEtaMean_lowPt_;
vector<float>  maxElePhiMean_lowPt_;
vector<float>  maxEleMVABWP_lowPt_;
vector<float>  maxEleMVAUnBWP_lowPt_;

vector<float> maxEleSvChi2_;
vector<float> maxEleSvNDOF_;
vector<float> maxEleSvProb_;
vector<float> maxEleSvX_;
vector<float> maxEleSvY_;
vector<float> maxEleSvZ_;
vector<float> maxEleSvXError_;
vector<float> maxEleSvYError_;
vector<float> maxEleSvZError_;
vector<float> maxEleSvMass_;
vector<float> maxEleSvCtxy_;
vector<float> maxEleSvCosAngle_;
vector<float> maxEleSvLxy_;
vector<float> maxEleSvLxyError_;

vector<int>    kaonMixCharge_lead_;
vector<float>  kaonMixD0_lead_;
vector<float>  kaonMixDz_lead_;
vector<float>  kaonMixD0Error_lead_;
vector<float>  kaonMixDzError_lead_;
vector<float>  kaonMixPt_lead_;
vector<float>  kaonMixEta_lead_;
vector<float>  kaonMixPhi_lead_;
vector<float>  kaonMixVx_lead_;
vector<float>  kaonMixVy_lead_;
vector<float>  kaonMixVz_lead_;
vector<float>  kaonMixTrkChi2_lead_;
vector<float>  kaonMixTrkNDOF_lead_;
vector<float>  kaonMixTrkNormChi2_lead_;
vector<float>  kaonMixJPsiMass_lead_;
vector<float>  kaonMixPhiMass_lead_;

vector<int>    kaonMixCharge_sublead_;
vector<float>  kaonMixD0_sublead_;
vector<float>  kaonMixDz_sublead_;
vector<float>  kaonMixD0Error_sublead_;
vector<float>  kaonMixDzError_sublead_;
vector<float>  kaonMixPt_sublead_;
vector<float>  kaonMixEta_sublead_;
vector<float>  kaonMixPhi_sublead_;
vector<float>  kaonMixVx_sublead_;
vector<float>  kaonMixVy_sublead_;
vector<float>  kaonMixVz_sublead_;
vector<float>  kaonMixTrkChi2_sublead_;
vector<float>  kaonMixTrkNDOF_sublead_;
vector<float>  kaonMixTrkNormChi2_sublead_;

vector<float>  bsMixdRele_;
vector<float>  bsMixdRkaon_;
vector<float>  bsMixdRJpsiPhi_;
vector<float>  bsMixJpsiMass_;
vector<float>  bsMixPhiMass_;
vector<float>  bsMixBsMass_;

void ggNtuplizer::branchesMixElectrons(TTree* tree) {

  tree->Branch("nMix",                    &nMix_);
  tree->Branch("maxEleCharge_pf",               &maxEleCharge_pf_);
  tree->Branch("maxEleD0_pf",                   &maxEleD0_pf_);
  tree->Branch("maxEleDz_pf",                   &maxEleDz_pf_);
  tree->Branch("maxEleD0Error_pf",              &maxEleD0Error_pf_);
  tree->Branch("maxEleDzError_pf",              &maxEleDzError_pf_);
  tree->Branch("maxElePt_pf",                   &maxElePt_pf_);
  tree->Branch("maxEleEta_pf",                  &maxEleEta_pf_);
  tree->Branch("maxElePhi_pf",                  &maxElePhi_pf_);
  tree->Branch("maxEleGsfPt_pf",                &maxEleGsfPt_pf_);
  tree->Branch("maxEleGsfEta_pf",               &maxEleGsfEta_pf_);
  tree->Branch("maxEleGsfPhi_pf",               &maxEleGsfPhi_pf_);

  tree->Branch("maxEleCharge_lowPt",               &maxEleCharge_lowPt_);
  tree->Branch("maxEleD0_lowPt",                   &maxEleD0_lowPt_);
  tree->Branch("maxEleDz_lowPt",                   &maxEleDz_lowPt_);
  tree->Branch("maxEleD0Error_lowPt",              &maxEleD0Error_lowPt_);
  tree->Branch("maxEleDzError_lowPt",              &maxEleDzError_lowPt_);
  tree->Branch("maxElePt_lowPt",                   &maxElePt_lowPt_);
  tree->Branch("maxEleEta_lowPt",                  &maxEleEta_lowPt_);
  tree->Branch("maxElePhi_lowPt",                  &maxElePhi_lowPt_);
  tree->Branch("maxElePtMean_lowPt",               &maxElePtMean_lowPt_);
  tree->Branch("maxEleEtaMean_lowPt",              &maxEleEtaMean_lowPt_);
  tree->Branch("maxElePhiMean_lowPt",              &maxElePhiMean_lowPt_);
  tree->Branch("maxEleMVABWP_lowPt",               &maxEleMVABWP_lowPt_);
  tree->Branch("maxEleMVAUnBWP_lowPt",             &maxEleMVAUnBWP_lowPt_);

  tree->Branch("maxEleSvChi2",                    &maxEleSvChi2_);
  tree->Branch("maxEleSvNDOF",                    &maxEleSvNDOF_);
  tree->Branch("maxEleSvProb",                    &maxEleSvProb_);
  tree->Branch("maxEleSvX",                       &maxEleSvX_);
  tree->Branch("maxEleSvY",                       &maxEleSvY_);
  tree->Branch("maxEleSvZ",                       &maxEleSvZ_);
  tree->Branch("maxEleSvXError",                  &maxEleSvXError_);
  tree->Branch("maxEleSvYError",                  &maxEleSvYError_);
  tree->Branch("maxEleSvZError",                  &maxEleSvZError_);
  tree->Branch("maxEleSvMass",                    &maxEleSvMass_);
  tree->Branch("maxEleSvCtxy",                    &maxEleSvCtxy_);
  tree->Branch("maxEleSvCosAngle",                    &maxEleSvCosAngle_);
  tree->Branch("maxEleSvLxy",                    	   &maxEleSvLxy_);
  tree->Branch("maxEleSvLxyError",                    &maxEleSvLxyError_);

  tree->Branch("kaonMixCharge_lead",               &kaonMixCharge_lead_);
  tree->Branch("kaonMixD0_lead",                   &kaonMixD0_lead_);
  tree->Branch("kaonMixDz_lead",                   &kaonMixDz_lead_);
  tree->Branch("kaonMixD0Error_lead",              &kaonMixD0Error_lead_);
  tree->Branch("kaonMixDzError_lead",              &kaonMixDzError_lead_);
  tree->Branch("kaonMixPt_lead",                   &kaonMixPt_lead_);
  tree->Branch("kaonMixEta_lead",                  &kaonMixEta_lead_);
  tree->Branch("kaonMixPhi_lead",                  &kaonMixPhi_lead_);
  tree->Branch("kaonMixVx_lead",                   &kaonMixVx_lead_);
  tree->Branch("kaonMixVy_lead",                   &kaonMixVy_lead_);
  tree->Branch("kaonMixVz_lead",                   &kaonMixVz_lead_);
  tree->Branch("kaonMixTrkChi2_lead",              &kaonMixTrkChi2_lead_);
  tree->Branch("kaonMixTrkNDOF_lead",              &kaonMixTrkNDOF_lead_);
  tree->Branch("kaonMixTrkNormChi2_lead",          &kaonMixTrkNormChi2_lead_);

  tree->Branch("kaonMixCharge_sublead",               &kaonMixCharge_sublead_);
  tree->Branch("kaonMixD0_sublead",                   &kaonMixD0_sublead_);
  tree->Branch("kaonMixDz_sublead",                   &kaonMixDz_sublead_);
  tree->Branch("kaonMixD0Error_sublead",              &kaonMixD0Error_sublead_);
  tree->Branch("kaonMixDzError_sublead",              &kaonMixDzError_sublead_);
  tree->Branch("kaonMixPt_sublead",                   &kaonMixPt_sublead_);
  tree->Branch("kaonMixEta_sublead",                  &kaonMixEta_sublead_);
  tree->Branch("kaonMixPhi_sublead",                  &kaonMixPhi_sublead_);
  tree->Branch("kaonMixVx_sublead",                   &kaonMixVx_sublead_);
  tree->Branch("kaonMixVy_sublead",                   &kaonMixVy_sublead_);
  tree->Branch("kaonMixVz_sublead",                   &kaonMixVz_sublead_);
  tree->Branch("kaonMixTrkChi2_sublead",              &kaonMixTrkChi2_sublead_);
  tree->Branch("kaonMixTrkNDOF_sublead",              &kaonMixTrkNDOF_sublead_);
  tree->Branch("kaonMixTrkNormChi2_sublead",          &kaonMixTrkNormChi2_sublead_);

  tree->Branch("bsMixdRele",                     &bsMixdRele_);
  tree->Branch("bsMixdRkaon",                    &bsMixdRkaon_);
  tree->Branch("bsMixdRJpsiPhi",                 &bsMixdRJpsiPhi_);
  tree->Branch("bsMixJpsiMass",                  &bsMixJpsiMass_);
  tree->Branch("bsMixPhiMass",                   &bsMixPhiMass_);
  tree->Branch("bsMixBsMass",                    &bsMixBsMass_);
  
}

void ggNtuplizer::fillMixElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  // cleanup from previous execution
  maxEleCharge_pf_                  .clear();
  maxEleD0_pf_                      .clear();
  maxEleDz_pf_                      .clear();
  maxEleD0Error_pf_                 .clear();
  maxEleDzError_pf_                 .clear();
  maxElePt_pf_                      .clear();
  maxEleEta_pf_                     .clear();
  maxElePhi_pf_                     .clear();
  maxEleGsfPt_pf_                   .clear();
  maxEleGsfEta_pf_                  .clear();
  maxEleGsfPhi_pf_                  .clear();

  maxEleCharge_lowPt_                  .clear();
  maxEleD0_lowPt_                      .clear();
  maxEleDz_lowPt_                      .clear();
  maxEleD0Error_lowPt_                 .clear();
  maxEleDzError_lowPt_                 .clear();
  maxElePt_lowPt_                      .clear();
  maxEleEta_lowPt_                     .clear();
  maxElePhi_lowPt_                     .clear();
  maxElePtMean_lowPt_                  .clear();
  maxEleEtaMean_lowPt_                 .clear();
  maxElePhiMean_lowPt_                 .clear();
  maxEleMVABWP_lowPt_                  .clear();
  maxEleMVAUnBWP_lowPt_                .clear();

  maxEleSvChi2_.clear();
  maxEleSvNDOF_.clear();
  maxEleSvProb_.clear();
  maxEleSvX_.clear();
  maxEleSvY_.clear();
  maxEleSvZ_.clear();
  maxEleSvXError_.clear();
  maxEleSvYError_.clear();
  maxEleSvZError_.clear();
  maxEleSvMass_.clear();
  maxEleSvCtxy_.clear();
  maxEleSvCosAngle_.clear();
  maxEleSvLxy_.clear();
  maxEleSvLxyError_.clear();

  kaonMixCharge_lead_                  .clear();
  kaonMixD0_lead_                      .clear();
  kaonMixDz_lead_                      .clear();
  kaonMixD0Error_lead_                 .clear();
  kaonMixDzError_lead_                 .clear();
  kaonMixPt_lead_                      .clear();
  kaonMixEta_lead_                     .clear();
  kaonMixPhi_lead_                     .clear();
  kaonMixVx_lead_                      .clear();
  kaonMixVy_lead_                      .clear();
  kaonMixVz_lead_                      .clear();
  kaonMixTrkChi2_lead_                 .clear();
  kaonMixTrkNDOF_lead_                 .clear();
  kaonMixTrkNormChi2_lead_             .clear();

  kaonMixCharge_sublead_                  .clear();
  kaonMixD0_sublead_                      .clear();
  kaonMixDz_sublead_                      .clear();
  kaonMixD0Error_sublead_                 .clear();
  kaonMixDzError_sublead_                 .clear();
  kaonMixPt_sublead_                      .clear();
  kaonMixEta_sublead_                     .clear();
  kaonMixPhi_sublead_                     .clear();
  kaonMixVx_sublead_                      .clear();
  kaonMixVy_sublead_                      .clear();
  kaonMixVz_sublead_                      .clear();
  kaonMixTrkChi2_sublead_                 .clear();
  kaonMixTrkNDOF_sublead_                 .clear();
  kaonMixTrkNormChi2_sublead_             .clear();

  bsMixdRele_		      .clear();
  bsMixdRkaon_                 .clear();
  bsMixdRJpsiPhi_              .clear();
  bsMixJpsiMass_		      .clear();
  bsMixPhiMass_		      .clear();
  bsMixBsMass_		      .clear();

  nMix_ = 0;

  if (isAOD_) {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<std::vector<reco::GsfElectron> > lowpTelectronHandle;
    e.getByToken(lowpTelectronlabel_, lowpTelectronHandle);

    edm::Handle<reco::TrackCollection> tracksHandle;
    e.getByToken(tracklabel_, tracksHandle);

    edm::Handle<edm::ValueMap<float> > ele_mva_wp_biased;
    e.getByToken( eleBWPToken_ ,ele_mva_wp_biased); 
    edm::Handle<edm::ValueMap<float> > ele_mva_wp_unbiased;
    e.getByToken( eleUnBWPToken_ ,ele_mva_wp_unbiased);

    edm::Handle<reco::ConversionCollection> conversions;
    e.getByToken(conversionsToken_, conversions);  

    if (!lowpTelectronHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no low pT electrons in event";
      return;
    }

    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      if (iEle->pt() < 0.5) continue;
      if (fabs(iEle->vz() - pv.z()) > 0.5) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      for (reco::GsfElectronCollection::const_iterator jEle = lowpTelectronHandle->begin(); jEle != lowpTelectronHandle->end(); ++jEle) {
	if (jEle->gsfTrack()->ptMode() < 0.5) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > 0.5) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->gsfTrack()->ptMode(), jEle->gsfTrack()->etaMode(), jEle->gsfTrack()->phiMode(), pmass);

	if ((iele_lv + jele_lv).M() < 2.6 || (iele_lv + jele_lv).M() > 3.6) continue;
	//if ((iele_lv + jele_lv).M() > 5.0) continue;

	//auto leadEle = iele_lv.Pt() > jele_lv.Pt() ? iEle : jEle;
	//auto subleadEle = iele_lv.Pt() > jele_lv.Pt() ? jEle : iEle;

	if ((*ele_mva_wp_unbiased)[jEle->gsfTrack()] < 0.19) continue;

	// remove duplicated electrons
	bool duplicatedEle = false;
	for (edm::View<pat::Electron>::const_iterator pfEle = electronHandle->begin(); pfEle != electronHandle->end(); ++pfEle) {
	  if (deltaR(jEle->eta(), jEle->phi(), pfEle->eta(), pfEle->phi()) < 0.001) {
	    duplicatedEle = true;
	    break;
	  }
	}

	if (duplicatedEle) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	//std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;

	//XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));


	//KinematicConstrainedVertexFitter kvFitter;
	//RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;

	for (reco::TrackCollection::const_iterator iHad = tracksHandle->begin(); iHad != tracksHandle->end(); ++iHad) {
	  if (iHad->pt() < 0.4) continue;
	  if (fabs(iHad->eta()) > 2.5) continue;
          if (fabs(iHad->vz() - pv.z()) > 0.5) continue;
	  if (iHad->normalizedChi2() < 0.0) continue;
	  if (iHad->normalizedChi2() > 20.0) continue;

	  for (reco::TrackCollection::const_iterator jHad = iHad+1; jHad != tracksHandle->end(); ++jHad) {
	    //if (iHad->charge()*jHad->charge() > 0.0) continue;
	    if (jHad->pt() <  0.4) continue;
	    if (fabs(jHad->eta()) > 2.5) continue;
	    if (fabs(jHad->vz() - pv.z()) > 0.5) continue;
	    if (jHad->normalizedChi2() < 0.0) continue;
	    if (jHad->normalizedChi2() > 20.0) continue;

	    // Phi mass window
	    float kpmass = 0.493677;
	    TLorentzVector iHad_lv, jHad_lv, bs_lv;
	    iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
	    jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);      
	    bs_lv = iele_lv + jele_lv + iHad_lv + jHad_lv;
	    if ((iHad_lv+jHad_lv).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.06) continue; 
	    if (bs_lv.M() < 4.5 || bs_lv.M() > 6.0) continue;

	    std::vector<RefCountedKinematicParticle> BsParticles;
	    float kpmasse = 1.e-6 * pmass;
	    //float bsM = 5.3663;

	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));


	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    //if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 20.0) continue;
	    if (TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()) < 0.001) continue;

	    // Accept these 4 tracks as a Bs candidate, fill ntuple

	    auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
	    auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

	    double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
	    
	    math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
	    math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	    math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	    double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	    if (cosAngle < 0.0) continue;

	    maxEleSvChi2_.push_back(DecayVtx->chiSquared());
	    maxEleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    maxEleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    maxEleSvX_.push_back(DecayVtx->position().x());
	    maxEleSvY_.push_back(DecayVtx->position().y());
	    maxEleSvZ_.push_back(DecayVtx->position().z());
	    maxEleSvXError_.push_back(DecayVtx->error().cxx());
	    maxEleSvYError_.push_back(DecayVtx->error().cyy());
	    maxEleSvZError_.push_back(DecayVtx->error().czz());
	    maxEleSvMass_.push_back(bs_lv.M());
	    maxEleSvCtxy_.push_back(ctxy);
	    maxEleSvCosAngle_.push_back(cosAngle);
	    maxEleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    maxEleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonMixCharge_lead_            .push_back(leadHad->charge());
	    kaonMixD0_lead_                .push_back(leadHad->dxy(pv));
	    kaonMixDz_lead_                .push_back(leadHad->dz(pv));
	    kaonMixD0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonMixDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonMixPt_lead_                .push_back(leadHad->pt());
	    kaonMixEta_lead_               .push_back(leadHad->eta());
	    kaonMixPhi_lead_               .push_back(leadHad->phi());
	    kaonMixVx_lead_ 		.push_back(leadHad->vx());
	    kaonMixVy_lead_ 		.push_back(leadHad->vy());
	    kaonMixVz_lead_ 		.push_back(leadHad->vz());
	    kaonMixTrkChi2_lead_ 		.push_back(leadHad->chi2());
	    kaonMixTrkNDOF_lead_ 		.push_back(leadHad->ndof());
	    kaonMixTrkNormChi2_lead_ 	.push_back(leadHad->normalizedChi2());

	    kaonMixCharge_sublead_         .push_back(subleadHad->charge());
	    kaonMixD0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonMixDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonMixD0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonMixDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonMixPt_sublead_             .push_back(subleadHad->pt());
	    kaonMixEta_sublead_            .push_back(subleadHad->eta());
	    kaonMixPhi_sublead_            .push_back(subleadHad->phi());
	    kaonMixVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonMixVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonMixVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonMixTrkChi2_sublead_ 	.push_back(subleadHad->chi2());
	    kaonMixTrkNDOF_sublead_ 	.push_back(subleadHad->ndof());
	    kaonMixTrkNormChi2_sublead_ 	.push_back(subleadHad->normalizedChi2());

	    bsMixdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsMixdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsMixdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsMixJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsMixPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsMixBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    maxEleCharge_pf_          .push_back(iEle->charge());
	    maxEleD0_pf_              .push_back(iEle->gsfTrack()->dxy(pv));
	    maxEleDz_pf_              .push_back(iEle->gsfTrack()->dz(pv));
	    maxEleD0Error_pf_         .push_back(iEle->gsfTrack()->dxyError());
	    maxEleDzError_pf_         .push_back(iEle->gsfTrack()->dzError());
	    maxElePt_pf_              .push_back(iEle->pt());
	    maxEleEta_pf_             .push_back(iEle->eta());
	    maxElePhi_pf_             .push_back(iEle->phi());
	    maxEleGsfPt_pf_           .push_back(iEle->gsfTrack()->ptMode());
	    maxEleGsfEta_pf_          .push_back(iEle->gsfTrack()->etaMode());
	    maxEleGsfPhi_pf_          .push_back(iEle->gsfTrack()->phiMode());

	    maxEleCharge_lowPt_          .push_back(jEle->gsfTrack()->charge());
	    maxEleD0_lowPt_              .push_back(jEle->gsfTrack()->dxy(pv));
	    maxEleDz_lowPt_              .push_back(jEle->gsfTrack()->dz(pv));
	    maxEleD0Error_lowPt_         .push_back(jEle->gsfTrack()->dxyError());
	    maxEleDzError_lowPt_         .push_back(jEle->gsfTrack()->dzError());
	    maxElePt_lowPt_              .push_back(jEle->gsfTrack()->ptMode());
	    maxEleEta_lowPt_             .push_back(jEle->gsfTrack()->etaMode());
	    maxElePhi_lowPt_             .push_back(jEle->gsfTrack()->phiMode());
	    maxElePtMean_lowPt_              .push_back(jEle->gsfTrack()->pt());
	    maxEleEtaMean_lowPt_             .push_back(jEle->gsfTrack()->eta());
	    maxElePhiMean_lowPt_             .push_back(jEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_lowPt = jEle->gsfTrack();
	    maxEleMVABWP_lowPt_          .push_back((*ele_mva_wp_biased)[mvaSeed_lowPt]);
	    maxEleMVAUnBWP_lowPt_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_lowPt]);

	    nMix_++;
	  }
	}
      }
    }
  } else {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<std::vector<reco::GsfElectron> > lowpTelectronHandle;
    e.getByToken(lowpTelectronlabel_, lowpTelectronHandle);

    edm::Handle<edm::ValueMap<float> > ele_mva_wp_biased;
    e.getByToken( eleBWPToken_ ,ele_mva_wp_biased); 
    edm::Handle<edm::ValueMap<float> > ele_mva_wp_unbiased;
    e.getByToken( eleUnBWPToken_ ,ele_mva_wp_unbiased);

    edm::Handle<reco::ConversionCollection> conversions;
    e.getByToken(conversionsToken_, conversions);  

    if (!lowpTelectronHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no low pT electrons in event";
      return;
    }

    edm::Handle<pat::PackedCandidateCollection> pfcands;
    e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<pat::PackedCandidateCollection> losttracks;
    e.getByToken(lostTracksLabel_, losttracks);

    std::vector<pat::PackedCandidate> alltracks;
    alltracks.reserve(pfcands->size() + losttracks->size());
    alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
    alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      if (iEle->pt() < 0.5) continue;
      if (fabs(iEle->vz() - pv.z()) > 0.5) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      for (reco::GsfElectronCollection::const_iterator jEle = lowpTelectronHandle->begin(); jEle != lowpTelectronHandle->end(); ++jEle) {
	//if (iEle->charge()*jEle->charge() > 0.0) continue;
	if (jEle->gsfTrack()->ptMode() < 0.5) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > 0.5) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	//if ((iele_lv + jele_lv).M() < 2.4 || (iele_lv + jele_lv).M() > 3.8) continue;
	if ((iele_lv + jele_lv).M() > 5.0) continue;

	//auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
	//auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;

	if ((*ele_mva_wp_unbiased)[jEle->gsfTrack()] < 0.19) continue;

	// remove duplicated electrons
	bool duplicatedEle = false;
	for (edm::View<pat::Electron>::const_iterator pfEle = electronHandle->begin(); pfEle != electronHandle->end(); ++pfEle) {
	  if (deltaR(jEle->eta(), jEle->phi(), pfEle->eta(), pfEle->phi()) < 0.001) {
	    duplicatedEle = true;
	    break;
	  }
	}

	if (duplicatedEle) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	//std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;

	//XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	//KinematicConstrainedVertexFitter kvFitter;
	//RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;

	for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {
	  if (iHad->pt() <= 0.4) continue;
          if (iHad->charge() == 0) continue;
          if (abs(iHad->pdgId()) != 211) continue;
          if (iHad->bestTrack() == nullptr) continue;
	  if (fabs(iHad->eta()) > 2.5) continue;
	  if (fabs(iHad->vz() - pv.z()) > 1.0) continue;
	  //if (iHad->normalizedChi2() < 0.0) continue;
	  //if (iHad->normalizedChi2() > 20.0) continue;

	  for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != alltracks.end(); ++jHad) {
	    if (jHad->pt() <= 0.4) continue;
            if (jHad->charge() == 0) continue;
            if (abs(jHad->pdgId()) != 211) continue;
            if (jHad->bestTrack() == nullptr) continue;
	    //if (iHad->charge()*jHad->charge() > 0.0) continue;
	    if (fabs(jHad->vz() - pv.z()) > 1.0) continue;
	    //if (jHad->normalizedChi2() < 0.0) continue;
	    //if (jHad->normalizedChi2() > 20) continue;

	    // Phi mass window
	    float kpmass = 0.493677;
	    TLorentzVector iHad_lv, jHad_lv, bs_lv;
	    iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
	    jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);      
	    bs_lv = iele_lv + jele_lv + iHad_lv + jHad_lv;
	    if ((iHad_lv+jHad_lv).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.10) continue; 
	    if ((iele_lv + jele_lv + iHad_lv + jHad_lv).M() < 4.0 || (iele_lv + jele_lv + iHad_lv + jHad_lv).M() > 6.0) continue;
	    if (fabs(jHad->eta()) > 2.5) continue;

	    std::vector<RefCountedKinematicParticle> BsParticles;
	    float kpmasse = 1.e-6 * pmass;
	    //float bsM = 5.3663;

	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    //if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 20.0) continue;
	    if (TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()) < 0.001) continue;

	    // Accept these 4 tracks as a Bs candidate, fill ntuple

	    auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
	    auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

	    double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
	    
	    math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
	    math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	    math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	    double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	    if (cosAngle < 0.0) continue;

	    maxEleSvChi2_.push_back(DecayVtx->chiSquared());
	    maxEleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    maxEleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    maxEleSvX_.push_back(DecayVtx->position().x());
	    maxEleSvY_.push_back(DecayVtx->position().y());
	    maxEleSvZ_.push_back(DecayVtx->position().z());
	    maxEleSvXError_.push_back(DecayVtx->error().cxx());
	    maxEleSvYError_.push_back(DecayVtx->error().cyy());
	    maxEleSvZError_.push_back(DecayVtx->error().czz());
	    maxEleSvMass_.push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());
	    maxEleSvCtxy_.push_back(ctxy);
	    maxEleSvCosAngle_.push_back(cosAngle);
	    maxEleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    maxEleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonMixCharge_lead_            .push_back(leadHad->charge());
	    kaonMixD0_lead_                .push_back(leadHad->dxy(pv));
	    kaonMixDz_lead_                .push_back(leadHad->dz(pv));
	    kaonMixD0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonMixDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonMixPt_lead_                .push_back(leadHad->pt());
	    kaonMixEta_lead_               .push_back(leadHad->eta());
	    kaonMixPhi_lead_               .push_back(leadHad->phi());
	    kaonMixVx_lead_ 		.push_back(leadHad->vx());
	    kaonMixVy_lead_ 		.push_back(leadHad->vy());
	    kaonMixVz_lead_ 		.push_back(leadHad->vz());
	    kaonMixTrkChi2_lead_ 		.push_back(leadHad->bestTrack()->chi2());
	    kaonMixTrkNDOF_lead_ 		.push_back(leadHad->bestTrack()->ndof());
	    kaonMixTrkNormChi2_lead_ 	.push_back(leadHad->bestTrack()->normalizedChi2());

	    kaonMixCharge_sublead_         .push_back(subleadHad->charge());
	    kaonMixD0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonMixDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonMixD0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonMixDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonMixPt_sublead_             .push_back(subleadHad->pt());
	    kaonMixEta_sublead_            .push_back(subleadHad->eta());
	    kaonMixPhi_sublead_            .push_back(subleadHad->phi());
	    kaonMixVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonMixVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonMixVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonMixTrkChi2_sublead_ 	.push_back(subleadHad->bestTrack()->chi2());
	    kaonMixTrkNDOF_sublead_ 	.push_back(subleadHad->bestTrack()->ndof());
	    kaonMixTrkNormChi2_sublead_ 	.push_back(subleadHad->bestTrack()->normalizedChi2());

	    bsMixdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsMixdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsMixdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsMixJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsMixPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsMixBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    maxEleCharge_pf_          .push_back(iEle->charge());
	    maxEleD0_pf_              .push_back(iEle->gsfTrack()->dxy(pv));
	    maxEleDz_pf_              .push_back(iEle->gsfTrack()->dz(pv));
	    maxEleD0Error_pf_         .push_back(iEle->gsfTrack()->dxyError());
	    maxEleDzError_pf_         .push_back(iEle->gsfTrack()->dzError());
	    maxElePt_pf_              .push_back(iEle->pt());
	    maxEleEta_pf_             .push_back(iEle->eta());
	    maxElePhi_pf_             .push_back(iEle->phi());
	    maxEleGsfPt_pf_           .push_back(iEle->gsfTrack()->ptMode());
	    maxEleGsfEta_pf_          .push_back(iEle->gsfTrack()->etaMode());
	    maxEleGsfPhi_pf_          .push_back(iEle->gsfTrack()->phiMode());

	    maxEleCharge_lowPt_          .push_back(jEle->gsfTrack()->charge());
	    maxEleD0_lowPt_              .push_back(jEle->gsfTrack()->dxy(pv));
	    maxEleDz_lowPt_              .push_back(jEle->gsfTrack()->dz(pv));
	    maxEleD0Error_lowPt_         .push_back(jEle->gsfTrack()->dxyError());
	    maxEleDzError_lowPt_         .push_back(jEle->gsfTrack()->dzError());
	    maxElePt_lowPt_              .push_back(jEle->gsfTrack()->ptMode());
	    maxEleEta_lowPt_             .push_back(jEle->gsfTrack()->etaMode());
	    maxElePhi_lowPt_             .push_back(jEle->gsfTrack()->phiMode());
	    maxElePtMean_lowPt_              .push_back(jEle->gsfTrack()->pt());
	    maxEleEtaMean_lowPt_             .push_back(jEle->gsfTrack()->eta());
	    maxElePhiMean_lowPt_             .push_back(jEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_lowPt = jEle->gsfTrack();
	    maxEleMVABWP_lowPt_          .push_back((*ele_mva_wp_biased)[mvaSeed_lowPt]);
	    maxEleMVAUnBWP_lowPt_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_lowPt]);


	    nMix_++;
	  }
	}
      }
    }

  }

}


