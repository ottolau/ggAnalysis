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
Int_t          nLowPt_;
vector<int>    lowPtCharge_lead_;
vector<float>  lowPtD0_lead_;
vector<float>  lowPtDz_lead_;
vector<float>  lowPtD0Error_lead_;
vector<float>  lowPtDzError_lead_;
vector<float>  lowPtPt_lead_;
vector<float>  lowPtEta_lead_;
vector<float>  lowPtPhi_lead_;
vector<float>  lowPtPtMean_lead_;
vector<float>  lowPtEtaMean_lead_;
vector<float>  lowPtPhiMean_lead_;
vector<float>  lowPtMVABWP_lead_;
vector<float>  lowPtMVAUnBWP_lead_;

vector<int>    lowPtCharge_sublead_;
vector<float>  lowPtD0_sublead_;
vector<float>  lowPtDz_sublead_;
vector<float>  lowPtD0Error_sublead_;
vector<float>  lowPtDzError_sublead_;
vector<float>  lowPtPt_sublead_;
vector<float>  lowPtEta_sublead_;
vector<float>  lowPtPhi_sublead_;
vector<float>  lowPtPtMean_sublead_;
vector<float>  lowPtEtaMean_sublead_;
vector<float>  lowPtPhiMean_sublead_;
vector<float>  lowPtMVABWP_sublead_;
vector<float>  lowPtMVAUnBWP_sublead_;

vector<float> lowPtSvChi2_;
vector<float> lowPtSvNDOF_;
vector<float> lowPtSvProb_;
vector<float> lowPtSvX_;
vector<float> lowPtSvY_;
vector<float> lowPtSvZ_;
vector<float> lowPtSvXError_;
vector<float> lowPtSvYError_;
vector<float> lowPtSvZError_;
vector<float> lowPtSvMass_;
vector<float> lowPtSvCtxy_;
vector<float> lowPtSvCosAngle_;
vector<float> lowPtSvLxy_;
vector<float> lowPtSvLxyError_;


void ggNtuplizer::branchesLowPtElectrons(TTree* tree) {

  tree->Branch("nLowPt",                    &nLowPt_);
  tree->Branch("lowPtCharge_lead",               &lowPtCharge_lead_);
  tree->Branch("lowPtD0_lead",                   &lowPtD0_lead_);
  tree->Branch("lowPtDz_lead",                   &lowPtDz_lead_);
  tree->Branch("lowPtD0Error_lead",              &lowPtD0Error_lead_);
  tree->Branch("lowPtDzError_lead",              &lowPtDzError_lead_);
  tree->Branch("lowPtPt_lead",                   &lowPtPt_lead_);
  tree->Branch("lowPtEta_lead",                  &lowPtEta_lead_);
  tree->Branch("lowPtPhi_lead",                  &lowPtPhi_lead_);
  tree->Branch("lowPtPtMean_lead",                   &lowPtPtMean_lead_);
  tree->Branch("lowPtEtaMean_lead",                  &lowPtEtaMean_lead_);
  tree->Branch("lowPtPhiMean_lead",                  &lowPtPhiMean_lead_);
  tree->Branch("lowPtMVABWP_lead",               &lowPtMVABWP_lead_);
  tree->Branch("lowPtMVAUnBWP_lead",             &lowPtMVAUnBWP_lead_);

  tree->Branch("lowPtCharge_sublead",               &lowPtCharge_sublead_);
  tree->Branch("lowPtD0_sublead",                   &lowPtD0_sublead_);
  tree->Branch("lowPtDz_sublead",                   &lowPtDz_sublead_);
  tree->Branch("lowPtD0Error_sublead",              &lowPtD0Error_sublead_);
  tree->Branch("lowPtDzError_sublead",              &lowPtDzError_sublead_);
  tree->Branch("lowPtPt_sublead",                   &lowPtPt_sublead_);
  tree->Branch("lowPtEta_sublead",                  &lowPtEta_sublead_);
  tree->Branch("lowPtPhi_sublead",                  &lowPtPhi_sublead_);
  tree->Branch("lowPtPtMean_sublead",                   &lowPtPtMean_sublead_);
  tree->Branch("lowPtEtaMean_sublead",                  &lowPtEtaMean_sublead_);
  tree->Branch("lowPtPhiMean_sublead",                  &lowPtPhiMean_sublead_);
  tree->Branch("lowPtMVABWP_sublead",               &lowPtMVABWP_sublead_);
  tree->Branch("lowPtMVAUnBWP_sublead",             &lowPtMVAUnBWP_sublead_);

  tree->Branch("lowPtSvChi2",                    &lowPtSvChi2_);
  tree->Branch("lowPtSvNDOF",                    &lowPtSvNDOF_);
  tree->Branch("lowPtSvProb",                    &lowPtSvProb_);
  tree->Branch("lowPtSvX",                       &lowPtSvX_);
  tree->Branch("lowPtSvY",                       &lowPtSvY_);
  tree->Branch("lowPtSvZ",                       &lowPtSvZ_);
  tree->Branch("lowPtSvXError",                  &lowPtSvXError_);
  tree->Branch("lowPtSvYError",                  &lowPtSvYError_);
  tree->Branch("lowPtSvZError",                  &lowPtSvZError_);
  tree->Branch("lowPtSvMass",                    &lowPtSvMass_);
  tree->Branch("lowPtSvCtxy",                    &lowPtSvCtxy_);
  tree->Branch("lowPtSvCosAngle",                    &lowPtSvCosAngle_);
  tree->Branch("lowPtSvLxy",                    	   &lowPtSvLxy_);
  tree->Branch("lowPtSvLxyError",                    &lowPtSvLxyError_);

  
}

void ggNtuplizer::fillLowPtElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  // cleanup from previous execution
  lowPtCharge_lead_                  .clear();
  lowPtD0_lead_                      .clear();
  lowPtDz_lead_                      .clear();
  lowPtD0Error_lead_                 .clear();
  lowPtDzError_lead_                 .clear();
  lowPtPt_lead_                      .clear();
  lowPtEta_lead_                     .clear();
  lowPtPhi_lead_                     .clear();
  lowPtPtMean_lead_                      .clear();
  lowPtEtaMean_lead_                     .clear();
  lowPtPhiMean_lead_                     .clear();
  lowPtMVABWP_lead_                  .clear();
  lowPtMVAUnBWP_lead_                .clear();

  lowPtCharge_sublead_                  .clear();
  lowPtD0_sublead_                      .clear();
  lowPtDz_sublead_                      .clear();
  lowPtD0Error_sublead_                 .clear();
  lowPtDzError_sublead_                 .clear();
  lowPtPt_sublead_                      .clear();
  lowPtEta_sublead_                     .clear();
  lowPtPhi_sublead_                     .clear();
  lowPtPtMean_sublead_                      .clear();
  lowPtEtaMean_sublead_                     .clear();
  lowPtPhiMean_sublead_                     .clear();
  lowPtMVABWP_sublead_                  .clear();
  lowPtMVAUnBWP_sublead_                .clear();

  lowPtSvChi2_.clear();
  lowPtSvNDOF_.clear();
  lowPtSvProb_.clear();
  lowPtSvX_.clear();
  lowPtSvY_.clear();
  lowPtSvZ_.clear();
  lowPtSvXError_.clear();
  lowPtSvYError_.clear();
  lowPtSvZError_.clear();
  lowPtSvMass_.clear();
  lowPtSvCtxy_.clear();
  lowPtSvCosAngle_.clear();
  lowPtSvLxy_.clear();
  lowPtSvLxyError_.clear();


  nLowPt_ = 0;

  if (isAOD_) {

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

    for (reco::GsfElectronCollection::const_iterator iEle = lowpTelectronHandle->begin(); iEle != lowpTelectronHandle->end(); ++iEle) {
      if (iEle->gsfTrack()->ptMode() < 0.5) continue;
      if (fabs(iEle->gsfTrack()->vz() - pv.z()) > 1.0) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      for (reco::GsfElectronCollection::const_iterator jEle = iEle+1; jEle != lowpTelectronHandle->end(); ++jEle) {
	if (jEle->pt() < 0.5) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > 1.0) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;
	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv, jpsi_lv;
	iele_lv.SetPtEtaPhiM(iEle->gsfTrack()->ptMode(), iEle->gsfTrack()->etaMode(), iEle->gsfTrack()->phiMode(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->gsfTrack()->ptMode(), jEle->gsfTrack()->etaMode(), jEle->gsfTrack()->phiMode(), pmass);
	//if ((iele_lv + jele_lv).M() < 2.4 || (iele_lv + jele_lv).M() > 3.8) continue;
	//if ((iele_lv + jele_lv).M() > 5.0) continue;
	jpsi_lv = iele_lv + jele_lv;

	auto leadEle = iEle->gsfTrack()->ptMode() > jEle->gsfTrack()->ptMode() ? iEle : jEle;
	auto subleadEle = iEle->gsfTrack()->ptMode() > jEle->gsfTrack()->ptMode() ? jEle : iEle;

	//if ((*ele_mva_wp_biased)[leadEle->gsfTrack()] < 4.6) continue;
	if ((*ele_mva_wp_unbiased)[leadEle->gsfTrack()] < 3.05) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;

	XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	KinematicConstrainedVertexFitter kvFitter;
	RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;

	RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();

	if (DecayVtx->chiSquared() < 0.0) continue;
	//if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 20.0) continue;

	// Accept these 4 tracks as a Bs candidate, fill ntuple

	double ctxy = ((DecayVtx->position().x() - pv.x())*jpsi_lv.Px() + (DecayVtx->position().y() - pv.y())*jpsi_lv.Py())/(pow(jpsi_lv.Pt(),2))*jpsi_lv.M();
	
	math::XYZVector perp(jpsi_lv.Px(), jpsi_lv.Py(), 0.);
	math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	lowPtSvChi2_.push_back(DecayVtx->chiSquared());
	lowPtSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	lowPtSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	lowPtSvX_.push_back(DecayVtx->position().x());
	lowPtSvY_.push_back(DecayVtx->position().y());
	lowPtSvZ_.push_back(DecayVtx->position().z());
	lowPtSvXError_.push_back(DecayVtx->error().cxx());
	lowPtSvYError_.push_back(DecayVtx->error().cyy());
	lowPtSvZError_.push_back(DecayVtx->error().czz());
	lowPtSvMass_.push_back(jpsi_lv.M());
	lowPtSvCtxy_.push_back(ctxy);
	lowPtSvCosAngle_.push_back(cosAngle);
	lowPtSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	lowPtSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	lowPtCharge_lead_          .push_back(leadEle->gsfTrack()->charge());
	lowPtD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
	lowPtDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
	lowPtD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
	lowPtDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
	lowPtPt_lead_              .push_back(leadEle->gsfTrack()->ptMode());
	lowPtEta_lead_             .push_back(leadEle->gsfTrack()->etaMode());
	lowPtPhi_lead_             .push_back(leadEle->gsfTrack()->phiMode());
	lowPtPtMean_lead_              .push_back(leadEle->gsfTrack()->pt());
	lowPtEtaMean_lead_             .push_back(leadEle->gsfTrack()->eta());
	lowPtPhiMean_lead_             .push_back(leadEle->gsfTrack()->phi());

	reco::GsfTrackRef mvaSeed_lead = leadEle->gsfTrack();
	lowPtMVABWP_lead_          .push_back((*ele_mva_wp_biased)[mvaSeed_lead]);
	lowPtMVAUnBWP_lead_     .push_back((*ele_mva_wp_unbiased)[mvaSeed_lead]);

	lowPtCharge_sublead_          .push_back(subleadEle->gsfTrack()->charge());
	lowPtD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
	lowPtDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
	lowPtD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
	lowPtDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
	lowPtPt_sublead_              .push_back(subleadEle->gsfTrack()->ptMode());
	lowPtEta_sublead_             .push_back(subleadEle->gsfTrack()->etaMode());
	lowPtPhi_sublead_             .push_back(subleadEle->gsfTrack()->phiMode());
	lowPtPtMean_sublead_              .push_back(subleadEle->gsfTrack()->pt());
	lowPtEtaMean_sublead_             .push_back(subleadEle->gsfTrack()->eta());
	lowPtPhiMean_sublead_             .push_back(subleadEle->gsfTrack()->phi());

	reco::GsfTrackRef mvaSeed_sublead = subleadEle->gsfTrack();
	lowPtMVABWP_sublead_          .push_back((*ele_mva_wp_biased)[mvaSeed_sublead]);
	lowPtMVAUnBWP_sublead_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_sublead]);

	nLowPt_++;
      
      }
    }

  }

}


