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
Int_t            nHads_;


vector<float> *HadCharge_ = 0;
vector<float> *HadD0_ = 0;
vector<float> *HadDz_ = 0;
vector<float> *HadD0Error_ = 0;
vector<float> *HadDzError_ = 0;
vector<float> *HadPt_ = 0;
vector<float> *HadEta_ = 0;
vector<float> *HadPhi_ = 0;
vector<float> *HadVx_ = 0;
vector<float> *HadVy_ = 0;
vector<float> *HadVz_ = 0;
vector<float> *HadTrkChi2_ = 0;
vector<float> *HadTrkNDOF_ = 0;
vector<float> *HadTrkNormChi2_ = 0;

void ggNtuplizer::branchesV1Tracks(TTree* tree) {

  tree->Branch("nHads",           &nHads_);
  
  tree->Branch("HadCharge", &HadCharge_);
  tree->Branch("HadD0", &HadD0_);
  tree->Branch("HadDz", &HadDz_);
  tree->Branch("HadD0Error", &HadD0Error_);
  tree->Branch("HadDzError", &HadDzError_);
  tree->Branch("HadPt", &HadPt_);
  tree->Branch("HadEta", &HadEta_);
  tree->Branch("HadPhi", &HadPhi_);
  tree->Branch("HadVx", &HadVx_);
  tree->Branch("HadVy", &HadVy_);
  tree->Branch("HadVz", &HadVz_);
  tree->Branch("HadTrkChi2", &HadTrkChi2_);
  tree->Branch("HadTrkNDOF", &HadTrkNDOF_);
  tree->Branch("HadTrkNormChi2", &HadTrkNormChi2_);

}

void ggNtuplizer::fillV1Tracks(const edm::Event& e, math::XYZPoint& pv, reco::Vertex vtx) {

  nHads_ = 0;

  HadCharge_->clear();
  HadD0_->clear();
  HadDz_->clear();
  HadD0Error_->clear();
  HadDzError_->clear();
  HadPt_->clear();
  HadEta_->clear();
  HadPhi_->clear();
  HadVx_->clear();
  HadVy_->clear();
  HadVz_->clear();
  HadTrkChi2_->clear();
  HadTrkNDOF_->clear();
  HadTrkNormChi2_->clear();
  

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  edm::Handle<pat::PackedCandidateCollection> losttracks;
  e.getByToken(lostTracksLabel_, losttracks);

  std::vector<pat::PackedCandidate> alltracks;
  alltracks.reserve(pfcands->size() + losttracks->size());
  alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
  alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

  if (!muonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
    //return;
  }

  vector<bool> muTrkMap;
  muTrkMap = muonTriggerMap(e);

  VertexDistanceXY vertTool;
      

  for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {
    if (iHad->pt() < 1.) continue;
		if (fabs(iHad->pdgId()) != 211) continue;
    if (iHad->bestTrack() == nullptr) continue;
    //if (iHad->bestTrack()->normalizedChi2() > 20.) continue; 
		if (fabs(iHad->eta()) > 2.4) continue;
	  if (fabs(iHad->dz(pv)) > .5) continue;
	  //if (fabs(iHad->bestTrack()->dxy(pv)/iHad->bestTrack()->dxyError()) < 1.) continue;

  	//if (fabs(iMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
  	//if (fabs(iMu->bestTrack()->vz() - iHad->vz()) > 0.5) continue;

    HadCharge_->push_back(iHad->charge());
    HadD0_->push_back(iHad->dxy(pv));
    HadDz_->push_back(iHad->dz(pv));
    HadD0Error_->push_back(iHad->dxyError());
    HadDzError_->push_back(iHad->dzError());
    HadPt_->push_back(iHad->pt());
    HadEta_->push_back(iHad->eta());
    HadPhi_->push_back(iHad->phi());
    HadVx_->push_back(iHad->vx());
    HadVy_->push_back(iHad->vy());
    HadVz_->push_back(iHad->vz());
    HadTrkChi2_->push_back(iHad->bestTrack()->chi2());
    HadTrkNDOF_->push_back(iHad->bestTrack()->ndof());
    HadTrkNormChi2_->push_back(iHad->bestTrack()->normalizedChi2());

    //prb_bestD0bar_->push_back(prb_D0bar);
    //prb_bestB_->push_back(prb_B);
    ++nHads_;
  }
}
