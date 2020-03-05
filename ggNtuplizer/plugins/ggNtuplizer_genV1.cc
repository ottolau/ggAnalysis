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
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"

using namespace std;

// (local) variables associated with tree branches
Int_t            GenV1nTracks_;
Int_t            GenV1nGenBu_;
Int_t            GenV1nGenD0bar_;
Int_t            GenV1nGenK_;
Int_t            GenV1nGenleadPhi_;
Int_t            GenV1nGensubleadPhi_;
Float_t          GenV1GenBuPt_;
Float_t          GenV1GenBuEta_;
Float_t          GenV1GenBuPhi_;

vector<float> *GenV1Chi2_ = 0;
vector<float> *GenV1NDOF_ = 0;
vector<float> *GenV1Prob_ = 0;
vector<float> *GenV1X_ = 0;
vector<float> *GenV1Y_ = 0;
vector<float> *GenV1Z_ = 0;
vector<float> *GenV1XError_ = 0;
vector<float> *GenV1YError_ = 0;
vector<float> *GenV1ZError_ = 0;
vector<float> *GenV1Mass_ = 0;
vector<float> *GenV1Ctxy_ = 0;
vector<float> *GenV1CosAngle_ = 0;
vector<float> *GenV1Dist_ = 0;
vector<float> *GenV1DistError_ = 0;

vector<float> *GenV1D0barChi2_ = 0;
vector<float> *GenV1D0barNDOF_ = 0;
vector<float> *GenV1D0barProb_ = 0;
vector<float> *GenV1D0barX_ = 0;
vector<float> *GenV1D0barY_ = 0;
vector<float> *GenV1D0barZ_ = 0;
vector<float> *GenV1D0barXError_ = 0;
vector<float> *GenV1D0barYError_ = 0;
vector<float> *GenV1D0barZError_ = 0;
vector<float> *GenV1D0barMass_ = 0;
vector<float> *GenV1D0barLxy_ = 0;
vector<float> *GenV1D0barLxyError_ = 0;


vector<float> *GenV1kaonCharge_ = 0;
vector<float> *GenV1kaonD0_ = 0;
vector<float> *GenV1kaonDz_ = 0;
vector<float> *GenV1kaonD0Error_ = 0;
vector<float> *GenV1kaonDzError_ = 0;
vector<float> *GenV1kaonPt_ = 0;
vector<float> *GenV1kaonEta_ = 0;
vector<float> *GenV1kaonPhi_ = 0;
vector<float> *GenV1kaonVx_ = 0;
vector<float> *GenV1kaonVy_ = 0;
vector<float> *GenV1kaonVz_ = 0;
vector<float> *GenV1kaonTrkChi2_ = 0;
vector<float> *GenV1kaonTrkNDOF_ = 0;
vector<float> *GenV1kaonTrkNormChi2_ = 0;

vector<float> *GenV1dRtrgMuB_ = 0;

vector<float> *GenV1phiCharge_lead_ = 0;
vector<float> *GenV1phiD0_lead_ = 0;
vector<float> *GenV1phiDz_lead_ = 0;
vector<float> *GenV1phiD0Error_lead_ = 0;
vector<float> *GenV1phiDzError_lead_ = 0;
vector<float> *GenV1phiPt_lead_ = 0;
vector<float> *GenV1phiEta_lead_ = 0;
vector<float> *GenV1phiPhi_lead_ = 0;
vector<float> *GenV1phiVx_lead_ = 0;
vector<float> *GenV1phiVy_lead_ = 0;
vector<float> *GenV1phiVz_lead_ = 0;
vector<float> *GenV1phiTrkChi2_lead_ = 0;
vector<float> *GenV1phiTrkNDOF_lead_ = 0;
vector<float> *GenV1phiTrkNormChi2_lead_ = 0;

vector<float> *GenV1phiCharge_sublead_ = 0;
vector<float> *GenV1phiD0_sublead_ = 0;
vector<float> *GenV1phiDz_sublead_ = 0;
vector<float> *GenV1phiD0Error_sublead_ = 0;
vector<float> *GenV1phiDzError_sublead_ = 0;
vector<float> *GenV1phiPt_sublead_ = 0;
vector<float> *GenV1phiEta_sublead_ = 0;
vector<float> *GenV1phiPhi_sublead_ = 0;
vector<float> *GenV1phiVx_sublead_ = 0;
vector<float> *GenV1phiVy_sublead_ = 0;
vector<float> *GenV1phiVz_sublead_ = 0;
vector<float> *GenV1phiTrkChi2_sublead_ = 0;
vector<float> *GenV1phiTrkNDOF_sublead_ = 0;
vector<float> *GenV1phiTrkNormChi2_sublead_ = 0;

void ggNtuplizer::branchesGenV1Tracks(TTree* tree) {

  tree->Branch("GenV1nTracks",           &GenV1nTracks_);
  tree->Branch("GenV1nGenBu",           &GenV1nGenBu_); 
  tree->Branch("GenV1nGenD0bar",           &GenV1nGenD0bar_); 
  tree->Branch("GenV1nGenK",           &GenV1nGenK_); 
  tree->Branch("GenV1nGenleadPhi",           &GenV1nGenleadPhi_); 
  tree->Branch("GenV1nGensubleadPhi",           &GenV1nGensubleadPhi_);
  tree->Branch("GenV1GenBuPt",           &GenV1GenBuPt_);
  tree->Branch("GenV1GenBuEta",           &GenV1GenBuEta_);
  tree->Branch("GenV1GenBuPhi",           &GenV1GenBuPhi_);


  tree->Branch("GenV1Chi2", &GenV1Chi2_);
  tree->Branch("GenV1NDOF", &GenV1NDOF_);
  tree->Branch("GenV1Prob", &GenV1Prob_);
  tree->Branch("GenV1X", &GenV1X_);
  tree->Branch("GenV1Y", &GenV1Y_);
  tree->Branch("GenV1Z", &GenV1Z_);
  tree->Branch("GenV1XError", &GenV1XError_);
  tree->Branch("GenV1YError", &GenV1YError_);
  tree->Branch("GenV1ZError", &GenV1ZError_);
  tree->Branch("GenV1Mass", &GenV1Mass_);
  tree->Branch("GenV1Ctxy", &GenV1Ctxy_);
  tree->Branch("GenV1CosAngle", &GenV1CosAngle_);
  tree->Branch("GenV1Dist", &GenV1Dist_);
  tree->Branch("GenV1DistError", &GenV1DistError_);

  tree->Branch("GenV1D0barChi2", &GenV1D0barChi2_);
  tree->Branch("GenV1D0barNDOF", &GenV1D0barNDOF_);
  tree->Branch("GenV1D0barProb", &GenV1D0barProb_);
  tree->Branch("GenV1D0barX", &GenV1D0barX_);
  tree->Branch("GenV1D0barY", &GenV1D0barY_);
  tree->Branch("GenV1D0barZ", &GenV1D0barZ_);
  tree->Branch("GenV1D0barXError", &GenV1D0barXError_);
  tree->Branch("GenV1D0barYError", &GenV1D0barYError_);
  tree->Branch("GenV1D0barZError", &GenV1D0barZError_);
  tree->Branch("GenV1D0barMass", &GenV1D0barMass_);
  tree->Branch("GenV1D0barLxy", &GenV1D0barLxy_);
  tree->Branch("GenV1D0barLxyError", &GenV1D0barLxyError_);


  tree->Branch("GenV1kaonCharge", &GenV1kaonCharge_);
  tree->Branch("GenV1kaonD0", &GenV1kaonD0_);
  tree->Branch("GenV1kaonDz", &GenV1kaonDz_);
  tree->Branch("GenV1kaonD0Error", &GenV1kaonD0Error_);
  tree->Branch("GenV1kaonDzError", &GenV1kaonDzError_);
  tree->Branch("GenV1kaonPt", &GenV1kaonPt_);
  tree->Branch("GenV1kaonEta", &GenV1kaonEta_);
  tree->Branch("GenV1kaonPhi", &GenV1kaonPhi_);
  tree->Branch("GenV1kaonVx", &GenV1kaonVx_);
  tree->Branch("GenV1kaonVy", &GenV1kaonVy_);
  tree->Branch("GenV1kaonVz", &GenV1kaonVz_);
  tree->Branch("GenV1kaonTrkChi2", &GenV1kaonTrkChi2_);
  tree->Branch("GenV1kaonTrkNDOF", &GenV1kaonTrkNDOF_);
  tree->Branch("GenV1kaonTrkNormChi2", &GenV1kaonTrkNormChi2_);


  tree->Branch("GenV1phiCharge_lead", &GenV1phiCharge_lead_);
  tree->Branch("GenV1phiD0_lead", &GenV1phiD0_lead_);
  tree->Branch("GenV1phiDz_lead", &GenV1phiDz_lead_);
  tree->Branch("GenV1phiD0Error_lead", &GenV1phiD0Error_lead_);
  tree->Branch("GenV1phiDzError_lead", &GenV1phiDzError_lead_);
  tree->Branch("GenV1phiPt_lead", &GenV1phiPt_lead_);
  tree->Branch("GenV1phiEta_lead", &GenV1phiEta_lead_);
  tree->Branch("GenV1phiPhi_lead", &GenV1phiPhi_lead_);
  tree->Branch("GenV1phiVx_lead", &GenV1phiVx_lead_);
  tree->Branch("GenV1phiVy_lead", &GenV1phiVy_lead_);
  tree->Branch("GenV1phiVz_lead", &GenV1phiVz_lead_);
  tree->Branch("GenV1phiTrkChi2_lead", &GenV1phiTrkChi2_lead_);
  tree->Branch("GenV1phiTrkNDOF_lead", &GenV1phiTrkNDOF_lead_);
  tree->Branch("GenV1phiTrkNormChi2_lead", &GenV1phiTrkNormChi2_lead_);

  tree->Branch("GenV1phiCharge_sublead", &GenV1phiCharge_sublead_);
  tree->Branch("GenV1phiD0_sublead", &GenV1phiD0_sublead_);
  tree->Branch("GenV1phiDz_sublead", &GenV1phiDz_sublead_);
  tree->Branch("GenV1phiD0Error_sublead", &GenV1phiD0Error_sublead_);
  tree->Branch("GenV1phiDzError_sublead", &GenV1phiDzError_sublead_);
  tree->Branch("GenV1phiPt_sublead", &GenV1phiPt_sublead_);
  tree->Branch("GenV1phiEta_sublead", &GenV1phiEta_sublead_);
  tree->Branch("GenV1phiPhi_sublead", &GenV1phiPhi_sublead_);
  tree->Branch("GenV1phiVx_sublead", &GenV1phiVx_sublead_);
  tree->Branch("GenV1phiVy_sublead", &GenV1phiVy_sublead_);
  tree->Branch("GenV1phiVz_sublead", &GenV1phiVz_sublead_);
  tree->Branch("GenV1phiTrkChi2_sublead", &GenV1phiTrkChi2_sublead_);
  tree->Branch("GenV1phiTrkNDOF_sublead", &GenV1phiTrkNDOF_sublead_);
  tree->Branch("GenV1phiTrkNormChi2_sublead", &GenV1phiTrkNormChi2_sublead_);

  tree->Branch("GenV1dRtrgMuB", &GenV1dRtrgMuB_);


}

void ggNtuplizer::fillGenV1Tracks(const edm::Event& e, math::XYZPoint& pv, reco::Vertex vtx) {

  GenV1nTracks_ = 0;
  GenV1nGenBu_ = 0;
  GenV1nGenD0bar_ = 0;
  GenV1nGenK_ = 0;
  GenV1nGenleadPhi_ = 0;
  GenV1nGensubleadPhi_ = 0;
  GenV1GenBuPt_ = 0;
  GenV1GenBuEta_ = 0;
  GenV1GenBuPhi_ = 0;

  GenV1Chi2_->clear();
  GenV1NDOF_->clear();
  GenV1Prob_->clear();
  GenV1X_->clear();
  GenV1Y_->clear();
  GenV1Z_->clear();
  GenV1XError_->clear();
  GenV1YError_->clear();
  GenV1ZError_->clear();
  GenV1Mass_->clear();
  GenV1Ctxy_->clear();
  GenV1CosAngle_->clear();
  GenV1Dist_->clear();
  GenV1DistError_->clear();

  GenV1D0barChi2_->clear();
  GenV1D0barNDOF_->clear();
  GenV1D0barProb_->clear();
  GenV1D0barX_->clear();
  GenV1D0barY_->clear();
  GenV1D0barZ_->clear();
  GenV1D0barXError_->clear();
  GenV1D0barYError_->clear();
  GenV1D0barZError_->clear();
  GenV1D0barMass_->clear();
  GenV1D0barLxy_->clear();
  GenV1D0barLxyError_->clear();


  GenV1kaonCharge_->clear();
  GenV1kaonD0_->clear();
  GenV1kaonDz_->clear();
  GenV1kaonD0Error_->clear();
  GenV1kaonDzError_->clear();
  GenV1kaonPt_->clear();
  GenV1kaonEta_->clear();
  GenV1kaonPhi_->clear();
  GenV1kaonVx_->clear();
  GenV1kaonVy_->clear();
  GenV1kaonVz_->clear();
  GenV1kaonTrkChi2_->clear();
  GenV1kaonTrkNDOF_->clear();
  GenV1kaonTrkNormChi2_->clear();
  
  GenV1dRtrgMuB_->clear();

  GenV1phiCharge_lead_->clear();
  GenV1phiD0_lead_->clear();
  GenV1phiDz_lead_->clear();
  GenV1phiD0Error_lead_->clear();
  GenV1phiDzError_lead_->clear();
  GenV1phiPt_lead_->clear();
  GenV1phiEta_lead_->clear();
  GenV1phiPhi_lead_->clear();
  GenV1phiVx_lead_->clear();
  GenV1phiVy_lead_->clear();
  GenV1phiVz_lead_->clear();
  GenV1phiTrkChi2_lead_->clear();
  GenV1phiTrkNDOF_lead_->clear();
  GenV1phiTrkNormChi2_lead_->clear();

  GenV1phiCharge_sublead_->clear();
  GenV1phiD0_sublead_->clear();
  GenV1phiDz_sublead_->clear();
  GenV1phiD0Error_sublead_->clear();
  GenV1phiDzError_sublead_->clear();
  GenV1phiPt_sublead_->clear();
  GenV1phiEta_sublead_->clear();
  GenV1phiPhi_sublead_->clear();
  GenV1phiVx_sublead_->clear();
  GenV1phiVy_sublead_->clear();
  GenV1phiVz_sublead_->clear();
  GenV1phiTrkChi2_sublead_->clear();
  GenV1phiTrkNDOF_sublead_->clear();
  GenV1phiTrkNormChi2_sublead_->clear();

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  edm::Handle<pat::PackedCandidateCollection> losttracks;
  e.getByToken(lostTracksLabel_, losttracks);

  std::vector<pat::PackedCandidate> alltracks;
  alltracks.reserve(pfcands->size() + losttracks->size());
  alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
  alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

  edm::Handle<vector<reco::GenParticle>> genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);


  TLorentzVector GenV1GenBu_lv, GenV1GenD0bar_lv, GenV1GenKaon_lv, GenV1GenleadPhi_lv, GenV1GensubleadPhi_lv;

  for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip !=genParticlesHandle->end(); ++ip){
    if (ip->pdgId() != 521) continue; 
    if (ip->numberOfDaughters() != 2) continue; 
    const Candidate *Bu = &(*ip);
    const Candidate *D0bar, *Kaon, *leadPhi, *subleadPhi;
    if ((ip->daughter(0)->pdgId() == -421) && (ip->daughter(1)->pdgId() == 211)){
      D0bar = ip->daughter(0);
      subleadPhi = ip->daughter(1);
    }else if ((ip->daughter(0)->pdgId() == 211) && (ip->daughter(1)->pdgId() == -421)){
      D0bar = ip->daughter(1);
      subleadPhi = ip->daughter(0);
    }else continue;
    if (D0bar->numberOfDaughters() != 2) continue;
    if ((D0bar->daughter(0)->pdgId() == 321) && (D0bar->daughter(1)->pdgId() == -211)){
      Kaon = D0bar->daughter(0);
      leadPhi = D0bar->daughter(1);
    }else if ((D0bar->daughter(0)->pdgId() == -211) && (D0bar->daughter(1)->pdgId() == 321)){
      Kaon = D0bar->daughter(1);
      leadPhi = D0bar->daughter(0);
    }else continue;
    //cout<<"gen particles?"<<endl; 
    GenV1GenBu_lv.SetPtEtaPhiM(Bu->pt(), Bu->eta(), Bu->phi(), Bu->mass());
    GenV1GenD0bar_lv.SetPtEtaPhiM(D0bar->pt(), D0bar->eta(), D0bar->phi(), D0bar->mass());
    GenV1GenKaon_lv.SetPtEtaPhiM(Kaon->pt(), Kaon->eta(), Kaon->phi(), Kaon->mass());
    GenV1GenleadPhi_lv.SetPtEtaPhiM(leadPhi->pt(), leadPhi->eta(), leadPhi->phi(), leadPhi->mass());
    GenV1GensubleadPhi_lv.SetPtEtaPhiM(subleadPhi->pt(), subleadPhi->eta(), subleadPhi->phi(), subleadPhi->mass()); 
   
    if ((fabs(Kaon->eta()) < 2.4) && (fabs(leadPhi->eta()) < 2.4) &&(fabs(subleadPhi->eta()) < 2.4)){
      GenV1GenBuPt_ = GenV1GenBu_lv.Pt();
      GenV1GenBuEta_ = GenV1GenBu_lv.Eta();
      GenV1GenBuPhi_ = GenV1GenBu_lv.Phi();
    }

  }

  /*vector<bool> IsGen;
  int numGen = 0; 
  for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {
    //if (iHad->pdgId() != 211) continue;
    double tempdR = 99999., bestdR = 99999.;
    if (fabs(iHad->pdgId()) > 211) cout<<iHad->pdgId()<<endl;
    for (int i = 0; i< (int) gen_had_Ind.size(); ++i){
      tempdR = deltaR(iHad->eta(), iHad->phi(), genParticlesHandle->at(gen_had_Ind.at(i)).eta(), genParticlesHandle->at(gen_had_Ind.at(i)).phi());
      if ((tempdR < bestdR) && (iHad->charge() == genParticlesHandle->at(gen_had_Ind.at(i)).charge())) bestdR = tempdR;
    }
    bool genmatch = ((fabs(iHad->pdgId()) == 211) && (bestdR < 0.15));
    IsGen.push_back(genmatch);
    if (genmatch) numGen++;
  }

  cout<<"Number of Cand gen matched : "<<numGen<<endl;
*/
  VertexDistanceXY vertTool;
      
  float pmass = 0.105658;  
  float pmasse = 1.e-8*pmass;  
  float kpmass = 0.493677;
  float kpmasse = 0.000016;
  float phipmass = 0.139570;
  float phipmasse = 0.00000024;
  float D0barpmass = 1.86484;
  //float D0barpmasse = 0.00017;
  //float Dspmass = 1.96847;
  //float Dspmasse = 0.00033;
  //float Dpmass = 1.86962;
  //float Dpmasse = 0.00020;
  //float K0starpmass = 0.89581;
  //float K0starpmasse = 0.00019;

	KinematicParticleFactoryFromTransientTrack pFactory;  	

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


  		for (pat::PackedCandidateCollection::const_iterator jHad = alltracks.begin(); jHad != alltracks.end(); ++jHad) {

        if (jHad == iHad) continue;
		    if (jHad->pt() < 1.) continue;
		    if (fabs(jHad->pdgId()) != 211) continue;
		    if (jHad->bestTrack() == nullptr) continue;
        //if (iHad->bestTrack()->normalizedChi2() > 20.) continue; 
		    if (fabs(jHad->eta()) > 2.4) continue;
  		  if (fabs(jHad->dz(pv)) > .5) continue;
    	  //if (fabs(jHad->bestTrack()->dxy(pv)/jHad->bestTrack()->dxyError()) < 1.) continue;

			 	//if (fabs(iMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
		  	//if (fabs(iMu->bestTrack()->vz() - jHad->vz()) > 0.5) continue;
	        //std::cout<<"all cut passed"<<std::endl;

        if (iHad->charge() != -jHad->charge()) continue;		  	

        pat::PackedCandidateCollection::const_iterator GenV1Kaon, GenV1leadPhi;
        GenV1Kaon = iHad;
        GenV1leadPhi = jHad;

        TLorentzVector GenV1Kaon_lv, GenV1leadPhi_lv;
        GenV1Kaon_lv.SetPtEtaPhiM(GenV1Kaon->pt(), GenV1Kaon->eta(), GenV1Kaon->phi(), kpmass);
        GenV1leadPhi_lv.SetPtEtaPhiM(GenV1leadPhi->pt(), GenV1leadPhi->eta(), GenV1leadPhi->phi(), phipmass);
 
        //GenV1D0bar_lv = GenV1Kaon_lv + GenV1leadPhi_lv;
        //if (GenV1D0bar_lv.Pt() < 5.) continue;

        for (pat::PackedCandidateCollection::const_iterator kHad = alltracks.begin(); kHad != alltracks.end(); ++kHad) {

	  	    if (fabs(kHad->pdgId()) != 211) continue;
          if ((kHad == iHad) | (kHad == jHad)) continue;
          if (kHad->pt() < 1.) continue;
          if (kHad->bestTrack() == nullptr) continue;
          if (fabs(kHad->eta()) > 2.4) continue;
       	  //if (fabs(kHad->bestTrack()->dxy(pv)/kHad->bestTrack()->dxyError()) < 1.) continue;
          //if (kHad->bestTrack()->normalizedChi2() > 20.) continue; 
          if (fabs(kHad->dz(pv)) > .5) continue;
          //if (kHad->dxy(pv)/kHad->dxyError() < 1.) continue;

          if (kHad->charge() != GenV1Kaon->charge()) continue;
          //cout<<"hadrons?"<<endl;
          pat::PackedCandidateCollection::const_iterator GenV1subleadPhi = kHad;

          TLorentzVector GenV1subleadPhi_lv;
          GenV1subleadPhi_lv.SetPtEtaPhiM(GenV1subleadPhi->pt(), GenV1subleadPhi->eta(), GenV1subleadPhi->phi(), phipmass);

          //GenV1b_lv = GenV1D0bar_lv + GenV1subleadPhi_lv;
          std::vector<RefCountedKinematicParticle> GenV1D0bar;
          
          GenV1D0bar.push_back(pFactory.particle(getTransientTrack( *(GenV1leadPhi->bestTrack()) ), phipmass, 0.0, 0, phipmasse));
          GenV1D0bar.push_back(pFactory.particle(getTransientTrack( *(GenV1Kaon->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
        
          KinematicConstrainedVertexFitter GenV1D0barKvFitter;
          RefCountedKinematicTree GenV1D0barKinVtx = GenV1D0barKvFitter.fit(GenV1D0bar); 

          if (!(GenV1D0barKinVtx->isValid())) break;

          GenV1D0barKinVtx->movePointerToTheTop();
          RefCountedKinematicParticle GenV1D0bar_part = GenV1D0barKinVtx->currentParticle(); 

          //AlgebraicVector7 temp1 = GenV1D0bar_part->currentState().kinematicParameters().vector();
          //TLorentzVector temp1_lv, temp2_lv;
          //temp1_lv.SetVectM(TVector3(temp1[3], temp1[4], temp1[5]), temp1[6]);
          //temp2_lv.SetPtEtaPhiM(GenV1subleadPhi->pt(), GenV1subleadPhi->eta(), GenV1subleadPhi->phi(), phipmass);

          RefCountedKinematicVertex GenV1D0barDecayVtx = GenV1D0barKinVtx->currentDecayVertex();

          float	prb_D0bar = TMath::Prob(GenV1D0barDecayVtx->chiSquared(), GenV1D0barDecayVtx->degreesOfFreedom());
    
          //if (GenV1D0barDecayVtx->chiSquared() < 0.0) break;
          //if (fabs(GenV1D0barDecayVtx->position().z() - pv.z()) > 1.) break;
          //if (prb_D0bar < 0.05) break;// || prb_D0bar > 0.975) break;
          AlgebraicVector7 D0bar_vec = GenV1D0bar_part->currentState().kinematicParameters().vector();
          TLorentzVector GenV1D0bar_lv;
          GenV1D0bar_lv.SetVectM(TVector3(D0bar_vec[3], D0bar_vec[4], D0bar_vec[5]), D0bar_vec[6]);
          //if (GenV1D0bar_lv.Pt() < 5.) break;

          std::vector<RefCountedKinematicParticle> GenV1BParticles;
          
          GenV1BParticles.push_back(pFactory.particle(getTransientTrack( *(GenV1subleadPhi->bestTrack()) ), phipmass, 0.0, 0, phipmasse));
          GenV1BParticles.push_back(GenV1D0bar_part);
          KinematicConstrainedVertexFitter GenV1BKvFitter;
          RefCountedKinematicTree GenV1BKinVtx = GenV1BKvFitter.fit(GenV1BParticles);
          if (!(GenV1BKinVtx->isValid())) continue;

          GenV1BKinVtx->movePointerToTheTop();
          RefCountedKinematicVertex GenV1DecayVtx = GenV1BKinVtx->currentDecayVertex();

          double prb_B = TMath::Prob(GenV1DecayVtx->chiSquared(), GenV1DecayVtx->degreesOfFreedom());

          AlgebraicVector7 b_vec = GenV1BKinVtx->currentParticle()->currentState().kinematicParameters().vector();
          TLorentzVector GenV1b_lv;
          GenV1b_lv.SetVectM(TVector3(b_vec[3], b_vec[4], b_vec[5]), b_vec[6]);

          //if (GenV1b_lv.DeltaR(GenV1GenBu_lv) > 0.01) continue;


          /*GenV1BKinVtx->movePointerToTheFirstChild();
          while (!(GenV1BKinVtx->isEmpty()) && (GenV1BKinVtx->currentParticle()->correspondingTree() == nullptr)){
            GenV1BKinVtx->movePointerToTheNextChild(); 
          }
          //cout<<!(GenV1BKinVtx->currentParticle()->correspondingTree() == nullptr)<<endl;
          GenV1D0barKinVtx = GenV1BKinVtx->currentParticle()->correspondingTree();
          GenV1D0barKinVtx->movePointerToTheTop();
          if (!(GenV1D0barKinVtx->isValid())) break;
          
          std::vector<RefCountedKinematicParticle> GenV1BParticles_refitted = GenV1BKinVtx->daughterParticles();
          RefCountedKinematicParticle GenV1D0bar_part_refitted = GenV1BParticles_refitted.at(1); 
          RefCountedKinematicTree temp = GenV1D0bar_part_refitted->correspondingTree();
          temp->movePointerToTheTop();
          RefCountedKinematicVertex GenV1D0barDecayVtx = temp->currentDecayVertex();

          float	prb_D0bar = TMath::Prob(GenV1D0barDecayVtx->chiSquared(), GenV1D0barDecayVtx->degreesOfFreedom());
    
          //if (GenV1D0barDecayVtx->chiSquared() < 0.0) break;
          if (fabs(GenV1D0barDecayVtx->position().z() - pv.z()) > 1.) break;
          if (prb_D0bar < 0.05) break;// || prb_D0bar > 0.975) break;
          AlgebraicVector7 D0bar_vec = GenV1D0bar_part_refitted->currentState().kinematicParameters().vector();
          TLorentzVector GenV1D0bar_lv;
          GenV1D0bar_lv.SetVectM(TVector3(D0bar_vec[3], D0bar_vec[4], D0bar_vec[5]), D0bar_vec[6]);
          if (GenV1D0bar_lv.Pt() < 5.) break;*/
          //if (prb_D0bar < prb_bestD0bar) break;
          //if (fabs(prb_D0bar > fabs(prb_bestD0bar - 0.5)) break;
          //cout<<"1: "<<temp1_lv.Pt()<<" "<<temp1_lv.Eta()<<" "<<temp1_lv.Phi()<<" "<<temp1_lv.M()<<endl;
          //cout<<"2: "<<temp2_lv.Pt()<<" "<<temp2_lv.Eta()<<" "<<temp2_lv.Phi()<<" "<<temp2_lv.M()<<endl;


          //cout<<"3: "<<GenV1D0bar_lv.Pt()<<" "<<GenV1D0bar_lv.Eta()<<" "<<GenV1D0bar_lv.Phi()<<" "<<GenV1D0bar_lv.M()<<endl;
          //cout<<endl;


          //if (fabs(GenV1DecayVtx->position().z() - pv.z()) > 1.) continue;
          //if (GenV1DecayVtx->chiSquared() < 0.0) continue;
          //if (GenV1DecayVtx->chiSquared()/GenV1DecayVtx->degreesOfFreedom() > 20.0) continue;
          //if (!(prb_B > 0.)) continue;
          //if (prb_B < prb_bestB) continue;
          //if (fabs(prb_B - 0.5) > fabs(prb_bestB -0.5)) continue;

          double GenV1ctxy = ((GenV1DecayVtx->position().x() - pv.x())*GenV1b_lv.Px() + (GenV1DecayVtx->position().y() - pv.y())*GenV1b_lv.Py())/(pow(GenV1b_lv.Pt(),2))*GenV1b_lv.M();
          math::XYZVector GenV1perp(GenV1b_lv.Px(), GenV1b_lv.Py(), 0.);
          math::XYZPoint GenV1dxybs(-1*(pv.x() - GenV1DecayVtx->position().x()), -1*(pv.y() - GenV1DecayVtx->position().y()), 0.);
          math::XYZVector GenV1vperp(GenV1dxybs.x(), GenV1dxybs.y(), 0.);
          double GenV1cosAngle = GenV1vperp.Dot(GenV1perp)/(GenV1vperp.R()*GenV1perp.R());
         
          //if (GenV1cosAngle < 0.999) continue;

          double GenV1Lxy= vertTool.distance(vtx, GenV1DecayVtx.get()->vertexState()).value();
          double GenV1LxyError= vertTool.distance(vtx, GenV1DecayVtx.get()->vertexState()).error();
    
          //if (GenV1Lxy/GenV1LxyError < 5.) continue;

          GenV1Chi2_->push_back(GenV1DecayVtx->chiSquared());
          GenV1NDOF_->push_back(GenV1DecayVtx->degreesOfFreedom());
          GenV1Prob_->push_back(TMath::Prob(GenV1DecayVtx->chiSquared(), GenV1DecayVtx->degreesOfFreedom()));
          GenV1X_->push_back(GenV1DecayVtx->position().x());
          GenV1Y_->push_back(GenV1DecayVtx->position().y());
          GenV1Z_->push_back(GenV1DecayVtx->position().z());
          GenV1XError_->push_back(GenV1DecayVtx->error().cxx());
          GenV1YError_->push_back(GenV1DecayVtx->error().cyy());
          GenV1ZError_->push_back(GenV1DecayVtx->error().czz());
          GenV1Mass_->push_back(GenV1b_lv.M());
          GenV1Ctxy_->push_back(GenV1ctxy);
          GenV1CosAngle_->push_back(GenV1cosAngle);
          GenV1Dist_->push_back(vertTool.distance(vtx, GenV1DecayVtx.get()->vertexState()).value());
          GenV1DistError_->push_back(vertTool.distance(vtx, GenV1DecayVtx.get()->vertexState()).error());

          GenV1D0barChi2_->push_back(GenV1D0barDecayVtx->chiSquared());
          GenV1D0barNDOF_->push_back(GenV1D0barDecayVtx->degreesOfFreedom());
          GenV1D0barProb_->push_back(TMath::Prob(GenV1D0barDecayVtx->chiSquared(), GenV1D0barDecayVtx->degreesOfFreedom()));
          GenV1D0barX_->push_back(GenV1D0barDecayVtx->position().x());
          GenV1D0barY_->push_back(GenV1D0barDecayVtx->position().y());
          GenV1D0barZ_->push_back(GenV1D0barDecayVtx->position().z());
          GenV1D0barXError_->push_back(GenV1D0barDecayVtx->error().cxx());
          GenV1D0barYError_->push_back(GenV1D0barDecayVtx->error().cyy());
          GenV1D0barZError_->push_back(GenV1D0barDecayVtx->error().czz());
          GenV1D0barMass_->push_back(GenV1D0bar_lv.M());
          GenV1D0barLxy_->push_back(vertTool.distance(vtx, GenV1D0barDecayVtx.get()->vertexState()).value());
          GenV1D0barLxyError_->push_back(vertTool.distance(vtx, GenV1D0barDecayVtx.get()->vertexState()).error());

          GenV1kaonCharge_->push_back(GenV1Kaon->charge());
          GenV1kaonD0_->push_back(GenV1Kaon->dxy(pv));
          GenV1kaonDz_->push_back(GenV1Kaon->dz(pv));
          GenV1kaonD0Error_->push_back(GenV1Kaon->dxyError());
          GenV1kaonDzError_->push_back(GenV1Kaon->dzError());
          GenV1kaonPt_->push_back(GenV1Kaon->pt());
          GenV1kaonEta_->push_back(GenV1Kaon->eta());
          GenV1kaonPhi_->push_back(GenV1Kaon->phi());
          GenV1kaonVx_->push_back(GenV1Kaon->vx());
          GenV1kaonVy_->push_back(GenV1Kaon->vy());
          GenV1kaonVz_->push_back(GenV1Kaon->vz());
          GenV1kaonTrkChi2_->push_back(GenV1Kaon->bestTrack()->chi2());
          GenV1kaonTrkNDOF_->push_back(GenV1Kaon->bestTrack()->ndof());
          GenV1kaonTrkNormChi2_->push_back(GenV1Kaon->bestTrack()->normalizedChi2());

          GenV1phiCharge_lead_->push_back(GenV1leadPhi->charge());
          GenV1phiD0_lead_->push_back(GenV1leadPhi->dxy(pv));
          GenV1phiDz_lead_->push_back(GenV1leadPhi->dz(pv));
          GenV1phiD0Error_lead_->push_back(GenV1leadPhi->dxyError());
          GenV1phiDzError_lead_->push_back(GenV1leadPhi->dzError());
          GenV1phiPt_lead_->push_back(GenV1leadPhi->pt());
          GenV1phiEta_lead_->push_back(GenV1leadPhi->eta());
          GenV1phiPhi_lead_->push_back(GenV1leadPhi->phi());
          GenV1phiVx_lead_->push_back(GenV1leadPhi->vx());
          GenV1phiVy_lead_->push_back(GenV1leadPhi->vy());
          GenV1phiVz_lead_->push_back(GenV1leadPhi->vz());
          GenV1phiTrkChi2_lead_->push_back(GenV1leadPhi->bestTrack()->chi2());
          GenV1phiTrkNDOF_lead_->push_back(GenV1leadPhi->bestTrack()->ndof());
          GenV1phiTrkNormChi2_lead_->push_back(GenV1leadPhi->bestTrack()->normalizedChi2());

          GenV1phiCharge_sublead_->push_back(GenV1subleadPhi->charge());
          GenV1phiD0_sublead_->push_back(GenV1subleadPhi->dxy(pv));
          GenV1phiDz_sublead_->push_back(GenV1subleadPhi->dz(pv));
          GenV1phiD0Error_sublead_->push_back(GenV1subleadPhi->dxyError());
          GenV1phiDzError_sublead_->push_back(GenV1subleadPhi->dzError());
          GenV1phiPt_sublead_->push_back(GenV1subleadPhi->pt());
          GenV1phiEta_sublead_->push_back(GenV1subleadPhi->eta());
          GenV1phiPhi_sublead_->push_back(GenV1subleadPhi->phi());
          GenV1phiVx_sublead_->push_back(GenV1subleadPhi->vx());
          GenV1phiVy_sublead_->push_back(GenV1subleadPhi->vy());
          GenV1phiVz_sublead_->push_back(GenV1subleadPhi->vz());
          GenV1phiTrkChi2_sublead_->push_back(GenV1subleadPhi->bestTrack()->chi2());
          GenV1phiTrkNDOF_sublead_->push_back(GenV1subleadPhi->bestTrack()->ndof());
          GenV1phiTrkNormChi2_sublead_->push_back(GenV1subleadPhi->bestTrack()->normalizedChi2());

          //bool Kaon_gen = IsGen.at(iHad - alltracks.begin());
          //bool leadPhi_gen = IsGen.at(jHad - alltracks.begin());
          //bool subleadPhi_gen = IsGen.at(kHad - alltracks.begin());

          if (GenV1Kaon_lv.DeltaR(GenV1GenKaon_lv) < 0.1)  ++GenV1nGenK_;
          if (GenV1leadPhi_lv.DeltaR(GenV1GenleadPhi_lv) < 0.1)  ++GenV1nGenleadPhi_;
          if (GenV1subleadPhi_lv.DeltaR(GenV1GensubleadPhi_lv) < 0.1)  ++GenV1nGensubleadPhi_;
          if (GenV1D0bar_lv.DeltaR(GenV1GenD0bar_lv) < 0.1) ++GenV1nGenD0bar_;
          if (GenV1b_lv.DeltaR(GenV1GenBu_lv) < 0.1) ++GenV1nGenBu_;

          //prb_bestD0bar_->push_back(prb_D0bar);
          //prb_bestB_->push_back(prb_B);
          ++GenV1nTracks_;
	      }
			}
  }
}
