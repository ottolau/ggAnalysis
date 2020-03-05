#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          npfPho_;
vector<float>  pfphoE_;
vector<float>  pfphoEt_;
vector<float>  pfphoEta_;
vector<float>  pfphoPhi_;

//Cut-based Selection variables

vector<float> pfphoLostInnerHits_;

void ggNtuplizer::branchesPFPhotons(TTree* tree) {
  
  tree->Branch("npfPho",                    &npfPho_);
  tree->Branch("pfphoE",                    &pfphoE_);
  tree->Branch("pfphoEt",                   &pfphoEt_);
  tree->Branch("pfphoEta",                  &pfphoEta_);
  tree->Branch("pfphoPhi",                  &pfphoPhi_);

  tree->Branch("pfphoLostInnerHits",&pfphoLostInnerHits_);
  

}

void ggNtuplizer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  
  pfphoE_                .clear();  
  pfphoEt_                .clear();
  pfphoEta_               .clear();
  pfphoPhi_               .clear();

  pfphoLostInnerHits_.clear();


  npfPho_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pckPFCandidates;
  e.getByToken(pckPFCandidateCollection_, pckPFCandidates);

  for (auto iPho = pckPFCandidates->begin(); iPho != pckPFCandidates->end(); ++iPho) {
   
    if (iPho->pdgId() != 22) continue;
    
    pfphoE_             .push_back(iPho->energy());
    pfphoEt_            .push_back(iPho->et());
    pfphoEta_           .push_back(iPho->eta());
    pfphoPhi_           .push_back(iPho->phi());
    //std::cout<<iPho->hcalOverEcal()<<std::endl;
    pfphoLostInnerHits_.push_back(iPho->lostInnerHits());
    
    npfPho_++;
  }

}

