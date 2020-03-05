#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          n_gen_pfPho_;
vector<float>  gen_pfphoE_;
vector<float>  gen_pfphoEt_;
vector<float>  gen_pfphoEta_;
vector<float>  gen_pfphoPhi_;

//Cut-based Selection variables

vector<float> gen_pfphoLostInnerHits_;


void ggNtuplizer::branchesGenPFPhotons(TTree* tree) {
  
  tree->Branch("n_gen_pfPho",                    &n_gen_pfPho_);
  tree->Branch("gen_pfphoE",                    &gen_pfphoE_);
  tree->Branch("gen_pfphoEt",                   &gen_pfphoEt_);
  tree->Branch("gen_pfphoEta",                  &gen_pfphoEta_);
  tree->Branch("gen_pfphoPhi",                  &gen_pfphoPhi_);

  tree->Branch("gen_pfphoLostInnerHits",        &gen_pfphoLostInnerHits_);
  
}

void ggNtuplizer::fillGenPFPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  
  gen_pfphoE_                .clear();  
  gen_pfphoEt_                .clear();
  gen_pfphoEta_               .clear();
  gen_pfphoPhi_               .clear();

  gen_pfphoLostInnerHits_.clear();

  n_gen_pfPho_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pckPFCandidates;
  e.getByToken(pckPFCandidateCollection_, pckPFCandidates);

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

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

  for (auto iPho = pckPFCandidates->begin(); iPho != pckPFCandidates->end(); ++iPho) {
    
    double temp_dR = 9999999999.0, bestdR = 9999999999.0;
    for (int i = 0; i < (int) pho_gen_ind.size(); ++i){
      temp_dR = deltaR(iPho->eta(), iPho->phi(), genParticlesHandle->at(pho_gen_ind.at(i)).eta(), genParticlesHandle->at(pho_gen_ind.at(i)).phi());
      if (temp_dR < bestdR) {
        bestdR = temp_dR;
      }
    }

    if (bestdR > 0.1) continue;
 
    if (iPho->pdgId() != 22) continue;

    gen_pfphoE_             .push_back(iPho->energy());
    gen_pfphoEt_            .push_back(iPho->et());
    gen_pfphoEta_           .push_back(iPho->eta());
    gen_pfphoPhi_           .push_back(iPho->phi());

    gen_pfphoLostInnerHits_.push_back(iPho->lostInnerHits());
       
    n_gen_pfPho_++;
  }
}
