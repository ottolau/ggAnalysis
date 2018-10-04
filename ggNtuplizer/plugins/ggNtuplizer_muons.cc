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
bool            matchedTrg_;
vector<pair<float,float>>    muPt_;
vector<pair<float,float>>    muEn_;
vector<pair<float,float>>    muEta_;
vector<pair<float,float>>    muPhi_;
vector<pair<int,int>>      muCharge_;
vector<pair<int,int>>      muType_;
vector<UShort_t> muIDbitfirst_;
vector<UShort_t> muIDbitsecond_;
vector<pair<float,float>>    muD0_;
vector<pair<float,float>>    muDz_;
vector<pair<float,float>>    muSIP_;
vector<pair<float,float>>    muD0Error_;
vector<pair<float,float>>    muDzError_;
vector<pair<float,float>>    muChi2NDF_;
vector<pair<float,float>>    muInnerD0_;
vector<pair<float,float>>    muInnerDz_;
vector<pair<int,int>>      muTrkLayers_;
vector<pair<int,int>>      muPixelLayers_;
vector<pair<int,int>>      muPixelHits_;
vector<pair<int,int>>      muMuonHits_;
vector<pair<int,int>>      muStations_;
vector<pair<int,int>>      muMatches_;
vector<pair<int,int>>      muTrkQuality_;
vector<pair<float,float>>    muIsoTrk_;
vector<pair<float,float>>    muPFChIso_;
vector<pair<float,float>>    muPFPhoIso_;
vector<pair<float,float>>    muPFNeuIso_;
vector<pair<float,float>>    muPFPUIso_;
vector<pair<float,float>>    muPFChIso03_;
vector<pair<float,float>>    muPFPhoIso03_;
vector<pair<float,float>>    muPFNeuIso03_;
vector<pair<float,float>>    muPFPUIso03_;
//vector<pair<float,float>>    muPFMiniIso_;
//vector<ULong64_t> muFiredTrgsfirst_;
vector<bool> muFiredTrgsfirst_;
//vector<ULong64_t> muFiredTrgssecond_;
vector<bool> muFiredTrgssecond_;
vector<ULong64_t> muFiredL1Trgsfirst_;
vector<ULong64_t> muFiredL1Trgssecond_;
vector<pair<float,float>>    muInnervalidFraction_;
vector<pair<float,float>>    musegmentCompatibility_;
vector<pair<float,float>>    muchi2LocalPosition_;
vector<pair<float,float>>    mutrkKink_;
vector<pair<float,float>>    muBestTrkPtError_;
vector<pair<float,float>>    muBestTrkPt_;
vector<pair<int,int>>      muBestTrkType_;

vector<float> muSvChi2_;
vector<float> muSvNDOF_;
vector<float> muSvX_;
vector<float> muSvY_;
vector<float> muSvZ_;
vector<float> muSvXError_;
vector<float> muSvYError_;
vector<float> muSvZError_;
vector<float> muSvMass_;
vector<float> muSvDxySig_;
vector<float> muSvCosAngle_;

vector<pair<int,int>>    uuHadCharge_;
vector<pair<float,float>>  uuHadD0_;
vector<pair<float,float>>  uuHadDz_;
vector<pair<float,float>>  uuHadD0Error_;
vector<pair<float,float>>  uuHadDzError_;
vector<pair<float,float>>  uuHadPt_;
vector<pair<float,float>>  uuHadEta_;
vector<pair<float,float>>  uuHadPhi_;
vector<pair<float,float>>  uuHadVx_;
vector<pair<float,float>>  uuHadVy_;
vector<pair<float,float>>  uuHadVz_;
vector<pair<float,float>>  uuHadEn_;
vector<pair<float,float>>  uuHadTrkChi2_;
vector<pair<float,float>>  uuHadTrkNDOF_;
vector<pair<float,float>>  uuHadTrkNormChi2_;
vector<float>  uuHadJPsiMass_;
vector<float>  uuHadPhiMass_;


void ggNtuplizer::branchesMuons(TTree* tree) {

  tree->Branch("nMu",           &nMu_);
  tree->Branch("matchedTrg",    &matchedTrg_);
  tree->Branch("muPt",          &muPt_);
  tree->Branch("muEn",          &muEn_);
  tree->Branch("muEta",         &muEta_);
  tree->Branch("muPhi",         &muPhi_);
  tree->Branch("muCharge",      &muCharge_);
  tree->Branch("muType",        &muType_);
  tree->Branch("muIDbitfirst",       &muIDbitfirst_);
  tree->Branch("muIDbitsecond",       &muIDbitsecond_);
  tree->Branch("muD0",          &muD0_);
  tree->Branch("muDz",          &muDz_);
  tree->Branch("muSIP",         &muSIP_);
  tree->Branch("muD0Error",          &muD0Error_);
  tree->Branch("muDzError",          &muDzError_);
  tree->Branch("muChi2NDF",     &muChi2NDF_);
  tree->Branch("muInnerD0",     &muInnerD0_);
  tree->Branch("muInnerDz",     &muInnerDz_);
  tree->Branch("muTrkLayers",   &muTrkLayers_);
  tree->Branch("muPixelLayers", &muPixelLayers_);
  tree->Branch("muPixelHits",   &muPixelHits_);
  tree->Branch("muMuonHits",    &muMuonHits_);
  tree->Branch("muStations",    &muStations_);
  tree->Branch("muMatches",     &muMatches_);
  tree->Branch("muTrkQuality",  &muTrkQuality_);
  tree->Branch("muIsoTrk",      &muIsoTrk_);
  tree->Branch("muPFChIso",     &muPFChIso_);
  tree->Branch("muPFPhoIso",    &muPFPhoIso_);
  tree->Branch("muPFNeuIso",    &muPFNeuIso_);
  tree->Branch("muPFPUIso",     &muPFPUIso_);
  tree->Branch("muPFChIso03",   &muPFChIso03_);
  tree->Branch("muPFPhoIso03",  &muPFPhoIso03_);
  tree->Branch("muPFNeuIso03",  &muPFNeuIso03_);
  tree->Branch("muPFPUIso03",   &muPFPUIso03_);
//  tree->Branch("muPFMiniIso",   &muPFMiniIso_);
  tree->Branch("muFiredTrgsfirst",   &muFiredTrgsfirst_);
  tree->Branch("muFiredTrgssecond",   &muFiredTrgssecond_);
  tree->Branch("muFiredL1Trgsfirst", &muFiredL1Trgsfirst_);
  tree->Branch("muFiredL1Trgssecond", &muFiredL1Trgssecond_);
  tree->Branch("muInnervalidFraction",   &muInnervalidFraction_);
  tree->Branch("musegmentCompatibility", &musegmentCompatibility_);
  tree->Branch("muchi2LocalPosition",    &muchi2LocalPosition_);
  tree->Branch("mutrkKink",              &mutrkKink_);
  tree->Branch("muBestTrkPtError",       &muBestTrkPtError_);
  tree->Branch("muBestTrkPt",            &muBestTrkPt_);
  tree->Branch("muBestTrkType",          &muBestTrkType_);
  tree->Branch("muSvChi2",                  &muSvChi2_);
  tree->Branch("muSvNDOF",                  &muSvNDOF_);
  tree->Branch("muSvX",                     &muSvX_);
  tree->Branch("muSvY",                     &muSvY_);
  tree->Branch("muSvZ",                     &muSvZ_);
  tree->Branch("muSvXError",                &muSvXError_);
  tree->Branch("muSvYError",                &muSvYError_);
  tree->Branch("muSvZError",                &muSvZError_);
  tree->Branch("muSvMass",                  &muSvMass_);
  tree->Branch("muSvDxySig",                  &muSvDxySig_);
  tree->Branch("muSvCosAngle",                  &muSvCosAngle_);
  if (!separateVtxFit_) {
    tree->Branch("uuHadCharge",               &uuHadCharge_);
    tree->Branch("uuHadD0",                   &uuHadD0_);
    tree->Branch("uuHadDz",                   &uuHadDz_);
    tree->Branch("uuHadD0Error",              &uuHadD0Error_);
    tree->Branch("uuHadDzError",              &uuHadDzError_);
    tree->Branch("uuHadPt",                   &uuHadPt_);
    tree->Branch("uuHadEta",                  &uuHadEta_);
    tree->Branch("uuHadPhi",                  &uuHadPhi_);
    tree->Branch("uuHadVx",                   &uuHadVx_);
    tree->Branch("uuHadVy",                   &uuHadVy_);
    tree->Branch("uuHadVz",                   &uuHadVz_);
    tree->Branch("uuHadEn",                   &uuHadEn_);
    tree->Branch("uuHadTrkChi2",              &uuHadTrkChi2_);
    tree->Branch("uuHadTrkNDOF",              &uuHadTrkNDOF_);
    tree->Branch("uuHadTrkNormChi2",          &uuHadTrkNormChi2_);
    tree->Branch("uuHadJPsiMass",                  &uuHadJPsiMass_);
    tree->Branch("uuHadPhiMass",                  &uuHadPhiMass_);
  

  }



}

void ggNtuplizer::fillMuons(const edm::Event& e, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  muPt_                  .clear();
  muEn_                  .clear();
  muEta_                 .clear();
  muPhi_                 .clear();
  muCharge_              .clear();
  muType_                .clear();
  muIDbitfirst_               .clear();
  muIDbitsecond_               .clear();
  muD0_                  .clear();
  muDz_                  .clear();
  muSIP_                 .clear();
  muD0Error_             .clear();
  muDzError_             .clear();
  muChi2NDF_             .clear();
  muInnerD0_             .clear();
  muInnerDz_             .clear();
  muTrkLayers_           .clear();
  muPixelLayers_         .clear();
  muPixelHits_           .clear();
  muMuonHits_            .clear();
  muStations_            .clear();
  muMatches_             .clear();
  muTrkQuality_          .clear();
  muIsoTrk_              .clear();
  muPFChIso_             .clear();
  muPFPhoIso_            .clear();
  muPFNeuIso_            .clear();
  muPFPUIso_             .clear();
  muPFChIso03_           .clear();
  muPFPhoIso03_          .clear();
  muPFNeuIso03_          .clear();
  muPFPUIso03_           .clear();
//  muPFMiniIso_           .clear();
  muFiredTrgsfirst_           .clear();
  muFiredTrgssecond_           .clear();
  muFiredL1Trgsfirst_         .clear();
  muFiredL1Trgssecond_         .clear();
  muInnervalidFraction_  .clear();
  musegmentCompatibility_.clear();
  muchi2LocalPosition_   .clear();
  mutrkKink_             .clear();
  muBestTrkPtError_      .clear();
  muBestTrkPt_           .clear();
  muBestTrkType_         .clear();
  muSvChi2_.clear();
  muSvNDOF_.clear();
  muSvX_.clear();
  muSvY_.clear();
  muSvZ_.clear();
  muSvXError_.clear();
  muSvYError_.clear();
  muSvZError_.clear();
  muSvMass_.clear();
  muSvDxySig_.clear();
  muSvCosAngle_.clear();

  uuHadCharge_                  .clear();
  uuHadD0_                      .clear();
  uuHadDz_                      .clear();
  uuHadD0Error_                 .clear();
  uuHadDzError_                 .clear();
  uuHadPt_                      .clear();
  uuHadEta_                     .clear();
  uuHadPhi_                     .clear();
  uuHadVx_                      .clear();
  uuHadVy_                      .clear();
  uuHadVz_                      .clear();
  uuHadEn_	              .clear();
  uuHadTrkChi2_                 .clear();
  uuHadTrkNDOF_                 .clear();
  uuHadTrkNormChi2_             .clear();
  uuHadJPsiMass_		      .clear();
  uuHadPhiMass_		      .clear();


  nMu_ = 0;
  matchedTrg_ = false;

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

  if (!(separateVtxFit_)) {



    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) == 1) matchedTrg_ = true;

      //if (iMu->pt() < 2) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
      if (fabs(iMu->eta()) > 2.5) continue;
      for (edm::View<pat::Muon>::const_iterator jMu = iMu+1; jMu != muonHandle->end(); ++jMu) {
        //if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) != 1 || matchMuonTriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()) != 1) continue;
        //if (jMu->pt() < 2) continue;
        if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;
        if (fabs(jMu->eta()) > 2.5) continue;
        if (iMu->charge() * jMu->charge() > 0.) continue;
        float pmass  = 0.1056583745;
        TLorentzVector iMu_lv, jMu_lv;
        iMu_lv.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), pmass);
        jMu_lv.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), pmass);      
        if (((iMu_lv+jMu_lv)).M() < 2.4 || (iMu_lv+jMu_lv).M() > 3.8) continue;

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

        if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 ||KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
        //KinVtx->movePointerToTheTop();
        //RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();

        //for (pat::PackedCandidateCollection::const_iterator iHad = pfcands->begin(); iHad != pfcands->end(); ++iHad) {
        for (reco::TrackCollection::const_iterator iHad = tracksHandle->begin(); iHad != tracksHandle->end(); ++iHad) {
          if (fabs(iHad->eta()) > 2.5) continue;
          if (fabs(iMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
          if (fabs(jMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(jMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(jMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
          if (fabs(iMu->bestTrack()->vz() - iHad->vz()) > 1) continue;
          if (iHad->normalizedChi2() < 0.0) continue;
          if (iHad->normalizedChi2() > 20.0) continue;

          //for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != pfcands->end(); ++jHad) {
          for (reco::TrackCollection::const_iterator jHad = iHad+1; jHad != tracksHandle->end(); ++jHad) {
            if (iHad->charge()*jHad->charge() > 0.0) continue;
            if (fabs(iMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
            if (fabs(jMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(jMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(jMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
            if (fabs(iMu->bestTrack()->vz() - jHad->vz()) > 1) continue;
            if (iHad->normalizedChi2() < 0.0) continue;
            if (iHad->normalizedChi2() > 20.0) continue;

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
            BsParticles.push_back(pFactory.particle(imuttk, pmass, 0.0, 0, pmasse));
            BsParticles.push_back(pFactory.particle(jmuttk, pmass, 0.0, 0, pmasse));

            KinematicConstrainedVertexFitter BsKvFitter;
            RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
            if (!(BsKinVtx->isValid())) continue;

            RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

            if (DecayVtx->chiSquared() < 0.0) continue;
            if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 30.0) continue;
            TVector3 vtxDisplace(DecayVtx->position().x()-pv.x(), DecayVtx->position().y()-pv.y(), DecayVtx->position().z()-pv.z());
            float cosAngle = vtxDisplace.Dot((iMu_lv + jMu_lv + iHad_lv + jHad_lv).Vect())/(vtxDisplace.Mag()*(iMu_lv + jMu_lv + iHad_lv + jHad_lv).Vect().Mag());
            //if (cosAngle < 0.0) continue;

            float dxy = TMath::Sqrt((DecayVtx->position().x()-pv.x())*(DecayVtx->position().x()-pv.x()) + (DecayVtx->position().y()-pv.y())*(DecayVtx->position().y()-pv.y()));
            float sigmadxy = TMath::Sqrt(DecayVtx->error().cxx()*DecayVtx->error().cxx() + DecayVtx->error().cyy()*DecayVtx->error().cyy());
            //if (dxy/sigmadxy < 2.0) continue;

            muSvChi2_.push_back(DecayVtx->chiSquared());
            muSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
            muSvX_.push_back(DecayVtx->position().x());
            muSvY_.push_back(DecayVtx->position().y());
            muSvZ_.push_back(DecayVtx->position().z());
            muSvXError_.push_back(DecayVtx->error().cxx());
            muSvYError_.push_back(DecayVtx->error().cyy());
            muSvZError_.push_back(DecayVtx->error().czz());
            muSvMass_.push_back((iMu_lv+jMu_lv+iHad_lv+jHad_lv).M());
            muSvDxySig_.push_back(dxy/sigmadxy);
            muSvCosAngle_.push_back(cosAngle);

            uuHadCharge_          .push_back(make_pair(iHad->charge(),jHad->charge()));
            uuHadD0_              .push_back(make_pair(iHad->dxy(pv),jHad->dxy(pv)));
            uuHadDz_              .push_back(make_pair(iHad->dz(pv),jHad->dz(pv)));
            uuHadD0Error_		.push_back(make_pair(iHad->dxyError(),jHad->dxyError()));
            uuHadDzError_		.push_back(make_pair(iHad->dzError(),jHad->dzError()));

            uuHadPt_              .push_back(make_pair(iHad->pt(),jHad->pt()));
            uuHadEta_             .push_back(make_pair(iHad->eta(),jHad->eta()));
            uuHadPhi_             .push_back(make_pair(iHad->phi(),jHad->phi()));
            uuHadVx_		.push_back(make_pair(iHad->vx(),jHad->vx()));
            uuHadVy_		.push_back(make_pair(iHad->vy(),jHad->vy()));
            uuHadVz_		.push_back(make_pair(iHad->vz(),jHad->vz()));

//            uuHadEn_ 	        .push_back(make_pair(iHad->energy(),jHad->energy()));
            uuHadTrkChi2_		.push_back(make_pair(iHad->chi2(),jHad->chi2()));
            uuHadTrkNDOF_		.push_back(make_pair(iHad->ndof(),jHad->ndof()));
            uuHadTrkNormChi2_	.push_back(make_pair(iHad->normalizedChi2(),jHad->normalizedChi2()));
            uuHadJPsiMass_        .push_back((iMu_lv+jMu_lv).M());
            uuHadPhiMass_         .push_back((iHad_lv+jHad_lv).M());

            muPt_    .push_back(make_pair(iMu->pt(),jMu->pt()));
            muEn_    .push_back(make_pair(iMu->energy(),jMu->energy()));
            muEta_   .push_back(make_pair(iMu->eta(),jMu->eta()));
            muPhi_   .push_back(make_pair(iMu->phi(),jMu->phi()));
            muCharge_.push_back(make_pair(iMu->charge(),jMu->charge()));
            muType_  .push_back(make_pair(iMu->type(),jMu->type()));
            muD0_    .push_back(make_pair(iMu->muonBestTrack()->dxy(pv),jMu->muonBestTrack()->dxy(pv)));
            muDz_    .push_back(make_pair(iMu->muonBestTrack()->dz(pv),jMu->muonBestTrack()->dz(pv)));
            muSIP_   .push_back(make_pair(fabs(iMu->dB(pat::Muon::PV3D))/iMu->edB(pat::Muon::PV3D),fabs(jMu->dB(pat::Muon::PV3D))/jMu->edB(pat::Muon::PV3D)));
            muD0Error_    .push_back(make_pair(iMu->muonBestTrack()->dxyError(),jMu->muonBestTrack()->dxyError()));
            muDzError_    .push_back(make_pair(iMu->muonBestTrack()->dzError(),jMu->muonBestTrack()->dzError()));

            UShort_t tmpmuIDbitfirst = 0;

            if (iMu->isLooseMuon())     setbit(tmpmuIDbitfirst, 0);
            if (iMu->isMediumMuon())    setbit(tmpmuIDbitfirst, 1);
            if (iMu->isTightMuon(vtx))  setbit(tmpmuIDbitfirst, 2);
            if (iMu->isSoftMuon(vtx))   setbit(tmpmuIDbitfirst, 3);
            if (iMu->isHighPtMuon(vtx)) setbit(tmpmuIDbitfirst, 4);

            UShort_t tmpmuIDbitsecond = 0;

            if (jMu->isLooseMuon())     setbit(tmpmuIDbitsecond, 0);
            if (jMu->isMediumMuon())    setbit(tmpmuIDbitsecond, 1);
            if (jMu->isTightMuon(vtx))  setbit(tmpmuIDbitsecond, 2);
            if (jMu->isSoftMuon(vtx))   setbit(tmpmuIDbitsecond, 3);
            if (jMu->isHighPtMuon(vtx)) setbit(tmpmuIDbitsecond, 4);

            muIDbitfirst_.push_back(tmpmuIDbitfirst);
            muIDbitsecond_.push_back(tmpmuIDbitsecond);

            //muFiredTrgsfirst_  .push_back(matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()));
            muFiredTrgsfirst_.push_back(muTrkMap[iMu - muonHandle->begin()]);
            muFiredL1Trgsfirst_.push_back(matchL1TriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()));
            //muFiredTrgssecond_  .push_back(matchMuonTriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()));
            muFiredTrgssecond_.push_back(muTrkMap[jMu - muonHandle->begin()]);
            muFiredL1Trgssecond_.push_back(matchL1TriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()));

            muBestTrkPtError_        .push_back(make_pair(iMu->muonBestTrack()->ptError(),jMu->muonBestTrack()->ptError()));
            muBestTrkPt_             .push_back(make_pair(iMu->muonBestTrack()->pt(),jMu->muonBestTrack()->pt()));
            muBestTrkType_           .push_back(make_pair(iMu->muonBestTrackType(),jMu->muonBestTrackType()));
            musegmentCompatibility_  .push_back(make_pair(iMu->segmentCompatibility(),jMu->segmentCompatibility()));
            muchi2LocalPosition_     .push_back(make_pair(iMu->combinedQuality().chi2LocalPosition,jMu->combinedQuality().chi2LocalPosition));
            mutrkKink_               .push_back(make_pair(iMu->combinedQuality().trkKink,jMu->combinedQuality().trkKink));

            const reco::TrackRef glbmufirst = iMu->globalTrack();
            const reco::TrackRef innmufirst = iMu->innerTrack();
            const reco::TrackRef glbmusecond = jMu->globalTrack();
            const reco::TrackRef innmusecond = jMu->innerTrack();

            if (glbmufirst.isNull() || glbmusecond.isNull()) {
              muChi2NDF_ .push_back(make_pair(-99.,-99.));
              muMuonHits_.push_back(make_pair(-99,-99.));
            } else {
              muChi2NDF_.push_back(make_pair(glbmufirst->normalizedChi2(),glbmusecond->normalizedChi2()));
              muMuonHits_.push_back(make_pair(glbmufirst->hitPattern().numberOfValidMuonHits(),glbmusecond->hitPattern().numberOfValidMuonHits()));
            }

            if (innmufirst.isNull() || innmusecond.isNull()) {
              muInnerD0_     .push_back(make_pair(-99.,-99.));
              muInnerDz_     .push_back(make_pair(-99.,-99.));
              muTrkLayers_   .push_back(make_pair(-99,-99.));
              muPixelLayers_ .push_back(make_pair(-99,-99.));
              muPixelHits_   .push_back(make_pair(-99,-99.));
              muTrkQuality_  .push_back(make_pair(-99,-99.));

              muInnervalidFraction_ .push_back(make_pair(-99,-99));
            } else {
              muInnerD0_     .push_back(make_pair(innmufirst->dxy(pv),innmusecond->dxy(pv)));
              muInnerDz_     .push_back(make_pair(innmufirst->dz(pv),innmusecond->dz(pv)));
              muTrkLayers_   .push_back(make_pair(innmufirst->hitPattern().trackerLayersWithMeasurement(),innmusecond->hitPattern().trackerLayersWithMeasurement()));
              muPixelLayers_ .push_back(make_pair(innmufirst->hitPattern().pixelLayersWithMeasurement(),innmusecond->hitPattern().pixelLayersWithMeasurement()));
              muPixelHits_   .push_back(make_pair(innmufirst->hitPattern().numberOfValidPixelHits(),innmusecond->hitPattern().numberOfValidPixelHits()));
              muTrkQuality_  .push_back(make_pair(innmufirst->quality(reco::TrackBase::highPurity),innmusecond->quality(reco::TrackBase::highPurity)));

              muInnervalidFraction_ .push_back(make_pair(innmufirst->validFraction(),innmusecond->validFraction()));
            }

            muStations_   .push_back(make_pair(iMu->numberOfMatchedStations(),jMu->numberOfMatchedStations()));
            muMatches_    .push_back(make_pair(iMu->numberOfMatches(),jMu->numberOfMatches()));
            muIsoTrk_     .push_back(make_pair(iMu->trackIso(),jMu->trackIso()));
            muPFChIso_    .push_back(make_pair(iMu->pfIsolationR04().sumChargedHadronPt,jMu->pfIsolationR04().sumChargedHadronPt));
            muPFPhoIso_   .push_back(make_pair(iMu->pfIsolationR04().sumPhotonEt,jMu->pfIsolationR04().sumPhotonEt));
            muPFNeuIso_   .push_back(make_pair(iMu->pfIsolationR04().sumNeutralHadronEt,jMu->pfIsolationR04().sumNeutralHadronEt));
            muPFPUIso_    .push_back(make_pair(iMu->pfIsolationR04().sumPUPt,jMu->pfIsolationR04().sumPUPt));
            muPFChIso03_  .push_back(make_pair(iMu->pfIsolationR03().sumChargedHadronPt,jMu->pfIsolationR03().sumChargedHadronPt));
            muPFPhoIso03_ .push_back(make_pair(iMu->pfIsolationR03().sumPhotonEt,jMu->pfIsolationR03().sumPhotonEt));
            muPFNeuIso03_ .push_back(make_pair(iMu->pfIsolationR03().sumNeutralHadronEt,jMu->pfIsolationR03().sumNeutralHadronEt));
            muPFPUIso03_  .push_back(make_pair(iMu->pfIsolationR03().sumPUPt,jMu->pfIsolationR03().sumPUPt));
            //muPFMiniIso_  .push_back(make_pair(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*iMu)), 0.05, 0.2, 10., false),getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*jMu)), 0.05, 0.2, 10., false)));

            nMu_++;
          }
        }
      }
    }
  }


  if (separateVtxFit_) {
    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) == 1) matchedTrg_ = true;

      //if (iMu->pt() < 2) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
      if (fabs(iMu->eta()) > 2.5) continue;

      for (edm::View<pat::Muon>::const_iterator jMu = iMu+1; jMu != muonHandle->end(); ++jMu) {
        //if (jMu->pt() < 2) continue;
        if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;
        if (fabs(jMu->eta()) > 2.5) continue;
        if (iMu->charge() * jMu->charge() > 0.) continue;
        TLorentzVector iMu_lv, jMu_lv;
        iMu_lv.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1056583745);
        jMu_lv.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1056583745);      
        if (((iMu_lv+jMu_lv)).M() < 2.4 || (iMu_lv+jMu_lv).M() > 3.8) continue;


        KinematicParticleFactoryFromTransientTrack pFactory;  
        std::vector<RefCountedKinematicParticle> XParticles;
        float pmass  = 0.1056583745;
        float pmasse = 1.e-6 * pmass;

        XParticles.push_back(pFactory.particle(getTransientTrack( *(iMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        XParticles.push_back(pFactory.particle(getTransientTrack( *(jMu->bestTrack()) ), pmass, 0.0, 0, pmasse));

        KinematicConstrainedVertexFitter kvFitter;
        RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

        if (KinVtx->isValid()) {
          RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();
          if (DecayVtx->chiSquared() < 0.0) continue;
          if (DecayVtx->chiSquared() > 10.0) continue;
          TVector3 vtxDisplace(DecayVtx->position().x()-pv.x(), DecayVtx->position().y()-pv.y(), DecayVtx->position().z()-pv.z());
          float cosAngle = vtxDisplace.Dot((iMu_lv+jMu_lv).Vect())/(vtxDisplace.Mag()*(iMu_lv+jMu_lv).Vect().Mag());
          if (cosAngle < 0.7) continue;

          float dxy = TMath::Sqrt((DecayVtx->position().x()-pv.x())*(DecayVtx->position().x()-pv.x()) + (DecayVtx->position().y()-pv.y())*(DecayVtx->position().y()-pv.y()));
          float sigmadxy = TMath::Sqrt(DecayVtx->error().cxx()*DecayVtx->error().cxx() + DecayVtx->error().cyy()*DecayVtx->error().cyy());
          //if (dxy/sigmadxy < 2.0) continue;

          muSvChi2_.push_back(DecayVtx->chiSquared());
          muSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
          muSvX_.push_back(DecayVtx->position().x());
          muSvY_.push_back(DecayVtx->position().y());
          muSvZ_.push_back(DecayVtx->position().z());
          muSvXError_.push_back(DecayVtx->error().cxx());
          muSvYError_.push_back(DecayVtx->error().cyy());
          muSvZError_.push_back(DecayVtx->error().czz());
          muSvMass_.push_back((iMu_lv+jMu_lv).M());
          muSvDxySig_.push_back(dxy/sigmadxy);
          muSvCosAngle_.push_back(cosAngle);

          muPt_    .push_back(make_pair(iMu->pt(),jMu->pt()));
          muEn_    .push_back(make_pair(iMu->energy(),jMu->energy()));
          muEta_   .push_back(make_pair(iMu->eta(),jMu->eta()));
          muPhi_   .push_back(make_pair(iMu->phi(),jMu->phi()));
          muCharge_.push_back(make_pair(iMu->charge(),jMu->charge()));
          muType_  .push_back(make_pair(iMu->type(),jMu->type()));
          muD0_    .push_back(make_pair(iMu->muonBestTrack()->dxy(pv),jMu->muonBestTrack()->dxy(pv)));
          muDz_    .push_back(make_pair(iMu->muonBestTrack()->dz(pv),jMu->muonBestTrack()->dz(pv)));
          muSIP_   .push_back(make_pair(fabs(iMu->dB(pat::Muon::PV3D))/iMu->edB(pat::Muon::PV3D),fabs(jMu->dB(pat::Muon::PV3D))/jMu->edB(pat::Muon::PV3D)));

          UShort_t tmpmuIDbitfirst = 0;

          if (iMu->isLooseMuon())     setbit(tmpmuIDbitfirst, 0);
          if (iMu->isMediumMuon())    setbit(tmpmuIDbitfirst, 1);
          if (iMu->isTightMuon(vtx))  setbit(tmpmuIDbitfirst, 2);
          if (iMu->isSoftMuon(vtx))   setbit(tmpmuIDbitfirst, 3);
          if (iMu->isHighPtMuon(vtx)) setbit(tmpmuIDbitfirst, 4);

          UShort_t tmpmuIDbitsecond = 0;

          if (jMu->isLooseMuon())     setbit(tmpmuIDbitsecond, 0);
          if (jMu->isMediumMuon())    setbit(tmpmuIDbitsecond, 1);
          if (jMu->isTightMuon(vtx))  setbit(tmpmuIDbitsecond, 2);
          if (jMu->isSoftMuon(vtx))   setbit(tmpmuIDbitsecond, 3);
          if (jMu->isHighPtMuon(vtx)) setbit(tmpmuIDbitsecond, 4);

          muIDbitfirst_.push_back(tmpmuIDbitfirst);
          muIDbitsecond_.push_back(tmpmuIDbitsecond);


          //muFiredTrgsfirst_  .push_back(matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()));
          muFiredTrgsfirst_.push_back(muTrkMap[iMu - muonHandle->begin()]);
          muFiredL1Trgsfirst_.push_back(matchL1TriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()));
          //muFiredTrgssecond_  .push_back(matchMuonTriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()));
          muFiredTrgssecond_.push_back(muTrkMap[jMu - muonHandle->begin()]);
          muFiredL1Trgssecond_.push_back(matchL1TriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()));


          muBestTrkPtError_        .push_back(make_pair(iMu->muonBestTrack()->ptError(),jMu->muonBestTrack()->ptError()));
          muBestTrkPt_             .push_back(make_pair(iMu->muonBestTrack()->pt(),jMu->muonBestTrack()->pt()));
          muBestTrkType_           .push_back(make_pair(iMu->muonBestTrackType(),jMu->muonBestTrackType()));
          musegmentCompatibility_  .push_back(make_pair(iMu->segmentCompatibility(),jMu->segmentCompatibility()));
          muchi2LocalPosition_     .push_back(make_pair(iMu->combinedQuality().chi2LocalPosition,jMu->combinedQuality().chi2LocalPosition));
          mutrkKink_               .push_back(make_pair(iMu->combinedQuality().trkKink,jMu->combinedQuality().trkKink));

          const reco::TrackRef glbmufirst = iMu->globalTrack();
          const reco::TrackRef innmufirst = iMu->innerTrack();
          const reco::TrackRef glbmusecond = jMu->globalTrack();
          const reco::TrackRef innmusecond = jMu->innerTrack();

          if (glbmufirst.isNull() || glbmusecond.isNull()) {
            muChi2NDF_ .push_back(make_pair(-99.,-99.));
            muMuonHits_.push_back(make_pair(-99,-99.));
          } else {
            muChi2NDF_.push_back(make_pair(glbmufirst->normalizedChi2(),glbmusecond->normalizedChi2()));
            muMuonHits_.push_back(make_pair(glbmufirst->hitPattern().numberOfValidMuonHits(),glbmusecond->hitPattern().numberOfValidMuonHits()));
          }

          if (innmufirst.isNull() || innmusecond.isNull()) {
            muInnerD0_     .push_back(make_pair(-99.,-99.));
            muInnerDz_     .push_back(make_pair(-99.,-99.));
            muTrkLayers_   .push_back(make_pair(-99,-99.));
            muPixelLayers_ .push_back(make_pair(-99,-99.));
            muPixelHits_   .push_back(make_pair(-99,-99.));
            muTrkQuality_  .push_back(make_pair(-99,-99.));

            muInnervalidFraction_ .push_back(make_pair(-99,-99));
          } else {
            muInnerD0_     .push_back(make_pair(innmufirst->dxy(pv),innmusecond->dxy(pv)));
            muInnerDz_     .push_back(make_pair(innmufirst->dz(pv),innmusecond->dz(pv)));
            muTrkLayers_   .push_back(make_pair(innmufirst->hitPattern().trackerLayersWithMeasurement(),innmusecond->hitPattern().trackerLayersWithMeasurement()));
            muPixelLayers_ .push_back(make_pair(innmufirst->hitPattern().pixelLayersWithMeasurement(),innmusecond->hitPattern().pixelLayersWithMeasurement()));
            muPixelHits_   .push_back(make_pair(innmufirst->hitPattern().numberOfValidPixelHits(),innmusecond->hitPattern().numberOfValidPixelHits()));
            muTrkQuality_  .push_back(make_pair(innmufirst->quality(reco::TrackBase::highPurity),innmusecond->quality(reco::TrackBase::highPurity)));

            muInnervalidFraction_ .push_back(make_pair(innmufirst->validFraction(),innmusecond->validFraction()));
          }

          muStations_   .push_back(make_pair(iMu->numberOfMatchedStations(),jMu->numberOfMatchedStations()));
          muMatches_    .push_back(make_pair(iMu->numberOfMatches(),jMu->numberOfMatches()));
          muIsoTrk_     .push_back(make_pair(iMu->trackIso(),jMu->trackIso()));
          muPFChIso_    .push_back(make_pair(iMu->pfIsolationR04().sumChargedHadronPt,jMu->pfIsolationR04().sumChargedHadronPt));
          muPFPhoIso_   .push_back(make_pair(iMu->pfIsolationR04().sumPhotonEt,jMu->pfIsolationR04().sumPhotonEt));
          muPFNeuIso_   .push_back(make_pair(iMu->pfIsolationR04().sumNeutralHadronEt,jMu->pfIsolationR04().sumNeutralHadronEt));
          muPFPUIso_    .push_back(make_pair(iMu->pfIsolationR04().sumPUPt,jMu->pfIsolationR04().sumPUPt));
          muPFChIso03_  .push_back(make_pair(iMu->pfIsolationR03().sumChargedHadronPt,jMu->pfIsolationR03().sumChargedHadronPt));
          muPFPhoIso03_ .push_back(make_pair(iMu->pfIsolationR03().sumPhotonEt,jMu->pfIsolationR03().sumPhotonEt));
          muPFNeuIso03_ .push_back(make_pair(iMu->pfIsolationR03().sumNeutralHadronEt,jMu->pfIsolationR03().sumNeutralHadronEt));
          muPFPUIso03_  .push_back(make_pair(iMu->pfIsolationR03().sumPUPt,jMu->pfIsolationR03().sumPUPt));
          //muPFMiniIso_  .push_back(make_pair(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*iMu)), 0.05, 0.2, 10., false),getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*jMu)), 0.05, 0.2, 10., false)));

          nMu_++;

        }
      }
    }
  }
}
