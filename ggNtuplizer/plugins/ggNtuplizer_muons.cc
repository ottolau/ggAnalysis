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

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

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
vector<float>    muVz_lead_;
vector<float>    muD0_lead_;
vector<float>    muDz_lead_;
vector<float>    muSIP_lead_;
vector<float>    muD0Error_lead_;
vector<float>    muDzError_lead_;
vector<float>    muChi2NDF_lead_;
vector<bool>     muFiredTrgs_lead_;
vector<ULong64_t> muFiredL1Trgs_lead_;

vector<float>    muPt_sublead_;
vector<float>    muEn_sublead_;
vector<float>    muEta_sublead_;
vector<float>    muPhi_sublead_;
vector<int>      muCharge_sublead_;
vector<int>      muType_sublead_;
vector<float>    muVz_sublead_;
vector<float>    muD0_sublead_;
vector<float>    muDz_sublead_;
vector<float>    muSIP_sublead_;
vector<float>    muD0Error_sublead_;
vector<float>    muDzError_sublead_;
vector<float>    muChi2NDF_sublead_;
vector<bool>     muFiredTrgs_sublead_;
vector<ULong64_t> muFiredL1Trgs_sublead_;

vector<float> muJpsiSvChi2_;
vector<float> muJpsiSvNDOF_;
vector<float> muJpsiSvProb_;
vector<float> muJpsiSvZ_;
vector<float> muJpsiSvCtxy_;
vector<float> muJpsiSvCosAngle_;
vector<float> muJpsiSvLxy_;
vector<float> muJpsiSvLxyError_;
vector<float> muLambdaSvChi2_;
vector<float> muLambdaSvNDOF_;
vector<float> muLambdaSvProb_;
vector<float> muLambdaSvCtxy_;
vector<float> muLambdaSvCosAngle_;
vector<float> muLambdaSvLxy_;
vector<float> muLambdaSvLxyError_;

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


void ggNtuplizer::branchesMuons(TTree* tree) {

  tree->Branch("nMu",           &nMu_);
  tree->Branch("matchedTrg",    &matchedTrg_);
  tree->Branch("muPt_lead",          &muPt_lead_);
  tree->Branch("muEn_lead",          &muEn_lead_);
  tree->Branch("muEta_lead",         &muEta_lead_);
  tree->Branch("muPhi_lead",         &muPhi_lead_);
  tree->Branch("muCharge_lead",      &muCharge_lead_);
  tree->Branch("muType_lead",        &muType_lead_);
  tree->Branch("muVz_lead",          &muVz_lead_);
  tree->Branch("muD0_lead",          &muD0_lead_);
  tree->Branch("muDz_lead",          &muDz_lead_);
  tree->Branch("muSIP_lead",         &muSIP_lead_);
  tree->Branch("muD0Error_lead",     &muD0Error_lead_);
  tree->Branch("muDzError_lead",     &muDzError_lead_);
  tree->Branch("muFiredTrgs_lead",   &muFiredTrgs_lead_);
  tree->Branch("muFiredL1Trgs_lead", &muFiredL1Trgs_lead_);

  tree->Branch("muPt_sublead",          &muPt_sublead_);
  tree->Branch("muEn_sublead",          &muEn_sublead_);
  tree->Branch("muEta_sublead",         &muEta_sublead_);
  tree->Branch("muPhi_sublead",         &muPhi_sublead_);
  tree->Branch("muCharge_sublead",      &muCharge_sublead_);
  tree->Branch("muType_sublead",        &muType_sublead_);
  tree->Branch("muVz_sublead",          &muVz_sublead_);
  tree->Branch("muD0_sublead",          &muD0_sublead_);
  tree->Branch("muDz_sublead",          &muDz_sublead_);
  tree->Branch("muSIP_sublead",         &muSIP_sublead_);
  tree->Branch("muD0Error_sublead",     &muD0Error_sublead_);
  tree->Branch("muDzError_sublead",     &muDzError_sublead_);
  tree->Branch("muFiredTrgs_sublead",   &muFiredTrgs_sublead_);
  tree->Branch("muFiredL1Trgs_sublead", &muFiredL1Trgs_sublead_);

  tree->Branch("muJpsiSvChi2",                  &muJpsiSvChi2_);
  tree->Branch("muJpsiSvNDOF",                  &muJpsiSvNDOF_);
  tree->Branch("muJpsiSvProb",                  &muJpsiSvProb_);
  tree->Branch("muJpsiSvZ",                     &muJpsiSvZ_);
  tree->Branch("muJpsiSvCtxy",                  &muJpsiSvCtxy_);
  tree->Branch("muJpsiSvCosAngle",              &muJpsiSvCosAngle_);
  tree->Branch("muJpsiSvLxy",                   &muJpsiSvLxy_);
  tree->Branch("muJpsiSvLxyError",              &muJpsiSvLxyError_);
  tree->Branch("muLambdaSvChi2",                  &muLambdaSvChi2_);
  tree->Branch("muLambdaSvNDOF",                  &muLambdaSvNDOF_);
  tree->Branch("muLambdaSvProb",                  &muLambdaSvProb_);
  tree->Branch("muLambdaSvCtxy",                  &muLambdaSvCtxy_);
  tree->Branch("muLambdaSvCosAngle",              &muLambdaSvCosAngle_);
  tree->Branch("muLambdaSvLxy",                   &muLambdaSvLxy_);
  tree->Branch("muLambdaSvLxyError",              &muLambdaSvLxyError_);

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


}

void ggNtuplizer::fillMuons(const edm::Event& e, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  muPt_lead_                  .clear();
  muEn_lead_                  .clear();
  muEta_lead_                 .clear();
  muPhi_lead_                 .clear();
  muCharge_lead_              .clear();
  muType_lead_                .clear();
  muVz_lead_                  .clear();
  muD0_lead_                  .clear();
  muDz_lead_                  .clear();
  muSIP_lead_                 .clear();
  muD0Error_lead_             .clear();
  muDzError_lead_             .clear();
  muFiredTrgs_lead_           .clear();
  muFiredL1Trgs_lead_         .clear();

  muPt_sublead_                  .clear();
  muEn_sublead_                  .clear();
  muEta_sublead_                 .clear();
  muPhi_sublead_                 .clear();
  muCharge_sublead_              .clear();
  muType_sublead_                .clear();
  muVz_sublead_                  .clear();
  muD0_sublead_                  .clear();
  muDz_sublead_                  .clear();
  muSIP_sublead_                 .clear();
  muD0Error_sublead_             .clear();
  muDzError_sublead_             .clear();
  muFiredTrgs_sublead_           .clear();
  muFiredL1Trgs_sublead_         .clear();

  muJpsiSvChi2_.clear();
  muJpsiSvNDOF_.clear();
  muJpsiSvProb_.clear();
  muJpsiSvZ_.clear();
  muJpsiSvCtxy_.clear();
  muJpsiSvCosAngle_.clear();
  muJpsiSvLxy_.clear();
  muJpsiSvLxyError_.clear();
  muLambdaSvChi2_.clear();
  muLambdaSvNDOF_.clear();
  muLambdaSvProb_.clear();
  muLambdaSvCtxy_.clear();
  muLambdaSvCosAngle_.clear();
  muLambdaSvLxy_.clear();
  muLambdaSvLxyError_.clear();

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

  kaonMMCharge_lead_                  .clear();
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

  nMu_ = 0;
  matchedTrg_ = false;

  if (isAOD_) {

    edm::Handle<edm::View<pat::Muon> > muonHandle;
    e.getByToken(muonCollection_, muonHandle);

  //  edm::Handle<pat::PackedCandidateCollection> pfcands;
  //  e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<reco::TrackCollection> tracksHandle;
    e.getByToken(tracklabel_, tracksHandle);

    edm::Handle<reco::VertexCompositeCandidateCollection> lambdasHandle;
    e.getByToken(lambdaLabel_, lambdasHandle);
    //e.getByLabel( edm::InputTag("generalV0Candidates:Lambda"), lambdasHandle );
    //e.getByLabel( edm::InputTag("generalV0Candidates","Lambda","RECO"), lambdasHandle );

    if (!muonHandle.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
      return;
    }

    vector<bool> muTrkMap;
    muTrkMap = muonTriggerMap(e);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) == 1) matchedTrg_ = true;
      //if (iMu->pt() < 2.0) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
      if (fabs(iMu->eta()) > 2.1) continue;
      if (fabs(iMu->vz() - pv.z()) > 1.0) continue;

      for (edm::View<pat::Muon>::const_iterator jMu = iMu+1; jMu != muonHandle->end(); ++jMu) {
        //if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) != 1 || matchMuonTriggerFilters(jMu->pt(), jMu->eta(), jMu->phi()) != 1) continue;
        //if (jMu->pt() < 2.0) continue;
        if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;
        if (fabs(jMu->eta()) > 2.1) continue;
        if (iMu->charge() * jMu->charge() > 0.0) continue;
        if (fabs(jMu->vz() - pv.z()) > 1.0) continue;

        float pmass  = 0.1056583745;
        TLorentzVector iMu_lv, jMu_lv;
        iMu_lv.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), pmass);
        jMu_lv.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), pmass);      
        if (((iMu_lv+jMu_lv)).M() < 2.6 || (iMu_lv+jMu_lv).M() > 3.6) continue;
        //if ((iMu_lv+jMu_lv).M() > 5.0) continue;

        KinematicParticleFactoryFromTransientTrack pFactory;  
        std::vector<RefCountedKinematicParticle> XParticles;
        float pmasse = 1.e-6 * pmass;
        //const reco::TransientTrack imuttk = getTransientTrack( *(iMu->bestTrack()) );
        //const reco::TransientTrack jmuttk = getTransientTrack( *(jMu->bestTrack()) );

        XParticles.push_back(pFactory.particle(getTransientTrack( *(iMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        XParticles.push_back(pFactory.particle(getTransientTrack( *(jMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        //XParticles.push_back(pFactory.particle(imuttk, pmass, 0.0, 0, pmasse));
        //XParticles.push_back(pFactory.particle(jmuttk, pmass, 0.0, 0, pmasse));

        KinematicConstrainedVertexFitter kvFitter;
        RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

        if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 || KinVtx->currentDecayVertex()->chiSquared() > 20.0) continue;
        //if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
        //KinVtx->movePointerToTheTop();
        //RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();
        
        RefCountedKinematicVertex JpsiDecayVtx = KinVtx->currentDecayVertex();

        for (reco::VertexCompositeCandidateCollection::const_iterator iLam = lambdasHandle->begin(); iLam != lambdasHandle->end(); ++iLam) {

          //if (fabs(iLam->vz() - pv.z()) > 1.0) continue;

          // Phi mass window
          //float kpmass = 0.493677;
          //float bsM = 5.3663;
          float protonM = 0.9382720813;
          float pionM = 0.13957018;
          float protonMe = 1.e-6 * protonM;
          float pionMe = 1.e-6 * pionM;

          auto leadMu = iMu->pt() > jMu->pt() ? iMu : jMu;
          auto subleadMu = iMu->pt() > jMu->pt() ? jMu : iMu;

          TLorentzVector jpsi_lv, phi_lv, bs_lv;
          phi_lv.SetPtEtaPhiM(iLam->pt(), iLam->eta(), iLam->phi(), iLam->mass());     
          std::cout<<iLam->mass()<<std::endl;
          jpsi_lv = iMu_lv + jMu_lv;
          bs_lv = jpsi_lv  +  phi_lv;
          //if (phi_lv.M() < 1.095 || phi_lv.M() > 1.135) continue; 
          if (bs_lv.M() < 4.8 || bs_lv.M() > 6.2) continue; 

          double Jpsictxy = ((JpsiDecayVtx->position().x() - pv.x())*jpsi_lv.Px() + (JpsiDecayVtx->position().y() - pv.y())*jpsi_lv.Py())/(pow(jpsi_lv.Pt(),2))*jpsi_lv.M();
          
          math::XYZVector Jpsiperp(jpsi_lv.Px(), jpsi_lv.Py(), 0.);
          math::XYZPoint Jpsidxybs(-1*(pv.x() - JpsiDecayVtx->position().x()), -1*(pv.y() - JpsiDecayVtx->position().y()), 0.);
          math::XYZVector Jpsivperp(Jpsidxybs.x(), Jpsidxybs.y(), 0.);
          double JpsicosAngle = Jpsivperp.Dot(Jpsiperp)/(Jpsivperp.R()*Jpsiperp.R());

          //if (JpsicosAngle < 0.0) continue;

          double Lambdactxy = ((iLam->vx() - pv.x())*phi_lv.Px() + (iLam->vy() - pv.y())*phi_lv.Py())/(pow(phi_lv.Pt(),2))*phi_lv.M();
          
          math::XYZVector Lambdaperp(phi_lv.Px(), phi_lv.Py(), 0.);
          math::XYZPoint Lambdadxybs(-1*(pv.x() - iLam->vx()), -1*(pv.y() - iLam->vy()), 0.);
          math::XYZVector Lambdavperp(Lambdadxybs.x(), Lambdadxybs.y(), 0.);
          double LambdacosAngle = Lambdavperp.Dot(Lambdaperp)/(Lambdavperp.R()*Lambdaperp.R());

          //if (LambdacosAngle < 0.0) continue;

          if (TMath::Prob(JpsiDecayVtx->chiSquared(), JpsiDecayVtx->degreesOfFreedom()) < 0.01) continue;
          //if (TMath::Prob(iLam->vertexChi2(), iLam->vertexNdof()) < 0.01) continue;

          muJpsiSvChi2_.push_back(JpsiDecayVtx->chiSquared());
          muJpsiSvNDOF_.push_back(JpsiDecayVtx->degreesOfFreedom());
          muJpsiSvProb_.push_back(TMath::Prob(JpsiDecayVtx->chiSquared(), JpsiDecayVtx->degreesOfFreedom()));
          muJpsiSvZ_.push_back(JpsiDecayVtx->position().z());
          muJpsiSvCtxy_.push_back(Jpsictxy);
          muJpsiSvCosAngle_.push_back(JpsicosAngle);
          muJpsiSvLxy_.push_back(vertTool.distance(vtx, JpsiDecayVtx.get()->vertexState()).value());
          muJpsiSvLxyError_.push_back(vertTool.distance(vtx, JpsiDecayVtx.get()->vertexState()).error());

          muLambdaSvChi2_.push_back(iLam->vertexChi2());
          muLambdaSvNDOF_.push_back(iLam->vertexNdof());
          muLambdaSvProb_.push_back(TMath::Prob(iLam->vertexChi2(), iLam->vertexNdof()));
          muLambdaSvCtxy_.push_back(Lambdactxy);
          muLambdaSvCosAngle_.push_back(LambdacosAngle);
          //muLambdaSvLxy_.push_back(vertTool.distance(vtx, iLam);
          //muLambdaSvLxyError_.push_back(vertTool.distance(vtx, iLam->vertexCovariance()));

          if (iLam->numberOfDaughters() != 2) {
            std::cout<<"numberOfDaughters!=2, "<<iLam->numberOfDaughters()<<std::endl;
          }

          kaonMMCharge_lead_          .push_back(iLam->daughter(0)->charge());
          //kaonMMD0_lead_              .push_back(iLam->daughter(0)->dxy(pv));
          //kaonMMDz_lead_              .push_back(iLam->daughter(0)->dz(pv));
          //kaonMMD0Error_lead_		  .push_back(iLam->daughter(0)->dxyError());
          //kaonMMDzError_lead_		  .push_back(iLam->daughter(0)->dzError());
          kaonMMPt_lead_              .push_back(iLam->daughter(0)->pt());
          kaonMMEta_lead_             .push_back(iLam->daughter(0)->eta());
          kaonMMPhi_lead_             .push_back(iLam->daughter(0)->phi());
          kaonMMVx_lead_		      .push_back(iLam->daughter(0)->vx());
          kaonMMVy_lead_		      .push_back(iLam->daughter(0)->vy());
          kaonMMVz_lead_		      .push_back(iLam->daughter(0)->vz());
//        kaonMMEn_lead_ 	          .push_back(iLam->daughter(0)->energy());
          //kaonMMTrkChi2_lead_		  .push_back(iLam->daughter(0)->chi2());
          //kaonMMTrkNDOF_lead_		  .push_back(iLam->daughter(0)->ndof());
          //kaonMMTrkNormChi2_lead_	  .push_back(iLam->daughter(0)->normalizedChi2());

          kaonMMCharge_sublead_       .push_back(iLam->daughter(1)->charge());
          //kaonMMD0_sublead_           .push_back(iLam->daughter(1)->dxy(pv));
          //kaonMMDz_sublead_           .push_back(iLam->daughter(1)->dz(pv));
          //kaonMMD0Error_sublead_      .push_back(iLam->daughter(1)->dxyError());
          //kaonMMDzError_sublead_	  .push_back(iLam->daughter(1)->dzError());
          kaonMMPt_sublead_           .push_back(iLam->daughter(1)->pt());
          kaonMMEta_sublead_          .push_back(iLam->daughter(1)->eta());
          kaonMMPhi_sublead_          .push_back(iLam->daughter(1)->phi());
          kaonMMVx_sublead_		      .push_back(iLam->daughter(1)->vx());
          kaonMMVy_sublead_		      .push_back(iLam->daughter(1)->vy());
          kaonMMVz_sublead_		      .push_back(iLam->daughter(1)->vz());
//        kaonMMEn_sublead_ 	      .push_back(iLam->daughter(1)->energy());
          //kaonMMTrkChi2_sublead_      .push_back(iLam->daughter(1)->chi2());
          //kaonMMTrkNDOF_sublead_      .push_back(iLam->daughter(1)->ndof());
          //kaonMMTrkNormChi2_sublead_  .push_back(iLam->daughter(1)->normalizedChi2());

          //bsMMdRmu_            .push_back(iMu_lv.DeltaR(jMu_lv));
          //bsMMdRkaon_          .push_back(iHad_lv.DeltaR(jHad_lv));
          //bsMMdRJpsiPhi_       .push_back((iMu_lv + jMu_lv).DeltaR(iHad_lv + jHad_lv));
          bsMMJpsiMass_        .push_back((jpsi_lv).M());
          bsMMPhiMass_         .push_back((phi_lv).M());
          bsMMBsMass_          .push_back(bs_lv.M());

          muPt_lead_    .push_back(leadMu->pt());
          muEn_lead_    .push_back(leadMu->energy());
          muEta_lead_   .push_back(leadMu->eta());
          muPhi_lead_   .push_back(leadMu->phi());
          muCharge_lead_.push_back(leadMu->charge());
          muType_lead_  .push_back(leadMu->type());
          muVz_lead_    .push_back(leadMu->vz());
          muD0_lead_    .push_back(leadMu->muonBestTrack()->dxy(pv));
          muDz_lead_    .push_back(leadMu->muonBestTrack()->dz(pv));
          muSIP_lead_   .push_back(leadMu->dB(pat::Muon::PV3D)/leadMu->edB(pat::Muon::PV3D));
          muD0Error_lead_    .push_back(leadMu->muonBestTrack()->dxyError());
          muDzError_lead_    .push_back(leadMu->muonBestTrack()->dzError());

          muFiredTrgs_lead_.push_back(muTrkMap[leadMu - muonHandle->begin()]);
          muFiredL1Trgs_lead_.push_back(matchL1TriggerFilters(leadMu->pt(), leadMu->eta(), leadMu->phi()));

          muPt_sublead_    .push_back(subleadMu->pt());
          muEn_sublead_    .push_back(subleadMu->energy());
          muEta_sublead_   .push_back(subleadMu->eta());
          muPhi_sublead_   .push_back(subleadMu->phi());
          muCharge_sublead_.push_back(subleadMu->charge());
          muType_sublead_  .push_back(subleadMu->type());
          muVz_sublead_    .push_back(subleadMu->vz());
          muD0_sublead_    .push_back(subleadMu->muonBestTrack()->dxy(pv));
          muDz_sublead_    .push_back(subleadMu->muonBestTrack()->dz(pv));
          muSIP_sublead_   .push_back(subleadMu->dB(pat::Muon::PV3D)/subleadMu->edB(pat::Muon::PV3D));
          muD0Error_sublead_    .push_back(subleadMu->muonBestTrack()->dxyError());
          muDzError_sublead_    .push_back(subleadMu->muonBestTrack()->dzError());

          muFiredTrgs_sublead_.push_back(muTrkMap[subleadMu - muonHandle->begin()]);
          muFiredL1Trgs_sublead_.push_back(matchL1TriggerFilters(subleadMu->pt(), subleadMu->eta(), subleadMu->phi()));


          nMu_++;
        }
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
        //if (jMu->pt() < 2) continue;
        if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;
        if (fabs(jMu->eta()) > 2.5) continue;
        //if (iMu->charge() * jMu->charge() > 0.) continue;
        if (fabs(jMu->vz() - pv.z()) > 0.5) continue;
        float pmass  = 0.1056583745;
        TLorentzVector iMu_lv, jMu_lv;
        iMu_lv.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), pmass);
        jMu_lv.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), pmass);      
        //if (((iMu_lv+jMu_lv)).M() < 2.4 || (iMu_lv+jMu_lv).M() > 3.8) continue;
        if ((iMu_lv+jMu_lv).M() > 5.0) continue;

        KinematicParticleFactoryFromTransientTrack pFactory;  
        //std::vector<RefCountedKinematicParticle> XParticles;
        float pmasse = 1.e-6 * pmass;
        //const reco::TransientTrack imuttk = getTransientTrack( *(iMu->bestTrack()) );
        //const reco::TransientTrack jmuttk = getTransientTrack( *(jMu->bestTrack()) );

        //XParticles.push_back(pFactory.particle(getTransientTrack( *(iMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        //XParticles.push_back(pFactory.particle(getTransientTrack( *(jMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
        //XParticles.push_back(pFactory.particle(imuttk, pmass, 0.0, 0, pmasse));
        //XParticles.push_back(pFactory.particle(jmuttk, pmass, 0.0, 0, pmasse));

        //KinematicConstrainedVertexFitter kvFitter;
        //RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

        //if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 || KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
        //if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
        //KinVtx->movePointerToTheTop();
        //RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();

        //for (pat::PackedCandidateCollection::const_iterator iHad = pfcands->begin(); iHad != pfcands->end(); ++iHad) {
        for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {

          if (iHad->pt() < 0.8) continue;
          if (iHad->charge() == 0) continue;
          if (abs(iHad->pdgId()) != 211) continue;
          if (iHad->bestTrack() == nullptr) continue;
          if (fabs(iHad->eta()) > 2.5) continue;
          //if (fabs(iMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
          //if (fabs(jMu->bestTrack()->eta() - iHad->eta()) < 0.001 && fabs(jMu->bestTrack()->phi() - iHad->phi()) < 0.001 && fabs(jMu->bestTrack()->pt() - iHad->pt()) < 0.001) continue;
          //if (fabs(iMu->bestTrack()->vz() - iHad->vz()) > 1) continue;
          if (fabs(iHad->vz() - pv.z()) > 0.5) continue;
          //if (iHad->normalizedChi2() < 0.0) continue;
          //if (iHad->normalizedChi2() > 20.0) continue;

          //for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != pfcands->end(); ++jHad) {
          for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != alltracks.end(); ++jHad) {

            if (jHad->pt() < 0.8) continue;
            if (jHad->charge() == 0) continue;
            if (abs(jHad->pdgId()) != 211) continue;
            if (jHad->bestTrack() == nullptr) continue;
            //if (iHad->charge()*jHad->charge() > 0.0) continue;
            //if (fabs(iMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(iMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(iMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
            //if (fabs(jMu->bestTrack()->eta() - jHad->eta()) < 0.001 && fabs(jMu->bestTrack()->phi() - jHad->phi()) < 0.001 && fabs(jMu->bestTrack()->pt() - jHad->pt()) < 0.001) continue;
            //if (fabs(iMu->bestTrack()->vz() - jHad->vz()) > 1) continue;
            if (fabs(jHad->vz() - pv.z()) > 0.5) continue;
            //if (jHad->normalizedChi2() < 0.0) continue;
            //if (jHad->normalizedChi2() > 20.0) continue;

            // Phi mass window
            float kpmass = 0.493677;
            //float bsM = 5.3663;

            TLorentzVector iHad_lv, jHad_lv, bs_lv;
            iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
            jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);     
            bs_lv = iMu_lv + jMu_lv + iHad_lv + jHad_lv;
            //if (((iHad_lv+jHad_lv)).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.06) continue; 
            if ((iHad_lv+jHad_lv).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.10) continue; 
            if ((iMu_lv + jMu_lv + iHad_lv + jHad_lv).M() < 4.5 || (iMu_lv + jMu_lv + iHad_lv + jHad_lv).M() > 6.0) continue;
            if (fabs(jHad->eta()) > 2.5) continue;

            std::vector<RefCountedKinematicParticle> BsParticles;
            float kpmasse = 1.e-6 * pmass;

            BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
            BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
            BsParticles.push_back(pFactory.particle(getTransientTrack( *(iMu->bestTrack()) ), pmass, 0.0, 0, pmasse));
            BsParticles.push_back(pFactory.particle(getTransientTrack( *(jMu->bestTrack()) ), pmass, 0.0, 0, pmasse));

            KinematicConstrainedVertexFitter BsKvFitter;
            RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
            if (!(BsKinVtx->isValid())) continue;

            RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

            if (DecayVtx->chiSquared() < 0.0) continue;
            //if (DecayVtx->chiSquared()/DecayVtx->degreesOfFreedom() > 20.0) continue;

            // Accept these 4 tracks as a Bs candidate, fill ntuple

            auto leadMu = iMu->pt() > jMu->pt() ? iMu : jMu;
            auto subleadMu = iMu->pt() > jMu->pt() ? jMu : iMu;
            auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
            auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

            double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
            
            math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
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
            muSvMass_.push_back((iMu_lv+jMu_lv+iHad_lv+jHad_lv).M());
            muSvCtxy_.push_back(ctxy);
            muSvCosAngle_.push_back(cosAngle);
            muSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
            muSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

            kaonMMCharge_lead_          .push_back(leadHad->charge());
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

            
            muFiredTrgs_lead_.push_back(muTrkMap[leadMu - muonHandle->begin()]);
            muFiredL1Trgs_lead_.push_back(matchL1TriggerFilters(leadMu->pt(), leadMu->eta(), leadMu->phi()));


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


            nMu_++;
          }
        }
        

      }
    }

  }
}
