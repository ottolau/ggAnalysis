#include "FWCore/MessageLogger/interface/MessageLogger.h"
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
Int_t          nHad_;
vector<pair<int,int>>    hadCharge_;
vector<pair<float,float>>  hadD0_;
vector<pair<float,float>>  hadDz_;
vector<pair<float,float>>  hadD0Error_;
vector<pair<float,float>>  hadDzError_;
vector<pair<float,float>>  hadPt_;
vector<pair<float,float>>  hadEta_;
vector<pair<float,float>>  hadPhi_;
vector<pair<float,float>>  hadVx_;
vector<pair<float,float>>  hadVy_;
vector<pair<float,float>>  hadVz_;
vector<pair<float,float>>  hadEn_;
vector<pair<float,float>>  hadTrkChi2_;
vector<pair<float,float>>  hadTrkNDOF_;
vector<pair<float,float>>  hadTrkNormChi2_;

vector<float> hadSvChi2_;
vector<float> hadSvNDOF_;
vector<float> hadSvX_;
vector<float> hadSvY_;
vector<float> hadSvZ_;
vector<float> hadSvXError_;
vector<float> hadSvYError_;
vector<float> hadSvZError_;
vector<float> hadSvMass_;
vector<float> hadSvDxySig_;
vector<float> hadSvCosAngle_;


void ggNtuplizer::branchesHadrons(TTree* tree) {

  tree->Branch("nHad",                    &nHad_);
  tree->Branch("hadCharge",               &hadCharge_);
  tree->Branch("hadD0",                   &hadD0_);
  tree->Branch("hadDz",                   &hadDz_);
  tree->Branch("hadD0Error",              &hadD0Error_);
  tree->Branch("hadDzError",              &hadDzError_);
  tree->Branch("hadPt",                   &hadPt_);
  tree->Branch("hadEta",                  &hadEta_);
  tree->Branch("hadPhi",                  &hadPhi_);
  tree->Branch("hadVx",                   &hadVx_);
  tree->Branch("hadVy",                   &hadVy_);
  tree->Branch("hadVz",                   &hadVz_);
  tree->Branch("hadEn",                   &hadEn_);
  tree->Branch("hadTrkChi2",              &hadTrkChi2_);
  tree->Branch("hadTrkNDOF",              &hadTrkNDOF_);
  tree->Branch("hadTrkNormChi2",          &hadTrkNormChi2_);
  tree->Branch("hadSvChi2",                  &hadSvChi2_);
  tree->Branch("hadSvNDOF",                  &hadSvNDOF_);
  tree->Branch("hadSvX",                     &hadSvX_);
  tree->Branch("hadSvY",                     &hadSvY_);
  tree->Branch("hadSvZ",                     &hadSvZ_);
  tree->Branch("hadSvXError",                &hadSvXError_);
  tree->Branch("hadSvYError",                &hadSvYError_);
  tree->Branch("hadSvZError",                &hadSvZError_);
  tree->Branch("hadSvMass",                  &hadSvMass_);
  tree->Branch("hadSvDxySig",                  &hadSvDxySig_);
  tree->Branch("hadSvCosAngle",                  &hadSvCosAngle_);


  
}

void ggNtuplizer::fillHadrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv) {
    
  // cleanup from previous execution
  hadCharge_                  .clear();
  hadD0_                      .clear();
  hadDz_                      .clear();
  hadD0Error_                 .clear();
  hadDzError_                 .clear();
  hadPt_                      .clear();
  hadEta_                     .clear();
  hadPhi_                     .clear();
  hadVx_                      .clear();
  hadVy_                      .clear();
  hadVz_                      .clear();
  hadEn_	              .clear();
  hadTrkChi2_                 .clear();
  hadTrkNDOF_                 .clear();
  hadTrkNormChi2_             .clear();

  hadSvChi2_.clear();
  hadSvNDOF_.clear();
  hadSvX_.clear();
  hadSvY_.clear();
  hadSvZ_.clear();
  hadSvXError_.clear();
  hadSvYError_.clear();
  hadSvZError_.clear();
  hadSvMass_.clear();
  hadSvDxySig_.clear();
  hadSvCosAngle_.clear();


  nHad_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  if (!pfcands.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no PF in event";
    return;
  }

  for (pat::PackedCandidateCollection::const_iterator iHad = pfcands->begin(); iHad != pfcands->end(); ++iHad) {
    if (abs(iHad->pdgId()) != 211) continue;
    if (iHad->bestTrack() == nullptr) continue;
    //if (iHad->bestTrack()->normalizedChi2() < 0.0) continue;
    //if (iHad->bestTrack()->normalizedChi2() >= 2.0) continue;
    if (fabs(iHad->eta()) > 2.5) continue;
    if (fabs(iHad->bestTrack()->dz(pv)) > 3.0) continue;
    //if (iHad->bestTrack()->normalizedChi2() == 0.0) cout<<"Chi2: "<<iHad->bestTrack()->chi2()<<" ndof: "<<iHad->bestTrack()->ndof()<<endl;

    for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != pfcands->end(); ++jHad) {
      if (abs(jHad->pdgId()) != 211) continue;
      if (iHad->charge()*jHad->charge() > 0.0) continue;
      // Phi mass window
      TLorentzVector iHad_lv, jHad_lv;
      iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), 0.493677);
      jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), 0.493677);      
      if (((iHad_lv+jHad_lv)).M() < 0.95 || (iHad_lv+jHad_lv).M() > 1.1) continue; 
      if (jHad->bestTrack() == nullptr) continue;
      //if (jHad->bestTrack()->normalizedChi2() < 0.0) continue;
      //if (jHad->bestTrack()->normalizedChi2() >= 2.0) continue;
      if (fabs(jHad->eta()) > 2.5) continue;
      if (fabs(jHad->bestTrack()->dz(pv)) > 3.0) continue;
      //if (iHad_lv.DeltaR(jHad_lv) > 3) continue;


      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;
      float pmass  = 0.493677;
      float pmasse = 1.e-6 * pmass;

      XParticles.push_back(pFactory.particle(getTransientTrack( *(iHad->bestTrack()) ), pmass, 0.0, 0, pmasse));
      XParticles.push_back(pFactory.particle(getTransientTrack( *(jHad->bestTrack()) ), pmass, 0.0, 0, pmasse));

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

      if (KinVtx->isValid()) {
	RefCountedKinematicVertex DecayVtx = KinVtx->currentDecayVertex();
	if (DecayVtx->chiSquared() < 0.0) continue;
	if (DecayVtx->chiSquared() > 10.0) continue;
	TVector3 vtxDisplace(DecayVtx->position().x()-pv.x(), DecayVtx->position().y()-pv.y(), DecayVtx->position().z()-pv.z());
	float cosAngle = vtxDisplace.Dot((iHad_lv+jHad_lv).Vect())/(vtxDisplace.Mag()*(iHad_lv+jHad_lv).Vect().Mag());
	if (cosAngle < 0.7) continue;

	float dxy = TMath::Sqrt((DecayVtx->position().x()-pv.x())*(DecayVtx->position().x()-pv.x()) + (DecayVtx->position().y()-pv.y())*(DecayVtx->position().y()-pv.y()));
	float sigmadxy = TMath::Sqrt(DecayVtx->error().cxx()*DecayVtx->error().cxx() + DecayVtx->error().cyy()*DecayVtx->error().cyy());
	if (dxy/sigmadxy < 3.0) continue;
	
	//cout<<"chiSquared: "<<DecayVtx->chiSquared()<<" dxy/sigmadxy "<<dxy/sigmadxy<<endl;

	hadSvChi2_.push_back(DecayVtx->chiSquared());
	hadSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	hadSvX_.push_back(DecayVtx->position().x());
	hadSvY_.push_back(DecayVtx->position().y());
	hadSvZ_.push_back(DecayVtx->position().z());
	hadSvXError_.push_back(DecayVtx->error().cxx());
	hadSvYError_.push_back(DecayVtx->error().cyy());
	hadSvZError_.push_back(DecayVtx->error().czz());
	hadSvMass_.push_back((iHad_lv+jHad_lv).M());
	hadSvDxySig_.push_back(dxy/sigmadxy);
	hadSvCosAngle_.push_back(cosAngle);

	hadCharge_          .push_back(make_pair(iHad->charge(),jHad->charge()));
	hadD0_              .push_back(make_pair(iHad->bestTrack()->dxy(pv),jHad->bestTrack()->dxy(pv)));
	hadDz_              .push_back(make_pair(iHad->bestTrack()->dz(pv),jHad->bestTrack()->dz(pv)));
	hadD0Error_		.push_back(make_pair(iHad->dxyError(),jHad->dxyError()));
	hadDzError_		.push_back(make_pair(iHad->dzError(),jHad->dzError()));

	hadPt_              .push_back(make_pair(iHad->pt(),jHad->pt()));
	hadEta_             .push_back(make_pair(iHad->eta(),jHad->eta()));
	hadPhi_             .push_back(make_pair(iHad->phi(),jHad->phi()));
	hadVx_		.push_back(make_pair(iHad->vx(),jHad->vx()));
	hadVy_		.push_back(make_pair(iHad->vy(),jHad->vy()));
	hadVz_		.push_back(make_pair(iHad->vz(),jHad->vz()));

	hadEn_ 	        .push_back(make_pair(iHad->energy(),jHad->energy()));
	hadTrkChi2_		.push_back(make_pair(iHad->bestTrack()->chi2(),jHad->bestTrack()->chi2()));
	hadTrkNDOF_		.push_back(make_pair(iHad->bestTrack()->ndof(),jHad->bestTrack()->ndof()));
	hadTrkNormChi2_	.push_back(make_pair(iHad->bestTrack()->normalizedChi2(),jHad->bestTrack()->normalizedChi2()));

	nHad_++;
      }
    }
  }

}

/*
reco::TransientTrack ggNtuplizer::getTransientTrack(const reco::Track& track) {   
    OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T"); 
    reco::TransientTrack transientTrack(track, paramField);
    return transientTrack;
}
*/

