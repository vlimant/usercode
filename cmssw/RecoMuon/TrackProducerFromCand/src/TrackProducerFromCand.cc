// -*- C++ -*-
//
// Package:    TrackProducerFromCand
// Class:      TrackProducerFromCand
// 
/**\class TrackProducerFromCand TrackProducerFromCand.cc RecoMuon/TrackProducerFromCand/src/TrackProducerFromCand.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Sat Oct 13 03:03:13 CEST 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMuon/TrackingTools/interface/MuonErrorMatrix.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"

#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

//
// class decleration
//

class TrackProducerFromCand : public edm::EDProducer {
   public:
      explicit TrackProducerFromCand(const edm::ParameterSet&);
      ~TrackProducerFromCand();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  edm::InputTag sourceTag_;
  std::string instanceLabel_;
  edm::ParameterSet assignerPset;
  MuonErrorMatrix * errorMatrix_Assigner_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TrackProducerFromCand::TrackProducerFromCand(const edm::ParameterSet& iConfig){
  sourceTag_=iConfig.getParameter<edm::InputTag>("sourceTag");
  instanceLabel_=iConfig.getParameter<std::string>("instanceName");

  assignerPset = iConfig.getParameter<edm::ParameterSet>("errorMatrix_Assigner");

  //very basic, do not produce rechits
  produces<reco::TrackCollection>(instanceLabel_);
}

TrackProducerFromCand::~TrackProducerFromCand(){}

void TrackProducerFromCand::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  //get the magnetic field
  edm::ESHandle<MagneticField> field;
  iSetup.get<IdealMagneticFieldRecord>().get(field);
  
  Handle<edm::View< reco::Candidate> > pCand;
  iEvent.getByLabel(sourceTag_ ,pCand);
  
  std::auto_ptr<reco::TrackCollection> pOut(new reco::TrackCollection());
  
  uint nC=pCand->size();
  uint iC=0;
  for (;iC!=nC;iC++){
    const reco::Candidate & c = (*pCand)[iC];

    reco::TrackBase::Point position(c.vx(), c.vy(), c.vz());
    reco::TrackBase::Vector momentum(c.px(), c.py(), c.pz());
    int charge = (int) c.charge();
    reco::TrackBase::CovarianceMatrix covariance_matrix;
    // do we have a clue what the error is ?
    if (errorMatrix_Assigner_){
      GlobalVector mom(c.px(), c.py(), c.pz());
      CurvilinearTrajectoryError cpe=errorMatrix_Assigner_->get(mom);
      const AlgebraicSymMatrix55 & mat= cpe.matrix();
      for (uint i=0;i!=5;i++){
	for (uint j=0;j!=5;j++){
	  covariance_matrix(i,j)=mat(i,j);
	}}
    }
    double chi2=-1;
    int ndof = -1;
    pOut->push_back( reco::Track(chi2, ndof, position, momentum, charge, covariance_matrix));
  }
  
  iEvent.put(pOut);
}

void TrackProducerFromCand::beginJob(const edm::EventSetup&){
  if (!assignerPset.empty()){
    errorMatrix_Assigner_ = new MuonErrorMatrix(assignerPset);}
}
void TrackProducerFromCand::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackProducerFromCand);
