// -*- C++ -*-
//
// Package:    NtuplerEDProducer
// Class:      NtuplerEDProducer
// 
/**\class NtuplerEDProducer NtuplerEDProducer.cc Workspace/NtuplerEDProducer/src/NtuplerEDProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Sun May 11 21:12:46 CEST 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"

//
// class decleration
//

class NtuplerEDProducer : public edm::EDProducer {
   public:
      explicit NtuplerEDProducer(const edm::ParameterSet&);
      ~NtuplerEDProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  NTupler * ntupler_;
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
NtuplerEDProducer::NtuplerEDProducer(const edm::ParameterSet& iConfig)
{
  //this Ntupler can work with the InputTagDistributor, but should not be configured as such.
  ntupler_ = new StringBasedNTupler(iConfig.getParameter<edm::ParameterSet>("Ntupler"));
  ntupler_->registerleaves(this);
  produces<double>("dummy");
}


NtuplerEDProducer::~NtuplerEDProducer(){}

// ------------ method called to produce the data  ------------
void
NtuplerEDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ntupler_->fill(iEvent);
  std::auto_ptr<double> v(new double(0));
  iEvent.put(v,"dummy");
}

// ------------ method called once each job just before starting event loop  ------------
void 
NtuplerEDProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtuplerEDProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerEDProducer);
