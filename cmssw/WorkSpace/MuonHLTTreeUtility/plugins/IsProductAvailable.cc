// -*- C++ -*-
//
// Package:    IsProductAvailable
// Class:      IsProductAvailable
// 
/**\class IsProductAvailable IsProductAvailable.cc Workspace/IsProductAvailable/src/IsProductAvailable.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu Mar  5 13:37:30 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/GenericHandle.h"

//
// class declaration
//

class IsProductAvailable : public edm::EDFilter {
   public:
      explicit IsProductAvailable(const edm::ParameterSet&);
      ~IsProductAvailable();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      
  std::string theClass;
  edm::InputTag theLabel;
  
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
IsProductAvailable::IsProductAvailable(const edm::ParameterSet& iConfig)
{
  theClass = iConfig.getParameter<std::string>("className");
  theLabel = iConfig.getParameter<edm::InputTag>("src");
}


IsProductAvailable::~IsProductAvailable()
{
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
IsProductAvailable::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::GenericHandle genericHandle(theClass);
  bool canGet = iEvent.getByLabel(theLabel, genericHandle);
  return canGet;
}


//define this as a plug-in
DEFINE_FWK_MODULE(IsProductAvailable);
