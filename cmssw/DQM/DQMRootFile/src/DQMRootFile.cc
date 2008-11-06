// -*- C++ -*-
//
// Package:    DQMRootFile
// Class:      DQMRootFile
// 
/**\class DQMRootFile DQMRootFile.cc DQM/DQMRootFile/src/DQMRootFile.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Fri Nov  2 20:13:09 CET 2007
// $Id: DQMRootFile.cc,v 1.1 2007/11/27 20:51:09 vlimant Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//
// class decleration
//

class DQMRootFile : public edm::EDAnalyzer {
   public:
      explicit DQMRootFile(const edm::ParameterSet&);
      ~DQMRootFile();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  DQMStore * theDBE;
  std::string theDQMRootFileName;
};

//
// constructors and destructor
//
DQMRootFile::DQMRootFile(const edm::ParameterSet& iConfig) :  theDBE(edm::Service<DQMStore>().operator->()), theDQMRootFileName(iConfig.getParameter<std::string>("DQMRootFileName")) {}



DQMRootFile::~DQMRootFile(){}

// ------------ method called to for each event  ------------
void
DQMRootFile::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //I am doing nothing but waiting for the job to end
}


// ------------ method called once each job just before starting event loop  ------------
void 
DQMRootFile::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DQMRootFile::endJob() {
  theDBE->showDirStructure();
  theDBE->save(theDQMRootFileName);  
}

//define this as a plug-in
DEFINE_FWK_MODULE(DQMRootFile);
