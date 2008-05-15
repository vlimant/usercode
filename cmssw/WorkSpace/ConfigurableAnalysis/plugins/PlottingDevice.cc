// -*- C++ -*-
//
// Package:    PlottingDevice
// Class:      PlottingDevice
// 
/**\class PlottingDevice PlottingDevice.cc Workspace/PlottingDevice/src/PlottingDevice.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu May 15 14:37:59 CEST 2008
// $Id$
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
//
// class decleration
//

#include "Workspace/ConfigurableAnalysis/interface/Plotter.h"
#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"

class PlottingDevice : public edm::EDAnalyzer {
   public:
      explicit PlottingDevice(const edm::ParameterSet&);
      ~PlottingDevice();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  std::string plotDirectoryName_;
  Plotter * plotter_;
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
PlottingDevice::PlottingDevice(const edm::ParameterSet& iConfig)
{
  //  plotDirectoryName_ = iConfig.getParameter<std::string>("@module_label");
  plotDirectoryName_="PlottingDevice";
  VariableHelperInstance::init(iConfig.getParameter<edm::ParameterSet>("Variables"));
  plotter_ = new Plotter(iConfig.getParameter<edm::ParameterSet>("Plotter"));
}


PlottingDevice::~PlottingDevice(){}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PlottingDevice::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  VariableHelperInstance::get().update(iEvent,iSetup);

  plotter_->setDir(plotDirectoryName_);

  plotter_->fill(plotDirectoryName_, iEvent);
}


void PlottingDevice::beginJob(const edm::EventSetup&){}
void PlottingDevice::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(PlottingDevice);
