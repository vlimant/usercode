// -*- C++ -*-
//
// Package:    ConfigurableAnalysis
// Class:      ConfigurableAnalysis
// 
/**\class ConfigurableAnalysis ConfigurableAnalysis.cc Workspace/ConfigurableAnalysis/src/ConfigurableAnalysis.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Mon Apr 14 11:39:51 CEST 2008
// $Id: ConfigurableAnalysis.cc,v 1.1 2008/05/11 21:24:40 vlimant Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Workspace/ConfigurableAnalysis/interface/Selections.h"
#include "Workspace/ConfigurableAnalysis/interface/Plotter.h"
#include "Workspace/ConfigurableAnalysis/interface/CombinedNTupler.h"

//
// class decleration
//

class ConfigurableAnalysis : public edm::EDFilter {
   public:
      explicit ConfigurableAnalysis(const edm::ParameterSet&);
      ~ConfigurableAnalysis();

   private:
      virtual void beginJob(const edm::EventSetup&);
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  Selections * selections_;
  Plotter * plotter_;
  NTupler * ntupler_;

  std::vector<std::string> flows_;
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
ConfigurableAnalysis::ConfigurableAnalysis(const edm::ParameterSet& iConfig) :
  selections_(0), plotter_(0), ntupler_(0)
{
  VariableHelperInstance::init(iConfig.getParameter<edm::ParameterSet>("Variables"));

  //configure the tools
  selections_ = new Selections(iConfig.getParameter<edm::ParameterSet>("Selections"));
  if (!iConfig.getParameter<edm::ParameterSet>("Plotter").empty())
    plotter_ = new Plotter(iConfig.getParameter<edm::ParameterSet>("Plotter"));

  if (!iConfig.getParameter<edm::ParameterSet>("Ntupler").empty())
    ntupler_ = new CombinedNTupler(iConfig.getParameter<edm::ParameterSet>("Ntupler"));
  
  flows_ = iConfig.getParameter<std::vector<std::string> >("flows");

  //vector of passed selections
  produces<std::vector<bool> >();
  ntupler_->registerleaves(this);
}

ConfigurableAnalysis::~ConfigurableAnalysis()
{
  
}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool ConfigurableAnalysis::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
//void ConfigurableAnalysis::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //will the filter pass or not.
  bool majorGlobalAccept=false;

  VariableHelperInstance::get().update(iEvent,iSetup);

  std::auto_ptr<std::vector<bool> > passedProduct(new std::vector<bool>(flows_.size(),false));
  bool filledOnce=false;  

  // loop the requested selections
  for (Selections::iterator selection=selections_->begin(); selection!=selections_->end();++selection){
    //was this flow of filter actually asked for
    bool skip=true;
    uint iFlow=0;
    for (;iFlow!=flows_.size();++iFlow){if (flows_[iFlow]==selection->name()){skip=false; break;}}
    if (skip) continue;

    //make a specific direction in the plotter
    if (plotter_) plotter_->setDir(selection->name());
    
    // apply individual filters on the event
    std::map<std::string, bool> accept=selection->accept(iEvent);
    
    bool globalAccept=true;
    std::string separator="";
    std::string cumulative="";
    std::string allButOne="allBut_";
    std::string fullAccept="fullAccept_";

    std::string fullContent="fullContent_";
    if (selection->makeContentPlots() && plotter_)
      plotter_->fill(fullContent,iEvent);

    //loop the filters to make cumulative and allButOne job
    for (Selection::iterator filterIt=selection->begin(); filterIt!=selection->end();++filterIt){
      Filter & filter=(**filterIt);
      //      bool lastCut=((filterIt+1)==selection->end());

      //increment the directory name
      cumulative+=separator+filter.name(); separator="_";

      if (accept[filter.name()]){
	//	if (globalAccept && selection->makeCumulativePlots() && !lastCut)
	if (globalAccept && selection->makeCumulativePlots() && plotter_)
	  plotter_->fill(cumulative,iEvent);
      }
      else{
	globalAccept=false;
	// did all the others filter fire
	bool goodForAllButThisOne=true;
	for (std::map<std::string,bool>::iterator decision=accept.begin(); decision!=accept.end();++decision){
	  if (decision->first==filter.name()) continue;
	  if (!decision->second) {
	    goodForAllButThisOne=false;
	    break;}
	}
	if (goodForAllButThisOne && selection->makeAllButOnePlots() && plotter_){
	  plotter_->fill(allButOne+filter.name(),iEvent);
	}
      }
      
    }// loop over the filters in this selection

    if (globalAccept){
      (*passedProduct)[iFlow]=true;
      majorGlobalAccept=true;
      //make final plots only if no cumulative plots
      if (selection->makeFinalPlots() && !selection->makeCumulativePlots() && plotter_)
	plotter_->fill(fullAccept,iEvent);

      //make the ntuple and put it in the event
      if (selection->ntuplize() && !filledOnce && ntupler_){
	ntupler_->fill(iEvent);
	filledOnce=true;}
    }
    
  }//loop the different filter order/number: loop the Selections

  iEvent.put(passedProduct);
  return majorGlobalAccept;
}
   

// ------------ method called once each job just before starting event loop  ------------
void 
ConfigurableAnalysis::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ConfigurableAnalysis::endJob() {
  //print summary tables
  selections_->print();
}


DEFINE_FWK_MODULE(ConfigurableAnalysis);
