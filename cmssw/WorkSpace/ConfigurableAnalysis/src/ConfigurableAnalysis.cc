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

#include "Workspace/ConfigurableAnalysis/interface/Selections.h"

//
// class decleration
//

class ConfigurableAnalysis : public edm::EDProducer {
   public:
      explicit ConfigurableAnalysis(const edm::ParameterSet&);
      ~ConfigurableAnalysis();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  Selections * selections_;
  Plotter * plotter_;
  Ntupler * ntupler_;

  //  std::map<std::string, std::map<std::string, uint> > counts_;
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
ConfigurableAnalysis::ConfigurableAnalysis(const edm::ParameterSet& iConfig)
{
  //configure the tools
  //  retriever_->configure(iConfig);
  selections_ = new Selections(iConfig.getParameter<edm::ParameterSet>("Selections"));
  plotter_ = new Plotter(iConfig.getParameter<edm::ParameterSet>("Plotter"));

  ntupler_ = new Ntupler(iConfig.getParameter<edm::ParameterSet>("Ntupler"));
  //register the configurable ntuple
  //  ntupler_->registerleaves(this);

  flows_ = iConfig.getParameter<std::vector<std::string> >("flows");
  
  //dummy output
  produces<double>();
}

ConfigurableAnalysis::~ConfigurableAnalysis()
{

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ConfigurableAnalysis::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // loop the requested selections
   for (Selections::iterator selection=selections_->begin(); selection!=selections_->end();++selection){
     //was this flow of filter actually asked for
     if (find(flows_.begin(), flows_.end(), selection->name())==flows_.end()) continue;

     //make a specific direction in the plotter
     plotter_->setDir(selection->name());

     // apply individual filters on the event
     std::map<std::string, bool> accept=selection->accept(iEvent);

     //loop the filters to make cumulative and allButOne job
     bool globalAccept=true;
     std::string separator="";
     std::string cumulative="";
     std::string allButOne="allBut_";
     for (Selection::iterator filterIt=selection->begin(); filterIt!=selection->end();++filterIt){
       Filter & filter=(**filterIt);
       bool lastCut=((filterIt+1)==selection->end());

       if (accept[filter.name()]){
	 //increment the directory name
	 cumulative+=separator+filter.name(); separator="_";
	 if ((selection->makeCumulativePlots() && !lastCut) || (selection->makeFinalPlots() && lastCut)){
	   plotter_->fill(cumulative,iEvent);
	 }
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
	 if (goodForAllButThisOne && selection->makeAllButOnePlots()){
	   plotter_->fill(allButOne+filter.name(),iEvent);
	 }
	 //do not global accept the event
	 globalAccept=false;
       }
     }
     if (globalAccept){
       //make the ntuple and put it in the event
       ntupler_->fill(selection->name(),iEvent);
     }
   }//loop the different filter order/number: loop the Selections

   
   //forget about this event.
   //   retriever_->clear();
   std::auto_ptr<double> dummy(new double(3));
   iEvent.put(dummy);
}
   

// ------------ method called once each job just before starting event loop  ------------
void 
ConfigurableAnalysis::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ConfigurableAnalysis::endJob() {
  //write out histograms
  plotter_->write();

  //print summary tables
  selections_->print();
}


DEFINE_FWK_MODULE(ConfigurableAnalysis);
