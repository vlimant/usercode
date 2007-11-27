// -*- C++ -*-
//
// Package:    BeamHaloMuonInTrackerFilter
// Class:      BeamHaloMuonInTrackerFilter
// 
/**\class BeamHaloMuonInTrackerFilter BeamHaloMuonInTrackerFilter.cc RecoTracker/BeamHaloMuonInTrackerFilter/src/BeamHaloMuonInTrackerFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Tue Apr  3 21:58:04 CDT 2007
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

#include <RecoMuon/L3MuonAnalyzer/interface/Analyzer.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

//
// class declaration
//

class BeamHaloMuonInTrackerFilter : public edm::EDFilter {
   public:
      explicit BeamHaloMuonInTrackerFilter(const edm::ParameterSet&);
      ~BeamHaloMuonInTrackerFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  std::string _category;
  IntrusiveAnalyzer  _analyzer;

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
BeamHaloMuonInTrackerFilter::BeamHaloMuonInTrackerFilter(const edm::ParameterSet& iConfig)
{
  _category = "BeamHaloMuonInTrackerFilter";
}


BeamHaloMuonInTrackerFilter::~BeamHaloMuonInTrackerFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
BeamHaloMuonInTrackerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //setup the analyzer
   _analyzer.setUp(iEvent,iSetup);

   //retreive sim tracks
   Handle<SimTrackContainer> simTracks;
   iEvent.getByType<SimTrackContainer>(simTracks);

   for (SimTrackContainer::const_iterator simIt = simTracks->begin(); simIt!=simTracks->end();++simIt){
     //only muons
     if (abs(simIt->type())!=13) continue;
     
     //number of tracker PSimHit in the tracker
     std::map< double ,TrackPSimHitRef > mc_hit=_analyzer.MapRefCount(simIt->trackId(),true,false,false,"HighTof");
     if (mc_hit.size()!=0) {edm::LogInfo(_category)<<mc_hit.size()<<" tracker HighTof PSimHits"; return true;}
     mc_hit=_analyzer.MapRefCount(simIt->trackId(),true,false,false,"LowTof");
     if (mc_hit.size()!=0) {edm::LogInfo(_category)<<mc_hit.size()<<" tracker LowTof PSimHits"; return true;}
   }

   //returns only true if the (anti-)muon has at least a simhit in the tracker (Hight/low time of flight)
   edm::LogInfo(_category)<<"event is skipped because no muon in the tracker";
   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
BeamHaloMuonInTrackerFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BeamHaloMuonInTrackerFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(BeamHaloMuonInTrackerFilter);
