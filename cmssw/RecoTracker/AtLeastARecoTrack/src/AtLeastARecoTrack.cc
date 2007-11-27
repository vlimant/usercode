// -*- C++ -*-
//
// Package:    AtLeastARecoTrack
// Class:      AtLeastARecoTrack
// 
/**\class AtLeastARecoTrack AtLeastARecoTrack.cc RecoTracker/AtLeastARecoTrack/src/AtLeastARecoTrack.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Mon Apr 16 12:17:49 CDT 2007
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

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"

//
// class declaration
//

class AtLeastARecoTrack : public edm::EDFilter {
   public:
      explicit AtLeastARecoTrack(const edm::ParameterSet&);
      ~AtLeastARecoTrack();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag _recoTrackCollectionLabel;
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
AtLeastARecoTrack::AtLeastARecoTrack(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  _recoTrackCollectionLabel = iConfig.getParameter<edm::InputTag>("trackCollectionLabel");
}


AtLeastARecoTrack::~AtLeastARecoTrack()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
AtLeastARecoTrack::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   const uint nMintrack=1;
   const bool inclusive=false;

   //retreive track collection
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(_recoTrackCollectionLabel,tracks);

   //default is reject the event
   bool answer=false;
   
   if (tracks->size()>0){

     //loop over tracks
     for (reco::TrackCollection::const_iterator Tit = tracks->begin();Tit!=tracks->end();++Tit){
       //are both ends of the track at opposit sign z ?
       if (( Tit->innerPosition().z() * Tit->outerPosition().z()) < 0)
	 {answer = true; break;}
       else {continue;}
     }//loop over tracks
   }//any track at all

   /*     if (tracks->size()==nMintrack){ return true;}
	  if (inclusive){if (tracks->size()>nMintrack) return true;}
	  else return false;*/

   edm::LogInfo("AtLeastARecoTrack")<<"filter "<<((answer)?"accepts ":"rejects ")<<"the event";
   return answer;
}

// ------------ method called once each job just before starting event loop  ------------
void 
AtLeastARecoTrack::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AtLeastARecoTrack::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(AtLeastARecoTrack);
