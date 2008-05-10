// -*- C++ -*-
//
// Package:    TrajectoryOrTrackCandidateTeller
// Class:      TrajectoryOrTrackCandidateTeller
// 
/**\class TrajectoryOrTrackCandidateTeller TrajectoryOrTrackCandidateTeller.cc RecoTracker/TrajectoryOrTrackCandidateTeller/src/TrajectoryOrTrackCandidateTeller.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Tue Jan 29 20:05:09 CET 2008
// $Id: TrajectoryOrTrackCandidateTeller.cc,v 1.1 2008/02/15 02:53:24 vlimant Exp $
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


#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"


class TrajectoryOrTrackCandidateTeller : public edm::EDAnalyzer {
   public:
      explicit TrajectoryOrTrackCandidateTeller(const edm::ParameterSet&);
      ~TrajectoryOrTrackCandidateTeller();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  
  edm::InputTag inputTag;
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
TrajectoryOrTrackCandidateTeller::TrajectoryOrTrackCandidateTeller(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  inputTag = iConfig.getParameter<edm::InputTag>("input");
}


TrajectoryOrTrackCandidateTeller::~TrajectoryOrTrackCandidateTeller()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrajectoryOrTrackCandidateTeller::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   const std::string category = "TrajectoryOrTrackCandidateTeller";

   edm::Handle<TrackCandidateCollection> TCH;
   iEvent.getByLabel(inputTag, TCH);

   edm::ESHandle<MagneticField> fieldH;
   iSetup.get<IdealMagneticFieldRecord>().get(fieldH);

   edm::ESHandle<TrackerGeometry> trackerGH;
   iSetup.get<TrackerDigiGeometryRecord>().get(trackerGH);
   
   if (!TCH.failedToGet()){
     edm::LogVerbatim(category)<<"getting: "<<TCH->size()<<" TrackCandidates from: "<<inputTag;
     TrajectoryStateTransform transform;
     uint iTc=0;
     for (TrackCandidateCollection::const_iterator tcIt = TCH->begin();tcIt!= TCH->end(); ++ tcIt){
       uint nRh = tcIt->recHits().second-tcIt->recHits().first;
       TrajectoryStateOnSurface state = transform.transientState(tcIt->trajectoryStateOnDet(),
								 &trackerGH->idToDet(DetId(tcIt->trajectoryStateOnDet().detId()))->surface(),
								 fieldH.product());
       
       edm::LogVerbatim(category)<<iTc++<<"] on "<<tcIt->trajectoryStateOnDet().detId()<<" with: "<<nRh<<" recHits.\n" <<state;
     }
   }
   else{edm::LogWarning(category)<<"did not get TrackCandidate from: "<<inputTag;}

   edm::Handle<std::vector<Trajectory> > TH;
   iEvent.getByLabel(inputTag, TH);
   if (!TH.failedToGet()){
     edm::LogVerbatim(category)<<"getting: "<<TH->size()<<" Trajectory from: "<<inputTag;
     uint iT=0;
     for (std::vector<Trajectory>::const_iterator tIt=TH->begin();tIt!=TH->end();++tIt){
       edm::LogVerbatim(category)<<iT++<<"] with: "<<tIt->foundHits()<<" found hits.";
     }
     
   }
   else{edm::LogWarning(category)<<"did not get Trajectory from: "<<inputTag;} 

   if (TCH.failedToGet() && TH.failedToGet())
     {edm::LogError(category)<<"failed to get either TrackCandidate or Trajectory from: "<< inputTag; return;}
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrajectoryOrTrackCandidateTeller::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrajectoryOrTrackCandidateTeller::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrajectoryOrTrackCandidateTeller);
