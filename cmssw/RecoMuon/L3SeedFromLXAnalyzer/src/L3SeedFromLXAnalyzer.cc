// -*- C++ -*-
//
// Package:    L3SeedFromLXAnalyzer
// Class:      L3SeedFromLXAnalyzer
// 
/**\class L3SeedFromLXAnalyzer L3SeedFromLXAnalyzer.cc RecoMuon/L3SeedFromLXAnalyzer/src/L3SeedFromLXAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu Jan 24 23:48:18 CET 2008
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

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

//
// class decleration
//

class L3SeedFromLXAnalyzer : public edm::EDAnalyzer {
   public:
      explicit L3SeedFromLXAnalyzer(const edm::ParameterSet&);
      ~L3SeedFromLXAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag source;
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
L3SeedFromLXAnalyzer::L3SeedFromLXAnalyzer(const edm::ParameterSet& iConfig)

{
  source = iConfig.getParameter<edm::InputTag>("l3SeedLabel");
}


L3SeedFromLXAnalyzer::~L3SeedFromLXAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L3SeedFromLXAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const std::string category = "L3SeedFromLXAnalyzer";

  edm::ESHandle<MagneticField> fieldH;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldH);

  edm::ESHandle<TrackerGeometry> trackerGH;
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGH);

  edm::Handle<L3MuonTrajectorySeedCollection> seedH;
  iEvent.getByLabel(source , seedH);

  TrajectoryStateTransform transformer;

  uint is=0;
  uint isMax = seedH->size();
  for (; is!=isMax;++is){
    const L3MuonTrajectorySeed & seed = (*seedH)[is];
    const reco::TrackRef & l2 =  seed.l2Track();
    const l1extra::L1MuonParticleRef & l1 = seed.l1Particle();
    if (l1.isNull()) {edm::LogError(category)<<"l1 reference is not valid.";}
    else{
      edm::LogWarning(category)<<" l1 initial state is: "
			       <<"x: "<<l1->vertex()
			       <<" p: "<<l1->momentum();
    }
    if (l2.isNull()){edm::LogError(category)<<"l2 reference is not valid.";}
    else{
      FreeTrajectoryState l2state =  transformer.initialFreeState(*l2, fieldH.product());
      edm::LogWarning(category)<<" l2 initial state is: "
			       <<"x: "<<l2state.position()
			       <<" p: "<<l2state.momentum();
    }

    //make the seed state a persistent state
    TrajectoryStateOnSurface seedstate = transformer.transientState(seed.startingState(),
								    &trackerGH->idToDet(DetId(seed.startingState().detId()))->surface(),
								    fieldH.product());
    edm::LogWarning(category)<<" the seed state is: "
			     <<"x: "<<seedstate.globalPosition()
			     <<" p: "<<seedstate.globalMomentum()
			     <<" with: "<<seed.nHits()<<" recHits.";
    /*
      TrajectorySeed::const_iterator rhi = seed.range().first;
      TrajectorySeed::const_iterator rhiEnd = seed.range().second;
      uint nRH=0;
      for (;rhi!=rhiEnd;++rhiEnd){
      nRH++;}
    */

    
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
L3SeedFromLXAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L3SeedFromLXAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L3SeedFromLXAnalyzer);
