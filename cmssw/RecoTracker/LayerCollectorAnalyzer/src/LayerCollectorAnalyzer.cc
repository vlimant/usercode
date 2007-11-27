// -*- C++ -*-
//
// Package:    LayerCollectorAnalyzer
// Class:      LayerCollectorAnalyzer
// 
/**\class LayerCollectorAnalyzer LayerCollectorAnalyzer.cc RecoTracker/LayerCollectorAnalyzer/src/LayerCollectorAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Tue Mar 20 19:54:05 CDT 2007
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

#include <RecoTracker/MeasurementDet/interface/MeasurementTracker.h>
#include <RecoTracker/Record/interface/CkfComponentsRecord.h>

#include "MagneticField/Engine/interface/MagneticField.h"
#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include <TrackingTools/Records/interface/TrackingComponentsRecord.h>

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "RecoTracker/TkNavigation/interface/StartingLayerFinder.h"
#include "RecoTracker/TkNavigation/interface/LayerCollector.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include <string>

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <FWCore/ParameterSet/interface/InputTag.h>

#include <RecoTracker/TkNavigation/interface/SimpleNavigationSchool.h>
#include <TrackingTools/DetLayers/interface/NavigationSetter.h>

//
// class decleration
//

class LayerCollectorAnalyzer : public edm::EDAnalyzer {
   public:
      explicit LayerCollectorAnalyzer(const edm::ParameterSet&);
      ~LayerCollectorAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag _STATrackCollectionLabel;
  bool _useCompatible;
  NavigationSchool * theNavigationSchool;
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
LayerCollectorAnalyzer::LayerCollectorAnalyzer(const edm::ParameterSet& iConfig)
  :  theNavigationSchool(NULL){
   //now do what ever initialization is needed
  _STATrackCollectionLabel = iConfig.getParameter<edm::InputTag>("src");
  _useCompatible = iConfig.getParameter<bool>("useCompatible");
}


LayerCollectorAnalyzer::~LayerCollectorAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
LayerCollectorAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get the measurementtracker
   edm::ESHandle<MeasurementTracker> _measurementTracker;
   iSetup.get<CkfComponentsRecord>().get(_measurementTracker);
     
   //get the magnetic field
   edm::ESHandle<MagneticField> _field;
   iSetup.get<IdealMagneticFieldRecord>().get(_field);
   
   //get the propagator
   edm::ESHandle<Propagator> _prop;
   //   iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlopng",_prop);
   iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong",_prop);

   //what are the tools we are testing
   Chi2MeasurementEstimator _estimator(9,3);
   StartingLayerFinder* _finder;
   StartingLayerFinder* _finderWithEstimator;
   /*
     if (_useCompatible)
     {     _finder =new StartingLayerFinder(_prop.product(),&_estimator,_measurementTracker.product()); }
     else
     {     _finder = new StartingLayerFinder(_prop.product()_measurementTracker.product());}*/
   _finder = new StartingLayerFinder(_prop.product(),NULL,_measurementTracker.product());
   _finderWithEstimator=new StartingLayerFinder(_prop.product(),&_estimator,_measurementTracker.product());

   //set a navigation school. otherwise ... no navigation ...
   if (!theNavigationSchool)
     theNavigationSchool   = new SimpleNavigationSchool(_measurementTracker->geometricSearchTracker(),_field.product());
   NavigationSetter setter(*theNavigationSchool);

   LayerCollector* _collector = new LayerCollector(_prop.product(),_finder,1,1);
   LayerCollector* _collectorWithEstimator = new LayerCollector(_prop.product(),_finderWithEstimator,1,1);
     
   //retreive a collection of Track
   Handle<reco::TrackCollection> _tracks;
   iEvent.getByLabel(_STATrackCollectionLabel,_tracks);

   //loop over then
   for (uint it=0;it!=_tracks->size();++it){
     const reco::Track & track = (*_tracks)[it];
       
     //get the state at impact
     TrajectoryStateTransform transform;
     FreeTrajectoryState cIPFTS=transform.initialFreeState(track,_field.product());
     
     if (cIPFTS.position().mag()==0) continue; //this is not very nice, people are aware of the problem.
     //https://hypernews.cern.ch/HyperNews/CMS/get/recoDevelopment/336.html
       
     //collect with old method
     std::vector<const DetLayer *> allOldLayers = _collector->allOldLayers(cIPFTS);

     //collect with new
     std::vector<const DetLayer *> allNewLayers = _collector->allLayers(cIPFTS);
     
     //collect with new and estimator
     std::vector<const DetLayer *> allNewLayersWithEstimator = _collectorWithEstimator->allLayers(cIPFTS);

     //crudely compare
     if (allOldLayers.size() != allNewLayers.size()){
       edm::LogError("LayerCollectorAnalyzer")<<"for state: \n"
					      <<cIPFTS
					      <<"(old) "<<allOldLayers.size()<<" != "<<allNewLayers.size()<<" (new)";}
     else{
       edm::LogInfo("LayerCollectorAnalyzer")<<"everything's fine at first order (old-new): ("<<allNewLayers.size()<<") layers collected";}

     if (allOldLayers.size() != allNewLayersWithEstimator.size()){
       edm::LogError("LayerCollectorAnalyzer")<<"for state: \n"
					      <<cIPFTS
					      <<"(old) "<<allOldLayers.size()<<" != "<<allNewLayersWithEstimator.size()<<" (est)";}
     else{
       edm::LogInfo("LayerCollectorAnalyzer")<<"everything's fine at first order (old-est): ("<<allNewLayersWithEstimator.size()<<") layers collected";}

     if (allNewLayers.size() != allNewLayersWithEstimator.size()){
       edm::LogError("LayerCollectorAnalyzer")<<"for state: \n"
					      <<cIPFTS
					      <<"(new) "<<allNewLayers.size()<<" != "<<allNewLayersWithEstimator.size()<<" (est)";}
     else{
       edm::LogInfo("LayerCollectorAnalyzer")<<"everything's fine at first order (new-est): ("<<allNewLayersWithEstimator.size()<<") layers collected";}
   }
   
   delete _finder;
   delete _collector;

}


// ------------ method called once each job just before starting event loop  ------------
void 
LayerCollectorAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LayerCollectorAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(LayerCollectorAnalyzer);
