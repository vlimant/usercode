#ifndef RecoTracker_TrackProducerFromSimTrack
#define  RecoTracker_TrackProducerFromSimTrack

// -*- C++ -*-
//
// Package:    TrackProducerFromSimTrack
// Class:      TrackProducerFromSimTrack
// 
/**\class TrackProducerFromSimTrack TrackProducerFromSimTrack.cc RecoTracker/TrackProducerFromSimTrack/src/TrackProducerFromSimTrack.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Fri Apr 13 17:42:20 CDT 2007
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

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackExtra.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <MagneticField/Engine/interface/MagneticField.h>

#include <TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h>
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>
#include <TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h>

#include <DataFormats/GeometrySurface/interface/Plane.h>
#include <TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToLocal.h>
#include <TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h>

#include "RecoTracker/ErrorMatrixAnalyzer/interface/ErrorMatrix.h"

#include <TrackingTools/GeomPropagators/interface/Propagator.h>
#include <TrackingTools/Records/interface/TrackingComponentsRecord.h>
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"


#include <DataFormats/GeometrySurface/interface/Plane.h>
#include <DataFormats/GeometrySurface/interface/Cylinder.h>

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>

#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>



#include "TMath.h"
#include "TRandom2.h"
#include "TString.h"


//
// class decleration
//

class TrackProducerFromSimTrack : public edm::EDProducer {
   public:
      explicit TrackProducerFromSimTrack(const edm::ParameterSet&);
      ~TrackProducerFromSimTrack();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
 
  bool selectSimTrack(const SimTrack & simtrack);
  reco::Track makeTrack(const SimTrack & simtrack);
  bool selectTrack(const reco::Track & recotrack);
  
  
  reco::TrackExtra * makeTrackExtra(const SimTrack &,
				    reco::Track &,
				    reco::TrackExtraCollection&);

  TH1 * get(const char * hname);
  bool attachRecHits(const SimTrack &,
		     reco::Track &,
		     reco::TrackExtra &,
		     TrackingRecHitCollection&);
  
  CurvilinearTrajectoryError smearingParameters(const FreeTrajectoryState & state);
  AlgebraicSymMatrix smearingParameters(const CurvilinearTrajectoryError & curv, const FreeTrajectoryState & state, const Plane &);
  AlgebraicVector smear(const AlgebraicVector & vector, const AlgebraicSymMatrix & matrix);
  int smear_charge(const FreeTrajectoryState &);

  
  // ----------member data ---------------------------
  
  std::string category;
  std::string instanceName;

  edm::Handle<edm::SimVertexContainer> simVtx;
  edm::Handle<edm::SimTrackContainer> simTracks;
  
  edm::ESHandle<MagneticField> field;
  edm::ESHandle<Propagator> propagatorAlong;

  //  std::string errorMatrixRootFile;
  edm::ParameterSet matrixProvider_pset;
  ErrorMatrix * matrixProvider;
  double errorMatrixOverEstimate;

  edm::RefProd<reco::TrackExtraCollection> refprodTE;
  edm::Ref<reco::TrackExtraCollection>::key_type TEi;

  edm::RefProd<TrackingRecHitCollection> refprodRH;
  edm::Ref<TrackingRecHitCollection>::key_type RHi;


  const edm::Event * ievent;
  const edm::EventSetup * isetup;

  TFile * plotFile;
  std::string plotFileName;

};

#endif
