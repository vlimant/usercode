#ifndef ANALYZER_H
#define ANALYZER_H


#include "FWCore/Framework/interface/Event.h"
//#include <Geometry/Vector/interface/LocalPoint.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>

#include <TrackingTools/TransientTrack/interface/TransientTrack.h>

#include <DataFormats/DetId/interface/DetId.h>

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include <Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h>
#include <Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h>
#include <RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h>

#include <RecoTracker/MeasurementDet/interface/MeasurementTracker.h>

#include <MagneticField/Engine/interface/MagneticField.h>

#include <TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h>

#include <TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h>
//#include <Geometry/Surface/interface/BoundPlane.h>
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include <TrackingTools/TransientTrack/interface/TransientTrackBuilder.h>
#include <TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h>

#include <TrackingTools/DetLayers/interface/NavigationSchool.h>

#include <SimDataFormats/Track/interface/SimTrackContainer.h>

#ifndef O
#define O(var) #var<<": "<<var<<"\n"
#endif

#ifndef Df
#define Df(var,os) os<<#var<<" "<<var<<"\n";
#endif
#ifndef Dfn
#define Dfn(var,os) os<<#var<<" "<<var<<" ";
#endif
#ifndef Dfl
#define Dfl(label,var,o) o<<label<<" "<<var<<"\n";
#endif
#ifndef Dfln
#define Dfln(label,var,o) o<<label<<" "<<var;
#endif


#ifndef D
#define D(var) std::cout<<#var<<" "<<var<<"\n";
#endif
#ifndef Dn
#define Dn(var) std::cout<<#var<<" "<<var<<" ";
#endif
#ifndef Dl
#define Dl(label,var) std::cout<<label<<" "<<var<<std::endl;
#endif
#ifndef Dln
#define Dln(label,var) std::cout<<label<<" "<<var;
#endif
#ifndef COMMENT
#define COMMENT(text) std::cout<<"### "<<text<<std::endl;
#endif
#ifndef COMMENTf
#define COMMENTf(text,o) o<<"### "<<text<<"\n";
#endif



using namespace std;

class IntrusiveAnalyzer{
 public:
  IntrusiveAnalyzer();
  ~IntrusiveAnalyzer();
  
  //will set every tools/data that are needed
  bool setUp(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void setDebug(bool d){_debug=d;};

  //recollect the proper topology of 'detid' and give the strip center (1D) of the given hit position (2D)
  LocalPoint stripcenter(LocalPoint striplocal, DetId detid);

  inline double delta_phi(double phi1,double phi2);
  inline double delta_R(double phi1,double eta1,double phi2,double eta2);
  
  //return the simulated track closest to the 'track_P' position when extrapolated to the given 'moduleId'
  //  const SimTrack * match(GlobalPoint track_P, const BoundPlane * surface,double dRmin=100000); //if (surface), then it is a position match on this surface
  SimTrackRef match(GlobalPoint track_P, const BoundPlane * surface,double dRmin=100000); //if (surface), then it is a position match on this surface
  //return the simulated track that match best the reconstructed track
  //first on position, and then on momentum if surface search fails
  //  const SimTrack * match(reco::TransientTrack & track,double dRmin=100000);
  SimTrackRef match(reco::TransientTrack & track,double dRmin=100000);


  //return a Trajectory state for a given simulated hit: no error matrix
  TrajectoryStateOnSurface convert(const PSimHit & simhit);//only from tracker simhits
  TrajectoryStateOnSurface convertFromMuonSystem(unsigned int ID);//same as above, from PSimHits in  the muon system
  //return a trajectory state of the simulated track: at vertex, no erro matrix
  FreeTrajectoryState convert(const SimTrack *);
  FreeTrajectoryState convert(SimTrackRef);
  FreeTrajectoryState convertFromSimTrack(uint trackId); //same as above, given the trackId
  
  //return the first psimhit on the 'detid', with the 'exact' 'trackID,' or not
  const PSimHit * simHitOnDetId(unsigned int trackID, DetId detid,bool exact = true);
  
  //count the number of hits on this track ID
  int Count(uint ID); 
  std::map< double, const PSimHit > MapCount(uint ID,bool silent=false,bool countStereo=false,bool processSelect=false,std::string suffix="LowTof",bool skipPixel=true);
  std::map< double ,TrackPSimHitRef > MapRefCount(uint ID,bool silent=true,bool countStereo=false,bool processSelect=false,std::string suffix="LowTof",bool skipPixel=true);
  std::map< double, DetId> DetCount(uint ID);
  
  //give the best simulated hit spatial match to the reconstructed hit
  //along with a code reflecting the process type, if there is a match, or not      
  enum RecHit_match_PSimHit { NOISE, OTHERTRACK, PRIMARY, SECONDARY };
  std::pair<RecHit_match_PSimHit, TrajectoryStateOnSurface> match(const TrackingRecHit & rechit,bool talk=false,unsigned int _trackID=0,bool Ltof=true);

  //give the best rechit spatial match of the simulated hit
  TransientTrackingRecHit::ConstRecHitPointer match(const PSimHit &  simhit); 
  
  inline SteppingHelixPropagator * prop(PropagationDirection d){ _prop->setPropagationDirection(d); return _prop;}
  
  std::string modulename(const DetId det);
  void show(const TrajectoryStateOnSurface & TSOS, char * label,ostream & o = cout);
  void show(const FreeTrajectoryState & FS,char * label,ostream & o = cout);
  void show(const DetId,bool showsim=false,ostream & o=std::cout);

  inline const TrackerGeometry * tracker(){return dynamic_cast<const TrackerGeometry*>( _measurementTracker->geomTracker());}
  inline const GeometricSearchTracker * geotracker(){return _measurementTracker->geometricSearchTracker();}
  inline edm::ESHandle<MeasurementTracker> & measurementTracker(){return _measurementTracker;}
  inline edm::ESHandle<MagneticField> & field(){return _field;}
  inline edm::ESHandle< GlobalTrackingGeometry > & globalTrackingGeometry(){ return _glbTrackinggeometry;}
  inline const TrajectoryStateTransform & transformer(){return _transformer;}
  inline const TransientTrackBuilder & TTbuilder() {return *_TTbuilder;}
  inline const TransientTrackingRecHitBuilder & TTRHbuilder() { return *_TTRHbuilder;}
  

  //  std::vector<const SimTrack*> particles(uint pdgid=13,int charge=0);
  std::vector<SimTrackRef> particles(uint pdgid=13,int charge=0);
  //const SimTrack * mother(uint trId);
  SimTrackRef mother(uint trId);
  //  TransientTrackingRecHit::ConstRecHitContainer roadContent(reco::TransientTrack & track , double chi2,bool withPxl=true, bool withStrip=true);
  
  typedef std::pair<TransientTrackingRecHit::ConstRecHitPointer, TrajectoryStateOnSurface> RecHitWithState;
  typedef std::vector<RecHitWithState> RecHitWithStateContainer;
  RecHitWithStateContainer roadContent(reco::TransientTrack & track , double chi2,bool withPxl=true, bool withStrip=true);

  std::vector<std::pair<SimTrackRef, reco::TrackRef> > match(std::vector<SimTrackRef> &, edm::Handle<reco::TrackCollection>&);
  //  std::vector<std::pair<SimTrackRef, reco::TrackRef> > match(std::vector<SimTrackRef> &, edm::Handle<reco::TrackCollection>&,const Plane & reference);
  //  std::vector<std::pair <int , int > > match(const SimTrackContainer, const reco::TrackCollection);

 private:

  const edm::Event * _theEvent;
  const edm::EventSetup * _theSetup;
 
  //  const TrackerGeometry * _tracker; 
  edm::ESHandle<MeasurementTracker>  _measurementTracker;
  edm::ESHandle< MagneticField> _field;
  edm::ESHandle<GlobalTrackingGeometry> _glbTrackinggeometry;
  TrajectoryStateTransform _transformer;

  TransientTrackBuilder * _TTbuilder;
  edm::ESHandle<TransientTrackingRecHitBuilder> _TTRHbuilder;

  SteppingHelixPropagator * _prop;
  bool _debug;

  NavigationSchool * _simpleNav;

  std::string _SimHitsSource;
  int _simIndexFix;

};


#endif
