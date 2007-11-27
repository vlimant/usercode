#include <RecoMuon/L3MuonAnalyzer/interface/Analyzer.h>


#include <RecoTracker/Record/interface/CkfComponentsRecord.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Geometry/CommonTopologies/interface/PixelTopology.h>
#include <Geometry/CommonTopologies/interface/StripTopology.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>
#include <Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h>


#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>

#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <DataFormats/MuonDetId/interface/DTWireId.h>

#include <DataFormats/SiPixelDetId/interface/PXBDetId.h>
#include <DataFormats/SiPixelDetId/interface/PXFDetId.h>
#include <DataFormats/SiStripDetId/interface/TIBDetId.h>
#include <DataFormats/SiStripDetId/interface/TOBDetId.h>
#include <DataFormats/SiStripDetId/interface/TIDDetId.h>
#include <DataFormats/SiStripDetId/interface/TECDetId.h>

#include <DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h>
#include <DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h>
#include <DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h>

//#include <RecoTracker/MeasurementDet/interface/MeasurementTracker.h>
#include <TrackingTools/MeasurementDet/interface/MeasurementDet.h>


#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "RecoTracker/TkNavigation/interface/StartingLayerFinder.h"
#include "RecoTracker/TkNavigation/interface/LayerCollector.h"

#include <TrackingTools/DetLayers/interface/NavigationSetter.h>
#include <RecoTracker/TkNavigation/interface/SimpleNavigationSchool.h>

IntrusiveAnalyzer::IntrusiveAnalyzer(): _TTbuilder(NULL),_prop(NULL),_simpleNav(NULL) 
{_debug=false;}

IntrusiveAnalyzer::~IntrusiveAnalyzer(){ 
  if (_simpleNav) delete _simpleNav;
  if (_prop) delete _prop;
}

bool IntrusiveAnalyzer::setUp(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //keep track of those to retreive collection and stuff
  _theEvent=&iEvent;_theSetup=&iSetup;
  
  iSetup.get<CkfComponentsRecord>().get( _measurementTracker);
  _measurementTracker->update(iEvent);

  //  _tracker = dynamic_cast<const TrackerGeometry*>( _measurementTracker->geomTracker());
  
  //get the global geometry
  iSetup.get<GlobalTrackingGeometryRecord>().get(_glbTrackinggeometry);

  //get the field
  iSetup.get<IdealMagneticFieldRecord>().get(_field);

  //get the propagator
  std::string _propagatorName = "SteppingHelixPropagatorAny";
  edm::ESHandle<Propagator> prop;
  iSetup.get<TrackingComponentsRecord>().get(_propagatorName,prop);

  if(_prop) delete _prop;
  _prop = dynamic_cast<SteppingHelixPropagator*>(prop->clone());

  if(!_TTbuilder)
    _TTbuilder = new  TransientTrackBuilder(field().product(),_glbTrackinggeometry);

  //setup the navigation in case it hasn't been set yet.
  if (!_simpleNav)
    _simpleNav= new SimpleNavigationSchool(geotracker(),_field.product());
  NavigationSetter setter(*_simpleNav);

  //get a tracking rechit builder
  const std::string builderName ="WithTrackAngle";
  iSetup.get<TransientRecHitRecord>().get(builderName,_TTRHbuilder);

  _SimHitsSource = "g4SimHits";
  _simIndexFix =0;
  
  return true;
}



LocalPoint IntrusiveAnalyzer::stripcenter(LocalPoint striplocal, DetId detid)
{
  //check if this is a glued thing
  const GeomDet *  gd = tracker()->idToDet(detid);
  if (gd->components().size()!=0) {/* this is a glued det: hence a 3d point */ return striplocal;}

  const GeomDetUnit * gdu = tracker()->idToDetUnit(detid); //this should not fail anymore

  if (!gdu){ edm::LogError("::stripcenter(...)")<<"no gdu for module in "<<modulename(detid)<<" ("<<detid.rawId()<<")";
    return striplocal;}

  const StripTopology * topology = dynamic_cast<const StripTopology *>(&gdu->topology());
  if (!topology)
    {
      if (dynamic_cast<const  PixelTopology*> (&gdu->topology())){
	return striplocal;//this is a pixel indeed...
      }
      else{
	edm::LogError("::stripcenter(...)")<<"topology does not cast ";
	COMMENT("::stripcenter(...) topology does not cast ");exit(1);}
    }

  float strip = topology->strip(striplocal);
  return topology->localPosition(strip);
}



inline double IntrusiveAnalyzer::delta_phi(double phi1,double phi2)
{  double dphi = phi1 -phi2;
  if (phi1<phi2) dphi = phi2-phi1;
  if (dphi>TMath::Pi()) 2*TMath::Pi() - dphi;
  return dphi; }
 
inline double IntrusiveAnalyzer::delta_R(double phi1,double eta1,double phi2,double eta2)
{  double deta=eta1-eta2;
  double dphi=delta_phi(phi1,phi2);
  return sqrt( deta*deta + dphi*dphi);}


//const SimTrack * IntrusiveAnalyzer::mother(uint trId){
SimTrackRef  IntrusiveAnalyzer::mother(uint trId){
  using  namespace edm;
  Handle<SimTrackContainer> simTracks;
  _theEvent->getByType<SimTrackContainer>(simTracks);
  Handle<SimVertexContainer> simVtx;
  _theEvent->getByType<SimVertexContainer>(simVtx);
  
  SimTrackContainer::const_iterator tracksCI = simTracks->begin();
  for(; tracksCI != simTracks->end(); tracksCI++) //loop over particles
    {
      if (trId!= tracksCI->trackId()) continue; //not the proper track
      if (tracksCI->vertIndex()==-1) continue; //no vertex for this guy

      const SimVertex & vtx= (*simVtx)[tracksCI->vertIndex()];

      if (vtx.parentIndex()!=-1){
	const SimTrack & mother_trk= (*simTracks)[vtx.parentIndex()];
	if (_debug) {Dn(trId);D(mother_trk.type());}
	//	return &mother_trk;
	return SimTrackRef(simTracks,vtx.parentIndex());
      }
    }
  if (_debug) {COMMENT("no mother to "<<trId);}
  //  return NULL;//by default
  return SimTrackRef();//by default
}


std::vector<std::pair<SimTrackRef,reco::TrackRef> > IntrusiveAnalyzer::match(std::vector<SimTrackRef> & sims, edm::Handle<reco::TrackCollection> & Tks){
  std::vector<std::pair<SimTrackRef,reco::TrackRef> > result;

  //loop over tracks and match them
  for (uint it=0;it!=Tks->size();++it){
    //build a transient track
    reco::TransientTrack TTk = TTbuilder().build(&(*Tks)[it]);
    
    SimTrackRef matchSim = match(TTk);
    
    if (matchSim.isNonnull()){ result.push_back(std::make_pair(matchSim,reco::TrackRef(Tks,it)));}
    else{result.push_back(std::make_pair(SimTrackRef(),reco::TrackRef(Tks,it)));}
  }
  //loop over the simtracks and get the non-matched ones
  
  for (uint is=0;is!=sims.size();++is){
    bool already=false;
    for (uint ip=0;ip!=result.size();++ip){
      if (result[ip].first==sims[is]) {already=true;break;}}

    if (!already){
      result.push_back(std::make_pair(sims[is],reco::TrackRef()));
    }
  }

  return result;
}



/*
std::vector<SimTrackRef> IntrusiveAnalyzer::particles(uint pdgid,int charge){
  using  namespace edm;
  Handle<SimTrackContainer> simTracks;
  _theEvent->getByType<SimTrackContainer>(simTracks);

  std::vector<SimTrackRef> result;
  for (uint iS=0;iS!=simTracks->size();++iS){
    const SimTrack * simT=(*simTracks)[iS];
    int trkPDG = simT->type();
    //select pdgid
    if (charge==0)
      {if( abs(trkPDG) != (int)pdgid) continue;}
    else
      {if (trkPDG!=(-1*charge)*((int)pdgid)) continue;}
    result.push_back(SimTrackRef(simTracks,iS));
  }
  return result;
}*/


  std::vector<SimTrackRef> IntrusiveAnalyzer::particles(uint pdgid,int charge)
  //std::vector<const SimTrack*> IntrusiveAnalyzer::particles(uint pdgid,int charge)
{
  using  namespace edm;
  Handle<SimTrackContainer> simTracks;
  _theEvent->getByType<SimTrackContainer>(simTracks);
  
  //  std::vector<const SimTrack*> result;
  std::vector<SimTrackRef> result;

  //  SimTrackContainer::const_iterator tracksCI = simTracks->begin();
  //  for(; tracksCI != simTracks->end(); tracksCI++){ //loop over particles
  for (uint is=0;is!=simTracks->size();++is)  {
    const SimTrack * tracksCI = &(*simTracks)[is];
      int trkPDG = tracksCI->type();
      //select pdgid
      if (charge==0)
	{if( abs(trkPDG) != (int)pdgid) continue;}
      else
	{if (trkPDG!=(-1*charge)*((int)pdgid)) continue;}
      //      result.push_back(&(*tracksCI)); 
      result.push_back(SimTrackRef(simTracks,is));
    }
  return result;
}

/*
  std::vector< std::pair <int , int > > IntrusiveAnalyzer::match(const SimTrackContainer, const reco::TrackCollection){
  std::map< uint , bool > simMatched;
  uint iSim=0;
  uint iTk=0;
  }*/



SimTrackRef IntrusiveAnalyzer::match(GlobalPoint track_P, const BoundPlane * surface,  double dRmin )
//const SimTrack *IntrusiveAnalyzer::match(GlobalPoint track_P, const BoundPlane * surface,  double dRmin )
{
  using namespace edm;

  //find a sim track according to the given track_P(osition) if surface is valid
  //otherwise do a momentum match

  GlobalPoint simtrack_P;
  SimTrackRef result;
  //  const SimTrack * result= NULL;  
  bool checkinside =(bool) (!surface);
  
  //generator level information
  //get the simulated particles list
  Handle<SimTrackContainer> simTracks;
  _theEvent->getByType<SimTrackContainer>(simTracks);
  if (!simTracks.isValid())/*abort*/
    {edm::LogError("IntrusiveAnalyzer::match()")<<"No Sim tracks found";return result;}

  double dRmin_init = dRmin;

  //loop over the particle list in the event
  //  SimTrackContainer::const_iterator tracksCI = simTracks->begin();
  //  for(; tracksCI != simTracks->end(); tracksCI++){ //loop over particles
  for (uint is=0;is!=simTracks->size();++is){
    const SimTrack * tracksCI = &(*simTracks)[is];
    int trkPDG = tracksCI->type();
    if (abs(trkPDG) != 13 )        continue; //I only do muons !

    FreeTrajectoryState state = convert(&(*tracksCI));
    
    if (checkinside)
      simtrack_P = GlobalPoint(tracksCI->momentum().px(),tracksCI->momentum().py(),tracksCI->momentum().pz());
    else
      {
	TrajectoryStateOnSurface state_at_mounsurface = prop(anyDirection)->propagate(state,*surface);
	if (!state_at_mounsurface.isValid())
	  { if (_debug){COMMENT("sim muon cannot be propagated to THIS muon system surface, take next one then !");} continue; }
	simtrack_P = state_at_mounsurface.globalPosition();
      }
    
    double deta = simtrack_P.eta() - track_P.eta();
    double dphi = delta_phi(simtrack_P.phi(),track_P.phi());
    double dR = sqrt(deta*deta + dphi*dphi);
    
    if(_debug){Dn(dR);D(simtrack_P);}
    
    if ( dR < dRmin ) //look for the closest dR match
      {
	dRmin = dR;
	//	  result = &(*tracksCI);
	result = SimTrackRef(simTracks,is);
      }
  }
  if (dRmin==dRmin_init) return SimTrackRef();
  else  return result;
}

//TransientTrackingRecHit::ConstRecHitContainer IntrusiveAnalyzer::roadContent(reco::TransientTrack & track , double chi2, bool withPxl,bool withStrip)
IntrusiveAnalyzer::RecHitWithStateContainer IntrusiveAnalyzer::roadContent(reco::TransientTrack & track , double chi2, bool withPxl,bool withStrip)
{
  //  TransientTrackingRecHit::ConstRecHitContainer result;
  RecHitWithStateContainer result;
  //create a chi2 estimator
  Chi2MeasurementEstimator estimator(chi2,sqrt(chi2));
  //create a starting layer finder
  StartingLayerFinder finder(prop(alongMomentum),&estimator,_measurementTracker.product());
  //  StartingLayerFinder finder(prop(alongMomentum),_measurementTracker.product());
  //create the layer collector
  LayerCollector collector(prop(alongMomentum),&finder,15,15); //10 cm around
  
  //  collect the layers with state
  std::vector<LayerCollector::LayerWithState> layersWithState = collector.allLayersWithState(track.initialFreeState());
  if(_debug) {D(layersWithState.size());}
  // find the compatible dets for each layer
  for (std::vector<LayerCollector::LayerWithState>::iterator lit=layersWithState.begin();lit!=layersWithState.end();++lit){
    if (!withPxl && (lit->first->subDetector()==GeomDetEnumerators::PixelBarrel  || lit->first->subDetector()==GeomDetEnumerators::PixelEndcap)) continue;
    if (!withStrip && (lit->first->subDetector()==GeomDetEnumerators::TIB ||
		       lit->first->subDetector()==GeomDetEnumerators::TOB ||
		       lit->first->subDetector()==GeomDetEnumerators::TID ||
		       lit->first->subDetector()==GeomDetEnumerators::TEC)) continue;
    std::vector<DetLayer::DetWithState> compatible =lit->first->compatibleDets(lit->second,*prop(anyDirection),estimator); 
    if (_debug){D(compatible.size());}
    for (std::vector<DetLayer::DetWithState>::iterator dit = compatible.begin();dit!=compatible.end();++dit){
      // collect the rechits for each  
      TransientTrackingRecHit::ConstRecHitContainer  thoseHits = _measurementTracker->idToDet(dit->first->geographicalId())->recHits(dit->second);
      if(_debug){D(thoseHits.size());}
      for (TransientTrackingRecHit::ConstRecHitContainer::const_iterator iTThit = thoseHits.begin();iTThit != thoseHits.end(); iTThit++){
	//estimate consistency
	MeasurementEstimator::HitReturnType est = estimator.estimate(dit->second,**iTThit);
	if (_debug) {Dn(est.first);D(est.second);}
	if (est.first){
	  //push them in the output
	  //	  result.push_back(*iTThit);
	  result.push_back(std::make_pair(*iTThit,dit->second));
	}
      }
    }
  }
  if (_debug){D(result.size());}
  return result;
}


SimTrackRef IntrusiveAnalyzer::match(reco::TransientTrack & track ,double dRmin)
//const SimTrack * IntrusiveAnalyzer::match(reco::TransientTrack & track ,double dRmin)
{
  //because of the poor resolution of the local muon reconstruction,  min dR match would fail
  //better to use DT,CSC, RPC local match ? Naaaa
  //better to swim the sim information to the outtermost measurement rather than swimming the reco track
  
  
  using namespace edm;
  //find a dR  match of this reco::track in the muon simulation information.
  //closest dR muon in simulation

  GlobalPoint track_P = track.outermostMeasurementState().globalPosition();

  const BoundPlane * surface = dynamic_cast<const BoundPlane* >(&track.outermostMeasurementState().surface());

  if (!surface)
    {edm::LogInfo(" IntrusiveAnalyzer::match(reco::TransientTrack...)")<<" outer most measurement surface does not cast to bound plane. check with momentum then";
      track_P = GlobalPoint(0,0,0)+track.trajectoryStateClosestToPoint(GlobalPoint(0,0,0)).momentum();}
  else
    {track_P = track.outermostMeasurementState().globalPosition();}

  return match(track_P,surface,dRmin);
}







TrajectoryStateOnSurface IntrusiveAnalyzer::convert(const PSimHit & simhit)
{
  const DetId did(simhit.detUnitId());
  if (!tracker()->idToDet(did)) { edm::LogError("IntrusiveAnalyzer::convert(const PSimHit & simhit)")<<did.rawId()<<" not recognized by the tracker geometry"; return TrajectoryStateOnSurface();}
  const Surface * surface = &tracker()->idToDet(did)->surface();
  if (!surface)
    {edm::LogError("IntrusiveAnalyzer::convert(const PSimHit & simhit)")<<" cannot find the surface of this simulated hit. most likely it is outside the tracker";return TrajectoryStateOnSurface(); }
  
  SurfaceSide surfaceside = atCenterOfSurface;
  GlobalPoint initialPoint=surface->toGlobal(simhit.localPosition());
  //  GlobalPoint initialPoint=surface->toGlobal(simhit.entryPoint());
  //  SurfaceSide surfaceside = beforeSurface;
  //  GlobalPoint initialPoint=surface->toGlobal(simhit.exitPoint());
  //  SurfaceSide surfaceside = afterSurface;

  GlobalVector initialMomentum=surface->toGlobal(simhit.momentumAtEntry());
  int initialCharge =  (simhit.particleType()>0) ? -1:1;
  CartesianTrajectoryError initialCartesianErrors(HepSymMatrix(6,0)); //no error at initial state
  const GlobalTrajectoryParameters initialParameters(initialPoint,initialMomentum,initialCharge,_field.product());
  return TrajectoryStateOnSurface(initialParameters,initialCartesianErrors,*surface,surfaceside);
}


TrajectoryStateOnSurface IntrusiveAnalyzer::convertFromMuonSystem(unsigned int trkInd)
{
  using namespace edm;
  //this will find the first sim track ID in the DT, CSC simhits to create a FTS seed

  PSimHitContainer::const_iterator muHits_CI;
  const GeomDetUnit * muHits_CI_layer=NULL;
  bool found=false;

  ESHandle<DTGeometry> dtGeomESH;
  Handle<PSimHitContainer> simHitsDT;

  ESHandle<CSCGeometry> cscGeomESH;
  Handle<PSimHitContainer> simHitsCSC;
  
  ESHandle<RPCGeometry> rpcGeomESH;
  Handle<PSimHitContainer> simHitsRPC;

  TrajectoryStateOnSurface empty;

  if(!found){
    //Muon DT geometry central
    _theSetup->get<MuonGeometryRecord>().get(dtGeomESH);
    
    //detector simulation level information
    //get list of hits in the DT
    _theEvent->getByLabel(_SimHitsSource, "MuonDTHits", simHitsDT);
    if (! simHitsDT.isValid() ){edm::LogError("IntrusiveAnalyzer::convertFromMuonSystem()")<<"No DT Hits found";return empty;}
    
    //loop over the DT hits to have an  initial measurement
    for ( muHits_CI = simHitsDT->begin(); muHits_CI != simHitsDT->end() ; muHits_CI++)       { 
      //      Dn(trkInd);D(muHits_CI->trackId());
      if (muHits_CI->trackId() != trkInd+_simIndexFix ) continue; 
      DTWireId wId(muHits_CI->detUnitId());
      muHits_CI_layer = dtGeomESH->layer(wId);
      if (muHits_CI_layer == 0){edm::LogError("IntrusiveAnalyzer::convertFromMuonSystem()")<<"Failed to get DT detector unit";continue;}
      found=true; break;
    }
  }

  if(!found){
    //Muon CSC geometry central
    _theSetup->get<MuonGeometryRecord>().get(cscGeomESH);
    
    //get list of hits in the CSC
    _theEvent->getByLabel(_SimHitsSource, "MuonCSCHits", simHitsCSC);
    if (! simHitsCSC.isValid() ){edm::LogError("IntrusiveAnalyzer::convertFromMuonSystem()")<<"No CSC Hits found";return empty;}

    //loop over the CSC hits to have an  initial measurement
    for ( muHits_CI = simHitsCSC->begin(); muHits_CI != simHitsCSC->end() ; muHits_CI++)       { 
      //      Dn(trkInd);D(muHits_CI->trackId());
      if (muHits_CI->trackId() != trkInd+_simIndexFix ) continue;
      CSCDetId cscId(muHits_CI->detUnitId());
      muHits_CI_layer = cscGeomESH->layer(cscId);
      if (muHits_CI_layer == 0){edm::LogError("IntrusiveAnalyzer::convertFromMuonSystem()")<<"Failed to get CSC detector unit";continue;}
      found=true; break;
    }
  }

  if(!found){
    //Muon RPC geometry central
    _theSetup->get<MuonGeometryRecord>().get(rpcGeomESH);
    
    //get list of hits in the RPC
    _theEvent->getByLabel(_SimHitsSource, "MuonRPCHits", simHitsRPC);
    if (! simHitsRPC.isValid() ){edm::LogError("IntrusiveAnalyzer::convertFromMuonSystem()")<<"No RPC Hits found";return empty;}
    
    //loop over the RPC hits to have an  initial measurement
    for (muHits_CI = simHitsRPC->begin(); muHits_CI != simHitsRPC->end() ; muHits_CI++)       { 
      //      Dn(trkInd);D(muHits_CI->trackId());
      if (muHits_CI->trackId() != trkInd+_simIndexFix ) continue;
      RPCDetId rpcId(muHits_CI->detUnitId());
      muHits_CI_layer = rpcGeomESH->roll(rpcId);
      if (muHits_CI_layer == 0){edm::LogError("IntrusiveAnalyzer::convertFromMuonSystem()")<<"Failed to get RPC detector unit";continue;}
      found=true; break;
    }
  }
  
  if (found)
    {
      if (_debug){      
	D(muHits_CI->localPosition());
	D(muHits_CI->momentumAtEntry());
	Dn(trkInd);D(muHits_CI->trackId());}

      const Surface& surf = muHits_CI_layer->surface();
      const GlobalPoint initialPoint(surf.toGlobal(muHits_CI->localPosition()));
      const GlobalVector initialMomentum(surf.toGlobal(muHits_CI->momentumAtEntry()));
      int initialCharge = (muHits_CI->particleType()>0) ? -1 : 1;
      const GlobalTrajectoryParameters initialParameters(initialPoint,initialMomentum,initialCharge,_field.product());
      //FIXME. could have a parametrized error matrix
      const CartesianTrajectoryError initialCartesianErrors(HepSymMatrix(6,0)); //no error at initial state
      return TrajectoryStateOnSurface(initialParameters,initialCartesianErrors,surf);
    }

  if (_debug) { COMMENT("no DT,RPC,CSC hit for track index:"<<trkInd);}
  return empty;
}
 
FreeTrajectoryState IntrusiveAnalyzer::convert(SimTrackRef simRef){
  if (simRef.isNonnull()){
    return convert(&(*simRef));}
  else { edm::LogError("IntrusiveAnalyzer::convert")<<"converting a null ref to simtrack into invalid FTS";return FreeTrajectoryState();}
};

FreeTrajectoryState IntrusiveAnalyzer::convert(const SimTrack *tracksCI)
{
  using namespace edm;
  GlobalPoint initialPoint(0,0,0);//na, this cannot be right !

  int vtx_i = tracksCI->vertIndex();
  if (vtx_i!=-1){
    Handle<SimVertexContainer> simVtx;
    _theEvent->getByType<SimVertexContainer>(simVtx);
    //get the proper vertex
    if ((uint)vtx_i<simVtx->size())
      {
	const HepLorentzVector& pos = (*simVtx)[vtx_i].position();
	initialPoint = GlobalPoint(pos.x(),pos.y(),pos.z());}}
    

  GlobalVector initialMomentum(tracksCI->momentum().px(),tracksCI->momentum().py(),tracksCI->momentum().pz());
  
  int initialCharge = (tracksCI->type()>0) ? -1 : 1;
  CartesianTrajectoryError initialCartesianErrors(HepSymMatrix(6,0)); //no error at initial state
  //FIXME. could have a parametrized error matrix
  const GlobalTrajectoryParameters initialParameters(initialPoint,initialMomentum,initialCharge,_field.product());
  return FreeTrajectoryState(initialParameters,initialCartesianErrors);
}

FreeTrajectoryState IntrusiveAnalyzer::convertFromSimTrack(uint trackId)
{
  using namespace edm;
  //convert the sim track with trackId to a free trajectory state on surface
  //generator level information
  //get the simulated particles list
  Handle<SimTrackContainer> simTracks;
  _theEvent->getByType<SimTrackContainer>(simTracks);
  if (! simTracks.isValid() )
    {edm::LogError("IntrusiveAnalyzer::convertFromSimTrack(uint trackId)")<<"No Sim tracks found";return FreeTrajectoryState();}
   
  //loop over the particle list in the event
  SimTrackContainer::const_iterator tracksCI = simTracks->begin();
  for(; tracksCI != simTracks->end(); tracksCI++) {if (tracksCI->trackId() ==trackId) return convert(&(*tracksCI)); }
  return FreeTrajectoryState();
}

const PSimHit * IntrusiveAnalyzer::simHitOnDetId(unsigned int ID, DetId det,bool exactId)
{
  //the ID is coming from SimTrack::trackId()
  using namespace edm;
  switch (det.subdetId())
    {
    case 1:{
      PXBDetId detid(det);
      //get list of sim hits in the PXB
      Handle<PSimHitContainer> simHitsPXB;
      _theEvent->getByLabel(_SimHitsSource, "TrackerHitsPixelBarrelLowTof", simHitsPXB);
      if (! simHitsPXB.isValid() ){edm::LogError("IntrusiveAnalyzer::simHitOnDetId(...)")<< "No PXB sim hits found";}
      else{
	for ( PSimHitContainer::const_iterator Hits_CI = simHitsPXB->begin();Hits_CI != simHitsPXB->end() ; Hits_CI++) 
	  {
	    PXBDetId id(Hits_CI->detUnitId());
	    //	      Dn(id.rawId());Dn(id.layer());D(Hits_CI->trackId());
	    if (id.layer() != detid.layer()) continue;
	    if ((id.rawId() == det.rawId()) && ((Hits_CI->trackId() == ID+ _simIndexFix) || !exactId ) && Hits_CI->localPosition().z()==0 ) return &(*Hits_CI);
	  }
      }
      break;}
    case 2:{
      PXFDetId detid(det);
      //get list of sim hits in the PXF
      Handle<PSimHitContainer> simHitsPXF;
      _theEvent->getByLabel(_SimHitsSource, "TrackerHitsPixelEndcapLowTof", simHitsPXF);
      if (! simHitsPXF.isValid() ){edm::LogError("IntrusiveAnalyzer::simHitOnDetId(...)")<< "No PXF sim hits found";}
      else{
	for ( PSimHitContainer::const_iterator Hits_CI = simHitsPXF->begin();Hits_CI != simHitsPXF->end() ; Hits_CI++) 
	  {
	    PXFDetId id(Hits_CI->detUnitId());
	    //	      Dn(id.rawId());Dn(id.layer());D(Hits_CI->trackId());
	    if (id.disk() != detid.disk()) continue;
	    if ((id.rawId() == det.rawId()) && ((Hits_CI->trackId() == ID+ _simIndexFix) || !exactId ) && Hits_CI->localPosition().z()==0 ) return &(*Hits_CI);
	  }
      }
      break;}
    case 3:{
      TIBDetId detid(det);
      //get list of sim hits in the TIB
      Handle<PSimHitContainer> simHitsTIB;
      _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTIBLowTof", simHitsTIB);
      if (! simHitsTIB.isValid() ){edm::LogError("IntrusiveAnalyzer::simHitOnDetId(...)")<< "No TIB sim hits found";}
      else{
	for ( PSimHitContainer::const_iterator Hits_CI = simHitsTIB->begin();Hits_CI != simHitsTIB->end() ; Hits_CI++) 
	  {
	    TIBDetId id(Hits_CI->detUnitId());
	    //	      Dn(id.rawId());Dn(id.layer());D(Hits_CI->trackId());
	    if (id.layer() != detid.layer()) continue;
	    if ((id.rawId() == det.rawId()) && ((Hits_CI->trackId() == ID+ _simIndexFix) || !exactId ) && Hits_CI->localPosition().z()==0 ) return &(*Hits_CI);
	  }
      }
      break;}
    case 4:{
      TIDDetId detid(det);
      //get list of sim hits in the TID
      Handle<PSimHitContainer> simHitsTID;
      _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTIDLowTof", simHitsTID);
      if (! simHitsTID.isValid() ){edm::LogError("IntrusiveAnalyzer::simHitOnDetId(...)")<< "No TID sim hits found";}
      else{
	for ( PSimHitContainer::const_iterator Hits_CI = simHitsTID->begin();
	      Hits_CI != simHitsTID->end() ; Hits_CI++) 
	  {
	    TIDDetId id(Hits_CI->detUnitId());
	    //	      Dn(id.subdetId());Dn(id.wheel());D(Hits_CI->trackId())
	    if (id.wheel() != detid.wheel()) continue;
	    if ((id.rawId() == det.rawId()) && ((Hits_CI->trackId() == ID+ _simIndexFix) || !exactId ) && Hits_CI->localPosition().z()==0 ) return &(*Hits_CI);
	    
	  }
      }
      break;}
    case 5:{     
      TOBDetId detid(det);
      //get list of sim hits in the TOB
      Handle<PSimHitContainer> simHitsTOB;
      _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTOBLowTof", simHitsTOB);
      if (! simHitsTOB.isValid() ){edm::LogError("IntrusiveAnalyzer::simHitOnDetId(...)")<< "No TOB sim hits found";}
      else{
	for ( PSimHitContainer::const_iterator Hits_CI = simHitsTOB->begin();
	      Hits_CI != simHitsTOB->end() ; Hits_CI++) 
	  {
	    TOBDetId id(Hits_CI->detUnitId());
	    //	      Dn(id.rawId());Dn(id.layer());D(Hits_CI->trackId());
	    if (id.layer() != detid.layer()) continue;
	    if ((id.rawId() == det.rawId()) && ((Hits_CI->trackId() == ID+ _simIndexFix) || !exactId ) && Hits_CI->localPosition().z()==0 ) return &(*Hits_CI);
	  }
      }
      break;}
    case 6:{
      TECDetId detid(det);
      //get list of sim hits in the TEC
      Handle<PSimHitContainer> simHitsTEC;
      _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTECLowTof", simHitsTEC);
      if (! simHitsTEC.isValid() ){edm::LogError("IntrusiveAnalyzer::simHitOnDetId(...)")<< "No TEC sim hits found";}
      else{
	for ( PSimHitContainer::const_iterator Hits_CI = simHitsTEC->begin();
	      Hits_CI != simHitsTEC->end() ; Hits_CI++) 
	  {
	    TECDetId id(Hits_CI->detUnitId());
	    //	      Dn(id.subdetId());Dn(id.wheel());D(Hits_CI->trackId())
	    if (id.wheel() != detid.wheel()) continue;
	    if ((id.rawId() == det.rawId()) && ((Hits_CI->trackId() == ID+ _simIndexFix) || !exactId ) && Hits_CI->localPosition().z()==0 ) return &(*Hits_CI);
	  }
      }
      break;}
    }
  return NULL;  //by default
}

int IntrusiveAnalyzer:: Count(uint ID)
{
  //will display the pair < distance , module ID of the "primary" PSimHits of ID
  using namespace std;
  map < double, const PSimHit > simmap = MapCount(ID,true);
  if (_debug){
    for ( map < double , const PSimHit >::const_iterator it =  simmap.begin();it!= simmap.end();it++)	
      {Dn(it->first); D(it->second.detUnitId());D(match(it->second));}
  }
  return simmap.size();
}


std::map< double, const PSimHit > IntrusiveAnalyzer::MapCount(uint ID ,bool silent,bool countStereo,bool processSelect, std::string suffix,bool skipPixel)
{
  edm::LogError("IntrusiveAnalyzer::MapCount")<<"you should not do that ! but that's alright";
  std::map< double, const PSimHit > result;
  std::map< double ,TrackPSimHitRef > temp = MapRefCount(ID,silent,countStereo,processSelect,suffix,skipPixel);
  for (std::map< double ,TrackPSimHitRef >::iterator it = temp.begin();it!=temp.end();++it){
    result.insert(make_pair(it->first,*it->second));  }
  return result;}

std::map< double ,TrackPSimHitRef > IntrusiveAnalyzer::MapRefCount(uint ID ,bool silent,bool countStereo,bool processSelect,std::string suffix,bool skipPixel)
{
  //ID is coming from SimTrack::trackId()
  using namespace edm;
  using namespace std;
  
  if(_debug){  Dn(ID);D(ID+_simIndexFix);}
  
  //  map< double, const PSimHit > hitlist;

  map< double ,TrackPSimHitRef > hitlist;

  const double magPlus=0.000001;
  double mag =0;
  map<uint32_t, bool> done;

  //find the number of simulated hits on the track ID
  int NhitsOnTrack=0;

  /* more generci way ? not straight forward because of stereo()*/
  string labels[6]={"TrackerHitsPixelBarrel",
		    "TrackerHitsPixelEndcap",
		    "TrackerHitsTIB",
		    "TrackerHitsTOB",
		    "TrackerHitsTID",
		    "TrackerHitsTEC"};

  Handle<PSimHitContainer> simHitsXXX;
  for (uint il=0;il!=6;++il){

    if (skipPixel){ 
      string::size_type loc =labels[il].find("Pixel");
      if (loc != string::npos ) {
	edm::LogInfo("IntrusiveAnalyzer::MapCount(...)")<<"skipping: "<<labels[il]<<" for PSimHit search";}
      }
    
    _theEvent->getByLabel(_SimHitsSource, labels[il]+suffix,simHitsXXX);
    if (! simHitsXXX.isValid() ){edm::LogError("IntrusiveAnalyzer::MapCount(...)")<< "no "<<labels[il]<<" sim hits found"; continue;}

    for (uint ips=0;ips!=simHitsXXX->size();++ips){
      const PSimHit * psh= &(*simHitsXXX)[ips];
      if (psh->trackId() == ID + _simIndexFix) {
	if (done.find(psh->detUnitId()) != done.end()) continue;
	if( processSelect && psh->processType() !=2) continue;
	DetId id(psh->detUnitId());
	bool jump=false;
	if (!countStereo){
	  switch(id.subdetId()){
	  case 1:
	  case 2:
	    jump=false;
	    break;
	  case 3:{ TIBDetId sid(id); if (sid.stereo()) jump=true; break;}
	  case 4:{ TIDDetId sid(id); if (sid.stereo()) jump=true; break;}
	  case 5:{ TOBDetId sid(id); if (sid.stereo()) jump=true; break;}
	  case 6:{ TECDetId sid(id); if (sid.stereo()) jump=true; break;}
	  }}
	if (jump) {/*std::cout<<"stereo det id "<<psh->detUnitId()<<"jumped"<<std::endl;*/continue;}
	NhitsOnTrack++;
	done[id.rawId()]=true;
	mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsXXX,ips)));
      }
    }
  }

  return hitlist;

  //get list of sim hits in the PXB
  Handle<PSimHitContainer> simHitsPXB;
  _theEvent->getByLabel(_SimHitsSource, "TrackerHitsPixelBarrelLowTof", simHitsPXB);
  if (! simHitsPXB.isValid() ){edm::LogError("IntrusiveAnalyzer::MapCount(...)")<< "No PXB sim hits found";  }
  else{
    //    for ( PSimHitContainer::const_iterator Hits_CI = simHitsPXB->begin(); Hits_CI != simHitsPXB->end() ; Hits_CI++) {
    for (uint ips = 0;ips!=simHitsPXB->size(); ++ips){ const PSimHit * Hits_CI=&(*simHitsPXB)[ips];
      if (Hits_CI->trackId() == ID + _simIndexFix)
	{
	    PXBDetId id(Hits_CI->detUnitId());
	    if (!silent && _debug){
	      Dn(id.rawId());Dn(id.subdetId());Dn(id.layer());Dn(id.module());D(Hits_CI->localPosition());
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(1,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,1,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,0,1)));
	      Dn(Hits_CI->particleType());Dn(Hits_CI->processType());D(Hits_CI->energyLoss());
	    }
	    //	    if (Hits_CI->localPosition().z()!=0) continue; //z!=0 seems to be a bug ...
	    if (done.find(id.rawId()) != done.end()) continue; //one hit per module. not more
	    //	    if ( !countStereo && id.stereo()) continue;
	    if( processSelect && Hits_CI->processType() !=2) continue;
	    NhitsOnTrack++;
	    done[id.rawId()]=true;
	    mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	    while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	    //	    hitlist.insert(make_pair(mag,*Hits_CI));
	    hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsPXB,ips)));
	  } 
      }
  }


  //get list of sim hits in the PXF
  Handle<PSimHitContainer> simHitsPXF;
  _theEvent->getByLabel(_SimHitsSource, "TrackerHitsPixelEndcapLowTof", simHitsPXF);
  if (! simHitsPXF.isValid() ){edm::LogError("IntrusiveAnalyzer::MapCount(...)")<< "No PXF sim hits found";  }
  else{
    //    for ( PSimHitContainer::const_iterator Hits_CI = simHitsPXF->begin(); Hits_CI != simHitsPXF->end() ; Hits_CI++) {
    for (uint ips = 0;ips!=simHitsPXF->size(); ++ips){ const PSimHit * Hits_CI=&(*simHitsPXF)[ips];
	if (Hits_CI->trackId() == ID + _simIndexFix)
	  {
	    PXFDetId id(Hits_CI->detUnitId());
	    if (!silent && _debug){
	      Dn(id.rawId());Dn(id.subdetId());Dn(id.disk());Dn(id.module());D(Hits_CI->localPosition());
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(1,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,1,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,0,1)));
	      Dn(Hits_CI->particleType());Dn(Hits_CI->processType());D(Hits_CI->energyLoss());
	    }
	    //	    if (Hits_CI->localPosition().z()!=0) continue; //z!=0 seems to be a bug ...
	    if (done.find(id.rawId()) != done.end()) continue; //one hit per module. not more
	    //	    if ( !countStereo && id.stereo()) continue;
	    if( processSelect && Hits_CI->processType() !=2) continue;
	    NhitsOnTrack++;
	    done[id.rawId()]=true;
	    mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	    while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	    //	    hitlist.insert(make_pair(mag,*Hits_CI));
	    hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsPXF,ips)));
	  } 
      }
  }

  //get list of sim hits in the TIB
  Handle<PSimHitContainer> simHitsTIB;
  _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTIBLowTof", simHitsTIB);
  if (! simHitsTIB.isValid() ){edm::LogError("IntrusiveAnalyzer::MapCount(...)")<< "No TIB sim hits found";  }
  else{
    //    for ( PSimHitContainer::const_iterator Hits_CI = simHitsTIB->begin(); Hits_CI != simHitsTIB->end() ; Hits_CI++)       {
    for (uint ips = 0;ips!=simHitsTIB->size(); ++ips){ const PSimHit * Hits_CI=&(*simHitsTIB)[ips]; 
	if (Hits_CI->trackId() == ID + _simIndexFix)
	  {
	    TIBDetId id(Hits_CI->detUnitId());
	    if (!silent && _debug){
	      Dn(id.rawId());Dn(id.subdetId());Dn(id.layer());Dn(id.module());Dn(id.stereo());Dn(id.glued());Dn(id.partnerDetId());D(Hits_CI->localPosition());
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(1,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,1,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,0,1)));
	      Dn(Hits_CI->particleType());Dn(Hits_CI->processType());D(Hits_CI->energyLoss());
	    }
	    //	    if (Hits_CI->localPosition().z()!=0) continue; //z!=0 seems to be a bug ...
	    if (done.find(id.rawId()) != done.end()) continue; //one hit per module. not more
	    if ( !countStereo && id.stereo()) continue;
	    if( processSelect && Hits_CI->processType() !=2) continue;
	    NhitsOnTrack++;
	    done[id.rawId()]=true;
	    mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	    while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	    //	    hitlist.insert(make_pair(mag,*Hits_CI));
	    hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsTIB,ips))); 
	  } 
      }
  }
  

  //get list of sim hits in the TID
  Handle<PSimHitContainer> simHitsTID;
  _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTIDLowTof", simHitsTID);
  if (! simHitsTID.isValid() ){edm::LogError("IntrusiveAnalyzer:::MapCount(...)")<< "No TID sim hits found"; }
  else{
    //    for ( PSimHitContainer::const_iterator Hits_CI = simHitsTID->begin();Hits_CI != simHitsTID->end() ; Hits_CI++)       {
    for (uint ips = 0;ips!=simHitsTID->size(); ++ips){ const PSimHit * Hits_CI=&(*simHitsTID)[ips];  
	if  (Hits_CI->trackId() == ID+ _simIndexFix) 
	  {
	    TIDDetId id(Hits_CI->detUnitId());
	    if (!silent &&_debug){
	      Dn(id.rawId());Dn(id.subdetId());Dn(id.wheel());Dn(id.ring());Dn(id.stereo());Dn(id.glued());Dn(id.partnerDetId());D(Hits_CI->localPosition());
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(1,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,1,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,0,1)));
	      Dn(Hits_CI->particleType());Dn(Hits_CI->processType());D(Hits_CI->energyLoss());
	    }
	    //	    if (Hits_CI->localPosition().z()!=0) continue; //z!=0 seems to be a bug ...
	    if (done.find(id.rawId()) != done.end()) continue; //one hit per module. not more
	    if ( !countStereo && id.stereo()) continue;
	    if( processSelect && Hits_CI->processType() !=2) continue;
	    NhitsOnTrack++;
	    done[id.rawId()]=true;
	    mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	    while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	    //	    hitlist.insert(make_pair(mag,*Hits_CI));
	    hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsTID,ips))); 
	  }
      }
  }
  
  //get list of sim hits in the TOB
  Handle<PSimHitContainer> simHitsTOB;
  _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTOBLowTof", simHitsTOB);
  if (! simHitsTOB.isValid() ){edm::LogError("IntrusiveAnalyzer:::MapCount(...)")<< "No TOB sim hits found"; }
  else{
    //    for ( PSimHitContainer::const_iterator Hits_CI = simHitsTOB->begin(); Hits_CI != simHitsTOB->end() ; Hits_CI++)      {
    for (uint ips = 0;ips!=simHitsTOB->size(); ++ips){ const PSimHit * Hits_CI=&(*simHitsTOB)[ips];  
	if (Hits_CI->trackId() == ID+ _simIndexFix) 
	  {
	    TOBDetId id(Hits_CI->detUnitId());
	    if (!silent && _debug){
	      Dn(id.rawId());Dn(id.subdetId());Dn(id.layer());Dn(id.module());Dn(id.stereo());Dn(id.glued());Dn(id.partnerDetId());D(Hits_CI->localPosition());
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(1,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,1,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,0,1)));
	      Dn(Hits_CI->particleType());Dn(Hits_CI->processType());D(Hits_CI->energyLoss());
	    }
	    //	    if (Hits_CI->localPosition().z()!=0) continue; //z!=0 seems to be a bug ...
	    if (done.find(id.rawId()) != done.end()) continue; //one hit per module. not more
	    if (!countStereo && id.stereo()) continue;
	    if( processSelect && Hits_CI->processType() !=2) continue;
	    NhitsOnTrack++;
	    done[id.rawId()]=true;
	    mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	    while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	    //	    hitlist.insert(make_pair(mag,*Hits_CI));
	    hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsTOB,ips))); 

	  }
      }
  }
  
  //get list of sim hits in the TEC
  Handle<PSimHitContainer> simHitsTEC;
  _theEvent->getByLabel(_SimHitsSource, "TrackerHitsTECLowTof", simHitsTEC);
  if (! simHitsTEC.isValid() ){edm::LogError("IntrusiveAnalyzer::find")<< "No TEC sim hits found";}
  else{
    //    for ( PSimHitContainer::const_iterator Hits_CI = simHitsTEC->begin();Hits_CI != simHitsTEC->end() ; Hits_CI++)       {
    for (uint ips = 0;ips!=simHitsTEC->size(); ++ips){ const PSimHit * Hits_CI=&(*simHitsTEC)[ips];  
	if (Hits_CI->trackId() == ID+ _simIndexFix)
	  { 
	    TECDetId id(Hits_CI->detUnitId());
	    if (!silent && _debug){
	      Dn(id.rawId());Dn(id.subdetId());Dn(id.wheel());Dn(id.ring());Dn(id.stereo());Dn(id.glued());Dn(id.partnerDetId());D(Hits_CI->localPosition());
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(1,0,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,1,0)));
	      D(tracker()->idToDet(id)->surface().toGlobal(LocalVector(0,0,1)));
	      Dn(Hits_CI->particleType());Dn(Hits_CI->processType());D(Hits_CI->energyLoss());
	    }
	    //	    if (Hits_CI->localPosition().z()!=0) continue; //z!=0 seems to be a bug ...
	    if (done.find(id.rawId()) != done.end()) continue; //one hit per module. not more
	    if (!countStereo &&id.stereo()) continue;
	    if( processSelect && Hits_CI->processType() !=2) continue;
	    NhitsOnTrack++;
	    done[id.rawId()]=true;
	    mag = tracker()->idToDet(id)->surface().toGlobal(LocalPoint(0,0,0)).mag();
	    while (hitlist.find(mag)!=hitlist.end())     mag+=magPlus;
	    //	    hitlist.insert(make_pair(mag,*Hits_CI));
	    hitlist.insert(make_pair(mag,TrackPSimHitRef(simHitsTEC,ips)));     
	  }
      }
  }
  return hitlist;
}


std::map< double, DetId> IntrusiveAnalyzer::DetCount(uint ID)
{
  using namespace std;
  using namespace edm;

  map < double, DetId> detlist;
  map < double, const PSimHit > simmap = MapCount(ID);
  for ( map < double , const PSimHit >::const_iterator it =  simmap.begin();it!= simmap.end();it++)
    { detlist[it->first] = DetId(it->second.detUnitId()); }
  return detlist;
}


std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> IntrusiveAnalyzer:: match(const  TrackingRecHit & rechit,bool talk,unsigned int  _trackID,bool Ltof)
{
  using namespace std;
  using namespace edm;

  const GeomDet * geodet = tracker()->idToDet(rechit.geographicalId());  
  if (!geodet){edm::LogError("RoadSearchMuonSeedeFinder::match(...)")
      <<"no detlayer for id "<<rechit.geographicalId().rawId(); return make_pair(NOISE,TrajectoryStateOnSurface());}
  const std::vector<const GeomDet * > & comp = geodet->components();
  if (_debug && talk) {for (std::vector<const GeomDet * >::const_iterator it = comp.begin(); it!=comp.end();it++) D((*it)->geographicalId ().rawId());}
  //the plane of the rechit
  const BoundPlane & plane = geodet->surface();
  
  //check if this is a 3D or a 2D measurement
  const SiStripMatchedRecHit2D * _3D=dynamic_cast<const SiStripMatchedRecHit2D *>(&rechit);
  const SiPixelRecHit * _PX = dynamic_cast<const SiPixelRecHit*>(&rechit);
  
  //purpose :
  //find the closest PSimHit on the module of the given rechit.
  //if none, it is noise
  //if one, it can be othertrack, secondary or primary
  //give a status about the corresponding PSimiHits within a local dX,dY cut
  const double dX_match = 0.05; //500 microns //that's all I can bear
  const double dY_match = 0.2; //2 millimeter in the y direction, because of the matcher reconstruction ?

  const PSimHit * bestmatch =NULL;
  TrajectoryStateOnSurface detTsos;
  double Xmin=dX_match;
  double Ymin=dY_match;

  //  what is the part of the tracker ?
  Handle<PSimHitContainer> simHitsXXX;
  int subdetid = rechit.geographicalId().subdetId();
  std::string suffix;
  if (Ltof) suffix= "LowTof";
  else suffix = "HighTof";

  switch (subdetid)
    {
    case 1:
      { _theEvent->getByLabel(_SimHitsSource, "TrackerHitsPixelBarrel"+suffix, simHitsXXX);      break; }
    case 2:
      { _theEvent->getByLabel(_SimHitsSource, "TrackerHitsPixelEndcap"+suffix, simHitsXXX);      break; }
    case 3:
      {	_theEvent->getByLabel(_SimHitsSource, "TrackerHitsTIB"+suffix, simHitsXXX);	break; }
    case 4 :
      {	_theEvent->getByLabel(_SimHitsSource, "TrackerHitsTID"+suffix, simHitsXXX);	break; }
    case 5 : 
      {	_theEvent->getByLabel(_SimHitsSource, "TrackerHitsTOB"+suffix, simHitsXXX);	break; }
    case 6:
      {	_theEvent->getByLabel(_SimHitsSource, "TrackerHitsTEC"+suffix, simHitsXXX);	break; }
    }

  if (! simHitsXXX.isValid() )
    {edm::LogError("IntrusiveAnalyzer::find")<< "No "<< modulename(rechit.geographicalId()) <<  "  sim hits found"; return make_pair(NOISE,TrajectoryStateOnSurface());}


  if (_debug && talk){
    Dn(rechit.geographicalId().rawId());
    Dn(modulename(rechit.geographicalId()));
    Dn(simHitsXXX->size());
    D(rechit.localPosition());}
  
  for ( PSimHitContainer::const_iterator Hits_CI = simHitsXXX->begin();
	Hits_CI != simHitsXXX->end() ; Hits_CI++) 
    {
      bool breakit= true;
      if (_debug && talk)  {Dn(Hits_CI->detUnitId());}
      //does it correspond to the det it or one of the components detid ?
      for (std::vector<const GeomDet * > ::const_iterator it = comp.begin();it!=comp.end();it++)
	{if ((*it)->geographicalId().rawId() == Hits_CI->detUnitId()){breakit=false; break;}}
      if (Hits_CI->detUnitId() == rechit.geographicalId().rawId()) {breakit = false ;}
      if (breakit) continue;


      //      if (Hits_CI->processType()!=2) continue; //this is very dangerous

      double abs_dx =0;
      double abs_dy =0;

      if (_3D || _PX)
	{//these are two 3D point that you are comparing ith: use transverse distance
	  // before that, propagate PSimHit to the surface

	  //get a tsos from PSmihit
	  TrajectoryStateOnSurface tsos = convert(*Hits_CI); 
	  if (tsos.globalMomentum().perp()<0.2) continue; //discard those with less than 200 MeV

	  if (!tsos.isValid())
	    {edm::LogError("IntrusiveAnalyzer::find")<<"converted PSimHit into TSOS is not valid"; continue;}

	  if (_debug && talk){show(tsos,"rsmsf::find  tsos:");}

	  if (Hits_CI->detUnitId() != rechit.geographicalId().rawId()) //propagate only if the surfaces are different
	    {
	      //the simhitpoint on the surface
	      detTsos = prop(anyDirection)->propagate(*tsos.freeState(),plane);
	      if (!detTsos.isValid())
		{ if (Hits_CI->processType()==2){edm::LogError("IntrusiveAnalyzer::find")<<" cannot propagate a simhit "<<tsos.globalMomentum()<<" from "<<Hits_CI->detUnitId()
												  <<" to the rechit surface id "<<rechit.geographicalId().rawId();
		    COMMENT(" cannot propagate a simhit");
		    /*
		      prop(anyDirection)->setDebug(true);
		      prop(anyDirection)->propagate(*tsos.freeState(),plane);
		      prop(anyDirection)->setDebug(false);
		    */
		  }
		  continue;}
	    }
	  else //surfaces are the same, this is it
	    detTsos = tsos;

	  if (_debug && talk) {show(detTsos,"rsmsf::find  detTsos:");}
	  LocalPoint simhitpoint= detTsos.localPosition();
	    
	  //	  D(detTsos.surfaceSide());

	  //use transverse distance between the points
	  abs_dx =abs((simhitpoint-rechit.localPosition()).x());
	  abs_dy =abs((simhitpoint-rechit.localPosition()).y());

	  if (_debug && talk) {
	    COMMENT(" matched rechit to be matched to a simulated hit, using distance in the transverse plane");
	    Dn(rechit.geographicalId().rawId());D(Hits_CI->detUnitId());
	    Dn(rechit.localPosition());Dn(simhitpoint);
	    D(plane.toGlobal(LocalPoint(0,0,0)));
		   
	    Dn(abs_dx);D((abs_dx < Xmin));
	    Dn(abs_dy);D((abs_dy < Ymin));	      }
	}
      else
	{
	  detTsos = convert(*Hits_CI);
	  LocalPoint Hits_CI_stripcenter = stripcenter (Hits_CI->entryPoint()
							,rechit.geographicalId());
	  // now that x coordinate is on the same footing.
	  abs_dx = fabs(Hits_CI_stripcenter.x() - rechit.localPosition().x());
	  if (_debug && talk){	Dn(Hits_CI->entryPoint());Dn(Hits_CI->exitPoint());Dn(Hits_CI->localDirection());Dn(Hits_CI->detUnitId());D(_trackID);
	    //	Dn(abs_dx);Dn(dX_match);Dn(hit_halfwidth);D((abs_dx > dX_match + hit_halfwidth ));
	    Dn(abs_dx);D((abs_dx < Xmin));
	  }
	}

      //check on Y only if 3D
      if ((abs_dx < Xmin) && ((abs_dy<Ymin) || (_3D==0 && _PX==0))  ) 
	{
	  Xmin = abs_dx;
	  Ymin = abs_dy;
	  bestmatch = &(*Hits_CI);
	  if (_debug && talk){
	    Dn(Hits_CI->detUnitId());
	    Dn(Hits_CI->trackId());
	    Dn(Hits_CI->processType());
	    Dn(Hits_CI->localPosition());
	    D(rechit.localPosition());
	  }

	}
    }

  if (!bestmatch) return make_pair(NOISE,TrajectoryStateOnSurface());

  if (bestmatch->trackId() == _trackID+ _simIndexFix)
    {
      //this is coming "from" the considered trackId
      unsigned short process = bestmatch->processType();
      if (process  == 2 )
	{ // this mean ionazition. straight through
	  return make_pair(PRIMARY,detTsos);}
      else
	{ // 9 would mean something
	  //this mean secondary hits, due to delta ray, etc
	  return make_pair(SECONDARY,detTsos);}
    }
  else
    { return make_pair(OTHERTRACK,detTsos); }

  //if nothing happened so far, then it is noise
  return make_pair(NOISE,TrajectoryStateOnSurface());
}


TransientTrackingRecHit::ConstRecHitPointer  IntrusiveAnalyzer:: match(const PSimHit & simhit)
{
  DetId did(simhit.detUnitId());
  DetId usethisId=did;

  //is it a mono layer, with a glued id too ?
  switch (did.subdetId()){
  case 3: { TIBDetId id(did); if (id.glued()!=0) { usethisId=DetId(id.glued());} break;}
  case 4: { TIDDetId id(did); if (id.glued()!=0) { usethisId=DetId(id.glued());} break;}
  case 5: { TOBDetId id(did); if (id.glued()!=0) { usethisId=DetId(id.glued());} break;}
  case 6: { TECDetId id(did); if (id.glued()!=0) { usethisId=DetId(id.glued());} break;}
  }
   
  //get the local center of the strip intercepted by the PSimHit
  LocalPoint center = stripcenter( simhit.entryPoint(), did);

  double dmin=10000;
  double abs_d=0;
  bool bestIs3D=false;
  TransientTrackingRecHit::ConstRecHitPointer best;
  // const    TrackingRecHit* best=NULL;
  //loop over the range of rechits on this module

  //convert the psimhit into a TSOS
  TrajectoryStateOnSurface TSOS = convert(simhit);
  if (!TSOS.isValid())
    { edm::LogError("IntrusiveAnalyzer:: match(const PSimHit & simhit)")<<" cannot build a TSOS from the simulted hit"; return best;     }

  //get the rechits for this detector id, glued if mono
  if (!_measurementTracker->idToDet(usethisId)){edm::LogError("IntrusiveAnalyzer:: match")<<usethisId.rawId()<<" is not recognized by the measurement tracker";}

  TransientTrackingRecHit::ConstRecHitContainer  thoseHits = _measurementTracker->idToDet(usethisId)->recHits(TSOS);
  for (TransientTrackingRecHit::ConstRecHitContainer::const_iterator iTThit = thoseHits.begin();
       iTThit != thoseHits.end(); iTThit++)
    {
      const TrackingRecHit * hit = (*iTThit)->hit();
      const SiStripMatchedRecHit2D * _3D = dynamic_cast<const SiStripMatchedRecHit2D *>(hit);
      const SiPixelRecHit * _PX = dynamic_cast<const SiPixelRecHit*>(hit);
      //what is the rechit type 
      if (_3D!=0 || _PX!=0)
	{//3D hits or pixel, compare distance in the module plane
	  const GeomDet * geodet = tracker()->idToDet(hit->geographicalId());
	  if (!geodet)
	    {edm::LogError("IntrusiveAnalyzer::match(...)")<<"no detlayer for id "<<hit->geographicalId().rawId()<<" from "<<did.rawId(); return NULL;}
	  const BoundPlane & plane = geodet->surface();
	   
	  TrajectoryStateOnSurface TSonsurface = prop(anyDirection)->propagate(*TSOS.freeState(),plane);

	  if (!TSonsurface.isValid())
	    {edm::LogError("IntrusiveAnalyzer::match(...)")<<"cannot propagate "<<TSOS.globalMomentum()<<" to "<<hit->geographicalId().rawId(); return NULL;}
	   
	  //abs_d = (simhit.entryPoint() - hit->localPosition()).perp();
	  abs_d = (TSonsurface.localPosition() - hit->localPosition()).perp();
	}
      else
	{//2D hit, stick with the transverse distance (local x)
	  abs_d = fabs(center.x() - hit->localPosition().x());}

      if (abs_d < dmin)
	{
	  bestIs3D=(_3D!=0 || _PX!=0);
	  dmin = abs_d;
	  //best = &(*hit);
	  best = *iTThit;
	}
    }
  if (best && _debug) {COMMENT("PSimHit with a spatial match to a tracking RecHit");Dn(bestIs3D);D(abs_d); }
  return best;
}


std::string  IntrusiveAnalyzer::modulename(const DetId det)
{
  std::stringstream ss;
  switch (det.subdetId())
    {
    case 1: {
      PXBDetId id(det);
      ss<<"PXB"<<id.layer();
      break; }
    case 2: {
      PXFDetId id(det);
      ss<<"PXF_"<<id.disk();
      break;}
    case 3: {
      TIBDetId id(det);
      ss<<"TIB_"<<id.layer();
      break;}
    case 5: {
      TOBDetId id(det);
      ss<<"TOB_"<<id.layer();
      break; }
    case 4:{
      TIDDetId id(det);
      ss<<"TID_"<<id.wheel();
      break; }
    case 6: {
      TECDetId id(det);
      ss<<"TEC_"<<id.wheel();
      break; }
    }
  return ss.str();
}

void IntrusiveAnalyzer::show(const DetId det,bool showsim,ostream & o){
  using namespace edm;
  if (det.det()!=1){    COMMENTf("not a tracker detid",o);    return; }
  
  std::string moduleName = modulename(det);
  Dfn(moduleName,o);
  switch (det.subdetId()) {
  case 1:
    { PXBDetId id(det);
      Dfn(id.rawId(),o);Dfn(id.subdetId(),o);Dfn(id.layer(),o);Df(id.module(),o);
      break;}
  case 2:
    { PXFDetId id(det);
      Dfn(id.rawId(),o);Dfn(id.subdetId(),o);Dfn(id.disk(),o);Df(id.module(),o);
      break;}
  case 3:
    {   TIBDetId id(det);
      Dfn(id.rawId(),o);Dfn(id.subdetId(),o);Dfn(id.layer(),o);Dfn(id.module(),o);Dfn(id.stereo(),o);Dfn(id.glued(),o);Df(id.partnerDetId(),o);
      break; }
  case 4:
    {   TIDDetId id(det);
      Dfn(id.rawId(),o);Dfn(id.subdetId(),o);Dfn(id.side(),o);Dfn(id.wheel(),o);Dfn(id.ring(),o);Dfn(id.stereo(),o);Dfn(id.glued(),o);Df(id.partnerDetId(),o);
      break; }
  case 5:
    {   TOBDetId id(det);
      Dfn(id.rawId(),o);Dfn(id.subdetId(),o);Dfn(id.layer(),o);Dfn(id.module(),o);Dfn(id.stereo(),o);Dfn(id.glued(),o);Df(id.partnerDetId(),o);
      break; }
  case 6:
    {   TECDetId id(det);
      Dfn(id.rawId(),o);Dfn(id.subdetId(),o);Dfn(id.side(),o);Dfn(id.wheel(),o);Dfn(id.ring(),o);Dfn(id.stereo(),o);Dfn(id.glued(),o);Df(id.partnerDetId(),o);
      break; }
  }


  const GeomDet * geomdet =tracker()->idToDet(det);
  if (!geomdet) {edm::LogError("Show(...)")<<"no geomdet for id "<<det.rawId(); return;}
  
  Df(geomdet->surface().toGlobal(LocalPoint(0,0,0)),o);
  Df(geomdet->surface().toGlobal(LocalVector(1,0,0)),o);
  Df(geomdet->surface().toGlobal(LocalVector(0,1,0)),o);
  Df(geomdet->surface().toGlobal(LocalVector(0,0,1)),o);
  
    if (showsim){
      const PSimHit * match_withoutId = simHitOnDetId(0,det,false);
      if (match_withoutId) {Dfn(match_withoutId->particleType(),o);Dfn(match_withoutId->trackId(),o);Df(match_withoutId->localPosition(),o);}
      else {COMMENTf(" not even a match without trackId requirement",o);} }
    return;
}



void IntrusiveAnalyzer::show(const TrajectoryStateOnSurface & TSOS, char * label,ostream & o)
{
  COMMENTf("--------------------------------------",o);
  COMMENTf("muon trajectory state on surface at "<<label,o);
  Dfn(TSOS.globalPosition(),o);Df(TSOS.globalPosition().mag(),o);
  Dfn(TSOS.globalMomentum(),o);Df(TSOS.globalMomentum().mag(),o);
  Dfn(TSOS.cartesianError().matrix(),o);
  COMMENTf("--------------------------------------",o);
}

void IntrusiveAnalyzer::show(const FreeTrajectoryState & FS,char * label,ostream & o)
{
  COMMENTf("--------------------------------------",o);
  COMMENTf(" muon free state at "<<label,o);
  Dfn(FS.position(),o);Df(FS.position().mag(),o);
  Dfn(FS.momentum(),o);Df(FS.momentum().mag(),o);
  //  Dn(FS.cartesianError().position().matrix());
  Dfn(FS.cartesianError().matrix(),o);
  COMMENTf("--------------------------------------",o);
}
