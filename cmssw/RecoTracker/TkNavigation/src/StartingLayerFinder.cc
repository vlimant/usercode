#include "RecoTracker/TkNavigation/interface/StartingLayerFinder.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <utility>

using namespace std;


std::vector<StartingLayerFinder::LayerWithState>
StartingLayerFinder::startingPixelLayerWithStates(const FTS& aFts, float dr, float dz)const{
  return  startingLayerWithStates(aFts,dr,dz,
			 firstPixelBarrelLayer(),
			 firstNegPixelFwdLayer(),
			 firstPosPixelFwdLayer());}


std::vector<StartingLayerFinder::LayerWithState>
StartingLayerFinder::startingStripLayerWithStates(const FTS& aFts, float dr, float dz)const{
  return  startingLayerWithStates(aFts,dr,dz,
			 firstStripBarrelLayer(),
			 firstNegStripFwdLayer(),
			 firstPosStripFwdLayer());}


std::vector<StartingLayerFinder::LayerWithState>
StartingLayerFinder::startingOuterStripLayerWithStates(const FTS& aFts, float dr, float dz)const{
  return  startingLayerWithStates(aFts,dr,dz,
				  lastStripBarrelLayer(),
				  lastNegStripFwdLayer(),
				  lastPosStripFwdLayer(),
				  false,false);}


std::vector<StartingLayerFinder::LayerWithState> 
StartingLayerFinder::startingLayerWithStates(const FTS& aFts, float dr, float dz,
					     const BarrelDetLayer* barrel,
					     const std::vector<ForwardDetLayer*> & nfwd,
					     const std::vector<ForwardDetLayer*> & pfwd,
					     bool stopAtMissedBarrel,
					     bool stopAtMissedDisk,
					     bool onlyOneLayer) const {

  LogDebug(_category)<<"StartingLayerFinder: \n"
		     <<"stopAtMissedBarrel("<<stopAtMissedBarrel<<")\n"
		     <<"stopAtMissedDisk("<<stopAtMissedDisk<<")\n"
		     <<"onlyOneLayer("<<onlyOneLayer<<")\n"
		     <<"dr="<<dr<<" dz="<<dz<<"\n"
		     <<"state "<<aFts; 
 
  std::vector<LayerWithState> mylayers; 
 
  FTS fastFts=aFts;
  //  FTS fastFts(aFts.parameters());
  
  TSOS pTsos;

  LogDebug(_category)<<"checking on barrel";
  //barrel detector
  if (inBarrel(fastFts,pTsos,barrel,dz)){
    mylayers.push_back(LayerWithState(barrel,pTsos));
    if (onlyOneLayer){LogDebug(_category)<<"only barrel layer"; return mylayers;}
    LogDebug(_category)<<"barrel layer is reached";}
  
  //something that does not reach the barrel is doomed later on
  if (stopAtMissedBarrel && !pTsos.isValid()) {
    LogDebug(_category)<<"state did not reached the barrrel surface.";return mylayers;}
  
  if (aFts.momentum().z()<0){
    LogDebug(_category)<<"checking on negative disk";
    
    //negative fwd detector
    for(std::vector<ForwardDetLayer*>::const_iterator infwd = nfwd.begin();
	infwd != nfwd.end(); infwd++) {

      /* implementation that would stop if your are below the disk inner radius of a reached disk
      diskIntersect inter=inDisk_(fastFts,pTsos,(*infwd),dr);
      if (inter==in){
	mylayers.push_back(LayerWithState(*infwd,pTsos));
        if (onlyOneLayer) {LogDebug(_category)<<"only negative forward layer";return mylayers;}
	//stop as soon as something is consistent
        break;}
      else if (inter=above){
	//you can still go to the next one man !
	continue;}
      else if(inter==below){
	//there is no hope anymore
	break;}
      else {
	//if the disk hasen't been reached at all. stop here
        if (!pTsos.isValid() && stopAtMissedDisk){break;}}*/

	
      
      if (inDisk(fastFts,pTsos,(*infwd),dr)){
	mylayers.push_back(LayerWithState(*infwd,pTsos));
	if (onlyOneLayer) {LogDebug(_category)<<"only negative forward layer";return mylayers;}
	LogDebug(_category)<<"a negative forward layer is reached";
	//stop as soon as something is consistent
	break;}
      else{
	//if the disk hasen't been reached at all. stop here
	if (!pTsos.isValid() && stopAtMissedDisk){break;}}
    }
  }

  if (aFts.momentum().z()>0){
    LogDebug(_category)<<"checking on positive disk";

    //positive fwd detector
    for(std::vector<ForwardDetLayer*>::const_iterator ipfwd = pfwd.begin();
	ipfwd != pfwd.end(); ipfwd++) {
      if (inDisk(fastFts,pTsos,(*ipfwd),dr)){
	mylayers.push_back(LayerWithState(*ipfwd,pTsos));
	if (onlyOneLayer){LogDebug(_category)<<"only positive forward layer"; return mylayers;}
	LogDebug(_category)<<"a positive forward layer is reached";
	//stop as soon as something is consistent
	break;}
      else{
	//if the disk hasen't been reached at all. stop here
	if (!pTsos.isValid() && stopAtMissedDisk){break;}}
    }
  }

  LogDebug(_category)<<" returning "<<mylayers.size()<<" starting layers";

  return mylayers;
}


StartingLayerFinder::LayerWithState 
StartingLayerFinder::startingPixelLayerWithState(const FTS& aFts, float dr, float dz)const{
  std::vector<LayerWithState> vect=startingLayerWithStates(aFts,dr,dz,
							   firstPixelBarrelLayer(),
							   firstNegPixelFwdLayer(),
							   firstPosPixelFwdLayer(),
							   true,true,false);
  if (vect.empty()) return LayerWithState(NULL,TrajectoryStateOnSurface()); else return vect.front();}

StartingLayerFinder::LayerWithState
StartingLayerFinder::startingStripLayerWithState(const FTS& aFts, float dr, float dz)const{
  std::vector<LayerWithState> vect=startingLayerWithStates(aFts,dr,dz,
							   firstStripBarrelLayer(),
							   firstNegStripFwdLayer(),
							   firstPosStripFwdLayer(),
							   true,true,false);
  if (vect.empty()) return LayerWithState(NULL,TrajectoryStateOnSurface()); else return vect.front();}

StartingLayerFinder::LayerWithState
StartingLayerFinder::startingOuterStripLayerWithState(const FTS& aFts, float dr, float dz)const{
  std::vector<LayerWithState> vect=startingLayerWithStates(aFts,dr,dz,
							   lastStripBarrelLayer(),
							   lastNegStripFwdLayer(),
							   lastPosStripFwdLayer(),
							   false,false,false);
  if (vect.empty()) return LayerWithState(NULL,TrajectoryStateOnSurface()); else return vect.front();}

bool StartingLayerFinder::
inLayer(const FTS& fastFts, TSOS & pTsos, const DetLayer * layer,float dr,float dz, bool onSurface) const{
  if (layer->location() == GeomDetEnumerators::barrel)
    {   return inBarrel(fastFts,pTsos,dynamic_cast<const BarrelDetLayer*>(layer),dz,onSurface);}
  else if (layer->location() ==  GeomDetEnumerators::endcap)
    {   return inDisk(fastFts,pTsos,dynamic_cast<const ForwardDetLayer*>(layer),dr,onSurface);}
  //by default
  return false;
}

bool StartingLayerFinder::
inDisk(const FTS& fastFts, TSOS & pTsos, const ForwardDetLayer * layer,float dr,bool onSurface) const{
  //  return (inDisk_(fastFts,pTsos,layer,dr,onSurface)==in);
  
  LogDebug(_category)<<"propagating to disk";
      if (!onSurface)
      pTsos = propagator()->propagate(fastFts, layer->specificSurface());
      
      if(pTsos.isValid()) {
      if (estimator()){
      LogDebug(_category)<<"using estimator";
      std::pair<bool,TrajectoryStateOnSurface> compatibility = layer->compatible(pTsos,*propagator(),*estimator());
      if (compatibility.first) return true;}
      else{
      
      Range nfwdRRange((layer)->specificSurface().innerRadius(),
      (layer)->specificSurface().outerRadius());       
      Range trajRRange(pTsos.globalPosition().perp() - dr,
      pTsos.globalPosition().perp() + dr);
      
      if(rangesIntersect(trajRRange, nfwdRRange)) {
      return true;}
      }}
      else{LogDebug(_category)<<"pTsos is not valid at disk surface";}
      return false;
}

/*
diskIntersect StartingLayerFinder::
inDisk_(const FTS& fastFts, TSOS & pTsos, const ForwardDetLayer * layer,float dr,bool onSurface) const{
  LogDebug(_category)<<"propagating to disk";
  if (!onSurface)
    pTsos = propagator()->propagate(fastFts, layer->specificSurface());

  if(pTsos.isValid()) {
    if (estimator()){
      LogDebug(_category)<<"using estimator";
      std::pair<bool,TrajectoryStateOnSurface> compatibility = layer->compatible(pTsos,*propagator(),*estimator());
      if (compatibility.first) return in;}
    else{

      Range nfwdRRange((layer)->specificSurface().innerRadius(),
                       (layer)->specificSurface().outerRadius());
      Range trajRRange(pTsos.globalPosition().perp() - dr,
                       pTsos.globalPosition().perp() + dr);

      //      diskIntercect inter=diskRangeIntercect(trajRRange, nfwdRRange);
      return diskRangeIntercect(trajRRange, nfwdRRange);
    }
  }
  else{LogDebug(_category)<<"pTsos is not valid at disk surface";}
  return out;
}
*/


bool StartingLayerFinder::
inBarrel(const FTS& fastFts, TSOS & pTsos, const BarrelDetLayer * layer,float dz,bool onSurface) const{
  LogDebug(_category)<<"propagating to barrel";
  if (!onSurface)
    pTsos = propagator()->propagate(fastFts, layer->specificSurface());

  if(pTsos.isValid()) {
    if (estimator()){
      LogDebug(_category)<<"using estimator";
      std::pair<bool,TrajectoryStateOnSurface> compatibility = layer->compatible(pTsos,*propagator(),*estimator());
      if (compatibility.first) return true;}
    else{
    
      Range barrZRange(layer->position().z() -
		       0.5*(layer->specificSurface().bounds().length()),
		       layer->position().z() +
		       0.5*(layer->specificSurface().bounds().length()));
      Range trajZRange(pTsos.globalPosition().z() - dz,
		       pTsos.globalPosition().z() + dz);

    if(rangesIntersect(trajZRange, barrZRange)) {
      return true;}
    }
  }
  else{LogDebug(_category)<<"pTsos is not valid at barrel surface";}
  return false;
}




//backward compatibility member functions:
//a different way of finding stating layers is implemented

std::vector<const DetLayer*> 
StartingLayerFinder::startingLayers(const FTS& aFts, float dr, float dz, Dest dest) const {
  switch (dest){
  case Pixel:
    return startingPixelLayers(aFts,dr,dz);
  case firstStrip:
    return startingStripLayers(aFts,dr,dz);
  case lastStrip:
    return startingOuterStripLayers(aFts,dr,dz);
  }
  return std::vector<const DetLayer*>();
}

std::vector<const DetLayer*>
StartingLayerFinder::startingPixelLayers(const FTS& aFts, float dr, float dz)const{
  return  startingLayers(aFts,dr,dz,
			 firstPixelBarrelLayer(),
			 firstNegPixelFwdLayer(),
			 firstPosPixelFwdLayer());}


std::vector<const DetLayer*>
StartingLayerFinder::startingStripLayers(const FTS& aFts, float dr, float dz)const{
  return  startingLayers(aFts,dr,dz,
			 firstStripBarrelLayer(),
			 firstNegStripFwdLayer(),
			 firstPosStripFwdLayer());}


std::vector<const DetLayer*>
StartingLayerFinder::startingOuterStripLayers(const FTS& aFts, float dr, float dz)const{
  return  startingLayers(aFts,dr,dz,
			 lastStripBarrelLayer(),
			 lastNegStripFwdLayer(),
			 lastPosStripFwdLayer());}



std::vector<const DetLayer*> StartingLayerFinder::startingLayers(const FTS& aFts, float dr, float dz,
								 const BarrelDetLayer* barrel,
								 const std::vector<ForwardDetLayer*> & nfwd_bla,
								 const std::vector<ForwardDetLayer*> & pfwd_bla) const {

  std::vector<const DetLayer*> mylayers; 
  //  mylayers.reserve(3);


  //  FTS fastFts(aFts.parameters());
  FTS fastFts=aFts;

  /*
  LogDebug(_category)<<"StartingLayerFinder: (orignal-modified) \n"
		     <<"dr="<<dr<<" dz="<<dz<<"\n"
		     <<"state "<<aFts; 
    
  */  
  /* more  way compact implementation
  TSOS pTsos;
  //barrel detector
  if (inBarrel(fastFts,pTsos,barrel,dz)){
    LogDebug(_category)<<"adding the barrel";
    mylayers.push_back(barrel);}


  //negative fwd detector
  for(std::vector<ForwardDetLayer*>::const_iterator infwd = nfwd.begin();
      infwd != nfwd.end(); infwd++) {
    if (inDisk(fastFts,pTsos,(*infwd),dr)){
      LogDebug(_category)<<"adding a negative disk";
      mylayers.push_back(*infwd);
    }
  }
     
  //positive fwd detector
  for(std::vector<ForwardDetLayer*>::const_iterator ipfwd = pfwd.begin();
      ipfwd != pfwd.end(); ipfwd++) {
    if (inDisk(fastFts,pTsos,(*ipfwd),dr)){
      LogDebug(_category)<<"adding a positive disk";
      mylayers.push_back(*ipfwd);
    }
  }
*/

  LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
		     <<"dr="<<dr<<" dz="<<dz<<"\n"
		     <<"state "<<aFts; 

  /* original implementation */
  //barrel pixel
  TSOS pTsos = 
    propagator()->propagate(fastFts, firstPixelBarrelLayer()->surface());
  
  if(pTsos.isValid()) {
    
    Range barrZRange(firstPixelBarrelLayer()->position().z() - 
		     0.5*(firstPixelBarrelLayer()->surface().bounds().length()),
		     firstPixelBarrelLayer()->position().z() + 
		     0.5*(firstPixelBarrelLayer()->surface().bounds().length()));
    Range trajZRange(pTsos.globalPosition().z() - dz,
		     pTsos.globalPosition().z() + dz);
    
    if(rangesIntersect(trajZRange, barrZRange)) {
      mylayers.push_back(firstPixelBarrelLayer());
	LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
			   <<"adding the barrel";
    }
  }
  else { LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
			    <<"state at barrel not valid";}
  
  
  //negative fwd pixel
  
  const vector<ForwardDetLayer*> nfwd = firstPosPixelFwdLayer();
  for(vector<ForwardDetLayer*>::const_iterator infwd = nfwd.begin();
      infwd != nfwd.end(); infwd++) {
    pTsos = propagator()->propagate(fastFts, (*infwd)->surface());  
    if(pTsos.isValid()) {
      Range nfwdRRange((*infwd)->specificSurface().innerRadius(),
		       (*infwd)->specificSurface().outerRadius());
      Range trajRRange(pTsos.globalPosition().perp() - dr,
		       pTsos.globalPosition().perp() + dr);
      if(rangesIntersect(trajRRange, nfwdRRange)) {
	LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
			   <<"adding a negative disk";
	mylayers.push_back(*infwd);
	
      }
    }
    else { LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
			      <<"state at negative disk not valid";}

  }
  
  //positive fwd pixel
  const vector<ForwardDetLayer*> pfwd = firstPosPixelFwdLayer();
  for(vector<ForwardDetLayer*>::const_iterator ipfwd = pfwd.begin();
      ipfwd != pfwd.end(); ipfwd++) {
    pTsos = propagator()->propagate(fastFts, (*ipfwd)->surface());
    if(pTsos.isValid()) {
      Range pfwdRRange((*ipfwd)->specificSurface().innerRadius(),
		       (*ipfwd)->specificSurface().outerRadius());
      Range trajRRange(pTsos.globalPosition().perp() - dr,
		       pTsos.globalPosition().perp() + dr);
      if(rangesIntersect(trajRRange, pfwdRRange)) {
	LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
			   <<"adding a positive disk";
	mylayers.push_back(*ipfwd);
	
      }
    }
    else { LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
			      <<"state at positive disk not valid";}

  }

  LogDebug(_category)<<"StartingLayerFinder: (orignal) \n"
		     <<"returning "<<mylayers.size()<<" layers";

  return mylayers;
}







std::vector<const DetLayer*> 
StartingLayerFinder::startingLayers(const TrajectorySeed& aSeed, Dest dest) const {



  float dr = 0., dz = 0.;


  if(propagator()->propagationDirection() != aSeed.direction())
    return std::vector<const DetLayer*>();

  if(aSeed.nHits() != 2) return std::vector<const DetLayer*>();
 

  TrackingRecHitCollection::const_iterator firstHit= aSeed.recHits().first;
  const TrackingRecHit* recHit1=&(*firstHit);
  const DetLayer* hit1Layer = theMeasurementTracker->geometricSearchTracker()->detLayer(recHit1->geographicalId());

  TrackingRecHitCollection::const_iterator secondHit= aSeed.recHits().second;
  const TrackingRecHit* recHit2=&(*secondHit);
  const DetLayer* hit2Layer = theMeasurementTracker->geometricSearchTracker()->detLayer(recHit2->geographicalId());

  
  GeomDetEnumerators::Location p1 =  hit1Layer->location();
  GeomDetEnumerators::Location p2 =  hit2Layer->location();

  if(p1 == GeomDetEnumerators::barrel && p2 == GeomDetEnumerators::barrel) {
    dr = 0.1; dz = 5.;
  } else if(p1 == GeomDetEnumerators::endcap && p2 == GeomDetEnumerators::endcap) {
    dr = 5.; dz = 0.1;
  } else {
    dr = 0.1; dz = 0.1;
  }


  
  const GeomDet* gdet = theMeasurementTracker->geomTracker()->idToDet( DetId( aSeed.startingState().detId()));
  
  TrajectoryStateTransform tsTransform;
  TrajectoryStateOnSurface tsos = tsTransform.transientState( aSeed.startingState(), &(gdet->surface()), 
							      thePropagator->magneticField());


  FreeTrajectoryState* fts=tsos.freeTrajectoryState();
  
  return startingLayers(*fts, dr, dz, dest);
}







//geometry member functions 
const BarrelDetLayer* StartingLayerFinder::firstPixelBarrelLayer() const {
  checkPixelLayers();  return theFirstPixelBarrelLayer;}

const std::vector<ForwardDetLayer*> & StartingLayerFinder::firstNegPixelFwdLayer() const {
  checkPixelLayers();  return theFirstNegPixelFwdLayer;}

const std::vector<ForwardDetLayer*> &  StartingLayerFinder::firstPosPixelFwdLayer() const {
  checkPixelLayers();  return theFirstPosPixelFwdLayer;}

void StartingLayerFinder::checkPixelLayers() const {


  if(!thePixelLayersValid) {
   
    const GeometricSearchTracker* theGeometricSearchTracker=theMeasurementTracker->geometricSearchTracker();

   
    theFirstPixelBarrelLayer = theGeometricSearchTracker->pixelBarrelLayers().front();
    theFirstNegPixelFwdLayer = theGeometricSearchTracker->negPixelForwardLayers();
    theFirstPosPixelFwdLayer = theGeometricSearchTracker->posPixelForwardLayers();
    thePixelLayersValid = true;
  }
}

 
const BarrelDetLayer* StartingLayerFinder::firstStripBarrelLayer() const {
  checkFirstStripLayers();  return theFirstStripBarrelLayer;}

const std::vector<ForwardDetLayer*> & StartingLayerFinder::firstNegStripFwdLayer() const {
  checkFirstStripLayers();  return theFirstNegStripFwdLayer;}

const std::vector<ForwardDetLayer*> & StartingLayerFinder::firstPosStripFwdLayer() const {
  checkFirstStripLayers();  return theFirstPosStripFwdLayer;}

void StartingLayerFinder::checkFirstStripLayers() const {


  if(!theFirstStripLayersValid) {
   
    const GeometricSearchTracker* theGeometricSearchTracker=theMeasurementTracker->geometricSearchTracker();

   
    theFirstStripBarrelLayer = theGeometricSearchTracker->tibLayers().front();
    std::vector<ForwardDetLayer*> tec=theGeometricSearchTracker->negTecLayers();
    theFirstNegStripFwdLayer = theGeometricSearchTracker->negTidLayers();
    theFirstNegStripFwdLayer.insert(theFirstNegStripFwdLayer.end(),tec.begin(),tec.end());
    tec=theGeometricSearchTracker->posTecLayers();
    theFirstPosStripFwdLayer = theGeometricSearchTracker->posTidLayers();
    theFirstPosStripFwdLayer.insert(theFirstPosStripFwdLayer.end(),tec.begin(),tec.end());
    theFirstStripLayersValid = true;
  }
}

 
const BarrelDetLayer* StartingLayerFinder::lastStripBarrelLayer() const {
  checkLastStripLayers();  return theLastStripBarrelLayer;}

const std::vector<ForwardDetLayer*> & StartingLayerFinder::lastNegStripFwdLayer() const {
  checkLastStripLayers();  return theLastNegStripFwdLayer;}

const std::vector<ForwardDetLayer*> & StartingLayerFinder::lastPosStripFwdLayer() const {
  checkLastStripLayers();  return theLastPosStripFwdLayer;}

void StartingLayerFinder::checkLastStripLayers() const {


  if(!theLastStripLayersValid) {

    //get first strip layers
    checkFirstStripLayers();

    const GeometricSearchTracker* theGeometricSearchTracker=theMeasurementTracker->geometricSearchTracker();

    //there is a possible crack in the tracker here between TOB6 and TEC1...
    theLastStripBarrelLayer = theGeometricSearchTracker->tobLayers().back();

    theLastNegStripFwdLayer =theFirstNegStripFwdLayer;
    std::reverse(theLastNegStripFwdLayer.begin(),theLastNegStripFwdLayer.end());
    theLastPosStripFwdLayer = theFirstPosStripFwdLayer;
    std::reverse(theLastPosStripFwdLayer.begin(),theLastPosStripFwdLayer.end());
    theLastStripLayersValid = true;
  }
}


