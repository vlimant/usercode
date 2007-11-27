#include "RecoTracker/TkNavigation/interface/LayerCollector.h"
#include "RecoTracker/TkNavigation/interface/StartingLayerFinder.h" 

#include <FWCore/MessageLogger/interface/MessageLogger.h>

using namespace std;


/* new implementation*/
vector<LayerCollector::LayerWithState> LayerCollector::allLayersWithState(const FTS& aFts) const {
  vector<LayerCollector::LayerWithState> result;

  FTS myFts = aFts;

  vector<const DetLayer*> nextLayers;
  vector<LayerCollector::LayerWithState> myLayersWithState = finder()->startingPixelLayerWithStates(myFts, deltaR(), deltaZ());
  for ( vector<StartingLayerFinder::LayerWithState>::iterator ilay = myLayersWithState.begin(); ilay != myLayersWithState.end(); ilay++) {
    nextLayers.push_back(ilay->first);}

  bool inside = true;
  TSOS pTsos;
  while(inside) {
    inside = false;//stop looking if none of the next layers are actually reached
    for(vector<const DetLayer*>::iterator ilay = nextLayers.begin(); ilay != nextLayers.end(); ilay++) {
      if (finder()->inLayer(myFts,pTsos,(*ilay),deltaR(),deltaZ())){
        //keep on looking for other layers
        inside = true;
        //gather this layer
	result.push_back(LayerCollector::LayerWithState(*ilay,pTsos));
	myFts = *pTsos.freeState();
	//get the next layers from this reached layer
        nextLayers = (**ilay).nextLayers(myFts,
                                         propagator()->propagationDirection());
	break;
      }
    }
  }

  LogDebug("TkNavigation")<<result.size()<<" layers collected";
  return result;
}

void LayerCollector::splittedLayersWithState(const FTS& aFts,
				    vector<LayerCollector::BarrelLayerWithState> &bLayers,
				    vector<LayerCollector::ForwardLayerWithState> &fLayers)const {
  vector<LayerCollector::LayerWithState> all = allLayersWithState(aFts);
  barrelLayersWithState(all,bLayers);
  forwardLayersWithState(all,fLayers);}

void LayerCollector::barrelLayersWithState(const vector<LayerCollector::LayerWithState> & all,
				  vector<LayerCollector::BarrelLayerWithState>& bLayers) const{
  for(vector<LayerCollector::LayerWithState>::const_iterator ilay = all.begin();ilay != all.end(); ilay++) {
    if (( ilay->first)->location()==GeomDetEnumerators::barrel)
      bLayers.push_back(LayerCollector::BarrelLayerWithState(dynamic_cast<const BarrelDetLayer*>(ilay->first),ilay->second)); 
}}

void LayerCollector::forwardLayersWithState(const vector<LayerCollector::LayerWithState> & all,
				   vector<LayerCollector::ForwardLayerWithState> &fLayers) const{
  for(vector<LayerCollector::LayerWithState>::const_iterator ilay = all.begin();ilay != all.end(); ilay++) {
    if ((ilay->first)->location()==GeomDetEnumerators::endcap)
      fLayers.push_back(LayerCollector::ForwardLayerWithState(dynamic_cast<const ForwardDetLayer*>(ilay->first),ilay->second)); 
}}

  
/* old implementation*/
vector<const DetLayer*> LayerCollector::allLayers(const FTS& aFts) const {
  
  vector<const DetLayer*> myLayers;
  
  //  FTS myFts(aFts.parameters());
  FTS myFts = aFts;

  //  vector<const DetLayer*> nextLayers = finder()->startingLayers(myFts, deltaR(), deltaZ());

  vector<const DetLayer*> nextLayers;
  vector<StartingLayerFinder::LayerWithState> nextLayersWithState = finder()->startingPixelLayerWithStates(myFts, deltaR(), deltaZ());
  for ( vector<StartingLayerFinder::LayerWithState>::iterator ilay = nextLayersWithState.begin(); ilay != nextLayersWithState.end(); ilay++) { 
    nextLayers.push_back(ilay->first);}
  
  bool inside = true;
  TSOS pTsos;
  while(inside) {
    inside = false;//stop looking if none of the next layers are actually reached
    for(vector<const DetLayer*>::iterator ilay = nextLayers.begin(); ilay != nextLayers.end(); ilay++) {
      if (finder()->inLayer(myFts,pTsos,(*ilay),deltaR(),deltaZ())){
	//keep on looking for other layers
	inside = true;
	//gather this layer
	myLayers.push_back(*ilay);
	//set the state on this layer as next starting point
	//	myFts = FTS(pTsos.globalParameters());
	myFts = *pTsos.freeState();
	//get the next layers from this reached layer
	nextLayers = (**ilay).nextLayers(myFts,
					 propagator()->propagationDirection());
	LogDebug(_category)<<nextLayers.size()<<" next layers."
			   <<"from :\n"
			   <<(*pTsos.freeState())
			   <<"going :"<<propagator()->propagationDirection();
	break;
      }
    }
  }
  return myLayers;
}


vector<const DetLayer*> LayerCollector::allOldLayers(const FTS& aFts) const {

  vector<const DetLayer*> myLayers;
  
  //  FTS myFts(aFts.parameters());
  FTS myFts = aFts;
  vector<const DetLayer*> nextLayers = finder()->startingLayers(myFts, deltaR(), deltaZ());
  
  bool inside = true;
  TSOS pTsos;
  while(inside) {
    inside = false;//stop looking if none of the next layers are actually reached
    for(vector<const DetLayer*>::iterator ilay = nextLayers.begin(); ilay != nextLayers.end(); ilay++) {

      TSOS pTsos = propagator()->propagate(myFts, (**ilay).surface());
      if(pTsos.isValid()) {
	inside = true;
	      
	if((**ilay).location() == GeomDetEnumerators::barrel) {
	      
	  Range barrZRange((**ilay).position().z() - 
			   0.5*((**ilay).surface().bounds().length()),
			   (**ilay).position().z() + 
			   0.5*((**ilay).surface().bounds().length()));
	  Range trajZRange(pTsos.globalPosition().z() - deltaZ(),
			   pTsos.globalPosition().z() + deltaZ());
	      
	  if(rangesIntersect(trajZRange, barrZRange)) 
	    myLayers.push_back(*ilay);
	      
	} else if((**ilay).location() == GeomDetEnumerators::endcap) {

	  const ForwardDetLayer* fwd = 
	    dynamic_cast<const ForwardDetLayer*>(*ilay);
	  Range fwdRRange((*fwd).specificSurface().innerRadius(),
			  (*fwd).specificSurface().outerRadius());
	  Range trajRRange(pTsos.globalPosition().perp() - deltaR(),
			   pTsos.globalPosition().perp() + deltaR());
	      
	  if(rangesIntersect(trajRRange, fwdRRange)) 
	    myLayers.push_back(*ilay);
	      
	}
	//	myFts = FTS(pTsos.globalParameters());
	myFts = *pTsos.freeState();
	nextLayers = (**ilay).nextLayers(*pTsos.freeState(), 
					 propagator()->propagationDirection());
	LogDebug(_category)<<nextLayers.size()<<" next layers (old)."
			   <<"from :\n"
			   <<(*pTsos.freeState())
			   <<"going :"<<propagator()->propagationDirection();
	break;
      }     
    }
  }
  
  return myLayers;
}



void LayerCollector::splittedLayers(const FTS& aFts,
				   vector<const BarrelDetLayer*> &bLayers,
				   vector<const ForwardDetLayer*> &fLayers)const {
  vector<const DetLayer*> all = allLayers(aFts);
  barrelLayers(all,bLayers);
  forwardLayers(all,fLayers);}

void LayerCollector::barrelLayers(const vector<const DetLayer*> & all,
				  vector<const BarrelDetLayer*>& bLayers) const{
  for(vector<const DetLayer*>::const_iterator ilay = all.begin();ilay != all.end(); ilay++) {
    if (( *ilay)->location()==GeomDetEnumerators::barrel)
      bLayers.push_back(dynamic_cast<const BarrelDetLayer*>(*ilay)); 
  }}

void LayerCollector::forwardLayers(const vector<const DetLayer*> & all,
				   vector<const ForwardDetLayer*> &fLayers) const{
  for(vector<const DetLayer*>::const_iterator ilay = all.begin();ilay != all.end(); ilay++) {
    if (( *ilay)->location()==GeomDetEnumerators::endcap)
      fLayers.push_back(dynamic_cast<const ForwardDetLayer*>(*ilay)); 
  }}


vector<const BarrelDetLayer*> LayerCollector::barrelLayers(const FTS& aFts) const {

  vector<const DetLayer*> all = allLayers(aFts);
  vector<const BarrelDetLayer*> bLayers;
  barrelLayers(all,bLayers);
  return bLayers;}
  
vector<const ForwardDetLayer*> LayerCollector::forwardLayers(const FTS& aFts) const {
  
  vector<const DetLayer*> all = allLayers(aFts);
  vector<const ForwardDetLayer*> fLayers;
  forwardLayers(all,fLayers);
  return fLayers;}


