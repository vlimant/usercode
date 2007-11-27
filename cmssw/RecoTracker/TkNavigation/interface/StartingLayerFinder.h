#ifndef TkNavigation_StartingLayerFinder_H_
#define TkNavigation_StartingLayerFinder_H_


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"


#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"

#include <vector>
#include <string>

/** Finds the nearest navigable layer.
 *  Needed to start trajectory building in case the seed does not 
 *  have a DetLayer
 */

class PTrajectoryStateOnDet;

class StartingLayerFinder {
private:

  typedef FreeTrajectoryState FTS;
  typedef TrajectoryStateOnSurface TSOS;
  typedef std::pair<float, float> Range;

public: 
  typedef std::pair< const DetLayer * , TrajectoryStateOnSurface > LayerWithState;


  StartingLayerFinder(const Propagator* aPropagator, const MeasurementTracker*  tracker ) :
  thePropagator(aPropagator),
    theEstimator(NULL),
    theMeasurementTracker(tracker),
       thePixelLayersValid(false),
    theFirstPixelBarrelLayer(0),
    theFirstNegPixelFwdLayer(0),
    theFirstPosPixelFwdLayer(0),
       theFirstStripLayersValid(false),
    theFirstStripBarrelLayer(0),
    theFirstNegStripFwdLayer(0),
    theFirstPosStripFwdLayer(0),
      theLastStripLayersValid(false),
    theLastStripBarrelLayer(0),
    theLastNegStripFwdLayer(0),
      theLastPosStripFwdLayer(0) 
      {_category="StartingLayerFinder";}
  
  StartingLayerFinder(const Propagator* aPropagator, const MeasurementEstimator* estimator, const MeasurementTracker*  tracker ) :
  thePropagator(aPropagator),
    theEstimator(estimator),
    theMeasurementTracker(tracker),
       thePixelLayersValid(false),
    theFirstPixelBarrelLayer(0),
    theFirstNegPixelFwdLayer(0),
    theFirstPosPixelFwdLayer(0),
       theFirstStripLayersValid(false),
    theFirstStripBarrelLayer(0),
    theFirstNegStripFwdLayer(0),
    theFirstPosStripFwdLayer(0),
      theLastStripLayersValid(false),
    theLastStripBarrelLayer(0),
    theLastNegStripFwdLayer(0),
      theLastPosStripFwdLayer(0) 
      {_category="StartingLayerFinder";}

  ~StartingLayerFinder() {}

  enum Dest { Pixel, firstStrip, lastStrip };

  //vector of LayerWithState
  std::vector<LayerWithState> startingPixelLayerWithStates(const FTS& aFts, float dr=1, float dz=1)const;
  std::vector<LayerWithState> startingStripLayerWithStates(const FTS& aFts, float dr=1, float dz=1)const;
  std::vector<LayerWithState> startingOuterStripLayerWithStates(const FTS& aFts, float dr=1, float dz=1)const;

  std::vector<LayerWithState> startingLayerWithStates(const FTS& aFts, float dr, float dz,
						      const BarrelDetLayer* barrel,
						      const std::vector<ForwardDetLayer*> & negDisk,
						      const std::vector<ForwardDetLayer*> & posDisk,
						      bool stopAtMissedBarrel=true,
						      bool stopAtMissedDisk=true,
						      bool onlyOneLayer=false) const;

  //explicit methods with only one layerwithstate (internaly use the vector output method)
  LayerWithState startingStripLayerWithState(const FTS& aFts, float dr=1, float dz=1)const;
  LayerWithState startingPixelLayerWithState(const FTS& aFts, float dr=1, float dz=1)const;
  LayerWithState startingOuterStripLayerWithState(const FTS& aFts, float dr=1, float dz=1)const;

  /* //output only the detlayer //useless
     const DetLayer* startingStripLayer(const FTS& aFts, float dr=1, float dz=1)const { return startingStripLayerWithState(aFts,dr,dz).first;}
     const DetLayer* startingPixelLayer(const FTS& aFts, float dr=1, float dz=1)const { return startingPixelLayerWithState(aFts,dr,dz).first;}
     const DetLayer* startingOuterStripLayer(const FTS& aFts, float dr=1, float dz=1)const  { return startingOuterStripLayerWithState(aFts,dr,dz).first;} */

  /* should be removed from the interface */
  //backward compatibility
  //search for layers is made differently here
  std::vector<const DetLayer*> startingLayers(const FTS& aFts, float dr, float dz,Dest=Pixel) const;
  std::vector<const DetLayer*> startingLayers(const TrajectorySeed& aSeed,Dest=Pixel) const;

  //vector of detlayer instead of LayerWithState
  std::vector<const DetLayer*> startingStripLayers(const FTS& aFts, float dr=1, float dz=1)const;
  std::vector<const DetLayer*> startingPixelLayers(const FTS& aFts, float dr=1, float dz=1)const;
  std::vector<const DetLayer*> startingOuterStripLayers(const FTS& aFts, float dr=1, float dz=1)const;
  
  //generic methods to output a vector of DetLayer
  std::vector<const DetLayer*> startingLayers(const FTS& aFts, float dr, float dz,
					      const BarrelDetLayer* barrel,
					      const std::vector<ForwardDetLayer*> & negDisk,
					      const std::vector<ForwardDetLayer*> & posDisk) const;
  /* until here */


  //geometry member functions
  //pixel layer methods
  const BarrelDetLayer* firstPixelBarrelLayer() const;
  const std::vector<ForwardDetLayer*> & firstNegPixelFwdLayer() const;
  const std::vector<ForwardDetLayer*> & firstPosPixelFwdLayer() const;

  //TIB/TID layer methods
  const BarrelDetLayer* firstStripBarrelLayer() const;
  const std::vector<ForwardDetLayer*> & firstNegStripFwdLayer() const;
  const std::vector<ForwardDetLayer*> & firstPosStripFwdLayer() const;

  //TOB/TEC layer methods
  const BarrelDetLayer* lastStripBarrelLayer() const;
  const std::vector<ForwardDetLayer*> & lastNegStripFwdLayer() const;
  const std::vector<ForwardDetLayer*> & lastPosStripFwdLayer() const;

  const Propagator* propagator() const {return thePropagator;}
  const MeasurementEstimator* estimator() const {return theEstimator;}

  bool inLayer(const FTS& fastFts, TSOS & pTsos, const DetLayer * layer,float dr,float dz, bool onSurface=false) const;
  bool inDisk(const FTS& fastFts, TSOS & pTsos, const ForwardDetLayer * layer,float dr, bool onSurface=false) const;
  bool inBarrel(const FTS& fastFts, TSOS & pTsos, const BarrelDetLayer * layer,float dz, bool onSurface=false) const;
  

private:

  std::string _category;

  const Propagator* thePropagator;
  const MeasurementEstimator* theEstimator;
  const MeasurementTracker* theMeasurementTracker;

  mutable bool thePixelLayersValid;
  mutable BarrelDetLayer* theFirstPixelBarrelLayer;
  mutable std::vector<ForwardDetLayer*> theFirstNegPixelFwdLayer;
  mutable std::vector<ForwardDetLayer*> theFirstPosPixelFwdLayer;

  mutable bool theFirstStripLayersValid;
  mutable BarrelDetLayer* theFirstStripBarrelLayer;
  mutable std::vector<ForwardDetLayer*> theFirstNegStripFwdLayer;
  mutable std::vector<ForwardDetLayer*> theFirstPosStripFwdLayer;

  mutable bool theLastStripLayersValid;
  mutable BarrelDetLayer* theLastStripBarrelLayer;
  mutable std::vector<ForwardDetLayer*> theLastNegStripFwdLayer;
  mutable std::vector<ForwardDetLayer*> theLastPosStripFwdLayer;


  void checkPixelLayers() const;
  void checkFirstStripLayers() const;
  void checkLastStripLayers() const;
  
  /*
    enum diskIntersect { below , in , above ,out };
    inline diskIntersect diskRangeIntersect(const Range& a, const Range& b) const {
    if (a.first > b.second) return below;
    if (b.first > a.second) return above;
    return in;
    }*/

  inline bool rangesIntersect( const Range& a, const Range& b) const {
    if ( a.first > b.second || b.first > a.second) return false;
    else return true;
  }

  /*  enum intersect { noneOut, noneIn, yes };
      inline intersect rangesIntersection( const Range& a, const Range& b) const {     
      if ( a.first > b.second || b.first > a.second) return yes; 
      else {
      if (a.second < b.first) return noneIn;
      else return noneOut;}
  }*/



};
#endif //TR_StartingLayerFinder_H_
