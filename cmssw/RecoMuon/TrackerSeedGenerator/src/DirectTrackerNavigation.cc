#include "RecoMuon/TrackerSeedGenerator/interface/DirectTrackerNavigation.h"

/** \file DirectTrackerNavigation
 *
 *  $Date: 2007/05/11 15:04:46 $
 *  $Revision: 1.1 $
 *  \author Chang Liu  -  Purdue University
 */

#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <algorithm>

using namespace std;

DirectTrackerNavigation::DirectTrackerNavigation(const edm::ESHandle<GeometricSearchTracker>& tkLayout) : theGeometricSearchTracker(tkLayout) {
   epsilon_ = 0.07; 

}

/* return compatible layers for given trajectory state  */ 
vector<const DetLayer*> 
DirectTrackerNavigation::compatibleLayers( const FreeTrajectoryState& fts,
                                        PropagationDirection dir ) const {

  bool inOut = outward(fts);

  vector<const DetLayer*> output;

  // check eta of DetLayers for compatibility

  if (inOut) { 

      inOutPx(fts,output);

      if ( fts.position().eta() > 1.55) inOutFPx(fts,output);
      else if ( fts.position().eta() < -1.55) inOutBPx(fts,output);

      if ( fabs(fts.position().eta()) < 1.67) inOutTIB(fts,output); 

      if ( fts.position().eta() > 1.17 ) inOutFTID(fts,output);
      else if ( fts.position().eta() < -1.17 ) inOutBTID(fts,output);

      if ( fabs(fts.position().eta()) < 1.35) inOutTOB(fts,output);

      if ( fts.position().eta() > 0.97 ) inOutFTEC(fts,output);
      else if ( fts.position().eta() < -0.97 )  inOutBTEC(fts,output);

   } else {
      LogTrace("Muon|RecoMuon|DirectionTrackerNavigation")<<"No implementation for inward state at this moment. ";
   }

  if ( dir == oppositeToMomentum ) std::reverse(output.begin(),output.end());

  return output;
}


void DirectTrackerNavigation::inOutPx(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<BarrelDetLayer*>& barrel = theGeometricSearchTracker->pixelBarrelLayers();

  for (vector<BarrelDetLayer*>::const_iterator iter_B = barrel.begin(); iter_B != barrel.end(); iter_B++){

      if ( checkCompatible(fts,(*iter_B))) {
      output.push_back((*iter_B));
      }
  }
}

void DirectTrackerNavigation::inOutTIB(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<BarrelDetLayer*>& barrel = theGeometricSearchTracker->tibLayers();

  for (vector<BarrelDetLayer*>::const_iterator iter_B = barrel.begin(); iter_B != barrel.end(); iter_B++){

      if ( checkCompatible(fts,(*iter_B))) {
      output.push_back((*iter_B));
      }
  }
}

void DirectTrackerNavigation::inOutTOB(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<BarrelDetLayer*>& barrel = theGeometricSearchTracker->tobLayers();

  for (vector<BarrelDetLayer*>::const_iterator iter_B = barrel.begin(); iter_B != barrel.end(); iter_B++){

      if ( checkCompatible(fts,(*iter_B))) {
      output.push_back((*iter_B));
      }
  }
}

void DirectTrackerNavigation::inOutFPx(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<ForwardDetLayer*>& forward = theGeometricSearchTracker->posPixelForwardLayers();
  for (vector<ForwardDetLayer*>::const_iterator iter_E = forward.begin(); iter_E != forward.end();
         iter_E++){
      if ( checkCompatible(fts,(*iter_E))) {
        output.push_back((*iter_E));
      }
    }
}

void DirectTrackerNavigation::inOutFTID(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<ForwardDetLayer*>& forward = theGeometricSearchTracker->posTidLayers();
  for (vector<ForwardDetLayer*>::const_iterator iter_E = forward.begin(); iter_E != forward.end();
         iter_E++){
      if ( checkCompatible(fts,(*iter_E))) {
        output.push_back((*iter_E));
      }
    }
}

void DirectTrackerNavigation::inOutFTEC(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<ForwardDetLayer*>& forward = theGeometricSearchTracker->posTecLayers();
  for (vector<ForwardDetLayer*>::const_iterator iter_E = forward.begin(); iter_E != forward.end();
         iter_E++){
      if ( checkCompatible(fts,(*iter_E))) {
        output.push_back((*iter_E));
      }
    }
}

void DirectTrackerNavigation::inOutBPx(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<ForwardDetLayer*>& forward = theGeometricSearchTracker->negPixelForwardLayers();
  for (vector<ForwardDetLayer*>::const_iterator iter_E = forward.begin(); iter_E != forward.end();
         iter_E++){
      if ( checkCompatible(fts,(*iter_E))) {
        output.push_back((*iter_E));
      }
    }
}

void DirectTrackerNavigation::inOutBTID(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<ForwardDetLayer*>& forward = theGeometricSearchTracker->negTidLayers();
  for (vector<ForwardDetLayer*>::const_iterator iter_E = forward.begin(); iter_E != forward.end();
         iter_E++){
      if ( checkCompatible(fts,(*iter_E))) {
        output.push_back((*iter_E));
      }
    }
}

void DirectTrackerNavigation::inOutBTEC(const FreeTrajectoryState& fts, vector<const DetLayer*>& output) const {

  const vector<ForwardDetLayer*>& forward = theGeometricSearchTracker->negTecLayers();
  for (vector<ForwardDetLayer*>::const_iterator iter_E = forward.begin(); iter_E != forward.end();
         iter_E++){
      if ( checkCompatible(fts,(*iter_E))) {
        output.push_back((*iter_E));
      }
    }
}

bool DirectTrackerNavigation::checkCompatible(const FreeTrajectoryState& fts,const BarrelDetLayer* dl) const {

  float eta0 = fts.position().eta();
//  float etam = fts.momentum().eta();

  const BoundCylinder& bc = dl->specificSurface();
  float radius = bc.radius();
  float length = bc.bounds().length()/2.;

  float eta = calculateEta(radius, length);

  return ( fabs(eta0) <= fabs(eta)+epsilon_ );

}

bool DirectTrackerNavigation::checkCompatible(const FreeTrajectoryState& fts,const ForwardDetLayer* dl) const {

  float eta0 = fts.position().eta();
//  float etam = fts.momentum().eta();

  const BoundDisk& bd = dl->specificSurface();

  float outRadius = bd.outerRadius();
  float inRadius = bd.innerRadius();
  float z = bd.position().z();

  float etaOut = calculateEta(outRadius, z);
  float etaIn = calculateEta(inRadius, z);
 
  if ( eta0 > 0 ) return ( eta0 > ( etaOut - epsilon_) && eta0 < (etaIn + epsilon_) );
  else return ( eta0 < (etaOut + epsilon_ ) && eta0 > ( etaIn - epsilon_ ) );

}

bool DirectTrackerNavigation::outward(const FreeTrajectoryState& fts) const {
 
  return (fts.position().basicVector().dot(fts.momentum().basicVector())>0);

}

/// calculate pseudorapidity from r and z
float DirectTrackerNavigation::calculateEta(float r, float z) const {

  if ( z > 0 ) return -log((tan(atan(r/z)/2.)));
  return log(-(tan(atan(r/z)/2.)));

}
