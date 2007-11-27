#include "RecoMuon/TrackerSeedGenerator/interface/DualByHitFractionTSG.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

DualByHitFractionTSG::DualByHitFractionTSG(const edm::ParameterSet &pset) : SeparatingTSG(pset){
  theCategory = "DualByHitFractionTSG";
  theFraction = pset.getParameter<double>("fraction");
  theRecHitAverage = pset.getParameter<double>("recHitAverage");
  if (theRecHitAverage==0)
    {edm::LogError(theCategory)<<"average hit is zero. there will be floating point exception.";}
  if (nTSGs()!=2)
    {edm::LogError(theCategory)<<"not two seed generators provided";}
}

uint DualByHitFractionTSG::selectTSG(const TrackCand & muonTrackCand, const TrackingRegion& region)
{
  double fraction = muonTrackCand.second-> recHitsSize()/theRecHitAverage;
  LogDebug(theCategory)<<"nRecHit: "<<muonTrackCand.second-> recHitsSize()
		       <<" average : "<<theRecHitAverage
		       <<" fraction: "<<fraction
		       <<" compared to: "<<theFraction;
return (fraction > theFraction);
}
    
