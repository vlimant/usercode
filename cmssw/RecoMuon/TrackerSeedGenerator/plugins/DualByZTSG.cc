#include "RecoMuon/TrackerSeedGenerator/interface/DualByZTSG.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

DualByZTSG::DualByZTSG(const edm::ParameterSet &pset) : SeparatingTSG(pset){
  theCategory ="DualByZTSG";
  theZSeparation = pset.getParameter<double>("zSeparation");
  if (nTSGs()!=2)
    {edm::LogError(theCategory)<<"not two seed generators provided";}
}

uint DualByZTSG::selectTSG(const TrackCand & muonTrackCand, const TrackingRegion& region)
{
  LogDebug(theCategory)<<"|vertex z|=|"<<muonTrackCand.second->vz()<<"|"
		       <<" compared to: "<<theZSeparation;
  return (fabs(muonTrackCand.second->vz()) > theZSeparation);
}
    
