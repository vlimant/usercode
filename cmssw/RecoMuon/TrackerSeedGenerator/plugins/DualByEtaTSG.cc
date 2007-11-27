#include "RecoMuon/TrackerSeedGenerator/interface/DualByEtaTSG.h"

DualByEtaTSG::DualByEtaTSG(const edm::ParameterSet &pset) : SeparatingTSG(pset){
  theCategory ="DualByEtaTSG";
  theEtaSeparation = pset.getParameter<double>("etaSeparation");
  if (nTSGs()!=2)
    {edm::LogError("DualTSG")<<"not two seed generators provided";}
}

uint DualByEtaTSG::selectTSG(const TrackCand & muonTrackCand, const TrackingRegion& region)
{
  LogDebug(theCategory)<<"|eta|=|"<<muonTrackCand.second->eta()<<"|"
		       <<" compared to: "<<theEtaSeparation;
  return (fabs(muonTrackCand.second->eta()) > theEtaSeparation);
}
    
