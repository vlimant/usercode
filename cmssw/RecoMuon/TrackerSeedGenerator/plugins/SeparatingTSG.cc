#include "RecoMuon/TrackerSeedGenerator/interface/SeparatingTSG.h"

SeparatingTSG::SeparatingTSG(const edm::ParameterSet &pset):CompositeTSG(pset){}

SeparatingTSG::~SeparatingTSG(){}

void SeparatingTSG::trackerSeeds(const TrackCand & muonTrackCand, const TrackingRegion& region, std::vector<TrajectorySeed> & result){
  uint sel = selectTSG(muonTrackCand,region);
  LogDebug(theCategory)<<"choosing: "<<theNames[sel]<<", at index ["<<sel<<"]";
  theTSGs[sel]->trackerSeeds(muonTrackCand,region,result);
}
