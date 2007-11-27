#ifndef RecoMuon_TrackerSeedGenerator_DualByHitFractionTSG_H
#define RecoMuon_TrackerSeedGenerator_DualByHitFractionTSG_H

#include "RecoMuon/TrackerSeedGenerator/interface/SeparatingTSG.h"

class DualByHitFractionTSG : public SeparatingTSG{
 public:
  DualByHitFractionTSG(const edm::ParameterSet &pset);

  uint selectTSG(const TrackCand&, const TrackingRegion&);

 private:
  std::string theCategory;
  double theFraction;
  double theRecHitAverage;
};

#endif
