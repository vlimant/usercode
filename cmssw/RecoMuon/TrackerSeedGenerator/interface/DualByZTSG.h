#ifndef RecoMuon_TrackerSeedGenerator_DualByZTSG_H
#define RecoMuon_TrackerSeedGenerator_DualByZTSG_H

#include "RecoMuon/TrackerSeedGenerator/interface/SeparatingTSG.h"


class DualByZTSG : public SeparatingTSG{
 public:
  DualByZTSG(const edm::ParameterSet &pset);

  uint selectTSG(const TrackCand&, const TrackingRegion&);
 private:
  std::string theCategory;
  double theZSeparation;
};

#endif
