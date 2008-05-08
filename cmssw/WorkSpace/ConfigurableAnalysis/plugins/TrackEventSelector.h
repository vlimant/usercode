#ifndef _TrackEventSelector_H
#define _TrackEventSelector_H

#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "PhysicsTools/Utilities/interface/StringCutObjectSelector.h"

class TrackEventSelector : public SusyEventSelector {
 public:
  TrackEventSelector(const edm::ParameterSet& pset) :
    SusyEventSelector(pset)
    {
      StringCutObjectSelector<reco::Track> f("pt > 15.0");
    }
 private:
  StringCutObjectSelector<reco::Track> f
};

#endif
